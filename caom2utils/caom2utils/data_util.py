# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2021.                            (c) 2021.
#  Government of Canada                 Gouvernement du Canada
#  National Research Council            Conseil national de recherches
#  Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
#  All rights reserved                  Tous droits réservés
#
#  NRC disclaims any warranties,        Le CNRC dénie toute garantie
#  expressed, implied, or               énoncée, implicite ou légale,
#  statutory, of any kind with          de quelque nature que ce
#  respect to the software,             soit, concernant le logiciel,
#  including without limitation         y compris sans restriction
#  any warranty of merchantability      toute garantie de valeur
#  or fitness for a particular          marchande ou de pertinence
#  purpose. NRC shall not be            pour un usage particulier.
#  liable in any event for any          Le CNRC ne pourra en aucun cas
#  damages, whether direct or           être tenu responsable de tout
#  indirect, special or general,        dommage, direct ou indirect,
#  consequential or incidental,         particulier ou général,
#  arising from the use of the          accessoire ou fortuit, résultant
#  software.  Neither the name          de l'utilisation du logiciel. Ni
#  of the National Research             le nom du Conseil National de
#  Council of Canada nor the            Recherches du Canada ni les noms
#  names of its contributors may        de ses  participants ne peuvent
#  be used to endorse or promote        être utilisés pour approuver ou
#  products derived from this           promouvoir les produits dérivés
#  software without specific prior      de ce logiciel sans autorisation
#  written permission.                  préalable et particulière
#                                       par écrit.
#
#  This file is part of the             Ce fichier fait partie du projet
#  OpenCADC project.                    OpenCADC.
#
#  OpenCADC is free software:           OpenCADC est un logiciel libre ;
#  you can redistribute it and/or       vous pouvez le redistribuer ou le
#  modify it under the terms of         modifier suivant les termes de
#  the GNU Affero General Public        la “GNU Affero General Public
#  License as published by the          License” telle que publiée
#  Free Software Foundation,            par la Free Software Foundation
#  either version 3 of the              : soit la version 3 de cette
#  License, or (at your option)         licence, soit (à votre gré)
#  any later version.                   toute version ultérieure.
#
#  OpenCADC is distributed in the       OpenCADC est distribué
#  hope that it will be useful,         dans l’espoir qu’il vous
#  but WITHOUT ANY WARRANTY;            sera utile, mais SANS AUCUNE
#  without even the implied             GARANTIE : sans même la garantie
#  warranty of MERCHANTABILITY          implicite de COMMERCIALISABILITÉ
#  or FITNESS FOR A PARTICULAR          ni d’ADÉQUATION À UN OBJECTIF
#  PURPOSE.  See the GNU Affero         PARTICULIER. Consultez la Licence
#  General Public License for           Générale Publique GNU Affero
#  more details.                        pour plus de détails.
#
#  You should have received             Vous devriez avoir reçu une
#  a copy of the GNU Affero             copie de la Licence Générale
#  General Public License along         Publique GNU Affero avec
#  with OpenCADC.  If not, see          OpenCADC ; si ce n’est
#  <http://www.gnu.org/licenses/>.      pas le cas, consultez :
#                                       <http://www.gnu.org/licenses/>.
#
#  : 4 $
#
# ***********************************************************************
#

import logging
import traceback
from datetime import datetime, timezone
from hashlib import md5
from io import BytesIO
from os import chdir, getcwd, path, stat

from astropy.io import fits
from urllib.parse import urlparse
from cadcdata import FileInfo, StorageInventoryClient
from cadcutils import exceptions


__all__ = [
    'get_file_encoding',
    'get_file_type',
    'get_local_file_headers',
    'get_local_file_info',
    'get_local_headers_from_fits',
    'make_headers_from_string',
    'StorageClientWrapper',
]


class StorageClientWrapper:
    """
    Wrap the metrics collection with StorageInventoryClient.
    """

    def __init__(self, subject, resource_id='ivo://cadc.nrc.ca/uvic/minoc', metrics=None):
        """
        :param subject: net.Subject instance for authentication and authorization
        :param resource_id: str identifies the StorageInventoryClient endpoint. Defaults to the installation closest to
            most of the current invocations.
        :param metrics: caom2pipe.manaage_composable.Metrics instance. If set, will track execution times, by action,
            from the beginning of the method invocation to the end of the method invocation, success or failure.
            Defaults to None, because fits2caom2 is a stand-alone application.
        """
        self._cadc_client = StorageInventoryClient(subject=subject, resource_id=resource_id)
        self._metrics = metrics
        self._logger = logging.getLogger(self.__class__.__name__)

    def _add_fail_metric(self, action, name):
        """Single location for the check for a self._metrics member in the failure case."""
        if self._metrics is not None:
            self._metrics.observe_failure(action, 'si', name)

    def _add_metric(self, action, name, start, value):
        """Single location for the check for a self._metrics member in the success case."""
        if self._metrics is not None:
            self._metrics.observe(start, StorageClientWrapper._current(), value, action, 'si', name)

    def get(self, working_directory, uri):
        """
        Retrieve data.
        :param working_directory: str where the file will be retrieved to.
            Assumes the same machine as this function is being called from.
        :param uri: str this is an Artifact URI, representing the file to
            be retrieved.
        """
        self._logger.debug(f'Begin get for {uri} in {working_directory}')
        start = StorageClientWrapper._current()
        try:
            archive, f_name = self._decompose(uri)
            fqn = path.join(working_directory, f_name)
            self._cadc_client.cadcget(uri, dest=fqn)
        except Exception as e:
            self._add_fail_metric('get', uri)
            self._logger.debug(traceback.format_exc())
            raise exceptions.UnexpectedException(
                f'Did not retrieve {uri} because {e}'
            )
        self._add_metric('get', uri, start, stat(fqn).st_size)
        self._logger.debug('End get')

    def get_head(self, uri):
        """
        Retrieve FITS file header data.
        :param uri: str that is an Artifact URI, representing the file for
            which to retrieve headers
        :return: list of fits.Header instances
        """
        self._logger.debug(f'Begin get_head for {uri}')
        start = StorageClientWrapper._current()
        try:
            b = BytesIO()
            b.name = uri
            self._cadc_client.cadcget(uri, b, fhead=True)
            fits_header = b.getvalue().decode('ascii')
            b.close()
            self._add_metric('get_head', uri, start, len(fits_header))
            temp = make_headers_from_string(fits_header)
            self._logger.debug('End get_head')
            return temp
        except Exception as e:
            self._add_fail_metric('get_header', uri)
            self._logger.debug(traceback.format_exc())
            self._logger.error(e)
            raise exceptions.UnexpectedException(
                f'Did not retrieve {uri} header because {e}'
            )

    def info(self, uri):
        """
        Retrieve the descriptive metadata associated with a file.
        :param uri: str that is an Artifact URI, representing the file for
            which to retrieve metadata
        :return: cadcdata.FileInfo instance, no scheme for md5sum
        """
        self._logger.debug(f'Begin info for {uri}')
        try:
            result = self._cadc_client.cadcinfo(uri)
            # make the result look like the other possible ways to
            # obtain metadata
            result.md5sum = result.md5sum.replace('md5:', '')
        except exceptions.NotFoundException:
            self._logger.info(f'cadcinfo:: {uri} not found')
            result = None
        self._logger.debug('End info')
        return result

    def put(self, working_directory, uri):
        """
        Store a file at CADC.
        :param working_directory: str fully-qualified name of where to find the file on the local machine
        :param uri: str that is an Artifact URI, representing the file to be stored at CADC.
        """
        self._logger.debug(f'Begin put for {uri} in {working_directory}')
        start = self._current()
        cwd = getcwd()
        archive, f_name = StorageClientWrapper._decompose(uri)
        fqn = path.join(working_directory, f_name)
        chdir(working_directory)
        try:
            local_meta = get_local_file_info(fqn)
            encoding = get_file_encoding(fqn)
            replace = True
            cadc_meta = self.info(uri)
            if cadc_meta is None:
                replace = False
            self._logger.debug(
                f'uri {uri} src {fqn} replace {replace} file_type {local_meta.file_type} encoding {encoding} '
                f'md5_checksum {local_meta.md5sum}'
            )
            self._cadc_client.cadcput(
                uri,
                src=fqn,
                replace=replace,
                file_type=local_meta.file_type,
                file_encoding=encoding,
                md5_checksum=local_meta.md5sum,
            )
            self._logger.info(f'Stored {fqn} at CADC.')
        except Exception as e:
            self._add_fail_metric('put', uri)
            self._logger.debug(traceback.format_exc())
            self._logger.error(e)
            raise exceptions.UnexpectedException(
                f'Failed to store data with {e}'
            )
        finally:
            chdir(cwd)
        self._add_metric('put', uri, start, local_meta.size)
        self._logger.debug('End put')

    def remove(self, uri):
        """
        Delete a file from CADC storage.
        :param uri: str that is an Artifact URI, representing the file to
            be removed from CADC.
        """
        self._logger.debug(f'Begin remove for {uri}')
        start = StorageClientWrapper._current()
        try:
            self._cadc_client.cadcremove(uri)
        except Exception as e:
            self._add_fail_metric('remove', uri)
            self._logger.debug(traceback.format_exc())
            self._logger.error(e)
            raise exceptions.UnexpectedException(
                f'Did not remove {uri} because {e}'
            )
        self._add_metric('remove', uri, start, value=None)
        self._logger.debug('End remove')

    @staticmethod
    def _current():
        """Encapsulate returning UTC now in microsecond resolution."""
        return datetime.now(tz=timezone.utc).timestamp()

    @staticmethod
    def _decompose(uri):
        temp = urlparse(uri)
        return path.dirname(temp.path), path.basename(temp.path)


def _clean_headers(fits_header):
    """
    Hopefully not Gemini specific.
    Remove invalid cards and add missing END cards after extensions.
    :param fits_header: fits_header a string of keyword/value pairs
    """
    new_header = []
    first_header_encountered = False
    for line in fits_header.split('\n'):
        if len(line.strip()) == 0:
            pass
        elif line.startswith('--- PHU ---'):
            first_header_encountered = True
        elif line.startswith('--- HDU 0'):
            if first_header_encountered:
                new_header.append('END\n')
            else:
                first_header_encountered = True
        elif line.startswith('--- HDU'):
            new_header.append('END\n')
        elif line.strip() == 'END':
            new_header.append('END\n')
        elif '=' not in line and not (line.startswith('COMMENT') or
                                      line.startswith('HISTORY')):
            pass
        else:
            new_header.append(f'{line}\n')
    new_header.append('END\n')
    return ''.join(new_header)


def get_local_headers_from_fits(fqn):
    """Create a list of fits.Header instances from a fits file.
    :param fqn str  fully-qualified name of the FITS file on disk
    :return list of fits.Header instances
    """
    hdulist = fits.open(fqn, memmap=True, lazy_load_hdus=True)
    hdulist.verify('fix')
    hdulist.close()
    headers = [h.header for h in hdulist]
    return headers


def get_local_file_headers(fqn):
    """
    Wrap two different attempts for header retrieval into a single
    function.
    :param fqn: str fully-qualified name of the FITS file on disk
    :return: list of fits.Header instances
    """
    file_uri = urlparse(fqn)
    try:
        fits_header = open(file_uri.path).read()
        headers = make_headers_from_string(fits_header)
    except UnicodeDecodeError:
        headers = get_local_headers_from_fits(file_uri.path)
    return headers


def get_local_file_info(fqn):
    """
    Gets descriptive metadata for a file on disk.
    :param fqn: Fully-qualified name of the file on disk.
    :return: FileInfo, no scheme on the md5sum value.
    """
    s = stat(fqn)
    # copy and paste from cadcdata/storageinventory.py
    hash_md5 = md5()
    with open(fqn, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b''):
            hash_md5.update(chunk)
    meta = FileInfo(
        id=path.basename(fqn),
        size=s.st_size,
        md5sum=hash_md5.hexdigest(),
        file_type=get_file_type(fqn),
    )
    return meta


def get_file_encoding(fqn):
    """Basic header extension to content_encoding lookup."""
    if fqn.endswith('.fits.fz'):
        return 'x-fits'
    elif fqn.endswith('.fits.gz'):
        return 'gzip'
    else:
        return None


def get_file_type(fqn):
    """Basic header extension to content_type lookup."""
    if (fqn.endswith('.header') or fqn.endswith('.txt') or
            fqn.endswith('.cat') or fqn.endswith('.dat')):
        return 'text/plain'
    elif fqn.endswith('.gif'):
        return 'image/gif'
    elif fqn.endswith('.png'):
        return 'image/png'
    elif fqn.endswith('.jpg'):
        return 'image/jpeg'
    elif fqn.endswith('.tar.gz'):
        return 'application/x-tar'
    elif fqn.endswith('.csv'):
        return 'text/csv'
    elif fqn.endswith('.hdf5') or fqn.endswith('.h5'):
        return 'application/x-hdf5'
    else:
        return 'application/fits'


def make_headers_from_string(fits_header):
    """Create a list of fits.Header instances from a string.
    ":param fits_header a string of keyword/value pairs"""
    fits_header = _clean_headers(fits_header)
    delim = 'END\n'
    extensions = \
        [e + delim for e in fits_header.split(delim) if e.strip()]
    headers = [fits.Header.fromstring(e, sep='\n') for e in extensions]
    return headers
