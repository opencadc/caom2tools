# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2018.                            (c) 2018.
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
#  $Revision: 4 $
#
# ***********************************************************************
#

from cadcdata import FileInfo
from caom2utils import legacy, caom2blueprint, data_util
from caom2 import ObservationReader, ObservationWriter
from caom2.diff import get_differences

import glob
import logging
import os
import sys
import tempfile
from unittest.mock import patch, Mock
from urllib.parse import urlparse

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TESTDATA_DIR = os.path.join(THIS_DIR, 'data')


def raise_exit_error():
    yield SystemExit("Raised")


@patch('sys.exit', Mock(side_effect=ImportError))
def test_differences(directory):
    """
    Note: This tests is parametrized from conftest.py file. Directories
    of the form TESTDATA_DIR/*/* become parameters to this test.
    This test assumes a directory contains the config, default,
    override, input FITS files (*.header) and expected observation (*.xml)
    files
    """
    expected_fname = _get_file('xml', directory)
    assert expected_fname
    assert len(expected_fname) == 1
    expected = _read_observation(expected_fname[0])  # expected observation
    assert len(expected.planes) == 1
    prod_id = [p.product_id for p in expected.planes.values()][0]
    product_id = f'--productID {prod_id}'
    collection_id = expected.collection
    data_files = _get_files(
        ['header', 'png', 'gif', 'cat', 'fits', 'h5'], directory)
    assert data_files

    file_meta = _get_uris(collection_id, data_files, expected)
    assert file_meta

    data_files_parameter = _get_data_files_parameter(data_files)

    config = _get_parameter('config', directory)
    if config is None:
        blueprints = _get_multi_parameter('blueprint', directory)
        assert blueprints
        module = _get_parameter('module', directory)
        cardinality = _get_cardinality(directory)
        if module is not None:
            inputs = f'{blueprints} {module}'
        else:
            inputs = blueprints
        application = '{} {} '.format('caom2gen', data_files_parameter)
        app_cmd = caom2blueprint.caom2gen
    else:
        defaults = _get_parameter('default', directory)
        assert defaults
        overrides = _get_parameter('override', directory)
        assert overrides
        inputs = f'{config} {defaults} {overrides}'
        application = '{} {}'.format('caom2blueprint', data_files_parameter)
        app_cmd = legacy.main_app
        temp = ' '.join(file_meta[0])
        cardinality = f'{product_id} {temp}'
        # return  # TODO shorter testing cycle

    with patch('caom2utils.data_util.StorageInventoryClient') as \
            swc_si_mock,\
            patch('cadcutils.net.ws.WsCapabilities.get_access_url',
                  autospec=True) as cap_mock,\
            patch('caom2utils.caom2blueprint.get_vos_headers') as gvh_mock, \
            patch('caom2utils.caom2blueprint._get_vos_meta') as gvm_mock, \
            patch('caom2utils.data_util.get_local_headers_from_fits') as \
            header_mock:
        def info_mock(uri):
            if uri.startswith('vos'):
                archive = uri.split('/')[-2]
            elif uri.startswith('/'):
                archive = uri.split('/')[-3].upper()
            else:
                archive = uri.split(':')[1].split('/')[0]
            file_id = uri.split('/')[-1]
            temp = file_meta[1][(archive, file_id)]
            logging.error(temp)
            return temp

        def _get_vos_headers(uri, subject=None):
            if uri.startswith('vos'):
                fname = data_files_parameter.split()[1].strip()
                fits_header = open(fname).read()
                return data_util.make_headers_from_string(fits_header)
            else:
                return None

        def _vos_client_meta(subject, uri):
            return FileInfo(id=uri,
                            md5sum='5b00b00d4b06aba986c3663d09aa581f',
                            size=682560,
                            file_type='application/fits')

        def _header(fqn):
            # during operation, want to use astropy on FITS files
            # but during testing want to use headers and built-in Python file
            # operations
            from urllib.parse import urlparse
            from astropy.io import fits
            file_uri = urlparse(fqn)
            try:
                fits_header = open(file_uri.path).read()
                headers = data_util.make_headers_from_string(fits_header)
            except UnicodeDecodeError:
                hdulist = fits.open(fqn, memmap=True, lazy_load_hdus=True)
                hdulist.verify('fix')
                hdulist.close()
                headers = [h.header for h in hdulist]
            return headers

        swc_si_mock.return_value.cadcinfo.side_effect = info_mock
        swc_si_mock.cadcget.return_value = []
        data_util.get_local_file_info.side_effect = info_mock
        gvh_mock.side_effect = _get_vos_headers
        gvm_mock.side_effect = _vos_client_meta
        cap_mock.return_value = 'https://localhost'
        header_mock.side_effect = _header

        temp = tempfile.NamedTemporaryFile()
        sys.argv = ('{} -o {} --no_validate --observation {} {} {} {} '
                    '--resource-id ivo://cadc.nrc.ca/test'.format(
                        application, temp.name,
                        expected.collection, expected.observation_id,
                        inputs, cardinality)).split()
        print(sys.argv)
        app_cmd()
    actual = _read_observation(temp.name)  # actual observation
    _write_observation(actual)
    _compare_observations(expected, actual, directory)


def _get_cardinality(directory):
    # TODO - read this from an aptly named file in the directory
    # The blueprints are named to reverse sort so that this
    # alignment of product id / artifact URI works
    if '/cfhtsg/' in directory:
        return '--lineage ' \
               'MegaPipe.080.156.Z.MP9801/ad:CFHTSG/' \
               'MegaPipe.080.156.Z.MP9801.weight.fits ' \
               'MegaPipe.080.156.Z.MP9801/ad:CFHTSG/' \
               'MegaPipe.080.156.Z.MP9801.fits ' \
               'MegaPipe.080.156.Z.MP9801/ad:CFHTSG/' \
               'MegaPipe.080.156.Z.MP9801.fits.gif'
    elif '/omm/' in directory:
        if 'SCIRED' in directory:
            return '--lineage Cdemo_ext2_SCIRED/ad:OMM/' \
                   'Cdemo_ext2_SCIRED.fits.gz'
        else:
            return '--lineage C190531_0432_SCI/ad:OMM/' \
                   'C190531_0432_SCI.fits.gz'
    elif 'apass/catalog' in directory:
        return '--lineage catalog/vos://cadc.nrc.ca!vospace/CAOMworkshop/' \
               'Examples/DAO/dao_c122_2016_012725.fits'
    elif 'taos_' in directory:
        if 'def' in directory:
            return '--lineage def/cadc:def/def.h5'
        else:
            return '--lineage star04239531/cadc:TAOSII/taos2_20220201T201317Z_star04239531.h5'
    else:
        return ''


def _get_common(fnames):
    common = os.path.basename(fnames[0])
    for jj in fnames:
        rhs = os.path.basename(jj)
        for kk in fnames:
            lhs = os.path.basename(kk)
            ii = 0
            while ii < len(rhs):
                if rhs[ii] == '.' or rhs[ii] != lhs[ii]:
                    if len(rhs[0:ii]) != 0 and len(rhs[0:ii]) < len(common):
                        common = rhs[0:ii]
                    break
                else:
                    ii += 1
    return common


def _get_subdirs(dir_name):
    return [name for name in os.listdir(dir_name) if
            os.path.isdir(os.path.join(dir_name, name))]


def _get_parameter(extension, dir_name):
    fnames = _get_file(extension, dir_name)
    if fnames:
        result = f'--{extension} {fnames[0]}'
        return result
    else:
        return None


def _get_multi_parameter(extension, dir_name):
    fnames = _get_file(extension, dir_name)
    if fnames:
        result = f'--{extension}'
        for fname in fnames:
            result = f'{result} {fname}'
        return result
    else:
        return None


def _get_data_files_parameter(fnames):
    if fnames:
        result = '--local'
        for fname in fnames:
            result = f'{result} {fname}'
        return result
    else:
        return None


def _get_uris(collection, fnames, obs):
    # NOTE: this function makes the assumption that the collection is
    # the same with the AD archive, which in many cases is not true (CFHT
    # for example has multiple collections)
    uris = []
    file_meta = {}
    if fnames:
        for fname in fnames:
            f = os.path.basename(fname).replace('.header', '')
            for p in obs.planes.values():
                for a in p.artifacts.values():
                    if (f'ad:{collection}/{f}' in a.uri or
                            (a.uri.startswith('vos') and f in a.uri)):
                        uris.append(a.uri)
                        meta = FileInfo(id=a.uri,
                                        file_type=a.content_type,
                                        size=a.content_length,
                                        md5sum=a.content_checksum.checksum)
                        file_url = urlparse(a.uri)
                        if file_url.scheme not in ['ad', 'vos', 'cadc']:
                            # TODO add hook to support other service providers
                            raise NotImplementedError(
                                'Only ad, vos type URIs supported')
                        archive, file_id = file_url.path.split('/')[-2:]
                        file_meta[(archive, file_id)] = meta
        return uris, file_meta
    else:
        return None


def _get_file(pattern, dir_name):
    files = glob.glob(f'{dir_name}/*.{pattern}')
    if files:
        files.sort(reverse=True)
        return files
    else:
        return None


def _get_files(patterns, dir_name):
    files = []
    for pattern in patterns:
        temp = _get_file(pattern, dir_name)
        if temp:
            files.extend(temp)
    return files


def _compare_observations(expected, actual, output_dir):

    result = get_differences(expected, actual, 'Observation')
    if result:
        tmp = '\n'.join([r for r in result])
        msg = f'Differences found observation {expected.observation_id} in ' \
              f'{output_dir}\n{tmp}'
        _write_observation(actual)
        raise AssertionError(msg)
    else:
        logging.info('Observation {} in {} match'.format(
            expected.observation_id, output_dir))


def _read_observation(fname):
    reader = ObservationReader(False)
    result = reader.read(fname)
    return result


def _write_observation(obs):
    writer = ObservationWriter(True, False, 'caom2',
                               'http://www.opencadc.org/caom2/xml/v2.4')
    writer.write(obs, './x.xml')
