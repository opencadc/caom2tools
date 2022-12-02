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

from pathlib import Path

from astropy.io import fits
from cadcdata import FileInfo
from cadcutils import exceptions
from caom2utils import data_util
from os.path import join

import pytest
from unittest.mock import ANY, Mock, patch
from . import test_fits2caom2


TEST_HEADERS = b"""SIMPLE  =                    T / Written by IDL:  Fri Oct  6 01:48:35 2017
BITPIX  =                  -32 / Bits per pixel
NAXIS   =                    2 / Number of dimensions
NAXIS1  =                 2048 /
NAXIS2  =                 2048 /
DATATYPE= 'REDUC   '           /Data type, SCIENCE/CALIB/REJECT/FOCUS/TEST
END"""


def test_get_file_type():
    vals = {
        'abc.cat': 'text/plain',
        'abc.gif': 'image/gif',
        'abc.png': 'image/png',
        'abc.jpg': 'image/jpeg',
        'abc.tar.gz': 'application/x-tar',
        'abc.csv': 'text/csv',
        'abc.hdf5': 'application/x-hdf5',
        'abc.fits': 'application/fits',
    }
    for key, value in vals.items():
        assert (
            data_util.get_file_type(key) == value
        ), f'wrong type {data_util.get_file_type(key)} for {key}'


@patch('caom2utils.data_util.StorageInventoryClient')
def test_storage_inventory_client(cadc_client_mock):
    test_subject = Mock(autospec=True)
    test_uri = 'cadc:TEST/test_file.fits'
    test_working_directory = Path(test_fits2caom2.TESTDATA_DIR)
    test_fqn = test_working_directory / 'test_file.fits'
    if test_fqn.exists():
        test_fqn.unlink()

    def info_si_mock(ignore):
        return FileInfo(id=test_uri, file_type='application/fits',
                        md5sum='abc', size=42)

    def get_si_mock(ignore2, dest, **kwargs):
        fhead = kwargs.get('fhead')
        if fhead:
            dest.write(TEST_HEADERS)
        else:
            test_fqn.write_text('StorageInventoryClient')

    cadc_client_mock.return_value.cadcinfo.side_effect = info_si_mock
    cadc_client_mock.return_value.cadcget.side_effect = get_si_mock
    cadc_client_mock.return_value.cadcput = Mock(autospec=True)
    cadc_client_mock.return_value.cadcremove = Mock(autospec=True)

    test_wrapper = data_util.StorageClientWrapper(subject=test_subject)
    assert test_wrapper is not None, 'ctor failure'

    # info
    test_result = test_wrapper.info(test_uri)
    _check_info_result(test_result)

    # get_head
    test_result = test_wrapper.get_head(test_uri)
    _check_header_result(test_result)

    # get
    test_wrapper.get(test_working_directory, test_uri)
    _check_get_result(test_fqn)

    # put
    test_wrapper.put(test_working_directory, test_uri)
    _check_put_result(cadc_client_mock.return_value.cadcput)

    # delete
    test_wrapper.remove(test_uri)
    assert cadc_client_mock.return_value.cadcremove.called, 'remove call'
    cadc_client_mock.return_value.cadcremove.assert_called_with(
        test_uri
    ), 'wrong remove args'

    cadc_client_mock.return_value.cadcinfo.side_effect = (
        exceptions.UnexpectedException('cadcinfo')
    )
    cadc_client_mock.return_value.cadcget.side_effect = (
        exceptions.UnexpectedException('cadcget')
    )
    cadc_client_mock.return_value.cadcput.side_effect = (
        exceptions.UnexpectedException('cadcput')
    )
    _fail_mock(test_wrapper, test_uri, test_working_directory)

    cadc_client_mock.return_value.cadcremove.side_effect = (
        exceptions.UnexpectedException('cadcremove')
    )
    with pytest.raises(exceptions.UnexpectedException):
        test_wrapper.remove(test_uri)

    cadc_client_mock.return_value.cadcinfo.side_effect = (
        exceptions.NotFoundException('cadcinfo')
    )
    test_result = test_wrapper.info(test_uri)
    assert test_result is None, 'expected when not found'


@patch('caom2utils.data_util.StorageInventoryClient')
def test_si_tracking(client_mock):
    test_subject = Mock(autospec=True)
    test_metrics = Mock(autospec=True)

    def _get(working_directory, uri):
        raise exceptions.UnexpectedException
    client_mock.return_value.cadcget.side_effect = _get
    client_mock.return_value.cadcremove.side_effect = Mock()

    test_wrapper = data_util.StorageClientWrapper(subject=test_subject, metrics=test_metrics)
    assert test_wrapper is not None, 'ctor failure'

    # test metrics failure
    with pytest.raises(exceptions.UnexpectedException):
        test_wrapper.get('/tmp', 'cadc:TEST/abc.fits')
    assert test_metrics.observe_failure.called, 'expect observe_failure call'
    test_metrics.observe_failure.assert_called_with(
        'get', 'si', 'cadc:TEST/abc.fits'
    )

    # test metrics success
    test_wrapper.remove('cadc:TEST/abc.fits')
    assert test_metrics.observe.called, 'expect observe call'
    test_metrics.observe.assert_called_with(
        ANY, ANY, None, 'remove', 'si', 'cadc:TEST/abc.fits'
    )


def test_clean_headers():
    test_input = """ilename: S20141130S0001.fits.bz2

AstroData Tags: {'CAL', 'RAW', 'AZEL_TARGET', 'DARK', 'F2', 'NON_SIDEREAL'}


--- PHU ---
SIMPLE  =                    T / file does conform to FITS standard
BITPIX  =                   16 / number of bits per data pixel
NAXIS   =                    0 / number of data axes
EXTEND  =                    T / FITS dataset may contain extensions
COMMENT   FITS (Flexible Image Transport System) format defined in 'Astronomy
COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H

--- HDU 0 ---
XTENSION= 'IMAGE   '           / IMAGE extension
BITPIX  =                   32 / number of bits per data pixel
NAXIS   =                    3 / number of data axes
NAXIS1  =                 2048 / length of data axis 1
NAXIS2  =                 2048 / length of data axis 2
NAXIS3  =                    1 / length of data axis 3
    """
    test_result = data_util.make_headers_from_string(test_input)
    assert test_result is not None, 'expect a result'
    assert len(test_result) == 2, 'expect two headers'


def test_unicode_decode_error():
    test_fqn = join(test_fits2caom2.TESTDATA_DIR, 'time_axes.fits')
    result = data_util.get_local_file_headers(test_fqn)
    assert result is not None, 'expect retry using a different method'


def test_get_file_encoding():
    test_subjects = {
        'abc.fits': None,
        'abc.fits.gz': 'gzip',
        'abc.fits.fz': 'x-fits'
    }
    for test_subject in test_subjects.keys():
        test_result = data_util.get_file_encoding(test_subject)
        assert (
            test_result == test_subjects.get(test_subject)
        ), f'got wrong extension {test_result} for {test_subject}'


def _check_get_result(test_fqn):
    assert test_fqn.exists(), 'expected file creation'


def _check_header_result(result):
    assert result is not None, 'expect a result'
    assert len(result) == 1, 'wrong header count'
    assert isinstance(result[0], fits.Header)
    assert result[0].get('DATATYPE') == 'REDUC', 'not quite correct'


def _check_info_result(result):
    assert result is not None, 'expect a result'
    assert result.size == 42, 'wrong size'
    assert result.file_type == 'application/fits', 'wrong file type'
    assert result.md5sum == 'abc', 'wrong md5sum'


def _check_put_result(client_mock):
    assert client_mock.called, 'expect put mock call'
    try:
        client_mock.assert_called_with(
            'TEST',
            'test_file.fits',
            archive_stream='default',
            mime_type='application/fits',
            mime_encoding=None,
            md5_check=True,
        ), 'wrong put args call'
    except AssertionError:
        client_mock.assert_called_with(
            'cadc:TEST/test_file.fits',
            src=f'{test_fits2caom2.TESTDATA_DIR}/test_file.fits',
            replace=True,
            file_type='application/fits',
            file_encoding=None,
            md5_checksum='3c66ee2cb6e0c2cfb5cd6824d353dc11',
        )


def _fail_mock(test_wrapper, test_uri, test_working_directory):
    with pytest.raises(exceptions.UnexpectedException):
        test_wrapper.info(test_uri)

    with pytest.raises(exceptions.UnexpectedException):
        test_wrapper.get_head(test_uri)

    with pytest.raises(exceptions.UnexpectedException):
        test_wrapper.get(test_working_directory, test_uri)

    with pytest.raises(exceptions.UnexpectedException):
        test_wrapper.put(test_working_directory, test_uri)
