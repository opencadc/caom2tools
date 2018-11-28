# -*- coding: utf-8 -*-
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

import os
import pytest
import sys

from mock import Mock, patch

from caom2 import ProductType, ReleaseType, Artifact, ChecksumURI
from caom2pipe import manage_composable as mc


THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TESTDATA_DIR = os.path.join(THIS_DIR, 'data')


@pytest.mark.skipif(not sys.version.startswith('3.6'),
                    reason='support 3.6 only')
def test_config_class():
    os.getcwd = Mock(return_value=TESTDATA_DIR)
    test_config = mc.Config()
    test_config.get_executors()
    assert test_config is not None
    assert test_config.work_file == 'todo.txt'
    assert test_config.features is not None
    assert test_config.features.supports_composite is False


@pytest.mark.skipif(not sys.version.startswith('3.6'),
                    reason='support 3.6 only')
def test_exec_cmd():
    test_cmd = 'ls /abc'
    with pytest.raises(mc.CadcException):
        mc.exec_cmd(test_cmd)


@pytest.mark.skipif(not sys.version.startswith('3.6'),
                    reason='support 3.6 only')
def test_exec_cmd_redirect():
    fqn = os.path.join(TESTDATA_DIR, 'exec_cmd_redirect.txt')
    if os.path.exists(fqn):
        os.remove(fqn)

    test_cmd = 'ls'
    mc.exec_cmd_redirect(test_cmd, fqn)
    assert os.path.exists(fqn)
    assert os.stat(fqn).st_size > 0


@pytest.mark.skipif(not sys.version.startswith('3.6'),
                    reason='support 3.6 only')
@patch('caom2utils.fits2caom2.CadcDataClient.get_file_info')
def test_compare_checksum(mock_get_file_info):

    # fail case - file doesn't exist
    test_file = os.path.join(TESTDATA_DIR, 'test_omm.fits.gz')
    test_netrc = os.path.join(TESTDATA_DIR, 'test_netrc')
    with pytest.raises(mc.CadcException):
        mc.compare_checksum(test_netrc, 'OMM', test_file)

    # fail case - file exists, different checksum - make a small test file
    test_file = os.path.join(TESTDATA_DIR, 'C111107_0694_SCI.fits')
    f = open(test_file, 'w')
    f.write('test')
    f.close()
    with pytest.raises(mc.CadcException):
        mc.compare_checksum(test_netrc, 'OMM', test_file)


@pytest.mark.skipif(not sys.version.startswith('3.6'),
                    reason='support 3.6 only')
def test_decompose_lineage():
    test_product_id = 'product_id'
    test_uri = 'ad:STARS/galaxies.fits.gz'
    test_lineage = '{}/{}'.format(test_product_id, test_uri)
    actual_product_id, actual_uri = mc.decompose_lineage(test_lineage)
    assert actual_product_id == test_product_id, 'expected {}'.format(
        test_product_id)
    assert actual_uri == test_uri, 'expected {}'.format(test_uri)

    with pytest.raises(mc.CadcException):
        mc.decompose_lineage('')


@pytest.mark.skipif(not sys.version.startswith('3.6'),
                    reason='support 3.6 only')
def test_read_csv_file():
    # bad read
    with pytest.raises(mc.CadcException):
        mc.read_csv_file(None)

    # good read
    test_file_name = os.path.join(TESTDATA_DIR, 'test_csv.csv')
    content = mc.read_csv_file(test_file_name)
    assert content is not None, 'empty results returned'
    assert len(content) == 1, 'missed the comment and the header'
    assert len(content[0]) == 24, 'missed the content'


@pytest.mark.skipif(not sys.version.startswith('3.6'),
                    reason='support 3.6 only')
def test_get_file_meta():
    # None
    with pytest.raises(mc.CadcException):
        mc.get_file_meta(None)

    # non-existent file
    fqn = os.path.join(TESTDATA_DIR, 'abc.txt')
    with pytest.raises(mc.CadcException):
        mc.get_file_meta(fqn)

    # empty file
    fqn = os.path.join(TESTDATA_DIR, 'todo.txt')
    result = mc.get_file_meta(fqn)
    assert result['size'] == 0, result['size']


@pytest.mark.skipif(not sys.version.startswith('3.6'),
                    reason='support 3.6 only')
@patch('cadcdata.core.net.BaseWsClient')
def test_read_file_list_from_archive(basews_mock):

    response = Mock()
    response.status_code.return_value = 200
    basews_mock.return_value.get.return_value = response
    test_config = mc.Config()
    result = mc.read_file_list_from_archive(test_config, 'test_app_name',
                                            '2018-11-18T22:39:56.186443+00:00',
                                            '2018-11-19T22:39:56.186443+00:00')
    assert result is not None
    assert type(result) is list
    assert len(result) == 0


@pytest.mark.skipif(not sys.version.startswith('3.6'),
                    reason='support 3.6 only')
def test_write_to_file():
    content = ['a.txt', 'b.jpg', 'c.fits.gz']
    test_fqn = '{}/test_out.txt'.format(TESTDATA_DIR)
    if os.path.exists(test_fqn):
        os.remove(test_fqn)

    mc.write_to_file(test_fqn, '\n'.join(content))
    assert os.path.exists(test_fqn)


@pytest.mark.skipif(not sys.version.startswith('3.6'),
                    reason='support 3.6 only')
def test_get_lineage():
    result = mc.get_lineage('TEST_COLLECTION', 'TEST_PRODUCT_ID',
                            'TEST_FILE_NAME.fits')
    assert result == 'TEST_PRODUCT_ID/ad:TEST_COLLECTION/TEST_FILE_NAME.fits'


@pytest.mark.skipif(not sys.version.startswith('3.6'),
                    reason='support 3.6 only')
def test_get_artifact_metadata():
    test_fqn = os.path.join(TESTDATA_DIR, 'config.yml')
    test_uri = 'ad:TEST/config.yml'

    # wrong command line parameters
    with pytest.raises(mc.CadcException):
        mc.get_artifact_metadata(test_fqn, ProductType.WEIGHT,
                                 ReleaseType.META)

    # create action
    result = mc.get_artifact_metadata(test_fqn, ProductType.WEIGHT,
                                 ReleaseType.META, uri=test_uri)
    assert result is not None, 'expect a result'
    assert isinstance(result, Artifact), 'expect an artifact'
    assert result.product_type == ProductType.WEIGHT, 'wrong product type'
    assert result.content_length == 255, 'wrong length'
    assert result.content_checksum.uri == \
           'md5:c649725745805d41fc1b601e85400e60', 'wrong checksum'

    # update action
    result.content_checksum = ChecksumURI('md5:abc')
    result = mc.get_artifact_metadata(test_fqn, ProductType.WEIGHT,
                                      ReleaseType.META, artifact=result)
    assert result is not None, 'expect a result'
    assert isinstance(result, Artifact), 'expect an artifact'
    assert result.content_checksum.uri == \
           'md5:c649725745805d41fc1b601e85400e60', 'wrong checksum'


# TODO understand python mocking .... :(
# @pytest.mark.skipif(not sys.version.startswith('3.6'),
#                     reason='support 3.6 only')
# @patch('cadcdata.core.net.BaseWsClient')
# @patch('cadcdata.core.TransferReader')
# @patch('cadcdata.core.CadcDataClient')
# def test_get_cadc_headers(basews_mock, trans_reader_mock, client_mock):
#     with pytest.raises(mc.CadcException):
#         mc.get_cadc_headers('file:GEM/TEST.fits')
#
#     # from cadcdata import transfer
#     # t = transfer.Transfer('ad:GEM/TEST.fits', 'pullFromVoSpace')
#     # p = transfer.Protocol
#     # p.endpoint = Mock()
#     # t.protocols = [p]
#     # trans_reader_mock.return_value.read.return_value = t
#     #
#     # file_content = 'ABCDEFGH12345'
#     # file_chunks = [file_content[i:i + 5].encode()
#     #                for i in range(0, len(file_content), 5)]
#     # response = Mock()
#     # response.headers.get.return_value = 'filename={}.gz'.format('TEST.fits')
#     # response.raw.read.side_effect = file_chunks
#     # response.history = []
#     # response.status_code = 200
#     # response.url = 'someurl'
#     # post_mock = Mock(return_value=response)
#     # basews_mock.return_value.post = post_mock
#     client_mock.return_value.get_file.return_value = ''
#     result = mc.get_cadc_headers('ad:GEM/TEST.fits')
#     assert result is not None
#     assert len(result) == 0
#
# a different approach
#     with patch('caom2utils.fits2caom2.CadcDataClient') as data_client_mock:
#         def get_file_info(archive, file_id):
#             if '_prev' in file_id:
#                 return {'size': 10290,
#                         'md5sum': md5('-37'.encode()).hexdigest(),
#                         'type': 'image/jpeg'}
#             else:
#                 return {'size': 37,
#                         'md5sum': md5('-37'.encode()).hexdigest(),
#                         'type': 'application/fits'}
#         data_client_mock.return_value.get_file_info.side_effect = \
#             get_file_info
