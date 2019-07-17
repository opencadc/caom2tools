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

import six

from datetime import datetime
from mock import Mock, patch

from caom2 import ProductType, ReleaseType, Artifact, ChecksumURI
from caom2 import SimpleObservation
if six.PY3:
    from caom2pipe import manage_composable as mc


PY_VERSION = '3.6'
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_DIR = os.path.join(THIS_DIR, 'data')
TEST_STATE_FILE = os.path.join(TEST_DATA_DIR, 'test_state.yml')
TEST_OBS_FILE = os.path.join(TEST_DATA_DIR, 'test_obs_id.fits.xml')
ISO8601_FORMAT = '%Y-%m-%dT%H:%M:%S.%f'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_read_obs():
    test_subject = mc.read_obs_from_file(TEST_OBS_FILE)
    assert test_subject is not None, 'expect a result'
    assert isinstance(test_subject, SimpleObservation), 'wrong read'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_read_from_file():
    test_subject = mc.read_from_file(TEST_OBS_FILE)
    assert test_subject is not None, 'expect a result'
    assert isinstance(test_subject, list), 'wrong type of result'
    assert len(test_subject) == 8, 'missed some content'
    assert test_subject[0].startswith('<?xml version'), 'read failed'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_build_uri():
    test_subject = mc.build_uri('archive', 'file_name.fits')
    assert test_subject is not None, 'expect a result'
    assert test_subject == 'ad:archive/file_name.fits', 'wrong result'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_query_endpoint():

    with patch('requests.Session.get') as session_get_mock:
        test_result = mc.query_endpoint('https://localhost', timeout=25)
        assert test_result is not None, 'expected result'
        assert session_get_mock.is_called, 'mock not called'
        session_get_mock.assert_called_with('https://localhost', timeout=25)


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_config_class():
    os.getcwd = Mock(return_value=TEST_DATA_DIR)
    mock_root = '/usr/src/app/omm2caom2/omm2caom2/tests/data'
    test_config = mc.Config()
    test_config.get_executors()
    assert test_config is not None
    assert test_config.work_file == 'todo.txt'
    assert test_config.features is not None
    assert test_config.features.supports_composite is False
    assert test_config.working_directory == mock_root, 'wrong dir'
    assert test_config.work_fqn == '{}/todo.txt'.format(mock_root), 'work_fqn'
    assert test_config.netrc_file == '.netrc', 'netrc'
    assert test_config.archive == 'TEST', 'archive'
    assert test_config.collection == 'TEST', 'collection'
    assert test_config.log_file_directory == mock_root, 'logging dir'
    assert test_config.success_fqn == '{}/success_log.txt'.format(mock_root), \
        'success fqn'
    assert test_config.success_log_file_name == 'success_log.txt', \
        'success file'
    assert test_config.failure_fqn == '{}/failure_log.txt'.format(mock_root), \
        'failure fqn'
    assert test_config.failure_log_file_name == 'failure_log.txt', \
        'failure file'
    assert test_config.retry_file_name == 'retries.txt', 'retry file'
    assert test_config.retry_fqn == '{}/retries.txt'.format(mock_root), \
        'retry fqn'
    assert test_config.proxy_file_name == 'test_proxy.pem', 'proxy file name'
    assert test_config.proxy_fqn == '{}/test_proxy.pem'.format(mock_root), \
        'proxy fqn'
    assert test_config.state_file_name == 'state.yml', 'state file name'
    assert test_config.state_fqn == '{}/state.yml'.format(mock_root), \
        'state fqn'
    assert test_config.rejected_directory == TEST_DATA_DIR, \
        'wrong rejected dir'
    assert test_config.rejected_file_name == 'rejected.yml', \
        'wrong rejected file'
    assert test_config.rejected_fqn == '{}/rejected.yml'.format(
        TEST_DATA_DIR), 'wrong rejected fqn'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_exec_cmd():
    test_cmd = 'ls /abc'
    with pytest.raises(mc.CadcException):
        mc.exec_cmd(test_cmd)


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support 3.6 only')
def test_exec_cmd_redirect():
    fqn = os.path.join(TEST_DATA_DIR, 'exec_cmd_redirect.txt')
    if os.path.exists(fqn):
        os.remove(fqn)

    test_cmd = 'ls'
    mc.exec_cmd_redirect(test_cmd, fqn)
    assert os.path.exists(fqn)
    assert os.stat(fqn).st_size > 0


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
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


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_read_csv_file():
    # bad read
    with pytest.raises(mc.CadcException):
        mc.read_csv_file(None)

    # good read
    test_file_name = os.path.join(TEST_DATA_DIR, 'test_csv.csv')
    content = mc.read_csv_file(test_file_name)
    assert content is not None, 'empty results returned'
    assert len(content) == 1, 'missed the comment and the header'
    assert len(content[0]) == 24, 'missed the content'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_get_file_meta():
    # None
    with pytest.raises(mc.CadcException):
        mc.get_file_meta(None)

    # non-existent file
    fqn = os.path.join(TEST_DATA_DIR, 'abc.txt')
    with pytest.raises(mc.CadcException):
        mc.get_file_meta(fqn)

    # empty file
    fqn = os.path.join(TEST_DATA_DIR, 'todo.txt')
    result = mc.get_file_meta(fqn)
    assert result['size'] == 0, result['size']


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
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


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_write_to_file():
    content = ['a.txt', 'b.jpg', 'c.fits.gz']
    test_fqn = '{}/test_out.txt'.format(TEST_DATA_DIR)
    if os.path.exists(test_fqn):
        os.remove(test_fqn)

    mc.write_to_file(test_fqn, '\n'.join(content))
    assert os.path.exists(test_fqn)


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support 3.6 only')
def test_get_lineage():
    result = mc.get_lineage('TEST_COLLECTION', 'TEST_PRODUCT_ID',
                            'TEST_FILE_NAME.fits')
    assert result == 'TEST_PRODUCT_ID/ad:TEST_COLLECTION/TEST_FILE_NAME.fits'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_get_artifact_metadata():
    test_fqn = os.path.join(TEST_DATA_DIR, 'config.yml')
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
    assert result.content_length == 314, 'wrong length'
    assert result.content_checksum.uri == \
        'md5:a75377d8d7cc55464944947c01cef816', 'wrong checksum'

    # update action
    result.content_checksum = ChecksumURI('md5:abc')
    result = mc.get_artifact_metadata(test_fqn, ProductType.WEIGHT,
                                      ReleaseType.META, artifact=result)
    assert result is not None, 'expect a result'
    assert isinstance(result, Artifact), 'expect an artifact'
    assert result.content_checksum.uri == \
        'md5:a75377d8d7cc55464944947c01cef816', 'wrong checksum'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('cadcdata.core.CadcDataClient')
@patch('caom2pipe.manage_composable.Metrics')
def test_data_put(mock_metrics, mock_client):
    mc.data_put(mock_client, TEST_DATA_DIR, 'TEST.fits', 'TEST', 'default',
                metrics=mock_metrics)
    mock_client.put_file.assert_called_with(
        'TEST', 'TEST.fits', archive_stream='default',
        md5_check=True, mime_encoding=None, mime_type=None), 'mock not called'
    assert mock_metrics.observe.is_called, 'mock not called'
    args, kwargs = mock_metrics.observe.call_args
    assert args[2] == 0, 'wrong size'
    assert args[3] == 'put', 'wrong endpoint'
    assert args[4] == 'data', 'wrong service'
    assert args[5] == 'TEST.fits', 'wrong id'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('cadcdata.core.CadcDataClient')
def test_data_get(mock_client):
    test_config = mc.Config()
    test_config.observe_execution = True
    test_metrics = mc.Metrics(test_config)
    with pytest.raises(mc.CadcException):
        mc.data_get(mock_client, TEST_DATA_DIR, 'TEST_get.fits', 'TEST',
                    test_metrics)
    assert len(test_metrics.failures) == 1, 'wrong failures'
    assert test_metrics.failures['data']['get']['TEST_get.fits'] == 1, 'count'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_state():
    test_start = os.path.getmtime(TEST_STATE_FILE)
    with pytest.raises(mc.CadcException):
        test_subject = mc.State('nonexistent')

    test_subject = mc.State(TEST_STATE_FILE)
    assert test_subject is not None, 'expect result'
    test_result = test_subject.get_bookmark('gemini_timestamp')
    assert test_result is not None, 'expect content'
    assert isinstance(test_result, datetime)

    test_subject.save_state('gemini_timestamp', test_result)
    test_end = os.path.getmtime(TEST_STATE_FILE)
    assert test_start != test_end, 'file should be modified'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_make_seconds():
    t1 = '2017-06-26T17:07:21.527+00'
    t1_dt = mc.make_seconds(t1)
    assert t1_dt is not None, 'expect a result'
    assert t1_dt == 1498496841.527, 'wrong result'

    t2 = '2017-07-26T17:07:21.527'
    t2_dt = mc.make_seconds(t2)
    assert t2_dt is not None, 'expect a result'
    assert t2_dt == 1501088841.527, 'wrong result'

    t3 = '16-Jul-2019 09:08'
    t3_dt = mc.make_seconds(t3)
    assert t3_dt is not None, 'expect a result'
    assert t3_dt == 1563268080.0, 'wrong result'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_increment_time():
    t1 = '2017-06-26T17:07:21.527'
    t1_dt = datetime.strptime(t1, mc.ISO_8601_FORMAT)
    result = mc.increment_time(t1_dt, 10)
    assert result is not None, 'expect a result'
    assert result == datetime(2017, 6, 26, 17, 17, 21, 527000),\
        'wrong result'

    t2 = '2017-07-26T17:07:21.527'
    t2_dt = datetime.strptime(t2, mc.ISO_8601_FORMAT)
    result = mc.increment_time(t2_dt, 5)
    assert result is not None, 'expect a result'
    assert result == datetime(2017, 7, 26, 17, 12, 21, 527000),\
        'wrong result'

    with pytest.raises(NotImplementedError):
        mc.increment_time(t2_dt, 23, '%f')


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('cadcdata.core.CadcDataClient')
@patch('caom2pipe.manage_composable.http_get')
def test_look_pull_and_put(http_mock, mock_client):
    stat_orig = os.stat
    os.stat = Mock()
    os.stat.return_value = Mock(st_size=1234)
    try:
        f_name = 'test_f_name.fits'
        url = 'https://localhost/{}'.format(f_name)
        test_config = mc.Config()
        test_config.observe_execution = True
        test_metrics = mc.Metrics(test_config)
        assert len(test_metrics.history) == 0, 'initial history conditions'
        assert len(test_metrics.failures) == 0, 'initial failure conditions'
        mc.look_pull_and_put(f_name, TEST_DATA_DIR, url, 'TEST', 'default',
                             'application/fits', mock_client, 'md5:01234',
                             test_metrics)
        mock_client.put_file.assert_called_with(
            'TEST', f_name, archive_stream='default', md5_check=True,
            mime_encoding=None,
            mime_type='application/fits'), 'mock not called'
        http_mock.assert_called_with(
            url, os.path.join(TEST_DATA_DIR, f_name)), 'http mock not called'
        assert len(test_metrics.history) == 1, 'history conditions'
        assert len(test_metrics.failures) == 0, 'failure conditions'
    finally:
        os.stat = stat_orig


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('caom2repo.core.CAOM2RepoClient')
def test_repo_create(mock_client):
    test_obs = mc.read_obs_from_file(TEST_OBS_FILE)
    test_config = mc.Config()
    test_config.observe_execution = True
    test_metrics = mc.Metrics(test_config)
    assert len(test_metrics.history) == 0, 'initial history conditions'
    assert len(test_metrics.failures) == 0, 'initial failure conditions'

    mc.repo_create(mock_client, test_obs, test_metrics)

    mock_client.create.assert_called_with(test_obs), 'mock not called'
    assert len(test_metrics.history) == 1, 'history conditions'
    assert len(test_metrics.failures) == 0, 'failure conditions'
    assert 'caom2' in test_metrics.history, 'history'
    assert 'create' in test_metrics.history['caom2'], 'create'
    assert 'test_obs_id' in test_metrics.history['caom2']['create'], 'obs id'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('caom2repo.core.CAOM2RepoClient')
def test_repo_get(mock_client):
    test_config = mc.Config()
    test_config.observe_execution = True
    test_metrics = mc.Metrics(test_config)
    assert len(test_metrics.history) == 0, 'initial history conditions'
    assert len(test_metrics.failures) == 0, 'initial failure conditions'

    mc.repo_get(mock_client, 'collection', 'test_obs_id', test_metrics)

    mock_client.read.assert_called_with('collection', 'test_obs_id'), \
        'mock not called'
    assert len(test_metrics.history) == 1, 'history conditions'
    assert len(test_metrics.failures) == 0, 'failure conditions'
    assert 'caom2' in test_metrics.history, 'history'
    assert 'read' in test_metrics.history['caom2'], 'create'
    assert 'test_obs_id' in test_metrics.history['caom2']['read'], 'obs id'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('caom2repo.core.CAOM2RepoClient')
def test_repo_update(mock_client):
    test_obs = mc.read_obs_from_file(TEST_OBS_FILE)
    test_config = mc.Config()
    test_config.observe_execution = True
    test_metrics = mc.Metrics(test_config)
    assert len(test_metrics.history) == 0, 'initial history conditions'
    assert len(test_metrics.failures) == 0, 'initial failure conditions'

    mc.repo_update(mock_client, test_obs, test_metrics)

    mock_client.update.assert_called_with(test_obs), 'mock not called'
    assert len(test_metrics.history) == 1, 'history conditions'
    assert len(test_metrics.failures) == 0, 'failure conditions'
    assert 'caom2' in test_metrics.history, 'history'
    assert 'update' in test_metrics.history['caom2'], 'update'
    assert 'test_obs_id' in test_metrics.history['caom2']['update'], 'obs id'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('caom2repo.core.CAOM2RepoClient')
def test_repo_delete(mock_client):
    test_config = mc.Config()
    test_config.observe_execution = True
    test_metrics = mc.Metrics(test_config)
    assert len(test_metrics.history) == 0, 'initial history conditions'
    assert len(test_metrics.failures) == 0, 'initial failure conditions'

    mc.repo_delete(mock_client, 'coll', 'test_id', test_metrics)

    mock_client.delete.assert_called_with('coll', 'test_id'), 'mock not called'
    assert len(test_metrics.history) == 1, 'history conditions'
    assert len(test_metrics.failures) == 0, 'failure conditions'
    assert 'caom2' in test_metrics.history, 'history'
    assert 'delete' in test_metrics.history['caom2'], 'delete'
    assert 'test_id' in test_metrics.history['caom2']['delete'], 'obs id'
