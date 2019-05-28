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

import distutils.sysconfig
import os
import pytest
import sys

import six

from mock import Mock, patch

from astropy.io import fits

from caom2 import SimpleObservation, Algorithm
from caom2repo import CAOM2RepoClient
from cadcdata import CadcDataClient


if six.PY3:
    from caom2pipe import CadcException
    from caom2pipe import execute_composable as ec
    from caom2pipe import manage_composable as mc

PY_VERSION = '3.6'
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_DIR = os.path.join(THIS_DIR, 'data')
TEST_APP = 'collection2caom2'


class MyExitError(Exception):
    pass


if six.PY3:
    class TestVisit:
        @staticmethod
        def visit(observation, **kwargs):
            x = kwargs['working_directory']
            assert x is not None, 'working directory'
            y = kwargs['science_file']
            assert y is not None, 'science file'
            z = kwargs['log_file_directory']
            assert z is not None, 'log file directory'
            assert observation is not None, 'undefined observation'

    class TestStorageName(ec.StorageName):
        def __init__(self, obs_id=None, file_name=None):
            super(TestStorageName, self).__init__(
                'test_obs_id', 'TEST', '*', 'test_file.fits.gz')
            self.url = 'https://test_url/'

        def is_valid(self):
            return True

    class TestChooser(ec.OrganizeChooser):
        def __init(self):
            super(TestChooser, self).__init__()

        def needs_delete(self, observation):
            return True

    def _init_config():
        test_config = mc.Config()
        test_config.working_directory = THIS_DIR
        test_config.collection = 'OMM'
        test_config.netrc_file = os.path.join(TEST_DATA_DIR, 'test_netrc')
        test_config.work_file = 'todo.txt'
        test_config.logging_level = 'DEBUG'
        test_config.log_file_directory = TEST_DATA_DIR
        test_config.failure_fqn = '{}/fail.txt'.format(TEST_DATA_DIR)
        test_config.retry_fqn = '{}/retry.txt'.format(TEST_DATA_DIR)
        test_config.success_fqn = '{}/good.txt'.format(TEST_DATA_DIR)
        test_config.resource_id = 'ivo://cadc.nrc.ca/sc2repo'
        test_config.features.run_in_airflow = False
        test_config.features.use_file_names = False
        test_config.stream = 'TEST'
        return test_config


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_meta_create_client_execute():
    test_config = _init_config()
    test_cred = None
    data_client_mock = Mock()
    data_client_mock.get_file_info.return_value = {'name': 'test_file.fits'}
    exec_cmd_orig = mc.exec_cmd
    mc.exec_cmd = Mock()
    repo_client_mock = Mock()
    mc.read_obs_from_file = Mock()
    mc.read_obs_from_file.return_value = _read_obs(None)

    test_executor = ec.MetaCreateClient(
        test_config, TestStorageName(), TEST_APP, test_cred,
        data_client_mock, repo_client_mock, meta_visitors=None)
    test_source = '{}/{}/{}.py'.format(distutils.sysconfig.get_python_lib(),
                                       TEST_APP, TEST_APP)
    try:
        test_executor.execute(None)
        assert mc.exec_cmd.called
        mc.exec_cmd.assert_called_with(
            '{} --debug None --observation OMM test_obs_id '
            '--out {}/test_obs_id/test_obs_id.fits.xml --plugin {} '
            '--module {} --lineage '
            'test_obs_id/ad:TEST/test_obs_id.fits.gz'.format(
                TEST_APP, THIS_DIR, test_source, test_source))
        assert repo_client_mock.create.is_called, 'create call missed'
        assert test_executor.url == 'https://test_url/', 'url'
    finally:
        mc.exec_cmd = exec_cmd_orig


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_meta_update_client_execute():
    test_config = _init_config()
    test_cred = None
    data_client_mock = Mock()
    data_client_mock.get_file_info.return_value = {'name': 'test_file.fits'}
    exec_cmd_orig = mc.exec_cmd
    mc.exec_cmd = Mock()
    repo_client_mock = Mock()
    test_executor = ec.MetaUpdateClient(
        test_config, TestStorageName(), TEST_APP, test_cred,
        data_client_mock, repo_client_mock, _read_obs(None),
        meta_visitors=None)
    test_source = '{}/{}/{}.py'.format(distutils.sysconfig.get_python_lib(),
                                       TEST_APP, TEST_APP)
    try:
        test_executor.execute(None)
        assert mc.exec_cmd.called
        mc.exec_cmd.assert_called_with(
            '{} --debug None --in {}/test_obs_id/test_obs_id.fits.xml '
            '--out {}/test_obs_id/test_obs_id.fits.xml --plugin {} '
            '--module {} --lineage '
            'test_obs_id/ad:TEST/test_obs_id.fits.gz'.format(
                TEST_APP, THIS_DIR, THIS_DIR, test_source, test_source))
        assert repo_client_mock.update.is_called, 'update call missed'
    finally:
        mc.exec_cmd = exec_cmd_orig


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_meta_delete_create_client_execute():
    test_config = _init_config()
    test_cred = None
    data_client_mock = Mock()
    data_client_mock.get_file_info.return_value = {'name': 'test_file.fits'}
    exec_cmd_orig = mc.exec_cmd
    mc.exec_cmd = Mock()
    repo_client_mock = Mock()
    test_executor = ec.MetaDeleteCreateClient(
        test_config, TestStorageName(), TEST_APP, test_cred,
        data_client_mock, repo_client_mock, _read_obs(None), None)
    test_source = '{}/{}/{}.py'.format(distutils.sysconfig.get_python_lib(),
                                       TEST_APP, TEST_APP)
    try:
        test_executor.execute(None)
        assert mc.exec_cmd.called
        mc.exec_cmd.assert_called_with(
            '{} --debug None --in {}/test_obs_id/test_obs_id.fits.xml '
            '--out {}/test_obs_id/test_obs_id.fits.xml --plugin {} '
            '--module {} --lineage '
            'test_obs_id/ad:TEST/test_obs_id.fits.gz'.format(
                TEST_APP, THIS_DIR, THIS_DIR, test_source, test_source))
        assert repo_client_mock.update.is_called, 'update call missed'
    finally:
        mc.exec_cmd = exec_cmd_orig


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_local_meta_create_client_execute():
    test_config = _init_config()
    test_cred = None
    data_client_mock = Mock()
    data_client_mock.get_file_info.return_value = {'name': 'test_file.fits'}
    exec_cmd_orig = mc.exec_cmd
    mc.exec_cmd = Mock()
    repo_client_mock = Mock()

    test_executor = ec.LocalMetaCreateClient(
        test_config, TestStorageName(), TEST_APP, test_cred,
        data_client_mock, repo_client_mock, meta_visitors=None)
    test_source = '{}/{}/{}.py'.format(distutils.sysconfig.get_python_lib(),
                                       TEST_APP, TEST_APP)
    try:
        test_executor.execute(None)
        assert mc.exec_cmd.called
        mc.exec_cmd.assert_called_with(
            '{} --debug None --observation OMM test_obs_id '
            '--local {}/test_file.fits --out {}/test_obs_id.fits.xml '
            '--plugin {} --module {} '
            '--lineage test_obs_id/ad:TEST/test_obs_id.fits.gz'.format(
                TEST_APP, THIS_DIR, THIS_DIR, test_source, test_source))
        assert repo_client_mock.create.is_called, 'create call missed'
    finally:
        mc.exec_cmd = exec_cmd_orig


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_local_meta_update_client_execute():
    test_config = _init_config()
    test_cred = None
    data_client_mock = Mock()
    data_client_mock.get_file_info.return_value = {'name': 'test_file.fits'}
    exec_cmd_orig = mc.exec_cmd
    mc.exec_cmd = Mock()
    repo_client_mock = Mock()
    test_executor = ec.LocalMetaUpdateClient(
        test_config, TestStorageName(), TEST_APP, test_cred,
        data_client_mock, repo_client_mock, _read_obs(None),
        meta_visitors=None)
    test_source = '{}/{}/{}.py'.format(distutils.sysconfig.get_python_lib(),
                                       TEST_APP, TEST_APP)
    try:
        test_executor.execute(None)
        assert mc.exec_cmd.called
        mc.exec_cmd.assert_called_with(
            '{} --debug None --in {}/test_obs_id.fits.xml '
            '--out {}/test_obs_id.fits.xml --local {}/test_file.fits '
            '--plugin {} --module {} '
            '--lineage test_obs_id/ad:TEST/test_obs_id.fits.gz'.format(
                TEST_APP, THIS_DIR, THIS_DIR, THIS_DIR, test_source,
                test_source))
        assert repo_client_mock.update.is_called, 'update call missed'
    finally:
        mc.exec_cmd = exec_cmd_orig


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_local_meta_delete_create_client_execute():
    test_config = _init_config()
    test_cred = None
    data_client_mock = Mock()
    data_client_mock.get_file_info.return_value = {'name': 'test_file.fits'}
    exec_cmd_orig = mc.exec_cmd
    mc.exec_cmd = Mock()
    repo_client_mock = Mock()
    test_executor = ec.LocalMetaDeleteCreateClient(
        test_config, TestStorageName(), TEST_APP, test_cred,
        data_client_mock, repo_client_mock, meta_visitors=None,
        observation=_read_obs(None))
    test_source = '{}/{}/{}.py'.format(distutils.sysconfig.get_python_lib(),
                                       TEST_APP, TEST_APP)
    try:
        test_executor.execute(None)
        assert mc.exec_cmd.called
        mc.exec_cmd.assert_called_with(
            '{} --debug None --in {}/test_obs_id.fits.xml '
            '--out {}/test_obs_id.fits.xml --local {}/test_file.fits '
            '--plugin {} --module {} '
            '--lineage test_obs_id/ad:TEST/test_obs_id.fits.gz'.format(
                TEST_APP, THIS_DIR, THIS_DIR, THIS_DIR, test_source,
                test_source))
        assert repo_client_mock.update.is_called, 'update call missed'
    finally:
        mc.exec_cmd = exec_cmd_orig


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_client_visit():
    test_config = _init_config()
    test_cred = None
    data_client_mock = Mock()
    repo_client_mock = Mock()

    test_executor = ec.ClientVisit(test_config,
                                   TestStorageName(), test_cred,
                                   data_client_mock,
                                   repo_client_mock,
                                   meta_visitors=None)

    test_executor.execute(None)
    assert repo_client_mock.read.is_called, 'read call missed'
    assert repo_client_mock.update.is_called, 'update call missed'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_data_execute():
    test_obs_id = 'TEST_OBS_ID'
    test_dir = os.path.join(THIS_DIR, test_obs_id)
    test_fits_fqn = os.path.join(test_dir,
                                 TestStorageName().file_name)
    if not os.path.exists(test_dir):
        os.mkdir(test_dir)
    precondition = open(test_fits_fqn, 'w')
    precondition.close()

    test_data_visitors = [TestVisit]
    os_path_exists_orig = os.path.exists
    os.path.exists = Mock(return_value=True)
    os_listdir_orig = os.listdir
    os.listdir = Mock(return_value=[])
    os_rmdir_orig = os.rmdir
    os.rmdir = Mock()
    test_config = _init_config()
    data_client_mock = Mock()
    data_client_mock.get_file_info.return_value = {'name': 'test_file.fits'}
    repo_client_mock = Mock()
    test_cred = None

    try:
        ec.CaomExecute._data_cmd_info = Mock(side_effect=_get_fname)

        # run the test
        test_executor = ec.DataClient(
            test_config, TestStorageName(), TEST_APP, test_cred,
            data_client_mock, repo_client_mock, test_data_visitors,
            mc.TaskType.MODIFY)
        try:
            test_executor.execute(None)
        except CadcException as e:
            assert False, e

        # check that things worked as expected
        assert data_client_mock.get_file_info.is_called, \
            'get_file_info call missed'
        assert data_client_mock.get_file.is_called, 'get_file call missed'
        assert repo_client_mock.read.is_called, 'read call missed'
        assert repo_client_mock.update.is_called, 'update call missed'

    finally:
        os.path.exists = os_path_exists_orig
        os.listdir = os_listdir_orig
        os.rmdir = os_rmdir_orig


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_data_local_execute():
    test_data_visitors = [TestVisit]

    data_client_mock = Mock()
    data_client_mock.get_file_info.return_value = {'name': 'test_file.fits'}
    repo_client_mock = Mock()
    repo_client_mock.read.return_value = _read_obs(None)
    test_cred = None

    test_config = _init_config()
    # run the test
    test_executor = ec.LocalDataClient(
        test_config, TestStorageName(), TEST_APP,
        test_cred, data_client_mock, repo_client_mock, test_data_visitors)
    try:
        test_executor.execute(None)
    except CadcException as e:
        assert False, e

    # check that things worked as expected - no cleanup
    assert data_client_mock.get_file_info.is_called, \
        'get_file_info call missed'
    assert data_client_mock.get_file.is_called, 'get_file call missed'
    assert repo_client_mock.read.is_called, 'read call missed'
    assert repo_client_mock.update.is_called, 'update call missed'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_data_store():
    test_config = _init_config()
    data_client_mock = Mock()
    repo_client_mock = Mock()
    test_executor = ec.StoreClient(
        test_config, TestStorageName(), 'command_name', '', data_client_mock,
        repo_client_mock)
    try:
        test_executor.execute(None)
    except CadcException as e:
        assert False, e

    # check that things worked as expected - no cleanup
    assert data_client_mock.put_file.is_called, 'put_file call missed'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_scrape():
    # clean up from previous tests
    if os.path.exists(TestStorageName().model_file_name):
        os.remove(TestStorageName().model_file_name)
    netrc = os.path.join(TEST_DATA_DIR, 'test_netrc')
    assert os.path.exists(netrc)

    test_config = _init_config()
    test_config.working_directory = TEST_DATA_DIR
    test_config.logging_level = 'INFO'
    exec_cmd_orig = mc.exec_cmd
    mc.exec_cmd = Mock()
    test_source = '{}/command_name/command_name.py'.format(
        distutils.sysconfig.get_python_lib())

    try:
        test_executor = ec.Scrape(
            test_config, TestStorageName(), 'command_name')
        try:
            test_executor.execute(None)
        except CadcException as e:
            assert False, e
        assert mc.exec_cmd.called
        mc.exec_cmd.assert_called_with(
            'command_name --verbose  '
            '--observation OMM test_obs_id --out {}/test_obs_id.fits.xml '
            '--plugin {} '
            '--module {} '
            '--local {}/test_file.fits.gz '
            '--lineage test_obs_id/ad:TEST/test_obs_id.fits.gz'.format(
                TEST_DATA_DIR, test_source, test_source, TEST_DATA_DIR))

    finally:
        mc.exec_cmd = exec_cmd_orig


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_data_scrape_execute():
    test_data_visitors = [TestVisit]
    read_orig = mc.read_obs_from_file
    mc.read_obs_from_file = Mock(side_effect=_read_obs)
    try:

        test_config = _init_config()

        # run the test
        test_executor = ec.DataScrape(
            test_config, TestStorageName(), TEST_APP,
            test_data_visitors)
        try:
            test_executor.execute(None)
        except CadcException as e:
            assert False, e

        assert mc.read_obs_from_file.is_called, 'read obs call missed'

    finally:
        mc.read_obs_from_file = read_orig


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_organize_executes_client():
    test_obs_id = TestStorageName()
    test_config = _init_config()
    test_config.use_local_files = True
    log_file_directory = os.path.join(THIS_DIR, 'logs')
    test_config.log_file_directory = log_file_directory
    success_log_file_name = 'success_log.txt'
    test_config.success_log_file_name = success_log_file_name
    failure_log_file_name = 'failure_log.txt'
    test_config.failure_log_file_name = failure_log_file_name
    test_config.features.use_clients = True
    retry_file_name = 'retries.txt'
    test_config.retry_file_name = retry_file_name
    exec_cmd_orig = mc.exec_cmd_info
    repo_cmd_orig = ec.CaomExecute.repo_cmd_get_client
    CadcDataClient.__init__ = Mock(return_value=None)
    CAOM2RepoClient.__init__ = Mock(return_value=None)

    try:
        ec.CaomExecute.repo_cmd_get_client = Mock(return_value=None)
        mc.exec_cmd_info = \
            Mock(return_value='INFO:cadc-data:info\n'
                              'File C170324_0054_SCI_prev.jpg:\n'
                              '    archive: OMM\n'
                              '   encoding: None\n'
                              '    lastmod: Mon, 25 Jun 2018 16:52:07 GMT\n'
                              '     md5sum: f37d21c53055498d1b5cb7753e1c6d6f\n'
                              '       name: C120902_sh2-132_J_old_'
                              'SCIRED.fits.gz\n'
                              '       size: 754408\n'
                              '       type: image/jpeg\n'
                              '    umd5sum: 704b494a972eed30b18b817e243ced7d\n'
                              '      usize: 754408\n'.encode('utf-8'))

        test_config.task_types = [mc.TaskType.SCRAPE]
        test_oe = ec.OrganizeExecutes(test_config)
        executors = test_oe.choose(test_obs_id, 'command_name', [], [])
        assert executors is not None
        assert len(executors) == 1
        assert isinstance(executors[0], ec.Scrape)

        test_config.task_types = [mc.TaskType.STORE,
                                  mc.TaskType.INGEST,
                                  mc.TaskType.MODIFY]
        test_oe = ec.OrganizeExecutes(test_config)
        executors = test_oe.choose(test_obs_id, 'command_name', [], [])
        assert executors is not None
        assert len(executors) == 4
        assert isinstance(executors[0], ec.StoreClient), \
            type(executors[0])
        assert isinstance(executors[1],
                          ec.LocalMetaCreateClient)
        assert isinstance(executors[2], ec.LocalDataClient)
        assert isinstance(
            executors[3], ec.CompareChecksumClient)
        assert CadcDataClient.__init__.is_called, 'mock not called'
        assert CAOM2RepoClient.__init__.is_called, 'mock not called'

        test_config.use_local_files = False
        test_config.task_types = [mc.TaskType.INGEST,
                                  mc.TaskType.MODIFY]
        test_oe = ec.OrganizeExecutes(test_config)
        executors = test_oe.choose(test_obs_id, 'command_name', [], [])
        assert executors is not None
        assert len(executors) == 2
        assert isinstance(executors[0], ec.MetaCreateClient)
        assert isinstance(executors[1], ec.DataClient)
        assert CadcDataClient.__init__.is_called, 'mock not called'
        assert CAOM2RepoClient.__init__.is_called, 'mock not called'

        test_config.use_local_files = True
        test_config.task_types = [mc.TaskType.INGEST,
                                  mc.TaskType.MODIFY]
        test_oe = ec.OrganizeExecutes(test_config)
        executors = test_oe.choose(test_obs_id, 'command_name', [], [])
        assert executors is not None
        assert len(executors) == 3
        assert isinstance(
            executors[0], ec.LocalMetaCreateClient)
        assert isinstance(executors[1], ec.LocalDataClient)
        assert isinstance(
            executors[2], ec.CompareChecksumClient)
        assert CadcDataClient.__init__.is_called, 'mock not called'
        assert CAOM2RepoClient.__init__.is_called, 'mock not called'

        test_config.task_types = [mc.TaskType.SCRAPE,
                                  mc.TaskType.MODIFY]
        test_config.use_local_files = True
        test_oe = ec.OrganizeExecutes(test_config)
        executors = test_oe.choose(test_obs_id, 'command_name', [], [])
        assert executors is not None
        assert len(executors) == 2
        assert isinstance(executors[0], ec.Scrape)
        assert isinstance(executors[1], ec.DataScrape)
        assert CadcDataClient.__init__.is_called, 'mock not called'
        assert CAOM2RepoClient.__init__.is_called, 'mock not called'

        test_config.task_types = [mc.TaskType.REMOTE]
        test_config.use_local_files = True
        test_oe = ec.OrganizeExecutes(test_config)
        executors = test_oe.choose(test_obs_id, 'command_name', [], [])
        assert executors is not None
        assert len(executors) == 1
        assert isinstance(executors[0], ec.LocalMetaCreateClientRemoteStorage)
        assert CadcDataClient.__init__.is_called, 'mock not called'
        assert CAOM2RepoClient.__init__.is_called, 'mock not called'

        test_config.task_types = [mc.TaskType.INGEST]
        test_config.use_local_files = False
        test_chooser = TestChooser()
        ec.CaomExecute.repo_cmd_get_client = Mock(return_value=_read_obs(None))
        test_oe = ec.OrganizeExecutes(test_config, test_chooser)
        executors = test_oe.choose(test_obs_id, 'command_name', [], [])
        assert executors is not None
        assert len(executors) == 1
        assert isinstance(executors[0],
                          ec.MetaDeleteCreateClient)
        assert CadcDataClient.__init__.is_called, 'mock not called'
        assert CAOM2RepoClient.__init__.is_called, 'mock not called'

        test_config.task_types = [mc.TaskType.PULL]
        test_config.use_local_files = False
        test_chooser = TestChooser()
        ec.CaomExecute.repo_cmd_get_client = Mock(return_value=_read_obs(None))
        test_oe = ec.OrganizeExecutes(test_config, test_chooser,
                                      '/tmp/todo.txt')
        executors = test_oe.choose(test_obs_id, 'command_name', [], [])
        assert executors is not None
        assert len(executors) == 1
        assert isinstance(executors[0], ec.PullClient)
        assert CadcDataClient.__init__.is_called, 'mock not called'
        assert CAOM2RepoClient.__init__.is_called, 'mock not called'
        assert executors[0].url == 'https://test_url/', 'url'
        assert executors[0].fname == 'test_obs_id.fits', 'file name'
        assert executors[0].stream == 'TEST', 'stream'
        assert executors[0].working_dir == '{}/test_obs_id'.format(THIS_DIR), \
            'working_dir'
        assert executors[0].local_fqn == \
            '{}/test_obs_id/test_obs_id.fits'.format(THIS_DIR), 'local_fqn'
        assert test_oe.success_fqn == \
            '{}/logs/todo_success_log.txt'.format(THIS_DIR), 'wrong success'
        assert test_oe.retry_fqn == \
            '{}/logs/todo_retries.txt'.format(THIS_DIR), 'wrong retry'
        assert test_oe.failure_fqn == \
            '{}/logs/todo_failure_log.txt'.format(THIS_DIR), 'wrong failure'
        assert test_oe.todo_fqn == '/tmp/todo.txt', 'wrong todo'
    finally:
        mc.exec_cmd_orig = exec_cmd_orig
        ec.CaomExecute.repo_cmd_get_client = repo_cmd_orig


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_organize_executes_chooser():
    test_obs_id = TestStorageName()
    test_config = _init_config()
    test_config.use_local_files = True
    log_file_directory = os.path.join(THIS_DIR, 'logs')
    test_config.log_file_directory = log_file_directory
    test_config.features.supports_composite = True
    exec_cmd_orig = mc.exec_cmd_info
    repo_cmd_orig = ec.CaomExecute.repo_cmd_get_client
    CadcDataClient.__init__ = Mock(return_value=None)
    CAOM2RepoClient.__init__ = Mock(return_value=None)

    try:
        ec.CaomExecute.repo_cmd_get_client = Mock(return_value=_read_obs(None))
        mc.exec_cmd_info = \
            Mock(return_value='INFO:cadc-data:info\n'
                              'File C170324_0054_SCI_prev.jpg:\n'
                              '    archive: OMM\n'
                              '   encoding: None\n'
                              '    lastmod: Mon, 25 Jun 2018 16:52:07 GMT\n'
                              '     md5sum: f37d21c53055498d1b5cb7753e1c6d6f\n'
                              '       name: C120902_sh2-132_J_old_'
                              'SCIRED.fits.gz\n'
                              '       size: 754408\n'
                              '       type: image/jpeg\n'
                              '    umd5sum: 704b494a972eed30b18b817e243ced7d\n'
                              '      usize: 754408\n'.encode('utf-8'))

        test_config.task_types = [mc.TaskType.INGEST]
        test_chooser = TestChooser()
        test_oe = ec.OrganizeExecutes(test_config, test_chooser)
        executors = test_oe.choose(test_obs_id, 'command_name', [], [])
        assert executors is not None
        assert len(executors) == 2
        assert isinstance(executors[0],
                          ec.LocalMetaDeleteCreateClient)
        assert executors[0].fname == 'test_obs_id.fits', 'file name'
        assert executors[0].stream == 'TEST', 'stream'
        assert executors[0].working_dir == THIS_DIR, 'working_dir'
        assert isinstance(executors[1],
                          ec.CompareChecksumClient)

        test_config.use_local_files = False
        test_config.task_types = [mc.TaskType.INGEST]
        test_oe = ec.OrganizeExecutes(test_config, test_chooser)
        executors = test_oe.choose(test_obs_id, 'command_name', [], [])
        assert executors is not None
        assert len(executors) == 1
        assert isinstance(executors[0],
                          ec.MetaDeleteCreateClient)
        assert CadcDataClient.__init__.is_called, 'mock not called'
        assert CAOM2RepoClient.__init__.is_called, 'mock not called'
    finally:
        mc.exec_cmd_orig = exec_cmd_orig
        ec.CaomExecute.repo_cmd_get_client = repo_cmd_orig


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_organize_executes_client_existing():
    test_obs_id = TestStorageName()
    test_config = _init_config()
    test_config.features.use_clients = True
    repo_cmd_orig = ec.CaomExecute.repo_cmd_get_client
    CadcDataClient.__init__ = Mock(return_value=None)
    CAOM2RepoClient.__init__ = Mock(return_value=None)
    try:

        ec.CaomExecute.repo_cmd_get_client = Mock(return_value=_read_obs(None))

        test_config.task_types = [mc.TaskType.INGEST]
        test_config.use_local_files = False
        test_oe = ec.OrganizeExecutes(test_config)
        executors = test_oe.choose(test_obs_id, 'command_name', [], [])
        assert executors is not None
        assert len(executors) == 1
        assert isinstance(executors[0], ec.MetaUpdateClient)
        assert CadcDataClient.__init__.is_called, 'mock not called'
        assert CAOM2RepoClient.__init__.is_called, 'mock not called'
    finally:
        ec.CaomExecute.repo_cmd_get_client = repo_cmd_orig


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_organize_executes_client_visit():
    test_obs_id = TestStorageName()
    test_config = _init_config()
    test_config.features.use_clients = True
    test_config.task_types = [mc.TaskType.VISIT]
    test_config.use_local_files = False
    test_oe = ec.OrganizeExecutes(test_config)
    CadcDataClient.__init__ = Mock(return_value=None)
    CAOM2RepoClient.__init__ = Mock(return_value=None)
    executors = test_oe.choose(test_obs_id, 'command_name', [], [])
    assert executors is not None
    assert len(executors) == 1
    assert isinstance(executors[0], ec.ClientVisit)
    assert CadcDataClient.__init__.is_called, 'mock not called'
    assert CAOM2RepoClient.__init__.is_called, 'mock not called'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_checksum_client():
    test_config = _init_config()
    test_executor = ec.CompareChecksumClient(
        test_config, TestStorageName(), 'test2caom2', None, None, None)
    compare_orig = mc.compare_checksum_client

    try:
        mc.compare_checksum_client = Mock()
        test_executor.execute(None)
        assert mc.compare_checksum_client.called
        assert test_executor.fname == 'test_file.fits.gz', 'fname'
        assert test_executor.working_dir == THIS_DIR, 'working dir'
        assert test_executor.model_fqn == os.path.join(
            THIS_DIR, 'test_obs_id.fits.xml'), 'model fqn'
        assert test_executor.url == 'https://test_url/', 'url'
    finally:
        mc.compare_checksum_client = compare_orig


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_capture_failure():
    test_obs_id = 'test_obs_id'
    test_config = _init_config()
    log_file_directory = os.path.join(THIS_DIR, 'logs')
    test_config.log_to_file = True
    test_config.log_file_directory = log_file_directory
    success_log_file_name = 'success_log.txt'
    test_config.success_log_file_name = success_log_file_name
    failure_log_file_name = 'failure_log.txt'
    test_config.failure_log_file_name = failure_log_file_name
    retry_file_name = 'retries.txt'
    test_config.retry_file_name = retry_file_name

    if not os.path.exists(log_file_directory):
        os.mkdir(log_file_directory)
    if os.path.exists(test_config.success_fqn):
        os.remove(test_config.success_fqn)
    if os.path.exists(test_config.failure_fqn):
        os.remove(test_config.failure_fqn)
    if os.path.exists(test_config.retry_fqn):
        os.remove(test_config.retry_fqn)

    # clean up from last execution

    test_oe = ec.OrganizeExecutes(test_config)
    test_oe.capture_failure(test_obs_id, None, 'exception text')
    test_oe.capture_success(test_obs_id, 'C121212_01234_CAL.fits.gz')

    assert os.path.exists(test_config.success_fqn)
    assert os.path.exists(test_config.failure_fqn)
    assert os.path.exists(test_config.retry_fqn)

    success_content = open(test_config.success_fqn).read()
    assert success_content.endswith('test_obs_id C121212_01234_CAL.fits.gz\n')
    retry_content = open(test_config.retry_fqn).read()
    assert retry_content == 'test_obs_id\n'
    failure_content = open(test_config.failure_fqn).read()
    assert failure_content.endswith('test_obs_id None exception text\n')


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('sys.exit', Mock(side_effect=MyExitError))
def test_run_by_file():
    try:
        os.getcwd = Mock(return_value=TEST_DATA_DIR)
        todo_file = os.path.join(os.getcwd(), 'todo.txt')
        f = open(todo_file, 'w')
        f.write('')
        f.close()
        ec.run_by_file(ec.StorageName, TEST_APP, 'collection',
                       proxy=None, meta_visitors=None, data_visitors=None)
    except mc.CadcException as e:
        assert False, 'but the work list is empty {}'.format(e)


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('sys.exit', Mock(side_effect=MyExitError))
def test_run_by_file_expects_retry():
    retry_dir = '{}_0'.format(TEST_DATA_DIR)
    if os.path.exists(retry_dir):
        f_log = '{}/failure_log.txt'.format(retry_dir)
        if os.path.exists(f_log):
            os.remove(f_log)
        r_txt = '{}/retries.txt'.format(retry_dir)
        if os.path.exists(r_txt):
            os.remove(r_txt)
        s_log = '{}/success_log.txt'.format(retry_dir)
        if os.path.exists(s_log):
            os.remove(s_log)
        t_log = '{}/test_obs_id.log'.format(retry_dir)
        if os.path.exists(t_log):
            os.remove(t_log)
        os.rmdir(retry_dir)

    test_config = _init_config()
    test_config.log_to_file = True
    test_config.features.expects_retry = True
    test_config.retry_failures = True
    test_config.retry_count = 1
    test_config.retry_file_name = 'retries.txt'
    test_config.success_log_file_name = 'success_log.txt'
    test_config.failure_log_file_name = 'failure_log.txt'
    test_retry_count = 0
    test_config.task_types = []
    assert test_config.log_file_directory == TEST_DATA_DIR
    assert test_config.work_file == 'todo.txt'

    assert test_config.need_to_retry(), 'should require retries'

    test_config.update_for_retry(test_retry_count)
    assert test_config.log_file_directory == '{}_{}'.format(TEST_DATA_DIR,
                                                            test_retry_count)
    assert test_config.work_file == 'retries.txt'
    assert test_config.work_fqn == os.path.join(TEST_DATA_DIR, 'retries.txt')
    try:
        ec._run_by_file(test_config, TestStorageName, TEST_APP,
                        proxy=None, meta_visitors=[], data_visitors=[])
    except mc.CadcException as e:
        assert False, 'but the work list is empty {}'.format(e)

    if TEST_DATA_DIR.startswith('/usr/src/app'):
        # these checks fail on travis ....
        assert os.path.exists('{}_0'.format(TEST_DATA_DIR))
        assert os.path.exists(test_config.success_fqn)
        assert os.path.exists(test_config.failure_fqn)
        assert os.path.exists(test_config.retry_fqn)


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_do_one():
    test_config = _init_config()
    test_config.task_types = []
    test_organizer = ec.OrganizeExecutes(test_config)
    # no client
    test_result = ec._do_one(config=test_config, organizer=test_organizer,
                             storage_name=TestStorageName(),
                             command_name='test2caom2',
                             meta_visitors=[], data_visitors=[])
    assert test_result is not None
    assert test_result == -1

    # client
    test_config.features.use_clients = True
    test_result = ec._do_one(config=test_config, organizer=test_organizer,
                             storage_name=TestStorageName(),
                             command_name='test2caom2',
                             meta_visitors=[], data_visitors=[])
    assert test_result is not None
    assert test_result == -1


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_storage_name():
    sn = ec.StorageName(obs_id='test_obs_id', collection='TEST',
                        collection_pattern='T[\\w+-]+')
    assert sn.file_uri == 'ad:TEST/test_obs_id.fits.gz'
    assert sn.file_name == 'test_obs_id.fits'
    assert sn.compressed_file_name == 'test_obs_id.fits.gz'
    assert sn.model_file_name == 'test_obs_id.fits.xml'
    assert sn.prev == 'test_obs_id_prev.jpg'
    assert sn.thumb == 'test_obs_id_prev_256.jpg'
    assert sn.prev_uri == 'ad:TEST/test_obs_id_prev.jpg'
    assert sn.thumb_uri == 'ad:TEST/test_obs_id_prev_256.jpg'
    assert sn.obs_id == 'test_obs_id'
    assert sn.log_file == 'test_obs_id.log'
    assert sn.product_id == 'test_obs_id'
    assert sn.fname_on_disk is None
    assert not sn.is_valid()
    sn = ec.StorageName(obs_id='Test_obs_id', collection='TEST',
                        collection_pattern='T[\\w+-]+')
    assert sn.is_valid()
    x = ec.StorageName.remove_extensions('test_obs_id.fits.header.gz')
    assert x == 'test_obs_id'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_caom_name():
    cn = ec.CaomName(uri='ad:TEST/test_obs_id.fits.gz')
    assert cn.file_id == 'test_obs_id'
    assert cn.file_name == 'test_obs_id.fits.gz'
    assert cn.uncomp_file_name == 'test_obs_id.fits'
    assert ec.CaomName.make_obs_uri_from_obs_id('TEST', 'test_obs_id') == \
        'caom:TEST/test_obs_id'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('sys.exit', Mock(side_effect=MyExitError))
def test_local_meta_create_client_remote_storage_execute():
    os_path_exists_orig = os.path.exists
    os.path.exists = Mock(return_value=True)
    os_listdir_orig = os.listdir
    os.listdir = Mock(return_value=[])
    os_rmdir_orig = os.rmdir
    os.rmdir = Mock()
    test_config = _init_config()
    data_client_mock = Mock()
    data_client_mock.get_file_info.return_value = {'name': 'test_file.fits'}
    repo_client_mock = Mock()
    test_cred = None
    exec_cmd_orig = mc.exec_cmd
    mc.exec_cmd = Mock()
    mc.read_obs_from_file = Mock()
    mc.read_obs_from_file.return_value = _read_obs(None)
    test_source = '{}/{}/{}.py'.format(distutils.sysconfig.get_python_lib(),
                                       TEST_APP, TEST_APP)
    test_local = '{}/test_obs_id.fits'.format(THIS_DIR)

    try:
        ec.CaomExecute._data_cmd_info = Mock(side_effect=_get_fname)

        # run the test
        test_executor = ec.LocalMetaCreateClientRemoteStorage(
            test_config, TestStorageName(), TEST_APP, test_cred,
            data_client_mock, repo_client_mock, None)
        try:
            test_executor.execute(None)
        except CadcException as e:
            assert False, e

        # check that things worked as expected
        assert repo_client_mock.create.is_called, 'create call missed'
        assert mc.exec_cmd.called
        mc.exec_cmd.assert_called_with(
            '{} --debug None --observation OMM test_obs_id --local {} '
            '--out {}/test_obs_id.fits.xml --plugin {} '
            '--module {} --lineage '
            'test_obs_id/ad:TEST/test_obs_id.fits.gz'.format(
                TEST_APP, test_local, THIS_DIR, test_source, test_source))
    finally:
        os.path.exists = os_path_exists_orig
        os.listdir = os_listdir_orig
        os.rmdir = os_rmdir_orig
        mc.exec_cmd = exec_cmd_orig


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('sys.exit', Mock(side_effect=MyExitError))
def test_local_meta_update_client_remote_storage_execute():
    os_path_exists_orig = os.path.exists
    os.path.exists = Mock(return_value=True)
    os_listdir_orig = os.listdir
    os.listdir = Mock(return_value=[])
    os_rmdir_orig = os.rmdir
    os.rmdir = Mock()
    test_config = _init_config()
    data_client_mock = Mock()
    data_client_mock.get_file_info.return_value = {'name': 'test_file.fits'}
    repo_client_mock = Mock()
    test_cred = None
    exec_cmd_orig = mc.exec_cmd
    mc.exec_cmd = Mock()
    mc.read_obs_from_file = Mock()
    mc.read_obs_from_file.return_value = _read_obs(None)
    test_source = '{}/{}/{}.py'.format(distutils.sysconfig.get_python_lib(),
                                       TEST_APP, TEST_APP)
    test_local = '{}/test_obs_id.fits'.format(THIS_DIR)

    try:
        ec.CaomExecute._data_cmd_info = Mock(side_effect=_get_fname)
        # run the test
        test_executor = ec.LocalMetaUpdateClientRemoteStorage(
            test_config, TestStorageName(), TEST_APP, test_cred,
            data_client_mock, repo_client_mock, _read_obs(None), None)
        try:
            test_executor.execute(None)
        except CadcException as e:
            assert False, e

        # check that things worked as expected
        assert repo_client_mock.read.is_called, 'read call missed'
        assert repo_client_mock.update.is_called, 'update call missed'
        assert mc.exec_cmd.called
        mc.exec_cmd.assert_called_with(
            '{} --debug None --in {}/test_obs_id.fits.xml '
            '--out {}/test_obs_id.fits.xml --local {} --plugin {} '
            '--module {} --lineage '
            'test_obs_id/ad:TEST/test_obs_id.fits.gz'.format(
                TEST_APP, THIS_DIR, THIS_DIR, test_local, test_source,
                test_source))
    finally:
        os.path.exists = os_path_exists_orig
        os.listdir = os_listdir_orig
        os.rmdir = os_rmdir_orig
        mc.exec_cmd = exec_cmd_orig


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_omm_name_dots():
    TEST_NAME = 'C121121_J024345.57-021326.4_K_SCIRED'
    TEST_URI = 'ad:OMM/{}.fits.gz'.format(TEST_NAME)
    test_file_id = ec.CaomName(TEST_URI).file_id
    assert TEST_NAME == test_file_id, 'dots messing with things'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('sys.exit', Mock(side_effect=MyExitError))
def test_pull_client():
    # Response mock
    class Object(object):
        pass

        def raise_for_status(self):
            pass

        def iter_content(self, chunk_size):
            return ['aaa'.encode(), 'bbb'.encode()]

        def __enter__(self):
            return self

        def __exit__(self, a, b, c):
            return None

    test_config = _init_config()
    data_client_mock = Mock()
    repo_client_mock = Mock()
    test_sn = TestStorageName()
    test_sn.url = 'file://{}/{}'.format(TEST_DATA_DIR, 'C111107_0694_SCI.fits')
    test_sn.fname_on_disk = '{}/{}'.format(TEST_DATA_DIR, 'x.fits')
    ec.CaomExecute._cleanup = Mock()
    with patch('requests.get') as get_mock:
        get_mock.return_value = Object()
        test_executor = ec.PullClient(test_config, test_sn, TEST_APP, None,
                                      data_client_mock, repo_client_mock)
        with pytest.raises(OSError):
            test_executor.execute(None)
            assert data_client_mock.put_file.is_called, 'call missed'
            assert ec.CaomExecute._cleanup.is_called, 'cleanup call missed'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('sys.exit', Mock(side_effect=MyExitError))
def test_choose_exceptions():
    test_config = _init_config()
    test_config.init_local_files = False
    test_config.task_types = [mc.TaskType.SCRAPE]
    with pytest.raises(mc.CadcException):
        test_organizer = ec.OrganizeExecutes(test_config)
        test_organizer.choose(TestStorageName(), 'command name', [], [])

    test_config.task_types = [mc.TaskType.STORE]
    with pytest.raises(mc.CadcException):
        test_organizer = ec.OrganizeExecutes(test_config)
        test_organizer.choose(TestStorageName(), 'command name', [], [])


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('sys.exit', Mock(side_effect=MyExitError))
def test_storage_name_failure():
    class TestStorageNameFails(TestStorageName):

        def __init__(self):
            super(TestStorageNameFails, self).__init__()

        def is_valid(self):
            return False
    test_config = _init_config()
    test_config.log_to_file = True
    good_start = os.path.getmtime(test_config.success_fqn)
    fail_start = os.path.getmtime(test_config.failure_fqn)
    retry_start = os.path.getmtime(test_config.retry_fqn)
    test_organizer = ec.OrganizeExecutes(test_config)
    test_organizer.choose(TestStorageNameFails(), 'command name', [], [])
    good_end = os.path.getmtime(test_config.success_fqn)
    fail_end = os.path.getmtime(test_config.failure_fqn)
    retry_end = os.path.getmtime(test_config.retry_fqn)
    assert good_end > good_start, 'good'
    assert retry_end > retry_start, 'retry'
    assert fail_end > fail_start, 'failure'


def _communicate():
    return ['return status', None]


def _get_headers(uri, subject):
    x = """SIMPLE  =                    T / Written by IDL:  Fri Oct  6 01:48:35 2017
BITPIX  =                  -32 / Bits per pixel
NAXIS   =                    2 / Number of dimensions
NAXIS1  =                 2048 /
NAXIS2  =                 2048 /
DATATYPE= 'REDUC   '           /Data type, SCIENCE/CALIB/REJECT/FOCUS/TEST
END
"""
    delim = '\nEND'
    extensions = \
        [e + delim for e in x.split(delim) if e.strip()]
    headers = [fits.Header.fromstring(e, sep='\n') for e in extensions]
    return headers


def _get_test_metadata(subject, path):
    return {'size': 37,
            'md5sum': 'e330482de75d5c4c88ce6f6ef99035ea',
            'type': 'applicaton/octect-stream'}


def _get_test_file_meta(path):
    return _get_test_metadata(None, None)


def _read_obs(arg1):
    return SimpleObservation(collection='test_collection',
                             observation_id='test_obs_id',
                             algorithm=Algorithm(str('exposure')))


def _get_file_headers(fname):
    return _get_headers(None, None)


def _get_fname():
    return 'TBD'


def _test_map_todo():
    """For a mock."""
    return ''


def _get_file_info():
    return {'fname': 'test_file.fits'}
