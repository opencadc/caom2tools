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

from mock import Mock

from astropy.io import fits

from caom2 import SimpleObservation, Algorithm
from caom2pipe import CadcException
from caom2pipe import execute_composable as ec
from caom2pipe import manage_composable as mc


THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TESTDATA_DIR = os.path.join(THIS_DIR, 'data')

TEST_OBS = SimpleObservation(collection='test_collection',
                             observation_id='test_obs_id',
                             algorithm=Algorithm(str('exposure')))


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

    def is_valid(self):
        return True


def _init_config():
    test_config = mc.Config()
    test_config.working_directory = THIS_DIR
    test_config.collection = 'OMM'
    test_config.netrc_file = os.path.join(TESTDATA_DIR, 'test_netrc')
    test_config.work_file = 'todo.txt'
    test_config.logging_level = 'DEBUG'
    test_config.log_file_directory = TESTDATA_DIR
    test_config.resource_id = 'ivo://cadc.nrc.ca/sc2repo'
    test_config.features.run_in_airflow = False
    test_config.features.use_clients = False
    test_config.features.use_file_names = False
    return test_config


def test_meta_execute():
    test_obs_id = 'test_obs_id'
    test_dir = os.path.join(THIS_DIR, test_obs_id)

    # clean up from previous tests
    if os.path.exists(test_dir):
        for ii in os.listdir(test_dir):
            os.remove(os.path.join(test_dir, ii))
        os.rmdir(test_dir)
    netrc = os.path.join(THIS_DIR, 'test_netrc')
    assert os.path.exists(netrc)

    # mocks for this test
    data_cmd_orig = ec.CaomExecute._data_cmd_info
    ec.CaomExecute._data_cmd_info = Mock(side_effect=_get_fname)
    exec_cmd_orig = mc.exec_cmd

    test_storage_name = ec.StorageName(
        'test_obs_id', 'OMM', collection_pattern=None)
    test_config = _init_config()
    try:
        # run the test
        mc.exec_cmd = Mock()
        test_executor = ec.Collection2CaomMeta(test_config, test_storage_name,
                                               'command_name')
        try:
            test_executor.execute(None)
        except CadcException as e:
            assert False, e

        # check that things worked as expected
        assert mc.exec_cmd.called
        mc.exec_cmd.assert_called_with(
            'caom2-repo create --debug --resource-id ivo://cadc.nrc.ca/sc2repo '
            '--netrc {}/data/test_netrc '
            '{}/test_obs_id/test_obs_id.fits.xml'.format(THIS_DIR, THIS_DIR))
        assert ec.CaomExecute._data_cmd_info.called
    finally:
        ec.CaomExecute._data_cmd_info = data_cmd_orig
        ec.exec_cmd = exec_cmd_orig


def test_meta_local_execute():
    # clean up from previous tests
    if os.path.exists(TestStorageName().get_model_file_name()):
        os.remove(TestStorageName().get_model_file_name())
    netrc = os.path.join(TESTDATA_DIR, 'test_netrc')
    assert os.path.exists(netrc)
    exec_cmd_orig = mc.exec_cmd
    mc.exec_cmd = Mock()

    test_config = _init_config()
    test_config.working_directory = TESTDATA_DIR
    test_config.logging_level = 'INFO'

    # run the test
    try:
        try:
            test_executor = ec.Collection2CaomLocalMeta(
                test_config, TestStorageName(), 'command_name')
            test_executor.execute(None)
        except CadcException as e:
            assert False, e
        assert mc.exec_cmd.called
        mc.exec_cmd.assert_called_with(
            'caom2-repo create --verbose --resource-id '
            'ivo://cadc.nrc.ca/sc2repo --netrc {}/test_netrc '
            '{}/test_obs_id.fits.xml'.format(TESTDATA_DIR, TESTDATA_DIR))
    finally:
        mc.exec_cmd = exec_cmd_orig


def test_data_execute():
    test_obs_id = 'TEST_OBS_ID'
    test_dir = os.path.join(THIS_DIR, test_obs_id)
    test_fits_fqn = os.path.join(test_dir,
                                 TestStorageName().get_file_name())
    if not os.path.exists(test_dir):
        os.mkdir(test_dir)
    precondition = open(test_fits_fqn, 'w')
    precondition.close()

    test_data_visitors = [TestVisit]
    read_orig = mc.read_obs_from_file
    mc.read_obs_from_file = Mock(side_effect=_read_obs)
    write_orig = mc.write_obs_to_file
    mc.write_obs_to_file = Mock()
    data_cmd_orig = ec.CaomExecute._data_cmd_info
    os_path_exists_orig = os.path.exists
    os.path.exists = Mock(return_value=True)
    os_listdir_orig = os.listdir
    os.listdir = Mock(return_value=[])
    os_rmdir_orig = os.rmdir
    os.rmdir = Mock()
    exec_cmd_orig = mc.exec_cmd
    mc.exec_cmd = Mock()
    test_config = _init_config()

    try:
        ec.CaomExecute._data_cmd_info = Mock(side_effect=_get_fname)

        # run the test
        test_executor = ec.Collection2CaomData(
            test_config, TestStorageName(), 'collection2caom2',
            test_data_visitors)
        try:
            test_executor.execute(None)
        except CadcException as e:
            assert False, e

        # check that things worked as expected
        assert mc.exec_cmd.called, 'exec cmd not called'
        assert mc.write_obs_to_file.called, 'write obs to file not called'
        mc.exec_cmd.assert_called_with(
            'caom2-repo update --debug --resource-id ivo://cadc.nrc.ca/sc2repo '
            '--netrc {}/test_netrc {}/test_obs_id/test_obs_id.fits.xml'.format(
                TESTDATA_DIR, THIS_DIR))

    finally:
        mc.read_obs_from_file = read_orig
        mc.write_obs_to_file = write_orig
        ec.CaomExecute._data_cmd_info = data_cmd_orig
        os.path.exists = os_path_exists_orig
        os.listdir = os_listdir_orig
        os.rmdir = os_rmdir_orig
        mc.exec_cmd = exec_cmd_orig


def test_data_local_execute():
    test_data_visitors = [TestVisit]

    read_orig = mc.read_obs_from_file
    mc.read_obs_from_file = Mock(side_effect=_read_obs)
    write_orig = mc.write_obs_to_file
    mc.write_obs_to_file = Mock()
    exec_cmd_orig = mc.exec_cmd
    mc.exec_cmd = Mock()

    try:
        test_config = _init_config()
        # run the test
        test_executor = ec.Collection2CaomLocalData(
            test_config, TestStorageName(), 'collection2caom2',
            test_data_visitors)
        try:
            test_executor.execute(None)
        except CadcException as e:
            assert False, e

        # check that things worked as expected - no cleanup
        assert mc.exec_cmd.called
        mc.exec_cmd.assert_called_with(
            'caom2-repo update --debug --resource-id ivo://cadc.nrc.ca/sc2repo '
            '--netrc {}/test_netrc {}/test_obs_id.fits.xml'.format(
                TESTDATA_DIR, THIS_DIR))
    finally:
        mc.read_obs_from_file = read_orig
        mc.write_obs_to_file = write_orig
        mc.exec_cmd = exec_cmd_orig


def test_data_store():
    test_config = _init_config()
    exec_cmd_orig = mc.exec_cmd
    mc.exec_cmd = Mock()
    try:
        # run the test
        test_executor = ec.Collection2CaomStore(test_config, TestStorageName(),
                                                'command_name')
        try:
            test_executor.execute(None)
        except CadcException as e:
            assert False, e
        assert mc.exec_cmd.called
        mc.exec_cmd.assert_called_with(
            'cadc-data put --debug -c --netrc '
            '{}/test_netrc OMM -s None test_file.fits.gz'.format(TESTDATA_DIR))

    finally:
        mc.exec_cmd = exec_cmd_orig


def test_scrape():
    # clean up from previous tests
    if os.path.exists(TestStorageName().get_model_file_name()):
        os.remove(TestStorageName().get_model_file_name())
    netrc = os.path.join(TESTDATA_DIR, 'test_netrc')
    assert os.path.exists(netrc)

    test_config = _init_config()
    test_config.working_directory = TESTDATA_DIR
    test_config.logging_level = 'INFO'
    exec_cmd_orig = mc.exec_cmd
    mc.exec_cmd = Mock()

    try:
        test_executor = ec.Collection2CaomScrape(
            test_config, TestStorageName(), 'command_name')
        try:
            test_executor.execute(None)
        except CadcException as e:
            assert False, e
        assert mc.exec_cmd.called
        mc.exec_cmd.assert_called_with(
            'command_name --verbose --netrc {}/test_netrc '
            '--observation OMM test_obs_id --out {}/test_obs_id.fits.xml '
            '--plugin '
            '/usr/local/lib/python3.6/site-packages/command_name/'
            'command_name.py '
            '--module '
            '/usr/local/lib/python3.6/site-packages/command_name/'
            'command_name.py '
            '--local {}/test_file.fits.gz '
            '--lineage test_obs_id/ad:TEST/test_obs_id.fits.gz'.format(
                TESTDATA_DIR, TESTDATA_DIR, TESTDATA_DIR))

    finally:
        mc.exec_cmd = exec_cmd_orig


def test_data_scrape_execute():
    test_data_visitors = [TestVisit]
    read_orig = mc.read_obs_from_file
    mc.read_obs_from_file = Mock(side_effect=_read_obs)
    try:

        test_config = _init_config()

        # run the test
        test_executor = ec.Collection2CaomDataScrape(
            test_config, TestStorageName(), 'collection2caom2',
            test_data_visitors)
        try:
            test_executor.execute(None)
        except CadcException as e:
            assert False, e

    finally:
        mc.read_obs_from_file = read_orig


def test_organize_executes():
    test_obs_id = TestStorageName()
    test_config = _init_config()
    test_config.use_local_files = True
    log_file_directory = os.path.join(THIS_DIR, 'logs')
    test_config.log_file_directory = log_file_directory
    success_log_file_name = 'success_log.txt'
    test_config.success_log_file_name = success_log_file_name
    failure_log_file_name = 'failure_log.txt'
    test_config.failure_log_file_name = failure_log_file_name
    retry_file_name = 'retries.txt'
    test_config.retry_file_name = retry_file_name
    exec_cmd_orig = mc.exec_cmd_info

    try:
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
        executors = test_oe.choose(test_obs_id, 'command_name')
        assert executors is not None
        assert len(executors) == 1
        assert isinstance(executors[0], ec.Collection2CaomScrape)

        test_config.task_types = [mc.TaskType.STORE,
                                  mc.TaskType.INGEST,
                                  mc.TaskType.MODIFY]
        test_oe = ec.OrganizeExecutes(test_config)
        executors = test_oe.choose(test_obs_id, 'command_name')
        assert executors is not None
        assert len(executors) == 4
        assert isinstance(executors[0], ec.Collection2CaomStore)
        assert isinstance(executors[1], ec.Collection2CaomLocalMeta)
        assert isinstance(executors[2], ec.Collection2CaomLocalData)
        assert isinstance(executors[3], ec.Collection2CaomCompareChecksum)

        test_config.use_local_files = False
        test_config.task_types = [mc.TaskType.INGEST,
                                  mc.TaskType.MODIFY]
        test_oe = ec.OrganizeExecutes(test_config)
        executors = test_oe.choose(test_obs_id, 'command_name')
        assert executors is not None
        assert len(executors) == 2
        assert isinstance(executors[0], ec.Collection2CaomMeta)
        assert isinstance(executors[1], ec.Collection2CaomData)

        test_config.task_types = [mc.TaskType.SCRAPE,
                                  mc.TaskType.MODIFY]
        test_config.use_local_files = True
        test_oe = ec.OrganizeExecutes(test_config)
        executors = test_oe.choose(test_obs_id, 'command_name')
        assert executors is not None
        assert len(executors) == 2
        assert isinstance(executors[0], ec.Collection2CaomScrape)
        assert isinstance(executors[1], ec.Collection2CaomDataScrape)
    finally:
        mc.exec_cmd_orig = exec_cmd_orig


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

    try:
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
        executors = test_oe._choose_client(test_obs_id, 'command_name')
        assert executors is not None
        assert len(executors) == 1
        assert isinstance(executors[0], ec.Collection2CaomScrape)

        test_config.task_types = [mc.TaskType.STORE,
                                  mc.TaskType.INGEST,
                                  mc.TaskType.MODIFY]
        test_oe = ec.OrganizeExecutes(test_config)
        executors = test_oe.choose(test_obs_id, 'command_name')
        assert executors is not None
        assert len(executors) == 4
        assert isinstance(executors[0], ec.Collection2CaomStoreClient), \
            type(executors[0])
        assert isinstance(executors[1], ec.Collection2CaomLocalMetaClient)
        assert isinstance(executors[2], ec.Collection2CaomLocalDataClient)
        assert isinstance(
            executors[3], ec.Collection2CaomCompareChecksumClient)

        test_config.use_local_files = False
        test_config.task_types = [mc.TaskType.INGEST,
                                  mc.TaskType.MODIFY]
        test_oe = ec.OrganizeExecutes(test_config)
        executors = test_oe.choose(test_obs_id, 'command_name')
        assert executors is not None
        assert len(executors) == 2
        assert isinstance(executors[0], ec.Collection2CaomMetaDeleteClient)
        assert isinstance(executors[1], ec.Collection2CaomDataClient)

        test_config.task_types = [mc.TaskType.SCRAPE,
                                  mc.TaskType.MODIFY]
        test_config.use_local_files = True
        test_oe = ec.OrganizeExecutes(test_config)
        executors = test_oe.choose(test_obs_id, 'command_name')
        assert executors is not None
        assert len(executors) == 2
        assert isinstance(executors[0], ec.Collection2CaomScrape)
        assert isinstance(executors[1], ec.Collection2CaomDataScrape)
    finally:
        mc.exec_cmd_orig = exec_cmd_orig


def test_organize_executes_client_augment():
    test_obs_id = TestStorageName()
    test_config = _init_config()
    test_config.features.use_clients = True
    repo_cmd_orig = ec.CaomExecute.repo_cmd_get_client
    try:

        ec.CaomExecute.repo_cmd_get_client = Mock(return_value=None)

        test_config.task_types = [mc.TaskType.AUGMENT]
        test_config.use_local_files = False
        test_oe = ec.OrganizeExecutes(test_config)
        executors = test_oe._choose_client(test_obs_id, 'command_name')
        assert executors is not None
        assert len(executors) == 1
        assert isinstance(executors[0], ec.Collection2CaomMetaCreateClient)

        ec.CaomExecute.repo_cmd_get_client = Mock(return_value=TEST_OBS)
        test_oe = ec.OrganizeExecutes(test_config)
        executors = test_oe.choose(test_obs_id, 'command_name')
        assert executors is not None
        assert len(executors) == 1
        assert isinstance(executors[0], ec.Collection2CaomMetaUpdateClient)
    finally:
        ec.CaomExecute.repo_cmd_get_client = repo_cmd_orig


def test_organize_executes_client_visit():
    test_obs_id = TestStorageName()
    test_config = _init_config()
    test_config.features.use_clients = True
    test_config.task_types = [mc.TaskType.VISIT]
    test_config.use_local_files = False
    test_oe = ec.OrganizeExecutes(test_config)
    executors = test_oe.choose(test_obs_id, 'command_name')
    assert executors is not None
    assert len(executors) == 1
    assert isinstance(executors[0], ec.Collection2CaomClientVisit)


def test_meta_client():
    test_config = _init_config()
    repo_client_mock = Mock()
    test_executor = ec.Collection2CaomLocalMetaClient(test_config,
                                                      TestStorageName(),
                                                      'test2caom2', None,
                                                      repo_client_mock, None)
    test_executor.execute(None)
    assert repo_client_mock.create.called


def test_checksum_client():
    test_config = _init_config()
    test_executor = ec.Collection2CaomCompareChecksumClient(
        test_config, TestStorageName(), 'test2caom2', None, None)
    compare_orig = mc.compare_checksum_client

    try:
        mc.compare_checksum_client = Mock()
        test_executor.execute(None)
        assert mc.compare_checksum_client.called
        assert test_executor.fname == 'test_file.fits.gz', 'fname'
        assert test_executor.working_dir == THIS_DIR, 'working dir'
        assert test_executor.model_fqn == os.path.join(
            THIS_DIR, 'test_obs_id.fits.xml'), 'model fqn'
    finally:
        mc.compare_checksum_client = compare_orig


def test_data_cmd_info():
    exec_cmd_orig = mc.exec_cmd_info
    try:
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
                              '      usize: 754408\n')
        test_config = _init_config()
        test_executor = ec.Collection2CaomMeta(test_config, TestStorageName(),
                                               'command_name')
        test_executor._find_file_name_storage()
        assert test_executor.fname is not None, test_executor.fname
        assert test_executor.fname == 'C120902_sh2-132_J_old_SCIRED.fits.gz', \
            test_executor.fname
        assert mc.exec_cmd_info.called
        mc.exec_cmd_info.assert_called_with(
            'cadc-data info --debug --netrc-file {}/test_netrc OMM '
            'test_obs_id'.format(TESTDATA_DIR))
    finally:
        mc.exec_cmd_orig = exec_cmd_orig


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


def test_run_by_file():
    try:
        os.getcwd = Mock(return_value=TESTDATA_DIR)
        todo_file = os.path.join(os.getcwd(), 'todo.txt')
        f = open(todo_file, 'w')
        f.write('')
        f.close()
        ec.run_by_file(ec.StorageName, 'collection2caom2', 'collection',
                       _test_map_todo)
    except mc.CadcException as e:
        assert False, 'but the work list is empty {}'.format(e)


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
    return TEST_OBS


def _get_file_headers(fname):
    return _get_headers(None, None)


def _get_fname():
    return 'TBD'


def _test_map_todo():
    """For a mock."""
    return ''


def _get_file_info():
    return {'fname': 'test_file.fits'}
