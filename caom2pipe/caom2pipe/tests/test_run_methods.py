# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2019.                            (c) 2019.
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
import six
import shutil
import sys

from datetime import datetime
from mock import patch, Mock

if six.PY3:
    from caom2pipe import execute_composable as ec
    from caom2pipe import manage_composable as mc

PY_VERSION = '3.6'
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_DIR = os.path.join(THIS_DIR, 'data')
COMMAND_NAME = 'test_command'
PROGRESS_FILE = os.path.join(TEST_DATA_DIR, 'progress.txt')
REJECTED_FILE = os.path.join(TEST_DATA_DIR, 'rejected.yml')
STATE_FILE = os.path.join(TEST_DATA_DIR, 'test_state.yml')
TEST_URL = 'http://localhost/vlass.fits'
TEST_OBS_ID = 'TEST_OBS_ID'
TEST_ENTRY = 'TEST_ENTRY'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('caom2pipe.execute_composable._do_one')
@patch('caom2pipe.manage_composable.Work')
def test_run_from_state(work_mock, do_mock, test_config):
    cleanup_log_txt(test_config)

    work_mock.todo.return_value = [TEST_ENTRY]
    work_mock.max_ts_s.return_value = datetime.utcnow().timestamp() + 6

    test_config.features.expects_retry = False
    test_config.progress_fqn = PROGRESS_FILE

    test_config.state_fqn = STATE_FILE
    test_config.interval = 5
    test_state = mc.State(test_config.state_fqn)
    test_state.save_state('gemini_timestamp', datetime.utcnow())

    do_mock.return_value = 0

    sys.argv = ['test_command']

    test_result = ec.run_from_state(
        test_config,
        sname=ec.StorageName,
        command_name=COMMAND_NAME,
        meta_visitors=None,
        data_visitors=None,
        bookmark_name='gemini_timestamp',
        work=work_mock)
    assert test_result is not None, 'expect a result'
    assert test_result == 0, 'wrong result'

    assert do_mock.called, 'do mock not called'
    assert do_mock.call_count == 1, do_mock.call_count
    args, kwargs = do_mock.call_args
    test_storage = args[2]
    assert isinstance(test_storage, ec.StorageName), type(test_storage)
    assert test_storage.obs_id is None, 'wrong obs id'
    assert test_storage.url == TEST_ENTRY, test_storage.url

    assert work_mock.todo.called, 'work todo mock not called'
    assert work_mock.initialize.called, 'work initialize mock not called'

    assert os.path.exists(REJECTED_FILE)


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('caom2pipe.execute_composable._do_one')
def test_run_single_from_state(do_mock, test_config):
    cleanup_log_txt(test_config)

    test_config.features.expects_retry = False
    test_config.progress_fqn = PROGRESS_FILE

    test_config.state_fqn = STATE_FILE
    test_config.interval = 5
    test_state = mc.State(test_config.state_fqn)
    test_state.save_state('gemini_timestamp', datetime.utcnow())

    do_mock.return_value = 0
    mock_organizer = Mock()

    test_url = 'http://localhost/test_url.fits'
    test_storage_name = ec.StorageName(url=test_url)

    test_result = ec.run_single_from_state(
        mock_organizer,
        test_config,
        test_storage_name,
        COMMAND_NAME,
        meta_visitors=None,
        data_visitors=None)
    assert test_result is not None, 'expect a result'
    assert test_result == 0, 'wrong result'

    assert do_mock.called, 'do mock not called'
    assert do_mock.call_count == 1, do_mock.call_count
    args, kwargs = do_mock.call_args
    test_storage = args[2]
    assert isinstance(test_storage, ec.StorageName), type(test_storage)
    assert test_storage.obs_id is None, 'wrong obs id'
    assert test_storage.url == test_url, test_storage.url

    assert mock_organizer.observable.rejected.persist_state().is_called, \
        'organizer should be called'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('caom2pipe.execute_composable._do_one')
def test_run_single(do_mock, test_config):
    cleanup_log_txt(test_config)

    test_config.features.expects_retry = False
    test_config.progress_fqn = PROGRESS_FILE

    test_config.state_fqn = STATE_FILE
    test_config.interval = 5
    test_state = mc.State(test_config.state_fqn)
    test_state.save_state('gemini_timestamp', datetime.utcnow())

    do_mock.return_value = -1

    test_url = 'http://localhost/test_url.fits'
    test_storage_name = ec.StorageName(url=test_url)

    test_result = ec.run_single(
        test_config,
        test_storage_name,
        COMMAND_NAME,
        meta_visitors=None,
        data_visitors=None)
    assert test_result is not None, 'expect a result'
    assert test_result == -1, 'wrong result'

    assert do_mock.called, 'do mock not called'
    assert do_mock.call_count == 1, do_mock.call_count
    args, kwargs = do_mock.call_args
    test_storage = args[2]
    assert isinstance(test_storage, ec.StorageName), type(test_storage)
    assert test_storage.obs_id is None, 'wrong obs id'
    assert test_storage.url == test_url, test_storage.url


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('caom2pipe.execute_composable._do_one')
def test_run_by_file(do_mock, test_config):
    cleanup_log_txt(test_config)

    test_config.features.expects_retry = True
    test_config.retry_failures = True
    test_config.retry_count = 1
    test_config.log_to_file = True
    test_config.features.use_urls = False
    test_config.features.use_file_names = False

    with open(test_config.work_fqn, 'w') as f:
        f.write('{}\n'.format(TEST_OBS_ID))

    do_mock.side_effect = _write_retry

    sys.argv = ['test_command']

    test_result = ec.run_by_file(test_config,
                                 storage_name=ec.StorageName,
                                 command_name=COMMAND_NAME,
                                 meta_visitors=None,
                                 data_visitors=None,
                                 chooser=None)
    assert test_result is not None, 'expect a result'
    assert test_result == -1, 'wrong result'

    assert do_mock.called, 'do mock not called'
    assert do_mock.call_count == 2, do_mock.call_count
    args, kwargs = do_mock.call_args
    test_storage = args[2]
    assert isinstance(test_storage, ec.StorageName), type(test_storage)
    assert test_storage.obs_id == TEST_ENTRY, 'wrong obs id'
    assert test_storage.url is None, test_storage.url

    assert os.path.exists(REJECTED_FILE)


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('caom2pipe.execute_composable._do_one')
def test_run_by_file_use_local_files(do_mock, test_config):
    cleanup_log_txt(test_config)

    test_config.use_local_files = True
    test_config.features.expects_retry = False
    test_config.log_to_file = False
    test_config.features.use_urls = False
    test_config.features.use_file_names = False
    test_config.working_directory = os.path.join(TEST_DATA_DIR, 'local_files')

    class TestStorageName(ec.StorageName):
        def __init__(self, file_name=None, fname_on_disk=None):
            super(TestStorageName, self).__init__()
            assert file_name in ['test_file.fits', 'test_file.fits.header'], \
                'wrong file name'

    test_result = ec.run_by_file(test_config,
                                 storage_name=TestStorageName,
                                 command_name=COMMAND_NAME,
                                 meta_visitors=None,
                                 data_visitors=None,
                                 chooser=None)
    assert test_result is not None, 'expect a result'
    assert test_result == 0, 'wrong result'

    # no local files, should not be called
    assert do_mock.called, 'do mock not called'
    assert do_mock.call_count == 2, do_mock.call_count
    assert os.path.exists(REJECTED_FILE)


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('caom2pipe.execute_composable._do_one')
def test_run_by_file_use_local_files_chooser(do_mock, test_config):
    cleanup_log_txt(test_config)

    test_config.use_local_files = True
    test_config.features.expects_retry = False
    test_config.log_to_file = False
    test_config.features.use_urls = False
    test_config.features.use_file_names = False
    test_config.working_directory = os.path.join(TEST_DATA_DIR, 'local_files')

    class TestStorageName(ec.StorageName):
        def __init__(self, file_name=None, fname_on_disk=None):
            super(TestStorageName, self).__init__()
            assert file_name in ['test_file.fits.gz',
                                 'test_file.fits.header'], 'wrong file name'

    class TestChooser(ec.OrganizeChooser):
        def __init__(self):
            super(TestChooser, self).__init__()

        def use_compressed(self):
            return True

    test_result = ec.run_by_file(test_config,
                                 storage_name=TestStorageName,
                                 command_name=COMMAND_NAME,
                                 meta_visitors=None,
                                 data_visitors=None,
                                 chooser=TestChooser())
    assert test_result is not None, 'expect a result'
    assert test_result == 0, 'wrong result'

    # no local files, should not be called
    assert do_mock.called, 'do mock not called'
    assert do_mock.call_count == 2, do_mock.call_count
    assert os.path.exists(REJECTED_FILE)


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_time_box(test_config):
    _write_state()
    test_bookmark = 'gemini_timestamp'
    test_config.state_fqn = STATE_FILE
    test_config.interval = 700

    class MakeWork(mc.Work):

        def __init__(self):
            super(MakeWork, self).__init__(
                mc.make_seconds('24-Jul-2019 09:20'))
            self.todo_call_count = 0
            self.zero_called = False
            self.one_called = False
            self.two_called = False

        def initialize(self):
            pass

        def todo(self, prev_exec_date, exec_date):
            if self.todo_call_count == 0:
                assert prev_exec_date == datetime(2019, 7, 23, 9, 51), \
                    'wrong prev'
                assert exec_date == datetime(2019, 7, 23, 21, 31), 'wrong exec'
                self.zero_called = True
            elif self.todo_call_count == 1:
                assert prev_exec_date == datetime(2019, 7, 23, 21, 31), \
                    'wrong prev'
                assert exec_date == datetime(2019, 7, 24, 9, 11), 'wrong exec'
                self.one_called = True
            elif self.todo_call_count == 2:
                assert prev_exec_date == datetime(2019, 7, 24, 9, 11), \
                    'wrong exec'
                assert exec_date == datetime(2019, 7, 24, 9, 20), 'wrong exec'
                self.two_called = True
            self.todo_call_count += 1
            assert self.todo_call_count <= 4, 'loop is written wrong'
            return []

    test_work = MakeWork()

    test_result = ec.run_from_state(test_config,
                                    sname=ec.StorageName,
                                    command_name=COMMAND_NAME,
                                    meta_visitors=None,
                                    data_visitors=None,
                                    bookmark_name=test_bookmark,
                                    work=test_work)
    assert test_result is not None, 'expect a result'

    test_state = mc.State(test_config.state_fqn)
    assert test_work.zero_called, 'missed zero'
    assert test_work.one_called, 'missed one'
    assert test_work.two_called, 'missed two'
    assert test_state.get_bookmark(test_bookmark) == \
        datetime(2019, 7, 24, 9, 20)
    assert test_work.todo_call_count == 3, 'wrong todo call count'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_time_box_equal(test_config):
    _write_state()
    test_bookmark = 'gemini_timestamp'
    test_config.state_fqn = STATE_FILE
    test_config.interval = 700

    class MakeWork(mc.Work):

        def __init__(self):
            super(MakeWork, self).__init__(
                mc.make_seconds('23-Jul-2019 09:51'))
            self.todo_call_count = 0
            self.zero_called = False
            self.one_called = False
            self.two_called = False

        def initialize(self):
            pass

        def todo(self, prev_exec_date, exec_date):
            self.todo_call_count += 1
            assert self.todo_call_count <= 4, 'loop is written wrong'
            return []

    test_work = MakeWork()

    test_result = ec.run_from_state(test_config,
                                    sname=ec.StorageName,
                                    command_name=COMMAND_NAME,
                                    meta_visitors=None,
                                    data_visitors=None,
                                    bookmark_name=test_bookmark,
                                    work=test_work)
    assert test_result is not None, 'expect a result'
    test_state = mc.State(test_config.state_fqn)
    assert test_state.get_bookmark(test_bookmark) == \
           datetime(2019, 7, 23, 9, 51)
    assert test_work.todo_call_count == 0, 'wrong todo call count'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
def test_time_box_once_through(test_config):
    _write_state()
    test_bookmark = 'gemini_timestamp'
    test_config.state_fqn = STATE_FILE
    test_config.interval = 700

    class MakeWork(mc.Work):

        def __init__(self):
            super(MakeWork, self).__init__(
                mc.make_seconds('23-Jul-2019 12:20'))
            self.todo_call_count = 0
            self.zero_called = False

        def initialize(self):
            pass

        def todo(self, prev_exec_date, exec_date):
            if self.todo_call_count == 0:
                assert prev_exec_date == datetime(2019, 7, 23, 9, 51), \
                    'wrong prev'
                assert exec_date == datetime(2019, 7, 23, 12, 20), 'wrong exec'
                self.zero_called = True
            self.todo_call_count += 1
            assert self.todo_call_count <= 4, 'loop is written wrong'
            return []

    test_work = MakeWork()

    test_result = ec.run_from_state(test_config,
                                    sname=ec.StorageName,
                                    command_name=COMMAND_NAME,
                                    meta_visitors=None,
                                    data_visitors=None,
                                    bookmark_name=test_bookmark,
                                    work=test_work)
    assert test_result is not None, 'expect a result'

    test_state = mc.State(test_config.state_fqn)
    assert test_work.zero_called, 'missed zero'
    assert test_state.get_bookmark(test_bookmark) == \
           datetime(2019, 7, 23, 12, 20)
    assert test_work.todo_call_count == 1, 'wrong todo call count'


def cleanup_log_txt(config):
    for fqn in [config.success_fqn, config.failure_fqn, config.retry_fqn,
                config.rejected_fqn, config.progress_fqn]:
        if os.path.exists(fqn):
            os.unlink(fqn)
    retry_dir = '{}_0'.format(TEST_DATA_DIR)
    if os.path.exists(retry_dir):
        shutil.rmtree(retry_dir)


def _write_retry(config, organizer, storage_name, command_name,
                 meta_visitors, data_visitors):
    with open(config.retry_fqn, 'w') as f:
        f.write('{}\n'.format(TEST_ENTRY))
    return -1


def _write_state():
    if os.path.exists(STATE_FILE):
        os.unlink(STATE_FILE)
    with open(STATE_FILE, 'w') as f:
        f.write(
            'bookmarks:\n  gemini_timestamp:\n    last_record: '
            '23-Jul-2019 09:51\n')
