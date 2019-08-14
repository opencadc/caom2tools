import os
import pytest
import six
import sys

if six.PY3:
    from caom2pipe import manage_composable as mc

PY_VERSION = '3.6'
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_DIR = os.path.join(THIS_DIR, 'data')
TEST_APP = 'collection2caom2'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@pytest.fixture(scope='function')
def test_config():
    test_config = mc.Config()
    test_config.working_directory = THIS_DIR
    test_config.collection = 'OMM'
    test_config.netrc_file = os.path.join(TEST_DATA_DIR, 'test_netrc')
    test_config.work_file = 'todo.txt'
    test_config.logging_level = 'DEBUG'
    test_config.log_file_directory = TEST_DATA_DIR
    test_config.failure_fqn = '{}/fail.txt'.format(TEST_DATA_DIR)
    test_config.failure_log_file_name = 'fail.txt'
    test_config.retry_fqn = '{}/retry.txt'.format(TEST_DATA_DIR)
    test_config.retry_file_name = 'retry.txt'
    test_config.success_fqn = '{}/good.txt'.format(TEST_DATA_DIR)
    test_config.success_log_file_name = 'good.txt'
    test_config.rejected_fqn = '{}/rejected.yml'.format(TEST_DATA_DIR)
    test_config.progress_fqn = '{}/progress.txt'.format(TEST_DATA_DIR)
    test_config.resource_id = 'ivo://cadc.nrc.ca/sc2repo'
    test_config.features.run_in_airflow = False
    test_config.features.use_file_names = False
    test_config.stream = 'TEST'
    for f_name in [test_config.failure_fqn, test_config.success_fqn,
                   test_config.retry_fqn]:
        if os.path.exists(f_name):
            os.unlink(f_name)
    return test_config
