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
import sys

from mock import patch

from caom2pipe import execute_composable as ec

PY_VERSION = '3.6'
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_DIR = os.path.join(THIS_DIR, 'data')
COMMAND_NAME = 'test_command'
REJECTED_FILE = os.path.join(TEST_DATA_DIR, 'rejected.yml')
TEST_URL = 'http://localhost/vlass.fits'
TEST_OBS_ID = 'TEST_OBS_ID'
TEST_ENTRY = 'TEST_ENTRY'


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('caom2pipe.execute_composable._do_one')
def test_run_from_state(do_mock, test_config):
    cleanup_log_txt(test_config)

    test_work = [TEST_URL]

    test_config.features.expects_retry = True
    test_config.retry_failures = True
    test_config.retry_count = 1
    test_config.log_to_file = True

    do_mock.side_effect = _write_retry
    do_mock.return_value = 0

    test_result = ec.run_from_state(test_config,
                                    sname=ec.StorageName,
                                    command_name=COMMAND_NAME,
                                    meta_visitors=None,
                                    data_visitors=None,
                                    todo=test_work)
    assert test_result is not None, 'expect a result'
    assert test_result == 0, 'wrong result'

    assert do_mock.called, 'do mock not called'
    assert do_mock.call_count == 2, do_mock.call_count
    args, kwargs = do_mock.call_args
    test_storage = args[2]
    assert isinstance(test_storage, ec.StorageName), type(test_storage)
    assert test_storage.obs_id is None, 'wrong obs id'
    assert test_storage.url == TEST_ENTRY, test_storage.url

    assert os.path.exists(REJECTED_FILE)


@pytest.mark.skipif(not sys.version.startswith(PY_VERSION),
                    reason='support one python version')
@patch('caom2pipe.execute_composable._do_one')
def test_run_by_file_prime(do_mock, test_config):
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
    do_mock.return_value = 0

    test_result = ec.run_by_file_prime(test_config,
                                       storage_name=ec.StorageName,
                                       command_name=COMMAND_NAME,
                                       meta_visitors=None,
                                       data_visitors=None,
                                       chooser=None)
    assert test_result is not None, 'expect a result'
    assert test_result == 0, 'wrong result'

    assert do_mock.called, 'do mock not called'
    assert do_mock.call_count == 2, do_mock.call_count
    args, kwargs = do_mock.call_args
    test_storage = args[2]
    assert isinstance(test_storage, ec.StorageName), type(test_storage)
    assert test_storage.obs_id == TEST_ENTRY, 'wrong obs id'
    assert test_storage.url is None, test_storage.url

    assert os.path.exists(REJECTED_FILE)


def cleanup_log_txt(config):
    for fqn in [config.success_fqn, config.failure_fqn, config.retry_fqn,
                config.rejected_fqn]:
        if os.path.exists(fqn):
            os.unlink(fqn)
    retry_dir = '{}_0'.format(TEST_DATA_DIR)
    if os.path.exists(retry_dir):
        for ii in os.listdir(retry_dir):
            os.unlink('{}/{}'.format(retry_dir, ii))
        os.rmdir(retry_dir)


def _write_retry(config, organizer, storage_name, command_name,
                 meta_visitors, data_visitors):
    with open(config.retry_fqn, 'w') as f:
        f.write('{}\n'.format(TEST_ENTRY))

