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

from mock import Mock, patch

from caom2pipe import execute_composable as ec
from caom2pipe import manage_composable as mc


THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TESTDATA_DIR = os.path.join(THIS_DIR, 'data')


def test_map_todo():
    """For a mock."""
    return ''


def test_config_class():
    os.getcwd = Mock(return_value=TESTDATA_DIR)
    test_config = mc.Config()
    test_config.get_executors()
    assert test_config is not None
    assert test_config.work_file == 'todo.txt'


def test_run_by_file():
    try:
        os.getcwd = Mock(return_value=TESTDATA_DIR)
        todo_file = os.path.join(os.getcwd(), 'todo.txt')
        f = open(todo_file, 'w')
        f.write('')
        f.close()
        ec.run_by_file(ec.StorageName, 'collection2caom2', 'collection',
                       test_map_todo)
    except mc.CadcException as e:
        assert False, 'but the work list is empty {}'.format(e)


def test_exec_cmd():
    test_cmd = 'ls'
    mc.exec_cmd(test_cmd)


def test_exec_cmd_redirect():
    fqn = os.path.join(TESTDATA_DIR, 'exec_cmd_redirect.txt')
    if os.path.exists(fqn):
        os.remove(fqn)

    test_cmd = 'ls'
    mc.exec_cmd_redirect(test_cmd, fqn)
    assert os.path.exists(fqn)
    assert os.stat(fqn).st_size > 0


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
