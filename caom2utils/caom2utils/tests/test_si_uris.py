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

import os
import sys
from cadcdata import FileInfo
from caom2 import obs_reader_writer
from caom2utils import caom2blueprint
from unittest.mock import patch
from . import test_collections as tc


@patch('cadcutils.net.ws.WsCapabilities.get_access_url', autospec=True)
@patch('caom2utils.data_util.StorageInventoryClient')
def test_cadc_uri(si_mock, ws_mock):
    def _get_mock(id_ignore, dest, fhead):
        dest.write(
            b"""SIMPLE  =                    T / Written by IDL:  Fri Oct  6 01:48:35 2017
BITPIX  =                  -32 / Bits per pixel
NAXIS   =                    2 / Number of dimensions
NAXIS1  =                 2048 /
NAXIS2  =                 2048 /
DATATYPE= 'REDUC   '           /Data type, SCIENCE/CALIB/REJECT/FOCUS/TEST
END"""
        )

    def _info_mock(uri):
        return FileInfo(
            id=uri, size=12, file_type='application/fits', md5sum='abc')

    si_mock.return_value.cadcget.side_effect = _get_mock
    si_mock.return_value.cadcinfo.side_effect = _info_mock

    working_dir = os.path.join(tc.TESTDATA_DIR, 'si')
    out_fqn = os.path.join(working_dir, 'test_out.xml')
    bp_fqn = os.path.join(working_dir, 'si.blueprint')

    if os.path.exists(out_fqn):
        os.unlink(out_fqn)

    sys.argv = ('caom2gen --debug -o {} --no_validate '
                '--resource-id ivo://cadc.nrc.ca/test '
                '--observation TEST_COLLECTION TEST_OBS_ID '
                '--lineage test_product_id/cadc:TEST/test_file.fits '
                '--blueprint {}'.format(out_fqn, bp_fqn)).split()
    caom2blueprint.caom2gen()

    assert os.path.exists(out_fqn), 'expect output file'
    obs_reader = obs_reader_writer.ObservationReader()
    obs = obs_reader.read(out_fqn)
    assert obs is not None, 'expect an Observation'
    assert obs.algorithm.name == 'exposure', 'wrong algorithm construction'

    if os.path.exists(out_fqn):
        os.unlink(out_fqn)
