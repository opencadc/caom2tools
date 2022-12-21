# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2022.                            (c) 2022.
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

""" Defines TestPlane class """

import hashlib
import os
import sys
from uuid import UUID
import logging

from builtins import int, str

from caom2 import obs_reader_writer, get_meta_checksum, get_acc_meta_checksum
from caom2 import update_meta_checksum
from caom2.caom_util import str2ivoa
from caom2.checksum import update_checksum, int_32, checksum_diff
import tempfile
from mock import patch
from io import StringIO

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA = 'data'


def test_primitive_checksum():
    md5 = hashlib.md5()
    # tests checksums of various primitives to match those in Java
    value = True
    update_checksum(md5, value, False)
    assert ('55a54008ad1ba589aa210d2629c1df41' == md5.hexdigest())
    md5 = hashlib.md5()
    value = False
    update_checksum(md5, value, False)
    assert ('93b885adfe0da089cdf634904fd59f71' == md5.hexdigest())
    md5 = hashlib.md5()
    value = 'hello'
    update_checksum(md5, value, False)
    assert ('5d41402abc4b2a76b9719d911017c592' == md5.hexdigest())
    md5 = hashlib.md5()
    value = int_32(3)
    update_checksum(md5, value, False)
    assert ('584a15a90f2f959d0703594ad447ae93' == md5.hexdigest())
    md5 = hashlib.md5()
    value = int(12345678910)
    update_checksum(md5, value, False)
    assert ('f61cbb413a37d320af998a215530bc78' == md5.hexdigest())
    # md5 = hashlib.md5()
    # value = common.float_32(1.1)
    # common.get_primitive_to_bytes(md5, value, False)
    # assert ('8ce670eb32869bc6b6109d970711f7c1' == md5.hexdigest())
    md5 = hashlib.md5()
    value = 2.2
    update_checksum(md5, value, False)
    assert ('0fec383169e99d1a6bebd89d1cd8fad9' == md5.hexdigest())
    md5 = hashlib.md5()
    value = str2ivoa('2012-07-11T13:26:37.123')
    update_checksum(md5, value, False)
    assert ('aedbcf5e27a17fc2daa5a0e0d7840009' == md5.hexdigest())
    # ensure that the milliseconds part is not part of checksum
    md5 = hashlib.md5()
    value = str2ivoa('2012-07-11T13:26:37.000')
    update_checksum(md5, value, False)
    assert ('aedbcf5e27a17fc2daa5a0e0d7840009' == md5.hexdigest())
    md5 = hashlib.md5()
    value = str('ad:file')
    update_checksum(md5, value, False)
    assert ('effad6d4f11ff5a2a8fdd4880b7f2081' == md5.hexdigest())
    md5 = hashlib.md5()
    value = UUID('00000000-0000-0000-9d25-b0383f3182a5')
    update_checksum(md5, value, False)
    assert ('5b71d023d4729575d550536dce8439e6' == md5.hexdigest())


def test_compatibility():
    # tests loads a previously generated observation and checks the checksums
    # against the previously calculated (in Java) checksums

    source_file_path = os.path.join(THIS_DIR, TEST_DATA,
                                    'SampleComposite-CAOM-2.3.xml')
    reader = obs_reader_writer.ObservationReader(True)
    with open(source_file_path, 'r'):
        obs = reader.read(source_file_path)

    writer = obs_reader_writer.ObservationWriter(
        True, namespace=obs_reader_writer.CAOM23_NAMESPACE)
    writer.write(obs, '/tmp/test.xml')
    _common_check(obs)

    # check observation
    assert obs.meta_checksum == get_meta_checksum(obs)
    assert obs.acc_meta_checksum == get_acc_meta_checksum(obs)

    # white spaces around strings should not affect the checksum
    obs.algorithm = ' {}\t\n'.format(obs.algorithm.name)
    assert obs.meta_checksum == get_meta_checksum(obs)

    # now change some attributes and see how the checksums start to diverge
    old_val = obs.collection
    obs.collection = 'OTHER'
    _common_check(obs)
    assert obs.meta_checksum != get_meta_checksum(obs)
    assert obs.acc_meta_checksum != get_acc_meta_checksum(obs)
    obs.collection = old_val

    # now change a plane
    aplane = list(obs.planes.values())[0]
    old_val = aplane.product_id
    aplane.product_id = 'TESTPRODID'
    for plane in obs.planes.values():
        for artifact in plane.artifacts.values():
            for part in artifact.parts.values():
                for chunk in part.chunks:
                    assert chunk.meta_checksum == get_meta_checksum(chunk)
                    assert chunk.acc_meta_checksum == get_acc_meta_checksum(
                        chunk)
                assert part.meta_checksum == get_meta_checksum(part)
                assert part.acc_meta_checksum == get_acc_meta_checksum(part)
            assert artifact.meta_checksum == get_meta_checksum(artifact)
            assert artifact.acc_meta_checksum == get_acc_meta_checksum(
                artifact)
        if plane._id == aplane._id:
            assert plane.meta_checksum != get_meta_checksum(plane)
            assert plane.acc_meta_checksum != get_acc_meta_checksum(plane)
        else:
            assert plane.meta_checksum == get_meta_checksum(plane)
            assert plane.acc_meta_checksum == get_acc_meta_checksum(plane)
    assert obs.meta_checksum == get_meta_checksum(obs)
    assert obs.acc_meta_checksum != get_acc_meta_checksum(obs)
    aplane.product_id = old_val

    # change an artifact
    anartifact = list(aplane.artifacts.values())[0]
    old_val = anartifact.content_length
    anartifact.content_length = 3344
    for plane in obs.planes.values():
        for artifact in plane.artifacts.values():
            for part in artifact.parts.values():
                for chunk in part.chunks:
                    assert chunk.meta_checksum == get_meta_checksum(chunk)
                    assert chunk.acc_meta_checksum == get_acc_meta_checksum(
                        chunk)
                assert part.meta_checksum == get_meta_checksum(part)
                assert part.acc_meta_checksum == get_acc_meta_checksum(part)
            if artifact._id == anartifact._id:
                assert artifact.meta_checksum != get_meta_checksum(artifact)
                assert artifact.acc_meta_checksum != get_acc_meta_checksum(
                    artifact)
            else:
                assert artifact.meta_checksum == get_meta_checksum(artifact)
                assert artifact.acc_meta_checksum == get_acc_meta_checksum(
                    artifact)
        assert plane.meta_checksum == get_meta_checksum(plane)
        if plane._id == aplane._id:
            assert plane.acc_meta_checksum != get_acc_meta_checksum(plane)
        else:
            assert plane.acc_meta_checksum == get_acc_meta_checksum(plane)
    assert obs.meta_checksum == get_meta_checksum(obs)
    assert obs.acc_meta_checksum != get_acc_meta_checksum(obs)
    anartifact.content_length = old_val

    apart = list(anartifact.parts.values())[0]
    old_val = apart.name
    apart.name = 'therealpart'
    for plane in obs.planes.values():
        for artifact in plane.artifacts.values():
            for part in artifact.parts.values():
                for chunk in part.chunks:
                    assert chunk.meta_checksum == get_meta_checksum(chunk)
                    assert chunk.acc_meta_checksum == get_acc_meta_checksum(
                        chunk)
                if part._id == apart._id:
                    assert part.meta_checksum != get_meta_checksum(part)
                    assert part.acc_meta_checksum != get_acc_meta_checksum(
                        part)
                else:
                    assert part.meta_checksum == get_meta_checksum(part)
                    assert part.acc_meta_checksum == get_acc_meta_checksum(
                        part)
            assert artifact.meta_checksum == get_meta_checksum(artifact)
            if artifact._id == anartifact._id:
                assert artifact.acc_meta_checksum != get_acc_meta_checksum(
                    artifact)
            else:
                assert artifact.acc_meta_checksum == get_acc_meta_checksum(
                    artifact)
        assert plane.meta_checksum == get_meta_checksum(plane)
        if plane._id == aplane._id:
            assert plane.acc_meta_checksum != get_acc_meta_checksum(plane)
        else:
            assert plane.acc_meta_checksum == get_acc_meta_checksum(plane)
    assert obs.meta_checksum == get_meta_checksum(obs)
    assert obs.acc_meta_checksum != get_acc_meta_checksum(obs)
    apart.name = old_val

    achunk = list(apart.chunks)[0]
    old_val = chunk.naxis
    if old_val == 5:
        achunk.naxis = 4
    else:
        achunk.naxis = old_val + 1
    for plane in obs.planes.values():
        for artifact in plane.artifacts.values():
            for part in artifact.parts.values():
                for chunk in part.chunks:
                    if chunk._id == achunk._id:
                        assert chunk.meta_checksum != get_meta_checksum(chunk)
                        assert chunk.acc_meta_checksum !=\
                            get_acc_meta_checksum(chunk)
                    else:
                        assert chunk.meta_checksum == get_meta_checksum(chunk)
                        assert chunk.acc_meta_checksum ==\
                            get_acc_meta_checksum(chunk)
                assert part.meta_checksum == get_meta_checksum(part)
                if part._id == apart._id:
                    assert part.acc_meta_checksum != get_acc_meta_checksum(
                        part)
                else:
                    assert part.acc_meta_checksum == get_acc_meta_checksum(
                        part)
            assert artifact.meta_checksum == get_meta_checksum(artifact)
            if artifact._id == anartifact._id:
                assert artifact.acc_meta_checksum != get_acc_meta_checksum(
                    artifact)
            else:
                assert artifact.acc_meta_checksum == get_acc_meta_checksum(
                    artifact)
        assert plane.meta_checksum == get_meta_checksum(plane)
        if plane._id == aplane._id:
            assert plane.acc_meta_checksum != get_acc_meta_checksum(plane)
        else:
            assert plane.acc_meta_checksum == get_acc_meta_checksum(plane)
    assert obs.meta_checksum == get_meta_checksum(obs)
    assert obs.acc_meta_checksum != get_acc_meta_checksum(obs)
    achunk.naxis = old_val

    # update the checksums and everything should match again
    update_meta_checksum(obs)
    _common_check(obs)


def test_compatibility_simple_obs():
    # tests loads a previously generated observation and checks the checksums
    # against the previously calculated (in Java) checksums
    logger = logging.getLogger('checksum')
    level = logger.getEffectiveLevel()
    logger.setLevel(logging.DEBUG)
    source_file_path = os.path.join(THIS_DIR, TEST_DATA,
                                    'SampleSimple-CAOM-2.3.xml')
    reader = obs_reader_writer.ObservationReader(True)
    with open(source_file_path, 'r'):
        obs = reader.read(source_file_path)

    writer = obs_reader_writer.ObservationWriter(
        True, namespace=obs_reader_writer.CAOM23_NAMESPACE)
    writer.write(obs, '/tmp/test.xml')
    _common_check(obs)

    # check observation
    assert obs.meta_checksum == get_meta_checksum(obs)
    assert obs.acc_meta_checksum == get_acc_meta_checksum(obs)
    logger.setLevel(level)


def test_round_trip():
    source_file_path = os.path.join(THIS_DIR, TEST_DATA,
                                    'SampleComposite-CAOM-2.3.xml')
    reader = obs_reader_writer.ObservationReader(True)
    with open(source_file_path, 'r'):
        obs = reader.read(source_file_path)

    filename = tempfile.TemporaryFile()
    writer = obs_reader_writer.ObservationWriter(
        True, namespace=obs_reader_writer.CAOM23_NAMESPACE)
    writer.write(obs, filename)

    # go back to the beginning of the file
    filename.seek(0)
    obs = reader.read(filename)
    _common_check(obs)

    # check observation
    assert obs.meta_checksum == get_meta_checksum(obs)
    assert obs.acc_meta_checksum == get_acc_meta_checksum(obs)


def test_checksum_diff():
    for source_file_path in \
            [os.path.join(THIS_DIR, TEST_DATA, x) for x in
             ['SampleDerived-CAOM-2.4.xml', 'SampleComposite-CAOM-2.3.xml']]:
        logging.debug(source_file_path)
        output_file = tempfile.NamedTemporaryFile()
        sys.argv = 'caom2_checksum -d -o {} {}'.format(
            output_file.name, source_file_path).split()
        with patch('sys.stdout', new_callable=StringIO) as stdout_mock:
            checksum_diff()
            output = stdout_mock.getvalue()
        assert 'mismatch' not in output, '{} should have correct checksum'.\
            format(source_file_path)
        assert 'chunk' in output
        assert 'part' in output
        assert 'artifact' in output
        assert 'plane' in output
        assert 'observation' in output

        # original observation and the one outputed should be identical
        reader = obs_reader_writer.ObservationReader()
        expected = reader.read(source_file_path)
        actual = reader.read(output_file.name)
        assert get_acc_meta_checksum(expected) == get_acc_meta_checksum(actual)


def _common_check(obs):
    for plane in obs.planes.values():
        for artifact in plane.artifacts.values():
            for part in artifact.parts.values():
                for chunk in part.chunks:
                    assert chunk.meta_checksum == get_meta_checksum(chunk)
                    assert chunk.acc_meta_checksum == get_acc_meta_checksum(
                        chunk)
                assert part.meta_checksum == get_meta_checksum(part)
                assert part.acc_meta_checksum == get_acc_meta_checksum(
                    part)
            assert artifact.meta_checksum == get_meta_checksum(artifact)
            assert artifact.acc_meta_checksum == get_acc_meta_checksum(
                artifact)
        assert plane.meta_checksum == get_meta_checksum(plane)
        assert plane.acc_meta_checksum == get_acc_meta_checksum(plane)
