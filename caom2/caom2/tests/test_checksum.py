# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2026.                            (c) 2026.
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
from caom2.shape import Circle, MultiShape, Point, Polygon
from caom2.caom_util import str2ivoa, TypedList, TypedSet
from caom2.checksum import update_checksum, int_32, checksum_diff
import tempfile
from unittest.mock import patch
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
    value = str2ivoa('2012-07-11T13:26:37.123200')
    update_checksum(md5, value, False)
    assert ('9f8af3a440b6e1c8e2a7ea86d90685ac' == md5.hexdigest())
    md5 = hashlib.md5()
    value = str2ivoa('2012-07-11T13:26:37.000000')
    update_checksum(md5, value, False)
    assert ('b35eea8d6e70a117ae7804f4e0f6cf58' == md5.hexdigest())
    md5 = hashlib.md5()
    value = str('ad:file')
    update_checksum(md5, value, False)
    assert ('effad6d4f11ff5a2a8fdd4880b7f2081' == md5.hexdigest())
    md5 = hashlib.md5()
    value = UUID('00000000-0000-0000-9d25-b0383f3182a5')
    update_checksum(md5, value, False)
    assert ('5b71d023d4729575d550536dce8439e6' == md5.hexdigest())


def _assert_stream_collision(left, right):
    """Same UTF-8 byte sequence; only boundaries between elements differ."""
    def stream(seq):
        return b''.join(
            x.encode('utf-8') if isinstance(x, str) else x for x in seq)

    assert stream(left) == stream(right), (
        'test data must use a collision pair (same bytes, different elements)')


def test_list_checksum_should_differ_when_element_boundaries_differ():
    """
    update_checksum feeds each list item with no terminator, so the MD5
    input is the concatenation of encodings only. Different lists whose
    items concatenate to the same bytes must still yield different checksums.
    """
    left = ['abcd', 'efgh']
    right = ['abcdef', 'gh']
    _assert_stream_collision(left, right)

    md5_left = hashlib.md5()
    update_checksum(md5_left, left, 'items')
    md5_right = hashlib.md5()
    update_checksum(md5_right, right, 'items')
    assert md5_left.hexdigest() != md5_right.hexdigest(), (
        'checksum must depend on list structure, not only on concatenated '
        'UTF-8 of elements')


def test_typed_list_checksum_should_differ_when_element_boundaries_differ():
    left = TypedList(str, 'abcd', 'efgh')
    right = TypedList(str, 'abcdef', 'gh')
    _assert_stream_collision(list(left), list(right))

    md5_left = hashlib.md5()
    update_checksum(md5_left, left, 'items')
    md5_right = hashlib.md5()
    update_checksum(md5_right, right, 'items')
    assert md5_left.hexdigest() != md5_right.hexdigest(), (
        'TypedList checksum must not match a different partition of the same '
        'byte sequence')


def test_set_checksum_should_differ_when_element_boundaries_differ():
    left = {'abcd', 'efgh'}
    right = {'abcdef', 'gh'}
    _assert_stream_collision(sorted(left), sorted(right))

    md5_left = hashlib.md5()
    update_checksum(md5_left, left, 'items')
    md5_right = hashlib.md5()
    update_checksum(md5_right, right, 'items')
    assert md5_left.hexdigest() != md5_right.hexdigest(), (
        'set checksum must not match a different partition of the same byte '
        'sequence (after sorted iteration)')


def test_typed_set_checksum_should_differ_when_element_boundaries_differ():
    left = TypedSet(str, 'abcd', 'efgh')
    right = TypedSet(str, 'abcdef', 'gh')
    _assert_stream_collision(sorted(left), sorted(right))

    md5_left = hashlib.md5()
    update_checksum(md5_left, left, 'items')
    md5_right = hashlib.md5()
    update_checksum(md5_right, right, 'items')
    assert md5_left.hexdigest() != md5_right.hexdigest(), (
        'TypedSet checksum must not match a different partition of the same '
        'byte sequence')


def _flatten_multishape_primitives(multishape):
    """Concatenate all shape primitive values, ignoring list boundaries."""
    out = []
    for shape_values in multishape.get_unwrapped_value():
        out.extend(shape_values)
    return out


def test_multishape_checksum_should_differ_when_shape_boundaries_differ():
    """
    MultiShape checksum must depend on shape boundaries (0x00 delimiters),
    not only on the contiguous primitive byte sequence.

    polygon(0,10;20,30) + circle(40,50,r=60)  ->  [[0,10,20,30],[40,50,60]]
    circle(0,10,r=20) + polygon(30,40;50,60)  ->  [[0,10,20],[30,40,50,60]]

    Both flatten to [0, 10, 20, 30, 40, 50, 60] but partition it differently.
    """
    polygon_quad = Polygon([Point(0.0, 10.0), Point(20.0, 30.0)])
    circle_tail = Circle(Point(40.0, 50.0), 60.0)
    circle_head = Circle(Point(0.0, 10.0), 20.0)
    polygon_tail = Polygon([Point(30.0, 40.0), Point(50.0, 60.0)])

    left = MultiShape([polygon_quad, circle_tail])
    right = MultiShape([circle_head, polygon_tail])

    assert left.get_unwrapped_value() == [[0.0, 10.0, 20.0, 30.0], [40.0, 50.0, 60.0]]
    assert right.get_unwrapped_value() == [[0.0, 10.0, 20.0], [30.0, 40.0, 50.0, 60.0]]
    assert _flatten_multishape_primitives(left) == _flatten_multishape_primitives(
        right)

    md5_left = hashlib.md5()
    update_checksum(md5_left, left, 'position.samples')
    md5_right = hashlib.md5()
    update_checksum(md5_right, right, 'position.samples')
    assert md5_left.hexdigest() != md5_right.hexdigest(), (
        'MultiShape checksum must depend on shape boundaries, not only on the '
        'contiguous list of primitive values')


def test_checksum_diff():
    for source_file_path in \
            [os.path.join(THIS_DIR, TEST_DATA, x) for x in ['SampleDerived-CAOM-2.5.xml']]:
        # TODO Maybe ['SampleDerived-CAOM-2.4.xml', 'SampleComposite-CAOM-2.3.xml']]:
        logging.debug(source_file_path)
        output_file = tempfile.NamedTemporaryFile()
        sys.argv = 'caom2_checksum -d -o {} {}'.format(
            output_file.name, source_file_path).split()
        with patch('sys.stdout', new_callable=StringIO) as stdout_mock:
            checksum_diff()
            output = stdout_mock.getvalue()
        assert 'mismatch' not in output, '{} should have correct checksum'.\
            format(source_file_path)
        # assert 'chunk' in output  - removed from sample file and maybe from 2.5?
        # assert 'part' in output  - removed from sample file and maybe from 2.5?
        assert 'artifact' in output
        assert 'plane' in output
        assert 'observation' in output

        # original observation and the one output should be identical
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
