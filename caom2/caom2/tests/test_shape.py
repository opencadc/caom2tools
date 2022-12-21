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

import math
import pytest
import unittest

from .. import shape


class TestEnums(unittest.TestCase):
    def test_all(self):
        # test for invalid key
        with self.assertRaises(KeyError):
            shape.SegmentType["foo"]
        with self.assertRaises(Exception):
            shape.SegmentType[None]
        with self.assertRaises(Exception):
            shape.SegmentType[999]
        # test for invalid value
        with self.assertRaises(ValueError):
            shape.SegmentType("foo")
        with self.assertRaises(ValueError):
            shape.SegmentType(None)
        with self.assertRaises(ValueError):
            shape.SegmentType(4)
        self.assertEqual(shape.SegmentType.CLOSE.value, 0)
        self.assertEqual(shape.SegmentType.LINE.value, 1)
        self.assertEqual(shape.SegmentType.MOVE.value, 2)


class TestBox(unittest.TestCase):
    def test_all(self):
        self.assertRaises(TypeError, shape.Box, None, None, None)
        self.assertRaises(TypeError, shape.Box, None, None, 1.0)
        self.assertRaises(TypeError, shape.Box, None, 1.0, None)
        self.assertRaises(TypeError, shape.Box, 1.0, None, None)
        self.assertRaises(TypeError, shape.Box, int(1), "string", int(1))
        self.assertRaises(TypeError, shape.Box, "string", int(1), "string")

        val1 = 1.0
        val2 = 2.0
        width = 3.0
        height = 4.0
        box = shape.Box(shape.Point(val1, val2), width, height)
        self.assertEqual(box.width, 3.0)
        self.assertEqual(box.height, 4.0)
        self.assertEqual(box.center.cval1, 1.0)
        self.assertEqual(box.center.cval2, 2.0)
        area = width * height
        self.assertEqual(box.get_area(), area)
        size = math.sqrt(width * width + height * height)
        self.assertEqual(box.get_size(), size)


class TestCircle(unittest.TestCase):
    def test_all(self):
        self.assertRaises(TypeError, shape.Circle, None, None)
        self.assertRaises(TypeError, shape.Circle, None, 1.0)
        self.assertRaises(TypeError, shape.Circle, 1.0, None)
        self.assertRaises(TypeError, shape.Circle, "string", int(1))
        self.assertRaises(TypeError, shape.Circle, int(1), "string")

        val1 = 1.0
        val2 = 2.0
        radius = 3.0
        circle = shape.Circle(shape.Point(val1, val2), radius)
        self.assertEqual(circle.radius, radius)
        self.assertEqual(circle.center.cval1, val1)
        self.assertEqual(circle.center.cval2, val2)
        area = math.pi * radius * radius
        self.assertEqual(circle.get_area(), area)
        self.assertEqual(circle.get_size(), 2.0 * radius)


class TestInterval(unittest.TestCase):
    def test_all(self):

        lower = 1.0
        upper = 2.0
        lower1 = 1.1
        upper1 = 2.1
        lower2 = 1.2
        upper2 = 2.2
        samples = [shape.SubInterval(lower, lower1),
                   shape.SubInterval(lower2, upper),
                   shape.SubInterval(upper1, upper2)]
        invalid_samples_lower_mismatch = [shape.SubInterval(lower, upper)]
        invalid_samples_upper_mismatch = [shape.SubInterval(lower, upper2)]
        invalid_samples_middle_bounds_overlap = [
            shape.SubInterval(lower, upper), shape.SubInterval(lower1, upper1)]

        self.assertRaises(TypeError, shape.Interval, None, None, None)
        self.assertRaises(TypeError, shape.Interval, None, None, 1.0)
        self.assertRaises(TypeError, shape.Interval, None, 1.0, None)
        self.assertRaises(TypeError, shape.Interval, 1.0, None, None)
        self.assertRaises(TypeError, shape.Interval, None, None, samples)
        self.assertRaises(TypeError, shape.Interval, None, int(1), samples)
        self.assertRaises(TypeError, shape.Interval, int(1), None, samples)
        self.assertRaises(TypeError, shape.Interval, None, "string", samples)
        self.assertRaises(TypeError, shape.Interval, "string", None, samples)
        self.assertRaises(TypeError, shape.Interval, "string1", "string2",
                          int(1))
        self.assertRaises(ValueError, shape.Interval, 2.0, 1.0, None)
        # validate errors
        self.assertRaises(ValueError, shape.Interval, lower, lower, [])
        self.assertRaises(ValueError, shape.Interval, lower1, upper,
                          invalid_samples_lower_mismatch)
        self.assertRaises(ValueError, shape.Interval, lower, upper,
                          invalid_samples_upper_mismatch)
        self.assertRaises(ValueError, shape.Interval, lower, upper2,
                          invalid_samples_middle_bounds_overlap)

        # test cannot set interval with upper < lower
        interval = shape.Interval(lower, upper2, samples)
        has_assertionError = False
        try:
            interval.upper = 0.5
        except ValueError:
            has_assertionError = True
        self.assertEqual(has_assertionError, True)

        # test intervals in samples
        actual_samples = interval.samples

        actual_subInterval = actual_samples[0]
        expected_subInterval = samples[0]
        actual_lower = actual_subInterval.lower
        actual_upper = actual_subInterval.upper
        expected_lower = expected_subInterval.lower
        expected_upper = expected_subInterval.upper
        self.assertEqual(actual_lower, expected_lower)
        self.assertEqual(actual_upper, expected_upper)

        actual_subInterval = actual_samples[1]
        expected_subInterval = samples[1]
        actual_lower = actual_subInterval.lower
        actual_upper = actual_subInterval.upper
        expected_lower = expected_subInterval.lower
        expected_upper = expected_subInterval.upper
        self.assertEqual(actual_lower, expected_lower)
        self.assertEqual(actual_upper, expected_upper)

        actual_subInterval = actual_samples[2]
        expected_subInterval = samples[2]
        actual_lower = actual_subInterval.lower
        actual_upper = actual_subInterval.upper
        expected_lower = expected_subInterval.lower
        expected_upper = expected_subInterval.upper
        self.assertEqual(actual_lower, expected_lower)
        self.assertEqual(actual_upper, expected_upper)

        # test instance methods
        i1 = shape.Interval(10.0, 15.0)
        self.assertEqual(i1.get_width(), 5)

        # test class methods
        i1 = shape.Interval(10.0, 15.0)
        i2 = shape.Interval(5.0, 8.0)
        intersec1 = shape.Interval.intersection(i1, i2)
        self.assertEqual(intersec1, None)
        intersec2 = shape.Interval.intersection(i2, i1)
        self.assertEqual(intersec2, None)
        i3 = shape.Interval(8.0, 12.0)
        lb = max(i1.lower, i3.lower)
        ub = min(i1.upper, i3.upper)
        intersec3 = shape.Interval.intersection(i1, i3)
        self.assertEqual(intersec3, shape.Interval(lb, ub))


class TestPoint(unittest.TestCase):
    def test_all(self):
        self.assertRaises(TypeError, shape.Point, None, None)
        self.assertRaises(TypeError, shape.Point, None, 1.0)
        self.assertRaises(TypeError, shape.Point, 1.0, None)
        self.assertRaises(TypeError, shape.Point, "string", int(1))
        self.assertRaises(TypeError, shape.Point, int(1), "string")

        point = shape.Point(1.0, 2.0)
        self.assertEqual(point.cval1, 1.0)
        self.assertEqual(point.cval2, 2.0)


class TestSubInterval(unittest.TestCase):
    def test_all(self):

        self.assertRaises(TypeError, shape.SubInterval, None, None)
        self.assertRaises(TypeError, shape.SubInterval, None, 1.0)
        self.assertRaises(TypeError, shape.SubInterval, 1.0, None)
        self.assertRaises(TypeError, shape.SubInterval, "string1", "string2")
        self.assertRaises(ValueError, shape.SubInterval, 2.0, 1.0)

        # test cannot set subInterval with upper < lower
        subInterval = shape.SubInterval(1.0, 2.0)
        has_assertionError = False
        try:
            subInterval.upper = 0.5
        except ValueError:
            has_assertionError = True
        self.assertEqual(has_assertionError, True)

        # test construction method
        shape.SubInterval(10.0, 15.0)


class TestVertex():
    def test_all(self):
        pytest.raises(TypeError, shape.Vertex, None, None, None)
        pytest.raises(TypeError, shape.Vertex, 1.0, 2.0, None)
        pytest.raises(TypeError, shape.Vertex, 1.0, 2.0, 1.0)
        pytest.raises(TypeError, shape.Vertex, None, None,
                      shape.SegmentType.LINE)
        pytest.raises(TypeError, shape.Vertex, None, 2.0,
                      shape.SegmentType.LINE)
        pytest.raises(TypeError, shape.Vertex, 1.0, None,
                      shape.SegmentType.LINE)
        pytest.raises(TypeError, shape.Vertex, None, "string",
                      shape.SegmentType.LINE)
        pytest.raises(TypeError, shape.Vertex, "string", None,
                      shape.SegmentType.LINE)
        pytest.raises(TypeError, shape.Vertex, None, int(1),
                      shape.SegmentType.LINE)
        pytest.raises(TypeError, shape.Vertex, int(1), None,
                      shape.SegmentType.LINE)

        vertex = shape.Vertex(1.0, 2.0, shape.SegmentType.LINE)
        assert(vertex.cval1 == 1.0)
        assert(vertex.cval2 == 2.0)
        assert(vertex.type == shape.SegmentType.LINE)
