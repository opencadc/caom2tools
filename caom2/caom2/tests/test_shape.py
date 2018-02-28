# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2016.                            (c) 2016.
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

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import math
import pytest
import unittest

from .. import shape


class TestEnums(unittest.TestCase):
    def test_all(self):
        # test for invalid value
        with self.assertRaises(KeyError):
            shape.SegmentType["foo"]
        with self.assertRaises(KeyError):
            shape.SegmentType[None]
        with self.assertRaises(KeyError):
            shape.SegmentType[999]

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
        self.assertRaises(AssertionError, shape.Interval, 2.0, 1.0, None)
        # validate errors
        self.assertRaises(AssertionError, shape.Interval, lower, lower, [])
        self.assertRaises(AssertionError, shape.Interval, lower1, upper,
                          invalid_samples_lower_mismatch)
        self.assertRaises(AssertionError, shape.Interval, lower, upper,
                          invalid_samples_upper_mismatch)
        self.assertRaises(AssertionError, shape.Interval, lower, upper2,
                          invalid_samples_middle_bounds_overlap)

        # test cannot set interval with upper < lower
        interval = shape.Interval(lower, upper2, samples)
        has_assertionError = False
        try:
            interval.upper = 0.5
        except AssertionError:
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
        self.assertRaises(AssertionError, shape.SubInterval, 2.0, 1.0)

        # test cannot set subInterval with upper < lower
        subInterval = shape.SubInterval(1.0, 2.0)
        has_assertionError = False
        try:
            subInterval.upper = 0.5
        except AssertionError:
            has_assertionError = True
        self.assertEqual(has_assertionError, True)

        # test construction method
        shape.SubInterval(10.0, 15.0)


class TestOpenPolygon():
    def test_all(self):
        p1 = shape.Point(-117.246094, 52.942018)
        p2 = shape.Point(-101.601563, 56.535258)
        p3 = shape.Point(-97.382813, 44.809122)
        p4 = shape.Point(-111.445313, 37.405074)
        # SphericalPolygon requires p1 == p5 for a closed polygon
        p5 = shape.Point(-117.246094, 52.942018)
        no_points = []
        too_few_points = [p1, p2]
        min_closed_points = [p1, p2, p3]
        closed_points = [p1, p2, p3, p4, p5]
        counter_clockwise_points = [p4, p3, p2, p1]

        # should detect that the polygons is not clockwise
        with pytest.raises(AssertionError) as ex:
            shape.Polygon(counter_clockwise_points)
        assert('not in clockwise direction' in str(ex.value))
        # should detect that polygon is requires a minimum of 4 points
        with pytest.raises(AssertionError) as ex:
            shape.Polygon(no_points)
        assert('invalid polygon: 0 points' in str(ex.value))
        with pytest.raises(AssertionError) as ex:
            shape.Polygon(too_few_points)
        assert('invalid polygon: 2 points' in str(ex.value))

        # polygon default constructor
        shape.Polygon()

        # should detect that polygon is closed
        shape.Polygon(min_closed_points)
        shape.Polygon(closed_points)

        # should detect that multipolygon is not closed
        v0 = shape.Vertex(-126.210938, 67.991108, shape.SegmentType.MOVE)
        v1 = shape.Vertex(-108.984375, 70.480896, shape.SegmentType.LINE)
        v2 = shape.Vertex(-98.789063, 66.912834, shape.SegmentType.LINE)
        v3 = shape.Vertex(-75.234375, 60.217991, shape.SegmentType.LINE)
        v4 = shape.Vertex(-87.890625, 52.241256, shape.SegmentType.LINE)
        v5 = shape.Vertex(-110.742188, 54.136696, shape.SegmentType.LINE)
        v6 = shape.Vertex(0.0, 0.0, shape.SegmentType.CLOSE)
        v7 = shape.Vertex(24.609375, 62.895218, shape.SegmentType.MOVE)
        v8 = shape.Vertex(43.593750, 67.322924, shape.SegmentType.LINE)
        v9 = shape.Vertex(55.898438, 62.734601, shape.SegmentType.LINE)
        v10 = shape.Vertex(46.757813, 56.145550, shape.SegmentType.LINE)
        v11 = shape.Vertex(26.015625, 55.354135, shape.SegmentType.LINE)
        v12 = shape.Vertex(0.0, 0.0, shape.SegmentType.CLOSE)
        closed_vertices = [
            v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12]

        # should detect that multipolygon is closed
        mp = shape.MultiPolygon(closed_vertices)

        # instantiated polygon should contain the same points
        p = shape.Polygon(points=closed_points, samples=mp)
        actual_points = p.points
        assert(actual_points[0].cval1 == closed_points[0].cval1)
        assert(actual_points[0].cval2 == closed_points[0].cval2)
        assert(actual_points[1].cval1 == closed_points[1].cval1)
        assert(actual_points[1].cval2 == closed_points[1].cval2)
        assert(actual_points[2].cval1 == closed_points[2].cval1)
        assert(actual_points[2].cval2 == closed_points[2].cval2)
        assert(actual_points[3].cval1 == closed_points[3].cval1)
        assert(actual_points[3].cval2 == closed_points[3].cval2)

        assert(mp is p.samples)


class TestOpenMultiPolygon():
    def test_all(self):
        # should detect that multipolygon is not closed
        v0 = shape.Vertex(-126.210938, 67.991108, shape.SegmentType.MOVE)
        v1 = shape.Vertex(-108.984375, 70.480896, shape.SegmentType.LINE)
        v2 = shape.Vertex(-98.789063, 66.912834, shape.SegmentType.LINE)
        v3 = shape.Vertex(-75.234375, 60.217991, shape.SegmentType.LINE)
        v4 = shape.Vertex(-87.890625, 52.241256, shape.SegmentType.LINE)
        v5 = shape.Vertex(-110.742188, 54.136696, shape.SegmentType.LINE)
        v6 = shape.Vertex(0.0, 0.0, shape.SegmentType.CLOSE)
        v7 = shape.Vertex(24.609375, 62.895218, shape.SegmentType.MOVE)
        v8 = shape.Vertex(43.593750, 67.322924, shape.SegmentType.LINE)
        v9 = shape.Vertex(55.898438, 62.734601, shape.SegmentType.LINE)
        v10 = shape.Vertex(46.757813, 56.145550, shape.SegmentType.LINE)
        v11 = shape.Vertex(26.015625, 55.354135, shape.SegmentType.LINE)
        v12 = shape.Vertex(0.0, 0.0, shape.SegmentType.CLOSE)

        no_vertices = []
        too_few_vertices = [v0, v1, v6]
        two_moves_vertices = [v0, v1, v7, v2, v3, v4, v5, v6]
        no_move_vertices = [v1, v2, v3, v4, v5, v6]
        two_closes_vertices = [
            v0, v1, v2, v3, v4, v5, v7, v8, v9, v10, v11, v12]
        no_close_vertices = [v0, v1, v2, v3, v4, v5]
        min_closed_vertices = [v0, v1, v2, v6]
        closed_vertices = [
            v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12]

        rv0 = shape.Vertex(26.015625, 55.354135, shape.SegmentType.MOVE)
        rv1 = shape.Vertex(46.757813, 56.145550, shape.SegmentType.LINE)
        rv2 = shape.Vertex(55.898438, 62.734601, shape.SegmentType.LINE)
        rv3 = shape.Vertex(43.593750, 67.322924, shape.SegmentType.LINE)
        rv4 = shape.Vertex(24.609375, 62.895218, shape.SegmentType.LINE)
        rv5 = shape.Vertex(0.0, 0.0, shape.SegmentType.CLOSE)
        rv6 = shape.Vertex(-110.742188, 54.136696, shape.SegmentType.MOVE)
        rv7 = shape.Vertex(-87.890625, 52.241256, shape.SegmentType.LINE)
        rv8 = shape.Vertex(-75.234375, 60.217991, shape.SegmentType.LINE)
        rv9 = shape.Vertex(-98.789063, 66.912834, shape.SegmentType.LINE)
        rv10 = shape.Vertex(-108.984375, 70.480896, shape.SegmentType.LINE)
        rv11 = shape.Vertex(-126.210938, 67.991108, shape.SegmentType.LINE)
        rv12 = shape.Vertex(0.0, 0.0, shape.SegmentType.CLOSE)
        counter_clockwise_vertices = [
            rv0, rv1, rv2, rv3, rv4, rv5, rv6, rv7, rv8, rv9, rv10, rv11, rv12]

        # should detect that the polygons is not clockwise
        with pytest.raises(AssertionError) as ex:
            shape.MultiPolygon(counter_clockwise_vertices)
        assert('not in clockwise direction' in str(ex.value))
        # should detect that there are not enough number of vertices to
        # produce a multipolygon
        with pytest.raises(AssertionError) as ex:
            shape.MultiPolygon(no_vertices)
        assert('invalid polygon: 0 vertices' in str(ex.value))
        with pytest.raises(AssertionError) as ex:
            shape.MultiPolygon(too_few_vertices)
        assert('invalid polygon: 3 vertices' in str(ex.value))
        # no close between two 'MOVE'
        with pytest.raises(AssertionError) as ex:
            shape.MultiPolygon(two_moves_vertices)
        assert(
            'invalid polygon: MOVE vertex when loop open' in str(ex.value))
        # no 'MOVE' before a 'CLOSE'
        with pytest.raises(AssertionError) as ex:
            shape.MultiPolygon(no_move_vertices)
        assert(
            'invalid polygon: first vertex is not a MOVE' in str(ex.value))
        # no 'MOVE' between two 'CLOSE'
        with pytest.raises(AssertionError) as ex:
            shape.MultiPolygon(two_closes_vertices)
        assert(
            'invalid polygon: MOVE vertex when loop open' in str(ex.value))
        # no 'CLOSE' after a 'MOVE'
        with pytest.raises(AssertionError) as ex:
            shape.MultiPolygon(no_close_vertices)
        assert(
            'invalid polygon: last vertex is not a CLOSE' in str(ex.value))

        # multipolygon default constructor
        shape.MultiPolygon(None)

        # should detect that multipolygon is closed
        shape.MultiPolygon(min_closed_vertices)

        # should detect that multipolygon is closed
        shape.MultiPolygon(closed_vertices)

        # instantiated multipolygon should contain the same vertices
        p = shape.MultiPolygon(vertices=closed_vertices)
        actual_vertices = p.vertices
        assert(actual_vertices[0].cval1 == closed_vertices[0].cval1)
        assert(actual_vertices[0].cval2 == closed_vertices[0].cval2)
        assert(actual_vertices[0].type == shape.SegmentType.MOVE)
        assert(actual_vertices[1].cval1 == closed_vertices[1].cval1)
        assert(actual_vertices[1].cval2 == closed_vertices[1].cval2)
        assert(actual_vertices[1].type == shape.SegmentType.LINE)
        assert(actual_vertices[2].cval1 == closed_vertices[2].cval1)
        assert(actual_vertices[2].cval2 == closed_vertices[2].cval2)
        assert(actual_vertices[2].type == shape.SegmentType.LINE)
        assert(actual_vertices[3].cval1 == closed_vertices[3].cval1)
        assert(actual_vertices[3].cval2 == closed_vertices[3].cval2)
        assert(actual_vertices[3].type == shape.SegmentType.LINE)
        assert(actual_vertices[4].cval1 == closed_vertices[4].cval1)
        assert(actual_vertices[4].cval2 == closed_vertices[4].cval2)
        assert(actual_vertices[4].type == shape.SegmentType.LINE)
        assert(actual_vertices[5].cval1 == closed_vertices[5].cval1)
        assert(actual_vertices[5].cval2 == closed_vertices[5].cval2)
        assert(actual_vertices[5].type == shape.SegmentType.LINE)
        assert(actual_vertices[6].cval1 == closed_vertices[6].cval1)
        assert(actual_vertices[6].cval2 == closed_vertices[6].cval2)
        assert(actual_vertices[6].type  ==shape.SegmentType.CLOSE)
        assert(actual_vertices[7].cval1 == closed_vertices[7].cval1)
        assert(actual_vertices[7].cval2 == closed_vertices[7].cval2)
        assert(actual_vertices[7].type == shape.SegmentType.MOVE)
        assert(actual_vertices[8].cval1 == closed_vertices[8].cval1)
        assert(actual_vertices[8].cval2 == closed_vertices[8].cval2)
        assert(actual_vertices[8].type == shape.SegmentType.LINE)
        assert(actual_vertices[9].cval1 == closed_vertices[9].cval1)
        assert(actual_vertices[9].cval2 == closed_vertices[9].cval2)
        assert(actual_vertices[9].type == shape.SegmentType.LINE)
        assert(actual_vertices[10].cval1 == closed_vertices[10].cval1)
        assert(actual_vertices[10].cval2 == closed_vertices[10].cval2)
        assert(actual_vertices[10].type == shape.SegmentType.LINE)
        assert(actual_vertices[11].cval1 == closed_vertices[11].cval1)
        assert(actual_vertices[11].cval2 == closed_vertices[11].cval2)
        assert(actual_vertices[11].type == shape.SegmentType.LINE)
        assert(actual_vertices[12].cval1 == closed_vertices[12].cval1)
        assert(actual_vertices[12].cval2 == closed_vertices[12].cval2)
        assert(actual_vertices[12].type == shape.SegmentType.CLOSE)


class TestSelfIntersectingPolygon():
    def test_all(self):
        # should detect self segment intersection of the polygon not near a
        # Pole
        p1 = shape.Point(-115.488281, 45.867063)
        p2 = shape.Point(-91.230469, 36.075742)
        p3 = shape.Point(-95.800781, 54.807017)
        p4 = shape.Point(-108.457031, 39.951859)
        p5 = shape.Point(0.0, 0.0)
        points_with_self_intersecting_segments = [p1, p2, p3, p4, p5]
        with pytest.raises(AssertionError) as ex:
            shape.Polygon(points_with_self_intersecting_segments)
        assert('self intersecting' in str(ex.value))

        # should detect self segment intersection of the polygon near the
        # South Pole, with the Pole outside the polygon
        p1 = shape.Point(0.6128286003, -89.8967940441)
        p2 = shape.Point(210.6391743183, -89.9073892376)
        p3 = shape.Point(90.6405151921, -89.8972874698)
        p4 = shape.Point(270.6114701911, -89.90689353)
        p5 = shape.Point(0.0, 0.0)
        points_with_self_intersecting_segments = [p1, p2, p3, p4, p5]
        with pytest.raises(AssertionError) as ex:
            shape.Polygon(points_with_self_intersecting_segments)
        assert('self intersecting' in str(ex.value))

        # should detect self segment intersection of the polygon near the
        # South Pole, with the Pole inside the polygon
        p1 = shape.Point(0.6128286003, -89.8967940441)
        p2 = shape.Point(130.6391743183, -89.9073892376)
        p3 = shape.Point(90.6405151921, -89.8972874698)
        p4 = shape.Point(270.6114701911, -89.90689353)
        p5 = shape.Point(0.0, 0.0)
        points_with_self_intersecting_segments = [p1, p2, p3, p4, p5]
        with pytest.raises(AssertionError) as ex:
            shape.Polygon(points_with_self_intersecting_segments)
        assert('self intersecting' in str(ex.value))

        # should detect self segment intersection of the polygon which
        # intersects with meridian = 0
        p1 = shape.Point(-7.910156, 13.293411)
        p2 = shape.Point(4.042969, 7.068185)
        p3 = shape.Point(4.746094, 18.030975)
        p4 = shape.Point(-6.855469, 6.369894)
        p5 = shape.Point(0.0, 0.0)
        points_with_self_intersecting_segments = [p1, p2, p3, p4, p5]
        with pytest.raises(AssertionError) as ex:
            shape.Polygon(points_with_self_intersecting_segments)
        assert('self intersecting' in str(ex.value))


class TestSelfIntersectingMultiPolygon():
    def test_all(self):
        # should detect self segment intersection of the multipolygon not
        # near a Pole
        v1 = shape.Vertex(-115.488281, 45.867063, shape.SegmentType.MOVE)
        v2 = shape.Vertex(-91.230469, 36.075742, shape.SegmentType.LINE)
        v3 = shape.Vertex(-95.800781, 54.807017, shape.SegmentType.LINE)
        v4 = shape.Vertex(-108.457031, 39.951859, shape.SegmentType.LINE)
        v5 = shape.Vertex(0.0, 0.0, shape.SegmentType.CLOSE)
        points_with_self_intersecting_segments = [v1, v2, v3, v4, v5]
        with pytest.raises(AssertionError) as ex:
            shape.MultiPolygon(points_with_self_intersecting_segments)
        assert('self intersecting' in str(ex.value))

        # should detect self segment intersection of the multipolygon near
        # the South Pole, with the Pole outside the multipolygon
        v1 = shape.Vertex(0.6128286003, -89.8967940441, shape.SegmentType.MOVE)
        v2 = shape.Vertex(
            210.6391743183, -89.9073892376, shape.SegmentType.LINE)
        v3 = shape.Vertex(
            90.6405151921, -89.8972874698, shape.SegmentType.LINE)
        v4 = shape.Vertex(270.6114701911, -89.90689353, shape.SegmentType.LINE)
        v5 = shape.Vertex(0.0, 0.0, shape.SegmentType.CLOSE)
        points_with_self_intersecting_segments = [v1, v2, v3, v4, v5]
        with pytest.raises(AssertionError) as ex:
            shape.Polygon(points_with_self_intersecting_segments)
        assert('self intersecting' in str(ex.value))

        # should detect self segment intersection of the multipolygon near the
        # South Pole, with the Pole inside the multipolygon
        v1 = shape.Vertex(0.6128286003, -89.8967940441, shape.SegmentType.MOVE)
        v2 = shape.Vertex(
            130.6391743183, -89.9073892376, shape.SegmentType.LINE)
        v3 = shape.Vertex(
            90.6405151921, -89.8972874698, shape.SegmentType.LINE)
        v4 = shape.Vertex(270.6114701911, -89.90689353, shape.SegmentType.LINE)
        v5 = shape.Vertex(0.0, 0.0, shape.SegmentType.CLOSE)
        points_with_self_intersecting_segments = [v1, v2, v3, v4, v5]
        with pytest.raises(AssertionError) as ex:
            shape.Polygon(points_with_self_intersecting_segments)
        assert('self intersecting' in str(ex.value))

        # should detect self segment intersection of the multipolygon which
        # intersects with meridian = 0
        v1 = shape.Vertex(-7.910156, 13.293411, shape.SegmentType.MOVE)
        v2 = shape.Vertex(4.042969, 7.068185, shape.SegmentType.LINE)
        v3 = shape.Vertex(4.746094, 18.030975, shape.SegmentType.LINE)
        v4 = shape.Vertex(-6.855469, 6.369894, shape.SegmentType.LINE)
        v5 = shape.Vertex(0.0, 0.0, shape.SegmentType.CLOSE)
        points_with_self_intersecting_segments = [v1, v2, v3, v4, v5]
        with pytest.raises(AssertionError) as ex:
            shape.Polygon(points_with_self_intersecting_segments)
        assert('self intersecting' in str(ex.value))


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
