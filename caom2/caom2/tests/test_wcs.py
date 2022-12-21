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

import unittest

from builtins import int

from .. import wcs


class TestAxis(unittest.TestCase):
    def test_init(self):
        self.assertRaises(TypeError, wcs.Axis, None, None)
        self.assertRaises(TypeError, wcs.Axis, None, "cunit")
        self.assertRaises(TypeError, wcs.Axis, "ctype", int(1))
        self.assertRaises(TypeError, wcs.Axis, int(1), "cunit")

        axis = wcs.Axis("ctype", "cunit")
        self.assertEqual(axis.ctype, "ctype")
        self.assertEqual(axis.cunit, "cunit")


class TestCoord2D(unittest.TestCase):
    def test_init(self):
        coord1 = wcs.RefCoord(float(1.0), float(2.0))
        coord2 = wcs.RefCoord(float(3.0), float(4.0))

        self.assertRaises(TypeError, wcs.Coord2D, None, None)
        self.assertRaises(TypeError, wcs.Coord2D, None, coord2)
        self.assertRaises(TypeError, wcs.Coord2D, coord1, None)
        self.assertRaises(TypeError, wcs.Coord2D, str("s"), coord2)
        self.assertRaises(TypeError, wcs.Coord2D, coord1, str("s"))

        coord_2d = wcs.Coord2D(coord1, coord2)
        self.assertEqual(coord_2d.coord1, coord1)
        self.assertEqual(coord_2d.coord2, coord2)


class TestCoordAxis1D(unittest.TestCase):
    def test_init(self):
        self.assertRaises(TypeError, wcs.CoordAxis1D, None)
        self.assertRaises(TypeError, wcs.CoordAxis1D, int(1))

        axis = wcs.Axis("ctype", "cunit")
        axis_1d = wcs.CoordAxis1D(axis)
        self.assertEqual(axis_1d.axis, axis)
        with self.assertRaises(TypeError):
            axis_1d.error = str("s")
            axis_1d.bounds = str("s")
            axis_1d.function = str("s")
            axis_1d.range = str("s")

        error = wcs.CoordError(float(1.0), float(2.0))
        axis_1d.error = error
        self.assertEqual(axis_1d.error, error)

        start = wcs.RefCoord(float(1.0), float(2.0))
        end = wcs.RefCoord(float(3.0), float(4.0))
        coord_range = wcs.CoordRange1D(start, end)
        axis_1d.range = coord_range
        self.assertEqual(axis_1d.range, coord_range)

        bounds = wcs.CoordBounds1D()
        axis_1d.bounds = bounds
        self.assertEqual(axis_1d.bounds, bounds)

        naxis = int(1)
        delta = float(2.5)
        ref_coord = wcs.RefCoord(float(1.0), float(2.0))
        function = wcs.CoordFunction1D(naxis, delta, ref_coord)
        axis_1d.function = function
        self.assertEqual(axis_1d.function, function)


class TestCoordAxis2D(unittest.TestCase):
    def test_init(self):
        self.assertRaises(TypeError, wcs.CoordAxis2D, None, None)
        self.assertRaises(TypeError, wcs.CoordAxis2D, None, int(1))
        self.assertRaises(TypeError, wcs.CoordAxis2D, int(1), None)

        axis1 = wcs.Axis("ctype1", "cunit1")
        axis2 = wcs.Axis("ctype2", "cunit2")
        axis_2d = wcs.CoordAxis2D(axis1, axis2)
        self.assertEqual(axis_2d.axis1, axis1)
        self.assertEqual(axis_2d.axis2, axis2)
        with self.assertRaises(TypeError):
            axis_2d.error1 = str("s")
            axis_2d.error2 = str("s")
            axis_2d.bounds = str("s")
            axis_2d.function = str("s")
            axis_2d.range = str("s")

        error1 = wcs.CoordError(float(1.0), float(2.0))
        axis_2d.error1 = error1
        self.assertEqual(axis_2d.error1, error1)

        error2 = wcs.CoordError(float(3.0), float(4.0))
        axis_2d.error2 = error2
        self.assertEqual(axis_2d.error2, error2)

        start = wcs.Coord2D(wcs.RefCoord(float(1.0), float(2.0)),
                            wcs.RefCoord(float(3.0), float(4.0)))
        end = wcs.Coord2D(wcs.RefCoord(float(5.0), float(6.0)),
                          wcs.RefCoord(float(7.0), float(8.0)))
        coord_range = wcs.CoordRange2D(start, end)
        axis_2d.range = coord_range
        self.assertEqual(axis_2d.range, coord_range)

        center = wcs.ValueCoord2D(float(1.0), float(2.0))
        radius = float(1.5)
        circle = wcs.CoordCircle2D(center, radius)
        axis_2d.bounds = circle
        self.assertEqual(axis_2d.bounds, circle)

        polygon = wcs.CoordPolygon2D()
        axis_2d.bounds = polygon
        self.assertEqual(axis_2d.bounds, polygon)

        dimension = wcs.Dimension2D(int(1), int(2))
        ref_coord = wcs.Coord2D(wcs.RefCoord(float(9.0), float(10.0)),
                                wcs.RefCoord(float(11.0), float(12.0)))
        cd11 = float(1.1)
        cd12 = float(1.2)
        cd21 = float(2.1)
        cd22 = float(2.2)
        function = wcs.CoordFunction2D(dimension, ref_coord,
                                       cd11, cd12, cd21, cd22)
        axis_2d.function = function
        self.assertEqual(axis_2d.function, function)


class TestCoordBounds1D(unittest.TestCase):
    def test_init(self):
        start = wcs.RefCoord(float(1.0), float(2.0))
        end = wcs.RefCoord(float(3.0), float(4.0))
        coord_range = wcs.CoordRange1D(start, end)

        bounds = wcs.CoordBounds1D()
        bounds.samples.append(coord_range)
        self.assertTrue(bounds.samples.count(coord_range) == 1)
        self.assertEqual(bounds.samples.pop(), coord_range)

        with self.assertRaises(TypeError):
            bounds.samples = [str("s")]


class TestCoordBounds2D(unittest.TestCase):
    def test_init(self):
        self.assertRaises(TypeError, wcs.CoordBounds2D, None)
        self.assertRaises(TypeError, wcs.CoordBounds2D, float(1.0))

        center = wcs.ValueCoord2D(float(1.0), float(2.0))
        radius = float(1.5)
        circle = wcs.CoordCircle2D(center, radius)

        polygon = wcs.CoordPolygon2D()
        polygon.vertices.append(wcs.ValueCoord2D(float(1.0), float(2.0)))

        bounds = wcs.CoordBounds2D(circle)
        self.assertEqual(bounds.bounds, circle)

        bounds = wcs.CoordBounds2D(polygon)
        self.assertEqual(bounds.bounds, polygon)


class TestCoordCircle2D(unittest.TestCase):
    def test_init(self):
        center = wcs.ValueCoord2D(float(1.0), float(2.0))
        radius = float(1.5)

        self.assertRaises(TypeError, wcs.CoordCircle2D, None, None)
        self.assertRaises(TypeError, wcs.CoordCircle2D, None, radius)
        self.assertRaises(TypeError, wcs.CoordCircle2D, center, None)
        self.assertRaises(TypeError, wcs.CoordCircle2D, int(1), radius)
        self.assertRaises(TypeError, wcs.CoordCircle2D, center, int(1))

        circle = wcs.CoordCircle2D(center, radius)
        self.assertEqual(circle.center, center)
        self.assertEqual(circle.radius, radius)


class TestCoordError(unittest.TestCase):
    def test_init(self):
        self.assertRaises(TypeError, wcs.CoordError, None, None)
        self.assertRaises(TypeError, wcs.CoordError, None, float(1.0))
        self.assertRaises(TypeError, wcs.CoordError, float(1.0), None)
        self.assertRaises(TypeError, wcs.CoordError, int(1), float(1.0))
        self.assertRaises(TypeError, wcs.CoordError, float(1.0), int(1))

        error = wcs.CoordError(float(1), float(2))
        self.assertIsNotNone(error)
        self.assertEqual(error.syser, float(1))
        self.assertEqual(error.rnder, float(2))


class TestCoordFunction1D(unittest.TestCase):
    def test_init(self):
        naxis = int(1)
        delta = float(2.5)
        ref_coord = wcs.RefCoord(float(1.0), float(2.0))

        self.assertRaises(TypeError, wcs.CoordFunction1D, None, None,
                          None)
        self.assertRaises(TypeError, wcs.CoordFunction1D, None, delta,
                          ref_coord)
        self.assertRaises(TypeError, wcs.CoordFunction1D, naxis, None,
                          ref_coord)
        self.assertRaises(TypeError, wcs.CoordFunction1D, naxis, delta,
                          None)
        self.assertRaises(TypeError, wcs.CoordFunction1D, 'w', delta,
                          ref_coord)
        self.assertRaises(TypeError, wcs.CoordFunction1D, naxis, 'w',
                          ref_coord)
        self.assertRaises(TypeError, wcs.CoordFunction1D, naxis, delta,
                          'w')

        function = wcs.CoordFunction1D(naxis, delta, ref_coord)
        self.assertEqual(function.naxis, naxis)
        self.assertEqual(function.delta, delta)
        self.assertEqual(function.ref_coord, ref_coord)


class TestCoordFunction2D(unittest.TestCase):
    def test_init(self):
        dimension = wcs.Dimension2D(int(1), int(2))
        ref_coord = wcs.Coord2D(wcs.RefCoord(float(9.0), float(10.0)),
                                wcs.RefCoord(float(11.0), float(12.0)))
        cd11 = float(1.1)
        cd12 = float(1.2)
        cd21 = float(2.1)
        cd22 = float(2.2)

        self.assertRaises(TypeError, wcs.CoordFunction2D, None,
                          ref_coord, cd11, cd12, cd21, cd22)
        self.assertRaises(TypeError, wcs.CoordFunction2D, dimension,
                          None, cd11, cd12, cd21, cd22)
        self.assertRaises(TypeError, wcs.CoordFunction2D, dimension,
                          ref_coord, None, cd12, cd21, cd22)
        self.assertRaises(TypeError, wcs.CoordFunction2D, dimension,
                          ref_coord, cd11, None, cd21, cd22)
        self.assertRaises(TypeError, wcs.CoordFunction2D, dimension,
                          ref_coord, cd11, cd12, None, cd22)
        self.assertRaises(TypeError, wcs.CoordFunction2D, dimension,
                          ref_coord, cd11, cd12, cd21, None)

        function = wcs.CoordFunction2D(dimension, ref_coord,
                                       cd11, cd12, cd21, cd22)
        self.assertEqual(function.dimension, dimension)
        self.assertEqual(function.ref_coord, ref_coord)
        self.assertEqual(function.cd11, cd11)
        self.assertEqual(function.cd12, cd12)
        self.assertEqual(function.cd21, cd21)
        self.assertEqual(function.cd22, cd22)


class TestCoordPolygon2D(unittest.TestCase):
    def test_init(self):
        value_coord2d = wcs.ValueCoord2D(float(1.0), float(2.0))

        polygon = wcs.CoordPolygon2D()
        polygon.vertices.append(value_coord2d)
        self.assertTrue(polygon.vertices.count(value_coord2d) == 1)
        self.assertEqual(polygon.vertices.pop(), value_coord2d)

        with self.assertRaises(TypeError):
            polygon.vertices = [str("s")]


class TestCoordRange1D(unittest.TestCase):
    def test_init(self):
        start = wcs.RefCoord(float(1.0), float(2.0))
        end = wcs.RefCoord(float(3.0), float(4.0))

        self.assertRaises(TypeError, wcs.CoordRange1D, None, None)
        self.assertRaises(TypeError, wcs.CoordRange1D, None, end)
        self.assertRaises(TypeError, wcs.CoordRange1D, start, None)
        self.assertRaises(TypeError, wcs.CoordRange1D, int(1), end)
        self.assertRaises(TypeError, wcs.CoordRange1D, start, int(1))

        coord_range = wcs.CoordRange1D(start, end)
        self.assertEqual(coord_range.start, start)
        self.assertEqual(coord_range.end, end)


class TestCoordRange2D(unittest.TestCase):
    def test_init(self):
        start = wcs.Coord2D(wcs.RefCoord(float(1.0), float(2.0)),
                            wcs.RefCoord(float(3.0), float(4.0)))
        end = wcs.Coord2D(wcs.RefCoord(float(5.0), float(6.0)),
                          wcs.RefCoord(float(7.0), float(8.0)))

        self.assertRaises(TypeError, wcs.CoordRange2D, None, None)
        self.assertRaises(TypeError, wcs.CoordRange2D, None, end)
        self.assertRaises(TypeError, wcs.CoordRange2D, start, None)
        self.assertRaises(TypeError, wcs.CoordRange2D, int(1), end)
        self.assertRaises(TypeError, wcs.CoordRange2D, start, int(1))

        coord_range = wcs.CoordRange2D(start, end)
        self.assertEqual(coord_range.start, start)
        self.assertEqual(coord_range.end, end)


class TestDimension2D(unittest.TestCase):
    def test_init(self):
        self.assertRaises(TypeError, wcs.Dimension2D, None, None)
        self.assertRaises(TypeError, wcs.Dimension2D, int(1), None)
        self.assertRaises(TypeError, wcs.Dimension2D, None, int(1))
        self.assertRaises(TypeError, wcs.Dimension2D, 'w', 1)
        self.assertRaises(TypeError, wcs.Dimension2D, 1, 'w')

        dimension = wcs.Dimension2D(int(1), int(2))
        self.assertEqual(dimension.naxis1, int(1))
        self.assertEqual(dimension.naxis2, int(2))


class TestRefCoord(unittest.TestCase):
    def test_init(self):
        self.assertRaises(TypeError, wcs.RefCoord, None, None)
        self.assertRaises(TypeError, wcs.RefCoord, None, float(1.0))
        self.assertRaises(TypeError, wcs.RefCoord, float(1.0), None)
        self.assertRaises(TypeError, wcs.RefCoord, int(1), float(1.0))
        self.assertRaises(TypeError, wcs.RefCoord, float(1.0), int(1))

        ref_coord = wcs.RefCoord(float(1), float(2))
        self.assertIsNotNone(ref_coord)
        self.assertEqual(ref_coord.pix, float(1))
        self.assertEqual(ref_coord.val, float(2))


class TestSlice(unittest.TestCase):
    def test_init(self):
        axis = wcs.Axis("ctype", "cunit")
        my_bin = int(1)

        self.assertRaises(TypeError, wcs.Slice, None, None)
        self.assertRaises(TypeError, wcs.Slice, None, my_bin)
        self.assertRaises(TypeError, wcs.Slice, axis, None)
        self.assertRaises(TypeError, wcs.Slice, str("s"), my_bin)
        self.assertRaises(TypeError, wcs.Slice, axis, 'a')

        my_slice = wcs.Slice(axis, my_bin)
        self.assertEqual(my_slice.axis, axis)
        self.assertEqual(my_slice.bin, int(1))


class TestValueCoord2d(unittest.TestCase):
    def test_init(self):
        self.assertRaises(TypeError, wcs.ValueCoord2D, None, None)
        self.assertRaises(TypeError, wcs.ValueCoord2D, None, float(1.0))
        self.assertRaises(TypeError, wcs.ValueCoord2D, float(1.0), None)
        self.assertRaises(TypeError, wcs.ValueCoord2D, 1, float(1.0))
        self.assertRaises(TypeError, wcs.ValueCoord2D, float(1.0), 1)

        value_coord2d = wcs.ValueCoord2D(float(1), float(2))
        self.assertIsNotNone(value_coord2d)
        self.assertEqual(value_coord2d.coord1, float(1))
        self.assertEqual(value_coord2d.coord2, float(2))
