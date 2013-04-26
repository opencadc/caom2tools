#!/usr/bin/env python
# -*- coding: utf-8 -*-
#***********************************************************************
#******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
#*************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2010.                            (c) 2010.
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
#***********************************************************************
#

""" Defines TestCoordAxis2D class """

from caom2.wcs.caom2_axis import Axis
from caom2.wcs.caom2_coord2d import Coord2D
from caom2.wcs.caom2_dimension2d import Dimension2D
from caom2.wcs.caom2_coord_axis2d import CoordAxis2D
from caom2.wcs.caom2_coord_circle2d import CoordCircle2D
from caom2.wcs.caom2_coord_polygon2d import CoordPolygon2D
from caom2.wcs.caom2_coord_error import CoordError
from caom2.wcs.caom2_coord_function2d import CoordFunction2D
from caom2.wcs.caom2_coord_range2d import CoordRange2D
from caom2.wcs.caom2_ref_coord import RefCoord
import os.path
import sys
import unittest

# put build at the start of the search path
sys.path.insert(0, os.path.abspath('../../lib.local/lib'))


class TestCoordAxis2D(unittest.TestCase):

    def testInit(self):

        self.assertRaises(TypeError, CoordAxis2D, None, None)
        self.assertRaises(TypeError, CoordAxis2D, None, int(1))
        self.assertRaises(TypeError, CoordAxis2D, int(1), None)

        axis1 = Axis("ctype1", "cunit1")
        axis2 = Axis("ctype2", "cunit2")
        axis_2d = CoordAxis2D(axis1, axis2)
        self.assertEqual(axis_2d.axis1, axis1)
        self.assertEqual(axis_2d.axis2, axis2)
        with self.assertRaises(TypeError):
            axis_2d.error1 = str("s")
            axis_2d.error2 = str("s")
            axis_2d.bounds = str("s")
            axis_2d.function = str("s")
            axis_2d.range = str("s")

        error1 = CoordError(float(1.0), float(2.0))
        axis_2d.error1 = error1
        self.assertEqual(axis_2d.error1, error1)

        error2 = CoordError(float(3.0), float(4.0))
        axis_2d.error2 = error2
        self.assertEqual(axis_2d.error2, error2)

        start = Coord2D(RefCoord(float(1.0), float(2.0)),
                        RefCoord(float(3.0), float(4.0)))
        end = Coord2D(RefCoord(float(5.0), float(6.0)),
                      RefCoord(float(7.0), float(8.0)))
        coordRange = CoordRange2D(start, end)
        axis_2d.range = coordRange
        self.assertEqual(axis_2d.range, coordRange)

        center = Coord2D(RefCoord(float(1.0), float(2.0)),
                         RefCoord(float(3.0), float(4.0)))
        radius = float(1.5)
        circle = CoordCircle2D(center, radius)
        axis_2d.bounds = circle
        self.assertEqual(axis_2d.bounds, circle)

        polygon = CoordPolygon2D()
        axis_2d.bounds = polygon
        self.assertEqual(axis_2d.bounds, polygon)

        dimension = Dimension2D(long(1), long(2))
        ref_coord = Coord2D(RefCoord(float(9.0), float(10.0)),
                            RefCoord(float(11.0), float(12.0)))
        cd11 = float(1.1)
        cd12 = float(1.2)
        cd21 = float(2.1)
        cd22 = float(2.2)
        function = CoordFunction2D(dimension, ref_coord,
                                   cd11, cd12, cd21, cd22)
        axis_2d.function = function
        self.assertEqual(axis_2d.function, function)

suite = unittest.TestLoader().loadTestsFromTestCase(TestCoordAxis2D)
unittest.TextTestRunner(verbosity=2).run(suite)
