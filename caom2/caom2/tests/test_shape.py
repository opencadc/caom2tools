# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2025.                            (c) 2025.
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
        assert (vertex.cval1 == 1.0)
        assert (vertex.cval2 == 2.0)
        assert (vertex.type == shape.SegmentType.LINE)
