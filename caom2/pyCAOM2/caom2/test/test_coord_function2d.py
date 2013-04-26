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

""" Defines TestCoordFunction2D class """

from caom2.wcs.caom2_coord_function2d import CoordFunction2D
from caom2.wcs.caom2_dimension2d import Dimension2D
from caom2.wcs.caom2_ref_coord import RefCoord
from caom2.wcs.caom2_coord2d import Coord2D
import os.path
import sys
import unittest

# put build at the start of the search path
sys.path.insert(0, os.path.abspath('../../lib.local/lib'))


class TestCoordFunction2D(unittest.TestCase):

    def testInit(self):

        dimension = Dimension2D(long(1), long(2))
        ref_coord = Coord2D(RefCoord(float(9.0), float(10.0)),
                            RefCoord(float(11.0), float(12.0)))
        cd11 = float(1.1)
        cd12 = float(1.2)
        cd21 = float(2.1)
        cd22 = float(2.2)

        self.assertRaises(TypeError, CoordFunction2D, None,
                          ref_coord, cd11, cd12, cd21, cd22)
        self.assertRaises(TypeError, CoordFunction2D, dimension,
                          None, cd11, cd12, cd21, cd22)
        self.assertRaises(TypeError, CoordFunction2D, dimension,
                          ref_coord, None, cd12, cd21, cd22)
        self.assertRaises(TypeError, CoordFunction2D, dimension,
                          ref_coord, cd11, None, cd21, cd22)
        self.assertRaises(TypeError, CoordFunction2D, dimension,
                          ref_coord, cd11, cd12, None, cd22)
        self.assertRaises(TypeError, CoordFunction2D, dimension,
                          ref_coord, cd11, cd12, cd21, None)

        function = CoordFunction2D(dimension, ref_coord,
                                   cd11, cd12, cd21, cd22)
        self.assertEqual(function.dimension, dimension)
        self.assertEqual(function.ref_coord, ref_coord)
        self.assertEqual(function.cd11, cd11)
        self.assertEqual(function.cd12, cd12)
        self.assertEqual(function.cd21, cd21)
        self.assertEqual(function.cd22, cd22)

suite = unittest.TestLoader().loadTestsFromTestCase(TestCoordFunction2D)
unittest.TextTestRunner(verbosity=2).run(suite)
