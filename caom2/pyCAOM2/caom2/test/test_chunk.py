#!/usr/bin/env python
# -*- coding: latin-1 -*-
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

""" Defines TestChunk class """

from caom2.caom2_enums import ProductType
from caom2.caom2_chunk import Chunk
from caom2.wcs.caom2_axis import Axis
from caom2.wcs.caom2_coord_axis1d import CoordAxis1D
from caom2.wcs.caom2_coord_axis2d import CoordAxis2D
from caom2.wcs.caom2_slice import Slice
from caom2.wcs.caom2_observable_axis import ObservableAxis
from caom2.wcs.caom2_spatial_wcs import SpatialWCS
from caom2.wcs.caom2_spectral_wcs import SpectralWCS
from caom2.wcs.caom2_temporal_wcs import TemporalWCS
from caom2.wcs.caom2_polarization_wcs import PolarizationWCS
import os.path
import sys
import unittest

# put build at the start of the search path
sys.path.insert(0, os.path.abspath('../../lib.local/lib'))


class TestChunk(unittest.TestCase):

    def testInit(self):

        chunk = Chunk()
        self.assertIsNone(chunk.product_type)
        self.assertIsNone(chunk.naxis)
        #self.assertIsNone(chunk.observable_axis)
        self.assertIsNone(chunk.position_axis_1)
        self.assertIsNone(chunk.position_axis_2)
        self.assertIsNone(chunk.energy_axis)
        self.assertIsNone(chunk.time_axis)
        self.assertIsNone(chunk.polarization_axis)
        self.assertIsNone(chunk.observable)
        self.assertIsNone(chunk.position)
        self.assertIsNone(chunk.energy)
        self.assertIsNone(chunk.time)
        self.assertIsNone(chunk.polarization)

    def testAttributes(self):

        chunk = Chunk()
        with self.assertRaises(TypeError):
            chunk.product_type = float(1.0)
            chunk.naxis = float(1.0)
            #chunk.observable_axis = float(1.0)
            chunk.position_axis_1 = float(1.0)
            chunk.position_axis_2 = float(1.0)
            chunk.energy_axis = float(1.0)
            chunk.time_axis = float(1.0)
            chunk.polarization_axis = float(1.0)
            chunk.observable = float(1.0)
            chunk.position = float(1.0)
            chunk.energy = float(1.0)
            chunk.time = float(1.0)
            chunk.polarization = float(1.0)

        chunk.product_type = ProductType.SCIENCE
        self.assertEqual(ProductType.SCIENCE, chunk.product_type)

        chunk.naxis = int(5)
        self.assertEqual(int(5), chunk.naxis)

        #chunk.observable_axis = int(2)
        #self.assertEqual(int(2), chunk.observable_axis)

        chunk.position_axis_1 = int(1)
        self.assertEqual(int(1), chunk.position_axis_1)

        chunk.position_axis_2 = int(2)
        self.assertEqual(int(2), chunk.position_axis_2)

        chunk.energy_axis = int(3)
        self.assertEqual(int(3), chunk.energy_axis)

        chunk.time_axis = int(4)
        self.assertEqual(int(4), chunk.time_axis)

        chunk.polarization_axis = int(5)
        self.assertEqual(int(5), chunk.polarization_axis)

        axis = Axis("ctype", "cunit")
        dependent = Slice(axis, long(1))
        observable = ObservableAxis(dependent)
        chunk.observable = observable
        self.assertEqual(observable, chunk.observable)

        axis1 = Axis("ctype1", "cunit1")
        axis2 = Axis("ctype2", "cunit2")
        axis_2d = CoordAxis2D(axis1, axis2)
        position = SpatialWCS(axis_2d)
        chunk.position = position
        self.assertEqual(position, chunk.position)

        axis_1d = CoordAxis1D(axis)
        energy = SpectralWCS(axis_1d, "specsys")
        chunk.energy = energy
        self.assertEqual(energy, chunk.energy)

        time = TemporalWCS(axis_1d)
        chunk.time = time
        self.assertEqual(time, chunk.time)

        polarization = PolarizationWCS(axis_1d)
        chunk.polarization = polarization
        self.assertEqual(polarization, chunk.polarization)

#    def testCompareTo(self):
#        # test for chunk1 == chunk2
#        chunk1 = Chunk()
#        chunk1.naxis = 1
#        chunk1.observableAxis = 2
#        chunk1.positionAxis1 = 3
#        chunk1.positionAxis2 = 4
#        chunk1.energyAxis = 5
#        chunk1.timeAxis = 6
#        chunk1.polarizationAxis = 7
#        chunk1.observable = ObservableAxis()
#        chunk1.position = SpatialWCS()
#        chunk1.energy = SpectralWCS()
#        chunk1.time = TemporalWCS()
#        chunk1.polarization = PolarizationWCS()

#        chunk2 = Chunk()
#        chunk2.naxis = 1
#        chunk2.observableAxis = 2
#        chunk2.positionAxis1 = 3
#        chunk2.positionAxis2 = 4
#        chunk2.energyAxis = 5
#        chunk2.timeAxis = 6
#        chunk2.polarizationAxis = 7
#        chunk2.observable = ObservableAxis()
#        chunk2.position = SpatialWCS()
#        chunk2.energy = SpectralWCS()
#        chunk2.time = TemporalWCS()
#        chunk2.polarization = PolarizationWCS()

        # test for chunk1 < chunk2
#        chunk1 = Chunk()
#        chunk1.naxis = 1
#        chunk1.observableAxis = 2
#        chunk1.positionAxis1 = 3
#        chunk1.positionAxis2 = 4
#        chunk1.energyAxis = 5
#        chunk1.timeAxis = 6
#        chunk1.polarizationAxis = 7
#        chunk1.observable = ObservableAxis()
#        chunk1.position = SpatialWCS()
#        chunk1.energy = SpectralWCS()
#        chunk1.time = TemporalWCS()
#        chunk1.polarization = PolarizationWCS()

#        chunk2 = Chunk()
#        chunk2.naxis = 2
#        chunk2.observableAxis = 2
#        chunk2.positionAxis1 = 3
#        chunk2.positionAxis2 = 4
#        chunk2.energyAxis = 5
#        chunk2.timeAxis = 6
#        chunk2.polarizationAxis = 7
#        chunk2.observable = ObservableAxis()
#        chunk2.position = SpatialWCS()
#        chunk2.energy = SpectralWCS()
#        chunk2.time = TemporalWCS()
#        chunk2.polarization = PolarizationWCS()
#
#        self.assertEqual(chunk1.compareTo(chunk2), -1,
#                         "compareTo equal failed")
#
        # test for chunk1 > chunk2
#        chunk1 = Chunk()
#        chunk1.naxis = 2
#        chunk1.observableAxis = 2
#        chunk1.positionAxis1 = 3
#        chunk1.positionAxis2 = 4
#        chunk1.energyAxis = 5
#        chunk1.timeAxis = 6
#        chunk1.polarizationAxis = 7
#        chunk1.observable = ObservableAxis()
#        chunk1.position = SpatialWCS()
#        chunk1.energy = SpectralWCS()
#        chunk1.time = TemporalWCS()
#        chunk1.polarization = PolarizationWCS()
#
#        chunk2 = Chunk()
#        chunk2.naxis = 1
#        chunk2.observableAxis = 2
#        chunk2.positionAxis1 = 3
#        chunk2.positionAxis2 = 4
#        chunk2.energyAxis = 5
#        chunk2.timeAxis = 6
#        chunk2.polarizationAxis = 7
#        chunk2.observable = ObservableAxis()
#        chunk2.position = SpatialWCS()
#        chunk2.energy = SpectralWCS()
#        chunk2.time = TemporalWCS()
#        chunk2.polarization = PolarizationWCS()

suite = unittest.TestLoader().loadTestsFromTestCase(TestChunk)
unittest.TextTestRunner(verbosity=2).run(suite)
