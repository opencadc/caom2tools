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

""" Defines TestChunk class """

import unittest

from .. import chunk
from .. import wcs


class TestEnums(unittest.TestCase):
    def test_all(self):
        # test for invalid value
        with self.assertRaises(ValueError):
            chunk.ProductType("no_such_string")
        with self.assertRaises(ValueError):
            chunk.ProductType(None)
        with self.assertRaises(ValueError):
            chunk.ProductType(1)

        # test that we can get the object for each enum by name
        self.assertEqual(chunk.ProductType.SCIENCE.name, "SCIENCE")
        self.assertEqual(chunk.ProductType[
                             chunk.ProductType.SCIENCE.name].name, "SCIENCE")
        self.assertEqual(chunk.ProductType['SCIENCE'].value, "science")
        self.assertEqual(chunk.ProductType[
                             chunk.ProductType.SCIENCE.name].value, "science")

        self.assertEqual(chunk.ProductType.SCIENCE.value, "science")
        self.assertEqual(chunk.ProductType.CALIBRATION.value, "calibration")
        self.assertEqual(chunk.ProductType.PREVIEW.value, "preview")
        self.assertEqual(chunk.ProductType.INFO.value, "info")
        self.assertEqual(chunk.ProductType.NOISE.value, "noise")
        self.assertEqual(chunk.ProductType.WEIGHT.value, "weight")
        self.assertEqual(chunk.ProductType.AUXILIARY.value, "auxiliary")
        self.assertEqual(chunk.ProductType.THUMBNAIL.value, "thumbnail")
        self.assertEqual(chunk.ProductType.BIAS.value, "bias")
        self.assertEqual(chunk.ProductType.DARK.value, "dark")
        self.assertEqual(chunk.ProductType.FLAT.value, "flat")
        self.assertEqual(chunk.ProductType.WAVECAL.value, "wavecal")


class TestChunk(unittest.TestCase):
    def test_init(self):
        test_chunk = chunk.Chunk()
        self.assertIsNone(test_chunk.product_type)
        self.assertIsNone(test_chunk.naxis)
        self.assertIsNone(test_chunk.position_axis_1)
        self.assertIsNone(test_chunk.position_axis_2)
        self.assertIsNone(test_chunk.energy_axis)
        self.assertIsNone(test_chunk.time_axis)
        self.assertIsNone(test_chunk.polarization_axis)
        self.assertIsNone(test_chunk.custom_axis)
        self.assertIsNone(test_chunk.observable)
        self.assertIsNone(test_chunk.position)
        self.assertIsNone(test_chunk.energy)
        self.assertIsNone(test_chunk.time)
        self.assertIsNone(test_chunk.polarization)
        self.assertIsNone(test_chunk.custom)

    def test_attributes(self):
        test_chunk = chunk.Chunk()
        with self.assertRaises(TypeError):
            test_chunk.product_type = float(1.0)
            test_chunk.naxis = float(1.0)
            test_chunk.position_axis_1 = float(1.0)
            test_chunk.position_axis_2 = float(1.0)
            test_chunk.energy_axis = float(1.0)
            test_chunk.time_axis = float(1.0)
            test_chunk.polarization_axis = float(1.0)
            test_chunk.custom_axis = float(1.0)
            test_chunk.observable = float(1.0)
            test_chunk.position = float(1.0)
            test_chunk.energy = float(1.0)
            test_chunk.time = float(1.0)
            test_chunk.polarization = float(1.0)
            test_chunk.custom = float(1.0)

        test_chunk.product_type = chunk.ProductType.SCIENCE
        self.assertEqual(chunk.ProductType.SCIENCE.name,
                         test_chunk.product_type.name)

        test_chunk.naxis = int(5)
        self.assertEqual(int(5), test_chunk.naxis)

        test_chunk.position_axis_1 = int(1)
        self.assertEqual(int(1), test_chunk.position_axis_1)

        test_chunk.position_axis_2 = int(2)
        self.assertEqual(int(2), test_chunk.position_axis_2)

        test_chunk.energy_axis = int(3)
        self.assertEqual(int(3), test_chunk.energy_axis)

        test_chunk.time_axis = int(4)
        self.assertEqual(int(4), test_chunk.time_axis)

        test_chunk.polarization_axis = int(5)
        self.assertEqual(int(5), test_chunk.polarization_axis)

        test_chunk.custom_axis = int(6)
        self.assertEqual(int(6), test_chunk.custom_axis)

        axis = wcs.Axis("ctype", "cunit")
        dependent = wcs.Slice(axis, 1)
        observable = chunk.ObservableAxis(dependent)
        test_chunk.observable = observable
        self.assertEqual(observable, test_chunk.observable)

        axis1 = wcs.Axis("ctype1", "cunit1")
        axis2 = wcs.Axis("ctype2", "cunit2")
        axis_2d = wcs.CoordAxis2D(axis1, axis2)
        position = chunk.SpatialWCS(axis_2d)
        test_chunk.position = position
        self.assertEqual(position, test_chunk.position)

        axis_1d = wcs.CoordAxis1D(axis)
        energy = chunk.SpectralWCS(axis_1d, "specsys")
        test_chunk.energy = energy
        self.assertEqual(energy, test_chunk.energy)

        time = chunk.TemporalWCS(axis_1d)
        test_chunk.time = time
        self.assertEqual(time, test_chunk.time)

        polarization = chunk.PolarizationWCS(
            wcs.CoordAxis1D(wcs.Axis('STOKES')))
        test_chunk.polarization = polarization
        self.assertEqual(polarization, test_chunk.polarization)

        custom = chunk.CustomWCS(wcs.CoordAxis1D(axis1))
        test_chunk.custom = custom
        self.assertEqual(custom, test_chunk.custom)


class TestObservableAxis(unittest.TestCase):
    def test_init(self):
        self.assertRaises(TypeError, chunk.ObservableAxis, None)
        self.assertRaises(TypeError, chunk.ObservableAxis, int(1))

        dependent = wcs.Slice(wcs.Axis("ctype1", "cunit1"), 1)
        independent = wcs.Slice(wcs.Axis("ctype2", "cunit2"), 2)

        observable = chunk.ObservableAxis(dependent)
        self.assertEqual(observable.dependent, dependent)

        observable.independent = independent
        self.assertEqual(observable.independent, independent)


class TestSpatialWCS(unittest.TestCase):
    def test_init(self):
        self.assertRaises(TypeError, chunk.SpatialWCS, None)
        self.assertRaises(TypeError, chunk.SpatialWCS, int(1))

        axis1 = wcs.Axis("ctype1", "cunit1")
        axis2 = wcs.Axis("ctype2", "cunit2")
        axis_2d = wcs.CoordAxis2D(axis1, axis2)
        position = chunk.SpatialWCS(axis_2d)
        self.assertEqual(position.axis, axis_2d)
        with self.assertRaises(TypeError):
            position.coordsys = float(1.0)
            position.bounds = str("s")
            position.function = str("s")

        position.coordsys = "coordsys"
        self.assertEqual(position.coordsys, "coordsys")

        with self.assertRaises(ValueError):
            position.equinox = float(1.0)
        position.equinox = float(2000.0)
        self.assertEqual(position.equinox, float(2000.0))

        position.resolution = float(2.0)
        self.assertEqual(position.resolution, float(2.0))


class TestSpectralWCS(unittest.TestCase):
    def test_init(self):
        axis = wcs.Axis("ctype", "cunit")
        axis_1d = wcs.CoordAxis1D(axis)

        self.assertRaises(TypeError, chunk.SpectralWCS, None, None)
        self.assertRaises(TypeError, chunk.SpectralWCS, None, str("s"))
        self.assertRaises(TypeError, chunk.SpectralWCS, axis_1d, None)
        self.assertRaises(TypeError, chunk.SpectralWCS, int(1), str("s"))
        self.assertRaises(TypeError, chunk.SpectralWCS, axis_1d, int(1))

        energy = chunk.SpectralWCS(axis_1d, "specsys")
        self.assertEqual(energy.axis, axis_1d)
        self.assertEqual(energy.specsys, "specsys")
        with self.assertRaises(TypeError):
            energy.ssysobs = int(1)
            energy.ssyssrc = int(1)
            energy.restfrq = int(1)
            energy.restwav = int(1)
            energy.velosys = int(1)
            energy.zsource = int(1)
            energy.velang = int(1)
            energy.bandpass_name = int(1)
            energy.transition = int(1)
            energy.resolving_power = int(1)

        with self.assertRaises(ValueError):
            energy.zsource = float(-1)
            energy.zsource = float(1201)
            energy.resolving_power = float(-1)
            energy.resolving_power = float(1.1e8)

        energy.ssysobs = "ssysobs"
        self.assertEqual(energy.ssysobs, "ssysobs")

        energy.ssyssrc = "ssyssrc"
        self.assertEqual(energy.ssyssrc, "ssyssrc")

        energy.restfrq = float(1.0)
        self.assertEqual(energy.restfrq, float(1.0))

        energy.restwav = float(2.0)
        self.assertEqual(energy.restwav, float(2.0))

        energy.velosys = float(3.0)
        self.assertEqual(energy.velosys, float(3.0))

        energy.zsource = float(4.0)
        self.assertEqual(energy.zsource, float(4.0))

        energy.velang = float(5.0)
        self.assertEqual(energy.velang, float(5.0))

        energy.bandpass_name = "bandpass_name"
        self.assertEqual(energy.bandpass_name, "bandpass_name")

        transition = wcs.EnergyTransition("species", "transition")
        energy.transition = transition
        self.assertEqual(energy.transition, transition)

        energy.resolving_power = float(6.0)
        self.assertEqual(energy.resolving_power, float(6.0))


class TestTemporalWCS(unittest.TestCase):
    def test_init(self):
        self.assertRaises(TypeError, chunk.TemporalWCS, None)
        self.assertRaises(TypeError, chunk.TemporalWCS, int(1))

        axis = wcs.Axis("ctype", "cunit")
        axis_1d = wcs.CoordAxis1D(axis)
        time = chunk.TemporalWCS(axis_1d)
        self.assertEqual(time.axis, axis_1d)
        with self.assertRaises(TypeError):
            time.exposure = str("s")
            time.resolution = str("s")

        time.exposure = float(1.0)
        self.assertEqual(time.exposure, float(1.0))

        time.exposure = 1E20
        self.assertEqual(time.exposure, 1E20)

        time.resolution = float(2.0)
        self.assertEqual(time.resolution, float(2.0))

        time.timesys = "timesys"
        self.assertEqual(time.timesys, "timesys")

        time.trefpos = "trefpos"
        self.assertEqual(time.trefpos, "trefpos")

        time.mjdref = float(3.0)
        self.assertEqual(time.mjdref, float(3.0))


class TestPolarizationWCS(unittest.TestCase):
    def test_init(self):
        self.assertRaises(TypeError, chunk.PolarizationWCS, None)
        self.assertRaises(TypeError, chunk.PolarizationWCS, int(1))

        axis = wcs.Axis('STOKES')
        axis_1d = wcs.CoordAxis1D(axis)
        polarization = chunk.PolarizationWCS(axis_1d)
        self.assertEqual(polarization.axis, axis_1d)


class TestCustomWCS(unittest.TestCase):
    def test_init(self):
        self.assertRaises(TypeError, chunk.CustomWCS, None)
        self.assertRaises(TypeError, chunk.CustomWCS, int(1))
        self.assertRaises(TypeError, chunk.CustomWCS, str("s"))

        axis = wcs.Axis("ctype", "cunit")
        axis_1d = wcs.CoordAxis1D(axis)
        custom = chunk.CustomWCS(axis_1d)
        self.assertEqual(custom.axis, axis_1d)
        with self.assertRaises(AttributeError):
            custom.axis = axis_1d
