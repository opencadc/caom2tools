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

import unittest
from datetime import datetime

from .. import artifact
from .. import chunk
from .. import observation
from .. import plane
from .. import shape
from .. import wcs


class TestEnums(unittest.TestCase):
    def test_all(self):
        # test for invalid value
        with self.assertRaises(ValueError):
            plane.CalibrationLevel("no_such_string")
        with self.assertRaises(ValueError):
            plane.CalibrationLevel(None)
        with self.assertRaises(ValueError):
            plane.CalibrationLevel(999)

        with self.assertRaises(ValueError):
            plane.DataProductType("no_such_string")
        with self.assertRaises(ValueError):
            plane.DataProductType(None)
        with self.assertRaises(ValueError):
            plane.DataProductType(1)

        with self.assertRaises(ValueError):
            plane.EnergyBand("no_such_string")
        with self.assertRaises(ValueError):
            plane.EnergyBand(None)
        with self.assertRaises(ValueError):
            plane.EnergyBand(1)

        with self.assertRaises(ValueError):
            plane.PolarizationState("no_such_string")
        with self.assertRaises(ValueError):
            plane.PolarizationState(None)
        with self.assertRaises(ValueError):
            plane.PolarizationState(1)

        with self.assertRaises(ValueError):
            plane.Quality("no_such_string")
        with self.assertRaises(ValueError):
            plane.Quality(None)
        with self.assertRaises(ValueError):
            plane.Quality(1)

        # test that we can get the object for each enum by name
        self.assertEqual(plane.CalibrationLevel.PLANNED.value, -1)
        self.assertEqual(plane.CalibrationLevel.RAW_INSTRUMENTAL.value, 0)
        self.assertEqual(plane.CalibrationLevel.RAW_STANDARD.value, 1)
        self.assertEqual(plane.CalibrationLevel.CALIBRATED.value, 2)
        self.assertEqual(plane.CalibrationLevel.PRODUCT.value, 3)

        self.assertEqual(plane.DataProductType.IMAGE.value, "image")
        self.assertEqual(plane.DataProductType.CUBE.value, "cube")
        self.assertEqual(plane.DataProductType.EVENTLIST.value, "eventlist")
        self.assertEqual(plane.DataProductType.SPECTRUM.value, "spectrum")
        self.assertEqual(plane.DataProductType.TIMESERIES.value, "timeseries")
        self.assertEqual(plane.DataProductType.VISIBILITY.value, "visibility")
        self.assertEqual(plane.DataProductType.MEASUREMENTS.value,
                         "measurements")
        self.assertEqual(
            plane.DataProductType.CATALOG.value,
            "http://www.opencadc.org/caom2/DataProductType#catalog")

        self.assertEqual(plane.EnergyBand['RADIO'].value, "Radio")
        self.assertEqual(plane.EnergyBand['MILLIMETER'].value, "Millimeter")
        self.assertEqual(plane.EnergyBand['INFRARED'].value, "Infrared")
        self.assertEqual(plane.EnergyBand['OPTICAL'].value, "Optical")
        self.assertEqual(plane.EnergyBand['UV'].value, "UV")
        self.assertEqual(plane.EnergyBand['EUV'].value, "EUV")
        self.assertEqual(plane.EnergyBand['XRAY'].value, "X-ray")
        self.assertEqual(plane.EnergyBand['GAMMARAY'].value, "Gamma-ray")

        self.assertEqual(plane.PolarizationState['I'].value, "I")
        self.assertEqual(plane.PolarizationState['Q'].value, "Q")
        self.assertEqual(plane.PolarizationState['U'].value, "U")
        self.assertEqual(plane.PolarizationState['V'].value, "V")
        self.assertEqual(plane.PolarizationState['LL'].value, "LL")
        self.assertEqual(plane.PolarizationState['LR'].value, "LR")
        self.assertEqual(plane.PolarizationState['RL'].value, "RL")
        self.assertEqual(plane.PolarizationState['RR'].value, "RR")
        self.assertEqual(plane.PolarizationState['XX'].value, "XX")
        self.assertEqual(plane.PolarizationState['XY'].value, "XY")
        self.assertEqual(plane.PolarizationState['YX'].value, "YX")
        self.assertEqual(plane.PolarizationState['YY'].value, "YY")

        self.assertEqual(plane.Quality['JUNK'].value, "junk")


class TestPlane(unittest.TestCase):
    def test_all(self):
        test_plane = plane.Plane("ProdID")
        self.assertEqual("ProdID", test_plane.product_id, "Product ID")
        self.assertEqual(None, test_plane.creator_id, "Creator ID")
        test_plane.creator_id = "ivo://cadc.nrc.ca/users?tester"
        self.assertEqual("ivo://cadc.nrc.ca/users?tester",
                         test_plane.creator_id, "Creator ID")
        self.assertEqual(0, len(test_plane.artifacts),
                         "Default number of artifacts")
        self.assertIsNone(test_plane.meta_release, "Default meta release date")
        date_now = datetime.now()
        test_plane.meta_release = date_now
        self.assertEqual(date_now, test_plane.meta_release,
                         "Metadata release date")
        self.assertIsNone(test_plane.data_release, "Default data release date")
        date_now = datetime.now()
        test_plane.data_release = date_now
        self.assertEqual(date_now, test_plane.data_release,
                         "Data release date")
        self.assertIsNone(test_plane.data_product_type,
                          "Default data product type")
        test_plane.data_product_type = plane.DataProductType.IMAGE
        self.assertEqual(plane.DataProductType.IMAGE,
                         test_plane.data_product_type,
                         "Data product type")
        plane.DataProductType.extend('http://www.myorg/std', 'mytype')
        test_plane.data_product_type = plane.DataProductType.MYTYPE
        self.assertEqual(plane.DataProductType.MYTYPE,
                         test_plane.data_product_type,
                         "Mytype product type")
        self.assertIsNone(test_plane.calibration_level,
                          "Default calibration level")
        test_plane.calibration_level = plane.CalibrationLevel.CALIBRATED
        self.assertEqual(plane.CalibrationLevel.CALIBRATED,
                         test_plane.calibration_level,
                         "plane.CalibrationLevel")
        self.assertIsNone(test_plane.quality,
                          "Default quality")
        quality = plane.DataQuality(plane.Quality.JUNK)
        test_plane.quality = quality
        self.assertEqual(quality,
                         test_plane.quality, "plane.Quality")
        self.assertIsNone(test_plane.provenance, "Default provenance")
        provenance = plane.Provenance("myProv")
        test_plane.provenance = provenance
        self.assertEqual("myProv", test_plane.provenance.name,
                         "Provenance - name")
        self.assertIsNone(test_plane.metrics, "Default metrics")
        metrics = plane.Metrics()
        test_plane.metrics = metrics
        self.assertEqual(metrics, test_plane.metrics, "Provenance - metrics")
        # self.assertIsNone(plane.observable, "Default observable")
        self.assertIsNone(test_plane.position, "Default position")
        self.assertIsNone(test_plane.energy, "Default energy")
        self.assertIsNone(test_plane.time, "Default time")
        self.assertIsNone(test_plane.polarization, "Default polarization")

        test_artifact1 = artifact.Artifact("caom:GEMINI/222/333",
                                           chunk.ProductType.SCIENCE,
                                           artifact.ReleaseType.DATA)
        test_plane.artifacts["caom:GEMINI/222/333"] = test_artifact1
        self.assertEqual(1, len(test_plane.artifacts), "Artifacts")
        self.assertTrue("caom:GEMINI/222/333" in test_plane.artifacts.keys())

        test_artifact2 = artifact.Artifact("caom:CFHT/55/66",
                                           chunk.ProductType.SCIENCE,
                                           artifact.ReleaseType.DATA)
        test_plane.artifacts["caom:CFHT/55/66"] = test_artifact2
        self.assertEqual(2, len(test_plane.artifacts), "Artifacts")
        self.assertTrue("caom:GEMINI/222/333" in test_plane.artifacts.keys())
        self.assertTrue("caom:CFHT/55/66" in test_plane.artifacts.keys())

        # try to append a duplicate artifact
        test_artifact3 = artifact.Artifact("caom:GEMINI/222/333",
                                           chunk.ProductType.SCIENCE,
                                           artifact.ReleaseType.DATA)
        test_plane.artifacts["caom:GEMINI/222/333"] = test_artifact3
        self.assertEqual(2, len(test_plane.artifacts), "Artifacts")
        self.assertTrue("caom:GEMINI/222/333" in test_plane.artifacts.keys())
        self.assertTrue("caom:CFHT/55/66" in test_plane.artifacts.keys())

        # Error cases
        exception = False
        try:
            test_plane = plane.Plane(None)
        except TypeError:
            exception = True
        self.assertTrue(exception, "Null argument in initialize")

        # exception = False
        # try:
        #     plane.compute_observable()
        # except TypeError:
        #     exception = True
        # self.assertTrue(exception,
        #                 "compute_observable implemented - Testing needed")

        # exception = False
        # try:
        #    plane.compute_position()
        # except TypeError:
        #     exception = True
        # self.assertTrue(exception,
        #                 "compute_position implemented - Testing needed")

        exception = False
        try:
            test_plane.compute_energy()
        except NotImplementedError:
            exception = True
        self.assertTrue(exception,
                        "compute_energy implemented - Testing needed")

        exception = False
        try:
            test_plane.compute_time()
        except NotImplementedError:
            exception = True
        self.assertTrue(exception, "compute_time implemented - Testing needed")

        exception = False
        try:
            test_plane.compute_polarization()
        except NotImplementedError:
            exception = True
        self.assertTrue(exception, "compute_polarization implemented"
                                   " - Testing needed")


class TestPlaneURI(unittest.TestCase):
    def test_all(self):
        plane_uri = plane.PlaneURI("caom:GEMINI/12345/3333")
        self.assertEqual("caom:GEMINI/12345/3333", plane_uri.uri,
                         "Plane URI")
        self.assertEqual("GEMINI", plane_uri.get_observation_uri().collection,
                         "Collection")
        self.assertEqual("12345",
                         plane_uri.get_observation_uri().observation_id,
                         "Observation ID")
        self.assertEqual("3333", plane_uri.get_product_id(), "Product ID")

        plane_uri = plane.PlaneURI.get_plane_uri(
            observation.ObservationURI("caom:CFHT/654321"),
            "555")
        self.assertEqual("caom:CFHT/654321/555", plane_uri.uri,
                         "Observation URI")
        self.assertEqual("CFHT", plane_uri.get_observation_uri().collection,
                         "Collection")
        self.assertEqual("654321",
                         plane_uri.get_observation_uri().observation_id,
                         "Observation ID")
        self.assertEqual("555", plane_uri.get_product_id(), "Product ID")

        exception = False
        try:
            plane_uri = plane.PlaneURI.get_plane_uri(None, "123")
        except TypeError:
            exception = True
        self.assertTrue(exception, "Missing exception")

        exception = False
        try:
            plane_uri = plane.PlaneURI.get_plane_uri("GEMINI", None)
        except TypeError:
            exception = True
        self.assertTrue(exception, "Missing exception")

        # wrong scheme
        exception = False
        try:
            plane_uri = plane.PlaneURI("somescheme:GEMINI/12345/3333")
        except ValueError:
            exception = True
        self.assertTrue(exception, "Missing exception")

        exception = False
        try:
            plane_uri = plane.PlaneURI("caom:GEMINI/12345")
        except ValueError:
            exception = True
        self.assertTrue(exception, "Missing exception")


class TestDataQuality(unittest.TestCase):
    def test_all(self):
        self.assertRaises(TypeError, plane.DataQuality, "string")
        quality = plane.DataQuality(plane.Quality.JUNK)
        self.assertEqual(plane.Quality.JUNK, quality.flag,
                         "DataQuality flag")


class TestMetrics(unittest.TestCase):
    def test_all(self):
        metrics = plane.Metrics()

        self.assertIsNone(metrics.source_number_density,
                          "Default source number density")
        metrics.source_number_density = 22.22
        self.assertEqual(22.22, metrics.source_number_density,
                         "Source number density")
        self.assertIsNone(metrics.background, "Default background")
        metrics.background = 12.34
        self.assertEqual(12.34, metrics.background, "Background")
        self.assertIsNone(metrics.background_std_dev,
                          "Default background standard deviation")
        metrics.background_std_dev = 34.34
        self.assertEqual(34.34, metrics.background_std_dev,
                         "Background standard deviation")
        self.assertIsNone(metrics.flux_density_limit,
                          "Default flux density limit")
        metrics.flux_density_limit = 55.55
        self.assertEqual(55.55, metrics.flux_density_limit,
                         "Flux density limit")
        self.assertIsNone(metrics.mag_limit, "Default mag limit")
        metrics.mag_limit = 20.08
        self.assertEqual(20.08, metrics.mag_limit, "Mag limit")


class TestProvenance(unittest.TestCase):
    def test_all(self):
        provenance = plane.Provenance("MyProvenance")
        self.assertEqual("MyProvenance", provenance.name, "Name")

        self.assertIsNone(provenance.version, "Default version")
        provenance.version = "XII"
        self.assertEqual("XII", provenance.version, "Version")
        self.assertIsNone(provenance.project, "Default project")
        provenance.project = "CFHTLS"
        self.assertEqual("CFHTLS", provenance.project, "Project")
        self.assertIsNone(provenance.producer, "Default producer")
        provenance.producer = "prod"
        self.assertEqual("prod", provenance.producer, "Producer")
        self.assertIsNone(provenance.run_id, "Default run ID")
        provenance.run_id = "A23"
        self.assertEqual("A23", provenance.run_id, "Run ID")
        self.assertIsNone(provenance.reference, "Default reference")

        self.assertEqual(0, len(provenance.inputs), "Default inputs")
        plane_uri1 = plane.PlaneURI("caom:HST/11/00")
        provenance.inputs.add(plane_uri1)
        self.assertEqual(1, len(provenance.inputs), "Default inputs")
        self.assertTrue(plane_uri1 in provenance.inputs)

        plane_uri2 = plane.PlaneURI("caom:HST/22/00")
        provenance.inputs.add(plane_uri2)
        self.assertEqual(2, len(provenance.inputs), "Default inputs")
        self.assertTrue(plane_uri1 in provenance.inputs)
        self.assertTrue(plane_uri2 in provenance.inputs)

        # testing duplicates
        plane_uri3 = plane.PlaneURI("caom:HST/22/00")
        provenance.inputs.add(plane_uri3)
        self.assertEqual(2, len(provenance.inputs), "Default inputs")
        self.assertTrue(plane_uri1 in provenance.inputs)
        self.assertTrue(plane_uri2 in provenance.inputs)

        self.assertIsNone(provenance.last_executed, "Default last executed")
        now_date = datetime.now()
        provenance.last_executed = now_date
        self.assertEqual(now_date, provenance.last_executed, "Last executed")

        self.assertEqual(0, len(provenance.keywords), "0 default keywords")
        provenance.keywords.add("keyword1")
        self.assertEqual(1, len(provenance.keywords), "1 keyword")
        self.assertTrue("keyword1" in provenance.keywords, "Keyword not found")

        provenance.keywords.add("keyword2")
        self.assertEqual(2, len(provenance.keywords), "2 keyword")
        self.assertTrue("keyword2" in provenance.keywords, "Keyword not found")

        # test the full constructor
        provenance = plane.Provenance("MyOtherProvenance",
                                      "Version2.0",
                                      "JCMT",
                                      "Mutt Lang",
                                      "b32",
                                      "caom:JCMT/33/00",
                                      now_date)

        self.assertIsNotNone(provenance.name)
        self.assertIsNotNone(provenance.version)
        self.assertIsNotNone(provenance.project)
        self.assertIsNotNone(provenance.producer)
        self.assertIsNotNone(provenance.run_id)
        self.assertIsNotNone(provenance.reference)
        self.assertIsNotNone(provenance.last_executed)

        self.assertEqual("MyOtherProvenance", provenance.name, "name")
        self.assertEqual("Version2.0", provenance.version, "version")
        self.assertEqual("JCMT", provenance.project, "project")
        self.assertEqual("Mutt Lang", provenance.producer, "producer")
        self.assertEqual("b32", provenance.run_id, "run_id")
        self.assertEqual("caom:JCMT/33/00", provenance.reference, "reference")
        self.assertEqual(now_date, provenance.last_executed, "last_executed")


class TestPosition(unittest.TestCase):
    def test_all(self):
        position = plane.Position()

        self.assertIsNone(position.bounds, "Default bounds")
        # position.bounds = 123
        # self.assertEqual(123, position.bounds, "Bounds")
        self.assertIsNone(position.dimension, "Default dimension")
        # position.dimension = 123
        # self.assertEqual(123, position.dimension, "Dimension")
        self.assertIsNone(position.resolution, "Default resolution")
        position.resolution = 123.321
        self.assertEqual(123.321, position.resolution, "Resolution")
        self.assertIsNone(position.sample_size, "Default sample size")
        position.sample_size = 321.123
        self.assertEqual(321.123, position.sample_size, "Sample size")
        self.assertFalse(position.time_dependent, "Default time dependent")
        position.time_dependent = True
        self.assertTrue(position.time_dependent, "Time dependent")


class TestEnergy(unittest.TestCase):
    def test_all(self):
        energy = plane.Energy()
        self.assertIsNone(energy.bounds, "Default energy bounds")
        energy.bounds = shape.Interval(1.0, 2.0)
        self.assertEqual(1.0, energy.bounds.lower, "Energy lower bounds")
        self.assertEqual(2.0, energy.bounds.upper, "Energy upper bounds")
        self.assertIsNone(energy.dimension, "Default energy dimension")
        energy.dimension = 1000
        self.assertEqual(1000, energy.dimension, "Energy dimension")
        self.assertIsNone(energy.resolving_power,
                          "Default energy resolving power")
        energy.resolving_power = 123.12
        self.assertEqual(123.12, energy.resolving_power,
                         "Energy resolving power")
        self.assertIsNone(energy.sample_size, "Default energy sample size")
        energy.sample_size = 123.321
        self.assertEqual(123.321, energy.sample_size, "Energy sample size")
        self.assertIsNone(energy.bandpass_name, "Default energy band pass")
        energy.bandpass_name = "EBN"
        self.assertEqual("EBN", energy.bandpass_name, "Energy bandpass name")
        self.assertIsNone(energy.em_band, "Default energy em band")
        self.assertTrue(0 == len(energy.energy_bands), "Default energy bands")
        energy.energy_bands.add(plane.EnergyBand.OPTICAL)
        self.assertEqual(1, len(energy.energy_bands), 'Energy bands')
        self.assertEqual(plane.EnergyBand.OPTICAL, energy.energy_bands.pop(),
                         "Energy band")
        self.assertIsNone(energy.transition, "Default energy transition")
        energy.transition = wcs.EnergyTransition("aSpecies", "aTransition")
        self.assertEqual("aSpecies", energy.transition.species,
                         "Energy transition species")
        self.assertEqual("aTransition", energy.transition.transition,
                         "Energy transition transition")


class TestEnergyTransition(unittest.TestCase):
    def test__init__(self):
        # test for invalid values
        self.assertRaises(TypeError, wcs.EnergyTransition, None, None)
        self.assertRaises(TypeError, wcs.EnergyTransition, 'aString', None)
        self.assertRaises(TypeError, wcs.EnergyTransition, None, 'aString')
        self.assertRaises(TypeError, wcs.EnergyTransition, 1, 'aString')
        self.assertRaises(TypeError, wcs.EnergyTransition, 'aString', 2)
        # test for happy path
        transition = wcs.EnergyTransition("aSpecies", "aTransition")
        self.assertEqual(transition._species, "aSpecies")
        self.assertEqual(transition._transition, "aTransition")

    def test_setters(self):
        # test that we cannot change the attribute values
        transition = wcs.EnergyTransition("aSpecies", "aTransition")
        try:
            transition.species = "newSpecies"
            transition.transition = "newTransition"
        except AttributeError:
            pass
        else:
            raise AttributeError("at least one attribute was changed")


class TestPolarizaton(unittest.TestCase):
    def test_all(self):
        polarization = plane.Polarization()

        self.assertIsNone(polarization.dimension,
                          "Default polarization dimension")

        # TODO add test for state


class TestTime(unittest.TestCase):
    def test_all(self):
        time = plane.Time()
        self.assertIsNone(time.bounds, "Default bounds")
        self.assertIsNone(time.dimension, "Default dimension")
        self.assertIsNone(time.resolution, "Default resolution")
        self.assertIsNone(time.sample_size, "Default sample size")
        self.assertIsNone(time.exposure, "Default exposure")

        time.dimension = 777
        self.assertEqual(777, time.dimension, "Dimension")
        time.resolution = 77.777
        self.assertEqual(77.777, time.resolution, "Resolution")
        time.sample_size = 12.34
        self.assertEqual(12.34, time.sample_size, "Sample size")
        time.exposure = 55.55
        self.assertEqual(55.55, time.exposure, "Exposure")


class TestCustomAxis(unittest.TestCase):
    def test_all(self):
        with self.assertRaises(AttributeError):
            plane.CustomAxis(None)
        my_axis = plane.CustomAxis('Foo')
        self.assertEqual('Foo', my_axis.ctype, 'CTYPE missmatch')
        self.assertIsNone(my_axis.bounds, "Default bounds")
        self.assertIsNone(my_axis.dimension, "Default dimension")

        my_axis.dimension = 777
        my_axis.bounds = shape.Interval(1.0, 2.0)
        self.assertEqual(777, my_axis.dimension, "Dimension")
        self.assertEqual(1.0, my_axis.bounds.lower, "Bounds mismatch")
        self.assertEqual(2.0, my_axis.bounds.upper, "Bounad mismatch")

        my_axis = plane.CustomAxis('Blah', bounds=shape.Interval(3.0, 4.0),
                                   dimension=33)
        self.assertEqual('Blah', my_axis.ctype, 'CTYPE missmatch')
        self.assertEqual(33, my_axis.dimension, 'Dimension missmatch')
        self.assertEqual(3.0, my_axis.bounds.lower, "Bounds mismatch")
        self.assertEqual(4.0, my_axis.bounds.upper, "Bounad mismatch")
