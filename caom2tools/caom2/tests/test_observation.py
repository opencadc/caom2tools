#!/usr/bin/env python2.7
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

""" Defines TestObservation class """

import unittest
from datetime import datetime

from caom2.shape import Point
from caom2.observation import Algorithm
from caom2.observation import CompositeObservation
from caom2.observation import Environment
from caom2.observation import Instrument
from caom2.observation import Observation
from caom2.observation import ObservationIntentType
from caom2.observation import ObservationURI
from caom2.observation import Proposal
from caom2.observation import Requirements
from caom2.observation import SimpleObservation
from caom2.observation import Status
from caom2.observation import Target
from caom2.observation import TargetPosition
from caom2.observation import TargetType
from caom2.observation import Telescope
from caom2.plane import Plane
from caom2.util import Util


class TestEnums(unittest.TestCase):

    def test_all(self):
        # test for invalid value
        self.assertEqual(ObservationIntentType.get("no_such_string"), None)
        self.assertRaises(AttributeError, ObservationIntentType.get, None)
        self.assertRaises(AttributeError, ObservationIntentType.get, 1)

        self.assertEqual(Status.get("no_such_string"), None)
        self.assertRaises(AttributeError, Status.get, None)
        self.assertRaises(AttributeError, Status.get, 1)

        self.assertEqual(TargetType.get("no_such_string"), None)
        self.assertRaises(AttributeError, TargetType.get, None)
        self.assertRaises(AttributeError, TargetType.get, 1)

        # test that we can get the object for each enum by name
        self.assertEqual(ObservationIntentType.CALIBRATION.value, "calibration")
        self.assertEqual(ObservationIntentType.SCIENCE.value, "science")

        self.assertEqual(Status.FAIL.value, "fail")

        self.assertEqual(TargetType.FIELD.value, "field")
        self.assertEqual(TargetType.OBJECT.value, "object")


class TestObservation(unittest.TestCase):

    def test_all(self):
        algorithm = Algorithm("myAlg")
        obs = Observation("GSA", "A12345", algorithm)
        self.assertEqual("GSA", obs.collection, "Collection")
        self.assertEqual("A12345", obs.observation_id, "Observation ID")
        self.assertEqual(algorithm, obs.algorithm, "Algorithm")

        new_algorithm = Algorithm("myNewAlg")
        obs.algorithm = new_algorithm
        self.assertEquals(new_algorithm, obs.algorithm, "New algorithm")

        self.assertIsNone(obs.intent, "Default intent")
        obs.intent = ObservationIntentType.CALIBRATION
        self.assertEqual(ObservationIntentType.CALIBRATION,
                         obs.intent, "Observation intent")

        self.assertIsNone(obs.obs_type, "Default obs_type")
        obs.obs_type = "obstype1"
        self.assertEqual("obstype1",
                         obs.obs_type, "obs type")

        self.assertIsNone(obs.proposal, "Default proposal")
        proposal = Proposal("ABC")
        obs.proposal = proposal
        self.assertEqual(proposal,
                         obs.proposal, "Proposal")

        self.assertIsNone(obs.telescope, "Default telescope")
        telescope = Telescope("GSAGN")
        obs.telescope = telescope
        self.assertEqual(telescope, obs.telescope, "Telescope")

        self.assertIsNone(obs.instrument, "Default instrument")
        instrument = Instrument("NIRI")
        obs.instrument = instrument
        self.assertEqual(instrument, obs.instrument, "Instrument")

        self.assertIsNone(obs.target, "Default target")
        target = Target("TGT")
        obs.target = target
        self.assertEqual(target, obs.target, "Target")

        self.assertIsNone(obs.target_position, "Default target position")
        target_position = TargetPosition(Point(1.0, 2.0), "coordsys")
        obs.target_position = target_position
        self.assertEqual(target_position,
                         obs.target_position, "TargetPosition")

        self.assertIsNone(obs.requirements, "Default requirements")
        requirements = Requirements(Status.FAIL)
        obs.requirements = requirements
        self.assertEqual(requirements,
                         obs.requirements, "Requirements")

        self.assertIsNone(obs.environment, "Default environment")
        environment = Environment()
        obs.environment = environment
        self.assertEqual(environment,
                         obs.environment, "Environment")

        self.assertIsNone(obs.meta_release, "Default metadata release")
        date_now = datetime.now()
        obs.meta_release = date_now
        self.assertEqual(date_now,
                         obs.meta_release, "Metadata release")

        self.assertEqual(0, len(obs.planes), "Default planes")
        plane1 = Plane("myPlaneID")
        obs.planes["myPlaneID"] = plane1
        self.assertEqual(1, len(obs.planes), "Planes")
        self.assertTrue("myPlaneID" in obs.planes.keys())

        plane2 = Plane("myPlaneID2")
        obs.planes["myPlaneID2"] = plane2
        self.assertEqual(2, len(obs.planes), "Planes")
        self.assertTrue("myPlaneID" in obs.planes)
        self.assertTrue("myPlaneID2" in obs.planes.keys())

        # test duplicates
        plane3 = Plane("myPlaneID2")
        obs.planes["myPlaneID2"] = plane3
        self.assertEqual(2, len(obs.planes), "Planes")
        self.assertTrue("myPlaneID" in obs.planes)
        self.assertTrue("myPlaneID2" in obs.planes.keys())

        obs2 = Observation(obs.collection,
                           obs.observation_id,
                           obs.algorithm,
                           planes=obs.planes,
                           sequence_number=obs.sequence_number,
                           intent=obs.intent,
                           obs_type=obs.obs_type,
                           proposal=obs.proposal,
                           telescope=obs.telescope,
                           instrument=obs.instrument,
                           target=obs.target,
                           meta_release=obs.meta_release,
                           environment=obs.environment,
                           target_position=obs.target_position)


class TestObservationURI(unittest.TestCase):

    def test_all(self):
        obsURI = ObservationURI("caom:GEMINI/12345")
        self.assertEqual("caom:GEMINI/12345", obsURI.uri, "Observation URI")
        self.assertEqual("GEMINI", obsURI.collection, "Collection")
        self.assertEqual("12345", obsURI.observation_id, "Observation ID")

        obsURI = ObservationURI.get_observation_uri("CFHT", "654321")
        self.assertEqual("caom:CFHT/654321", obsURI.uri, "Observation URI")
        self.assertEqual("CFHT", obsURI.collection, "Collection")
        self.assertEqual("654321", obsURI.observation_id, "Observation ID")

        exception = False
        try:
            obsURI = ObservationURI.get_observation_uri(None, "123")
        except TypeError:
            exception = True
        self.assertTrue(exception, "Missing exception")

        exception = False
        try:
            obsURI = ObservationURI.get_observation_uri("GEMINI", None)
        except TypeError:
            exception = True
        self.assertTrue(exception, "Missing exception")


class TestSimpleObservation(unittest.TestCase):

    def test_all(self):
        algorithm = SimpleObservation._ALGORITHM
        obs = SimpleObservation("GSA", "A12345")
        self.assertEqual("GSA", obs.collection, "Collection")
        self.assertEqual("A12345", obs.observation_id, "Observation ID")

        self.assertEqual(algorithm, obs.algorithm, "Algorithm")
        obs.algorithm = algorithm
        self.assertEqual(algorithm, obs.algorithm, "Algorithm")

        # try to set algorithm
        exception = False
        try:
            obs.algorithm = Algorithm("myAlg")
        except ValueError:
            exception = True
        self.assertTrue(exception, "Missing exception")

        # run the rest of the Observation tests
        self.assertIsNone(obs.intent, "Default intent")
        obs.intent = ObservationIntentType.CALIBRATION
        self.assertEqual(ObservationIntentType.CALIBRATION,
                         obs.intent, "Observation intent")

        self.assertIsNone(obs.obs_type, "Default obs_type")
        obs.obs_type = "obstype1"
        self.assertEqual("obstype1",
                         obs.obs_type, "obs type")

        self.assertIsNone(obs.proposal, "Default proposal")
        proposal = Proposal("ABC")
        obs.proposal = proposal
        self.assertEqual(proposal,
                         obs.proposal, "Proposal")

        self.assertIsNone(obs.telescope, "Default telescope")
        telescope = Telescope("GSAGN")
        obs.telescope = telescope
        self.assertEqual(telescope,
                         obs.telescope, "Telescope")

        self.assertIsNone(obs.instrument, "Default instrument")
        instrument = Instrument("NIRI")
        obs.instrument = instrument
        self.assertEqual(instrument,
                         obs.instrument, "Instrument")

        self.assertIsNone(obs.target, "Default target")
        target = Target("TGT")
        obs.target = target
        self.assertEqual(target,
                         obs.target, "Target")

        self.assertIsNone(obs.environment, "Default environment")
        environment = Environment()
        obs.environment = environment
        self.assertEqual(environment,
                         obs.environment, "Environment")

        self.assertIsNone(obs.target_position, "Default target position")
        target_position = TargetPosition(Point(1.0, 2.0), "coordsys")
        obs.target_position = target_position
        self.assertEqual(target_position,
                         obs.target_position, "TargetPosition")

        self.assertIsNone(obs.requirements, "Default requirements")
        requirements = Requirements(Status.FAIL)
        obs.requirements = requirements
        self.assertEquals(requirements, obs.requirements, "requirements")

        self.assertIsNone(obs.meta_release, "Default metadata release")
        date_now = datetime.now()
        obs.meta_release = date_now
        self.assertEqual(date_now,
                         obs.meta_release, "Metadata release")

    # Test the complete constructor
    def test_complete_init(self):
        collection = str("CFHT")
        observationID = str("543210")
        algorithm = SimpleObservation._ALGORITHM
        sequence_number = int(3)
        intent = ObservationIntentType.SCIENCE
        obs_type = str("foo")
        proposal = Proposal("123")
        telescope = Telescope("TEL")
        instrument = Instrument("INST")
        target = Target("LMC")
        meta_release = datetime.now()
        planes = Util.TypedOrderedDict((Plane),)
        environment = Environment()

        obs = SimpleObservation(collection,
                                observationID,
                                algorithm,
                                sequence_number,
                                intent,
                                obs_type,
                                proposal,
                                telescope,
                                instrument,
                                target,
                                meta_release,
                                planes,
                                environment)

        self.assertIsNotNone(obs.collection, "Collection")
        self.assertEqual(collection, obs.collection, "Collection")

        self.assertIsNotNone(obs.observation_id, "Observation ID")
        self.assertEqual(observationID, obs.observation_id, "Observation ID")

        self.assertIsNotNone(obs.algorithm, "Algorithm")
        self.assertEqual(algorithm, obs.algorithm, "Algorithm")

        self.assertIsNotNone(obs.intent, "Observation intent")
        self.assertEqual(intent, obs.intent, "Observation intent")

        self.assertIsNotNone(obs.obs_type, "obs type")
        self.assertEqual(obs_type, obs.obs_type, "obs type")

        self.assertIsNotNone(obs.proposal, "Proposal")
        self.assertEqual(proposal, obs.proposal, "Proposal")

        self.assertIsNotNone(obs.telescope, "Telescope")
        self.assertEqual(telescope, obs.telescope, "Telescope")

        self.assertIsNotNone(obs.instrument, "Instrument")
        self.assertEqual(instrument, obs.instrument, "Instrument")

        self.assertIsNotNone(obs.target, "Target")
        self.assertEqual(target, obs.target, "Target")

        self.assertIsNotNone(obs.meta_release, "Metadata release")
        self.assertEqual(meta_release, obs.meta_release, "Metadata release")

        self.assertIsNotNone(obs.planes, "Planes")
        self.assertEqual(planes, obs.planes, "Planes")

        self.assertIsNotNone(obs.environment, "Environment")
        self.assertEqual(environment, obs.environment, "Environment")


class TestCompositeObservation(unittest.TestCase):

    def test_all(self):
        algorithm = Algorithm("mozaic")
        obs = CompositeObservation("GSA", "A12345", algorithm)
        self.assertEqual("GSA", obs.collection, "Collection")
        self.assertEqual("A12345", obs.observation_id, "Observation ID")
        self.assertEqual(algorithm, obs.algorithm, "Algorithm")
        obs.algorithm = algorithm
        self.assertEqual(algorithm, obs.algorithm, "Algorithm")

        # try to set algorithm to an invalid value
        exception = False
        try:
            obs.algorithm = SimpleObservation._ALGORITHM
        except ValueError:
            exception = True
        self.assertTrue(exception, "Missing exception")

        # try to set algorithm to None
        exception = False
        try:
            obs.algorithm = None
        except ValueError:
            exception = True
        self.assertTrue(exception, "Missing exception")

        self.assertEqual(0, len(obs.members), "Members")
        observationURI1 = ObservationURI("caom:collection/obsID")
        obs.members.add(observationURI1)
        self.assertEqual(1, len(obs.members), "Members")
        self.assertTrue(observationURI1 in obs.members)

        observationURI2 = ObservationURI("caom:collection/obsID2")
        obs.members.add(observationURI2)
        self.assertEqual(2, len(obs.members), "Members")
        self.assertTrue(observationURI1 in obs.members)
        self.assertTrue(observationURI2 in obs.members)

        #duplicates
        observationURI3 = ObservationURI("caom:collection/obsID")
        obs.members.add(observationURI3)
        self.assertEqual(2, len(obs.members), "Members")
        self.assertTrue(observationURI1 in obs.members)
        self.assertTrue(observationURI2 in obs.members)

        # run the rest of the Observation tests
        self.assertIsNone(obs.intent, "Default intent")
        obs.intent = ObservationIntentType.CALIBRATION
        self.assertEqual(ObservationIntentType.CALIBRATION,
                         obs.intent, "Observation intent")

        self.assertIsNone(obs.obs_type, "Default obs_type")
        obs.obs_type = "obstype1"
        self.assertEqual("obstype1",
                         obs.obs_type, "obs type")

        self.assertIsNone(obs.proposal, "Default proposal")
        proposal = Proposal("ABC")
        obs.proposal = proposal
        self.assertEqual(proposal,
                         obs.proposal, "Proposal")

        self.assertIsNone(obs.telescope, "Default telescope")
        telescope = Telescope("GSAGN")
        obs.telescope = telescope
        self.assertEqual(telescope,
                         obs.telescope, "Telescope")

        self.assertIsNone(obs.instrument, "Default instrument")
        instrument = Instrument("NIRI")
        obs.instrument = instrument
        self.assertEqual(instrument,
                         obs.instrument, "Instrument")

        self.assertIsNone(obs.target, "Default target")
        target = Target("TGT")
        obs.target = target
        self.assertEqual(target,
                         obs.target, "Target")

        self.assertIsNone(obs.environment, "Default environment")
        environment = Environment()
        obs.environment = environment
        self.assertEqual(environment,
                         obs.environment, "Environment")

        self.assertIsNone(obs.target_position, "Default target position")
        target_position = TargetPosition(Point(1.0, 2.0), "coordsys")
        obs.target_position = target_position
        self.assertEqual(target_position,
                         obs.target_position, "TargetPosition")

        self.assertIsNone(obs.requirements, "Default requirements")
        requirements = Requirements(Status.FAIL)
        obs.requirements = requirements
        self.assertEquals(requirements, obs.requirements, "requirements")

        self.assertIsNone(obs.meta_release, "Default metadata release")
        date_now = datetime.now()
        obs.meta_release = date_now
        self.assertEqual(date_now,
                         obs.meta_release, "Metadata release")

    # Test the complete constructor
    def test_complete_init(self):
        collection = str("CFHT")
        observationID = str("543210")
        algorithm = str("algo")
        sequence_number = int(3)
        intent = ObservationIntentType.SCIENCE
        obs_type = str("foo")
        proposal = Proposal("123")
        telescope = Telescope("TEL")
        instrument = Instrument("INST")
        target = Target("LMC")
        meta_release = datetime.now()
        planes = Util.TypedOrderedDict((Plane),)
        environment = Environment()
        target_position = TargetPosition(Point(1.0, 2.0), "coordsys")

        obs = CompositeObservation(collection,
                                   observationID,
                                   algorithm,
                                   sequence_number,
                                   intent,
                                   obs_type,
                                   proposal,
                                   telescope,
                                   instrument,
                                   target,
                                   meta_release,
                                   planes,
                                   environment,
                                   target_position)

        self.assertIsNotNone(obs.collection, "Collection")
        self.assertEqual(collection, obs.collection, "Collection")

        self.assertIsNotNone(obs.observation_id, "Observation ID")
        self.assertEqual(observationID, obs.observation_id, "Observation ID")

        self.assertIsNotNone(obs.algorithm, "Algorithm")
        self.assertEqual(algorithm, obs.algorithm, "Algorithm")

        self.assertIsNotNone(obs.intent, "Observation intent")
        self.assertEqual(intent, obs.intent, "Observation intent")

        self.assertIsNotNone(obs.obs_type, "obs type")
        self.assertEqual(obs_type, obs.obs_type, "obs type")

        self.assertIsNotNone(obs.proposal, "Proposal")
        self.assertEqual(proposal, obs.proposal, "Proposal")

        self.assertIsNotNone(obs.telescope, "Telescope")
        self.assertEqual(telescope, obs.telescope, "Telescope")

        self.assertIsNotNone(obs.instrument, "Instrument")
        self.assertEqual(instrument, obs.instrument, "Instrument")

        self.assertIsNotNone(obs.target, "Target")
        self.assertEqual(target, obs.target, "Target")

        self.assertIsNotNone(obs.meta_release, "Metadata release")
        self.assertEqual(meta_release, obs.meta_release, "Metadata release")

        self.assertIsNotNone(obs.planes, "Planes")
        self.assertEqual(planes, obs.planes, "Planes")

        self.assertIsNotNone(obs.environment, "Environment")
        self.assertEqual(environment, obs.environment, "Environment")

        self.assertIsNotNone(obs.target_position, "TargetPosition")
        self.assertEqual(target_position, obs.target_position,
                         "TargetPosition")

        # Try changing the algorithm
        algorithm2 = str("new algo")
        obs.algorithm = algorithm2
        self.assertIsNotNone(obs.algorithm, "Algorithm")
        self.assertNotEqual(algorithm, obs.algorithm, "Algorithm")
        self.assertEqual(algorithm2, obs.algorithm, "Algorithm")


class TestAlgorithm(unittest.TestCase):

    def test_all(self):
        algorithm = Algorithm("myAlgorithm")
        self.assertEqual("myAlgorithm", algorithm.name, "Algorithm name")


class TestEnvironment(unittest.TestCase):

    def test_all(self):
        environment = Environment()

        self.assertIsNone(environment.seeing, "Default seeing")
        environment.seeing = 123.321
        self.assertEqual(123.321, environment.seeing, "Seeing")
        self.assertIsNone(environment.humidity, "Default humidity")
        environment.humidity = 0.333
        self.assertEqual(0.333, environment.humidity, "humidity")
        self.assertIsNone(environment.elevation, "Default elevation")
        environment.elevation = 12.12
        self.assertEqual(12.12, environment.elevation, "Elevation")
        self.assertIsNone(environment.tau, "Default tau")
        environment.tau = 0.456
        self.assertEqual(0.456, environment.tau, "Tau")
        self.assertIsNone(environment.wavelength_tau, "Default wavelength tau")
        environment.wavelength_tau = 200.02
        self.assertEqual(200.02, environment.wavelength_tau, "Wavelength tau")
        self.assertIsNone(environment.ambient_temp,
                          "Default ambient temperature")
        environment.ambient_temp = 12.44
        self.assertEqual(12.44, environment.ambient_temp,
                         "Ambient temperature")
        self.assertIsNone(environment.photometric, "Default photometric")
        environment.photometric = True
        self.assertTrue(environment.photometric, "Photometric")


class TestIntrument(unittest.TestCase):

    def test_all(self):
        instrument = Instrument("myInstrument")
        self.assertEqual("myInstrument", instrument.name, "Instrument name")
        self.assertEqual(0, len(instrument.keywords), "Default number of keywords")

        instrument.keywords.add("optical")
        self.assertEqual(1, len(instrument.keywords), "Number of keywords")
        self.assertTrue("optical" in instrument.keywords, "Keyword not found")

        instrument.keywords.add("radio")
        self.assertEqual(2, len(instrument.keywords), "Number of keywords")
        self.assertTrue("radio" in instrument.keywords, "Keyword not found")


class TestProposal(unittest.TestCase):

    def test_all(self):
        proposal = Proposal("myProposal")
        self.assertEqual("myProposal", proposal.proposal_id, "Proposal ID")
        self.assertEqual(0, len(proposal.keywords), "Default number of keywords")
        proposal.keywords.add("optical")
        self.assertEqual(1, len(proposal.keywords), "Number of keywords")
        self.assertTrue("optical" in proposal.keywords, "Keyword not found")
        self.assertIsNone(proposal.pi_name, "Default PI")
        proposal.pi_name = "John Doe"
        self.assertEqual("John Doe", proposal.pi_name, "PI")
        self.assertIsNone(proposal.project, "Default PI")
        proposal.project = "Project A"
        self.assertEqual("Project A", proposal.project, "Project")
        self.assertIsNone(proposal.title, "Default title")
        proposal.title = "Something Interesting"
        self.assertEqual("Something Interesting", proposal.title, "Title")


class TestRequirements(unittest.TestCase):

    def test_all(self):
        self.assertRaises(TypeError, Requirements, "string")
        requirements = Requirements(Status.FAIL)
        self.assertEqual(Status.FAIL, requirements.flag,
                         "Requirements flag")


class TestTarget(unittest.TestCase):

    def test_all(self):
        target = Target("myTarget")
        self.assertEqual("myTarget", target.name, "target name")

        target.target_type = TargetType.FIELD
        self.assertEqual(TargetType.FIELD, target.target_type, "target type")

        self.assertEqual(0, len(target.keywords), "Default number of keywords")
        target.keywords.add("optical")
        self.assertEqual(1, len(target.keywords), "Number of keywords")
        self.assertTrue("optical" in target.keywords, "Keyword not found")

        self.assertIsNone(target.redshift, "Default redshift")
        target.redshift = 123.321
        self.assertEqual(123.321, target.redshift, "Redshift")

        self.assertIsNone(target.standard, "Default standard")
        target.standard = True
        self.assertTrue(target.standard, "Standard")

        self.assertIsNone(target.moving, "Default moving")
        target.moving = True
        self.assertTrue(target.moving, "Moving")

        target = Target("myOtherTarget", TargetType.OBJECT, False, 1.2, {"radio"}, False)
        self.assertEquals("myOtherTarget", target.name, "target name")
        self.assertEquals(TargetType.OBJECT, target.target_type, "target type")
        self.assertFalse(target.standard, "Standard")
        self.assertEquals(1.2, target.redshift, "Redshift")
        self.assertEquals(1, len(target.keywords), "Keywords")
        self.assertTrue("radio" in target.keywords, "Keywords")
        self.assertFalse(target.moving, "Moving")


class TestTargetPosition(unittest.TestCase):

    def test_all(self):
        self.assertRaises(TypeError, TargetPosition, "string")
        point = Point(1.0, 2.0)
        target_position = TargetPosition(point, "coordsys")
        self.assertIsNotNone(target_position.coordinates,
                             "target position coordinates")
        self.assertEqual(point.cval1, target_position.coordinates.cval1,
                         "coordinates cval1")
        self.assertEqual(point.cval2, target_position.coordinates.cval2,
                         "coordinates cval2")
        self.assertIsNotNone(target_position.coordsys,
                             "target position coordsys")
        self.assertEqual("coordsys", target_position.coordsys, "coordsys")
        self.assertIsNone(target_position.equinox,
                          "target position equinox")

        target_position = TargetPosition(point, "coordsys", 1.0)
        self.assertIsNotNone(target_position.equinox,
                             "target position equinox")
        self.assertEqual(1.0, target_position.equinox,
                         "equinox")


class TestTelescope(unittest.TestCase):

    def test_all(self):
        telescope = Telescope("myTelescope")
        self.assertEqual("myTelescope", telescope.name, "telescope name")
        self.assertEqual(0, len(telescope.keywords), "Default number of keywords")

        telescope.keywords.add("optical")
        self.assertEqual(1, len(telescope.keywords), "Number of keywords")
        self.assertTrue("optical" in telescope.keywords, "Keyword not found")

        self.assertIsNone(telescope.geo_location_x, "Default geo location x")
        telescope.geo_location_x = 123.321
        self.assertEqual(123.321, telescope.geo_location_x, "Geo location x")

        self.assertIsNone(telescope.geo_location_y, "Default geo location y")
        telescope.geo_location_y = 333.33
        self.assertEqual(333.33, telescope.geo_location_y, "Geo location y")

        self.assertIsNone(telescope.geo_location_z, "Default geo location z")
        telescope.geo_location_z = 12.12
        self.assertEqual(12.12, telescope.geo_location_z, "Geo location z")



if __name__ == '__main__':
    unittest.main()
