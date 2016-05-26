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

""" Defines TestSimpleObservation class """

import unittest
from datetime import datetime

from caom2.caom2_algorithm import Algorithm
from caom2.caom2_enums import ObservationIntentType
from caom2.caom2_enums import Status
from caom2.caom2_environment import Environment
from caom2.caom2_instrument import Instrument
from caom2.caom2_plane import Plane
from caom2.caom2_proposal import Proposal
from caom2.caom2_requirements import Requirements
from caom2.caom2_simple_observation import SimpleObservation
from caom2.caom2_target import Target
from caom2.caom2_target_position import TargetPosition
from caom2.caom2_telescope import Telescope
from caom2.types.caom2_point import Point
from caom2.util.caom2_util import TypedOrderedDict


class TestSimpleObservation(unittest.TestCase):

    def testAll(self):
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
    def testCompleteInit(self):
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
        planes = TypedOrderedDict((Plane),)
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


if __name__ == '__main__':
    unittest.main()

