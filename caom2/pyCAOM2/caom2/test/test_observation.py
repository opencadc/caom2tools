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

from caom2.caom2_observation import Observation
from caom2.caom2_algorithm import Algorithm
from caom2.caom2_proposal import Proposal
from caom2.caom2_telescope import Telescope
from caom2.caom2_instrument import Instrument
from caom2.caom2_target import Target
from caom2.caom2_environment import Environment
from caom2.caom2_plane import Plane
from caom2.caom2_enums import ObservationIntentType
import os
import sys
import unittest

from datetime import datetime

# put build at the start of the search path
sys.path.insert(0, os.path.abspath('../../lib.local/lib'))


class TestObservation(unittest.TestCase):

    def testAll(self):
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
                           environment=obs.environment)


suite = unittest.TestLoader().loadTestsFromTestCase(TestObservation)
unittest.TextTestRunner(verbosity=2).run(suite)
