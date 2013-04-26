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

""" Defines TestPlane class """


from caom2.caom2_plane import Plane
from caom2.caom2_artifact import Artifact
from caom2.caom2_provenance import Provenance
from caom2.caom2_metrics import Metrics
from caom2.caom2_enums import DataProductType
from caom2.caom2_enums import CalibrationLevel
import os
import sys
import unittest
from datetime import datetime

# put build at the start of the search path
sys.path.insert(0, os.path.abspath('../../lib.local/lib'))


class TestPlane(unittest.TestCase):

    def testAll(self):
        plane = Plane("ProdID")
        self.assertEqual("ProdID", plane.product_id, "Product ID")
        self.assertEqual(0, len(plane.artifacts),
                         "Default number of artifacts")
        self.assertIsNone(plane.meta_release, "Default meta release date")
        date_now = datetime.now()
        plane.meta_release = date_now
        self.assertEqual(date_now, plane.meta_release, "Metadata release date")
        self.assertIsNone(plane.data_release, "Default data release date")
        date_now = datetime.now()
        plane.data_release = date_now
        self.assertEqual(date_now, plane.data_release, "Data release date")
        self.assertIsNone(plane.data_product_type, "Default data product type")
        plane.data_product_type = DataProductType.IMAGE
        self.assertEqual(DataProductType.IMAGE, plane.data_product_type,
                         "Data product type")
        self.assertIsNone(plane.calibration_level,
                          "Default calibration level")
        plane.calibration_level = CalibrationLevel.CALIBRATED
        self.assertEqual(CalibrationLevel.CALIBRATED,
                         plane.calibration_level, "CalibrationLevel")
        self.assertIsNone(plane.provenance, "Default provenance")
        provenance = Provenance("myProv")
        plane.provenance = provenance
        self.assertEqual("myProv", plane.provenance.name, "Provenance - name")
        self.assertIsNone(plane.metrics, "Default metrics")
        metrics = Metrics()
        plane.metrics = metrics
        self.assertEqual(metrics, plane.metrics, "Provenance - metrics")
        #self.assertIsNone(plane.observable, "Default observable")
        self.assertIsNone(plane.position, "Default position")
        self.assertIsNone(plane.energy, "Default energy")
        self.assertIsNone(plane.time, "Default time")
        self.assertIsNone(plane.polarization, "Default polarization")

        artifact1 = Artifact("caom:GEMINI/222/333")
        plane.artifacts["caom:GEMINI/222/333"] = artifact1
        self.assertEquals(1, len(plane.artifacts), "Artifacts")
        self.assertTrue("caom:GEMINI/222/333" in plane.artifacts.keys())

        artifact2 = Artifact("caom:CFHT/55/66")
        plane.artifacts["caom:CFHT/55/66"] = artifact2
        self.assertEquals(2, len(plane.artifacts), "Artifacts")
        self.assertTrue("caom:GEMINI/222/333" in plane.artifacts.keys())
        self.assertTrue("caom:CFHT/55/66" in plane.artifacts.keys())

        #try to append a duplicate artifact
        artifact3 = Artifact("caom:GEMINI/222/333")
        plane.artifacts["caom:GEMINI/222/333"] = artifact3
        self.assertEquals(2, len(plane.artifacts), "Artifacts")
        self.assertTrue("caom:GEMINI/222/333" in plane.artifacts.keys())
        self.assertTrue("caom:CFHT/55/66" in plane.artifacts.keys())

        #Error cases
        exception = False
        try:
            plane = Plane(None)
        except TypeError:
            exception = True
        self.assertTrue(exception, "Null argument in initialize")

        #exception = False
        #try:
        #    plane.compute_observable()
        #except TypeError:
        #    exception = True
        #self.assertTrue(exception,
        #                "compute_observable implemented - Testing needed")

        #exception = False
        #try:
        #    plane.compute_position()
        #except TypeError:
        #    exception = True
        #self.assertTrue(exception,
        #                "compute_position implemented - Testing needed")

        exception = False
        try:
            plane.compute_energy()
        except NotImplementedError:
            exception = True
        self.assertTrue(exception,
                        "compute_energy implemented - Testing needed")

        exception = False
        try:
            plane.compute_time()
        except NotImplementedError:
            exception = True
        self.assertTrue(exception, "compute_time implemented - Testing needed")

        exception = False
        try:
            plane.compute_polarization()
        except NotImplementedError:
            exception = True
        self.assertTrue(exception, "compute_polarization implemented"
                                    " - Testing needed")

suite = unittest.TestLoader().loadTestsFromTestCase(TestPlane)
unittest.TextTestRunner(verbosity=2).run(suite)
