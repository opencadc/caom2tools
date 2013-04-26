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

""" Defines TestCaom2Enums class """

from caom2.caom2_enums import EnergyBand, CalibrationLevel
import os
import sys
import unittest


# put build at the start of the search path
sys.path.insert(0, os.path.abspath('../../lib.local/lib'))


class TestCaom2Enums(unittest.TestCase):

    # EnergyBand and CalibrationLevel are used to test the open source
    # Enum package that the Caom2Enums use. Caom2Enums only use a subset
    # of the features provided by the Enum package. Only these subset
    # features are tested.
    def testToOBject(self):
        # test for invalid value
        self.assertEqual(EnergyBand.get("no_such_string"), None)
        self.assertEqual(EnergyBand.get("g"), None)
        self.assertEqual(EnergyBand.get("no_such_string"), None)
        self.assertRaises(AttributeError, EnergyBand.get, None)
        self.assertRaises(AttributeError, EnergyBand.get, 1)
        # test that we can get the object for each enum by name
        self.assertEqual(CalibrationLevel.CALIBRATED.name, "CALIBRATED",
                         "incorrect enum name")
        self.assertEqual(CalibrationLevel.get(CalibrationLevel.CALIBRATED.name)
                         .name, "CALIBRATED")
        self.assertEqual(CalibrationLevel.get('CALIBRATED').value, 2)
        self.assertEqual(CalibrationLevel.get(CalibrationLevel.CALIBRATED.name)
                         .value, 2)
        self.assertEqual(CalibrationLevel.
                         getByValue(CalibrationLevel.CALIBRATED.value)
                         .value, 2)
        self.assertEqual(CalibrationLevel.
                         getByValue(CalibrationLevel.CALIBRATED.value).name,
                         "CALIBRATED")
        self.assertEqual(EnergyBand.get('RADIO').value, "Radio")
        self.assertEqual(EnergyBand.get('MILLIMETER').value, "Millimeter")
        self.assertEqual(EnergyBand.get('INFRARED').value, "Infrared")
        self.assertEqual(EnergyBand.get('OPTICAL').value, "Optical")
        self.assertEqual(EnergyBand.get('UV').value, "UV")
        self.assertEqual(EnergyBand.get('EUV').value, "EUV")
        self.assertEqual(EnergyBand.get('XRAY').value, "X-ray")
        self.assertEqual(EnergyBand.get('GAMMARAY').value, "Gamma-ray")
        # test that we can get the object for each enum by value
        self.assertEqual(EnergyBand.getByValue(EnergyBand.RADIO.value).value,
                         "Radio")
        self.assertEqual(EnergyBand.getByValue(EnergyBand.MILLIMETER.value)
                         .value, "Millimeter")
        self.assertEqual(EnergyBand.getByValue(EnergyBand.INFRARED.value)
                         .value, "Infrared")
        self.assertEqual(EnergyBand.getByValue(EnergyBand.OPTICAL.value)
                         .value, "Optical")
        self.assertEqual(EnergyBand.getByValue(EnergyBand.UV.value).value,
                         "UV")
        self.assertEqual(EnergyBand.getByValue(EnergyBand.EUV.value).value,
                         "EUV")
        self.assertEqual(EnergyBand.getByValue(EnergyBand.XRAY.value).value,
                         "X-ray")
        self.assertEqual(EnergyBand.getByValue(EnergyBand.GAMMARAY.value)
                         .value, "Gamma-ray")

    def testGetValue(self):
        # test that we can get each enum value
        self.assertEqual(EnergyBand.RADIO.value, "Radio")
        self.assertEqual(EnergyBand.MILLIMETER.value, "Millimeter")
        self.assertEqual(EnergyBand.INFRARED.value, "Infrared")
        self.assertEqual(EnergyBand.OPTICAL.value, "Optical")
        self.assertEqual(EnergyBand.UV.value, "UV")
        self.assertEqual(EnergyBand.EUV.value, "EUV")
        self.assertEqual(EnergyBand.XRAY.value, "X-ray")
        self.assertEqual(EnergyBand.GAMMARAY.value, "Gamma-ray")

    def testSetValue(self):
        # test that we cannot change an enum value
        self.assertRaises(TypeError, EnergyBand.RADIO.value, "InvalidValue")
        self.assertRaises(TypeError, EnergyBand.RADIO.name, "InvalidName")

suite = unittest.TestLoader().loadTestsFromTestCase(TestCaom2Enums)
unittest.TextTestRunner(verbosity=2).run(suite)
