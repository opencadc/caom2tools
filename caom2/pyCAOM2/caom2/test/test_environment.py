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

""" Defines TestEnvironment class """

from caom2.caom2_environment import Environment
import os
import sys
import unittest

# put build at the start of the search path
sys.path.insert(0, os.path.abspath('../../lib.local/lib'))


class TestEnvironment(unittest.TestCase):

    def testAll(self):
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

suite = unittest.TestLoader().loadTestsFromTestCase(TestEnvironment)
unittest.TextTestRunner(verbosity=2).run(suite)
