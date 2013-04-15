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

""" Defines TestSpectralWCS class """


from caom2.wcs.caom2_spectral_wcs import SpectralWCS
from caom2.caom2_energy_transition import EnergyTransition
from caom2.wcs.caom2_axis import Axis
from caom2.wcs.caom2_coord_axis1d import CoordAxis1D
import os.path
import sys
import unittest

# put build at the start of the search path
sys.path.insert(0, os.path.abspath('../../lib.local/lib'))


class TestSpectralWCS(unittest.TestCase):

    def testInit(self):

        axis = Axis("ctype", "cunit")
        axis_1d = CoordAxis1D(axis)

        self.assertRaises(TypeError, SpectralWCS, None, None)
        self.assertRaises(TypeError, SpectralWCS, None, str("s"))
        self.assertRaises(TypeError, SpectralWCS, axis_1d, None)
        self.assertRaises(TypeError, SpectralWCS, int(1), str("s"))
        self.assertRaises(TypeError, SpectralWCS, axis_1d, int(1))

        energy = SpectralWCS(axis_1d, "specsys")
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
            energy.bandpassName = int(1)
            energy.transition = int(1)
            energy.resolvingPower = int(1)

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

        energy.bandpassName = "bandpassName"
        self.assertEqual(energy.bandpassName, "bandpassName")

        transition = EnergyTransition("species", "transition")
        energy.transition = transition
        self.assertEqual(energy.transition, transition)

        energy.resolvingPower = float(6.0)
        self.assertEqual(energy.resolvingPower, float(6.0))

suite = unittest.TestLoader().loadTestsFromTestCase(TestSpectralWCS)
unittest.TextTestRunner(verbosity=2).run(suite)
