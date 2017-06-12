# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2017.                            (c) 2017.
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

""" Defines TestCaom2IdGenerator class """

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import unittest
import binascii
import os
import sys

from cadcutils.net import auth
from cadcutils.util import utils
from caom2repo import CAOM2RepoClient
from io import BytesIO
from caom2 import common
from caom2 import chunk
from caom2 import part
from caom2 import plane
from caom2 import artifact
from caom2 import observation
from caom2 import ObservationReader
from datetime import datetime

THIS_DIR = os.path.dirname(os.path.realpath(__file__))

class TestCaom2Integration(unittest.TestCase):
            
    def test_create_and_visit(self):
        
        start = datetime.now()
        name = obs_name = 'caom2pyinttest{}'.format(start.microsecond)
        
        try:
            env_a = os.environ['A']
            cert_file = env_a + '/test-certificates/x509_CADCAuthtest1.pem'
            subject = auth.Subject(certificate=cert_file)
            client = CAOM2RepoClient(subject)
            
            # create one observation for today
            algorithm = observation.SimpleObservation._ALGORITHM
            
            print("test obs name {}".format(name))
            obs = observation.SimpleObservation("TEST", obs_name)
            obs.algorithm = algorithm
            client.put_observation(obs)
            
            plugin = os.path.join(THIS_DIR, 'visitor-plugin.py')
            (visited, updated, skipped, failed) = client.visit(plugin, 'TEST', start=start, halt_on_error=True)
            print("observations visited: {}".format(len(visited)))
            
            self.assertGreater(visited, 0, msg="No Observations Visited")
        
        finally:
            try:
                client.delete_observation("TEST", name)
            except:
                print('Failed to delete test observation, continuing')
        

if __name__ == '__main__':
    unittest.main()
