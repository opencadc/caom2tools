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

""" Defines TestRoundTrip class """

from caom2.xml.caom2_observation_reader import ObservationReader
from caom2.xml.caom2_observation_writer import ObservationWriter
import filecmp
import errno
import glob
import os
import shutil
import sys
import unittest

# put build at the start of the search path
sys.path.insert(0, os.path.abspath('../../lib.local/lib'))


class TestRoundTrip(unittest.TestCase):

    TEST_FILE_SOURCE_DIR = '../caom2/test/data'
    XML_FILE_SOURCE_DIR = '/tmp/caom2-round-trip-test'
    XML_FILE_DEST_DIR = '/tmp/caom2-round-trip-test/pyCAOM2'

    def makeTestDir(self):
        try:
            os.makedirs(TestRoundTrip.XML_FILE_SOURCE_DIR)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    def copyFiles(self):
        for filename in glob.glob(os.path.join(TestRoundTrip.TEST_FILE_SOURCE_DIR,
                                               '*.xml')):
            shutil.copy(filename, TestRoundTrip.XML_FILE_SOURCE_DIR)

    def init(self):
        self.makeTestDir()
        self.copyFiles()

    def cleanup(self):
        shutil.rmtree(TestRoundTrip.XML_FILE_SOURCE_DIR)

    def get_file_list(self):
        return [f for f in os.listdir(TestRoundTrip.XML_FILE_SOURCE_DIR)\
                if f.endswith('.xml')]

    # This test reads each file in XML_FILE_SOURCE_DIR, creates the CAOM2
    # objects and writes a file in XML_FILE_DEST_DIR based on the CAOM2
    # objects. The two XML files are then compared to ensure that they
    # are the same. The test fails if the files are not the same. The
    # test/data/*.xml files can be used in this test.

    def test_round_trip(self):
        print "Test Round Trip"

        try:
            self.init()
            files = self.get_file_list()
            self.assertTrue(len(files) > 0, 'No XML files in ' +
                            TestRoundTrip.XML_FILE_SOURCE_DIR)
            try:
                os.makedirs(TestRoundTrip.XML_FILE_DEST_DIR)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise

            reader = ObservationReader(True)
            writer = ObservationWriter(True, False)
            for filename in files:
                sourceFilePath = (TestRoundTrip.XML_FILE_SOURCE_DIR +
                                  '/' + filename)
                print "test file: " + sourceFilePath
                sourceXMLfp = open(sourceFilePath, 'r')
                observation = reader.read(sourceFilePath)
                sourceXMLfp.close()
                destFilePath = TestRoundTrip.XML_FILE_DEST_DIR + '/' + filename
                destXMLfp = open(destFilePath, 'w')
                writer.write(observation, destXMLfp)
                destXMLfp.close()
                self.assertTrue(filecmp.cmp(sourceFilePath, destFilePath),
                                'files are different, ' +\
                                'file from Java was: ' + sourceFilePath + ' '\
                                'file from Python was: ' + destFilePath)
        #finally:
        #    self.cleanup()
        except Exception as e:
            #if e.errno != errno.EEXIST:
            raise

suite = unittest.TestLoader().loadTestsFromTestCase(TestRoundTrip)
unittest.TextTestRunner(verbosity=2).run(suite)
