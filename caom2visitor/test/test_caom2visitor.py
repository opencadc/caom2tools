## -*- coding: utf-8 -*-
#***********************************************************************
#******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
#*************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#                                                                                                                                                          
#  (c) 2016.                            (c) 2016.                                                                                                          
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

import unittest
import imp
import sys
from datetime import datetime
from os import path
from StringIO import StringIO
from mock import Mock, patch, MagicMock
from caom2visitor import CAOM2Visitor


THIS_DIR = path.dirname(path.realpath(__file__))

class TestCAOM2Visitor(unittest.TestCase):

    """Test the Caom2Visitor class"""

    @patch('sys.stdout', new_callable=StringIO)
    @patch('caom2visitor.CAOM2Visitor')
    def test_application(self, mock_visitor, mock_print):
        # invokes the caom2visitor application with different arguments
        # and tests the arguments that are being passed to the constructor
        # of CAOM2Visitor class
        # Because caom2visitor is not part of a package, the code
        # is dynamically loaded and invoked directly from the file.
        scripts_dir = path.join(path.dirname(THIS_DIR), 'scripts')
        
        plugin_file = path.join(path.dirname(THIS_DIR), 
                                'caom2visitor', 'test', 'passplugin.py')
        sys.argv = ['caom2visitor', '--plugin', plugin_file, 'ad']
        # dynamically load the source code
        py_mod = imp.load_source('caom2visitor', path.join(scripts_dir, 'caom2visitor'))
        self.assertTrue(hasattr(py_mod, 'run'))
        # this is how to call the run function
        getattr(py_mod, 'run')()
        mock_visitor.assert_called_with(
            plugin=plugin_file, start=None, end=None, collection='ad', retries=None, server=None)
        
        
        # test with start time argument        
        mock_visitor.reset_mock()
        start_date = '2001-01-01'
        sys.argv = ['caom2visitor', '--plugin', plugin_file, '--start', start_date, 'ad']
        getattr(py_mod, 'run')()
        mock_visitor.assert_called_with(
            plugin=plugin_file, start=datetime.strptime(start_date, '%Y-%m-%d'), 
            end=None, collection='ad', retries=None, server=None)
        
        
        # test with end time argument        
        mock_visitor.reset_mock()
        end_date = '2011-01-01'
        sys.argv = ['caom2visitor', '--plugin', plugin_file, '--end', end_date, 'ad']
        getattr(py_mod, 'run')()
        mock_visitor.assert_called_with(
            plugin=plugin_file, end=datetime.strptime(end_date, '%Y-%m-%d'), 
            start=None, collection='ad', retries=None, server=None)
        
        
        # test with start and end time arguments      
        mock_visitor.reset_mock()
        sys.argv = ['caom2visitor', '--plugin', plugin_file, '--end', end_date, 
                    '--start', start_date, 'ad']
        getattr(py_mod, 'run')()
        mock_visitor.assert_called_with(
            plugin=plugin_file, end=datetime.strptime(end_date, '%Y-%m-%d'), 
            start=datetime.strptime(start_date, '%Y-%m-%d'), collection='ad', 
            retries=None, server=None)
        

        # test with retries
        mock_visitor.reset_mock()
        sys.argv = ['caom2visitor', '--plugin', plugin_file, '--end', end_date, 
                    '--start', start_date, '--retries' , '7', '--server', 'www.test.net', 'ad']
        getattr(py_mod, 'run')()
        mock_visitor.assert_called_with(
            plugin=plugin_file, end=datetime.strptime(end_date, '%Y-%m-%d'), 
            start=datetime.strptime(start_date, '%Y-%m-%d'), collection='ad', 
            retries=7, server='www.test.net')
        
        
        # test the help option 
        mock_visitor.reset_mock()
        mock_print.seek(0)
        sys.argv = ['caom2visitor', '--help']
        with self.assertRaises(SystemExit):
            getattr(py_mod, 'run')()
            
            
        expected_out=\
"""usage: caom2visitor [-h] [--cert <CertFile>] [--token <TokenString>]
                    [--version] [-d] [-v] [-w] --plugin <pluginClassFile>
                    [--start <datetime start point>]
                    [--end <datetime end point>]
                    [--retries <number of retries>] [-s <CAOM2 service URL>]
                    <datacollection>

    Visitor over a data collection in CAOM2 observation repository. The
    code provided in the plugin is for updating each visited observation.

positional arguments:
  <datacollection>      data collection in CAOM2 repo

optional arguments:
  -h, --help            show this help message and exit
  --cert <CertFile>     location of your CADC security certificate file
  --token <TokenString>
                        token string (alternative to certfile)
  --version             show program's version number and exit
  -d, --debug           Print debug level log messages
  -v, --verbose         Print verbose level log messages
  -w, --warning         Print warning level log messages
  --plugin <pluginClassFile>
                        Pluging class to update each observation
  --start <datetime start point>
                        oldest dataset to visit (UTC %Y-%m-%d format)
  --end <datetime end point>
                        earliest dataset to visit (UTC %Y-%m-%d format)
  --retries <number of retries>
                        number of tries with transient server errors
  -s <CAOM2 service URL>, --server <CAOM2 service URL>
                        URL of the CAOM2 repo server

Environment:
       CAOM2_VISITOR_BATCH_SIZE: the size of a get batch

Minimum plugin file format:
----
   from caom2.caom2_observation import Observation

   class ObservationUpdater:

    def update(self, observation):
        assert isinstance(observation, Observation), (
            'observation {} is not an Observation'.format(observation))
        # custom code to update the observation
----
"""
        self.assertEquals(expected_out, mock_print.getvalue())
        
        
        
        # test incorrect pluging file
        mock_visitor.reset_mock()
        mock_print.seek(0)
        sys.argv = ['caom2visitor', '--plugin', 'blah', 'add']
        with self.assertRaises(SystemExit):
            getattr(py_mod, 'run')()
            
            
            
def run():
    suite1 = unittest.TestLoader().loadTestsFromTestCase(TestCAOM2Visitor)
    all_tests = unittest.TestSuite([suite1])
    return unittest.TextTestRunner(verbosity=2).run(all_tests)


if __name__ == "__main__":
    run()
