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
import copy
import sys
import os
from mock import Mock, patch, MagicMock
#TODO to be changed to io.StringIO when caom2 is prepared for python3
from StringIO import StringIO
from datetime import datetime

from caom2visitor import CAOM2Visitor
from caom2visitor.visitor import DATE_FORMAT
from caom2.caom2_simple_observation import SimpleObservation
from caom2.xml.caom2_observation_writer import ObservationWriter
from caom2.xml.caom2_observation_reader import ObservationReader

THIS_DIR = os.path.dirname(os.path.realpath(__file__))

class TestCAOM2Visitor(unittest.TestCase):

    """Test the Caom2Visitor class"""


    def test_plugin_class(self):
        
        
        
        # plugin class does not change the observation
        visitor = CAOM2Visitor(os.path.join(
                THIS_DIR, 'passplugin.py'), 'cfht')
        obs = SimpleObservation('cfht', '7000000o')
        expect_obs = copy.deepcopy(obs)
        visitor.plugin.update(obs)
        self.assertEquals(expect_obs, obs)
        
        # plugin class adds a plane to the observation
        visitor = CAOM2Visitor(os.path.join(
                THIS_DIR, 'addplaneplugin.py'), 'cfht')
        obs = SimpleObservation('cfht', '7000000o')
        expect_obs = copy.deepcopy(obs)
        visitor.plugin.update(obs)
        self.assertNotEquals(expect_obs, obs)
        self.assertEquals(len(expect_obs.planes) + 1, len(obs.planes))
        
        # non-existent the plugin file
        with self.assertRaises(Exception):
            visitor = CAOM2Visitor('blah.py', 'cfht')
        
        # non-existent ObservationUpdater class in the plugin file
        with self.assertRaises(Exception):
            visitor = CAOM2Visitor(os.path.join(
                THIS_DIR, 'test_visitor.py'), 'cfht')
            
        # non-existent update method in ObservationUpdater class
        with self.assertRaises(Exception):    
            visitor = CAOM2Visitor(os.path.join(
                THIS_DIR, 'noupdateplugin.py'), 'cfht')
   
   
    # patch sleep to stop the test from sleeping and slowing down execution
    @patch('caom2repoClient.caom2repoClient.time.sleep', MagicMock(), create=True)
    @patch('caom2repoClient.caom2repoClient.open', MagicMock(), create=True)
    @patch('caom2repoClient.caom2repoClient.HTTPSConnection')
    def test_get_observation(self, mock_conn):
                
        obs = SimpleObservation('cfht', '7000000o')
        writer = ObservationWriter()
        ibuffer = StringIO()
        writer.write(obs, ibuffer)
        response = MagicMock()
        response.status = 200
        response.read.return_value = ibuffer.getvalue()
        conn = mock_conn.return_value
        conn.getresponse.return_value = response
        ibuffer.seek(0) # reposition the buffer for reading
        
        visitor = CAOM2Visitor(os.path.join(
                THIS_DIR, 'passplugin.py'), 'cfht')
        
        self.assertEquals(
            'ivo://cadc.nrc.ca/caom2repo/cfht', visitor.collection_uri)
        
        self.assertEquals(obs, visitor._get_observation('700000o'))
        conn.request.assert_called_once_with(
            'GET', '/caom2repo/pub/cfht/700000o', '',
            {})
        
        # signal some problems
        response.status = 404
        with self.assertRaises(IOError):
            visitor._get_observation('700000o')
            
        response.status = 405
        with self.assertRaises(IOError):
            visitor._get_observation('700000o')
        
        # transient problem no retries
        response.status = 503
        with self.assertRaises(IOError):
            visitor._get_observation('700000o')
            
        # permanent transient problem with 2 retries
        response.status = 503
        visitor = CAOM2Visitor(os.path.join(
                THIS_DIR, 'passplugin.py'), 'cfht', retries=2)
        with self.assertRaises(IOError):
            visitor._get_observation('700000o')
            
        # temporary transient problem. Succeeds after 2 retries
        response.status = 200
        ibuffer.seek(0) # reposition the buffer for reading
        transient_response = MagicMock()
        transient_response.status = 503
        conn = MagicMock()
        conn.getresponse.side_effect = [transient_response, transient_response, response]
        mock_conn.return_value = conn 
        visitor = CAOM2Visitor(os.path.join(
                THIS_DIR, 'passplugin.py'), 'cfht', retries=2)
        self.assertEquals(obs, visitor._get_observation('700000o'))


    # patch sleep to stop the test from sleeping and slowing down execution
    @patch('caom2repoClient.caom2repoClient.time.sleep', MagicMock(), create=True)
    @patch('caom2repoClient.caom2repoClient.open', MagicMock(), create=True)
    @patch('caom2visitor.visitor.CAOM2RepoClient')
    def test_get_observations(self, mock_repo):
        # This is almost similar to the previous test except that it gets
        # observations matching a collection and start/end criteria
        # Also, patch the CAOM2RepoClient now.

        response = MagicMock()
        response.status = 200
        last_datetime = '2000-10-10T12:30:00.333'
        response.read.return_value = '700000o,2000-10-10T12:20:11.123\n700001o,' +\
            last_datetime
        repo = mock_repo.return_value
        repo.send_request.return_value = response
        
        visitor = CAOM2Visitor(os.path.join(
                THIS_DIR, 'passplugin.py'), 'cfht')
        
        self.assertEquals(None, visitor.start)
        self.assertEquals(None, visitor.end)
        self.assertEquals(None, visitor.current_start)
        end_date = datetime.strptime(last_datetime, DATE_FORMAT)
        
        expect_observations = ['700000o', '700001o']
        self.assertEquals(expect_observations, visitor._get_observations())
        self.assertEquals(end_date, visitor.current_start)
        repo.send_request.assert_called_once_with('GET', 
            '/cfht',
            {"Content-type": "application/x-www-form-urlencoded"}, 
            'MAXREC=10000')
        
        # check the generated urls for different parameters
        repo.reset_mock()
        #set the size of a batch
        os.environ['CAOM2_VISITOR_BATCH_SIZE'] = '100'
        visitor = CAOM2Visitor(os.path.join(
                THIS_DIR, 'passplugin.py'), 'cfht', 
                start=datetime.strptime('2000-11-11', '%Y-%m-%d'))
        visitor._get_observations()
        repo.send_request.assert_called_once_with('GET', 
            '/cfht',
            {"Content-type": "application/x-www-form-urlencoded"}, 
            'START=2000-11-11T00%3A00%3A00.000000&MAXREC=100')
        del os.environ['CAOM2_VISITOR_BATCH_SIZE']
        
        repo.reset_mock()
        visitor = CAOM2Visitor(os.path.join(
                THIS_DIR, 'passplugin.py'), 'cfht', 
                end=datetime.strptime('2000-11-11', '%Y-%m-%d'))
        visitor._get_observations()
        repo.send_request.assert_called_once_with('GET', 
            '/cfht',
            {"Content-type": "application/x-www-form-urlencoded"}, 
            'END=2000-11-11T00%3A00%3A00.000000&MAXREC=10000')
        
        repo.reset_mock()
        visitor = CAOM2Visitor(os.path.join(
                THIS_DIR, 'passplugin.py'), 'cfht', 
                start=datetime.strptime('2000-11-11', '%Y-%m-%d'),
                end=datetime.strptime('2000-11-12', '%Y-%m-%d'))
        visitor._get_observations()
        repo.send_request.assert_called_once_with('GET', 
            '/cfht',
            {"Content-type": "application/x-www-form-urlencoded"}, 
            'START=2000-11-11T00%3A00%3A00.000000&' +\
            'END=2000-11-12T00%3A00%3A00.000000&MAXREC=10000')
        
        # signal some problems
        response.status = 404
        with self.assertRaises(IOError):
            visitor._get_observations()
            
        response.status = 405
        with self.assertRaises(IOError):
            visitor._get_observations()
        
        # transient problem no retries
        response.status = 503
        with self.assertRaises(IOError):
            visitor._get_observations()
            
        # permanent transient problem with 2 retries
        response.status = 503
        visitor = CAOM2Visitor(os.path.join(
                THIS_DIR, 'passplugin.py'), 'cfht', retries=2)
        repo.retry.return_value = response
        with self.assertRaises(IOError):
            visitor._get_observations()
            
        # temporary transient problem. Succeeds after 2 retries
        repo.reset_mock()
        response.status = 200
        transient_response = MagicMock()
        transient_response.status = 503
        repo.send_request.return_value = transient_response
        repo.retry.return_value = response
        visitor = CAOM2Visitor(os.path.join(
                THIS_DIR, 'passplugin.py'), 'cfht', retries=2)
        self.assertEquals(expect_observations, visitor._get_observations())
        self.assertEquals(end_date, visitor.current_start)
        repo.send_request.assert_called_once_with('GET', 
            '/cfht',
            {"Content-type": "application/x-www-form-urlencoded"}, 
            'MAXREC=10000')
        repo.retry.assert_called_once_with('GET', 
            '/cfht',
            {"Content-type": "application/x-www-form-urlencoded"}, 
            'MAXREC=10000')
                
    
    # patch sleep to stop the test from sleeping and slowing down execution
    @patch('caom2repoClient.caom2repoClient.time.sleep', MagicMock(), create=True)
    @patch('caom2repoClient.caom2repoClient.open', MagicMock(), create=True)
    @patch('caom2repoClient.caom2repoClient.HTTPSConnection')
    def test_persist_observation(self, mock_conn):
        obs = SimpleObservation('cfht', '7000000o')
        visitor = CAOM2Visitor(os.path.join(
                THIS_DIR, 'passplugin.py'), 'cfht')
        response = MagicMock()
        response.status = 200
        conn = mock_conn.return_value
        conn.getresponse.return_value = response
        #mock_conn.return_value = conn        
        iobuffer = StringIO()
        ObservationWriter().write(obs, iobuffer)
        obsxml = iobuffer.getvalue()
        
        visitor._persist_observation('700000o', obs)
        
        conn.request.assert_called_once_with(
            'POST', '/caom2repo/pub/cfht/700000o', obsxml,
            {'Content-Type': 'text/xml'})
        
        conn.reset_mock()
        
        # signal some problems
        response.status = 404
        with self.assertRaises(IOError):
            visitor._persist_observation('700000o', obs)
            
        response.status = 405
        with self.assertRaises(IOError):
            visitor._persist_observation('700000o', obs)
        
        # transient problem no retries
        response.status = 503
        with self.assertRaises(IOError):
            visitor._persist_observation('700000o', obs)
            
        # permanent transient problem with 2 retries
        response.status = 503
        visitor = CAOM2Visitor(os.path.join(
                THIS_DIR, 'passplugin.py'), 'cfht', retries=2)
        with self.assertRaises(IOError):
            visitor._persist_observation('700000o', obs)
            
        # temporary transient problem. Succeeds after 2 retries
        response.status = 200
        transient_response = MagicMock()
        transient_response.status = 503
        conn = MagicMock()
        conn.getresponse.side_effect = [transient_response, transient_response, response]
        mock_conn.return_value = conn 
        visitor = CAOM2Visitor(os.path.join(
                THIS_DIR, 'passplugin.py'), 'cfht', retries=2)
        visitor._persist_observation('700000o', obs) 
        conn.request.assert_called_with(
            'POST', '/caom2repo/pub/cfht/700000o', obsxml,
            {'Content-Type': 'text/xml'})
        

    @patch('test_visitor.CAOM2Visitor._get_observation', MagicMock(), create = True)
    @patch('test_visitor.CAOM2Visitor._persist_observation', MagicMock(), create = True)
    @patch('test_visitor.CAOM2Visitor._get_observations')
    def test_process(self, mock_obs):
        os.environ['CAOM2_VISITOR_BATCH_SIZE'] = '3' # size of the batch is 3
        mock_obs.side_effect = [['a', 'b', 'c'], ['d'], []]
        visitor = CAOM2Visitor(os.path.join(
                THIS_DIR, 'passplugin.py'), 'cfht')
        visitor.plugin = MagicMock()
        self.assertEquals(4, visitor.process())

        mock_obs.side_effect = [['a', 'b', 'c'], ['d', 'e', 'f'], []]
        visitor = CAOM2Visitor(os.path.join(
                THIS_DIR, 'passplugin.py'), 'cfht')
        visitor.plugin = MagicMock()
        self.assertEquals(6, visitor.process())
