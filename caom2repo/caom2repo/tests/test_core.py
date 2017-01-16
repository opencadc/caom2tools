# # -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
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
# ***********************************************************************
#
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import copy
import os
import sys
import unittest
# TODO to be changed to io.StringIO when caom2 is prepared for python3
from StringIO import StringIO
from datetime import datetime

import requests
from cadcutils import util, exceptions
from cadcutils.net import auth
from caom2.obs_reader_writer import ObservationWriter
from caom2.observation import SimpleObservation
from mock import Mock, patch, MagicMock, ANY

from caom2repo import core
from caom2repo.core import CAOM2RepoClient, DATE_FORMAT

THIS_DIR = os.path.dirname(os.path.realpath(__file__))


class MyExitError(Exception):
    pass


class TestCAOM2Repo(unittest.TestCase):

    """Test the Caom2Visitor class"""

    def test_plugin_class(self):
        # plugin class does not change the observation
        collection = 'cfht'
        observation_id = '7000000o'
        visitor = CAOM2RepoClient(auth.Subject())
        obs = SimpleObservation(collection, observation_id)
        expect_obs = copy.deepcopy(obs)
        visitor._load_plugin_class(os.path.join(THIS_DIR, 'passplugin.py'))
        visitor.plugin.update(obs)
        self.assertEquals(expect_obs, obs)
        
        # plugin class adds a plane to the observation
        visitor = CAOM2RepoClient(auth.Subject())
        obs = SimpleObservation('cfht', '7000000o')
        expect_obs = copy.deepcopy(obs)
        visitor._load_plugin_class(os.path.join(THIS_DIR, 'addplaneplugin.py'))
        visitor.plugin.update(obs)
        self.assertNotEquals(expect_obs, obs)
        self.assertEquals(len(expect_obs.planes) + 1, len(obs.planes))
        
        # non-existent the plugin file
        with self.assertRaises(Exception):
            visitor._load_plugin_class(os.path.join(THIS_DIR, 'blah.py'))
        
        # non-existent ObservationUpdater class in the plugin file
        with self.assertRaises(Exception):
            visitor._load_plugin_class(os.path.join(THIS_DIR, 'test_visitor.py'))
            
        # non-existent update method in ObservationUpdater class
        with self.assertRaises(Exception):    
            visitor._load_plugin_class(os.path.join(THIS_DIR, 'noupdateplugin.py'))

    # patch sleep to stop the test from sleeping and slowing down execution
    @patch('cadcutils.net.ws.time.sleep', MagicMock(), create=True)
    @patch('cadcutils.net.ws.open', MagicMock(), create=True)
    @patch('cadcutils.net.ws.Session.send')
    def test_get_observation(self, mock_get):
        collection = 'cfht'
        observation_id = '7000000o'
        service_url = 'www.cadc.nrc.ca/caom2repo'
        obs = SimpleObservation(collection, observation_id)
        writer = ObservationWriter()
        ibuffer = StringIO()
        writer.write(obs, ibuffer)
        response = MagicMock()
        response.status_code = 200
        response.content = ibuffer.getvalue()
        mock_get.return_value = response
        ibuffer.seek(0)  # reposition the buffer for reading
        visitor = CAOM2RepoClient(auth.Subject(), host=service_url)
        self.assertEquals(obs, visitor.get_observation(collection, observation_id))
        
        # signal problems
        http_error = requests.HTTPError()
        response.status_code = 500
        http_error.response = response
        response.raise_for_status.side_effect = [http_error]
        with self.assertRaises(exceptions.InternalServerException):
            visitor.get_observation(collection, observation_id)
            
        # temporary transient errors
        http_error = requests.HTTPError()
        response.status_code = 503
        http_error.response = response
        response.raise_for_status.side_effect = [http_error, None]
        visitor.get_observation(collection, observation_id)

        # permanent transient errors
        http_error = requests.HTTPError()
        response.status_code = 503
        http_error.response = response

        def raise_error(): raise http_error
        response.raise_for_status.side_effect = raise_error
        with self.assertRaises(exceptions.HttpException):
            visitor.get_observation(collection, observation_id)

    # patch sleep to stop the test from sleeping and slowing down execution
    @patch('cadcutils.net.ws.time.sleep', MagicMock(), create=True)
    @patch('cadcutils.net.ws.open', MagicMock(), create=True)
    @patch('caom2repo.core.net.BaseWsClient.get')
    def test_get_observations(self, mock_get):
        # This is almost similar to the previous test except that it gets
        # observations matching a collection and start/end criteria
        # Also, patch the CAOM2RepoClient now.

        response = MagicMock()
        response.status_code = 200
        last_datetime = '2000-10-10T12:30:00.333'
        response.content = '700000o,2000-10-10T12:20:11.123\n700001o,' +\
            last_datetime
        mock_get.return_value = response
        
        visitor = CAOM2RepoClient(auth.Subject())
        end_date = datetime.strptime(last_datetime, DATE_FORMAT)
        
        expect_observations = ['700000o', '700001o']
        self.assertEquals(expect_observations, visitor._get_observations('cfht'))
        self.assertEquals(end_date, visitor._start)
        mock_get.assert_called_once_with('cfht', params={'MAXREC': core.BATCH_SIZE})

        mock_get.reset_mock()
        visitor._get_observations('cfht', end=datetime.strptime('2000-11-11', '%Y-%m-%d'))
        mock_get.assert_called_once_with('cfht', params={'END': '2000-11-11T00:00:00.000000',
                                                         'MAXREC': core.BATCH_SIZE})

        mock_get.reset_mock()
        visitor._get_observations('cfht',
                                  start=datetime.strptime('2000-11-11', '%Y-%m-%d'),
                                  end=datetime.strptime('2000-11-12', '%Y-%m-%d'))
        mock_get.assert_called_once_with('cfht', params={'START': '2000-11-11T00:00:00.000000',
                                                         'END': '2000-11-12T00:00:00.000000',
                                                         'MAXREC': core.BATCH_SIZE})

    # patch sleep to stop the test from sleeping and slowing down execution
    @patch('cadcutils.net.ws.time.sleep', MagicMock(), create=True)
    @patch('cadcutils.net.auth.netrclib', MagicMock(), create=True)
    @patch('cadcutils.net.ws.Session.send')
    def test_post_observation(self, mock_conn):
        collection = 'cfht'
        observation_id = '7000000o'
        service = 'caom2repo'
        service_url = 'www.cadc.nrc.ca'

        obs = SimpleObservation(collection, observation_id)
        visitor = CAOM2RepoClient(auth.Subject(netrc='somenetrc'), host=service_url)
        response = MagicMock()
        response.status = 200
        mock_conn.return_value = response
        iobuffer = StringIO()
        ObservationWriter().write(obs, iobuffer)
        obsxml = iobuffer.getvalue()
        response.content = obsxml

        visitor.post_observation(obs)
        self.assertEqual('POST', mock_conn.call_args[0][0].method)
        self.assertEqual('/{}/auth/{}/{}'.format(service, collection, observation_id),
                         mock_conn.call_args[0][0].path_url)
        self.assertEqual('application/xml', mock_conn.call_args[0][0].headers['Content-Type'])
        self.assertEqual(obsxml, mock_conn.call_args[0][0].body)

        # signal problems
        http_error = requests.HTTPError()
        response.status_code = 500
        http_error.response = response
        response.raise_for_status.side_effect = [http_error]
        with self.assertRaises(exceptions.InternalServerException):
            visitor.post_observation(obs)

        # temporary transient errors
        http_error = requests.HTTPError()
        response.status_code = 503
        http_error.response = response
        response.raise_for_status.side_effect = [http_error, None]
        visitor.post_observation(obs)

        # permanent transient errors
        http_error = requests.HTTPError()
        response.status_code = 503
        http_error.response = response

        def raise_error(): raise http_error
        response.raise_for_status.side_effect = raise_error
        with self.assertRaises(exceptions.HttpException):
            visitor.post_observation(obs)

    # patch sleep to stop the test from sleeping and slowing down execution
    @patch('cadcutils.net.ws.time.sleep', MagicMock(), create=True)
    @patch('cadcutils.net.auth.os', MagicMock(), create=True)
    @patch('cadcutils.net.ws.Session.send')
    def test_put_observation(self, mock_conn):
        collection = 'cfht'
        observation_id = '7000000o'
        service = 'caom2repo'
        service_url = 'www.cadc.nrc.ca/'

        obs = SimpleObservation(collection, observation_id)
        subject = auth.Subject(certificate='somefile.pem')
        visitor = CAOM2RepoClient(subject, host=service_url)
        response = MagicMock()
        response.status = 200
        mock_conn.return_value = response
        iobuffer = StringIO()
        ObservationWriter().write(obs, iobuffer)
        obsxml = iobuffer.getvalue()
        response.content = obsxml

        visitor.put_observation(obs)
        self.assertEqual('PUT', mock_conn.call_args[0][0].method)
        self.assertEqual('/{}/pub/{}/{}'.format(service, collection, observation_id),
                         mock_conn.call_args[0][0].path_url)
        self.assertEqual('application/xml', mock_conn.call_args[0][0].headers['Content-Type'])
        self.assertEqual(obsxml, mock_conn.call_args[0][0].body)

        # signal problems
        http_error = requests.HTTPError()
        response.status_code = 500
        http_error.response = response
        response.raise_for_status.side_effect = [http_error]
        with self.assertRaises(exceptions.InternalServerException):
            visitor.put_observation(obs)

        # temporary transient errors
        http_error = requests.HTTPError()
        response.status_code = 503
        http_error.response = response
        response.raise_for_status.side_effect = [http_error, None]
        visitor.put_observation(obs)

        # permanent transient errors
        http_error = requests.HTTPError()
        response.status_code = 503
        http_error.response = response

        def raise_error(): raise http_error

        response.raise_for_status.side_effect = raise_error
        with self.assertRaises(exceptions.HttpException):
            visitor.put_observation(obs)

    # patch sleep to stop the test from sleeping and slowing down execution
    @patch('cadcutils.net.ws.time.sleep', MagicMock(), create=True)
    @patch('cadcutils.net.ws.Session.send')
    def test_delete_observation(self, mock_conn):
        collection = 'cfht'
        observation_id = '7000000o'
        service_url = 'www.cadc.nrc.ca/caom2repo'

        obs = SimpleObservation(collection, observation_id)
        visitor = CAOM2RepoClient(auth.Subject(), host=service_url)
        response = MagicMock()
        response.status = 200
        mock_conn.return_value = response

        visitor.delete_observation(collection, observation_id)
        self.assertEqual('DELETE', mock_conn.call_args[0][0].method)

        # signal problems
        http_error = requests.HTTPError()
        response.status_code = 500
        http_error.response = response
        response.raise_for_status.side_effect = [http_error]
        with self.assertRaises(exceptions.InternalServerException):
            visitor.delete_observation(collection, observation_id)

        # temporary transient errors
        http_error = requests.HTTPError()
        response.status_code = 503
        http_error.response = response
        response.raise_for_status.side_effect = [http_error, None]
        visitor.delete_observation(collection, observation_id)

        # permanent transient errors
        http_error = requests.HTTPError()
        response.status_code = 503
        http_error.response = response

        def raise_error(): raise http_error

        response.raise_for_status.side_effect = raise_error
        with self.assertRaises(exceptions.HttpException):
            visitor.delete_observation(collection, observation_id)

    def test_process(self):
        core.BATCH_SIZE = 3  # size of the batch is 3
        obs = [['a', 'b', 'c'], ['d'], []]
        visitor = CAOM2RepoClient(auth.Subject())
        visitor.get_observation = MagicMock(return_value=MagicMock(spec=SimpleObservation))
        visitor.post_observation = MagicMock()
        visitor._get_observations = MagicMock(side_effect=obs)

        self.assertEquals(4, visitor.visit(os.path.join(
                THIS_DIR, 'passplugin.py'), 'cfht'))

        obs = [['a', 'b', 'c'], ['d', 'e', 'f'], []]
        visitor._get_observations = MagicMock(side_effect=obs)
        self.assertEquals(6, visitor.visit(os.path.join(
                THIS_DIR, 'passplugin.py'), 'cfht'))

    @patch('caom2repo.core.CAOM2RepoClient')
    def test_main_app(self, client_mock):
        collection = 'cfht'
        observation_id = '7000000o'
        service = 'caom2repo'
        ifile = '/tmp/inputobs'

        obs = SimpleObservation(collection, observation_id)

        # test create
        with open(ifile, 'w') as infile:
            ObservationWriter().write(obs, infile)
        sys.argv = ["caom2tools", "create", '--resourceID', 'ivo://ca.nrc.ca/resource', ifile]
        core.main_app()
        client_mock.return_value.put_observation.assert_called_with(obs)

        # test update
        sys.argv = ["caom2tools", "update", '--resourceID', 'ivo://ca.nrc.ca/resource', ifile]
        core.main_app()
        client_mock.return_value.post_observation.assert_called_with(obs)

        # test read
        sys.argv = ["caom2tools", "read", '--resourceID', 'ivo://ca.nrc.ca/resource',
                    "--collection", collection, observation_id]
        client_mock.return_value.get_observation.return_value = obs
        core.main_app()
        client_mock.return_value.get_observation.assert_called_with(collection, observation_id)
        # repeat with output argument
        sys.argv = ["caom2tools", "read", '--resourceID', 'ivo://ca.nrc.ca/resource',
                    "--collection", collection, "--output", ifile, observation_id]
        client_mock.return_value.get_observation.return_value = obs
        core.main_app()
        client_mock.return_value.get_observation.assert_called_with(collection, observation_id)
        os.remove(ifile)

        # test delete
        sys.argv = ["caom2tools", "delete", '--resourceID', 'ivo://ca.nrc.ca/resource',
                    "--collection", collection, observation_id]
        core.main_app()
        client_mock.return_value.delete_observation.assert_called_with(collection=collection,
                                                                       observation_id=observation_id)

        # test visit
        # get the absolute path to be able to run the tests with the astropy frameworks
        plugin_file = THIS_DIR + "/passplugin.py"
        sys.argv = ["caom2tools", "visit", '--resourceID', 'ivo://ca.nrc.ca/resource',
                    "--plugin", plugin_file, "--start", "2012-01-01T11:22:33.44",
                    "--end", "2013-01-01T11:33:22.443", collection]
        with open(plugin_file, 'r') as infile:
            core.main_app()
            client_mock.return_value.visit.assert_called_with(
                ANY, collection,
                start=util.str2ivoa("2012-01-01T11:22:33.44"),
                end=util.str2ivoa("2013-01-01T11:33:22.443"))

    @patch('sys.exit', Mock(side_effect=[MyExitError, MyExitError, MyExitError,
                                         MyExitError, MyExitError, MyExitError]))
    def test_help(self):
        """ Tests the helper displays for commands and subcommands in main"""

        # expected helper messages
        usage =\
"""usage: caom2-client [-h] {create,read,update,delete,visit} ...

Client for a CAOM2 repo. In addition to CRUD (Create, Read, Update and Delete) operations it also implements a visitor operation that allows for updating multiple observations in a collection

positional arguments:
  {create,read,update,delete,visit}
    read                Read an existing observation
    update              Update an existing observation
    delete              Delete an existing observation
    visit               Visit observations in a collection

optional arguments:
  -h, --help            show this help message and exit
"""

        create_usage =\
"""usage: caom2-client create [-h] [-V]
                           [--cert CERT | -n | --netrc-file NETRC_FILE | -u USER]
                           [--host HOST] [--resourceID RESOURCEID]
                           [-d | -q | -v]
                           <new observation file in XML format>

Create a new observation

positional arguments:
  <new observation file in XML format>

optional arguments:
  --cert CERT           location of your X509 certificate to use for
                        authentication (unencrypted, in PEM format)
  -d, --debug           debug messages
  -h, --help            show this help message and exit
  --host HOST           Base hostname for services - used mainly for testing
                        (default: www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca)
  -n                    Use .netrc in $HOME for authentication
  --netrc-file NETRC_FILE
                        netrc file to use for authentication
  -q, --quiet           run quietly
  --resourceID RESOURCEID
                        resource identifier (default
                        ivo://cadc.nrc.ca/caom2repo)
  -u, --user USER       Name of user to authenticate. Note: application
                        prompts for the corresponding password!
  -v, --verbose         verbose messages
  -V, --version         show program's version number and exit
"""

        read_usage =\
"""usage: caom2-client read [-h] [-V]
                         [--cert CERT | -n | --netrc-file NETRC_FILE | -u USER]
                         [--host HOST] [--resourceID RESOURCEID]
                         [-d | -q | -v] --collection <collection>
                         [--output <destination file>]
                         <observation>

Read an existing observation

positional arguments:
  <observation>

optional arguments:
  --cert CERT           location of your X509 certificate to use for
                        authentication (unencrypted, in PEM format)
  --collection <collection>
  -d, --debug           debug messages
  -h, --help            show this help message and exit
  --host HOST           Base hostname for services - used mainly for testing
                        (default: www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca)
  -n                    Use .netrc in $HOME for authentication
  --netrc-file NETRC_FILE
                        netrc file to use for authentication
  --output, -o <destination file>
  -q, --quiet           run quietly
  --resourceID RESOURCEID
                        resource identifier (default
                        ivo://cadc.nrc.ca/caom2repo)
  -u, --user USER       Name of user to authenticate. Note: application
                        prompts for the corresponding password!
  -v, --verbose         verbose messages
  -V, --version         show program's version number and exit
"""

        update_usage =\
"""usage: caom2-client update [-h] [-V]
                           [--cert CERT | -n | --netrc-file NETRC_FILE | -u USER]
                           [--host HOST] [--resourceID RESOURCEID]
                           [-d | -q | -v]
                           <observation file>

Update an existing observation

positional arguments:
  <observation file>

optional arguments:
  --cert CERT           location of your X509 certificate to use for
                        authentication (unencrypted, in PEM format)
  -d, --debug           debug messages
  -h, --help            show this help message and exit
  --host HOST           Base hostname for services - used mainly for testing
                        (default: www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca)
  -n                    Use .netrc in $HOME for authentication
  --netrc-file NETRC_FILE
                        netrc file to use for authentication
  -q, --quiet           run quietly
  --resourceID RESOURCEID
                        resource identifier (default
                        ivo://cadc.nrc.ca/caom2repo)
  -u, --user USER       Name of user to authenticate. Note: application
                        prompts for the corresponding password!
  -v, --verbose         verbose messages
  -V, --version         show program's version number and exit
"""

        delete_usage =\
"""usage: caom2-client delete [-h] [-V]
                           [--cert CERT | -n | --netrc-file NETRC_FILE | -u USER]
                           [--host HOST] [--resourceID RESOURCEID]
                           [-d | -q | -v] --collection <collection>
                           <ID of observation>

Delete an existing observation

positional arguments:
  <ID of observation>

optional arguments:
  --cert CERT           location of your X509 certificate to use for
                        authentication (unencrypted, in PEM format)
  --collection <collection>
  -d, --debug           debug messages
  -h, --help            show this help message and exit
  --host HOST           Base hostname for services - used mainly for testing
                        (default: www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca)
  -n                    Use .netrc in $HOME for authentication
  --netrc-file NETRC_FILE
                        netrc file to use for authentication
  -q, --quiet           run quietly
  --resourceID RESOURCEID
                        resource identifier (default
                        ivo://cadc.nrc.ca/caom2repo)
  -u, --user USER       Name of user to authenticate. Note: application
                        prompts for the corresponding password!
  -v, --verbose         verbose messages
  -V, --version         show program's version number and exit
"""

        visit_usage =\
"""usage: caom2-client visit [-h] [-V]
                          [--cert CERT | -n | --netrc-file NETRC_FILE | -u USER]
                          [--host HOST] [--resourceID RESOURCEID]
                          [-d | -q | -v] --plugin <pluginClassFile>
                          [--start <datetime start point>]
                          [--end <datetime end point>]
                          [-s <CAOM2 service URL>]
                          <datacollection>

Visit observations in a collection

positional arguments:
  <datacollection>      data collection in CAOM2 repo

optional arguments:
  --cert CERT           location of your X509 certificate to use for
                        authentication (unencrypted, in PEM format)
  -d, --debug           debug messages
  --end <datetime end point>
                        earliest dataset to visit (UTC IVOA format: YYYY-mm-
                        ddTH:M:S.f)
  -h, --help            show this help message and exit
  --host HOST           Base hostname for services - used mainly for testing
                        (default: www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca)
  -n                    Use .netrc in $HOME for authentication
  --netrc-file NETRC_FILE
                        netrc file to use for authentication
  --plugin <pluginClassFile>
                        Plugin class to update each observation
  -q, --quiet           run quietly
  --resourceID RESOURCEID
                        resource identifier (default
                        ivo://cadc.nrc.ca/caom2repo)
  -s, --server <CAOM2 service URL>
                        URL of the CAOM2 repo server
  --start <datetime start point>
                        oldest dataset to visit (UTC IVOA format: YYYY-mm-
                        ddTH:M:S.f)
  -u, --user USER       Name of user to authenticate. Note: application
                        prompts for the corresponding password!
  -v, --verbose         verbose messages
  -V, --version         show program's version number and exit

Minimum plugin file format:
----
   from caom2 import Observation

   class ObservationUpdater:

    def update(self, observation):
        assert isinstance(observation, Observation), (
            'observation {} is not an Observation'.format(observation))
        # custom code to update the observation
----
"""

        self.maxDiff = None  # Display the entire difference
        # --help
        with patch('sys.stdout', new_callable=StringIO) as stdout_mock:
            sys.argv = ["caom2-client", "--help"]
            with self.assertRaises(MyExitError):
                core.main_app()
            self.assertEqual(usage, stdout_mock.getvalue())

        # create --help
        with patch('sys.stdout', new_callable=StringIO) as stdout_mock:
            sys.argv = ["caom2-client", "create", "--help"]
            with self.assertRaises(MyExitError):
                core.main_app()
            self.assertEqual(create_usage, stdout_mock.getvalue())
        #print(stdout_mock.getvalue())

        # read --help
        with patch('sys.stdout', new_callable=StringIO) as stdout_mock:
            sys.argv = ["caom2-client", "read", "--help"]
            with self.assertRaises(MyExitError):
                core.main_app()
            self.assertEqual(read_usage, stdout_mock.getvalue())
        #print(stdout_mock.getvalue())

        # update --help
        with patch('sys.stdout', new_callable=StringIO) as stdout_mock:
            sys.argv = ["caom2-client", "update", "--help"]
            with self.assertRaises(MyExitError):
                core.main_app()
            self.assertEqual(update_usage, stdout_mock.getvalue())
        #print(stdout_mock.getvalue())

        # delete --help
        with patch('sys.stdout', new_callable=StringIO) as stdout_mock:
            sys.argv = ["caom2-client", "delete", "--help"]
            with self.assertRaises(MyExitError):
                core.main_app()
            self.assertEqual(delete_usage, stdout_mock.getvalue())
        #print(stdout_mock.getvalue())

        # visit --help
        with patch('sys.stdout', new_callable=StringIO) as stdout_mock:
            sys.argv = ["caom2-client", "visit", "--help"]
            with self.assertRaises(MyExitError):
                core.main_app()
            self.assertEqual(visit_usage, stdout_mock.getvalue())
        #print(stdout_mock.getvalue())
