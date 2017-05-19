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
from six import StringIO
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

    @patch('cadcutils.net.ws.WsCapabilities')
    def test_plugin_class(self, caps_mock):
        caps_mock.get_service_host.return_value = 'some.host.com'
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
    @patch('cadcutils.net.ws.WsCapabilities')
    @patch('cadcutils.net.ws.time.sleep', MagicMock(), create=True)
    @patch('cadcutils.net.ws.open', MagicMock(), create=True)
    @patch('cadcutils.net.ws.Session.send')
    def test_get_observation(self, mock_get, caps_mock):
        caps_mock.get_service_host.return_value = 'some.host.com'
        caps_mock.return_value.get_access_url.return_value = 'http://serviceurl/caom2repo/pub'
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
    @patch('cadcutils.net.ws.WsCapabilities')
    @patch('cadcutils.net.ws.time.sleep', MagicMock(), create=True)
    @patch('cadcutils.net.ws.open', MagicMock(), create=True)
    @patch('caom2repo.core.net.BaseWsClient.get')
    def test_get_observations(self, mock_get, caps_mock):
        # This is almost similar to the previous test except that it gets
        # observations matching a collection and start/end criteria
        # Also, patch the CAOM2RepoClient now.
        caps_mock.get_service_host.return_value = 'some.host.com'
        caps_mock.return_value.get_access_url.return_value = 'http://serviceurl/caom2repo/pub'
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
        mock_get.assert_called_once_with((
            'vos://cadc.nrc.ca~vospace/CADC/std/CAOM2Repository#obs-1.0', 'cfht'),
            params={'MAXREC': core.BATCH_SIZE})

        mock_get.reset_mock()
        visitor._get_observations('cfht', end=datetime.strptime('2000-11-11', '%Y-%m-%d'))
        mock_get.assert_called_once_with((
            'vos://cadc.nrc.ca~vospace/CADC/std/CAOM2Repository#obs-1.0', 'cfht'),
            params={'END': '2000-11-11T00:00:00.000000', 'MAXREC': core.BATCH_SIZE})

        mock_get.reset_mock()
        visitor._get_observations('cfht',
                                  start=datetime.strptime('2000-11-11', '%Y-%m-%d'),
                                  end=datetime.strptime('2000-11-12', '%Y-%m-%d'))
        mock_get.assert_called_once_with((
            'vos://cadc.nrc.ca~vospace/CADC/std/CAOM2Repository#obs-1.0', 'cfht')
            , params={'START': '2000-11-11T00:00:00.000000',
            'END': '2000-11-12T00:00:00.000000', 'MAXREC': core.BATCH_SIZE})

    # patch sleep to stop the test from sleeping and slowing down execution
    @patch('cadcutils.net.ws.WsCapabilities')
    @patch('cadcutils.net.ws.time.sleep', MagicMock(), create=True)
    @patch('cadcutils.net.auth.netrclib', MagicMock(), create=True)
    @patch('cadcutils.net.ws.Session.send')
    def test_post_observation(self, mock_conn, caps_mock):
        caps_mock.get_service_host.return_value = 'some.host.com'
        caps_mock.return_value.get_access_url.return_value = 'http://serviceurl/caom2repo/auth'
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
    @patch('cadcutils.net.ws.WsCapabilities')
    @patch('cadcutils.net.ws.time.sleep', MagicMock(), create=True)
    @patch('cadcutils.net.auth.os', MagicMock(), create=True)
    @patch('cadcutils.net.ws.Session.send')
    def test_put_observation(self, mock_conn, caps_mock):
        caps_mock.get_service_host.return_value = 'some.host.com'
        caps_mock.return_value.get_access_url.return_value = 'http://serviceurl/caom2repo/pub'
        collection = 'cfht'
        observation_id = '7000000o'
        service = 'caom2repo'
        service_url = 'www.cadc.nrc.ca'

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
    @patch('cadcutils.net.ws.WsCapabilities')
    @patch('cadcutils.net.ws.time.sleep', MagicMock(), create=True)
    @patch('cadcutils.net.ws.Session.send')
    def test_delete_observation(self, mock_conn, caps_mock):
        caps_mock.get_service_host.return_value = 'some.host.com'
        caps_mock.return_value.get_access_url.return_value = 'http://serviceurl/caom2repo/pub'
        collection = 'cfht'
        observation_id = '7000000o'
        service_url = 'www.cadc.nrc.ca'

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


    @patch('cadcutils.net.ws.WsCapabilities')
    def test_process(self, caps_mock):
        caps_mock.get_service_host.return_value = 'some.host.com'
        core.BATCH_SIZE = 3  # size of the batch is 3
        obs = [['a', 'b', 'c'], ['d'], []]
        visitor = CAOM2RepoClient(auth.Subject())
        visitor.get_observation = MagicMock(return_value=MagicMock(spec=SimpleObservation))
        visitor.post_observation = MagicMock()
        visitor._get_observations = MagicMock(side_effect=obs)

        (visited, updated, skipped, failed) = visitor.visit(os.path.join(THIS_DIR, 'passplugin.py'), 'cfht')
        self.assertEqual(4, len(visited))
        self.assertEqual(4, len(updated))
        self.assertEqual(0, len(skipped))
        self.assertEqual(0, len(failed))

        obs = [['a', 'b', 'c'], ['d', 'e', 'f'], []]
        visitor._get_observations = MagicMock(side_effect=obs)
        (visited, updated, skipped, failed) = visitor.visit(os.path.join(THIS_DIR, 'passplugin.py'), 'cfht')
        self.assertEqual(6, len(visited))
        self.assertEqual(6, len(updated))
        self.assertEqual(0, len(skipped))
        self.assertEqual(0, len(failed))

        # make it return different status. errorplugin returns according to the
        # id of the observation: True for 'UPDATE', False for 'SKIP' and
        # raises exception for 'ERROR'
        obs_ids = [['UPDATE', 'SKIP', 'ERROR'], []]
        obs = [SimpleObservation(collection='TEST', observation_id='UPDATE'),
                SimpleObservation(collection='TEST', observation_id='SKIP'),
                SimpleObservation(collection='TEST', observation_id='ERROR')]
        visitor._get_observations = MagicMock(side_effect=obs_ids)
        visitor.get_observation = MagicMock(side_effect=obs)
        (visited, updated, skipped, failed) = visitor.visit(os.path.join(THIS_DIR, 'errorplugin.py'), 'cfht')
        self.assertEqual(3, len(visited))
        self.assertEqual(1, len(updated))
        self.assertEqual(1, len(skipped))
        self.assertEqual(1, len(failed))

        # repeat with other obs
        obs_ids = [['UPDATE', 'SKIP', 'ERROR'], ['UPDATE', 'SKIP']]
        obs = [SimpleObservation(collection='TEST', observation_id='UPDATE'),
               SimpleObservation(collection='TEST', observation_id='SKIP'),
               SimpleObservation(collection='TEST', observation_id='ERROR'),
               SimpleObservation(collection='TEST', observation_id='UPDATE'),
               SimpleObservation(collection='TEST', observation_id='SKIP')]
        visitor._get_observations = MagicMock(side_effect=obs_ids)
        visitor.get_observation = MagicMock(side_effect=obs)
        (visited, updated, skipped, failed) = visitor.visit(os.path.join(THIS_DIR, 'errorplugin.py'), 'cfht')
        self.assertEqual(5, len(visited))
        self.assertEqual(2, len(updated))
        self.assertEqual(2, len(skipped))
        self.assertEqual(1, len(failed))

        # repeat but halt on first ERROR -> process only 3 observations
        obs_ids = [['UPDATE', 'SKIP', 'ERROR'], ['UPDATE', 'SKIP']]
        obs = [SimpleObservation(collection='TEST', observation_id='UPDATE'),
               SimpleObservation(collection='TEST', observation_id='SKIP'),
               SimpleObservation(collection='TEST', observation_id='ERROR'),
               SimpleObservation(collection='TEST', observation_id='UPDATE'),
               SimpleObservation(collection='TEST', observation_id='SKIP')]
        visitor._get_observations = MagicMock(side_effect=obs_ids)
        visitor.get_observation = MagicMock(side_effect=obs)
        with self.assertRaises(SystemError):
            visitor.visit(os.path.join(
                THIS_DIR, 'errorplugin.py'), 'cfht', halt_on_error=True)


    def test_shortcuts(self):
        target = CAOM2RepoClient(auth.Subject())
        obs = SimpleObservation('CFHT', 'abc')

        target.put_observation = Mock()
        target.create(obs)
        target.put_observation.assert_called_with(obs)

        target.get_observation = Mock()
        target.read('CFHT', 'abc')
        target.get_observation.assert_called_with('CFHT', 'abc')

        target.post_observation = Mock()
        target.update(obs)
        target.post_observation.assert_called_with(obs)

        target.delete_observation = Mock()
        target.delete('CFHT', 'abc')
        target.delete_observation.assert_called_with('CFHT', 'abc')


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
        sys.argv = ["caom2tools", "create", '--resource-id', 'ivo://ca.nrc.ca/resource', ifile]
        core.main_app()
        client_mock.return_value.put_observation.assert_called_with(obs)

        # test update
        sys.argv = ["caom2tools", "update", '--resource-id', 'ivo://ca.nrc.ca/resource', ifile]
        core.main_app()
        client_mock.return_value.post_observation.assert_called_with(obs)

        # test read
        sys.argv = ["caom2tools", "read", '--resource-id', 'ivo://ca.nrc.ca/resource',
                    collection, observation_id]
        client_mock.return_value.get_observation.return_value = obs
        core.main_app()
        client_mock.return_value.get_observation.assert_called_with(collection, observation_id)
        # repeat with output argument
        sys.argv = ["caom2tools", "read", '--resource-id', 'ivo://ca.nrc.ca/resource',
                    "--output", ifile, collection, observation_id]
        client_mock.return_value.get_observation.return_value = obs
        core.main_app()
        client_mock.return_value.get_observation.assert_called_with(collection, observation_id)
        os.remove(ifile)

        # test delete
        sys.argv = ["caom2tools", "delete", '--resource-id', 'ivo://ca.nrc.ca/resource',
                    collection, observation_id]
        core.main_app()
        client_mock.return_value.delete_observation.assert_called_with(collection=collection,
                                                                       observation_id=observation_id)

        # test visit
        # get the absolute path to be able to run the tests with the astropy frameworks
        plugin_file = THIS_DIR + "/passplugin.py"
        sys.argv = ["caom2tools", "visit", '--resource-id', 'ivo://ca.nrc.ca/resource',
                    "--plugin", plugin_file, "--start", "2012-01-01T11:22:33",
                    "--end", "2013-01-01T11:33:22", collection]
        client_mock.return_value.visit.return_value = ['1'], ['1'], [], []
        with open(plugin_file, 'r') as infile:
            core.main_app()
            client_mock.return_value.visit.assert_called_with(
                ANY, collection, halt_on_error=False,
                start=core.str2date("2012-01-01T11:22:33"),
                end=core.str2date("2013-01-01T11:33:22"))

        # repeat visit test with halt-on-error
        sys.argv = ["caom2tools", "visit", '--resource-id', 'ivo://ca.nrc.ca/resource',
                    "--plugin", plugin_file,'--halt-on-error',
                    "--start", "2012-01-01T11:22:33",
                    "--end", "2013-01-01T11:33:22", collection]
        client_mock.return_value.visit.return_value = ['1'], ['1'], [], []
        with open(plugin_file, 'r') as infile:
            core.main_app()
            client_mock.return_value.visit.assert_called_with(
                ANY, collection, halt_on_error=True,
                start=core.str2date("2012-01-01T11:22:33"),
                end=core.str2date("2013-01-01T11:33:22"))

    @patch('sys.exit', Mock(side_effect=[MyExitError, MyExitError, MyExitError,
                                         MyExitError, MyExitError, MyExitError]))
    def test_help(self):
        """ Tests the helper displays for commands and subcommands in main"""

        # expected helper messages
        usage =\
"""usage: caom2-client [-h] [-V] {create,read,update,delete,visit} ...

Client for a CAOM2 repo. In addition to CRUD (Create, Read, Update and Delete) operations it also implements a visitor operation that allows for updating multiple observations in a collection

positional arguments:
  {create,read,update,delete,visit}
    create              Create a new observation
    read                Read an existing observation
    update              Update an existing observation
    delete              Delete an existing observation
    visit               Visit observations in a collection

optional arguments:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit
"""

        create_usage =\
"""usage: caom2-client create [-h]
                           [--cert CERT | -n | --netrc-file NETRC_FILE | -u USER]
                           [--host HOST] [--resource-id RESOURCE_ID]
                           [-d | -q | -v]
                           observation

Create a new observation

positional arguments:
  observation           XML file containing the observation

optional arguments:
  --cert CERT           location of your X509 certificate to use for
                        authentication (unencrypted, in PEM format)
  -d, --debug           debug messages
  -h, --help            show this help message and exit
  --host HOST           base hostname for services - used mainly for testing
                        (default: www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca)
  -n                    use .netrc in $HOME for authentication
  --netrc-file NETRC_FILE
                        netrc file to use for authentication
  -q, --quiet           run quietly
  --resource-id RESOURCE_ID
                        resource identifier (default
                        ivo://cadc.nrc.ca/caom2repo)
  -u, --user USER       name of user to authenticate. Note: application
                        prompts for the corresponding password!
  -v, --verbose         verbose messages
"""

        read_usage =\
"""usage: caom2-client read [-h]
                         [--cert CERT | -n | --netrc-file NETRC_FILE | -u USER]
                         [--host HOST] [--resource-id RESOURCE_ID]
                         [-d | -q | -v] [--output OUTPUT]
                         collection observationID

Read an existing observation

positional arguments:
  collection            collection name in CAOM2 repo
  observationID         observation identifier

optional arguments:
  --cert CERT           location of your X509 certificate to use for
                        authentication (unencrypted, in PEM format)
  -d, --debug           debug messages
  -h, --help            show this help message and exit
  --host HOST           base hostname for services - used mainly for testing
                        (default: www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca)
  -n                    use .netrc in $HOME for authentication
  --netrc-file NETRC_FILE
                        netrc file to use for authentication
  --output, -o OUTPUT   destination file
  -q, --quiet           run quietly
  --resource-id RESOURCE_ID
                        resource identifier (default
                        ivo://cadc.nrc.ca/caom2repo)
  -u, --user USER       name of user to authenticate. Note: application
                        prompts for the corresponding password!
  -v, --verbose         verbose messages
"""

        update_usage =\
"""usage: caom2-client update [-h]
                           [--cert CERT | -n | --netrc-file NETRC_FILE | -u USER]
                           [--host HOST] [--resource-id RESOURCE_ID]
                           [-d | -q | -v]
                           observation

Update an existing observation

positional arguments:
  observation           XML file containing the observation

optional arguments:
  --cert CERT           location of your X509 certificate to use for
                        authentication (unencrypted, in PEM format)
  -d, --debug           debug messages
  -h, --help            show this help message and exit
  --host HOST           base hostname for services - used mainly for testing
                        (default: www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca)
  -n                    use .netrc in $HOME for authentication
  --netrc-file NETRC_FILE
                        netrc file to use for authentication
  -q, --quiet           run quietly
  --resource-id RESOURCE_ID
                        resource identifier (default
                        ivo://cadc.nrc.ca/caom2repo)
  -u, --user USER       name of user to authenticate. Note: application
                        prompts for the corresponding password!
  -v, --verbose         verbose messages
"""

        delete_usage =\
"""usage: caom2-client delete [-h]
                           [--cert CERT | -n | --netrc-file NETRC_FILE | -u USER]
                           [--host HOST] [--resource-id RESOURCE_ID]
                           [-d | -q | -v]
                           collection observationID

Delete an existing observation

positional arguments:
  collection            collection name in CAOM2 repo
  observationID         observation identifier

optional arguments:
  --cert CERT           location of your X509 certificate to use for
                        authentication (unencrypted, in PEM format)
  -d, --debug           debug messages
  -h, --help            show this help message and exit
  --host HOST           base hostname for services - used mainly for testing
                        (default: www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca)
  -n                    use .netrc in $HOME for authentication
  --netrc-file NETRC_FILE
                        netrc file to use for authentication
  -q, --quiet           run quietly
  --resource-id RESOURCE_ID
                        resource identifier (default
                        ivo://cadc.nrc.ca/caom2repo)
  -u, --user USER       name of user to authenticate. Note: application
                        prompts for the corresponding password!
  -v, --verbose         verbose messages
"""

        visit_usage =\
"""usage: caom2-client visit [-h]
                          [--cert CERT | -n | --netrc-file NETRC_FILE | -u USER]
                          [--host HOST] [--resource-id RESOURCE_ID]
                          [-d | -q | -v] --plugin PLUGIN [--start START]
                          [--end END] [--halt-on-error] [-s SERVER]
                          collection

Visit observations in a collection

positional arguments:
  collection            data collection in CAOM2 repo

optional arguments:
  --cert CERT           location of your X509 certificate to use for
                        authentication (unencrypted, in PEM format)
  -d, --debug           debug messages
  --end END             latest observation to visit (UTC IVOA format: YYYY-mm-
                        ddTH:M:S)
  --halt-on-error       stop visitor on first update exception raised by
                        plugin
  -h, --help            show this help message and exit
  --host HOST           base hostname for services - used mainly for testing
                        (default: www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca)
  -n                    use .netrc in $HOME for authentication
  --netrc-file NETRC_FILE
                        netrc file to use for authentication
  --plugin PLUGIN       plugin class to update each observation
  -q, --quiet           run quietly
  --resource-id RESOURCE_ID
                        resource identifier (default
                        ivo://cadc.nrc.ca/caom2repo)
  -s, --server SERVER   URL of the CAOM2 repo server
  --start START         earliest observation to visit (UTC IVOA format: YYYY-
                        mm-ddTH:M:S)
  -u, --user USER       name of user to authenticate. Note: application
                        prompts for the corresponding password!
  -v, --verbose         verbose messages

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
