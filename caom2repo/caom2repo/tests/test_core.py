# # -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2022.                            (c) 2022.
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

import copy
import logging
import os
import sys
import unittest
from datetime import datetime

import requests
from cadcutils import util, exceptions
from cadcutils.net import auth
from caom2.obs_reader_writer import ObservationWriter
from caom2 import obs_reader_writer, ChecksumURI
from caom2.observation import SimpleObservation
from mock import Mock, patch, MagicMock, ANY, call
# TODO to be changed to io.BytesIO when caom2 is prepared for python3
from io import BytesIO, StringIO

from caom2repo import core
from caom2repo.core import CAOM2RepoClient

# The following is a temporary workaround for Python issue 25532
# (https://bugs.python.org/issue25532)
call.__wrapped__ = None


THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TESTDATA_DIR = os.path.join(THIS_DIR, 'data')


class MyExitError(Exception):
    pass


class PickableMagicMock(MagicMock):
    def __reduce__(self):
        return MagicMock, ()


class TestCAOM2Repo(unittest.TestCase):
    """Test the Caom2Visitor class"""

    def test_get_obs_from_file(self):
        level = logging.DEBUG
        visitor = CAOM2RepoClient(auth.Subject(), level)

        # no start or end
        with open(os.path.join(THIS_DIR, 'data/obs_id.txt')) as obs_file:
            obs_id_list = visitor._get_obs_from_file(obs_file, None, None,
                                                     False)
            self.assertEqual('obs_id_1', obs_id_list[0])
            self.assertEqual('obs_id_2', obs_id_list[1])
            self.assertEqual('obs_id_3', obs_id_list[2])

        # last_modified_date is earlier than start
        with open(os.path.join(THIS_DIR, 'data/obs_id.txt')) as obs_file:
            obs_id_list = visitor._get_obs_from_file(obs_file, util.str2ivoa(
                '2000-10-11T12:30:00.333'), None, False)
            self.assertEqual('obs_id_1', obs_id_list[0])

        # last_modified_date is between start and end
        with open(os.path.join(THIS_DIR, 'data/obs_id.txt')) as obs_file:
            obs_id_list = visitor._get_obs_from_file(
                obs_file, util.str2ivoa('2000-10-9T12:30:00.333'),
                util.str2ivoa('2016-10-11T12:30:00.333'), False)
            self.assertEqual('obs_id_1', obs_id_list[0])
            self.assertEqual('obs_id_2', obs_id_list[1])

        # last_modified_date is after end
        with open(os.path.join(THIS_DIR, 'data/obs_id.txt')) as obs_file:
            obs_id_list = visitor._get_obs_from_file(
                obs_file, util.str2ivoa('2000-10-9T12:30:00.333'),
                util.str2ivoa('2017-10-11T12:30:00.333'), False)
            self.assertEqual('obs_id_1', obs_id_list[0])
            self.assertEqual('obs_id_2', obs_id_list[1])
            self.assertEqual('obs_id_3', obs_id_list[2])

        # error in file
        with open(os.path.join(THIS_DIR, 'data/obs_id_error.txt')) as obs_file:
            with self.assertRaises(Exception):
                obs_id_list = visitor._get_obs_from_file(
                    obs_file, util.str2ivoa('2000-10-9T12:30:00.333'),
                    util.str2ivoa('2016-10-11T12:30:00.333'), True)

    @patch('caom2repo.core.net.BaseWsClient', Mock())
    def test_plugin_class(self):
        # plugin class does not change the observation
        collection = 'cfht'
        observation_id = '7000000o'
        level = logging.DEBUG
        visitor = CAOM2RepoClient(auth.Subject(), level)
        obs = SimpleObservation(collection, observation_id)
        expect_obs = copy.deepcopy(obs)
        visitor._load_plugin_class(os.path.join(THIS_DIR, 'passplugin.py'))
        visitor.plugin.update(obs)
        self.assertEqual(expect_obs, obs)

        # plugin class adds a plane to the observation
        visitor = CAOM2RepoClient(auth.Subject(), level)
        obs = SimpleObservation('cfht', '7000000o')
        expect_obs = copy.deepcopy(obs)
        visitor._load_plugin_class(os.path.join(THIS_DIR, 'addplaneplugin.py'))
        visitor.plugin.update(obs)
        self.assertNotEqual(expect_obs, obs)
        self.assertEqual(len(expect_obs.planes) + 1, len(obs.planes))

        # non-existent the plugin file
        with self.assertRaises(Exception):
            visitor._load_plugin_class(os.path.join(THIS_DIR, 'blah.py'))

        # non-existent ObservationUpdater class in the plugin file
        with self.assertRaises(Exception):
            visitor._load_plugin_class(
                os.path.join(THIS_DIR, 'test_visitor.py'))

        # non-existent update method in ObservationUpdater class
        with self.assertRaises(Exception):
            visitor._load_plugin_class(
                os.path.join(THIS_DIR, 'noupdateplugin.py'))

    # patch sleep to stop the test from sleeping and slowing down execution
    @patch('cadcutils.net.ws.WsCapabilities')
    @patch('cadcutils.net.ws.time.sleep', MagicMock(), create=True)
    @patch('cadcutils.net.ws.open', MagicMock(), create=True)
    @patch('cadcutils.net.ws.Session.send')
    def test_get_observation(self, mock_get, caps_mock):
        caps_mock.get_service_host.return_value = 'some.host.com'
        caps_mock.return_value.get_access_url.return_value =\
            'http://serviceurl/caom2repo/pub'
        collection = 'cfht'
        observation_id = '7000000o'
        service_url = 'www.cadc.nrc.ca/caom2repo'
        obs = SimpleObservation(collection, observation_id)
        writer = ObservationWriter()
        ibuffer = BytesIO()
        writer.write(obs, ibuffer)
        response = MagicMock()
        response.status_code = 200
        response.content = ibuffer.getvalue()
        mock_get.return_value = response
        ibuffer.seek(0)  # reposition the buffer for reading
        level = logging.DEBUG
        visitor = CAOM2RepoClient(auth.Subject(), level, host=service_url)
        self.assertEqual(obs, visitor.get_observation(
            collection, observation_id))

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
        visitor.read(collection, observation_id)

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
        caps_mock.return_value.get_access_url.return_value =\
            'http://serviceurl/caom2repo/pub'
        response = MagicMock()
        response.status_code = 200
        last_datetime = '2000-10-10T12:30:00.333'
        response.text = \
            ('CFHT\t700000o\t2000-10-10T12:20:11.123\t'
             '3e00ca6129dc8358315015204ab9fe15\nCFHT\t700001o\t' +
             last_datetime + '\t3e00ca6129dc8358315015204ab9fe15')
        mock_get.return_value = response

        level = logging.DEBUG
        visitor = CAOM2RepoClient(auth.Subject(), level)
        end_date = util.utils.str2ivoa(last_datetime)

        expect_observations = ['700000o', '700001o']
        self.assertEqual(expect_observations,
                         visitor._get_observations('cfht'))
        self.assertEqual(end_date, visitor._start)
        mock_get.assert_called_once_with((
            'vos://cadc.nrc.ca~vospace/CADC/std/CAOM2Repository#obs-1.2',
            'cfht'),
            params={'MAXREC': core.BATCH_SIZE})

        mock_get.reset_mock()
        visitor._get_observations('cfht', end=datetime.strptime('2000-11-11',
                                                                '%Y-%m-%d'))
        mock_get.assert_called_once_with((
            'vos://cadc.nrc.ca~vospace/CADC/std/CAOM2Repository#obs-1.2',
            'cfht'),
            params={'END': '2000-11-11T00:00:00.000',
                    'MAXREC': core.BATCH_SIZE})

        mock_get.reset_mock()
        visitor._get_observations('cfht',
                                  start=datetime.strptime('2000-11-11',
                                                          '%Y-%m-%d'),
                                  end=datetime.strptime('2000-11-12',
                                                        '%Y-%m-%d'))
        mock_get.assert_called_once_with((
            'vos://cadc.nrc.ca~vospace/CADC/std/CAOM2Repository#obs-1.2',
            'cfht'), params={'START': '2000-11-11T00:00:00.000',
                             'END': '2000-11-12T00:00:00.000',
                             'MAXREC': core.BATCH_SIZE})

    # patch sleep to stop the test from sleeping and slowing down execution
    @patch('cadcutils.net.ws.WsCapabilities')
    @patch('cadcutils.net.ws.time.sleep', MagicMock(), create=True)
    @patch('cadcutils.net.auth.netrclib', MagicMock(), create=True)
    @patch('cadcutils.net.ws.Session.send')
    def test_post_observation(self, mock_conn, caps_mock):
        caps_mock.get_service_host.return_value = 'some.host.com'
        caps_mock.return_value.get_access_url.return_value =\
            'http://serviceurl/caom2repo/auth'
        collection = 'cfht'
        observation_id = '7000000o'
        service = 'caom2repo'
        service_url = 'www.cadc.nrc.ca'

        obs = SimpleObservation(collection, observation_id)
        level = logging.DEBUG
        visitor = CAOM2RepoClient(auth.Subject(netrc='somenetrc'), level,
                                  host=service_url)
        response = MagicMock()
        response.status = 200
        mock_conn.return_value = response
        iobuffer = BytesIO()
        ObservationWriter().write(obs, iobuffer)
        obsxml = iobuffer.getvalue()
        response.content = obsxml

        visitor.post_observation(obs)
        self.assertEqual('POST', mock_conn.call_args[0][0].method)
        self.assertEqual(
            '/{}/auth/{}/{}'.format(service, collection, observation_id),
            mock_conn.call_args[0][0].path_url)
        self.assertEqual('application/xml',
                         mock_conn.call_args[0][0].headers['Content-Type'])
        self.assertEqual(obsxml, mock_conn.call_args[0][0].body)

        # signal problems
        http_error = requests.HTTPError()
        response.status_code = 500
        http_error.response = response
        response.raise_for_status.side_effect = [http_error]
        with self.assertRaises(exceptions.InternalServerException):
            visitor.update(obs)

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
        caps_mock.return_value.get_access_url.return_value =\
            'http://serviceurl/caom2repo/pub'
        collection = 'cfht'
        observation_id = '7000000o'
        service = 'caom2repo'
        service_url = 'www.cadc.nrc.ca'

        obs = SimpleObservation(collection, observation_id)
        subject = auth.Subject(certificate='somefile.pem')
        level = logging.DEBUG
        visitor = CAOM2RepoClient(subject, level, host=service_url)
        response = MagicMock()
        response.status = 200
        mock_conn.return_value = response
        iobuffer = BytesIO()
        ObservationWriter().write(obs, iobuffer)
        obsxml = iobuffer.getvalue()
        response.content = obsxml

        visitor.put_observation(obs)
        self.assertEqual('PUT', mock_conn.call_args[0][0].method)
        self.assertEqual(
            '/{}/pub/{}/{}'.format(service, collection, observation_id),
            mock_conn.call_args[0][0].path_url)
        self.assertEqual('application/xml',
                         mock_conn.call_args[0][0].headers['Content-Type'])
        self.assertEqual(obsxml, mock_conn.call_args[0][0].body)

        # signal problems
        http_error = requests.HTTPError()
        response.status_code = 500
        http_error.response = response
        response.raise_for_status.side_effect = [http_error]
        with self.assertRaises(exceptions.InternalServerException):
            visitor.create(obs)

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
        caps_mock.return_value.get_access_url.return_value =\
            'http://serviceurl/caom2repo/pub'
        collection = 'cfht'
        observation_id = '7000000o'
        service_url = 'www.cadc.nrc.ca'
        level = logging.DEBUG

        visitor = CAOM2RepoClient(auth.Subject(), level, host=service_url)
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
            visitor.delete(collection, observation_id)

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

    @patch('caom2repo.core.net.BaseWsClient', Mock())
    def test_process(self):
        core.BATCH_SIZE = 3  # size of the batch is 3
        obs = [['a', 'b', 'c'], ['d'], []]
        level = logging.DEBUG
        visitor = CAOM2RepoClient(auth.Subject(), level)
        visitor.get_observation = MagicMock(
            return_value=MagicMock(spec=SimpleObservation))
        visitor.post_observation = MagicMock()
        visitor._get_observations = MagicMock(side_effect=obs)

        (visited, updated, skipped, failed) = visitor.visit(
            os.path.join(THIS_DIR, 'passplugin.py'), 'cfht')
        self.assertEqual(4, len(visited))
        self.assertEqual(4, len(updated))
        self.assertEqual(0, len(skipped))
        self.assertEqual(0, len(failed))

        obs = [['a', 'b', 'c'], ['d', 'e', 'f'], []]
        visitor._get_observations = MagicMock(side_effect=obs)
        (visited, updated, skipped, failed) = visitor.visit(
            os.path.join(THIS_DIR, 'passplugin.py'), 'cfht')
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
        (visited, updated, skipped, failed) = visitor.visit(
            os.path.join(THIS_DIR, 'errorplugin.py'), 'cfht')
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
        (visited, updated, skipped, failed) = visitor.visit(
            os.path.join(THIS_DIR, 'errorplugin.py'), 'cfht')
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

        # test with time boundaries
        core.BATCH_SIZE = 3  # size of the batch is 3
        response = MagicMock()
        response.text = """ARCHIVE\ta\t2011-01-01T11:00:00.000
                           ARCHIVE\tb\t211-01-01T11:00:10.000
                           ARCHIVE\tc\t2011-01-01T12:00:00.000"""
        response2 = MagicMock()
        response2.text = """ARCHIVE\td\t2011-02-02T11:00:00.000"""
        level = logging.DEBUG
        visitor = CAOM2RepoClient(auth.Subject(), level)
        visitor.get_observation = MagicMock(
            return_value=MagicMock(spec=SimpleObservation))
        visitor.post_observation = MagicMock()
        visitor._repo_client.get = MagicMock(side_effect=[response, response2])

        start = '2010-10-10T12:00:00.000'
        end = '2012-12-12T11:11:11.000'
        (visited, updated, skipped, failed) = visitor.visit(
            os.path.join(THIS_DIR, 'passplugin.py'), 'cfht',
            start=util.str2ivoa(start),
            end=util.str2ivoa(end))

        self.assertEqual(4, len(visited))
        self.assertEqual(4, len(updated))
        self.assertEqual(0, len(skipped))
        self.assertEqual(0, len(failed))
        calls = [call((core.CURRENT_CAOM2REPO_OBS_CAPABILITY_ID, 'cfht'),
                      params={'START': start, 'END': end, 'MAXREC': 3}),
                 call((core.CURRENT_CAOM2REPO_OBS_CAPABILITY_ID, 'cfht'),
                      params={'START': '2011-01-01T12:00:00.000',
                              # datetime of the last record in the batch
                              'END': end,
                              'MAXREC': 3})]
        visitor._repo_client.get.assert_has_calls(calls)

    @patch('caom2repo.core.net.BaseWsClient', Mock())
    def test_visit_retry_on_412(self):
        # observation changed on server while visited
        core.BATCH_SIZE = 3  # size of the batch is 3
        obs = [['a'], []]
        level = logging.DEBUG
        visitor = CAOM2RepoClient(auth.Subject(), level)
        observation = SimpleObservation('cfht', 'a')
        observation.acc_meta_checksum = ChecksumURI('md5:abc')
        visitor.get_observation = MagicMock(side_effect=[observation,
                                                         observation])

        exception_412 = exceptions.UnexpectedException()
        exception_412.orig_exception = Mock()
        exception_412.orig_exception.response = Mock(status_code=412)
        visitor.post_observation = MagicMock(side_effect=[exception_412, None])
        visitor._get_observations = MagicMock(side_effect=obs)

        (visited, updated, skipped, failed) = visitor.visit(
            os.path.join(THIS_DIR, 'passplugin.py'), 'cfht')
        self.assertEqual(1, len(visited))
        self.assertEqual(1, len(updated))
        self.assertEqual(0, len(skipped))
        self.assertEqual(0, len(failed))
        # get and post called twice to recover from error HTTP status 412 -
        # precondition
        self.assertEqual(2, visitor.get_observation.call_count)
        self.assertEqual(2, visitor.post_observation.call_count)
        visitor.post_observation.assert_called_with(
            observation, observation.acc_meta_checksum.uri)

    def mock_get_observation(self, collection, observationID):
        return SimpleObservation(collection, observationID)

    def mock_get_observation_with_expected_type_error(self, collection,
                                                      observationID):
        raise TypeError("unexpected keyword argument")

    def mock_get_observation_with_unexpected_type_error(self, collection,
                                                        observationID):
        raise TypeError("unexpected TypeError")

    def mock_post_observation_with_exception(self, observation,
                                             obs_checksum=None):
        raise Exception("exception with observation")

    @patch('caom2repo.core.CAOM2RepoClient.get_observation')
    @patch('caom2repo.core.CAOM2RepoClient.post_observation')
    def test_multiprocess_with_exception(self, get_mock, post_mock):
        core.BATCH_SIZE = 3  # size of the batch is 3
        obs_ids = [['a', 'b', 'c'], ['d'], []]
        get_mock.side_effect = \
            self.mock_get_observation
        post_mock.side_effect = \
            self.mock_post_observation_with_exception
        level = logging.DEBUG
        visitor = CAOM2RepoClient(auth.Subject(), level)
        visitor.get_observation = PickableMagicMock(
            return_value=PickableMagicMock(spec=SimpleObservation))
        visitor.post_observation = PickableMagicMock()
        visitor._get_observations = PickableMagicMock(side_effect=obs_ids)

        try:
            (visited, updated, skipped, failed) = visitor.visit(
                os.path.join(THIS_DIR, 'passplugin.py'), 'cfht', start=None,
                end=None, obs_file=None, nthreads=3,
                halt_on_error=True)
        except Exception as e:
            self.assertTrue("exception with observation" in str(e))
        finally:
            logging.info("DONE")

    @patch('caom2repo.core.CAOM2RepoClient.get_observation')
    @patch('caom2repo.core.CAOM2RepoClient.post_observation', Mock())
    def test_multiprocess_with_expected_type_error(self, get_mock):
        core.BATCH_SIZE = 3  # size of the batch is 3
        obs_ids = [['a', 'b', 'c'], ['d'], []]
        get_mock.side_effect = \
            self.mock_get_observation_with_expected_type_error
        level = logging.DEBUG
        visitor = CAOM2RepoClient(auth.Subject(), level)
        visitor.get_observation = PickableMagicMock(
            return_value=PickableMagicMock(spec=SimpleObservation))
        visitor.post_observation = PickableMagicMock()
        visitor._get_observations = PickableMagicMock(side_effect=obs_ids)

        try:
            (visited, updated, skipped, failed) = visitor.visit(
                os.path.join(THIS_DIR, 'passplugin.py'), 'cfht', start=None,
                end=None, obs_file=None, nthreads=3,
                halt_on_error=True)
        except RuntimeError as e:
            self.assertTrue("To fix the problem" in str(e))
        finally:
            logging.info("DONE")

    @patch('caom2repo.core.CAOM2RepoClient.get_observation')
    @patch('caom2repo.core.CAOM2RepoClient.post_observation', Mock())
    def test_multiprocess_with_unexpected_type_error(self, get_mock):
        core.BATCH_SIZE = 3  # size of the batch is 3
        obs_ids = [['a', 'b', 'c'], ['d'], []]
        get_mock.side_effect = \
            self.mock_get_observation_with_unexpected_type_error
        level = logging.DEBUG
        visitor = CAOM2RepoClient(auth.Subject(), level)
        visitor.get_observation = PickableMagicMock(
            return_value=PickableMagicMock(spec=SimpleObservation))
        visitor.post_observation = PickableMagicMock()
        visitor._get_observations = PickableMagicMock(side_effect=obs_ids)

        try:
            (visited, updated, skipped, failed) = visitor.visit(
                os.path.join(THIS_DIR, 'passplugin.py'), 'cfht', start=None,
                end=None, obs_file=None, nthreads=3,
                halt_on_error=True)
        except TypeError as e:
            self.assertTrue("unexpected TypeError" in str(e))
        finally:
            logging.info("DONE")

    @patch('caom2repo.core.CAOM2RepoClient.get_observation')
    @patch('caom2repo.core.CAOM2RepoClient.post_observation', Mock())
    def test_multiprocess_with_obs_id(self, get_mock):
        core.BATCH_SIZE = 3  # size of the batch is 3
        obs_ids = [['a', 'b', 'c'], ['d'], []]
        get_mock.side_effect = self.mock_get_observation
        level = logging.DEBUG
        visitor = CAOM2RepoClient(auth.Subject(), level)
        visitor.get_observation = PickableMagicMock(
            return_value=PickableMagicMock(spec=SimpleObservation))
        visitor.post_observation = PickableMagicMock()
        visitor._get_observations = PickableMagicMock(side_effect=obs_ids)

        (visited, updated, skipped, failed) = visitor.visit(
            os.path.join(THIS_DIR, 'passplugin.py'), 'cfht', start=None,
            end=None, obs_file=None, nthreads=3)

        try:
            self.assertEqual(4, len(visited))
            self.assertEqual(4, len(updated))
            self.assertEqual(0, len(skipped))
            self.assertEqual(0, len(failed))
            self.assertTrue('a' in visited)
            self.assertTrue('b' in visited)
            self.assertTrue('c' in visited)
            self.assertTrue('d' in visited)
            self.assertFalse('e' in visited)
            self.assertTrue('a' in updated)
            self.assertTrue('b' in updated)
            self.assertTrue('c' in updated)
            self.assertTrue('d' in updated)
            self.assertFalse('e' in updated)
        finally:
            # lp.join()
            logging.info("DONE")

    @patch('caom2repo.core.CAOM2RepoClient.get_observation')
    @patch('caom2repo.core.CAOM2RepoClient.post_observation', Mock())
    def test_multiprocess_with_more_obs_id(self, get_mock):
        core.BATCH_SIZE = 3  # size of the batch is 3
        obs_ids = [['a', 'b', 'c'], ['d', 'e', 'f'], []]
        get_mock.side_effect = self.mock_get_observation
        level = logging.DEBUG
        visitor = CAOM2RepoClient(auth.Subject(), level)
        visitor.get_observation = PickableMagicMock(
            return_value=PickableMagicMock(spec=SimpleObservation))
        visitor.post_observation = PickableMagicMock()
        visitor._get_observations = PickableMagicMock(side_effect=obs_ids)

        (visited, updated, skipped, failed) = visitor.visit(
            os.path.join(THIS_DIR, 'passplugin.py'), 'cfht', start=None,
            end=None, obs_file=None, nthreads=3)

        try:
            self.assertEqual(6, len(visited))
            self.assertEqual(6, len(updated))
            self.assertEqual(0, len(skipped))
            self.assertEqual(0, len(failed))
            self.assertTrue('a' in visited)
            self.assertTrue('b' in visited)
            self.assertTrue('c' in visited)
            self.assertTrue('d' in visited)
            self.assertTrue('e' in visited)
            self.assertTrue('f' in visited)
            self.assertFalse('g' in visited)
            self.assertTrue('a' in updated)
            self.assertTrue('b' in updated)
            self.assertTrue('c' in updated)
            self.assertTrue('d' in updated)
            self.assertTrue('e' in updated)
            self.assertTrue('f' in updated)
            self.assertFalse('g' in updated)
        finally:
            # lp.join()
            logging.info("DONE")

    @patch('caom2repo.core.CAOM2RepoClient.get_observation')
    @patch('caom2repo.core.CAOM2RepoClient.post_observation', Mock())
    def test_multiprocess_with_different_statuses(self, get_mock):
        core.BATCH_SIZE = 3  # size of the batch is 3
        # make it return different status. errorplugin returns according to the
        # id of the observation: True for 'UPDATE', False for 'SKIP' and
        # raises exception for 'ERROR'
        obs_ids = [['UPDATE', 'SKIP', 'ERROR'], []]
        get_mock.side_effect = self.mock_get_observation
        level = logging.DEBUG
        visitor = CAOM2RepoClient(auth.Subject(), level)
        visitor.get_observation = PickableMagicMock(
            return_value=PickableMagicMock(spec=SimpleObservation))
        visitor.post_observation = PickableMagicMock()
        visitor._get_observations = PickableMagicMock(side_effect=obs_ids)

        (visited, updated, skipped, failed) = visitor.visit(
            os.path.join(THIS_DIR, 'errorplugin.py'), 'cfht', start=None,
            end=None, obs_file=None, nthreads=3)

        try:
            self.assertEqual(3, len(visited))
            self.assertEqual(1, len(updated))
            self.assertEqual(1, len(skipped))
            self.assertEqual(1, len(failed))
        finally:
            # lp.join()
            logging.info("DONE")

    @patch('caom2repo.core.CAOM2RepoClient.get_observation')
    @patch('caom2repo.core.CAOM2RepoClient.post_observation', Mock())
    def test_multiprocess_with_more_different_statuses(self, get_mock):
        core.BATCH_SIZE = 3  # size of the batch is 3
        # make it return different status. errorplugin returns according to the
        # id of the observation: True for 'UPDATE', False for 'SKIP' and
        # raises exception for 'ERROR'
        obs_ids = [['UPDATE', 'SKIP', 'ERROR'], ['UPDATE', 'SKIP']]
        get_mock.side_effect = self.mock_get_observation
        level = logging.DEBUG
        visitor = CAOM2RepoClient(auth.Subject(), level)
        visitor.get_observation = PickableMagicMock(
            return_value=PickableMagicMock(spec=SimpleObservation))
        visitor.post_observation = PickableMagicMock()
        visitor._get_observations = PickableMagicMock(side_effect=obs_ids)

        (visited, updated, skipped, failed) = visitor.visit(
            os.path.join(THIS_DIR, 'errorplugin.py'), 'cfht', start=None,
            end=None, obs_file=None, nthreads=3)

        try:
            self.assertEqual(5, len(visited))
            self.assertEqual(2, len(updated))
            self.assertEqual(2, len(skipped))
            self.assertEqual(1, len(failed))
        finally:
            # lp.join()
            logging.info("DONE")

    def test_shortcuts(self):
        level = logging.DEBUG
        target = CAOM2RepoClient(auth.Subject(), level)
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
        ifile = '/tmp/inputobs'

        obs = SimpleObservation(collection, observation_id)

        # test create
        with open(ifile, 'wb') as infile:
            ObservationWriter().write(obs, infile)
        sys.argv = ["caom2tools", "create", '--resource-id',
                    'ivo://ca.nrc.ca/resource', ifile]
        core.main_app()
        client_mock.return_value.put_observation.assert_called_with(obs)

        # test update
        sys.argv = ["caom2tools", "update", '--resource-id',
                    'ivo://ca.nrc.ca/resource', ifile]
        core.main_app()
        client_mock.return_value.post_observation.assert_called_with(obs)

        # test read
        sys.argv = ["caom2tools", "read", '--resource-id',
                    'ivo://ca.nrc.ca/resource',
                    collection, observation_id]
        client_mock.return_value.get_observation.return_value = obs
        client_mock.return_value.namespace = obs_reader_writer.CAOM24_NAMESPACE
        core.main_app()
        client_mock.return_value.get_observation.\
            assert_called_with(collection, observation_id)
        # repeat with output argument
        sys.argv = ["caom2tools", "read", '--resource-id',
                    'ivo://ca.nrc.ca/resource',
                    "--output", ifile, collection, observation_id]
        client_mock.return_value.get_observation.return_value = obs
        core.main_app()
        client_mock.return_value.get_observation.\
            assert_called_with(collection, observation_id)
        os.remove(ifile)

        # test delete
        sys.argv = ["caom2tools", "delete", '--resource-id',
                    'ivo://ca.nrc.ca/resource',
                    collection, observation_id]
        core.main_app()
        client_mock.return_value.delete_observation.assert_called_with(
            collection=collection,
            observation_id=observation_id)

        # test visit
        # get the absolute path to be able to run the tests with the
        # astropy frameworks
        plugin_file = THIS_DIR + "/passplugin.py"
        sys.argv = ["caom2tools", "visit", '--resource-id',
                    'ivo://ca.nrc.ca/resource',
                    "--plugin", plugin_file, "--start", "2012-01-01T11:22:33",
                    "--end", "2013-01-01T11:33:22", collection]
        client_mock.return_value.visit.return_value = ['1'], ['1'], [], []
        with open(plugin_file, 'r') as infile:
            core.main_app()
            client_mock.return_value.visit.assert_called_with(
                ANY, collection, halt_on_error=False, nthreads=None,
                obs_file=None,
                start=core.str2date("2012-01-01T11:22:33"),
                end=core.str2date("2013-01-01T11:33:22"))

        # repeat visit test with halt-on-error
        sys.argv = ["caom2tools", "visit", '--resource-id',
                    'ivo://ca.nrc.ca/resource',
                    "--plugin", plugin_file, '--halt-on-error',
                    "--start", "2012-01-01T11:22:33",
                    "--end", "2013-01-01T11:33:22", collection]
        client_mock.return_value.visit.return_value = ['1'], ['1'], [], []
        with open(plugin_file, 'r') as infile:
            core.main_app()
            client_mock.return_value.visit.assert_called_with(
                ANY, collection, halt_on_error=True, nthreads=None,
                obs_file=None,
                start=core.str2date("2012-01-01T11:22:33"),
                end=core.str2date("2013-01-01T11:33:22"))

    @patch('sys.exit', Mock(side_effect=[MyExitError, MyExitError, MyExitError,
                                         MyExitError, MyExitError, MyExitError,
                                         MyExitError, MyExitError]))
    def test_help(self):
        """ Tests the helper displays for commands and subcommands in main"""

        # expected helper messages
        with open(os.path.join(TESTDATA_DIR, 'help.txt'), 'r') as myfile:
            usage = myfile.read()
        with open(
                os.path.join(TESTDATA_DIR, 'create_help.txt'), 'r') as myfile:
            create_usage = myfile.read()
        with open(os.path.join(TESTDATA_DIR, 'read_help.txt'), 'r') as myfile:
            read_usage = myfile.read()
        with open(
                os.path.join(TESTDATA_DIR, 'update_help.txt'), 'r') as myfile:
            update_usage = myfile.read()
        with open(
                os.path.join(TESTDATA_DIR, 'delete_help.txt'), 'r') as myfile:
            delete_usage = myfile.read()
        with open(os.path.join(TESTDATA_DIR, 'visit_help.txt'), 'r') as myfile:
            visit_usage = myfile.read()

        self.maxDiff = None  # Display the entire difference
        # --help
        with patch('sys.stdout', new_callable=StringIO) as stdout_mock:
            sys.argv = ["caom2-repo", "--help"]
            with self.assertRaises(MyExitError):
                core.main_app()
            # Python 3.10 difference in titles
            actual = stdout_mock.getvalue().replace(
                'options:', 'optional arguments:').strip('\n')
            assert usage.strip('\n') == actual

        # create --help
        with patch('sys.stdout', new_callable=StringIO) as stdout_mock:
            sys.argv = ["caom2-repo", "create", "--help"]
            with self.assertRaises(MyExitError):
                core.main_app()
            actual = stdout_mock.getvalue().replace(
                'options:', 'optional arguments:').strip('\n')
            assert create_usage.strip('\n') == actual

        # read --help
        with patch('sys.stdout', new_callable=StringIO) as stdout_mock:
            sys.argv = ["caom2-repo", "read", "--help"]
            with self.assertRaises(MyExitError):
                core.main_app()
            actual = stdout_mock.getvalue().replace(
                'options:', 'optional arguments:').strip('\n')
            assert read_usage.strip('\n') == actual

        # update --help
        with patch('sys.stdout', new_callable=StringIO) as stdout_mock:
            sys.argv = ["caom2-repo", "update", "--help"]
            with self.assertRaises(MyExitError):
                core.main_app()
            actual = stdout_mock.getvalue().replace(
                'options:', 'optional arguments:').strip('\n')
            assert update_usage.strip('\n') == actual

        # delete --help
        with patch('sys.stdout', new_callable=StringIO) as stdout_mock:
            sys.argv = ["caom2-repo", "delete", "--help"]
            with self.assertRaises(MyExitError):
                core.main_app()
            actual = stdout_mock.getvalue().replace(
                'options:', 'optional arguments:').strip('\n')
            assert delete_usage.strip('\n') == actual

        # visit --help
        with patch('sys.stdout', new_callable=StringIO) as stdout_mock:
            sys.argv = ["caom2-repo", "visit", "--help"]
            with self.assertRaises(MyExitError):
                core.main_app()
            actual = stdout_mock.getvalue().replace(
                'options:', 'optional arguments:').strip('\n')
            assert visit_usage.strip('\n') == actual

        # visit too few number of threads
        with patch('sys.stderr', new_callable=StringIO) as stderr_mock:
            sys.argv = ["caom2-repo", "visit", "--threads", "1",
                        "--obs_file",
                        os.path.join(THIS_DIR, 'data/obs_id.txt'),
                        "--plugin", os.path.join(THIS_DIR, 'passplugin.py'),
                        "TEST"]
            with self.assertRaises(MyExitError):
                core.main_app()
            self.assertTrue('error: argument --threads: invalid choice' in
                            stderr_mock.getvalue())

        # visit too many number of threads
        with patch('sys.stderr', new_callable=StringIO) as stderr_mock:
            sys.argv = ["caom2-repo", "visit", "--threads", "11",
                        "--obs_file",
                        os.path.join(THIS_DIR, 'data/obs_id.txt'),
                        "--plugin",
                        os.path.join(THIS_DIR, 'passplugin.py'), "TEST"]
            with self.assertRaises(MyExitError):
                core.main_app()
            self.assertTrue('error: argument --threads: invalid choice' in
                            stderr_mock.getvalue())
