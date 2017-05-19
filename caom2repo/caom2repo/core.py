#!/usr/bin/env python2.7
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

import argparse
import imp
import logging
import os
import os.path
import sys
from six import StringIO
from datetime import datetime

from cadcutils import net
from cadcutils import util

from caom2repo import version

from caom2.obs_reader_writer import ObservationReader, ObservationWriter
from caom2.version import version as caom2_version
from six.moves.urllib.parse import urlparse

# from . import version as caom2repo_version
from caom2repo import version

__all__ = ['CAOM2RepoClient']

BATCH_SIZE = int(10000)
DATE_FORMAT = "%Y-%m-%dT%H:%M:%S.%f" #IVOA dateformat

CAOM2REPO_OBS_CAPABILITY_ID = 'vos://cadc.nrc.ca~vospace/CADC/std/CAOM2Repository#obs-1.0'
# IVOA dateformat
DATE_FORMAT = "%Y-%m-%dT%H:%M:%S.%f"
# resource ID for info
DEFAULT_RESOURCE_ID = 'ivo://cadc.nrc.ca/caom2repo'
APP_NAME = 'caom2repo'


class CAOM2RepoClient(object):

    """Class to do CRUD + visitor actions on a CAOM2 collection repo."""

    logger = logging.getLogger('CAOM2RepoClient')

    def __init__(self, subject, resource_id=DEFAULT_RESOURCE_ID, host=None):
        """
        Instance of a CAOM2RepoClient
        :param subject: the subject performing the action
        :type cadcutils.auth.Subject
        :param server: Host server for the caom2repo service
        """
        agent = '{}/{}'.format(APP_NAME, version.version)
        self.host = host
        self._repo_client = net.BaseWsClient(resource_id, subject, agent, retry=True, host=host)

        agent = "caom2-repo-client/{} caom2/{}".format(version.version, caom2_version)

        self._repo_client = net.BaseWsClient(resource_id, subject,
                                             agent, retry=True, host=self.host)

    # shortcuts for the CRUD operations
    def create(self, observation):
        """
        Creates an observation in the repo. 
        :param observation: Observation to create
        :return: Created observation
        """
        self.put_observation(observation)

    def read(self, collection, observation_id):
        """
        Read an observation from the repo
        :param collection: Name of the collection
        :param observation_id: Observation identifier
        :return: Observation
        """
        self.get_observation(collection, observation_id)

    def update(self, observation):
        """
        Update an observation in the repo
        :param observation: Observation to update
        :return: Updated observation
        """
        self.post_observation(observation)

    def delete(self, collection, observation_id):
        """
        Delete an observation from the repo
        :param collection: Name of the collection
        :param observation_id: Observation identifier
        """
        self.delete_observation(collection, observation_id)


    def visit(self, plugin, collection, start=None, end=None, halt_on_error=False):
        """
        Main processing function that iterates through the observations of
        the collection and updates them according to the algorithm
        of the plugin function
        :param plugin: path to python file that contains the algorithm to be applied to visited
                        observations
        :param collection: name of the CAOM2 collection
        :param start: optional earliest date-time of the targeted observation set
        :param end: optional latest date-time of the targeted observation set
        :param halt_on_error if True halts the execution on the first exception raised by
               the plugin update function otherwise logs the error and continues
        :return: tuple (list of visited observations, list of updated observation,
                list of skipped observations, list of failured observation)
        """
        if not os.path.isfile(plugin):
            raise Exception('Cannot find plugin file ' + plugin)
        assert collection is not None
        if start is not None:
            assert type(start) is datetime
        if end is not None:
            assert type(end) is datetime
        self._load_plugin_class(plugin)

        # this is updated by _get_observations with the timestamp of last observation in the batch
        self._start = start
        visited = []
        failed = []
        updated = []
        skipped = []
        observations = self._get_observations(collection, self._start, end)
        while len(observations) > 0:
            for observationID in observations:
                self.logger.info('Process observation: ' + observationID)
                observation = self.get_observation(collection, observationID)
                try:
                    if self.plugin.update(observation) is False:
                        self.logger.info('SKIP {}'.format(observation.observation_id))
                        skipped.append(observation.observation_id)
                    else:
                        self.post_observation(observation)
                        self.logger.debug('UPDATED {}'.format(observation.observation_id))
                        updated.append(observation.observation_id)
                except Exception as e:
                    failed.append(observation.observation_id)
                    self.logger.error('FAILED {} - Reason: '.
                                  format(observation.observation_id), e)
                    if halt_on_error:
                        raise e
                visited.append(observation.observation_id)
            if len(observations) == BATCH_SIZE:
                observations = self._get_observations(collection)
            else:
                # the last batch was smaller so it must have been the last
                break
        return visited, updated, skipped, failed

    def _get_observations(self, collection, start=None, end=None):
        """
        Returns a list of observations from the collection
        :param collection: name of the collection
        :param start: earliest observation
        :param end: latest observation
        :return: list of observation ids
        """
        assert collection is not None
        observations = []
        params = {'MAXREC': BATCH_SIZE}
        if start is not None:
            params['START'] = start.strftime(DATE_FORMAT)
        if end is not None:
            params['END'] = end.strftime(DATE_FORMAT)

        response = self._repo_client.get((CAOM2REPO_OBS_CAPABILITY_ID, collection),
                                         params=params)
        last_datetime = None
        for line in response.content.splitlines():
            (obs, last_datetime) = line.split(',')
            observations.append(obs)
        if last_datetime is not None:
            self._start = datetime.strptime(last_datetime, DATE_FORMAT)
        return observations

    def _load_plugin_class(self, filepath):
        """
        Loads the plugin method and sets the self.plugin to refer to it.
        :param filepath: path to the file containing the python function
        """
        expected_class = 'ObservationUpdater'
    
        mod_name, file_ext = os.path.splitext(os.path.split(filepath)[-1])

        if file_ext.lower() == '.pyc':
            py_mod = imp.load_compiled(mod_name, filepath)
        else:
            py_mod = imp.load_source(mod_name, filepath)
    
        if hasattr(py_mod, expected_class):
            self.plugin = getattr(py_mod, expected_class)()
        else:
            raise Exception(
                'Cannot find ObservationUpdater class in pluging file ' + filepath)
        
        if not hasattr(self.plugin, 'update'):
            raise Exception('Cannot find update method in plugin class ' + filepath)

    def get_observation(self, collection, observation_id):
        """
        Get an observation from the CAOM2 repo
        :param collection: name of the collection
        :param observation_id: the ID of the observation
        :return: the caom2.observation.Observation object
        """
        assert collection is not None
        assert observation_id is not None
        path = '/{}/{}'.format(collection, observation_id)
        logging.debug('GET '.format(path))

        response = self._repo_client.get((CAOM2REPO_OBS_CAPABILITY_ID, path))
        obs_reader = ObservationReader()
        content = response.content
        if len(content) == 0:
            logging.error(response.status_code)
            response.close()
            raise Exception('Got empty response for resource: {}'.format(path))
        return obs_reader.read(StringIO(content))

    def post_observation(self, observation):
        """
        Updates an observation in the CAOM2 repo
        :param observation: observation to update
        :return: updated observation
        """
        assert observation.collection is not None
        assert observation.observation_id is not None
        path = '/{}/{}'.format(observation.collection, observation.observation_id)
        logging.debug('POST {}'.format(path))

        ibuffer = StringIO()
        ObservationWriter().write(observation, ibuffer)
        obs_xml = ibuffer.getvalue()
        headers = {'Content-Type': 'application/xml'}
        response = self._repo_client.post(
            (CAOM2REPO_OBS_CAPABILITY_ID, path), headers=headers, data=obs_xml)
        logging.debug('Successfully updated Observation\n')

    def put_observation(self, observation):
        """
        Add an observation to the CAOM2 repo
        :param observation: observation to add to the CAOM2 repo
        :return: Added observation
        """
        assert observation.collection is not None
        assert observation.observation_id is not None
        path = '/{}/{}'.format(observation.collection, observation.observation_id)
        logging.debug('PUT {}'.format(path))

        ibuffer = StringIO()
        ObservationWriter().write(observation, ibuffer)
        obs_xml = ibuffer.getvalue()
        headers = {'Content-Type': 'application/xml'}
        response = self._repo_client.put(
            (CAOM2REPO_OBS_CAPABILITY_ID, path), headers=headers, data=obs_xml)
        logging.debug('Successfully put Observation\n')

    def delete_observation(self, collection, observation_id):
        """
        Delete an observation from the CAOM2 repo
        :param collection: Name of the collection
        :param observation_id: ID of the observation
        """
        assert observation_id is not None
        path = '/{}/{}'.format(collection, observation_id)
        logging.debug('DELETE {}'.format(path))
        response = self._repo_client.delete((CAOM2REPO_OBS_CAPABILITY_ID, path))
        logging.info('Successfully deleted Observation {}\n')


def str2date(s):
    """
    Takes a date formatted string and returns a datetime.

    """
    date_format = '%Y-%m-%dT%H:%M:%S'
    if s is None:
        return None
    return datetime.strptime(s, date_format)

def main_app():


    parser = util.get_base_parser(version=version.version, default_resource_id=DEFAULT_RESOURCE_ID)

    parser.description = ('Client for a CAOM2 repo. In addition to CRUD (Create, Read, Update and Delete) '
                          'operations it also implements a visitor operation that allows for updating '
                          'multiple observations in a collection')
    parser.formatter_class = argparse.RawTextHelpFormatter

    subparsers = parser.add_subparsers(dest='cmd')
    create_parser = subparsers.add_parser('create', description='Create a new observation',
                                          help='Create a new observation')
    create_parser.add_argument('observation', help='XML file containing the observation',
                               type=argparse.FileType('r'))

    read_parser = subparsers.add_parser('read',
                                        description='Read an existing observation',
                                        help='Read an existing observation')
    read_parser.add_argument('--output', '-o', help='destination file', required=False)
    read_parser.add_argument('collection', help='collection name in CAOM2 repo')
    read_parser.add_argument('observationID', help='observation identifier')

    update_parser = subparsers.add_parser('update',
                                          description='Update an existing observation',
                                          help='Update an existing observation')
    update_parser.add_argument('observation', help='XML file containing the observation',
                               type=argparse.FileType('r'))

    delete_parser = subparsers.add_parser('delete',
                                          description='Delete an existing observation',
                                          help='Delete an existing observation')
    delete_parser.add_argument('collection', help='collection name in CAOM2 repo')
    delete_parser.add_argument('observationID', help='observation identifier')

    # Note: RawTextHelpFormatter allows for the use of newline in epilog
    visit_parser = subparsers.add_parser('visit',
                                         formatter_class=argparse.RawTextHelpFormatter,
                                         description='Visit observations in a collection',
                                         help='Visit observations in a collection')
    visit_parser.add_argument('--plugin', required=True, type=argparse.FileType('r'),
                              help='plugin class to update each observation')
    visit_parser.add_argument('--start', type=str2date,

                        help='earliest observation to visit (UTC IVOA format: YYYY-mm-ddTH:M:S)')
    visit_parser.add_argument('--end', type=str2date,
                        help='latest observation to visit (UTC IVOA format: YYYY-mm-ddTH:M:S)')
    visit_parser.add_argument('--halt-on-error', action='store_true',
                              help='stop visitor on first update exception raised by plugin')
    visit_parser.add_argument("-s", "--server", help='URL of the CAOM2 repo server')
    visit_parser.add_argument('collection', help='data collection in CAOM2 repo')
    
    visit_parser.epilog =\
"""
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
    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(level=logging.INFO, stream=sys.stdout)
    elif args.debug:
        logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)
    else:
        logging.basicConfig(level=logging.WARN, stream=sys.stdout)

    subject = net.Subject.from_cmd_line_args(args)
    client = CAOM2RepoClient(subject, args.resource_id, host=args.host)
    if args.cmd == 'visit':
        print ("Visit")
        logging.debug("Call visitor with plugin={}, start={}, end={}, collection={}".
                      format(args.plugin.name, args.start, args.end, args.collection))
        (visited, updated, skipped, failed) = \
            client.visit(args.plugin.name, args.collection, start=args.start, end=args.end,
                         halt_on_error=args.halt_on_error)
        logging.info('Visitor stats: visited/updated/skipped/errors: {}/{}/{}/{}'.format(
                     len(visited), len(updated), len(skipped), len(failed)))

    elif args.cmd == 'create':
        logging.info("Create")
        obs_reader = ObservationReader()
        client.put_observation(obs_reader.read(args.observation))
    elif args.cmd == 'read':
        logging.info("Read")
        observation = client.get_observation(args.collection, args.observationID)
        observation_writer = ObservationWriter()
        if args.output:
            with open(args.output, 'w') as obsfile:
                observation_writer.write(observation, obsfile)
        else:
            observation_writer.write(observation, sys.stdout)
    elif args.cmd == 'update':
        logging.info("Update")
        obs_reader = ObservationReader()
        # TODO not sure if need to read in string first
        client.post_observation(obs_reader.read(args.observation))
    else:
        logging.info("Delete")
        client.delete_observation(collection=args.collection, observation_id=args.observationID)

    logging.info("DONE")

if __name__ == '__main__':
    main_app()
