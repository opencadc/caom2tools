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
import sys
print('\n'.join(sys.path))
from datetime import datetime
import errno
import logging
import httplib
import imp
import os
import os.path
#TODO to be changed to io.StringIO when caom2 is prepared for python3
from StringIO import StringIO 
import urllib

from cadctools.net import ws

from ..caom2 import ObservationReader, ObservationWriter

BATCH_SIZE = int(10000)
SERVICE_URI = 'ivo://cadc.nrc.ca/caom2repo'
DATE_FORMAT = "%Y-%m-%dT%H:%M:%S.%f" #IVOA dateformat

class CAOM2RepoClient:

    """Class to do CRUD + visitor actions on a CAOM2 collection repo."""

    def __init__(self, plugin, collection=None, start=None, end=None, retries=0, server=None):

        if retries is not None:
            self.retries = int(retries)
        else:
            self.retries = 0

        self.collection_uri = SERVICE_URI + '/' +\
            self.collection
        # repo client to use
        self._repo_client = ws.BaseWsServer(server)
        logging.info('Service URL: {}'.format(self._repo_client.SERVICE_URL))


    def visit(self, plugin, collection=None, start=None, end=None):
        """
        Main processing function that iterates through the observations of
        the collection and updates them according to the algorithm 
        of the plugin function
        """
        if not os.path.isfile(plugin):
            raise Exception('Cannot find plugin file ' + plugin)
        if start is not None:
            assert type(start) is datetime
        if end is not None:
            assert type(end) is datetime
        self._load_plugin_class(plugin)
        count = 0
        observations = self._get_observations(collection, start, end)
        while len(observations) > 0:
            for observationID in observations:
                observation = self._get_observation(observationID)
                logging.info("Process observation: " + observation.observation_id)
                self.plugin.update(observation)
                self._persist_observation(observationID, observation)
                count = count + 1
            if len(observations) == BATCH_SIZE:
                observations = self._get_observations()
            else:
                # the last batch was smaller so it must have been the last
                break;
        return count

    def _get_observations(self, collection, start=None, end=None):
        """
        Returns a list of datasets from the collection
        """
        assert collection is not None
        observations = []
        params = {'MAXREC':BATCH_SIZE}
        if start:
            params['START'] = start.strftime(DATE_FORMAT)
        if end:
            params['END'] = end.strftime(DATE_FORMAT)

        #self._repo_client.get(collection)
        response = self._repo_client.send_request("GET", "/{}".format(self.collection),
                {"Content-type": "application/x-www-form-urlencoded"},
            urllib.urlencode(params)) #TODO use BaseWs
        if response.status != 200:
            logging.error('Cannot find the %s collection at %s. \n' %
                          (self.collection, self._repo_client.SERVICE_URL))
            raise IOError(errno.ENOENT)
        elif response.status >= 400:
            logging.error(
                'Unable to retrieve list of observations from \'%s\'.\n' % 
                self.collection_uri)
            logging.error('Server Returned: ' +\
                          httplib.responses[response.status] + ' ('
                          + str(response.status) + ')\n' + response.read())
            raise IOError(errno.ENOEXEC)
        else:
            data = response.read()
            last_datetime = None
            for line in data.splitlines():
                (obs, last_datetime) = line.split(',')
                observations.append(obs)
        if last_datetime is not None:
            self.current_start = datetime.strptime(last_datetime, DATE_FORMAT)
        return observations
        
         
    def _load_plugin_class(self, filepath):
        expected_class = 'ObservationUpdater'
    
        mod_name,file_ext = os.path.splitext(os.path.split(filepath)[-1])
    
        if file_ext.lower() == '.py':
            py_mod = imp.load_source(mod_name, filepath)
    
        elif file_ext.lower() == '.pyc':
            py_mod = imp.load_compiled(mod_name, filepath)
    
        if hasattr(py_mod, expected_class):
            self.plugin = getattr(py_mod, expected_class)()
        else:
            raise Exception(
                'Cannot find ObservationUpdater class in pluging file ' +\
                filepath)
        
        if not hasattr(self.plugin, 'update'):
            raise Exception('Cannot find update method in plugin class ' +\
                filepath)
            
    
    def _get_observation(self, observationID):
        logging.debug("GET " + observationID)
        uri = '/' + self.collection + '/' + observationID
        response = self._repo_client.send_request("GET", uri, {}, '')

        logging.debug("GET status: " + str(response.status))
        if response.status == 503 and (self.retries > 0):
            response = self._repo_client.retry("GET", uri, {}, '')

        if response.status == 404:
            logging.error('No such Observation found with URI \'%s\'.\n' % uri)
            raise IOError(errno.ENOENT)
        elif response.status >= 400:
            logging.error('Unable to retrieve Observation with URI \'%s\'.\n' % uri)
            logging.error('Server Returned: ' +\
                          httplib.responses[response.status] + ' ('
                          + str(response.status) + ')\n' + response.read())
            raise IOError(errno.ENOEXEC)
        else:
            obs_reader = ObservationReader()
            content = response.read()
            if len(content) == 0:
                logging.error(response.status)
                response.close()
                raise Exception("Got empty response for uri: " + uri)
            return obs_reader.read(StringIO(content))
        
    def _persist_observation(self, observationID, observation):
        uri = '/' + self.collection + '/' + observationID
        ibuffer = StringIO()
        ObservationWriter().write(observation, ibuffer)
        obs_xml = ibuffer.getvalue()
        response = self._repo_client.send_request(
            "POST", uri, {'Content-Type': 'text/xml'}, obs_xml)

        if response.status == 503 and (self.retries > 0):
            response = self._repo_client.retry(
            "POST", uri, {'Content-Type': 'text/xml'}, obs_xml)

        if response.status == 404:
            logging.error('Observation with URI \'%s\' does not exist.\n' % uri)
            raise IOError(errno.ENOENT)
        elif response.status >= 400:
            msg = ''
            for hmsg in response.msg.headers:
                msg = msg + hmsg
            logging.error('Unable to update Observation ' + str(observationID)
                          + '\nServer Returned: ' +\
                        httplib.responses[response.status] + ' ('
                          + str(response.status) + ')\n' + msg + response.read())
            raise IOError(errno.ENOEXEC)
        else:
            logging.debug('Successfully updated Observation\n')


def run():

    parser = CommonParser(description=description,
                          formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--plugin', required=True,
                        metavar=('<pluginClassFile>'),
                        help='Pluging class to update each observation')
    parser.add_argument('--start', metavar='<datetime start point>',
                        type=parse_date,
                        help='oldest dataset to visit (UTC %%Y-%%m-%%d format)');
    parser.add_argument('--end', metavar='<datetime end point>',
                        type=parse_date,
                        help='earliest dataset to visit (UTC %%Y-%%m-%%d format)');
    parser.add_argument('--retries', metavar='<number of retries>',
                        type=int,
                        help='number of tries with transient server errors');
    parser.add_argument("-s", "--server", metavar=('<CAOM2 service URL>'),
                      help="URL of the CAOM2 repo server")

    parser.add_argument('collection', metavar='<datacollection>', type=str,
                help='data collection in CAOM2 repo')
    parser.epilog =\
"""
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

    args = parser.parse_args()


    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    plugin = args.plugin
    if not os.path.isfile(args.plugin):
        logging.error('Cannot find plugin file ' + plugin)
        sys.exit(-1)

    start = args.start
    end = args.end
    retries = args.retries
    collection = args.collection

    logging.debug("Call visitor with plugin={}, start={}, end={}, dataset={}".
           format(plugin, start, end, collection, retries))
    client = CAOM2Visitor(plugin=plugin, start=start, end=end,
                          collection=collection, retries=retries, server=args.server)
    client.process()
    logging.debug("debug")

    print("DONE")


