#
# -*- coding: utf-8 -*-
"""Client script to access the CAOM-2 repository Observations."""

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

__author__ = 'jenkinsd'

from util.caom2_properties import RegistryProperties
#from caom2.caom2_exceptions import ObservationParsingException
from urlparse import urlparse
from httplib import HTTPSConnection, HTTPConnection
from argparse import ArgumentParser
import httplib
import netrc
import sys
import base64
import os
import errno
import logging

class CAOM2RepoClient:

    CAOM2REPO_REGISTRY_PROP_DIR_VAR = "CADC_ROOT"

    #
    # Main function for execution.  This function will delegate to proper functions.
    #
    def main(self):
        # Some constants.
        parser = ArgumentParser(description="CAOM-2 Observation Repository client.")

        parser.add_argument('--version', action='version', version='%(prog)s 1.0')
        parser.add_argument('-v', '--verbose', required=False, action="store_true",
                            dest="verbose")
        parser.add_argument('-d', '--debug', required=False, action="store_true",
                            dest="debug")
        parser.add_argument('-g', '--get', required=False, dest="get_action",
                            nargs=2, metavar=("<observationURI>", "<filename>"),
                            help="Get observation <observationURI> (in the form caom:<collection>/<observationID>) from the repository and write to <filename>")
        parser.add_argument('-p', '--put', required=False, dest='create_action',
                            nargs=2, metavar=("<observationURI>", "<filename>"),
                            help="Create observation <observationURI> (in the form caom:<collection>/<observationID>) in the repository using <filename>")
        parser.add_argument('-u', '--update', required=False, dest='update_action',
                            nargs=2, metavar=("<observationURI>", "<filename>"),
                            help="Update observation <observationURI> (in the form caom:<collection>/<observationID>) in the repository using <filename>")
        parser.add_argument('-r', '--remove', required=False, dest='delete_action',
                            nargs=1, metavar="<observationURI>",
                            help="Remove observation <observationURI> (in the form caom:<collection>/<observationID>) from the repository")
#        parser.add_argument('observationURI', metavar='<observationURI>', nargs=1, help="A CAOM-2.0 Observation URI of the form caom:<collection>/<observationID>")
#        parser.add_argument('filename', metavar='<filename>', nargs=1, help="CAOM-2.0 XML document (required for get, put, update)")

        arguments = parser.parse_args(sys.argv[1:])

        if arguments.verbose:
            logging.basicConfig(level=logging.INFO)

        if arguments.debug:
            logging.basicConfig(level=logging.DEBUG)

        # CAOM-2 Repository registry lookup.
        self.CAOM2_SERVICE_URI = 'ivo://cadc.nrc.ca/caom2repo'
        if CAOM2RepoClient.CAOM2REPO_REGISTRY_PROP_DIR_VAR in os.environ \
           and os.path.isfile(os.environ[CAOM2RepoClient.CAOM2REPO_REGISTRY_PROP_DIR_VAR]
                              + '/config/' + RegistryProperties.FILE_NAME):
            registryClientPropertiesDir= \
                os.environ[CAOM2RepoClient.CAOM2REPO_REGISTRY_PROP_DIR_VAR] \
                + '/config'
        else:
            # Otherwise, read from the EGG file.
            registryClientPropertiesDir = None


        self.SERVICE_PROTOCOL = 'https'
        self.REGISTRY = RegistryProperties(registryClientPropertiesDir)
        self.SERVICE_URL = self.REGISTRY.get_service_url(self.CAOM2_SERVICE_URI,
                                                         self.SERVICE_PROTOCOL)

        logging.info("Found Service URL: '%s'" % self.SERVICE_URL)

        if arguments.get_action:
            logging.info("GET ACTION")
            self.get(arguments.get_action[0], arguments.get_action[1])
        elif arguments.create_action:
            logging.info("PUT ACTION")
            self.put(arguments.create_action[0], arguments.create_action[1])
        elif arguments.update_action:
            logging.info("POST ACTION")
            self.update(arguments.update_action[0], arguments.update_action[1])
        elif arguments.delete_action:
            logging.info("DELETE ACTION")
            self.remove(arguments.delete_action[0])
        else:
            parser.print_help()


    #
    # Obtain the CAOM-2 Observation instance from the CAOM-2 service.
    #
    # @param    observationURI - The URI of the CAOM-2 Observation
    # @param    filename - The full path to the XML File to construct an
    #                      Observation.
    #
    def get(self, observationURI, filename):
        logging.info("GETting " + observationURI)
        observationResponse = self.send_request("GET", observationURI, {}, '')

        status = observationResponse.status

        if status == 404:
            logging.error('No such Observation found with URI \'%s\'.\n' % observationURI)
            sys.exit(errno.ENOENT)
        elif status > 404:
            logging.error('Unable to retrieve Observation with URI \'%s\'.\n' % observationURI)
            logging.error('Server Returned: ' + httplib.responses[status] + ' ('\
                          + str(status) + ')\n' + observationResponse.read())
            sys.exit(errno.ENOEXEC)
        else:
            f = open(filename, 'w')
            f.write(observationResponse.read())
            f.close()
            logging.info("Successfully saved Observation at '%s'" % filename)
    #
    # Create a new Observation resource.
    #
    # @param    observationURI - The URI of the CAOM-2 Observation
    # @param    filename - The full path to the XML File to construct an
    #                      Observation.
    #
    def put(self, observationURI, filename):
        logging.debug("PUTting " + filename)
        xmlfile = None

        try:

            xmlfile = open(filename)
            response = self.send_request("PUT", observationURI,
                                         {'Content-Type': 'text/xml'},
                                         xmlfile.read())
            status = response.status

            if status == 404:
                logging.error('No such collection found for URI \'%s\'' % observationURI)
                sys.exit(errno.ENOENT)
            elif status >= 300:
                msg = ''
                for hmsg in response.msg.headers:
                    msg = msg + hmsg
                logging.error('Unable to create Observation from file ' + filename \
                             + '\nServer Returned: ' + httplib.responses[status] + ' (' \
                             + str(status) + ')\n' + msg + response.read())
                sys.exit(errno.ENOEXEC)
            else:
                logging.info('Successfully created Observation\n')

        #except ObservationParsingException, err:
        #    logging.error('\nAborting due to error!\nUnable to parse an Observation from '\
        #                  + filename + '\n' + str(err) + '\n')
        #    sys.exit(errno.EIO)
        except IOError as ioerr:
            logging.error('\nAborting due to error!\nUnable to read file '\
                         + filename + '\n' + str(ioerr) + '\n')
            sys.exit(errno.EIO)
        finally:
            if not xmlfile is None:
                xmlfile.close()

    #
    # Update an existing CAOM-2 Observation.
    #
    # @param    observationURI  - The URI of the Observation to update.
    # @param    filename        - The full path to the XML File of the Observation.
    #
    def update(self, observationURI, filename):
        logging.debug("POSTing " + observationURI + " with filename " + filename)
        xmlfile = None

        try:
            xmlfile = open(filename)
            response = self.send_request("POST", observationURI,
                                         {'Content-Type': 'text/xml'},
                                         xmlfile.read())
            status = response.status

            if status == 404:
                logging.error('Observation with URI \'%s\' does not exist.\n' % observationURI)
                sys.exit(errno.ENOENT)
            elif status >= 300:
                logging.error('Unable to update Observation from file ' + filename\
                              + '\nServer Returned: ' + httplib.responses[status] + ' ('\
                              + str(status) + ')\n')
                sys.exit(errno.ENOEXEC)
            else:
                logging.info('Successfully updated Observation\n')

        #except ObservationParsingException, err:
        #    logging.error('Aborting due to error!\nUnable to parse an Observation from '\
        #                  + filename + '\n' + str(err) + '\n')
        #    sys.exit(errno.EIO)
        except IOError as ioerr:
            logging.error('Aborting due to error!\nUnable to read file '\
                          + filename + '\n' + str(ioerr) + '\n')
            sys.exit(errno.EIO)
        finally:
            if not xmlfile is None:
                xmlfile.close()
    #
    # Permanently remove an Observation resource.
    #
    # @param    observationURI - The URI of the Observation to delete.
    #
    def remove(self, observationURI):
        logging.debug("DELETing " + observationURI)
        response = self.send_request("DELETE", observationURI, {}, '')
        status = response.status

        if status == 404:
            logging.error('No such Observation found with URI \'%s\'.\n' % observationURI)
            sys.exit(errno.ENOENT)
        elif status >= 300:
            logging.error('Unable to remove Observation with URI \'' + observationURI \
                  + '\'.\n\n' + httplib.responses[status] + ' (' + str(status) \
                  + ')\n' + response.read())
            sys.exit(errno.ENOEXEC)
        else:
            logging.info('Successfully removed Observation %s' % observationURI + '\n')

    #
    # Send the HTTP(S) request.
    #
    # @param    method      - The HTTP Method to use (String).
    # @param    path        - The path of the URL (After the host).  Should begin
    #                         with a slash ('/') (String).
    # @param    headers     - Any custom headers.  This is a dictionary.
    # @param    payload     - The payload to send for a write operation (String).
    #
    # @return      The httplib.HTTPResponse object.
    #
    def send_request(self, method, observationURI, headers, payload):
        parseResult = urlparse(observationURI)
        serviceURLResult = urlparse(self.SERVICE_URL)
        path = parseResult.path
        logging.debug("Found path: " + path)

        if self.SERVICE_PROTOCOL == 'https':
            try:
                with open(os.path.join(os.environ['HOME'], '.ssl/cadcproxy.pem')) \
                                                                            as certfile:
                    conn = HTTPSConnection(serviceURLResult.hostname, 443,
                                           None, certfile.name)
                    conn.request(method, serviceURLResult.path + '/pub/' + path,
                                 payload, headers)
                    logging.debug("Making request to " + self.SERVICE_URL
                                  + '/pub/' + path)
                    return conn.getresponse()
            except IOError as e:
                logging.error('No usable credentials to connect to '
                              + self.SERVICE_URL + '/' + path + '\n')
                logging.error(str(e) + "\n")
                sys.exit(errno.EACCES)
        elif self.SERVICE_PROTOCOL == 'http':
            try:
                netrcfile = netrc.netrc()
                auth = netrcfile.authenticators(serviceURLResult.hostname)
            except netrc.NetrcParseError as err:
                logging.error('Unable to read netrc file.\n')
                logging.error(str(err) + "\n")
                sys.exit(errno.EIO)

            conn = HTTPConnection(serviceURLResult.hostname)
            username = auth[0]
            password = auth[2]
            base64string = base64.encodestring('%s:%s' % (username, password)).replace('\n', '')

            if not headers or not len(headers):
                headers = {'Authorization': 'Basic %s' % base64string}
            elif not headers.has_key('Authorization'):
                headers.update({'Authorization': 'Basic %s' % base64string})

            conn.request(method, serviceURLResult.path + '/auth/' + path, payload, headers)
            logging.debug('Making request to ' + self.SERVICE_URL + '/auth/'
                          + path)
            return conn.getresponse()
        else:
            logging.error("Unknown protocol '%s'" % self.SERVICE_PROTOCOL)

