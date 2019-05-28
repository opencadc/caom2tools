# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2018.                            (c) 2018.
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

import csv
import logging
import os
import requests
import subprocess
import yaml

from datetime import datetime
from enum import Enum
from hashlib import md5
from io import BytesIO
from os import stat
from requests.adapters import HTTPAdapter
from urllib import parse as parse
from urllib3 import Retry

from cadcutils import net
from cadcdata import CadcDataClient
from caom2 import ObservationWriter, ObservationReader, Artifact
from caom2 import ChecksumURI


__all__ = ['CadcException', 'Config', 'State', 'to_float', 'TaskType',
           'exec_cmd', 'exec_cmd_redirect', 'exec_cmd_info',
           'get_cadc_meta', 'get_file_meta', 'compare_checksum',
           'decompose_lineage', 'check_param', 'read_csv_file',
           'write_obs_to_file', 'read_obs_from_file',
           'compare_checksum_client', 'Features', 'write_to_file',
           'read_from_file', 'read_file_list_from_archive', 'update_typed_set',
           'get_cadc_headers', 'get_lineage', 'get_artifact_metadata',
           'data_put', 'data_get', 'build_uri']


class CadcException(Exception):
    """Generic exception raised by failure cases within the caom2pipe
    module."""
    pass


class Features(object):
    """Boolean feature flag implementation."""

    def __init__(self):
        self.use_file_names = True
        self.run_in_airflow = True
        self.supports_composite = True
        self.supports_catalog = True
        self.expects_retry = True

    @property
    def use_file_names(self):
        """If true, the lists of work to be done are expected to
        identify file names. If false, they are expected to identify
        observation IDs."""
        return self._use_file_names

    @use_file_names.setter
    def use_file_names(self, value):
        self._use_file_names = value

    @property
    def run_in_airflow(self):
        """If true, will treat command-line arguments as if the application
        is running in airflow."""
        return self._run_in_airflow

    @run_in_airflow.setter
    def run_in_airflow(self, value):
        self._run_in_airflow = value

    @property
    def supports_composite(self):
        """If true, will execute any specific code for composite observation
        definition."""
        return self._supports_composite

    @supports_composite.setter
    def supports_composite(self, value):
        self._supports_composite = value

    @property
    def supports_catalog(self):
        """If true, will execute any specific code for catalog handling
        when creating a CAOM instance."""
        return self._supports_catalog

    @supports_catalog.setter
    def supports_catalog(self, value):
        self._supports_catalog = value

    @property
    def expects_retry(self):
        """If true, will execute any specific code for running retries
        based on retries_log.txt content."""
        return self._expects_retry

    @expects_retry.setter
    def expects_retry(self, value):
        self._expects_retry = value

    def __str__(self):
        return ' '.join(
            '{} {}'.format(ii, getattr(self, ii)) for ii in vars(self))


class TaskType(Enum):
    """The possible steps in a Collection pipeline. A short-hand, user-facing
    way to identify the work to be done by a pipeline."""
    STORE = 'store'  # store a local file to ad
    SCRAPE = 'scrape'  # local CAOM instance creation, no network required
    INGEST = 'ingest'  # create a CAOM instance from metadata only
    MODIFY = 'modify'  # modify a CAOM instance from data
    CHECKSUM = 'checksum'  # is the checksum on local disk the same as in ad?
    VISIT = 'visit'    # visit an observation
    # remote file storage, create CAOM instance via local metadata
    REMOTE = 'remote'
    # retrieve file via HTTP to local temp storage, store to ad
    PULL = 'pull'


class State(object):
    """Persist information between pipeline invocations.

    Currently the State class persists the concept of a bookmark, which is the
    place in the flow of data that was last processed. This 'place' may be a
    timestamp, or an id. That value is up to clients of this class.
    """

    def __init__(self, fqn):
        self.fqn = fqn
        self.bookmarks = {}
        self.logger = logging.getLogger('State')
        result = read_as_yaml(self.fqn)
        if result is None:
            raise CadcException('Could not load state from {}'.format(fqn))
        else:
            self.bookmarks = result.get('bookmarks')
            self.content = result

    def get_bookmark(self, key):
        """Lookup for last_record key."""
        result = None
        if key in self.bookmarks:
            if 'last_record' in self.bookmarks[key]:
                result = self.bookmarks[key]['last_record']
            else:
                self.logger.warning('No record found for {}'.format(key))
        else:
            self.logger.warning('No bookmarks found for {}'.format(key))
        return result

    def save_state(self, key, value):
        """Write the current state as a YAML file.
        :param key which record is being updated
        :param value the value to update the record with
        """
        if key in self.bookmarks:
            if 'last_record' in self.bookmarks[key]:
                self.bookmarks[key]['last_record'] = value
                write_as_yaml(self.content, self.fqn)
            else:
                self.logger.warning('No record found for {}'.format(key))
        else:
            self.logger.warning('No bookmarks found for {}'.format(key))


class Config(object):
    """Configuration information that remains the same for all steps and all
     work in a pipeline execution."""

    def __init__(self):
        self.working_directory = None
        self.work_file = None
        # the fully qualified name for the work file
        self.work_fqn = None
        self.netrc_file = None
        self.archive = None
        self.collection = None
        self.use_local_files = False
        self.resource_id = None
        self.tap_id = None
        self.logging_level = None
        self.log_to_file = False
        self.log_file_directory = None
        self.stream = None
        self.storage_host = None
        self.task_types = None
        self.success_log_file_name = None
        # the fully qualified name for the file
        self.success_fqn = None
        self.failure_log_file_name = None
        # the fully qualified name for the file
        self.failure_fqn = None
        self.retry_file_name = None
        # the fully qualified name for the file
        self.retry_fqn = None
        self.retry_failures = False
        self.retry_count = 1
        self.proxy_file_name = None
        # the fully qualified name for the file
        self.proxy_fqn = None
        self.state_file_name = None
        # the fully qualified name for the file
        self.state_fqn = None
        self.features = Features()

    @property
    def working_directory(self):
        """the root directory for all executor operations"""
        return self._working_directory

    @working_directory.setter
    def working_directory(self, value):
        self._working_directory = value

    @property
    def work_file(self):
        """ the file that contains the list of work to be passed through the
        pipeline"""
        return self._work_file

    @work_file.setter
    def work_file(self, value):
        self._work_file = value
        if self.working_directory is not None:
            self.work_fqn = os.path.join(
                self.working_directory, self.work_file)

    @property
    def netrc_file(self):
        """credentials for any service calls"""
        return self._netrc_file

    @netrc_file.setter
    def netrc_file(self, value):
        self._netrc_file = value

    @property
    def collection(self):
        """which collection is addressed by the pipeline"""
        return self._collection

    @collection.setter
    def collection(self, value):
        self._collection = value

    @property
    def archive(self):
        """which archive is addressed by the pipeline"""
        return self._archive

    @archive.setter
    def archive(self, value):
        self._archive = value

    @property
    def use_local_files(self):
        """changes expectations of the executors for handling files on disk"""
        return self._use_local_files

    @use_local_files.setter
    def use_local_files(self, value):
        self._use_local_files = value

    @property
    def resource_id(self):
        """which service instance to use"""
        return self._resource_id

    @resource_id.setter
    def resource_id(self, value):
        self._resource_id = value

    @property
    def tap_id(self):
        """which tap service instance to use"""
        return self._tap_id

    @tap_id.setter
    def tap_id(self, value):
        self._tap_id = value

    @property
    def log_to_file(self):
        """boolean - write the log to a file?"""
        return self._log_to_file

    @log_to_file.setter
    def log_to_file(self, value):
        self._log_to_file = value

    @property
    def log_file_directory(self):
        """where log files are written to - defaults to working_directory"""
        return self._log_file_directory

    @log_file_directory.setter
    def log_file_directory(self, value):
        self._log_file_directory = value

    @property
    def logging_level(self):
        """the logging level - enforced throughout the pipeline"""
        return self._logging_level

    @logging_level.setter
    def logging_level(self, value):
        lookup = {'DEBUG': logging.DEBUG,
                  'INFO': logging.INFO,
                  'WARNING': logging.WARNING,
                  'ERROR': logging.ERROR}
        if value in lookup:
            self._logging_level = lookup[value]

    @property
    def stream(self):
        """the ad 'stream' that goes with the archive - use when storing
        files"""
        return self._stream

    @stream.setter
    def stream(self, value):
        self._stream = value

    @property
    def storage_host(self):
        """the ad 'host' to store files to - used for testing cadc-data put
        commands only, should usually be None"""
        return self._storage_host

    @storage_host.setter
    def storage_host(self, value):
        self._storage_host = value

    @property
    def task_type(self):
        """the way to control which steps get executed"""
        return self._task_type

    @task_type.setter
    def task_type(self, value):
        self._task_type = value

    @property
    def success_log_file_name(self):
        """the filename where success logs are written, this will be created
        in log_file_directory"""
        return self._success_log_file_name

    @success_log_file_name.setter
    def success_log_file_name(self, value):
        self._success_log_file_name = value
        if self.log_file_directory is not None:
            self.success_fqn = os.path.join(
                self.log_file_directory, self.success_log_file_name)

    @property
    def failure_log_file_name(self):
        """the filename where failure logs are written this will be created
        in log_file_directory"""
        return self._failure_log_file_name

    @failure_log_file_name.setter
    def failure_log_file_name(self, value):
        self._failure_log_file_name = value
        if self.log_file_directory is not None:
            self.failure_fqn = os.path.join(
                self.log_file_directory, self.failure_log_file_name)

    @property
    def retry_file_name(self):
        """the filename where retry entries are written this will be created
        in log_file_directory"""
        return self._retry_file_name

    @retry_file_name.setter
    def retry_file_name(self, value):
        self._retry_file_name = value
        if self.log_file_directory is not None:
            self.retry_fqn = os.path.join(
                self.log_file_directory, self.retry_file_name)

    @property
    def retry_failures(self):
        """Will the application retry the entries in the
        retries.txt file? If True, the application will attempt to re-run
        the work to do for each entry in the retries.txt file. If False,
        it will do nothing."""
        return self._retry_failures

    @retry_failures.setter
    def retry_failures(self, value):
        self._retry_failures = value

    @property
    def retry_count(self):
        """how many times the application will retry the entries in the
        retries.txt file."""
        return self._retry_count

    @retry_count.setter
    def retry_count(self, value):
        self._retry_count = value

    @property
    def proxy_file_name(self):
        """If using a proxy certificate for authentication, identify the
        fully-qualified pathname here."""
        return self._proxy_file_name

    @proxy_file_name.setter
    def proxy_file_name(self, value):
        self._proxy_file_name = value
        if (self.working_directory is not None and
                self.proxy_file_name is not None):
            self.proxy_fqn = os.path.join(
                self.working_directory, self.proxy_file_name)

    @property
    def state_file_name(self):
        """If using a state file to communicate persistent information between
        invocations, identify the fully-qualified pathname here."""
        return self._state_file_name

    @state_file_name.setter
    def state_file_name(self, value):
        self._state_file_name = value
        if (self.working_directory is not None and
                self.state_file_name is not None):
            self.state_fqn = os.path.join(
                self.working_directory, self.state_file_name)

    @property
    def features(self):
        """Feature flag setting access."""
        return self._features

    @features.setter
    def features(self, value):
        self._features = value

    @staticmethod
    def _lookup(config, lookup, default):
        if lookup in config:
            result = config[lookup]
        else:
            result = default
        return result

    def __str__(self):
        return 'working_directory:: \'{}\' ' \
               'work_fqn:: \'{}\' ' \
               'netrc_file:: \'{}\' ' \
               'archive:: \'{}\' ' \
               'collection:: \'{}\' ' \
               'task_types:: \'{}\' ' \
               'stream:: \'{}\' ' \
               'resource_id:: \'{}\' ' \
               'tap_id:: \'{}\' ' \
               'use_local_files:: \'{}\' ' \
               'log_to_file:: \'{}\' ' \
               'log_file_directory:: \'{}\' ' \
               'success_log_file_name:: \'{}\' ' \
               'success_fqn:: \'{}\' ' \
               'failure_log_file_name:: \'{}\' ' \
               'failure_fqn:: \'{}\' ' \
               'retry_file_name:: \'{}\' ' \
               'retry_fqn:: \'{}\' ' \
               'retry_failures:: \'{}\' ' \
               'retry_count:: \'{}\' ' \
               'proxy_file:: \'{}\' ' \
               'state_fqn:: \'{}\' ' \
               'features:: \'{}\' ' \
               'logging_level:: \'{}\''.format(
                self.working_directory, self.work_fqn, self.netrc_file,
                self.archive, self.collection, self.task_types, self.stream,
                self.resource_id,
                self.tap_id, self.use_local_files, self.log_to_file,
                self.log_file_directory, self.success_log_file_name,
                self.success_fqn, self.failure_log_file_name,
                self.failure_fqn, self.retry_file_name, self.retry_fqn,
                self.retry_failures, self.retry_count, self.proxy_fqn,
                self.state_fqn, self.features, self.logging_level)

    @staticmethod
    def _obtain_task_types(config, default=None):
        """Make the configuration file entries into the Enum."""
        task_types = []
        if 'task_types' in config:
            for ii in config['task_types']:
                task_types.append(TaskType(ii))
            return task_types
        else:
            return default

    @staticmethod
    def _obtain_features(config):
        """Make the configuration file entries into the class members."""
        feature_flags = Features()
        if 'features' in config:
            for ii in config['features']:
                if not config['features'][ii]:
                    getattr(feature_flags, ii)
                    setattr(feature_flags, ii, False)
        return feature_flags

    def get(self):
        """Look up the configuration values in the data structure extracted
        from the configuration file."""
        return self.get_executors()

    def get_executors(self):
        """Look up the configuration values in the data structure extracted
        from the configuration file.

        Consider this deprecated - use get instead, because the name is
        non-representative of the work being done.
        """
        try:
            config = self.get_config()
            self.working_directory = \
                self._lookup(config, 'working_directory', os.getcwd())
            self.work_file = self._lookup(config, 'todo_file_name', 'todo.txt')
            self.netrc_file = \
                self._lookup(config, 'netrc_filename', 'test_netrc')
            self.resource_id = self._lookup(
                config, 'resource_id', 'ivo://cadc.nrc.ca/sc2repo')
            self.tap_id = self._lookup(
                config, 'tap_id', 'ivo://cadc.nrc.ca/sc2tap')
            self.use_local_files = bool(
                self._lookup(config, 'use_local_files', False))
            self.logging_level = self._lookup(config, 'logging_level', 'DEBUG')
            self.log_to_file = self._lookup(config, 'log_to_file', False)
            self.log_file_directory = self._lookup(
                config, 'log_file_directory', self.working_directory)
            self.stream = self._lookup(config, 'stream', 'raw')
            self.task_types = self._obtain_task_types(
                config, [TaskType.SCRAPE])
            self.collection = self._lookup(config, 'collection', 'TEST')
            self.archive = self._lookup(config, 'archive', self.collection)
            self.success_log_file_name = self._lookup(config,
                                                      'success_log_file_name',
                                                      'success_log.txt')
            self.failure_log_file_name = self._lookup(config,
                                                      'failure_log_file_name',
                                                      'failure_log.txt')
            self.retry_file_name = self._lookup(config, 'retry_file_name',
                                                'retries.txt')
            self.retry_failures = self._lookup(config, 'retry_failures', False)
            self.retry_count = self._lookup(config, 'retry_count', 1)
            self.features = self._obtain_features(config)
            self.proxy_file_name = self._lookup(
                config, 'proxy_file_name', None)
            self.state_file_name = self._lookup(
                config, 'state_file_name', None)
        except KeyError as e:
            raise CadcException(
                'Error in config file {}'.format(e))

    def get_config(self):
        """Return a configuration dictionary. Assumes a file named config.yml
        in the current working directory."""
        config_fqn = os.path.join(os.getcwd(), 'config.yml')
        config = self.load_config(config_fqn)
        if config is None:
            raise CadcException(
                'Could not find the file {}'.format(config_fqn))
        return config

    def need_to_retry(self):
        """Evaluate the need to have the pipeline try to re-execute for any
         files/observations that have been logged as failures.

        If log_to_file is not set to True, there is no retry file content
        to retry on.

         :param config does the configuration identify retry information?
         :return True if the configuration and logging information indicate a
            need to attempt to retry the pipeline execution for any entries.
         """
        result = True
        if (self.features is not None and self.features.expects_retry and
                self.retry_failures and self.log_to_file):
            meta = get_file_meta(self.retry_fqn)
            if meta['size'] == 0:
                logging.info('Checked the retry file {}. There are no logged '
                             'failures.'.format(self.retry_fqn))
                result = False
        else:
            result = False
        return result

    def update_for_retry(self, count):
        """
        When retrying, the application will:

        - use the retries.txt file as the todo list
        - retry as many times as the 'retry count' in the config.yml file.
        - make a new log directory, in the working directory, with the name
            logs_{retry_count}. Any failures for the retry execution that
            need to be logged will be logged here.
        - in the new log directory, make a new .xml file for the
            output, with the name {obs_id}.xml

        :param count the current retry iteration
        """
        self.work_file = '{}'.format(self.retry_file_name)
        self.work_fqn = self.retry_fqn
        if '_' in self.log_file_directory:
            temp = self.log_file_directory.split('_')[0]
            self.log_file_directory = '{}_{}'.format(temp, count)
        else:
            self.log_file_directory = '{}_{}'.format(
                self.log_file_directory, count)
        # reset the location of the log file names
        self.success_log_file_name = self.success_log_file_name
        self.failure_log_file_name = self.failure_log_file_name
        self.retry_file_name = self.retry_file_name

        logging.info('Retry work file is {}'.format(self.work_fqn))

    @staticmethod
    def load_config(config_fqn):
        """Read a configuration as a YAML file.
        :param config_fqn the fully qualified name for the configuration
            file.
        """
        try:
            logging.debug('Begin load_config.')
            with open(config_fqn) as f:
                data_map = yaml.safe_load(f)
                logging.debug('End load_config.')
                return data_map
        except (yaml.scanner.ScannerError, FileNotFoundError) as e:
            logging.error(e)
            return None


def to_float(value):
    """Cast to float, without throwing an exception."""
    return float(value) if value is not None else None


def exec_cmd(cmd):
    """
    This does command execution as a subprocess call.

    :param cmd the text version of the command being executed
    :return None
    """
    logging.debug(cmd)
    cmd_array = cmd.split()
    try:
        child = subprocess.Popen(cmd_array, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
        output, outerr = child.communicate()
        logging.debug('stdout {}'.format(output.decode('utf-8')))
        logging.debug('stderr {}'.format(outerr.decode('utf-8')))
        if child.returncode != 0:
            logging.debug('Command {} failed.'.format(cmd))
            raise CadcException(
                'Command {} had stdout{} stderr {}'.format(
                    cmd, output.decode('utf-8'), outerr.decode('utf-8')))
    except Exception as e:
        logging.debug('Error with command {}:: {}'.format(cmd, e))
        raise CadcException('Could not execute cmd {}. '
                            'Exception {}'.format(cmd, e))


def exec_cmd_info(cmd):
    """
    This does command execution as a subprocess call.

    :param cmd the text version of the command being executed
    :return The text from stdout.
    """
    logging.debug(cmd)
    cmd_array = cmd.split()
    try:
        output, outerr = subprocess.Popen(cmd_array, stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE).communicate()
        if outerr is not None and len(outerr) > 0 and outerr[0] is not None:
            raise CadcException('Command {} had stderr {}'.format(
                cmd, outerr.decode('utf-8')))
        if output is not None and len(output) > 0:
            return output.decode('utf-8')
    except Exception as e:
        logging.debug('Error with command {}:: {}'.format(cmd, e))
        raise CadcException('Could not execute cmd {}.'
                            'Exception {}'.format(cmd, e))


def exec_cmd_redirect(cmd, fqn):
    """
    This does command execution as a subprocess call. It redirects stdout
    to fqn, and assumes binary output for the re-direct.

    :param cmd the text version of the command being executed
    :param fqn the fully-qualified name of the file to which stdout is
        re-directed
    :return None
    """
    logging.debug(cmd)
    cmd_array = cmd.split()
    try:
        with open(fqn, 'wb') as outfile:
            outerr = subprocess.Popen(
                cmd_array, stdout=outfile,
                stderr=subprocess.PIPE).communicate()
            if (outerr is not None and len(outerr) > 0 and
                    outerr[0] is not None):
                logging.debug('Command {} had stderr {}'.format(
                    cmd, outerr.decode('utf-8')))
                raise CadcException(
                    'Command {} had outerr {}'.format(
                        cmd, outerr.decode('utf-8')))
    except Exception as e:
        logging.debug('Error with command {}:: {}'.format(cmd, e))
        raise CadcException('Could not execute cmd {}.'
                            'Exception {}'.format(cmd, e))


def get_cadc_headers(uri):
    """
    Creates the FITS headers object by fetching the FITS headers of a CADC
    file. The function takes advantage of the fhead feature of the CADC
    storage service and retrieves just the headers and no data, minimizing
    the transfer time.

    The file must be public, because the header retrieval is done as an
    anonymous user.

    :param uri: CADC URI
    :return: a string of keyword/value pairs.
    """
    file_url = parse.urlparse(uri)
    # create possible types of subjects
    subject = net.Subject()
    client = CadcDataClient(subject)
    # do a fhead on the file
    archive, file_id = file_url.path.split('/')
    b = BytesIO()
    b.name = uri
    client.get_file(archive, file_id, b, fhead=True)
    fits_header = b.getvalue().decode('ascii')
    b.close()
    return fits_header


def get_cadc_meta(netrc_fqn, archive, fname):
    """
    Gets contentType, contentLength and contentChecksum of a CADC artifact
    :param netrc_fqn: user credentials
    :param archive: archive file has been stored to
    :param fname: name of file in the archive
    :return:
    """
    subject = net.Subject(username=None, certificate=None, netrc=netrc_fqn)
    client = CadcDataClient(subject)
    return client.get_file_info(archive, fname)


def get_file_meta(fqn):
    """
    Gets contentType, contentLength and contentChecksum of an artifact on disk.

    :param fqn: Fully-qualified name of the file for which to get the metadata.
    :return:
    """
    if fqn is None or not os.path.exists(fqn):
        raise CadcException('Could not find {} in get_file_meta'.format(fqn))
    meta = {}
    s = stat(fqn)
    meta['size'] = s.st_size
    meta['md5sum'] = md5(open(fqn, 'rb').read()).hexdigest()
    if fqn.endswith('.header') or fqn.endswith('.txt'):
        meta['type'] = 'text/plain'
    elif fqn.endswith('.csv'):
        meta['type'] = 'text/csv'
    elif fqn.endswith('.gif'):
        meta['type'] = 'image/gif'
    elif fqn.endswith('.png'):
        meta['type'] = 'image/png'
    elif fqn.endswith('.jpg'):
        meta['type'] = 'image/jpeg'
    elif fqn.endswith('tar.gz'):
        meta['type'] = 'application/gzip'
    else:
        meta['type'] = 'application/fits'
    logging.debug(meta)
    return meta


def _check_checksums(fqn, archive, local_meta, ad_meta):
    """Raise CadcException if the checksum of a file in ad is not the same as
    the checksum of a file on disk.

    :param fqn: Fully-qualified name of file for which to compare metadata.
    :param archive: archive file has been stored to
    :param local_meta: md5 checksum for the file on disk
    :param ad_meta: md5 checksum for the file in ad storage
    """

    if ((fqn.endswith('.gz') and local_meta['md5sum'] !=
         ad_meta['md5sum']) or (
            not fqn.endswith('.gz') and local_meta['md5sum'] !=
            ad_meta['umd5sum'])):
        raise CadcException(
            '{} md5sum not the same as the one in the ad '
            '{} archive.'.format(fqn, archive))


def compare_checksum(netrc_fqn, archive, fqn):
    """
    Raise CadcException if the checksum of a file in ad is not the same as
    the checksum of a file on disk.

    :param netrc_fqn: fully-qualified file name for the netrc file
    :param archive: archive file has been stored to
    :param fqn: Fully-qualified name of the file for which to get the metadata.
    """
    fname = os.path.basename(fqn)
    try:
        local_meta = get_file_meta(fqn)
        ad_meta = get_cadc_meta(netrc_fqn, archive, fname)
    except Exception as e:
        raise CadcException('Could not find md5 checksum for {} in the ad {} '
                            'archive. {}'.format(fqn, archive, e))
    _check_checksums(fqn, archive, local_meta, ad_meta)


def compare_checksum_client(client, archive, fqn):
    """
    Raise CadcException if the checksum of a file in ad is not the same as
    the checksum of a file on disk.

    :param client: access to CADC data service
    :param archive: archive file has been stored to
    :param fqn: Fully-qualified name of the file for which to get the metadata.
    """
    fname = os.path.basename(fqn)
    try:
        local_meta = get_file_meta(fqn)
        ad_meta = client.get_file_info(archive, fname)
    except Exception as e:
        raise CadcException('Could not find md5 checksum for {} in the ad {} '
                            'archive. {}'.format(fqn, archive, e))
    _check_checksums(fqn, archive, local_meta, ad_meta)


def create_dir(dir_name):
    """Create the working area if it does not already exist."""
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
        if not os.path.exists(dir_name):
            raise CadcException(
                'Could not mkdir {}'.format(dir_name))
        if not os.access(dir_name, os.W_OK | os.X_OK):
            raise CadcException(
                '{} is not writeable.'.format(dir_name))


def decompose_lineage(lineage):
    """Returns a product id and an artifact uri from the command line."""
    try:
        result = lineage.split('/', 1)
        return result[0], result[1]
    except Exception as e:
        logging.debug('Lineage {} caused error {}. Expected '
                      'product_id/ad:ARCHIVE/FILE_NAME'.format(
                        lineage, e))
        raise CadcException('Expected product_id/ad:ARCHIVE/FILE_NAME')


def check_param(param, param_type):
    """Generic code to check if a parameter is not None, and is of the
    expected type."""
    if param is None or not isinstance(param, param_type):
        raise CadcException(
            'Parameter {} failed check for {}'.format(param, param_type))


def read_csv_file(fqn):
    """Read a csv file.

    :returns a list of lists."""
    results = []
    try:
        with open(fqn) as csv_file:
            reader = csv.reader(csv_file)
            for row in reader:
                if row[0].startswith('#'):
                    continue
                results.append(row)
    except Exception as e:
        logging.error('Could not read from csv file {}'.format(fqn))
        raise CadcException(e)
    return results


def write_obs_to_file(obs, fqn):
    """Common code to write a CAOM Observation to a file."""
    ow = ObservationWriter()
    ow.write(obs, fqn)


def read_obs_from_file(fqn):
    """Common code to read a CAOM Observation from a file."""
    if not os.path.exists(fqn):
        raise CadcException('Could not find {}'.format(fqn))
    reader = ObservationReader(False)
    return reader.read(fqn)


def read_from_file(fqn):
    """Common code to read from a text file. Mostly to make it easy to
    mock."""
    if not os.path.exists(fqn):
        raise CadcException('Could not find {}'.format(fqn))
    with open(fqn, 'r') as f:
        return f.readlines()


def write_to_file(fqn, content):
    """Common code to write to a fully-qualified file name. Mostly to make
    it easy to mock."""
    try:
        with open(fqn, 'w') as f:
            f.write(content)
    except Exception:
        logging.error('Could not write file {}'.format(fqn))
        raise CadcException('Could not write file {}'.format(fqn))


def update_typed_set(typed_set, new_set):
    """Common code to remove all the entries from an existing set, and
    then replace those entries with a new set."""
    # remove the previous values
    while len(typed_set) > 0:
        typed_set.pop()
    typed_set.update(new_set)


def format_time_for_query(from_time):
    length = len(datetime.now().strftime('%Y-%m-%dT%H:%M:%S'))
    return datetime.strptime(from_time[:length], '%Y-%m-%dT%H:%M:%S')


def read_file_list_from_archive(config, app_name, prev_exec_date, exec_date):
    """Code to execute a time-boxed query for files that have arrived in ad.

    :param config Config instance
    :param app_name Information used in the http connection for tracing
        queries.
    :param prev_exec_date Timestamp that indicates the beginning of the
        chunk of time. Results will be > than this time.
    :param exec_date Timestamp. that indicates the end of the chunk of time.
        Results will be <= this time.
    """
    start_time = format_time_for_query(prev_exec_date)
    end_time = format_time_for_query(exec_date)
    ad_resource_id = 'ivo://cadc.nrc.ca/ad'
    agent = '{}/{}'.format(app_name, '1.0')
    subject = net.Subject(certificate=config.proxy_fqn)
    client = net.BaseWsClient(resource_id=ad_resource_id,
                              subject=subject, agent=agent, retry=True)
    query_meta = "SELECT fileName FROM archive_files WHERE " \
                 "archiveName = '{}' AND ingestDate > '{}' and " \
                 "ingestDate <= '{}' ORDER BY ingestDate".format(
                    config.archive, start_time, end_time)
    data = {'QUERY': query_meta, 'LANG': 'ADQL', 'FORMAT': 'csv'}
    logging.debug('Query is {}'.format(query_meta))
    try:
        response = client.get('https://{}/ad/sync?{}'.format(
            client.host, parse.urlencode(data)), cert=config.proxy_fqn)
        if response.status_code == 200:
            # ignore the column name as the first part of the response
            artifact_files_list = response.text.split()[1:]
            return artifact_files_list
        else:
            logging.warning('No work to do. Query failure {!r}'.format(
                response))
            return []
    except Exception as e:
        raise CadcException('Failed ad content query: {}'.format(e))


def get_lineage(archive, product_id, file_name, scheme='ad'):
    """Construct an instance of the caom2gen lineage parameter.
    :param archive archive name at CADC.
    :param product_id CAOM2 Plane unique identifier.
    :param file_name String representation of the file name.
    :param scheme Usually 'ad', otherwise an indication of external storage.
    :return str understood by the caom2gen application, lineage parameter
        value"""
    return '{}/{}:{}/{}'.format(product_id, scheme, archive, file_name)


def get_artifact_metadata(fqn, product_type, release_type, uri=None,
                          artifact=None):
    """
    Build or update artifact content metadata using the CAOM2 objects, and
    with access to a file on disk.

    :param fqn: The fully-qualified name of the file on disk, for which an
        Artifact is being created or updated.
    :param product_type: which ProductType enumeration value
    :param release_type: which ReleaseType enumeration value
    :param uri: mandatory if creating an Artifact, a URI of the form
        scheme:ARCHIVE/file_name
    :param artifact: use when updating an existing Artifact instance

    :return: the created or updated Artifact instance, with the
        content_* elements filled in.
    """
    local_meta = get_file_meta(fqn)
    md5uri = ChecksumURI('md5:{}'.format(local_meta['md5sum']))
    if artifact is None:
        if uri is None:
            raise CadcException('Cannot build an Artifact without a URI.')
        return Artifact(uri, product_type, release_type, local_meta['type'],
                        local_meta['size'], md5uri)
    else:
        artifact.product_type = product_type
        artifact.content_type = local_meta['type']
        artifact.content_length = local_meta['size']
        artifact.content_checksum = md5uri
        return artifact


def data_put(client, working_directory, file_name, archive, stream='raw',
             mime_type=None):
    """
    Make a copy of a locally available file by writing it to CADC. Assumes
    file and directory locations are correct. Does a checksum comparison to
    test whether the file made it to storage as it exists on disk.

    :param client: The CadcDataClient for write access to CADC storage.
    :param working_directory: Where 'file_name' exists locally.
    :param file_name: What to copy to CADC storage.
    :param archive: Which archive to associate the file with.
    :param stream: Defaults to raw - use is deprecated, however necessary it
        may be at the current moment to the 'put_file' call.
    :param mime_type: Because libmagic can't see inside a zipped fits file.
    """
    cwd = os.getcwd()
    try:
        os.chdir(working_directory)
        client.put_file(archive, file_name, archive_stream=stream,
                        mime_type=mime_type)
    except Exception as e:
        raise CadcException('Failed to store data with {}'.format(e))
    finally:
        os.chdir(cwd)
    compare_checksum_client(client, archive,
                            os.path.join(working_directory, file_name))


def data_get(client, working_directory, file_name, archive):
    """
    Retrieve a local copy of a file available from CADC. Assumes the working
    directory location exists and is writeable.

    :param client: The CadcDataClient for read access to CADC storage.
    :param working_directory: Where 'file_name' will be written.
    :param file_name: What to copy from CADC storage.
    :param archive: Which archive to retrieve the file from.
    """
    fqn = os.path.join(working_directory, file_name)
    try:
        client.get_file(archive, file_name, destination=fqn)
        if not os.path.exists(fqn):
            raise CadcException(
                'Retrieve failed. {} does not exist.'.format(fqn))
    except Exception as e:
        raise CadcException('Did not retrieve {} because {}'.format(
            fqn, e))


def build_uri(archive, file_name, scheme='ad'):
    """One location to keep the syntax for an Artifact URI."""
    return '{}:{}/{}'.format(scheme, archive, file_name)


def query_endpoint(url, timeout=20):
    """Return a response for an endpoint. Caller needs to call 'close'
    on the response.
    """

    # Open the URL and fetch the JSON document for the observation
    session = requests.Session()
    retries = 10
    retry = Retry(total=retries, read=retries, connect=retries,
                  backoff_factor=0.5)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    try:
        response = session.get(url, timeout=timeout)
        response.raise_for_status()
        return response
    except Exception as e:
        raise CadcException('Endpoint {} failure {}'.format(url, str(e)))


def read_as_yaml(fqn):
    """Read and return YAML content of 'fqn'."""
    try:
        logging.debug('Begin read_as_yaml for {}.'.format(fqn))
        with open(fqn) as f:
            data_map = yaml.safe_load(f)
            logging.debug('End read_as_yaml.')
            return data_map
    except (yaml.scanner.ScannerError, FileNotFoundError) as e:
        logging.error(e)
        return None


def write_as_yaml(content, fqn):
    """Write 'content' to 'fqn' as YAML."""
    try:
        logging.debug('Begin write_as_yaml for {}.'.format(fqn))
        with open(fqn, 'w') as f:
            yaml.dump(content, f, default_flow_style=False)
            logging.debug('End write_as_yaml.')
    except Exception as e:
        logging.error(e)
