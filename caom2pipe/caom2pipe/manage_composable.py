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
import io
import logging
import os
import requests
import subprocess
import sys
import yaml

from datetime import datetime
from enum import Enum
from ftplib import FTP
from ftputil import FTPHost
from hashlib import md5
from io import BytesIO
from pytz import timezone
from requests.adapters import HTTPAdapter
from urllib import parse as parse
from urllib3 import Retry

from astropy.io.votable import parse_single_table
from astropy.table import Table

from cadcutils import net, exceptions
from cadcdata import CadcDataClient
from cadctap import CadcTapClient
from caom2 import ObservationWriter, ObservationReader, Artifact
from caom2 import ChecksumURI


__all__ = ['CadcException', 'Config', 'State', 'to_float', 'TaskType',
           'exec_cmd', 'exec_cmd_redirect', 'exec_cmd_info',
           'get_cadc_meta', 'get_file_meta',
           'decompose_lineage', 'check_param', 'read_csv_file',
           'write_obs_to_file', 'read_obs_from_file',
           'Features', 'write_to_file',
           'read_from_file', 'read_file_list_from_archive', 'update_typed_set',
           'get_cadc_headers', 'get_lineage', 'get_artifact_metadata',
           'data_put', 'data_get', 'build_uri', 'make_seconds',
           'increment_time', 'ISO_8601_FORMAT', 'http_get', 'Rejected',
           'record_progress', 'Builder', 'Work', 'look_pull_and_put',
           'Observable', 'Metrics', 'repo_create', 'repo_delete', 'repo_get',
           'repo_update', 'ftp_get', 'ftp_get_timeout', 'VALIDATE_OUTPUT']

ISO_8601_FORMAT = '%Y-%m-%dT%H:%M:%S.%f'
READ_BLOCK_SIZE = 8 * 1024
VALIDATE_OUTPUT = 'validated.yml'


class CadcException(Exception):
    """Generic exception raised by failure cases within the caom2pipe
    module."""
    pass


class Features(object):
    """Boolean feature flag implementation."""

    def __init__(self):
        self._use_file_names = True
        self._use_urls = True
        self._run_in_airflow = True
        self._supports_composite = True
        self._supports_catalog = True
        self._expects_retry = True

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
    def use_urls(self):
        """If true, the lists of work to be done are expected to
        identify URLs. If false, they are expected to identify
        observation IDs. This and use_file_names cannot both be True."""
        return self._use_urls

    @use_urls.setter
    def use_urls(self, value):
        self._use_urls = value

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
                logging.debug('Saving timestamp {} {}'.format(value, self.fqn))
                write_as_yaml(self.content, self.fqn)
            else:
                self.logger.warning('No record found for {}'.format(key))
        else:
            self.logger.warning('No bookmarks found for {}'.format(key))


class Builder(object):
    """"Abstract-like class that defines the operations used to build
    StorageName instances"""

    def __init__(self, config):
        """
        Ctor
        :param config: the source of all knowledge.
        """
        self._config = config
        self._todo_list = None

    @property
    def todo_list(self):
        return self._todo_list

    @todo_list.setter
    def todo_list(self, to_list):
        """Set the list of work to be done. It probably contains useful
        information like timestamps, and relationships between obs ids
        and file names, hence make it possible to set it independently."""
        raise NotImplementedError

    def build(self, entry):
        """Returns an instance of StorageName, based on all the configuration
         information, and the knowledge about the relationship between
         file names and CAOM2 Observations that is specific to a Collection."""
        raise NotImplementedError


class Work(object):
    """"Abstract-like class that defines the operations used to chunk work when
    controlling execution by State."""

    def __init__(self, max_ts_s):
        """
        Ctor
        :param max_ts_s: the maximum timestamp in seconds, for work that
        is chunked by time-boxes.
        """
        self._max_ts_s = max_ts_s

    @property
    def max_ts_s(self):
        # State execution is currently chunked only by time-boxing, so
        # set the maximum time-box ending for the work to be done.
        return self._max_ts_s

    def initialize(self):
        """Anything necessary to make todo work."""
        raise NotImplementedError

    def todo(self, prev_exec_date, exec_date):
        """Returns a list of entries for processing by Execute.
        :param prev_exec_date when chunking by time-boxes, the start time
        :param exec_date when chunking by time-boxes, the end time"""
        raise NotImplementedError


class Rejected(object):
    """Persist information between pipeline invocations about the observation
    IDs that will fail a particular TaskType.
    """
    NO_REASON = ''
    BAD_METADATA = 'bad_metadata'
    NO_PREVIEW = 'no_preview'

    # A map to the logging message string representing acknowledged rejections
    reasons = {BAD_METADATA: 'Cannot build an observation',
               NO_PREVIEW: '404 Client Error: Not Found for url'}

    def __init__(self, fqn):
        """
        Ctor
        :param fqn: fully-qualified name for reading/writing rejected
        information records. If the file exists, initialize the in-memory
        records with the content of the file. Otherwise initialize the
        in-memory records to be empty.
        """
        self.fqn = fqn
        if os.path.exists(fqn):
            try:
                self.content = read_as_yaml(fqn)
            except yaml.constructor.ConstructorError:
                logging.error('ConstructorError reading {}'.format(self.fqn))
                self.content = {}
        else:
            dir_name = os.path.dirname(fqn)
            create_dir(dir_name)
            self.content = {}
        for reason in Rejected.reasons.keys():
            if reason not in self.content:
                self.content[reason] = []

    def check_and_record(self, message, obs_id):
        """Keep track of an additional entry.
        :returns boolean True if the failure is known and tracked."""
        for reason, reason_str in Rejected.reasons.items():
            if reason_str in message:
                self.record(reason, obs_id)
                return True
        return False

    def record(self, reason, entry):
        """Keep track of an additional entry."""
        self.content[reason].append(entry)

    def persist_state(self):
        """Write the current state as a YAML file."""
        for key, value in self.content.items():
            # ensure unique entries
            self.content[key] = list(set(value))
        write_as_yaml(self.content, self.fqn)

    def is_bad_metadata(self, obs_id):
        return obs_id in self.content[Rejected.BAD_METADATA]

    def is_no_preview(self, file_name):
        return file_name in self.content[Rejected.NO_PREVIEW]

    @staticmethod
    def known_failure(message):
        """Returns the REASON for the failure, or an empty string if
        the failure is of an unexpected type.
        """
        result = Rejected.NO_REASON
        for reason, text in Rejected.reasons.items():
            if text in message:
                result = reason
        return result


class Metrics(object):
    """
    A class to capture execution metrics.
    """

    def __init__(self, config):
        self.enabled = config.observe_execution
        if self.enabled:
            self.history = {}
            self.failures = {}
            self.observable_dir = config.observable_directory

    def observe(self, start, stop, size, action, service, label):
        if self.enabled:
            elapsed = round(stop - start, 3)
            rate = round(size / (stop - start), 3)
            if service not in self.history:
                self.history[service] = {}
            if action not in self.history[service]:
                self.history[service][action] = {}
            self.history[service][action][label] = [elapsed, rate, start]

    def observe_failure(self, action, service, label):
        if self.enabled:
            if service not in self.failures:
                self.failures[service] = {}
            if action not in self.failures[service]:
                self.failures[service][action] = {}
            if label not in self.failures[service][action]:
                self.failures[service][action][label] = 1
            else:
                self.failures[service][action][label] += 1

    def capture(self):
        if self.enabled:
            create_dir(self.observable_dir)
            now = datetime.utcnow().timestamp()
            for service in self.history.keys():
                fqn = '{}/{}.{}.yml'.format(self.observable_dir, now, service)
                write_as_yaml(self.history[service], fqn)

            fqn = '{}/{}.fail.yml'.format(self.observable_dir, now)
            # if os.path.exists(fqn):
            #     temp = read_as_yaml(fqn)
            write_as_yaml(self.failures, fqn)


class Observable(object):
    """
    A class to contain all the classes that maintain information between
    pipeline execution instances.
    """

    def __init__(self, rejected, metrics):
        self._rejected = rejected
        self._metrics = metrics

    @property
    def rejected(self):
        return self._rejected

    @property
    def metrics(self):
        return self._metrics


class Config(object):
    """Configuration information that remains the same for all steps and all
     work in a pipeline execution."""

    def __init__(self):
        self._working_directory = None
        self._work_file = None
        # the fully qualified name for the work file
        self.work_fqn = None
        self._netrc_file = None
        self._archive = None
        self._collection = None
        self._use_local_files = False
        self._resource_id = None
        self._tap_id = None
        self._logging_level = None
        self._log_to_file = False
        self._log_file_directory = None
        self._stream = None
        self._storage_host = None
        self._task_types = None
        self._success_log_file_name = None
        # the fully qualified name for the file
        self.success_fqn = None
        self._failure_log_file_name = None
        # the fully qualified name for the file
        self.failure_fqn = None
        self._retry_file_name = None
        # the fully qualified name for the file
        self.retry_fqn = None
        self._retry_failures = False
        self._retry_count = 1
        self._proxy_file_name = None
        # the fully qualified name for the file
        self.proxy_fqn = None
        self._state_file_name = None
        # the fully qualified name for the file
        self.state_fqn = None
        self._rejected_file_name = None
        self._rejected_directory = None
        # the fully qualified name for the file
        self.rejected_fqn = None
        self._progress_file_name = None
        self.progress_fqn = None
        self._interval = None
        self._observe_execution = False
        self._observable_directory = None
        self._source_host = None
        self._features = Features()

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
        if self._working_directory is not None:
            self.work_fqn = os.path.join(
                self._working_directory, self._work_file)

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
    def task_types(self):
        """the way to control which steps get executed"""
        return self._task_types

    @task_types.setter
    def task_types(self, value):
        self._task_types = value

    @property
    def success_log_file_name(self):
        """the filename where success logs are written, this will be created
        in log_file_directory"""
        return self._success_log_file_name

    @success_log_file_name.setter
    def success_log_file_name(self, value):
        self._success_log_file_name = value
        if self._log_file_directory is not None:
            self.success_fqn = os.path.join(
                self._log_file_directory, self._success_log_file_name)

    @property
    def failure_log_file_name(self):
        """the filename where failure logs are written this will be created
        in log_file_directory"""
        return self._failure_log_file_name

    @failure_log_file_name.setter
    def failure_log_file_name(self, value):
        self._failure_log_file_name = value
        if self._log_file_directory is not None:
            self.failure_fqn = os.path.join(
                self._log_file_directory, self._failure_log_file_name)

    @property
    def retry_file_name(self):
        """the filename where retry entries are written this will be created
        in log_file_directory"""
        return self._retry_file_name

    @retry_file_name.setter
    def retry_file_name(self, value):
        self._retry_file_name = value
        if self._log_file_directory is not None:
            self.retry_fqn = os.path.join(
                self._log_file_directory, self._retry_file_name)

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
    def rejected_directory(self):
        """the directory where rejected entry files are written, by default
        will be the log file directory"""
        return self._rejected_directory

    @rejected_directory.setter
    def rejected_directory(self, value):
        self._rejected_directory = value
        if self._rejected_directory is not None:
            self.rejected_fqn = os.path.join(
                self._rejected_directory, self._rejected_file_name)
        elif self._log_file_directory is not None:
            self.rejected_fqn = os.path.join(
                self._log_file_directory, self._rejected_file_name)

    @property
    def rejected_file_name(self):
        """the filename where rejected entries are written, this will be created
        in log_file_directory"""
        return self._rejected_file_name

    @rejected_file_name.setter
    def rejected_file_name(self, value):
        self._rejected_file_name = value
        if self._log_file_directory is not None:
            self.rejected_fqn = os.path.join(
                self._log_file_directory, self._rejected_file_name)

    @property
    def progress_file_name(self):
        """the filename where pipeline progress is written, this will be created
        in log_file_directory. Useful when using timestamp windows for
        execution. """
        return self._progress_file_name

    @progress_file_name.setter
    def progress_file_name(self, value):
        self._progress_file_name = value
        if self._log_file_directory is not None:
            self.progress_fqn = os.path.join(
                self._log_file_directory, self._progress_file_name)

    @property
    def proxy_file_name(self):
        """If using a proxy certificate for authentication, identify the
        fully-qualified pathname here."""
        return self._proxy_file_name

    @proxy_file_name.setter
    def proxy_file_name(self, value):
        self._proxy_file_name = value
        if (self._working_directory is not None and
                self._proxy_file_name is not None):
            self.proxy_fqn = os.path.join(
                self._working_directory, self._proxy_file_name)

    @property
    def state_file_name(self):
        """If using a state file to communicate persistent information between
        invocations, identify the fully-qualified pathname here."""
        return self._state_file_name

    @state_file_name.setter
    def state_file_name(self, value):
        self._state_file_name = value
        if (self._working_directory is not None and
                self._state_file_name is not None):
            self.state_fqn = os.path.join(
                self._working_directory, self._state_file_name)

    @property
    def features(self):
        """Feature flag setting access."""
        return self._features

    @features.setter
    def features(self, value):
        self._features = value

    @property
    def interval(self):
        """Interval - used for setting timestamp chunks, for now."""
        return self._interval

    @interval.setter
    def interval(self, value):
        self._interval = value

    @property
    def observe_execution(self):
        """If true, time and track the CADC service invocations."""
        return self._observe_execution

    @observe_execution.setter
    def observe_execution(self, value):
        self._observe_execution = value

    @property
    def observable_directory(self):
        """Directory where to track CADC service invocation information."""
        return self._observable_directory

    @observable_directory.setter
    def observable_directory(self, value):
        self._observable_directory = value

    @property
    def source_host(self):
        """Host that is the source of something. Initial use case as ftp
        host name."""
        return self._source_host

    @source_host.setter
    def source_host(self, value):
        self._source_host = value

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
               'rejected_directory:: \'{}\' ' \
               'rejected_file_name:: \'{}\' ' \
               'rejected_fqn:: \'{}\' ' \
               'progress_file_name:: \'{}\' ' \
               'progress_fqn:: \'{}\' ' \
               'proxy_file:: \'{}\' ' \
               'state_fqn:: \'{}\' ' \
               'features:: \'{}\' ' \
               'interval:: \'{}\' ' \
               'observe_execution:: \'{}\' ' \
               'observable_directory:: \'{}\' ' \
               'source_host:: \'{}\' ' \
               'logging_level:: \'{}\''.format(
                self.working_directory, self.work_fqn, self.netrc_file,
                self.archive, self.collection, self.task_types, self.stream,
                self.resource_id,
                self.tap_id, self.use_local_files, self.log_to_file,
                self.log_file_directory, self.success_log_file_name,
                self.success_fqn, self.failure_log_file_name,
                self.failure_fqn, self.retry_file_name, self.retry_fqn,
                self.retry_failures, self.retry_count, self.rejected_file_name,
                self.rejected_directory, self.rejected_fqn,
                self.progress_file_name,
                self.progress_fqn, self.proxy_fqn, self.state_fqn,
                self.features, self.interval, self.observe_execution,
                self.observable_directory, self.source_host,
                self.logging_level)

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
            self.working_directory = config.get('working_directory',
                                                os.getcwd())
            self.work_file = config.get('todo_file_name', 'todo.txt')
            self.netrc_file = config.get('netrc_filename', 'test_netrc')
            self.resource_id = config.get('resource_id',
                                          'ivo://cadc.nrc.ca/sc2repo')
            self.tap_id = config.get('tap_id', 'ivo://cadc.nrc.ca/sc2tap')
            self.use_local_files = bool(config.get('use_local_files', False))
            self.logging_level = config.get('logging_level', 'DEBUG')
            self.log_to_file = config.get('log_to_file', False)
            self.log_file_directory = config.get('log_file_directory',
                                                 self.working_directory)
            self.stream = config.get('stream', 'raw')
            self.task_types = self._obtain_task_types(
                config, [TaskType.SCRAPE])
            self.collection = config.get('collection', 'TEST')
            self.archive = config.get('archive', self.collection)
            self.success_log_file_name = config.get('success_log_file_name',
                                                    'success_log.txt')
            self.failure_log_file_name = config.get('failure_log_file_name',
                                                    'failure_log.txt')
            self.retry_file_name = config.get('retry_file_name',
                                              'retries.txt')
            self.retry_failures = config.get('retry_failures', False)
            self.retry_count = config.get('retry_count', 1)
            self.rejected_file_name = config.get('rejected_file_name',
                                                 'rejected.yml')
            self.rejected_directory = config.get('rejected_directory',
                                                 os.getcwd())
            self.progress_file_name = config.get('progress_file_name',
                                                 'progress.txt')
            self.interval = config.get('interval', 10)
            self.features = self._obtain_features(config)
            self.proxy_file_name = config.get('proxy_file_name', None)
            self.state_file_name = config.get('state_file_name', None)
            self.observe_execution = config.get('observe_execution', False)
            self.observable_directory = config.get(
                'observable_directory', None)
            self.source_host = config.get('source_host', None)
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
            logging.debug('Begin load_config from {}.'.format(config_fqn))
            with open(config_fqn) as f:
                data_map = yaml.safe_load(f)
                logging.debug('End load_config.')
                return data_map
        except (yaml.scanner.ScannerError, FileNotFoundError) as e:
            logging.error(e)
            return None


class Validator(object):
    """
    Compares the state of CAOM entries at CADC with the state of the source
    of the data. Identifies files that are in CAOM entries that do not exist
    at the source, and files at the source that are not represented in CAOM.
    Checks that the timestamp associated with the file at the source is less
    than the ad timestamp.

    CADC is the destination, where the data and metadata originate from
    is the source.
    """
    def __init__(self, source_name, scheme='ad', preview_suffix='jpg',
                 source_tz='UTC'):
        """

        :param source_name: String value used for logging
        :param scheme: string which encapsulates scheme as used in CAOM
            Artifact URIs. The default of 'ad' means the canonical version
            of the file is stored at CADC.
        :param preview_suffix String value that is excluded from queries,
            because it's produced at CADC, so something like '256.jpg' should
            work for Gemini.
        :param source_tz String representation of timezone name, as understood
            by pytz.
        """
        self._config = Config()
        self._config.get_executors()
        logger = logging.getLogger()
        logger.setLevel(self._config.logging_level)
        logging.debug(self._config)
        self._source = None
        self._destination_data = None
        self._destination_meta = None
        self._source_name = source_name
        self._scheme = scheme
        self._preview_suffix = preview_suffix
        self._source_tz = timezone(source_tz)

    def _find_unaligned_dates(self, source, meta, data):
        result = []
        for f_name in meta:
            if f_name in source and f_name in data:
                source_dt = datetime.strptime(source[f_name], ISO_8601_FORMAT)
                source_utc = source_dt.astimezone(timezone(self._source_tz))
                dest_dt = datetime.strptime(data[f_name], ISO_8601_FORMAT)
                dest_utc = dest_dt.astimezone(timezone('US/Pacific'))
                if dest_utc < source_utc:
                    result.append(f_name)
        return list(set(result))

    def _read_list_from_destination_data(self):
        """Code to execute a query for files and the arrival time, that are in
        CADC storage.
        """
        ad_resource_id = 'ivo://cadc.nrc.ca/ad'
        agent = '{}/{}'.format(self._source_name, '1.0')
        subject = define_subject(self._config)
        client = net.BaseWsClient(resource_id=ad_resource_id,
                                  subject=subject, agent=agent, retry=True)
        query_meta = f"SELECT fileName, ingestDate FROM archive_files WHERE " \
                     f"archiveName = '{self._config.archive}' " \
                     f"AND fileName not like '%{self._preview_suffix} " \
                     f"ORDER BY fileName"
        data = {'QUERY': query_meta, 'LANG': 'ADQL', 'FORMAT': 'csv'}
        logging.debug('Query is {}'.format(query_meta))
        try:
            response = client.get(
                f'https://{client.host}/ad/sync?{parse.urlencode(data)}')
            if response.status_code == 200:
                table = Table.read(response.text, format='csv')
                return table
            else:
                logging.warning('No work to do. Query failure {!r}'.format(
                    response))
                return Table()
        except Exception as e:
            raise CadcException('Failed ad content query: {}'.format(e))

    def _read_list_from_destination_meta(self):
        query = f"SELECT A.uri FROM caom2.Observation AS O " \
                f"JOIN caom2.Plane AS P ON O.obsID = P.obsID " \
                f"JOIN caom2.Artifact AS A ON P.planeID = A.planeID " \
                f"WHERE O.collection='{self._config.collection}' " \
                f"AND A.uri not like '%{self._preview_suffix}'"
        subject = net.Subject(certificate=self._config.proxy_fqn)
        tap_client = CadcTapClient(subject, resource_id=self._config.tap_id)
        buffer = io.BytesIO()
        tap_client.query(query, output_file=buffer, data_only=True)
        temp = parse_single_table(buffer).to_table()
        return [ii.decode().replace(
            f'{self._scheme}:{self._config.archive}/', '').strip()
                for ii in temp['uri']]

    def read_from_source(self):
        """Read the entire source site listing. This function is expected to
        return a dict of all the file names available from the source, where
        the keys are the file names, and the values are the timestamps at the
        source."""
        raise NotImplementedError()

    def validate(self):
        logging.info('Query destination metadata.')
        dest_meta_temp = self._read_list_from_destination_meta()

        logging.info('Query source metadata.')
        source_temp = self.read_from_source()

        logging.info('Find files that are missing from CADC.')
        self._destination_meta = find_missing(
            dest_meta_temp, source_temp.keys())

        logging.info(f'Find files that do not appear at {self._source_name}.')
        self._source = find_missing(source_temp.keys(), dest_meta_temp)

        logging.info('Query destination data.')
        dest_data_temp = self._read_list_from_destination_data()

        logging.info(f'Find files that are newer at {self._source_name} '
                     f'than at CADC.')
        self._destination_data = self._find_unaligned_dates(
            source_temp, dest_meta_temp, dest_data_temp)

        logging.info('Log the results.')
        result = {f'{self._source_name}': self._source,
                  'cadc': self._destination_meta,
                  'timestamps': self._destination_data}
        result_fqn = os.path.join(
            self._config.working_directory, VALIDATE_OUTPUT)
        write_as_yaml(result, result_fqn)

        logging.info(f'There are {len(self._source)} files at '
                     f'{self._source_name} that are not represented at CADC, '
                     f'and {len(self._destination_meta)} CAOM entries at '
                     f'CADC that are not available from {self._source_name}.\n'
                     f'There are {len(self._destination_data)} files that are '
                     f'newer at {self._source_name}.')
        return self._source, self._destination_meta, self._destination_data

    def write_todo(self):
        """Write a todo.txt file, given the list of entries available from
        the source, that are not currently at the destination (CADC)."""
        raise NotImplementedError()


def to_float(value):
    """Cast to float, without throwing an exception."""
    return float(value) if value is not None else None


def to_int(value):
    """Cast to int, without throwing an exception."""
    return int(value) if value is not None else None


def define_subject(config):
    """Common code to figure out which credentials to use based on the
    content of a Config instance."""
    subject = None
    if config.proxy_fqn is not None and os.path.exists(config.proxy_fqn):
        logging.debug(
            f'Using proxy certificate {config.proxy_fqn} for credentials.')
        subject = net.Subject(username=None,
                              certificate=config.proxy_fqn)
    elif config.netrc_file is not None:
        netrc_fqn = os.path.join(config.working_directory, config.netrc_file)
        if os.path.exists(netrc_fqn):
            logging.error(
                f'Using netrc file {netrc_fqn} for credentials.')
            subject = net.Subject(username=None, certificate=None,
                                  netrc=netrc_fqn)
        else:
            logging.warning(f'Cannot find netrc file {netrc_fqn}')
    else:
        logging.warning(
            'No credentials provided (proxy certificate or netrc file).')
        logging.warning(
            'Proxy certificate is {}, netrc file is {}.'.format(
                config.proxy_fqn, config.netrc_file))
    return subject


def exec_cmd(cmd, log_leval_as=logging.debug):
    """
    This does command execution as a subprocess call.

    :param cmd the text version of the command being executed
    :param log_leval_as control the logging level from the exec call
    :return None
    """
    logging.debug(cmd)
    cmd_array = cmd.split()
    try:
        child = subprocess.Popen(cmd_array, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
        output, outerr = child.communicate()
        if len(output) > 0:
            log_leval_as('stdout {}'.format(output.decode('utf-8')))
        if len(outerr) > 0:
            logging.error('stderr {}'.format(outerr.decode('utf-8')))
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


def ftp_get(ftp_host_name, source_fqn, dest_fqn):
    """
    :param ftp_host_name name from which to originate the FTP transfer. Assume
        anonymous login.
    :param source_fqn fully-qualified name on the FTP host, of the file to be
        transferred.
    :param dest_fqn fully-qualified name, locally valid, for where to transfer
        the file to.

    Uses ftputil, which always transfers files in binary mode.
    """
    try:
        with FTPHost(ftp_host_name, 'anonymous', '@anonymous') as ftp_host:
            ftp_host.download(source_fqn, dest_fqn)
            source_stats = ftp_host.stat(source_fqn)
            ftp_host.close()
            dest_meta = get_file_meta(dest_fqn)
            if source_stats.st_size == dest_meta.get('size'):
                logging.info('Downloaded {} from {}'.format(
                    source_fqn, ftp_host_name))
            else:
                raise CadcException(
                    'File size error when transferring {} from {}'.format(
                        source_fqn, ftp_host_name))
    except Exception as e:
        logging.error(e)
        raise CadcException('Could not transfer {} from {}'.format(
            source_fqn, ftp_host_name))


def ftp_get_timeout(ftp_host_name, source_fqn, dest_fqn, timeout=20):
    """
    :param ftp_host_name name from which to originate the FTP transfer. Assume
        anonymous login.
    :param source_fqn fully-qualified name on the FTP host, of the file to be
        transferred.
    :param dest_fqn fully-qualified name, locally valid, for where to transfer
        the file to.
    :param timeout in seconds for blocking operations

    Uses ftplib, which supports specifying timeouts in the connection.
    """
    try:
        with FTP(ftp_host_name, timeout=timeout) as ftp_host:
            ftp_host.login()
            with open(dest_fqn, 'wb') as fp:
                ftp_host.retrbinary(f'RETR {source_fqn}', fp.write)
            ftp_host.voidcmd('TYPE I')
            source_size = ftp_host.size(source_fqn)
            ftp_host.quit()
            dest_meta = get_file_meta(dest_fqn)
            if source_size == dest_meta.get('size'):
                logging.info('Downloaded {} from {}'.format(
                    source_fqn, ftp_host_name))
            else:
                raise CadcException(
                    'File size error when transferring {} from {}'.format(
                        source_fqn, ftp_host_name))
    except Exception as e:
        logging.error(e)
        raise CadcException('Could not transfer {} from {}'.format(
            source_fqn, ftp_host_name))


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
    s = os.stat(fqn)
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
        meta['type'] = 'application/tar+gzip'
    elif fqn.endswith('h5'):
        meta['type'] = 'application/x-hdf5'
    else:
        meta['type'] = 'application/fits'
    logging.debug(meta)
    return meta


def create_dir(dir_name):
    """Create the working area if it does not already exist."""
    if os.path.exists(dir_name):
        if os.path.isfile(dir_name):
            raise CadcException(f'{dir_name} already exists as a file.')
    else:
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
    expected type.
    """
    if param is None or not isinstance(param, param_type):
        raise CadcException(
            'Parameter {} failed check for {}'.format(param, param_type))


def read_csv_file(fqn):
    """Read a csv file.

    :returns a list of lists.
    """
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


def record_progress(config, application, count, cumulative, start_time):
    """Common code to write the number of entries processed to a file."""
    with open(config.progress_fqn, 'a') as progress:
        progress.write(
            '{} {} current:: {} since {}:: {}\n'.format(
                datetime.now(), application, count, start_time, cumulative))


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
    mock.
    """
    if not os.path.exists(fqn):
        raise CadcException('Could not find {}'.format(fqn))
    with open(fqn, 'r') as f:
        return f.readlines()


def write_to_file(fqn, content):
    """Common code to write to a fully-qualified file name. Mostly to make
    it easy to mock.
    """
    try:
        with open(fqn, 'w') as f:
            f.write(content)
    except Exception:
        logging.error('Could not write file {}'.format(fqn))
        raise CadcException('Could not write file {}'.format(fqn))


def update_typed_set(typed_set, new_set):
    """Common code to remove all the entries from an existing set, and
    then replace those entries with a new set.
    """
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
        value
    """
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


def current():
    """Encapsulate returning UTC now in microsecond resolution."""
    return datetime.utcnow().timestamp()


def sizeof(x):
    """Encapsulate returning the memory size in bytes."""
    return sys.getsizeof(x)


def data_put(client, working_directory, file_name, archive, stream='raw',
             mime_type=None, mime_encoding=None, metrics=None):
    """
    Make a copy of a locally available file by writing it to CADC. Assumes
    file and directory locations are correct. Requires a checksum comparison
    by the client.

    :param client: The CadcDataClient for write access to CADC storage.
    :param working_directory: Where 'file_name' exists locally.
    :param file_name: What to copy to CADC storage.
    :param archive: Which archive to associate the file with.
    :param stream: Defaults to raw - use is deprecated, however necessary it
        may be at the current moment to the 'put_file' call.
    :param mime_type: Because libmagic can't see inside a zipped fits file.
    :param mime_encoding: Also because libmagic can't see inside a zipped
        fits file.
    :param metrics: Tracking success execution times, and failure counts.
    """
    start = current()
    cwd = os.getcwd()
    try:
        os.chdir(working_directory)
        client.put_file(archive, file_name, archive_stream=stream,
                        mime_type=mime_type, mime_encoding=mime_encoding,
                        md5_check=True)
        file_size = os.stat(file_name).st_size
    except Exception as e:
        metrics.observe_failure('get', 'data', file_name)
        raise CadcException('Failed to store data with {}'.format(e))
    finally:
        os.chdir(cwd)
    end = current()
    metrics.observe(start, end, file_size, 'put', 'data', file_name)


def data_get(client, working_directory, file_name, archive, metrics):
    """
    Retrieve a local copy of a file available from CADC. Assumes the working
    directory location exists and is writeable.

    :param client: The CadcDataClient for read access to CADC storage.
    :param working_directory: Where 'file_name' will be written.
    :param file_name: What to copy from CADC storage.
    :param archive: Which archive to retrieve the file from.
    """
    start = current()
    fqn = os.path.join(working_directory, file_name)
    try:
        client.get_file(archive, file_name, destination=fqn)
        if not os.path.exists(fqn):
            raise CadcException(
                'ad retrieve failed. {} does not exist.'.format(fqn))
    except Exception as e:
        metrics.observe_failure('get', 'data', file_name)
        raise CadcException('Did not retrieve {} because {}'.format(
            fqn, e))
    end = current()
    file_size = os.stat(fqn).st_size
    metrics.observe(start, end, file_size, 'get', 'data', file_name)


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


def make_seconds(from_time):
    """Deal with different time formats to get the number of
    seconds since the epoch.

    Timezone information may be present as +00, strip that for returned
    results.
    :param from_time a string representing some time
    :return the time as a timestamp, so seconds.microseconds
    """
    try:
        index = from_time.index('+00')
    except ValueError:
        index = len(from_time)

    for fmt in [ISO_8601_FORMAT, '%Y-%m-%dT%H:%M:%S', '%d-%b-%Y %H:%M',
                '%b %d %Y', '%b %d %H:%M', '%Y%m%d-%H%M%S']:
        try:
            seconds_since_epoch = datetime.strptime(
                from_time[:index], fmt).timestamp()
            if fmt == '%b %d %H:%M':
                # the format '%b %d %H:%M' results in a timestamp based on
                # 1900, so need to set it to 'this' year
                year = datetime.utcnow().year
                dt = '{} {}'.format(from_time[:index], year)
                dt_format = '{} %Y'.format(fmt)
                seconds_since_epoch = datetime.strptime(
                    dt, dt_format).timestamp()
            break
        except ValueError:
            seconds_since_epoch = None
    if seconds_since_epoch is None:
        raise CadcException('Could not make seconds from {}'.format(from_time))
    return seconds_since_epoch


def increment_time(this_ts, by_interval, unit='%M'):
    """
    Increment time by an interval. Times should be in datetime format, but
    a modest attempt is made to check for otherwise.

    :param this_ts: datetime
    :param by_interval: integer - e.g. 10, for a 10 minute increment
    :param unit: the formatting string, default is minutes
    :return: this_ts incremented by interval amount
    """
    if isinstance(this_ts, datetime):
        time_s = this_ts.timestamp()
    elif isinstance(this_ts, str):
        time_s = make_seconds(this_ts)
    else:
        time_s = this_ts
    if unit == '%M':
        factor = 60
    else:
        raise NotImplementedError('Unexpected unit {}'.format(unit))
    interval = by_interval * factor
    temp = time_s + interval
    return datetime.fromtimestamp(temp)


def http_get(url, local_fqn):
    """Retrieve a file via http.

    :param url where the file can be found.
    :param local_fqn fully qualified name for where to file the file
        locally.
    """
    try:
        with requests.get(url, stream=True, timeout=10) as r:
            r.raise_for_status()
            with open(local_fqn, 'wb') as f:
                for chunk in r.iter_content(chunk_size=READ_BLOCK_SIZE):
                    f.write(chunk)
            length = to_int(r.headers.get('Content-Length'))
            if length is not None:
                file_meta = get_file_meta(local_fqn)
                if file_meta['size'] != length:
                    raise CadcException(
                        'Could not retrieve {} from {}. File size '
                        'error.'.format(local_fqn, url))
            checksum = r.headers.get('Content-Checksum')
            if checksum is not None:
                file_meta = get_file_meta(local_fqn)
                if file_meta['md5sum'] != checksum:
                    raise CadcException(
                        'Could not retrieve {} from {}. File checksum '
                        'error.'.format(local_fqn, url))
        if not os.path.exists(local_fqn):
            raise CadcException(
                'Retrieve failed. {} does not exist.'.format(local_fqn))
    except exceptions.HttpException as e:
        raise CadcException(
            'Could not retrieve {} from {}. Failed with {}'.format(
                local_fqn, url, e))


def look_pull_and_put(f_name, working_dir, url, archive, stream, mime_type,
                      cadc_client, checksum, metrics):
    """Checks to see if a file exists in ad. If yes, stop. If no,
    pull via https to local storage, then put to ad.

    TODO - stream

    :param f_name file name on disk for caching between the
        pull and the put
    :param working_dir together with f_name, location for caching
    :param url for retrieving the file externally, if it does not exist
    :param archive for storing in ad
    :param stream for storing in ad
    :param mime_type because libmagic is not always available
    :param cadc_client access to the data web service
    :param checksum what the CAOM observation says the checksum should be -
        just the checksum part of ChecksumURI please, or the comparison will
        always fail.
    :param metrics track how long operations take
    """
    retrieve = False
    try:
        meta = cadc_client.get_file_info(archive, f_name)
        if checksum is not None and meta['md5sum'] != checksum:
            logging.debug('Different checksums: CADC {} Source {}'.format(
                meta['md5sum'], checksum))
            retrieve = True
        else:
            logging.info('{} already exists at CADC/{}'.format(
                f_name, archive))
    except exceptions.NotFoundException:
        retrieve = True

    if retrieve:
        logging.info('Retrieving {} for {}'.format(f_name, archive))
        fqn = os.path.join(working_dir, f_name)
        http_get(url, fqn)
        data_put(cadc_client, working_dir, f_name, archive, stream,
                 mime_type, mime_encoding=None, metrics=metrics)


def repo_create(client, observation, metrics):
    start = current()
    try:
        client.create(observation)
    except Exception as e:
        metrics.observe_failure('create', 'caom2', observation.observation_id)
        raise CadcException(
            'Could not create an observation record for {}. '
            '{}'.format(observation.observation_id, e))
    end = current()
    metrics.observe(start, end, sizeof(observation), 'create', 'caom2',
                    observation.observation_id)


def repo_delete(client, collection, obs_id, metrics):
    start = current()
    try:
        client.delete(collection, obs_id)
    except Exception as e:
        metrics.observe_failure('delete', 'caom2', obs_id)
        raise CadcException(
            'Could not delete the observation record for {}. '
            '{}'.format(obs_id, e))
    end = current()
    metrics.observe(start, end, 0, 'delete', 'caom2', obs_id)


def repo_get(client, collection, obs_id, metrics):
    start = current()
    try:
        observation = client.read(collection, obs_id)
    except exceptions.NotFoundException:
        observation = None
    except Exception:
        metrics.observe_failure('read', 'caom2', obs_id)
        raise CadcException(
            'Could not retrieve an observation record for {}.'.format(obs_id))
    end = current()
    metrics.observe(
        start, end, sizeof(observation), 'read', 'caom2', obs_id)
    return observation


def repo_update(client, observation, metrics):
    start = current()
    try:
        client.update(observation)
    except Exception as e:
        metrics.observe_failure('update', 'caom2', observation.observation_id)
        raise CadcException(
            'Could not update an observation record for {}. '
            '{}'.format(observation.observation_id, e))
    end = current()
    metrics.observe(start, end, sizeof(observation), 'update', 'caom2',
                    observation.observation_id)


def find_missing(compare_this, to_this):
    return [ii for ii in compare_this if ii not in to_this]
