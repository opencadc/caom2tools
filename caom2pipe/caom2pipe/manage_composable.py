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

import logging
import os
import subprocess
import yaml

from aenum import Enum
from hashlib import md5
from os import stat

from cadcutils import net
from cadcdata import CadcDataClient


__all__ = ['CadcException', 'Config', 'to_float', 'TaskType',
           'exec_cmd', 'exec_cmd_redirect', 'exec_cmd_info',
           'get_cadc_meta', 'get_file_meta', 'compare_checksum']


class CadcException(Exception):
    pass


class TaskType(Enum):
    """The possible steps in the OMM pipeline."""
    STORE = 'store'
    SCRAPE = 'scrape'
    INGEST = 'ingest'
    MODIFY = 'modify'
    UNKNOWN = 'unknown'


class Config(object):
    """Configuration information that remains the same for all steps and all
     work in a pipeline execution."""

    def __init__(self):
        # the root directory for all executor operations
        self.working_directory = None

        # the file that contains the list of work to be passed through the
        # pipeline
        self.work_file = None
        # the fully qualified name for the work file
        self.work_fqn = None

        # credentials for any service calls
        self.netrc_file = None

        # which collection is addressed by the pipeline
        self.collection = None

        # changes expectations of the executors for handling files on disk
        self.use_local_files = False

        # which service instance to use
        self.resource_id = None

        # the logging level - enforced throughout the pipeline
        self.logging_level = None

        # write the log to a file?
        self.log_to_file = False

        # where log files are written to - defaults to working_directory
        self.log_file_directory = None

        # the ad 'stream' that goes with the collection - use when storing
        # files
        self.stream = None

        # the ad 'host' to store files to - used for testing cadc-data put
        # commands only, should usually be None
        self.storage_host = None

        # the way to control which steps get executed
        self.task_types = None

        # the filename where success logs are written
        # this will be created in log_file_directory
        self.success_log_file_name = None
        # the fully qualified name for the file
        self.success_fqn = None

        # the filename where failure logs are written
        # this will be created in log_file_directory
        self.failure_log_file_name = None
        # the fully qualified name for the file
        self.failure_fqn = None

        # the filename where retry entries are writter
        # this will be created in log_file_directory
        self.retry_file_name = None
        # the fully qualified name for the file
        self.retry_fqn = None

    @property
    def working_directory(self):
        return self._working_directory

    @working_directory.setter
    def working_directory(self, value):
        self._working_directory = value

    @property
    def work_file(self):
        return self._work_file

    @work_file.setter
    def work_file(self, value):
        self._work_file = value
        if self.working_directory is not None:
            self.work_fqn = os.path.join(
                self.working_directory, self.work_file)

    @property
    def netrc_file(self):
        return self._netrc_file

    @netrc_file.setter
    def netrc_file(self, value):
        self._netrc_file = value

    @property
    def collection(self):
        return self._collection

    @collection.setter
    def collection(self, value):
        self._collection = value

    @property
    def use_local_files(self):
        return self._use_local_files

    @use_local_files.setter
    def use_local_files(self, value):
        self._use_local_files = value

    @property
    def resource_id(self):
        return self._resource_id

    @resource_id.setter
    def resource_id(self, value):
        self._resource_id = value

    @property
    def log_to_file(self):
        return self._log_to_file

    @log_to_file.setter
    def log_to_file(self, value):
        self._log_to_file = value

    @property
    def log_file_directory(self):
        return self._log_file_directory

    @log_file_directory.setter
    def log_file_directory(self, value):
        self._log_file_directory = value

    @property
    def logging_level(self):
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
        return self._stream

    @stream.setter
    def stream(self, value):
        self._stream = value

    @property
    def storage_host(self):
        return self._storage_host

    @storage_host.setter
    def storage_host(self, value):
        self._storage_host = value

    @property
    def task_type(self):
        return self._task_type

    @task_type.setter
    def task_type(self, value):
        self._task_type = value

    @property
    def success_log_file_name(self):
        return self._success_log_file_name

    @success_log_file_name.setter
    def success_log_file_name(self, value):
        self._success_log_file_name = value
        if self.log_file_directory is not None:
            self.success_fqn = os.path.join(
                self.log_file_directory, self.success_log_file_name)

    @property
    def failure_log_file_name(self):
        return self._failure_log_file_name

    @failure_log_file_name.setter
    def failure_log_file_name(self, value):
        self._failure_log_file_name = value
        if self.log_file_directory is not None:
            self.failure_fqn = os.path.join(
                self.log_file_directory, self.failure_log_file_name)

    @property
    def retry_file_name(self):
        return self._retry_file_name

    @retry_file_name.setter
    def retry_file_name(self, value):
        self._retry_file_name = value
        if self.log_file_directory is not None:
            self.retry_fqn = os.path.join(
                self.log_file_directory, self.retry_file_name)

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
               'collection:: \'{}\' ' \
               'task_types:: \'{}\' ' \
               'stream:: \'{}\' ' \
               'resource_id:: \'{}\' ' \
               'use_local_files:: \'{}\' ' \
               'log_to_file:: \'{}\' ' \
               'log_file_directory:: \'{}\' ' \
               'success_log_file_name:: \'{}\' ' \
               'success_fqn:: \'{}\' ' \
               'failure_log_file_name:: \'{}\' ' \
               'failure_fqn:: \'{}\' ' \
               'retry_file_name:: \'{}\' ' \
               'retry_fqn:: \'{}\' ' \
               'retry_count:: \'{}\' ' \
               'logging_level:: \'{}\''.format(
                self.working_directory, self.work_fqn, self.netrc_file,
                self.collection, self.task_types, self.stream,
                self.resource_id, self.use_local_files, self.log_to_file,
                self.log_file_directory, self.success_log_file_name,
                self.success_fqn, self.failure_log_file_name,
                self.failure_fqn, self.retry_file_name, self.retry_fqn,
                self.retry_count, self.logging_level)

    @staticmethod
    def _set_task_types(config, default=None):
        task_types = []
        if 'task_types' in config:
            for ii in config['task_types']:
                task_types.append(TaskType(ii))
            return task_types
        else:
            return default

    def get_executors(self):
        """Look up the configuration values in the data structure extracted
        from the configuration file."""
        try:
            config = self.get_config()
            self.working_directory = \
                self._lookup(config, 'working_directory', os.getcwd())
            self.work_file = self._lookup(config, 'todo_file_name', 'todo.txt')
            self.netrc_file = \
                self._lookup(config, 'netrc_filename', 'test_netrc')
            self.resource_id = self._lookup(
                config, 'resource_id', 'ivo://cadc.nrc.ca/sc2repo')
            self.use_local_files = bool(
                self._lookup(config, 'use_local_files', False))
            self.logging_level = self._lookup(config, 'logging_level', 'DEBUG')
            self.log_to_file = self._lookup(config, 'log_to_file', False)
            self.log_file_directory = self._lookup(config, 'log_file_directory',
                                                   self.working_directory)
            self.stream = self._lookup(config, 'stream', 'raw')
            self.task_types = self._set_task_types(config, [TaskType.SCRAPE])
            self.collection = self._lookup(config, 'collection', 'TEST')
            self.success_log_file_name = self._lookup(config,
                                                      'success_log_file_name',
                                                      'success_log.txt')
            self.failure_log_file_name = self._lookup(config,
                                                      'failure_log_file_name',
                                                      'failure_log.txt')
            self.retry_file_name = self._lookup(config, 'retry_file_name',
                                                'retries.txt')
            self.retry_count = self._lookup(config, 'retry_count', 0)
        except KeyError as e:
            raise CadcException(
                'Error in config file {}'.format(e))

    def get_config(self):
        """Return a configuration dictionary. Assumes a file named config.yml
        in the current working directory."""
        config_fqn = os.path.join(os.getcwd(), 'config.yml')
        config = self.load_config(config_fqn)
        if config is None:
            raise CadcException('Could not find the file {}'.format(config_fqn))
        return config

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
        raise CadcException('Could not execute cmd {}'.format(cmd))


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
        raise CadcException('Could not execute cmd {}'.format(cmd))


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
                cmd_array, stdout=outfile, stderr=subprocess.PIPE).communicate()
            if outerr is not None and len(outerr) > 0 and outerr[0] is not None:
                logging.debug('Command {} had stderr {}'.format(
                    cmd, outerr.decode('utf-8')))
                raise CadcException(
                    'Command {} had outerr {}'.format(
                        cmd, outerr.decode('utf-8')))
    except Exception as e:
        logging.debug('Error with command {}:: {}'.format(cmd, e))
        raise CadcException('Could not execute cmd {}'.format(cmd))


def get_cadc_meta(netrc_fqn, collection, fname):
    """
    Gets contentType, contentLength and contentChecksum of a CADC artifact
    :param netrc_fqn: user credentials
    :param collection: archive file has been stored to
    :param fname: name of file in the archive
    :return:
    """
    subject = net.Subject(username=None, certificate=None, netrc=netrc_fqn)
    client = CadcDataClient(subject)
    return client.get_file_info(collection, fname)


def get_file_meta(fqn):
    """
    Gets contentType, contentLength and contentChecksum of an artifact on disk.

    :param fqn: Fully-qualified name of the file for which to get the metadata.
    :return:
    """
    meta = {}
    s = stat(fqn)
    meta['size'] = s.st_size
    meta['md5sum'] = md5(open(fqn, 'rb').read()).hexdigest()
    if fqn.endswith('.header') or fqn.endswith('.txt'):
        meta['type'] = 'text/plain'
    elif fqn.endswith('.gif'):
        meta['type'] = 'image/gif'
    elif fqn.endswith('.png'):
        meta['type'] = 'image/png'
    else:
        meta['type'] = 'application/octet-stream'
    return meta


def compare_checksum(netrc_fqn, collection, fqn):
    """
    Raise CadcException if the checksum of a file in ad is not the same as
    the checksum of a file on disk.

    :param netrc_fqn: fully-qualified file name for the netrc file
    :param collection: archive file has been stored to
    :param fqn: Fully-qualified name of the file for which to get the metadata.
    """
    fname = os.path.basename(fqn)
    try:
        local_meta = get_file_meta(fqn)
        ad_meta = get_cadc_meta(netrc_fqn, collection, fname)
    except Exception as e:
        raise CadcException('Could not find md5 checksum for {} in the ad {} '
                            'collection. {}'.format(fqn, collection, e))

    if ((fqn.endswith('.gz') and local_meta['md5sum'] !=
         ad_meta['md5sum']) or (
            not fqn.endswith('.gz') and local_meta['md5sum'] !=
            ad_meta['umd5sum'])):
        raise CadcException(
            '{} md5sum not the same as the one in the ad '
            '{} collection.'.format(fqn, collection))
