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

"""
This module contains pipeline execution classes. Each execution class
corresponds to a single task type, correlated with specific configuration
or implementation assumptions.

The execute methods in each of the class definitions require no if statements.
All the if statements are limited to the choose* methods in the
OrganizeExecutes class. If you find yourself adding an if statement to an
execute method, create a new *Execute class instead. The result is execute
methods that are composable into complex and varied pipelines, while
remaining easily tested.

Raise the CadcException upon encountering an error. There is no recovery
effort as part of a failure. Log the error and stop the pipeline
execution for an Observation.

The correlations that currently exist:
- use_local_data: True => classes have "Local" in their name
- uses the CadcDataClient, the Caom2RepoClient => classes have "Client"
    in their name
- requires metadata access only => classes have "Meta" in their name
- requires data access => classes have "Data" in their name

"""

import distutils.sysconfig
import logging
import os
import re
import sys
import traceback

from argparse import ArgumentParser
from datetime import datetime

from cadcutils import net, exceptions
from cadcdata import CadcDataClient
from caom2repo import CAOM2RepoClient
from caom2pipe import manage_composable as mc


__all__ = ['OrganizeExecutes', 'StorageName', 'CaomName']


class StorageName(object):
    """Naming rules for a collection:
    - support mixed-case file name storage
    - support gzipped and not zipped file names

    This class assumes the obs_id is part of the file name. This assumption
    may be broken in the future, in which case lots of CaomExecute
    implementations will need to be re-addressed somehow.

    This class assumes the file name in storage, and the file name on disk
    are not necessarily the same thing.
    """

    def __init__(self, obs_id, collection, collection_pattern,
                 fname_on_disk=None):
        """

        :param obs_id: string value for Observation.observationID
        :param collection: string value for Observation.collection
        :param collection_pattern: regular expression that can be used to
            determine if a file name or observation id meets particular
            patterns.
        :param fname_on_disk: string value for the name of a file on disk,
            which is not necessarily the same thing as the name of the file
            in storage (i.e. extensions may exist in one location that do
            not exist in another.
        """
        self.obs_id = obs_id
        self.collection = collection
        self.collection_pattern = collection_pattern
        self.fname_on_disk = fname_on_disk

    @property
    def file_uri(self):
        """The ad URI for the file. Assumes compression."""
        return 'ad:{}/{}.gz'.format(self.collection, self.file_name)

    @property
    def file_name(self):
        """The file name."""
        return '{}.fits'.format(self.obs_id)

    @property
    def compressed_file_name(self):
        """The compressed file name - adds the .gz extension."""
        return '{}.fits.gz'.format(self.obs_id)

    @property
    def model_file_name(self):
        """The file name used on local disk that holds the CAOM2 Observation
        XML."""
        return '{}.fits.xml'.format(self.obs_id)

    @property
    def prev(self):
        """The preview file name for the file."""
        return '{}_prev.jpg'.format(self.obs_id)

    @property
    def thumb(self):
        """The thumbnail file name for the file."""
        return '{}_prev_256.jpg'.format(self.obs_id)

    @property
    def prev_uri(self):
        """The ad preview URI."""
        return self._get_uri(self.prev)

    @property
    def thumb_uri(self):
        """The ad thumbnail URI."""
        return self._get_uri(self.thumb)

    @property
    def obs_id(self):
        """The observation ID associated with the file name."""
        return self._obs_id

    @obs_id.setter
    def obs_id(self, value):
        self._obs_id = value

    @property
    def log_file(self):
        """The log file name used when running any of the 'execute' steps."""
        return '{}.log'.format(self.obs_id)

    @property
    def product_id(self):
        """The relationship between the observation ID of an observation, and
        the product ID of a plane."""
        return self.obs_id

    @property
    def fname_on_disk(self):
        """The file name on disk, which is not necessarily the same as the
        file name in ad."""
        return self._fname_on_disk

    @fname_on_disk.setter
    def fname_on_disk(self, value):
        self._fname_on_disk = value

    def is_valid(self):
        """:return True if the observation ID conforms to naming rules."""
        pattern = re.compile(self.collection_pattern)
        return pattern.match(self.obs_id)

    def _get_uri(self, fname):
        """The ad URI for a file, without consideration for compression."""
        return 'ad:{}/{}'.format(self.collection, fname)

    @staticmethod
    def remove_extensions(name):
        """How to get the file_id from a file_name."""
        return name.replace('.fits', '').replace('.gz', '').replace('.header',
                                                                    '')


class CaomName(object):
    """The naming rules for making and decomposing CAOM URIs (i.e. Observation
    URIs, Plane URIs, and archive URIs, all isolated in one class. There are
    probably OMM assumptions built in, but those will slowly go away :). """

    def __init__(self, uri):
        self.uri = uri

    @property
    def file_id(self):
        """

        :return: Extracted from an Artifact URI, the file_id is the file
        name portion of the URI with all file type and compression type
        extensions removed.
        """
        return self.uri.split('/')[1].split('.')[0]

    @property
    def file_name(self):
        """:return The file name extracted from an Artifact URI."""
        return self.uri.split('/')[1]

    @property
    def uncomp_file_name(self):
        """:return The file name extracted from an Artifact URI, without
        the compression extension."""
        return self.file_name.replace('.gz', '')

    @staticmethod
    def make_obs_uri_from_obs_id(collection, obs_id):
        """:return A string that conforms to the Observation URI
        specification from CAOM."""
        return 'caom:{}/{}'.format(collection, obs_id)


class CaomExecute(object):
    """Abstract class that defines the operations common to all Execute
    classes."""

    def __init__(self, config, task_type, storage_name, command_name,
                 cred_param, cadc_data_client, caom_repo_client,
                 meta_visitors):
        """

        :param config: Configurable parts of execution, as stored in
            manage_composable.Config.
        :param task_type: manage_composable.TaskType enumeration - identifies
            the work to do, in words that are user-facing. Used in logging
            messages.
        :param storage_name: An instance of StorageName.
        :param command_name: The collection-specific application to apply a
            blueprint. May be 'fits2caom2'.
        :param cred_param: either --netrc <value> or --cert <value>,
            depending on which credentials have been supplied to the
            process.
        :param cadc_data_client: Instance of CadcDataClient. Used for data
            service access.
        :param caom_repo_client: Instance of CAOM2Repo client. Used for
            caom2 repository service access.
        :param meta_visitors: List of classes with a
            'visit(observation, **kwargs)' method signature. Requires access
            to metadata only.
        """
        self.logger = logging.getLogger()
        self.logger.setLevel(config.logging_level)
        formatter = logging.Formatter(
            '%(asctime)s:%(levelname)s:%(name)-12s:%(lineno)d:%(message)s')
        for handler in self.logger.handlers:
            handler.setLevel(config.logging_level)
            handler.setFormatter(formatter)
        self.logging_level_param = self._set_logging_level_param(
            config.logging_level)
        self.obs_id = storage_name.obs_id
        self.product_id = storage_name.product_id
        self.uri = storage_name.file_uri
        self.fname = storage_name.file_name
        self.command_name = command_name
        self.root_dir = config.working_directory
        self.collection = config.collection
        self.working_dir = os.path.join(self.root_dir, self.obs_id)
        if config.log_to_file:
            self.model_fqn = os.path.join(config.log_file_directory,
                                          storage_name.model_file_name)
        else:
            self.model_fqn = os.path.join(self.working_dir,
                                          storage_name.model_file_name)
        self.resource_id = config.resource_id
        self.cadc_data_client = cadc_data_client
        self.caom_repo_client = caom_repo_client
        self.stream = config.stream
        self.meta_visitors = meta_visitors
        self.task_type = task_type
        self.cred_param = cred_param

    def _cleanup(self):
        """Remove a directory and all its contents."""
        if os.path.exists(self.working_dir):
            for ii in os.listdir(self.working_dir):
                os.remove(os.path.join(self.working_dir, ii))
            os.rmdir(self.working_dir)

    def _create_dir(self):
        """Create the working area if it does not already exist."""
        mc.create_dir(self.working_dir)

    def _define_local_dirs(self, storage_name):
        """when files are on disk don't worry about a separate directory
        per observation"""
        self.working_dir = self.root_dir
        self.model_fqn = os.path.join(self.working_dir,
                                      storage_name.model_file_name)

    def _find_fits2caom2_plugin(self):
        """Find the code that is passed as the --plugin parameter to
        fits2caom2.

        This code makes the assumption that execution always occurs within
        the context of a Docker container, and therefore the
        get_python_lib call will always have the appropriately-named
        module installed in a site package location.
        """
        packages = distutils.sysconfig.get_python_lib()
        return os.path.join(packages, '{}/{}.py'.format(self.command_name,
                                                        self.command_name))

    def _fits2caom2_cmd_local(self):
        """Execute fits2caom with a --local parameter."""
        fqn = os.path.join(self.working_dir, self.fname)
        plugin = self._find_fits2caom2_plugin()
        # so far, the plugin is also the module :)
        cmd = '{} {} {} --observation {} {} --out {} ' \
              '--plugin {} --module {} --local {} --lineage {}/{}'.format(
                self.command_name,
                self.logging_level_param, self.cred_param, self.collection,
                self.obs_id, self.model_fqn, plugin, plugin, fqn, self.obs_id,
                self.uri)
        mc.exec_cmd(cmd)

    def _fits2caom2_cmd_client(self):
        """Execute fits2caom with a --cert parameter."""
        plugin = self._find_fits2caom2_plugin()
        # so far, the plugin is also the module :)
        cmd = '{} {} {} --observation {} {} --out {} ' \
              '--plugin {} --module {} --lineage {}/{}'.format(
                self.command_name, self.logging_level_param, self.cred_param,
                self.collection, self.obs_id, self.model_fqn, plugin, plugin,
                self.product_id, self.uri)
        mc.exec_cmd(cmd)

    def _fits2caom2_cmd_in_out_client(self):
        """Execute fits2caom with a --local and a --cert parameter."""
        plugin = self._find_fits2caom2_plugin()
        # so far, the plugin is also the module :)
        # TODO add an input parameter
        cmd = '{} {} {} --in {} --out {} ' \
              '--plugin {} --module {} --lineage {}/{}'.format(
                self.command_name, self.logging_level_param, self.cred_param,
                self.model_fqn, self.model_fqn, plugin, plugin,
                self.product_id, self.uri)
        mc.exec_cmd(cmd)

    def _compare_checksums_client(self, fname):
        """Compare the checksum of a file on disk with a file in ad,
        using the client instance from this class."""
        fqn = os.path.join(self.working_dir, fname)
        mc.compare_checksum_client(
            self.cadc_data_client, self.collection, fqn)

    def _repo_cmd_create_client(self, observation):
        """Create an observation instance from the input parameter."""
        try:
            self.caom_repo_client.create(observation)
        except Exception as e:
            raise mc.CadcException(
                'Could not create an observation record for {} in {}. '
                '{}'.format(self.obs_id, self.resource_id, e))

    def _repo_cmd_update_client(self, observation):
        """Update an existing observation instance.  Assumes the obs_id
        values are set correctly."""
        try:
            self.caom_repo_client.update(observation)
        except Exception as e:
            raise mc.CadcException(
                'Could not update an observation record for {} in {}. '
                '{}'.format(self.obs_id, self.resource_id, e))

    def _repo_cmd_read_client(self):
        """Retrieve the existing observation model metadata."""
        try:
            return self.caom_repo_client.read(self.collection, self.obs_id)
        except Exception as e:
            raise mc.CadcException(
                'Could not read observation record for {} in {}. {}'.format(
                    self.obs_id, self.resource_id, e))

    def _cadc_data_put_client(self, fname):
        """Store a collection file."""
        try:
            self.cadc_data_client.put_file(self.collection, fname, self.stream)
        except Exception as e:
            raise mc.CadcException(
                'Did not store {} with {}'.format(fname, e))

    def _cadc_data_get_client(self):
        """Retrieve a collection file, even if it already exists. This might
        ensure that the latest version of the file is retrieved from
        storage."""

        # all the gets are supposed to be unzipped, so that the
        # footprintfinder is more efficient, so make sure the fully
        # qualified output name isn't the gzip'd version

        if self.fname.endswith('.gz'):
            fqn = os.path.join(self.working_dir, self.fname.replace('.gz', ''))
        else:
            fqn = os.path.join(self.working_dir, self.fname)

        try:
            self.cadc_data_client.get_file(
                self.collection, self.fname, destination=fqn, decompress=True)
            if not os.path.exists(fqn):
                raise mc.CadcException(
                    '{} does not exist.'.format(fqn))
        except Exception:
            raise mc.CadcException(
                'Did not retrieve {}'.format(fqn))

    def _cadc_data_info_file_name_client(self):
        """Execute CadcDataClient.get_file_info with the client instance from
        this class."""
        file_info = self.cadc_data_client.get_file_info(
            self.collection, self.fname)
        self.fname = file_info['name']

    def _read_model(self):
        """Read an observation into memory from an XML file on disk."""
        return mc.read_obs_from_file(self.model_fqn)

    def _write_model(self, observation):
        """Write an observation to disk from memory, represented in XML."""
        mc.write_obs_to_file(observation, self.model_fqn)

    def _visit_meta(self, observation):
        """Execute metadata-only visitors on an Observation in
        memory."""
        if self.meta_visitors is not None and len(self.meta_visitors) > 0:
            kwargs = {'working_directory': self.working_dir}
            for visitor in self.meta_visitors:
                try:
                    self.logger.debug('Visit for {}'.format(visitor))
                    visitor.visit(observation, **kwargs)
                except Exception as e:
                    raise mc.CadcException(e)

    @staticmethod
    def _set_logging_level_param(logging_level):
        """Make a configured logging level into command-line parameters."""
        lookup = {logging.DEBUG: '--debug',
                  logging.INFO: '--verbose',
                  logging.WARNING: '',
                  logging.ERROR: '--quiet'}
        if logging_level in lookup:
            result = lookup[logging_level]
        else:
            result = ''
        return result

    @staticmethod
    def repo_cmd_get_client(caom_repo_client, collection, observation_id):
        """Execute the CAOM2Repo 'read' operation using the client instance
        from this class.
        :return an Observation instance, or None, if the observation id
        does not exist."""
        try:
            observation = caom_repo_client.read(collection, observation_id)
            return observation
        except exceptions.NotFoundException:
            return None
        except Exception:
            raise mc.CadcException(
                'Could not retrieve an observation record for {}.'.format(
                    observation_id))


class Collection2CaomMetaCreateClient(CaomExecute):
    """Defines the pipeline step for Collection ingestion of metadata into CAOM.
    This requires access to only header information.

    This pipeline step will execute a caom2-repo create."""

    def __init__(self, config, storage_name, command_name,
                 cred_param, cadc_data_client, caom_repo_client,
                 meta_visitors):
        super(Collection2CaomMetaCreateClient, self).__init__(
            config, mc.TaskType.INGEST, storage_name, command_name,
            cred_param, cadc_data_client, caom_repo_client, meta_visitors)

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('Find the file name as stored.')
        self._cadc_data_info_file_name_client()

        self.logger.debug('create the work space, if it does not exist')
        self._create_dir()

        self.logger.debug('the observation does not exist, so go '
                          'straight to generating the xml, as the main_app '
                          'will retrieve the headers')
        self._fits2caom2_cmd_client()

        self.logger.debug('read the xml into memory from the file')
        observation = self._read_model()

        self.logger.debug('the metadata visitors')
        self._visit_meta(observation)

        self.logger.debug('store the xml')
        self._repo_cmd_create_client(observation)

        self.logger.debug('clean up the workspace')
        self._cleanup()

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomMetaUpdateClient(CaomExecute):
    """Defines the pipeline step for Collection ingestion of metadata into CAOM.
    This requires access to only header information.

    This pipeline step will execute a caom2-repo update."""

    def __init__(self, config, storage_name, command_name, cred_param,
                 cadc_data_client, caom_repo_client, observation,
                 meta_visitors):
        super(Collection2CaomMetaUpdateClient, self).__init__(
            config, mc.TaskType.INGEST, storage_name, command_name, cred_param,
            cadc_data_client, caom_repo_client, meta_visitors)
        self.observation = observation

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('Find the file name as stored.')
        self._cadc_data_info_file_name_client()

        self.logger.debug('create the work space, if it does not exist')
        self._create_dir()

        self.logger.debug('write the observation to disk for next step')
        self._write_model(self.observation)

        self.logger.debug('generate the xml, as the main_app will retrieve '
                          'the headers')
        self._fits2caom2_cmd_in_out_client()

        self.logger.debug('read the xml from disk')
        self.observation = self._read_model()

        self.logger.debug('the metadata visitors')
        self._visit_meta(self.observation)

        self.logger.debug('store the xml')
        self._repo_cmd_update_client(self.observation)

        self.logger.debug('clean up the workspace')
        self._cleanup()

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomLocalMetaCreateClient(CaomExecute):
    """Defines the pipeline step for Collection ingestion of metadata into CAOM.
    This requires access to only header information.

    This pipeline step will execute a caom2-repo create."""

    def __init__(self, config, storage_name, command_name, cred_param,
                 cadc_data_client, caom_repo_client, meta_visitors):
        super(Collection2CaomLocalMetaCreateClient, self).__init__(
            config, mc.TaskType.INGEST, storage_name, command_name, cred_param,
            cadc_data_client, caom_repo_client, meta_visitors)
        self._define_local_dirs(storage_name)

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('Find the file name as stored.')
        self._cadc_data_info_file_name_client()

        self.logger.debug('the observation does not exist, so go '
                          'straight to generating the xml, as the main_app '
                          'will retrieve the headers')
        self._fits2caom2_cmd_client()

        self.logger.debug('read the xml from disk')
        observation = self._read_model()

        self.logger.debug('the metadata visitors')
        self._visit_meta(observation)

        self.logger.debug('store the xml')
        self._repo_cmd_create_client(observation)

        self.logger.debug('write the updated xml to disk for debugging')
        self._write_model(observation)

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomLocalMetaUpdateClient(CaomExecute):
    """Defines the pipeline step for Collection ingestion of metadata into CAOM.
    This requires access to only header information.

    This pipeline step will execute a caom2-repo update."""

    def __init__(self, config, storage_name, command_name, cred_param,
                 cadc_data_client, caom_repo_client, observation,
                 meta_visitors):
        super(Collection2CaomLocalMetaUpdateClient, self).__init__(
            config, mc.TaskType.INGEST, storage_name, command_name, cred_param,
            cadc_data_client, caom_repo_client, meta_visitors)
        self._define_local_dirs(storage_name)
        self.observation = observation

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('Find the file name as stored.')
        self._cadc_data_info_file_name_client()

        self.logger.debug('write the observation to disk for next step')
        self._write_model(self.observation)

        self.logger.debug('generate the xml, as the main_app will retrieve '
                          'the headers')
        self._fits2caom2_cmd_in_out_client()

        self.logger.debug('read the xml from disk')
        self.observation = self._read_model()

        self.logger.debug('the metadata visitors')
        self._visit_meta(self.observation)

        self.logger.debug('store the xml')
        self._repo_cmd_update_client(self.observation)

        self.logger.debug('write the updated xml to disk for debugging')
        self._write_model(self.observation)

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomClientVisit(CaomExecute):
    """Defines the pipeline step for Collection augmentation by a visitor
     of metadata into CAOM. This assumes a record already exists in CAOM,
     and the update DOES NOT require access to either the header or the data.

    This pipeline step will execute a caom2-repo update."""

    def __init__(self, config, storage_name, cred_param,
                 cadc_data_client, caom_repo_client, meta_visitors):
        super(Collection2CaomClientVisit, self).__init__(
            config, mc.TaskType.VISIT, storage_name, command_name=None,
            cred_param=cred_param, cadc_data_client=cadc_data_client,
            caom_repo_client=caom_repo_client,
            meta_visitors=meta_visitors)
        self.fname = None

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')

        # TODO - run a test to see if this is necessary
        # self.logger.debug('Find the file name as stored.')
        # self._find_file_name_storage_client()

        self.logger.debug('retrieve the existing observation, if it exists')
        observation = self._repo_cmd_read_client()

        self.logger.debug('the metadata visitors')
        self._visit_meta(observation)

        self.logger.debug('store the xml')
        self._repo_cmd_update_client(observation)

        self.logger.debug('clean up the workspace')
        self._cleanup()

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomDataClient(CaomExecute):
    """Defines the pipeline step for all the operations that
    require access to the file on disk, not just the header data. """

    def __init__(self, config, storage_name, command_name, cred_param,
                 cadc_data_client, caom_repo_client, data_visitors,
                 task_type):
        super(Collection2CaomDataClient, self).__init__(
            config, task_type, storage_name, command_name, cred_param,
            cadc_data_client, caom_repo_client, meta_visitors=None)
        self.log_file_directory = config.log_file_directory
        self.data_visitors = data_visitors
        self.prev_fname = storage_name.prev
        self.thumb_fname = storage_name.thumb

    def execute(self, context):
        self.logger.debug('Begin execute for {} Data'.format(__name__))

        self.logger.debug('Find the file name as stored.')
        self._cadc_data_info_file_name_client()

        self.logger.debug('create the work space, if it does not exist')
        self._create_dir()

        self.logger.debug('get the input file')
        self._cadc_data_get_client()

        self.logger.debug('get the observation for the existing model')
        observation = self._repo_cmd_read_client()

        self.logger.debug('execute the data visitors')
        self._visit_data(observation)

        self.logger.debug('store the updated xml')
        self._repo_cmd_update_client(observation)

        self.logger.debug('clean up the workspace')
        self._cleanup()

        self.logger.debug('End execute for {}'.format(__name__))

    def _visit_data(self, observation):
        """Execute the visitors that require access to the full data content
        of a file."""
        kwargs = {'working_directory': self.working_dir,
                  'science_file': self.fname,
                  'log_file_directory': self.log_file_directory,
                  'cadc_client': self.cadc_data_client}
        for visitor in self.data_visitors:
            try:
                self.logger.debug('Visit for {}'.format(visitor))
                visitor.visit(observation, **kwargs)
            except Exception as e:
                raise mc.CadcException(e)


class Collection2CaomLocalDataClient(Collection2CaomDataClient):
    """Defines the pipeline step for all the operations that
    require access to the file on disk. This class assumes it has access to
    the files on disk."""

    def __init__(self, config, storage_name, command_name, cred_param,
                 cadc_data_client, caom_repo_client, data_visitors):
        super(Collection2CaomLocalDataClient, self).__init__(
            config, storage_name, command_name, cred_param,
            cadc_data_client=cadc_data_client,
            caom_repo_client=caom_repo_client, data_visitors=data_visitors,
            task_type=mc.TaskType.MODIFY)
        self._define_local_dirs(storage_name)
        self.fname = storage_name.fname_on_disk

    def execute(self, context):
        self.logger.debug('Begin execute for {} Data'.format(__name__))

        self.logger.debug('get the observation for the existing model')
        observation = self._repo_cmd_read_client()

        self.logger.debug('execute the data visitors')
        self._visit_data(observation)

        self.logger.debug('store the updated xml')
        self._repo_cmd_update_client(observation)

        self.logger.debug('write the updated xml to disk for debugging')
        self._write_model(observation)

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomStoreClient(CaomExecute):
    """Defines the pipeline step for Collection storage of a file. This requires
    access to the file on disk."""

    def __init__(self, config, storage_name, command_name, cred_param,
                 cadc_data_client, caom_repo_client):
        super(Collection2CaomStoreClient, self).__init__(
            config, mc.TaskType.STORE, storage_name, command_name, cred_param,
            cadc_data_client, caom_repo_client, meta_visitors=None)
        # when files are on disk don't worry about a separate directory
        # per observation
        self.working_dir = self.root_dir
        self.storage_host = config.storage_host
        self.stream = config.stream
        self.fname = storage_name.fname_on_disk

    def execute(self, context):
        self.logger.debug('Begin execute for {} Data'.format(__name__))

        self.logger.debug('store the input file to ad')
        self._cadc_data_put_client(self.fname)

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomScrape(CaomExecute):
    """Defines the pipeline step for Collection creation of a CAOM model
    observation. The file containing the metadata is located on disk.
    No record is written to a web service."""

    def __init__(self, config, storage_name, command_name):
        super(Collection2CaomScrape, self).__init__(
            config, mc.TaskType.SCRAPE, storage_name, command_name,
            cred_param='', cadc_data_client=None, caom_repo_client=None,
            meta_visitors=None)
        self._define_local_dirs(storage_name)
        self.fname = storage_name.fname_on_disk

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('generate the xml from the file on disk')
        self._fits2caom2_cmd_local()

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomDataScrape(Collection2CaomDataClient):
    """Defines the pipeline step for Collection generation and ingestion of
    operations that require access to the file on disk, with no update to the
    service at the end. This class assumes it has access to the files on disk.
    The organization of this class assumes the 'Scrape' task has been done
    previously, so the model instance exists on disk."""

    def __init__(self, config, storage_name, command_name, data_visitors):
        super(Collection2CaomDataScrape, self).__init__(
            config, storage_name, command_name, cred_param='',
            cadc_data_client=None, caom_repo_client=None,
            data_visitors=data_visitors, task_type=mc.TaskType.SCRAPE)
        self._define_local_dirs(storage_name)
        self.fname = storage_name.fname_on_disk
        self.log_file_directory = config.log_file_directory
        self.data_visitors = data_visitors
        self.prev_fname = storage_name.prev
        self.thumb_fname = storage_name.thumb

    def execute(self, context):
        self.logger.debug('Begin execute for {} Data'.format(__name__))

        self.logger.debug('get observation for the existing model from disk')
        observation = self._read_model()

        self.logger.debug('execute the data visitors')
        self._visit_data(observation)

        self.logger.debug('output the updated xml')
        self._write_model(observation)

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomCompareChecksumClient(CaomExecute):
    """Defines the pipeline step for comparing the checksum of a file on disk
    with the checksum of the supposedly-the-same file stored at CADC.

    This step should be invoked with any other task type that relies on
    files on local disk.
    """

    def __init__(self, config, storage_name, command_name, cred_param,
                 cadc_data_client, caom_repo_client):
        super(Collection2CaomCompareChecksumClient, self).__init__(
            config, mc.TaskType.CHECKSUM, storage_name, command_name,
            cred_param, cadc_data_client, caom_repo_client,
            meta_visitors=None)
        self._define_local_dirs(storage_name)
        self.fname = storage_name.fname_on_disk

    def execute(self, context):
        self.logger.debug('Begin execute for {} '
                          'CompareChecksum'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('generate the xml from the file on disk')
        self._compare_checksums_client(self.fname)

        self.logger.debug('End execute for {}'.format(__name__))


class OrganizeExecutes(object):
    """How to turn on/off various task types in the a CaomExecute pipeline."""
    def __init__(self, config, todo_file=None):
        self.config = config
        self.task_types = config.task_types
        self.logger = logging.getLogger()
        self.logger.setLevel(config.logging_level)
        if todo_file is not None:
            self.todo_fqn = todo_file
            todo_name = os.path.basename(todo_file).split('.')[0]
            self.success_fqn = os.path.join(
                self.config.log_file_directory,
                '{}_success_log.txt'.format(todo_name))
            config.success_fqn = self.success_fqn
            self.failure_fqn = os.path.join(
                self.config.log_file_directory,
                '{}_failure_log.txt'.format(todo_name))
            config.failure_fqn = self.failure_fqn
            self.retry_fqn = os.path.join(
                self.config.log_file_directory,
                '{}_retries.txt'.format(todo_name))
            config.retry_fqn = self.retry_fqn
        else:
            self.todo_fqn = config.work_fqn
            self.success_fqn = config.success_fqn
            self.failure_fqn = config.failure_fqn
            self.retry_fqn = config.retry_fqn

        if self.config.log_to_file:
            mc.create_dir(self.config.log_file_directory)
            failure = open(self.failure_fqn, 'w')
            failure.close()
            retry = open(self.retry_fqn, 'w')
            retry.close()
            success = open(self.success_fqn, 'w')
            success.close()
        self.success_count = 0
        self.complete_record_count = 0

    @property
    def complete_record_count(self):
        """:return integer indicating how many inputs (files or observations,
         depending on the configuration) have been processed."""
        return self._complete_record_count

    @complete_record_count.setter
    def complete_record_count(self, value):
        self._complete_record_count = value

    def choose(self, storage_name, command_name, meta_visitors, data_visitors):
        """The logic that decides which descendants of CaomExecute to
        instantiate. This is based on the content of the config.yml file
        for an application.
        :storage_name StorageName extension that handles the naming rules for
            a file in ad.
        :command_name Extension of fits2caom2 (or fits2caom2) that is executed
            for blueprint handling.
        :meta_visitors List of methods that implement the
            visit(observation, **kwargs) signature that require metadata
            access.
        :data_visitors List of methods that implement the
            visit(observation, **kwargs) signature that require data access."""
        executors = []
        if storage_name.is_valid():
            subject, cred_param = self._define_subject()
            cadc_data_client = CadcDataClient(subject)
            caom_repo_client = CAOM2RepoClient(
                subject, self.config.logging_level, self.config.resource_id)
            for task_type in self.task_types:
                self.logger.debug(task_type)
                if task_type == mc.TaskType.SCRAPE:
                    if self.config.use_local_files:
                        executors.append(
                            Collection2CaomScrape(self.config, storage_name,
                                                  command_name))
                    else:
                        raise mc.CadcException(
                            'use_local_files must be True with '
                            'Task Type "SCRAPE"')
                elif task_type == mc.TaskType.STORE:
                    if self.config.use_local_files:
                        executors.append(
                            Collection2CaomStoreClient(
                                self.config, storage_name, command_name,
                                cred_param, cadc_data_client,
                                caom_repo_client))
                    else:
                        raise mc.CadcException(
                            'use_local_files must be True with '
                            'Task Type "STORE"')
                elif task_type == mc.TaskType.INGEST:
                    observation = CaomExecute.repo_cmd_get_client(
                        caom_repo_client, self.config.collection,
                        storage_name.obs_id)
                    if observation is None:
                        if self.config.use_local_files:
                            executors.append(
                                Collection2CaomLocalMetaCreateClient(
                                    self.config, storage_name, command_name,
                                    cred_param, cadc_data_client,
                                    caom_repo_client, meta_visitors))
                        else:
                            executors.append(Collection2CaomMetaCreateClient(
                                self.config, storage_name, command_name,
                                cred_param, cadc_data_client, caom_repo_client,
                                meta_visitors))
                    else:
                        if self.config.use_local_files:
                            executors.append(
                                Collection2CaomLocalMetaUpdateClient(
                                    self.config, storage_name, command_name,
                                    cred_param, cadc_data_client,
                                    caom_repo_client, observation,
                                    meta_visitors))
                        else:
                            executors.append(Collection2CaomMetaUpdateClient(
                                self.config, storage_name, command_name,
                                cred_param, cadc_data_client, caom_repo_client,
                                observation, meta_visitors))
                elif task_type == mc.TaskType.MODIFY:
                    if self.config.use_local_files:
                        if (executors is not None and len(executors) > 0 and
                                isinstance(
                                    executors[0], Collection2CaomScrape)):
                            executors.append(
                                Collection2CaomDataScrape(self.config,
                                                          storage_name,
                                                          command_name,
                                                          data_visitors))
                        else:
                            executors.append(
                                Collection2CaomLocalDataClient(
                                    self.config, storage_name, command_name,
                                    cred_param, cadc_data_client,
                                    caom_repo_client, data_visitors))
                    else:
                        executors.append(Collection2CaomDataClient(
                            self.config, storage_name,
                            command_name, cred_param,
                            cadc_data_client, caom_repo_client, data_visitors,
                            mc.TaskType.MODIFY))
                elif task_type == mc.TaskType.VISIT:
                    executors.append(Collection2CaomClientVisit(
                        self.config, storage_name, cred_param,
                        cadc_data_client, caom_repo_client, meta_visitors))
                else:
                    raise mc.CadcException(
                        'Do not understand task type {}'.format(task_type))
            if (self.config.use_local_files and
                    mc.TaskType.SCRAPE not in self.task_types):
                executors.append(
                    Collection2CaomCompareChecksumClient(
                        self.config, storage_name, command_name,
                        cred_param, cadc_data_client, caom_repo_client))
        else:
            logging.error('{} failed naming validation check.'.format(
                storage_name.obs_id))
            self.capture_failure(storage_name.obs_id,
                                 storage_name.file_name,
                                 'Invalid observation ID')
        return executors

    def capture_failure(self, obs_id, file_name, e):
        """Log an error message to the failure file.
        :obs_id observation ID being processed
        :file_name file name being processed
        :e Exception to log"""
        if self.config.log_to_file:
            failure = open(self.failure_fqn, 'a')
            try:
                min_error = self._minimize_error_message(e)
                failure.write(
                    '{} {} {} {}\n'.format(datetime.now(), obs_id, file_name,
                                           min_error))
            finally:
                failure.close()

            retry = open(self.retry_fqn, 'a')
            try:
                if self.config.features.use_file_names:
                    retry.write('{}\n'.format(file_name))
                else:
                    retry.write('{}\n'.format(obs_id))
            finally:
                retry.close()

    def capture_success(self, obs_id, file_name):
        """Capture, with a timestamp, the successful observations/file names
        that have been processed.
        :obs_id observation ID being processed
        :file_name file name being processed"""
        self.success_count += 1
        if self.config.log_to_file:
            success = open(self.success_fqn, 'a')
            try:
                success.write(
                    '{} {} {}\n'.format(datetime.now(), obs_id, file_name))
                logging.info('Progress - processed {} of {} records.'.format(
                    self.success_count, self.complete_record_count))
            finally:
                success.close()

    def _define_subject(self):
        """Common code to figure out which credentials to use when
        creating an instance of the CadcDataClient and the CAOM2Repo client."""
        if (self.config.proxy_fqn is not None and os.path.exists(
                self.config.proxy_fqn)):
            subject = net.Subject(username=None,
                                  certificate=self.config.proxy_fqn)
            cred_param = '--cert {}'.format(self.config.proxy_fqn)
        elif (self.config.netrc_file is not None and os.path.exists(
                self.config.netrc_file)):
            subject = net.Subject(username=None, certificate=None,
                                  netrc=self.config.netrc_file)
            cred_param = '--netrc {}'.format(self.config.netrc_file)
        else:
            subject = None
            cred_param = ''
            logging.warning(
                'No credentials provided (proxy certificate or netrc file).')
        return subject, cred_param

    @staticmethod
    def _minimize_error_message(e):
        """Turn the long-winded stack trace into something minimal that lends
        itself to awk."""
        if 'Read timed out' in e:
            return 'Read timed out'
        elif 'failed to load external entity' in e:
            return 'caom2repo xml error'
        elif 'Did not retrieve' in e:
            return 'cadc-data get error'
        elif 'NAXES was not set' in e:
            return 'NAXES was not set'
        elif 'Invalid SpatialWCS' in e:
            return 'Invalid SpatialWCS'
        elif 'getProxyCertficate failed' in e:
            return 'getProxyCertificate failed'
        elif 'AlreadyExistsException' in e:
            return 'already exists'
        elif 'Could not find the file' in e:
            return 'cadc-data info failed'
        elif 'md5sum not the same' in e:
            return 'md5sum not the same'
        elif 'Start tag expected' in e:
            return 'XML Syntax Exception'
        elif 'failed to compute metadata' in e:
            return 'Failed to compute metadata'
        else:
            return str(e)


def _set_up_file_logging(config, storage_name):
    """Configure logging to a separate file for each item being processed."""
    log_h = None
    if config.log_to_file:
        log_fqn = os.path.join(config.working_directory,
                               storage_name.log_file)
        if config.log_file_directory is not None:
            log_fqn = os.path.join(config.log_file_directory,
                                   storage_name.log_file)
        log_h = logging.FileHandler(log_fqn)
        formatter = logging.Formatter(
            '%(asctime)s:%(levelname)s:%(name)-12s:%(lineno)d:%(message)s')
        log_h.setLevel(config.logging_level)
        log_h.setFormatter(formatter)
        logging.getLogger().addHandler(log_h)
    return log_h


def _unset_file_logging(config, log_h):
    """Turn off the logging to the separate file for each item being
    processed."""
    if config.log_to_file:
        logging.getLogger().removeHandler(log_h)


def _do_one(config, organizer, storage_name, command_name, meta_visitors,
            data_visitors):
    """Process one item.
    :config mc.Config
    :organizer instance of OrganizeExecutes - for calling the choose method.
    :storage_name instance of StorageName for the collection
    :command_name extension of fits2caom2 for the collection
    :meta_visitors List of metadata visit methods.
    :data_visitors List of data visit methods.
    """
    log_h = _set_up_file_logging(config, storage_name)
    try:
        executors = organizer.choose(storage_name, command_name,
                                     meta_visitors, data_visitors)
        for executor in executors:
            logging.info('Step {} for {}'.format(
                executor.task_type, storage_name.obs_id))
            executor.execute(context=None)
        if len(executors) > 0:
            organizer.capture_success(storage_name.obs_id,
                                      storage_name.file_name)
            return 0
        else:
            logging.info('No executors for {}'.format(
                storage_name.obs_id))
            return -1  # cover the case where file name validation fails
    except Exception as e:
        organizer.capture_failure(storage_name.obs_id,
                                  storage_name.file_name,
                                  e=traceback.format_exc())
        logging.info('Execution failed for {} with {}'.format(
            storage_name.obs_id, e))
        logging.error(traceback.format_exc())
        return -1
    finally:
        _unset_file_logging(config, log_h)


def _run_by_file_list(config, organizer, sname, command_name, proxy,
                      meta_visitors, data_visitors, entry):
    """Process an item from a list of files. Creates the correct instance
    of the StorageName extension, based on Config values.

    :config mc.Config
    :organizer instance of OrganizeExecutes - for calling the choose method.
    :sname which extension of StorageName to instantiate for the collection
    :command_name extension of fits2caom2 for the collection
    :proxy Certificate proxy.
    :meta_visitors List of metadata visit methods.
    :data_visitors List of data visit methods.
    """
    if config.features.use_file_names:
        if config.use_local_files:
            storage_name = sname(file_name=entry, fname_on_disk=entry)
        else:
            storage_name = sname(file_name=entry)
    else:
        if config.use_local_files:
            storage_name = sname(file_name=entry, fname_on_disk=entry)
        else:
            storage_name = sname(obs_id=entry)
    logging.info('Process observation id {}'.format(storage_name.obs_id))
    config.proxy_fqn = proxy
    _do_one(config, organizer, storage_name, command_name,
            meta_visitors, data_visitors)


def _run_todo_file(config, organizer, sname, command_name, proxy,
                   meta_visitors, data_visitors):
    """Process all entries listed in a file.

    :config mc.Config
    :organizer instance of OrganizeExecutes - for calling the choose method.
    :sname which extension of StorageName to instantiate for the collection
    :command_name extension of fits2caom2 for the collection
    :proxy Certificate proxy.
    :meta_visitors List of metadata visit methods.
    :data_visitors List of data visit methods.
    """
    with open(organizer.todo_fqn) as f:
        todo_list_length = sum(1 for _ in f)
    organizer.complete_record_count = todo_list_length
    with open(organizer.todo_fqn) as f:
        for line in f:
            _run_by_file_list(config, organizer, sname, command_name,
                              proxy, meta_visitors, data_visitors,
                              line.strip())


def _run_local_files(config, organizer, sname, command_name, proxy,
                     meta_visitors, data_visitors):
    """Process all entries located in the current working directory.

    :config mc.Config
    :organizer instance of OrganizeExecutes - for calling the choose method.
    :sname which extension of StorageName to instantiate for the collection
    :command_name extension of fits2caom2 for the collection
    :proxy Certificate proxy.
    :meta_visitors List of metadata visit methods.
    :data_visitors List of data visit methods.
    """
    file_list = os.listdir(config.working_directory)
    todo_list = []
    for f in file_list:
        if f.endswith('.fits') or f.endswith('.fits.gz'):
            todo_list.append(f)

    organizer.complete_record_count = len(todo_list)
    for do_file in todo_list:
        _run_by_file_list(config, organizer, sname, command_name,
                          proxy, meta_visitors, data_visitors, do_file)


def _run_by_file(config, storage_name, command_name, proxy, meta_visitors,
                 data_visitors):
    """Process all entries by file name. The file names may be obtained
    from the Config todo entry, from the --todo parameter, or from listing
    files on local disk.

    :param config configures the execution of the application
    :param storage_name which extension of StorageName to instantiate for the
        collection
    :param command_name extension of fits2caom2 for the collection
    :param proxy Certificate proxy.
    :param meta_visitors List of metadata visit methods.
    :param data_visitors List of data visit methods.
    """
    try:
        if config.use_local_files:
            logging.debug(
                'Using files from {}'.format(config.working_directory))
            organize = OrganizeExecutes(config)
            _run_local_files(config, organize, storage_name, command_name,
                             proxy, meta_visitors, data_visitors)
        else:
            parser = ArgumentParser()
            parser.add_argument('--todo',
                                help='Fully-qualified todo file name.')
            args = parser.parse_args()
            if args.todo is not None:
                logging.debug('Using entries from file {}'.format(
                    args.todo))
                organize = OrganizeExecutes(config, args.todo)
            else:
                logging.debug('Using entries from file {}'.format(
                    config.work_file))
                organize = OrganizeExecutes(config)
            _run_todo_file(
                config, organize, storage_name, command_name,
                proxy, meta_visitors, data_visitors)
        logging.info('Done, processed {} of {} correctly.'.format(
            organize.success_count, organize.complete_record_count))
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.error(tb)


def _need_to_retry(config):
    """Evaluate the need to have the pipeline try to re-execute for any
     files/observations that have been logged as failures.

    If log_to_file is not set to True, there is no retry file content
    to retry on..

     :param config does the configuration identify retry information?
     :return True if the configuration and logging information indicate a
        need to attempt to retry the pipeline execution for any entries.
     """
    result = True
    if (config is not None and config.features is not None and
            config.features.expects_retry and config.retry_failures and
            config.log_to_file):
        meta = mc.get_file_meta(config.retry_fqn)
        if meta['size'] == 0:
            logging.info('Checked the retry file {}. There are no logged '
                         'failures.'.format(config.retry_fqn))
            result = False
    else:
        result = False
    return result


def _update_config_for_retry(config, count):
    """
    :param config how the retry information is identified to the pipeline
    :param count the current retry iteration
    :return: the Config instance, updated with locations that reflect
        retry information. The logging should not over-write on a retry
        attempt.
    """
    config.work_file = '{}'.format(config.retry_fqn)
    config.log_file_directory = '{}_{}'.format(
        config.log_file_directory, count)
    # reset the location of the log file names
    config.success_log_file_name = config.success_log_file_name
    config.failure_log_file_name = config.failure_log_file_name
    config.retry_file_name = config.retry_file_name
    return config


def run_by_file(storage_name, command_name, collection, proxy, meta_visitors,
                data_visitors):
    """Process all entries by file name. The file names may be obtained
    from the Config todo entry, from the --todo parameter, or from listing
    files on local disk.

    :param storage_name which extension of StorageName to instantiate for the
        collection
    :param command_name extension of fits2caom2 for the collection
    :param collection string which indicates which collection CAOM instances
        are being created for
    :param proxy Certificate proxy.
    :param meta_visitors List of metadata visit methods.
    :param data_visitors List of data visit methods.
    """
    try:
        config = mc.Config()
        config.get_executors()
        config.collection = collection
        logging.debug(config)
        logger = logging.getLogger()
        logger.setLevel(config.logging_level)
        config.features.supports_composite = False
        _run_by_file(config, storage_name, command_name, proxy, meta_visitors,
                     data_visitors)
        if _need_to_retry(config):
            for count in range(0, config.retry_count):
                _update_config_for_retry(config, count)
                _run_by_file(config, storage_name, command_name, proxy,
                             meta_visitors, data_visitors)
                if not _need_to_retry(config):
                    break
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.error(tb)


def run_single(config, storage_name, command_name, meta_visitors,
               data_visitors):
    """Process a single entry by StorageName detail.

    :param config mc.Config
    :param storage_name instance of StorageName for the collection
    :param command_name extension of fits2caom2 for the collection
    :param meta_visitors List of metadata visit methods.
    :param data_visitors List of data visit methods.
    """
    organizer = OrganizeExecutes(config)
    result = _do_one(config, organizer, storage_name,
                     command_name, meta_visitors, data_visitors)
    sys.exit(result)
