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

import logging
import os
import re
import site
import traceback

from argparse import ArgumentParser
from datetime import datetime

from cadcutils import net, exceptions
from cadcdata import CadcDataClient
from caom2repo import CAOM2RepoClient
from caom2pipe import manage_composable as mc


__all__ = ['OrganizeExecutes', 'StorageName', 'CaomName']


class StorageName(object):
    """Naming rules:
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
        self.obs_id = obs_id
        self.collection = collection
        self.collection_pattern = collection_pattern
        self.fname_on_disk = fname_on_disk

    def get_file_uri(self):
        return 'ad:{}/{}.gz'.format(self.collection, self.get_file_name())

    def get_file_name(self):
        return '{}.fits'.format(self.obs_id)

    def get_compressed_file_name(self):
        return '{}.fits.gz'.format(self.obs_id)

    def get_model_file_name(self):
        return '{}.fits.xml'.format(self.obs_id)

    def get_prev(self):
        return '{}_prev.jpg'.format(self.obs_id)

    def get_thumb(self):
        return '{}_prev_256.jpg'.format(self.obs_id)

    def get_prev_uri(self):
        return self._get_uri(self.get_prev())

    def get_thumb_uri(self):
        return self._get_uri(self.get_thumb())

    def get_obs_id(self):
        return self.obs_id

    def get_log_file(self):
        return '{}.log'.format(self.obs_id)

    def get_product_id(self):
        return self.obs_id

    def get_fname_on_disk(self):
        return self.fname_on_disk

    def is_valid(self):
        pattern = re.compile(self.collection_pattern)
        return pattern.match(self.obs_id)

    def _get_uri(self, fname):
        return 'ad:{}/{}'.format(self.collection, fname)

    @staticmethod
    def remove_extensions(name):
        return name.replace('.fits', '').replace('.gz', '').replace('.header',
                                                                    '')


class CaomName(object):
    """The naming rules for making and decomposing CAOM URIs, all isolated in
    one class. There are probably OMM assumptions built in, but those will
    slowly go away :). """

    def __init__(self, uri):
        self.uri = uri

    def get_file_id(self):
        return self.uri.split('/')[1].split('.')[0]

    def get_file_name(self):
        return self.uri.split('/')[1]

    def get_uncomp_file_name(self):
        return self.get_file_name().replace('.gz', '')


class CaomExecute(object):
    """Abstract class that defines the operations common to all Execute
    classes."""

    def __init__(self, config, task_type, storage_name, command_name,
                 cadc_data_client=None, caom_repo_client=None, cert=None,
                 meta_visitors=None):
        self.logger = logging.getLogger()
        self.logger.setLevel(config.logging_level)
        formatter = logging.Formatter(
            '%(asctime)s:%(levelname)s:%(name)-12s:%(lineno)d:%(message)s')
        for handler in self.logger.handlers:
            handler.setLevel(config.logging_level)
            handler.setFormatter(formatter)
        self.logging_level_param = self._set_logging_level_param(
            config.logging_level)
        self.obs_id = storage_name.get_obs_id()
        self.product_id = storage_name.get_product_id()
        self.uri = storage_name.get_file_uri()
        self.fname = storage_name.get_file_name()
        self.command_name = command_name
        self.root_dir = config.working_directory
        self.collection = config.collection
        self.working_dir = os.path.join(self.root_dir, self.obs_id)
        self.model_fqn = os.path.join(self.working_dir,
                                      storage_name.get_model_file_name())
        self.netrc_fqn = config.netrc_file
        self.resource_id = config.resource_id
        self.task_type = task_type
        self.cadc_data_client = cadc_data_client
        self.caom_repo_client = caom_repo_client
        self.stream = config.stream
        self.cert = cert
        self.meta_visitors = meta_visitors

    def _create_dir(self):
        """Create the working area if it does not already exist."""
        mc.create_dir(self.working_dir)

    def _cleanup(self):
        """Remove a directory and all its contents."""
        if os.path.exists(self.working_dir):
            for ii in os.listdir(self.working_dir):
                os.remove(os.path.join(self.working_dir, ii))
            os.rmdir(self.working_dir)

    def _check_credentials_exist(self):
        """Ensure named credentials exist in this environment."""
        if not os.path.exists(self.netrc_fqn):
            raise mc.CadcException(
                'Credentials do not exist {}.'.format(self.netrc_fqn))

    def _repo_cmd_read(self):
        """Retrieve the existing observaton model metadata."""
        repo_cmd = 'caom2-repo read {} --resource-id {} --netrc {} ' \
                   '{} {} -o {}'.format(self.logging_level_param,
                                        self.resource_id, self.netrc_fqn,
                                        self.collection,
                                        self.obs_id, self.model_fqn)
        mc.exec_cmd(repo_cmd)

    def _repo_cmd_delete(self):
        """Retrieve the existing observaton model metadata."""
        repo_cmd = 'caom2-repo delete {} --resource-id {} --netrc {} ' \
                   '{} {}'.format(self.logging_level_param,
                                  self.resource_id, self.netrc_fqn,
                                  self.collection,
                                  self.obs_id)
        try:
            mc.exec_cmd(repo_cmd)
            if os.path.exists(self.model_fqn):
                os.remove(self.model_fqn)
        except mc.CadcException as e:
            pass
        # TODO - how to tell the difference between 'it doesn't exist', and
        # there's a real failure to pay attention to?
        # raise CadcException('Could not delete the observation in {}'.format(
        #     self.model_fqn))

    def _repo_cmd(self, operation):
        """This repo operation will work for either create or update."""
        repo_cmd = 'caom2-repo {} {} --resource-id {} --netrc ' \
                   '{} {}'.format(operation, self.logging_level_param,
                                  self.resource_id, self.netrc_fqn,
                                  self.model_fqn)
        mc.exec_cmd(repo_cmd)

    def _define_local_dirs(self, storage_name):
        """when files are on disk don't worry about a separate directory
        per observation"""
        self.working_dir = self.root_dir
        self.model_fqn = os.path.join(
            self.working_dir, storage_name.get_model_file_name())

    def _find_file_name_storage(self):
        ad_lookup = self._data_cmd_info()
        if ad_lookup is not None:
            self.fname = ad_lookup
        else:
            raise mc.CadcException(
                'Could not find the file {} in storage as expected.'.format(
                    self.obs_id))

    def _data_cmd_info(self):
        cmd = 'cadc-data info {} --netrc-file {} {} {}'.format(
            self.logging_level_param, self.netrc_fqn, self.collection,
            self.obs_id)
        try:
            result = mc.exec_cmd_info(cmd)
            looked_up_name = None
            if result is not None and 'name: ' in result:
                looked_up_name = result.split('name: ')[1].split()[0]
                self.logger.debug(
                    'Found file name in storage {}'.format(looked_up_name))
            return looked_up_name
        except Exception as e:
            self.logger.debug(e)
            raise mc.CadcException(
                'Failed to execute {} with {}.'.format(cmd, e))

    def _data_cmd_put(self, fname):
        """Store a collection file."""
        data_cmd = 'cadc-data put {} -c --netrc {} {} -s {} {}'.format(
            self.logging_level_param, self.netrc_fqn, self.collection,
            self.stream, fname)
        mc.exec_cmd(data_cmd)

    def _data_cmd_put_client(self, fname):
        """Store a collection file."""
        try:
            self.cadc_data_client.put_file(self.collection, fname, self.stream)
        except Exception as e:
            raise mc.CadcException(
                'Did not store {} with {}'.format(fname, e))

    def _fits2caom2_cmd_local(self):
        fqn = os.path.join(self.working_dir, self.fname)
        plugin = self._find_fits2caom2_plugin()
        cmd = '{} {} --netrc {} --observation {} {} --out {} ' \
              '--plugin {} --local {} --lineage {}/{}'.format(
                self.command_name,
                self.logging_level_param, self.netrc_fqn, self.collection,
                self.obs_id, self.model_fqn, plugin, fqn, self.obs_id,
                self.uri)
        mc.exec_cmd(cmd)

    def _fits2caom2_cmd_local_client(self):
        fqn = os.path.join(self.working_dir, self.fname)
        plugin = self._find_fits2caom2_plugin()
        cmd = '{} {} --cert {} --observation {} {} --out {} ' \
              '--plugin {} --local {} --lineage {}/{}'.format(
                self.command_name, self.logging_level_param, self.cert,
                self.collection,
                self.obs_id, self.model_fqn, plugin, fqn, self.obs_id,
                self.uri)
        mc.exec_cmd(cmd)

    def _compare_checksums(self, fname):
        fqn = os.path.join(self.working_dir, fname)
        mc.compare_checksum(
            self.netrc_fqn, self.collection, fqn)

    def _compare_checksums_client(self, fname):
        fqn = os.path.join(self.working_dir, fname)
        mc.compare_checksum_client(
            self.cadc_data_client, self.collection, fqn)

    def _read_model(self):
        return mc.read_obs_from_file(self.model_fqn)

    def _find_file_name_storage_client(self):
        file_info = self.cadc_data_client.get_file_info(
            self.collection, self.fname)
        self.fname = file_info['name']

    def _repo_cmd_delete_client(self):
        """Retrieve the existing observaton model metadata."""
        try:
            self.caom_repo_client.delete(self.collection, self.obs_id)
        except Exception as e:
            pass  # TODO for now

    def _repo_cmd_create_client(self, observation):
        try:
            self.caom_repo_client.create(observation)
        except Exception as e:
            raise mc.CadcException(
                'Could not create an observation record for {} in {}'.format(
                    self.obs_id, self.resource_id))

    def _repo_cmd_update_client(self, observation):
        try:
            self.caom_repo_client.update(observation)
        except Exception as e:
            raise mc.CadcException(
                'Could not update an observation record for {} in {}'.format(
                    self.obs_id, self.resource_id))

    def _repo_cmd_read_client(self):
        """Retrieve the existing observaton model metadata."""
        try:
            return self.caom_repo_client.read(self.collection, self.obs_id)
        except Exception as e:
            raise mc.CadcException(
                'Could not read observation record for {} in {}'.format(
                    self.obs_id, self.resource_id))

    def _fits2caom2_cmd(self):
        plugin = self._find_fits2caom2_plugin()
        cmd = '{} {} --netrc {} --observation {} {} --out {} ' \
              '--plugin {} --lineage {}/{}'.format(self.command_name,
                                                   self.logging_level_param,
                                                   self.netrc_fqn,
                                                   self.collection,
                                                   self.obs_id, self.model_fqn,
                                                   plugin, self.product_id,
                                                   self.uri)
        mc.exec_cmd(cmd)

    def _fits2caom2_cmd_client(self):
        plugin = self._find_fits2caom2_plugin()
        cmd = '{} {} --cert {} --observation {} {} --out {} ' \
              '--plugin {} --lineage {}/{}'.format(self.command_name,
                                                   self.logging_level_param,
                                                   self.cert, self.collection,
                                                   self.obs_id, self.model_fqn,
                                                   plugin, self.product_id,
                                                   self.uri)
        mc.exec_cmd(cmd)

    def _fits2caom2_cmd_in_out_client(self):
        plugin = self._find_fits2caom2_plugin()
        # TODO add an input parameter
        cmd = '{} {} --cert {} --in {} --out {} ' \
              '--plugin {} --lineage {}/{}'.format(self.command_name,
                                                   self.logging_level_param,
                                                   self.cert, self.model_fqn,
                                                   self.model_fqn,
                                                   plugin, self.product_id,
                                                   self.uri)
        mc.exec_cmd(cmd)

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
        except Exception as e:
            raise mc.CadcException(
                'Did not retrieve {}'.format(fqn))

    def _write_model(self, observation):
        mc.write_obs_to_file(observation, self.model_fqn)

    def _find_fits2caom2_plugin(self):
        packages = site.getsitepackages()
        return os.path.join(packages[0], '{}/{}.py'.format(self.command_name,
                                                           self.command_name))

    def _visit_meta(self, observation):
        kwargs = {'working_directory': self.working_dir}
        for visitor in self.meta_visitors:
            try:
                self.logger.debug('Visit for {}'.format(visitor))
                visitor.visit(observation, **kwargs)
            except Exception as e:
                raise mc.CadcException(e)

    @staticmethod
    def _set_logging_level_param(logging_level):
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
        try:
            observation = caom_repo_client.read(collection, observation_id)
            return observation
        except exceptions.NotFoundException as e:
            return None
        except Exception as e2:
            raise mc.CadcException(
                'Could not retrieve an observation record for {}.'.format(
                    observation_id))


class Collection2CaomMeta(CaomExecute):
    """Defines the pipeline step for OMM ingestion of metadata into CAOM2.
    This requires access to only header information."""

    def __init__(self, config, storage_name, command_name):
        super(Collection2CaomMeta, self).__init__(
            config, mc.TaskType.INGEST, storage_name, command_name)

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')
        self.logger.debug('make sure named credentials exist')
        self._check_credentials_exist()

        self.logger.debug('Find the file name as stored.')
        self._find_file_name_storage()

        self.logger.debug('create the work space, if it does not exist')
        self._create_dir()

        self.logger.debug('remove the existing observation, if it exists, '
                          'because metadata generation is less repeatable '
                          'for updates than for creates.')
        self._repo_cmd_delete()

        self.logger.debug('generate the xml, as the main_app will retrieve '
                          'the headers')
        self._fits2caom2_cmd()

        self.logger.debug('store the xml')
        self._repo_cmd('create')

        self.logger.debug('clean up the workspace')
        self._cleanup()

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomMetaCreateClient(CaomExecute):
    """Defines the pipeline step for Collection ingestion of metadata into CAOM.
    This requires access to only header information.

    This pipeline step will execute a caom2-repo create."""

    def __init__(self, config, storage_name, command_name,
                 cadc_data_client, caom_repo_client, meta_visitors=None):
        super(Collection2CaomMetaCreateClient, self).__init__(
            config, mc.TaskType.AUGMENT, storage_name, command_name,
            cadc_data_client, caom_repo_client, config.proxy, meta_visitors)

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('Find the file name as stored.')
        self._find_file_name_storage_client()

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

    def __init__(self, config, storage_name, command_name,
                 cadc_data_client, caom_repo_client, meta_visitors=None):
        super(Collection2CaomMetaUpdateClient, self).__init__(
            config, mc.TaskType.AUGMENT, storage_name, command_name,
            cadc_data_client, caom_repo_client, config.proxy, meta_visitors)

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('Find the file name as stored.')
        self._find_file_name_storage_client()

        self.logger.debug('create the work space, if it does not exist')
        self._create_dir()

        self.logger.debug('retrieve the existing observation, if it exists')
        observation = CaomExecute.repo_cmd_get_client(
            self.caom_repo_client, self.collection, self.obs_id)

        self.logger.debug('write the observation to disk for next step')
        self._write_model(observation)

        self.logger.debug('generate the xml, as the main_app will retrieve '
                          'the headers')
        self._fits2caom2_cmd_in_out_client()

        self.logger.debug('read the xml from disk')
        observation = self._read_model()

        self.logger.debug('the metadata visitors')
        self._visit_meta(observation)

        self.logger.debug('store the xml')
        self._repo_cmd_update_client(observation)

        self.logger.debug('clean up the workspace')
        self._cleanup()

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomVisit(CaomExecute):
    """Defines the pipeline step for Collection augmentation by a visitor
     of metadata into CAOM. This assumes a record already exists in CAOM,
     and the update DOES NOT require access to either the header or the data.

    This pipeline step will execute a caom2-repo update."""

    def __init__(self, config, storage_name,
                 cadc_data_client, caom_repo_client, meta_visitors=None):
        super(Collection2CaomVisit, self).__init__(
            config, mc.TaskType.VISIT, storage_name, command_name=None,
            cadc_data_client, caom_repo_client, proxy=None,
            meta_visitors=meta_visitors)

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('Find the file name as stored.')
        self._find_file_name_storage_client()

        self.logger.debug('retrieve the existing observation, if it exists')
        observation = CaomExecute.repo_cmd_get_client(
            self.caom_repo_client, self.collection, self.obs_id)

        self.logger.debug('the metadata visitors')
        self._visit_meta(observation)

        self.logger.debug('store the xml')
        self._repo_cmd_update_client(observation)

        self.logger.debug('clean up the workspace')
        self._cleanup()

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomMetaDeleteClient(CaomExecute):
    """Defines the pipeline step for Collection ingestion of metadata into CAOM.
    This requires access to only header information.

    This pipeline step will always delete the record from caom2repo before
    creating a new version."""

    def __init__(self, config, storage_name, command_name,
                 cadc_data_client, caom_repo_client, cert):
        super(Collection2CaomMetaDeleteClient, self).__init__(
            config, mc.TaskType.INGEST, storage_name, command_name,
            cadc_data_client, caom_repo_client, cert)

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('Find the file name as stored.')
        self._find_file_name_storage_client()

        self.logger.debug('create the work space, if it does not exist')
        self._create_dir()

        self.logger.debug('remove the existing observation, if it exists, '
                          'because metadata generation is less repeatable '
                          'for updates than for creates.')
        self._repo_cmd_delete_client()

        self.logger.debug('generate the xml, as the main_app will retrieve '
                          'the headers')
        self._fits2caom2_cmd_client()

        self.logger.debug('read the xml into memory from the file')
        observation = self._read_model()

        self.logger.debug('store the xml')
        self._repo_cmd_create_client(observation)

        self.logger.debug('clean up the workspace')
        self._cleanup()

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomLocalMeta(CaomExecute):
    """Defines the pipeline step for ingestion of metadata into CAOM.
    The file containing the metadata is located on disk."""

    def __init__(self, config, storage_name, command_name):
        super(Collection2CaomLocalMeta, self).__init__(
            config, mc.TaskType.INGEST, storage_name, command_name)
        self._define_local_dirs(storage_name)
        self.fname = storage_name.get_fname_on_disk()

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')
        self.logger.debug('make sure named credentials exist')
        self._check_credentials_exist()

        self.logger.debug('remove the existing observation, if it exists, '
                          'because metadata generation is less repeatable '
                          'for updates than for creates.')
        self._repo_cmd_delete()

        self.logger.debug('generate the xml from the file on disk')
        self._fits2caom2_cmd_local()

        self.logger.debug('store the xml')
        self._repo_cmd('create')

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomLocalMetaClient(CaomExecute):
    """Defines the pipeline step for Collection ingestion of metadata into
    CAOM. The file containing the metadata is located on disk."""

    def __init__(self, config, storage_name, command_name,
                 cadc_data_client, caom_repo_client, cert):
        super(Collection2CaomLocalMetaClient, self).__init__(
            config, mc.TaskType.INGEST, storage_name, command_name,
            cadc_data_client, caom_repo_client, cert)
        self._define_local_dirs(storage_name)
        self.fname = storage_name.get_fname_on_disk()

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('remove the existing observation, if it exists, '
                          'because metadata generation is less repeatable '
                          'for updates than for creates.')
        self._repo_cmd_delete_client()

        self.logger.debug('generate the xml from the file on disk')
        self._fits2caom2_cmd_local_client()

        self.logger.debug('read the xml into memory from the file')
        observation = self._read_model()

        self.logger.debug('store the xml')
        self._repo_cmd_create_client(observation)

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomData(CaomExecute):
    """Defines the pipeline step for all the operations that require
    access to the file on disk, not just the header data. """

    def __init__(self, config, storage_name, command_name, data_visitors=None):
        super(Collection2CaomData, self).__init__(
            config, mc.TaskType.MODIFY, storage_name, command_name)
        self.log_file_directory = config.log_file_directory
        self.prev_fname = storage_name.get_prev()
        self.thumb_fname = storage_name.get_thumb()
        self.data_visitors = data_visitors

    def execute(self, context):
        self.logger.debug('Begin execute for {} Data'.format(__name__))
        self.logger.debug('make sure named credentials exist')
        self._check_credentials_exist()

        self.logger.debug('Find the file name as stored.')
        self._find_file_name_storage()

        self.logger.debug('create the work space, if it does not exist')
        self._create_dir()

        self.logger.debug('get the input file')
        self._cadc_data_get()

        self.logger.debug('get the observation for the existing model')
        self._repo_cmd_read()
        observation = self._read_model()

        self.logger.debug('execute the data visitors')
        self._visit_data(observation)

        self.logger.debug('output the updated xml')
        self._write_model(observation)

        self.logger.debug('store the updated xml')
        self._repo_cmd('update')

        self.logger.debug('clean up the workspace')
        self._cleanup()

        self.logger.debug('End execute for {}'.format(__name__))

    def _cadc_data_get(self):
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

        data_cmd = 'cadc-data get {} -z --netrc {} {} {} -o {}'.format(
            self.logging_level_param, self.netrc_fqn, self.collection,
            self.obs_id, fqn)
        mc.exec_cmd(data_cmd)
        if not os.path.exists(fqn):
            raise mc.CadcException(
                'Failed {}. Did not retrieve {}'.format(data_cmd, fqn))

    def _find_file_name(self):
        self._find_file_name_storage()

    def _visit_data(self, observation):
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


class Collection2CaomDataClient(CaomExecute):
    """Defines the pipeline step for all the operations that
    require access to the file on disk, not just the header data. """

    def __init__(self, config, storage_name, command_name,
                 cadc_data_client, caom_repo_client, data_visitors=None):
        super(Collection2CaomDataClient, self).__init__(
            config, mc.TaskType.MODIFY, storage_name, command_name,
            cadc_data_client, caom_repo_client)
        self.log_file_directory = config.log_file_directory
        self.data_visitors = data_visitors

    def execute(self, context):
        self.logger.debug('Begin execute for {} Data'.format(__name__))

        self.logger.debug('Find the file name as stored.')
        self._find_file_name_storage_client()

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


class Collection2CaomLocalData(Collection2CaomData):
    """Defines the pipeline step for all the operations that require
    access to the file on disk. This class assumes it has access to the
    files on disk."""

    def __init__(self, config, storage_name, command_name, data_visitors=None):
        super(Collection2CaomLocalData, self).__init__(config, storage_name,
                                                       command_name,
                                                       data_visitors)
        self._define_local_dirs(storage_name)
        self.fname = storage_name.get_fname_on_disk()

    def execute(self, context):
        self.logger.debug('Begin execute for {} Data'.format(__name__))
        self.logger.debug('make sure named credentials exist')
        self._check_credentials_exist()

        self.logger.debug('get the observation for the existing model')
        self._repo_cmd_read()
        observation = self._read_model()

        self.logger.debug('execute the data visitors')
        self._visit_data(observation)

        self.logger.debug('output the updated xml')
        self._write_model(observation)

        self.logger.debug('store the updated xml')
        self._repo_cmd('update')

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomLocalDataClient(Collection2CaomDataClient):
    """Defines the pipeline step for  all the operations that
    require access to the file on disk. This class assumes it has access to
    the files on disk."""

    def __init__(self, config, storage_name, command_name,
                 cadc_data_client, caom_repo_client, data_visitors=None):
        super(Collection2CaomLocalDataClient, self).__init__(
            config, storage_name, command_name,
            cadc_data_client=cadc_data_client,
            caom_repo_client=caom_repo_client, data_visitors=data_visitors)
        self._define_local_dirs(storage_name)
        self.fname = storage_name.get_fname_on_disk()

    def execute(self, context):
        self.logger.debug('Begin execute for {} Data'.format(__name__))

        self.logger.debug('get the observation for the existing model')
        observation = self._repo_cmd_read_client()

        self.logger.debug('execute the data visitors')
        self._visit_data(observation)

        self.logger.debug('store the updated xml')
        self._repo_cmd_update_client(observation)

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomStore(CaomExecute):
    """Defines the pipeline step for OMM storage of a file. This requires
    access to the file on disk. It will gzip compress the file."""

    def __init__(self, config, storage_name, command_name):
        super(Collection2CaomStore, self).__init__(
            config, mc.TaskType.STORE, storage_name, command_name)
        # when files are on disk don't worry about a separate directory
        # per observation
        self.working_dir = self.root_dir
        self.storage_host = config.storage_host
        self.stream = config.stream
        self.fname = storage_name.get_file_name()

    def execute(self, context):
        self.logger.debug('Begin execute for {} Data'.format(__name__))
        self.logger.debug('make sure named credentials exist')
        self._check_credentials_exist()

        self.logger.debug('store the input file to ad')
        self._cadc_data_put()

        self.logger.debug('End execute for {}'.format(__name__))

    def _cadc_data_put(self):
        """Store a collection file."""
        data_cmd = 'cadc-data put {} -c --netrc {} {} -s {} {}'.format(
            self.logging_level_param, self.netrc_fqn, self.collection,
            self.stream, self.fname)
        mc.exec_cmd(data_cmd)


class Collection2CaomStoreClient(CaomExecute):
    """Defines the pipeline step for Collection storage of a file. This requires
    access to the file on disk."""

    def __init__(self, config, storage_name, command_name,
                 cadc_data_client, caom_repo_client):
        super(Collection2CaomStoreClient, self).__init__(
            config, mc.TaskType.STORE, storage_name, command_name,
            cadc_data_client, caom_repo_client)
        # when files are on disk don't worry about a separate directory
        # per observation
        self.working_dir = self.root_dir
        self.storage_host = config.storage_host
        self.stream = config.stream
        self.fname = storage_name.get_file_name()

    def execute(self, context):
        self.logger.debug('Begin execute for {} Data'.format(__name__))

        self.logger.debug('store the input file to ad')
        self._data_cmd_put_client(self.fname)

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomScrape(CaomExecute):
    """Defines the pipeline step for Collection creation of a CAOM model
    observation. The file containing the metadata is located on disk.
    No record is written to a web service."""

    def __init__(self, config, storage_name, command_name):
        super(Collection2CaomScrape, self).__init__(
            config, mc.TaskType.SCRAPE, storage_name, command_name)
        self._define_local_dirs(storage_name)
        self.fname = storage_name.get_fname_on_disk()

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('generate the xml from the file on disk')
        self._fits2caom2_cmd_local()

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomDataScrape(Collection2CaomLocalData):
    """Defines the pipeline step for Collection generation and ingestion of
    operations that require access to the file on disk, with no update to the
    service at the end. This class assumes it has access to the files on disk.
    The organization of this class assumes the 'Scrape' task has been done
    previously, so the model instance exists on disk."""

    def __init__(self, config, storage_name, command_name, data_visitors=None):
        super(Collection2CaomDataScrape, self).__init__(
            config, storage_name, command_name, data_visitors)

    def execute(self, context):
        self.logger.debug('Begin execute for {} Data'.format(__name__))

        self.logger.debug('get observation for the existing model from disk')
        observation = self._read_model()

        self.logger.debug('execute the data visitors')
        self._visit_data(observation)

        self.logger.debug('output the updated xml')
        self._write_model(observation)

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomCompareChecksum(CaomExecute):
    """Defines the pipeline step for comparing the checksum of a file on disk
    with the checksum of the supposedly-the-same file stored at CADC.

    This step should be invoked with any other task type that relies on
    files on local disk.
    """

    def __init__(self, config, storage_name, command_name):
        super(Collection2CaomCompareChecksum, self).__init__(
            config, mc.TaskType.CHECKSUM, storage_name, command_name)
        self._define_local_dirs(storage_name)
        self.fname = storage_name.get_fname_on_disk()

    def execute(self, context):
        self.logger.debug('Begin execute for {} '
                          'CompareChecksum'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('generate the xml from the file on disk')
        self._compare_checksums(self.fname)

        self.logger.debug('End execute for {}'.format(__name__))


class Collection2CaomCompareChecksumClient(CaomExecute):
    """Defines the pipeline step for comparing the checksum of a file on disk
    with the checksum of the supposedly-the-same file stored at CADC.

    This step should be invoked with any other task type that relies on
    files on local disk.
    """

    def __init__(self, config, storage_name, command_name, cadc_data_client,
                 caom_repo_client):
        super(Collection2CaomCompareChecksumClient, self).__init__(
            config, mc.TaskType.CHECKSUM, storage_name, command_name,
            cadc_data_client, caom_repo_client)
        self._define_local_dirs(storage_name)
        self.fname = storage_name.get_fname_on_disk()

    def execute(self, context):
        self.logger.debug('Begin execute for {} '
                          'CompareChecksum'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('generate the xml from the file on disk')
        self._compare_checksums_client(self.fname)

        self.logger.debug('End execute for {}'.format(__name__))


class OrganizeExecutes(object):
    """How to turn on/off various steps in the OMM pipeline."""
    def __init__(self, config, todo_file=None):
        self.config = config
        self.task_types = config.task_types
        self.logger = logging.getLogger()
        self.logger.setLevel(config.logging_level)
        if todo_file is not None:
            self.todo_fqn = todo_file
            todo_name = os.path.basename(todo_file).split('.')[0]
            self.success_fqn = os.path.join(self.config.log_file_directory,
                                            '{}_success.log'.format(todo_name))
            self.failure_fqn = os.path.join(self.config.log_file_directory,
                                            '{}_failure.log'.format(todo_name))
            self.retry_fqn = os.path.join(self.config.log_file_directory,
                                          '{}_retry.log'.format(todo_name))
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
        return self._complete_record_count

    @complete_record_count.setter
    def complete_record_count(self, value):
        self._complete_record_count = value

    def choose(self, storage_name, command_name, meta_visitors=None,
               data_visitors=None):
        if self.config.features.use_clients:
            return self._choose_client(storage_name, command_name, meta_visitors,
                                      data_visitors)
        else:
            return self._choose_no_client(storage_name, command_name,
                                         meta_visitors, data_visitors)

    def _choose_no_client(self, storage_name, command_name, meta_visitors=None,
                          data_visitors=None):
        executors = []
        if storage_name.is_valid():
            logging.error('what is the task type {}'.format(self.task_types))
            for task_type in self.task_types:
                if task_type == mc.TaskType.SCRAPE:
                    if self.config.use_local_files:
                        executors.append(
                            Collection2CaomScrape(
                                self.config, storage_name, command_name))
                    else:
                        raise mc.CadcException(
                            'use_local_files must be True with '
                            'Task Type "SCRAPE"')
                elif task_type == mc.TaskType.STORE:
                    if self.config.use_local_files:
                        executors.append(
                            Collection2CaomStore(
                                self.config, storage_name, command_name))
                    else:
                        raise mc.CadcException(
                            'use_local_files must be True with '
                            'Task Type "STORE"')
                elif task_type == mc.TaskType.INGEST:
                    if self.config.use_local_files:
                        executors.append(
                            Collection2CaomLocalMeta(
                                self.config, storage_name, command_name))
                    else:
                        executors.append(
                            Collection2CaomMeta(self.config, storage_name,
                                                command_name))
                elif task_type == mc.TaskType.MODIFY:
                    if self.config.use_local_files:
                        if isinstance(executors[0], Collection2CaomScrape):
                            executors.append(
                                Collection2CaomDataScrape(self.config,
                                                          storage_name,
                                                          command_name,
                                                          data_visitors))
                        else:
                            executors.append(
                                Collection2CaomLocalData(self.config,
                                                         storage_name,
                                                         command_name,
                                                         data_visitors))
                    else:
                        executors.append(
                            Collection2CaomData(self.config, storage_name,
                                                command_name, data_visitors))
                else:
                    raise mc.CadcException(
                        'Do not understand task type {}'.format(task_type))
            if (self.config.use_local_files and
                    mc.TaskType.SCRAPE not in self.task_types):
                executors.append(
                    Collection2CaomCompareChecksum(
                        self.config, storage_name, command_name))
        else:
            logging.error('{} failed naming validation check.'.format(
                storage_name.get_obs_id()))
            self.capture_failure(storage_name.get_obs_id(),
                                 storage_name.get_file_name(),
                                 'Invalid observation ID')
        return executors

    def _choose_client(self, storage_name, command_name, meta_visitors=None,
                       data_visitors=None):
        executors = []
        if storage_name.is_valid():
            subject = net.Subject(username=None, certificate=self.config.proxy)
            cadc_data_client = None
            caom_repo_client = None
            # cadc_data_client = CadcDataClient(subject)
            # caom_repo_client = CAOM2RepoClient(
            #     subject, self.config.logging_level, self.config.resource_id)
            for task_type in self.task_types:
                self.logger.debug(task_type)
                if task_type == mc.TaskType.SCRAPE:
                    if self.config.use_local_files:
                        executors.append(
                            Collection2CaomScrape(
                                self.config, storage_name, command_name))
                    else:
                        raise mc.CadcException(
                            'use_local_files must be True with '
                            'Task Type "SCRAPE"')
                elif task_type == mc.TaskType.STORE:
                    if self.config.use_local_files:
                        executors.append(
                            Collection2CaomStoreClient(
                                self.config, storage_name, command_name,
                                cadc_data_client, caom_repo_client))
                    else:
                        raise mc.CadcException(
                            'use_local_files must be True with '
                            'Task Type "STORE"')
                elif task_type == mc.TaskType.INGEST:
                    if self.config.use_local_files:
                        executors.append(
                            Collection2CaomLocalMetaClient(
                                self.config, storage_name, command_name,
                                cadc_data_client, caom_repo_client,
                                self.config.proxy))
                    else:
                        executors.append(Collection2CaomMetaDeleteClient(
                            self.config, storage_name, command_name,
                            cadc_data_client, caom_repo_client,
                            self.config.proxy))
                elif task_type == mc.TaskType.MODIFY:
                    if self.config.use_local_files:
                        if isinstance(executors[0], Collection2CaomScrape):
                            executors.append(
                                Collection2CaomDataScrape(self.config,
                                                          storage_name,
                                                          command_name,
                                                          data_visitors))
                        else:
                            executors.append(
                                Collection2CaomLocalDataClient(
                                    self.config, storage_name, command_name,
                                    cadc_data_client, caom_repo_client,
                                    data_visitors))
                    else:
                        executors.append(Collection2CaomDataClient(
                            self.config, storage_name, command_name,
                            cadc_data_client, caom_repo_client, data_visitors))
                elif task_type == mc.TaskType.AUGMENT:
                    observation = CaomExecute.repo_cmd_get_client(
                        caom_repo_client, self.config.collection,
                        storage_name.get_obs_id())
                    if observation is None:
                        executors.append(Collection2CaomMetaCreateClient(
                            self.config, storage_name, command_name,
                            cadc_data_client,
                            caom_repo_client))
                    else:
                        executors.append(Collection2CaomMetaUpdateClient(
                            self.config, storage_name, command_name,
                            cadc_data_client,
                            caom_repo_client, meta_visitors))
                elif task_type == mc.TaskType.VISIT:
                    executors.append(Collection2CaomVisit(self.config,
                                                          storage_name,
                                                          caom_repo_client))
                else:
                    raise mc.CadcException(
                        'Do not understand task type {}'.format(task_type))
            if (self.config.use_local_files and
                    mc.TaskType.SCRAPE not in self.task_types):
                executors.append(
                    Collection2CaomCompareChecksumClient(
                        self.config, storage_name, command_name,
                        cadc_data_client, caom_repo_client))
        else:
            logging.error('{} failed naming validation check.'.format(
                storage_name.get_obs_id()))
            self.capture_failure(storage_name.get_obs_id(),
                                 storage_name.get_file_name(),
                                 'Invalid observation ID')
        return executors

    def capture_failure(self, obs_id, file_name, e):
        if self.config.log_to_file:
            failure = open(self.failure_fqn, 'a')
            try:
                min_error = self._minimize_error_message(e)
                failure.write(
                    '{} {} {} {}\n'.format(datetime.now(), obs_id, file_name,
                                           min_error))
            finally:
                failure.close()

            if not self.config.use_local_files:
                retry = open(self.retry_fqn, 'a')
                try:
                    retry.write('{}\n'.format(obs_id))
                finally:
                    retry.close()

    def capture_success(self, obs_id, file_name):
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

    @staticmethod
    def _minimize_error_message(e):
        """Turn the long-winded stack trace into something minimal that lends
        itself to awk."""
        if 'IllegalPolygonException' in e:
            return 'IllegalPolygonException'
        elif 'Read timed out' in e:
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
        else:
            return e


def _set_up_file_logging(config, storage_name):
    log_h = None
    if config.log_to_file:
        log_fqn = os.path.join(config.working_directory,
                               storage_name.get_log_file())
        if config.log_file_directory is not None:
            log_fqn = os.path.join(config.log_file_directory,
                                   storage_name.get_log_file())
        log_h = logging.FileHandler(log_fqn)
        formatter = logging.Formatter(
            '%(asctime)s:%(levelname)s:%(name)-12s:%(lineno)d:%(message)s')
        log_h.setLevel(config.logging_level)
        log_h.setFormatter(formatter)
        logging.getLogger().addHandler(log_h)
    return log_h


def _unset_file_logging(config, log_h):
    if config.log_to_file:
        logging.getLogger().removeHandler(log_h)


def _do_one(config, organizer, storage_name, command_name,
            meta_visitors=None, data_visitors=None):
    log_h = _set_up_file_logging(config, storage_name)
    try:
        executors = organizer.choose(storage_name, command_name,
                                     meta_visitors, data_visitors)
        for executor in executors:
            logging.info('Step {} for {}'.format(
                executor.task_type, storage_name.get_obs_id()))
            executor.execute(context=None)
        if len(executors) > 0:
            organizer.capture_success(storage_name.get_obs_id(),
                                      storage_name.get_file_name())
            return 0
        else:
            logging.info('No executors for {}'.format(
                storage_name.get_obs_id()))
            return -1  # cover the case where file name validation fails
    except Exception as e:
        organizer.capture_failure(storage_name.get_obs_id(),
                                  storage_name.get_file_name(),
                                  e=traceback.format_exc())
        logging.info('Execution failed for {} with {}'.format(
            storage_name.get_obs_id(), e))
        logging.error(traceback.format_exc())
        return -1
    finally:
        _unset_file_logging(config, log_h)


def _run_by_file_list(config, organizer, sname, command_name, proxy, meta_visitors,
                      data_visitors, entry):
    if config.features.use_file_names:
        storage_name = sname(file_id=entry)
    else:
        storage_name = sname(obs_id=entry)
    logging.info('Process observation id {}'.format(
        storage_name.get_obs_id()))
    if config.features.use_clients:
        config.proxy = proxy
    _do_one(config, organizer, storage_name, command_name,
            meta_visitors, data_visitors)


def _run_todo_file(config, organizer, sname, command_name,
                   proxy=None, meta_visitors=None,
                   data_visitors=None):
    with open(organizer.todo_fqn) as f:
        todo_list_length = sum(1 for _ in f)
    organizer.complete_record_count = todo_list_length
    with open(organizer.todo_fqn) as f:
        for line in f:
            _run_by_file_list(config, organizer, sname, command_name,
                              proxy, meta_visitors, data_visitors, line)


def _run_local_files(config, organizer, sname, command_name,
                     proxy=None, meta_visitors=None,
                     data_visitors=None):
    file_list = os.listdir(config.working_directory)
    todo_list = []
    for f in file_list:
        if f.endswith('.fits') or f.endswith('.fits.gz'):
            todo_list.append(f)

    organizer.complete_record_count = len(todo_list)
    for do_file in todo_list:
        _run_by_file_list(config, organizer, sname, command_name,
                          proxy, meta_visitors, data_visitors, do_file)


def run_by_file(storage_name, command_name, collection, proxy=None,
                meta_visitors=None, data_visitors=None):
    try:
        config = mc.Config()
        config.get_executors()
        config.collection = collection
        logging.debug(config)
        logger = logging.getLogger()
        logger.setLevel(config.logging_level)
        if config.use_local_files:
            organize = OrganizeExecutes(config)
            _run_local_files(config, organize, storage_name, command_name,
                             proxy, meta_visitors, data_visitors)
        else:
            parser = ArgumentParser()
            parser.add_argument('--todo',
                                help='Fully-qualified todo file name.')
            args = parser.parse_args()
            if args.todo is not None:
                organize = OrganizeExecutes(config, args.todo)
            else:
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


def run_single(config, storage_name, command_name,
               meta_visitors=None, data_visitors=None):
    import sys
    organizer = OrganizeExecutes(config)
    result = _do_one(config, organizer, storage_name,
                     command_name, meta_visitors, data_visitors)
    sys.exit(result)
