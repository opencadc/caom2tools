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
remaining easily tested. The execute methods do conform to an Airflow API
for operator extension, but please, please, please, do not ever import an
Airflow class here.

Alternatively, if you find yourself adding an if statement here, evaluate
whether the need could be solved with a collection-specific visit method
instead.

Raise the CadcException upon encountering an error. There is no recovery
effort as part of a failure. Log the error and stop the pipeline
execution for an Observation.

The correlations that currently exist:
- use_local_data: True => classes have "Local" in their name
- uses the CadcDataClient, the Caom2RepoClient => classes have "Client"
    in their name
- requires metadata access only => classes have "Meta" in their name
- requires data access => classes have "Data" in their name


On the structure and responsibility of the 'run' methods:
- config file reads are in the 'composable' modules of the individual
  pipelines, so that collection-specific changes will not be surprises
  during execution
- exception handling is also in composable, except for retry loops
-
"""

import distutils.sysconfig
import logging
import os
import re
import traceback

from argparse import ArgumentParser
from astropy.io import fits
from datetime import datetime
from shutil import move

from cadcdata import CadcDataClient
from caom2repo import CAOM2RepoClient
from caom2pipe import manage_composable as mc

__all__ = ['OrganizeExecutes', 'StorageName', 'CaomName', 'OrganizeChooser',
           'run_single', 'run_by_file', 'run_single_from_state',
           'run_from_state', 'run_from_storage_name_instance',
           'run_by_file_storage_name']


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

    def __init__(self, obs_id=None, collection=None, collection_pattern=None,
                 fname_on_disk=None, scheme='ad', archive=None, url=None,
                 mime_encoding=None, mime_type='application/fits'):
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
        :param scheme: string value for the scheme of the file URI.
        :param archive: ad storage unit, defaults to value of
            'collection'
        :param url: if the metadata/data is externally available via http,
            the url for retrieval
        """
        self.obs_id = obs_id
        self.collection = collection
        self.collection_pattern = collection_pattern
        self.scheme = scheme
        self.fname_on_disk = fname_on_disk
        if archive is not None:
            self.archive = archive
        else:
            self.archive = collection
        self._url = url
        self._mime_encoding = mime_encoding
        self._mime_type = mime_type

    @property
    def file_uri(self):
        """The ad URI for the file. Assumes compression."""
        return '{}:{}/{}.gz'.format(
            self.scheme, self.archive, self.file_name)

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
        """The preview URI."""
        return self._get_uri(self.prev)

    @property
    def thumb_uri(self):
        """The thumbnail URI."""
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

    @property
    def url(self):
        """A URL from where a file can be retrieved."""
        return self._url

    @url.setter
    def url(self, value):
        self._url = value

    @property
    def mime_encoding(self):
        """The mime encoding for the file, defaults to None."""
        return self._mime_encoding

    @mime_encoding.setter
    def mime_encoding(self, value):
        self._mime_encoding = value

    @property
    def mime_type(self):
        """The mime type for the file, defaults to application/fits."""
        return self._mime_type

    @mime_type.setter
    def mime_type(self, value):
        self._mime_type = value

    @property
    def lineage(self):
        """The value provided to the --lineage parameter for
        fits2caom2 extensions."""
        return '{}/{}'.format(self.product_id, self.file_uri)

    @property
    def external_urls(self):
        return None

    @fname_on_disk.setter
    def fname_on_disk(self, value):
        self._fname_on_disk = value

    def is_valid(self):
        """:return True if the observation ID conforms to naming rules."""
        pattern = re.compile(self.collection_pattern)
        return pattern.match(self.obs_id)

    def _get_uri(self, fname):
        """The ad URI for a file, without consideration for compression."""
        return '{}:{}/{}'.format(self.scheme, self.archive, fname)

    @staticmethod
    def remove_extensions(name):
        """How to get the file_id from a file_name."""
        return name.replace('.fits', '').replace('.gz', '').replace('.header',
                                                                    '')

    @staticmethod
    def is_preview(entry):
        return '.jpg' in entry


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
        return StorageName.remove_extensions(self.uri.split('/')[1])

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
                 meta_visitors, observable):
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
        :param observable: things that last longer than a pipeline execution
        """
        self.logger = logging.getLogger()
        self.logger.setLevel(config.logging_level)
        formatter = logging.Formatter(
            '%(asctime)s:%(levelname)s:%(name)-12s:%(lineno)d:%(message)s')
        for handler in self.logger.handlers:
            handler.setLevel(config.logging_level)
            handler.setFormatter(formatter)
        self.logging_level_param, self.log_level_as = \
            self._specify_logging_level_param(config.logging_level)
        self.obs_id = storage_name.obs_id
        self.uri = storage_name.file_uri
        self.fname = storage_name.file_name
        self.command_name = command_name
        self.root_dir = config.working_directory
        self.collection = config.collection
        self.archive = config.archive
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
        self.url = storage_name.url
        self.lineage = storage_name.lineage
        self.external_urls_param = self._specify_external_urls_param(
            storage_name.external_urls)
        self.observable = observable
        self.mime_encoding = storage_name.mime_encoding
        self.mime_type = storage_name.mime_type

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
              '--plugin {} --module {} --local {} --lineage {}'.format(
                self.command_name,
                self.logging_level_param, self.cred_param, self.collection,
                self.obs_id, self.model_fqn, plugin, plugin, fqn, self.lineage)
        mc.exec_cmd(cmd, self.log_level_as)

    def _fits2caom2_cmd_client(self):
        """Execute fits2caom with a --cert parameter."""
        plugin = self._find_fits2caom2_plugin()
        # so far, the plugin is also the module :)
        cmd = '{} {} {} --observation {} {} --out {} {} ' \
              '--plugin {} --module {} --lineage {}'.format(
                self.command_name, self.logging_level_param, self.cred_param,
                self.collection, self.obs_id, self.model_fqn,
                self.external_urls_param, plugin, plugin, self.lineage)
        mc.exec_cmd(cmd, self.log_level_as)

    def _fits2caom2_cmd_client_local(self):
        """
        Execute fits2caom with a --cert parameter and a --local parameter.
        """
        plugin = self._find_fits2caom2_plugin()
        # so far, the plugin is also the module :)
        local_fqn = os.path.join(self.working_dir, self.fname)
        cmd = '{} {} {} --observation {} {} --local {} --out {} ' \
              '--plugin {} --module {} --lineage {}'.format(
                self.command_name, self.logging_level_param, self.cred_param,
                self.collection, self.obs_id, local_fqn, self.model_fqn,
                plugin, plugin, self.lineage)
        mc.exec_cmd(cmd, self.log_level_as)

    def _fits2caom2_cmd_in_out_client(self):
        """Execute fits2caom with a --in, a --external_url and a --cert
        parameter."""
        plugin = self._find_fits2caom2_plugin()
        # so far, the plugin is also the module :)
        cmd = '{} {} {} --in {} --out {} {} ' \
              '--plugin {} --module {} --lineage {}'.format(
                self.command_name, self.logging_level_param,
                self.cred_param, self.model_fqn, self.model_fqn,
                self.external_urls_param, plugin, plugin, self.lineage)
        mc.exec_cmd(cmd, self.log_level_as)

    def _fits2caom2_cmd_in_out_local_client(self):
        """Execute fits2caom with a --in, --local and a --cert parameter."""
        plugin = self._find_fits2caom2_plugin()
        # so far, the plugin is also the module :)
        local_fqn = os.path.join(self.working_dir, self.fname)
        cmd = '{} {} {} --in {} --out {} --local {} ' \
              '--plugin {} --module {} --lineage {}'.format(
                self.command_name, self.logging_level_param, self.cred_param,
                self.model_fqn, self.model_fqn, local_fqn, plugin, plugin,
                self.lineage)
        mc.exec_cmd(cmd, self.log_level_as)

    def _repo_cmd_create_client(self, observation):
        """Create an observation instance from the input parameter."""
        mc.repo_create(self.caom_repo_client, observation,
                       self.observable.metrics)

    def _repo_cmd_update_client(self, observation):
        """Update an existing observation instance.  Assumes the obs_id
        values are set correctly."""
        mc.repo_update(self.caom_repo_client, observation,
                       self.observable.metrics)

    def _repo_cmd_read_client(self):
        """Retrieve the existing observation model metadata."""
        return mc.repo_get(self.caom_repo_client, self.collection,
                           self.obs_id, self.observable.metrics)

    def _repo_cmd_delete_client(self, observation):
        """Delete an observation instance based on an input parameter."""
        mc.repo_delete(self.caom_repo_client, observation.collection,
                       observation.observation_id, self.observable.metrics)

    def _cadc_data_put_client(self):
        """Store a collection file."""
        mc.data_put(self.cadc_data_client, self.working_dir,
                    self.fname, self.archive, self.stream, self.mime_type,
                    self.mime_encoding, metrics=self.observable.metrics)

    def _cadc_data_get_client(self):
        """Retrieve an archive file, even if it already exists. This might
        ensure that the latest version of the file is retrieved from
        storage."""
        mc.data_get(self.cadc_data_client, self.working_dir,
                    self.fname, self.archive, self.observable.metrics)

    def _cadc_data_info_file_name_client(self):
        """Execute CadcDataClient.get_file_info with the client instance from
        this class."""
        if self.fname is not None:
            file_info = self.cadc_data_client.get_file_info(
                self.archive, self.fname)
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
            kwargs = {'working_directory': self.working_dir,
                      'cadc_client': self.cadc_data_client,
                      'stream': self.stream,
                      'url': self.url,
                      'observable': self.observable}
            for visitor in self.meta_visitors:
                try:
                    self.logger.debug('Visit for {}'.format(visitor))
                    visitor.visit(observation, **kwargs)
                except Exception as e:
                    raise mc.CadcException(e)

    @staticmethod
    def _specify_external_urls_param(external_urls):
        """Make a list of external urls into a command-line parameter."""
        if external_urls is None:
            result = ''
        else:
            result = '--external_url {}'.format(external_urls)
        return result

    @staticmethod
    def _specify_logging_level_param(logging_level):
        """Make a configured logging level into command-line parameters."""
        lookup = {logging.DEBUG: ('--debug', logging.debug),
                  logging.INFO: ('--verbose', logging.info),
                  logging.WARNING: ('', logging.warning),
                  logging.ERROR: ('--quiet', logging.error)}
        return lookup.get(logging_level, ('', logging.info))

    @staticmethod
    def repo_cmd_get_client(caom_repo_client, collection, observation_id,
                            metrics):
        """Execute the CAOM2Repo 'read' operation using the client instance
        from this class.
        :return an Observation instance, or None, if the observation id
        does not exist."""
        return mc.repo_get(caom_repo_client, collection, observation_id,
                           metrics)


class MetaCreateClient(CaomExecute):
    """Defines the pipeline step for Collection ingestion of metadata into CAOM.
    This requires access to only header information.

    This pipeline step will execute a caom2-repo create."""

    def __init__(self, config, storage_name, command_name,
                 cred_param, cadc_data_client, caom_repo_client,
                 meta_visitors, observable):
        super(MetaCreateClient, self).__init__(
            config, mc.TaskType.INGEST, storage_name, command_name,
            cred_param, cadc_data_client, caom_repo_client, meta_visitors,
            observable)

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

        self.logger.debug('create the observation with xml')
        self._repo_cmd_create_client(observation)

        self.logger.debug('clean up the workspace')
        self._cleanup()

        self.logger.debug('End execute for {}'.format(__name__))


class MetaUpdateClient(CaomExecute):
    """Defines the pipeline step for Collection ingestion of metadata into CAOM.
    This requires access to only header information.

    This pipeline step will execute a caom2-repo update."""

    def __init__(self, config, storage_name, command_name, cred_param,
                 cadc_data_client, caom_repo_client, observation,
                 meta_visitors, observable):
        super(MetaUpdateClient, self).__init__(
            config, mc.TaskType.INGEST, storage_name, command_name, cred_param,
            cadc_data_client, caom_repo_client, meta_visitors, observable)
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

        self.logger.debug('write the observation to disk for next step')
        self._write_model(self.observation)

        self.logger.debug('update the observation with xml')
        self._repo_cmd_update_client(self.observation)

        self.logger.debug('clean up the workspace')
        self._cleanup()

        self.logger.debug('End execute for {}'.format(__name__))


class MetaDeleteCreateClient(CaomExecute):
    """Defines the pipeline step for Collection ingestion of metadata into CAOM.
    This requires access to only header information.

    This pipeline step will execute a caom2-repo delete followed by
    a create, because an update will not support a Simple->Composite
    or Composite->Simple type change for the Observation
    structure."""

    def __init__(self, config, storage_name, command_name,
                 cred_param, cadc_data_client, caom_repo_client,
                 observation, meta_visitors, observable):
        super(MetaDeleteCreateClient, self).__init__(
            config, mc.TaskType.INGEST, storage_name, command_name,
            cred_param, cadc_data_client, caom_repo_client, meta_visitors,
            observable)
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

        self.logger.debug('make a new observation from an existing '
                          'observation')
        self._fits2caom2_cmd_in_out_client()

        self.logger.debug('read the xml into memory from the file')
        self.observation = self._read_model()

        self.logger.debug('the metadata visitors')
        self._visit_meta(self.observation)

        self.logger.debug('the observation exists, delete it')
        self._repo_cmd_delete_client(self.observation)

        self.logger.debug('store the xml')
        self._repo_cmd_create_client(self.observation)

        self.logger.debug('clean up the workspace')
        self._cleanup()

        self.logger.debug('End execute for {}'.format(__name__))


class LocalMetaCreateClient(CaomExecute):
    """Defines the pipeline step for Collection ingestion of metadata into CAOM.
    This requires access to only header information.

    This pipeline step will execute a caom2-repo create."""

    def __init__(self, config, storage_name, command_name, cred_param,
                 cadc_data_client, caom_repo_client, meta_visitors,
                 observable):
        super(LocalMetaCreateClient, self).__init__(
            config, mc.TaskType.INGEST, storage_name, command_name, cred_param,
            cadc_data_client, caom_repo_client, meta_visitors, observable)
        self._define_local_dirs(storage_name)
        self.fname = storage_name.fname_on_disk

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('Find the file name as stored.')
        self._cadc_data_info_file_name_client()

        self.logger.debug('the observation does not exist, so go '
                          'straight to generating the xml, as the main_app '
                          'will retrieve the headers')
        self._fits2caom2_cmd_client_local()

        self.logger.debug('read the xml from disk')
        observation = self._read_model()

        self.logger.debug('the metadata visitors')
        self._visit_meta(observation)

        self.logger.debug('store the xml')
        self._repo_cmd_create_client(observation)

        self.logger.debug('write the updated xml to disk for debugging')
        self._write_model(observation)

        self.logger.debug('End execute for {}'.format(__name__))


class LocalMetaDeleteCreateClient(CaomExecute):
    """Defines the pipeline step for Collection ingestion of metadata into CAOM.
    This requires access to only header information.

    This pipeline step will execute a caom2-repo delete followed by
    a create, because an update will not support a Simple->Composite
    or Composite->Simple type change for the Observation
    structure."""

    def __init__(self, config, storage_name, command_name, cred_param,
                 cadc_data_client, caom_repo_client, observation,
                 meta_visitors, observable):
        super(LocalMetaDeleteCreateClient, self).__init__(
            config, mc.TaskType.INGEST, storage_name, command_name, cred_param,
            cadc_data_client, caom_repo_client, meta_visitors, observable)
        self._define_local_dirs(storage_name)
        self.observation = observation

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('Find the file name as stored.')
        self._cadc_data_info_file_name_client()

        self.logger.debug('write the observation to disk for next step')
        self._write_model(self.observation)

        self.logger.debug('make a new observation from an existing '
                          'observation')
        self._fits2caom2_cmd_in_out_local_client()

        self.logger.debug('read the xml from disk')
        observation = self._read_model()

        self.logger.debug('the metadata visitors')
        self._visit_meta(observation)

        self.logger.debug('the observation exists, delete it')
        self._repo_cmd_delete_client(self.observation)

        self.logger.debug('store the xml')
        self._repo_cmd_create_client(observation)

        self.logger.debug('write the updated xml to disk for debugging')
        self._write_model(observation)

        self.logger.debug('End execute for {}'.format(__name__))


class LocalMetaUpdateClient(CaomExecute):
    """Defines the pipeline step for Collection ingestion of metadata into CAOM.
    This requires access to only header information.

    This pipeline step will execute a caom2-repo update."""

    def __init__(self, config, storage_name, command_name, cred_param,
                 cadc_data_client, caom_repo_client, observation,
                 meta_visitors, observable):
        super(LocalMetaUpdateClient, self).__init__(
            config, mc.TaskType.INGEST, storage_name, command_name, cred_param,
            cadc_data_client, caom_repo_client, meta_visitors, observable)
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
        self._fits2caom2_cmd_in_out_local_client()

        self.logger.debug('read the xml from disk')
        self.observation = self._read_model()

        self.logger.debug('the metadata visitors')
        self._visit_meta(self.observation)

        self.logger.debug('store the xml')
        self._repo_cmd_update_client(self.observation)

        self.logger.debug('write the updated xml to disk for debugging')
        self._write_model(self.observation)

        self.logger.debug('End execute for {}'.format(__name__))


class ClientVisit(CaomExecute):
    """Defines the pipeline step for Collection augmentation by a visitor
     of metadata into CAOM. This assumes a record already exists in CAOM,
     and the update DOES NOT require access to either the header or the data.

    This pipeline step will execute a caom2-repo update."""

    def __init__(self, config, storage_name, cred_param,
                 cadc_data_client, caom_repo_client, meta_visitors,
                 observable):
        super(ClientVisit, self).__init__(
            config, mc.TaskType.VISIT, storage_name, command_name=None,
            cred_param=cred_param, cadc_data_client=cadc_data_client,
            caom_repo_client=caom_repo_client,
            meta_visitors=meta_visitors, observable=observable)
        self.fname = None

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('create the work space, if it does not exist')
        self._create_dir()

        # TODO - run a test to see if this is necessary
        # self.logger.debug('Find the file name as stored.')
        # self._find_file_name_storage_client()

        self.logger.debug('retrieve the existing observation, if it exists')
        observation = self._repo_cmd_read_client()

        self.logger.debug('the metadata visitors')
        self._visit_meta(observation)

        self.logger.debug('store the xml')
        self._repo_cmd_update_client(observation)

        self.logger.debug('write the updated xml to disk for debugging')
        self._write_model(observation)

        self.logger.debug('clean up the workspace')
        self._cleanup()

        self.logger.debug('End execute for {}'.format(__name__))


class DataClient(CaomExecute):
    """Defines the pipeline step for all the operations that
    require access to the file on disk, not just the header data. """

    def __init__(self, config, storage_name, command_name, cred_param,
                 cadc_data_client, caom_repo_client, data_visitors,
                 task_type, observable):
        super(DataClient, self).__init__(
            config, task_type, storage_name, command_name, cred_param,
            cadc_data_client, caom_repo_client, meta_visitors=None,
            observable=observable)
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
                  'cadc_client': self.cadc_data_client,
                  'stream': self.stream,
                  'observable': self.observable}
        for visitor in self.data_visitors:
            try:
                self.logger.debug('Visit for {}'.format(visitor))
                visitor.visit(observation, **kwargs)
            except Exception as e:
                raise mc.CadcException(e)


class LocalDataClient(DataClient):
    """Defines the pipeline step for all the operations that
    require access to the file on disk. This class assumes it has access to
    the files on disk - i.e. there is not need to retrieve the files from
    the CADC storage system."""

    def __init__(self, config, storage_name, command_name, cred_param,
                 cadc_data_client, caom_repo_client, data_visitors,
                 observable):
        super(LocalDataClient, self).__init__(
            config, storage_name, command_name, cred_param,
            cadc_data_client=cadc_data_client,
            caom_repo_client=caom_repo_client, data_visitors=data_visitors,
            task_type=mc.TaskType.MODIFY, observable=observable)
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


class PullClient(CaomExecute):
    """Defines the pipeline step for Collection storage of a file that
    is retrieved via http. The file will be temporarily stored on disk,
    because the cadc-data client doesn't support streaming (yet)."""

    def __init__(self, config, storage_name, command_name, cred_param,
                 cadc_data_client, caom_repo_client, observable):
        super(PullClient, self).__init__(
            config, mc.TaskType.PULL, storage_name, command_name, cred_param,
            cadc_data_client, caom_repo_client, meta_visitors=None,
            observable=observable)
        self.stream = config.stream
        self.fname = storage_name.file_name
        self.local_fqn = os.path.join(self.working_dir, self.fname)

    def execute(self, context):
        self.logger.debug('Begin execute for {} Data'.format(__name__))

        self.logger.debug('create the work space, if it does not exist')
        self._create_dir()

        self.logger.debug('get the input file')
        self._transfer_get()

        self.logger.debug(
            'store the input file {} to ad'.format(self.local_fqn))
        self._cadc_data_put_client()

        self.logger.debug('clean up the workspace')
        self._cleanup()

        self.logger.debug('End execute for {}'.format(__name__))

    def _transfer_get(self):
        """Retrieve a file via http to temporary local storage."""
        self.logger.debug('retrieve {} from {}'.format(self.fname, self.url))
        mc.http_get(self.url, self.local_fqn)
        self._transfer_check()
        self.logger.debug('Successfully retrieved {}'.format(self.fname))

    def _transfer_check(self):
        try:
            hdulist = fits.open(self.local_fqn, memmap=True,
                                lazy_load_hdus=False)
            hdulist.verify('warn')
            for h in hdulist:
                h.verify('warn')
            hdulist.close()
        except (fits.VerifyError, OSError) as e:
            raise mc.CadcException(
                f'astropy verify error {self.local_fqn} when reading {e}')
        # a second check that fails for some NEOSSat cases - if this works,
        # the file might have been correctly retrieved
        try:
            # ignore the return value - if the file is corrupted, the getdata
            # fails, which is the only interesting behaviour here
            fits.getdata(self.local_fqn, ext=0)
        except TypeError as e:
            raise mc.CadcException(
                f'astropy getdata error {self.local_fqn} when reading {e}')


class FtpPullClient(PullClient):
    """Defines the pipeline step for Collection storage of a file that
    is retrieved via ftp. The file will be temporarily stored on disk,
    because the cadc-data client doesn't support streaming (yet)."""

    def __init__(self, config, storage_name, command_name, cred_param,
                 cadc_data_client, caom_repo_client, observable):
        super(FtpPullClient, self).__init__(
            config, storage_name, command_name, cred_param,
            cadc_data_client, caom_repo_client, observable)
        self.stream = config.stream
        self.fname = storage_name.file_name
        self.local_fqn = os.path.join(self.working_dir, self.fname)
        self.ftp_host = config.source_host  # TODO TODO TODO
        self.ftp_fqn = storage_name._ftp_fqn

    def _transfer_get(self):
        """Retrieve a file via ftp to temporary local storage."""
        # Right now, using the complete file name as the entry in the
        # todo.txt file, because ftp clients have been beyond me. Later,
        # that might change, but I'm working on that, so deal with
        # the info that I have for now.
        self.logger.debug('retrieve {} from {}'.format(
            self.ftp_fqn, self.ftp_host))
        mc.ftp_get_timeout(self.ftp_host, self.ftp_fqn, self.local_fqn)
        self._transfer_check()
        self.logger.debug('Successfully retrieved {}'.format(self.fname))


class StoreClient(CaomExecute):
    """Defines the pipeline step for Collection storage of a file. This
    requires access to the file on disk."""

    def __init__(self, config, storage_name, command_name, cred_param,
                 cadc_data_client, caom_repo_client, observable):
        super(StoreClient, self).__init__(
            config, mc.TaskType.STORE, storage_name, command_name, cred_param,
            cadc_data_client, caom_repo_client, meta_visitors=None,
            observable=observable)
        # when files are on disk don't worry about a separate directory
        # per observation
        self.working_dir = self.root_dir
        self.stream = config.stream
        self.fname = storage_name.fname_on_disk

    def execute(self, context):
        self.logger.debug('Begin execute for {} Data'.format(__name__))

        self.logger.debug('store the input file {} to ad'.format(self.fname))
        self._cadc_data_put_client()

        self.logger.debug('End execute for {}'.format(__name__))


class Scrape(CaomExecute):
    """Defines the pipeline step for Collection creation of a CAOM model
    observation. The file containing the metadata is located on disk.
    No record is written to a web service."""

    def __init__(self, config, storage_name, command_name, observable):
        super(Scrape, self).__init__(
            config, mc.TaskType.SCRAPE, storage_name, command_name,
            cred_param='', cadc_data_client=None, caom_repo_client=None,
            meta_visitors=None, observable=observable)
        self._define_local_dirs(storage_name)
        self.fname = storage_name.fname_on_disk
        if self.fname is None:
            self.fname = storage_name.file_name

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('generate the xml from the file on disk')
        self._fits2caom2_cmd_local()

        self.logger.debug('End execute for {}'.format(__name__))


class DataScrape(DataClient):
    """Defines the pipeline step for Collection generation and ingestion of
    operations that require access to the file on disk, with no update to the
    service at the end. This class assumes it has access to the files on disk.
    The organization of this class assumes the 'Scrape' task has been done
    previously, so the model instance exists on disk."""

    def __init__(self, config, storage_name, command_name, data_visitors,
                 observable):
        super(DataScrape, self).__init__(
            config, storage_name, command_name, cred_param='',
            cadc_data_client=None, caom_repo_client=None,
            data_visitors=data_visitors, task_type=mc.TaskType.SCRAPE,
            observable=observable)
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


class LocalMetaCreateClientRemoteStorage(CaomExecute):
    """Defines the pipeline step for Collection ingestion of metadata into
    CAOM. This requires access to only header information.

    The file that contains the metadata is available locally, but this file
    is not, nor will it, be stored in CADC.

    This pipeline step will execute a caom2-repo create."""

    def __init__(self, config, storage_name, command_name,
                 cred_param, cadc_data_client, caom_repo_client,
                 meta_visitors, observable):
        super(LocalMetaCreateClientRemoteStorage, self).__init__(
            config, mc.TaskType.REMOTE, storage_name, command_name,
            cred_param, cadc_data_client, caom_repo_client, meta_visitors,
            observable)
        self._define_local_dirs(storage_name)

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('create the work space, if it does not exist')
        self._create_dir()

        self.logger.debug('the observation does not exist, so go '
                          'straight to generating the xml, as the main_app '
                          'will retrieve the headers')
        self._fits2caom2_cmd_client_local()

        self.logger.debug('read the xml into memory from the file')
        observation = self._read_model()

        self.logger.debug('the metadata visitors')
        self._visit_meta(observation)

        self.logger.debug('store the xml')
        self._repo_cmd_create_client(observation)

        self.logger.debug('clean up the workspace')
        self._cleanup()

        self.logger.debug('End execute for {}'.format(__name__))


class LocalMetaUpdateClientRemoteStorage(CaomExecute):
    """Defines the pipeline step for Collection ingestion of metadata into
    CAOM. This requires access to only header information.

    The file that contains the metadata is available locally, but this file
    is not, nor will it, be stored in CADC.

    This pipeline step will execute a caom2-repo update."""

    def __init__(self, config, storage_name, command_name, cred_param,
                 cadc_data_client, caom_repo_client, observation,
                 meta_visitors, observable):
        super(LocalMetaUpdateClientRemoteStorage, self).__init__(
            config, mc.TaskType.INGEST, storage_name, command_name, cred_param,
            cadc_data_client, caom_repo_client, meta_visitors, observable)
        self._define_local_dirs(storage_name)
        self.observation = observation

    def execute(self, context):
        self.logger.debug('Begin execute for {} Meta'.format(__name__))
        self.logger.debug('the steps:')

        self.logger.debug('write the observation to disk for next step')
        self._write_model(self.observation)

        self.logger.debug('generate the xml, as the main_app will retrieve '
                          'the headers')
        self._fits2caom2_cmd_in_out_local_client()

        self.logger.debug('read the xml from disk')
        self.observation = self._read_model()

        self.logger.debug('the metadata visitors')
        self._visit_meta(self.observation)

        self.logger.debug('store the xml')
        self._repo_cmd_update_client(self.observation)

        self.logger.debug('write the updated xml to disk for debugging')
        self._write_model(self.observation)

        self.logger.debug('End execute for {}'.format(__name__))


class OrganizeChooser(object):
    """Extend this class to provide a way to make collection-specific
    complex conditions available within the OrganizeExecute class."""
    def __init__(self):
        pass

    def needs_delete(self, observation):
        return False

    def use_compressed(self):
        return False


class OrganizeExecutes(object):
    """How to turn on/off various task types in a CaomExecute pipeline."""
    def __init__(self, config, chooser=None, todo_file=None):
        self.config = config
        self.chooser = chooser
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
            self.rejected_fqn = os.path.join(
                self.config.log_file_directory,
                '{}_rejected.yml'.format(todo_name))
            config.rejected_fqn = self.rejected_fqn
        else:
            self.todo_fqn = config.work_fqn
            self.success_fqn = config.success_fqn
            self.failure_fqn = config.failure_fqn
            self.retry_fqn = config.retry_fqn
            self.rejected_fqn = config.rejected_fqn

        mc.create_dir(self.config.log_file_directory)
        now_s = datetime.utcnow().timestamp()
        # if log_to_file is False, should leave no logging trace behind
        if self.config.log_to_file:
            for fqn in [self.failure_fqn, self.retry_fqn, self.success_fqn]:
                back_fqn = '{}.{}.txt'.format(fqn.replace('.txt', ''), now_s)
                OrganizeExecutes._init_log_file(fqn, back_fqn)
        self._success_count = 0
        self._complete_record_count = 0
        self.observable = mc.Observable(mc.Rejected(self.rejected_fqn),
                                        mc.Metrics(config))

    @property
    def success_count(self):
        """:return integer indicating how many inputs (files or observations,
         depending on the configuration) have been successfully processed."""
        return self._success_count

    @success_count.setter
    def success_count(self, value):
        self._success_count = value

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
            if mc.TaskType.SCRAPE in self.task_types:
                cred_param = None
                cadc_data_client = None
                caom_repo_client = None
            else:
                subject, cred_param = self._define_subject()
                cadc_data_client = CadcDataClient(subject)
                caom_repo_client = CAOM2RepoClient(
                    subject, self.config.logging_level,
                    self.config.resource_id)
            for task_type in self.task_types:
                self.logger.debug(task_type)
                if task_type == mc.TaskType.SCRAPE:
                    if self.config.use_local_files:
                        executors.append(
                            Scrape(self.config, storage_name, command_name,
                                   self.observable))
                    else:
                        raise mc.CadcException(
                            'use_local_files must be True with '
                            'Task Type "SCRAPE"')
                elif task_type == mc.TaskType.STORE:
                    if self.config.use_local_files:
                        executors.append(
                            StoreClient(
                                self.config, storage_name, command_name,
                                cred_param, cadc_data_client,
                                caom_repo_client, self.observable))
                    else:
                        raise mc.CadcException(
                            'use_local_files must be True with '
                            'Task Type "STORE"')
                elif task_type == mc.TaskType.INGEST:
                    observation = CaomExecute.repo_cmd_get_client(
                        caom_repo_client, self.config.collection,
                        storage_name.obs_id, self.observable.metrics)
                    if observation is None:
                        if self.config.use_local_files:
                            executors.append(
                                LocalMetaCreateClient(
                                    self.config, storage_name, command_name,
                                    cred_param, cadc_data_client,
                                    caom_repo_client, meta_visitors,
                                    self.observable))
                        else:
                            executors.append(MetaCreateClient(
                                self.config, storage_name, command_name,
                                cred_param, cadc_data_client, caom_repo_client,
                                meta_visitors, self.observable))
                    else:
                        if self.config.use_local_files:
                            if (self.chooser is not None and
                                    self.chooser.needs_delete(observation)):
                                executors.append(
                                    LocalMetaDeleteCreateClient(
                                        self.config, storage_name,
                                        command_name,
                                        cred_param, cadc_data_client,
                                        caom_repo_client, observation,
                                        meta_visitors, self.observable))
                            else:
                                executors.append(
                                    LocalMetaUpdateClient(
                                        self.config, storage_name,
                                        command_name,
                                        cred_param, cadc_data_client,
                                        caom_repo_client, observation,
                                        meta_visitors, self.observable))
                        else:
                            if (self.chooser is not None and
                                    self.chooser.needs_delete(observation)):
                                executors.append(
                                    MetaDeleteCreateClient(
                                        self.config, storage_name,
                                        command_name,
                                        cred_param, cadc_data_client,
                                        caom_repo_client, observation,
                                        meta_visitors, self.observable))
                            else:
                                executors.append(
                                    MetaUpdateClient(
                                        self.config, storage_name,
                                        command_name, cred_param,
                                        cadc_data_client, caom_repo_client,
                                        observation, meta_visitors,
                                        self.observable))
                elif task_type == mc.TaskType.MODIFY:
                    if self.config.use_local_files:
                        if (executors is not None and len(executors) > 0 and
                                isinstance(
                                    executors[0], Scrape)):
                            executors.append(
                                DataScrape(self.config,
                                           storage_name,
                                           command_name,
                                           data_visitors, self.observable))
                        else:
                            executors.append(
                                LocalDataClient(
                                    self.config, storage_name, command_name,
                                    cred_param, cadc_data_client,
                                    caom_repo_client, data_visitors,
                                    self.observable))
                    else:
                        executors.append(DataClient(
                            self.config, storage_name,
                            command_name, cred_param,
                            cadc_data_client, caom_repo_client, data_visitors,
                            mc.TaskType.MODIFY, self.observable))
                elif task_type == mc.TaskType.VISIT:
                    executors.append(ClientVisit(
                        self.config, storage_name, cred_param,
                        cadc_data_client, caom_repo_client, meta_visitors,
                        self.observable))
                elif task_type == mc.TaskType.REMOTE:
                    observation = CaomExecute.repo_cmd_get_client(
                        caom_repo_client, self.config.collection,
                        storage_name.obs_id, self.observable.metrics)
                    if observation is None:
                        if self.config.use_local_files:
                            executors.append(
                                LocalMetaCreateClientRemoteStorage(
                                    self.config, storage_name, command_name,
                                    cred_param, cadc_data_client,
                                    caom_repo_client, meta_visitors,
                                    self.observable))
                        else:
                            raise mc.CadcException(
                                'use_local_files must be True with '
                                'Task Type "REMOTE"')
                    else:
                        if self.config.use_local_files:
                            executors.append(
                                LocalMetaUpdateClientRemoteStorage(
                                    self.config, storage_name, command_name,
                                    cred_param, cadc_data_client,
                                    caom_repo_client, observation,
                                    meta_visitors, self.observable))
                        else:
                            raise mc.CadcException(
                                'use_local_files must be True with '
                                'Task Type "REMOTE"')
                elif task_type == mc.TaskType.PULL:
                    if self.config.features.use_urls:
                        executors.append(
                            PullClient(self.config, storage_name, command_name,
                                       cred_param, cadc_data_client,
                                       caom_repo_client, self.observable))
                    else:
                        executors.append(
                            FtpPullClient(self.config, storage_name,
                                          command_name, cred_param,
                                          cadc_data_client, caom_repo_client,
                                          self.observable))
                elif task_type == mc.TaskType.DEFAULT:
                    pass
                else:
                    raise mc.CadcException(
                        'Do not understand task type {}'.format(task_type))
        else:
            logging.error('{} failed naming validation check.'.format(
                storage_name.obs_id))
            self.capture_failure(storage_name.obs_id,
                                 storage_name.file_name,
                                 'Invalid observation ID')
        return executors

    def capture_failure(self, obs_id, file_name, e):
        """Log an error message to the failure file.

        If the failure is of a known type, also capture it to the rejected
        list. The rejected list will be saved to disk when the execute method
        completes.

        If log_to_file is set to False, don't track the failure log file or
        the retry log file entries, because the application should leave no
        logging trace.

        :obs_id observation ID being processed
        :file_name file name being processed
        :e Exception to log - the entire stack trace, which, if logging
            level is not set to debug, will be lost for debugging purposes.
        """
        if self.config.log_to_file:
            with open(self.failure_fqn, 'a') as failure:
                min_error = self._minimize_error_message(e)
                failure.write(
                    '{} {} {} {}\n'.format(datetime.now(), obs_id, file_name,
                                           min_error))

        # only retry entries that are not permanently marked as rejected
        reason = mc.Rejected.known_failure(e)
        if reason == mc.Rejected.NO_REASON:
            if self.config.log_to_file:
                with open(self.retry_fqn, 'a') as retry:
                    if (self.config.features.use_file_names or
                            self.config.use_local_files):
                        retry.write('{}\n'.format(file_name))
                    else:
                        retry.write('{}\n'.format(obs_id))
        else:
            self.observable.rejected.record(reason, obs_id)

    def capture_success(self, obs_id, file_name, start_time):
        """Capture, with a timestamp, the successful observations/file names
        that have been processed.
        :obs_id observation ID being processed
        :file_name file name being processed
        :start_time seconds since beginning of execution.
        """
        self.success_count += 1
        execution_s = datetime.utcnow().timestamp() - start_time
        # if log_to_file is set to False, leave no logging trace
        if self.config.log_to_file:
            success = open(self.success_fqn, 'a')
            try:
                success.write(
                    '{} {} {} {:.2f}\n'.format(
                        datetime.now(), obs_id, file_name, execution_s))
            finally:
                success.close()
        logging.info('Progress - record {} of {} records processed in '
                     '{:.2f} s.'.format(self.success_count,
                                        self.complete_record_count,
                                        execution_s))

    def is_rejected(self, storage_name):
        """Common code to use the appropriate identifier when checking for
        rejected entries."""
        if self.config.features.use_urls:
            result = self.observable.rejected.is_bad_metadata(
                storage_name.url)
        elif self.config.features.use_file_names:
            result = self.observable.rejected.is_bad_metadata(
                storage_name.file_name)
        else:
            result = self.observable.rejected.is_bad_metadata(
                storage_name.obs_id)
        if result:
            logging.info(
                'Rejected observation {} because of bad metadata'.format(
                    storage_name.obs_id))
        return result

    def _define_subject(self):
        """Common code to figure out which credentials to use when
        creating an instance of the CadcDataClient and the CAOM2Repo client."""
        subject = mc.define_subject(self.config)
        if (self.config.proxy_fqn is not None and os.path.exists(
                self.config.proxy_fqn)):
            logging.debug('Using proxy certificate {} for credentials.'.format(
                self.config.proxy_fqn))
            cred_param = '--cert {}'.format(self.config.proxy_fqn)
        elif (self.config.netrc_file is not None and os.path.exists(
                self.config.netrc_file)):
            logging.debug('Using netrc file {} for credentials.'.format(
                self.config.netrc_file))
            cred_param = '--netrc {}'.format(self.config.netrc_file)
        else:
            cred_param = ''
            logging.warning(
                'No credentials provided (proxy certificate or netrc file).')
            logging.warning(
                'Proxy certificate is {}, netrc file is {}.'.format(
                    self.config.proxy_fqn, self.config.netrc_file))
        return subject, cred_param

    @staticmethod
    def _init_log_file(log_fqn, log_fqn_backup):
        """Keep old versions of the progress files."""
        if os.path.exists(log_fqn) and os.path.getsize(log_fqn) != 0:
            move(log_fqn, log_fqn_backup)
        f_handle = open(log_fqn, 'w')
        f_handle.close()

    @staticmethod
    def _minimize_error_message(e):
        """Turn the long-winded stack trace into something minimal that lends
        itself to awk."""
        if e is None:
            return 'None'
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
        elif 'failed to compute metadata' in e:
            return 'Failed to compute metadata'
        elif 'reset by peer' in e:
            return 'Connection reset by peer'
        elif 'ConnectTimeoutError' in e:
            return 'Connection to host timed out'
        elif 'FileNotFoundError' in e:
            return 'No such file or directory'
        elif 'Must set a value of' in e:
            return 'Value Error'
        elif 'This does not look like a FITS file' in e:
            return 'Not a FITS file'
        elif 'invalid Polygon: segment intersect' in e:
            return 'Segment intersect in polygon'
        elif 'Could not read observation record' in e:
            return 'Observation not found'
        elif 'Broken pipe' in e:
            return 'Broken pipe'
        else:
            return str(e)


def _set_up_file_logging(config, storage_name):
    """Configure logging to a separate file for each entry being processed."""
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
    """Turn off the logging to the separate file for each entry being
    processed."""
    if config.log_to_file:
        logging.getLogger().removeHandler(log_h)


def _finish_run(organizer, config):
    """Code common to the end of all _run_<something> methods."""
    mc.create_dir(config.log_file_directory)
    organizer.observable.rejected.persist_state()
    organizer.observable.metrics.capture()


def _do_one(config, organizer, storage_name, command_name, meta_visitors,
            data_visitors):
    """Process one entry.
    :param config mc.Config
    :param organizer instance of OrganizeExecutes - for calling the choose
        method.
    :param storage_name instance of StorageName for the collection
    :param command_name extension of fits2caom2 for the collection
    :param meta_visitors List of metadata visit methods.
    :param data_visitors List of data visit methods.
    """
    log_h = _set_up_file_logging(config, storage_name)
    start_s = datetime.utcnow().timestamp()
    try:
        if organizer.is_rejected(storage_name):
            organizer.capture_failure(storage_name.obs_id,
                                      storage_name.file_name, e='Rejected')
            # successful rejection of the execution case
            return 0
        executors = organizer.choose(storage_name, command_name,
                                     meta_visitors, data_visitors)
        for executor in executors:
            logging.info('Step {} for {}'.format(
                executor.task_type, storage_name.obs_id))
            executor.execute(context=None)
        if len(executors) > 0:
            organizer.capture_success(storage_name.obs_id,
                                      storage_name.file_name,
                                      start_s)
            return 0
        else:
            logging.info('No executors for {}'.format(
                storage_name))
            return -1  # cover the case where file name validation fails
    except Exception as e:
        organizer.capture_failure(storage_name.obs_id,
                                  storage_name.file_name,
                                  e=traceback.format_exc())
        logging.info('Execution failed for {} with {}'.format(
            storage_name.obs_id, e))
        logging.debug(traceback.format_exc())
        return -1
    finally:
        _unset_file_logging(config, log_h)


def _run_by_file_list(config, organizer, sname, command_name,
                      meta_visitors, data_visitors, entry):
    """Process an entry from a list. The list may contain obs ids,
    file names, or URLs. Creates the expected instance of the StorageName
    extension, based on Config values.

    :param config mc.Config
    :param organizer instance of OrganizeExecutes - for calling the choose
        method.
    :param sname which extension of StorageName to instantiate for the
        collection
    :param command_name extension of fits2caom2 for the collection
    :param meta_visitors List of metadata visit methods.
    :param data_visitors List of data visit methods.
    :param entry what is being processed.
    """
    if config.features.use_file_names:
        if config.use_local_files:
            storage_name = sname(file_name=entry, fname_on_disk=entry)
        else:
            storage_name = sname(file_name=entry)
    elif config.features.use_urls:
        if config.use_local_files:
            raise mc.CadcException('Incompatible configuration. '
                                   'use_local_files does not coexist with'
                                   'feature use_urls.')
        storage_name = sname(url=entry)
    else:
        if config.use_local_files:
            storage_name = sname(file_name=entry, fname_on_disk=entry)
        else:
            storage_name = sname(obs_id=entry)
    return _do_one(config, organizer, storage_name, command_name,
                   meta_visitors, data_visitors)


def _run_todo_file(config, organizer, sname, command_name,
                   meta_visitors, data_visitors):
    """Process all entries listed in a file.

    :param config mc.Config
    :param organizer instance of OrganizeExecutes - for calling the choose
        method.
    :param sname which extension of StorageName to instantiate for the
        collection
    :param command_name extension of fits2caom2 for the collection
    :param meta_visitors List of metadata visit methods.
    :param data_visitors List of data visit methods.
    """
    with open(organizer.todo_fqn) as f:
        todo_list_length = sum(1 for _ in f)
    organizer.complete_record_count = todo_list_length
    result = 0
    with open(organizer.todo_fqn) as f:
        for line in f:
            try:
                entry = line.strip()
                result |= _run_by_file_list(config, organizer, sname,
                                            command_name, meta_visitors,
                                            data_visitors, entry)
            except Exception as e:
                organizer.capture_failure(
                    entry, entry, e=traceback.format_exc())
                logging.info('Execution failed for {} with {}'.format(
                    entry, e))
                logging.debug(traceback.format_exc())
                # then keep processing the rest of the lines in the file
                result = -1
    _finish_run(organizer, config)
    return result


def _run_local_files(config, organizer, sname, command_name,
                     meta_visitors, data_visitors, chooser):
    """Process all entries located in the current working directory.

    :param config mc.Config
    :param organizer instance of OrganizeExecutes - for calling the choose
        method.
    :param sname which extension of StorageName to instantiate for the
        collection
    :param command_name extension of fits2caom2 for the collection
    :param meta_visitors List of metadata visit methods.
    :param data_visitors List of data visit methods.
    :param chooser OrganizeChooser access to collection-specific rules
    """
    file_list = os.listdir(config.working_directory)
    temp_list = []
    for f in file_list:
        if f.endswith('.fits') or f.endswith('.fits.gz'):
            if chooser is not None and chooser.use_compressed():
                if f.endswith('.fits'):
                    temp_list.append('{}.gz'.format(f))
                else:
                    temp_list.append(f)
            else:
                if f.endswith('.fits.gz'):
                    temp_list.append(f.replace('.gz', ''))
                else:
                    temp_list.append(f)
        elif f.endswith('.header'):
            temp_list.append(f)
        elif f.endswith('.gz'):
            temp_list.append(f)

    # make the entries unique
    todo_list = list(set(temp_list))
    organizer.complete_record_count = len(todo_list)
    for do_file in todo_list:
        _run_by_file_list(config, organizer, sname, command_name,
                          meta_visitors, data_visitors, do_file)
    _finish_run(organizer, config)

    if config.need_to_retry():
        for count in range(0, config.retry_count):
            logging.warning(
                'Beginning retry {} in {}'.format(count + 1, os.getcwd()))
            config.update_for_retry(count)

            # make another file list
            temp_list = mc.read_from_file(config.work_fqn)
            todo_list = []
            for ii in temp_list:
                # because the entries in retry aren't compressed names
                todo_list.append('{}.gz'.format(ii.strip()))
            organizer = OrganizeExecutes(config, chooser)
            organizer.complete_record_count = len(todo_list)
            logging.info('Retry {} entries'.format(
                organizer.complete_record_count))
            for redo_file in todo_list:
                try:
                    _run_by_file_list(config, organizer, sname, command_name,
                                      meta_visitors, data_visitors,
                                      redo_file.strip())
                except Exception as e:
                    logging.error(e)
            _finish_run(organizer, config)
            if not config.need_to_retry():
                break
        logging.warning('Done retry attempts.')


def run_by_file(config, storage_name, command_name, meta_visitors,
                data_visitors, chooser=None):
    """Process all entries by file name. The file names may be obtained
    from the Config todo entry, from the --todo parameter, or from listing
    files on local disk.

    :param config configures the execution of the application
    :param storage_name which extension of StorageName to instantiate for the
        collection
    :param command_name extension of fits2caom2 for the collection
    :param meta_visitors List of metadata visit methods.
    :param data_visitors List of data visit methods.
    :param chooser OrganizeChooser instance for detailed CaomExecute
        descendant choices
    """
    logger = logging.getLogger()
    logger.setLevel(config.logging_level)
    logging.debug(config)

    result = 0

    if config.use_local_files:
        logging.debug(
            'Using files from {}'.format(config.working_directory))
        organize = OrganizeExecutes(config, chooser)
        _run_local_files(config, organize, storage_name, command_name,
                         meta_visitors, data_visitors, chooser)
    else:
        parser = ArgumentParser()
        parser.add_argument('--todo',
                            help='Fully-qualified todo file name.')
        args = parser.parse_args()
        if args.todo is not None:
            logging.debug('Using entries from todo file {}'.format(
                args.todo))
            organize = OrganizeExecutes(config, chooser, args.todo)
        else:
            logging.debug('Using entries from file {}'.format(
                config.work_fqn))
            organize = OrganizeExecutes(config, chooser)
        result |= _run_todo_file(config, organize, storage_name, command_name,
                                 meta_visitors, data_visitors)
        if config.need_to_retry():
            for count in range(0, config.retry_count):
                logging.warning('Beginning retry {}'.format(count + 1))
                config.update_for_retry(count)
                organize = OrganizeExecutes(config, chooser)
                try:
                    _run_todo_file(config, organize, storage_name,
                                   command_name, meta_visitors, data_visitors)
                except Exception as e:
                    logging.error(e)
                    result = -1
                if not config.need_to_retry():
                    break
            logging.warning('Done retry attempts.')

    logging.info('Done, processed {} of {} correctly.'.format(
        organize.success_count, organize.complete_record_count))
    return result


def run_by_file_storage_name(config, command_name, meta_visitors,
                             data_visitors, name_builder, chooser=None):
    """Process all entries by a list of work provided by the 'work'
    parameter.

    :param config configures the execution of the application
    :param command_name extension of fits2caom2 for the collection
    :param name_builder builds StorageName instances cleverly
    :param meta_visitors List of metadata visit methods.
    :param data_visitors List of data visit methods.
    :param chooser OrganizeChooser instance for detailed CaomExecute
        descendant choices
    """
    logger = logging.getLogger()
    logger.setLevel(config.logging_level)
    logging.debug(config)

    result = 0

    organizer = OrganizeExecutes(config, chooser)

    with open(organizer.todo_fqn) as f:
        todo_list_length = sum(1 for _ in f)
    organizer.complete_record_count = todo_list_length

    with open(organizer.todo_fqn) as f:
        for line in f:
            try:
                entry = line.strip()
                storage_name = name_builder.build(entry)
                result |= _do_one(config, organizer, storage_name,
                                  command_name, meta_visitors,
                                  data_visitors)
            except Exception as e:
                organizer.capture_failure(
                    entry, entry, e=traceback.format_exc())
                logging.info('Execution failed for {} with {}'.format(
                    entry, e))
                logging.debug(traceback.format_exc())
                # then keep processing the rest of the entries
                result = -1
    _finish_run(organizer, config)

    # if config.need_to_retry():  # TODO
    #     for count in range(0, config.retry_count):
    #         logging.warning('Beginning retry {}'.format(count + 1))
    #         config.update_for_retry(count)
    #         organize = OrganizeExecutes(config, chooser)
    #         try:
    #             _run_todo_file(config, organize, storage_name,
    #                            command_name, meta_visitors,
    #                            data_visitors)
    #         except Exception as e:
    #             logging.error(e)
    #             result = -1
    #         if not config.need_to_retry():
    #             break
    #     logging.warning('Done retry attempts.')

    logging.info('Done, processed {} of {} correctly.'.format(
        organizer.success_count, organizer.complete_record_count))
    return result


def run_by_builder(config, command_name, bookmark_name, meta_visitors,
                   data_visitors, work, name_builder, chooser=None):
    """Process all entries by a list of work provided by the 'work'
    parameter.

    :param config configures the execution of the application
    :param command_name extension of fits2caom2 for the collection
    :param bookmark_name
    :param meta_visitors List of metadata visit methods.
    :param data_visitors List of data visit methods.
    :param work
    :param name_builder builds StorageName instances cleverly
    :param chooser OrganizeChooser instance for detailed CaomExecute
        descendant choices
    """
    return _common_state(config, command_name, meta_visitors, data_visitors,
                         bookmark_name, work, middle=None, sname=name_builder,
                         through=_for_loop_through)


def run_single(config, storage_name, command_name, meta_visitors,
               data_visitors, chooser=None):
    """Process a single entry by StorageName detail.

    :param config mc.Config
    :param storage_name instance of StorageName for the collection
    :param command_name extension of fits2caom2 for the collection
    :param meta_visitors List of metadata visit methods.
    :param data_visitors List of data visit methods.
    :param chooser OrganizeChooser instance for detailed CaomExecute
        descendant choices
    """
    organizer = OrganizeExecutes(config, chooser)
    organizer.complete_record_count = 1
    result = _do_one(config, organizer, storage_name,
                     command_name, meta_visitors, data_visitors)
    logging.debug('run_single result is {}'.format(result))
    _finish_run(organizer, config)
    return result


def run_single_from_state(organizer, config, storage_name, command_name,
                          meta_visitors, data_visitors):
    """Process a single entry by StorageName detail. No sys.exit call.

    :param config mc.Config
    :param storage_name instance of StorageName for the collection
    :param command_name extension of fits2caom2 for the collection
    :param meta_visitors List of metadata visit methods.
    :param data_visitors List of data visit methods.
    :param organizer single organizer instance, maintains log records.
    """
    result = _do_one(config, organizer, storage_name,
                     command_name, meta_visitors, data_visitors)
    logging.info('Result is {} for {}'.format(result, storage_name.file_name))
    _finish_run(organizer, config)
    return result


def run_from_state(config, sname, command_name, meta_visitors,
                   data_visitors, bookmark_name, work):
    """Process a list of entries by queries that are bounded and controlled
    by the content of a state file. No sys.exit call.

    :param config mc.Config
    :param sname StorageName extension for the collection
    :param command_name extension of fits2caom2 for the collection
    :param meta_visitors List of metadata visit methods
    :param data_visitors List of data visit methods
    :param bookmark_name
    :param work Work extension for the collection, query returns a list of
        entries to be processed, init does anything required for work.query
        to make sense.
    """
    return _common_state(config, command_name, meta_visitors, data_visitors,
                         bookmark_name, work, _state_middle, sname,
                         _timebox_through)


def run_from_storage_name_instance(config, command_name, meta_visitors,
                                   data_visitors, bookmark_name, work):
    """Process a list of StorageName instances by queries that are bounded and
    controlled by the content of a state file. No sys.exit call. Avoid an
    interim todo.txt file.

    :param config mc.Config
    :param command_name extension of fits2caom2 for the collection
    :param meta_visitors List of metadata visit methods
    :param data_visitors List of data visit methods
    :param bookmark_name
    :param work Work extension for the collection, query returns a list of
        entries to be processed, init does anything required for work.query
        to make sense.
    """
    return _common_state(config, command_name, meta_visitors, data_visitors,
                         bookmark_name, work, _storage_name_middle, sname=None,
                         through=_timebox_through)


def _storage_name_middle(organizer, entries, config, command_name,
                         meta_visitors, data_visitors, result, sname=None):
    organizer.complete_record_count = len(entries)
    organizer.success_count = 0
    for storage_name in entries:
        try:
            result |= _do_one(config, organizer, storage_name,
                              command_name, meta_visitors,
                              data_visitors)
        except Exception as e:
            organizer.capture_failure(
                storage_name.obs_id, storage_name.file_name,
                e=traceback.format_exc())
            logging.info('Execution failed for {} with {}'.format(
                storage_name.file_name, e))
            logging.debug(traceback.format_exc())
            # then keep processing the rest of the lines in the file
            result |= -1
    _finish_run(organizer, config)

    # TODO - figure this out, do it without using a todo file
    # if config.need_to_retry():
    #     for count in range(0, config.retry_count):
    #         logging.warning('Beginning retry {}'.format(
    # count + 1))
    #         config.update_for_retry(count)
    #         organize = OrganizeExecutes(config, chooser)
    #         try:
    #             _run_todo_file(config, organize, storage_name,
    #                            command_name, meta_visitors,
    # data_visitors)
    #         except Exception as e:
    #             logging.error(e)
    #             result = -1
    #         if not config.need_to_retry():
    #             break
    #     logging.warning('Done retry attempts.')

    return result


def _state_middle(organizer=None, entries=None, config=None,
                  command_name=None, meta_visitors=None, data_visitors=None,
                  result=None, sname=None):
    mc.write_to_file(config.work_fqn, '\n'.join(entries))
    result |= run_by_file(config, sname, command_name,
                          meta_visitors, data_visitors)
    return result


def _common_state(config, command_name, meta_visitors,
                  data_visitors, bookmark_name, work, middle, sname=None,
                  through=None):
    """Process a list of StorageName instances by queries that are bounded and
    controlled by the content of a state file. No sys.exit call. Avoid an
    interim todo.txt file.

    :param config mc.Config
    :param command_name extension of fits2caom2 for the collection
    :param meta_visitors List of metadata visit methods
    :param data_visitors List of data visit methods
    :param bookmark_name
    :param work Work extension for the collection, query returns a list of
        entries to be processed, init does anything required for work.query
        to make sense.
    """
    if not os.path.exists(os.path.dirname(config.progress_fqn)):
        os.makedirs(os.path.dirname(config.progress_fqn))

    logger = logging.getLogger()
    logger.setLevel(config.logging_level)
    logging.debug(config)

    state = mc.State(config.state_fqn)
    start_time = state.get_bookmark(bookmark_name)
    end_time = datetime.fromtimestamp(work.max_ts_s)

    # make sure prev_exec_time is type datetime
    prev_exec_time = mc.increment_time(start_time, 0)
    exec_time = min(
        mc.increment_time(prev_exec_time, config.interval), end_time)

    logging.debug('Starting at {}, ending at {}'.format(start_time, end_time))
    result = 0
    if prev_exec_time == end_time:
        logging.info('Start time is the same as end time {}, stopping.'.format(
            start_time))
        exec_time = prev_exec_time
    else:
        temp_result, exec_time = through(
            config, state, work, middle, command_name, bookmark_name,
            meta_visitors, data_visitors, sname, exec_time, end_time,
            prev_exec_time, start_time)
        result |= temp_result
    state.save_state(bookmark_name, exec_time)
    logging.info(
        'Done {}, saved state is {}'.format(command_name, exec_time))
    return result


def _timebox_through(config, state, work, middle, command_name, bookmark_name,
                     meta_visitors, data_visitors, sname, exec_time, end_time,
                     prev_exec_time, start_time):
    cumulative = 0
    result = 0
    organizer = OrganizeExecutes(config, chooser=None)
    while exec_time <= end_time:
        logging.info(
            'Processing from {} to {}'.format(prev_exec_time, exec_time))
        entries = work.todo(prev_exec_time, exec_time)
        if len(entries) > 0:
            work.initialize()
            logging.info('Processing {} entries.'.format(len(entries)))
            result |= middle(organizer, entries, config, command_name,
                             meta_visitors, data_visitors, result, sname)
        cumulative += len(entries)
        mc.record_progress(
            config, command_name, len(entries), cumulative, start_time)

        state.save_state(bookmark_name, exec_time)

        if exec_time == end_time:
            # the last interval will always have the exec time
            # equal to the end time, which will fail the while check
            # so leave after the last interval has been processed
            #
            # but the while <= check is required so that an interval
            # smaller than exec_time -> end_time will get executed,
            # so don't get rid of the '=' in the while loop
            # comparison, just because this one exists
            break

        prev_exec_time = exec_time
        exec_time = min(
            mc.increment_time(prev_exec_time, config.interval),
            end_time)

    return result, exec_time


def _for_loop_through(config, state, work, middle, command_name, bookmark_name,
                      meta_visitors, data_visitors, name_builder, exec_time,
                      end_time, prev_exec_time, start_time):
    result = 0
    organizer = OrganizeExecutes(config, chooser=None)
    todo_list = work.todo()
    name_builder.todo_list = todo_list
    organizer.complete_record_count = len(todo_list)
    for entry in todo_list.values():
        logging.info(f'Processing {entry}')
        try:
            storage_name = name_builder.build(entry)
            result |= _do_one(config, organizer, storage_name,
                              command_name, meta_visitors,
                              data_visitors)
            # the todo list is sorted by last modified timestamps
            exec_time = storage_name.last_modified_s
        except Exception as e:
            organizer.capture_failure(
                entry, entry, e=traceback.format_exc())
            logging.info(f'Execution failed for {entry} with {e}.')
            logging.debug(traceback.format_exc())
            # then keep processing the rest of the entries
            result = -1

        # interim save - see how this works for now
        if organizer.success_count % 100 == 0:
            if isinstance(exec_time, float):
                exec_time = datetime.fromtimestamp(exec_time)
            logging.debug(f'Saving interim state of {exec_time}')
            state.save_state(bookmark_name, exec_time)
            # note that this will save timestamps past failures, but
            # the failures are recorded in the failure log, so they're
            # known

    _finish_run(organizer, config)

    # if config.need_to_retry():  # TODO
    #     for count in range(0, config.retry_count):
    #         logging.warning('Beginning retry {}'.format(count + 1))
    #         config.update_for_retry(count)
    #         organize = OrganizeExecutes(config, chooser)
    #         try:
    #             _run_todo_file(config, organize, storage_name,
    #                            command_name, meta_visitors,
    #                            data_visitors)
    #         except Exception as e:
    #             logging.error(e)
    #             result = -1
    #         if not config.need_to_retry():
    #             break
    #     logging.warning('Done retry attempts.')

    logging.info(f'Done, processed {organizer.success_count} of '
                 f'{organizer.complete_record_count} correctly.')

    if isinstance(exec_time, float):
        exec_time = datetime.fromtimestamp(exec_time)
    return result, exec_time
