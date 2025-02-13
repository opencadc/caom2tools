# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2025.                            (c) 2025.
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
#  Revision: 4
#
# ***********************************************************************
#

import logging
import re
import sys
import traceback

from astropy.io import fits
from astropy.time import Time
from collections import defaultdict
from datetime import datetime

import caom2

from caom2utils import data_util
from caom2utils.blueprints import ObsBlueprint, _to_int, _to_int_32, _to_float, _to_str
from caom2utils.wcs_parsers import FitsWcsParser, Hdf5WcsParser, WcsParser


class Caom2Exception(Exception):
    """Exception raised when an attempt to create or update a CAOM2 record fails for some reason."""

    pass


class BlueprintParser:
    """
    Extract CAOM2 metadata from files with no WCS information.
    """

    def __init__(self, obs_blueprint=None, uri=None):
        if obs_blueprint:
            self._blueprint = obs_blueprint
        else:
            self._blueprint = ObsBlueprint()
        self._errors = []
        self.logger = logging.getLogger(__name__)
        self.uri = uri
        self.apply_blueprint()

    @property
    def blueprint(self):
        return self._blueprint

    @blueprint.setter
    def blueprint(self, value):
        self._blueprint = value
        self.apply_blueprint()

    def apply_blueprint(self):
        plan = self.blueprint._plan

        #  first apply the functions
        if self.blueprint._module is not None or self.blueprint._module_instance is not None:
            for key, value in plan.items():
                if ObsBlueprint.is_function(value):
                    if self._blueprint._module_instance is None:
                        plan[key] = self._execute_external(value, key, 0)
                    else:
                        plan[key] = self._execute_external_instance(value, key, 0)

        # apply defaults
        for key, value in plan.items():
            if ObsBlueprint.has_default_value(value):
                # there is a default value set
                if key in plan:
                    plan[key] = value[1]

    def augment_observation(self, observation, artifact_uri, product_id=None):
        """
        Augments a given observation with plane structure only.
        :param observation: existing CAOM2 observation to be augmented.
        :param artifact_uri: the key for finding the artifact to augment
        :param product_id: the key for finding for the plane to augment
        """
        self.logger.debug(f'Begin CAOM2 observation augmentation for URI {artifact_uri}.')
        if observation is None or not isinstance(observation, caom2.Observation):
            raise ValueError(f'Observation type mis-match for {observation}.')

        observation.meta_release = self._get_datetime(
            self._get_from_list('Observation.metaRelease', index=0, current=observation.meta_release)
        )
        observation.meta_read_groups = self._get_from_list(
            'Observation.metaReadGroups', index=0, current=observation.meta_read_groups
        )
        observation.meta_producer = self._get_from_list(
            'Observation.metaProducer', index=0, current=observation.meta_producer
        )

        plane = None
        if not product_id:
            product_id = self._get_from_list('Plane.productID', index=0)
        if product_id is None:
            raise ValueError('product ID required')

        for ii in observation.planes:
            if observation.planes[ii].product_id == product_id:
                plane = observation.planes[product_id]
                break
        if plane is None:
            plane = caom2.Plane(product_id=product_id)
            observation.planes[product_id] = plane
        self.augment_plane(plane, artifact_uri)
        self.logger.debug(f'End CAOM2 observation augmentation for {artifact_uri}.')

    def augment_plane(self, plane, artifact_uri):
        """
        Augments a given plane with artifact structure only.
        :param plane: existing CAOM2 plane to be augmented.
        :param artifact_uri:
        """
        self.logger.debug(f'Begin CAOM2 plane augmentation for {artifact_uri}.')
        if plane is None or not isinstance(plane, caom2.Plane):
            raise ValueError(f'Plane type mis-match for {plane}')

        plane.meta_release = self._get_datetime(
            self._get_from_list('Plane.metaRelease', index=0, current=plane.meta_release)
        )
        plane.data_release = self._get_datetime(
            self._get_from_list('Plane.dataRelease', index=0, current=plane.data_release)
        )
        plane.data_product_type = self._to_data_product_type(
            self._get_from_list('Plane.dataProductType', index=0, current=plane.data_product_type)
        )
        plane.calibration_level = self._to_calibration_level(
            _to_int_32(self._get_from_list('Plane.calibrationLevel', index=0, current=plane.calibration_level))
        )
        plane.meta_producer = self._get_from_list('Plane.metaProducer', index=0, current=plane.meta_producer)
        plane.observable = self._get_observable(current=plane.observable)
        plane.provenance = self._get_provenance(plane.provenance)
        plane.metrics = self._get_metrics(current=plane.metrics)
        plane.quality = self._get_quality(current=plane.quality)

        artifact = None
        for ii in plane.artifacts:
            artifact = plane.artifacts[ii]
            if artifact.uri == artifact_uri:
                break
        if artifact is None or artifact.uri != artifact_uri:
            artifact = caom2.Artifact(
                artifact_uri,
                self._to_product_type(self._get_from_list('Artifact.productType', index=0)),
                self._to_release_type(self._get_from_list('Artifact.releaseType', index=0)),
            )
            plane.artifacts[artifact_uri] = artifact
        self.augment_artifact(artifact)
        self.logger.debug(f'End CAOM2 plane augmentation for {artifact_uri}.')

    def augment_artifact(self, artifact):
        """
        Augments a given CAOM2 artifact with available information
        :param artifact: existing CAOM2 artifact to be augmented
        :param index: int Part name, used in specializing classes
        """
        self.logger.debug(f'Begin CAOM2 artifact augmentation for {self.uri}.')
        if artifact is None or not isinstance(artifact, caom2.Artifact):
            raise ValueError(f'Artifact type mis-match for {artifact}')

        artifact.product_type = self._to_product_type(
            self._get_from_list('Artifact.productType', index=0, current=artifact.product_type)
        )
        artifact.release_type = self._to_release_type(
            self._get_from_list('Artifact.releaseType', index=0, current=artifact.release_type)
        )
        artifact.content_type = self._get_from_list('Artifact.contentType', index=0, current=artifact.content_type)
        artifact.content_length = self._get_from_list(
            'Artifact.contentLength', index=0, current=artifact.content_length
        )
        artifact.content_checksum = _to_checksum_uri(
            self._get_from_list('Artifact.contentChecksum', index=0, current=artifact.content_checksum)
        )
        artifact.content_release = self._get_from_list(
            'Artifact.contentRelease', index=0, current=artifact.content_release
        )
        artifact.content_read_groups = self._get_from_list(
            'Artifact.contentReadGroups', index=0, current=artifact.content_read_groups
        )
        artifact.meta_producer = self._get_from_list('Artifact.metaProducer', index=0, current=artifact.meta_producer)
        self.logger.debug(f'End CAOM2 artifact augmentation for {self.uri}.')

    def _get_from_list(self, lookup, index, current=None):
        value = None
        try:
            keywords = self.blueprint._get(lookup, index)
        except KeyError:
            self.add_error(lookup, sys.exc_info()[1])
            self.logger.debug(f'Could not find {lookup} in configuration.')
            if current:
                self.logger.debug(f'{lookup}: using current value of {current!r}.')
                value = current
            return value
        if (
            keywords is not None
            and not ObsBlueprint.needs_lookup(keywords)
            and not ObsBlueprint.is_function(keywords)
        ):
            value = keywords
        elif self._blueprint.update:
            # The first clause: boolean attributes are used to represent three different values: True, False, and
            # unknown. For boolean attributes _only_ assessed that the risk of setting to None accidentally was
            # better than being unable to set a value of 'unknown'.
            #
            # The second clause: the default value for the current parameter in the method signature is 'None', so
            # do not want to inadvertently assign the default value.
            if isinstance(value, bool) or current is not None:
                value = current

        self.logger.debug(f'{lookup}: value is {value}')
        return value

    def _get_set_from_list(self, lookup, index):
        value = None
        keywords = None
        try:
            keywords = self.blueprint._get(lookup)
        except KeyError:
            self.add_error(lookup, sys.exc_info()[1])
            self.logger.debug(f'Could not find \'{lookup}\' in caom2blueprint ' f'configuration.')

        # if there's something useful as a value in the keywords, extract it
        if keywords is not None and any(keywords):
            if ObsBlueprint.needs_lookup(keywords):
                # if there's a default value use it
                if keywords[1]:
                    value = keywords[1]
                    self.logger.debug(f'{lookup}: assigned default value {value}.')
            elif not ObsBlueprint.is_function(keywords):
                value = keywords
                self.logger.debug(f'{lookup}: assigned value {value}.')
        return value

    def add_error(self, key, message):
        self._errors.append('{} {} {}'.format(datetime.now().strftime('%Y-%m-%dT%H:%M:%S'), key, message))

    def _to_data_product_type(self, value):
        return self._to_enum_type(value, caom2.DataProductType)

    def _to_calibration_level(self, value):
        return self._to_enum_type(value, caom2.CalibrationLevel)

    def _to_product_type(self, value):
        return self._to_enum_type(value, caom2.ProductType)

    def _to_release_type(self, value):
        return self._to_enum_type(value, caom2.ReleaseType)

    def _to_enum_type(self, value, to_enum_type):
        if value is None:
            raise ValueError(f'Must set a value of {to_enum_type.__name__} for ' f'{self.uri}.')
        elif isinstance(value, to_enum_type):
            return value
        else:
            return to_enum_type(value)

    def _execute_external(self, value, key, extension):
        """Execute a function supplied by a user, assign a value to a blueprint entry. The input parameters passed
        to the function are the headers as read in by astropy, or the artifact uri.

        :param value the name of the function to apply.
        :param key:
        :param extension: the current extension name or number.
        """
        # determine which of the possible values for parameter the user is hoping for
        if 'uri' in value:
            parameter = self.uri
        elif 'header' in value and isinstance(self, FitsParser):
            parameter = self._headers[extension]
        elif isinstance(self, FitsParser):
            parameter = {'uri': self.uri, 'header': self._headers[extension]}
        else:
            if hasattr(self, '_file'):
                parameter = {'base': self._file}
            else:
                parameter = {'uri': self.uri, 'header': None}

        result = ''
        execute = None
        try:
            execute = getattr(self.blueprint._module, value.split('(')[0])
        except Exception as e:
            msg = 'Failed to find {}.{} for {}'.format(self.blueprint._module.__name__, value.split('(')[0], key)
            self.logger.error(msg)
            self._errors.append(msg)
            tb = traceback.format_exc()
            self.logger.debug(tb)
            self.logger.error(e)
            return result
        if execute:
            try:
                result = execute(parameter)
                self.logger.debug(f'Key {key} calculated value of {result} using {value} type {type(result)}')
            except Exception as e:
                msg = f'Failed to execute {execute.__name__} for {key} in {self.uri}'
                self.logger.error(msg)
                self.logger.debug(f'Input parameter was {parameter}, value was {value}')
                self._errors.append(msg)
                tb = traceback.format_exc()
                self.logger.debug(tb)
                self.logger.error(e)
        return result

    def _execute_external_instance(self, value, key, extension):
        """Execute a function supplied by a user, assign a value to a blueprint entry. The input parameters passed
        to the function are the headers as read in by astropy, or the artifact uri.

        :param value the name of the function to apply.
        :param key:
        :param extension: the current extension name or number.
        :raise Caom2Exception exception raised when there is a recognizable error in the information being used to
            create a CAOM2 record. A correct and consistent CAOM2 record cannot be created from the input metadata.
            The client should treat the Observation instance under construction as invalid.
        """
        result = ''
        try:
            execute = getattr(self.blueprint._module_instance, value.split('(')[0])
        except Exception as e:
            msg = 'Failed to find {}.{} for {}'.format(
                self.blueprint._module_instance.__class__.__name__, value.split('(')[0], key
            )
            self.logger.error(msg)
            self._errors.append(msg)
            tb = traceback.format_exc()
            self.logger.debug(tb)
            self.logger.error(e)
            return result
        try:
            result = execute(extension)
            self.logger.debug(f'Key {key} calculated value of {result} using {value}')
        except ValueError as e2:
            # DB 23-03-22
            # Anything that you can do to make the CAOM2 record creation fail in this case of bad WCS metadata
            # would be useful. Use ValueError because that happens to be what astropy is throwing for a SkyCoord
            # construction failure.
            raise Caom2Exception(e2)
        except Exception as e:
            msg = 'Failed to execute {} for {} in {}'.format(execute, key, self.uri)
            self.logger.error(msg)
            self.logger.debug('Input value was {}'.format(value))
            self._errors.append(msg)
            tb = traceback.format_exc()
            self.logger.debug(tb)
            self.logger.error(e)
        return result

    def _get_datetime(self, from_value):
        """
        Ensure datetime values are in MJD. Really. Just not yet.
        :param from_value:
        :return:
        """
        if from_value:
            if isinstance(from_value, datetime):
                return from_value
            elif isinstance(from_value, Time):
                return from_value.datetime
            else:
                result = None
                # CFHT 2003/03/29,01:34:54
                # CFHT 2003/03/29
                # DDO 12/02/95
                # TAOSII 2024-01-26T14:52:49Z
                for dt_format in [
                    '%Y-%m-%dT%H:%M:%S',
                    '%Y-%m-%dT%H:%M:%S.%f',
                    '%Y-%m-%d %H:%M:%S.%f',
                    '%Y-%m-%d',
                    '%Y/%m/%d %H:%M:%S',
                    '%Y-%m-%d %H:%M:%S',
                    '%Y/%m/%d,%H:%M:%S',
                    '%Y/%m/%d',
                    '%d/%m/%y',
                    '%d/%m/%y %H:%M:%S',
                    '%d-%m-%Y',
                    '%Y-%m-%dT%H:%M:%SZ',
                ]:
                    try:
                        result = datetime.strptime(from_value, dt_format)
                    except ValueError:
                        pass

                if result is None:
                    self.logger.error('Cannot parse datetime {}'.format(from_value))
                    self.add_error('get_datetime', sys.exc_info()[1])
                return result
        else:
            return None

    def _get_metrics(self, current):
        """
        Create a Metrics instance populated with available content information.
        :return: Metrics
        """
        self.logger.debug('Begin Metrics augmentation.')
        source_number_density = self._get_from_list(
            'Plane.metrics.sourceNumberDensity',
            index=0,
            current=None if current is None else current.source_number_density,
        )
        background = self._get_from_list(
            'Plane.metrics.background', index=0, current=None if current is None else current.background
        )
        background_stddev = self._get_from_list(
            'Plane.metrics.backgroundStddev', index=0, current=None if current is None else current.background_std_dev
        )
        flux_density_limit = self._get_from_list(
            'Plane.metrics.fluxDensityLimit', index=0, current=None if current is None else current.flux_density_limit
        )
        mag_limit = self._get_from_list(
            'Plane.metrics.magLimit', index=0, current=None if current is None else current.mag_limit
        )
        sample_snr = self._get_from_list(
            'Plane.metrics.sampleSNR', index=0, current=None if current is None else current.sample_snr
        )

        metrics = None
        if source_number_density or background or background_stddev or flux_density_limit or mag_limit or sample_snr:
            metrics = caom2.Metrics()
            metrics.source_number_density = source_number_density
            metrics.background = background
            metrics.background_std_dev = background_stddev
            metrics.flux_density_limit = flux_density_limit
            metrics.mag_limit = mag_limit
            metrics.sample_snr = sample_snr
        self.logger.debug('End Metrics augmentation.')
        return metrics

    def _get_observable(self, current):
        """
        Create a Observable instance populated with available content
        information.
        :return: Observable
        """
        self.logger.debug('Begin Observable augmentation.')
        ucd = self._get_from_list('Plane.observable.ucd', index=0, current=None if current is None else current.ucd)
        observable = caom2.Observable(ucd) if ucd else None
        self.logger.debug('End Observable augmentation.')
        return observable

    def _get_provenance(self, current):
        """
        Create a Provenance instance populated with available Content
        information.
        :return: Provenance
        """
        self.logger.debug('Begin Provenance augmentation.')
        name = _to_str(
            self._get_from_list('Plane.provenance.name', index=0, current=None if current is None else current.name)
        )
        p_version = _to_str(
            self._get_from_list(
                'Plane.provenance.version', index=0, current=None if current is None else current.version
            )
        )
        project = _to_str(
            self._get_from_list(
                'Plane.provenance.project', index=0, current=None if current is None else current.project
            )
        )
        producer = _to_str(
            self._get_from_list(
                'Plane.provenance.producer', index=0, current=None if current is None else current.producer
            )
        )
        run_id = _to_str(
            self._get_from_list(
                'Plane.provenance.runID', index=0, current=None if current is None else current.run_id
            )
        )
        reference = _to_str(
            self._get_from_list(
                'Plane.provenance.reference', index=0, current=None if current is None else current.reference
            )
        )
        last_executed = self._get_datetime(
            self._get_from_list(
                'Plane.provenance.lastExecuted', index=0, current=None if current is None else current.last_executed
            )
        )
        keywords = self._get_set_from_list('Plane.provenance.keywords', index=0)
        inputs = self._get_set_from_list('Plane.provenance.inputs', index=0)
        prov = None
        if name:
            prov = caom2.Provenance(name, p_version, project, producer, run_id, reference, last_executed)
            ContentParser._add_keywords(keywords, current, prov)
            if inputs is not None and any(inputs):
                if isinstance(inputs, caom2.TypedSet):
                    for i in inputs:
                        prov.inputs.add(i)
                else:
                    if isinstance(inputs, str):
                        for i in inputs.split():
                            prov.inputs.add(caom2.PlaneURI(str(i)))
                    else:
                        for i in inputs:
                            prov.inputs.add(caom2.PlaneURI(str(i)))
            else:
                if current is not None and len(current.inputs) > 0:
                    # preserve the original value
                    prov.inputs.update(current.inputs)
        self.logger.debug('End Provenance augmentation.')
        return prov

    def _get_quality(self, current):
        """
        Create a Quality instance populated with available content information.
        :return: Quality
        """
        self.logger.debug('Begin Quality augmentation.')
        flag = self._get_from_list('Plane.dataQuality', index=0, current=None if current is None else current.flag)
        quality = caom2.DataQuality(flag) if flag else None
        self.logger.debug('End Quality augmentation.')
        return quality


class ContentParser(BlueprintParser):
    def __init__(self, obs_blueprint=None, uri=None, extension_start_index=0, extension_end_index=None):
        super().__init__(obs_blueprint, uri)
        # for those cases where the extensions of interest are not all the extensions in the original file
        self._extension_start_index = extension_start_index
        self._extension_end_index = extension_end_index if extension_end_index else self._get_num_parts()
        self._wcs_parsers = {}
        self._set_wcs_parsers(obs_blueprint)

    def _get_chunk_naxis(self, chunk, index):
        chunk.naxis = self._get_from_list('Chunk.naxis', index, self.blueprint.get_configed_axes_count())

    def _get_num_parts(self):
        """return the number of Parts to create for a CAOM record
        """
        return len(self._blueprint._extensions) + 1

    def _set_wcs_parsers(self, obs_blueprint):
        self._wcs_parsers[0] = WcsParser(obs_blueprint, extension=self._extension_start_index)

    def augment_artifact(self, artifact):
        """
        Augments a given CAOM2 artifact with available content information
        :param artifact: existing CAOM2 artifact to be augmented
        :param index: int Part name
        """
        super().augment_artifact(artifact)
        self.logger.debug(f'Begin content artifact augmentation for {artifact.uri}')

        if self.blueprint.get_configed_axes_count() == 0:
            raise TypeError(f'No WCS Data. End content artifact augmentation for ' f'{artifact.uri}.')

        for index in range(self._extension_start_index, self._extension_end_index):
            if self.add_parts(artifact, index):
                part = artifact.parts[str(index)]
                part.product_type = self._get_from_list('Part.productType', index)
                part.meta_producer = self._get_from_list(
                    'Part.metaProducer', index=self._extension_start_index, current=part.meta_producer
                )

                # each Part has one Chunk, if it's not an empty part as determined just previously
                if not part.chunks:
                    part.chunks.append(caom2.Chunk())
                chunk = part.chunks[0]
                chunk.meta_producer = self._get_from_list(
                    'Chunk.metaProducer', index=self._extension_start_index, current=chunk.meta_producer
                )

                self._get_chunk_naxis(chunk, index)

                # order by which the blueprint is used to set WCS information:
                # 1 - try to construct the information for an axis from WCS information
                # 2 - if the WCS information is insufficient, try to construct the information from the blueprint
                # 3 - Always try to fill the range metadata from the blueprint.
                if self.blueprint._pos_axes_configed:
                    self._wcs_parsers[index].augment_position(chunk)
                    self._try_position_with_blueprint(chunk, index)

                if self.blueprint._energy_axis_configed:
                    self._wcs_parsers[index].augment_energy(chunk)
                    self._try_energy_with_blueprint(chunk, index)

                if self.blueprint._time_axis_configed:
                    self._wcs_parsers[index].augment_temporal(chunk)
                    self._try_time_with_blueprint(chunk, index)

                if self.blueprint._polarization_axis_configed:
                    self._wcs_parsers[index].augment_polarization(chunk)
                    self._try_polarization_with_blueprint(chunk, index)

                if self.blueprint._obs_axis_configed:
                    self._wcs_parsers[index].augment_observable(chunk)
                    self._try_observable_with_blueprint(chunk, index)

                if self.blueprint._custom_axis_configed:
                    self._wcs_parsers[index].augment_custom(chunk)
                    self._try_custom_with_blueprint(chunk, index)

        self.logger.debug(f'End content artifact augmentation for {artifact.uri}.')

    def augment_observation(self, observation, artifact_uri, product_id=None):
        """
        Augments a given observation with available content information.
        :param observation: existing CAOM2 observation to be augmented.
        :param artifact_uri: the key for finding the artifact to augment
        :param product_id: the key for finding for the plane to augment
        """
        super().augment_observation(observation, artifact_uri, product_id)
        self.logger.debug(f'Begin content observation augmentation for URI {artifact_uri}.')
        members = self._get_members(observation)
        if members:
            if isinstance(members, caom2.TypedSet):
                for m in members:
                    observation.members.add(m)
            else:
                for m in members.split():
                    observation.members.add(caom2.ObservationURI(m))
        observation.algorithm = self._get_algorithm(observation)

        observation.sequence_number = _to_int(self._get_from_list('Observation.sequenceNumber', index=0))
        observation.intent = self._get_from_list(
            'Observation.intent',
            0,
            (caom2.ObservationIntentType.SCIENCE if observation.intent is None else observation.intent),
        )
        observation.type = self._get_from_list('Observation.type', 0, current=observation.type)
        observation.meta_release = self._get_datetime(
            self._get_from_list('Observation.metaRelease', 0, current=observation.meta_release)
        )
        observation.meta_read_groups = self._get_from_list('Observation.metaReadGroups', 0)
        observation.meta_producer = self._get_from_list(
            'Observation.metaProducer', 0, current=observation.meta_producer
        )
        observation.requirements = self._get_requirements(observation.requirements)
        observation.instrument = self._get_instrument(observation.instrument)
        observation.proposal = self._get_proposal(observation.proposal)
        observation.target = self._get_target(observation.target)
        observation.target_position = self._get_target_position(observation.target_position)
        observation.telescope = self._get_telescope(observation.telescope)
        observation.environment = self._get_environment(observation.environment)
        self.logger.debug('End content observation augmentation.')

    def _get_algorithm(self, obs):
        """
        Create an Algorithm instance populated with available content
        information.
        :return: Algorithm
        """
        self.logger.debug('Begin Algorithm augmentation.')
        # TODO DEFAULT VALUE
        name = self._get_from_list('Observation.algorithm.name', index=0, current=obs.algorithm.name)
        if name is not None and name == 'exposure' and isinstance(obs, caom2.DerivedObservation):
            # stop the raising of a ValueError when adding a Plane representing a SimpleObservation to a
            # DerivedObservation under construction. It results in attempting to change Algorithm.name value to
            # 'exposure' otherwise.
            result = obs.algorithm
        else:
            result = caom2.Algorithm(str(name)) if name else None
        self.logger.debug('End Algorithm augmentation.')
        return result

    def _get_energy_transition(self, current):
        """
        Create an EnergyTransition instance populated with available content information.
        :return: EnergyTransition
        """
        self.logger.debug('Begin EnergyTransition augmentation.')
        species = self._get_from_list(
            'Chunk.energy.transition.species', index=0, current=None if current is None else current.species
        )
        transition = self._get_from_list(
            'Chunk.energy.transition.transition', index=0, current=None if current is None else current.transition
        )
        result = None
        if species is not None and transition is not None:
            result = caom2.EnergyTransition(species, transition)
        self.logger.debug('End EnergyTransition augmentation.')
        return result

    def _get_environment(self, current):
        """
        Create an Environment instance populated with available content information.
        :current Environment instance, if one already exists in the Observation
        :return: Environment
        """
        self.logger.debug('Begin Environment augmentation.')
        seeing = self._get_from_list(
            'Observation.environment.seeing', index=0, current=None if current is None else current.seeing
        )
        humidity = _to_float(
            self._get_from_list(
                'Observation.environment.humidity', index=0, current=None if current is None else current.humidity
            )
        )
        elevation = self._get_from_list(
            'Observation.environment.elevation', index=0, current=None if current is None else current.elevation
        )
        tau = self._get_from_list(
            'Observation.environment.tau', index=0, current=None if current is None else current.tau
        )
        wavelength_tau = self._get_from_list(
            'Observation.environment.wavelengthTau',
            index=0,
            current=None if current is None else current.wavelength_tau,
        )
        ambient = _to_float(
            self._get_from_list(
                'Observation.environment.ambientTemp',
                index=0,
                current=None if current is None else current.ambient_temp,
            )
        )
        photometric = self._cast_as_bool(
            self._get_from_list(
                'Observation.environment.photometric',
                index=0,
                current=None if current is None else current.photometric,
            )
        )
        enviro = None
        if seeing or humidity or elevation or tau or wavelength_tau or ambient:
            enviro = caom2.Environment()
            enviro.seeing = seeing
            enviro.humidity = humidity
            enviro.elevation = elevation
            enviro.tau = tau
            enviro.wavelength_tau = wavelength_tau
            enviro.ambient_temp = ambient
            enviro.photometric = photometric
        self.logger.debug('End Environment augmentation.')
        return enviro

    def _get_instrument(self, current):
        """
        Create an Instrument instance populated with available content information.
        :return: Instrument
        """
        self.logger.debug('Begin Instrument augmentation.')
        name = self._get_from_list(
            'Observation.instrument.name', index=0, current=None if current is None else current.name
        )
        keywords = self._get_set_from_list('Observation.instrument.keywords', index=0)
        instr = None
        if name:
            instr = caom2.Instrument(str(name))
            ContentParser._add_keywords(keywords, current, instr)
        self.logger.debug('End Instrument augmentation.')
        return instr

    def _get_members(self, obs):
        """
        Returns the members of a derived observation (if specified)
        :param obs: observation to augment
        :return: members value
        """
        members = None
        self.logger.debug('Begin Members augmentation.')
        if isinstance(obs, caom2.SimpleObservation) and (
            self.blueprint._get('DerivedObservation.members') or self.blueprint._get('CompositeObservation.members')
        ):
            raise TypeError('Cannot apply blueprint for DerivedObservation to a ' 'simple observation')
        elif isinstance(obs, caom2.DerivedObservation):
            lookup = self.blueprint._get('DerivedObservation.members', extension=1)
            if ObsBlueprint.is_table(lookup):
                *_, extension = lookup
                member_list = self._get_from_table('DerivedObservation.members', int(extension))
                # ensure the members are good little ObservationURIs
                if member_list.startswith('caom:'):
                    members = member_list
                else:
                    members = ' '.join(
                        [
                            'caom:{}/{}'.format(obs.collection, i) if not i.startswith('caom') else i
                            for i in member_list.split()
                        ]
                    )
            else:
                if obs.members is None:
                    members = self._get_from_list('DerivedObservation.members', index=0)
                else:
                    members = self._get_from_list('DerivedObservation.members', index=0, current=obs.members)
        elif isinstance(obs, caom2.CompositeObservation):
            lookup = self.blueprint._get('CompositeObservation.members', extension=1)
            if ObsBlueprint.is_table(lookup):
                *_, extension = lookup
                member_list = self._get_from_table('CompositeObservation.members', int(extension))
                # ensure the members are good little ObservationURIs
                if member_list.startswith('caom:'):
                    members = member_list
                else:
                    members = ' '.join(
                        [
                            'caom:{}/{}'.format(obs.collection, i) if not i.startswith('caom') else i
                            for i in member_list.split()
                        ]
                    )
            else:
                if obs.members is None:
                    members = self._get_from_list('CompositeObservation.members', index=0)
                else:
                    members = self._get_from_list('CompositeObservation.members', index=0, current=obs.members)
        self.logger.debug('End Members augmentation.')
        return members

    def _get_axis_wcs(self, label, wcs, index):
        """Helper function to construct a CoordAxis1D instance, with all
        it's members, from the blueprint.

        :param label: axis name - must be one of 'custom', 'energy', 'time', or 'polarization', as it's used for the
            blueprint lookup.
        :param index: which blueprint index to find a value in
        :return an instance of CoordAxis1D
        """
        self.logger.debug(f'Begin {label} axis construction from blueprint.')

        aug_axis = None
        aug_error = None
        aug_axis_ctype = self._get_from_list(f'Chunk.{label}.axis.axis.ctype', index)
        aug_axis_cunit = self._get_from_list(f'Chunk.{label}.axis.axis.cunit', index)
        if aug_axis_ctype is not None:
            aug_axis = caom2.Axis(aug_axis_ctype, aug_axis_cunit)
            self.logger.debug(f'Creating {label} Axis for {self.uri} from blueprint')

        aug_error = self._two_param_constructor(
            f'Chunk.{label}.axis.error.syser', f'Chunk.{label}.axis.error.rnder', index, _to_float, caom2.CoordError,
        )

        aug_naxis = None
        aug_range = self._try_range(index, label)
        aug_naxis_index = None
        if aug_axis is not None:
            if aug_range is None:
                self.logger.debug(f'Try Function construction since Range construction failed for {label}.')
                if wcs is None or wcs.axis is None or wcs.axis.function is None:
                    aug_ref_coord = self._two_param_constructor(
                        f'Chunk.{label}.axis.function.refCoord.pix',
                        f'Chunk.{label}.axis.function.refCoord.val',
                        index,
                        _to_float,
                        caom2.RefCoord,
                    )
                    aug_delta = _to_float(self._get_from_list(f'Chunk.{label}.axis.function.delta', index))
                    aug_length = _to_int(self._get_from_list(f'Chunk.{label}.axis.function.naxis', index))
                    aug_function = None
                    if aug_length is not None and aug_delta is not None and aug_ref_coord is not None:
                        aug_function = caom2.CoordFunction1D(aug_length, aug_delta, aug_ref_coord)
                    aug_naxis = caom2.CoordAxis1D(aug_axis, aug_error, None, None, aug_function)
                    if aug_function is not None:
                        # if the WCS is described with a Function, cutouts can be supported, so specify an axis
                        aug_naxis_index = _to_int(self._get_from_list(f'Chunk.{label}Axis', index))
                    self.logger.debug(f'Creating function {label} CoordAxis1D for {self.uri} from blueprint')
            else:
                aug_naxis = caom2.CoordAxis1D(axis=aug_axis, error=aug_error, range=aug_range)
                self.logger.debug(f'Creating range {label} CoordAxis1D for {self.uri} from blueprint')

        self.logger.debug(f'End {label} axis construction from blueprint.')
        return aug_naxis, aug_naxis_index

    def _get_proposal(self, current):
        """
        Create a Proposal instance populated with available content
        information.
        :return: Proposal
        """
        self.logger.debug('Begin Proposal augmentation.')
        prop_id = self._get_from_list(
            'Observation.proposal.id', index=0, current=None if current is None else current.id
        )
        pi = self._get_from_list(
            'Observation.proposal.pi', index=0, current=None if current is None else current.pi_name
        )
        project = self._get_from_list(
            'Observation.proposal.project', index=0, current=None if current is None else current.project
        )
        title = self._get_from_list(
            'Observation.proposal.title', index=0, current=None if current is None else current.title
        )
        keywords = self._get_set_from_list('Observation.proposal.keywords', index=0)
        proposal = current
        if prop_id:
            proposal = caom2.Proposal(str(prop_id), pi, project, title)
            ContentParser._add_keywords(keywords, current, proposal)
        self.logger.debug(f'End Proposal augmentation {prop_id}.')
        return proposal

    def _get_requirements(self, current):
        """
        Create a Requirements instance populated with available content
        information.
        :return: Requirements
        """
        self.logger.debug('Begin Requirement augmentation.')
        flag = self._get_from_list(
            'Observation.requirements.flag', index=0, current=None if current is None else current.flag
        )
        reqts = caom2.Requirements(flag) if flag else None
        self.logger.debug('End Requirement augmentation.')
        return reqts

    def _get_target(self, current):
        """
        Create a Target instance populated with available content information.
        :return: Target
        """
        self.logger.debug('Begin Target augmentation.')
        name = self._get_from_list(
            'Observation.target.name', index=0, current=None if current is None else current.name
        )
        target_type = self._get_from_list(
            'Observation.target.type', index=0, current=None if current is None else current.target_type
        )
        standard = self._cast_as_bool(
            self._get_from_list(
                'Observation.target.standard', index=0, current=None if current is None else current.standard
            )
        )
        redshift = self._get_from_list(
            'Observation.target.redshift', index=0, current=None if current is None else current.redshift
        )
        keywords = self._get_set_from_list('Observation.target.keywords', index=0)
        moving = self._cast_as_bool(
            self._get_from_list(
                'Observation.target.moving', index=0, current=None if current is None else current.moving
            )
        )
        target_id = _to_str(
            self._get_from_list(
                'Observation.target.targetID', index=0, current=None if current is None else current.target_id
            )
        )
        target = None
        if name:
            target = caom2.Target(str(name), target_type, standard, redshift, moving=moving, target_id=target_id)
            ContentParser._add_keywords(keywords, current, target)
        self.logger.debug('End Target augmentation.')
        return target

    def _get_target_position(self, current):
        """
        Create a Target Position instance populated with available content information.
        :return: Target Position
        """
        self.logger.debug('Begin CAOM2 TargetPosition augmentation.')
        x = self._get_from_list(
            'Observation.target_position.point.cval1',
            index=0,
            current=None if current is None else current.coordinates.cval1,
        )
        y = self._get_from_list(
            'Observation.target_position.point.cval2',
            index=0,
            current=None if current is None else current.coordinates.cval2,
        )
        coordsys = self._get_from_list(
            'Observation.target_position.coordsys', index=0, current=None if current is None else current.coordsys
        )
        equinox = self._get_from_list(
            'Observation.target_position.equinox', index=0, current=None if current is None else current.equinox
        )
        aug_target_position = None
        if x and y:
            aug_point = caom2.Point(x, y)
            aug_target_position = caom2.TargetPosition(aug_point, coordsys)
            aug_target_position.equinox = _to_float(equinox)
        self.logger.debug('End CAOM2 TargetPosition augmentation.')
        return aug_target_position

    def _get_telescope(self, current):
        """
        Create a Telescope instance populated with available content information.
        :return: Telescope
        """
        self.logger.debug('Begin Telescope augmentation.')
        name = self._get_from_list(
            'Observation.telescope.name', index=0, current=None if current is None else current.name
        )
        geo_x = _to_float(
            self._get_from_list(
                'Observation.telescope.geoLocationX',
                index=0,
                current=None if current is None else current.geo_location_x,
            )
        )
        geo_y = _to_float(
            self._get_from_list(
                'Observation.telescope.geoLocationY',
                index=0,
                current=None if current is None else current.geo_location_y,
            )
        )
        geo_z = _to_float(
            self._get_from_list(
                'Observation.telescope.geoLocationZ',
                index=0,
                current=None if current is None else current.geo_location_z,
            )
        )
        keywords = self._get_set_from_list('Observation.telescope.keywords', index=0)
        aug_tel = None
        if name:
            aug_tel = caom2.Telescope(str(name), geo_x, geo_y, geo_z)
            ContentParser._add_keywords(keywords, current, aug_tel)
        self.logger.debug('End Telescope augmentation.')
        return aug_tel

    def _cast_as_bool(self, from_value):
        """
        Make lower case Java booleans into capitalized python booleans.
        :param from_value: Something that represents a boolean value
        :return: a python boolean value
        """
        if isinstance(from_value, bool):
            return from_value
        result = None
        # so far, these are the only options that are coming in from the config files - may need to add more as
        # more types are experienced
        if from_value == 'false':
            result = False
        elif from_value == 'true':
            result = True
        return result

    def _try_custom_with_blueprint(self, chunk, index):
        """
        A mechanism to augment the Custom WCS completely from the blueprint. Do nothing if the WCS information cannot
        be correctly created.

        :param chunk: The chunk to modify with the addition of custom information.
        :param index: The index in the blueprint for looking up plan information.
        """
        self.logger.debug('Begin augmentation with blueprint for custom.')
        aug_naxis, aug_naxis_index = self._get_axis_wcs('custom', chunk.custom, index)
        if aug_naxis is None:
            self.logger.debug('No blueprint custom information.')
        else:
            # always create a new CustomWCS instance because there's no setter for 'axis' parameter
            chunk.custom = caom2.CustomWCS(aug_naxis)
            chunk.custom_axis = aug_naxis_index
            self.logger.debug(f'Updating CustomWCS for {self.uri}.')
        self.logger.debug('End augmentation with blueprint for custom.')

    def _try_energy_with_blueprint(self, chunk, index):
        """
        A mechanism to augment the Energy WCS completely from the blueprint. Do nothing if the WCS information
        cannot be correctly created.

        :param chunk: The chunk to modify with the addition of energy information.
        :param index: The index in the blueprint for looking up plan information.
        """
        self.logger.debug(f'Begin augmentation with blueprint for energy with index {index}.')
        aug_axis, aug_naxis_index = self._get_axis_wcs('energy', chunk.energy, index)
        specsys = _to_str(self._get_from_list('Chunk.energy.specsys', index))
        if aug_axis is None:
            self.logger.debug('No blueprint energy information.')
        else:
            if chunk.energy:
                chunk.energy.axis = aug_axis
                chunk.energy.specsys = specsys
            else:
                chunk.energy = caom2.SpectralWCS(aug_axis, specsys)
                self.logger.debug(f'Creating SpectralWCS for {self.uri} from blueprint')
            chunk.energy_axis = aug_naxis_index

        if chunk.energy:
            chunk.energy.ssysobs = self._get_from_list('Chunk.energy.ssysobs', index, chunk.energy.ssysobs)
            chunk.energy.restfrq = self._get_from_list('Chunk.energy.restfrq', index, chunk.energy.restfrq)
            chunk.energy.restwav = self._get_from_list('Chunk.energy.restwav', index, chunk.energy.restwav)
            chunk.energy.velosys = self._get_from_list('Chunk.energy.velosys', index, chunk.energy.velosys)
            chunk.energy.zsource = self._get_from_list('Chunk.energy.zsource', index, chunk.energy.zsource)
            chunk.energy.ssyssrc = self._get_from_list('Chunk.energy.ssyssrc', index, chunk.energy.ssyssrc)
            chunk.energy.velang = self._get_from_list('Chunk.energy.velang', index, chunk.energy.velang)
            chunk.energy.bandpass_name = self._get_from_list(
                'Chunk.energy.bandpassName', index, chunk.energy.bandpass_name
            )
            chunk.energy.transition = self._get_energy_transition(chunk.energy.transition)
            chunk.energy.resolving_power = _to_float(
                self._get_from_list('Chunk.energy.resolvingPower', index, chunk.energy.resolving_power)
            )
        self.logger.debug('End augmentation with blueprint for energy.')

    def _try_observable_with_blueprint(self, chunk, index):
        """
        A mechanism to augment the Observable WCS completely from the blueprint. Do nothing if the WCS information
        cannot be correctly created.

        :param chunk: The chunk to modify with the addition of observable information.
        :param index: The index in the blueprint for looking up plan information.
        """
        self.logger.debug('Begin augmentation with blueprint for ' 'observable.')
        aug_axis = self._two_param_constructor(
            'Chunk.observable.dependent.axis.ctype',
            'Chunk.observable.dependent.axis.cunit',
            index,
            _to_str,
            caom2.Axis,
        )
        aug_bin = _to_int(self._get_from_list('Chunk.observable.dependent.bin', index))
        if aug_axis is not None and aug_bin is not None:
            chunk.observable = caom2.ObservableAxis(caom2.Slice(aug_axis, aug_bin))
            chunk.observable_axis = _to_int(self._get_from_list('Chunk.observableAxis', index))
        self.logger.debug('End augmentation with blueprint for polarization.')

    def _try_polarization_with_blueprint(self, chunk, index):
        """
        A mechanism to augment the Polarization WCS completely from the blueprint. Do nothing if the WCS information
        cannot be correctly created.

        :param chunk: The chunk to modify with the addition of polarization information.
        :param index: The index in the blueprint for looking up plan information.
        """
        self.logger.debug('Begin augmentation with blueprint for ' 'polarization.')
        aug_axis, aug_naxis_index = self._get_axis_wcs('polarization', chunk.polarization, index)
        if aug_axis is not None:
            if chunk.polarization:
                chunk.polarization.axis = aug_axis
            else:
                chunk.polarization = caom2.PolarizationWCS(aug_axis)
                self.logger.debug(f'Creating PolarizationWCS for {self.uri} from blueprint')
            chunk.polarization_axis = aug_naxis_index

        self.logger.debug('End augmentation with blueprint for polarization.')

    def _try_position_range(self, index):
        self.logger.debug('Try to set the range for position from blueprint, since there is no function')
        aug_range = None
        aug_range_c1_start = self._two_param_constructor(
            'Chunk.position.axis.range.start.coord1.pix',
            'Chunk.position.axis.range.start.coord1.val',
            index,
            _to_float,
            caom2.RefCoord,
        )
        aug_range_c1_end = self._two_param_constructor(
            'Chunk.position.axis.range.end.coord1.pix',
            'Chunk.position.axis.range.end.coord1.val',
            index,
            _to_float,
            caom2.RefCoord,
        )
        aug_range_c2_start = self._two_param_constructor(
            'Chunk.position.axis.range.start.coord2.pix',
            'Chunk.position.axis.range.start.coord2.val',
            index,
            _to_float,
            caom2.RefCoord,
        )
        aug_range_c2_end = self._two_param_constructor(
            'Chunk.position.axis.range.end.coord2.pix',
            'Chunk.position.axis.range.end.coord2.val',
            index,
            _to_float,
            caom2.RefCoord,
        )
        if aug_range_c1_start and aug_range_c1_end and aug_range_c2_start and aug_range_c2_end:
            aug_range = caom2.CoordRange2D(
                caom2.Coord2D(aug_range_c1_start, aug_range_c1_end),
                caom2.Coord2D(aug_range_c2_start, aug_range_c2_end),
            )
            self.logger.debug('Completed setting range for position')
        return aug_range

    def _try_position_with_blueprint(self, chunk, index):
        """
        A mechanism to augment the Position WCS completely from the blueprint. Do nothing if the WCS information
        cannot be correctly created.

        :param chunk: The chunk to modify with the addition of position information.
        :param index: The index in the blueprint for looking up plan information.
        """
        self.logger.debug('Begin augmentation with blueprint for position.')
        aug_axis = None
        if (
            chunk.position is not None
            and chunk.position.axis is not None
            and chunk.position.axis.axis1 is not None
            and chunk.position.axis.axis2 is not None
        ):
            # preserve the values obtained from file data
            aug_x_axis = chunk.position.axis.axis1
            aug_y_axis = chunk.position.axis.axis2
            aug_x_error = chunk.position.axis.error1
            aug_y_error = chunk.position.axis.error2
        else:
            aug_x_axis = self._two_param_constructor(
                'Chunk.position.axis.axis1.ctype', 'Chunk.position.axis.axis1.cunit', index, _to_str, caom2.Axis
            )
            aug_y_axis = self._two_param_constructor(
                'Chunk.position.axis.axis2.ctype', 'Chunk.position.axis.axis2.cunit', index, _to_str, caom2.Axis
            )
            aug_x_error = self._two_param_constructor(
                'Chunk.position.axis.error1.syser',
                'Chunk.position.axis.error1.rnder',
                index,
                _to_float,
                caom2.CoordError,
            )
            aug_y_error = self._two_param_constructor(
                'Chunk.position.axis.error2.syser',
                'Chunk.position.axis.error2.rnder',
                index,
                _to_float,
                caom2.CoordError,
            )
        aug_range = self._try_position_range(index)
        if aug_range is None:
            if chunk.position is None or chunk.position.axis is None or chunk.position.axis.function is None:
                aug_dimension = self._two_param_constructor(
                    'Chunk.position.axis.function.dimension.naxis1',
                    'Chunk.position.axis.function.dimension.naxis2',
                    index,
                    _to_int,
                    caom2.Dimension2D,
                )
                aug_x_ref_coord = self._two_param_constructor(
                    'Chunk.position.axis.function.refCoord.coord1.pix',
                    'Chunk.position.axis.function.refCoord.coord1.val',
                    index,
                    _to_float,
                    caom2.RefCoord,
                )
                aug_y_ref_coord = self._two_param_constructor(
                    'Chunk.position.axis.function.refCoord.coord2.pix',
                    'Chunk.position.axis.function.refCoord.coord2.val',
                    index,
                    _to_float,
                    caom2.RefCoord,
                )
                aug_cd11 = _to_float(self._get_from_list('Chunk.position.axis.function.cd11', index))
                aug_cd12 = _to_float(self._get_from_list('Chunk.position.axis.function.cd12', index))
                aug_cd21 = _to_float(self._get_from_list('Chunk.position.axis.function.cd21', index))
                aug_cd22 = _to_float(self._get_from_list('Chunk.position.axis.function.cd22', index))

                aug_ref_coord = None
                if aug_x_ref_coord is not None and aug_y_ref_coord is not None:
                    aug_ref_coord = caom2.Coord2D(aug_x_ref_coord, aug_y_ref_coord)
                    self.logger.debug(f'Creating position Coord2D for {self.uri}')

                aug_function = None
                if (
                    aug_dimension is not None
                    and aug_ref_coord is not None
                    and aug_cd11 is not None
                    and aug_cd12 is not None
                    and aug_cd21 is not None
                    and aug_cd22 is not None
                ):
                    aug_function = caom2.CoordFunction2D(
                        aug_dimension, aug_ref_coord, aug_cd11, aug_cd12, aug_cd21, aug_cd22
                    )
                    self.logger.debug(f'Creating position CoordFunction2D for {self.uri}')

                if aug_x_axis is not None and aug_y_axis is not None and aug_function is not None:
                    aug_axis = caom2.CoordAxis2D(
                        aug_x_axis, aug_y_axis, aug_x_error, aug_y_error, None, None, aug_function
                    )
                    self.logger.debug(f'Creating position CoordAxis2D for {self.uri}')

                    chunk.position_axis_1 = _to_int(self._get_from_list('Chunk.positionAxis1', index))
                    chunk.position_axis_2 = _to_int(self._get_from_list('Chunk.positionAxis2', index))
        else:
            aug_axis = caom2.CoordAxis2D(aug_x_axis, aug_y_axis, aug_x_error, aug_y_error, range=aug_range)

        if aug_axis is not None:
            if chunk.position:
                chunk.position.axis = aug_axis
            else:
                chunk.position = caom2.SpatialWCS(aug_axis)
                self.logger.debug(f'Creating SpatialWCS for {self.uri} from blueprint')

        if chunk.position:
            chunk.position.coordsys = self._get_from_list('Chunk.position.coordsys', index, chunk.position.coordsys)
            chunk.position.equinox = _to_float(
                self._get_from_list('Chunk.position.equinox', index, chunk.position.equinox)
            )
            chunk.position.resolution = self._get_from_list(
                'Chunk.position.resolution', index, chunk.position.resolution
            )
        self.logger.debug('End augmentation with blueprint for position.')

    def _try_range(self, index, lookup):
        self.logger.debug(f'Try to set the range for {lookup}')
        result = None
        aug_range_start = self._two_param_constructor(
            f'Chunk.{lookup}.axis.range.start.pix',
            f'Chunk.{lookup}.axis.range.start.val',
            index,
            _to_float,
            caom2.RefCoord,
        )
        aug_range_end = self._two_param_constructor(
            f'Chunk.{lookup}.axis.range.end.pix',
            f'Chunk.{lookup}.axis.range.end.val',
            index,
            _to_float,
            caom2.RefCoord,
        )
        if aug_range_start and aug_range_end:
            result = caom2.CoordRange1D(aug_range_start, aug_range_end)
            self.logger.debug(f'Completed setting range with return for {lookup}')
        return result

    def _try_time_with_blueprint(self, chunk, index):
        """
        A mechanism to augment the Time WCS completely from the blueprint. Do nothing if the WCS information cannot
        be correctly created.

        :param chunk: The chunk to modify with the addition of time information.
        :param index: The index in the blueprint for looking up plan information.
        """
        self.logger.debug('Begin augmentation with blueprint for temporal.')

        aug_axis, aug_axis_index = self._get_axis_wcs('time', chunk.time, index)
        if aug_axis is not None:
            if chunk.time:
                chunk.time.axis = aug_axis
            else:
                chunk.time = caom2.TemporalWCS(aug_axis)
                self.logger.debug(f'Creating TemporalWCS for {self.uri} from blueprint')
            chunk.time_axis = aug_axis_index

        if chunk.time:
            chunk.time.exposure = _to_float(self._get_from_list('Chunk.time.exposure', index, chunk.time.exposure))
            chunk.time.resolution = _to_float(
                self._get_from_list('Chunk.time.resolution', index, chunk.time.resolution)
            )
            chunk.time.timesys = _to_str(self._get_from_list('Chunk.time.timesys', index, chunk.time.timesys))
            chunk.time.trefpos = self._get_from_list('Chunk.time.trefpos', index, chunk.time.trefpos)
            chunk.time.mjdref = self._get_from_list('Chunk.time.mjdref', index, chunk.time.mjdref)

        self.logger.debug('End augmentation with blueprint for temporal.')

    def _two_param_constructor(self, lookup1, lookup2, index, to_type, ctor):
        """
        Helper function to build from the blueprint, a CAOM2 entity that has two required parameters.

        :param lookup1: Blueprint lookup text for the first constructor parameter.
        :param lookup2: Blueprint lookup text for the second constructor parameter.
        :param index:  Which index in the blueprint to do the lookup on.
        :param to_type: Function to cast the blueprint value to a particular type.
        :param ctor: The constructor that has two parameters to build.
        :return: The instance returned by the constructor, or None if any of the values are undefined.
        """
        param1 = to_type(self._get_from_list(lookup1, index))
        param2 = to_type(self._get_from_list(lookup2, index))
        new_object = None
        if param1 is not None and param2 is not None:
            new_object = ctor(param1, param2)
        return new_object

    # TODO - is this the right implementation?
    def add_parts(self, artifact, index=0):
        result = False
        if self.blueprint.has_chunk(index):
            artifact.parts.add(caom2.Part(str(index)))
            result = True
        return result

    @staticmethod
    def _add_keywords(keywords, current, to_set):
        """
        Common code for adding keywords to a CAOM2 entity, capturing all the weird metadata cases that happen at CADC.

        :param keywords: Keywords to add to a CAOM2 set.
        :param current: Existing CAOM2 entity with a keywords attribute.
        :param to_set: A CAOM2 entity with a keywords attribute.
        """
        if keywords:
            if isinstance(keywords, set):
                to_set.keywords.update(keywords)
            else:
                for k in keywords.split():
                    to_set.keywords.add(k)
        else:
            if current is not None:
                # preserve the original value
                to_set.keywords.update(current.keywords)
        if to_set.keywords is not None and None in to_set.keywords:
            to_set.keywords.remove(None)
        if to_set.keywords is not None and 'none' in to_set.keywords:
            to_set.keywords.remove('none')


class FitsParser(ContentParser):
    """
    Parses a FITS file and extracts the CAOM2 related information which can be used to augment an existing CAOM2
    observation, plane or artifact. The constructor takes either a FITS file as argument or a list of dictionaries
    (FITS keyword=value) corresponding to each extension.

    The WCS-related keywords of the FITS file are consumed by the astropy.wcs package which might display warnings
    with regards to compliance.

    Example 1:
    parser = FitsParser(input = '/staging/700000o.fits.gz')
    ...
    # customize parser.headers by deleting, changing or adding attributes

    obs = Observation(collection='TEST', observation_id='700000', algorithm='exposure')
    plane = Plane(plane_id='700000-1')
    obs.plane.add(plane)

    artifact = Artifact(uri='ad:CFHT/700000o.fits.gz', product_type='science', release_type='data')
    plane.artifacts.add(artifact)

    parser.augment_observation(obs)

    # further update obs


    Example 2:

    headers = [] # list of dictionaries headers
    # populate headers
    parser = FitsParser(input=headers)

    parser.augment_observation(obs)
    ...

    """

    def __init__(self, src, obs_blueprint=None, uri=None, extension_start_index=0, extension_end_index=None):
        """
        Ctor
        :param src: List of headers (dictionary of FITS keywords:value) with one header for each extension or a FITS
            input file.
        :param obs_blueprint: externally provided blueprint
        :param uri: which artifact augmentation is based on
        """
        self.logger = logging.getLogger(__name__)
        self._wcs_parsers = {}
        self._headers = []
        self.parts = 0
        self.file = ''
        if isinstance(src, list):
            # assume this is the list of headers
            self._headers = src
        else:
            # assume file
            self.file = src
            self._headers = data_util.get_local_headers_from_fits(self.file)
        if obs_blueprint:
            self._blueprint = obs_blueprint
        else:
            self._blueprint = ObsBlueprint()
        self._errors = []
        # for command-line parameter to module execution
        self.uri = uri
        self._extension_start_index = extension_start_index
        self._extension_end_index = extension_end_index if extension_end_index is not None else self._get_num_parts()
        self.apply_blueprint()

    def _get_num_parts(self):
        """return the number of Parts to create for a CAOM record
        """
        return len(self._headers)

    @property
    def headers(self):
        """
        List of headers where each header should allow dictionary like access to the FITS attribute in that header
        :return:
        """
        return self._headers

    def add_parts(self, artifact, index):
        # there is one Part per extension, the name is the extension number
        if FitsParser._has_data_array(self._headers[index]) and self.blueprint.has_chunk(index):
            if str(index) not in artifact.parts.keys():
                # TODO use extension name?
                artifact.parts.add(caom2.Part(str(index)))
                self.logger.debug(f'Part created for HDU {index}.')
            result = True
        else:
            artifact.parts.add(caom2.Part(str(index)))
            self.logger.debug(f'Create empty part for HDU {index}')
            result = False
        return result

    def apply_blueprint(self):
        self.logger.debug(f'Begin apply_blueprint {self.uri}')
        # pointers that are short to type
        exts = self.blueprint._extensions
        wcs_std = self.blueprint._wcs_std
        plan = self.blueprint._plan

        # firstly, apply the functions
        if self.blueprint._module is not None or self.blueprint._module_instance is not None:
            for key, value in plan.items():
                if ObsBlueprint.is_function(value):
                    if self._blueprint._module_instance is None:
                        plan[key] = self._execute_external(value, key, 0)
                    else:
                        plan[key] = self._execute_external_instance(value, key, 0)
            for extension in exts:
                for key, value in exts[extension].items():
                    if ObsBlueprint.is_function(value):
                        if self._blueprint._module_instance is None:
                            exts[extension][key] = self._execute_external(value, key, extension)
                        else:
                            exts[extension][key] = self._execute_external_instance(value, key, extension)

        # apply overrides from blueprint to all extensions
        for key, value in plan.items():
            if key in wcs_std:
                if ObsBlueprint.needs_lookup(value):
                    # alternative attributes provided for standard wcs attrib.
                    for header in self.headers:
                        for v in value[0]:
                            if v in header and v not in wcs_std[key].split(','):
                                keywords = wcs_std[key].split(',')
                                for keyword in keywords:
                                    _set_by_type(header, keyword, str(header[v]))
                elif ObsBlueprint.is_function(value):
                    continue
                elif ObsBlueprint.has_no_value(value):
                    continue
                else:
                    # value provided for standard wcs attribute
                    if ObsBlueprint.needs_lookup(wcs_std[key]):
                        keywords = wcs_std[key][0]
                    elif ObsBlueprint.is_function(wcs_std[key]):
                        continue
                    else:
                        keywords = wcs_std[key].split(',')
                    for keyword in keywords:
                        for header in self.headers:
                            _set_by_type(header, keyword, str(value))

        # apply overrides to the remaining extensions
        for extension in exts:
            if extension >= len(self.headers):
                logging.error('More extensions configured {} than headers ' '{}'.format(extension, len(self.headers)))
                continue
            hdr = self.headers[extension]
            for key, value in exts[extension].items():
                if ObsBlueprint.needs_lookup(value):
                    # alternative attributes provided for standard wcs attrib.
                    for v in value[0]:
                        if v in hdr and v not in wcs_std[key].split(','):
                            keywords = wcs_std[key].split(',')
                            for keyword in keywords:
                                _set_by_type(hdr, keyword, str(hdr[v]))
                elif ObsBlueprint.is_table(value):
                    continue
                elif ObsBlueprint.has_no_value(value):
                    continue
                else:
                    if key in wcs_std.keys():
                        keywords = wcs_std[key].split(',')
                        for keyword in keywords:
                            _set_by_type(hdr, keyword, value)
                            logging.debug(f'{keyword}: set to {value} in extension {extension}')
                    else:
                        exts[extension][key] = value
        # apply defaults to all extensions
        for key, value in plan.items():
            if ObsBlueprint.has_default_value(value):
                for index, header in enumerate(self.headers):
                    for keywords in value[0]:
                        for keyword in keywords.split(','):
                            if (
                                not header.get(keyword.strip())
                                and keyword == keywords
                                and keywords == value[0][-1]  # checking a string
                            ):  # last item
                                # apply a default if a value does not already exist, and all possible values of
                                # keywords have been checked
                                _set_by_type(header, keyword.strip(), value[1])
                                logging.debug(
                                    '{}: set default value of {} in HDU {}.'.format(keyword, value[1], index)
                                )

        # TODO wcs in astropy ignores cdelt attributes when it finds a cd attribute even if it's in a different axis
        for header in self.headers:
            cd_present = False
            for i in range(1, 6):
                if 'CD{0}_{0}'.format(i) in header:
                    cd_present = True
                    break
            if cd_present:
                for i in range(1, 6):
                    if f'CDELT{i}' in header and 'CD{0}_{0}'.format(i) not in header:
                        header['CD{0}_{0}'.format(i)] = header[f'CDELT{i}']

        # TODO When a projection is specified, wcslib expects corresponding DP arguments with NAXES attributes.
        # Normally, omitting the attribute signals no distortion which is the assumption in caom2blueprint for
        # energy and polarization axes. Following is a workaround for SIP projections. For more details see:
        # http://www.atnf.csiro.au/people/mcalabre/WCS/dcs_20040422.pdf
        for header in self.headers:
            sip = False
            for i in range(1, 6):
                if (
                    (f'CTYPE{i}' in header)
                    and isinstance(header[f'CTYPE{i}'], str)
                    and ('-SIP' in header[f'CTYPE{i}'])
                ):
                    sip = True
                    break
            if sip:
                for i in range(1, 6):
                    if (f'CTYPE{i}' in header) and ('-SIP' not in header[f'CTYPE{i}']) and (f'DP{i}' not in header):
                        header[f'DP{i}'] = 'NAXES: 1'
        return

    def augment_artifact(self, artifact):
        """
        Augments a given CAOM2 artifact with available FITS information
        :param artifact: existing CAOM2 artifact to be augmented
        """
        self.logger.debug(f'Begin artifact augmentation for {artifact.uri} with {len(self.headers)} HDUs.')

        if self.blueprint.get_configed_axes_count() == 0:
            raise TypeError('No WCS Data. End artifact augmentation for {}.'.format(artifact.uri))

        for i, header in enumerate(self.headers):
            if not self.add_parts(artifact, i):
                # artifact-level attributes still require updating
                BlueprintParser.augment_artifact(self, artifact)
                continue
            self._wcs_parsers[i] = FitsWcsParser(header, self.file, str(i))
        super().augment_artifact(artifact)

        self.logger.debug(f'End artifact augmentation for {artifact.uri}.')

    def _get_chunk_naxis(self, chunk, index=None):
        # NOTE: astropy.wcs does not distinguished between WCS axes and data array axes. naxis in astropy.wcs
        # represents in fact the number of WCS axes, whereas chunk.axis represents the naxis of the data array.
        # Solution is to determine it directly from the header
        if 'ZNAXIS' in self._headers[index]:
            chunk.naxis = _to_int(self._headers[index]['ZNAXIS'])
        elif 'NAXIS' in self._headers[index]:
            chunk.naxis = _to_int(self._headers[index]['NAXIS'])
        else:
            super()._get_chunk_naxis(chunk)

    def _get_from_list(self, lookup, index, current=None):
        value = None
        try:
            keys = self.blueprint._get(lookup, index)
        except KeyError:
            self.add_error(lookup, sys.exc_info()[1])
            self.logger.debug(f'Could not find {lookup!r} in caom2blueprint configuration.')
            if current:
                self.logger.debug(f'{lookup}: using current value of {current!r}.')
                value = current
            return value

        if ObsBlueprint.needs_lookup(keys):
            for ii in keys[0]:
                try:
                    value = self.headers[index].get(ii)
                    if value:
                        self.logger.debug(f'{lookup}: assigned value {value} based on ' f'keyword {ii}.')
                        break
                except (KeyError, IndexError):
                    if keys[0].index(ii) == len(keys[0]) - 1:
                        self.add_error(lookup, sys.exc_info()[1])
                    # assign a default value, if one exists
                    if keys[1]:
                        if current is None:
                            value = keys[1]
                            self.logger.debug(f'{lookup}: assigned default value {value}.')
                        else:
                            value = current
            if value is None:
                # checking current does not work in the general case, because current might legitimately be 'None'
                if self._blueprint.update:
                    if current is not None or (current is None and isinstance(value, bool)):
                        value = current
                        self.logger.debug(f'{lookup}: used current value {value}.')
                else:
                    # assign a default value, if one exists
                    if keys[1]:
                        if current is None:
                            value = keys[1]
                            self.logger.debug(f'{lookup}: assigned default value {value}.')
                        else:
                            value = current

        elif (keys is not None) and (keys != ''):
            if keys == 'None':
                value = None
            else:
                value = keys
        elif current:
            value = current

        self.logger.debug(f'{lookup}: value is {value}')
        return value

    def _get_from_table(self, lookup, extension):
        """
        Return a space-delimited list of all the row values from a column.

        This is a straight FITS BINTABLE lookup. There is no support for default values. Unless someone provides a
        compelling use case.

        :param lookup: where to find the column name
        :param extension: which extension
        :return: A string, which is a space-delimited list of all the values.
        """
        value = ''
        try:
            keywords = self.blueprint._get(lookup, extension)
        except KeyError as e:
            self.add_error(lookup, sys.exc_info()[1])
            self.logger.debug('Could not find {!r} in fits2caom2 configuration.'.format(lookup))
            raise e

        if isinstance(keywords, tuple) and keywords[0] == 'BINTABLE':
            # BINTABLE, so need to retrieve the data from the file
            if self.file is not None and self.file != '':
                with fits.open(self.file) as fits_data:
                    if fits_data[extension].header['XTENSION'] != 'BINTABLE':
                        raise ValueError(
                            f'Got {fits_data[extension].header["XTENSION"]} when looking for a BINTABLE extension.'
                        )
                    for ii in fits_data[extension].data[keywords[1]]:
                        value = f'{ii} {value}'

        self.logger.debug(f'{lookup}: value is {value}')
        return value

    def _get_set_from_list(self, lookup, index):
        value = None
        keywords = None
        try:
            keywords = self.blueprint._get(lookup)
        except KeyError:
            self.add_error(lookup, sys.exc_info()[1])
            self.logger.debug(f'Could not find \'{lookup}\' in caom2blueprint ' f'configuration.')

        if isinstance(keywords, tuple):
            for ii in keywords[0]:
                try:
                    value = self.headers[index].get(ii)
                    break
                except KeyError:
                    self.add_error(lookup, sys.exc_info()[1])
                    if keywords[1]:
                        value = keywords[1]
                        self.logger.debug('{}: assigned default value {}.'.format(lookup, value))
        elif keywords:
            value = keywords
            self.logger.debug(f'{lookup}: assigned value {value}.')

        return value

    @staticmethod
    def _has_data_array(header):
        """

        :param header:
        :return:
        """
        naxis = 0
        if 'ZNAXIS' in header:
            naxis = _to_int(header['ZNAXIS'])
        elif 'NAXIS' in header:
            naxis = _to_int(header['NAXIS'])
        if not naxis:
            return False

        data_axes = 0
        for i in range(1, naxis + 1):
            axis = f'NAXIS{i}'
            if axis in header:
                data_axis = _to_int(header[axis])
                if not data_axes:
                    data_axes = data_axis
                else:
                    data_axes = data_axes * data_axis
        if not data_axes:
            return False

        bitpix = 0
        if 'BITPIX' in header:
            bitpix = _to_int(header['BITPIX'])
        if not bitpix:
            return False
        return True


class Hdf5Parser(ContentParser):
    """
    Parses an HDF5 file and extracts the CAOM2 related information which can be used to augment an existing CAOM2
    observation, plane, or artifact.

    If there is per-Chunk metadata in the file, the constructor parameter 'find_roots_here' is the address location
    in the file where the N Chunk metadata starts.

    The WCS-related keywords of the HDF5 files are used to create instances of astropy.wcs.WCS so that verify might
    be called.

    There is no CADC support for the equivalent of the FITS --fhead parameter for HDF5 files, which is why the name
    of the file on a local disk is required.

    How the classes work together for HDF5 files:
    - build an HDF5ObsBlueprint, with _CAOM2_ELEMENT keys, and HDF5 metadata path names as keys
    - cache the metadata from an HDF5 file in the HDF5ObsBlueprint. This caching is done in the
        "apply_blueprint_from_file" method in the Hdf5Parser class, and replaces the path names in the blueprint
        with the values from the HDF5 file. The caching is done so that all HDF5 file access is isolated to one point
        in time.
    - use the cached metadata to build astropy.wcs instances for verification in Hdf5WcsParser.
    - use the astropy.wcs instance and other blueprint metadata to fill the CAOM2 record.
    """

    def __init__(self, obs_blueprint, uri, h5_file, extension_names=None, extension_start_index=0,
                 extension_end_index=None):
        """
        :param obs_blueprint: Hdf5ObsBlueprint instance
        :param uri: which artifact augmentation is based on
        :param h5_file: h5py file handle
        :param extension_names: list of str where Chunk metadata starts. There is one Part/Chunk per list entry
        """
        self._file = h5_file
        # the length of the array is the number of Parts in an HDF5 file,
        # and the values are HDF5 lookup path names.
        self._extension_names = extension_names
        super().__init__(obs_blueprint, uri, extension_start_index, extension_end_index)

    def _get_num_parts(self):
        """return the number of Parts to create for a CAOM record
        """
        result = len(self._blueprint._extensions)
        if result == 0:
            # for HDF5 files, cutouts should be supported in the future, so the minimum is one Part/Chunk construction
            result = 1
        return result

    def _set_wcs_parsers(self, obs_blueprint):
        # used to set the astropy wcs info, resulting in a validated wcs that can be used to construct a valid CAOM2
        # record
        # This method call is over-writing the default behaviour in the ContentParser class. The default behaviour
        # uses the obs_blueprint. This method is called in the ContentParser constructor.
        self._wcs_parsers = {}

    def apply_blueprint_from_file(self):
        """
        Retrieve metadata from file, cache in the blueprint.
        """
        self.logger.debug('Begin apply_blueprint_from_file')
        # h5py is an extra in this package since most collections do not require it
        import h5py

        individual, multi, attributes, candidate_extensions = self._extract_path_names_from_blueprint()
        if self._extension_names is None and len(candidate_extensions) > 0:
            self._find_extension_names(candidate_extensions)
            for index, _ in enumerate(self._extension_names):
                self._blueprint._extensions[index] = {}
        else:
            self._blueprint._extensions[0] = {}
        filtered_individual = [ii for ii in individual.keys() if '(' in ii]

        def _extract_from_item(name, object):
            """
            Function signature dictated by h5py visititems implementation. Executed for each dataset/group in an
            HDF5 file.

            :param name: fully-qualified HDF5 path name
            :param object: what the HDF5 path name points to
            """
            # If it's the Part/Chunk metadata, capture it to extensions.
            # Syntax of the keys described in Hdf5ObsBlueprint class.
            for part_index, part_name in enumerate(self._extension_names):
                if name.startswith(part_name) and isinstance(object, h5py.Dataset) and object.dtype.names is not None:
                    for d_name in object.dtype.names:
                        temp_path = f'{name.replace(part_name, "")}/{d_name}'
                        for path_name in multi.keys():
                            if path_name == temp_path:
                                for jj in multi.get(path_name):
                                    self._blueprint.set(jj, object[d_name], part_index)
                            elif path_name.startswith(temp_path) and '(' in path_name:
                                z = path_name.split('(')
                                if ':' in z[1]:
                                    a = z[1].split(')')[0].split(':')
                                    if len(a) > 2:
                                        raise NotImplementedError
                                    for jj in multi.get(path_name):
                                        self._blueprint.set(
                                            jj,
                                            object[d_name][int(a[0])][int(a[1])],
                                            part_index,
                                        )
                                else:
                                    index = int(z[1].split(')')[0])
                                    for jj in multi.get(path_name):
                                        self._blueprint.set(
                                            jj,
                                            object[d_name][index],
                                            part_index,
                                        )

            # if it's Observation/Plane/Artifact metadata, capture it to the base blueprint
            if isinstance(object, h5py.Dataset):
                if object.dtype.names is not None:
                    for d_name in object.dtype.names:
                        temp = f'//{name}/{d_name}'
                        if temp in individual.keys():
                            for jj in individual.get(temp):
                                self._blueprint.set(jj, object[d_name], 0)
                        else:
                            for ind_path in filtered_individual:
                                if ind_path.startswith(temp):
                                    z = ind_path.split('(')
                                    index = int(z[1].split(')')[0])
                                    for jj in individual.get(ind_path):
                                        self._blueprint.set(jj, object[d_name][index], 0)

        if len(individual) == 0 and len(multi) == 0:
            # CFHT SITELLE
            self.logger.debug(f'attrs for {self.uri}')
            self._extract_from_attrs(attributes)
        else:
            # TAOSII
            self.logger.debug(f'visititems for {self.uri}')
            self._file.visititems(_extract_from_item)
        self.logger.debug('Done apply_blueprint_from_file')

    def _extract_from_attrs(self, attributes):
        # I don't currently see any way to have more than one Part, if relying on attrs for metadata
        part_index = 0
        # v == list of blueprint keys
        for k, v in attributes.items():
            if k in self._file.attrs:
                value = self._file.attrs[k]
                for entry in v:
                    self._blueprint.set(entry, value, part_index)

    def _extract_path_names_from_blueprint(self):
        """
        :return: individual - a dictionary of lists, keys are unique path names for finding metadata once per file.
            Values are _CAOM2_ELEMENT strings.
            multiple - a dictionary of lists, keys are unique path names for finding metadata N times per file. Values
            are _CAOM2_ELEMENT strings.
            attributes - a dictionary of lists, keys reference expected content from the h5py.File().attrs data
            structure and its keys.
            extensions - a list of prefixes for identifying extensions
        """
        individual = defaultdict(list)
        multi = defaultdict(list)
        attributes = defaultdict(list)
        extensions = []
        for key, value in self._blueprint._plan.items():
            if ObsBlueprint.needs_lookup(value):
                for ii in value[0]:
                    if ii.startswith('//'):
                        individual[ii].append(key)
                    elif ii.startswith('/'):
                        if '{}' in ii:
                            bits = ii.split('{}')
                            extensions.append(bits[0])
                            multi[bits[1]].append(key)
                        else:
                            multi[ii].append(key)
                    else:
                        attributes[ii].append(key)

        return individual, multi, attributes, list(set(extensions))

    def _find_extension_names(self, candidates):
        """ if the HDF5 file has a structure where-by more than one Chunk (the equivalent of a FITS HDU extension)
        is defined, try to guess that structure
        """
        candidate_extension_names = []

        def _extract_extension_prefixes(name, object):
            """
            Function signature dictated by h5py visititems implementation. Executed for each dataset/group in an
            HDF5 file.

            :param name: fully-qualified HDF5 path name
            :param object: what the HDF5 path name points to
            """
            for part_name in candidates:
                y = part_name.replace('/', '', 1)
                if name.startswith(y):
                    x = name.split(y)[1].split('/')
                    temp = f'{y}{x[0]}'
                    candidate_extension_names.append(temp)
            self._extension_names = list(sorted(set(candidate_extension_names)))

        self._file.visititems(_extract_extension_prefixes)
        msg = '\n'.join(ii for ii in self._extension_names)
        self.logger.info(f'Found extension_names:\n{msg}')

    def apply_blueprint(self):
        self.logger.debug('Begin apply_blueprint')
        self.apply_blueprint_from_file()

        # after the apply_blueprint_from_file call, all the metadata from the file has been applied to the blueprint,
        # so now do the bits that require no access to file content

        # pointers that are short to type
        exts = self._blueprint._extensions
        plan = self._blueprint._plan

        # apply the functions
        if self._blueprint._module is not None or self._blueprint._module_instance is not None:
            for key, value in plan.items():
                if ObsBlueprint.is_function(value):
                    if self._blueprint._module_instance is None:
                        plan[key] = self._execute_external(value, key, 0)
                    else:
                        plan[key] = self._execute_external_instance(value, key, 0)
            for extension in exts:
                for key, value in exts[extension].items():
                    if ObsBlueprint.is_function(value):
                        if self._blueprint._module_instance is None:
                            exts[extension][key] = self._execute_external(value, key, extension)
                        else:
                            exts[extension][key] = self._execute_external_instance(value, key, extension)

        # apply overrides
        # blueprint already contains all the overrides, only need to make sure the overrides get applied to all the
        # extensions
        for extension in exts:
            for key, value in exts[extension].items():
                if (
                    ObsBlueprint.is_table(value)
                    # already been looked up
                    or ObsBlueprint.needs_lookup(value)
                    # already been executed
                    or ObsBlueprint.is_function(value)
                    # nothing to assign
                    or ObsBlueprint.has_no_value(value)
                ):
                    continue
                exts[extension][key] = value
                self.logger.debug(f'{key}: set to {value} in extension {extension}')

        # apply defaults
        # if no values have been set by file lookups, function execution, or applying overrides, apply defaults,
        # including to all extensions
        for key, value in plan.items():
            if ObsBlueprint.needs_lookup(value) and value[1]:
                # there is a default value in the blueprint that can be used
                for extension in exts:
                    q = exts[extension].get(key)
                    if q is None:
                        exts[extension][key] = value[1]
                        self.logger.debug(
                            f'Add {key} and assign default value of ' f'{value[1]} in extension {extension}.'
                        )
                    elif ObsBlueprint.needs_lookup(value):
                        exts[extension][key] = value[1]
                        self.logger.debug(f'{key}: set value to default of {value[1]} in ' f'extension {extension}.')
                plan[key] = value[1]
                self.logger.debug(f'{key}: set value to default of {value[1]}')

        self.logger.debug('Done apply_blueprint')
        return

    def augment_artifact(self, artifact):
        for ii in range(self._extension_start_index, self._extension_end_index):
            # one WCS parser per Part/Chunk
            self._wcs_parsers[ii] = Hdf5WcsParser(self._blueprint, ii)
        super().augment_artifact(artifact)

    def _get_chunk_naxis(self, chunk, index):
        chunk.naxis = self._get_from_list('Chunk.naxis', index, chunk.naxis)

    def add_parts(self, artifact, index=0):
        artifact.parts.add(caom2.Part(str(index)))
        return True


def _set_by_type(header, keyword, value):
    """astropy documentation says that the type of the second parameter in the 'set' call is 'str', and then warns
    of expectations for floating-point values when the code does that, so make float values into floats, and int
    values into ints."""
    float_value = None
    int_value = None
    if value is not None:
        try:
            float_value = float(value)
        except ValueError:
            pass

        try:
            int_value = int(value)
        except ValueError:
            pass
    if float_value and not str(value).isdecimal() or re.match(r'0\.0*', str(value)):
        header.set(keyword, float_value)
    elif int_value:
        header.set(keyword, int_value)
    else:
        header.set(keyword, value)


def _to_checksum_uri(value):
    if value is None:
        return None
    elif isinstance(value, caom2.ChecksumURI):
        return value
    else:
        return caom2.ChecksumURI(value)
