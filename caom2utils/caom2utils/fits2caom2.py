# -*- coding: utf-8 -*-
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
from builtins import str
from datetime import datetime

import math
from astropy.wcs import WCS
from astropy.io import fits
from cadcutils import version
from caom2 import Artifact, Part, Chunk, Plane, Observation, CoordError
from caom2 import SpectralWCS, CoordAxis1D, Axis, CoordFunction1D, RefCoord
from caom2 import SpatialWCS, Dimension2D, Coord2D, CoordFunction2D
from caom2 import CoordAxis2D, PolarizationWCS, TemporalWCS
from caom2 import ObservationReader, ObservationWriter, Algorithm
from caom2 import ReleaseType, ProductType, ObservationIntentType
from caom2 import DataProductType, TargetType, Telescope, Environment
from caom2 import Instrument, Proposal, Target, Provenance, Metrics, Quality
from caom2 import CalibrationLevel, Status, Requirements, DataQuality
import logging
import sys
from six.moves.urllib.parse import urlparse
from cadcutils import net
from cadcdata import CadcDataClient
from io import BytesIO

APP_NAME = 'fits2caom2'

__all__ = ['FitsParser', 'WcsParser', 'get_cadc_headers', 'main_app']

ENERGY_CTYPES = [
    'FREQ',
    'ENER',
    'WAVN',
    'VRAD',
    'WAVE',
    'VOPT',
    'ZOPT',
    'AWAV',
    'VELO',
    'BETA']

# From http://hea-www.cfa.harvard.edu/~arots/TimeWCS/
TIME_KEYWORDS = [
    'TIME',
    'TAI',
    'TT',
    'TDT',
    'ET',
    'IAT',
    'UT1',
    'UTC',
    'GMT',
    'GPS',
    'TCG',
    'TCB',
    'TDB',
    'LOCAL']

POLARIZATION_CTYPES = ['STOKES']


class LoggingFilter(logging.Filter):
    """Add the HDU number to logging messages as a default."""

    def __init__(self):
        super(LoggingFilter, self).__init__()
        self._extension = -1

    def filter(self, record):
        record.hdu = self._extension
        return True

    def extension(self, value):
        self._extension = value


class FitsParser(object):
    """
    Parses a FITS file and extracts the CAOM2 related information which can
    be used to augment an existing CAOM2 observation, plane or artifact. If
    the FITS file is not provided, a list of dictionaries (FITS keyword=value)
    with a dictionary for each "extension" need to be provided instead.

    The WCS-related keywords of the FITS file are consumed by the astropy.wcs
    package which might display warnings with regards to compliance.

    Example 1:
    parser = FitsParser(file = '/staging/700000o.fits.gz')
    ...
    # customize parser.headers by deleting, changing or adding attributes

    obs = Observation(collection='TEST', observation_id='700000',
                      algorithm='exposure')
    plane = Plane(plane_id='700000-1')
    obs.plane.add(plane)

    artifact = Artifact(uri='ad:CFHT/700000o.fits.gz', product_type='science',
                        release_type='data')
    plane.artifacts.add(artifact)

    parser.augment_observation(obs)

    # further update obs


    Example 2:
    parser = FitsParser()

    headers = [] # list of dictionaries headers
    # populate headers

    parser.headers = headers

    parser.augment_observation(obs)
    ...

    """

    CONFIG = {'Observation.meta_release':
              ['DATE', 'DATE-OBS', 'UTCOBS', 'UTCDATE',
               'UTC-DATE', 'MJDOBS', 'MJD_OBS'],
              'Observation.instrument.name': ['INSTRUME'],
              'Observation.target.name': ['OBJECT'],
              'Observation.type': ['OBSTYPE'],
              'Observation.telescope.name': ['INSTRUME'],
              'Observation.environment.ambientTemp': ['TEMPERAT'],
              'Observation.algorithm.name': ['PROCNAME'],
              'Observation.intent': [],
              'Observation.sequenceNumber': [],
              'Observation.instrument.keywords': ['INSTMODE'],
              'Observation.proposal.id': ['RUNID'],
              'Observation.proposal.pi': [],
              'Observation.proposal.project': [],
              'Observation.proposal.title': [],
              'Observation.proposal.keywords': [],
              'Observation.target.type': [],
              'Observation.target.standard': [],
              'Observation.target.redshift': [],
              'Observation.target.keywords': [],
              'Observation.telescope.geo_location_x': ['OBSGEO-X'],
              'Observation.telescope.geo_location_y': ['OBSGEO-Y'],
              'Observation.telescope.geo_location_z': ['OBSGEO-Z'],
              'Observation.telescope.keywords': [],
              'Observation.environment.seeing': [],
              'Observation.environment.humidity': [],
              'Observation.environment.elevation': [],
              'Observation.environment.tau': [],
              'Observation.environment.wavelengthTau': [],
              'Observation.environment.photometric': [],
              'Observation.observation_id': ['OBSID'],
              'Plane.meta_release': ['RELEASE', 'REL_DATE'],
              'Plane.data_release': ['RELEASE', 'REL_DATE'],
              'Plane.dataProductType': [],
              'Plane.product_id': ['RUNID'],
              'Plane.calibrationLevel': [],
              'Plane.provenance.name': ['XPRVNAME'],
              'Plane.provenance.version': [],
              'Plane.provenance.project': ['ADC_ARCH'],
              'Plane.provenance.producer': ['ORIGIN'],
              'Plane.provenance.runID': [],
              'Plane.provenance.reference': ['XREFER'],
              'Plane.provenance.lastExecuted': ['DATE-FTS'],
              'Plane.provenance.keywords': [],
              'Plane.provenance.inputs': [],
              'Plane.metrics.sourceNumberDensity': [],
              'Plane.metrics.background': [],
              'Plane.metrics.backgroundStddev': [],
              'Plane.metrics.fluxDensityLimit': [],
              'Plane.metrics.magLimit': []
              }

    def __init__(self,
                 file=None):
        """
        Ctor
        :param file: FITS file
        """
        self.logger = logging.getLogger(__name__)
        self._headers = []
        self.parts = 0
        self.file = ''
        if file:
            self.file = file
            hdulist = fits.open(file, memmap=True, lazy_load_hdus=False)
            hdulist.close()
            self._headers = [h.header for h in hdulist]

    @property
    def headers(self):
        """
        List of headers where each header should allow dictionary like
        access to the FITS attribute in that header
        :return:
        """
        return self._headers

    @headers.setter
    def headers(self, headers):
        self._headers = headers

    def augment_artifact(self, artifact):
        """
        Augments a given CAOM2 artifact with available FITS information
        :param artifact: existing CAOM2 artifact to be augmented
        """
        assert artifact
        assert isinstance(artifact, Artifact)

        self.logger.debug(
            'Begin CAOM2 artifact augmentation for {} with {} HDUs.'.format(
                artifact.uri, len(self.headers)))

        for i, header in enumerate(self.headers):
            ii = str(i)

            # there is one Part per extension, the name is the extension number
            # Assumption:
            #    Only primary headers for 1 extension files or the extensions
            # for multiple extension files can have data and therefore
            # corresponding parts
            if (i > 0) or (len(self.headers) == 1) :
                if ii not in artifact.parts.keys():
                    artifact.parts.add(Part(ii))  # TODO use extension name?
                    self.logger.debug('Part created for HDU {}.'.format(ii))
            else:
                artifact.parts.add(Part(ii))
                self.logger.debug('Create empty part for HDU {}'.format(ii))
                continue

            part = artifact.parts[ii]

            # each Part has one Chunk
            if not part.chunks:
                part.chunks.append(Chunk())
            chunk = part.chunks[0]

            wcs_parser = WcsParser(header, self.file, i)
            wcs_parser.augment_position(chunk)
            wcs_parser.augment_energy(chunk)
            wcs_parser.augment_temporal(chunk)
            wcs_parser.augment_polarization(chunk)

        self.logger.debug(
            'End CAOM2 artifact augmentation for {}.'.format(artifact.uri))

    def augment_observation(self, observation, artifact_uri, product_id=None):
        """
        Augments a given observation with available FITS information.
        :param observation: existing CAOM2 observation to be augmented.
        :param artifact_uri:
        """
        self.logger.debug(
            'Begin CAOM2 observation augmentation for URI {}.'.format(
                artifact_uri))

        assert observation
        assert isinstance(observation, Observation)

        # TODO default value
        observation.collection = self._get_from_list('Observation.collection',
                                                     0,
                                                     'UNKNOWN')
        observation.observation_id = str(
            self._get_from_list('Observation.observation_id', 0))
        observation.algorithm = self._get_algorithm()

        # TODO default values for the following fields
        observation.sequence_number = self._get_from_list(
            'Observation.sequence_number', 0, -1)
        observation.intent = self._get_from_list('Observation.intent', 0,
                                                 ObservationIntentType.SCIENCE)
        observation.type = self._get_from_list('Observation.type', 0)
        observation.meta_release = self._get_datetime(
            self._get_from_list('Observation.meta_release', 0, datetime.now()))
        observation.requirements = self._get_requirements()
        observation.instrument = self._get_instrument()
        observation.proposal = self._get_proposal()
        observation.target = self._get_target()
        observation.target_position = self._get_target_position()
        observation.telescope = self._get_telescope()
        observation.environment = self._get_environment()

        plane = None
        if not product_id:
            product_id = self._get_from_list('Plane.product_id', index=0,
                                             default=None)
        assert product_id, 'product ID required'

        for ii in observation.planes:
            if observation.planes[ii].product_id == product_id:
                plane = observation.planes[product_id]
                break
        if plane is None:
            plane = Plane(str(product_id))
            observation.planes[product_id] = plane

        self.augment_plane(plane, artifact_uri)
        self.logger.debug(
            'End CAOM2 observation augmentation for {}.'.format(artifact_uri))

    def augment_plane(self, plane, artifact_uri):
        """
        Augments a given plane with available FITS information.
        :param plane: existing CAOM2 plane to be augmented.
        :param artifact_uri:
        """
        self.logger.debug(
            'Begin CAOM2 plane augmentation for {}.'.format(artifact_uri))

        assert plane
        assert isinstance(plane, Plane)

        plane.creator = self._get_from_list('Plane.creator_id', index=0,
                                            default='UNKNOWN')
        plane.meta_release = self._get_from_list('Plane.meta_release', index=0,
                                                 default=None)
        plane.data_release = self._get_from_list('Plane.data_release', index=0,
                                                 default=None)
        plane.data_product_type = \
            self._get_from_list('Plane.data_product_type',
                                index=0,
                                default=DataProductType.CUBE)
        plane.calibration_level = \
            self._get_from_list('Plane.calibration_level',
                                index=0,
                                default=CalibrationLevel.CALIBRATED)
        plane.provenance = self._get_provenance()
        plane.metrics = self._get_metrics()
        plane.quality = self._get_quality()

        artifact = None
        for ii in plane.artifacts:
            artifact = plane.artifacts[ii]
            if artifact.uri == artifact_uri:
                break

        if artifact is None:
            artifact = Artifact(artifact_uri, ProductType.SCIENCE,
                                ReleaseType.DATA)  # TODO
            plane.artifacts[artifact_uri] = artifact

        self.augment_artifact(artifact)

        self.logger.debug(
            'End CAOM2 plane augmentation for {}.'.format(artifact_uri))

    def _get_algorithm(self):
        """
        Create an Algorithm instance populated with available FITS information.
        :return: Algorithm
        """
        self.logger.debug('Begin CAOM2 Algorithm augmentation.')
        name = self._get_from_list('Observation.algorithm.name', index=0,
                                   default='DEFAULT')  # TODO DEFAULT VALUE
        self.logger.debug('End CAOM2 Algorithm augmentation.')
        if name:
            return Algorithm(str(name))
        else:
            return None

    def _get_instrument(self):
        """
        Create an Instrument instance populated with available FITS
        information.
        :return: Instrument
        """
        self.logger.debug('Begin CAOM2 Instrument augmentation.')
        name = self._get_from_list('Observation.instrument.name', index=0,
                                   default='UNKNOWN')  # TODO DEFAULT VALUE
        keywords = self._get_from_list('Observation.instrument.keywords',
                                       index=0, default=['UNKNOWN'])  # TODO
        self.logger.debug('End CAOM2 Instrument augmentation.')
        if name:
            instr = Instrument(str(name))
            instr.keywords.union(keywords)
            return instr
        else:
            return None

    def _get_proposal(self):
        """
        Create a Proposal instance populated with available FITS information.
        :return: Proposal
        """
        self.logger.debug('Begin CAOM2 Proposal augmentation.')
        id = self._get_from_list('Observation.proposal.id', index=0)  # TODO
        pi = self._get_from_list('Observation.proposal.pi_name',
                                 index=0)  # TODO
        project = self._get_from_list('Observation.proposal.project',
                                      index=0)  # TODO
        title = self._get_from_list('Observation.proposal.title',
                                    index=0)  # TODO
        self.logger.debug('End CAOM2 Proposal augmentation.')
        if id:
            return Proposal(str(id), pi, project, title)
        else:
            return None

    def _get_target(self):
        """
        Create a Target instance populated with available FITS information.
        :return: Target
        """
        self.logger.debug('Begin CAOM2 Target augmentation.')
        name = self._get_from_list('Observation.target.name', index=0,
                                   default='UNKNOWN')  # TODO
        target_type = self._get_from_list('Observation.target.target_type',
                                          index=0, default=TargetType.FIELD)
        standard = self._get_from_list('Observation.target.standard', index=0,
                                       default=False)  # TODO
        redshift = self._get_from_list('Observation.target.redshift', index=0)
        keywords = self._get_set_from_list('Observation.target.keywords',
                                           index=0)  # TODO
        moving = self._get_from_list('Observation.target.moving', index=0,
                                     default=False)  # TODO
        self.logger.debug('End CAOM2 Target augmentation.')
        if name:
            return Target(str(name), target_type, standard, redshift, keywords,
                          moving)
        else:
            return None

    def _get_target_position(self):
        """
        Create a Target Position instance populated with available FITS
        information.
        :return: Target Position
        """
        self.logger.debug('Begin CAOM2 TargetPosition augmentation.')
        # TODO don't know what to do here, since config file says this is
        # Chunk-level metadata
        self.logger.debug('End CAOM2 TargetPosition augmentation.')
        return None

    def _get_telescope(self):
        """
        Create a Telescope instance populated with available FITS information.
        :return: Telescope
        """
        self.logger.debug('Begin CAOM2 Telescope augmentation.')
        name = self._get_from_list('Observation.telescope.name', index=0)
        geo_x = self._get_from_list('Observation.telescope.geo_location_x',
                                    index=0)
        geo_y = self._get_from_list('Observation.telescope.geo_location_y',
                                    index=0)
        geo_z = self._get_from_list('Observation.telescope.geo_location_z',
                                    index=0)
        keywords = self._get_set_from_list('Observation.telescope.keywords',
                                           index=0)  # TODO
        if name:
            self.logger.debug('End CAOM2 Telescope augmentation.')
            return Telescope(str(name), geo_x, geo_y, geo_z, keywords)
        else:
            return None

    def _get_environment(self):
        """
        Create an Environment instance populated with available FITS
        information.
        :return: Environment
        """
        self.logger.debug('Begin CAOM2 Environment augmentation.')
        seeing = self._get_from_list('Observation.environment.seeing', index=0,
                                     default=None)  # TODO
        humidity = self._get_from_list('Observation.environment.humidity',
                                       index=0, default=None)  # TODO
        elevation = self._get_from_list('Observation.environment.elevation',
                                        index=0, default=None)  # TODO
        tau = self._get_from_list('Observation.environment.tau', index=0,
                                  default=None)  # TODO
        wavelength_tau = self._get_from_list(
            'Observation.environment.wavelengthTau', index=0,
            default=None)  # TODO
        ambient = self._get_from_list('Observation.environment.ambientTemp',
                                      index=0, default=None)  # TODO
        photometric = \
            self._get_from_list('Observation.environment.photometric',
                                index=0, default=None)  # TODO
        enviro = Environment()
        enviro.seeing = seeing
        enviro.humidity = humidity
        enviro.elevation = elevation
        enviro.tau = tau
        enviro.wavelength_tau = wavelength_tau
        enviro.ambient_temp = ambient
        enviro.photometric = photometric
        self.logger.debug('End CAOM2 Environment augmentation.')
        return enviro

    def _get_requirements(self):
        """
        Create a Requirements instance populated with available FITS
        information.
        :return: Requirements
        """
        self.logger.debug('Begin CAOM2 Requirement augmentation.')
        flag = self._get_from_list('Observation.requirements.flag', index=0,
                                   default=Status.FAIL)  # TODO DEFAULT VALUE
        self.logger.debug('End CAOM2 Requirement augmentation.')
        if flag:
            return Requirements(flag)
        else:
            return None

    def _get_from_list(self, lookup, index, default=None):
        value = default
        try:
            keywords = self.CONFIG[lookup]
        except KeyError:
            self.logger.debug(
                'Could not find \'{}\' in fits2caom2 configuration.'.format(
                    lookup))
            return value

        for ii in keywords:
            value = self.headers[index].get(ii, default)
            self.logger.debug(
                'Assigned value {} based on keyword {}'.format(value, ii))
            if value is not default:
                break

        return value

    def _get_set_from_list(self, lookup, index, default=None):
        value = default
        try:
            keywords = self.CONFIG[lookup]
        except KeyError:
            self.logger.debug(
                'Could not find \'{}\' in fits2caom2 configuration.'.format(
                    lookup))
            return value

        for ii in keywords:
            temp = self.headers[index].get(ii, default)
            self.logger.debug(
                'Assigned value {} based on keyword {}'.format(temp, ii))
            if temp is not default:
                value = set()
                for jj in temp.split(','):
                    value.add(jj)
                    print('do i get here?')
                break

        return value

    def _get_provenance(self):
        """
        Create a Provenance instance populated with available FITS information.
        :return: Provenance
        """
        self.logger.debug('Begin CAOM2 Provenance augmentation.')
        name = self._get_from_list('Plane.provenance.name', index=0)
        version = self._get_from_list('Plane.provenance.version',
                                      index=0)  # TODO DEFAULT VALUE
        project = self._get_from_list('Plane.provenance.project', index=0)
        producer = self._get_from_list('Plane.provenance.producer', index=0)
        run_id = self._get_from_list('Plane.provenance.runID', index=0)
        reference = self._get_from_list('Plane.provenance.reference', index=0)
        last_executed = self._get_datetime(
            self._get_from_list('Plane.provenance.lastExecuted',
                                index=0))  # TODO DEFAULT VALUE
        keywords = self._get_from_list('Plane.provenance.keywords', index=0,
                                       default='DEFAULT')  # TODO DEFAULT VALUE
        # inputs = self._get_from_list('Plane.provenance.inputs', index=0,
        # default=set(PlaneURI('caom:UNKNOWN/UNKNOWN/UNKNOWN')))
        #  TODO DEFAULT VALUE
        inputs = self._get_from_list('Plane.provenance.inputs', index=0,
                                     default=None)  # TODO DEFAULT VALUE
        self.logger.debug('End CAOM2 Provenance augmentation.')
        if name:
            prov = Provenance(str(name), str(version), str(project),
                              str(producer), run_id, str(reference),
                              last_executed)
            prov.keywords.union(keywords)
            if inputs:
                prov.inputs.add(inputs)
            return prov
        else:
            return None

    def _get_metrics(self):
        """
        Create a Metrics instance populated with available FITS information.
        :return: Metrics
        """
        self.logger.debug('Begin CAOM2 Metrics augmentation.')
        source_number_density = self._get_from_list(
            'Plane.metrics.sourceNumberDensity', index=0)  # TODO DEFAULT VALUE
        background = self._get_from_list('Plane.metrics.background',
                                         index=0)  # TODO DEFAULT VALUE
        background_stddev = self._get_from_list(
            'Plane.metrics.backgroundStddev', index=0)  # TODO DEFAULT VALUE
        flux_density_limit = self._get_from_list(
            'Plane.metrics.fluxDensityLimit', index=0)  # TODO DEFAULT VALUE
        mag_limit = self._get_from_list('Plane.metrics.magLimit',
                                        index=0)  # TODO DEFAULT VALUE
        metrics = Metrics()
        metrics.source_number_density = source_number_density
        metrics.background = background
        metrics.background_std_dev = background_stddev
        metrics.flux_density_limit = flux_density_limit
        metrics.mag_limit = mag_limit
        self.logger.debug('End CAOM2 Metrics augmentation.')
        return metrics

    def _get_quality(self):
        """
        Create a Quality instance populated with available FITS information.
        :return: Quality
        """
        self.logger.debug('Begin CAOM2 Quality augmentation.')
        flag = self._get_from_list('Plane.dataQuality', index=0,
                                   default=Quality.JUNK)  # TODO DEFAULT VALUE
        self.logger.debug('End CAOM2 Quality augmentation.')
        if flag:
            return DataQuality(flag)
        else:
            return None

    def _get_datetime(self, from_value):
        """
        Ensure datetime values are in MJD. Really. Just not yet.
        :param from_value:
        :return:
        """
        units = self.headers[0].get('TIMEUNIT', 'd')
        if units == 'd':
            return datetime.strptime(from_value, '%Y-%m-%d')
        else:
            return datetime(1990, 1, 1, 12, 12, 12)


class WcsParser(object):
    """
    Parser to augment chunks with positional, temporal, energy and polarization
    information based on the WCS keywords in an extension of a FITS header.

    Note: Under the hood, this class uses the astropy.wcs package to parse the
    header and any inconsistencies or missing keywords are reported back as
    warnings.
    """

    ENERGY_AXIS = 'energy'
    POLARIZATION_AXIS = 'polarization'
    TIME_AXIS = 'time'

    def __init__(self, header, file, extension):
        """

        :param header: FITS extension header
        :param file: name of FITS file
        """

        # add the HDU extension to every logging message
        self.logger = logging.getLogger(__name__ + '.fitsparser')
        self.log_filter = LoggingFilter()
        self.logger.addFilter(self.log_filter)
        self.log_filter.extension(extension)
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            '%(levelname)s:%(name)-12s:HDU:%(hdu)-2d:%(message)s')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)
        self.logger.propagate = False

        self.wcs = WCS(header)
        self.header = header
        self.file = file

    def augment_energy(self, chunk):
        """
        Augments the energy information in a chunk
        :param chunk:
        """
        assert chunk
        assert isinstance(chunk, Chunk)

        # get the energy axis
        energy_axis = self._get_axis_index(ENERGY_CTYPES)

        if energy_axis is None:
            self.logger.debug('No WCS Energy info.')
            return

        chunk.energy_axis = energy_axis
        naxis = CoordAxis1D(self._get_axis(energy_axis))

        # Note - could not avoid using _naxis private attributes...
        naxis.function = CoordFunction1D(
            self._sanitize(self.wcs._naxis[energy_axis]),
            self._sanitize(self.wcs.wcs.cdelt[energy_axis]),
            RefCoord(self._sanitize(self.wcs.wcs.crpix[energy_axis]),
                     self._sanitize(self.wcs.wcs.crval[energy_axis])))
        specsys = str(self.wcs.wcs.specsys)
        if not chunk.energy:
            chunk.energy = SpectralWCS(naxis, specsys)
        else:
            chunk.energy.naxis = naxis
            chunk.energy.specsys = specsys

        chunk.energy.ssysobs = self._sanitize(self.wcs.wcs.ssysobs)
        chunk.energy.restfrq = self._sanitize(self.wcs.wcs.restfrq)
        chunk.energy.restwav = self._sanitize(self.wcs.wcs.restwav)
        chunk.energy.velosys = self._sanitize(self.wcs.wcs.velosys)
        chunk.energy.zsource = self._sanitize(self.wcs.wcs.zsource)
        chunk.energy.ssyssrc = self._sanitize(self.wcs.wcs.ssyssrc)
        chunk.energy.velang = self._sanitize(self.wcs.wcs.velangl)

    def augment_position(self, chunk):
        """
        Augments a chunk with spatial WCS information
        :param chunk:
        :return:
        """
        self.logger.debug('Begin Spatial WCS augmentation.')

        assert chunk
        assert isinstance(chunk, Chunk)

        if self.wcs.has_celestial:
            chunk.positionAxis1, chunk.positionAxis2 = \
                self._get_position_axis()
            axis = self._get_spatial_axis(None, chunk.positionAxis1 - 1,
                                          chunk.positionAxis2 - 1)
            if not chunk.position:
                chunk.position = SpatialWCS(axis)
            else:
                chunk.position.axis = axis

            radesys = self._sanitize(self.wcs.celestial.wcs.radesys)
            chunk.position.coordsys = None if radesys is None else str(radesys)
            chunk.position.equinox = \
                self._sanitize(self.wcs.celestial.wcs.equinox)
            self.logger.debug('End Spatial WCS augmentation.')
        else:
            self.logger.warning('No celestial metadata')

    def augment_temporal(self, chunk):
        """
        Augments a chunk with temporal WCS information

        The expected caom2 - FITS keywords mapping is:

            time.exposure = EXPTIME
            time.resolution = TIMEDEL
            time.timesys = TIMESYS default UTC
            time.trefpos = TREFPOS
            time.mjdref = MJDREF | MJDDATE

        :param chunk:
        :return:
        """
        self.logger.debug('Begin TemporalWCS augmentation.')
        assert chunk
        assert isinstance(chunk, Chunk)

        time_axis = self._get_axis_index(TIME_KEYWORDS)

        if time_axis is None:
            self.logger.warning('No WCS Time info.')
            return

        chunk.timeAxis = time_axis
        # set chunk.time
        self.logger.debug('Begin temporal axis augmentation.')

        aug_naxis = self._get_axis(time_axis)
        aug_error = self._get_coord_error(None, time_axis)
        aug_ref_coord = self._get_ref_coord(None, time_axis)
        aug_function = CoordFunction1D(self.wcs._naxis[time_axis],
                                       self.wcs.wcs.cdelt[time_axis],
                                       aug_ref_coord)
        naxis = CoordAxis1D(aug_naxis, aug_error, None, None, aug_function)
        if not chunk.time:
            chunk.time = TemporalWCS(naxis)
        else:
            chunk.time.naxis = naxis

        chunk.time.exposure = self.header.get('EXPTIME')
        chunk.time.resolution = self.header.get('TIMEDEL')
        chunk.time.timesys = str(self.header.get('TIMESYS', 'UTC'))
        chunk.time.trefpos = self.header.get('TREFPOS', None)
        chunk.time.mjdref = self.header.get('MJDREF',
                                            self.header.get('MJDDATE'))
        self.logger.debug('End TemporalWCS augmentation.')

    def augment_polarization(self, chunk):
        """
        Augments a chunk with polarization WCS information
        :param chunk:
        :return:
        """
        self.logger.debug('Begin TemporalWCS augmentation.')
        assert chunk
        assert isinstance(chunk, Chunk)

        polarization_axis = self._get_axis_index(POLARIZATION_CTYPES)
        if polarization_axis is None:
            self.logger.debug('No WCS Polarization info')
            return

        chunk.polarization_axis = polarization_axis

        naxis = CoordAxis1D(Axis(str(self.wcs.wcs.ctype[polarization_axis]),
                                 str(self.wcs.wcs.cunit[polarization_axis])))
        naxis.function = CoordFunction1D(
            self._sanitize(self.wcs._naxis[polarization_axis]),
            self._sanitize(self.wcs.wcs.cdelt[polarization_axis]),
            RefCoord(self._sanitize(self.wcs.wcs.crpix[polarization_axis]),
                     self._sanitize(self.wcs.wcs.crval[polarization_axis])))
        if not chunk.polarization:
            chunk.polarization = PolarizationWCS(naxis)
        else:
            chunk.polarization.naxis = naxis

    def _get_axis_index(self, keywords):
        """
        Return the index of a specific axis type or None of it doesn't exist
        :param keywords:
        :return:
        """
        axis = None
        for i, elem in enumerate(self.wcs.axis_type_names):
            if elem in keywords:
                axis = i
                break
        return axis

    def _get_axis(self, index, over_ctype=None, over_cunit=None):
        """ Assemble a generic axis """
        aug_ctype = over_ctype if over_ctype is not None \
            else str(self.wcs.wcs.ctype[index])
        aug_cunit = over_cunit if over_cunit is not None \
            else str(self.wcs.wcs.cunit[index])
        aug_axis = Axis(aug_ctype, aug_cunit)
        return aug_axis

    def _get_spatial_axis(self, aug_axis, xindex, yindex):
        """Assemble the bits to make the axis parameter needed for
        SpatialWCS construction."""

        if aug_axis:
            raise NotImplementedError
        else:

            aug_dimension = self._get_dimension(None, xindex, yindex)

            aug_ref_coord = Coord2D(self._get_ref_coord(None, xindex),
                                    self._get_ref_coord(None, yindex))

            aug_cd11, aug_cd12, aug_cd21, aug_cd22 = \
                self._get_cd(xindex, yindex)

            aug_function = CoordFunction2D(aug_dimension, aug_ref_coord,
                                           aug_cd11, aug_cd12,
                                           aug_cd21, aug_cd22)

            aug_axis = CoordAxis2D(self._get_axis(xindex),
                                   self._get_axis(yindex),
                                   self._get_coord_error(None, xindex),
                                   self._get_coord_error(None, yindex),
                                   None, None, aug_function)
        return aug_axis

    def _get_cd(self, x_index, y_index):
        """ returns cd info"""

        try:
            if self.wcs.wcs.has_cd():
                cd11 = self.wcs.wcs.cd[x_index][x_index]
                cd12 = self.wcs.wcs.cd[x_index][y_index]
                cd21 = self.wcs.wcs.cd[y_index][x_index]
                cd22 = self.wcs.wcs.cd[y_index][y_index]
            else:
                cd11 = self.wcs.wcs.cdelt[x_index]
                cd12 = self.wcs.wcs.crota[x_index]
                cd21 = self.wcs.wcs.crota[y_index]
                cd22 = self.wcs.wcs.cdelt[y_index]
        except AttributeError:
            self.logger.warning(
                'Error searching for CD* values {}'.format(sys.exc_info()[1]))
            cd11 = -1.0  # TODO what if neither of these are defined?
            cd12 = -1.0
            cd21 = -1.0
            cd22 = -1.0

        return cd11, cd12, cd21, cd22

    def _get_coord_error(self, aug_coord_error, index, over_csyer=None,
                         over_crder=None):
        if aug_coord_error:
            raise NotImplementedError
        else:
            aug_csyer = over_csyer if over_csyer is not None \
                else self._sanitize(self.wcs.wcs.csyer[index])
            aug_crder = over_crder if over_crder is not None \
                else self._sanitize(self.wcs.wcs.crder[index])

            if aug_csyer and aug_crder:
                aug_coord_error = CoordError(aug_csyer, aug_crder)

        return aug_coord_error

    def _get_dimension(self, aug_dimension, xindex, yindex):

        if aug_dimension:
            raise NotImplementedError
        else:
            aug_dimension = Dimension2D(self.wcs._naxis[xindex],
                                        self.wcs._naxis[yindex])

        return aug_dimension

    def _get_position_axis(self):

        # there are two celestial axes, get the applicable indices from
        # the axis_types
        xindex = None
        yindex = None
        axis_types = self.wcs.get_axis_types()

        for ii in axis_types:
            if ii['coordinate_type'] == 'celestial':
                if xindex is None:
                    xindex = axis_types.index(ii)
                else:
                    yindex = axis_types.index(ii)

        # TODO determine what value to return if there is no index for an axis
        xaxis = -1 if xindex is None else int(axis_types[xindex]['number']) + 1
        yaxis = -1 if yindex is None else int(axis_types[yindex]['number']) + 1

        self.logger.debug(
            'Setting positionAxis1 to {}, positionAxis2 to {}'.format(xaxis,
                                                                      yaxis))

        return xaxis, yaxis

    def _get_ref_coord(self, aug_ref_coord, index, over_crpix=None,
                       over_crval=None):
        if aug_ref_coord:
            raise NotImplementedError
        else:
            aug_crpix = over_crpix if over_crpix is not None \
                else self.wcs.wcs.crpix[index]
            aug_crval = over_crval if over_crval is not None \
                else self.wcs.wcs.crval[index]
            aug_ref_coord = RefCoord(aug_crpix, aug_crval)
        return aug_ref_coord

    def _sanitize(self, value):
        """
        Sanitazes values from FITS to caom2
        :param value:
        :return:
        """
        if isinstance(value, float) and math.isnan(value):
            return None
        elif not str(value):
            return None  # empty string
        else:
            return value


def get_cadc_headers(uri, cert=None):
    """
    Fetches the FITS headers of a CADC file. The function takes advantage
    of the fhead feature of the CADC storage service and retrieves just the
    headers and no data, minimizing the transfer time.
    :param uri: CADC (AD like) file URI
    :param cert: X509 certificate for accessing proprietary files
    :return: List of headers corresponding to each extension. Each header is
    of astropy.wcs.Header type - essentially a dictionary of FITS keywords.
    """
    file_url = urlparse(uri)
    if file_url.scheme != 'ad':
        # TODO add hook to support other service providers
        raise NotImplementedError('Only ad type URIs supported')
    # create possible types of subjects
    subject = net.Subject(cert)
    client = CadcDataClient(subject)
    # do a fhead on the file
    archive, file_id = file_url.path.split('/')
    b = BytesIO()
    b.name = uri
    client.get_file(archive, file_id, b, fhead=True)
    fits_header = b.getvalue().decode('ascii')
    b.close()
    delim = '\nEND'
    extensions = \
        [e + delim for e in fits_header.split(delim) if e.strip()]
    headers = [fits.Header.fromstring(e, sep='\n') for e in extensions]
    return headers


def main_app():
    parser = argparse.ArgumentParser()

    parser.description = (
        'Augments an observation with information in one or more fits files.')

    if version.version is not None:
        parser.add_argument('-V', '--version', action='version',
                            version=version)

    log_group = parser.add_mutually_exclusive_group()
    log_group.add_argument('-d', '--debug', action='store_true',
                           help='debug messages')
    log_group.add_argument('-q', '--quiet', action='store_true',
                           help='run quietly')
    log_group.add_argument('-v', '--verbose', action='store_true',
                           help='verbose messages')

    parser.add_argument('--dumpconfig', action='store_true',
                        help=('output the utype to keyword mapping to '
                              'the console'))

    parser.add_argument('--ignorePartialWCS', action='store_true',
                        help='do not stop and exit upon finding partial WCS')

    parser.add_argument('-o', '--out', dest='out_obs_xml',
                        help='output of augmented observation in XML',
                        required=False)

    in_group = parser.add_mutually_exclusive_group(required=True)
    in_group.add_argument('-i', '--in', dest='in_obs_xml',
                          type=argparse.FileType('r'),
                          help='input of observation to be augmented in XML')
    in_group.add_argument('--observation', nargs=2,
                          help='observation in a collection',
                          metavar=('collection', 'observationID'))

    parser.add_argument('--config', required=False,
                        help=('optional CAOM2 utype to keyword config file to '
                              'merge with the internal configuration'))

    parser.add_argument('--default',
                        help='file with default values for keywords')
    parser.add_argument('--override',
                        help='file with override values for keywords')
    parser.add_argument('--local', nargs='+',
                        help=('list of files in local filesystem (same order '
                              'as uri)'))
    parser.add_argument('--log', help='log file name > (instead of console)')
    parser.add_argument('--keep', action='store_true',
                        help='keep the locally stored files after ingestion')
    parser.add_argument('--test', action='store_true',
                        help='test mode, do not persist to database')
    parser.add_argument('--cert', help='Proxy Cert&Key PEM file')

    parser.add_argument('productID',
                        help='product ID of the plane in the observation')
    parser.add_argument('fileURI', help='URI of a fits file', nargs='+')

    if len(sys.argv) < 2:
        # correct error message when running python3
        parser.print_usage(file=sys.stderr)
        sys.stderr.write("{}: error: too few arguments\n".format(APP_NAME))
        sys.exit(-1)

    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(level=logging.INFO, stream=sys.stdout)
    if args.debug:
        logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)
    else:
        logging.basicConfig(level=logging.WARN, stream=sys.stdout)

    if args.local and (len(args.local) != len(args.fileURI)):
        sys.stderr.write(('number of local arguments not the same with file '
                          'URIs ({} vs {})').format(len(args.local),
                                                    args.fileURI))
        sys.exit(-1)

    # invoke the appropriate function based on the inputs
    if args.in_obs_xml:
        # append to existing observation
        reader = ObservationReader(validate=True)
        obs = reader.read(args.in_obs_xml)
    else:
        obs = Observation(collection=args.observation[0],
                          observation_id=args.observation[1],
                          algorithm=Algorithm('blah'))  # TODO

    if args.productID not in obs.planes.keys():
        obs.planes.add(Plane(product_id=args.productID))

    plane = obs.planes[args.productID]
    for i, uri in enumerate(args.fileURI):
        if args.local:
            file = args.local[i]
            if uri not in plane.artifacts.keys():
                plane.artifacts.add(
                    Artifact(uri=uri,
                             product_type=ProductType.SCIENCE,
                             release_type=ReleaseType.DATA))
            artifact = plane.artifacts[uri]
            parser = FitsParser(file)
        else:
            headers = get_cadc_headers(uri, args.cert)

            if uri not in plane.artifacts.keys():
                plane.artifacts.add(
                    Artifact(uri=uri,
                             product_type=ProductType.SCIENCE,
                             release_type=ReleaseType.DATA))
            parser = FitsParser()
            parser.headers = headers

        parser.augment_observation(observation=obs, artifact_uri=uri,
                                   product_id=plane.product_id)

    writer = ObservationWriter()
    if args.out_obs_xml:
        writer.write(obs, args.out_obs_xml)
    else:
        writer.write(obs, sys.stdout)

    logging.info("DONE")

