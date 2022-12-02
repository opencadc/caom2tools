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

import argparse
from datetime import datetime
from logging.handlers import TimedRotatingFileHandler

import math
from astropy.wcs import Wcsprm, WCS
from astropy.io import fits
from astropy.time import Time
from cadcutils import version
from caom2.caom_util import int_32
from caom2 import (
    Artifact, Part, Chunk, Plane, Observation, CoordError,
    SpectralWCS, CoordAxis1D, Axis, CoordFunction1D, RefCoord,
    SpatialWCS, Dimension2D, Coord2D, CoordFunction2D,
    CoordAxis2D, CoordRange1D, PolarizationWCS, TemporalWCS,
    ObservationReader, ObservationWriter, Algorithm,
    ReleaseType, ProductType, ObservationIntentType,
    DataProductType, Telescope, Environment,
    Instrument, Proposal, Target, Provenance, Metrics,
    CalibrationLevel, Requirements, DataQuality, PlaneURI,
    SimpleObservation, DerivedObservation, ChecksumURI,
    ObservationURI, ObservableAxis, Slice, Point, TargetPosition,
    CoordRange2D, TypedSet, CustomWCS, Observable,
    CompositeObservation, EnergyTransition
)
from caom2utils import data_util
from caom2utils.caomvalidator import validate
from caom2utils.wcsvalidator import InvalidWCSError
import importlib
import logging
import os
import re
import requests
import sys
import tempfile
import traceback
from collections import defaultdict
from urllib.parse import urlparse
from cadcutils import net, util
from cadcdata import FileInfo
from vos import Client
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

APP_NAME = 'caom2gen'

__all__ = ['Caom2Exception', 'ContentParser',  'FitsParser', 'FitsWcsParser',
           'DispatchingFormatter', 'ObsBlueprint', 'get_arg_parser', 'proc',
           'POLARIZATION_CTYPES', 'gen_proc', 'get_gen_proc_arg_parser',
           'BlueprintParser', 'augment', 'get_vos_headers',
           'get_external_headers', 'Hdf5Parser', 'Hdf5ObsBlueprint',
           'Hdf5WcsParser', 'update_artifact_meta']

CUSTOM_CTYPES = [
    'RM',
    'FDEP'
]

POSITION_CTYPES = [
    ['RA',
     'GLON',
     'ELON',
     'HLON',
     'SLON'],
    ['DEC',
     'GLAT',
     'ELAT',
     'HLAT',
     'SLAT']
]

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

OBSERVABLE_CTYPES = [
    'observable',
    'FLUX']


class Caom2Exception(Exception):
    """Exception raised when an attempt to create or update a CAOM2 record
    fails for some reason."""
    pass


class HDULoggingFilter(logging.Filter):
    """Add the HDU number to logging messages as a default."""

    def __init__(self):
        super().__init__()
        self._extension = -1

    def filter(self, record):
        record.hdu = self._extension
        return True

    def extension(self, value):
        self._extension = value


class classproperty:
    """
    Class property used for CAOM2_ELEMENTS in ObsBleprint
    """
    def __init__(self, f):
        self.f = f

    def __get__(self, obj, owner):
        return self.f(owner)


class DispatchingFormatter:
    """Dispatch formatter for logger and it's sub-logger, so there can
    be multiple formatters."""

    def __init__(self, formatters, default_formatter):
        self._formatters = formatters
        self._default_formatter = default_formatter

    def format(self, record):
        logger = logging.getLogger(record.name)
        while logger:
            # check if suitable formatter for current logger exists
            if logger.name in self._formatters:
                formatter = self._formatters[logger.name]
                break
            else:
                logger = logger.parent
        else:
            # if no formatter found, just use the default
            formatter = self._default_formatter
        return formatter.format(record)


class ObsBlueprint:
    """
    Class that represents the blueprint of a CAOM2 Observation that can be
    used to build an observation.

    The following CAOM2 elements can be specified in the blueprint:
    _CAOM2_ELEMENTS

    The blueprint designates the source of each of these attributes as either
    FITS keywords with possible default values or sets the actual values.
    The blueprint can be checked by simply displaying it.

    For example:

    # display the default blueprint when WCS axes are not specified
    print(ObsBlueprint())

    # display the default blueprint when WCS axes are specified
    print(ObsBlueprint(position_axis=(1, 2), energy_axis=3,
                       polarization_axis=4, time_axis=5))

    # create a blueprint and customize it
    ob = ObsBlueprint(position_axis=(1, 2), energy_axis=3,
                      polarization_axis=4, time_axis=5))
    ob.set('Observation.algorithm.name', 'exposure')
    ob.add_attribute('Chunk.energy.axis.axis.ctype', ['MYCTYPE'],
                      extension=1)
    ob.add_attribute('Chunk.energy.axis.axis.ctype', 'MYCTYPE2',
                          extension=1)
    ob.set('Chunk.energy.velang', 33, extension=1)
    ob.set_default('Chunk.position.coordsys', 'RA-DEC', extension=1)

    ob.set('Chunk.energy.velang', 44, extension=2)
    print(ob)

    """
    _CAOM2_ELEMENTS = [
        'CompositeObservation.members',
        'DerivedObservation.members',
        'Observation.observationID',
        'Observation.type',
        'Observation.intent',
        'Observation.sequenceNumber',
        'Observation.metaRelease',
        'Observation.metaReadGroups',
        'Observation.metaProducer',
        'Observation.requirements.flag',

        'Observation.algorithm.name',

        'Observation.instrument.name',
        'Observation.instrument.keywords',

        'Observation.proposal.id',
        'Observation.proposal.pi',
        'Observation.proposal.project',
        'Observation.proposal.title',
        'Observation.proposal.keywords',

        'Observation.target.name',
        'Observation.target.type',
        'Observation.target.standard',
        'Observation.target.redshift',
        'Observation.target.keywords',
        'Observation.target.moving',
        'Observation.target.targetID',

        'Observation.target_position.point.cval1',
        'Observation.target_position.point.cval2',
        'Observation.target_position.coordsys',
        'Observation.target_position.equinox',

        'Observation.telescope.name',
        'Observation.telescope.geoLocationX',
        'Observation.telescope.geoLocationY',
        'Observation.telescope.geoLocationZ',
        'Observation.telescope.keywords',

        'Observation.environment.seeing',
        'Observation.environment.humidity',
        'Observation.environment.elevation',
        'Observation.environment.tau',
        'Observation.environment.wavelengthTau',
        'Observation.environment.ambientTemp',
        'Observation.environment.photometric',

        'Plane.productID',
        'Plane.metaRelease',
        'Plane.dataRelease',
        'Plane.dataProductType',
        'Plane.calibrationLevel',
        'Plane.dataQuality',
        'Plane.metaReadGroups',
        'Plane.dataReadGroups',
        'Plane.metaProducer',

        'Plane.provenance.name',
        'Plane.provenance.version',
        'Plane.provenance.project',
        'Plane.provenance.producer',
        'Plane.provenance.runID',
        'Plane.provenance.reference',
        'Plane.provenance.lastExecuted',
        'Plane.provenance.keywords',
        'Plane.provenance.inputs',

        'Plane.metrics.sourceNumberDensity',
        'Plane.metrics.background',
        'Plane.metrics.backgroundStddev',
        'Plane.metrics.fluxDensityLimit',
        'Plane.metrics.magLimit',
        'Plane.metrics.sampleSNR',

        'Plane.observable.ucd',

        'Artifact.productType',
        'Artifact.releaseType',
        'Artifact.contentChecksum',
        'Artifact.contentLength',
        'Artifact.contentType',
        'Artifact.contentRelease',
        'Artifact.contentReadGroups',
        'Artifact.uri',
        'Artifact.metaProducer',

        'Part.name',
        'Part.productType',
        'Part.metaProducer',

        'Chunk',
        'Chunk.naxis',
        'Chunk.observableAxis',
        'Chunk.positionAxis1',
        'Chunk.positionAxis2',
        'Chunk.energyAxis',
        'Chunk.timeAxis',
        'Chunk.polarizationAxis',
        'Chunk.metaProducer',

        'Chunk.observable.dependent.bin',
        'Chunk.observable.dependent.axis.ctype',
        'Chunk.observable.dependent.axis.cunit',
        'Chunk.observable.independent.bin',
        'Chunk.observable.independent.axis.ctype',
        'Chunk.observable.independent.axis.cunit',

        'Chunk.position.coordsys',
        'Chunk.position.equinox',
        'Chunk.position.resolution',
        'Chunk.position.axis.axis1.ctype',
        'Chunk.position.axis.axis1.cunit',
        'Chunk.position.axis.axis2.ctype',
        'Chunk.position.axis.axis2.cunit',
        'Chunk.position.axis.error1.syser',
        'Chunk.position.axis.error1.rnder',
        'Chunk.position.axis.error2.syser',
        'Chunk.position.axis.error2.rnder',
        'Chunk.position.axis.function.cd11',
        'Chunk.position.axis.function.cd12',
        'Chunk.position.axis.function.cd21',
        'Chunk.position.axis.function.cd22',
        'Chunk.position.axis.function.dimension.naxis1',
        'Chunk.position.axis.function.dimension.naxis2',
        'Chunk.position.axis.function.refCoord.coord1.pix',
        'Chunk.position.axis.function.refCoord.coord1.val',
        'Chunk.position.axis.function.refCoord.coord2.pix',
        'Chunk.position.axis.function.refCoord.coord2.val',
        'Chunk.position.axis.range.start.coord1.pix',
        'Chunk.position.axis.range.start.coord1.val',
        'Chunk.position.axis.range.start.coord2.pix',
        'Chunk.position.axis.range.start.coord2.val',
        'Chunk.position.axis.range.end.coord1.pix',
        'Chunk.position.axis.range.end.coord1.val',
        'Chunk.position.axis.range.end.coord2.pix',
        'Chunk.position.axis.range.end.coord2.val',

        'Chunk.energy.specsys',
        'Chunk.energy.ssysobs',
        'Chunk.energy.restfrq',
        'Chunk.energy.restwav',
        'Chunk.energy.velosys',
        'Chunk.energy.zsource',
        'Chunk.energy.ssyssrc',
        'Chunk.energy.velang',
        'Chunk.energy.bandpassName',
        'Chunk.energy.resolvingPower',
        'Chunk.energy.transition',
        'Chunk.energy.transition.species',
        'Chunk.energy.transition.transition',
        'Chunk.energy.axis.axis.ctype',
        'Chunk.energy.axis.axis.cunit',
        'Chunk.energy.axis.bounds.samples',
        'Chunk.energy.axis.error.syser',
        'Chunk.energy.axis.error.rnder',
        'Chunk.energy.axis.function.naxis',
        'Chunk.energy.axis.function.delta',
        'Chunk.energy.axis.function.refCoord.pix',
        'Chunk.energy.axis.function.refCoord.val',
        'Chunk.energy.axis.range.start.pix',
        'Chunk.energy.axis.range.start.val',
        'Chunk.energy.axis.range.end.pix',
        'Chunk.energy.axis.range.end.val',

        'Chunk.polarization.axis.axis.ctype',
        'Chunk.polarization.axis.axis.cunit',
        'Chunk.polarization.axis.bounds.samples',
        'Chunk.polarization.axis.error.syser',
        'Chunk.polarization.axis.error.rnder',
        'Chunk.polarization.axis.function.naxis',
        'Chunk.polarization.axis.function.delta',
        'Chunk.polarization.axis.function.refCoord.pix',
        'Chunk.polarization.axis.function.refCoord.val',
        'Chunk.polarization.axis.range.start.pix',
        'Chunk.polarization.axis.range.start.val',
        'Chunk.polarization.axis.range.end.pix',
        'Chunk.polarization.axis.range.end.val',

        'Chunk.time.exposure',
        'Chunk.time.resolution',
        'Chunk.time.timesys',
        'Chunk.time.trefpos',
        'Chunk.time.mjdref',
        'Chunk.time.axis.axis.ctype',
        'Chunk.time.axis.axis.cunit',
        'Chunk.time.axis.bounds.samples',
        'Chunk.time.axis.error.syser',
        'Chunk.time.axis.error.rnder',
        'Chunk.time.axis.function.naxis',
        'Chunk.time.axis.function.delta',
        'Chunk.time.axis.function.refCoord.pix',
        'Chunk.time.axis.function.refCoord.val',
        'Chunk.time.axis.range.start.pix',
        'Chunk.time.axis.range.start.val',
        'Chunk.time.axis.range.end.pix',
        'Chunk.time.axis.range.end.val',

        'Chunk.observable.axis.axis.ctype',
        'Chunk.observable.axis.axis.cunit',
        'Chunk.observable.axis.function.refCoord.pix',

        'Chunk.custom.axis.axis.ctype',
        'Chunk.custom.axis.axis.cunit',
        'Chunk.custom.axis.bounds.samples',
        'Chunk.custom.axis.error.syser',
        'Chunk.custom.axis.error.rnder',
        'Chunk.custom.axis.function.naxis',
        'Chunk.custom.axis.function.delta',
        'Chunk.custom.axis.function.refCoord.pix',
        'Chunk.custom.axis.function.refCoord.val',
        'Chunk.custom.axis.range.start.pix',
        'Chunk.custom.axis.range.start.val',
        'Chunk.custom.axis.range.end.pix',
        'Chunk.custom.axis.range.end.val'
        ]

    # replace _CAOM2_ELEMENTS in __doc__ with the real elements
    __doc__ = __doc__.replace('_CAOM2_ELEMENTS', '\n'.join(['\t\t{}'.format(
        elem) for elem in _CAOM2_ELEMENTS]))

    def __init__(self, position_axes=None, energy_axis=None,
                 polarization_axis=None, time_axis=None,
                 obs_axis=None, custom_axis=None, module=None,
                 update=True, instantiated_class=None):
        """
        Ctor
        :param position_axes: tuple of form (int, int) indicating the indexes
        of position axis
        :param energy_axis: index of energy axis (int)
        :param polarization_axis: index of polarization axis (int)
        :param time_axis: index of time axis (int)
        :param obs_axis: index of observable axis (int)
        :param custom_axis: index of custom axis (int)
        :param module: user-provided code, will be loaded with
            importlib.import_module if a value is provided.
        """

        if position_axes and isinstance(position_axes, tuple) and\
                (len(position_axes) != 2):
            raise ValueError(
                'Invalid position axis: {}. Must be tuple with 2 elements'.
                format(str(position_axes)))

        self.logger = logging.getLogger(__name__)

        # this is the default blueprint
        self._plan = {}
        tmp = {'Observation.metaRelease':
               (['DATE', 'DATE-OBS', 'UTCOBS', 'UTCDATE',
                 'UTC-DATE', 'MJDOBS', 'MJD_OBS'], None),
               'Observation.instrument.name': (['INSTRUME'], None),
               'Observation.type': (['OBSTYPE'], None),
               'Observation.environment.ambientTemp': (['TEMPERAT'],
                                                       None),
               # set the default for SimpleObservation construction
               'Observation.algorithm.name': (['PROCNAME'], 'exposure'),
               'Observation.instrument.keywords': (['INSTMODE'], None),
               'Observation.proposal.id': (['RUNID'], None),
               'Observation.target.name': (['OBJECT'], None),
               'Observation.telescope.name': (['TELESCOP'], None),
               'Observation.telescope.geoLocationX': (['OBSGEO-X'],
                                                      None),
               'Observation.telescope.geoLocationY': (['OBSGEO-Y'],
                                                      None),
               'Observation.telescope.geoLocationZ': (['OBSGEO-Z'],
                                                      None),
               'Observation.observationID': (['OBSID'], None),
               'Plane.calibrationLevel': ([], CalibrationLevel.RAW_STANDARD),
               'Plane.dataProductType': ([], DataProductType.IMAGE),
               'Plane.metaRelease': (['RELEASE', 'REL_DATE'], None),
               'Plane.dataRelease': (['RELEASE', 'REL_DATE'], None),
               'Plane.productID': (['RUNID'], None),
               'Plane.provenance.name': (['XPRVNAME'], None),
               'Plane.provenance.project': (['ADC_ARCH'], None),
               'Plane.provenance.producer': (['ORIGIN'], None),
               'Plane.provenance.reference': (['XREFER'], None),
               'Plane.provenance.lastExecuted': (['DATE-FTS'], None),
               'Artifact.releaseType': ([], ReleaseType.DATA),
               'Chunk': 'include'
               }
        # using the tmp to make sure that the keywords are valid
        for key in tmp:
            self.set(key, tmp[key])

        self._extensions = {}

        # contains the standard WCS keywords in the FITS file expected by the
        # astropy.WCS package.
        self._wcs_std = {
            'Chunk.naxis': 'ZNAXIS,NAXIS'
        }
        self._pos_axes_configed = False
        self._energy_axis_configed = False
        self._time_axis_configed = False
        self._polarization_axis_configed = False
        self._obs_axis_configed = False
        self._custom_axis_configed = False
        if position_axes:
            self.configure_position_axes(position_axes)

        if energy_axis:
            self.configure_energy_axis(energy_axis)

        if polarization_axis:
            self.configure_polarization_axis(polarization_axis)

        if time_axis:
            self.configure_time_axis(time_axis)

        if obs_axis:
            self.configure_observable_axis(obs_axis)

        if custom_axis:
            self.configure_custom_axis(custom_axis)

        if module:
            self._module = module
        else:
            self._module = None
        self._module_instance = instantiated_class
        # if True, existing values are used instead of defaults
        self._update = update
        # a data structure to carry around twelve bits of data at a time:
        # the first item in the set is the ctype index, and the second is
        # whether or not the index means anything, resulting in a
        # call to the blueprint configure_* methods if it's True.
        self._axis_info = {
            'custom': (0, False),
            'dec': (0, False),
            'energy': (0, False),
            'obs': (0, False),
            'polarization': (0, False),
            'ra': (0, False),
            'time': (0, False)}

    def configure_custom_axis(self, axis, override=True):
        """
        Set the expected FITS custom keywords by index in the blueprint
        and the wcs_std lookup.

        :param axis: The index expected for the custom axis.
        :param override: Set to False when reading from a file.
        :return:
        """
        if self._custom_axis_configed:
            self.logger.debug(
                'Attempt to configure already-configured custom axis.')
            return

        if override:
            self.set('Chunk.custom.axis.axis.ctype',
                     ([f'CTYPE{axis}'], None))
            self.set('Chunk.custom.axis.axis.cunit',
                     ([f'CUNIT{axis}'], None))
            self.set('Chunk.custom.axis.function.naxis',
                     ([f'NAXIS{axis}'], None))
            self.set('Chunk.custom.axis.function.delta',
                     ([f'CDELT{axis}'], None))
            self.set('Chunk.custom.axis.function.refCoord.pix',
                     ([f'CRPIX{axis}'], None))
            self.set('Chunk.custom.axis.function.refCoord.val',
                     ([f'CRVAL{axis}'], None))

        self._wcs_std['Chunk.custom.axis.axis.ctype'] = f'CTYPE{axis}'
        self._wcs_std['Chunk.custom.axis.axis.cunit'] = f'CUNIT{axis}'
        self._wcs_std['Chunk.custom.axis.function.naxis'] = f'NAXIS{axis}'
        self._wcs_std['Chunk.custom.axis.function.delta'] = f'CDELT{axis}'
        self._wcs_std['Chunk.custom.axis.function.refCoord.pix'] = \
            f'CRPIX{axis}'
        self._wcs_std['Chunk.custom.axis.function.refCoord.val'] = \
            f'CRVAL{axis}'

        self._custom_axis_configed = True

    def configure_position_axes(self, axes, override=True):
        """
        Set the expected FITS spatial keywords by indices in the blueprint and
        the wcs_std lookup.

        :param axes: The index expected for the position axes.
        :return:
        """
        if self._pos_axes_configed:
            self.logger.debug(
                'Attempt to configure already-configured position axes.')
            return

        if override:
            self.set('Chunk.position.coordsys', (['RADESYS'], None))
            self.set('Chunk.position.equinox', (['EQUINOX', 'EPOCH'], None))
            self.set('Chunk.position.axis.axis1.ctype',
                     ([f'CTYPE{axes[0]}'], None))
            self.set('Chunk.position.axis.axis1.cunit',
                     ([f'CUNIT{axes[0]}'], None))
            self.set('Chunk.position.axis.axis2.ctype',
                     ([f'CTYPE{axes[1]}'], None))
            self.set('Chunk.position.axis.axis2.cunit',
                     ([f'CUNIT{axes[1]}'], None))
            self.set('Chunk.position.axis.error1.syser',
                     ([f'CSYER{axes[0]}'], None))
            self.set('Chunk.position.axis.error1.rnder',
                     ([f'CRDER{axes[0]}'], None))
            self.set('Chunk.position.axis.error2.syser',
                     ([f'CSYER{axes[1]}'], None))
            self.set('Chunk.position.axis.error2.rnder',
                     ([f'CRDER{axes[1]}'], None))
            self.set('Chunk.position.axis.function.cd11',
                     ([f'CD{axes[0]}_{axes[0]}'], None))
            self.set('Chunk.position.axis.function.cd12',
                     ([f'CD{axes[0]}_{axes[1]}'], None))
            self.set('Chunk.position.axis.function.cd21',
                     ([f'CD{axes[1]}_{axes[0]}'], None))
            self.set('Chunk.position.axis.function.cd22',
                     ([f'CD{axes[1]}_{axes[1]}'], None))
            self.set('Chunk.position.axis.function.dimension.naxis1',
                     ([f'ZNAXIS{axes[0]}',
                       f'NAXIS{axes[0]}'], None))
            self.set('Chunk.position.axis.function.dimension.naxis2',
                     ([f'ZNAXIS{axes[1]}',
                       f'NAXIS{axes[1]}'], None))
            self.set('Chunk.position.axis.function.refCoord.coord1.pix',
                     ([f'CRPIX{axes[0]}'], None))
            self.set('Chunk.position.axis.function.refCoord.coord1.val',
                     ([f'CRVAL{axes[0]}'], None))
            self.set('Chunk.position.axis.function.refCoord.coord2.pix',
                     ([f'CRPIX{axes[1]}'], None))
            self.set('Chunk.position.axis.function.refCoord.coord2.val',
                     ([f'CRVAL{axes[1]}'], None))

        self._wcs_std['Chunk.position.coordsys'] = 'RADESYS'
        self._wcs_std['Chunk.position.equinox'] = 'EQUINOX'

        self._wcs_std['Chunk.position.axis.axis1.ctype'] = \
            f'CTYPE{axes[0]}'
        self._wcs_std['Chunk.position.axis.axis1.cunit'] = \
            f'CUNIT{axes[0]}'
        self._wcs_std['Chunk.position.axis.axis2.ctype'] = \
            f'CTYPE{axes[1]}'
        self._wcs_std['Chunk.position.axis.axis2.cunit'] = \
            f'CUNIT{axes[1]}'
        self._wcs_std['Chunk.position.axis.error1.syser'] = \
            f'CSYER{axes[0]}'
        self._wcs_std['Chunk.position.axis.error1.rnder'] = \
            f'CRDER{axes[0]}'
        self._wcs_std['Chunk.position.axis.error2.syser'] = \
            f'CSYER{axes[1]}'
        self._wcs_std['Chunk.position.axis.error2.rnder'] = \
            f'CRDER{axes[1]}'
        self._wcs_std['Chunk.position.axis.function.cd11'] = \
            f'CD{axes[0]}_{axes[0]}'
        self._wcs_std['Chunk.position.axis.function.cd12'] = \
            f'CD{axes[0]}_{axes[1]}'
        self._wcs_std['Chunk.position.axis.function.cd21'] = \
            f'CD{axes[1]}_{axes[0]}'
        self._wcs_std['Chunk.position.axis.function.cd22'] = \
            f'CD{axes[1]}_{axes[1]}'
        self._wcs_std['Chunk.position.axis.function.dimension.naxis1'] = \
            f'NAXIS{axes[0]}'
        self._wcs_std['Chunk.position.axis.function.dimension.naxis2'] = \
            f'NAXIS{axes[1]}'
        self._wcs_std['Chunk.position.axis.function.refCoord.coord1.pix'] \
            = f'CRPIX{axes[0]}'
        self._wcs_std['Chunk.position.axis.function.refCoord.coord1.val'] \
            = f'CRVAL{axes[0]}'
        self._wcs_std['Chunk.position.axis.function.refCoord.coord2.pix'] \
            = f'CRPIX{axes[1]}'
        self._wcs_std['Chunk.position.axis.function.refCoord.coord2.val'] \
            = f'CRVAL{axes[1]}'

        self._pos_axes_configed = True

    def configure_energy_axis(self, axis, override=True):
        """
        Set the expected FITS energy keywords by index in the blueprint and
        the wcs_std lookup.

        :param axis: The index expected for the energy axis.
        :param override: Set to False when reading from a file.
        :return:
        """
        if self._energy_axis_configed:
            self.logger.debug(
                'Attempt to configure already-configured energy axis.')
            return

        if override:
            self.set('Chunk.energy.specsys', (['SPECSYS'], None))
            self.set('Chunk.energy.ssysobs', (['SSYSOBS'], None))
            self.set('Chunk.energy.restfrq', (['RESTFRQ'], None))
            self.set('Chunk.energy.restwav', (['RESTWAV'], None))
            self.set('Chunk.energy.velosys', (['VELOSYS'], None))
            self.set('Chunk.energy.zsource', (['ZSOURCE'], None))
            self.set('Chunk.energy.ssyssrc', (['SSYSSRC'], None))
            self.set('Chunk.energy.velang', (['VELANG'], None))

            self.set('Chunk.energy.bandpassName', ([], None))
            self.set('Chunk.energy.resolvingPower', ([], None))

            self.set('Chunk.energy.axis.axis.ctype',
                     ([f'CTYPE{axis}'], None))
            self.set('Chunk.energy.axis.axis.cunit',
                     ([f'CUNIT{axis}'], None))
            self.set('Chunk.energy.axis.error.syser',
                     ([f'CSYER{axis}'], None))
            self.set('Chunk.energy.axis.error.rnder',
                     ([f'CRDER{axis}'], None))
            self.set('Chunk.energy.axis.function.naxis',
                     ([f'NAXIS{axis}'], None))
            self.set('Chunk.energy.axis.function.delta',
                     ([f'CDELT{axis}'], None))
            self.set('Chunk.energy.axis.function.refCoord.pix',
                     ([f'CRPIX{axis}'], None))
            self.set('Chunk.energy.axis.function.refCoord.val',
                     ([f'CRVAL{axis}'], None))

        self._wcs_std['Chunk.energy.specsys'] = 'SPECSYS'
        self._wcs_std['Chunk.energy.ssysobs'] = 'SSYSOBS'
        self._wcs_std['Chunk.energy.restfrq'] = 'RESTFRQ'
        self._wcs_std['Chunk.energy.restwav'] = 'RESTWAV'
        self._wcs_std['Chunk.energy.velosys'] = 'VELOSYS'
        self._wcs_std['Chunk.energy.zsource'] = 'ZSOURCE'
        self._wcs_std['Chunk.energy.ssyssrc'] = 'SSYSSRC'
        self._wcs_std['Chunk.energy.velang'] = 'VELANG'

        self._wcs_std['Chunk.energy.axis.axis.ctype'] = \
            f'CTYPE{axis}'
        self._wcs_std['Chunk.energy.axis.axis.cunit'] = \
            f'CUNIT{axis}'
        self._wcs_std['Chunk.energy.axis.error.syser'] = \
            f'CSYER{axis}'
        self._wcs_std['Chunk.energy.axis.error.rnder'] = \
            f'CRDER{axis}'
        self._wcs_std['Chunk.energy.axis.function.naxis'] = \
            f'NAXIS{axis}'
        self._wcs_std['Chunk.energy.axis.function.delta'] = \
            f'CDELT{axis}'
        self._wcs_std['Chunk.energy.axis.function.refCoord.pix'] = \
            f'CRPIX{axis}'
        self._wcs_std['Chunk.energy.axis.function.refCoord.val'] = \
            f'CRVAL{axis}'

        self._energy_axis_configed = True

    def configure_polarization_axis(self, axis, override=True):
        """
        Set the expected FITS polarization keywords by index in the blueprint
        and the wcs_std lookup.

        :param axis: The index expected for the polarization axis.
        :param override: Set to False when reading from a file.
        :return:
        """
        if self._polarization_axis_configed:
            self.logger.debug(
                'Attempt to configure already-configured polarization axis.')
            return

        if override:
            self.set('Chunk.polarization.axis.axis.ctype',
                     ([f'CTYPE{axis}'], None))
            self.set('Chunk.polarization.axis.axis.cunit',
                     ([f'CUNIT{axis}'], None))
            self.set('Chunk.polarization.axis.function.naxis',
                     ([f'NAXIS{axis}'], None))
            self.set('Chunk.polarization.axis.function.delta',
                     ([f'CDELT{axis}'], None))
            self.set('Chunk.polarization.axis.function.refCoord.pix',
                     ([f'CRPIX{axis}'], None))
            self.set('Chunk.polarization.axis.function.refCoord.val',
                     ([f'CRVAL{axis}'], None))

        self._wcs_std['Chunk.polarization.axis.axis.ctype'] = \
            f'CTYPE{axis}'
        self._wcs_std['Chunk.polarization.axis.axis.cunit'] = \
            f'CUNIT{axis}'
        self._wcs_std['Chunk.polarization.axis.function.naxis'] = \
            f'NAXIS{axis}'
        self._wcs_std['Chunk.polarization.axis.function.delta'] = \
            f'CDELT{axis}'
        self._wcs_std['Chunk.polarization.axis.function.refCoord.pix'] = \
            f'CRPIX{axis}'
        self._wcs_std['Chunk.polarization.axis.function.refCoord.val'] = \
            f'CRVAL{axis}'

        self._polarization_axis_configed = True

    def configure_observable_axis(self, axis, override=True):
        """
        Set the expected FITS observable keywords by index in the blueprint
        and the wcs_std lookup.
        Note: observable axis is not a standard WCS and it's not used by
        astropy.wcs so, arguably, it can be removed. It is here for now for
        consistency purposes.
        :param axis: The index expected for the observable axis.
        :param override: Set to False when reading from a file.
        :return:
        """
        if self._obs_axis_configed:
            self.logger.debug(
                'Attempt to configure already-configured observable axis.')
            return

        if override:
            self.set('Chunk.observable.axis.axis.ctype',
                     ([f'CTYPE{axis}'], None))
            self.set('Chunk.observable.axis.axis.cunit',
                     ([f'CUNIT{axis}'], None))
            self.set('Chunk.observable.axis.function.refCoord.pix',
                     ([f'CRPIX{axis}'], None))

        self._wcs_std['Chunk.observable.axis.axis.ctype'] = \
            f'CTYPE{axis}'
        self._wcs_std['Chunk.observable.axis.axis.cunit'] = \
            f'CUNIT{axis}'
        self._wcs_std['Chunk.observable.axis.function.refCoord.pix'] = \
            f'CRPIX{axis}'

        self._obs_axis_configed = True

    def configure_time_axis(self, axis, override=True):
        """
        Set the expected FITS time keywords by index in the blueprint and
        the wcs_std lookup.

        :param axis: The index expected for the time axis.
        :param override: Set to False when reading from a file.
        :return:
        """
        if self._time_axis_configed:
            self.logger.debug(
                'Attempt to configure already-configured time axis.')
            return

        if override:
            self.set('Chunk.time.exposure', (['EXPTIME', 'INTTIME'], None))
            self.set('Chunk.time.timesys', (['TIMESYS'], None))
            self.set('Chunk.time.trefpos', (['TREFPOS'], None))
            self.set('Chunk.time.mjdref', (['MJDREF'], None))
            self.set('Chunk.time.resolution', (['TIMEDEL'], None))
            self.set('Chunk.time.axis.axis.ctype',
                     ([f'CTYPE{axis}'], None))
            self.set('Chunk.time.axis.axis.cunit',
                     ([f'CUNIT{axis}'], None))
            self.set('Chunk.time.axis.error.syser',
                     ([f'CSYER{axis}'], None))
            self.set('Chunk.time.axis.error.rnder',
                     ([f'CRDER{axis}'], None))
            self.set('Chunk.time.axis.function.naxis',
                     ([f'NAXIS{axis}'], None))
            self.set('Chunk.time.axis.function.delta',
                     ([f'CDELT{axis}'], None))
            self.set('Chunk.time.axis.function.refCoord.pix',
                     ([f'CRPIX{axis}'], None))
            self.set('Chunk.time.axis.function.refCoord.val',
                     ([f'CRVAL{axis}'], None))

        self._wcs_std['Chunk.time.exposure'] = 'EXPTIME'
        self._wcs_std['Chunk.time.resolution'] = 'TIMEDEL'
        self._wcs_std['Chunk.time.timesys'] = 'TIMESYS'
        self._wcs_std['Chunk.time.trefpos'] = 'TREFPOS'
        self._wcs_std['Chunk.time.mjdref'] = 'MJDREF'

        self._wcs_std['Chunk.time.axis.axis.ctype'] = \
            f'CTYPE{axis}'
        self._wcs_std['Chunk.time.axis.axis.cunit'] = \
            f'CUNIT{axis}'
        self._wcs_std['Chunk.time.axis.error.syser'] = \
            f'CSYER{axis}'
        self._wcs_std['Chunk.time.axis.error.rnder'] = \
            f'CRDER{axis}'
        self._wcs_std['Chunk.time.axis.function.naxis'] = \
            f'NAXIS{axis}'
        self._wcs_std['Chunk.time.axis.function.delta'] = \
            f'CDELT{axis}'
        self._wcs_std['Chunk.time.axis.function.refCoord.pix'] = \
            f'CRPIX{axis}'
        self._wcs_std['Chunk.time.axis.function.refCoord.val'] = \
            f'CRVAL{axis}'

        self._time_axis_configed = True

    def _guess_axis_info(self):
        """Look for info regarding axis types in the blueprint wcs_std.
        Configure the blueprint according to the guesses.
        """
        for ii in self._plan:
            if isinstance(self._plan[ii], tuple):
                for value in self._plan[ii][0]:
                    if (value.startswith('CTYPE')) and value[-1].isdigit():
                        value = value.split('-')[0]
                        self._guess_axis_info_from_ctypes(ii, int(value[-1]))
            else:
                value = self._plan[ii]
                if value is None:
                    continue
                if (value.startswith('CTYPE')) and value[-1].isdigit():
                    value = value.split('-')[0]
                    self._guess_axis_info_from_ctypes(ii, int(value[-1]))

        self._guess_axis_info_from_plan()

    def _guess_axis_info_from_plan(self):
        for ii in self._plan:
            if ii.startswith('Chunk.position') and ii.endswith('axis1.ctype') \
                    and not self._axis_info['ra'][1]:
                configured_index = self._get_configured_index(
                    self._axis_info, 'ra')
                self._axis_info['ra'] = (configured_index, True)
            elif ii.startswith('Chunk.position') and \
                    ii.endswith('axis2.ctype') and not \
                    self._axis_info['dec'][1]:
                configured_index = self._get_configured_index(self._axis_info,
                                                              'dec')
                self._axis_info['dec'] = (configured_index, True)
            elif ii.startswith('Chunk.energy') and not \
                    self._axis_info['energy'][1]:
                configured_index = self._get_configured_index(self._axis_info,
                                                              'energy')
                self._axis_info['energy'] = (configured_index, True)
            elif ii.startswith('Chunk.time') and not \
                    self._axis_info['time'][1]:
                configured_index = self._get_configured_index(self._axis_info,
                                                              'time')
                self._axis_info['time'] = (configured_index, True)
            elif ii.startswith('Chunk.polarization') \
                    and not self._axis_info['polarization'][1]:
                configured_index = self._get_configured_index(self._axis_info,
                                                              'polarization')
                self._axis_info['polarization'] = (configured_index, True)
            elif ii.startswith('Chunk.observable') and not \
                    self._axis_info['obs'][1]:
                configured_index = self._get_configured_index(self._axis_info,
                                                              'obs')
                self._axis_info['obs'] = (configured_index, True)
            elif ii.startswith('Chunk.custom') and not \
                    self._axis_info['custom'][1]:
                configured_index = self._get_configured_index(self._axis_info,
                                                              'custom')
                self._axis_info['custom'] = (configured_index, True)

        if self._axis_info['ra'][1] and self._axis_info['dec'][1]:
            self.configure_position_axes(
                (self._axis_info['ra'][0], self._axis_info['dec'][0]), False)
        elif self._axis_info['ra'][1] or self._axis_info['dec'][1]:
            raise ValueError('Only one positional axis found '
                             '(ra/dec): {}/{}'.
                             format(self._axis_info['ra'][0],
                                    self._axis_info['dec'][0]))
        else:
            # assume that positional axis are 1 and 2 by default
            if (self._axis_info['time'][0] in [1, 2] or
                self._axis_info['energy'][0] in [1, 2] or
                self._axis_info['polarization'][0] in [1, 2] or
                self._axis_info['obs'][0] in [1, 2] or
                    self._axis_info['custom'][0] in [1, 2]):
                raise ValueError('Cannot determine the positional axis')
            else:
                self.configure_position_axes((1, 2), False)

        if self._axis_info['time'][1]:
            self.configure_time_axis(self._axis_info['time'][0], False)
        if self._axis_info['energy'][1]:
            self.configure_energy_axis(self._axis_info['energy'][0], False)
        if self._axis_info['polarization'][1]:
            self.configure_polarization_axis(
                self._axis_info['polarization'][0], False)
        if self._axis_info['obs'][1]:
            self.configure_observable_axis(self._axis_info['obs'][0], False)
        if self._axis_info['custom'][1]:
            self.configure_custom_axis(self._axis_info['custom'][0], False)

    def _guess_axis_info_from_ctypes(self, lookup, counter):
        """
        Check for the presence of blueprint keys in the plan, and whether or
        not they indicate an index in their configuration.

        :param lookup: Blueprint plan key.
        :param counter: Value to set the index to for an axis.
        :param axis_info: local data structure to pass around what is
            configured, and what is it's value.
        """
        if lookup.startswith('Chunk.energy'):
            self._axis_info['energy'] = (counter, True)
        elif lookup.startswith('Chunk.polarization'):
            self._axis_info['polarization'] = (counter, True)
        elif lookup.startswith('Chunk.time'):
            self._axis_info['time'] = (counter, True)
        elif lookup.startswith('Chunk.position') and lookup.endswith(
                'axis1.ctype'):
            self._axis_info['ra'] = (counter, True)
        elif lookup.startswith('Chunk.position') and lookup.endswith(
                'axis2.ctype'):
            self._axis_info['dec'] = (counter, True)
        elif lookup.startswith('Chunk.observable'):
            self._axis_info['obs'] = (counter, True)
        elif lookup.startswith('Chunk.custom'):
            self._axis_info['custom'] = (counter, True)
        else:
            raise ValueError(
                f'Unrecognized axis type: {lookup}')

    def _get_configured_index(self, axis_info, lookup):
        """Find the next available index value among those that are not set.

        :param axis_info: local data structure to pass around what is
            configured, and what is it's value."""
        DEFAULT_INDICES = {'ra': 1,
                           'dec': 2,
                           'energy': 3,
                           'time': 4,
                           'polarization': 5,
                           'obs': 6,
                           'custom': 7}

        # the logic - if the default index is already used, assign the lowest
        # index that is unused, otherwise use the default index

        max_index = 0
        min_index = 7
        default_index = DEFAULT_INDICES[lookup]
        default_used = False
        for axis in axis_info:
            # do two unrelated things in this for loop
            # 1. determine where to start counting
            if axis_info[axis][1]:
                max_index = max(max_index, axis_info[axis][0])
                min_index = min(min_index, axis_info[axis][0])
            # 2. determine if the default is used
            if axis_info[axis][1] and default_index == axis_info[axis][0]:
                default_used = True

        configured_index = 0
        if default_used:
            if min_index == 1:
                configured_index = max_index + 1
            else:
                configured_index = min(1, min_index)
        else:
            configured_index = default_index
        return configured_index

    def load_from_file(self, file_name):
        """
        Load a blueprint from a file. The expected input format is the same
        as is output by _serialize. This means there's lots of stripping of
        extra spaces, equals signs, and the word default. Also manage
        square brackets as list construction.

        Accept comments that start with '#'.

        :param file_name: The fully-qualified pathname for the blueprint
        file on disk.
        """
        with open(file_name) as file:
            for line in file:
                if '=' in line:
                    if '#' in line:
                        if line.find('#') == 0:
                            # ignore lines starting with a comment
                            continue
                        line = line.split('#')[0]
                    key, value = line.split('=', 1)
                    if 'default' in value:
                        temp = value.replace('default', ''). \
                            replace('=', '').strip('\n').strip()
                        default = temp.rsplit(',')[1]
                        temp_list = temp.rsplit(',')[0].replace('[', ''). \
                            replace(']', '').replace('\'', '').split(',')
                        if 'None' in default:
                            default = None
                        else:
                            default = default.strip()
                        cleaned_up_value = (temp_list, default)
                    else:
                        if '[' in value:
                            temp_list = value.replace('[', ''). \
                                replace(']', '').replace('\'', '').split(',')
                            temp_list_2 = []
                            for ii in temp_list:
                                temp_list_2.append(ii.strip().strip('\n'))
                            cleaned_up_value = (temp_list_2, None)
                        else:
                            cleaned_up_value = value.strip('\n').strip()
                            if cleaned_up_value == 'None':
                                cleaned_up_value = None
                    self.set(key.strip(), cleaned_up_value)
        self._guess_axis_info()

    @classproperty
    def CAOM2_ELEMENTS(cls):
        """
        List of valid names of CAOM2 elements.
        :return:
        """
        return list(ObsBlueprint._CAOM2_ELEMENTS)  # return a copy

    @classmethod
    def check_caom2_element(cls, caom2_element):
        """
        Checks that an element is a valid caom2_element in the blueprint. It
        checks that it's part of the ObsBlueprint._CAOM2_ELEMENTS
        :param caom2_element: name CAOM2 element to check
        :raises KeyError
        """
        if caom2_element not in cls._CAOM2_ELEMENTS:
            raise KeyError(
                '{} not a valid CAOM2 element name (mispelling?).'.
                format(caom2_element))

    @staticmethod
    def check_chunk(caom2_element):
        """
        Checks that an element is a valid Chunk-type caom2_element
        :param caom2_element: name CAOM2 element to check
        :raises ValueError
        """
        if not caom2_element.startswith('Chunk'):
            raise ValueError(
                "Extension number refers to Chunk elements only")

    @staticmethod
    def check_extension(extension):
        if extension is not None and extension < 0:
            raise ValueError(
                f'Extension count failure. {extension} should be >= 0')

    def __str__(self):
        plan = self._serialize(self._plan)

        extensions = ''
        if self._extensions:
            for key in sorted(self._extensions):
                extensions = extensions + f'\nextension {key}:\n' +\
                    self._serialize(self._extensions[key])
        return plan + extensions

    def _serialize(self, src):
        return '\n'.join(
            ['{} = {}'.format(key, '{}, default = {}'.format(src[key][0],
                                                             src[key][1])
             if isinstance(src[key], tuple)
             else src[key])
             for key in ObsBlueprint._CAOM2_ELEMENTS
             if key in src])

    def set(self, caom2_element, value, extension=0):
        """
        Sets the value associated with an element in the CAOM2 model. Value
        cannot be a tuple.
        :param caom2_element: name CAOM2 element (as in
        ObsBlueprint.CAOM2_ELEMEMTS)
        :param value: new value of the CAOM2 element
        :param extension: extension number (used only for Chunk elements)
        """
        ObsBlueprint.check_caom2_element(caom2_element)
        ObsBlueprint.check_extension(extension)
        if extension:
            ObsBlueprint.check_chunk(caom2_element)
            if extension not in self._extensions:
                self._extensions[extension] = {}
            self._extensions[extension][caom2_element] = value
        else:
            self._plan[caom2_element] = value

    def add_attribute(self, caom2_element, attribute, extension=0):
        """
        Adds an attribute in the list of other attributes associated
        with an caom2 element.
        :param caom2_element: name CAOM2 element (as in
        ObsBlueprint.CAOM2_ELEMEMTS)
        :param attribute: name of attribute the element is mapped to
        :param extension: extension number (used only for Chunk elements)
        :raises AttributeError if the caom2 element has already an associated
        value or KeyError if the caom2 element does not exists.
        """
        ObsBlueprint.check_caom2_element(caom2_element)
        ObsBlueprint.check_extension(extension)
        if extension:
            ObsBlueprint.check_chunk(caom2_element)
            if extension not in self._extensions:
                raise AttributeError(
                    f'No extension {extension} in the blueprint')
            else:
                if caom2_element in self._extensions[extension]:
                    if (isinstance(self._extensions[extension][caom2_element],
                                   tuple)):
                        if (attribute not in
                                self._extensions[extension][caom2_element][0]):
                            self._extensions[extension][caom2_element][0].\
                                insert(0, attribute)
                    else:
                        raise AttributeError(
                            (f'No attributes in extension {extension} '
                             f'associated with keyword {caom2_element}'))
                else:
                    self._extensions[extension][caom2_element] = \
                        ([attribute], None)
        else:
            if caom2_element in self._plan:
                if isinstance(self._plan[caom2_element], tuple):
                    if attribute not in self._plan[caom2_element][0]:
                        self._plan[caom2_element][0].insert(0, attribute)
                else:
                    raise AttributeError(f'No attributes associated with '
                                         f'keyword {caom2_element}')
            else:
                self._plan[caom2_element] = ([attribute], None)

    def add_table_attribute(self, caom2_element, ttype_attribute, extension=0,
                            index=0):
        """
        Adds a FITS BINTABLE TTYPE* lookup, to a list of other FITS attributes
        associated with an caom2 element. This does not co-exist with
        non-table attributes.

        There is no support for default values for table attributes.

        :param caom2_element: name CAOM2 element (as in
        ObsBlueprint.CAOM2_ELEMEMTS)
        :param ttype_attribute: name of TTYPE attribute element is mapped to
        :param extension: extension number (used only for Chunk elements)
        :param index: which row values to return. If index is None, all row
            values will be returned as a comma-separated list.
        :raises AttributeError if the caom2 element has already an associated
        value or KeyError if the caom2 element does not exists.
        """
        ObsBlueprint.check_caom2_element(caom2_element)
        ObsBlueprint.check_extension(extension)
        if extension:
            if extension in self._extensions:
                if caom2_element in self._extensions[extension]:
                    if (ObsBlueprint.is_table(
                            self._extensions[extension][caom2_element])):
                        if (ttype_attribute not in
                                self._extensions[extension][caom2_element][1]):
                            self._extensions[extension][caom2_element][1]. \
                                insert(0, ttype_attribute)
                    else:
                        raise AttributeError(
                            ('No TTYPE attributes in extension {} associated '
                             'with keyword {}').format(extension,
                                                       caom2_element))
                else:
                    self._extensions[extension][caom2_element] = \
                        ('BINTABLE', [ttype_attribute], index)
            else:
                self._extensions[extension] = {}
                self._extensions[extension][caom2_element] = \
                    ('BINTABLE', [ttype_attribute], index)
        else:
            if caom2_element in self._plan:
                if ObsBlueprint.is_table(self._plan[caom2_element]):
                    if ttype_attribute not in self._plan[caom2_element][1]:
                        self._plan[caom2_element][1].insert(0, ttype_attribute)
                else:
                    raise AttributeError('No TTYPE attributes associated '
                                         'with keyword {}'.format(
                                            caom2_element))
            else:
                self._plan[caom2_element] = (
                    'BINTABLE', [ttype_attribute], None)

    def set_default(self, caom2_element, default, extension=0):
        """
        Sets the default value of a caom2 element that is associated with
        attributes. If the element does not exist or does not have a list of
        associated attributes, default is set as the associated value
        of the element.

        If set is called for the same caom2_element after this, the default
        value will be reset to None.

        :param caom2_element: name CAOM2 element (as in
        ObsBlueprint.CAOM2_ELEMEMTS)
        :param default: default value
        :param extension: extension number (used only for Chunk elements)
        """
        ObsBlueprint.check_caom2_element(caom2_element)
        ObsBlueprint.check_extension(extension)
        if extension:
            ObsBlueprint.check_chunk(caom2_element)
            if extension not in self._extensions:
                self._extensions[extension] = {}
            if caom2_element in self._extensions[extension] and \
               isinstance(self._extensions[extension][caom2_element], tuple):
                self._extensions[extension][caom2_element] = \
                    (self._extensions[extension][caom2_element][0], default)
            else:
                # default is the only value
                self._extensions[extension][caom2_element] = default
        else:
            if (caom2_element in self._plan) and \
                    isinstance(self._plan[caom2_element], tuple):
                self._plan[caom2_element] = (self._plan[caom2_element][0],
                                             default)
            else:
                # override the value
                self._plan[caom2_element] = default

    def delete(self, caom2_element, extension=0):
        """
        Deletes an element from the blueprint
        :param caom2_element: name CAOM2 element (as in
        ObsBlueprint.CAOM2_ELEMEMTS)
        :param extension: extension number
        :raises exceptions if the element or extension not found
        """
        ObsBlueprint.check_caom2_element(caom2_element)
        ObsBlueprint.check_extension(extension)
        if extension:
            ObsBlueprint.check_chunk(caom2_element)
            if extension not in self._extensions:
                raise ValueError('Extension {} not configured in blueprint'.
                                 format(extension))
            if caom2_element in self._extensions[extension]:
                del self._extensions[extension][caom2_element]
            if len(self._extensions[extension]) == 0:
                del self._extensions[extension]
        else:
            if caom2_element in self._plan:
                del self._plan[caom2_element]

    def clear(self, caom2_element, extension=0):
        """
        Clears the value for an element in the blueprint by resetting it to an
        empty list with no default.

        :param caom2_element: name CAOM2 element (as in
        ObsBlueprint.CAOM2_ELEMEMTS)
        :param extension: extension number
        :raises exceptions if the element or extension not found
        """
        ObsBlueprint.check_caom2_element(caom2_element)
        ObsBlueprint.check_extension(extension)
        if extension:
            ObsBlueprint.check_chunk(caom2_element)
            if extension not in self._extensions:
                raise ValueError('Extension {} not configured in blueprint'.
                                 format(extension))
            if caom2_element in self._extensions[extension]:
                self._extensions[extension][caom2_element] = ([], None)
        else:
            if caom2_element in self._plan:
                self._plan[caom2_element] = ([], None)

    def _get(self, caom2_element, extension=0):
        """
        Returns the source associated with a CAOM2 element
        :param caom2_element: name CAOM2 element (as in
        ObsBlueprint.CAOM2_ELEMEMTS)
        :param extension: extension number
        :return: Tuple of the form (list_of_associated_attributes,
        default_value) OR the actual value associated with the CAOM2 element
        """
        ObsBlueprint.check_caom2_element(caom2_element)
        ObsBlueprint.check_extension(extension)
        if extension:
            if (extension in self._extensions) and \
                    (caom2_element in self._extensions[extension]):
                return self._extensions[extension][caom2_element]

        # look in the minimal plan
        if caom2_element not in self._plan:
            return None
        else:
            return self._plan[caom2_element]

    def has_chunk(self, extension):
        """What does the plan say about creating chunks for an
        extension?

        :return True if there should be a chunk to go along with a part
        """
        value = ''
        if extension is not None and extension in self._extensions:
            if 'Chunk' in self._extensions[extension]:
                value = self._extensions[extension]['Chunk']
        elif 'Chunk' in self._plan:
            if ((extension is not None and extension == 0) or (
                    extension is None)):
                value = self._plan['Chunk']
        return not value == '{ignore}'

    @staticmethod
    def is_table(value):
        """Hide the blueprint structure from clients - they shouldn't need
        to know that a value of type tuple requires special processing."""
        return ObsBlueprint.needs_lookup(value) and value[0] == 'BINTABLE'

    @staticmethod
    def is_function(value):
        """
        Check if a blueprint value has Python 'function' syntax. The
        "'/' not in value" clause excludes strings with syntax that enables
        addressing HDF5 arrays.

        :return: True if the value is the name of a function to be executed,
            False, otherwise
        """
        return (not ObsBlueprint.needs_lookup(value) and isinstance(value, str)
                and isinstance(value, str) and '(' in value and ')' in value
                and '/' not in value)

    @staticmethod
    def has_default_value(value):
        """"""
        return isinstance(value, tuple) and value[1]

    @staticmethod
    def has_no_value(value):
        """If functions return None, try not to update the WCS with this
        value."""
        return value is None or (
                isinstance(value, str) and 'None' in value.strip())

    @staticmethod
    def needs_lookup(value):
        """Hide the blueprint structure from clients - they shouldn't need
        to know that a value of type tuple requires special processing."""
        return isinstance(value, tuple)

    def get_configed_axes_count(self):
        """:return how many axes have been configured to read from WCS"""
        configed_axes = 0
        if self._pos_axes_configed:
            configed_axes += 2
        if self._energy_axis_configed:
            configed_axes += 1
        if self._time_axis_configed:
            configed_axes += 1
        if self._polarization_axis_configed:
            configed_axes += 1
        if self._obs_axis_configed:
            configed_axes += 1
        if self._custom_axis_configed:
            configed_axes += 1
        return configed_axes

    @property
    def update(self):
        return self._update

    @update.setter
    def update(self, value):
        self._update = value


class Hdf5ObsBlueprint(ObsBlueprint):
    """
    Class that specializes the CAOM2 Observation construction based on HDF5
    file content.

    The blueprint designates the source of each of these attributes as either
    HDF5 Dataset or Group values. Specific or default values may also be
    indicated in the same fashion os for an ObsBlueprint. The blueprint can
    be checked by simply displaying it.

    HDF5-specific example:
    # create a blueprint and customize it
    ob = Hdf5ObsBlueprint(position_axes=(1, 2)

    # lookup value starting with // means rooted at base of the hdf5 file
    ob.add_attribute('Observation.target.name', '//header/object/obj_id')

    # lookup value starting with / means rooted at the base of the
    # "find_roots_here" parameter for Hdf5Parser
    #
    # (integer) means return only the value with the index of "integer"
    #  from a list
    ob.add_attribute(
        'Chunk.position.axis.function.refCoord.coord1.pix',
        '/header/wcs/crpix(0)')

    # (integer:integer) means return only the value with the index of
    #   "integer" from a list, followed by "integer" from the list in the
    #   list
    ob.add_attribute(
        'Chunk.position.axis.function.cd11', '/header/wcs/cd(0:0)')
    print(ob)

    """
    def __init__(self, position_axes=None, energy_axis=None,
                 polarization_axis=None, time_axis=None,
                 obs_axis=None, custom_axis=None, module=None,
                 update=True, instantiated_class=None):
        """
        There are no sensible/known HDF5 defaults for WCS construction, so
        default to ensuring the blueprint executes with mostly values of None.

        Use the attribute _wcs_std, so that the list of WCS keywords used
        as input is known.
        """
        super().__init__(
            position_axes,
            energy_axis,
            polarization_axis,
            time_axis,
            obs_axis,
            custom_axis,
            module,
            update,
            instantiated_class,
        )
        tmp = {
               'Observation.algorithm.name': ([], 'exposure'),
               'Plane.calibrationLevel': ([], CalibrationLevel.RAW_STANDARD),
               'Plane.dataProductType': ([], DataProductType.IMAGE),
               'Artifact.releaseType': ([], ReleaseType.DATA),
               'Chunk': 'include'
        }
        # using the tmp to make sure that the keywords are valid
        for key in tmp:
            self.set(key, tmp[key])

    def configure_custom_axis(self, axis, override=True):
        """
        Set the expected custom keywords by index in the blueprint
        and the wcs_std lookup.

        :param axis: The index expected for the custom axis.
        :param override: Set to False when reading from a file.
        """
        if self._custom_axis_configed:
            self.logger.debug(
                'Attempt to configure already-configured custom axis.')
            return

        if override:
            self.set('Chunk.custom.axis.axis.ctype', ([], None))
            self.set('Chunk.custom.axis.axis.cunit', ([], None))
            self.set('Chunk.custom.axis.function.naxis', ([], 1))
            self.set('Chunk.custom.axis.function.delta', ([], None))
            self.set('Chunk.custom.axis.function.refCoord.pix', ([], None))
            self.set('Chunk.custom.axis.function.refCoord.val', ([], None))

        self._wcs_std['Chunk.custom.axis.axis.ctype'] = ''
        self._wcs_std['Chunk.custom.axis.axis.cunit'] = ''
        self._wcs_std['Chunk.custom.axis.function.naxis'] = ''
        self._wcs_std['Chunk.custom.axis.function.delta'] = ''
        self._wcs_std['Chunk.custom.axis.function.refCoord.pix'] = ''
        self._wcs_std['Chunk.custom.axis.function.refCoord.val'] = ''
        self._custom_axis_configed = True

    def configure_position_axes(self, axes, override=True):
        """
        Set the expected spatial keywords by indices in the blueprint and
        the wcs_std lookup.

        :param axes: The index expected for the position axes.
        :param override: Set to False when reading from a file.
        """
        if self._pos_axes_configed:
            self.logger.debug(
                'Attempt to configure already-configured position axes.')
            return

        if override:
            self.set('Chunk.position.coordsys', ([], None))
            self.set('Chunk.position.equinox', ([], None))
            self.set('Chunk.position.axis.axis1.ctype', ([], None))
            self.set('Chunk.position.axis.axis1.cunit', ([], None))
            self.set('Chunk.position.axis.axis2.ctype', ([], None))
            self.set('Chunk.position.axis.axis2.cunit', ([], None))
            self.set('Chunk.position.axis.error1.syser', ([], None))
            self.set('Chunk.position.axis.error1.rnder', ([], None))
            self.set('Chunk.position.axis.error2.syser', ([], None))
            self.set('Chunk.position.axis.error2.rnder', ([], None))
            self.set('Chunk.position.axis.function.cd11', ([], None))
            self.set('Chunk.position.axis.function.cd12', ([], None))
            self.set('Chunk.position.axis.function.cd21', ([], None))
            self.set('Chunk.position.axis.function.cd22', ([], None))
            self.set('Chunk.position.axis.function.dimension.naxis1',
                     ([], 1))
            self.set('Chunk.position.axis.function.dimension.naxis2',
                     ([], 1))
            self.set('Chunk.position.axis.function.refCoord.coord1.pix',
                     ([], None))
            self.set('Chunk.position.axis.function.refCoord.coord1.val',
                     ([], None))
            self.set('Chunk.position.axis.function.refCoord.coord2.pix',
                     ([], None))
            self.set('Chunk.position.axis.function.refCoord.coord2.val',
                     ([], None))

        self._wcs_std['Chunk.position.coordsys'] = ''
        self._wcs_std['Chunk.position.equinox'] = ''

        self._wcs_std['Chunk.position.axis.axis1.ctype'] = ''
        self._wcs_std['Chunk.position.axis.axis1.cunit'] = ''
        self._wcs_std['Chunk.position.axis.axis2.ctype'] = ''
        self._wcs_std['Chunk.position.axis.axis2.cunit'] = ''
        self._wcs_std['Chunk.position.axis.error1.syser'] = ''
        self._wcs_std['Chunk.position.axis.error1.rnder'] = ''
        self._wcs_std['Chunk.position.axis.error2.syser'] = ''
        self._wcs_std['Chunk.position.axis.error2.rnder'] = ''
        self._wcs_std['Chunk.position.axis.function.cd11'] = ''
        self._wcs_std['Chunk.position.axis.function.cd12'] = ''
        self._wcs_std['Chunk.position.axis.function.cd21'] = ''
        self._wcs_std['Chunk.position.axis.function.cd22'] = ''
        self._wcs_std['Chunk.position.axis.function.dimension.naxis1'] = ''
        self._wcs_std['Chunk.position.axis.function.dimension.naxis2'] = ''
        self._wcs_std['Chunk.position.axis.function.refCoord.coord1.pix'] = ''
        self._wcs_std['Chunk.position.axis.function.refCoord.coord1.val'] = ''
        self._wcs_std['Chunk.position.axis.function.refCoord.coord2.pix'] = ''
        self._wcs_std['Chunk.position.axis.function.refCoord.coord2.val'] = ''

        self._pos_axes_configed = True

    def configure_energy_axis(self, axis, override=True):
        """
        :param axis: The index expected for the energy axis.
        :param override: Set to False when reading from a file.
        :return:
        """
        if self._energy_axis_configed:
            self.logger.debug(
                'Attempt to configure already-configured energy axis.')
            return

        if override:
            self.set('Chunk.energy.specsys', ([], None))
            self.set('Chunk.energy.ssysobs', ([], None))
            self.set('Chunk.energy.restfrq', ([], None))
            self.set('Chunk.energy.restwav', ([], None))
            self.set('Chunk.energy.velosys', ([], None))
            self.set('Chunk.energy.zsource', ([], None))
            self.set('Chunk.energy.ssyssrc', ([], None))
            self.set('Chunk.energy.velang', ([], None))

            self.set('Chunk.energy.bandpassName', ([], None))
            self.set('Chunk.energy.resolvingPower', ([], None))

            self.set('Chunk.energy.axis.axis.ctype', ([], None))
            self.set('Chunk.energy.axis.axis.cunit', ([], None))
            self.set('Chunk.energy.axis.error.syser', ([], None))
            self.set('Chunk.energy.axis.error.rnder', ([], None))
            self.set('Chunk.energy.axis.function.naxis', ([], 1))
            self.set('Chunk.energy.axis.function.delta', ([], None))
            self.set('Chunk.energy.axis.function.refCoord.pix', ([], None))
            self.set('Chunk.energy.axis.function.refCoord.val', ([], None))

        self._wcs_std['Chunk.energy.specsys'] = ''
        self._wcs_std['Chunk.energy.ssysobs'] = ''
        self._wcs_std['Chunk.energy.restfrq'] = ''
        self._wcs_std['Chunk.energy.restwav'] = ''
        self._wcs_std['Chunk.energy.velosys'] = ''
        self._wcs_std['Chunk.energy.zsource'] = ''
        self._wcs_std['Chunk.energy.ssyssrc'] = ''
        self._wcs_std['Chunk.energy.velang'] = ''

        self._wcs_std['Chunk.energy.axis.axis.ctype'] = ''
        self._wcs_std['Chunk.energy.axis.axis.cunit'] = ''
        self._wcs_std['Chunk.energy.axis.error.syser'] = ''
        self._wcs_std['Chunk.energy.axis.error.rnder'] = ''
        self._wcs_std['Chunk.energy.axis.function.naxis'] = ''
        self._wcs_std['Chunk.energy.axis.function.delta'] = ''
        self._wcs_std['Chunk.energy.axis.function.refCoord.pix'] = ''
        self._wcs_std['Chunk.energy.axis.function.refCoord.val'] = ''
        self._energy_axis_configed = True

    def configure_polarization_axis(self, axis, override=True):
        """
        Set the expected polarization keywords by index in the blueprint
        and the wcs_std lookup.

        :param axis: The index expected for the polarization axis.
        :param override: Set to False when reading from a file.
        :return:
        """
        if self._polarization_axis_configed:
            self.logger.debug(
                'Attempt to configure already-configured polarization axis.')
            return

        if override:
            # STOKES is the only value allowed for PolarizationWCS ctype.
            self.set('Chunk.polarization.axis.axis.ctype', ([], 'STOKES'))
            self.set('Chunk.polarization.axis.axis.cunit', ([], None))
            self.set('Chunk.polarization.axis.function.naxis', ([], 1))
            self.set('Chunk.polarization.axis.function.delta', ([], None))
            self.set('Chunk.polarization.axis.function.refCoord.pix',
                     ([], None))
            self.set('Chunk.polarization.axis.function.refCoord.val',
                     ([], None))

        self._wcs_std['Chunk.polarization.axis.axis.ctype'] = ''
        self._wcs_std['Chunk.polarization.axis.axis.cunit'] = ''
        self._wcs_std['Chunk.polarization.axis.function.naxis'] = ''
        self._wcs_std['Chunk.polarization.axis.function.delta'] = ''
        self._wcs_std['Chunk.polarization.axis.function.refCoord.pix'] = ''
        self._wcs_std['Chunk.polarization.axis.function.refCoord.val'] = ''

        self._polarization_axis_configed = True

    def configure_observable_axis(self, axis, override=True):
        """
        Set the expected observable keywords by index in the blueprint
        and the wcs_std lookup.
        Note: observable axis is not a standard WCS and it's not used by
        astropy.wcs so, arguably, it can be removed. It is here for now for
        consistency purposes.
        :param axis: The index expected for the observable axis.
        :param override: Set to False when reading from a file.
        :return:
        """
        if self._obs_axis_configed:
            self.logger.debug(
                'Attempt to configure already-configured observable axis.')
            return

        if override:
            self.set('Chunk.observable.axis.axis.ctype', ([], None))
            self.set('Chunk.observable.axis.axis.cunit', ([], None))
            self.set('Chunk.observable.axis.function.refCoord.pix', ([], None))

        self._wcs_std['Chunk.observable.axis.axis.ctype'] = ''
        self._wcs_std['Chunk.observable.axis.axis.cunit'] = ''
        self._wcs_std['Chunk.observable.axis.function.refCoord.pix'] = ''

        self._obs_axis_configed = True

    def configure_time_axis(self, axis, override=True):
        """
        Set the expected time keywords by index in the blueprint and
        the wcs_std lookup.

        :param axis: The index expected for the time axis.
        :param override: Set to False when reading from a file.
        :return:
        """
        if self._time_axis_configed:
            self.logger.debug(
                'Attempt to configure already-configured time axis.')
            return

        if override:
            self.set('Chunk.time.exposure', ([], None))
            self.set('Chunk.time.timesys', ([], None))
            self.set('Chunk.time.trefpos', ([], None))
            self.set('Chunk.time.mjdref', ([], None))
            self.set('Chunk.time.resolution', ([], None))
            self.set('Chunk.time.axis.axis.ctype', ([], None))
            self.set('Chunk.time.axis.axis.cunit', ([], None))
            self.set('Chunk.time.axis.error.syser', ([], None))
            self.set('Chunk.time.axis.error.rnder', ([], None))
            self.set('Chunk.time.axis.function.naxis', ([], 1))
            self.set('Chunk.time.axis.function.delta', ([], None))
            self.set('Chunk.time.axis.function.refCoord.pix', ([], None))
            self.set('Chunk.time.axis.function.refCoord.val', ([], None))

        self._wcs_std['Chunk.time.exposure'] = ''
        self._wcs_std['Chunk.time.resolution'] = ''
        self._wcs_std['Chunk.time.timesys'] = ''
        self._wcs_std['Chunk.time.trefpos'] = ''
        self._wcs_std['Chunk.time.mjdref'] = ''

        self._wcs_std['Chunk.time.axis.axis.ctype'] = ''
        self._wcs_std['Chunk.time.axis.axis.cunit'] = ''
        self._wcs_std['Chunk.time.axis.error.syser'] = ''
        self._wcs_std['Chunk.time.axis.error.rnder'] = ''
        self._wcs_std['Chunk.time.axis.function.naxis'] = ''
        self._wcs_std['Chunk.time.axis.function.delta'] = ''
        self._wcs_std['Chunk.time.axis.function.refCoord.pix'] = ''
        self._wcs_std['Chunk.time.axis.function.refCoord.val'] = ''

        self._time_axis_configed = True

    def set(self, caom2_element, value, extension=0):
        """
        Sets the value associated with an element in the CAOM2 model. Value
        cannot be a tuple.
        :param caom2_element: name CAOM2 element (as in
        ObsBlueprint.CAOM2_ELEMEMTS)
        :param value: new value of the CAOM2 element
        :param extension: extension number (used only for Chunk elements)
        """
        if hasattr(value, 'decode'):
            value = value.decode('utf-8')
        super().set(caom2_element, value, extension)

    def _guess_axis_info(self):
        self._guess_axis_info_from_plan()


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
        if (self.blueprint._module is not None or
                self.blueprint._module_instance is not None):
            for key, value in plan.items():
                if ObsBlueprint.is_function(value):
                    if self._blueprint._module_instance is None:
                        plan[key] = self._execute_external(value, key, 0)
                    else:
                        plan[key] = self._execute_external_instance(
                            value, key, 0)

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
        self.logger.debug(
            f'Begin CAOM2 observation augmentation for URI {artifact_uri}.')
        if observation is None or not isinstance(observation, Observation):
            raise ValueError(
                f'Observation type mis-match for {observation}.')

        observation.meta_release = self._get_datetime(self._get_from_list(
            'Observation.metaRelease', index=0,
            current=observation.meta_release))
        observation.meta_read_groups = self._get_from_list(
            'Observation.metaReadGroups', index=0,
            current=observation.meta_read_groups)
        observation.meta_producer = self._get_from_list(
            'Observation.metaProducer', index=0,
            current=observation.meta_producer)

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
            plane = Plane(product_id=product_id)
            observation.planes[product_id] = plane
        self.augment_plane(plane, artifact_uri)
        self.logger.debug(
            f'End CAOM2 observation augmentation for {artifact_uri}.')

    def augment_plane(self, plane, artifact_uri):
        """
        Augments a given plane with artifact structure only.
        :param plane: existing CAOM2 plane to be augmented.
        :param artifact_uri:
        """
        self.logger.debug(
            f'Begin CAOM2 plane augmentation for {artifact_uri}.')
        if plane is None or not isinstance(plane, Plane):
            raise ValueError(f'Plane type mis-match for {plane}')

        plane.meta_release = self._get_datetime(self._get_from_list(
            'Plane.metaRelease', index=0, current=plane.meta_release))
        plane.data_release = self._get_datetime(self._get_from_list(
            'Plane.dataRelease', index=0, current=plane.data_release))
        plane.data_product_type = self._to_data_product_type(
            self._get_from_list('Plane.dataProductType', index=0,
                                current=plane.data_product_type))
        plane.calibration_level = self._to_calibration_level(_to_int_32(
            self._get_from_list('Plane.calibrationLevel', index=0,
                                current=plane.calibration_level)))
        plane.meta_producer = self._get_from_list(
            'Plane.metaProducer', index=0, current=plane.meta_producer)

        artifact = None
        for ii in plane.artifacts:
            artifact = plane.artifacts[ii]
            if artifact.uri == artifact_uri:
                break
        if artifact is None or artifact.uri != artifact_uri:
            artifact = Artifact(artifact_uri, self._to_product_type(
                self._get_from_list('Artifact.productType', index=0)),
                                self._to_release_type(self._get_from_list(
                                    'Artifact.releaseType', index=0)))
            plane.artifacts[artifact_uri] = artifact
        self.augment_artifact(artifact, 0)
        self.logger.debug(
            f'End CAOM2 plane augmentation for {artifact_uri}.')

    def augment_artifact(self, artifact, index):
        """
        Augments a given CAOM2 artifact with available information
        :param artifact: existing CAOM2 artifact to be augmented
        :param index: int Part name, used in specializing classes
        """
        self.logger.debug(f'Begin CAOM2 artifact augmentation for {self.uri}.')
        if artifact is None or not isinstance(artifact, Artifact):
            raise ValueError(
                f'Artifact type mis-match for {artifact}')

        artifact.product_type = self._to_product_type(self._get_from_list(
            'Artifact.productType', index=0, current=artifact.product_type))
        artifact.release_type = self._to_release_type(self._get_from_list(
            'Artifact.releaseType', index=0, current=artifact.release_type))
        artifact.content_type = self._get_from_list(
            'Artifact.contentType', index=0, current=artifact.content_type)
        artifact.content_length = self._get_from_list(
            'Artifact.contentLength', index=0, current=artifact.content_length)
        artifact.content_checksum = _to_checksum_uri(self._get_from_list(
            'Artifact.contentChecksum', index=0,
            current=artifact.content_checksum))
        artifact.content_release = self._get_from_list(
            'Artifact.contentRelease', index=0,
            current=artifact.content_release)
        artifact.content_read_groups = self._get_from_list(
            'Artifact.contentReadGroups', index=0,
            current=artifact.content_read_groups)
        artifact.meta_producer = self._get_from_list(
            'Artifact.metaProducer', index=0, current=artifact.meta_producer)
        self.logger.debug(f'End CAOM2 artifact augmentation for {self.uri}.')

    def _get_from_list(self, lookup, index, current=None):
        value = None
        try:
            keywords = self.blueprint._get(lookup)
        except KeyError:
            self.add_error(lookup, sys.exc_info()[1])
            self.logger.debug(
                f'Could not find {lookup} in configuration.')
            if current:
                self.logger.debug(
                    f'{lookup}: using current value of {current!r}.')
                value = current
            return value
        if (keywords and not ObsBlueprint.needs_lookup(keywords)
                and not ObsBlueprint.is_function(keywords)):
            value = keywords
        elif self._blueprint.update:
            # The first clause: boolean attributes are used to represent
            # three different values: True, False, and unknown. For boolean
            # attributes _only_ assessed that the risk of setting to None
            # accidentally was better than being unable to set a value of
            # 'unknown'.
            #
            # The second clause: the default value for the current parameter
            # in the method signature is 'None', so do not want to
            # inadvertently assign the default value.
            #
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
            self.logger.debug(f'Could not find \'{lookup}\' in caom2blueprint '
                              f'configuration.')

        # if there's something useful as a value in the keywords,
        # extract it
        if keywords:
            if ObsBlueprint.needs_lookup(keywords):
                # if there's a default value use it
                if keywords[1]:
                    value = keywords[1]
                    self.logger.debug(
                        f'{lookup}: assigned default value {value}.')
            elif not ObsBlueprint.is_function(keywords):
                value = keywords
                self.logger.debug(f'{lookup}: assigned value {value}.')
        return value

    def add_error(self, key, message):
        self._errors.append('{} {} {}'.format(
            datetime.now().strftime('%Y-%m-%dT%H:%M:%S'), key, message))

    def _to_data_product_type(self, value):
        return self._to_enum_type(value, DataProductType)

    def _to_calibration_level(self, value):
        return self._to_enum_type(value, CalibrationLevel)

    def _to_product_type(self, value):
        return self._to_enum_type(value, ProductType)

    def _to_release_type(self, value):
        return self._to_enum_type(value, ReleaseType)

    def _to_enum_type(self, value, to_enum_type):
        if value is None:
            raise ValueError(
                f'Must set a value of {to_enum_type.__name__} for '
                f'{self.uri}.')
        elif isinstance(value, to_enum_type):
            return value
        else:
            return to_enum_type(value)

    def _execute_external(self, value, key, extension):
        """Execute a function supplied by a user, assign a value to a
        blueprint entry. The input parameters passed to the function are the
        headers as read in by astropy, or the artifact uri.

        :param value the name of the function to apply.
        :param key:
        :param extension: the current extension name or number.
        """
        # determine which of the possible values for parameter the user
        # is hoping for
        if 'uri' in value:
            parameter = self.uri
        elif 'header' in value and isinstance(self, FitsParser):
            parameter = self._headers[extension]
        elif isinstance(self, FitsParser):
            parameter = {'uri': self.uri,
                         'header': self._headers[extension]}
        else:
            if hasattr(self, '_file'):
                parameter = {'base': self._file}
            else:
                parameter = {'uri': self.uri,
                             'header': None}

        result = ''
        execute = None
        try:
            execute = getattr(self.blueprint._module, value.split('(')[0])
        except Exception as e:
            msg = 'Failed to find {}.{} for {}'.format(
                self.blueprint._module.__name__, value.split('(')[0], key)
            logging.error(msg)
            self._errors.append(msg)
            tb = traceback.format_exc()
            logging.debug(tb)
            logging.error(e)
        try:
            result = execute(parameter)
            logging.debug(
                f'Key {key} calculated value of {result} using {value} type '
                f'{type(result)}')
        except Exception as e:
            msg = 'Failed to execute {} for {} in {}'.format(
                execute.__name__, key, self.uri)
            logging.error(msg)
            logging.debug('Input parameter was {}, value was {}'.format(
                parameter, value))
            self._errors.append(msg)
            tb = traceback.format_exc()
            logging.debug(tb)
            logging.error(e)
        return result

    def _execute_external_instance(self, value, key, extension):
        """Execute a function supplied by a user, assign a value to a
        blueprint entry. The input parameters passed to the function are the
        headers as read in by astropy, or the artifact uri.

        :param value the name of the function to apply.
        :param key:
        :param extension: the current extension name or number.
        :raise Caom2Exception exception raised when there is a recognizable
            error in the information being used to create a CAOM2 record. A
            correct and consistent CAOM2 record cannot be created from the
            input metadata. The client should treat the Observation instance
            under construction as invalid.
        """
        result = ''
        try:
            execute = getattr(
                self.blueprint._module_instance, value.split('(')[0])
        except Exception as e:
            msg = 'Failed to find {}.{} for {}'.format(
                self.blueprint._module_instance.__class__.__name__,
                value.split('(')[0], key)
            logging.error(msg)
            self._errors.append(msg)
            tb = traceback.format_exc()
            logging.debug(tb)
            logging.error(e)
            return result
        try:
            result = execute(extension)
            logging.debug(
                'Key {} calculated value of {} using {}'.format(
                    key, result, value))
        except ValueError as e2:
            # DB 23-03-22
            # Anything that you can do to make the CAOM2 record creation fail
            # in this case of bad WCS metadata would be useful. Use
            # ValueError because that happens to be what astropy is throwing
            # for a SkyCoord construction failure.
            raise Caom2Exception(e2)
        except Exception as e:
            msg = 'Failed to execute {} for {} in {}'.format(
                execute, key, self.uri)
            logging.error(msg)
            logging.debug('Input value was {}'.format(value))
            self._errors.append(msg)
            tb = traceback.format_exc()
            logging.debug(tb)
            logging.error(e)
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
                for dt_format in ['%Y-%m-%dT%H:%M:%S', '%Y-%m-%dT%H:%M:%S.%f',
                                  '%Y-%m-%d %H:%M:%S.%f', '%Y-%m-%d',
                                  '%Y/%m/%d %H:%M:%S', '%Y-%m-%d %H:%M:%S',
                                  '%Y/%m/%d,%H:%M:%S', '%Y/%m/%d',
                                  '%d/%m/%y', '%d/%m/%y %H:%M:%S', '%d-%m-%Y']:
                    try:
                        result = datetime.strptime(from_value, dt_format)
                    except ValueError:
                        pass

                if result is None:
                    self.logger.error('Cannot parse datetime {}'.format(
                            from_value))
                    self.add_error('get_datetime', sys.exc_info()[1])
                return result
        else:
            return None


class ContentParser(BlueprintParser):

    def __init__(self, obs_blueprint=None, uri=None):
        super().__init__(obs_blueprint, uri)
        self._wcs_parser = WcsParser()

    def _get_chunk_naxis(self, chunk, index):
        chunk.naxis = self._get_from_list(
            'Chunk.naxis', index, self._wcs_parser.wcs.wcs.naxis)

    def augment_artifact(self, artifact, index):
        """
        Augments a given CAOM2 artifact with available content information
        :param artifact: existing CAOM2 artifact to be augmented
        :param index: int Part name
        """
        super().augment_artifact(artifact, index)

        self.logger.debug(
            f'Begin content artifact augmentation for {artifact.uri}')

        if self.blueprint.get_configed_axes_count() == 0:
            raise TypeError(
                f'No WCS Data. End content artifact augmentation for '
                f'{artifact.uri}.')

        if self.ignore_chunks(artifact, index):
            return

        part = artifact.parts[str(index)]
        part.product_type = self._get_from_list('Part.productType', index)
        part.meta_producer = self._get_from_list(
            'Part.metaProducer', index=0, current=part.meta_producer)

        # each Part has one Chunk, if it's not an empty part as determined
        # just previously
        if not part.chunks:
            part.chunks.append(Chunk())
        chunk = part.chunks[0]
        chunk.meta_producer = self._get_from_list(
            'Chunk.metaProducer', index=0, current=chunk.meta_producer)

        self._get_chunk_naxis(chunk, index)
        if self.blueprint._pos_axes_configed:
            self._wcs_parser.augment_position(chunk)
            if chunk.position is None:
                self._try_position_with_blueprint(chunk, index)
        if chunk.position:
            chunk.position.resolution = self._get_from_list(
                'Chunk.position.resolution', index=index)
        if self.blueprint._energy_axis_configed:
            self._wcs_parser.augment_energy(chunk)
        if chunk.energy:
            chunk.energy.bandpass_name = self._get_from_list(
                'Chunk.energy.bandpassName', index=index)
            chunk.energy.transition = self._get_energy_transition(
                chunk.energy.transition)
            chunk.energy.resolving_power = _to_float(self._get_from_list(
                'Chunk.energy.resolvingPower', index=index))
        else:
            if self.blueprint._energy_axis_configed:
                self._try_energy_with_blueprint(chunk, index)
        if self.blueprint._time_axis_configed:
            self._wcs_parser.augment_temporal(chunk)
            if chunk.time is None:
                self._try_time_with_blueprint(chunk, index)
        if self.blueprint._polarization_axis_configed:
            self._wcs_parser.augment_polarization(chunk)
            if chunk.polarization is None:
                self._try_polarization_with_blueprint(chunk, index)
        if self.blueprint._obs_axis_configed:
            self._wcs_parser.augment_observable(chunk)
            if chunk.observable is None and chunk.observable_axis is None:
                self._try_observable_with_blueprint(chunk, index)
        if self.blueprint._custom_axis_configed:
            self._wcs_parser.augment_custom(chunk)

        # try to set smaller bits of the chunk WCS elements from the
        # blueprint
        self._try_range_with_blueprint(chunk, index)

        self.logger.debug(
            f'End content artifact augmentation for {artifact.uri}.')

    def augment_observation(self, observation, artifact_uri, product_id=None):
        """
        Augments a given observation with available content information.
        :param observation: existing CAOM2 observation to be augmented.
        :param artifact_uri: the key for finding the artifact to augment
        :param product_id: the key for finding for the plane to augment
        """
        super().augment_observation(observation, artifact_uri, product_id)
        self.logger.debug(
            f'Begin content observation augmentation for URI {artifact_uri}.')
        members = self._get_members(observation)
        if members:
            if isinstance(members, TypedSet):
                for m in members:
                    observation.members.add(m)
            else:
                for m in members.split():
                    observation.members.add(ObservationURI(m))
        observation.algorithm = self._get_algorithm(observation)

        observation.sequence_number = _to_int(self._get_from_list(
            'Observation.sequenceNumber', index=0))
        observation.intent = self._get_from_list(
            'Observation.intent', 0, (ObservationIntentType.SCIENCE if
                                      observation.intent is None else
                                      observation.intent))
        observation.type = self._get_from_list('Observation.type', 0,
                                               current=observation.type)
        observation.meta_release = self._get_datetime(
            self._get_from_list('Observation.metaRelease', 0,
                                current=observation.meta_release))
        observation.meta_read_groups = self._get_from_list(
            'Observation.metaReadGroups', 0)
        observation.meta_producer = self._get_from_list(
            'Observation.metaProducer', 0, current=observation.meta_producer)
        observation.requirements = self._get_requirements(
            observation.requirements)
        observation.instrument = self._get_instrument(observation.instrument)
        observation.proposal = self._get_proposal(observation.proposal)
        observation.target = self._get_target(observation.target)
        observation.target_position = self._get_target_position(
            observation.target_position)
        observation.telescope = self._get_telescope(observation.telescope)
        observation.environment = self._get_environment(
            observation.environment)
        self.logger.debug(
            f'End content observation augmentation for {artifact_uri}.')

    def augment_plane(self, plane, artifact_uri):
        """
        Augments a given plane with available content information.
        :param plane: existing CAOM2 plane to be augmented.
        :param artifact_uri:
        """
        super().augment_plane(plane, artifact_uri)
        self.logger.debug(
            f'Begin content plane augmentation for {artifact_uri}.')

        plane.meta_release = self._get_datetime(self._get_from_list(
            'Plane.metaRelease', index=0, current=plane.meta_release))
        plane.data_release = self._get_datetime(self._get_from_list(
            'Plane.dataRelease', index=0))
        plane.data_product_type = self._to_data_product_type(
            self._get_from_list('Plane.dataProductType', index=0,
                                current=plane.data_product_type))
        plane.calibration_level = self._to_calibration_level(_to_int_32(
            self._get_from_list('Plane.calibrationLevel', index=0,
                                current=plane.calibration_level)))
        plane.meta_producer = self._get_from_list(
            'Plane.metaProducer', index=0, current=plane.meta_producer)
        plane.observable = self._get_observable(current=plane.observable)
        plane.provenance = self._get_provenance(plane.provenance)
        plane.metrics = self._get_metrics(current=plane.metrics)
        plane.quality = self._get_quality(current=plane.quality)

        self.logger.debug(
            f'End content plane augmentation for {artifact_uri}.')

    def _get_algorithm(self, obs):
        """
        Create an Algorithm instance populated with available content
        information.
        :return: Algorithm
        """
        self.logger.debug('Begin Algorithm augmentation.')
        # TODO DEFAULT VALUE
        name = self._get_from_list('Observation.algorithm.name', index=0,
                                   current=obs.algorithm.name)
        result = Algorithm(str(name)) if name else None
        self.logger.debug('End Algorithm augmentation.')
        return result

    def _get_energy_transition(self, current):
        """
        Create an EnergyTransition instance populated with available content
        information.
        :return: EnergyTransition
        """
        self.logger.debug('Begin EnergyTransition augmentation.')
        species = self._get_from_list(
            'Chunk.energy.transition.species', index=0,
            current=None if current is None else current.species)
        transition = self._get_from_list(
            'Chunk.energy.transition.transition', index=0,
            current=None if current is None else current.transition)
        result = None
        if species is not None and transition is not None:
            result = EnergyTransition(species, transition)
        self.logger.debug('End EnergyTransition augmentation.')
        return result

    def _get_environment(self, current):
        """
        Create an Environment instance populated with available content
        information.
        :current Environment instance, if one already exists in the
            Observation
        :return: Environment
        """
        self.logger.debug('Begin Environment augmentation.')
        seeing = self._get_from_list(
            'Observation.environment.seeing', index=0,
            current=None if current is None else current.seeing)
        humidity = _to_float(
            self._get_from_list(
                'Observation.environment.humidity', index=0,
                current=None if current is None else current.humidity))
        elevation = self._get_from_list(
            'Observation.environment.elevation', index=0,
            current=None if current is None else current.elevation)
        tau = self._get_from_list(
            'Observation.environment.tau', index=0,
            current=None if current is None else current.tau)
        wavelength_tau = self._get_from_list(
            'Observation.environment.wavelengthTau', index=0,
            current=None if current is None else current.wavelength_tau)
        ambient = _to_float(
            self._get_from_list(
                'Observation.environment.ambientTemp', index=0,
                current=None if current is None else current.ambient_temp))
        photometric = self._cast_as_bool(self._get_from_list(
            'Observation.environment.photometric', index=0,
            current=None if current is None else current.photometric))
        enviro = None
        if seeing or humidity or elevation or tau or wavelength_tau or ambient:
            enviro = Environment()
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
        Create an Instrument instance populated with available content
        information.
        :return: Instrument
        """
        self.logger.debug('Begin Instrument augmentation.')
        name = self._get_from_list(
            'Observation.instrument.name', index=0,
            current=None if current is None else current.name)
        keywords = self._get_set_from_list(
            'Observation.instrument.keywords', index=0)
        instr = None
        if name:
            instr = Instrument(str(name))
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
        if (isinstance(obs, SimpleObservation) and
            (self.blueprint._get('DerivedObservation.members') or
             self.blueprint._get('CompositeObservation.members'))):
            raise TypeError(
                'Cannot apply blueprint for DerivedObservation to a '
                'simple observation')
        elif isinstance(obs, DerivedObservation):
            lookup = self.blueprint._get('DerivedObservation.members',
                                         extension=1)
            if ObsBlueprint.is_table(lookup) and len(self.headers) > 1:
                member_list = self._get_from_table(
                    'DerivedObservation.members', 1)
                # ensure the members are good little ObservationURIs
                if member_list.startswith('caom:'):
                    members = member_list
                else:
                    members = ' '.join(['caom:{}/{}'.format(
                        obs.collection, i) if not i.startswith('caom') else i
                                        for i in member_list.split()])
            else:
                if obs.members is None:
                    members = self._get_from_list(
                        'DerivedObservation.members', index=0)
                else:
                    members = self._get_from_list(
                        'DerivedObservation.members', index=0,
                        current=obs.members)
        elif isinstance(obs, CompositeObservation):
            lookup = self.blueprint._get('CompositeObservation.members',
                                         extension=1)
            if ObsBlueprint.is_table(lookup) and len(self.headers) > 1:
                member_list = self._get_from_table(
                    'CompositeObservation.members', 1)
                # ensure the members are good little ObservationURIs
                if member_list.startswith('caom:'):
                    members = member_list
                else:
                    members = ' '.join(['caom:{}/{}'.format(
                        obs.collection, i) if not i.startswith('caom') else i
                                        for i in member_list.split()])
            else:
                if obs.members is None:
                    members = self._get_from_list(
                        'CompositeObservation.members', index=0)
                else:
                    members = self._get_from_list(
                        'CompositeObservation.members', index=0,
                        current=obs.members)
        self.logger.debug('End Members augmentation.')
        return members

    def _get_metrics(self, current):
        """
        Create a Metrics instance populated with available content information.
        :return: Metrics
        """
        self.logger.debug('Begin Metrics augmentation.')
        source_number_density = self._get_from_list(
            'Plane.metrics.sourceNumberDensity', index=0,
            current=None if current is None else current.source_number_density)
        background = self._get_from_list(
            'Plane.metrics.background', index=0,
            current=None if current is None else current.background)
        background_stddev = self._get_from_list(
            'Plane.metrics.backgroundStddev', index=0,
            current=None if current is None else current.background_std_dev)
        flux_density_limit = self._get_from_list(
            'Plane.metrics.fluxDensityLimit', index=0,
            current=None if current is None else current.flux_density_limit)
        mag_limit = self._get_from_list(
            'Plane.metrics.magLimit', index=0,
            current=None if current is None else current.mag_limit)
        sample_snr = self._get_from_list(
            'Plane.metrics.sampleSNR', index=0,
            current=None if current is None else current.sample_snr)

        metrics = None
        if (source_number_density or background or background_stddev or
                flux_density_limit or mag_limit or sample_snr):
            metrics = Metrics()
            metrics.source_number_density = source_number_density
            metrics.background = background
            metrics.background_std_dev = background_stddev
            metrics.flux_density_limit = flux_density_limit
            metrics.mag_limit = mag_limit
            metrics.sample_snr = sample_snr
        self.logger.debug('End Metrics augmentation.')
        return metrics

    def _get_naxis(self, label, index):
        """Helper function to construct a CoordAxis1D instance, with all
        it's members, from the blueprint.

        :param label: axis name - must be one of 'energy', 'time', or
        'polarization', as it's used for the blueprint lookup.
        :param index: which blueprint index to find a value in
        :return an instance of CoordAxis1D
        """
        self.logger.debug(
            f'Begin {label} naxis construction from blueprint.')

        aug_axis_ctype = self._get_from_list(
            f'Chunk.{label}.axis.axis.ctype', index)
        aug_axis_cunit = self._get_from_list(
            f'Chunk.{label}.axis.axis.cunit', index)
        aug_axis = None
        if aug_axis_ctype is not None:
            aug_axis = Axis(aug_axis_ctype, aug_axis_cunit)
            self.logger.debug(
                f'Creating {label} Axis for {self.uri} from blueprint')

        aug_error = self._two_param_constructor(
            f'Chunk.{label}.axis.error.syser',
            f'Chunk.{label}.axis.error.rnder',
            index, _to_float, CoordError)
        aug_ref_coord = self._two_param_constructor(
            f'Chunk.{label}.axis.function.refCoord.pix',
            f'Chunk.{label}.axis.function.refCoord.val',
            index, _to_float, RefCoord)
        aug_delta = _to_float(
            self._get_from_list(f'Chunk.{label}.axis.function.delta',
                                index))
        aug_length = _to_int(
            self._get_from_list(f'Chunk.{label}.axis.function.naxis',
                                index))

        aug_function = None
        if (aug_length is not None and aug_delta is not None and
                aug_ref_coord is not None):
            aug_function = \
                CoordFunction1D(aug_length, aug_delta, aug_ref_coord)
            self.logger.debug(
                f'Creating {label} function for {self.uri} from blueprint')

        aug_naxis = None
        if aug_function is None:
            aug_range = self._try_range_return(index, label)
            if aug_axis is not None and aug_range is not None:
                aug_naxis = CoordAxis1D(
                    axis=aug_axis, error=aug_error, range=aug_range)
                self.logger.debug(
                    f'Creating range {label} CoordAxis1D for {self.uri} from '
                    f'blueprint')
        else:
            if aug_axis is not None and aug_function is not None:
                aug_naxis = CoordAxis1D(aug_axis, aug_error, None, None,
                                        aug_function)
                self.logger.debug(
                    f'Creating function {label} CoordAxis1D for {self.uri} '
                    f'from blueprint')
        self.logger.debug(
            f'End {label} naxis construction from blueprint.')
        return aug_naxis

    def _get_observable(self, current):
        """
        Create a Observable instance populated with available content
        information.
        :return: Observable
        """
        self.logger.debug('Begin Observable augmentation.')
        ucd = self._get_from_list(
            'Plane.observable.ucd', index=0,
            current=None if current is None else current.ucd)
        observable = Observable(ucd) if ucd else None
        self.logger.debug('End Observable augmentation.')
        return observable

    def _get_proposal(self, current):
        """
        Create a Proposal instance populated with available content
        information.
        :return: Proposal
        """
        self.logger.debug('Begin Proposal augmentation.')
        prop_id = self._get_from_list(
            'Observation.proposal.id', index=0,
            current=None if current is None else current.id)
        pi = self._get_from_list(
            'Observation.proposal.pi', index=0,
            current=None if current is None else current.pi_name)
        project = self._get_from_list(
            'Observation.proposal.project', index=0,
            current=None if current is None else current.project)
        title = self._get_from_list(
            'Observation.proposal.title', index=0,
            current=None if current is None else current.title)
        keywords = self._get_set_from_list(
            'Observation.proposal.keywords', index=0)
        proposal = current
        if prop_id:
            proposal = Proposal(str(prop_id), pi, project, title)
            ContentParser._add_keywords(keywords, current, proposal)
        self.logger.debug(f'End Proposal augmentation {prop_id}.')
        return proposal

    def _get_provenance(self, current):
        """
        Create a Provenance instance populated with available Content
        information.
        :return: Provenance
        """
        self.logger.debug('Begin Provenance augmentation.')
        name = _to_str(
            self._get_from_list(
                'Plane.provenance.name', index=0,
                current=None if current is None else current.name))
        p_version = _to_str(self._get_from_list(
            'Plane.provenance.version', index=0,
            current=None if current is None else current.version))
        project = _to_str(
            self._get_from_list(
                'Plane.provenance.project', index=0,
                current=None if current is None else current.project))
        producer = _to_str(
            self._get_from_list(
                'Plane.provenance.producer', index=0,
                current=None if current is None else current.producer))
        run_id = _to_str(
            self._get_from_list(
                'Plane.provenance.runID', index=0,
                current=None if current is None else current.run_id))
        reference = _to_str(
            self._get_from_list(
                'Plane.provenance.reference', index=0,
                current=None if current is None else current.reference))
        last_executed = self._get_datetime(
            self._get_from_list(
                'Plane.provenance.lastExecuted', index=0,
                current=None if current is None else current.last_executed))
        keywords = self._get_set_from_list(
            'Plane.provenance.keywords', index=0)
        inputs = self._get_set_from_list('Plane.provenance.inputs', index=0)
        prov = None
        if name:
            prov = Provenance(name, p_version, project, producer, run_id,
                              reference, last_executed)
            ContentParser._add_keywords(keywords, current, prov)
            if inputs:
                if isinstance(inputs, TypedSet):
                    for i in inputs:
                        prov.inputs.add(i)
                else:
                    for i in inputs.split():
                        prov.inputs.add(PlaneURI(str(i)))
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
        flag = self._get_from_list(
            'Plane.dataQuality', index=0,
            current=None if current is None else current.flag)
        quality = DataQuality(flag) if flag else None
        self.logger.debug('End Quality augmentation.')
        return quality

    def _get_requirements(self, current):
        """
        Create a Requirements instance populated with available content
        information.
        :return: Requirements
        """
        self.logger.debug('Begin Requirement augmentation.')
        flag = self._get_from_list(
            'Observation.requirements.flag', index=0,
            current=None if current is None else current.flag)
        reqts = Requirements(flag) if flag else None
        self.logger.debug('End Requirement augmentation.')
        return reqts

    def _get_target(self, current):
        """
        Create a Target instance populated with available content information.
        :return: Target
        """
        self.logger.debug('Begin Target augmentation.')
        name = self._get_from_list(
            'Observation.target.name', index=0,
            current=None if current is None else current.name)
        target_type = self._get_from_list(
            'Observation.target.type', index=0,
            current=None if current is None else current.target_type)
        standard = self._cast_as_bool(self._get_from_list(
            'Observation.target.standard', index=0,
            current=None if current is None else current.standard))
        redshift = self._get_from_list(
            'Observation.target.redshift', index=0,
            current=None if current is None else current.redshift)
        keywords = self._get_set_from_list(
            'Observation.target.keywords', index=0)
        moving = self._cast_as_bool(
            self._get_from_list(
                'Observation.target.moving', index=0,
                current=None if current is None else current.moving))
        target_id = self._get_from_list(
            'Observation.target.targetID', index=0,
            current=None if current is None else current.target_id)
        target = None
        if name:
            target = Target(str(name), target_type, standard, redshift,
                            moving=moving, target_id=target_id)
            ContentParser._add_keywords(keywords, current, target)
        self.logger.debug('End Target augmentation.')
        return target

    def _get_target_position(self, current):
        """
        Create a Target Position instance populated with available content
        information.
        :return: Target Position
        """
        self.logger.debug('Begin CAOM2 TargetPosition augmentation.')
        x = self._get_from_list(
            'Observation.target_position.point.cval1', index=0,
            current=None if current is None else current.coordinates.cval1)
        y = self._get_from_list(
            'Observation.target_position.point.cval2', index=0,
            current=None if current is None else current.coordinates.cval2)
        coordsys = self._get_from_list(
            'Observation.target_position.coordsys', index=0,
            current=None if current is None else current.coordsys)
        equinox = self._get_from_list(
            'Observation.target_position.equinox', index=0,
            current=None if current is None else current.equinox)
        aug_target_position = None
        if x and y:
            aug_point = Point(x, y)
            aug_target_position = TargetPosition(aug_point, coordsys)
            aug_target_position.equinox = _to_float(equinox)
        self.logger.debug('End CAOM2 TargetPosition augmentation.')
        return aug_target_position

    def _get_telescope(self, current):
        """
        Create a Telescope instance populated with available content
        information.
        :return: Telescope
        """
        self.logger.debug('Begin Telescope augmentation.')
        name = self._get_from_list(
            'Observation.telescope.name', index=0,
            current=None if current is None else current.name)
        geo_x = _to_float(
            self._get_from_list(
                'Observation.telescope.geoLocationX', index=0,
                current=None if current is None else current.geo_location_x))
        geo_y = _to_float(
            self._get_from_list(
                'Observation.telescope.geoLocationY', index=0,
                current=None if current is None else current.geo_location_y))
        geo_z = _to_float(
            self._get_from_list(
                'Observation.telescope.geoLocationZ', index=0,
                current=None if current is None else current.geo_location_z))
        keywords = self._get_set_from_list(
            'Observation.telescope.keywords', index=0)
        aug_tel = None
        if name:
            aug_tel = Telescope(str(name), geo_x, geo_y, geo_z)
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
        # so far, these are the only options that are coming in from the
        # config files - may need to add more as more types are experienced
        if from_value == 'false':
            result = False
        elif from_value == 'true':
            result = True
        return result

    def _try_energy_with_blueprint(self, chunk, index):
        """
        A mechanism to augment the Energy WCS completely from the blueprint.
        Do nothing if the WCS information cannot be correctly created.

        :param chunk: The chunk to modify with the addition of energy
            information.
        :param index: The index in the blueprint for looking up plan
            information.
        """
        self.logger.debug('Begin augmentation with blueprint for energy.')
        aug_naxis = self._get_naxis('energy', index)

        specsys = _to_str(self._get_from_list('Chunk.energy.specsys', index))
        if aug_naxis is None:
            self.logger.debug('No blueprint energy information.')
        else:
            if not chunk.energy:
                chunk.energy = SpectralWCS(aug_naxis, specsys)
            else:
                chunk.energy.naxis = aug_naxis
                chunk.energy.specsys = specsys

            if chunk.energy is not None:
                chunk.energy.ssysobs = self._get_from_list(
                    'Chunk.energy.ssysobs', index)
                chunk.energy.restfrq = self._get_from_list(
                    'Chunk.energy.restfrq', index)
                chunk.energy.restwav = self._get_from_list(
                    'Chunk.energy.restwav', index)
                chunk.energy.velosys = self._get_from_list(
                    'Chunk.energy.velosys', index)
                chunk.energy.zsource = self._get_from_list(
                    'Chunk.energy.zsource', index)
                chunk.energy.ssyssrc = self._get_from_list(
                    'Chunk.energy.ssyssrc', index)
                chunk.energy.velang = self._get_from_list(
                    'Chunk.energy.velang', index)
                chunk.energy.bandpass_name = self._get_from_list(
                    'Chunk.energy.bandpassName', index)
                chunk.energy.transition = self._get_from_list(
                    'Chunk.energy.transition', index)
                chunk.energy.resolving_power = _to_float(self._get_from_list(
                    'Chunk.energy.resolvingPower', index))
        self.logger.debug('End augmentation with blueprint for energy.')

    def _try_observable_with_blueprint(self, chunk, index):
        """
        A mechanism to augment the Observable WCS completely from the
        blueprint. Do nothing if the WCS information cannot be correctly
        created.

        :param chunk: The chunk to modify with the addition of observable
            information.
        :param index: The index in the blueprint for looking up plan
            information.
        """
        self.logger.debug('Begin augmentation with blueprint for '
                          'observable.')
        chunk.observable_axis = _to_int(
            self._get_from_list('Chunk.observableAxis', index))
        aug_axis = self._two_param_constructor(
            'Chunk.observable.dependent.axis.ctype',
            'Chunk.observable.dependent.axis.cunit', index, _to_str, Axis)
        aug_bin = _to_int(
            self._get_from_list('Chunk.observable.dependent.bin', index))
        if aug_axis is not None and aug_bin is not None:
            chunk.observable = ObservableAxis(Slice(aug_axis, aug_bin))
        self.logger.debug('End augmentation with blueprint for polarization.')

    def _try_polarization_with_blueprint(self, chunk, index):
        """
        A mechanism to augment the Polarization WCS completely from the
        blueprint. Do nothing if the WCS information cannot be correctly
        created.

        :param chunk: The chunk to modify with the addition of polarization
            information.
        :param index: The index in the blueprint for looking up plan
            information.
        """
        self.logger.debug('Begin augmentation with blueprint for '
                          'polarization.')
        chunk.polarization_axis = _to_int(
            self._get_from_list('Chunk.polarizationAxis', index))
        aug_naxis = self._get_naxis('polarization', index)
        if aug_naxis is not None:
            if chunk.polarization:
                chunk.polarization.naxis = aug_naxis
            else:
                chunk.polarization = PolarizationWCS(aug_naxis)
                self.logger.debug(
                    f'Creating PolarizationWCS for {self.uri} from blueprint')

        self.logger.debug('End augmentation with blueprint for polarization.')

    def _try_position_range(self, chunk, index):
        self.logger.debug('Try to set the range for position from blueprint')
        if (self.blueprint._pos_axes_configed and chunk.position is not None
                and chunk.position.axis is not None):
            aug_range_c1_start = self._two_param_constructor(
                'Chunk.position.axis.range.start.coord1.pix',
                'Chunk.position.axis.range.start.coord1.val',
                index, _to_float, RefCoord)
            aug_range_c1_end = self._two_param_constructor(
                'Chunk.position.axis.range.end.coord1.pix',
                'Chunk.position.axis.range.end.coord1.val',
                index, _to_float, RefCoord)
            aug_range_c2_start = self._two_param_constructor(
                'Chunk.position.axis.range.start.coord2.pix',
                'Chunk.position.axis.range.start.coord2.val',
                index, _to_float, RefCoord)
            aug_range_c2_end = self._two_param_constructor(
                'Chunk.position.axis.range.end.coord2.pix',
                'Chunk.position.axis.range.end.coord2.val',
                index, _to_float, RefCoord)
            if (aug_range_c1_start and aug_range_c1_end and aug_range_c2_start
                    and aug_range_c2_end):
                chunk.position.axis.range = CoordRange2D(
                    Coord2D(aug_range_c1_start, aug_range_c1_end),
                    Coord2D(aug_range_c2_start, aug_range_c2_end))
                self.logger.debug('Completed setting range for position')

    def _try_position_with_blueprint(self, chunk, index):
        """
        A mechanism to augment the Position WCS completely from the blueprint.
        Do nothing if the WCS information cannot be correctly created.

        :param chunk: The chunk to modify with the addition of position
            information.
        :param index: The index in the blueprint for looking up plan
            information.
        """
        self.logger.debug('Begin augmentation with blueprint for position.')

        aug_x_axis = self._two_param_constructor(
            'Chunk.position.axis.axis1.ctype',
            'Chunk.position.axis.axis1.cunit', index, _to_str, Axis)
        aug_y_axis = self._two_param_constructor(
            'Chunk.position.axis.axis2.ctype',
            'Chunk.position.axis.axis2.cunit', index, _to_str, Axis)
        aug_x_error = self._two_param_constructor(
            'Chunk.position.axis.error1.syser',
            'Chunk.position.axis.error1.rnder', index, _to_float, CoordError)
        aug_y_error = self._two_param_constructor(
            'Chunk.position.axis.error2.syser',
            'Chunk.position.axis.error2.rnder', index, _to_float, CoordError)
        aug_dimension = self._two_param_constructor(
            'Chunk.position.axis.function.dimension.naxis1',
            'Chunk.position.axis.function.dimension.naxis2',
            index, _to_int, Dimension2D)
        aug_x_ref_coord = self._two_param_constructor(
            'Chunk.position.axis.function.refCoord.coord1.pix',
            'Chunk.position.axis.function.refCoord.coord1.val',
            index, _to_float, RefCoord)
        aug_y_ref_coord = self._two_param_constructor(
            'Chunk.position.axis.function.refCoord.coord2.pix',
            'Chunk.position.axis.function.refCoord.coord2.val',
            index, _to_float, RefCoord)
        aug_cd11 = _to_float(self._get_from_list(
            'Chunk.position.axis.function.cd11', index))
        aug_cd12 = _to_float(self._get_from_list(
            'Chunk.position.axis.function.cd12', index))
        aug_cd21 = _to_float(self._get_from_list(
            'Chunk.position.axis.function.cd21', index))
        aug_cd22 = _to_float(self._get_from_list(
            'Chunk.position.axis.function.cd22', index))

        aug_ref_coord = None
        if aug_x_ref_coord is not None and aug_y_ref_coord is not None:
            aug_ref_coord = Coord2D(aug_x_ref_coord, aug_y_ref_coord)
            self.logger.debug(
                f'Creating position Coord2D for {self.uri}')

        aug_function = None
        if (aug_dimension is not None and aug_ref_coord is not None and
            aug_cd11 is not None and aug_cd12 is not None and
                aug_cd21 is not None and aug_cd22 is not None):
            aug_function = CoordFunction2D(aug_dimension, aug_ref_coord,
                                           aug_cd11, aug_cd12, aug_cd21,
                                           aug_cd22)
            self.logger.debug(
                f'Creating position CoordFunction2D for {self.uri}')

        aug_axis = None
        if (aug_x_axis is not None and aug_y_axis is not None and
                aug_function is not None):
            aug_axis = CoordAxis2D(aug_x_axis, aug_y_axis, aug_x_error,
                                   aug_y_error, None, None, aug_function)
            self.logger.debug(
                f'Creating position CoordAxis2D for {self.uri}')

        if aug_axis is not None:
            if chunk.position:
                chunk.position.axis = aug_axis
            else:
                chunk.position = SpatialWCS(aug_axis)

        if chunk.position:
            chunk.position.coordsys = self._get_from_list(
                'Chunk.position.coordsys', index)
            chunk.position.equinox = _to_float(self._get_from_list(
                'Chunk.position.equinox', index))
            chunk.position.resolution = self._get_from_list(
                'Chunk.position.resolution', index)
        self.logger.debug('End augmentation with blueprint for position.')

    def _try_range(self, wcs, index, lookup):
        self.logger.debug(f'Try to set the range for {lookup}')
        aug_range_start = self._two_param_constructor(
            f'Chunk.{lookup}.axis.range.start.pix',
            f'Chunk.{lookup}.axis.range.start.val',
            index, _to_float, RefCoord)
        aug_range_end = self._two_param_constructor(
            f'Chunk.{lookup}.axis.range.end.pix',
            f'Chunk.{lookup}.axis.range.end.val',
            index, _to_float, RefCoord)
        if aug_range_start and aug_range_end:
            wcs.axis.range = CoordRange1D(aug_range_start, aug_range_end)
            self.logger.debug(f'Completed setting range for {lookup}')

    def _try_range_return(self, index, lookup):
        self.logger.debug(f'Try to set the range for {lookup}')
        range = None
        aug_range_start = self._two_param_constructor(
            f'Chunk.{lookup}.axis.range.start.pix',
            f'Chunk.{lookup}.axis.range.start.val',
            index, _to_float, RefCoord)
        aug_range_end = self._two_param_constructor(
            f'Chunk.{lookup}.axis.range.end.pix',
            f'Chunk.{lookup}.axis.range.end.val',
            index, _to_float, RefCoord)
        if aug_range_start and aug_range_end:
            range = CoordRange1D(aug_range_start, aug_range_end)
        self.logger.debug(f'Completed setting range for {lookup}')
        return range

    def _try_range_with_blueprint(self, chunk, index):
        """Use the blueprint to set elements and attributes that
        are not in the scope of astropy and files content, and therefore are
        not covered by the *WcsParser classes. Per PD 19/04/18, bounds and
        range are not covered by WCS keywords."""

        for i in ['energy', 'time', 'polarization']:
            axis_configed = getattr(self.blueprint,
                                    f'_{i}_axis_configed')
            if axis_configed:
                wcs = getattr(chunk, i)
                if wcs is not None and wcs.axis is not None:
                    if wcs.axis.range is None:
                        self._try_range(wcs, index, i)
        self._try_position_range(chunk, index)

    def _try_time_with_blueprint(self, chunk, index):
        """
        A mechanism to augment the Time WCS completely from the blueprint.
        Do nothing if the WCS information cannot be correctly created.

        :param chunk: The chunk to modify with the addition of time
            information.
        :param index: The index in the blueprint for looking up plan
            information.
        """
        self.logger.debug('Begin augmentation with blueprint for temporal.')

        chunk.time_axis = _to_int(self._get_from_list('Chunk.timeAxis', index))
        aug_naxis = self._get_naxis('time', index)
        if aug_naxis is not None:
            if chunk.time:
                chunk.time.naxis = aug_naxis
            else:
                chunk.time = TemporalWCS(aug_naxis)
                self.logger.debug('Creating TemporalWCS for {} from blueprint'.
                                  format(self.uri))
        if chunk.time is not None:
            chunk.time.exposure = _to_float(
                self._get_from_list('Chunk.time.exposure', index))
            chunk.time.resolution = _to_float(
                self._get_from_list('Chunk.time.resolution', index))
            chunk.time.timesys = _to_str(
                self._get_from_list('Chunk.time.timesys', index))
            chunk.time.trefpos = self._get_from_list('Chunk.time.trefpos',
                                                     index)
            chunk.time.mjdref = self._get_from_list('Chunk.time.mjdref', index)

        self.logger.debug('End augmentation with blueprint for temporal.')

    def _two_param_constructor(self, lookup1, lookup2, index, to_type, ctor):
        """
        Helper function to build from the blueprint, a CAOM2 entity that
        has two required parameters.

        :param lookup1: Blueprint lookup text for the first constructor
            parameter.
        :param lookup2: Blueprint lookup text for the second constructor
            parameter.
        :param index:  Which index in the blueprint to do the lookup on.
        :param to_type: Function to cast the blueprint value to a particular
            type.
        :param ctor: The constructor that has two parameters to build.
        :return: The instance returned by the constructor, or None if any of
            the values are undefined.
        """
        param1 = to_type(self._get_from_list(lookup1, index))
        param2 = to_type(self._get_from_list(lookup2, index))
        new_object = None
        if param1 is not None and param2 is not None:
            new_object = ctor(param1, param2)
        return new_object

    @staticmethod
    def _add_keywords(keywords, current, to_set):
        """
        Common code for adding keywords to a CAOM2 entity, capturing all
        the weird metadata cases that happen at CADC.

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
    Parses a FITS file and extracts the CAOM2 related information which can
    be used to augment an existing CAOM2 observation, plane or artifact. The
    constructor takes either a FITS file as argument or a list of dictionaries
    (FITS keyword=value) corresponding to each extension.

    The WCS-related keywords of the FITS file are consumed by the astropy.wcs
    package which might display warnings with regards to compliance.

    Example 1:
    parser = FitsParser(input = '/staging/700000o.fits.gz')
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

    headers = [] # list of dictionaries headers
    # populate headers
    parser = FitsParser(input=headers)

    parser.augment_observation(obs)
    ...

    """

    def __init__(self, src, obs_blueprint=None, uri=None):
        """
        Ctor
        :param src: List of headers (dictionary of FITS keywords:value) with
        one header for each extension or a FITS input file.
        :param obs_blueprint: externally provided blueprint
        :param uri: which artifact augmentation is based on
        """
        self.logger = logging.getLogger(__name__)
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
        self.apply_blueprint()

    @property
    def headers(self):
        """
        List of headers where each header should allow dictionary like
        access to the FITS attribute in that header
        :return:
        """
        return self._headers

    def ignore_chunks(self, artifact, index):
        # there is one Part per extension, the name is the extension number
        if (
            FitsParser._has_data_array(self._headers[index])
            and self.blueprint.has_chunk(index)
        ):
            if str(index) not in artifact.parts.keys():
                # TODO use extension name?
                artifact.parts.add(Part(str(index)))
                self.logger.debug(f'Part created for HDU {index}.')
            result = False
        else:
            artifact.parts.add(Part(str(index)))
            self.logger.debug(f'Create empty part for HDU {index}')
            result = True
        return result

    def apply_blueprint(self):

        # pointers that are short to type
        exts = self.blueprint._extensions
        wcs_std = self.blueprint._wcs_std
        plan = self.blueprint._plan

        # firstly, apply the functions
        if (self.blueprint._module is not None or
                self.blueprint._module_instance is not None):
            for key, value in plan.items():
                if ObsBlueprint.is_function(value):
                    if self._blueprint._module_instance is None:
                        plan[key] = self._execute_external(value, key, 0)
                    else:
                        plan[key] = self._execute_external_instance(
                            value, key, 0)
            for extension in exts:
                for key, value in exts[extension].items():
                    if ObsBlueprint.is_function(value):
                        if self._blueprint._module_instance is None:
                            exts[extension][key] = self._execute_external(
                                value, key, extension)
                        else:
                            exts[extension][key] = \
                                self._execute_external_instance(
                                    value, key, extension)

        # apply overrides from blueprint to all extensions
        for key, value in plan.items():
            if key in wcs_std:
                if ObsBlueprint.needs_lookup(value):
                    # alternative attributes provided for standard wcs attrib.
                    for header in self.headers:
                        for v in value[0]:
                            if v in header and \
                                    v not in wcs_std[key].split(','):
                                keywords = wcs_std[key].split(',')
                                for keyword in keywords:
                                    _set_by_type(header, keyword,
                                                 str(header[v]))
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
                logging.error('More extensions configured {} than headers '
                              '{}'.format(extension, len(self.headers)))
                continue
            hdr = self.headers[extension]
            for key, value in exts[extension].items():
                if ObsBlueprint.is_table(value):
                    continue
                keywords = wcs_std[key].split(',')
                for keyword in keywords:
                    _set_by_type(hdr, keyword, value)
                    logging.debug(
                        '{}: set to {} in extension {}'.format(keyword, value,
                                                               extension))
        # apply defaults to all extensions
        for key, value in plan.items():
            if ObsBlueprint.has_default_value(value):
                for index, header in enumerate(self.headers):
                    for keywords in value[0]:
                        for keyword in keywords.split(','):
                            if (not header.get(keyword.strip()) and
                                keyword == keywords and  # checking a string
                                    keywords == value[0][-1]):  # last item
                                # apply a default if a value does not already
                                # exist, and all possible values of
                                # keywords have been checked
                                _set_by_type(header, keyword.strip(), value[1])
                                logging.debug(
                                    '{}: set default value of {} in HDU {}.'.
                                    format(keyword, value[1], index))

        # TODO wcs in astropy ignores cdelt attributes when it finds a cd
        # attribute even if it's in a different axis
        for header in self.headers:
            cd_present = False
            for i in range(1, 6):
                if 'CD{0}_{0}'.format(i) in header:
                    cd_present = True
                    break
            if cd_present:
                for i in range(1, 6):
                    if f'CDELT{i}' in header and \
                            'CD{0}_{0}'.format(i) not in header:
                        header['CD{0}_{0}'.format(i)] = \
                            header[f'CDELT{i}']

        # TODO When a projection is specified, wcslib expects corresponding
        # DP arguments with NAXES attributes. Normally, omitting the attribute
        # signals no distortion which is the assumption in caom2blueprint for
        # energy and polarization axes. Following is a workaround for
        # SIP projections.
        # For more details see:
        # http://www.atnf.csiro.au/people/mcalabre/WCS/dcs_20040422.pdf
        for header in self.headers:
            sip = False
            for i in range(1, 6):
                if ((f'CTYPE{i}' in header) and
                        isinstance(header[f'CTYPE{i}'], str) and
                        ('-SIP' in header[f'CTYPE{i}'])):
                    sip = True
                    break
            if sip:
                for i in range(1, 6):
                    if (f'CTYPE{i}' in header) and \
                            ('-SIP' not in header[f'CTYPE{i}']) and \
                            (f'DP{i}' not in header):
                        header[f'DP{i}'] = 'NAXES: 1'

        return

    def augment_artifact(self, artifact, index=0):
        """
        Augments a given CAOM2 artifact with available FITS information
        :param artifact: existing CAOM2 artifact to be augmented
        """
        self.logger.debug(
            'Begin artifact augmentation for {} with {} HDUs.'.format(
                artifact.uri, len(self.headers)))

        if self.blueprint.get_configed_axes_count() == 0:
            raise TypeError(
                'No WCS Data. End artifact augmentation for {}.'.format(
                    artifact.uri))

        for i, header in enumerate(self.headers):
            if self.ignore_chunks(artifact, i):
                # artifact-level attributes still require updating
                BlueprintParser.augment_artifact(self, artifact, 0)
                continue
            self._wcs_parser = FitsWcsParser(header, self.file, str(i))
            super().augment_artifact(artifact, i)

        self.logger.debug(
            f'End artifact augmentation for {artifact.uri}.')

    def _get_chunk_naxis(self, chunk, index=None):
        # NOTE: astropy.wcs does not distinguished between WCS axes and
        # data array axes. naxis in astropy.wcs represents in fact the
        # number of WCS axes, whereas chunk.axis represents the naxis
        # of the data array. Solution is to determine it directly from
        # the header
        if 'ZNAXIS' in self._headers[index]:
            chunk.naxis = _to_int(self._headers[index]['ZNAXIS'])
        elif 'NAXIS' in self._headers[index]:
            chunk.naxis = _to_int(self._headers[index]['NAXIS'])
        else:
            super()._get_chunk_naxis(chunk)

    def _get_from_list(self, lookup, index, current=None):
        value = None
        try:
            keys = self.blueprint._get(lookup)
        except KeyError:
            self.add_error(lookup, sys.exc_info()[1])
            self.logger.debug(
                f'Could not find {lookup!r} in caom2blueprint configuration.')
            if current:
                self.logger.debug(
                    f'{lookup}: using current value of {current!r}.')
                value = current
            return value

        if ObsBlueprint.needs_lookup(keys):
            for ii in keys[0]:
                try:
                    value = self.headers[index].get(ii)
                    if value:
                        self.logger.debug(
                            f'{lookup}: assigned value {value} based on '
                            f'keyword {ii}.')
                        break
                except (KeyError, IndexError):
                    if keys[0].index(ii) == len(keys[0]) - 1:
                        self.add_error(lookup, sys.exc_info()[1])
                    # assign a default value, if one exists
                    if keys[1]:
                        if current is None:
                            value = keys[1]
                            self.logger.debug(
                                f'{lookup}: assigned default value {value}.')
                        else:
                            value = current
            if value is None:
                # checking current does not work in the general case,
                # because current might legitimately be 'None'
                if self._blueprint.update:
                    if (
                        current is not None
                        or (current is None and isinstance(value, bool))
                    ):
                        value = current
                        self.logger.debug(
                            f'{lookup}: used current value {value}.')
                else:
                    # assign a default value, if one exists
                    if keys[1]:
                        if current is None:
                            value = keys[1]
                            self.logger.debug(
                                f'{lookup}: assigned default value {value}.')
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

        This is a straight FITS BINTABLE lookup. There is no support for
        default values. Unless someone provides a compelling use case.

        :param lookup: where to find the column name
        :param extension: which extension
        :return: A string, which is a space-delimited list of all the values.
        """
        value = ''
        try:
            keywords = self.blueprint._get(lookup, extension)
        except KeyError as e:
            self.add_error(lookup, sys.exc_info()[1])
            self.logger.debug(
                'Could not find {!r} in fits2caom2 configuration.'.format(
                    lookup))
            raise e

        if isinstance(keywords, tuple) and keywords[0] == 'BINTABLE':

            # BINTABLE, so need to retrieve the data from the file
            if self.file is not None and self.file != '':
                with fits.open(self.file) as fits_data:
                    if fits_data[extension].header['XTENSION'] != 'BINTABLE':
                        raise ValueError(
                            'Got {} when looking for a BINTABLE '
                            'extension.'.format(
                                fits_data[extension].header['XTENSION']))
                    for ii in keywords[1]:
                        for jj in fits_data[extension].data[keywords[2]][ii]:
                            value = f'{jj} {value}'

        self.logger.debug(f'{lookup}: value is {value}')
        return value

    def _get_set_from_list(self, lookup, index):
        value = None
        keywords = None
        try:
            keywords = self.blueprint._get(lookup)
        except KeyError:
            self.add_error(lookup, sys.exc_info()[1])
            self.logger.debug(f'Could not find \'{lookup}\' in caom2blueprint '
                              f'configuration.')

        if isinstance(keywords, tuple):
            for ii in keywords[0]:
                try:
                    value = self.headers[index].get(ii)
                    break
                except KeyError:
                    self.add_error(lookup, sys.exc_info()[1])
                    if keywords[1]:
                        value = keywords[1]
                        self.logger.debug(
                            '{}: assigned default value {}.'.format(lookup,
                                                                    value))
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
    Parses an HDF5 file and extracts the CAOM2 related information which
    can be used to augment an existing CAOM2 observation, plane, or artifact.

    If there is per-Chunk metadata in the file, the constructor parameter
    'find_roots_here' is the address location in the file where the N Chunk
    metadata starts.

    The WCS-related keywords of the HDF5 files are used to create instances of
    astropy.wcs.WCS so that verify might be called.

    There is no CADC support for the equivalent of the FITS --fhead parameter
    for HDF5 files, which is why the name of the file on a local disk is
    required.

    How the classes work together for HDF5 files:
    - build an HDF5ObsBlueprint, with _CAOM2_ELEMENT keys, and HDF5 metadata
        path names as keys
    - cache the metadata from an HDF5 file in the HDF5ObsBlueprint. This
        caching is done in the "apply_blueprint_from_file" method in the
        Hdf5Parser class, and replaces the path names in the blueprint with
        the values from the HDF5 file. The caching is done so that all HDF5
        file access is isolated to one point in time.
    - use the cached metadata to build astropy.wcs instances for verification
        in Hdf5WcsParser.
    - use the astropy.wcs instance and other blueprint metadata to fill the
        CAOM2 record.
    """

    def __init__(
        self, obs_blueprint, uri, h5_file, find_roots_here='sitedata'
    ):
        """
        :param obs_blueprint: Hdf5ObsBlueprint instance
        :param uri: which artifact augmentation is based on
        :param h5_file: h5py file handle
        :param find_roots_here: str location where Chunk metadata starts
        """
        self._file = h5_file
        # where N Chunk metadata starts
        self._find_roots_here = find_roots_here
        # the length of the array is the number of Parts in an HDF5 file,
        # and the values are HDF5 lookup path names.
        self._extension_names = []
        super().__init__(obs_blueprint, uri)
        # used to set the astropy wcs info, resulting in a validated wcs
        # that can be used to construct a valid CAOM2 record
        self._wcs_parser = None

    def apply_blueprint_from_file(self):
        """
        Retrieve metadata from file, cache in the blueprint.
        """
        self.logger.debug('Begin apply_blueprint_from_file')
        # h5py is an extra in this package since most collections do not
        # require it
        import h5py
        individual, multi, attributes = self._extract_path_names_from_blueprint()

        def _extract_from_item(name, object):
            """
            Function signature dictated by h5py visititems implementation.
            Executed for each dataset/group in an HDF5 file.

            :param name: fully-qualified HDF5 path name
            :param object: what the HDF5 path name points to
            """
            if name == self._find_roots_here:
                for ii, path_name in enumerate(object.keys()):
                    # store the names and locations of the Part/Chunk metadata
                    temp = f'{name}/{path_name}'
                    self.logger.debug(f'Adding extension {temp}')
                    self._extension_names.append(temp)
                    self._blueprint._extensions[ii] = {}

            # If it's the Part/Chunk metadata, capture it to extensions.
            # Syntax of the keys described in Hdf5ObsBlueprint class.
            for part_index, part_name in enumerate(self._extension_names):
                if (
                    name.startswith(part_name)
                    and isinstance(object, h5py.Dataset)
                    and object.dtype.names is not None
                ):
                    for d_name in object.dtype.names:
                        temp_path = f'{name.replace(part_name, "")}/{d_name}'
                        for path_name in multi.keys():
                            if path_name == temp_path:
                                for jj in multi.get(path_name):
                                    self._blueprint.set(
                                        jj, object[d_name], part_index
                                    )
                            elif (path_name.startswith(temp_path)
                                  and '(' in path_name):
                                z = path_name.split('(')
                                if ':' in z[1]:
                                    a = z[1].split(')')[0].split(':')
                                    if len(a) > 2:
                                        raise NotImplementedError
                                    for jj in multi.get(path_name):
                                        self._blueprint.set(
                                            jj,
                                            object[d_name][int(a[0])][
                                                int(a[1])],
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

            # if it's Observation/Plane/Artifact metadata, capture it to
            # the base blueprint
            if isinstance(object, h5py.Dataset):
                if object.dtype.names is not None:
                    for d_name in object.dtype.names:
                        temp = f'//{name}/{d_name}'
                        if temp in individual.keys():
                            for jj in individual.get(temp):
                                self._blueprint.set(jj, object[d_name], 0)

        if len(individual) == 0 and len(multi) == 0:
            self._extract_from_attrs(attributes)
        else:
            self._file.visititems(_extract_from_item)
        self.logger.debug('Done apply_blueprint_from_file')

    def _extract_from_attrs(self, attributes):
        # I don't currently see any way to have more than one Part, if relying on
        # attrs for metadata
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
        """
        individual = defaultdict(list)
        multi = defaultdict(list)
        attributes = defaultdict(list)
        for key, value in self._blueprint._plan.items():
            if ObsBlueprint.needs_lookup(value):
                for ii in value[0]:
                    if ii.startswith('//'):
                        individual[ii].append(key)
                    elif ii.startswith('/'):
                        multi[ii].append(key)
                    else:
                        attributes[ii].append(key)
        return individual, multi, attributes

    def apply_blueprint(self):
        self.logger.debug('Begin apply_blueprint')
        self.apply_blueprint_from_file()

        # after the apply_blueprint_from_file call, all the metadata from the
        # file has been applied to the blueprint, so now do the bits that
        # require no access to file content

        # pointers that are short to type
        exts = self._blueprint._extensions
        plan = self._blueprint._plan

        # apply the functions
        if (self._blueprint._module is not None or
                self._blueprint._module_instance is not None):
            for key, value in plan.items():
                if ObsBlueprint.is_function(value):
                    if self._blueprint._module_instance is None:
                        plan[key] = self._execute_external(value, key, 0)
                    else:
                        plan[key] = self._execute_external_instance(
                            value, key, 0)
            for extension in exts:
                for key, value in exts[extension].items():
                    if ObsBlueprint.is_function(value):
                        if self._blueprint._module_instance is None:
                            exts[extension][key] = self._execute_external(
                                value, key, extension)
                        else:
                            exts[extension][key] = \
                                self._execute_external_instance(
                                    value, key, extension)

        # blueprint already contains all the overrides, only need to make
        # sure the overrides get applied to all the extensions
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
                self.logger.debug(
                    f'{key}: set to {value} in extension {extension}')

        # if no values have been set by file lookups, function execution,
        # or applying overrides, apply defaults, including to all extensions
        for key, value in plan.items():
            if ObsBlueprint.needs_lookup(value) and value[1]:
                # there is a default value in the blueprint that can be used
                for extension in exts:
                    q = exts[extension].get(key)
                    if q is None:
                        exts[extension][key] = value[1]
                        self.logger.debug(
                            f'Add {key} and assign default value of '
                            f'{value[1]} in extension {extension}.')
                    elif ObsBlueprint.needs_lookup(value):
                        exts[extension][key] = value[1]
                        self.logger.debug(
                            f'{key}: set value to default of {value[1]} in '
                            f'extension {extension}.')
                plan[key] = value[1]
                self.logger.debug(f'{key}: set value to default of {value[1]}')

        self._file.close()
        self.logger.debug('Done apply_blueprint')
        return

    def augment_artifact(self, artifact, index=0):
        self._wcs_parser = Hdf5WcsParser(self._blueprint, 0)
        super().augment_artifact(artifact, 0)
        for ii in range(1, len(self._blueprint._extensions)):
            self._wcs_parser = Hdf5WcsParser(self._blueprint, ii)
            super().augment_artifact(artifact, ii)

    def _get_chunk_naxis(self, chunk, index):
        chunk.naxis = self._get_from_list('Chunk.naxis', index, chunk.naxis)

    def ignore_chunks(self, artifact, index=0):
        artifact.parts.add(Part(str(index)))
        return False


class WcsParser:
    """
    WCS axes methods.
    """

    ENERGY_AXIS = 'energy'
    POLARIZATION_AXIS = 'polarization'
    TIME_AXIS = 'time'

    def __init__(self):
        self.logger = logging.getLogger(self.__class__.__name__)
        self.wcs = None

    def augment_custom(self, chunk):
        """
        Augments a chunk with custom WCS information
        :param chunk:
        :return:
        """
        self.logger.debug('Begin Custom WCS augmentation.')
        if chunk is None or not isinstance(chunk, Chunk):
            raise ValueError(f'Chunk type mis-match for {chunk}.')

        custom_axis_index = self._get_axis_index(CUSTOM_CTYPES)
        if custom_axis_index is None:
            self.logger.debug('No WCS Custom info')
            return

        chunk.custom_axis = custom_axis_index + 1

        naxis = CoordAxis1D(self._get_axis(custom_axis_index))
        if self.wcs.has_cd():
            delta = self.wcs.cd[custom_axis_index][
                custom_axis_index]
        else:
            delta = self.wcs.cdelt[custom_axis_index]
        naxis.function = CoordFunction1D(
            self._get_axis_length(custom_axis_index + 1),
            delta,
            self._get_ref_coord(custom_axis_index))
        if not chunk.custom:
            chunk.custom = CustomWCS(naxis)
        else:
            chunk.custom.axis = naxis

        self.logger.debug('End Custom WCS augmentation.')

    def augment_energy(self, chunk):
        """
        Augments the energy information in a chunk
        :param chunk:
        """
        self.logger.debug('Begin Energy WCS augmentation.')
        if chunk is None or not isinstance(chunk, Chunk):
            raise ValueError(f'Chunk type mis-match for {chunk}.')

        # get the energy axis
        energy_axis_index = self._get_axis_index(ENERGY_CTYPES)

        if energy_axis_index is None:
            self.logger.debug('No WCS Energy info.')
            return

        chunk.energy_axis = energy_axis_index + 1
        naxis = CoordAxis1D(self._get_axis(energy_axis_index))
        naxis.error = self._get_coord_error(energy_axis_index)
        if self.wcs.has_cd():
            delta = self.wcs.cd[energy_axis_index][energy_axis_index]
        else:
            delta = self.wcs.cdelt[energy_axis_index]
        naxis.function = CoordFunction1D(
            self._get_axis_length(energy_axis_index + 1), delta,
            self._get_ref_coord(energy_axis_index))

        specsys = _to_str(self.wcs.specsys)
        if not chunk.energy:
            chunk.energy = SpectralWCS(naxis, specsys)
        else:
            chunk.energy.axis = naxis
            chunk.energy.specsys = specsys

        chunk.energy.ssysobs = _to_str(self._sanitize(self.wcs.ssysobs))
        # wcs returns 0.0 by default
        if self._sanitize(self.wcs.restfrq) != 0:
            chunk.energy.restfrq = self._sanitize(self.wcs.restfrq)
        if self._sanitize(self.wcs.restwav) != 0:
            chunk.energy.restwav = self._sanitize(self.wcs.restwav)
        chunk.energy.velosys = self._sanitize(self.wcs.velosys)
        chunk.energy.zsource = self._sanitize(self.wcs.zsource)
        chunk.energy.ssyssrc = _to_str(self._sanitize(self.wcs.ssyssrc))
        chunk.energy.velang = self._sanitize(self.wcs.velangl)
        self.logger.debug('End Energy WCS augmentation.')

    def augment_position(self, chunk):
        """
        Augments a chunk with spatial WCS information
        :param chunk:
        :return:
        """
        self.logger.debug('Begin Spatial WCS augmentation.')
        if chunk is None or not isinstance(chunk, Chunk):
            raise ValueError(f'Chunk type mis-match for {chunk}.')

        position_axes_indices = self._get_position_axis()
        if not position_axes_indices:
            self.logger.debug('No Spatial WCS found')
            return

        chunk.position_axis_1 = position_axes_indices[0]
        chunk.position_axis_2 = position_axes_indices[1]
        axis = self._get_spatial_axis(chunk.position_axis_1 - 1,
                                      chunk.position_axis_2 - 1)
        if chunk.position:
            chunk.position.axis = axis
        else:
            chunk.position = SpatialWCS(axis)

        chunk.position.coordsys = _to_str(self._sanitize(self.wcs.radesys))
        temp = self._sanitize(self.wcs.equinox)
        if (temp is not None and 1800.0 <= temp <= 2500) or temp is None:
            chunk.position.equinox = temp
        self.logger.debug('End Spatial WCS augmentation.')

    def augment_temporal(self, chunk):
        """
        Augments a chunk with temporal WCS information

        :param chunk:
        :return:
        """
        self.logger.debug('Begin TemporalWCS augmentation.')
        if chunk is None or not isinstance(chunk, Chunk):
            raise ValueError(f'Chunk type mis-match for {chunk}.')

        time_axis_index = self._get_axis_index(TIME_KEYWORDS)

        if time_axis_index is None:
            self.logger.debug('No WCS Time info.')
            return

        chunk.time_axis = time_axis_index + 1
        # set chunk.time
        self.logger.debug('Begin temporal axis augmentation.')

        aug_naxis = self._get_axis(time_axis_index)
        aug_error = self._get_coord_error(time_axis_index)
        aug_ref_coord = self._get_ref_coord(time_axis_index)
        if self.wcs.has_cd():
            delta = self.wcs.cd[time_axis_index][time_axis_index]
        else:
            delta = self.wcs.cdelt[time_axis_index]
        axis_length = self._get_axis_length(time_axis_index + 1)
        if aug_ref_coord is not None and axis_length is not None:
            aug_function = CoordFunction1D(axis_length, delta, aug_ref_coord)
            naxis = CoordAxis1D(aug_naxis, aug_error, None, None, aug_function)
            if not chunk.time:
                chunk.time = TemporalWCS(naxis)
            else:
                chunk.time.axis = naxis

            self._finish_chunk_time(chunk)
        self.logger.debug('End TemporalWCS augmentation.')

    def augment_polarization(self, chunk):
        """
        Augments a chunk with polarization WCS information
        :param chunk:
        :return:
        """
        self.logger.debug('Begin Polarization WCS augmentation.')
        if chunk is None or not isinstance(chunk, Chunk):
            raise ValueError(f'Chunk type mis-match for {chunk}.')

        polarization_axis_index = self._get_axis_index(POLARIZATION_CTYPES)
        if polarization_axis_index is None:
            self.logger.debug('No WCS Polarization info')
            return

        chunk.polarization_axis = polarization_axis_index + 1

        naxis = CoordAxis1D(self._get_axis(polarization_axis_index))
        if self.wcs.has_cd():
            delta = self.wcs.cd[polarization_axis_index][
                polarization_axis_index]
        else:
            delta = self.wcs.cdelt[polarization_axis_index]
        naxis.function = CoordFunction1D(
            self._get_axis_length(polarization_axis_index + 1),
            delta,
            self._get_ref_coord(polarization_axis_index))
        if not chunk.polarization:
            chunk.polarization = PolarizationWCS(naxis)
        else:
            chunk.polarization.axis = naxis

        self.logger.debug('End Polarization WCS augmentation.')

    def augment_observable(self, chunk):
        """
        Augments a chunk with an observable axis.

        :param chunk:
        :return:
        """
        self.logger.debug('Begin Observable WCS augmentation.')
        if chunk is None or not isinstance(chunk, Chunk):
            raise ValueError(f'Chunk type mis-match for {chunk}.')

        observable_axis_index = self._get_axis_index(OBSERVABLE_CTYPES)
        if observable_axis_index is None:
            self.logger.debug('No Observable axis info')
            return

        chunk.observable_axis = observable_axis_index + 1
        self._finish_chunk_observable(chunk)
        self.logger.debug('End Observable WCS augmentation.')

    def _finish_chunk_time(self, chunk):
        raise NotImplementedError

    def _get_axis(self, index, over_ctype=None, over_cunit=None):
        """ Assemble a generic axis """
        aug_ctype = str(self.wcs.ctype[index]) if over_ctype is None \
            else over_ctype
        aug_cunit = str(self.wcs.cunit[index]) if over_cunit is None \
            else over_cunit
        if aug_cunit is not None and len(aug_cunit) == 0:
            aug_cunit = None
        aug_axis = Axis(aug_ctype, aug_cunit)
        return aug_axis

    def _get_axis_index(self, keywords):
        """
        Return the index of a specific axis type or None of it doesn't exist
        :param keywords:
        :return:
        """
        axis = None
        for i, elem in enumerate(self.wcs.ctype):
            elem = elem.split('-')[0]
            if elem in keywords:
                axis = i
                break
            elif len(elem) == 0:
                check = self.wcs.ctype[i]
                if check in keywords:
                    axis = i
                    break
        return axis

    def _get_axis_length(self, index):
        raise NotImplementedError

    def _get_cd(self, x_index, y_index):
        """ returns cd info"""

        try:
            if self.wcs.has_cd():
                cd11 = self.wcs.cd[x_index][x_index]
                cd12 = self.wcs.cd[x_index][y_index]
                cd21 = self.wcs.cd[y_index][x_index]
                cd22 = self.wcs.cd[y_index][y_index]
            else:
                cd11 = self.wcs.cdelt[x_index]
                cd12 = self.wcs.crota[x_index]
                cd21 = self.wcs.crota[y_index]
                cd22 = self.wcs.cdelt[y_index]
        except AttributeError:
            self.logger.debug(
                f'Error searching for CD* values {sys.exc_info()[1]}')
            cd11 = None
            cd12 = None
            cd21 = None
            cd22 = None

        return cd11, cd12, cd21, cd22

    def _get_coord_error(self, index):
        aug_coord_error = None
        aug_csyer = self._sanitize(self.wcs.csyer[index])
        aug_crder = self._sanitize(self.wcs.crder[index])
        if aug_csyer is not None and aug_crder is not None:
            aug_coord_error = CoordError(aug_csyer, aug_crder)
        return aug_coord_error

    def _get_dimension(self, xindex, yindex):
        aug_dimension = None
        aug_dim1 = _to_int(self._get_axis_length(xindex + 1))
        aug_dim2 = _to_int(self._get_axis_length(yindex + 1))
        if aug_dim1 and aug_dim2:
            aug_dimension = Dimension2D(aug_dim1, aug_dim2)
            self.logger.debug('End 2D dimension augmentation.')
        return aug_dimension

    def _get_position_axis(self):
        # there are two celestial axes, get the applicable indices from
        # the axis_types
        xindex = self._get_axis_index(POSITION_CTYPES[0])
        yindex = self._get_axis_index(POSITION_CTYPES[1])

        if (xindex is not None) and (yindex is not None):
            return xindex + 1, yindex + 1
        elif (xindex is None) and (yindex is None):
            return None
        else:
            raise ValueError('Found only one position axis ra/dec: {}/{} in '
                             '{}'.
                             format(xindex, yindex, self.file))

    def _get_ref_coord(self, index):
        aug_crpix = _to_float(self._sanitize(self.wcs.crpix[index]))
        aug_crval = _to_float(self._sanitize(self.wcs.crval[index]))
        aug_ref_coord = None
        if aug_crpix is not None and aug_crval is not None:
            aug_ref_coord = RefCoord(aug_crpix, aug_crval)
        return aug_ref_coord

    def _get_spatial_axis(self, xindex, yindex):
        """Assemble the bits to make the axis parameter needed for
        SpatialWCS construction."""
        aug_dimension = self._get_dimension(xindex, yindex)

        aug_ref_coord = Coord2D(self._get_ref_coord(xindex),
                                self._get_ref_coord(yindex))

        aug_cd11, aug_cd12, aug_cd21, aug_cd22 = \
            self._get_cd(xindex, yindex)

        if aug_dimension is not None and \
            aug_ref_coord is not None and \
            aug_cd11 is not None and \
            aug_cd12 is not None and \
            aug_cd21 is not None and \
                aug_cd22 is not None:
            aug_function = CoordFunction2D(aug_dimension, aug_ref_coord,
                                           aug_cd11, aug_cd12,
                                           aug_cd21, aug_cd22)
            self.logger.debug('End CoordFunction2D augmentation.')
        else:
            aug_function = None

        aug_axis = CoordAxis2D(self._get_axis(xindex),
                               self._get_axis(yindex),
                               self._get_coord_error(xindex),
                               self._get_coord_error(yindex),
                               None, None, aug_function)
        self.logger.debug('End CoordAxis2D augmentation.')
        return aug_axis

    def _sanitize(self, value):
        """
        Sanitizes values from content to caom2
        :param value:
        :return:
        """
        if value is None:
            return None
        elif isinstance(value, float) and math.isnan(value):
            return None
        elif not str(value):
            return None  # empty string
        else:
            return value


class FitsWcsParser(WcsParser):
    """
    Parser to augment chunks with positional, temporal, energy and polarization
    information based on the WCS keywords in an extension of a FITS header.

    Note: Under the hood, this class uses the astropy.wcs package to parse the
    header and any inconsistencies or missing keywords are reported back as
    warnings.
    """

    def __init__(self, header, file, extension):
        """

        :param header: FITS extension header
        :param file: name of FITS file
        :param extension: which HDU
        WCS axes methods of this class.
        """
        super().__init__()
        self.log_filter = HDULoggingFilter()
        self.log_filter.extension(extension)
        self.logger.addFilter(self.log_filter)
        logastro = logging.getLogger('astropy')
        logastro.addFilter(self.log_filter)
        logastro.propagate = False
        header_string = header.tostring().rstrip()
        header_string = header_string.replace('END' + ' ' * 77, '')
        self.wcs = Wcsprm(header_string.encode('ascii'))
        self.wcs.fix()
        self.header = header
        self.file = file
        self.extension = extension

    def _finish_chunk_observable(self, chunk):
        ctype = self.header.get(f'CTYPE{chunk.observable_axis}')
        cunit = self.header.get(f'CUNIT{chunk.observable_axis}')
        pix_bin = self.header.get(f'CRPIX{chunk.observable_axis}')
        if ctype is not None and cunit is not None and pix_bin is not None:
            chunk.observable = ObservableAxis(
                Slice(self._get_axis(0, ctype, cunit), pix_bin))

    def _finish_chunk_time(self, chunk):
        """
        The expected caom2 - FITS keywords mapping is:

        time.exposure = EXPTIME
        time.resolution = TIMEDEL
        time.timesys = TIMESYS default UTC
        time.trefpos = TREFPOS
        time.mjdref = MJDREF | MJDDATE
        """
        chunk.time.exposure = _to_float(self.header.get('EXPTIME'))
        chunk.time.resolution = _to_float(self.header.get('TIMEDEL'))
        chunk.time.timesys = str(self.header.get('TIMESYS', 'UTC'))
        chunk.time.trefpos = self.header.get('TREFPOS', None)
        chunk.time.mjdref = self.header.get('MJDREF',
                                            self.header.get('MJDDATE'))

    def _get_axis_length(self, for_axis):
        # try ZNAXIS first in order to get the size of the original
        # image in case it was FITS compressed
        result = _to_int(self._sanitize(
            self.header.get(f'ZNAXIS{for_axis}')))
        if result is None:
            result = _to_int(self._sanitize(
                self.header.get(f'NAXIS{for_axis}')))
        if result is None:
            msg = f'Could not find axis length for axis {for_axis}'
            raise ValueError(msg)
        return result


class Hdf5WcsParser(WcsParser):
    """
    This class initializes an astropy.wcs instance with metadata from an
    Hdf5ObsBlueprint populated using an Hdf5Parser.
    """

    def __init__(self, blueprint, extension):
        """
        :param blueprint: ObsBlueprint
        """
        super().__init__()
        self._wcs = None
        self._axes = {
            'ra': [0, False],
            'dec': [0, False],
            'time': [0, False],
            'energy': [0, False],
            'polarization': [0, False],
            'observable': [0, False],
            'custom': [0, False],
        }
        self._blueprint = blueprint
        # int - index into blueprint._plan extensions
        self._extension = extension
        self._set_wcs()

    @property
    def wcs(self):
        return self._wcs.wcs

    @wcs.setter
    def wcs(self, value):
        self._wcs = value

    def _get_axis_index(self, keywords):
        result = self._axes['custom'][0]
        if 'RA' in keywords:
            result = self._axes['ra'][0]
        elif 'DEC' in keywords:
            result = self._axes['dec'][0]
        elif 'TIME' in keywords:
            result = self._axes['time'][0]
        elif 'FREQ' in keywords:
            result = self._axes['energy'][0]
        elif 'STOKES' in keywords:
            result = self._axes['polarization'][0]
        elif 'FLUX' in keywords:
            result = self._axes['observable'][0]
        return result

    def _get_axis_length(self, for_axis):
        if self._wcs.array_shape is None:
            # TODO I think this is wrong
            return 1
        else:
            if len(self._wcs.array_shape) == 1:
                result = self._wcs.array_shape[0]
            else:
                result = self._wcs.array_shape[for_axis-1]
            return _to_int(result)

    def assign_sanitize(self, assignee, index, key, sanitize=True):
        """
        Do not want to blindly assign None to astropy.wcs attributes, so
        use this method for conditional assignment.

        The current implementation is that ff there is a legitimate need to
        assign None to a value, either use 'set' in the Hdf5ObsBlueprint, and
        specifically assign None, or execute a function to set it to None
        conditionally. There will be no support for a Default value of None
        with HDF5 files.
        """
        x = self._blueprint._get(key, self._extension)
        if sanitize:
            x = self._sanitize(x)
        if x is not None and not ObsBlueprint.needs_lookup(x):
            assignee[index] = x

    def _set_wcs(self):
        self._wcs = WCS(naxis=self._blueprint.get_configed_axes_count())
        array_shape = [0] * self._blueprint.get_configed_axes_count()
        count = 0
        if self._blueprint._pos_axes_configed:
            self._axes['ra'][1] = True
            self._axes['dec'][1] = True
            self._axes['ra'][0] = count
            self._axes['dec'][0] = count + 1
            temp = [0] * self._blueprint.get_configed_axes_count()
            cd = [temp.copy()
                  for ii in range(self._blueprint.get_configed_axes_count())]
            self.assign_sanitize(self._wcs.wcs.ctype, count,
                                 'Chunk.position.axis.axis1.ctype')
            self.assign_sanitize(self._wcs.wcs.ctype, count + 1,
                                 'Chunk.position.axis.axis2.ctype')
            self.assign_sanitize(self._wcs.wcs.cunit, count,
                                 'Chunk.position.axis.axis1.cunit')
            self.assign_sanitize(self._wcs.wcs.cunit, count + 1,
                                 'Chunk.position.axis.axis2.cunit')
            array_shape[count] = self._blueprint._get(
                'Chunk.position.axis.function.dimension.naxis1')
            array_shape[count + 1] = self._blueprint._get(
                'Chunk.position.axis.function.dimension.naxis2')
            self.assign_sanitize(
                self._wcs.wcs.crpix, count,
                'Chunk.position.axis.function.refCoord.coord1.pix')
            self.assign_sanitize(
                self._wcs.wcs.crpix, count + 1,
                'Chunk.position.axis.function.refCoord.coord2.pix')
            self.assign_sanitize(
                self._wcs.wcs.crval, count,
                'Chunk.position.axis.function.refCoord.coord1.val')
            self.assign_sanitize(
                self._wcs.wcs.crval, count + 1,
                'Chunk.position.axis.function.refCoord.coord2.val')
            x = self._blueprint._get('Chunk.position.axis.function.cd11',
                                     self._extension)
            if x is not None and not ObsBlueprint.needs_lookup(x):
                cd[count][0] = x
            x = self._blueprint._get('Chunk.position.axis.function.cd12',
                                     self._extension)
            if x is not None and not ObsBlueprint.needs_lookup(x):
                cd[count][1] = x
            x = self._blueprint._get('Chunk.position.axis.function.cd21',
                                     self._extension)
            if x is not None and not ObsBlueprint.needs_lookup(x):
                cd[count + 1][0] = x
            x = self._blueprint._get('Chunk.position.axis.function.cd22',
                                     self._extension)
            if x is not None and not ObsBlueprint.needs_lookup(x):
                cd[count + 1][1] = x
            self.assign_sanitize(self._wcs.wcs.crder, count,
                                 'Chunk.position.axis.error1.rnder')
            self.assign_sanitize(self._wcs.wcs.crder, count + 1,
                                 'Chunk.position.axis.error2.rnder')
            self.assign_sanitize(self._wcs.wcs.csyer, count,
                                 'Chunk.position.axis.error1.syser')
            self.assign_sanitize(self._wcs.wcs.csyer, count + 1,
                                 'Chunk.position.axis.error2.syser')
            self._finish_position()
            self._wcs.wcs.cd = cd
            count += 2
        if self._blueprint._time_axis_configed:
            self._axes['time'][1] = True
            self._axes['time'][0] = count
            self.assign_sanitize(self._wcs.wcs.ctype, count,
                                 'Chunk.time.axis.axis.ctype', False)
            self.assign_sanitize(self._wcs.wcs.cunit, count,
                                 'Chunk.time.axis.axis.cunit', False)
            array_shape[count] = self._blueprint._get(
                'Chunk.time.axis.function.naxis', self._extension)
            self.assign_sanitize(
                self._wcs.wcs.crpix, count,
                'Chunk.time.axis.function.refCoord.pix', False)
            self.assign_sanitize(
                self._wcs.wcs.crval, count,
                'Chunk.time.axis.function.refCoord.val', False)
            self.assign_sanitize(self._wcs.wcs.crder, count,
                                 'Chunk.time.axis.error.rnder')
            self.assign_sanitize(self._wcs.wcs.csyer, count,
                                 'Chunk.time.axis.error.syser')
            self._finish_time()
            count += 1
        if self._blueprint._energy_axis_configed:
            self._axes['energy'][1] = True
            self._axes['energy'][0] = count
            self.assign_sanitize(self._wcs.wcs.ctype, count,
                                 'Chunk.energy.axis.axis.ctype', False)
            self.assign_sanitize(self._wcs.wcs.cunit, count,
                                 'Chunk.energy.axis.axis.cunit', False)
            array_shape[count] = self._blueprint._get(
                'Chunk.energy.axis.function.naxis', self._extension)
            self.assign_sanitize(
                self._wcs.wcs.crpix, count,
                'Chunk.energy.axis.function.refCoord.pix', False)
            self.assign_sanitize(
                self._wcs.wcs.crval, count,
                'Chunk.energy.axis.function.refCoord.val', False)
            self.assign_sanitize(
                self._wcs.wcs.crder, count, 'Chunk.energy.axis.error.rnder')
            self.assign_sanitize(
                self._wcs.wcs.csyer, count, 'Chunk.energy.axis.error.syser')
            self._finish_energy()
            count += 1
        if self._blueprint._polarization_axis_configed:
            self._axes['polarization'][1] = True
            self._axes['polarization'][0] = count
            self.assign_sanitize(self._wcs.wcs.ctype, count,
                                 'Chunk.polarization.axis.axis.ctype', False)
            self.assign_sanitize(self._wcs.wcs.cunit, count,
                                 'Chunk.polarization.axis.axis.cunit', False)
            array_shape[count] = self._blueprint._get(
                'Chunk.polarization.axis.function.naxis', self._extension)
            self.assign_sanitize(
                self._wcs.wcs.crpix, count,
                'Chunk.polarization.axis.function.refCoord.pix', False)
            self.assign_sanitize(
                self._wcs.wcs.crval, count,
                'Chunk.polarization.axis.function.refCoord.val', False)
            count += 1
            # TODO - where's the delta?
        if self._blueprint._obs_axis_configed:
            self._axes['observable'][1] = True
            self._axes['observable'][0] = count
            self.assign_sanitize(self._wcs.wcs.ctype, count,
                                 'Chunk.observable.axis.axis.ctype', False)
            self.assign_sanitize(self._wcs.wcs.cunit, count,
                                 'Chunk.observable.axis.axis.cunit', False)
            array_shape[count] = 1.0
            self.assign_sanitize(self._wcs.wcs.crpix, count,
                                 'Chunk.observable.axis.function.refCoord.pix',
                                 False)
            self._wcs.wcs.crval[count] = 0.0
            count += 1
        if self._blueprint._custom_axis_configed:
            self._axes['custom'][1] = True
            self._axes['custom'][0] = count
            self.assign_sanitize(self._wcs.wcs.ctype, count,
                                 'Chunk.custom.axis.axis.ctype', False)
            self.assign_sanitize(self._wcs.wcs.cunit, count,
                                 'Chunk.custom.axis.axis.cunit', False)
            array_shape[count] = self._blueprint._get(
                'Chunk.custom.axis.function.naxis', self._extension)
            # TODO delta
            self.assign_sanitize(self._wcs.wcs.crpix, count,
                                 'Chunk.custom.axis.function.refCoord.pix',
                                 False)
            self.assign_sanitize(self._wcs.wcs.crval, count,
                                 'Chunk.custom.axis.function.refCoord.val',
                                 False)
            count += 1

        self._wcs.array_shape = array_shape

    def _finish_chunk_observable(self, chunk):
        ctype = self._wcs.wcs.ctype[chunk.observable_axis-1]
        cunit = self._wcs.wcs.ctype[chunk.observable_axis-1]
        pix_bin = _to_int(self._wcs.wcs.crpix[chunk.observable_axis-1])
        if ctype is not None and cunit is not None and pix_bin is not None:
            chunk.observable = ObservableAxis(
                Slice(self._get_axis(0, ctype, cunit), pix_bin))

    def _finish_chunk_time(self, chunk):
        if not math.isnan(self._wcs.wcs.xposure):
            chunk.time.exposure = self._wcs.wcs.xposure
        chunk.time.timesys = self._wcs.wcs.timesys
        if self._wcs.wcs.trefpos is not None and self._wcs.wcs.trefpos != '':
            chunk.time.trefpos = self._wcs.wcs.trefpos
        # convert from the numpy array length 2 of self._wcs.wcs.mjdref
        # to a single value
        # TODO chunk.time.mjdref = self._wcs.to_header().get('MJDREF')

    def _finish_energy(self):
        x = self._blueprint._get('Chunk.energy.specsys', self._extension)
        if x and not ObsBlueprint.needs_lookup(x):
            self._wcs.wcs.specsys = x
        x = self._blueprint._get('Chunk.energy.ssysobs', self._extension)
        if x and not ObsBlueprint.needs_lookup(x):
            self._wcs.wcs.ssysobs = x
        x = self._blueprint._get('Chunk.energy.restfrq', self._extension)
        if x and not ObsBlueprint.needs_lookup(x):
            self._wcs.wcs.restfrq = _to_float(x)
        x = self._blueprint._get('Chunk.energy.restwav', self._extension)
        if x and not ObsBlueprint.needs_lookup(x):
            self._wcs.wcs.restwav = x
        x = self._blueprint._get('Chunk.energy.velosys', self._extension)
        if x and not ObsBlueprint.needs_lookup(x):
            self._wcs.wcs.velosys = x
        x = self._blueprint._get('Chunk.energy.zsource', self._extension)
        if x and not ObsBlueprint.needs_lookup(x):
            self._wcs.wcs.zsource = x
        x = self._blueprint._get('Chunk.energy.ssyssrc', self._extension)
        if x and not ObsBlueprint.needs_lookup(x):
            self._wcs.wcs.ssyssrc = x
        x = self._blueprint._get('Chunk.energy.velang', self._extension)
        if x and not ObsBlueprint.needs_lookup(x):
            self._wcs.wcs.velangl = x
        return

    def _finish_position(self):
        x = self._blueprint._get('Chunk.position.coordsys', self._extension)
        if x and not ObsBlueprint.needs_lookup(x):
            self._wcs.wcs.radesys = x
        x = self._blueprint._get('Chunk.position.equinox', self._extension)
        if x and not ObsBlueprint.needs_lookup(x):
            self._wcs.wcs.equinox = _to_float(x)

    def _finish_time(self):
        x = self._blueprint._get('Chunk.time.exposure', self._extension)
        if x and not ObsBlueprint.needs_lookup(x):
            self._wcs.wcs.xposure = x
        x = self._blueprint._get('Chunk.time.timesys', self._extension)
        if x and not ObsBlueprint.needs_lookup(x):
            self._wcs.wcs.timesys = x
        x = self._blueprint._get('Chunk.time.trefpos', self._extension)
        if x and not ObsBlueprint.needs_lookup(x):
            self._wcs.wcs.trefpos = x
        x = self._blueprint._get('Chunk.time.mjdref', self._extension)
        if x and not ObsBlueprint.needs_lookup(x):
            self._wcs.wcs.mjdref = x


def _to_str(value):
    return str(value).strip() if value is not None else None


def _to_float(value):
    return float(value) if value is not None else None


def _to_int(value):
    return int(value) if value is not None else None


def _to_int_32(value):
    if value is None:
        return None
    elif isinstance(value, str):
        return int_32(value)
    else:
        return value


def _to_checksum_uri(value):
    if value is None:
        return None
    elif isinstance(value, ChecksumURI):
        return value
    else:
        return ChecksumURI(value)


def _set_by_type(header, keyword, value):
    """astropy documentation says that the type of the second
    parameter in the 'set' call is 'str', and then warns of expectations
    for floating-point values when the code does that, so make float values
    into floats, and int values into ints."""
    float_value = None
    int_value = None

    try:
        float_value = float(value)
    except ValueError:
        pass

    try:
        int_value = int(value)
    except ValueError:
        pass

    if (float_value and not str(value).isdecimal() or
            re.match(r'0\.0*', str(value))):
        header.set(keyword, float_value)
    elif int_value:
        header.set(keyword, int_value)
    else:
        header.set(keyword, value)


def get_external_headers(external_url):
    try:
        session = requests.Session()
        retries = 10
        retry = Retry(total=retries, read=retries, connect=retries,
                      backoff_factor=0.5)
        adapter = HTTPAdapter(max_retries=retry)
        session.mount('http://', adapter)
        session.mount('https://', adapter)
        r = session.get(external_url, timeout=20)
        if r.status_code == requests.codes.ok:
            headers = data_util.make_headers_from_string(r.text)
        else:
            headers = None
            logging.warning('Error {} when retrieving {} headers.'.format(
                r.status_code, external_url))
        r.close()
        return headers
    except Exception as e:
        logging.error(f'Connection failed to {external_url}.\n{e}')
        raise RuntimeError(e)


def get_vos_headers(uri, subject=None):
    """
    Creates the FITS headers object from a vospace file.
    :param uri: vos URI
    :param subject: user credentials. Anonymous if subject is None
    :return: List of headers corresponding to each extension. Each header is
    of astropy.wcs.Header type - essentially a dictionary of FITS keywords.
    """
    if uri.startswith('vos'):
        if subject is not None and subject.certificate is not None:
            client = Client(subject.certificate)
        else:
            client = Client()

        temp_filename = tempfile.NamedTemporaryFile()
        client.copy(uri, temp_filename.name, head=True)
        return data_util.get_local_file_headers(
            f'file://{temp_filename.name}'
        )
    else:
        # this should be a programming error by now
        raise NotImplementedError('Only vos type URIs supported')


def _get_and_update_artifact_meta(uri, artifact, subject=None, connected=True,
                                  client=None):
    """
    Updates contentType, contentLength and contentChecksum of an artifact
    :param artifact:
    :param subject: User credentials
    :param connected: True if there's a network connection
    :param client: connection to CADC storage
    :return:
    """
    logging.debug(f'Begin _get_and_update_artifact_meta for {uri}')
    file_url = urlparse(uri)
    if file_url.scheme == 'gemini' and '.jpg' not in file_url.path:
        # will get file metadata from Gemini JSON summary for fits,
        # because the metadata is available long before the data
        # will be stored at CADC
        return
    elif file_url.scheme == 'vos':
        metadata = _get_vos_meta(subject, uri)
    elif file_url.scheme == 'file':
        if (file_url.path.endswith('.header') and subject is not None and
                connected):
            if artifact.uri.startswith('vos'):
                metadata = _get_vos_meta(subject, artifact.uri)
            else:
                # if header is on disk, get the content_* from CADC
                metadata = client.info(artifact.uri)
                if metadata is None:
                    logging.info(
                        'Could not find {} at CADC. No Artifact '
                        'metadata.'.format(artifact.uri))
                    return
        else:
            metadata = data_util.get_local_file_info(file_url.path)
    else:
        metadata = client.info(uri)
        if metadata is None:
            logging.info('Could not find {} at CADC. No Artifact '
                         'metadata.'.format(artifact.uri))
            return

    update_artifact_meta(artifact, metadata)


def update_artifact_meta(artifact, file_info):
    """
    Updates contentType, contentLength and contentChecksum of an artifact
    :param artifact:
    :param file_info
    :return:
    """
    logging.debug('old artifact metadata - '
                  'uri({}), encoding({}), size({}), type({})'.
                  format(artifact.uri,
                         artifact.content_checksum,
                         artifact.content_length,
                         artifact.content_type))
    if file_info.md5sum is not None:
        if file_info.md5sum.startswith('md5:'):
            checksum = ChecksumURI(file_info.md5sum)
        else:
            checksum = ChecksumURI(f'md5:{file_info.md5sum}')
        artifact.content_checksum = checksum
    artifact.content_length = _to_int(file_info.size)
    artifact.content_type = _to_str(file_info.file_type)
    logging.debug('updated artifact metadata - '
                  'uri({}), encoding({}), size({}), type({})'.
                  format(artifact.uri,
                         artifact.content_checksum,
                         artifact.content_length,
                         artifact.content_type))


def _get_vos_meta(subject, uri):
    """
    Gets contentType, contentLength and contentChecksum of a VOS artifact
    :param subject: user credentials
    :param uri:
    :return:
    """
    if subject is not None and subject.certificate is not None:
        client = Client(subject.certificate)
    else:
        client = Client()
    node = client.get_node(uri, limit=None, force=False)
    return FileInfo(id=uri, size=node.props['length'],
                    md5sum=node.props['MD5'],
                    file_type=data_util.get_file_type(uri))


def _lookup_blueprint(blueprints, uri):
    """
    Blueprint handling may be one-per-observation, or one-per-URI. Find
    the correct one here.
    :param blueprints: The collection of blueprints provided by the user.
    :param uri: Which blueprint to look for
    :return: the blueprint to apply to Observation creation.
    """
    if len(blueprints) == 1:
        key, value = blueprints.popitem()
        return value
    else:
        return blueprints[uri]


def _extract_ids(cardinality):
    """
    Localize cardinality structure knowledge.

    :param cardinality:
    :return: product_id, artifact URI
    """
    return cardinality.split('/', 1)


def _augment(obs, product_id, uri, blueprint, subject, dumpconfig=False,
             validate_wcs=True, plugin=None, local=None,
             external_url=None, connected=True, use_blueprint_parser=False,
             client=None, **kwargs):
    """
    Find or construct a plane and an artifact to go with the observation
    under augmentation.

    :param obs: Observation - target of CAOM2 model augmentation
    :param product_id: Unique identifier for a plane in an Observation
    :param uri: Unique identifier for an artifact in a plane
    :param blueprint: Which blueprint to use when mapping from a telescope
        data model to CAOM2
    :param subject: authorization for any metdata access
    :param dumpconfig: print the blueprint to stdout
    :param validate_wcs: if true, call the validate method on the constructed
        observation, which checks that the WCS in the CAOM model is valid,
    :param plugin: what code to use for modifying a CAOM instance
    :param local: the input is the name of a file on disk
    :param external_url: if header information should be retrieved
        externally, this is where to find it
    :param client: StorageClientWrapper
    :return: an updated Observation
    """
    if dumpconfig:
        print(f'Blueprint for {uri}: {blueprint}')

    if product_id not in obs.planes.keys():
        obs.planes.add(Plane(product_id=str(product_id)))

    plane = obs.planes[product_id]

    if uri not in plane.artifacts.keys():
        plane.artifacts.add(
            Artifact(uri=str(uri),
                     product_type=ProductType.SCIENCE,
                     release_type=ReleaseType.DATA))

    meta_uri = uri
    visit_local = None
    if use_blueprint_parser:
        logging.debug(
            f'Using a BlueprintParser as requested for {uri}')
        parser = BlueprintParser(blueprint, uri=uri)
    elif local:
        if uri.startswith('vos'):
            if '.fits' in local or '.fits.gz' in local:
                meta_uri = f'file://{local}'
                logging.debug(
                    f'Using a FitsParser for vos local {local}')
                headers = data_util.get_local_file_headers(local)
                parser = FitsParser(headers, blueprint, uri=uri)
            elif '.csv' in local:
                logging.debug(
                    f'Using a BlueprintParser for vos local {local}')
                parser = BlueprintParser(blueprint, uri=uri)
            else:
                raise ValueError(f'Unexpected file type {local}')
        else:
            meta_uri = f'file://{local}'
            visit_local = local
            if ('.header' in local or
                data_util.get_file_type(local) ==
                    'application/fits'):
                logging.debug(
                    f'Using a FitsParser for local file {local}')
                parser = FitsParser(local, blueprint, uri=uri)
            elif '.h5' in local:
                logging.debug(
                    f'Using an Hdf5Parser for local file {local}')
                # h5py is an extra in this package since most collections do
                # not require it
                import h5py
                temp = h5py.File(local)
                parser = Hdf5Parser(blueprint, uri, temp)
            else:
                # explicitly ignore headers for txt and image files
                logging.debug(f'Using a BlueprintParser for {local}')
                parser = BlueprintParser(blueprint, uri=uri)
    elif external_url:
        headers = get_external_headers(external_url)
        if headers is None:
            logging.debug(
                'Using a BlueprintParser for un-retrievable remote headers '
                '{}'.format(uri)
            )
            parser = BlueprintParser(blueprint, uri=uri)
        else:
            logging.debug(
                f'Using a FitsParser for remote headers {uri}')
            parser = FitsParser(headers, blueprint, uri=uri)
    else:
        if '.fits' in uri:
            if uri.startswith('vos'):
                headers = get_vos_headers(uri, subject)
            elif uri.startswith('file'):
                headers = data_util.get_local_headers_from_fits(uri)
            else:
                headers = client.get_head(uri)
            logging.debug(f'Using a FitsParser for remote file {uri}')
            parser = FitsParser(headers, blueprint, uri=uri)
        else:
            # explicitly ignore headers for txt and image files
            logging.debug(
                f'Using a BlueprintParser for remote file {uri}')
            parser = BlueprintParser(blueprint, uri=uri)

    if parser is None:
        result = None
    else:
        _get_and_update_artifact_meta(
            meta_uri, plane.artifacts[uri], subject, connected, client)
        parser.augment_observation(observation=obs, artifact_uri=uri,
                                   product_id=plane.product_id)

        result = _visit(plugin, parser, obs, visit_local, product_id, uri,
                        subject, **kwargs)

        if result is not None:
            if validate_wcs:
                try:
                    validate(obs)
                except InvalidWCSError as e:
                    logging.error(e)
                    tb = traceback.format_exc()
                    logging.debug(tb)
                    raise e

        if len(parser._errors) > 0:
            logging.debug(
                '{} errors encountered while processing {!r}.'.format(
                    len(parser._errors), uri))
            logging.debug(f'{parser._errors}')

    return result


def _load_module(module):
    """If a user provides code for execution during blueprint configuration,
    add that code to the execution environment of the interpreter here.

    :param module the fully-qualified path name to the source code from a
        user.
    """
    mname = os.path.basename(module)
    if '.' in mname:
        # remove extension from the provided name
        mname = mname.split('.')[0]
    pname = os.path.dirname(module)
    sys.path.append(pname)
    try:
        return importlib.import_module(mname)
    except ImportError as e:
        logging.debug(f'Looking for {mname!r} in {pname!r}')
        raise e


def caom2gen():
    parser = get_gen_proc_arg_parser()
    parser.add_argument('--blueprint', nargs='+', required=True,
                        help=('list of files with blueprints for CAOM2 '
                              'construction, in serialized format. If the '
                              'list is of length 1, the same blueprint will '
                              'be applied to all lineage entries. Otherwise, '
                              'there must be a blueprint file per lineage '
                              'entry.'))

    if len(sys.argv) < 2:
        parser.print_usage(file=sys.stderr)
        sys.stderr.write(f'{APP_NAME}: error: too few arguments\n')
        sys.exit(-1)

    args = parser.parse_args()
    _set_logging(args.verbose, args.debug, args.quiet)

    module = None
    if args.module:
        module = _load_module(args.module)

    blueprints = {}
    if len(args.blueprint) == 1:
        # one blueprint to rule them all
        temp = ' '.join(ii for ii in args.lineage)
        if '.h5' in temp:
            blueprint = Hdf5ObsBlueprint(module=module)
        else:
            blueprint = ObsBlueprint(module=module)
        blueprint.load_from_file(args.blueprint[0])
        for i, cardinality in enumerate(args.lineage):
            product_id, uri = _extract_ids(cardinality)
            blueprints[uri] = blueprint
    else:
        # there needs to be the same number of blueprints as plane/artifact
        # identifiers
        if len(args.lineage) != len(args.blueprint):
            logging.debug(f'Lineage: {args.lineage}')
            logging.debug(f'Blueprints: {args.blueprint}')
            sys.stderr.write(
                '{}: error: different number of blueprints '
                '{}  and files {}.'.format(APP_NAME, len(args.blueprint),
                                           len(args.lineage)))
            sys.exit(-1)

        for i, cardinality in enumerate(args.lineage):
            product_id, uri = _extract_ids(cardinality)
            logging.debug('Loading blueprint for {} from {}'.format(
                uri, args.blueprint[i]))
            if '.h5' in uri:
                blueprint = Hdf5ObsBlueprint(module=module)
            else:
                blueprint = ObsBlueprint(module=module)
            blueprint.load_from_file(args.blueprint[i])
            blueprints[uri] = blueprint

    try:
        gen_proc(args, blueprints)
    except Exception as e:
        logging.error('Failed caom2gen execution.')
        logging.error(e)
        tb = traceback.format_exc()
        logging.error(tb)
        sys.exit(-1)

    logging.debug(f'Done {APP_NAME} processing.')


def _gen_obs(obs_blueprints, in_obs_xml, collection=None, obs_id=None):
    """
    Determine whether to create a Simple or Derived Observation, or to
    read an existing Observation from an input file.

    :param obs_blueprints: Collection of blueprints provided to application.
    :param in_obs_xml: Existing observation information, contains the
        collection and obs_id values.
    :param collection: This plus the obs_id is a unique key for an
        observation.
    :param obs_id: This plus the collection is a unique key for an
        observation.
    :return: Initially constructed Observation.
    """
    obs = None
    if in_obs_xml:
        # append to existing observation
        reader = ObservationReader(validate=True)
        obs = reader.read(in_obs_xml)
    else:
        # determine the type of observation to create by looking for the
        # the DerivedObservation.members in the blueprints. If present
        # in any of it assume derived
        for bp in obs_blueprints.values():
            if bp._get('DerivedObservation.members') is not None:
                logging.debug('Build a DerivedObservation')
                obs = DerivedObservation(
                    collection=collection,
                    observation_id=obs_id,
                    algorithm=Algorithm('composite'))
                break
            elif bp._get('CompositeObservation.members') is not None:
                logging.debug(
                    'Build a CompositeObservation with obs_id {}'.format(
                        obs_id))
                obs = CompositeObservation(
                    collection=collection, observation_id=obs_id,
                    algorithm=Algorithm('composite'))
                break
    if not obs:
        # build a simple observation
        logging.debug(
            f'Build a SimpleObservation with obs_id {obs_id}')
        obs = SimpleObservation(collection=collection,
                                observation_id=obs_id,
                                algorithm=Algorithm('exposure'))
    return obs


def _set_logging(verbose, debug, quiet):
    logger = logging.getLogger()
    # replace the StreamHandler with one that has custom formatters
    if logger.handlers:
        for handler in logger.handlers:
            if not isinstance(handler, TimedRotatingFileHandler):
                logger.removeHandler(handler)

    handler = logging.StreamHandler()
    handler.setFormatter(DispatchingFormatter({
        'caom2utils.fits2caom2.FitsWcsParser': logging.Formatter(
            '%(asctime)s:%(levelname)s:%(name)-12s:HDU:%(hdu)-2s:'
            '%(lineno)d:%(message)s'),
        'astropy': logging.Formatter(
            '%(asctime)s:%(levelname)s:%(name)-12s:HDU:%(hdu)-2s:'
            '%(lineno)d:%(message)s')
    },
        logging.Formatter('%(asctime)s:%(levelname)s:%(name)-12s:'
                          '%(lineno)d:%(message)s')
    ))
    logger.addHandler(handler)
    if verbose:
        logger.setLevel(logging.INFO)
        handler.setLevel(logging.INFO)
    elif debug:
        logger.setLevel(logging.DEBUG)
        handler.setLevel(logging.DEBUG)
    elif quiet:
        logger.setLevel(logging.ERROR)
        handler.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.WARN)
        handler.setLevel(logging.WARN)


def _get_common_arg_parser():
    """
    Returns the arg parser with common arguments between
    fits2caom2 and caom2gen
    :return: args parser
    """
    resource_id = "ivo://cadc.nrc.ca/fits2caom2"
    parser = util.get_base_parser(subparsers=False,
                                  version=version.version,
                                  default_resource_id=resource_id)

    parser.description = (
        'Augments an observation with information in one or more fits files.')

    parser.add_argument('--dumpconfig', action='store_true',
                        help=('output the utype to keyword mapping to '
                              'the console'))

    parser.add_argument('--not_connected', action='store_true',
                        help=('if set, there is no internet connection, so '
                              'skip service invocations.'))

    parser.add_argument('--no_validate', action='store_true',
                        help=('by default, the application will validate the '
                              'WCS information for an observation. '
                              'Specifying this flag skips that step.'))

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
    parser.add_argument('--local', nargs='+',
                        help=('list of files in local filesystem (same order '
                              'as uri)'))
    return parser


def get_arg_parser():
    """
    Returns the arg parser with minimum arguments required to run
    fits2caom2
    :return: args parser
    """
    parser = _get_common_arg_parser()
    parser.add_argument('--productID',
                        help='product ID of the plane in the observation',
                        required=False)
    parser.add_argument('fileURI', help='URI of a fits file', nargs='+')
    return parser


def proc(args, obs_blueprints):
    """
    Function to process an observation according to command line arguments
    and a dictionary of blueprints.

    This implementation mirrors the Java implementation of fits2caom2, and
    the command line arguments it handles are productID and fileURI or
    local.

    There is no support for plugin execution to modify the blueprint with
    this access point.

    :param args: argparse args object containing the user supplied arguments.
    Arguments correspond to the parser returned by the get_arg_parser function
    :param obs_blueprints: dictionary of blueprints reguired to process the
    observation. The fileURIs represent the keys in this dictionary. Every
    fileURI in args.fileURI should have a corresponding blueprint.
    :return:
    """

    _set_logging(args.verbose, args.debug, args.quiet)

    if args.local and (len(args.local) != len(args.fileURI)):
        msg = ('number of local arguments not the same with file '
               'URIs ({} vs {})').format(len(args.local), args.fileURI)
        raise RuntimeError(msg)

    if args.in_obs_xml:
        obs = _gen_obs(obs_blueprints, args.in_obs_xml)
    else:
        obs = _gen_obs(obs_blueprints, None, args.observation[0],
                       args.observation[1])

    if args.in_obs_xml and len(obs.planes) != 1:
        if not args.productID:
            msg = '{}{}{}'.format(
                'A productID parameter is required if ',
                'there are zero or more than one planes ',
                'in the input observation.')
            raise RuntimeError(msg)

    subject = net.Subject.from_cmd_line_args(args)
    client = data_util.StorageClientWrapper(subject, resource_id=args.resource_id)
    validate_wcs = True
    if args.no_validate:
        validate_wcs = False

    for i, uri in enumerate(args.fileURI):
        blueprint = obs_blueprints[uri]
        # override the command-line argument for the plane product ID value
        product_id = blueprint._get('Plane.productID')
        if product_id is None or isinstance(product_id, tuple):
            if args.productID:
                product_id = args.productID
            else:
                msg = '{}{}'.format(
                    'A productID parameter is required if one is not ',
                    'identified in the blueprint.')
                raise RuntimeError(msg)

        file_name = None
        if args.local:
            file_name = args.local[i]

        obs = _augment(obs, product_id, uri, blueprint, subject,
                       args.dumpconfig, validate_wcs, plugin=None,
                       local=file_name, client=client)

    _write_observation(obs, args)


def _load_plugin(plugin_name):
    plgin = _load_module(plugin_name)
    if hasattr(plgin, 'update'):
        pass
    elif hasattr(plgin, 'ObservationUpdater'):
        # for backwards compatibility with caom2repo
        plgin = getattr(plgin, 'ObservationUpdater')()

    if not hasattr(plgin, 'update'):
        msg = 'The plugin {} is not correct.  It must provide one ' \
                      'of:\n' \
                      '1 - a function named update, or\n' \
                      '2 - a class ObservationUpdater with a function named ' \
                      'update.\n In either case, the update signature needs ' \
                      'to be (Observation, **kwargs).'.format(plugin_name)
        raise ImportError(msg)
    return plgin


def _visit(plugin_name, parser, obs, visit_local, product_id=None, uri=None,
           subject=None, **kwargs):
    result = obs
    if plugin_name is not None and len(plugin_name) > 0:
        # TODO make a check that's necessary under both calling conditions here
        logging.debug(
            'Begin plugin execution {!r} update method on '
            'observation {!r}'.format(plugin_name, obs.observation_id))
        plgin = _load_plugin(plugin_name)
        if isinstance(parser, FitsParser):
            kwargs['headers'] = parser.headers
        if visit_local is not None:
            kwargs['fqn'] = visit_local
        if product_id is not None:
            kwargs['product_id'] = product_id
        if uri is not None:
            kwargs['uri'] = uri
        if subject is not None:
            kwargs['subject'] = subject
        try:
            result = plgin.update(observation=obs, **kwargs)
            if result is not None:
                logging.debug(
                    'Finished executing plugin {!r} update '
                    'method on observation {!r}'.format(
                        plugin_name, obs.observation_id))
        except Exception as e:
            logging.error(e)
            tb = traceback.format_exc()
            logging.debug(tb)
            raise e
    return result


def _write_observation(obs, args):
    writer = ObservationWriter()
    if args.out_obs_xml:
        writer.write(obs, args.out_obs_xml)
    else:
        sys.stdout.flush()
        writer.write(obs, sys.stdout)


def gen_proc(args, blueprints, **kwargs):
    """The implementation that expects a product ID to be provided as
    part of the lineage parameter, and blueprints as input parameters,
    and a plugin parameter, that supports external programmatic blueprint
    modification."""
    _set_logging(args.verbose, args.debug, args.quiet)
    result = 0

    if args.in_obs_xml:
        obs = _gen_obs(blueprints, args.in_obs_xml)
    else:
        obs = _gen_obs(blueprints, None, args.observation[0],
                       args.observation[1])

    validate_wcs = True
    if args.no_validate:
        validate_wcs = False
    connected = True
    client = None
    subject = None
    if args.not_connected:
        connected = False
    else:
        subject = net.Subject.from_cmd_line_args(args)
        if args.resource_id is None:
            # if the resource_id is Undefined, using CadcDataClient
            client = data_util.StorageClientWrapper(
                subject, using_storage_inventory=False)
        else:
            # if the resource_id is defined, assume that the caller intends to
            # use the Storage Inventory system, as it's the CADC storage
            # client that depends on a resource_id
            client = data_util.StorageClientWrapper(
                subject, resource_id=args.resource_id)

    for ii, cardinality in enumerate(args.lineage):
        product_id, uri = _extract_ids(cardinality)
        blueprint = _lookup_blueprint(blueprints, uri)
        logging.debug(
            'Begin augmentation for product_id {}, uri {}'.format(product_id,
                                                                  uri))

        file_name = None
        if args.local:
            file_name = args.local[ii]

        external_url = None
        if args.external_url:
            external_url = args.external_url[ii]

        use_blueprint_parser = False
        if args.use_blueprint_parser:
            use_blueprint_parser = uri in args.use_blueprint_parser

        obs = _augment(obs, product_id, uri, blueprint, subject,
                       args.dumpconfig, validate_wcs, args.plugin, file_name,
                       external_url, connected, use_blueprint_parser, client,
                       **kwargs)

        if obs is None:
            logging.warning('No observation. Stop processing.')
            break

    if obs is None:
        if args.in_obs_xml:
            log_id = args.lineage
        else:
            log_id = args.observation
        logging.warning(f'No Observation generated for {log_id}')
        result = -1
    else:
        _write_observation(obs, args)
    return result


def get_gen_proc_arg_parser():
    """
    Returns the arg parser with minimum arguments required to run
    caom2gen
    :return: args parser
    """
    parser = _get_common_arg_parser()
    parser.add_argument('--external_url',  nargs='+',
                        help=('service endpoint(s) that '
                              'return(s) a string that can be '
                              'made into FITS headers. Cardinality should'
                              'be consistent with lineage.'))
    parser.add_argument('--module', help=('if the blueprint contains function '
                                          'calls, call '
                                          'importlib.import_module '
                                          'for the named module. Provide a '
                                          'fully qualified name. Parameter '
                                          'choices are the artifact URI (uri) '
                                          'or a list of astropy Header '
                                          'instances (header). This will '
                                          'allow the update of a single '
                                          'blueprint entry with a single '
                                          'call.'))
    parser.add_argument('--plugin', help=('if this parameter is specified, '
                                          'call importlib.import_module '
                                          'for the named module. Then '
                                          'execute the method "update", '
                                          'with the signature '
                                          '(Observation, **kwargs). '
                                          'This will allow '
                                          'for the update of multiple '
                                          'observation data members with one '
                                          'call.'))
    parser.add_argument('--lineage', nargs='+',
                        help=('productID/artifactURI. List of plane/artifact '
                              'identifiers that will be'
                              'created for the identified observation.'))
    parser.add_argument('--use_blueprint_parser', nargs='+',
                        help=('productID/artifactURI. List of lineage entries '
                              'that will be processed with a BlueprintParser. '
                              'Good for files with no metadata in the '
                              'content.'))
    return parser


def augment(blueprints, no_validate=False, dump_config=False, plugin=None,
            out_obs_xml=None, in_obs_xml=None, collection=None,
            observation=None, product_id=None, uri=None, netrc=False,
            file_name=None, verbose=False, debug=False, quiet=False, **kwargs):
    _set_logging(verbose, debug, quiet)
    logging.debug(
        'Begin augmentation for product_id {}, uri {}'.format(product_id,
                                                              uri))

    # The 'visit_args' are a dictionary within the 'params' dictionary.
    # They are set by the collection-specific implementation, as they are
    # dependent on that collection-specific implementation. The args to the
    # visit function are not set in fits2caom2.

    params = kwargs.get('params')
    kwargs = {}
    if params is not None:
        kwargs = params.get('visit_args')

    obs = _gen_obs(blueprints, in_obs_xml, collection, observation)
    subject = net.Subject(username=None, certificate=None, netrc=netrc)
    validate_wcs = True
    if no_validate is not None:
        validate_wcs = not no_validate

    for ii in blueprints:
        obs = _augment(obs, product_id, uri, blueprints[ii], subject,
                       dump_config, validate_wcs, plugin, file_name, **kwargs)

    writer = ObservationWriter()
    if out_obs_xml:
        writer.write(obs, out_obs_xml)
    else:
        sys.stdout.flush()
        writer.write(obs, sys.stdout)
    logging.info('Done augment.')


augment.__doc__ = get_gen_proc_arg_parser().format_help()
