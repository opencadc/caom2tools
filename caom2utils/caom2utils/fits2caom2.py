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
from astropy.wcs import Wcsprm
from astropy.io import fits
from cadcutils import version
from caom2.caom_util import int_32
from caom2 import Artifact, Part, Chunk, Plane, Observation, CoordError
from caom2 import SpectralWCS, CoordAxis1D, Axis, CoordFunction1D, RefCoord
from caom2 import SpatialWCS, Dimension2D, Coord2D, CoordFunction2D
from caom2 import CoordAxis2D, PolarizationWCS, TemporalWCS
from caom2 import ObservationReader, ObservationWriter, Algorithm
from caom2 import ReleaseType, ProductType, ObservationIntentType
from caom2 import DataProductType, Telescope, Environment
from caom2 import Instrument, Proposal, Target, Provenance, Metrics
from caom2 import CalibrationLevel, Requirements, DataQuality, PlaneURI
from caom2 import SimpleObservation, CompositeObservation, ChecksumURI
from caom2 import ObservationURI, ObservableAxis, Slice
from caom2utils.caomvalidator import validate
import importlib
import logging
import os
import re
import sys
import traceback
from abc import ABCMeta
from six import add_metaclass
from hashlib import md5
from os import stat
from six.moves.urllib.parse import urlparse
from cadcutils import net, util
from cadcdata import CadcDataClient
from io import BytesIO
import six

APP_NAME = 'caom2gen'

__all__ = ['FitsParser', 'WcsParser', 'DispatchingFormatter',
           'ObsBlueprint', 'get_cadc_headers', 'get_arg_parser', 'proc',
           'POLARIZATION_CTYPES']

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


class HDULoggingFilter(logging.Filter):
    """Add the HDU number to logging messages as a default."""

    def __init__(self):
        super(HDULoggingFilter, self).__init__()
        self._extension = -1

    def filter(self, record):
        record.hdu = self._extension
        return True

    def extension(self, value):
        self._extension = value


class classproperty(object):
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


class ObsBlueprint(object):
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
    ob.set_fits_attribute('Chunk.energy.axis.axis.ctype', ['MYCTYPE'],
                          extension=1)
    ob.add_fits_attribute('Chunk.energy.axis.axis.ctype', 'MYCTYPE2',
                          extension=1)
    ob.set('Chunk.energy.velang', 33, extension=1)
    ob.set_default('Chunk.position.coordsys', 'RA-DEC', extension=1)

    ob.set('Chunk.energy.velang', 44, extension=2)
    print(ob)

    """
    _CAOM2_ELEMENTS = [
        'CompositeObservation.members',
        'Observation.observationID',
        'Observation.type',
        'Observation.intent',
        'Observation.sequenceNumber',
        'Observation.metaRelease',
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

        'Artifact.productType',
        'Artifact.releaseType',
        'Artifact.contentChecksum',
        'Artifact.contentLength',
        'Artifact.contentType',
        'Artifact.uri',

        'Part.name',
        'Part.productType',

        'Chunk.naxis',
        'Chunk.observableAxis',
        'Chunk.positionAxis1',
        'Chunk.positionAxis2',
        'Chunk.energyAxis',
        'Chunk.timeAxis',
        'Chunk.polarizationAxis',

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
        'Chunk.observable.axis.function.refCoord.pix'
        ]

    # replace _CAOM2_ELEMENTS in __doc__ with the real elements
    __doc__ = __doc__.replace('_CAOM2_ELEMENTS', '\n'.join(['\t\t{}'.format(
        elem) for elem in _CAOM2_ELEMENTS]))

    def __init__(self, position_axes=None, energy_axis=None,
                 polarization_axis=None, time_axis=None,
                 obs_axis=None, module=None):
        """
        Ctor
        :param position_axes: tuple of form (int, int) indicating the indexes
        of position axis
        :param energy_axis: index of energy axis (int)
        :param polarization_axis: index of polarization axis (int)
        :param time_axis: index of time axis (int)
        :param obs_axis: index of observable axis (int)
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
               'Observation.algorithm.name': (['PROCNAME'], None),
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
               'Plane.metaRelease': (['RELEASE', 'REL_DATE'], None),
               'Plane.dataRelease': (['RELEASE', 'REL_DATE'], None),
               'Plane.productID': (['RUNID'], None),
               'Plane.provenance.name': (['XPRVNAME'], None),
               'Plane.provenance.project': (['ADC_ARCH'], None),
               'Plane.provenance.producer': (['ORIGIN'], None),
               'Plane.provenance.reference': (['XREFER'], None),
               'Plane.provenance.lastExecuted': (['DATE-FTS'], None)
               }
        # using the tmp to make sure that the keywords are valid
        for key in tmp:
            self.set(key, tmp[key])

        self._extensions = {}

        # contains the standard WCS keywords in the FITS file expected by the
        # astropy.WCS package.
        self._wcs_std = {
            # 'Chunk.naxis': (['ZNAXIS', 'NAXIS'], None)
            'Chunk.naxis': 'ZNAXIS,NAXIS'
        }
        self._pos_axes_configed = False
        self._energy_axis_configed = False
        self._time_axis_configed = False
        self._pol_axis_configed = False
        self._obs_axis_configed = False
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

        if module:
            self._module = module
        else:
            self._module = None

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
            self.set('Chunk.position.coordsys', (['RADECSYS', 'RADESYS'],
                                                 None))
            self.set('Chunk.position.equinox', (['EQUINOX', 'EPOCH'], None))
            self.set('Chunk.position.axis.axis1.ctype',
                     (['CTYPE{}'.format(axes[0])], None))
            self.set('Chunk.position.axis.axis1.cunit',
                     (['CUNIT{}'.format(axes[0])], None))
            self.set('Chunk.position.axis.axis2.ctype',
                     (['CTYPE{}'.format(axes[1])], None))
            self.set('Chunk.position.axis.axis2.cunit',
                     (['CUNIT{}'.format(axes[1])], None))
            self.set('Chunk.position.axis.error1.syser',
                     (['CSYER{}'.format(axes[0])], None))
            self.set('Chunk.position.axis.error1.rnder',
                     (['CRDER{}'.format(axes[0])], None))
            self.set('Chunk.position.axis.error2.syser',
                     (['CSYER{}'.format(axes[1])], None))
            self.set('Chunk.position.axis.error2.rnder',
                     (['CRDER{}'.format(axes[1])], None))
            self.set('Chunk.position.axis.function.cd11',
                     (['CD{}_{}'.format(axes[0], axes[0])], None))
            self.set('Chunk.position.axis.function.cd12',
                     (['CD{}_{}'.format(axes[0], axes[1])], None))
            self.set('Chunk.position.axis.function.cd21',
                     (['CD{}_{}'.format(axes[1], axes[0])], None))
            self.set('Chunk.position.axis.function.cd22',
                     (['CD{}_{}'.format(axes[1], axes[1])], None))
            self.set('Chunk.position.axis.function.dimension.naxis1',
                     (['ZNAXIS{}'.format(axes[0]),
                       'NAXIS{}'.format(axes[0])], None))
            self.set('Chunk.position.axis.function.dimension.naxis2',
                     (['ZNAXIS{}'.format(axes[1]),
                       'NAXIS{}'.format(axes[1])], None))
            self.set('Chunk.position.axis.function.refCoord.coord1.pix',
                     (['CRPIX{}'.format(axes[0])], None))
            self.set('Chunk.position.axis.function.refCoord.coord1.val',
                     (['CRVAL{}'.format(axes[0])], None))
            self.set('Chunk.position.axis.function.refCoord.coord2.pix',
                     (['CRPIX{}'.format(axes[1])], None))
            self.set('Chunk.position.axis.function.refCoord.coord2.val',
                     (['CRVAL{}'.format(axes[1])], None))

        self._wcs_std['Chunk.position.coordsys'] = 'RADECSYS'
        self._wcs_std['Chunk.position.equinox'] = 'EQUINOX'

        self._wcs_std['Chunk.position.axis.axis1.ctype'] = \
            'CTYPE{}'.format(axes[0])
        self._wcs_std['Chunk.position.axis.axis1.cunit'] = \
            'CUNIT{}'.format(axes[0])
        self._wcs_std['Chunk.position.axis.axis2.ctype'] = \
            'CTYPE{}'.format(axes[1])
        self._wcs_std['Chunk.position.axis.axis2.cunit'] = \
            'CUNIT{}'.format(axes[1])
        self._wcs_std['Chunk.position.axis.error1.syser'] = \
            'CSYER{}'.format(axes[0])
        self._wcs_std['Chunk.position.axis.error1.rnder'] = \
            'CRDER{}'.format(axes[0])
        self._wcs_std['Chunk.position.axis.error2.syser'] = \
            'CSYER{}'.format(axes[1])
        self._wcs_std['Chunk.position.axis.error2.rnder'] = \
            'CRDER{}'.format(axes[1])
        self._wcs_std['Chunk.position.axis.function.cd11'] = \
            'CD{}_{}'.format(axes[0], axes[0])
        self._wcs_std['Chunk.position.axis.function.cd12'] = \
            'CD{}_{}'.format(axes[0], axes[1])
        self._wcs_std['Chunk.position.axis.function.cd21'] = \
            'CD{}_{}'.format(axes[1], axes[0])
        self._wcs_std['Chunk.position.axis.function.cd22'] = \
            'CD{}_{}'.format(axes[1], axes[1])
        self._wcs_std['Chunk.position.axis.function.dimension.naxis1'] = \
            'NAXIS{}'.format(axes[0])
        self._wcs_std['Chunk.position.axis.function.dimension.naxis2'] = \
            'NAXIS{}'.format(axes[1])
        self._wcs_std['Chunk.position.axis.function.refCoord.coord1.pix'] \
            = 'CRPIX{}'.format(axes[0])
        self._wcs_std['Chunk.position.axis.function.refCoord.coord1.val'] \
            = 'CRVAL{}'.format(axes[0])
        self._wcs_std['Chunk.position.axis.function.refCoord.coord2.pix'] \
            = 'CRPIX{}'.format(axes[1])
        self._wcs_std['Chunk.position.axis.function.refCoord.coord2.val'] \
            = 'CRVAL{}'.format(axes[1])

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
                     (['CTYPE{}'.format(axis)], None))
            self.set('Chunk.energy.axis.axis.cunit',
                     (['CUNIT{}'.format(axis)], None))
            self.set('Chunk.energy.axis.error.syser',
                     (['CSYER{}'.format(axis)], None))
            self.set('Chunk.energy.axis.error.rnder',
                     (['CRDER{}'.format(axis)], None))
            self.set('Chunk.energy.axis.function.naxis',
                     (['NAXIS{}'.format(axis)], None))
            self.set('Chunk.energy.axis.function.delta',
                     (['CDELT{}'.format(axis)], None))
            self.set('Chunk.energy.axis.function.refCoord.pix',
                     (['CRPIX{}'.format(axis)], None))
            self.set('Chunk.energy.axis.function.refCoord.val',
                     (['CRVAL{}'.format(axis)], None))

        self._wcs_std['Chunk.energy.specsys'] = 'SPECSYS'
        self._wcs_std['Chunk.energy.ssysobs'] = 'SSYSOBS'
        self._wcs_std['Chunk.energy.restfrq'] = 'RESTFRQ'
        self._wcs_std['Chunk.energy.restwav'] = 'RESTWAV'
        self._wcs_std['Chunk.energy.velosys'] = 'VELOSYS'
        self._wcs_std['Chunk.energy.zsource'] = 'ZSOURCE'
        self._wcs_std['Chunk.energy.ssyssrc'] = 'SSYSSRC'
        self._wcs_std['Chunk.energy.velang'] = 'VELANG'

        self._wcs_std['Chunk.energy.axis.axis.ctype'] = \
            'CTYPE{}'.format(axis)
        self._wcs_std['Chunk.energy.axis.axis.cunit'] = \
            'CUNIT{}'.format(axis)
        self._wcs_std['Chunk.energy.axis.error.syser'] = \
            'CSYER{}'.format(axis)
        self._wcs_std['Chunk.energy.axis.error.rnder'] = \
            'CRDER{}'.format(axis)
        self._wcs_std['Chunk.energy.axis.function.naxis'] = \
            'NAXIS{}'.format(axis)
        self._wcs_std['Chunk.energy.axis.function.delta'] = \
            'CDELT{}'.format(axis)
        self._wcs_std['Chunk.energy.axis.function.refCoord.pix'] = \
            'CRPIX{}'.format(axis)
        self._wcs_std['Chunk.energy.axis.function.refCoord.val'] = \
            'CRVAL{}'.format(axis)

        self._energy_axis_configed = True

    def configure_polarization_axis(self, axis, override=True):
        """
        Set the expected FITS polarization keywords by index in the blueprint
        and the wcs_std lookup.

        :param axis: The index expected for the polarization axis.
        :param override: Set to False when reading from a file.
        :return:
        """
        if self._pol_axis_configed:
            self.logger.debug(
                'Attempt to configure already-configured polarization axis.')
            return

        if override:
            self.set('Chunk.polarization.axis.axis.ctype',
                     (['CTYPE{}'.format(axis)], None))
            self.set('Chunk.polarization.axis.axis.cunit',
                     (['CUNIT{}'.format(axis)], None))
            self.set('Chunk.polarization.axis.function.naxis',
                     (['NAXIS{}'.format(axis)], None))
            self.set('Chunk.polarization.axis.function.delta',
                     (['CDELT{}'.format(axis)], None))
            self.set('Chunk.polarization.axis.function.refCoord.pix',
                     (['CRPIX{}'.format(axis)], None))
            self.set('Chunk.polarization.axis.function.refCoord.val',
                     (['CRVAL{}'.format(axis)], None))

        self._wcs_std['Chunk.polarization.axis.axis.ctype'] = \
            'CTYPE{}'.format(axis)
        self._wcs_std['Chunk.polarization.axis.axis.cunit'] = \
            'CUNIT{}'.format(axis)
        self._wcs_std['Chunk.polarization.axis.function.naxis'] = \
            'NAXIS{}'.format(axis)
        self._wcs_std['Chunk.polarization.axis.function.delta'] = \
            'CDELT{}'.format(axis)
        self._wcs_std['Chunk.polarization.axis.function.refCoord.pix'] = \
            'CRPIX{}'.format(axis)
        self._wcs_std['Chunk.polarization.axis.function.refCoord.val'] = \
            'CRVAL{}'.format(axis)

        self._pol_axis_configed = True

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
                     (['CTYPE{}'.format(axis)], None))
            self.set('Chunk.observable.axis.axis.cunit',
                     (['CUNIT{}'.format(axis)], None))
            self.set('Chunk.observable.axis.function.refCoord.pix',
                     (['CRPIX{}'.format(axis)], None))

        self._wcs_std['Chunk.observable.axis.axis.ctype'] = \
            'CTYPE{}'.format(axis)
        self._wcs_std['Chunk.observable.axis.axis.cunit'] = \
            'CUNIT{}'.format(axis)
        self._wcs_std['Chunk.observable.axis.function.refCoord.pix'] = \
            'CRPIX{}'.format(axis)

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
                     (['CTYPE{}'.format(axis)], None))
            self.set('Chunk.time.axis.axis.cunit',
                     (['CUNIT{}'.format(axis)], None))
            self.set('Chunk.time.axis.error.syser',
                     (['CSYER{}'.format(axis)], None))
            self.set('Chunk.time.axis.error.rnder',
                     (['CRDER{}'.format(axis)], None))
            self.set('Chunk.time.axis.function.naxis',
                     (['NAXIS{}'.format(axis)], None))
            self.set('Chunk.time.axis.function.delta',
                     (['CDELT{}'.format(axis)], None))
            self.set('Chunk.time.axis.function.refCoord.pix',
                     (['CRPIX{}'.format(axis)], None))
            self.set('Chunk.time.axis.function.refCoord.val',
                     (['CRVAL{}'.format(axis)], None))

        self._wcs_std['Chunk.time.exposure'] = 'EXPTIME'
        self._wcs_std['Chunk.time.resolution'] = 'TIMEDEL'
        self._wcs_std['Chunk.time.timesys'] = 'TIMESYS'
        self._wcs_std['Chunk.time.trefpos'] = 'TREFPOS'
        self._wcs_std['Chunk.time.mjdref'] = 'MJDREF'

        self._wcs_std['Chunk.time.axis.axis.ctype'] = \
            'CTYPE{}'.format(axis)
        self._wcs_std['Chunk.time.axis.axis.cunit'] = \
            'CUNIT{}'.format(axis)
        self._wcs_std['Chunk.time.axis.error.syser'] = \
            'CSYER{}'.format(axis)
        self._wcs_std['Chunk.time.axis.error.rnder'] = \
            'CRDER{}'.format(axis)
        self._wcs_std['Chunk.time.axis.function.naxis'] = \
            'NAXIS{}'.format(axis)
        self._wcs_std['Chunk.time.axis.function.delta'] = \
            'CDELT{}'.format(axis)
        self._wcs_std['Chunk.time.axis.function.refCoord.pix'] = \
            'CRPIX{}'.format(axis)
        self._wcs_std['Chunk.time.axis.function.refCoord.val'] = \
            'CRVAL{}'.format(axis)

        self._time_axis_configed = True

    def _guess_axis_info_from_plan(self):
        """Look for info regarding axis types in the blueprint wcs_std.
        Configure the blueprint according to the guesses.
        """
        # a data structure to carry around twelve bits of data at a time:
        # the first item in the set is the ctype index, and the second is
        # whether or not the index means anything, resulting in a
        # call to the blueprint configure_* methods if it's True.
        axis_info = {
            'dec': (0, False),
            'energy': (0, False),
            'obs': (0, False),
            'polarization': (0, False),
            'ra': (0, False),
            'time': (0, False)}

        for ii in self._plan:
            if isinstance(self._plan[ii], tuple):
                for value in self._plan[ii][0]:
                    if (value.startswith('CTYPE')) and value[-1].isdigit():
                        value = value.split('-')[0]
                        self._guess_axis_info_from_ctypes(ii, int(value[-1]),
                                                          axis_info)
            else:
                value = self._plan[ii]
                if (value.startswith('CTYPE')) and value[-1].isdigit():
                    value = value.split('-')[0]
                    self._guess_axis_info_from_ctypes(ii, int(value[-1]),
                                                      axis_info)

        configured_index = 0
        for ii in self._plan:
            if ii.startswith('Chunk.position') and ii.endswith('axis1.ctype') \
                    and not axis_info['ra'][1]:
                configured_index = self._get_configured_index(axis_info, 'ra')
                axis_info['ra'] = (configured_index, True)
            elif ii.startswith('Chunk.position') and \
                    ii.endswith('axis2.ctype') and not axis_info['dec'][1]:
                configured_index = self._get_configured_index(axis_info,
                                                              'dec')
                axis_info['dec'] = (configured_index, True)
            elif ii.startswith('Chunk.energy') and not axis_info['energy'][1]:
                configured_index = self._get_configured_index(axis_info,
                                                              'energy')
                axis_info['energy'] = (configured_index, True)
            elif ii.startswith('Chunk.time') and not axis_info['time'][1]:
                configured_index = self._get_configured_index(axis_info,
                                                              'time')
                axis_info['time'] = (configured_index, True)
            elif ii.startswith('Chunk.polarization') \
                    and not axis_info['polarization'][1]:
                configured_index = self._get_configured_index(axis_info,
                                                              'polarization')
                axis_info['polarization'] = (configured_index, True)
            elif ii.startswith('Chunk.observable') and not axis_info['obs'][1]:
                configured_index = self._get_configured_index(axis_info,
                                                              'obs')
                axis_info['obs'] = (configured_index, True)

        if axis_info['ra'][1] and axis_info['dec'][1]:
            self.configure_position_axes(
                (axis_info['ra'][0], axis_info['dec'][0]), False)
        elif axis_info['ra'][1] or axis_info['dec'][1]:
            raise ValueError('Only one positional axis found '
                             '(ra/dec): {}/{}'.
                             format(axis_info['ra'][0], axis_info['dec'][0]))
        else:
            # assume that positional axis are 1 and 2 by default
            if axis_info['time'][0] in [1, 2] or \
                    axis_info['energy'][0] in [1, 2] or \
                    axis_info['polarization'][0] in [1, 2] or \
                    axis_info['obs'][0] in [1, 2]:
                raise ValueError('Cannot determine the positional axis')
            else:
                self.configure_position_axes((1, 2), False)

        if axis_info['time'][1]:
            self.configure_time_axis(axis_info['time'][0], False)
        if axis_info['energy'][1]:
            self.configure_energy_axis(axis_info['energy'][0], False)
        if axis_info['polarization'][1]:
            self.configure_polarization_axis(axis_info['polarization'][0],
                                             False)
        if axis_info['obs'][1]:
            self.configure_observable_axis(axis_info['obs'][0], False)

    def _guess_axis_info_from_ctypes(self, lookup, counter, axis_info):
        """
        Check for the presence of blueprint keys in the plan, and whether or
        not they indicate an index in their configuration.

        :param lookup: Blueprint plan key.
        :param counter: Value to set the index to for an axis.
        :param axis_info: local data structure to pass around what is
            configured, and what is it's value.
        """
        if lookup.startswith('Chunk.energy'):
            axis_info['energy'] = (counter, True)
        elif lookup.startswith('Chunk.polarization'):
            axis_info['polarization'] = (counter, True)
        elif lookup.startswith('Chunk.time'):
            axis_info['time'] = (counter, True)
        elif lookup.startswith('Chunk.position') and lookup.endswith(
                'axis1.ctype'):
            axis_info['ra'] = (counter, True)
        elif lookup.startswith('Chunk.position') and lookup.endswith(
                'axis2.ctype'):
            axis_info['dec'] = (counter, True)
        elif lookup.startswith('Chunk.observable'):
            axis_info['obs'] = (counter, True)
        else:
            raise ValueError(
                'Unrecognized axis type: {}'.format(lookup))

    def _get_configured_index(self, axis_info, lookup):
        """Find the next available index value among those that are not set.

        :param axis_info: local data structure to pass around what is
            configured, and what is it's value."""
        DEFAULT_INDICES = {'ra': 1,
                           'dec': 2,
                           'energy': 3,
                           'time': 4,
                           'polarization': 5,
                           'obs': 6}

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
                    self.set(key.strip(), cleaned_up_value)
        self._guess_axis_info_from_plan()

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

    def __str__(self):
        plan = self._serialize(self._plan)

        extensions = ''
        if self._extensions:
            for key in sorted(self._extensions):
                extensions = extensions + '\nextension {}:\n'.format(key) +\
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
        assert extension >= 0, 'Extension count failure when setting a value.'
        if extension:
            ObsBlueprint.check_chunk(caom2_element)
            if extension not in self._extensions:
                self._extensions[extension] = {}
            self._extensions[extension][caom2_element] = value
        else:
            self._plan[caom2_element] = value

    def set_fits_attribute(self, caom2_element, fits_attribute_list,
                           extension=0):
        """
        Associates a CAOM2 element with a FITS attribute
        :param caom2_element: name CAOM2 element (as in
        ObsBlueprint.CAOM2_ELEMEMTS)
        :param fits_attribute_list: list of FITS attributes the element is
        potentially mapped to.
        :param extension: extension number (used only for Chunk elements)
        """
        ObsBlueprint.check_caom2_element(caom2_element)
        assert isinstance(fits_attribute_list,
                          list), 'Type mis-match for fits_attribute_list'
        assert extension >= 0, 'Extension count failure when setting a value'
        if extension:
            ObsBlueprint.check_chunk(caom2_element)
            if extension not in self._extensions:
                self._extensions[extension] = {}
            self._extensions[extension][caom2_element] = (fits_attribute_list,
                                                          None)
        else:
            self._plan[caom2_element] = (fits_attribute_list, None)

    def add_fits_attribute(self, caom2_element, fits_attribute, extension=0):
        """
        Adds a FITS attribute in the list of other FITS attributes associated
        with an caom2 element.
        :param caom2_element: name CAOM2 element (as in
        ObsBlueprint.CAOM2_ELEMEMTS)
        :param fits_attribute: name of FITS attribute the element is mapped to
        :param extension: extension number (used only for Chunk elements)
        :raises AttributeError if the caom2 element has already an associated
        value or KeyError if the caom2 element does not exists.
        """
        ObsBlueprint.check_caom2_element(caom2_element)
        assert extension >= 0, 'Extension failure when adding FITS attribute.'
        if extension:
            ObsBlueprint.check_chunk(caom2_element)
            if extension not in self._extensions:
                raise AttributeError(
                    'No extension {} in the blueprint'.format(extension))
            else:
                if caom2_element in self._extensions[extension]:
                    if isinstance(self._extensions[extension][caom2_element],
                                  tuple):
                        self._extensions[extension][caom2_element][0].\
                            insert(0, fits_attribute)
                    else:
                        raise AttributeError(
                            ('No FITS attributes in extension {} associated '
                             'with keyword {}').format(extension,
                                                       caom2_element))
                else:
                    raise KeyError(
                        ('Keyword {} not found in the extension {} of '
                         'the blueprint').format(caom2_element, extension))
        else:
            if caom2_element in self._plan:
                if isinstance(self._plan[caom2_element], tuple):
                    self._plan[caom2_element][0].insert(0, fits_attribute)
                else:
                    raise AttributeError(
                        'No FITS attributes associated with keyword {}'.
                        format(caom2_element))
            else:
                raise KeyError(
                    'Keyword {} not found in the blueprint'.
                    format(caom2_element))

    def set_default(self, caom2_element, default, extension=0):
        """
        Sets the default value of a caom2 element that is associated with FITS
        attributes. If the element does not exist or does not have a list of
        associated FITS attributes, default is set as the associated value
        of the element.

        If set_fits_attribute is called for the same caom2_element after this,
        the default value will be reset to None.

        :param caom2_element: name CAOM2 element (as in
        ObsBlueprint.CAOM2_ELEMEMTS)
        :param default: default value
        :param extension: extension number (used only for Chunk elements)
        """
        ObsBlueprint.check_caom2_element(caom2_element)
        assert extension >= 0, 'Extension failure when setting default.'
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
        assert extension >= 0, 'Extension failure when deleting.'
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

    def _get(self, caom2_element, extension=0):
        """
        Returns the source associated with a CAOM2 element
        :param caom2_element: name CAOM2 element (as in
        ObsBlueprint.CAOM2_ELEMEMTS)
        :param extension: extension number (used only for Chunk elements)
        :return: Tuple of the form (list_of_associated_fits_attributes,
        default_value) OR the actual value associated with the CAOM2 element
        """
        ObsBlueprint.check_caom2_element(caom2_element)
        if extension:
            ObsBlueprint.check_chunk(caom2_element)
            if (extension in self._extensions) and \
                    (caom2_element in self._extensions[extension]):
                return self._extensions[extension][caom2_element]

        # look in the generic plan
        if caom2_element not in self._plan:
            return None
        else:
            return self._plan[caom2_element]

    def get_configed_axes_count(self):
        configed_axes = 0
        if self._pos_axes_configed:
            configed_axes += 2
        if self._energy_axis_configed:
            configed_axes += 1
        if self._time_axis_configed:
            configed_axes += 1
        if self._pol_axis_configed:
            configed_axes += 1
        if self._obs_axis_configed:
            configed_axes += 1
        return configed_axes

    @staticmethod
    def is_fits(value):
        """Hide the blueprint structure from clients - they shouldn't need
        to know that a value of type tuple requires special processing."""
        return isinstance(value, tuple)

    @staticmethod
    def is_function(value):
        return (not ObsBlueprint.is_fits(value)
                and isinstance(value, str) and '(' in value and ')' in value)


@add_metaclass(ABCMeta)
class GenericParser:
    """
    Extract CAOM2 metadata from files with no WCS information.
    """
    def __init__(self, obs_blueprint=None, logging_name=None, uri=None):
        if obs_blueprint:
            self._blueprint = obs_blueprint
        else:
            self._blueprint = ObsBlueprint()
        self._errors = []
        self.logging_name = logging_name
        self.logger = logging.getLogger(__name__)
        self.uri = uri

    @property
    def blueprint(self):
        return self._blueprint

    @blueprint.setter
    def blueprint(self, value):
        self._blueprint = value

    def augment_observation(self, observation, artifact_uri, product_id=None):
        """
        Augments a given observation with plane structure only.
        :param observation: existing CAOM2 observation to be augmented.
        :param artifact_uri: the key for finding the artifact to augment
        :param product_id: the key for finding for the plane to augment
        """
        self.logger.debug(
            'Begin generic CAOM2 observation augmentation for URI {}.'.format(
                artifact_uri))
        assert observation, 'observation instance required.'
        assert isinstance(observation,
                          Observation), 'Observation type mis-match.'

        plane = None
        if not product_id:
            product_id = self._get_from_list('Plane.productID', index=0)
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
            'End generic CAOM2 observation augmentation for {}.'.format(
                artifact_uri))

    def augment_plane(self, plane, artifact_uri):
        """
        Augments a given plane with artifact structure only.
        :param plane: existing CAOM2 plane to be augmented.
        :param artifact_uri:
        """
        self.logger.debug(
            'Begin generic CAOM2 plane augmentation for {}.'.format(
                artifact_uri))
        assert plane, 'plane instance required.'
        assert isinstance(plane, Plane), 'Plane type mis-match.'
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
        self.augment_artifact(artifact)
        self.logger.debug(
            'End generic CAOM2 plane augmentation for {}.'.format(
                artifact_uri))

    def augment_artifact(self, artifact):
        """
        Augments a given CAOM2 artifact with available FITS information
        :param artifact: existing CAOM2 artifact to be augmented
        """
        self.logger.debug(
            'Begin generic CAOM2 artifact augmentation for {}.'.format(
                self.logging_name))
        assert artifact, 'artifact instance required.'
        assert isinstance(artifact, Artifact), 'Artifact type mis-match'

        artifact.uri = self._get_from_list('Artifact.uri', index=0,
                                           current=artifact.uri)
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
        self.logger.debug(
            'End generic CAOM2 artifact augmentation for {}.'.format(
                self.logging_name))

    def _get_set_from_list(self, lookup, index):
        value = None
        keywords = None
        try:
            keywords = self.blueprint._get(lookup)
        except KeyError:
            self.add_error(lookup, sys.exc_info()[1])
            self.logger.debug(
                'Could not find \'{}\' in fits2caom2 configuration.'.format(
                    lookup))
        if keywords:
            value = keywords
            self.logger.debug('{}: assigned value {}.'.format(lookup, value))

        return value

    def _get_from_list(self, lookup, index, current=None):
        value = None
        try:
            keywords = self.blueprint._get(lookup)
        except KeyError:
            self.add_error(lookup, sys.exc_info()[1])
            self.logger.debug(
                'Could not find {!r} in fits2caom2 configuration.'.format(
                    lookup))
            if current:
                self.logger.debug(
                    '{}: using current value of {!r}.'.format(lookup, current))
                value = current
            return value
        if keywords:
            value = keywords
        elif current:
            value = current

        self.logger.debug('{}: value is {}'.format(lookup, value))
        return value

    def add_error(self, key, message):
        self._errors.append(('{} {} {}'.format(
            datetime.now().strftime('%Y-%m-%dT%H:%M:%S'), key, message)))

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
                'Must set a value of {} for {}.'.format(to_enum_type.__name__,
                                                        self.logging_name))
        elif isinstance(value, to_enum_type):
            return value
        else:
            return to_enum_type(value)


class FitsParser(GenericParser):
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
            self._headers = _get_headers_from_fits(self.file)
        super(FitsParser, self).__init__(obs_blueprint, self.file, uri)
        # for command-line parameter to module execution
        self.uri = uri
        self.apply_blueprint_to_fits()

    @property
    def headers(self):
        """
        List of headers where each header should allow dictionary like
        access to the FITS attribute in that header
        :return:
        """
        return self._headers

    @property
    def blueprint(self):
        return self._blueprint

    @blueprint.setter
    def blueprint(self, value):
        self._blueprint = value
        self.apply_blueprint_to_fits()

    def augment_artifact(self, artifact):
        """
        Augments a given CAOM2 artifact with available FITS information
        :param artifact: existing CAOM2 artifact to be augmented
        """
        super(FitsParser, self).augment_artifact(artifact)

        self.logger.debug(
            'Begin artifact augmentation for {} with {} HDUs.'.format(
                artifact.uri, len(self.headers)))

        if self.blueprint.get_configed_axes_count() == 0:
            self.logger.debug(
                'No WCS Data. End artifact augmentation for {}.'.format(
                    artifact.uri))
            return

        for i, header in enumerate(self.headers):
            ii = str(i)

            # there is one Part per extension, the name is the extension number
            # Assumption:
            #    Only primary headers for 1 extension files or the extensions
            # for multiple extension files can have data and therefore
            # corresponding parts
            if (i > 0) or (len(self.headers) == 1):
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

            wcs_parser = WcsParser(header, self.file, ii)
            # NOTE: astropy.wcs does not distinguished between WCS axes and
            # data array axes. naxis in astropy.wcs represents in fact the
            # number of WCS axes, whereas chunk.axis represents the naxis
            # of the data array. Solution is to determine it directly from
            # the header
            if 'ZNAXIS' in header:
                chunk.naxis = _to_int(header['ZNAXIS'])
            elif 'NAXIS' in header:
                chunk.naxis = _to_int(header['NAXIS'])
            else:
                chunk.naxis = self._get_from_list('Chunk.naxis', 0,
                                                  wcs_parser.wcs.wcs.naxis)
            if self.blueprint._pos_axes_configed:
                wcs_parser.augment_position(chunk)
                if chunk.position is None:
                    self._try_position_with_blueprint(chunk, i)
            if self.blueprint._energy_axis_configed:
                wcs_parser.augment_energy(chunk)
            if chunk.energy:
                chunk.energy.bandpass_name = self._get_from_list(
                    'Chunk.energy.bandpassName', index=i)
                chunk.energy.transition = self._get_from_list(
                    'Chunk.energy.transition', index=i)
                chunk.energy.resolving_power = _to_float(self._get_from_list(
                    'Chunk.energy.resolvingPower', index=i))
            else:
                if self.blueprint._energy_axis_configed:
                    self._try_energy_with_blueprint(chunk, i)
            if self.blueprint._time_axis_configed:
                wcs_parser.augment_temporal(chunk)
                if chunk.time is None:
                    self._try_time_with_blueprint(chunk, i)
            if self.blueprint._pol_axis_configed:
                wcs_parser.augment_polarization(chunk)
                if chunk.polarization is None:
                    self._try_polarization_with_blueprint(chunk, i)
            if self.blueprint._obs_axis_configed:
                wcs_parser.augment_observable(chunk)

        self.logger.debug(
            'End artifact augmentation for {}.'.format(artifact.uri))

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
            'Chunk.position.axis.axis1.csyer',
            'Chunk.position.axis.axis1.crder', index, _to_str, CoordError)
        aug_y_error = self._two_param_constructor(
            'Chunk.position.axis.axis2.csyer',
            'Chunk.position.axis.axis2.crder', index, _to_str, CoordError)
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
                'Creating position Coord2D for {}'.format(self.uri))

        aug_function = None
        if (aug_dimension is not None and aug_ref_coord is not None and
                aug_cd11 is not None and aug_cd12 is not None and
                aug_cd21 is not None and aug_cd22 is not None):
            aug_function = CoordFunction2D(aug_dimension, aug_ref_coord,
                                           aug_cd11, aug_cd12, aug_cd21,
                                           aug_cd22)
            self.logger.debug(
                'Creating position CoordFunction2D for {}'.format(self.uri))

        aug_axis = None
        if (aug_x_axis is not None and aug_y_axis is not None and
                aug_function is not None):
            aug_axis = CoordAxis2D(aug_x_axis, aug_y_axis, aug_x_error,
                                   aug_y_error, None, None, aug_function)
            self.logger.debug(
                'Creating position CoordAxis2D for {}'.format(self.uri))

        if aug_axis is not None:
            if chunk.position:
                chunk.position.axis = aug_axis
            else:
                chunk.position = SpatialWCS(aug_axis)

        if chunk.position:
            chunk.position.coordsys = self._get_from_list(
                'Chunk.position.coordsys', index)
            chunk.position.equinox = self._get_from_list(
                'Chunk.position.equinox', index)
            chunk.position.resolution = self._get_from_list(
                'Chunk.position.resolution', index)
        self.logger.debug('End augmentation with blueprint for position.')

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

        chunk.time_axis = self._get_from_list('Chunk.energyAxis', index)

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
                    'Creating PolarizationWCS for {} from blueprint'.
                    format(self.uri))

        self.logger.debug('End augmentation with blueprint for polarization.')

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
        if not chunk.energy:
            chunk.energy = SpectralWCS(aug_naxis, specsys)
        else:
            chunk.energy.naxis = aug_naxis
            chunk.energy.specsys = specsys

        if chunk.energy is not None:
            chunk.energy.ssysobs = self._get_from_list('Chunk.energy.ssysobs',
                                                       index)
            chunk.energy.restfrq = self._get_from_list('Chunk.energy.restfrq',
                                                       index)
            chunk.energy.restwav = self._get_from_list('Chunk.energy.restwav',
                                                       index)
            chunk.energy.velosys = self._get_from_list('Chunk.energy.velosys',
                                                       index)
            chunk.energy.zsource = self._get_from_list('Chunk.energy.zsource',
                                                       index)
            chunk.energy.ssyssrc = self._get_from_list('Chunk.energy.ssyssrc',
                                                       index)
            chunk.energy.velang = self._get_from_list('Chunk.energy.velang',
                                                      index)
            chunk.energy.bandpass_name = self._get_from_list(
                'Chunk.energy.bandpassName', index)
            chunk.energy.transition = self._get_from_list(
                'Chunk.energy.transition', index)
            chunk.energy.resolving_power = _to_float(self._get_from_list(
                'Chunk.energy.resolvingPower', index))
        self.logger.debug('End augmentation with blueprint for energy.')

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

    def _get_naxis(self, label, index):
        """Helper function to construct a CoordAxis1D instance, with all
        it's members, from the blueprint.

        :param label: axis name - must be one of 'energy', 'time', or
        'polarization', as it's used for the blueprint lookup.
        :param index: which blueprint index to find a value in
        :return an instance of CoordAxis1D
        """
        self.logger.debug(
            'Begin {} naxis construction from blueprint.'.format(label))

        aug_axis_ctype = self._get_from_list(
            'Chunk.{}.axis.axis.ctype'.format(label), index)
        aug_axis_cunit = self._get_from_list(
            'Chunk.{}.axis.axis.cunit'.format(label), index)
        aug_axis = None
        if aug_axis_ctype is not None:
            aug_axis = Axis(aug_axis_ctype, aug_axis_cunit)
            self.logger.debug(
                'Creating polarization Axis for {} from blueprint'.
                format(self.uri))

        aug_error = self._two_param_constructor(
            'Chunk.{}.axis.error.syser'.format(label),
            'Chunk.{}.axis.error.rnder'.format(label),
            index, _to_float, CoordError)
        aug_error = self._two_param_constructor(
            'Chunk.{}.axis.error.syser'.format(label),
            'Chunk.{}.axis.error.rnder'.format(label),
            index, _to_float, CoordError)
        aug_ref_coord = self._two_param_constructor(
            'Chunk.{}.axis.function.refCoord.pix'.format(label),
            'Chunk.{}.axis.function.refCoord.val'.format(label),
            index, _to_float, RefCoord)
        aug_delta = _to_float(
            self._get_from_list('Chunk.{}.axis.function.delta'.format(label),
                                index))
        aug_length = _to_int(
            self._get_from_list('Chunk.{}.axis.function.naxis'.format(label),
                                index))

        aug_function = None
        if (aug_length is not None and aug_delta is not None and
                aug_ref_coord is not None):
            aug_function = \
                CoordFunction1D(aug_length, aug_delta, aug_ref_coord)
            self.logger.debug(
                'Creating {} function for {} from blueprint'.
                format(label, self.uri))

        aug_naxis = None
        if aug_axis is not None and aug_function is not None:
            aug_naxis = CoordAxis1D(aug_axis, aug_error, None, None,
                                    aug_function)
            self.logger.debug(
                'Creating {} CoordAxis1D for {} from blueprint'.
                format(label, self.uri))
        self.logger.debug(
            'End {} naxis construction from blueprint.'.format(label))
        return aug_naxis

    def augment_observation(self, observation, artifact_uri, product_id=None):
        """
        Augments a given observation with available FITS information.
        :param observation: existing CAOM2 observation to be augmented.
        :param artifact_uri: the key for finding the artifact to augment
        :param product_id: the key for finding for the plane to augment
        """
        super(FitsParser, self).augment_observation(observation, artifact_uri,
                                                    product_id)
        self.logger.debug(
            'Begin observation augmentation for URI {}.'.format(
                artifact_uri))
        members = self._get_members(observation)
        if members:
            for m in members.split():
                observation.members.add(ObservationURI(m))
        observation.algorithm = self._get_algorithm(observation)

        observation.sequence_number = _to_int(self._get_from_list(
            'Observation.sequenceNumber', index=0))
        observation.intent = self._get_from_list('Observation.intent', 0,
                                                 ObservationIntentType.SCIENCE)
        observation.type = self._get_from_list('Observation.type', 0)
        observation.meta_release = self._get_datetime(
            self._get_from_list('Observation.metaRelease', 0))
        observation.requirements = self._get_requirements()
        observation.instrument = self._get_instrument()
        observation.proposal = self._get_proposal()
        observation.target = self._get_target()
        observation.target_position = self._get_target_position()
        observation.telescope = self._get_telescope()
        observation.environment = self._get_environment()
        self.logger.debug(
            'End observation augmentation for {}.'.format(artifact_uri))

    def augment_plane(self, plane, artifact_uri):
        """
        Augments a given plane with available FITS information.
        :param plane: existing CAOM2 plane to be augmented.
        :param artifact_uri:
        """
        super(FitsParser, self).augment_plane(plane, artifact_uri)
        self.logger.debug(
            'Begin plane augmentation for {}.'.format(artifact_uri))

        plane.meta_release = self._get_datetime(self._get_from_list(
            'Plane.metaRelease', index=0))
        plane.data_release = self._get_datetime(self._get_from_list(
            'Plane.dataRelease', index=0))
        plane.data_product_type = self._to_data_product_type(
            self._get_from_list('Plane.dataProductType', index=0,
                                current=plane.data_product_type))
        plane.calibration_level = self._to_calibration_level(_to_int_32(
            self._get_from_list('Plane.calibrationLevel', index=0)))
        plane.provenance = self._get_provenance()
        plane.metrics = self._get_metrics()
        plane.quality = self._get_quality()

        self.logger.debug(
            'End plane augmentation for {}.'.format(artifact_uri))

    def apply_blueprint_to_fits(self):

        # pointers that are short to type
        exts = self.blueprint._extensions
        wcs_std = self.blueprint._wcs_std
        plan = self.blueprint._plan

        # firstly, apply the functions
        if self.blueprint._module is not None:
            for key, value in plan.items():
                if ObsBlueprint.is_function(value):
                    plan[key] = self._execute_external(value, key)
            for extension in exts:
                for key, value in exts[extension].items():
                    if ObsBlueprint.is_function(value):
                        exts[extension][key] = self._execute_external(value,
                                                                      key)

        # apply overrides from blueprint to all extensions
        for key, value in plan.items():
            if key in wcs_std:
                if ObsBlueprint.is_fits(value):
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
                else:
                    # value provided for standard wcs attribute
                    if ObsBlueprint.is_fits(wcs_std[key]):
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
            hdr = self.headers[extension]
            for key, value in exts[extension].items():
                keywords = wcs_std[key].split(',')
                for keyword in keywords:
                    _set_by_type(hdr, keyword, value)
                    logging.debug(
                        '{}: set to {} in extension {}'.format(keyword, value,
                                                               extension))
        # apply defaults to all extensions
        for key, value in plan.items():
            if ObsBlueprint.is_fits(value) and value[1]:
                # there is a default value set
                for index, header in enumerate(self.headers):
                    for keyword in value[0]:
                        if not header.get(keyword):
                            # apply a default if a value does not already exist
                            _set_by_type(header, keyword, value[1])
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
                    if 'CDELT{}'.format(i) in header and \
                            'CD{0}_{0}'.format(i) not in header:
                        header['CD{0}_{0}'.format(i)] = \
                            header['CDELT{}'.format(i)]

        # TODO When a projection is specified, wcslib expects corresponding
        # DP arguments with NAXES attributes. Normally, omitting the attribute
        # signals no distortion which is the assumption in fits2caom2 for
        # energy and polarization axes. Following is a workaround for
        # SIP projections.
        # For more details see:
        # http://www.atnf.csiro.au/people/mcalabre/WCS/dcs_20040422.pdf
        for header in self.headers:
            sip = False
            for i in range(1, 6):
                if ('CTYPE{}'.format(i) in header) and \
                        ('-SIP' in header['CTYPE{}'.format(i)]):
                    sip = True
                    break
            if sip:
                for i in range(1, 6):
                    if ('CTYPE{}'.format(i) in header) and \
                            ('-SIP' not in header['CTYPE{}'.format(i)]) and \
                            ('DP{}'.format(i) not in header):
                        header['DP{}'.format(i)] = 'NAXES: 1'

        return

    def _execute_external(self, value, key):
        """Execute a function supplied by a user, assign a value to a
        blueprint entry. The input parameters passed to the function are the
        headers as read in by astropy, or the artifact uri.

        :param value the name of the function to apply.
        """
        # determine which of the two possible values for parameter the user
        # is hoping for
        parameter = ''
        if 'uri' in value:
            parameter = self.uri
        elif 'header' in value:
            parameter = self._headers

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
            logging.error(tb)
            logging.error(e)
        try:
            result = execute(parameter)
            logging.debug(
                'Key {} calculated value of {} using {}'.format(
                    key, result, value))
        except Exception as e:
            msg = 'Failed to execute {} for {}'.format(execute.__name__, key)
            logging.error(msg)
            logging.debug('Input parameter was {}'.format(parameter))
            self._errors.append(msg)
            tb = traceback.format_exc()
            logging.error(tb)
            logging.error(e)
        return result

    def _get_members(self, obs):
        """
        Returns the members of a composite observation (if specified)
        :param obs: observation to augment
        :return: members value
        """
        members = None
        self.logger.debug('Begin Members augmentation.')
        if isinstance(obs, SimpleObservation) and \
           self.blueprint._get('CompositeObservation.members'):
            raise TypeError(
                ('Cannot apply blueprint for CompositeObservation to a '
                 'simple observation'))
        elif isinstance(obs, CompositeObservation):
            members = self._get_from_list('CompositeObservation.members',
                                          index=0)
        self.logger.debug('End Members augmentation.')
        return members

    def _get_algorithm(self, obs):
        """
        Create an Algorithm instance populated with available FITS information.
        :return: Algorithm
        """
        self.logger.debug('Begin Algorithm augmentation.')
        # TODO DEFAULT VALUE
        name = self._get_from_list('Observation.algorithm.name', index=0,
                                   current=obs.algorithm.name)
        self.logger.debug('End Algorithm augmentation.')
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
        self.logger.debug('Begin Instrument augmentation.')
        name = self._get_from_list('Observation.instrument.name', index=0)
        keywords = self._get_from_list('Observation.instrument.keywords',
                                       index=0)
        self.logger.debug('End Instrument augmentation.')
        if name:
            instr = Instrument(str(name))
            if keywords:
                for k in keywords.split():
                    instr.keywords.add(k)
            return instr
        else:
            return None

    def _get_proposal(self):
        """
        Create a Proposal instance populated with available FITS information.
        :return: Proposal
        """
        self.logger.debug('Begin Proposal augmentation.')
        prop_id = self._get_from_list('Observation.proposal.id', index=0)
        pi = self._get_from_list('Observation.proposal.pi', index=0)
        project = self._get_from_list('Observation.proposal.project', index=0)
        title = self._get_from_list('Observation.proposal.title', index=0)
        self.logger.debug('End Proposal augmentation.')
        if prop_id:
            return Proposal(str(prop_id), pi, project, title)
        else:
            return None

    def _get_target(self):
        """
        Create a Target instance populated with available FITS information.
        :return: Target
        """
        self.logger.debug('Begin Target augmentation.')
        name = self._get_from_list('Observation.target.name', index=0)
        target_type = self._get_from_list('Observation.target.type',
                                          index=0)
        standard = self._cast_as_bool(self._get_from_list(
            'Observation.target.standard', index=0))
        redshift = self._get_from_list('Observation.target.redshift', index=0)
        keywords = self._get_set_from_list('Observation.target.keywords',
                                           index=0)  # TODO
        moving = self._cast_as_bool(
            self._get_from_list('Observation.target.moving', index=0))
        self.logger.debug('End Target augmentation.')
        if name:
            return Target(str(name), target_type, standard, redshift,
                          keywords, moving)
        else:
            return None

    def _get_target_position(self):
        """
        Create a Target Position instance populated with available FITS
        information.
        :return: Target Position
        """
        self.logger.debug('Begin TargetPosition augmentation.')
        # TODO don't know what to do here, since config file says this is
        # Chunk-level metadata
        self.logger.debug('End TargetPosition augmentation.')
        return None

    def _get_telescope(self):
        """
        Create a Telescope instance populated with available FITS information.
        :return: Telescope
        """
        self.logger.debug('Begin Telescope augmentation.')
        name = self._get_from_list('Observation.telescope.name', index=0)
        geo_x = _to_float(
            self._get_from_list('Observation.telescope.geoLocationX', index=0))
        geo_y = _to_float(
            self._get_from_list('Observation.telescope.geoLocationY', index=0))
        geo_z = _to_float(
            self._get_from_list('Observation.telescope.geoLocationZ', index=0))
        keywords = self._get_set_from_list('Observation.telescope.keywords',
                                           index=0)  # TODO
        self.logger.debug('End Telescope augmentation.')
        if name:
            self.logger.debug('name is {}'.format(name))
            return Telescope(str(name), geo_x, geo_y, geo_z, keywords)
        else:
            return None

    def _get_environment(self):
        """
        Create an Environment instance populated with available FITS
        information.
        :return: Environment
        """
        self.logger.debug('Begin Environment augmentation.')
        seeing = self._get_from_list('Observation.environment.seeing', index=0)
        humidity = _to_float(
            self._get_from_list('Observation.environment.humidity', index=0))
        elevation = self._get_from_list('Observation.environment.elevation',
                                        index=0)
        tau = self._get_from_list('Observation.environment.tau', index=0)
        wavelength_tau = self._get_from_list(
            'Observation.environment.wavelengthTau', index=0)
        ambient = _to_float(
            self._get_from_list('Observation.environment.ambientTemp',
                                index=0))
        photometric = self._cast_as_bool(self._get_from_list(
            'Observation.environment.photometric', index=0))

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
        else:
            return None

    def _get_requirements(self):
        """
        Create a Requirements instance populated with available FITS
        information.
        :return: Requirements
        """
        self.logger.debug('Begin Requirement augmentation.')
        flag = self._get_from_list('Observation.requirements.flag', index=0)
        self.logger.debug('End Requirement augmentation.')
        if flag:
            return Requirements(flag)
        else:
            return None

    def _get_from_list(self, lookup, index, current=None):
        value = None
        try:
            keywords = self.blueprint._get(lookup)
        except KeyError:
            self.add_error(lookup, sys.exc_info()[1])
            self.logger.debug(
                'Could not find {!r} in fits2caom2 configuration.'.format(
                    lookup))
            if current:
                self.logger.debug(
                    '{}: using current value of {!r}.'.format(lookup, current))
                value = current
            return value

        if isinstance(keywords, tuple):
            for ii in keywords[0]:
                try:
                    value = self.headers[index].get(ii)
                    if value:
                        self.logger.debug(
                            '{}: assigned value {} based on keyword {}.'.
                            format(lookup, value, ii))
                        break
                except (KeyError, IndexError):
                    if keywords[0].index(ii) == len(keywords[0]) - 1:
                        self.add_error(lookup, sys.exc_info()[1])
                    # assign a default value, if one exists
                    if keywords[1]:
                        value = keywords[1]
                        self.logger.debug(
                            '{}: assigned default value {}.'.format(lookup,
                                                                    value))
            if value is None:
                if current:
                    value = current
                    self.logger.debug(
                        '{}: used current value {!r}.'.format(lookup, value))
                else:
                    # assign a default value, if one exists
                    if keywords[1]:
                        value = keywords[1]
                        self.logger.debug(
                            '{}: assigned default value {}.'.format(lookup,
                                                                    value))

        elif (keywords is not None) and (keywords != ''):
            value = keywords
        elif current:
            value = current

        self.logger.debug('{}: value is {}'.format(lookup, value))
        return value

    def _get_set_from_list(self, lookup, index):
        value = None
        keywords = None
        try:
            keywords = self.blueprint._get(lookup)
        except KeyError:
            self.add_error(lookup, sys.exc_info()[1])
            self.logger.debug(
                'Could not find \'{}\' in fits2caom2 configuration.'.format(
                    lookup))

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
            self.logger.debug('{}: assigned value {}.'.format(lookup, value))

        return value

    def _get_provenance(self):
        """
        Create a Provenance instance populated with available FITS information.
        :return: Provenance
        """
        self.logger.debug('Begin Provenance augmentation.')
        name = _to_str(
            self._get_from_list('Plane.provenance.name', index=0))
        p_version = _to_str(self._get_from_list('Plane.provenance.version',
                                                index=0))
        project = _to_str(
            self._get_from_list('Plane.provenance.project', index=0))
        producer = _to_str(
            self._get_from_list('Plane.provenance.producer', index=0))
        run_id = _to_str(
            self._get_from_list('Plane.provenance.runID', index=0))
        reference = _to_str(
            self._get_from_list('Plane.provenance.reference', index=0))
        last_executed = self._get_datetime(
            self._get_from_list('Plane.provenance.lastExecuted', index=0))
        keywords = self._get_from_list('Plane.provenance.keywords', index=0)
        inputs = self._get_from_list('Plane.provenance.inputs', index=0)
        if name:
            prov = Provenance(name, p_version, project, producer, run_id,
                              reference, last_executed)
            if keywords:
                for k in keywords.split():
                    prov.keywords.add(k)
            if inputs:
                for i in inputs.split():
                    prov.inputs.add(PlaneURI(str(i)))
            self.logger.debug('End Provenance augmentation.')
            return prov
        else:
            self.logger.debug(
                'End Provenance augmentation - no provenance information.')
            return None

    def _get_metrics(self):
        """
        Create a Metrics instance populated with available FITS information.
        :return: Metrics
        """
        self.logger.debug('Begin Metrics augmentation.')
        source_number_density = self._get_from_list(
            'Plane.metrics.sourceNumberDensity', index=0)
        background = self._get_from_list('Plane.metrics.background', index=0)
        background_stddev = self._get_from_list(
            'Plane.metrics.backgroundStddev', index=0)
        flux_density_limit = self._get_from_list(
            'Plane.metrics.fluxDensityLimit', index=0)
        mag_limit = self._get_from_list('Plane.metrics.magLimit', index=0)

        if source_number_density or background or background_stddev or \
                flux_density_limit or mag_limit:
            metrics = Metrics()
            metrics.source_number_density = source_number_density
            metrics.background = background
            metrics.background_std_dev = background_stddev
            metrics.flux_density_limit = flux_density_limit
            metrics.mag_limit = mag_limit
        else:
            metrics = None
        self.logger.debug('End Metrics augmentation.')
        return metrics

    def _get_quality(self):
        """
        Create a Quality instance populated with available FITS information.
        :return: Quality
        """
        self.logger.debug('Begin Quality augmentation.')
        flag = self._get_from_list('Plane.dataQuality', index=0)
        self.logger.debug('End Quality augmentation.')
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

        if from_value:
            try:
                return datetime.strptime(from_value, '%Y-%m-%dT%H:%M:%S')
            except ValueError:
                try:
                    return datetime.strptime(from_value,
                                             '%Y-%m-%d %H:%M:%S.%f')
                except ValueError:
                    try:
                        return datetime.strptime(from_value, '%Y-%m-%d')
                    except ValueError:
                        self.logger.error(
                            'Cannot parse datetime {}'.format(from_value))
                        self.add_error('get_datetime', sys.exc_info()[1])
                        return None
        else:
            return None

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

        # add the HDU extension to logging messages from this class
        self.logger = logging.getLogger(__name__ + '.WcsParser')
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

    def augment_energy(self, chunk):
        """
        Augments the energy information in a chunk
        :param chunk:
        """
        assert chunk, 'chunk instance required.'
        assert isinstance(chunk, Chunk), 'Chunk type mis-match.'

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

        specsys = str(self.wcs.specsys)
        if not chunk.energy:
            chunk.energy = SpectralWCS(naxis, specsys)
        else:
            chunk.energy.naxis = naxis
            chunk.energy.specsys = specsys

        chunk.energy.ssysobs = _to_str(self._sanitize(self.wcs.ssysobs))
        # TODO not sure why, but wcs returns 0.0 when the FITS keywords for the
        # following two keywords are actually not present in the header
        if self._sanitize(self.wcs.restfrq) != 0:
            chunk.energy.restfrq = self._sanitize(self.wcs.restfrq)
        if self._sanitize(self.wcs.restwav) != 0:
            chunk.energy.restwav = self._sanitize(self.wcs.restwav)
        chunk.energy.velosys = self._sanitize(self.wcs.velosys)
        chunk.energy.zsource = self._sanitize(self.wcs.zsource)
        chunk.energy.ssyssrc = _to_str(self._sanitize(self.wcs.ssyssrc))
        chunk.energy.velang = self._sanitize(self.wcs.velangl)

    def augment_position(self, chunk):
        """
        Augments a chunk with spatial WCS information
        :param chunk:
        :return:
        """
        self.logger.debug('Begin Spatial WCS augmentation.')

        assert chunk, 'chunk instance required.'
        assert isinstance(chunk, Chunk), 'Chunk type mis-match.'

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
        chunk.position.equinox = self._sanitize(self.wcs.equinox)
        self.logger.debug('End Spatial WCS augmentation.')

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
        assert chunk, 'chunk instance required.'
        assert isinstance(chunk, Chunk), 'Chunk type mis-match.'

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
        aug_function = CoordFunction1D(
            self._get_axis_length(time_axis_index + 1), delta, aug_ref_coord)
        naxis = CoordAxis1D(aug_naxis, aug_error, None, None, aug_function)
        if not chunk.time:
            chunk.time = TemporalWCS(naxis)
        else:
            chunk.time.naxis = naxis

        chunk.time.exposure = _to_float(self.header.get('EXPTIME'))
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
        self.logger.debug('Begin Polarization WCS augmentation.')
        assert chunk, 'chunk instance required.'
        assert isinstance(chunk, Chunk), 'Chunk type mis-match.'

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
            chunk.polarization.naxis = naxis
        self.logger.debug('End Polarization WCS augmentation.')

    def augment_observable(self, chunk):
        """
        Augments a chunk with an observable axis
        :param chunk:
        :return:
        """
        self.logger.debug('Begin Observable WCS augmentation.')
        assert chunk, 'chunk instance required.'
        assert isinstance(chunk, Chunk), 'Chunk type mis-match.'

        observable_axis_index = self._get_axis_index(OBSERVABLE_CTYPES)
        if observable_axis_index is None:
            self.logger.debug('No Observable axis info')
            return

        chunk.observable_axis = observable_axis_index + 1
        ctype = self.header.get('CTYPE{}'.format(chunk.observable_axis))
        cunit = self.header.get('CUNIT{}'.format(chunk.observable_axis))
        pix_bin = self.header.get('CRPIX{}'.format(chunk.observable_axis))
        if ctype is not None and cunit is not None and pix_bin is not None:
            chunk.observable = ObservableAxis(
                Slice(self._get_axis(0, ctype, cunit), pix_bin))
        self.logger.debug('End Observable WCS augmentation.')

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
        else:
            aug_function = None

        aug_axis = CoordAxis2D(self._get_axis(xindex),
                               self._get_axis(yindex),
                               self._get_coord_error(xindex),
                               self._get_coord_error(yindex),
                               None, None, aug_function)
        return aug_axis

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
                'Error searching for CD* values {}'.format(sys.exc_info()[1]))
            cd11 = -1.0  # TODO what if neither of these are defined?
            cd12 = -1.0
            cd21 = -1.0
            cd22 = -1.0

        return cd11, cd12, cd21, cd22

    def _get_coord_error(self, index):
        aug_coord_error = None
        aug_csyer = self._sanitize(self.wcs.csyer[index])
        aug_crder = self._sanitize(self.wcs.crder[index])
        if aug_csyer and aug_crder:
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
        aug_ref_coord = RefCoord(aug_crpix, aug_crval)
        return aug_ref_coord

    def _get_axis_length(self, for_axis):
        result = -1
        try:
            # try ZNAXIS first in order to get the size of the original
            # image in case it was FITS compressed
            result = int(self._sanitize(
                self.header.get('ZNAXIS{}'.format(for_axis))))
            return result
        except TypeError:
            try:
                result = _to_int(self._sanitize(
                    self.header.get('NAXIS{}'.format(for_axis))))
            except ValueError:
                self.logger.debug(
                    'Could not find axis length for axis {}'.format(
                        for_axis))
        if isinstance(result, tuple):
            result = result[0]
        return result

    def _sanitize(self, value):
        """
        Sanitizes values from FITS to caom2
        :param value:
        :return:
        """
        if isinstance(value, float) and math.isnan(value):
            return None
        elif not str(value):
            return None  # empty string
        else:
            return value


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
    """astropy documentations says that the type of the second
    parameter in the 'set' call is 'str', and then warns of expectations
    for floating-point values."""
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

    if float_value and not value.isdecimal() or re.match('0\.0*', value):
        header.set(keyword, float_value)
    elif int_value:
        header.set(keyword, int_value)
    else:
        header.set(keyword, value)


def get_cadc_headers(uri, subject=None):
    """
    Creates the FITS headers object from a either a local file or it
    fetches the FITS headers of a CADC file. The function takes advantage
    of the fhead feature of the CADC storage service and retrieves just the
    headers and no data, minimizing the transfer time.
    :param uri: CADC ('ad:') or local file ('file:') URI
    :param subject: user credentials. Anonymous if subject is None
    :return: List of headers corresponding to each extension. Each header is
    of astropy.wcs.Header type - essentially a dictionary of FITS keywords.
    """
    file_url = urlparse(uri)
    if file_url.scheme == 'ad':
        # create possible types of subjects
        if not subject:
            subject = net.Subject()
        client = CadcDataClient(subject)
        # do a fhead on the file
        archive, file_id = file_url.path.split('/')
        b = BytesIO()
        b.name = uri
        client.get_file(archive, file_id, b, fhead=True)
        fits_header = b.getvalue().decode('ascii')
        b.close()
        headers = _make_headers_from_string(fits_header)
    elif file_url.scheme == 'file':
        try:
            fits_header = open(file_url.path).read()
            headers = _make_headers_from_string(fits_header)
        except UnicodeDecodeError:
            headers = _get_headers_from_fits(file_url.path)
    else:
        # TODO add hook to support other service providers
        raise NotImplementedError('Only ad type URIs supported')
    return headers


def _make_headers_from_string(fits_header):
    """Create a list of fits.Header instances from a string.
    ":param fits_header a string of keyword/value pairs"""
    delim = '\nEND'
    extensions = \
        [e + delim for e in fits_header.split(delim) if e.strip()]
    headers = [fits.Header.fromstring(e, sep='\n') for e in extensions]
    return headers


def _get_headers_from_fits(path):
    """Create a list of fits.Header instances from a fits file.
    :param path where the FITS files resides on disk."""
    hdulist = fits.open(path, memmap=True, lazy_load_hdus=False)
    hdulist.close()
    headers = [h.header for h in hdulist]
    return headers


def _update_artifact_meta(artifact, subject=None):
    """
    Updates contentType, contentLength and contentChecksum of an artifact
    :param artifact:
    :param subject: User credentials
    :return:
    """
    file_url = urlparse(artifact.uri)
    if file_url.scheme == 'ad':
        metadata = _get_cadc_meta(subject, file_url.path)
    elif file_url.scheme == 'file':
        metadata = _get_file_meta(file_url.path)
    else:
        # TODO add hook to support other service providers
        raise NotImplementedError('Only ad type URIs supported')

    checksum = ChecksumURI('md5:{}'.format(metadata['md5sum']))
    logging.debug('old artifact metadata - '
                  'uri({}), encoding({}), size({}), type({})'.
                  format(artifact.uri,
                         artifact.content_checksum,
                         artifact.content_length,
                         artifact.content_type))
    artifact.content_checksum = checksum
    artifact.content_length = int(metadata['size'])
    artifact.content_type = str(metadata['type'])
    logging.debug('updated artifact metadata - '
                  'uri({}), encoding({}), size({}), type({})'.
                  format(artifact.uri,
                         artifact.content_checksum,
                         artifact.content_length,
                         artifact.content_type))


def _get_cadc_meta(subject, path):
    """
    Gets contentType, contentLength and contentChecksum of a CADC artifact
    :param subject: user credentials
    :param path:
    :return:
    """
    client = CadcDataClient(subject)
    archive, file_id = path.split('/')
    return client.get_file_info(archive, file_id)


def _get_file_meta(path):
    """
    Gets contentType, contentLength and contentChecksum of an artifact on disk.
    :param path:
    :return:
    """
    meta = {}
    s = stat(path)
    meta['size'] = s.st_size
    meta['md5sum'] = md5(open(path, 'rb').read()).hexdigest()
    if path.endswith('.header') or path.endswith('.txt'):
        meta['type'] = 'text/plain'
    else:
        meta['type'] = 'application/octet-stream'
    return meta


def _lookup_blueprint(blueprints, uri):
    """
    Blueprint handling may be one-per-observation, or one-per-URI. Find
    the correct one here.
    :param blueprints: The collection of blueprints provided by the user.
    :param uri: Which blueprint to look for
    :return: the blueprint to apply to Observation creation.
    """
    if len(blueprints) == 1:
        return six.next(six.itervalues(blueprints))
    else:
        return blueprints[uri]


def _lookup_blueprint_name(index, blueprint_names):
    """
    Blueprint handling may be one-per-observation, or one-per-URI. Reference
    the correct name of the file here.
    :param blueprint_names: The collection of blueprint names provided by the
        user.
    :param index: Which blueprint to look for
    :return: the blueprint to apply to Observation creation.
    """
    if len(blueprint_names) == 1:
        return 'The One.'
    else:
        return blueprint_names[index]


def _extract_ids(cardinality):
    """
    Localize cardinality structure knowledge.

    :param cardinality:
    :return: product_id, artifact URI
    """
    return cardinality.split('/', 1)


def _augment(obs, product_id, uri, args, blueprint, index):
    """
    Find or construct a plane and an artifact to go with the observation
    under augmentation.

    :param obs: Observation - target of CAOM2 model augmentation
    :param product_id: Unique identifier for a plane in an Observation
    :param args: Command line arguments.
    :param uri: Unique identifier for an artifact in a plane
    :param blueprint: Which blueprint to use when mapping from a telescope
        data model to CAOM2
    :param index: How to find the file in the input parameter local that is
        the metadata augmentation source.
    :return:
    """
    if product_id not in obs.planes.keys():
        obs.planes.add(Plane(product_id=str(product_id)))

    plane = obs.planes[product_id]

    subject = net.Subject.from_cmd_line_args(args)
    if uri not in plane.artifacts.keys():
        plane.artifacts.add(
            Artifact(uri=str(uri),
                     product_type=ProductType.SCIENCE,
                     release_type=ReleaseType.DATA))
    if args.local:
        file = args.local[index]
        if file.endswith('.fits'):
            logging.debug('Using a FitsParser for {}'.format(file))
            parser = FitsParser(file, blueprint, uri=uri)
        elif '.header' in file:
            logging.debug('Using a FitsParser for {}'.format(file))
            parser = FitsParser(get_cadc_headers('file://{}'.format(file)),
                                blueprint, uri=uri)
        else:
            # explicitly ignore headers for txt and image files
            logging.debug('Using a GenericParser for {}'.format(file))
            parser = GenericParser(blueprint, uri=uri)
    else:
        if uri.endswith('.fits'):
            logging.debug('Using a FitsParser for {}'.format(uri))
            headers = get_cadc_headers(uri, subject)
            parser = FitsParser(headers, blueprint, uri=uri)
        else:
            # explicitly ignore headers for txt and image files
            logging.debug('Using a GenericParser for {}'.format(uri))
            parser = GenericParser(blueprint, uri=uri)

    _update_artifact_meta(plane.artifacts[uri], subject)

    if args.dumpconfig:
        print('Blueprint for {}: {}'.format(uri, blueprint))

    parser.augment_observation(observation=obs, artifact_uri=uri,
                               product_id=plane.product_id)

    if len(parser._errors) > 0:
        logging.debug(
            '{} errors encountered while processing {!r}.'.format(
                len(parser._errors), uri))
        logging.debug('{}'.format(parser._errors))

    if not args.no_validate:
        validate(obs)


def _load_module(module):
    """If a user provides code for execution during blueprint configuration,
    add that code to the execution environment of the interpreter here.

    :param module the fully-qualified path name to the source code from a
        user.
    """
    mname = os.path.basename(module)
    if '.' in mname:
        mname = mname.split('.')[0]
    pname = os.path.dirname(module)
    sys.path.append(pname)
    try:
        return importlib.import_module(mname)
    except ImportError as e:
        logging.error('Looking for {} in {}'.format(mname, pname))
        raise e


def caom2gen():
    parser = _get_common_arg_parser()

    parser.add_argument('--module', help=('if the blueprint contains function '
                                          'calls, call '
                                          'importlib.import_module '
                                          'for the named module. Provide a '
                                          'fully qualified name. Parameter '
                                          'choices are the artifact URI (uri) '
                                          'or a list of astropy Header '
                                          'instances (header).'))
    parser.add_argument('--blueprint', nargs='+', required=True,
                        help=('list of files with blueprints for CAOM2 '
                              'construction, in serialized format. If the '
                              'list is of length 1, the same blueprint will '
                              'be applied to all lineage entries. Otherwise, '
                              'there must be a blueprint file per lineage '
                              'entry.'))
    parser.add_argument('--lineage', nargs='+',
                        help=('productID/artifactURI. List of plane/artifact '
                              'identifiers that will be'
                              'created for the identified observation.'))

    if len(sys.argv) < 2:
        parser.print_usage(file=sys.stderr)
        sys.stderr.write('{}: error: too few arguments\n'.format(APP_NAME))
        sys.exit(-1)

    args = parser.parse_args()
    _set_arg_parser_logging(args)

    module = None
    if args.module:
        module = _load_module(args.module)

    blueprints = {}
    if len(args.blueprint) == 1:
        # one blueprint to rule them all
        blueprint = ObsBlueprint(module=module)
        blueprint.load_from_file(args.blueprint[0])
        for i, cardinality in enumerate(args.lineage):
            product_id, uri = _extract_ids(cardinality)
            blueprints[uri] = blueprint
    else:
        # there needs to be the same number of blueprints as plane/artifact
        # identifiers
        if len(args.lineage) != len(args.blueprint):
            logging.debug('Lineage: {}'.format(args.lineage))
            logging.debug('Blueprints: {}'.format(args.blueprint))
            sys.stderr.write(
                '{}: error: different number of blueprints '
                '{}  and files {}.'.format(APP_NAME, len(args.blueprint),
                                           len(args.lineage)))
            sys.exit(-1)

        for i, cardinality in enumerate(args.lineage):
            product_id, uri = _extract_ids(cardinality)
            logging.debug('Loading blueprint for {} from {}'.format(
                uri, args.blueprint[i]))
            blueprint = ObsBlueprint(module=module)
            blueprint.load_from_file(args.blueprint[i])
            blueprints[uri] = blueprint

    try:
        obs = _set_obs(args, blueprints)

        for ii, cardinality in enumerate(args.lineage):
            product_id, uri = _extract_ids(cardinality)
            bp_name = _lookup_blueprint_name(ii, args.blueprint)
            blueprint = _lookup_blueprint(blueprints, uri)
            logging.debug('Begin augmentation for product_id {}, uri {}, '
                          'with blueprint {}'.format(product_id, uri, bp_name))
            _augment(obs, product_id, uri, args, blueprint, ii)

        writer = ObservationWriter()
        if args.out_obs_xml:
            writer.write(obs, args.out_obs_xml)
        else:
            sys.stdout.flush()
            writer.write(obs, sys.stdout)

    except Exception as e:
        logging.error('Failed caom2gen execution.')
        logging.error(e)
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)

    logging.debug(
        'Done {} processing for {}'.format(APP_NAME, args.observation[1]))


def _set_obs(args, obs_blueprints):
    """
    Determine whether to create a Simple or Composite Observation.

    :param args: Command-line parameters.
    :param obs_blueprints: Collection of blueprints provided to application.
    :return: Initially constructed Observation.
    """
    obs = None
    if args.in_obs_xml:
        # append to existing observation
        reader = ObservationReader(validate=True)
        obs = reader.read(args.in_obs_xml)
    else:
        # determine the type of observation to create by looking for the
        # the CompositeObservation.members in the blueprints. If present
        # in any of it assume composite
        for bp in obs_blueprints.values():
            if bp._get('CompositeObservation.members'):
                obs = CompositeObservation(
                    collection=str(args.observation[0]),
                    observation_id=str(args.observation[1]),
                    algorithm=Algorithm(str('composite')))
                break
    if not obs:
        # build a simple observation
        obs = SimpleObservation(collection=str(args.observation[0]),
                                observation_id=str(args.observation[1]),
                                algorithm=Algorithm(str('exposure')))
    return obs


def _set_arg_parser_logging(args):
    logger = logging.getLogger()
    if args.verbose:
        logger.setLevel(logging.INFO)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    elif args.quiet:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.WARN)

    if logger.handlers:
        handler = logger.handlers[0]
        logger.removeHandler(handler)
    handler = logging.StreamHandler()
    handler.setFormatter(DispatchingFormatter({
        'caom2utils.fits2caom2.WcsParser': logging.Formatter(
            '%(levelname)s:%(name)-12s:HDU:%(hdu)-2s:%(lineno)d:%(message)s'),
        'astropy': logging.Formatter(
            '%(levelname)s:%(name)-12s:HDU:%(hdu)-2s:%(lineno)d:%(message)s')
    },
        logging.Formatter('%(levelname)s:%(name)-12s:%(lineno)d:%(message)s')
    ))
    logger.addHandler(handler)


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

    parser.add_argument('--no_validate', action='store_true',
                        help=('by default, the application will validate the '
                              'WCS information for an observation. '
                              'Specifying this flag skips that step.'))

    parser.add_argument('--ignorePartialWCS', action='store_true',
                        help='do not stop and exit upon finding partial WCS')

    parser.add_argument('-o', '--out', dest='out_obs_xml',
                        type=argparse.FileType('wb'),
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
    parser.add_argument('--keep', action='store_true',
                        help='keep the locally stored files after ingestion')
    parser.add_argument('--test', action='store_true',
                        help='test mode, do not persist to database')
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
    and a dictionary of blueprints
    :param args: argparse args object containing the user supplied arguments.
    Arguments correspond to the parser returned by the get_arg_parser function
    :param obs_blueprints: dictionary of blueprints reguired to process the
    observation. The fileURIs represent the keys in this dictionary. Every
    fileURI in args.fileURI should have a corresponding blueprint.
    :return:
    """

    _set_arg_parser_logging(args)

    if args.local and (len(args.local) != len(args.fileURI)):
        msg = ('number of local arguments not the same with file '
               'URIs ({} vs {})').format(len(args.local), args.fileURI)
        raise RuntimeError(msg)

    obs = _set_obs(args, obs_blueprints)

    if args.in_obs_xml and len(obs.planes) != 1:
        if not args.productID:
            msg = '{}{}{}'.format(
                'A productID parameter is required if ',
                'there are zero or more than one planes ',
                'in the input observation.')
            raise RuntimeError(msg)

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

        _augment(obs, product_id, uri, args, blueprint, i)

    writer = ObservationWriter()
    if args.out_obs_xml:
        writer.write(obs, args.out_obs_xml)
    else:
        sys.stdout.flush()
        writer.write(obs, sys.stdout)
