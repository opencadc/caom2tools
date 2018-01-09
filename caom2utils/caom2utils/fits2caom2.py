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
from caom2.caom_util import int_32
from caom2 import Artifact, Part, Chunk, Plane, Observation, CoordError
from caom2 import SpectralWCS, CoordAxis1D, Axis, CoordFunction1D, RefCoord
from caom2 import SpatialWCS, Dimension2D, Coord2D, CoordFunction2D
from caom2 import CoordAxis2D, PolarizationWCS, TemporalWCS
from caom2 import ObservationReader, ObservationWriter, Algorithm
from caom2 import ReleaseType, ProductType, ObservationIntentType
from caom2 import DataProductType, Telescope, Environment
from caom2 import Instrument, Proposal, Target, Provenance, Metrics
from caom2 import CalibrationLevel
from caom2 import SimpleObservation
import logging
import sys
from six.moves.urllib.parse import urlparse
from cadcutils import net
from cadcdata import CadcDataClient
from io import BytesIO

APP_NAME = 'fits2caom2'

__all__ = ['FitsParser', 'WcsParser', 'DispatchingFormatter',
           'ObservationBlueprint', 'get_cadc_headers', 'main_app',
           'update_fits_headers', 'load_config']

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


class ObservationBlueprint(object):
    """
    Class that captures the blueprint of building a CAOM2 Observation
    """
    _CAOM2_ELEMENTS = [
        #TODO not used for now 'CompositeObservation.members',
        'Observation.type',
        'Observation.intent',
        'Observation.sequenceNumber',
        'Observation.metaRelease',

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

        'Plane.metaRelease',
        'Plane.dataRelease',
        'Plane.dataProductType',
        'Plane.calibrationLevel',

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
        'Chunk.time.axis.range.end.val'
        ]

    def __init__(self, position_axis=None, energy_axis=None,
                 polarization_axis=None, time_axis=None,
                 user_supplied_config=None):
        """
        Ctor
        """

        if position_axis and isinstance(position_axis, tuple) and\
                (len(position_axis) != 2):
            raise ValueError(
                'Invalid position axis: {}. Must be tuple with 2 elements'.
                    format(str(position_axis)))

        # this is the default blueprint
        self._plan = {'Observation.metaRelease':
                      (['DATE', 'DATE-OBS', 'UTCOBS', 'UTCDATE',
                        'UTC-DATE', 'MJDOBS', 'MJD_OBS'], None),
                      'Observation.instrument.name': (['INSTRUME'], None),
                      'Observation.type': (['OBSTYPE'], None),
                      'Observation.environment.ambientTemp': (['TEMPERAT'], None),
                      'Observation.algorithm.name': (['PROCNAME'], None),
                      'Observation.instrument.keywords': (['INSTMODE'], None),
                      'Observation.proposal.id': (['RUNID'], None),
                      'Observation.target.name': (['OBJECT'], None),
                      'Observation.telescope.name': (['INSTRUME'], None),
                      'Observation.telescope.geoLocationX': (['OBSGEO-X'], None),
                      'Observation.telescope.geoLocationY': (['OBSGEO-Y'], None),
                      'Observation.telescope.geoLocationZ': (['OBSGEO-Z'], None),
                      'Observation.observation_id': (['OBSID'], None),
                      'Plane.metaRelease': (['RELEASE', 'REL_DATE'], None),
                      'Plane.dataRelease': (['RELEASE', 'REL_DATE'], None),
                      'Plane.product_id': (['RUNID'], None),
                      'Plane.provenance.name': (['XPRVNAME'], None),
                      'Plane.provenance.project': (['ADC_ARCH'], None),
                      'Plane.provenance.producer': (['ORIGIN'], None),
                      'Plane.provenance.reference': (['XREFER'], None),
                      'Plane.provenance.lastExecuted': (['DATE-FTS'], None),
                     }
        # TODO position axis, time axis, polarization axis

        if position_axis:
            self._plan['Chunk.position.coordsys'] = (['RADECSYS', 'RADESYS'], None)
            self._plan['Chunk.position.equinox'] = (['EQUINOX' ,'EPOCH'], None)
            self._plan['Chunk.position.axis.axis1.ctype'] = (['CTYPE{}'.format(position_axis[0])], None)
            self._plan['Chunk.position.axis.axis1.cunit'] = (['CUNIT{}'.format(position_axis[0])], None)
            self._plan['Chunk.position.axis.axis2.ctype'] = (['CTYPE{}'.format(position_axis[1])], None)
            self._plan['Chunk.position.axis.axis2.cunit'] = (['CUNIT{}'.format(position_axis[1])], None)
            self._plan['Chunk.position.axis.error1.syser'] = (['CSYER{}'.format(position_axis[0])], None)
            self._plan['Chunk.position.axis.error1.rnder'] = (['CRDER{}'.format(position_axis[0])], None)
            self._plan['Chunk.position.axis.error2.syser'] = (['CSYER{}'.format(position_axis[1])], None)
            self._plan['Chunk.position.axis.error2.rnder'] = (['CRDER{}'.format(position_axis[1])], None)
            self._plan['Chunk.position.axis.function.cd11'] = (['CD{}_{}'.format(position_axis[0], position_axis[0])], None)
            self._plan['Chunk.position.axis.function.cd12'] = (['CD{}_{}'.format(position_axis[0], position_axis[1])], None)
            self._plan['Chunk.position.axis.function.cd21'] = (['CD{}_{}'.format(position_axis[1], position_axis[0])], None)
            self._plan['Chunk.position.axis.function.cd22'] = (['CD{}_{}'.format(position_axis[1], position_axis[1])], None)
            self._plan['Chunk.position.axis.function.dimension.naxis1'] = (['ZNAXIS{}', 'NAXIS{}'.format(position_axis[0])], None)
            self._plan['Chunk.position.axis.function.dimension.naxis2'] = (['ZNAXIS{}' ,'NAXIS{}'.format(position_axis[1])], None)
            self._plan['Chunk.position.axis.function.refCoord.coord1.pix'] = (['CRPIX{}'.format(position_axis[0])], None)
            self._plan['Chunk.position.axis.function.refCoord.coord1.val'] = (['CRVAL{}'.format(position_axis[0])], None)
            self._plan['Chunk.position.axis.function.refCoord.coord2.pix'] = (['CRPIX{}'.format(position_axis[1])], None)
            self._plan['Chunk.position.axis.function.refCoord.coord2.val'] = (['CRVAL{}'.format(position_axis[1])], None)

        if energy_axis:
            self._plan['Chunk.energy.specsys'] = (['SPECSYS'], None,)
            self._plan['Chunk.energy.ssysobs'] = (['SSYSOBS'], None)
            self._plan['Chunk.energy.restfrq'] = (['RESTFRQ'], None)
            self._plan['Chunk.energy.restwav'] = (['RESTWAV'], None)
            self._plan['Chunk.energy.velosys'] = (['VELOSYS'], None)
            self._plan['Chunk.energy.zsource'] = (['ZSOURCE'], None)
            self._plan['Chunk.energy.ssyssrc'] = (['SSYSSRC'], None)
            self._plan['Chunk.energy.velang'] = (['VELANG'], None)
            self._plan['Chunk.energy.axis.axis.ctype'] = (['CTYPE{}'.format(energy_axis)], None)
            self._plan['Chunk.energy.axis.axis.cunit'] = (['CUNIT{}'.format(energy_axis)], None)
            self._plan['Chunk.energy.axis.error.syser'] = (['CSYER{}'.format(energy_axis)], None)
            self._plan['Chunk.energy.axis.error.rnder'] = (['CRDER{}'.format(energy_axis)], None)
            self._plan['Chunk.energy.axis.function.naxis'] = (['NAXIS{}'.format(energy_axis)], None)
            self._plan['Chunk.energy.axis.function.delta'] = (['CDELT{}'.format(energy_axis)], None)
            self._plan['Chunk.energy.axis.function.refCoord.pix'] = (['CRPIX{}'.format(energy_axis)], None)
            self._plan['Chunk.energy.axis.function.refCoord.val'] = (['CRVAL{}'.format(energy_axis)], None)

        if polarization_axis:
            self._plan['Chunk.polarization.axis.axis.ctype'] = (['CTYPE{}'.format(polarization_axis)], None)
            self._plan['Chunk.polarization.axis.axis.cunit'] = (['CUNIT{}'.format(polarization_axis)], None)
            self._plan['Chunk.polarization.axis.function.naxis'] = (['NAXIS{}'.format(polarization_axis)], None)
            self._plan['Chunk.polarization.axis.function.delta'] = (['CDELT{}'.format(polarization_axis)], None)
            self._plan['Chunk.polarization.axis.function.refCoord.pix'] = (['CRPIX{}'.format(polarization_axis)], None)
            self._plan['Chunk.polarization.axis.function.refCoord.val'] = (['CRVAL{}'.format(polarization_axis)], None)

        if time_axis:
            self._plan['Chunk.time.exposure'] = (['EXPTIME', 'INTTIME'], None)
            self._plan['Chunk.time.timesys'] = (['TIMESYS'], None)
            self._plan['Chunk.time.trefpos'] = (['TREFPOS'], None)
            self._plan['Chunk.time.mjdref'] = (['MJDREF'], None)
            self._plan['Chunk.time.axis.axis.ctype'] = (['CTYPE{}'.format(time_axis)], None)
            self._plan['Chunk.time.axis.axis.cunit'] = (['CUNIT{}'.format(time_axis)], None)
            self._plan['Chunk.time.axis.error.syser'] = (['CSYER{}'.format(time_axis)], None)
            self._plan['Chunk.time.axis.error.rnder'] = (['CRDER{}'.format(time_axis)], None)
            self._plan['Chunk.time.axis.function.naxis'] = (['NAXIS{}'.format(time_axis)], None)
            self._plan['Chunk.time.axis.function.delta'] = (['CDELT{}'.format(time_axis)], None)
            self._plan['Chunk.time.axis.function.refCoord.pix'] = (['CRPIX{}'.format(time_axis)], None)
            self._plan['Chunk.time.axis.function.refCoord.val'] = (['CRVAL{}'.format(time_axis)], None)

        self._extensions = {}

        # contains the standard WCS keywords in the FITS file expected by the
        # astropy.WCS package.
        self._wcs_std = {
            'Chunk.naxis': 'ZNAXIS,NAXIS',
            'Chunk.position.coordsys': 'RADECSYS,RADESYS',
            'Chunk.position.equinox': 'EQUINOX,EPOCH',
            'Chunk.energy.specsys': 'SPECSYS',
            'Chunk.energy.ssysobs': 'SSYSOBS',
            'Chunk.energy.restfrq': 'RESTFRQ',
            'Chunk.energy.restwav': 'RESTWAV',
            'Chunk.energy.velosys': 'VELOSYS',
            'Chunk.energy.zsource': 'ZSOURCE',
            'Chunk.energy.ssyssrc': 'SSYSSRC',
            'Chunk.energy.velang': 'VELANG',
            'Chunk.time.exposure': 'EXPTIME,INTTIME',
            'Chunk.time.timesys': 'TIMESYS',
            'Chunk.time.trefpos': 'TREFPOS',
            'Chunk.time.mjdref': 'MJDREF'
        }

        if position_axis:
            self._wcs_std['Chunk.position.axis.axis1.ctype'] = 'CTYPE{}'.format(position_axis[0])
            self._wcs_std['Chunk.position.axis.axis1.cunit'] = 'CUNIT{}'.format(position_axis[0])
            self._wcs_std['Chunk.position.axis.axis2.ctype'] = 'CTYPE{}'.format(position_axis[1])
            self._wcs_std['Chunk.position.axis.axis2.cunit'] = 'CUNIT{}'.format(position_axis[1])
            self._wcs_std['Chunk.position.axis.error1.syser'] = 'CSYER{}'.format(position_axis[0])
            self._wcs_std['Chunk.position.axis.error1.rnder'] = 'CRDER{}'.format(position_axis[0])
            self._wcs_std['Chunk.position.axis.error2.syser'] = 'CSYER{}'.format(position_axis[1])
            self._wcs_std['Chunk.position.axis.error2.rnder'] = 'CRDER{}'.format(position_axis[1])
            self._wcs_std['Chunk.position.axis.function.cd11'] = 'CD{}_{}'.format(position_axis[0], position_axis[0])
            self._wcs_std['Chunk.position.axis.function.cd12'] = 'CD{}_{}'.format(position_axis[0], position_axis[1])
            self._wcs_std['Chunk.position.axis.function.cd21'] = 'CD{}_{}'.format(position_axis[1], position_axis[0])
            self._wcs_std['Chunk.position.axis.function.cd22'] = 'CD{}_{}'.format(position_axis[1], position_axis[1])
            self._wcs_std['Chunk.position.axis.function.dimension.naxis1'] = 'NAXIS{}'.format(position_axis[0])
            self._wcs_std['Chunk.position.axis.function.dimension.naxis2'] = 'NAXIS{}'.format(position_axis[1])
            self._wcs_std['Chunk.position.axis.function.refCoord.coord1.pix'] = 'CRPIX{}'.format(position_axis[0])
            self._wcs_std['Chunk.position.axis.function.refCoord.coord1.val'] = 'CRVAL{}'.format(position_axis[0])
            self._wcs_std['Chunk.position.axis.function.refCoord.coord2.pix'] = 'CRPIX{}'.format(position_axis[1])
            self._wcs_std['Chunk.position.axis.function.refCoord.coord2.val'] = 'CRVAL{}'.format(position_axis[1])

        if energy_axis:
            self._wcs_std['Chunk.energy.axis.axis.ctype'] = 'CTYPE{}'.format(energy_axis)
            self._wcs_std['Chunk.energy.axis.axis.cunit'] = 'CUNIT{}'.format(energy_axis)
            self._wcs_std['Chunk.energy.axis.error.syser'] = 'CSYER{}'.format(energy_axis)
            self._wcs_std['Chunk.energy.axis.error.rnder'] = 'CRDER{}'.format(energy_axis)
            self._wcs_std['Chunk.energy.axis.function.naxis'] = 'NAXIS{}'.format(energy_axis)
            self._wcs_std['Chunk.energy.axis.function.delta'] = 'CDELT{}'.format(energy_axis)
            self._wcs_std['Chunk.energy.axis.function.refCoord.pix'] = 'CRPIX{}'.format(energy_axis)
            self._wcs_std['Chunk.energy.axis.function.refCoord.val'] = 'CRVAL{}'.format(energy_axis)

        if position_axis:
            self._wcs_std['Chunk.polarization.axis.axis.ctype'] = 'CTYPE{}'.format(polarization_axis)
            self._wcs_std['Chunk.polarization.axis.axis.cunit'] = 'CUNIT{}'.format(polarization_axis)
            self._wcs_std['Chunk.polarization.axis.function.naxis'] = 'NAXIS{}'.format(polarization_axis)
            self._wcs_std['Chunk.polarization.axis.function.delta'] = 'CDELT{}'.format(polarization_axis)
            self._wcs_std['Chunk.polarization.axis.function.refCoord.pix'] = 'CRPIX{}'.format(polarization_axis)
            self._wcs_std['Chunk.polarization.axis.function.refCoord.val'] = 'CRVAL{}'.format(polarization_axis)

        if time_axis:
            self._wcs_std['Chunk.time.axis.axis.ctype'] = 'CTYPE{}'.format(time_axis)
            self._wcs_std['Chunk.time.axis.axis.cunit'] = 'CUNIT{}'.format(time_axis)
            self._wcs_std['Chunk.time.axis.error.syser'] = 'CSYER{}'.format(time_axis)
            self._wcs_std['Chunk.time.axis.error.rnder'] = 'CRDER{}'.format(time_axis)
            self._wcs_std['Chunk.time.axis.function.naxis'] = 'NAXIS{}'.format(time_axis)
            self._wcs_std['Chunk.time.axis.function.delta'] = 'CDELT{}'.format(time_axis)
            self._wcs_std['Chunk.time.axis.function.refCoord.pix'] = 'CRPIX{}'.format(time_axis)
            self._wcs_std['Chunk.time.axis.function.refCoord.val'] = 'CRVAL{}'.format(time_axis)

        # for a quick lookup of keywords referenced by the plan
        self._inverse_plan = {}
        for key, value in self._plan.items():
            if isinstance(value, tuple):
                for ii in value[0]:
                    self._inverse_plan[ii] = key

        # for a quick lookup of config reference values
        if user_supplied_config:
            self._inverse_user_supplied_config = \
                {v: k for k, v in user_supplied_config.items()}

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
            ['{} = {}'.format(key, '{}, default = {}'.format(
                src[key][0], src[key][1])
            if isinstance(src[key], tuple)
            else src[key])
             for key in ObservationBlueprint._CAOM2_ELEMENTS
             if key in src])

    def set(self, caom2_element, value, extension=None):
        """
        Sets the value associated with an element in the CAOM2 model. Value
        cannot be a tuple.
        :param caom2_element: CAOM2 element
        :param value: new value of the CAOM2 element
        :param extension: extension number (used only for Chunk elements)
        """
        if caom2_element not in ObservationBlueprint._CAOM2_ELEMENTS:
            raise ValueError(
                '{} not a caom2 element (spelling?).'.format(caom2_element))
        if extension:
            if not caom2_element.startswith('Chunk'):
                raise ValueError(
                    "Extension number refers to Chunk elements only")
            if extension not in self._extensions:
                self._extensions[extension] = {}
            self._extensions[extension][caom2_element] = value
        else:
            self._plan[caom2_element] = value

    def set_fits_attribute(self, caom2_element, fits_attribute_list, extension=None):
        """
        Associates a CAOM2 element with a FITS attribute
        :param caom2_element:
        :param fits_attribute_list:
        :param extension: extension number (used only for Chunk elements)
        """
        if caom2_element not in ObservationBlueprint._CAOM2_ELEMENTS:
            raise ValueError(
                '{} not a caom2 element (spelling?).'.format(caom2_element))
        if extension:
            if not caom2_element.startswith('Chunk'):
                raise ValueError(
                    "Extension number refers to Chunk elements only")
            if extension not in self._extensions:
                self._extensions[extension] = {}
            self._extensions[extension][caom2_element] = (fits_attribute_list, None)
        else:
            self._plan[caom2_element] = (fits_attribute_list, None)

    def add_fits_attribute(self, caom2_element, fits_attribute, extension=None):
        """
        Adds a FITS attribute in the list of other FITS attributes associated
        with an caom2 element.
        :param caom2_element:
        :param fits_attribute:
        :param extension: extension number (used only for Chunk elements)
        :raises AttributeError if the caom2 element has already an associated
        value or KeyError if the caom2 element does not exists.
        """
        if caom2_element not in ObservationBlueprint._CAOM2_ELEMENTS:
            raise ValueError(
                '{} not a caom2 element (spelling?).'.format(caom2_element))
        if extension:
            if not caom2_element.startswith('Chunk'):
                raise ValueError(
                    "Extension number refers to Chunk elements only")
            if extension not in self._extensions:
                raise AttributeError(
                    'No extension {} in the blueprint'.format(extension))
            else:
                if caom2_element in self._extensions[extension]:
                    if isinstance(self._extensions[extension][caom2_element], tuple):
                        self._extensions[extension][caom2_element][0].insert(0, fits_attribute)
                    else:
                        raise AttributeError(
                            'No FITS attributes in extension {} associated with keyword {}'.
                                format(extension, caom2_element))
                else:
                    raise KeyError(
                        'Keyword {} not found in the extension {} of the blueprint'.
                            format(caom2_element, extension))
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

    def set_default(self, caom2_element, default, extension=None):
        """
        Sets the default value of a caom2 element that is associated with FITS
        attributes. If the element does not exist or does not have a list of
        associated FITS attributes, default is set as the associated value
        of the element
        :param caom2_element:
        :param default: default value
        :param extension: extension number (used only for Chunk elements)
        """
        element = self._lookup_element(caom2_element)
        if extension:
            if not element.startswith('Chunk'):
                raise ValueError(
                    "Extension number refers to Chunk elements only")
            if extension not in self._extensions:
                self._extensions[extension] = {}
            if element in self._extensions[extension] and\
                isinstance(self._extensions[extension][element], tuple):
                self._extensions[extension][element] = (self._extensions[extension][element][0], default)
            else:
                # default is the only value
                self._extensions[extension][element] = default
        else:
            if (element in self._plan) and \
                    isinstance(self._plan[element], tuple):
                self._plan[element] = (self._plan[element][0], default)
            else:
                # override the value
                self._plan[element] = default

    def _get(self, caom2_element, extension=None):
        """
        Returns the value associated with a CAOM2 element
        :param caom2_element:
        :param extension: extension number (used only for Chunk elements)
        :return: Tuple of the form (list_of_associated_fits_attributes,
        default_value) OR the actual associated value of the CAOM2 element
        """
        element = self._lookup_element(caom2_element)
        if extension:
            if not element.startswith('Chunk'):
                raise ValueError(
                    "Extension number refers to Chunk elements only")
            if (extension in self._extensions) and (element in self._extensions[extension]):
                return self._extensions[extension][element]

        # look in the generic plan
        if element not in self._plan:
            return None
        else:
            return self._plan[element]

    def _lookup_element(self, caom2_element):
        if caom2_element in ObservationBlueprint._CAOM2_ELEMENTS:
            return caom2_element
        elif caom2_element in self._inverse_plan.keys():
            return self._inverse_plan[caom2_element]
        elif caom2_element in self._inverse_user_supplied_config.keys():
            return self._inverse_user_supplied_config[caom2_element]
        else:
            raise ValueError(
                '{} caom2 element not found in the plan (spelling?).'.
                format(caom2_element))


class FitsParser(object):
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

    def __init__(self, src, obs_blueprint=None):
        """
        Ctor
        :param src: List of headers (dictionary of FITS keywords:value) with
        one header for each extension or a FITS input file.
        """
        self.logger = logging.getLogger(__name__)
        if obs_blueprint:
            self.blueprint = obs_blueprint
        else:
            self.blueprint = ObservationBlueprint()
        self._headers = []
        self.parts = 0
        self.file = ''
        if isinstance(src, list):
            # assume this is the list of headers
            self._headers = src
        else:
            # assume file
            self.file = src
            hdulist = fits.open(self.file, memmap=True, lazy_load_hdus=False)
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

            wcs_parser = WcsParser(header, self.file, ii)
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
        :param artifact_uri: the key for finding the artifact to augment
        :param product_id: the key for finding for the plane to augment
        """
        self.logger.debug(
            'Begin CAOM2 observation augmentation for URI {}.'.format(
                artifact_uri))

        assert observation
        assert isinstance(observation, Observation)

        observation.algorithm = self._get_algorithm(observation)

        observation.sequence_number = int(self._get_from_list(
            'Observation.sequenceNumber', index=0, default=-1))
        observation.intent = self._get_from_list('Observation.intent', 0,
                                                 ObservationIntentType.SCIENCE)
        observation.type = self._get_from_list('Observation.type', 0)
        observation.meta_release = self._get_datetime(

            self._get_from_list('Observation.metaRelease', 0))
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

        plane.meta_release = self._get_datetime(self._get_from_list(
            'Plane.metaRelease', index=0, default=None), 's')
        plane.data_release = self._get_datetime(self._get_from_list(
            'Plane.dataRelease', index=0, default=None), 's')
        plane.data_product_type = \
            DataProductType(self._get_from_list('Plane.dataProductType',
                                                index=0,
                                                default=DataProductType.CUBE))
        plane.calibration_level = \
            CalibrationLevel(int_32(
                self._get_from_list('Plane.calibrationLevel',
                                    index=0,
                                    default=CalibrationLevel.CALIBRATED.value)))
        plane.provenance = self._get_provenance()
        plane.metrics = self._get_metrics()

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

    def _get_algorithm(self, obs):
        """
        Create an Algorithm instance populated with available FITS information.
        :return: Algorithm
        """
        self.logger.debug('Begin CAOM2 Algorithm augmentation.')
        # TODO DEFAULT VALUE
        name = self._get_from_list('Observation.algorithm.name', index=0,
                                   current=obs.algorithm.name)
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
            if keywords:
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
        prop_id = self._get_from_list('Observation.proposal.id', index=0)
        pi = self._get_from_list('Observation.proposal.pi', index=0)
        project = self._get_from_list(
            'Observation.proposal.project', index=0)
        title = self._get_from_list('Observation.proposal.title', index=0)
        self.logger.debug('End CAOM2 Proposal augmentation.')
        if prop_id:
            return Proposal(str(prop_id), pi, project, title)
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
        target_type = self._get_from_list('Observation.target.type',
                                          index=0)
        standard = self._cast_as_bool(self._get_from_list(
            'Observation.target.standard', index=0))
        redshift = self._get_from_list('Observation.target.redshift', index=0)
        keywords = self._get_set_from_list('Observation.target.keywords',
                                           index=0)  # TODO
        self.logger.debug('End CAOM2 Target augmentation.')
        if name:
            return Target(str(name), target_type, standard, redshift, keywords)
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
        answer = self._get_from_list('Observation.telescope.geoLocationX',
                                     index=0)
        geo_x = float(answer) if answer else None
        answer = self._get_from_list('Observation.telescope.geoLocationY',
                                     index=0)
        geo_y = float(answer) if answer else None
        answer = self._get_from_list('Observation.telescope.geoLocationZ',
                                     index=0)
        geo_z = float(answer) if answer else None
        keywords = self._get_set_from_list('Observation.telescope.keywords',
                                           index=0)  # TODO
        self.logger.debug('End CAOM2 Telescope augmentation.')
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
                                index=0, default=False)  # TODO

        if seeing or humidity or elevation or tau or wavelength_tau or ambient:
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
        else:
            return None

    def _get_from_list(self, lookup, index, default=None, current=None):
        value = default
        try:
            keywords = self.blueprint._get(lookup)
        except KeyError:
            self.logger.warning(
                'Could not find {!r} in fits2caom2 configuration.'.format(
                    lookup))
            if current:
                self.logger.warning(
                    '{}: using current value of {!r}.'.format(lookup, current))
                value = current
            return value

        if type(keywords) == list:
            for ii in keywords:
                if isinstance(ii, dict):
                    value = ii['default']
                else:
                    value = self.headers[index].get(ii, default)
                if value is None and current:
                    value = current
                    self.logger.debug(
                        '{}: used current value {!r}.'.format(
                            lookup, value))
                else:
                    self.logger.debug(
                        '{}: assigned value {} based on keyword {}.'.format(
                            lookup, value, ii))
                if value is not default:
                    break
            # TODO set the default
        elif type(keywords) == tuple:
            for ii in keywords[0]:
                try:
                    value = self.headers[index].get(ii)
                    self.logger.debug(
                        '{}: assigned value {} based on keyword {}.'.format(
                            lookup, value, ii))
                    break
                except KeyError:
                    pass

            if value is None and current:
                value = current
                self.logger.debug(
                    '{}: used current value {!r}.'.format(
                        lookup, value))

        elif keywords:
            value = keywords
        elif current:
            value = current

        self.logger.debug('{}: value is {}'.format(lookup, value))
        return value

    def _get_set_from_list(self, lookup, index, default=None):
        value = default
        try:
            keywords = self.blueprint._get(lookup)
        except KeyError:
            self.logger.debug(
                'Could not find \'{}\' in fits2caom2 configuration.'.format(
                    lookup))
            return value

        if isinstance(keywords, list):
            for ii in keywords:
                temp = self.headers[index].get(ii, default)
                self.logger.debug(
                    'Assigned value {} based on keyword {}'.format(temp, ii))
                if temp is not default:
                    value = set()
                    for jj in temp.split(','):
                        value.add(jj)
                    break
        elif isinstance(keywords, tuple):
            for ii in keywords[0]:
                if isinstance(ii, dict):
                    value = ii['default']
                elif isinstance(ii, list):
                    for jj in ii:
                        try:
                            value = self.headers[index].get(jj)
                            break
                        except KeyError:
                            pass
                if value is not default:
                    break

        return value

    def _get_provenance(self):
        """
        Create a Provenance instance populated with available FITS information.
        :return: Provenance
        """
        self.logger.debug('Begin CAOM2 Provenance augmentation.')
        name = self._to_str(
            self._get_from_list('Plane.provenance.name', index=0))
        p_version = self._to_str(self._get_from_list('Plane.provenance.version',
                                                     index=0))  # TODO DEFAULT VALUE
        project = self._to_str(
            self._get_from_list('Plane.provenance.project', index=0))
        producer = self._to_str(
            self._get_from_list('Plane.provenance.producer', index=0))
        run_id = self._to_str(
            self._get_from_list('Plane.provenance.runID', index=0))
        reference = self._to_str(
            self._get_from_list('Plane.provenance.reference', index=0))
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
            prov = Provenance(name, p_version, project, producer, run_id,
                              reference, last_executed)
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
        self.logger.debug('End CAOM2 Metrics augmentation.')
        return metrics

    def _get_datetime(self, from_value, default_units=None):
        """
        Ensure datetime values are in MJD. Really. Just not yet.
        :param from_value:
        :return:
        """

        if default_units:
            units = default_units
        else:
            # FITS Time Paper is the source for defaults
            # http://hea-www.cfa.harvard.edu/~arots/TimeWCS/WCSPaperV0.90.pdf
            units = self.headers[0].get('TIMEUNIT', 's')
        try:
            if from_value:
                if units == 'd':
                    return datetime.strptime(from_value, '%Y-%m-%d')
                elif units == 's':
                    return datetime.strptime(from_value, '%Y-%m-%dT%H:%M:%S')
                else:
                    return datetime(1999, 1, 1, 0, 0, 0)  # TODO better
            else:
                # return datetime(1999, 1, 1, 0, 0, 0)  # TODO better
                return None
        except ValueError:
            self.logger.warning('{}'.format(sys.exc_info()[1]))
            # return datetime(1999, 1, 1, 0, 0, 0)  # TODO better
            return None

    def _cast_as_bool(self, from_value):
        """
        Make lower case Java booleans into capitalized python booleans.
        :param from_value: Something that represents a boolean value
        :return: a python boolean value
        """
        result = False
        # so far, these are the only options that are coming in from the
        # config files - may need to add more as more types are experienced
        if from_value == 'false':
            result = False
        elif from_value == 'true':
            result = True
        return result

    def _to_str(self, value):
        return str(value) if value else None


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
        self.logger.debug('Begin Polarization WCS augmentation.')
        assert chunk
        assert isinstance(chunk, Chunk)

        polarization_axis = self._get_axis_index(POLARIZATION_CTYPES)
        if polarization_axis is None:
            self.logger.debug('No WCS Polarization info')
            return

        chunk.polarization_axis = polarization_axis

        naxis = CoordAxis1D(self._get_axis(polarization_axis))
        naxis.function = CoordFunction1D(
            self._sanitize(self.wcs._naxis[polarization_axis]),
            self._sanitize(self.wcs.wcs.cdelt[polarization_axis]),
            RefCoord(self._sanitize(self.wcs.wcs.crpix[polarization_axis]),
                     self._sanitize(self.wcs.wcs.crval[polarization_axis])))
        if not chunk.polarization:
            chunk.polarization = PolarizationWCS(naxis)
        else:
            chunk.polarization.naxis = naxis
        self.logger.debug('End Polarization WCS augmentation.')

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
            elif len(elem) == 0:
                check = self.wcs.wcs.ctype[i]
                if check in keywords:
                    axis = i
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


def load_config(file_name):
    """
    Override CONFIG with externally-supplied values.

    The override file can contain information for more than one input file,
    as well as providing information for different HDUs.

    :param file_name Name of the configuration file to load.
    :return: dict representation of file content.
    """
    d = {}
    extension = 0
    artifact = ''
    ptr = d
    with open(file_name) as file:
        for line in file:
            if line.startswith('?'):
                # file, extension-specific content
                if line.find('#[') == -1:
                    artifact = line.split('?')[1].strip()
                    extension = 0
                    if 'artifacts' not in d.keys():
                        d['artifacts'] = {}
                else:
                    artifact = line.split('?')[1].split('#[')[0].strip()
                    extension = \
                        int(line.split('#[')[1].split(']')[0].strip())
                    logging.debug(
                        'Adding overrides for artifact {} in extension {}'.
                        format(artifact, extension))
                if artifact not in d['artifacts'].keys():
                    d['artifacts'][artifact] = {}
                if extension not in d['artifacts'][artifact].keys():
                    d['artifacts'][artifact][extension] = {}
                ptr = d['artifacts'][artifact][extension]
            elif line.find('=') == -1 or line.startswith('#'):
                continue
            else:
                key, value = line.split('=')
                ptr[key.strip()] = value.strip()
    return d


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


def update_fits_headers(parser, artifact_uri=None, config=None, defaults=None,
                        overrides=None):
    """
    Update the in-memory representation of FITS headers according to defaults
    and/or overrides as configured by the user.
    :param parser: access to FITS keywords as type astropy.wcs.Header - a
    dictionary of FITS keywords, as well as the default parser CONFIG
    :param artifact_uri: Where the overrides come from, and where
    to apply them.
    :param config: Input configuration in a dict.
    :param defaults: FITS header and configuration default values in a dict.
    :param overrides: FITS header keyword and configuration default overrides
    in a dict.
    :return: updated headers
    """

    # this_config = _merge_configs(parser.CONFIG, config)
    # parser.set_config(this_config)

    # for lack of better criteria, anything that's all upper case, on the
    # left-hand side of an '=' is assumed to be a FITS keyword. On the right
    # hand side of an '=', it is assumed to be a value. I'm sure that will
    # bite somewhere in the future :)

    if defaults:
        logging.debug('Setting defaults for {}'.format(artifact_uri))
        for key, value in defaults.items():
            parser.blueprint.set_default(key, value)
            logging.debug('{} setting default value to {}'.format(key, value))

    logging.debug('Defaults set for {}. Start overrides.'.format(artifact_uri))

    # if overrides:
    #     _set_overrides(overrides, this_config, parser.headers, 0)
    #     _set_overrides_for_artifacts(
    #         overrides, parser, artifact_uri, this_config)
    #
    # logging.debug('Overrides set for {}.'.format(artifact_uri))
    # return parser


def _set_default_keyword_value_in_blueprint(blueprint, key, value):
    # TODO
    return None


def _set_default_value_in_config(_config, _key, _value):
    key_found = False
    for k, v in _config.items():
        if isinstance(v, list):
            for jj in v:
                logging.warning('list entry is {}, _key is {}'.format(jj, _key))
                if jj == _key:
                    index = v.index(jj)
                    v[index] = {'default': _value}
                    logging.debug('{}: Set default value of {}'.format(
                        k, _value))
                    key_found = True
                    break
        elif v == _key:
            index = _config.index(v)
            _config[index] = {'default': _value}
            logging.debug('{}: Set to default value of {}'.format(k, _value))
            key_found = True
            break

        if key_found:
            break

    if not key_found:
        # one last-ditch attempt? for the case where there's a FITS keyword
        # lookup value, that might not exist in the source file?
        if not _last_ditch_find(_config, _key, _value):
            msg = 'Could not find configuration key for default {!r}'.format(_key)
            logging.warning(msg)
            raise KeyError(msg)


def _set_overrides(_overrides, _config, _headers, _index):
    for ii in _overrides.keys():
        if ii == 'artifacts':
            # the overrides that are part of the 'artifacts' dict entry are
            # handled in a separate function
            continue

        if ii.isupper() and ii.find('.') == -1:
            _set_override_keyword_value_in_header(_headers, [ii],
                                                  _overrides[ii], _index)

        else:
            # config key, add as a default if it's not found in the
            # configuration
            keywords = _find_fits_keyword_in_config(_config, ii)
            if len(keywords) > 0:

                for jj, value in enumerate(keywords):
                    if value.find('.') == -1 and value.isupper():
                        _set_override_keyword_value_in_header(
                            _headers, keywords, _overrides[ii], _index)
                    else:
                        _set_override_value_in_config(
                            _config, ii, _overrides[ii])
                    break
            else:
                _set_override_value_in_config(_config, ii,
                                              _overrides[ii])


def _set_overrides_for_artifacts(_overrides, _parser, _uri, _config):
    if 'artifacts' in _overrides.keys() and _uri in _overrides['artifacts']:
        logging.debug(
            'Found extension overrides for URI {}. Update headers accordingly.'.
            format(_uri))
        if len(_overrides['artifacts'][_uri]) > len(_parser.headers):
            msg = 'Disconnect in {} and overrides.'.format(_uri)
            logging.error(msg)
            raise IndexError(msg)

        for ii in _overrides['artifacts'][_uri]:
            _set_overrides(
                _overrides['artifacts'][_uri][ii], _config, _parser.headers, ii)
    return


def _set_override_value_in_config(_config, _key, _value):
    key_not_found = True
    for k, v in _config.items():
        if isinstance(v, list):
            for jj in v:
                if jj == _key:
                    _config[k] = _set_value(k, _value)
                    key_not_found = False
                    break
        elif v == _key:
            _config[k] = _set_value(k, _value)
            key_not_found = False
            break
    # if a value does not already exist in a configuration, add that
    # value to the configuration, giving over-ride values the effect
    # of appending undefined values
    if key_not_found:
        _config[_key] = _value
        logging.debug('{}: Add override value of {} to configuration.'.format(
            _key, _config[_key]))
    return


def _set_value(_here, _to):
    logging.debug('{}: Set override value of {}'.format(
        _here, _to))
    if len(_to) == 0:
        return None
    else:
        return _to


def _find_fits_keyword_in_config(_config, _key):
    keywords = []
    for k, v in _config.items():
        if _key == k:
            if len(v) > 0:
                keywords = v
            break
    return keywords


def _last_ditch_find(_config, _key, _value):
    found = False
    for ii in _config.keys():
        # the second clause is for the case where plane.dataProductType
        # is part of the defaults
        if ii.endswith(_key) or ii.lower() == _key.lower():
            if isinstance(_config[ii], list):
                _config[ii].append({'default': _value})
            else:
                _config[ii] = [{'default': _value}]
            found = True
            break

    return found

def _set_default_keyword_value_in_header(_headers, _keys, _value):
    """
    The configuration provides a FITS keyword, so modify the headers to have
    the default value for the keyword.
    :param _headers: Metadata from the FITS files.
    :param _keys: The FITS keyword list.
    :param _value: What value to give to the keyword.
    :return:
    """

    # set will append if the keyword doesn't exist, or update an existing value

    if _value.find('{') == -1:
        # the default value does not contain index markup, add the
        # default value only to the first header
        logging.debug(
            'Set header {} to default value of {} in extension 0.'.format(
                _keys[0], _value))
        _headers[0].set(_keys[0], _value, 'fits2caom2 set value')
    else:
        # the default value contains index markup, add the default value
        # to all the headers
        for ii, header in enumerate(_headers):
            logging.debug(
                'Set header {} to default value of {} in extension {}'.format(
                    _keys[0], _value, ii))
            header.set(_keys[0], _value, 'fits2caom2 set value')


def _set_override_keyword_value_in_header(_headers, _keys, _value, _index):
    for key in _keys:
        if _value.find('{') == -1:
            # the default value does not contain index markup
            # if key in _headers[_index].keys():
            logging.debug(
                'Set {} to override value of {} in HDU {}.'.format(
                    key, _value, _index))
            _headers[_index].set(
                key, _value, 'Updated HDU {} value'.format(_index))
            break
        else:
            # the default value contains index markup, check all the headers
            for ii, header in enumerate(_headers):
                # if key in header.keys():
                logging.debug(
                    'Set {} to override value of {} in HDU {}'.format(
                        key, _value, ii))
                header.set(key, _value, 'Updated HDU {} value'.format(ii))
            break


def _merge_configs(default_config, input_config):
    merged = None
    if input_config and default_config:
        logging.debug('Merging two separate configurations into single map.')
        merged = default_config.copy()
        # if there was only a need to support python 3.5+, do it like this:
        # return {**x, **y}

        # input_config takes precedence

        # make the input config values into lists, because that's how the
        # default config is done, that's why and also because of the default
        # value of 'ZNAXIS, NAXIS' for Chunk.naxis

        for ii in input_config.keys():
            if not isinstance(input_config[ii], list):
                values = input_config[ii].split(',')
                input_config[ii] = values

        merged.update(input_config)
    elif input_config and not default_config:
        merged = input_config.copy()
    elif default_config and not input_config:
        merged = default_config.copy()
    return merged


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
                        type=argparse.FileType('w'),
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
    elif args.debug:
        logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)
    elif args.quiet:
        logging.basicConfig(level=logging.ERROR, stream=sys.stdout)
    else:
        logging.basicConfig(level=logging.WARN, stream=sys.stdout)

    handler = logging.StreamHandler()
    handler.setFormatter(DispatchingFormatter({
        'caom2utils.fits2caom2.WcsParser': logging.Formatter(
            '%(levelname)s:%(name)-12s:HDU:%(hdu)-2s:%(message)s'),
        'astropy': logging.Formatter(
            '%(levelname)s:%(name)-12s:HDU:%(hdu)-2s:%(message)s')
    },
        logging.Formatter('%(levelname)s:%(name)-12s:%(message)s')
    ))
    logging.getLogger().addHandler(handler)

    if args.local and (len(args.local) != len(args.fileURI)):
        sys.stderr.write(('number of local arguments not the same with file '
                          'URIs ({} vs {})').format(len(args.local),
                                                    args.fileURI))
        sys.exit(-1)

    config = None
    if args.config:
        config = load_config(args.config)
        logging.debug('Apply configuration from {}.'.format(args.config))

    defaults = None
    if args.default:
        defaults = load_config(args.default)
        logging.debug('Apply defaults from {}.'.format(args.default))

    overrides = None
    if args.override:
        overrides = load_config(args.override)
        logging.debug('Apply overrides from {}.'.format(args.override))

    # invoke the appropriate function based on the inputs
    if args.in_obs_xml:
        # append to existing observation
        reader = ObservationReader(validate=True)
    else:
        obs = SimpleObservation(collection=args.observation[0],
                                observation_id=args.observation[1],
                                algorithm=Algorithm('EXPOSURE'))  # TODO

    if args.productID not in obs.planes.keys():
        obs.planes.add(Plane(product_id=str(args.productID)))

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
                    Artifact(uri=str(uri),
                             product_type=ProductType.SCIENCE,
                             release_type=ReleaseType.DATA))
            parser = FitsParser(headers)

        # TODO update_fits_headers(parser, uri, config, defaults, overrides)
        parser.augment_observation(observation=obs, artifact_uri=uri,
                                   product_id=plane.product_id)

    writer = ObservationWriter()
    if args.out_obs_xml:
        writer.write(obs, args.out_obs_xml)
    else:
        sys.stdout.flush()
        writer.write(obs, sys.stdout)

    logging.info("DONE")

