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

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from astropy.wcs import Wcsprm
from caom2utils.wcs_util import TimeUtil, EnergyUtil, ORIGIN
from . import wcs_util
from .wcs_util import PolarizationWcsUtil
from caom2 import Artifact, Chunk, Interval, Observation, Part, Plane, \
    PolarizationState
import numpy as np
import logging


APP_NAME = 'wcsvalidator'
logger = logging.getLogger(APP_NAME)

__all__ = ['validate_wcs', 'InvalidWCSError']


class InvalidWCSError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def validate_wcs(caom2_entity):
    """
    validate_wcs: validates WCS coordinates in the CAOM2 entity passed in
    :param caom2_entity: a caom2 entity such as an Observation, a Plane
    """
    if caom2_entity is not None:
        if isinstance(caom2_entity, Observation):
            _validate_observation(caom2_entity)
            logger.debug('observation validation succeeded.')
        elif isinstance(caom2_entity, Plane):
            _validate_plane(caom2_entity)
            logger.debug('plane validation succeeded.')
        elif isinstance(caom2_entity, Artifact):
            _validate_artifact(caom2_entity)
            logger.debug('artifact validation succeeded.')
        elif isinstance(caom2_entity, Part):
            _validate_part(caom2_entity)
            logger.debug('part validation succeeded.')
        elif isinstance(caom2_entity, Chunk):
            _validate_chunk(caom2_entity)
            logger.debug('chunk validation succeeded.')
        else:
            raise InvalidWCSError("Not a CAOM2 entity")


def _validate_observation(obs):
    if obs is not None and obs.planes is not None:
        for pkey in obs.planes.keys():
            p = obs.planes[pkey]
            _validate_plane(p)


def _validate_plane(plane):
    if plane is not None and plane.artifacts is not None:
        for akey in plane.artifacts.keys():
            a = plane.artifacts[akey]
            _validate_artifact(a)


def _validate_artifact(artifact):
    if artifact is not None:
        for pkey in artifact.parts.keys():
            _validate_part(artifact.parts[pkey])


def _validate_part(part):
    if part is not None:
        for c in part.chunks:
            _validate_chunk(c)


def _validate_chunk(chunk):
    """
    Validate all WCS in this chunk individually
    """
    _validate_spatial_wcs(chunk.position)
    _validate_spectral_wcs(chunk.energy)
    _validate_temporal_wcs(chunk.time)
    _validate_polarization_wcs(chunk.polarization)


def _validate_spatial_wcs(position):
    # position is a SpatialWCS
    error_string = ""
    if position is not None and position.axis is not None:
        try:
            # There's not much that can be validated about range & bounds
            if position.axis.function is not None:
                fn2D = position.axis.function
                _check_transform(
                    Interval(float(fn2D.dimension.naxis1/2),
                             float(fn2D.dimension.naxis2/2)))
                logger.debug('position_axis.function succeeded.')
        except Exception as e:
            error_string = repr(e)

        if len(error_string) > 0:
            raise InvalidWCSError(
                "Invalid SpatialWCS: {}: {}".format(
                    error_string, str(position)))


def _check_transform(coords):
    # Coords is a shape.Subinterval
    wcsprm = Wcsprm()
    coord_array = np.array([[coords.lower, coords.upper]])
    sky_transform = wcsprm.p2s(coord_array, ORIGIN)
    pix_transform = wcsprm.s2p(sky_transform['world'], ORIGIN)
    transformed_coords = pix_transform['pixcrd']

    if not (transformed_coords[0][0] == coords.lower
            and transformed_coords[0][1] == coords.upper):
        raise ValueError(
            "Could not transform coordinates pixel to sky, sky to pixel")


def _validate_spectral_wcs(energy):
    error_msg = ""
    if energy is not None:
        try:
            energy_axis = energy.axis
            si = None

            if energy_axis.range is not None:
                si = EnergyUtil.range1d_to_interval(energy_axis.range)
                _check_transform(si)
                logger.debug('time_axis.range succeeded.')

            if energy_axis.bounds is not None:
                for tile in energy_axis.bounds.samples:
                    si = EnergyUtil.range1d_to_interval(tile)
                    _check_transform(si)
                logger.debug('time_axis.bounds succeeded.')

            if energy_axis.function is not None:
                si = EnergyUtil.function1d_to_interval(energy)
                _check_transform(si)
                logger.debug('time_axis.function succeeded.')

        except Exception as ex:
            error_msg = repr(ex)

        if len(error_msg) > 0:
            raise InvalidWCSError(
                "Invalid Spectral WCS: {}: {}".format(
                    error_msg, str(energy)))


def _validate_temporal_wcs(time):
    error_msg = ""
    if time is not None:
        try:
            time_axis = time.axis

            if time_axis.range is not None:
                logger.debug('time_axis.range to interval validation.')
                TimeUtil.range1d_to_interval(time, time_axis.range)
                logger.debug('time_axis.range to interval succeeded.')

            if time_axis.bounds is not None:
                logger.debug('time_axis.bounds to interval validation.')
                for cr in time_axis.bounds.samples:
                    TimeUtil.range1d_to_interval(time, cr)
                logger.debug('time_axis.bounds to interval succeeded.')

            if time_axis.function is not None:
                logger.debug('time_axis.function to interval validation.')
                TimeUtil.function1d_to_interval(time, time_axis.function)
                logger.debug('time_axis.function to interval succeeded.')

        except Exception as e:
            error_msg = repr(e)

        if len(error_msg) > 0:
            raise InvalidWCSError(
                "Invalid Temporal WCS: {}: {}".format(
                    error_msg, str(time)))


def _validate_range(a_range):
    keys = PolarizationWcsUtil.get_keys(a_range)
    if keys is not None:
        for key in keys:
            WcsPolarizationState.to_value(key)


def _validate_bounds(bounds):
    sample_ranges = PolarizationWcsUtil.get_ranges_from_bounds(bounds)
    if len(sample_ranges) > 0:
        for srange in sample_ranges:
            for key in srange:
                WcsPolarizationState.to_value(key)


def _validate_function(a_function):
    naxis_range = \
        PolarizationWcsUtil.get_range_from_function(a_function)
    if naxis_range is not None:
        for pix in naxis_range:
            WcsPolarizationState.to_value(
                int(round(wcs_util.pix2val(a_function, pix))))


def _validate_polarization_wcs(polarization_wcs):
    """
    Validates the PolarizationWCS.
    :param polarization_wcs: PolarizationWCS to be validated

    An InvalidWCSError is thrown if the PolarizationWCS is determined
    to be invalid.
    """
    if polarization_wcs is not None:
        try:
            axis = polarization_wcs.axis
            _validate_range(axis.range)
            logger.debug('polarization_axis.range succeeded.')
            _validate_bounds(axis.bounds)
            logger.debug('polarization_axis.bounds succeeded.')
            _validate_function(axis.function)
            logger.debug('polarization_axis.function succeeded.')
        except Exception as e:
            raise InvalidWCSError(
                "Invalid Polarization WCS: {}".format(str(e)))


class WcsPolarizationState():
    """
    A dictionary which maps an integer to a PolarizationState value.
    """
    MAP = {
        1: PolarizationState.I, 2: PolarizationState.Q,
        3: PolarizationState.U, 4: PolarizationState.V,
        5: PolarizationState.POLI, 6: PolarizationState.FPOLI,
        7: PolarizationState.POLA, 8: PolarizationState.EPOLI,
        9: PolarizationState.CPOLI, 10: PolarizationState.NPOLI,
        -1: PolarizationState.RR, -2: PolarizationState.LL,
        -3: PolarizationState.RL, -4: PolarizationState.LR,
        -5: PolarizationState.XX, -6: PolarizationState.YY,
        -7: PolarizationState.XY, -8: PolarizationState.YX}

    @staticmethod
    def to_value(key):
        return WcsPolarizationState.MAP[key]
