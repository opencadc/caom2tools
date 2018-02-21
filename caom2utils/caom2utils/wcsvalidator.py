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

# from astropy.wcs import WCS
# from cadcutils import version
from caom2 import Artifact, Part, Chunk, Plane, Observation, CoordError
# from caom2 import SpectralWCS,CoordAxis1D, Axis, CoordFunction1D, RefCoord
# from caom2 import shape
from astropy.wcs import Wcsprm
from caom2utils import TimeUtil, EnergyUtil, ORIGIN
import numpy as np


# from caom2util import TimeUtil

APP_NAME = 'wcsvalidator'

__all__ = ['WcsValidator', 'InvalidWCSError']


class InvalidWCSError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


# WcsValidator class and functions
class WcsValidator():
    """
    WcsValidator: validates WCS coordinates in the artifact passed in

    Example:
        ...

    """
    def __init__(self):
        pass

    def __str__(self):
        pass

    @staticmethod
    def validate_artifact(self, artifact):
        if artifact is not None:
            for p in artifact.parts:
                if p is not None:
                    for c in p.chunks:
                        context = artifact.uri + "[" + p.name + "]:" + c.id + " "
                        WcsValidator.validate_chunk(context, c)

    @staticmethod
    def validate_chunk(context, chunk):
        """
        Validate all WCS in this chunk individually
        """
        WcsValidator.validate_spatial_wcs(chunk.position)
        WcsValidator.validate_spectral_wcs(chunk.energy)
        WcsValidator.validate_temporal_wcs(chunk.time)
        WcsValidator.validate_polarization_wcs(chunk.polarization)



    @staticmethod
    def validate_spatial_wcs(position):
        # position is a SpatialWCS
        error_string = ""
        if position is not None and position.axis is not None:
            try:
                # There's not much that can be validated about range & bounds
                if position.axis.function is not None:
                    fn2D = position.axis.function
                    naxis1_half = float(fn2D.dimension.naxis1/2)
                    naxis2_half = float(fn2D.dimension.naxis2/2)

                    wcsprm = Wcsprm()
                    coord_array = np.array([[naxis1_half, naxis2_half]])
                    sky_transform = wcsprm.p2s(coord_array, ORIGIN)
                    pix_transform = wcsprm.s2p(sky_transform['world'], ORIGIN)

                    transformed_coords = pix_transform['pixcrd']

                    if not (transformed_coords[0][0] == naxis1_half and transformed_coords[0][1] == naxis2_half):
                        error_string = "Could not transform centre coordinate"

            except Exception as e:
                error_string = repr(e)

            if len(error_string) > 0:
                raise InvalidWCSError("Invalid SpatialWCS: " + error_string + ": " + str(position))

    @staticmethod
    def check_transform(coords):
        # Coords is a shape.Subinterval
        wcsprm = Wcsprm()
        coord_array = np.array([[coords.lower, coords.upper]])
        sky_transform = wcsprm.p2s(coord_array, ORIGIN)
        pix_transform = wcsprm.s2p(sky_transform['world'], ORIGIN)
        transformed_coords = pix_transform['pixcrd']

        if not (transformed_coords[0][0] == coords.lower and transformed_coords[0][1] == coords.upper):
            raise ValueError("Could not transform coordinates pixel to sky, sky to pixel")

    @staticmethod
    def validate_spectral_wcs(energy):
        error_msg = ""
        if energy is not None:
            try:
                energyAxis = energy.axis
                si = None
                energy_util = EnergyUtil()

                if energyAxis.range is not None:
                    si = energy_util.range1d_to_interval(energy, energyAxis.range)
                    WcsValidator.check_transform(si)

                if energyAxis.bounds is not None:
                    for tile in energyAxis.bounds.samples:
                        si = energy_util.range1d_to_interval(energy, tile)
                        WcsValidator.check_transform(si)

                if energyAxis.function is not None:
                    si = energy_util.function1d_to_interval(energy)
                    WcsValidator.check_transform(si)

            except Exception as ex:
                error_msg = repr(ex)

            if len(error_msg) > 0:
                raise InvalidWCSError("Invalid Spectral WCS: " + error_msg + ": " + str(energy))

    @staticmethod
    def validate_temporal_wcs(time):
        error_msg = ""
        if time is not None:
            subinterval = None
            try:
                timeutil = TimeUtil()
                time_axis = time.axis

                if time_axis.range is not None:
                    subinterval = timeutil.range1d_to_interval(time, time_axis.range)

                if time_axis.bounds is not None:
                    for cr in time_axis.bounds.samples:
                        subinterval = timeutil.range1d_to_interval(time, cr)

                if time_axis.function is not None:
                        subinterval = timeutil.function1d_to_interval(time, time_axis.function)

            except Exception as e:
                error_msg = repr(e)

            if len(error_msg) > 0:
                raise InvalidWCSError("Invalid Temporal WCS: " + error_msg + ": " + str(time))


    @staticmethod
    def validate_polarization_wcs(polarization):
        return True




