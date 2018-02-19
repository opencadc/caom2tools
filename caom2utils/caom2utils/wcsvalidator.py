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

__all__ = ['WcsValidator', 'validate_polarization_wcs', 'validate_chunk', 'validate_temporal_wcs', 'validate_spatial_wcs', 'validate_spectral_wcs']


# WcsValidator class and functions
class WcsValidator(Artifact):
    """
    WcsValidator: validates WCS coordinates in the artifact passed in

    Example:
        ...

    """
    def __init__(self):
        pass

    def validate(self, artifact):
        # for i, uri in enumerate(args.fileURI):
        if artifact is not None:
            for p in artifact.parts:
                if p is not None:
                    for c in p.chunks:
                        context = artifact.uri + "[" + p.name + "]:" + c.id + " "
                        self.validate_chunk(context, c)

    def __str__(self):
        pass


def validate_chunk(context, chunk):
    """
    Validate all WCS in this chunk individually
    """
    validate_spatial_wcs(chunk.position)
    validate_spectral_wcs(chunk.energy)
    validate_temporal_wcs(chunk.time)
    validate_polarization_wcs(chunk.polarization)


def validate_spatial_wcs(position):
    # position is a SpatialWCS
    retval = True
    if position is not None and position.axis is not None:
        if position.axis.function is not None:
            fn2D = position.axis.function
            naxis1_half = fn2D.dimension.naxis1/2
            naxis2_half = fn2D.dimension.naxis2/2
        #     use this as centre point for p2s translation
        #  now into wcsprm...

        # // need 2D array of doubles to pass in to p2s
        wcsprm = Wcsprm()
        coord_array = np.array([[naxis1_half, naxis2_half]])
        # [[] * n for x in xrange(n)]

        sky_transform = wcsprm.p2s(coord_array, ORIGIN)

        pix_transform = wcsprm.s2p(sky_transform['world'], ORIGIN)

        transformed_coords = pix_transform['pixcrd']
        # TODO: not sure this is the right thing to have returned?
        retval = transformed_coords[0][0] == naxis1_half and transformed_coords[0][1] == naxis2_half

        # Errors from Java code that need to be handled
        # Transform.Result tr = transform.sky2pix(coords)
        #     }
        # } catch (NoSuchKeywordException ex) {
        #     throw new IllegalArgumentException(SPATIAL_WCS_VALIDATION_ERROR + ex.getMessage() + " in " + context, ex);
        # } catch (WCSLibRuntimeException ex) {
        #     throw new IllegalArgumentException(SPATIAL_WCS_VALIDATION_ERROR + ex.getMessage() + " in " + context , ex);

    return retval


def validate_spectral_wcs(energy):
    # Not complete. TODO
    retval = True
    if energy is not None:
        energyAxis = energy.axis
        si = None
        energy_util = EnergyUtil()

        if energyAxis.range is not None:
            si = energy_util.range1d_to_interval(energy, energyAxis.range)

        if energyAxis.bounds is not None:
            for tile in energyAxis.bounds.samples:
                si = energy_util.range1d_to_interval(energy, tile)

        if energyAxis.function is not None:
            print("wcsvalidator: energyAxis has a function")
            si = energy_util.function1d_to_interval(energy)

            wcsprm = Wcsprm()
            print(si)
            coord_array = np.array([[si.lower, si.upper]])

            sky_transform = wcsprm.p2s(coord_array, ORIGIN)
            print(sky_transform)
            pix_transform = wcsprm.s2p(sky_transform['world'], ORIGIN)
            print(pix_transform)

            transformed_coords = pix_transform['pixcrd']
            print(transformed_coords)
            print(si)

            # TODO: not sure this is the right thing to have returned?
            retval = transformed_coords[0][0] == si.lower and transformed_coords[0][1] == si.upper

            # Exceptions from the Java code: does this get handled somehow here or not necessary?
            # } catch (NoSuchKeywordException ex) {
            # throw new IllegalArgumentException(SPECTRAL_WCS_VALIDATION_ERROR + ex.getMessage() + " in " + context, ex);

    return retval


def validate_temporal_wcs(time):
    # TODO: what to report if there are no range, bounds, function?
    subinterval = None
    if time is not None:
        timeutil = TimeUtil()
        time_axis = time.axis

        if time_axis.range is not None:
            subinterval = timeutil.range1d_to_interval(time, time_axis.range)

        if time_axis.bounds is not None:
            for cr in time_axis.bounds.samples:
                subinterval = timeutil.range1d_to_interval(time, cr)

        if time_axis.function is not None:
                subinterval = timeutil.function1d_to_interval(time, time_axis.function)

    return subinterval


def validate_polarization_wcs(polarization):
    pass






