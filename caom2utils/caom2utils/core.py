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

import math
from astropy.wcs import WCS
from astropy.io import fits
from caom2 import Artifact, Part, ProductType, ReleaseType, Chunk, CoordError
from caom2 import SpectralWCS, CoordAxis1D, Axis, CoordFunction1D, RefCoord
from caom2 import SpatialWCS, Dimension2D, Coord2D, CoordFunction2D
from caom2 import CoordAxis2D
import logging

ENERGY_KEYWORDS = [
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

logger = logging.getLogger()


def augment_artifact(artifact, file, collection=None):
    hdulist = fits.open(file, memmap=True)
    hdulist.close()
    parts = len(hdulist)

    if not artifact:
        assert not collection
        artifact = Artifact('ad:{}/{}'.format(collection, file),
                            ProductType.SCIENCE, ReleaseType.DATA)  # TODO

    # there is one part per extension, the name is the extension number,
    # and each part has one chunk

    for i in range(parts):
        hdu = hdulist[i]
        if hdu.size:
            if str(i) not in artifact.parts.keys():
                artifact.parts.add(Part(str(i))) #TODO use extension name
        part = artifact.parts[str(i)]
        if not part.chunks:
            part.chunks.append(Chunk())

        chunk = part.chunks[i]
        header = hdulist[i].header
        header['RESTFRQ'] = header['OBSFREQ'] #TODO remove
        header['VELREF'] = 256 #TODO remove
        wcs = WCS(header)
        augment_position(chunk, wcs, file)
        augment_energy(chunk, wcs, file)
    return artifact


def augment_energy(chunk, wcs, file):
    # get the energy axis
    energy_axis = None
    for i, elem in enumerate(wcs.axis_type_names):
        if elem in ENERGY_KEYWORDS:
            energy_axis = i
            break

    if energy_axis is None:
        logger.debug('No WCS Energy info for {}'.format(file))
        return

    chunk.energy_axis = energy_axis

    naxis = CoordAxis1D(Axis(wcs.wcs.ctype[energy_axis],
                             wcs.wcs.cunit[energy_axis].to_string()))
    naxis.function = \
        CoordFunction1D(fix_value(wcs._naxis[energy_axis]), #TODO
                        fix_value(wcs.wcs.cdelt[energy_axis]),
                        RefCoord(fix_value(wcs.wcs.crpix[energy_axis]),
                                 fix_value(wcs.wcs.crval[energy_axis])))
    specsys = wcs.wcs.specsys
    if not chunk.energy:
        chunk.energy = SpectralWCS(naxis, specsys)
    else:
        chunk.energy.naxis = naxis
        chunk.energy.specsys = specsys

    chunk.energy.ssysobs = fix_value(wcs.wcs.ssysobs)
    chunk.energy.restfrq = fix_value(wcs.wcs.restfrq)
    chunk.energy.restwav = fix_value(wcs.wcs.restwav)
    chunk.energy.velosys = fix_value(wcs.wcs.velosys)
    chunk.energy.zsource = fix_value(wcs.wcs.zsource)
    chunk.energy.ssyssrc = fix_value(wcs.wcs.ssyssrc)
    chunk.energy.velang = fix_value(wcs.wcs.velangl)


def augment_position(chunk, wcs, file):
    if wcs.has_celestial:
        chunk.positionAxis1, chunk.positionAxis2 = get_position_axis(wcs)
        axis = get_axis(None, wcs.celestial, chunk.positionAxis1 - 1,
            chunk.positionAxis2 - 1)
        # Chunk.position.coordsys = RADECSYS,RADESYS
        # Chunk.position.equinox = EQUINOX,EPOCH
        # Chunk.position.resolution = position.resolution

        if not chunk.position:
            chunk.position = SpatialWCS(axis)
        else:
            chunk.position.axis = axis

        chunk.position.coordsys = fix_value(wcs.celestial.wcs.radesys)
        chunk.position.equinox = fix_value(wcs.celestial.wcs.equinox)

    else:
        logger.debug('No celestial metadata for {}'.format(file))


def get_axis(aug_axis, wcs, xindex, yindex):
    """Assemble the bits to make the axis parameter needed for SpatialWCS construction."""

    if aug_axis:
        raise NotImplementedError

    else:

        # Chunk.position.axis.axis1.ctype = CTYPE{positionAxis1}
        # Chunk.position.axis.axis1.cunit = CUNIT{positionAxis1}
        # Chunk.position.axis.axis2.ctype = CTYPE{positionAxis2}
        # Chunk.position.axis.axis2.cunit = CUNIT{positionAxis2}

        aug_axis1 = Axis(getattr(wcs.wcs, 'ctype')[xindex], getattr(wcs.wcs, 'cunit')[xindex].name)
        aug_axis2 = Axis(getattr(wcs.wcs, 'ctype')[yindex], getattr(wcs.wcs, 'cunit')[yindex].name)

        aug_error1 = get_coord_error(None, wcs.wcs, xindex)
        aug_error2 = get_coord_error(None, wcs.wcs, yindex)

        # Chunk.position.axis.function.dimension.naxis1 = ZNAXIS{positionAxis1},NAXIS{positionAxis1}
        # Chunk.position.axis.function.dimension.naxis2 = ZNAXIS{positionAxis2},NAXIS{positionAxis2}

        aug_dimension = Dimension2D(getattr(wcs, '_naxis{}'.format(xindex + 1)), getattr(wcs, '_naxis{}'.format(yindex + 1)))

        aug_ref_coord = Coord2D(get_ref_coord(None, wcs.wcs, xindex), get_ref_coord(None, wcs.wcs, yindex))

        aug_cd11, aug_cd12, aug_cd21, aug_cd22 = get_cd(wcs.wcs, xindex, yindex)

        aug_function = CoordFunction2D(aug_dimension, aug_ref_coord, aug_cd11, aug_cd12, aug_cd21, aug_cd22)

        aug_axis = CoordAxis2D(aug_axis1, aug_axis2, aug_error1, aug_error2, None, None, aug_function)

    return aug_axis


def get_cd(wcs, x_index, y_index):

    # Chunk.position.axis.function.cd11 = CD{positionAxis1}_{positionAxis1}
    # Chunk.position.axis.function.cd12 = CD{positionAxis1}_{positionAxis2}
    # Chunk.position.axis.function.cd21 = CD{positionAxis2}_{positionAxis1}
    # Chunk.position.axis.function.cd22 = CD{positionAxis2}_{positionAxis2}

    if wcs.has_cd():
        cd11 = getattr(wcs, 'cd')[x_index][x_index]
        cd12 = getattr(wcs, 'cd')[x_index][y_index]
        cd21 = getattr(wcs, 'cd')[y_index][x_index]
        cd22 = getattr(wcs, 'cd')[y_index][y_index]
    else:
        cd11 = getattr(wcs, 'cdelt')[x_index]
        cd12 = getattr(wcs, 'crota')[x_index]
        cd21 = getattr(wcs, 'crota')[y_index]
        cd22 = getattr(wcs, 'cdelt')[y_index]
    return cd11, cd12, cd21, cd22


def get_coord_error(aug_coord_error, wcs, index):
    if aug_coord_error:
        raise NotImplementedError

    else:
        # Chunk.position.axis.error1.syser = CSYER{positionAxis1}
        # Chunk.position.axis.error1.rnder = CRDER{positionAxis1}
        # Chunk.position.axis.error2.syser = CSYER{positionAxis2}
        # Chunk.position.axis.error2.rnder = CRDER{positionAxis2}

        aug_csyer = fix_value(wcs.csyer[index])
        aug_crder = fix_value(wcs.crder[index])

        if aug_csyer and aug_crder:
            aug_coord_error = CoordError(aug_csyer, aug_crder)

    return aug_coord_error


def get_position_axis(wcs):
    axis_types = wcs.get_axis_types()
    return int(axis_types[0]['number']) + 1, int(axis_types[1]['number']) + 1


def get_ref_coord(aug_ref_coord, wcs, index):
    if aug_ref_coord:
        raise NotImplementedError
    else:
        # Chunk.position.axis.function.refCoord.coord1.pix = CRPIX{positionAxis1}
        # Chunk.position.axis.function.refCoord.coord1.val = CRVAL{positionAxis1}
        # Chunk.position.axis.function.refCoord.coord2.pix = CRPIX{positionAxis2}
        # Chunk.position.axis.function.refCoord.coord2.val = CRVAL{positionAxis2}

        aug_ref_coord = RefCoord(getattr(wcs, 'crpix')[index], getattr(wcs, 'crval')[index])
    return aug_ref_coord



def fix_value(value):
    if isinstance(value, float) and math.isnan(value):
        return None
    elif not str(value):
        return None # empyt string
    else:
        return value
