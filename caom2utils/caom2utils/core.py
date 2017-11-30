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

""" Defines utilities for working with fits files """

# TODO

# defaults are based on test observation
# http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ui/view/CGPS/MA1_DRAO-ST
#
# plus this:
#
# from cgps.config in git repo wcaom2archive
#
# DATE-OBS                             = NO_RELIABLE_DATE_IN_HEADER
# Plane.provenance.lastExecuted        = DATE-FTS
# Chunk.polarization.axis.axis.cunit   = THIS_HAS_NO_VALID_VALUE
#
# from cgps.default
#
# target.type             = field
# provenance.project      = CGPS
# CUNIT1                  = deg
# CUNIT2                  = deg
# obs.intent              = science
#
# from cgps.defaults
#
# target.classification   = FIELD
# process.version         = 1
# process.out.version     = 1
#
BOUNDS_DEFAULT = None
COORDSYS_DEFAULT = None
CTYPE_DEFAULT = None
CUNIT_DEFAULT = 'deg'
EQUINOX_DEFAULT = None
RANGE_DEFAULT = None
RESOLUTION_DEFAULT = None

# from ca.nrc.cadc.caom2.fits.Ctypes.java
#
POSITION_CTYPES = [
    'RA-',
    'DEC-',
    'GLON-',
    'GLAT-',
    'ELON-',
    'ELAT-',
    'HLON-',
    'HLAT-',
    'SLON-',
    'SLAT-',
    'CUBEFACE-']

POSITION_NAXES = [
    'ZAXIS',
    'NAXIS']

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

import logging

from astropy.io import fits
from astropy.wcs import WCS as awcs
from caom2 import Artifact, Part, ProductType, ReleaseType, Chunk
from caom2 import Axis, SpatialWCS, SpectralWCS, CoordAxis1D, CoordAxis2D, Coord2D, CoordError, CoordFunction2D, Dimension2D
from caom2 import RefCoord

logger = logging.getLogger()

def augment_artifact(artifact, file, collection=None):
    hdulist = fits.open(file)
    hdulist.close()
    parts = len(hdulist)

    if not artifact:
        assert not collection
        artifact = Artifact('ad:{}/{}'.format(collection, file),
                            ProductType.SCIENCE, ReleaseType.DATA) #TODO

    # there is one part per extension, the name is the extension number,
    # and each part has one chunk

    for p in range(parts):
        if str(p) not in artifact.parts.keys():
            artifact.parts.add(Part(str(p)))
        part = artifact.parts[str(p)]
        if not part.chunks:
            part.chunks.append(Chunk())
        chunk = part.chunks[p]
        #augment_energy(chunk, hdulist[0].header)
        augment_position(chunk, hdulist[p].header, awcs(hdulist[p]))
        #print('*******Chunk - {}'.format(chunk))

    return artifact


def augment_position(aug_chunk, header, w):

    if aug_chunk.position:
        raise NotImplementedError

    else:
        # Chunk.position.coordsys = RADECSYS,RADESYS
        # Chunk.position.equinox = EQUINOX,EPOCH
        # Chunk.position.resolution = position.resolution

        print(w.wcs.axis_types)

        #naxis = get_axis_value(header, 'naxis')
        #xtension = get_axis_value(header, 'xtension')
        if((naxis >=2) and (xtension == None)):

            logger.debug('Assuming position WCS exists in this HDU {}', 0)

            aug_chunk.positionAxis1 = 1
            aug_chunk.positionAxis2 = 2
            aug_chunk.position = SpatialWCS(get_axis(None, header, aug_chunk.positionAxis1, aug_chunk.positionAxis2),
                                            get_choice_value(header, ['RADECSYS', 'RADESYS'], default = COORDSYS_DEFAULT),
                                            get_choice_value(header, ['EQUINOX', 'EPOCH'], default = EQUINOX_DEFAULT),
                                            RESOLUTION_DEFAULT)
        else:
            raise NotImplementedError

    return

# from https://github.com/opencadc/caom2/blob/master/fits2caom2/src/main/resources/fits2caom2.config
#
# {positionAxis1} is the index of the positional axis
#
# Chunk.position.axis.range.start.coord1.pix = position.range.start.coord1.pix
# Chunk.position.axis.range.start.coord1.val = position.range.start.coord1.val
# Chunk.position.axis.range.start.coord2.pix = position.range.start.coord2.pix
# Chunk.position.axis.range.start.coord2.val = position.range.start.coord2.val
# Chunk.position.axis.range.end.coord1.pix = position.range.end.coord1.pix
# Chunk.position.axis.range.end.coord1.val = position.range.end.coord1.val
# Chunk.position.axis.range.end.coord2.pix = position.range.end.coord2.pix
# Chunk.position.axis.range.end.coord2.val = position.range.end.coord2.val


def get_axis(aug_axis, header, x_index, y_index):

    if aug_axis:
        raise NotImplementedError

    else:

        # Chunk.position.axis.axis1.ctype = CTYPE{positionAxis1}
        # Chunk.position.axis.axis1.cunit = CUNIT{positionAxis1}
        # Chunk.position.axis.axis2.ctype = CTYPE{positionAxis2}
        # Chunk.position.axis.axis2.cunit = CUNIT{positionAxis2}

        aug_axis1 = Axis(get_axis_value(header, 'ctype', x_index),
                         get_axis_value(header, 'cunit', x_index, CUNIT_DEFAULT))
        aug_axis2 = Axis(get_axis_value(header, 'ctype', y_index),
                         get_axis_value(header, 'cunit', y_index, CUNIT_DEFAULT))

        aug_error1 = get_coord_error(None, header, x_index)
        aug_error2 = get_coord_error(None, header, y_index)
        aug_range = RANGE_DEFAULT
        aug_bounds = BOUNDS_DEFAULT

        # Chunk.position.axis.function.dimension.naxis1 = ZNAXIS{positionAxis1},NAXIS{positionAxis1}
        # Chunk.position.axis.function.dimension.naxis2 = ZNAXIS{positionAxis2},NAXIS{positionAxis2}

        aug_dimension = Dimension2D(get_choice_value(header, POSITION_NAXES, x_index),
                                    get_choice_value(header, POSITION_NAXES, y_index))

        aug_ref_coord = Coord2D(get_ref_coord(None, header, x_index),
                                get_ref_coord(None, header, y_index))

        # Chunk.position.axis.function.cd11 = CD{positionAxis1}_{positionAxis1}
        # Chunk.position.axis.function.cd12 = CD{positionAxis1}_{positionAxis2}
        # Chunk.position.axis.function.cd21 = CD{positionAxis2}_{positionAxis1}
        # Chunk.position.axis.function.cd22 = CD{positionAxis2}_{positionAxis2}

        aug_function = CoordFunction2D(aug_dimension, aug_ref_coord,
                                       get_axis_value(header, 'cdelt', x_index),
                                       get_axis_value(header, 'crota', x_index),
                                       get_axis_value(header, 'crota', y_index),
                                       get_axis_value(header, 'cdelt', y_index))

        aug_axis = CoordAxis2D(aug_axis1, aug_axis2, aug_error1, aug_error2, aug_range, aug_bounds, aug_function)

    return aug_axis


def get_choice_value(header, choices, axis = '', default = None):
    lookup = default

    for ii in choices:
        try:
            lookup = header[ii + str(axis)]
            logger.debug('Found {} with value {}', ii + str(axis), lookup)
            break

        except KeyError:
            pass

    return lookup


def get_coord_error(aug_coord_error, header, axis):
    if aug_coord_error:
        raise NotImplementedError

    else:
        # Chunk.position.axis.error1.syser = CSYER{positionAxis1}
        # Chunk.position.axis.error1.rnder = CRDER{positionAxis1}
        # Chunk.position.axis.error2.syser = CSYER{positionAxis2}
        # Chunk.position.axis.error2.rnder = CRDER{positionAxis2}

        aug_csyer = get_axis_value(header, 'csyer', axis)
        aug_crder = get_axis_value(header, 'crder', axis)
        if aug_csyer and aug_crder:
            aug_coord_error = CoordError(aug_csyer, aug_crder)

    return aug_coord_error


def get_axis_value(header, keyword, axis = '', default = None):
    lookup = default

    try:
        lookup = header[keyword + str(axis)]
        logger.debug('Found {} with value {}', keyword, lookup)
    except KeyError:
        logger.debug('Could not find a value for {}. Using the default of {}', keyword, default)
        lookup = default

    return lookup


def get_ref_coord(aug_ref_coord, header, axis):
    if aug_ref_coord:
        raise NotImplementedError
    else:
        # Chunk.position.axis.function.refCoord.coord1.pix = CRPIX{positionAxis1}
        # Chunk.position.axis.function.refCoord.coord1.val = CRVAL{positionAxis1}
        # Chunk.position.axis.function.refCoord.coord2.pix = CRPIX{positionAxis2}
        # Chunk.position.axis.function.refCoord.coord2.val = CRVAL{positionAxis2}

        aug_ref_coord = RefCoord(get_axis_value(header, 'crpix', axis),
                                 get_axis_value(header, 'crval', axis))
    return aug_ref_coord


def augment_energy(chunk, header):
    specsys = header.get('SPECSYS', None)
    naxis = header.get('NAXIS', 0) or header.get('ZAXIS', 0)
    if not chunk.energy:
        chunk.energy = SpectralWCS(naxis, specsys)


    # determine the index of energy axis
    index = None
    for i in range(1, naxis + 1):
        if header['CTYPE{}'.format(i)] is not None:
           if header['CTYPE{}'.format(i)].split('-')[0].strip() in ENERGY_KEYWORDS:
               index = i
               break

    if not index:
        return

    chunk.ctype = header['CTYPE{}'.format(index).split('-')[0].strip()]
    #chunk.naxis = header['NAXIS{}'.format(index)]


    # Artifact.productType = artifact.productType
    # Artifact.releaseType = artifact.releaseType
    #
    # Part.name = part.name
    # Part.productType = part.productType
    #
    # Chunk.naxis = ZNAXIS, NAXIS
    # Chunk.observableAxis = chunk.observableAxis
    # Chunk.positionAxis1 = getPositionAxis()
    # Chunk.positionAxis2 = getPositionAxis()
    # Chunk.energyAxis = getEnergyAxis()
    # Chunk.timeAxis = getTimeAxis()
    # Chunk.polarizationAxis = getPolarizationAxis()
    #
    # Chunk.observable.dependent.bin = observable.dependent.bin
    # Chunk.observable.dependent.axis.ctype = observable.dependent.ctype
    # Chunk.observable.dependent.axis.cunit = observable.dependent.cunit
    # Chunk.observable.independent.bin = observable.independent.bin
    # Chunk.observable.independent.axis.ctype = observable.independent.ctype
    # Chunk.observable.independent.axis.cunit = observable.independent.cunit
    # Chunk.energy.specsys = SPECSYS
    # Chunk.energy.ssysobs = SSYSOBS
    # Chunk.energy.restfrq = RESTFRQ
    # Chunk.energy.restwav = RESTWAV
    # Chunk.energy.velosys = VELOSYS
    # Chunk.energy.zsource = ZSOURCE
    # Chunk.energy.ssyssrc = SSYSSRC
    # Chunk.energy.velang = VELANG
    # Chunk.energy.bandpassName = bandpassName
    # Chunk.energy.resolvingPower = resolvingPower
    # Chunk.energy.transition.species = energy.transition.species
    # Chunk.energy.transition.transition = energy.transition.transition
    # Chunk.energy.axis.axis.ctype = CTYPE{energyAxis}
    # Chunk.energy.axis.axis.cunit = CUNIT{energyAxis}
    # Chunk.energy.axis.bounds.samples = energy.samples
    # Chunk.energy.axis.error.syser = CSYER{energyAxis}
    # Chunk.energy.axis.error.rnder = CRDER{energyAxis}
    # Chunk.energy.axis.function.naxis = NAXIS{energyAxis}
    # Chunk.energy.axis.function.delta = CDELT{energyAxis}
    # Chunk.energy.axis.function.refCoord.pix = CRPIX{energyAxis}
    # Chunk.energy.axis.function.refCoord.val = CRVAL{energyAxis}
    # Chunk.energy.axis.range.start.pix = energy.range.start.pix
    # Chunk.energy.axis.range.start.val = energy.range.start.val
    # Chunk.energy.axis.range.end.pix = energy.range.end.pix
    # Chunk.energy.axis.range.end.val = energy.range.end.val