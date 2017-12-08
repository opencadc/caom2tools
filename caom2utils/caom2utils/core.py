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

from argparse import ArgumentParser
import logging
import sys

from astropy.io import fits
from astropy.wcs import WCS as AWCS
from cadcutils import util, version
from caom2 import Artifact, Part, ProductType, ReleaseType, Chunk
from caom2 import Axis, SpatialWCS, SpectralWCS, CoordAxis2D, Coord2D, CoordError, CoordFunction2D, Dimension2D
from caom2 import RefCoord

# TODO - check defaults
#
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
APP_NAME = 'fits2caom2'
BOUNDS_DEFAULT = None
COORDSYS_DEFAULT = None
EQUINOX_DEFAULT = None
RANGE_DEFAULT = None
RESOLUTION_DEFAULT = None

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
    hdulist = fits.open(file)
    hdulist.close()
    parts = len(hdulist)

    if not artifact:
        assert not collection
        artifact = Artifact('ad:{}/{}'.format(collection, file),
                            ProductType.SCIENCE, ReleaseType.DATA)  # TODO

    # there is one part per extension, the name is the extension number,
    # and each part has one chunk

    for p in range(parts):
        if str(p) not in artifact.parts.keys():
            artifact.parts.add(Part(str(p)))
        part = artifact.parts[str(p)]
        if not part.chunks:
            part.chunks.append(Chunk())
        chunk = part.chunks[p]
        # augment_energy(chunk, hdulist[0].header)
        augment_position(chunk, hdulist[p].header, AWCS(hdulist[p]), file)
        # print('*******Chunk - {}'.format(chunk))

    return artifact


def augment_position(aug_chunk, header, w, file):

    if aug_chunk.position:
        raise NotImplementedError

    else:

        if w.has_celestial:

            aug_chunk.positionAxis1, aug_chunk.positionAxis2 = get_position_axis(w)

            # Chunk.position.coordsys = RADECSYS,RADESYS
            # Chunk.position.equinox = EQUINOX,EPOCH
            # Chunk.position.resolution = position.resolution

            aug_chunk.position = SpatialWCS(get_axis(None,
                                                     w.celestial,
                                                     aug_chunk.positionAxis1 - 1,
                                                     aug_chunk.positionAxis2 - 1),
                                            avoid_nan(w.celestial.wcs, 'radsys', COORDSYS_DEFAULT),
                                            avoid_nan(w.celestial.wcs, 'equinox', EQUINOX_DEFAULT),
                                            RESOLUTION_DEFAULT)
        else:
            logger.error('No celestial metadata for %s', file)
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


def avoid_nan(wcs, keyword, default=None):
    """astropy sets values to 'nan' if they're undefined in the fits file. caom2 types don't understand this."""
    x = getattr(wcs, keyword, default)
    if str(x) == 'nan':
        x = default
    return x


def avoid_nan_axis(wcs, keyword, index, default=None):
    """astropy sets values to 'nan' if they're undefined in the fits file. caom2 types don't understand this,
    so make the value into None when looking up by array index."""
    x = getattr(wcs, keyword, default)[index]
    if str(x) == 'nan':
        x = default
    return x


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
        aug_range = RANGE_DEFAULT
        aug_bounds = BOUNDS_DEFAULT

        # Chunk.position.axis.function.dimension.naxis1 = ZNAXIS{positionAxis1},NAXIS{positionAxis1}
        # Chunk.position.axis.function.dimension.naxis2 = ZNAXIS{positionAxis2},NAXIS{positionAxis2}

        aug_dimension = Dimension2D(getattr(wcs, '_naxis' + str(xindex + 1)), getattr(wcs, '_naxis' + str(yindex + 1)))

        aug_ref_coord = Coord2D(get_ref_coord(None, wcs.wcs, xindex), get_ref_coord(None, wcs.wcs, yindex))

        aug_cd11, aug_cd12, aug_cd21, aug_cd22 = get_cd(wcs.wcs, xindex, yindex);

        aug_function = CoordFunction2D(aug_dimension, aug_ref_coord, aug_cd11, aug_cd12, aug_cd21, aug_cd22)

        aug_axis = CoordAxis2D(aug_axis1, aug_axis2, aug_error1, aug_error2, aug_range, aug_bounds, aug_function)

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

        aug_csyer = avoid_nan_axis(wcs, 'csyer', index)
        aug_crder = avoid_nan_axis(wcs, 'crder', index)

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


def main_app():
    parser = ArgumentParser()

    parser.description = (
        'Augments an observation with information in one or more fits files.')

    if version.version is not None:
        parser.add_argument('-V', '--version', action='version', version=version)

    log_group = parser.add_mutually_exclusive_group()
    log_group.add_argument('-d', '--debug', action='store_true',
                           help='debug messages')
    log_group.add_argument('-q', '--quiet', action='store_true',
                           help='run quietly')
    log_group.add_argument('-v', '--verbose', action='store_true',
                           help='verbose messages')

    parser.add_argument('-o', '--out', dest='out_obs_xml', help='output of augmented observation in XML', required=False)
    parser.add_argument('productID', help='product ID of the plane in the observation')
    parser.add_argument('fileURI', help='URI of a fits file', nargs='+')

    in_group = parser.add_mutually_exclusive_group(required=True)
    in_group.add_argument('-i', '--in', dest='in_obs_xml', help='input of observation to be augmented in XML')
    in_group.add_argument('--observation', nargs=2, help='observation in a collection',
                          metavar=('collection', 'observationID'))

    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(level=logging.INFO, stream=sys.stdout)
    if args.debug:
        logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)
    else:
        logging.basicConfig(level=logging.WARN, stream=sys.stdout)

    if args.out_obs_xml:
        outObsXml = args.out_obs_xml

    if args.in_obs_xml:
        inObsXml = args.in_obs_xml
    else:
        collection = args.observation[0]
        observationID = args.observation[1]

    fileURIs = args.fileURI

    # invoke the appropriate function based on the inputs


    logging.info("DONE")

if __name__ == '__main__':
    main_app()