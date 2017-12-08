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

from builtins import str

import math
from astropy.wcs import WCS
from astropy.io import fits
from caom2 import Artifact, Part, ProductType, ReleaseType, Chunk, CoordError
from caom2 import SpectralWCS, CoordAxis1D, Axis, CoordFunction1D, RefCoord
from caom2 import SpatialWCS, Dimension2D, Coord2D, CoordFunction2D
from caom2 import CoordAxis2D, PolarizationWCS
import logging

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

POLARIZATION_CTYPES = ['STOKES']

class FitsParser(object):
    """
    Augments a CAOM2 Artifact object with the content of a FITS file.

    Assumes an existing observation + plane construct.

    May add to an existing artifact or create a new artifact from the FITS file.

    May override artifact content with default values.
    """

    ENERGY_AXIS = 'energy'
    POLARIZATION_AXIS = 'polarization'
    TIME_AXIS = 'time'

    def __init__(self, filename,
                 artifact=None,
                 defaults=None,
                 collection=None):
        self.filename = filename
        self.defaults = defaults
        self.artifact = artifact
        self.collection = collection
        self.wcs = None
        self.chunk = None
        self.header_hdu = None
        self.logger = logging.getLogger()

        self.hdulist = fits.open(filename, memmap=True)
        self.hdulist.close()
        self.parts = len(self.hdulist)

        if not self.artifact:
            assert not self.collection
            self.artifact = Artifact('ad:{}/{}'.format(collection, self.filename),
                                ProductType.SCIENCE, ReleaseType.DATA)  # TODO


    def augment_artifact(self):

        # there is one part per extension, the name is the extension number,
        # and each part has one chunk

        for i in range(self.parts):
            hdu = self.hdulist[i]
            ii = str(i)
            if hdu.size:
                if ii not in self.artifact.parts.keys():
                    self.artifact.parts.add(Part(ii))  # TODO use extension name
            part = self.artifact.parts[ii]
            if not part.chunks:
                part.chunks.append(Chunk())

            self.chunk = part.chunks[i]
            header = self.hdulist[i].header
            self.wcs = WCS(header)
            self.augment_position()
            self.augment_energy()
            self.augment_polarization()
        return self.artifact

    def _get_axis(self, keywords):
        axis = None
        for i, elem in enumerate(self.wcs.axis_type_names):
            if elem in keywords:
                axis = i
                break
        return axis

    def augment_energy(self):
        # get the energy axis
        energy_axis = self._get_axis(ENERGY_CTYPES)

        if energy_axis is None:
            self.logger.debug('No WCS Energy info for {}'.format(self.filename))
            return

        self.chunk.energy_axis = energy_axis

        naxis = CoordAxis1D(Axis(str(self.wcs.wcs.ctype[energy_axis]),
                                 str(self.wcs.wcs.cunit[energy_axis])))
        naxis.function = \
            CoordFunction1D(self.fix_value(self.wcs._naxis[energy_axis]),  # TODO
                            self.fix_value(self.wcs.wcs.cdelt[energy_axis]),
                            RefCoord(self.fix_value(self.wcs.wcs.crpix[energy_axis]),
                                     self.fix_value(self.wcs.wcs.crval[energy_axis])))
        specsys = str(self.wcs.wcs.specsys)
        if not self.chunk.energy:
            self.chunk.energy = SpectralWCS(naxis, specsys)
        else:
            self.chunk.energy.naxis = naxis
            self.chunk.energy.specsys = specsys

        self.chunk.energy.ssysobs = self.fix_value(self.wcs.wcs.ssysobs)
        self.chunk.energy.restfrq = self.fix_value(self.wcs.wcs.restfrq)
        self.chunk.energy.restwav = self.fix_value(self.wcs.wcs.restwav)
        self.chunk.energy.velosys = self.fix_value(self.wcs.wcs.velosys)
        self.chunk.energy.zsource = self.fix_value(self.wcs.wcs.zsource)
        self.chunk.energy.ssyssrc = self.fix_value(self.wcs.wcs.ssyssrc)
        self.chunk.energy.velang = self.fix_value(self.wcs.wcs.velangl)


    def augment_position(self):
        if self.wcs.has_celestial:
            self.chunk.positionAxis1, self.chunk.positionAxis2 = self.get_position_axis()
            axis = self.get_axis(None, self.chunk.positionAxis1 - 1, self.chunk.positionAxis2 - 1)

            # Chunk.position.coordsys = RADECSYS,RADESYS
            # Chunk.position.equinox = EQUINOX,EPOCH
            # Chunk.position.resolution = position.resolution

            if not self.chunk.position:
                self.chunk.position = SpatialWCS(axis)
            else:
                self.chunk.position.axis = axis

            self.chunk.position.coordsys = self.fix_value(self.wcs.celestial.wcs.radesys)
            self.chunk.position.equinox = self.fix_value(self.wcs.celestial.wcs.equinox)

        else:
            self.logger.debug('No celestial metadata for {}'.format(self.filename))


    def augment_polarization(self):
        polarization_axis = self._get_axis(POLARIZATION_CTYPES)
        if polarization_axis is None:
            self.logger.debug('No WCS Polarization info for {}'.format(self.filename))
            return

        self.chunk.polarization_axis = polarization_axis

        naxis = CoordAxis1D(Axis(str(self.wcs.wcs.ctype[polarization_axis]),
                                 str(self.wcs.wcs.cunit[polarization_axis])))
        naxis.function = \
            CoordFunction1D(self.fix_value(self.wcs._naxis[polarization_axis]),  # TODO
                            self.fix_value(self.wcs.wcs.cdelt[polarization_axis]),
                            RefCoord(self.fix_value(self.wcs.wcs.crpix[polarization_axis]),
                                     self.fix_value(self.wcs.wcs.crval[polarization_axis])))

        if not self.chunk.polarization:
            self.chunk.polarization = PolarizationWCS(naxis)
        else:
            self.chunk.polarization.naxis = naxis


    def get_axis(self, aug_axis, xindex, yindex):
        """Assemble the bits to make the axis parameter needed for SpatialWCS construction."""

        if aug_axis:
            raise NotImplementedError

        else:

            # Chunk.position.axis.axis1.ctype = CTYPE{positionAxis1}
            # Chunk.position.axis.axis1.cunit = CUNIT{positionAxis1}
            # Chunk.position.axis.axis2.ctype = CTYPE{positionAxis2}
            # Chunk.position.axis.axis2.cunit = CUNIT{positionAxis2}

            aug_axis1 = Axis(str(self.wcs.wcs.ctype[xindex]), str(self.wcs.wcs.cunit[xindex].name))
            aug_axis2 = Axis(str(self.wcs.wcs.ctype[yindex]), str(self.wcs.wcs.cunit[yindex].name))

            aug_error1 = self.get_coord_error(None, xindex)
            aug_error2 = self.get_coord_error(None, yindex)

            # Chunk.position.axis.function.dimension.naxis1 = ZNAXIS{positionAxis1},NAXIS{positionAxis1}
            # Chunk.position.axis.function.dimension.naxis2 = ZNAXIS{positionAxis2},NAXIS{positionAxis2}

            aug_dimension = Dimension2D(self.wcs._naxis[xindex + 1], self.wcs._naxis[yindex + 1])

            aug_ref_coord = Coord2D(self.get_ref_coord(None, xindex), self.get_ref_coord(None, yindex))

            aug_cd11, aug_cd12, aug_cd21, aug_cd22 = self.get_cd(xindex, yindex)

            aug_function = CoordFunction2D(aug_dimension, aug_ref_coord, aug_cd11, aug_cd12, aug_cd21, aug_cd22)

            aug_axis = CoordAxis2D(aug_axis1, aug_axis2, aug_error1, aug_error2, None, None, aug_function)

        return aug_axis


    def get_cd(self, x_index, y_index):

        # Chunk.position.axis.function.cd11 = CD{positionAxis1}_{positionAxis1}
        # Chunk.position.axis.function.cd12 = CD{positionAxis1}_{positionAxis2}
        # Chunk.position.axis.function.cd21 = CD{positionAxis2}_{positionAxis1}
        # Chunk.position.axis.function.cd22 = CD{positionAxis2}_{positionAxis2}

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
        return cd11, cd12, cd21, cd22


    def get_coord_error(self, aug_coord_error, index):
        if aug_coord_error:
            raise NotImplementedError

        else:
            # Chunk.position.axis.error1.syser = CSYER{positionAxis1}
            # Chunk.position.axis.error1.rnder = CRDER{positionAxis1}
            # Chunk.position.axis.error2.syser = CSYER{positionAxis2}
            # Chunk.position.axis.error2.rnder = CRDER{positionAxis2}

            aug_csyer = self.fix_value(self.wcs.wcs.csyer[index])
            aug_crder = self.fix_value(self.wcs.wcs.crder[index])

            if aug_csyer and aug_crder:
                aug_coord_error = CoordError(aug_csyer, aug_crder)

        return aug_coord_error


    def get_position_axis(self):

        # there are two celestial axes, get the applicable indices from the axis_types
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

        return xaxis, yaxis


    def get_ref_coord(self, aug_ref_coord, index):
        if aug_ref_coord:
            raise NotImplementedError
        else:
            # Chunk.position.axis.function.refCoord.coord1.pix = CRPIX{positionAxis1}
            # Chunk.position.axis.function.refCoord.coord1.val = CRVAL{positionAxis1}
            # Chunk.position.axis.function.refCoord.coord2.pix = CRPIX{positionAxis2}
            # Chunk.position.axis.function.refCoord.coord2.val = CRVAL{positionAxis2}

            aug_ref_coord = RefCoord(self.wcs.wcs.crpix[index], self.wcs.wcs.crval[index])
        return aug_ref_coord


    def fix_value(self, value):
        if isinstance(value, float) and math.isnan(value):
            return None
        elif not str(value):
            return None  # empty string
        else:
            return value
