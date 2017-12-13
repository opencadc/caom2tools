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

from argparse import ArgumentParser
from builtins import str

import math
from astropy.wcs import WCS
from astropy.io import fits
from cadcutils import version
from caom2 import Artifact, Part, Chunk, CoordError
from caom2 import SpectralWCS, CoordAxis1D, Axis, CoordFunction1D, RefCoord
from caom2 import SpatialWCS, Dimension2D, Coord2D, CoordFunction2D
from caom2 import CoordAxis2D, PolarizationWCS, TemporalWCS
import logging
import os
import sys

APP_NAME = 'fits2caom2'

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


class LoggingFilter(logging.Filter):
    """Add the HDU number to logging messages as a default."""

    def __init__(self):
        super(LoggingFilter, self).__init__()
        self._hdu_for_log = -1

    def filter(self, record):
        record.hdu = self._hdu_for_log
        return True

    def hdu_for_log(self, value):
        self._hdu_for_log = value


class FitsParser(object):
    """
    Parses a FITS file and extracts the CAOM2 related information which can
    be used to augment an existing CAOM2 observation, plane or artifact. If
    the FITS file is not provided, a list of dictionaries (FITS keyword=value)
    with a dictionary for each "extension" need to be provided instead.

    The WCS-related keywords of the FITS file are consumed by the astropy.wcs
    package which might display warnings with regards to compliance.

    Example 1:
    parser = FitsParser(file = '/staging/700000o.fits.gz')
    ...
    # customize parser.hdulist by deleting, changing or adding attributes

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
    parser = FitsParser()

    headers = [] # list of dictionaries headers
    # populate headers

    parser.hdulist = headers

    parser.augment_observation(obs)
    ...

    """

    def __init__(self,
                 file):
        """
        Ctor
        :param file: FITS file
        """

        _configure_logging(self)

        self.file = file
        self.filename = os.path.basename(file)
        self.hdulist = fits.open(file, memmap=True, lazy_load_hdus=False)
        self.hdulist.close()
        self.parts = len(self.hdulist)

    def augment_artifact(self, artifact):
        """
        Augments a given CAOM2 artifact with available FITS information
        :param artifact: existing CAOM2 artifact to be augmented
        """
        self.logger.debug(
            'Begin CAOM2 artifact augmentation for {} with {} HDUs.'.format(
                self.file, self.parts))

        assert artifact
        assert isinstance(artifact, Artifact)

        for i in range(self.parts):
            hdu = self.hdulist[i]
            ii = str(i)
            self.log_filter.hdu_for_log(i)

            # there is one Part per extension, the name is the extension number
            if hdu.size:
                if ii not in artifact.parts.keys():
                    artifact.parts.add(Part(ii))  # TODO use extension name
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

            wcs_parser = WcsParser(hdu.header, self.filename)
            wcs_parser.augment_position(chunk)
            wcs_parser.augment_energy(chunk)
            wcs_parser.augment_temporal(chunk)
            wcs_parser.augment_polarization(chunk)

        self.logger.debug(
            'End CAOM2 artifact augmentation for {}'.format(self.filename))


class WcsParser(object):
    """
    Parser to augment chunks with positional, temporal, energy and polarization
    information based on the WCS keywords in an extension of a FITS header
    """

    ENERGY_AXIS = 'energy'
    POLARIZATION_AXIS = 'polarization'
    TIME_AXIS = 'time'

    def __init__(self, header, filename):
        """

        :param header: FITS extension header
        :param filename: name of FITS file
        """
        _configure_logging(self)
        self.wcs = WCS(header)
        self.header = header
        self.filename = filename

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
            self.logger.debug(
                'No WCS Energy info for {}'.format(self.filename))
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
            self.logger.warning(
                'No celestial metadata for {}'.format(self.filename))

    def augment_temporal(self, chunk):
        """
        Augments a chunk with temporal WCS information
        :param chunk:
        :return:
        """
        self.logger.debug('Begin TemporalWCS augmentation.')
        assert chunk
        assert isinstance(chunk, Chunk)

        time_axis = self._get_axis_index(TIME_KEYWORDS)

        if time_axis is None:
            self.logger.warning('No WCS Time for {}'.format(self.filename))
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

        # TODO need to set the values for the keywords in the test file headers
        chunk.time.exposure = self.header.get('EXPTIME', 0.02)
        chunk.time.resolution = self.header.get('TODO', 0.02)
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
        self.logger.debug('Begin TemporalWCS augmentation.')
        assert chunk
        assert isinstance(chunk, Chunk)

        polarization_axis = self._get_axis_index(POLARIZATION_CTYPES)
        if polarization_axis is None:
            self.logger.debug(
                'No WCS Polarization info for {}'.format(self.filename))
            return

        chunk.polarization_axis = polarization_axis

        naxis = CoordAxis1D(Axis(str(self.wcs.wcs.ctype[polarization_axis]),
                                 str(self.wcs.wcs.cunit[polarization_axis])))
        naxis.function = CoordFunction1D(
            self._sanitize(self.wcs._naxis[polarization_axis]),
            self._sanitize(self.wcs.wcs.cdelt[polarization_axis]),
            RefCoord(self._sanitize(self.wcs.wcs.crpix[polarization_axis]),
                     self._sanitize(self.wcs.wcs.crval[polarization_axis])))
        if not chunk.polarization:
            chunk.polarization = PolarizationWCS(naxis)
        else:
            chunk.polarization.naxis = naxis

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
        Sanitazes values from FITS to caom2
        :param value:
        :return:
        """
        if isinstance(value, float) and math.isnan(value):
            return None
        elif not str(value):
            return None  # empty string
        else:
            return value


def _configure_logging(obj):
    obj.logger = logging.getLogger(__name__)
    if not len(obj.logger.handlers):
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            ('%(asctime)s %(name)-12s '
             '%(levelname)-8s HDU:%(hdu)-2d %(message)s'))
        handler.setFormatter(formatter)
        obj.logger.addHandler(handler)

    obj.logger.setLevel(logging.DEBUG)
    obj.log_filter = LoggingFilter()
    obj.logger.addFilter(obj.log_filter)
    logging.getLogger('astropy').addFilter(obj.log_filter)


def main_app():
    parser = ArgumentParser()

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

    parser.add_argument('-o', '--out', dest='out_obs_xml',
                        help='output of augmented observation in XML',
                        required=False)
    parser.add_argument('productID',
                        help='product ID of the plane in the observation')
    parser.add_argument('fileURI', help='URI of a fits file', nargs='+')

    in_group = parser.add_mutually_exclusive_group(required=True)
    in_group.add_argument('-i', '--in', dest='in_obs_xml',
                          help='input of observation to be augmented in XML')
    in_group.add_argument('--observation', nargs=2,
                          help='observation in a collection',
                          metavar=('collection', 'observationID'))

    if len(sys.argv) < 2:
        parser.print_usage(file=sys.stderr)
        sys.stderr.write("{}: error: too few arguments\n".format(APP_NAME))
        sys.exit(-1)
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
