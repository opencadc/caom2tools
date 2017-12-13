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
from cadcutils import util, version
from caom2 import Artifact, Part, ProductType, ReleaseType, Chunk, CoordError
from caom2 import SpectralWCS, CoordAxis1D, Axis, CoordFunction1D, RefCoord
from caom2 import SpatialWCS, Dimension2D, Coord2D, CoordFunction2D
from caom2 import CoordAxis2D, PolarizationWCS, TemporalWCS
from caom2 import TemporalWCS
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
    Creates a CAOM2 Artifact object based on the content of the in-coming map,
    which is currently assumed to be the hdulist returned by astropy.

    Assumes an existing observation + plane construct.
    """

    ENERGY_AXIS = 'energy'
    POLARIZATION_AXIS = 'polarization'
    TIME_AXIS = 'time'

    def __init__(self,
                 defaults=None,
                 artifact=None,
                 collection=None):
        self.defaults = defaults
        self.artifact = artifact
        self.collection = collection
        self.wcs = None
        self.chunk = None
        self.header_hdu = None

        # configure logging

        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            '%(asctime)s %(name)-12s %(levelname)-8s HDU:%(hdu)-2d %(message)s')
        handler.setFormatter(formatter)
        self.logger = logging.getLogger()
        self.logger.addHandler(handler)
        self.logger.setLevel(logging.DEBUG)

        self.logfilter = LoggingFilter()
        self.logger.addFilter(self.logfilter)
        logging.getLogger('astropy').addFilter(self.logfilter)

        self.filename = None
        self.parts = None
        self.hdulist = None
        self.header = None
        self.artifact = None

    def augment_artifact(self, file):
        self.filename = os.path.basename(file)
        self.hdulist = fits.open(file, memmap=True, lazy_load_hdus=False)
        self.hdulist.close()
        self.parts = len(self.hdulist)

        self.logger.debug('Begin CAOM2 artifact augmentation for {} with {} HDUs.'.format(file, self.parts))

        if not self.artifact:
            assert not self.collection
            self.artifact = Artifact('ad:{}/{}'.format(self.collection, self.filename),
                                     ProductType.SCIENCE, ReleaseType.DATA)  # TODO

        for i in range(self.parts):
            hdu = self.hdulist[i]
            # print(repr(hdu.header))
            ii = str(i)
            self.logfilter.hdu_for_log(i)

            # there is one Part per extension, the name is the extension number
            if hdu.size:
                if ii not in self.artifact.parts.keys():
                    self.artifact.parts.add(Part(ii))  # TODO use extension name
                    self.logger.debug('Part created for HDU {}.'.format(ii))
            else:
                self.artifact.parts.add(Part(ii))
                self.logger.debug('Create empty part for HDU {}'.format(ii))
                continue
            part = self.artifact.parts[ii]

            # each Part has one Chunk
            if not part.chunks:
                part.chunks.append(Chunk())
            self.chunk = part.chunks[0]

            self.header = hdu.header
            self.wcs = WCS(self.header)
            # print(repr(header))
            # self.wcs.wcs.print_contents()
            self.augment_position()
            self.augment_energy()
            self.augment_temporal()
            self.augment_polarization()

        self.logger.debug('End CAOM2 artifact augmentation for {}'.format(self.filename))
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
        naxis = CoordAxis1D(self.get_axis(energy_axis))
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
        self.logger.debug('Begin Spatial WCS augmentation.')

        if self.wcs.has_celestial:
            self.chunk.positionAxis1, self.chunk.positionAxis2 = self.get_position_axis()
            axis = self.get_spatial_axis(None, self.chunk.positionAxis1 - 1, self.chunk.positionAxis2 - 1)

            # Chunk.position.coordsys = RADECSYS,RADESYS
            # Chunk.position.equinox = EQUINOX,EPOCH
            # Chunk.position.resolution = position.resolution

            if not self.chunk.position:
                self.chunk.position = SpatialWCS(axis)
            else:
                self.chunk.position.axis = axis

            radesys = self.fix_value(self.wcs.celestial.wcs.radesys)
            self.chunk.position.coordsys = None if radesys is None else str(radesys)
            self.chunk.position.equinox = self.fix_value(self.wcs.celestial.wcs.equinox)

            self.logger.debug('End Spatial WCS augmentation.')
        else:
            self.logger.warning('No celestial metadata for {}'.format(self.filename))

    def augment_temporal(self):
        if self.chunk.time:
            raise NotImplementedError

        else:
            self.logger.debug('Begin TemporalWCS augmentation.')
            time_axis = self.get_time_axis()

            if time_axis is None:
                self.logger.warning('No WCS Time for {}'.format(self.filename))
                return

            self.chunk.timeAxis = time_axis
            # set self.chunk.time
            self.get_temporal_axis(time_axis)

            # Chunk.time.exposure = EXPTIME,INTTIME
            # Chunk.time.resolution = time.resolution
            # Chunk.time.timesys = TIMESYS
            # Chunk.time.trefpos = TREFPOS
            # Chunk.time.mjdref = MJDREF

            # TODO need to set the values for the keywords in the test file headers
            self.chunk.time.exposure = self.header.get('EXPTIME', 0.02)
            self.chunk.time.resolution = self.header.get('TODO', 0.02)
            self.chunk.time.timesys = str(self.header.get('TIMESYS', 'UTC'))
            self.chunk.time.trefpos = self.header.get('TREFPOS', None)
            self.chunk.time.mjdref = self.header.get('MJDREF', self.header.get('MJDDATE'))

            self.logger.debug('End TemporalWCS augmentation.')

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

    def get_axis(self, index, over_ctype=None, over_cunit=None):
        aug_ctype = over_ctype if over_ctype is not None else str(self.wcs.wcs.ctype[index])
        aug_cunit = over_cunit if over_cunit is not None else str(self.wcs.wcs.cunit[index])
        aug_axis = Axis(aug_ctype, aug_cunit)
        return aug_axis

    def get_spatial_axis(self, aug_axis, xindex, yindex):
        """Assemble the bits to make the axis parameter needed for SpatialWCS construction."""

        if aug_axis:
            raise NotImplementedError

        else:

            # Chunk.position.axis.axis1.ctype = CTYPE{positionAxis1}
            # Chunk.position.axis.axis1.cunit = CUNIT{positionAxis1}
            # Chunk.position.axis.axis2.ctype = CTYPE{positionAxis2}
            # Chunk.position.axis.axis2.cunit = CUNIT{positionAxis2}
            # Chunk.position.axis.function.dimension.naxis1 = ZNAXIS{positionAxis1},NAXIS{positionAxis1}
            # Chunk.position.axis.function.dimension.naxis2 = ZNAXIS{positionAxis2},NAXIS{positionAxis2}

            aug_dimension = self.get_dimension(None, xindex, yindex)

            aug_ref_coord = Coord2D(self.get_ref_coord(None, xindex), self.get_ref_coord(None, yindex))

            aug_cd11, aug_cd12, aug_cd21, aug_cd22 = self.get_cd(xindex, yindex)

            aug_function = CoordFunction2D(aug_dimension, aug_ref_coord, aug_cd11, aug_cd12, aug_cd21, aug_cd22)

            aug_axis = CoordAxis2D(self.get_axis(xindex),
                                   self.get_axis(yindex),
                                   self.get_coord_error(None, xindex),
                                   self.get_coord_error(None, yindex),
                                   None, None, aug_function)

        return aug_axis


    def get_cd(self, x_index, y_index):

        # Chunk.position.axis.function.cd11 = CD{positionAxis1}_{positionAxis1}
        # Chunk.position.axis.function.cd12 = CD{positionAxis1}_{positionAxis2}
        # Chunk.position.axis.function.cd21 = CD{positionAxis2}_{positionAxis1}
        # Chunk.position.axis.function.cd22 = CD{positionAxis2}_{positionAxis2}

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
            self.logger.warning('Error searching for CD* values {}'.format(sys.exc_info()[1]))
            cd11 = -1.0  # TODO what if neither of these are defined?
            cd12 = -1.0
            cd21 = -1.0
            cd22 = -1.0

        return cd11, cd12, cd21, cd22

    def get_coord_error(self, aug_coord_error, index, over_csyer=None, over_crder=None):
        if aug_coord_error:
            raise NotImplementedError

        else:
            # Chunk.position.axis.error1.syser = CSYER{positionAxis1}
            # Chunk.position.axis.error1.rnder = CRDER{positionAxis1}
            # Chunk.position.axis.error2.syser = CSYER{positionAxis2}
            # Chunk.position.axis.error2.rnder = CRDER{positionAxis2}

            aug_csyer = over_csyer if over_csyer is not None else self.fix_value(self.wcs.wcs.csyer[index])
            aug_crder = over_crder if over_crder is not None else self.fix_value(self.wcs.wcs.crder[index])

            if aug_csyer and aug_crder:
                aug_coord_error = CoordError(aug_csyer, aug_crder)

        return aug_coord_error

    def get_dimension(self, aug_dimension, xindex, yindex):

        if aug_dimension:
            raise NotImplementedError
        else:
            aug_dimension = Dimension2D(self.wcs._naxis[xindex], self.wcs._naxis[yindex])

        return aug_dimension

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

        self.logger.debug('Setting positionAxis1 to {}, positionAxis2 to {}'.format(xaxis, yaxis))

        return xaxis, yaxis

    def get_temporal_axis(self, index):
        self.logger.debug('Begin temporal axis augmentation.')

        # Chunk.time.axis.axis.ctype = CTYPE{timeAxis}
        # Chunk.time.axis.axis.cunit = CUNIT{timeAxis}

        # TODO - need to change headers in test file
        aug_naxis = self.get_axis(index, over_ctype='TIME', over_cunit='d')

        # Chunk.time.axis.bounds.samples = time.samples

        # Chunk.time.axis.error.syser = CSYER{timeAxis}
        # Chunk.time.axis.error.rnder = CRDER{timeAxis}
        # TODO - need to change headers in test file
        aug_error = self.get_coord_error(None, index, over_csyer=1e-07, over_crder=1e-07)

        # Chunk.time.axis.function.naxis = NAXIS{timeAxis}
        # Chunk.time.axis.function.delta = CDELT{timeAxis}
        # Chunk.time.axis.function.refCoord.pix = CRPIX{timeAxis}
        # Chunk.time.axis.function.refCoord.val = CRVAL{timeAxis}
        # TODO - need to change headers in test file
        over_cdelt = 2.31481e-07
        aug_ref_coord = self.get_ref_coord(None, index, over_crpix=0.5, over_crval=56789.4298069)
        aug_cdelt = over_cdelt if not None else self.wcs.wcs.cdelt[index]
        aug_function = CoordFunction1D(self.wcs._naxis[index], aug_cdelt, aug_ref_coord)

        # Chunk.time.axis.range.start.pix = time.range.start.pix
        # Chunk.time.axis.range.start.val = time.range.start.val
        # Chunk.time.axis.range.end.pix = time.range.end.pix
        # Chunk.time.axis.range.end.val = time.range.end.val

        self.chunk.time = TemporalWCS(CoordAxis1D(aug_naxis, aug_error, None, None, aug_function))
        self.logger.debug('End temporal axis augmentation.')


    def get_time_axis(self):
        axis_types = self.wcs.axis_type_names
        if(len(axis_types) >= 3) and (axis_types[2] == ''):
            # TODO - set the metadata so axis_types[2] == 'TIME', in this case
            # what axis_types looks like in this case: ['RA', 'DEC', '']
            return 2
        else:
            index = None
            for i, elem in enumerate(self.wcs.axis_type_names):
                if elem in TIME_KEYWORDS:
                    index = i
                    break
            return index

    def get_ref_coord(self, aug_ref_coord, index, over_crpix=None, over_crval=None):
        if aug_ref_coord:
            raise NotImplementedError
        else:
            # Chunk.position.axis.function.refCoord.coord1.pix = CRPIX{positionAxis1}
            # Chunk.position.axis.function.refCoord.coord1.val = CRVAL{positionAxis1}
            # Chunk.position.axis.function.refCoord.coord2.pix = CRPIX{positionAxis2}
            # Chunk.position.axis.function.refCoord.coord2.val = CRVAL{positionAxis2}

            aug_crpix = over_crpix if over_crpix is not None else self.wcs.wcs.crpix[index]
            aug_crval = over_crval if over_crval is not None else self.wcs.wcs.crval[index]
            aug_ref_coord = RefCoord(aug_crpix, aug_crval)
        return aug_ref_coord

    def fix_value(self, value):
        if isinstance(value, float) and math.isnan(value):
            return None
        elif not str(value):
            return None  # empty string
        else:
            return value

    def convert_argv(self):
        convertedArgv = []
        for arg in sys.argv:
            tempArg = arg.split("=")

    def main_app(self):
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

        parser.add_argument('--dumpconfig', help='output the utype to keyword mapping to the console',
                            action='store_true')
        parser.add_argument('--ignorePartialWCS', help='do not stop and exit upon finding partial WCS',
                            action='store_true')

        parser.add_argument('-o', '--out', dest='out_obs_xml', help='output of augmented observation in XML',
                            required=False)
        in_group = parser.add_mutually_exclusive_group(required=True)
        in_group.add_argument('-i', '--in', dest='in_obs_xml', help='input of observation to be augmented in XML')
        in_group.add_argument('--observation', nargs=2, help='observation in a collection',
                              metavar=('collection', 'observationID'))

        parser.add_argument('--config',
                            help='optional CAOM2 utype to keyword config file to merge with the internal configuration',
                            required=False)
        parser.add_argument('--default', help='file with default values for keywords')
        parser.add_argument('--override', help='file with override values for keywords')
        parser.add_argument('--local', help='list of files in local filesystem (same order as uri)', nargs='+')
        parser.add_argument('--log', help='log file name > (instead of console)')
        parser.add_argument('--keep', help='keep the locally stored files after ingestion', action='store_true')
        parser.add_argument('--test', help='test mode, do not persist to database', action='store_true')
        parser.add_argument('--cert', help='Cert File or Proxy Cert&Key PEM file')

        parser.add_argument('productID', help='product ID of the plane in the observation')
        parser.add_argument('fileURI', help='URI of a fits file', nargs='+')

        myArgv = self.convert_argv()
        args = parser.parse_args(myArgv)
        if len(sys.argv) < 2:
            # correct error message when running python3
            parser.print_usage(file=sys.stderr)
            sys.stderr.write("{}: error: too few arguments\n".format(APP_NAME))
            sys.exit(-1)
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


