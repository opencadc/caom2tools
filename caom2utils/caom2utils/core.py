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

from datetime import datetime

import math
from astropy.wcs import WCS
from astropy.io import fits
from cadcutils import util, version
from caom2 import Artifact, Part, ProductType, ReleaseType, Chunk, CoordError
from caom2 import SpectralWCS, CoordAxis1D, Axis, CoordFunction1D, RefCoord
from caom2 import SpatialWCS, Dimension2D, Coord2D, CoordFunction2D
from caom2 import CoordAxis2D, PolarizationWCS
from caom2 import TemporalWCS
from caom2 import Observation, Plane, Proposal, Telescope, Instrument, Target, TargetPosition, Environment, Algorithm
from caom2 import Requirements, ObservationIntentType, Status, TargetType, Provenance, Metrics, Quality
from caom2 import DataProductType, CalibrationLevel, PlaneURI, DataQuality
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
                 file,
                 defaults=None,
                 artifact=None,
                 collection=None):
        self.defaults = defaults
        self.artifact = artifact
        self.collection = collection
        self.wcs = None
        self.chunk = None
        self.header = None

        self._configure_logging()

        self.file = file
        self.filename = os.path.basename(file)
        self.hdulist = fits.open(file, memmap=True, lazy_load_hdus=False)
        self.hdulist.close()
        self.parts = len(self.hdulist)

    def augment_artifact(self, artifact_uri):
        self.logger.debug('Begin CAOM2 artifact augmentation for {} with {} HDUs.'.format(self.file, self.parts))

        if not self.artifact:
            assert not self.collection
            self.artifact = Artifact(artifact_uri,
                                     ProductType.SCIENCE, ReleaseType.DATA)  # TODO
            # self.artifact = Artifact('ad:{}/{}'.format(self.collection, self.filename),
            #                          ProductType.SCIENCE, ReleaseType.DATA)  # TODO

        for i in range(self.parts):
            hdu = self.hdulist[i]
            # print(repr(hdu.header))
            ii = str(i)
            self.log_filter.hdu_for_log(i)

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
        aug_naxis = self.get_axis(index)

        # Chunk.time.axis.bounds.samples = time.samples

        # Chunk.time.axis.error.syser = CSYER{timeAxis}
        # Chunk.time.axis.error.rnder = CRDER{timeAxis}
        aug_error = self.get_coord_error(None, index)

        # Chunk.time.axis.function.naxis = NAXIS{timeAxis}
        # Chunk.time.axis.function.delta = CDELT{timeAxis}
        # Chunk.time.axis.function.refCoord.pix = CRPIX{timeAxis}
        # Chunk.time.axis.function.refCoord.val = CRVAL{timeAxis}
        aug_ref_coord = self.get_ref_coord(None, index)
        aug_function = CoordFunction1D(self.wcs._naxis[index], self.wcs.wcs.cdelt[index], aug_ref_coord)

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

        parser.add_argument('-o', '--out', dest='out_obs_xml', help='output of augmented observation in XML',
                            required=False)
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

    def _configure_logging(self):
        self.logger = logging.getLogger(__name__)
        if not len(self.logger.handlers):
            handler = logging.StreamHandler()
            formatter = logging.Formatter(
                '%(asctime)s %(name)-12s %(levelname)-8s HDU:%(hdu)-2d %(message)s')
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)

        self.logger.setLevel(logging.DEBUG)
        self.log_filter = LoggingFilter()
        self.logger.addFilter(self.log_filter)
        logging.getLogger('astropy').addFilter(self.log_filter)

    def augment_observation(self, aug_obs, artifact_uri):
        self.logger.debug('Begin CAOM2 observation augmentation for {} based on URI {}.'.format(self.file, artifact_uri))
        if aug_obs:
            raise NotImplementedError
        else:
            aug_collection = self.get_from_list('Observation.collection', 0, 'UNKNOWN')  # TODO default value
            aug_observation_id = self.get_from_list('Observation.observation_id', 0, 'UNKNOWN')  # TODO default value
            aug_algorithm = self.get_algorithm()
            aug_obs = Observation(aug_collection, aug_observation_id, aug_algorithm)

            # TODO default values for the following fields
            aug_obs.sequence_number = self.get_from_list('Observation.sequence_number', 0, -1)
            aug_obs.intent = self.get_from_list('Observation.intent', 0, ObservationIntentType.SCIENCE)
            aug_obs.type = self.get_from_list('Observation.type', 0, 'UNKNOWN')
            aug_obs.meta_release = self.get_datetime(self.get_from_list('Observation.meta_release', 0, datetime.now()))
            aug_obs.requirements = self.get_requirements()
            aug_obs.instrument = self.get_instrument()
            aug_obs.proposal = self.get_proposal()
            aug_obs.target = self.get_target()
            aug_obs.target_position = self.get_target_position()
            aug_obs.telescope = self.get_telescope()
            aug_obs.environment = self.get_environment()

            if aug_obs.planes:
                raise NotImplementedError
            else:
                prod_id = self.get_from_list('Plane.product_id', index=0, default='None')
                aug_obs.planes[prod_id] = Plane(prod_id)
                aug_plane = aug_obs.planes[prod_id]

            self.augment_plane(aug_plane, artifact_uri)

        self.logger.debug('End CAOM2 observation augmentation for {}.'.format(self.file))
        return aug_obs

    def augment_plane(self, aug_plane, artifact_uri):
        self.logger.debug('Begin CAOM2 plane augmentation.')
        if aug_plane:
            aug_plane.creator = self.get_from_list('Plane.creator_id', index=0, default='UNKNOWN')
            aug_plane.meta_release = self.get_from_list('Plane.meta_release', index=0, default=None)
            aug_plane.data_release = self.get_from_list('Plane.data_release', index=0, default=None)
            aug_plane.data_product_type = self.get_from_list('Plane.data_product_type', index=0,
                                                             default=DataProductType.IMAGE)
            aug_plane.calibration_level = self.get_from_list('Plane.calibration_level', index=0,
                                                             default=CalibrationLevel.PLANNED)
            aug_plane.provenance = self.get_provenance()
            aug_plane.metrics = self.get_metrics()
            aug_plane.quality = self.get_quality()

            for ii in aug_plane.artifacts:
                aug_artifact = aug_plane[ii]
                if aug_artifact.uri == artifact_uri:
                    self.artifact = aug_artifact
                    break

            self.augment_artifact(artifact_uri)
            aug_plane.artifacts[artifact_uri] = self.artifact

        else:
            raise NotImplementedError

        self.logger.debug('End CAOM2 plane augmentation.')

    def get_algorithm(self):
        self.logger.debug('Begin CAOM2 Algorithm augmentation.')
        name = self.get_from_list('Observation.algorithm.name', index=0, default='DEFAULT')  # TODO DEFAULT VALUE
        self.logger.debug('End CAOM2 Algorithm augmentation.')
        if name:
            return Algorithm(name)
        else:
            return None

    def get_instrument(self):
        self.logger.debug('Begin CAOM2 Instrument augmentation.')
        name = self.get_from_list('Observation.instrument.name', index=0, default='UNKNOWN')  # TODO DEFAULT VALUE
        keywords = self.get_from_list('Observation.instrument.keywords', index=0, default=['UNKNOWN'])  # TODO
        self.logger.debug('End CAOM2 Instrument augmentation.')
        if name:
            instr = Instrument(str(name))
            instr.keywords.union(keywords)
            return instr
        else:
            return None

    def get_proposal(self):
        self.logger.debug('Begin CAOM2 Proposal augmentation.')
        id = self.get_from_list('Observation.proposal.id', index=0, default='UNKNOWN')  # TODO
        pi = self.get_from_list('Observation.proposal.pi_name', index=0, default='UNKNOWN')  # TODO
        project = self.get_from_list('Observation.proposal.project', index=0, default='UNKNOWN')  # TODO
        title = self.get_from_list('Observation.proposal.title', index=0, default='UNKNOWN')  # TODO
        self.logger.debug('End CAOM2 Proposal augmentation.')
        if id:
            return Proposal(id, pi, project, title)
        else:
            return None

    def get_target(self):
        self.logger.debug('Begin CAOM2 Target augmentation.')
        name = self.get_from_list('Observation.target.name', index=0, default='UNKNOWN')  # TODO
        target_type = self.get_from_list('Observation.target.target_type', index=0, default=TargetType.OBJECT)  # TODO
        standard = self.get_from_list('Observation.target.standard', index=0, default=False)  # TODO
        redshift = self.get_from_list('Observation.target.redshift', index=0, default=-0.5)  # TODO
        keywords = self.get_from_list('Observation.target.keywords', index=0, default=set('UNKNOWN'))  # TODO
        moving = self.get_from_list('Observation.target.moving', index=0, default=False)  # TODO
        self.logger.debug('End CAOM2 Target augmentation.')
        if name:
            return Target(str(name), target_type, standard, redshift, keywords, moving)
        else:
            return None

    def get_target_position(self):
        self.logger.debug('Begin CAOM2 TargetPosition augmentation.')
        # TODO don't know what to do here, since config file says this is Chunk-level metadata
        self.logger.debug('End CAOM2 TargetPosition augmentation.')
        return None

    def get_telescope(self):
        self.logger.debug('Begin CAOM2 Telescope augmentation.')
        name = self.get_from_list('Observation.telescope.name', index=0, default='UNKNOWN')  # TODO
        geo_x = self.get_from_list('Observation.telescope.geo_location_x', index=0, default=-1.0)  # TODO
        geo_y = self.get_from_list('Observation.telescope.geo_location_y', index=0, default=-1.0)  # TODO
        geo_z = self.get_from_list('Observation.telescope.geo_location_z', index=0, default=-1.0)  # TODO
        keywords = self.get_from_list('Observation.telescope.keywords', index=0, default=set('UNKNOWN'))  # TODO
        self.logger.debug('End CAOM2 Telescope augmentation.')
        if name:
            return Telescope(name, geo_x, geo_y, geo_z, keywords)
        else:
            return None

    def get_environment(self):
        self.logger.debug('Begin CAOM2 Environment augmentation.')
        seeing = self.get_from_list('Observation.environment.seeing', index=0, default=None)  # TODO
        humidity = self.get_from_list('Observation.environment.humidity', index=0, default=None)  # TODO
        elevation = self.get_from_list('Observation.environment.elevation', index=0, default=None)  # TODO
        tau = self.get_from_list('Observation.environment.tau', index=0, default=None)  # TODO
        wavelength_tau = self.get_from_list('Observation.environment.wavelengthTau', index=0, default=None)  # TODO
        ambient = self.get_from_list('Observation.environment.ambientTemp', index=0, default=None)  # TODO
        photometric = self.get_from_list('Observation.environment.photometric', index=0, default=None)  # TODO
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

    def get_requirements(self):
        self.logger.debug('Begin CAOM2 Requirement augmentation.')
        flag = self.get_from_list('Observation.requirements.flag', index=0, default=Status.FAIL)  # TODO DEFAULT VALUE
        self.logger.debug('End CAOM2 Requirement augmentation.')
        if flag:
            return Requirements(flag)
        else:
            return None

    def get_from_list(self, lookup, index, default=None):
        value = default
        try:
            keywords = CONFIG[lookup]
        except KeyError:
            self.logger.debug('Could not find lookup value \'{}\' in fits2caom2 configuration.'.format(lookup))
            return value

        for ii in keywords:
            value = self.hdulist[index].header.get(ii, default)
            self.logger.debug('Assigned value {} based on keyword {}'.format(value, ii))
            if value is not default:
                break

        return value

    def get_provenance(self):
        self.logger.debug('Begin CAOM2 Provenance augmentation.')
        name = self.get_from_list('Plane.provenance.name', index=0, default='DEFAULT')  # TODO DEFAULT VALUE
        version = self.get_from_list('Plane.provenance.version', index=0, default='DEFAULT')  # TODO DEFAULT VALUE
        project = self.get_from_list('Plane.provenance.project', index=0, default='DEFAULT')  # TODO DEFAULT VALUE
        producer = self.get_from_list('Plane.provenance.producer', index=0, default='DEFAULT')  # TODO DEFAULT VALUE
        run_id = self.get_from_list('Plane.provenance.runID', index=0, default='DEFAULT')  # TODO DEFAULT VALUE
        reference = self.get_from_list('Plane.provenance.reference', index=0, default='DEFAULT')  # TODO DEFAULT VALUE
        last_executed = self.get_datetime(self.get_from_list('Plane.provenance.lastExecuted', index=0))  # TODO DEFAULT VALUE
        keywords = self.get_from_list('Plane.provenance.keywords', index=0, default='DEFAULT')  # TODO DEFAULT VALUE
        #inputs = self.get_from_list('Plane.provenance.inputs', index=0, default=set(PlaneURI('caom:UNKNOWN/UNKNOWN/UNKNOWN')))  # TODO DEFAULT VALUE
        inputs = self.get_from_list('Plane.provenance.inputs', index=0, default=None)  # TODO DEFAULT VALUE
        self.logger.debug('End CAOM2 Provenance augmentation.')
        if name:
            prov = Provenance(name, version, project, producer, run_id, reference, last_executed)
            prov.keywords.union(keywords)
            if inputs:
                prov.inputs.add(inputs)
            return prov
        else:
            return None

    def get_metrics(self):
        self.logger.debug('Begin CAOM2 Metrics augmentation.')
        source_number_density = self.get_from_list('Plane.metrics.sourceNumberDensity', index=0, default=1.0)  # TODO DEFAULT VALUE
        background = self.get_from_list('Plane.metrics.background', index=0, default=1.0)  # TODO DEFAULT VALUE
        background_stddev = self.get_from_list('Plane.metrics.backgroundStddev', index=0, default=1.0)  # TODO DEFAULT VALUE
        flux_density_limit = self.get_from_list('Plane.metrics.fluxDensityLimit', index=0, default=1.0)  # TODO DEFAULT VALUE
        mag_limit = self.get_from_list('Plane.metrics.magLimit', index=0, default=1.0)  # TODO DEFAULT VALUE
        metrics = Metrics()
        metrics.source_number_density = source_number_density
        metrics.background = background
        metrics.background_std_dev = background_stddev
        metrics.flux_density_limit = flux_density_limit
        metrics.mag_limit = mag_limit
        self.logger.debug('End CAOM2 Metrics augmentation.')
        return metrics

    def get_quality(self):
        self.logger.debug('Begin CAOM2 Quality augmentation.')
        flag = self.get_from_list('Plane.dataQuality', index=0, default=Quality.JUNK)  # TODO DEFAULT VALUE
        self.logger.debug('End CAOM2 Quality augmentation.')
        if flag:
            return DataQuality(flag)
        else:
            return None

    def get_datetime(self, from_value):
        return datetime(1990, 1, 1, 12, 12, 12)


CONFIG = {'Observation.meta_release': ['DATE', 'DATE-OBS', 'UTCOBS', 'UTCDATE', 'UTC-DATE', 'MJDOBS', 'MJD_OBS'],
          'Observation.instrument.name': ['INSTRUME'],
          'Observation.target.name': ['OBJECT'],
          'Observation.type': ['OBSTYPE'],
          'Observation.telescope.name': ['TELESCOP'],
          'Observation.environment.ambientTemp': ['TEMPERAT'],
          'Observation.algorithm.name': [],  # TODO
          'Observation.intent': [],
          'Observation.sequenceNumber': [],
          'Observation.instrument.keywords': ['INSTMODE'],
          'Observation.proposal.id': ['RUNID'],
          'Observation.proposal.pi': [],
          'Observation.proposal.project': [],
          'Observation.proposal.title': [],
          'Observation.proposal.keywords': [],
          'Observation.target.type': [],
          'Observation.target.standard': [],
          'Observation.target.redshift': [],
          'Observation.target.keywords': [],
          'Observation.telescope.geoLocationX': [],
          'Observation.telescope.geoLocationY': [],
          'Observation.telescope.geoLocationZ': [],
          'Observation.telescope.keywords': [],
          'Observation.environment.seeing': [],
          'Observation.environment.humidity': [],
          'Observation.environment.elevation': [],
          'Observation.environment.tau': [],
          'Observation.environment.wavelengthTau': [],
          'Observation.environment.photometric': [],
          'Plane.meta_release': ['RELEASE', 'REL_DATE'],
          'Plane.data_release': ['RELEASE', 'REL_DATE'],
          'Plane.dataProductType': [],
          'Plane.calibrationLevel': [],
          'Plane.provenance.name': [],
          'Plane.provenance.version': [],
          'Plane.provenance.project': [],
          'Plane.provenance.producer': [],
          'Plane.provenance.runID': [],
          'Plane.provenance.reference': [],
          'Plane.provenance.lastExecuted': [],
          'Plane.provenance.keywords': [],
          'Plane.provenance.inputs': [],
          'Plane.metrics.sourceNumberDensity': [],
          'Plane.metrics.background': [],
          'Plane.metrics.backgroundStddev': [],
          'Plane.metrics.fluxDensityLimit': [],
          'Plane.metrics.magLimit': []
          }
#          '': [],
