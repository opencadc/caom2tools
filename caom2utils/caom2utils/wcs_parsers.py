# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2024.                            (c) 2024.
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
#  Revision: 4
#
# ***********************************************************************
#

import logging
import math
import sys

from astropy.wcs import SingularMatrixError, utils, Wcsprm, WCS
from caom2 import (
    Axis,
    Chunk,
    Coord2D,
    CoordAxis1D,
    CoordAxis2D,
    CoordError,
    CoordFunction1D,
    CoordFunction2D,
    CustomWCS,
    Dimension2D,
    ObservableAxis,
    PolarizationWCS,
    RefCoord,
    Slice,
    SpatialWCS,
    SpectralWCS,
    TemporalWCS,
)
from caom2utils.blueprints import ObsBlueprint, _to_float, _to_int, _to_str


CUSTOM_CTYPES = ['RM', 'FDEP']

POSITION_CTYPES = [['RA', 'GLON', 'ELON', 'HLON', 'SLON'], ['DEC', 'GLAT', 'ELAT', 'HLAT', 'SLAT']]

ENERGY_CTYPES = ['FREQ', 'ENER', 'WAVN', 'VRAD', 'WAVE', 'VOPT', 'ZOPT', 'AWAV', 'VELO', 'BETA']

# From http://hea-www.cfa.harvard.edu/~arots/TimeWCS/
TIME_KEYWORDS = ['TIME', 'TAI', 'TT', 'TDT', 'ET', 'IAT', 'UT1', 'UTC', 'GMT', 'GPS', 'TCG', 'TCB', 'TDB', 'LOCAL']

POLARIZATION_CTYPES = ['STOKES']

OBSERVABLE_CTYPES = ['observable', 'FLUX']


class HDULoggingFilter(logging.Filter):
    """Add the HDU number to logging messages as a default."""

    def __init__(self):
        super().__init__()
        self._extension = -1

    def filter(self, record):
        record.hdu = self._extension
        return True

    def extension(self, value):
        self._extension = value


class WcsParser:
    """
    WCS axes methods.
    """

    ENERGY_AXIS = 'energy'
    POLARIZATION_AXIS = 'polarization'
    TIME_AXIS = 'time'

    def __init__(self, blueprint, extension):
        self._wcs = None
        self.wcs = None
        self._blueprint = blueprint
        self._axes = {
            'ra': [0, False],
            'dec': [0, False],
            'time': [0, False],
            'energy': [0, False],
            'polarization': [0, False],
            'observable': [0, False],
            'custom': [0, False],
        }
        # int - index into blueprint._plan extensions
        self._extension = extension
        self.logger = logging.getLogger(self.__class__.__name__)
        self._set_wcs()

    def _assign_cd(self, key, cd, count):
        x = self._blueprint._get(key, self._extension)
        if x is not None:
            if ObsBlueprint.needs_lookup(x):
                cd[count][count] = 1.0
            else:
                cd[count][count] = x

    def assign_sanitize(self, assignee, index, key, sanitize=True):
        """
        Do not want to blindly assign None to astropy.wcs attributes, so use this method for conditional assignment.

        The current implementation is that if there is a legitimate need to assign None to a value, either use 'set'
        in the Hdf5ObsBlueprint, and specifically assign None, or execute a function to set it to None conditionally.
        There will be no support for a Default value of None with HDF5 files.

        By the time this method is called, if the value still passes the "ObsBlueprint.needs_lookup" check, the
        value should be ignored for fulfilling the WCS needs of the record under construction.
        """
        x = self._blueprint._get(key, self._extension)
        if sanitize:
            x = self._sanitize(x)
        if x is not None and not ObsBlueprint.needs_lookup(x):
            assignee[index] = x

    def _set_wcs(self):
        num_axes = self._blueprint.get_configed_axes_count()
        self._wcs = WCS(naxis=num_axes)
        self.wcs = self._wcs.wcs
        array_shape, crder, crpix, crval, csyer, ctype, cunit, temp = [[0] * num_axes for _ in range(8)]
        cd = [temp.copy() for _ in range(num_axes)]
        count = 0
        if self._blueprint._pos_axes_configed:
            self._axes['ra'][1] = True
            self._axes['dec'][1] = True
            self._axes['ra'][0] = count
            self._axes['dec'][0] = count + 1
            self.assign_sanitize(ctype, count, 'Chunk.position.axis.axis1.ctype')
            self.assign_sanitize(ctype, count + 1, 'Chunk.position.axis.axis2.ctype')
            self.assign_sanitize(cunit, count, 'Chunk.position.axis.axis1.cunit')
            self.assign_sanitize(cunit, count + 1, 'Chunk.position.axis.axis2.cunit')
            self.assign_sanitize(array_shape, count, 'Chunk.position.axis.function.dimension.naxis1')
            self.assign_sanitize(array_shape, count + 1, 'Chunk.position.axis.function.dimension.naxis2')
            self.assign_sanitize(crpix, count, 'Chunk.position.axis.function.refCoord.coord1.pix')
            self.assign_sanitize(crpix, count + 1, 'Chunk.position.axis.function.refCoord.coord2.pix')
            self.assign_sanitize(crval, count, 'Chunk.position.axis.function.refCoord.coord1.val')
            self.assign_sanitize(crval, count + 1, 'Chunk.position.axis.function.refCoord.coord2.val')
            x = self._blueprint._get('Chunk.position.axis.function.cd11', self._extension)
            if x is not None and not ObsBlueprint.needs_lookup(x):
                cd[count][0] = x
            x = self._blueprint._get('Chunk.position.axis.function.cd12', self._extension)
            if x is not None and not ObsBlueprint.needs_lookup(x):
                cd[count][1] = x
            x = self._blueprint._get('Chunk.position.axis.function.cd21', self._extension)
            if x is not None and not ObsBlueprint.needs_lookup(x):
                cd[count + 1][0] = x
            x = self._blueprint._get('Chunk.position.axis.function.cd22', self._extension)
            if x is not None and not ObsBlueprint.needs_lookup(x):
                cd[count + 1][1] = x
            self.assign_sanitize(crder, count, 'Chunk.position.axis.error1.rnder')
            self.assign_sanitize(crder, count + 1, 'Chunk.position.axis.error2.rnder')
            self.assign_sanitize(csyer, count, 'Chunk.position.axis.error1.syser')
            self.assign_sanitize(csyer, count + 1, 'Chunk.position.axis.error2.syser')
            count += 2
        if self._blueprint._time_axis_configed:
            self._axes['time'][1] = True
            self._axes['time'][0] = count
            self.assign_sanitize(ctype, count, 'Chunk.time.axis.axis.ctype', False)
            self.assign_sanitize(cunit, count, 'Chunk.time.axis.axis.cunit', False)
            self.assign_sanitize(array_shape, count, 'Chunk.time.axis.function.naxis', False)
            self.assign_sanitize(crpix, count, 'Chunk.time.axis.function.refCoord.pix', False)
            self.assign_sanitize(crval, count, 'Chunk.time.axis.function.refCoord.val', False)
            self.assign_sanitize(crder, count, 'Chunk.time.axis.error.rnder')
            self.assign_sanitize(csyer, count, 'Chunk.time.axis.error.syser')
            self._assign_cd('Chunk.time.axis.function.delta', cd, count)
            count += 1
        if self._blueprint._energy_axis_configed:
            self._axes['energy'][1] = True
            self._axes['energy'][0] = count
            self.assign_sanitize(ctype, count, 'Chunk.energy.axis.axis.ctype', False)
            self.assign_sanitize(cunit, count, 'Chunk.energy.axis.axis.cunit', False)
            self.assign_sanitize(array_shape, count, 'Chunk.energy.axis.function.naxis', False)
            self.assign_sanitize(crpix, count, 'Chunk.energy.axis.function.refCoord.pix', False)
            self.assign_sanitize(crval, count, 'Chunk.energy.axis.function.refCoord.val', False)
            self.assign_sanitize(crder, count, 'Chunk.energy.axis.error.rnder')
            self.assign_sanitize(csyer, count, 'Chunk.energy.axis.error.syser')
            self._assign_cd('Chunk.energy.axis.function.delta', cd, count)
            count += 1
        if self._blueprint._polarization_axis_configed:
            self._axes['polarization'][1] = True
            self._axes['polarization'][0] = count
            self.assign_sanitize(ctype, count, 'Chunk.polarization.axis.axis.ctype', False)
            self.assign_sanitize(cunit, count, 'Chunk.polarization.axis.axis.cunit', False)
            self.assign_sanitize(array_shape, count, 'Chunk.polarization.axis.function.naxis', False)
            self.assign_sanitize(crpix, count, 'Chunk.polarization.axis.function.refCoord.pix', False)
            self.assign_sanitize(crval, count, 'Chunk.polarization.axis.function.refCoord.val', False)
            self._assign_cd('Chunk.polarization.axis.function.delta', cd, count)
            count += 1
        if self._blueprint._obs_axis_configed:
            self._axes['observable'][1] = True
            self._axes['observable'][0] = count
            self.assign_sanitize(ctype, count, 'Chunk.observable.axis.axis.ctype', False)
            self.assign_sanitize(cunit, count, 'Chunk.observable.axis.axis.cunit', False)
            array_shape[count] = 1.0
            self.assign_sanitize(crpix, count, 'Chunk.observable.axis.function.refCoord.pix', False)
            crval[count] = 0.0
            cd[count][count] = 1.0
            count += 1
        if self._blueprint._custom_axis_configed:
            self._axes['custom'][1] = True
            self._axes['custom'][0] = count
            self.assign_sanitize(ctype, count, 'Chunk.custom.axis.axis.ctype', False)
            self.assign_sanitize(cunit, count, 'Chunk.custom.axis.axis.cunit', False)
            self.assign_sanitize(array_shape, count, 'Chunk.custom.axis.function.naxis', False)
            self.assign_sanitize(crpix, count, 'Chunk.custom.axis.function.refCoord.pix', False)
            self.assign_sanitize(crval, count, 'Chunk.custom.axis.function.refCoord.val', False)
            self._assign_cd('Chunk.custom.axis.function.delta', cd, count)
            count += 1

        if not all(val == 0 for val in array_shape):
            self._wcs.array_shape = array_shape
        if not all(val == 0 for val in cunit):
            self._wcs.wcs.cunit = cunit
        if not all(val == 0 for val in ctype):
            self._wcs.wcs.ctype = ctype
        if not all(val == 0 for val in crpix):
            self._wcs.wcs.crpix = crpix
        if not all(val == 0 for val in crval):
            self._wcs.wcs.crval = crval
        if not all(val == 0 for val in crder):
            self._wcs.wcs.crder = crder
        if not all(val == 0 for val in csyer):
            self._wcs.wcs.csyer = csyer
        self._wcs.wcs.cd = cd
        self._finish_position()
        self._finish_time()
        self._finish_energy()

    def augment_custom(self, chunk):
        """
        Augments a chunk with custom WCS information
        :param chunk:
        :return:
        """
        self.logger.debug('Begin Custom WCS augmentation.')
        if chunk is None or not isinstance(chunk, Chunk):
            raise ValueError(f'Chunk type mis-match for {chunk}.')

        custom_axis_index = self._get_axis_index(CUSTOM_CTYPES)
        if custom_axis_index is None:
            self.logger.debug('No WCS Custom info')
            return
        try:
            custom_axis_length = self._get_axis_length(custom_axis_index + 1)
        except ValueError:
            self.logger.debug('No WCS Custom axis.function')
            return

        if custom_axis_length:
            chunk.custom_axis = custom_axis_index + 1
            naxis = CoordAxis1D(self._get_axis(custom_axis_index))
            if self.wcs.has_cd():
                delta = self.wcs.cd[custom_axis_index][custom_axis_index]
            else:
                delta = self.wcs.cdelt[custom_axis_index]
            ref_coord = self._get_ref_coord(custom_axis_index)
            if delta and ref_coord:
                naxis.function = CoordFunction1D(custom_axis_length, delta, ref_coord)
            chunk.custom = CustomWCS(naxis)

        self.logger.debug('End Custom WCS augmentation.')

    def augment_energy(self, chunk):
        """
        Augments the energy information in a chunk
        :param chunk:
        """
        self.logger.debug('Begin Energy WCS augmentation.')
        if chunk is None or not isinstance(chunk, Chunk):
            raise ValueError(f'Chunk type mis-match for {chunk}.')

        # get the energy axis
        energy_axis_index = self._get_axis_index(ENERGY_CTYPES)

        if energy_axis_index is None:
            self.logger.debug('No WCS Energy info.')
            return
        try:
            energy_axis_length = self._get_axis_length(energy_axis_index + 1)
        except ValueError:
            self.logger.debug('No WCS Energy axis.function')
            return

        if energy_axis_length:
            chunk.energy_axis = energy_axis_index + 1
            naxis = CoordAxis1D(self._get_axis(energy_axis_index))
            naxis.error = self._get_coord_error(energy_axis_index)
            if self.wcs.has_cd():
                delta = self.wcs.cd[energy_axis_index][energy_axis_index]
            else:
                delta = self.wcs.cdelt[energy_axis_index]
            ref_coord = self._get_ref_coord(energy_axis_index)
            if delta and ref_coord:
                naxis.function = CoordFunction1D(energy_axis_length, delta, ref_coord)

            specsys = _to_str(self.wcs.specsys) if self.wcs.specsys else ''
            if not chunk.energy:
                chunk.energy = SpectralWCS(naxis, specsys)
            else:
                chunk.energy.axis = naxis
                chunk.energy.specsys = specsys

            chunk.energy.ssysobs = _to_str(self._sanitize(self.wcs.ssysobs))
            # wcs returns 0.0 by default
            if self._sanitize(self.wcs.restfrq) != 0:
                chunk.energy.restfrq = self._sanitize(self.wcs.restfrq)
            if self._sanitize(self.wcs.restwav) != 0:
                chunk.energy.restwav = self._sanitize(self.wcs.restwav)
            chunk.energy.velosys = self._sanitize(self.wcs.velosys)
            chunk.energy.zsource = self._sanitize(self.wcs.zsource)
            chunk.energy.ssyssrc = _to_str(self._sanitize(self.wcs.ssyssrc))
            chunk.energy.velang = self._sanitize(self.wcs.velangl)
        self.logger.debug('End Energy WCS augmentation.')

    def augment_position(self, chunk):
        """
        Augments a chunk with spatial WCS information
        :param chunk:
        :return:
        """
        self.logger.debug('Begin Spatial WCS augmentation.')
        if chunk is None or not isinstance(chunk, Chunk):
            raise ValueError(f'Chunk type mis-match for {chunk}.')

        position_axes_indices = self._get_position_axis()
        if not position_axes_indices:
            self.logger.debug('No Spatial WCS found')
            return

        chunk.position_axis_1 = position_axes_indices[0]
        chunk.position_axis_2 = position_axes_indices[1]
        axis = self._get_spatial_axis(chunk.position_axis_1 - 1, chunk.position_axis_2 - 1)

        if axis is None:
            self.logger.debug('No WCS Position axis.function')
            return

        if chunk.position:
            chunk.position.axis = axis
        else:
            chunk.position = SpatialWCS(axis)

        chunk.position.coordsys = _to_str(self._sanitize(self.wcs.radesys))
        temp = self._sanitize(self.wcs.equinox)
        if (temp is not None and 1800.0 <= temp <= 2500) or temp is None:
            chunk.position.equinox = temp

        self._finish_chunk_position(chunk)
        self.logger.debug('End Spatial WCS augmentation.')

    def augment_temporal(self, chunk):
        """
        Augments a chunk with temporal WCS information

        :param chunk:
        :return:
        """
        self.logger.debug('Begin TemporalWCS augmentation.')
        if chunk is None or not isinstance(chunk, Chunk):
            raise ValueError(f'Chunk type mis-match for {chunk}.')

        time_axis_index = self._get_axis_index(TIME_KEYWORDS)

        if time_axis_index is None:
            self.logger.debug('No WCS Time info.')
            return

        chunk.time_axis = time_axis_index + 1
        # set chunk.time
        self.logger.debug('Begin temporal axis augmentation.')

        try:
            axis_length = self._get_axis_length(time_axis_index + 1)
        except ValueError:
            self.logger.debug('No WCS Temporal axis.function')
            return

        if axis_length:
            aug_naxis = self._get_axis(time_axis_index)
            aug_error = self._get_coord_error(time_axis_index)
            aug_ref_coord = self._get_ref_coord(time_axis_index)
            if self.wcs.has_cd():
                delta = self.wcs.cd[time_axis_index][time_axis_index]
            else:
                delta = self.wcs.cdelt[time_axis_index]
            if aug_ref_coord is not None:
                aug_function = CoordFunction1D(axis_length, delta, aug_ref_coord)
                naxis = CoordAxis1D(aug_naxis, aug_error, None, None, aug_function)
                if not chunk.time:
                    chunk.time = TemporalWCS(naxis)
                else:
                    chunk.time.axis = naxis

                self._finish_chunk_time(chunk)
        self.logger.debug('End TemporalWCS augmentation.')

    def augment_polarization(self, chunk):
        """
        Augments a chunk with polarization WCS information
        :param chunk:
        :return:
        """
        self.logger.debug('Begin Polarization WCS augmentation.')
        if chunk is None or not isinstance(chunk, Chunk):
            raise ValueError(f'Chunk type mis-match for {chunk}.')

        polarization_axis_index = self._get_axis_index(POLARIZATION_CTYPES)
        if polarization_axis_index is None:
            self.logger.debug('No WCS Polarization info')
            return

        try:
            axis_length = self._get_axis_length(polarization_axis_index + 1)
        except ValueError:
            self.logger.debug('No WCS Polarization axis.function')
            return

        if axis_length:
            chunk.polarization_axis = polarization_axis_index + 1

            naxis = CoordAxis1D(self._get_axis(polarization_axis_index))
            if self.wcs.has_cd():
                delta = self.wcs.cd[polarization_axis_index][polarization_axis_index]
            else:
                delta = self.wcs.cdelt[polarization_axis_index]
            ref_coord = self._get_ref_coord(polarization_axis_index)
            if delta and ref_coord:
                naxis.function = CoordFunction1D(axis_length, delta, ref_coord)
            if not chunk.polarization:
                chunk.polarization = PolarizationWCS(naxis)
            else:
                chunk.polarization.axis = naxis

        self.logger.debug('End Polarization WCS augmentation.')

    def augment_observable(self, chunk):
        """
        Augments a chunk with an observable axis.

        :param chunk:
        :return:
        """
        self.logger.debug('Begin Observable WCS augmentation.')
        if chunk is None or not isinstance(chunk, Chunk):
            raise ValueError(f'Chunk type mis-match for {chunk}.')

        observable_axis_index = self._get_axis_index(OBSERVABLE_CTYPES)
        if observable_axis_index is None:
            self.logger.debug('No Observable axis info')
            return

        chunk.observable_axis = observable_axis_index + 1
        self._finish_chunk_observable(chunk)
        self.logger.debug('End Observable WCS augmentation.')

    def _finish_chunk_observable(self, chunk):
        self.logger.debug('Begin _finish_chunk_observable')
        ctype = self._wcs.wcs.ctype[chunk.observable_axis - 1]
        cunit = self._wcs.wcs.ctype[chunk.observable_axis - 1]
        pix_bin = _to_int(self._wcs.wcs.crpix[chunk.observable_axis - 1])
        if ctype is not None and cunit is not None and pix_bin is not None:
            chunk.observable = ObservableAxis(Slice(self._get_axis(0, ctype, cunit), pix_bin))
        self.logger.debug('End _finish_chunk_observable')

    def _finish_chunk_position(self, chunk):
        self.logger.debug('Begin _finish_chunk_position')
        if chunk.position.resolution is None:
            try:
                # JJK 30-01-23
                # In a spatial data chunk the resolution is 2 times the pixel size.  We can get the pixel size from
                # the wcs
                temp = utils.proj_plane_pixel_scales(self._wcs)
                chunk.position.resolution = temp[0]
            except SingularMatrixError as e:
                # cannot calculate position.resolution, ignore and continue on
                self.logger.warning(f'Not calculating resolution due to {e}')
        self.logger.debug('End _finish_chunk_position')

    def _finish_chunk_time(self, chunk):
        self.logger.debug('Begin _finish_chunk_time')
        if not math.isnan(self._wcs.wcs.xposure):
            chunk.time.exposure = self._wcs.wcs.xposure
        if self._wcs.wcs.timesys is not None and self._wcs.wcs.timesys != '':
            chunk.time.timesys = self._wcs.wcs.timesys
        if self._wcs.wcs.trefpos is not None and self._wcs.wcs.trefpos != '':
            chunk.time.trefpos = self._wcs.wcs.trefpos
        if self._wcs.wcs.mjdref is not None and self._wcs.wcs.mjdref[0] != '' and self._wcs.wcs.mjdref[0] != 0.0:
            # the astropy value is an array of length 2, use the first value
            chunk.time.mjdref = self._wcs.wcs.mjdref[0]
        self.logger.debug('End _finish_chunk_time')

    def _finish_energy(self):
        self.logger.debug('Begin _finish_energy')
        if self._blueprint._energy_axis_configed:
            x = self._blueprint._get('Chunk.energy.specsys', self._extension)
            if x and not ObsBlueprint.needs_lookup(x):
                self._wcs.wcs.specsys = x
            x = self._blueprint._get('Chunk.energy.ssysobs', self._extension)
            if x and not ObsBlueprint.needs_lookup(x):
                self._wcs.wcs.ssysobs = x
            x = self._blueprint._get('Chunk.energy.restfrq', self._extension)
            if x and not ObsBlueprint.needs_lookup(x):
                self._wcs.wcs.restfrq = _to_float(x)
            x = self._blueprint._get('Chunk.energy.restwav', self._extension)
            if x and not ObsBlueprint.needs_lookup(x):
                self._wcs.wcs.restwav = x
            x = self._blueprint._get('Chunk.energy.velosys', self._extension)
            if x and not ObsBlueprint.needs_lookup(x):
                self._wcs.wcs.velosys = x
            x = self._blueprint._get('Chunk.energy.zsource', self._extension)
            if x and not ObsBlueprint.needs_lookup(x):
                self._wcs.wcs.zsource = x
            x = self._blueprint._get('Chunk.energy.ssyssrc', self._extension)
            if x and not ObsBlueprint.needs_lookup(x):
                self._wcs.wcs.ssyssrc = x
            x = self._blueprint._get('Chunk.energy.velang', self._extension)
            if x and not ObsBlueprint.needs_lookup(x):
                self._wcs.wcs.velangl = x
        self.logger.debug('End _finish_energy')

    def _finish_position(self):
        self.logger.debug('Begin _finish_position')
        if self._blueprint._pos_axes_configed:
            x = self._blueprint._get('Chunk.position.coordsys', self._extension)
            if x and not ObsBlueprint.needs_lookup(x):
                self._wcs.wcs.radesys = x
            x = self._blueprint._get('Chunk.position.equinox', self._extension)
            if x and not ObsBlueprint.needs_lookup(x):
                self._wcs.wcs.equinox = _to_float(x)
        self.logger.debug('End _finish_position')

    def _finish_time(self):
        self.logger.debug('Begin _finish_time')
        if self._blueprint._time_axis_configed:
            x = self._blueprint._get('Chunk.time.exposure', self._extension)
            if x and not ObsBlueprint.needs_lookup(x):
                self._wcs.wcs.xposure = _to_float(x)
            x = self._blueprint._get('Chunk.time.timesys', self._extension)
            if x and not ObsBlueprint.needs_lookup(x):
                self._wcs.wcs.timesys = x
            x = self._blueprint._get('Chunk.time.trefpos', self._extension)
            if x and not ObsBlueprint.needs_lookup(x):
                self._wcs.wcs.trefpos = x
            x = self._blueprint._get('Chunk.time.mjdref', self._extension)
            if x and not ObsBlueprint.needs_lookup(x):
                self._wcs.wcs.mjdref = [x, x]
        self.logger.debug('End _finish_time')

    def _get_axis(self, index, over_ctype=None, over_cunit=None):
        """Assemble a generic axis"""
        aug_ctype = str(self.wcs.ctype[index]) if over_ctype is None else over_ctype
        aug_cunit = str(self.wcs.cunit[index]) if over_cunit is None else over_cunit
        if aug_cunit is not None and len(aug_cunit) == 0:
            aug_cunit = None
        aug_axis = Axis(aug_ctype, aug_cunit)
        return aug_axis

    def _get_axis_index(self, keywords):
        """
        Return the index of a specific axis type or None of it doesn't exist
        :param keywords:
        :return:
        """
        axis = None
        for i, elem in enumerate(self.wcs.ctype):
            elem = elem.split('-')[0]
            if elem in keywords:
                axis = i
                break
            elif len(elem) == 0:
                check = self.wcs.ctype[i]
                if check in keywords:
                    axis = i
                    break
        return axis

    def _get_axis_length(self, for_axis):
        if self._wcs.array_shape is None:
            return 0
        else:
            if len(self._wcs.array_shape) == 1:
                result = self._wcs.array_shape[0]
            else:
                result = self._wcs.array_shape[for_axis - 1]
            if isinstance(result, tuple):
                raise ValueError(
                    f'Could not find axis length for axis {for_axis}. The blueprint is incompletely configured.'
                )
            return _to_int(result)

    def _get_cd(self, x_index, y_index):
        """returns cd info"""

        try:
            if self.wcs.has_cd():
                cd11 = self.wcs.cd[x_index][x_index]
                cd12 = self.wcs.cd[x_index][y_index]
                cd21 = self.wcs.cd[y_index][x_index]
                cd22 = self.wcs.cd[y_index][y_index]
            else:
                cd11 = self.wcs.cdelt[x_index]
                cd12 = self.wcs.crota[x_index]
                cd21 = self.wcs.crota[y_index]
                cd22 = self.wcs.cdelt[y_index]
        except AttributeError:
            self.logger.debug(f'Error searching for CD* values {sys.exc_info()[1]}')
            cd11 = None
            cd12 = None
            cd21 = None
            cd22 = None

        return cd11, cd12, cd21, cd22

    def _get_coord_error(self, wcs_index):
        aug_coord_error = None
        aug_csyer = self._sanitize(self.wcs.csyer[wcs_index])
        aug_crder = self._sanitize(self.wcs.crder[wcs_index])
        if aug_csyer is not None and aug_crder is not None:
            aug_coord_error = CoordError(aug_csyer, aug_crder)
        return aug_coord_error

    def _get_dimension(self, xindex, yindex):
        aug_dimension = None
        try:
            xindex_axis_length = self._get_axis_length(xindex + 1)
            yindex_axis_length = self._get_axis_length(yindex + 1)
        except ValueError:
            self.logger.debug('No WCS Energy axis.function')
            return None

        if xindex_axis_length > 0 and yindex_axis_length > 0:
            aug_dim1 = _to_int(xindex_axis_length)
            aug_dim2 = _to_int(yindex_axis_length)
            if aug_dim1 and aug_dim2:
                aug_dimension = Dimension2D(aug_dim1, aug_dim2)
                self.logger.debug('End 2D dimension augmentation.')
        return aug_dimension

    def _get_position_axis(self):
        # there are two celestial axes, get the applicable indices from the axis_types
        xindex = self._get_axis_index(POSITION_CTYPES[0])
        yindex = self._get_axis_index(POSITION_CTYPES[1])

        if (xindex is not None) and (yindex is not None):
            return xindex + 1, yindex + 1
        elif (xindex is None) and (yindex is None):
            return None
        else:
            raise ValueError('Found only one position axis ra/dec: {}/{} in ' '{}'.format(xindex, yindex, self.file))

    def _get_ref_coord(self, index):
        aug_crpix = _to_float(self._sanitize(self.wcs.crpix[index]))
        aug_crval = _to_float(self._sanitize(self.wcs.crval[index]))
        aug_ref_coord = None
        if aug_crpix is not None and aug_crval is not None:
            aug_ref_coord = RefCoord(aug_crpix, aug_crval)
        return aug_ref_coord

    def _get_spatial_axis(self, xindex, yindex):
        """Assemble the bits to make the axis parameter needed for SpatialWCS construction."""
        aug_dimension = self._get_dimension(xindex, yindex)
        if aug_dimension is None:
            return None

        x_ref_coord = self._get_ref_coord(xindex)
        y_ref_coord = self._get_ref_coord(yindex)
        aug_ref_coord = None
        if x_ref_coord and y_ref_coord:
            aug_ref_coord = Coord2D(x_ref_coord, y_ref_coord)

        aug_cd11, aug_cd12, aug_cd21, aug_cd22 = self._get_cd(xindex, yindex)

        if (
            aug_dimension is not None
            and aug_ref_coord is not None
            and aug_cd11 is not None
            and aug_cd12 is not None
            and aug_cd21 is not None
            and aug_cd22 is not None
        ):
            aug_function = CoordFunction2D(aug_dimension, aug_ref_coord, aug_cd11, aug_cd12, aug_cd21, aug_cd22)
            self.logger.debug('End CoordFunction2D augmentation.')
        else:
            aug_function = None

        aug_axis = CoordAxis2D(
            self._get_axis(xindex),
            self._get_axis(yindex),
            self._get_coord_error(xindex),
            self._get_coord_error(yindex),
            None,
            None,
            aug_function,
        )
        self.logger.debug('End CoordAxis2D augmentation.')
        return aug_axis

    def _sanitize(self, value):
        """
        Sanitizes values from content to caom2
        :param value:
        :return:
        """
        if value is None:
            return None
        elif isinstance(value, float) and math.isnan(value):
            return None
        elif not str(value):
            return None  # empty string
        else:
            return value


class FitsWcsParser(WcsParser):
    """
    Parser to augment chunks with positional, temporal, energy and polarization information based on the WCS keywords
    in an extension of a FITS header.

    Note: Under the hood, this class uses the astropy.wcs package to parse the header and any inconsistencies or
    missing keywords are reported back as warnings.
    """

    def __init__(self, header, file, extension):
        """
        :param header: FITS extension header
        :param file: name of FITS file
        :param extension: which HDU
        """
        self.logger = logging.getLogger(self.__class__.__name__)
        self.log_filter = HDULoggingFilter()
        self.log_filter.extension(extension)
        self.logger.addFilter(self.log_filter)
        logastro = logging.getLogger('astropy')
        logastro.addFilter(self.log_filter)
        logastro.propagate = False
        header_string = header.tostring().rstrip()
        header_string = header_string.replace('END' + ' ' * 77, '')
        self.wcs = Wcsprm(header_string.encode('ascii'))
        self.wcs.fix()
        self.header = header
        self.file = file
        self.extension = extension

    def _finish_chunk_observable(self, chunk):
        self.logger.debug('Begin _finish_chunk_observable')
        ctype = self.header.get(f'CTYPE{chunk.observable_axis}')
        cunit = self.header.get(f'CUNIT{chunk.observable_axis}')
        pix_bin = self.header.get(f'CRPIX{chunk.observable_axis}')
        if ctype is not None and cunit is not None and pix_bin is not None:
            chunk.observable = ObservableAxis(Slice(self._get_axis(0, ctype, cunit), pix_bin))
        self.logger.debug('End _finish_chunk_observable')

    def _finish_chunk_position(self, chunk):
        pass

    def _finish_chunk_time(self, chunk):
        """
        The expected caom2 - FITS keywords mapping is:

        time.exposure = EXPTIME
        time.resolution = TIMEDEL
        time.timesys = TIMESYS default UTC
        time.trefpos = TREFPOS
        time.mjdref = MJDREF | MJDDATE
        """
        self.logger.debug('Begin _finish_chunk_time')
        chunk.time.exposure = _to_float(self.header.get('EXPTIME'))
        chunk.time.resolution = _to_float(self.header.get('TIMEDEL'))
        chunk.time.timesys = str(self.header.get('TIMESYS', 'UTC'))
        chunk.time.trefpos = self.header.get('TREFPOS', None)
        chunk.time.mjdref = self.header.get('MJDREF', self.header.get('MJDDATE'))
        self.logger.debug('End _finish_chunk_time')

    def _get_axis_length(self, for_axis):
        # try ZNAXIS first in order to get the size of the original image in case it was FITS compressed
        result = _to_int(self._sanitize(self.header.get(f'ZNAXIS{for_axis}')))
        if result is None:
            result = _to_int(self._sanitize(self.header.get(f'NAXIS{for_axis}')))
        if result is None:
            msg = f'Could not find axis length for axis {for_axis}'
            raise ValueError(msg)
        return result


class Hdf5WcsParser(WcsParser):
    """
    This class initializes an astropy.wcs instance with metadata from an Hdf5ObsBlueprint populated using an
    Hdf5Parser.
    """

    def __init__(self, blueprint, extension):
        """
        :param blueprint: ObsBlueprint
        """
        super().__init__(blueprint, extension)

    def _get_axis_index(self, keywords):
        result = self._axes['custom'][0]
        if 'RA' in keywords:
            result = self._axes['ra'][0]
        elif 'DEC' in keywords:
            result = self._axes['dec'][0]
        elif 'TIME' in keywords:
            result = self._axes['time'][0]
        elif 'FREQ' in keywords:
            result = self._axes['energy'][0]
        elif 'STOKES' in keywords:
            result = self._axes['polarization'][0]
        elif 'FLUX' in keywords:
            result = self._axes['observable'][0]
        return result
