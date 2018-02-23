# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2018.                            (c) 2016.
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

"""
WCS Validation Utilities

"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from caom2 import shape
from astropy.wcs import Wcsprm
import numpy as np

__all__ = ['TimeUtil', 'EnergyUtil', 'ORIGIN']

# TODO both these are very bad, implement more sensibly
TARGET_TIMESYS = "UTC"
TARGET_CTYPE = "TIME"
TARGET_CUNIT = "d"
ORIGIN = 0


# CoordFunction1D, double/float
def pix2val(function, pix):
    ref_pix = float(function.ref_coord.pix)
    return function.ref_coord.val + function.delta * (pix - ref_pix)


class TimeUtil:
    def __init__(self):
        pass

    def range1d_to_interval(self, temporal_wcs, axis_1d):
        """
        """
        self.validate_wcs(temporal_wcs)

        # TODO: (comment pulled from Java code):
        # if mjdref has a value then the units of axis values could be any time
        # units, like days, hours, minutes, seconds, and smaller
        # since they are offsets from mjdref
        a = axis_1d.getStart().val
        b = axis_1d.getEnd().val
        if temporal_wcs.mjdref is not None:
            a += float(temporal_wcs.mjdref)
            b += float(temporal_wcs.mjdref)

        return shape.SubInterval(min(a,b),max(a,b))

    def function1d_to_interval(self, temporal_wcs, function_1d):
        """
        """
        try:
            self.validate_wcs(temporal_wcs)
            # // TODO: (comment pulled from Java code):
            # if mjdref has a value then the units of axis values could be any time
            # // units, like days, hours, minutes, seconds, and smaller since they are offsets
            # // from mjdref

            p1 = float(0.5)
            p2 = float(function_1d.naxis + 0.5)
            a = pix2val(function_1d, p1)
            b = pix2val(function_1d, p2)

            if temporal_wcs.mjdref is not None:
                a += float(temporal_wcs.mjdref)
                b += float(temporal_wcs.mjdref)

            return shape.SubInterval(min(a, b), max(a, b))

        except Exception:
            raise ValueError("Invalid function in Temporal WCS")

    def validate_wcs(self,  temporal_wcs):
        """
        """
        ctype = temporal_wcs.axis.axis.ctype
        sb = ""
        if ctype == TARGET_CTYPE and (temporal_wcs.timesys is None or temporal_wcs.timesys == TARGET_TIMESYS):
            pass
        elif ctype == TARGET_TIMESYS and temporal_wcs.timesys is None:
            pass
        else:
            sb = sb + "unexpected TIMESYS, CTYPE: " + temporal_wcs.timesys + "," + ctype

        cunit = temporal_wcs.axis.axis.cunit
        if TARGET_CUNIT != cunit:
            sb = sb + "unexpected CUNIT: " + cunit

        if len(sb) > 0:
            raise ValueError(sb)

# #  TODO: Under Construction
# #  Q: astropy units & Quanitities could be used here?
# #  On first pass duplicating what is in the Java code so the
# #  validations can be as close as possible mathematically
# #  Plus it's not clear to me what advantage the Quantities/Units
# #  have (other than handling conversions differently, as there needs
# #  to be similar somersaulting in order to determine what units
# #  need to be applied to the values in a particular spectralWCS. -
# #  much of the work in the freqUnits & related arrays would still need
# #  to be done.
# class EnergyConverter():
#
#     def __init__(self):
#         self._CORE_SPECSYS = "BARYCENT"
#         self._CORE_CTYPE = "WAVE"
#         self._CORE_CUNIT = "m"
#         self.BASE_UNIT_FREQ = "Hz"
#
#         self.c = float(2.9979250e8) # m / sec
#         self.h = float(6.62620e-27) # erg / sec
#         self.eV = float(1.602192e-12) # erg
#
#         self.allUnits = [];
#
#         self.freqUnits = [ "Hz", "kHz", "MHz", "GHz" ]
#         self.freqMult = [ 1.0, 1.0e3, 1.0e6, 1.0e9 ]
#         self.enUnits = [ "eV", "keV", "MeV", "GeV" ]
#         self.enMult = [ 1.0, 1.0e3, 1.0e6, 1.0e9 ]
#
#         self.waveUnits = [ "m", "cm", "mm","um", "µm", "nm", "A" ]
#         self.waveMult = [ 1.0, 1.0e-2, 1.0e-3, 1.0e-6, 1.0e-6, 1.0e-9, 1.0e-10 ]
#         # todo: see if astropy units & quantities have representations of
#         #  each of the waveUnits
#         # self.bla = [u.m, u.cm, u.mm, u.]
#
#     # Properties
#     @property
#     def CORE_CUNIT(self):
#         """
#         """
#         return self._CORE_CUNIT
#
#     @property
#     def CORE_CTYPE(self):
#         """
#         """
#         return self._CORE_CTYPE
#
#     @property
#     def CORE_SPECSYS(self):
#         """
#         """
#         return self._CORE_SPECSYS
#
#     # // Lay out the actual units only once, then coalesce them.
#     # static {
#     #     final List<String> allUnitList = new ArrayList<String>(
#     #             Arrays.asList(freqUnits));
#     #     allUnitList.addAll(Arrays.asList(enUnits));
#     #     allUnitList.addAll(Arrays.asList(waveUnits));
#     #
#     #     allUnits = allUnitList.toArray(new String[allUnitList.size()]);
#
#     def getSupportedUnits(self):
#         return self.allUnits;
#
#     def convert(self, value, ctype, cunit):
#         # TODO: check ctype instead of just relying on units
#         return self.to_meters(value, cunit)
#
#     def convert_specsys(self, value, specsys):
#         return value # noop
#
#     # /**
#     #  * Convert the energy value d from the specified units to wavelength in
#     #  * meters.
#     #  *
#     #  * @param d
#     #  * @param units
#     #  * @return wavelength in meters
#     #  */
#     def to_meters(self, d, units):
#         #  if a list of units is allowed,
#         #  This wil turn into: discovering which array the wcs value unit
#         #  is in, assigning that astropy unit to a variable.
#         #  multiplying the wcs value (d) against the astropy unit (makng a quantity,)
#         #  then applying and equivalence in the .to() function to do the
#         #  conversion.
#         #  It just seems like more work than using what is in the
#         #  java code. :(
#         # try:
#         # i = ArrayUtil.matches("^" + units + "$", self.freqUnits, True)
#         if units in self.freqUnits:
#             i = self.freqUnits.index(units)
#             return self.freq_to_meters(d, i)
#
#         # i = ArrayUtil.matches("^" + units + "$", self.enUnits, True)
#         if units in self.enUnits:
#             i = self.enUnits.index(units)
#             return self.energy_to_meters(d, i)
#
#         # i = ArrayUtil.matches("^" + units + "$", self.waveUnits, True)
#         if units in self.waveUnits:
#             i = self.waveUnits.index(units)
#             return self.wavelength_to_meters(d, i)
#         # except ValueError:
#         #     pass
#
#
#         # throw new IllegalArgumentException("Unknown units: " + units);
#
#     # /**
#     #  * Convert the energy value d from the specified units to frequency in Hz.
#     #  *
#     #  * @param d
#     #  * @param units
#     #  * @return frequency in Hz
#     #  */
#     def to_hz(self, d, units):
#         # i = ArrayUtil.matches("^" + units + "$", self.freqUnits, True)
#         try:
#             i = self.freqUnits.index(units)
#             if units in self.freqUnits:
#                 return self.freq_to_hz(d, i)
#
#             # i = ArrayUtil.matches("^" + units + "$", self.enUnits, True)
#             i = self.enUnits.index(units)
#             if units in self.enUnits:
#                 return self.energy_to_hz(d, i)
#
#             # i = ArrayUtil.matches("^" + units + "$", self.waveUnits, True)
#             i = self.waveUnits.index(units)
#             if units in self.waveUnits:
#                 return self.wavelength_to_hz(d, i)
#         except ValueError:
#             pass
#             #  unknown units
#
#
#
#     def freq_to_meters(self, d, i):
#         nu = float(d * self.freqMult[i])
#         return self.c / nu
#
#     def energy_to_meters(self, d, i):
#         e = float(self.eV * d * self.enMult[i])
#         return float(self.c * self.h / e)
#
#     def wavelength_to_meters(self, d, i):
#         return float(d * self.waveMult[i])
#
#     def freq_to_hz(self, d, i):
#         return float(self.d * self.freqMult[i])
#
#     def energy_to_hz(self, d, i):
#         w = self.energy_to_meters(d, i)
#         return float(self.c / w)
#
#     def wavelength_to_hz(self, d, i):
#         w = d * self.waveMult[i]
#         return float(self.c / w)


class EnergyUtil():
    def __init__(self):
        pass

    def range1d_to_interval(self, temporal_wcs, range_1d):
        """
        """
        a = float(range_1d.start.val)
        b = float(range_1d.end.val)
        #  The energy converter work done in the Java code is skipped here.
        #  Doing it here introduced some precision errors which lead to false invalids.
        #  Ignoring the units for validation sounds like it's ok, as long as the
        #  same values come out of the p2s, s2p calculations in the main validator
        #  code. TODO: is this sufficient for validation?
        # conv = EnergyConverter()
        #
        # ctype = temporal_wcs.axis.axis.ctype
        # cunit = temporal_wcs.axis.axis.cunit
        # if not ctype.startswith(str(EnergyConverter.CORE_CTYPE)) or EnergyConverter.CORE_CUNIT != cunit:
        #     a = conv.convert(a, conv.CORE_CTYPE, cunit)
        #     b = conv.convert(b, conv.CORE_CTYPE, cunit)

        return shape.SubInterval(min(a, b), max(a, b))

    def function1d_to_interval(self, temporal_wcs):
        """
        """
        naxis = temporal_wcs.axis.function.naxis
        p1 = 0.5
        p2 = naxis + 0.5
        return shape.SubInterval(p1, p2)


class PolarizationWcsUtil():
    def _get_xrange(self, range):
        if range is not None:
            lb = int(round(range.start.val))
            ub = int(round(range.end.val))
            return xrange(lb, ub+1)
        return None

    @staticmethod
    def get_keys(range):
        """
        Examines the lower bound (lb) and upper bound (ub) of
        PolarizationWCS.axis.range and returns range(lb, ub+1) if a range is
        defined, else returns None. Since Python range iterates over [lb,ub),
        the returned range is ub+1 to ensure that ub is included in the
        range iteration.
        """
        return PolarizationWcsUtil()._get_xrange(range)

    @staticmethod
    def get_ranges_from_bounds(bounds):
        """
        Examines the ranges in PolarizationWCS.axis.bounds and returns the
        list of ranges in the bounds if the bounds is defined, else returns
        an empty list. The upper bound of each range is incremented by 1
        (refer to comments for get_range above)
        """
        ranges = []
        if bounds is not None:
            samples = bounds.samples
            if samples is not None:
                for sample in samples:
                    ranges.append(PolarizationWcsUtil()._get_xrange(sample))
        return ranges

    @staticmethod
    def get_range_from_function(function):
        """
        Examines the ranges in PolarizationWCS.axis.function and returns a
        range from 1 to Naxis+1 if the function is defined, else returns
        None. The upper bound of the range is incremented by 1 (refer to
        comments for get_range above)
        """
        if function is not None:
            if function.naxis >= 1:
                return xrange(1, function.naxis + 1)
            else:
                raise ValueError(
                    'Invalid naxis value: {}'.format(function.naxis))
        return None
