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

__all__ = ['TimeUtil', 'EnergyUtil', 'EnergyConverter', 'ORIGIN']

# TODO both these are very bad, implement more sensibly
TARGET_TIMESYS = "UTC"
TARGET_CTYPE = "TIME"
TARGET_CUNIT = "d"
ORIGIN = 0


# CoordFunction1D, double/float
def pix2val(function, pix):
    refPix = float(function.ref_coord.pix)
    return function.ref_coord.val + function.delta * (pix - refPix)


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


#  TODO: Under Construction
#  Q: astropy units & Quanitities could be used here?
#  On first pass duplicating what is in the Java code so the
#  validations can be as close as possible mathematically
#  Plus it's not clear to me what advantage the Quantities/Units
#  have (other than handling conversions differently, as there needs
#  to be similar somersaulting in order to determine what units
#  need to be applied to the values in a particular spectralWCS. -
#  much of the work in the freqUnits & related arrays would still need
#  to be done.
class EnergyConverter():

    def __init__(self):
        self._CORE_SPECSYS = "BARYCENT"
        self._CORE_CTYPE = "WAVE"
        self._CORE_CUNIT = "m"
        self.BASE_UNIT_FREQ = "Hz"

        self.c = float(2.9979250e8) # m / sec
        self.h = float(6.62620e-27) # erg / sec
        self.eV = float(1.602192e-12) # erg

        self.allUnits = [];

        self.freqUnits = [ "Hz", "kHz", "MHz", "GHz" ]
        self.freqMult = [ 1.0, 1.0e3, 1.0e6, 1.0e9 ]
        self.enUnits = [ "eV", "keV", "MeV", "GeV" ]
        self.enMult = [ 1.0, 1.0e3, 1.0e6, 1.0e9 ]

        self.waveUnits = [ "m", "cm", "mm","um", "µm", "nm", "A" ]
        self.waveMult = [ 1.0, 1.0e-2, 1.0e-3, 1.0e-6, 1.0e-6, 1.0e-9, 1.0e-10 ]
        # todo: see if astropy units & quantities have representations of
        #  each of the waveUnits
        # self.bla = [u.m, u.cm, u.mm, u.]

    # Properties
    @property
    def CORE_CUNIT(self):
        """
        """
        return self._CORE_CUNIT

    @property
    def CORE_CTYPE(self):
        """
        """
        return self._CORE_CTYPE

        # // Lay out the actual units only once, then coalesce them.
        # static {
        #     final List<String> allUnitList = new ArrayList<String>(
        #             Arrays.asList(freqUnits));
        #     allUnitList.addAll(Arrays.asList(enUnits));
        #     allUnitList.addAll(Arrays.asList(waveUnits));
        #
        #     allUnits = allUnitList.toArray(new String[allUnitList.size()]);

        def getSupportedUnits(self):
            return self.allUnits;

        def convert(value, ctype, cunit):
            # TODO: check ctype instead of just relying on units
            return to_meters(value, cunit)

        def convertSpecsys(value, specsys):
            return value # noop

        # /**
        #  * Convert the energy value d from the specified units to wavelength in
        #  * meters.
        #  *
        #  * @param d
        #  * @param units
        #  * @return wavelength in meters
        #  */
        def to_meters(d, units):
            #  if a list of units is allowed,
            #  This wil turn into: discovering which array the wcs value unit
            #  is in, assigning that astropy unit to a variable.
            #  multiplying the wcs value (d) against the astropy unit (makng a quantity,)
            #  then applying and equivalence in the .to() function to do the
            #  conversion.
            #  It just seems like more work than using what is in the
            #  java code. :(
            try:
                # i = ArrayUtil.matches("^" + units + "$", self.freqUnits, True)
                if units in self.freqUnits:
                    i = self.freqUnits.index(units)
                    return freq_to_meters(d, i)

                # i = ArrayUtil.matches("^" + units + "$", self.enUnits, True)
                if units in self.enUnits:
                    i = self.enUnits.index(units)
                    return energy_to_meters(d, i)

                # i = ArrayUtil.matches("^" + units + "$", self.waveUnits, True)
                if units in self.waveUnits:
                    i = self.waveUnits.index(units)
                    return wavelength_to_meters(d, i)
            except ValueError:
                pass
                # TODO: what is the reporting structure for errors?

            # throw new IllegalArgumentException("Unknown units: " + units);

        # /**
        #  * Convert the energy value d from the specified units to frequency in Hz.
        #  *
        #  * @param d
        #  * @param units
        #  * @return frequency in Hz
        #  */
        def to_hz(self, d, units):
            # i = ArrayUtil.matches("^" + units + "$", self.freqUnits, True)
            try:
                i = self.freqUnits.index(units)
                if units in self.freqUnits:
                    return freq_to_hz(d, i)

                # i = ArrayUtil.matches("^" + units + "$", self.enUnits, True)
                i = self.enUnits.index(units)
                if units in self.enUnits:
                    return energy_to_hz(d, i)

                # i = ArrayUtil.matches("^" + units + "$", self.waveUnits, True)
                i = self.waveUnits.index(units)
                if units in self.waveUnits:
                    return wavelength_to_hz(d, i)
            except ValueError:
                pass
                # // unknown units

        # throw new IllegalArgumentException("unknown units: " + units);

        # /**
        #  * Compute the range of energy values to a wavelength width in meters.
        #  *
        #  * @param d1
        #  * @param d2
        #  * @param units
        #  * @return delta lambda in meters
        #  */
        # def to_delta_meters(self, d1, d2, units):
        #     w1 = to_meters(d1, units)
        #     w2 = to_meters(d2, units)
        #     return abs(w2 - w1)

        # /**
        #  * Compute the range of energy values to a frequency width in Hz.
        #  *
        #  * @param d1
        #  * @param d2
        #  * @param units
        #  * @return delta nu in Hz
        #  */
        # def to_delta_hz(self, d1, d2, units):
        #     f1 = to_hz(d1, units)
        #     f2 = to_hz(d2, units)
        #     return abs(f2 - f1)

        def freq_to_meters(self, d, i):
            nu = float(d * self.freqMult[i])
            return self.c / nu

        def energy_to_meters(self, d, i):
            e = float(self.eV * d * self.enMult[i])
            return float(self.c * self.h / e)

        def wavelength_to_meters(self, d, i):
            return float(d * self.waveMult[i])

        def freq_to_hz(self, d, i):
            return float(self.d * self.freqMult[i])

        def energy_to_hz(self, d, i):
            w = energy_to_meters(d, i)
            return float(self.c / w)

        def wavelength_to_hz(self, d, i):
            w = d * self.waveMult[i]
            return float(self.c / w)

#  TODO: Under Construction
class EnergyUtil():
    def __init__(self):
        pass

    def range1d_to_interval(self, temporal_wcs, axis_1d):
        """
        """
        a = float(axis_1d.start.val)
        b = float(axis_1d.end.val)
        conv = EnergyConverter()

        # String specsys = wcs.getSpecsys();
        # if (!EnergyConverter.CORE_SPECSYS.equals(specsys)) {
        #     a = conv.convertSpecsys(a, specsys);
        #     b = conv.convertSpecsys(b, specsys);
        # }

        ctype = temporal_wcs.axis.axis.ctype
        cunit = temporal_wcs.axis.axis.cunit
        # TODO: energy converter convert function not accessible?
        # if not ctype.startswith(str(EnergyConverter.CORE_CTYPE)) or EnergyConverter.CORE_CUNIT is not cunit:
        #         # log.debug("toInterval: converting " + a + cunit);
        #     a = conv.convert(a, EnergyConverter.CORE_CTYPE, cunit)
        #     # log.debug("toInterval: converting " + b + cunit);
        #     b = conv.convert(b, EnergyConverter.CORE_CTYPE, cunit)

        return shape.SubInterval(min(a, b), max(a, b))

    def function1d_to_interval(self, temporal_wcs):
        """
            needs Util.pix2val equivalent in here.
        """
        # throws
        # NoSuchKeywordException, WCSLibRuntimeException

        #  TODO: translate function needed?
        # ctype = wcs.axis.axis.ctype
        # if not ctype.startsWith(EnergyConverter.CORE_CTYPE):
        # # log.debug("toInterval: transform from " + ctype + " to " + EnergyConverter.CORE_CTYPE + "-???");
        #     kw = trans.translate(EnergyConverter.CORE_CTYPE + "-???"); // any linearization algorithm
        #     trans = new Transform(kw);

        # naxis = kw.getDoubleValue("NAXIS1"); // axis set to 1 above
        # function associated with naxis 1 ? TODO: is this right?
        naxis = temporal_wcs.axis.function.naxis
        p1 = 0.5
        p2 = naxis + 0.5

        wcsprm = Wcsprm()
        start_coord = np.array([[p1, p1]])
        end_coord = np.array([[p2, p2]])

        start = wcsprm.p2s(start_coord, ORIGIN)
        end = wcsprm.p2s(end_coord, ORIGIN)

        a = start['world'][0][0]
        b = end['world'][0][0]
        # log.debug("toInterval: wcslib returned " + a + start.units[0] + "," + b + end.units[0]);

        # String specsys = wcs.getSpecsys();
        # if not EnergyConverter.CORE_SPECSYS.equals(specsys)):
        # a = conv.convertSpecsys(a, specsys);
        # b = conv.convertSpecsys(b, specsys);
        # }

        # TODO: what needs to happen here?
        # where can units be pulled from?
        # // wcslib convert to WAVE-??? but units might be a multiple of EnergyConverter.CORE_UNIT
        # cunit = start.units[0]; # assume same as end.units[0]
        # if (!EnergyConverter.CORE_UNIT.equals(cunit)) {
        # log.debug("toInterval: converting " + a + " " + cunit);
        # a = conv.convert(a, EnergyConverter.CORE_CTYPE, cunit);
        # log.debug("toInterval: converting " + b + " " + cunit);
        # b = conv.convert(b, EnergyConverter.CORE_CTYPE, cunit);
        # }

        return shape.SubInterval(min(a, b), max(a, b))



