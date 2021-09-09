# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2019.                            (c) 2019.
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

from caom2 import shape, chunk, plane
import logging
import sys


logger = logging.getLogger(__name__)

__all__ = ['TimeUtil', 'CustomAxisUtil', 'EnergyUtil', 'ORIGIN']

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

    @staticmethod
    def range1d_to_interval(temporal_wcs, range_1d):
        TimeUtil.validate_wcs(temporal_wcs)

        # TODO: (comment pulled from Java code):
        # if mjdref has a value then the units of axis values could be any time
        # units, like days, hours, minutes, seconds, and smaller
        # since they are offsets from mjdref
        a = range_1d.start.val
        b = range_1d.end.val
        if b < a:
            raise ValueError("range.end not >= range.start in Temporal WCS")

        if temporal_wcs.mjdref is not None:
            a += float(temporal_wcs.mjdref)
            b += float(temporal_wcs.mjdref)

        return shape.SubInterval(min(a, b), max(a, b))

    @staticmethod
    def function1d_to_interval(temporal_wcs, function_1d):
        try:
            TimeUtil.validate_wcs(temporal_wcs)
            # // TODO: (comment pulled from Java code):
            # if mjdref has a value then the units of axis values could be
            # any time
            # // units, like days, hours, minutes, seconds, and smaller
            # // since they are offsets from mjdref

            # PD - 16-04-20
            # technically there is nothing wrong with a WCS axis that
            # decreases in coord values while increasing in pixel values
            # (it's just a line with negative slope).
            #
            # the computation of the time bounds interval sorts out the
            # min/max so it also doesn't care about the "direction"

            p1 = float(0.5)
            p2 = float(function_1d.naxis + 0.5)
            a = pix2val(function_1d, p1)
            b = pix2val(function_1d, p2)

            if temporal_wcs.mjdref is not None:
                a += float(temporal_wcs.mjdref)
                b += float(temporal_wcs.mjdref)

            return shape.SubInterval(min(a, b), max(a, b))

        except Exception as ex:
            raise ValueError(
                f"Invalid function in Temporal WCS: {repr(ex)}")

    @staticmethod
    def validate_wcs(temporal_wcs):
        ctype = temporal_wcs.axis.axis.ctype
        sb = ""
        if ctype == TARGET_CTYPE \
                and (temporal_wcs.timesys is None
                     or temporal_wcs.timesys == TARGET_TIMESYS):
            pass
        elif ctype == TARGET_TIMESYS and temporal_wcs.timesys is None:
            pass
        else:
            sb = "unexpected TIMESYS, CTYPE: {},{}".format(
                temporal_wcs.timesys, ctype)

        cunit = temporal_wcs.axis.axis.cunit
        if TARGET_CUNIT != cunit:
            sb = sb + f"unexpected CUNIT: {cunit}"

        if len(sb) > 0:
            raise ValueError(sb)


class CustomAxisUtil:
    """
    Utility class for Custom calculations. Ported from the Java version
    of CustomAxisUtil. Additional functions were ported from
    caom2/compute/Util.java
    """

    ctype_cunit_map = {"FARADAY": "rad/m**2", "RM": "rad/m**2"}

    def __init__(self):
        pass

    @staticmethod
    def _merge_into_list(si, samples, union_scale):
        """Merge an Interval into a list of Intervals

        :param si An Interval
        :param samples A list of Intervals
        :param union_scale A float
        """
        snew = si
        if len(samples) > 0:
            f = union_scale * (si.upper - si.lower)
            a = si.lower - f
            b = si.upper + f
            tmp = []
            # Find intervals that overlap the new one, move from samples -> tmp
            for sample in samples:
                f = union_scale * (sample.upper - sample.lower)
                c = sample.lower - f
                d = sample.upper + f
                # [a,b] U [c,d]
                if b >= c and d >= a:
                    # there is overlap
                    tmp.append(sample)
                    samples.pop(0)

            # Merge all overlapping sub-intervals
            # compute the outer bounds of the sub-intervals
            if len(tmp) > 0:
                lb = si.lower
                ub = si.upper
                for s in tmp:
                    if lb > s.lower:
                        lb = s.lower

                    if ub < s.upper:
                        ub = s.upper

                snew = shape.Interval(lb, ub)

        # Insert new sub to preserve order
        added = False
        for i in range(len(samples)):
            sample = samples[i]
            if snew.lower < sample.lower:
                samples.insert(i, snew)
                added = True
                break

        if not added:
            samples.append(snew)

    @staticmethod
    def _get_num_pixels(axis, use_func):
        """
        Count the number of pixels
        :param axis: A CoordAxis1D
        :param use_func: A boolean
        :return: A float
        """
        range = axis.range
        bounds = axis.bounds
        function = axis.function
        if range is not None:
            return abs(range.end.pix - range.start.pix)

        if bounds is not None:
            # Count number of distinct bins
            bins = []
            for cr in bounds.samples:
                si = shape.Interval(cr.start.pix, cr.end.pix)
                CustomAxisUtil._merge_into_list(si, bins, float(0.0))

            ret = float(0.0)
            for si in bins:
                ret += abs(si.upper - si.lower)

            return ret

        if use_func and function is not None:
            return function.naxis

        return float(0.0)

    @staticmethod
    def _use_chunk(atype, ptype, ctype, matches):
        """
        Determines whether to use the Chunk.
        :param atype: Product_type in Artifact
        :param ptype: Product_type in Plane
        :param ctype: Product_type in Chunk
        :param matches: Product_type to match
        :return: True if use Chunk, False otherwise
        """
        if matches is None:
            return False

        if ctype is not None and ctype != matches:
            logger.debug(
                f"use_chunk=False: Chunk.product_type={ctype}")
            return False

        if ptype is not None and ptype != matches:
            logger.debug(
                f"use_chunk=False: Part.product_type={ptype}")
            return False

        if atype == matches:
            logger.debug(
                f"use_chunk=True: Artifact.product_type={atype}")
            return True

        logger.debug("use_chunk=False: product_type={},{},{}".
                     format(atype, ptype, ctype))
        return False

    @staticmethod
    def val2pix(wcs, func, val):
        """ val2pix calculates pixel from value

        :param wcs A CustomWCS
        :param func A CoordFunction1D
        :param val A float
        :return A float
        """
        CustomAxisUtil.validate_wcs(wcs)
        ref_val = func.ref_coord.val
        return func.ref_coord.pix + (val - ref_val) / func.delta

    @staticmethod
    def function1d_to_interval(wcs, func):
        """ function1d_to_interval calculates interval for CoordFunction1D

        :param wcs A CustomWCS
        :param r A CoordFunction1D
        :return An Interval
        """
        CustomAxisUtil.validate_wcs(wcs)
        if func.delta == 0.0 and func.naxis > 1:
            raise ValueError(
                "Invalid CoordFunction1D: found {} pixels and delta = 0.0".
                format(func.naxis))

        p1 = 0.5
        p2 = float(func.naxis) + 0.5
        a = CustomAxisUtil.val2pix(wcs, func, p1)
        b = CustomAxisUtil.val2pix(wcs, func, p2)

        return shape.Interval(min(a, b), max(a, b))

    @staticmethod
    def range1d_to_interval(wcs, r):
        """ range1d_to_interval calculates interval for CoordRange1D

        :param wcs A CustomWCS
        :param r A CoordRange1D
        :return An Interval
        """
        CustomAxisUtil.validate_wcs(wcs)
        np = abs(r.start.pix - r.end.pix)
        a = r.start.val
        b = r.end.val
        delta = abs(b - a)
        if delta == 0.0 and np > 1.0:
            raise ValueError(
                "Invalid CoordRange1D: found {} + pixels and delta = 0.0 in \
                [{},{}]".format(np, a, b)
            )

        return shape.Interval(min(a, b), max(a, b))

    @staticmethod
    def _chose_product_type(artifacts):
        ret = None
        for a in artifacts:
            if chunk.ProductType.SCIENCE == a.product_type:
                return chunk.ProductType.SCIENCE

            if chunk.ProductType.CALIBRATION == a.product_type:
                return chunk.ProductType.CALIBRATION

            for p_key in a.parts:
                p = a.parts[p_key]
                if chunk.ProductType.SCIENCE == p.product_type:
                    return chunk.ProductType.SCIENCE

                if chunk.ProductType.CALIBRATION == p.product_type:
                    return chunk.ProductType.CALIBRATION

                for c in p.chunks:
                    if chunk.ProductType.SCIENCE == c.product_type:
                        return chunk.ProductType.SCIENCE

                    if chunk.ProductType.CALIBRATION == c.product_type:
                        ret = chunk.ProductType.CALIBRATION

        return ret

    @staticmethod
    def _get_ctype(artifacts, product_type):
        first_ctype = None
        for a in artifacts:
            for p_key in a.parts:
                p = a.parts[p_key]
                for c in p.chunks:
                    if c.custom is not None and \
                            CustomAxisUtil._use_chunk(
                                a.product_type, p.product_type,
                                c.product_type, product_type):
                        current_ctype = c.custom.axis.axis.ctype
                        if first_ctype is None:
                            if current_ctype in CustomAxisUtil.ctype_cunit_map:
                                first_ctype = current_ctype
                            else:
                                raise ValueError("Unsupported CTYPE: {}".
                                                 format(current_ctype))

                        if current_ctype != first_ctype:
                            raise ValueError(
                                "CTYPE must be the same across all Artifacts. \
                                Found: {} and {}".format(current_ctype,
                                                         first_ctype))

        return first_ctype

    @staticmethod
    def compute(artifacts):
        product_type = CustomAxisUtil._chose_product_type(artifacts)
        axis_ctype = CustomAxisUtil._get_ctype(artifacts, product_type)
        if axis_ctype is not None:
            c = plane.CustomAxis(axis_ctype)
            if product_type is not None:
                c.bounds = CustomAxisUtil.compute_bounds(
                    artifacts, product_type, axis_ctype)
                if c.dimension is None:
                    c.dimension = CustomAxisUtil.compute_dimension_from_wcs(
                        c.bounds, artifacts, product_type, axis_ctype)
            return c
        else:
            # No ctype found for chosen product type
            return None

    @staticmethod
    def compute_bounds(artifacts, product_type, expected_ctype):
        """ Compute bounds.

        :param artifacts List of Artifacts
        :param product_type A Product_type
        :param expected_ctype Expected CTYPE
        :return An Interval
        """
        union_scale = float(0.02)
        subs = []
        for a in artifacts:
            for p_key in a.parts:
                p = a.parts[p_key]
                for c in p.chunks:
                    if c is not None and c.custom is not None and \
                            CustomAxisUtil._use_chunk(
                                a.product_type, p.product_type,
                                c.product_type, product_type):
                        current_ctype = c.custom.axis.axis.ctype
                        if current_ctype is None or \
                                current_ctype != expected_ctype:
                            raise ValueError(
                                "CTYPE must be the same across all Artifacts. \
                                Found: {}. Expected: {}".format(
                                    current_ctype, expected_ctype))
                        else:
                            range = c.custom.axis.range
                            bounds = c.custom.axis.bounds
                            function = c.custom.axis.function
                            if range is not None:
                                s = CustomAxisUtil.range1d_to_interval(
                                    c.custom, range)
                                logger.debug(
                                    "[compute_bounds] range -> sub: {}".
                                    format(s))
                                CustomAxisUtil._merge_into_list(
                                    s, subs, union_scale)
                            elif bounds is not None:
                                for cr in bounds.samples:
                                    s = CustomAxisUtil.range1d_to_interval(
                                        c.custom, cr)
                                    logger.debug(
                                        "[compute_bounds] bounds -> sub: {}".
                                        format(s))
                                    CustomAxisUtil._merge_into_list(
                                        s, subs, union_scale)
                            elif function is not None:
                                s = CustomAxisUtil.function1d_to_interval(
                                    c.custom, function)
                                logger.debug(
                                    "[compute_bounds] function -> sub: {}".
                                    format(s))
                                CustomAxisUtil._merge_into_list(
                                    s, subs, union_scale)

        if len(subs) == 0:
            return None

        # Compute the outer bounds of the sub intervals
        lb = sys.float_info.max
        ub = sys.float_info.min
        for sub in subs:
            lb = min(lb, sub.lower)
            ub = max(ub, sub.upper)

        return shape.Interval(lb, ub, subs)

    @staticmethod
    def compute_dimension_from_wcs(bounds, artifacts,
                                   product_type, expected_ctype):
        """ Compute dimensionality (number of pixels).

        :param bounds A sampled interval
        :param artifacts List of Artifacts
        :param product_type A Product_type
        :param expected_ctype Expected CTYPE
        :return long
        """
        logger.debug(f"compute_dimension_from_wcs: {bounds}")
        if bounds is None:
            return None

        # Pick the WCS with the largest pixel size
        sw = None
        scale = float(0.0)
        num = 0
        for a in artifacts:
            for p_key in a.parts:
                p = a.parts[p_key]
                for c in p.chunks:
                    if c is not None and c.custom is not None and \
                            CustomAxisUtil._use_chunk(
                                a.product_type, p.product_type,
                                c.product_type, product_type):
                        current_ctype = c.custom.axis.axis.ctype
                        if current_ctype is None or \
                                current_ctype != expected_ctype:
                            raise ValueError(
                                "CTYPE must be the same across all Artifacts. \
                                Found: {}. Expected: {}".format(
                                    current_ctype, expected_ctype))
                        else:
                            num += 1
                            ss = abs(c.custom.axis.function.delta)
                            if ss >= scale:
                                scale = ss
                                sw = c.custom

        if sw is None:
            return None

        if sw.axis.function is None:
            return None

        if num == 1:
            return sw.axis.function.naxis

        x1 = CustomAxisUtil.val2pix(sw, sw.axis.function, bounds.lower)
        x2 = CustomAxisUtil.val2pix(sw, sw.axis.function, bounds.upper)

        logger.debug(f"compute_dimension_from_wcs: {x1},{x2}")
        return int(round(abs(x2 - x1)))

    @staticmethod
    def compute_dimension_from_range_bounds(
            artifacts, product_type, expected_ctype):
        """ Compute dimensionality (number of pixels).

        :param artifacts List of Artifacts
        :param product_type A Product_type
        :param expected_ctype Expected CTYPE
        :return long
        """
        # Assumption: all ...pixels are distinct so just add up the
        # number of pixels
        num_pixels = float(0.0)
        for a in artifacts:
            for p_key in a.parts:
                p = a.parts[p_key]
                for c in p.chunks:
                    if CustomAxisUtil._use_chunk(
                            a.product_type, p.product_type,
                            c.product_type, product_type):
                        current_ctype = c.custom.axis.axis.ctype
                        if current_ctype is None or \
                                current_ctype != expected_ctype:
                            raise ValueError(
                                "CTYPE must be the same across all Artifacts. \
                                Found: {}. Expected: {}".format(
                                    current_ctype, expected_ctype))
                        else:
                            n = CustomAxisUtil._get_num_pixels(
                                c.custom.axis, False)
                            num_pixels += n

        if num_pixels > 0.0:
            logger.debug("compute_dimension_from_range_bounds: {}".format(
                num_pixels))
            return int(num_pixels)

        logger.debug("compute_dimension_from_range_bounds: None")
        return None

    @staticmethod
    def validate_wcs(custom_wcs):
        ctype = custom_wcs.axis.axis.ctype
        raw_cunit = custom_wcs.axis.axis.cunit
        cunit = CustomAxisUtil()._normalize_unit(raw_cunit)
        map_cunit = None
        if ctype in CustomAxisUtil.ctype_cunit_map:
            map_cunit = CustomAxisUtil.ctype_cunit_map[ctype]

        if map_cunit is None:
            raise ValueError(
               f"Invalid CTYPE: {ctype}")

        if map_cunit != cunit:
            raise ValueError(
                "Invalid CUNIT for CTYPE: {}. Expected: {}. Found {} \
                (normalized, raw: {})".format(
                    ctype, map_cunit, cunit, raw_cunit))

    @staticmethod
    def _normalize_unit(raw_cunit):
        normalized_unit = raw_cunit
        if "^" in raw_cunit:
            normalized_unit = raw_cunit.replace("^", "**")
            logger.debug("normalized unit: {} to {}".format(
                raw_cunit, normalized_unit))

        return normalized_unit


class EnergyUtil:
    def __init__(self):
        pass

    @staticmethod
    def range1d_to_interval(range_1d):
        a = float(range_1d.start.val)
        b = float(range_1d.end.val)
        #  The energy converter work done in the Java code is skipped here.
        #  Doing it here introduced some precision errors which lead to false
        #  invalids. Ignoring the units for validation sounds like it's ok,
        #  as long as the same values come out of the p2s, s2p calculations
        #  in the main validator code. It's assumed that doing the
        #  conversions in the native units is sufficient, as long as the
        #  values are the same after p2s -> s2p is done.

        return shape.SubInterval(min(a, b), max(a, b))

    @staticmethod
    def function1d_to_interval(temporal_wcs):
        naxis = temporal_wcs.axis.function.naxis
        p1 = 0.5
        p2 = naxis + 0.5
        return shape.SubInterval(p1, p2)


class PolarizationWcsUtil:
    def _get_range(self, from_range):
        if from_range is not None:
            lb = int(round(from_range.start.val))
            ub = int(round(from_range.end.val))
            return range(lb, ub+1)
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
        return PolarizationWcsUtil()._get_range(range)

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
                    ranges.append(PolarizationWcsUtil()._get_range(sample))
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
                return range(1, function.naxis + 1)
            else:
                raise ValueError(
                    f'Invalid naxis value: {function.naxis}')
        return None
