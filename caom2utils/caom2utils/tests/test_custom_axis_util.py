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
#  General Public License for           Générale Publique GNU AfferoF
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

from caom2utils import wcs_util
from caom2 import ReleaseType, Artifact, Part, Chunk, plane, caom_util, \
    chunk, CoordAxis1D, CoordBounds1D, CoordFunction1D, CoordRange1D, \
    Interval, RefCoord, wcs
from caom2.caom_util import TypedList, TypedOrderedDict
import pytest
import unittest

single_test = False


# Custom WCS Util tests
@pytest.mark.skipif(single_test, reason='Single test mode')
class CustomAxisUtilTests(unittest.TestCase):
    def test_function1d_to_interval(self):
        # validate wcs
        wcs = CustomTestUtil.good_wcs()
        wcs_util.CustomAxisUtil.validate_wcs(wcs)
        # bad ctype
        wcs = CustomTestUtil.bad_ctype_wcs()
        with pytest.raises(ValueError) as ex:
            wcs_util.CustomAxisUtil.validate_wcs(wcs)
        assert ('Invalid CTYPE:' in str(ex.value))
        # bad cunit
        wcs = CustomTestUtil.bad_cunit_wcs()
        with pytest.raises(ValueError) as ex:
            wcs_util.CustomAxisUtil.validate_wcs(wcs)
        assert ('Invalid CUNIT for CTYPE:' in str(ex.value))

    def test_val2pix(self):
        # happy path
        wcs = CustomTestUtil.good_wcs()
        naxis = int(100)
        delta = -0.01
        ref_coord = RefCoord(0.0, 0.0)
        func = CoordFunction1D(naxis, delta, ref_coord)
        val = 0.1
        pix = wcs_util.CustomAxisUtil.val2pix(wcs, func, val)
        expected_pix = -10.0
        self.assertEqual(pix, expected_pix)

    def test_function1d_to_interval_happy_path(self):
        # happy path
        wcs = CustomTestUtil.good_wcs()
        naxis = int(100)
        delta = -0.2
        ref_coord = RefCoord(0.0, 0.0)
        function_1d = CoordFunction1D(naxis, delta, ref_coord)
        actual_interval = wcs_util.CustomAxisUtil.function1d_to_interval(
            wcs, function_1d)
        expected_interval = Interval(-502.5, -2.5)
        self.assertEqual(expected_interval.lower, actual_interval.lower)
        self.assertEqual(expected_interval.upper, actual_interval.upper)
        self.assertEqual(None, actual_interval.samples)
        # function_1d.delta == 0.0 && function_1d.naxis > 1
        naxis = int(100)
        delta = 0.0
        ref_coord = RefCoord(0.0, 0.0)
        function_1d = CoordFunction1D(naxis, delta, ref_coord)
        with pytest.raises(ValueError) as ex:
            wcs_util.CustomAxisUtil.function1d_to_interval(wcs, function_1d)
        assert ('Invalid CoordFunction1D:' in str(ex.value))

    def test_range1d_to_interval(self):
        # happy path
        wcs = CustomTestUtil.good_wcs()
        start = RefCoord(float(0.9), float(1.1))
        end = RefCoord(float(10.9), float(11.1))
        range_1d = CoordRange1D(start, end)
        actual_interval = wcs_util.CustomAxisUtil.range1d_to_interval(
            wcs, range_1d)
        expected_interval = Interval(1.1, 11.1)
        self.assertEqual(expected_interval.lower, actual_interval.lower)
        self.assertEqual(expected_interval.upper, actual_interval.upper)
        self.assertEqual(None, actual_interval.samples)
        # function_1d.delta == 0.0 && function_1d.naxis > 1
        start = RefCoord(float(0.9), float(1.1))
        end = RefCoord(float(10.9), float(1.1))
        range_1d = CoordRange1D(start, end)
        with pytest.raises(ValueError) as ex:
            wcs_util.CustomAxisUtil.range1d_to_interval(wcs, range_1d)
        assert ('Invalid CoordRange1D:' in str(ex.value))

    def test_compute_dimension_from_range_bounds(self):
        # user_chunk = False, matches is None
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.SCIENCE
        test_chunk.custom = CustomTestUtil.good_wcs()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.SCIENCE
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.SCIENCE
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        product_type = None
        expected_ctype = "RM"
        actual_num_pixels = \
            wcs_util.CustomAxisUtil.compute_dimension_from_range_bounds(
                artifacts, product_type, expected_ctype)
        expected_num_pixels = None
        self.assertEqual(expected_num_pixels, actual_num_pixels)
        # user_chunk = False, ctype not match
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.SCIENCE
        test_chunk.custom = CustomTestUtil.good_wcs()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.SCIENCE
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.SCIENCE
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        product_type = chunk.ProductType.CALIBRATION
        expected_ctype = "RM"
        actual_num_pixels = \
            wcs_util.CustomAxisUtil.compute_dimension_from_range_bounds(
                artifacts, product_type, expected_ctype)
        expected_num_pixels = None
        self.assertEqual(expected_num_pixels, actual_num_pixels)
        # user_chunk = False, ptype not match
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.CALIBRATION
        test_chunk.custom = CustomTestUtil.good_wcs()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.SCIENCE
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.SCIENCE
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        product_type = chunk.ProductType.CALIBRATION
        expected_ctype = "RM"
        actual_num_pixels = \
            wcs_util.CustomAxisUtil.compute_dimension_from_range_bounds(
                artifacts, product_type, expected_ctype)
        expected_num_pixels = None
        self.assertEqual(expected_num_pixels, actual_num_pixels)
        # user_chunk = False, atype not match
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.CALIBRATION
        test_chunk.custom = CustomTestUtil.good_wcs()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.CALIBRATION
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.SCIENCE
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        product_type = chunk.ProductType.CALIBRATION
        expected_ctype = "RM"
        actual_num_pixels = \
            wcs_util.CustomAxisUtil.compute_dimension_from_range_bounds(
                artifacts, product_type, expected_ctype)
        expected_num_pixels = None
        self.assertEqual(expected_num_pixels, actual_num_pixels)
        # user_chunk = True, current_type != expected_ctype
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.CALIBRATION
        test_chunk.custom = CustomTestUtil.good_wcs_with_function()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.CALIBRATION
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.CALIBRATION
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        product_type = chunk.ProductType.CALIBRATION
        expected_ctype = "FARADAY"
        with pytest.raises(ValueError) as ex:
            wcs_util.CustomAxisUtil.compute_dimension_from_range_bounds(
                artifacts, product_type, expected_ctype)
        assert ('CTYPE must be the same across all Artifacts' in str(ex.value))
        # user_chunk = True, get_num_pixels: range is not None
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.CALIBRATION
        test_chunk.custom = CustomTestUtil.good_wcs_with_range()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.CALIBRATION
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.CALIBRATION
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        product_type = chunk.ProductType.CALIBRATION
        expected_ctype = "RM"
        actual_num_pixels = \
            wcs_util.CustomAxisUtil.compute_dimension_from_range_bounds(
                artifacts, product_type, expected_ctype)
        expected_num_pixels = 10
        self.assertEqual(expected_num_pixels, actual_num_pixels)
        # user_chunk = True, get_num_pixels: bounds with 3 samples that
        # traverses _merge_into_list completely
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.CALIBRATION
        test_chunk.custom = CustomTestUtil.good_wcs_with_bounds_3_samples()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.CALIBRATION
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.CALIBRATION
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        product_type = chunk.ProductType.CALIBRATION
        expected_ctype = "RM"
        actual_num_pixels = \
            wcs_util.CustomAxisUtil.compute_dimension_from_range_bounds(
                artifacts, product_type, expected_ctype)
        expected_num_pixels = 11
        self.assertEqual(expected_num_pixels, actual_num_pixels)
        # user_chunk = True, range = None, bounds = None, use_func and
        # function = None
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.CALIBRATION
        test_chunk.custom = CustomTestUtil.good_wcs()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.CALIBRATION
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.CALIBRATION
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        product_type = chunk.ProductType.CALIBRATION
        expected_ctype = "RM"
        actual_num_pixels = \
            wcs_util.CustomAxisUtil.compute_dimension_from_range_bounds(
                artifacts, product_type, expected_ctype)
        expected_num_pixels = None
        self.assertEqual(expected_num_pixels, actual_num_pixels)

    def test_compute_dimension_from_wcs(self):
        # bounds is None
        bounds = None
        artifacts = None
        product_type = None
        expected_ctype = None
        actual_dimension = wcs_util.CustomAxisUtil.compute_dimension_from_wcs(
            bounds, artifacts, product_type, expected_ctype)
        expected_dimension = None
        self.assertEqual(expected_dimension, actual_dimension)
        # bounds is not None, user_chunk = True, current_type != expected_ctype
        bounds = Interval(1.1, 11.1)
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.CALIBRATION
        test_chunk.custom = CustomTestUtil.good_wcs_with_function()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.CALIBRATION
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.CALIBRATION
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        product_type = chunk.ProductType.CALIBRATION
        expected_ctype = "FARADAY"
        with pytest.raises(ValueError) as ex:
            wcs_util.CustomAxisUtil.compute_dimension_from_wcs(
                bounds, artifacts, product_type, expected_ctype)
        assert ('CTYPE must be the same across all Artifacts' in str(ex.value))
        # bounds is not None, user_chunk = True, current_type is not None and
        # current_type == expected_ctype, ss >= scale, num = 1
        bounds = Interval(1.1, 11.1)
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.CALIBRATION
        test_chunk.custom = CustomTestUtil.good_wcs_with_negative_delta()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.CALIBRATION
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.CALIBRATION
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        product_type = chunk.ProductType.CALIBRATION
        expected_ctype = "RM"
        actual_dimension = wcs_util.CustomAxisUtil.compute_dimension_from_wcs(
            bounds, artifacts, product_type, expected_ctype)
        expected_dimension = 200
        self.assertEqual(expected_dimension, actual_dimension)
        # bounds is not None, user_chunk = False, sw = None
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.SCIENCE
        test_chunk.custom = CustomTestUtil.good_wcs()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.SCIENCE
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.SCIENCE
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        product_type = None
        expected_ctype = "RM"
        actual_dimension = wcs_util.CustomAxisUtil.compute_dimension_from_wcs(
            bounds, artifacts, product_type, expected_ctype)
        expected_dimension = None
        self.assertEqual(expected_dimension, actual_dimension)
        # bounds is not None, user_chunk = True, current_type is not None and
        # current_type == expected_ctype, ss >= scale, num = 2
        bounds = Interval(1.1, 11.1)
        test_chunk1 = Chunk()
        test_chunk1.product_type = chunk.ProductType.CALIBRATION
        test_chunk1.custom = CustomTestUtil.good_wcs_with_negative_delta()
        test_chunk2 = Chunk()
        test_chunk2.product_type = chunk.ProductType.CALIBRATION
        test_chunk2.custom = CustomTestUtil.good_wcs()
        test_chunks = TypedList(Chunk, test_chunk1, test_chunk2)
        part_name = "test_part"
        part_product_type = chunk.ProductType.CALIBRATION
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.CALIBRATION
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        product_type = chunk.ProductType.CALIBRATION
        expected_ctype = "RM"
        actual_dimension = wcs_util.CustomAxisUtil.compute_dimension_from_wcs(
            bounds, artifacts, product_type, expected_ctype)
        expected_dimension = 500
        self.assertEqual(expected_dimension, actual_dimension)

    def test_compute_bounds(self):
        # user_chunk = False
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.SCIENCE
        test_chunk.custom = CustomTestUtil.good_wcs()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.SCIENCE
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.SCIENCE
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        product_type = None
        expected_ctype = "RM"
        actual_bounds = wcs_util.CustomAxisUtil.compute_bounds(
            artifacts, product_type, expected_ctype)
        expected_bounds = None
        self.assertEqual(expected_bounds, actual_bounds)
        # user_chunk = True, current_type != expected_ctype
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.CALIBRATION
        test_chunk.custom = CustomTestUtil.good_wcs_with_function()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.CALIBRATION
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.CALIBRATION
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        product_type = chunk.ProductType.CALIBRATION
        expected_ctype = "FARADAY"
        with pytest.raises(ValueError) as ex:
            wcs_util.CustomAxisUtil.compute_bounds(
                artifacts, product_type, expected_ctype)
        assert ('CTYPE must be the same across all Artifacts' in str(ex.value))
        # user_chunk = True, range is not None
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.CALIBRATION
        test_chunk.custom = CustomTestUtil.good_wcs_with_range()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.CALIBRATION
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.CALIBRATION
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        product_type = chunk.ProductType.CALIBRATION
        expected_ctype = "RM"
        actual_interval = wcs_util.CustomAxisUtil.compute_bounds(
            artifacts, product_type, expected_ctype)
        expected_interval = Interval(1.1, 11.1)
        self.assertEqual(expected_interval.lower, actual_interval.lower)
        self.assertEqual(expected_interval.upper, actual_interval.upper)
        # user_chunk = True, get_num_pixels: bounds with 3 samples that
        # traverses _merge_into_list completely
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.CALIBRATION
        test_chunk.custom = CustomTestUtil.good_wcs_with_bounds_3_samples()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.CALIBRATION
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.CALIBRATION
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        product_type = chunk.ProductType.CALIBRATION
        expected_ctype = "RM"
        actual_interval = wcs_util.CustomAxisUtil.compute_bounds(
            artifacts, product_type, expected_ctype)
        expected_interval = Interval(-1.2, 11.2)
        self.assertEqual(expected_interval.lower, actual_interval.lower)
        self.assertEqual(expected_interval.upper, actual_interval.upper)
        # user_chunk = True, function is not None
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.CALIBRATION
        test_chunk.custom = CustomTestUtil.good_wcs_with_function()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.CALIBRATION
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.CALIBRATION
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        product_type = chunk.ProductType.CALIBRATION
        expected_ctype = "RM"
        actual_interval = wcs_util.CustomAxisUtil.compute_bounds(
            artifacts, product_type, expected_ctype)
        expected_interval = Interval(-49.5, 19950.5)
        self.assertEqual(expected_interval.lower, actual_interval.lower)
        self.assertEqual(expected_interval.upper, actual_interval.upper)

    def test_compute(self):
        # _choose_product returns Artifact.product (SCIENCE),
        # user_chunk = False
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.CALIBRATION
        test_chunk.custom = CustomTestUtil.good_wcs()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.CALIBRATION
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.SCIENCE
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        actual_axis = wcs_util.CustomAxisUtil.compute(artifacts)
        expected_axis = None
        self.assertEqual(expected_axis, actual_axis)
        # _choose_product returns Artifact.product (CALIBRATION),
        # user_chunk = False
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.SCIENCE
        test_chunk.custom = CustomTestUtil.good_wcs()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.SCIENCE
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.CALIBRATION
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        actual_axis = wcs_util.CustomAxisUtil.compute(artifacts)
        expected_axis = None
        self.assertEqual(expected_axis, actual_axis)
        # _choose_product returns Part.product (SCIENCE),
        # user_chunk = False
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.SCIENCE
        test_chunk.custom = CustomTestUtil.good_wcs()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.SCIENCE
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.PREVIEW
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        actual_axis = wcs_util.CustomAxisUtil.compute(artifacts)
        expected_axis = None
        self.assertEqual(expected_axis, actual_axis)
        # _choose_product returns Part.product (CALIBRATION),
        # user_chunk = False
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.SCIENCE
        test_chunk.custom = CustomTestUtil.good_wcs()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.CALIBRATION
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.PREVIEW
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        actual_axis = wcs_util.CustomAxisUtil.compute(artifacts)
        expected_axis = None
        self.assertEqual(expected_axis, actual_axis)
        # _choose_product returns Chunk.product (SCIENCE), user_chunk = False
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.SCIENCE
        test_chunk.custom = CustomTestUtil.good_wcs()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.PREVIEW
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.PREVIEW
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        actual_axis = wcs_util.CustomAxisUtil.compute(artifacts)
        expected_axis = None
        self.assertEqual(expected_axis, actual_axis)
        # _choose_product returns Chunk.product (CALIBRATION),
        # user_chunk = False
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.CALIBRATION
        test_chunk.custom = CustomTestUtil.good_wcs()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.PREVIEW
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.PREVIEW
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        actual_axis = wcs_util.CustomAxisUtil.compute(artifacts)
        expected_axis = None
        self.assertEqual(expected_axis, actual_axis)
        # _choose_product returns None, user_chunk = False
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.PREVIEW
        test_chunk.custom = CustomTestUtil.good_wcs()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.PREVIEW
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.PREVIEW
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        actual_axis = wcs_util.CustomAxisUtil.compute(artifacts)
        expected_axis = None
        self.assertEqual(expected_axis, actual_axis)
        # _choose_product returns Artifact.product (SCIENCE),
        # user_chunk = True, Chunk.custom is None
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.SCIENCE
        test_chunk.custom = None
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.SCIENCE
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.SCIENCE
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        expected_axis = None
        actual_axis = wcs_util.CustomAxisUtil.compute(artifacts)
        self.assertEqual(expected_axis, actual_axis)
        # _choose_product returns Artifact.product (SCIENCE),
        # user_chunk = True, Chunk.custom is not None
        # bad Chunk.custom.axis.axis.ctype
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.SCIENCE
        test_chunk.custom = CustomTestUtil.bad_ctype_wcs()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.SCIENCE
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.SCIENCE
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        with pytest.raises(ValueError) as ex:
            actual_axis = wcs_util.CustomAxisUtil.compute(artifacts)
        assert ('Unsupported CTYPE:' in str(ex.value))
        # _choose_product returns Artifact.product (SCIENCE),
        # user_chunk = True, Chunk.custom is not None
        # first_ctype == Chunk.custom.axis.axis.ctype
        test_chunk = Chunk()
        test_chunk.product_type = chunk.ProductType.SCIENCE
        test_chunk.custom = CustomTestUtil.good_wcs_with_function()
        test_chunks = TypedList(Chunk, test_chunk)
        part_name = "test_part"
        part_product_type = chunk.ProductType.SCIENCE
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.SCIENCE
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        expected_ctype = "RM"
        expected_sample = Interval(-49.5, 19950.5)
        expected_samples = [expected_sample]
        expected_bounds = Interval(-49.5, 19950.5, expected_samples)
        expected_dimension = 200
        expected_axis = plane.CustomAxis(expected_ctype, expected_bounds,
                                         expected_dimension)
        actual_axis = wcs_util.CustomAxisUtil.compute(artifacts)
        self.assertEqual(expected_axis.ctype, actual_axis.ctype)
        self.assertEqual(expected_axis.bounds.lower, actual_axis.bounds.lower)
        self.assertEqual(expected_axis.bounds.upper, actual_axis.bounds.upper)
        self.assertEqual(expected_axis.bounds.samples[0].lower,
                         actual_axis.bounds.samples[0].lower)
        self.assertEqual(expected_axis.bounds.samples[0].upper,
                         actual_axis.bounds.samples[0].upper)
        self.assertEqual(expected_axis.dimension, actual_axis.dimension)
        # _choose_product returns Artifact.product (SCIENCE),
        # user_chunk = True, Chunk.custom is not None
        # first_ctype == Chunk.custom.axis.axis.ctype
        test_chunk_1 = Chunk()
        test_chunk_1.product_type = chunk.ProductType.SCIENCE
        test_chunk_1.custom = CustomTestUtil.good_wcs_with_function()
        test_chunk_2 = Chunk()
        test_chunk_2.product_type = chunk.ProductType.SCIENCE
        test_chunk_2.custom = CustomTestUtil.good_wcs_with_function()
        test_chunk_2.custom.axis.axis.ctype = "FARADAY"
        test_chunks = TypedList(Chunk, test_chunk_1, test_chunk_2)
        part_name = "test_part"
        part_product_type = chunk.ProductType.SCIENCE
        part = Part(part_name, part_product_type, test_chunks)
        uri = 'mast:HST/product/test_file.jpg'
        artifact_product_type = chunk.ProductType.SCIENCE
        release_type = ReleaseType.DATA
        artifact = Artifact(uri, artifact_product_type, release_type)
        artifact.parts = TypedOrderedDict((Part), (part_name, part))
        artifacts = TypedList(Artifact, artifact)
        expected_ctype = "RM"
        expected_sample = Interval(-49.5, 19950.5)
        expected_samples = [expected_sample]
        expected_bounds = Interval(-49.5, 19950.5, expected_samples)
        expected_dimension = 200
        expected_axis = plane.CustomAxis(expected_ctype, expected_bounds,
                                         expected_dimension)
        with pytest.raises(ValueError) as ex:
            actual_axis = wcs_util.CustomAxisUtil.compute(artifacts)
        assert ('CTYPE must be the same across all Artifacts' in str(ex.value))


# Supporting Classes for generating test data
class CustomTestUtil:
    def __init__(self):
        pass

    @staticmethod
    def good_wcs():
        ctype = "RM"
        unit = "rad/m^2"
        px = float(0.5)
        sx = float(54321.0)
        nx = 200
        ds = float(0.01)
        return CustomTestUtil.get_test_function(ctype, unit, px, sx, nx, ds)

    @staticmethod
    def good_wcs_with_negative_delta():
        ctype = "RM"
        unit = "rad/m^2"
        px = float(0.5)
        sx = float(54321.0)
        nx = 200
        ds = float(-0.02)
        return CustomTestUtil.get_test_function(ctype, unit, px, sx, nx, ds)

    @staticmethod
    def good_wcs_with_range():
        ctype = "RM"
        unit = "rad/m^2"
        px = float(0.5)
        sx = float(54321.0)
        nx = 200
        ds = float(0.01)
        return CustomTestUtil.get_test_function_with_range(
            ctype, unit, px, sx, nx, ds)

    @staticmethod
    def good_wcs_with_bounds_3_samples():
        ctype = "RM"
        unit = "rad/m^2"
        px = float(0.5)
        sx = float(54321.0)
        nx = 200
        ds = float(0.01)
        return CustomTestUtil.get_test_function_with_bounds_3_samples(
            ctype, unit, px, sx, nx, ds)

    @staticmethod
    def good_wcs_with_function():
        ctype = "RM"
        unit = "rad/m^2"
        px = float(0.5)
        sx = float(1.0)
        nx = 200
        ds = float(0.01)
        return CustomTestUtil.get_test_function_with_function(
            ctype, unit, px, sx, nx, ds)

    @staticmethod
    def bad_ctype_wcs():
        #  Should fail on the ctype
        ctype = "bla"
        unit = "rad/m^2"
        px = float(0.5)
        sx = float(54321.0)
        nx = 200
        ds = float(0.01)
        return CustomTestUtil.get_test_function(ctype, unit, px, sx, nx, ds)

    @staticmethod
    def bad_cunit_wcs():
        ctype = "RM"
        unit = "foo"
        px = float(0.5)
        sx = float(54321.0)
        nx = 200
        ds = float(0.01)
        return CustomTestUtil.get_test_function(ctype, unit, px, sx, nx, ds)

    @staticmethod
    def bad_delta():
        axis_1d = CoordAxis1D(wcs.Axis("RM", "rad/m**2"))
        # delta < 0.0 is bad
        ref_coord = RefCoord(float(1.0), float(2.0))
        axis_1d.function = CoordFunction1D(
            int(100), -0.01, ref_coord)

        return chunk.CustomWCS(axis_1d)

    @staticmethod
    def bad_range_wcs():
        px = float(0.5)
        sx = float(54321.0)
        nx = 200
        ds = float(0.01)
        axis_1d = wcs.CoordAxis1D(wcs.Axis("RM", "rad/m**2"))

        # divide into 2 samples with a gap between
        c1 = RefCoord(px, sx)
        c2 = RefCoord(0.0, 0.0)
        c3 = RefCoord(px + nx * 0.66, sx + nx * ds * 0.66)
        c4 = RefCoord(px + nx, sx + nx * ds)
        axis_1d.bounds = CoordBounds1D()
        axis_1d.bounds.samples.append(CoordRange1D(c1, c3))
        axis_1d.bounds.samples.append(CoordRange1D(c4, c2))

        return chunk.CustomWCS(axis_1d)

    @staticmethod
    def get_test_function(ctype, unit, px, sx, nx, ds):
        axis_1d = CoordAxis1D(wcs.Axis(ctype, unit))
        ref_coord = RefCoord(px, sx)
        axis_1d.function = CoordFunction1D(nx, ds, ref_coord)
        custom_wcs = chunk.CustomWCS(axis_1d)
        return custom_wcs

    @staticmethod
    def get_test_function_with_range(ctype, unit, px, sx, nx, ds):
        error = None
        start = RefCoord(float(0.9), float(1.1))
        end = RefCoord(float(10.9), float(11.1))
        range = CoordRange1D(start, end)
        axis_1d = CoordAxis1D(wcs.Axis(ctype, unit), error, range)
        ref_coord = RefCoord(px, sx)
        axis_1d.function = CoordFunction1D(nx, ds, ref_coord)
        custom_wcs = chunk.CustomWCS(axis_1d)
        return custom_wcs

    @staticmethod
    def get_test_function_with_bounds_3_samples(ctype, unit, px, sx, nx, ds):
        error = None
        range = None
        start = RefCoord(float(0.8), float(1.1))
        end = RefCoord(float(10.8), float(11.1))
        b_range_1 = CoordRange1D(start, end)
        start = RefCoord(float(0.9), float(1.2))
        end = RefCoord(float(10.9), float(11.2))
        b_range_2 = CoordRange1D(start, end)
        start = RefCoord(float(-0.9), float(-1.2))
        end = RefCoord(float(0.6), float(0.2))
        b_range_3 = CoordRange1D(start, end)
        samples = caom_util.TypedList(CoordRange1D, b_range_1, b_range_2,
                                      b_range_3)
        bounds = CoordBounds1D(samples)
        axis_1d = wcs.CoordAxis1D(wcs.Axis(ctype, unit), error, range, bounds)
        ref_coord = wcs.RefCoord(px, sx)
        axis_1d.function = wcs.CoordFunction1D(nx, ds, ref_coord)
        custom_wcs = chunk.CustomWCS(axis_1d)
        return custom_wcs

    @staticmethod
    def get_test_function_with_function(ctype, unit, px, sx, nx, ds):
        error = None
        range = None
        bounds = None
        naxis = int(1)
        delta = float(2.5)
        ref_coord = wcs.RefCoord(float(1.0), float(2.0))
        function = CoordFunction1D(naxis, delta, ref_coord)
        axis_1d = CoordAxis1D(wcs.Axis(ctype, unit), error, range, bounds,
                              function)
        ref_coord = RefCoord(px, sx)
        axis_1d.function = CoordFunction1D(nx, ds, ref_coord)
        custom_wcs = chunk.CustomWCS(axis_1d)
        return custom_wcs
