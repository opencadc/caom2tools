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

from caom2utils import wcsvalidator, validate_wcs, InvalidWCSError
from caom2 import artifact, observation, part, plane, caom_util, Axis, \
    chunk, CoordAxis1D, CoordBounds1D, CoordFunction1D, CoordRange1D, \
    PolarizationWCS, RefCoord, wcs
from caom2.caom_util import TypedList, TypedOrderedDict
from ..wcsvalidator import WcsPolarizationState
import pytest
import unittest

single_test = False


# TemporalWCS validator tests
@pytest.mark.skipif(single_test, reason='Single test mode')
class TemporalWCSValidatorTests(unittest.TestCase):
    def test_temporalwcs_validator(self):
        good_temporal_wcs = TimeTestUtil.good_wcs()
        assert(good_temporal_wcs.axis.function is not None)

        wcsvalidator._validate_temporal_wcs(good_temporal_wcs)
        wcsvalidator._validate_temporal_wcs(None)

    def test_bad_temporalwcs(self):
        bad_temporal_wcs = TimeTestUtil.bad_ctype_wcs()
        with self.assertRaisesRegex(
                InvalidWCSError, 'unexpected TIMESYS, CTYPE'):
            wcsvalidator._validate_temporal_wcs(bad_temporal_wcs)

        bad_temporal_wcs = TimeTestUtil.bad_cunit_wcs()
        with self.assertRaisesRegex(InvalidWCSError, 'unexpected CUNIT'):
            wcsvalidator._validate_temporal_wcs(bad_temporal_wcs)

        bad_temporal_wcs = TimeTestUtil.bad_range_wcs()
        with self.assertRaisesRegex(
                InvalidWCSError, 'range.end not >= range.start'):
            wcsvalidator._validate_temporal_wcs(bad_temporal_wcs)


# CustomWCS validator tests
@pytest.mark.skipif(single_test, reason='Single test mode')
class CustomWCSValidatorTests(unittest.TestCase):
    def test_customwcs_validator(self):
        good_custom_wcs = CustomTestUtil.good_wcs()
        assert(good_custom_wcs.axis is not None)

        wcsvalidator._validate_custom_wcs(good_custom_wcs)
        wcsvalidator._validate_custom_wcs(None)

    def test_bad_customwcs(self):
        bad_custom_wcs = CustomTestUtil.bad_ctype_wcs()
        with self.assertRaisesRegex(
                InvalidWCSError, 'CUSTOM_WCS_VALIDATION_ERROR:'):
            wcsvalidator._validate_custom_wcs(bad_custom_wcs)

        bad_custom_wcs = CustomTestUtil.bad_cunit_wcs()
        with self.assertRaisesRegex(
                InvalidWCSError, 'CUSTOM_WCS_VALIDATION_ERROR:'):
            wcsvalidator._validate_custom_wcs(bad_custom_wcs)

        bad_custom_wcs = CustomTestUtil.bad_range_wcs()
        with self.assertRaisesRegex(
                InvalidWCSError, 'CUSTOM_WCS_VALIDATION_ERROR:'):
            wcsvalidator._validate_custom_wcs(bad_custom_wcs)

        bad_custom_wcs = CustomTestUtil.bad_bounds_wcs()
        with self.assertRaisesRegex(
                InvalidWCSError, 'CUSTOM_WCS_VALIDATION_ERROR:'):
            wcsvalidator._validate_custom_wcs(bad_custom_wcs)

        bad_custom_wcs = CustomTestUtil.bad_function_wcs()
        with self.assertRaisesRegex(
                InvalidWCSError, 'CUSTOM_WCS_VALIDATION_ERROR:'):
            wcsvalidator._validate_custom_wcs(bad_custom_wcs)


@pytest.mark.skipif(single_test, reason='Single test mode')
class SpatialWCSValidatorTests(unittest.TestCase):
    def test_spatialwcs_validator(self):

        spatialtest = SpatialTestUtil()
        good_spatial_wcs = spatialtest.good_wcs()
        assert(good_spatial_wcs.axis.function is not None)
        wcsvalidator._validate_spatial_wcs(good_spatial_wcs)
        # None is valid
        wcsvalidator._validate_spatial_wcs(None)

    # Without the toPolygon function moved over from Java, there's not
    # alot of way to test an invalid spatial wcs.
    # def test_invalid_spatial_wcs(self):
    #     spatialtest = SpatialTestUtil()
    #     position = spatialtest.bad_wcs()
    #     with self.assertRaises(InvalidWCSError):
    #         WcsValidator.validate_spatial_wcs(position)


@pytest.mark.skipif(single_test, reason='Single test mode')
class SpectralWCSValidatorTests(unittest.TestCase):
    def test_spectralwcs_validator(self):
        energyTest = EnergyTestUtil()
        good_spectral_wcs = energyTest.good_wcs()
        assert(good_spectral_wcs.axis.function is not None)

        wcsvalidator._validate_spectral_wcs(good_spectral_wcs)
        wcsvalidator._validate_spectral_wcs(None)

    def test_invalid_spectral_wcs(self):
        energytest = EnergyTestUtil()
        energy_wcs = energytest.bad_wcs()
        with pytest.raises(InvalidWCSError):
            wcsvalidator._validate_spectral_wcs(energy_wcs)


# validate_wcs tests
@pytest.mark.skipif(single_test, reason='Single test mode')
class ValidateWCSTests(unittest.TestCase):
    def test_None(self):
        # CAOM2 entity is None
        validate_wcs(None)

    def test_observation(self):
        # CAOM2 entity is an Observation
        obs = ObservationTestUtil.get_test_observation()
        validate_wcs(obs)

    def test_plane(self):
        # CAOM2 entity is a Plane
        plane = PlaneTestUtil.get_test_plane('plane1')
        validate_wcs(plane)

    def test_part(self):
        # CAOM2 entity is a Part
        pname = "part1"
        product_type = chunk.ProductType.SCIENCE
        part = PartTestUtil.get_test_part(pname, product_type)
        validate_wcs(part)

    def test_artifact(self):
        # CAOM2 entity is an Artifact
        auri = "uri:foo/bar"
        product_type = chunk.ProductType.SCIENCE
        # with valid wcs
        artifact = ArtifactTestUtil.get_test_artifact(auri, product_type)
        validate_wcs(artifact)

    def test_artifact_with_null_wcs(self):
        # with null wcs
        auri = "uri:foo/bar"
        product_type = chunk.ProductType.SCIENCE
        artifact = ArtifactTestUtil.get_test_artifact(auri, product_type)
        validate_wcs(artifact)
        c = artifact.parts['test_part'].chunks[0]
        c.position_axis_1 = 1
        c.position_axis_2 = 2
        c.energy_axis = 3
        c.time_axis = 4
        c.custom_axis = 5

        c.naxis = 6
        with pytest.raises(InvalidWCSError):
            validate_wcs(c)

        # Not probably reasonable Chunks, but should still be valid
        # Different combinations of this will be represented in
        # different data sets
        c.position = None
        c.position_axis_1 = None
        c.position_axis_2 = None
        c.energy_axis = 1
        c.time_axis = 2
        c.custom_axis = 3
        c.naxis = 3
        validate_wcs(artifact)

        c.position = SpatialTestUtil.good_wcs()
        c.energy = None
        c.position_axis_1 = 1
        c.position_axis_2 = 2
        c.energy_axis = None
        c.time_axis = 3
        c.custom_axis = 4
        c.naxis = 4
        validate_wcs(artifact)

        c.energy = EnergyTestUtil.good_wcs()
        c.time = None
        c.energy_axis = 3
        c.time_axis = None
        validate_wcs(artifact)

        c.time = TimeTestUtil.good_wcs()
        c.time_axis = 5
        c.naxis = 5
        c.polarization = None
        validate_wcs(artifact)

        c.energy = None
        c.energy_axis = None
        c.naxis = 2
        validate_wcs(artifact)

        c.time = None
        c.time_axis = None
        validate_wcs(artifact)

        # Assert: all WCS should be null at this step
        c.position = None
        c.position_axis_1 = None
        c.position_axis_2 = None
        c.naxis = None
        validate_wcs(artifact)

    def test_failure(self):
        test_object = type('', (), {})()
        with pytest.raises(InvalidWCSError):
            validate_wcs(test_object)


# Supporting Classes for generating test data
class TimeTestUtil:
    def __init__(self):
        pass

    @staticmethod
    def good_wcs():
        px = float(0.5)
        sx = float(54321.0)
        nx = 200
        ds = float(0.01)
        goodwcs = TimeTestUtil.get_test_function(True, px, sx*nx*ds, nx, ds)
        return goodwcs

    @staticmethod
    def bad_ctype_wcs():
        px = float(0.5)
        sx = float(54321.0)
        nx = 200
        ds = float(0.01)

        badwcs = TimeTestUtil.get_test_function(True, px, sx*nx*ds, nx, ds)
        #  Should fail on the ctype
        badwcs.axis.axis.ctype = "bla"
        return badwcs

    @staticmethod
    def bad_cunit_wcs():
        badcunit = TimeTestUtil.good_wcs()
        badcunit.axis.axis.cunit = "foo"
        return badcunit

    @staticmethod
    def bad_range_wcs():
        px = float(0.5)
        sx = float(54321.0)
        nx = 200
        ds = float(0.01)
        axis_1d = wcs.CoordAxis1D(wcs.Axis("UTC", "d"))
        temporal_wcs = chunk.TemporalWCS(axis_1d)
        temporal_wcs.exposure = 300.0
        temporal_wcs.resolution = 0.1

        # divide into 2 samples with a gap between
        c1 = wcs.RefCoord(px, sx)
        c2 = wcs.RefCoord(0.0, 0.0)
        c3 = wcs.RefCoord(px + nx * 0.66, sx + nx * ds * 0.66)
        c4 = wcs.RefCoord(px + nx, sx + nx * ds)
        temporal_wcs.axis.bounds = wcs.CoordBounds1D()
        temporal_wcs.axis.bounds.samples.append(wcs.CoordRange1D(c1, c3))
        temporal_wcs.axis.bounds.samples.append(wcs.CoordRange1D(c4, c2))

        return temporal_wcs

    @staticmethod
    def get_test_function(complete, px, sx, nx, ds):
        axis_1d = wcs.CoordAxis1D(wcs.Axis("UTC", "d"))

        if complete:
            wcs.exposure = 300.0
            wcs.resolution = 0.1

        temporal_wcs = chunk.TemporalWCS(axis_1d)
        ref_coord = wcs.RefCoord(px, sx)
        temporal_wcs.axis.function = wcs.CoordFunction1D(nx, ds, ref_coord)
        return temporal_wcs


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
        sx = float(1.0)
        nx = 200
        ds = float(0.01)
        return CustomTestUtil.get_test_function(ctype, unit, px, sx, nx, ds)

    @staticmethod
    def bad_range_wcs():
        ctype = "RM"
        unit = "rad/m^2"
        error = None
        start = RefCoord(float(0.9), float(1.1))
        end = RefCoord(float(10.9), float(1.1))
        range = CoordRange1D(start, end)
        axis_1d = CoordAxis1D(wcs.Axis(ctype, unit), error, range)
        return chunk.CustomWCS(axis_1d)

    @staticmethod
    def bad_bounds_wcs():
        ctype = "RM"
        unit = "rad/m^2"
        error = None
        range = None
        c1 = RefCoord(float(0.9), float(1.1))
        c2 = RefCoord(float(10.9), float(1.1))
        bounds = CoordBounds1D()
        bounds.samples.append(CoordRange1D(c1, c2))
        axis_1d = CoordAxis1D(wcs.Axis(ctype, unit), error, range, bounds)
        return chunk.CustomWCS(axis_1d)

    @staticmethod
    def bad_function_wcs():
        ctype = "RM"
        unit = "rad/m^2"
        error = None
        range = None
        c1 = RefCoord(float(0.9), float(1.1))
        c2 = RefCoord(float(10.9), float(1.1))
        bounds = CoordBounds1D()
        bounds.samples.append(CoordRange1D(c1, c2))
        naxis = 1
        delta = 0.0
        ref_coord = RefCoord(float(0.9), float(1.1))
        func = CoordFunction1D(naxis, delta, ref_coord)
        axis_1d = CoordAxis1D(wcs.Axis(ctype, unit), error, range, bounds,
                              func)
        return chunk.CustomWCS(axis_1d)

    @staticmethod
    def get_test_function(ctype, unit, px, sx, nx, ds):
        axis_1d = wcs.CoordAxis1D(wcs.Axis(ctype, unit))
        ref_coord = wcs.RefCoord(px, sx)
        axis_1d.function = wcs.CoordFunction1D(nx, ds, ref_coord)
        custom_wcs = chunk.CustomWCS(axis_1d)
        return custom_wcs


class SpatialTestUtil:
    def __init__(self):
        pass

    @staticmethod
    def good_wcs():
        px = float(0.5)
        py = float(0.5)
        sx = float(20.0)
        sy = float(10.0)
        # not used in java code although declared?
        # double dp = 1000.0;
        # double ds = 1.0;
        return SpatialTestUtil.get_test_function(px, py, sx, sy, False)

    # This is in the java code, but without the toPolygon function in python
    # is possibly irrelevant here. With the basic validator this will produce
    # a good WCS value.
    @staticmethod
    def bad_wcs():
        axis1 = wcs.Axis("RA---TAN", "deg")
        axis2 = wcs.Axis("DEC--TAN", "deg")
        axis = wcs.CoordAxis2D(axis1, axis2)
        spatial_wcs = chunk.SpatialWCS(axis)
        spatial_wcs.equinox = None
        dim = wcs.Dimension2D(1024, 1024)
        ref = wcs.Coord2D(wcs.RefCoord(
            512.0, 10.0),  wcs.RefCoord(512.0, 20.0))
        #  Create Invalid function
        axis.function = wcs.CoordFunction2D(
            dim, ref, 1.0e-3, 0.0, 0.0, 0.0)  # singular CD matrix
        return spatial_wcs

    @staticmethod
    def get_test_function(px, py, sx, sy, gal):
        axis1 = wcs.Axis("RA", "deg")
        axis2 = wcs.Axis("DEC", "deg")

        if gal:
            axis1 = wcs.Axis("GLON", "deg")
            axis2 = wcs.Axis("GLAT", "deg")

        axis_2d = wcs.CoordAxis2D(axis1, axis2)

        spatial_wcs = chunk.SpatialWCS(axis_2d)
        spatial_wcs.coordsys = "ICRS"

        wcs.coordsys = "ICRS"
        if gal:
            spatial_wcs.coordsys = None

        wcs.equinox = None

        # Simple frame set: 1000x1000 pixels, 1 pixel = 1.0e-3 deg
        dim = wcs.Dimension2D(1000, 1000)
        ref = wcs.Coord2D(wcs.RefCoord(px, sx), wcs.RefCoord(py, sy))
        axis_2d.function = wcs.CoordFunction2D(
            dim, ref, 1.e-3, 0.0, 0.0, 1.0e-3)
        return spatial_wcs


#  TODO: add functions to create test data here
class PolarizationTestUtil:
    def __init__(self):
        pass

    @staticmethod
    def good_wcs(self):
        return None

    @staticmethod
    def bad_wcs(self):
        return None


# There's got to be a better way to do this, but of the 2 ways I tried,
# none is really great. Leaving these here for now.
BANDPASS_NAME = "H-Alpha-narrow"
TRANSITION = wcs.EnergyTransition("H", "alpha")


class EnergyTestUtil:
    def __init__(self):
        pass

    @staticmethod
    def good_wcs():
        px = float(0.5)
        sx = float(400.0)
        nx = float(200.0)
        ds = float(1.0)
        # SpectralWCS
        energy = EnergyTestUtil.getTestRange(True, px, sx * nx * ds, nx, ds)

        c1 = wcs.RefCoord(0.5, 2000.0)
        energy.axis.function = wcs.CoordFunction1D(100, 10.0, c1)
        return energy

    @staticmethod
    def bad_wcs():
        px = float(0.5)
        sx = float(400.0)
        nx = 200
        ds = float(1.0)
        bad_energy = EnergyTestUtil.getTestFunction(
            True, px, sx * nx * ds, nx, ds)
        # Make function invalid
        c1 = wcs.RefCoord(0.5, 2000.0)
        bad_energy.axis.function = wcs.CoordFunction1D(100, 10.0, c1)
        return wcs

    @staticmethod
    def getTestRange(complete, px, sx, nx, ds):
        axis = wcs.CoordAxis1D(wcs.Axis("WAVE", "nm"))
        spectral_wcs = chunk.SpectralWCS(axis, "TOPOCENT")
        if complete:
            spectral_wcs.bandpassName = BANDPASS_NAME
            spectral_wcs.restwav = 6563.0e-10  # meters
            spectral_wcs.resolvingPower = 33000.0
            spectral_wcs.transition = TRANSITION

        c1 = wcs.RefCoord(px, sx)
        c2 = wcs.RefCoord(px + nx, sx + nx * ds)
        spectral_wcs.axis.range = wcs.CoordRange1D(c1, c2)
        # log.debug("test range: " + axis.range)
        return spectral_wcs

    @staticmethod
    def getTestFunction(complete, px, sx, nx, ds):
        axis = wcs.CoordAxis1D(wcs.Axis("WAVE", "nm"))
        spectral_wcs = chunk.SpectralWCS(axis, "TOPOCENT")
        if complete:
            spectral_wcs.bandpassName = BANDPASS_NAME
            spectral_wcs.restwav = 6563.0e-10  # meters
            spectral_wcs.resolvingPower = 33000.0
            spectral_wcs.transition = TRANSITION

        c1 = wcs.RefCoord(px, sx)
        spectral_wcs.axis.function = wcs.CoordFunction1D(nx, ds, c1)
        return spectral_wcs


class ObservationTestUtil():
    def __init__(self):
        pass

    @staticmethod
    def get_test_observation():
        obs = observation.Observation(
           "collection", "obsID", observation.Algorithm("algo"))
        p1 = PlaneTestUtil.get_test_plane("planeID1")
        p2 = PlaneTestUtil.get_test_plane("planeID2")
        p3 = PlaneTestUtil.get_test_plane("planeID3")
        obs.planes["planeID1"] = p1
        obs.planes["planeID2"] = p2
        obs.planes["planeID3"] = p3
        return obs


class PlaneTestUtil():
    def __init__(self):
        pass

    @staticmethod
    def get_test_plane(planeID):
        test_plane = plane.Plane(planeID)
        uri1 = 'uri:foo/bar1'
        uri2 = 'uri:foo/bar2'
        uri3 = 'uri:foo/bar3'
        product_type = chunk.ProductType.SCIENCE
        a1 = ArtifactTestUtil.get_test_artifact(uri1, product_type)
        a2 = ArtifactTestUtil.get_test_artifact(uri2, product_type)
        a3 = ArtifactTestUtil.get_test_artifact(uri3, product_type)
        test_plane.artifacts[uri1] = a1
        test_plane.artifacts[uri2] = a2
        test_plane.artifacts[uri3] = a3
        return test_plane


class ArtifactTestUtil():
    def __init__(self):
        pass

    @staticmethod
    def get_good_test_chunk(ptype):
        test_chunk = chunk.Chunk()
        test_chunk.position = SpatialTestUtil.good_wcs()
        test_chunk.energy = EnergyTestUtil.good_wcs()
        test_chunk.time = TimeTestUtil.good_wcs()
        # test_chunk.polarization = PolarizationTestUtil.good_wcs()
        test_chunk.custom_axis = 1
        test_chunk.custom = CustomTestUtil.good_wcs()

        return test_chunk

    @staticmethod
    def get_test_artifact(uri, ptype):
        # chunk.ProductType.SCIENCE is a common type
        if ptype is None:
            ptype = chunk.ProductType.SCIENCE
        test_artifact = artifact.Artifact(
            uri, ptype, artifact.ReleaseType.DATA)
        chunks = TypedList(chunk.Chunk)
        chunks.append(ArtifactTestUtil.get_good_test_chunk(ptype))

        pname = "test_part"
        test_part = PartTestUtil.get_test_part(pname, ptype)
        test_artifact.parts = TypedOrderedDict(part.Part)
        test_artifact.parts[pname] = test_part
        return test_artifact


class PartTestUtil():
    def __init__(self):
        pass

    @staticmethod
    def get_test_part(pname, ptype):
        # chunk.ProductType.SCIENCE is a common type
        chunks = TypedList(chunk.Chunk)
        chunks.append(ArtifactTestUtil.get_good_test_chunk(ptype))

        test_part = part.Part(pname, ptype, chunks)
        return test_part


@pytest.mark.skipif(single_test, reason='Single test mode')
class TestValidatePolarizationWcs():
    def test_none_polarization_wcs(self):
        # Polarization is None, should not produce an error
        wcsvalidator._validate_polarization_wcs(None)

    def test_range(self):
        # Polarization range is None, should not produce an error
        axis = Axis("STOKES", "cunit")
        axis_1d = CoordAxis1D(axis)
        polarization = PolarizationWCS(axis_1d)
        wcsvalidator._validate_polarization_wcs(polarization)

        # Polarization axis range contains valid positive values
        start = RefCoord(float(0.9), float(1.1))
        end = RefCoord(float(9.9), float(10.1))
        p_range = CoordRange1D(start, end)
        axis_1d.range = p_range
        polarization = PolarizationWCS(axis_1d)
        wcsvalidator._validate_polarization_wcs(polarization)

        # Polarization axis range contains valid negative values
        start = RefCoord(float(-8.1), float(-7.9))
        end = RefCoord(float(-1.1), float(-0.9))
        n_range = CoordRange1D(start, end)
        axis_1d.range = n_range
        polarization = PolarizationWCS(axis_1d)
        wcsvalidator._validate_polarization_wcs(polarization)

        # Polarization axis range contains invalid positive values
        start = RefCoord(float(0.9), float(1.1))
        end = RefCoord(float(10.9), float(11.1))
        p_range = CoordRange1D(start, end)
        axis_1d.range = p_range
        polarization = PolarizationWCS(axis_1d)
        with pytest.raises(InvalidWCSError) as ex:
            wcsvalidator._validate_polarization_wcs(polarization)
        assert('Invalid Polarization WCS' in str(ex.value))
        assert('11' in str(ex.value))

        # Polarization axis range contains invalid negative values
        start = RefCoord(float(-9.1), float(-8.9))
        end = RefCoord(float(-1.1), float(-0.9))
        n_range = CoordRange1D(start, end)
        axis_1d.range = n_range
        polarization = PolarizationWCS(axis_1d)
        with pytest.raises(InvalidWCSError) as ex:
            wcsvalidator._validate_polarization_wcs(polarization)
        assert('Invalid Polarization WCS' in str(ex.value))
        assert('-9' in str(ex.value))

        # Polarization axis range contains an invalid value (0) within a range
        start = RefCoord(float(-8.1), float(-7.9))
        end = RefCoord(float(9.9), float(10.1))
        range = CoordRange1D(start, end)
        axis_1d.range = range
        polarization = PolarizationWCS(axis_1d)
        with pytest.raises(InvalidWCSError) as ex:
            wcsvalidator._validate_polarization_wcs(polarization)
        assert('Invalid Polarization WCS' in str(ex.value))
        assert('0' in str(ex.value))

    def test_bounds(self):
        # Polarization bounds is None, should not produce an error
        axis = Axis("STOKES", "cunit")
        axis_1d = CoordAxis1D(axis)
        polarization = PolarizationWCS(axis_1d)
        wcsvalidator._validate_polarization_wcs(polarization)

        # Polarization axis bounds contains one valid range
        start = RefCoord(float(0.9), float(1.1))
        end = RefCoord(float(9.9), float(10.1))
        p_range = CoordRange1D(start, end)
        samples = caom_util.TypedList(CoordRange1D, p_range)
        axis_1d.bounds = CoordBounds1D(samples)
        polarization = PolarizationWCS(axis_1d)
        wcsvalidator._validate_polarization_wcs(polarization)

        # Polarization axis bounds contains more than one valid range
        start = RefCoord(float(0.9), float(1.1))
        end = RefCoord(float(9.9), float(10.1))
        p_range = CoordRange1D(start, end)
        start = RefCoord(float(-8.1), float(-7.9))
        end = RefCoord(float(-1.1), float(-0.9))
        n_range = CoordRange1D(start, end)
        samples = caom_util.TypedList(CoordRange1D, p_range, n_range)
        axis_1d.bounds = CoordBounds1D(samples)
        polarization = PolarizationWCS(axis_1d)
        wcsvalidator._validate_polarization_wcs(polarization)

        # Polarization axis bounds contains one invalid range
        start = RefCoord(float(0.9), float(1.1))
        end = RefCoord(float(10.9), float(11.1))
        p_range = CoordRange1D(start, end)
        samples = caom_util.TypedList(CoordRange1D, p_range)
        axis_1d.bounds = CoordBounds1D(samples)
        polarization = PolarizationWCS(axis_1d)
        with pytest.raises(InvalidWCSError) as ex:
            wcsvalidator._validate_polarization_wcs(polarization)
        assert('Invalid Polarization WCS' in str(ex.value))
        assert('11' in str(ex.value))

        # Polarization axis bounds contains more than one invalid range
        start = RefCoord(float(0.9), float(1.1))
        end = RefCoord(float(9.9), float(10.1))
        p_range = CoordRange1D(start, end)
        start = RefCoord(float(-9.1), float(-8.9))
        end = RefCoord(float(-1.1), float(-0.9))
        n_range = CoordRange1D(start, end)
        samples = caom_util.TypedList(CoordRange1D, p_range, n_range)
        axis_1d.bounds = CoordBounds1D(samples)
        polarization = PolarizationWCS(axis_1d)
        with pytest.raises(InvalidWCSError) as ex:
            wcsvalidator._validate_polarization_wcs(polarization)
        assert('Invalid Polarization WCS' in str(ex.value))
        assert('-9' in str(ex.value))

    def test_function(self):
        # Polarization function is None, should not produce an error
        axis = Axis("STOKES", "cunit")
        axis_1d = CoordAxis1D(axis)
        polarization = PolarizationWCS(axis_1d)
        wcsvalidator._validate_polarization_wcs(polarization)

        # Polarization axis function with naxis=1
        naxis = int(1)
        delta = float(2.5)
        ref_coord = wcs.RefCoord(float(1.0), float(2.0))
        axis_1d.function = CoordFunction1D(naxis, delta, ref_coord)
        polarization = PolarizationWCS(axis_1d)
        wcsvalidator._validate_polarization_wcs(polarization)

        # Polarization axis function with naxis>1
        naxis = int(3)
        delta = float(2.5)
        ref_coord = wcs.RefCoord(float(1.0), float(2.0))
        axis_1d.function = CoordFunction1D(naxis, delta, ref_coord)
        polarization = PolarizationWCS(axis_1d)
        wcsvalidator._validate_polarization_wcs(polarization)

        # Polarization axis function with invalid naxis=0
        naxis = int(0)
        delta = float(2.5)
        ref_coord = wcs.RefCoord(float(1.0), float(2.0))
        axis_1d.function = CoordFunction1D(naxis, delta, ref_coord)
        polarization = PolarizationWCS(axis_1d)
        with pytest.raises(InvalidWCSError) as ex:
            wcsvalidator._validate_polarization_wcs(polarization)
        assert('Invalid Polarization WCS' in str(ex.value))
        assert('Invalid naxis value' in str(ex.value))


@pytest.mark.skipif(single_test, reason='Single test mode')
class TestWcsPolarizationState():
    def test_all(self):
        # valid keys
        for i in range(1, 11):
            WcsPolarizationState.to_value(i)
        for i in range(-8, 0):
            WcsPolarizationState.to_value(i)
        # invalid keys
        with pytest.raises(KeyError):
            WcsPolarizationState.to_value(int(-9))
        with pytest.raises(KeyError):
            WcsPolarizationState.to_value(int(0))
        with pytest.raises(KeyError):
            WcsPolarizationState.to_value(int(11))
