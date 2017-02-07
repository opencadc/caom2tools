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

""" Defines TestObservationReaderWriter class """

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from six import StringIO
import os
import unittest

from lxml import etree
from .xml_compare import xml_compare

from . import caom_test_instances
from .. import obs_reader_writer
from .. import observation
from .. import wcs

THIS_DIR = os.path.dirname(os.path.realpath(__file__))


def minimal_simple(depth, bounds_is_circle, version):
    instances = caom_test_instances.Caom2TestInstances()
    instances.complete = False
    instances.depth = depth
    instances.bounds_is_circle = bounds_is_circle
    instances.caom_version = version
    return instances.get_simple_observation()


def complete_simple(depth, bounds_is_circle, version):
    instances = caom_test_instances.Caom2TestInstances()
    instances.complete = True
    instances.depth = depth
    instances.bounds_is_circle = bounds_is_circle
    instances.caom_version = version
    return instances.get_simple_observation()


def minimal_composite(depth, bounds_is_circle, version):
    instances = caom_test_instances.Caom2TestInstances()
    instances.complete = False
    instances.depth = depth
    instances.bounds_is_circle = bounds_is_circle
    instances.caom_version = version
    return instances.get_composite_observation()


def complete_composite(depth, bounds_is_circle, version):
    instances = caom_test_instances.Caom2TestInstances()
    instances.complete = True
    instances.depth = depth
    instances.bounds_is_circle = bounds_is_circle
    instances.caom_version = version
    return instances.get_composite_observation()


class TestObservationReaderWriter(unittest.TestCase):

    def test_invalid_long_id(self):
        simple_observation = minimal_simple(1, False, 20)
        writer = obs_reader_writer.ObservationWriter(
            False, False, "caom2", obs_reader_writer.CAOM20_NAMESPACE)
        output = StringIO()
        writer.write(simple_observation, output)
        xml = output.getvalue()
        output.close()
        xml = xml.replace("caom2:id=\"", "caom2:id=\"x")
        f = open('/tmp/test.xml', 'w')
        f.write(xml)
        f.close()
        reader = obs_reader_writer.ObservationReader(False)
        try:
            reader.read('/tmp/test.xml')
            self.fail("invalid long id should throw ValueError")
        except ValueError:
            pass

    def test_invalid_uuid(self):
        simple_observation = minimal_simple(1, False, 21)
        writer = obs_reader_writer.ObservationWriter(False, False)  # default writer is 2.1
        output = StringIO()
        writer.write(simple_observation, output)
        xml = output.getvalue()
        output.close()
        xml = xml.replace("0000", "xxxx", 1)
        f = open('/tmp/test.xml', 'w')
        f.write(xml)
        f.close()
        reader = obs_reader_writer.ObservationReader(False)
        try:
            reader.read('/tmp/test.xml')
            self.fail("invalid uuid id should throw ValueError")
        except ValueError:
            pass

    def test_minimal_simple(self):
        for version in (20, 21):
            for i in range(1, 6):
                print("Test Minimal Simple {} version {}".format(i, version))
                # CoordBounds2D as CoordCircle2D
                simple_observation = minimal_simple(i, True, version)
                # write empty elements
                self.observation_test(simple_observation, True, True, version)
                # do not write empty elements
                self.observation_test(simple_observation, True, False, version)
                # CoordBounds2D as CoordPolygon2D
                simple_observation = minimal_simple(i, False, version)
                # write empty elements
                self.observation_test(simple_observation, True, True, version)
                # do not write empty elements
                self.observation_test(simple_observation, True, False, version)

    def test_complete_simple(self):
        for version in (20, 21, 22):
            for i in range(1, 6):
                print("Test Complete Simple {} version {}".format(i, version))
                # CoordBounds2D as CoordCircle2D
                simple_observation = complete_simple(i, True, version)
                # write empty elements
                self.observation_test(simple_observation, True, True, version)
                # do not write empty elements
                self.observation_test(simple_observation, True, False, version)
                # CoordBounds2D as CoordPolygon2D
                simple_observation = complete_simple(i, False, version)
                # write empty elements
                self.observation_test(simple_observation, True, True, version)
                # do not write empty elements
                self.observation_test(simple_observation, True, False, version)

    def test_minimal_composite(self):
        for version in (20, 21, 22):
            for i in range(1, 6):
                print("Test Minimal Composite {} version {}".format(i, version))
                # CoordBounds2D as CoordCircle2D
                composite_observation = minimal_composite(i, True, version)
                # write empty elements
                self.observation_test(composite_observation, True, True, version)
                # do not write empty elements
                self.observation_test(composite_observation, True, False, version)
                # CoordBounds2D as CoordPolygon2D
                composite_observation = minimal_composite(i, False, version)
                # write empty elements
                self.observation_test(composite_observation, True, True, version)
                # do not write empty elements
                self.observation_test(composite_observation, True, False, version)

    def test_complete_composite(self):
        for version in (20, 21, 22):
            for i in range(1, 6):
                print("Test Complete Composite {} version {}".format(i, version))
                # CoordBounds2D as CoordCircle2D
                composite_observation = complete_composite(i, True, version)
                # write empty elements
                self.observation_test(composite_observation, True, True, version)
                # do not write empty elements
                self.observation_test(composite_observation, True, False, version)
                # CoordBounds2D as CoordPolygon2D
                composite_observation = complete_composite(i, False, version)
                # write empty elements
                self.observation_test(composite_observation, True, True, version)
                # do not write empty elements
                self.observation_test(composite_observation, True, False, version)

    def test_versions(self):
        composite_observation = complete_composite(6, True, 20)
        self.observation_test(composite_observation, True, True, 20)
        self.observation_test(composite_observation, True, True, 21)
        self.observation_test(composite_observation, True, True, 22)

        composite_observation = complete_composite(6, True, 21)
        self.observation_test(composite_observation, True, True, 20)
        self.observation_test(composite_observation, True, True, 21)
        self.observation_test(composite_observation, True, True, 22)

        composite_observation = complete_composite(6, True, 22)
        self.observation_test(composite_observation, True, True, 20)
        self.observation_test(composite_observation, True, True, 21)
        self.observation_test(composite_observation, True, True, 22)

    def observation_test(self, obs, validate, write_empty_collections, version):
        if version == 20:
            writer = obs_reader_writer.ObservationWriter(
                validate, write_empty_collections, "caom2",
                obs_reader_writer.CAOM20_NAMESPACE)
        elif version == 21:
            writer = obs_reader_writer.ObservationWriter(
                validate, write_empty_collections, "caom2",
                obs_reader_writer.CAOM21_NAMESPACE)
        else:
            writer = obs_reader_writer.ObservationWriter(
                validate, write_empty_collections)
        xml_file = open('/tmp/test.xml', 'w')
        writer.write(obs, xml_file)
        xml_file.close()
        reader = obs_reader_writer.ObservationReader(True)
        returned = reader.read('/tmp/test.xml')
        self.compare_observations(obs, returned, version)

    def compare_observations(self, expected, actual, version):

        assert ((isinstance(expected, observation.SimpleObservation) and
                isinstance(actual, observation.SimpleObservation)) or
                (isinstance(expected, observation.CompositeObservation) and
                isinstance(actual, observation.CompositeObservation))), (
                    "Observation types do not match 0 vs 1".
                    format(expected.__class__.__name__,
                           actual.__class__.__name__))

        self.assertIsNotNone(expected.collection)
        self.assertIsNotNone(actual.collection)
        self.assertEqual(expected.collection, actual.collection)

        self.assertIsNotNone(expected.observation_id)
        self.assertIsNotNone(actual.observation_id)
        self.assertEqual(expected.observation_id, actual.observation_id)

        self.assertIsNotNone(expected._id)
        self.assertIsNotNone(actual._id)
        self.assertEqual(expected._id, actual._id)

        self.assertIsNotNone(expected._last_modified)
        self.assertIsNotNone(actual._last_modified)
        self.assertEqual(expected._last_modified, actual._last_modified)

        self.assertIsNotNone(expected.algorithm)
        self.assertIsNotNone(actual.algorithm)
        self.assertEqual(expected.algorithm.name, actual.algorithm.name)

        self.assertEqual(expected.sequence_number, actual.sequence_number)
        self.assertEqual(expected.intent, actual.intent)
        self.assertEqual(expected.meta_release, actual.meta_release)
        self.compare_proposal(expected.proposal, actual.proposal)
        self.compare_target(expected.target, actual.target)
        self.compare_target_position(expected.target_position,
                                     actual.target_position)
        self.compare_telescope(expected.telescope, actual.telescope)
        self.compare_instrument(expected.instrument, actual.instrument)
        self.compare_environment(expected.environment, actual.environment)
        if version == 21:
            self.compare_requirements(expected.requirements, actual.requirements)

        self.compare_planes(expected.planes, actual.planes, version)

        if (isinstance(expected, observation.CompositeObservation) and
                isinstance(actual, observation.CompositeObservation)):
            self.compare_members(expected.members, actual.members)

    def compare_proposal(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.assertEqual(expected.proposal_id, actual.proposal_id)
        self.assertEqual(expected.pi_name, actual.pi_name)
        self.assertEqual(expected.project, actual.project)
        self.assertEqual(expected.title, actual.title)
        self.assertEqual(len(expected.keywords), len(actual.keywords))
        for keyword in expected.keywords:
            self.assertTrue(keyword in actual.keywords)

    def compare_target(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.assertEqual(expected.name, actual.name)
        self.assertEqual(expected.target_type, actual.target_type)
        self.assertEqual(expected.redshift, actual.redshift)
        for keyword in expected.keywords:
            self.assertTrue(keyword in actual.keywords)

    def compare_target_position(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.assertIsNotNone(actual.coordinates)
        self.assertIsNotNone(actual.coordsys)
        self.compare_point(expected.coordinates, actual.coordinates)
        self.assertEqual(expected.coordsys, actual.coordsys)
        self.assertEqual(expected.equinox, actual.equinox)

    def compare_telescope(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.assertEqual(expected.name, actual.name)
        self.assertEqual(expected.geo_location_x, actual.geo_location_x)
        self.assertEqual(expected.geo_location_y, actual.geo_location_y)
        self.assertEqual(expected.geo_location_z, actual.geo_location_z)
        for keyword in expected.keywords:
            self.assertTrue(keyword in actual.keywords)

    def compare_instrument(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.assertEqual(expected.name, actual.name)
        for keyword in expected.keywords:
            self.assertTrue(keyword in actual.keywords)

    def compare_environment(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.assertEqual(expected.seeing, actual.seeing)
        self.assertEqual(expected.humidity, actual.humidity)
        self.assertEqual(expected.elevation, actual.elevation)
        self.assertEqual(expected.tau, actual.tau)
        self.assertEqual(expected.wavelength_tau, actual.wavelength_tau)
        self.assertEqual(expected.ambient_temp, actual.ambient_temp)
        self.assertEqual(expected.photometric, actual.photometric)

    def compare_members(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.assertEqual(len(expected), len(actual))
        for expected_member, actual_member in zip(expected, actual):
            self.compare_observation_uri(expected_member, actual_member)

    def compare_observation_uri(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.assertEquals(expected.uri, actual.uri)
        self.assertEquals(expected.collection, actual.collection)
        self.assertEquals(expected.observation_id, actual.observation_id)

    def compare_requirements(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.assertEquals(expected.flag, actual.flag)

    def compare_planes(self, expected, actual, version):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.assertEqual(len(expected), len(actual))
        for key in expected:
            self.assertTrue(key in actual)
            expected_plane = expected[key]
            actual_plane = actual[key]
            self.assertIsNotNone(expected_plane)
            self.assertIsNotNone(actual_plane)
            self.assertEqual(expected_plane.product_id, actual_plane.product_id)
            self.assertIsNotNone(expected_plane._id)
            self.assertIsNotNone(actual_plane._id)
            self.assertEqual(expected_plane._id, actual_plane._id)
            self.assertIsNotNone(expected_plane._last_modified)
            self.assertIsNotNone(actual_plane._last_modified)
            self.assertEqual(expected_plane._last_modified,
                             actual_plane._last_modified)
            self.assertEqual(expected_plane.meta_release,
                             actual_plane.meta_release)
            self.assertEqual(expected_plane.data_release,
                             actual_plane.data_release)
            self.assertEqual(expected_plane.data_product_type,
                             actual_plane.data_product_type)
            self.assertEqual(expected_plane.calibration_level,
                             actual_plane.calibration_level)
            self.compare_provenance(expected_plane.provenance,
                                    actual_plane.provenance)
            self.compare_metrics(expected_plane.metrics, actual_plane.metrics)
            if version == 21:
                self.compare_quality(expected_plane.quality,
                                     actual_plane.quality)

            self.compare_artifacts(expected_plane.artifacts,
                                   actual_plane.artifacts, version)

    def compare_provenance(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.assertEqual(expected.version, actual.version)
        self.assertEqual(expected.project, actual.project)
        self.assertEqual(expected.producer, actual.producer)
        self.assertEqual(expected.run_id, actual.run_id)
        self.assertEqual(expected.reference, actual.reference)
        self.assertEqual(expected.last_executed, actual.last_executed)
        self.compare_inputs(expected.inputs, actual.inputs)

    def compare_metrics(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.assertEqual(expected.source_number_density,
                         actual.source_number_density)
        self.assertEqual(expected.background, actual.background)
        self.assertEqual(expected.background_std_dev, actual.background_std_dev)
        self.assertEqual(expected.flux_density_limit, actual.flux_density_limit)
        self.assertEqual(expected.mag_limit, actual.mag_limit)

    def compare_quality(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.assertEqual(expected.flag, actual.flag)

    def compare_inputs(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.assertEqual(len(expected), len(actual))
        for expected_plane_uri, actual_plane_uri in zip(expected, actual):
            self.assertEqual(expected_plane_uri, actual_plane_uri)

    def compare_artifacts(self, expected, actual, version):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.assertEqual(len(expected), len(actual))
        for key in expected:
            self.assertTrue(key in actual)
            expected_artifact = expected[key]
            actual_artifact = actual[key]
            self.assertIsNotNone(expected_artifact)
            self.assertIsNotNone(actual_artifact)
            self.assertIsNotNone(expected_artifact._id)
            self.assertIsNotNone(actual_artifact._id)
            self.assertEqual(expected_artifact._id, actual_artifact._id)
            self.assertIsNotNone(expected_artifact._last_modified)
            self.assertIsNotNone(actual_artifact._last_modified)
            self.assertEqual(expected_artifact._last_modified,
                             actual_artifact._last_modified)
            self.assertEqual(expected_artifact.uri, actual_artifact.uri)
            self.assertEqual(expected_artifact.content_type,
                             actual_artifact.content_type)
            self.assertEqual(expected_artifact.content_length,
                             actual_artifact.content_length)
            self.assertEqual(expected_artifact.product_type,
                             actual_artifact.product_type)
            if version > 21:
                self.assertEqual(expected_artifact.release_type,
                                 actual_artifact.release_type)
            self.compare_parts(expected_artifact.parts,
                               actual_artifact.parts, version)

    def compare_parts(self, expected, actual, version):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.assertEqual(len(expected), len(actual))
        for key in expected:
            self.assertTrue(key in actual)
            expected_part = expected[key]
            actual_part = actual[key]
            self.assertIsNotNone(expected_part)
            self.assertIsNotNone(actual_part)
            self.assertIsNotNone(expected_part._id)
            self.assertIsNotNone(actual_part._id)
            self.assertEqual(expected_part._id, actual_part._id)
            self.assertIsNotNone(expected_part._last_modified)
            self.assertIsNotNone(actual_part._last_modified)
            self.assertEqual(expected_part._last_modified,
                             actual_part._last_modified)
            self.assertEqual(expected_part.name, actual_part.name)
            self.assertEqual(expected_part.product_type, actual_part.product_type)
            self.compare_chunks(expected_part.chunks, actual_part.chunks)

    def compare_chunks(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.assertEqual(len(expected), len(actual))
        for expected_chunk, actual_chunk in zip(expected, actual):
            self.assertIsNotNone(expected_chunk)
            self.assertIsNotNone(actual_chunk)
            self.assertIsNotNone(expected_chunk._id)
            self.assertIsNotNone(actual_chunk._id)
            self.assertEqual(expected_chunk._id, actual_chunk._id)
            self.assertIsNotNone(expected_chunk._last_modified)
            self.assertIsNotNone(actual_chunk._last_modified)
            self.assertEqual(expected_chunk._last_modified,
                             actual_chunk._last_modified)
            self.assertEqual(expected_chunk.product_type,
                             actual_chunk.product_type)
            self.assertEqual(expected_chunk.naxis, actual_chunk.naxis)
            self.assertEqual(expected_chunk.observable_axis,
                             actual_chunk.observable_axis)
            self.assertEqual(expected_chunk.position_axis_1,
                             actual_chunk.position_axis_1)
            self.assertEqual(expected_chunk.position_axis_2,
                             actual_chunk.position_axis_2)
            self.assertEqual(expected_chunk.energy_axis, actual_chunk.energy_axis)
            self.assertEqual(expected_chunk.time_axis, actual_chunk.time_axis)
            self.assertEqual(expected_chunk.polarization_axis,
                             actual_chunk.polarization_axis)
            self.compare_observable_axis(expected_chunk.observable,
                                         actual_chunk.observable)
            self.compare_spatial_wcs(expected_chunk.position,
                                     actual_chunk.position)
            self.compare_spectral_wcs(expected_chunk.energy, actual_chunk.energy)
            self.compare_temporal_wcs(expected_chunk.time, actual_chunk.time)
            self.compare_polarization_wcs(expected_chunk.polarization,
                                          actual_chunk.polarization)

    def compare_observable_axis(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.compare_slice(expected.dependent, actual.dependent)
        self.compare_slice(expected.independent, actual.independent)

    def compare_spatial_wcs(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.compare_coord_axis2d(expected.axis, actual.axis)
        self.assertEqual(expected.coordsys, actual.coordsys)
        self.assertEqual(expected.equinox, actual.equinox)
        self.assertEqual(expected.resolution, actual.resolution)

    def compare_spectral_wcs(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.compare_coord_axis1d(expected.axis, actual.axis)
        self.assertEqual(expected.bandpass_name, actual.bandpass_name)
        self.assertEqual(expected.resolving_power, actual.resolving_power)
        self.assertEqual(expected.restfrq, actual.restfrq)
        self.assertEqual(expected.restwav, actual.restwav)
        self.assertEqual(expected.specsys, actual.specsys)
        self.assertEqual(expected.ssysobs, actual.ssysobs)
        self.assertEqual(expected.ssyssrc, actual.ssyssrc)
        self.assertEqual(expected.velang, actual.velang)
        self.assertEqual(expected.velosys, actual.velosys)
        self.assertEqual(expected.zsource, actual.zsource)

    def compare_temporal_wcs(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.compare_coord_axis1d(expected.axis, actual.axis)
        self.assertEqual(expected.exposure, actual.exposure)
        self.assertEqual(expected.resolution, actual.resolution)
        self.assertEqual(expected.timesys, actual.timesys)
        self.assertEqual(expected.trefpos, actual.trefpos)
        self.assertEqual(expected.mjdref, actual.mjdref)

    def compare_polarization_wcs(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.compare_coord_axis1d(expected.axis, actual.axis)

    def compare_axis(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        self.assertIsNotNone(actual.ctype)
        self.assertIsNotNone(actual.cunit)
        self.assertEqual(expected.ctype, actual.ctype)
        self.assertEqual(expected.cunit, actual.cunit)

    def compare_coord2d(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        self.compare_ref_coord(expected.coord1, actual.coord1)
        self.compare_ref_coord(expected.coord2, actual.coord2)

    def compare_value_coord2d(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        self.assertEqual(expected.coord1, actual.coord1)
        self.assertEqual(expected.coord2, actual.coord2)

    def compare_coord_axis1d(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        self.compare_coord_error(expected.error, actual.error)
        self.compare_coord_range1d(expected.range, actual.range)
        self.compare_coord_bounds1d(expected.bounds, actual.bounds)
        self.compare_coord_function1d(expected.function, actual.function)

    def compare_coord_axis2d(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        self.assertIsNotNone(actual.axis1)
        self.assertIsNotNone(actual.axis2)
        self.compare_axis(expected.axis1, actual.axis1)
        self.compare_axis(expected.axis2, actual.axis2)
        self.compare_coord_error(expected.error1, actual.error1)
        self.compare_coord_error(expected.error2, actual.error2)
        self.compare_coord_range2d(expected.range, actual.range)
        self.compare_coord_bounds2d(expected.bounds, actual.bounds)
        self.compare_coord_function2d(expected.function, actual.function)

    def compare_coord_bounds1d(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        self.assertIsNotNone(expected.samples)
        self.assertIsNotNone(actual.samples)
        self.assertEqual(len(expected.samples), len(actual.samples))
        for expected_range, actual_range in zip(expected.samples, actual.samples):
            self.compare_coord_range1d(expected_range, actual_range)

    def compare_coord_bounds2d(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        if (isinstance(expected, wcs.CoordCircle2D) and
                isinstance(actual, wcs.CoordCircle2D)):
            self.compare_coord_circle2d(expected, actual)
        elif (isinstance(expected, wcs.CoordPolygon2D) and
              isinstance(actual, wcs.CoordPolygon2D)):
            self.compare_coord_polygon2d(expected, actual)
        else:
            self.fail("CoordBounds2D expected and actual are different types.")

    def compare_coord_circle2d(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        self.assertIsNotNone(actual.center)
        self.assertIsNotNone(actual.radius)
        self.compare_value_coord2d(expected.center, actual.center)
        self.assertEqual(expected.radius, actual.radius)

    def compare_coord_error(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual)
            return

        self.assertIsNotNone(actual)
        if expected.syser:
            self.assertIsNotNone(actual.syser)
            self.assertEqual(expected.syser, actual.syser)
        if expected.rnder:
            self.assertIsNotNone(actual.rnder)
            self.assertEqual(expected.rnder, actual.rnder)

    def compare_coord_function1d(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        self.assertEqual(expected.naxis, actual.naxis)
        self.assertEqual(expected.delta, actual.delta)
        self.compare_ref_coord(expected.ref_coord, actual.ref_coord)

    def compare_coord_function2d(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        self.assertIsNotNone(actual.dimension)
        self.assertIsNotNone(actual.ref_coord)
        self.assertIsNotNone(actual.cd11)
        self.assertIsNotNone(actual.cd12)
        self.assertIsNotNone(actual.cd21)
        self.assertIsNotNone(actual.cd22)
        self.compare_dimension2d(expected.dimension, actual.dimension)
        self.compare_coord2d(expected.ref_coord, actual.ref_coord)
        self.assertEqual(expected.cd11, actual.cd11, 0.0)
        self.assertEqual(expected.cd12, actual.cd12, 0.0)
        self.assertEqual(expected.cd21, actual.cd21, 0.0)
        self.assertEqual(expected.cd22, actual.cd22, 0.0)

    def compare_coord_polygon2d(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        self.assertIsNotNone(expected.vertices)
        self.assertIsNotNone(actual.vertices)
        self.assertEqual(len(expected.vertices), len(actual.vertices))
        for expected_coord_2d, actual_coord_2d in zip(expected.vertices,
                                                      actual.vertices):
            self.compare_value_coord2d(expected_coord_2d, actual_coord_2d)

    def compare_coord_range1d(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        self.compare_ref_coord(expected.start, actual.start)
        self.compare_ref_coord(expected.end, actual.end)

    def compare_coord_range2d(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        self.assertIsNotNone(actual.start)
        self.assertIsNotNone(actual.end)
        self.compare_coord2d(expected.start, actual.start)
        self.compare_coord2d(expected.end, actual.end)

    def compare_dimension2d(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        self.assertEqual(expected.naxis1, actual.naxis1)
        self.assertEqual(expected.naxis2, actual.naxis2)

    def compare_ref_coord(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        self.assertEqual(expected.pix, actual.pix)
        self.assertEqual(expected.val, actual.val)

    def compare_slice(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        self.assertIsNotNone(actual.bin)
        self.assertIsNotNone(actual.axis)
        self.assertEqual(expected.bin, actual.bin)
        self.compare_axis(expected.axis, actual.axis)

    def compare_point(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        self.assertIsNotNone(actual.cval1)
        self.assertIsNotNone(actual.cval2)
        self.assertEqual(expected.cval1, actual.cval1)
        self.assertEqual(expected.cval2, actual.cval2)


class TestRoundTrip(unittest.TestCase):

    TEST_DATA = 'data'

    def init(self):
        pass

    def get_file_list(self):
        return [f for f in os.listdir(THIS_DIR + "/" + TestRoundTrip.TEST_DATA)
                if f.endswith('.xml')]

    def do_test(self, reader, writer, filename):
        source_file_path = os.path.join(THIS_DIR + "/" +
                                        TestRoundTrip.TEST_DATA + "/" + filename)
        source_xml_fp = open(source_file_path, 'r')
        obs = reader.read(source_file_path)
        source_xml_fp.close()
        dest_file = StringIO()
        writer.write(obs, dest_file)

        source_dom = etree.parse(source_file_path).getroot()
        dest_dom = etree.fromstring(dest_file.getvalue())
        self.assertTrue(xml_compare(source_dom, dest_dom, reporter=print),
                        'files are different')

    # This test reads each file in XML_FILE_SOURCE_DIR, creates the CAOM2
    # objects and writes a file in XML_FILE_DEST_DIR based on the CAOM2
    # objects. The two XML files are then compared to ensure that they
    # are the same. The test fails if the files are not the same. The
    # test/data/*.xml files can be used in this test.

    def test_round_trip(self):
        print("Test Round Trip")

        try:
            self.init()
            files = self.get_file_list()
            self.assertTrue(len(files) > 0, 'No XML files in test data directory')

            reader = obs_reader_writer.ObservationReader(True)
            writer20 = obs_reader_writer.ObservationWriter(
                True, False, "caom2", obs_reader_writer.CAOM20_NAMESPACE)
            writer21 = obs_reader_writer.ObservationWriter(
                True, False, "caom2", obs_reader_writer.CAOM21_NAMESPACE)
            writer22 = obs_reader_writer.ObservationWriter(True, False, "caom2")
            for filename in files:
                if filename.endswith("CAOM-2.2.xml"):
                    self.do_test(reader, writer22, filename)
                elif filename.endswith("CAOM-2.1.xml"):
                    self.do_test(reader, writer21, filename)
                else:
                    self.do_test(reader, writer20, filename)

        except Exception:
            raise
