from caom2 import SimpleObservation, CompositeObservation, Circle, Polygon
from caom2 import CoordCircle2D, CoordPolygon2D

import unittest


class ObsCompare(unittest.TestCase):

    def __init__(self, _ignored):
        self._testMethodName = 'compare_observations'
        super(ObsCompare, self).__init__()

    def runTest(self):
        pass

    def compare_observations(self, expected, actual, version, compare_ids=False):

        assert ((isinstance(expected, SimpleObservation) and
                 isinstance(actual, SimpleObservation)) or
                (isinstance(expected, CompositeObservation) and
                 isinstance(actual, CompositeObservation))), \
            ("Observation types do not match 0 vs 1".
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
        if compare_ids:
            self.assertEqual(expected._id, actual._id)

        self.compare_entity_attributes(expected, actual)

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
            self.compare_requirements(expected.requirements,
                                      actual.requirements)

        self.compare_planes(expected.planes, actual.planes, version, compare_ids)

        if (isinstance(expected, CompositeObservation) and
                isinstance(actual, CompositeObservation)):
            self.compare_members(expected.members, actual.members)

    def compare_proposal(self, expected, actual):
        if expected is None and actual is None:
            return
        self.assertIsNotNone(expected)
        self.assertIsNotNone(actual)
        self.assertEqual(expected.id, actual.id)
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

    def compare_planes(self, expected, actual, version, compare_ids=False):
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
            self.assertEqual(expected_plane.product_id,
                             actual_plane.product_id)
            self.assertIsNotNone(expected_plane._id)
            self.assertIsNotNone(actual_plane._id)
            if compare_ids:
                self.assertEqual(expected_plane._id, actual_plane._id)

            if version >= 23:
                self.assertEqual(
                    expected_plane.creator_id, actual_plane.creator_id,
                    "creator_id")

            self.compare_entity_attributes(expected_plane, actual_plane)

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

            if version >= 22:
                self.compare_position(expected_plane.position,
                                      actual_plane.position)
                self.compare_energy(expected_plane.energy, actual_plane.energy)
                print('comparing time')
                self.compare_time(expected_plane.time, actual_plane.time)
                print('compared time')
                self.compare_polarization(expected_plane.polarization,
                                          actual_plane.polarization)

            self.compare_artifacts(expected_plane.artifacts,
                                   actual_plane.artifacts, version)

    def compare_position(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual, "position")
        else:
            self.compare_shape(expected.bounds, actual.bounds)
            self.compare_dimension2d(expected.dimension, actual.dimension)
            self.assertEqual(expected.resolution, actual.resolution,
                             "resolution")
            self.assertEqual(expected.sample_size, actual.sample_size,
                             "sample_size")
            self.assertEqual(expected.time_dependent, actual.time_dependent,
                             "time_dependent")

    def compare_energy(self, expected, actual):
        print("comparing energy")
        if expected is None:
            self.assertIsNone(actual, "energy")
        else:
            self.compare_interval(expected.bounds, actual.bounds)
            self.assertEqual(expected.dimension, actual.dimension, "dimension")
            self.assertEqual(expected.resolving_power, actual.resolving_power,
                             "resolving_power")
            self.assertEqual(expected.sample_size, actual.sample_size,
                             "sample_size")
            self.assertEqual(expected.bandpass_name, actual.bandpass_name,
                             "bandpass_name")
            self.assertEqual(expected.em_band, actual.em_band, "em_band")
            self.assertEqual(expected.bandpass_name, actual.bandpass_name,
                             "bandpass_name")
            self.compare_wcs_energy_transition(expected.transition,
                                               actual.transition)

    def compare_time(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual, "time")
        else:
            self.compare_interval(expected.bounds, actual.bounds)
            self.assertEqual(expected.dimension, actual.dimension, "dimension")
            self.assertEqual(expected.resolution, actual.resolution,
                             "resolution")
            self.assertEqual(expected.sample_size, actual.sample_size,
                             "sample_size")
            self.assertEqual(expected.exposure, actual.exposure, "exposure")

    def compare_polarization(self, expected, actual):
        if expected is None:
            self.assertIsNone(expected, "polarization")
        else:
            self.assertEqual(expected.dimension, actual.dimension, "dimension")
            if expected.polarization_states is None:
                self.assertIsNone(actual.polarization_states,
                                  "polarization_states")
            else:
                self.assertEqual(len(expected.polarization_states),
                                 len(actual.polarization_states),
                                 "different number of polarization_states")
                for index, state in enumerate(expected.polarization_states):
                    self.assertEqual(state, actual.polarization_states[index],
                                     "polarization_state")

    def compare_shape(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual, "shape")
        else:
            if isinstance(expected, Polygon):
                self.assertIsNotNone(actual, "shape is None")
                self.assertTrue(isinstance(actual, Polygon),
                                "mismatched shapes" +
                                actual.__class__.__name__)
                expected_points = expected.points
                actual_points = actual.points
                self.assertEqual(len(expected_points), len(actual_points),
                                 "different number of points")
                for index, point in enumerate(expected_points):
                    self.compare_point(point, actual_points[index])
                actual_samples = actual.samples
                self.assertIsNotNone(actual_samples, "shape is None")
                # TODO - figure out the import for MultiPolygon
                # self.assertTrue(isinstance(actual_samples, MultiPolygon),
                #                 "mismatched shapes" +
                #                 actual.__class__.__name__)
                expected_samples = expected.samples
                expected_vertices = expected_samples.vertices
                actual_vertices = actual_samples.vertices
                self.assertEqual(len(expected_vertices), len(actual_vertices),
                                 "different number of vertices")
                for index, vertex in enumerate(expected_vertices):
                    self.compare_vertices(vertex, actual_vertices[index])
            elif isinstance(expected, Circle):
                self.assertIsNotNone(actual, "shape is None")
                self.assertTrue(isinstance(actual, Circle),
                                "mismatched shapes" +
                                actual.__class__.__name__)
                self.assertEqual(expected.center.cval1, actual.center.cval1,
                                 "mismatched centers (cval1)")
                self.assertEqual(expected.center.cval2, actual.center.cval2,
                                 "mismatched centers (cval2)")
                self.assertEqual(expected.radius, actual.radius,
                                 "mismatched radius")
            else:
                raise TypeError("Unsupported shape type "
                                + expected.__class__.__name__)

    def compare_interval(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual, "interval")
        else:
            self.assertEqual(expected.lower, actual.lower, "lower")
            self.assertEqual(expected.upper, actual.upper, "upper")
            if expected.samples is None:
                self.assertIsNone(actual.samples, "samples")
            else:
                self.assertEqual(len(actual.samples), len(expected.samples),
                                 "samples")
                for index, sample in enumerate(expected.samples):
                    self.compare_sub_interval(sample, actual.samples[index])

    def compare_sub_interval(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual, "sub_interval")
        else:
            self.assertEqual(expected.lower, actual.lower, "lower")
            self.assertEqual(expected.upper, actual.upper, "upper")

    def compare_wcs_energy_transition(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual, "wcs_energy_transition")
        else:
            self.assertEqual(expected.species, actual.species, "species")
            self.assertEqual(expected.transition, actual.transition,
                             "transition")

    def compare_vertices(self, expected, actual):
        self.assertEqual(expected.cval1, actual.cval1)
        self.assertEqual(expected.cval2, actual.cval2)
        self.assertEqual(expected.type, actual.type)

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
        self.assertEqual(expected.background_std_dev,
                         actual.background_std_dev)
        self.assertEqual(expected.flux_density_limit,
                         actual.flux_density_limit)
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
        for expected_plane_uri in expected:
            self.assertTrue(expected_plane_uri in actual)

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

            self.compare_entity_attributes(expected_artifact, actual_artifact)

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

            self.compare_entity_attributes(expected_part, actual_part)

            self.assertEqual(expected_part.name, actual_part.name)
            self.assertEqual(expected_part.product_type,
                             actual_part.product_type)
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

            self.compare_entity_attributes(expected_chunk, actual_chunk)

            self.assertEqual(expected_chunk.product_type,
                             actual_chunk.product_type)
            self.assertEqual(expected_chunk.naxis, actual_chunk.naxis)
            self.assertEqual(expected_chunk.observable_axis,
                             actual_chunk.observable_axis)
            self.assertEqual(expected_chunk.position_axis_1,
                             actual_chunk.position_axis_1)
            self.assertEqual(expected_chunk.position_axis_2,
                             actual_chunk.position_axis_2)
            self.assertEqual(expected_chunk.energy_axis,
                             actual_chunk.energy_axis)
            self.assertEqual(expected_chunk.time_axis, actual_chunk.time_axis)
            self.assertEqual(expected_chunk.polarization_axis,
                             actual_chunk.polarization_axis)
            self.compare_observable_axis(expected_chunk.observable,
                                         actual_chunk.observable)
            self.compare_spatial_wcs(expected_chunk.position,
                                     actual_chunk.position)
            self.compare_spectral_wcs(expected_chunk.energy,
                                      actual_chunk.energy)
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
        for expected_range, actual_range in zip(expected.samples,
                                                actual.samples):
            self.compare_coord_range1d(expected_range, actual_range)

    def compare_coord_bounds2d(self, expected, actual):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        if (isinstance(expected, CoordCircle2D) and
                isinstance(actual, CoordCircle2D)):
            self.compare_coord_circle2d(expected, actual)
        elif (isinstance(expected, CoordPolygon2D) and
              isinstance(actual, CoordPolygon2D)):
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

    def compare_entity_attributes(self, expected, actual):
        if expected.last_modified is not None and \
                actual.last_modified is not None:
            self.assertEqual(expected.last_modified, actual.last_modified,
                             "last_modified")

        if expected.max_last_modified is not None and \
                actual.max_last_modified is not None:
            self.assertEqual(expected.max_last_modified,
                             actual.max_last_modified, "max_last_modified")

        if expected.meta_checksum is not None and \
                actual.meta_checksum is not None:
            self.assertEqual(expected.meta_checksum, actual.meta_checksum,
                             "meta_checksum")

        if expected.acc_meta_checksum is not None and \
                actual.acc_meta_checksum is not None:
            self.assertEqual(expected.acc_meta_checksum,
                             actual.acc_meta_checksum, "acc_meta_checksum")

