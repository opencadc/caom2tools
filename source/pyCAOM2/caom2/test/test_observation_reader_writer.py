#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
#***********************************************************************
#******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
#*************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2010.                            (c) 2010.
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
#***********************************************************************
#

""" Defines TestObservationReaderWriter class """

from caom2.caom2_simple_observation import SimpleObservation
from caom2.caom2_composite_observation import CompositeObservation
from caom2.wcs.caom2_coord_circle2d import CoordCircle2D
from caom2.wcs.caom2_coord_polygon2d import CoordPolygon2D
from caom2.xml.caom2_observation_writer import ObservationWriter
from caom2.xml.caom2_observation_reader import ObservationReader
from caom2testinstances import Caom2TestInstances
import os
import sys
import unittest

# put build at the start of the search path
sys.path.insert(0, os.path.abspath('../../lib.local/lib'))


class TestObservationReaderWriter(unittest.TestCase):

    def test_minimal_simple(self):

        for i in range(1, 6):
            print "Test Minimal Simple ", i
            # CoordBounds2D as CoordCircle2D
            observation = minimal_simple(i, True)
            # write empty elements
            test_observation(self, observation, True, True)
            # do not write empty elements
            test_observation(self, observation, True, False)
            # CoordBounds2D as CoordPolygon2D
            observation = minimal_simple(i, False)
            # write empty elements
            test_observation(self, observation, True, True)
            # do not write empty elements
            test_observation(self, observation, True, False)

    def test_complete_simple(self):

        for i in range(1, 6):
            print "Test Complete Simple ", i
            # CoordBounds2D as CoordCircle2D
            observation = complete_simple(i, True)
            # write empty elements
            test_observation(self, observation, True, True)
            # do not write empty elements
            test_observation(self, observation, True, False)
            # CoordBounds2D as CoordPolygon2D
            observation = complete_simple(i, False)
            # write empty elements
            test_observation(self, observation, True, True)
            # do not write empty elements
            test_observation(self, observation, True, False)

    def test_minimal_composite(self):

        for i in range(1, 6):
            print "Test Minimal Composite ", i
            # CoordBounds2D as CoordCircle2D
            observation = minimal_composite(i, True)
            # write empty elements
            test_observation(self, observation, True, True)
            # do not write empty elements
            test_observation(self, observation, True, False)
            # CoordBounds2D as CoordPolygon2D
            observation = minimal_composite(i, False)
            # write empty elements
            test_observation(self, observation, True, True)
            # do not write empty elements
            test_observation(self, observation, True, False)

    def test_complete_composite(self):

        for i in range(1, 6):
            print "Test Complete Composite ", i
            # CoordBounds2D as CoordCircle2D
            observation = complete_composite(i, True)
            # write empty elements
            test_observation(self, observation, True, True)
            # do not write empty elements
            test_observation(self, observation, True, False)
            # CoordBounds2D as CoordPolygon2D
            observation = complete_composite(i, False)
            # write empty elements
            test_observation(self, observation, True, True)
            # do not write empty elements
            test_observation(self, observation, True, False)


def minimal_simple(depth, bounds_is_circle):
    instances = Caom2TestInstances()
    instances.complete = False
    instances.depth = depth
    instances.bounds_is_circle = bounds_is_circle
    return instances.get_simple_observation()


def complete_simple(depth, bounds_is_circle):
    instances = Caom2TestInstances()
    instances.complete = True
    instances.depth = depth
    instances.bounds_is_circle = bounds_is_circle
    return instances.get_simple_observation()


def minimal_composite(depth, bounds_is_circle):
    instances = Caom2TestInstances()
    instances.complete = False
    instances.depth = depth
    instances.bounds_is_circle = bounds_is_circle
    return instances.get_composite_observation()


def complete_composite(depth, bounds_is_circle):
    instances = Caom2TestInstances()
    instances.complete = True
    instances.depth = depth
    instances.bounds_is_circle = bounds_is_circle
    return instances.get_composite_observation()


def test_observation(self, observation, validate, write_empty_collections):
    writer = ObservationWriter(validate, write_empty_collections)
    xmlfile = open('/tmp/test.xml', 'w')
    writer.write(observation, xmlfile)
    xmlfile.close()
    reader = ObservationReader(True)
    returned = reader.read('/tmp/test.xml')
    compareObservations(self, observation, returned)


def compareObservations(self, expected, actual):

    assert ((isinstance(expected, SimpleObservation) and
            isinstance(actual, SimpleObservation)) or
            (isinstance(expected, CompositeObservation) and
            isinstance(actual, CompositeObservation))), (
                "Observation types do not match 0 vs 1".
                format(expected.__class__.__name__, actual.__class__.__name__))

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
    compareProposal(self, expected.proposal, actual.proposal)
    compareTarget(self, expected.target, actual.target)
    compareTelescope(self, expected.telescope, actual.telescope)
    compareInstrument(self, expected.instrument, actual.instrument)

    comparePlanes(self, expected.planes, actual.planes)

    if (isinstance(expected, CompositeObservation) and
        isinstance(actual, CompositeObservation)):
        compareMembers(self, expected.members, actual.members)


def compareProposal(self, expected, actual):
    if (expected == None and actual == None):
        return
    self.assertIsNotNone(expected)
    self.assertIsNotNone(actual)
    self.assertEqual(expected.proposal_id, actual.proposal_id)
    self.assertEqual(expected.pi_name, actual.pi_name)
    self.assertEqual(expected.project, actual.project)
    self.assertEqual(expected.title, actual.title)
    self.assertEqual(len(expected.keywords), len(actual.keywords))
    for i in range(len(expected.keywords)):
        self.assertEquals(expected.keywords[i], actual.keywords[i])


def compareTarget(self, expected, actual):
    if (expected == None and actual == None):
        return
    self.assertIsNotNone(expected)
    self.assertIsNotNone(actual)
    self.assertEqual(expected.name, actual.name)
    self.assertEqual(expected.target_type, actual.target_type)
    self.assertEqual(expected.redshift, actual.redshift)
    for i in range(len(expected.keywords)):
        self.assertEquals(expected.keywords[i], actual.keywords[i])


def compareTargetPosition(self, expected, actual):
    if (expected == None and actual == None):
        return
    self.assertIsNotNone(expected)
    self.assertIsNotNone(actual)
    self.assertNotNone(actual.coordinates)
    self.assertNotNone(actual.coordsys)
    self.comparePoint(expected.coordinates, actual.coordinates)
    self.assertEqual(expected.coordsys, actual.coordsys)
    self.assertEqual(expected.equinox, actual.equinox)


def compareTelescope(self, expected, actual):
    if (expected == None and actual == None):
        return
    self.assertIsNotNone(expected)
    self.assertIsNotNone(actual)
    self.assertEqual(expected.name, actual.name)
    self.assertEqual(expected.geo_location_x, actual.geo_location_x)
    self.assertEqual(expected.geo_location_y, actual.geo_location_y)
    self.assertEqual(expected.geo_location_z, actual.geo_location_z)
    for i in range(len(expected.keywords)):
        self.assertEquals(expected.keywords[i], actual.keywords[i])


def compareInstrument(self, expected, actual):
    if (expected == None and actual == None):
        return
    self.assertIsNotNone(expected)
    self.assertIsNotNone(actual)
    self.assertEqual(expected.name, actual.name)
    for i in range(len(expected.keywords)):
        self.assertEquals(expected.keywords[i], actual.keywords[i])


def compareMembers(self, expected, actual):
    if (expected == None and actual == None):
        return
    self.assertIsNotNone(expected)
    self.assertIsNotNone(actual)
    self.assertEqual(len(expected), len(actual))
    for expected_member, actual_member in zip(expected, actual):
        compareObservationURI(self, expected_member, actual_member)


def compareObservationURI(self, expected, actual):
    if (expected == None and actual == None):
        return
    self.assertIsNotNone(expected)
    self.assertIsNotNone(actual)
    self.assertEquals(expected.uri, actual.uri)
    self.assertEquals(expected.collection, actual.collection)
    self.assertEquals(expected.observation_id, actual.observation_id)


def comparePlanes(self, expected, actual):
    if (expected == None and actual == None):
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
        compareProvenance(self, expected_plane.provenance,
                          actual_plane.provenance)
        compareArtifacts(self, expected_plane.artifacts,
                         actual_plane.artifacts)


def compareProvenance(self, expected, actual):
    if (expected == None and actual == None):
        return
    self.assertIsNotNone(expected)
    self.assertIsNotNone(actual)
    self.assertEqual(expected.version, actual.version)
    self.assertEqual(expected.project, actual.project)
    self.assertEqual(expected.producer, actual.producer)
    self.assertEqual(expected.run_id, actual.run_id)
    self.assertEqual(expected.reference, actual.reference)
    self.assertEqual(expected.last_executed, actual.last_executed)
    compareInputs(self, expected.inputs, actual.inputs)


def compareInputs(self, expected, actual):
    if (expected == None and actual == None):
        return
    self.assertIsNotNone(expected)
    self.assertIsNotNone(actual)
    self.assertEqual(len(expected), len(actual))
    for expected_plane_uri, actual_plane_uri in zip(expected, actual):
        self.assertEqual(expected_plane_uri, actual_plane_uri)


def compareArtifacts(self, expected, actual):
    if (expected == None and actual == None):
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
        self.assertEqual(expected_artifact.alternative,
                         actual_artifact.alternative)
        compareParts(self, expected_artifact.parts, actual_artifact.parts)


def compareParts(self, expected, actual):
    if (expected == None and actual == None):
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
        compareChunks(self, expected_part.chunks, actual_part.chunks)


def compareChunks(self, expected, actual):
    if (expected == None and actual == None):
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
        compareObservableAxis(self, expected_chunk.observable,
                              actual_chunk.observable)
        compareSpatialWCS(self, expected_chunk.position, actual_chunk.position)
        compareSpectralWCS(self, expected_chunk.energy, actual_chunk.energy)
        compareTemporalWCS(self, expected_chunk.time, actual_chunk.time)
        comparePolarizationWCS(self, expected_chunk.polarization,
                               actual_chunk.polarization)


def compareObservableAxis(self, expected, actual):
    if (expected == None and actual == None):
        return
    self.assertIsNotNone(expected)
    self.assertIsNotNone(actual)
    compareSlice(self, expected.dependent, actual.dependent)
    compareSlice(self, expected.independent, actual.independent)


def compareSpatialWCS(self, expected, actual):
    if (expected == None and actual == None):
        return
    self.assertIsNotNone(expected)
    self.assertIsNotNone(actual)
    compareCoordAxis2D(self, expected.axis, actual.axis)
    self.assertEqual(expected.coordsys, actual.coordsys)
    self.assertEqual(expected.equinox, actual.equinox)
    self.assertEqual(expected.resolution, actual.resolution)


def compareSpectralWCS(self, expected, actual):
    if (expected == None and actual == None):
        return
    self.assertIsNotNone(expected)
    self.assertIsNotNone(actual)
    compareCoordAxis1D(self, expected.axis, actual.axis)
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


def compareTemporalWCS(self, expected, actual):
    if (expected == None and actual == None):
        return
    self.assertIsNotNone(expected)
    self.assertIsNotNone(actual)
    compareCoordAxis1D(self, expected.axis, actual.axis)
    self.assertEqual(expected.exposure, actual.exposure)
    self.assertEqual(expected.resolution, actual.resolution)
    self.assertEqual(expected.timesys, actual.timesys)
    self.assertEqual(expected.trefpos, actual.trefpos)
    self.assertEqual(expected.mjdref, actual.mjdref)


def comparePolarizationWCS(self, expected, actual):
    if (expected == None and actual == None):
        return
    self.assertIsNotNone(expected)
    self.assertIsNotNone(actual)
    compareCoordAxis1D(self, expected.axis, actual.axis)


def compareAxis(self, expected, actual):
    if (expected == None):
        self.assertIsNone(actual)
        return
    self.assertIsNotNone(actual)
    self.assertIsNotNone(actual.ctype)
    self.assertIsNotNone(actual.cunit)
    self.assertEqual(expected.ctype, actual.ctype)
    self.assertEqual(expected.cunit, actual.cunit)


def compareCoord2D(self, expected, actual):
    if (expected == None):
        self.assertIsNone(actual)
        return
    self.assertIsNotNone(actual)
    compareRefCoord(self, expected.coord1, actual.coord1)
    compareRefCoord(self, expected.coord2, actual.coord2)


def compareValueCoord2D(self, expected, actual):
    if (expected == None):
        self.assertIsNone(actual)
        return
    self.assertIsNotNone(actual)
    self.assertEqual(expected.coord1, actual.coord1)
    self.assertEqual(expected.coord2, actual.coord2)


def compareCoordAxis1D(self, expected, actual):
    if (expected == None):
        self.assertIsNone(actual)
        return
    self.assertIsNotNone(actual)
    compareCoordError(self, expected.error, actual.error)
    compareCoordRange1D(self, expected.range, actual.range)
    compareCoordBounds1D(self, expected.bounds, actual.bounds)
    compareCoordFunction1D(self, expected.function, actual.function)


def compareCoordAxis2D(self, expected, actual):
    if (expected == None):
        self.assertIsNone(actual)
        return
    self.assertIsNotNone(actual)
    self.assertIsNotNone(actual.axis1)
    self.assertIsNotNone(actual.axis2)
    compareAxis(self, expected.axis1, actual.axis1)
    compareAxis(self, expected.axis2, actual.axis2)
    compareCoordError(self, expected.error1, actual.error1)
    compareCoordError(self, expected.error2, actual.error2)
    compareCoordRange2D(self, expected.range, actual.range)
    compareCoordBounds2D(self, expected.bounds, actual.bounds)
    compareCoordFunction2D(self, expected.function, actual.function)


def compareCoordBounds1D(self, expected, actual):
    if (expected == None):
        self.assertIsNone(actual)
        return
    self.assertIsNotNone(actual)
    self.assertIsNotNone(expected.samples)
    self.assertIsNotNone(actual.samples)
    self.assertEqual(len(expected.samples), len(actual.samples))
    for expected_range, actual_range in zip(expected.samples, actual.samples):
        compareCoordRange1D(self, expected_range, actual_range)


def compareCoordBounds2D(self, expected, actual):
    if (expected == None):
        self.assertIsNone(actual)
        return
    self.assertIsNotNone(actual)
    if (isinstance(expected, CoordCircle2D) and
        isinstance(actual, CoordCircle2D)):
        compareCoordCircle2D(self, expected, actual)
    elif (isinstance(expected, CoordPolygon2D) and
          isinstance(actual, CoordPolygon2D)):
        compareCoordPolygon2D(self, expected, actual)
    else:
        self.fail("CoordBounds2D expected and actual are different types.")


def compareCoordCircle2D(self, expected, actual):
    if (expected == None):
        self.assertIsNone(actual)
        return
    self.assertIsNotNone(actual)
    self.assertIsNotNone(actual.center)
    self.assertIsNotNone(actual.radius)
    compareValueCoord2D(self, expected.center, actual.center)
    self.assertEqual(expected.radius, actual.radius)


def compareCoordError(self, expected, actual):
    if (expected == None):
        self.assertIsNone(actual)
        return

    self.assertIsNotNone(actual)
    if (expected.syser):
        self.assertIsNotNone(actual.syser)
        self.assertEqual(expected.syser, actual.syser)
    if (expected.rnder):
        self.assertIsNotNone(actual.rnder)
        self.assertEqual(expected.rnder, actual.rnder)


def compareCoordFunction1D(self, expected, actual):
    if (expected == None):
        self.assertIsNone(actual)
        return
    self.assertIsNotNone(actual)
    self.assertEqual(expected.naxis, actual.naxis)
    self.assertEqual(expected.delta, actual.delta)
    compareRefCoord(self, expected.ref_coord, actual.ref_coord)


def compareCoordFunction2D(self, expected, actual):
    if (expected == None):
        self.assertIsNone(actual)
        return
    self.assertIsNotNone(actual)
    self.assertIsNotNone(actual.dimension)
    self.assertIsNotNone(actual.ref_coord)
    self.assertIsNotNone(actual.cd11)
    self.assertIsNotNone(actual.cd12)
    self.assertIsNotNone(actual.cd21)
    self.assertIsNotNone(actual.cd22)
    compareDimension2D(self, expected.dimension, actual.dimension)
    compareCoord2D(self, expected.ref_coord, actual.ref_coord)
    self.assertEqual(expected.cd11, actual.cd11, 0.0)
    self.assertEqual(expected.cd12, actual.cd12, 0.0)
    self.assertEqual(expected.cd21, actual.cd21, 0.0)
    self.assertEqual(expected.cd22, actual.cd22, 0.0)


def compareCoordPolygon2D(self, expected, actual):
    if (expected == None):
        self.assertIsNone(actual)
        return
    self.assertIsNotNone(actual)
    self.assertIsNotNone(expected.vertices)
    self.assertIsNotNone(actual.vertices)
    self.assertEqual(len(expected.vertices), len(actual.vertices))
    for expected_coord_2d, actual_coord_2d in zip(expected.vertices,
                                                  actual.vertices):
        compareValueCoord2D(self, expected_coord_2d, actual_coord_2d)


def compareCoordRange1D(self, expected, actual):
    if (expected == None):
        self.assertIsNone(actual)
        return
    self.assertIsNotNone(actual)
    compareRefCoord(self, expected.start, actual.start)
    compareRefCoord(self, expected.end, actual.end)


def compareCoordRange2D(self, expected, actual):
    if (expected == None):
        self.assertIsNone(actual)
        return
    self.assertIsNotNone(actual)
    self.assertIsNotNone(actual.start)
    self.assertIsNotNone(actual.end)
    compareCoord2D(self, expected.start, actual.start)
    compareCoord2D(self, expected.end, actual.end)


def compareDimension2D(self, expected, actual):
    if (expected == None):
        self.assertIsNone(actual)
        return
    self.assertIsNotNone(actual)
    self.assertEqual(expected.naxis1, actual.naxis1)
    self.assertEqual(expected.naxis2, actual.naxis2)


def compareRefCoord(self, expected, actual):
    if (expected == None):
        self.assertIsNone(actual)
        return
    self.assertIsNotNone(actual)
    self.assertEqual(expected.pix, actual.pix)
    self.assertEqual(expected.val, actual.val)


def compareSlice(self, expected, actual):
    if (expected == None):
        self.assertIsNone(actual)
        return
    self.assertIsNotNone(actual)
    self.assertIsNotNone(actual.bin)
    self.assertIsNotNone(actual.axis)
    self.assertEqual(expected.bin, actual.bin)
    compareAxis(self, expected.axis, actual.axis)


def comparePoint(self, expected, actual):
    if (expected == None):
        self.assertIsNone(actual)
        return
    self.assertIsNotNone(actual)
    self.assertIsNotNone(actual.cval1)
    self.assertIsNotNone(actual.cval2)
    self.assertEqual(expected.cval1, actual.cval1)
    self.assertEqual(expected.cval2, actual.cval2)


if __name__ == '__main__':
    unittest.main()
