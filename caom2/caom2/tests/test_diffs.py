# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2022.                            (c) 2022.
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

import os
import unittest

from caom2 import Point, shape, Vertex, SegmentType, Position

from .. import diff
from .. import observation
from .. import obs_reader_writer
from . import caom_test_instances


class TestCaomUtil(unittest.TestCase):
    def test_get_differences(self):
        expected_simple = observation.SimpleObservation(
            collection='test_collection',
            observation_id='test_observation_id',
            algorithm=observation.Algorithm('EXPOSURE'))
        # meta_producer field is ignored by get_difference
        expected_simple.meta_producer = 'meta_producer/0.1'
        report = diff.get_differences(expected_simple, expected_simple,
                                      'obs')
        self.assertTrue(report is None, repr(report))

        actual_simple = observation.SimpleObservation(
            collection='test_collection',
            observation_id='test_observation_id',
            algorithm=observation.Algorithm('EXPOSURE'))
        report = diff.get_differences(expected_simple, actual_simple,
                                      'obs')
        self.assertTrue(report is None, repr(report))

        act_plane = observation.Plane(product_id='test_product_id1')
        actual_simple.planes['test_product_id1'] = act_plane

        report = diff.get_differences(expected_simple, actual_simple,
                                      'obs')
        self.assertTrue(report is not None, repr(report))
        self.assertTrue(len(report) == 1, repr(report))

        ex_plane = observation.Plane(product_id='test_product_id2')
        expected_simple.planes['test_product_id2'] = ex_plane
        report = diff.get_differences(expected_simple, actual_simple,
                                      'obs')
        self.assertTrue(report is not None, repr(report))
        self.assertTrue(len(report) == 2, repr(report))

        instances = caom_test_instances.Caom2TestInstances()
        instances.complete = True
        obs1 = instances.get_composite_observation()
        obs2 = instances.get_composite_observation()

        report = diff.get_differences(obs1, obs2, 'caom_test_instances')
        assert report is None

        obs3 = instances.get_simple_observation()

        report = diff.get_differences(obs1, obs3, 'caom_test_instances')
        assert len(report) == 1

    def test_samples(self):
        instances = caom_test_instances.Caom2TestInstances()
        seq1 = instances.get_coord_axis1d()
        seq2 = instances.get_coord_axis1d()

        report = diff.get_differences(seq1, seq2, 'samples')
        assert report is None

        for k, v in enumerate(seq2.bounds.samples):
            seq2.bounds.samples[k].pix = 0.0
        report = diff.get_differences(seq1, seq2, 'samples')
        assert report is not None
        assert len(report) == 2

    def test_chunks(self):
        instances = caom_test_instances.Caom2TestInstances()
        seq1 = instances.get_chunks()
        seq2 = instances.get_chunks()

        report = diff.get_differences(seq1, seq2, 'chunks')
        assert report is None

        seq2[0].observable.independent.bin = 0
        report = diff.get_differences(seq1, seq2, 'chunks')
        assert report is not None
        assert len(report) == 1

    def test_plane_level_position(self):
        # special handling, because nan == nan is False
        p1 = [Point(cval1=float('nan'), cval2=float('nan')),
              Point(cval1=100.25, cval2=-30.5),
              Point(cval1=100.25, cval2=30.0),
              Point(cval1=float('nan'), cval2=float('nan'))]
        p2 = [Point(cval1=float('nan'), cval2=float('nan')),
              Point(cval1=100.25, cval2=-30.5),
              Point(cval1=100.25, cval2=30.0),
              Point(cval1=float('nan'), cval2=float('nan'))]

        v1 = [Vertex(p1[0].cval1, p1[0].cval2, SegmentType.MOVE),
              Vertex(p1[1].cval1, p1[1].cval2, SegmentType.LINE),
              Vertex(p1[2].cval1, p1[2].cval2, SegmentType.LINE),
              Vertex(p1[3].cval1, p1[3].cval2, SegmentType.LINE),
              Vertex(p1[0].cval1, p1[0].cval2, SegmentType.CLOSE)]
        v2 = [Vertex(p2[0].cval1, p2[0].cval2, SegmentType.MOVE),
              Vertex(p2[1].cval1, p2[1].cval2, SegmentType.LINE),
              Vertex(p2[2].cval1, p2[2].cval2, SegmentType.LINE),
              Vertex(p2[3].cval1, p2[3].cval2, SegmentType.LINE),
              Vertex(p2[0].cval1, p2[0].cval2, SegmentType.CLOSE)]

        poly1 = shape.Polygon(points=p1, samples=shape.MultiPolygon(v1))
        poly2 = shape.Polygon(points=p2, samples=shape.MultiPolygon(v2))

        o1 = Position(time_dependent=False, bounds=poly1)
        o2 = Position(time_dependent=False, bounds=poly2)

        report = diff.get_differences(o1, o2, 'caom test instances')
        assert report is None, 'NaN comparison failure'

    def test_xml(self):
        test_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                'data')
        actual_fqn = os.path.join(test_dir, 'diff-actual-CAOM-2.3.xml')
        expected_fqn = os.path.join(test_dir, 'diff-expected-CAOM-2.3.xml')
        reader = obs_reader_writer.ObservationReader(False)
        actual_obs = reader.read(actual_fqn)
        expected_obs = reader.read(expected_fqn)
        report = diff.get_differences(expected_obs, actual_obs)
        assert report is None, report

    def test_xml_sequence(self):
        test_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                'data')
        actual_fqn = os.path.join(test_dir, 'diff-actual-CAOM-2.4.xml')
        expected_fqn = os.path.join(test_dir, 'diff-expected-CAOM-2.4.xml')
        reader = obs_reader_writer.ObservationReader(False)
        actual_obs = reader.read(actual_fqn)
        expected_obs = reader.read(expected_fqn)
        report = diff.get_differences(expected_obs, actual_obs)
        assert report is None, report
