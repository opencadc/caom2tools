# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2018.                            (c) 2018.
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
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import unittest

from .. import diff
from .. import observation
from . import caom_test_instances


class TestCaomUtil(unittest.TestCase):
    def test_get_differences(self):
        expected_simple = observation.SimpleObservation(
            collection='test_collection',
            observation_id='test_observation_id',
            algorithm=observation.Algorithm('EXPOSURE'))
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
        self.assertTrue(len(report) == 2, repr(report))

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
