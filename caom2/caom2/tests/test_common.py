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

""" Defines TestCaom2IdGenerator class """

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import unittest

from .. import common
from .. import chunk
from .. import part
from .. import plane
from .. import artifact
from .. import observation


class TestCaom2IdGenerator(unittest.TestCase):

    def test_all(self):
        # Not much for now. Just to make sure that all the clients work
        test_entity = common.AbstractCaomEntity()
        print(test_entity._id, test_entity._last_modified)
        test_artifact = artifact.Artifact("caom2:/blah/blah",
                                          chunk.ProductType.SCIENCE,
                                          artifact.ReleaseType.DATA)
        print(test_artifact._id, test_artifact._last_modified)

        test_chunk = chunk.Chunk()
        print(test_chunk._id, test_chunk._last_modified)

        algorithm = observation.Algorithm("myAlg")
        test_observation = observation.Observation("colect", "obs", algorithm)
        print(test_observation._id, test_observation._last_modified)

        test_part = part.Part("part")
        print(test_part._id, test_part._last_modified)

        test_plane = plane.Plane("prodid")
        print(test_plane._id, test_plane._last_modified)


class TestObservationURI(unittest.TestCase):

    def test_all(self):
        obs_uri = observation.ObservationURI("caom:GEMINI/12345")
        self.assertEqual("caom:GEMINI/12345", obs_uri.uri, "Observation URI")
        self.assertEqual("GEMINI", obs_uri.collection, "Collection")
        self.assertEqual("12345", obs_uri.observation_id, "Observation ID")

        obs_uri = observation.ObservationURI.get_observation_uri("CFHT", "654321")
        self.assertEqual("caom:CFHT/654321", obs_uri.uri, "Observation URI")
        self.assertEqual("CFHT", obs_uri.collection, "Collection")
        self.assertEqual("654321", obs_uri.observation_id, "Observation ID")

        exception = False
        try:
            obs_uri = observation.ObservationURI.get_observation_uri(None, "123")
        except TypeError:
            exception = True
        self.assertTrue(exception, "Missing exception")

        exception = False
        try:
            obs_uri = observation.ObservationURI.get_observation_uri("GEMINI", None)
        except TypeError:
            exception = True
        self.assertTrue(exception, "Missing exception")
