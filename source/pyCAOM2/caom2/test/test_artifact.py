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

""" Defines TestArtifact class """

import unittest
from urlparse import urlparse

from caom2.caom2_artifact import Artifact
from caom2.caom2_enums import ProductType, ReleaseType
from caom2.caom2_part import Part


class TestArtifact(unittest.TestCase):

    def test_all(self):
        with self.assertRaises(TypeError):
            artifact = Artifact("caom:GEMINI/12345")
        with self.assertRaises(TypeError):
            artifact = Artifact("caom:GEMINI/12345", ReleaseType('META'), ProductType('THUMBNAIL'))
        with self.assertRaises(TypeError):
            artifact = Artifact("caom:GEMINI/12345", ProductType('THUMBNAIL'), None)
        with self.assertRaises(TypeError):
            artifact = Artifact("caom:GEMINI/12345", None, ReleaseType('META'))

        artifact = Artifact("caom:GEMINI/12345", ProductType('THUMBNAIL'), ReleaseType('META'))
        urlparse("caom:GEMINI/12345")
        self.assertEqual("caom:GEMINI/12345", artifact.uri, "Artifact URI")
        self.assertEqual(ProductType('THUMBNAIL'), artifact.product_type, "Artifact ProductType")
        self.assertEqual(ReleaseType('META'), artifact.release_type, "Artifact ReleaseType")

        self.assertIsNone(artifact.content_type, "Default content type")
        artifact.content_type = "FITS"
        self.assertEquals("FITS", artifact.content_type, "Content type")
        self.assertIsNone(artifact.content_length, "Default content length")
        artifact.content_length = 23L
        self.assertEquals(23L, artifact.content_length, "Content length")
        artifact.product_type = ProductType.PREVIEW
        self.assertEquals(ProductType.PREVIEW, artifact.product_type, "Product type")
        self.assertEquals(0, len(artifact.parts), "Default parts")
        part1 = Part("1")
        artifact.parts["1"] = part1
        self.assertEquals(1, len(artifact.parts), "Parts")
        self.assertTrue("1" in artifact.parts.keys())
        #add same part again
        part2 = Part("2")
        artifact.parts["2"] = part2
        self.assertEquals(2, len(artifact.parts), "Parts")
        self.assertTrue("1" in artifact.parts.keys())
        self.assertTrue("2" in artifact.parts.keys())

        # try to add duplicates
        part3 = part1
        artifact.parts["1"] = part3
        self.assertEquals(2, len(artifact.parts), "Parts")
        self.assertTrue("1" in artifact.parts.keys())
        self.assertTrue("2" in artifact.parts.keys())

        part4 = Part("1")
        artifact.parts["1"] = part4
        self.assertEquals(2, len(artifact.parts), "Parts")
        self.assertTrue("1" in artifact.parts.keys())
        self.assertTrue("2" in artifact.parts.keys())

        #incorrect URI
        exception = False
        try:
            artifact = Artifact("caom://#observation://? something#//", ReleaseType('META'), ProductType('THUMBNAIL'))
            print artifact.uri
        except ValueError:
            exception = True
        self.assertTrue(exception, "Missing exception")


if __name__ == '__main__':
    unittest.main()

