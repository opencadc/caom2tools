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

from .. import artifact
from .. import part


class TestEnums(unittest.TestCase):

    def test_all(self):
        # test for invalid value
        self.assertEqual(artifact.ProductType.get("no_such_string"), None)
        self.assertRaises(AttributeError, artifact.ProductType.get, None)
        self.assertRaises(AttributeError, artifact.ProductType.get, 1)

        self.assertEqual(artifact.ReleaseType.get("no_such_string"), None)
        self.assertRaises(AttributeError, artifact.ReleaseType.get, None)
        self.assertRaises(AttributeError, artifact.ReleaseType.get, 1)

        # test that we can get the object for each enum by name
        self.assertEqual(artifact.ProductType.SCIENCE.name, "SCIENCE")
        self.assertEqual(artifact.ProductType.get(
            artifact.ProductType.SCIENCE.name).name, "SCIENCE")
        self.assertEqual(artifact.ProductType.get(
            'SCIENCE').value, "science")
        self.assertEqual(artifact.ProductType.get(
            artifact.ProductType.SCIENCE.name).value, "science")
        self.assertEqual(artifact.ProductType.getByValue(
            artifact.ProductType.SCIENCE.value).value, "science")
        self.assertEqual(artifact.ProductType.getByValue(
            artifact.ProductType.SCIENCE.value).name, "SCIENCE")

        self.assertEqual(artifact.ProductType.SCIENCE.value, "science")
        self.assertEqual(artifact.ProductType.CALIBRATION.value, "calibration")
        self.assertEqual(artifact.ProductType.PREVIEW.value, "preview")
        self.assertEqual(artifact.ProductType.INFO.value, "info")
        self.assertEqual(artifact.ProductType.NOISE.value, "noise")
        self.assertEqual(artifact.ProductType.WEIGHT.value, "weight")
        self.assertEqual(artifact.ProductType.AUXILIARY.value, "auxiliary")
        self.assertEqual(artifact.ProductType.THUMBNAIL.value, "thumbnail")

        self.assertEqual(artifact.ReleaseType.DATA.value, "data")
        self.assertEqual(artifact.ReleaseType.META.value, "meta")


class TestArtifact(unittest.TestCase):

    def test_all(self):
        with self.assertRaises(TypeError):
            test_artifact = artifact.Artifact("caom:GEMINI/12345")
        with self.assertRaises(TypeError):
            test_artifact = artifact.Artifact("caom:GEMINI/12345",
                                              artifact.ReleaseType('META'),
                                              artifact.ProductType('THUMBNAIL'))
        with self.assertRaises(TypeError):
            test_artifact = artifact.Artifact("caom:GEMINI/12345",
                                              artifact.ProductType('THUMBNAIL'),
                                              None)
        with self.assertRaises(TypeError):
            test_artifact = artifact.Artifact("caom:GEMINI/12345",
                                              None,
                                              artifact.ReleaseType('META'))

        test_artifact = artifact.Artifact("caom:GEMINI/12345",
                                          artifact.ProductType('THUMBNAIL'),
                                          artifact.ReleaseType('META'))
        urlparse("caom:GEMINI/12345")
        self.assertEqual("caom:GEMINI/12345",
                         test_artifact.uri,
                         "Artifact URI")
        self.assertEqual(artifact.ProductType('THUMBNAIL'),
                         test_artifact.product_type,
                         "Artifact ProductType")
        self.assertEqual(artifact.ReleaseType('META'),
                         test_artifact.release_type,
                         "Artifact ReleaseType")

        self.assertIsNone(test_artifact.content_type, "Default content type")
        test_artifact.content_type = "FITS"
        self.assertEquals("FITS", test_artifact.content_type, "Content type")
        self.assertIsNone(test_artifact.content_length,
                          "Default content length")
        test_artifact.content_length = 23L
        self.assertEquals(23L, test_artifact.content_length, "Content length")
        test_artifact.product_type = artifact.ProductType.PREVIEW
        self.assertEquals(artifact.ProductType.PREVIEW,
                          test_artifact.product_type,
                          "Product type")
        self.assertEquals(0, len(test_artifact.parts), "Default parts")
        part1 = part.Part("1")
        test_artifact.parts["1"] = part1
        self.assertEquals(1, len(test_artifact.parts), "Parts")
        self.assertTrue("1" in test_artifact.parts.keys())
        #add same part again
        part2 = part.Part("2")
        test_artifact.parts["2"] = part2
        self.assertEquals(2, len(test_artifact.parts), "Parts")
        self.assertTrue("1" in test_artifact.parts.keys())
        self.assertTrue("2" in test_artifact.parts.keys())

        # try to add duplicates
        part3 = part1
        test_artifact.parts["1"] = part3
        self.assertEquals(2, len(test_artifact.parts), "Parts")
        self.assertTrue("1" in test_artifact.parts.keys())
        self.assertTrue("2" in test_artifact.parts.keys())

        part4 = part.Part("1")
        test_artifact.parts["1"] = part4
        self.assertEquals(2, len(test_artifact.parts), "Parts")
        self.assertTrue("1" in test_artifact.parts.keys())
        self.assertTrue("2" in test_artifact.parts.keys())

        #incorrect URI
        exception = False
        try:
            test_artifact = artifact.Artifact(
                "caom://#observation://? something#//",
                artifact.ReleaseType('META'),
                artifact.ProductType('THUMBNAIL'))
            print artifact.uri
        except ValueError:
            exception = True
        self.assertTrue(exception, "Missing exception")