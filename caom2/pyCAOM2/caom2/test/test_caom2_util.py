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

""" Defines TestCaomUtil class """

from caom2.util.caom2_util import TypedList
from caom2.util.caom2_util import TypedOrderedDict
from caom2.util.caom2_util import TypedSet
from caom2.util.caom2_util import validate_path_component
from caom2.caom2_artifact import Artifact
from caom2.caom2_energy import Energy
from caom2.caom2_part import Part
from caom2.caom2_plane import Plane
from caom2.caom2_plane_uri import PlaneURI

import os
import sys
import unittest

# put build at the start of the search path
sys.path.insert(0, os.path.abspath('../../lib.local/lib'))


class TestCaomUtil(unittest.TestCase):

    def testTypedList(self):
        mylist1 = TypedList((str), "Test1")
        self.assertEquals(1, len(mylist1), "list1 length")
        self.assertEqual("Test1", mylist1[0], "Non matching elements")

        mylist1.append("Test2")
        self.assertEquals(2, len(mylist1), "list1 length")
        self.assertEqual("Test1", mylist1[0], "Non matching elements")
        self.assertEqual("Test2", mylist1[1], "Non matching elements")

        # try to add the wrong type
        exception = False
        try:
            mylist1.append(3)
        except AssertionError:
            exception = True
        self.assertTrue(exception, "Exception thrown")

        exception = False
        try:
            mylist1.extend([2, 4])
        except AssertionError:
            exception = True
        self.assertTrue(exception, "Exception thrown")

        exception = False
        try:
            mylist2 = TypedList((int), 2)
            mylist1.extend(mylist2)
        except AssertionError:
            exception = True
        self.assertTrue(exception, "Exception thrown")

        exception = False
        try:
            mylist1.insert(1, 3)
        except AssertionError:
            exception = True
        self.assertTrue(exception, "Exception thrown")

        exception = False
        try:
            mylist2 = TypedList((str), 1, 3)
        except AssertionError:
            exception = True
        self.assertTrue(exception, "Exception thrown")

        self.assertEquals(2, len(mylist1), "list1 length")
        self.assertEqual("Test1", mylist1[0], "Non matching elements")
        self.assertEqual("Test2", mylist1[1], "Non matching elements")

        mylist2 = TypedList((str), "Test3", "Test4")
        mylist1.extend(mylist2)
        self.assertEquals(4, len(mylist1), "list1 length")
        self.assertEqual("Test1", mylist1[0], "Non matching elements")
        self.assertEqual("Test2", mylist1[1], "Non matching elements")
        self.assertEqual("Test3", mylist1[2], "Non matching elements")
        self.assertEqual("Test4", mylist1[3], "Non matching elements")

        mylist1.insert(0, "Test0")
        self.assertEquals(5, len(mylist1), "list1 length")
        self.assertEqual("Test0", mylist1[0], "Non matching elements")
        self.assertEqual("Test1", mylist1[1], "Non matching elements")
        self.assertEqual("Test2", mylist1[2], "Non matching elements")
        self.assertEqual("Test3", mylist1[3], "Non matching elements")
        self.assertEqual("Test4", mylist1[4], "Non matching elements")

        mylist2 = TypedList((Energy),)
        self.assertEquals(0, len(mylist2), "list2 length")

    def test_validate_path_component(self):
        energy = Energy()
        validate_path_component(energy, "something", "some:test\\path")

        exception = False
        try:
            validate_path_component(energy, "energyfield", "some:test path")
        except AssertionError:
            exception = True
        self.assertTrue(exception, "Missing exception")

        exception = False
        try:
            validate_path_component(energy, "energyfield", "some:test/path")
        except AssertionError:
            exception = True
        self.assertTrue(exception, "Missing exception")

        exception = False
        try:
            validate_path_component(energy, "energyfield", "some:test||path")
        except AssertionError:
            exception = True
        self.assertTrue(exception, "Missing exception")

        exception = False
        try:
            validate_path_component(energy, "energyfield", "some:test %path")
        except AssertionError:
            exception = True
        self.assertTrue(exception, "Missing exception")

    def testTypedSet(self):

        myset = TypedSet((str),)
        with self.assertRaises(AssertionError):
            myset.add(float(1.0))
            myset.add(long(1))
            myset.add(bool(1))

        self.assertRaises(AssertionError, TypedSet, (str), float(1.0))

        myset = TypedSet((str), "Test1")
        self.assertEquals(1, len(myset))
        self.assertEqual("Test1", myset.pop())

        myset.add("Test1")
        myset.add("Test2")
        self.assertEquals(2, len(myset), "set length")
        with self.assertRaises(AssertionError):
            myset.add(float(1.0))
            myset.add(long(1))
            myset.add(bool(1))
            myset.add("Test1")

        myset = TypedSet((str),)
        myset.add("Test1")
        myset.add("Test1")
        self.assertTrue(len(myset) == 1)

    def testTypedOrderedDict(self):

        # test validation and constructor with an empty dictionary
        testPlane10 = Plane('key10')
        testArtifact66 = Artifact('caom:CFHT/55/66')
        testPart10 = Part("10")
        testPlaneURI = PlaneURI('caom:CFHT/55/66')
        myDictPlane = TypedOrderedDict((Plane),)
        with self.assertRaises(ValueError):
            myDictPlane['key11'] = testPlane10
        myDictArtifact = TypedOrderedDict((Artifact),)
        with self.assertRaises(ValueError):
            myDictArtifact['caom:CFHT/55/6'] = testArtifact66
        myDictPart = TypedOrderedDict((Part),)
        with self.assertRaises(ValueError):
            myDictPart['11'] = testPart10
        myDictWrongType = TypedOrderedDict((PlaneURI),)
        with self.assertRaises(AttributeError):
            myDictWrongType['caom:CFHT/55/66'] = testPlaneURI
        with self.assertRaises(TypeError):
            myDictPlane['key2'] = 'value2'
        with self.assertRaises(TypeError):
            myDictPlane['key1'] = float(2.0)
        # test assignment
        myDict = TypedOrderedDict((Plane),)
        testPlane2 = Plane('key2')
        testPlane1 = Plane('key1')
        myDict['key2'] = testPlane2
        myDict['key1'] = testPlane1
        self.assertEqual(2, len(myDict),
                         'mismatch in the number of entries in dictionary.')
        self.assertEqual('key2', myDict.keys()[0],
                         'key mismatch for 1st key')
        self.assertEqual('key1', myDict.keys()[1],
                         'key mismatch for 2nd key')
        self.assertEqual(testPlane2, myDict.values()[0],
                         'value mismatch for 1st value')
        self.assertEqual(testPlane1, myDict.values()[1],
                         'value mismatch for 2nd value')
        # test constructor with non-empty dictionary
        testPlane1 = Plane('key1')
        testPlane2 = Plane('key2')
        myDict1 = TypedOrderedDict((Plane), ('key1', testPlane1),
                                            ('key2', testPlane2))
        self.assertEqual(2, len(myDict1),
                         'mismatch in the number of entries in dictionary.')
        # test assignment via setdefault
        self.assertRaises(TypeError, myDict1.setdefault,
                          'key3', 'wrong value')
        myDict1.setdefault('key3', Plane('key3'))
        self.assertEqual(3, len(myDict1),
                         'mismatch in the number of entries in dictionary.')
        # test assignment via update
        myDict1.update(myDict)
        self.assertEqual(3, len(myDict1),
                         'mismatch in the number of entries in dictionary.')
        self.assertEqual('key2', myDict.keys()[0],
                         'key mismatch for 1st key')
        self.assertEqual('key1', myDict.keys()[1],
                         'key mismatch for 2nd key')
        #self.assertEqual(testPlane1, myDict.values()[0],
        #                 'value mismatch for 1st value')
        #self.assertEqual(testPlane2, myDict.values()[1],
        #                 'value mismatch for 2nd value')

        # test add function
        myDict1.add(Plane('key4'))
        self.assertEqual(4, len(myDict1),
                         'mismatch in the number of entries in dictionary.')
        self.assertEqual('key1', myDict1.keys()[0],
                         'key mismatch for 1st key')
        self.assertEqual('key2', myDict1.keys()[1],
                         'key mismatch for 2nd key')
        self.assertEqual('key3', myDict1.keys()[2],
                         'key mismatch for 3rd key')
        self.assertEqual('key4', myDict1.keys()[3],
                         'key mismatch for 4th key')
        # self.assertEqual(testPlane1, myDict1.values()[0],
        #                  'value mismatch for 1st value')
        #self.assertEqual(testPlane2, myDict1.values()[1],
        #                 'value mismatch for 2nd value')
        plane5 = Plane("key5")
        myDict1[plane5.key] = plane5

        with self.assertRaises(AttributeError):
            myDict1.add(testPlaneURI)



if __name__ == '__main__':
    unittest.main()

