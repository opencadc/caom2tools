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

import unittest
import pytest

from builtins import str, int

from .. import artifact
from .. import caom_util
from .. import chunk
from .. import part
from .. import plane


class TestCaomUtil(unittest.TestCase):
    def test_typed_list(self):
        my_list1 = caom_util.TypedList(str, "Test1")
        self.assertEqual(1, len(my_list1), "list1 length")
        self.assertEqual("Test1", my_list1[0], "Non matching elements")

        my_list1.append("Test2")
        self.assertEqual(2, len(my_list1), "list1 length")
        self.assertEqual("Test1", my_list1[0], "Non matching elements")
        self.assertEqual("Test2", my_list1[1], "Non matching elements")

        # try to add the wrong type
        with pytest.raises(TypeError):
            my_list1.append(3)

        with pytest.raises(TypeError):
            my_list1.extend([2, 4])

        with pytest.raises(TypeError):
            my_list2 = caom_util.TypedList(int, 2)
            my_list1.extend(my_list2)

        with pytest.raises(TypeError):
            my_list1.insert(1, 3)

        with pytest.raises(TypeError):
            my_list2 = caom_util.TypedList(str, 1, 3)

        self.assertEqual(2, len(my_list1), "list1 length")
        self.assertEqual("Test1", my_list1[0], "Non matching elements")
        self.assertEqual("Test2", my_list1[1], "Non matching elements")

        my_list2 = caom_util.TypedList(str, "Test3", "Test4")
        my_list1.extend(my_list2)
        self.assertEqual(4, len(my_list1), "list1 length")
        self.assertEqual("Test1", my_list1[0], "Non matching elements")
        self.assertEqual("Test2", my_list1[1], "Non matching elements")
        self.assertEqual("Test3", my_list1[2], "Non matching elements")
        self.assertEqual("Test4", my_list1[3], "Non matching elements")

        my_list1.insert(0, "Test0")
        self.assertEqual(5, len(my_list1), "list1 length")
        self.assertEqual("Test0", my_list1[0], "Non matching elements")
        self.assertEqual("Test1", my_list1[1], "Non matching elements")
        self.assertEqual("Test2", my_list1[2], "Non matching elements")
        self.assertEqual("Test3", my_list1[3], "Non matching elements")
        self.assertEqual("Test4", my_list1[4], "Non matching elements")

        my_list2 = caom_util.TypedList(plane.Energy, )
        self.assertEqual(0, len(my_list2), "list2 length")

    def test_validate_path_component(self):
        energy = plane.Energy()
        caom_util.validate_path_component(energy, "something",
                                          "some:test\\path")

        with pytest.raises(ValueError):
            caom_util.validate_path_component(energy, "energyfield",
                                              "some:test path")

        with pytest.raises(ValueError):
            caom_util.validate_path_component(energy, "energyfield",
                                              "some:test||path")

        with pytest.raises(ValueError):
            caom_util.validate_path_component(energy, "energyfield",
                                              "some:test %path")

    def test_typed_set(self):

        my_set = caom_util.TypedSet(str, )
        with self.assertRaises(TypeError):
            my_set.add(float(1.0))
            my_set.add(int(1))
            my_set.add(bool(1))

        self.assertRaises(TypeError, caom_util.TypedSet, str, float(1.0))

        my_set = caom_util.TypedSet(str, "Test1")
        self.assertEqual(1, len(my_set))
        self.assertEqual("Test1", my_set.pop())

        my_set.add("Test1")
        my_set.add("Test2")
        self.assertEqual(2, len(my_set), "set length")
        with self.assertRaises(TypeError):
            my_set.add(float(1.0))
            my_set.add(int(1))
            my_set.add(bool(1))
            my_set.add("Test1")

        my_set = caom_util.TypedSet(str, )
        my_set.add("Test1")
        my_set.add("Test1")
        self.assertTrue(len(my_set) == 1)

    def test_typed_ordered_dict(self):

        # test validation and constructor with an empty dictionary
        test_plane10 = plane.Plane('key10')
        test_artifact66 = artifact.Artifact("caom:CFHT/55/66",
                                            chunk.ProductType.SCIENCE,
                                            artifact.ReleaseType.DATA)
        test_part10 = part.Part("10")
        test_plane_uri = plane.PlaneURI('caom:CFHT/55/66')
        my_dict_plane = caom_util.TypedOrderedDict(plane.Plane, )
        with self.assertRaises(ValueError):
            my_dict_plane['key11'] = test_plane10
        my_dict_artifact = caom_util.TypedOrderedDict(artifact.Artifact, )
        with self.assertRaises(ValueError):
            my_dict_artifact['caom:CFHT/55/6'] = test_artifact66
        my_dict_part = caom_util.TypedOrderedDict(part.Part, )
        with self.assertRaises(ValueError):
            my_dict_part['11'] = test_part10
        my_dict_wrong_type = caom_util.TypedOrderedDict(plane.PlaneURI, )
        with self.assertRaises(ValueError):
            my_dict_wrong_type['caom:CFHT/55/67'] = test_plane_uri
        with self.assertRaises(TypeError):
            my_dict_plane['key2'] = 'value2'
        with self.assertRaises(TypeError):
            my_dict_plane['key1'] = float(2.0)
        # test assignment
        my_dict = caom_util.TypedOrderedDict(plane.Plane, )
        test_plane2 = plane.Plane('key2')
        test_plane1 = plane.Plane('key1')
        my_dict['key2'] = test_plane2
        my_dict['key1'] = test_plane1
        # need to cast to list in order to make it work with both python
        # 2 and 3
        self.assertEqual(2, len(my_dict),
                         'mismatch in the number of entries in dictionary.')
        self.assertEqual('key2', list(my_dict.keys())[0],
                         'key mismatch for 1st key')
        self.assertEqual('key1', list(my_dict.keys())[1],
                         'key mismatch for 2nd key')
        self.assertEqual(test_plane2, list(my_dict.values())[0],
                         'value mismatch for 1st value')
        self.assertEqual(test_plane1, list(my_dict.values())[1],
                         'value mismatch for 2nd value')
        # test constructor with non-empty dictionary
        test_plane1 = plane.Plane('key1')
        test_plane2 = plane.Plane('key2')
        my_dict1 = caom_util.TypedOrderedDict(plane.Plane,
                                              ('key1', test_plane1),
                                              ('key2', test_plane2))
        self.assertEqual(2, len(my_dict1),
                         'mismatch in the number of entries in dictionary.')
        # test assignment via setdefault
        self.assertRaises(TypeError, my_dict1.setdefault,
                          'key3', 'wrong value')
        my_dict1.setdefault('key3', plane.Plane('key3'))
        self.assertEqual(3, len(my_dict1),
                         'mismatch in the number of entries in dictionary.')
        # test assignment via update
        my_dict1.update(my_dict)
        self.assertEqual(3, len(my_dict1),
                         'mismatch in the number of entries in dictionary.')
        self.assertEqual('key2', list(my_dict.keys())[0],
                         'key mismatch for 1st key')
        self.assertEqual('key1', list(my_dict.keys())[1],
                         'key mismatch for 2nd key')

        # test add function
        my_dict1.add(plane.Plane('key4'))
        self.assertEqual(4, len(my_dict1),
                         'mismatch in the number of entries in dictionary.')
        self.assertEqual('key1', list(my_dict1.keys())[0],
                         'key mismatch for 1st key')
        self.assertEqual('key2', list(my_dict1.keys())[1],
                         'key mismatch for 2nd key')
        self.assertEqual('key3', list(my_dict1.keys())[2],
                         'key mismatch for 3rd key')
        self.assertEqual('key4', list(my_dict1.keys())[3],
                         'key mismatch for 4th key')

        plane5 = plane.Plane("key5")
        my_dict1[plane5._key()] = plane5

        with self.assertRaises(TypeError):
            my_dict1.add(test_plane_uri)

        # test pop function
        self.assertEqual(5, len(my_dict1),
                         'mismatch in the number of entries in dictionary.')
        my_value = my_dict1.pop('key4')
        self.assertEqual('key4', my_value._key(),
                         'popped the wrong entry from dictionary.')
        self.assertEqual(4, len(my_dict1),
                         'mismatch in the number of entries in dictionary.')
        my_value = my_dict1.pop('key5')
        self.assertEqual('key5', my_value._key(),
                         'popped the wrong entry from dictionary.')
        self.assertEqual(3, len(my_dict1),
                         'mismatch in the number of entries in dictionary.')
        my_value = my_dict1.pop('key3')
        self.assertEqual('key3', my_value._key(),
                         'popped the wrong entry from dictionary.')
        self.assertEqual(2, len(my_dict1),
                         'mismatch in the number of entries in dictionary.')
        my_value = my_dict1.pop('key2')
        self.assertEqual('key2', my_value._key(),
                         'popped the wrong entry from dictionary.')
        self.assertEqual(1, len(my_dict1),
                         'mismatch in the number of entries in dictionary.')
        my_value = my_dict1.pop('key1')
        self.assertEqual('key1', my_value._key(),
                         'popped the wrong entry from dictionary.')
