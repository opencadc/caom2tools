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

from caom2utils import fits2caom2
from caom2 import ObservationReader
from caom2.diff import get_differences

import glob
import logging
import os
import sys
import tempfile
import traceback

import pytest


THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TESTDATA_DIR = os.path.join(THIS_DIR, 'data')


def test_differences():
    """This test assumes a directory structure of
    tests/data/'collection'/'test case'/files. """

    collections = _get_subdirs(TESTDATA_DIR)
    for collection in collections:
        test_cases = _get_subdirs(os.path.join(TESTDATA_DIR, collection))
        for test_case in test_cases:
            location = os.path.join(os.path.join(TESTDATA_DIR, collection),
                                    test_case)
            product_id = test_case
            collection_id = collection.upper()
            config = _get_parameter('config', location)
            defaults = _get_parameter('default', location)
            overrides = _get_parameter('override', location)
            data_files = _get_file('header', location)
            data_files_parameter = _get_data_files_parameter(data_files)
            expected_fname = _get_file('xml', location)
            obs_id = _get_common(data_files)
            uris = _get_uris(collection_id, data_files)
            temp = tempfile.NamedTemporaryFile()
            sys.argv = ('fits2caom2 --debug --dumpconfig '
                        '{} -o {} '
                        '--observation {} {} {} {} {} {} {} '.format(
                    data_files_parameter, temp.name, collection_id, obs_id,
                        config, defaults, overrides, product_id, uris)).split()
            try:
                print(sys.argv)
                fits2caom2.main_app()
                _compare_observations(expected_fname, temp.name, location,
                                      obs_id)
            except:
                print(
                    'Execution for collection \'{}\' test case \'{}\' failed with {}'.format(
                        collection_id, obs_id, sys.exc_info()[1]))
                logging.exception('test_differences')


def _get_common(fnames):
    common = os.path.basename(fnames[0])
    for jj in fnames:
        rhs = os.path.basename(jj)
        for kk in fnames:
            lhs = os.path.basename(kk)
            ii = 0
            while ii < len(rhs):
                if rhs[ii] == '.' or rhs[ii] != lhs[ii]:
                    if len(rhs[0:ii]) != 0 and len(rhs[0:ii]) < len(common):
                        common = rhs[0:ii]
                    break
                else:
                    ii += 1
    return common


def _get_subdirs(dir_name):
    return [name for name in os.listdir(dir_name) if
            os.path.isdir(os.path.join(dir_name, name))]


def _get_parameter(extension, dir_name):
    fnames = _get_file(extension, dir_name)
    if fnames:
        result = '--{} {}'.format(extension, fnames[0])
        return result
    else:
        return None


def _get_data_files_parameter(fnames):
    if fnames:
        result = '--local'
        for fname in fnames:
            result = '{} {}'.format(result, fname)
        return result
    else:
        return None


def _get_uris(collection, fnames):
    if fnames:
        result = ''
        for fname in fnames:
            result = '{} ad:{}/{}'.format(result, collection,
                                          os.path.basename(fname))
        return result
    else:
        return None


def _get_file(pattern, dir_name):
    files = glob.glob('{}/*.{}'.format(dir_name, pattern))
    if files:
        return files
    else:
        return None


def _compare_observations(expected_fname, actual_fname, output_dir, obsid):
    expected = _read_observation(expected_fname)
    actual = _read_observation(actual_fname)
    result = get_differences(expected, actual, 'Observation')
    result_fname = os.path.join(output_dir, '{}.compare'.format(obsid))
    print(
        '{} differences between expected and actual, in {}'.format(
            len(result), result_fname))
    f = open(result_fname, 'w')
    f.write('\n'.join(str(p) for p in result))
    f.write('\n')
    f.close()


def _read_observation(fname):
    reader = ObservationReader(False)
    result = reader.read(fname)
    return result
