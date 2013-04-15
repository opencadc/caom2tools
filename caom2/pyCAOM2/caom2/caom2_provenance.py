# -*- coding: latin-1 -*-
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

"""defines the Provenance class"""


from caom2_plane_uri import PlaneURI
from util.caom2_util import TypedList
from util.caom2_util import TypedSet
from datetime import datetime
from urlparse import urlsplit
import util.caom2_util as util


class Provenance(object):
    """ Provenance """

    def __init__(self, name):
        """
        Initializes a Provenance instance

        Arguments:
        name - name of the provenance
        """

        assert name is not None, "No name provided"
        assert isinstance(name, str), "name is not a str: {0}".format(name)
        self._name = name

        self._version = None
        self._project = None
        self._producer = None
        self._run_id = None
        self._reference = None
        self._last_executed = None

        self._keywords = TypedList((str),)
        self._inputs = TypedSet((PlaneURI),)

    # Properties

    @property
    def name(self):
        """ Name """
        return self._name

    @property
    def version(self):
        """ Version """
        return self._version

    @version.setter
    def version(self, value):
        if value is not None:
            assert isinstance(value, str), (
                    "version is not a str: {0}".format(value))
        self._version = value

    @property
    def project(self):
        """ Project """
        return self._project

    @project.setter
    def project(self, value):
        if value is not None:
            assert isinstance(value, str), (
                    "project is not a str: {0}".format(value))
        self._project = value

    @property
    def producer(self):
        """ Producer """
        return self._producer

    @producer.setter
    def producer(self, value):
        if value is not None:
            assert isinstance(value, str), (
                    "producer is not a str: {0}".format(value))
        self._producer = value

    @property
    def run_id(self):
        """ Run ID """
        return self._run_id

    @run_id.setter
    def run_id(self, value):
        if value is not None:
            assert isinstance(value, str), (
                    "Run ID is not a str: {0}".format(value))
        self._run_id = value

    @property
    def reference(self):
        """ Reference """
        return self._reference

    @reference.setter
    def reference(self, value):
        if value is not None:
            assert isinstance(value, str), (
                    "reference is not a URI: {0}".format(value))
        tmp = urlsplit(value)
        assert tmp.geturl() == value, "Invalid URI: " + value
        self._reference = value

    @property
    def last_executed(self):
        """ Version """
        return self._last_executed

    @last_executed.setter
    def last_executed(self, value):
        if value is not None:
            assert isinstance(value, datetime), (
                    "sample size is not a datetime: {0}".format(value))
        self._last_executed = value

    @property
    def keywords(self):
        """ List of keywords as str"""
        return self._keywords

    @property
    def inputs(self):
        """ Set of inputs as PlaneURI"""
        return self._inputs
