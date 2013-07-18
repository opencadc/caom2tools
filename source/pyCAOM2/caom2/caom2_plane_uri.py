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

"""defines the caom2.PlaneURI class"""

from caom2_observation_uri import ObservationURI
from util.caom2_util import validate_path_component
from urlparse import urlsplit
from urlparse import SplitResult
from caom2_object import Caom2Object
import util.caom2_util as util


class PlaneURI(Caom2Object):
    """ Plane URI """

    def __init__(self, uri):
        """
        Initializes an Plane instance

        Arguments:
        uri : URI corresponding to the plane

        Throws:
        TypeError  : if uri is not a string
        ValueError : if uri is invalid
        ValueError : if the uri is valid but does not contain the expected
        fields (collection, observation_id and product_id)
        """

        self.uri = uri

    def _key(self):
        return (self.uri)

    def __hash__(self):
        return hash(self._key())

    @classmethod
    def get_plane_uri(cls, observation_uri, product_id):
        """
        Initializes an Plane URI instance

        Arguments:
        observation_uri : the uri of the observation
        product_id : ID of the product
        """
        util.typeCheck(observation_uri, ObservationURI, "observation_uri",
                       override=False)
        util.typeCheck(product_id, str, "observation_uri", override=False)
        validate_path_component(cls, "product_id", product_id)

        path = urlsplit(observation_uri.uri).path
        uri = SplitResult(ObservationURI._SCHEME, "", path + "/" +
                          product_id, "", "").geturl()
        return cls(uri)

    # Properties
    @property
    def uri(self):
        """A uri that locates the plane object inside caom"""
        return self._uri

    @uri.setter
    def uri(self, value):

        util.typeCheck(value, str, "uri", override=False)
        tmp = urlsplit(value)

        if tmp.scheme != ObservationURI._SCHEME:
            raise ValueError("{} doesn't have an allowed scheme".format(value))
        if tmp.geturl() != value:
            raise ValueError("Failed to parse uri correctly: {}".format(value))

        (collection, observation_id, product_id) = tmp.path.split("/")

        if product_id is None:
            raise ValueError("Faield to get product ID from uri: {}"
                             .format(value))

        self._product_id = product_id
        self._observation_uri = \
            ObservationURI.get_observation_uri(collection,
                                               observation_id)
        self._uri = value

    @property
    def product_id(self):
        """the product_id associated with this plane"""
        return self._product_id

    @property
    def observation_uri(self):
        """The uri that can be used to find the caom2 observation object that
        this plane belongs to"""
        return self._observation_uri
