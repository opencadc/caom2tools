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

"""Defines the caom2.Artifact class.

"""

from caom2_entity import AbstractCaom2Entity
from caom2_enums import ProductType
from caom2_part import Part
from util.caom2_util import TypedOrderedDict
from urlparse import urlparse
import util.caom2_util as util

class Artifact(AbstractCaom2Entity):
    """Contains the meta data assocaited with a file. 

    - location of the file (uri)
    - the http content-type
    - size of the file (content-lenght)

    As well as a pointer (parts) to content of the file. 

    eg:  Artificat('ad:CFHT/1234567o')
    where 'ad:CFHT/1234567o' is a uri that refernce the file... 

    """

    def __init__(self, 
                 uri,
                 content_type=None,
                 content_length=None,
                 product_type=None,
                 alternative=False,
                 parts=None
                 ):
        """
        Initialize a Artifact instance.

        Arguments: uri of the artifact.  eg:
           vos://cadc.nrc.ca!vospace/APASS/apass_north/proc/100605/n100605.0384.new.fz
           ad:CFHT/123456p
        
        """
        super(Artifact, self).__init__()
        self.uri = uri
        self.content_type = content_type
        self.content_length = content_length
        self.product_type = product_type
        self.alternative = alternative
        if parts is None:
            parts = TypedOrderedDict((Part),)
        self.parts = parts

    def _key(self):
        return (self.uri)

    def __hash__(self):
        return hash(self._key())

    @property
    def content_type(self):
        """content-type (http header style) of the content represented by this
        artifact.

        """
        return self._content_type

    @property
    def key(self):
        """Dictionary key for artifact is its URI 

        """
        return self._uri

    @content_type.setter
    def content_type(self, value):
        util.typeCheck(value, str, "content_type" ) 
        self._content_type = value

    @property
    def content_length(self):
        """size of the artifact.

        Unit: byte
        type: int

        """
        return self._content_length

    @content_length.setter
    def content_length(self, value):
        util.typeCheck(value, long, "content_length") 
        util.valueCheck(value, 0, 1E10, "content_length")
        self._content_length = value

    @property
    def uri(self):
        """The URI corresponding to the Artifact.

        should be a well formed uri thatenables the archive system to
        find the file.

        eg.: 
           ad:CFHT/1234567p
           vos://cadc.nrc.ca!vospace/APASS/apass_north/proc/100605/n100605.0384.new.fz

        """
        return self._uri

    @uri.setter
    def uri(self, value):
        util.typeCheck(value, str, 'uri')
        uri = urlparse(value).geturl()
        util.valueCheck(value, None, None, 'uri', override=uri)
        self._uri = uri
        

    @property
    def product_type(self):
        """The product type associated with the Artifact. 

        type:  caom2.ProductType
        restricted to caom2.ProductType.names()

        eg.  Artifcat.product_type = caom2.ProductType('SCIENCE')

        """
        
        if self._product_type == None:
            return self._compute_product_type()
        return self._product_type

    def _compute_product_type(self):
        """Loop over the parts and if they all have the same product_type then
        return that, else None

        """
        common = None
        for part in self.parts.itervalues():
            if common is None:
                common = part.product_type
            else:
                if common != part.product_type:
                    common = None
                    break

        return common

    @product_type.setter
    def product_type(self, value):
        if value is not None:
            assert isinstance(value, ProductType), (
                    "product type is not a ProductType: {0}".format(value))
        self._product_type = value

    @property
    def alternative(self):
        """This artifact is an alternative representation.

        type: bool  (True/False)

        """
        return self._alternative

    @alternative.setter
    def alternative(self, value):
        if not isinstance(value,bool):
            raise TypeError(
                "alternative must be boolean. received: {0}".format(value))
        self._alternative = value

    @property
    def parts(self):
        """A list of caom2.Part objects.

        This is currently implemented as caom2.util.TypedOrderedDict
        which isn't really a list or a dictionary.

        see caom2.Part for help on making a Part object to add to the
        list.

        Usage example:

        a = Artifact('ad:CFHT/123456p')
        a.parts.add(Part('partName')) # to add a new part to an
                                      # existing artifact.

        print a.parts.items()   #  a list of part names
        for part in parts:  # to iterate over all parts in an Artifact

        """
        return self._parts

    @parts.setter
    def parts(self, value):
        util.typeCheck(value, 
                        TypedOrderedDict, 
                        'parts', 
                        override=False)
        self._parts = value
