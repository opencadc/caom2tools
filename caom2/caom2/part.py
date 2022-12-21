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

"""Defines the caom2.Part class which describes
the caom2_Observation_Plane_Artifact_Part object."""

from builtins import str

from . import caom_util
from .chunk import Chunk
from .chunk import ProductType
from .common import AbstractCaomEntity

__all__ = ['Part']


class Part(AbstractCaomEntity):
    """A qualitative subsection of an artifact.
       eg: a extension of a FITS file.


       This object should contain the product_tpye attribute
       and the list of chunks.
    """

    def __init__(self, name, product_type=None, chunks=None):
        super(Part, self).__init__()
        self.name = name
        self.product_type = product_type
        if chunks is None:
            chunks = caom_util.TypedList(Chunk, )
        self.chunks = chunks

    def _key(self):
        return self.name

    def __eq__(self, y):
        if isinstance(y, Part):
            return self._key() == y._key()
        return False

    def __hash__(self):
        return hash(self._key())

    @property
    def product_type(self):
        """The type of data product referred to by this part.

        Must be one of the allowed data product types:
        caom2.ProductType"""
        return self._product_type

    @product_type.setter
    def product_type(self, value):
        caom_util.type_check(value, ProductType, "product_type")
        self._product_type = value

    @property
    def name(self):
        """The name of this part, normally the FITS extension.

        This values is also used as the key to find the part in the
        Artifact.parts dictionary"""
        return self._name

    @name.setter
    def name(self, value):
        caom_util.type_check(value, str, 'name', override=False)
        self._name = value

    @property
    def chunks(self):
        """A list of chunks that this part contains"""
        return self._chunks

    @chunks.setter
    def chunks(self, value):
        caom_util.type_check(value, caom_util.TypedList, 'chunks',
                             override=False)
        self._chunks = value
