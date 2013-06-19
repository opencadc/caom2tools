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

"""contains the definition of the caom2.Instrument object"""

from caom2_object import Caom2Object
from util.caom2_util import TypedList
import util.caom2_util as util


class Instrument(Caom2Object):
    """The telescopic instrument that recorded a given observation.

    Each observation should have an associated instrument.

    eg:
    inst=caom2.Observation('CFHT','123456p','exposure').Instrument("MEGAPRIME")

    Each instrument can have list of keywords.
    eg.
    inst.keywords.append("shutter=closed")
    """

    def __init__(self, name):
        """
        Initializes a Instrument instance

        Arguments:
        name - name of the instrument
        """
        util.typeCheck(name, str, 'name', override='none')
        self._name = name
        self._keywords = TypedList((str),)

    # Properties
    @property
    def name(self):
        """The name of the instrument.

        type: str
        """
        return self._name

    @property
    def keywords(self):
        """A list of strings that are keywords associated with the instrument.

        eg.  keywords.append(("ccd=off","hot","shutter broken"))

        Keywords are stored as text in the database and are searched as text,
        so in the above example you can search for ccd=off but that is just
        matching the text 'ccd=off' and you can not do 'ccd!=on' to find that
        key/value pair.

        Also, the list is concated together and inserted into the
        model as a single text field, so keywords like 'shutter
        broken' will match 'shutter' and 'broken'.  If these keywords
        appear in a pick list 'shutter' and 'broken' will be seperate
        entries.  So shutter_broken is a better choice.
        """
        return self._keywords
