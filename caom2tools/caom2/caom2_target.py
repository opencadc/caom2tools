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

"""defines the caom2.Target class"""

import util.caom2_util as util
from caom2_enums import TargetType
from caom2_object import Caom2Object


class Target(Caom2Object):
    """ Target """

    def __init__(self, name,
                 target_type=None,
                 standard=None,
                 redshift=None,
                 keywords=None,
                 moving=None):
        """
        Initializes a Target instance

        Arguments:
        name : name of the target
        type : type of the target
        """

        self.name = name
        self.target_type = target_type
        self.standard = standard
        self.redshift = redshift
        if keywords is None:
            keywords = set()
        self.keywords = keywords
        self.moving = moving

    # Properties

    @property
    def name(self):
        """A name for the target

        eg 'NGC 3115'
        type: str

        """
        return self._name

    @name.setter
    def name(self, value):
        util.typeCheck(value, str, "name", override=False)
        self._name = value

    @property
    def target_type(self):
        """A keyword describing the type of target.

        must be from the list
        """ + ",".join(TargetType.names()) + """
        type: TargetType

        """
        return self._type

    @target_type.setter
    def target_type(self, value):
        if isinstance(value, str):
            value = TargetType(value)
        util.typeCheck(value, TargetType, "target_type")
        self._type = value

    @property
    def keywords(self):
        """A set of keywords associated with this target.

        eg. keywords.add('galaxy')
        type: set

        """
        return self._keywords

    @keywords.setter
    def keywords(self, value):
        util.typeCheck(value, set, 'keywords', override=False)
        self._keywords = value

    @property
    def standard(self):
        """Is this a standard field?

        eg True
        type: bool

        """
        return self._standard

    @standard.setter
    def standard(self, value):
        util.typeCheck(value, bool, 'standard')
        self._standard = value

    @property
    def redshift(self):
        """The redshift of the observed target.

        eg 1.2  (can be None)
        type: float

        """
        return self._redshift

    @redshift.setter
    def redshift(self, value):
        util.typeCheck(value, float, 'redshift')
        util.valueCheck(value, -0.5, 1200, 'redshift')
        self._redshift = value

    @property
    def moving(self):
        """Is this a moving target?

        eg True
        type: bool

        """
        return self._moving

    @moving.setter
    def moving(self, value):
        util.typeCheck(value, bool, 'moving')
        self._moving = value
