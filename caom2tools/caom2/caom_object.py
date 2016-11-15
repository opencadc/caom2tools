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

import inspect
import random
import time
import uuid
from datetime import datetime

from caom_util import Util


class CaomObject(object):
    """
    setup all objects with the same generic equality, str and repr methods
    """

    def __init__(self):
        pass

    def __str__(self):
        args = inspect.getargspec(self.__init__).args[1:]
        class_name = self.__class__.__name__
        return "\n".join(["{}.{} : {}".
                         format(class_name, arg, getattr(self, arg, None))
                         for arg in args])

    def __eq__(self, other):
        if type(other) == type(self):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __repr__(self):
        args = inspect.getargspec(self.__init__).args[1:]
        class_name = ""
        if self.__class__.__module__ != '__main__':
            class_name += self.__class__.__module__ + "."
        class_name += self.__class__.__name__
        pading = " " * (len(class_name) + 1)
        return class_name + "(" + (
            ",\n" + pading).join(
            ["%s=%r" % (arg, getattr(self, arg, None)
                        ) for arg in args]) + ")"


class AbstractCaomEntity(CaomObject):
    """Class that defines the persistence unique ID and last mod date """

    def __init__(self, fulluuid=False):
        self._id = AbstractCaomEntity._gen_id(fulluuid)
        self._last_modified = AbstractCaomEntity._gen_last_modified()

    @classmethod
    def _gen_id(cls, fulluuid=False):
        """Generate a 128 but UUID by default. For backwards compatibility
        allow creation of a 64 bit UUID using a rand number for the
        lower 64 bits. First two bytes of the random number are generated
        with the random and the last 6 bytes from the current time
        in microseconds.

        return: UUID
        """

        if fulluuid:
            return uuid.uuid4()
        else:
            vmrandom = random.randint(-int(0x7fff), int(0x7fff)) << 8 * 6
            randtime = int(round(time.time() * 1000000))
            randtime = randtime & 0xffffffffffff
            rand = vmrandom | randtime
            if rand & 0x8000000000000000:
                rand = 0x1000000000000000 + rand
            return Util.long2uuid(rand)

    @classmethod
    def _gen_last_modified(cls):
        """Generate a datetime with 3 digit microsecond precision.

        return: datatime
            IVOA date format to millisecond precision.
        """
        now = datetime.now()
        return datetime(now.year, now.month, now.day, now.hour, now.minute, \
                        now.second, long(str(now.microsecond)[:-3] + '000'))
