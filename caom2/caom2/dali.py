# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2025.                            (c) 2025.
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

from . import common
from . import caom_util

__all__ = ['Interval']


class Interval(common.CaomObject):
    def __init__(self, lower, upper):

        super().__init__()
        self.lower = lower
        self.upper = upper

    def get_width(self):
        return self._upper - self._lower

    @classmethod
    def intersection(cls, i1, i2):
        if i1.lower > i2.upper or i1.upper < i2.lower:
            return None

        lb = max(i1.lower, i2.lower)
        ub = min(i1.upper, i2.upper)
        return cls(lb, ub)

    # Properties

    @property
    def lower(self):
        """
        type: float
        """
        return self._lower

    @lower.setter
    def lower(self, value):
        caom_util.type_check(value, float, 'lower', override=False)
        has_upper = True
        try:
            self._upper
        except AttributeError:
            has_upper = False
        if has_upper and self._upper < value:
            raise ValueError("Interval: attempt to set upper < lower "
                             "for {}, {}".format(self._upper, value))
        self._lower = value

    @property
    def upper(self):
        """
        type: float
        """
        return self._upper

    @upper.setter
    def upper(self, value):
        caom_util.type_check(value, float, 'upper', override=False)
        has_lower = True
        try:
            self._lower
        except AttributeError:
            has_lower = False
        if has_lower and value < self._lower:
            raise ValueError("Interval: attempt to set upper < lower "
                             "for {}, {}".format(value, self._lower))
        self._upper = value

    # def validate(self):
    #     """
    #     Performs a validation of the current object.
    #
    #     An AssertionError is thrown if the object does not represent an
    #     Interval
    #     """
    #     if self._samples is not None:
    #
    #         if len(self._samples) == 0:
    #             raise ValueError(
    #                 'invalid interval (samples cannot be empty)')
    #
    #         prev = None
    #         for sample in self._samples:
    #             if sample.lower < self._lower:
    #                 raise ValueError(
    #                     'invalid interval: sample extends below lower bound: '
    #                     '{} vs {}'.format(sample, self._lower))
    #             if sample.upper > self._upper:
    #                 raise ValueError(
    #                     'invalid interval: sample extends above upper bound: '
    #                     '{} vs {}'.format(sample, self._upper))
    #             if prev is not None:
    #                 if sample.lower <= prev.upper:
    #                     raise ValueError(
    #                         'invalid interval: sample overlaps previous '
    #                         'sample:\n{}\nvs\n{}'.format(sample, prev))
    #             prev = sample
