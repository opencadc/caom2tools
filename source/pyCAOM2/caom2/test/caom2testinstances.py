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

""" Defines Caom2TestInstances class """

import collections
from datetime import datetime

from caom2.caom2_algorithm import Algorithm
from caom2.caom2_artifact import Artifact
from caom2.caom2_chunk import Chunk
from caom2.caom2_composite_observation import CompositeObservation
from caom2.caom2_data_quality import DataQuality
from caom2.caom2_energy_transition import EnergyTransition
from caom2.caom2_enums import CalibrationLevel
from caom2.caom2_enums import DataProductType
from caom2.caom2_enums import ObservationIntentType
from caom2.caom2_enums import ProductType
from caom2.caom2_enums import Quality
from caom2.caom2_enums import Status
from caom2.caom2_enums import TargetType
from caom2.caom2_enums import ReleaseType
from caom2.caom2_environment import Environment
from caom2.caom2_instrument import Instrument
from caom2.caom2_metrics import Metrics
from caom2.caom2_observation_uri import ObservationURI
from caom2.caom2_part import Part
from caom2.caom2_plane import Plane
from caom2.caom2_plane_uri import PlaneURI
from caom2.caom2_proposal import Proposal
from caom2.caom2_provenance import Provenance
from caom2.caom2_requirements import Requirements
from caom2.caom2_simple_observation import SimpleObservation
from caom2.caom2_target import Target
from caom2.caom2_target_position import TargetPosition
from caom2.caom2_telescope import Telescope
from caom2.types.caom2_point import Point
from caom2.util.caom2_util import TypedList, TypedSet
from caom2.wcs.caom2_axis import Axis
from caom2.wcs.caom2_coord2d import Coord2D
from caom2.wcs.caom2_coord_axis1d import CoordAxis1D
from caom2.wcs.caom2_coord_axis2d import CoordAxis2D
from caom2.wcs.caom2_coord_bounds1d import CoordBounds1D
from caom2.wcs.caom2_coord_circle2d import CoordCircle2D
from caom2.wcs.caom2_coord_error import CoordError
from caom2.wcs.caom2_coord_function1d import CoordFunction1D
from caom2.wcs.caom2_coord_function2d import CoordFunction2D
from caom2.wcs.caom2_coord_polygon2d import CoordPolygon2D
from caom2.wcs.caom2_coord_range1d import CoordRange1D
from caom2.wcs.caom2_coord_range2d import CoordRange2D
from caom2.wcs.caom2_dimension2d import Dimension2D
from caom2.wcs.caom2_observable_axis import ObservableAxis
from caom2.wcs.caom2_polarization_wcs import PolarizationWCS
from caom2.wcs.caom2_ref_coord import RefCoord
from caom2.wcs.caom2_slice import Slice
from caom2.wcs.caom2_spatial_wcs import SpatialWCS
from caom2.wcs.caom2_spectral_wcs import SpectralWCS
from caom2.wcs.caom2_temporal_wcs import TemporalWCS
from caom2.wcs.caom2_value_coord2d import ValueCoord2D


class Caom2TestInstances(object):

    _collection = "collection"
    _observation_id = "observationID"
    _product_id = "productId"
    _keywords = TypedList(str, "keyword1", "keyword2")
    _ivoa_date = datetime(2012, 07, 11, 13, 26, 37, 0)

    def __init__(self):
        self.depth = 5
        self.complete = True
        self.bounds_is_circle = True
        self.caom_version = 20

    @property
    def depth(self):
        return self._depth

    @depth.setter
    def depth(self, v):
        self._depth = v

    @property
    def complete(self):
        return self._complete

    @complete.setter
    def complete(self, v):
        self._complete = v

    @property
    def bounds_is_circle(self):
        return self._bounds_is_circle

    @bounds_is_circle.setter
    def bounds_is_circle(self, v):
        self._bounds_is_circle = v

    @property
    def caom_version(self):
        return self._caom_version

    @caom_version.setter
    def caom_version(self, v):
        self._caom_version = v

    def get_simple_observation(self):
        observation = SimpleObservation(Caom2TestInstances._collection,
                                        Caom2TestInstances._observation_id)
        if self.complete:
            observation.sequence_number = int(5)
            observation.obs_type = "flat"
            observation.intent = (
                ObservationIntentType.CALIBRATION)
            observation.meta_release = Caom2TestInstances._ivoa_date
            observation.proposal = self.get_proposal()
            observation.target = self.get_target()
            observation.target_position = self.get_target_position()
            if self.caom_version == 21:
                observation.requirements = self.get_requirements()
            observation.telescope = self.get_telescope()
            observation.instrument = self.get_instrument()
            observation.environment = self.get_environment()
        if self.depth > 1:
            observation.planes.update(self.get_planes())
        return observation

    def get_composite_observation(self):
        observation = CompositeObservation(Caom2TestInstances._collection,
                                           Caom2TestInstances._observation_id,
                                        self.get_algorithm())
        if self.complete:
            observation.sequence_number = int(10)
            observation.obs_type = "filed"
            observation.intent = (
                ObservationIntentType.SCIENCE)
            observation.meta_release = Caom2TestInstances._ivoa_date
            observation.proposal = self.get_proposal()
            observation.target = self.get_target()
            observation.target_position = self.get_target_position()
            if self.caom_version == 21:
                observation.requirements = self.get_requirements()
            observation.telescope = self.get_telescope()
            observation.instrument = self.get_instrument()
            observation.environment = self.get_environment()
        if self.depth > 1:
            observation.planes.update(self.get_planes())
            observation.members.update(self.get_members())
        return observation

    def get_algorithm(self):
        return Algorithm("algorithmName")

    def get_proposal(self):
        proposal = Proposal("proposalId")
        proposal.pi_name = "proposalPi"
        proposal.project = "proposalProject"
        proposal.title = "proposalTitle"
        proposal.keywords.extend(Caom2TestInstances._keywords)
        return proposal

    def get_target(self):
        target = Target("targetName")
        target.target_type = TargetType.OBJECT
        target.standard = False
        target.redshift = 1.5
        target.keywords.extend(Caom2TestInstances._keywords)
        return target

    def get_target_position(self):
        point = Point(1.0, 2.0)
        target_position = TargetPosition(point, "coordsys")
        target_position.equinox = 3.0
        return target_position

    def get_requirements(self):
        return Requirements(Status.FAIL)

    def get_telescope(self):
        telescope = Telescope("telescopeName")
        telescope.geo_location_x = 1.0
        telescope.geo_location_y = 2.0
        telescope.geo_location_z = 3.0
        telescope.keywords.extend(Caom2TestInstances._keywords)
        return telescope

    def get_instrument(self):
        instrument = Instrument("instrumentName")
        instrument.keywords.extend(Caom2TestInstances._keywords)
        return instrument

    def get_environment(self):
        env = Environment()
        env.seeing = 0.08
        env.humidity = 0.35
        env.elevation = 2.7
        env.tau = 0.7
        env.wavelength_tau = 450e-6
        env.ambient_temp = 20.0
        env.photometric = True
        return env

    def get_members(self):
        members = TypedSet(
            ObservationURI, ObservationURI("caom:foo/bar"))
        return members

    def get_planes(self):
        planes = collections.OrderedDict()
        plane = Plane("productID")
        if self.complete:
            plane.meta_release = Caom2TestInstances._ivoa_date
            plane.data_release = Caom2TestInstances._ivoa_date
            plane.data_product_type = DataProductType.IMAGE
            plane.calibration_level = CalibrationLevel.PRODUCT
            plane.provenance = self.get_provenance()
            plane.metrics = self.get_metrics()
            if self.caom_version == 21:
                plane.quality = self.get_quality()

        if self.depth > 2:
            for k, v in self.get_artifacts().iteritems():
                plane.artifacts[k] = v
        planes["productID"] = plane
        return planes

    def get_provenance(self):
        provenance = Provenance("name")
        provenance.version = "version"
        provenance.product = "product"
        provenance.producer = "producer"
        provenance.run_id = "run_id"
        provenance.reference = "http://foo/bar"
        provenance.last_executed = Caom2TestInstances._ivoa_date
        provenance.keywords.extend(Caom2TestInstances._keywords)
        provenance.inputs.update(self.get_inputs())
        return provenance

    def get_inputs(self):
        return TypedSet(PlaneURI, PlaneURI("caom:foo/bar/plane1"),
                        PlaneURI("caom:foo/bar/plane2"))

    def get_metrics(self):
        metrics = Metrics()
        metrics.source_number_density = float(1.0)
        metrics.background = float(2.0)
        metrics.background_std_dev = float(3.0)
        metrics.flux_density_limit = float(4.0)
        metrics.mag_limit = float(5.0)
        return metrics

    def get_quality(self):
        return DataQuality(Quality.JUNK)

    def get_artifacts(self):
        artifacts = collections.OrderedDict()
        artifact = Artifact("ad:foo/bar1", ProductType.SCIENCE, ReleaseType.META)
        if self.complete:
            artifact.content_type = "application/fits"
            artifact.content_length = 12345L
        if self.depth > 3:
            for k, v in self.get_parts().iteritems():
                artifact.parts[k] = v
        artifacts["ad:foo/bar1"] = artifact
        return artifacts

    def get_parts(self):
        parts = collections.OrderedDict()
        part = Part("x")
        if self.complete:
            part.product_type = ProductType.SCIENCE
        if self.depth > 4:
            part.chunk = self.get_chunk()
        parts["x"] = part
        return parts

    def get_chunk(self):
        chunk = Chunk()
        if self.complete:
            chunk.naxis = 5
            chunk.observable_axis = 1
            chunk.position_axis_1 = 1
            chunk.position_axis_2 = 2
            chunk.energy_axis = 3
            chunk.time_axis = 4
            chunk.polarization_axis = 5
            chunk.observable = self.get_observable_axis()
            chunk.position = self.get_spatial_wcs()
            chunk.energy = self.get_spectral_wcs()
            chunk.time = self.get_temporal_wcs()
            chunk.polarization = self.get_polarization_wcs()
        return chunk

    def get_observable_axis(self):
        observable = ObservableAxis(self.get_slice())
        if self.complete:
            observable.independent = self.get_slice()
        return observable

    def get_spatial_wcs(self):
        coord_axis_2d = self.get_coord_axis_2d()
        position = SpatialWCS(coord_axis_2d)
        if self.complete:
            position.coordsys = "position coordsys"
            position.equinox = 2000.0
            position.resolution = 0.5
        return position

    def get_spectral_wcs(self):
        axis = self.get_coord_axis_1d()
        energy = SpectralWCS(axis, "energy specsys")
        if self.complete:
            energy.ssysobs = "energy ssysobs"
            energy.ssyssrc = "energy ssyssrc"
            energy.restfrq = 1.0
            energy.restwav = 2.0
            energy.velosys = 3.0
            energy.zsource = 4.0
            energy.velang = 5.0
            energy.bandpassName = "energy bandpassName"
            energy.resolvingPower = 6.0
            energy.transition = EnergyTransition("H", "21cm")
        return energy

    def get_temporal_wcs(self):
        axis = self.get_coord_axis_1d()
        time = TemporalWCS(axis)
        if self.complete:
            time.exposure = 1.0
            time.resolution = 2.0
            time.timesys = "UTC"
            time.trefpos = "TOPOCENTER"
            time.mjdref = 3.0
        return time

    def get_polarization_wcs(self):
        axis = Axis('STOKES')
        axis_1d = CoordAxis1D(axis)
        #IQUV
        axis_1d.function = CoordFunction1D(4L, 1.0, RefCoord(1.0, 1.0))
        pol = PolarizationWCS(axis_1d)
        return pol

    def get_slice(self):
        return Slice(Axis("sliceCtype", "sliceCunit"), 1L)

    def get_coord_axis_1d(self):
        coord_axis_1d = CoordAxis1D(Axis("axisCtype", "axisCunit"))
        if self.complete:
            coord_axis_1d.error = CoordError(1.0, 1.5)
            coord_axis_1d.range = CoordRange1D(RefCoord(2.0, 2.5),
                                               RefCoord(3.0, 3.5))
            coord_axis_1d.function = (
                CoordFunction1D(4L, 4.5, RefCoord(5.0, 5.5)))
            bounds = CoordBounds1D()
            bounds.samples.append(CoordRange1D(RefCoord(6.0, 6.5),
                                               RefCoord(7.0, 7.5)))
            bounds.samples.append(CoordRange1D(RefCoord(8.0, 8.5),
                                               RefCoord(9.0, 9.5)))
            coord_axis_1d.bounds = bounds
        return coord_axis_1d

    def get_coord_axis_2d(self):
        axis1 = Axis("axis1Ctype", "axis1Cunit")
        axis2 = Axis("axis2Ctype", "axis2Cunit")
        coord_axis_2d = CoordAxis2D(axis1, axis2)
        if self.complete:
            coord_axis_2d.error1 = CoordError(1.0, 1.5)
            coord_axis_2d.error2 = CoordError(2.0, 2.5)
            start = Coord2D(RefCoord(3.0, 3.5), RefCoord(4.0, 4.5))
            end = Coord2D(RefCoord(5.0, 5.5), RefCoord(6.0, 6.5))
            coord_axis_2d.range = CoordRange2D(start, end)
            dimension = Dimension2D(7L, 8L)
            ref_coord = Coord2D(RefCoord(9.0, 9.5), RefCoord(10.0, 10.5))
            coord_axis_2d.function = (CoordFunction2D(dimension, ref_coord,
                                      11.0, 12.0, 13.0, 14.0))
            if self.bounds_is_circle:
                center = ValueCoord2D(15.0, 16.0)
                coord_axis_2d.bounds = CoordCircle2D(center, 17.0)
            else:
                polygon = CoordPolygon2D()
                polygon.vertices.append(ValueCoord2D(15.0, 16.0))
                polygon.vertices.append(ValueCoord2D(17.0, 18.0))
                polygon.vertices.append(ValueCoord2D(19.0, 20.0))
                coord_axis_2d.bounds = polygon
        return coord_axis_2d
