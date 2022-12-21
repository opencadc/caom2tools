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

""" Defines Caom2TestInstances class """

import collections
from datetime import datetime
import uuid
from builtins import int

from caom2 import artifact
from caom2 import caom_util
from caom2 import chunk
from caom2 import common
from caom2 import observation
from caom2 import part
from caom2 import plane
from caom2 import shape
from caom2 import wcs


def fix_ids(observation):
    # temporary function to make the internal ids 64 bit UUID instead of
    # the default 128. This is only required for caom2.0
    def get_64bit_uuid(id):
        return uuid.UUID(fields=(0x00000000, 0x0000, 0x0000,
                         id.clock_seq_hi_variant,
                         id.clock_seq_low, id.node))

    observation._id = get_64bit_uuid(observation._id)
    for p in observation.planes.values():
        p._id = get_64bit_uuid(p._id)
        for a in p.artifacts.values():
            a._id = get_64bit_uuid(a._id)
            for pa in a.parts.values():
                pa._id = get_64bit_uuid(pa._id)
                for c in pa.chunks:
                    c._id = get_64bit_uuid(c._id)


class Caom2TestInstances(object):
    _collection = "collection"
    _observation_id = "observationID"
    _product_id = "productId"
    _keywords = {"keyword1", "keyword2"}
    _ivoa_date = datetime(2012, 7, 11, 13, 26, 37, 0)

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

    def get_simple_observation(self, short_uuid=False):
        simple_observation = \
            observation.SimpleObservation(Caom2TestInstances._collection,
                                          Caom2TestInstances._observation_id)
        if self.complete:
            simple_observation.sequence_number = int(5)
            simple_observation.obs_type = "flat"
            simple_observation.intent =\
                observation.ObservationIntentType.CALIBRATION
            simple_observation.meta_release = Caom2TestInstances._ivoa_date
            simple_observation.proposal = self.get_proposal()
            simple_observation.target = self.get_target()
            simple_observation.target_position = self.get_target_position()
            simple_observation.telescope = self.get_telescope()
            simple_observation.instrument = self.get_instrument()
            simple_observation.environment = self.get_environment()
            simple_observation.last_modified = common.get_current_ivoa_time()
            if self.caom_version >= 23:
                simple_observation.max_last_modified =\
                    common.get_current_ivoa_time()
                simple_observation.meta_checksum = common.ChecksumURI(
                    "md5:9882dbbf9cadc221019b712fd402bcbd")
                simple_observation.acc_meta_checksum = common.ChecksumURI(
                    "md5:844ce247db0844ad9f721430c80e7a21")
            if self.caom_version >= 24:
                simple_observation.meta_read_groups.add(
                    "ivo://cadc.nrc.ca/groups?A")
                simple_observation.meta_read_groups.add(
                    "ivo://cadc.nrc.ca/groups?B")
        if self.depth > 1:
            simple_observation.planes.update(self.get_planes())
        if self.caom_version == 20 or short_uuid:
            fix_ids(simple_observation)
        return simple_observation

    def get_composite_observation(self, short_uuid=False):
        composite_observation = \
            observation.CompositeObservation(
                Caom2TestInstances._collection,
                Caom2TestInstances._observation_id,
                self.get_algorithm())
        print("Creating test composite observation of version " + str(
            self.caom_version))
        if self.complete:
            composite_observation.sequence_number = int(10)
            composite_observation.obs_type = "filed"
            composite_observation.intent =\
                observation.ObservationIntentType.SCIENCE
            composite_observation.meta_release = Caom2TestInstances._ivoa_date
            composite_observation.proposal = self.get_proposal()
            composite_observation.target = self.get_target()
            composite_observation.target_position = self.get_target_position()
            composite_observation.telescope = self.get_telescope()
            composite_observation.instrument = self.get_instrument()
            composite_observation.environment = self.get_environment()
            composite_observation.last_modified =\
                common.get_current_ivoa_time()
            if self.caom_version >= 23:
                composite_observation.max_last_modified = \
                    common.get_current_ivoa_time()
                composite_observation.meta_checksum = common.ChecksumURI(
                    "md5:9882dbbf9cadc221019b712fd402bcbd")
                composite_observation.acc_meta_checksum = common.ChecksumURI(
                    "md5:844ce247db0844ad9f721430c80e7a21")
        if self.depth > 1:
            composite_observation.planes.update(self.get_planes())
            composite_observation.members.update(self.get_members())
        return composite_observation

    def get_derived_observation(self, short_uuid=False):
        derived_observation = \
            observation.DerivedObservation(
                Caom2TestInstances._collection,
                Caom2TestInstances._observation_id,
                self.get_algorithm())
        print("Creating test composite observation of version " + str(
            self.caom_version))
        if self.complete:
            derived_observation.sequence_number = int(10)
            derived_observation.obs_type = "filed"
            derived_observation.intent =\
                observation.ObservationIntentType.SCIENCE
            derived_observation.meta_release = Caom2TestInstances._ivoa_date
            derived_observation.proposal = self.get_proposal()
            derived_observation.target = self.get_target()
            derived_observation.target_position = self.get_target_position()
            derived_observation.telescope = self.get_telescope()
            derived_observation.instrument = self.get_instrument()
            derived_observation.environment = self.get_environment()
            derived_observation.last_modified =\
                common.get_current_ivoa_time()
            derived_observation.max_last_modified = \
                common.get_current_ivoa_time()
            derived_observation.meta_checksum = common.ChecksumURI(
                "md5:9882dbbf9cadc221019b712fd402bcbd")
            derived_observation.acc_meta_checksum = common.ChecksumURI(
                "md5:844ce247db0844ad9f721430c80e7a21")
            derived_observation.meta_read_groups.add(
                "ivo://cadc.nrc.ca/groups?A")
            derived_observation.meta_read_groups.add(
                "ivo://cadc.nrc.ca/groups?B")
        if self.depth > 1:
            derived_observation.planes.update(self.get_planes())
            derived_observation.members.update(self.get_members())
        return derived_observation

    def get_algorithm(self):
        return observation.Algorithm("algorithmName")

    def get_proposal(self):
        proposal = observation.Proposal("proposalId")
        proposal.pi_name = "proposalPi"
        proposal.project = "proposalProject"
        proposal.title = "proposalTitle"
        proposal.keywords.update(Caom2TestInstances._keywords)
        return proposal

    def get_target(self):
        target = observation.Target("targetName")
        if self.caom_version >= 24:
            target.target_id = 'someTargetID'
        target.target_type = observation.TargetType.OBJECT
        target.standard = False
        target.redshift = 1.5
        target.keywords.update(Caom2TestInstances._keywords)
        return target

    def get_target_position(self):
        point = shape.Point(1.0, 2.0)
        target_position = observation.TargetPosition(point, "coordsys")
        target_position.equinox = 3.0
        return target_position

    def get_requirements(self):
        return observation.Requirements(observation.Status.FAIL)

    def get_telescope(self):
        telescope = observation.Telescope("telescopeName")
        telescope.geo_location_x = 1.0
        telescope.geo_location_y = 2.0
        telescope.geo_location_z = 3.0
        telescope.keywords.update(Caom2TestInstances._keywords)
        return telescope

    def get_instrument(self):
        instrument = observation.Instrument("instrumentName")
        instrument.keywords.update(Caom2TestInstances._keywords)
        return instrument

    def get_environment(self):
        env = observation.Environment()
        env.seeing = 0.08
        env.humidity = 0.35
        env.elevation = 2.7
        env.tau = 0.7
        env.wavelength_tau = 450e-6
        env.ambient_temp = 20.0
        env.photometric = True
        return env

    def get_members(self):
        members = caom_util.TypedSet(
            observation.ObservationURI,
            observation.ObservationURI("caom:foo/bar"))
        return members

    def get_planes(self):
        planes = collections.OrderedDict()
        if self.caom_version < 22:
            shapes = ['']
        elif self.caom_version >= 23:
            shapes = ['polygon', 'circle']
        else:
            shapes = ['polygon']

        for s in shapes:
            prod_id = "productID{}".format(s)
            _plane = plane.Plane(prod_id)
            if self.complete:
                _plane.meta_release = Caom2TestInstances._ivoa_date
                _plane.data_release = Caom2TestInstances._ivoa_date
                _plane.data_product_type = plane.DataProductType.IMAGE
                _plane.calibration_level = plane.CalibrationLevel.PRODUCT
                _plane.provenance = self.get_provenance()
                _plane.metrics = self.get_metrics()
                _plane.last_modified = common.get_current_ivoa_time()
                if self.caom_version >= 23:
                    _plane.creator_id = "ivo://cadc.nrc.ca?testuser"
                    _plane.max_last_modified = common.get_current_ivoa_time()
                    _plane.meta_checksum = common.ChecksumURI(
                        "md5:9882dbbf9cadc221019b712fd402bcbd")
                    _plane.acc_meta_checksum = common.ChecksumURI(
                        "md5:844ce247db0844ad9f721430c80e7a21")
                if s == 'polygon':
                    _plane.position = self.get_poly_position()
                if s == 'circle':
                    _plane.position = self.get_circle_position()
                _plane.energy = self.get_energy()
                _plane.time = self.get_time()
                _plane.polarization = self.get_polarization()
                if self.caom_version >= 24:
                    _plane.custom_axis = self.get_custom()
                    _plane.meta_read_groups.add('ivo://cadc.nrc.ca/groups?A')
                    _plane.meta_read_groups.add('ivo://cadc.nrc.ca/groups?D')
                    _plane.data_read_groups.add('ivo://cadc.nrc.ca/groups?B')
                    _plane.data_read_groups.add('ivo://cadc.nrc.ca/groups?C')

            if self.depth > 2:
                for k, v in self.get_artifacts().items():
                    _plane.artifacts[k] = v
            planes[prod_id] = _plane
        return planes

    def get_poly_position(self):
        position = plane.Position()

        if self.caom_version >= 23:
            v0 = shape.Vertex(0.0, 0.0, shape.SegmentType.MOVE)
            v1 = shape.Vertex(3.0, 4.0, shape.SegmentType.LINE)
            v2 = shape.Vertex(2.0, 3.0, shape.SegmentType.LINE)
            v3 = shape.Vertex(1.0, 2.0, shape.SegmentType.LINE)
            v4 = shape.Vertex(0.0, 0.0, shape.SegmentType.CLOSE)
            vl = [v0, v1, v2, v3, v4]

            samples = shape.MultiPolygon(vertices=vl)

            p1 = shape.Point(0.0, 0.0)
            p2 = shape.Point(3.0, 4.0)
            p3 = shape.Point(2.0, 3.0)
            p4 = shape.Point(1.0, 2.0)
            p = [p1, p2, p3, p4]
            polygon = shape.Polygon(points=p, samples=samples)

            position.bounds = polygon

        position.dimension = wcs.Dimension2D(10, 20)
        position.resolution = 0.5
        position.sample_size = 1.1
        position.time_dependent = False
        if self.caom_version >= 24:
            position.resolution_bounds = shape.Interval(1.0, 2.0)

        return position

    def get_circle_position(self):
        position = plane.Position()
        position.bounds = shape.Circle(shape.Point(1.1, 2.2), 3.0)
        position.dimension = wcs.Dimension2D(10, 20)
        position.resolution = 0.5
        position.sample_size = 1.1
        position.time_dependent = False
        if self.caom_version >= 24:
            position.resolution_bounds = shape.Interval(1.0, 2.0)
        return position

    def get_energy(self):
        energy = plane.Energy()

        lower = 1.0
        upper = 2.0
        lower1 = 1.1
        upper1 = 2.1
        lower2 = 1.2
        upper2 = 2.2
        samples = [shape.SubInterval(lower, lower1),
                   shape.SubInterval(lower2, upper),
                   shape.SubInterval(upper1, upper2)]

        interval = shape.Interval(lower, upper2, samples)

        energy.bounds = interval
        energy.dimension = 100
        energy.resolving_power = 2.0
        energy.sample_size = 1.1
        energy.bandpass_name = "e"
        if self.caom_version >= 24:
            energy.energy_bands.add(plane.EnergyBand.GAMMARAY)
            energy.energy_bands.add(plane.EnergyBand.OPTICAL)
        energy.transition = wcs.EnergyTransition("species", "transition")

        return energy

    def get_time(self):
        time = plane.Time()

        lower = 1.0
        upper = 2.0
        lower1 = 1.1
        upper1 = 2.1
        lower2 = 1.2
        upper2 = 2.2
        samples = [shape.SubInterval(lower, lower1),
                   shape.SubInterval(lower2, upper),
                   shape.SubInterval(upper1, upper2)]

        interval = shape.Interval(lower, upper2, samples)

        time.bounds = interval
        if self.caom_version >= 24:
            time.resolution_bounds = shape.Interval(22.2, 33.3)
        time.dimension = 1
        time.resolution = 2.1
        time.sample_size = 3.0
        time.exposure = 10.3

        return time

    def get_custom(self):
        custom = plane.CustomAxis('MyAxis')
        custom.bounds = shape.Interval(2.2, 3.3)
        custom.dimension = 1

    def get_polarization(self):
        polarization = plane.Polarization()

        p_states = [plane.PolarizationState.LL, plane.PolarizationState.XY]

        polarization.dimension = 2
        polarization.polarization_states = p_states

        return polarization

    def get_provenance(self):
        provenance = plane.Provenance("name")
        provenance.version = "version"
        provenance.product = "product"
        provenance.producer = "producer"
        provenance.run_id = "run_id"
        provenance.reference = "http://foo/bar"
        provenance.last_executed = Caom2TestInstances._ivoa_date
        provenance.keywords.update(Caom2TestInstances._keywords)
        provenance.inputs.update(self.get_inputs())
        return provenance

    def get_inputs(self):
        return caom_util.TypedSet(plane.PlaneURI,
                                  plane.PlaneURI("caom:foo/bar/plane1"),
                                  plane.PlaneURI("caom:foo/bar/plane2"))

    def get_metrics(self):
        metrics = plane.Metrics()
        metrics.source_number_density = float(1.0)
        metrics.background = float(2.0)
        metrics.background_std_dev = float(3.0)
        metrics.flux_density_limit = float(4.0)
        metrics.mag_limit = float(5.0)
        if self.caom_version >= 24:
            metrics.sample_snr = 33.3
        return metrics

    def get_quality(self):
        return plane.DataQuality(plane.Quality.JUNK)

    def get_artifacts(self):
        artifacts = collections.OrderedDict()
        _artifact = artifact.Artifact("ad:foo/bar1",
                                      chunk.ProductType.SCIENCE,
                                      artifact.ReleaseType.META)
        if self.complete:
            _artifact.content_type = "application/fits"
            _artifact.content_length = int(12345)
            _artifact.last_modified = common.get_current_ivoa_time()
            if self.caom_version >= 23:
                _artifact.max_last_modified = common.get_current_ivoa_time()
                _artifact.meta_checksum = common.ChecksumURI(
                    "md5:9882dbbf9cadc221019b712fd402bcbd")
                _artifact.acc_meta_checksum = common.ChecksumURI(
                    "md5:844ce247db0844ad9f721430c80e7a21")
        if self.depth > 3:
            for k, v in self.get_parts().items():
                _artifact.parts[k] = v
        artifacts["ad:foo/bar1"] = _artifact
        if self.caom_version >= 24:
            _artifact.content_release = \
                caom_util.str2ivoa("2050-01-11T00:00:00.000")
            _artifact.content_read_groups.add("ivo://cadc.nrc.ca/gms?B")
            _artifact.content_read_groups.add("ivo://cadc.nrc.ca/gms?A")
        return artifacts

    def get_parts(self):
        parts = collections.OrderedDict()
        _part = part.Part("x")
        if self.complete:
            _part.product_type = chunk.ProductType.SCIENCE
        if self.depth > 4:
            for _chunk in self.get_chunks():
                _part.chunks.append(_chunk)
        parts["x"] = _part
        return parts

    def get_chunks(self):
        chunks = caom_util.TypedList(chunk.Chunk, )
        _chunk = chunk.Chunk()
        if self.complete:
            _chunk.product_type = chunk.ProductType.SCIENCE
            _chunk.naxis = 5
            _chunk.observable_axis = 1
            _chunk.position_axis_1 = 1
            _chunk.position_axis_2 = 2
            _chunk.energy_axis = 3
            _chunk.time_axis = 4
            _chunk.polarization_axis = 5
            _chunk.observable = self.get_observable_axis()
            _chunk.position = self.get_spatial_wcs()
            _chunk.energy = self.get_spectral_wcs()
            _chunk.time = self.get_temporal_wcs()
            _chunk.polarization = self.get_polarization_wcs()
            _chunk.last_modified = common.get_current_ivoa_time()
            if self.caom_version >= 23:
                _chunk.max_last_modified = common.get_current_ivoa_time()
                _chunk.meta_checksum = common.ChecksumURI(
                    "md5:9882dbbf9cadc221019b712fd402bcbd")
                _chunk.acc_meta_checksum = common.ChecksumURI(
                    "md5:844ce247db0844ad9f721430c80e7a21")
            if self.caom_version >= 24:
                _chunk.custom_axis = 3
                _chunk.custom = self.get_custom_wcs()
        chunks.append(_chunk)
        return chunks

    def get_observable_axis(self):
        observable = chunk.ObservableAxis(self.get_slice())
        if self.complete:
            observable.independent = self.get_slice()
        return observable

    def get_spatial_wcs(self):
        coord_axis2d = self.get_coord_axis2d()
        position = chunk.SpatialWCS(coord_axis2d)
        if self.complete:
            position.coordsys = "position coordsys"
            position.equinox = 2000.0
            position.resolution = 0.5
        return position

    def get_spectral_wcs(self):
        axis = self.get_coord_axis1d()
        energy = chunk.SpectralWCS(axis, "energy specsys")
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
            energy.transition = wcs.EnergyTransition("H", "21cm")
        return energy

    def get_temporal_wcs(self):
        axis = self.get_coord_axis1d()
        time = chunk.TemporalWCS(axis)
        if self.complete:
            time.exposure = 1.0
            time.resolution = 2.0
            time.timesys = "UTC"
            time.trefpos = "TOPOCENTER"
            time.mjdref = 3.0
        return time

    def get_polarization_wcs(self):
        axis = wcs.Axis('STOKES')
        axis1d = wcs.CoordAxis1D(axis)
        # IQUV
        axis1d.function = wcs.CoordFunction1D(int(4), 1.0,
                                              wcs.RefCoord(1.0, 1.0))
        pol = chunk.PolarizationWCS(axis1d)
        return pol

    def get_custom_wcs(self):
        axis = wcs.Axis('Foo')
        axis1d = wcs.CoordAxis1D(axis)
        axis1d.function = wcs.CoordFunction1D(int(5), 1.0,
                                              wcs.RefCoord(1.0, 1.0))
        custom = chunk.CustomWCS(axis1d)
        return custom

    def get_slice(self):
        return wcs.Slice(wcs.Axis("sliceCtype", "sliceCunit"), int(1))

    def get_coord_axis1d(self):
        coord_axis1d = wcs.CoordAxis1D(wcs.Axis("axisCtype", "axisCunit"))
        if self.complete:
            coord_axis1d.error = wcs.CoordError(1.0, 1.5)
            coord_axis1d.range = wcs.CoordRange1D(wcs.RefCoord(2.0, 2.5),
                                                  wcs.RefCoord(3.0, 3.5))
            coord_axis1d.function = (
                wcs.CoordFunction1D(4, 4.5, wcs.RefCoord(5.0, 5.5)))
            bounds = wcs.CoordBounds1D()
            bounds.samples.append(wcs.CoordRange1D(wcs.RefCoord(6.0, 6.5),
                                                   wcs.RefCoord(7.0, 7.5)))
            bounds.samples.append(wcs.CoordRange1D(wcs.RefCoord(8.0, 8.5),
                                                   wcs.RefCoord(9.0, 9.5)))
            coord_axis1d.bounds = bounds
        return coord_axis1d

    def get_coord_axis2d(self):
        axis1 = wcs.Axis("axis1Ctype", "axis1Cunit")
        axis2 = wcs.Axis("axis2Ctype", "axis2Cunit")
        coord_axis2d = wcs.CoordAxis2D(axis1, axis2)
        if self.complete:
            coord_axis2d.error1 = wcs.CoordError(1.0, 1.5)
            coord_axis2d.error2 = wcs.CoordError(2.0, 2.5)
            start = wcs.Coord2D(wcs.RefCoord(3.0, 3.5),
                                wcs.RefCoord(4.0, 4.5))
            end = wcs.Coord2D(wcs.RefCoord(5.0, 5.5),
                              wcs.RefCoord(6.0, 6.5))
            coord_axis2d.range = wcs.CoordRange2D(start, end)
            dimension = wcs.Dimension2D(7, 8)
            ref_coord = wcs.Coord2D(wcs.RefCoord(9.0, 9.5),
                                    wcs.RefCoord(10.0, 10.5))
            coord_axis2d.function = (
                wcs.CoordFunction2D(dimension, ref_coord,
                                    11.0, 12.0, 13.0, 14.0))
            if self.bounds_is_circle:
                center = wcs.ValueCoord2D(15.0, 16.0)
                coord_axis2d.bounds = wcs.CoordCircle2D(center, 17.0)
            else:
                polygon = wcs.CoordPolygon2D()
                polygon.vertices.append(wcs.ValueCoord2D(15.0, 16.0))
                polygon.vertices.append(wcs.ValueCoord2D(17.0, 18.0))
                polygon.vertices.append(wcs.ValueCoord2D(19.0, 20.0))
                coord_axis2d.bounds = polygon
        return coord_axis2d
