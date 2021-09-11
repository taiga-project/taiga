import numpy


def set_beamlet(z, tor):
    beamlet_geometry = BeamletGeometry()
    beamlet_geometry.rad = numpy.linspace(0.78, 0.35, 200)
    beamlet_geometry.set_with_value(z, 'z', 'rad')
    beamlet_geometry.set_with_value(tor, 'tor', 'rad')
    return beamlet_geometry


class BeamletGeometry:
    def __init__(self):
        self.rad = []
        self.z = []
        self.tor = []
        self.set_default_values()

    def set_with_value(self, value, field, reference):
        value_array = numpy.full_like(getattr(self, reference), value)
        setattr(self, field, value_array)

    def set_default_values(self):
        self.rad = numpy.linspace(0.78, 0.35, 200)
        self.set_with_value(0, 'z', 'rad')
        self.set_with_value(0, 'tor', 'rad')

    def get_distance(self):
        return self.rad[0]-self.rad
