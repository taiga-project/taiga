import pandas
import lxml.etree
import sys

from profiles import *

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             os.pardir, os.pardir, 'ext', 'renate_od')))
from crm_solver.beamlet import Beamlet


class RenateODManager:
    def __init__(self, beamlet_geometry, shot_number, time, species, energy, scenario='default'):
        self.beamlet_geometry = beamlet_geometry
        self.shot_number = str(shot_number)
        self.time = str(time)
        self.species = species
        self.energy = str(energy)
        self.scenario = scenario

        p = Profiles(beamlet_geometry=self.beamlet_geometry, shot_number=self.shot_number, time=self.time)
        self.distance = p.get_distance()
        self.density = p.get_density()
        self.temperature = p.get_temperature()

        self.beamlet = self.get_beamlet()
        self.relative_attenuation = self.get_relative_attenuation()

    def get_electron_density(self):
        if self.scenario in ['default', 'just electron']:
            return self.density
        else:
            return numpy.zeros_like(self.density)

    def get_electron_temperature(self):
        if self.scenario in ['default', 'just electron']:
            return self.temperature
        else:
            return numpy.zeros_like(self.temperature)

    def get_ion_density(self):
        if self.scenario in ['default', 'just ion']:
            return self.density
        else:
            return numpy.zeros_like(self.density)

    def get_ion_temperature(self):
        if self.scenario in ['default', 'just ion']:
            return self.temperature
        else:
            return numpy.zeros_like(self.temperature)

    def get_param(self):
        xml_content = '<xml lang="en"><head><id>taiga beamlet</id></head><body>' \
                      '<beamlet_energy unit = "keV">' + self.energy + '</beamlet_energy>' \
                                                                      '<beamlet_species unit = "">' + self.species + '</beamlet_species>' \
                                                                                                                     '<beamlet_current unit = "A">0.001</beamlet_current>' \
                                                                                                                     '</body></xml>'
        xml_root = lxml.etree.XML(xml_content)
        return lxml.etree.ElementTree(xml_root)

    @staticmethod
    def get_components():
        return pandas.DataFrame(
            {'q': [-1, 1], 'Z': [0, 1], 'A': [0, 2]},
            index=['electron', 'ion'])

    def get_profiles(self):
        tuples = [('beamlet grid', 'distance', 'm'),
                  ('electron', 'density', 'm-3'),
                  ('electron', 'temperature', 'eV'),
                  ('ion1', 'density', 'm-3'),
                  ('ion1', 'temperature', 'eV')]

        header = pandas.MultiIndex.from_tuples(tuples, names=['type', 'property', 'unit'])
        profiles_data = numpy.transpose(numpy.array([self.distance,
                                                     self.get_electron_density(), self.get_electron_temperature(),
                                                     self.get_ion_density(), self.get_ion_temperature()]))
        return pandas.DataFrame(profiles_data, columns=header)

    def get_relative_attenuation(self):
        self.beamlet.compute_linear_density_attenuation()
        absolute_attenuation = self.beamlet.profiles['linear_density_attenuation']
        return absolute_attenuation / absolute_attenuation.max()

    def get_beamlet(self):
        param = self.get_param()
        profiles = self.get_profiles()
        components = self.get_components()
        return Beamlet(param=param, profiles=profiles, components=components)

    def get_attenuation_profile(self):
        radial_coordinate = pandas.DataFrame(self.beamlet_geometry.rad)
        return radial_coordinate, self.relative_attenuation
