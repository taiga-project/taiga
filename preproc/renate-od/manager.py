import pandas
import lxml

from crm_solver.beamlet import Beamlet
from load_ts_profile import *


class RenateODManager:
    def __init__(self, beamlet_geometry, shot_number, time, species, energy):
        self.beamlet_geometry = beamlet_geometry
        self.shot_number = shot_number
        self.time = time
        self.species = species
        self.energy = energy

    @staticmethod
    def set_param(self):
        xml_content = '<xml lang="en"><head><id>taiga beamlet</id></head><body>' \
                      '<beamlet_energy unit = "keV">' + self.energy + '</beamlet_energy>' \
                      '<beamlet_species unit = "">' + self.species + '</beamlet_species>' \
                      '<beamlet_current unit = "A">0.001</beamlet_current>' \
                      '</body></xml>'
        xml_root = lxml.etree.XML(xml_content)
        self.param = lxml.etree.ElementTree(xml_root)

    @staticmethod
    def set_components(self):
        self.components = pandas.DataFrame(
            {'q': [-1, 1], 'Z': [0, 1], 'A': [0, 2]},
            index=['electron', 'ion'])

    @staticmethod
    def set_profiles(self):
        tuples = [('beamlet grid', 'distance', 'm'),
                  ('electron', 'density', 'm-3'),
                  ('electron', 'temperature', 'eV'),
                  ('ion1', 'density', 'm-3'),
                  ('ion1', 'temperature', 'eV')]

        header = pandas.MultiIndex.from_tuples(tuples, names=['type', 'property', 'unit'])

        p = Profiles(beamlet_geometry=self.beamlet_geometry, shot_number=self.shot_number, time=self.time)
        distance = p.get_distance()
        density = p.get_density()
        temperature = p.get_temperature()

        profiles_data = numpy.transpose(numpy.array([distance, density, temperature, density, temperature]))
        self.profiles = pandas.DataFrame(profiles_data, columns=header)

    @staticmethod
    def calculate_relative_attenuation(self):
        if not hasattr(self.__class__, 'relative_attenuation'):
            self.beamlet.compute_linear_density_attenuation()
            absolute_attenuation = self.beamlet.profiles['linear_density_attenuation']
            self.relative_attenuation = absolute_attenuation / absolute_attenuation.max()

    @staticmethod
    def calculate_beamlet(self):
        self.set_param()
        self.set_profiles()
        self.set_components()
        self.beamlet = Beamlet(param=self.param, profiles=self.profiles, components=self.components)
        self.calculate_relative_attenuation()

    @staticmethod
    def get_attenuation_profile(self):
        radial_coordinate = pandas.DataFrame(self.beamlet_geometry.rad)
        if not hasattr(self.__class__, 'beamlet'):
            self.calculate_beamlet()
        if not hasattr(self.__class__, 'relative_attenuation'):
            self.calculate_relative_attenuation()
        return radial_coordinate, self.relative_attenuation
