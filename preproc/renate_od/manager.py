import pandas
import lxml.etree
import sys

from profiles import *

sys.path.append('../../ext/renate_od')
from crm_solver.beamlet import Beamlet


class RenateODManager:
    def __init__(self, beamlet_geometry, shot_number, time, species, energy):
        self.beamlet_geometry = beamlet_geometry
        self.shot_number = shot_number
        self.time = time
        self.species = species
        self.energy = energy
        self.beamlet = self.get_beamlet()
        self.relative_attenuation = self.get_relative_attenuation()

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

        p = Profiles(beamlet_geometry=self.beamlet_geometry, shot_number=self.shot_number, time=self.time)
        distance = p.get_distance()
        density = p.get_density()
        temperature = p.get_temperature()

        profiles_data = numpy.transpose(numpy.array([distance, density, temperature, density, temperature]))
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
