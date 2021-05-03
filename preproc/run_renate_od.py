from crm_solver.beamlet import Beamlet
import pandas
import numpy
import lxml
import os


def get_param():
    xml_content = '<xml lang="en"><head><id>taiga beamlet</id></head><body>' \
                  '<beamlet_energy unit = "keV">60</beamlet_energy>' \
                  '<beamlet_species unit = "">Na</beamlet_species>' \
                  '<beamlet_current unit = "A">0.001</beamlet_current>' \
                  '</body></xml>'
    xml_root = lxml.etree.XML(xml_content)
    return lxml.etree.ElementTree(xml_root)


def get_components():
    components = pandas.DataFrame(
        {'q': [-1, 1], 'Z': [0, 1], 'A': [0, 2]},
        index=['electron', 'ion'])
    return components


def load_profiles():
    distance = [0, 0.1, 0.2]
    density = [1e19, 1.2e19, 1.5e19]
    temperature = [1000, 1500, 2500]
    return distance, density, temperature


def get_profiles():
    tuples = [('beamlet grid', 'distance', 'm'),
              ('electron', 'density', 'm-3'),
              ('electron', 'temperature', 'eV'),
              ('ion1', 'density', 'm-3'),
              ('ion1', 'temperature', 'eV')]

    header = pandas.MultiIndex.from_tuples(tuples, names=['type', 'property', 'unit'])
    distance, density, temperature = load_profiles()

    profiles_data = numpy.transpose(numpy.array([distance, density, temperature, density, temperature]))
    profiles = pandas.DataFrame(profiles_data, columns=header)

    return profiles


def export_beamlet_profile():
    b = Beamlet(param=get_param(), profiles=get_profiles(), components=get_components())
    pops = b.profiles.filter(like='level', axis=1)
    reference_pop = pops.iat[0, 0]
    rates = pops / reference_pop
    ionisation_degree = 1 - rates.sum(1)
    print(ionisation_degree)

    b.profiles['beamlet grid'].to_csv('data/output/matyi/rad.txt', index=False, header=False)
    ionisation_degree.to_csv('data/output/matyi/degree.txt', index=False, header=False)

    print(b.profiles)


if __name__ == "__main__":
    export_beamlet_profile()
