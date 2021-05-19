import pandas
import lxml

from crm_solver.beamlet import Beamlet
from load_ts_profile import *


def get_param():
    xml_content = '<xml lang="en"><head><id>taiga beamlet</id></head><body>' \
                  '<beamlet_energy unit = "keV">40</beamlet_energy>' \
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


def get_profiles(beamlet_geometry):
    tuples = [('beamlet grid', 'distance', 'm'),
              ('electron', 'density', 'm-3'),
              ('electron', 'temperature', 'eV'),
              ('ion1', 'density', 'm-3'),
              ('ion1', 'temperature', 'eV')]

    header = pandas.MultiIndex.from_tuples(tuples, names=['type', 'property', 'unit'])

    TestProfilePlot()
    p = MockedProfiles()
    p = Profiles(beamlet_geometry=beamlet_geometry)
    distance = p.get_distance()
    density = p.get_density()
    temperature = p.get_temperature()

    profiles_data = numpy.transpose(numpy.array([distance, density, temperature, density, temperature]))
    profiles = pandas.DataFrame(profiles_data, columns=header)

    return profiles


def export_beamlet_profile(export_directory='data/output/matyi'):

    beamlet_geometry = BeamletGeometry()
    beamlet_geometry.rad = numpy.linspace(0.78, 0.35, 200)
    beamlet_geometry.set_with_value(0, 'z', 'rad')
    beamlet_geometry.set_with_value(0, 'tor', 'rad')

    beamlet = Beamlet(param=get_param(), profiles=get_profiles(beamlet_geometry=beamlet_geometry),
                      components=get_components())
    pops = beamlet.profiles.filter(like='level', axis=1)
    reference_pop = pops.iat[0, 0]
    rates = pops / reference_pop
    ionisation_degree = 1 - rates.sum(1)
    radial_coordinate = pandas.DataFrame(beamlet_geometry.rad)

    radial_coordinate.to_csv(export_directory+'/rad.dat', index=False, header=False)
    rates.sum(1).to_csv(export_directory+'/ionyeald.dat', index=False, header=False)

    plot_ionisation_profile(radial_coordinate, ionisation_degree)


def plot_ionisation_profile(radial_coordinate, ionisation_degree):
    fig, ax = matplotlib.pyplot.subplots()
    ax.plot(radial_coordinate, ionisation_degree, '.')
    matplotlib.pyplot.show()
    print(ionisation_degree)


if __name__ == "__main__":
    export_beamlet_profile()
