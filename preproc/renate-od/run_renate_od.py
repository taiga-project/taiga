import pandas
import lxml

from crm_solver.beamlet import Beamlet
from load_ts_profile import *


def get_param():
    xml_content = '<xml lang="en"><head><id>taiga beamlet</id></head><body>' \
                  '<beamlet_energy unit = "keV">80</beamlet_energy>' \
                  '<beamlet_species unit = "">Li</beamlet_species>' \
                  '<beamlet_current unit = "A">0.001</beamlet_current>' \
                  '</body></xml>'
    xml_root = lxml.etree.XML(xml_content)
    return lxml.etree.ElementTree(xml_root)


def get_components():
    components = pandas.DataFrame(
        {'q': [-1, 1], 'Z': [0, 1], 'A': [0, 2]},
        index=['electron', 'ion'])
    return components


def get_profiles(beamlet_geometry, shot_number, time):
    tuples = [('beamlet grid', 'distance', 'm'),
              ('electron', 'density', 'm-3'),
              ('electron', 'temperature', 'eV'),
              ('ion1', 'density', 'm-3'),
              ('ion1', 'temperature', 'eV')]

    header = pandas.MultiIndex.from_tuples(tuples, names=['type', 'property', 'unit'])

    #TestProfilePlot()
    #p = MockedProfiles()
    p = Profiles(beamlet_geometry=beamlet_geometry, shot_number=shot_number, time=time)
    distance = p.get_distance()
    density = p.get_density()
    temperature = p.get_temperature()

    profiles_data = numpy.transpose(numpy.array([distance, density, temperature, density, temperature]))
    profiles = pandas.DataFrame(profiles_data, columns=header)

    return profiles


def calculate_ionisation_degree(beamlet):
    pops = beamlet.profiles.filter(like='level', axis=1)
    reference_pop = pops.iat[0, 0]
    rates = pops / reference_pop
    unionisation_degree = rates.sum(1)
    ionisation_degree = 1 - unionisation_degree
    return ionisation_degree, unionisation_degree


def set_beamlet(z, tor):
    beamlet_geometry = BeamletGeometry()
    beamlet_geometry.rad = numpy.linspace(0.78, 0.35, 200)
    beamlet_geometry.set_with_value(z, 'z', 'rad')
    beamlet_geometry.set_with_value(tor, 'tor', 'rad')
    return beamlet_geometry


def calculate_beamlet(beamlet_geometry, shot_number, time):
    beamlet = Beamlet(param=get_param(),
                      profiles=get_profiles(beamlet_geometry=beamlet_geometry, shot_number=shot_number, time=time),
                      components=get_components())
    ionisation_degree, attenuation_degree = calculate_ionisation_degree(beamlet)
    radial_coordinate = pandas.DataFrame(beamlet_geometry.rad)
    return radial_coordinate, ionisation_degree, attenuation_degree


def export_beamlet_profile(export_directory='data/output/matyi'):
    shot_number = '17178'
    time = '1097'
    z = 0
    tor = 0
    beamlet_geometry = set_beamlet(z, tor)
    radial_coordinate, ionisation_degree, attenuation_degree = calculate_beamlet(beamlet_geometry, shot_number, time)
    radial_coordinate.to_csv(export_directory+'/rad.dat', index=False, header=False)
    attenuation_degree.to_csv(export_directory+'/ionyeald.dat', index=False, header=False)
    plot_ionisation_profile(radial_coordinate, ionisation_degree)


def plot_ionisation_profile(radial_coordinate, ionisation_degree):
    fig, ax = matplotlib.pyplot.subplots()
    ax.plot(radial_coordinate, ionisation_degree, '.')
    matplotlib.pyplot.show()
    print(ionisation_degree)


def mock_beam(diameter=5e-3, z_length=3, tor_length=3):
    shot_number = '17178'
    time = '1097'
    fig, ax = matplotlib.pyplot.subplots()
    for z in numpy.linspace(-diameter/2, diameter/2, z_length):
        if z_length == 1:
            z = 0
        for tor in numpy.linspace(-diameter/2, diameter/2, tor_length):
            if tor_length == 1:
                tor = 0
            beamlet_geometry = set_beamlet(z, tor)
            radial_coordinate, ionisation_degree,  attenuation_degree = \
                calculate_beamlet(beamlet_geometry, shot_number, time)
            ax.plot(radial_coordinate, attenuation_degree, '-')

    ax.set_xlabel('R [m]')
    ax.set_ylabel('normalised linear density attenuation')
    ax.set_title('COMPASS #'+shot_number+' ('+time+' ms) ')
    matplotlib.pyplot.xlim(0.6, 0.75)
    matplotlib.pyplot.savefig('attenuation_'+shot_number+'_'+time+'.svg')
    matplotlib.pyplot.show()


if __name__ == "__main__":
     export_beamlet_profile()
    #mock_beam()
