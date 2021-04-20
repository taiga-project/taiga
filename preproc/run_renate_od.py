from crm_solver.beamlet import Beamlet
from utility.writedata import WriteData
from visualization.profiles import BeamletProfiles
import pandas
import numpy

# 54b5f90

components = pandas.DataFrame(
    {'q': [-1, 1], 'Z': [0, 1], 'A': [0, 2]},
    index=['electron', 'ion'])


tuples = [('beamlet grid', 'distance', 'm'),
          ('electron', 'density', 'm-3'),
          ('electron', 'temperature', 'eV'),
          ('ion1', 'density', 'm-3'),
          ('ion1', 'temperature', 'eV')]

header = pandas.MultiIndex.from_tuples(tuples, names=['type', 'property', 'unit'])

distance = [0, 0.1, 0.2]
density = [1e19, 1.2e19, 1.5e19]
temperature = [1000, 1500, 2500]
profiles_data = numpy.transpose(numpy.array([distance, density, temperature, density, temperature]))
profiles = pandas.DataFrame(profiles_data, columns=header)

# run
b = Beamlet(profiles=profiles, components=components, data_path="beamlet/taiga_beamlet.xml")
pops = b.profiles.filter(like='level', axis=1)
reference_pop = pops.iat[0, 0]
rates = pops/reference_pop
ionisation_degree = 1-rates.sum(1)
print(ionisation_degree)

#b = Beamlet(data_path="beamlet/H_plasma_Na_beam@60keV.xml")
if False:
    ionisation_degree = 1
    for i in range(7):
        ionisation_degree -= b.profiles['RENATE level '+str(i)]
    b.profiles['Degree of ionisation'] = ionisation_degree
    print(b.profiles)
    print(b.profiles['beamlet grid'].to_numpy)
    print(ionisation_degree.to_numpy)
    print(b.profiles['Degree of ionisation'])
b.profiles.to_html('data/output/matyi/temp.html')
b.profiles.to_csv('data/output/matyi/temp.txt')
b.profiles['beamlet grid'].to_csv('data/output/matyi/rad.txt', index=False, header=False)
ionisation_degree.to_csv('data/output/matyi/degree.txt', index=False, header=False)

#b.compute_relative_populations(reference_level='level 0')
if False:
    w = WriteData()
    w.write_beamlet_profiles(b.param, b.profiles)
    p = BeamletProfiles('output/beamlet/H_plasma_Na_beam@60keV.xml')
    p.plot_all_profiles()
