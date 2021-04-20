from crm_solver.beamlet import Beamlet
import pandas
import numpy


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

b = Beamlet(profiles=profiles, components=components, data_path="beamlet/taiga_beamlet.xml")
pops = b.profiles.filter(like='level', axis=1)
reference_pop = pops.iat[0, 0]
rates = pops/reference_pop
ionisation_degree = 1-rates.sum(1)
print(ionisation_degree)

b.profiles.to_html('data/output/matyi/temp.html')
b.profiles.to_csv('data/output/matyi/temp.txt')
b.profiles['beamlet grid'].to_csv('data/output/matyi/rad.txt', index=False, header=False)
ionisation_degree.to_csv('data/output/matyi/degree.txt', index=False, header=False)
