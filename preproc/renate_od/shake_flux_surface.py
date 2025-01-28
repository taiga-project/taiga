import matplotlib.pyplot
import pandas
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

from preproc.renate_od.interface import get_lcfs_radial_coordinate
from preproc.renate_od.manager import RenateODManager
from preproc.renate_od.beamlet import set_beamlet
from preproc.renate_od.efit import EFITManager
from preproc.renate_od.profiles import *
from preproc.renate_od.utils import *


class ProfilesMock(Profiles):
    def __init__(self, shot_number='17178', time='1097', shift=0,
                 beamlet_geometry=BeamletGeometry(),
                 efit_reconstruction_id=1, thomson_reconstruction_id=1,
                 database_directory='input/cdb', thomson_subdir='THOMSON', efit_subdir='EFITXX'):
        self.data_directory = self.get_data_directory(database_directory, shot_number)
        self.export_directory = self.get_export_directory(shot_number, time)
        thomson_directory, efit_file = self.set_path(efit_subdir, thomson_subdir, efit_reconstruction_id)
        self.thomson_profiles = ThomsonProfiles(thomson_directory, self.export_directory, shot_number, time,
                                                thomson_reconstruction_id)
        self.efit = EFITManager(efit_file, time)
        beamlet_normalised_poloidal_flux = self.efit.get_normalised_poloidal_flux(beamlet_geometry) + shift

        self.distance = beamlet_geometry.get_distance()
        self.density = self.thomson_profiles.density.get_value(beamlet_normalised_poloidal_flux) \
                       / (1 + 3 * shift)
        self.temperature = self.thomson_profiles.temperature.get_value(beamlet_normalised_poloidal_flux)


class RenateODManagerMock(RenateODManager):
    def __init__(self, beamlet_geometry, shot_number, time, species, energy, shift=0, scenario='default'):
        self.beamlet_geometry = beamlet_geometry
        self.shot_number = str(int(shot_number))
        self.time = str(int(time))
        self.species = species
        self.energy = str(int(energy))
        self.scenario = scenario

        p = ProfilesMock(beamlet_geometry=self.beamlet_geometry, shot_number=self.shot_number, time=self.time, shift=shift)
        self.distance = p.get_distance()
        self.density = p.get_density()
        self.temperature = p.get_temperature()

        self.beamlet = self.get_beamlet()
        self.relative_attenuation = self.get_relative_attenuation()


def shake_flux_surface(species, shot_number, time, energy):
    export_directory='.'
    z = 0
    tor = 0
    beamlet_geometry = set_beamlet(z, tor)

    fig, ax = matplotlib.pyplot.subplots()
    fig.set_size_inches(5, 2)
    matplotlib.pyplot.minorticks_on()
    matplotlib.pyplot.grid(which='major')
    matplotlib.pyplot.xlabel('$R$ [m]', labelpad=-10.5, loc='right')
    matplotlib.pyplot.ylabel('normalized\nneutral beam density')
    matplotlib.pyplot.title('reference discharge, ' + species + ' beam, ' + str(energy) + ' keV')
    R_LCFS = get_lcfs_radial_coordinate(shot_number=shot_number, time=time)
    matplotlib.pyplot.axvline(R_LCFS, c='red', ls='--')
    matplotlib.pyplot.text(R_LCFS+0.005, 0.56, 'reference', c='red', fontsize=6)
    matplotlib.pyplot.text(R_LCFS+0.005, 0.45, 'LCFS', c='red', fontsize=12)

    export_directory = get_home_directory() + '/input/ionProf/' + shot_number + '_' + time
    for shift in numpy.linspace(-0.02, 0.02, 5):
        r = RenateODManagerMock(beamlet_geometry, shot_number, time, species, energy, shift)
        radial_coordinate, relative_attenuation = r.get_attenuation_profile()
        label = r'with %3.f%% $\rho$' % (100*(1+shift))
        if shift == 0:
            label = r'reference  $\rho$'
        ax.plot(radial_coordinate, relative_attenuation, '-', linewidth=2, label=label,
                color=(max(abs(shift*10), shift*50), abs(shift*30), max(abs(shift*10), -shift*50)))
        relative_attenuation.fillna(0).to_csv(export_directory + '/ionyeald' + str(int(100*(1+shift)))+'.dat', index=False, header=False)
    ax.legend(labelspacing=0.3, borderpad=0)
    matplotlib.pyplot.xlim(0.6, 0.75)
    matplotlib.pyplot.subplots_adjust(left=0.18, right=0.95, top=0.85, bottom=0.15)
    matplotlib.pyplot.savefig(export_directory+'/attenuation_shaked.pdf')
    matplotlib.pyplot.savefig(export_directory+'/attenuation_shaked.svg')


def shake_flux_surface_silent(species, shot_number, time, energy):
    beamlet_geometry = BeamletGeometry()
    beamlet_geometry.rad = numpy.linspace(0.72, 0.6, 1000)
    beamlet_geometry.set_with_value(0, 'z', 'rad')
    beamlet_geometry.set_with_value(0, 'tor', 'rad')

    ref_r = RenateODManagerMock(beamlet_geometry, shot_number, time, species, energy, 0.00)
    reference_radial_coordinate, reference_relative_attenuation = ref_r.get_attenuation_profile()

    shift = 0.01
    r = RenateODManagerMock(beamlet_geometry, shot_number, time, species, energy, shift)
    radial_coordinate, relative_attenuation = r.get_attenuation_profile()
    d = relative_attenuation.diff()
    diff = numpy.abs(relative_attenuation.diff() - reference_relative_attenuation.diff()) / \
        (beamlet_geometry.rad[1]-beamlet_geometry.rad[2])
    return numpy.max(diff)


def test_multi_energy():
    a_species = 'Li'
    a_shot_number = '17178'
    a_time = '1097'
    N = 6
    energies = numpy.linspace(50, 100, N)
    max_diff = numpy.empty(N)
    for i in range(N):
        max_diff[i] = shake_flux_surface_silent(a_species, a_shot_number, a_time, energies[i])
    export_directory = get_home_directory() + '/input/ionProf/' + a_shot_number + '_' + a_time

    pandas.DataFrame(max_diff).fillna(0).to_csv(export_directory + '/maxdiff.dat',
                                          index=False, header=False)

    fig, ax = matplotlib.pyplot.subplots()
    fig.set_size_inches(5, 2)
    matplotlib.pyplot.xlabel('$E$ [keV]', labelpad=-10.5, loc='right')
    matplotlib.pyplot.ylabel('maximal difference in beam attenuation [1/m]')
    ax.plot(energies, max_diff, '-', linewidth=2)
    matplotlib.pyplot.savefig(export_directory+'/maxdiff.pdf')
    matplotlib.pyplot.savefig(export_directory+'/maxdiff.svg')


def test_single_energy():
    a_species = 'Li'
    a_shot_number = '17178'
    a_time = '1097'
    an_energy = 80
    shake_flux_surface(a_species, a_shot_number, a_time, an_energy)


def test_single_energy2():
    a_species = 'Na'
    a_shot_number = '17178'
    a_time = '1097'
    an_energy = 50
    shake_flux_surface(a_species, a_shot_number, a_time, an_energy)


if __name__ == "__main__":
    test_single_energy2()
