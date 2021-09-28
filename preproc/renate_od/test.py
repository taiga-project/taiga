import numpy
import matplotlib.pyplot
import unittest

from manager import RenateODManager
from profiles import Profiles
from beamlet import BeamletGeometry, set_beamlet
from thomson import ThomsonProfiles
from efit import EFITManager
from utils import *


class TestROD(unittest.TestCase):
    @staticmethod
    def test_normalised_poloidal_flux_profile(shot_number='17178', time='1097'):
        thomson_directory = get_home_directory() + '/input/cdb/' + str(shot_number) + '/THOMSON'
        efit_file = get_home_directory() + '/input/cdb/' + str(shot_number) + '/EFITXX/EFITXX.1.h5'

        ts = ThomsonProfiles(thomson_directory, shot_number, time)
        ts_flux = ts.normalised_poloidal_flux_profile

        beamlet = BeamletGeometry()
        beamlet.z = ts.get_dataset('TS_z_axis')
        beamlet.set_with_value(0, 'tor', 'z')

        efit = EFITManager(efit_file, time)

        R1 = 0.5565
        R2 = 0.5575
        R_TS = 0.557

        beamlet.set_with_value(R1, 'rad', 'z')
        efit_flux1 = efit.get_normalised_poloidal_flux(beamlet)
        beamlet.set_with_value(R2, 'rad', 'z')
        efit_flux2 = efit.get_normalised_poloidal_flux(beamlet)

        R = (R1 * (efit_flux2 - ts_flux) + R2 * (ts_flux - efit_flux1)) / (efit_flux2 - efit_flux1)

        numpy.testing.assert_allclose(R, R_TS, atol=1e-4)

    @staticmethod
    def test_lcfs(shot_number='17178', time='1097'):
        efit_file = get_home_directory() + '/input/cdb/' + str(shot_number) + '/EFITXX/EFITXX.1.h5'
        efit = EFITManager(efit_file, time)
        beamlet = BeamletGeometry()
        R_reference = efit.get_time_sliced_data('output/separatrixGeometry/rmidplaneOut')
        beamlet.rad = numpy.linspace(0.6, 0.75, 100)
        beamlet.set_with_value(0, 'z', 'rad')
        beamlet.set_with_value(0, 'tor', 'rad')

        flux = efit.get_normalised_poloidal_flux(beamlet)
        R_LCFS = numpy.interp(1.0, flux, beamlet.rad)

        numpy.testing.assert_allclose(R_LCFS, R_reference, atol=1e-5)


def mock_beam(shot_number='17178', time='1097', diameter=5e-3, z_length=3, tor_length=3):
    fig, ax = matplotlib.pyplot.subplots()
    species = 'Li'
    energy = '80'
    for z in numpy.linspace(-diameter/2, diameter/2, z_length):
        if z_length == 1:
            z = 0
        for tor in numpy.linspace(-diameter/2, diameter/2, tor_length):
            if tor_length == 1:
                tor = 0
            beamlet_geometry = set_beamlet(z, tor)
            r = RenateODManager(beamlet_geometry, shot_number, time, species, energy)
            radial_coordinate, relative_attenuation = r.get_attenuation_profile()
            ax.plot(radial_coordinate, relative_attenuation, '-')

    ax.set_xlabel('R [m]')
    ax.set_ylabel('normalised linear density attenuation')
    ax.set_title('COMPASS #'+shot_number+' ('+time+' ms) ')
    matplotlib.pyplot.xlim(0.6, 0.75)
    matplotlib.pyplot.savefig('attenuation_'+shot_number+'_'+time+'.svg')
    matplotlib.pyplot.show()


# noinspection PyMissingConstructor
class MockedProfiles(Profiles):
    def __init__(self):
        self.distance = [0, 0.1, 0.2]
        self.density = [1e19, 1.2e19, 1.5e19]
        self.temperature = [1000, 1500, 2500]


class TestProfilePlot(Profiles):
    reference_beamlet = BeamletGeometry()

    def __init__(self, shot_number='17178', time='1097',
                 beamlet_geometry=reference_beamlet, ts_time_source='EFIT',
                 efit_reconstruction_id=1, thomson_reconstruction_id=1, database_directory='input/cdb',
                 thomson_subdir='THOMSON', efit_subdir='EFITXX'):
        super().__init__(shot_number, time, beamlet_geometry, ts_time_source, efit_reconstruction_id,
                         thomson_reconstruction_id, database_directory, thomson_subdir, efit_subdir)
        self.plot()

    def plot(self):
        ts_normalised_poloidal_flux = self.thomson_profiles.normalised_poloidal_flux_profile

        ts_beamlet = BeamletGeometry()
        ts_beamlet.z = self.thomson_profiles.z_axis
        ts_beamlet.rad = numpy.linspace(0.556, 0.558, 11)
        ts_beamlet.set_with_value(0, 'tor', 'rad')

        efit_normalised_poloidal_flux = self.efit.get_normalised_poloidal_flux(ts_beamlet, grid=True)

        fig, ax = matplotlib.pyplot.subplots()
        ax.plot(efit_normalised_poloidal_flux.T, (ts_normalised_poloidal_flux - efit_normalised_poloidal_flux).T, '.')

        matplotlib.pyplot.legend(ts_beamlet.rad.round(4), bbox_to_anchor=(0.92, 1), loc='upper left')
        matplotlib.pyplot.show()


if __name__ == "__main__":
    a_shot_number = '17178'
    a_time = '1097'
    #mock_beam(shot_number=a_shot_number, time=a_time)
    unittest.main()
