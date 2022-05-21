import unittest

from ..beamlet import BeamletGeometry
from ..thomson import ThomsonProfiles
from ..efit import EFITManager
from ..utils import *


class TestROD(unittest.TestCase):
    @staticmethod
    def test_normalised_poloidal_flux_profile(shot_number='17178', time='1097'):
        thomson_directory = get_home_directory() + '/input/cdb/' + str(shot_number) + '/THOMSON'
        efit_file = get_home_directory() + '/input/cdb/' + str(shot_number) + '/EFITXX/EFITXX.1.h5'

        ts = ThomsonProfiles(thomson_directory, '/tmp', shot_number, time)
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


if __name__ == "__main__":
    unittest.main()
