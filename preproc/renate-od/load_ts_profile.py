from statsmodels.nonparametric.smoothers_lowess import lowess

from beamlet import BeamletGeometry
from thomson import ThomsonProfiles
from efit import EFITManager
from utils import *


class Profiles:
    def __init__(self, shot_number='17178', time='1097',
                 beamlet_geometry=BeamletGeometry(),
                 ts_time_source='EFIT',
                 efit_reconstruction_id=1, thomson_reconstruction_id=1,
                 database_directory='input/cdb', thomson_subdir='THOMSON', efit_subdir='EFITXX'):
        self.set_data_directory(database_directory, shot_number)
        thomson_directory, efit_file = self.set_path(efit_subdir, thomson_subdir, efit_reconstruction_id)
        self.thomson_profiles = ThomsonProfiles(thomson_directory, shot_number, time, thomson_reconstruction_id,
                                                ts_time_source)
        self.efit = EFITManager(efit_file, time)
        beamlet_normalised_poloidal_flux = self.efit.get_normalised_poloidal_flux(beamlet_geometry)

        self.distance = beamlet_geometry.get_distance()
        self.density = self.thomson_profiles.density.get_value(beamlet_normalised_poloidal_flux)
        self.temperature = self.thomson_profiles.temperature.get_value(beamlet_normalised_poloidal_flux)
        self.export_profiles()
        self.thomson_profiles.plot_profiles()

    @staticmethod
    def set_data_directory(self, database_directory, shot_number):
        self.data_directory = get_home_directory() + '/' + database_directory + '/' + str(shot_number)

    def set_path(self, efit_subdir, thomson_subdir, efit_reconstruction_id):
        thomson_directory = self.data_directory + '/' + thomson_subdir
        efit_file = self.data_directory + '/' + efit_subdir + '/' + 'EFITXX.' + str(efit_reconstruction_id) + '.h5'
        return thomson_directory, efit_file

    def get_distance(self):
        return self.distance

    def get_density(self):
        return self.density

    def get_temperature(self):
        return self.temperature

    def export_profiles(self, path=get_home_directory() + '/input/tsProf'):
        self.thomson_profiles.export_profiles(path)
