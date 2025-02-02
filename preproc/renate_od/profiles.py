from preproc.renate_od.beamlet import BeamletGeometry
from preproc.renate_od.thomson import ThomsonProfiles
from preproc.renate_od.efit import EFITManager
from preproc.renate_od.utils import *


class Profiles:
    def __init__(self, shot_number='17178', time='1097',
                 beamlet_geometry=BeamletGeometry(),
                 efit_reconstruction_id=1, thomson_reconstruction_id=1,
                 database_directory='input/cdb', thomson_subdir='THOMSON', efit_subdir='EFITXX',
                 is_export=True, is_plot=True):
        self.data_directory = self.get_data_directory(database_directory, shot_number)
        self.export_directory = self.get_export_directory(shot_number, time)
        thomson_directory, efit_file = self.set_path(efit_subdir, thomson_subdir, efit_reconstruction_id)
        self.thomson_profiles = ThomsonProfiles(thomson_directory, self.export_directory, shot_number, time, thomson_reconstruction_id)
        self.efit = EFITManager(efit_file, time)
        beamlet_normalised_poloidal_flux = self.efit.get_normalised_poloidal_flux(beamlet_geometry)

        self.distance = beamlet_geometry.get_distance()
        self.density = self.thomson_profiles.density.get_value(beamlet_normalised_poloidal_flux)
        self.temperature = self.thomson_profiles.temperature.get_value(beamlet_normalised_poloidal_flux)
        if is_export:
            self.thomson_profiles.export_profiles()
        if is_plot:
            self.thomson_profiles.plot_profiles()

    @staticmethod
    def get_data_directory(database_directory, shot_number):
        return get_home_directory() + '/' + database_directory + '/' + str(shot_number)

    @staticmethod
    def get_export_directory(shot_number, time):
        return get_home_directory() + '/' + '/input/tsProf' '/' + str(shot_number) + '_' + str(time)

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


if __name__ == "__main__":
    Profiles(is_export=False, is_plot=True)
