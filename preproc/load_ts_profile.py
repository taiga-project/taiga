import os
import h5py
import matplotlib
import numpy
import scipy
from statsmodels.nonparametric.smoothers_lowess import lowess

from stefanikova2016 import *


class Profiles:
    def __init__(self, shot_number='17178', time='1097', ts_time_source='EFIT',
                 efit_reconstruction_id=1, thomson_reconstruction_id=2,
                 r_min=0.35, r_max=0.78, r_points=200,
                 database_directory='input/cdb', thomson_subdir='THOMSON', efit_subdir='EFITXX'):
        self.data_directory = None
        self.set_data_directory(database_directory, shot_number)
        thomson_directory, efit_file = self.set_path(efit_subdir, thomson_subdir, efit_reconstruction_id)
        thomson_profiles = ThomsonProfiles(thomson_directory, time, thomson_reconstruction_id, ts_time_source)
        efit = EFITManager(efit_file, time)
        beamlet_r = numpy.linspace(r_max, r_min, r_points)
        beamlet_z = numpy.full_like(beamlet_r, 0)
        beamlet_tor = numpy.full_like(beamlet_r, 0)
        beamlet_normalised_poloidal_flux = efit.get_normalised_poloidal_flux(beamlet_r, beamlet_z, beamlet_tor)

        self.distance = beamlet_r[0] - beamlet_r
        self.density = thomson_profiles.density.get_value(beamlet_normalised_poloidal_flux)
        self.temperature = thomson_profiles.temperature.get_value(beamlet_normalised_poloidal_flux)

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


class MockedProfiles(Profiles):
    def __init__(self):
        self.distance = [0, 0.1, 0.2]
        self.density = [1e19, 1.2e19, 1.5e19]
        self.temperature = [1000, 1500, 2500]


def get_home_directory():
    try:
        return os.environ['TAIGA_HOME']
    except KeyError:
        return '/home/matyi/work/taiga_local'  # '.'


class TestProfilePlot(Profiles):
    def __init__(self, shot_number='17178', time='1097', ts_time_source='EFIT',
                 efit_reconstruction_id=1, thomson_reconstruction_id=2,
                 database_directory='input/cdb', thomson_subdir='THOMSON', efit_subdir='EFITXX'):
        self.data_directory = None
        self.set_data_directory(database_directory, shot_number)
        thomson_directory, efit_file = self.set_path(efit_subdir, thomson_subdir, efit_reconstruction_id)
        thomson_profiles = ThomsonProfiles(thomson_directory, time, thomson_reconstruction_id, ts_time_source)
        thomson_profiles.plot_profiles()

        ts_normalised_poloidal_flux = thomson_profiles.normalised_poloidal_flux_profile
        ts_r = numpy.linspace(0.556, 0.558, 11)
        ts_z = thomson_profiles.z_axis
        efit = EFITManager(efit_file, time)
        efit_normalised_poloidal_flux = efit.get_normalised_poloidal_flux(ts_r, ts_z, grid=True)

        fig, ax = matplotlib.pyplot.subplots()
        ax.plot(efit_normalised_poloidal_flux.T, (ts_normalised_poloidal_flux - efit_normalised_poloidal_flux).T, '.')

        matplotlib.pyplot.legend(ts_r.round(4), bbox_to_anchor=(0.92, 1), loc='upper left')
        matplotlib.pyplot.show()


class EFITManager:
    def __init__(self, efit_file, time):
        self.efit_file = efit_file
        self.time = time
        self.time_index = self.get_time_index()
        print(self.time_index)
        print(time)
        self.r = self.get_time_sliced_data('output/profiles2D/r')
        self.z = self.get_time_sliced_data('output/profiles2D/z')
        self.poloidal_flux_profile = self.get_time_sliced_data('output/profiles2D/poloidalFlux')
        self.normalise_poloidal_flux = None
        self.get_poloidal_flux = None
        self.set_normalised_poloidal_flux()
        self.set_poloidal_flux()

    def get_time_index(self):
        time_dataset = self.get_data('time')
        return (numpy.abs(time_dataset - int(self.time) / 1000)).argmin()

    def get_data(self, field):
        file = h5py.File(self.efit_file)
        return file[field][()]

    def get_time_sliced_data(self, field):
        file = h5py.File(self.efit_file)
        return file[field][self.time_index]

    def set_normalised_poloidal_flux(self):
        normalised_poloidal_flux = self.get_data('output/fluxFunctionProfiles/normalizedPoloidalFlux')
        poloidal_flux = self.get_time_sliced_data('output/fluxFunctionProfiles/poloidalFlux')
        self.normalise_poloidal_flux = \
            scipy.interpolate.UnivariateSpline(poloidal_flux, normalised_poloidal_flux)

    def set_poloidal_flux(self):
        print(self.r.shape)
        print(self.z.shape)
        print(self.poloidal_flux_profile.shape)
        self.get_poloidal_flux = scipy.interpolate.RectBivariateSpline(self.r, self.z, self.poloidal_flux_profile)

    def get_normalised_poloidal_flux(self, r, z, tor=numpy.array(0), grid=False):
        r_poloidal_cross_section = self.transform_to_poloidal_cross_section(r, tor)
        poloidal_flux = self.get_poloidal_flux(r_poloidal_cross_section, z, grid=grid)
        return self.normalise_poloidal_flux(poloidal_flux)

    def transform_to_poloidal_cross_section(self, r, tor):
        return numpy.hypot(r, tor)


class ThomsonProfiles:
    def __init__(self, thomson_directory, time, reconstruction_id, ts_time_source):
        self.thomson_directory = thomson_directory
        self.reconstruction_id = reconstruction_id
        self.ts_time_source = ts_time_source
        self.time = time

        self.time_index = []
        self.z_axis = []
        self.temperature_profile = []
        self.temperature_error_profile = []
        self.density_profile = []
        self.density_error_profile = []
        self.normalised_poloidal_flux_profile = []

        self.read_thomson_database()
        self.density = ProfileManager(self.normalised_poloidal_flux_profile,
                                      self.density_profile, self.density_error_profile)
        self.temperature = ProfileManager(self.normalised_poloidal_flux_profile,
                                          self.temperature_profile, self.temperature_error_profile)

    def read_thomson_database(self):
        try:
            time_dataset = self.get_dataset('TS_' + self.ts_time_source + '_time', reconstruction_id=1)
            self.z_axis = self.get_dataset('TS_z_axis', self.reconstruction_id)
            print('Thomson scattering time and geometry files read successfully from: ' + self.thomson_directory)
        except OSError:
            raise OSError('Invalid Thomson scattering file structure! '
                          '\nPlease check the directory tree here:\n' + self.thomson_directory)
        except KeyError:
            raise KeyError(
                'Invalid Thomson scattering data structure!\nExample for a correct structure:\nTe.1.h5\n\tTe')

        self.time_index = self.get_time_index(time_dataset)

        try:
            self.temperature_profile = self.get_profile('Te', reconstruction_id=self.reconstruction_id)
            self.temperature_error_profile = self.get_profile('Te_err', reconstruction_id=self.reconstruction_id)
            self.density_profile = self.get_profile('ne', reconstruction_id=self.reconstruction_id)
            self.density_error_profile = self.get_profile('ne_err', reconstruction_id=self.reconstruction_id)
            self.normalised_poloidal_flux_profile = self.get_profile('psi_n', reconstruction_id=1)
            print('Thomson scattering profile files read successfully from: ' + self.thomson_directory)
        except OSError:
            raise OSError('Invalid Thomson scattering file structure! '
                          '\nPlease check the directory tree here:\n' + self.thomson_directory)
        except KeyError:
            raise KeyError(
                'Invalid Thomson scattering data structure!\nExample for a correct structure:\nTe.1.h5\n\tTe')

    def plot_profiles(self):
        self.density.plot_profile()
        self.temperature.plot_profile()

    def get_time_index(self, time_dataset):
        return (numpy.abs(time_dataset - int(self.time))).argmin()

    def get_dataset(self, field, reconstruction_id):
        file = h5py.File(self.thomson_directory + '/' + field + '.' + str(reconstruction_id) + '.h5')
        return file[field][()]

    def get_profile(self, field, reconstruction_id):
        return self.get_dataset(field, reconstruction_id)[self.time_index]


def set_negatives_to_zero(values):
    return values * (values > 0)


class ProfileManager:
    def __init__(self, x, y, y_error):
        self.x_raw = x
        self.y_raw = y
        self.y_error_raw = y_error
        self.x = x
        self.y = y
        self.y_error = y_error
        self.y_smooth = []
        self.f_smooth = ValueError

        self.x_fine = []
        self.y_fine = []

        self.get_valid_profile()
        self.get_smoothed_profile()

    def get_valid_profile(self):
        valid_indices = self.get_valid_profile_indices()
        self.x = self.x_raw[valid_indices]
        self.y = self.y_raw[valid_indices]
        self.y_error = self.y_error_raw[valid_indices]

    def get_valid_profile_indices(self):
        return numpy.where(~numpy.isnan(self.y) & (self.x > 0.0) & (self.x < 1.2) &
                           (self.y_error < self.y) & self.filter_sol_outliers())

    def filter_sol_outliers(self):
        return (numpy.fmin.accumulate(self.y) == self.y) | (self.x < 1.02)

    def get_smoothed_profile(self):
        self.smooth_profile()
        self.refine_smoothed_profile()

    def smooth_profile(self):
        smooth = lowess(self.y, self.x, is_sorted=True, frac=0.3, it=0)
        y_smooth = smooth[:, 1]
        self.y_smooth = set_negatives_to_zero(y_smooth)

    def refine_smoothed_profile(self):
        self.interpolate_smoothed_profile()
        self.x_fine = numpy.linspace(self.x[0], 1.2, 1000)
        self.y_fine = self.f_smooth(self.x_fine)

    def interpolate_smoothed_profile(self):
        x_ext = numpy.append(self.x, [self.x[-1] + 0.001, 10])
        y_ext = numpy.append(self.y_smooth, [0., 0.])
        self.f_smooth = scipy.interpolate.interp1d(x_ext, y_ext, kind='linear', bounds_error=False)

    def plot_profile(self):
        fig, ax = matplotlib.pyplot.subplots()
        ax.plot(self.x_raw, self.y_raw, '.')
        ax.plot(self.x, self.y, '.')
        ax.fill_between(self.x, self.y - 1 * self.y_error, self.y + 1 * self.y_error, color='gray', alpha=0.4)
        ax.fill_between(self.x, self.y - 3 * self.y_error, self.y + 3 * self.y_error, color='gray', alpha=0.3)
        ax.fill_between(self.x, self.y - 5 * self.y_error, self.y + 5 * self.y_error, color='gray', alpha=0.2)
        ax.plot(self.x_fine, self.y_fine)
        matplotlib.pyplot.show()

    def get_value(self, normalised_poloidal_flux):
        return self.f_smooth(normalised_poloidal_flux)
