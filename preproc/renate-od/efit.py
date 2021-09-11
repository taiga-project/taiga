import scipy
import h5py

from utils import *


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

    def get_normalised_poloidal_flux(self, beamlet_geometry, grid=False):
        r_poloidal_cross_section = transform_to_poloidal_cross_section(beamlet_geometry.rad, beamlet_geometry.tor)
        poloidal_flux = self.get_poloidal_flux(r_poloidal_cross_section, beamlet_geometry.z, grid=grid)
        return self.normalise_poloidal_flux(poloidal_flux)
