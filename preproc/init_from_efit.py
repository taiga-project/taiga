import os
import h5py
import numpy
import scipy.interpolate
from h5py import File
import matplotlib.pyplot as plt


def get_home_directory():
    try:
        return os.environ['TAIGA_HOME']
    except KeyError:
        return '.'


class EFITManager:
    file: File

    def __init__(self, shot_number, time, reconstruction_id=1):
        self.shot_number = shot_number
        self.time = time
        self.reconstruction_id = reconstruction_id
        self.file_path = ''
        self.init_efit_file()
        self.time_index = self.get_time_index()

    def set_efit_file(self, file_path):
        self.file_path = file_path
        self.open_efit_file()

    def open_efit_file(self):
        try:
            self.file = h5py.File(self.file_path, 'r')
            print('Reading data from: ' + self.file_path)
        except OSError:
            raise

    def init_efit_file(self):
        self.set_efit_file(get_home_directory() + '/input/cdb/' + str(self.shot_number) +
                           '/EFITXX/EFITXX.' + str(self.reconstruction_id) + '.h5')

    def get_time_index(self):
        time_dataset = self.get_data('time')
        return (numpy.abs(time_dataset - int(self.time) / 1000)).argmin()

    def get_data(self, field):
        try:
            return self.file[field][()]
        except KeyError:
            raise KeyError('Invalid hdf5 format.')

    def get_time_sliced_data(self, field):
        try:
            return self.file[field][self.time_index]
        except KeyError:
            raise KeyError('Invalid hdf5 format.')


class CDBReader:
    B_R: scipy.interpolate.RectBivariateSpline
    B_Z: scipy.interpolate.RectBivariateSpline
    B_phi: scipy.interpolate.RectBivariateSpline
    psi: scipy.interpolate.RectBivariateSpline
    normalise_poloidal_flux: scipy.interpolate.UnivariateSpline

    def __init__(self, shot_number, time):
        self.shot_number = shot_number
        self.time = time
        self.efit = EFITManager(self.shot_number, self.time)

        self.R = self.efit.get_time_sliced_data('/output/profiles2D/r')
        self.Z = self.efit.get_time_sliced_data('/output/profiles2D/z')

        self.poloidal_flux = []
        self.get_poloidal_flux()

        self.normalised_poloidal_flux = []
        self.get_normalised_poloidal_flux()

        self.set_psi()
        self.set_B_R()
        self.set_B_Z()
        self.set_B_phi()

    def set_interpolation_for_poloidal_flux(self, y):
        poloidal_flux = self.efit.get_time_sliced_data('output/fluxFunctionProfiles/poloidalFlux')
        return scipy.interpolate.UnivariateSpline(poloidal_flux, y)

    def set_normalised_poloidal_flux(self):
        normalised_poloidal_flux = self.efit.get_data('output/fluxFunctionProfiles/normalizedPoloidalFlux')
        self.normalise_poloidal_flux = self.set_interpolation_for_poloidal_flux(normalised_poloidal_flux)

    def get_normalised_poloidal_flux(self):
        self.set_normalised_poloidal_flux()
        if len(self.poloidal_flux) == 0:
            self.get_poloidal_flux()
        self.normalised_poloidal_flux = self.normalise_poloidal_flux(self.poloidal_flux)

    def get_poloidal_flux(self):
        self.poloidal_flux = numpy.transpose(self.efit.get_time_sliced_data('/output/profiles2D/poloidalFlux'))

    def set_psi(self):
        self.psi = scipy.interpolate.RectBivariateSpline(self.Z, self.R, self.poloidal_flux, kx=5, ky=5)

    def set_B_R(self):
        B_R_2d = -self.psi(self.Z, self.R, dx=1, dy=0, grid=True)/self.R
        self.B_R = scipy.interpolate.RectBivariateSpline(self.Z, self.R, B_R_2d)

    def set_B_Z(self):
        B_Z_2d = self.psi(self.Z, self.R, dx=0, dy=1, grid=True)/self.R
        self.B_Z = scipy.interpolate.RectBivariateSpline(self.Z, self.R, B_Z_2d)

    def set_B_phi(self):
        R_B_phi_1d = self.efit.get_time_sliced_data('/output/fluxFunctionProfiles/rBphi')
        get_R_B_phi = self.set_interpolation_for_poloidal_flux(R_B_phi_1d)
        if len(self.poloidal_flux) == 0:
            self.get_poloidal_flux()
        R_B_phi_2d = get_R_B_phi(self.poloidal_flux)
        B_phi_2d = -R_B_phi_2d/self.R
        self.B_phi = scipy.interpolate.RectBivariateSpline(self.Z, self.R, B_phi_2d)


class PlotCDB:
    def __init__(self, cdb):
        self.cdb = cdb
        self.value_range = [0, 1]
        self.plot('B_R')
        self.plot('B_Z')
        self.plot('B_phi')
        self.plot('psi')

    def plot(self, attribute):
        self.get_range(attribute)
        f = getattr(self.cdb, attribute)

        fig, (ax_taiga, ax_reference) = plt.subplots(1, 2, sharex='all', sharey='all',
                                                     subplot_kw=dict(aspect='equal'))

        ax_reference.tick_params(direction='out', left=True, right=True, labelleft=True)
        value = f(self.cdb.Z, self.cdb.R, grid=True)
        c = ax_reference.contourf(self.cdb.R, self.cdb.Z, value, levels=self.get_levels())
        ax_reference.set_xlabel('R [m]')
        ax_reference.tick_params(direction='out', left=True, right=True)
        ax_reference.set_title('CDB')
        plt.suptitle(attribute + ' @ COMPASS #' + self.cdb.shot_number + ' (' + self.cdb.time + ' ms) ')
        if attribute == 'psi':
            ax_reference.tick_params(direction='out', left=True, right=True, labelright=True)
        else:
            cbar = fig.colorbar(c, ax=(ax_taiga, ax_reference), pad=0.02, fraction=0.034)
            cbar.set_ticks(self.get_tick_levels())

        fig.show()

    def get_levels(self):
        return numpy.linspace(self.value_range[0], self.value_range[1], 21)

    def get_tick_levels(self):
        return numpy.linspace(self.value_range[0], self.value_range[1], 9)

    def get_range(self, attribute):
        if attribute == 'B_R':
            self.value_range = [-0.2, 0.2]
        elif attribute == 'B_Z':
            self.value_range = [-0.4, 0.4]
        elif attribute == 'B_phi':
            self.value_range = [0, 2]
        elif attribute == 'psi':
            self.value_range = [-0.02, 0.02]
        elif attribute == 'normalised_poloidal_flux':
            self.value_range = [0, 1]


class ParseToTaiga:
    def __init__(self, cdb):
        export_dir = get_home_directory() + '/input/fieldSpl/' + str(cdb.shot_number) + '_' + str(cdb.time)
        print('Save B-spline coefficients to: ' + export_dir)
        numpy.savetxt(export_dir + '/z.bspl', cdb.B_R.tck[0])
        numpy.savetxt(export_dir + '/r.bspl', cdb.B_R.tck[1])
        numpy.savetxt(export_dir + '/brad.bspl', cdb.B_R.tck[2])
        numpy.savetxt(export_dir + '/bz.bspl', cdb.B_Z.tck[2])
        numpy.savetxt(export_dir + '/btor.bspl', cdb.B_phi.tck[2])
        numpy.savetxt(export_dir + '/polflux.bspl', cdb.psi.tck[2])


if __name__ == "__main__":
    cr = CDBReader(17178, 1097)
    #PlotCDB(cr)
    ParseToTaiga(cr)
