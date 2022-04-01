import os
import h5py
import numpy
import scipy.interpolate
from h5py import File


def get_home_directory():
    try:
        return os.environ['TAIGA_HOME']
    except KeyError:
        return '.'


class EFITDataReader:
    file: File

    def __init__(self, shot_number, reconstruction_id=1):
        self.shot_number = shot_number
        self.reconstruction_id = reconstruction_id
        self.file_path = ''
        self.init_efit_file()

    def get_data(self, field):
        try:
            return self.file[field][()]
        except KeyError:
            raise KeyError('Invalid hdf5 format.')

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


class EFITManager(EFITDataReader):
    def __init__(self, shot_number, time, reconstruction_id=1):
        super().__init__(shot_number, reconstruction_id)
        self.time = time
        self.time_index = self.get_time_index()

    def get_time_index(self):
        time_dataset = self.get_data('time')
        return (numpy.abs(time_dataset - int(self.time) / 1000)).argmin()

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
    psi_n: scipy.interpolate.RectBivariateSpline
    normalise_poloidal_flux: scipy.interpolate.UnivariateSpline

    def __init__(self, shot_number, time):
        self.shot_number = shot_number
        self.time = time
        self.efit = EFITManager(self.shot_number, self.time)

        self.R = self.efit.get_time_sliced_data('/output/profiles2D/r')
        self.one_over_R = numpy.diag(1/self.R)
        self.Z = self.efit.get_time_sliced_data('/output/profiles2D/z')

        self.poloidal_flux = []
        self.get_poloidal_flux()

        self.normalised_poloidal_flux = []
        self.get_normalised_poloidal_flux()

        self.set_psi()
        self.set_psi_n()
        self.set_B_R()
        self.set_B_Z()
        self.set_B_phi()

    def set_interpolation_for_poloidal_flux(self, y):
        poloidal_flux = self.efit.get_time_sliced_data('output/fluxFunctionProfiles/poloidalFlux')
        return scipy.interpolate.UnivariateSpline(poloidal_flux, y, k=3, s=0)

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

    def set_psi_n(self):
        self.psi_n = scipy.interpolate.RectBivariateSpline(self.Z, self.R, self.normalised_poloidal_flux)

    def set_B_R(self):
        B_R_2d = numpy.dot(-self.psi(self.Z, self.R, dx=1, dy=0, grid=True), self.one_over_R)
        self.B_R = scipy.interpolate.RectBivariateSpline(self.Z, self.R, B_R_2d)

    def set_B_Z(self):
        B_Z_2d = numpy.dot(self.psi(self.Z, self.R, dx=0, dy=1, grid=True), self.one_over_R)
        self.B_Z = scipy.interpolate.RectBivariateSpline(self.Z, self.R, B_Z_2d)

    def set_B_phi(self):
        R_B_phi_1d = self.efit.get_time_sliced_data('/output/fluxFunctionProfiles/rBphi')
        get_R_B_phi = self.set_interpolation_for_poloidal_flux(R_B_phi_1d)
        if len(self.poloidal_flux) == 0:
            self.get_poloidal_flux()
        R_B_phi_2d = get_R_B_phi(self.poloidal_flux)
        B_phi_2d = numpy.dot(-R_B_phi_2d, self.one_over_R)
        self.B_phi = scipy.interpolate.RectBivariateSpline(self.Z, self.R, B_phi_2d)


class ParseToTaiga:
    def __init__(self, cdb):
        export_dir = get_home_directory() + '/input/fieldSpl/' + str(cdb.shot_number) + '_' + str(cdb.time)
        print('Save B-spline coefficients to: ' + export_dir)
        try:
            os.mkdir(export_dir)
            print('Create directory and write data to ' + export_dir)
        except FileExistsError:
            print('Write data to ' + export_dir)
        else:
            pass
        numpy.savetxt(export_dir + '/z.bspl', cdb.B_R.tck[0])
        numpy.savetxt(export_dir + '/r.bspl', cdb.B_R.tck[1])
        numpy.savetxt(export_dir + '/brad.bspl', cdb.B_R.tck[2])
        numpy.savetxt(export_dir + '/bz.bspl', cdb.B_Z.tck[2])
        numpy.savetxt(export_dir + '/btor.bspl', cdb.B_phi.tck[2])
        numpy.savetxt(export_dir + '/psi_n.bspl', cdb.psi_n.tck[2])


if __name__ == "__main__":
    cr = CDBReader(17178, 1097)
    ParseToTaiga(cr)

