import os

import matplotlib.pyplot
import numpy
import scipy.interpolate


def set_negatives_to_zero(values):
    return values * (values > 0)


def transform_to_poloidal_cross_section(r, tor):
    return numpy.hypot(r, tor)


def get_home_directory():
    try:
        return os.environ['TAIGA_HOME']
    except KeyError:
        return '.'


class ProfileManager:
    def __init__(self, export_directory, x, y, y_error, field):
        self.export_directory = export_directory
        self.field = field
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
        from statsmodels.nonparametric.smoothers_lowess import lowess
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

    def plot_profile(self, shot_number, time, y_scale=1):
        fig, ax = matplotlib.pyplot.subplots()
        fig.set_size_inches(3.5, 3.5)
        fig.subplots_adjust(left=0.2, bottom=0.15)
        ax.plot(self.x_raw, self.y_raw*y_scale, '.')
        ax.plot(self.x, self.y*y_scale, '.')
        ax.fill_between(self.x, (self.y - 1 * self.y_error)*y_scale, (self.y + 1 * self.y_error)*y_scale, color='gray', alpha=0.4)
        ax.fill_between(self.x, (self.y - 3 * self.y_error)*y_scale, (self.y + 3 * self.y_error)*y_scale, color='gray', alpha=0.3)
        ax.fill_between(self.x, (self.y - 5 * self.y_error)*y_scale, (self.y + 5 * self.y_error)*y_scale, color='gray', alpha=0.2)
        ax.plot(self.x_fine, self.y_fine*y_scale)
        ax.set_xlabel('normalised poloidal flux')
        ax.set_ylabel(self.get_y_label(y_scale))
        ax.set_title('COMPASS #'+shot_number+' ('+time+' ms) ')
        matplotlib.pyplot.minorticks_on()
        matplotlib.pyplot.grid(which='both')
        matplotlib.pyplot.axvline(1, c='red', ls='--')
        matplotlib.pyplot.text(1.005, max(self.y), 'LCFS', c='red', fontsize=12)
        matplotlib.pyplot.savefig(self.export_directory + '/' + self.field + '.pdf')
        matplotlib.pyplot.show()

    def get_y_label(self, y_scale=1):
        log_y_scale = int(numpy.log10(y_scale))
        scale_text = '' if y_scale == 1 else '10 ^ {' + str(log_y_scale) + '} ~'
        if self.field == 'density':
            return r'$n_e~[' + scale_text + '\mathrm{m}^{-3}]$'
        elif self.field == 'temperature':
            return r'$T_e~[' + scale_text + '\mathrm{eV}]$'
        else:
            raise ValueError

    def get_value(self, normalised_poloidal_flux):
        return self.f_smooth(normalised_poloidal_flux)

    def set_path(self, path):
        try:
            os.mkdir(path)
            print('Create directory and write data to ' + path)
        except FileExistsError:
            print('Write data to ' + path)
        else:
            pass

    def export_profile(self):
        flux_file = open(self.export_directory + '/flux.prof', 'w')
        # noinspection PyTypeChecker
        numpy.savetxt(flux_file, self.x_fine)
        flux_file.close()

        data_file = open(self.export_directory + '/' + self.field + '.prof', 'w')
        # noinspection PyTypeChecker
        numpy.savetxt(data_file, self.y_fine)
        data_file.close()
