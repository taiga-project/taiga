import numpy
import scipy
# from scipy.optimize import curve_fit
from scipy import odr


def get_fitted_profile(self):
    fit = scipy.optimize.curve_fit(stefanikova_ped_old, self.x, self.y,
                                   p0=[self.y[0], self.y[-1], 1, 0.02, self.y[0] / 1000])
    p = fit[0]
    self.y_fit = stefanikova_ped(p, self.x)

def stefanikova_core(a, r):
    a_height, a_width, a_exp = a
    return a_height * (numpy.exp(-r ** 2 / a_width) ** a_exp)


def stefanikova_core_old(r, a_height, a_width, a_exp):
    return a_height * (numpy.exp(-r ** 2 / a_width) ** a_exp)


def stefanikova_ped_old(r, b_height, b_sol, b_pos, b_width, b_slope):
    return (b_height - b_sol) / 2 * (mtanh((b_pos - r) / (2 * b_width), b_slope) + 1) + b_sol


def stefanikova_ped(b, r):
    b_height, b_sol, b_pos, b_width, b_slope = b
    return (b_height - b_sol) / 2 * (mtanh((b_pos - r) / (2 * b_width), b_slope) + 1) + b_sol


def stefanikova_full(ab, r):
    a_height, a_width, a_exp, b_height, b_sol, b_pos, b_width, b_slope = ab
    F_ped = (b_height - b_sol) / 2 * (mtanh((b_pos - r) / (2 * b_width), b_slope) + 1) + b_sol
    return F_ped + (a_height - F_ped) * numpy.exp((-r / a_width) ** a_exp)


def mtanh(x, b_slope):
    return ((1 + b_slope * x) * scipy.exp(x) - scipy.exp(-x)) / (scipy.exp(x) + scipy.exp(-x))


def get_fitted_profile(x, y, yerr):
    data = scipy.odr.RealData(x, y, sy=yerr)

    ped_index = numpy.where((x > 0.7) & (x < 1.04))
    core_index = numpy.where((x > 0) & (x < 0.3))
    core_data = scipy.odr.RealData(x[core_index], y[core_index], sy=yerr)
    ped_data = scipy.odr.RealData(x[ped_index], y[ped_index], sy=yerr)

    stefanikova_core_model = scipy.odr.Model(stefanikova_core)
    stefanikova_core_odr = scipy.odr.ODR(core_data, stefanikova_core_model,
                                         beta0=[y[0], 0.5, 4])
    stefanikova_a = stefanikova_core_odr.run().beta

    stefanikova_ped_model = scipy.odr.Model(stefanikova_ped)
    stefanikova_ped_odr = scipy.odr.ODR(ped_data, stefanikova_ped_model,
                                        beta0=[y[0]/2, y[0] / 50, 1., 0.02, y[0] / 1000],
                                        maxit=10000)
    stefanikova_b = stefanikova_ped_odr.run().beta

    stefanikova_a = [y[0], 0.5, 2.]
    stefanikova_full_model = scipy.odr.Model(stefanikova_full)
    stefanikova_full_odr = scipy.odr.ODR(data, stefanikova_full_model,
                                         beta0=numpy.concatenate((stefanikova_a, stefanikova_b)))
    stefanikova_ab = stefanikova_full_odr.run().beta

    return stefanikova_ped(stefanikova_b, x)
