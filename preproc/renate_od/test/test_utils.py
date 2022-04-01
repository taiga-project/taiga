import unittest

from taiga.preproc.renate_od.utils import ProfileManager


class MockProfileManagerNoInit(ProfileManager):
    def __init__(self, x, y, y_error):
        pass


class TestRODUtils(unittest.TestCase):
    def test_get_y_label_density(self):
        x = [0, 1]
        y = [0, 1]
        yerr = [0, 0]
        x = None
        y = None
        yerr = None
        p = MockProfileManagerNoInit(x, y, yerr)
        p.field = 'density'
        reference = '$n_e~[10 ^ {-19} ~\mathrm{m}^{-3}]$'
        self.assertEqual(reference, p.get_y_label(1e-19))

    def test_get_y_label_temperature(self):
        x = [0, 1]
        y = [0, 1]
        yerr = [0, 0]
        x = None
        y = None
        yerr = None
        p = MockProfileManagerNoInit(x, y, yerr)
        p.field = 'temperature'
        reference = '$T_e~[\mathrm{eV}]$'
        self.assertEqual(reference, p.get_y_label())


if __name__ == '__main__':
    unittest.main()
