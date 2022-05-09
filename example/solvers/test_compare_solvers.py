import unittest
from taiga.example.solvers.compare_solvers import SolverData


class TestCompareSolvers(unittest):
    def test_timestep(self):
        s = SolverData('rk4')
        reference = '$10^{-9}$ s'
        self.assertEqual(reference, s.get_timestep())

    def test_timestep_1p5em10(self):
        s = SolverData('rk4')
        s.timestep = 1.5e-10
        reference = '$1.5 \cdot 10^{-10}$ s'
        self.assertEqual(reference, s.get_timestep())
