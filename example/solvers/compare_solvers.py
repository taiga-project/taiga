import numpy
import matplotlib.pyplot


class SolverData:
    directory = "bin/example"

    def __init__(self, solver_name):
        self.solver_name = solver_name
        self.data = self.get_data()
        self.x = self.get_field(0)
        self.y = self.get_field(1)
        self.v = numpy.sqrt(self.get_field(3)**2 + self.get_field(4)**2 + self.get_field(5)**2)

    def get_data(self):
        return numpy.loadtxt(self.directory + '/' + self.solver_name + '.dat')

    def get_field(self, identifier):
        return self.data[:, identifier]

    def plot(self):
        matplotlib.pyplot.plot(self.x, self.y, '.', label=self.solver_name)
        matplotlib.pyplot.show()

    def plot_v(self):
        matplotlib.pyplot.plot(self.v, '.', label=self.solver_name)
        matplotlib.pyplot.show()


def compare_solvers():
    SolverData('rk4').plot_v()
    SolverData('rkn').plot_v()
    SolverData('verlet').plot_v()
    SolverData('yoshida').plot_v()


if __name__ == "__main__":
    compare_solvers()
