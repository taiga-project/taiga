import numpy
import matplotlib.pyplot


class SolverData:
    directory = "bin/example"
    reference_velocity = 400000
    number_of_cyclotron_periods = 10000000
    frequency_of_export = 100000
    time_step = 1e-9
    field_scenario = 'B = 1'

    def __init__(self, solver_id, scenario='default'):
        self.solver_file = solver_id
        self.solver_name = self.get_solver_name(solver_id)
        self.set_scenario(scenario)
        self.data = self.get_data()
        self.x = self.get_field(0)
        self.y = self.get_field(1)
        self.v = numpy.sqrt(self.get_field(3)**2 + self.get_field(4)**2 + self.get_field(5)**2)
        self.v_relative = self.v / self.reference_velocity
        self.x_axis = numpy.arange(0, self.number_of_cyclotron_periods * 2 + 1, self.frequency_of_export)
        self.line_style = self.get_line_style(solver_id)

    def set_scenario(self, scenario):
        if scenario != 'default':
            self.solver_file += '_' + scenario

        if scenario == 'start10':
            self.time_step = 1e-10

        if scenario in {'start8', 'gradb_start8'}:
            self.time_step = 1e-8

        if scenario in {'gradb', 'gradb_start'}:
            self.field_scenario = 'B_z = 1 + 0.01 y'

        if scenario in {'eparb', 'eparb_start'}:
            self.field_scenario = 'B_z = 1, E_z = 0.01'

        if scenario in {'b__r'}:
            self.field_scenario = '$B_z = \sqrt{x^2+y^2}$'#, E_z = 0.01'

        if scenario in {'start', 'gradb_start', 'b__r_start'}:
            self.reference_velocity = 400000
            self.number_of_cyclotron_periods = 100#000
            self.frequency_of_export = 1#000

        if scenario in {'start8', 'gradb_start8'}:
            self.reference_velocity = 400000
            self.number_of_cyclotron_periods = 1000  # 000
            self.frequency_of_export = 1  # 000

        if scenario in {'b__r'}:
            self.reference_velocity = 400000
            self.number_of_cyclotron_periods = 10000000
            self.frequency_of_export = 100

    @staticmethod
    def get_solver_name(solver_id):
        if solver_id == 'rk4':
            return 'Runge--Kutta 4/5'
        elif solver_id == 'rkn':
            return 'Runge--Kutta--Nyström 4/5'
        elif solver_id == 'verlet':
            return 'Verlet--Boris 2/3'
        elif solver_id == 'yoshida':
            return 'Yoshida--Boris 4/5'
        else:
            return 'unknown'

    @staticmethod
    def get_line_style(solver_id):
        if solver_id == 'rk4':
            return '-'
        elif solver_id == 'rkn':
            return '--'
        elif solver_id == 'verlet':
            return '-'
        elif solver_id == 'yoshida':
            return '--'
        else:
            return ':'

    def get_data(self):
        return numpy.loadtxt(self.directory + '/' + self.solver_file + '.dat')

    def get_field(self, identifier):
        return self.data[:, identifier]

    def plot_xy(self):
        matplotlib.pyplot.plot(self.x, self.y, self.line_style, markersize=1.5, label=self.solver_name)
        matplotlib.pyplot.title('timestep: ' + str(self.time_step) + ', ' + self.field_scenario)

    def plot_v(self):
        matplotlib.pyplot.plot(self.x_axis, self.v_relative, self.line_style, linewidth=2.5, label=self.solver_name)
        matplotlib.pyplot.title('timestep: ' + str(self.time_step) + ', ' + self.field_scenario)


def visualise():
    matplotlib.pyplot.xlabel('steps')
    matplotlib.pyplot.ylabel('relative velocity')
    matplotlib.pyplot.legend()
    matplotlib.pyplot.tight_layout()


def visualise_trajectory():
    matplotlib.pyplot.xlabel('x [m]')
    matplotlib.pyplot.ylabel('y [m]')
    matplotlib.pyplot.axis('square')
    matplotlib.pyplot.legend()
    matplotlib.pyplot.tight_layout()


def compare_solvers(scenario='default'):
    SolverData('rk4', scenario).plot_v()
    SolverData('rkn', scenario).plot_v()
    SolverData('verlet', scenario).plot_v()
    SolverData('yoshida', scenario).plot_v()
    visualise()
    matplotlib.pyplot.savefig('bin/example/' + scenario + '.pdf')
    matplotlib.pyplot.savefig('bin/example/' + scenario + '.png')
    matplotlib.pyplot.show()


def compare_solvers_trajectory(scenario='default'):
    SolverData('rk4', scenario).plot_xy()
    SolverData('rkn', scenario).plot_xy()
    SolverData('verlet', scenario).plot_xy()
    SolverData('yoshida', scenario).plot_xy()
    visualise_trajectory()
    matplotlib.pyplot.savefig('bin/example/' + scenario + '_traj.pdf')
    matplotlib.pyplot.savefig('bin/example/' + scenario + '_traj.png')
    matplotlib.pyplot.show()


if __name__ == "__main__":
    compare_solvers('b__r')
    compare_solvers_trajectory('b__r')
    compare_solvers()
    compare_solvers_trajectory()
    compare_solvers_trajectory('start')
    compare_solvers('gradb')
    compare_solvers_trajectory('gradb')
    compare_solvers('eparb_start')
    compare_solvers('eparb')
    compare_solvers('gradb_start8')
def aa():
    compare_solvers('eparb_start')
    compare_solvers('start10')
    compare_solvers('start8')
    compare_solvers_trajectory('gradb')
    compare_solvers_trajectory('br')
    compare_solvers('br')
