import elastica
import numba
import numpy as np
from elastica.utils import _bspline
import copy
from matrix_operators import _batch_product_i_k_to_ik, _batch_matvec,inplace_addition,inplace_substraction
from elastica.external_forces import NoForces
from numba import njit, float64

class StretchAndTwitch(NoForces):
    def __init__(self, b_coeff, rest_length,period,wave_length,direction,percent_crawling = 0):
        self.b_coeff = np.array(b_coeff)
        self.total_force = 0
        self.rest_length = copy.deepcopy(rest_length)
        self.wave_length = wave_length
        self.period = period
        self.direction = direction
        self.angular_frequency = 2.0 * np.pi / period
        self.wave_number = 2.0 * np.pi / wave_length
        self.percent_crawling = percent_crawling
        self.s = np.cumsum(self.rest_length)
        self.s /= self.s[-1]
        spline, _, _ = _bspline(self.b_coeff)
        self.spline = spline(self.s)
        # self.node_s = np.insert(self.s,0,0)
        # self.node_spline = spline(self.node_s)

    def calc_total_force(self, time: np.float64 = 0.0):
        self.s = np.cumsum(self.rest_length)
        self.s /= self.s[-1]
        spline,_,_ = _bspline(self.b_coeff)
        self.spline = spline(self.s)
        # calculate total force of the snake
        self.total_force = self.spline * np.sin(-self.angular_frequency * time + self.wave_number * self.s)


    def apply_torques(self,system,time: np.float64 = 0.0):
        calc_twist(
            self.s,
            self.spline,
            self.angular_frequency,
            time,
            self.wave_number,
            self.percent_crawling,
            self.direction,
            system.radius,
            system.external_torques,
            system.director_collection
        )

    def apply_forces(self,system,time: np.float64 = 0.0):
        calc_stretch(
                self.s,
                self.spline,
                self.angular_frequency,
                time,
                self.wave_number,
                self.percent_crawling,
                system.n_elems,
                system.tangents,
                system.external_forces,
                )
@staticmethod
@njit(cache=True)
def calc_stretch(
            s,
            spline,
            angular_frequency,
            time,
            wave_number,
            percent_crawling,
            n_elems,
            tangents,
            external_forces,
            ):
    #calculate total force
    total_force = spline * np.sin(angular_frequency * time + wave_number * s)

    force_mag = total_force * percent_crawling
    force = np.zeros((3, n_elems),dtype=numba.float64)
    # we need to find the tangent of the snake at each nodes and multiply them with the stretch force
    for i in range(3):
        for element in range(n_elems):
            force[i, element] = force_mag[element] * tangents[i, element]

    inplace_addition(external_forces[:, :-1], force)
    inplace_substraction(external_forces[:, 1:], force)

@staticmethod
@njit(cache=True)
def calc_twist(
        s,
        spline,
        angular_frequency,
        time,
        wave_number,
        percent_crawling,
        radius,
        direction,
        external_torques,
        director_collection
        ):
    # calculate total force
    total_force = spline * np.sin(-angular_frequency * time + wave_number * s)
    force_mag = total_force * (1-percent_crawling)
    torque_mag = radius * force_mag
    torque = _batch_product_i_k_to_ik(direction, torque_mag[::-1])
    inplace_addition(
        external_torques[..., 1:],
        _batch_matvec(director_collection, torque)[..., 1:],
    )
    inplace_substraction(
        external_torques[..., :-1],
        _batch_matvec(director_collection[..., :-1], torque[..., 1:]),
    )

# b_coeff = [0,23,23,23,23,0]
# rest_length = np.array([0.1,0.1,0.1,0.1,0.1])
# period = 1
# wave_length = 1
# direction = np.array([[0,0,1]])
# time = 2
# test_case = StretchAndTwitch(b_coeff, rest_length, period, wave_length, direction)
# print(test_case.calc_total_force(time))