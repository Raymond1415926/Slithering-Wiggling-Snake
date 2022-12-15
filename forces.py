import elastica
import numpy as np
from elastica.utils import _bspline
import copy
from matrix_operators import _batch_product_i_k_to_ik, _batch_matvec,inplace_addition,inplace_substraction
elastica.CosseratRod

class StretchAndTwitch():
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

    def calc_total_force(self, time: np.float64 = 0.0):
        self.s = np.cumsum(self.rest_length)
        self.s /= self.s[-1]
        spline,_,_ = _bspline(self.b_coeff)
        self.spline = spline(self.s)
        # calculate total force of the snake
        self.total_force = self.spline * np.sin(-self.angular_frequency * time + self.wave_number * self.s)


    def apply_torques(self,system,time: np.float64 = 0.0):
        #calculate the twitching torque
        self.calc_total_force(time)
        force_mag = self.total_force * (1-self.percent_crawling)
        torque_mag = system.base_radius * force_mag
        torque = _batch_product_i_k_to_ik(self.direction, torque_mag[::-1])
        inplace_addition(
            system.external_torques[..., 1:],
            _batch_matvec(system.director_collection, torque)[..., 1:],
        )
        inplace_substraction(
            system.external_torques[..., :-1],
            _batch_matvec(system.director_collection[..., :-1], torque[..., 1:]),
        )

    def apply_forces(self,system,time: np.float64 = 0.0):
        #calculate the stretching force
        self.calc_total_force(time)
        force_mag = self.total_force * self.percent_crawling
        #we need to find the tangent of the snake at each nodes and multiply them with the stretch force
        inplace_addition(system.external_forces[:,1:], 0.5 * np.multiply(system.tangents, force_mag[1:]))
        inplace_addition(system.external_forces[:, :-1], 0.5 * np.multiply(system.tangents, force_mag[:-1]))

# b_coeff = [0,23,23,23,23,0]
# rest_length = np.array([0.1,0.1,0.1,0.1,0.1])
# period = 1
# wave_length = 1
# direction = np.array([[0,0,1]])
# time = 2
# test_case = StretchAndTwitch(b_coeff, rest_length, period, wave_length, direction)
# print(test_case.calc_total_force(time))