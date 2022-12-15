import numpy as np
from elastica import external_forces
from elastica.utils import _bspline
import copy

class StretchAndTwitch():
    def __init__(self, b_coeff, rest_length,period,wave_length,direction):
        self.b_coeff = b_coeff
        self.total_force = 0
        self.rest_length = rest_length.copy()
        self.wave_length = wave_length
        self.period = period
        self.direction = direction
        self.angular_frequency = 2.0 * np.pi / period
        self.wave_number = 2.0 * np.pi / wave_length

    def calc_total_force(self, time: np.float64 = 0.0):
        self.s = np.cumsum(self.rest_length)
        self.s /= self.s[-1]
        spline = _bspline(self.b_coeff)
        self.spline = spline(self.s)

        # calculate torque magnitude
        magnitude = self.spline * np.sin(-self.angular_frequency * time + self.wave_number * self.s)

        return magnitude

b_coeff = [0,23,23,23,23,0]
rest_length = np.array([0.1,0.1,0.1,0.1,0.1])
period = 1
wave_length = 1
direction = np.array([[0,0,1]])
time = 2
test_case = StretchAndTwitch(b_coeff, rest_length, period, wave_length, direction)
print(test_case.calc_total_force(time))