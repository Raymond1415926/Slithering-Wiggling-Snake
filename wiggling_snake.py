__doc__ = """Snake friction case from X. Zhang et. al. Nat. Comm. 2021"""

import os
import numpy as np
from elastica import *
from plotting import (
    plot_snake_velocity,
    plot_video,
    compute_projected_velocity,
    plot_curvature,
)
from forces import*
from math import isnan



class SnakeSimulator(BaseSystemCollection, Constraints, Forcing, Damping, CallBacks):
    pass


def run_snake_pure_wiggle(
    b_coeff, wave_length, run_time=1, n_elements=10,PLOT_FIGURE=False, SAVE_FIGURE=False, SAVE_VIDEO=False,\
        xlim = (0,5)
):
    # Initialize the simulation class
    snake_sim = SnakeSimulator()

    # Simulation parameters
    period = 1
    final_time = run_time

    # setting up test params
    n_elem = n_elements
    start = np.array([0.0, 0.0, 0.0])
    direction = np.array([0.0, 0.0, 1.0])
    normal = np.array([0.0, 1.0, 0.0])
    base_length = 1
    base_radius = 0.025
    base_area = np.pi * base_radius ** 2
    density = 1000
    E = 1e6
    poisson_ratio = 0.5
    shear_modulus = E / (poisson_ratio + 1.0)

    shearable_rod = CosseratRod.straight_rod(
        n_elem,
        start,
        direction,
        normal,
        base_length,
        base_radius,
        density,
        nu=0.0,
        youngs_modulus=E,
        shear_modulus=shear_modulus,
    )

    snake_sim.append(shearable_rod)

    # Add gravitational forces
    gravitational_acc = -9.80665
    snake_sim.add_forcing_to(shearable_rod).using(
        GravityForces, acc_gravity=np.array([0.0, gravitational_acc, 0.0])
    )

    # Add muscle torques
    # snake_sim.add_forcing_to(shearable_rod).using(
    #     MuscleTorques,
    #     base_length=base_length,
    #     b_coeff= np.array(b_coeff),
    #     period=period,
    #     wave_number=2.0 * np.pi / (wave_length),
    #     phase_shift=0.0,
    #     rest_lengths=shearable_rod.rest_lengths,
    #     ramp_up_time=period,
    #     direction=normal,
    #     with_spline=True,
    # )


    #add wiggling
    """
    add wiggling/stretching
    """
    snake_sim.add_forcing_to(shearable_rod).using(
        StretchAndTwitch,
        b_coeff=b_coeff,
        rest_length=shearable_rod.rest_lengths,
        period=period,
        wave_length=wave_length,
        direction=normal,
        percent_crawling=1
    )
    # Add friction forces
    origin_plane = np.array([0.0, -base_radius, 0.0])
    normal_plane = normal
    slip_velocity_tol = 1e-8
    mu = 1.1019368
    kinetic_mu_array = np.array(
        [mu, 1.5*mu, 2*mu]
    )  # [forward, backward, sideways]
    static_mu_array = 2 * kinetic_mu_array
    snake_sim.add_forcing_to(shearable_rod).using(
        AnisotropicFrictionalPlane,
        k=1.0,
        nu=1e-6,
        plane_origin=origin_plane,
        plane_normal=normal_plane,
        slip_velocity_tol=slip_velocity_tol,
        static_mu_array=static_mu_array,
        kinetic_mu_array=kinetic_mu_array,
    )

    damping_constant = 1.9
    time_step = 1e-4
    snake_sim.dampen(shearable_rod).using(
        AnalyticalLinearDamper,
        damping_constant=damping_constant,
        time_step=time_step,
    )

    total_steps = int(final_time / time_step)
    rendering_fps = 60
    step_skip = int(1.0 / (rendering_fps * time_step))

    # Add call backs
    class ContinuumSnakeCallBack(CallBackBaseClass):
        """
        Call back function for continuum snake
        """

        def __init__(self, step_skip: int, callback_params: dict):
            CallBackBaseClass.__init__(self)
            self.every = step_skip
            self.callback_params = callback_params

        def make_callback(self, system, time, current_step: int):

            if current_step % self.every == 0:

                self.callback_params["time"].append(time)
                self.callback_params["step"].append(current_step)
                self.callback_params["position"].append(
                    system.position_collection.copy()
                )
                self.callback_params["velocity"].append(
                    system.velocity_collection.copy()
                )
                self.callback_params["avg_velocity"].append(
                    system.compute_velocity_center_of_mass()
                )

                self.callback_params["center_of_mass"].append(
                    system.compute_position_center_of_mass()
                )
                self.callback_params["curvature"].append(system.kappa.copy())

                return

    pp_list = defaultdict(list)
    snake_sim.collect_diagnostics(shearable_rod).using(
        ContinuumSnakeCallBack, step_skip=step_skip, callback_params=pp_list
    )

    snake_sim.finalize()

    timestepper = PositionVerlet()
    integrate(timestepper, snake_sim, final_time, total_steps)

    if PLOT_FIGURE:
        filename_plot = "continuum_snake_velocity.png"
        plot_snake_velocity(pp_list, period, filename_plot, SAVE_FIGURE)
        plot_curvature(pp_list, shearable_rod.rest_lengths, period, SAVE_FIGURE)

    if SAVE_VIDEO:
        filename_video = "continuum_snake.mp4"
        plot_video(
            pp_list,
            video_name=filename_video,
            fps=rendering_fps,
            xlim=xlim,
            ylim=(-1, 1),
        )


    distance_traveled = (pp_list["center_of_mass"][-1] - pp_list["center_of_mass"][0])[2]
    print(distance_traveled)
    return distance_traveled

def run_snake_pure_twitching(
    b_coeff, wave_length, run_time=1, n_elements=10,PLOT_FIGURE=False, SAVE_FIGURE=False, SAVE_VIDEO=False,\
        xlim = (0,5)
):
    # Initialize the simulation class
    snake_sim = SnakeSimulator()

    # Simulation parameters
    period = 1
    final_time = run_time

    # setting up test params
    n_elem = n_elements
    start = np.array([0.0, 0.0, 0.0])
    direction = np.array([0.0, 0.0, 1.0])
    normal = np.array([0.0, 1.0, 0.0])
    base_length = 1
    base_radius = 0.025
    base_area = np.pi * base_radius ** 2
    density = 1000
    E = 1e6
    poisson_ratio = 0.5
    shear_modulus = E / (poisson_ratio + 1.0)

    shearable_rod = CosseratRod.straight_rod(
        n_elem,
        start,
        direction,
        normal,
        base_length,
        base_radius,
        density,
        nu=0.0,
        youngs_modulus=E,
        shear_modulus=shear_modulus,
    )

    snake_sim.append(shearable_rod)

    # Add gravitational forces
    gravitational_acc = -9.80665
    snake_sim.add_forcing_to(shearable_rod).using(
        GravityForces, acc_gravity=np.array([0.0, gravitational_acc, 0.0])
    )

    #Add muscle torques
    # snake_sim.add_forcing_to(shearable_rod).using(
    #     MuscleTorques,
    #     base_length=base_length,
    #     b_coeff= np.array(b_coeff),
    #     period=period,
    #     wave_number=2.0 * np.pi / (wave_length),
    #     phase_shift=0.0,
    #     rest_lengths=shearable_rod.rest_lengths,
    #     ramp_up_time=period,
    #     direction=normal,
    #     with_spline=True,
    # )

    # Add friction forces
    """
    add twitching
    """
    snake_sim.add_forcing_to(shearable_rod).using(
        StretchAndTwitch,
        b_coeff=b_coeff,
        rest_length=shearable_rod.rest_lengths,
        period=period,
        wave_length=wave_length,
        direction=normal,
        percent_crawling=0
    )
    origin_plane = np.array([0.0, -base_radius, 0.0])
    normal_plane = normal
    slip_velocity_tol = 1e-8
    mu = 1.1019368
    kinetic_mu_array = np.array(
        [mu, 1.5*mu, 2*mu]
    )  # [forward, backward, sideways]
    static_mu_array = 2 * kinetic_mu_array
    snake_sim.add_forcing_to(shearable_rod).using(
        AnisotropicFrictionalPlane,
        k=1.0,
        nu=1e-6,
        plane_origin=origin_plane,
        plane_normal=normal_plane,
        slip_velocity_tol=slip_velocity_tol,
        static_mu_array=static_mu_array,
        kinetic_mu_array=kinetic_mu_array,
    )

    damping_constant = 1.9
    time_step = 1e-4
    snake_sim.dampen(shearable_rod).using(
        AnalyticalLinearDamper,
        damping_constant=damping_constant,
        time_step=time_step,
    )

    total_steps = int(final_time / time_step)
    rendering_fps = 60
    step_skip = int(1.0 / (rendering_fps * time_step))

    # Add call backs
    class ContinuumSnakeCallBack(CallBackBaseClass):
        """
        Call back function for continuum snake
        """

        def __init__(self, step_skip: int, callback_params: dict):
            CallBackBaseClass.__init__(self)
            self.every = step_skip
            self.callback_params = callback_params

        def make_callback(self, system, time, current_step: int):

            if current_step % self.every == 0:

                self.callback_params["time"].append(time)
                self.callback_params["step"].append(current_step)
                self.callback_params["position"].append(
                    system.position_collection.copy()
                )
                self.callback_params["velocity"].append(
                    system.velocity_collection.copy()
                )
                self.callback_params["avg_velocity"].append(
                    system.compute_velocity_center_of_mass()
                )

                self.callback_params["center_of_mass"].append(
                    system.compute_position_center_of_mass()
                )
                self.callback_params["curvature"].append(system.kappa.copy())

                return

    pp_list = defaultdict(list)
    snake_sim.collect_diagnostics(shearable_rod).using(
        ContinuumSnakeCallBack, step_skip=step_skip, callback_params=pp_list
    )

    snake_sim.finalize()

    timestepper = PositionVerlet()
    integrate(timestepper, snake_sim, final_time, total_steps)

    if PLOT_FIGURE:
        filename_plot = "continuum_snake_velocity.png"
        plot_snake_velocity(pp_list, period, filename_plot, SAVE_FIGURE)
        plot_curvature(pp_list, shearable_rod.rest_lengths, period, SAVE_FIGURE)

    if SAVE_VIDEO:
        filename_video = "continuum_snake.mp4"
        plot_video(
            pp_list,
            video_name=filename_video,
            fps=rendering_fps,
            xlim=xlim,
            ylim=(-1, 1),
        )


    distance_traveled = (pp_list["center_of_mass"][-1] - pp_list["center_of_mass"][0])[2]
    print(distance_traveled)
    return distance_traveled

def run_snake(
    b_coeff, wave_length, percent_crawling=0.5, run_time=1, n_elements=10,PLOT_FIGURE=False, SAVE_FIGURE=False, SAVE_VIDEO=False,\
        xlim = (0,5)
):
    # Initialize the simulation class
    snake_sim = SnakeSimulator()

    # Simulation parameters
    period = 1
    final_time = run_time

    # setting up test params
    n_elem = n_elements
    start = np.array([0.0, 0.0, 0.0])
    direction = np.array([0.0, 0.0, 1.0])
    normal = np.array([0.0, 1.0, 0.0])
    base_length = 1
    base_radius = 0.025
    base_area = np.pi * base_radius ** 2
    density = 1000
    E = 1e6
    poisson_ratio = 0.5
    shear_modulus = E / (poisson_ratio + 1.0)

    shearable_rod = CosseratRod.straight_rod(
        n_elem,
        start,
        direction,
        normal,
        base_length,
        base_radius,
        density,
        nu=0.0,
        youngs_modulus=E,
        shear_modulus=shear_modulus,
    )

    snake_sim.append(shearable_rod)

    # Add gravitational forces
    gravitational_acc = -9.80665
    snake_sim.add_forcing_to(shearable_rod).using(
        GravityForces, acc_gravity=np.array([0.0, gravitational_acc, 0.0])
    )

    percent_crawling = percent_crawling
    # Add muscle torques
    # snake_sim.add_forcing_to(shearable_rod).using(
    #     MuscleTorques,
    #     base_length=base_length,
    #     b_coeff= np.array(b_coeff) * percent_crawling,
    #     period=period,
    #     wave_number=2.0 * np.pi / (wave_length),
    #     phase_shift=0.0,
    #     rest_lengths=shearable_rod.rest_lengths,
    #     ramp_up_time=period,
    #     direction=normal,
    #     with_spline=True,
    # )
    #Add wiggling

    """
    Add both stretching and twitching
    """
    snake_sim.add_forcing_to(shearable_rod).using(
        StretchAndTwitch,
        b_coeff=b_coeff,
        rest_length=shearable_rod.rest_lengths,
        period=period,
        wave_length=wave_length,
        direction=normal,
        percent_crawling=percent_crawling
    )

    # Add friction forces
    origin_plane = np.array([0.0, -base_radius, 0.0])
    normal_plane = normal
    slip_velocity_tol = 1e-8
    mu = 1.1019368
    kinetic_mu_array = np.array(
        [mu, 1.5*mu, 2*mu]
    )  # [forward, backward, sideways]
    # static_mu_array = np.zeros(kinetic_mu_array.shape)
    static_mu_array = 2 * kinetic_mu_array
    snake_sim.add_forcing_to(shearable_rod).using(
        AnisotropicFrictionalPlane,
        k=1.0,
        nu=1e-6,
        plane_origin=origin_plane,
        plane_normal=normal_plane,
        slip_velocity_tol=slip_velocity_tol,
        static_mu_array=static_mu_array,
        kinetic_mu_array=kinetic_mu_array,
    )

    damping_constant = 1.9
    time_step = 1e-4
    snake_sim.dampen(shearable_rod).using(
        AnalyticalLinearDamper,
        damping_constant=damping_constant,
        time_step=time_step,
    )

    total_steps = int(final_time / time_step)
    rendering_fps = 60
    step_skip = int(1.0 / (rendering_fps * time_step))

    # Add call backs
    class ContinuumSnakeCallBack(CallBackBaseClass):
        """
        Call back function for continuum snake
        """

        def __init__(self, step_skip: int, callback_params: dict):
            CallBackBaseClass.__init__(self)
            self.every = step_skip
            self.callback_params = callback_params

        def make_callback(self, system, time, current_step: int):

            if current_step % self.every == 0:

                self.callback_params["time"].append(time)
                self.callback_params["step"].append(current_step)
                self.callback_params["position"].append(
                    system.position_collection.copy()
                )
                self.callback_params["velocity"].append(
                    system.velocity_collection.copy()
                )
                self.callback_params["avg_velocity"].append(
                    system.compute_velocity_center_of_mass()
                )

                self.callback_params["center_of_mass"].append(
                    system.compute_position_center_of_mass()
                )
                self.callback_params["curvature"].append(system.kappa.copy())

                return

    pp_list = defaultdict(list)
    snake_sim.collect_diagnostics(shearable_rod).using(
        ContinuumSnakeCallBack, step_skip=step_skip, callback_params=pp_list
    )

    snake_sim.finalize()

    timestepper = PositionVerlet()
    integrate(timestepper, snake_sim, final_time, total_steps, progress_bar=False)

    if PLOT_FIGURE:
        filename_plot = "continuum_snake_velocity.png"
        plot_snake_velocity(pp_list, period, filename_plot, SAVE_FIGURE)
        plot_curvature(pp_list, shearable_rod.rest_lengths, period, SAVE_FIGURE)

    if SAVE_VIDEO:
        filename_video = "continuum_snake.mp4"
        plot_video(
            pp_list,
            video_name=filename_video,
            fps=rendering_fps,
            xlim=xlim,
            ylim=(-3, 3),
        )


    distance_traveled = ((pp_list["center_of_mass"][-1] - pp_list["center_of_mass"][0])[2]**2 \
                         + (pp_list["center_of_mass"][-1] - pp_list["center_of_mass"][0])[0]**2)**0.5
    if isnan(distance_traveled):
        distance_traveled = -99999999
    print(distance_traveled)
    return distance_traveled


