from wiggling_snake import*

b_coeff = [ 0,8.245666,   10, 43.68705795, 49.09816949,0]
wave_length = 1.47437266

# b_coeff = [ 0,-5039.3269229,    8290.45412652,   2938.70286774,   3291.909462,0]
#
# wave_length = -12766.60197089
distance_traveled = run_snake_pure_wiggle(b_coeff=b_coeff,wave_length=wave_length,n_elements=20,run_time=5,SAVE_VIDEO=True)
boundary = 300
b_coeff = np.array(b_coeff)
fitness = distance_traveled
abs_coeff = np.abs(b_coeff)
penalty = 10
out_of_bound = abs_coeff > boundary
if np.any(out_of_bound):
    for idx in np.where(abs_coeff > boundary)[0]:
        fitness -= abs_coeff[idx] * penalty
if wave_length < 0:
        fitness -= abs(wave_length) * penalty * 30
