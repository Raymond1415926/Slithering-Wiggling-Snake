from wiggling_snake import*

b_coeff = [ 0,300,   300, 200, 100,0]
wave_length = 1.47437266

distance_traveled = run_snake(b_coeff=b_coeff,wave_length=wave_length,n_elements=20,run_time=5,SAVE_VIDEO=True)
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