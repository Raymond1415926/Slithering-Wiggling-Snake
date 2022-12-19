from wiggling_snake import*

pure_wiggling = False
pure_twitching = False
combined_case = False
rattle_snake = True


n_elements = 50
runtime = 1
dt = 1e-4


#pure wiggling 95 generations
#{'MaxIter': False, 'EqualFunVals': False, 'ConditionCov': False, 'NoEffectCoor': True, 'Stagnation': False, 'TolXUp': False}
if pure_wiggling:
    b_coeff = [ 0, -150.19929554, 499.68705001, 499.708699,  489.26091413, 0 ] #in N
    wave_length = 285.0082#in cm
    percent_crawling = 100
    percent_crawling /= 100 #obtain input for run fuction in digital form instead of percentage
    wave_length /= 100 #obtain wave_length in meters
    no_friction = False
    period = 1

# pure twitching
# {'MaxIter': False, 'EqualFunVals': False, 'ConditionCov': False, 'NoEffectCoor': True, 'Stagnation': False, 'TolXUp': False}
if pure_twitching:
    b_coeff_and_lambda = [135.65041025, 211.66918051, 499.5960613, 498.99765331, 147.94281119]
    b_coeff = [0,0,0,0,0,0] #in N
    b_coeff[1:5] = b_coeff_and_lambda[0:4]
    wave_length = b_coeff_and_lambda[-1] #in cm
    percent_crawling = 0
    percent_crawling /= 100 #obtain input for run fuction in digital form instead of percentage
    wave_length /= 100 #obtain wave_length in meters
    no_friction = False
    period = 1

#combined case
#{'MaxIter': False, 'EqualFunVals': False, 'ConditionCov': True, 'NoEffectCoor': False, 'Stagnation': False, 'TolXUp': False}
if combined_case:
    b_coeff_and_lambda_percentage = [1.17293510e+02,2.86312534e+02,4.30534094e+02,4.85093523e+02,1.49971895e+02,2.74094314e-01]
    b_coeff = [0,0,0,0,0,0] #in N
    b_coeff[1:5] = b_coeff_and_lambda_percentage[0:4]
    wave_length = b_coeff_and_lambda_percentage[4] #in cm
    percent_crawling = b_coeff_and_lambda_percentage[-1] / 100
    percent_crawling /= 100 #obtain input for run fuction in digital form instead of percentage
    wave_length /= 100 #obtain wave_length in meters
    no_friction = False
    period = 1

#rattle snake
#{'MaxIter': False, 'EqualFunVals': False, 'ConditionCov': True, 'NoEffectCoor': False, 'Stagnation': False, 'TolXUp': False}
if rattle_snake:
    b_coeff_and_lambda_percentage_period = np.array([-4.89530774e+02, -4.92543301e+02, -4.93365239e+02, -2.11807898e+02,
    1.19699537e+02,  1.94942169e-01, 2.87774780e+02])
    b_coeff = [0,0,0,0,0,0] #in N
    b_coeff[1:5] = b_coeff_and_lambda_percentage_period[0:4]
    wave_length = b_coeff_and_lambda_percentage_period[4] #in cm
    percent_crawling = b_coeff_and_lambda_percentage_period[5] #in %
    period = b_coeff_and_lambda_percentage_period[6] #in ms
    period /= 1000
    percent_crawling /= 100 #obtain input for run fuction in digital form instead of percentage
    wave_length /= 100 #obtain wave_length in meters
    no_friction = True

distance_traveled = run_snake(period = period, no_fwd_fric=no_friction,b_coeff=b_coeff,wave_length=wave_length, dt = dt,\
percent_crawling=percent_crawling,SAVE_VIDEO=False, xlim=(0,6), ylim=(-1,1), n_elements=n_elements,run_time=runtime)

