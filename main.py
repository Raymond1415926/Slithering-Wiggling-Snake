from wiggling_snake import*

n_elements = 50

#pure wiggling 95 generations
#{'MaxIter': False, 'EqualFunVals': False, 'ConditionCov': False, 'NoEffectCoor': True, 'Stagnation': False, 'TolXUp': False}
#[-150.19929554, 499.68705001, 499.708699,  489.26091413, 285.00816832]
# b_coeff = [ 0, -150.19929554, 499.68705001, 499.708699,  489.26091413, 0 ] #in N
# wave_length = 285.0082#in cm
# percent_crawling = 100
# percent_crawling /= 100 #obtain input for run fuction in digital form instead of percentage
# wave_length /= 100 #obtain wave_length in meters

#50% case
#[-47.09704241,490.06000744,499.88922264,499.98137414,82.77147074,50]
# b_coeff = [0,-47.09704241,490.06000744,499.88922264,499.98137414,0] #in N
# wave_length = 82.77147074 #in cm
# percent_crawling = 50 #in %
# percent_crawling /= 100 #obtain input for run fuction in digital form instead of percentage
# wave_length /= 100 #obtain wave_length in meters

# pure twitching
b_coeff_and_lambda = [219.25880586, 265.95335158, 270.85175921, 344.27440363,  139.08981327]
b_coeff = [0,0,0,0,0,0] #in N
b_coeff[1:5] = b_coeff_and_lambda[0:4]
wave_length = b_coeff_and_lambda[-1] #in cm
percent_crawling = 0
percent_crawling /= 100 #obtain input for run fuction in digital form instead of percentage
wave_length /= 100 #obtain wave_length in meters

#combined case
# b_coeff_and_lambda_percentage = [-66.08819673,367.36506461,493.36972054,488.02881696,78.20414871,24.81136202]
# b_coeff = [-66.08819673,367.36506461,493.36972054,488.02881696] #in N
# wave_length = 78.20414871 #in cm
# percent_crawling = 24.81136202 #in %
# percent_crawling /= 100 #obtain input for run fuction in digital form instead of percentage
# wave_length /= 100 #obtain wave_length in meters

distance_traveled = run_snake(b_coeff=b_coeff,wave_length=wave_length,percent_crawling=percent_crawling,\
                                  SAVE_VIDEO=True, xlim=(0,5), n_elements=20,run_time=1)

