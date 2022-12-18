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
# {'MaxIter': False, 'EqualFunVals': False, 'ConditionCov': False, 'NoEffectCoor': True, 'Stagnation': False, 'TolXUp': False}
# b_coeff_and_lambda = [135.65041025, 211.66918051, 499.5960613, 498.99765331, 147.94281119]
# b_coeff = [0,0,0,0,0,0] #in N
# b_coeff[1:5] = b_coeff_and_lambda[0:4]
# wave_length = b_coeff_and_lambda[-1] #in cm
# percent_crawling = 0
# percent_crawling /= 100 #obtain input for run fuction in digital form instead of percentage
# wave_length /= 100 #obtain wave_length in meters

#combined case
#{'MaxIter': False, 'EqualFunVals': False, 'ConditionCov': True, 'NoEffectCoor': False, 'Stagnation': False, 'TolXUp': False}
# b_coeff_and_lambda_percentage = [1.17293510e+02,2.86312534e+02,4.30534094e+02,4.85093523e+02,1.49971895e+02,2.74094314e-01]
# b_coeff = [0,0,0,0,0,0] #in N
# b_coeff[1:5] = b_coeff_and_lambda_percentage[0:4]
# wave_length = b_coeff_and_lambda_percentage[4] #in cm
# percent_crawling = b_coeff_and_lambda_percentage[-1] / 100
# percent_crawling /= 100 #obtain input for run fuction in digital form instead of percentage
# wave_length /= 100 #obtain wave_length in meters

#test around
b_coeff_and_lambda_percentage = [174.72321114,254.19679463,411.91380427,292.05419199,141.12135051,3.87812994]
b_coeff = [0,0,0,0,0,0] #in N
b_coeff[1:5] = b_coeff_and_lambda_percentage[0:4]
wave_length = b_coeff_and_lambda_percentage[4] #in cm
percent_crawling = b_coeff_and_lambda_percentage[-1] / 100
percent_crawling /= 100 #obtain input for run fuction in digital form instead of percentage
wave_length /= 100 #obtain wave_length in meters

distance_traveled = run_snake(b_coeff=b_coeff,wave_length=wave_length,percent_crawling=percent_crawling,\
                                  SAVE_VIDEO=False, xlim=(0,5), n_elements=20,run_time=1)

