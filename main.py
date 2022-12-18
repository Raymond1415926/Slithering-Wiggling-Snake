from wiggling_snake import*

n_elements = 50

#pure wiggling
#[496.07627175, 499.68782902, 499.87980135, 474.23855433, 255.0977935 ]
b_coeff = [ 0,496.07627175, 499.68782902, 499.87980135, 474.23855433,0 ] #in N
wave_length = 255.0977935#in cm
percent_crawling = 100
percent_crawling /= 100 #obtain input for run fuction in digital form instead of percentage
wave_length /= 100 #obtain wave_length in meters

#50% case
#[-47.09704241,490.06000744,499.88922264,499.98137414,82.77147074,50]
# b_coeff = [0,-47.09704241,490.06000744,499.88922264,499.98137414,0] #in N
# wave_length = 82.77147074 #in cm
# percent_crawling = 50 #in %
# percent_crawling /= 100 #obtain input for run fuction in digital form instead of percentage
# wave_length /= 100 #obtain wave_length in meters

#pure twitching
# b_coeff = [0,-8.10107625,325.16388691,499.58505016,301.97029055,0] #in N
# wave_length = 149.56854175 #in cm
# percent_crawling = 0
# percent_crawling /= 100 #obtain input for run fuction in digital form instead of percentage
# wave_length /= 100 #obtain wave_length in meters

#combined case
# b_coeff_and_lambda_percentage = [-66.08819673,367.36506461,493.36972054,488.02881696,78.20414871,24.81136202]
# b_coeff = [-66.08819673,367.36506461,493.36972054,488.02881696] #in N
# wave_length = 78.20414871 #in cm
# percent_crawling = 24.81136202 #in %
# percent_crawling /= 100 #obtain input for run fuction in digital form instead of percentage
# wave_length /= 100 #obtain wave_length in meters

# b_coeff_and_lambda_percentage = [496.07627175, 499.68782902, 499.87980135, 474.23855433, 255.0977935, 100]
# def calc_fitness_combined(*b_coeff_and_lambda_percentage):
#     b_coeff_and_lambda_percentage = np.array(b_coeff_and_lambda_percentage)
#     if (len(b_coeff_and_lambda_percentage) == 6):
#         b_coeff_and_lambda = b_coeff_and_lambda_percentage.flatten()
#     else:
#         assert False, "wrong number of input parameters"
#     b_coeff = b_coeff_and_lambda_percentage[0:4]
#     lambda_m = b_coeff_and_lambda_percentage[4] / 100
#     crawling_percentage = b_coeff_and_lambda_percentage[-1] / 100
#     b_coeffs = np.zeros(6)
#     b_coeffs[1:5] = b_coeff
#     print(b_coeffs)
#     print(lambda_m)
#     print(crawling_percentage)
#     distance_traveled = run_snake(b_coeff=b_coeffs,wave_length=lambda_m,percent_crawling=crawling_percentage,\
#                                   SAVE_VIDEO=True, xlim=(-3,3), n_elements=20,run_time=1)
#
#     boundary = 500
#     b_coeff = np.array(b_coeff)
#     fitness = distance_traveled
#     abs_coeff = np.abs(b_coeff)
#     penalty = 100.0
#     out_of_bound = abs_coeff > boundary
#     if np.any(out_of_bound):
#         for idx in np.where(abs_coeff > boundary)[0]:
#             fitness -= abs_coeff[idx] * penalty
#     if lambda_m < 0:
#         fitness -= abs(lambda_m) * penalty * 30
#     if crawling_percentage < 0 or crawling_percentage > 1:
#         fitness -= abs(crawling_percentage) * penalty * 100
#     return fitness
# print(calc_fitness_combined(*b_coeff_and_lambda_percentage))

distance_traveled = run_snake(b_coeff=b_coeff,wave_length=wave_length,percent_crawling=percent_crawling,\
                                  SAVE_VIDEO=True, xlim=(-3,3), n_elements=50,run_time=5)