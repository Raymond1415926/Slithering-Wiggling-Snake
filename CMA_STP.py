import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm 
from mpl_toolkits.mplot3d import axes3d
from numpy.random import multivariate_normal
import copy
from matplotlib import animation
from numpy.random import multivariate_normal
from wiggling_snake import*
from tqdm import tqdm
import copy
import sys

    
class CMAES:
    """Naive CMA implementation"""
    
    def __init__(self, initial_mean, sigma, popsize, **kwargs):
        """Please do all the initialization. The reserve space and 
        code for collecting the statistics are already provided."""

        # Things that evolve : centroid, sigma, paths etc.
        self.centroid = np.asarray(initial_mean).copy()
        self.sigma = sigma
        self.pc = np.zeros_like(initial_mean)
        self.ps = np.zeros_like(initial_mean)        
        self.C = np.eye(initial_mean.shape[0])
        self.B = np.eye(self.C.shape[0])
        self.diagD = np.ones(initial_mean.shape[0])
        
        self.stopTolXUP = 1e4 * self.sigma
        self.stopping = {"MaxIter": False, "EqualFunVals": False, "ConditionCov": False, "NoEffectCoor": False, "Stagnation": False, "TolXUp": False}

        # Optimal popsize 
        self.popsize = popsize
        self.mu = popsize // 2
        
        # Update weights later on
        # Constant weight policy
        # self.weights = np.ones((self.mu, )) / self.mu

        # Decreasing weight policy
        self.weights = np.arange(self.mu, 0.0, -1.0)
        self.weights /= np.sum(self.weights)
        
        # Negative, Positive weight policy
        # unscaled_weights = np.arange(1.0 ,  1.0 + popsize)
        # unscaled_weights = np.log(0.5 * (popsize + 1.0) / unscaled_weights)

        # Utility variables
        self.dim = initial_mean.shape[0]

        # Expectation of a normal distribution
        self.chiN = np.sqrt(self.dim) * (1.0 - 0.25 / self.dim + 1.0/(21.0 * self.dim**2))
        self.mueff = 1.0 / np.linalg.norm(self.weights, 2)**2
        self.generations = 0
        
        # Options
 
        # Sigma adaptation
        # cs is short for c_sigma
        self.cs = kwargs.get("cs", (2.0 + self.mueff) / (self.dim + self.mueff + 5.0))
        # ds is short for d_sigma
        self.ds = 1.0 + 2.0 * max(0.0, np.sqrt((self.mueff - 1.0)/ (self.dim + 1.0)) - 1.0) + self.cs
        
        # Covariance adaptation
        self.cc = kwargs.get("cc", (4.0 + self.mueff/self.dim) / (self.dim + 4.0 + 2.0 * self.mueff/self.dim))
        self.ccov = 0.0
        # If implementing the latest version of CMA according to the tutorial, 
        # these parameters can be useful
        self.ccov1 = 2.0 / ((self.dim + 1.3)**2 + self.mueff)
        self.ccovmu = min(1.0 - self.ccov1, 2.0 * (self.mueff - 2.0 + 1.0/self.mueff)/((self.dim + 2.0)**2 + self.mueff))
        self.cond = 0
        
        self.reverse = kwargs.get("reverse", False)
        self.numgens = kwargs.get("generations", 200)

        self.stats_centroids = []
        self.stats_new_centroids = []
        self.stats_covs = []
        self.stats_new_covs = []
        self.stats_offspring = []
        self.stats_offspring_weights = []
        self.stats_ps = []
        self.stats_maxfitness = []
        self.stats_fitness = []
        self.stats_medianfitness = []
        
    def update(self, problem, population):
        """Update the current covariance matrix strategy from the
        *population*.
        
        :param population: A list of individuals from which to update the
                           parameters.
        """
        # -- store current state of the algorithm
        self.stats_centroids.append(copy.deepcopy(self.centroid))
        self.stats_covs.append(copy.deepcopy(self.C))

        
        

        
        
        population.sort(key=lambda ind: problem(*ind), reverse=self.reverse)
        # population.sort(key=lambda ind: problem(ind[0], ind[1]))
        # population.sort(key=problem)
        
        
        
        
        # -- store sorted offspring
        self.stats_offspring.append(copy.deepcopy(population))
        
        old_centroid = self.centroid
        # Note : the following does m <- <x>_w
        # Note : this is equivalent to doing m <- m + sigma * <z>_w
        # as x = m + sigma * z provided the weights sum to 1.0 which it
        # does
        self.centroid = np.dot(self.weights, population[0:self.mu])
        
        # -- store new centroid
        self.stats_new_centroids.append(copy.deepcopy(self.centroid))
        
        c_diff = self.centroid - old_centroid
        
        # Cumulation : update evolution path
        # Equivalent to in-class definition
        self.ps = (1 - self.cs) * self.ps \
             + np.sqrt(self.cs * (2 - self.cs) * self.mueff) / self.sigma \
             * np.dot(self.B, (1. / self.diagD) * np.dot(self.B.T, c_diff))
        
        # -- store new evol path
        self.stats_ps.append(copy.deepcopy(self.ps))
        
        hsig = float((np.linalg.norm(self.ps) / 
                np.sqrt(1. - (1. - self.cs)**(2. * (self.generations + 1.))) / self.chiN
                < (1.4 + 2. / (self.dim + 1.))))
        
        self.pc = (1 - self.cc) * self.pc + hsig \
                  * np.sqrt(self.cc * (2 - self.cc) * self.mueff) / self.sigma \
                  * c_diff
        
        # Update covariance matrix
        artmp = population[0:self.mu] - old_centroid
        self.C = (1 - self.ccov1 - self.ccovmu + (1 - hsig) \
                   * self.ccov1 * self.cc * (2 - self.cc)) * self.C \
                + self.ccov1 * np.outer(self.pc, self.pc) \
                + self.ccovmu * np.dot((self.weights * artmp.T), artmp) \
                / self.sigma**2
        
        # -- store new covs
        self.stats_new_covs.append(copy.deepcopy(self.C))
        
        self.sigma *= np.exp((np.linalg.norm(self.ps) / self.chiN - 1.) \
                                * self.cs / self.ds)
        
        self.diagD, self.B = np.linalg.eigh(self.C)
        indx = np.argsort(self.diagD)
        
        self.cond = self.diagD[indx[-1]]/self.diagD[indx[0]]
        
        self.diagD = self.diagD[indx]**0.5
        self.B = self.B[:, indx]
        self.BD = self.B * self.diagD
        
    def run(self, problem):
        # At the start, clear all stored cache and start a new campaign
        self.reset()
        while any(self.stopping.values()) == False:
            # Sample the population here!
            population = list(multivariate_normal(self.centroid, self.sigma**2 * self.C, self.popsize))
            

            # Pass the population to update, which computes all new parameters 
            self.update(problem, population)
            # print(np.array(population).shape)
            self.generations += 1
            
            self.stats_fitness.append(population[0])
            self.stats_maxfitness.append(problem(*(population[0].tolist())))
            self.stats_medianfitness.append(problem(*(population[self.popsize // 2 - 1].tolist())))
            
            self.stconds()
            print("\n *******************************************")
            print(f"The {self.generations} generation")
            print(self.stats_fitness[-1])
        else:
            return population[0]
        
    
    def stconds(self):
        if (self.generations > self.numgens):
            self.stopping["MaxIter"] = True
        last = int(10 + np.ceil(30*self.dim / self.popsize))
        if (last < len(self.stats_medianfitness) and np.isclose(np.max(self.stats_maxfitness[-last:]) - np.min(self.stats_maxfitness[-last:]), 0)):
            self.stopping["EqualFunVals"] = True
        if (self.cond > 10**4):
            self.stopping["ConditionCov"] = True
        if all(np.isclose(self.centroid, self.centroid + 0.2 * self.sigma * np.diag(self.C))):
            self.stopping["NoEffectCoor"] = True
        stagsize = int(len(self.stats_medianfitness) // 5)
        minstagsize = int(120 + np.ceil(30*self.dim / self.popsize))
        minstagsize = max(stagsize, minstagsize)
        if (minstagsize < len(self.stats_medianfitness) - 1):
            minmedian = self.stats_medianfitness[-minstagsize:]
            minbest = self.stats_maxfitness[-minstagsize:]
            l = int(np.floor(len(minbest) / 3))
            if (len(self.stats_medianfitness) > minstagsize):
                self.stopping["Stagnation"] = (np.median(minbest[:l]) < np.median(minbest[-l:]) and np.median(minmedian[:l]) < np.median(minmedian[-l:]))
        if (self.sigma * np.max(self.diagD) > self.stopTolXUP):
            self.stopping["TolXUp"] = True
            
    def reset(self):
        # Clears everything to rerun the problem
        self.stats_centroids = []
        self.stats_new_centroids = []
        self.stats_covs = []
        self.stats_new_covs = []
        self.stats_offspring = []
        self.stats_offspring_weights = []
        self.stats_ps = []
        self.stats_fitness = []
        self.stats_maxfitness = []
        self.stats_medianfitness = []

def calc_fitness_pure(*b_coeff_and_lambda):
    b_coeff_and_lambda = np.array(b_coeff_and_lambda)
    if (len(b_coeff_and_lambda) == 5):
        b_coeff_and_lambda = b_coeff_and_lambda.flatten()
    b_coeff = b_coeff_and_lambda[0:4]
    lambda_m = b_coeff_and_lambda[4]
    b_coeffs = np.zeros(6)
    b_coeffs[1:5] = b_coeff

    original_stdout = sys.stdout
    sys.stdout = open(os.devnull,"w")
    distance_traveled = run_snake_pure_wiggle(b_coeff=b_coeff,wave_length=lambda_m,n_elements=20,run_time=1)
    sys.stdout = original_stdout

    boundary = 500
    b_coeff = np.array(b_coeff)
    fitness = distance_traveled
    abs_coeff = np.abs(b_coeff)
    penalty = 100.0
    out_of_bound = abs_coeff > boundary
    if np.any(out_of_bound):
        for idx in np.where(abs_coeff > boundary)[0]:
            fitness -= abs_coeff[idx] * penalty
    if lambda_m < 0:
            fitness -= abs(lambda_m) * penalty * 30
    return fitness

def calc_fitness_combined(*b_coeff_and_lambda_percentage):
    b_coeff_and_lambda_percentage = np.array(b_coeff_and_lambda_percentage)
    if (len(b_coeff_and_lambda_percentage) == 6):
        b_coeff_and_lambda = b_coeff_and_lambda_percentage.flatten()
    else:
        assert len(b_coeff_and_lambda_percentage) == 6, "wrong number of input parameters"
    b_coeff = b_coeff_and_lambda_percentage[0:4]
    lambda_m = b_coeff_and_lambda_percentage[4]
    crawling_percentage = b_coeff_and_lambda_percentage[-1]
    b_coeffs = np.zeros(6)
    b_coeffs[1:5] = b_coeff

    original_stdout = sys.stdout
    sys.stdout = open(os.devnull,"w")
    distance_traveled = run_snake(b_coeff=b_coeff,wave_length=lambda_m,percent_crawling=crawling_percentage, n_elements=20,run_time=1)
    sys.stdout = original_stdout

    boundary = 500
    b_coeff = np.array(b_coeff)
    fitness = distance_traveled
    abs_coeff = np.abs(b_coeff)
    penalty = 100.0
    out_of_bound = abs_coeff > boundary
    if np.any(out_of_bound):
        for idx in np.where(abs_coeff > boundary)[0]:
            fitness -= abs_coeff[idx] * penalty
    if lambda_m < 0:
        fitness -= abs(lambda_m) * penalty * 30
    if crawling_percentage < 0 or crawling_percentage > 1:
        fitness -= abs(crawling_percentage) * penalty * 50
    return fitness

initial_mean = np.array([150,150,150,150,1])
sigma = 15
pop_size = 30
snake_optimization = CMAES(initial_mean=initial_mean,sigma=sigma,popsize=pop_size,generations=200, reverse=True)
answer = snake_optimization.run(calc_fitness_pure)
print(answer)
plt.plot(snake_optimization.stats_maxfitness)
print(snake_optimization.stopping)
plt.ylabel("Fitness")
plt.xlabel("Generations")
plt.show()