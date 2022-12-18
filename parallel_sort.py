import multiprocessing
import psutil

def parallel_sort(population, problem, reverse=False, cores="physical", sort_problem=False):
    """
    Perform parallel sorting on a given population.
    
    Parameters:
    population (list): A list of elements to be sorted.
    problem (function): A function that takes an element from the population as input and returns a value that will be used to determine the sort order.
    reverse (bool, optional): A boolean that determines whether the list should be sorted in ascending or descending order. If set to True, the list will be sorted in descending order. Default is False.
    cores (str, optional): A string that determines the number of cores to use for parallel sorting. Valid values are "physical" and "logical". Default is "physical".
    sort_problem (bool, optional): A boolean that determines whether the problem values should also be returned along with the sorted population. If set to True, a tuple containing the sorted problem values and the sorted population is returned. If set to False (default), only the sorted population is returned.
    
    Returns:
    tuple or list: If sort_problem is set to True, returns a tuple containing the sorted problem values and the sorted population. If set to False (default), only the sorted population is returned.
    """
    
    if cores == "physical":
        num_cores = psutil.cpu_count(logical=False)
    elif cores == "logical":
        num_cores = psutil.cpu_count(logical=True)
    else:
        raise ValueError("Invalid value for cores. Valid values are 'physical' and 'logical'.")

    with multiprocessing.Pool(processes=num_cores) as pool:
        # Use the starmap method if the problem function takes more than one argument
        results = pool.map(problem, population)

        # Sort the population list directly based on the corresponding element in the results list
        sorted_population = sorted(population, key=lambda x: results[population.index(x)], reverse=reverse)

        # Close the pool and wait for all worker processes to terminate
        pool.close()
        pool.join()

    if sort_problem:
        # Extract the corresponding element in the results list for each element in the sorted population list
        sorted_output = [results[population.index(x)] for x in sorted_population]
        return sorted_output, sorted_population
    else:
        return sorted_population