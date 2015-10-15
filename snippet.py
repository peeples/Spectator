import multiprocessing as mp

pool = mp.Pool(processes=10)
pool.map(modify_data, raw_files)

Where modify_data is a function that takes one input variable.

def modify_data(item):
