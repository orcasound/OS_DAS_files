import numpy as np
import pandas as pd
from scipy.optimize import minimize
from sklearn.metrics import r2_score
import pickle
import matplotlib.pyplot as plt
from DAS_classes import DAS_line
from scipy.optimize import approx_fprime
import time
import DAS_classes

"""
Starting from time of arrival differences on all the hydrophones,
     determine the location of the sound source that created these differences
"""
##########################################
def loss(source_location):   #  Note: we are using the fiber's optimized hydrophone positions
    global ITER              #        to calculate the predicted arrival time differences for
                             #        this particular source_location during the search
    predicted_diffs = fiber.predict_time_diffs(fiber.optimized_positions, [source_location])
    time_diff_mse = np.mean((predicted_diffs - observed_time_diffs) ** 2)
    if ITER % 10 == 0: print(ITER, time_diff_mse)
    ITER += 1
    return time_diff_mse

##########################################

fiberfile = "fiber_files/fiber_NH_100_L_1000_c_1485_NS_20_-5__descend_-30__-100__sinWiggle_15__5__shape_L_x-axis__0.5__North.pickle"
fiberfile = "fiber_files/fiber_NH_100_L_1000_c_1485_NS_10_-5__descend_-30__-100__sinWiggle_15__5__shape_L_x-axis__0.5__North.pickle"
fiber = pickle.load(open(fiberfile, "rb"))

S_source = [400, -100, -50]

DEBUG = 0
method = 'L-BFGS-B'
method = 'BFGS'
tol = 1e-9
ITER = 0
###########  Model a source, since we don't have any real data right now

true_hydrophone_positions = fiber.xyzsAlongFiber
signal_locations = [S_source]
observed_time_diffs = fiber.predict_time_diffs(true_hydrophone_positions, signal_locations)
print("predicted time diffs for assumed source", observed_time_diffs)   # These are a model of what will be returned from DAS fiber

##################
## Find location of the 'unknown' signal source which produced these time differences
## use fiber's optimized hydrophone locations

source_location = [10, 10, 10]  # first guess for fitting to 'observed' time differences

x = loss(source_location)
print("initial loss is ", x)

result = minimize(loss, source_location, method=method,  tol=tol)
optimized_source = result.x.reshape(-1, 3)
print("Calculated Source location is ", optimized_source)
print("'Actual' source location is ", S_source)





