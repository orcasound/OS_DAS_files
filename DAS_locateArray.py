import numpy as np
from scipy.optimize import minimize
from sklearn.metrics import r2_score
import pickle
import matplotlib.pyplot as plt
from DAS_classes import DAS_line
from scipy.optimize import approx_fprime
import time
import DAS_classes
#################################################
# This program is for finding the local positions (xyz) of the hydrophone array
# Send a number of signals from a variety of locations and use the
# recorded (or modeled) arrival time differences to find the actual location of the hydrophones


#################################################


def loss(hydrophone_positions):
    global ITER
    distance_mse = 0
    predicted_diffs = fiber.predict_time_diffs(hydrophone_positions, signal_locations)
    time_diff_mse = np.mean((predicted_diffs - observed_time_diffs) ** 2)

    if B!= 0:
        # Calculate MSE of distances between successive hydrophones
        distances = np.linalg.norm(np.diff(hydrophone_positions.reshape(-1, 3), axis=0), axis=1)
        distancesFiber = fiber.specsGapsAlongFiber
        distance_mse = np.mean((distances-distancesFiber) ** 2)  # Assuming the ideal distances are 0

    loss =  fiber.c * A * time_diff_mse + B * distance_mse
    if ITER % 5000 == 0: print(ITER, "loss values", loss, fiber.c * A * time_diff_mse, B * distance_mse)
    ITER += 1
    return loss

def gradient_of_loss(hydrophone_positions):
    return approx_fprime(hydrophone_positions, loss, epsilon=1e-8)  # Adjust epsilon if needed

#################################################
fiberFile = "fiber_files/fiber_N_100_L_1000_c_1485_initialXYZ_0__0__-5__descend_-5__-50__sinWiggle_5__5__shape_straight_x-axis.pickle"
fiberFile = "fiber_files/fiber_N_100_L_1000_c_1485_initialXYZ_0__0__-5__descend_-30__-100__sinWiggle_15__5__shape_L_x-axis__0.5__North.pickle"
fiberFile = "fiber_files/fiber_N_10_L_1000_c_1485_initialXYZ_0__0__-5__descend_-30__-100__sinWiggle_15__5__shape_L_x-axis__0.5__North.pickle"
fiber = pickle.load(open(fiberFile, "rb"))
DEBUG = 0
ITER = 0
# Constants (Adjust as needed)
NUM_HYDROPHONES = fiber.Nhydros
NUM_SIGNAL_SOURCES = 5
A = 1.0  # Weight for time difference MSE
B = 0.0  # Weight for cable configuration MSE  THIS DOES NOT SEEM TO HELP - IT HURTS THE FITTING

assert isinstance(fiber.xyzsAlongFiber, object)
true_hydrophone_positions = fiber.xyzsAlongFiber
        # In the real world these (x,y,z)'s are only known approximately at first

# setup signal locations at the surface to localize the hydrophones
signal_locations = np.random.rand(NUM_SIGNAL_SOURCES, 3)*fiber.lengthMeters
signal_locations[:,2] = 0   # delete third dimension for surface sources only
#df_signal_locations = pd.DataFrame(signal_locations, columns=['x', 'y', 'z'])

observed_time_diffs = fiber.calculateTimeDiffs(true_hydrophone_positions, signal_locations, True)
# Store in dataframes - Replace with your loading logic
#df_observed_diffs = pd.DataFrame(observed_time_diffs)

print("true_hydrophone_positions")
print(true_hydrophone_positions)
# print("df_signal_locations")
# print(df_signal_locations)
# print("df_observed_diffs")
# print(df_observed_diffs)
print("True loss: ", loss(true_hydrophone_positions))
# Initial guess for hydrophone positions (linear along x-axis)
#initial_positions = np.linspace(0, 1, NUM_HYDROPHONES)[:, np.newaxis] * np.array([[1, 0, 0]])
initial_positions = fiber.initial_xyzsAlongFiber
#initial_positions = true_hydrophone_positions * 0.1
print("initial_positions")
print(initial_positions)
DAS_classes.quadPlot(true_hydrophone_positions, initial_positions, "True hydrophone (red) positions vs initial assumptions (blue)")

#minmization method
method = 'BFGS'  # or 'L-BFGS-B'
method = 'L-BFGS-B'
method = 'BFGS'
tol = 1e-6
ITER = 0

print("Starting minimization --------------------")
start_time = time.time()  # Record the start time
result = minimize(loss, initial_positions.flatten(), method=method, jac=gradient_of_loss, tol=tol)
optimized_positions = result.x.reshape(-1, 3)
fiber.optimized_positions = optimized_positions

print("True vs Optimized losses: ",loss(true_hydrophone_positions), loss(optimized_positions))
if DEBUG == 0:
    print("        x         xmod         y         ymod         z         zmod")
    for i in range(len(true_hydrophone_positions)):
        print("{:d} {:10.2f} {:10.2f}   {:10.2f}{:10.4f}   {:10.2f} {:10.2f}".format(i,\
             true_hydrophone_positions[i][0], optimized_positions[i][0], \
             true_hydrophone_positions[i][1], optimized_positions[i][1],\
             true_hydrophone_positions[i][2], optimized_positions[i][2]))

    print("         dx          dy         dz  errors in hydro positions")
    for i in range(len(true_hydrophone_positions)):
       # print("{:d} {:10.2f}  {:10.2f}  {:10.2f}".format(i,true_hydrophone_positions[i][0] - optimized_positions[i][0], true_hydrophone_positions[i][1] - optimized_positions[i][1]))
       print("{:d} {:10.2f} {:10.2f} {:10.4f}".format(i, \
             true_hydrophone_positions[i][0] - optimized_positions[i][0], \
             true_hydrophone_positions[i][1] - optimized_positions[i][1],\
             true_hydrophone_positions[i][2] - optimized_positions[i][2]))

DAS_classes.quadPlot(true_hydrophone_positions, optimized_positions, "True Hydrophone Positions (red) vs Optimized Positions (blue)")

# Get predicted time differences for the optimized solution
print(".........Calculating predicted time diffs using the optimized_postions for plot and R**2")
predicted_diffs_optimized = fiber.predict_time_diffs(optimized_positions, signal_locations)
plt.scatter(observed_time_diffs.flatten(), predicted_diffs_optimized.flatten(), marker='.')
plt.title("Predicted time differences vs observed time differences for {} pairs".format(len(predicted_diffs_optimized.flatten())))
plt.show()

# Calculate R-squared using sklearn
r2 = r2_score(observed_time_diffs.flatten(), predicted_diffs_optimized.flatten())
print("R-squared:", r2)

elapsed_time = time.time() - start_time
print("With {} hydrophones------------Elapsed time {:.2f} seconds using {} sources".format(NUM_HYDROPHONES ,elapsed_time, NUM_SIGNAL_SOURCES))

pickle.dump(fiber, open(fiberFile, "wb"))
print("saved optimized positions in fiber as {}".format(fiberFile))