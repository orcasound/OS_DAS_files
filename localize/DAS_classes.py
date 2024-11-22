import numpy as np
import math
import textwrap
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Setting our seed value so that repeated runs use the same random numnbers
np.random.seed(444)
DEBUG = 0
class DAS_line:
    def __init__(self, lengthMeters, N, fiberShape, c = 1500, locationUncertainty=5):
        self.lengthMeters = lengthMeters
        self.Nhydros = N
        self.originX = fiberShape['initialXYZ'][0]
        self.originY = fiberShape['initialXYZ'][1]
        self.originZ = fiberShape['initialXYZ'][2]
        self.fiberShape = fiberShape
        self.sinShape = fiberShape['sinWiggle']  #[ amplitude of sin shift from straight line, Num of periods along line]
        self.c = c
        self.locationUncertainty = locationUncertainty
           # "specs" refers to optical scattering elements
        specsDistanceAlongFiber = np.random.random(N)*self.lengthMeters   # Model: place scattering centers randomly
        specsDistanceAlongFiber.sort()
        self.specsDistanceAlongFiber = specsDistanceAlongFiber
        self.specsGapsAlongFiber = np.diff(self.specsDistanceAlongFiber)
        self.construct_xyzsOffXaxis(N)   # This is a modeled wiggly fiber to be able to simulate real arrival times
        self.optimizedPositions = []
        self.num_sources = 0   # number of sources used to localize and find optimizedPositions of hydrophones                #
        if DEBUG >= 10:
            print("  hydros  X,      Y,     Z")
            for i in range(N):
                print(i, self.xyzsAlongFiber[i][0], self.xyzsAlongFiber[i][1], self.xyzsAlongFiber[i][2])
    def construct_xyzsOffXaxis(self, N):
        x = self.originX
        z = self.originZ
        dd = 0.01   # step along an axis in meters for integrating along curve to find hydro locs
        l=0
        iSpec = 0   # index to step along fiber
        self.xyzsAlongFiber = []
        self.initial_xyzsAlongFiber = []
        dz = dd*math.tan(self.fiberShape['descend'][0]*math.pi/180.0)  # this is delta z as we are decending
        xL = -99
        while iSpec in range(len(self.specsDistanceAlongFiber)) and l < self.lengthMeters:
            try:
                if self.fiberShape['shape']["straight"] == "x-axis":
                    y = self.originY + self.sinShape[0]*math.sin(2*3.14*self.sinShape[1]*x/self.lengthMeters)
                    dl = math.sqrt((self.sinShape[0]*math.cos(2*3.14*self.sinShape[1]*x/self.lengthMeters)*2*3.14*self.sinShape[1]/self.lengthMeters )**2 + 1) * dd
                    if DEBUG == 1:
                        print(iSpec, l, self.specsDistanceAlongFiber[iSpec], l + dl , x, y)
                    if self.specsDistanceAlongFiber[iSpec] > l and self.specsDistanceAlongFiber[iSpec] <= l + dl:
                        self.xyzsAlongFiber.append([x, y, z])
                        # N.B.  The initial locations are 'near' the actual (but unknown in the real world) locations
                        #       The uncertainty in z is 1/10 of the x and y uncertainties as the bottom is the bottom!
                        self.initial_xyzsAlongFiber.append([x+self.locationUncertainty*(np.random.rand()-0.5), (self.originY+self.locationUncertainty*(np.random.rand()-0.5)), z+self.locationUncertainty*(np.random.rand()-0.5)/10.0])
                        iSpec += 1
                    x += dd
                    l += dl
                    z += dz
                    if z < self.fiberShape['descend'][1]:
                        dz = 0
            except:
                try:
                    if self.fiberShape['shape']['L']:   # the fiber will descend and then turn 90 deg at specified fraction of the length of the fiber
                        if xL == -99:
                            y = self.originY + self.sinShape[0] * math.sin(
                                2 * 3.14 * self.sinShape[1] * l / self.lengthMeters)
                            dl = math.sqrt((self.sinShape[0] * math.cos(
                                2 * 3.14 * self.sinShape[1] * l / self.lengthMeters) * 2 * 3.14 * self.sinShape[
                                                1] / self.lengthMeters) ** 2 + 1) * dd
                            if DEBUG == 1:
                                print(iSpec, l, self.specsDistanceAlongFiber[iSpec], l + dl, x, y)
                            if self.specsDistanceAlongFiber[iSpec] > l and self.specsDistanceAlongFiber[iSpec] <= l + dl:
                                self.xyzsAlongFiber.append([x, y, z])
                                #  set z component to have 1/10 error range as do x and y
                                approxLoc = [x+self.locationUncertainty*(np.random.rand()-0.5), (y+self.locationUncertainty*(np.random.rand()-0.5)), z+self.locationUncertainty*(np.random.rand()-0.5)/10.0]
                                self.initial_xyzsAlongFiber.append(approxLoc)
                                iSpec += 1
                            x += dd
                            l += dl
                            z += dz
                        if l > self.lengthMeters * self.fiberShape['shape']['L'][1]:
                            if xL == -99: xL = x
                            x = xL + self.sinShape[0] * math.sin(
                                2 * 3.14 * self.sinShape[1] * l / self.lengthMeters)
                            dl = math.sqrt((self.sinShape[0] * math.cos(
                                2 * 3.14 * self.sinShape[1] * l / self.lengthMeters) * 2 * 3.14 * self.sinShape[
                                                1] / self.lengthMeters) ** 2 + 1) * dd
                            if self.specsDistanceAlongFiber[iSpec] > l and self.specsDistanceAlongFiber[
                                iSpec] <= l + dl:
                                self.xyzsAlongFiber.append([x, y, z])
                                self.initial_xyzsAlongFiber.append(
                                    [xL + self.locationUncertainty * (np.random.rand() - 0.5), (y + self.locationUncertainty * (np.random.rand() - 0.5)),
                                     z + self.locationUncertainty * (np.random.rand() - 0.5)/10.0])
                                iSpec += 1
                            if self.fiberShape['shape']['L'][2] == 'North':
                                y -= dd
                            else:
                                y += dd
                            l += dl
                            z += dz
                        if z < self.fiberShape['descend'][1]:
                            dz = 0
                except:
                    print("No shape matched")
                    return
        self.initial_xyzsAlongFiber = np.array(self.initial_xyzsAlongFiber)
        self.xyzsAlongFiber = np.array(self.xyzsAlongFiber)    # actual xyz's of hydrophones  -- modeled here
        quadPlot(self.xyzsAlongFiber, self.initial_xyzsAlongFiber, None, "'Actual-red' vs 'Assumed-blue' hydrophone loations")

    def calculateSpeckSourceArrivals(self, Ssource, xyzsAlongFiber, DEBUG=0):
        Sx, Sy, Sz = Ssource
        rangesToSource = np.zeros(len(self.specsDistanceAlongFiber))
        arrivalTimesAtSpecs =  np.zeros(len(self.specsDistanceAlongFiber))
        shoreSuccessiveDeltaTs = np.zeros(len(self.specsDistanceAlongFiber))
        minArrivalTime = 9e9
        x, y, z = zip(*xyzsAlongFiber)
        for i in range(len(self.specsDistanceAlongFiber)):
            rangesToSource[i] = math.sqrt((x[i] - Sx)**2 +(y[i] - Sy)**2 +(z[i] - Sz)**2 )
            arrivalTimesAtSpecs[i] = rangesToSource[i] / self.c
            if i > 0:
                shoreSuccessiveDeltaTs[i-1] = arrivalTimesAtSpecs[i]  - arrivalTimesAtSpecs[i-1]
            if (arrivalTimesAtSpecs[i] < minArrivalTime):
                minArrivalTime = arrivalTimesAtSpecs[i]
        if DEBUG == 1:
            print(" i  fiber dist   deltaTs(ms) ")
            for i in range(len(self.specsDistanceAlongFiber)):
                ##print(i, self.specsAlongFiber[i], self.shoreArrivalTimes[i]*1000, self.shoreSuccessiveDeltaTs[i]*1000)
                print("{} {:12.2f} {:12.4f}".format(i, self.specsDistanceAlongFiber[i], shoreSuccessiveDeltaTs[i]*1000))

        return shoreSuccessiveDeltaTs

    def calculateTimeDiffs(self, hydrophone_positions, signal_locations, DEBUG=False):
        # Reshape into an array of (x, y, z) triads
        positions = hydrophone_positions.reshape(-1, 3)

        predicted_diffs = []
        for signal_loc in signal_locations:
            distances = np.linalg.norm(positions - signal_loc, axis=1)

            diffs = (distances[:, np.newaxis] - distances)
            ###
            # Each row i of diffs represents the pairwise distance differences between the ith signal source and all combinations of hydrophone pairs.
            # For example, diffs[0, j, k] would represent the difference between the distance of the 1st signal source to the jth hydrophone and the distance of the same signal source to the kth hydrophone.
            ###
            # Extract upper triangle elements
            upper_tri_indices = np.triu_indices(self.Nhydros, k=1)
            predicted_diffs.append(diffs[upper_tri_indices])

        return np.array(predicted_diffs)/self.c   # convert to time

    def predict_time_diffs(self, hydrophone_positions, signal_locations):
            # return upper triangle time differences between i,j pairs of hydrophones
        # Reshape into an array of (x, y, z) triplets
        positions = hydrophone_positions.reshape(-1, 3)
        predicted_diffs = []
        for signal_loc in signal_locations:
            distances = np.linalg.norm(positions - signal_loc, axis=1)
            diffs = distances[:, np.newaxis] - distances
            ###
            # Each row i of diffs represents the pairwise distance differences between the ith signal source and all combinations of hydrophone pairs.
            # For example, diffs[0, j, k] would represent the difference between the distance of the 1st signal source to the jth hydrophone and the distance of the same signal source to the kth hydrophone.
            ###
            # Extract upper triangle elements
            upper_tri_indices = np.triu_indices(self.Nhydros, k=1)
            predicted_diffs.append(diffs[upper_tri_indices])
        predicted_diffs = np.array(predicted_diffs) / self.c #convert to time
        return predicted_diffs  # converted to time
def quadPlot(true_hydrophone_positions, calculated_positions, signal_positions, supTitle):
    fig = plt.figure(figsize=(10, 8))  # Create the figure
    fig.suptitle(supTitle, fontsize=14)
    # 2D subplots
    axs = [fig.add_subplot(2, 2, i) for i in range(1, 4)]  # Create the first 3 axes
    # 3D subplot
    ax3d = fig.add_subplot(2, 2, 4, projection='3d')  # Create the 3D axis for the 4th subplot
    # Top view (xy-plane)

    axs[0].plot(calculated_positions[:, 0], calculated_positions[:, 1], c='b', marker='o')
    axs[0].plot(true_hydrophone_positions[:, 0], true_hydrophone_positions[:, 1], c='r', marker='+')
    if isinstance(signal_positions,np.ndarray):
        axs[0].scatter(signal_positions[:, 0], signal_positions[:, 1], c='g', marker='o')
    axs[0].set_title('Top View (Y vs X)')
    axs[0].set_xlabel('X-West', labelpad=10)
    axs[0].set_ylabel('Y-South', labelpad=10)
    # Side view (xz-plane)

    axs[1].plot(calculated_positions[:, 0], calculated_positions[:, 2], c='b', marker='o')
    axs[1].plot(true_hydrophone_positions[:, 0], true_hydrophone_positions[:, 2], c='r', marker='+')
    if isinstance(signal_positions,np.ndarray):
        axs[1].scatter(signal_positions[:, 0], signal_positions[:, 2], c='g', marker='o')
    axs[1].set_title('North View (Z vs X)')
    axs[1].set_xlabel('X-West', labelpad=10)
    axs[1].set_ylabel('Z-Vertical', labelpad=10)
    # Front view (yz-plane)

    axs[2].plot(calculated_positions[:, 1], calculated_positions[:, 2], c='b', marker='o')
    axs[2].plot(true_hydrophone_positions[:, 1], true_hydrophone_positions[:, 2], c='r', marker='+')
    if isinstance(signal_positions,np.ndarray):
        axs[2].scatter(signal_positions[:, 1], signal_positions[:, 2], c='g', marker='o')
    axs[2].set_title('West View (Z vs Y)')
    axs[2].set_xlabel('Y=South', labelpad=10)
    axs[2].set_ylabel('Z-Vertical', labelpad=10)
    # Plot 3D data in the 4th subplot

    ax3d.plot(calculated_positions[:, 0], calculated_positions[:, 1], calculated_positions[:, 2], c='blue', marker='o')
    ax3d.plot(true_hydrophone_positions[:, 0], true_hydrophone_positions[:, 1], true_hydrophone_positions[:, 2], c='red', marker='+')
    if isinstance(signal_positions,np.ndarray):
        ax3d.scatter(signal_positions[:, 0], signal_positions[:, 1], signal_positions[:, 2], c='green', marker='o')
    ax3d.set_title('3D Perspective')
    ax3d.view_init(30, 45)
    # Adjust spacing and labels
    plt.subplots_adjust(hspace=0.4, wspace=0.4)
    ax3d.set_xlabel('X', labelpad=10)
    ax3d.set_ylabel('Y', labelpad=10)
    ax3d.set_zlabel('Z', labelpad=10)
    plt.show()