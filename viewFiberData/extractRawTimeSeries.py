import pickle
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
from librosa.core.audio import samplerate
from scipy.signal import correlate
from scipy.signal import correlation_lags
import math
import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import h5py  # DAS usualy comes in HDF5 format

##########################################
############# Methods ######################
def print_info(name, obj):
    if isinstance(obj, h5py.Group):
        print(f"Group: {name}")
        print(f"  Number of members: {len(obj)}")
    elif isinstance(obj, h5py.Dataset):
        print(f"Dataset: {name}")
        print(f"  Shape: {obj.shape}")
        print(f"  Dtype: {obj.dtype}")
    else:
        print(f"Other object: {name}, Type: {type(obj)}")

    # Print attributes for both Groups and Datasets
    if hasattr(obj, 'attrs'):
        print("  Attributes:")
        for key, value in obj.attrs.items():
            print(f"    {key}: {value}")

def get_metadata(filename):
    """
    Extracts metadata from an HDF5 file and returns it as a dictionary.
    Args:
    filename: The path to the HDF5 file.
    Returns:
        metadata: A dictionary containing the metadata.
    """
    metadata = {}

    def get_info(name, obj):
        """
        Callback function for h5py.visititems to collect metadata.
        """
        metadata[name] = {
            'type': type(obj),
            'shape': obj.shape if hasattr(obj, 'shape') else None,
            'dtype': obj.dtype if hasattr(obj, 'dtype') else None,
            'attrs': dict(obj.attrs) if hasattr(obj, 'attrs') else None,
        }
        return

filename = "data/decimator_2024-11-06_23.43.00_UTC_023887.h5"

plotOne = False

#  if plotOne == True then choose a single channel to graph
desiredChannel = 491

#  if plotOne == False  select the channels for istart to istop and plot all these
istart = 480   # these are channel numbers in the fiber dataset
istop  = 496
nrows = 4
ncols = 4

#  Choose start_time and stop_time for the time series of raw fiber stress observations
time_start = 10  # desired times in seconds for plots
time_stop = 20

multiChannelList = np.linspace(istart, istop, istop - istart).astype('int')

x_range = [time_start, time_stop]  # choose x_range and y_range for plots
y_range = [-10, 20]

with h5py.File(filename, 'r') as f:
    f.visititems(print_info)  # print metadata for h5 file
    # metadata = get_metadata(filename)
    # print(metadata.keys())
    for key, value in f['Acquisition'].attrs.items():
        print(f"{key}: {value}")
        if key == 'GaugeLength':
            guagelength = value
        if key == 'PulseRate':
            pulserate = value
        if key == 'MeasurementStartTime':
            starttimestr = value
        if key == 'NumberOfLoci':
            numberofloci = value
    for key, value in f['Acquisition/Raw[0]/RawDataTime'].attrs.items():
        print(f"{key}: {value}")
        if key == 'Count':
            count = value
    dataset = f['Acquisition/Raw[0]/RawData']
 #   pulse_rate = energyData.metadata['Acquisition/Raw[0]']['attrs']['OutputDataRate']


    idxTime1 = int(time_start * pulserate)
    idxTime2 = int(min(time_stop * pulserate, count))

    if plotOne:
        tsTimeseries = dataset[idxTime1:idxTime2, desiredChannel]
        print("tsTimeseries.shape", tsTimeseries.shape)

        fig, ax1 = plt.subplots()
        n_points = len(tsTimeseries)
        timeaxis = np.linspace(time_start, time_stop, n_points)
      #  ax1.set_xticks(np.linspace(idxTime1, idxTime2, num=5).astype(int))
      #  ax1.set_xticklabels(np.linspace(time_start, time_stop, num=5).astype(int))
        ax1.plot(timeaxis, tsTimeseries)
        fig.suptitle(f"Stress vs Seconds from time {time_start}s to {time_stop}s \n Hydrophone is number {desiredChannel} at {(desiredChannel*guagelength)/1000:0.2f} km from beginning of fiber")
        fig.supxlabel('Time (s)', fontsize=12)
        plt.show()

#######################  plot a grid of timeseries
    else:
        # Create the grid of subplots  Choose nrows and ncols for the len(multiChannelList) graphs
        fig, axes = plt.subplots(nrows, ncols, figsize=(12, 12))
        n_plots = len(multiChannelList)
        for i in range(n_plots):
            tsTimeseries = dataset[idxTime1:idxTime2, multiChannelList[i]]
            indexaxis = np.linspace(idxTime1, idxTime2, idxTime2 - idxTime1).astype('int')
            row = i // ncols
            col = i % nrows
            ax = axes[row, col]
            ax.plot(indexaxis, tsTimeseries)
            ax.set_title(f"Channel {multiChannelList[i]}")
        fig.suptitle(f"Stress vs Time: index1={idxTime1} to index2={idxTime2}, Seconds {time_start}s to {time_stop}s ", fontsize=16)
        # Adjust layout
        plt.tight_layout()
        plt.show()



# # Plot each time series
# for i in range(data.shape[1]):
#     row = i // 3  # Calculate subplot row index
#     col = i % 3   # Calculate subplot column index
#     ax = axes[row, col]
#     ax.plot(time_values, data[i])
#     ax.set_title(f'Time Series {i+1}')
#     ax.set_xlabel('Time')
#     ax.set_ylabel('Value')
#
# # Adjust layout
# plt.tight_layout()
# plt.show()
