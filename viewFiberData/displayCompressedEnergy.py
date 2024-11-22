import pickle
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
# from librosa.core.audio import samplerate
# from scipy.signal import correlate
# from scipy.signal import correlation_lags
# import math
import os

####from ValFiles.extractRawTimeSeries import pulserate

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import h5py  # DAS usualy comes in HDF5 format

############# Methods ######################
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

  with h5py.File(filename, 'r') as f:
    f.visititems(get_info)

  return metadata

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
def loadVariables(metadata):
    gaugelength = pulserate = n_channels = n_samples = None
    for key in metadata.keys():
        group = metadata[key]
        attributes = group['attrs']
        # Iterate through attributes and print their names and values
        # metadata
        for key1, value in group['attrs'].items():
            print(f"key: {key} Attribute: {key1}, Value: {value}")
            if key1 == "GaugeLength":
                gaugelength = value  # meters per channel
            if key1 == "OutputDataRate":
                pulserate = int(value)
                # direct access to an attribute example: pulse_rate = energyData.metadata['Acquisition/Raw[0]']['attrs']['OutputDataRate']
            if key1 == "NumberOfLoci":
                n_channels = int(value)
            if key1 == "Count":
                n_samples = int(value)
    return gaugelength, pulserate, n_channels, n_samples
def getHighestPeaks(islice, ary, N_top=5):
    abs_data = np.abs(ary)
    indices = np.argpartition(abs_data, -N_top)[-N_top:]  # Efficiently find top 10 indices
    # Sort the indices by their corresponding absolute values in descending order
    indices = indices[np.argsort(abs_data[indices])][::-1]
    peaks = [[ary[i], i, islice] for i in indices]
    return peaks


def compress_2d_array(large_array, output_shape, channel_range, NstdMax, transpose=False):
    arr = np.array(large_array[:, channel_range[0]:channel_range[1]]) # DAS files are (seconds, channels)
    rows, cols = arr.shape
    out_rows, out_cols = output_shape

    # Calculate compression factors
    row_factor = rows // out_rows
    col_factor = cols // out_cols

    # Trim the array if necessary
    trimmed_arr = arr[:out_rows * row_factor, :out_cols * col_factor]

    # Reshape and sum
    reshaped = trimmed_arr.reshape(out_rows, row_factor, out_cols, col_factor)
    compressed = reshaped.sum(axis=(1, 3))
    if transpose:
        compressed = np.transpose(compressed)   # we are switching to the y-axis in the first position, x-axis in the second
    mean = np.mean(compressed,)
    std = np.std(compressed)
    max = np.max(compressed)

    if max > mean + NstdMax * std:
        compressed = np.where(compressed > mean+NstdMax*std, 1, compressed)
#        compressed = np.where(compressed > mean + Ncut * std, mean + Ncut * std, compressed)
    return compressed

def compress_array(Ary2, N1, M1):
    """
    Compresses Ary2 into Ary1 by summing values from Ary2 and placing the
    sums in the corresponding elements of Ary1.

    Args:
      Ary1: The first array (N1 by M1) where the sums will be stored.
      Ary2: The second array (N2 by M2) to be compressed.

    Returns:
      Ary1: The modified first array with the compressed sums.
    """
            # e.g.      shape is (120000, 2500)  time, distance
    N2, M2 = Ary2.shape   # This is input array DAS array M2 is distance and N2 is seconds
    Ary1 = np.zeros([N1, M1])  # Similarly, N1 is distance, M1 is seconds

    if N2 <= N1 or M2 <= M1:
        print("can't compress array smaller in a dimension into the larger one")
        return None

    # Compress rows
    row_scale = N2 // N1
    temp_array = np.zeros((N1, M2))  # Initialize with zeros
    for i in range(N1):  # sum over time in Ary2  (first index) and store in bin in temp_array
        row_start = i * row_scale
        row_end = min((i + 1) * row_scale, N2)  # Ensure row_end doesn't exceed N2
        temp_array[i, :] = np.sum(Ary2[row_start:row_end, :], axis=0) # axis=0 -> rows

    # Compress columns
    col_scale = M2 // M1
    for j in range(M1):  # sum over distance in Ary2 (second index)
        col_start = j * col_scale
        col_end = min((j + 1) * col_scale, M2)  # Ensure col_end doesn't exceed M2
        Ary1[:, j] = np.sum(temp_array[:, col_start:col_end], axis=1) # axis=1 -> columns

    return Ary1
############# Classes ######################


class selectedDASenergydata:
    def __init__(self, filename, metadata, energy_raw, energy_filtered, channel_start, lowcut, highcut, order):
        self.filename = filename
        self.metadata = metadata
        self.energy_raw = energy_raw
        self.energy_filtered = energy_filtered
        self.channel_start = channel_start
        self.lowcut = lowcut
        self.highcut = highcut
        self.order = order

############# Parameters ###################
print("This python program is: displayCompressedEnergy.py")

filename = "/media/bigbox/3053-CEE8/decimator_2024-10-09_19.09.00_UTC_000951.h5"
filename = "/home/bigbox/PycharmProjects/DAS_Files/OOI_DAS_2024/ValFiles/data/decimator_2024-11-06_23.43.00_UTC_023887.h5"
loadFrom_h5 = False   # False dignifies use pickled binary file of previous selection

# laser pulses per sec from the print_info method

start_time = "2024-10-09T19:09:00.000714+00:00"
v_sound = 1490

n_Xslices = 50 #10    # sum the square of the filtered timeseries**2 into this number of bins (distance axis)

# these parameters determine how many channel's timeseries are extracted from the h5py database
channel_start =  450# 0# 450#1875
channel_stop  = 500# #500 # 2200   # this can't be larger than the second axis of the input array which is n_channels!!!!

deltaTslices = 60  # sum the samples axis into bins of this width in seconds
# Bandpass filter parameters
lowcut = 500
highcut = 960
order = 4  # Filter order
#  plot controls
plotLogs = True  # If True, display log10 of filtered, compressed stress in heatmap
NstdMax = 4 # in heatmap, 'remove' (set energy = 1) any points greater than NstdMax * std above mean
nrows = 3   # multiple rows should spread the time over subtime intervals of the dataset's time axis max

############# Run analysis #################
    ########## Load fiber variables ################
if loadFrom_h5:
    metadata = get_metadata(filename)
    gaugelength, pulserate, n_channels, n_samples = loadVariables(metadata)
else:   # load from pickle files
    pickleFile = "data/energyDataWhid_1.pkl"
    with open(pickleFile, 'rb') as f:
        energyData= pickle.load(f)
    print("loaded pickle file: ")
    gaugelength, pulserate, n_channels, n_samples = loadVariables(energyData.metadata)

    # Initialize some parameters
    ############# Variables ####################
    n_channelsPerXSlice = (channel_stop - channel_start) // n_Xslices
    n_Tslices = deltaTslices * n_samples // pulserate

    ## Endfire sends a pulse along the fiber and the distance between slices / v_sound is deltaT
    ##     multiply this deltaT by the pulserate to get the maximum lag between slices in samples
    maximumLagBetweenSlices = int((n_channelsPerXSlice * gaugelength / v_sound) * pulserate)

    ## Butterworth filter
    nyquist = 0.5 * pulserate
    low = lowcut / nyquist
    high = highcut / nyquist
    b, a = signal.butter(order, [low, high], btype='band')

    ## event and lag lists
    events = []  # Event will be a list of lags linked by close association with adjacent lags
    lags = []  # each item in lags will be [stackIdx, lag_maxCorr, idx_hmaxCorr, maxCorr

if loadFrom_h5:
    with h5py.File(filename, 'r') as f:
        slices_raw = np.zeros([n_Tslices, n_Xslices])  # placeholders for the slices  N.B. [samples, hydro]
        slices_filtered = np.zeros([n_Tslices, n_Xslices])
        channel_stop = min(channel_stop, n_channelsPerXSlice * n_Xslices)
        array_raw = np.zeros([n_samples, channel_stop])  # placeholders   N. B. DON'T load the whole fiber !!
        array_filtered = np.zeros([n_samples, channel_stop])
        energy_raw = np.zeros([n_samples, channel_stop])
        energy_filtered = np.zeros([n_samples, channel_stop])
        # Access a specific dataset
        dataset = f['Acquisition/Raw[0]/RawData']
        # print("dataset shape", dataset.shape)

        for i in range(0, channel_stop - channel_start):
            try:
                array_raw[:, i] = dataset[:, channel_start + i ] # ??????+ channel_start]
                array_filtered[:, i] = signal.filtfilt(b, a, array_raw[:, i])
            except:
                print('i, channel_stop, channel_start, dataset.shape, array_raw.shape, array_filtered.shape \n ',i, channel_stop, channel_start, dataset.shape, array_raw.shape, array_filtered.shape)
        # Now square the data to get energy
        energy_raw = np.square(array_raw)
        energy_filtered = np.square(array_filtered)

        # save slices_filtered and array_raw in pkl file
        print("starting to save pickle file")
        energyData = selectedDASenergydata(filename, metadata, energy_raw, energy_filtered, channel_start, lowcut, highcut, order)
        pickleFile = "data/energyDataWhid_1.pkl"
        # Check if the file exists and remove it if it does
        if os.path.exists(pickleFile):
            os.remove(pickleFile)
        # Write the new content to the file
        with open(pickleFile, 'ab') as f:
            pickle.dump(energyData, f)


# direct access to an attribute:: example: pulse_rate = energyData.metadata['Acquisition/Raw[0]']['attrs']['OutputDataRate']



############# Variables ####################
n_channelsPerXSlice = (channel_stop - channel_start) // n_Xslices
n_Tslices = deltaTslices * n_samples // pulserate

    ## Endfire sends a pulse along the fiber and the distance between slices / v_sound is deltaT
    ##     multiply this deltaT by the pulserate to get the maximum possible real lag between slices in samples
maximumLagBetweenSlices = int((n_channelsPerXSlice * gaugelength / v_sound) * pulserate)

    ## Butterworth filter
nyquist = 0.5 * pulserate
low = lowcut / nyquist
high = highcut / nyquist
b, a = signal.butter(order, [low, high], btype='band')

    ## event and lag lists
events = []  # Event will be a list of lags linked by close association with adjacent lags
lags = []    # each item in lags will be [stackIdx, lag_maxCorr, idx_hmaxCorr, maxCorr

# slice down to n_Xslices x n_Tslices  --  compress data
# sliced_energy_raw = compress_array(energyData.energy_raw, n_Tslices, n_Xslices)
# sliced_energy_filtered = compress_array(energyData.energy_filtered, n_Tslices, n_Xslices)

channel_stop = min(energyData.energy_raw.shape[1], channel_stop) # number of channels in dataset: energyData.energy_raw.shape[1]
if channel_start >= channel_stop:
    print("!!!!!!!!!!!!!!!!!!!!!channel_start > channel_stop", channel_start, channel_stop)
    quit()
sliced_energy_raw = compress_2d_array(energyData.energy_raw, (n_Tslices, n_Xslices), (channel_start, channel_stop), NstdMax, transpose=True)
sliced_energy_filtered = compress_2d_array(energyData.energy_filtered, (n_Tslices, n_Xslices), (channel_start, channel_stop), NstdMax, transpose=True)


print(f"Input array shape: {energyData.energy_raw.shape} [seconds, channels]")
print(f"Output array shape: {sliced_energy_raw.shape}")
plt.hist(np.log10(sliced_energy_filtered.flatten()), bins=100)
plt.title("Histogram of the log10 of the energy values")
plt.show()
##
## plot arrays

y0 = channel_start * gaugelength / 1000  # convert to km
y1 = channel_stop * gaugelength / 1000
y_range = [y0, y1]
print("y range", y_range, y0, y1, "data shape", sliced_energy_raw.shape, "channel_start", channel_start, "channel_stop", channel_stop, f"N channels/hydro {n_channelsPerXSlice}")
if nrows > 1:
    fig, axes = plt.subplots(nrows, 1, figsize=(15, 15))  # Adjust figsize as needed
    fig.suptitle(f"Compressed Energy vs Time: shape {sliced_energy_raw.shape}\channel start {channel_start} -> channel stop {channel_stop}", fontsize=20)
    # Flatten the axes array for easier iteration

    deltaSecs = 60 // nrows
    i1 = 0
    plotwidth = sliced_energy_raw.shape[1] // nrows
    i2 = plotwidth
    tvalues = np.linspace(0, 60, 5)
    for i in range(nrows):
        ax = axes[i]
        x_range = [i*deltaSecs, (i+1)*deltaSecs]  # For example, x goes from 0 to 30
        theSlice = sliced_energy_filtered[:,i1:i2]
        if plotLogs:
            ax.imshow(np.log10(theSlice), extent = [*x_range, *y_range],  aspect='auto', origin='lower')
            fig.suptitle(
                f"Log 10 Compressed Energy vs Time: shape {sliced_energy_raw.shape}\channel start {channel_start} -> channel stop {channel_stop}",
                fontsize=20)
        else:
            ax.imshow(theSlice, extent = [*x_range, *y_range],  aspect='auto', origin='lower')
            fig.suptitle(
                f"Compressed Energy vs Time: shape {sliced_energy_raw.shape}\channel start {channel_start} -> channel stop {channel_stop}",
                fontsize=20)

        i1 = i2 + 1
        i2 = i2 + plotwidth

        # Set custom tick locationsd
        x_ticks = np.linspace(x_range[0], x_range[1], 5)  # 6 ticks from 0 to 100
        y_ticks = np.linspace(y_range[0], y_range[1], 6)
        ax.set_xticks(x_ticks)
        ax.set_yticks(y_ticks)
        ax.set_xticklabels([f'{x:.1f}' for x in x_ticks], fontsize=16)
        ax.set_yticklabels([f'{x:.1f}' for x in y_ticks], fontsize=16)

        #ax.set_xticklabels(np.linspace(0, 30, 5))
    fig.supxlabel('Time (s)', fontsize=20)
    fig.supylabel("Distance in km", fontsize=20)
    # Adjust the layout to make room for the common x-axis label
    plt.tight_layout()
    # Add some extra space at the bottom for the label
    #plt.subplots_adjust(bottom=0.1)
    plt.show()
else:
    x_range = [0, 60]

    fig, ax1 = plt.subplots()

    if plotLogs:
        im = ax1.imshow(np.log10(sliced_energy_filtered), extent=[*x_range, *y_range], aspect='auto', origin='lower')
        plt.title(
            f"Log10 Compressed Energy vs Time: shape {sliced_energy_raw.shape} \n channel start {channel_start} -> channel stop {channel_stop}, band pass {lowcut}->{highcut} \n number of channels per hydrophone {n_channelsPerXSlice}",
            fontsize=10)
    else:
        im = ax1.imshow(sliced_energy_filtered, extent=[*x_range, *y_range], aspect='auto', origin='lower')
        plt.title(
            f"Compressed Energy vs Time: shape {sliced_energy_raw.shape} \n channel start {channel_start} -> channel stop {channel_stop}, band pass {lowcut}->{highcut} \n number of channels per hydrophone {n_channelsPerXSlice}",
            fontsize=10)

    fig.supxlabel('Time (s)', fontsize=12)
    fig.supylabel("Distance from start of fiber in km", fontsize=12)
    # Adjust the layout to make room for the common x-axis label
    #plt.tight_layout()

    # Create a second y-axis
    ax2 = ax1.twinx()
    chan_ticks = np.linspace(channel_start, channel_stop, 5).astype(int)
    print(chan_ticks)
    # Set labels and ticks for the second axis
    ax2.set_ylabel('Channel number')
    ax2.set_yticks(np.linspace(0, 100, 5).astype(int))#chan_ticks)  # Example tick positions
    ax2.set_yticklabels(chan_ticks)  # Example tick labels
    plt.show()



print("All Done")
# plt.hist(sliced_energy_filtered.flatten())
# plt.show()
# # plot the ten highest slice values for the timeseries
# peakAmpsList = []
# N_top = 1
# for i in range(n_Xslices):
# #for slice in slices_filtered:
#     slice = slices_filtered[:, i]
#     # Get the peakAmps and indices of the N_top largest absolute values of the stacked signals
#     peakAmps = getHighestPeaks(i, slice, N_top)
#     peakAmpsList.append(peakAmps)
#
# #peakAmpsList_sorted = sorted(peakAmpsList, key=lambda x: x[1]) # sort by second element, the peak idx putting together peakAmps with similar channel numbers
# # Flatten the list of lists of lists into a single list of lists
# flat_list = [sublist for sublist_of_sublists in peakAmpsList for sublist in sublist_of_sublists]
# # Sort the flattened list by the second element of each sublist
# sorted_list = sorted(flat_list, key=lambda sublist: sublist[1])
#
#
# # Jitter amount
# jitter_amount = 0.1
# x = []
# y = []
# for i in range(len(peakAmpsList)):
#     #plt.scatter(np.array(peakAmpsList[i])[:,1], np.array(peakAmpsList[i])[:,0])
#     for j in range(N_top):
#         x = np.array(peakAmpsList[i])[j,1]
#         y = np.array(peakAmpsList[i])[j,0]
#         plt.text(x+ 50000*np.random.uniform(-jitter_amount, jitter_amount, size=1)[0],y+ np.random.uniform(-jitter_amount, jitter_amount, size=1)[0], str(f"{i}"), ha='center', va='center', fontsize=10)
#         plt.scatter(x, y, s=3)
#
# plt.title(f"10 highest amplitudes for {len(peakAmpsList)} slices\n "
#           f"5,6 6,7, 7,8 are sequential and close in amplitude")
# plt.xlabel('channel')
# plt.ylabel(' highest absolute amplitdes')
# plt.show()
#
# # choose a first slice and then
# #      with adjacent slices, calculate lag about some of the highest amplitude points in the first slice
#
# i_0 = 6
# maxcorr = peakAmpsList[0][0][0]
# idxmaxcorr = peakAmpsList[0][0][1]
# delta_t = pulserate # this is plus/minus a second
#
# slice0 = slices_filtered[idxmaxcorr-delta_t: idxmaxcorr+delta_t, i_0]
# slice1 = slices_filtered[idxmaxcorr-delta_t: idxmaxcorr+delta_t, i_0+1]
#
# corr = correlate(slice0, slice1, mode='full')
# lags = correlation_lags(len(slice0), len(slice1), mode='full') # returns list of lag indicies
#
# peakCorrs = getHighestPeaks(i_0, corr, 20)
#
# pkIdxsCorr = [itm[1] for itm in peakCorrs]
# lags_at_highest_corr = lags[pkIdxsCorr]
# print(f"Maximum lag (endfire) on the two slices is {maximumLagBetweenSlices}")
# print(f"lags   corr Peaks  peak Indices (slices {i_0}-{i_0+1})")
# for i in range(len(lags_at_highest_corr)):
#     print(f"{lags_at_highest_corr[i]} \t {peakCorrs[i][0]:0.1f} \t   {peakCorrs[i][1]}")

# above is list of high peakCorrs and corr lags in number of samples
# each slice has n_channelsPerXSlice (eg 150) stacked

# work through some of the indexes of the peakCorrs for this correlation pair
#    now examining i_0+1 correlated with i_0+2  and later with - steps
#    build linked lists of 'nearby' channels with 'similar' lags

