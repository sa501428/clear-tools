
import matplotlib
import matplotlib.pyplot as plt
import sys
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
matplotlib.use('Agg')


color_dict = {'red': ((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)), 'green': ((0.0, 1.0, 1.0), (1.0, 0.0, 0.0)),
              'blue': ((0.0, 1.0, 1.0), (1.0, 0.0, 0.0))}
red_map = LinearSegmentedColormap('Redmap', color_dict)


def plot_map(matrix, colormap, clim1, clim2, score, ax, x_label):
    img_plot = ax.matshow(matrix)
    img_plot.set_cmap(colormap)
    img_plot.set_clim(clim1, clim2)
    ax.set_frame_on(True)
    ax.xaxis.set_ticklabels(np.arange(0, 1), visible=False)
    ax.yaxis.set_ticklabels(np.arange(0, 1), visible=False)
    ax.xaxis.set_tick_params(length=0, labelsize=5, labeltop='off', labelbottom='on')
    ax.yaxis.set_tick_params(length=0, labelsize=5)
    ax.set_title("APA score: " + '{:.2f}'.format(score))
    ax.set_xlabel(x_label)


def plot_from_npy(m_name, ax, p, is_intra=True):
    matrix = np.load(m_name)
    n = np.shape(matrix)[0]
    color_lim1 = 0 # np.min(matrix)/2
    color_lim2 = max(3 * np.mean(matrix[:n // 4, -n // 4:]), 1)
    score = matrix[n // 2, n // 2] / np.mean(matrix[-n // 4:, :n // 4])
    if is_intra:
        xdist0 = 2 ** (p - 1)
        xdist1 = 2 ** p
        x_label_string = str(xdist0) + "MB-" + str(xdist1) + "MB"
    else:
        x_label_string = "INTER"
    plot_map(matrix, red_map, color_lim1, color_lim2, score, ax, x_label_string)


def get_intra_name(stem, k):
    return stem + 'intra_' + str(k) + '_apa.npy'


def get_inter_name(stem):
    return stem + 'inter_apa.npy'


def plot_series(folder_name, out_name):
    figure, axis = plt.subplots(1, 9, figsize=(16, 2))
    for z in range(8):
        f_name = get_intra_name(folder_name, z)
        plot_from_npy(f_name, axis[z], z)
    plot_from_npy(get_inter_name(folder_name), axis[8], 0, False)
    plt.savefig(out_name)
    plt.close()


plot_series(sys.argv[1], sys.argv[2])


