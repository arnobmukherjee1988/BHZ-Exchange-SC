import numpy as np
import math
import pprint
import inspect
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import os
import matplotlib.cm as cm
from matplotlib import rc
#print(plt.style.available)
plt.style.use("seaborn-v0_8-paper")
rc('axes', edgecolor='k')

nice_fonts = {
        # Use LaTeX to write all text
        "text.usetex": True,
        "font.family": "serif",
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 10,
        "font.size": 10,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        #'text.latex.preamble' : [r'\usepackage{amsmath}'],
        'mathtext.fontset' : 'stix',
        'mathtext.rm' : 'serif'
}
mpl.rcParams.update(nice_fonts)

def visualize_matrices(matrices, titles, multiplot=False, cmap='gnuplot', norm=None, aspect='auto', interpolation='none', vmin=None, vmax=None, extent=None, SaveAs=None):
    """
    Visualizes matrices with optional customization for single or multiplot.

    Parameters:
    - matrices (list): List of matrices to be visualized.
    - titles (list): List of titles corresponding to each matrix in 'matrices'.
    - multiplot (bool, optional): If True, creates a multiplot for each matrix. Default is False.

    Additional Arguments:
    - cmap (str, optional): Colormap for mapping the data values to colors. Default is 'gnuplot'.
    - norm (Normalize, optional): An instance of Normalize to map values to the interval [0, 1]. Default is None.
    - aspect (str, {'auto', 'equal'}, optional): Aspect ratio of the plot. Default is 'auto'.
    - interpolation (str, optional): Interpolation method for displaying the matrix. Default is 'none'.
    - vmin, vmax (float, optional): The range of values to be displayed. Default is the full range of the data.
    - extent (float, optional): The [left, right, bottom, top] extent of the plot. Default is determined by the shape of the input matrix.

    Examples:
    - Single plot: visualize_matrices([Ham_ref], ['Ham_ref'])
    - Multiplot: visualize_matrices([Ham_ref, Ham_Matrix], ['Ham_ref', 'Ham_Matrix'], multiplot=True)

    """
    if multiplot:
        num_plots = len(matrices)
        rows = 1
        cols = num_plots
        fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 6 * rows))

        for i in range(num_plots):
            axes[i].imshow(matrices[i].real, cmap=cmap, norm=norm, aspect=aspect, interpolation=interpolation, vmin=vmin, vmax=vmax, extent=extent)
            axes[i].set_title(titles[i])

            # Create an axes on the right side of the plot and shrink it to match the height
            divider = make_axes_locatable(axes[i])
            cax = divider.append_axes("right", size="5%", pad=0.1)

            # Add colorbar with the modified axes
            cbar = plt.colorbar(axes[i].images[0], cax=cax)

        plt.tight_layout()  # Adjust layout to prevent overlapping
        if SaveAs == 'None':
            plt.show()
        else:
            plt.savefig(SaveAs, bbox_inches='tight')
    else:
        for matrix, title in zip(matrices, titles):
            plt.figure()
            plt.imshow(matrix.real, cmap=cmap, norm=norm, aspect=aspect, interpolation=interpolation, vmin=vmin, vmax=vmax, extent=extent)
            plt.title(title)
            plt.colorbar()
            if SaveAs == 'None':
                plt.show()
            else:
                plt.savefig(SaveAs, bbox_inches='tight')



def plot_data(data, xlim=None, ylim=None, linetype='-', pointtype='o', pointsize=5, color='blue', xlabel='', ylabel='', title='', multiplot=False, SaveAs=None):
    """
    Plots data with optional customization for single or multiplot.

    Parameters:
    - data (list): List of datasets to be plotted. Each dataset is represented as a pair [X_arr, Y_arr].
    - xlim (tuple, optional): Tuple specifying the X-axis limits. Default is None.
    - ylim (tuple, optional): Tuple specifying the Y-axis limits. Default is None.
    - linetype (str, optional): Line style for the plot. Default is '-' (solid line).
    - pointtype (str, optional): Marker style for data points. Default is 'o' (circle).
    - pointsize (int, optional): Size of markers for data points. Default is 5.
    - color (str, optional): Color of the plot. Default is 'blue'.
    - xlabel (str, optional): Label for the X-axis. Default is an empty string.
    - ylabel (str, optional): Label for the Y-axis. Default is an empty string.
    - title (str, optional): Title of the plot. Default is an empty string.
    - multiplot (bool, optional): If True, creates a multiplot for each dataset. Default is False.

    Examples:
    - Single line plot: plot_data([[X1_arr, Y1_arr]], xlim=(0, 10), ylim=(0, 20), linetype='-', pointtype='o', 
                                  pointsize=5, color='blue', xlabel='X-axis', ylabel='Y-axis', title='Plot 1')

    - Multiplot: plot_data([[X1_arr, Y1_arr], [X2_arr, Y2_arr]], xlim=[(0, 10), (0, 15)], ylim=[(0, 20), (0, 25)], 
                            linetype=['-', '--'], pointtype=['o', '^'], pointsize=[5, 7], color=['blue', 'red'], 
                            xlabel=['X-axis', 'X-axis'], ylabel=['Y-axis', 'Y-axis'], title=['Plot 1', 'Plot 2'], multiplot=True)

    - Single scatter plot: plot_data([[X1_arr, Y1_arr]], xlim=(0, 10), ylim=(0, 20), linetype='', pointtype='o', 
                                      pointsize=5, color='blue', xlabel='X-axis', ylabel='Y-axis', title='Scatter Plot')
    
    """
    if multiplot:
        num_plots = len(data)
        rows = 1
        cols = num_plots
        fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 6 * rows))

        for i in range(num_plots):
            x_data, y_data = data[i]
            axes[i].plot(x_data, y_data, linetype, marker=pointtype, markersize=pointsize, color=color, zorder=3)
            axes[i].set_title(title[i] if len(title) > i else '')
            axes[i].set_xlabel(xlabel[i] if len(xlabel) > i else '')
            axes[i].set_ylabel(ylabel[i] if len(ylabel) > i else '')
            axes[i].set_xlim(xlim[i] if xlim and len(xlim) > i else None)
            axes[i].set_ylim(ylim[i] if ylim and len(ylim) > i else None)

        plt.tight_layout()  # Adjust layout to prevent overlapping
        plt.grid()
        if SaveAs == 'None':
            plt.show()
        else:
            plt.savefig(SaveAs, bbox_inches='tight')
    else:
        for dataset in data:
            x_data, y_data = dataset
            plt.figure()
            if linetype:
                plt.plot(x_data, y_data, linetype, marker=pointtype, markersize=pointsize, color=color, zorder=3)
            else:
                plt.scatter(x_data, y_data, marker=pointtype, s=pointsize, color=color, zorder=3)
            plt.title(title)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.xlim(xlim)
            plt.ylim(ylim)
            plt.grid()
            if SaveAs == 'None':
                plt.show()
            else:
                plt.savefig(SaveAs, bbox_inches='tight')



def SpinPlot(datafile=None, Mat1=None, Mat2=None, Mat3=None, arrow_length=1.5, SaveAs='None'):
    """
    Visualize spin configurations using matplotlib.

    Parameters:
    - datafile (str, optional): Path to a data file containing spin configuration. Default is None.
    - sx (ndarray, optional): Matrix representing x-component of spins. Default is None.
    - sy (ndarray, optional): Matrix representing y-component of spins. Default is None.
    - sz (ndarray, optional): Matrix representing z-component of spins. Default is None.
    - separation (float, optional): Distance between neighboring spins in the plot. Default is 4.
    - arrow_length (float, optional): Length of the arrows representing spins. Default is 1.5.

    Notes:
    - Either provide a datafile or all three matrices (sx, sy, sz).
    - If using a datafile, the file should contain columns: x, y, sx, sy, sz.
    - The function creates a 2D plot with arrows representing spin configurations.

    Example usage:
    - Using a datafile: SpinPlot(datafile="your_data_file.txt", separation=4, arrow_length=1.5)
    - Using matrices: SpinPlot(sx=sx_matrix, sy=sy_matrix, sz=sz_matrix, separation=4, arrow_length=1.5)

    """

    if datafile is not None:
        # Load data from file
        x, y, Mat1, Mat2, Mat3 = np.genfromtxt(datafile, usecols=(0, 1, 3, 4, 5), unpack=True)
        num_of_x_data = 1
        for i in range(len(y)):
            if np.abs(y[i + 1] - y[i]) <= 10e-8:
                num_of_x_data = num_of_x_data + 1
            else:
                break
        num_of_y_data = np.int64(np.size(y) / num_of_x_data)
        nx = num_of_x_data
        ny = num_of_y_data

        Mat1 = np.reshape(Mat1, (nx, ny), order='C')
        Mat2 = np.reshape(Mat2, (nx, ny), order='C')
        Mat3 = np.reshape(Mat3, (nx, ny), order='C')

    elif Mat1 is not None and Mat2 is not None and Mat3 is not None:
        # Ensure all matrices are provided
        if Mat1.shape != Mat2.shape or Mat1.shape != Mat3.shape:
            raise ValueError("Input matrices must have the same shape.")
        nx, ny = Mat1.shape

    else:
        raise ValueError("Either provide a datafile or all three matrices (Mat1, Mat2, Mat3).")

    # Create the plot
    fig = plt.figure(figsize=(6, 6))
    mag_sx_sy = np.sqrt(Mat1 ** 2 + Mat2 ** 2)

    ax = plt.axes([0.01, 0.03, 0.85, 0.85])  # [left, bottom, width, height]
    fig.gca().set_aspect('equal', adjustable='box')
    separation = arrow_length * 4.0
    ax.set_xlim(-separation + separation * 0, separation * nx)
    ax.set_ylim(-separation + separation * 0, separation * ny)
    ax.set_xticks([])
    ax.set_yticks([])

    norm = mpl.colors.Normalize(vmin=-1.0, vmax=1.0)
    cmap = plt.cm.gnuplot

    for ix in range(nx):
        for iy in range(ny):
            arrow_width = 0.60 * mag_sx_sy[ix, iy] / np.max(mag_sx_sy)
            arrow_head_width = 3.0 * arrow_width
            arrow_head_length = 1.5 * arrow_head_width
            x_pos = separation * ix - arrow_length * Mat1[ix, iy]
            y_pos = separation * iy - arrow_length * Mat2[ix, iy]
            dx = arrow_length * Mat1[ix, iy]
            dy = arrow_length * Mat2[ix, iy]
            arrow_color = cmap(norm(Mat3[ix, iy]))
            ax.arrow(
                x_pos,
                y_pos,
                dx,
                dy,
                width=arrow_width,
                ec=arrow_color,
                fc=arrow_color,
                head_length=arrow_head_length,
                head_width=arrow_head_width,
            )

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cax = plt.axes([0.87, 0.03, 0.03, 0.85])
    cax.tick_params(labelsize=15)
    fig.colorbar(sm, cax=cax, ticks=[-1.0, -0.5, 0.0, 0.5, 1.0])
    cax.set_yticklabels([r'$-1.0$', r'$-0.5$', r'$0.0$', r'$0.5$', r'$1.0$'])

    if SaveAs == 'None':
        plt.show()
    else:
        plt.savefig(SaveAs, bbox_inches='tight')
        
os.system ('mkdir -p Plot')
# plot the spin distribution    
# Check if the file exists
if os.path.isfile("spin_config.txt"):
    SpinPlot(datafile="spin_config.txt", arrow_length=1.5, SaveAs='Plot/spin.pdf')
else:
    print("Spin config data does not exist.")


# plot the eigen spectrum
if os.path.isfile("eigenvals.txt"):
    X, Y = np.genfromtxt ('eigenvals.txt', usecols = (0,1), unpack = 'True')
    plot_data([[X,Y]], xlim=None, ylim=None, linetype='', pointtype='.', 
                    pointsize=20, color='red', xlabel='Eigen Index', ylabel='Eigen Values',
                    title='Eigen Spectrum', SaveAs='Plot/eigenvalue.pdf')
else:
    print("Eigenvalue data file does not exist.")

# plot the dos
if os.path.isfile ("dos.txt"):
    X, Y = np.genfromtxt ('dos.txt', usecols = (0,1), unpack = 'True')
    plot_data([[X,Y]], xlim=None, ylim=None, linetype='', pointtype='.', 
                    pointsize=20, color='red', xlabel='E', ylabel='g(E)',
                    title='DOS', SaveAs='Plot/dos.pdf')
else:
    print("DOS data file does not exist.")

# plot the ldos at fixed site
if os.path.isfile ("Site_265_ldos_SiteFixed.txt"):
    X, Y = np.genfromtxt ('Site_265_ldos_SiteFixed.txt', usecols = (0,1), unpack = 'True')
    plot_data([[X,Y]], xlim=None, ylim=None, linetype='-', pointtype='', 
                    pointsize=20, color='blue', xlabel='E', ylabel='g(E)',
                    title='LDOS at Fixed Site', SaveAs='Plot/ldos_SiteFixed.pdf')
else:
    print("LDOS data file does not exist.")

# plot the ldos distribution
if os.path.isfile ("E_0.0000_ldos_EnergyFixed.txt"):
    (x,y,z) = np.genfromtxt ('E_0.0000_ldos_EnergyFixed.txt', usecols=(0,1,2), unpack=True)
    xmin = np.min(x) ; xmax = np.max(x)
    ymin = np.min(y) ; ymax = np.max(y)
    zmin = np.min(z) ; zmax = np.max(z)
    num_of_x_data = 1
    for i in range(len(y)):
        if(np.abs(y[i+1]-y[i])<=10e-8):
            num_of_x_data = num_of_x_data+1
        else:
            break
    num_of_y_data = np.int64(np.size(y)/num_of_x_data)
    Z = np.transpose(np.reshape(z,(num_of_x_data,num_of_y_data),order="C"))
    Z1 = np.zeros((num_of_x_data,num_of_y_data))
    for iy in range (num_of_y_data):
        for ix in range (num_of_x_data):
            jy = ix
            jx = num_of_y_data - iy - 1
            Z1[jx,jy] = Z[ix,iy]

    visualize_matrices([Z1], ['LDOS at E=0.0'], extent=[xmin,xmax,ymin,ymax], 
                    vmin = zmin, vmax = zmax, interpolation='None',
                    aspect='equal', SaveAs='Plot/ldos_at_E_0.0.pdf')
else:
    print("LDOS distribution data file does not exist")
