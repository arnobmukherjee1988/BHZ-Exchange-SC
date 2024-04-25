import numpy as np
import scipy as sp
from scipy.linalg import expm
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

# ### Defining aliases for different mathematical fucntions
atan = np.arctan
sin = np.sin
cos = np.cos
exp = np.exp
sqrt = np.sqrt
ln = np.log
Im = np.imag
PI = np.pi
iota = complex(0.0, 1.0)
conj = np.conj
PauliX = np.array([[0, 1], [1, 0]])
PauliY = np.array([[0, -1j], [1j, 0]])
PauliZ = np.array([[1, 0], [0, -1]])

# ### Defining various generic functions
# function for easy matrix print
def PrintMatrix (matrix):
  # Loop over each row
  for i in range(matrix.shape[0]):
    # Loop over each column in the current row
    for j in range(matrix.shape[1]):
      # Print element at row i, column j
      print(f"{i:<6} {j:<10} [{i}][{j}]: {matrix[i][j]:<4}")
      # print(f"i: {i:<2}, j: {j:<2}, matrix[{i}][{j}]: {matrix[i][j]:<4}")

# function to calculate distance between two 2-D points
def Distance (x1, y1, x2, y2):
  val = np.sqrt ( (x2-x1)**2 + (y2-y1)**2 )
  return val

# Kronecker delta function.
def delta(i, j):
  return 1 if i == j else 0

# Theta function
def ThetaFunc(x):
  if (x>0.0):
      val=1.0
  else:
      val=0.0
  return val

def save_matrices_to_file(datafile, *data):
  if not data:
    raise ValueError("At least one array or matrix must be provided.")

  # Determine if the data is 1D or 2D
  if isinstance(data[0], np.ndarray) and data[0].ndim == 1:
    NumElements = len(data[0])

    with open(datafile, "w") as file:
      # Print headers
      file.write("#\t\t{:<6}\t{:<6}\t".format("i", "Value"))
      file.write("\n")

      # Iterate through array and print values
      for i in range(NumElements):
        line = "\t\t{:<6}\t{:<15.6f}\t".format(i+1, float(data[0][i]))
        line += "\n"
        file.write(line)
  elif isinstance(data[0], np.ndarray) and data[0].ndim == 2:
    NumRow, NumColumn = data[0].shape

    with open(datafile, "w") as file:
      # Print headers
      file.write("#\t\t{:<6}\t{:<6}\t{:<6}\t".format("i", "j", "Site"))

      # Print matrix names based on variable names
      for matrix in data:
        matrix_name = [name for name, obj in inspect.currentframe().f_back.f_locals.items() if obj is matrix][0]
        file.write("{:<15}\t".format(matrix_name))

      file.write("\n")

      # Iterate through matrices and print values
      for i in range(NumRow):
        for j in range(NumColumn):
          line = "\t\t{:<6}\t{:<6}\t{:<6}\t".format(i+1, j+1, (i*NumRow+j)+1)

          # Print values for each matrix
          for k in range(len(data)):
            line += "{:<15.6f}\t".format(data[k][i][j])

          line += "\n"
          file.write(line)
  else:
    raise ValueError("Input must be a 1D or 2D NumPy array.")
