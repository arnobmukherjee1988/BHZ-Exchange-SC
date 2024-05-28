import numpy as np
import scipy as sp
from scipy.linalg import expm
import time
from typing import Optional, Tuple
from numba import jit
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

# Defining aliases for different mathematical fucntions
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
PauliX = np.array([[0, 1], [1, 0]], dtype=np.float32)
PauliY = np.array([[0, -1j], [1j, 0]], dtype=np.complex64)
PauliZ = np.array([[1, 0], [0, -1]], dtype=np.float32)

# function for easy matrix print
def PrintMatrix(matrix: np.ndarray):
  for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):
      print(f"{i:<6} {j:<10} [{i}][{j}]: {matrix[i][j]:<4}")

# function to calculate distance between two 2-D points
def Distance(x1: float, y1: float, x2: float, y2: float) -> float:
  return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

# Kronecker delta function
def delta(i: int, j: int) -> int:
  return 1 if i == j else 0
  
# Theta function
def ThetaFunc(x: float) -> float:
  return 1.0 if x > 0.0 else 0.0

# function for save matrices in datafile
def save_matrices_to_file(datafile: str, *data: np.ndarray):
  if not data:
    raise ValueError("At least one array or matrix must be provided.")

  with open(datafile, "w") as file:
    if data[0].ndim == 1:
      NumElements = len(data[0])
      file.write("#\t\t{:<6}\t".format("i"))

      for matrix in data:
        matrix_name = next(name for name, obj in inspect.currentframe().f_back.f_locals.items() if obj is matrix)
        file.write("{:<15}\t".format(matrix_name))
      file.write("\n")

      for i in range(NumElements):
        line = "\t\t{:<6}\t".format(i + 1)
        for k in range(len(data)):
            line += "{:<15.6f}\t".format(data[k][i])
        line += "\n"
        file.write(line)

    elif data[0].ndim == 2:
      NumRow, NumColumn = data[0].shape
      file.write("#\t\t{:<6}\t{:<6}\t{:<6}\t".format("i", "j", "Site"))
      for matrix in data:
        matrix_name = next(name for name, obj in inspect.currentframe().f_back.f_locals.items() if obj is matrix)
        file.write("{:<15}\t".format(matrix_name))
      file.write("\n")

      for i in range(NumRow):
        for j in range(NumColumn):
          line = "\t\t{:<6}\t{:<6}\t{:<6}\t".format(i + 1, j + 1, (i * NumRow + j) + 1)
          for k in range(len(data)):
            line += "{:<15.6f}\t".format(data[k][i][j])
          line += "\n"
          file.write(line)
        file.write("\n")
    else:
        raise ValueError("Input must be a 1D or 2D NumPy array.")

# Lorentzian function
@jit(nopython=True)
def Lorentzian(BroadFac: float, x: float) -> float:
  return (1.0 / PI) * (BroadFac / (BroadFac ** 2 + x ** 2))

# Function to log the runtime to a file with formatted time
def log_runtime(message: str, duration: float):
  def format_time(duration: float) -> str:
    days, remainder = divmod(duration, 86400)
    hours, remainder = divmod(remainder, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{int(days):02d}-{int(hours):02d}:{int(minutes):02d}:{seconds:.6f}"

  formatted_duration = format_time(duration)
  with open("runtimes.log", "a") as file:
    file.write(f"{message}: {formatted_duration}\n")
