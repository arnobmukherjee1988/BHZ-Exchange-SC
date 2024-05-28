from config import *
from utils import *

# Function to calculate Skyrmion number
import numpy as np

def CalculateSkyrmionNumber(Sx: np.ndarray, Sy: np.ndarray, Sz: np.ndarray) -> float:
  chi: float = 0.0
  # Pre-calculate some constant values
  factor: float = 1.0 / (8.0 * PI)
  # Initialize arrays to store temporary results
  chi1 = np.zeros((Lx, Ly), dtype=float)
  chi2 = np.zeros((Lx, Ly), dtype=float)
  # Loop through the arrays
  for ix in range(Lx):
    for iy in range(Ly):
      jx1 = (ix + 1) % Lx
      jx2 = (ix - 1) % Lx
      jy1 = (iy + 1) % Ly
      jy2 = (iy - 1) % Ly

      chi1[ix, iy] = factor * (
          Sx[ix, iy] * (Sy[jx1, iy] * Sz[ix, jy1] - Sz[jx1, iy] * Sy[ix, jy1]) +
          Sy[ix, iy] * (Sz[jx1, iy] * Sx[ix, jy1] - Sx[jx1, iy] * Sz[ix, jy1]) +
          Sz[ix, iy] * (Sx[jx1, iy] * Sy[ix, jy1] - Sy[jx1, iy] * Sx[ix, jy1])
      )

      chi2[ix, iy] = factor * (
          Sx[ix, iy] * (Sy[jx2, iy] * Sz[ix, jy2] - Sz[jx2, iy] * Sy[ix, jy2]) +
          Sy[ix, iy] * (Sz[jx2, iy] * Sx[ix, jy2] - Sx[jx2, iy] * Sz[ix, jy2]) +
          Sz[ix, iy] * (Sx[jx2, iy] * Sy[ix, jy2] - Sy[jx2, iy] * Sx[ix, jy2])
      )
  # Summing up the results
  chi = np.sum(chi1) + np.sum(chi2)
  return chi
      
# Different way to skyrmion number calculation
def calculate_skyrmion_number(Sx: np.ndarray, Sy: np.ndarray, Sz: np.ndarray) -> float:
  skyrmion_number: float = 0.0
  dx: float = 1.0 / Lx  # Adjust for the number of points along x
  dy: float = 1.0 / Ly  # Adjust for the number of points along y

  for i in range(Lx):
    for j in range(Ly):
      # Calculate derivatives using finite differences
      dSx_dx: float = (Sx[(i + 1) % Lx][j] - Sx[(i - 1) % Lx][j]) / (2 * dx)
      dSx_dy: float = (Sx[i][(j + 1) % Ly] - Sx[i][(j - 1) % Ly]) / (2 * dy)

      dSy_dx: float = (Sy[(i + 1) % Lx][j] - Sy[(i - 1) % Lx][j]) / (2 * dx)
      dSy_dy: float = (Sy[i][(j + 1) % Ly] - Sy[i][(j - 1) % Ly]) / (2 * dy)

      dSz_dx: float = (Sz[(i + 1) % Lx][j] - Sz[(i - 1) % Lx][j]) / (2 * dx)
      dSz_dy: float = (Sz[i][(j + 1) % Ly] - Sz[i][(j - 1) % Ly]) / (2 * dy)
      
      # Calculate cross product
      cross_product: np.ndarray = np.array([dSy_dx * dSz_dy - dSz_dx * dSy_dy,
                                            dSz_dx * dSx_dy - dSx_dx * dSz_dy,
                                            dSx_dx * dSy_dy - dSy_dx * dSx_dy])

      # Calculate dot product
      dot_product: float = Sx[i][j] * cross_product[0] + Sy[i][j] * cross_product[1] + Sz[i][j] * cross_product[2]
      # Summing up Skyrmion number
      skyrmion_number += dot_product

  # Normalize Skyrmion number
  skyrmion_number /= (4 * np.pi * Nsites)
  return skyrmion_number

# Function to create a skyrmion
def NeelSkyrmion(Beta: float) -> tuple:
  Theta = np.zeros((Lx, Ly), dtype=float)
  Phi = np.zeros((Lx, Ly), dtype=float)
  Sx = np.zeros((Lx, Ly), dtype=float)
  Sy = np.zeros((Lx, Ly), dtype=float)
  Sz = np.zeros((Lx, Ly), dtype=float)

  SkyrRadius: int = int(SkyrDiameter / 2.0)
  NumSkAlongX: int = int(Lx / SkyrDiameter)
  NumSkAlongY: int = int(Ly / SkyrDiameter)

  for Skyr_ix in range(NumSkAlongX):
    for Skyr_iy in range(NumSkAlongY):
      SkyrCenterX: float = SkyrDiameter * Skyr_ix + SkyrRadius
      SkyrCenterY: float = SkyrDiameter * Skyr_iy + SkyrRadius
      for ix in range(Lx):
        for iy in range(Ly):
          if BraviasLattice == "SquareLattice":
            dis_x: float = abs(ix * 1.0 - SkyrCenterX)
            dis_y: float = abs(iy * 1.0 - SkyrCenterY)
          elif BraviasLattice == "TriangularLattice":
            dis_x: float = abs((ix * 1.0 - SkyrCenterX) + ((iy * 1.0 - SkyrCenterY) * 0.5))
            dis_y: float = abs((sqrt(3.0) / 2.0) * (iy * 1.0 - SkyrCenterY))
          dis: float = Distance(dis_x, dis_y, 0, 0)

          Theta[ix][iy] += 2.0 * atan(SkyrRadius / (dis + offset)) * exp(-Beta * dis) * ThetaFunc(SkyrRadius - dis)
          if iy >= SkyrCenterY and ix >= SkyrCenterX:
            Phi[ix][iy] += atan((dis_y + offset) / (dis_x + offset)) * ThetaFunc(SkyrRadius - dis)
          elif iy >= SkyrCenterY and ix < SkyrCenterX:
            Phi[ix][iy] += (pi - atan((abs(dis_y) + offset) / (abs(dis_x) + offset))) * ThetaFunc(SkyrRadius - dis)
          elif iy < SkyrCenterY and ix <= SkyrCenterX:
            Phi[ix][iy] += (pi + atan((abs(dis_y) + offset) / (abs(dis_x) + offset))) * ThetaFunc(SkyrRadius - dis)
          elif iy < SkyrCenterY and ix > SkyrCenterX:
            Phi[ix][iy] += (2.0 * pi - atan((abs(dis_y) + offset) / (abs(dis_x) + offset))) * ThetaFunc(SkyrRadius - dis)

  Sx = np.sin(Theta) * np.cos(Phi)
  Sy = np.sin(Theta) * np.sin(Phi)
  Sz = np.cos(Theta)
  
  return Sx, Sy, Sz

def BlochSkyrmion(Beta: float) -> tuple:
  Theta = np.zeros((Lx, Ly), dtype=float)
  Phi = np.zeros((Lx, Ly), dtype=float)
  Sx = np.zeros((Lx, Ly), dtype=float)
  Sy = np.zeros((Lx, Ly), dtype=float)
  Sz = np.zeros((Lx, Ly), dtype=float)

  SkyrRadius: int = int(SkyrDiameter / 2.0)
  NumSkAlongX: int = int(Lx / SkyrDiameter)
  NumSkAlongY: int = int(Ly / SkyrDiameter)

  for Skyr_ix in range(NumSkAlongX):
    for Skyr_iy in range(NumSkAlongY):
      SkyrCenterX: float = SkyrDiameter * Skyr_ix + SkyrRadius
      SkyrCenterY: float = SkyrDiameter * Skyr_iy + SkyrRadius
      for ix in range(Lx):
        for iy in range(Ly):
          if BraviasLattice == "SquareLattice":
            dis_x = abs(ix * 1.0 - SkyrCenterX)
            dis_y = abs(iy * 1.0 - SkyrCenterY)
          elif BraviasLattice == "TriangularLattice":
            dis_x = abs((ix * 1.0 - SkyrCenterX) + ((iy * 1.0 - SkyrCenterY) * 0.5))
            dis_y = abs((sqrt(3.0) / 2.0) * (iy * 1.0 - SkyrCenterY))
          dis = Distance(dis_x, dis_y, 0, 0)
          Theta[ix][iy] += 2.0 * atan(SkyrRadius / (dis + offset)) * exp(-Beta * dis) * ThetaFunc(SkyrRadius - dis)
          if iy >= SkyrCenterY and ix >= SkyrCenterX:
            Phi[ix][iy] += (atan((dis_y + offset) / (dis_x + offset)) + pi / 2.0) * ThetaFunc(SkyrRadius - dis)
          elif iy >= SkyrCenterY and ix < SkyrCenterX:
            Phi[ix][iy] += ((pi - atan((abs(dis_y) + offset) / (abs(dis_x) + offset))) + pi / 2.0) * ThetaFunc(SkyrRadius - dis)
          elif iy < SkyrCenterY and ix <= SkyrCenterX:
            Phi[ix][iy] += ((pi + atan((abs(dis_y) + offset) / (abs(dis_x) + offset))) + pi / 2.0) * ThetaFunc(SkyrRadius - dis)
          elif iy < SkyrCenterY and ix > SkyrCenterX:
            Phi[ix][iy] += ((2.0 * pi - atan((abs(dis_y) + offset) / (abs(dis_x) + offset))) + pi / 2.0) * ThetaFunc(SkyrRadius - dis)

  Sx = np.sin(Theta) * np.cos(Phi)
  Sy = np.sin(Theta) * np.sin(Phi)
  Sz = np.cos(Theta)

  return Sx, Sy, Sz


def AntiSkyrmion(Beta: float) -> tuple:
  Theta = np.zeros((Lx, Ly), dtype=float)
  Phi = np.zeros((Lx, Ly), dtype=float)
  Sx = np.zeros((Lx, Ly), dtype=float)
  Sy = np.zeros((Lx, Ly), dtype=float)
  Sz = np.zeros((Lx, Ly), dtype=float)

  SkyrRadius: int = int(SkyrDiameter / 2.0)
  NumSkAlongX: int = int(Lx / SkyrDiameter)
  NumSkAlongY: int = int(Ly / SkyrDiameter)

  for Skyr_ix in range(NumSkAlongX):
    for Skyr_iy in range(NumSkAlongY):
      SkyrCenterX: float = SkyrDiameter * Skyr_ix + SkyrRadius
      SkyrCenterY: float = SkyrDiameter * Skyr_iy + SkyrRadius
      for ix in range(Lx):
        for iy in range(Ly):
          if BraviasLattice == "SquareLattice":
            dis_x: float = abs(ix * 1.0 - SkyrCenterX)
            dis_y: float = abs(iy * 1.0 - SkyrCenterY)
          elif BraviasLattice == "TriangularLattice":
            dis_x: float = abs((ix * 1.0 - SkyrCenterX) + ((iy * 1.0 - SkyrCenterY) * 0.5))
            dis_y: float = abs((sqrt(3.0) / 2.0) * (iy * 1.0 - SkyrCenterY))
          dis: float = Distance(dis_x, dis_y, 0, 0)
          Theta[ix, iy] += 2.0 * atan(SkyrRadius / (dis + offset)) * exp(-Beta * dis) * ThetaFunc(SkyrRadius - dis)
          if iy >= SkyrCenterY and ix >= SkyrCenterX:
              Phi[ix, iy] += atan((dis_x + offset) / (dis_y + offset)) * ThetaFunc(SkyrRadius - dis)
          elif iy < SkyrCenterY and ix >= SkyrCenterX:
              Phi[ix, iy] += (pi - atan((abs(dis_x) + offset) / (abs(dis_y) + offset))) * ThetaFunc(SkyrRadius - dis)
          elif iy <= SkyrCenterY and ix < SkyrCenterX:
              Phi[ix, iy] += (pi + atan((abs(dis_x) + offset) / (abs(dis_y) + offset))) * ThetaFunc(SkyrRadius - dis)
          elif iy > SkyrCenterY and ix < SkyrCenterX:
              Phi[ix, iy] += (2 * pi - atan((abs(dis_x) + offset) / (abs(dis_y) + offset))) * ThetaFunc(SkyrRadius - dis)

  Sx = np.sin(Theta) * np.cos(Phi)
  Sy = np.sin(Theta) * np.sin(Phi)
  Sz = np.cos(Theta)

  return Sx, Sy, Sz

def SpinSpiral() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
  ix_array, iy_array = np.meshgrid(np.arange(Lx), np.arange(Ly), indexing='ij')
  Theta = qx * ix_array + qy * iy_array
  Phi = Theta  # Since Phi has the same calculation as Theta

  Sx = np.sin(Theta) * np.cos(Phi)
  Sy = np.sin(Theta) * np.sin(Phi)
  Sz = np.cos(Theta)

  return Sx, Sy, Sz

def PRB_109_L041409() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
  ix_array, iy_array = np.meshgrid(np.arange(Lx), np.arange(Ly), indexing='ij')
  Theta = np.full((Lx, Ly), PI / 2.0)
  Phi = qx * ix_array + qy * iy_array

  Sx = np.sin(Theta) * np.cos(Phi)
  Sy = np.sin(Theta) * np.sin(Phi)
  Sz = np.cos(Theta)

  return Sx, Sy, Sz

def PRR_4_013225_SWCB4() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
  ix_array, iy_array = np.meshgrid(np.arange(Lx), np.arange(Ly), indexing='ij')

  Sx = Mx * np.sin(qx * ix_array)
  Sy = My * np.sin(qy * iy_array)
  Sz = Mz * (np.cos(qx * ix_array) + np.cos(qy * iy_array)) + Bz

  return Sx, Sy, Sz

def Ferro() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
  Theta = np.zeros((Lx, Ly), dtype=float)
  Phi = np.zeros((Lx, Ly), dtype=float)

  Sx = np.sin(Theta) * np.cos(Phi)
  Sy = np.sin(Theta) * np.sin(Phi)
  Sz = np.cos(Theta)

  return Sx, Sy, Sz

def AntiFerro() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
  ix_array, iy_array = np.meshgrid(np.arange(Lx), np.arange(Ly), indexing='ij')
  Theta = np.full((Lx, Ly), pi / 2.0)
  Phi = np.where((ix_array + iy_array) % 2 == 0, pi / 2.0, 3.0 * pi / 2.0)

  Sx = np.sin(Theta) * np.cos(Phi)
  Sy = np.sin(Theta) * np.sin(Phi)
  Sz = np.cos(Theta)

  return Sx, Sy, Sz

def CustomSkyrmion(Beta: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
  Theta: np.ndarray = np.zeros((Lx, Ly))
  Phi: np.ndarray = np.zeros((Lx, Ly))
  Sx: np.ndarray = np.zeros((Lx, Ly))
  Sy: np.ndarray = np.zeros((Lx, Ly))
  Sz: np.ndarray = np.zeros((Lx, Ly))
  
  SkyrRadius: int = int(SkyrDiameter / 2.0)
  
  cx: int = int(2 * SkyrDiameter)
  cy: int = int(2 * SkyrDiameter)
  ix_o: float = cx / 2
  iy_o: float = cy / 2
  
  x_grid, y_grid = np.meshgrid(np.arange(cx), np.arange(cy))
  r_grid = Distance(x_grid, ix_o, y_grid, iy_o)
  
  # Handle division by zero for r_grid
  with np.errstate(divide='ignore', invalid='ignore'):
    phi_c = np.arccos(vorticity * (x_grid - ix_o) / r_grid) - helicity
    mask = y_grid - iy_o < 0.0
    phi_c[mask] = -helicity - np.arccos(vorticity * (x_grid[mask] - ix_o) / r_grid[mask])
    phi_c[(x_grid - ix_o) == 0.0 & (y_grid - iy_o) == 0.0] = helicity
    
    num = SkyrRadius * exp(Beta * (SkyrRadius - r_grid))
    theta_c = 2.0 * arctan(num / r_grid) + arccos(-polarity)
    
    corner_mask = ((x_grid == 0) | (x_grid == cx - 1)) & ((y_grid == 0) | (y_grid == cy - 1))
    theta_c[corner_mask] = arccos(polarity)
  
  for j in range(Ly // cy):
    for i in range(Lx // cx):
      Phi[i * cx:(i + 1) * cx, j * cy:(j + 1) * cy] = phi_c
      Theta[i * cx:(i + 1) * cx, j * cy:(j + 1) * cy] = theta_c

  Sx = sin(Theta) * cos(Phi)
  Sy = sin(Theta) * sin(Phi)
  Sz = cos(Theta)
  
  return Sx, Sy, Sz

# Function to optimize Skyrmion configuration
def get_spin() -> (np.ndarray, np.ndarray, np.ndarray, float):
  prev_skyr_number: float = None
  beta_critical: float = None
  for i_Beta in range(no_Beta):
    beta_value: float = Beta_min + ((Beta_max - Beta_min) / no_Beta) * i_Beta
    if SpinConfig_Type == "NeelSkyrmion":
      Sx, Sy, Sz = NeelSkyrmion(beta_value)
    elif SpinConfig_Type == "BlochSkyrmion":
      Sx, Sy, Sz = BlochSkyrmion(beta_value)
    elif SpinConfig_Type == "AntiSkyrmion":
      Sx, Sy, Sz = AntiSkyrmion(beta_value)
    elif SpinConfig_Type == "CustomSkyrmion":
      Sx, Sy, Sz = CustomSkyrmion(beta_value)
    elif SpinConfig_Type == "SpinSpiral":
      Sx, Sy, Sz = SpinSpiral()
    elif SpinConfig_Type == "Ferro":
      Sx, Sy, Sz = Ferro()
    elif SpinConfig_Type == "AntiFerro":
      Sx, Sy, Sz = AntiFerro()
    elif SpinConfig_Type == 'PRR_4_013225_SWCB4':
      Sx, Sy, Sz = PRR_4_013225_SWCB4()
    elif SpinConfig_Type == 'PRB_109_L041409':
      Sx, Sy, Sz = PRB_109_L041409()

    skyr_number: float = calculate_skyrmion_number(Sx, Sy, Sz)

    if i_Beta > 0 and np.abs(skyr_number) <= np.abs(prev_skyr_number):
      beta_critical = Beta_min + ((Beta_max - Beta_min) / no_Beta) * (i_Beta - 1)
      break
    prev_skyr_number = skyr_number
  
  if beta_critical is not None:
    if SpinConfig_Type == "NeelSkyrmion":
      Sx, Sy, Sz = NeelSkyrmion(beta_critical)
    elif SpinConfig_Type == "BlochSkyrmion":
      Sx, Sy, Sz = BlochSkyrmion(beta_critical)
    elif SpinConfig_Type == "AntiSkyrmion":
      Sx, Sy, Sz = AntiSkyrmion(beta_critical)
    elif SpinConfig_Type == "CustomSkyrmion":
      Sx, Sy, Sz = CustomSkyrmion(beta_critical)
    elif SpinConfig_Type == "SpinSpiral":
      Sx, Sy, Sz = SpinSpiral()
    elif SpinConfig_Type == "Ferro":
      Sx, Sy, Sz = Ferro()
    elif SpinConfig_Type == "AntiFerro":
      Sx, Sy, Sz = AntiFerro()
    elif SpinConfig_Type == 'PRR_4_013225_SWCB4':
      Sx, Sy, Sz = PRR_4_013225_SWCB4()
    elif SpinConfig_Type == 'PRB_109_L041409':
      Sx, Sy, Sz = PRB_109_L041409()
    return Sx, Sy, Sz, beta_critical
  else:
    return None

