from config import *
from utils import *

# Function to calculate Skyrmion number
def CalculateSkyrmionNumber(Sx, Sy, Sz):
  SkyrRadius = int(SkyrDiameter/2.0)
  NumSkAlongX = int(Lx/SkyrDiameter)
  NumSkAlongY = int(Ly/SkyrDiameter)
  TotalSkyrmions = NumSkAlongX * NumSkAlongY
  chi = 0.0
  for ix in range(Lx):
    for iy in range(Ly):
      chi1 = 0.0
      chi2 = 0.0
      jy = iy
      jx1 = ix + 1
      if jx1 >= Lx:
        jx1 = 0
      jy = iy
      jx2 = ix - 1
      if jx2 < 0:
        jx2 = Lx - 1
      jy1 = iy + 1
      jx = ix
      if jy1 >= Ly:
        jy1 = 0
      jy2 = iy - 1
      jx = ix
      if jy2 < 0:
        jy2 = Ly - 1

      chi1 = (1 / (8.0 * PI)) * ( (Sx[ix][iy] * (Sy[jx1][jy] * Sz[jx][jy1] - Sz[jx1][jy] * Sy[jx][jy1])) +
                                  (Sy[ix][iy] * (Sz[jx1][jy] * Sx[jx][jy1] - Sx[jx1][jy] * Sz[jx][jy1])) +
                                  (Sz[ix][iy] * (Sx[jx1][jy] * Sy[jx][jy1] - Sy[jx1][jy] * Sx[jx][jy1])) )
      chi2 = (1 / (8.0 * PI)) * ( (Sx[ix][iy] * (Sy[jx2][jy] * Sz[jx][jy2] - Sz[jx2][jy] * Sy[jx][jy2])) +
                                  (Sy[ix][iy] * (Sz[jx2][jy] * Sx[jx][jy2] - Sx[jx2][jy] * Sz[jx][jy2])) +
                                  (Sz[ix][iy] * (Sx[jx2][jy] * Sy[jx][jy2] - Sy[jx2][jy] * Sx[jx][jy2])) )
      chi = chi + chi1 + chi2
  return chi
      
# Different way to skyrmion number calculation
def calculate_skyrmion_number(Sx, Sy, Sz):
  SkyrRadius = int(SkyrDiameter/2.0)
  NumSkAlongX = int(Lx/SkyrDiameter)
  NumSkAlongY = int(Ly/SkyrDiameter)
  TotalSkyrmions = NumSkAlongX * NumSkAlongY
  # Initialize variables
  skyrmion_number = 0.0
  dx = 1.0 / (Lx)  # Adjust for the number of points along x
  dy = 1.0 / (Ly)  # Adjust for the number of points along y
  
  for i in range(Lx):
    for j in range(Ly):
      # Calculate derivatives using finite differences
      dSx_dx = (Sx[(i+1)%Lx][j] - Sx[(i-1)%Lx][j]) / (2 * dx)
      dSx_dy = (Sx[i][(j+1)%Ly] - Sx[i][(j-1)%Ly]) / (2 * dy)
      
      dSy_dx = (Sy[(i+1)%Lx][j] - Sy[(i-1)%Lx][j]) / (2 * dx)
      dSy_dy = (Sy[i][(j+1)%Ly] - Sy[i][(j-1)%Ly]) / (2 * dy)
      
      dSz_dx = (Sz[(i+1)%Lx][j] - Sz[(i-1)%Lx][j]) / (2 * dx)
      dSz_dy = (Sz[i][(j+1)%Ly] - Sz[i][(j-1)%Ly]) / (2 * dy)
      # Calculate cross product
      cross_product = ( dSy_dx * dSz_dy - dSz_dx * dSy_dy,
                        dSz_dx * dSx_dy - dSx_dx * dSz_dy,
                        dSx_dx * dSy_dy - dSy_dx * dSx_dy )
      # Calculate dot product
      dot_product = Sx[i][j] * cross_product[0] + Sy[i][j] * cross_product[1] + Sz[i][j] * cross_product[2]
      # Summing up Skyrmion number
      skyrmion_number += dot_product    
  # Normalize Skyrmion number
  skyrmion_number /= (4 * PI * Nsites)
  return skyrmion_number

# Function to create a skyrmion
def NeelSkyrmion(Beta):
  Theta = np.zeros ((Lx,Ly)) ; Phi = np.zeros ((Lx,Ly))
  Sx = np.zeros ((Lx,Ly)) ; Sy = np.zeros ((Lx,Ly)) ; Sz = np.zeros ((Lx,Ly))
  
  SkyrRadius = int(SkyrDiameter/2.0)
  NumSkAlongX = int(Lx/SkyrDiameter)
  NumSkAlongY = int(Ly/SkyrDiameter)
  TotalSkyrmions = NumSkAlongX * NumSkAlongY
  
  for ix in range (Lx):
    for iy in range (Ly):
      for Skyr_ix in range (NumSkAlongX):
        for Skyr_iy in range (NumSkAlongY):
          SkyrCenterX = SkyrDiameter*Skyr_ix  + SkyrRadius
          SkyrCenterY = SkyrDiameter*Skyr_iy  + SkyrRadius
          if (BraviasLattice=="SquareLattice"):
            dis_x = abs(ix*1.0 - SkyrCenterX)
            dis_y = abs(iy*1.0 - SkyrCenterY)
          if(BraviasLattice=="TriangularLattice"):
            dis_x = abs((ix*1.0 - SkyrCenterX) + ((iy*1.0 - SkyrCenterY)*0.5))
            dis_y = abs((sqrt(3.0)/2.0)*(iy*1.0 - SkyrCenterY))
          dis = Distance(dis_x, dis_y, 0, 0)
          
          Theta[ix][iy] += 2.0*atan(SkyrRadius/(dis+offset)) * exp(-Beta*dis) * ThetaFunc(SkyrRadius-dis)
          if (iy>=SkyrCenterY and ix>=SkyrCenterX):
            Phi[ix][iy] += atan((dis_y+offset)/(dis_x+offset)) * ThetaFunc(SkyrRadius-dis)
          elif (iy>=SkyrCenterY and ix<SkyrCenterX):
            Phi[ix][iy] += (PI - atan((abs(dis_y)+offset)/(abs(dis_x)+offset))) * ThetaFunc(SkyrRadius-dis)
          elif (iy<SkyrCenterY and ix<=SkyrCenterX):
            Phi[ix][iy] += (PI + atan((abs(dis_y)+offset)/(abs(dis_x)+offset))) * ThetaFunc(SkyrRadius-dis)
          elif (iy<SkyrCenterY and ix>SkyrCenterX):
            Phi[ix][iy] += (2.0*PI - atan((abs(dis_y)+offset)/(abs(dis_x)+offset))) * ThetaFunc(SkyrRadius-dis)
          
      Sx[ix][iy] = sin(Theta[ix][iy]) * cos (Phi[ix][iy])
      Sy[ix][iy] = sin(Theta[ix][iy]) * sin (Phi[ix][iy])
      Sz[ix][iy] = cos(Theta[ix][iy])
        
  return Sx, Sy, Sz


def BlochSkyrmion(Beta):
  Theta = np.zeros ((Lx,Ly)) ; Phi = np.zeros ((Lx,Ly))
  Sx = np.zeros ((Lx,Ly)) ; Sy = np.zeros ((Lx,Ly)) ; Sz = np.zeros ((Lx,Ly))
  
  SkyrRadius = int(SkyrDiameter/2.0)
  NumSkAlongX = int(Lx/SkyrDiameter)
  NumSkAlongY = int(Ly/SkyrDiameter)
  TotalSkyrmions = NumSkAlongX * NumSkAlongY
  
  for ix in range (Lx):
    for iy in range (Ly):
      for Skyr_ix in range (NumSkAlongX):
        for Skyr_iy in range (NumSkAlongY):
          SkyrCenterX = SkyrDiameter*Skyr_ix  + SkyrRadius
          SkyrCenterY = SkyrDiameter*Skyr_iy  + SkyrRadius
          if (BraviasLattice=="SquareLattice"):
            dis_x = abs(ix*1.0 - SkyrCenterX)
            dis_y = abs(iy*1.0 - SkyrCenterY)
          if(BraviasLattice=="TriangularLattice"):
            dis_x = abs((ix*1.0 - SkyrCenterX) + ((iy*1.0 - SkyrCenterY)*0.5))
            dis_y = abs((sqrt(3.0)/2.0)*(iy*1.0 - SkyrCenterY))
          dis = Distance(dis_x, dis_y, 0, 0)
          
          Theta[ix][iy] += 2.0*atan(SkyrRadius/(dis+offset)) * exp(-Beta*dis) * ThetaFunc(SkyrRadius-dis)
          if (iy>=SkyrCenterY and ix>=SkyrCenterX):
            Phi[ix][iy] += (atan((dis_y+offset)/(dis_x+offset)) + PI/2.0) * ThetaFunc(SkyrRadius-dis)
          elif (iy>=SkyrCenterY and ix<SkyrCenterX):
            Phi[ix][iy] += ((PI - atan((abs(dis_y)+offset)/(abs(dis_x)+offset))) + PI/2.0) * ThetaFunc(SkyrRadius-dis)
          elif (iy<SkyrCenterY and ix<=SkyrCenterX):
            Phi[ix][iy] += ((PI + atan((abs(dis_y)+offset)/(abs(dis_x)+offset))) + PI/2.0) * ThetaFunc(SkyrRadius-dis)
          elif (iy<SkyrCenterY and ix>SkyrCenterX):
            Phi[ix][iy] += ((2.0*PI - atan((abs(dis_y)+offset)/(abs(dis_x)+offset))) + PI/2.0) * ThetaFunc(SkyrRadius-dis)
          
      Sx[ix][iy] = sin(Theta[ix][iy]) * cos (Phi[ix][iy])
      Sy[ix][iy] = sin(Theta[ix][iy]) * sin (Phi[ix][iy])
      Sz[ix][iy] = cos(Theta[ix][iy])

  return Sx, Sy, Sz


def AntiSkyrmion(Beta):
  Theta = np.zeros ((Lx,Ly)) ; Phi = np.zeros ((Lx,Ly))
  Sx = np.zeros ((Lx,Ly)) ; Sy = np.zeros ((Lx,Ly)) ; Sz = np.zeros ((Lx,Ly))
  
  SkyrRadius = int(SkyrDiameter/2.0)
  NumSkAlongX = int(Lx/SkyrDiameter)
  NumSkAlongY = int(Ly/SkyrDiameter)
  TotalSkyrmions = NumSkAlongX * NumSkAlongY
  
  for ix in range (Lx):
    for iy in range (Ly):
      for Skyr_ix in range (NumSkAlongX):
        for Skyr_iy in range (NumSkAlongY):
          SkyrCenterX = SkyrDiameter*Skyr_ix  + SkyrRadius
          SkyrCenterY = SkyrDiameter*Skyr_iy  + SkyrRadius
          if (BraviasLattice=="SquareLattice"):
            dis_x = abs(ix*1.0 - SkyrCenterX)
            dis_y = abs(iy*1.0 - SkyrCenterY)
          if(BraviasLattice=="TriangularLattice"):
            dis_x = abs((ix*1.0 - SkyrCenterX) + ((iy*1.0 - SkyrCenterY)*0.5))
            dis_y = abs((sqrt(3.0)/2.0)*(iy*1.0 - SkyrCenterY))
          dis = Distance(dis_x, dis_y, 0, 0)
          
          Theta[ix][iy] += 2.0*atan(SkyrRadius/(dis+offset)) * exp(-Beta*dis) * ThetaFunc(SkyrRadius-dis)
          if( iy>=SkyrCenterY and ix>=SkyrCenterX ):
            Phi[ix][iy] += (atan((dis_x+offset)/(dis_y+offset)) ) * ThetaFunc(SkyrRadius-dis)
          elif (iy<SkyrCenterY and ix>=SkyrCenterX):
            Phi[ix][iy] += (PI - atan(  (abs(dis_x)+offset)/(abs(dis_y)+offset)  ) ) * ThetaFunc(SkyrRadius-dis)
          elif (iy<=SkyrCenterY and ix<SkyrCenterX):
            Phi[ix][iy] += (PI + atan(  (abs(dis_x)+offset)/(abs(dis_y)+offset)  ) ) * ThetaFunc(SkyrRadius-dis)
          elif (iy>SkyrCenterY and ix<SkyrCenterX):
            Phi[ix][iy] += (2*PI - atan(  (abs(dis_x)+offset)/(abs(dis_y)+offset)  ) ) * ThetaFunc(SkyrRadius-dis)
          
      Sx[ix][iy] = sin(Theta[ix][iy]) * cos (Phi[ix][iy])
      Sy[ix][iy] = sin(Theta[ix][iy]) * sin (Phi[ix][iy])
      Sz[ix][iy] = cos(Theta[ix][iy])
      
  return Sx, Sy, Sz

def SpinSpiral():
  Theta = np.zeros ((Lx,Ly)) ; Phi = np.zeros ((Lx,Ly))
  Sx = np.zeros ((Lx,Ly)) ; Sy = np.zeros ((Lx,Ly)) ; Sz = np.zeros ((Lx,Ly))
  for ix in range (Lx):
    for iy in range (Ly):
      Theta[ix, iy] = (qx * ix) + (qy * iy) ; Phi[ix, iy] = (qx * ix) + (qy * iy)
      Sx[ix][iy] = sin(Theta[ix][iy]) * cos (Phi[ix][iy])
      Sy[ix][iy] = sin(Theta[ix][iy]) * sin (Phi[ix][iy])
      Sz[ix][iy] = cos(Theta[ix][iy])  
  return Sx, Sy, Sz

def PRB_109_L041409():
  Theta = np.zeros ((Lx,Ly)) ; Phi = np.zeros ((Lx,Ly))
  Sx = np.zeros ((Lx,Ly)) ; Sy = np.zeros ((Lx,Ly)) ; Sz = np.zeros ((Lx,Ly))
  for ix in range (Lx):
    for iy in range (Ly):
      Theta[ix, iy] = PI/2.0 ; Phi[ix, iy] = (qx * ix) + (qy * iy)
      Sx[ix][iy] = sin(Theta[ix][iy]) * cos (Phi[ix][iy])
      Sy[ix][iy] = sin(Theta[ix][iy]) * sin (Phi[ix][iy])
      Sz[ix][iy] = cos(Theta[ix][iy])  
  return Sx, Sy, Sz

def PRR_4_013225_SWCB4():
  Sx = np.zeros ((Lx,Ly)) ; Sy = np.zeros ((Lx,Ly)) ; Sz = np.zeros ((Lx,Ly))
  for ix in range (Lx):
    for iy in range (Ly):
      Sx[ix][iy] = Mx * sin(qx*ix)
      Sy[ix][iy] = My * sin(qy*iy)
      Sz[ix][iy] = Mz * (cos(qx*ix) + cos(qx*ix)) + Bz
  return Sx, Sy, Sz

def Ferro():
  Theta = np.zeros ((Lx,Ly)) ; Phi = np.zeros ((Lx,Ly))
  Sx = np.zeros ((Lx,Ly)) ; Sy = np.zeros ((Lx,Ly)) ; Sz = np.zeros ((Lx,Ly))  
  for ix in range (Lx):
    for iy in range (Ly):
      Theta[ix, iy] = 0.0 ; Phi[ix, iy] = 0.0
      Sx[ix][iy] = sin(Theta[ix][iy]) * cos (Phi[ix][iy])
      Sy[ix][iy] = sin(Theta[ix][iy]) * sin (Phi[ix][iy])
      Sz[ix][iy] = cos(Theta[ix][iy])
  return Sx, Sy, Sz

def AntiFerro():
  Theta = np.zeros ((Lx,Ly)) ; Phi = np.zeros ((Lx,Ly))
  Sx = np.zeros ((Lx,Ly)) ; Sy = np.zeros ((Lx,Ly)) ; Sz = np.zeros ((Lx,Ly))
  for ix in range (Lx):
    for iy in range (Ly):
      if ( (ix+iy)%2 == 0 ):
        Theta[ix, iy] = PI/2.0 ; Phi[ix, iy] = PI/2.0
      else:
        Theta[ix, iy] = PI/2.0 ; Phi[ix, iy] = 3.0*PI/2.0
      Sx[ix][iy] = sin(Theta[ix][iy]) * cos (Phi[ix][iy])
      Sy[ix][iy] = sin(Theta[ix][iy]) * sin (Phi[ix][iy])
      Sz[ix][iy] = cos(Theta[ix][iy])
  return Sx, Sy, Sz


def CustomSkyrmion(Beta):
  Theta = np.zeros ((Lx,Ly)) ; Phi = np.zeros ((Lx,Ly))
  Sx = np.zeros ((Lx,Ly)) ; Sy = np.zeros ((Lx,Ly)) ; Sz = np.zeros ((Lx,Ly))
  SkyrRadius = int(SkyrDiameter/2.0)
  NumSkAlongX = int(Lx/SkyrDiameter)
  NumSkAlongY = int(Ly/SkyrDiameter)
  TotalSkyrmions = NumSkAlongX * NumSkAlongY
  
  cx, cy = 2*SkyrDiameter, 2*SkyrDiameter ; area_c = cx*cy
  theta_c = np.zeros(area_c) ; phi_c = np.zeros(area_c)
  ix_o, iy_o = cx / 2, cy / 2
  for iy in range(cy):
    for ix in range(cx):
      i = iy * cx + ix
      r = Distance (ix, ix_o, iy, iy_o)
      # Handle division by zero
      if r != 0:
        phi_c[i] = np.arccos(vorticity * (ix - ix_o) / r) - helicity
        if iy - iy_o < 0.0:
          phi_c[i] = -helicity - np.arccos(vorticity * (ix - ix_o) / r)
        if (ix - ix_o) == 0.0 and (iy - iy_o) == 0.0:
          phi_c[i] = 0.0 + helicity
        num = SkyrRadius * np.exp(Beta * (SkyrRadius - r))
        den = r
        if den != 0:
          theta_c[i] = 2.0 * np.arctan(num / den) + np.arccos(-polarity)
          if (ix == 0 or ix == cx - 1) and (iy == 0 or iy == cy - 1):
            theta_c[i] = 0.0 + np.arccos(polarity)
  for j in range(Ly // cy):
    for i in range(Lx // cx):
      for j2 in range(cy):
        for i1 in range(cx):
          ii = i * cx + i1
          jj = j * cy + j2
          iii = jj * Lx + ii
          Phi[ii,jj] = phi_c[j2 * cx + i1]
          Theta[ii,jj] = theta_c[j2 * cx + i1]
  for ix in range (Lx):
    for iy in range (Ly):     
      Sx[ix][iy] = sin(Theta[ix][iy]) * cos (Phi[ix][iy])
      Sy[ix][iy] = sin(Theta[ix][iy]) * sin (Phi[ix][iy])
      Sz[ix][iy] = cos(Theta[ix][iy])
          
  return Sx, Sy, Sz

# Function to optimize Skyrmion configuration
def get_spin():
  prev_skyr_number = None
  beta_critical = None
  for i_Beta in range(no_Beta):
    beta_value = Beta_min + ((Beta_max - Beta_min) / no_Beta) * i_Beta
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
    skyr_number = calculate_skyrmion_number(Sx, Sy, Sz)
    # print(i_Beta, beta_value, skyr_number)
    if i_Beta > 0 and np.abs(skyr_number) <= np.abs(prev_skyr_number):
      beta_critical = beta_value = Beta_min + ((Beta_max - Beta_min) / no_Beta) * (i_Beta-1)
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
    # Handle the case when beta_critical is not found
    return None