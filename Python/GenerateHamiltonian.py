from config import *
from utils import *

# Function to get nearest neighbor arrays
def neighbor_array_square_lattice():
  J1_arr = np.empty(Lx * Ly, dtype=int)
  J2_arr = np.empty(Lx * Ly, dtype=int)
  for iy in range(Ly):
    for ix in range(Lx):
      i = (iy * Lx) + ix
      if BoundaryCondition == 'PBCx + PBCy':
        jx = ix + 1 if ix + 1 < Lx else 0
        jy = iy + 1 if iy + 1 < Ly else 0
        J1_arr[i] = (iy * Lx) + jx
        J2_arr[i] = (jy * Lx) + ix
      elif BoundaryCondition == 'OBCx + PBCy':
        jx = ix + 1 if ix + 1 < Lx else 0
        jy = iy + 1 if iy + 1 < Ly else 0
        J1_arr[i] = (iy * Lx) + jx
        J2_arr[i] = (jy * Lx) + ix
        if ix == Lx - 1:
          J1_arr[i] = -i
      elif BoundaryCondition == 'PBCx + OBCy':
        jx = ix + 1 if ix + 1 < Lx else 0
        jy = iy + 1 if iy + 1 < Ly else 0
        J1_arr[i] = (iy * Lx) + jx
        J2_arr[i] = (jy * Lx) + ix
        if iy == Ly - 1:
          J2_arr[i] = -i
      elif BoundaryCondition == 'OBCx + OBCy':
        jx = ix + 1 if ix + 1 < Lx else 0
        jy = iy + 1 if iy + 1 < Ly else 0
        J1_arr[i] = (iy * Lx) + jx
        J2_arr[i] = (jy * Lx) + ix
        if ix == Lx - 1:
          J1_arr[i] = -i
        if iy == Ly - 1:
          J2_arr[i] = -i
      else:
        raise ValueError(BC_Error)
  return J1_arr, J2_arr
  
  # Function to get nearest neighbor arrays
def neighbor_array_triangular_lattice():
  J1_arr = np.empty(Lx * Ly, dtype=int)
  J2_arr = np.empty(Lx * Ly, dtype=int)
  J3_arr = np.empty(Lx * Ly, dtype=int)
  for iy in range(Ly):
    for ix in range(Lx):
      i = (iy * Lx) + ix
      if BoundaryCondition == 'PBCx + PBCy':
        jx = ix + 1 if ix + 1 < Lx else 0
        jy = iy + 1 if iy + 1 < Ly else 0
        J1_arr[i] = (iy * Lx) + jx
        J2_arr[i] = (jy * Lx) + ix
        J3_arr[i] = (jy * Lx) + jx
      elif BoundaryCondition == 'OBCx + PBCy':
        jx = ix + 1 if ix + 1 < Lx else 0
        jy = iy + 1 if iy + 1 < Ly else 0
        J1_arr[i] = (iy * Lx) + jx
        J2_arr[i] = (jy * Lx) + ix
        J3_arr[i] = (jy * Lx) + ix
        if ix == Lx - 1:
          J1_arr[i] = -i
      elif BoundaryCondition == 'PBCx + OBCy':
        jx = ix + 1 if ix + 1 < Lx else 0
        jy = iy + 1 if iy + 1 < Ly else 0
        J1_arr[i] = (iy * Lx) + jx
        J2_arr[i] = (jy * Lx) + ix
        J3_arr[i] = (jy * Lx) + jx
      elif BoundaryCondition == 'OBCx + OBCy':
        jx = ix + 1 if ix + 1 < Lx else 0
        jy = iy + 1 if iy + 1 < Ly else 0
        J1_arr[i] = (iy * Lx) + jx
        J2_arr[i] = (jy * Lx) + ix
        J3_arr[i] = (jy * Lx) + jx
        if ix == Lx - 1:
          J1_arr[i] = -i
        if iy == Ly - 1:
          J2_arr[i] = -i
        if (iy == Ly - 1) or (ix == Lx - 1):
          J3_arr[i] = -i
      else:
        raise ValueError(BC_Error)
  return J1_arr, J2_arr, J3_arr
  
# ### Hamiltonian Matrix construction:
# 
# \begin{eqnarray}
# H & = & H_{\text{TB}} + H_{\text{Ex}} + H_{\text{sSc}} + H_{\text{pSc}} + H_{\text{dSc}} + H_{\text{ChemPot}} \nonumber \\
# & = & - \sum_{\langle ij \rangle,\sigma} t_{ij} (c^\dagger_{i\sigma} c^{}_{j\sigma} + {\textrm H.c.}) - J_{\text{H}} \sum_{i} {\bf S}_i \cdot {\bf s}_i + \sum_{i} \Delta^{s} (c^{\dagger}_{i \uparrow} c^{\dagger}_{i \downarrow} + h.c.) \nonumber \\
# & & + \sum_{ij, \sigma} \Delta^{p}_{ij} (c^{\dagger}_{i \sigma} c^{\dagger}_{j \sigma} + h.c.) + \sum_{ij} \Delta^{d}_{ij} (c^{\dagger}_{i \uparrow} c^{\dagger}_{j \downarrow} - c^{\dagger}_{i \downarrow} c^{\dagger}_{j \uparrow} + h.c.) + \mu \sum_{i, \sigma} c^{\dagger}_{i \sigma} c^{}_{i \sigma} \nonumber \\
# \end{eqnarray}
# 
# where,
# 
# $t_{i,i+\hat{x}} = t_x ; t_{i,i+\hat{y}} = t_y ; tx = ty$
# 
# $\Delta^{p}_{i,i+\hat{x}} = \Delta^{p}_{x} ; \Delta^{p}_{i,i+\hat{y}} = \Delta^{p}_{y} ; \Delta^{p}_{x} = \Delta^{p}_{y}$
# 
# $\Delta^{d}_{i,i+\hat{x}} = \Delta^{d}_{x} ; \Delta^{d}_{i,i+\hat{y}} = \Delta^{d}_{y} ; \Delta^{d}_{x} = -\Delta^{d}_{y}$


def hamiltonian(Sx, Sy, Sz):
  # Initialize arrays
  H = np.zeros((n, n), dtype=np.complex128)
  # setting up lattice geometry
  # J11, J22 = neighbor_array_square_lattice()
  J11, J22, J33 = neighbor_array_triangular_lattice()
  Shift = SDOF * Nsites
  
  # nearest neighbour tight-binding hopping along X (tx), and Y (ty)
  # TB = \sum_{ij, \sigma} -t_{ij} (c^{\dagger}_{i \sigma} c^{}_{j \sigma} + h.c.)
  for spindex2 in range (SDOF):
    for spindex1 in range (SDOF):
      for iy in range (Ly):
        for ix in range (Lx):
          if spindex1 == spindex2:
            ii = iy * Lx + ix
            i = (iy * Lx * SDOF) + (ix * SDOF) + spindex1
            if J11[ii] >= 0 and J11[ii] != ii:
              j1 = J11[ii]*SDOF + spindex2
              term = -tx * PHFactor
              H[i, j1] += term
              H[j1, i] += conj(H[i, j1])
              if PHDOF == 2:
                H[i + Shift, j1 + Shift] += -conj(term)
                H[j1 + Shift, i + Shift] += conj(H[i + Shift, j1 + Shift])
            if J22[ii] >= 0 and J22[ii] != ii:
              j2 = J22[ii]*SDOF + spindex2
              term = -ty * PHFactor
              H[i, j2] += term
              H[j2, i] += conj(H[i, j2])
              if PHDOF == 2:
                H[i + Shift, j2 + Shift] += -conj(term)
                H[j2 + Shift, i + Shift] += conj(H[i + Shift, j2 + Shift])
            if J33[ii] >= 0 and J33[ii] != ii:
              j3 = J33[ii]*SDOF + spindex2
              term = -ty * PHFactor
              H[i, j3] += term
              H[j3, i] += conj(H[i, j3])
              if PHDOF == 2:
                H[i + Shift, j3 + Shift] += -conj(term)
                H[j3 + Shift, i + Shift] += conj(H[i + Shift, j3 + Shift])

  # nearest neighbour p-wave SC pairing along X (Delta_px), and Y (Delta_py)
  # pSC = \sum_{ij, \sigma} \Delta^{p}{ij} (c^{\dagger}_{i \sigma} c^{\dagger}_{j \sigma} + h.c.)
  for spindex2 in range (SDOF):
    for spindex1 in range (SDOF):
      for iy in range (Ly):
        for ix in range (Lx):
          if spindex1 == spindex2:
            ii = iy * Lx + ix
            i = (iy * Lx * SDOF) + (ix * SDOF) + spindex1
            if J11[ii] >= 0 and J11[ii] != ii:
              j1 = J11[ii]*SDOF + spindex2
              term = Delta_px * PHFactor
              if PHDOF == 2:
                H[i, j1 + Shift] += term
                H[j1 + Shift, i] += conj(H[i, j1 + Shift])
                H[i + Shift, j1] += -conj(term)
                H[j1, i + Shift] += conj(H[i + Shift, j1])
            if J22[ii] >= 0 and J22[ii] != ii:
              j2 = J22[ii]*SDOF + spindex2
              term = Delta_py * PHFactor
              if PHDOF == 2:
                H[i, j2 + Shift] += term
                H[j2 + Shift, i] += conj(H[i, j2 + Shift])
                H[i + Shift, j2] += -conj(term)
                H[j2, i + Shift] += conj(H[i + Shift, j2])
            if J33[ii] >= 0 and J33[ii] != ii:
              j3 = J33[ii]*SDOF + spindex2
              term = Delta_py * PHFactor
              if PHDOF == 2:
                H[i, j3 + Shift] += term
                H[j3 + Shift, i] += conj(H[i, j3 + Shift])
                H[i + Shift, j3] += -conj(term)
                H[j3, i + Shift] += conj(H[i + Shift, j3])

  # chemical potential (\mu): \sum_{i, \sigma} c^{\dagger}_{i \sigma} c^{}_{i \sigma}
  for spindex2 in range (SDOF):
    for spindex1 in range (SDOF):
      for iy in range (Ly):
        for ix in range (Lx):
          if spindex1 == spindex2:
            i = (iy * Lx * SDOF) + (ix * SDOF) + spindex1
            j = (iy * Lx * SDOF) + (ix * SDOF) + spindex2
            term = -mu * PHFactor
            H[i, j] += term
            if PHDOF == 2:
              H[i + Shift, j + Shift] += -conj(term)

  # onsite s-wave SC pairing (Delta_s): \sum_{i} \Delta^{s} (c^{\dagger}_{i \uparrow} c^{\dagger}_{i \downarrow} + h.c.)
  for spindex2 in range (SDOF):
    for spindex1 in range (SDOF):
      for iy in range (Ly):
        for ix in range (Lx):
          if spindex1 == 0 and spindex2 == 1:
            i = (iy * Lx * SDOF) + (ix * SDOF) + spindex1
            j = (iy * Lx * SDOF) + (ix * SDOF) + spindex2
            term = PHFactor * Delta_s
            if PHDOF == 2:
              H[i, j + Shift] += term
              H[j + Shift, i] += conj(H[i, j + Shift])
              H[i + Shift, j] += -conj(term)
              H[j, i + Shift] += conj(H[i + Shift, j])

  # onsite exchange coupling between classical spins dof and itenary electron spin dof
  # $\text{Exchange Coupliong} = \frac{J_{\text{Hund}}}{2} \sum_{i} {\bf S}_i \cdot {\bf s}_i$
  for spindex2 in range (SDOF):
    for spindex1 in range (SDOF):
      for iy in range (Ly):
        for ix in range (Lx):
          i = (iy * Lx * SDOF) + (ix * SDOF) + spindex1
          j = (iy * Lx * SDOF) + (ix * SDOF) + spindex2
          # $ \mathbf{S}_i \times \mathbf{\sigma}_i  $
          term = - PHFactor * (J_Hundx * Sx[ix, iy] * PauliX[spindex1, spindex2] +
                               J_Hundy * Sy[ix, iy] * PauliY[spindex1, spindex2] +
                               J_Hundz * Sz[ix, iy] * PauliZ[spindex1, spindex2])
          H[i, j] += term
          if PHDOF == 2:
            H[i + Shift, j + Shift] += -conj(term)

  # nearest neighbour d-wave SC pairing:
  # $\text{dSC} = \sum_{ij} \Delta^{d}_{ij} (c^{\dagger}_{i \uparrow} c^{\dagger}_{j \downarrow} - 
  #                                          c^{\dagger}_{i \downarrow} c^{\dagger}_{j \uparrow} + h.c.) $
  for spindex2 in range (SDOF):
    for spindex1 in range (SDOF):
      for iy in range (Ly):
        for ix in range (Lx):
          ii = iy * Lx + ix
          if spindex1 == 0 and spindex2 == 1:
            i = (iy * Lx * SDOF) + (ix * SDOF) + spindex1
            if J11[ii] >= 0 and J11[ii] != ii:
              j1 = J11[ii]*SDOF + spindex2
              term = Delta_dx * PHFactor
              if PHDOF == 2:
                H[i, j1 + Shift] += term
                H[j1 + Shift, i] += conj(H[i, j1 + Shift])
                H[i + Shift, j1] += -conj(term)
                H[j1, i + Shift] += conj(H[i + Shift, j1])
            if J22[ii] >= 0 and J22[ii] != ii:
              j2 = J22[ii]*SDOF + spindex2
              term = Delta_dy * PHFactor
              if PHDOF == 2:
                H[i, j2 + Shift] += term
                H[j2 + Shift, i] += conj(H[i, j2 + Shift])
                H[i + Shift, j2] += -conj(term)
                H[j2, i + Shift] += conj(H[i + Shift, j2])
            if J33[ii] >= 0 and J33[ii] != ii:
              j3 = J33[ii]*SDOF + spindex2
              term = Delta_dy * PHFactor
              if PHDOF == 2:
                H[i, j3 + Shift] += term
                H[j3 + Shift, i] += conj(H[i, j3 + Shift])
                H[i + Shift, j3] += -conj(term)
                H[j3, i + Shift] += conj(H[i + Shift, j3])
          if spindex1 == 1 and spindex2 == 0:
            i = (iy * Lx * SDOF) + (ix * SDOF) + spindex1
            if J11[ii] >= 0 and J11[ii] != ii:
              j1 = J11[ii]*SDOF + spindex2
              term = -Delta_dx * PHFactor
              if PHDOF == 2:
                H[i, j1 + Shift] += term
                H[j1 + Shift, i] += conj(H[i, j1 + Shift])
                H[i + Shift, j1] += -conj(term)
                H[j1, i + Shift] += conj(H[i + Shift, j1])
            if J22[ii] >= 0 and J22[ii] != ii:
              j2 = J22[ii]*SDOF + spindex2
              term = -Delta_dy * PHFactor
              if PHDOF == 2:
                H[i, j2 + Shift] += term
                H[j2 + Shift, i] += conj(H[i, j2 + Shift])
                H[i + Shift, j2] += -conj(term)
                H[j2, i + Shift] += conj(H[i + Shift, j2])
            if J33[ii] >= 0 and J33[ii] != ii:
              j3 = J33[ii]*SDOF + spindex2
              term = -Delta_dy * PHFactor
              if PHDOF == 2:
                H[i, j3 + Shift] += term
                H[j3 + Shift, i] += conj(H[i, j3 + Shift])
                H[i + Shift, j3] += -conj(term)
                H[j3, i + Shift] += conj(H[i + Shift, j3])
  
  return H


def nambu_hamiltonian(Sx, Sy, Sz):
  # Initialize arrays
  H = np.zeros((n, n), dtype=np.complex128)
  # setting up lattice geometry
  # J11, J22 = neighbor_array_square_lattice()
  J11, J22, J33 = neighbor_array_triangular_lattice()
  # Pauli matrices for spin degrees of freedom
  sigma_x = np.array([[0, 1], [1, 0]])
  sigma_y = np.array([[0, -1j], [1j, 0]])
  sigma_z = np.array([[1, 0], [0, -1]])
  sigma_0 = np.identity(2)
  # Pauli matrices for particle-hole (PH) degrees of freedom
  tau_x = np.array([[0, 1], [1, 0]])
  tau_y = np.array([[0, -1j], [1j, 0]])
  tau_z = np.array([[1, 0], [0, -1]])
  tau_0 = np.identity(2)
  # Kronecker Product for PH and Spin space, [PH ⊗ Spin]
  Gamma1 = np.kron(tau_z, sigma_0)
  Gamma2 = np.kron(tau_x, sigma_0)
  Gamma3 = np.kron(tau_0, sigma_x)
  Gamma4 = np.kron(tau_0, sigma_y)
  Gamma5 = np.kron(tau_0, sigma_z)
  
  # The Hamiltonian is given in Ref. PRB 109, L041409 (2024)
  # \begin{eqnarray}
  #   H & = & \sum_{i,j} c^{\dagger}_{i,j}[\{\mu \Gamma_1 + \Delta_0 \Gamma_2 + J_{\text{Hund}} (S_x \Gamma_3 + S_y \Gamma_4 + S_z \Gamma_5) \}c^{}_{i,j} \nonumber \\
  #                                         & & - t \Gamma_1c^{}_{i+1,j} - t \Gamma_1c^{}_{i,j+1} ] + H.c.
  # \end{eqnarray}
  
  for iy in range (Ly):
    for ix in range (Lx):
        i = iy * Lx + ix
        for c in range (TDOF):
          for r in range (TDOF):
            # hopping along x
            if J11[i] >= 0 and J11[i] != i:
              j1 = J11[i]
              row = TDOF*i + r ; col = TDOF*j1 + c
              H[row,col] = (tx * Gamma1[r,c])
              H[col,row] = conj(H[row,col])
            # hopping along y
            if J22[i] >= 0 and J22[i] != i:
              j2 = J22[i]
              row = TDOF*i + r ; col = TDOF*j2 + c
              H[row,col] = (ty * Gamma1[r,c])
              H[col,row] = conj(H[row,col])
            # hopping along diagonal
            if J33[i] >= 0 and J33[i] != i:
              j3 = J33[i]
              row = TDOF*i + r ; col = TDOF*j3 + c
              H[row,col] = (td * Gamma1[r,c])
              H[col,row] = conj(H[row,col])
            # onsite terms, Chemical potential, s-Wave SC, exchange J_Hund
            row = TDOF*i + r ; col = TDOF*i + c
            H[row,col] = ((mu * Gamma1[r,c]) +
                          (Delta_s * Gamma2[r,c]) + 
                          (J_Hundx * Sx[ix, iy] * Gamma3[r,c] + 
                           J_Hundy * Sy[ix, iy] * Gamma4[r,c] + 
                           J_Hundz * Sz[ix, iy] * Gamma5[r,c]) )
            H[col,row] = conj(H[row,col])
          
  return H
  
  
def Bott_Index (Z):
  Eth = np.zeros((n, n), dtype=complex)
  Eph = np.zeros((n, n), dtype=complex)
  IdnttyMat = np.identity(4)
  
  # lattice coordinate projection onto spherical torus
  for iy in range (Ly):
    for ix in range (Lx):
      i = iy * Lx + ix
      for c in range (TDOF):
        for r in range (TDOF):
          row = TDOF*i + r ; col = TDOF*i + c
          Eth[row,col] = exp ( iota*(2.0*PI/Lx)*ix ) * IdnttyMat[r,c]
          Eph[row,col] = exp ( iota*(2.0*PI/Ly)*iy ) * IdnttyMat[r,c]
  
  # define projector
  Prj = np.zeros ((n,n), dtype=complex)  
  # for k in range (n//2):
  #   for i in range(n):
  #     for j in range(n):
  #       Prj[i, j] += Z[i, k] * conj(Z[j, k])
  
  # for k in range(n//2):
  #   # Select the k-th column and its conjugate transpose
  #   Z_k = Z[:, k]
  #   Z_k_conj = np.conj(Z_k)
    # Compute outer product and accumulate
    # Prj += np.outer(Z_k, Z_k_conj)
  Prj = np.dot(Z[:, :n//2], np.conj(Z[:, :n//2].T))
  
  # define position projected operators      
  U = np.dot (Prj, np.dot(Eth,Prj))
  V = np.dot (Prj, np.dot(Eph,Prj))
  
  # singular value decomposition of UMatrix and VMatrix
  
  # Bott matrix(W), VUV^{\dagger}U^{\dagger}
  W = np.dot(V, np.dot(U, np.dot(conj(V.T), conj(U.T))))
  
  # Compute the expression (1 / (2 * np.pi)) * Im[tr[log W]]
  BottIndex = (1 / (2.0 * PI)) * np.imag(np.trace(np.log(W)))
  
  return BottIndex


def BHZ_hamiltonian(Sx, Sy, Sz):
  # Initialize arrays
  H = np.zeros((n, n), dtype=np.complex128)
  # setting up lattice geometry
  J11, J22 = neighbor_array_square_lattice()
  # J11, J22, J33 = neighbor_array_triangular_lattice()  
  # Pauli matrices for orbital degrees of freedom
  sigma_x = np.array([[0, 1], [1, 0]])
  sigma_y = np.array([[0, -1j], [1j, 0]])
  sigma_z = np.array([[1, 0], [0, -1]])
  sigma_0 = np.identity(2)
  # Pauli matrices for spin degrees of freedom
  s_x = np.array([[0, 1], [1, 0]])
  s_y = np.array([[0, -1j], [1j, 0]])
  s_z = np.array([[1, 0], [0, -1]])
  s_0 = np.identity(2)
  # Pauli matrices for particle-hole degrees of freedom
  tau_x = np.array([[0, 1], [1, 0]])
  tau_y = np.array([[0, -1j], [1j, 0]])
  tau_z = np.array([[1, 0], [0, -1]])
  tau_0 = np.identity(2)
  # Kronecker Product for PH and Spin space, [PH ⊗ Spin]
  Gamma1 = np.kron(tau_z, np.kron(sigma_z, s_0))
  Gamma2 = np.kron(tau_x, np.kron(sigma_0, s_0))
  Gamma3 = np.kron(tau_0, np.kron(sigma_0, s_x))
  Gamma4 = np.kron(tau_0, np.kron(sigma_0, s_y))
  Gamma5 = np.kron(tau_0, np.kron(sigma_0, s_z))
  Gamma6 = np.kron(tau_z, np.kron(sigma_x, s_z))
  Gamma7 = np.kron(tau_z, np.kron(sigma_y, s_0))
  Gamma8 = np.kron(tau_z, np.kron(sigma_x, s_x))
  
  # The Hamiltonian is given in Ref. PRB 109, L041409 (2024)
  # \begin{eqnarray}
  #   H & = & \sum_{i,j} c^{\dagger}_{i,j}[\{\mu \Gamma_1 + \Delta_0 \Gamma_2 + J_{\text{Hund}} (S_x \Gamma_3 + S_y \Gamma_4 + S_z \Gamma_5) \}c^{}_{i,j} \nonumber \\
  #                                         & & - t \Gamma_1c^{}_{i+1,j} - t \Gamma_1c^{}_{i,j+1} ] + H.c.
  # \end{eqnarray}
  
  for iy in range (Ly):
    for ix in range (Lx):
      i = iy * Lx + ix
      for c in range (TDOF):
        for r in range (TDOF):
          # hopping along x
          if J11[i] >= 0 and J11[i] != i:
            j1 = J11[i]
            row = TDOF*i + r ; col = TDOF*j1 + c
            H[row,col] = ( - tx * Gamma1[r,c]
                           - 1j * Lambda_SOC_x * Gamma6[r,c]
                           + Lambda_WD_x * Gamma8[r,c])
            H[col,row] = conj(H[row,col])
          # hopping along y
          if J22[i] >= 0 and J22[i] != i:
            j2 = J22[i]
            row = TDOF*i + r ; col = TDOF*j2 + c
            H[row,col] = ( - ty * Gamma1[r,c]
                           - 1j * Lambda_SOC_y * Gamma7[r,c]
                           - Lambda_WD_y * Gamma8[r,c])
            H[col,row] = conj(H[row,col])
          # hopping along diagonal
          # if J33[i] >= 0 and J33[i] != i:
          #   j3 = J33[i]
          #   row = TDOF*i + r ; col = TDOF*j3 + c
          #   H[row,col] = ( -0.5 * td * Gamma1[r,c] -
          #                   1j * 0.5 * Lambda_SOC_d * Gamma4[r,c] -
          #                   0.5 * Lambda_WD_d * Gamma3[r,c] )
          #   H[col,row] = conj(H[row,col])
          # onsite terms, Chemical potential, s-Wave SC, exchange J_Hund
          row = TDOF*i + r ; col = TDOF*i + c
          H[row,col] = (epsilon0 * Gamma1[r,c] + 
                        Delta_s * Gamma2[r,c] + 
                        (J_Hundx * Sx[ix, iy] * Gamma3[r,c] + 
                         J_Hundy * Sy[ix, iy] * Gamma4[r,c] + 
                         J_Hundz * Sz[ix, iy] * Gamma5[r,c]))
          H[col,row] = conj(H[row,col])
          
  return H
  
def QuadrupoleMoment (Z):
  Idntty = np.identity(TDOF)
  
  U = np.zeros ((n, int(n/2)), dtype=complex)
  for j in range (int(n/2)):
    for i in range (n):
      U[i,j] = Z[i,j]
  
  q = np.zeros ((n, n), dtype=complex)
  for iy in range (Ly):
    for ix in range (Lx):
      i = iy * Lx + ix
      for c in range (TDOF):
        for r in range (TDOF):
          row = TDOF*i + r ; col = TDOF*i + c
          q[row, col] = (ix*iy*Idntty[r,c])/Nsites
            
  W = expm(1j*2.0*PI*q)
      
  V = np.dot(conj(U).T, np.dot(W, U))
  c1 = sp.linalg.det(V)
  c2 = exp(-1j*PI*np.trace(q))
  qxy = (1.0/(2.0*PI))*(ln(c1*c2))

  return Im(qxy)

  
  
  