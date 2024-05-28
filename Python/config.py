from utils import *

# Lattice Type and lattice parameters
Lx: int = 20
Ly: int = 20
Nsites: int = Lx * Ly
BoundaryConditionList: list[str] = ['PBCx + PBCy', 'OBCx + PBCy', 'PBCx + OBCy', 'OBCx + OBCy']
BoundaryCondition: str = BoundaryConditionList[3]
BraviasLattice: str = "SquareLattice"

# Hamiltonian coupling parameters
# hopping parameters
tx: float = 1.0
ty: float = 1.0
td: float = 0.0
# superconducting order parameters
Delta_s: float = 0.4
Delta_px: float = 0.0
Delta_py: float = 0.0
Delta_dx: float = 0.0
Delta_dy: float = 0.0
# Hunds coupling parameters
J_Hundx: float = 0.8
J_Hundy: float = 0.8
J_Hundz: float = 0.8
# chemical potential & crystal field splitting
mu: float = 0.0
epsilon0: float = 1.0
# the SOC strength
Lambda_SOC_x: float = 0.5
Lambda_SOC_y: float = 0.5
Lambda_SOC_d: float = 0.0
Lambda_WD_x: float = 0.0
Lambda_WD_y: float = 0.0

# Hamiltonian size parameters
# degrees of freedom, spin (SDOF), particle-hole (PHDOF), orbital (ODOF)
SDOF: int = 2
ODOF: int = 2
PHDOF: int = 2
TDOF: int = PHDOF * ODOF * SDOF
PHFactor: float = 1.0 if PHDOF == 2 else 1.0
n: int = Nsites * TDOF

# Spin Configuration Parameters
SpinConfigTypeList: list[str] = [
    "NeelSkyrmion",  # 0
    "BlochSkyrmion",  # 1
    "AntiSkyrmion",  # 2
    "SpinSpiral",  # 3
    "Ferro",  # 4
    "AntiFerro",  # 5
    "CustomSkyrmion",  # 6
    "PRB_109_L041409",  # 7
    "PRR_4_013225_SWCB4"  # 8
]
SpinConfig_Type: str = SpinConfigTypeList[7]
SkyrDiameter: int = 20
vorticity: float = 1.0
helicity: float = 0.0
polarity: float = -1.0

qx: float = PI/5.0
qy: float = PI/5.0
Mx: float = 1.0
My: float = 1.0
Mz: float = 1.0
Bz: float = 0.5
# Skyrmion optimization parameters
no_Beta: int = 1000
Beta_min: float = 0.05
Beta_max: float = 0.3

# generic offset
offset: float = 0.00001

BC_Error: str = 'Invalid BoundaryCondition! You have entered ' + BoundaryCondition + '\n' + 'Valid options are: PBCx + PBCy, OBCx + PBCy, PBCx + OBCy, OBCx + OBCy'
line1: str = 'Invalid SpinConfig_Type! You have entered ' + SpinConfig_Type + '\n'
SpinConfig_Error: str = line1

def print_parameters(print_info: bool = True):
  if not print_info:
    return
  print("-" * 88)
  print("-" * 88)
  print("Parameters:")
  print(f"  Bravais Lattice: {BraviasLattice}")
  print(f"  Lattice Size: Lx = {Lx}, Ly = {Ly}")
  print(f"  Boundary Condition: {BoundaryCondition}")
  print(f"  Tight-Binding Parameters: tx = {tx}, ty = {ty}")
  print(f"  Chemical Potential: mu = {mu}")
  print(f"  S-wave pairing: Delta_s = {Delta_s}")
  print(f"  P-wave pairing: Delta_px = {Delta_px}, Delta_py = {Delta_py}")
  print(f"  D-Wave pairing: Delta_dx = {Delta_dx}, Delta_dy = {Delta_dy}")
  print(f"  Spin Configuration Type: {SpinConfig_Type}")
  if SpinConfig_Type == 'SpinSpiral':
      print(f"  Wave vector: qx = {qx}, qy = {qy}")
  elif 'Skyrmion' in SpinConfig_Type:
      print(f"  Skyrmion Parameters: SkyrDiameter = {SkyrDiameter}, SkyrRadius = {int(SkyrDiameter/2.0)}")
      print(f"  Number of Skyrmions Along X, Y: {int(Lx/SkyrDiameter)}, {int(Ly/SkyrDiameter)}")
      print(f"  Total Skyrmions = {int(Lx/SkyrDiameter) * int(Ly/SkyrDiameter)}")
  print("-" * 88)
  print("-" * 88)

energy_of_interest: float = 0.0
site_of_interest: int = 1
eta: float = 0.01  # 1e-2
DOS_NE: int = 200
LDOS_NE: int = 200

filename_spin: str = 'spin_config.txt'
filename_eigenval: str = 'eigenvals.txt'
filename_dos: str = 'dos.txt'
filename_ldos_E_fix: str = 'ldos_EnergyFixed.txt'
filename_ldos_Site_fix: str = 'ldos_SiteFixed.txt'
