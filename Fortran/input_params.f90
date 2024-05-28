MODULE input_parameters
  IMPLICIT NONE

  ! Lattice Type and lattice parameters
  INTEGER, PARAMETER :: Lx = 10 , Ly = 10
  INTEGER, PARAMETER :: Nsites = Lx * Ly
  CHARACTER(LEN=12), PARAMETER :: BoundaryConditionList(4) = ['PBCx + PBCy', 'OBCx + PBCy', 'PBCx + OBCy', 'OBCx + OBCy']
  CHARACTER(LEN=12), PARAMETER :: BoundaryCondition = BoundaryConditionList(4)

  ! Hamiltonian coupling parameters
  REAL, PARAMETER :: tx = 1.0 , ty = 1.0
  REAL, PARAMETER :: Delta_s = 0.4
  REAL, PARAMETER :: Delta_px = 0.0 , Delta_py = 0.0
  REAL, PARAMETER :: Delta_dx = 0.0 , Delta_dy = 0.0
  REAL, PARAMETER :: J_Hundx = 0.8 , J_Hundy = 0.8 , J_Hundz = 0.8
  REAL, PARAMETER :: mu = 0.0

  ! BHZ hamiltonian parameters
  REAL, PARAMETER :: epsilon_0 = 1.0
  REAL, PARAMETER :: Lambda_SOC_x = 0.5 , Lambda_SOC_y = 0.5
  REAL, PARAMETER :: Lambda_WD_x = 0.0 , Lambda_WD_y = 0.0

  ! Hamiltonian size parameters
  INTEGER, PARAMETER :: SDOF = 2 , ODOF = 2 , PHDOF = 2
  INTEGER, PARAMETER :: TDOF = PHDOF*ODOF*SDOF
  INTEGER, PARAMETER :: n = Nsites * ODOF * SDOF * PHDOF

  ! Spin Configuration Parameters
  COMPLEX, PARAMETER :: iota = CMPLX(0.0, 1.0)
  REAL, PARAMETER :: PI = 2.0 * ASIN(1.0)
  REAL, PARAMETER :: qx = 0.2 , qy = 0.2
  ! CHARACTER(LEN=20), PARAMETER :: SpinConfigType(7) = ["Skyrmion", "Spiral", "PRB_109_L041409", &
  !                                                      "AKG_code", "Ferro", "AntiFerro"]
  CHARACTER(LEN=*), PARAMETER :: SpinConfig = "Ferro"
  REAL, PARAMETER :: spin = 1.0
  INTEGER, PARAMETER :: SkyrRadius = 6
  INTEGER, PARAMETER :: NumSkAlongX = 1 , NumSkAlongY = 1

  REAL, PARAMETER :: Betaval = 0.03, vorticity = -1.0 , helicity = PI/2.0 , polarity = -1.0

  

  ! Skyrmion optimization paramters
  INTEGER, PARAMETER :: no_beta = 1000
  REAL, PARAMETER :: beta_min = 0.05 , beta_max = 0.3

  ! generic offset
  REAL, PARAMETER :: offset = 0.00001

  COMPLEX, PARAMETER :: PauliX(2,2) = RESHAPE([CMPLX(0.0,0.0), CMPLX(1.0,0.0), CMPLX(1.0,0.0), CMPLX(0.0,0.0)], [2, 2])
  COMPLEX, PARAMETER :: PauliY(2,2) = RESHAPE([CMPLX(0.0,0.0), CMPLX(0.0,1.0), CMPLX(0.0,-1.0), CMPLX(0.0,0.0)], [2, 2])
  COMPLEX, PARAMETER :: PauliZ(2,2) = RESHAPE([CMPLX(1.0,0.0), CMPLX(0.0,0.0), CMPLX(0.0,0.0), CMPLX(-1.0,0.0)], [2, 2])
  COMPLEX, PARAMETER :: Pauli0(2,2) = RESHAPE([CMPLX(1.0,0.0), CMPLX(0.0,0.0), CMPLX(0.0,0.0), CMPLX(1.0,0.0)], [2, 2])

  REAL, PARAMETER :: energy_of_interest = 0.0
  INTEGER, PARAMETER :: site_of_interest = 1
  REAL, PARAMETER :: eta = 1e-2
  INTEGER, PARAMETER :: DOS_NE = 5000, LDOS_NE = 5000

  CHARACTER(LEN=100), PARAMETER :: filename_spin = 'spin_config.txt'
  CHARACTER(LEN=100), PARAMETER :: filename_eigenval = 'eigenvals.txt'
  CHARACTER(LEN=100), PARAMETER :: filename_dos = 'dos.txt'
  CHARACTER(LEN=100), PARAMETER :: filename_ldos_E_fix = 'ldos_EnergyFixed.txt'
  CHARACTER(LEN=100), PARAMETER :: filename_ldos_Site_fix = 'ldos_SiteFixed.txt'

! BC_Error = 'Invalid BoundaryCondition! You have entered ' + BoundaryCondition + '\n' + 'Valid options are: PBCx + PBCy, OBCx + PBCy, PBCx + OBCy, OBCx + OBCy'
! line1 = 'Invalid SpinConfig_Type! You have entered ' + SpinConfig_Type + '\n'
! line2 = 'Valid options are: "NeelSkyrmion", "BlochSkyrmion", "AntiSkyrmion", "SpinSpiral", "Ferro", "AntiFerro"'
! SpinConfig_Error = line1 + line2

! CHARACTER (LEN=11), PARAMETER :: BoundaryCondition = 'PBCx + PBCy'
CHARACTER (LEN=200), PARAMETER :: BC_Error = 'Invalid BoundaryCondition! You have entered '&
                                              & //BoundaryCondition//NEW_LINE('A')// &
                                              & 'Valid options are: PBCx + PBCy, OBCx + PBCy, PBCx + OBCy, OBCx + OBCy'


END MODULE input_parameters