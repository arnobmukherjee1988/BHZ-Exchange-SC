INCLUDE 'input_params.f90'
INCLUDE 'GenerateCoordinate.f90'
INCLUDE 'GenerateSpinConfig.f90'
INCLUDE 'GenerateHamiltonian.f90'
INCLUDE 'subroutines.f90'
INCLUDE 'Observables.f90'
! INCLUDE 'bott_index_svd.f90'
!*********************************************************************************************************************************
PROGRAM topoSC
  USE input_parameters
  USE gen_coordinate
  IMPLICIT NONE

  INTEGER :: ix, iy
  REAL, DIMENSION(lx,ly) :: Sx, Sy, Sz, Theta, Phi
  REAL :: skyr_number
  COMPLEX :: H (n,n), Z (n,n)
  REAL :: W(n)

  COMPLEX :: qxy

  REAL :: ldos_all_sites(lx,ly), ldos_all_energies(0:LDOS_NE)
  REAL :: EnergyMin, EnergyMax
  ! REAL :: DOSGreen
  ! COMPLEX :: GreenFuncMat(n, n), TRACE
  ! REAL :: BottIndex

  ! setting up skyrmion spin configuration
  CALL BuildSpinConfig (Sx, Sy, Sz, Theta, Phi)
  ! Sx(:,:) = 1.0 ; Sy(:,:) = 1.0 ; Sz(:,:) = 1.0
  DO iy = 0, Ly-1
    DO ix = 0, Lx-1
      Sx(ix+1, iy+1) = spin*SIN(PI/2.0)*COS(qx*ix + qy*iy)
      Sy(ix+1, iy+1) = spin*SIN(PI/2.0)*SIN(qx*ix + qy*iy)
      Sz(ix+1, iy+1) = spin*COS(PI/2.0)*0.0
      ! print*, ix, iy, Sx(ix+1, iy+1), Sy(ix+1, iy+1), Sz(ix+1, iy+1)
    ENDDO
  ENDDO
  CALL calculate_skyrmion_number(Sx, Sy, Sz, skyr_number)
  CALL print_SpinVectors (lx, ly, Sx, Sy, Sz, Theta, Phi, filename_spin)

  ! subroutines to construct the hamiltonian matrix and then diagonalize it to get eigenvalues W, eigenvectors Z
  ! Ev_Yes (Ev_No) ---> to calculate (not calculate) the eigenvectors
  CALL BHZ_Jex_sSC (Sx, Sy, Sz, H)
  ! CALL BHZ_hamiltonian (H)
  CALL diagonalization (n, H, 'Ev_Yes', W, Z)

  ! printing total eigen spectrum
  CALL print_eigenvals (n, W, TRIM(ADJUSTL(filename_eigenval)))

  ! Quadrupole Moment calculation for topological characterization
  ! CALL QuadrupoleMoment (Z, qxy)
  ! print*, qxy

  ! calculating and printing total DOS data
  ! CALL dos_cal (W, TRIM(ADJUSTL(filename_dos)))

  ! Calculate LDOS for all sites at the energy of interest and Print
  ! CALL LDOS_ForAllSites_AtFixedEnergy(energy_of_interest, W, Z, ldos_all_sites, TRIM(ADJUSTL(filename_ldos_E_fix)))

  ! Calculate LDOS for all sites at the energy of interest and Print
  ! EnergyMin = W(1) - 1.0 ; EnergyMax = W(n) + 1.0
  ! CALL LDOS_ForAllEnergies_AtFixedSite(site_of_interest, EnergyMin, EnergyMax, W, Z, &
  !                                      ldos_all_energies, TRIM(ADJUSTL(filename_ldos_Site_fix)))
  

  ! OPEN (UNIT = 20, FILE = 'dos_GreenFunction.txt', STATUS = 'unknown')
  ! EnergyMin = W(1) - 1.0 ; EnergyMax = W(n) + 1.0
  ! ! print*, EnergyMin, EnergyMax
  ! DO E_index = 0, DOS_NE
  !   EnergyVal = EnergyMin + ((EnergyMax - EnergyMin)/DOS_NE) * E_index
  !   CALL calculate_green_function (EnergyVal, W, Z, GreenFuncMat)
  !   ! DOSGreen = 0.0
  !   ! GreenFuncMat_Trace = TRACE (n, GreenFuncMat)
  !   DOSGreen = - AIMAG(TRACE (n, GreenFuncMat))/pi
  !   WRITE(20,*) EnergyVal, DOSGreen
  ! ENDDO



END PROGRAM topoSC
!*********************************************************************************************************************************


