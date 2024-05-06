SUBROUTINE dos_cal (E, filename)    !sys_size=system size(n*pi), n=energy dimension
  USE input_parameters
  IMPLICIT NONE

  REAL, INTENT(IN) :: E(n)
  CHARACTER(LEN=*), INTENT(IN) :: filename
  REAL :: dos, Energy, E_min, E_max
  INTEGER :: i, j

  OPEN (UNIT = 2, FILE = TRIM(ADJUSTL(filename)), STATUS = 'unknown')
  E_min = E(1) - 1.0
  E_max = E(n) + 1.0
  DO i = 0, DOS_NE
    dos = 0.0
    Energy = ((E_max - E_min)/DOS_NE) * i + E_min
    DO j = 1, n
      dos = dos + (1.0/(pi*n)) * ( eta/((Energy-E(j))**2 + eta**2) )
    ENDDO
    WRITE(2,*) Energy, dos
  ENDDO

END SUBROUTINE dos_cal


! SUBROUTINE LDOS_ForAllSites_AtFixedEnergy(E, W, Z, ldos, filename)
!   USE input_parameters
!   IMPLICIT NONE

!   REAL, INTENT(IN) :: E, W(n)
!   COMPLEX, INTENT(IN) :: Z(n, n)
!   REAL, INTENT(OUT) :: ldos(lx, ly)
!   CHARACTER (LEN=*), INTENT(IN) :: filename
!   INTEGER :: i,j, ix, iy, k, spindex
!   CHARACTER (LEN=20) :: strE
!   REAL :: Lorentzian

!   WRITE (strE,'(F10.4)') E
!   OPEN (UNIT = 3, FILE = 'E_'//TRIM(ADJUSTL(strE))//'_'//TRIM(ADJUSTL(filename)), STATUS = 'unknown') 
!   ldos(:,:) = 0.0
!   DO iy = 1, ly
!     DO ix = 1, lx
!       DO spindex = 1, SDOF
!         i = (iy-1)*lx*SDOF + (ix-1)*SDOF + spindex
!         j = i + Nsites*SDOF
!         DO k = 1, n
!           ldos(ix,iy) = ldos(ix,iy) + Lorentzian(eta, E-W(k)) * ABS( Z(i,k) * CONJG(Z(i,k)) + Z(j,k) * CONJG(Z(j,k)) )
!         ENDDO
!       ENDDO
!     ENDDO
!   ENDDO

!   DO iy = 1, ly
!     DO ix = 1, lx
!       WRITE(3,*) ix, iy, ldos(ix, iy)/ MAXVAL(ldos)
!     ENDDO
!     WRITE(3,*)
!   ENDDO 

! END SUBROUTINE LDOS_ForAllSites_AtFixedEnergy



SUBROUTINE LDOS_ForAllSites_AtFixedEnergy(E, W, Z, ldos, filename)
  USE input_parameters
  IMPLICIT NONE

  REAL, INTENT(IN) :: E, W(n)
  COMPLEX, INTENT(IN) :: Z(n, n)
  REAL, INTENT(OUT) :: ldos(lx, ly)
  CHARACTER (LEN=*), INTENT(IN) :: filename
  INTEGER :: i, ix, iy, k, r, c, row, col
  CHARACTER (LEN=20) :: strE
  REAL :: Lorentzian

  WRITE (strE,'(F10.4)') E
  OPEN (UNIT = 3, FILE = 'E_'//TRIM(ADJUSTL(strE))//'_'//TRIM(ADJUSTL(filename)), STATUS = 'unknown') 
  ldos(:,:) = 0.0
  DO k = 1, n
    DO iy = 1, Ly
      DO ix = 1, Lx
        i = (iy-1) * Lx + ix
        DO c = 1, TDOF
          DO r = 1, TDOF
            row = TDOF*(i-1) + r ; col = TDOF*(i-1) + c
            ldos(ix,iy) = ldos(ix,iy) + Lorentzian(eta, E-W(k)) * ABS( Z(row,k) * CONJG(Z(col,k)) )
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO iy = 1, ly
    DO ix = 1, lx
      WRITE(3,*) ix, iy, ldos(ix, iy)/ MAXVAL(ldos)
    ENDDO
    WRITE(3,*)
  ENDDO 

END SUBROUTINE LDOS_ForAllSites_AtFixedEnergy

! SUBROUTINE LDOS_ForAllEnergies_AtFixedSite(site, EnergyMin, EnergyMax, W, Z, ldos, filename)
!   USE input_parameters
!   USE gen_coordinate
!   IMPLICIT NONE

!   INTEGER, INTENT (IN) :: site
!   REAL, INTENT(IN) :: EnergyMax, EnergyMin, W(n)
!   COMPLEX, INTENT(IN) :: Z(n, n)
!   REAL, INTENT(OUT) :: ldos(0:LDOS_NE)
!   CHARACTER (LEN=*), INTENT(IN) :: filename
!   REAL :: E, Lorentzian
!   INTEGER :: i, j, ix, iy, k, spindex, E_index
!   CHARACTER (LEN=20) :: strSITE

!   CALL GetSiteCoordinates(site, ix, iy)
!   ldos(:) = 0.0
!   DO E_index = 0, LDOS_NE
!     E = EnergyMin + ((EnergyMax - EnergyMin)/LDOS_NE) * E_index
!     DO spindex = 1, SDOF
!       i = (iy-1)*lx*SDOF + (ix-1)*SDOF + spindex
!       j = i + Nsites*SDOF
!       DO k = 1, n
!         ldos(E_index) = ldos(E_index) + Lorentzian(eta, E-W(k)) * ABS( Z(i,k) * CONJG(Z(i,k)) + Z(j,k) * CONJG(Z(j,k)) )
!       ENDDO
!     ENDDO
!   ENDDO

!   WRITE (strSITE,'(I6)') site_of_interest
!   OPEN (UNIT = 4, FILE = 'Site_'//TRIM(ADJUSTL(strSITE))//'_'//TRIM(ADJUSTL(filename)), STATUS = 'unknown')
!   DO E_index = 0, LDOS_NE
!     E = EnergyMin + ((EnergyMax - EnergyMin)/LDOS_NE) * E_index
!     WRITE(4,*) E, ldos(E_index)/ MAXVAL(ldos)
!   ENDDO

! END SUBROUTINE LDOS_ForAllEnergies_AtFixedSite

SUBROUTINE LDOS_ForAllEnergies_AtFixedSite(site, EnergyMin, EnergyMax, W, Z, ldos, filename)
  USE input_parameters
  USE gen_coordinate
  IMPLICIT NONE

  INTEGER, INTENT (IN) :: site
  REAL, INTENT(IN) :: EnergyMax, EnergyMin, W(n)
  COMPLEX, INTENT(IN) :: Z(n, n)
  REAL, INTENT(OUT) :: ldos(0:LDOS_NE)
  CHARACTER (LEN=*), INTENT(IN) :: filename
  REAL :: E, Lorentzian
  INTEGER :: i, k, c, r, col, row, E_index
  CHARACTER (LEN=20) :: strSITE

  i = site
  ldos(:) = 0.0
  DO E_index = 0, LDOS_NE
    E = EnergyMin + ((EnergyMax - EnergyMin)/LDOS_NE) * E_index
    DO k = 1, n
      DO c = 1, TDOF
        DO r = 1, TDOF
          row = TDOF*(i-1) + r ; col = TDOF*(i-1) + c
          ldos(E_index) = ldos(E_index) + Lorentzian(eta, E-W(k)) * ABS( Z(row,k) * CONJG(Z(col,k)) )
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  WRITE (strSITE,'(I6)') site_of_interest
  OPEN (UNIT = 4, FILE = 'Site_'//TRIM(ADJUSTL(strSITE))//'_'//TRIM(ADJUSTL(filename)), STATUS = 'unknown')
  DO E_index = 0, LDOS_NE
    E = EnergyMin + ((EnergyMax - EnergyMin)/LDOS_NE) * E_index
    ! WRITE(4,*) E, ldos(E_index)/ MAXVAL(ldos)
    WRITE(4,*) E, ldos(E_index)/ (Nsites)
  ENDDO

END SUBROUTINE LDOS_ForAllEnergies_AtFixedSite

SUBROUTINE calculate_green_function (E, W, Z, GreenFuncMat)
  USE input_parameters
  IMPLICIT NONE

  REAL, INTENT(IN) :: E
  REAL, INTENT(IN) :: W(n)
  COMPLEX, INTENT(IN) :: Z(n, n)
  COMPLEX, INTENT(OUT) :: GreenFuncMat(n, n)
  COMPLEX, SAVE :: diag_inv_denominator(n,n)
  INTEGER :: i, j, k

  diag_inv_denominator(:,:) = CMPLX (0.0, 0.0)
  DO i = 1, n
    diag_inv_denominator (i,i) = 1.0 / (E - W(i) + CMPLX (0.0, eta))
  ENDDO

  GreenFuncMat(:,:) = CMPLX (0.0, 0.0)
  DO i = 1, n
    DO j = 1, n
      DO k = 1, n
        GreenFuncMat(i,j) = GreenFuncMat(i,j) + Z(i, k) * diag_inv_denominator(k, k) * CONJG(Z(j, k))
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE calculate_green_function

! SUBROUTINE QuadrupoleMoment(Z, qxy)
!   complex(8), dimension(:,:), intent(in) :: Z
!   real(8), intent(out) :: qxy
  
!   complex(8) :: Idntty(TDOF,TDOF)
!   complex(8) :: U(n, n/2), q(n, n)
!   complex(8) :: W(TDOF*Nsites, TDOF*Nsites)
!   complex(8) :: V(n/2, n/2)
!   complex(8) :: c1, c2
  
!   integer :: iy, ix, i, c, r, row, col
  
!   ! Identity matrix
!   Idntty = 0.0
!   do i = 1, TDOF
!       Idntty(i, i) = 1.0
!   end do
  
!   ! Populate U matrix
!   U = 0.0
!   do j = 1, n/2
!       do i = 1, n
!           U(i, j) = Z(i, j)
!       end do
!   end do
  
!   ! Populate q matrix
!   q = 0.0
!   do iy = 1, Ly
!       do ix = 1, Lx
!           i = (iy - 1) * Lx + ix
!           do c = 1, TDOF
!               do r = 1, TDOF
!                   row = TDOF * (i - 1) + r
!                   col = TDOF * (i - 1) + c
!                   q(row, col) = (ix * iy * Idntty(r, c)) / Nsites
!               end do
!           end do
!       end do
!   end do
  
!   ! Calculate W matrix
!   call expMatrix(1j * 2.0 * PI * q, W)
  
!   ! Calculate V matrix
!   call matmul(conjg(transpose(U)), matmul(W, U), V)
  
!   ! Calculate c1 and c2
!   c1 = det(V)
!   c2 = exp(-1j * PI * trace(q))
  
!   ! Calculate qxy
!   qxy = (1.0 / (2.0 * PI)) * aimag(log(c1 * c2))
  
! contains

!   subroutine expMatrix(A, ExpA)
!       complex(8), dimension(:,:), intent(in) :: A
!       complex(8), dimension(size(A,1), size(A,2)), intent(out) :: ExpA
      
!       ! Your implementation of matrix exponential here
!       ! For simplicity, you can assume ExpA = exp(A)
      
!   end subroutine expMatrix

! end subroutine QuadrupoleMoment


SUBROUTINE QuadrupoleMoment (Z, qxy)
  USE input_parameters
  IMPLICIT NONE

  COMPLEX, INTENT(IN) :: Z(n,n)
  COMPLEX, INTENT(OUT) :: qxy

  INTEGER :: i, j, ix, iy, c, r, row, col
  COMPLEX :: Idntty(TDOF,TDOF)
  COMPLEX :: U(n, n/2), q(n, n), W(n, n)
  COMPLEX :: conj_U_transpose(n/2, n), V(n/2, n/2)
  COMPLEX :: c1, c2, trace_q
  COMPLEX :: Z1 (n/2,n/2)
  REAL :: W1(n/2)

  Idntty (:,:) = CMPLX (0.0, 0.0) 
  DO i = 1, TDOF
    Idntty (i,i) = CMPLX (1.0, 0.0)
  ENDDO
  
  DO j = 1, INT(n/2)
    DO i = 1, n
      U(i,j) = Z(i,j)
    ENDDO
  ENDDO
  
  DO iy = 1, Ly
    DO ix = 1, Lx
      i = (iy-1) * Lx + ix
      DO c = 1, TDOF
        DO r = 1, TDOF
          row = TDOF*(i-1) + r ; col = TDOF*(i-1) + c
          q(row, col) = ((ix)*(iy)*Idntty(r,c))/Nsites
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  CALL compute_exp_matrix(n, q, W)
  DO j = 1, n
    DO i = 1, n
      W(i,j) = iota * 2.0 * PI * W(i,j)
    ENDDO
  ENDDO

  ! W = expm(1j*2.0*PI*q)
  ! V = np.dot(conj(U).T, np.dot(W, U))
  ! c1 = sp.linalg.det(V)
  ! c2 = exp(-1j*PI*np.trace(q))
  ! qxy = (1.0/(2.0*PI))*(ln(c1*c2))

  DO i = 1, n
    DO j = 1, INT(n/2)
      conj_U_transpose(j,i) = CONJG(U(i,j))
    ENDDO
  ENDDO
  V = MATMUL (conj_U_transpose, MATMUL (W, U))
  
  CALL diagonalization (INT(n/2), V, 'Ev_No', W1, Z1)
  c1 = 1.0
  DO i = 1, n/2
    c1 = c1 * W1(i)
  ENDDO
  
  trace_q = 0.0
  DO i = 1, n
    trace_q = trace_q + q(i,i)
  ENDDO
  c2 = EXP (-iota * PI * trace_q)

  qxy = (1.0/(2.0*PI)) * LOG (c1*c2)

END SUBROUTINE QuadrupoleMoment