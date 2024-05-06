! ### Hamiltonian Matrix construction:

! \begin{eqnarray}
! H & = & H_{\text{TB}} + H_{\text{Ex}} + H_{\text{sSc}} + H_{\text{pSc}} + H_{\text{dSc}} + H_{\text{ChemPot}} \nonumber \\
! & = & - \sum_{\langle ij \rangle,\sigma} t_{ij} (c^\dagger_{i\sigma} c^{}_{j\sigma} + {\textrm H.c.}) - J_{\text{H}} \sum_{i} {\bf S}_i \cdot {\bf s}_i + \sum_{i} \Delta^{s} (c^{\dagger}_{i \uparrow} c^{\dagger}_{i \downarrow} + h.c.) \nonumber \\
! & & + \sum_{ij, \sigma} \Delta^{p}_{ij} (c^{\dagger}_{i \sigma} c^{\dagger}_{j \sigma} + h.c.) + \sum_{ij} \Delta^{d}_{ij} (c^{\dagger}_{i \uparrow} c^{\dagger}_{j \downarrow} - c^{\dagger}_{i \downarrow} c^{\dagger}_{j \uparrow} + h.c.) + \mu \sum_{i, \sigma} c^{\dagger}_{i \sigma} c^{}_{i \sigma} \nonumber \\
! \end{eqnarray}

! where,

! $t_{i,i+\hat{x}} = t_x ; t_{i,i+\hat{y}} = t_y ; tx = ty$

! $\Delta^{p}_{i,i+\hat{x}} = \Delta^{p}_{x} ; \Delta^{p}_{i,i+\hat{y}} = \Delta^{p}_{y} ; \Delta^{p}_{x} = \Delta^{p}_{y}$

! $\Delta^{d}_{i,i+\hat{x}} = \Delta^{d}_{x} ; \Delta^{d}_{i,i+\hat{y}} = \Delta^{d}_{y} ; \Delta^{d}_{x} = -\Delta^{d}_{y}$


SUBROUTINE BuildHamiltonian (Sx, Sy, Sz, H)
  USE input_parameters
  USE gen_coordinate
  IMPLICIT NONE

  REAL, DIMENSION(Lx,Ly), INTENT(IN) :: Sx, Sy, Sz
  COMPLEX, INTENT(OUT) :: H (n,n)
  
  INTEGER :: i, ii, j, ix, iy, spindex1, spindex2, Shift
  INTEGER, DIMENSION(Lx*Ly) :: J11, J22
  REAL :: PHFactor
  COMPLEX :: term

  IF (PHDOF .EQ. 2) THEN
    PHFactor = 1.0
  ELSEIF (PHDOF .EQ. 1) THEN
    PHFactor = 1.0
  ENDIF

  CALL NeighborArray (J11, J22)
  Shift = SDOF * Nsites

  H(:,:) = CMPLX(0.0,0.0)

  ! nearest neighbour tight-binding hopping along X (tx), and Y (ty)
  ! TB = \sum_{ij, \sigma} -t_{ij} (c^{\dagger}_{i \sigma} c^{}_{j \sigma} + h.c.)
  DO spindex2 = 1, SDOF
    DO spindex1 = 1, SDOF
      DO iy = 1, Ly
        DO ix = 1, Lx

          ii = (iy-1) * Lx + ix
          IF (spindex1 .EQ. spindex2) THEN
            
            i = (iy-1)*Lx*SDOF + (ix-1)*SDOF + spindex1
            IF ( (J11(ii) .GT. 0) .AND. (J11(ii) .NE. ii) ) THEN
              j = (J11(ii)-1)*SDOF + spindex2 ; term = -tx * PHFactor * CMPLX(1.0, 0.0)
              H(i,j) = H(i,j) + term
              H(j,i) = H(j,i) + CONJG(H(i,j))
              IF (PHDOF .EQ. 2) THEN
                H(i+Shift, j+Shift) = H(i+Shift, j+Shift) - CONJG(term)
                H(j+Shift, i+Shift) = H(j+Shift, i+Shift) + CONJG(H(i+Shift, j+Shift))
              ENDIF
            ENDIF

            IF ( (J22(ii) .GT. 0) .AND. (J22(ii) .NE. ii) ) THEN
              j = (J22(ii)-1)*SDOF + spindex2 ; term = -ty * PHFactor * CMPLX(1.0, 0.0)
              H(i,j) = H(i,j) + term
              H(j,i) = H(j,i) + CONJG(H(i,j))
              IF (PHDOF .EQ. 2) THEN
                H(i+Shift, j+Shift) = H(i+Shift, j+Shift) - CONJG(term)
                H(j+Shift, i+Shift) = H(j+Shift, i+Shift) + CONJG(H(i+Shift, j+Shift))
              ENDIF
            ENDIF

          ENDIF

        ENDDO
      ENDDO  
    ENDDO
  ENDDO

  ! nearest neighbour p-wave SC pairing along X (Delta_px), and Y (Delta_py)
  ! pSC = \sum_{ij, \sigma} \Delta^{p}{ij} (c^{\dagger}_{i \sigma} c^{\dagger}_{j \sigma} + h.c.)
  DO spindex2 = 1, SDOF
    DO spindex1 = 1, SDOF
      DO iy = 1, Ly
        DO ix = 1, Lx

          ii = (iy-1) * Lx + ix
          IF (spindex1 .EQ. spindex2) THEN
            i = (iy-1)*Lx*SDOF + (ix-1)*SDOF + spindex1
            
            IF ( (J11(ii) .GT. 0) .AND. (J11(ii) .NE. ii) ) THEN
              j = (J11(ii)-1)*SDOF + spindex2 ; term = Delta_px * PHFactor * CMPLX(1.0, 0.0)
              IF (PHDOF .EQ. 2) THEN
                H(i, j+Shift) = H(i, j+Shift) + term
                H(j+Shift, i) = H(j+Shift, i) + CONJG(H(i, j+Shift))
                H(i+Shift, j) = H(i+Shift, j) - CONJG(term)
                H(j, i+Shift) = H(j, i+Shift) + CONJG(H(i+Shift, j))
              ENDIF
            ENDIF
            
            IF ( (J22(ii) .GT. 0) .AND. (J22(ii) .NE. ii) ) THEN
              j = (J22(ii)-1)*SDOF + spindex2 ; term = Delta_py * PHFactor * CMPLX(1.0, 0.0)
              IF (PHDOF .EQ. 2) THEN
                H(i, j+Shift) = H(i, j+Shift) + term
                H(j+Shift, i) = H(j+Shift, i) + CONJG(H(i, j+Shift))
                H(i+Shift, j) = H(i+Shift, j) - CONJG(term)
                H(j, i+Shift) = H(j, i+Shift) + CONJG(H(i+Shift, j))
              ENDIF
            ENDIF

          ENDIF

        ENDDO
      ENDDO  
    ENDDO
  ENDDO

  ! chemical potential (\mu): \sum_{i, \sigma} c^{\dagger}_{i \sigma} c^{}_{i \sigma}
  DO spindex2 = 1, SDOF
    DO spindex1 = 1, SDOF
      DO iy = 1, Ly
        DO ix = 1, Lx
          
          IF (spindex1 .EQ. spindex2) THEN
            i = (iy-1)*Lx*SDOF + (ix-1)*SDOF + spindex1
            j = (iy-1)*Lx*SDOF + (ix-1)*SDOF + spindex2
            term = -mu * PHFactor * CMPLX(1.0, 0.0)
            H(i, j) = H(i, j) + term
            IF (PHDOF .EQ. 2) H(i+Shift, j+Shift) = H(i+Shift, j+Shift) - CONJG(term)
          ENDIF

        ENDDO
      ENDDO  
    ENDDO
  ENDDO

  ! onsite s-wave SC pairing (Delta_s): \sum_{i} \Delta^{s} (c^{\dagger}_{i \uparrow} c^{\dagger}_{i \downarrow} + h.c.)
  DO spindex2 = 1, SDOF
    DO spindex1 = 1, SDOF
      DO iy = 1, Ly
        DO ix = 1, Lx
          
          IF ( (spindex1 .EQ. 1) .AND. (spindex2 .EQ. 2) )THEN
            i = (iy-1)*Lx*SDOF + (ix-1)*SDOF + spindex1
            j = (iy-1)*Lx*SDOF + (ix-1)*SDOF + spindex2
            term = Delta_s * PHFactor * CMPLX(1.0, 0.0)
            IF (PHDOF .EQ. 2) THEN
              H(i, j+Shift) = H(i, j+Shift) + term
              H(j+Shift, i) = H(j+Shift, i) + CONJG(H(i, j+Shift))
              H(i+Shift, j) = H(i+Shift, j) - CONJG(term)
              H(j, i+Shift) = H(j, i+Shift) + CONJG(H(i+Shift, j))
            ENDIF
          ENDIF

        ENDDO
      ENDDO  
    ENDDO
  ENDDO

  ! onsite exchange coupling between classical spins dof and itenary electron spin dof
  ! $\text{Exchange Coupliong} = \frac{J_{\text{Hund}}}{2} \sum_{i} {\bf S}_i \cdot {\bf s}_i$
  DO spindex2 = 1, SDOF
    DO spindex1 = 1, SDOF
      DO iy = 1, Ly
        DO ix = 1, Lx
            i = (iy-1)*Lx*SDOF + (ix-1)*SDOF + spindex1
            j = (iy-1)*Lx*SDOF + (ix-1)*SDOF + spindex2
            term = - PHFactor * ( J_Hundx * Sx(ix, iy) * PauliX(spindex1, spindex2) + &
                                  J_Hundy * Sy(ix, iy) * PauliY(spindex1, spindex2) + &
                                  J_Hundz * Sz(ix, iy) * PauliZ(spindex1, spindex2) )
            H(i, j) = H(i, j) + term
            IF (PHDOF .EQ. 2) H(i+Shift, j+Shift) = H(i+Shift, j+Shift) - CONJG(term)
        ENDDO
      ENDDO  
    ENDDO
  ENDDO

  ! nearest neighbour d-wave SC pairing:
  ! $\text{dSC} = \sum_{ij} \Delta^{d}_{ij} (c^{\dagger}_{i \uparrow} c^{\dagger}_{j \downarrow} - 
  !                                          c^{\dagger}_{i \downarrow} c^{\dagger}_{j \uparrow} + h.c.) $
  DO spindex2 = 1, SDOF
    DO spindex1 = 1, SDOF
      DO iy = 1, Ly
        DO ix = 1, Lx

          ii = (iy-1) * Lx + ix
          IF ( (spindex1 .EQ. 1) .AND. (spindex2 .EQ. 2) ) THEN
            i = (iy-1)*Lx*SDOF + (ix-1)*SDOF + spindex1
            
            IF ( (J11(ii) .GT. 0) .AND. (J11(ii) .NE. ii) ) THEN
              j = (J11(ii)-1)*SDOF + spindex2 ; term = Delta_dx * PHFactor * CMPLX(1.0, 0.0)
              IF (PHDOF .EQ. 2) THEN
                H(i, j+Shift) = H(i, j+Shift) + term
                H(j+Shift, i) = H(j+Shift, i) + CONJG(H(i, j+Shift))
                H(i+Shift, j) = H(i+Shift, j) - CONJG(term)
                H(j, i+Shift) = H(j, i+Shift) + CONJG(H(i+Shift, j))
              ENDIF
            ENDIF
            
            IF ( (J22(ii) .GT. 0) .AND. (J22(ii) .NE. ii) ) THEN
              j = (J22(ii)-1)*SDOF + spindex2 ; term = Delta_py * PHFactor * CMPLX(1.0, 0.0)
              IF (PHDOF .EQ. 2) THEN
                H(i, j+Shift) = H(i, j+Shift) + term
                H(j+Shift, i) = H(j+Shift, i) + CONJG(H(i, j+Shift))
                H(i+Shift, j) = H(i+Shift, j) - CONJG(term)
                H(j, i+Shift) = H(j, i+Shift) + CONJG(H(i+Shift, j))
              ENDIF
            ENDIF

          ENDIF

          IF ( (spindex1 .EQ. 2) .AND. (spindex2 .EQ. 1) ) THEN
            i = (iy-1)*Lx*SDOF + (ix-1)*SDOF + spindex1
            
            IF ( (J11(ii) .GT. 0) .AND. (J11(ii) .NE. ii) ) THEN
              j = (J11(ii)-1)*SDOF + spindex2 ; term = -Delta_dx * PHFactor * CMPLX(1.0, 0.0)
              IF (PHDOF .EQ. 2) THEN
                H(i, j+Shift) = H(i, j+Shift) + term
                H(j+Shift, i) = H(j+Shift, i) + CONJG(H(i, j+Shift))
                H(i+Shift, j) = H(i+Shift, j) - CONJG(term)
                H(j, i+Shift) = H(j, i+Shift) + CONJG(H(i+Shift, j))
              ENDIF
            ENDIF
            
            IF ( (J22(ii) .GT. 0) .AND. (J22(ii) .NE. ii) ) THEN
              j = (J22(ii)-1)*SDOF + spindex2 ; term = -Delta_py * PHFactor * CMPLX(1.0, 0.0)
              IF (PHDOF .EQ. 2) THEN
                H(i, j+Shift) = H(i, j+Shift) + term
                H(j+Shift, i) = H(j+Shift, i) + CONJG(H(i, j+Shift))
                H(i+Shift, j) = H(i+Shift, j) - CONJG(term)
                H(j, i+Shift) = H(j, i+Shift) + CONJG(H(i+Shift, j))
              ENDIF
            ENDIF

          ENDIF

        ENDDO
      ENDDO  
    ENDDO
  ENDDO

END SUBROUTINE BuildHamiltonian


SUBROUTINE BuildHamiltonian_In_NumbuSpinor (Sx, Sy, Sz, H)
  USE input_parameters
  USE gen_coordinate
  IMPLICIT NONE

  REAL, DIMENSION(Lx,Ly), INTENT(IN) :: Sx, Sy, Sz
  COMPLEX, INTENT(OUT) :: H (n,n)
  
  INTEGER :: i, j1, j2, ix, iy, c, r, row, col
  INTEGER, DIMENSION(Lx*Ly) :: J11, J22
  COMPLEX, DIMENSION(SDOF, SDOF) :: sigma_0, sigma_x, sigma_y, sigma_z
  COMPLEX, DIMENSION(PHDOF, PHDOF) :: tau_0, tau_x, tau_y, tau_z
  COMPLEX, DIMENSION(SDOF*PHDOF,SDOF*PHDOF) :: Gamma1, Gamma2, Gamma3, Gamma4, Gamma5

  CALL NeighborArray (J11, J22)
  H(:,:) = CMPLX(0.0,0.0)
  ! Pauli matrices for spin (sigma) and particle-hole (tau) degrees of freedom
  sigma_x = PauliX ; tau_x = PauliX
  sigma_y = PauliY ; tau_y = PauliY
  sigma_z = PauliZ ; tau_z = PauliZ
  sigma_0 = Pauli0 ;tau_0 = Pauli0  
  ! Kronecker Product for PH and Spin space, [PH ⊗ Spin]
  CALL kronecker_product2(PHDOF, SDOF, tau_z, sigma_0, Gamma1)
  CALL kronecker_product2(PHDOF, SDOF, tau_x, sigma_0, Gamma2)
  CALL kronecker_product2(PHDOF, SDOF, tau_0, sigma_x, Gamma3)
  CALL kronecker_product2(PHDOF, SDOF, tau_0, sigma_y, Gamma4)
  CALL kronecker_product2(PHDOF, SDOF, tau_0, sigma_z, Gamma5)

  ! The Hamiltonian is given in Ref. PRB 109, L041409 (2024)
  ! \begin{eqnarray}
  !   H & = & \sum_{i,j} c^{\dagger}_{i,j}[\{\mu \Gamma_1 + \Delta_0 \Gamma_2 + J_{\text{Hund}} (S_x \Gamma_3 + S_y \Gamma_4 + S_z \Gamma_5) \}c^{}_{i,j} \nonumber \\
  !                                         & & - t \Gamma_1c^{}_{i+1,j} - t \Gamma_1c^{}_{i,j+1} ] + H.c.
  ! \end{eqnarray}

  DO iy = 1, Ly
    DO ix = 1, Lx
      i = (iy-1) * Lx + ix
      DO c = 1, SDOF*PHDOF
        DO r = 1, SDOF*PHDOF
          ! hopping along x
          IF ( (J11(i) .GT. 0) .AND. (J11(i) .NE. i) ) THEN
            j1 = J11(i)
            row = SDOF*PHDOF*(i-1) + r ; col = SDOF*PHDOF*(j1-1) + c
            H(row,col) = (tx * Gamma1(r,c))
            H(col,row) = CONJG(H(row,col))
          ENDIF
          ! hopping along y
          IF ( (J22(i) .GT. 0) .AND. (J22(i) .NE. i) ) THEN
            j2 = J22(i)
            row = SDOF*PHDOF*(i-1) + r ; col = SDOF*PHDOF*(j2-1) + c
            H(row,col) = (ty * Gamma1(r,c))
            H(col,row) = CONJG(H(row,col))
          ENDIF
          ! onsite terms, Chemical potential, s-Wave SC, exchange J_Hund
          row = SDOF*PHDOF*(i-1) + r ; col = SDOF*PHDOF*(i-1) + c
          H(row,col) = ( (mu * Gamma1(r,c)) + &
                         (Delta_s * Gamma2(r,c)) + &
                         (J_Hundx * Sx(ix, iy) * Gamma3(r,c) + &
                          J_Hundy * Sy(ix, iy) * Gamma4(r,c) + &
                          J_Hundz * Sz(ix, iy) * Gamma5(r,c)) )
          H(col,row) = CONJG(H(row,col))
        ENDDO
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE BuildHamiltonian_In_NumbuSpinor

SUBROUTINE BHZ_hamiltonian (H)
  USE input_parameters
  USE gen_coordinate
  IMPLICIT NONE

  COMPLEX, INTENT(OUT) :: H (n,n)
  
  INTEGER :: i, j1, j2, ix, iy, c, r, row, col
  INTEGER, DIMENSION(Lx*Ly) :: J11, J22
  COMPLEX, DIMENSION(SDOF, SDOF) :: sigma_0, sigma_x, sigma_y, sigma_z
  COMPLEX, DIMENSION(PHDOF, PHDOF) :: tau_0, tau_x, tau_y, tau_z
  COMPLEX, DIMENSION(SDOF*PHDOF,SDOF*PHDOF) :: Gamma1, Gamma2, Gamma3, Gamma4

  CALL NeighborArray (J11, J22)
  H(:,:) = CMPLX(0.0,0.0)
  ! Pauli matrices for spin (sigma) and particle-hole (tau) degrees of freedom
  sigma_x = PauliX ; tau_x = PauliX
  sigma_y = PauliY ; tau_y = PauliY
  sigma_z = PauliZ ; tau_z = PauliZ
  sigma_0 = Pauli0 ; tau_0 = Pauli0  
  ! Kronecker Product for PH and Spin space, [PH ⊗ Spin]
  CALL kronecker_product2(PHDOF, SDOF, tau_z, sigma_0, Gamma1)
  CALL kronecker_product2(PHDOF, SDOF, tau_x, sigma_z, Gamma2)
  CALL kronecker_product2(PHDOF, SDOF, tau_x, sigma_x, Gamma3)
  CALL kronecker_product2(PHDOF, SDOF, tau_y, sigma_0, Gamma4)

  ! The Hamiltonian is given in Ref. PRB 109, L041409 (2024)
  ! \begin{eqnarray}
  !   H & = & \sum_{i,j} c^{\dagger}_{i,j}[\{\mu \Gamma_1 + \Delta_0 \Gamma_2 + J_{\text{Hund}} (S_x \Gamma_3 + S_y \Gamma_4 + S_z \Gamma_5) \}c^{}_{i,j} \nonumber \\
  !                                         & & - t \Gamma_1c^{}_{i+1,j} - t \Gamma_1c^{}_{i,j+1} ] + H.c.
  ! \end{eqnarray}

  DO iy = 1, Ly
    DO ix = 1, Lx
      i = (iy-1) * Lx + ix
      DO c = 1, SDOF*PHDOF
        DO r = 1, SDOF*PHDOF
          ! hopping along x
          IF ( (J11(i) .GT. 0) .AND. (J11(i) .NE. i) ) THEN
            j1 = J11(i)
            row = SDOF*PHDOF*(i-1) + r ; col = SDOF*PHDOF*(j1-1) + c            
            H(row,col) = - 0.5 * tx * Gamma1(r,c) &
                         - iota * 0.5 * Lambda_SOC_x * Gamma2(r,c) &
                         + 0.5 * Lambda_WD_x * Gamma3(r,c)
            H(col,row) = CONJG(H(row,col))
          ENDIF
          ! hopping along y
          IF ( (J22(i) .GT. 0) .AND. (J22(i) .NE. i) ) THEN
            j2 = J22(i)
            row = SDOF*PHDOF*(i-1) + r ; col = SDOF*PHDOF*(j2-1) + c
            H(row,col) = - 0.5 * ty * Gamma1(r,c) &
                         - iota * 0.5 * Lambda_SOC_y * Gamma4(r,c) &
                         - 0.5 * Lambda_WD_y * Gamma3(r,c)
            H(col,row) = CONJG(H(row,col))
          ENDIF
          ! onsite terms, Chemical potential, s-Wave SC, exchange J_Hund
          row = SDOF*PHDOF*(i-1) + r ; col = SDOF*PHDOF*(i-1) + c
          H(row,col) = epsilon_0 * Gamma1(r,c)
          H(col,row) = CONJG(H(row,col))
        ENDDO
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE BHZ_hamiltonian


SUBROUTINE BHZ_Jex_sSC (Sx, Sy, Sz, H)
  USE input_parameters
  USE gen_coordinate
  IMPLICIT NONE

  REAL, DIMENSION(Lx,Ly), INTENT(IN) :: Sx, Sy, Sz
  COMPLEX, INTENT(OUT) :: H (n,n)
  
  INTEGER :: i, j1, j2, ix, iy, c, r, row, col
  INTEGER, DIMENSION(Lx*Ly) :: J11, J22
  COMPLEX, DIMENSION(PHDOF, PHDOF) :: tau_0, tau_x, tau_y, tau_z
  COMPLEX, DIMENSION(ODOF, ODOF) :: sigma_0, sigma_x, sigma_y, sigma_z
  COMPLEX, DIMENSION(SDOF, SDOF) :: s_0, s_x, s_y, s_z
  COMPLEX, DIMENSION(TDOF, TDOF) :: Gamma1, Gamma2, Gamma3, Gamma4, Gamma5, Gamma6, Gamma7, Gamma8

  CALL NeighborArray (J11, J22)
  H(:,:) = CMPLX(0.0,0.0)
  ! Pauli matrices for spin (sigma) and particle-hole (tau) degrees of freedom
  tau_x = PauliX ; sigma_x = PauliX ; s_x = PauliX
  tau_y = PauliY ; sigma_y = PauliY ; s_y = PauliY
  tau_z = PauliZ ; sigma_z = PauliZ ; s_z = PauliZ
  tau_0 = Pauli0 ; sigma_0 = Pauli0 ; s_0 = Pauli0
  ! Kronecker Product for PH and Spin space, [PH ⊗ Spin]
  CALL kronecker_product3(PHDOF, ODOF, SDOF, tau_z, sigma_z, s_0, Gamma1)
  CALL kronecker_product3(PHDOF, ODOF, SDOF, tau_x, sigma_0, s_0, Gamma2)
  CALL kronecker_product3(PHDOF, ODOF, SDOF, tau_0, sigma_0, s_x, Gamma3)
  CALL kronecker_product3(PHDOF, ODOF, SDOF, tau_0, sigma_0, s_y, Gamma4)
  CALL kronecker_product3(PHDOF, ODOF, SDOF, tau_0, sigma_0, s_z, Gamma5)
  CALL kronecker_product3(PHDOF, ODOF, SDOF, tau_z, sigma_x, s_z, Gamma6)
  CALL kronecker_product3(PHDOF, ODOF, SDOF, tau_z, sigma_y, s_0, Gamma7)
  CALL kronecker_product3(PHDOF, ODOF, SDOF, tau_z, sigma_x, s_x, Gamma8)

  ! The Hamiltonian is given in Ref. PRB 109, L041409 (2024)
  ! \begin{eqnarray}
  !   H & = & \sum_{i,j} c^{\dagger}_{i,j}[\{\mu \Gamma_1 + \Delta_0 \Gamma_2 + J_{\text{Hund}} (S_x \Gamma_3 + S_y \Gamma_4 + S_z \Gamma_5) \}c^{}_{i,j} \nonumber \\
  !                                         & & - t \Gamma_1c^{}_{i+1,j} - t \Gamma_1c^{}_{i,j+1} ] + H.c.
  ! \end{eqnarray}

  DO iy = 1, Ly
    DO ix = 1, Lx
      i = (iy-1) * Lx + ix
      DO c = 1, TDOF
        DO r = 1, TDOF
          ! hopping along x
          IF ( (J11(i) .GT. 0) .AND. (J11(i) .NE. i) ) THEN
            j1 = J11(i)
            row = TDOF*(i-1) + r ; col = TDOF*(j1-1) + c            
            H(row,col) = - tx * Gamma1(r,c) &
                         - iota * Lambda_SOC_x * Gamma6(r,c) &
                         + Lambda_WD_x * Gamma8(r,c)
            H(col,row) = CONJG(H(row,col))
          ENDIF
          ! hopping along y
          IF ( (J22(i) .GT. 0) .AND. (J22(i) .NE. i) ) THEN
            j2 = J22(i)
            row = TDOF*(i-1) + r ; col = TDOF*(j2-1) + c
            H(row,col) = - ty * Gamma1(r,c) &
                         - iota * Lambda_SOC_y * Gamma7(r,c) &
                         - Lambda_WD_y * Gamma8(r,c)
            H(col,row) = CONJG(H(row,col))
          ENDIF
          ! onsite terms, Chemical potential, s-Wave SC, exchange J_Hund
          row = TDOF*(i-1) + r ; col = TDOF*(i-1) + c
          H(row,col) = epsilon_0 * Gamma1(r,c) + &
                       Delta_s * Gamma2(r,c) + &
                       (J_Hundx * Sx(ix, iy) * Gamma3(r,c) + &
                        J_Hundy * Sy(ix, iy) * Gamma4(r,c) + &
                        J_Hundz * Sz(ix, iy) * Gamma5(r,c))
          H(col,row) = CONJG(H(row,col))
        ENDDO
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE BHZ_Jex_sSC


