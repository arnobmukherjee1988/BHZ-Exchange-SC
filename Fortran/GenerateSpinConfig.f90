SUBROUTINE DummySkyrmion (beta, Sx, Sy, Sz, Theta, Phi)
  USE input_parameters
  IMPLICIT NONE

  REAL, INTENT (IN) :: beta
  REAL, DIMENSION (lx,ly), INTENT(OUT) :: Sx, Sy, Sz, Theta, Phi

  REAL, DIMENSION(INT(Lx/NumSkAlongX), INT(Ly/NumSkAlongY)) :: theta_c, phi_c
  INTEGER :: ix_o, iy_o, cx, cy
  INTEGER :: ix, iy, i, j, ii, jj, i1, j1
  REAL :: r, num, den

  theta_c(:,:) = 0.0 ; phi_c(:,:) = 0.0
  sx(:,:) = 0.0 ; sy(:,:) = 0.0 ; sz(:,:) = 0.0
  Theta(:,:) = 0.0 ; Phi(:,:) = 0.0
  
  cx = INT(Lx/NumSkAlongX) ; cy = INT(Ly/NumSkAlongY)
  ix_o = cx/2 ; iy_o = cy/2
  DO iy = 1,cy
    DO ix = 1,cx
      i = (iy-1)*cx + ix
      r = SQRT((REAL(ix-ix_o))**2 + (REAL(iy-iy_o))**2)
      phi_c(ix,iy) =  ACOS(vorticity*(ix-ix_o)/r) - helicity
      IF ((iy-iy_o) .LT. 0) phi_c(ix,iy) = - helicity - ACOS(vorticity*(ix-ix_o)/r)
      IF (((ix-ix_o) .EQ. 0) .AND. ((iy-iy_o) .EQ. 0)) phi_c(ix,iy) = 0.0 + helicity
      num = SkyrRadius * EXP(beta*(SkyrRadius-r)) ; den = r
      theta_c(ix,iy) = 2.0 * ATAN(num/den) + ACOS(-polarity)
      IF((ix.EQ.1 .OR. ix.EQ.cx) .AND. (iy.EQ.1 .OR. iy.EQ.cy)) theta_c(ix,iy) = 0.0 + ACOS(polarity)
    ENDDO
  ENDDO

  DO j = 1, Ly/cy
    DO i = 1, Lx/cx
      DO j1 = 1, cy
        DO i1 = 1, cx
          ii = ((i-1) * cx) + i1
          jj = ((j-1) * cy) + j1
          Phi(ii,jj) = phi_c(i1, j1)
          Theta(ii,jj) = theta_c(i1, j1)

          Sx(ii,jj) = spin*SIN(Theta(ii,jj))*COS(Phi(ii,jj))
          Sy(ii,jj) = spin*SIN(Theta(ii,jj))*SIN(Phi(ii,jj))
          Sz(ii,jj) = spin*COS(Theta(ii,jj))
        ENDDO
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE DummySkyrmion


SUBROUTINE calculate_skyrmion_number(Sx, Sy, Sz, skyrmion_number)
  USE input_parameters
  IMPLICIT NONE
  
  REAL, DIMENSION(Lx, Ly), INTENT(IN) :: Sx, Sy, Sz
  REAL :: dx, dy
  REAL :: skyrmion_number
  INTEGER :: i, j
  REAL :: dSx_dx, dSx_dy, dSy_dx, dSy_dy, dSz_dx, dSz_dy
  REAL :: cross_product(3), dot_product
  
  dx = 1.0 / Lx
  dy = 1.0 / Ly
  skyrmion_number = 0.0
  DO i = 2, Lx
    DO j = 2, Ly
      ! Calculate derivatives using finite differences
      dSx_dx = (Sx(MOD(i+1, Lx), j) - Sx(MOD(i-1, Lx), j)) / (2.0 * dx)
      dSx_dy = (Sx(i, MOD(j+1, Ly)) - Sx(i, MOD(j-1, Ly))) / (2.0 * dy)
      
      dSy_dx = (Sy(MOD(i+1, Lx), j) - Sy(MOD(i-1, Lx), j)) / (2.0 * dx)
      dSy_dy = (Sy(i, MOD(j+1, Ly)) - Sy(i, MOD(j-1, Ly))) / (2.0 * dy)
      
      dSz_dx = (Sz(MOD(i+1, Lx), j) - Sz(MOD(i-1, Lx), j)) / (2.0 * dx)
      dSz_dy = (Sz(i, MOD(j+1, Ly)) - Sz(i, MOD(j-1, Ly))) / (2.0 * dy)
      
      ! Calculate cross product
      cross_product(1) = dSy_dx * dSz_dy - dSz_dx * dSy_dy
      cross_product(2) = dSz_dx * dSx_dy - dSx_dx * dSz_dy
      cross_product(3) = dSx_dx * dSy_dy - dSy_dx * dSx_dy
      
      ! Calculate dot product
      dot_product = Sx(i, j) * cross_product(1) + Sy(i, j) * cross_product(2) + Sz(i, j) * cross_product(3)
      
      ! Summing up Skyrmion number
      skyrmion_number = skyrmion_number + dot_product
    ENDDO
  ENDDO
  
  ! Normalize Skyrmion number
  skyrmion_number = skyrmion_number / (4.0 * PI * REAL(Lx) * REAL(Ly))

END SUBROUTINE calculate_skyrmion_number


SUBROUTINE Skyrmion (Sx, Sy, Sz, Theta, Phi)
  USE input_parameters
  IMPLICIT NONE
  
  REAL, DIMENSION(Lx,Ly), INTENT(OUT) :: Sx, Sy, Sz, Theta, Phi
  REAL :: beta_value, prev_skyr_number, skyr_number, beta_critical
  INTEGER :: i_beta
  
  prev_skyr_number = 0.0
  beta_critical = 0.0
  DO i_beta = 1, no_beta
    beta_value = beta_min + ((beta_max - beta_min) / real(no_beta)) * real(i_beta)
    CALL DummySkyrmion(beta_value, Sx, Sy, Sz, Theta, Phi)
    CALL calculate_skyrmion_number(Sx, Sy, Sz, skyr_number)
    IF ((i_Beta .GT. 1) .AND. (ABS(skyr_number) <= ABS(prev_skyr_number))) THEN
      beta_critical = beta_min + ((beta_max - beta_min) / real(no_beta)) * real(i_beta-1)
      exit
    ENDIF
    prev_skyr_number = skyr_number
  ENDDO
  IF ( ABS(beta_critical) > epsilon(1.0) ) THEN
    CALL DummySkyrmion(beta_critical, Sx, Sy, Sz, Theta, Phi)
  ELSE
    ! Handle the case when beta_critical is not found
    Sx(:,:) = 0.0
    Sy(:,:) = 0.0
    Sz(:,:) = 0.0
    Theta(:,:) = 0.0
    Phi(:,:) = 0.0
  ENDIF
END SUBROUTINE Skyrmion


SUBROUTINE Spiral (Sx, Sy, Sz, Theta, Phi)
  USE input_parameters
  IMPLICIT NONE

  REAL, DIMENSION (Lx,Ly), INTENT(OUT) :: Sx, Sy, Sz, Theta, Phi
  INTEGER :: ix, iy, i

  Sx(:,:) = 0.0 ; Sy(:,:) = 0.0 ; Sz(:,:) = 0.0
  Theta(:,:) = 0.0 ; Phi(:,:) = 0.0
  DO iy = 1, Ly
    DO ix = 1, Lx
      i = (iy-1)*Lx + ix
      Theta(ix, iy) = qx*ix + qy*iy ; Phi(ix, iy) = qx*ix + qy*iy
      Sx(ix, iy) = spin*SIN(Theta(ix, iy))*COS(Phi(ix, iy))
      Sy(ix, iy) = spin*SIN(Theta(ix, iy))*SIN(Phi(ix, iy))
      Sz(ix, iy) = spin*COS(Theta(ix, iy))
    ENDDO
  ENDDO
END SUBROUTINE Spiral


SUBROUTINE PRB_109_L041409 (Sx, Sy, Sz, Theta, Phi)
  USE input_parameters
  IMPLICIT NONE

  REAL, DIMENSION (Lx,Ly), INTENT(OUT) :: Sx, Sy, Sz, Theta, Phi
  INTEGER :: ix, iy, i

  Sx(:,:) = 0.0 ; Sy(:,:) = 0.0 ; Sz(:,:) = 0.0
  Theta(:,:) = 0.0 ; Phi(:,:) = 0.0
  DO iy = 1, Ly
    DO ix = 1, Lx
      i = (iy-1)*Lx + ix
      Theta(ix, iy) = PI/2.0 ; Phi(ix, iy) = qx*ix + qy*iy
      Sx(ix, iy) = spin*SIN(Theta(ix, iy))*COS(Phi(ix, iy))
      Sy(ix, iy) = spin*SIN(Theta(ix, iy))*SIN(Phi(ix, iy))
      Sz(ix, iy) = spin*COS(Theta(ix, iy))
    ENDDO
  ENDDO
END SUBROUTINE PRB_109_L041409


SUBROUTINE AKG_code (Sx, Sy, Sz, Theta, Phi)
  USE input_parameters
  IMPLICIT NONE

  REAL, DIMENSION (Lx,Ly), INTENT(OUT) :: Sx, Sy, Sz, Theta, Phi
  INTEGER :: ix, iy, i

  Sx(:,:) = 0.0 ; Sy(:,:) = 0.0 ; Sz(:,:) = 0.0
  Theta(:,:) = 0.0 ; Phi(:,:) = 0.0
  DO iy = 1, Ly
    DO ix = 1, Lx
      i = (iy-1)*Lx + ix
      Theta(ix, iy) = PI/2.0 ; Phi(ix, iy) = qx*ix + qy*iy
      Sx(ix, iy) = spin*SIN(Theta(ix, iy))*COS(Phi(ix, iy))
      Sy(ix, iy) = spin*SIN(Theta(ix, iy))*COS(Phi(ix, iy))
      Sz(ix, iy) = spin*COS(Theta(ix, iy))
    ENDDO
  ENDDO
END SUBROUTINE AKG_code


SUBROUTINE Ferro (Sx, Sy, Sz, Theta, Phi)
  USE input_parameters
  IMPLICIT NONE

  REAL, DIMENSION (Lx,Ly), INTENT(OUT) :: sx, sy, sz, theta, phi
  INTEGER :: ix, iy, i

  Sx(:,:) = 0.0 ; Sy(:,:) = 0.0 ; Sz(:,:) = 0.0
  Theta(:,:) = 0.0; Phi(:,:) = 0.0
  DO iy = 1, Ly
    DO ix = 1, Lx
      i = (iy-1)*Lx + ix
      Theta(ix, iy) = 0.0 ; Phi(ix, iy) = qx*ix + qy*iy
      Sx(ix, iy) = spin*SIN(Theta(ix, iy))*COS(Phi(ix, iy))
      Sy(ix, iy) = spin*SIN(Theta(ix, iy))*SIN(Phi(ix, iy))
      Sz(ix, iy) = spin*COS(Theta(ix, iy))
    ENDDO
  ENDDO

END SUBROUTINE Ferro


SUBROUTINE AntiFerro (Sx, Sy, Sz, Theta, Phi)
  USE input_parameters
  IMPLICIT NONE

  REAL, DIMENSION (Lx,Ly), INTENT(OUT) :: Sx, Sy, Sz, Theta, Phi
  INTEGER :: ix, iy, i

  Sx(:,:) = 0.0 ; Sy(:,:) = 0.0 ; Sz(:,:) = 1.0
  Theta(:,:) = 0.0; Phi(:,:) = 0.0
  DO iy = 1, Ly
    DO ix = 1, Lx
      i = (iy-1)*Lx + ix
      Phi(ix, iy) = 0.0
      IF (MOD(ix+iy, 2) .EQ. 0) THEN
        Theta(ix, iy) = 0.0
      ELSE
        Theta(ix, iy) = -PI
      ENDIF
      Sx(ix, iy) = spin*SIN(Theta(ix, iy))*COS(Phi(ix, iy))
      Sy(ix, iy) = spin*SIN(Theta(ix, iy))*SIN(Phi(ix, iy))
      Sz(ix, iy) = spin*COS(Theta(ix, iy))
    ENDDO
  ENDDO
END SUBROUTINE AntiFerro

SUBROUTINE BuildSpinConfig (Sx, Sy, Sz, Theta, Phi)
  USE input_parameters
  IMPLICIT NONE

  REAL, DIMENSION (Lx,Ly), INTENT(OUT) :: Sx, Sy, Sz, Theta, Phi

  SELECT CASE (SpinConfig)
    CASE ('Skyrmion')  ! Periodic boundary conditions along x and y
      CALL Skyrmion (Sx, Sy, Sz, Theta, Phi)
    CASE ('Spiral')  ! Open boundary condition along x and periodic along y
      CALL Spiral (Sx, Sy, Sz, Theta, Phi)
    CASE ('PRB_109_L041409')  ! Periodic boundary condition along x and open along y
      CALL PRB_109_L041409 (Sx, Sy, Sz, Theta, Phi)
    CASE ('AKG_code')  ! Open boundary conditions along x and y
      CALL AKG_code (Sx, Sy, Sz, Theta, Phi)
    CASE ('Ferro')  ! Open boundary conditions along x and y
      CALL Ferro (Sx, Sy, Sz, Theta, Phi)
    CASE ('AntiFerro')  ! Open boundary conditions along x and y
      CALL AntiFerro (Sx, Sy, Sz, Theta, Phi)
    CASE DEFAULT
      PRINT *, BC_Error
      STOP
  END SELECT

END SUBROUTINE BuildSpinConfig


