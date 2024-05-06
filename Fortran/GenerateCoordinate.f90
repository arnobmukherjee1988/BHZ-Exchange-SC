MODULE gen_coordinate
  USE input_parameters
  IMPLICIT NONE

  CONTAINS
  SUBROUTINE NeighborArray(J1_arr, J2_arr)
  USE input_parameters
  IMPLICIT NONE

  INTEGER, DIMENSION(Lx*Ly), INTENT(OUT) :: J1_arr, J2_arr
  INTEGER :: ix, iy, jx, jy, i

  DO iy = 1, ly
    DO ix = 1, lx
      i = (iy-1)*Lx + ix
      
      SELECT CASE (BoundaryCondition)
      CASE ('PBCx + PBCy')  ! Periodic boundary conditions along x and y
        jx = ix + 1 ; IF (jx .GT. lx) jx = 1
        jy = iy + 1 ; IF (jy .GT. ly) jy = 1
        J1_arr(i) =  (iy-1)*Lx + jx
        J2_arr(i) =  (jy-1)*Lx + ix
      CASE ('OBCx + PBCy')  ! Open boundary condition along x and periodic along y
        jx = ix + 1 ; IF (jx .GT. lx) jx = 1
        jy = iy + 1 ; IF (jy .GT. ly) jy = 1
        J1_arr(i) =  (iy-1)*Lx + jx
        J2_arr(i) =  (jy-1)*Lx + ix
        IF (ix .EQ. lx) THEN
          J1_arr(i) =  -i
        ENDIF
      CASE ('PBCx + OBCy')  ! Periodic boundary condition along x and open along y
        jx = ix + 1 ; IF (jx .GT. lx) jx = 1
        jy = iy + 1 ; IF (jy .GT. ly) jy = 1
        J1_arr(i) =  (iy-1)*Lx + jx
        J2_arr(i) =  (jy-1)*Lx + ix
        IF (iy .EQ. ly) THEN
          J2_arr(i) =  -i
        ENDIF
      CASE ('OBCx + OBCy')  ! Open boundary conditions along x and y
        jx = ix + 1 ; IF (jx .GT. lx) jx = 1
        jy = iy + 1 ; IF (jy .GT. ly) jy = 1
        J1_arr(i) =  (iy-1)*Lx + jx
        J2_arr(i) =  (jy-1)*Lx + ix
        IF (ix .EQ. lx) THEN
          J1_arr(i) =  -i
        ENDIF
        IF (iy .EQ. ly) THEN
          J2_arr(i) =  -i
        ENDIF
      CASE DEFAULT
        PRINT *, BC_Error
        STOP
      END SELECT

    ENDDO
  ENDDO

  END SUBROUTINE NeighborArray

  SUBROUTINE GetSiteCoordinates(i, ix, iy)
  USE input_parameters
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i
  INTEGER, INTENT(OUT) :: ix, iy

  iy = (i - 1) / lx + 1
  ix = MOD(i - 1, lx) + 1
  END SUBROUTINE GetSiteCoordinates

END MODULE gen_coordinate


