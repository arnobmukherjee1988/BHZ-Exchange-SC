FUNCTION TRACE (MatrixSize, Matrix)
  USE input_parameters
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: MatrixSize
  COMPLEX, INTENT(IN) :: Matrix(MatrixSize,MatrixSize)
  COMPLEX :: TRACE
  INTEGER :: i

  TRACE = CMPLX (0.0, 0.0)
  DO i = 1, MatrixSize
    TRACE = TRACE + Matrix(i, i)
  END DO
END FUNCTION TRACE

SUBROUTINE print_eigenvals (NumEV, W, filename)
  USE input_parameters
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NumEV
  REAL, INTENT(IN) :: W(NumEV)
  CHARACTER(LEN=*), INTENT(IN) :: filename
  INTEGER :: i

  OPEN (UNIT = 1, FILE = TRIM(ADJUSTL(filename)), STATUS = 'unknown')
  DO i = 1, NumEV
    WRITE(1,*) i, W(i)
  ENDDO
END SUBROUTINE print_eigenvals
  

SUBROUTINE print_SpinVectors (NumRow, NumCol, Sx, Sy, Sz, Theta, Phi, filename)
  USE input_parameters
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NumRow, NumCol
  REAL, DIMENSION (NumRow, NumCol), INTENT(IN) :: Sx, Sy, Sz, Theta, Phi
  CHARACTER(LEN=*), INTENT(IN) :: filename
  INTEGER :: ix, iy, i

  OPEN (UNIT = 100, FILE = TRIM(ADJUSTL(filename)), STATUS = 'unknown')
  DO iy = 1, NumCol
    DO ix = 1, NumRow
      i = (iy-1)*NumRow + ix
      WRITE (100, *) ix, iy, i, Sx(ix,iy), Sy(ix,iy), Sz(ix,iy), Theta(ix,iy), Phi(ix,iy)
    ENDDO
  ENDDO
END SUBROUTINE print_SpinVectors

SUBROUTINE diagonalization (size, H, CalEV, W, Z)
  USE input_parameters
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: size
  COMPLEX, INTENT(IN) :: H (size, size)
  CHARACTER (LEN=*), INTENT(IN) :: CalEV
  COMPLEX, INTENT(OUT) :: Z (size, size)
  REAL, INTENT(OUT) :: W(size)
  
  INTEGER :: IFAIL(size), IWORK(5*size)
  REAL :: RWORK(7*size)
  COMPLEX, DIMENSION (2*size) :: WORK
  INTEGER :: IL, INFO, IU, LDA, LDZ, LWORK, M
  REAL :: ABSTOL, SLAMCH, VL, VU

  LWORK = 2*size
  LDA = size
  LDZ = size
  ABSTOL = 2*SLAMCH('S') + 10.0**(-8.0)

  ! diagonalization subroutine
  IF (CalEV .EQ. 'Ev_Yes') THEN
    CALL cheevx('V','A','U',size,H,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)
  ELSEIF (CalEV .EQ. 'Ev_No') THEN
    CALL cheevx('N','A','U',size,H,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)
  ELSE
    print*, "Diagonalization error: Put 'Ev_Yes' or 'Ev_No' for calculating or not of eigenvectors"
  ENDIF
  END SUBROUTINE diagonalization

  FUNCTION fermi (x)
  USE input_parameters
  IMPLICIT NONE
  REAL, INTENT(IN) :: x
  REAL :: fermi, width_parameter
  width_parameter = 8.0
  IF (ABS(x) .LE. width_parameter) THEN
    fermi = ( 1.0 + EXP(x) ) ** (-1.0)
  ELSEIF (x .GT. width_parameter) THEN
    fermi = 0.0
  ELSE
    fermi = 1.0
  ENDIF
END FUNCTION fermi

REAL FUNCTION Distance(x1, y1, x2, y2)
  USE input_parameters
  IMPLICIT NONE

  REAL, INTENT(IN) :: x1, y1, x2, y2
  REAL :: dx, dy, dist

  ! Use KIND of x1 for calculations to ensure consistent precision
  dist = REAL(0.0, kind=KIND(x1))

  ! Convert integer arguments to real, if necessary
  dx = REAL(x2 - x1, kind=KIND(x1))
  dy = REAL(y2 - y1, kind=KIND(x1))

  ! Calculate squared distance
  dist = dx*dx + dy*dy

  ! Square root the distance, handling potential negative inputs
  IF (dist .gt. 0.0) THEN
    dist = SQRT(dist)
  ELSE
    dist = 0.0
  ENDIF

  distance = dist
END FUNCTION Distance


! Kronecker delta function.
FUNCTION delta(i, j) result(result)
  USE input_parameters
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i, j
  INTEGER :: result
  IF (i == j) THEN
    result = 1
  ELSE
    result = 0
  ENDIF
ENDFUNCTION delta


FUNCTION Lorentzian (BroadFac, x) 
  USE input_parameters
  REAL, INTENT(IN) :: BroadFac, x
  REAL :: Lorentzian

  Lorentzian = (1.0/pi) * (BroadFac/( BroadFac**2 + x**2 ))
END FUNCTION Lorentzian

! Kronecker tensor product
SUBROUTINE kronecker_product2(sizeA, sizeB, matrixA, matrixB, result)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: sizeA, sizeB
  COMPLEX, INTENT(IN) :: matrixA(sizeA, sizeA), matrixB(sizeB, sizeB)
  COMPLEX, INTENT(OUT) :: result(sizeA*sizeB, sizeA*sizeB)
  INTEGER :: i, j, k, l

  ! Perform Kronecker product
  DO i = 1, sizeA
    DO j = 1, sizeA
      DO k = 1, sizeB
        DO l = 1, sizeB
          result((i-1)*sizeB + k, (j-1)*sizeB + l) = matrixA(i, j) * matrixB(k, l)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
end subroutine kronecker_product2

! Kronecker tensor product for three matrices
SUBROUTINE kronecker_product3(sizeA, sizeB, sizeC, matrixA, matrixB, matrixC, result)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: sizeA, sizeB, sizeC
  COMPLEX, INTENT(IN) :: matrixA(sizeA, sizeA), matrixB(sizeB, sizeB), matrixC(sizeC, sizeC)
  COMPLEX, INTENT(OUT) :: result(sizeA*sizeB*sizeC, sizeA*sizeB*sizeC)
  INTEGER :: i, j, k, l, m, n

  ! Perform Kronecker product
  DO i = 1, sizeA
    DO j = 1, sizeA
      DO k = 1, sizeB
        DO l = 1, sizeB
          DO m = 1, sizeC
            DO n = 1, sizeC
              result((i-1)*sizeB*sizeC + (k-1)*sizeC + m, &
                     (j-1)*sizeB*sizeC + (l-1)*sizeC + n) = &
                     matrixA(i, j) * matrixB(k, l) * matrixC(m, n)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE kronecker_product3

subroutine compute_exp_matrix(size, A, expA)
  implicit none
  
  ! Inputs
  integer, intent(in) :: size
  complex, dimension(size,size), intent(in) :: A
  
  ! Outputs
  complex, dimension(size,size), intent(out) :: expA
  
  ! Declarations
  integer :: info ! Error flag
  integer :: lda ! Leading dimension of array A
  integer :: lwork ! Size of workspace array
  complex, allocatable :: work(:) ! Workspace array
  integer, allocatable :: ipiv(:) ! Pivot indices array
  integer :: i, j
  
  ! Perform LU factorization
  lda = size
  allocate(ipiv(size), work(1))
  call zgetrf(size, size, A, lda, ipiv, info)
  if (info /= 0) then
      print *, "Error: zgetrf failed with error code", info
      stop
  end if
  
  ! Compute matrix exponential
  lwork = size * size * 4
  allocate(work(lwork))
  call zgetri(size, A, lda, ipiv, work, lwork, info)
  if (info /= 0) then
      print *, "Error: zgetri failed with error code", info
      stop
  end if
  
  ! Store the result in expA
  expA = A
  
  deallocate(ipiv, work)
    
end subroutine compute_exp_matrix


