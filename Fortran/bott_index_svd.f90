SUBROUTINE bott_index_svd (Z, bott_indx)
  USE input_parameters
  IMPLICIT NONE
  
  COMPLEX, DIMENSION(n,n),INTENT(IN) :: Z
  REAL, INTENT(OUT) :: bott_indx
  
  INTEGER :: i,j,k,ix,iy,spindex,j1,j2,ii,jj, Shift
  COMPLEX, DIMENSION(n,n) :: Prj, UMatrix, VMatrix, Udmatrix, Vdmatrix, Bott_Matrix
  COMPLEX, DIMENSION(n,n) :: UTildaMatrix, VTildaMatrix, UTildaMatrix_d, VTildaMatrix_d
  COMPLEX, DIMENSION(2*n) :: Wb
  COMPLEX :: B_ind, Det

  complex, dimension(n, n) :: Eth, Eph
  
  ! variable declarations for use in CGEES (LAPACK subroutine)
  CHARACTER*1 :: JOBVS,SORT
  INTEGER :: LWORK,INFO,SDIM,LDVS
  INTEGER :: IFAIL(n/2),LDZ,label,LDA
  COMPLEX :: WORK(n), VS(1,n)
  REAL :: RWORK(n)
  LOGICAL :: BWORK(n), SELECT
  
  !---------- SVD subroutine variables CGESVD ---------------------
  INTEGER, PARAMETER :: LDA_svd = n, LDU_svd = n, LDVT_svd = n
  INTEGER, PARAMETER :: LWORK_svd = 5*n
  INTEGER :: INFO_svd
  REAL :: S_svd(n), RWORK_svd(5*n)
  COMPLEX :: U_svd(n,n), VT_svd(n,n), WORK_svd(LWORK_svd)  
  !----------------------------------------------------------------
  
  !  For Bott Index
  Shift = SDOF * Nsites
  DO spindex = 1, SDOF
    DO iy = 1, Ly
      DO ix = 1, Lx
        i = (iy-1)*Lx*SDOF + (ix-1)*SDOF + spindex
        Eth (i,i) = cexp(cmplx(0.,1.)*2.*Pi/real(Lx)*real(ix-1))
        Eph (i,i) = cexp(cmplx(0.,1.)*2.*Pi/real(Lx)*real(ix-1))
        Eth (i+Shift,i+Shift) = cexp(cmplx(0.,1.)*2.*Pi/real(Lx)*real(ix-1))
        Eph (i+Shift,i+Shift) = cexp(cmplx(0.,1.)*2.*Pi/real(Lx)*real(ix-1))
      ENDDO
    ENDDO
  ENDDO

  !Define Projector
  Prj(:,:) = CMPLX(0.0,0.0)
  DO k = 1, n
    DO i = 1, n
      DO j = 1, n
        Prj(i,j) = Prj(i,j) + Z(i,k) * conjg(Z(j,k))
      ENDDO
    ENDDO
  ENDDO

  ! ! Sum over all eigenstates
  ! DO k = 1, Nel
  !   ! Get the k-th eigenvector
  !   DO i = 1, n
  !     psi_k(i) = Z(i, k)
  !   END DO
      
  !   ! Calculate the outer product and add it to the projection operator
  !   DO i = 1, n
  !     psi_k_i = psi_k(i)
  !     DO j = 1, N
  !       psi_k_conj_j = CONJG(Z(j, k))
  !       Prj(i, j) = Prj(i, j) + psi_k_i * psi_k_conj_j
  !     END DO
  !   END DO
  ! END DO
  
  UMatrix = MATMUL( Prj,MATMUL(Eth,Prj) )
  VMatrix = MATMUL( Prj,MATMUL(Eph,Prj) )
  
  DO ii = 1,N
    DO jj = 1,N
      IF (ii.EQ.jj) THEN
        UMatrix(ii,jj) = UMatrix(ii,jj) + CMPLX(1.0,0.0) - Prj(ii,jj)
        VMatrix(ii,jj) = VMatrix(ii,jj) + CMPLX(1.0,0.0) - Prj(ii,jj)
      ELSE
        UMatrix(ii,jj) = UMatrix(ii,jj) + CMPLX(0.0,0.0) - Prj(ii,jj)
        VMatrix(ii,jj) = VMatrix(ii,jj) + CMPLX(0.0,0.0) - Prj(ii,jj)
      ENDIF      
    ENDDO
  ENDDO
  
  !******************** SVD of UMatrix and VMatrix *******************************************************
  
  ! Compute SVD for UMatrix
  CALL CGESVD( 'All', 'All', n, n, UMatrix, LDA_svd, S_svd, U_svd, LDU_svd, VT_svd, LDVT_svd, WORK_svd, LWORK_svd, RWORK_svd, &
                                                                                                                     & INFO_svd )
  
  ! Check for convergence
  IF( INFO_svd .GT. 0 ) THEN
      WRITE(*,*)'The algorithm computing SVD failed to converge for UMatrix.'
      STOP
  END IF
  
  UTildaMatrix = MATMUL (U_svd,VT_svd)

  ! Compute SVD
  CALL CGESVD( 'All', 'All', n, n, VMatrix, LDA_svd, S_svd, U_svd, LDU_svd, VT_svd, LDVT_svd, WORK_svd, LWORK_svd, RWORK_svd, &
                                                                                                                     & INFO_svd )
  
  ! Check for convergence
  IF( INFO_svd .GT. 0 ) THEN
      WRITE(*,*)'The algorithm computing SVD failed to converge for VMatrix.'
      STOP
  END IF
  
  VTildaMatrix = MATMUL (U_svd,VT_svd)
  
  !******************** SVD of UMatrix and VMatrix *******************************************************
  
    DO ii = 1,N
      DO jj = 1,N
        UTildaMatrix_d(ii,jj) = CONJG(UTildaMatrix(jj,ii))
        VTildaMatrix_d(ii,jj) = CONJG(VTildaMatrix(jj,ii))
      ENDDO
    ENDDO   
  
  Bott_Matrix = MATMUL( MATMUL(VTildaMatrix,UTildaMatrix),MATMUL(VTildaMatrix_d,UTildaMatrix_d) )
  
  JOBVS='N' ; SORT = 'N' ; LDVS = 1 ;  LWORK = 2*N
  LDA = N  
  LDZ = N
  
  CALL CGEES (JOBVS, SORT, SELECT, N, Bott_Matrix, LDA, SDIM, Wb, VS, LDVS, WORK, LWORK, RWORK, BWORK, INFO)

  B_ind = CMPLX(0.0,0.0)
  Det = CMPLX(0.0,0.0)
  
  DO i = 1,N
    B_ind = B_ind + CLOG(Wb(i))
    Det = Det*Wb(i)
  ENDDO
  
  bott_indx = AIMAG(B_ind)/(2.*Pi)

END SUBROUTINE bott_index_svd
