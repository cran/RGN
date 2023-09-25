SUBROUTINE Qr_f90(nRow,nCol,a,q,r)
    ! Modified Gram-Schmidt process

   use constantsMod, only: rk

    IMPLICIT NONE

    INTERFACE

    FUNCTION Norm(x)
    USE constantsMod, ONLY: rk
    IMPLICIT NONE
    REAL(rk), DIMENSION(:), INTENT(in) :: x
    REAL(rk) :: Norm
    END FUNCTION Norm

    END INTERFACE

    integer, intent(in) :: nRow, nCol
    REAL(rk), DIMENSION(nRow,nCol), INTENT(in) :: a
    REAL(rk), DIMENSION(nRow,nCol), INTENT(out) :: q
    REAL(rk), DIMENSION(nCol,nCol), INTENT(out) :: r
    REAL(rk), DIMENSION(nRow,nCol) :: a0
    integer :: k,n
    n = nCol
    q(:,:) = 0.0_rk
    r(:,:) = 0.0_rk
    a0(:,:) = a(:,:)

    DO k = 1, n
       r(k,k) = Norm(a0(:,k))
!       q(:,k) = a0(:,k) / r(k,k)
       where(a0(:,k)/=0._rk)  !FIXME_DK: rough fix will need checks+refinement
       q(:,k) = a0(:,k) / r(k,k)
       elsewhere
        q(:,k) = 0._rk
       endwhere
       r(k,k+1:n) = MATMUL(q(:,k), a0(:,k+1:n))
       a0(:,k+1:n) = a0(:,k+1:n) - MATMUL(q(:,k:k), r(k:k,k+1:n))
    END DO

END SUBROUTINE Qr_f90

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 PURE FUNCTION Norm(x) RESULT(valu)
    ! L2-norm
    USE constantsMod, ONLY: rk
    IMPLICIT NONE
    REAL(rk), DIMENSION(:), INTENT(in) :: x
    REAL(rk) :: valu
    valu = SQRT(SUM(x**2))
  END FUNCTION Norm


