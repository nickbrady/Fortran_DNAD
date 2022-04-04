! MATINV computes 𝐁⁻¹𝐃 using Gaussian elimination.
! The matrices 𝐁 and 𝐃 are received from subroutine BAND(J)
! 𝐁 is transformed into the identity matrix 𝐈 using Guassian elimination
! The row operations that are performed on 𝐁 are simultaneously performed on 𝐃
! with the end result being that 𝐃 is transformed into 𝐁⁻¹𝐃 (i.e. 𝐃 <-- 𝐁⁻¹𝐃)

SUBROUTINE MATINV(N, M, DETERM)
  use variables, only: B, D, ID ! A imported but not used

  implicit none
  integer, intent(in) :: N, M
  ! integer :: NN, I, J, JC, jcol
  ! real    :: BMAX, BNEXT, BTRY,
  real    :: F, DETERM

  integer :: i_max, irow, icol, ii

  DETERM = 1.0

  ! The opposite of >= is <, but strange behavior can happen when encountering NaN's
  !(.NOT.(BNEXT >= BMAX*BTRY))

  do irow = 1,N
        icol = irow

        i_max = MAXLOC( ABS(B(irow:,icol)), 1 )
        i_max = (i_max - 1) + irow

        if (B(i_max,icol) == 0) then
            !/* No pivot in this column, pass to next column */
            if (icol == N) then
                DETERM = 0
            end if
            cycle
        end if

        ! else:
        ! swap rows(irow, i_max)
        if (i_max /= irow) then
            call swap(B(IROW,:), B(i_max,:))         ! row swap
            call swap(D(IROW,:), D(i_max,:))
        end if

        ! make the leading values 1
        F = B(IROW, ICOL)
        B(IROW,:)  = B(IROW,:)  / F
        D(IROW,:M) = D(IROW,:M) / F

        do ii = 1, N
            if (ii /= irow) then
                F = B(ii,icol)
                B(ii,:)   = B(ii,:)  - F * B(irow,:)
                D(ii,:M)  = D(ii,:M) - F * D(irow,:M)
            end if
        end do

  end do  ! NN = 1,N

  return

contains
  elemental subroutine swap(a, b)
     ! row exchange: switch row A and row B
     real, intent(in out) :: a, b
     real :: save_val

     save_val = a
     a = b
     b = save_val
  end subroutine swap

end subroutine MATINV
