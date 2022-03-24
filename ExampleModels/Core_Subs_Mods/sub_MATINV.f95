! MATINV computes 𝐁⁻¹𝐃 using Gaussian elimination.
! The matrices 𝐁 and 𝐃 are received from subroutine BAND(J)
! 𝐁 is transformed into the identity matrix 𝐈 using Guassian elimination
! The row operations that are performed on 𝐁 are simultaneously performed on 𝐃
! with the end result being that 𝐃 is transformed into 𝐁⁻¹𝐃 (i.e. 𝐃 <-- 𝐁⁻¹𝐃)


SUBROUTINE MATINV(N, M, DETERM)
  use variables, only: B, D, ID

  implicit none
  integer, intent(in) :: N, M
  integer :: NN, I, J, JC, IROW, JCOL
  real    :: BMAX, BNEXT, BTRY, F, DETERM

  !     BTRY    - largest absolute value in a row
  !     BNEXT   - second largest absolute value in a row
  !     BMAX    - ratio of the two largest magnitude values in a row (BNEXT / BTRY)
  !     ID      - keeps track of available pivot rows
  !
  ! pivot row    - row with the smallest BMAX value (excluding already used pivot columns)
  ! pivot column - max magnitude value in pivot row

  ! The opposite of >= is <, but strange behavior can happen when encountering NaN's
  ! (.NOT.(BNEXT >= BMAX*BTRY))


  DETERM = 1.0
  ID = 0.0

  do NN = 1,N
      BMAX = 1.1

      do I = 1,N
          if (ID(I) == 0) then
              BNEXT = 0.0
              BTRY = 0.0

              do J = 1,N
                  if (ID(J) == 0) then

                      if (.NOT.(ABS(B(I,J)) <= BNEXT)) then
                          BNEXT = ABS(B(I,J))

                          if (.NOT.(BNEXT <= BTRY)) then
                              BNEXT = BTRY
                              BTRY = ABS(B(I,J))
                              JC = J
                          end if
                      end if
                  end if
              end do                                          ! J=1,N

              if (.NOT.(BNEXT >= BMAX*BTRY)) then             ! BMAX <= BTRY / BNEXT
                  BMAX = BNEXT/BTRY
                  IROW = I
                  JCOL = JC
              end if
          end if                                              ! (ID(I) == 0)
      end do                                                  ! I = 1,N

      if (ID(JC) /= 0) then
        DETERM = 0.0
        return
      end if


      ID(JCOL) = 1
      if (JCOL /= IROW) then
        call swap(B(IROW,:), B(JCOL,:))         ! row swap
        call swap(D(IROW,:), D(JCOL,:))
      end if


      F = 1.0/B(JCOL,JCOL)
      B(JCOL,:)   = B(JCOL,:)   * F
      D(JCOL,1:M) = D(JCOL,1:M) * F

      do I = 1,N
        if (I /= JCOL) then
          F        = B(I,JCOL)
          B(i,:)   = B(i,:)   - F*B(JCOL,:)
          D(i,1:M) = D(i,1:M) - F*D(JCOL,1:M)
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
