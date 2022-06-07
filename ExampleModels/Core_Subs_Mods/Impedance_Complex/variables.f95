! the module variables contains the variables that are subsequently used in
! the modules and subroutines ABDGXY, BAND, MATINV, and GOV_EQNS
! The size of the arrays created depend only on the number of dependent variables (N),
! and the number of mesh points (NJ); both of which are provided in user_input

! edits: June 3, 2022
! real --> complex: cprev, delC
!                   A, B, D, G, X, Y
!                   dW, dE, fW, fE, rj, smG
!                   E, ID

module variables
  use user_input, only: N, NJ

  implicit none

  complex, dimension(N, NJ) :: cprev, delC
  real, dimension(NJ)       :: xx, delX

  ! ABDGXY variables
  complex, dimension(N, N)     :: A, B, X, Y
  complex, dimension(N, 2*N+1) :: D
  complex, dimension(N)        :: G

  ! BAND and MATINV variables
  complex, dimension(N, N+1, NJ) :: E
  complex, dimension(N)          :: ID

  ! auto_fill variables
  real                   :: alphaW, alphaE, betaW, betaE
  complex, dimension(N,N):: dW, dE, fW, fE, rj
  complex, dimension(N)  :: smG
  real, dimension(NJ)    :: Cntrl_Vol       ! Control Volume
  real, dimension(2, NJ) :: Crx_Area        ! dimesion = (2,NJ) because east-west
                                            ! cross-sectional areaCrx_Area

contains

  ! ****************************************************************************
  ! ** functions to define the node control volumes and cross-sectional area ***
  ! ----------------------------------------------------------------------------
  function Control_Volume(Geometry) result(ctrl_vol)
    use user_input, only: PI
    real, dimension(NJ) :: ctrl_vol
    character(len=:), allocatable, intent(in) :: Geometry
    integer :: j

    ctrl_vol(1) = 0.0
    ctrl_vol(NJ) = 0.0

    if (Geometry == 'Rectangular') then
      do j = 2, NJ-1
        ctrl_vol(j) = delX(j)
      end do

    else if (Geometry == 'Cylindrical') then
      do j = 2, NJ-1
        ctrl_vol(j) = PI*( (xx(j) + delX(j)/2.0)**2 - (xx(j) - delX(j)/2.0)**2 )
      end do

    else if (Geometry == 'Spherical') then
      do j = 2, NJ-1
        ctrl_vol(j) = 4.0/3.0*PI*( (xx(j) + delX(j)/2.0)**3 - (xx(j) - delX(j)/2.0)**3 )
      end do

    end if

  end function Control_Volume


  function Cross_Sectional_Area(Geometry) result(Cross_Area)
    use user_input, only: PI
    real, dimension(2, NJ) :: Cross_Area      ! 1 - west side, 2 - east side
    character(len=:), allocatable, intent(in) :: Geometry
    integer :: j

    if (Geometry == 'Rectangular') then
      Cross_Area = 1.0

    else if (Geometry == 'Cylindrical') then
      do j = 1, NJ
        Cross_Area(1,j) = 2.0 * PI * (xx(j) - delX(j)/2.0)
        Cross_Area(2,j) = 2.0 * PI * (xx(j) + delX(j)/2.0)
      end do

    else if (Geometry == 'Spherical') then
      do j = 1, NJ
        Cross_Area(1,j) = 4.0 * PI * (xx(j) - delX(j)/2.0)**2.0
        Cross_Area(2,j) = 4.0 * PI * (xx(j) + delX(j)/2.0)**2.0
      end do

    end if

  end function Cross_Sectional_Area
! ******************************************************************************

end module variables
