! * Written by Nicholas Brady August 10, 2020
! Uses John Newman's Band algorithm to solve coupled non-linear partial
! differential equations
! Incorporates DNAD which anables to use of automatic Differentiation to
! linearize the PDE's
! variables or names that begin and end with '_', i.e. _variable_ are changed by the python program: RunFortran.py

module number_of_variables
  implicit none
  integer, parameter :: N = 3*2+1
end module number_of_variables

module user_input
  use number_of_variables
  use dnadmod
  implicit none

  integer, parameter :: NJ = 1001                          ! Number of mesh points
  integer, parameter :: Numbertimesteps = 30 * 3.6e3     ! Number of time steps
  integer, parameter :: maxIterations = 1e6
  real               :: delT = 1.0e1                       ! size of timestep [s]
  real               :: time                             ! [s]
  real, parameter    :: xmax = 4.0                      ! 100 um

  real, parameter :: Rigc   = 8.314             ! Ideal gas constant [J/(mol*K)]
  real, parameter :: Fconst = 96485             ! Faraday's Constant [C/mol]
  real, parameter :: Temp   = 298
  real, parameter :: PI = 4.0 * ATAN(1.0)       ! pi - Geometric constant

  real, parameter :: c_total = 55.5e-3
  ! real, parameter :: xbulk_O2 = 2e-7 / c_total
  ! real, parameter :: xbulk_OH = 0.0
  real, parameter :: xbulk_Na = 10.0e-3 / c_total
  real, parameter :: xbulk_Cl = 10.0e-3 / c_total
  real, parameter :: xbulk_H2O = 1.0 - xbulk_Na - xbulk_Cl

  ! real, parameter :: D_O2 = 2.00000e-5
  ! real, parameter :: D_OH = 5.26000e-5
  real, parameter :: D_Na = 1.33400e-5
  real, parameter :: D_Cl = 2.03200e-5

  ! real, parameter :: z_O2 =  0.0
  ! real, parameter :: z_OH = 0.0!-1.0
  real, parameter :: z_Na = +1.0
  real, parameter :: z_Cl = -1.0
  real, parameter :: z_H2O = 0.0

  ! real, parameter :: s_O2 = +1.0
  ! real, parameter :: s_OH = -4.0
  real, parameter :: s_Na =  1.0
  real, parameter :: s_Cl =  0.0
  real, parameter :: s_H2O = 0.0!+2.0

  character(len=65) :: geometry = 'Rectangular'   ! system Geometry: Rectangular, Cylindrical, Spherical

end module user_input


module variables
  use user_input
  implicit none

  real, dimension(N, NJ)  :: cprev, delC
  real, dimension(NJ)     :: xx, delX
  real                    :: h

  ! ABDGXY_VARS
  real, dimension(N, N)     :: A, B, X, Y
  real, dimension(N, 2*N+1) :: D
  real, dimension(N)        :: G

  ! BAND and MATINV variables
  real, dimension(N, N+1, NJ) :: E
  real, dimension(N)          :: ID
  integer :: NP1

end module variables

! ******************************************************************************
module write_data_mod
  use user_input
  use variables, only: cprev, delC, xx

  implicit none

  character(len=65) :: header_fmt = '(1(A12,2X),   20(A15,2X)) '
  character(len=65) :: data_fmt   = '(1(F12.5,2X), 20(ES15.8,2X))'

  real :: t_write

  real :: conc_1, conc_NJ, vel_1, vel_NJ
  real :: c_O2, c_OH, c_Na, c_Cl, Phi
  real :: conc, vel, density_

contains

  subroutine t_write__Equiv__Li_Bal_Assignment
    t_write     = time / 3600.0

  end subroutine t_write__Equiv__Li_Bal_Assignment


  subroutine write_condition(it)
    integer :: it
    real :: last_write_time = 0.0
    real :: write_every_x_sec = 1.0e-1           ! 3600 s = 1 hour

    call t_write__Equiv__Li_Bal_Assignment

    call write_to_screen(it)
    call write_positional_information_to_file(it)

  end subroutine write_condition


  subroutine write_to_screen(it)
    integer :: it, j
    real :: xi, zeta

    if (it == 0) then       ! write the headers on the first entrance into write all voltage
      write(*, header_fmt) 'Time',  'Position', 'x_Na', 'x_Cl', 'x_H2O', 'Flux_Na', 'Flux_Cl', 'Flux_H2O', 'Phi'
      write(*, header_fmt) 'hours', 'ξ'       , ''    , ''    , ''     , ''       , ''       , ''        , 'Volts'
    end if

    do j = 1, NJ                                    !
                                                          !
      write(*, data_fmt) t_write, xx(j), cprev(1,j)*55.5, cprev(2,j)*55.5, cprev(3,j)*55.5, cprev(4,j), cprev(5,j), &
                          & cprev(6,j), cprev(7,j)

    end do

    ! write(*, data_fmt) 0.0, xx(1), cprev(6,1), cprev(7,1), cprev(8,1), cprev(9,1), cprev(10,1), cprev(11,1)
    ! write(*, data_fmt) 0.0, xx(NJ), cprev(6,NJ), cprev(7,NJ), cprev(8,NJ), cprev(9,NJ), cprev(10,NJ), cprev(11,NJ)

  end subroutine write_to_screen


  subroutine write_positional_information_to_file(it)
    use variables, only: xx
    integer :: it, j
    real :: xi, zeta

    open(56, file = 'Time_Conc_Position.txt', status = 'unknown')

    if (it == 0) then
!           write the headers on the first entrance into write all voltage
      write(56, header_fmt) 'Time',  'Position', 'x_Na', 'x_Cl', 'x_H2O', 'Flux_Na', 'Flux_Cl', 'Flux_H2O', 'Phi'
      write(56, header_fmt) 'hours', 'ξ'       , ''    , ''    , ''     , ''       , ''       , ''        , 'Volts'
                                                      !
    end if                                              !


    do j = 1, NJ                                    !
                                                          !
      write(56, data_fmt) t_write, xx(j), cprev(1,j)*55.5, cprev(2,j)*55.5, cprev(3,j)*55.5, cprev(4,j), cprev(5,j), &
                        & cprev(6,j), cprev(7,j)

    end do

  end subroutine write_positional_information_to_file

end module write_data_mod


module Solution_Properties_Mod
  use user_input
  use dnadmod
  implicit none

contains
  function d_chemical_potential_dx(c_vars, dcdx_vars)   result(c_dmudx)
    type(dual), dimension(N)              :: c_dmudx
    type(dual), dimension(N), intent(in)  :: c_vars, dcdx_vars
    type(dual)  :: c_O2_, dc_O2_dx_

    c_O2_ = c_vars(1)
    dc_O2_dx_ = dcdx_vars(1)

    c_dmudx(1)    = Rigc*Temp * dc_O2_dx_
    c_dmudx(2:N)  = 0.0

  end function d_chemical_potential_dx


end module Solution_Properties_Mod



module GOV_EQNS
  use user_input
  use variables
  use dnadmod
  use Solution_Properties_Mod

  implicit none

  ! initialize the fill_mat variables - shared with auto_fill and ABDGXY
  real                  :: alphaW, alphaE, betaW, betaE
  real, dimension(N,N)  :: dW, dE, fW, fE, rj
  real, dimension(N,N)  :: smP, smE, smA, smB, smD
  real, dimension(N)    :: smG, smF
  real, dimension(NJ)   :: Cntrl_Vol
  real, dimension(2, NJ):: Crx_Area        ! dimesion = (2,NJ) because east-west
                                           ! cross-sectional area
  ! type(dual), dimension(N) :: c_vars_dual, dcdx_vars_dual

  type(dual) :: x_O2, dcdx_O2, d2cdx2_O2
  type(dual) :: x_OH, dcdx_OH, d2cdx2_OH
  type(dual) :: x_Na, dcdx_Na, d2cdx2_Na
  type(dual) :: x_Cl, dcdx_Cl, d2cdx2_Cl
  type(dual) :: x_H2O, dcdx_H2O, d2cdx2_H2O
  type(dual) :: Phi, dPhidx, d2Phidx2
  type(dual) :: N_O2, dNdx_O2     ! Fluxes
  type(dual) :: N_OH, dNdx_OH     ! Fluxes
  type(dual) :: N_Na, dNdx_Na     ! Fluxes
  type(dual) :: N_Cl, dNdx_Cl     ! Fluxes
  type(dual) :: N_H2O, dNdx_H2O     ! Fluxes
  real       :: position


contains

! ******************************************************************************
! **************************** GOVERNING EQUATIONS *****************************
! **************** Accumulation = FluxIn - FluxOut + Generation ****************
! ******************************************************************************
! can define intermediate variables with the functions to improve readability
! i.e c0 = c_vars_dual(1), dPhi2_dx = dcdx_vars_dual(2)

  function GOV_RHS(c_vars_dual, dcdx_vars_dual, d2cdx2_vars_dual, position)  result(GOV_)
    type(dual), dimension(N)              :: GOV_
    type(dual), dimension(N), intent(in)  :: c_vars_dual, dcdx_vars_dual, d2cdx2_vars_dual
    real, intent(in)                      :: position

    x_Na       = c_vars_dual(1)
    dcdx_Na    = dcdx_vars_dual(1)
    d2cdx2_Na  = d2cdx2_vars_dual(1)

    x_Cl       = c_vars_dual(2)
    dcdx_Cl    = dcdx_vars_dual(2)
    d2cdx2_Cl  = d2cdx2_vars_dual(2)

    x_H2O      = c_vars_dual(3)
    dcdx_H2O   = dcdx_vars_dual(3)
    d2cdx2_H2O = d2cdx2_vars_dual(3)

    N_Na       = c_vars_dual(4)
    dNdx_Na    = dcdx_vars_dual(4)

    N_Cl       = c_vars_dual(5)
    dNdx_Cl    = dcdx_vars_dual(5)

    N_H2O      = c_vars_dual(6)
    dNdx_H2O   = dcdx_vars_dual(6)

    Phi        = c_vars_dual(7)
    dPhidx     = dcdx_vars_dual(7)
    d2Phidx2   = d2cdx2_vars_dual(7)

    !             dN/dη = 2η dc/dη
    GOV_(1)    = 2*position*dcdx_Na - dNdx_Na
    GOV_(2)    = 2*position*dcdx_Cl - dNdx_Cl
    GOV_(3)    = 2*position*dcdx_H2O - dNdx_H2O

    !             Σ x_i = 1
    GOV_(4)    = x_Na + x_Cl + x_H2O - 1.0

    !           dx_i/dη + z_i F x_i / RT * dΦ/dη = Σ D_OR / D_ij * (x_i N_j - x_j N_i)
    !           only the diffusion coefficients with H2O are significant
    GOV_(5)    = dcdx_Na + z_Na*Fconst*x_Na/(Rigc*Temp)*dPhidx - D_Na/D_Na*(x_Na * N_H2O - x_H2O * N_Na)
    GOV_(6)   = dcdx_Cl + z_Cl*Fconst*x_Cl/(Rigc*Temp)*dPhidx - D_Na/D_Cl*(x_Cl * N_H2O - x_H2O * N_Cl)

    !             Σ c_i z_i = 0 = Σ x_i z_i
    GOV_(7)   = z_Na*x_Na + z_Cl*x_Cl + z_H2O*x_H2O - 0.0
    ! GOV_(11) = Phi - 0.0

  end function GOV_RHS

  ! time-varying portion of governing equation
  function ACCUM(c_vars_dual, position)                           result(Accum_)
    type(dual), dimension(N)              :: Accum_
    type(dual), dimension(N), intent(in)  :: c_vars_dual
    real, intent(in)                      :: position

    Accum_    = 0.0
    Accum_%x  = 0.0

  end function ACCUM
! ******************************************************************************


! ******************************************************************************
! **************************** BOUNDARY CONDITIONS *****************************
! ******************************************************************************
! boundary conditions need to be written as homogeneous equations
! i.e.    c = 1.0   -->      c - 1.0 = 0.0    -->   BC_WEST_(N) = c - 1.0
! i.e. flux = 1.0   -->   flux - 1.0 = 0.0    -->   BC_WEST_(N) = flux - 1.0
! ********** special case **********
! i.e. dc/dt = rxn  -->   dc/dt - rxn = 0.0
! dc_n   = c_vars_dual(n)
! dc_n%x = 0.0      -->   dc        -->           BC_WEST_(2) = dc_n/delT - rxn
  function Boundary_WEST (c_vars_dual, dcdx_vars_dual, position) result(BC_WEST_)
    type(dual), dimension(N)               :: BC_WEST_
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual
    real, intent(in)                       :: position

    x_Na       = c_vars_dual(1)
    dcdx_Na    = dcdx_vars_dual(1)

    x_Cl       = c_vars_dual(2)
    dcdx_Cl    = dcdx_vars_dual(2)

    x_H2O      = c_vars_dual(3)
    dcdx_H2O   = dcdx_vars_dual(3)

    N_Na       = c_vars_dual(4)
    dNdx_Na    = dcdx_vars_dual(4)

    N_Cl       = c_vars_dual(5)
    dNdx_Cl    = dcdx_vars_dual(5)

    N_H2O      = c_vars_dual(6)
    dNdx_H2O   = dcdx_vars_dual(6)

    Phi        = c_vars_dual(7)
    dPhidx     = dcdx_vars_dual(7)

    ! O2 + 4e- + 2*H2O --> 4 OH-
    ! x_O2 = 0.0
    BC_WEST_(1) = x_Na - xbulk_Na*0.5
    BC_WEST_(2) = N_Cl - 0.0
    BC_WEST_(3) = N_H2O - 0.0

    !             Σ x_i = 1
    BC_WEST_(4)    = x_Na + x_Cl + x_H2O - 1.0

    !           dx_i/dη + z_i F x_i / RT * dΦ/dη = Σ D_OR / D_ij * (x_i N_j - x_j N_i)
    BC_WEST_(5)   = dcdx_Na + z_Na*Fconst*x_Na/(Rigc*Temp)*dPhidx - D_Na/D_Na*(x_Na * N_H2O - x_H2O * N_Na)
    BC_WEST_(6)   = dcdx_Cl + z_Cl*Fconst*x_Cl/(Rigc*Temp)*dPhidx - D_Na/D_Cl*(x_Cl * N_H2O - x_H2O * N_Cl)

    !             Σ c_i z_i = 0 = Σ x_i z_i
    BC_WEST_(7)   = z_Na*x_Na + z_Cl*x_Cl + z_H2O*x_H2O - 0.0
    ! BC_WEST_(11) = Phi - 0.0

  end function Boundary_WEST

  function Boundary_EAST (c_vars_dual, dcdx_vars_dual, position) result(BC_EAST_)
    type(dual), dimension(N)               :: BC_EAST_
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual
    real, intent(in)                       :: position
    type(dual) :: c_O2_accum, c_H2O_accum

    x_Na       = c_vars_dual(1)
    dcdx_Na    = dcdx_vars_dual(1)

    x_Cl       = c_vars_dual(2)
    dcdx_Cl    = dcdx_vars_dual(2)

    x_H2O      = c_vars_dual(3)
    dcdx_H2O   = dcdx_vars_dual(3)

    N_Na       = c_vars_dual(4)
    dNdx_Na    = dcdx_vars_dual(4)

    N_Cl       = c_vars_dual(5)
    dNdx_Cl    = dcdx_vars_dual(5)

    N_H2O      = c_vars_dual(6)
    dNdx_H2O   = dcdx_vars_dual(6)

    Phi        = c_vars_dual(7)
    dPhidx     = dcdx_vars_dual(7)

    !             dN/dη = 2η dc/dη
    BC_EAST_(1) = 2*position*dcdx_Na - dNdx_Na
    BC_EAST_(2) = 2*position*dcdx_Cl - dNdx_Cl
    BC_EAST_(3) = 2*position*dcdx_H2O - dNdx_H2O

    !             Σ x_i = 1
    BC_EAST_(4)    = x_Na + x_Cl + x_H2O - 1.0

    ! x_i = x_i_bulk
    ! Φ = 0.0
    BC_EAST_(5)   = x_Na - xbulk_Na
    BC_EAST_(6)   = x_Cl - xbulk_Cl
    BC_EAST_(7)   = Phi - 0.0

  end function Boundary_EAST
! ******************************************************************************


  ! ****************************************************************************
  ! ******* functions to convert concentration and dcdx variables to dual ******
  ! ----------------------------------------------------------------------------
  function c_to_dual(c_vars)          result(c_dual)
    type(dual), dimension(N)          :: c_dual
    real, dimension(N), intent(in)    :: c_vars
    real, dimension(ndv)              :: dx_array
    INTEGER :: ic

    do ic = 1, N
      dx_array = 0.0              ! set the dx_array to zero (all elements)
      dx_array(ic) = 1.0

      c_dual(ic) = dual(c_vars(ic), dx_array)
    end do
  end function c_to_dual


  function dcdx_to_dual(dcdx)         result(dcdx_dual)
    type(dual), dimension(N)          :: dcdx_dual
    real, dimension(N), intent(in)    :: dcdx
    real, dimension(ndv)              :: dx_array
    INTEGER :: ic

    do ic = 1, N
      dx_array = 0.0              ! set the dx_array to zero (all elements)
      dx_array(N+ic) = 1.0

      dcdx_dual(ic) = dual(dcdx(ic), dx_array)
    end do

  end function dcdx_to_dual

  function d2cdx2_to_dual(d2cdx2)         result(d2cdx2_dual)
    type(dual), dimension(N)          :: d2cdx2_dual
    real, dimension(N), intent(in)    :: d2cdx2
    real, dimension(ndv)              :: dx_array
    INTEGER :: ic

    do ic = 1, N
      dx_array = 0.0              ! set the dx_array to zero (all elements)
      dx_array(2*N+ic) = 1.0

      d2cdx2_dual(ic) = dual(d2cdx2(ic), dx_array)
    end do

  end function d2cdx2_to_dual
! ******************************************************************************


  ! ****************************************************************************
  ! ** functions to define the node control volumes and cross-sectional area ***
  ! ----------------------------------------------------------------------------
  function Control_Volume(Geometry) result(ctrl_vol)
    use variables, only: xx, delX

    real, dimension(NJ) :: ctrl_vol
    character(len=:), allocatable, intent(in) :: Geometry
    integer :: j

    print*, Geometry
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
    use variables, only: xx, delX
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

end module GOV_EQNS


! ******************************************************************************
! ******************************** MAIN PROGRAM ********************************
! ******************************************************************************

program unsteady
  use user_input
  use variables, only : cprev, delC
  use write_data_mod

  implicit none
  integer :: t1, t2, clock_rate, clock_max
  integer :: it, j
  real :: tol = 1e-4, convergence_error = 1e3

  it = 0

  call system_clock(t1,clock_rate,clock_max)
  call initial_condition()
  call write_condition(it) ! writes the headers
  it = 1

  do while ((convergence_error > tol).AND.(it <= maxIterations)) !it = 1, 360*2 !Numbertimesteps

      ! **********************************************************************
      ! ******************************** BOUND VAL ***************************
      ! **********************************************************************
        do j = 1, NJ
          call auto_fill(j)
          call ABDGXY(j)
          call BAND(j)
        end do

        cprev = cprev + delC      ! Update the dependent variables
        time  = time + delT
        ! ********************************************************************
        ! call write_condition(it)

        convergence_error = MAXVAL(ABS(delC))
        it = it + 1

  end do

  call write_condition(it)



  call system_clock(t2,clock_rate,clock_max)
  write ( *, * ) 'Elapsed real time =', real(t2-t1)/real(clock_rate )

end program unsteady

! ------------------------------------------------------------------------------
! ***************************** end MAIN PROGRAM *******************************
! ------------------------------------------------------------------------------


subroutine initial_condition()
  use user_input
  use variables
  use GOV_EQNS, only: Control_Volume, Cross_Sectional_Area,  &    ! functions
                      Cntrl_Vol, Crx_Area                         ! variables

  implicit none
  integer :: j

  character(len=:), allocatable :: control_volume_input

  h = xmax/float(nj-1)
  do j=1,NJ

    if (j.EQ.1) then
      xx(j) = 0.0
    else if (j.EQ.NJ) then
      xx(NJ) = 0.0 + xmax
    else
      xx(j) = xx(1) + h*float(j-1)
    end if

  end do

  do j=2, NJ-1
     delx(j) = h
  end do
     delx(1) = h/2                      ! it is common in the finite volume
     delx(NJ) = h/2                     ! algorithm to set the control
                                        ! volumes at the boundaries to zero
  cprev(1,:) = xbulk_Na
  cprev(2,:) = xbulk_Cl
  cprev(3,:) = xbulk_H2O

  cprev(4,:) = 0.0
  cprev(5,:) = 0.0
  cprev(6,:) = 0.0

  cprev(7,:) = 0.0

  control_volume_input = trim(geometry)                   ! Define the control volume
  Cntrl_Vol = Control_Volume(control_volume_input)        ! size based on the
  Crx_Area = Cross_Sectional_Area(control_volume_input)   ! specified system geometry

  return                                                  ! is this necessary?
end subroutine initial_condition



!------------------------------------------------------------------------------!
!********************************* AUTO_FILL **********************************!
!------------------------------------------------------------------------------!
subroutine auto_fill(j)
  use user_input, only : N, NJ
  use GOV_EQNS
  use dnadmod
  use variables

  implicit none

  integer :: j
  integer :: ic

  ! variables and their derivatives at the control volume interfaces
  real, dimension(N) :: dcdx, d2cdx2

  ! dual variables
  type(dual), dimension(N) :: cj_dual, dcdx_dual, d2cdx2_dual
  type(dual), dimension(N) :: GOV_dual, accumulation_dual
  type(dual), dimension(N) :: boundary_conditionW, boundary_conditionE

  ! set all matrix variables (dW, dE, fW, fE, rj, smG) to 0.0
  smA = 0.0
  smB = 0.0
  smD = 0.0
  smG = 0.0

  smP = 0.0
  smE = 0.0
  smF = 0.0

  !-----------------------------------------------------------------------------
  ! Boundary Conditions
  !-----------------------------------------------------------------------------
  if (j == 1) then                              ! West Side Boundary Conditions
    dcdx = (-3.0*cprev(:,j) + 4.0*cprev(:,j+1) - cprev(:,j+2)) /(2.0 * h)

    cj_dual     = c_to_dual(cprev(:,j))
    dcdx_dual   = dcdx_to_dual(dcdx)

    boundary_conditionW = Boundary_WEST(cj_dual, dcdx_dual, xx(j))

    do ic = 1,N
        smP(ic, :) = boundary_conditionW(ic)%dx(N+1:2*N)
        smE(ic, :) = boundary_conditionW(ic)%dx(1:N)

        smF(ic)    = -boundary_conditionW(ic)%x
    end do

                                                !-------------------------------
  else if (j == NJ) then                        ! East Side Boundary Conditions
    dcdx = (3.0*cprev(:,j) - 4.0*cprev(:,j-1) + cprev(:,j-2)) /(2.0 * h)

    cj_dual     = c_to_dual(cprev(:,j))
    dcdx_dual   = dcdx_to_dual(dcdx)

    boundary_conditionE = Boundary_EAST(cj_dual, dcdx_dual, xx(j))

    do ic = 1,N
        smP(ic, :) = boundary_conditionE(ic)%dx(N+1:2*N)
        smE(ic, :) = boundary_conditionE(ic)%dx(1:N)

        smF(ic)    = -boundary_conditionE(ic)%x
    end do

  !-----------------------------------------------------------------------------
  ! Governing Equations
  !-----------------------------------------------------------------------------
  else
    dcdx   = (cprev(:,j+1) - cprev(:,j-1)) /(2.0 * h)
    d2cdx2 = (cprev(:,j+1) - 2.0*cprev(:,j) + cprev(:,j-1))/h**2

    cj_dual     = c_to_dual(cprev(:,j))
    dcdx_dual   = dcdx_to_dual(dcdx)
    d2cdx2_dual = d2cdx2_to_dual(d2cdx2)

    GOV_dual = GOV_RHS(cj_dual, dcdx_dual, d2cdx2_dual, xx(j))
    accumulation_dual = ACCUM(cj_dual, xx(j))

    do ic = 1,N                                   ! ic - equation number
        smA(ic,:) = GOV_dual(ic)%dx(2*N+1:3*N)    ! smA * d2cdx2 +
        smB(ic,:) = GOV_dual(ic)%dx(N+1:2*N)      ! smB * dcdx   +
        smD(ic,:) = GOV_dual(ic)%dx(1:N)          ! smD * c      +

        smG(ic)   = -GOV_dual(ic)%x                ! smG = 0

    end do

  end if

end subroutine auto_fill

!------------------------------------------------------------------------------!
!*********************************** ABDGXY ***********************************!
!------------------------------------------------------------------------------!
subroutine ABDGXY(j)
  ! ABDGXY equates the large coefficents based on the small coefficents.
  ! The coefficents A, B, D, G, X, Y can be found in Newman appendix C.
      use user_input, only: N, NJ
      use GOV_EQNS
      use variables, only: A, B, D, G, X, Y
      implicit none
      integer :: j, i


      ! if (j.eq.1) then
      !   X = 0.0
      !   B = rj - (1.0 - alphaE)*fE + betaE*dE
      !   D(1:N,1:N) = -alphaE*fE - betaE*dE
      !   G = smG
      !
      !   return
      ! end if
      !
      ! if (j.eq.NJ) then
      !   Y = 0.0
      !   A = (1.d0 - alphaW)*fW - betaW*dW
      !   B = rj + betaW*dW + alphaW*fW
      !   G = smG
      !
      !   return
      ! end if
      !
      ! A = (1.d0 - alphaW)*fW - betaW*dW
      ! B = rj + betaW*dW + alphaW*fW - (1.0 - alphaE)*fE + betaE*dE
      ! D(1:N,1:N) = -alphaE*fE - betaE*dE
      ! G = smG
      !
      ! return

      if(j.eq.1) then
        B = h*smE - 1.5*smP
        D(1:N,1:N) = 2.0*smP
        X = -0.5*smP
        G = h*smF
        return
      end if

      if (j.eq.NJ) then
        A = -2.0*smP
        B = h*smE + 1.5*smP
        Y = 0.5*smP
        G = h*smF
        return
      end if

      A = smA - h/2.0*smB
      B = -2.0*smA + h**2*smD
      D(1:N,1:N) = smA + h/2.0*smB
      G = h**2*smG

      return

end subroutine ABDGXY

!------------------------------------------------------------------------------!
!*********************************** MATINV ***********************************!
!------------------------------------------------------------------------------!
SUBROUTINE MATINV(N, M, DETERM)
  use variables, only: A, B, delC, D, ID ! A imported but not used
  implicit double precision (A-H,O-Z) ! implicits are not good coding practice
 ! use variables, only: delC ! A imported but not used
 ! use ABDGXY_VARS, only: A, B, D
 ! implicit double precision (A-H,O-Z)
 ! real, dimension(N) :: ID

      DETERM=1.0
      ! ID = 0.0
      DO 1 I=1,N
1       ID(I)=0
      DO 18 NN=1,N
        BMAX=1.1
        DO 6 I=1,N
          IF (ID(I).NE.0) GOTO 6
          BNEXT=0.0
          BTRY=0.0
          DO 5 J=1,N
            IF (ID(J).NE.0) GOTO 5
            IF (ABS(B(I,J)).LE.BNEXT) GOTO 5
            BNEXT=ABS(B(I,J))
            IF (BNEXT.LE.BTRY) GOTO 5
            BNEXT=BTRY
            BTRY=ABS(B(I,J))
            JC=J
5           CONTINUE
          IF (BNEXT.GE.BMAX*BTRY) GOTO 6
          BMAX=BNEXT/BTRY
          IROW=I
          JCOL=JC
6       CONTINUE
        IF (ID(JC).EQ.0) GOTO 8
        DETERM=0.0
        RETURN
8       ID(JCOL)=1
        IF (JCOL.EQ.IROW) GOTO 12
9       DO 10 J=1,N
          SAVE=B(IROW,J)
          B(IROW,J)=B(JCOL,J)
10        B(JCOL,J)=SAVE
        DO 11 K=1,M
          SAVE=D(IROW,K)
          D(IROW,K)=D(JCOL,K)
11        D(JCOL,K)=SAVE
12      F=1.0/B(JCOL,JCOL)
        DO 13 J=1,N
13        B(JCOL,J)=B(JCOL,J)*F
        DO 14 K=1,M
14        D(JCOL,K)=D(JCOL,K)*F
        DO 18 I=1,N
          IF (I.EQ.JCOL) GOTO 18
          F=B(I,JCOL)
          DO 16 J=1,N
16          B(I,J)=B(I,J)-F*B(JCOL,J)
          DO 17 K=1,M
17          D(I,K)=D(I,K)-F*D(JCOL,K)
18    CONTINUE
      RETURN
      end

!------------------------------------------------------------------------------!
!************************************ BAND ************************************!
!------------------------------------------------------------------------------!
SUBROUTINE BAND(J)
! BAND(J) computes delC and calls MATINV to solve the problem using gaussian elimination.
use variables, only: A, B, delC, D, G, X, Y, NP1, E
! use variables, only: delC
use user_input, only: N, NJ
! use ABDGXY_VARS
! use BAND_J_VARS
implicit double precision (A-H,O-Z)
! integer :: NP1

101   FORMAT(15H DETERM=0 AT J=,I4)
      IF (J-2) 1,6,8
1     NP1=N+1
      DO 2 I=1,N
        D(I,2*N+1)=G(I)
        DO 2 L=1,N
          LPN=L+N
2         D(I,LPN)=X(I,L)
      CALL MATINV(N,2*N+1,DETERM)
      IF (DETERM) 4,3,4
3     PRINT 101,J
4     DO 5 K=1,N
        E(K,NP1,1)=D(K,2*N+1)
        DO 5 L=1,N
          E(K,L,1)=-D(K,L)
          LPN=L+N
5         X(K,L)=-D(K,LPN)
      RETURN
6     DO 7 I=1,N
        DO 7 K=1,N
          DO 7 L=1,N
7           D(I,K)=D(I,K)+A(I,L)*X(L,K)
8     IF (J-NJ) 11,9,9
9     DO 10 I=1,N
        DO 10 L=1,N
          G(I)=G(I)-Y(I,L)*E(L,NP1,J-2)
          DO 10 M=1,N
10          A(I,L)=A(I,L) + Y(I,M)*E(M,L,J-2)
11    DO 12 I=1,N
        D(I,NP1)=-G(I)
        DO 12 L=1,N
          D(I,NP1)=D(I,NP1)+A(I,L)*E(L,NP1,J-1)
          DO 12 K=1,N
12          B(I,K)=B(I,K) + A(I,L)*E(L,K,J-1)
      CALL MATINV(N,NP1,DETERM)
      IF (DETERM) 14,13,14
13    PRINT 101,J
14    DO 15 K=1,N
        DO 15 M=1,NP1
15        E(K,M,J)=-D(K,M)
      IF (J-NJ) 20,16,16
16    DO 17 K=1,N
17      delC(K,J)=E(K,NP1,J)
      DO 18 JJ=2,NJ
        M=NJ-JJ+1
        DO 18 K=1,N
          delC(K,M)=E(K,NP1,M)
          DO 18 L=1,N
18          delC(K,M)=delC(K,M) +E(K,L,M)*delC(L,M+1)
      DO 19 L=1,N
        DO 19 K=1,N
19        delC(K,1)=delC(K,1)+X(K,L)*delC(L,3)
20    RETURN
      end



!******************************************************************************
!* dual Number Automatic Differentiation (DNAD) of Fortran Codes
!*-----------------------------------------------------------------------------
!* COPYRIGHT (c) Joshua Hodson, All rights reserved, you are free to copy,
!* modify, or translate this code to other languages such as c/c++. This is a
!* fork of the original Fortran DNAD module developed by Dr. Wenbin Yu. See
!* original copyright information below. You can download the original version
!* at https://cdmhub.org/resources/374
!*
!* COPYRIGHT (c) Wenbin Yu, All rights reserved, you are free to copy,
!* modify or translate this code to other languages such as c/c++. If
!* you find a bug please let me know through wenbinyu.heaven@gmail.com. If
!* you added new functions and want to share with others, please let me know
!* too. You are welcome to share your successful stories with us through
!* http://groups.google.com/group/hifi-comp.
!******************************************************************************
!* Acknowledgements
!*-----------------------------------------------------------------------------
!* The development of DNAD is supported, in part, by the Chief Scientist
!* Innovative Research Fund at AFRL/RB WPAFB, and by Department of Army
!* SBIR (Topic A08-022) through Advanced Dynamics Inc. The views and
!* conclusions contained herein are those of the authors and should not be
!* interpreted as necessarily representing the official policies or
!* endorsement, either expressed or implied, of the funding agency.
!*
!* Additional development of DNAD has been supported under a Department of
!* Energy (DOE) Nuclear Energy University Program (NEUP) Graduate Fellowship.
!* Any opinions, findings, conclusions or recommendations expressed in this
!* publication are those of the authors and do not necessarily reflect the
!* views of the Department of Energy Office of Nuclear Energy.
!******************************************************************************
!* Citation
!*-----------------------------------------------------------------------------
!* Your citation of the following two papers is appreciated:
!* Yu, W. and Blair, M.: "DNAD, a Simple Tool for Automatic Differentiation of
!* Fortran Codes Using dual Numbers," Computer Physics Communications, vol.
!* 184, 2013, pp. 1446-1452.
!*
!* Spall, R. and Yu, W.: "Imbedded dual-Number Automatic Differentiation for
!* CFD Sensitivity Analysis," Journal of Fluids Engineering, vol. 135, 2013,
!* 014501.
!******************************************************************************
!* Quick Start Guide
!*-----------------------------------------------------------------------------
!* To integrate DNAD into an existing Fortran program, do the following:
!*
!*   1. Include the DNAD module in the source files by adding "use dnadmod" to
!*      the beginning of all modules, global functions, and global subroutines
!*      that include definitions of floating-point variables.
!*   2. Redefine all floating-point variables as type(dual). This can be done
!*      using precompiler directives so that the integration can be turned on
!*      or off at compile-time, eliminating the need for maintaining two
!*      separate code bases for the same project.
!*   3. All I/O involving floating-point variables will need to be examined.
!*      A method will need to be determined for inputting and outputting
!*      derivative values. This customization is typically unique for each
!*      piece of software and needs to be determined on a case-by-case basis.
!*   4. When compiling DNAD, use the compiler option "-Dndv=#", where # is the
!*      number of design variables desired. This sizes the derivative array
!*      that is stored with each floating point number.
!*   5. When compiling DNAD, use compiler options to specify precision. If no
!*      compiler options are specified, DNAD will default to single-precision
!*      floating-point arithmetic. Most popular Fortran compilers provide
!*      options for specifying precision at compile-time so that it does not
!*      have to be hard-coded into the source code. For example, use the
!*      "-fdefault-real-8" compiler in gfortran or the "-r8" compiler option
!*      with Intel Fortran to compile DNAD as double-precision.
!*   6. Modify the compilation process for the target software to include the
!*      DNAD module in the resulting executable or library.
!******************************************************************************
!* Change Log
!*-----------------------------------------------------------------------------
!*
!*  2016-04-29  Joshua Hodson
!*  - Updated copyright, acknowledgments, and quick start guide.
!*  - Removed overloads for single-precision reals.
!*  - Added tan, dtan, atan, and atan2 intrinsic function overloads.
!*  - Removed macro for precision and defined all floating-point variables as
!*    default real. Compiler options can now be used to set precision.
!*  - Added checks for undefined derivatives when only constants are used in
!*    the calculation (i.e. all partial derivatives are zero). This limits the
!*    perpetuation of NaN values in the code.
!*  - Combined the header and source files into a single file.
!*
!*  2015-07-29  Joshua Hodson
!*  - Added maxloc intrinsic function overload.
!*  - Converted UPPERCASE to lowercase for readability.
!*  - Added macros for defining precision and number of design variables.
!*  - Renamed module from dual_Num_Auto_Diff to dnadmod
!*  - Renamed dual number type from dual_NUM to dual
!*  - Renamed components of dual number type from (xp_ad_, xp_ad_) to (x, dx)
!*
!*  2014-06-05  Wenbin Yu
!*  - Forked from original DNAD repository, see https://cdmhub.org/resources/374
!*
!******************************************************************************

! Number of design variables (default = 1)
! #ifndef ndv
! #define ndv 1
! #endif

module dnadmod

    ! use user_input, only: N
  use number_of_variables
    implicit none

    integer, PARAMETER, public :: ndv = N*3   ! cprev_vars, dcdx_vars, d2cdx2_vars

    private

    real :: negative_one = -1.0
    type,public :: dual  ! make this private will create difficulty to use the
                        ! original write/read commands, hence x and dx are
                        ! variables which can be accessed using D%x and D%dx in
                        ! other units using this module in which D is defined
                        ! as type(dual).
        sequence
        real :: x  ! functional value
        real :: dx(ndv)  ! derivative
    end type dual


!******** Interfaces for operator overloading
    public assignment (=)
    interface assignment (=)
        module procedure assign_di  ! dual=integer, elemental
        module procedure assign_dr  ! dual=real, elemental
        module procedure assign_id  ! integer=dual, elemental
    end interface


    public operator (+)
    interface operator (+)
        module procedure add_d   ! +dual number, elemental
        module procedure add_dd  ! dual + dual, elemental
        module procedure add_di  ! dual + integer, elemental
        module procedure add_dr  ! dual + real, elemental
        module procedure add_id  ! integer + dual, elemental
        module procedure add_rd  ! real + dual, elemental
    end interface

    public operator (-)
    interface operator (-)
        module procedure minus_d   ! negate a dual number,elemental
        module procedure minus_dd  ! dual -dual,elemental
        module procedure minus_di  ! dual-integer,elemental
        module procedure minus_dr  ! dual-real,elemental
        module procedure minus_id  ! integer-dual,elemental
        module procedure minus_rd  ! real-dual,elemental
    end interface

    public operator (*)
    interface operator (*)
        module procedure mult_dd    ! dual*dual, elemental
        module procedure mult_di    ! dual*integer,elemental
        module procedure mult_dr    ! dual*real,elemental
        module procedure mult_id    ! integer*dual,elemental
        module procedure mult_rd    ! real*dual,elemental
    end interface

    public operator (/)
    interface operator (/)
        module procedure div_dd ! dual/dual,elemental
        module procedure div_di ! dual/integer, elemental
        module procedure div_dr ! dual/real,emental
        module procedure div_id ! integer/dual, elemental
        module procedure div_rd ! real/dual, elemental
    end interface

    public operator (**)
    interface operator (**)
        module procedure pow_i ! dual number to an integer power,elemental
        module procedure pow_r ! dual number to a real power, elemental
        module procedure pow_d ! dual number to a dual power, elemental
    end interface

    public operator (==)
    interface operator (==)
        module procedure eq_dd ! compare two dual numbers, elemental
        module procedure eq_di ! compare a dual and an integer, elemental
        module procedure eq_dr ! compare a dual and a real, elemental
        module procedure eq_id ! compare integer with a dual number, elemental
        module procedure eq_rd ! compare a real with a dual number, elemental
    end interface

    public operator (<=)
    interface operator (<=)
        module procedure le_dd  ! compare two dual numbers, elemental
        module procedure le_di  ! compare a dual and an integer, elemental
        module procedure le_dr  ! compare a dual and a real,elemental
        module procedure le_id ! compare integer with a dual number, elemental
        module procedure le_rd ! compare a real with a dual number, elemental
    end interface

    public operator (<)
    interface operator (<)
        module procedure lt_dd  !compare two dual numbers, elemental
        module procedure lt_di  !compare a dual and an integer, elemental
        module procedure lt_dr  !compare dual with a real, elemental
        module procedure lt_id ! compare integer with a dual number, elemental
        module procedure lt_rd ! compare a real with a dual number, elemental
    end interface

    public operator (>=)
    interface operator (>=)
        module procedure ge_dd ! compare two dual numbers, elemental
        module procedure ge_di ! compare dual with integer, elemental
        module procedure ge_dr ! compare dual with a real number, elemental
        module procedure ge_id ! compare integer with a dual number, elemental
        module procedure ge_rd ! compare a real with a dual number, elemental
    end interface

    public operator (>)
    interface operator (>)
        module procedure gt_dd  !compare two dual numbers, elemental
        module procedure gt_di  !compare a dual and an integer, elemental
        module procedure gt_dr  !compare dual with a real, elemental
        module procedure gt_id ! compare integer with a dual number, elemental
        module procedure gt_rd ! compare a real with a dual number, elemental
    end interface

    public operator (/=)
    interface operator (/=)
        module procedure ne_dd  !compare two dual numbers, elemental
        module procedure ne_di  !compare a dual and an integer, elemental
        module procedure ne_dr  !compare dual with a real, elemental
        module procedure ne_id ! compare integer with a dual number, elemental
        module procedure ne_rd ! compare a real with a dual number, elemental
    end interface


!------------------------------------------------
! Interfaces for intrinsic functions overloading
!------------------------------------------------
    public abs
    interface abs
        module procedure abs_d  ! absolute value of a dual number, elemental
    end interface

    public dabs
    interface dabs
        module procedure abs_d ! same as abs, used for some old fortran commands
    end interface

    public acos
    interface acos
        module procedure acos_d ! arccosine of a dual number, elemental
    end interface

    public asin
    interface asin
        module procedure asin_d ! arcsine of a dual number, elemental
    end interface

    public atan
    interface atan
        module procedure atan_d ! arctan of a dual number, elemental
    end interface

    public atan2
    interface atan2
        module procedure atan2_d ! arctan of a dual number, elemental
    end interface

    public cos
    interface cos
        module procedure cos_d ! cosine of a dual number, elemental
    end interface

    public dcos
    interface dcos
        module procedure cos_d ! cosine of a dual number, elemental
    end interface

    public dot_product
    interface dot_product
        module procedure dot_product_dd ! dot product two dual number vectors
    end interface

    public exp
    interface exp
        module procedure exp_d ! exponential of a dual number, elemental
    end interface

    public int
    interface int
        module procedure int_d ! integer part of a dual number, elemental
    end interface

    public log
    interface log
        module procedure log_d ! log of a dual number, elemental
    end interface

    public log10
    interface log10
        module procedure log10_d ! log of a dual number, elemental
    end interface

    public matmul
    interface matmul
        module procedure matmul_dd ! multiply two dual matrices
        module procedure matmul_dv ! multiply a dual matrix with a dual vector
        module procedure matmul_vd ! multiply a dual vector with a dual matrix
    end interface


    public max
    interface max
        module procedure max_dd ! max of from two to four dual numbers, elemental
        module procedure max_di ! max of a dual number and an integer, elemental
        module procedure max_dr ! max of a dual number and a real, elemental
        module procedure max_rd ! max of a real,and a dual number,  elemental
    end interface

    public dmax1
    interface dmax1
        module procedure max_dd ! max of from two to four dual numbers, elemental
    end interface

    public maxval
    interface maxval
        module procedure maxval_d ! maxval of a dual number vector
    end interface

    public min
    interface min
        module procedure min_dd ! min of from two to four dual numbers, elemental
        module procedure min_dr ! min of a dual and a real, elemental
    end interface

    public dmin1
    interface dmin1
        module procedure min_dd ! min of from two to four dual numbers, elemental
    end interface

    public minval
    interface minval
        module procedure minval_d ! obtain the maxval  of a dual number vectgor
    end interface

    public nint
    interface nint
        module procedure nint_d ! nearest integer to the argument, elemental
    end interface

    public sign
    interface  sign
      module procedure  sign_dd ! sign(a,b) with two dual numbers, elemental
      module procedure  sign_rd ! sign(a,b) with a real and a dual, elemental
    end interface

    public sin
    interface sin
        module procedure sin_d ! obtain sine of a dual number, elemental
    end interface

    public dsin
    interface dsin
        module procedure sin_d ! obtain sine of a dual number, elemental
    end interface

    public tan
    interface tan
        module procedure tan_d ! obtain sine of a dual number, elemental
    end interface

    public dtan
    interface dtan
        module procedure tan_d ! obtain sine of a dual number, elemental
    end interface

    public sqrt
    interface sqrt
        module procedure sqrt_d ! obtain the sqrt of a dual number, elemental
    end interface

    public sum
    interface sum
        module procedure sum_d ! sum a dual array
    end interface

    public maxloc
    interface maxloc
        module procedure maxloc_d ! location of max in a dual array
    end interface

    public sinh
    interface sinh
        module procedure sinh_d ! obtain sinh of a dual number, elemental
    end interface

    public cosh
    interface cosh
        module procedure cosh_d ! obtain cosh of a dual number, elemental
    end interface

    public tanh
    interface tanh
        module procedure tanh_d ! obtain tanh of a dual number, elemental
    end interface

    public asinh
    interface asinh
        module procedure asinh_d ! obtain asinh of a dual number, elemental
    end interface

    public acosh
    interface acosh
        module procedure acosh_d ! obtain acosh of a dual number, elemental
    end interface

    public atanh
    interface atanh
        module procedure atanh_d ! obtain atanh of a dual number, elemental
    end interface

contains

!*********Begin: functions/subroutines for overloading operators

!******* Begin: (=)
!---------------------

    !-----------------------------------------
    ! dual = integer
    ! <u, du> = <i, 0>
    !-----------------------------------------
    elemental subroutine assign_di(u, i)
         type(dual), intent(out) :: u
         integer, intent(in) :: i

         u%x = real(i)  ! This is faster than direct assignment
         u%dx = 0.0

    end subroutine assign_di


    !-----------------------------------------
    ! dual = real(double)
    ! <u, du> = <r, 0>
    !-----------------------------------------
    elemental subroutine assign_dr(u, r)
        type(dual), intent(out) :: u
        real, intent(in) :: r

        u%x = r
        u%dx = 0.0

    end subroutine assign_dr


    !-----------------------------------------
    ! integer = dual
    ! i = <u, du>
    !-----------------------------------------
    elemental subroutine assign_id(i, v)
         type(dual), intent(in) :: v
         integer, intent(out) :: i

         i = int(v%x)

    end subroutine assign_id

!******* end: (=)
!---------------------


!******* Begin: (+)
!---------------------

    !-----------------------------------------
    ! Unary positive
    ! <res, dres> = +<u, du>
    !-----------------------------------------
    elemental function add_d(u) result(res)
         type(dual), intent(in) :: u
         type(dual) :: res

         res = u  ! Faster than assigning component wise

    end function add_d


    !-----------------------------------------
    ! dual + dual
    ! <res, dres> = <u, du> + <v, dv> = <u + v, du + dv>
    !-----------------------------------------
    elemental function add_dd(u, v) result(res)
         type(dual), intent(in) :: u, v
         type(dual) :: res

         res%x = u%x + v%x
         res%dx = u%dx + v%dx

    end function add_dd


    !-----------------------------------------
    ! dual + integer
    ! <res, dres> = <u, du> + i = <u + i, du>
    !-----------------------------------------
    elemental function add_di(u, i) result(res)
         type(dual), intent(in) :: u
         integer, intent(in) :: i
         type(dual) :: res

         res%x = real(i) + u%x
         res%dx = u%dx

    end function add_di


    !-----------------------------------------
    ! dual + double
    ! <res, dres> = <u, du> + <r, 0> = <u + r, du>
    !-----------------------------------------
    elemental function add_dr(u, r) result(res)
        type(dual), intent(in) :: u
        real, intent(in) :: r
        type(dual) :: res

        res%x = r + u%x
        res%dx = u%dx

    end function add_dr


    !-----------------------------------------
    ! integer + dual
    ! <res, dres> = <i, 0> + <v, dv> = <i + v, dv>
    !-----------------------------------------
    elemental function add_id(i, v) result(res)
        integer, intent(in) :: i
        type(dual), intent(in) :: v
        type(dual) :: res

        res%x = real(i) + v%x
        res%dx = v%dx

    end function add_id


    !-----------------------------------------
    ! double + dual
    ! <res, dres> = <r, 0> + <v, dv> = <r + v, dv>
    !-----------------------------------------
    elemental function add_rd(r, v) result(res)
        real, intent(in) :: r
        type(dual), intent(in) :: v
        type(dual) :: res

        res%x = r + v%x
        res%dx = v%dx

    end function add_rd

!******* end: (+)
!---------------------


!******* Begin: (-)
!---------------------

    !-------------------------------------------------
    ! negate a dual
    ! <res, dres> = -<u, du>
    !-------------------------------------------------
    elemental function minus_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = -u%x
        res%dx = -u%dx

    end function minus_d


    !-------------------------------------------------
    ! dual - dual
    ! <res, dres> = <u, du> - <v, dv> = <u - v, du - dv>
    !-------------------------------------------------
    elemental function minus_dd(u, v) result(res)
        type(dual), intent(in) :: u, v
        type(dual) :: res

        res%x = u%x - v%x
        res%dx = u%dx - v%dx

    end function minus_dd

    !-------------------------------------------------
    ! dual - integer
    ! <res, dres> = <u, du> - i = <u - i, du>
    !-------------------------------------------------
    elemental function minus_di(u, i) result(res)
        type(dual), intent(in) :: u
        integer, intent(in) :: i
        type(dual) :: res

        res%x = u%x - real(i)
        res%dx = u%dx

    end function minus_di


    !-------------------------------------------------
    ! dual - double
    ! <res, dres> = <u, du> - r = <u - r, du>
    !-------------------------------------------------
    elemental function minus_dr(u, r) result(res)
        type(dual), intent(in) :: u
        real,intent(in) :: r
        type(dual) :: res

        res%x = u%x - r
        res%dx = u%dx

    end function minus_dr


    !-------------------------------------------------
    ! integer - dual
    ! <res, dres> = i - <v, dv> = <i - v, -dv>
    !-------------------------------------------------
    elemental function minus_id(i, v) result(res)
        integer, intent(in) :: i
        type(dual), intent(in) :: v
        type(dual) :: res

        res%x = real(i) - v%x
        res%dx = -v%dx

    end function minus_id


    !-------------------------------------------------
    ! double - dual
    ! <res, dres> = r - <v, dv> = <r - v, -dv>
    !-------------------------------------------------
    elemental function minus_rd(r, v) result(res)
         real, intent(in) :: r
         type(dual), intent(in) :: v
         type(dual) :: res

        res%x = r - v%x
        res%dx = -v%dx

    end function minus_rd

!******* end: (-)
!---------------------


!******* BEGIN: (*)
!---------------------

    !----------------------------------------
    ! dual * dual
    ! <res, dres> = <u, du> * <v, dv> = <u * v, u * dv + v * du>
    !----------------------------------------
    elemental function mult_dd(u, v) result(res)
        type(dual), intent(in) :: u, v
        type(dual) :: res

        res%x = u%x * v%x
        res%dx = u%x * v%dx + v%x * u%dx

    end function mult_dd


    !-----------------------------------------
    ! dual * integer
    ! <res, dres> = <u, du> * i = <u * i, du * i>
    !-----------------------------------------
    elemental function mult_di(u, i) result(res)
        type(dual), intent(in) :: u
        integer, intent(in) :: i
        type(dual) :: res

        real :: r

        r = real(i)
        res%x = r * u%x
        res%dx = r * u%dx

    end function mult_di

    !-----------------------------------------
    ! dual * double
    ! <res, dres> = <u, du> * r = <u * r, du * r>
    !----------------------------------------
    elemental function mult_dr(u, r) result(res)
        type(dual), intent(in) :: u
        real, intent(in) :: r
        type(dual) :: res

        res%x = u%x * r
        res%dx = u%dx * r

    end function mult_dr


    !-----------------------------------------
    ! integer * dual
    ! <res, dres> = i * <v, dv> = <i * v, i * dv>
    !-----------------------------------------
    elemental function mult_id(i, v) result(res)
        integer, intent(in) :: i
        type(dual), intent(in) :: v
        type(dual) :: res

        real :: r

        r = real(i)
        res%x = r * v%x
        res%dx = r * v%dx

    end function mult_id


    !-----------------------------------------
    ! double * dual
    ! <res, dres> = r * <v, dv> = <r * v, r * dv>
    !-----------------------------------------
    elemental function mult_rd(r, v) result(res)
        real, intent(in) :: r
        type(dual), intent(in) :: v
        type(dual) :: res

        res%x = r * v%x
        res%dx = r * v%dx

    end function mult_rd

!******* end: (*)
!---------------------


!******* BEGIN: (/)
!---------------------

    !-----------------------------------------
    ! dual / dual
    ! <res, dres> = <u, du> / <v, dv> = <u / v, du / v - u * dv / v^2>
    !-----------------------------------------
    elemental function div_dd(u, v) result(res)
        type(dual), intent(in) :: u, v
        type(dual) :: res

        real :: inv

        inv = 1.0 / v%x
        res%x = u%x * inv
        res%dx = (u%dx - res%x * v%dx) * inv

    end function div_dd


    !-----------------------------------------
    ! dual / integer
    ! <res, dres> = <u, du> / i = <u / i, du / i>
    !-----------------------------------------
    elemental function div_di(u, i) result(res)
        type(dual), intent(in) :: u
        integer, intent(in) :: i
        type(dual) :: res

        real :: inv

        inv = 1.0 / real(i)
        res%x = u%x * inv
        res%dx = u%dx * inv

    end function div_di


    !-----------------------------------------
    ! dual / double
    ! <res, dres> = <u, du> / r = <u / r, du / r>
    !----------------------------------------
    elemental function div_dr(u, r) result(res)
        type(dual), intent(in) :: u
        real, intent(in) :: r
        type(dual):: res

        real :: inv

        inv = 1.0 / r
        res%x = u%x * inv
        res%dx = u%dx * inv

    end function div_dr


    !-----------------------------------------
    ! integer / dual
    ! <res, dres> = i / <v, dv> = <i / v, -i / v^2 * du>
    !-----------------------------------------
    elemental function div_id(i, v) result(res)
        integer, intent(in) :: i
        type(dual), intent(in) :: v
        type(dual) :: res

        real :: inv

        inv = 1.0 / v%x
        res%x = real(i) * inv
        res%dx = -res%x * inv * v%dx

    end function div_id


    !-----------------------------------------
    ! double / dual
    ! <res, dres> = r / <u, du> = <r / u, -r / u^2 * du>
    !-----------------------------------------
    elemental function div_rd(r, v) result(res)
        real, intent(in) :: r
        type(dual), intent(in) :: v
        type(dual) :: res

        real :: inv

        inv = 1.0 / v%x
        res%x = r * inv
        res%dx = -res%x * inv * v%dx

    end function div_rd

!******* end: (/)
!---------------------

!******* BEGIN: (**)
!---------------------

    !-----------------------------------------
    ! power(dual, integer)
    ! <res, dres> = <u, du> ^ i = <u ^ i, i * u ^ (i - 1) * du>
    !-----------------------------------------
    elemental function pow_i(u, i) result(res)
        type(dual), intent(in) :: u
        integer, intent(in) :: i
        type(dual) :: res

        real :: pow_x

        pow_x = u%x ** (i - 1)
        res%x = u%x * pow_x
        res%dx = real(i) * pow_x * u%dx

    end function pow_i

    !-----------------------------------------
    ! power(dual, double)
    ! <res, dres> = <u, du> ^ r = <u ^ r, r * u ^ (r - 1) * du>
    !-----------------------------------------
    elemental function pow_r(u, r) result(res)
        type(dual), intent(in) :: u
        real, intent(in) :: r
        type(dual) :: res

        real :: pow_x

        pow_x = u%x ** (r - 1.0)
        res%x = u%x * pow_x
        res%dx = r * pow_x * u%dx

    end function pow_r

    !-----------------------------------------
    ! POWER dual numbers to a dual power
    ! <res, dres> = <u, du> ^ <v, dv>
    !     = <u ^ v, u ^ v * (v / u * du + Log(u) * dv)>
    !-----------------------------------------
    elemental function pow_d(u, v) result(res)
        type(dual), intent(in)::u, v
        type(dual) :: res

        res%x = u%x ** v%x
        res%dx = res%x * (v%x / u%x * u%dx + log(u%x) * v%dx)

    end function pow_d

!******* end: (**)
!---------------------


!******* BEGIN: (==)
!---------------------
    !-----------------------------------------
    ! compare two dual numbers,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function eq_dd(lhs, rhs) result(res)
         type(dual), intent(in) :: lhs, rhs
         logical :: res

         res = (lhs%x == rhs%x)

    end function eq_dd


    !-----------------------------------------
    ! compare a dual with an integer,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function eq_di(lhs, rhs) result(res)
         type(dual), intent(in) :: lhs
         integer, intent(in) :: rhs
         logical :: res

         res = (lhs%x == real(rhs))

    end function eq_di


    !-----------------------------------------
    ! compare a dual number with a real number,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function eq_dr(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        real, intent(in) :: rhs
        logical::res

        res = (lhs%x == rhs)

    end function eq_dr


    !-----------------------------------------
    ! compare an integer with a dual,
    ! simply compare the functional value.
    !----------------------------------------
    elemental function eq_id(lhs, rhs) result(res)
         integer, intent(in) :: lhs
         type(dual), intent(in) :: rhs
         logical :: res

         res = (lhs == rhs%x)

    end function eq_id


    !-----------------------------------------
    ! compare a real with a dual,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function eq_rd(lhs, rhs) result(res)
         real, intent(in) :: lhs
         type(dual), intent(in) :: rhs
         logical :: res

         res = (lhs == rhs%x)

    end function eq_rd

!******* end: (==)
!---------------------


!******* BEGIN: (<=)
!---------------------
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function le_dd(lhs, rhs) result(res)
         type(dual), intent(in) :: lhs, rhs
         logical :: res

         res = (lhs%x <= rhs%x)

    end function le_dd


    !-----------------------------------------
    ! compare a dual with an integer,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function le_di(lhs, rhs) result(res)
         type(dual), intent(in) :: lhs
         integer, intent(in) :: rhs
         logical :: res

         res = (lhs%x <= rhs)

    end function le_di


    !-----------------------------------------
    ! compare a dual number with a real number,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function le_dr(lhs, rhs) result(res)
         type(dual), intent(in) :: lhs
         real, intent(in) :: rhs
         logical :: res

         res = (lhs%x <= rhs)

    end function le_dr


    !-----------------------------------------
    ! compare a dual number with an integer,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function le_id(i, rhs) result(res)
         integer, intent(in) :: i
         type(dual), intent(in) :: rhs
         logical :: res

         res = (i <= rhs%x)

    end function le_id


    !-----------------------------------------
    ! compare a real with a dual,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function le_rd(lhs, rhs) result(res)
         real, intent(in) :: lhs
         type(dual), intent(in) :: rhs
         logical :: res

         res = (lhs <= rhs%x)

    end function le_rd

!******* end: (<=)
!---------------------

!******* BEGIN: (<)
!---------------------
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function lt_dd(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs, rhs
        logical :: res

        res = (lhs%x < rhs%x)

    end function lt_dd

    !-----------------------------------------
    ! compare a dual with an integer,
    ! simply compare the functional value.
    !-----------------------------------------
    elemental function lt_di(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        integer, intent(in) :: rhs
        logical :: res

        res = (lhs%x < rhs)

    end function lt_di


    !-----------------------------------------
    ! compare a dual number with a real number, simply compare
    ! the functional value.
    !----------------------------------------
    elemental function lt_dr(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        real, intent(in) :: rhs
        logical :: res

        res = (lhs%x < rhs)

    end function lt_dr


    !-----------------------------------------
    ! compare a dual number with an integer
    !-----------------------------------------
    elemental function lt_id(i, rhs) result(res)
         integer, intent(in) :: i
         type(dual), intent(in) :: rhs
         logical :: res

         res = (i < rhs%x)

    end function lt_id


    !-----------------------------------------
    ! compare a real with a dual
    !----------------------------------------
    elemental function lt_rd(lhs, rhs) result(res)
         real, intent(in) :: lhs
         type(dual), intent(in) :: rhs
         logical :: res

         res = (lhs < rhs%x)

    end function lt_rd

!******* end: (<)
!---------------------

!******* BEGIN: (>=)
!---------------------
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !----------------------------------------
    elemental function ge_dd(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs, rhs
        logical :: res

        res = (lhs%x >= rhs%x)

    end function ge_dd


    !-----------------------------------------
    ! compare a dual with an integer
    !-----------------------------------------
    elemental function ge_di(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        integer, intent(in) :: rhs
        logical :: res

        res = (lhs%x >= rhs)

    end function ge_di


    !-----------------------------------------
    ! compare a dual number with a real number, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function ge_dr(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        real, intent(in) :: rhs
        logical :: res

        res = (lhs%x >= rhs)

    end function ge_dr


    !-----------------------------------------
    ! compare a dual number with an integer
    !-----------------------------------------
    elemental function ge_id(i, rhs) result(res)
        integer, intent(in) :: i
        type(dual), intent(in) :: rhs
        logical :: res

        res = (i >= rhs%x)

    end function ge_id


    !-----------------------------------------
    ! compare a real with a dual
    !-----------------------------------------
    elemental function ge_rd(lhs, rhs) result(res)
         real, intent(in) :: lhs
         type(dual), intent(in) :: rhs
         logical :: res

         res = (lhs >= rhs%x)

    end function ge_rd

!******* end: (>=)
!---------------------

!******* BEGIN: (>)
!---------------------
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function gt_dd(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs, rhs
        logical :: res

        res = (lhs%x > rhs%x)

    end function gt_dd


    !-----------------------------------------
    ! compare a dual with an integer
    !-----------------------------------------
    elemental function gt_di(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        integer, intent(in) :: rhs
        logical :: res

        res = (lhs%x > rhs)

    end function gt_di


    !-----------------------------------------
    ! compare a dual number with a real number, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function gt_dr(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        real, intent(in) :: rhs
        logical :: res

        res = (lhs%x > rhs)

    end function gt_dr


    !-----------------------------------------
    ! compare a dual number with an integer
    !-----------------------------------------
    elemental function gt_id(i, rhs) result(res)
        integer, intent(in) :: i
        type(dual), intent(in) :: rhs
        logical :: res

        res = (i > rhs%x)

    end function gt_id


    !-----------------------------------------
    ! compare a real with a dual
    !-----------------------------------------
    elemental function gt_rd(lhs, rhs) result(res)
         real, intent(in) :: lhs
         type(dual), intent(in) :: rhs
         logical :: res

         res = (lhs > rhs%x)

    end function gt_rd

!******* end: (>)
!---------------------

!******* BEGIN: (/=)
!---------------------
    !-----------------------------------------
    ! compare two dual numbers, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function ne_dd(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs, rhs
        logical :: res

        res = (lhs%x /= rhs%x)

    end function ne_dd


    !-----------------------------------------
    ! compare a dual with an integer
    !-----------------------------------------
    elemental function ne_di(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        integer, intent(in) :: rhs
        logical :: res

        res = (lhs%x /= rhs)

    end function ne_di


    !-----------------------------------------
    ! compare a dual number with a real number, simply compare
    ! the functional value.
    !-----------------------------------------
    elemental function ne_dr(lhs, rhs) result(res)
        type(dual), intent(in) :: lhs
        real, intent(in) :: rhs
        logical :: res

        res = (lhs%x /= rhs)

    end function ne_dr


    !-----------------------------------------
    ! compare a dual number with an integer
    !-----------------------------------------
    elemental function ne_id(i, rhs) result(res)
        integer, intent(in) :: i
        type(dual), intent(in) :: rhs
        logical :: res

        res = (i /= rhs%x)

    end function ne_id


    !-----------------------------------------
    ! compare a real with a dual
    !-----------------------------------------
    elemental function ne_rd(lhs, rhs) result(res)
        real, intent(in) :: lhs
        type(dual), intent(in) :: rhs
        logical :: res

        res = (lhs /= rhs%x)

    end function ne_rd

!******* end: (/=)
!---------------------

    !---------------------------------------------------
    ! Absolute value of dual numbers
    ! <res, dres> = abs(<u, du>) = <abs(u), du * sign(u)>
    !---------------------------------------------------
    elemental function abs_d(u) result(res)
         type(dual), intent(in) :: u
         type(dual) :: res
         integer :: i

         if(u%x > 0) then
            res%x = u%x
            res%dx = u%dx
         else if (u%x < 0) then
            res%x = -u%x
            res%dx = -u%dx
         else
            res%x = 0.0
            do i = 1, ndv
                if (u%dx(i) .eq. 0.0) then
                    res%dx(i) = 0.0
                else
                    res%dx(i) = set_NaN()
                end if
            end do
         endif

    end function abs_d


    !-----------------------------------------
    ! ACOS of dual numbers
    ! <res, dres> = acos(<u, du>) = <acos(u), -du / sqrt(1 - u^2)>
    !----------------------------------------
    elemental function acos_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = acos(u%x)
        if (u%x == 1.0 .or. u%x == -1.0) then
            res%dx = set_Nan()  ! Undefined derivative
        else
            res%dx = -u%dx / sqrt(1.0 - u%x**2)
        end if

    end function acos_d


    !-----------------------------------------
    ! ASIN of dual numbers
    ! <res, dres> = asin(<u, du>) = <asin(u), du / sqrt(1 - u^2)>
    !----------------------------------------
    elemental function asin_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = asin(u%x)
        if (u%x == 1.0 .or. u%x == -1.0) then
            res%dx = set_NaN()  ! Undefined derivative
        else
            res%dx = u%dx / sqrt(1.0 - u%x**2)
        end if

    end function asin_d


    !-----------------------------------------
    ! ATAN of dual numbers
    ! <res, dres> = atan(<u, du>) = <atan(u), du / (1 + u^2)>
    !----------------------------------------
    elemental function atan_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = atan(u%x)
        res%dx = u%dx / (1.0 + u%x**2)

    end function atan_d


    !-----------------------------------------
    ! ATAN2 of dual numbers
    ! <res, dres> = atan2(<u, du>, <v, dv>)
    !             = <atan2(u, v), v / (u^2 + v^2) * du - u / (u^2 + v^2) * dv>
    !----------------------------------------
    elemental function atan2_d(u, v) result(res)
        type(dual), intent(in) :: u, v
        type(dual) :: res

        real :: usq_plus_vsq

        res%x = atan2(u%x, v%x)

        usq_plus_vsq = u%x**2 + v%x**2
        res%dx = v%x / usq_plus_vsq * u%dx - u%x / usq_plus_vsq * v%dx

    end function atan2_d


    !-----------------------------------------
    ! COS of dual numbers
    ! <res, dres> = cos(<u, du>) = <cos(u), -sin(u) * du>
    !----------------------------------------
    elemental function cos_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = cos(u%x)
        res%dx = -sin(u%x) * u%dx

    end function cos_d


    !-----------------------------------------
    ! DOT PRODUCT two dual number vectors
    ! <res, dres> = <u, du> . <v, dv> = <u . v, u . dv + v . du>
    !-----------------------------------------
    function dot_product_dd(u, v) result(res)
        type(dual), intent(in) :: u(:), v(:)
        type(dual) :: res

        integer :: i

        res%x = dot_product(u%x, v%x)
        do i = 1, ndv
            res%dx(i) = dot_product(u%x, v%dx(i)) + dot_product(v%x, u%dx(i))
        end do

    end function dot_product_dd


    !-----------------------------------------
    ! EXPONENTIAL OF dual numbers
    ! <res, dres> = exp(<u, du>) = <exp(u), exp(u) * du>
    !-----------------------------------------
    elemental function exp_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        real :: exp_x

        exp_x = exp(u%x)
        res%x = exp_x
        res%dx = u%dx * exp_x

    end function exp_d


    !-----------------------------------------
    ! Convert dual to integer
    ! i = int(<u, du>) = int(u)
    !----------------------------------------
    elemental function int_d(u) result(res)
         type(dual), intent(in) :: u
         integer :: res

         res = int(u%x)

    end function int_d


    !-----------------------------------------
    ! LOG OF dual numbers,defined for u%x>0 only
    ! the error control should be done in the original code
    ! in other words, if u%x<=0, it is not possible to obtain LOG.
    ! <res, dres> = log(<u, du>) = <log(u), du / u>
    !----------------------------------------
    elemental function log_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        real :: inv

        inv = 1.0 / u%x
        res%x = log(u%x)
        res%dx = u%dx * inv

    end function log_d


    !-----------------------------------------
    ! LOG10 OF dual numbers,defined for u%x>0 only
    ! the error control should be done in the original code
    ! in other words, if u%x<=0, it is not possible to obtain LOG.
    ! <res, dres> = log10(<u, du>) = <log10(u), du / (u * log(10))>
    ! LOG<u,up>=<LOG(u),up/u>
    !----------------------------------------
    elemental function log10_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        real :: inv

        inv = 1.0 / (u%x * log(10.0))
        res%x = log10(u%x)
        res%dx = u%dx * inv

    end function log10_d


    !-----------------------------------------
    ! MULTIPLY two dual number matrices
    ! <res, dres> = <u, du> . <v, dv> = <u . v, du . v + u . dv>
    !----------------------------------------
    function matmul_dd(u,v) result(res)
        type(dual), intent(in) :: u(:,:), v(:,:)
        type(dual) :: res(size(u,1), size(v,2))

        integer :: i

        res%x = matmul(u%x, v%x)
        do i = 1, ndv
            res%dx(i) = matmul(u%dx(i), v%x) + matmul(u%x, v%dx(i))
        end do

    end function matmul_dd


    !-----------------------------------------
    ! MULTIPLY a dual number matrix with a dual number
    ! vector
    !
    ! <u,up>.<v,vp>=<u.v,up.v+u.vp>
    !----------------------------------------
    function matmul_dv(u, v) result(res)
        type(dual), intent(in) :: u(:,:), v(:)
        type(dual) :: res(size(u,1))
        integer :: i

        res%x = matmul(u%x, v%x)
        do i = 1, ndv
            res%dx(i) = matmul(u%dx(i), v%x) + matmul(u%x, v%dx(i))
        end do

    end function matmul_dv


    !-----------------------------------------
    ! MULTIPLY a dual vector with a  dual matrix
    !
    ! <u,up>.<v,vp>=<u.v,up.v+u.vp>
    !----------------------------------------
    function matmul_vd(u, v) result(res)
        type(dual), intent(in) :: u(:), v(:,:)
        type(dual) :: res(size(v, 2))
        integer::i

        res%x = matmul(u%x, v%x)
        do i = 1, ndv
            res%dx(i) = matmul(u%dx(i), v%x) + matmul(u%x, v%dx(i))
        end do

    end function matmul_vd

    !-----------------------------------------
    ! Obtain the max of 2 to 5 dual numbers
    !----------------------------------------
    elemental function max_dd(val1, val2, val3, val4,val5) result(res)
        type(dual), intent(in) :: val1, val2
        type(dual), intent(in), optional :: val3, val4,val5
        type(dual) :: res

        if (val1%x > val2%x) then
            res = val1
        else
            res = val2
        endif
        if(present(val3))then
           if(res%x < val3%x) res = val3
        endif
        if(present(val4))then
           if(res%x < val4%x) res = val4
        endif
        if(present(val5))then
           if(res%x < val5%x) res = val5
        endif

    end function max_dd


    !-----------------------------------------
    ! Obtain the max of a dual number and an integer
    !----------------------------------------
    elemental function max_di(u, i) result(res)
        type(dual), intent(in) :: u
        integer, intent(in) :: i
        type(dual) :: res

        if (u%x > i) then
            res = u
        else
            res = i
        endif

    end function max_di

    !-----------------------------------------
    ! Obtain the max of a dual number and a real number
    !----------------------------------------
    elemental function max_dr(u, r) result(res)
        type(dual), intent(in) :: u
        real, intent(in) :: r
        type(dual) :: res

        if (u%x > r) then
            res = u
        else
            res = r
        endif

    end function max_dr


    !---------------------------------------------------
    ! Obtain the max of a real and a dual
    !---------------------------------------------------
     elemental function max_rd(n, u) result(res)
        real, intent(in) :: n
        type(dual), intent(in) :: u
        type(dual) :: res

        if (u%x > n) then
            res = u
        else
            res = n
        endif

    end function max_rd


    !-----------------------------------------
    ! Obtain the max value of vector u
    !----------------------------------------
    function maxval_d(u) result(res)
        type(dual), intent(in) :: u(:)
        integer :: iloc(1)
        type(dual) :: res

        iloc=maxloc(u%x)
        res=u(iloc(1))

    end function maxval_d


    !-----------------------------------------
    ! Obtain the min of 2 to 4 dual numbers
    !----------------------------------------
    elemental function min_dd(val1, val2, val3, val4) result(res)
        type(dual), intent(in) :: val1, val2
        type(dual), intent(in), optional :: val3, val4
        type(dual) :: res

        if (val1%x < val2%x) then
            res = val1
        else
            res = val2
        endif
        if(present(val3))then
           if(res%x > val3%x) res = val3
        endif
        if(present(val4))then
           if(res%x > val4%x) res = val4
        endif

    end function min_dd


    !-----------------------------------------
    ! Obtain the min of a dual and a double
    !----------------------------------------
    elemental function min_dr(u, r) result(res)
        type(dual), intent(in) :: u
        real, intent(in) :: r
        type(dual) :: res

        if (u%x < r) then
            res = u
        else
            res = r
        endif

    end function min_dr


  !-----------------------------------------
    ! Obtain the min value of vector u
    !----------------------------------------
    function minval_d(u) result(res)
        type(dual), intent(in) :: u(:)
        integer :: iloc(1)
        type(dual) :: res

        iloc=minloc(u%x)
        res=u(iloc(1))

    end function minval_d


    !------------------------------------------------------
    !Returns the nearest integer to u%x, ELEMENTAL
    !------------------------------------------------------
    elemental function nint_d(u) result(res)
        type(dual), intent(in) :: u
        integer :: res

        res=nint(u%x)

    end function nint_d


    !----------------------------------------------------------------
    ! SIGN(a,b) with two dual numbers as inputs,
    ! the result will be |a| if b%x>=0, -|a| if b%x<0,ELEMENTAL
    !----------------------------------------------------------------
    elemental function sign_dd(val1, val2) result(res)
        type(dual), intent(in) :: val1, val2
        type(dual) :: res

        if (val2%x < 0.0) then
            res = -abs(val1)
        else
            res =  abs(val1)
        endif

     end function sign_dd


    !----------------------------------------------------------------
    ! SIGN(a,b) with one real and one dual number as inputs,
    ! the result will be |a| if b%x>=0, -|a| if b%x<0,ELEMENTAL
    !----------------------------------------------------------------
    elemental function sign_rd(val1, val2) result(res)
        real, intent(in) :: val1
        type(dual), intent(in) :: val2
        type(dual) :: res

        if (val2%x < 0.0) then
            res = -abs(val1)
        else
            res = abs(val1)
        endif

     end function sign_rd


    !-----------------------------------------
    ! SIN of dual numbers
    ! <res, dres> = sin(<u, du>) = <sin(u), cos(u) * du>
    !----------------------------------------
    elemental function sin_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = sin(u%x)
        res%dx = cos(u%x) * u%dx

    end function sin_d


    !-----------------------------------------
    ! TAN of dual numbers
    ! <res, dres> = tan(<u, du>) = <tan(u), du / cos(u)^2>
    !----------------------------------------
    elemental function tan_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = tan(u%x)
        res%dx = u%dx / cos(u%x)**2

    end function tan_d


    !-----------------------------------------
    ! SQRT of dual numbers
    ! <res, dres> = sqrt(<u, du>) = <sqrt(u), du / (2 * sqrt(u))>
    !----------------------------------------
    elemental function sqrt_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res
        integer :: i

        res%x = sqrt(u%x)

        if (res%x .ne. 0.0) then
            res%dx = 0.5 * u%dx / res%x
        else
            do i = 1, ndv
                if (u%dx(i) .eq. 0.0) then
                    res%dx(i) = 0.0
                else
                    res%dx(i) = set_NaN()
                end if
            end do
        end if

    end function sqrt_d


    !-----------------------------------------
    ! Sum of a dual array
    !-----------------------------------------
    function sum_d(u) result(res)
        type(dual), intent(in) :: u(:)
        type(dual) :: res
        integer :: i

        res%x = sum(u%x)
        do i = 1, ndv
            res%dx(i) = sum(u%dx(i))
        end do

    end function sum_d


    !-----------------------------------------
    ! Find the location of the max value in an
    ! array of dual numbers
    !-----------------------------------------
    function maxloc_d(array) result(ind)
        type(dual), intent(in) :: array(:)
        integer :: ind(1)

        ind = maxloc(array%x)

    end function maxloc_d


    elemental function set_NaN() result(res)
        real :: res

        res = sqrt(negative_one)

    end function set_NaN


    !-----------------------------------------
    ! Hyperbolic functions: sinh, cosh, tanh
    ! and their inverses: asinh, acosh, atanh
    !-----------------------------------------
    !-----------------------------------------
    ! SINH OF dual numbers
    ! <res, dres> = sinh(<u, du>) = <sinh(u), cosh(u) * du>
    !-----------------------------------------
    elemental function sinh_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = sinh(u%x)
        res%dx = u%dx * cosh(u%x)

    end function sinh_d

    !-----------------------------------------
    ! COSH OF dual numbers
    ! <res, dres> = cosh(<u, du>) = <cosh(u), sinh(u) * du>
    !-----------------------------------------
    elemental function cosh_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = cosh(u%x)
        res%dx = u%dx * sinh(u%x)

    end function cosh_d

    !-----------------------------------------
    ! TANH OF dual numbers
    ! <res, dres> = tanh(<u, du>) = <tanh(u), 1.0/cosh(u)**2 * du>
    !-----------------------------------------
    elemental function tanh_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = tanh(u%x)
        res%dx = u%dx * 1.0/cosh(u%x)**2

    end function tanh_d

    !-----------------------------------------
    ! ASINH OF dual numbers
    ! <res, dres> = asinh(<u, du>) = <asinh(u), 1/sqrt(u**2 + 1) * du>
    !-----------------------------------------
    elemental function asinh_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = asinh(u%x)
        res%dx = u%dx * 1.0/sqrt(u%x**2 + 1.0)

    end function asinh_d

    !-----------------------------------------
    ! ACOSH OF dual numbers
    ! <res, dres> = acosh(<u, du>) = <acosh(u), 1/sqrt(u**2 - 1) * du>
    !-----------------------------------------
    elemental function acosh_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = acosh(u%x)
        if (u%x <= 1.0) then
            res%dx = set_Nan()  ! Undefined derivative
        else
            res%dx = u%dx * 1.0/sqrt(u%x**2 - 1.0)
        end if

    end function acosh_d

    !-----------------------------------------
    ! ATAHN OF dual numbers
    ! <res, dres> = atanh(<u, du>) = <atanh(u), 1/(1 - u**2) * du>
    !-----------------------------------------
    elemental function atanh_d(u) result(res)
        type(dual), intent(in) :: u
        type(dual) :: res

        res%x = atanh(u%x)
        if (abs(u%x) >= 1.0) then
            res%dx = set_Nan()  ! Undefined derivative
        else
            res%dx = u%dx * 1.0/(1.0 - u%x**2)
        end if

    end function atanh_d


end module dnadmod
! ******************************************************************************
! ******************************************************************************
