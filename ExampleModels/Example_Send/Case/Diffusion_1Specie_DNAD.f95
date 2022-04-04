! * Written by Nicholas August 10, 2020
! Uses John Newman's Band algorithm to solve coupled non-linear partial
! differential equations
! Incorporates DNAD which anables to use of automatic Differentiation to
! linearize the PDE's
! variables or names that begin and end with '_', i.e. _variable_ are changed by the python program: RunFortran.py

module number_of_variables
  implicit none
  integer, parameter :: N = 1
  integer, parameter :: NJ = 22                          ! Number of mesh points
end module number_of_variables


! include '../Core_Subs_Mods/dnadmod.f95'
! include '../Core_Subs_Mods/variables.f95'


module user_input
  use number_of_variables
  implicit none

  integer, parameter :: Numbertimesteps = 1e4     ! Number of time steps
  real               :: delT = 1.0                       ! size of timestep [s]
  real               :: time                             ! [s]
  real               :: xmax = 1.0                      ! 500 um is 500e-4 cm

  real, parameter :: PI = 4.0 * ATAN(1.0)       ! pi - Geometric constant

  real :: diff  = 1e-1                       ! diffusion coefficient [cm^2/s]
  real :: cbulk = 0.001                       ! concentration [mol/cm3]

  character(len=65) :: geometry = 'Rectangular'   ! system Geometry: Rectangular, Cylindrical, Spherical

end module user_input


! ******************************************************************************
module write_data_mod
  use user_input
  use variables, only: cprev, delC, xx

  implicit none

  real :: c0

  character(len=65) :: header_fmt = '(2(A12,1X),   20(A15,1X)) '
  character(len=65) :: data_fmt   = '(2(F12.5,1X), 20(ES15.5,1X))'

  real :: t_write

contains

  subroutine t_write__Equiv__Li_Bal_Assignment
    t_write     = time / 3600.0

  end subroutine t_write__Equiv__Li_Bal_Assignment


  subroutine write_condition(it)
    integer :: it
    real :: last_write_time = 0.0
    real :: write_every_x_sec = 3600.0/10.0           ! 3600 s = 1 hour

    call t_write__Equiv__Li_Bal_Assignment()
    call write_to_screen(it)

    ! if (it == 1) then
    !   call write_to_screen(it)
    !   call write_positional_information_to_file(it)
    !   last_write_time = time
    ! else if ( (time - last_write_time).GE.write_every_x_sec ) then
    !   call write_to_screen(it)
    !   call write_positional_information_to_file(it)
    !   last_write_time = time
    ! end if
  end subroutine write_condition


  subroutine write_to_screen(it)
    integer :: it

    if (it.EQ.1) then       ! write the headers on the first entrance into write all voltage
      write(*, header_fmt) 'Time',  'Conc'
      write(*, header_fmt) 'hours', 'mol/L'
    end if

    c0 = cprev(1,NJ-1)

    write(*, data_fmt) t_write,    c0

  end subroutine write_to_screen


  subroutine write_positional_information_to_file(it)
    use variables, only: xx
    integer :: it, j

    open(56, file = 'Time_Conc_Position.txt', status = 'unknown')

    if (it.EQ.1) then
!           write the headers on the first entrance into write all voltage
      write(56, header_fmt) 'Time',  'Position', 'Conc'
      write(56, header_fmt) 'hours', 'um'      , 'mol/L'
                                                    !
    end if                                          !
                                                    !
    do j = 1, NJ                                    !
      c0 = cprev(1,j) * 1e3                         !
                                                    !
      write(56, data_fmt) t_write,    xx(j)*1e4,    c0

    end do

  end subroutine write_positional_information_to_file


end module write_data_mod



module GOV_EQNS
  use user_input
  use variables, only: xx, delX
  use dnadmod

  implicit none

  real, dimension(NJ)    :: Cntrl_Vol       ! Control Volume
  real, dimension(2, NJ) :: Crx_Area        ! dimesion = (2,NJ) because east-west
                                           ! cross-sectional area
contains

! ******************************************************************************
! **************************** GOVERNING EQUATIONS *****************************
! **************** Accumulation = FluxIn - FluxOut + Generation ****************
! ******************************************************************************
! can define intermediate variables with the functions to improve readability
! i.e c0 = c_vars_dual(1), dPhi2_dx = dcdx_vars_dual(2)

  function FLUX(c_vars_dual, dcdx_vars_dual) result(Flux_)
    ! (1)   N_0 = -D dc/dx
    type(dual), dimension(N)              :: Flux_
    type(dual), dimension(N), intent(in)  :: c_vars_dual, dcdx_vars_dual
    type(dual)                            :: diff_

    type(dual) :: c0
    type(dual) :: dc0dx

    c0       = c_vars_dual(1)
    dc0dx    = dcdx_vars_dual(1)

    Flux_(1) = -diff * dc0dx

  end function FLUX

  function RXN(c_vars_dual) result(Rxn_)
    ! (1)  no reaction
    type(dual), dimension(N)               :: Rxn_
    type(dual), dimension(N), intent(in)   :: c_vars_dual

    type(dual) :: c0

    c0    = c_vars_dual(1)

    Rxn_(1) = 0.0

  end function RXN

  function ACCUM(c_vars_dual) result(Accum_)
    ! (1) dc0/dt
    type(dual), dimension(N)              :: Accum_
    type(dual), dimension(N), intent(in)  :: c_vars_dual

    type(dual) :: c0

    c0    = c_vars_dual(1)

    Accum_(1) = c0/delT

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
  function Boundary_WEST (c_vars_dual, dcdx_vars_dual) result(BC_WEST_)
    ! (1)   c0 = cbulk
    type(dual), dimension(N)               :: BC_WEST_
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual

    type(dual) :: c0

    c0      = c_vars_dual(1)

    BC_WEST_(1) = c0 - 1.0

  end function Boundary_WEST

  function Boundary_EAST (c_vars_dual, dcdx_vars_dual) result(BC_EAST_)
    ! (1)   c0 = 0.1 * cbulk
    type(dual), dimension(N)               :: BC_EAST_
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual

    type(dual) :: c0

    c0    = c_vars_dual(1)

    BC_EAST_(1) = c0 - 0.0                                ! c0 = 0.1*cbulk

  end function Boundary_EAST
! ******************************************************************************


  ! ****************************************************************************
  ! ******* functions to convert concentration and dcdx variables to dual ******
  ! ----------------------------------------------------------------------------
  function c_to_dual(c_vars) result(c_dual)
    type(dual), dimension(N)       :: c_dual
    real, dimension(N), intent(in) :: c_vars
    real, dimension(N*2)            :: dx_array
    INTEGER :: ic

    do ic = 1, N
      dx_array = 0.0              ! set the dx_array to zero (all elements)
      dx_array(ic) = 1.0

      c_dual(ic) = dual(c_vars(ic), dx_array)
    end do
  end function c_to_dual


  function dcdx_to_dual(dcdx) result(dcdx_dual)
    type(dual), dimension(N)       :: dcdx_dual
    real, dimension(N), intent(in) :: dcdx
    real, dimension(N*2)            :: dx_array
    INTEGER :: ic

    do ic = 1, N
      dx_array = 0.0              ! set the dx_array to zero (all elements)
      dx_array(N+ic) = 1.0

      dcdx_dual(ic) = dual(dcdx(ic), dx_array)
    end do

  end function dcdx_to_dual
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

! include '../Core_Subs_Mods/BAND_MOD.f95'

! ******************************************************************************
! ******************************** MAIN PROGRAM ********************************
! ******************************************************************************

program unsteady
  use user_input
  use variables, only : cprev, delC
  use write_data_mod
  use BAND_mod

  implicit none
  integer :: t1, t2, clock_rate, clock_max
  integer :: it, j

  call system_clock(t1,clock_rate,clock_max)
  call initial_condition()

  do it = 1, N
      print '(22(F8.3,1X))', cprev(it, :)
  end do

  do it = 1, Numbertimesteps

    ! call write_condition(it)
                                        ! function to dynamically and integlligently change delT
    ! **************************************************************************
    ! ******************************** BOUND VAL *******************************
    ! **************************************************************************
      do j = 1, NJ
        call auto_fill(j)               ! These can be changed to functions
        call ABDGXY(j)                  ! because they each have a set of inputs
        call BAND(j)                    ! and a set of outputs

      end do

      cprev = cprev + delC              ! Update the dependent variables
      time  = time  + delT              ! Update the time

  ! ****************************************************************************

  end do

  do it = 1, N
    print '(22(F8.3,1X))', cprev(it, :)
  end do

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
  real    :: h
  integer :: j

  character(len=:), allocatable :: control_volume_input

  h = xmax/float(nj-2)
  do j=1,NJ

    if (j == 1) then
      xx(j) = 0.0
    else if (j == NJ) then
      xx(NJ) = xmax
    else
      xx(j) = xx(1) + h*float(j-1) - h/2.0
    end if

  end do

  delX(2:NJ-1) = h
  delx(1) = 0.0                      ! it is common in the finite volume
  delx(NJ) = 0.0                     ! algorithm to set the control
                                        ! volumes at the boundaries to zero
  cprev(1,:) = 1.0

  control_volume_input = trim(geometry)                   ! Define the control volume
  Cntrl_Vol = Control_Volume(control_volume_input)        ! size based on the
  Crx_Area = Cross_Sectional_Area(control_volume_input)   ! specified system geometry

  return                                                  ! is this necessary?
end subroutine initial_condition
