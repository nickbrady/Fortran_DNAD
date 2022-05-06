module user_input
  implicit none

  integer, parameter :: N = 2
  integer, parameter :: NJ = 1002                        ! Number of mesh points
  integer, parameter :: Numbertimesteps = 100 * 60 * 1e3       ! Number of time steps
  integer, parameter :: maxIterations = 1e6
  real               :: delT = 1.0e-1                   ! size of timestep [s]
  real               :: time = 0.0                      ! [s]
  logical            :: UPWIND = .FALSE.
  character(len=65)  :: direction = 'WestToEast'

  real, parameter    :: xmax = 1.0                      ! [cm] 1 cm

  real, parameter :: Rigc   = 8.314                     ! Ideal gas constant [J/(mol*K)]
  real, parameter :: Fconst = 96485                     ! Faraday's Constant [C/mol]
  real, parameter :: Temp   = 298
  real, parameter :: PI = 4.0 * ATAN(1.0)               ! pi - Geometric constant

  real, parameter :: cbulk_Li = 1e-3
  ! real, parameter :: trans_Li = 0.5

  real, parameter :: D_Li     = 2e-6
  real, parameter :: u_Li     = D_Li/(Rigc*Temp)
  real, parameter :: z_Li     = +1.0
  real, parameter :: nu_Li    = 1.0
  real, parameter :: s_Li     = +1.0

  real, parameter :: D_PF6     = 1e-6
  real, parameter :: u_PF6     = D_PF6/(Rigc*Temp)
  real, parameter :: z_PF6     = -1.0
  real, parameter :: nu_PF6    = 1.0
  real, parameter :: s_PF6     = 0.0

  real :: i_applied_cm2            ! 1 μA/cm2

  character(len=65) :: geometry = 'Rectangular'   ! system Geometry: Rectangular, Cylindrical, Spherical

end module user_input

!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/variables.f95'
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/dnadmod.f95'
!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*

module Electrochemical_Program_Mod
  use user_input
  implicit none

  real    :: rest_time = 30*60
  real    :: pulse_time = 60
  real    :: end_time_prev_cyc = 0.0
  integer :: cycle_n = 1
  integer :: cycle_loops = 10
  real    :: experiment_time = 100*3600.

contains
  ! Current Profile
  subroutine current_profile()
    experiment_time = cycle_loops*(rest_time+pulse_time)

    if ((time - end_time_prev_cyc) >= (rest_time + pulse_time)) then
      end_time_prev_cyc = time
      cycle_n = cycle_n + 1
      ! call initial_condition() ! reset the concentration profiles (mimic waiting a long time)
    end if

    if ((time - end_time_prev_cyc) < pulse_time) then
      i_applied_cm2 = 1.0e-3 * (cycle_n)*2

    else
      i_applied_cm2 = 0.0

    end if

  end subroutine current_profile
end module Electrochemical_Program_Mod

module write_data_mod
  use user_input
  use variables, only: cprev, delC, xx
  use Electrochemical_Program_Mod

  implicit none

  character(len=65) :: header_fmt = '(1(A12,1X),   1(A10,1X), 20(A15,1X)) '
  character(len=65) :: data_fmt   = '(1(F12.5,1X), 1(I10,1X), 20(ES20.10,1X))'

  real :: t_write

  real :: conc_Li, Phi_2, Delta_Phi_2

contains

  subroutine t_write__Equiv__Li_Bal_Assignment
    t_write     = time

  end subroutine t_write__Equiv__Li_Bal_Assignment


  subroutine write_condition(it)
    integer :: it
    real :: last_write_time = -1e6
    real :: write_every_x_sec = 1.0e0           ! 3600 s = 1 hour

    call t_write__Equiv__Li_Bal_Assignment

    ! if (time <= 360) then
    !   if ( (time - last_write_time) >= write_every_x_sec) then
    !     call write_to_screen(it)
    !     call write_positional_information_to_file(it)
    !
    !     last_write_time = time - delT*0.5       ! subtract half a time-step for numerical reasons to ensure data is written at
    !   end if                                    ! even intervals
    !
    ! else
      if ( ((time - end_time_prev_cyc) >= pulse_time).AND. &
          & ((time - end_time_prev_cyc) <= (1.0 + 0.5*delT + pulse_time)) ) then
        ! write every time step for first 1 second of relaxation
        call write_positional_information_to_file(it)

        last_write_time = time! - delT*0.5
      end if

      if ( (time - last_write_time) >= write_every_x_sec) then
        call write_to_screen(it)
        call write_positional_information_to_file(it)

        last_write_time = time - delT*0.5
      end if
    ! end if

  end subroutine write_condition


  subroutine write_to_screen(it)
    integer :: it
    real    :: conc_Li_east, conc_Li_west

    if (it == 0) then       ! write the headers on the first entrance into write all voltage
      write(*, header_fmt)    'Time', 'Cycle#', 'Current',  'c_west',  'c_east', 'ΔΦ_2'
      write(*, header_fmt) 'seconds',      '#',  'mA/cm2', 'mol/cm3', 'mol/cm3', 'V'
    end if

      conc_Li_west  = cprev(1,1)
      conc_Li_east  = cprev(1,NJ)
      Delta_Phi_2   = cprev(2,NJ) - cprev(2,1)

      write(*, data_fmt) t_write, cycle_n, i_applied_cm2, conc_Li_west, conc_Li_east, Delta_Phi_2


  end subroutine write_to_screen


  subroutine write_positional_information_to_file(it)
    use variables, only: xx
    integer :: it, j

    open(56, file = 'Time_Conc_Position.txt', status = 'unknown')

    if (it == 0) then
!           write the headers on the first entrance into write all voltage
      write(56, header_fmt) 'Time', 'Cycle#', 'Current', 'Position',    'Conc', 'Potential'
      write(56, header_fmt) 'seconds',   '#',  'mA/cm2',       'cm', 'mol/cm3', 'V'
                                                      !
    end if                                              !


    do j = 1, NJ, int(NJ/10)                                    !
      conc_Li    = cprev(1,j)
      Phi_2      = cprev(2,j)
                                                    !
      write(56, data_fmt) t_write, cycle_n, i_applied_cm2, xx(j), conc_Li, Phi_2

    end do

  end subroutine write_positional_information_to_file

end module write_data_mod
!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*

module GOV_EQNS
  use user_input
  use variables
  use dnadmod

  implicit none

  type(dual) :: c_Li, dcdx_Li
  type(dual) :: c_PF6, dcdx_PF6
  type(dual) :: Phi_2, dPhi2_dx
  real       :: position

contains

! ******************************************************************************
! **************************** GOVERNING EQUATIONS *****************************
! **************** Accumulation = FluxIn - FluxOut + Generation ****************
! ******************************************************************************
! can define intermediate variables with the functions to improve readability
! i.e c0 = c_vars_dual(1), dPhi2_dx = dcdx_vars_dual(2)

  function FLUX(c_vars_dual, dcdx_vars_dual)  result(Flux_)
    type(dual), dimension(N)              :: Flux_
    type(dual), dimension(N), intent(in)  :: c_vars_dual, dcdx_vars_dual

    c_Li        = c_vars_dual(1)
    dcdx_Li     = dcdx_vars_dual(1)
    c_PF6       = c_Li
    dcdx_PF6    = dcdx_Li

    Phi_2       = c_vars_dual(2)
    dPhi2_dx    = dcdx_vars_dual(2)

    Flux_(1) = -D_Li  * dcdx_Li  - z_Li  * u_Li  * c_Li  * Fconst * dPhi2_dx
    Flux_(2) = -D_PF6 * dcdx_PF6 - z_PF6 * u_PF6 * c_PF6 * Fconst * dPhi2_dx


  end function FLUX

  function RXN(c_vars_dual)                   result(Rxn_)
    type(dual), dimension(N)               :: Rxn_
    type(dual), dimension(N), intent(in)   :: c_vars_dual

    Rxn_(1) = 0.0
    Rxn_(2) = 0.0

  end function RXN

  function ACCUM(c_vars_dual)                           result(Accum_)
    type(dual), dimension(N)              :: Accum_
    type(dual), dimension(N), intent(in)  :: c_vars_dual

    c_Li        = c_vars_dual(1)
    c_PF6       = c_Li

    Accum_(1) = c_Li/delT
    Accum_(2) = c_PF6/delT

    Accum_%x = 0.0

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
    type(dual), dimension(N)               :: BC_WEST_, flux_temp, rxn_temp
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual

    c_Li        = c_vars_dual(1)
    dcdx_Li     = dcdx_vars_dual(1)

    Phi_2       = c_vars_dual(2)

    flux_temp = FLUX(c_vars_dual, dcdx_vars_dual)
    rxn_temp  = RXN(c_vars_dual)

    BC_WEST_(1) = z_Li*nu_Li*flux_temp(1) - i_applied_cm2 / Fconst
    BC_WEST_(2) = Phi_2 - 0.0

  end function Boundary_WEST

  function Boundary_EAST (c_vars_dual, dcdx_vars_dual) result(BC_EAST_)
    type(dual), dimension(N)               :: BC_EAST_, flux_temp, rxn_temp
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual

    c_Li        = c_vars_dual(1)
    dcdx_Li     = dcdx_vars_dual(1)

    flux_temp = FLUX(c_vars_dual, dcdx_vars_dual)
    rxn_temp  = RXN(c_vars_dual)

    BC_EAST_(1) = z_Li*nu_Li*flux_temp(1) - i_applied_cm2 / Fconst
    BC_EAST_(2) = flux_temp(2) - 0.0


  end function Boundary_EAST

! ******************************************************************************

end module GOV_EQNS


! ******************************************************************************
! ******************************** MAIN PROGRAM ********************************
! ******************************************************************************

program unsteady
  use user_input
  use variables, only : cprev, delC
  use write_data_mod
  use Electrochemical_Program_Mod

  implicit none
  integer :: t1, t2, clock_rate, clock_max
  integer :: it, j
  real :: tol = 1e-6, convergence_error = 1e3

  it = 0

  call system_clock(t1,clock_rate,clock_max)
  call initial_condition()
  call current_profile()
  call write_condition(it) ! writes the headers

  do it = 1, Numbertimesteps

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
        call write_condition(it)
        call current_profile()

        if (time >= experiment_time) then
          call write_condition(it)
          exit
        end if


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

  implicit none
  real    :: h
  integer :: j

  character(len=:), allocatable :: control_volume_input

  h = xmax/float(nj-2)
  do j=1,NJ

    if (j.EQ.1) then
      xx(j) = 0.0
    else if (j.EQ.NJ) then
      xx(NJ) = xmax
    else
      xx(j) = xx(1) + h*float(j-1) - h/2.0
    end if

  end do

  do j=2, NJ-1
     delx(j) = h
  end do
     delx(1) = 0.0                      ! it is common in the finite volume
     delx(NJ) = 0.0                     ! algorithm to set the control
                                        ! volumes at the boundaries to zero
  do j = 1, NJ
    cprev(1,j) = cbulk_Li
    cprev(2,j) = 0.0
  end do

  control_volume_input = trim(geometry)                   ! Define the control volume
  Cntrl_Vol = Control_Volume(control_volume_input)        ! size based on the
  Crx_Area = Cross_Sectional_Area(control_volume_input)   ! specified system geometry

  return                                                  ! is this necessary?
end subroutine initial_condition

include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/sub_auto_fill.f95'
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/sub_ABDGXY.f95'
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/sub_MATINV.f95'
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/sub_BAND.f95'
