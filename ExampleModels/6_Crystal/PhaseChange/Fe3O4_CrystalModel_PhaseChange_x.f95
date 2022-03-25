! * Written by Nicholas Brady August 10, 2020
! Uses John Newman's Band algorithm to solve coupled non-linear partial
! differential equations
! Incorporates DNAD which anables to use of automatic Differentiation to
! linearize the PDE's
! variables or names that begin and end with '_', i.e. _variable_ are changed by the python program: RunFortran.py

module user_input
  implicit none

  integer, parameter :: N = 2
  integer, parameter :: NJ = 42                          ! Number of mesh points
  integer, parameter :: Numbertimesteps = 3.6d3*1000     ! Number of time steps
  real               :: delT = 1.0                       ! size of timestep [s]
  real               :: time                             ! [s]

  ! **************************** Physical Constants ****************************
  real, parameter :: Rigc   = 8.314             ! Ideal gas constant [J/(mol*K)]
  real, parameter :: Temp   = 298               ! Temperature [K]
  real, parameter :: Fconst = 96485             ! Faraday's Constant [C/mol]
  real, parameter :: PI = 4.0 * ATAN(1.0)       ! pi - Geometric constant

  real, parameter :: density_Fe3O4  = 5.175            ! [g/cm3]
  real, parameter :: max_mAhg       = 926.0
  real, parameter :: mAhg_to_conc   = 3.6 * density_Fe3O4 / Fconst  ! mAh/g (x C/mAh * g/cm3 * mol/C --> mol/cm3)
  real, parameter :: max_Li_conc = max_mAhg * mAhg_to_conc

  ! **************************** Physical Parameters ***************************

  real :: diff_Li = 1e-18                       ! diffusion coefficient [cm^2/s]

  real :: cbulk = 0.001

  real :: xmax  = 10e-7         ! Crystal radius (cm) [10 nm]

  real :: c_x_init   = 1.0e-7

  real :: k_beta = 1e-4
  real :: c_beta = max_Li_conc
  real :: c_alpha_sat = max_Li_conc / 8.0

  real :: applied_specific_current             ! [mA/g]  926 mAh/g / xxx hours
  real :: applied_current_A
  real :: mAhg = 0.0                           ! Cummulative mAh/g

  integer :: cycle_number = 1                   ! cycle number counter


  character(len=3) :: state = 'D'               ! Discharge - D, Charge - C, Recovery - R

  character(len=65) :: geometry = 'Spherical'   ! system Geometry: Rectangular, Cylindrical, Spherical

end module user_input

! ******************************************************************************
! dnadmod and variables
! ******************************************************************************
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/dnadmod.f95'
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/variables.f95'
! ******************************************************************************

! ******************************************************************************
module echem_mod
  use user_input, only: N, Rigc, Fconst, Temp, max_Li_conc, cbulk
  use dnadmod

  implicit none

  real :: alpha_a = 0.5
  real :: alpha_c = 0.5
  real :: k_rxn   = 1.0e-8

  interface OCP
      module procedure OCP_real ! obtain OCP using real numbers
      module procedure OCP_dual ! obtain OCP using dual numbers
  end interface

contains


  function OCP_dual(c00, cii) result(OpenCircuitVoltage)
    type(dual)                :: OpenCircuitVoltage
    type(dual), intent(in)    :: c00, cii
    type(dual)                :: theta
    type(dual)                :: Vint
    real                      :: U_ref
    real, dimension(0:7)      :: AK
    integer                   :: kk


    theta = cii / max_Li_conc

    U_ref =  1.5617                                     ! U_ref  =  1.561677

    AK = (/ -6.5811e-1, &                               ! RK_A0  = -6.581135e-1
             6.5863e-3, &                               ! RK_A1  =  6.586353e-3
             1.2249e-1, &                               ! RK_A2  =  1.224915e-1
             2.7651e-1, &                               ! RK_A3  =  2.765079e-1
            -5.1470e-1, &                               ! RK_A4  = -5.146964e-1
            -1.2049e-4, &                               ! RK_A5  = -1.204923e-4
            -4.3649e-8, &                               ! RK_A6  = -4.364893e-8
             1.1099e-1 /)                               ! RK_A7  =  1.109934e-1

      Vint = 0.0
      if (theta == 0.5) then ! problems at theta = 0.5 b/c divide by 0
        ! ((2.0*theta-1.0)**(kk+1) = 0 for all k
        ! (2.0*theta*kk*(1.0-theta)) * (2.0*theta-1.0)**(kk-1)) = 0   for k != 1
        !                                                       = -0.5 for k = 1
        ! derivative also is "undefined" at theta = 0.5
        ! taking the limit as theta --> 0.5, we can simplify to:
        Vint%x  = -0.5 * AK(1)
        Vint%dx =  2.0 * AK(0) - 2.0 * AK(2)
      else
        do kk = 0, SIZE(AK)-1
            Vint = Vint + AK(kk)*( (2.0*theta-1.0)**(kk+1) - (2.0*theta*kk*(1.0-theta)) * (2.0*theta-1.0)**(kk-1) )
        end do
      end if

      OpenCircuitVoltage = U_ref + Rigc*Temp/Fconst * LOG( c00/cbulk * (1.0 - theta)/theta ) + Vint

  end function OCP_dual

  function OCP_real(c00, cii) result(OpenCircuitVoltage)
    real                      :: OpenCircuitVoltage
    real, intent(in)          :: c00, cii
    real                      :: theta
    real                      :: Vint
    real                      :: U_ref
    real, dimension(0:7)      :: AK
    integer                   :: kk


    theta = cii / max_Li_conc

    U_ref =  1.5617                                     ! U_ref  =  1.561677

    AK = (/ -6.5811e-1, &                               ! RK_A0  = -6.581135e-1
             6.5863e-3, &                               ! RK_A1  =  6.586353e-3
             1.2249e-1, &                               ! RK_A2  =  1.224915e-1
             2.7651e-1, &                               ! RK_A3  =  2.765079e-1
            -5.1470e-1, &                               ! RK_A4  = -5.146964e-1
            -1.2049e-4, &                               ! RK_A5  = -1.204923e-4
            -4.3649e-8, &                               ! RK_A6  = -4.364893e-8
             1.1099e-1 /)                               ! RK_A7  =  1.109934e-1

      Vint = 0.0
      if (theta == 0.5) then ! problems at theta = 0.5 b/c divide by 0
        ! ((2.0*theta-1.0)**(kk+1) = 0 for all k
        ! (2.0*theta*kk*(1.0-theta)) * (2.0*theta-1.0)**(kk-1)) = 0   for k != 1
        !                                                       = -0.5 for k = 1
        ! derivative also is "undefined" at theta = 0.5
        ! taking the limit as theta --> 0.5, we can simplify to:
        Vint  = -0.5 * AK(1)
      else
        do kk = 0, SIZE(AK)-1
            Vint = Vint + AK(kk)*( (2.0*theta-1.0)**(kk+1) - (2.0*theta*kk*(1.0-theta)) * (2.0*theta-1.0)**(kk-1) )
        end do
      end if

      OpenCircuitVoltage = U_ref + Rigc*Temp/Fconst * LOG( c00/cbulk * (1.0 - theta)/theta ) + Vint

  end function OCP_real



  function exchange_current_density(c0, c_x)  result(ex_curr)
  type(dual)              :: ex_curr
  type(dual), intent(in)  :: c0, c_x

  ex_curr = Fconst * k_rxn * c0**alpha_a * c_x**alpha_c * (max_Li_conc - c_x)**alpha_a

  end function exchange_current_density


  function Butler_Volmer_Rxn(c0, c_x, Phi_1, Phi_2)    result(i_BV)
    ! positive when eta > 0; negative when eta < 0
    ! eta positive when Phi_1 > U_ocp; eta negative when Phi_1 < U_ocp
    ! Phi_1 < U_ocp is lithiation (eta < 0, i_rxn < 0)
    ! Phi_1 > U_ocp is delithiation (eta > 0, i_rxn > 0)
    type(dual)              :: i_BV
    type(dual), intent(in)  :: c0, c_x, Phi_1, Phi_2
    type(dual)              :: eta, U_ocp
    type(dual)              :: i_ex_curr

    i_ex_curr = exchange_current_density(c0, c_x)

    U_ocp = OCP(c0, c_x)
    eta = Phi_1 - Phi_2 - U_ocp

    i_BV = i_ex_curr * ( EXP(alpha_a * Fconst * eta / (Rigc * Temp)) &
                        - EXP(-alpha_c * Fconst * eta / (Rigc * Temp)) )

  end function Butler_Volmer_Rxn


end module echem_mod


module write_data_mod
  use user_input
  use variables, only: cprev, delC, xx
  use echem_mod

  implicit none

  real :: c_x_Li, Phi_1, theta_beta

  character(len=65) :: header_fmt = '(1(A5,1X), 1(A5,1X), 3(A12,1X), 20(A15,1X)) '
  character(len=65) :: data_fmt   = '(1(I5,1X), 1(A5,1X), 3(F12.5,1X), 20(ES15.5,1X))'


  real :: t_write
  real :: Equiv
  real :: Li_balance


contains

  subroutine t_write__Equiv__Li_Bal_Assignment
    t_write     = time / 3600.0
    Equiv       = mAhg / max_mAhg * 8.

  end subroutine t_write__Equiv__Li_Bal_Assignment


  subroutine write_condition(it)
    integer :: it
    real :: last_write_time = 0.0
    real :: write_every_x_sec = 3600.0/10.0           ! 3600 s = 1 hour

    call t_write__Equiv__Li_Bal_Assignment()

    if (it == 1) then
      call write_to_screen(it)
      call write_positional_information_to_file(it)
      last_write_time = time
    else if ( (time - last_write_time).GE.write_every_x_sec ) then
      call write_to_screen(it)
      call write_positional_information_to_file(it)
      last_write_time = time
    end if
  end subroutine write_condition


  subroutine write_to_screen(it)
    integer :: it

    if (it.EQ.1) then       ! write the headers on the first entrance into write all voltage
      write(*, header_fmt) 'Cycle', 'State', 'Time', 'Equiv',    'Potential', 'conc_x', 'Theta_beta'
      write(*, header_fmt) '#', 'CDR',  'hours', 'LixFe3O4', 'Volts'  ,      'LixFe3O4', 'fraction'
    end if


    c_x_Li = cprev(1,NJ) / max_Li_conc * 8.0        ! Li_xFe_3O_4
    Phi_1  = OCP(cbulk, cprev(1,NJ))
    theta_beta = cprev(2,NJ)

    write(*, data_fmt) cycle_number, state, t_write, Equiv, Phi_1, c_x_Li, theta_beta

  end subroutine write_to_screen



  subroutine write_positional_information_to_file(it)
    use variables, only: xx
    integer :: it, j

    open(56, file = 'Time_Voltage_Position.txt', status = 'unknown')

    if (it.EQ.1) then
!           write the headers on the first entrance into write all voltage
      write(56, header_fmt) 'Cycle', 'State', 'Time', 'Equivalence', 'Voltage', 'Position', 'Solid_Conc', 'Theta_beta'
      write(56, header_fmt) '#',     'CDR',  'hours', 'LixFe3O4' , 'Volts'    , 'nm'      , 'LixFe3O4', 'fraction'
                                                    !           !
    end if                                          !           !

    Phi_1  = OCP(cbulk, cprev(1,NJ))
    do j = 1, NJ                                    !           !
      c_x_Li = cprev(1,j) / max_Li_conc * 8.0        ! Li_xFe_3O_4
                                                    !           !
      write(56, data_fmt) cycle_number, state,    t_write,    Equiv,       Phi_1,     xx(j)*1e7,    c_x_Li,   cprev(2,j)

    end do

  end subroutine write_positional_information_to_file


end module write_data_mod



module GOV_EQNS
  use user_input
  use dnadmod
  use echem_mod
  implicit none

contains

! ******************************************************************************
! **************************** GOVERNING EQUATIONS *****************************
! **************** Accumulation = FluxIn - FluxOut + Generation ****************
! ******************************************************************************
! can define intermediate variables with the functions to improve readability
! i.e c0 = c_vars_dual(1), dPhi2_dx = dcdx_vars_dual(2)

  ! (1) d(c_a * 0_a + c_b * 0_b)/dt = D d^2c_a/dx^2
  ! (2) c_b d(0_b)/dt = k_b (c_a - c_{a,sat}) * (1 - 0_b)

  function FLUX(c_vars_dual, dcdx_vars_dual)                      result(Flux_)
    ! (1)   N = -D * dc(1)/dx
    ! (2)   N = 0.0
    type(dual), dimension(N)              :: Flux_
    type(dual), dimension(N), intent(in)  :: c_vars_dual, dcdx_vars_dual

    type(dual) :: c_a, theta_b
    type(dual) :: dc_adx

    c_a     = c_vars_dual(1)
    dc_adx  = dcdx_vars_dual(1)

    Flux_(1) = -diff_Li * dc_adx
    Flux_(2) = 0.0

  end function FLUX

  function RXN(c_vars_dual)                                       result(Rxn_)
    ! (1)  0
    ! (2) k_b (c_a - c_{a,sat}) (1 - 0_b)
    type(dual), dimension(N)               :: Rxn_
    type(dual), dimension(N), intent(in)   :: c_vars_dual
    type(dual) :: c_a, theta_b

    c_a       = c_vars_dual(1)
    theta_b   = c_vars_dual(2)

    Rxn_(1) = 0

    if (c_a%x >= c_alpha_sat) then
      Rxn_(2) = k_beta * (c_a - c_alpha_sat) * (1.0 - theta_b)
    else if ( (c_a%x < c_alpha_sat).AND.(theta_b%x > 0) ) then
      Rxn_(2) = k_beta * (c_a - c_alpha_sat) * theta_b
    else
      Rxn_(2) = 0
    end if

  end function RXN

  function ACCUM(c_vars_dual)                                    result(Accum_)
    ! (1) d(c_a * 0_a + c_b * 0_b)/dt
    ! (2) d(c_b * 0_b)/dt
    type(dual), dimension(N)              :: Accum_
    type(dual), dimension(N), intent(in)  :: c_vars_dual
    type(dual) :: c_a, theta_b, theta_a

    c_a         = c_vars_dual(1)
    theta_b     = c_vars_dual(2)
    theta_a     = 1.0 - theta_b


    Accum_(1) = (c_a * theta_a + c_beta * theta_b)/delT
    Accum_(2) = (c_beta * theta_b)/delT                 ! d(c_b * 0_b)/dt

    ! accumulation terms are implicitly d/dt, so the function value is unphysical; only the derivatives are physically relevant
    Accum_(:)%x = 0.0
    ! setting Accum_(:)%x = 0 makes it easier to use ACCUM in writing boundary condtions

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
  function Boundary_WEST (c_vars_dual, dcdx_vars_dual)         result(BC_WEST_)
    ! (1)   N_x = 0
    ! (2) d(c_b * 0_b)/dt = k_b (c_a - c_{a,sat}) (1 - 0_b)
    type(dual), dimension(N)               :: BC_WEST_
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual
    type(dual), dimension(N)               :: flux_temp, rxn_temp, accum_temp

    type(dual) :: c, theta_b

    c       = c_vars_dual(1)
    theta_b = c_vars_dual(2)

    flux_temp   = FLUX(c_vars_dual, dcdx_vars_dual)
    rxn_temp    = RXN(c_vars_dual)
    accum_temp  = ACCUM(c_vars_dual)

    BC_WEST_(1) = flux_temp(1) - 0.0
    BC_WEST_(2) = accum_temp(2) - rxn_temp(2)

  end function Boundary_WEST

  function Boundary_EAST (c_vars_dual, dcdx_vars_dual)         result(BC_EAST_)
    ! (1)   N_x * (0_a) = i_app / F
    ! (2) d(c_b * 0_b)/dt = k_b (c_a - c_{a,sat}) (1 - 0_b)
    type(dual), dimension(N)               :: BC_EAST_
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual
    type(dual), dimension(N)               :: flux_temp, rxn_temp, accum_temp

    type(dual) :: c, theta_b, theta_a

    c       = c_vars_dual(1)
    theta_b = c_vars_dual(2)
    theta_a = 1. - theta_b

    flux_temp   = FLUX(c_vars_dual, dcdx_vars_dual)
    rxn_temp    = RXN(c_vars_dual)
    accum_temp  = ACCUM(c_vars_dual)

    BC_EAST_(1) = flux_temp(1) * theta_a - applied_current_A / Fconst
    BC_EAST_(2) = accum_temp(2) - rxn_temp(2)

  end function Boundary_EAST

end module GOV_EQNS





module experiment_mod
  use user_input
  use variables, only: cprev, delX
  use echem_mod
  use write_data_mod
  implicit none

  real    :: c_rate
  real    :: discharge_time
  real    :: discharge_time_0
  real    :: charge_time, rest_time, CCV_time
  character(len=:), allocatable :: experimental_program

contains

  subroutine experimental_conditions(it)
  integer :: it

  experimental_program = 'C_200_Discharge_Rest_Charge'
  ! _expt_program_ changed by python
  ! 'C_200_Discharge_Rest_Charge'
  ! 'C_200_Discharge_Charge_Rest'
  ! 'C_100_Discharge_Charge_cycle'
  ! 'C_10_Discharge'

! *****************************************************************************
  if (experimental_program == 'C_200_Discharge_Rest_Charge') then

    c_rate = 1.0/200.0

    discharge_time = 1.0/c_rate * ( 8.0 / 8.0) * 3600.0 !
    rest_time      = 20 * 3600                  ! 200 hours of recovery time

    if (time <= discharge_time) then
      state = 'D'
      applied_specific_current = +max_mAhg * c_rate

    else if (time <= (discharge_time + rest_time)) then
      state = 'R'
      applied_specific_current = 0.0

    else if (time > (discharge_time + rest_time)) then
      state = 'C'
      applied_specific_current = -max_mAhg * c_rate

    end if

  end if

! *****************************************************************************
  mAhg = mAhg + applied_specific_current * delT/3600.0
  if (trim(geometry) == 'Rectangular') then
    applied_current_A = -applied_specific_current/1000 * density_Fe3O4 * xmax

  else if (trim(geometry) == 'Cylindrical') then
    applied_current_A = -applied_specific_current/1000 * density_Fe3O4 * xmax/2.0

  else if (trim(geometry) == 'Spherical') then
    applied_current_A = -applied_specific_current/1000 * density_Fe3O4 * xmax/3.0
  end if



  if (time > 1000 * 3600) then
    print*, 'Time > 1000 hours'
    stop
  else if (cprev(1,NJ) < 0.0) then
    call t_write__Equiv__Li_Bal_Assignment
    call write_to_screen(0)
    print*, 'c0 < 0.0'
    stop
  else if (any(isnan(cprev))) then
    call t_write__Equiv__Li_Bal_Assignment
    call write_to_screen(0)
    print*, 'Dependent Variable is NaN'
    stop
  end if

  end subroutine experimental_conditions


end module experiment_mod

! ******************************************************************************
! auto_fill, ABDGXY, BAND, MATINV
! ******************************************************************************
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/BAND_MOD.f95'
! ******************************************************************************

! ******************************************************************************
! ******************************** MAIN PROGRAM ********************************
! ******************************************************************************
program unsteady
  use user_input
  use variables, only : cprev, delC
  use write_data_mod
  use echem_mod
  use experiment_mod
  use BAND_mod

  implicit none
  integer :: t1, t2, clock_rate, clock_max, it
  integer :: i, k, j
  real :: crystal_volume, crystal_surface_area

  call system_clock(t1,clock_rate,clock_max)
  call initial_condition()

  do it = 1, Numbertimesteps

    call experimental_conditions(it)
    call write_condition(it)

    do j = 1, NJ
      call auto_fill(j)               ! These can be changed to functions
      call ABDGXY(j)                  ! because they each have a set of inputs
      call BAND(j)                    ! and a set of outputs
    end do

    cprev = cprev + delC              ! Update the dependent variables
    time  = time  + delT              ! Update the time

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
      xx(j) = h*float(j-1) - h/2.0
    end if

  end do

  do j=2, NJ-1
     delx(j) = h
  end do
     delx(1) = 0.0                      ! it is common in the finite volume
     delx(NJ) = 0.0                     ! algorithm to set the control
                                        ! volumes at the boundaries to zero
  do j = 1, NJ
    cprev(1,j) = c_x_init
    cprev(2,j) = 0.0
  end do

  control_volume_input = trim(geometry)                   ! Define the control volume
  Cntrl_Vol = Control_Volume(control_volume_input)        ! size based on the
  Crx_Area = Cross_Sectional_Area(control_volume_input)   ! specified system geometry

  return
end subroutine initial_condition
