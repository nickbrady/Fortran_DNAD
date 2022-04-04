! * Written by Nicholas Brady August 10, 2020
! Uses John Newman's Band algorithm to solve coupled non-linear partial
! differential equations
! Incorporates DNAD which anables to use of automatic Differentiation to
! linearize the PDE's
! variables or names that begin and end with '_', i.e. _variable_ are changed by the python program: RunFortran.py

module user_input
  implicit none

  integer, parameter :: N = 4
  integer, parameter :: NJ = 102                          ! Number of mesh points

  integer, parameter :: Numbertimesteps = 3.6d3*1000     ! Number of time steps
  real               :: delT = 1.0                       ! size of timestep [s]
  real               :: time                             ! [s]
  logical            :: UPWIND = .FALSE.
  character(len=65)  :: direction = ''

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

  real :: diff_Li = 1e-13                       ! diffusion coefficient [cm^2/s]
  real :: diff_PF6 = 1e-13                      ! diff coeff PF6
  real :: z_1 = 1.0                             ! cation charge z_Li
  real :: z_2 = -1.0                            ! anion charge, z_PF6

  real :: cbulk = 0.001
  real :: agg_porosity = 0.26
  real :: volumetric_surface_area

  real :: xmax          = 1e-4              ! 500 um is 500e-4 cm
  real :: crystal_xmax  = 10e-7         ! Crystal radius (cm) [10 nm]

  ! real :: c_0_init  = 1.0e-7              ! initial c0 lithium conc  [mol/cm3]
  real :: c_x_init   = 1.0e-7
  real :: Phi_1_init = 3.3
  real :: Phi_2_init = 0.0
  real :: sigma      = 1.0e1

  real :: CCV_Voltage = 3.0
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

contains


  function OCP(c00, cii) result(OpenCircuitVoltage)
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

  end function OCP


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




! ******************************************************************************
module write_data_mod
  use user_input
  use variables, only: cprev, delC, xx
  use echem_mod

  implicit none

  real :: c_0_Li, c_x_Li, Phi_1, Phi_2

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
      write(*, header_fmt) 'Cycle', 'State', 'Time', 'Equiv',    'Voltage', 'Solution_Pot', 'c_Li_0', 'conc_x'
      write(*, header_fmt) '#', 'CDR',  'hours', 'LixFe3O4', 'Volts'  ,      'Volts',        'mol/L',  'LixFe3O4'
    end if


    c_0_Li = cprev(1,1)*1e3
    c_x_Li = cprev(2,NJ) / max_Li_conc * 8.0        ! Li_xFe_3O_4
    Phi_1  = cprev(3,NJ)
    Phi_2  = cprev(4,1)

    write(*, data_fmt) cycle_number, state, t_write, Equiv, Phi_1, Phi_2, c_0_Li, c_x_Li

  end subroutine write_to_screen



  subroutine write_edge_information_to_file(it)
    integer :: it




  end subroutine



  subroutine write_positional_information_to_file(it)
    use variables, only: xx
    integer :: it, j

    open(56, file = 'Time_Voltage_Position.txt', status = 'unknown')

    if (it.EQ.1) then
!           write the headers on the first entrance into write all voltage
      write(56, header_fmt) 'Cycle', 'State', 'Time', 'Equivalence', 'Position', 'Voltage', 'Soln_Conc', 'Solid_Conc', &
      & 'Solution_Pot'
      write(56, header_fmt) '#',     'CDR',  'hours', 'LixFe3O4' , 'um'      , 'Volts'    , 'mol/cm3', 'LixFe3O4', &
      & 'Volts'
                                                    !           !
    end if                                          !           !

    do j = 1, NJ                                    !           !
      c_0_Li = cprev(1,j) * 1e3
      c_x_Li = cprev(2,j) / max_Li_conc * 8.0        ! Li_xFe_3O_4
      Phi_1  = cprev(3,j)
      Phi_2  = cprev(4,j)
                                                    !           !
      write(56, data_fmt) cycle_number, state,    t_write,    Equiv,       xx(j)*1e4,    Phi_1,      c_0_Li,  c_x_Li, Phi_2

    end do

  end subroutine write_positional_information_to_file


end module write_data_mod



module GOV_EQNS
  use user_input
  ! use variables, only: xx, delX
  use dnadmod
  use echem_mod
  implicit none

  type(dual) :: c0, c_x, Phi_1, Phi_2
  type(dual) :: dc0dx, dc_xdx, dPhi_1dx, dPhi_2dx
  type(dual) :: i_rxn
contains

! ******************************************************************************
! **************************** GOVERNING EQUATIONS *****************************
! **************** Accumulation = FluxIn - FluxOut + Generation ****************
! ******************************************************************************
! can define intermediate variables with the functions to improve readability
! i.e c0 = c_vars_dual(1), dPhi2_dx = dcdx_vars_dual(2)

  function FLUX(c_vars_dual, dcdx_vars_dual) result(Flux_)
    ! (1)   N_0 = -eps * (D_Li * dc0/dx + z_Li * u_Li * c0 * F * dPhi_2/dx)
    ! (2)   N_x = 0
    ! (3)   i_1 = -(1-eps) * sigma * dPhi_1/dx
    ! (4)   i_2/F = -eps * [ (z_1 * diff_1 + z_2 * diff_2)*dc0/dx
    !                       (z_1**2 * u_1 + z_2**2 * u_2) * F * c0 * dPhi_2/dx]
    !       i.e. i_2 / F = sum_i (z_i * N_i)
    type(dual), dimension(N)              :: Flux_
    type(dual), dimension(N), intent(in)  :: c_vars_dual, dcdx_vars_dual
    type(dual)                            :: diff_
    real                                  :: u_1, u_2    ! mobility

    c0    = c_vars_dual(1)
    c_x   = c_vars_dual(2)
    Phi_1 = c_vars_dual(3)
    Phi_2 = c_vars_dual(4)

    dc0dx    = dcdx_vars_dual(1)
    dc_xdx   = dcdx_vars_dual(2)
    dPhi_1dx = dcdx_vars_dual(3)
    dPhi_2dx = dcdx_vars_dual(4)

    ! Nernst-Einstein Relationship
    u_1 = diff_Li / (Rigc * Temp)
    u_2 = diff_PF6 / (Rigc * Temp)

    Flux_(1) = -agg_porosity * (diff_Li * dc0dx + z_1 * u_1 * c0 * Fconst * dPhi_2dx)
    Flux_(2) = 0.0
    Flux_(3) = -(1 - agg_porosity) * sigma * dPhi_1dx
    Flux_(4) = -agg_porosity *( Fconst*(z_1 * diff_Li + z_2 * diff_PF6) * dc0dx &
        + (z_1**2 * u_1 + z_2**2 * u_2) * Fconst**2 * c0 * dPhi_2dx )

  end function FLUX

  function RXN(c_vars_dual) result(Rxn_)
    ! (1)  a * i_rxn / F
    ! (2) -a * i_rxn / F
    ! (3) -a * i_rxn
    ! (4)  a * i_rxn
    type(dual), dimension(N)               :: Rxn_
    type(dual), dimension(N), intent(in)   :: c_vars_dual

    c0    = c_vars_dual(1)
    c_x   = c_vars_dual(2)
    Phi_1 = c_vars_dual(3)
    Phi_2 = c_vars_dual(4)

    i_rxn = Butler_Volmer_Rxn(c0, c_x, Phi_1, Phi_2)

    Rxn_(1) = +volumetric_surface_area * i_rxn / Fconst     ! c_0
    Rxn_(2) = -volumetric_surface_area * i_rxn / Fconst     ! c_x
    Rxn_(3) = -volumetric_surface_area * i_rxn              ! i_1
    Rxn_(4) = +volumetric_surface_area * i_rxn              ! i_2

  end function RXN

  function ACCUM(c_vars_dual) result(Accum_)
    ! (1) eps * dc0/dt
    ! (2) (1-eps) * dc_x/dt
    ! (3) 0                       (no accumulation of e-)
    ! (4) 0                       (electroneutrality - ions)
    type(dual), dimension(N)              :: Accum_
    type(dual), dimension(N), intent(in)  :: c_vars_dual

    c0    = c_vars_dual(1)
    c_x   = c_vars_dual(2)
    Phi_1 = c_vars_dual(3)
    Phi_2 = c_vars_dual(4)

    Accum_(1) = (agg_porosity) * c0/delT
    Accum_(2) = (1.0 - agg_porosity) * c_x/delT
    Accum_(3) = 0.0                  ! no accumulation of e-
    Accum_(4) = 0.0                  ! electroneutrality

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
    ! (1)   N_0 = 0
    ! (2)   (1-eps) dc/dt = -a * i_rxn/F
    ! (3)   i_1 = 0                             (flux_3 = 0)
    ! (4)   i_2 = 0                             (flux_4 = 0)
    type(dual), dimension(N)               :: BC_WEST_
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual
    type(dual), dimension(N)               :: flux_temp
    type(dual)                             :: dc_x!, i_rxn

    c0    = c_vars_dual(1)
    c_x   = c_vars_dual(2)
    Phi_1 = c_vars_dual(3)
    Phi_2 = c_vars_dual(4)

    i_rxn = Butler_Volmer_Rxn(c0, c_x, Phi_1, Phi_2)

    flux_temp = FLUX(c_vars_dual, dcdx_vars_dual)

    dc_x   = c_vars_dual(2)
    dc_x%x = 0.0             ! --> d(c_x)  -- dc_x/delT = d(c_x)/dt

    BC_WEST_(1) = flux_temp(1) - 0.0
    BC_WEST_(2) = (1 - agg_porosity) * dc_x/delT + &  ! (1-eps)dc/dt = -ai_rxn/F
                & volumetric_surface_area * i_rxn/Fconst
    BC_WEST_(3) = flux_temp(3) - 0.0                              ! i_1 = 0.0
    BC_WEST_(4) = flux_temp(4) - 0.0                              ! i_2 = 0.0

  end function Boundary_WEST

  function Boundary_EAST (c_vars_dual, dcdx_vars_dual) result(BC_EAST_)
    ! (1)   c_0 = cbulk
    ! (2)   (1-eps) dc/dt = -a * i_rxn/F
    ! (3)   i_1   = i_applied     / OR /     Phi_1 = const_V
    ! (4)   Phi_2 = 0                                   (arbitrary ref Voltage)
    type(dual), dimension(N)               :: BC_EAST_
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual
    type(dual), dimension(N)               :: flux_temp
    type(dual)                             :: dc_x!, i_rxn

    c0    = c_vars_dual(1)
    c_x   = c_vars_dual(2)
    Phi_1 = c_vars_dual(3)
    Phi_2 = c_vars_dual(4)

    i_rxn = Butler_Volmer_Rxn(c0, c_x, Phi_1, Phi_2)

    flux_temp = FLUX(c_vars_dual, dcdx_vars_dual)

    dc_x   = c_vars_dual(2)
    dc_x%x = 0.0                    ! --> d(c_x)  -- dc_x/delT = d(c_x)/dt

    if (state == 'CCV') then
      ! constant voltage boundary condition
      BC_EAST_(3) = c_vars_dual(3) - CCV_Voltage
      ! applied_specific_current =
    else
      ! constant current boundary condition
      BC_EAST_(3) = flux_temp(3) - applied_current_A               ! i_1 = i_app
    end if

    BC_EAST_(1) = c0 - cbulk                                      ! c0 = cbulk
    BC_EAST_(2) = (1 - agg_porosity) * dc_x/delT + &  ! (1-eps)dc/dt = -ai_rxn/F
                & volumetric_surface_area * i_rxn/Fconst
    BC_EAST_(4) = c_vars_dual(4) - 0.0                            ! Phi_2 = 0.0

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

    discharge_time = 1.0/c_rate * ( 0.5 / 8.0) * 3600.0 !
    rest_time      = 0.0!20 * 3600                  ! 200 hours of recovery time

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
  applied_current_A = applied_specific_current/1000 * (1 - agg_porosity) * density_Fe3O4 * xmax/3.0

  ! if (time > 1000 * 3600) then
  !   print*, 'Time > 1000 hours'
  !   exit
  ! else if (cprev(1,NJ) < 0.0) then
  !   call t_write__Equiv__Li_Bal_Assignment
  !   call write_to_screen(0)
  !   print*, 'c0 < 0.0'
  !   exit
  ! else if (any(isnan(cprev))) then
  !   call t_write__Equiv__Li_Bal_Assignment
  !   call write_to_screen(0)
  !   print*, 'Dependent Variable is NaN'
  !   exit
  ! end if

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

  crystal_surface_area    = 4.0 * PI * crystal_xmax**2
  crystal_volume          = 4/3 * PI * crystal_xmax**3
  volumetric_surface_area = (1 - agg_porosity)/crystal_volume * crystal_surface_area

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

! *********** STOPPING CONDITIONS ************* !
    if (time > 1000 * 3600) then
      print*, 'Time > 1000 hours'
      exit
    else if (cprev(1,NJ) < 0.0) then
      call t_write__Equiv__Li_Bal_Assignment
      call write_to_screen(0)
      print*, 'c0 < 0.0'
      exit
    else if (any(isnan(cprev))) then
      call t_write__Equiv__Li_Bal_Assignment
      call write_to_screen(0)
      print*, 'Dependent Variable is NaN'
      exit
    end if

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
    cprev(1,j) = cbulk
    cprev(2,j) = c_x_init
    cprev(3,j) = Phi_1_init
    cprev(4,j) = Phi_2_init
  end do

  control_volume_input = trim(geometry)                   ! Define the control volume
  Cntrl_Vol = Control_Volume(control_volume_input)        ! size based on the
  Crx_Area = Cross_Sectional_Area(control_volume_input)   ! specified system geometry

  return
end subroutine initial_condition
