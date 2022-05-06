module user_input
    implicit none
  
    integer, parameter :: N = 2
    integer, parameter :: NJ = 42                          ! Number of mesh points
    integer, parameter :: Numbertimesteps = 30 * 3.6e3     ! Number of time steps
    real               :: delT = 1.0                       ! size of timestep [s]
    real               :: time                             ! [s]
    real               :: xmax = 1e-0                      ! 1 cm
    logical            :: UPWIND = .FALSE.
    character(len=65)  :: direction = 'WestToEast'
  
    real, parameter :: Rigc   = 8.314             ! Ideal gas constant [J/(mol*K)]
    real, parameter :: Temp   = 298               ! Temperature [K]
    real, parameter :: Fconst = 96485             ! Faraday's Constant [C/mol]
    real, parameter :: PI = 4.0 * ATAN(1.0)       ! pi - Geometric constant
  
    real, parameter :: diff_cat  = 1e-5                       ! diffusion coefficient [cm^2/s]
    real, parameter :: diff_an   = diff_cat/2
    real :: u_cat  = diff_cat / (Rigc*Temp)
    real :: u_an   = diff_an / (Rigc*Temp)
    real, parameter :: nu_cat = 1.0
    real, parameter :: nu_an  = 1.0
    real, parameter :: z_cat  = +1.0
    real, parameter :: z_an   = -1.0
    real :: cbulk = 0.001                       ! concentration [mol/cm3]
  
    real :: i_applied_cm2 = -2e-3
    ! real :: Delta_V = 600e-3
  
    character(len=65) :: geometry = 'Rectangular'   ! system Geometry: Rectangular, Cylindrical, Spherical
  
  end module user_input
  !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*
  include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/variables.f95'
  include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/dnadmod.f95'
  !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*
  
  module echem_transport_mod
    use dnadmod
    use user_input
  
    implicit none
  
    real :: alpha_a = 0.5
    real :: alpha_c = 0.5
    real :: k_rxn   = 1.0e-7
  
    interface Soln_Pot
      module procedure Soln_Pot_dual
      module procedure Soln_Pot_real
    end interface Soln_Pot
  
    interface exchange_current_density
      module procedure exchange_current_density_dual
      module procedure exchange_current_density_real
    end interface exchange_current_density
  
    interface Butler_Volmer_Rxn
      module procedure Butler_Volmer_Rxn_dual
      module procedure Butler_Volmer_Rxn_real
      module procedure Butler_Volmer_Rxn_Phi_dual
    end interface Butler_Volmer_Rxn
  
    contains
  
      function Diffusion_Coef(c_vars_dual)    Result(diff_coef_)
        type(dual), dimension(N)              :: diff_coef_
        type(dual), dimension(N), intent(in)  :: c_vars_dual
  
        type(dual) :: c0
  
        real :: D_Li_0 = 1e-5
        real :: D_Li_min = 1e-7
        real :: c0_0 = 9.0e-4
  
        c0 = c_vars_dual(1)
  
        diff_coef_ = D_Li_0 * exp(-c0 / c0_0) + D_Li_min
  
      end function Diffusion_Coef
  
  
  !*******************************************************************************
      function Soln_Pot_dual(c0)                   result(SolutionPotential)
        type(dual)                :: SolutionPotential
        type(dual), intent(in)    :: c0
        real                      :: U_0
  
        U_0 =  0.0
  
        SolutionPotential = U_0 + Rigc*Temp/Fconst * LOG( c0/cbulk )
  
      end function Soln_Pot_dual
  
  
      function Soln_Pot_real(c0)                   result(SolutionPotential)
        real                :: SolutionPotential
        real, intent(in)    :: c0
        real                :: U_0
  
        U_0 =  0.0
  
        SolutionPotential = U_0 + Rigc*Temp/Fconst * LOG( c0/cbulk )
  
      end function Soln_Pot_real
  !*******************************************************************************
  
  !*******************************************************************************
      function exchange_current_density_dual(c0)  result(ex_curr)
        type(dual)              :: ex_curr
        type(dual), intent(in)  :: c0
  
        ex_curr = Fconst * k_rxn * c0**alpha_a
  
      end function exchange_current_density_dual
  
      function exchange_current_density_real(c0)  result(ex_curr)
        real              :: ex_curr
        real, intent(in)  :: c0
  
        ex_curr = Fconst * k_rxn * c0**alpha_a
  
      end function exchange_current_density_real
  !*******************************************************************************
  
  !*******************************************************************************
      function Butler_Volmer_Rxn_dual(c0, Phi_1)    result(i_BV)
        ! positive when eta > 0; negative when eta < 0
        ! eta positive when Phi_1 > U_ocp; eta negative when Phi_1 < U_ocp
        ! Phi_1 < U_ocp is lithiation (eta < 0, i_rxn < 0)
        ! Phi_1 > U_ocp is delithiation (eta > 0, i_rxn > 0)
        type(dual)              :: i_BV
        type(dual), intent(in)  :: c0, Phi_1
        type(dual)              :: eta, U_ocp
        type(dual)              :: i_ex_curr
  
        i_ex_curr = exchange_current_density(c0)
  
        ! U_ocp = Soln_Pot(c0)
        ! eta = Phi_1 - Phi_2 - U_ocp
        eta = Phi_1
  
        i_BV = i_ex_curr * ( EXP(alpha_a * Fconst * eta / (Rigc * Temp)) &
                            - EXP(-alpha_c * Fconst * eta / (Rigc * Temp)) )
  
      end function Butler_Volmer_Rxn_dual
  
      function Butler_Volmer_Rxn_Phi_dual(c0, Phi_1)    result(i_BV)
        ! positive when eta > 0; negative when eta < 0
        ! eta positive when Phi_1 > U_ocp; eta negative when Phi_1 < U_ocp
        ! Phi_1 < U_ocp is lithiation (eta < 0, i_rxn < 0)
        ! Phi_1 > U_ocp is delithiation (eta > 0, i_rxn > 0)
        type(dual)              :: i_BV
        real, intent(in)        :: c0
        type(dual), intent(in)  :: Phi_1
        type(dual)              :: eta
        real                    :: i_ex_curr, U_ocp
  
        i_ex_curr = exchange_current_density(c0)
  
        ! U_ocp = Soln_Pot(c0)
        ! eta = Phi_1 - Phi_2 - U_ocp
        eta = Phi_1
  
        i_BV = i_ex_curr * ( EXP(alpha_a * Fconst * eta / (Rigc * Temp)) &
                            - EXP(-alpha_c * Fconst * eta / (Rigc * Temp)) )
  
      end function Butler_Volmer_Rxn_Phi_dual
  
      function Butler_Volmer_Rxn_real(c0, Phi_1)    result(i_BV)
        ! positive when eta > 0; negative when eta < 0
        ! eta positive when Phi_1 > U_ocp; eta negative when Phi_1 < U_ocp
        ! Phi_1 < U_ocp is lithiation (eta < 0, i_rxn < 0)
        ! Phi_1 > U_ocp is delithiation (eta > 0, i_rxn > 0)
        real              :: i_BV
        real, intent(in)  :: c0, Phi_1
        real              :: eta, U_ocp
        real              :: i_ex_curr
  
        i_ex_curr = exchange_current_density(c0)
  
        ! U_ocp = Soln_Pot(c0)
        ! eta = Phi_1 - Phi_2 - U_ocp
        eta = Phi_1
  
        i_BV = i_ex_curr * ( EXP(alpha_a * Fconst * eta / (Rigc * Temp)) &
                            - EXP(-alpha_c * Fconst * eta / (Rigc * Temp)) )
  
      end function Butler_Volmer_Rxn_real
  !*******************************************************************************
  
  end module echem_transport_mod
  !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*
  
  module write_data_mod
    use user_input
    use variables, only: cprev, delC, xx
    use echem_transport_mod
  
    implicit none
  
    real :: c0, Phi_2
    real :: Phi_2_1, Phi_2_NJ
    real :: i_rxn, i_rxn_1, i_rxn_NJ
  
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
      integer :: it, j
      real :: Delta_c, Delta_x, dcdx, Delta_Phi, D_Phi
      real :: c_cat, t_cat, c_an
  
      if (it.EQ.1) then       ! write the headers on the first entrance into write all voltage
        write(*, header_fmt) 'Time',  'Conc', 'ΔΦ_2', 'Δc0', 'ΔΦ_2'
        write(*, header_fmt) 'hours', 'mol/L', 'Volts', 'mol/cm3', 'Volts'
      end if
  
  
      Phi_2_1 = cprev(2,1)
      Phi_2_NJ = cprev(2,NJ)
  
  
      Delta_Phi = 0.0
      do j = 2,NJ
        c0 = cprev(1,j)
        c_cat = c0*nu_cat
        c_an  = c0*nu_an
        t_cat = z_cat**2 * u_cat * c_cat/(z_cat**2 * u_cat * c_cat + z_an**2 * u_an * c_an)
        Delta_c = c0 - cprev(1,j-1)
  
        Delta_Phi = Delta_Phi - Rigc*Temp/Fconst * 1/(z_cat**2*diff_cat) * t_cat/c_cat*(z_cat*diff_cat + z_an*diff_an)*Delta_c
      end do
      c0 = cprev(1,NJ)*1e3
      Delta_c  = cprev(1,NJ) - cprev(1,1)
  
      D_Phi = -Rigc*Temp/Fconst * 1/(z_cat*diff_cat) * t_cat*(z_cat*diff_cat + z_an*diff_an)*LOG(cprev(1,NJ)/cprev(1,1))
  
      write(*, data_fmt) t_write,    c0, Phi_2_NJ - Phi_2_1, Delta_c, Delta_Phi, D_Phi
  
    end subroutine write_to_screen
  
  
    subroutine write_positional_information_to_file(it)
      use variables, only: xx
      integer :: it, j
  
      open(56, file = 'Time_Conc_Position.txt', status = 'unknown')
  
      if (it.EQ.1) then
  !           write the headers on the first entrance into write all voltage
        write(56, header_fmt) 'Time',  'Position', 'Conc', 'Phi_1', 'i_rxn'
        write(56, header_fmt) 'hours', 'cm'      , 'mol/L',  'Volts', 'A/cm2'
                                                      !
      end if                                          !
  
  
      ! Phi_2 = Delta_V - Phi_1
                                                      !
      do j = 1, NJ                                    !
        c0 = cprev(1,j) * 1e3                         !
        Phi_2 = cprev(2,j)
                                                      !
        write(56, data_fmt) t_write,    xx(j),    c0, Phi_2
  
      end do
  
    end subroutine write_positional_information_to_file
  
  
  end module write_data_mod
  !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*
  
  module GOV_EQNS
    use user_input
    use variables
    use echem_transport_mod
    use dnadmod
  
    implicit none
  
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
  
      type(dual) :: c0, c_cat, c_an
      type(dual) :: dc0dx, dc_cat_dx, dc_an_dx
      type(dual) :: dPhi2_dx
  
      c0        = c_vars_dual(1)
      dc0dx     = dcdx_vars_dual(1)
  
      dPhi2_dx  = dcdx_vars_dual(2)
  
      c_cat     = nu_cat*c0
      dc_cat_dx = nu_cat*dc0dx
  
      c_an      = nu_an*c0
      dc_an_dx  = nu_an*dc0dx
  
      Flux_(1) = -diff_cat * dc_cat_dx - z_cat*u_cat*Fconst*c_cat*dPhi2_dx
      Flux_(2) = -diff_an  * dc_an_dx  - z_an*u_an*Fconst*c_an*dPhi2_dx
  
    end function FLUX
  
    function RXN(c_vars_dual) result(Rxn_)
      ! (1)  no reaction
      type(dual), dimension(N)               :: Rxn_
      type(dual), dimension(N), intent(in)   :: c_vars_dual
  
      type(dual) :: c0, Phi_1
  
      c0    = c_vars_dual(1)
      Phi_1 = c_vars_dual(2)
  
      Rxn_(1) = 0.0
      Rxn_(2) = 0.0
  
    end function RXN
  
    function ACCUM(c_vars_dual) result(Accum_)
      ! (1) dc0/dt
      type(dual), dimension(N)              :: Accum_
      type(dual), dimension(N), intent(in)  :: c_vars_dual
  
      type(dual) :: c0, c_cat, c_an
  
      c0    = c_vars_dual(1)
      c_cat     = nu_cat*c0
      c_an      = nu_an*c0
  
      Accum_(1) = c_cat/delT
      Accum_(2) = c_an/delT
  
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
      ! (1)   c0 = cbulk
      type(dual), dimension(N)               :: BC_WEST_, flux_temp
      type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual
      type(dual) :: Phi_2
  
      flux_temp = FLUX(c_vars_dual, dcdx_vars_dual)
      Phi_2 = c_vars_dual(2)
  
  
      BC_WEST_(1) = z_cat*Fconst*flux_temp(1) - i_applied_cm2
      BC_WEST_(2) = Phi_2 - 0.0
  
    end function Boundary_WEST
  
    function Boundary_EAST (c_vars_dual, dcdx_vars_dual) result(BC_EAST_)
      ! (1)   c0 = 0.1 * cbulk
      type(dual), dimension(N)               :: BC_EAST_, flux_temp
      type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual
  
      flux_temp = FLUX(c_vars_dual, dcdx_vars_dual)
  
      BC_EAST_(1) = z_cat*Fconst*flux_temp(1) - i_applied_cm2
      BC_EAST_(2) = z_an*Fconst*flux_temp(2) - 0.0
  
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
  
    implicit none
    integer :: t1, t2, clock_rate, clock_max
    integer :: it, j
  
    call system_clock(t1,clock_rate,clock_max)
    call initial_condition()
  
    do it = 1, Numbertimesteps
  
      call write_condition(it)
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
  
      if (time >= 10*3600) then
        i_applied_cm2 = 0.0
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
    use GOV_EQNS, only: Control_Volume, Cross_Sectional_Area,  &    ! functions
                        Cntrl_Vol, Crx_Area                         ! variables
  
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