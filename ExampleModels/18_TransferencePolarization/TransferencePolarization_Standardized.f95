module user_input
  implicit none
  
  integer, parameter :: N = 1                        ! Number of mesh points
  integer, parameter :: NJ = 202                        ! Number of mesh points
  integer, parameter :: Numbertimesteps = 100 * 60 * 1e3       ! Number of time steps
  integer, parameter :: maxIterations = 1e6
  real               :: delT = 1.0e-3                   ! size of timestep [s]
  real               :: time                            ! [s]
  real, parameter    :: xmax = 2.0                      ! [cm] 1 cm
  logical            :: UPWIND = .FALSE.
  character(len=65)  :: direction = 'WestToEast'

  real, parameter :: Rigc   = 8.314                     ! Ideal gas constant [J/(mol*K)]
  real, parameter :: Fconst = 96485                     ! Faraday's Constant [C/mol]
  real, parameter :: Temp   = 298
  real, parameter :: PI = 4.0 * ATAN(1.0)               ! pi - Geometric constant

  real, parameter :: cbulk_Li = 1e-3
  real, parameter :: trans_Li = 0.5

  real, parameter :: D_bulk   = 1e-5
  real, parameter :: z_Li     = +1.0
  real, parameter :: nu_Li    = 1.0
  real, parameter :: s_Li     = +1.0

  real :: i_applied_cm2 = 10e-6           ! 1 μA/cm2

  character(len=65) :: geometry = 'Rectangular'   ! system Geometry: Rectangular, Cylindrical, Spherical

end module user_input

!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/variables.f95'
include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/dnadmod.f95'
!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*


module write_data_mod
  use user_input
  use variables, only: cprev, delC, xx

  implicit none

  character(len=65) :: header_fmt = '(1(A12,1X),   20(A15,1X)) '
  character(len=65) :: data_fmt   = '(1(F12.5,1X), 20(ES15.5,1X))'

  real :: t_write

  real :: conc_Li

contains

  subroutine t_write__Equiv__Li_Bal_Assignment
    t_write    = time

  end subroutine t_write__Equiv__Li_Bal_Assignment


  subroutine write_condition(it)
    integer :: it
    real :: last_write_time = -1e6
    real :: write_every_x_sec = 1.0e-1           ! 3600 s = 1 hour

    call t_write__Equiv__Li_Bal_Assignment

    if (time <= 360) then
      if ( (time - last_write_time) >= write_every_x_sec - delT*0.5) then
        call write_to_screen(it)
        call write_positional_information_to_file(it)

        last_write_time = time
      end if
      
    else
      if ( (time - last_write_time) >= 1.0 - delT*0.5) then
        call write_to_screen(it)
        call write_positional_information_to_file(it)

        last_write_time = time
      end if
    end if

  end subroutine write_condition


  subroutine write_to_screen(it)
    integer :: it

    if (it == 0) then       ! write the headers on the first entrance into write all voltage
      write(*, header_fmt) 'Time',  'conc_1',  'conc_NJ', 'vel_1', 'vel_NJ', 'ξ_max', 'ζ_max'
      write(*, header_fmt) 'seconds', 'mol/cm3', 'mol/cm3', 'cm/s' , 'cm/s'  , ''     , ''
    end if

      conc_Li   = cprev(1,1)

      write(*, data_fmt) t_write, conc_Li-cbulk_Li, cprev(1,2)-cbulk_Li, cprev(1,3)-cbulk_Li, cprev(1,4)-cbulk_Li, &
                        & cprev(1,5)-cbulk_Li


  end subroutine write_to_screen


  subroutine write_positional_information_to_file(it)
    use variables, only: xx
    integer :: it, j

    open(56, file = 'Time_Conc_Position.txt', status = 'unknown')

    if (it == 0) then
!           write the headers on the first entrance into write all voltage
      write(56, header_fmt) 'Time', 'Position', 'Delta_Conc',   'Velocity', 'Xi', 'Zeta'
      write(56, header_fmt) 'seconds', 'cm',      'mol/cm3', 'cm/s'   , ''  , ''
                                                      !
    end if                                              !


    do j = 1, NJ                                    !
      conc_Li    = cprev(1,j)
                                                    !
      write(56, data_fmt) t_write, xx(j), conc_Li - cbulk_Li

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

    Flux_(1) = -D_bulk * dcdx_Li


  end function FLUX

  function RXN(c_vars_dual)                   result(Rxn_)
    type(dual), dimension(N)               :: Rxn_
    type(dual), dimension(N), intent(in)   :: c_vars_dual

    Rxn_(1) = 0.0

  end function RXN

  function ACCUM(c_vars_dual)                           result(Accum_)
    type(dual), dimension(N)              :: Accum_
    type(dual), dimension(N), intent(in)  :: c_vars_dual

    Accum_(1) = c_vars_dual(1)/delT

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
    ! real :: position = 0.0

    c_Li        = c_vars_dual(1)
    dcdx_Li     = dcdx_vars_dual(1)

    flux_temp = FLUX(c_vars_dual, dcdx_vars_dual)
    rxn_temp  = RXN(c_vars_dual)

    BC_WEST_(1) = z_Li*nu_Li*flux_temp(1) - i_applied_cm2 * (1 - trans_Li)/Fconst

  end function Boundary_WEST

  function Boundary_EAST (c_vars_dual, dcdx_vars_dual) result(BC_EAST_)
    type(dual), dimension(N)               :: BC_EAST_, flux_temp
    type(dual), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual

    c_Li        = c_vars_dual(1)
    dcdx_Li     = dcdx_vars_dual(1)


    BC_EAST_(1) = c_Li - cbulk_Li


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
  real :: tol = 1e-6, convergence_error = 1e3

  it = 0

  call system_clock(t1,clock_rate,clock_max)
  call initial_condition()
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

        if (time >= 3*60) then
          i_applied_cm2 = 0.0
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
      xx(NJ) = 0.0 + xmax
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