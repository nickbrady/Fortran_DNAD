! * Written by Nicholas Brady June 3, 2022
! Uses John Newman's Band algorithm to solve coupled non-linear partial
! differential equations
! Incorporates DNAD which anables to use of automatic Differentiation to
! linearize the PDE's
! variables or names that begin and end with '_', i.e. _variable_ are changed by the python program: RunFortran.py


module user_input
  implicit none

  integer, parameter :: N = 1
  integer, parameter :: NJ = 42                   ! Number of mesh points

  integer, parameter :: Numbertimesteps = 1e2     ! Number of time steps
  real               :: delT = 1.0                ! size of timestep [s]
  real               :: time                      ! [s]
  logical            :: UPWIND = .FALSE.
  character(len=65)  :: direction = ''

  real               :: xmax = 1.0                ! 500 um is 500e-4 cm

  real, parameter :: PI = 4.0 * ATAN(1.0)         ! pi - Geometric constant

  real :: diff  = 1e-2                            ! diffusion coefficient [cm^2/s]
  real :: cbulk = 1.0                             ! concentration [mol/cm3]
  real :: freq  = 1e2

  character(len=65) :: geometry = 'Rectangular'   ! system Geometry: Rectangular, Cylindrical, Spherical

end module user_input

! ******************************************************************************
! dnadmod and variables
! ******************************************************************************
! include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/Impedance_Complex/dnadmod.f95'
! include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/Impedance_Complex/variables.f95'
include 'dnadmod.f95'
include 'variables.f95'
! ******************************************************************************

! ******************************************************************************
module write_data_mod
  use user_input
  use variables, only: cprev, delC, xx

  implicit none

  complex :: c0

  character(len=65) :: header_fmt = '(20(A15,1X)) '
  character(len=65) :: data_fmt   = '(20(ES15.5E3,1X))'

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
    ! call write_to_screen(it)

    if (it == 1) then
      ! call write_to_screen(it)
      call write_positional_information_to_file(it)
      last_write_time = time

    else if ( (time - last_write_time).GE.write_every_x_sec ) then
      ! call write_to_screen(it)
      call write_positional_information_to_file(it)
      last_write_time = time

    end if
  end subroutine write_condition


  subroutine write_to_screen(it)
    integer :: it
    logical :: first_write = .TRUE.

    if (first_write) then       ! write the headers on the first entrance into write all voltage
      first_write = .FALSE.
      write(*, header_fmt) 'Diff_Coef', 'Frequency', 'Real(c)', 'Imag(c)'
      write(*, header_fmt) 'cm2/s'    , 'Hz'       , ''       , ''
    end if

    c0 = cprev(1,1)

    write(*, data_fmt) diff, freq, realpart(c0), imagpart(c0)

  end subroutine write_to_screen

  subroutine write_to_file(it)
    integer :: it
    logical :: first_write = .TRUE.

    open(60, file = 'Impedance_Data_.txt', status = 'unknown')

    if (first_write) then       ! write the headers on the first entrance into write all voltage
      first_write = .FALSE.
      write(60, header_fmt) 'Diff_Coef', 'Frequency', 'Real(c)', 'Imag(c)'
      write(60, header_fmt) 'cm2/s'    , 'Hz'       , ''       , ''
    end if

    c0 = cprev(1,1)

    write(60, data_fmt) diff, freq, realpart(c0), imagpart(c0)

  end subroutine write_to_file


  subroutine write_positional_information_to_file(it)
    use variables, only: xx
    integer :: it, j
    logical :: first_write = .TRUE.

    open(56, file = 'Time_Conc_Position_Transmissive.txt', status = 'unknown')

    if (first_write) then       ! write the headers on the first entrance into write all voltage
      first_write = .FALSE.
!           write the headers on the first entrance into write all voltage
      write(56, header_fmt) 'Frequency',  'Position', 'Real(c)', 'Imag(c)'
      write(56, header_fmt) 'Hz'       ,  'um'      , 'mol/L'  , 'mol/L'
                                                    !
    end if                                          !
                                                    !
    do j = 1, NJ                                    !
      c0 = cprev(1,j)                               !
                                                    !
      write(56, data_fmt) freq         ,    xx(j),  realpart(c0), imagpart(c0)

    end do

  end subroutine write_positional_information_to_file


end module write_data_mod



module GOV_EQNS
  use user_input
  use dnadmod

  implicit none

contains
! General Equations
! Accumulation = FLUX_in - FLUX_out + Reaction
! Accum = FLUX_west - FLUX_east + Rxn

! Analytic:
!     ∂cᵢ/∂t = -∇⋅𝐍ᵢ + Rᵢ
! Finite Volume (Control Volume):
!     ΔV ∂c/∂t = (Aₓᵢ⋅𝐍ᵢ - Aₓₒ⋅𝐍ₒ) + ΔV ⋅ Rⱼ
!
! frequently, the control volume (ΔV) and cross-sectional area (Aₓ) are given
! by the system geometry
! 𝐍ᵢ is the flux of specie i
! Rᵢ is the rate of generation (reaction rate) of specie i

! (1) ∂c/∂t = D∇²c
! 𝐍 = -D ∇c
! R = 0

! ******************************************************************************
! **************************** GOVERNING EQUATIONS *****************************
! ******************************************************************************
! **************** Accumulation = FluxIn - FluxOut + Generation ****************
! ******************************************************************************
! can define intermediate variables within the functions to improve readability
! i.e c0 = c_vars_dual(1), dPhi2_dx = dcdx_vars_dual(2)

  function FLUX(c_vars_dual, dcdx_vars_dual) result(Flux_)
    ! (1)   N = -D ∇c
    type(dual_complex), dimension(N)              :: Flux_
    type(dual_complex), dimension(N), intent(in)  :: c_vars_dual, dcdx_vars_dual
    type(dual_complex)                            :: diff_

    type(dual_complex) :: c0
    type(dual_complex) :: dc0dx

    c0       = c_vars_dual(1)
    dc0dx    = dcdx_vars_dual(1)

    Flux_(1) = -diff * dc0dx

  end function FLUX

  function RXN(c_vars_dual) result(Rxn_)
    ! (1)  no reaction
    type(dual_complex), dimension(N)               :: Rxn_
    type(dual_complex), dimension(N), intent(in)   :: c_vars_dual

    type(dual_complex) :: c0
    complex :: i_imag = complex(0, 1)

    c0    = c_vars_dual(1)

    Rxn_(1) = 0.0 - i_imag * freq * c0

  end function RXN

  function ACCUM(c_vars_dual) result(Accum_)
    ! (1) dc0/dt
    type(dual_complex), dimension(N)              :: Accum_
    type(dual_complex), dimension(N), intent(in)  :: c_vars_dual

    type(dual_complex) :: c0

    c0    = c_vars_dual(1)

    Accum_(1) = c0/delT
    Accum_ = 0.0        ! steady_state no accumulation

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
    type(dual_complex), dimension(N)               :: BC_WEST_
    type(dual_complex), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual
    type(dual_complex), dimension(N)               :: flux_temp, rxn_temp

    type(dual_complex) :: c0

    flux_temp = FLUX(c_vars_dual, dcdx_vars_dual)
    c0        = c_vars_dual(1)

    BC_WEST_ (1) = flux_temp(1) - 1e0!0.0

  end function Boundary_WEST

  function Boundary_EAST (c_vars_dual, dcdx_vars_dual) result(BC_EAST_)
    type(dual_complex), dimension(N)               :: BC_EAST_
    type(dual_complex), dimension(N), intent(in)   :: c_vars_dual, dcdx_vars_dual
    type(dual_complex), dimension(N)               :: flux_temp, rxn_temp
    type(dual_complex) :: c0

    flux_temp = FLUX(c_vars_dual, dcdx_vars_dual)
    c0        = c_vars_dual(1)

    ! BC_EAST_(1) = flux_temp(1) - 0.0                                ! Reflective
    ! BC_EAST_(1) = flux_temp(1) - 1e0                                ! Symmetric
    BC_EAST_(1) = c0 - 0.0                                          ! Transmissive

  end function Boundary_EAST
! ******************************************************************************
end module GOV_EQNS

! ******************************************************************************
! ******************************** MAIN PROGRAM ********************************
! ******************************************************************************

program steady
  use user_input
  use variables, only : cprev, delC
  use write_data_mod

  implicit none
  integer :: t1, t2, clock_rate, clock_max
  integer :: it, j, ii, dd
  integer :: maxiter = 1e2
  real    :: tol = 5e-16, error

  call system_clock(t1,clock_rate,clock_max)
  call initial_condition()

  do dd = 1, 1!6
    diff = 10**(-float(dd))

    do ii = 1, 50!-50, 50
      if (ii == 1) then
        freq = 0.0!10**(-5)!10**(float(ii)/10.)
      else
        freq = 10**(float(ii)/10.)
      end if

      error = 1e2
      it = 0

      do while ( (it < maxiter).AND.(error > tol) )

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
        error = maxval(ABS(delC / cprev))
        it = it + 1
      ! ****************************************************************************
      end do

      call write_to_screen(it)
      call write_to_file(it)
      call write_positional_information_to_file(it)
    end do
  end do


  call system_clock(t2,clock_rate,clock_max)
  write ( *, * ) 'Elapsed real time =', real(t2-t1)/real(clock_rate )

end program steady

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

  h = xmax/float(NJ-2)

  ! LINEAR MESH
    delx(1) = 0.0               ! it is common in the finite volume
    delx(2:NJ-1) = h            ! algorithm to set the control
    delx(NJ) = 0.0              ! volumes at the boundaries to zero

! !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
!   ! LOGARITHMIC MESH
!     ! calculate exponential / logarithmic mesh using secant-method
!     ! delX_max is the maximum desired spacing between node points
!     ! produces a grid symmetric about xmax/2
!     delX_max = h*4
!     ! use small values of x0 and x1 (1e-300, 1e-299)
!     ! plotting the function reveals that convergence is more likely when
!     ! the initial guess values are too small as opposed to too big
!                               !  x0,    x1,             tol, max_it
!     delX_init = secant_method(1e-30, 1e-29, delX_max, 1e-15, 100)
!     delX(:NJ/2) = delX_mesh(delX_init, delX_max, NJ/2, xmax/2)
!     delX(NJ/2) = xmax/2.0 - sum(delX(:NJ/2-1))  ! set delX(NJ/2) to exact value necessary
!     delX(NJ/2+1:NJ-1) = delX(NJ/2:2:-1)
! !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!

  ! xx positions are calculated from delX
  xx(1) = 0.0
  do j = 2, NJ
    xx(j) = xx(j-1) + (delX(j) + delX(j-1))/2
  end do

  cprev(1,:) = cbulk

  control_volume_input = trim(geometry)                   ! Define the control volume
  Cntrl_Vol = Control_Volume(control_volume_input)        ! size based on the
  Crx_Area = Cross_Sectional_Area(control_volume_input)   ! specified system geometry

  return                                                  ! is this necessary?

contains

  function delX_mesh(initial_delX, delX_max, mesh_pts, xmax)    result(delX)
    ! defines an exponential grid mesh
    implicit none
    integer, intent(in) :: mesh_pts
    real, intent(in)    :: initial_delX, xmax, delX_max
    real, dimension(mesh_pts) :: delX
    integer             :: j
    real                :: h_log

    h_log = (log10(delX_max) - log10(initial_delX)) / float(mesh_pts-2)

    do j = 1, mesh_pts
      delX(j+1) = log10(initial_delX) + h_log * (j-1)
    end do
    delX(1) = 0.0
    delX(2:) = 10**(delX(2:))

    return
  end function delX_mesh

  function secant_method(x0_init, x1_init, delX_max, tolerance, max_iterations)    result(x2)
    ! """Return the root calculated using the secant method."""
    ! finds the exponential mesh size that fits the system boundaries
    use user_input, only: NJ, xmax
    implicit none
    integer, intent(in) :: max_iterations
    real, intent(in)    :: tolerance
    real, intent(in)    :: x0_init, x1_init, delX_max
    real                :: x0, x1, x2
    integer             :: i
    real                :: x1_x0_diff, f_x0, f_x1

    x0 = x0_init
    x1 = x1_init

    i = 0
    x1_x0_diff = 1e30


    do while ((abs(x1_x0_diff) > tolerance)             .AND. (i < max_iterations)                  .AND. (x0 /= x1)  )

      f_x0 = sum(delX_mesh(x0, delX_max, NJ/2, xmax/2)) - xmax/2
      f_x1 = sum(delX_mesh(x1, delX_max, NJ/2, xmax/2)) - xmax/2
      x1_x0_diff = f_x1 - f_x0

      if (x1_x0_diff == 0.0) return ! x1_x0_diff = 0 produces divide by zero error

      x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0)
      x0 = x1
      x1 = x2

      i = i + 1
    end do

    return

  end function secant_method

end subroutine initial_condition


! include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/Impedance_Complex/sub_auto_fill.f95'
! include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/Impedance_Complex/sub_ABDGXY.f95'
! include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/Impedance_Complex/sub_MATINV.f95'
! include '/Users/nicholasbrady/Documents/Post-Doc/Projects/Fortran/Fortran_DNAD/ExampleModels/Core_Subs_Mods/Impedance_Complex/sub_BAND.f95'
include 'sub_auto_fill.f95'
include 'sub_ABDGXY.f95'
include 'sub_MATINV.f95'
include 'sub_BAND.f95'
