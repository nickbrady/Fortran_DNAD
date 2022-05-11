! Implemented a method for a log-space grid mesh
! This is useful because frequently the change profiles have large
! gradients near interfaces and this method allows us to capture 
! these sharp gradients without adding excessive node points, which
! would hinder computational performance

subroutine initial_condition()
    use user_input
    use variables
  
    implicit none
    real    :: h
    integer :: j
    real    :: delX_init, delX_max
  
    character(len=:), allocatable :: control_volume_input
  
    h = xmax/float(nj-2)
    
  ! LINEAR MESH
    delx(1) = 0.0               ! it is common in the finite volume
    delx(2:NJ-1) = h            ! algorithm to set the control
    delx(NJ) = 0.0              ! volumes at the boundaries to zero
    ! do j = 1, NJ
    !   print*, j, delX(j)
    ! end do
  
  !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
  ! LOGARITHMIC MESH
    ! calculate exponential / logarithmic mesh using secant-method
    ! delX_max is the maximum desired spacing between node points
    ! produces a grid symmetric about xmax/2
    delX_max = h*4
    ! use small values of x0 and x1 (1e-300, 1e-299)
    ! plotting the function reveals that convergence is more likely when
    ! the initial guess values are too small as opposed to too big
                              !  x0,    x1,             tol, max_it
    delX_init = secant_method(1e-300, 1e-299, delX_max, 1e-15, 100)
    delX(:NJ/2) = delX_mesh(delX_init, delX_max, NJ/2, xmax/2)
    delX(NJ/2) = xmax/2.0 - sum(delX(:NJ/2-1))  ! set delX(NJ/2) to exact value necessary
    delX(NJ/2+1:NJ-1) = delX(NJ/2:2:-1)
  !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
  
    ! xx positions are calculated from delX
    xx(1) = 0.0
    do j = 2, NJ
      xx(j) = xx(j-1) + (delX(j) + delX(j-1))/2
    end do
  
    do j = 1, NJ
      cprev(1,j) = c_LiTFSI_bulk
      cprev(2,j) = vel_0
      cprev(3,j) = 0.0
    end do
  
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
  
  
      do while ((abs(x1_x0_diff) > tolerance) &
        &       .AND. (i < max_iterations)    &
                .AND. (x0 /= x1)  )
        
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