! Dependencies
! user_input > really_long_module > MOD_Echem > CircleArea

module user_input
  integer, PARAMETER :: N = 2
  real :: PI = 4.0D0*ATAN(1.0D0) ! pi - Geometric constant

end module user_input


! module really_long_module
!     use user_input, only: N
!     implicit none
!
! 		real :: test_number = N
!
! end module really_long_module


module useful_module
  use user_input
  use really_long_module
  SAVE

  real :: useful_output

end module useful_module


program main
  use user_input
  use really_long_module
  use useful_module

	useful_output = 2 * N

	print*, useful_output

end program main

!### Would like to put really_long_module at the end of file
module really_long_module
    use user_input, only: N
    implicit none

		real :: test_number = 0.5

end module really_long_module


module A
    use user_input
    implicit none

end module A

module B
    use H
		use F
    implicit none

end module B

module C
    use B
		use D
    implicit none

end module C

module D
    use A
		use F
    implicit none

end module D

module E
    use F
    implicit none

end module E

module F
    use user_input
    implicit none

end module F

module G
    use D
    implicit none

end module G

module H
    use E
    implicit none

end module H
