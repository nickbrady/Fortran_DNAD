! TEST MATINV

module number_of_variables
  implicit none
  integer, parameter :: N = 3
  integer, parameter :: NJ = 42
end module number_of_variables

include '../CoreFortran_Subroutines_Modules/variables.f95'
include '../CoreFortran_Subroutines_Modules/subroutine_MATINV.f95'


program test
  use number_of_variables
  use variables

  implicit none

  integer :: nn
  real :: DETERM

  B = reshape( (/ 0.7866495,  0.59600742, 0.67644582, &
                    0.38397028, 0.06693267, 0.10238761, &
                    0.62778101, 0.29118293, 0.0223777 /), &
                  shape(B), order=(/2, 1/) )

  D(:,1:N) = reshape( (/ 0.33298009, 0.3387619,  0.95379304, &
                         0.63159401, 0.56107099, 0.12250367, &
                         0.59474099, 0.87774137, 0.35123293 /), &
                      shape(D(:,1:N)), order=(/2, 1/) )

  do nn = 1, N
    print*, B(nn,:)
  end do

  call MATINV(N, N, DETERM)

  do nn = 1, N
    print*, B(nn,:)
  end do

  print*,

  do nn = 1, N
    print*, D(nn,1:N)
  end do

end program test
