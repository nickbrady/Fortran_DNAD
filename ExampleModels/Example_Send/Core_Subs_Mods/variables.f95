module variables
  ! use user_input
  use number_of_variables
  implicit none

  real, dimension(N, NJ)  :: cprev, delC
  real, dimension(NJ)     :: xx, delX

  ! ABDGXY_VARS
  real, dimension(N, N)     :: A, B, X, Y
  real, dimension(N, 2*N+1) :: D
  real, dimension(N)        :: G

  ! BAND and MATINV variables
  real, dimension(N, N+1, NJ) :: E
  real, dimension(N)          :: ID
  integer                     :: NP1      ! this can be removed

  ! GOV_EQNS variables
  real                   :: alphaW, alphaE, betaW, betaE
  real, dimension(N,N)   :: dW, dE, fW, fE, rj
  real, dimension(N)     :: smG

end module variables
