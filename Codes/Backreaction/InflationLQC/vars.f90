module vars
  
  use kinds

  implicit none

  real(wp), parameter :: kconst = 2.0_wp &
                                 /(3.0_wp*sqrt(3.0_wp*sqrt(3.0_wp))), gamma=0.2375_wp
  real(wp), parameter :: half = 0.5_wp
  real(wp), parameter :: eps = 1.0d-16 ! eps= 1.0e-10
  real(wp), parameter :: third = 1.0_wp/3.0_wp
  real(wp), parameter :: pi = acos(-1.0_wp)
  real(wp), parameter :: lambdaabs = sqrt(4.0_wp*sqrt(3.0_wp)*pi*gamma)
  real(wp), parameter :: twelvepi = 12.0_wp*pi
  real(wp), parameter :: sqsixteenpi3 = sqrt(16.0_wp*pi/3.0_wp)
  real(wp), parameter :: rhomax = 3.0_wp/(8.0_wp*pi*gamma**2*lambdaabs**2)
  real(wp), parameter :: factor = (1.0_wp/kconst) &
                                  *(8.0_wp*pi*gamma/6.0_wp)**(3.0_wp/2.0_wp)
  complex(wp), parameter :: zi = (0.0_wp, 1.0_wp), zzero = (0.0_wp, 0.0_wp)
  integer(ip) ::  n, sgn, i, j, rhs_type, numk, nsteps, outevery, iter, maxiter, &
                    outbgrndevery, inidatatype
 
  real(wp) :: mass, vol0, rhoratio, phi0, pphi0, tmax, tmin, direc, H0, cini, &
                        lambda, mu, rhofactor=1.0_wp, beta0, kmin, kmax
  character(6) :: pert_type
  character(11) :: potential_type
  real(wp) :: time, dtime, dk

  complex(wp), dimension(:), allocatable :: yback, dydtback, y, yout, dydt, ystart
  real(wp), dimension(:), allocatable :: k, rhoqsum, pqsum, pqsumder, aa, HH, HHd, & 
                                         bckgrnd, aa1, aa2, &
                                         aa3,  aa4, phiarray, phidarray, &
                                           rqsumtmp, pqsumtmp, rtmp, ptmp
  real(wp), dimension(:, :), allocatable :: rhoqtmp, pqtmp
  complex(wp), dimension(:), allocatable :: q, qdot, ysol, dydtsol

  real(wp) :: vol, pressure, rho, Vphi1, Vphi2

  namelist /input/ mass, vol0, H0, cini, rhoratio, rhofactor, dk, maxiter, phi0, pphi0, pert_type, potential_type, tmax, tmin, numk, outevery, outbgrndevery, nsteps, direc, mu, inidatatype, beta0, kmin, kmax
  
end module vars

