program lqcpower
  use kinds
  use vars
  use rk
  use omp_lib
  use rhs

  implicit none
   integer(ip) :: cntr =1, tmp, ii, idx, idx2
   character(len=3) :: iname,jname 
   real(wp), dimension(4) :: dydttmp 
   real(wp) :: hmin, x1, x2, h1
   real(wp) :: at,at1,at2,at3,at4,at3tmp,phi,phid,phi2d,rhot,hub, &
               hubd,hub2d,rhoq,rhoq1,rhoq2, pressq
   call system("mkdir -p data")
   call system("rm ./data/*")

   call readparams ()

   allocate(k(numk), q(numk), qdot(numk), rhoqsum(nsteps), pqsum(nsteps), &
            rhoqtmp(numk,nsteps), pqtmp(numk,nsteps), &
            y(2), dydt(2), yback(4), dydtback(4), bckgrnd(4), aa(nsteps), &
            HH(nsteps), HHd(nsteps))

   call initialize ()

   print*, 'Initialization is done!'


   hmin=0.0_wp  
 

!   call OMP_set_num_threads(1)


   print*, 'Going to start iteration!'

   rhoqsum(1:nsteps) = 0.0_wp 
   pqsum(1:nsteps) = 0.0_wp 

do iter=1, maxiter

  print*, 'initial backreaction=', rhoqsum(1)

   write(jname, '(i3.3)') iter
   print*, 'iteration=', iter

     bckgrnd(2) = sqrt(8.0_wp*pi*third*(bckgrnd(4)**2/(2.0_wp*bckgrnd(1)**6) + &
                         potential(bckgrnd(3))  + rhoqsum(1)))


   print*, bckgrnd

     x1=tmin
     x2=x1
     h1=dtime/5.0_wp
     hmin=0.0_wp

     yback(1) = bckgrnd(1)
     yback(2) = bckgrnd(2)
     yback(3) = bckgrnd(3)
     yback(4) = bckgrnd(4)
     dydttmp(:) = 0.0_wp

     call eqnback(yback, dydttmp, 0.0_wp, x1, rhoqsum(1), pqsum(1))


     dtime=(tmax-tmin)/nsteps
     idx = 1
     aa(idx) = yback(1)
     HH(idx) = yback(2)
     HHd(idx) = dydttmp(2)

     print*, "Opening the file to print background data."

     open(iter, file='./data/background_iter_'//jname//'.dat', status='replace',action='write')
       write(iter, 20) x2, yback, rhoqsum(idx),pqsum(idx)

     do while (x2 <= tmax)
       call rk4(yback, dydtback, x2, dtime, k(i),rhoqsum(idx),pqsum(idx))
 
       idx = idx +1 
       aa(idx) = yback(1)
       HH(idx) = yback(2)
       HHd(idx) = dydttmp(2)

       write(iter, 20) x2, yback, rhoqsum(idx),pqsum(idx)
     end do
  
     close(iter)

     print*, 'Background solved for iter=', iter

     rhoqtmp(:,:) = 0.0_wp
     pqtmp(:,:) = 0.0_wp

!$OMP PARALLEL DO private(i,y,dydt,time,tmp,cntr,x1,x2,h1, &
!$OMP                  at,phi,phid,phi2d,rhot,hub,hubd,hub2d,at1,at2, &
!$OMP                  at3,at4,at3tmp,rhoq,rhoq1,rhoq2,iname,idx2) &
!$OMP             shared(tmax,k,vol0,nsteps,numk,tmin,q,qdot,bckgrnd,jname,iter, &
!$OMP                  rhoqsum, pqsum, rhoqtmp, pqtmp,dtime,aa,HH,HHd)
   do i =1, numk

       write(iname, '(i3.3)') i
       open(i, file='./data/powerspec_iter_'//jname//'-kid'//iname//'.dat', status='replace',action='write')

     dydt(:) = (0.0_wp, 0.0_wp)   
 
     time=tmin
    
     x1=tmin
     x2=x1+dtime
     h1=dtime/5.0_wp
     hmin=0.0_wp
      
     idx2 = 1
 
     y(1) = q(i)
     y(2) = qdot(i)

       at = aa(idx2)
       hub= HH(idx2)
       hubd = HHd(idx2)
!       phi = y(3)
!       phid = y(4)/at**3
!       rhot =  half*phid**2+potential(phi)
!
!       hub = y(2)
!       hubd = -4.0_wp*pi*(phid**2)
!       hub2d =8.0_wp*pi*(3.0_wp*hub*phid**2+derpotential(phi)*phid) 
!
!       at1 = at*hub
!       at2 = at*(hubd+hub**2)
!       at3 = at1*(hubd+hub**2) + at*(hub2d+2*hub*hubd)
!
!       phi2d = -3*hub*phid - derpotential(phi)

!       rhoq = k(i)**2/(2.0_wp*pi**2) * & 
!              (abs(y(6))**2+(effpotscalar(at**2,phi,phid)+ &
!                    k(i)**2)*abs(y(5))**2/at**2 - &
!              2.0_wp*renorm (k(i),at,at1,at2,at3,phi,phid,phi2d))

       rhoq = k(i)**2/(4.0_wp*pi**2) * & 
              (abs(y(2))**2+(& !effpotscalar(at**2,phi,phid)+ &
                    k(i)**2)*abs(y(1))**2/at**2 )!- &
              !2.0_wp*renorm (k(i),at,at1,at2,at3,phi,phid,phi2d))

       pressq = k(i)**2/(4.0_wp*pi**2) * & 
              (abs(y(2))**2-( & !effpotscalar(at**2,phi,phid)+ &
                    k(i)**2)*abs(y(1))**2/(3.0_wp*at**2) )!- &
              !2.0_wp*renormP(k(i),at,at1,at2,at3,at4,phi,phid,phi2d))
    
       rhoq = rhoq * rhofactor
       pressq = pressq * rhofactor

         write(i, 10), x1, k(i), y(1), y(2), &
                 2*aimag(y(1)*conjg(y(2))*at**3), rhoq, pressq

       rhoqtmp(i,idx2) = rhoq
       pqtmp(i,idx2) = pressq
 
     cntr=1
!     do while (x2 <= tmax)
     do while (x2 <= tmax)
       call odeint(y,x1,x2,eps,h1,hmin,k(i),aa(idx2),HH(idx2))
 
!     do while (time <= tmax)
!       call rk4(y, dydt, x2, dtime, k(i),rhoqsum(idx),pqsum(idx))

       idx2 = idx2 + 1

       x1=x2
       x2=x1+dtime

       at = aa(idx2)
       hub= HH(idx2)
       hubd = HHd(idx2)
!       phi = y(3)
!       phid = y(4)/at**3
!       rhot =  half*phid**2+potential(phi)
!
!       hub = y(2)
!       hubd = -4.0_wp*pi*(phid**2)
!       hub2d =8.0_wp*pi*(3.0_wp*hub*phid**2+derpotential(phi)*phid) 
!
!       at1 = at*hub
!       at2 = at*(hubd+hub**2)
!       at3 = at1*(hubd+hub**2) + at*(hub2d+2*hub*hubd)
!
!       phi2d = -3*hub*phid - derpotential(phi)

!       rhoq = k(i)**2/(2.0_wp*pi**2) * & 
!              (abs(y(6))**2+(effpotscalar(at**2,phi,phid)+ &
!                    k(i)**2)*abs(y(5))**2/at**2 - &
!              2.0_wp*renorm (k(i),at,at1,at2,at3,phi,phid,phi2d))

       rhoq = k(i)**2/(4.0_wp*pi**2) * & 
              (abs(y(2))**2+(& !effpotscalar(at**2,phi,phid)+ &
                    k(i)**2)*abs(y(1))**2/at**2) !- &
              !2.0_wp*renorm (k(i),at,at1,at2,at3,phi,phid,phi2d))

       pressq = k(i)**2/(4.0_wp*pi**2) * & 
              (abs(y(2))**2-( & !effpotscalar(at**2,phi,phid)+ &
                    k(i)**2)*abs(y(1))**2/(3.0_wp*at**2) ) !- &
              !2.0_wp*renormP(k(i),at,at1,at2,at3,at4,phi,phid,phi2d))

       rhoq = rhoq * rhofactor!* 1.0d-4
       pressq = pressq * rhofactor!* 1.0d-4

       rhoqtmp(i,idx2) = rhoq
       pqtmp(i,idx2) = pressq
 
    
!       tmp = ceiling(real(10**floor(log10(real(cntr,wp)))/4.0_wp,wp))
        tmp = 100
     
       if(mod(cntr, tmp) == 0 ) then 
         write(i, 10), x1, k(i), y(1), y(2), &
                 2*aimag(y(1)*conjg(y(2))*at**3), rhoq, pressq
       end if
       cntr= cntr+1
     end do
  
   end do
  
!$OMP END PARALLEL DO

   do ii=1,numk
     close(ii)
   end do


   do ii=1,nsteps
     rhoqsum(ii) = sum(rhoqtmp(:,ii))
     pqsum(ii) = sum(pqtmp(:,ii))
   end do 



end do 

   123 format(6a20)
   10 format(15es30.15e3)
   20 format(7es30.15e3)


   deallocate (k, q, qdot, rhoqsum, pqsum, rhoqtmp, pqtmp, y, dydt, bckgrnd, &
                aa, HH, HHd, yback, dydtback)

  contains
 
  subroutine readparams ()
    use vars
    implicit none

! Read the input parameters
    open(1,file = 'input.dat', status='old', form = 'formatted', &
                             action = 'read' )
    read(1, nml=input)
    close (1)
  end subroutine readparams


  subroutine initialize ()
    use kinds
    use vars
    
    implicit none

    real(wp) :: tmp1, tmp2, tmp3, tmp4, tmp5, rhotmp
    integer(ip) :: i

    lambda = lambdaabs*sqrt(rhoratio)

! Initilization for the background:
! bckgrnd(1) is the scale factor 'a'
    bckgrnd(1) = (vol0)**(1.0_wp/3.0_wp)

! bckgrnd(3) is the scalar field phi. 
    bckgrnd(3) = phi0

! Now we compute the value of pphi by imposing that the total energy density at
! the bounce is rhomax/rhoratio, and save it in bckgrnd(4).
! bckgrnd(4) = sqrt(2.0_wp*(rhomax/rhoratio - potential(real(bckgrnd(3)))  &
!                      *bckgrnd(1)**(3.0)))

! In general when initial conditions are not given at the bounce then:
    bckgrnd(4) = pphi0

! Check whether the potential is correctly selected:
!    print*, potential(3.0_wp)
!    print*, derpotential(3.0_wp)
!    print*, der2potential(3.0_wp)
     
    print*, 'Potential type selected: ', potential_type


! Now using the initial conditions on p, phi and pphi compute the initial value
! of the connection 'c', say bckgrnd(2).
    tmp1 = bckgrnd(4)**2/(2*bckgrnd(1)**3) + potential(bckgrnd(3))
    tmp2 = asin(abs(sqrt(tmp1*8.0_wp*pi*gamma**2*lambda**2/(3.0_wp))))

    rhotmp = sqrt(8.0_wp*pi*third*(pphi0**2/(2*vol0**2)+potential(phi0)))

    bckgrnd(2) = H0


! Output the initial data on the background:
    print*, 'gamma=', gamma, 'lambda=', lambda
    print*, 'vol0=', vol0
    print*, 'factor=', factor
    print*, 'pini=', bckgrnd(1), 'H0=',bckgrnd(2)
    print*, 'phi0=', bckgrnd(3), 'p_phi=', bckgrnd(4)
    print*, 'H0 from rho =', rhotmp, 'numk=', numk

! Read the initial data for the mode that is generated using the mathematica
! script for various values of wavenumber 'k'.
    open(2, file = './initialdata.dat', status = 'old', &                
                    form ='formatted', action = 'read' )
    do i =1, numk
      read(2, *) tmp1, tmp2, tmp3, tmp4, tmp5
      k(i) = tmp1
      q(i) = tmp2 + zi*tmp3
      qdot(i) = tmp4 + zi*tmp5
    end do

    print*, 'Initial data for the modes read!'

    if (pert_type == 'scalar') then
      print*, 'Scalar perturbation selected!'
    else if (pert_type == 'tensor') then
      print*, 'Tensor perturbation selected!'
    else 
      print*, 'The perturbation can be either "tensor" or "scalar". &
             Check pert_type in the input parameter file.'
    end if

    close(2)

  end subroutine initialize

end program lqcpower
