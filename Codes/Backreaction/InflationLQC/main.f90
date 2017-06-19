program lqcpower
  use kinds
  use vars
  use rk
  use omp_lib
  use rhs

  implicit none
   integer(ip) :: cntr =1, tmp, ii, idx, idx2, ll, iii
   character(len=3) :: iname,jname 
   complex(wp), dimension(4) :: dydttmp
   complex(wp) :: qvac, qdvac, rhochk, betai
   real(wp) :: hmin, x1, x2, h1, rhoqc, pressqc, wq, vq1, vq2
   real(wp) :: at,at1,at2,at3,at4,at3tmp,phi,phid,phi2d,rhot,hub, &
               hubd,hub2d,rhoq,rhoq1,rhoq2, pressq
   real(wp), dimension(:), allocatable :: allofthem
   call system("mkdir -p data")
   call system("rm ./data/*")

   call readparams ()

   allocate(k(numk), q(numk), qdot(numk), rhoqsum(nsteps), pqsum(nsteps), &
            pqsumder(nsteps), rhoqtmp(numk,nsteps), pqtmp(numk,nsteps), &
            allofthem(numk), &
            y(2), dydt(2), yback(4), dydtback(4), bckgrnd(4), aa(nsteps), &
            ysol(6), dydtsol(6), &
            HH(nsteps), HHd(nsteps), aa1(nsteps), aa2(nsteps), aa3(nsteps), &
             aa4(nsteps),phiarray(nsteps),phidarray(nsteps),rqsumtmp(nsteps), pqsumtmp(nsteps), &
             rtmp(numk), ptmp(numk))

   call initialize ()

   print*, 'Initialization is done!'


   hmin=0.0_wp  
 

   call OMP_set_num_threads(16)


   print*, 'Going to start iteration!'

   rhoqsum(1:nsteps) = 0.0_wp 
   pqsum(1:nsteps) = 0.0_wp 
   pqsumder(1:nsteps) = 0.0_wp 

do iter=1, maxiter
 
!   if (iter == 2) then
!     rhofactor = real(rhofactor, wp)
!   else
!     rhofactor = 1.0_wp
!   end if 


  print*, 'initial backreaction=', rhoqsum(1)

   write(jname, '(i3.3)') iter
   print*, 'iteration=', iter

     bckgrnd(2) = (H0/dabs(H0))*bckgrnd(1)*sqrt(8.0_wp*pi*third*(bckgrnd(4)**2/(2.0_wp*bckgrnd(1)**6) + &
                         potential(bckgrnd(3))  + (rhoqsum(1))))


   print*, bckgrnd
   print*, 'H0 =', (bckgrnd(2)/bckgrnd(1))
   print*, '(H0**2) 3/(8 Pi)=', (bckgrnd(2)/bckgrnd(1))**2*3.0_wp/(8.0_wp*pi)

     dtime=(tmax-tmin)/nsteps

     x1=tmin
     x2=x1
     h1=dtime/5.0_wp
     hmin=0.0_wp

     yback(1) = bckgrnd(1)
     yback(2) = bckgrnd(2)
     yback(3) = bckgrnd(3)
     yback(4) = bckgrnd(4)
     dydttmp(:) = 0.0_wp

!     call eqn(yback, dydttmp, 0.0_wp, x1, real(yback(1), wp), real(yback(2), wp)/real(yback(1), wp), &
!                   rhoqsum(1), pqsum(1), 0.0_wp)

     call eqn(yback, dydttmp, 0.0_wp, x1, real(yback(1), wp), real(yback(2), wp)/real(yback(1), wp), &
                   rhoqsum(1), pqsumder(1), 0.0_wp)

     idx = 1
     aa(idx) = real(yback(1), wp)
     aa1(idx) = real(yback(2), wp)
     HH(idx) = aa1(idx)/aa(idx)
     aa2(idx) = real(dydttmp(2), wp)
     phiarray(idx) = real(yback(3), wp)
     phidarray(idx) = real(yback(4), wp)/aa(idx)**3

       rhot = phidarray(idx)**2/2.0_wp + potential(phiarray(idx))
     print*, "Opening the file to print background data."

     open(iter, file='./data/background_iter_'//jname//'.dat', status='replace',action='write')
       write(iter, 20) x2, real(yback, wp), real(dydttmp,wp), rhoqsum(idx),pqsum(idx),pqsumder(idx), rhot

     do while (x2 <= tmax)
!       call odeint(yback,x1,x2,eps,h1,hmin,k(i),aa(idx),HH(idx), rhoqsum(idx), pqsum(idx),0.0_wp)
       call odeint(yback,x1,x2,eps,h1,hmin,k(i),aa(idx),HH(idx), rhoqsum(idx), pqsum(idx),0.0_wp)
!       call rk4(yback, dydtback, x2, dtime, 0.0_wp,rhoqsum(idx),pqsum(idx))

!       call eqn(yback, dydttmp, 0.0_wp, x1, aa(idx),HH(idx),rhoqsum(idx), pqsum(idx), 0.0_wp)
       call eqn(yback, dydttmp, 0.0_wp, x1, aa(idx),HH(idx),rhoqsum(idx), pqsum(idx), 0.0_wp)

       x1=x2
       x2=x1+dtime
 
       idx = idx +1 
       aa(idx) = real(yback(1), wp)
!       aa1(idx) = real(dydttmp(1), wp)
       aa1(idx) = real(yback(2), wp)
       HH(idx) = aa1(idx)/aa(idx)
       aa2(idx) = real(dydttmp(2), wp)
       phiarray(idx) = real(yback(3), wp)
       phidarray(idx) = real(yback(4), wp)/aa(idx)**3


       rhot = phidarray(idx)**2/2.0_wp + potential(phiarray(idx))

       if (mod(idx,outbgrndevery) == 0) then
         write(iter, 20) x2, real(yback, wp), real(dydttmp, wp), rhoqsum(idx),pqsum(idx),pqsumder(idx), rhot
       end if
     end do
  
     close(iter)


!     aa3(1) = (aa2(2)-aa2(1))/dtime
!     aa3(2) = (aa2(3)-aa2(1))/(2.0_wp*dtime)
!     aa3(nsteps) = (aa2(nsteps)-aa2(nsteps-1))/dtime
!     aa3(nsteps-1) = (aa2(nsteps)-aa2(nsteps-2))/(2.0_wp*dtime)
!
!     aa4(1) = (aa3(2)-aa3(1))/dtime
!     aa4(nsteps) = (aa3(nsteps)-aa3(nsteps-1))/dtime

!     aa3(1) = (-aa2(3)+4.0_wp*aa2(2)-3.0_wp*aa2(1))/(2.0_wp*dtime)
!     aa3(2) = (-aa2(4)+4.0_wp*aa2(3)-3.0_wp*aa2(2))/(2.0_wp*dtime)
!     aa3(nsteps) = (3.0_wp*aa2(nsteps)-4.0_wp*aa2(nsteps-1)+aa2(nsteps-2))/(2.0_wp*dtime)
!     aa3(nsteps-1) = (3.0_wp*aa2(nsteps-1)-4.0_wp*aa2(nsteps-2)+aa2(nsteps-3))/(2.0_wp*dtime)

!     aa2(1) = (-aa(4)+4.0_wp*aa(3)-5.0_wp*aa(2)+2.0_wp*aa(1))/(dtime**2)
!     aa2(2) = (-aa(5)+4.0_wp*aa(4)-5.0_wp*aa(3)+2.0_wp*aa(2))/(dtime**2)
!     aa2(nsteps) = (2.0_wp*aa(nsteps)-5.0_wp*aa(nsteps-1)+4.0_wp*aa(nsteps-2)-aa(nsteps-3))/(dtime**2)
!     aa2(nsteps-1) = (2.0_wp*aa(nsteps-1)-5.0_wp*aa(nsteps-2)+4.0_wp*aa(nsteps-3)-aa(nsteps-4))/(dtime**2)

     aa3(1) = (-aa1(4)+4.0_wp*aa1(3)-5.0_wp*aa1(2)+2.0_wp*aa1(1))/(dtime**2)
     aa3(2) = (-aa1(5)+4.0_wp*aa1(4)-5.0_wp*aa1(3)+2.0_wp*aa1(2))/(dtime**2)
     aa3(nsteps) = (2.0_wp*aa1(nsteps)-5.0_wp*aa1(nsteps-1)+4.0_wp*aa1(nsteps-2)-aa1(nsteps-3))/(dtime**2)
     aa3(nsteps-1) = (2.0_wp*aa1(nsteps-1)-5.0_wp*aa1(nsteps-2)+4.0_wp*aa1(nsteps-3)-aa1(nsteps-4))/(dtime**2)

     aa4(1) = (-aa2(4)+4.0_wp*aa2(3)-5.0_wp*aa2(2)+2.0_wp*aa2(1))/(dtime**2)
     aa4(2) = (-aa2(5)+4.0_wp*aa2(4)-5.0_wp*aa2(3)+2.0_wp*aa2(2))/(dtime**2)
     aa4(nsteps) = (2.0_wp*aa2(nsteps)-5.0_wp*aa2(nsteps-1)+4.0_wp*aa2(nsteps-2)-aa2(nsteps-3))/(dtime**2)
     aa4(nsteps-1) = (2.0_wp*aa2(nsteps-1)-5.0_wp*aa2(nsteps-2)+4.0_wp*aa2(nsteps-3)-aa2(nsteps-4))/(dtime**2)

!     do ii = 3, nsteps-2
!       aa2(ii) = (-aa(ii+2)+16.0_wp*aa(ii+1)-30.0_wp*aa(ii)+16.0_wp*aa(ii-1)-aa(ii-2))/(12.0_wp*dtime**2)
!     end do

     do ii = 3, nsteps-2
!       aa3(ii) = (aa2(ii+1)-aa2(ii-1))/(2.0_wp*dtime)
!       aa4(ii) = (aa2(ii+1)-2.0_wp*aa2(ii)+aa2(ii-1))/(dtime**2)
       aa3(ii) = (-aa1(ii+2)+16.0_wp*aa1(ii+1)-30.0_wp*aa1(ii)+16.0_wp*aa1(ii-1)-aa1(ii-2))/(12.0_wp*dtime**2)
       aa4(ii) = (-aa2(ii+2)+16.0_wp*aa2(ii+1)-30.0_wp*aa2(ii)+16.0_wp*aa2(ii-1)-aa2(ii-2))/(12.0_wp*dtime**2)
     end do

     print*, 'Background solved for iter=', iter

!     rhoqtmp(:,:) = 0.0_wp
!     pqtmp(:,:) = 0.0_wp

!$OMP PARALLEL DO private(i,y,dydt,time,tmp,cntr,x1,x2,h1,ysol,dydtsol, &
!$OMP                  at,phi,phid,phi2d,rhot,hub,hubd,hub2d,at1,at2, &
!$OMP                  at3,at4,at3tmp,rhoq,rhoq1,rhoq2,iname,idx2,pressq, &
!$OMP                  rhoqC, pressqC, qvac, qdvac, wq, vq1, vq2, rhochk, betai) &
!$OMP             shared(tmax,k,vol0,nsteps,numk,tmin,q,qdot,bckgrnd,jname,iter, &
!$OMP                  rhoqsum, pqsum, rhoqtmp, pqtmp,dtime,aa,HH,HHd,aa1,aa2)
   do i =1, numk

       write(iname, '(i3.3)') i
       open(i, file='./data/powerspec_iter_'//jname//'-kid'//iname//'.dat', status='replace',action='write')

     dydtsol(:) = (0.0_wp, 0.0_wp)   
 
     time=real(tmin, wp)
    
     x1=real(tmin, wp)
     x2=x1+dtime
     h1=dtime/5.0_wp
     hmin=0.0_wp
      
     idx2 = 1
 
!     y(1) = q(i)
!     y(2) = qdot(i)

       at = aa(idx2)
       at1 = aa1(idx2)
       at2 = aa2(idx2)
       at3 = aa3(idx2)!thirdder(aa, idx2)
       at4 = aa4(idx2)!fourthder(aa, idx2)


!       rhoq = k(i)**2/(2.0_wp*pi**2) * & 
!              (abs(y(6))**2+(effpotscalar(at**2,phi,phid)+ &
!                    k(i)**2)*abs(y(5))**2/at**2 - &
!              2.0_wp*renorm (k(i),at,at1,at2,at3,phi,phid,phi2d))

       rhoqC=  renorm (k(i),at,at1,at2,at3)!,phi,phid,phi2d)

       if (idx2 < 20) then
         pressqC = third*rhoqC
       else 
         pressqC = renormP(k(i),at,at1,at2,at3,at4)!,phi,phid,phi2d)
       end if
       
!       print*, "counter terms:", rhoqC, pressqC


       wq = k(i) !-(k(i)**2/(3.0_wp*at**4*(pressqC-rhoqC)))
       qvac=1.0_wp/(at*sqrt(2.0_wp*wq))
       qdvac=(-zi*wq - at1)*qvac/at

!       if (at1 > 0) then 
!         vq1= -2.0_wp*(sqrt(abs(-(k(i)**2-4.0_wp*at**4*rhoqC*wq+wq**2)))-at1)
!         qdvac=(-zi*wq + vq1/2.0_wp - at1)*qvac/at
!       else 
!         vq1= -2.0_wp*(-sqrt(abs(-(k(i)**2-4.0_wp*at**4*rhoqC*wq+wq**2)))-at1)
!         qdvac=(-zi*wq + vq1/2.0_wp - at1)*qvac/at
! 
!       end if
 
       betai = beta0*exp(-(k(i)/10.0_wp)**16)!*exp(zi*2.0_wp*pi*k(i))
!       betai = beta0*exp(-(k(i)/10.0_wp)**4)
!       betai = beta0*exp(-(k(i)/1.5_wp))
!       betai = beta0/k(i)**4.5_wp!*exp(-(k(i)/10.0_wp)**16)

!       y(1) = sqrt(1.0_wp+betai**2)*qvac + conjg(betai*qvac)
!       y(2) = sqrt(1.0_wp+betai**2)*qdvac + conjg(betai*qdvac)

       ysol(1) = bckgrnd(1)
       ysol(2) = bckgrnd(2)
       ysol(3) = bckgrnd(3)
       ysol(4) = bckgrnd(4)
       ysol(5) = sqrt(1.0_wp+betai**2)*qvac + conjg(betai*qvac)
       ysol(6) = sqrt(1.0_wp+betai**2)*qdvac + conjg(betai*qdvac)
!       ysol(5) = q(i)!sqrt(1.0_wp+betai**2)*qvac + conjg(betai*qvac)
!       ysol(6) = qdot(i)!sqrt(1.0_wp+betai**2)*qdvac + conjg(betai*qdvac)

       if(i == 1) then
         print*, 'ysol initial=', ysol
         print*, 'wq:', wq
         print*, 'vq1:', vq1
         print*, 'vq1 theory:', -2.0_wp*(-sqrt(abs(-(k(i)**2-4.0_wp*at**4*rhoqC*wq+wq**2)))-at1)
         print*, 'rhoC', rhoqC
         print*, 'at1:', at1
         print*, 'k(i):', k(i)
         print*, " "
       end if


!       print*, 'freq:', wq, 'at:', at, 'at1:', at1, 'vq1:', vq1
!       print*, 'at2:', at2, 'at3:', at3, 'at4:', at4
!
!STOP
!
!       rhochk =  k(i)**2/(4.0_wp*pi**2) * &
!              (abs(y(2))**2-abs(qdvac)**2+(& !effpotscalar(at**2,phi,phid)+ &
!                    k(i)**2)*(abs(y(1))**2-abs(qvac)**2)/at**2)

       rhochk = -zi * at**3 * ( ysol(6)*qvac - ysol(5)*qdvac)

       rhoq = k(i)**2/(4.0_wp*pi**2) * & 
              (abs(ysol(6))**2+(& !effpotscalar(at**2,phi,phid)+ &
                    k(i)**2)*abs(ysol(5))**2/at**2 - 2.0_wp*rhoqC)

       pressq = k(i)**2/(4.0_wp*pi**2) * & 
              (abs(ysol(6))**2-( & !effpotscalar(at**2,phi,phid)+ &
                    k(i)**2)*abs(ysol(5))**2/(3.0_wp*at**2) - 2.0_wp*pressqC)

!       rhoq = k(i)**2/(4.0_wp*pi**2) * & 
!              ((abs(ysol(6))**2-abs(qdvac)**2)+(& !effpotscalar(at**2,phi,phid)+ &
!                    k(i)**2)*(abs(ysol(5))**2-abs(qvac)**2)/at**2) !- 2.0_wp*rhoqC)
!
!       pressq = k(i)**2/(4.0_wp*pi**2) * & 
!              ((abs(ysol(6))**2-abs(qdvac)**2)-( & !effpotscalar(at**2,phi,phid)+ &
!                    k(i)**2)*(abs(ysol(5))**2-abs(qvac)**2)/(3.0_wp*at**2)) !- 2.0_wp*pressqC)
!       rhoq = k(i)**2/(4.0_wp*pi**2) * & 
!              ((abs(y(2))**2-abs(qdvac)**2)+(& !effpotscalar(at**2,phi,phid)+ &
!                    k(i)**2)*(abs(y(1))**2-abs(qvac)**2)/at**2)! - rhoqC)
!
!       pressq = k(i)**2/(4.0_wp*pi**2) * & 
!              ((abs(y(2))**2-abs(qdvac)**2)-( & !effpotscalar(at**2,phi,phid)+ &
!                    k(i)**2)*(abs(y(1))**2-abs(qvac)**2)/(3.0_wp*at**2))! - pressqC)

       rhoq = rhoq * rhofactor
       pressq = pressq * rhofactor

       rhoqtmp(i,idx2) = rhoq
       pqtmp(i,idx2) = pressq

 
         write(i, 10), x1, k(i), ysol(5), ysol(6), at, at1, at2, at3, at4, &
              2*aimag(ysol(5)*conjg(ysol(6))*at**3), rhoq, 2.0_wp*rhoqC, pressq, 2.0_wp*pressqC, &
               qvac, qdvac, rhochk

     cntr=1
!     do while (x2 <= tmax)

     do while (x2 <= tmax)
!       call odeint(ysol,x1,x2,eps,h1,hmin,k(i),aa(idx2),HH(idx2),rhoqsum(idx),pqsum(idx),1.0_wp)
       call odeint(ysol,x1,x2,eps,h1,hmin,k(i),aa(idx2),HH(idx2),rhoqsum(idx),pqsumder(idx),1.0_wp)
 
!     do while (time <= tmax)
!       call rk4(y, dydt, x2, dtime, k(i),rhoqsum(idx),pqsum(idx))

       idx2 = idx2 + 1

       x1=x2
       x2=x1+dtime

       at = aa(idx2)
       at1 = aa1(idx2)
       at2 = aa2(idx2)
       at3 = aa3(idx2)!thirdder(aa, idx2)
       at4 = aa4(idx2)!fourthder(aa, idx2)

!       rhoq = k(i)**2/(2.0_wp*pi**2) * & 
!              (abs(y(6))**2+(effpotscalar(at**2,phi,phid)+ &
!                    k(i)**2)*abs(y(5))**2/at**2 - &
!              2.0_wp*renorm (k(i),at,at1,at2,at3,phi,phid,phi2d))

       rhoqC =  renorm (k(i),at,at1,at2,at3)!,phi,phid,phi2d)


       pressqC = renormP(k(i),at,at1,at2,at3,at4)!,phi,phid,phi2d)


!       rhoqtmp(i,idx2) = rhoq
!       pqtmp(i,idx2) = pressq
 
!       wq = -k(i)**2/(3.0_wp*at**4*(pressqC-rhoqC))
!       vq1= -2.0_wp*(sqrt(-(k(i)**2-4.0_wp*at**4*rhoqC*wq+wq**2))-at1)
!       qvac=1.0_wp/(at*sqrt(2.0_wp*wq))
!       qdvac=(-zi*wq + vq1/2.0_wp - at1)*qvac/at

       wq = k(i) !-(k(i)**2/(3.0_wp*at**4*(pressqC-rhoqC)))
       qvac=1.0_wp/(at*sqrt(2.0_wp*wq))
       qdvac=(-zi*wq - at1)*qvac/at
    
!       rhochk =  k(i)**2/(4.0_wp*pi**2) * &
!              (abs(y(2))**2-abs(qdvac)**2+(& !effpotscalar(at**2,phi,phid)+ &
!                    k(i)**2)*(abs(y(1))**2-abs(qvac)**2)/at**2)

        rhochk = -zi * at**3 * ( y(6)*qvac - y(5)*qdvac)

       rhoq = k(i)**2/(4.0_wp*pi**2) * & 
              (abs(ysol(6))**2+(& !effpotscalar(at**2,phi,phid)+ &
                    k(i)**2)*abs(ysol(5))**2/at**2 - 2.0_wp*rhoqC)

       pressq = k(i)**2/(4.0_wp*pi**2) * & 
              (abs(ysol(6))**2-( & !effpotscalar(at**2,phi,phid)+ &
                    k(i)**2)*abs(ysol(5))**2/(3.0_wp*at**2) - 2.0_wp*pressqC)

!       rhoq = k(i)**2/(4.0_wp*pi**2) * & 
!              ((abs(y(2))**2-abs(qdvac)**2)+(& !effpotscalar(at**2,phi,phid)+ &
!                    k(i)**2)*(abs(y(1))**2-abs(qvac)**2)/at**2)! - rhoqC)
!
!       pressq = k(i)**2/(4.0_wp*pi**2) * & 
!              ((abs(y(2))**2-abs(qdvac)**2)-( & !effpotscalar(at**2,phi,phid)+ &
!                    k(i)**2)*(abs(y(1))**2-abs(qvac)**2)/(3.0_wp*at**2))! - pressqC)

       rhoq = rhoq * rhofactor
       pressq = pressq * rhofactor

       rhoqtmp(i,idx2) = rhoq
       pqtmp(i,idx2) = pressq
!       tmp = ceiling(real(10**floor(log10(real(cntr,wp)))/4.0_wp,wp))
        tmp = 2000
     
       if(mod(cntr, outevery) == 0 ) then 
         write(i, 10), x1, k(i), ysol(5), ysol(6), at, at1, at2, at3, at4, &
              2*aimag(ysol(5)*conjg(ysol(6))*at**3), rhoq, 2.0_wp*rhoqC, pressq, 2.0_wp*pressqC, &
              qvac, qdvac, rhochk
       end if
       cntr= cntr+1
     end do
  
   end do
  
!$OMP END PARALLEL DO
   
   Print*, "Begin summing rho and p for iteration=", iter


!   do ii=1,numk
!     close(ii)
!   end do
!
!   do ii =1, numk
!       write(iname, '(i3.3)') ii
!       open(ii, file='./data/powerspec_iter_'//jname//'-kid'//iname//'.dat', &
!                status='old', form='formatted',action='read')
!   end do

!   do ll =1, nsteps
!     do i =1, numk
!       read(i, *), allofthem(:14)
!       rtmp(i) = allofthem(12)
!       ptmp(i) = allofthem(14)
!     end do
!      
!     rhoqsum(ll) = dk*sum(rtmp(:))
!     pqsum(ll) = dk*sum(ptmp(:))
!   end do
!
!   do ii=1,numk
!     close(ii)
!   end do
!
!   Print*, "Done summing rho and p for iteration=", iter
!     rhoqsum = dk*sum(rhoqtmp, dim=1)
!     pqsum = dk*sum(pqtmp, dim=1)


   do ii=1,nsteps
!     rqsumtmp(:) = rhoqtmp(:,ii)
!     pqsumtmp(:) = pqtmp(:,ii)
     rhoqsum(ii) = dk*sum(rhoqtmp(2:numk,ii)+rhoqtmp(1:numk-1,ii))*half
!     pqsum(ii) = rhoqsum(ii)*third!dabs(dk*sum(pqtmp(2:numk,ii)+pqtmp(1:numk-1,ii))*half)
     pqsum(ii) = dabs(dk*sum(pqtmp(2:numk,ii)+pqtmp(1:numk-1,ii))*half)
!     pqsum(ii) = dk*sum(pqsumtmp)
   end do 

   do ii=1,nsteps
     if ((ii >= 3) .and. (ii <=nsteps-2)) then
       pqsumder(ii) = -rhoqsum(ii) -(-rhoqsum(ii+2)+8.0_wp*rhoqsum(ii+1) - &
                      8.0_wp*rhoqsum(ii-1) + rhoqsum(ii-2))/(12.0_wp*dtime)/&
                      (3.0_wp*aa1(ii)/aa(ii))
      elseif (ii<3) then
       pqsumder(ii) = -rhoqsum(ii) &
                      -(-3.0_wp*rhoqsum(ii)+4.0_wp*rhoqsum(ii+1)-rhoqsum(ii+2))/&
                      (2.0_wp*dtime)/(3.0_wp*aa1(ii)/aa(ii))
      elseif (ii>nsteps-2) then
       pqsumder(ii) = -rhoqsum(ii) &
                      -(3.0_wp*rhoqsum(ii)-4.0_wp*rhoqsum(ii-1)+rhoqsum(ii-2))/&
                      (2.0_wp*dtime)/(3.0_wp*aa1(ii)/aa(ii))
      end if
   end do

!   ii=1
!
!   do while (ii < nsteps-30)
!     if (pqsum(ii) < 0.0_wp) then
!      do iii=ii,ii+100
!        pqsum(iii) = sum(pqsum(iii-25:iii+25))/50.0_wp
!      end do
!      ii=ii+100
!     else  
!      ii=ii+1
!     end if
!   end do


end do 

   123 format(6a20)
   10 format(22es30.15e3)
   20 format(13es30.15e3)


   deallocate (k, q, qdot, rhoqsum, pqsum, pqsumder, rhoqtmp, pqtmp, y, dydt, bckgrnd, &
                aa, HH, HHd, yback,dydtback,phiarray,phidarray,rqsumtmp,pqsumtmp, &
                 rtmp, ptmp)

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

    bckgrnd(2) = bckgrnd(1)*H0


! Output the initial data on the background:
    print*, 'gamma=', gamma, 'lambda=', lambda
    print*, 'vol0=', vol0
    print*, 'factor=', factor
    print*, 'aini=', bckgrnd(1), 'adini=',bckgrnd(2)
    print*, 'phi0=', bckgrnd(3), 'p_phi=', bckgrnd(4)
    print*, 'H0 from adot and a =', bckgrnd(2)/bckgrnd(1) 
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


  real(wp) function thirdder(array, id)
    use vars

    implicit none
    real(wp), dimension(:), intent(in) :: array
    integer(ip), intent(in) :: id


    if (id < 4) then
      thirdder = half/(dtime**3)* (-3.0_wp*array(id+4) + 14.0_wp*array(id+3) &
                 -24.0_wp*array(id+2) + 18.0_wp*array(id+1) - 5.0_wp*array(id))
    elseif (id > nsteps-4) then
      thirdder = -half/(dtime**3)* (-3.0_wp*array(id-4) + 14.0_wp*array(id-3) &
                 -24.0_wp*array(id-2) + 18.0_wp*array(id-1) - 5.0_wp*array(id))
    else
      thirdder = 1.0_wp/(8.0_wp*dtime**3) *(-array(id+3) + 8.0_wp*array(id+2) &
                   -13.0_wp*array(id+1) + 13.0_wp*array(id-1) &
                   -8.0_wp*array(id-2) + array(id-3) )
    end if
     
    
  end function thirdder


  real(wp) function fourthder(array, id)
    use vars

    implicit none
    real(wp), dimension(:), intent(in) :: array
    integer(ip), intent(in) :: id


    if (id < 4) then
      fourthder = 1.0_wp/(dtime**4)* (-2.0_wp*array(id+5)+11.0_wp*array(id+4)-24.0_wp*array(id+3) &
                 +26.0_wp*array(id+2) - 14.0_wp*array(id+1) + 3.0_wp*array(id))
    elseif (id > nsteps-4) then
      fourthder = 1.0_wp/(dtime**4)* (-2.0_wp*array(id-5)+11.0_wp*array(id-4)-24.0_wp*array(id-3) &
                 +26.0_wp*array(id-2) - 14.0_wp*array(id-1) + 3.0_wp*array(id))
    else
      fourthder = 1.0_wp/(6.0_wp*dtime**4) *(-array(id+3) + 12.0_wp*array(id+2) &
                   -39.0_wp*array(id+1) + 56.0_wp*array(id) - 39.0_wp*array(id-1) &
                   -12.0_wp*array(id-2) - array(id-3) )    
    end if
     
    
  end function fourthder



  subroutine gaulegf(x1, x2, x, w, n)
    implicit none
    integer, intent(in) :: n
    real(wp), intent(in) :: x1, x2
    real(wp), dimension(n), intent(out) :: x, w
    integer :: i, j, m
    real(wp) :: p1, p2, p3, pp, xl, xm, z, z1
    real(wp), parameter :: eps=3.d-14
        
    m = (n+1)/2
    xm = 0.5d0*(x2+x1)
    xl = 0.5d0*(x2-x1)
    do i=1,m
      z = cos(3.141592654d0*(i-0.25d0)/(n+0.5d0))
      z1 = 0.0
      do while(abs(z-z1) .gt. eps)
        p1 = 1.0d0
        p2 = 0.0d0
        do j=1,n
          p3 = p2
          p2 = p1
          p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
        end do
        pp = n*(z*p1-p2)/(z*z-1.0d0)
        z1 = z
        z = z1 - p1/pp
      end do
      x(i) = xm - xl*z
      x(n+1-i) = xm + xl*z
      w(i) = (2.0d0*xl)/((1.0d0-z*z)*pp*pp)
      w(n+1-i) = w(i)
    end do
  end subroutine gaulegf

end program lqcpower
