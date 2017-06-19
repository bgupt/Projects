module rk
  use kinds
  use rhs
 
  implicit none

!  procedure, private :: rhsroutine => eqn

 contains
!********************************************
!**************RK4 Subroutine******************
!********************************************

  subroutine rk4 ( u1, u2, time, dtime, localk, rqloc, pqloc)
    real(wp), dimension(:), intent(inout) :: u1, u2
    real(wp), intent(inout) :: time
    real(wp), intent(in) :: dtime, localk, rqloc, pqloc

!    interface
!      subroutine rhsroutine(u1, u2, localk, localtime)
!        use kinds
!        complex(wp), dimension(:) :: u1, u2
!        real(wp) :: localtime, localk
!      end subroutine rhsroutine
!    end interface

    integer(ip), parameter :: nsteps = 5
    real(wp), dimension(nsteps), parameter :: &
      rk4a = (/ 0.0_wp, &
                -567301805773.0_wp/1357537059087.0_wp, &
                -2404267990393.0_wp/2016746695238.0_wp, &
                -3550918686646.0_wp/2091501179385.0_wp, &
                -1275806237668.0_wp/842570457699.0_wp /)
    real(wp), dimension(nsteps), parameter :: &
      rk4b = (/ 1432997174477.0_wp/9575080441755.0_wp, &
                5161836677717.0_wp/13612068292357.0_wp, &
                1720146321549.0_wp/2090206949498.0_wp, &
                3134564353537.0_wp/4481467310338.0_wp, &
                2277821191437.0_wp/14882151754819.0_wp /)
    real(wp), dimension(nsteps), parameter :: &
      rk4c = (/ 0.0_wp, &
                1432997174477.0_wp/9575080441755.0_wp, &
                2526269341429.0_wp/6820363962896.0_wp, &
                2006345519317.0_wp/3224310063776.0_wp, &
                2802321613138.0_wp/2924317926251.0_wp /)
    complex(wp), dimension(:), allocatable :: resu1, resu2
    real(wp) :: localtime
    integer(ip) :: i, j, n

    n = size(u1)
    allocate ( resu1(n), resu2(n) )

!!$OMP PARALLEL DO private(i) shared(resu1,resu2)
    do i=1, n
      resu1(i) = 0.0_wp
      resu2(i) = 0.0_wp
    end do
!!$OMP END PARALLEL DO

    do j = 1, nsteps
      localtime = time + rk4c(j)*dtime
!      print*,'localtime = ', localtime
      call eqnback(u1, u2, localk, localtime, rqloc, pqloc)
!!$OMP PARALLEL DO private(i) shared(resu1,resu2,dtime,j)
      do i = 1, n
        resu1(i) = rk4a(j) * resu1(i) + dtime * u2(i)
      end do
!!$OMP END PARALLEL DO

!!$OMP PARALLEL DO private(i) shared(u1,u2,resu1,resu2,j)
      do i = 1, n
        u1(i) = u1(i) + rk4b(j) * resu1(i)
      end do
!!$OMP END PARALLEL DO

    end do

    time = time + dtime

    deallocate ( resu1, resu2 )

  end subroutine rk4



!********************************************
!****ODEINT***************************
!********************************************
SUBROUTINE odeint(ystart,x1,x2,eps,h1,hmin,kloc,rqloc,pqloc,phloc,phdloc,bckprt)
  use kinds
integer(ip) :: nbad,nok,nvar
real(wp) :: eps,h1,hmin,x1,x2, kloc, rqloc, pqloc,phloc,phdloc,bckprt
complex(wp), DIMENSION(:):: ystart

!    interface
!      subroutine rhsroutine(u1, u2, localk, localtime)
!        use kinds
!        complex(wp), dimension(:) :: u1, u2
!        real(wp) :: localtime, localk
!      end subroutine rhsroutine
!
!      subroutine rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,kloc)
!         use kinds
!         INTEGER :: n
!         double precision :: eps,hdid,hnext,htry,x,kloc
!         complex(wp), DIMENSION(:)::dydx,y,yscal
!      end subroutine rkqs
!    end interface

integer(ip), PARAMETER :: MAXSTP=100000000,NMAX=50,KMAXX=200
real(wp), parameter :: TINY=1.0d-40

INTEGER :: i,kmax,kount,nstp
real(wp) :: dxsav,h,hdid,hnext,x,xsav
complex(wp), dimension(1:NMAX)::dydx,y,yscal
real(wp), dimension(1:NMAX)::xp
complex(wp), dimension(1:NMAX,1:KMAXX) :: yp
COMMON /path/ kmax,kount,dxsav,xp,yp

x=x1
h=sign(h1,x2-x1)
nok=0
nbad=0
kount=0
nvar=size(ystart)
do i=1,nvar
y(i)=ystart(i)
end do

if (kmax.gt.0) xsav=x-2.0_wp*dxsav

do nstp=1,MAXSTP

        call eqn(y,dydx,kloc,x,rqloc,pqloc,phloc,phdloc,bckprt)
        do i=1,nvar
        yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
	end do

	if(kmax.gt.0)then
	if(abs(x-xsav).gt.abs(dxsav)) then

	if(kount.lt.kmax-1)then
	kount=kount+1
	xp(kount)=x
	
	do i=1,nvar
	yp(i,kount)=y(i)
	end do
	xsav=x
	end if
	end if
	end if
	if((x+h-x2)*(x+h-x1).gt.0.0_wp) h=x2-x

	call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,kloc, rqloc, pqloc,phloc,phdloc,bckprt)
	if(hdid.eq.h)then
	nok=nok+1
	else
	nbad=nbad+1
	end if
	if((x-x2)*(x2-x1).ge.0.0_wp)then

	do i=1,nvar
	ystart(i)=y(i)
	end do
	if(kmax.ne.0)then
	kount=kount+1

	xp(kount)=x
	do i=1,nvar
	yp(i,kount)=y(i)
	end do
	end if
	return

	end if
	if(abs(hnext).lt.hmin) pause 'stepsize smaller than minimum'
	h=hnext
end do
	pause 'too many steps '
	return
END SUBROUTINE odeint



!*******************************************
!************RKQS***************************
!*******************************************
SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,kloc, rqloc, &
                 pqloc,phloc,phdloc,bckprt)
INTEGER :: n
real(wp) :: eps,hdid,hnext,htry,x,kloc, rqloc, pqloc,phloc,phdloc,bckprt
complex(wp), DIMENSION(1:n)::dydx,y,yscal

!INTERFACE
!      subroutine rhsroutine(u1, u2, localk, localtime)
!        use kinds
!        complex(wp), dimension(:) :: u1, u2
!        real(wp) :: localtime, localk
!      end subroutine rhsroutine
!
! subroutine rkck(y,dydx,n,x,h,yout,yerr,kloc)
!  use kinds
! INTEGER(ip) :: n
! complex(wp), DIMENSION(1:n) :: dydx,y,yerr,yout
! DOUBLE PRECISION :: h,x,kloc
! end subroutine rkck
!END INTERFACE

!EXTERNAL eqn
integer, PARAMETER ::NMAX=50


INTEGER :: i
real(wp):: errmax,h,htemp,xnew
complex(wp), dimension(1:NMAX) :: yerr,ytemp

real(wp), PARAMETER :: SAFETY=0.9_wp,PGROW=-.2_wp,PSHRNK=-.25_wp,ERRCON=1.89_wp/10000.0_wp

h=htry

1 call rkck(y,dydx,n,x,h,ytemp,yerr,kloc, rqloc, pqloc,phloc,phdloc, bckprt)

errmax=0.0_wp

do i=1,n
errmax=max(errmax,abs(yerr(i)/yscal(i)))
end do
errmax=errmax/eps

if(errmax.gt.1.0_wp)then

htemp=SAFETY*h*(errmax**PSHRNK)
h=sign(max(abs(htemp),0.1_wp*abs(h)),h)

xnew=x+h
if(xnew.eq.x) then 
print*, 'stepsize and htry =', h, htry
pause 'stepsize underflow rkqs'
end if

goto 1

else

if(errmax.gt.ERRCON)then
hnext=SAFETY*h*(errmax**PGROW)
else

hnext=5.0_wp*h
end if
hdid=h
x=x+h
       do i=1,n
       y(i)=ytemp(i)
        end do
return
end if
END SUBROUTINE rkqs


!**************************************************
!***************RKCK*******************************
!**************************************************
SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,kloc, rqloc, pqloc,phloc,phdloc,bckprt)
INTEGER(ip) :: n
complex(wp),DIMENSION(1:n) :: dydx,y,yerr,yout
real(wp) :: h,x, kloc, rqloc, pqloc,phloc,phdloc, bckprt

!    interface
!      subroutine rhsroutine(u1, u2, localk, localtime)
!        use kinds
!        complex(wp), dimension(:) :: u1, u2
!        real(wp) :: localtime, localk
!      end subroutine rhsroutine
!    end interface


!EXTERNAL eqn
integer, PARAMETER :: NMAX=50
 
INTEGER :: i
complex(wp), dimension(1:nmax) :: ak2,ak3,ak4,ak5,ak6,ytemp
!double precision ::A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
real(wp), PARAMETER :: A2=.2_wp,A3=.3_wp,A4=.6_wp,A5=1.0_wp,A6=.875_wp,B21=.2_wp,B31=3.0_wp/40.0_wp  
real(wp), PARAMETER :: B32=9.0_wp/40.0_wp,B41=.30_wp,B42=-.90_wp,B43=1.20_wp,B51=-11.0_wp/54.0_wp,B52=2.5_wp
real(wp), PARAMETER :: B53=-70.0_wp/27.0_wp,B54=35.0_wp/27.0_wp,B61=1631.0_wp/55296.0_wp,B62=175.0_wp/512.0_wp
real(wp), PARAMETER :: B63=575.0_wp/13824.0_wp,B64=44275.0_wp/110592.0_wp,B65=253.0_wp/4096.0_wp
real(wp), PARAMETER :: C1=37.0_wp/378.0_wp,C3=250.0_wp/621.0_wp,C4=125.0_wp/594.0_wp,C6=512.0_wp/1771.0_wp
real(wp), PARAMETER :: DC1=C1-2825.0_wp/27648.0_wp,DC3=C3-18575.0_wp/48384.0_wp
real(wp), PARAMETER :: DC4=C4-13525.0_wp/55296.0_wp,DC5=-277.0_wp/14336.0_wp,DC6=C6-.250_wp

do i=1,n

ytemp(i)=y(i)+B21*h*dydx(i)
end do
 call eqn(ytemp,ak2,kloc,x+A2*h, rqloc, pqloc,phloc,phdloc,bckprt) 

	do i=1,n
	ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
	end do
 call eqn(ytemp,ak3,kloc,x+A3*h, rqloc, pqloc,phloc,phdloc,bckprt) 

	do i=1,n
	ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))

	end do
 call eqn(ytemp,ak4,kloc,x+A4*h, rqloc, pqloc,phloc,phdloc,bckprt) 

	do i=1,n
	ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
	end do
 call eqn(ytemp,ak5,kloc,x+A5*h, rqloc, pqloc,phloc,phdloc,bckprt) 

	do i=1,n
	ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
	end do 
 call eqn(ytemp,ak6,kloc,x+A6*h, rqloc, pqloc,phloc,phdloc,bckprt) 

	do i=1,n

	yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
	end do
	do i=1,n

	yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
	end do 
return
END subroutine rkck



end module rk
