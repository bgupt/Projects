module rhs
  use vars
  use kinds

 contains

  real(wp) function potential (x)
    use kinds
    use vars

    implicit none
    real(wp), intent(in) :: x

!    print*, potential_type

    if (potential_type == 'massive') then
      potential = massive (x)
    else if (potential_type == 'starobinsky') then
      potential = starobinsky (x)
    else 
    print*, 'Wrong potential type! Aborting evolution...'
    STOP
    end if

  end function potential


  real(wp) function derpotential (x)
    use kinds
    use vars

    implicit none
    real(wp), intent(in) :: x

!    print*, potential_type

    if (potential_type == 'massive') then
      derpotential = dermassive (x)
    else if (potential_type == 'starobinsky') then
      derpotential = derstarobinsky (x)
    else 
    print*, 'Wrong potential type! Aborting evolution...'
    STOP
    end if

  end function derpotential


  real(wp) function der2potential (x)
    use kinds
    use vars

    implicit none
    real(wp), intent(in) :: x

!    print*, potential_type

    if (potential_type == 'massive') then
      der2potential = der2massive (x)
    else if (potential_type == 'starobinsky') then
      der2potential = der2starobinsky (x)
    else 
    print*, 'Wrong potential type! Aborting evolution...'
    STOP
    end if

  end function der2potential

  real(wp) function massive (x)
    use kinds
    use vars
  
    implicit none
    real(wp), intent(in) :: x

    massive = half*(mass*x)**2
  end function massive

  real(wp) function starobinsky (x)
    use kinds
    use vars
  
    implicit none
    real(wp), intent(in) :: x

    starobinsky = 3.0_wp*mass**2/(32.0_wp*pi) * &
                  (1.0_wp - exp(-sqsixteenpi3*x))**2
  end function starobinsky


  real(wp) function dermassive (x)
    use kinds
    use vars
  
    implicit none
    real(wp), intent(in) :: x

    dermassive = x*mass**2
  end function dermassive

  real(wp) function derstarobinsky (x)
    use kinds
    use vars
  
    implicit none
    real(wp), intent(in) :: x

    derstarobinsky = exp(-2.0_wp*sqsixteenpi3*x) * &
                     (exp(sqsixteenpi3*x)-1.0_wp)*mass**2/sqsixteenpi3
                  
  end function derstarobinsky

  real(wp) function der2massive (x)
    use kinds
    use vars
  
    implicit none
    real(wp), intent(in) :: x

    der2massive = mass**2
  end function der2massive

  real(wp) function der2starobinsky (x)
    use kinds
    use vars
  
    implicit none
    real(wp), intent(in) :: x

    der2starobinsky = exp(-2.0_wp*sqsixteenpi3*x) - &
                      exp(-sqsixteenpi3*x)*(1.0_wp-exp(-sqsixteenpi3*x))
                     
    der2starobinsky = der2starobinsky*mass**2 
                  
  end function der2starobinsky


  subroutine eqnback(yin, yderiv, klocal, tphi, rhoqlocal, pqlocal)
    use kinds
    use vars
    implicit none

    real(wp) :: y1, y2, y3, y4, phidloc, rhoqlocal, pqlocal
    complex(wp) :: y5, y6 
    real(wp), intent(in) :: klocal, tphi
    real(wp), dimension(:) :: yin, yderiv

! y1: scale factor 'a'
    y1= real(yin(1))
! y2: hubble parameter 'h'
    y2= real(yin(2))
! y3: scalar field 'phi'
    y3= real(yin(3))
! y4: momentum of the scalar field 'pphi'
    y4= real(yin(4))
   
   phidloc = y4/y1**3
 
! adot
    yderiv(1) = y2 !* sqrt(8.0_wp*pi*third*(phidloc**2/2.0_wp &
                      !    +potential(y3)+rhoqlocal))

! hdot
    yderiv(2) = -4.0_wp*y1*pi*third*(2.0_wp*phidloc**2 - 2.0_wp*potential(y3) &
                                   +rhoqlocal+3.0_wp*pqlocal) 

! phidot 
    yderiv(3) = y4/y1**3

! pphidot
    yderiv(4) = -derpotential(y3)*y1**3
    
  end subroutine eqnback



  subroutine eqn(yin, yderiv, klocal, tphi, aalocal, HHlocal, philoc, phidloc, backpert)
    use kinds
    use vars
    implicit none

    real(wp) :: y1, y2, y3, y4, phidloctmp, aalocal, HHlocal, philoc, phidloc, &  
                 backpert, rhoqlocal, pqlocal
    complex(wp) :: y5, y6 
    real(wp), intent(in) :: klocal, tphi
    complex(wp), dimension(:) :: yin, yderiv
   
    if (backpert == 0.0_wp) then 

       ! y1: scale factor 'a'
           y1= real(yin(1))
       ! y2: hubble parameter 'h'
           y2= real(yin(2))
       ! y3: scalar field 'phi'
           y3= real(yin(3))
       ! y4: momentum of the scalar field 'pphi'
           y4= real(yin(4))
       
          phidloctmp = y4/y1**3
      
          rhoqlocal = philoc
          pqlocal = phidloc
 
       ! adot
           yderiv(1) = y2 !* sqrt(8.0_wp*pi*third*(phidloc**2/2.0_wp &
!           yderiv(1) = y1*sqrt(8.0_wp*pi*third*(phidloctmp**2/2.0_wp &
!                                 +potential(y3)+rhoqlocal))
       
       ! hdot
!           yderiv(2) = -4.0_wp*y1*pi*third*(2.0_wp*phidloctmp**2 - &
!                             2.0_wp*potential(y3) +rhoqlocal+3.0_wp*pqlocal) 

           yderiv(2) = -4.0_wp*pi*y1*(phidloctmp**2+rhoqlocal+pqlocal) &
                         + 8.0_wp*pi*third*y1*(phidloctmp**2/2.0_wp &
                                 +potential(y3)+rhoqlocal)
       
       ! phidot 
           yderiv(3) = y4/y1**3
       
       ! pphidot
           yderiv(4) = -derpotential(y3)*y1**3

     elseif (backpert == 1.0_wp) then

       ! y1: scale factor 'a'
           y1= real(yin(1))
       ! y2: hubble parameter 'h'
           y2= real(yin(2))
       ! y3: scalar field 'phi'
           y3= real(yin(3))
       ! y4: momentum of the scalar field 'pphi'
           y4= real(yin(4))
       
          phidloctmp = y4/y1**3
      
          rhoqlocal = philoc
!          pqlocal = rhoqlocal*third !phidloc
          pqlocal = phidloc
 
       ! adot
           yderiv(1) = y2 !* sqrt(8.0_wp*pi*third*(phidloc**2/2.0_wp &
!           yderiv(1) = y1 * sqrt(8.0_wp*pi*third*(phidloctmp**2/2.0_wp &
!                                 +potential(y3)+rhoqlocal))
       
       ! hdot
!           yderiv(2) = -4.0_wp*y1*pi*third*(2.0_wp*phidloctmp**2 - &
!                             2.0_wp*potential(y3) +rhoqlocal+3.0_wp*pqlocal) 

            yderiv(2) = -4.0_wp*pi*y1*(phidloctmp**2+rhoqlocal+pqlocal) &
                         + 8.0_wp*pi*third*y1*(phidloctmp**2/2.0_wp &
                                 +potential(y3)+rhoqlocal)
       
       ! phidot 
           yderiv(3) = y4/y1**3
       
       ! pphidot
           yderiv(4) = -derpotential(y3)*y1**3


       ! y5: the mode function 'q'
           y5= yin(5)
       ! y6: derivative of the mode function 'qdot'
           y6= yin(6)
          
!    print*, 'yin=', yin
!STOP
       
       ! qdot  
           yderiv(5) = y6
       
       ! qdot_dot 
!           yderiv(6) = -3.0_wp*(HHlocal)*y6 - (klocal**2 + &
!                       effpotscalar(aalocal, philoc, phidloc) ) &
!                         *y5/aalocal**2
!           yderiv(6) = -3.0_wp*(y2/y1)*y6 - (klocal**2 + &
!                       effpotscalar(aalocal, philoc, phidloc) ) &
!                         *y5/y1**2
           yderiv(6) = -3.0_wp*(yderiv(1)/y1)*y6 - (klocal**2 + &
                       effpotscalar(aalocal, philoc, phidloc) ) &
                         *y5/y1**2

!       print*, 'yderive q:', yderiv(5), yderiv(6)
!STOP

      else 
        print*, 'Wrong backpert while calling Odeint!'
      end if
    
  end subroutine eqn

!  subroutine eqnpert(yin, yderiv, klocal, tphi, aalocal, HHlocal, backpert)
!    use kinds
!    use vars
!    implicit none
!
!    real(wp) :: y1, y2, y3, y4, phidloc, aalocal, HHlocal, backpert
!    complex(wp) :: y5, y6 
!    real(wp), intent(in) :: klocal, tphi
!    complex(wp), dimension(:) :: yin, yderiv
!
!! y1: Scale factor
!    y1 = aalocal
!
!! y2: Hubble rate
!    y2 = HHlocal
!
!! y5: the mode function 'q'
!    y5= yin(1)
!! y6: derivative of the mode function 'qdot'
!    y6= yin(2)
!   
!   phidloc = y4/y1**3
! 
!
!! qdot  
!    yderiv(1) = y6
!
!! qdot_dot 
!    yderiv(2) = -3.0_wp*(y2)*y6 - (klocal**2)*y5/y1**2
!    
!  end subroutine eqnpert



  real(wp) function effpotscalar (a, x, xd)
    use vars
    use kinds
   
    implicit none

    real(wp), intent(in) :: x, xd, a
    real(wp) :: rhotmp, rtmp, a2
  
    a2 = a**2
 
    rhotmp = xd**2/2.0_wp + potential(x)
    rtmp = 3.0_wp*8.0_wp*pi*xd**2/rhotmp 

    if (pert_type == 'scalar') then
    effpotscalar = potential(x)*rtmp - 2.0_wp*derpotential(x)*sqrt(rtmp) & 
                                       + der2potential(x)
    effpotscalar = effpotscalar*a2
    else if (pert_type == 'tensor') then
    effpotscalar = 0.0_wp 
    else 
    print*, 'wrong perturbation type! aborting the evolution...'
    stop
    end if

  end function effpotscalar


  real(wp) function renorm (klocal, at, at1, at2, at3)!, phit, phit1, phit2)
   use vars
   use kinds 

   implicit none

   real(wp), intent(in) :: klocal, at, at1, at2, at3!, phit, phit1, phit2

   renorm =   (0.5*klocal)/at**4 + &
               (0.25*at1**2)/(at**4*klocal) + &
              (0.015625*at**4*at1**4*klocal**2*mu**8)/&
               (klocal**2 + at**2*mu**2)**6.5 + &
              (0.015625*at**6*at1**4*mu**10)/&
               (klocal**2 + at**2*mu**2)**6.5 - &
              (0.03125*at**2*at1**4*klocal**2*mu**6)/&
               (klocal**2 + at**2*mu**2)**5.5 - &
              (0.03125*at**2*at1**2*(at1**2 + at*at2)*klocal**2*&
                 mu**6)/(klocal**2 + at**2*mu**2)**5.5 - &
              (0.45703125*at**4*at1**4*mu**8)/&
               (klocal**2 + at**2*mu**2)**5.5 - &
              (0.03125*at**4*at1**2*(at1**2 + at*at2)*mu**8)/&
               (klocal**2 + at**2*mu**2)**5.5 + &
              (0.015625*at1**4*klocal**2*mu**4)/&
               (klocal**2 + at**2*mu**2)**4.5 + &
              (0.03125*at1**2*(at1**2 + at*at2)*klocal**2*mu**4)/&
               (klocal**2 + at**2*mu**2)**4.5 + &
              (0.015625*(at1**2 + at*at2)**2*klocal**2*mu**4)/&
               (klocal**2 + at**2*mu**2)**4.5 - &
              (0.828125*at**2*at1**4*mu**6)/&
               (klocal**2 + at**2*mu**2)**4.5 + &
              (0.28125*at**2*at1**2*(at1**2 + at*at2)*mu**6)/&
               (klocal**2 + at**2*mu**2)**4.5 + &
              (0.015625*at**2*(at1**2 + at*at2)**2*mu**6)/&
               (klocal**2 + at**2*mu**2)**4.5 + &
              (0.25*at1**2*(at1**2 + at*at2)*klocal**4)/&
               (at**4*(klocal**2 + at**2*mu**2)**3.5) + &
              (0.0625*(at1**2 + at*at2)**2*klocal**4)/&
               (at**4*(klocal**2 + at**2*mu**2)**3.5) - &
              (0.125*at1*(at1**3 + 4.*at*at1*at2 + at**2*at3)*&
                 klocal**4)/(at**4*(klocal**2 + at**2*mu**2)**3.5) &
               + (0.9375*at1**2*(at1**2 + at*at2)*klocal**2*mu**2)/&
               (at**2*(klocal**2 + at**2*mu**2)**3.5) + &
              (0.125*(at1**2 + at*at2)**2*klocal**2*mu**2)/&
               (at**2*(klocal**2 + at**2*mu**2)**3.5) - &
              (0.3125*at1*(at1**3 + 4.*at*at1*at2 + at**2*at3)*&
                 klocal**2*mu**2)/&
               (at**2*(klocal**2 + at**2*mu**2)**3.5) + &
              (0.46875*at1**4*mu**4)/&
               (klocal**2 + at**2*mu**2)**3.5 + &
              (1.21875*at1**2*(at1**2 + at*at2)*mu**4)/&
               (klocal**2 + at**2*mu**2)**3.5 + &
              (0.0625*(at1**2 + at*at2)**2*mu**4)/&
               (klocal**2 + at**2*mu**2)**3.5 - &
              (0.21875*at1*(at1**3 + 4.*at*at1*at2 + at**2*at3)*&
                 mu**4)/(klocal**2 + at**2*mu**2)**3.5 + &
              (0.0625*at1**4*mu**2)/&
               (at**2*(klocal**2 + at**2*mu**2)**2.5) - &
              (0.0625*at1**2*(at1**2 + at*at2)*mu**2)/&
               (at**2*(klocal**2 + at**2*mu**2)**2.5) + &
              (0.0625*(at1**2 + at*at2)**2*mu**2)/&
               (at**2*(klocal**2 + at**2*mu**2)**2.5) - &
              (0.0625*at1*(at1**3 + 4.*at*at1*at2 + at**2*at3)*&
                 mu**2)/(at**2*(klocal**2 + at**2*mu**2)**2.5)

  end function renorm



  real(wp) function renormp (klocal, at, at1, at2, at3, at4)!, phit, phit1, phit2)
   use vars
   use kinds 

   implicit none

   real(wp), intent(in) :: klocal, at, at1, at2, at3, at4!, phit, phit1, phit2

   
   renormp =  (0.083333333333333333333_wp*at1**2)/(at**2*klocal) - &
            (0.16666666666666666667_wp*(at1**2 + at*at2))/(at**2*klocal) + &
            (0.16666666666666666667_wp*klocal)/at**4 - &
            (0.140625_wp*at**8*at1**4*klocal**4*mu**8)/&
             (klocal**2 + at**2*mu**2)**7.5_wp - &
            (0.359375_wp*at**10*at1**4*klocal**2*mu**10)/&
             (klocal**2 + at**2*mu**2)**7.5_wp - &
            (0.21875_wp*at**12*at1**4*mu**12)/(klocal**2 + at**2*mu**2)**7.5_wp + &
            (0.5625_wp*at**6*at1**4*klocal**4*mu**6)/&
             (klocal**2 + at**2*mu**2)**6.5_wp + &
            (0.28125_wp*at**6*at1**2*(at1**2 + at*at2)*klocal**4*mu**6)/&
             (klocal**2 + at**2*mu**2)**6.5_wp + &
            (1.4375_wp*at**8*at1**4*klocal**2*mu**8)/&
             (klocal**2 + at**2*mu**2)**6.5_wp + &
            (0.71875_wp*at**8*at1**2*(at1**2 + at*at2)*klocal**2*mu**8)/&
             (klocal**2 + at**2*mu**2)**6.5_wp - &
            (0.55078125_wp*at**10*at1**4*mu**10)/&
             (klocal**2 + at**2*mu**2)**6.5_wp + &
            (0.4375_wp*at**10*at1**2*(at1**2 + at*at2)*mu**10)/&
             (klocal**2 + at**2*mu**2)**6.5_wp - &
            (0.5625_wp*at**4*at1**4*klocal**4*mu**4)/&
             (klocal**2 + at**2*mu**2)**5.5_wp - &
            (0.5625_wp*at**4*at1**2*(at1**2 + at*at2)*klocal**4*mu**4)/&
             (klocal**2 + at**2*mu**2)**5.5_wp - &
            (0.140625_wp*at**4*(at1**2 + at*at2)**2*klocal**4*mu**4)/&
             (klocal**2 + at**2*mu**2)**5.5_wp - &
            (1.4375_wp*at**6*at1**4*klocal**2*mu**6)/&
             (klocal**2 + at**2*mu**2)**5.5_wp - &
            (1.4375_wp*at**6*at1**2*(at1**2 + at*at2)*klocal**2*mu**6)/&
             (klocal**2 + at**2*mu**2)**5.5_wp - &
            (0.359375_wp*at**6*(at1**2 + at*at2)**2*klocal**2*mu**6)/&
             (klocal**2 + at**2*mu**2)**5.5_wp - &
            (2.05859375_wp*at**8*at1**4*mu**8)/(klocal**2 + at**2*mu**2)**5.5_wp + &
            (0.171875_wp*at**8*at1**2*(at1**2 + at*at2)*mu**8)/&
             (klocal**2 + at**2*mu**2)**5.5_wp - &
            (0.21875*at**8*(at1**2 + at*at2)**2*mu**8)/&
             (klocal**2 + at**2*mu**2)**5.5_wp + &
            (0.0625_wp*at1**4*klocal**6)/(klocal**2 + at**2*mu**2)**4.5_wp - &
            (0.20833333333333333333_wp*at1**2*(at1**2 + at*at2)*klocal**6)/&
             (klocal**2 + at**2*mu**2)**4.5_wp + &
            (0.0625_wp*(at1**2 + at*at2)**2*klocal**6)/&
             (klocal**2 + at**2*mu**2)**4.5_wp + &
            (0.083333333333333333333_wp*at1*(at1**3 + 4*at*at1*at2 + at**2*at3)*&
               klocal**6)/(klocal**2 + at**2*mu**2)**4.5_wp + &
            (0.041666666666666666667_wp*&
               (at1**4 + 12*at*at1**2*at2 + 4*at**2*at2**2 + &
                 7*at**2*at1*at3 + at**3*at4)*klocal**6)/&
             (klocal**2 + at**2*mu**2)**4.5_wp + &
            (0.5_wp*at**2*at1**4*klocal**4*mu**2)/&
             (klocal**2 + at**2*mu**2)**4.5_wp - &
            (1.1041666666666666667_wp*at**2*at1**2*(at1**2 + at*at2)*klocal**4*&
               mu**2)/(klocal**2 + at**2*mu**2)**4.5_wp + &
            (0.1875_wp*at**2*(at1**2 + at*at2)**2*klocal**4*mu**2)/&
             (klocal**2 + at**2*mu**2)**4.5 + &
            (0.083333333333333333333_wp*at**2*at1*&
               (at1**3 + 4*at*at1*at2 + at**2*at3)*klocal**4*mu**2)/&
             (klocal**2 + at**2*mu**2)**4.5_wp + &
            (0.14583333333333333333_wp*at**2*&
               (at1**4 + 12*at*at1**2*at2 + 4*at**2*at2**2 + &
                 7*at**2*at1*at3 + at**3*at4)*klocal**4*mu**2)/&
             (klocal**2 + at**2*mu**2)**4.5_wp + &
            (1.3020833333333333333_wp*at**4*at1**4*klocal**2*mu**4)/&
             (klocal**2 + at**2*mu**2)**4.5_wp - &
            (1.3125_wp*at**4*at1**2*(at1**2 + at*at2)*klocal**2*mu**4)/&
             (klocal**2 + at**2*mu**2)**4.5_wp + &
            (0.1875_wp*at**4*(at1**2 + at*at2)**2*klocal**2*mu**4)/&
             (klocal**2 + at**2*mu**2)**4.5_wp - &
            (0.1875_wp*at**4*at1*(at1**3 + 4*at*at1*at2 + at**2*at3)*klocal**2*&
               mu**4)/(klocal**2 + at**2*mu**2)**4.5_wp + &
            (0.16666666666666666667_wp*at**4*&
               (at1**4 + 12*at*at1**2*at2 + 4*at**2*at2**2 + &
                 7*at**2*at1*at3 + at**3*at4)*klocal**2*mu**4)/&
             (klocal**2 + at**2*mu**2)**4.5_wp + &
            (4.1145833333333333333_wp*at**6*at1**4*mu**6)/&
             (klocal**2 + at**2*mu**2)**4.5_wp + &
            (1.1458333333333333333_wp*at**6*at1**2*(at1**2 + at*at2)*mu**6)/&
             (klocal**2 + at**2*mu**2)**4.5_wp + &
            (0.03125_wp*at**6*(at1**2 + at*at2)**2*mu**6)/&
             (klocal**2 + at**2*mu**2)**4.5_wp - &
            (0.33333333333333333333_wp*at**6*at1*&
               (at1**3 + 4*at*at1*at2 + at**2*at3)*mu**6)/&
             (klocal**2 + at**2*mu**2)**4.5_wp + &
            (0.0625_wp*at**6*(at1**4 + 12*at*at1**2*at2 + 4*at**2*at2**2 + &
                 7*at**2*at1*at3 + at**3*at4)*mu**6)/&
             (klocal**2 + at**2*mu**2)**4.5_wp - &
            (0.33333333333333333333_wp*at**4*at1**4*mu**4)/&
             (klocal**2 + at**2*mu**2)**3.5_wp - &
            (1.6979166666666666667_wp*at**4*at1**2*(at1**2 + at*at2)*mu**4)/&
             (klocal**2 + at**2*mu**2)**3.5_wp - &
            (0.09375_wp*at**4*(at1**2 + at*at2)**2*mu**4)/&
             (klocal**2 + at**2*mu**2)**3.5_wp - &
            (0.20833333333333333333_wp*at**4*at1*&
               (at1**3 + 4*at*at1*at2 + at**2*at3)*mu**4)/&
             (klocal**2 + at**2*mu**2)**3.5_wp + &
            (0.010416666666666666667_wp*at**4*&
               (at1**4 + 12*at*at1**2*at2 + 4*at**2*at2**2 + &
                 7*at**2*at1*at3 + at**3*at4)*mu**4)/&
             (klocal**2 + at**2*mu**2)**3.5_wp - &
            (0.25_wp*at**2*at1**4*mu**2)/(klocal**2 + at**2*mu**2)**2.5_wp - &
            (0.10416666666666666667_wp*at**2*at1**2*(at1**2 + at*at2)*mu**2)/&
             (klocal**2 + at**2*mu**2)**2.5_wp + &
            (0.16666666666666666667_wp*at**2*at1*&
               (at1**3 + 4*at*at1*at2 + at**2*at3)*mu**2)/&
             (klocal**2 + at**2*mu**2)**2.5_wp + &
            (0.020833333333333333333_wp*at**2*&
               (at1**4 + 12*at*at1**2*at2 + 4*at**2*at2**2 + &
                 7*at**2*at1*at3 + at**3*at4)*mu**2)/&
             (klocal**2 + at**2*mu**2)**2.5_wp  

!   renormp = (0.25*at1**2)/(at**4*klocal) - &
!             (0.16666666666666666*(at1**2 + at*at2))/(at**4*klocal) + &
!             (0.16666666666666666*klocal)/at**4 - &
!             (0.140625*at**4*at1**4*klocal**4*mu**8)/&
!              (klocal**2 + at**2*mu**2)**7.5 - &
!             (0.359375*at**6*at1**4*klocal**2*mu**10)/&
!              (klocal**2 + at**2*mu**2)**7.5 - &
!             (0.21875*at**8*at1**4*mu**12)/&
!              (klocal**2 + at**2*mu**2)**7.5 + &
!             (0.28125*at**2*at1**4*klocal**4*mu**6)/&
!              (klocal**2 + at**2*mu**2)**6.5 + &
!             (0.28125*at**2*at1**2*(at1**2 + at*at2)*klocal**4*mu**6)/&
!              (klocal**2 + at**2*mu**2)**6.5 + &
!             (0.71875*at**4*at1**4*klocal**2*mu**8)/&
!              (klocal**2 + at**2*mu**2)**6.5 + &
!             (0.71875*at**4*at1**2*(at1**2 + at*at2)*klocal**2*mu**8)/&
!              (klocal**2 + at**2*mu**2)**6.5 - &
!             (0.98828125*at**6*at1**4*mu**10)/&
!              (klocal**2 + at**2*mu**2)**6.5 + &
!             (0.4375*at**6*at1**2*(at1**2 + at*at2)*mu**10)/&
!              (klocal**2 + at**2*mu**2)**6.5 - &
!             (0.140625*at1**4*klocal**4*mu**4)/&
!              (klocal**2 + at**2*mu**2)**5.5 - &
!             (0.28125*at1**2*(at1**2 + at*at2)*klocal**4*mu**4)/&
!              (klocal**2 + at**2*mu**2)**5.5 - &
!             (0.140625*(at1**2 + at*at2)**2*klocal**4*mu**4)/&
!              (klocal**2 + at**2*mu**2)**5.5 - &
!             (0.359375*at**2*at1**4*klocal**2*mu**6)/&
!              (klocal**2 + at**2*mu**2)**5.5 - &
!             (0.71875*at**2*at1**2*(at1**2 + at*at2)*klocal**2*mu**6)/&
!              (klocal**2 + at**2*mu**2)**5.5 - &
!             (0.359375*at**2*(at1**2 + at*at2)**2*klocal**2*mu**6)/&
!              (klocal**2 + at**2*mu**2)**5.5 - &
!             (2.44921875*at**4*at1**4*mu**8)/&
!              (klocal**2 + at**2*mu**2)**5.5 + &
!             (0.609375*at**4*at1**2*(at1**2 + at*at2)*mu**8)/&
!              (klocal**2 + at**2*mu**2)**5.5 - &
!             (0.21875*at**4*(at1**2 + at*at2)**2*mu**8)/&
!              (klocal**2 + at**2*mu**2)**5.5 + &
!             (0.3333333333333333*at1**2*(at1**2 + at*at2)*klocal**6)/&
!              (at**4*(klocal**2 + at**2*mu**2)**4.5) - &
!             (0.10416666666666667*(at1**2 + at*at2)**2*klocal**6)/&
!              (at**4*(klocal**2 + at**2*mu**2)**4.5) - &
!             (0.20833333333333334*at1*&
!                (at1**3 + 4.*at*at1*at2 + at**2*at3)*klocal**6)/&
!              (at**4*(klocal**2 + at**2*mu**2)**4.5) + &
!             (0.041666666666666664*&
!                (at1**4 + 12.*at*at1**2*at2 + 4.*at**2*at2**2 + &
!                  7.*at**2*at1*at3 + at**3*at4)*klocal**6)/&
!              (at**4*(klocal**2 + at**2*mu**2)**4.5) + &
!             (1.6875*at1**2*(at1**2 + at*at2)*klocal**4*mu**2)/&
!              (at**2*(klocal**2 + at**2*mu**2)**4.5) - &
!             (0.3958333333333333*(at1**2 + at*at2)**2*klocal**4*mu**2)/&
!              (at**2*(klocal**2 + at**2*mu**2)**4.5) - &
!             (0.9375*at1*(at1**3 + 4.*at*at1*at2 + at**2*at3)*&
!                klocal**4*mu**2)/(at**2*(klocal**2 + at**2*mu**2)**4.5)&
!               + (0.14583333333333334*&
!                (at1**4 + 12.*at*at1**2*at2 + 4.*at**2*at2**2 + &
!                  7.*at**2*at1*at3 + at**3*at4)*klocal**4*mu**2)/&
!              (at**2*(klocal**2 + at**2*mu**2)**4.5) - &
!             (0.09375*at1**4*klocal**2*mu**4)/&
!              (klocal**2 + at**2*mu**2)**4.5 + &
!             (3.0625*at1**2*(at1**2 + at*at2)*klocal**2*mu**4)/&
!              (klocal**2 + at**2*mu**2)**4.5 - &
!             (0.4791666666666667*(at1**2 + at*at2)**2*klocal**2*mu**4)/&
!              (klocal**2 + at**2*mu**2)**4.5 - &
!             (1.3541666666666667*at1*&
!                (at1**3 + 4.*at*at1*at2 + at**2*at3)*klocal**2*mu**4)/&
!              (klocal**2 + at**2*mu**2)**4.5 + &
!             (0.16666666666666666*&
!                (at1**4 + 12.*at*at1**2*at2 + 4.*at**2*at2**2 + &
!                  7.*at**2*at1*at3 + at**3*at4)*klocal**2*mu**4)/&
!              (klocal**2 + at**2*mu**2)**4.5 + &
!             (1.125*at**2*at1**4*mu**6)/&
!              (klocal**2 + at**2*mu**2)**4.5 + &
!             (3.9166666666666665*at**2*at1**2*(at1**2 + at*at2)*mu**6)/&
!              (klocal**2 + at**2*mu**2)**4.5 - &
!             (0.21875*at**2*(at1**2 + at*at2)**2*mu**6)/&
!              (klocal**2 + at**2*mu**2)**4.5 - &
!             (0.7708333333333334*at**2*at1*&
!                (at1**3 + 4.*at*at1*at2 + at**2*at3)*mu**6)/&
!              (klocal**2 + at**2*mu**2)**4.5 + &
!             (0.0625*at**2*(at1**4 + 12.*at*at1**2*at2 + &
!                  4.*at**2*at2**2 + 7.*at**2*at1*at3 + at**3*at4)*mu**6&
!                )/(klocal**2 + at**2*mu**2)**4.5 + &
!             (0.5*at1**4*mu**4)/(klocal**2 + at**2*mu**2)**3.5 - &
!             (0.4270833333333333*at1**2*(at1**2 + at*at2)*mu**4)/&
!              (klocal**2 + at**2*mu**2)**3.5 - &
!             (0.13541666666666666*(at1**2 + at*at2)**2*mu**4)/&
!              (klocal**2 + at**2*mu**2)**3.5 - &
!             (0.28125*at1*(at1**3 + 4.*at*at1*at2 + at**2*at3)*mu**4)/&
!              (klocal**2 + at**2*mu**2)**3.5 + &
!             (0.010416666666666666*&
!                (at1**4 + 12.*at*at1**2*at2 + 4.*at**2*at2**2 + &
!                  7.*at**2*at1*at3 + at**3*at4)*mu**4)/&
!              (klocal**2 + at**2*mu**2)**3.5 + &
!             (0.0625*at1**4*mu**2)/&
!              (at**2*(klocal**2 + at**2*mu**2)**2.5) - &
!             (0.2708333333333333*at1**2*(at1**2 + at*at2)*mu**2)/&
!              (at**2*(klocal**2 + at**2*mu**2)**2.5) - &
!             (0.08333333333333333*(at1**2 + at*at2)**2*mu**2)/&
!              (at**2*(klocal**2 + at**2*mu**2)**2.5) + &
!             (0.020833333333333332*at1*&
!                (at1**3 + 4.*at*at1*at2 + at**2*at3)*mu**2)/&
!              (at**2*(klocal**2 + at**2*mu**2)**2.5) + &
!             (0.020833333333333332*&
!                (at1**4 + 12.*at*at1**2*at2 + 4.*at**2*at2**2 + &
!                  7.*at**2*at1*at3 + at**3*at4)*mu**2)/&
!              (at**2*(klocal**2 + at**2*mu**2)**2.5)

  end function renormp

end module rhs
