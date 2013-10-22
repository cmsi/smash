!---------------------------------------------------------------------------
  subroutine funcsdlda(excor,energy,rhoa,rhob,rhoa13,rhob13,weight,csdlda)
!---------------------------------------------------------------------------
!
! Calculate Slater-Dirac LDA exchange functional
!
      implicit none
      real(8),parameter :: facp=-1.240700981798800D+00 !-(6/pi)**(1/3)
      real(8),parameter :: face=-1.861051472698200D+00 !-((3/2)*(3/(4*pi))**(1/3))*2
      real(8),intent(in) :: rhoa, rhob, rhoa13, rhob13, weight, csdlda
      real(8),intent(inout) :: excor(4), energy
      real(8) :: rhoa43
!
      rhoa43= rhoa13*rhoa
      energy= energy+csdlda*face*rhoa43*weight
      excor(1)= excor(1)+csdlda*facp*rhoa13
      return
end


!---------------------------------------------------------------------------------------
  subroutine funcbecke88(excor,energy,rhoa,rhob,grhoa,grhob,rhoa13,rhob13,weight,cb88)
!---------------------------------------------------------------------------------------
!
! Calculate Becke88 exchange functional
!
      implicit none
      real(8),parameter :: prmb=0.0042D+00
      real(8),parameter :: face=-0.9305257363491000D+00 !-(3/2)*(3/(4*pi))**(1/3)
      real(8),parameter :: one=1.0D+00, two=2.0D+00, six=6.0D+00, fourthird=1.333333333333333D+00
      real(8),intent(in) :: rhoa, rhob, grhoa(3), grhob(3), rhoa13, rhob13, weight, cb88
      real(8),intent(inout) :: excor(4), energy
      real(8) :: gradaa, rhoa43, xa, sqrt1xa, sinhinv, denom, gxa, gxap, tmp
!
      gradaa= sqrt(grhoa(1)*grhoa(1)+grhoa(2)*grhoa(2)+grhoa(3)*grhoa(3))
      rhoa43= rhoa13*rhoa
      xa= gradaa/rhoa43
!
      sqrt1xa= sqrt(one+xa*xa)
      sinhinv= log(xa+sqrt1xa)
      denom= one/(one+six*prmb*xa*sinhinv)
      gxa= face-prmb*xa*xa*denom
      gxap=(six*prmb*prmb*xa*xa*(xa/sqrt1xa-sinhinv)-two*prmb*xa)*(denom*denom)
!
      tmp= cb88*gxap/gradaa
      energy= energy+cb88*gxa*rhoa43*weight*two
      excor(1)= excor(1)+cb88*fourthird*rhoa13*(gxa-xa*gxap)
      excor(2)= excor(2)+tmp*grhoa(1)
      excor(3)= excor(3)+tmp*grhoa(2)
      excor(4)= excor(4)+tmp*grhoa(3)
      return
end


!---------------------------------------------------------------------------------------
  subroutine funcvwn5(excor,energy,rhoa,rhob,rhoa13,rhob13,weight,cvwn)
!---------------------------------------------------------------------------------------
!
! Calculate VWN (Formula V) correlation functional
!
      implicit none
      real(8),parameter :: facx=0.7016926042943222D+00 ! (3/(4*2*pi))**(1/6)  
      real(8),parameter :: prmqp=6.151990819759080D+00 ! Q=sqrt(4*c-b*b)
      real(8),parameter :: prmqpinv=0.1625490071910025D+00 ! 1/Q
      real(8),parameter :: prmxx0pinv=7.965008666058621D-02! 1/Xx0=1/(x0*x0+b*x0+c)
      real(8),parameter :: prmap=3.10907D-02, prmbp=3.72744D+00, prmcp=1.29352D+01
      real(8),parameter :: prmx0p=-0.10498D+00
      real(8),parameter :: prmp1=9.690227711544374D-04, prmp2=3.878329487811301D-02
      real(8),parameter :: one=1.0D+00, two=2.0D+00, onethird=0.3333333333333333D+00
      real(8),intent(in) :: rhoa, rhob, rhoa13, rhob13, weight, cvwn
      real(8),intent(inout) :: excor(4), energy
      real(8) :: rhot, x, xx, xxinv, tmp1, tmp2, tmp3, eps, deps
!
      rhot= rhoa+rhob
      x= facx/sqrt(rhoa13)
      xx= x*x+prmbp*x+prmcp
      xxinv=one/xx
!
      tmp1= log(x*x*xxinv)
      tmp2= log((x-prmx0p)*(x-prmx0p)*xxinv)
      tmp3= atan(prmqp/(two*x+prmbp))
      eps= prmap*tmp1+prmp1*tmp2+prmp2*tmp3
      deps= onethird*prmap*(one/x-(x*xxinv)*(one+prmbp/(x-prmx0p)))
      energy= energy+cvwn*rhot*eps*weight
      excor(1)= excor(1)+cvwn*(eps-x*deps)
      return
end


!---------------------------------------------------------------------------------------
  subroutine funclyp(excor,energy,rhoa,rhob,grhoa,grhob,rhoa13,rhob13,weight,clyp)
!---------------------------------------------------------------------------------------
!
! Calculate LYP correlation functional
!
      implicit none
      real(8),parameter :: face=3.646239897876478D+01 ! 2**(11/3)*3/10*(3*pi*pi)**(2/3)
      real(8),parameter :: fac13=1.259921049894873D+00 ! 2**(1/3)
      real(8),parameter :: prma=0.04918D+00, prmb=0.132D+00, prmc=0.2533D+00, prmd=0.349D+00
      real(8),parameter :: half=0.5D+00, one=1.0D+00, two=2.0D+00, three=3.0D+00
      real(8),parameter :: four=4.0D+00, seven=7.0D+00, eleven=11.0D+00, p47=47.0D+00
      real(8),parameter :: onethird=0.3333333333333333D+00, oneninth=0.1111111111111111D+00
      real(8),parameter :: fourthird=1.333333333333333D+00, eightthird=2.6666666666666666D+00
      real(8),parameter :: p14third=4.666666666666666D+00
      real(8),intent(in) :: rhoa, rhob, grhoa(3), grhob(3), rhoa13, rhob13, weight, clyp
      real(8),intent(inout) :: excor(4), energy
      real(8) :: rhot, rhoab, rhot13, rhotinv, rhot13inv, rhot43inv, rhot53inv, rhoa83, rhob83
      real(8) :: rhoatinv, denom, omega, delta, delta11, omegap, deltap
      real(8) :: gradaa, gradab, gradbb, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
      real(8) :: dlypaa, dlypab, dlypbb, dlypaaa, dlypaab, dlypabb, dlyp
!
      rhot= rhoa+rhob
      rhoab= rhoa*rhob
      rhot13= fac13*rhoa13
      rhotinv= one/rhot
      rhot13inv= rhot13*rhot13*rhotinv
      rhot43inv= rhotinv*rhot13inv
      rhot53inv= rhotinv*rhot13inv*rhot13inv
      rhoa83= rhoa*rhoa*rhoa13*rhoa13
      rhob83= rhoa83
      rhoatinv= rhoa*rhotinv
      denom= one/(one+prmd*rhot13inv)
      omega=exp(-prmc*rhot13inv)*denom*rhot53inv*rhot53inv*rhot13inv*prma*prmb
      delta=prmc*rhot13inv+prmd*rhot13inv*denom
      delta11=(delta-eleven)*rhoatinv
      gradaa= grhoa(1)*grhoa(1)+grhoa(2)*grhoa(2)+grhoa(3)*grhoa(3)
      gradab= gradaa
      gradbb= gradaa
!
      tmp1= oneninth*rhoab
      tmp2= four*prma*denom*rhob*rhoatinv
      dlypaa=-omega*(tmp1*(one-three*delta-delta11)-rhob*rhob)
      dlypab=-omega*(tmp1*(p47-seven*delta)-fourthird*rhot*rhot)
      dlypbb= dlypaa
      tmp3= dlypaa*gradaa+dlypab*gradab+dlypbb*gradbb
      energy= energy+(-tmp2-face*omega*rhoab*(rhoa83+rhob83)+tmp3)*weight*clyp
!
      omegap=-onethird*rhot43inv*(eleven*rhot13-prmc-prmd*denom)
      deltap=onethird*(prmd*prmd*rhot53inv*denom*denom-rhotinv*delta)*oneninth*rhoab
      tmp4= oneninth*rhob*(one-three*delta-delta11)-(three+rhoatinv)*deltap
      tmp5= oneninth*rhob*delta11*rhoatinv
      dlypaaa= omegap*dlypaa-omega*(tmp4-tmp5)
      dlypaab= omegap*dlypab-omega*(oneninth*rhob*(p47-seven*delta)-seven*deltap-eightthird*rhot)
      dlypabb= omegap*dlypbb-omega*(tmp4+tmp5-two*rhoa)
      tmp6=-four*prma*denom*rhoatinv*(onethird*prmd*rhot43inv*denom*rhoa+one-rhoatinv)
      tmp7=-face*omega*rhob*rhoa83*(omegap*rhoa*two+p14third)
      tmp8= dlypaaa*gradaa+dlypaab*gradab+dlypabb*gradbb
      tmp9= clyp*(two*dlypaa+dlypab)
      dlyp= tmp6+tmp7+tmp8
      excor(1)= excor(1)+clyp*dlyp
      excor(2)= excor(2)+tmp9*grhoa(1)
      excor(3)= excor(3)+tmp9*grhoa(2)
      excor(4)= excor(4)+tmp9*grhoa(3)
      return
end
