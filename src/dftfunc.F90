! Copyright 2014-2017  Kazuya Ishimura
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!-----------------------------------------------------------------------------------------
  subroutine funcsdlda(excora,excorb,energy,rhoa,rhob,rhoa13,rhob13,weight,csdlda,itype)
!-----------------------------------------------------------------------------------------
!
! Calculate Slater-Dirac LDA exchange functional
!
      implicit none
      integer,intent(in) :: itype
      real(8),parameter :: facp=-1.240700981798800D+00  !-(6/pi)**(1/3)
      real(8),parameter :: face1=-1.861051472698200D+00 !-((3/2)*(3/(4*pi))**(1/3))*2
      real(8),parameter :: face2=-9.305257363491000D-01 !-((3/4)*(3/(4*pi))**(1/3))*2
      real(8),intent(in) :: rhoa, rhob, rhoa13, rhob13, weight, csdlda
      real(8),intent(inout) :: excora(4), excorb(4), energy
      real(8) :: rhoa43, rhob43
!
      if(itype == 1) then
        rhoa43= rhoa13*rhoa
        energy= energy+csdlda*face1*rhoa43*weight
        excora(1)= excora(1)+csdlda*facp*rhoa13
      elseif(itype == 2) then
        rhoa43= rhoa13*rhoa
        rhob43= rhob13*rhob
        energy= energy+csdlda*face2*(rhoa43+rhob43)*weight
        excora(1)= excora(1)+csdlda*facp*rhoa13
        excorb(1)= excorb(1)+csdlda*facp*rhob13
      endif
      return
end


!-------------------------------------------------------------------------------------
  subroutine funcbecke88(excora,excorb,energy,rhoa,rhob,grhoa,grhob,rhoa13,rhob13, &
&                        weight,cb88,itype)
!-------------------------------------------------------------------------------------
!
! Calculate Becke88 exchange functional
!
      implicit none
      integer,intent(in) :: itype
      real(8),parameter :: small=1.0D-15
      real(8),parameter :: prmb=0.0042D+00
      real(8),parameter :: face=-0.9305257363491000D+00 !-(3/2)*(3/(4*pi))**(1/3)
      real(8),parameter :: one=1.0D+00, two=2.0D+00, six=6.0D+00, fourthird=1.333333333333333D+00
      real(8),intent(in) :: rhoa, rhob, grhoa(3), grhob(3), rhoa13, rhob13, weight, cb88
      real(8),intent(inout) :: excora(4), excorb(4), energy
      real(8) :: gradaa, gradbb, rhoa43, rhob43, xa, xb, sqrt1xa, sqrt1xb, sinhinva
      real(8) :: sinhinvb, denoma, denomb, gxa, gxb, gxap, gxbp, tmp
!
      if(itype == 1) then
        gradaa= sqrt(grhoa(1)*grhoa(1)+grhoa(2)*grhoa(2)+grhoa(3)*grhoa(3))
        rhoa43= rhoa13*rhoa
        xa= gradaa/rhoa43
!
        sqrt1xa= sqrt(one+xa*xa)
        sinhinva= log(xa+sqrt1xa)
        denoma= one/(one+six*prmb*xa*sinhinva)
        gxa= face-prmb*xa*xa*denoma
        gxap=(six*prmb*prmb*xa*xa*(xa/sqrt1xa-sinhinva)-two*prmb*xa)*(denoma*denoma)
!
        tmp= cb88*gxap/gradaa
        energy= energy+cb88*gxa*rhoa43*weight*two
        excora(1)= excora(1)+cb88*fourthird*rhoa13*(gxa-xa*gxap)
        excora(2)= excora(2)+tmp*grhoa(1)
        excora(3)= excora(3)+tmp*grhoa(2)
        excora(4)= excora(4)+tmp*grhoa(3)
      elseif(itype == 2) then
        gradaa= sqrt(grhoa(1)*grhoa(1)+grhoa(2)*grhoa(2)+grhoa(3)*grhoa(3))
        rhoa43= rhoa13*rhoa
        xa= gradaa/rhoa43
        sqrt1xa= sqrt(one+xa*xa)
        sinhinva= log(xa+sqrt1xa)
        denoma= one/(one+six*prmb*xa*sinhinva)
        gxa= face-prmb*xa*xa*denoma
        gxap=(six*prmb*prmb*xa*xa*(xa/sqrt1xa-sinhinva)-two*prmb*xa)*(denoma*denoma)
!
        tmp= cb88*gxap/gradaa
        energy= energy+cb88*gxa*rhoa43*weight
        excora(1)= excora(1)+cb88*fourthird*rhoa13*(gxa-xa*gxap)
        excora(2)= excora(2)+tmp*grhoa(1)
        excora(3)= excora(3)+tmp*grhoa(2)
        excora(4)= excora(4)+tmp*grhoa(3)
!
        if(rhob >= small) then
          gradbb= sqrt(grhob(1)*grhob(1)+grhob(2)*grhob(2)+grhob(3)*grhob(3))
          rhob43= rhob13*rhob
          xb= gradbb/rhob43
          sqrt1xb= sqrt(one+xb*xb)
          sinhinvb= log(xb+sqrt1xb)
          denomb= one/(one+six*prmb*xb*sinhinvb)
          gxb= face-prmb*xb*xb*denomb
          gxbp=(six*prmb*prmb*xb*xb*(xb/sqrt1xb-sinhinvb)-two*prmb*xb)*(denomb*denomb)
!  
          tmp= cb88*gxbp/gradbb
          energy= energy+cb88*gxb*rhob43*weight
          excorb(1)= excorb(1)+cb88*fourthird*rhob13*(gxb-xb*gxbp)
          excorb(2)= excorb(2)+tmp*grhob(1)
          excorb(3)= excorb(3)+tmp*grhob(2)
          excorb(4)= excorb(4)+tmp*grhob(3)
        endif
      endif
      return
end


!------------------------------------------------------------------------
  subroutine funcvwn5(excora,excorb,energy,rhoa,rhob,weight,cvwn,itype)
!------------------------------------------------------------------------
!
! Calculate VWN (Formula V) correlation functional
!
      implicit none
      integer,intent(in):: itype
      real(8),parameter :: facx=0.7876233178997432D+00 ! (3/(4*pi))**(1/6)
      real(8),parameter :: fach=1.7099209341613656D+00 ! 4/(9*(2**(1/3)-1))
      real(8),parameter :: prmqp=6.151990819759080D+00 ! Qp=sqrt(4*cp-bp*bp)
      real(8),parameter :: prmqf=4.730926909560112D+00 ! Qf=sqrt(4*cf-bf*bf)
      real(8),parameter :: prmqa=7.123108917818117D+00 ! Qa=sqrt(4*ca-ba*ba)
      real(8),parameter :: prmqpinv=0.1625490071910025D+00 ! 1/Qp
      real(8),parameter :: prmqfinv=0.2113750685894619D+00 ! 1/Qf
      real(8),parameter :: prmqainv=0.1403881383167604D+00 ! 1/Qa
      real(8),parameter :: prmxx0pinv=7.965008666058621D-02 ! 1/Xx0p=1/(x0p*x0p+bp*x0p+cp)
      real(8),parameter :: prmxx0finv=6.301678291320096D-02 ! 1/Xx0f=1/(x0f*x0f+bf*x0f+cf)
      real(8),parameter :: prmxx0ainv=7.692816270970167D-02 ! 1/Xx0a=1/(x0a*x0a+ba*x0a+ca)
      real(8),parameter :: prmap=3.10907D-02, prmbp=3.72744D+00, prmcp=1.29352D+01
      real(8),parameter :: prmaf=1.554535D-02, prmbf=7.06042D+00, prmcf=1.80578D+01
      real(8),parameter :: prmaa=-1.688686D-02, prmba=1.13107D+00, prmca=1.30045D+01
      real(8),parameter :: prmx0p=-0.10498D+00, prmx0f=-0.325D+00, prmx0a=-0.0047584D+00
      real(8),parameter :: prmp2=9.690227711544374D-04  ! -ap*bp*x0p/Xx0p
      real(8),parameter :: prmf2=2.247867095542611D-03  ! -af*bf*x0f/Xx0f
      real(8),parameter :: prma2=-6.991730719309994D-06 ! -aa*ba*x0a/Xx0a
      real(8),parameter :: prmp3=3.878329487811301D-02  ! ap*2*bp/Qp*(cp-x0p*x0p)/Xx0p
      real(8),parameter :: prmf3=5.249139316978093D-02  ! af*2*bf/Qf*(cf-x0f*x0f)/Xx0f
      real(8),parameter :: prma3=-5.365090596990275D-03 ! aa*2*ba/Qa*(ca-x0a*x0a)/Xx0a
      real(8),parameter :: one=1.0D+00, two=2.0D+00, onethird=0.3333333333333333D+00
      real(8),parameter :: onesixth=0.1666666666666666D+00, four=4.0D+00
      real(8),parameter :: fourthird=1.333333333333333D+00, nineeighth=1.125D+00
      real(8),parameter :: onehalf=1.5D+00
      real(8),intent(in) :: rhoa, rhob, weight, cvwn
      real(8),intent(inout) :: excora(4), excorb(4), energy
      real(8) :: rhot, x, xxp, xxf, xxa, xxpinv, xxfinv, xxainv, tmp1p, tmp2p, tmp3p
      real(8) :: tmp1f, tmp2f, tmp3f, tmp1a, tmp2a, tmp3a, epsp, epsf, epsa, xbq2p, xbq2f, xbq2a
      real(8) :: depsp, depsf, depsa, xbp, xbf, xba, bracketp, bracketf, bracketa, hx, dhx
      real(8) :: zeta, zeta2, zeta3, zeta4, gzeta, dgzeta, vwnpot, dvwnpot1, dvwnpot2
!
      if(itype == 1) then
        rhot= rhoa*two
        x= facx*(rhot**(-onesixth))
        xxp= x*x+prmbp*x+prmcp
        xxpinv=one/xxp
        xbp= two*x+prmbp
!  
        tmp1p= log(x*x*xxpinv)
        tmp2p= log((x-prmx0p)*(x-prmx0p)*xxpinv)
        tmp3p= atan(prmqp/xbp)
        epsp= prmap*tmp1p+prmp2*tmp2p+prmp3*tmp3p
        depsp= onethird*prmap*(one/x-(x*xxpinv)*(one+prmbp/(x-prmx0p)))
        energy= energy+cvwn*rhot*epsp*weight
        excora(1)= excora(1)+cvwn*(epsp-x*depsp)
      elseif(itype == 2) then
        rhot= rhoa+rhob
        x= facx*(rhot**(-onesixth))
        xxp= x*x+prmbp*x+prmcp
        xxf= x*x+prmbf*x+prmcf
        xxa= x*x+prmba*x+prmca
        xxpinv=one/xxp
        xxfinv=one/xxf
        xxainv=one/xxa
        xbp= two*x+prmbp
        xbf= two*x+prmbf
        xba= two*x+prmba
!
        tmp1p= log(x*x*xxpinv)
        tmp1f= log(x*x*xxfinv)
        tmp1a= log(x*x*xxainv)
        tmp2p= log((x-prmx0p)*(x-prmx0p)*xxpinv)
        tmp2f= log((x-prmx0f)*(x-prmx0f)*xxfinv)
        tmp2a= log((x-prmx0a)*(x-prmx0a)*xxainv)
        tmp3p= atan(prmqp/xbp)
        tmp3f= atan(prmqf/xbf)
        tmp3a= atan(prmqa/xba)
        epsp= prmap*tmp1p+prmp2*tmp2p+prmp3*tmp3p
        epsf= prmaf*tmp1f+prmf2*tmp2f+prmf3*tmp3f
        epsa= prmaa*tmp1a+prma2*tmp2a+prma3*tmp3a
!
        xbq2p= one/(xbp*xbp+prmqp*prmqp)
        xbq2f= one/(xbf*xbf+prmqf*prmqf)
        xbq2a= one/(xba*xba+prmqa*prmqa)
        bracketp= two/(x-prmx0p)-xbp*xxpinv-four*(two*prmx0p+prmbp)*xbq2p
        bracketf= two/(x-prmx0f)-xbf*xxfinv-four*(two*prmx0f+prmbf)*xbq2f
        bracketa= two/(x-prmx0a)-xba*xxainv-four*(two*prmx0a+prmba)*xbq2a
        depsp=prmap*(two/x-xbp*xxpinv-four*prmbp*xbq2p-prmbp*prmx0p*prmxx0pinv*bracketp)
        depsf=prmaf*(two/x-xbf*xxfinv-four*prmbf*xbq2f-prmbf*prmx0f*prmxx0finv*bracketf)
        depsa=prmaa*(two/x-xba*xxainv-four*prmba*xbq2a-prmba*prmx0a*prmxx0ainv*bracketa)
!
        hx= fach*(epsf-epsp)/epsa-one
        dhx= fach*(depsf-depsp-(epsf-epsp)*depsa/epsa)/epsa
!
        zeta=(rhoa-rhob)/rhot
        zeta2= zeta *zeta
        zeta3= zeta2*zeta
        zeta4= zeta2*zeta2
        gzeta= nineeighth*((one+zeta)**fourthird+(one-zeta)**fourthird-two)
        dgzeta= onehalf*((one+zeta)**onethird-(one-zeta)**onethird)
        vwnpot= epsp+epsa*gzeta*(one+hx*zeta4)
        energy= energy+cvwn*rhot*vwnpot*weight
        dvwnpot1=-x*onesixth*(depsp+depsa*gzeta*(one+hx*zeta4)+epsa*gzeta*dhx*zeta4)
        dvwnpot2= epsa*(dgzeta*(one+hx*zeta4)+four*gzeta*hx*zeta3)
        excora(1)= excora(1)+cvwn*(vwnpot+dvwnpot1+dvwnpot2*(one-zeta))
        excorb(1)= excorb(1)+cvwn*(vwnpot+dvwnpot1-dvwnpot2*(one+zeta))
      endif
      return
end


!------------------------------------------------------------------------
  subroutine funcvwn3(excora,excorb,energy,rhoa,rhob,weight,cvwn,itype)
!------------------------------------------------------------------------
!
! Calculate VWN (Formula III) correlation functional
!
      implicit none
      integer,intent(in):: itype
      real(8),parameter :: facx=0.7876233178997432D+00 ! (3/(4*pi))**(1/6)
      real(8),parameter :: fach=1.7099209341613656D+00 ! 4/(9*(2**(1/3)-1))
      real(8),parameter :: prmqp= 4.489988864157680D-02 ! Qp=sqrt(4*cp-bp*bp)
      real(8),parameter :: prmqf= 1.171685277708971D+00 ! Qf=sqrt(4*cf-bf*bf)
      real(8),parameter :: prmqa= 6.692072046645942D+00 ! Qa=sqrt(4*ca-ba*ba)
      real(8),parameter :: prmqpinv= 2.227177015922510D+01 ! 1/Qp
      real(8),parameter :: prmqfinv= 8.534715072594644D-01 ! 1/Qf
      real(8),parameter :: prmqainv= 1.494305490182520D-01 ! 1/Qa
      real(8),parameter :: prmxx0pinv= 2.664029033699615D-02 ! 1/Xx0p=1/(x0p*x0p+bp*x0p+cp)
      real(8),parameter :: prmxx0finv= 2.664029033699615D-02 ! 1/Xx0f=1/(x0f*x0f+bf*x0f+cf)
      real(8),parameter :: prmxx0ainv= 2.664029033699615D-02 ! 1/Xx0a=1/(x0a*x0a+ba*x0a+ca)
      real(8),parameter :: prmap=3.10907D-02, prmbp=1.30720D+01, prmcp=4.27198D+01
      real(8),parameter :: prmaf=1.554535D-02, prmbf=2.01231D+01, prmcf=1.01578D+02
      real(8),parameter :: prmaa=-1.688686D-02, prmba=1.06835D+00, prmca=1.14813D+01
      real(8),parameter :: prmx0p=-4.09286D-01, prmx0f=-7.43294D-01, prmx0a=-2.28344D-01
      real(8),parameter :: prmp2= 4.431373767749538D-03 !-ap*bp*x0p/Xx0p
      real(8),parameter :: prmf2= 2.667310007273315D-03 !-af*bf*x0f/Xx0f
      real(8),parameter :: prma2=-3.649032666450384D-04 !-aa*ba*x0a/Xx0a
      real(8),parameter :: prmp3= 2.052197293770518D+01 ! ap*2*bp/Qp*(cp-x0p*x0p)/Xx0p
      real(8),parameter :: prmf3= 6.188180297906176D-01 ! af*2*bf/Qf*(cf-x0f*x0f)/Xx0f
      real(8),parameter :: prma3=-5.458481084953849D-03 ! aa*2*ba/Qa*(ca-x0a*x0a)/Xx0a
      real(8),parameter :: one=1.0D+00, two=2.0D+00, onethird=0.3333333333333333D+00
      real(8),parameter :: onesixth=0.1666666666666666D+00, four=4.0D+00
      real(8),parameter :: fourthird=1.333333333333333D+00, nineeighth=1.125D+00
      real(8),parameter :: onehalf=1.5D+00
      real(8),intent(in) :: rhoa, rhob, weight, cvwn
      real(8),intent(inout) :: excora(4), excorb(4), energy
      real(8) :: rhot, x, xxp, xxf, xxa, xxpinv, xxfinv, xxainv, tmp1p, tmp2p, tmp3p
      real(8) :: tmp1f, tmp2f, tmp3f, tmp1a, tmp2a, tmp3a, epsp, epsf, epsa, xbq2p, xbq2f, xbq2a
      real(8) :: depsp, depsf, depsa, xbp, xbf, xba, bracketp, bracketf, bracketa, hx, dhx
      real(8) :: zeta, zeta2, zeta3, zeta4, gzeta, dgzeta, vwnpot, dvwnpot1, dvwnpot2
!
      if(itype == 1) then
        rhot= rhoa*two
        x= facx*(rhot**(-onesixth))
        xxp= x*x+prmbp*x+prmcp
        xxpinv=one/xxp
        xbp= two*x+prmbp
!  
        tmp1p= log(x*x*xxpinv)
        tmp2p= log((x-prmx0p)*(x-prmx0p)*xxpinv)
        tmp3p= atan(prmqp/xbp)
        epsp= prmap*tmp1p+prmp2*tmp2p+prmp3*tmp3p
        depsp= onethird*prmap*(one/x-(x*xxpinv)*(one+prmbp/(x-prmx0p)))
        energy= energy+cvwn*rhot*epsp*weight
        excora(1)= excora(1)+cvwn*(epsp-x*depsp)
      elseif(itype == 2) then
        rhot= rhoa+rhob
        x= facx*(rhot**(-onesixth))
        xxp= x*x+prmbp*x+prmcp
        xxf= x*x+prmbf*x+prmcf
        xxa= x*x+prmba*x+prmca
        xxpinv=one/xxp
        xxfinv=one/xxf
        xxainv=one/xxa
        xbp= two*x+prmbp
        xbf= two*x+prmbf
        xba= two*x+prmba
!
        tmp1p= log(x*x*xxpinv)
        tmp1f= log(x*x*xxfinv)
        tmp1a= log(x*x*xxainv)
        tmp2p= log((x-prmx0p)*(x-prmx0p)*xxpinv)
        tmp2f= log((x-prmx0f)*(x-prmx0f)*xxfinv)
        tmp2a= log((x-prmx0a)*(x-prmx0a)*xxainv)
        tmp3p= atan(prmqp/xbp)
        tmp3f= atan(prmqf/xbf)
        tmp3a= atan(prmqa/xba)
        epsp= prmap*tmp1p+prmp2*tmp2p+prmp3*tmp3p
        epsf= prmaf*tmp1f+prmf2*tmp2f+prmf3*tmp3f
        epsa= prmaa*tmp1a+prma2*tmp2a+prma3*tmp3a
!
        xbq2p= one/(xbp*xbp+prmqp*prmqp)
        xbq2f= one/(xbf*xbf+prmqf*prmqf)
        xbq2a= one/(xba*xba+prmqa*prmqa)
        bracketp= two/(x-prmx0p)-xbp*xxpinv-four*(two*prmx0p+prmbp)*xbq2p
        bracketf= two/(x-prmx0f)-xbf*xxfinv-four*(two*prmx0f+prmbf)*xbq2f
        bracketa= two/(x-prmx0a)-xba*xxainv-four*(two*prmx0a+prmba)*xbq2a
        depsp=prmap*(two/x-xbp*xxpinv-four*prmbp*xbq2p-prmbp*prmx0p*prmxx0pinv*bracketp)
        depsf=prmaf*(two/x-xbf*xxfinv-four*prmbf*xbq2f-prmbf*prmx0f*prmxx0finv*bracketf)
        depsa=prmaa*(two/x-xba*xxainv-four*prmba*xbq2a-prmba*prmx0a*prmxx0ainv*bracketa)
!
        hx= fach*(epsf-epsp)/epsa-one
        dhx= fach*(depsf-depsp-(epsf-epsp)*depsa/epsa)/epsa
!
        zeta=(rhoa-rhob)/rhot
        zeta2= zeta *zeta
        zeta3= zeta2*zeta
        zeta4= zeta2*zeta2
        gzeta= nineeighth*((one+zeta)**fourthird+(one-zeta)**fourthird-two)
        dgzeta= onehalf*((one+zeta)**onethird-(one-zeta)**onethird)
        vwnpot= epsp+epsa*gzeta*(one+hx*zeta4)
        energy= energy+cvwn*rhot*vwnpot*weight
        dvwnpot1=-x*onesixth*(depsp+depsa*gzeta*(one+hx*zeta4)+epsa*gzeta*dhx*zeta4)
        dvwnpot2= epsa*(dgzeta*(one+hx*zeta4)+four*gzeta*hx*zeta3)
        excora(1)= excora(1)+cvwn*(vwnpot+dvwnpot1+dvwnpot2*(one-zeta))
        excorb(1)= excorb(1)+cvwn*(vwnpot+dvwnpot1-dvwnpot2*(one+zeta))
      endif
      return
end


!-------------------------------------------------------------------------------------------------
  subroutine funclyp(excora,excorb,energy,rhoa,rhob,grhoa,grhob,rhoa13,rhob13,weight,clyp,itype)
!-------------------------------------------------------------------------------------------------
!
! Calculate LYP correlation functional
!
      implicit none
      integer,intent(in) :: itype
      real(8),parameter :: face=3.646239897876478D+01 ! 2**(11/3)*3/10*(3*pi*pi)**(2/3)
      real(8),parameter :: prma=0.04918D+00, prmb=0.132D+00, prmc=0.2533D+00, prmd=0.349D+00
      real(8),parameter :: half=0.5D+00, one=1.0D+00, two=2.0D+00, three=3.0D+00, four=4.0D+00
      real(8),parameter :: five=5.0D+00, seven=7.0D+00, eleven=11.0D+00, p47=47.0D+00
      real(8),parameter :: onethird=0.3333333333333333D+00, oneninth=0.1111111111111111D+00
      real(8),parameter :: fourthird=1.333333333333333D+00, eightthird=2.6666666666666666D+00
      real(8),parameter :: p11third=3.666666666666666D+00
      real(8),intent(in) :: rhoa, rhob, grhoa(3), grhob(3), rhoa13, rhob13, weight, clyp
      real(8),intent(inout) :: excora(4), excorb(4), energy
      real(8) :: rhot, rhoab, rhot13, rhotinv, rhot13inv, rhot43inv, rhot53inv, rhoa83, rhob83
      real(8) :: rhoatinv, rhobtinv, denom, omega, delta, delta11a, delta11b, omegap, deltap
      real(8) :: gradaa, gradab, gradbb, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
      real(8) :: tmp10, tmp11, tmp12
      real(8) :: dlypaa, dlypab, dlypbb, dlypaaa, dlypaab, dlypabb, dlypbaa, dlypbab, dlypbbb
!
      rhot= rhoa+rhob
      rhoab= rhoa*rhob
      rhot13= rhot**onethird
      rhotinv= one/rhot
      rhot13inv= rhot13*rhot13*rhotinv
      rhot43inv= rhotinv*rhot13inv
      rhot53inv= rhotinv*rhot13inv*rhot13inv
      rhoa83= rhoa*rhoa*rhoa13*rhoa13
      rhob83= rhob*rhob*rhob13*rhob13
      rhoatinv= rhoa*rhotinv
      rhobtinv= rhob*rhotinv
      denom= one/(one+prmd*rhot13inv)
      omega=exp(-prmc*rhot13inv)*denom*rhot53inv*rhot53inv*rhot13inv*prma*prmb
      delta=prmc*rhot13inv+prmd*rhot13inv*denom
      delta11a=(delta-eleven)*rhoatinv
      delta11b=(delta-eleven)*rhobtinv
      gradaa= grhoa(1)*grhoa(1)+grhoa(2)*grhoa(2)+grhoa(3)*grhoa(3)
      gradab= grhoa(1)*grhob(1)+grhoa(2)*grhob(2)+grhoa(3)*grhob(3)
      gradbb= grhob(1)*grhob(1)+grhob(2)*grhob(2)+grhob(3)*grhob(3)
!
      tmp1= oneninth*rhoab
      tmp2= four*prma*denom*rhob*rhoatinv
      dlypaa=-omega*(tmp1*(one-three*delta-delta11a)-rhob*rhob)
      dlypab=-omega*(tmp1*(p47-seven*delta)-fourthird*rhot*rhot)
      dlypbb=-omega*(tmp1*(one-three*delta-delta11b)-rhoa*rhoa)
      tmp3= dlypaa*gradaa+dlypab*gradab+dlypbb*gradbb
      energy= energy+clyp*(-tmp2-face*omega*rhoab*(rhoa83+rhob83)+tmp3)*weight
!
      omegap=-onethird*rhot43inv*(eleven*rhot13-prmc-prmd*denom)
      deltap=onethird*(prmd*prmd*rhot53inv*denom*denom-rhotinv*delta)*oneninth*rhoab
      tmp4= oneninth*rhob*(one-three*delta-delta11a)-(three+rhoatinv)*deltap
      tmp5= oneninth*rhob*(one-three*delta-delta11b)-(three+rhobtinv)*deltap
      dlypaaa= omegap*dlypaa-omega*(tmp4-oneninth*rhoab*delta11b*rhotinv)
      dlypaab= omegap*dlypab-omega*(oneninth*rhob*(p47-seven*delta)-seven*deltap-eightthird*rhot)
      dlypabb= omegap*dlypbb-omega*(tmp5+oneninth*rhoab*delta11b*rhotinv-two*rhoa)
      tmp6=-four*prma*denom*rhobtinv*(onethird*prmd*rhot43inv*denom*rhoa+one-rhoatinv)
      tmp7=-face*omega*rhob*(omegap*rhoa*(rhoa83+rhob83)+(p11third*rhoa83+rhob83))
      tmp8= dlypaaa*gradaa+dlypaab*gradab+dlypabb*gradbb
      tmp9= clyp*(dlypaa+dlypab+dlypbb)
!
      if(itype == 1) then
        excora(1)= excora(1)+clyp*(tmp6+tmp7+tmp8)
        excora(2)= excora(2)+tmp9*grhoa(1)
        excora(3)= excora(3)+tmp9*grhoa(2)
        excora(4)= excora(4)+tmp9*grhoa(3)
      elseif(itype == 2) then
        tmp4= oneninth*rhoa*(one-three*delta-delta11a)-(three+rhoatinv)*deltap
        tmp5= oneninth*rhoa*(one-three*delta-delta11b)-(three+rhobtinv)*deltap
        dlypbaa= omegap*dlypaa-omega*(tmp4+oneninth*rhoab*delta11a*rhotinv-two*rhob)
        dlypbab= omegap*dlypab-omega*(oneninth*rhoa*(p47-seven*delta)-seven*deltap-eightthird*rhot)
        dlypbbb= omegap*dlypbb-omega*(tmp5-oneninth*rhoab*delta11a*rhotinv)
        tmp10=-four*prma*denom*rhoatinv*(onethird*prmd*rhot43inv*denom*rhob+one-rhobtinv)
        tmp11=-face*omega*rhoa*(omegap*rhob*(rhoa83+rhob83)+(rhoa83+p11third*rhob83))
        tmp12= dlypbaa*gradaa+dlypbab*gradab+dlypbbb*gradbb
        excora(1)= excora(1)+clyp*(tmp6+tmp7+tmp8)
        excora(2)= excora(2)+clyp*(two*dlypaa*grhoa(1)+dlypab*grhob(1))
        excora(3)= excora(3)+clyp*(two*dlypaa*grhoa(2)+dlypab*grhob(2))
        excora(4)= excora(4)+clyp*(two*dlypaa*grhoa(3)+dlypab*grhob(3))
        excorb(1)= excorb(1)+clyp*(tmp10+tmp11+tmp12)
        excorb(2)= excorb(2)+clyp*(two*dlypbb*grhob(1)+dlypab*grhoa(1))
        excorb(3)= excorb(3)+clyp*(two*dlypbb*grhob(2)+dlypab*grhoa(2))
        excorb(4)= excorb(4)+clyp*(two*dlypbb*grhob(3)+dlypab*grhoa(3))
      endif
      return
end
