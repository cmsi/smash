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
!--------------------------------------------------------------------------
  subroutine int2dddd(phmdint,exfac1,exfac2,xyziq,xzkl,rot,nijkl,nbfijkl)
!--------------------------------------------------------------------------
!
! Calculate (dd|dp) integrals
!
! In  : exfac12, exfac34 (exponents and coefficients of primitive pair functions)
!       xyziq    (x,y,z elements of i and q)
!       xzkl     (x,z elements of k and l)
!       nij, nkl (number of primitive pair functions)
! Out : phmdint  (two-electron integral)
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nijkl(2), nbfijkl(4)
      integer :: ij, kl, igrid, i, j, k, l, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, sqrtpi4=0.8862269254527580D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00, four=4.0D+00, five=5.0D+00
      real(8),parameter :: six=6.0D+00, seven=7.0D+00, eight=8.0D+00, nine=9.0D+00, ten=1.0D+01
      real(8),parameter :: p11=1.1D+01, p12=1.2D+01, p13=1.3D+01, p15=1.5D+01, p16=1.6D+01
      real(8),parameter :: p18=1.8D+01, p20=2.0D+01, p21=2.1D+01, p24=2.4D+01, p28=2.8D+01
      real(8),parameter :: p30=3.0D+01, p36=3.6D+01, p45=4.5D+01, p105=1.05D+2, p210=2.1D+02
      real(8),parameter :: p420=4.2D+02, sqrt3=1.73205080756888D+00, sqrt3h=8.660254037844386D-01
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*), xzkl(2), rot(3,3)
      real(8),intent(out) :: phmdint(6,6,6,6)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:8)
      real(8) :: f0(5), f1(2,13), f2(3,16), f3(4,16), f4(5,16), f5(6,11), f6(7,7), f7(8,3)
      real(8) :: f8(9), ftw(8,16), r0(25), r1(3,40), r2(6,56), r3(10,52), r4(15,42), r5(21,24)
      real(8) :: r6(28,12), r7(36,4), r8(45)
      real(8) :: ex12, ex34, ex43, ex14, ex41, expq, expq2, expq3, expq4, ex3q, ex4q, c12, c34
      real(8) :: zip, zip2, ziq, xiq, yiq, xiq2, yiq2, xyiq, xiq4, yiq4, xyiq2, xiq6, yiq6
      real(8) :: x4y2, x2y4, xypq2, xiq8, yiq8, x6y2, x4y4, x2y6
      real(8) :: zpq, zpq2, fac, ex33q, ex34q, ex44q, zjp, zjp2
      real(8) :: pmd, pmd2, pmd3, qmd, qmd2, qmd3, qmd4, qmd4x, qmd4y, qmd4xy, qx, qz
      real(8) :: eri(6,6,6,6), work(15), rot2(6,6), rot3(6,5)
      real(8) :: f1w(3,9), f2w(6,13), f3w(10,15), f4w(15,16), f5w(21,11), f6w(28,7), f7w(36,3)
!
! Zero-clear
!
      r0(1:25)     = zero
      r1(1:3 ,1:40)= zero
      r2(1:6 ,1:56)= zero
      r3(1:10,1:52)= zero
      r4(1:15,1:42)= zero
      r5(1:21,1:24)= zero
      r6(1:28,1:12)= zero
      r7(1:36,1:4) = zero
      r8(1:45)     = zero
!
      do kl= 1,nijkl(2)
        ex34= exfac2(1,kl)
        ex43= exfac2(2,kl)
        ex3q= exfac2(3,kl)
        ex4q= exfac2(4,kl)
        c34 = exfac2(5,kl)
        xiq = xyziq(1,kl)
        yiq = xyziq(2,kl)
        ziq = xyziq(3,kl)
        xiq2= xiq*xiq
        yiq2= yiq*yiq
        xyiq= xiq*yiq
        xiq4= xiq2*xiq2
        yiq4= yiq2*yiq2
        xyiq2= xiq2*yiq2
        xiq6 = xiq4*xiq2
        yiq6 = yiq4*yiq2
        xiq8 = xiq4*xiq4
        yiq8 = yiq4*yiq4
        x4y2 = xiq4*yiq2
        x2y4 = xiq2*yiq4
        x6y2 = xiq6*yiq2
        x4y4 = xiq4*yiq4
        x2y6 = xiq2*yiq6
        xypq2= xiq2+yiq2
        f0(1:5)     = zero
        f1(1:2,1:13)= zero
        f2(1:3,1:16)= zero
        f3(1:4,1:16)= zero
        f4(1:5,1:16)= zero
        f5(1:6,1:11)= zero
        f6(1:7,1:7) = zero
        f7(1:8,1:3) = zero
        f8(1:9)     = zero
        do ij= 1,nijkl(1)
          ex12= exfac1(1,ij)
          pmd = exfac1(2,ij)
          zjp = exfac1(3,ij)
          zip = exfac1(4,ij)
          c12 = exfac1(5,ij)
          ex14= ex12+ex34
          zpq = ziq-zip
          zpq2= zpq*zpq
          expq= ex12*ex34
          tval=(xypq2+zpq2)*expq
!
! Calculate Fm(T)
!
          if(tval >= threshtval*ex14) then
            tinv= one/sqrt(tval)
            ft(0)= c12*sqrtpi4*tinv
            expq= expq*tinv*tinv
            ft(1)= ft(0)*expq
            ft(2)= ft(1)*expq*three
            ft(3)= ft(2)*expq*five
            ft(4)= ft(3)*expq*seven
            ft(5)= ft(4)*expq*nine
            ft(6)= ft(5)*expq*p11
            ft(7)= ft(6)*expq*p13
            ft(8)= ft(7)*expq*p15
          else
            ex41= one/sqrt(ex14)
            tval= tval*ex41*ex41
            igrid= int(tval)
            tval2= tval *tval
            tval3= tval2*tval
            tval4= tval2*tval2
            tval5= tval2*tval3
            tval6= tval3*tval3
            tval7= tval4*tval3
            tval8= tval4*tval4
            tval9= tval4*tval5
            tval10=tval5*tval5
            do ii= 0,8
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
            fac= c12*ex41
            expq= expq*two*ex41*ex41
            expq2= expq*expq
            expq3= expq2*expq
            expq4= expq2*expq2
            ft(0)= ft(0)*fac
            ft(1)= ft(1)*fac*expq
            ft(2)= ft(2)*fac*expq2
            ft(3)= ft(3)*fac*expq3
            ft(4)= ft(4)*fac*expq4
            ft(5)= ft(5)*fac*expq4*expq
            ft(6)= ft(6)*fac*expq4*expq2
            ft(7)= ft(7)*fac*expq4*expq3
            ft(8)= ft(8)*fac*expq4*expq4
          endif
          work(2)= zpq
          work(3)= zpq2
          work(4)= zpq2*zpq
          work(5)= zpq2*zpq2
          work(6)= work(5)*zpq
          work(7)= work(5)*zpq2
          work(8)= work(6)*zpq2
          pmd2= pmd*pmd
          pmd3= pmd*pmd*pmd
          zip2=zip*zip
          zjp2=zjp*zjp
          ftw(1, 1)= zjp2*zip2
          ftw(1, 2)= pmd*zjp2
          ftw(1, 3)= pmd*zjp*zip
          ftw(1, 4)= pmd*zip2
          ftw(1, 5)= pmd2
          ftw(1, 6)= pmd*zjp2*zip
          ftw(1, 7)= pmd*zip2*zjp
          ftw(1, 8)= pmd2*zjp
          ftw(1, 9)= pmd2*zip
          ftw(1,10)= pmd2*zjp2
          ftw(1,11)= pmd2*zjp*zip
          ftw(1,12)= pmd2*zip2
          ftw(1,13)= pmd3
          ftw(1,14)= pmd3*zjp
          ftw(1,15)= pmd3*zip
          ftw(1,16)= pmd3*pmd
          do i= 1,5
            do j= 2,5
              ftw(j,i)= ftw(1,i)*work(j)
            enddo
          enddo
          do i= 6,9
            do j= 2,6
              ftw(j,i)= ftw(1,i)*work(j)
            enddo
          enddo
          do i= 10,13
            do j= 2,7
              ftw(j,i)= ftw(1,i)*work(j)
            enddo
          enddo
          do i= 14,16
            do j= 2,8
              ftw(j,i)= ftw(1,i)*work(j)
            enddo
          enddo
          f0(1)= f0(1)+ft(0)*ftw(1,1)
          f0(2)= f0(2)+ft(0)*ftw(1,2)
          f0(3)= f0(3)+ft(0)*ftw(1,3)
          f0(4)= f0(4)+ft(0)*ftw(1,4)
          f0(5)= f0(5)+ft(0)*ftw(1,5)
          do i= 1,13
            f1(1,i)= f1(1,i)-ft(1)*ftw(1,i)
            f1(2,i)= f1(2,i)-ft(1)*ftw(2,i)
          enddo
          do i= 1,16
            f2(1,i)= f2(1,i)+ft(2)*ftw(1,i)
            f2(2,i)= f2(2,i)+ft(2)*ftw(2,i)
            f2(3,i)= f2(3,i)+ft(2)*ftw(3,i)
          enddo
          do i= 1,16
            do j= 1,4
              f3(j,i)= f3(j,i)-ft(3)*ftw(j,i)
            enddo
          enddo
          do i= 1,16
            do j= 1,5
              f4(j,i)= f4(j,i)+ft(4)*ftw(j,i)
            enddo
          enddo
          do i= 1,11
            do j= 1,6
              f5(j,i)= f5(j,i)-ft(5)*ftw(j,i+5)
            enddo
          enddo
          do i= 1,7
            do j= 1,7
              f6(j,i)= f6(j,i)+ft(6)*ftw(j,i+9)
            enddo
          enddo
          do i= 1,3
            do j= 1,8
              f7(j,i)= f7(j,i)-ft(7)*ftw(j,i+13)
            enddo
          enddo
          do j= 1,8
            f8(j)= f8(j)+ft(8)*ftw(j,16)
          enddo
          f8(9)= f8(9)+ft(8)*ftw(8,16)*zpq
        enddo
!
        qmd = ex43*c34
        qmd2= qmd*ex43
        qmd3= qmd*ex43*ex43
        qmd4= qmd*ex43*ex43*ex43
        ex33q= ex3q*ex3q
        ex34q= ex3q*ex4q
        ex44q= ex4q*ex4q
        qmd4x= qmd4*xiq
        qmd4y= qmd4*yiq
        qmd4xy=qmd4*xiq*yiq
!
        work( 1)= qmd2
        work( 2)= qmd*ex33q
        work( 3)= qmd*ex34q
        work( 4)= qmd*ex44q
        work( 5)= ex33q*ex44q*c34
        work( 6)= qmd2*ex3q
        work( 7)= qmd2*ex4q
        work( 8)= qmd*ex34q*ex3q
        work( 9)= qmd*ex34q*ex4q
        work(10)= qmd3
        work(11)= qmd2*ex33q
        work(12)= qmd2*ex34q
        work(13)= qmd2*ex44q
        work(14)= qmd3*ex3q
        work(15)= qmd3*ex4q
!
        do i= 1,5
          r0(i   )= r0(i   )+f0(1)*work(i)
          r0(i+5 )= r0(i+5 )+f0(2)*work(i)
          r0(i+10)= r0(i+10)+f0(3)*work(i)
          r0(i+15)= r0(i+15)+f0(4)*work(i)
          r0(i+20)= r0(i+20)+f0(5)*work(i)
        enddo
!
        do i= 1,9
          f1w(1,i)= f1(1,i)*xiq
          f1w(2,i)= f1(1,i)*yiq
          f1w(3,i)= f1(2,i)
        enddo
        do i= 1,4
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,1)*work(i+5)
          enddo
        enddo
        do i= 5,8
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,2)*work(i+1)
          enddo
        enddo
        do i= 9,12
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,3)*work(i-3)
          enddo
        enddo
        do i= 13,16
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,4)*work(i-7)
          enddo
        enddo
        do i= 17,20
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,5)*work(i-11)
          enddo
        enddo
        do i= 21,25
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,6)*work(i-20)
          enddo
        enddo
        do i= 26,30
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,7)*work(i-25)
          enddo
        enddo
        do i= 31,35
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,8)*work(i-30)
          enddo
        enddo
        do i= 36,40
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,9)*work(i-35)
          enddo
        enddo
!
        do i= 1,13
          f2w(1,i)= f2(1,i)*xiq2+f1(1,i)
          f2w(2,i)= f2(1,i)*yiq2+f1(1,i)
          f2w(3,i)= f2(3,i)     +f1(1,i)
          f2w(4,i)= f2(1,i)*xyiq
          f2w(5,i)= f2(2,i)*xiq
          f2w(6,i)= f2(2,i)*yiq
        enddo
        do i= 1,4
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,1)*work(i+9)
          enddo
        enddo
        do i= 5,8
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,2)*work(i+5)
          enddo
        enddo
        do i= 9,12
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,3)*work(i+1)
          enddo
        enddo
        do i= 13,16
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,4)*work(i-3)
          enddo
        enddo
        do i= 17,20
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,5)*work(i-7)
          enddo
        enddo
        do i= 21,24
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,6)*work(i-15)
          enddo
        enddo
        do i= 25,28
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,7)*work(i-19)
          enddo
        enddo
        do i= 29,32
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,8)*work(i-23)
          enddo
        enddo
        do i= 33,36
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,9)*work(i-27)
          enddo
        enddo
        do i= 37,41
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,10)*work(i-36)
          enddo
        enddo
        do i= 42,46
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,11)*work(i-41)
          enddo
        enddo
        do i= 47,51
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,12)*work(i-46)
          enddo
        enddo
        do i= 52,56
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,13)*work(i-51)
          enddo
        enddo
! 
        do i= 1,15
          f3w( 1,i)=(f3(1,i)*xiq2+f2(1,i)*three)*xiq
          f3w( 2,i)=(f3(1,i)*xiq2+f2(1,i)      )*yiq
          f3w( 3,i)=(f3(2,i)*xiq2+f2(2,i)      )
          f3w( 4,i)=(f3(1,i)*yiq2+f2(1,i)      )*xiq
          f3w( 5,i)=(f3(2,i)*xyiq              )
          f3w( 6,i)=(f3(3,i)     +f2(1,i)      )*xiq
          f3w( 7,i)=(f3(1,i)*yiq2+f2(1,i)*three)*yiq
          f3w( 8,i)=(f3(2,i)*yiq2+f2(2,i)      )
          f3w( 9,i)=(f3(3,i)     +f2(1,i)      )*yiq
          f3w(10,i)=(f3(4,i)     +f2(2,i)*three)
        enddo
        do i= 1,2
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,1)*work(i+13)
          enddo
        enddo
        do i= 3,4
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,2)*work(i+11)
          enddo
        enddo
        do i= 5,6
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,3)*work(i+9)
          enddo
        enddo
        do i= 7,8
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,4)*work(i+7)
          enddo
        enddo
        do i= 9,10
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,5)*work(i+5)
          enddo
        enddo
        do i= 11,14
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,6)*work(i-1)
          enddo
        enddo
        do i= 15,18
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,7)*work(i-5)
          enddo
        enddo
        do i= 19,22
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,8)*work(i-9)
          enddo
        enddo
        do i= 23,26
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,9)*work(i-13)
          enddo
        enddo
        do i= 27,30
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,10)*work(i-21)
          enddo
        enddo
        do i= 31,34
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,11)*work(i-25)
          enddo
        enddo
        do i= 35,38
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,12)*work(i-29)
          enddo
        enddo
        do i= 39,42
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,13)*work(i-33)
          enddo
        enddo
        do i= 43,47
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,14)*work(i-42)
          enddo
        enddo
        do i= 48,52
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,15)*work(i-47)
          enddo
        enddo
!
        do i= 1,16
          f4w( 1,i)=(f4(1,i)*xiq4 +f3(1,i)*xiq2*six         +f2(1,i)*three)
          f4w( 2,i)=(f4(1,i)*xiq2 +f3(1,i)*three                          )*xyiq
          f4w( 3,i)=(f4(2,i)*xiq2 +f3(2,i)*three                          )*xiq
          f4w( 4,i)=(f4(1,i)*xyiq2+f3(1,i)*xiq2+f3(1,i)*yiq2+f2(1,i)      )
          f4w( 5,i)=(f4(2,i)*xiq2 +f3(2,i)                                )*yiq
          f4w( 6,i)=(f4(3,i)*xiq2 +f3(1,i)*xiq2+f3(3,i)     +f2(1,i)      )
          f4w( 7,i)=(f4(1,i)*yiq2 +f3(1,i)*three                          )*xyiq
          f4w( 8,i)=(f4(2,i)*yiq2 +f3(2,i)                                )*xiq
          f4w( 9,i)=(f4(3,i)      +f3(1,i)                                )*xyiq
          f4w(10,i)=(f4(4,i)      +f3(2,i)*three                          )*xiq
          f4w(11,i)=(f4(1,i)*yiq4 +f3(1,i)*yiq2*six         +f2(1,i)*three)
          f4w(12,i)=(f4(2,i)*yiq2 +f3(2,i)*three                          )*yiq
          f4w(13,i)=(f4(3,i)*yiq2 +f3(1,i)*yiq2+f3(3,i)     +f2(1,i)      )
          f4w(14,i)=(f4(4,i)      +f3(2,i)*three                          )*yiq
          f4w(15,i)=(f4(5,i)      +f3(3,i)*six              +f2(1,i)*three)
        enddo
        do i= 1,5
          do j= 1,15
            r4(j,i)= r4(j,i)+f4w(j,i)*qmd4
          enddo
        enddo
        do i= 6,7
          do j= 1,15
            r4(j,i)= r4(j,i)+f4w(j,6)*work(i+8)
          enddo
        enddo
        do i= 8,9
          do j= 1,15
            r4(j,i)= r4(j,i)+f4w(j,7)*work(i+6)
          enddo
        enddo
        do i= 10,11
          do j= 1,15
            r4(j,i)= r4(j,i)+f4w(j,8)*work(i+4)
          enddo
        enddo
        do i= 12,13
          do j= 1,15
            r4(j,i)= r4(j,i)+f4w(j,9)*work(i+2)
          enddo
        enddo
        do i= 14,17
          do j= 1,15
            r4(j,i)= r4(j,i)+f4w(j,10)*work(i-4)
          enddo
        enddo
        do i= 18,21
          do j= 1,15
            r4(j,i)= r4(j,i)+f4w(j,11)*work(i-8)
          enddo
        enddo
        do i= 22,25
          do j= 1,15
            r4(j,i)= r4(j,i)+f4w(j,12)*work(i-12)
          enddo
        enddo
        do i= 26,29
          do j= 1,15
            r4(j,i)= r4(j,i)+f4w(j,13)*work(i-16)
          enddo
        enddo
        do i= 30,33
          do j= 1,15
            r4(j,i)= r4(j,i)+f4w(j,14)*work(i-24)
          enddo
        enddo
        do i= 34,37
          do j= 1,15
            r4(j,i)= r4(j,i)+f4w(j,15)*work(i-28)
          enddo
        enddo
        do i= 38,42
          do j= 1,15
            r4(j,i)= r4(j,i)+f4w(j,16)*work(i-37)
          enddo
        enddo
!
        do i= 1,11
          f5w( 1,i)=(f5(1,i)*xiq4 +f4(1,i+5)*xiq2*ten                 +f3(1,i+5)*p15  )*xiq
          f5w( 2,i)=(f5(1,i)*xiq4 +f4(1,i+5)*xiq2*six                 +f3(1,i+5)*three)*yiq
          f5w( 3,i)=(f5(2,i)*xiq4 +f4(2,i+5)*xiq2*six                 +f3(2,i+5)*three)
          f5w( 4,i)=(f5(1,i)*xyiq2+f4(1,i+5)*xiq2+f4(1,i+5)*yiq2*three+f3(1,i+5)*three)*xiq
          f5w( 5,i)=(f5(2,i)*xiq2 +f4(2,i+5)*three                                    )*xyiq
          f5w( 6,i)=(f5(3,i)*xiq2 +f4(1,i+5)*xiq2+f4(3,i+5)*three     +f3(1,i+5)*three)*xiq
          f5w( 7,i)=(f5(1,i)*xyiq2+f4(1,i+5)*xiq2*three+f4(1,i+5)*yiq2+f3(1,i+5)*three)*yiq
          f5w( 8,i)=(f5(2,i)*xyiq2+f4(2,i+5)*xiq2+f4(2,i+5)*yiq2      +f3(2,i+5)      )
          f5w( 9,i)=(f5(3,i)*xiq2 +f4(1,i+5)*xiq2+f4(3,i+5)           +f3(1,i+5)      )*yiq
          f5w(10,i)=(f5(4,i)*xiq2 +f4(2,i+5)*xiq2*three+f4(4,i+5)     +f3(2,i+5)*three)
          f5w(11,i)=(f5(1,i)*yiq4 +f4(1,i+5)*yiq2*six                 +f3(1,i+5)*three)*xiq
          f5w(12,i)=(f5(2,i)*yiq2 +f4(2,i+5)*three                                    )*xyiq
          f5w(13,i)=(f5(3,i)*yiq2 +f4(1,i+5)*yiq2+f4(3,i+5)           +f3(1,i+5)      )*xiq
          f5w(14,i)=(f5(4,i)      +f4(2,i+5)*three                                    )*xyiq
          f5w(15,i)=(f5(5,i)      +f4(3,i+5)*six                      +f3(1,i+5)*three)*xiq
          f5w(16,i)=(f5(1,i)*yiq4 +f4(1,i+5)*yiq2*ten                 +f3(1,i+5)*p15  )*yiq
          f5w(17,i)=(f5(2,i)*yiq4 +f4(2,i+5)*yiq2*six                 +f3(2,i+5)*three)
          f5w(18,i)=(f5(3,i)*yiq2 +f4(1,i+5)*yiq2+f4(3,i+5)*three     +f3(1,i+5)*three)*yiq
          f5w(19,i)=(f5(4,i)*yiq2 +f4(2,i+5)*yiq2*three+f4(4,i+5)     +f3(2,i+5)*three)
          f5w(20,i)=(f5(5,i)      +f4(3,i+5)*six                      +f3(1,i+5)*three)*yiq
          f5w(21,i)=(f5(6,i)      +f4(4,i+5)*ten                      +f3(2,i+5)*p15  )
        enddo
        do i= 1,4
          do j= 1,21
            r5(j,i)= r5(j,i)+f5w(j,i)*qmd4
          enddo
        enddo
        do i= 5,6
          do j= 1,21
            r5(j,i)= r5(j,i)+f5w(j,5)*work(i+9)
          enddo
        enddo
        do i= 7,8
          do j= 1,21
            r5(j,i)= r5(j,i)+f5w(j,6)*work(i+7)
          enddo
        enddo
        do i= 9,10
          do j= 1,21
            r5(j,i)= r5(j,i)+f5w(j,7)*work(i+5)
          enddo
        enddo
        do i= 11,12
          do j= 1,21
            r5(j,i)= r5(j,i)+f5w(j,8)*work(i+3)
          enddo
        enddo
        do i= 13,16
          do j= 1,21
            r5(j,i)= r5(j,i)+f5w(j,9)*work(i-3)
          enddo
        enddo
        do i= 17,20
          do j= 1,21
            r5(j,i)= r5(j,i)+f5w(j,10)*work(i-7)
          enddo
        enddo
        do i= 21,24
          do j= 1,21
            r5(j,i)= r5(j,i)+f5w(j,11)*work(i-15)
          enddo
        enddo
!
        do i= 1,7
          f6w( 1,i)=(f6(1,i)*xiq6+f5(1,i+4)*xiq4*p15+f4(1,i+9)*xiq2*p45 &
&                   +f3(1,i+9)*p15)
          f6w( 2,i)=(f6(1,i)*xiq4+f5(1,i+4)*xiq2*ten+f4(1,i+9)*p15)*xyiq
          f6w( 3,i)=(f6(2,i)*xiq4+f5(2,i+4)*xiq2*ten+f4(2,i+9)*p15)*xiq
          f6w( 4,i)=(f6(1,i)*x4y2+f5(1,i+4)*xiq4+f5(1,i+4)*xyiq2*six &
&                   +f4(1,i+9)*xiq2*six+f4(1,i+9)*yiq2*three+f3(1,i+9)*three)
          f6w( 5,i)=(f6(2,i)*xiq4+f5(2,i+4)*xiq2*six+f4(2,i+9)*three)*yiq
          f6w( 6,i)=(f6(3,i)*xiq4+f5(1,i+4)*xiq4+f5(3,i+4)*xiq2*six &
&                   +f4(1,i+9)*xiq2*six+f4(3,i+9)*three+f3(1,i+9)*three)
          f6w( 7,i)=(f6(1,i)*xyiq2+f5(1,i+4)*xiq2*three+f5(1,i+4)*yiq2*three &
&                   +f4(1,i+9)*nine)*xyiq
          f6w( 8,i)=(f6(2,i)*xyiq2+f5(2,i+4)*xiq2+f5(2,i+4)*yiq2*three &
&                   +f4(2,i+9)*three)*xiq
          f6w( 9,i)=(f6(3,i)*xiq2+f5(1,i+4)*xiq2+f5(3,i+4)*three &
&                   +f4(1,i+9)*three)*xyiq
          f6w(10,i)=(f6(4,i)*xiq2+f5(2,i+4)*xiq2*three+f5(4,i+4)*three &
&                   +f4(2,i+9)*nine)*xiq
          f6w(11,i)=(f6(1,i)*x2y4+f5(1,i+4)*xyiq2*six+f5(1,i+4)*yiq4 &
&                   +f4(1,i+9)*xiq2*three+f4(1,i+9)*yiq2*six+f3(1,i+9)*three)
          f6w(12,i)=(f6(2,i)*xyiq2+f5(2,i+4)*xiq2*three+f5(2,i+4)*yiq2 &
&                   +f4(2,i+9)*three)*yiq
          f6w(13,i)=(f6(3,i)*xyiq2+f5(1,i+4)*xyiq2+f5(3,i+4)*xiq2+f5(3,i+4)*yiq2 &
&                   +f4(1,i+9)*xiq2+f4(1,i+9)*yiq2+f4(3,i+9)+f3(1,i+9))
          f6w(14,i)=(f6(4,i)*xiq2+f5(2,i+4)*xiq2*three+f5(4,i+4) &
&                   +f4(2,i+9)*three)*yiq
          f6w(15,i)=(f6(5,i)*xiq2+f5(3,i+4)*xiq2*six+f5(5,i+4) &
&                   +f4(1,i+9)*xiq2*three+f4(3,i+9)*six+f3(1,i+9)*three)
          f6w(16,i)=(f6(1,i)*yiq4+f5(1,i+4)*yiq2*ten+f4(1,i+9)*p15)*xyiq
          f6w(17,i)=(f6(2,i)*yiq4+f5(2,i+4)*yiq2*six+f4(2,i+9)*three)*xiq
          f6w(18,i)=(f6(3,i)*yiq2+f5(1,i+4)*yiq2+f5(3,i+4)*three &
&                   +f4(1,i+9)*three)*xyiq
          f6w(19,i)=(f6(4,i)*yiq2+f5(2,i+4)*yiq2*three+f5(4,i+4) &
&                   +f4(2,i+9)*three)*xiq
          f6w(20,i)=(f6(5,i)+f5(3,i+4)*six+f4(1,i+9)*three)*xyiq
          f6w(21,i)=(f6(6,i)+f5(4,i+4)*ten+f4(2,i+9)*p15)*xiq
          f6w(22,i)=(f6(1,i)*yiq6+f5(1,i+4)*yiq4*p15+f4(1,i+9)*yiq2*p45 &
&                   +f3(1,i+9)*p15)
          f6w(23,i)=(f6(2,i)*yiq4+f5(2,i+4)*yiq2*ten+f4(2,i+9)*p15)*yiq
          f6w(24,i)=(f6(3,i)*yiq4+f5(1,i+4)*yiq4+f5(3,i+4)*yiq2*six &
&                   +f4(1,i+9)*yiq2*six+f4(3,i+9)*three+f3(1,i+9)*three)
          f6w(25,i)=(f6(4,i)*yiq2+f5(2,i+4)*yiq2*three+f5(4,i+4)*three &
&                   +f4(2,i+9)*nine)*yiq
          f6w(26,i)=(f6(5,i)*yiq2+f5(3,i+4)*yiq2*six+f5(5,i+4) &
&                   +f4(1,i+9)*yiq2*three+f4(3,i+9)*six+f3(1,i+9)*three)
          f6w(27,i)=(f6(6,i)+f5(4,i+4)*ten+f4(2,i+9)*p15)*yiq
          f6w(28,i)=(f6(7,i)+f5(5,i+4)*p15+f4(3,i+9)*p45+f3(1,i+9)*p15)
        enddo
        do i= 1,4
          do j= 1,28
            r6(j,i)= r6(j,i)+f6w(j,i)*qmd4
          enddo
        enddo
        do i= 5,6
          do j= 1,28
            r6(j,i)= r6(j,i)+f6w(j,5)*work(i+9)
          enddo
        enddo
        do i= 7,8
          do j= 1,28
            r6(j,i)= r6(j,i)+f6w(j,6)*work(i+7)
          enddo
        enddo
        do i= 9,12
          do j= 1,28
            r6(j,i)= r6(j,i)+f6w(j,7)*work(i+1)
          enddo
        enddo
!
        do i=1,3
          f7w( 1,i)=(f7(1,i)*xiq6+f6(1,i+4)*xiq4*p21+f5(1,i+8)*xiq2*p105+f4(1,i+13)*p105)*xiq
          f7w( 2,i)=(f7(1,i)*xiq6+f6(1,i+4)*xiq4*p15+f5(1,i+8)*xiq2*p45+f4(1,i+13)*p15)*yiq
          f7w( 3,i)=(f7(2,i)*xiq6+f6(2,i+4)*xiq4*p15+f5(2,i+8)*xiq2*p45+f4(2,i+13)*p15)
          f7w( 4,i)=(f7(1,i)*x4y2+f6(1,i+4)*xiq4+f6(1,i+4)*xyiq2*ten+f5(1,i+8)*xiq2*ten &
&                   +f5(1,i+8)*yiq2*p15+f4(1,i+13)*p15)*xiq
          f7w( 5,i)=(f7(2,i)*xiq4+f6(2,i+4)*xiq2*ten+f5(2,i+8)*p15)*xyiq
          f7w( 6,i)=(f7(3,i)*xiq4+f6(1,i+4)*xiq4+f6(3,i+4)*xiq2*ten+f5(1,i+8)*xiq2*ten &
&                   +f5(3,i+8)*p15+f4(1,i+13)*p15)*xiq
          f7w( 7,i)=(f7(1,i)*x4y2+f6(1,i+4)*xiq4*three+f6(1,i+4)*xyiq2*six+f5(1,i+8)*xiq2*p18 &
&                   +f5(1,i+8)*yiq2*three+f4(1,i+13)*nine)*yiq
          f7w( 8,i)=(f7(2,i)*x4y2+f6(2,i+4)*xiq4+f6(2,i+4)*xyiq2*six+f5(2,i+8)*xiq2*six &
&                   +f5(2,i+8)*yiq2*three+f4(2,i+13)*three)
          f7w( 9,i)=(f7(3,i)*xiq4+f6(1,i+4)*xiq4+f6(3,i+4)*xiq2*six+f5(1,i+8)*xiq2*six &
&                   +f5(3,i+8)*three+f4(1,i+13)*three)*yiq
          f7w(10,i)=(f7(4,i)*xiq4+f6(2,i+4)*xiq4*three+f6(4,i+4)*xiq2*six+f5(2,i+8)*xiq2*p18 &
&                   +f5(4,i+8)*three+f4(2,i+13)*nine)
          f7w(11,i)=(f7(1,i)*x2y4+f6(1,i+4)*xyiq2*six+f6(1,i+4)*yiq4*three+f5(1,i+8)*xiq2*three &
&                   +f5(1,i+8)*yiq2*p18+f4(1,i+13)*nine)*xiq
          f7w(12,i)=(f7(2,i)*xyiq2+f6(2,i+4)*xiq2*three+f6(2,i+4)*yiq2*three+f5(2,i+8)*nine)*xyiq
          f7w(13,i)=(f7(3,i)*xyiq2+f6(1,i+4)*xyiq2+f6(3,i+4)*xiq2+f6(3,i+4)*yiq2*three &
&                   +f5(1,i+8)*xiq2+f5(1,i+8)*yiq2*three+f5(3,i+8)*three+f4(1,i+13)*three)*xiq
          f7w(14,i)=(f7(4,i)*xiq2+f6(2,i+4)*xiq2*three+f6(4,i+4)*three+f5(2,i+8)*nine)*xyiq
          f7w(15,i)=(f7(5,i)*xiq2+f6(3,i+4)*xiq2*six+f6(5,i+4)*three+f5(1,i+8)*xiq2*three &
&                   +f5(3,i+8)*p18+f4(1,i+13)*nine)*xiq
          f7w(16,i)=(f7(1,i)*x2y4+f6(1,i+4)*xyiq2*ten+f6(1,i+4)*yiq4+f5(1,i+8)*xiq2*p15 &
&                   +f5(1,i+8)*yiq2*ten+f4(1,i+13)*p15)*yiq
          f7w(17,i)=(f7(2,i)*x2y4+f6(2,i+4)*xyiq2*six+f6(2,i+4)*yiq4+f5(2,i+8)*xiq2*three &
&                   +f5(2,i+8)*yiq2*six+f4(2,i+13)*three)
          f7w(18,i)=(f7(3,i)*xyiq2+f6(1,i+4)*xyiq2+f6(3,i+4)*xiq2*three+f6(3,i+4)*yiq2 &
&                   +f5(1,i+8)*xiq2*three+f5(1,i+8)*yiq2+f5(3,i+8)*three+f4(1,i+13)*three)*yiq
          f7w(19,i)=(f7(4,i)*xyiq2+f6(2,i+4)*xyiq2*three+f6(4,i+4)*xiq2+f6(4,i+4)*yiq2 &
&                   +f5(2,i+8)*xiq2*three+f5(2,i+8)*yiq2*three+f5(4,i+8)+f4(2,i+13)*three)
          f7w(20,i)=(f7(5,i)*xiq2+f6(3,i+4)*xiq2*six+f6(5,i+4)+f5(1,i+8)*xiq2*three+f5(3,i+8)*six &
&                   +f4(1,i+13)*three)*yiq
          f7w(21,i)=(f7(6,i)*xiq2+f6(4,i+4)*xiq2*ten+f6(6,i+4)+f5(2,i+8)*xiq2*p15+f5(4,i+8)*ten &
&                   +f4(2,i+13)*p15)
          f7w(22,i)=(f7(1,i)*yiq6+f6(1,i+4)*yiq4*p15+f5(1,i+8)*yiq2*p45+f4(1,i+13)*p15)*xiq
          f7w(23,i)=(f7(2,i)*yiq4+f6(2,i+4)*yiq2*ten+f5(2,i+8)*p15)*xyiq
          f7w(24,i)=(f7(3,i)*yiq4+f6(1,i+4)*yiq4+f6(3,i+4)*yiq2*six+f5(1,i+8)*yiq2*six &
&                   +f5(3,i+8)*three+f4(1,i+13)*three)*xiq
          f7w(25,i)=(f7(4,i)*yiq2+f6(2,i+4)*yiq2*three+f6(4,i+4)*three+f5(2,i+8)*nine)*xyiq
          f7w(26,i)=(f7(5,i)*yiq2+f6(3,i+4)*yiq2*six+f6(5,i+4)+f5(1,i+8)*yiq2*three+f5(3,i+8)*six &
&                   +f4(1,i+13)*three)*xiq
          f7w(27,i)=(f7(6,i)+f6(4,i+4)*ten+f5(2,i+8)*p15)*xyiq
          f7w(28,i)=(f7(7,i)+f6(5,i+4)*p15+f5(3,i+8)*p45+f4(1,i+13)*p15)*xiq
          f7w(29,i)=(f7(1,i)*yiq6+f6(1,i+4)*yiq4*p21+f5(1,i+8)*yiq2*p105+f4(1,i+13)*p105)*yiq
          f7w(30,i)=(f7(2,i)*yiq6+f6(2,i+4)*yiq4*p15+f5(2,i+8)*yiq2*p45+f4(2,i+13)*p15)
          f7w(31,i)=(f7(3,i)*yiq4+f6(1,i+4)*yiq4+f6(3,i+4)*yiq2*ten+f5(1,i+8)*yiq2*ten &
&                   +f5(3,i+8)*p15+f4(1,i+13)*p15)*yiq
          f7w(32,i)=(f7(4,i)*yiq4+f6(2,i+4)*yiq4*three+f6(4,i+4)*yiq2*six+f5(2,i+8)*yiq2*p18 &
&                   +f5(4,i+8)*three+f4(2,i+13)*nine)
          f7w(33,i)=(f7(5,i)*yiq2+f6(3,i+4)*yiq2*six+f6(5,i+4)*three+f5(1,i+8)*yiq2*three &
&                   +f5(3,i+8)*p18+f4(1,i+13)*nine)*yiq
          f7w(34,i)=(f7(6,i)*yiq2+f6(4,i+4)*yiq2*ten+f6(6,i+4)+f5(2,i+8)*yiq2*p15+f5(4,i+8)*ten &
&                   +f4(2,i+13)*p15)
          f7w(35,i)=(f7(7,i)+f6(5,i+4)*p15+f5(3,i+8)*p45+f4(1,i+13)*p15)*yiq
          f7w(36,i)=(f7(8,i)+f6(6,i+4)*p21+f5(4,i+8)*p105+f4(2,i+13)*p105)
        enddo
        do j= 1,36
          r7(j,1)= r7(j,1)+f7w(j,1)*qmd4
          r7(j,2)= r7(j,2)+f7w(j,2)*qmd4
          r7(j,3)= r7(j,3)+f7w(j,3)*work(14)
          r7(j,4)= r7(j,4)+f7w(j,3)*work(15)
        enddo
!
        r8( 1)= r8( 1)+(f8(1)*xiq8+f7(1,3)*xiq6*p28+f6(1,7)*xiq4*p210+f5(1,11)*xiq2*p420 &
&                     +f4(1,16)*p105)*qmd4
        r8( 2)= r8( 2)+(f8(1)*xiq6+f7(1,3)*xiq4*p21+f6(1,7)*xiq2*p105+f5(1,11)*p105)*qmd4xy
        r8( 3)= r8( 3)+(f8(2)*xiq6+f7(2,3)*xiq4*p21+f6(2,7)*xiq2*p105+f5(2,11)*p105)*qmd4x
        r8( 4)= r8( 4)+(f8(1)*x6y2+f7(1,3)*xiq6+f7(1,3)*x4y2*p15+f6(1,7)*xiq4*p15 &
&                    +f6(1,7)*xyiq2*p45+f5(1,11)*xiq2*p45+f5(1,11)*yiq2*p15+f4(1,16)*p15)*qmd4
        r8( 5)= r8( 5)+(f8(2)*xiq6+f7(2,3)*xiq4*p15+f6(2,7)*xiq2*p45+f5(2,11)*p15)*qmd4y
        r8( 6)= r8( 6)+(f8(3)*xiq6+f7(1,3)*xiq6+f7(3,3)*xiq4*p15+f6(1,7)*xiq4*p15 &
&                    +f6(3,7)*xiq2*p45+f5(1,11)*xiq2*p45+f5(3,11)*p15+f4(1,16)*p15)*qmd4
        r8( 7)= r8( 7)+(f8(1)*x4y2+f7(1,3)*xiq4*three+f7(1,3)*xyiq2*ten+f6(1,7)*xiq2*p30 &
&                    +f6(1,7)*yiq2*p15+f5(1,11)*p45)*qmd4xy
        r8( 8)= r8( 8)+(f8(2)*x4y2+f7(2,3)*xiq4+f7(2,3)*xyiq2*ten+f6(2,7)*xiq2*ten &
&                    +f6(2,7)*yiq2*p15+f5(2,11)*p15)*qmd4x
        r8( 9)= r8( 9)+(f8(3)*xiq4+f7(1,3)*xiq4+f7(3,3)*xiq2*ten+f6(1,7)*xiq2*ten &
&                    +f6(3,7)*p15+f5(1,11)*p15)*qmd4xy
        r8(10)= r8(10)+(f8(4)*xiq4+f7(2,3)*xiq4*three+f7(4,3)*xiq2*ten+f6(2,7)*xiq2*p30 &
&                    +f6(4,7)*p15+f5(2,11)*p45)*qmd4x
        r8(11)= r8(11)+(f8(1)*x4y4+f7(1,3)*x4y2*six+f7(1,3)*x2y4*six+f6(1,7)*xiq4*three &
&                    +f6(1,7)*xyiq2*p36+f6(1,7)*yiq4*three+f5(1,11)*xiq2*p18 &
&                    +f5(1,11)*yiq2*p18+f4(1,16)*nine)*qmd4
        r8(12)= r8(12)+(f8(2)*x4y2+f7(2,3)*xiq4*three+f7(2,3)*xyiq2*six+f6(2,7)*xiq2*p18 &
&                    +f6(2,7)*yiq2*three+f5(2,11)*nine)*qmd4y
        r8(13)= r8(13)+(f8(3)*x4y2+f7(1,3)*x4y2+f7(3,3)*xiq4+f7(3,3)*xyiq2*six+f6(1,7)*xiq4 &
&                    +f6(1,7)*xyiq2*six+f6(3,7)*xiq2*six+f6(3,7)*yiq2*three+f5(1,11)*xiq2*six &
&                    +f5(1,11)*yiq2*three+f5(3,11)*three+f4(1,16)*three)*qmd4
        r8(14)= r8(14)+(f8(4)*xiq4+f7(2,3)*xiq4*three+f7(4,3)*xiq2*six+f6(2,7)*xiq2*p18 &
&                    +f6(4,7)*three+f5(2,11)*nine)*qmd4y
        r8(15)= r8(15)+(f8(5)*xiq4+f7(3,3)*xiq4*six+f7(5,3)*xiq2*six+f6(1,7)*xiq4*three &
&                    +f6(3,7)*xiq2*p36+f6(5,7)*three+f5(1,11)*xiq2*p18+f5(3,11)*p18 &
&                    +f4(1,16)*nine)*qmd4
        r8(16)= r8(16)+(f8(1)*x2y4+f7(1,3)*xyiq2*ten+f7(1,3)*yiq4*three+f6(1,7)*xiq2*p15 &
&                    +f6(1,7)*yiq2*p30+f5(1,11)*p45)*qmd4xy
        r8(17)= r8(17)+(f8(2)*x2y4+f7(2,3)*xyiq2*six+f7(2,3)*yiq4*three+f6(2,7)*xiq2*three &
&                    +f6(2,7)*yiq2*p18+f5(2,11)*nine)*qmd4x
        r8(18)= r8(18)+(f8(3)*xyiq2+f7(1,3)*xyiq2+f7(3,3)*xiq2*three+f7(3,3)*yiq2*three &
&                    +f6(1,7)*xiq2*three+f6(1,7)*yiq2*three+f6(3,7)*nine+f5(1,11)*nine)*qmd4xy
        r8(19)= r8(19)+(f8(4)*xyiq2+f7(2,3)*xyiq2*three+f7(4,3)*xiq2+f7(4,3)*yiq2*three &
&                    +f6(2,7)*xiq2*three+f6(2,7)*yiq2*nine+f6(4,7)*three+f5(2,11)*nine)*qmd4x
        r8(20)= r8(20)+(f8(5)*xiq2+f7(3,3)*xiq2*six+f7(5,3)*three+f6(1,7)*xiq2*three &
&                    +f6(3,7)*p18+f5(1,11)*nine)*qmd4xy
        r8(21)= r8(21)+(f8(6)*xiq2+f7(4,3)*xiq2*ten+f7(6,3)*three+f6(2,7)*xiq2*p15 &
&                    +f6(4,7)*p30+f5(2,11)*p45)*qmd4x
        r8(22)= r8(22)+(f8(1)*x2y6+f7(1,3)*x2y4*p15+f7(1,3)*yiq6+f6(1,7)*xyiq2*p45 &
&                    +f6(1,7)*yiq4*p15+f5(1,11)*xiq2*p15+f5(1,11)*yiq2*p45+f4(1,16)*p15)*qmd4
        r8(23)= r8(23)+(f8(2)*x2y4+f7(2,3)*xyiq2*ten+f7(2,3)*yiq4+f6(2,7)*xiq2*p15 &
&                    +f6(2,7)*yiq2*ten+f5(2,11)*p15)*qmd4y
        r8(24)= r8(24)+(f8(3)*x2y4+f7(1,3)*x2y4+f7(3,3)*xyiq2*six+f7(3,3)*yiq4 &
&                    +f6(1,7)*xyiq2*six+f6(1,7)*yiq4+f6(3,7)*xiq2*three+f6(3,7)*yiq2*six &
&                    +f5(1,11)*xiq2*three+f5(1,11)*yiq2*six+f5(3,11)*three+f4(1,16)*three)*qmd4
        r8(25)= r8(25)+(f8(4)*xyiq2+f7(2,3)*xyiq2*three+f7(4,3)*xiq2*three+f7(4,3)*yiq2 &
&                    +f6(2,7)*xiq2*nine+f6(2,7)*yiq2*three+f6(4,7)*three+f5(2,11)*nine)*qmd4y
        r8(26)= r8(26)+(f8(5)*xyiq2+f7(3,3)*xyiq2*six+f7(5,3)*xiq2+f7(5,3)*yiq2 &
&                    +f6(1,7)*xyiq2*three+f6(3,7)*xiq2*six+f6(3,7)*yiq2*six+f6(5,7) &
&                    +f5(1,11)*xiq2*three+f5(1,11)*yiq2*three+f5(3,11)*six+f4(1,16)*three)*qmd4
        r8(27)= r8(27)+(f8(6)*xiq2+f7(4,3)*xiq2*ten+f7(6,3)+f6(2,7)*xiq2*p15+f6(4,7)*ten &
&                    +f5(2,11)*p15)*qmd4y
        r8(28)= r8(28)+(f8(7)*xiq2+f7(5,3)*xiq2*p15+f7(7,3)+f6(3,7)*xiq2*p45+f6(5,7)*p15 &
&                    +f5(1,11)*xiq2*p15+f5(3,11)*p45+f4(1,16)*p15)*qmd4
        r8(29)= r8(29)+(f8(1)*yiq6+f7(1,3)*yiq4*p21+f6(1,7)*yiq2*p105+f5(1,11)*p105)*qmd4xy
        r8(30)= r8(30)+(f8(2)*yiq6+f7(2,3)*yiq4*p15+f6(2,7)*yiq2*p45+f5(2,11)*p15)*qmd4x
        r8(31)= r8(31)+(f8(3)*yiq4+f7(1,3)*yiq4+f7(3,3)*yiq2*ten+f6(1,7)*yiq2*ten &
&                    +f6(3,7)*p15+f5(1,11)*p15)*qmd4xy
        r8(32)= r8(32)+(f8(4)*yiq4+f7(2,3)*yiq4*three+f7(4,3)*yiq2*six+f6(2,7)*yiq2*p18 &
&                    +f6(4,7)*three+f5(2,11)*nine)*qmd4x
        r8(33)= r8(33)+(f8(5)*yiq2+f7(3,3)*yiq2*six+f7(5,3)*three+f6(1,7)*yiq2*three &
&                    +f6(3,7)*p18+f5(1,11)*nine)*qmd4xy
        r8(34)= r8(34)+(f8(6)*yiq2+f7(4,3)*yiq2*ten+f7(6,3)+f6(2,7)*yiq2*p15+f6(4,7)*ten &
&                    +f5(2,11)*p15)*qmd4x
        r8(35)= r8(35)+(f8(7)+f7(5,3)*p15+f6(3,7)*p45+f5(1,11)*p15)*qmd4xy
        r8(36)= r8(36)+(f8(8)+f7(6,3)*p21+f6(4,7)*p105+f5(2,11)*p105)*qmd4x
        r8(37)= r8(37)+(f8(1)*yiq8+f7(1,3)*yiq6*p28+f6(1,7)*yiq4*p210+f5(1,11)*yiq2*p420 &
&                    +f4(1,16)*p105)*qmd4
        r8(38)= r8(38)+(f8(2)*yiq6+f7(2,3)*yiq4*p21+f6(2,7)*yiq2*p105+f5(2,11)*p105)*qmd4y
        r8(39)= r8(39)+(f8(3)*yiq6+f7(1,3)*yiq6+f7(3,3)*yiq4*p15+f6(1,7)*yiq4*p15 &
&                    +f6(3,7)*yiq2*p45+f5(1,11)*yiq2*p45+f5(3,11)*p15+f4(1,16)*p15)*qmd4
        r8(40)= r8(40)+(f8(4)*yiq4+f7(2,3)*yiq4*three+f7(4,3)*yiq2*ten+f6(2,7)*yiq2*p30 &
&                    +f6(4,7)*p15+f5(2,11)*p45)*qmd4y
        r8(41)= r8(41)+(f8(5)*yiq4+f7(3,3)*yiq4*six+f7(5,3)*yiq2*six+f6(1,7)*yiq4*three &
&                    +f6(3,7)*yiq2*p36+f6(5,7)*three+f5(1,11)*yiq2*p18+f5(3,11)*p18 &
&                    +f4(1,16)*nine)*qmd4
        r8(42)= r8(42)+(f8(6)*yiq2+f7(4,3)*yiq2*ten+f7(6,3)*three+f6(2,7)*yiq2*p15 &
&                    +f6(4,7)*p30+f5(2,11)*p45)*qmd4y
        r8(43)= r8(43)+(f8(7)*yiq2+f7(5,3)*yiq2*p15+f7(7,3)+f6(3,7)*yiq2*p45+f6(5,7)*p15 &
&                    +f5(1,11)*yiq2*p15+f5(3,11)*p45+f4(1,16)*p15)*qmd4
        r8(44)= r8(44)+(f8(8)+f7(6,3)*p21+f6(4,7)*p105+f5(2,11)*p105)*qmd4y
        r8(45)= r8(45)+(f8(9)+f7(7,3)*p28+f6(5,7)*p210+f5(3,11)*p420+f4(1,16)*p105)*qmd4
      enddo
!
      qx= xzkl(1)
      qz= xzkl(2)
      call int2dddd1(eri,r0,r1,r2,r3,r4,r5,r6,r7,r8,qx,qz)
      call int2dddd2(eri,r0,r1,r2,r3,r4,r5,r6,r7,r8,qx,qz)
      call int2dddd3(eri,r0,r1,r2,r3,r4,r5,r6,r7,r8,qx,qz)
      call int2dddd4(eri,r0,r1,r2,r3,r4,r5,r6,r7,r8,qx,qz)
      call int2dddd5(eri,r0,r1,r2,r3,r4,r5,r6,r7,r8,qx,qz)
!
      rot2(1,1)= rot(1,1)*rot(1,1)
      rot2(2,1)= rot(2,1)*rot(2,1)
      rot2(3,1)= rot(3,1)*rot(3,1)
      rot2(4,1)= rot(1,1)*rot(2,1)*two
      rot2(5,1)= rot(1,1)*rot(3,1)*two
      rot2(6,1)= rot(2,1)*rot(3,1)*two
      do l= 2,3
        rot2(1,l)= rot(1,1)*rot(1,l)*sqrt3
        rot2(2,l)= rot(2,1)*rot(2,l)*sqrt3
        rot2(3,l)= rot(3,1)*rot(3,l)*sqrt3
        rot2(4,l)=(rot(1,1)*rot(2,l)+rot(2,1)*rot(1,l))*sqrt3
        rot2(5,l)=(rot(1,1)*rot(3,l)+rot(3,1)*rot(1,l))*sqrt3
        rot2(6,l)=(rot(2,1)*rot(3,l)+rot(3,1)*rot(2,l))*sqrt3
      enddo
      do l= 2,3
        rot2(1,l*2)= rot(1,l)*rot(1,l)
        rot2(2,l*2)= rot(2,l)*rot(2,l)
        rot2(3,l*2)= rot(3,l)*rot(3,l)
        rot2(4,l*2)= rot(1,l)*rot(2,l)*two
        rot2(5,l*2)= rot(1,l)*rot(3,l)*two
        rot2(6,l*2)= rot(2,l)*rot(3,l)*two
      enddo
      rot2(1,5)= rot(1,2)*rot(1,3)*sqrt3
      rot2(2,5)= rot(2,2)*rot(2,3)*sqrt3
      rot2(3,5)= rot(3,2)*rot(3,3)*sqrt3
      rot2(4,5)=(rot(1,2)*rot(2,3)+rot(2,2)*rot(1,3))*sqrt3
      rot2(5,5)=(rot(1,2)*rot(3,3)+rot(3,2)*rot(1,3))*sqrt3
      rot2(6,5)=(rot(2,2)*rot(3,3)+rot(3,2)*rot(2,3))*sqrt3
!
      do l= 1,6
        rot3(l,1)= rot2(l,2)
        rot3(l,2)= rot2(l,5)
        rot3(l,3)= rot2(l,6)-(rot2(l,1)+rot2(l,4))*half
        rot3(l,4)= rot2(l,3)
        rot3(l,5)=(rot2(l,1)-rot2(l,4))*sqrt3h
      enddo
!
      if(nbfijkl(1) == 6)then
        do j= 1,6
          do k= 1,6
            do l= 1,6
              do i= 1,6
                work(i)= eri(l,k,j,i)
              enddo
              do i= 1,6
                eri(l,k,j,i)= work(1)*rot2(1,i)+work(2)*rot2(2,i)+work(3)*rot2(3,i) &
&                            +work(4)*rot2(4,i)+work(5)*rot2(5,i)+work(6)*rot2(6,i)
              enddo
            enddo
          enddo
        enddo
      else
        do j= 1,6
          do k= 1,6
            do l= 1,6
              do i= 1,6
                work(i)= eri(l,k,j,i)
              enddo
              do i= 1,5
                eri(l,k,j,i)= work(1)*rot3(1,i)+work(2)*rot3(2,i)+work(3)*rot3(3,i) &
&                            +work(4)*rot3(4,i)+work(5)*rot3(5,i)+work(6)*rot3(6,i)
              enddo
            enddo
          enddo
        enddo
      endif
!
      if(nbfijkl(2) == 6)then
        do i= 1,nbfijkl(1)
          do k= 1,6
            do l= 1,6
              do j= 1,6
                work(j)= eri(l,k,j,i)
              enddo
              do j= 1,6
                eri(l,k,j,i)= work(1)*rot2(1,j)+work(2)*rot2(2,j)+work(3)*rot2(3,j) &
&                            +work(4)*rot2(4,j)+work(5)*rot2(5,j)+work(6)*rot2(6,j)
              enddo
            enddo
          enddo
        enddo
      else
        do i= 1,nbfijkl(1)
          do k= 1,6
            do l= 1,6
              do j= 1,6
                work(j)= eri(l,k,j,i)
              enddo
              do j= 1,5
                eri(l,k,j,i)= work(1)*rot3(1,j)+work(2)*rot3(2,j)+work(3)*rot3(3,j) &
&                            +work(4)*rot3(4,j)+work(5)*rot3(5,j)+work(6)*rot3(6,j)
              enddo
            enddo
          enddo
        enddo
      endif
!
      if(nbfijkl(3) == 6)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do l= 1,6
              do k= 1,6
                work(k)= eri(l,k,j,i)
              enddo
              do k= 1,6
                eri(l,k,j,i)= work(1)*rot2(1,k)+work(2)*rot2(2,k)+work(3)*rot2(3,k) &
&                            +work(4)*rot2(4,k)+work(5)*rot2(5,k)+work(6)*rot2(6,k)
              enddo
            enddo
          enddo
        enddo
      else
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do l= 1,6
              do k= 1,6
                work(k)= eri(l,k,j,i)
              enddo
              do k= 1,5
                eri(l,k,j,i)= work(1)*rot3(1,k)+work(2)*rot3(2,k)+work(3)*rot3(3,k) &
&                            +work(4)*rot3(4,k)+work(5)*rot3(5,k)+work(6)*rot3(6,k)
              enddo
            enddo
          enddo
        enddo
      endif
!
      if(nbfijkl(4) == 6)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do k= 1,nbfijkl(3)
              do l= 1,6
                phmdint(l,k,j,i)= &
&                   eri(1,k,j,i)*rot2(1,l)+eri(2,k,j,i)*rot2(2,l)+eri(3,k,j,i)*rot2(3,l) &
&                  +eri(4,k,j,i)*rot2(4,l)+eri(5,k,j,i)*rot2(5,l)+eri(6,k,j,i)*rot2(6,l) 
              enddo
            enddo
          enddo
        enddo
      else
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do k= 1,nbfijkl(3)
              do l= 1,5
                phmdint(l,k,j,i)= &
&                   eri(1,k,j,i)*rot3(1,l)+eri(2,k,j,i)*rot3(2,l)+eri(3,k,j,i)*rot3(3,l) &
&                  +eri(4,k,j,i)*rot3(4,l)+eri(5,k,j,i)*rot3(5,l)+eri(6,k,j,i)*rot3(6,l)
              enddo
            enddo
          enddo
        enddo
      endif
!
      return
end
!-------------------------------------------------------------
  subroutine int2dddd1(eri,r0,r1,r2,r3,r4,r5,r6,r7,r8,qx,qz)
!-------------------------------------------------------------
!
      implicit none
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00, four=4.0D+00, five=5.0D+00
      real(8),parameter :: six=6.0D+00, seven=7.0D+00, eight=8.0D+00, nine=9.0D+00, ten=1.0D+01
      real(8),parameter :: p11=1.1D+01, p12=1.2D+01, p13=1.3D+01, p15=1.5D+01, p16=1.6D+01
      real(8),parameter :: p18=1.8D+01, p20=2.0D+01, p21=2.1D+01, p24=2.4D+01, p28=2.8D+01
      real(8),parameter :: p30=3.0D+01, p36=3.6D+01, p45=4.5D+01, p105=1.05D+2, p210=2.1D+02
      real(8),parameter :: p420=4.2D+02
      real(8),intent(in) :: r0(25), r1(3,40), r2(6,56), r3(10,52), r4(15,42), r5(21,24)
      real(8),intent(in) :: r6(28,12), r7(36,4), r8(45), qx, qz
      real(8),intent(inout) :: eri(6,6,6,6)
      real(8) :: r400, r310, r301, r220, r211, r202, r130, r121, r112, r103, r040, r031
      real(8) :: r022, r013, r004, rxyz(20), xx, xz, zz, xxx, xxz, xzz, zzz
      real(8) :: xxxx, xxxz, xxzz, xzzz, zzzz
!
      xx= qx*qx
      xz= qx*qz
      zz= qz*qz
      xxx= qx*qx*qx
      xxz= qx*qx*qz
      xzz= qx*qz*qz
      zzz= qz*qz*qz
      xxxx= xx*xx
      xxxz= xx*xz
      xxzz= xx*zz
      xzzz= xz*zz
      zzzz= zz*zz
!
      r400= r8(1)+r6(1,4)*six+r6(1,9)*six+r4(1,5)*three+r4(1,26)*p36+r4(1,38)*three &
&          +r2(1,17)*p18+r2(1,52)*p18+r0(21)*nine
      r310= r8(2)+r6(2,4)*six+r6(2,9)*three+r4(2,5)*three+r4(2,26)*p18+r2(4,17)*nine
      r301= r8(3)+r6(3,4)*six+r6(3,9)*three+r4(3,5)*three+r4(3,26)*p18+r2(5,17)*nine
      r220= r8(4)+r6(4,4)*six+r6(1,9)+r6(4,9)+r4(4,5)*three+r4(1,26)*six+r4(4,26)*six &
&          +r4(1,38)+r2(1,17)*three+r2(2,17)*three+r2(1,52)*six+r0(21)*three
      r211= r8(5)+r6(5,4)*six+r6(5,9)+r4(5,5)*three+r4(5,26)*six+r2(6,17)*three
      r202= r8(6)+r6(6,4)*six+r6(1,9)+r6(6,9)+r4(6,5)*three+r4(1,26)*six+r4(6,26)*six &
&          +r4(1,38)+r2(1,17)*three+r2(3,17)*three+r2(1,52)*six+r0(21)*three
      r130= r8(7)+r6(7,4)*six+r6(2,9)*three+r4(7,5)*three+r4(2,26)*p18+r2(4,17)*nine
      r121= r8(8)+r6(8,4)*six+r6(3,9)+r4(8,5)*three+r4(3,26)*six+r2(5,17)*three
      r112= r8(9)+r6(9,4)*six+r6(2,9)+r4(9,5)*three+r4(2,26)*six+r2(4,17)*three
      r103= r8(10)+r6(10,4)*six+r6(3,9)*three+r4(10,5)*three+r4(3,26)*p18+r2(5,17)*nine
      r040= r8(11)+r6(11,4)*six+r6(4,9)*six+r4(11,5)*three+r4(4,26)*p36+r4(1,38)*three &
&          +r2(2,17)*p18+r2(1,52)*p18+r0(21)*nine
      r031= r8(12)+r6(12,4)*six+r6(5,9)*three+r4(12,5)*three+r4(5,26)*p18+r2(6,17)*nine
      r022= r8(13)+r6(13,4)*six+r6(4,9)+r6(6,9)+r4(13,5)*three+r4(4,26)*six+r4(6,26)*six &
&          +r4(1,38)+r2(2,17)*three+r2(3,17)*three+r2(1,52)*six+r0(21)*three
      r013= r8(14)+r6(14,4)*six+r6(5,9)*three+r4(14,5)*three+r4(5,26)*p18+r2(6,17)*nine
      r004= r8(15)+r6(15,4)*six+r6(6,9)*six+r4(15,5)*three+r4(6,26)*p36+r4(1,38)*three &
&          +r2(3,17)*p18+r2(1,52)*p18+r0(21)*nine
      rxyz(1)=+r4(1,42)+r2(1,56)*six+r0(25)*three
      rxyz(2)=+r5(2,24)+r3(2,42)*six+r1(2,20)*three
      rxyz(3)=+r5(2,23)+r3(2,41)*six+r1(2,19)*three
      rxyz(4)=+r6(4,10)+r4(4,27)*six+r4(1,39)+r2(2,18)*three+r2(1,53)*six+r0(22)*three
      rxyz(5)=+r6(4,12)+r4(4,29)*six+r4(1,41)+r2(2,20)*three+r2(1,55)*six+r0(24)*three
      rxyz(6)=+r6(4,11)+r4(4,28)*six+r4(1,40)+r2(2,19)*three+r2(1,54)*six+r0(23)*three
      rxyz(7)=+r7(7,3)+r5(7,11)*six+r5(2,21)*three+r3(7,9)*three+r3(2,39)*p18 &
&             +r1(2,17)*nine
      rxyz(8)=+r7(7,4)+r5(7,12)*six+r5(2,22)*three+r3(7,10)*three+r3(2,40)*p18 &
&             +r1(2,18)*nine
      rxyz(9)=+r7(5,3)+r7(5,4)+r5(5,11)*six+r5(5,12)*six+r3(5,9)*three+r3(5,10)*three
      rxyz(10)=+r6(5,11)+r4(5,28)*six+r2(6,19)*three
      rxyz(11)=+r6(2,11)+r4(2,28)*six+r2(4,19)*three
      rxyz(12)=+r7(8,3)+r5(8,11)*six+r5(3,21)+r3(8,9)*three+r3(3,39)*six+r1(3,17)*three
      rxyz(13)=+r7(8,4)+r5(8,12)*six+r5(3,22)+r3(8,10)*three+r3(3,40)*six+r1(3,18)*three
      rxyz(14)=+r7(4,4)+r5(4,12)*six+r5(1,22)+r3(4,10)*three+r3(1,40)*six+r1(1,18)*three
      rxyz(15)=+r7(4,3)+r5(4,11)*six+r5(1,21)+r3(4,9)*three+r3(1,39)*six+r1(1,17)*three
      rxyz(16)=+r7(2,4)+r5(2,12)*six+r5(2,22)+r3(2,10)*three+r3(2,40)*six+r1(2,18)*three
      rxyz(17)=+r7(2,3)+r5(2,11)*six+r5(2,21)+r3(2,9)*three+r3(2,39)*six+r1(2,17)*three
      rxyz(18)=+r7(9,3)+r5(9,11)*six+r5(2,21)+r3(9,9)*three+r3(2,39)*six+r1(2,17)*three
      rxyz(19)=+r7(9,4)+r5(9,12)*six+r5(2,22)+r3(9,10)*three+r3(2,40)*six+r1(2,18)*three
      rxyz(20)=+r6(3,11)*four+r4(3,28)*p24+r2(5,19)*p12
      eri(1,1,1,1)=r400+(+r7(1,3)*two+r7(1,4)*two+r5(1,11)*p12+r5(1,12)*p12+r5(1,21)*six &
&                  +r5(1,22)*six+r3(1,9)*six+r3(1,10)*six+r3(1,39)*p36+r3(1,40)*p36 &
&                  +r1(1,17)*p18+r1(1,18)*p18)*qx+(+r6(1,10)+r6(1,11)*four+r6(1,12) &
&                  +r4(1,27)*six+r4(1,28)*p24+r4(1,29)*six+r4(1,39)+r4(1,40)*four+r4(1,41) &
&                  +r2(1,18)*three+r2(1,19)*p12+r2(1,20)*three+r2(1,53)*six+r2(1,54)*p24 &
&                  +r2(1,55)*six+r0(22)*three+r0(23)*p12+r0(24)*three)*xx+(+r5(1,23)*two &
&                  +r5(1,24)*two+r3(1,41)*p12+r3(1,42)*p12+r1(1,19)*six+r1(1,20)*six)*xxx &
&                  +rxyz(1)*xxxx
      eri(2,1,1,1)=r220+(+r7(4,4)*two+r5(4,12)*p12+r5(1,22)*two+r3(4,10)*six+r3(1,40)*p12 &
&                  +r1(1,18)*six)*qx+rxyz(5)*xx
      eri(3,1,1,1)=r202+(+r7(6,4)*two+r5(6,12)*p12+r5(1,22)*two+r3(6,10)*six+r3(1,40)*p12 &
&                  +r1(1,18)*six)*qx+(+r7(3,3)*two+r5(3,11)*p12+r5(3,21)*two+r3(3,9)*six &
&                  +r3(3,39)*p12+r1(3,17)*six)*qz+(+r6(6,12)+r4(6,29)*six+r4(1,41) &
&                  +r2(3,20)*three+r2(1,55)*six+r0(24)*three)*xx+rxyz(20)*xz+(+r6(1,10) &
&                  +r4(1,27)*six+r4(1,39)+r2(1,18)*three+r2(1,53)*six+r0(22)*three)*zz+( &
&                  +r5(3,24)*two+r3(3,42)*p12+r1(3,20)*six)*xxz+(+r5(1,23)*two+r3(1,41)*p12 &
&                  +r1(1,19)*six)*xzz+rxyz(1)*xxzz
      eri(4,1,1,1)=r310+(+r7(2,3)+r7(2,4)*two+r5(2,11)*six+r5(2,12)*p12+r5(2,21) &
&                  +r5(2,22)*two+r3(2,9)*three+r3(2,10)*six+r3(2,39)*six+r3(2,40)*p12 &
&                  +r1(2,17)*three+r1(2,18)*six)*qx+(+r6(2,11)*two+r6(2,12)+r4(2,28)*p12 &
&                  +r4(2,29)*six+r2(4,19)*six+r2(4,20)*three)*xx+rxyz(2)*xxx
      eri(5,1,1,1)=r301+(+r7(3,3)+r7(3,4)*two+r5(3,11)*six+r5(3,12)*p12+r5(3,21) &
&                  +r5(3,22)*two+r3(3,9)*three+r3(3,10)*six+r3(3,39)*six+r3(3,40)*p12 &
&                  +r1(3,17)*three+r1(3,18)*six)*qx+(+r7(1,3)+r5(1,11)*six+r5(1,21)*three &
&                  +r3(1,9)*three+r3(1,39)*p18+r1(1,17)*nine)*qz+(+r6(3,11)*two+r6(3,12) &
&                  +r4(3,28)*p12+r4(3,29)*six+r2(5,19)*six+r2(5,20)*three)*xx+(+r6(1,10) &
&                  +r6(1,11)*two+r4(1,27)*six+r4(1,28)*p12+r4(1,39)+r4(1,40)*two &
&                  +r2(1,18)*three+r2(1,19)*six+r2(1,53)*six+r2(1,54)*p12+r0(22)*three &
&                  +r0(23)*six)*xz+(+r5(3,24)+r3(3,42)*six+r1(3,20)*three)*xxx+(+r5(1,23)*two &
&                  +r5(1,24)+r3(1,41)*p12+r3(1,42)*six+r1(1,19)*six+r1(1,20)*three)*xxz &
&                  +rxyz(1)*xxxz
      eri(6,1,1,1)=r211+(+r7(5,4)*two+r5(5,12)*p12+r3(5,10)*six)*qx+rxyz(17)*qz+( &
&                  +r6(5,12)+r4(5,29)*six+r2(6,20)*three)*xx+(+r6(2,11)*two+r4(2,28)*p12 &
&                  +r2(4,19)*six)*xz+rxyz(2)*xxz
      eri(1,2,1,1)=r220+(+r7(4,3)*two+r5(4,11)*p12+r5(1,21)*two+r3(4,9)*six+r3(1,39)*p12 &
&                  +r1(1,17)*six)*qx+rxyz(4)*xx
      eri(2,2,1,1)=r040
      eri(3,2,1,1)=r022+(+r7(8,3)*two+r5(8,11)*p12+r5(3,21)*two+r3(8,9)*six+r3(3,39)*p12 &
&                  +r1(3,17)*six)*qz+rxyz(4)*zz
      eri(4,2,1,1)=r130+rxyz(7)*qx
      eri(5,2,1,1)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,1,1)=r031+rxyz(7)*qz
      eri(1,3,1,1)=r202+(+r7(6,3)*two+r5(6,11)*p12+r5(1,21)*two+r3(6,9)*six+r3(1,39)*p12 &
&                  +r1(1,17)*six)*qx+(+r7(3,4)*two+r5(3,12)*p12+r5(3,22)*two+r3(3,10)*six &
&                  +r3(3,40)*p12+r1(3,18)*six)*qz+(+r6(6,10)+r4(6,27)*six+r4(1,39) &
&                  +r2(3,18)*three+r2(1,53)*six+r0(22)*three)*xx+rxyz(20)*xz+(+r6(1,12) &
&                  +r4(1,29)*six+r4(1,41)+r2(1,20)*three+r2(1,55)*six+r0(24)*three)*zz+( &
&                  +r5(3,23)*two+r3(3,41)*p12+r1(3,19)*six)*xxz+(+r5(1,24)*two+r3(1,42)*p12 &
&                  +r1(1,20)*six)*xzz+rxyz(1)*xxzz
      eri(2,3,1,1)=r022+(+r7(8,4)*two+r5(8,12)*p12+r5(3,22)*two+r3(8,10)*six+r3(3,40)*p12 &
&                  +r1(3,18)*six)*qz+rxyz(5)*zz
      eri(3,3,1,1)=r004+(+r7(10,3)*two+r7(10,4)*two+r5(10,11)*p12+r5(10,12)*p12 &
&                  +r5(3,21)*six+r5(3,22)*six+r3(10,9)*six+r3(10,10)*six+r3(3,39)*p36 &
&                  +r3(3,40)*p36+r1(3,17)*p18+r1(3,18)*p18)*qz+(+r6(6,10)+r6(6,11)*four &
&                  +r6(6,12)+r4(6,27)*six+r4(6,28)*p24+r4(6,29)*six+r4(1,39)+r4(1,40)*four &
&                  +r4(1,41)+r2(3,18)*three+r2(3,19)*p12+r2(3,20)*three+r2(1,53)*six &
&                  +r2(1,54)*p24+r2(1,55)*six+r0(22)*three+r0(23)*p12+r0(24)*three)*zz+( &
&                  +r5(3,23)*two+r5(3,24)*two+r3(3,41)*p12+r3(3,42)*p12+r1(3,19)*six &
&                  +r1(3,20)*six)*zzz+rxyz(1)*zzzz
      eri(4,3,1,1)=r112+rxyz(18)*qx+(+r7(5,4)*two+r5(5,12)*p12+r3(5,10)*six)*qz+( &
&                  +r6(5,11)*two+r4(5,28)*p12+r2(6,19)*six)*xz+(+r6(2,12)+r4(2,29)*six &
&                  +r2(4,20)*three)*zz+rxyz(2)*xzz
      eri(5,3,1,1)=r103+(+r7(10,3)+r5(10,11)*six+r5(3,21)*three+r3(10,9)*three &
&                  +r3(3,39)*p18+r1(3,17)*nine)*qx+(+r7(6,3)+r7(6,4)*two+r5(6,11)*six &
&                  +r5(6,12)*p12+r5(1,21)+r5(1,22)*two+r3(6,9)*three+r3(6,10)*six &
&                  +r3(1,39)*six+r3(1,40)*p12+r1(1,17)*three+r1(1,18)*six)*qz+(+r6(6,10) &
&                  +r6(6,11)*two+r4(6,27)*six+r4(6,28)*p12+r4(1,39)+r4(1,40)*two &
&                  +r2(3,18)*three+r2(3,19)*six+r2(1,53)*six+r2(1,54)*p12+r0(22)*three &
&                  +r0(23)*six)*xz+(+r6(3,11)*two+r6(3,12)+r4(3,28)*p12+r4(3,29)*six &
&                  +r2(5,19)*six+r2(5,20)*three)*zz+(+r5(3,23)*two+r5(3,24)+r3(3,41)*p12 &
&                  +r3(3,42)*six+r1(3,19)*six+r1(3,20)*three)*xzz+(+r5(1,24)+r3(1,42)*six &
&                  +r1(1,20)*three)*zzz+rxyz(1)*xzzz
      eri(6,3,1,1)=r013+(+r7(9,3)+r7(9,4)*two+r5(9,11)*six+r5(9,12)*p12+r5(2,21) &
&                  +r5(2,22)*two+r3(9,9)*three+r3(9,10)*six+r3(2,39)*six+r3(2,40)*p12 &
&                  +r1(2,17)*three+r1(2,18)*six)*qz+(+r6(5,11)*two+r6(5,12)+r4(5,28)*p12 &
&                  +r4(5,29)*six+r2(6,19)*six+r2(6,20)*three)*zz+rxyz(2)*zzz
      eri(1,4,1,1)=r310+(+r7(2,3)*two+r7(2,4)+r5(2,11)*p12+r5(2,12)*six+r5(2,21)*two &
&                  +r5(2,22)+r3(2,9)*six+r3(2,10)*three+r3(2,39)*p12+r3(2,40)*six &
&                  +r1(2,17)*six+r1(2,18)*three)*qx+(+r6(2,10)+r6(2,11)*two+r4(2,27)*six &
&                  +r4(2,28)*p12+r2(4,18)*three+r2(4,19)*six)*xx+rxyz(3)*xxx
      eri(2,4,1,1)=r130+rxyz(8)*qx
      eri(3,4,1,1)=r112+rxyz(19)*qx+(+r7(5,3)*two+r5(5,11)*p12+r3(5,9)*six)*qz+( &
&                  +r6(5,11)*two+r4(5,28)*p12+r2(6,19)*six)*xz+(+r6(2,10)+r4(2,27)*six &
&                  +r2(4,18)*three)*zz+rxyz(3)*xzz
      eri(4,4,1,1)=r220+(+r7(4,3)+r7(4,4)+r5(4,11)*six+r5(4,12)*six+r5(1,21)+r5(1,22) &
&                  +r3(4,9)*three+r3(4,10)*three+r3(1,39)*six+r3(1,40)*six+r1(1,17)*three &
&                  +r1(1,18)*three)*qx+rxyz(6)*xx
      eri(5,4,1,1)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(2,10)+r6(2,11) &
&                  +r4(2,27)*six+r4(2,28)*six+r2(4,18)*three+r2(4,19)*three)*xz+rxyz(3)*xxz
      eri(6,4,1,1)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,1,1)=r301+(+r7(3,3)*two+r7(3,4)+r5(3,11)*p12+r5(3,12)*six+r5(3,21)*two &
&                  +r5(3,22)+r3(3,9)*six+r3(3,10)*three+r3(3,39)*p12+r3(3,40)*six &
&                  +r1(3,17)*six+r1(3,18)*three)*qx+(+r7(1,4)+r5(1,12)*six+r5(1,22)*three &
&                  +r3(1,10)*three+r3(1,40)*p18+r1(1,18)*nine)*qz+(+r6(3,10)+r6(3,11)*two &
&                  +r4(3,27)*six+r4(3,28)*p12+r2(5,18)*three+r2(5,19)*six)*xx+(+r6(1,11)*two &
&                  +r6(1,12)+r4(1,28)*p12+r4(1,29)*six+r4(1,40)*two+r4(1,41)+r2(1,19)*six &
&                  +r2(1,20)*three+r2(1,54)*p12+r2(1,55)*six+r0(23)*six+r0(24)*three)*xz+( &
&                  +r5(3,23)+r3(3,41)*six+r1(3,19)*three)*xxx+(+r5(1,23)+r5(1,24)*two &
&                  +r3(1,41)*six+r3(1,42)*p12+r1(1,19)*three+r1(1,20)*six)*xxz+rxyz(1)*xxxz
      eri(2,5,1,1)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,1,1)=r103+(+r7(10,4)+r5(10,12)*six+r5(3,22)*three+r3(10,10)*three &
&                  +r3(3,40)*p18+r1(3,18)*nine)*qx+(+r7(6,3)*two+r7(6,4)+r5(6,11)*p12 &
&                  +r5(6,12)*six+r5(1,21)*two+r5(1,22)+r3(6,9)*six+r3(6,10)*three &
&                  +r3(1,39)*p12+r3(1,40)*six+r1(1,17)*six+r1(1,18)*three)*qz+(+r6(6,11)*two &
&                  +r6(6,12)+r4(6,28)*p12+r4(6,29)*six+r4(1,40)*two+r4(1,41)+r2(3,19)*six &
&                  +r2(3,20)*three+r2(1,54)*p12+r2(1,55)*six+r0(23)*six+r0(24)*three)*xz+( &
&                  +r6(3,10)+r6(3,11)*two+r4(3,27)*six+r4(3,28)*p12+r2(5,18)*three &
&                  +r2(5,19)*six)*zz+(+r5(3,23)+r5(3,24)*two+r3(3,41)*six+r3(3,42)*p12 &
&                  +r1(3,19)*three+r1(3,20)*six)*xzz+(+r5(1,23)+r3(1,41)*six+r1(1,19)*three) &
&                  *zzz+rxyz(1)*xzzz
      eri(4,5,1,1)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(2,11)+r6(2,12) &
&                  +r4(2,28)*six+r4(2,29)*six+r2(4,19)*three+r2(4,20)*three)*xz+rxyz(2)*xxz
      eri(5,5,1,1)=r202+(+r7(6,3)+r7(6,4)+r5(6,11)*six+r5(6,12)*six+r5(1,21)+r5(1,22) &
&                  +r3(6,9)*three+r3(6,10)*three+r3(1,39)*six+r3(1,40)*six+r1(1,17)*three &
&                  +r1(1,18)*three)*qx+(+r7(3,3)+r7(3,4)+r5(3,11)*six+r5(3,12)*six+r5(3,21) &
&                  +r5(3,22)+r3(3,9)*three+r3(3,10)*three+r3(3,39)*six+r3(3,40)*six &
&                  +r1(3,17)*three+r1(3,18)*three)*qz+(+r6(6,11)+r4(6,28)*six+r4(1,40) &
&                  +r2(3,19)*three+r2(1,54)*six+r0(23)*three)*xx+(+r6(3,10)+r6(3,11)*two &
&                  +r6(3,12)+r4(3,27)*six+r4(3,28)*p12+r4(3,29)*six+r2(5,18)*three &
&                  +r2(5,19)*six+r2(5,20)*three)*xz+(+r6(1,11)+r4(1,28)*six+r4(1,40) &
&                  +r2(1,19)*three+r2(1,54)*six+r0(23)*three)*zz+(+r5(3,23)+r5(3,24) &
&                  +r3(3,41)*six+r3(3,42)*six+r1(3,19)*three+r1(3,20)*three)*xxz+(+r5(1,23) &
&                  +r5(1,24)+r3(1,41)*six+r3(1,42)*six+r1(1,19)*three+r1(1,20)*three)*xzz &
&                  +rxyz(1)*xxzz
      eri(6,5,1,1)=r112+rxyz(19)*qx+(+r7(5,3)+r7(5,4)+r5(5,11)*six+r5(5,12)*six &
&                  +r3(5,9)*three+r3(5,10)*three)*qz+(+r6(5,11)+r6(5,12)+r4(5,28)*six &
&                  +r4(5,29)*six+r2(6,19)*three+r2(6,20)*three)*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,1,1)=r211+(+r7(5,3)*two+r5(5,11)*p12+r3(5,9)*six)*qx+rxyz(16)*qz+(+r6(5,10) &
&                  +r4(5,27)*six+r2(6,18)*three)*xx+(+r6(2,11)*two+r4(2,28)*p12+r2(4,19)*six) &
&                  *xz+rxyz(3)*xxz
      eri(2,6,1,1)=r031+rxyz(8)*qz
      eri(3,6,1,1)=r013+(+r7(9,3)*two+r7(9,4)+r5(9,11)*p12+r5(9,12)*six+r5(2,21)*two &
&                  +r5(2,22)+r3(9,9)*six+r3(9,10)*three+r3(2,39)*p12+r3(2,40)*six &
&                  +r1(2,17)*six+r1(2,18)*three)*qz+(+r6(5,10)+r6(5,11)*two+r4(5,27)*six &
&                  +r4(5,28)*p12+r2(6,18)*three+r2(6,19)*six)*zz+rxyz(3)*zzz
      eri(4,6,1,1)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,1,1)=r112+rxyz(18)*qx+(+r7(5,3)+r7(5,4)+r5(5,11)*six+r5(5,12)*six &
&                  +r3(5,9)*three+r3(5,10)*three)*qz+(+r6(5,10)+r6(5,11)+r4(5,27)*six &
&                  +r4(5,28)*six+r2(6,18)*three+r2(6,19)*three)*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,1,1)=r022+(+r7(8,3)+r7(8,4)+r5(8,11)*six+r5(8,12)*six+r5(3,21)+r5(3,22) &
&                  +r3(8,9)*three+r3(8,10)*three+r3(3,39)*six+r3(3,40)*six+r1(3,17)*three &
&                  +r1(3,18)*three)*qz+rxyz(6)*zz
!
      r400= r8(4)+r6(1,4)+r6(4,4)+r6(4,9)*six+r4(1,5)+r4(1,26)*six+r4(4,26)*six &
&          +r4(4,38)*three+r2(1,17)*six+r2(1,52)*three+r2(2,52)*three+r0(21)*three
      r310= r8(7)+r6(2,4)+r6(7,4)+r6(7,9)*three+r4(2,5)+r4(2,26)*three+r4(7,26)*three &
&          +r2(4,17)*three
      r301= r8(8)+r6(3,4)+r6(8,4)+r6(8,9)*three+r4(3,5)+r4(3,26)*three+r4(8,26)*three &
&          +r2(5,17)*three
      r220= r8(11)+r6(4,4)+r6(11,4)+r6(4,9)+r6(11,9)+r4(4,5)+r4(1,26)+r4(4,26)*two &
&          +r4(11,26)+r4(4,38)+r2(1,17)+r2(2,17)+r2(1,52)+r2(2,52)+r0(21)
      r211= r8(12)+r6(5,4)+r6(12,4)+r6(12,9)+r4(5,5)+r4(5,26)+r4(12,26)+r2(6,17)
      r202= r8(13)+r6(6,4)+r6(13,4)+r6(4,9)+r6(13,9)+r4(6,5)+r4(1,26)+r4(4,26)+r4(6,26) &
&          +r4(13,26)+r4(4,38)+r2(1,17)+r2(3,17)+r2(1,52)+r2(2,52)+r0(21)
      r130= r8(16)+r6(7,4)+r6(16,4)+r6(7,9)*three+r4(7,5)+r4(2,26)*three+r4(7,26)*three &
&          +r2(4,17)*three
      r121= r8(17)+r6(8,4)+r6(17,4)+r6(8,9)+r4(8,5)+r4(3,26)+r4(8,26)+r2(5,17)
      r112= r8(18)+r6(9,4)+r6(18,4)+r6(7,9)+r4(9,5)+r4(2,26)+r4(7,26)+r2(4,17)
      r103= r8(19)+r6(10,4)+r6(19,4)+r6(8,9)*three+r4(10,5)+r4(3,26)*three+r4(8,26)*three &
&          +r2(5,17)*three
      r040= r8(22)+r6(11,4)+r6(22,4)+r6(11,9)*six+r4(11,5)+r4(4,26)*six+r4(11,26)*six &
&          +r4(4,38)*three+r2(2,17)*six+r2(1,52)*three+r2(2,52)*three+r0(21)*three
      r031= r8(23)+r6(12,4)+r6(23,4)+r6(12,9)*three+r4(12,5)+r4(5,26)*three+r4(12,26)*three &
&          +r2(6,17)*three
      r022= r8(24)+r6(13,4)+r6(24,4)+r6(11,9)+r6(13,9)+r4(13,5)+r4(4,26)+r4(6,26)+r4(11,26) &
&          +r4(13,26)+r4(4,38)+r2(2,17)+r2(3,17)+r2(1,52)+r2(2,52)+r0(21)
      r013= r8(25)+r6(14,4)+r6(25,4)+r6(12,9)*three+r4(14,5)+r4(5,26)*three+r4(12,26)*three &
&          +r2(6,17)*three
      r004= r8(26)+r6(15,4)+r6(26,4)+r6(13,9)*six+r4(15,5)+r4(6,26)*six+r4(13,26)*six &
&          +r4(4,38)*three+r2(3,17)*six+r2(1,52)*three+r2(2,52)*three+r0(21)*three
      rxyz(1)=+r4(4,42)+r2(1,56)+r2(2,56)+r0(25)
      rxyz(2)=+r5(7,24)+r3(2,42)+r3(7,42)+r1(2,20)
      rxyz(3)=+r5(7,23)+r3(2,41)+r3(7,41)+r1(2,19)
      rxyz(4)=+r6(11,10)+r4(4,27)+r4(11,27)+r4(4,39)+r2(2,18)+r2(1,53)+r2(2,53)+r0(22)
      rxyz(5)=+r6(11,12)+r4(4,29)+r4(11,29)+r4(4,41)+r2(2,20)+r2(1,55)+r2(2,55)+r0(24)
      rxyz(6)=+r6(11,11)+r4(4,28)+r4(11,28)+r4(4,40)+r2(2,19)+r2(1,54)+r2(2,54)+r0(23)
      rxyz(7)=+r7(16,3)+r5(7,11)+r5(16,11)+r5(7,21)*three+r3(7,9)+r3(2,39)*three &
&             +r3(7,39)*three+r1(2,17)*three
      rxyz(8)=+r7(16,4)+r5(7,12)+r5(16,12)+r5(7,22)*three+r3(7,10)+r3(2,40)*three &
&             +r3(7,40)*three+r1(2,18)*three
      rxyz(9)=+r7(12,3)+r7(12,4)+r5(5,11)+r5(12,11)+r5(5,12)+r5(12,12)+r3(5,9)+r3(5,10)
      rxyz(10)=+r6(12,11)+r4(5,28)+r4(12,28)+r2(6,19)
      rxyz(11)=+r6(7,11)+r4(2,28)+r4(7,28)+r2(4,19)
      rxyz(12)=+r7(17,3)+r5(8,11)+r5(17,11)+r5(8,21)+r3(8,9)+r3(3,39)+r3(8,39)+r1(3,17)
      rxyz(13)=+r7(17,4)+r5(8,12)+r5(17,12)+r5(8,22)+r3(8,10)+r3(3,40)+r3(8,40)+r1(3,18)
      rxyz(14)=+r7(11,4)+r5(4,12)+r5(11,12)+r5(4,22)+r3(4,10)+r3(1,40)+r3(4,40)+r1(1,18)
      rxyz(15)=+r7(11,3)+r5(4,11)+r5(11,11)+r5(4,21)+r3(4,9)+r3(1,39)+r3(4,39)+r1(1,17)
      rxyz(16)=+r7(7,4)+r5(2,12)+r5(7,12)+r5(7,22)+r3(2,10)+r3(2,40)+r3(7,40)+r1(2,18)
      rxyz(17)=+r7(7,3)+r5(2,11)+r5(7,11)+r5(7,21)+r3(2,9)+r3(2,39)+r3(7,39)+r1(2,17)
      rxyz(18)=+r7(18,3)+r5(9,11)+r5(18,11)+r5(7,21)+r3(9,9)+r3(2,39)+r3(7,39)+r1(2,17)
      rxyz(19)=+r7(18,4)+r5(9,12)+r5(18,12)+r5(7,22)+r3(9,10)+r3(2,40)+r3(7,40)+r1(2,18)
      rxyz(20)=+r6(8,11)*four+r4(3,28)*four+r4(8,28)*four+r2(5,19)*four
      eri(1,1,2,1)=r400+(+r7(4,3)*two+r7(4,4)*two+r5(1,11)*two+r5(4,11)*two+r5(1,12)*two &
&                  +r5(4,12)*two+r5(4,21)*six+r5(4,22)*six+r3(1,9)*two+r3(1,10)*two &
&                  +r3(1,39)*six+r3(4,39)*six+r3(1,40)*six+r3(4,40)*six+r1(1,17)*six &
&                  +r1(1,18)*six)*qx+(+r6(4,10)+r6(4,11)*four+r6(4,12)+r4(1,27)+r4(4,27) &
&                  +r4(1,28)*four+r4(4,28)*four+r4(1,29)+r4(4,29)+r4(4,39)+r4(4,40)*four &
&                  +r4(4,41)+r2(1,18)+r2(1,19)*four+r2(1,20)+r2(1,53)+r2(2,53)+r2(1,54)*four &
&                  +r2(2,54)*four+r2(1,55)+r2(2,55)+r0(22)+r0(23)*four+r0(24))*xx+( &
&                  +r5(4,23)*two+r5(4,24)*two+r3(1,41)*two+r3(4,41)*two+r3(1,42)*two &
&                  +r3(4,42)*two+r1(1,19)*two+r1(1,20)*two)*xxx+rxyz(1)*xxxx
      eri(2,1,2,1)=r220+(+r7(11,4)*two+r5(4,12)*two+r5(11,12)*two+r5(4,22)*two &
&                  +r3(4,10)*two+r3(1,40)*two+r3(4,40)*two+r1(1,18)*two)*qx+rxyz(5)*xx
      eri(3,1,2,1)=r202+(+r7(13,4)*two+r5(6,12)*two+r5(13,12)*two+r5(4,22)*two &
&                  +r3(6,10)*two+r3(1,40)*two+r3(4,40)*two+r1(1,18)*two)*qx+(+r7(8,3)*two &
&                  +r5(3,11)*two+r5(8,11)*two+r5(8,21)*two+r3(3,9)*two+r3(3,39)*two &
&                  +r3(8,39)*two+r1(3,17)*two)*qz+(+r6(13,12)+r4(6,29)+r4(13,29)+r4(4,41) &
&                  +r2(3,20)+r2(1,55)+r2(2,55)+r0(24))*xx+rxyz(20)*xz+(+r6(4,10)+r4(1,27) &
&                  +r4(4,27)+r4(4,39)+r2(1,18)+r2(1,53)+r2(2,53)+r0(22))*zz+(+r5(8,24)*two &
&                  +r3(3,42)*two+r3(8,42)*two+r1(3,20)*two)*xxz+(+r5(4,23)*two+r3(1,41)*two &
&                  +r3(4,41)*two+r1(1,19)*two)*xzz+rxyz(1)*xxzz
      eri(4,1,2,1)=r310+(+r7(7,3)+r7(7,4)*two+r5(2,11)+r5(7,11)+r5(2,12)*two+r5(7,12)*two &
&                  +r5(7,21)+r5(7,22)*two+r3(2,9)+r3(2,10)*two+r3(2,39)+r3(7,39)+r3(2,40)*two &
&                  +r3(7,40)*two+r1(2,17)+r1(2,18)*two)*qx+(+r6(7,11)*two+r6(7,12) &
&                  +r4(2,28)*two+r4(7,28)*two+r4(2,29)+r4(7,29)+r2(4,19)*two+r2(4,20))*xx &
&                  +rxyz(2)*xxx
      eri(5,1,2,1)=r301+(+r7(8,3)+r7(8,4)*two+r5(3,11)+r5(8,11)+r5(3,12)*two+r5(8,12)*two &
&                  +r5(8,21)+r5(8,22)*two+r3(3,9)+r3(3,10)*two+r3(3,39)+r3(8,39)+r3(3,40)*two &
&                  +r3(8,40)*two+r1(3,17)+r1(3,18)*two)*qx+(+r7(4,3)+r5(1,11)+r5(4,11) &
&                  +r5(4,21)*three+r3(1,9)+r3(1,39)*three+r3(4,39)*three+r1(1,17)*three)*qz+( &
&                  +r6(8,11)*two+r6(8,12)+r4(3,28)*two+r4(8,28)*two+r4(3,29)+r4(8,29) &
&                  +r2(5,19)*two+r2(5,20))*xx+(+r6(4,10)+r6(4,11)*two+r4(1,27)+r4(4,27) &
&                  +r4(1,28)*two+r4(4,28)*two+r4(4,39)+r4(4,40)*two+r2(1,18)+r2(1,19)*two &
&                  +r2(1,53)+r2(2,53)+r2(1,54)*two+r2(2,54)*two+r0(22)+r0(23)*two)*xz+( &
&                  +r5(8,24)+r3(3,42)+r3(8,42)+r1(3,20))*xxx+(+r5(4,23)*two+r5(4,24) &
&                  +r3(1,41)*two+r3(4,41)*two+r3(1,42)+r3(4,42)+r1(1,19)*two+r1(1,20))*xxz &
&                  +rxyz(1)*xxxz
      eri(6,1,2,1)=r211+(+r7(12,4)*two+r5(5,12)*two+r5(12,12)*two+r3(5,10)*two)*qx &
&                  +rxyz(17)*qz+(+r6(12,12)+r4(5,29)+r4(12,29)+r2(6,20))*xx+(+r6(7,11)*two &
&                  +r4(2,28)*two+r4(7,28)*two+r2(4,19)*two)*xz+rxyz(2)*xxz
      eri(1,2,2,1)=r220+(+r7(11,3)*two+r5(4,11)*two+r5(11,11)*two+r5(4,21)*two &
&                  +r3(4,9)*two+r3(1,39)*two+r3(4,39)*two+r1(1,17)*two)*qx+rxyz(4)*xx
      eri(2,2,2,1)=r040
      eri(3,2,2,1)=r022+(+r7(17,3)*two+r5(8,11)*two+r5(17,11)*two+r5(8,21)*two &
&                  +r3(8,9)*two+r3(3,39)*two+r3(8,39)*two+r1(3,17)*two)*qz+rxyz(4)*zz
      eri(4,2,2,1)=r130+rxyz(7)*qx
      eri(5,2,2,1)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,2,1)=r031+rxyz(7)*qz
      eri(1,3,2,1)=r202+(+r7(13,3)*two+r5(6,11)*two+r5(13,11)*two+r5(4,21)*two &
&                  +r3(6,9)*two+r3(1,39)*two+r3(4,39)*two+r1(1,17)*two)*qx+(+r7(8,4)*two &
&                  +r5(3,12)*two+r5(8,12)*two+r5(8,22)*two+r3(3,10)*two+r3(3,40)*two &
&                  +r3(8,40)*two+r1(3,18)*two)*qz+(+r6(13,10)+r4(6,27)+r4(13,27)+r4(4,39) &
&                  +r2(3,18)+r2(1,53)+r2(2,53)+r0(22))*xx+rxyz(20)*xz+(+r6(4,12)+r4(1,29) &
&                  +r4(4,29)+r4(4,41)+r2(1,20)+r2(1,55)+r2(2,55)+r0(24))*zz+(+r5(8,23)*two &
&                  +r3(3,41)*two+r3(8,41)*two+r1(3,19)*two)*xxz+(+r5(4,24)*two+r3(1,42)*two &
&                  +r3(4,42)*two+r1(1,20)*two)*xzz+rxyz(1)*xxzz
      eri(2,3,2,1)=r022+(+r7(17,4)*two+r5(8,12)*two+r5(17,12)*two+r5(8,22)*two &
&                  +r3(8,10)*two+r3(3,40)*two+r3(8,40)*two+r1(3,18)*two)*qz+rxyz(5)*zz
      eri(3,3,2,1)=r004+(+r7(19,3)*two+r7(19,4)*two+r5(10,11)*two+r5(19,11)*two &
&                  +r5(10,12)*two+r5(19,12)*two+r5(8,21)*six+r5(8,22)*six+r3(10,9)*two &
&                  +r3(10,10)*two+r3(3,39)*six+r3(8,39)*six+r3(3,40)*six+r3(8,40)*six &
&                  +r1(3,17)*six+r1(3,18)*six)*qz+(+r6(13,10)+r6(13,11)*four+r6(13,12) &
&                  +r4(6,27)+r4(13,27)+r4(6,28)*four+r4(13,28)*four+r4(6,29)+r4(13,29) &
&                  +r4(4,39)+r4(4,40)*four+r4(4,41)+r2(3,18)+r2(3,19)*four+r2(3,20)+r2(1,53) &
&                  +r2(2,53)+r2(1,54)*four+r2(2,54)*four+r2(1,55)+r2(2,55)+r0(22)+r0(23)*four &
&                  +r0(24))*zz+(+r5(8,23)*two+r5(8,24)*two+r3(3,41)*two+r3(8,41)*two &
&                  +r3(3,42)*two+r3(8,42)*two+r1(3,19)*two+r1(3,20)*two)*zzz+rxyz(1)*zzzz
      eri(4,3,2,1)=r112+rxyz(18)*qx+(+r7(12,4)*two+r5(5,12)*two+r5(12,12)*two &
&                  +r3(5,10)*two)*qz+(+r6(12,11)*two+r4(5,28)*two+r4(12,28)*two+r2(6,19)*two) &
&                  *xz+(+r6(7,12)+r4(2,29)+r4(7,29)+r2(4,20))*zz+rxyz(2)*xzz
      eri(5,3,2,1)=r103+(+r7(19,3)+r5(10,11)+r5(19,11)+r5(8,21)*three+r3(10,9) &
&                  +r3(3,39)*three+r3(8,39)*three+r1(3,17)*three)*qx+(+r7(13,3)+r7(13,4)*two &
&                  +r5(6,11)+r5(13,11)+r5(6,12)*two+r5(13,12)*two+r5(4,21)+r5(4,22)*two &
&                  +r3(6,9)+r3(6,10)*two+r3(1,39)+r3(4,39)+r3(1,40)*two+r3(4,40)*two+r1(1,17) &
&                  +r1(1,18)*two)*qz+(+r6(13,10)+r6(13,11)*two+r4(6,27)+r4(13,27) &
&                  +r4(6,28)*two+r4(13,28)*two+r4(4,39)+r4(4,40)*two+r2(3,18)+r2(3,19)*two &
&                  +r2(1,53)+r2(2,53)+r2(1,54)*two+r2(2,54)*two+r0(22)+r0(23)*two)*xz+( &
&                  +r6(8,11)*two+r6(8,12)+r4(3,28)*two+r4(8,28)*two+r4(3,29)+r4(8,29) &
&                  +r2(5,19)*two+r2(5,20))*zz+(+r5(8,23)*two+r5(8,24)+r3(3,41)*two &
&                  +r3(8,41)*two+r3(3,42)+r3(8,42)+r1(3,19)*two+r1(3,20))*xzz+(+r5(4,24) &
&                  +r3(1,42)+r3(4,42)+r1(1,20))*zzz+rxyz(1)*xzzz
      eri(6,3,2,1)=r013+(+r7(18,3)+r7(18,4)*two+r5(9,11)+r5(18,11)+r5(9,12)*two &
&                  +r5(18,12)*two+r5(7,21)+r5(7,22)*two+r3(9,9)+r3(9,10)*two+r3(2,39) &
&                  +r3(7,39)+r3(2,40)*two+r3(7,40)*two+r1(2,17)+r1(2,18)*two)*qz+( &
&                  +r6(12,11)*two+r6(12,12)+r4(5,28)*two+r4(12,28)*two+r4(5,29)+r4(12,29) &
&                  +r2(6,19)*two+r2(6,20))*zz+rxyz(2)*zzz
      eri(1,4,2,1)=r310+(+r7(7,3)*two+r7(7,4)+r5(2,11)*two+r5(7,11)*two+r5(2,12)+r5(7,12) &
&                  +r5(7,21)*two+r5(7,22)+r3(2,9)*two+r3(2,10)+r3(2,39)*two+r3(7,39)*two &
&                  +r3(2,40)+r3(7,40)+r1(2,17)*two+r1(2,18))*qx+(+r6(7,10)+r6(7,11)*two &
&                  +r4(2,27)+r4(7,27)+r4(2,28)*two+r4(7,28)*two+r2(4,18)+r2(4,19)*two)*xx &
&                  +rxyz(3)*xxx
      eri(2,4,2,1)=r130+rxyz(8)*qx
      eri(3,4,2,1)=r112+rxyz(19)*qx+(+r7(12,3)*two+r5(5,11)*two+r5(12,11)*two+r3(5,9)*two &
&                  )*qz+(+r6(12,11)*two+r4(5,28)*two+r4(12,28)*two+r2(6,19)*two)*xz+( &
&                  +r6(7,10)+r4(2,27)+r4(7,27)+r2(4,18))*zz+rxyz(3)*xzz
      eri(4,4,2,1)=r220+(+r7(11,3)+r7(11,4)+r5(4,11)+r5(11,11)+r5(4,12)+r5(11,12) &
&                  +r5(4,21)+r5(4,22)+r3(4,9)+r3(4,10)+r3(1,39)+r3(4,39)+r3(1,40)+r3(4,40) &
&                  +r1(1,17)+r1(1,18))*qx+rxyz(6)*xx
      eri(5,4,2,1)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(7,10)+r6(7,11)+r4(2,27) &
&                  +r4(7,27)+r4(2,28)+r4(7,28)+r2(4,18)+r2(4,19))*xz+rxyz(3)*xxz
      eri(6,4,2,1)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,2,1)=r301+(+r7(8,3)*two+r7(8,4)+r5(3,11)*two+r5(8,11)*two+r5(3,12)+r5(8,12) &
&                  +r5(8,21)*two+r5(8,22)+r3(3,9)*two+r3(3,10)+r3(3,39)*two+r3(8,39)*two &
&                  +r3(3,40)+r3(8,40)+r1(3,17)*two+r1(3,18))*qx+(+r7(4,4)+r5(1,12)+r5(4,12) &
&                  +r5(4,22)*three+r3(1,10)+r3(1,40)*three+r3(4,40)*three+r1(1,18)*three)*qz &
&                  +(+r6(8,10)+r6(8,11)*two+r4(3,27)+r4(8,27)+r4(3,28)*two+r4(8,28)*two &
&                  +r2(5,18)+r2(5,19)*two)*xx+(+r6(4,11)*two+r6(4,12)+r4(1,28)*two &
&                  +r4(4,28)*two+r4(1,29)+r4(4,29)+r4(4,40)*two+r4(4,41)+r2(1,19)*two &
&                  +r2(1,20)+r2(1,54)*two+r2(2,54)*two+r2(1,55)+r2(2,55)+r0(23)*two+r0(24)) &
&                  *xz+(+r5(8,23)+r3(3,41)+r3(8,41)+r1(3,19))*xxx+(+r5(4,23)+r5(4,24)*two &
&                  +r3(1,41)+r3(4,41)+r3(1,42)*two+r3(4,42)*two+r1(1,19)+r1(1,20)*two)*xxz &
&                  +rxyz(1)*xxxz
      eri(2,5,2,1)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,2,1)=r103+(+r7(19,4)+r5(10,12)+r5(19,12)+r5(8,22)*three+r3(10,10) &
&                  +r3(3,40)*three+r3(8,40)*three+r1(3,18)*three)*qx+(+r7(13,3)*two+r7(13,4) &
&                  +r5(6,11)*two+r5(13,11)*two+r5(6,12)+r5(13,12)+r5(4,21)*two+r5(4,22) &
&                  +r3(6,9)*two+r3(6,10)+r3(1,39)*two+r3(4,39)*two+r3(1,40)+r3(4,40) &
&                  +r1(1,17)*two+r1(1,18))*qz+(+r6(13,11)*two+r6(13,12)+r4(6,28)*two &
&                  +r4(13,28)*two+r4(6,29)+r4(13,29)+r4(4,40)*two+r4(4,41)+r2(3,19)*two &
&                  +r2(3,20)+r2(1,54)*two+r2(2,54)*two+r2(1,55)+r2(2,55)+r0(23)*two+r0(24)) &
&                  *xz+(+r6(8,10)+r6(8,11)*two+r4(3,27)+r4(8,27)+r4(3,28)*two+r4(8,28)*two &
&                  +r2(5,18)+r2(5,19)*two)*zz+(+r5(8,23)+r5(8,24)*two+r3(3,41)+r3(8,41) &
&                  +r3(3,42)*two+r3(8,42)*two+r1(3,19)+r1(3,20)*two)*xzz+(+r5(4,23)+r3(1,41) &
&                  +r3(4,41)+r1(1,19))*zzz+rxyz(1)*xzzz
      eri(4,5,2,1)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(7,11)+r6(7,12)+r4(2,28) &
&                  +r4(7,28)+r4(2,29)+r4(7,29)+r2(4,19)+r2(4,20))*xz+rxyz(2)*xxz
      eri(5,5,2,1)=r202+(+r7(13,3)+r7(13,4)+r5(6,11)+r5(13,11)+r5(6,12)+r5(13,12) &
&                  +r5(4,21)+r5(4,22)+r3(6,9)+r3(6,10)+r3(1,39)+r3(4,39)+r3(1,40)+r3(4,40) &
&                  +r1(1,17)+r1(1,18))*qx+(+r7(8,3)+r7(8,4)+r5(3,11)+r5(8,11)+r5(3,12) &
&                  +r5(8,12)+r5(8,21)+r5(8,22)+r3(3,9)+r3(3,10)+r3(3,39)+r3(8,39)+r3(3,40) &
&                  +r3(8,40)+r1(3,17)+r1(3,18))*qz+(+r6(13,11)+r4(6,28)+r4(13,28)+r4(4,40) &
&                  +r2(3,19)+r2(1,54)+r2(2,54)+r0(23))*xx+(+r6(8,10)+r6(8,11)*two+r6(8,12) &
&                  +r4(3,27)+r4(8,27)+r4(3,28)*two+r4(8,28)*two+r4(3,29)+r4(8,29)+r2(5,18) &
&                  +r2(5,19)*two+r2(5,20))*xz+(+r6(4,11)+r4(1,28)+r4(4,28)+r4(4,40)+r2(1,19) &
&                  +r2(1,54)+r2(2,54)+r0(23))*zz+(+r5(8,23)+r5(8,24)+r3(3,41)+r3(8,41) &
&                  +r3(3,42)+r3(8,42)+r1(3,19)+r1(3,20))*xxz+(+r5(4,23)+r5(4,24)+r3(1,41) &
&                  +r3(4,41)+r3(1,42)+r3(4,42)+r1(1,19)+r1(1,20))*xzz+rxyz(1)*xxzz
      eri(6,5,2,1)=r112+rxyz(19)*qx+(+r7(12,3)+r7(12,4)+r5(5,11)+r5(12,11)+r5(5,12) &
&                  +r5(12,12)+r3(5,9)+r3(5,10))*qz+(+r6(12,11)+r6(12,12)+r4(5,28)+r4(12,28) &
&                  +r4(5,29)+r4(12,29)+r2(6,19)+r2(6,20))*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,2,1)=r211+(+r7(12,3)*two+r5(5,11)*two+r5(12,11)*two+r3(5,9)*two)*qx &
&                  +rxyz(16)*qz+(+r6(12,10)+r4(5,27)+r4(12,27)+r2(6,18))*xx+(+r6(7,11)*two &
&                  +r4(2,28)*two+r4(7,28)*two+r2(4,19)*two)*xz+rxyz(3)*xxz
      eri(2,6,2,1)=r031+rxyz(8)*qz
      eri(3,6,2,1)=r013+(+r7(18,3)*two+r7(18,4)+r5(9,11)*two+r5(18,11)*two+r5(9,12) &
&                  +r5(18,12)+r5(7,21)*two+r5(7,22)+r3(9,9)*two+r3(9,10)+r3(2,39)*two &
&                  +r3(7,39)*two+r3(2,40)+r3(7,40)+r1(2,17)*two+r1(2,18))*qz+(+r6(12,10) &
&                  +r6(12,11)*two+r4(5,27)+r4(12,27)+r4(5,28)*two+r4(12,28)*two+r2(6,18) &
&                  +r2(6,19)*two)*zz+rxyz(3)*zzz
      eri(4,6,2,1)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,2,1)=r112+rxyz(18)*qx+(+r7(12,3)+r7(12,4)+r5(5,11)+r5(12,11)+r5(5,12) &
&                  +r5(12,12)+r3(5,9)+r3(5,10))*qz+(+r6(12,10)+r6(12,11)+r4(5,27)+r4(12,27) &
&                  +r4(5,28)+r4(12,28)+r2(6,18)+r2(6,19))*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,2,1)=r022+(+r7(17,3)+r7(17,4)+r5(8,11)+r5(17,11)+r5(8,12)+r5(17,12) &
&                  +r5(8,21)+r5(8,22)+r3(8,9)+r3(8,10)+r3(3,39)+r3(8,39)+r3(3,40)+r3(8,40) &
&                  +r1(3,17)+r1(3,18))*qz+rxyz(6)*zz
!
      r400= r8(6)-r7(3,1)*two+r6(1,1)+r6(1,4)+r6(6,4)+r6(6,9)*six-r5(3,3)*two-r5(3,13)*p12 &
&          +r4(1,2)+r4(1,5)+r4(1,14)*six+r4(1,26)*six+r4(6,26)*six+r4(6,38)*three &
&          -r3(3,19)*p12-r3(3,43)*six+r2(1,5)*six+r2(1,17)*six+r2(1,37)*three &
&          +r2(1,52)*three+r2(3,52)*three-r1(3,31)*six+r0(6)*three+r0(21)*three
      r310= r8(9)-r7(5,1)*two+r6(2,1)+r6(2,4)+r6(9,4)+r6(9,9)*three-r5(5,3)*two &
&          -r5(5,13)*six+r4(2,2)+r4(2,5)+r4(2,14)*three+r4(2,26)*three+r4(9,26)*three &
&          -r3(5,19)*six+r2(4,5)*three+r2(4,17)*three
      r301= r8(10)-r7(6,1)*two+r6(3,1)+r6(3,4)+r6(10,4)+r6(10,9)*three-r5(6,3)*two &
&          -r5(6,13)*six+r4(3,2)+r4(3,5)+r4(3,14)*three+r4(3,26)*three+r4(10,26)*three &
&          -r3(6,19)*six+r2(5,5)*three+r2(5,17)*three
      r220= r8(13)-r7(8,1)*two+r6(4,1)+r6(4,4)+r6(13,4)+r6(6,9)+r6(13,9)-r5(8,3)*two &
&          -r5(3,13)*two-r5(8,13)*two+r4(4,2)+r4(4,5)+r4(1,14)+r4(4,14)+r4(1,26)+r4(4,26) &
&          +r4(6,26)+r4(13,26)+r4(6,38)-r3(3,19)*two-r3(8,19)*two-r3(3,43)*two+r2(1,5) &
&          +r2(2,5)+r2(1,17)+r2(2,17)+r2(1,37)+r2(1,52)+r2(3,52)-r1(3,31)*two+r0(6)+r0(21)
      r211= r8(14)-r7(9,1)*two+r6(5,1)+r6(5,4)+r6(14,4)+r6(14,9)-r5(9,3)*two-r5(9,13)*two &
&          +r4(5,2)+r4(5,5)+r4(5,14)+r4(5,26)+r4(14,26)-r3(9,19)*two+r2(6,5)+r2(6,17)
      r202= r8(15)-r7(10,1)*two+r6(6,1)+r6(6,4)+r6(15,4)+r6(6,9)+r6(15,9)-r5(10,3)*two &
&          -r5(3,13)*two-r5(10,13)*two+r4(6,2)+r4(6,5)+r4(1,14)+r4(6,14)+r4(1,26) &
&          +r4(6,26)*two+r4(15,26)+r4(6,38)-r3(3,19)*two-r3(10,19)*two-r3(3,43)*two &
&          +r2(1,5)+r2(3,5)+r2(1,17)+r2(3,17)+r2(1,37)+r2(1,52)+r2(3,52)-r1(3,31)*two &
&          +r0(6)+r0(21)
      r130= r8(18)-r7(12,1)*two+r6(7,1)+r6(7,4)+r6(18,4)+r6(9,9)*three-r5(12,3)*two &
&          -r5(5,13)*six+r4(7,2)+r4(7,5)+r4(2,14)*three+r4(2,26)*three+r4(9,26)*three &
&          -r3(5,19)*six+r2(4,5)*three+r2(4,17)*three
      r121= r8(19)-r7(13,1)*two+r6(8,1)+r6(8,4)+r6(19,4)+r6(10,9)-r5(13,3)*two-r5(6,13)*two &
&          +r4(8,2)+r4(8,5)+r4(3,14)+r4(3,26)+r4(10,26)-r3(6,19)*two+r2(5,5)+r2(5,17)
      r112= r8(20)-r7(14,1)*two+r6(9,1)+r6(9,4)+r6(20,4)+r6(9,9)-r5(14,3)*two-r5(5,13)*two &
&          +r4(9,2)+r4(9,5)+r4(2,14)+r4(2,26)+r4(9,26)-r3(5,19)*two+r2(4,5)+r2(4,17)
      r103= r8(21)-r7(15,1)*two+r6(10,1)+r6(10,4)+r6(21,4)+r6(10,9)*three-r5(15,3)*two &
&          -r5(6,13)*six+r4(10,2)+r4(10,5)+r4(3,14)*three+r4(3,26)*three+r4(10,26)*three &
&          -r3(6,19)*six+r2(5,5)*three+r2(5,17)*three
      r040= r8(24)-r7(17,1)*two+r6(11,1)+r6(11,4)+r6(24,4)+r6(13,9)*six-r5(17,3)*two &
&          -r5(8,13)*p12+r4(11,2)+r4(11,5)+r4(4,14)*six+r4(4,26)*six+r4(13,26)*six &
&          +r4(6,38)*three-r3(8,19)*p12-r3(3,43)*six+r2(2,5)*six+r2(2,17)*six &
&          +r2(1,37)*three+r2(1,52)*three+r2(3,52)*three-r1(3,31)*six+r0(6)*three &
&          +r0(21)*three
      r031= r8(25)-r7(18,1)*two+r6(12,1)+r6(12,4)+r6(25,4)+r6(14,9)*three-r5(18,3)*two &
&          -r5(9,13)*six+r4(12,2)+r4(12,5)+r4(5,14)*three+r4(5,26)*three+r4(14,26)*three &
&          -r3(9,19)*six+r2(6,5)*three+r2(6,17)*three
      r022= r8(26)-r7(19,1)*two+r6(13,1)+r6(13,4)+r6(26,4)+r6(13,9)+r6(15,9)-r5(19,3)*two &
&          -r5(8,13)*two-r5(10,13)*two+r4(13,2)+r4(13,5)+r4(4,14)+r4(6,14)+r4(4,26) &
&          +r4(6,26)+r4(13,26)+r4(15,26)+r4(6,38)-r3(8,19)*two-r3(10,19)*two-r3(3,43)*two &
&          +r2(2,5)+r2(3,5)+r2(2,17)+r2(3,17)+r2(1,37)+r2(1,52)+r2(3,52)-r1(3,31)*two &
&          +r0(6)+r0(21)
      r013= r8(27)-r7(20,1)*two+r6(14,1)+r6(14,4)+r6(27,4)+r6(14,9)*three-r5(20,3)*two &
&          -r5(9,13)*six+r4(14,2)+r4(14,5)+r4(5,14)*three+r4(5,26)*three+r4(14,26)*three &
&          -r3(9,19)*six+r2(6,5)*three+r2(6,17)*three
      r004= r8(28)-r7(21,1)*two+r6(15,1)+r6(15,4)+r6(28,4)+r6(15,9)*six-r5(21,3)*two &
&          -r5(10,13)*p12+r4(15,2)+r4(15,5)+r4(6,14)*six+r4(6,26)*six+r4(15,26)*six &
&          +r4(6,38)*three-r3(10,19)*p12-r3(3,43)*six+r2(3,5)*six+r2(3,17)*six &
&          +r2(1,37)*three+r2(1,52)*three+r2(3,52)*three-r1(3,31)*six+r0(6)*three &
&          +r0(21)*three
      rxyz(1)=+r4(6,42)-r3(3,47)*two+r2(1,41)+r2(1,56)+r2(3,56)-r1(3,35)*two+r0(10)+r0(25)
      rxyz(2)=+r5(9,24)-r4(5,33)*two+r3(2,30)+r3(2,42)+r3(9,42)-r2(6,32)*two+r1(2,8) &
&             +r1(2,20)
      rxyz(3)=+r5(9,23)-r4(5,32)*two+r3(2,29)+r3(2,41)+r3(9,41)-r2(6,31)*two+r1(2,7) &
&             +r1(2,19)
      rxyz(4)=+r6(13,10)-r5(8,14)*two+r4(4,15)+r4(4,27)+r4(13,27)+r4(6,39)-r3(8,20)*two &
&             -r3(3,44)*two+r2(2,6)+r2(2,18)+r2(1,38)+r2(1,53)+r2(3,53)-r1(3,32)*two+r0(7) &
&             +r0(22)
      rxyz(5)=+r6(13,12)-r5(8,16)*two+r4(4,17)+r4(4,29)+r4(13,29)+r4(6,41)-r3(8,22)*two &
&             -r3(3,46)*two+r2(2,8)+r2(2,20)+r2(1,40)+r2(1,55)+r2(3,55)-r1(3,34)*two+r0(9) &
&             +r0(24)
      rxyz(6)=+r6(13,11)-r5(8,15)*two+r4(4,16)+r4(4,28)+r4(13,28)+r4(6,40)-r3(8,21)*two &
&             -r3(3,45)*two+r2(2,7)+r2(2,19)+r2(1,39)+r2(1,54)+r2(3,54)-r1(3,33)*two+r0(8) &
&             +r0(23)
      rxyz(7)=+r7(18,3)-r6(12,5)*two+r5(7,5)+r5(7,11)+r5(18,11)+r5(9,21)*three &
&             -r4(12,10)*two-r4(5,30)*six+r3(7,3)+r3(7,9)+r3(2,27)*three+r3(2,39)*three &
&             +r3(9,39)*three-r2(6,29)*six+r1(2,5)*three+r1(2,17)*three
      rxyz(8)=+r7(18,4)-r6(12,6)*two+r5(7,6)+r5(7,12)+r5(18,12)+r5(9,22)*three &
&             -r4(12,11)*two-r4(5,31)*six+r3(7,4)+r3(7,10)+r3(2,28)*three+r3(2,40)*three &
&             +r3(9,40)*three-r2(6,30)*six+r1(2,6)*three+r1(2,18)*three
      rxyz(9)=+r7(14,3)+r7(14,4)-r6(9,5)*two-r6(9,6)*two+r5(5,5)+r5(5,6)+r5(5,11) &
&             +r5(14,11)+r5(5,12)+r5(14,12)-r4(9,10)*two-r4(9,11)*two+r3(5,3)+r3(5,4) &
&             +r3(5,9)+r3(5,10)
      rxyz(10)=+r6(14,11)-r5(9,15)*two+r4(5,16)+r4(5,28)+r4(14,28)-r3(9,21)*two+r2(6,7) &
&             +r2(6,19)
      rxyz(11)=+r6(9,11)-r5(5,15)*two+r4(2,16)+r4(2,28)+r4(9,28)-r3(5,21)*two+r2(4,7) &
&             +r2(4,19)
      rxyz(12)=+r7(19,3)-r6(13,5)*two+r5(8,5)+r5(8,11)+r5(19,11)+r5(10,21)-r4(13,10)*two &
&             -r4(6,30)*two+r3(8,3)+r3(8,9)+r3(3,27)+r3(3,39)+r3(10,39)-r2(3,29)*two &
&             +r1(3,5)+r1(3,17)
      rxyz(13)=+r7(19,4)-r6(13,6)*two+r5(8,6)+r5(8,12)+r5(19,12)+r5(10,22)-r4(13,11)*two &
&             -r4(6,31)*two+r3(8,4)+r3(8,10)+r3(3,28)+r3(3,40)+r3(10,40)-r2(3,30)*two &
&             +r1(3,6)+r1(3,18)
      rxyz(14)=+r7(13,4)-r6(8,6)*two+r5(4,6)+r5(4,12)+r5(13,12)+r5(6,22)-r4(8,11)*two &
&             -r4(3,31)*two+r3(4,4)+r3(4,10)+r3(1,28)+r3(1,40)+r3(6,40)-r2(5,30)*two &
&             +r1(1,6)+r1(1,18)
      rxyz(15)=+r7(13,3)-r6(8,5)*two+r5(4,5)+r5(4,11)+r5(13,11)+r5(6,21)-r4(8,10)*two &
&             -r4(3,30)*two+r3(4,3)+r3(4,9)+r3(1,27)+r3(1,39)+r3(6,39)-r2(5,29)*two &
&             +r1(1,5)+r1(1,17)
      rxyz(16)=+r7(9,4)-r6(5,6)*two+r5(2,6)+r5(2,12)+r5(9,12)+r5(9,22)-r4(5,11)*two &
&             -r4(5,31)*two+r3(2,4)+r3(2,10)+r3(2,28)+r3(2,40)+r3(9,40)-r2(6,30)*two &
&             +r1(2,6)+r1(2,18)
      rxyz(17)=+r7(9,3)-r6(5,5)*two+r5(2,5)+r5(2,11)+r5(9,11)+r5(9,21)-r4(5,10)*two &
&             -r4(5,30)*two+r3(2,3)+r3(2,9)+r3(2,27)+r3(2,39)+r3(9,39)-r2(6,29)*two &
&             +r1(2,5)+r1(2,17)
      rxyz(18)=+r7(20,3)-r6(14,5)*two+r5(9,5)+r5(9,11)+r5(20,11)+r5(9,21)-r4(14,10)*two &
&             -r4(5,30)*two+r3(9,3)+r3(9,9)+r3(2,27)+r3(2,39)+r3(9,39)-r2(6,29)*two &
&             +r1(2,5)+r1(2,17)
      rxyz(19)=+r7(20,4)-r6(14,6)*two+r5(9,6)+r5(9,12)+r5(20,12)+r5(9,22)-r4(14,11)*two &
&             -r4(5,31)*two+r3(9,4)+r3(9,10)+r3(2,28)+r3(2,40)+r3(9,40)-r2(6,30)*two &
&             +r1(2,6)+r1(2,18)
      rxyz(20)=+r6(10,11)*four-r5(6,15)*eight+r4(3,16)*four+r4(3,28)*four+r4(10,28)*four &
&             -r3(6,21)*eight+r2(5,7)*four+r2(5,19)*four
      eri(1,1,3,1)=r400+(+r7(6,3)*two+r7(6,4)*two-r6(3,5)*four-r6(3,6)*four+r5(1,5)*two &
&                  +r5(1,6)*two+r5(1,11)*two+r5(6,11)*two+r5(1,12)*two+r5(6,12)*two &
&                  +r5(6,21)*six+r5(6,22)*six-r4(3,10)*four-r4(3,11)*four-r4(3,30)*p12 &
&                  -r4(3,31)*p12+r3(1,3)*two+r3(1,4)*two+r3(1,9)*two+r3(1,10)*two &
&                  +r3(1,27)*six+r3(1,28)*six+r3(1,39)*six+r3(6,39)*six+r3(1,40)*six &
&                  +r3(6,40)*six-r2(5,29)*p12-r2(5,30)*p12+r1(1,5)*six+r1(1,6)*six &
&                  +r1(1,17)*six+r1(1,18)*six)*qx+(+r6(6,10)+r6(6,11)*four+r6(6,12) &
&                  -r5(3,14)*two-r5(3,15)*eight-r5(3,16)*two+r4(1,15)+r4(1,16)*four+r4(1,17) &
&                  +r4(1,27)+r4(6,27)+r4(1,28)*four+r4(6,28)*four+r4(1,29)+r4(6,29)+r4(6,39) &
&                  +r4(6,40)*four+r4(6,41)-r3(3,20)*two-r3(3,21)*eight-r3(3,22)*two &
&                  -r3(3,44)*two-r3(3,45)*eight-r3(3,46)*two+r2(1,6)+r2(1,7)*four+r2(1,8) &
&                  +r2(1,18)+r2(1,19)*four+r2(1,20)+r2(1,38)+r2(1,39)*four+r2(1,40)+r2(1,53) &
&                  +r2(3,53)+r2(1,54)*four+r2(3,54)*four+r2(1,55)+r2(3,55)-r1(3,32)*two &
&                  -r1(3,33)*eight-r1(3,34)*two+r0(7)+r0(8)*four+r0(9)+r0(22)+r0(23)*four &
&                  +r0(24))*xx+(+r5(6,23)*two+r5(6,24)*two-r4(3,32)*four-r4(3,33)*four &
&                  +r3(1,29)*two+r3(1,30)*two+r3(1,41)*two+r3(6,41)*two+r3(1,42)*two &
&                  +r3(6,42)*two-r2(5,31)*four-r2(5,32)*four+r1(1,7)*two+r1(1,8)*two &
&                  +r1(1,19)*two+r1(1,20)*two)*xxx+rxyz(1)*xxxx
      eri(2,1,3,1)=r220+(+r7(13,4)*two-r6(8,6)*four+r5(4,6)*two+r5(4,12)*two &
&                  +r5(13,12)*two+r5(6,22)*two-r4(8,11)*four-r4(3,31)*four+r3(4,4)*two &
&                  +r3(4,10)*two+r3(1,28)*two+r3(1,40)*two+r3(6,40)*two-r2(5,30)*four &
&                  +r1(1,6)*two+r1(1,18)*two)*qx+rxyz(5)*xx
      eri(3,1,3,1)=r202+(+r7(15,4)*two-r6(10,6)*four+r5(6,6)*two+r5(6,12)*two &
&                  +r5(15,12)*two+r5(6,22)*two-r4(10,11)*four-r4(3,31)*four+r3(6,4)*two &
&                  +r3(6,10)*two+r3(1,28)*two+r3(1,40)*two+r3(6,40)*two-r2(5,30)*four &
&                  +r1(1,6)*two+r1(1,18)*two)*qx+(+r7(10,3)*two-r6(6,5)*four+r5(3,5)*two &
&                  +r5(3,11)*two+r5(10,11)*two+r5(10,21)*two-r4(6,10)*four-r4(6,30)*four &
&                  +r3(3,3)*two+r3(3,9)*two+r3(3,27)*two+r3(3,39)*two+r3(10,39)*two &
&                  -r2(3,29)*four+r1(3,5)*two+r1(3,17)*two)*qz+(+r6(15,12)-r5(10,16)*two &
&                  +r4(6,17)+r4(6,29)+r4(15,29)+r4(6,41)-r3(10,22)*two-r3(3,46)*two+r2(3,8) &
&                  +r2(3,20)+r2(1,40)+r2(1,55)+r2(3,55)-r1(3,34)*two+r0(9)+r0(24))*xx &
&                  +rxyz(20)*xz+(+r6(6,10)-r5(3,14)*two+r4(1,15)+r4(1,27)+r4(6,27)+r4(6,39) &
&                  -r3(3,20)*two-r3(3,44)*two+r2(1,6)+r2(1,18)+r2(1,38)+r2(1,53)+r2(3,53) &
&                  -r1(3,32)*two+r0(7)+r0(22))*zz+(+r5(10,24)*two-r4(6,33)*four+r3(3,30)*two &
&                  +r3(3,42)*two+r3(10,42)*two-r2(3,32)*four+r1(3,8)*two+r1(3,20)*two)*xxz+( &
&                  +r5(6,23)*two-r4(3,32)*four+r3(1,29)*two+r3(1,41)*two+r3(6,41)*two &
&                  -r2(5,31)*four+r1(1,7)*two+r1(1,19)*two)*xzz+rxyz(1)*xxzz
      eri(4,1,3,1)=r310+(+r7(9,3)+r7(9,4)*two-r6(5,5)*two-r6(5,6)*four+r5(2,5) &
&                  +r5(2,6)*two+r5(2,11)+r5(9,11)+r5(2,12)*two+r5(9,12)*two+r5(9,21) &
&                  +r5(9,22)*two-r4(5,10)*two-r4(5,11)*four-r4(5,30)*two-r4(5,31)*four &
&                  +r3(2,3)+r3(2,4)*two+r3(2,9)+r3(2,10)*two+r3(2,27)+r3(2,28)*two+r3(2,39) &
&                  +r3(9,39)+r3(2,40)*two+r3(9,40)*two-r2(6,29)*two-r2(6,30)*four+r1(2,5) &
&                  +r1(2,6)*two+r1(2,17)+r1(2,18)*two)*qx+(+r6(9,11)*two+r6(9,12) &
&                  -r5(5,15)*four-r5(5,16)*two+r4(2,16)*two+r4(2,17)+r4(2,28)*two &
&                  +r4(9,28)*two+r4(2,29)+r4(9,29)-r3(5,21)*four-r3(5,22)*two+r2(4,7)*two &
&                  +r2(4,8)+r2(4,19)*two+r2(4,20))*xx+rxyz(2)*xxx
      eri(5,1,3,1)=r301+(+r7(10,3)+r7(10,4)*two-r6(6,5)*two-r6(6,6)*four+r5(3,5) &
&                  +r5(3,6)*two+r5(3,11)+r5(10,11)+r5(3,12)*two+r5(10,12)*two+r5(10,21) &
&                  +r5(10,22)*two-r4(6,10)*two-r4(6,11)*four-r4(6,30)*two-r4(6,31)*four &
&                  +r3(3,3)+r3(3,4)*two+r3(3,9)+r3(3,10)*two+r3(3,27)+r3(3,28)*two+r3(3,39) &
&                  +r3(10,39)+r3(3,40)*two+r3(10,40)*two-r2(3,29)*two-r2(3,30)*four+r1(3,5) &
&                  +r1(3,6)*two+r1(3,17)+r1(3,18)*two)*qx+(+r7(6,3)-r6(3,5)*two+r5(1,5) &
&                  +r5(1,11)+r5(6,11)+r5(6,21)*three-r4(3,10)*two-r4(3,30)*six+r3(1,3) &
&                  +r3(1,9)+r3(1,27)*three+r3(1,39)*three+r3(6,39)*three-r2(5,29)*six &
&                  +r1(1,5)*three+r1(1,17)*three)*qz+(+r6(10,11)*two+r6(10,12)-r5(6,15)*four &
&                  -r5(6,16)*two+r4(3,16)*two+r4(3,17)+r4(3,28)*two+r4(10,28)*two+r4(3,29) &
&                  +r4(10,29)-r3(6,21)*four-r3(6,22)*two+r2(5,7)*two+r2(5,8)+r2(5,19)*two &
&                  +r2(5,20))*xx+(+r6(6,10)+r6(6,11)*two-r5(3,14)*two-r5(3,15)*four+r4(1,15) &
&                  +r4(1,16)*two+r4(1,27)+r4(6,27)+r4(1,28)*two+r4(6,28)*two+r4(6,39) &
&                  +r4(6,40)*two-r3(3,20)*two-r3(3,21)*four-r3(3,44)*two-r3(3,45)*four &
&                  +r2(1,6)+r2(1,7)*two+r2(1,18)+r2(1,19)*two+r2(1,38)+r2(1,39)*two+r2(1,53) &
&                  +r2(3,53)+r2(1,54)*two+r2(3,54)*two-r1(3,32)*two-r1(3,33)*four+r0(7) &
&                  +r0(8)*two+r0(22)+r0(23)*two)*xz+(+r5(10,24)-r4(6,33)*two+r3(3,30) &
&                  +r3(3,42)+r3(10,42)-r2(3,32)*two+r1(3,8)+r1(3,20))*xxx+(+r5(6,23)*two &
&                  +r5(6,24)-r4(3,32)*four-r4(3,33)*two+r3(1,29)*two+r3(1,30)+r3(1,41)*two &
&                  +r3(6,41)*two+r3(1,42)+r3(6,42)-r2(5,31)*four-r2(5,32)*two+r1(1,7)*two &
&                  +r1(1,8)+r1(1,19)*two+r1(1,20))*xxz+rxyz(1)*xxxz
      eri(6,1,3,1)=r211+(+r7(14,4)*two-r6(9,6)*four+r5(5,6)*two+r5(5,12)*two &
&                  +r5(14,12)*two-r4(9,11)*four+r3(5,4)*two+r3(5,10)*two)*qx+rxyz(17)*qz+( &
&                  +r6(14,12)-r5(9,16)*two+r4(5,17)+r4(5,29)+r4(14,29)-r3(9,22)*two+r2(6,8) &
&                  +r2(6,20))*xx+(+r6(9,11)*two-r5(5,15)*four+r4(2,16)*two+r4(2,28)*two &
&                  +r4(9,28)*two-r3(5,21)*four+r2(4,7)*two+r2(4,19)*two)*xz+rxyz(2)*xxz
      eri(1,2,3,1)=r220+(+r7(13,3)*two-r6(8,5)*four+r5(4,5)*two+r5(4,11)*two &
&                  +r5(13,11)*two+r5(6,21)*two-r4(8,10)*four-r4(3,30)*four+r3(4,3)*two &
&                  +r3(4,9)*two+r3(1,27)*two+r3(1,39)*two+r3(6,39)*two-r2(5,29)*four &
&                  +r1(1,5)*two+r1(1,17)*two)*qx+rxyz(4)*xx
      eri(2,2,3,1)=r040
      eri(3,2,3,1)=r022+(+r7(19,3)*two-r6(13,5)*four+r5(8,5)*two+r5(8,11)*two &
&                  +r5(19,11)*two+r5(10,21)*two-r4(13,10)*four-r4(6,30)*four+r3(8,3)*two &
&                  +r3(8,9)*two+r3(3,27)*two+r3(3,39)*two+r3(10,39)*two-r2(3,29)*four &
&                  +r1(3,5)*two+r1(3,17)*two)*qz+rxyz(4)*zz
      eri(4,2,3,1)=r130+rxyz(7)*qx
      eri(5,2,3,1)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,3,1)=r031+rxyz(7)*qz
      eri(1,3,3,1)=r202+(+r7(15,3)*two-r6(10,5)*four+r5(6,5)*two+r5(6,11)*two &
&                  +r5(15,11)*two+r5(6,21)*two-r4(10,10)*four-r4(3,30)*four+r3(6,3)*two &
&                  +r3(6,9)*two+r3(1,27)*two+r3(1,39)*two+r3(6,39)*two-r2(5,29)*four &
&                  +r1(1,5)*two+r1(1,17)*two)*qx+(+r7(10,4)*two-r6(6,6)*four+r5(3,6)*two &
&                  +r5(3,12)*two+r5(10,12)*two+r5(10,22)*two-r4(6,11)*four-r4(6,31)*four &
&                  +r3(3,4)*two+r3(3,10)*two+r3(3,28)*two+r3(3,40)*two+r3(10,40)*two &
&                  -r2(3,30)*four+r1(3,6)*two+r1(3,18)*two)*qz+(+r6(15,10)-r5(10,14)*two &
&                  +r4(6,15)+r4(6,27)+r4(15,27)+r4(6,39)-r3(10,20)*two-r3(3,44)*two+r2(3,6) &
&                  +r2(3,18)+r2(1,38)+r2(1,53)+r2(3,53)-r1(3,32)*two+r0(7)+r0(22))*xx &
&                  +rxyz(20)*xz+(+r6(6,12)-r5(3,16)*two+r4(1,17)+r4(1,29)+r4(6,29)+r4(6,41) &
&                  -r3(3,22)*two-r3(3,46)*two+r2(1,8)+r2(1,20)+r2(1,40)+r2(1,55)+r2(3,55) &
&                  -r1(3,34)*two+r0(9)+r0(24))*zz+(+r5(10,23)*two-r4(6,32)*four+r3(3,29)*two &
&                  +r3(3,41)*two+r3(10,41)*two-r2(3,31)*four+r1(3,7)*two+r1(3,19)*two)*xxz+( &
&                  +r5(6,24)*two-r4(3,33)*four+r3(1,30)*two+r3(1,42)*two+r3(6,42)*two &
&                  -r2(5,32)*four+r1(1,8)*two+r1(1,20)*two)*xzz+rxyz(1)*xxzz
      eri(2,3,3,1)=r022+(+r7(19,4)*two-r6(13,6)*four+r5(8,6)*two+r5(8,12)*two &
&                  +r5(19,12)*two+r5(10,22)*two-r4(13,11)*four-r4(6,31)*four+r3(8,4)*two &
&                  +r3(8,10)*two+r3(3,28)*two+r3(3,40)*two+r3(10,40)*two-r2(3,30)*four &
&                  +r1(3,6)*two+r1(3,18)*two)*qz+rxyz(5)*zz
      eri(3,3,3,1)=r004+(+r7(21,3)*two+r7(21,4)*two-r6(15,5)*four-r6(15,6)*four &
&                  +r5(10,5)*two+r5(10,6)*two+r5(10,11)*two+r5(21,11)*two+r5(10,12)*two &
&                  +r5(21,12)*two+r5(10,21)*six+r5(10,22)*six-r4(15,10)*four-r4(15,11)*four &
&                  -r4(6,30)*p12-r4(6,31)*p12+r3(10,3)*two+r3(10,4)*two+r3(10,9)*two &
&                  +r3(10,10)*two+r3(3,27)*six+r3(3,28)*six+r3(3,39)*six+r3(10,39)*six &
&                  +r3(3,40)*six+r3(10,40)*six-r2(3,29)*p12-r2(3,30)*p12+r1(3,5)*six &
&                  +r1(3,6)*six+r1(3,17)*six+r1(3,18)*six)*qz+(+r6(15,10)+r6(15,11)*four &
&                  +r6(15,12)-r5(10,14)*two-r5(10,15)*eight-r5(10,16)*two+r4(6,15) &
&                  +r4(6,16)*four+r4(6,17)+r4(6,27)+r4(15,27)+r4(6,28)*four+r4(15,28)*four &
&                  +r4(6,29)+r4(15,29)+r4(6,39)+r4(6,40)*four+r4(6,41)-r3(10,20)*two &
&                  -r3(10,21)*eight-r3(10,22)*two-r3(3,44)*two-r3(3,45)*eight-r3(3,46)*two &
&                  +r2(3,6)+r2(3,7)*four+r2(3,8)+r2(3,18)+r2(3,19)*four+r2(3,20)+r2(1,38) &
&                  +r2(1,39)*four+r2(1,40)+r2(1,53)+r2(3,53)+r2(1,54)*four+r2(3,54)*four &
&                  +r2(1,55)+r2(3,55)-r1(3,32)*two-r1(3,33)*eight-r1(3,34)*two+r0(7) &
&                  +r0(8)*four+r0(9)+r0(22)+r0(23)*four+r0(24))*zz+(+r5(10,23)*two &
&                  +r5(10,24)*two-r4(6,32)*four-r4(6,33)*four+r3(3,29)*two+r3(3,30)*two &
&                  +r3(3,41)*two+r3(10,41)*two+r3(3,42)*two+r3(10,42)*two-r2(3,31)*four &
&                  -r2(3,32)*four+r1(3,7)*two+r1(3,8)*two+r1(3,19)*two+r1(3,20)*two)*zzz &
&                  +rxyz(1)*zzzz
      eri(4,3,3,1)=r112+rxyz(18)*qx+(+r7(14,4)*two-r6(9,6)*four+r5(5,6)*two+r5(5,12)*two &
&                  +r5(14,12)*two-r4(9,11)*four+r3(5,4)*two+r3(5,10)*two)*qz+(+r6(14,11)*two &
&                  -r5(9,15)*four+r4(5,16)*two+r4(5,28)*two+r4(14,28)*two-r3(9,21)*four &
&                  +r2(6,7)*two+r2(6,19)*two)*xz+(+r6(9,12)-r5(5,16)*two+r4(2,17)+r4(2,29) &
&                  +r4(9,29)-r3(5,22)*two+r2(4,8)+r2(4,20))*zz+rxyz(2)*xzz
      eri(5,3,3,1)=r103+(+r7(21,3)-r6(15,5)*two+r5(10,5)+r5(10,11)+r5(21,11) &
&                  +r5(10,21)*three-r4(15,10)*two-r4(6,30)*six+r3(10,3)+r3(10,9) &
&                  +r3(3,27)*three+r3(3,39)*three+r3(10,39)*three-r2(3,29)*six+r1(3,5)*three &
&                  +r1(3,17)*three)*qx+(+r7(15,3)+r7(15,4)*two-r6(10,5)*two-r6(10,6)*four &
&                  +r5(6,5)+r5(6,6)*two+r5(6,11)+r5(15,11)+r5(6,12)*two+r5(15,12)*two &
&                  +r5(6,21)+r5(6,22)*two-r4(10,10)*two-r4(10,11)*four-r4(3,30)*two &
&                  -r4(3,31)*four+r3(6,3)+r3(6,4)*two+r3(6,9)+r3(6,10)*two+r3(1,27) &
&                  +r3(1,28)*two+r3(1,39)+r3(6,39)+r3(1,40)*two+r3(6,40)*two-r2(5,29)*two &
&                  -r2(5,30)*four+r1(1,5)+r1(1,6)*two+r1(1,17)+r1(1,18)*two)*qz+(+r6(15,10) &
&                  +r6(15,11)*two-r5(10,14)*two-r5(10,15)*four+r4(6,15)+r4(6,16)*two+r4(6,27) &
&                  +r4(15,27)+r4(6,28)*two+r4(15,28)*two+r4(6,39)+r4(6,40)*two-r3(10,20)*two &
&                  -r3(10,21)*four-r3(3,44)*two-r3(3,45)*four+r2(3,6)+r2(3,7)*two+r2(3,18) &
&                  +r2(3,19)*two+r2(1,38)+r2(1,39)*two+r2(1,53)+r2(3,53)+r2(1,54)*two &
&                  +r2(3,54)*two-r1(3,32)*two-r1(3,33)*four+r0(7)+r0(8)*two+r0(22)+r0(23)*two &
&                  )*xz+(+r6(10,11)*two+r6(10,12)-r5(6,15)*four-r5(6,16)*two+r4(3,16)*two &
&                  +r4(3,17)+r4(3,28)*two+r4(10,28)*two+r4(3,29)+r4(10,29)-r3(6,21)*four &
&                  -r3(6,22)*two+r2(5,7)*two+r2(5,8)+r2(5,19)*two+r2(5,20))*zz+( &
&                  +r5(10,23)*two+r5(10,24)-r4(6,32)*four-r4(6,33)*two+r3(3,29)*two+r3(3,30) &
&                  +r3(3,41)*two+r3(10,41)*two+r3(3,42)+r3(10,42)-r2(3,31)*four-r2(3,32)*two &
&                  +r1(3,7)*two+r1(3,8)+r1(3,19)*two+r1(3,20))*xzz+(+r5(6,24)-r4(3,33)*two &
&                  +r3(1,30)+r3(1,42)+r3(6,42)-r2(5,32)*two+r1(1,8)+r1(1,20))*zzz+rxyz(1) &
&                  *xzzz
      eri(6,3,3,1)=r013+(+r7(20,3)+r7(20,4)*two-r6(14,5)*two-r6(14,6)*four+r5(9,5) &
&                  +r5(9,6)*two+r5(9,11)+r5(20,11)+r5(9,12)*two+r5(20,12)*two+r5(9,21) &
&                  +r5(9,22)*two-r4(14,10)*two-r4(14,11)*four-r4(5,30)*two-r4(5,31)*four &
&                  +r3(9,3)+r3(9,4)*two+r3(9,9)+r3(9,10)*two+r3(2,27)+r3(2,28)*two+r3(2,39) &
&                  +r3(9,39)+r3(2,40)*two+r3(9,40)*two-r2(6,29)*two-r2(6,30)*four+r1(2,5) &
&                  +r1(2,6)*two+r1(2,17)+r1(2,18)*two)*qz+(+r6(14,11)*two+r6(14,12) &
&                  -r5(9,15)*four-r5(9,16)*two+r4(5,16)*two+r4(5,17)+r4(5,28)*two &
&                  +r4(14,28)*two+r4(5,29)+r4(14,29)-r3(9,21)*four-r3(9,22)*two+r2(6,7)*two &
&                  +r2(6,8)+r2(6,19)*two+r2(6,20))*zz+rxyz(2)*zzz
      eri(1,4,3,1)=r310+(+r7(9,3)*two+r7(9,4)-r6(5,5)*four-r6(5,6)*two+r5(2,5)*two &
&                  +r5(2,6)+r5(2,11)*two+r5(9,11)*two+r5(2,12)+r5(9,12)+r5(9,21)*two+r5(9,22) &
&                  -r4(5,10)*four-r4(5,11)*two-r4(5,30)*four-r4(5,31)*two+r3(2,3)*two+r3(2,4) &
&                  +r3(2,9)*two+r3(2,10)+r3(2,27)*two+r3(2,28)+r3(2,39)*two+r3(9,39)*two &
&                  +r3(2,40)+r3(9,40)-r2(6,29)*four-r2(6,30)*two+r1(2,5)*two+r1(2,6) &
&                  +r1(2,17)*two+r1(2,18))*qx+(+r6(9,10)+r6(9,11)*two-r5(5,14)*two &
&                  -r5(5,15)*four+r4(2,15)+r4(2,16)*two+r4(2,27)+r4(9,27)+r4(2,28)*two &
&                  +r4(9,28)*two-r3(5,20)*two-r3(5,21)*four+r2(4,6)+r2(4,7)*two+r2(4,18) &
&                  +r2(4,19)*two)*xx+rxyz(3)*xxx
      eri(2,4,3,1)=r130+rxyz(8)*qx
      eri(3,4,3,1)=r112+rxyz(19)*qx+(+r7(14,3)*two-r6(9,5)*four+r5(5,5)*two+r5(5,11)*two &
&                  +r5(14,11)*two-r4(9,10)*four+r3(5,3)*two+r3(5,9)*two)*qz+(+r6(14,11)*two &
&                  -r5(9,15)*four+r4(5,16)*two+r4(5,28)*two+r4(14,28)*two-r3(9,21)*four &
&                  +r2(6,7)*two+r2(6,19)*two)*xz+(+r6(9,10)-r5(5,14)*two+r4(2,15)+r4(2,27) &
&                  +r4(9,27)-r3(5,20)*two+r2(4,6)+r2(4,18))*zz+rxyz(3)*xzz
      eri(4,4,3,1)=r220+(+r7(13,3)+r7(13,4)-r6(8,5)*two-r6(8,6)*two+r5(4,5)+r5(4,6) &
&                  +r5(4,11)+r5(13,11)+r5(4,12)+r5(13,12)+r5(6,21)+r5(6,22)-r4(8,10)*two &
&                  -r4(8,11)*two-r4(3,30)*two-r4(3,31)*two+r3(4,3)+r3(4,4)+r3(4,9)+r3(4,10) &
&                  +r3(1,27)+r3(1,28)+r3(1,39)+r3(6,39)+r3(1,40)+r3(6,40)-r2(5,29)*two &
&                  -r2(5,30)*two+r1(1,5)+r1(1,6)+r1(1,17)+r1(1,18))*qx+rxyz(6)*xx
      eri(5,4,3,1)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(9,10)+r6(9,11) &
&                  -r5(5,14)*two-r5(5,15)*two+r4(2,15)+r4(2,16)+r4(2,27)+r4(9,27)+r4(2,28) &
&                  +r4(9,28)-r3(5,20)*two-r3(5,21)*two+r2(4,6)+r2(4,7)+r2(4,18)+r2(4,19))*xz &
&                  +rxyz(3)*xxz
      eri(6,4,3,1)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,3,1)=r301+(+r7(10,3)*two+r7(10,4)-r6(6,5)*four-r6(6,6)*two+r5(3,5)*two &
&                  +r5(3,6)+r5(3,11)*two+r5(10,11)*two+r5(3,12)+r5(10,12)+r5(10,21)*two &
&                  +r5(10,22)-r4(6,10)*four-r4(6,11)*two-r4(6,30)*four-r4(6,31)*two &
&                  +r3(3,3)*two+r3(3,4)+r3(3,9)*two+r3(3,10)+r3(3,27)*two+r3(3,28) &
&                  +r3(3,39)*two+r3(10,39)*two+r3(3,40)+r3(10,40)-r2(3,29)*four-r2(3,30)*two &
&                  +r1(3,5)*two+r1(3,6)+r1(3,17)*two+r1(3,18))*qx+(+r7(6,4)-r6(3,6)*two &
&                  +r5(1,6)+r5(1,12)+r5(6,12)+r5(6,22)*three-r4(3,11)*two-r4(3,31)*six &
&                  +r3(1,4)+r3(1,10)+r3(1,28)*three+r3(1,40)*three+r3(6,40)*three &
&                  -r2(5,30)*six+r1(1,6)*three+r1(1,18)*three)*qz+(+r6(10,10)+r6(10,11)*two &
&                  -r5(6,14)*two-r5(6,15)*four+r4(3,15)+r4(3,16)*two+r4(3,27)+r4(10,27) &
&                  +r4(3,28)*two+r4(10,28)*two-r3(6,20)*two-r3(6,21)*four+r2(5,6)+r2(5,7)*two &
&                  +r2(5,18)+r2(5,19)*two)*xx+(+r6(6,11)*two+r6(6,12)-r5(3,15)*four &
&                  -r5(3,16)*two+r4(1,16)*two+r4(1,17)+r4(1,28)*two+r4(6,28)*two+r4(1,29) &
&                  +r4(6,29)+r4(6,40)*two+r4(6,41)-r3(3,21)*four-r3(3,22)*two-r3(3,45)*four &
&                  -r3(3,46)*two+r2(1,7)*two+r2(1,8)+r2(1,19)*two+r2(1,20)+r2(1,39)*two &
&                  +r2(1,40)+r2(1,54)*two+r2(3,54)*two+r2(1,55)+r2(3,55)-r1(3,33)*four &
&                  -r1(3,34)*two+r0(8)*two+r0(9)+r0(23)*two+r0(24))*xz+(+r5(10,23) &
&                  -r4(6,32)*two+r3(3,29)+r3(3,41)+r3(10,41)-r2(3,31)*two+r1(3,7)+r1(3,19)) &
&                  *xxx+(+r5(6,23)+r5(6,24)*two-r4(3,32)*two-r4(3,33)*four+r3(1,29) &
&                  +r3(1,30)*two+r3(1,41)+r3(6,41)+r3(1,42)*two+r3(6,42)*two-r2(5,31)*two &
&                  -r2(5,32)*four+r1(1,7)+r1(1,8)*two+r1(1,19)+r1(1,20)*two)*xxz+rxyz(1)*xxxz
      eri(2,5,3,1)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,3,1)=r103+(+r7(21,4)-r6(15,6)*two+r5(10,6)+r5(10,12)+r5(21,12) &
&                  +r5(10,22)*three-r4(15,11)*two-r4(6,31)*six+r3(10,4)+r3(10,10) &
&                  +r3(3,28)*three+r3(3,40)*three+r3(10,40)*three-r2(3,30)*six+r1(3,6)*three &
&                  +r1(3,18)*three)*qx+(+r7(15,3)*two+r7(15,4)-r6(10,5)*four-r6(10,6)*two &
&                  +r5(6,5)*two+r5(6,6)+r5(6,11)*two+r5(15,11)*two+r5(6,12)+r5(15,12) &
&                  +r5(6,21)*two+r5(6,22)-r4(10,10)*four-r4(10,11)*two-r4(3,30)*four &
&                  -r4(3,31)*two+r3(6,3)*two+r3(6,4)+r3(6,9)*two+r3(6,10)+r3(1,27)*two &
&                  +r3(1,28)+r3(1,39)*two+r3(6,39)*two+r3(1,40)+r3(6,40)-r2(5,29)*four &
&                  -r2(5,30)*two+r1(1,5)*two+r1(1,6)+r1(1,17)*two+r1(1,18))*qz+( &
&                  +r6(15,11)*two+r6(15,12)-r5(10,15)*four-r5(10,16)*two+r4(6,16)*two &
&                  +r4(6,17)+r4(6,28)*two+r4(15,28)*two+r4(6,29)+r4(15,29)+r4(6,40)*two &
&                  +r4(6,41)-r3(10,21)*four-r3(10,22)*two-r3(3,45)*four-r3(3,46)*two &
&                  +r2(3,7)*two+r2(3,8)+r2(3,19)*two+r2(3,20)+r2(1,39)*two+r2(1,40) &
&                  +r2(1,54)*two+r2(3,54)*two+r2(1,55)+r2(3,55)-r1(3,33)*four-r1(3,34)*two &
&                  +r0(8)*two+r0(9)+r0(23)*two+r0(24))*xz+(+r6(10,10)+r6(10,11)*two &
&                  -r5(6,14)*two-r5(6,15)*four+r4(3,15)+r4(3,16)*two+r4(3,27)+r4(10,27) &
&                  +r4(3,28)*two+r4(10,28)*two-r3(6,20)*two-r3(6,21)*four+r2(5,6)+r2(5,7)*two &
&                  +r2(5,18)+r2(5,19)*two)*zz+(+r5(10,23)+r5(10,24)*two-r4(6,32)*two &
&                  -r4(6,33)*four+r3(3,29)+r3(3,30)*two+r3(3,41)+r3(10,41)+r3(3,42)*two &
&                  +r3(10,42)*two-r2(3,31)*two-r2(3,32)*four+r1(3,7)+r1(3,8)*two+r1(3,19) &
&                  +r1(3,20)*two)*xzz+(+r5(6,23)-r4(3,32)*two+r3(1,29)+r3(1,41)+r3(6,41) &
&                  -r2(5,31)*two+r1(1,7)+r1(1,19))*zzz+rxyz(1)*xzzz
      eri(4,5,3,1)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(9,11)+r6(9,12) &
&                  -r5(5,15)*two-r5(5,16)*two+r4(2,16)+r4(2,17)+r4(2,28)+r4(9,28)+r4(2,29) &
&                  +r4(9,29)-r3(5,21)*two-r3(5,22)*two+r2(4,7)+r2(4,8)+r2(4,19)+r2(4,20))*xz &
&                  +rxyz(2)*xxz
      eri(5,5,3,1)=r202+(+r7(15,3)+r7(15,4)-r6(10,5)*two-r6(10,6)*two+r5(6,5)+r5(6,6) &
&                  +r5(6,11)+r5(15,11)+r5(6,12)+r5(15,12)+r5(6,21)+r5(6,22)-r4(10,10)*two &
&                  -r4(10,11)*two-r4(3,30)*two-r4(3,31)*two+r3(6,3)+r3(6,4)+r3(6,9)+r3(6,10) &
&                  +r3(1,27)+r3(1,28)+r3(1,39)+r3(6,39)+r3(1,40)+r3(6,40)-r2(5,29)*two &
&                  -r2(5,30)*two+r1(1,5)+r1(1,6)+r1(1,17)+r1(1,18))*qx+(+r7(10,3)+r7(10,4) &
&                  -r6(6,5)*two-r6(6,6)*two+r5(3,5)+r5(3,6)+r5(3,11)+r5(10,11)+r5(3,12) &
&                  +r5(10,12)+r5(10,21)+r5(10,22)-r4(6,10)*two-r4(6,11)*two-r4(6,30)*two &
&                  -r4(6,31)*two+r3(3,3)+r3(3,4)+r3(3,9)+r3(3,10)+r3(3,27)+r3(3,28)+r3(3,39) &
&                  +r3(10,39)+r3(3,40)+r3(10,40)-r2(3,29)*two-r2(3,30)*two+r1(3,5)+r1(3,6) &
&                  +r1(3,17)+r1(3,18))*qz+(+r6(15,11)-r5(10,15)*two+r4(6,16)+r4(6,28) &
&                  +r4(15,28)+r4(6,40)-r3(10,21)*two-r3(3,45)*two+r2(3,7)+r2(3,19)+r2(1,39) &
&                  +r2(1,54)+r2(3,54)-r1(3,33)*two+r0(8)+r0(23))*xx+(+r6(10,10)+r6(10,11)*two &
&                  +r6(10,12)-r5(6,14)*two-r5(6,15)*four-r5(6,16)*two+r4(3,15)+r4(3,16)*two &
&                  +r4(3,17)+r4(3,27)+r4(10,27)+r4(3,28)*two+r4(10,28)*two+r4(3,29)+r4(10,29) &
&                  -r3(6,20)*two-r3(6,21)*four-r3(6,22)*two+r2(5,6)+r2(5,7)*two+r2(5,8) &
&                  +r2(5,18)+r2(5,19)*two+r2(5,20))*xz+(+r6(6,11)-r5(3,15)*two+r4(1,16) &
&                  +r4(1,28)+r4(6,28)+r4(6,40)-r3(3,21)*two-r3(3,45)*two+r2(1,7)+r2(1,19) &
&                  +r2(1,39)+r2(1,54)+r2(3,54)-r1(3,33)*two+r0(8)+r0(23))*zz+(+r5(10,23) &
&                  +r5(10,24)-r4(6,32)*two-r4(6,33)*two+r3(3,29)+r3(3,30)+r3(3,41)+r3(10,41) &
&                  +r3(3,42)+r3(10,42)-r2(3,31)*two-r2(3,32)*two+r1(3,7)+r1(3,8)+r1(3,19) &
&                  +r1(3,20))*xxz+(+r5(6,23)+r5(6,24)-r4(3,32)*two-r4(3,33)*two+r3(1,29) &
&                  +r3(1,30)+r3(1,41)+r3(6,41)+r3(1,42)+r3(6,42)-r2(5,31)*two-r2(5,32)*two &
&                  +r1(1,7)+r1(1,8)+r1(1,19)+r1(1,20))*xzz+rxyz(1)*xxzz
      eri(6,5,3,1)=r112+rxyz(19)*qx+(+r7(14,3)+r7(14,4)-r6(9,5)*two-r6(9,6)*two+r5(5,5) &
&                  +r5(5,6)+r5(5,11)+r5(14,11)+r5(5,12)+r5(14,12)-r4(9,10)*two-r4(9,11)*two &
&                  +r3(5,3)+r3(5,4)+r3(5,9)+r3(5,10))*qz+(+r6(14,11)+r6(14,12)-r5(9,15)*two &
&                  -r5(9,16)*two+r4(5,16)+r4(5,17)+r4(5,28)+r4(14,28)+r4(5,29)+r4(14,29) &
&                  -r3(9,21)*two-r3(9,22)*two+r2(6,7)+r2(6,8)+r2(6,19)+r2(6,20))*xz+rxyz(11) &
&                  *zz+rxyz(2)*xzz
      eri(1,6,3,1)=r211+(+r7(14,3)*two-r6(9,5)*four+r5(5,5)*two+r5(5,11)*two &
&                  +r5(14,11)*two-r4(9,10)*four+r3(5,3)*two+r3(5,9)*two)*qx+rxyz(16)*qz+( &
&                  +r6(14,10)-r5(9,14)*two+r4(5,15)+r4(5,27)+r4(14,27)-r3(9,20)*two+r2(6,6) &
&                  +r2(6,18))*xx+(+r6(9,11)*two-r5(5,15)*four+r4(2,16)*two+r4(2,28)*two &
&                  +r4(9,28)*two-r3(5,21)*four+r2(4,7)*two+r2(4,19)*two)*xz+rxyz(3)*xxz
      eri(2,6,3,1)=r031+rxyz(8)*qz
      eri(3,6,3,1)=r013+(+r7(20,3)*two+r7(20,4)-r6(14,5)*four-r6(14,6)*two+r5(9,5)*two &
&                  +r5(9,6)+r5(9,11)*two+r5(20,11)*two+r5(9,12)+r5(20,12)+r5(9,21)*two &
&                  +r5(9,22)-r4(14,10)*four-r4(14,11)*two-r4(5,30)*four-r4(5,31)*two &
&                  +r3(9,3)*two+r3(9,4)+r3(9,9)*two+r3(9,10)+r3(2,27)*two+r3(2,28) &
&                  +r3(2,39)*two+r3(9,39)*two+r3(2,40)+r3(9,40)-r2(6,29)*four-r2(6,30)*two &
&                  +r1(2,5)*two+r1(2,6)+r1(2,17)*two+r1(2,18))*qz+(+r6(14,10)+r6(14,11)*two &
&                  -r5(9,14)*two-r5(9,15)*four+r4(5,15)+r4(5,16)*two+r4(5,27)+r4(14,27) &
&                  +r4(5,28)*two+r4(14,28)*two-r3(9,20)*two-r3(9,21)*four+r2(6,6)+r2(6,7)*two &
&                  +r2(6,18)+r2(6,19)*two)*zz+rxyz(3)*zzz
      eri(4,6,3,1)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,3,1)=r112+rxyz(18)*qx+(+r7(14,3)+r7(14,4)-r6(9,5)*two-r6(9,6)*two+r5(5,5) &
&                  +r5(5,6)+r5(5,11)+r5(14,11)+r5(5,12)+r5(14,12)-r4(9,10)*two-r4(9,11)*two &
&                  +r3(5,3)+r3(5,4)+r3(5,9)+r3(5,10))*qz+(+r6(14,10)+r6(14,11)-r5(9,14)*two &
&                  -r5(9,15)*two+r4(5,15)+r4(5,16)+r4(5,27)+r4(14,27)+r4(5,28)+r4(14,28) &
&                  -r3(9,20)*two-r3(9,21)*two+r2(6,6)+r2(6,7)+r2(6,18)+r2(6,19))*xz+rxyz(11) &
&                  *zz+rxyz(3)*xzz
      eri(6,6,3,1)=r022+(+r7(19,3)+r7(19,4)-r6(13,5)*two-r6(13,6)*two+r5(8,5)+r5(8,6) &
&                  +r5(8,11)+r5(19,11)+r5(8,12)+r5(19,12)+r5(10,21)+r5(10,22)-r4(13,10)*two &
&                  -r4(13,11)*two-r4(6,30)*two-r4(6,31)*two+r3(8,3)+r3(8,4)+r3(8,9)+r3(8,10) &
&                  +r3(3,27)+r3(3,28)+r3(3,39)+r3(10,39)+r3(3,40)+r3(10,40)-r2(3,29)*two &
&                  -r2(3,30)*two+r1(3,5)+r1(3,6)+r1(3,17)+r1(3,18))*qz+rxyz(6)*zz
!
      r400= r8(2)+r6(2,4)*three+r6(2,9)*six+r4(2,26)*p18+r4(2,38)*three+r2(4,52)*nine
      r310= r8(4)+r6(4,4)*three+r6(4,9)*three+r4(4,26)*nine
      r301= r8(5)+r6(5,4)*three+r6(5,9)*three+r4(5,26)*nine
      r220= r8(7)+r6(7,4)*three+r6(2,9)+r6(7,9)+r4(2,26)*three+r4(7,26)*three+r4(2,38) &
&          +r2(4,52)*three
      r211= r8(8)+r6(8,4)*three+r6(8,9)+r4(8,26)*three
      r202= r8(9)+r6(9,4)*three+r6(2,9)+r6(9,9)+r4(2,26)*three+r4(9,26)*three+r4(2,38) &
&          +r2(4,52)*three
      r130= r8(11)+r6(11,4)*three+r6(4,9)*three+r4(4,26)*nine
      r121= r8(12)+r6(12,4)*three+r6(5,9)+r4(5,26)*three
      r112= r8(13)+r6(13,4)*three+r6(4,9)+r4(4,26)*three
      r103= r8(14)+r6(14,4)*three+r6(5,9)*three+r4(5,26)*nine
      r040= r8(16)+r6(16,4)*three+r6(7,9)*six+r4(7,26)*p18+r4(2,38)*three+r2(4,52)*nine
      r031= r8(17)+r6(17,4)*three+r6(8,9)*three+r4(8,26)*nine
      r022= r8(18)+r6(18,4)*three+r6(7,9)+r6(9,9)+r4(7,26)*three+r4(9,26)*three+r4(2,38) &
&          +r2(4,52)*three
      r013= r8(19)+r6(19,4)*three+r6(8,9)*three+r4(8,26)*nine
      r004= r8(20)+r6(20,4)*three+r6(9,9)*six+r4(9,26)*p18+r4(2,38)*three+r2(4,52)*nine
      rxyz(1)=+r4(2,42)+r2(4,56)*three
      rxyz(2)=+r5(4,24)+r3(4,42)*three
      rxyz(3)=+r5(4,23)+r3(4,41)*three
      rxyz(4)=+r6(7,10)+r4(7,27)*three+r4(2,39)+r2(4,53)*three
      rxyz(5)=+r6(7,12)+r4(7,29)*three+r4(2,41)+r2(4,55)*three
      rxyz(6)=+r6(7,11)+r4(7,28)*three+r4(2,40)+r2(4,54)*three
      rxyz(7)=+r7(11,3)+r5(11,11)*three+r5(4,21)*three+r3(4,39)*nine
      rxyz(8)=+r7(11,4)+r5(11,12)*three+r5(4,22)*three+r3(4,40)*nine
      rxyz(9)=+r7(8,3)+r7(8,4)+r5(8,11)*three+r5(8,12)*three
      rxyz(10)=+r6(8,11)+r4(8,28)*three
      rxyz(11)=+r6(4,11)+r4(4,28)*three
      rxyz(12)=+r7(12,3)+r5(12,11)*three+r5(5,21)+r3(5,39)*three
      rxyz(13)=+r7(12,4)+r5(12,12)*three+r5(5,22)+r3(5,40)*three
      rxyz(14)=+r7(7,4)+r5(7,12)*three+r5(2,22)+r3(2,40)*three
      rxyz(15)=+r7(7,3)+r5(7,11)*three+r5(2,21)+r3(2,39)*three
      rxyz(16)=+r7(4,4)+r5(4,12)*three+r5(4,22)+r3(4,40)*three
      rxyz(17)=+r7(4,3)+r5(4,11)*three+r5(4,21)+r3(4,39)*three
      rxyz(18)=+r7(13,3)+r5(13,11)*three+r5(4,21)+r3(4,39)*three
      rxyz(19)=+r7(13,4)+r5(13,12)*three+r5(4,22)+r3(4,40)*three
      rxyz(20)=+r6(5,11)*four+r4(5,28)*p12
      eri(1,1,4,1)=r400+(+r7(2,3)*two+r7(2,4)*two+r5(2,11)*six+r5(2,12)*six+r5(2,21)*six &
&                  +r5(2,22)*six+r3(2,39)*p18+r3(2,40)*p18)*qx+(+r6(2,10)+r6(2,11)*four &
&                  +r6(2,12)+r4(2,27)*three+r4(2,28)*p12+r4(2,29)*three+r4(2,39) &
&                  +r4(2,40)*four+r4(2,41)+r2(4,53)*three+r2(4,54)*p12+r2(4,55)*three)*xx+( &
&                  +r5(2,23)*two+r5(2,24)*two+r3(2,41)*six+r3(2,42)*six)*xxx+rxyz(1)*xxxx
      eri(2,1,4,1)=r220+(+r7(7,4)*two+r5(7,12)*six+r5(2,22)*two+r3(2,40)*six)*qx+rxyz(5) &
&                  *xx
      eri(3,1,4,1)=r202+(+r7(9,4)*two+r5(9,12)*six+r5(2,22)*two+r3(2,40)*six)*qx+( &
&                  +r7(5,3)*two+r5(5,11)*six+r5(5,21)*two+r3(5,39)*six)*qz+(+r6(9,12) &
&                  +r4(9,29)*three+r4(2,41)+r2(4,55)*three)*xx+rxyz(20)*xz+(+r6(2,10) &
&                  +r4(2,27)*three+r4(2,39)+r2(4,53)*three)*zz+(+r5(5,24)*two+r3(5,42)*six) &
&                  *xxz+(+r5(2,23)*two+r3(2,41)*six)*xzz+rxyz(1)*xxzz
      eri(4,1,4,1)=r310+(+r7(4,3)+r7(4,4)*two+r5(4,11)*three+r5(4,12)*six+r5(4,21) &
&                  +r5(4,22)*two+r3(4,39)*three+r3(4,40)*six)*qx+(+r6(4,11)*two+r6(4,12) &
&                  +r4(4,28)*six+r4(4,29)*three)*xx+rxyz(2)*xxx
      eri(5,1,4,1)=r301+(+r7(5,3)+r7(5,4)*two+r5(5,11)*three+r5(5,12)*six+r5(5,21) &
&                  +r5(5,22)*two+r3(5,39)*three+r3(5,40)*six)*qx+(+r7(2,3)+r5(2,11)*three &
&                  +r5(2,21)*three+r3(2,39)*nine)*qz+(+r6(5,11)*two+r6(5,12)+r4(5,28)*six &
&                  +r4(5,29)*three)*xx+(+r6(2,10)+r6(2,11)*two+r4(2,27)*three+r4(2,28)*six &
&                  +r4(2,39)+r4(2,40)*two+r2(4,53)*three+r2(4,54)*six)*xz+(+r5(5,24) &
&                  +r3(5,42)*three)*xxx+(+r5(2,23)*two+r5(2,24)+r3(2,41)*six+r3(2,42)*three) &
&                  *xxz+rxyz(1)*xxxz
      eri(6,1,4,1)=r211+(+r7(8,4)*two+r5(8,12)*six)*qx+rxyz(17)*qz+(+r6(8,12) &
&                  +r4(8,29)*three)*xx+(+r6(4,11)*two+r4(4,28)*six)*xz+rxyz(2)*xxz
      eri(1,2,4,1)=r220+(+r7(7,3)*two+r5(7,11)*six+r5(2,21)*two+r3(2,39)*six)*qx+rxyz(4) &
&                  *xx
      eri(2,2,4,1)=r040
      eri(3,2,4,1)=r022+(+r7(12,3)*two+r5(12,11)*six+r5(5,21)*two+r3(5,39)*six)*qz &
&                  +rxyz(4)*zz
      eri(4,2,4,1)=r130+rxyz(7)*qx
      eri(5,2,4,1)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,4,1)=r031+rxyz(7)*qz
      eri(1,3,4,1)=r202+(+r7(9,3)*two+r5(9,11)*six+r5(2,21)*two+r3(2,39)*six)*qx+( &
&                  +r7(5,4)*two+r5(5,12)*six+r5(5,22)*two+r3(5,40)*six)*qz+(+r6(9,10) &
&                  +r4(9,27)*three+r4(2,39)+r2(4,53)*three)*xx+rxyz(20)*xz+(+r6(2,12) &
&                  +r4(2,29)*three+r4(2,41)+r2(4,55)*three)*zz+(+r5(5,23)*two+r3(5,41)*six) &
&                  *xxz+(+r5(2,24)*two+r3(2,42)*six)*xzz+rxyz(1)*xxzz
      eri(2,3,4,1)=r022+(+r7(12,4)*two+r5(12,12)*six+r5(5,22)*two+r3(5,40)*six)*qz &
&                  +rxyz(5)*zz
      eri(3,3,4,1)=r004+(+r7(14,3)*two+r7(14,4)*two+r5(14,11)*six+r5(14,12)*six &
&                  +r5(5,21)*six+r5(5,22)*six+r3(5,39)*p18+r3(5,40)*p18)*qz+(+r6(9,10) &
&                  +r6(9,11)*four+r6(9,12)+r4(9,27)*three+r4(9,28)*p12+r4(9,29)*three &
&                  +r4(2,39)+r4(2,40)*four+r4(2,41)+r2(4,53)*three+r2(4,54)*p12 &
&                  +r2(4,55)*three)*zz+(+r5(5,23)*two+r5(5,24)*two+r3(5,41)*six+r3(5,42)*six) &
&                  *zzz+rxyz(1)*zzzz
      eri(4,3,4,1)=r112+rxyz(18)*qx+(+r7(8,4)*two+r5(8,12)*six)*qz+(+r6(8,11)*two &
&                  +r4(8,28)*six)*xz+(+r6(4,12)+r4(4,29)*three)*zz+rxyz(2)*xzz
      eri(5,3,4,1)=r103+(+r7(14,3)+r5(14,11)*three+r5(5,21)*three+r3(5,39)*nine)*qx+( &
&                  +r7(9,3)+r7(9,4)*two+r5(9,11)*three+r5(9,12)*six+r5(2,21)+r5(2,22)*two &
&                  +r3(2,39)*three+r3(2,40)*six)*qz+(+r6(9,10)+r6(9,11)*two+r4(9,27)*three &
&                  +r4(9,28)*six+r4(2,39)+r4(2,40)*two+r2(4,53)*three+r2(4,54)*six)*xz+( &
&                  +r6(5,11)*two+r6(5,12)+r4(5,28)*six+r4(5,29)*three)*zz+(+r5(5,23)*two &
&                  +r5(5,24)+r3(5,41)*six+r3(5,42)*three)*xzz+(+r5(2,24)+r3(2,42)*three)*zzz &
&                  +rxyz(1)*xzzz
      eri(6,3,4,1)=r013+(+r7(13,3)+r7(13,4)*two+r5(13,11)*three+r5(13,12)*six+r5(4,21) &
&                  +r5(4,22)*two+r3(4,39)*three+r3(4,40)*six)*qz+(+r6(8,11)*two+r6(8,12) &
&                  +r4(8,28)*six+r4(8,29)*three)*zz+rxyz(2)*zzz
      eri(1,4,4,1)=r310+(+r7(4,3)*two+r7(4,4)+r5(4,11)*six+r5(4,12)*three+r5(4,21)*two &
&                  +r5(4,22)+r3(4,39)*six+r3(4,40)*three)*qx+(+r6(4,10)+r6(4,11)*two &
&                  +r4(4,27)*three+r4(4,28)*six)*xx+rxyz(3)*xxx
      eri(2,4,4,1)=r130+rxyz(8)*qx
      eri(3,4,4,1)=r112+rxyz(19)*qx+(+r7(8,3)*two+r5(8,11)*six)*qz+(+r6(8,11)*two &
&                  +r4(8,28)*six)*xz+(+r6(4,10)+r4(4,27)*three)*zz+rxyz(3)*xzz
      eri(4,4,4,1)=r220+(+r7(7,3)+r7(7,4)+r5(7,11)*three+r5(7,12)*three+r5(2,21)+r5(2,22) &
&                  +r3(2,39)*three+r3(2,40)*three)*qx+rxyz(6)*xx
      eri(5,4,4,1)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(4,10)+r6(4,11) &
&                  +r4(4,27)*three+r4(4,28)*three)*xz+rxyz(3)*xxz
      eri(6,4,4,1)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,4,1)=r301+(+r7(5,3)*two+r7(5,4)+r5(5,11)*six+r5(5,12)*three+r5(5,21)*two &
&                  +r5(5,22)+r3(5,39)*six+r3(5,40)*three)*qx+(+r7(2,4)+r5(2,12)*three &
&                  +r5(2,22)*three+r3(2,40)*nine)*qz+(+r6(5,10)+r6(5,11)*two+r4(5,27)*three &
&                  +r4(5,28)*six)*xx+(+r6(2,11)*two+r6(2,12)+r4(2,28)*six+r4(2,29)*three &
&                  +r4(2,40)*two+r4(2,41)+r2(4,54)*six+r2(4,55)*three)*xz+(+r5(5,23) &
&                  +r3(5,41)*three)*xxx+(+r5(2,23)+r5(2,24)*two+r3(2,41)*three+r3(2,42)*six) &
&                  *xxz+rxyz(1)*xxxz
      eri(2,5,4,1)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,4,1)=r103+(+r7(14,4)+r5(14,12)*three+r5(5,22)*three+r3(5,40)*nine)*qx+( &
&                  +r7(9,3)*two+r7(9,4)+r5(9,11)*six+r5(9,12)*three+r5(2,21)*two+r5(2,22) &
&                  +r3(2,39)*six+r3(2,40)*three)*qz+(+r6(9,11)*two+r6(9,12)+r4(9,28)*six &
&                  +r4(9,29)*three+r4(2,40)*two+r4(2,41)+r2(4,54)*six+r2(4,55)*three)*xz+( &
&                  +r6(5,10)+r6(5,11)*two+r4(5,27)*three+r4(5,28)*six)*zz+(+r5(5,23) &
&                  +r5(5,24)*two+r3(5,41)*three+r3(5,42)*six)*xzz+(+r5(2,23)+r3(2,41)*three) &
&                  *zzz+rxyz(1)*xzzz
      eri(4,5,4,1)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(4,11)+r6(4,12) &
&                  +r4(4,28)*three+r4(4,29)*three)*xz+rxyz(2)*xxz
      eri(5,5,4,1)=r202+(+r7(9,3)+r7(9,4)+r5(9,11)*three+r5(9,12)*three+r5(2,21)+r5(2,22) &
&                  +r3(2,39)*three+r3(2,40)*three)*qx+(+r7(5,3)+r7(5,4)+r5(5,11)*three &
&                  +r5(5,12)*three+r5(5,21)+r5(5,22)+r3(5,39)*three+r3(5,40)*three)*qz+( &
&                  +r6(9,11)+r4(9,28)*three+r4(2,40)+r2(4,54)*three)*xx+(+r6(5,10) &
&                  +r6(5,11)*two+r6(5,12)+r4(5,27)*three+r4(5,28)*six+r4(5,29)*three)*xz+( &
&                  +r6(2,11)+r4(2,28)*three+r4(2,40)+r2(4,54)*three)*zz+(+r5(5,23)+r5(5,24) &
&                  +r3(5,41)*three+r3(5,42)*three)*xxz+(+r5(2,23)+r5(2,24)+r3(2,41)*three &
&                  +r3(2,42)*three)*xzz+rxyz(1)*xxzz
      eri(6,5,4,1)=r112+rxyz(19)*qx+(+r7(8,3)+r7(8,4)+r5(8,11)*three+r5(8,12)*three)*qz+( &
&                  +r6(8,11)+r6(8,12)+r4(8,28)*three+r4(8,29)*three)*xz+rxyz(11)*zz+rxyz(2) &
&                  *xzz
      eri(1,6,4,1)=r211+(+r7(8,3)*two+r5(8,11)*six)*qx+rxyz(16)*qz+(+r6(8,10) &
&                  +r4(8,27)*three)*xx+(+r6(4,11)*two+r4(4,28)*six)*xz+rxyz(3)*xxz
      eri(2,6,4,1)=r031+rxyz(8)*qz
      eri(3,6,4,1)=r013+(+r7(13,3)*two+r7(13,4)+r5(13,11)*six+r5(13,12)*three &
&                  +r5(4,21)*two+r5(4,22)+r3(4,39)*six+r3(4,40)*three)*qz+(+r6(8,10) &
&                  +r6(8,11)*two+r4(8,27)*three+r4(8,28)*six)*zz+rxyz(3)*zzz
      eri(4,6,4,1)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,4,1)=r112+rxyz(18)*qx+(+r7(8,3)+r7(8,4)+r5(8,11)*three+r5(8,12)*three)*qz+( &
&                  +r6(8,10)+r6(8,11)+r4(8,27)*three+r4(8,28)*three)*xz+rxyz(11)*zz+rxyz(3) &
&                  *xzz
      eri(6,6,4,1)=r022+(+r7(12,3)+r7(12,4)+r5(12,11)*three+r5(12,12)*three+r5(5,21) &
&                  +r5(5,22)+r3(5,39)*three+r3(5,40)*three)*qz+rxyz(6)*zz
!
      r400= r8(3)-r7(1,1)+r6(3,4)*three+r6(3,9)*six-r5(1,3)*three-r5(1,13)*six+r4(3,26)*p18 &
&          +r4(3,38)*three-r3(1,19)*p18-r3(1,43)*three+r2(5,52)*nine-r1(1,31)*nine
      r310= r8(5)-r7(2,1)+r6(5,4)*three+r6(5,9)*three-r5(2,3)*three-r5(2,13)*three &
&          +r4(5,26)*nine-r3(2,19)*nine
      r301= r8(6)-r7(3,1)+r6(6,4)*three+r6(6,9)*three-r5(3,3)*three-r5(3,13)*three &
&          +r4(6,26)*nine-r3(3,19)*nine
      r220= r8(8)-r7(4,1)+r6(8,4)*three+r6(3,9)+r6(8,9)-r5(4,3)*three-r5(1,13)-r5(4,13) &
&          +r4(3,26)*three+r4(8,26)*three+r4(3,38)-r3(1,19)*three-r3(4,19)*three-r3(1,43) &
&          +r2(5,52)*three-r1(1,31)*three
      r211= r8(9)-r7(5,1)+r6(9,4)*three+r6(9,9)-r5(5,3)*three-r5(5,13)+r4(9,26)*three &
&          -r3(5,19)*three
      r202= r8(10)-r7(6,1)+r6(10,4)*three+r6(3,9)+r6(10,9)-r5(6,3)*three-r5(1,13)-r5(6,13) &
&          +r4(3,26)*three+r4(10,26)*three+r4(3,38)-r3(1,19)*three-r3(6,19)*three-r3(1,43) &
&          +r2(5,52)*three-r1(1,31)*three
      r130= r8(12)-r7(7,1)+r6(12,4)*three+r6(5,9)*three-r5(7,3)*three-r5(2,13)*three &
&          +r4(5,26)*nine-r3(2,19)*nine
      r121= r8(13)-r7(8,1)+r6(13,4)*three+r6(6,9)-r5(8,3)*three-r5(3,13)+r4(6,26)*three &
&          -r3(3,19)*three
      r112= r8(14)-r7(9,1)+r6(14,4)*three+r6(5,9)-r5(9,3)*three-r5(2,13)+r4(5,26)*three &
&          -r3(2,19)*three
      r103= r8(15)-r7(10,1)+r6(15,4)*three+r6(6,9)*three-r5(10,3)*three-r5(3,13)*three &
&          +r4(6,26)*nine-r3(3,19)*nine
      r040= r8(17)-r7(11,1)+r6(17,4)*three+r6(8,9)*six-r5(11,3)*three-r5(4,13)*six &
&          +r4(8,26)*p18+r4(3,38)*three-r3(4,19)*p18-r3(1,43)*three+r2(5,52)*nine &
&          -r1(1,31)*nine
      r031= r8(18)-r7(12,1)+r6(18,4)*three+r6(9,9)*three-r5(12,3)*three-r5(5,13)*three &
&          +r4(9,26)*nine-r3(5,19)*nine
      r022= r8(19)-r7(13,1)+r6(19,4)*three+r6(8,9)+r6(10,9)-r5(13,3)*three-r5(4,13) &
&          -r5(6,13)+r4(8,26)*three+r4(10,26)*three+r4(3,38)-r3(4,19)*three-r3(6,19)*three &
&          -r3(1,43)+r2(5,52)*three-r1(1,31)*three
      r013= r8(20)-r7(14,1)+r6(20,4)*three+r6(9,9)*three-r5(14,3)*three-r5(5,13)*three &
&          +r4(9,26)*nine-r3(5,19)*nine
      r004= r8(21)-r7(15,1)+r6(21,4)*three+r6(10,9)*six-r5(15,3)*three-r5(6,13)*six &
&          +r4(10,26)*p18+r4(3,38)*three-r3(6,19)*p18-r3(1,43)*three+r2(5,52)*nine &
&          -r1(1,31)*nine
      rxyz(1)=+r4(3,42)-r3(1,47)+r2(5,56)*three-r1(1,35)*three
      rxyz(2)=+r5(5,24)-r4(2,33)+r3(5,42)*three-r2(4,32)*three
      rxyz(3)=+r5(5,23)-r4(2,32)+r3(5,41)*three-r2(4,31)*three
      rxyz(4)=+r6(8,10)-r5(4,14)+r4(8,27)*three+r4(3,39)-r3(4,20)*three-r3(1,44) &
&             +r2(5,53)*three-r1(1,32)*three
      rxyz(5)=+r6(8,12)-r5(4,16)+r4(8,29)*three+r4(3,41)-r3(4,22)*three-r3(1,46) &
&             +r2(5,55)*three-r1(1,34)*three
      rxyz(6)=+r6(8,11)-r5(4,15)+r4(8,28)*three+r4(3,40)-r3(4,21)*three-r3(1,45) &
&             +r2(5,54)*three-r1(1,33)*three
      rxyz(7)=+r7(12,3)-r6(7,5)+r5(12,11)*three+r5(5,21)*three-r4(7,10)*three &
&             -r4(2,30)*three+r3(5,39)*nine-r2(4,29)*nine
      rxyz(8)=+r7(12,4)-r6(7,6)+r5(12,12)*three+r5(5,22)*three-r4(7,11)*three &
&             -r4(2,31)*three+r3(5,40)*nine-r2(4,30)*nine
      rxyz(9)=+r7(9,3)+r7(9,4)-r6(5,5)-r6(5,6)+r5(9,11)*three+r5(9,12)*three &
&             -r4(5,10)*three-r4(5,11)*three
      rxyz(10)=+r6(9,11)-r5(5,15)+r4(9,28)*three-r3(5,21)*three
      rxyz(11)=+r6(5,11)-r5(2,15)+r4(5,28)*three-r3(2,21)*three
      rxyz(12)=+r7(13,3)-r6(8,5)+r5(13,11)*three+r5(6,21)-r4(8,10)*three-r4(3,30) &
&             +r3(6,39)*three-r2(5,29)*three
      rxyz(13)=+r7(13,4)-r6(8,6)+r5(13,12)*three+r5(6,22)-r4(8,11)*three-r4(3,31) &
&             +r3(6,40)*three-r2(5,30)*three
      rxyz(14)=+r7(8,4)-r6(4,6)+r5(8,12)*three+r5(3,22)-r4(4,11)*three-r4(1,31) &
&             +r3(3,40)*three-r2(1,30)*three
      rxyz(15)=+r7(8,3)-r6(4,5)+r5(8,11)*three+r5(3,21)-r4(4,10)*three-r4(1,30) &
&             +r3(3,39)*three-r2(1,29)*three
      rxyz(16)=+r7(5,4)-r6(2,6)+r5(5,12)*three+r5(5,22)-r4(2,11)*three-r4(2,31) &
&             +r3(5,40)*three-r2(4,30)*three
      rxyz(17)=+r7(5,3)-r6(2,5)+r5(5,11)*three+r5(5,21)-r4(2,10)*three-r4(2,30) &
&             +r3(5,39)*three-r2(4,29)*three
      rxyz(18)=+r7(14,3)-r6(9,5)+r5(14,11)*three+r5(5,21)-r4(9,10)*three-r4(2,30) &
&             +r3(5,39)*three-r2(4,29)*three
      rxyz(19)=+r7(14,4)-r6(9,6)+r5(14,12)*three+r5(5,22)-r4(9,11)*three-r4(2,31) &
&             +r3(5,40)*three-r2(4,30)*three
      rxyz(20)=+r6(6,11)*four-r5(3,15)*four+r4(6,28)*p12-r3(3,21)*p12
      eri(1,1,5,1)=r400+(+r7(3,3)*two+r7(3,4)*two-r6(1,5)*two-r6(1,6)*two+r5(3,11)*six &
&                  +r5(3,12)*six+r5(3,21)*six+r5(3,22)*six-r4(1,10)*six-r4(1,11)*six &
&                  -r4(1,30)*six-r4(1,31)*six+r3(3,39)*p18+r3(3,40)*p18-r2(1,29)*p18 &
&                  -r2(1,30)*p18)*qx+(+r6(3,10)+r6(3,11)*four+r6(3,12)-r5(1,14)-r5(1,15)*four &
&                  -r5(1,16)+r4(3,27)*three+r4(3,28)*p12+r4(3,29)*three+r4(3,39) &
&                  +r4(3,40)*four+r4(3,41)-r3(1,20)*three-r3(1,21)*p12-r3(1,22)*three &
&                  -r3(1,44)-r3(1,45)*four-r3(1,46)+r2(5,53)*three+r2(5,54)*p12 &
&                  +r2(5,55)*three-r1(1,32)*three-r1(1,33)*p12-r1(1,34)*three)*xx+( &
&                  +r5(3,23)*two+r5(3,24)*two-r4(1,32)*two-r4(1,33)*two+r3(3,41)*six &
&                  +r3(3,42)*six-r2(1,31)*six-r2(1,32)*six)*xxx+rxyz(1)*xxxx
      eri(2,1,5,1)=r220+(+r7(8,4)*two-r6(4,6)*two+r5(8,12)*six+r5(3,22)*two-r4(4,11)*six &
&                  -r4(1,31)*two+r3(3,40)*six-r2(1,30)*six)*qx+rxyz(5)*xx
      eri(3,1,5,1)=r202+(+r7(10,4)*two-r6(6,6)*two+r5(10,12)*six+r5(3,22)*two &
&                  -r4(6,11)*six-r4(1,31)*two+r3(3,40)*six-r2(1,30)*six)*qx+(+r7(6,3)*two &
&                  -r6(3,5)*two+r5(6,11)*six+r5(6,21)*two-r4(3,10)*six-r4(3,30)*two &
&                  +r3(6,39)*six-r2(5,29)*six)*qz+(+r6(10,12)-r5(6,16)+r4(10,29)*three &
&                  +r4(3,41)-r3(6,22)*three-r3(1,46)+r2(5,55)*three-r1(1,34)*three)*xx &
&                  +rxyz(20)*xz+(+r6(3,10)-r5(1,14)+r4(3,27)*three+r4(3,39)-r3(1,20)*three &
&                  -r3(1,44)+r2(5,53)*three-r1(1,32)*three)*zz+(+r5(6,24)*two-r4(3,33)*two &
&                  +r3(6,42)*six-r2(5,32)*six)*xxz+(+r5(3,23)*two-r4(1,32)*two+r3(3,41)*six &
&                  -r2(1,31)*six)*xzz+rxyz(1)*xxzz
      eri(4,1,5,1)=r310+(+r7(5,3)+r7(5,4)*two-r6(2,5)-r6(2,6)*two+r5(5,11)*three &
&                  +r5(5,12)*six+r5(5,21)+r5(5,22)*two-r4(2,10)*three-r4(2,11)*six-r4(2,30) &
&                  -r4(2,31)*two+r3(5,39)*three+r3(5,40)*six-r2(4,29)*three-r2(4,30)*six)*qx &
&                  +(+r6(5,11)*two+r6(5,12)-r5(2,15)*two-r5(2,16)+r4(5,28)*six+r4(5,29)*three &
&                  -r3(2,21)*six-r3(2,22)*three)*xx+rxyz(2)*xxx
      eri(5,1,5,1)=r301+(+r7(6,3)+r7(6,4)*two-r6(3,5)-r6(3,6)*two+r5(6,11)*three &
&                  +r5(6,12)*six+r5(6,21)+r5(6,22)*two-r4(3,10)*three-r4(3,11)*six-r4(3,30) &
&                  -r4(3,31)*two+r3(6,39)*three+r3(6,40)*six-r2(5,29)*three-r2(5,30)*six)*qx &
&                  +(+r7(3,3)-r6(1,5)+r5(3,11)*three+r5(3,21)*three-r4(1,10)*three &
&                  -r4(1,30)*three+r3(3,39)*nine-r2(1,29)*nine)*qz+(+r6(6,11)*two+r6(6,12) &
&                  -r5(3,15)*two-r5(3,16)+r4(6,28)*six+r4(6,29)*three-r3(3,21)*six &
&                  -r3(3,22)*three)*xx+(+r6(3,10)+r6(3,11)*two-r5(1,14)-r5(1,15)*two &
&                  +r4(3,27)*three+r4(3,28)*six+r4(3,39)+r4(3,40)*two-r3(1,20)*three &
&                  -r3(1,21)*six-r3(1,44)-r3(1,45)*two+r2(5,53)*three+r2(5,54)*six &
&                  -r1(1,32)*three-r1(1,33)*six)*xz+(+r5(6,24)-r4(3,33)+r3(6,42)*three &
&                  -r2(5,32)*three)*xxx+(+r5(3,23)*two+r5(3,24)-r4(1,32)*two-r4(1,33) &
&                  +r3(3,41)*six+r3(3,42)*three-r2(1,31)*six-r2(1,32)*three)*xxz+rxyz(1)*xxxz
      eri(6,1,5,1)=r211+(+r7(9,4)*two-r6(5,6)*two+r5(9,12)*six-r4(5,11)*six)*qx+rxyz(17) &
&                  *qz+(+r6(9,12)-r5(5,16)+r4(9,29)*three-r3(5,22)*three)*xx+(+r6(5,11)*two &
&                  -r5(2,15)*two+r4(5,28)*six-r3(2,21)*six)*xz+rxyz(2)*xxz
      eri(1,2,5,1)=r220+(+r7(8,3)*two-r6(4,5)*two+r5(8,11)*six+r5(3,21)*two-r4(4,10)*six &
&                  -r4(1,30)*two+r3(3,39)*six-r2(1,29)*six)*qx+rxyz(4)*xx
      eri(2,2,5,1)=r040
      eri(3,2,5,1)=r022+(+r7(13,3)*two-r6(8,5)*two+r5(13,11)*six+r5(6,21)*two &
&                  -r4(8,10)*six-r4(3,30)*two+r3(6,39)*six-r2(5,29)*six)*qz+rxyz(4)*zz
      eri(4,2,5,1)=r130+rxyz(7)*qx
      eri(5,2,5,1)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,5,1)=r031+rxyz(7)*qz
      eri(1,3,5,1)=r202+(+r7(10,3)*two-r6(6,5)*two+r5(10,11)*six+r5(3,21)*two &
&                  -r4(6,10)*six-r4(1,30)*two+r3(3,39)*six-r2(1,29)*six)*qx+(+r7(6,4)*two &
&                  -r6(3,6)*two+r5(6,12)*six+r5(6,22)*two-r4(3,11)*six-r4(3,31)*two &
&                  +r3(6,40)*six-r2(5,30)*six)*qz+(+r6(10,10)-r5(6,14)+r4(10,27)*three &
&                  +r4(3,39)-r3(6,20)*three-r3(1,44)+r2(5,53)*three-r1(1,32)*three)*xx &
&                  +rxyz(20)*xz+(+r6(3,12)-r5(1,16)+r4(3,29)*three+r4(3,41)-r3(1,22)*three &
&                  -r3(1,46)+r2(5,55)*three-r1(1,34)*three)*zz+(+r5(6,23)*two-r4(3,32)*two &
&                  +r3(6,41)*six-r2(5,31)*six)*xxz+(+r5(3,24)*two-r4(1,33)*two+r3(3,42)*six &
&                  -r2(1,32)*six)*xzz+rxyz(1)*xxzz
      eri(2,3,5,1)=r022+(+r7(13,4)*two-r6(8,6)*two+r5(13,12)*six+r5(6,22)*two &
&                  -r4(8,11)*six-r4(3,31)*two+r3(6,40)*six-r2(5,30)*six)*qz+rxyz(5)*zz
      eri(3,3,5,1)=r004+(+r7(15,3)*two+r7(15,4)*two-r6(10,5)*two-r6(10,6)*two &
&                  +r5(15,11)*six+r5(15,12)*six+r5(6,21)*six+r5(6,22)*six-r4(10,10)*six &
&                  -r4(10,11)*six-r4(3,30)*six-r4(3,31)*six+r3(6,39)*p18+r3(6,40)*p18 &
&                  -r2(5,29)*p18-r2(5,30)*p18)*qz+(+r6(10,10)+r6(10,11)*four+r6(10,12) &
&                  -r5(6,14)-r5(6,15)*four-r5(6,16)+r4(10,27)*three+r4(10,28)*p12 &
&                  +r4(10,29)*three+r4(3,39)+r4(3,40)*four+r4(3,41)-r3(6,20)*three &
&                  -r3(6,21)*p12-r3(6,22)*three-r3(1,44)-r3(1,45)*four-r3(1,46) &
&                  +r2(5,53)*three+r2(5,54)*p12+r2(5,55)*three-r1(1,32)*three-r1(1,33)*p12 &
&                  -r1(1,34)*three)*zz+(+r5(6,23)*two+r5(6,24)*two-r4(3,32)*two-r4(3,33)*two &
&                  +r3(6,41)*six+r3(6,42)*six-r2(5,31)*six-r2(5,32)*six)*zzz+rxyz(1)*zzzz
      eri(4,3,5,1)=r112+rxyz(18)*qx+(+r7(9,4)*two-r6(5,6)*two+r5(9,12)*six-r4(5,11)*six) &
&                  *qz+(+r6(9,11)*two-r5(5,15)*two+r4(9,28)*six-r3(5,21)*six)*xz+(+r6(5,12) &
&                  -r5(2,16)+r4(5,29)*three-r3(2,22)*three)*zz+rxyz(2)*xzz
      eri(5,3,5,1)=r103+(+r7(15,3)-r6(10,5)+r5(15,11)*three+r5(6,21)*three &
&                  -r4(10,10)*three-r4(3,30)*three+r3(6,39)*nine-r2(5,29)*nine)*qx+(+r7(10,3) &
&                  +r7(10,4)*two-r6(6,5)-r6(6,6)*two+r5(10,11)*three+r5(10,12)*six+r5(3,21) &
&                  +r5(3,22)*two-r4(6,10)*three-r4(6,11)*six-r4(1,30)-r4(1,31)*two &
&                  +r3(3,39)*three+r3(3,40)*six-r2(1,29)*three-r2(1,30)*six)*qz+(+r6(10,10) &
&                  +r6(10,11)*two-r5(6,14)-r5(6,15)*two+r4(10,27)*three+r4(10,28)*six &
&                  +r4(3,39)+r4(3,40)*two-r3(6,20)*three-r3(6,21)*six-r3(1,44)-r3(1,45)*two &
&                  +r2(5,53)*three+r2(5,54)*six-r1(1,32)*three-r1(1,33)*six)*xz+( &
&                  +r6(6,11)*two+r6(6,12)-r5(3,15)*two-r5(3,16)+r4(6,28)*six+r4(6,29)*three &
&                  -r3(3,21)*six-r3(3,22)*three)*zz+(+r5(6,23)*two+r5(6,24)-r4(3,32)*two &
&                  -r4(3,33)+r3(6,41)*six+r3(6,42)*three-r2(5,31)*six-r2(5,32)*three)*xzz+( &
&                  +r5(3,24)-r4(1,33)+r3(3,42)*three-r2(1,32)*three)*zzz+rxyz(1)*xzzz
      eri(6,3,5,1)=r013+(+r7(14,3)+r7(14,4)*two-r6(9,5)-r6(9,6)*two+r5(14,11)*three &
&                  +r5(14,12)*six+r5(5,21)+r5(5,22)*two-r4(9,10)*three-r4(9,11)*six-r4(2,30) &
&                  -r4(2,31)*two+r3(5,39)*three+r3(5,40)*six-r2(4,29)*three-r2(4,30)*six)*qz &
&                  +(+r6(9,11)*two+r6(9,12)-r5(5,15)*two-r5(5,16)+r4(9,28)*six+r4(9,29)*three &
&                  -r3(5,21)*six-r3(5,22)*three)*zz+rxyz(2)*zzz
      eri(1,4,5,1)=r310+(+r7(5,3)*two+r7(5,4)-r6(2,5)*two-r6(2,6)+r5(5,11)*six &
&                  +r5(5,12)*three+r5(5,21)*two+r5(5,22)-r4(2,10)*six-r4(2,11)*three &
&                  -r4(2,30)*two-r4(2,31)+r3(5,39)*six+r3(5,40)*three-r2(4,29)*six &
&                  -r2(4,30)*three)*qx+(+r6(5,10)+r6(5,11)*two-r5(2,14)-r5(2,15)*two &
&                  +r4(5,27)*three+r4(5,28)*six-r3(2,20)*three-r3(2,21)*six)*xx+rxyz(3)*xxx
      eri(2,4,5,1)=r130+rxyz(8)*qx
      eri(3,4,5,1)=r112+rxyz(19)*qx+(+r7(9,3)*two-r6(5,5)*two+r5(9,11)*six-r4(5,10)*six) &
&                  *qz+(+r6(9,11)*two-r5(5,15)*two+r4(9,28)*six-r3(5,21)*six)*xz+(+r6(5,10) &
&                  -r5(2,14)+r4(5,27)*three-r3(2,20)*three)*zz+rxyz(3)*xzz
      eri(4,4,5,1)=r220+(+r7(8,3)+r7(8,4)-r6(4,5)-r6(4,6)+r5(8,11)*three+r5(8,12)*three &
&                  +r5(3,21)+r5(3,22)-r4(4,10)*three-r4(4,11)*three-r4(1,30)-r4(1,31) &
&                  +r3(3,39)*three+r3(3,40)*three-r2(1,29)*three-r2(1,30)*three)*qx+rxyz(6) &
&                  *xx
      eri(5,4,5,1)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(5,10)+r6(5,11)-r5(2,14) &
&                  -r5(2,15)+r4(5,27)*three+r4(5,28)*three-r3(2,20)*three-r3(2,21)*three)*xz &
&                  +rxyz(3)*xxz
      eri(6,4,5,1)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,5,1)=r301+(+r7(6,3)*two+r7(6,4)-r6(3,5)*two-r6(3,6)+r5(6,11)*six &
&                  +r5(6,12)*three+r5(6,21)*two+r5(6,22)-r4(3,10)*six-r4(3,11)*three &
&                  -r4(3,30)*two-r4(3,31)+r3(6,39)*six+r3(6,40)*three-r2(5,29)*six &
&                  -r2(5,30)*three)*qx+(+r7(3,4)-r6(1,6)+r5(3,12)*three+r5(3,22)*three &
&                  -r4(1,11)*three-r4(1,31)*three+r3(3,40)*nine-r2(1,30)*nine)*qz+(+r6(6,10) &
&                  +r6(6,11)*two-r5(3,14)-r5(3,15)*two+r4(6,27)*three+r4(6,28)*six &
&                  -r3(3,20)*three-r3(3,21)*six)*xx+(+r6(3,11)*two+r6(3,12)-r5(1,15)*two &
&                  -r5(1,16)+r4(3,28)*six+r4(3,29)*three+r4(3,40)*two+r4(3,41)-r3(1,21)*six &
&                  -r3(1,22)*three-r3(1,45)*two-r3(1,46)+r2(5,54)*six+r2(5,55)*three &
&                  -r1(1,33)*six-r1(1,34)*three)*xz+(+r5(6,23)-r4(3,32)+r3(6,41)*three &
&                  -r2(5,31)*three)*xxx+(+r5(3,23)+r5(3,24)*two-r4(1,32)-r4(1,33)*two &
&                  +r3(3,41)*three+r3(3,42)*six-r2(1,31)*three-r2(1,32)*six)*xxz+rxyz(1)*xxxz
      eri(2,5,5,1)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,5,1)=r103+(+r7(15,4)-r6(10,6)+r5(15,12)*three+r5(6,22)*three &
&                  -r4(10,11)*three-r4(3,31)*three+r3(6,40)*nine-r2(5,30)*nine)*qx+( &
&                  +r7(10,3)*two+r7(10,4)-r6(6,5)*two-r6(6,6)+r5(10,11)*six+r5(10,12)*three &
&                  +r5(3,21)*two+r5(3,22)-r4(6,10)*six-r4(6,11)*three-r4(1,30)*two-r4(1,31) &
&                  +r3(3,39)*six+r3(3,40)*three-r2(1,29)*six-r2(1,30)*three)*qz+( &
&                  +r6(10,11)*two+r6(10,12)-r5(6,15)*two-r5(6,16)+r4(10,28)*six &
&                  +r4(10,29)*three+r4(3,40)*two+r4(3,41)-r3(6,21)*six-r3(6,22)*three &
&                  -r3(1,45)*two-r3(1,46)+r2(5,54)*six+r2(5,55)*three-r1(1,33)*six &
&                  -r1(1,34)*three)*xz+(+r6(6,10)+r6(6,11)*two-r5(3,14)-r5(3,15)*two &
&                  +r4(6,27)*three+r4(6,28)*six-r3(3,20)*three-r3(3,21)*six)*zz+(+r5(6,23) &
&                  +r5(6,24)*two-r4(3,32)-r4(3,33)*two+r3(6,41)*three+r3(6,42)*six &
&                  -r2(5,31)*three-r2(5,32)*six)*xzz+(+r5(3,23)-r4(1,32)+r3(3,41)*three &
&                  -r2(1,31)*three)*zzz+rxyz(1)*xzzz
      eri(4,5,5,1)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(5,11)+r6(5,12)-r5(2,15) &
&                  -r5(2,16)+r4(5,28)*three+r4(5,29)*three-r3(2,21)*three-r3(2,22)*three)*xz &
&                  +rxyz(2)*xxz
      eri(5,5,5,1)=r202+(+r7(10,3)+r7(10,4)-r6(6,5)-r6(6,6)+r5(10,11)*three &
&                  +r5(10,12)*three+r5(3,21)+r5(3,22)-r4(6,10)*three-r4(6,11)*three-r4(1,30) &
&                  -r4(1,31)+r3(3,39)*three+r3(3,40)*three-r2(1,29)*three-r2(1,30)*three)*qx &
&                  +(+r7(6,3)+r7(6,4)-r6(3,5)-r6(3,6)+r5(6,11)*three+r5(6,12)*three+r5(6,21) &
&                  +r5(6,22)-r4(3,10)*three-r4(3,11)*three-r4(3,30)-r4(3,31)+r3(6,39)*three &
&                  +r3(6,40)*three-r2(5,29)*three-r2(5,30)*three)*qz+(+r6(10,11)-r5(6,15) &
&                  +r4(10,28)*three+r4(3,40)-r3(6,21)*three-r3(1,45)+r2(5,54)*three &
&                  -r1(1,33)*three)*xx+(+r6(6,10)+r6(6,11)*two+r6(6,12)-r5(3,14)-r5(3,15)*two &
&                  -r5(3,16)+r4(6,27)*three+r4(6,28)*six+r4(6,29)*three-r3(3,20)*three &
&                  -r3(3,21)*six-r3(3,22)*three)*xz+(+r6(3,11)-r5(1,15)+r4(3,28)*three &
&                  +r4(3,40)-r3(1,21)*three-r3(1,45)+r2(5,54)*three-r1(1,33)*three)*zz+( &
&                  +r5(6,23)+r5(6,24)-r4(3,32)-r4(3,33)+r3(6,41)*three+r3(6,42)*three &
&                  -r2(5,31)*three-r2(5,32)*three)*xxz+(+r5(3,23)+r5(3,24)-r4(1,32)-r4(1,33) &
&                  +r3(3,41)*three+r3(3,42)*three-r2(1,31)*three-r2(1,32)*three)*xzz+rxyz(1) &
&                  *xxzz
      eri(6,5,5,1)=r112+rxyz(19)*qx+(+r7(9,3)+r7(9,4)-r6(5,5)-r6(5,6)+r5(9,11)*three &
&                  +r5(9,12)*three-r4(5,10)*three-r4(5,11)*three)*qz+(+r6(9,11)+r6(9,12) &
&                  -r5(5,15)-r5(5,16)+r4(9,28)*three+r4(9,29)*three-r3(5,21)*three &
&                  -r3(5,22)*three)*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,5,1)=r211+(+r7(9,3)*two-r6(5,5)*two+r5(9,11)*six-r4(5,10)*six)*qx+rxyz(16) &
&                  *qz+(+r6(9,10)-r5(5,14)+r4(9,27)*three-r3(5,20)*three)*xx+(+r6(5,11)*two &
&                  -r5(2,15)*two+r4(5,28)*six-r3(2,21)*six)*xz+rxyz(3)*xxz
      eri(2,6,5,1)=r031+rxyz(8)*qz
      eri(3,6,5,1)=r013+(+r7(14,3)*two+r7(14,4)-r6(9,5)*two-r6(9,6)+r5(14,11)*six &
&                  +r5(14,12)*three+r5(5,21)*two+r5(5,22)-r4(9,10)*six-r4(9,11)*three &
&                  -r4(2,30)*two-r4(2,31)+r3(5,39)*six+r3(5,40)*three-r2(4,29)*six &
&                  -r2(4,30)*three)*qz+(+r6(9,10)+r6(9,11)*two-r5(5,14)-r5(5,15)*two &
&                  +r4(9,27)*three+r4(9,28)*six-r3(5,20)*three-r3(5,21)*six)*zz+rxyz(3)*zzz
      eri(4,6,5,1)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,5,1)=r112+rxyz(18)*qx+(+r7(9,3)+r7(9,4)-r6(5,5)-r6(5,6)+r5(9,11)*three &
&                  +r5(9,12)*three-r4(5,10)*three-r4(5,11)*three)*qz+(+r6(9,10)+r6(9,11) &
&                  -r5(5,14)-r5(5,15)+r4(9,27)*three+r4(9,28)*three-r3(5,20)*three &
&                  -r3(5,21)*three)*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,5,1)=r022+(+r7(13,3)+r7(13,4)-r6(8,5)-r6(8,6)+r5(13,11)*three &
&                  +r5(13,12)*three+r5(6,21)+r5(6,22)-r4(8,10)*three-r4(8,11)*three-r4(3,30) &
&                  -r4(3,31)+r3(6,39)*three+r3(6,40)*three-r2(5,29)*three-r2(5,30)*three)*qz &
&                  +rxyz(6)*zz
!
      r400= r8(5)-r7(2,1)+r6(5,4)+r6(5,9)*six-r5(2,3)-r5(2,13)*six+r4(5,26)*six &
&          +r4(5,38)*three-r3(2,19)*six-r3(2,43)*three+r2(6,52)*three-r1(2,31)*three
      r310= r8(8)-r7(4,1)+r6(8,4)+r6(8,9)*three-r5(4,3)-r5(4,13)*three+r4(8,26)*three &
&          -r3(4,19)*three
      r301= r8(9)-r7(5,1)+r6(9,4)+r6(9,9)*three-r5(5,3)-r5(5,13)*three+r4(9,26)*three &
&          -r3(5,19)*three
      r220= r8(12)-r7(7,1)+r6(12,4)+r6(5,9)+r6(12,9)-r5(7,3)-r5(2,13)-r5(7,13)+r4(5,26) &
&          +r4(12,26)+r4(5,38)-r3(2,19)-r3(7,19)-r3(2,43)+r2(6,52)-r1(2,31)
      r211= r8(13)-r7(8,1)+r6(13,4)+r6(13,9)-r5(8,3)-r5(8,13)+r4(13,26)-r3(8,19)
      r202= r8(14)-r7(9,1)+r6(14,4)+r6(5,9)+r6(14,9)-r5(9,3)-r5(2,13)-r5(9,13)+r4(5,26) &
&          +r4(14,26)+r4(5,38)-r3(2,19)-r3(9,19)-r3(2,43)+r2(6,52)-r1(2,31)
      r130= r8(17)-r7(11,1)+r6(17,4)+r6(8,9)*three-r5(11,3)-r5(4,13)*three+r4(8,26)*three &
&          -r3(4,19)*three
      r121= r8(18)-r7(12,1)+r6(18,4)+r6(9,9)-r5(12,3)-r5(5,13)+r4(9,26)-r3(5,19)
      r112= r8(19)-r7(13,1)+r6(19,4)+r6(8,9)-r5(13,3)-r5(4,13)+r4(8,26)-r3(4,19)
      r103= r8(20)-r7(14,1)+r6(20,4)+r6(9,9)*three-r5(14,3)-r5(5,13)*three+r4(9,26)*three &
&          -r3(5,19)*three
      r040= r8(23)-r7(16,1)+r6(23,4)+r6(12,9)*six-r5(16,3)-r5(7,13)*six+r4(12,26)*six &
&          +r4(5,38)*three-r3(7,19)*six-r3(2,43)*three+r2(6,52)*three-r1(2,31)*three
      r031= r8(24)-r7(17,1)+r6(24,4)+r6(13,9)*three-r5(17,3)-r5(8,13)*three+r4(13,26)*three &
&          -r3(8,19)*three
      r022= r8(25)-r7(18,1)+r6(25,4)+r6(12,9)+r6(14,9)-r5(18,3)-r5(7,13)-r5(9,13)+r4(12,26) &
&          +r4(14,26)+r4(5,38)-r3(7,19)-r3(9,19)-r3(2,43)+r2(6,52)-r1(2,31)
      r013= r8(26)-r7(19,1)+r6(26,4)+r6(13,9)*three-r5(19,3)-r5(8,13)*three+r4(13,26)*three &
&          -r3(8,19)*three
      r004= r8(27)-r7(20,1)+r6(27,4)+r6(14,9)*six-r5(20,3)-r5(9,13)*six+r4(14,26)*six &
&          +r4(5,38)*three-r3(9,19)*six-r3(2,43)*three+r2(6,52)*three-r1(2,31)*three
      rxyz(1)=+r4(5,42)-r3(2,47)+r2(6,56)-r1(2,35)
      rxyz(2)=+r5(8,24)-r4(4,33)+r3(8,42)-r2(2,32)
      rxyz(3)=+r5(8,23)-r4(4,32)+r3(8,41)-r2(2,31)
      rxyz(4)=+r6(12,10)-r5(7,14)+r4(12,27)+r4(5,39)-r3(7,20)-r3(2,44)+r2(6,53)-r1(2,32)
      rxyz(5)=+r6(12,12)-r5(7,16)+r4(12,29)+r4(5,41)-r3(7,22)-r3(2,46)+r2(6,55)-r1(2,34)
      rxyz(6)=+r6(12,11)-r5(7,15)+r4(12,28)+r4(5,40)-r3(7,21)-r3(2,45)+r2(6,54)-r1(2,33)
      rxyz(7)=+r7(17,3)-r6(11,5)+r5(17,11)+r5(8,21)*three-r4(11,10)-r4(4,30)*three &
&             +r3(8,39)*three-r2(2,29)*three
      rxyz(8)=+r7(17,4)-r6(11,6)+r5(17,12)+r5(8,22)*three-r4(11,11)-r4(4,31)*three &
&             +r3(8,40)*three-r2(2,30)*three
      rxyz(9)=+r7(13,3)+r7(13,4)-r6(8,5)-r6(8,6)+r5(13,11)+r5(13,12)-r4(8,10)-r4(8,11)
      rxyz(10)=+r6(13,11)-r5(8,15)+r4(13,28)-r3(8,21)
      rxyz(11)=+r6(8,11)-r5(4,15)+r4(8,28)-r3(4,21)
      rxyz(12)=+r7(18,3)-r6(12,5)+r5(18,11)+r5(9,21)-r4(12,10)-r4(5,30)+r3(9,39)-r2(6,29)
      rxyz(13)=+r7(18,4)-r6(12,6)+r5(18,12)+r5(9,22)-r4(12,11)-r4(5,31)+r3(9,40)-r2(6,30)
      rxyz(14)=+r7(12,4)-r6(7,6)+r5(12,12)+r5(5,22)-r4(7,11)-r4(2,31)+r3(5,40)-r2(4,30)
      rxyz(15)=+r7(12,3)-r6(7,5)+r5(12,11)+r5(5,21)-r4(7,10)-r4(2,30)+r3(5,39)-r2(4,29)
      rxyz(16)=+r7(8,4)-r6(4,6)+r5(8,12)+r5(8,22)-r4(4,11)-r4(4,31)+r3(8,40)-r2(2,30)
      rxyz(17)=+r7(8,3)-r6(4,5)+r5(8,11)+r5(8,21)-r4(4,10)-r4(4,30)+r3(8,39)-r2(2,29)
      rxyz(18)=+r7(19,3)-r6(13,5)+r5(19,11)+r5(8,21)-r4(13,10)-r4(4,30)+r3(8,39)-r2(2,29)
      rxyz(19)=+r7(19,4)-r6(13,6)+r5(19,12)+r5(8,22)-r4(13,11)-r4(4,31)+r3(8,40)-r2(2,30)
      rxyz(20)=+r6(9,11)*four-r5(5,15)*four+r4(9,28)*four-r3(5,21)*four
      eri(1,1,6,1)=r400+(+r7(5,3)*two+r7(5,4)*two-r6(2,5)*two-r6(2,6)*two+r5(5,11)*two &
&                  +r5(5,12)*two+r5(5,21)*six+r5(5,22)*six-r4(2,10)*two-r4(2,11)*two &
&                  -r4(2,30)*six-r4(2,31)*six+r3(5,39)*six+r3(5,40)*six-r2(4,29)*six &
&                  -r2(4,30)*six)*qx+(+r6(5,10)+r6(5,11)*four+r6(5,12)-r5(2,14)-r5(2,15)*four &
&                  -r5(2,16)+r4(5,27)+r4(5,28)*four+r4(5,29)+r4(5,39)+r4(5,40)*four+r4(5,41) &
&                  -r3(2,20)-r3(2,21)*four-r3(2,22)-r3(2,44)-r3(2,45)*four-r3(2,46)+r2(6,53) &
&                  +r2(6,54)*four+r2(6,55)-r1(2,32)-r1(2,33)*four-r1(2,34))*xx+(+r5(5,23)*two &
&                  +r5(5,24)*two-r4(2,32)*two-r4(2,33)*two+r3(5,41)*two+r3(5,42)*two &
&                  -r2(4,31)*two-r2(4,32)*two)*xxx+rxyz(1)*xxxx
      eri(2,1,6,1)=r220+(+r7(12,4)*two-r6(7,6)*two+r5(12,12)*two+r5(5,22)*two &
&                  -r4(7,11)*two-r4(2,31)*two+r3(5,40)*two-r2(4,30)*two)*qx+rxyz(5)*xx
      eri(3,1,6,1)=r202+(+r7(14,4)*two-r6(9,6)*two+r5(14,12)*two+r5(5,22)*two &
&                  -r4(9,11)*two-r4(2,31)*two+r3(5,40)*two-r2(4,30)*two)*qx+(+r7(9,3)*two &
&                  -r6(5,5)*two+r5(9,11)*two+r5(9,21)*two-r4(5,10)*two-r4(5,30)*two &
&                  +r3(9,39)*two-r2(6,29)*two)*qz+(+r6(14,12)-r5(9,16)+r4(14,29)+r4(5,41) &
&                  -r3(9,22)-r3(2,46)+r2(6,55)-r1(2,34))*xx+rxyz(20)*xz+(+r6(5,10)-r5(2,14) &
&                  +r4(5,27)+r4(5,39)-r3(2,20)-r3(2,44)+r2(6,53)-r1(2,32))*zz+(+r5(9,24)*two &
&                  -r4(5,33)*two+r3(9,42)*two-r2(6,32)*two)*xxz+(+r5(5,23)*two-r4(2,32)*two &
&                  +r3(5,41)*two-r2(4,31)*two)*xzz+rxyz(1)*xxzz
      eri(4,1,6,1)=r310+(+r7(8,3)+r7(8,4)*two-r6(4,5)-r6(4,6)*two+r5(8,11)+r5(8,12)*two &
&                  +r5(8,21)+r5(8,22)*two-r4(4,10)-r4(4,11)*two-r4(4,30)-r4(4,31)*two &
&                  +r3(8,39)+r3(8,40)*two-r2(2,29)-r2(2,30)*two)*qx+(+r6(8,11)*two+r6(8,12) &
&                  -r5(4,15)*two-r5(4,16)+r4(8,28)*two+r4(8,29)-r3(4,21)*two-r3(4,22))*xx &
&                  +rxyz(2)*xxx
      eri(5,1,6,1)=r301+(+r7(9,3)+r7(9,4)*two-r6(5,5)-r6(5,6)*two+r5(9,11)+r5(9,12)*two &
&                  +r5(9,21)+r5(9,22)*two-r4(5,10)-r4(5,11)*two-r4(5,30)-r4(5,31)*two &
&                  +r3(9,39)+r3(9,40)*two-r2(6,29)-r2(6,30)*two)*qx+(+r7(5,3)-r6(2,5) &
&                  +r5(5,11)+r5(5,21)*three-r4(2,10)-r4(2,30)*three+r3(5,39)*three &
&                  -r2(4,29)*three)*qz+(+r6(9,11)*two+r6(9,12)-r5(5,15)*two-r5(5,16) &
&                  +r4(9,28)*two+r4(9,29)-r3(5,21)*two-r3(5,22))*xx+(+r6(5,10)+r6(5,11)*two &
&                  -r5(2,14)-r5(2,15)*two+r4(5,27)+r4(5,28)*two+r4(5,39)+r4(5,40)*two &
&                  -r3(2,20)-r3(2,21)*two-r3(2,44)-r3(2,45)*two+r2(6,53)+r2(6,54)*two &
&                  -r1(2,32)-r1(2,33)*two)*xz+(+r5(9,24)-r4(5,33)+r3(9,42)-r2(6,32))*xxx+( &
&                  +r5(5,23)*two+r5(5,24)-r4(2,32)*two-r4(2,33)+r3(5,41)*two+r3(5,42) &
&                  -r2(4,31)*two-r2(4,32))*xxz+rxyz(1)*xxxz
      eri(6,1,6,1)=r211+(+r7(13,4)*two-r6(8,6)*two+r5(13,12)*two-r4(8,11)*two)*qx &
&                  +rxyz(17)*qz+(+r6(13,12)-r5(8,16)+r4(13,29)-r3(8,22))*xx+(+r6(8,11)*two &
&                  -r5(4,15)*two+r4(8,28)*two-r3(4,21)*two)*xz+rxyz(2)*xxz
      eri(1,2,6,1)=r220+(+r7(12,3)*two-r6(7,5)*two+r5(12,11)*two+r5(5,21)*two &
&                  -r4(7,10)*two-r4(2,30)*two+r3(5,39)*two-r2(4,29)*two)*qx+rxyz(4)*xx
      eri(2,2,6,1)=r040
      eri(3,2,6,1)=r022+(+r7(18,3)*two-r6(12,5)*two+r5(18,11)*two+r5(9,21)*two &
&                  -r4(12,10)*two-r4(5,30)*two+r3(9,39)*two-r2(6,29)*two)*qz+rxyz(4)*zz
      eri(4,2,6,1)=r130+rxyz(7)*qx
      eri(5,2,6,1)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,6,1)=r031+rxyz(7)*qz
      eri(1,3,6,1)=r202+(+r7(14,3)*two-r6(9,5)*two+r5(14,11)*two+r5(5,21)*two &
&                  -r4(9,10)*two-r4(2,30)*two+r3(5,39)*two-r2(4,29)*two)*qx+(+r7(9,4)*two &
&                  -r6(5,6)*two+r5(9,12)*two+r5(9,22)*two-r4(5,11)*two-r4(5,31)*two &
&                  +r3(9,40)*two-r2(6,30)*two)*qz+(+r6(14,10)-r5(9,14)+r4(14,27)+r4(5,39) &
&                  -r3(9,20)-r3(2,44)+r2(6,53)-r1(2,32))*xx+rxyz(20)*xz+(+r6(5,12)-r5(2,16) &
&                  +r4(5,29)+r4(5,41)-r3(2,22)-r3(2,46)+r2(6,55)-r1(2,34))*zz+(+r5(9,23)*two &
&                  -r4(5,32)*two+r3(9,41)*two-r2(6,31)*two)*xxz+(+r5(5,24)*two-r4(2,33)*two &
&                  +r3(5,42)*two-r2(4,32)*two)*xzz+rxyz(1)*xxzz
      eri(2,3,6,1)=r022+(+r7(18,4)*two-r6(12,6)*two+r5(18,12)*two+r5(9,22)*two &
&                  -r4(12,11)*two-r4(5,31)*two+r3(9,40)*two-r2(6,30)*two)*qz+rxyz(5)*zz
      eri(3,3,6,1)=r004+(+r7(20,3)*two+r7(20,4)*two-r6(14,5)*two-r6(14,6)*two &
&                  +r5(20,11)*two+r5(20,12)*two+r5(9,21)*six+r5(9,22)*six-r4(14,10)*two &
&                  -r4(14,11)*two-r4(5,30)*six-r4(5,31)*six+r3(9,39)*six+r3(9,40)*six &
&                  -r2(6,29)*six-r2(6,30)*six)*qz+(+r6(14,10)+r6(14,11)*four+r6(14,12) &
&                  -r5(9,14)-r5(9,15)*four-r5(9,16)+r4(14,27)+r4(14,28)*four+r4(14,29) &
&                  +r4(5,39)+r4(5,40)*four+r4(5,41)-r3(9,20)-r3(9,21)*four-r3(9,22)-r3(2,44) &
&                  -r3(2,45)*four-r3(2,46)+r2(6,53)+r2(6,54)*four+r2(6,55)-r1(2,32) &
&                  -r1(2,33)*four-r1(2,34))*zz+(+r5(9,23)*two+r5(9,24)*two-r4(5,32)*two &
&                  -r4(5,33)*two+r3(9,41)*two+r3(9,42)*two-r2(6,31)*two-r2(6,32)*two)*zzz &
&                  +rxyz(1)*zzzz
      eri(4,3,6,1)=r112+rxyz(18)*qx+(+r7(13,4)*two-r6(8,6)*two+r5(13,12)*two-r4(8,11)*two &
&                  )*qz+(+r6(13,11)*two-r5(8,15)*two+r4(13,28)*two-r3(8,21)*two)*xz+( &
&                  +r6(8,12)-r5(4,16)+r4(8,29)-r3(4,22))*zz+rxyz(2)*xzz
      eri(5,3,6,1)=r103+(+r7(20,3)-r6(14,5)+r5(20,11)+r5(9,21)*three-r4(14,10) &
&                  -r4(5,30)*three+r3(9,39)*three-r2(6,29)*three)*qx+(+r7(14,3)+r7(14,4)*two &
&                  -r6(9,5)-r6(9,6)*two+r5(14,11)+r5(14,12)*two+r5(5,21)+r5(5,22)*two &
&                  -r4(9,10)-r4(9,11)*two-r4(2,30)-r4(2,31)*two+r3(5,39)+r3(5,40)*two &
&                  -r2(4,29)-r2(4,30)*two)*qz+(+r6(14,10)+r6(14,11)*two-r5(9,14)-r5(9,15)*two &
&                  +r4(14,27)+r4(14,28)*two+r4(5,39)+r4(5,40)*two-r3(9,20)-r3(9,21)*two &
&                  -r3(2,44)-r3(2,45)*two+r2(6,53)+r2(6,54)*two-r1(2,32)-r1(2,33)*two)*xz+( &
&                  +r6(9,11)*two+r6(9,12)-r5(5,15)*two-r5(5,16)+r4(9,28)*two+r4(9,29) &
&                  -r3(5,21)*two-r3(5,22))*zz+(+r5(9,23)*two+r5(9,24)-r4(5,32)*two-r4(5,33) &
&                  +r3(9,41)*two+r3(9,42)-r2(6,31)*two-r2(6,32))*xzz+(+r5(5,24)-r4(2,33) &
&                  +r3(5,42)-r2(4,32))*zzz+rxyz(1)*xzzz
      eri(6,3,6,1)=r013+(+r7(19,3)+r7(19,4)*two-r6(13,5)-r6(13,6)*two+r5(19,11) &
&                  +r5(19,12)*two+r5(8,21)+r5(8,22)*two-r4(13,10)-r4(13,11)*two-r4(4,30) &
&                  -r4(4,31)*two+r3(8,39)+r3(8,40)*two-r2(2,29)-r2(2,30)*two)*qz+( &
&                  +r6(13,11)*two+r6(13,12)-r5(8,15)*two-r5(8,16)+r4(13,28)*two+r4(13,29) &
&                  -r3(8,21)*two-r3(8,22))*zz+rxyz(2)*zzz
      eri(1,4,6,1)=r310+(+r7(8,3)*two+r7(8,4)-r6(4,5)*two-r6(4,6)+r5(8,11)*two+r5(8,12) &
&                  +r5(8,21)*two+r5(8,22)-r4(4,10)*two-r4(4,11)-r4(4,30)*two-r4(4,31) &
&                  +r3(8,39)*two+r3(8,40)-r2(2,29)*two-r2(2,30))*qx+(+r6(8,10)+r6(8,11)*two &
&                  -r5(4,14)-r5(4,15)*two+r4(8,27)+r4(8,28)*two-r3(4,20)-r3(4,21)*two)*xx &
&                  +rxyz(3)*xxx
      eri(2,4,6,1)=r130+rxyz(8)*qx
      eri(3,4,6,1)=r112+rxyz(19)*qx+(+r7(13,3)*two-r6(8,5)*two+r5(13,11)*two-r4(8,10)*two &
&                  )*qz+(+r6(13,11)*two-r5(8,15)*two+r4(13,28)*two-r3(8,21)*two)*xz+( &
&                  +r6(8,10)-r5(4,14)+r4(8,27)-r3(4,20))*zz+rxyz(3)*xzz
      eri(4,4,6,1)=r220+(+r7(12,3)+r7(12,4)-r6(7,5)-r6(7,6)+r5(12,11)+r5(12,12)+r5(5,21) &
&                  +r5(5,22)-r4(7,10)-r4(7,11)-r4(2,30)-r4(2,31)+r3(5,39)+r3(5,40)-r2(4,29) &
&                  -r2(4,30))*qx+rxyz(6)*xx
      eri(5,4,6,1)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(8,10)+r6(8,11)-r5(4,14) &
&                  -r5(4,15)+r4(8,27)+r4(8,28)-r3(4,20)-r3(4,21))*xz+rxyz(3)*xxz
      eri(6,4,6,1)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,6,1)=r301+(+r7(9,3)*two+r7(9,4)-r6(5,5)*two-r6(5,6)+r5(9,11)*two+r5(9,12) &
&                  +r5(9,21)*two+r5(9,22)-r4(5,10)*two-r4(5,11)-r4(5,30)*two-r4(5,31) &
&                  +r3(9,39)*two+r3(9,40)-r2(6,29)*two-r2(6,30))*qx+(+r7(5,4)-r6(2,6) &
&                  +r5(5,12)+r5(5,22)*three-r4(2,11)-r4(2,31)*three+r3(5,40)*three &
&                  -r2(4,30)*three)*qz+(+r6(9,10)+r6(9,11)*two-r5(5,14)-r5(5,15)*two+r4(9,27) &
&                  +r4(9,28)*two-r3(5,20)-r3(5,21)*two)*xx+(+r6(5,11)*two+r6(5,12) &
&                  -r5(2,15)*two-r5(2,16)+r4(5,28)*two+r4(5,29)+r4(5,40)*two+r4(5,41) &
&                  -r3(2,21)*two-r3(2,22)-r3(2,45)*two-r3(2,46)+r2(6,54)*two+r2(6,55) &
&                  -r1(2,33)*two-r1(2,34))*xz+(+r5(9,23)-r4(5,32)+r3(9,41)-r2(6,31))*xxx+( &
&                  +r5(5,23)+r5(5,24)*two-r4(2,32)-r4(2,33)*two+r3(5,41)+r3(5,42)*two &
&                  -r2(4,31)-r2(4,32)*two)*xxz+rxyz(1)*xxxz
      eri(2,5,6,1)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,6,1)=r103+(+r7(20,4)-r6(14,6)+r5(20,12)+r5(9,22)*three-r4(14,11) &
&                  -r4(5,31)*three+r3(9,40)*three-r2(6,30)*three)*qx+(+r7(14,3)*two+r7(14,4) &
&                  -r6(9,5)*two-r6(9,6)+r5(14,11)*two+r5(14,12)+r5(5,21)*two+r5(5,22) &
&                  -r4(9,10)*two-r4(9,11)-r4(2,30)*two-r4(2,31)+r3(5,39)*two+r3(5,40) &
&                  -r2(4,29)*two-r2(4,30))*qz+(+r6(14,11)*two+r6(14,12)-r5(9,15)*two-r5(9,16) &
&                  +r4(14,28)*two+r4(14,29)+r4(5,40)*two+r4(5,41)-r3(9,21)*two-r3(9,22) &
&                  -r3(2,45)*two-r3(2,46)+r2(6,54)*two+r2(6,55)-r1(2,33)*two-r1(2,34))*xz+( &
&                  +r6(9,10)+r6(9,11)*two-r5(5,14)-r5(5,15)*two+r4(9,27)+r4(9,28)*two &
&                  -r3(5,20)-r3(5,21)*two)*zz+(+r5(9,23)+r5(9,24)*two-r4(5,32)-r4(5,33)*two &
&                  +r3(9,41)+r3(9,42)*two-r2(6,31)-r2(6,32)*two)*xzz+(+r5(5,23)-r4(2,32) &
&                  +r3(5,41)-r2(4,31))*zzz+rxyz(1)*xzzz
      eri(4,5,6,1)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(8,11)+r6(8,12)-r5(4,15) &
&                  -r5(4,16)+r4(8,28)+r4(8,29)-r3(4,21)-r3(4,22))*xz+rxyz(2)*xxz
      eri(5,5,6,1)=r202+(+r7(14,3)+r7(14,4)-r6(9,5)-r6(9,6)+r5(14,11)+r5(14,12)+r5(5,21) &
&                  +r5(5,22)-r4(9,10)-r4(9,11)-r4(2,30)-r4(2,31)+r3(5,39)+r3(5,40)-r2(4,29) &
&                  -r2(4,30))*qx+(+r7(9,3)+r7(9,4)-r6(5,5)-r6(5,6)+r5(9,11)+r5(9,12)+r5(9,21) &
&                  +r5(9,22)-r4(5,10)-r4(5,11)-r4(5,30)-r4(5,31)+r3(9,39)+r3(9,40)-r2(6,29) &
&                  -r2(6,30))*qz+(+r6(14,11)-r5(9,15)+r4(14,28)+r4(5,40)-r3(9,21)-r3(2,45) &
&                  +r2(6,54)-r1(2,33))*xx+(+r6(9,10)+r6(9,11)*two+r6(9,12)-r5(5,14) &
&                  -r5(5,15)*two-r5(5,16)+r4(9,27)+r4(9,28)*two+r4(9,29)-r3(5,20) &
&                  -r3(5,21)*two-r3(5,22))*xz+(+r6(5,11)-r5(2,15)+r4(5,28)+r4(5,40)-r3(2,21) &
&                  -r3(2,45)+r2(6,54)-r1(2,33))*zz+(+r5(9,23)+r5(9,24)-r4(5,32)-r4(5,33) &
&                  +r3(9,41)+r3(9,42)-r2(6,31)-r2(6,32))*xxz+(+r5(5,23)+r5(5,24)-r4(2,32) &
&                  -r4(2,33)+r3(5,41)+r3(5,42)-r2(4,31)-r2(4,32))*xzz+rxyz(1)*xxzz
      eri(6,5,6,1)=r112+rxyz(19)*qx+(+r7(13,3)+r7(13,4)-r6(8,5)-r6(8,6)+r5(13,11) &
&                  +r5(13,12)-r4(8,10)-r4(8,11))*qz+(+r6(13,11)+r6(13,12)-r5(8,15)-r5(8,16) &
&                  +r4(13,28)+r4(13,29)-r3(8,21)-r3(8,22))*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,6,1)=r211+(+r7(13,3)*two-r6(8,5)*two+r5(13,11)*two-r4(8,10)*two)*qx &
&                  +rxyz(16)*qz+(+r6(13,10)-r5(8,14)+r4(13,27)-r3(8,20))*xx+(+r6(8,11)*two &
&                  -r5(4,15)*two+r4(8,28)*two-r3(4,21)*two)*xz+rxyz(3)*xxz
      eri(2,6,6,1)=r031+rxyz(8)*qz
      eri(3,6,6,1)=r013+(+r7(19,3)*two+r7(19,4)-r6(13,5)*two-r6(13,6)+r5(19,11)*two &
&                  +r5(19,12)+r5(8,21)*two+r5(8,22)-r4(13,10)*two-r4(13,11)-r4(4,30)*two &
&                  -r4(4,31)+r3(8,39)*two+r3(8,40)-r2(2,29)*two-r2(2,30))*qz+(+r6(13,10) &
&                  +r6(13,11)*two-r5(8,14)-r5(8,15)*two+r4(13,27)+r4(13,28)*two-r3(8,20) &
&                  -r3(8,21)*two)*zz+rxyz(3)*zzz
      eri(4,6,6,1)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,6,1)=r112+rxyz(18)*qx+(+r7(13,3)+r7(13,4)-r6(8,5)-r6(8,6)+r5(13,11) &
&                  +r5(13,12)-r4(8,10)-r4(8,11))*qz+(+r6(13,10)+r6(13,11)-r5(8,14)-r5(8,15) &
&                  +r4(13,27)+r4(13,28)-r3(8,20)-r3(8,21))*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,6,1)=r022+(+r7(18,3)+r7(18,4)-r6(12,5)-r6(12,6)+r5(18,11)+r5(18,12) &
&                  +r5(9,21)+r5(9,22)-r4(12,10)-r4(12,11)-r4(5,30)-r4(5,31)+r3(9,39)+r3(9,40) &
&                  -r2(6,29)-r2(6,30))*qz+rxyz(6)*zz
      return
end

!-------------------------------------------------------------
  subroutine int2dddd2(eri,r0,r1,r2,r3,r4,r5,r6,r7,r8,qx,qz)
!-------------------------------------------------------------
!
      implicit none
      integer :: i,j
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00, four=4.0D+00, five=5.0D+00
      real(8),parameter :: six=6.0D+00, seven=7.0D+00, eight=8.0D+00, nine=9.0D+00, ten=1.0D+01
      real(8),parameter :: p11=1.1D+01, p12=1.2D+01, p13=1.3D+01, p15=1.5D+01, p16=1.6D+01
      real(8),parameter :: p18=1.8D+01, p20=2.0D+01, p21=2.1D+01, p24=2.4D+01, p28=2.8D+01
      real(8),parameter :: p30=3.0D+01, p36=3.6D+01, p45=4.5D+01, p105=1.05D+2, p210=2.1D+02
      real(8),parameter :: p420=4.2D+02
      real(8),intent(in) :: r0(25), r1(3,40), r2(6,56), r3(10,52), r4(15,42), r5(21,24)
      real(8),intent(in) :: r6(28,12), r7(36,4), r8(45), qx, qz
      real(8),intent(inout) :: eri(6,6,6,6)
      real(8) :: r400, r310, r301, r220, r211, r202, r130, r121, r112, r103, r040, r031
      real(8) :: r022, r013, r004, rxyz(20), xx, xz, zz, xxx, xxz, xzz, zzz
      real(8) :: xxxx, xxxz, xxzz, xzzz, zzzz
!
      xx= qx*qx
      xz= qx*qz
      zz= qz*qz
      xxx= qx*qx*qx
      xxz= qx*qx*qz
      xzz= qx*qz*qz
      zzz= qz*qz*qz
      xxxx= xx*xx
      xxxz= xx*xz
      xxzz= xx*zz
      xzzz= xz*zz
      zzzz= zz*zz
!
!
      do i= 1,6
        do j= 1,6
          eri(j,i,1,2)= eri(j,i,2,1)
        enddo
      enddo
!
      r400= r8(11)+r6(4,4)*six+r6(11,9)*six+r4(1,5)*three+r4(4,26)*p36+r4(11,38)*three &
&          +r2(1,17)*p18+r2(2,52)*p18+r0(21)*nine
      r310= r8(16)+r6(7,4)*six+r6(16,9)*three+r4(2,5)*three+r4(7,26)*p18+r2(4,17)*nine
      r301= r8(17)+r6(8,4)*six+r6(17,9)*three+r4(3,5)*three+r4(8,26)*p18+r2(5,17)*nine
      r220= r8(22)+r6(11,4)*six+r6(11,9)+r6(22,9)+r4(4,5)*three+r4(4,26)*six+r4(11,26)*six &
&          +r4(11,38)+r2(1,17)*three+r2(2,17)*three+r2(2,52)*six+r0(21)*three
      r211= r8(23)+r6(12,4)*six+r6(23,9)+r4(5,5)*three+r4(12,26)*six+r2(6,17)*three
      r202= r8(24)+r6(13,4)*six+r6(11,9)+r6(24,9)+r4(6,5)*three+r4(4,26)*six+r4(13,26)*six &
&          +r4(11,38)+r2(1,17)*three+r2(3,17)*three+r2(2,52)*six+r0(21)*three
      r130= r8(29)+r6(16,4)*six+r6(16,9)*three+r4(7,5)*three+r4(7,26)*p18+r2(4,17)*nine
      r121= r8(30)+r6(17,4)*six+r6(17,9)+r4(8,5)*three+r4(8,26)*six+r2(5,17)*three
      r112= r8(31)+r6(18,4)*six+r6(16,9)+r4(9,5)*three+r4(7,26)*six+r2(4,17)*three
      r103= r8(32)+r6(19,4)*six+r6(17,9)*three+r4(10,5)*three+r4(8,26)*p18+r2(5,17)*nine
      r040= r8(37)+r6(22,4)*six+r6(22,9)*six+r4(11,5)*three+r4(11,26)*p36+r4(11,38)*three &
&          +r2(2,17)*p18+r2(2,52)*p18+r0(21)*nine
      r031= r8(38)+r6(23,4)*six+r6(23,9)*three+r4(12,5)*three+r4(12,26)*p18+r2(6,17)*nine
      r022= r8(39)+r6(24,4)*six+r6(22,9)+r6(24,9)+r4(13,5)*three+r4(11,26)*six &
&          +r4(13,26)*six+r4(11,38)+r2(2,17)*three+r2(3,17)*three+r2(2,52)*six &
&          +r0(21)*three
      r013= r8(40)+r6(25,4)*six+r6(23,9)*three+r4(14,5)*three+r4(12,26)*p18+r2(6,17)*nine
      r004= r8(41)+r6(26,4)*six+r6(24,9)*six+r4(15,5)*three+r4(13,26)*p36+r4(11,38)*three &
&          +r2(3,17)*p18+r2(2,52)*p18+r0(21)*nine
      rxyz(1)=+r4(11,42)+r2(2,56)*six+r0(25)*three
      rxyz(2)=+r5(16,24)+r3(7,42)*six+r1(2,20)*three
      rxyz(3)=+r5(16,23)+r3(7,41)*six+r1(2,19)*three
      rxyz(4)=+r6(22,10)+r4(11,27)*six+r4(11,39)+r2(2,18)*three+r2(2,53)*six+r0(22)*three
      rxyz(5)=+r6(22,12)+r4(11,29)*six+r4(11,41)+r2(2,20)*three+r2(2,55)*six+r0(24)*three
      rxyz(6)=+r6(22,11)+r4(11,28)*six+r4(11,40)+r2(2,19)*three+r2(2,54)*six+r0(23)*three
      rxyz(7)=+r7(29,3)+r5(16,11)*six+r5(16,21)*three+r3(7,9)*three+r3(7,39)*p18 &
&             +r1(2,17)*nine
      rxyz(8)=+r7(29,4)+r5(16,12)*six+r5(16,22)*three+r3(7,10)*three+r3(7,40)*p18 &
&             +r1(2,18)*nine
      rxyz(9)=+r7(23,3)+r7(23,4)+r5(12,11)*six+r5(12,12)*six+r3(5,9)*three+r3(5,10)*three
      rxyz(10)=+r6(23,11)+r4(12,28)*six+r2(6,19)*three
      rxyz(11)=+r6(16,11)+r4(7,28)*six+r2(4,19)*three
      rxyz(12)=+r7(30,3)+r5(17,11)*six+r5(17,21)+r3(8,9)*three+r3(8,39)*six+r1(3,17)*three
      rxyz(13)=+r7(30,4)+r5(17,12)*six+r5(17,22)+r3(8,10)*three+r3(8,40)*six &
&             +r1(3,18)*three
      rxyz(14)=+r7(22,4)+r5(11,12)*six+r5(11,22)+r3(4,10)*three+r3(4,40)*six &
&             +r1(1,18)*three
      rxyz(15)=+r7(22,3)+r5(11,11)*six+r5(11,21)+r3(4,9)*three+r3(4,39)*six+r1(1,17)*three
      rxyz(16)=+r7(16,4)+r5(7,12)*six+r5(16,22)+r3(2,10)*three+r3(7,40)*six+r1(2,18)*three
      rxyz(17)=+r7(16,3)+r5(7,11)*six+r5(16,21)+r3(2,9)*three+r3(7,39)*six+r1(2,17)*three
      rxyz(18)=+r7(31,3)+r5(18,11)*six+r5(16,21)+r3(9,9)*three+r3(7,39)*six+r1(2,17)*three
      rxyz(19)=+r7(31,4)+r5(18,12)*six+r5(16,22)+r3(9,10)*three+r3(7,40)*six &
&             +r1(2,18)*three
      rxyz(20)=+r6(17,11)*four+r4(8,28)*p24+r2(5,19)*p12
      eri(1,1,2,2)=r400+(+r7(11,3)*two+r7(11,4)*two+r5(4,11)*p12+r5(4,12)*p12 &
&                  +r5(11,21)*six+r5(11,22)*six+r3(1,9)*six+r3(1,10)*six+r3(4,39)*p36 &
&                  +r3(4,40)*p36+r1(1,17)*p18+r1(1,18)*p18)*qx+(+r6(11,10)+r6(11,11)*four &
&                  +r6(11,12)+r4(4,27)*six+r4(4,28)*p24+r4(4,29)*six+r4(11,39)+r4(11,40)*four &
&                  +r4(11,41)+r2(1,18)*three+r2(1,19)*p12+r2(1,20)*three+r2(2,53)*six &
&                  +r2(2,54)*p24+r2(2,55)*six+r0(22)*three+r0(23)*p12+r0(24)*three)*xx+( &
&                  +r5(11,23)*two+r5(11,24)*two+r3(4,41)*p12+r3(4,42)*p12+r1(1,19)*six &
&                  +r1(1,20)*six)*xxx+rxyz(1)*xxxx
      eri(2,1,2,2)=r220+(+r7(22,4)*two+r5(11,12)*p12+r5(11,22)*two+r3(4,10)*six &
&                  +r3(4,40)*p12+r1(1,18)*six)*qx+rxyz(5)*xx
      eri(3,1,2,2)=r202+(+r7(24,4)*two+r5(13,12)*p12+r5(11,22)*two+r3(6,10)*six &
&                  +r3(4,40)*p12+r1(1,18)*six)*qx+(+r7(17,3)*two+r5(8,11)*p12+r5(17,21)*two &
&                  +r3(3,9)*six+r3(8,39)*p12+r1(3,17)*six)*qz+(+r6(24,12)+r4(13,29)*six &
&                  +r4(11,41)+r2(3,20)*three+r2(2,55)*six+r0(24)*three)*xx+rxyz(20)*xz+( &
&                  +r6(11,10)+r4(4,27)*six+r4(11,39)+r2(1,18)*three+r2(2,53)*six+r0(22)*three &
&                  )*zz+(+r5(17,24)*two+r3(8,42)*p12+r1(3,20)*six)*xxz+(+r5(11,23)*two &
&                  +r3(4,41)*p12+r1(1,19)*six)*xzz+rxyz(1)*xxzz
      eri(4,1,2,2)=r310+(+r7(16,3)+r7(16,4)*two+r5(7,11)*six+r5(7,12)*p12+r5(16,21) &
&                  +r5(16,22)*two+r3(2,9)*three+r3(2,10)*six+r3(7,39)*six+r3(7,40)*p12 &
&                  +r1(2,17)*three+r1(2,18)*six)*qx+(+r6(16,11)*two+r6(16,12)+r4(7,28)*p12 &
&                  +r4(7,29)*six+r2(4,19)*six+r2(4,20)*three)*xx+rxyz(2)*xxx
      eri(5,1,2,2)=r301+(+r7(17,3)+r7(17,4)*two+r5(8,11)*six+r5(8,12)*p12+r5(17,21) &
&                  +r5(17,22)*two+r3(3,9)*three+r3(3,10)*six+r3(8,39)*six+r3(8,40)*p12 &
&                  +r1(3,17)*three+r1(3,18)*six)*qx+(+r7(11,3)+r5(4,11)*six+r5(11,21)*three &
&                  +r3(1,9)*three+r3(4,39)*p18+r1(1,17)*nine)*qz+(+r6(17,11)*two+r6(17,12) &
&                  +r4(8,28)*p12+r4(8,29)*six+r2(5,19)*six+r2(5,20)*three)*xx+(+r6(11,10) &
&                  +r6(11,11)*two+r4(4,27)*six+r4(4,28)*p12+r4(11,39)+r4(11,40)*two &
&                  +r2(1,18)*three+r2(1,19)*six+r2(2,53)*six+r2(2,54)*p12+r0(22)*three &
&                  +r0(23)*six)*xz+(+r5(17,24)+r3(8,42)*six+r1(3,20)*three)*xxx+( &
&                  +r5(11,23)*two+r5(11,24)+r3(4,41)*p12+r3(4,42)*six+r1(1,19)*six &
&                  +r1(1,20)*three)*xxz+rxyz(1)*xxxz
      eri(6,1,2,2)=r211+(+r7(23,4)*two+r5(12,12)*p12+r3(5,10)*six)*qx+rxyz(17)*qz+( &
&                  +r6(23,12)+r4(12,29)*six+r2(6,20)*three)*xx+(+r6(16,11)*two+r4(7,28)*p12 &
&                  +r2(4,19)*six)*xz+rxyz(2)*xxz
      eri(1,2,2,2)=r220+(+r7(22,3)*two+r5(11,11)*p12+r5(11,21)*two+r3(4,9)*six &
&                  +r3(4,39)*p12+r1(1,17)*six)*qx+rxyz(4)*xx
      eri(2,2,2,2)=r040
      eri(3,2,2,2)=r022+(+r7(30,3)*two+r5(17,11)*p12+r5(17,21)*two+r3(8,9)*six &
&                  +r3(8,39)*p12+r1(3,17)*six)*qz+rxyz(4)*zz
      eri(4,2,2,2)=r130+rxyz(7)*qx
      eri(5,2,2,2)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,2,2)=r031+rxyz(7)*qz
      eri(1,3,2,2)=r202+(+r7(24,3)*two+r5(13,11)*p12+r5(11,21)*two+r3(6,9)*six &
&                  +r3(4,39)*p12+r1(1,17)*six)*qx+(+r7(17,4)*two+r5(8,12)*p12+r5(17,22)*two &
&                  +r3(3,10)*six+r3(8,40)*p12+r1(3,18)*six)*qz+(+r6(24,10)+r4(13,27)*six &
&                  +r4(11,39)+r2(3,18)*three+r2(2,53)*six+r0(22)*three)*xx+rxyz(20)*xz+( &
&                  +r6(11,12)+r4(4,29)*six+r4(11,41)+r2(1,20)*three+r2(2,55)*six+r0(24)*three &
&                  )*zz+(+r5(17,23)*two+r3(8,41)*p12+r1(3,19)*six)*xxz+(+r5(11,24)*two &
&                  +r3(4,42)*p12+r1(1,20)*six)*xzz+rxyz(1)*xxzz
      eri(2,3,2,2)=r022+(+r7(30,4)*two+r5(17,12)*p12+r5(17,22)*two+r3(8,10)*six &
&                  +r3(8,40)*p12+r1(3,18)*six)*qz+rxyz(5)*zz
      eri(3,3,2,2)=r004+(+r7(32,3)*two+r7(32,4)*two+r5(19,11)*p12+r5(19,12)*p12 &
&                  +r5(17,21)*six+r5(17,22)*six+r3(10,9)*six+r3(10,10)*six+r3(8,39)*p36 &
&                  +r3(8,40)*p36+r1(3,17)*p18+r1(3,18)*p18)*qz+(+r6(24,10)+r6(24,11)*four &
&                  +r6(24,12)+r4(13,27)*six+r4(13,28)*p24+r4(13,29)*six+r4(11,39) &
&                  +r4(11,40)*four+r4(11,41)+r2(3,18)*three+r2(3,19)*p12+r2(3,20)*three &
&                  +r2(2,53)*six+r2(2,54)*p24+r2(2,55)*six+r0(22)*three+r0(23)*p12 &
&                  +r0(24)*three)*zz+(+r5(17,23)*two+r5(17,24)*two+r3(8,41)*p12+r3(8,42)*p12 &
&                  +r1(3,19)*six+r1(3,20)*six)*zzz+rxyz(1)*zzzz
      eri(4,3,2,2)=r112+rxyz(18)*qx+(+r7(23,4)*two+r5(12,12)*p12+r3(5,10)*six)*qz+( &
&                  +r6(23,11)*two+r4(12,28)*p12+r2(6,19)*six)*xz+(+r6(16,12)+r4(7,29)*six &
&                  +r2(4,20)*three)*zz+rxyz(2)*xzz
      eri(5,3,2,2)=r103+(+r7(32,3)+r5(19,11)*six+r5(17,21)*three+r3(10,9)*three &
&                  +r3(8,39)*p18+r1(3,17)*nine)*qx+(+r7(24,3)+r7(24,4)*two+r5(13,11)*six &
&                  +r5(13,12)*p12+r5(11,21)+r5(11,22)*two+r3(6,9)*three+r3(6,10)*six &
&                  +r3(4,39)*six+r3(4,40)*p12+r1(1,17)*three+r1(1,18)*six)*qz+(+r6(24,10) &
&                  +r6(24,11)*two+r4(13,27)*six+r4(13,28)*p12+r4(11,39)+r4(11,40)*two &
&                  +r2(3,18)*three+r2(3,19)*six+r2(2,53)*six+r2(2,54)*p12+r0(22)*three &
&                  +r0(23)*six)*xz+(+r6(17,11)*two+r6(17,12)+r4(8,28)*p12+r4(8,29)*six &
&                  +r2(5,19)*six+r2(5,20)*three)*zz+(+r5(17,23)*two+r5(17,24)+r3(8,41)*p12 &
&                  +r3(8,42)*six+r1(3,19)*six+r1(3,20)*three)*xzz+(+r5(11,24)+r3(4,42)*six &
&                  +r1(1,20)*three)*zzz+rxyz(1)*xzzz
      eri(6,3,2,2)=r013+(+r7(31,3)+r7(31,4)*two+r5(18,11)*six+r5(18,12)*p12+r5(16,21) &
&                  +r5(16,22)*two+r3(9,9)*three+r3(9,10)*six+r3(7,39)*six+r3(7,40)*p12 &
&                  +r1(2,17)*three+r1(2,18)*six)*qz+(+r6(23,11)*two+r6(23,12)+r4(12,28)*p12 &
&                  +r4(12,29)*six+r2(6,19)*six+r2(6,20)*three)*zz+rxyz(2)*zzz
      eri(1,4,2,2)=r310+(+r7(16,3)*two+r7(16,4)+r5(7,11)*p12+r5(7,12)*six+r5(16,21)*two &
&                  +r5(16,22)+r3(2,9)*six+r3(2,10)*three+r3(7,39)*p12+r3(7,40)*six &
&                  +r1(2,17)*six+r1(2,18)*three)*qx+(+r6(16,10)+r6(16,11)*two+r4(7,27)*six &
&                  +r4(7,28)*p12+r2(4,18)*three+r2(4,19)*six)*xx+rxyz(3)*xxx
      eri(2,4,2,2)=r130+rxyz(8)*qx
      eri(3,4,2,2)=r112+rxyz(19)*qx+(+r7(23,3)*two+r5(12,11)*p12+r3(5,9)*six)*qz+( &
&                  +r6(23,11)*two+r4(12,28)*p12+r2(6,19)*six)*xz+(+r6(16,10)+r4(7,27)*six &
&                  +r2(4,18)*three)*zz+rxyz(3)*xzz
      eri(4,4,2,2)=r220+(+r7(22,3)+r7(22,4)+r5(11,11)*six+r5(11,12)*six+r5(11,21) &
&                  +r5(11,22)+r3(4,9)*three+r3(4,10)*three+r3(4,39)*six+r3(4,40)*six &
&                  +r1(1,17)*three+r1(1,18)*three)*qx+rxyz(6)*xx
      eri(5,4,2,2)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(16,10)+r6(16,11) &
&                  +r4(7,27)*six+r4(7,28)*six+r2(4,18)*three+r2(4,19)*three)*xz+rxyz(3)*xxz
      eri(6,4,2,2)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,2,2)=r301+(+r7(17,3)*two+r7(17,4)+r5(8,11)*p12+r5(8,12)*six+r5(17,21)*two &
&                  +r5(17,22)+r3(3,9)*six+r3(3,10)*three+r3(8,39)*p12+r3(8,40)*six &
&                  +r1(3,17)*six+r1(3,18)*three)*qx+(+r7(11,4)+r5(4,12)*six+r5(11,22)*three &
&                  +r3(1,10)*three+r3(4,40)*p18+r1(1,18)*nine)*qz+(+r6(17,10)+r6(17,11)*two &
&                  +r4(8,27)*six+r4(8,28)*p12+r2(5,18)*three+r2(5,19)*six)*xx+(+r6(11,11)*two &
&                  +r6(11,12)+r4(4,28)*p12+r4(4,29)*six+r4(11,40)*two+r4(11,41)+r2(1,19)*six &
&                  +r2(1,20)*three+r2(2,54)*p12+r2(2,55)*six+r0(23)*six+r0(24)*three)*xz+( &
&                  +r5(17,23)+r3(8,41)*six+r1(3,19)*three)*xxx+(+r5(11,23)+r5(11,24)*two &
&                  +r3(4,41)*six+r3(4,42)*p12+r1(1,19)*three+r1(1,20)*six)*xxz+rxyz(1)*xxxz
      eri(2,5,2,2)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,2,2)=r103+(+r7(32,4)+r5(19,12)*six+r5(17,22)*three+r3(10,10)*three &
&                  +r3(8,40)*p18+r1(3,18)*nine)*qx+(+r7(24,3)*two+r7(24,4)+r5(13,11)*p12 &
&                  +r5(13,12)*six+r5(11,21)*two+r5(11,22)+r3(6,9)*six+r3(6,10)*three &
&                  +r3(4,39)*p12+r3(4,40)*six+r1(1,17)*six+r1(1,18)*three)*qz+(+r6(24,11)*two &
&                  +r6(24,12)+r4(13,28)*p12+r4(13,29)*six+r4(11,40)*two+r4(11,41) &
&                  +r2(3,19)*six+r2(3,20)*three+r2(2,54)*p12+r2(2,55)*six+r0(23)*six &
&                  +r0(24)*three)*xz+(+r6(17,10)+r6(17,11)*two+r4(8,27)*six+r4(8,28)*p12 &
&                  +r2(5,18)*three+r2(5,19)*six)*zz+(+r5(17,23)+r5(17,24)*two+r3(8,41)*six &
&                  +r3(8,42)*p12+r1(3,19)*three+r1(3,20)*six)*xzz+(+r5(11,23)+r3(4,41)*six &
&                  +r1(1,19)*three)*zzz+rxyz(1)*xzzz
      eri(4,5,2,2)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(16,11)+r6(16,12) &
&                  +r4(7,28)*six+r4(7,29)*six+r2(4,19)*three+r2(4,20)*three)*xz+rxyz(2)*xxz
      eri(5,5,2,2)=r202+(+r7(24,3)+r7(24,4)+r5(13,11)*six+r5(13,12)*six+r5(11,21) &
&                  +r5(11,22)+r3(6,9)*three+r3(6,10)*three+r3(4,39)*six+r3(4,40)*six &
&                  +r1(1,17)*three+r1(1,18)*three)*qx+(+r7(17,3)+r7(17,4)+r5(8,11)*six &
&                  +r5(8,12)*six+r5(17,21)+r5(17,22)+r3(3,9)*three+r3(3,10)*three &
&                  +r3(8,39)*six+r3(8,40)*six+r1(3,17)*three+r1(3,18)*three)*qz+(+r6(24,11) &
&                  +r4(13,28)*six+r4(11,40)+r2(3,19)*three+r2(2,54)*six+r0(23)*three)*xx+( &
&                  +r6(17,10)+r6(17,11)*two+r6(17,12)+r4(8,27)*six+r4(8,28)*p12+r4(8,29)*six &
&                  +r2(5,18)*three+r2(5,19)*six+r2(5,20)*three)*xz+(+r6(11,11)+r4(4,28)*six &
&                  +r4(11,40)+r2(1,19)*three+r2(2,54)*six+r0(23)*three)*zz+(+r5(17,23) &
&                  +r5(17,24)+r3(8,41)*six+r3(8,42)*six+r1(3,19)*three+r1(3,20)*three)*xxz+( &
&                  +r5(11,23)+r5(11,24)+r3(4,41)*six+r3(4,42)*six+r1(1,19)*three &
&                  +r1(1,20)*three)*xzz+rxyz(1)*xxzz
      eri(6,5,2,2)=r112+rxyz(19)*qx+(+r7(23,3)+r7(23,4)+r5(12,11)*six+r5(12,12)*six &
&                  +r3(5,9)*three+r3(5,10)*three)*qz+(+r6(23,11)+r6(23,12)+r4(12,28)*six &
&                  +r4(12,29)*six+r2(6,19)*three+r2(6,20)*three)*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,2,2)=r211+(+r7(23,3)*two+r5(12,11)*p12+r3(5,9)*six)*qx+rxyz(16)*qz+( &
&                  +r6(23,10)+r4(12,27)*six+r2(6,18)*three)*xx+(+r6(16,11)*two+r4(7,28)*p12 &
&                  +r2(4,19)*six)*xz+rxyz(3)*xxz
      eri(2,6,2,2)=r031+rxyz(8)*qz
      eri(3,6,2,2)=r013+(+r7(31,3)*two+r7(31,4)+r5(18,11)*p12+r5(18,12)*six+r5(16,21)*two &
&                  +r5(16,22)+r3(9,9)*six+r3(9,10)*three+r3(7,39)*p12+r3(7,40)*six &
&                  +r1(2,17)*six+r1(2,18)*three)*qz+(+r6(23,10)+r6(23,11)*two+r4(12,27)*six &
&                  +r4(12,28)*p12+r2(6,18)*three+r2(6,19)*six)*zz+rxyz(3)*zzz
      eri(4,6,2,2)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,2,2)=r112+rxyz(18)*qx+(+r7(23,3)+r7(23,4)+r5(12,11)*six+r5(12,12)*six &
&                  +r3(5,9)*three+r3(5,10)*three)*qz+(+r6(23,10)+r6(23,11)+r4(12,27)*six &
&                  +r4(12,28)*six+r2(6,18)*three+r2(6,19)*three)*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,2,2)=r022+(+r7(30,3)+r7(30,4)+r5(17,11)*six+r5(17,12)*six+r5(17,21) &
&                  +r5(17,22)+r3(8,9)*three+r3(8,10)*three+r3(8,39)*six+r3(8,40)*six &
&                  +r1(3,17)*three+r1(3,18)*three)*qz+rxyz(6)*zz
!
      r400= r8(13)-r7(8,1)*two+r6(4,1)+r6(4,4)+r6(6,4)+r6(13,9)*six-r5(3,3)*two &
&          -r5(8,13)*p12+r4(1,2)+r4(1,5)+r4(4,14)*six+r4(4,26)*six+r4(6,26)*six &
&          +r4(13,38)*three-r3(3,19)*p12-r3(8,43)*six+r2(1,5)*six+r2(1,17)*six &
&          +r2(2,37)*three+r2(2,52)*three+r2(3,52)*three-r1(3,31)*six+r0(6)*three &
&          +r0(21)*three
      r310= r8(18)-r7(12,1)*two+r6(7,1)+r6(7,4)+r6(9,4)+r6(18,9)*three-r5(5,3)*two &
&          -r5(12,13)*six+r4(2,2)+r4(2,5)+r4(7,14)*three+r4(7,26)*three+r4(9,26)*three &
&          -r3(5,19)*six+r2(4,5)*three+r2(4,17)*three
      r301= r8(19)-r7(13,1)*two+r6(8,1)+r6(8,4)+r6(10,4)+r6(19,9)*three-r5(6,3)*two &
&          -r5(13,13)*six+r4(3,2)+r4(3,5)+r4(8,14)*three+r4(8,26)*three+r4(10,26)*three &
&          -r3(6,19)*six+r2(5,5)*three+r2(5,17)*three
      r220= r8(24)-r7(17,1)*two+r6(11,1)+r6(11,4)+r6(13,4)+r6(13,9)+r6(24,9)-r5(8,3)*two &
&          -r5(8,13)*two-r5(17,13)*two+r4(4,2)+r4(4,5)+r4(4,14)+r4(11,14)+r4(4,26) &
&          +r4(6,26)+r4(11,26)+r4(13,26)+r4(13,38)-r3(3,19)*two-r3(8,19)*two-r3(8,43)*two &
&          +r2(1,5)+r2(2,5)+r2(1,17)+r2(2,17)+r2(2,37)+r2(2,52)+r2(3,52)-r1(3,31)*two &
&          +r0(6)+r0(21)
      r211= r8(25)-r7(18,1)*two+r6(12,1)+r6(12,4)+r6(14,4)+r6(25,9)-r5(9,3)*two &
&          -r5(18,13)*two+r4(5,2)+r4(5,5)+r4(12,14)+r4(12,26)+r4(14,26)-r3(9,19)*two &
&          +r2(6,5)+r2(6,17)
      r202= r8(26)-r7(19,1)*two+r6(13,1)+r6(13,4)+r6(15,4)+r6(13,9)+r6(26,9)-r5(10,3)*two &
&          -r5(8,13)*two-r5(19,13)*two+r4(6,2)+r4(6,5)+r4(4,14)+r4(13,14)+r4(4,26) &
&          +r4(6,26)+r4(13,26)+r4(15,26)+r4(13,38)-r3(3,19)*two-r3(10,19)*two-r3(8,43)*two &
&          +r2(1,5)+r2(3,5)+r2(1,17)+r2(3,17)+r2(2,37)+r2(2,52)+r2(3,52)-r1(3,31)*two &
&          +r0(6)+r0(21)
      r130= r8(31)-r7(23,1)*two+r6(16,1)+r6(16,4)+r6(18,4)+r6(18,9)*three-r5(12,3)*two &
&          -r5(12,13)*six+r4(7,2)+r4(7,5)+r4(7,14)*three+r4(7,26)*three+r4(9,26)*three &
&          -r3(5,19)*six+r2(4,5)*three+r2(4,17)*three
      r121= r8(32)-r7(24,1)*two+r6(17,1)+r6(17,4)+r6(19,4)+r6(19,9)-r5(13,3)*two &
&          -r5(13,13)*two+r4(8,2)+r4(8,5)+r4(8,14)+r4(8,26)+r4(10,26)-r3(6,19)*two+r2(5,5) &
&          +r2(5,17)
      r112= r8(33)-r7(25,1)*two+r6(18,1)+r6(18,4)+r6(20,4)+r6(18,9)-r5(14,3)*two &
&          -r5(12,13)*two+r4(9,2)+r4(9,5)+r4(7,14)+r4(7,26)+r4(9,26)-r3(5,19)*two+r2(4,5) &
&          +r2(4,17)
      r103= r8(34)-r7(26,1)*two+r6(19,1)+r6(19,4)+r6(21,4)+r6(19,9)*three-r5(15,3)*two &
&          -r5(13,13)*six+r4(10,2)+r4(10,5)+r4(8,14)*three+r4(8,26)*three+r4(10,26)*three &
&          -r3(6,19)*six+r2(5,5)*three+r2(5,17)*three
      r040= r8(39)-r7(30,1)*two+r6(22,1)+r6(22,4)+r6(24,4)+r6(24,9)*six-r5(17,3)*two &
&          -r5(17,13)*p12+r4(11,2)+r4(11,5)+r4(11,14)*six+r4(11,26)*six+r4(13,26)*six &
&          +r4(13,38)*three-r3(8,19)*p12-r3(8,43)*six+r2(2,5)*six+r2(2,17)*six &
&          +r2(2,37)*three+r2(2,52)*three+r2(3,52)*three-r1(3,31)*six+r0(6)*three &
&          +r0(21)*three
      r031= r8(40)-r7(31,1)*two+r6(23,1)+r6(23,4)+r6(25,4)+r6(25,9)*three-r5(18,3)*two &
&          -r5(18,13)*six+r4(12,2)+r4(12,5)+r4(12,14)*three+r4(12,26)*three &
&          +r4(14,26)*three-r3(9,19)*six+r2(6,5)*three+r2(6,17)*three
      r022= r8(41)-r7(32,1)*two+r6(24,1)+r6(24,4)+r6(26,4)+r6(24,9)+r6(26,9)-r5(19,3)*two &
&          -r5(17,13)*two-r5(19,13)*two+r4(13,2)+r4(13,5)+r4(11,14)+r4(13,14)+r4(11,26) &
&          +r4(13,26)*two+r4(15,26)+r4(13,38)-r3(8,19)*two-r3(10,19)*two-r3(8,43)*two &
&          +r2(2,5)+r2(3,5)+r2(2,17)+r2(3,17)+r2(2,37)+r2(2,52)+r2(3,52)-r1(3,31)*two &
&          +r0(6)+r0(21)
      r013= r8(42)-r7(33,1)*two+r6(25,1)+r6(25,4)+r6(27,4)+r6(25,9)*three-r5(20,3)*two &
&          -r5(18,13)*six+r4(14,2)+r4(14,5)+r4(12,14)*three+r4(12,26)*three &
&          +r4(14,26)*three-r3(9,19)*six+r2(6,5)*three+r2(6,17)*three
      r004= r8(43)-r7(34,1)*two+r6(26,1)+r6(26,4)+r6(28,4)+r6(26,9)*six-r5(21,3)*two &
&          -r5(19,13)*p12+r4(15,2)+r4(15,5)+r4(13,14)*six+r4(13,26)*six+r4(15,26)*six &
&          +r4(13,38)*three-r3(10,19)*p12-r3(8,43)*six+r2(3,5)*six+r2(3,17)*six &
&          +r2(2,37)*three+r2(2,52)*three+r2(3,52)*three-r1(3,31)*six+r0(6)*three &
&          +r0(21)*three
      rxyz(1)=+r4(13,42)-r3(8,47)*two+r2(2,41)+r2(2,56)+r2(3,56)-r1(3,35)*two+r0(10) &
&             +r0(25)
      rxyz(2)=+r5(18,24)-r4(12,33)*two+r3(7,30)+r3(7,42)+r3(9,42)-r2(6,32)*two+r1(2,8) &
&             +r1(2,20)
      rxyz(3)=+r5(18,23)-r4(12,32)*two+r3(7,29)+r3(7,41)+r3(9,41)-r2(6,31)*two+r1(2,7) &
&             +r1(2,19)
      rxyz(4)=+r6(24,10)-r5(17,14)*two+r4(11,15)+r4(11,27)+r4(13,27)+r4(13,39) &
&             -r3(8,20)*two-r3(8,44)*two+r2(2,6)+r2(2,18)+r2(2,38)+r2(2,53)+r2(3,53) &
&             -r1(3,32)*two+r0(7)+r0(22)
      rxyz(5)=+r6(24,12)-r5(17,16)*two+r4(11,17)+r4(11,29)+r4(13,29)+r4(13,41) &
&             -r3(8,22)*two-r3(8,46)*two+r2(2,8)+r2(2,20)+r2(2,40)+r2(2,55)+r2(3,55) &
&             -r1(3,34)*two+r0(9)+r0(24)
      rxyz(6)=+r6(24,11)-r5(17,15)*two+r4(11,16)+r4(11,28)+r4(13,28)+r4(13,40) &
&             -r3(8,21)*two-r3(8,45)*two+r2(2,7)+r2(2,19)+r2(2,39)+r2(2,54)+r2(3,54) &
&             -r1(3,33)*two+r0(8)+r0(23)
      rxyz(7)=+r7(31,3)-r6(23,5)*two+r5(16,5)+r5(16,11)+r5(18,11)+r5(18,21)*three &
&             -r4(12,10)*two-r4(12,30)*six+r3(7,3)+r3(7,9)+r3(7,27)*three+r3(7,39)*three &
&             +r3(9,39)*three-r2(6,29)*six+r1(2,5)*three+r1(2,17)*three
      rxyz(8)=+r7(31,4)-r6(23,6)*two+r5(16,6)+r5(16,12)+r5(18,12)+r5(18,22)*three &
&             -r4(12,11)*two-r4(12,31)*six+r3(7,4)+r3(7,10)+r3(7,28)*three+r3(7,40)*three &
&             +r3(9,40)*three-r2(6,30)*six+r1(2,6)*three+r1(2,18)*three
      rxyz(9)=+r7(25,3)+r7(25,4)-r6(18,5)*two-r6(18,6)*two+r5(12,5)+r5(12,6)+r5(12,11) &
&             +r5(14,11)+r5(12,12)+r5(14,12)-r4(9,10)*two-r4(9,11)*two+r3(5,3)+r3(5,4) &
&             +r3(5,9)+r3(5,10)
      rxyz(10)=+r6(25,11)-r5(18,15)*two+r4(12,16)+r4(12,28)+r4(14,28)-r3(9,21)*two+r2(6,7) &
&             +r2(6,19)
      rxyz(11)=+r6(18,11)-r5(12,15)*two+r4(7,16)+r4(7,28)+r4(9,28)-r3(5,21)*two+r2(4,7) &
&             +r2(4,19)
      rxyz(12)=+r7(32,3)-r6(24,5)*two+r5(17,5)+r5(17,11)+r5(19,11)+r5(19,21)-r4(13,10)*two &
&             -r4(13,30)*two+r3(8,3)+r3(8,9)+r3(8,27)+r3(8,39)+r3(10,39)-r2(3,29)*two &
&             +r1(3,5)+r1(3,17)
      rxyz(13)=+r7(32,4)-r6(24,6)*two+r5(17,6)+r5(17,12)+r5(19,12)+r5(19,22)-r4(13,11)*two &
&             -r4(13,31)*two+r3(8,4)+r3(8,10)+r3(8,28)+r3(8,40)+r3(10,40)-r2(3,30)*two &
&             +r1(3,6)+r1(3,18)
      rxyz(14)=+r7(24,4)-r6(17,6)*two+r5(11,6)+r5(11,12)+r5(13,12)+r5(13,22)-r4(8,11)*two &
&             -r4(8,31)*two+r3(4,4)+r3(4,10)+r3(4,28)+r3(4,40)+r3(6,40)-r2(5,30)*two &
&             +r1(1,6)+r1(1,18)
      rxyz(15)=+r7(24,3)-r6(17,5)*two+r5(11,5)+r5(11,11)+r5(13,11)+r5(13,21)-r4(8,10)*two &
&             -r4(8,30)*two+r3(4,3)+r3(4,9)+r3(4,27)+r3(4,39)+r3(6,39)-r2(5,29)*two &
&             +r1(1,5)+r1(1,17)
      rxyz(16)=+r7(18,4)-r6(12,6)*two+r5(7,6)+r5(7,12)+r5(9,12)+r5(18,22)-r4(5,11)*two &
&             -r4(12,31)*two+r3(2,4)+r3(2,10)+r3(7,28)+r3(7,40)+r3(9,40)-r2(6,30)*two &
&             +r1(2,6)+r1(2,18)
      rxyz(17)=+r7(18,3)-r6(12,5)*two+r5(7,5)+r5(7,11)+r5(9,11)+r5(18,21)-r4(5,10)*two &
&             -r4(12,30)*two+r3(2,3)+r3(2,9)+r3(7,27)+r3(7,39)+r3(9,39)-r2(6,29)*two &
&             +r1(2,5)+r1(2,17)
      rxyz(18)=+r7(33,3)-r6(25,5)*two+r5(18,5)+r5(18,11)+r5(20,11)+r5(18,21)-r4(14,10)*two &
&             -r4(12,30)*two+r3(9,3)+r3(9,9)+r3(7,27)+r3(7,39)+r3(9,39)-r2(6,29)*two &
&             +r1(2,5)+r1(2,17)
      rxyz(19)=+r7(33,4)-r6(25,6)*two+r5(18,6)+r5(18,12)+r5(20,12)+r5(18,22)-r4(14,11)*two &
&             -r4(12,31)*two+r3(9,4)+r3(9,10)+r3(7,28)+r3(7,40)+r3(9,40)-r2(6,30)*two &
&             +r1(2,6)+r1(2,18)
      rxyz(20)=+r6(19,11)*four-r5(13,15)*eight+r4(8,16)*four+r4(8,28)*four+r4(10,28)*four &
&             -r3(6,21)*eight+r2(5,7)*four+r2(5,19)*four
      eri(1,1,3,2)=r400+(+r7(13,3)*two+r7(13,4)*two-r6(8,5)*four-r6(8,6)*four+r5(4,5)*two &
&                  +r5(4,6)*two+r5(4,11)*two+r5(6,11)*two+r5(4,12)*two+r5(6,12)*two &
&                  +r5(13,21)*six+r5(13,22)*six-r4(3,10)*four-r4(3,11)*four-r4(8,30)*p12 &
&                  -r4(8,31)*p12+r3(1,3)*two+r3(1,4)*two+r3(1,9)*two+r3(1,10)*two &
&                  +r3(4,27)*six+r3(4,28)*six+r3(4,39)*six+r3(6,39)*six+r3(4,40)*six &
&                  +r3(6,40)*six-r2(5,29)*p12-r2(5,30)*p12+r1(1,5)*six+r1(1,6)*six &
&                  +r1(1,17)*six+r1(1,18)*six)*qx+(+r6(13,10)+r6(13,11)*four+r6(13,12) &
&                  -r5(8,14)*two-r5(8,15)*eight-r5(8,16)*two+r4(4,15)+r4(4,16)*four+r4(4,17) &
&                  +r4(4,27)+r4(6,27)+r4(4,28)*four+r4(6,28)*four+r4(4,29)+r4(6,29)+r4(13,39) &
&                  +r4(13,40)*four+r4(13,41)-r3(3,20)*two-r3(3,21)*eight-r3(3,22)*two &
&                  -r3(8,44)*two-r3(8,45)*eight-r3(8,46)*two+r2(1,6)+r2(1,7)*four+r2(1,8) &
&                  +r2(1,18)+r2(1,19)*four+r2(1,20)+r2(2,38)+r2(2,39)*four+r2(2,40)+r2(2,53) &
&                  +r2(3,53)+r2(2,54)*four+r2(3,54)*four+r2(2,55)+r2(3,55)-r1(3,32)*two &
&                  -r1(3,33)*eight-r1(3,34)*two+r0(7)+r0(8)*four+r0(9)+r0(22)+r0(23)*four &
&                  +r0(24))*xx+(+r5(13,23)*two+r5(13,24)*two-r4(8,32)*four-r4(8,33)*four &
&                  +r3(4,29)*two+r3(4,30)*two+r3(4,41)*two+r3(6,41)*two+r3(4,42)*two &
&                  +r3(6,42)*two-r2(5,31)*four-r2(5,32)*four+r1(1,7)*two+r1(1,8)*two &
&                  +r1(1,19)*two+r1(1,20)*two)*xxx+rxyz(1)*xxxx
      eri(2,1,3,2)=r220+(+r7(24,4)*two-r6(17,6)*four+r5(11,6)*two+r5(11,12)*two &
&                  +r5(13,12)*two+r5(13,22)*two-r4(8,11)*four-r4(8,31)*four+r3(4,4)*two &
&                  +r3(4,10)*two+r3(4,28)*two+r3(4,40)*two+r3(6,40)*two-r2(5,30)*four &
&                  +r1(1,6)*two+r1(1,18)*two)*qx+rxyz(5)*xx
      eri(3,1,3,2)=r202+(+r7(26,4)*two-r6(19,6)*four+r5(13,6)*two+r5(13,12)*two &
&                  +r5(15,12)*two+r5(13,22)*two-r4(10,11)*four-r4(8,31)*four+r3(6,4)*two &
&                  +r3(6,10)*two+r3(4,28)*two+r3(4,40)*two+r3(6,40)*two-r2(5,30)*four &
&                  +r1(1,6)*two+r1(1,18)*two)*qx+(+r7(19,3)*two-r6(13,5)*four+r5(8,5)*two &
&                  +r5(8,11)*two+r5(10,11)*two+r5(19,21)*two-r4(6,10)*four-r4(13,30)*four &
&                  +r3(3,3)*two+r3(3,9)*two+r3(8,27)*two+r3(8,39)*two+r3(10,39)*two &
&                  -r2(3,29)*four+r1(3,5)*two+r1(3,17)*two)*qz+(+r6(26,12)-r5(19,16)*two &
&                  +r4(13,17)+r4(13,29)+r4(15,29)+r4(13,41)-r3(10,22)*two-r3(8,46)*two &
&                  +r2(3,8)+r2(3,20)+r2(2,40)+r2(2,55)+r2(3,55)-r1(3,34)*two+r0(9)+r0(24))*xx &
&                  +rxyz(20)*xz+(+r6(13,10)-r5(8,14)*two+r4(4,15)+r4(4,27)+r4(6,27)+r4(13,39) &
&                  -r3(3,20)*two-r3(8,44)*two+r2(1,6)+r2(1,18)+r2(2,38)+r2(2,53)+r2(3,53) &
&                  -r1(3,32)*two+r0(7)+r0(22))*zz+(+r5(19,24)*two-r4(13,33)*four+r3(8,30)*two &
&                  +r3(8,42)*two+r3(10,42)*two-r2(3,32)*four+r1(3,8)*two+r1(3,20)*two)*xxz+( &
&                  +r5(13,23)*two-r4(8,32)*four+r3(4,29)*two+r3(4,41)*two+r3(6,41)*two &
&                  -r2(5,31)*four+r1(1,7)*two+r1(1,19)*two)*xzz+rxyz(1)*xxzz
      eri(4,1,3,2)=r310+(+r7(18,3)+r7(18,4)*two-r6(12,5)*two-r6(12,6)*four+r5(7,5) &
&                  +r5(7,6)*two+r5(7,11)+r5(9,11)+r5(7,12)*two+r5(9,12)*two+r5(18,21) &
&                  +r5(18,22)*two-r4(5,10)*two-r4(5,11)*four-r4(12,30)*two-r4(12,31)*four &
&                  +r3(2,3)+r3(2,4)*two+r3(2,9)+r3(2,10)*two+r3(7,27)+r3(7,28)*two+r3(7,39) &
&                  +r3(9,39)+r3(7,40)*two+r3(9,40)*two-r2(6,29)*two-r2(6,30)*four+r1(2,5) &
&                  +r1(2,6)*two+r1(2,17)+r1(2,18)*two)*qx+(+r6(18,11)*two+r6(18,12) &
&                  -r5(12,15)*four-r5(12,16)*two+r4(7,16)*two+r4(7,17)+r4(7,28)*two &
&                  +r4(9,28)*two+r4(7,29)+r4(9,29)-r3(5,21)*four-r3(5,22)*two+r2(4,7)*two &
&                  +r2(4,8)+r2(4,19)*two+r2(4,20))*xx+rxyz(2)*xxx
      eri(5,1,3,2)=r301+(+r7(19,3)+r7(19,4)*two-r6(13,5)*two-r6(13,6)*four+r5(8,5) &
&                  +r5(8,6)*two+r5(8,11)+r5(10,11)+r5(8,12)*two+r5(10,12)*two+r5(19,21) &
&                  +r5(19,22)*two-r4(6,10)*two-r4(6,11)*four-r4(13,30)*two-r4(13,31)*four &
&                  +r3(3,3)+r3(3,4)*two+r3(3,9)+r3(3,10)*two+r3(8,27)+r3(8,28)*two+r3(8,39) &
&                  +r3(10,39)+r3(8,40)*two+r3(10,40)*two-r2(3,29)*two-r2(3,30)*four+r1(3,5) &
&                  +r1(3,6)*two+r1(3,17)+r1(3,18)*two)*qx+(+r7(13,3)-r6(8,5)*two+r5(4,5) &
&                  +r5(4,11)+r5(6,11)+r5(13,21)*three-r4(3,10)*two-r4(8,30)*six+r3(1,3) &
&                  +r3(1,9)+r3(4,27)*three+r3(4,39)*three+r3(6,39)*three-r2(5,29)*six &
&                  +r1(1,5)*three+r1(1,17)*three)*qz+(+r6(19,11)*two+r6(19,12)-r5(13,15)*four &
&                  -r5(13,16)*two+r4(8,16)*two+r4(8,17)+r4(8,28)*two+r4(10,28)*two+r4(8,29) &
&                  +r4(10,29)-r3(6,21)*four-r3(6,22)*two+r2(5,7)*two+r2(5,8)+r2(5,19)*two &
&                  +r2(5,20))*xx+(+r6(13,10)+r6(13,11)*two-r5(8,14)*two-r5(8,15)*four &
&                  +r4(4,15)+r4(4,16)*two+r4(4,27)+r4(6,27)+r4(4,28)*two+r4(6,28)*two &
&                  +r4(13,39)+r4(13,40)*two-r3(3,20)*two-r3(3,21)*four-r3(8,44)*two &
&                  -r3(8,45)*four+r2(1,6)+r2(1,7)*two+r2(1,18)+r2(1,19)*two+r2(2,38) &
&                  +r2(2,39)*two+r2(2,53)+r2(3,53)+r2(2,54)*two+r2(3,54)*two-r1(3,32)*two &
&                  -r1(3,33)*four+r0(7)+r0(8)*two+r0(22)+r0(23)*two)*xz+(+r5(19,24) &
&                  -r4(13,33)*two+r3(8,30)+r3(8,42)+r3(10,42)-r2(3,32)*two+r1(3,8)+r1(3,20)) &
&                  *xxx+(+r5(13,23)*two+r5(13,24)-r4(8,32)*four-r4(8,33)*two+r3(4,29)*two &
&                  +r3(4,30)+r3(4,41)*two+r3(6,41)*two+r3(4,42)+r3(6,42)-r2(5,31)*four &
&                  -r2(5,32)*two+r1(1,7)*two+r1(1,8)+r1(1,19)*two+r1(1,20))*xxz+rxyz(1)*xxxz
      eri(6,1,3,2)=r211+(+r7(25,4)*two-r6(18,6)*four+r5(12,6)*two+r5(12,12)*two &
&                  +r5(14,12)*two-r4(9,11)*four+r3(5,4)*two+r3(5,10)*two)*qx+rxyz(17)*qz+( &
&                  +r6(25,12)-r5(18,16)*two+r4(12,17)+r4(12,29)+r4(14,29)-r3(9,22)*two &
&                  +r2(6,8)+r2(6,20))*xx+(+r6(18,11)*two-r5(12,15)*four+r4(7,16)*two &
&                  +r4(7,28)*two+r4(9,28)*two-r3(5,21)*four+r2(4,7)*two+r2(4,19)*two)*xz &
&                  +rxyz(2)*xxz
      eri(1,2,3,2)=r220+(+r7(24,3)*two-r6(17,5)*four+r5(11,5)*two+r5(11,11)*two &
&                  +r5(13,11)*two+r5(13,21)*two-r4(8,10)*four-r4(8,30)*four+r3(4,3)*two &
&                  +r3(4,9)*two+r3(4,27)*two+r3(4,39)*two+r3(6,39)*two-r2(5,29)*four &
&                  +r1(1,5)*two+r1(1,17)*two)*qx+rxyz(4)*xx
      eri(2,2,3,2)=r040
      eri(3,2,3,2)=r022+(+r7(32,3)*two-r6(24,5)*four+r5(17,5)*two+r5(17,11)*two &
&                  +r5(19,11)*two+r5(19,21)*two-r4(13,10)*four-r4(13,30)*four+r3(8,3)*two &
&                  +r3(8,9)*two+r3(8,27)*two+r3(8,39)*two+r3(10,39)*two-r2(3,29)*four &
&                  +r1(3,5)*two+r1(3,17)*two)*qz+rxyz(4)*zz
      eri(4,2,3,2)=r130+rxyz(7)*qx
      eri(5,2,3,2)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,3,2)=r031+rxyz(7)*qz
      eri(1,3,3,2)=r202+(+r7(26,3)*two-r6(19,5)*four+r5(13,5)*two+r5(13,11)*two &
&                  +r5(15,11)*two+r5(13,21)*two-r4(10,10)*four-r4(8,30)*four+r3(6,3)*two &
&                  +r3(6,9)*two+r3(4,27)*two+r3(4,39)*two+r3(6,39)*two-r2(5,29)*four &
&                  +r1(1,5)*two+r1(1,17)*two)*qx+(+r7(19,4)*two-r6(13,6)*four+r5(8,6)*two &
&                  +r5(8,12)*two+r5(10,12)*two+r5(19,22)*two-r4(6,11)*four-r4(13,31)*four &
&                  +r3(3,4)*two+r3(3,10)*two+r3(8,28)*two+r3(8,40)*two+r3(10,40)*two &
&                  -r2(3,30)*four+r1(3,6)*two+r1(3,18)*two)*qz+(+r6(26,10)-r5(19,14)*two &
&                  +r4(13,15)+r4(13,27)+r4(15,27)+r4(13,39)-r3(10,20)*two-r3(8,44)*two &
&                  +r2(3,6)+r2(3,18)+r2(2,38)+r2(2,53)+r2(3,53)-r1(3,32)*two+r0(7)+r0(22))*xx &
&                  +rxyz(20)*xz+(+r6(13,12)-r5(8,16)*two+r4(4,17)+r4(4,29)+r4(6,29)+r4(13,41) &
&                  -r3(3,22)*two-r3(8,46)*two+r2(1,8)+r2(1,20)+r2(2,40)+r2(2,55)+r2(3,55) &
&                  -r1(3,34)*two+r0(9)+r0(24))*zz+(+r5(19,23)*two-r4(13,32)*four+r3(8,29)*two &
&                  +r3(8,41)*two+r3(10,41)*two-r2(3,31)*four+r1(3,7)*two+r1(3,19)*two)*xxz+( &
&                  +r5(13,24)*two-r4(8,33)*four+r3(4,30)*two+r3(4,42)*two+r3(6,42)*two &
&                  -r2(5,32)*four+r1(1,8)*two+r1(1,20)*two)*xzz+rxyz(1)*xxzz
      eri(2,3,3,2)=r022+(+r7(32,4)*two-r6(24,6)*four+r5(17,6)*two+r5(17,12)*two &
&                  +r5(19,12)*two+r5(19,22)*two-r4(13,11)*four-r4(13,31)*four+r3(8,4)*two &
&                  +r3(8,10)*two+r3(8,28)*two+r3(8,40)*two+r3(10,40)*two-r2(3,30)*four &
&                  +r1(3,6)*two+r1(3,18)*two)*qz+rxyz(5)*zz
      eri(3,3,3,2)=r004+(+r7(34,3)*two+r7(34,4)*two-r6(26,5)*four-r6(26,6)*four &
&                  +r5(19,5)*two+r5(19,6)*two+r5(19,11)*two+r5(21,11)*two+r5(19,12)*two &
&                  +r5(21,12)*two+r5(19,21)*six+r5(19,22)*six-r4(15,10)*four-r4(15,11)*four &
&                  -r4(13,30)*p12-r4(13,31)*p12+r3(10,3)*two+r3(10,4)*two+r3(10,9)*two &
&                  +r3(10,10)*two+r3(8,27)*six+r3(8,28)*six+r3(8,39)*six+r3(10,39)*six &
&                  +r3(8,40)*six+r3(10,40)*six-r2(3,29)*p12-r2(3,30)*p12+r1(3,5)*six &
&                  +r1(3,6)*six+r1(3,17)*six+r1(3,18)*six)*qz+(+r6(26,10)+r6(26,11)*four &
&                  +r6(26,12)-r5(19,14)*two-r5(19,15)*eight-r5(19,16)*two+r4(13,15) &
&                  +r4(13,16)*four+r4(13,17)+r4(13,27)+r4(15,27)+r4(13,28)*four &
&                  +r4(15,28)*four+r4(13,29)+r4(15,29)+r4(13,39)+r4(13,40)*four+r4(13,41) &
&                  -r3(10,20)*two-r3(10,21)*eight-r3(10,22)*two-r3(8,44)*two-r3(8,45)*eight &
&                  -r3(8,46)*two+r2(3,6)+r2(3,7)*four+r2(3,8)+r2(3,18)+r2(3,19)*four+r2(3,20) &
&                  +r2(2,38)+r2(2,39)*four+r2(2,40)+r2(2,53)+r2(3,53)+r2(2,54)*four &
&                  +r2(3,54)*four+r2(2,55)+r2(3,55)-r1(3,32)*two-r1(3,33)*eight-r1(3,34)*two &
&                  +r0(7)+r0(8)*four+r0(9)+r0(22)+r0(23)*four+r0(24))*zz+(+r5(19,23)*two &
&                  +r5(19,24)*two-r4(13,32)*four-r4(13,33)*four+r3(8,29)*two+r3(8,30)*two &
&                  +r3(8,41)*two+r3(10,41)*two+r3(8,42)*two+r3(10,42)*two-r2(3,31)*four &
&                  -r2(3,32)*four+r1(3,7)*two+r1(3,8)*two+r1(3,19)*two+r1(3,20)*two)*zzz &
&                  +rxyz(1)*zzzz
      eri(4,3,3,2)=r112+rxyz(18)*qx+(+r7(25,4)*two-r6(18,6)*four+r5(12,6)*two &
&                  +r5(12,12)*two+r5(14,12)*two-r4(9,11)*four+r3(5,4)*two+r3(5,10)*two)*qz+( &
&                  +r6(25,11)*two-r5(18,15)*four+r4(12,16)*two+r4(12,28)*two+r4(14,28)*two &
&                  -r3(9,21)*four+r2(6,7)*two+r2(6,19)*two)*xz+(+r6(18,12)-r5(12,16)*two &
&                  +r4(7,17)+r4(7,29)+r4(9,29)-r3(5,22)*two+r2(4,8)+r2(4,20))*zz+rxyz(2)*xzz
      eri(5,3,3,2)=r103+(+r7(34,3)-r6(26,5)*two+r5(19,5)+r5(19,11)+r5(21,11) &
&                  +r5(19,21)*three-r4(15,10)*two-r4(13,30)*six+r3(10,3)+r3(10,9) &
&                  +r3(8,27)*three+r3(8,39)*three+r3(10,39)*three-r2(3,29)*six+r1(3,5)*three &
&                  +r1(3,17)*three)*qx+(+r7(26,3)+r7(26,4)*two-r6(19,5)*two-r6(19,6)*four &
&                  +r5(13,5)+r5(13,6)*two+r5(13,11)+r5(15,11)+r5(13,12)*two+r5(15,12)*two &
&                  +r5(13,21)+r5(13,22)*two-r4(10,10)*two-r4(10,11)*four-r4(8,30)*two &
&                  -r4(8,31)*four+r3(6,3)+r3(6,4)*two+r3(6,9)+r3(6,10)*two+r3(4,27) &
&                  +r3(4,28)*two+r3(4,39)+r3(6,39)+r3(4,40)*two+r3(6,40)*two-r2(5,29)*two &
&                  -r2(5,30)*four+r1(1,5)+r1(1,6)*two+r1(1,17)+r1(1,18)*two)*qz+(+r6(26,10) &
&                  +r6(26,11)*two-r5(19,14)*two-r5(19,15)*four+r4(13,15)+r4(13,16)*two &
&                  +r4(13,27)+r4(15,27)+r4(13,28)*two+r4(15,28)*two+r4(13,39)+r4(13,40)*two &
&                  -r3(10,20)*two-r3(10,21)*four-r3(8,44)*two-r3(8,45)*four+r2(3,6) &
&                  +r2(3,7)*two+r2(3,18)+r2(3,19)*two+r2(2,38)+r2(2,39)*two+r2(2,53)+r2(3,53) &
&                  +r2(2,54)*two+r2(3,54)*two-r1(3,32)*two-r1(3,33)*four+r0(7)+r0(8)*two &
&                  +r0(22)+r0(23)*two)*xz+(+r6(19,11)*two+r6(19,12)-r5(13,15)*four &
&                  -r5(13,16)*two+r4(8,16)*two+r4(8,17)+r4(8,28)*two+r4(10,28)*two+r4(8,29) &
&                  +r4(10,29)-r3(6,21)*four-r3(6,22)*two+r2(5,7)*two+r2(5,8)+r2(5,19)*two &
&                  +r2(5,20))*zz+(+r5(19,23)*two+r5(19,24)-r4(13,32)*four-r4(13,33)*two &
&                  +r3(8,29)*two+r3(8,30)+r3(8,41)*two+r3(10,41)*two+r3(8,42)+r3(10,42) &
&                  -r2(3,31)*four-r2(3,32)*two+r1(3,7)*two+r1(3,8)+r1(3,19)*two+r1(3,20))*xzz &
&                  +(+r5(13,24)-r4(8,33)*two+r3(4,30)+r3(4,42)+r3(6,42)-r2(5,32)*two+r1(1,8) &
&                  +r1(1,20))*zzz+rxyz(1)*xzzz
      eri(6,3,3,2)=r013+(+r7(33,3)+r7(33,4)*two-r6(25,5)*two-r6(25,6)*four+r5(18,5) &
&                  +r5(18,6)*two+r5(18,11)+r5(20,11)+r5(18,12)*two+r5(20,12)*two+r5(18,21) &
&                  +r5(18,22)*two-r4(14,10)*two-r4(14,11)*four-r4(12,30)*two-r4(12,31)*four &
&                  +r3(9,3)+r3(9,4)*two+r3(9,9)+r3(9,10)*two+r3(7,27)+r3(7,28)*two+r3(7,39) &
&                  +r3(9,39)+r3(7,40)*two+r3(9,40)*two-r2(6,29)*two-r2(6,30)*four+r1(2,5) &
&                  +r1(2,6)*two+r1(2,17)+r1(2,18)*two)*qz+(+r6(25,11)*two+r6(25,12) &
&                  -r5(18,15)*four-r5(18,16)*two+r4(12,16)*two+r4(12,17)+r4(12,28)*two &
&                  +r4(14,28)*two+r4(12,29)+r4(14,29)-r3(9,21)*four-r3(9,22)*two+r2(6,7)*two &
&                  +r2(6,8)+r2(6,19)*two+r2(6,20))*zz+rxyz(2)*zzz
      eri(1,4,3,2)=r310+(+r7(18,3)*two+r7(18,4)-r6(12,5)*four-r6(12,6)*two+r5(7,5)*two &
&                  +r5(7,6)+r5(7,11)*two+r5(9,11)*two+r5(7,12)+r5(9,12)+r5(18,21)*two &
&                  +r5(18,22)-r4(5,10)*four-r4(5,11)*two-r4(12,30)*four-r4(12,31)*two &
&                  +r3(2,3)*two+r3(2,4)+r3(2,9)*two+r3(2,10)+r3(7,27)*two+r3(7,28) &
&                  +r3(7,39)*two+r3(9,39)*two+r3(7,40)+r3(9,40)-r2(6,29)*four-r2(6,30)*two &
&                  +r1(2,5)*two+r1(2,6)+r1(2,17)*two+r1(2,18))*qx+(+r6(18,10)+r6(18,11)*two &
&                  -r5(12,14)*two-r5(12,15)*four+r4(7,15)+r4(7,16)*two+r4(7,27)+r4(9,27) &
&                  +r4(7,28)*two+r4(9,28)*two-r3(5,20)*two-r3(5,21)*four+r2(4,6)+r2(4,7)*two &
&                  +r2(4,18)+r2(4,19)*two)*xx+rxyz(3)*xxx
      eri(2,4,3,2)=r130+rxyz(8)*qx
      eri(3,4,3,2)=r112+rxyz(19)*qx+(+r7(25,3)*two-r6(18,5)*four+r5(12,5)*two &
&                  +r5(12,11)*two+r5(14,11)*two-r4(9,10)*four+r3(5,3)*two+r3(5,9)*two)*qz+( &
&                  +r6(25,11)*two-r5(18,15)*four+r4(12,16)*two+r4(12,28)*two+r4(14,28)*two &
&                  -r3(9,21)*four+r2(6,7)*two+r2(6,19)*two)*xz+(+r6(18,10)-r5(12,14)*two &
&                  +r4(7,15)+r4(7,27)+r4(9,27)-r3(5,20)*two+r2(4,6)+r2(4,18))*zz+rxyz(3)*xzz
      eri(4,4,3,2)=r220+(+r7(24,3)+r7(24,4)-r6(17,5)*two-r6(17,6)*two+r5(11,5)+r5(11,6) &
&                  +r5(11,11)+r5(13,11)+r5(11,12)+r5(13,12)+r5(13,21)+r5(13,22)-r4(8,10)*two &
&                  -r4(8,11)*two-r4(8,30)*two-r4(8,31)*two+r3(4,3)+r3(4,4)+r3(4,9)+r3(4,10) &
&                  +r3(4,27)+r3(4,28)+r3(4,39)+r3(6,39)+r3(4,40)+r3(6,40)-r2(5,29)*two &
&                  -r2(5,30)*two+r1(1,5)+r1(1,6)+r1(1,17)+r1(1,18))*qx+rxyz(6)*xx
      eri(5,4,3,2)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(18,10)+r6(18,11) &
&                  -r5(12,14)*two-r5(12,15)*two+r4(7,15)+r4(7,16)+r4(7,27)+r4(9,27)+r4(7,28) &
&                  +r4(9,28)-r3(5,20)*two-r3(5,21)*two+r2(4,6)+r2(4,7)+r2(4,18)+r2(4,19))*xz &
&                  +rxyz(3)*xxz
      eri(6,4,3,2)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,3,2)=r301+(+r7(19,3)*two+r7(19,4)-r6(13,5)*four-r6(13,6)*two+r5(8,5)*two &
&                  +r5(8,6)+r5(8,11)*two+r5(10,11)*two+r5(8,12)+r5(10,12)+r5(19,21)*two &
&                  +r5(19,22)-r4(6,10)*four-r4(6,11)*two-r4(13,30)*four-r4(13,31)*two &
&                  +r3(3,3)*two+r3(3,4)+r3(3,9)*two+r3(3,10)+r3(8,27)*two+r3(8,28) &
&                  +r3(8,39)*two+r3(10,39)*two+r3(8,40)+r3(10,40)-r2(3,29)*four-r2(3,30)*two &
&                  +r1(3,5)*two+r1(3,6)+r1(3,17)*two+r1(3,18))*qx+(+r7(13,4)-r6(8,6)*two &
&                  +r5(4,6)+r5(4,12)+r5(6,12)+r5(13,22)*three-r4(3,11)*two-r4(8,31)*six &
&                  +r3(1,4)+r3(1,10)+r3(4,28)*three+r3(4,40)*three+r3(6,40)*three &
&                  -r2(5,30)*six+r1(1,6)*three+r1(1,18)*three)*qz+(+r6(19,10)+r6(19,11)*two &
&                  -r5(13,14)*two-r5(13,15)*four+r4(8,15)+r4(8,16)*two+r4(8,27)+r4(10,27) &
&                  +r4(8,28)*two+r4(10,28)*two-r3(6,20)*two-r3(6,21)*four+r2(5,6)+r2(5,7)*two &
&                  +r2(5,18)+r2(5,19)*two)*xx+(+r6(13,11)*two+r6(13,12)-r5(8,15)*four &
&                  -r5(8,16)*two+r4(4,16)*two+r4(4,17)+r4(4,28)*two+r4(6,28)*two+r4(4,29) &
&                  +r4(6,29)+r4(13,40)*two+r4(13,41)-r3(3,21)*four-r3(3,22)*two-r3(8,45)*four &
&                  -r3(8,46)*two+r2(1,7)*two+r2(1,8)+r2(1,19)*two+r2(1,20)+r2(2,39)*two &
&                  +r2(2,40)+r2(2,54)*two+r2(3,54)*two+r2(2,55)+r2(3,55)-r1(3,33)*four &
&                  -r1(3,34)*two+r0(8)*two+r0(9)+r0(23)*two+r0(24))*xz+(+r5(19,23) &
&                  -r4(13,32)*two+r3(8,29)+r3(8,41)+r3(10,41)-r2(3,31)*two+r1(3,7)+r1(3,19)) &
&                  *xxx+(+r5(13,23)+r5(13,24)*two-r4(8,32)*two-r4(8,33)*four+r3(4,29) &
&                  +r3(4,30)*two+r3(4,41)+r3(6,41)+r3(4,42)*two+r3(6,42)*two-r2(5,31)*two &
&                  -r2(5,32)*four+r1(1,7)+r1(1,8)*two+r1(1,19)+r1(1,20)*two)*xxz+rxyz(1)*xxxz
      eri(2,5,3,2)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,3,2)=r103+(+r7(34,4)-r6(26,6)*two+r5(19,6)+r5(19,12)+r5(21,12) &
&                  +r5(19,22)*three-r4(15,11)*two-r4(13,31)*six+r3(10,4)+r3(10,10) &
&                  +r3(8,28)*three+r3(8,40)*three+r3(10,40)*three-r2(3,30)*six+r1(3,6)*three &
&                  +r1(3,18)*three)*qx+(+r7(26,3)*two+r7(26,4)-r6(19,5)*four-r6(19,6)*two &
&                  +r5(13,5)*two+r5(13,6)+r5(13,11)*two+r5(15,11)*two+r5(13,12)+r5(15,12) &
&                  +r5(13,21)*two+r5(13,22)-r4(10,10)*four-r4(10,11)*two-r4(8,30)*four &
&                  -r4(8,31)*two+r3(6,3)*two+r3(6,4)+r3(6,9)*two+r3(6,10)+r3(4,27)*two &
&                  +r3(4,28)+r3(4,39)*two+r3(6,39)*two+r3(4,40)+r3(6,40)-r2(5,29)*four &
&                  -r2(5,30)*two+r1(1,5)*two+r1(1,6)+r1(1,17)*two+r1(1,18))*qz+( &
&                  +r6(26,11)*two+r6(26,12)-r5(19,15)*four-r5(19,16)*two+r4(13,16)*two &
&                  +r4(13,17)+r4(13,28)*two+r4(15,28)*two+r4(13,29)+r4(15,29)+r4(13,40)*two &
&                  +r4(13,41)-r3(10,21)*four-r3(10,22)*two-r3(8,45)*four-r3(8,46)*two &
&                  +r2(3,7)*two+r2(3,8)+r2(3,19)*two+r2(3,20)+r2(2,39)*two+r2(2,40) &
&                  +r2(2,54)*two+r2(3,54)*two+r2(2,55)+r2(3,55)-r1(3,33)*four-r1(3,34)*two &
&                  +r0(8)*two+r0(9)+r0(23)*two+r0(24))*xz+(+r6(19,10)+r6(19,11)*two &
&                  -r5(13,14)*two-r5(13,15)*four+r4(8,15)+r4(8,16)*two+r4(8,27)+r4(10,27) &
&                  +r4(8,28)*two+r4(10,28)*two-r3(6,20)*two-r3(6,21)*four+r2(5,6)+r2(5,7)*two &
&                  +r2(5,18)+r2(5,19)*two)*zz+(+r5(19,23)+r5(19,24)*two-r4(13,32)*two &
&                  -r4(13,33)*four+r3(8,29)+r3(8,30)*two+r3(8,41)+r3(10,41)+r3(8,42)*two &
&                  +r3(10,42)*two-r2(3,31)*two-r2(3,32)*four+r1(3,7)+r1(3,8)*two+r1(3,19) &
&                  +r1(3,20)*two)*xzz+(+r5(13,23)-r4(8,32)*two+r3(4,29)+r3(4,41)+r3(6,41) &
&                  -r2(5,31)*two+r1(1,7)+r1(1,19))*zzz+rxyz(1)*xzzz
      eri(4,5,3,2)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(18,11)+r6(18,12) &
&                  -r5(12,15)*two-r5(12,16)*two+r4(7,16)+r4(7,17)+r4(7,28)+r4(9,28)+r4(7,29) &
&                  +r4(9,29)-r3(5,21)*two-r3(5,22)*two+r2(4,7)+r2(4,8)+r2(4,19)+r2(4,20))*xz &
&                  +rxyz(2)*xxz
      eri(5,5,3,2)=r202+(+r7(26,3)+r7(26,4)-r6(19,5)*two-r6(19,6)*two+r5(13,5)+r5(13,6) &
&                  +r5(13,11)+r5(15,11)+r5(13,12)+r5(15,12)+r5(13,21)+r5(13,22)-r4(10,10)*two &
&                  -r4(10,11)*two-r4(8,30)*two-r4(8,31)*two+r3(6,3)+r3(6,4)+r3(6,9)+r3(6,10) &
&                  +r3(4,27)+r3(4,28)+r3(4,39)+r3(6,39)+r3(4,40)+r3(6,40)-r2(5,29)*two &
&                  -r2(5,30)*two+r1(1,5)+r1(1,6)+r1(1,17)+r1(1,18))*qx+(+r7(19,3)+r7(19,4) &
&                  -r6(13,5)*two-r6(13,6)*two+r5(8,5)+r5(8,6)+r5(8,11)+r5(10,11)+r5(8,12) &
&                  +r5(10,12)+r5(19,21)+r5(19,22)-r4(6,10)*two-r4(6,11)*two-r4(13,30)*two &
&                  -r4(13,31)*two+r3(3,3)+r3(3,4)+r3(3,9)+r3(3,10)+r3(8,27)+r3(8,28)+r3(8,39) &
&                  +r3(10,39)+r3(8,40)+r3(10,40)-r2(3,29)*two-r2(3,30)*two+r1(3,5)+r1(3,6) &
&                  +r1(3,17)+r1(3,18))*qz+(+r6(26,11)-r5(19,15)*two+r4(13,16)+r4(13,28) &
&                  +r4(15,28)+r4(13,40)-r3(10,21)*two-r3(8,45)*two+r2(3,7)+r2(3,19)+r2(2,39) &
&                  +r2(2,54)+r2(3,54)-r1(3,33)*two+r0(8)+r0(23))*xx+(+r6(19,10)+r6(19,11)*two &
&                  +r6(19,12)-r5(13,14)*two-r5(13,15)*four-r5(13,16)*two+r4(8,15) &
&                  +r4(8,16)*two+r4(8,17)+r4(8,27)+r4(10,27)+r4(8,28)*two+r4(10,28)*two &
&                  +r4(8,29)+r4(10,29)-r3(6,20)*two-r3(6,21)*four-r3(6,22)*two+r2(5,6) &
&                  +r2(5,7)*two+r2(5,8)+r2(5,18)+r2(5,19)*two+r2(5,20))*xz+(+r6(13,11) &
&                  -r5(8,15)*two+r4(4,16)+r4(4,28)+r4(6,28)+r4(13,40)-r3(3,21)*two &
&                  -r3(8,45)*two+r2(1,7)+r2(1,19)+r2(2,39)+r2(2,54)+r2(3,54)-r1(3,33)*two &
&                  +r0(8)+r0(23))*zz+(+r5(19,23)+r5(19,24)-r4(13,32)*two-r4(13,33)*two &
&                  +r3(8,29)+r3(8,30)+r3(8,41)+r3(10,41)+r3(8,42)+r3(10,42)-r2(3,31)*two &
&                  -r2(3,32)*two+r1(3,7)+r1(3,8)+r1(3,19)+r1(3,20))*xxz+(+r5(13,23)+r5(13,24) &
&                  -r4(8,32)*two-r4(8,33)*two+r3(4,29)+r3(4,30)+r3(4,41)+r3(6,41)+r3(4,42) &
&                  +r3(6,42)-r2(5,31)*two-r2(5,32)*two+r1(1,7)+r1(1,8)+r1(1,19)+r1(1,20))*xzz &
&                  +rxyz(1)*xxzz
      eri(6,5,3,2)=r112+rxyz(19)*qx+(+r7(25,3)+r7(25,4)-r6(18,5)*two-r6(18,6)*two &
&                  +r5(12,5)+r5(12,6)+r5(12,11)+r5(14,11)+r5(12,12)+r5(14,12)-r4(9,10)*two &
&                  -r4(9,11)*two+r3(5,3)+r3(5,4)+r3(5,9)+r3(5,10))*qz+(+r6(25,11)+r6(25,12) &
&                  -r5(18,15)*two-r5(18,16)*two+r4(12,16)+r4(12,17)+r4(12,28)+r4(14,28) &
&                  +r4(12,29)+r4(14,29)-r3(9,21)*two-r3(9,22)*two+r2(6,7)+r2(6,8)+r2(6,19) &
&                  +r2(6,20))*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,3,2)=r211+(+r7(25,3)*two-r6(18,5)*four+r5(12,5)*two+r5(12,11)*two &
&                  +r5(14,11)*two-r4(9,10)*four+r3(5,3)*two+r3(5,9)*two)*qx+rxyz(16)*qz+( &
&                  +r6(25,10)-r5(18,14)*two+r4(12,15)+r4(12,27)+r4(14,27)-r3(9,20)*two &
&                  +r2(6,6)+r2(6,18))*xx+(+r6(18,11)*two-r5(12,15)*four+r4(7,16)*two &
&                  +r4(7,28)*two+r4(9,28)*two-r3(5,21)*four+r2(4,7)*two+r2(4,19)*two)*xz &
&                  +rxyz(3)*xxz
      eri(2,6,3,2)=r031+rxyz(8)*qz
      eri(3,6,3,2)=r013+(+r7(33,3)*two+r7(33,4)-r6(25,5)*four-r6(25,6)*two+r5(18,5)*two &
&                  +r5(18,6)+r5(18,11)*two+r5(20,11)*two+r5(18,12)+r5(20,12)+r5(18,21)*two &
&                  +r5(18,22)-r4(14,10)*four-r4(14,11)*two-r4(12,30)*four-r4(12,31)*two &
&                  +r3(9,3)*two+r3(9,4)+r3(9,9)*two+r3(9,10)+r3(7,27)*two+r3(7,28) &
&                  +r3(7,39)*two+r3(9,39)*two+r3(7,40)+r3(9,40)-r2(6,29)*four-r2(6,30)*two &
&                  +r1(2,5)*two+r1(2,6)+r1(2,17)*two+r1(2,18))*qz+(+r6(25,10)+r6(25,11)*two &
&                  -r5(18,14)*two-r5(18,15)*four+r4(12,15)+r4(12,16)*two+r4(12,27)+r4(14,27) &
&                  +r4(12,28)*two+r4(14,28)*two-r3(9,20)*two-r3(9,21)*four+r2(6,6) &
&                  +r2(6,7)*two+r2(6,18)+r2(6,19)*two)*zz+rxyz(3)*zzz
      eri(4,6,3,2)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,3,2)=r112+rxyz(18)*qx+(+r7(25,3)+r7(25,4)-r6(18,5)*two-r6(18,6)*two &
&                  +r5(12,5)+r5(12,6)+r5(12,11)+r5(14,11)+r5(12,12)+r5(14,12)-r4(9,10)*two &
&                  -r4(9,11)*two+r3(5,3)+r3(5,4)+r3(5,9)+r3(5,10))*qz+(+r6(25,10)+r6(25,11) &
&                  -r5(18,14)*two-r5(18,15)*two+r4(12,15)+r4(12,16)+r4(12,27)+r4(14,27) &
&                  +r4(12,28)+r4(14,28)-r3(9,20)*two-r3(9,21)*two+r2(6,6)+r2(6,7)+r2(6,18) &
&                  +r2(6,19))*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,3,2)=r022+(+r7(32,3)+r7(32,4)-r6(24,5)*two-r6(24,6)*two+r5(17,5)+r5(17,6) &
&                  +r5(17,11)+r5(19,11)+r5(17,12)+r5(19,12)+r5(19,21)+r5(19,22)-r4(13,10)*two &
&                  -r4(13,11)*two-r4(13,30)*two-r4(13,31)*two+r3(8,3)+r3(8,4)+r3(8,9) &
&                  +r3(8,10)+r3(8,27)+r3(8,28)+r3(8,39)+r3(10,39)+r3(8,40)+r3(10,40) &
&                  -r2(3,29)*two-r2(3,30)*two+r1(3,5)+r1(3,6)+r1(3,17)+r1(3,18))*qz+rxyz(6) &
&                  *zz
!
      r400= r8(7)+r6(2,4)*three+r6(7,9)*six+r4(2,26)*p18+r4(7,38)*three+r2(4,52)*nine
      r310= r8(11)+r6(4,4)*three+r6(11,9)*three+r4(4,26)*nine
      r301= r8(12)+r6(5,4)*three+r6(12,9)*three+r4(5,26)*nine
      r220= r8(16)+r6(7,4)*three+r6(7,9)+r6(16,9)+r4(2,26)*three+r4(7,26)*three+r4(7,38) &
&          +r2(4,52)*three
      r211= r8(17)+r6(8,4)*three+r6(17,9)+r4(8,26)*three
      r202= r8(18)+r6(9,4)*three+r6(7,9)+r6(18,9)+r4(2,26)*three+r4(9,26)*three+r4(7,38) &
&          +r2(4,52)*three
      r130= r8(22)+r6(11,4)*three+r6(11,9)*three+r4(4,26)*nine
      r121= r8(23)+r6(12,4)*three+r6(12,9)+r4(5,26)*three
      r112= r8(24)+r6(13,4)*three+r6(11,9)+r4(4,26)*three
      r103= r8(25)+r6(14,4)*three+r6(12,9)*three+r4(5,26)*nine
      r040= r8(29)+r6(16,4)*three+r6(16,9)*six+r4(7,26)*p18+r4(7,38)*three+r2(4,52)*nine
      r031= r8(30)+r6(17,4)*three+r6(17,9)*three+r4(8,26)*nine
      r022= r8(31)+r6(18,4)*three+r6(16,9)+r6(18,9)+r4(7,26)*three+r4(9,26)*three+r4(7,38) &
&          +r2(4,52)*three
      r013= r8(32)+r6(19,4)*three+r6(17,9)*three+r4(8,26)*nine
      r004= r8(33)+r6(20,4)*three+r6(18,9)*six+r4(9,26)*p18+r4(7,38)*three+r2(4,52)*nine
      rxyz(1)=+r4(7,42)+r2(4,56)*three
      rxyz(2)=+r5(11,24)+r3(4,42)*three
      rxyz(3)=+r5(11,23)+r3(4,41)*three
      rxyz(4)=+r6(16,10)+r4(7,27)*three+r4(7,39)+r2(4,53)*three
      rxyz(5)=+r6(16,12)+r4(7,29)*three+r4(7,41)+r2(4,55)*three
      rxyz(6)=+r6(16,11)+r4(7,28)*three+r4(7,40)+r2(4,54)*three
      rxyz(7)=+r7(22,3)+r5(11,11)*three+r5(11,21)*three+r3(4,39)*nine
      rxyz(8)=+r7(22,4)+r5(11,12)*three+r5(11,22)*three+r3(4,40)*nine
      rxyz(9)=+r7(17,3)+r7(17,4)+r5(8,11)*three+r5(8,12)*three
      rxyz(10)=+r6(17,11)+r4(8,28)*three
      rxyz(11)=+r6(11,11)+r4(4,28)*three
      rxyz(12)=+r7(23,3)+r5(12,11)*three+r5(12,21)+r3(5,39)*three
      rxyz(13)=+r7(23,4)+r5(12,12)*three+r5(12,22)+r3(5,40)*three
      rxyz(14)=+r7(16,4)+r5(7,12)*three+r5(7,22)+r3(2,40)*three
      rxyz(15)=+r7(16,3)+r5(7,11)*three+r5(7,21)+r3(2,39)*three
      rxyz(16)=+r7(11,4)+r5(4,12)*three+r5(11,22)+r3(4,40)*three
      rxyz(17)=+r7(11,3)+r5(4,11)*three+r5(11,21)+r3(4,39)*three
      rxyz(18)=+r7(24,3)+r5(13,11)*three+r5(11,21)+r3(4,39)*three
      rxyz(19)=+r7(24,4)+r5(13,12)*three+r5(11,22)+r3(4,40)*three
      rxyz(20)=+r6(12,11)*four+r4(5,28)*p12
      eri(1,1,4,2)=r400+(+r7(7,3)*two+r7(7,4)*two+r5(2,11)*six+r5(2,12)*six+r5(7,21)*six &
&                  +r5(7,22)*six+r3(2,39)*p18+r3(2,40)*p18)*qx+(+r6(7,10)+r6(7,11)*four &
&                  +r6(7,12)+r4(2,27)*three+r4(2,28)*p12+r4(2,29)*three+r4(7,39) &
&                  +r4(7,40)*four+r4(7,41)+r2(4,53)*three+r2(4,54)*p12+r2(4,55)*three)*xx+( &
&                  +r5(7,23)*two+r5(7,24)*two+r3(2,41)*six+r3(2,42)*six)*xxx+rxyz(1)*xxxx
      eri(2,1,4,2)=r220+(+r7(16,4)*two+r5(7,12)*six+r5(7,22)*two+r3(2,40)*six)*qx+rxyz(5) &
&                  *xx
      eri(3,1,4,2)=r202+(+r7(18,4)*two+r5(9,12)*six+r5(7,22)*two+r3(2,40)*six)*qx+( &
&                  +r7(12,3)*two+r5(5,11)*six+r5(12,21)*two+r3(5,39)*six)*qz+(+r6(18,12) &
&                  +r4(9,29)*three+r4(7,41)+r2(4,55)*three)*xx+rxyz(20)*xz+(+r6(7,10) &
&                  +r4(2,27)*three+r4(7,39)+r2(4,53)*three)*zz+(+r5(12,24)*two+r3(5,42)*six) &
&                  *xxz+(+r5(7,23)*two+r3(2,41)*six)*xzz+rxyz(1)*xxzz
      eri(4,1,4,2)=r310+(+r7(11,3)+r7(11,4)*two+r5(4,11)*three+r5(4,12)*six+r5(11,21) &
&                  +r5(11,22)*two+r3(4,39)*three+r3(4,40)*six)*qx+(+r6(11,11)*two+r6(11,12) &
&                  +r4(4,28)*six+r4(4,29)*three)*xx+rxyz(2)*xxx
      eri(5,1,4,2)=r301+(+r7(12,3)+r7(12,4)*two+r5(5,11)*three+r5(5,12)*six+r5(12,21) &
&                  +r5(12,22)*two+r3(5,39)*three+r3(5,40)*six)*qx+(+r7(7,3)+r5(2,11)*three &
&                  +r5(7,21)*three+r3(2,39)*nine)*qz+(+r6(12,11)*two+r6(12,12)+r4(5,28)*six &
&                  +r4(5,29)*three)*xx+(+r6(7,10)+r6(7,11)*two+r4(2,27)*three+r4(2,28)*six &
&                  +r4(7,39)+r4(7,40)*two+r2(4,53)*three+r2(4,54)*six)*xz+(+r5(12,24) &
&                  +r3(5,42)*three)*xxx+(+r5(7,23)*two+r5(7,24)+r3(2,41)*six+r3(2,42)*three) &
&                  *xxz+rxyz(1)*xxxz
      eri(6,1,4,2)=r211+(+r7(17,4)*two+r5(8,12)*six)*qx+rxyz(17)*qz+(+r6(17,12) &
&                  +r4(8,29)*three)*xx+(+r6(11,11)*two+r4(4,28)*six)*xz+rxyz(2)*xxz
      eri(1,2,4,2)=r220+(+r7(16,3)*two+r5(7,11)*six+r5(7,21)*two+r3(2,39)*six)*qx+rxyz(4) &
&                  *xx
      eri(2,2,4,2)=r040
      eri(3,2,4,2)=r022+(+r7(23,3)*two+r5(12,11)*six+r5(12,21)*two+r3(5,39)*six)*qz &
&                  +rxyz(4)*zz
      eri(4,2,4,2)=r130+rxyz(7)*qx
      eri(5,2,4,2)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,4,2)=r031+rxyz(7)*qz
      eri(1,3,4,2)=r202+(+r7(18,3)*two+r5(9,11)*six+r5(7,21)*two+r3(2,39)*six)*qx+( &
&                  +r7(12,4)*two+r5(5,12)*six+r5(12,22)*two+r3(5,40)*six)*qz+(+r6(18,10) &
&                  +r4(9,27)*three+r4(7,39)+r2(4,53)*three)*xx+rxyz(20)*xz+(+r6(7,12) &
&                  +r4(2,29)*three+r4(7,41)+r2(4,55)*three)*zz+(+r5(12,23)*two+r3(5,41)*six) &
&                  *xxz+(+r5(7,24)*two+r3(2,42)*six)*xzz+rxyz(1)*xxzz
      eri(2,3,4,2)=r022+(+r7(23,4)*two+r5(12,12)*six+r5(12,22)*two+r3(5,40)*six)*qz &
&                  +rxyz(5)*zz
      eri(3,3,4,2)=r004+(+r7(25,3)*two+r7(25,4)*two+r5(14,11)*six+r5(14,12)*six &
&                  +r5(12,21)*six+r5(12,22)*six+r3(5,39)*p18+r3(5,40)*p18)*qz+(+r6(18,10) &
&                  +r6(18,11)*four+r6(18,12)+r4(9,27)*three+r4(9,28)*p12+r4(9,29)*three &
&                  +r4(7,39)+r4(7,40)*four+r4(7,41)+r2(4,53)*three+r2(4,54)*p12 &
&                  +r2(4,55)*three)*zz+(+r5(12,23)*two+r5(12,24)*two+r3(5,41)*six &
&                  +r3(5,42)*six)*zzz+rxyz(1)*zzzz
      eri(4,3,4,2)=r112+rxyz(18)*qx+(+r7(17,4)*two+r5(8,12)*six)*qz+(+r6(17,11)*two &
&                  +r4(8,28)*six)*xz+(+r6(11,12)+r4(4,29)*three)*zz+rxyz(2)*xzz
      eri(5,3,4,2)=r103+(+r7(25,3)+r5(14,11)*three+r5(12,21)*three+r3(5,39)*nine)*qx+( &
&                  +r7(18,3)+r7(18,4)*two+r5(9,11)*three+r5(9,12)*six+r5(7,21)+r5(7,22)*two &
&                  +r3(2,39)*three+r3(2,40)*six)*qz+(+r6(18,10)+r6(18,11)*two+r4(9,27)*three &
&                  +r4(9,28)*six+r4(7,39)+r4(7,40)*two+r2(4,53)*three+r2(4,54)*six)*xz+( &
&                  +r6(12,11)*two+r6(12,12)+r4(5,28)*six+r4(5,29)*three)*zz+(+r5(12,23)*two &
&                  +r5(12,24)+r3(5,41)*six+r3(5,42)*three)*xzz+(+r5(7,24)+r3(2,42)*three)*zzz &
&                  +rxyz(1)*xzzz
      eri(6,3,4,2)=r013+(+r7(24,3)+r7(24,4)*two+r5(13,11)*three+r5(13,12)*six+r5(11,21) &
&                  +r5(11,22)*two+r3(4,39)*three+r3(4,40)*six)*qz+(+r6(17,11)*two+r6(17,12) &
&                  +r4(8,28)*six+r4(8,29)*three)*zz+rxyz(2)*zzz
      eri(1,4,4,2)=r310+(+r7(11,3)*two+r7(11,4)+r5(4,11)*six+r5(4,12)*three+r5(11,21)*two &
&                  +r5(11,22)+r3(4,39)*six+r3(4,40)*three)*qx+(+r6(11,10)+r6(11,11)*two &
&                  +r4(4,27)*three+r4(4,28)*six)*xx+rxyz(3)*xxx
      eri(2,4,4,2)=r130+rxyz(8)*qx
      eri(3,4,4,2)=r112+rxyz(19)*qx+(+r7(17,3)*two+r5(8,11)*six)*qz+(+r6(17,11)*two &
&                  +r4(8,28)*six)*xz+(+r6(11,10)+r4(4,27)*three)*zz+rxyz(3)*xzz
      eri(4,4,4,2)=r220+(+r7(16,3)+r7(16,4)+r5(7,11)*three+r5(7,12)*three+r5(7,21) &
&                  +r5(7,22)+r3(2,39)*three+r3(2,40)*three)*qx+rxyz(6)*xx
      eri(5,4,4,2)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(11,10)+r6(11,11) &
&                  +r4(4,27)*three+r4(4,28)*three)*xz+rxyz(3)*xxz
      eri(6,4,4,2)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,4,2)=r301+(+r7(12,3)*two+r7(12,4)+r5(5,11)*six+r5(5,12)*three+r5(12,21)*two &
&                  +r5(12,22)+r3(5,39)*six+r3(5,40)*three)*qx+(+r7(7,4)+r5(2,12)*three &
&                  +r5(7,22)*three+r3(2,40)*nine)*qz+(+r6(12,10)+r6(12,11)*two+r4(5,27)*three &
&                  +r4(5,28)*six)*xx+(+r6(7,11)*two+r6(7,12)+r4(2,28)*six+r4(2,29)*three &
&                  +r4(7,40)*two+r4(7,41)+r2(4,54)*six+r2(4,55)*three)*xz+(+r5(12,23) &
&                  +r3(5,41)*three)*xxx+(+r5(7,23)+r5(7,24)*two+r3(2,41)*three+r3(2,42)*six) &
&                  *xxz+rxyz(1)*xxxz
      eri(2,5,4,2)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,4,2)=r103+(+r7(25,4)+r5(14,12)*three+r5(12,22)*three+r3(5,40)*nine)*qx+( &
&                  +r7(18,3)*two+r7(18,4)+r5(9,11)*six+r5(9,12)*three+r5(7,21)*two+r5(7,22) &
&                  +r3(2,39)*six+r3(2,40)*three)*qz+(+r6(18,11)*two+r6(18,12)+r4(9,28)*six &
&                  +r4(9,29)*three+r4(7,40)*two+r4(7,41)+r2(4,54)*six+r2(4,55)*three)*xz+( &
&                  +r6(12,10)+r6(12,11)*two+r4(5,27)*three+r4(5,28)*six)*zz+(+r5(12,23) &
&                  +r5(12,24)*two+r3(5,41)*three+r3(5,42)*six)*xzz+(+r5(7,23)+r3(2,41)*three) &
&                  *zzz+rxyz(1)*xzzz
      eri(4,5,4,2)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(11,11)+r6(11,12) &
&                  +r4(4,28)*three+r4(4,29)*three)*xz+rxyz(2)*xxz
      eri(5,5,4,2)=r202+(+r7(18,3)+r7(18,4)+r5(9,11)*three+r5(9,12)*three+r5(7,21) &
&                  +r5(7,22)+r3(2,39)*three+r3(2,40)*three)*qx+(+r7(12,3)+r7(12,4) &
&                  +r5(5,11)*three+r5(5,12)*three+r5(12,21)+r5(12,22)+r3(5,39)*three &
&                  +r3(5,40)*three)*qz+(+r6(18,11)+r4(9,28)*three+r4(7,40)+r2(4,54)*three)*xx &
&                  +(+r6(12,10)+r6(12,11)*two+r6(12,12)+r4(5,27)*three+r4(5,28)*six &
&                  +r4(5,29)*three)*xz+(+r6(7,11)+r4(2,28)*three+r4(7,40)+r2(4,54)*three)*zz &
&                  +(+r5(12,23)+r5(12,24)+r3(5,41)*three+r3(5,42)*three)*xxz+(+r5(7,23) &
&                  +r5(7,24)+r3(2,41)*three+r3(2,42)*three)*xzz+rxyz(1)*xxzz
      eri(6,5,4,2)=r112+rxyz(19)*qx+(+r7(17,3)+r7(17,4)+r5(8,11)*three+r5(8,12)*three)*qz &
&                  +(+r6(17,11)+r6(17,12)+r4(8,28)*three+r4(8,29)*three)*xz+rxyz(11)*zz &
&                  +rxyz(2)*xzz
      eri(1,6,4,2)=r211+(+r7(17,3)*two+r5(8,11)*six)*qx+rxyz(16)*qz+(+r6(17,10) &
&                  +r4(8,27)*three)*xx+(+r6(11,11)*two+r4(4,28)*six)*xz+rxyz(3)*xxz
      eri(2,6,4,2)=r031+rxyz(8)*qz
      eri(3,6,4,2)=r013+(+r7(24,3)*two+r7(24,4)+r5(13,11)*six+r5(13,12)*three &
&                  +r5(11,21)*two+r5(11,22)+r3(4,39)*six+r3(4,40)*three)*qz+(+r6(17,10) &
&                  +r6(17,11)*two+r4(8,27)*three+r4(8,28)*six)*zz+rxyz(3)*zzz
      eri(4,6,4,2)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,4,2)=r112+rxyz(18)*qx+(+r7(17,3)+r7(17,4)+r5(8,11)*three+r5(8,12)*three)*qz &
&                  +(+r6(17,10)+r6(17,11)+r4(8,27)*three+r4(8,28)*three)*xz+rxyz(11)*zz &
&                  +rxyz(3)*xzz
      eri(6,6,4,2)=r022+(+r7(23,3)+r7(23,4)+r5(12,11)*three+r5(12,12)*three+r5(12,21) &
&                  +r5(12,22)+r3(5,39)*three+r3(5,40)*three)*qz+rxyz(6)*zz
!
      r400= r8(8)-r7(4,1)+r6(3,4)+r6(8,9)*six-r5(1,3)-r5(4,13)*six+r4(3,26)*six &
&          +r4(8,38)*three-r3(1,19)*six-r3(4,43)*three+r2(5,52)*three-r1(1,31)*three
      r310= r8(12)-r7(7,1)+r6(5,4)+r6(12,9)*three-r5(2,3)-r5(7,13)*three+r4(5,26)*three &
&          -r3(2,19)*three
      r301= r8(13)-r7(8,1)+r6(6,4)+r6(13,9)*three-r5(3,3)-r5(8,13)*three+r4(6,26)*three &
&          -r3(3,19)*three
      r220= r8(17)-r7(11,1)+r6(8,4)+r6(8,9)+r6(17,9)-r5(4,3)-r5(4,13)-r5(11,13)+r4(3,26) &
&          +r4(8,26)+r4(8,38)-r3(1,19)-r3(4,19)-r3(4,43)+r2(5,52)-r1(1,31)
      r211= r8(18)-r7(12,1)+r6(9,4)+r6(18,9)-r5(5,3)-r5(12,13)+r4(9,26)-r3(5,19)
      r202= r8(19)-r7(13,1)+r6(10,4)+r6(8,9)+r6(19,9)-r5(6,3)-r5(4,13)-r5(13,13)+r4(3,26) &
&          +r4(10,26)+r4(8,38)-r3(1,19)-r3(6,19)-r3(4,43)+r2(5,52)-r1(1,31)
      r130= r8(23)-r7(16,1)+r6(12,4)+r6(12,9)*three-r5(7,3)-r5(7,13)*three+r4(5,26)*three &
&          -r3(2,19)*three
      r121= r8(24)-r7(17,1)+r6(13,4)+r6(13,9)-r5(8,3)-r5(8,13)+r4(6,26)-r3(3,19)
      r112= r8(25)-r7(18,1)+r6(14,4)+r6(12,9)-r5(9,3)-r5(7,13)+r4(5,26)-r3(2,19)
      r103= r8(26)-r7(19,1)+r6(15,4)+r6(13,9)*three-r5(10,3)-r5(8,13)*three+r4(6,26)*three &
&          -r3(3,19)*three
      r040= r8(30)-r7(22,1)+r6(17,4)+r6(17,9)*six-r5(11,3)-r5(11,13)*six+r4(8,26)*six &
&          +r4(8,38)*three-r3(4,19)*six-r3(4,43)*three+r2(5,52)*three-r1(1,31)*three
      r031= r8(31)-r7(23,1)+r6(18,4)+r6(18,9)*three-r5(12,3)-r5(12,13)*three+r4(9,26)*three &
&          -r3(5,19)*three
      r022= r8(32)-r7(24,1)+r6(19,4)+r6(17,9)+r6(19,9)-r5(13,3)-r5(11,13)-r5(13,13) &
&          +r4(8,26)+r4(10,26)+r4(8,38)-r3(4,19)-r3(6,19)-r3(4,43)+r2(5,52)-r1(1,31)
      r013= r8(33)-r7(25,1)+r6(20,4)+r6(18,9)*three-r5(14,3)-r5(12,13)*three+r4(9,26)*three &
&          -r3(5,19)*three
      r004= r8(34)-r7(26,1)+r6(21,4)+r6(19,9)*six-r5(15,3)-r5(13,13)*six+r4(10,26)*six &
&          +r4(8,38)*three-r3(6,19)*six-r3(4,43)*three+r2(5,52)*three-r1(1,31)*three
      rxyz(1)=+r4(8,42)-r3(4,47)+r2(5,56)-r1(1,35)
      rxyz(2)=+r5(12,24)-r4(7,33)+r3(5,42)-r2(4,32)
      rxyz(3)=+r5(12,23)-r4(7,32)+r3(5,41)-r2(4,31)
      rxyz(4)=+r6(17,10)-r5(11,14)+r4(8,27)+r4(8,39)-r3(4,20)-r3(4,44)+r2(5,53)-r1(1,32)
      rxyz(5)=+r6(17,12)-r5(11,16)+r4(8,29)+r4(8,41)-r3(4,22)-r3(4,46)+r2(5,55)-r1(1,34)
      rxyz(6)=+r6(17,11)-r5(11,15)+r4(8,28)+r4(8,40)-r3(4,21)-r3(4,45)+r2(5,54)-r1(1,33)
      rxyz(7)=+r7(23,3)-r6(16,5)+r5(12,11)+r5(12,21)*three-r4(7,10)-r4(7,30)*three &
&             +r3(5,39)*three-r2(4,29)*three
      rxyz(8)=+r7(23,4)-r6(16,6)+r5(12,12)+r5(12,22)*three-r4(7,11)-r4(7,31)*three &
&             +r3(5,40)*three-r2(4,30)*three
      rxyz(9)=+r7(18,3)+r7(18,4)-r6(12,5)-r6(12,6)+r5(9,11)+r5(9,12)-r4(5,10)-r4(5,11)
      rxyz(10)=+r6(18,11)-r5(12,15)+r4(9,28)-r3(5,21)
      rxyz(11)=+r6(12,11)-r5(7,15)+r4(5,28)-r3(2,21)
      rxyz(12)=+r7(24,3)-r6(17,5)+r5(13,11)+r5(13,21)-r4(8,10)-r4(8,30)+r3(6,39)-r2(5,29)
      rxyz(13)=+r7(24,4)-r6(17,6)+r5(13,12)+r5(13,22)-r4(8,11)-r4(8,31)+r3(6,40)-r2(5,30)
      rxyz(14)=+r7(17,4)-r6(11,6)+r5(8,12)+r5(8,22)-r4(4,11)-r4(4,31)+r3(3,40)-r2(1,30)
      rxyz(15)=+r7(17,3)-r6(11,5)+r5(8,11)+r5(8,21)-r4(4,10)-r4(4,30)+r3(3,39)-r2(1,29)
      rxyz(16)=+r7(12,4)-r6(7,6)+r5(5,12)+r5(12,22)-r4(2,11)-r4(7,31)+r3(5,40)-r2(4,30)
      rxyz(17)=+r7(12,3)-r6(7,5)+r5(5,11)+r5(12,21)-r4(2,10)-r4(7,30)+r3(5,39)-r2(4,29)
      rxyz(18)=+r7(25,3)-r6(18,5)+r5(14,11)+r5(12,21)-r4(9,10)-r4(7,30)+r3(5,39)-r2(4,29)
      rxyz(19)=+r7(25,4)-r6(18,6)+r5(14,12)+r5(12,22)-r4(9,11)-r4(7,31)+r3(5,40)-r2(4,30)
      rxyz(20)=+r6(13,11)*four-r5(8,15)*four+r4(6,28)*four-r3(3,21)*four
      eri(1,1,5,2)=r400+(+r7(8,3)*two+r7(8,4)*two-r6(4,5)*two-r6(4,6)*two+r5(3,11)*two &
&                  +r5(3,12)*two+r5(8,21)*six+r5(8,22)*six-r4(1,10)*two-r4(1,11)*two &
&                  -r4(4,30)*six-r4(4,31)*six+r3(3,39)*six+r3(3,40)*six-r2(1,29)*six &
&                  -r2(1,30)*six)*qx+(+r6(8,10)+r6(8,11)*four+r6(8,12)-r5(4,14)-r5(4,15)*four &
&                  -r5(4,16)+r4(3,27)+r4(3,28)*four+r4(3,29)+r4(8,39)+r4(8,40)*four+r4(8,41) &
&                  -r3(1,20)-r3(1,21)*four-r3(1,22)-r3(4,44)-r3(4,45)*four-r3(4,46)+r2(5,53) &
&                  +r2(5,54)*four+r2(5,55)-r1(1,32)-r1(1,33)*four-r1(1,34))*xx+(+r5(8,23)*two &
&                  +r5(8,24)*two-r4(4,32)*two-r4(4,33)*two+r3(3,41)*two+r3(3,42)*two &
&                  -r2(1,31)*two-r2(1,32)*two)*xxx+rxyz(1)*xxxx
      eri(2,1,5,2)=r220+(+r7(17,4)*two-r6(11,6)*two+r5(8,12)*two+r5(8,22)*two &
&                  -r4(4,11)*two-r4(4,31)*two+r3(3,40)*two-r2(1,30)*two)*qx+rxyz(5)*xx
      eri(3,1,5,2)=r202+(+r7(19,4)*two-r6(13,6)*two+r5(10,12)*two+r5(8,22)*two &
&                  -r4(6,11)*two-r4(4,31)*two+r3(3,40)*two-r2(1,30)*two)*qx+(+r7(13,3)*two &
&                  -r6(8,5)*two+r5(6,11)*two+r5(13,21)*two-r4(3,10)*two-r4(8,30)*two &
&                  +r3(6,39)*two-r2(5,29)*two)*qz+(+r6(19,12)-r5(13,16)+r4(10,29)+r4(8,41) &
&                  -r3(6,22)-r3(4,46)+r2(5,55)-r1(1,34))*xx+rxyz(20)*xz+(+r6(8,10)-r5(4,14) &
&                  +r4(3,27)+r4(8,39)-r3(1,20)-r3(4,44)+r2(5,53)-r1(1,32))*zz+(+r5(13,24)*two &
&                  -r4(8,33)*two+r3(6,42)*two-r2(5,32)*two)*xxz+(+r5(8,23)*two-r4(4,32)*two &
&                  +r3(3,41)*two-r2(1,31)*two)*xzz+rxyz(1)*xxzz
      eri(4,1,5,2)=r310+(+r7(12,3)+r7(12,4)*two-r6(7,5)-r6(7,6)*two+r5(5,11)+r5(5,12)*two &
&                  +r5(12,21)+r5(12,22)*two-r4(2,10)-r4(2,11)*two-r4(7,30)-r4(7,31)*two &
&                  +r3(5,39)+r3(5,40)*two-r2(4,29)-r2(4,30)*two)*qx+(+r6(12,11)*two+r6(12,12) &
&                  -r5(7,15)*two-r5(7,16)+r4(5,28)*two+r4(5,29)-r3(2,21)*two-r3(2,22))*xx &
&                  +rxyz(2)*xxx
      eri(5,1,5,2)=r301+(+r7(13,3)+r7(13,4)*two-r6(8,5)-r6(8,6)*two+r5(6,11)+r5(6,12)*two &
&                  +r5(13,21)+r5(13,22)*two-r4(3,10)-r4(3,11)*two-r4(8,30)-r4(8,31)*two &
&                  +r3(6,39)+r3(6,40)*two-r2(5,29)-r2(5,30)*two)*qx+(+r7(8,3)-r6(4,5) &
&                  +r5(3,11)+r5(8,21)*three-r4(1,10)-r4(4,30)*three+r3(3,39)*three &
&                  -r2(1,29)*three)*qz+(+r6(13,11)*two+r6(13,12)-r5(8,15)*two-r5(8,16) &
&                  +r4(6,28)*two+r4(6,29)-r3(3,21)*two-r3(3,22))*xx+(+r6(8,10)+r6(8,11)*two &
&                  -r5(4,14)-r5(4,15)*two+r4(3,27)+r4(3,28)*two+r4(8,39)+r4(8,40)*two &
&                  -r3(1,20)-r3(1,21)*two-r3(4,44)-r3(4,45)*two+r2(5,53)+r2(5,54)*two &
&                  -r1(1,32)-r1(1,33)*two)*xz+(+r5(13,24)-r4(8,33)+r3(6,42)-r2(5,32))*xxx+( &
&                  +r5(8,23)*two+r5(8,24)-r4(4,32)*two-r4(4,33)+r3(3,41)*two+r3(3,42) &
&                  -r2(1,31)*two-r2(1,32))*xxz+rxyz(1)*xxxz
      eri(6,1,5,2)=r211+(+r7(18,4)*two-r6(12,6)*two+r5(9,12)*two-r4(5,11)*two)*qx &
&                  +rxyz(17)*qz+(+r6(18,12)-r5(12,16)+r4(9,29)-r3(5,22))*xx+(+r6(12,11)*two &
&                  -r5(7,15)*two+r4(5,28)*two-r3(2,21)*two)*xz+rxyz(2)*xxz
      eri(1,2,5,2)=r220+(+r7(17,3)*two-r6(11,5)*two+r5(8,11)*two+r5(8,21)*two &
&                  -r4(4,10)*two-r4(4,30)*two+r3(3,39)*two-r2(1,29)*two)*qx+rxyz(4)*xx
      eri(2,2,5,2)=r040
      eri(3,2,5,2)=r022+(+r7(24,3)*two-r6(17,5)*two+r5(13,11)*two+r5(13,21)*two &
&                  -r4(8,10)*two-r4(8,30)*two+r3(6,39)*two-r2(5,29)*two)*qz+rxyz(4)*zz
      eri(4,2,5,2)=r130+rxyz(7)*qx
      eri(5,2,5,2)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,5,2)=r031+rxyz(7)*qz
      eri(1,3,5,2)=r202+(+r7(19,3)*two-r6(13,5)*two+r5(10,11)*two+r5(8,21)*two &
&                  -r4(6,10)*two-r4(4,30)*two+r3(3,39)*two-r2(1,29)*two)*qx+(+r7(13,4)*two &
&                  -r6(8,6)*two+r5(6,12)*two+r5(13,22)*two-r4(3,11)*two-r4(8,31)*two &
&                  +r3(6,40)*two-r2(5,30)*two)*qz+(+r6(19,10)-r5(13,14)+r4(10,27)+r4(8,39) &
&                  -r3(6,20)-r3(4,44)+r2(5,53)-r1(1,32))*xx+rxyz(20)*xz+(+r6(8,12)-r5(4,16) &
&                  +r4(3,29)+r4(8,41)-r3(1,22)-r3(4,46)+r2(5,55)-r1(1,34))*zz+(+r5(13,23)*two &
&                  -r4(8,32)*two+r3(6,41)*two-r2(5,31)*two)*xxz+(+r5(8,24)*two-r4(4,33)*two &
&                  +r3(3,42)*two-r2(1,32)*two)*xzz+rxyz(1)*xxzz
      eri(2,3,5,2)=r022+(+r7(24,4)*two-r6(17,6)*two+r5(13,12)*two+r5(13,22)*two &
&                  -r4(8,11)*two-r4(8,31)*two+r3(6,40)*two-r2(5,30)*two)*qz+rxyz(5)*zz
      eri(3,3,5,2)=r004+(+r7(26,3)*two+r7(26,4)*two-r6(19,5)*two-r6(19,6)*two &
&                  +r5(15,11)*two+r5(15,12)*two+r5(13,21)*six+r5(13,22)*six-r4(10,10)*two &
&                  -r4(10,11)*two-r4(8,30)*six-r4(8,31)*six+r3(6,39)*six+r3(6,40)*six &
&                  -r2(5,29)*six-r2(5,30)*six)*qz+(+r6(19,10)+r6(19,11)*four+r6(19,12) &
&                  -r5(13,14)-r5(13,15)*four-r5(13,16)+r4(10,27)+r4(10,28)*four+r4(10,29) &
&                  +r4(8,39)+r4(8,40)*four+r4(8,41)-r3(6,20)-r3(6,21)*four-r3(6,22)-r3(4,44) &
&                  -r3(4,45)*four-r3(4,46)+r2(5,53)+r2(5,54)*four+r2(5,55)-r1(1,32) &
&                  -r1(1,33)*four-r1(1,34))*zz+(+r5(13,23)*two+r5(13,24)*two-r4(8,32)*two &
&                  -r4(8,33)*two+r3(6,41)*two+r3(6,42)*two-r2(5,31)*two-r2(5,32)*two)*zzz &
&                  +rxyz(1)*zzzz
      eri(4,3,5,2)=r112+rxyz(18)*qx+(+r7(18,4)*two-r6(12,6)*two+r5(9,12)*two-r4(5,11)*two &
&                  )*qz+(+r6(18,11)*two-r5(12,15)*two+r4(9,28)*two-r3(5,21)*two)*xz+( &
&                  +r6(12,12)-r5(7,16)+r4(5,29)-r3(2,22))*zz+rxyz(2)*xzz
      eri(5,3,5,2)=r103+(+r7(26,3)-r6(19,5)+r5(15,11)+r5(13,21)*three-r4(10,10) &
&                  -r4(8,30)*three+r3(6,39)*three-r2(5,29)*three)*qx+(+r7(19,3)+r7(19,4)*two &
&                  -r6(13,5)-r6(13,6)*two+r5(10,11)+r5(10,12)*two+r5(8,21)+r5(8,22)*two &
&                  -r4(6,10)-r4(6,11)*two-r4(4,30)-r4(4,31)*two+r3(3,39)+r3(3,40)*two &
&                  -r2(1,29)-r2(1,30)*two)*qz+(+r6(19,10)+r6(19,11)*two-r5(13,14) &
&                  -r5(13,15)*two+r4(10,27)+r4(10,28)*two+r4(8,39)+r4(8,40)*two-r3(6,20) &
&                  -r3(6,21)*two-r3(4,44)-r3(4,45)*two+r2(5,53)+r2(5,54)*two-r1(1,32) &
&                  -r1(1,33)*two)*xz+(+r6(13,11)*two+r6(13,12)-r5(8,15)*two-r5(8,16) &
&                  +r4(6,28)*two+r4(6,29)-r3(3,21)*two-r3(3,22))*zz+(+r5(13,23)*two+r5(13,24) &
&                  -r4(8,32)*two-r4(8,33)+r3(6,41)*two+r3(6,42)-r2(5,31)*two-r2(5,32))*xzz+( &
&                  +r5(8,24)-r4(4,33)+r3(3,42)-r2(1,32))*zzz+rxyz(1)*xzzz
      eri(6,3,5,2)=r013+(+r7(25,3)+r7(25,4)*two-r6(18,5)-r6(18,6)*two+r5(14,11) &
&                  +r5(14,12)*two+r5(12,21)+r5(12,22)*two-r4(9,10)-r4(9,11)*two-r4(7,30) &
&                  -r4(7,31)*two+r3(5,39)+r3(5,40)*two-r2(4,29)-r2(4,30)*two)*qz+( &
&                  +r6(18,11)*two+r6(18,12)-r5(12,15)*two-r5(12,16)+r4(9,28)*two+r4(9,29) &
&                  -r3(5,21)*two-r3(5,22))*zz+rxyz(2)*zzz
      eri(1,4,5,2)=r310+(+r7(12,3)*two+r7(12,4)-r6(7,5)*two-r6(7,6)+r5(5,11)*two+r5(5,12) &
&                  +r5(12,21)*two+r5(12,22)-r4(2,10)*two-r4(2,11)-r4(7,30)*two-r4(7,31) &
&                  +r3(5,39)*two+r3(5,40)-r2(4,29)*two-r2(4,30))*qx+(+r6(12,10)+r6(12,11)*two &
&                  -r5(7,14)-r5(7,15)*two+r4(5,27)+r4(5,28)*two-r3(2,20)-r3(2,21)*two)*xx &
&                  +rxyz(3)*xxx
      eri(2,4,5,2)=r130+rxyz(8)*qx
      eri(3,4,5,2)=r112+rxyz(19)*qx+(+r7(18,3)*two-r6(12,5)*two+r5(9,11)*two-r4(5,10)*two &
&                  )*qz+(+r6(18,11)*two-r5(12,15)*two+r4(9,28)*two-r3(5,21)*two)*xz+( &
&                  +r6(12,10)-r5(7,14)+r4(5,27)-r3(2,20))*zz+rxyz(3)*xzz
      eri(4,4,5,2)=r220+(+r7(17,3)+r7(17,4)-r6(11,5)-r6(11,6)+r5(8,11)+r5(8,12)+r5(8,21) &
&                  +r5(8,22)-r4(4,10)-r4(4,11)-r4(4,30)-r4(4,31)+r3(3,39)+r3(3,40)-r2(1,29) &
&                  -r2(1,30))*qx+rxyz(6)*xx
      eri(5,4,5,2)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(12,10)+r6(12,11)-r5(7,14) &
&                  -r5(7,15)+r4(5,27)+r4(5,28)-r3(2,20)-r3(2,21))*xz+rxyz(3)*xxz
      eri(6,4,5,2)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,5,2)=r301+(+r7(13,3)*two+r7(13,4)-r6(8,5)*two-r6(8,6)+r5(6,11)*two+r5(6,12) &
&                  +r5(13,21)*two+r5(13,22)-r4(3,10)*two-r4(3,11)-r4(8,30)*two-r4(8,31) &
&                  +r3(6,39)*two+r3(6,40)-r2(5,29)*two-r2(5,30))*qx+(+r7(8,4)-r6(4,6) &
&                  +r5(3,12)+r5(8,22)*three-r4(1,11)-r4(4,31)*three+r3(3,40)*three &
&                  -r2(1,30)*three)*qz+(+r6(13,10)+r6(13,11)*two-r5(8,14)-r5(8,15)*two &
&                  +r4(6,27)+r4(6,28)*two-r3(3,20)-r3(3,21)*two)*xx+(+r6(8,11)*two+r6(8,12) &
&                  -r5(4,15)*two-r5(4,16)+r4(3,28)*two+r4(3,29)+r4(8,40)*two+r4(8,41) &
&                  -r3(1,21)*two-r3(1,22)-r3(4,45)*two-r3(4,46)+r2(5,54)*two+r2(5,55) &
&                  -r1(1,33)*two-r1(1,34))*xz+(+r5(13,23)-r4(8,32)+r3(6,41)-r2(5,31))*xxx+( &
&                  +r5(8,23)+r5(8,24)*two-r4(4,32)-r4(4,33)*two+r3(3,41)+r3(3,42)*two &
&                  -r2(1,31)-r2(1,32)*two)*xxz+rxyz(1)*xxxz
      eri(2,5,5,2)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,5,2)=r103+(+r7(26,4)-r6(19,6)+r5(15,12)+r5(13,22)*three-r4(10,11) &
&                  -r4(8,31)*three+r3(6,40)*three-r2(5,30)*three)*qx+(+r7(19,3)*two+r7(19,4) &
&                  -r6(13,5)*two-r6(13,6)+r5(10,11)*two+r5(10,12)+r5(8,21)*two+r5(8,22) &
&                  -r4(6,10)*two-r4(6,11)-r4(4,30)*two-r4(4,31)+r3(3,39)*two+r3(3,40) &
&                  -r2(1,29)*two-r2(1,30))*qz+(+r6(19,11)*two+r6(19,12)-r5(13,15)*two &
&                  -r5(13,16)+r4(10,28)*two+r4(10,29)+r4(8,40)*two+r4(8,41)-r3(6,21)*two &
&                  -r3(6,22)-r3(4,45)*two-r3(4,46)+r2(5,54)*two+r2(5,55)-r1(1,33)*two &
&                  -r1(1,34))*xz+(+r6(13,10)+r6(13,11)*two-r5(8,14)-r5(8,15)*two+r4(6,27) &
&                  +r4(6,28)*two-r3(3,20)-r3(3,21)*two)*zz+(+r5(13,23)+r5(13,24)*two-r4(8,32) &
&                  -r4(8,33)*two+r3(6,41)+r3(6,42)*two-r2(5,31)-r2(5,32)*two)*xzz+(+r5(8,23) &
&                  -r4(4,32)+r3(3,41)-r2(1,31))*zzz+rxyz(1)*xzzz
      eri(4,5,5,2)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(12,11)+r6(12,12)-r5(7,15) &
&                  -r5(7,16)+r4(5,28)+r4(5,29)-r3(2,21)-r3(2,22))*xz+rxyz(2)*xxz
      eri(5,5,5,2)=r202+(+r7(19,3)+r7(19,4)-r6(13,5)-r6(13,6)+r5(10,11)+r5(10,12) &
&                  +r5(8,21)+r5(8,22)-r4(6,10)-r4(6,11)-r4(4,30)-r4(4,31)+r3(3,39)+r3(3,40) &
&                  -r2(1,29)-r2(1,30))*qx+(+r7(13,3)+r7(13,4)-r6(8,5)-r6(8,6)+r5(6,11) &
&                  +r5(6,12)+r5(13,21)+r5(13,22)-r4(3,10)-r4(3,11)-r4(8,30)-r4(8,31)+r3(6,39) &
&                  +r3(6,40)-r2(5,29)-r2(5,30))*qz+(+r6(19,11)-r5(13,15)+r4(10,28)+r4(8,40) &
&                  -r3(6,21)-r3(4,45)+r2(5,54)-r1(1,33))*xx+(+r6(13,10)+r6(13,11)*two &
&                  +r6(13,12)-r5(8,14)-r5(8,15)*two-r5(8,16)+r4(6,27)+r4(6,28)*two+r4(6,29) &
&                  -r3(3,20)-r3(3,21)*two-r3(3,22))*xz+(+r6(8,11)-r5(4,15)+r4(3,28)+r4(8,40) &
&                  -r3(1,21)-r3(4,45)+r2(5,54)-r1(1,33))*zz+(+r5(13,23)+r5(13,24)-r4(8,32) &
&                  -r4(8,33)+r3(6,41)+r3(6,42)-r2(5,31)-r2(5,32))*xxz+(+r5(8,23)+r5(8,24) &
&                  -r4(4,32)-r4(4,33)+r3(3,41)+r3(3,42)-r2(1,31)-r2(1,32))*xzz+rxyz(1)*xxzz
      eri(6,5,5,2)=r112+rxyz(19)*qx+(+r7(18,3)+r7(18,4)-r6(12,5)-r6(12,6)+r5(9,11) &
&                  +r5(9,12)-r4(5,10)-r4(5,11))*qz+(+r6(18,11)+r6(18,12)-r5(12,15)-r5(12,16) &
&                  +r4(9,28)+r4(9,29)-r3(5,21)-r3(5,22))*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,5,2)=r211+(+r7(18,3)*two-r6(12,5)*two+r5(9,11)*two-r4(5,10)*two)*qx &
&                  +rxyz(16)*qz+(+r6(18,10)-r5(12,14)+r4(9,27)-r3(5,20))*xx+(+r6(12,11)*two &
&                  -r5(7,15)*two+r4(5,28)*two-r3(2,21)*two)*xz+rxyz(3)*xxz
      eri(2,6,5,2)=r031+rxyz(8)*qz
      eri(3,6,5,2)=r013+(+r7(25,3)*two+r7(25,4)-r6(18,5)*two-r6(18,6)+r5(14,11)*two &
&                  +r5(14,12)+r5(12,21)*two+r5(12,22)-r4(9,10)*two-r4(9,11)-r4(7,30)*two &
&                  -r4(7,31)+r3(5,39)*two+r3(5,40)-r2(4,29)*two-r2(4,30))*qz+(+r6(18,10) &
&                  +r6(18,11)*two-r5(12,14)-r5(12,15)*two+r4(9,27)+r4(9,28)*two-r3(5,20) &
&                  -r3(5,21)*two)*zz+rxyz(3)*zzz
      eri(4,6,5,2)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,5,2)=r112+rxyz(18)*qx+(+r7(18,3)+r7(18,4)-r6(12,5)-r6(12,6)+r5(9,11) &
&                  +r5(9,12)-r4(5,10)-r4(5,11))*qz+(+r6(18,10)+r6(18,11)-r5(12,14)-r5(12,15) &
&                  +r4(9,27)+r4(9,28)-r3(5,20)-r3(5,21))*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,5,2)=r022+(+r7(24,3)+r7(24,4)-r6(17,5)-r6(17,6)+r5(13,11)+r5(13,12) &
&                  +r5(13,21)+r5(13,22)-r4(8,10)-r4(8,11)-r4(8,30)-r4(8,31)+r3(6,39)+r3(6,40) &
&                  -r2(5,29)-r2(5,30))*qz+rxyz(6)*zz
!
      r400= r8(12)-r7(7,1)+r6(5,4)*three+r6(12,9)*six-r5(2,3)*three-r5(7,13)*six &
&          +r4(5,26)*p18+r4(12,38)*three-r3(2,19)*p18-r3(7,43)*three+r2(6,52)*nine &
&          -r1(2,31)*nine
      r310= r8(17)-r7(11,1)+r6(8,4)*three+r6(17,9)*three-r5(4,3)*three-r5(11,13)*three &
&          +r4(8,26)*nine-r3(4,19)*nine
      r301= r8(18)-r7(12,1)+r6(9,4)*three+r6(18,9)*three-r5(5,3)*three-r5(12,13)*three &
&          +r4(9,26)*nine-r3(5,19)*nine
      r220= r8(23)-r7(16,1)+r6(12,4)*three+r6(12,9)+r6(23,9)-r5(7,3)*three-r5(7,13) &
&          -r5(16,13)+r4(5,26)*three+r4(12,26)*three+r4(12,38)-r3(2,19)*three &
&          -r3(7,19)*three-r3(7,43)+r2(6,52)*three-r1(2,31)*three
      r211= r8(24)-r7(17,1)+r6(13,4)*three+r6(24,9)-r5(8,3)*three-r5(17,13)+r4(13,26)*three &
&          -r3(8,19)*three
      r202= r8(25)-r7(18,1)+r6(14,4)*three+r6(12,9)+r6(25,9)-r5(9,3)*three-r5(7,13) &
&          -r5(18,13)+r4(5,26)*three+r4(14,26)*three+r4(12,38)-r3(2,19)*three &
&          -r3(9,19)*three-r3(7,43)+r2(6,52)*three-r1(2,31)*three
      r130= r8(30)-r7(22,1)+r6(17,4)*three+r6(17,9)*three-r5(11,3)*three-r5(11,13)*three &
&          +r4(8,26)*nine-r3(4,19)*nine
      r121= r8(31)-r7(23,1)+r6(18,4)*three+r6(18,9)-r5(12,3)*three-r5(12,13)+r4(9,26)*three &
&          -r3(5,19)*three
      r112= r8(32)-r7(24,1)+r6(19,4)*three+r6(17,9)-r5(13,3)*three-r5(11,13)+r4(8,26)*three &
&          -r3(4,19)*three
      r103= r8(33)-r7(25,1)+r6(20,4)*three+r6(18,9)*three-r5(14,3)*three-r5(12,13)*three &
&          +r4(9,26)*nine-r3(5,19)*nine
      r040= r8(38)-r7(29,1)+r6(23,4)*three+r6(23,9)*six-r5(16,3)*three-r5(16,13)*six &
&          +r4(12,26)*p18+r4(12,38)*three-r3(7,19)*p18-r3(7,43)*three+r2(6,52)*nine &
&          -r1(2,31)*nine
      r031= r8(39)-r7(30,1)+r6(24,4)*three+r6(24,9)*three-r5(17,3)*three-r5(17,13)*three &
&          +r4(13,26)*nine-r3(8,19)*nine
      r022= r8(40)-r7(31,1)+r6(25,4)*three+r6(23,9)+r6(25,9)-r5(18,3)*three-r5(16,13) &
&          -r5(18,13)+r4(12,26)*three+r4(14,26)*three+r4(12,38)-r3(7,19)*three &
&          -r3(9,19)*three-r3(7,43)+r2(6,52)*three-r1(2,31)*three
      r013= r8(41)-r7(32,1)+r6(26,4)*three+r6(24,9)*three-r5(19,3)*three-r5(17,13)*three &
&          +r4(13,26)*nine-r3(8,19)*nine
      r004= r8(42)-r7(33,1)+r6(27,4)*three+r6(25,9)*six-r5(20,3)*three-r5(18,13)*six &
&          +r4(14,26)*p18+r4(12,38)*three-r3(9,19)*p18-r3(7,43)*three+r2(6,52)*nine &
&          -r1(2,31)*nine
      rxyz(1)=+r4(12,42)-r3(7,47)+r2(6,56)*three-r1(2,35)*three
      rxyz(2)=+r5(17,24)-r4(11,33)+r3(8,42)*three-r2(2,32)*three
      rxyz(3)=+r5(17,23)-r4(11,32)+r3(8,41)*three-r2(2,31)*three
      rxyz(4)=+r6(23,10)-r5(16,14)+r4(12,27)*three+r4(12,39)-r3(7,20)*three-r3(7,44) &
&             +r2(6,53)*three-r1(2,32)*three
      rxyz(5)=+r6(23,12)-r5(16,16)+r4(12,29)*three+r4(12,41)-r3(7,22)*three-r3(7,46) &
&             +r2(6,55)*three-r1(2,34)*three
      rxyz(6)=+r6(23,11)-r5(16,15)+r4(12,28)*three+r4(12,40)-r3(7,21)*three-r3(7,45) &
&             +r2(6,54)*three-r1(2,33)*three
      rxyz(7)=+r7(30,3)-r6(22,5)+r5(17,11)*three+r5(17,21)*three-r4(11,10)*three &
&             -r4(11,30)*three+r3(8,39)*nine-r2(2,29)*nine
      rxyz(8)=+r7(30,4)-r6(22,6)+r5(17,12)*three+r5(17,22)*three-r4(11,11)*three &
&             -r4(11,31)*three+r3(8,40)*nine-r2(2,30)*nine
      rxyz(9)=+r7(24,3)+r7(24,4)-r6(17,5)-r6(17,6)+r5(13,11)*three+r5(13,12)*three &
&             -r4(8,10)*three-r4(8,11)*three
      rxyz(10)=+r6(24,11)-r5(17,15)+r4(13,28)*three-r3(8,21)*three
      rxyz(11)=+r6(17,11)-r5(11,15)+r4(8,28)*three-r3(4,21)*three
      rxyz(12)=+r7(31,3)-r6(23,5)+r5(18,11)*three+r5(18,21)-r4(12,10)*three-r4(12,30) &
&             +r3(9,39)*three-r2(6,29)*three
      rxyz(13)=+r7(31,4)-r6(23,6)+r5(18,12)*three+r5(18,22)-r4(12,11)*three-r4(12,31) &
&             +r3(9,40)*three-r2(6,30)*three
      rxyz(14)=+r7(23,4)-r6(16,6)+r5(12,12)*three+r5(12,22)-r4(7,11)*three-r4(7,31) &
&             +r3(5,40)*three-r2(4,30)*three
      rxyz(15)=+r7(23,3)-r6(16,5)+r5(12,11)*three+r5(12,21)-r4(7,10)*three-r4(7,30) &
&             +r3(5,39)*three-r2(4,29)*three
      rxyz(16)=+r7(17,4)-r6(11,6)+r5(8,12)*three+r5(17,22)-r4(4,11)*three-r4(11,31) &
&             +r3(8,40)*three-r2(2,30)*three
      rxyz(17)=+r7(17,3)-r6(11,5)+r5(8,11)*three+r5(17,21)-r4(4,10)*three-r4(11,30) &
&             +r3(8,39)*three-r2(2,29)*three
      rxyz(18)=+r7(32,3)-r6(24,5)+r5(19,11)*three+r5(17,21)-r4(13,10)*three-r4(11,30) &
&             +r3(8,39)*three-r2(2,29)*three
      rxyz(19)=+r7(32,4)-r6(24,6)+r5(19,12)*three+r5(17,22)-r4(13,11)*three-r4(11,31) &
&             +r3(8,40)*three-r2(2,30)*three
      rxyz(20)=+r6(18,11)*four-r5(12,15)*four+r4(9,28)*p12-r3(5,21)*p12
      eri(1,1,6,2)=r400+(+r7(12,3)*two+r7(12,4)*two-r6(7,5)*two-r6(7,6)*two+r5(5,11)*six &
&                  +r5(5,12)*six+r5(12,21)*six+r5(12,22)*six-r4(2,10)*six-r4(2,11)*six &
&                  -r4(7,30)*six-r4(7,31)*six+r3(5,39)*p18+r3(5,40)*p18-r2(4,29)*p18 &
&                  -r2(4,30)*p18)*qx+(+r6(12,10)+r6(12,11)*four+r6(12,12)-r5(7,14) &
&                  -r5(7,15)*four-r5(7,16)+r4(5,27)*three+r4(5,28)*p12+r4(5,29)*three &
&                  +r4(12,39)+r4(12,40)*four+r4(12,41)-r3(2,20)*three-r3(2,21)*p12 &
&                  -r3(2,22)*three-r3(7,44)-r3(7,45)*four-r3(7,46)+r2(6,53)*three &
&                  +r2(6,54)*p12+r2(6,55)*three-r1(2,32)*three-r1(2,33)*p12-r1(2,34)*three) &
&                  *xx+(+r5(12,23)*two+r5(12,24)*two-r4(7,32)*two-r4(7,33)*two+r3(5,41)*six &
&                  +r3(5,42)*six-r2(4,31)*six-r2(4,32)*six)*xxx+rxyz(1)*xxxx
      eri(2,1,6,2)=r220+(+r7(23,4)*two-r6(16,6)*two+r5(12,12)*six+r5(12,22)*two &
&                  -r4(7,11)*six-r4(7,31)*two+r3(5,40)*six-r2(4,30)*six)*qx+rxyz(5)*xx
      eri(3,1,6,2)=r202+(+r7(25,4)*two-r6(18,6)*two+r5(14,12)*six+r5(12,22)*two &
&                  -r4(9,11)*six-r4(7,31)*two+r3(5,40)*six-r2(4,30)*six)*qx+(+r7(18,3)*two &
&                  -r6(12,5)*two+r5(9,11)*six+r5(18,21)*two-r4(5,10)*six-r4(12,30)*two &
&                  +r3(9,39)*six-r2(6,29)*six)*qz+(+r6(25,12)-r5(18,16)+r4(14,29)*three &
&                  +r4(12,41)-r3(9,22)*three-r3(7,46)+r2(6,55)*three-r1(2,34)*three)*xx &
&                  +rxyz(20)*xz+(+r6(12,10)-r5(7,14)+r4(5,27)*three+r4(12,39)-r3(2,20)*three &
&                  -r3(7,44)+r2(6,53)*three-r1(2,32)*three)*zz+(+r5(18,24)*two-r4(12,33)*two &
&                  +r3(9,42)*six-r2(6,32)*six)*xxz+(+r5(12,23)*two-r4(7,32)*two+r3(5,41)*six &
&                  -r2(4,31)*six)*xzz+rxyz(1)*xxzz
      eri(4,1,6,2)=r310+(+r7(17,3)+r7(17,4)*two-r6(11,5)-r6(11,6)*two+r5(8,11)*three &
&                  +r5(8,12)*six+r5(17,21)+r5(17,22)*two-r4(4,10)*three-r4(4,11)*six &
&                  -r4(11,30)-r4(11,31)*two+r3(8,39)*three+r3(8,40)*six-r2(2,29)*three &
&                  -r2(2,30)*six)*qx+(+r6(17,11)*two+r6(17,12)-r5(11,15)*two-r5(11,16) &
&                  +r4(8,28)*six+r4(8,29)*three-r3(4,21)*six-r3(4,22)*three)*xx+rxyz(2)*xxx
      eri(5,1,6,2)=r301+(+r7(18,3)+r7(18,4)*two-r6(12,5)-r6(12,6)*two+r5(9,11)*three &
&                  +r5(9,12)*six+r5(18,21)+r5(18,22)*two-r4(5,10)*three-r4(5,11)*six &
&                  -r4(12,30)-r4(12,31)*two+r3(9,39)*three+r3(9,40)*six-r2(6,29)*three &
&                  -r2(6,30)*six)*qx+(+r7(12,3)-r6(7,5)+r5(5,11)*three+r5(12,21)*three &
&                  -r4(2,10)*three-r4(7,30)*three+r3(5,39)*nine-r2(4,29)*nine)*qz+( &
&                  +r6(18,11)*two+r6(18,12)-r5(12,15)*two-r5(12,16)+r4(9,28)*six &
&                  +r4(9,29)*three-r3(5,21)*six-r3(5,22)*three)*xx+(+r6(12,10)+r6(12,11)*two &
&                  -r5(7,14)-r5(7,15)*two+r4(5,27)*three+r4(5,28)*six+r4(12,39)+r4(12,40)*two &
&                  -r3(2,20)*three-r3(2,21)*six-r3(7,44)-r3(7,45)*two+r2(6,53)*three &
&                  +r2(6,54)*six-r1(2,32)*three-r1(2,33)*six)*xz+(+r5(18,24)-r4(12,33) &
&                  +r3(9,42)*three-r2(6,32)*three)*xxx+(+r5(12,23)*two+r5(12,24)-r4(7,32)*two &
&                  -r4(7,33)+r3(5,41)*six+r3(5,42)*three-r2(4,31)*six-r2(4,32)*three)*xxz &
&                  +rxyz(1)*xxxz
      eri(6,1,6,2)=r211+(+r7(24,4)*two-r6(17,6)*two+r5(13,12)*six-r4(8,11)*six)*qx &
&                  +rxyz(17)*qz+(+r6(24,12)-r5(17,16)+r4(13,29)*three-r3(8,22)*three)*xx+( &
&                  +r6(17,11)*two-r5(11,15)*two+r4(8,28)*six-r3(4,21)*six)*xz+rxyz(2)*xxz
      eri(1,2,6,2)=r220+(+r7(23,3)*two-r6(16,5)*two+r5(12,11)*six+r5(12,21)*two &
&                  -r4(7,10)*six-r4(7,30)*two+r3(5,39)*six-r2(4,29)*six)*qx+rxyz(4)*xx
      eri(2,2,6,2)=r040
      eri(3,2,6,2)=r022+(+r7(31,3)*two-r6(23,5)*two+r5(18,11)*six+r5(18,21)*two &
&                  -r4(12,10)*six-r4(12,30)*two+r3(9,39)*six-r2(6,29)*six)*qz+rxyz(4)*zz
      eri(4,2,6,2)=r130+rxyz(7)*qx
      eri(5,2,6,2)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,6,2)=r031+rxyz(7)*qz
      eri(1,3,6,2)=r202+(+r7(25,3)*two-r6(18,5)*two+r5(14,11)*six+r5(12,21)*two &
&                  -r4(9,10)*six-r4(7,30)*two+r3(5,39)*six-r2(4,29)*six)*qx+(+r7(18,4)*two &
&                  -r6(12,6)*two+r5(9,12)*six+r5(18,22)*two-r4(5,11)*six-r4(12,31)*two &
&                  +r3(9,40)*six-r2(6,30)*six)*qz+(+r6(25,10)-r5(18,14)+r4(14,27)*three &
&                  +r4(12,39)-r3(9,20)*three-r3(7,44)+r2(6,53)*three-r1(2,32)*three)*xx &
&                  +rxyz(20)*xz+(+r6(12,12)-r5(7,16)+r4(5,29)*three+r4(12,41)-r3(2,22)*three &
&                  -r3(7,46)+r2(6,55)*three-r1(2,34)*three)*zz+(+r5(18,23)*two-r4(12,32)*two &
&                  +r3(9,41)*six-r2(6,31)*six)*xxz+(+r5(12,24)*two-r4(7,33)*two+r3(5,42)*six &
&                  -r2(4,32)*six)*xzz+rxyz(1)*xxzz
      eri(2,3,6,2)=r022+(+r7(31,4)*two-r6(23,6)*two+r5(18,12)*six+r5(18,22)*two &
&                  -r4(12,11)*six-r4(12,31)*two+r3(9,40)*six-r2(6,30)*six)*qz+rxyz(5)*zz
      eri(3,3,6,2)=r004+(+r7(33,3)*two+r7(33,4)*two-r6(25,5)*two-r6(25,6)*two &
&                  +r5(20,11)*six+r5(20,12)*six+r5(18,21)*six+r5(18,22)*six-r4(14,10)*six &
&                  -r4(14,11)*six-r4(12,30)*six-r4(12,31)*six+r3(9,39)*p18+r3(9,40)*p18 &
&                  -r2(6,29)*p18-r2(6,30)*p18)*qz+(+r6(25,10)+r6(25,11)*four+r6(25,12) &
&                  -r5(18,14)-r5(18,15)*four-r5(18,16)+r4(14,27)*three+r4(14,28)*p12 &
&                  +r4(14,29)*three+r4(12,39)+r4(12,40)*four+r4(12,41)-r3(9,20)*three &
&                  -r3(9,21)*p12-r3(9,22)*three-r3(7,44)-r3(7,45)*four-r3(7,46) &
&                  +r2(6,53)*three+r2(6,54)*p12+r2(6,55)*three-r1(2,32)*three-r1(2,33)*p12 &
&                  -r1(2,34)*three)*zz+(+r5(18,23)*two+r5(18,24)*two-r4(12,32)*two &
&                  -r4(12,33)*two+r3(9,41)*six+r3(9,42)*six-r2(6,31)*six-r2(6,32)*six)*zzz &
&                  +rxyz(1)*zzzz
      eri(4,3,6,2)=r112+rxyz(18)*qx+(+r7(24,4)*two-r6(17,6)*two+r5(13,12)*six &
&                  -r4(8,11)*six)*qz+(+r6(24,11)*two-r5(17,15)*two+r4(13,28)*six-r3(8,21)*six &
&                  )*xz+(+r6(17,12)-r5(11,16)+r4(8,29)*three-r3(4,22)*three)*zz+rxyz(2)*xzz
      eri(5,3,6,2)=r103+(+r7(33,3)-r6(25,5)+r5(20,11)*three+r5(18,21)*three &
&                  -r4(14,10)*three-r4(12,30)*three+r3(9,39)*nine-r2(6,29)*nine)*qx+( &
&                  +r7(25,3)+r7(25,4)*two-r6(18,5)-r6(18,6)*two+r5(14,11)*three+r5(14,12)*six &
&                  +r5(12,21)+r5(12,22)*two-r4(9,10)*three-r4(9,11)*six-r4(7,30)-r4(7,31)*two &
&                  +r3(5,39)*three+r3(5,40)*six-r2(4,29)*three-r2(4,30)*six)*qz+(+r6(25,10) &
&                  +r6(25,11)*two-r5(18,14)-r5(18,15)*two+r4(14,27)*three+r4(14,28)*six &
&                  +r4(12,39)+r4(12,40)*two-r3(9,20)*three-r3(9,21)*six-r3(7,44)-r3(7,45)*two &
&                  +r2(6,53)*three+r2(6,54)*six-r1(2,32)*three-r1(2,33)*six)*xz+( &
&                  +r6(18,11)*two+r6(18,12)-r5(12,15)*two-r5(12,16)+r4(9,28)*six &
&                  +r4(9,29)*three-r3(5,21)*six-r3(5,22)*three)*zz+(+r5(18,23)*two+r5(18,24) &
&                  -r4(12,32)*two-r4(12,33)+r3(9,41)*six+r3(9,42)*three-r2(6,31)*six &
&                  -r2(6,32)*three)*xzz+(+r5(12,24)-r4(7,33)+r3(5,42)*three-r2(4,32)*three) &
&                  *zzz+rxyz(1)*xzzz
      eri(6,3,6,2)=r013+(+r7(32,3)+r7(32,4)*two-r6(24,5)-r6(24,6)*two+r5(19,11)*three &
&                  +r5(19,12)*six+r5(17,21)+r5(17,22)*two-r4(13,10)*three-r4(13,11)*six &
&                  -r4(11,30)-r4(11,31)*two+r3(8,39)*three+r3(8,40)*six-r2(2,29)*three &
&                  -r2(2,30)*six)*qz+(+r6(24,11)*two+r6(24,12)-r5(17,15)*two-r5(17,16) &
&                  +r4(13,28)*six+r4(13,29)*three-r3(8,21)*six-r3(8,22)*three)*zz+rxyz(2)*zzz
      eri(1,4,6,2)=r310+(+r7(17,3)*two+r7(17,4)-r6(11,5)*two-r6(11,6)+r5(8,11)*six &
&                  +r5(8,12)*three+r5(17,21)*two+r5(17,22)-r4(4,10)*six-r4(4,11)*three &
&                  -r4(11,30)*two-r4(11,31)+r3(8,39)*six+r3(8,40)*three-r2(2,29)*six &
&                  -r2(2,30)*three)*qx+(+r6(17,10)+r6(17,11)*two-r5(11,14)-r5(11,15)*two &
&                  +r4(8,27)*three+r4(8,28)*six-r3(4,20)*three-r3(4,21)*six)*xx+rxyz(3)*xxx
      eri(2,4,6,2)=r130+rxyz(8)*qx
      eri(3,4,6,2)=r112+rxyz(19)*qx+(+r7(24,3)*two-r6(17,5)*two+r5(13,11)*six &
&                  -r4(8,10)*six)*qz+(+r6(24,11)*two-r5(17,15)*two+r4(13,28)*six-r3(8,21)*six &
&                  )*xz+(+r6(17,10)-r5(11,14)+r4(8,27)*three-r3(4,20)*three)*zz+rxyz(3)*xzz
      eri(4,4,6,2)=r220+(+r7(23,3)+r7(23,4)-r6(16,5)-r6(16,6)+r5(12,11)*three &
&                  +r5(12,12)*three+r5(12,21)+r5(12,22)-r4(7,10)*three-r4(7,11)*three &
&                  -r4(7,30)-r4(7,31)+r3(5,39)*three+r3(5,40)*three-r2(4,29)*three &
&                  -r2(4,30)*three)*qx+rxyz(6)*xx
      eri(5,4,6,2)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(17,10)+r6(17,11) &
&                  -r5(11,14)-r5(11,15)+r4(8,27)*three+r4(8,28)*three-r3(4,20)*three &
&                  -r3(4,21)*three)*xz+rxyz(3)*xxz
      eri(6,4,6,2)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,6,2)=r301+(+r7(18,3)*two+r7(18,4)-r6(12,5)*two-r6(12,6)+r5(9,11)*six &
&                  +r5(9,12)*three+r5(18,21)*two+r5(18,22)-r4(5,10)*six-r4(5,11)*three &
&                  -r4(12,30)*two-r4(12,31)+r3(9,39)*six+r3(9,40)*three-r2(6,29)*six &
&                  -r2(6,30)*three)*qx+(+r7(12,4)-r6(7,6)+r5(5,12)*three+r5(12,22)*three &
&                  -r4(2,11)*three-r4(7,31)*three+r3(5,40)*nine-r2(4,30)*nine)*qz+(+r6(18,10) &
&                  +r6(18,11)*two-r5(12,14)-r5(12,15)*two+r4(9,27)*three+r4(9,28)*six &
&                  -r3(5,20)*three-r3(5,21)*six)*xx+(+r6(12,11)*two+r6(12,12)-r5(7,15)*two &
&                  -r5(7,16)+r4(5,28)*six+r4(5,29)*three+r4(12,40)*two+r4(12,41)-r3(2,21)*six &
&                  -r3(2,22)*three-r3(7,45)*two-r3(7,46)+r2(6,54)*six+r2(6,55)*three &
&                  -r1(2,33)*six-r1(2,34)*three)*xz+(+r5(18,23)-r4(12,32)+r3(9,41)*three &
&                  -r2(6,31)*three)*xxx+(+r5(12,23)+r5(12,24)*two-r4(7,32)-r4(7,33)*two &
&                  +r3(5,41)*three+r3(5,42)*six-r2(4,31)*three-r2(4,32)*six)*xxz+rxyz(1)*xxxz
      eri(2,5,6,2)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,6,2)=r103+(+r7(33,4)-r6(25,6)+r5(20,12)*three+r5(18,22)*three &
&                  -r4(14,11)*three-r4(12,31)*three+r3(9,40)*nine-r2(6,30)*nine)*qx+( &
&                  +r7(25,3)*two+r7(25,4)-r6(18,5)*two-r6(18,6)+r5(14,11)*six+r5(14,12)*three &
&                  +r5(12,21)*two+r5(12,22)-r4(9,10)*six-r4(9,11)*three-r4(7,30)*two-r4(7,31) &
&                  +r3(5,39)*six+r3(5,40)*three-r2(4,29)*six-r2(4,30)*three)*qz+( &
&                  +r6(25,11)*two+r6(25,12)-r5(18,15)*two-r5(18,16)+r4(14,28)*six &
&                  +r4(14,29)*three+r4(12,40)*two+r4(12,41)-r3(9,21)*six-r3(9,22)*three &
&                  -r3(7,45)*two-r3(7,46)+r2(6,54)*six+r2(6,55)*three-r1(2,33)*six &
&                  -r1(2,34)*three)*xz+(+r6(18,10)+r6(18,11)*two-r5(12,14)-r5(12,15)*two &
&                  +r4(9,27)*three+r4(9,28)*six-r3(5,20)*three-r3(5,21)*six)*zz+(+r5(18,23) &
&                  +r5(18,24)*two-r4(12,32)-r4(12,33)*two+r3(9,41)*three+r3(9,42)*six &
&                  -r2(6,31)*three-r2(6,32)*six)*xzz+(+r5(12,23)-r4(7,32)+r3(5,41)*three &
&                  -r2(4,31)*three)*zzz+rxyz(1)*xzzz
      eri(4,5,6,2)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(17,11)+r6(17,12) &
&                  -r5(11,15)-r5(11,16)+r4(8,28)*three+r4(8,29)*three-r3(4,21)*three &
&                  -r3(4,22)*three)*xz+rxyz(2)*xxz
      eri(5,5,6,2)=r202+(+r7(25,3)+r7(25,4)-r6(18,5)-r6(18,6)+r5(14,11)*three &
&                  +r5(14,12)*three+r5(12,21)+r5(12,22)-r4(9,10)*three-r4(9,11)*three &
&                  -r4(7,30)-r4(7,31)+r3(5,39)*three+r3(5,40)*three-r2(4,29)*three &
&                  -r2(4,30)*three)*qx+(+r7(18,3)+r7(18,4)-r6(12,5)-r6(12,6)+r5(9,11)*three &
&                  +r5(9,12)*three+r5(18,21)+r5(18,22)-r4(5,10)*three-r4(5,11)*three &
&                  -r4(12,30)-r4(12,31)+r3(9,39)*three+r3(9,40)*three-r2(6,29)*three &
&                  -r2(6,30)*three)*qz+(+r6(25,11)-r5(18,15)+r4(14,28)*three+r4(12,40) &
&                  -r3(9,21)*three-r3(7,45)+r2(6,54)*three-r1(2,33)*three)*xx+(+r6(18,10) &
&                  +r6(18,11)*two+r6(18,12)-r5(12,14)-r5(12,15)*two-r5(12,16)+r4(9,27)*three &
&                  +r4(9,28)*six+r4(9,29)*three-r3(5,20)*three-r3(5,21)*six-r3(5,22)*three) &
&                  *xz+(+r6(12,11)-r5(7,15)+r4(5,28)*three+r4(12,40)-r3(2,21)*three-r3(7,45) &
&                  +r2(6,54)*three-r1(2,33)*three)*zz+(+r5(18,23)+r5(18,24)-r4(12,32) &
&                  -r4(12,33)+r3(9,41)*three+r3(9,42)*three-r2(6,31)*three-r2(6,32)*three) &
&                  *xxz+(+r5(12,23)+r5(12,24)-r4(7,32)-r4(7,33)+r3(5,41)*three+r3(5,42)*three &
&                  -r2(4,31)*three-r2(4,32)*three)*xzz+rxyz(1)*xxzz
      eri(6,5,6,2)=r112+rxyz(19)*qx+(+r7(24,3)+r7(24,4)-r6(17,5)-r6(17,6)+r5(13,11)*three &
&                  +r5(13,12)*three-r4(8,10)*three-r4(8,11)*three)*qz+(+r6(24,11)+r6(24,12) &
&                  -r5(17,15)-r5(17,16)+r4(13,28)*three+r4(13,29)*three-r3(8,21)*three &
&                  -r3(8,22)*three)*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,6,2)=r211+(+r7(24,3)*two-r6(17,5)*two+r5(13,11)*six-r4(8,10)*six)*qx &
&                  +rxyz(16)*qz+(+r6(24,10)-r5(17,14)+r4(13,27)*three-r3(8,20)*three)*xx+( &
&                  +r6(17,11)*two-r5(11,15)*two+r4(8,28)*six-r3(4,21)*six)*xz+rxyz(3)*xxz
      eri(2,6,6,2)=r031+rxyz(8)*qz
      eri(3,6,6,2)=r013+(+r7(32,3)*two+r7(32,4)-r6(24,5)*two-r6(24,6)+r5(19,11)*six &
&                  +r5(19,12)*three+r5(17,21)*two+r5(17,22)-r4(13,10)*six-r4(13,11)*three &
&                  -r4(11,30)*two-r4(11,31)+r3(8,39)*six+r3(8,40)*three-r2(2,29)*six &
&                  -r2(2,30)*three)*qz+(+r6(24,10)+r6(24,11)*two-r5(17,14)-r5(17,15)*two &
&                  +r4(13,27)*three+r4(13,28)*six-r3(8,20)*three-r3(8,21)*six)*zz+rxyz(3)*zzz
      eri(4,6,6,2)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,6,2)=r112+rxyz(18)*qx+(+r7(24,3)+r7(24,4)-r6(17,5)-r6(17,6)+r5(13,11)*three &
&                  +r5(13,12)*three-r4(8,10)*three-r4(8,11)*three)*qz+(+r6(24,10)+r6(24,11) &
&                  -r5(17,14)-r5(17,15)+r4(13,27)*three+r4(13,28)*three-r3(8,20)*three &
&                  -r3(8,21)*three)*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,6,2)=r022+(+r7(31,3)+r7(31,4)-r6(23,5)-r6(23,6)+r5(18,11)*three &
&                  +r5(18,12)*three+r5(18,21)+r5(18,22)-r4(12,10)*three-r4(12,11)*three &
&                  -r4(12,30)-r4(12,31)+r3(9,39)*three+r3(9,40)*three-r2(6,29)*three &
&                  -r2(6,30)*three)*qz+rxyz(6)*zz
      return
end

!-------------------------------------------------------------
  subroutine int2dddd3(eri,r0,r1,r2,r3,r4,r5,r6,r7,r8,qx,qz)
!-------------------------------------------------------------
!
      implicit none
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00, four=4.0D+00, five=5.0D+00
      real(8),parameter :: six=6.0D+00, seven=7.0D+00, eight=8.0D+00, nine=9.0D+00, ten=1.0D+01
      real(8),parameter :: p11=1.1D+01, p12=1.2D+01, p13=1.3D+01, p15=1.5D+01, p16=1.6D+01
      real(8),parameter :: p18=1.8D+01, p20=2.0D+01, p21=2.1D+01, p24=2.4D+01, p28=2.8D+01
      real(8),parameter :: p30=3.0D+01, p36=3.6D+01, p45=4.5D+01, p105=1.05D+2, p210=2.1D+02
      real(8),parameter :: p420=4.2D+02
      real(8),intent(in) :: r0(25), r1(3,40), r2(6,56), r3(10,52), r4(15,42), r5(21,24)
      real(8),intent(in) :: r6(28,12), r7(36,4), r8(45), qx, qz
      real(8),intent(inout) :: eri(6,6,6,6)
      real(8) :: r400, r310, r301, r220, r211, r202, r130, r121, r112, r103, r040, r031
      real(8) :: r022, r013, r004, rxyz(20), xx, xz, zz, xxx, xxz, xzz, zzz
      real(8) :: xxxx, xxxz, xxzz, xzzz, zzzz
!
      xx= qx*qx
      xz= qx*qz
      zz= qz*qz
      xxx= qx*qx*qx
      xxz= qx*qx*qz
      xzz= qx*qz*qz
      zzz= qz*qz*qz
      xxxx= xx*xx
      xxxz= xx*xz
      xxzz= xx*zz
      xzzz= xz*zz
      zzzz= zz*zz
!
      r400= r8(6)-r7(3,2)*two+r6(1,3)+r6(1,4)+r6(6,4)+r6(6,9)*six-r5(3,4)*two-r5(3,17)*p12 &
&          +r4(1,4)+r4(1,5)+r4(1,22)*six+r4(1,26)*six+r4(6,26)*six+r4(6,38)*three &
&          -r3(3,23)*p12-r3(3,48)*six+r2(1,13)*six+r2(1,17)*six+r2(1,47)*three &
&          +r2(1,52)*three+r2(3,52)*three-r1(3,36)*six+r0(16)*three+r0(21)*three
      r310= r8(9)-r7(5,2)*two+r6(2,3)+r6(2,4)+r6(9,4)+r6(9,9)*three-r5(5,4)*two &
&          -r5(5,17)*six+r4(2,4)+r4(2,5)+r4(2,22)*three+r4(2,26)*three+r4(9,26)*three &
&          -r3(5,23)*six+r2(4,13)*three+r2(4,17)*three
      r301= r8(10)-r7(6,2)*two+r6(3,3)+r6(3,4)+r6(10,4)+r6(10,9)*three-r5(6,4)*two &
&          -r5(6,17)*six+r4(3,4)+r4(3,5)+r4(3,22)*three+r4(3,26)*three+r4(10,26)*three &
&          -r3(6,23)*six+r2(5,13)*three+r2(5,17)*three
      r220= r8(13)-r7(8,2)*two+r6(4,3)+r6(4,4)+r6(13,4)+r6(6,9)+r6(13,9)-r5(8,4)*two &
&          -r5(3,17)*two-r5(8,17)*two+r4(4,4)+r4(4,5)+r4(1,22)+r4(4,22)+r4(1,26)+r4(4,26) &
&          +r4(6,26)+r4(13,26)+r4(6,38)-r3(3,23)*two-r3(8,23)*two-r3(3,48)*two+r2(1,13) &
&          +r2(2,13)+r2(1,17)+r2(2,17)+r2(1,47)+r2(1,52)+r2(3,52)-r1(3,36)*two+r0(16) &
&          +r0(21)
      r211= r8(14)-r7(9,2)*two+r6(5,3)+r6(5,4)+r6(14,4)+r6(14,9)-r5(9,4)*two-r5(9,17)*two &
&          +r4(5,4)+r4(5,5)+r4(5,22)+r4(5,26)+r4(14,26)-r3(9,23)*two+r2(6,13)+r2(6,17)
      r202= r8(15)-r7(10,2)*two+r6(6,3)+r6(6,4)+r6(15,4)+r6(6,9)+r6(15,9)-r5(10,4)*two &
&          -r5(3,17)*two-r5(10,17)*two+r4(6,4)+r4(6,5)+r4(1,22)+r4(6,22)+r4(1,26) &
&          +r4(6,26)*two+r4(15,26)+r4(6,38)-r3(3,23)*two-r3(10,23)*two-r3(3,48)*two &
&          +r2(1,13)+r2(3,13)+r2(1,17)+r2(3,17)+r2(1,47)+r2(1,52)+r2(3,52)-r1(3,36)*two &
&          +r0(16)+r0(21)
      r130= r8(18)-r7(12,2)*two+r6(7,3)+r6(7,4)+r6(18,4)+r6(9,9)*three-r5(12,4)*two &
&          -r5(5,17)*six+r4(7,4)+r4(7,5)+r4(2,22)*three+r4(2,26)*three+r4(9,26)*three &
&          -r3(5,23)*six+r2(4,13)*three+r2(4,17)*three
      r121= r8(19)-r7(13,2)*two+r6(8,3)+r6(8,4)+r6(19,4)+r6(10,9)-r5(13,4)*two-r5(6,17)*two &
&          +r4(8,4)+r4(8,5)+r4(3,22)+r4(3,26)+r4(10,26)-r3(6,23)*two+r2(5,13)+r2(5,17)
      r112= r8(20)-r7(14,2)*two+r6(9,3)+r6(9,4)+r6(20,4)+r6(9,9)-r5(14,4)*two-r5(5,17)*two &
&          +r4(9,4)+r4(9,5)+r4(2,22)+r4(2,26)+r4(9,26)-r3(5,23)*two+r2(4,13)+r2(4,17)
      r103= r8(21)-r7(15,2)*two+r6(10,3)+r6(10,4)+r6(21,4)+r6(10,9)*three-r5(15,4)*two &
&          -r5(6,17)*six+r4(10,4)+r4(10,5)+r4(3,22)*three+r4(3,26)*three+r4(10,26)*three &
&          -r3(6,23)*six+r2(5,13)*three+r2(5,17)*three
      r040= r8(24)-r7(17,2)*two+r6(11,3)+r6(11,4)+r6(24,4)+r6(13,9)*six-r5(17,4)*two &
&          -r5(8,17)*p12+r4(11,4)+r4(11,5)+r4(4,22)*six+r4(4,26)*six+r4(13,26)*six &
&          +r4(6,38)*three-r3(8,23)*p12-r3(3,48)*six+r2(2,13)*six+r2(2,17)*six &
&          +r2(1,47)*three+r2(1,52)*three+r2(3,52)*three-r1(3,36)*six+r0(16)*three &
&          +r0(21)*three
      r031= r8(25)-r7(18,2)*two+r6(12,3)+r6(12,4)+r6(25,4)+r6(14,9)*three-r5(18,4)*two &
&          -r5(9,17)*six+r4(12,4)+r4(12,5)+r4(5,22)*three+r4(5,26)*three+r4(14,26)*three &
&          -r3(9,23)*six+r2(6,13)*three+r2(6,17)*three
      r022= r8(26)-r7(19,2)*two+r6(13,3)+r6(13,4)+r6(26,4)+r6(13,9)+r6(15,9)-r5(19,4)*two &
&          -r5(8,17)*two-r5(10,17)*two+r4(13,4)+r4(13,5)+r4(4,22)+r4(6,22)+r4(4,26) &
&          +r4(6,26)+r4(13,26)+r4(15,26)+r4(6,38)-r3(8,23)*two-r3(10,23)*two-r3(3,48)*two &
&          +r2(2,13)+r2(3,13)+r2(2,17)+r2(3,17)+r2(1,47)+r2(1,52)+r2(3,52)-r1(3,36)*two &
&          +r0(16)+r0(21)
      r013= r8(27)-r7(20,2)*two+r6(14,3)+r6(14,4)+r6(27,4)+r6(14,9)*three-r5(20,4)*two &
&          -r5(9,17)*six+r4(14,4)+r4(14,5)+r4(5,22)*three+r4(5,26)*three+r4(14,26)*three &
&          -r3(9,23)*six+r2(6,13)*three+r2(6,17)*three
      r004= r8(28)-r7(21,2)*two+r6(15,3)+r6(15,4)+r6(28,4)+r6(15,9)*six-r5(21,4)*two &
&          -r5(10,17)*p12+r4(15,4)+r4(15,5)+r4(6,22)*six+r4(6,26)*six+r4(15,26)*six &
&          +r4(6,38)*three-r3(10,23)*p12-r3(3,48)*six+r2(3,13)*six+r2(3,17)*six &
&          +r2(1,47)*three+r2(1,52)*three+r2(3,52)*three-r1(3,36)*six+r0(16)*three &
&          +r0(21)*three
      rxyz(1)=+r4(6,42)-r3(3,52)*two+r2(1,51)+r2(1,56)+r2(3,56)-r1(3,40)*two+r0(20)+r0(25)
      rxyz(2)=+r5(9,24)-r4(5,37)*two+r3(2,38)+r3(2,42)+r3(9,42)-r2(6,36)*two+r1(2,16) &
&             +r1(2,20)
      rxyz(3)=+r5(9,23)-r4(5,36)*two+r3(2,37)+r3(2,41)+r3(9,41)-r2(6,35)*two+r1(2,15) &
&             +r1(2,19)
      rxyz(4)=+r6(13,10)-r5(8,18)*two+r4(4,23)+r4(4,27)+r4(13,27)+r4(6,39)-r3(8,24)*two &
&             -r3(3,49)*two+r2(2,14)+r2(2,18)+r2(1,48)+r2(1,53)+r2(3,53)-r1(3,37)*two &
&             +r0(17)+r0(22)
      rxyz(5)=+r6(13,12)-r5(8,20)*two+r4(4,25)+r4(4,29)+r4(13,29)+r4(6,41)-r3(8,26)*two &
&             -r3(3,51)*two+r2(2,16)+r2(2,20)+r2(1,50)+r2(1,55)+r2(3,55)-r1(3,39)*two &
&             +r0(19)+r0(24)
      rxyz(6)=+r6(13,11)-r5(8,19)*two+r4(4,24)+r4(4,28)+r4(13,28)+r4(6,40)-r3(8,25)*two &
&             -r3(3,50)*two+r2(2,15)+r2(2,19)+r2(1,49)+r2(1,54)+r2(3,54)-r1(3,38)*two &
&             +r0(18)+r0(23)
      rxyz(7)=+r7(18,3)-r6(12,7)*two+r5(7,9)+r5(7,11)+r5(18,11)+r5(9,21)*three &
&             -r4(12,12)*two-r4(5,34)*six+r3(7,7)+r3(7,9)+r3(2,35)*three+r3(2,39)*three &
&             +r3(9,39)*three-r2(6,33)*six+r1(2,13)*three+r1(2,17)*three
      rxyz(8)=+r7(18,4)-r6(12,8)*two+r5(7,10)+r5(7,12)+r5(18,12)+r5(9,22)*three &
&             -r4(12,13)*two-r4(5,35)*six+r3(7,8)+r3(7,10)+r3(2,36)*three+r3(2,40)*three &
&             +r3(9,40)*three-r2(6,34)*six+r1(2,14)*three+r1(2,18)*three
      rxyz(9)=+r7(14,3)+r7(14,4)-r6(9,7)*two-r6(9,8)*two+r5(5,9)+r5(5,10)+r5(5,11) &
&             +r5(14,11)+r5(5,12)+r5(14,12)-r4(9,12)*two-r4(9,13)*two+r3(5,7)+r3(5,8) &
&             +r3(5,9)+r3(5,10)
      rxyz(10)=+r6(14,11)-r5(9,19)*two+r4(5,24)+r4(5,28)+r4(14,28)-r3(9,25)*two+r2(6,15) &
&             +r2(6,19)
      rxyz(11)=+r6(9,11)-r5(5,19)*two+r4(2,24)+r4(2,28)+r4(9,28)-r3(5,25)*two+r2(4,15) &
&             +r2(4,19)
      rxyz(12)=+r7(19,3)-r6(13,7)*two+r5(8,9)+r5(8,11)+r5(19,11)+r5(10,21)-r4(13,12)*two &
&             -r4(6,34)*two+r3(8,7)+r3(8,9)+r3(3,35)+r3(3,39)+r3(10,39)-r2(3,33)*two &
&             +r1(3,13)+r1(3,17)
      rxyz(13)=+r7(19,4)-r6(13,8)*two+r5(8,10)+r5(8,12)+r5(19,12)+r5(10,22)-r4(13,13)*two &
&             -r4(6,35)*two+r3(8,8)+r3(8,10)+r3(3,36)+r3(3,40)+r3(10,40)-r2(3,34)*two &
&             +r1(3,14)+r1(3,18)
      rxyz(14)=+r7(13,4)-r6(8,8)*two+r5(4,10)+r5(4,12)+r5(13,12)+r5(6,22)-r4(8,13)*two &
&             -r4(3,35)*two+r3(4,8)+r3(4,10)+r3(1,36)+r3(1,40)+r3(6,40)-r2(5,34)*two &
&             +r1(1,14)+r1(1,18)
      rxyz(15)=+r7(13,3)-r6(8,7)*two+r5(4,9)+r5(4,11)+r5(13,11)+r5(6,21)-r4(8,12)*two &
&             -r4(3,34)*two+r3(4,7)+r3(4,9)+r3(1,35)+r3(1,39)+r3(6,39)-r2(5,33)*two &
&             +r1(1,13)+r1(1,17)
      rxyz(16)=+r7(9,4)-r6(5,8)*two+r5(2,10)+r5(2,12)+r5(9,12)+r5(9,22)-r4(5,13)*two &
&             -r4(5,35)*two+r3(2,8)+r3(2,10)+r3(2,36)+r3(2,40)+r3(9,40)-r2(6,34)*two &
&             +r1(2,14)+r1(2,18)
      rxyz(17)=+r7(9,3)-r6(5,7)*two+r5(2,9)+r5(2,11)+r5(9,11)+r5(9,21)-r4(5,12)*two &
&             -r4(5,34)*two+r3(2,7)+r3(2,9)+r3(2,35)+r3(2,39)+r3(9,39)-r2(6,33)*two &
&             +r1(2,13)+r1(2,17)
      rxyz(18)=+r7(20,3)-r6(14,7)*two+r5(9,9)+r5(9,11)+r5(20,11)+r5(9,21)-r4(14,12)*two &
&             -r4(5,34)*two+r3(9,7)+r3(9,9)+r3(2,35)+r3(2,39)+r3(9,39)-r2(6,33)*two &
&             +r1(2,13)+r1(2,17)
      rxyz(19)=+r7(20,4)-r6(14,8)*two+r5(9,10)+r5(9,12)+r5(20,12)+r5(9,22)-r4(14,13)*two &
&             -r4(5,35)*two+r3(9,8)+r3(9,10)+r3(2,36)+r3(2,40)+r3(9,40)-r2(6,34)*two &
&             +r1(2,14)+r1(2,18)
      rxyz(20)=+r6(10,11)*four-r5(6,19)*eight+r4(3,24)*four+r4(3,28)*four+r4(10,28)*four &
&             -r3(6,25)*eight+r2(5,15)*four+r2(5,19)*four
      eri(1,1,1,3)=r400+(+r7(6,3)*two+r7(6,4)*two-r6(3,7)*four-r6(3,8)*four+r5(1,9)*two &
&                  +r5(1,10)*two+r5(1,11)*two+r5(6,11)*two+r5(1,12)*two+r5(6,12)*two &
&                  +r5(6,21)*six+r5(6,22)*six-r4(3,12)*four-r4(3,13)*four-r4(3,34)*p12 &
&                  -r4(3,35)*p12+r3(1,7)*two+r3(1,8)*two+r3(1,9)*two+r3(1,10)*two &
&                  +r3(1,35)*six+r3(1,36)*six+r3(1,39)*six+r3(6,39)*six+r3(1,40)*six &
&                  +r3(6,40)*six-r2(5,33)*p12-r2(5,34)*p12+r1(1,13)*six+r1(1,14)*six &
&                  +r1(1,17)*six+r1(1,18)*six)*qx+(+r6(6,10)+r6(6,11)*four+r6(6,12) &
&                  -r5(3,18)*two-r5(3,19)*eight-r5(3,20)*two+r4(1,23)+r4(1,24)*four+r4(1,25) &
&                  +r4(1,27)+r4(6,27)+r4(1,28)*four+r4(6,28)*four+r4(1,29)+r4(6,29)+r4(6,39) &
&                  +r4(6,40)*four+r4(6,41)-r3(3,24)*two-r3(3,25)*eight-r3(3,26)*two &
&                  -r3(3,49)*two-r3(3,50)*eight-r3(3,51)*two+r2(1,14)+r2(1,15)*four+r2(1,16) &
&                  +r2(1,18)+r2(1,19)*four+r2(1,20)+r2(1,48)+r2(1,49)*four+r2(1,50)+r2(1,53) &
&                  +r2(3,53)+r2(1,54)*four+r2(3,54)*four+r2(1,55)+r2(3,55)-r1(3,37)*two &
&                  -r1(3,38)*eight-r1(3,39)*two+r0(17)+r0(18)*four+r0(19)+r0(22)+r0(23)*four &
&                  +r0(24))*xx+(+r5(6,23)*two+r5(6,24)*two-r4(3,36)*four-r4(3,37)*four &
&                  +r3(1,37)*two+r3(1,38)*two+r3(1,41)*two+r3(6,41)*two+r3(1,42)*two &
&                  +r3(6,42)*two-r2(5,35)*four-r2(5,36)*four+r1(1,15)*two+r1(1,16)*two &
&                  +r1(1,19)*two+r1(1,20)*two)*xxx+rxyz(1)*xxxx
      eri(2,1,1,3)=r220+(+r7(13,4)*two-r6(8,8)*four+r5(4,10)*two+r5(4,12)*two &
&                  +r5(13,12)*two+r5(6,22)*two-r4(8,13)*four-r4(3,35)*four+r3(4,8)*two &
&                  +r3(4,10)*two+r3(1,36)*two+r3(1,40)*two+r3(6,40)*two-r2(5,34)*four &
&                  +r1(1,14)*two+r1(1,18)*two)*qx+rxyz(5)*xx
      eri(3,1,1,3)=r202+(+r7(15,4)*two-r6(10,8)*four+r5(6,10)*two+r5(6,12)*two &
&                  +r5(15,12)*two+r5(6,22)*two-r4(10,13)*four-r4(3,35)*four+r3(6,8)*two &
&                  +r3(6,10)*two+r3(1,36)*two+r3(1,40)*two+r3(6,40)*two-r2(5,34)*four &
&                  +r1(1,14)*two+r1(1,18)*two)*qx+(+r7(10,3)*two-r6(6,7)*four+r5(3,9)*two &
&                  +r5(3,11)*two+r5(10,11)*two+r5(10,21)*two-r4(6,12)*four-r4(6,34)*four &
&                  +r3(3,7)*two+r3(3,9)*two+r3(3,35)*two+r3(3,39)*two+r3(10,39)*two &
&                  -r2(3,33)*four+r1(3,13)*two+r1(3,17)*two)*qz+(+r6(15,12)-r5(10,20)*two &
&                  +r4(6,25)+r4(6,29)+r4(15,29)+r4(6,41)-r3(10,26)*two-r3(3,51)*two+r2(3,16) &
&                  +r2(3,20)+r2(1,50)+r2(1,55)+r2(3,55)-r1(3,39)*two+r0(19)+r0(24))*xx &
&                  +rxyz(20)*xz+(+r6(6,10)-r5(3,18)*two+r4(1,23)+r4(1,27)+r4(6,27)+r4(6,39) &
&                  -r3(3,24)*two-r3(3,49)*two+r2(1,14)+r2(1,18)+r2(1,48)+r2(1,53)+r2(3,53) &
&                  -r1(3,37)*two+r0(17)+r0(22))*zz+(+r5(10,24)*two-r4(6,37)*four+r3(3,38)*two &
&                  +r3(3,42)*two+r3(10,42)*two-r2(3,36)*four+r1(3,16)*two+r1(3,20)*two)*xxz+( &
&                  +r5(6,23)*two-r4(3,36)*four+r3(1,37)*two+r3(1,41)*two+r3(6,41)*two &
&                  -r2(5,35)*four+r1(1,15)*two+r1(1,19)*two)*xzz+rxyz(1)*xxzz
      eri(4,1,1,3)=r310+(+r7(9,3)+r7(9,4)*two-r6(5,7)*two-r6(5,8)*four+r5(2,9) &
&                  +r5(2,10)*two+r5(2,11)+r5(9,11)+r5(2,12)*two+r5(9,12)*two+r5(9,21) &
&                  +r5(9,22)*two-r4(5,12)*two-r4(5,13)*four-r4(5,34)*two-r4(5,35)*four &
&                  +r3(2,7)+r3(2,8)*two+r3(2,9)+r3(2,10)*two+r3(2,35)+r3(2,36)*two+r3(2,39) &
&                  +r3(9,39)+r3(2,40)*two+r3(9,40)*two-r2(6,33)*two-r2(6,34)*four+r1(2,13) &
&                  +r1(2,14)*two+r1(2,17)+r1(2,18)*two)*qx+(+r6(9,11)*two+r6(9,12) &
&                  -r5(5,19)*four-r5(5,20)*two+r4(2,24)*two+r4(2,25)+r4(2,28)*two &
&                  +r4(9,28)*two+r4(2,29)+r4(9,29)-r3(5,25)*four-r3(5,26)*two+r2(4,15)*two &
&                  +r2(4,16)+r2(4,19)*two+r2(4,20))*xx+rxyz(2)*xxx
      eri(5,1,1,3)=r301+(+r7(10,3)+r7(10,4)*two-r6(6,7)*two-r6(6,8)*four+r5(3,9) &
&                  +r5(3,10)*two+r5(3,11)+r5(10,11)+r5(3,12)*two+r5(10,12)*two+r5(10,21) &
&                  +r5(10,22)*two-r4(6,12)*two-r4(6,13)*four-r4(6,34)*two-r4(6,35)*four &
&                  +r3(3,7)+r3(3,8)*two+r3(3,9)+r3(3,10)*two+r3(3,35)+r3(3,36)*two+r3(3,39) &
&                  +r3(10,39)+r3(3,40)*two+r3(10,40)*two-r2(3,33)*two-r2(3,34)*four+r1(3,13) &
&                  +r1(3,14)*two+r1(3,17)+r1(3,18)*two)*qx+(+r7(6,3)-r6(3,7)*two+r5(1,9) &
&                  +r5(1,11)+r5(6,11)+r5(6,21)*three-r4(3,12)*two-r4(3,34)*six+r3(1,7) &
&                  +r3(1,9)+r3(1,35)*three+r3(1,39)*three+r3(6,39)*three-r2(5,33)*six &
&                  +r1(1,13)*three+r1(1,17)*three)*qz+(+r6(10,11)*two+r6(10,12)-r5(6,19)*four &
&                  -r5(6,20)*two+r4(3,24)*two+r4(3,25)+r4(3,28)*two+r4(10,28)*two+r4(3,29) &
&                  +r4(10,29)-r3(6,25)*four-r3(6,26)*two+r2(5,15)*two+r2(5,16)+r2(5,19)*two &
&                  +r2(5,20))*xx+(+r6(6,10)+r6(6,11)*two-r5(3,18)*two-r5(3,19)*four+r4(1,23) &
&                  +r4(1,24)*two+r4(1,27)+r4(6,27)+r4(1,28)*two+r4(6,28)*two+r4(6,39) &
&                  +r4(6,40)*two-r3(3,24)*two-r3(3,25)*four-r3(3,49)*two-r3(3,50)*four &
&                  +r2(1,14)+r2(1,15)*two+r2(1,18)+r2(1,19)*two+r2(1,48)+r2(1,49)*two &
&                  +r2(1,53)+r2(3,53)+r2(1,54)*two+r2(3,54)*two-r1(3,37)*two-r1(3,38)*four &
&                  +r0(17)+r0(18)*two+r0(22)+r0(23)*two)*xz+(+r5(10,24)-r4(6,37)*two+r3(3,38) &
&                  +r3(3,42)+r3(10,42)-r2(3,36)*two+r1(3,16)+r1(3,20))*xxx+(+r5(6,23)*two &
&                  +r5(6,24)-r4(3,36)*four-r4(3,37)*two+r3(1,37)*two+r3(1,38)+r3(1,41)*two &
&                  +r3(6,41)*two+r3(1,42)+r3(6,42)-r2(5,35)*four-r2(5,36)*two+r1(1,15)*two &
&                  +r1(1,16)+r1(1,19)*two+r1(1,20))*xxz+rxyz(1)*xxxz
      eri(6,1,1,3)=r211+(+r7(14,4)*two-r6(9,8)*four+r5(5,10)*two+r5(5,12)*two &
&                  +r5(14,12)*two-r4(9,13)*four+r3(5,8)*two+r3(5,10)*two)*qx+rxyz(17)*qz+( &
&                  +r6(14,12)-r5(9,20)*two+r4(5,25)+r4(5,29)+r4(14,29)-r3(9,26)*two+r2(6,16) &
&                  +r2(6,20))*xx+(+r6(9,11)*two-r5(5,19)*four+r4(2,24)*two+r4(2,28)*two &
&                  +r4(9,28)*two-r3(5,25)*four+r2(4,15)*two+r2(4,19)*two)*xz+rxyz(2)*xxz
      eri(1,2,1,3)=r220+(+r7(13,3)*two-r6(8,7)*four+r5(4,9)*two+r5(4,11)*two &
&                  +r5(13,11)*two+r5(6,21)*two-r4(8,12)*four-r4(3,34)*four+r3(4,7)*two &
&                  +r3(4,9)*two+r3(1,35)*two+r3(1,39)*two+r3(6,39)*two-r2(5,33)*four &
&                  +r1(1,13)*two+r1(1,17)*two)*qx+rxyz(4)*xx
      eri(2,2,1,3)=r040
      eri(3,2,1,3)=r022+(+r7(19,3)*two-r6(13,7)*four+r5(8,9)*two+r5(8,11)*two &
&                  +r5(19,11)*two+r5(10,21)*two-r4(13,12)*four-r4(6,34)*four+r3(8,7)*two &
&                  +r3(8,9)*two+r3(3,35)*two+r3(3,39)*two+r3(10,39)*two-r2(3,33)*four &
&                  +r1(3,13)*two+r1(3,17)*two)*qz+rxyz(4)*zz
      eri(4,2,1,3)=r130+rxyz(7)*qx
      eri(5,2,1,3)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,1,3)=r031+rxyz(7)*qz
      eri(1,3,1,3)=r202+(+r7(15,3)*two-r6(10,7)*four+r5(6,9)*two+r5(6,11)*two &
&                  +r5(15,11)*two+r5(6,21)*two-r4(10,12)*four-r4(3,34)*four+r3(6,7)*two &
&                  +r3(6,9)*two+r3(1,35)*two+r3(1,39)*two+r3(6,39)*two-r2(5,33)*four &
&                  +r1(1,13)*two+r1(1,17)*two)*qx+(+r7(10,4)*two-r6(6,8)*four+r5(3,10)*two &
&                  +r5(3,12)*two+r5(10,12)*two+r5(10,22)*two-r4(6,13)*four-r4(6,35)*four &
&                  +r3(3,8)*two+r3(3,10)*two+r3(3,36)*two+r3(3,40)*two+r3(10,40)*two &
&                  -r2(3,34)*four+r1(3,14)*two+r1(3,18)*two)*qz+(+r6(15,10)-r5(10,18)*two &
&                  +r4(6,23)+r4(6,27)+r4(15,27)+r4(6,39)-r3(10,24)*two-r3(3,49)*two+r2(3,14) &
&                  +r2(3,18)+r2(1,48)+r2(1,53)+r2(3,53)-r1(3,37)*two+r0(17)+r0(22))*xx &
&                  +rxyz(20)*xz+(+r6(6,12)-r5(3,20)*two+r4(1,25)+r4(1,29)+r4(6,29)+r4(6,41) &
&                  -r3(3,26)*two-r3(3,51)*two+r2(1,16)+r2(1,20)+r2(1,50)+r2(1,55)+r2(3,55) &
&                  -r1(3,39)*two+r0(19)+r0(24))*zz+(+r5(10,23)*two-r4(6,36)*four+r3(3,37)*two &
&                  +r3(3,41)*two+r3(10,41)*two-r2(3,35)*four+r1(3,15)*two+r1(3,19)*two)*xxz+( &
&                  +r5(6,24)*two-r4(3,37)*four+r3(1,38)*two+r3(1,42)*two+r3(6,42)*two &
&                  -r2(5,36)*four+r1(1,16)*two+r1(1,20)*two)*xzz+rxyz(1)*xxzz
      eri(2,3,1,3)=r022+(+r7(19,4)*two-r6(13,8)*four+r5(8,10)*two+r5(8,12)*two &
&                  +r5(19,12)*two+r5(10,22)*two-r4(13,13)*four-r4(6,35)*four+r3(8,8)*two &
&                  +r3(8,10)*two+r3(3,36)*two+r3(3,40)*two+r3(10,40)*two-r2(3,34)*four &
&                  +r1(3,14)*two+r1(3,18)*two)*qz+rxyz(5)*zz
      eri(3,3,1,3)=r004+(+r7(21,3)*two+r7(21,4)*two-r6(15,7)*four-r6(15,8)*four &
&                  +r5(10,9)*two+r5(10,10)*two+r5(10,11)*two+r5(21,11)*two+r5(10,12)*two &
&                  +r5(21,12)*two+r5(10,21)*six+r5(10,22)*six-r4(15,12)*four-r4(15,13)*four &
&                  -r4(6,34)*p12-r4(6,35)*p12+r3(10,7)*two+r3(10,8)*two+r3(10,9)*two &
&                  +r3(10,10)*two+r3(3,35)*six+r3(3,36)*six+r3(3,39)*six+r3(10,39)*six &
&                  +r3(3,40)*six+r3(10,40)*six-r2(3,33)*p12-r2(3,34)*p12+r1(3,13)*six &
&                  +r1(3,14)*six+r1(3,17)*six+r1(3,18)*six)*qz+(+r6(15,10)+r6(15,11)*four &
&                  +r6(15,12)-r5(10,18)*two-r5(10,19)*eight-r5(10,20)*two+r4(6,23) &
&                  +r4(6,24)*four+r4(6,25)+r4(6,27)+r4(15,27)+r4(6,28)*four+r4(15,28)*four &
&                  +r4(6,29)+r4(15,29)+r4(6,39)+r4(6,40)*four+r4(6,41)-r3(10,24)*two &
&                  -r3(10,25)*eight-r3(10,26)*two-r3(3,49)*two-r3(3,50)*eight-r3(3,51)*two &
&                  +r2(3,14)+r2(3,15)*four+r2(3,16)+r2(3,18)+r2(3,19)*four+r2(3,20)+r2(1,48) &
&                  +r2(1,49)*four+r2(1,50)+r2(1,53)+r2(3,53)+r2(1,54)*four+r2(3,54)*four &
&                  +r2(1,55)+r2(3,55)-r1(3,37)*two-r1(3,38)*eight-r1(3,39)*two+r0(17) &
&                  +r0(18)*four+r0(19)+r0(22)+r0(23)*four+r0(24))*zz+(+r5(10,23)*two &
&                  +r5(10,24)*two-r4(6,36)*four-r4(6,37)*four+r3(3,37)*two+r3(3,38)*two &
&                  +r3(3,41)*two+r3(10,41)*two+r3(3,42)*two+r3(10,42)*two-r2(3,35)*four &
&                  -r2(3,36)*four+r1(3,15)*two+r1(3,16)*two+r1(3,19)*two+r1(3,20)*two)*zzz &
&                  +rxyz(1)*zzzz
      eri(4,3,1,3)=r112+rxyz(18)*qx+(+r7(14,4)*two-r6(9,8)*four+r5(5,10)*two+r5(5,12)*two &
&                  +r5(14,12)*two-r4(9,13)*four+r3(5,8)*two+r3(5,10)*two)*qz+(+r6(14,11)*two &
&                  -r5(9,19)*four+r4(5,24)*two+r4(5,28)*two+r4(14,28)*two-r3(9,25)*four &
&                  +r2(6,15)*two+r2(6,19)*two)*xz+(+r6(9,12)-r5(5,20)*two+r4(2,25)+r4(2,29) &
&                  +r4(9,29)-r3(5,26)*two+r2(4,16)+r2(4,20))*zz+rxyz(2)*xzz
      eri(5,3,1,3)=r103+(+r7(21,3)-r6(15,7)*two+r5(10,9)+r5(10,11)+r5(21,11) &
&                  +r5(10,21)*three-r4(15,12)*two-r4(6,34)*six+r3(10,7)+r3(10,9) &
&                  +r3(3,35)*three+r3(3,39)*three+r3(10,39)*three-r2(3,33)*six+r1(3,13)*three &
&                  +r1(3,17)*three)*qx+(+r7(15,3)+r7(15,4)*two-r6(10,7)*two-r6(10,8)*four &
&                  +r5(6,9)+r5(6,10)*two+r5(6,11)+r5(15,11)+r5(6,12)*two+r5(15,12)*two &
&                  +r5(6,21)+r5(6,22)*two-r4(10,12)*two-r4(10,13)*four-r4(3,34)*two &
&                  -r4(3,35)*four+r3(6,7)+r3(6,8)*two+r3(6,9)+r3(6,10)*two+r3(1,35) &
&                  +r3(1,36)*two+r3(1,39)+r3(6,39)+r3(1,40)*two+r3(6,40)*two-r2(5,33)*two &
&                  -r2(5,34)*four+r1(1,13)+r1(1,14)*two+r1(1,17)+r1(1,18)*two)*qz+(+r6(15,10) &
&                  +r6(15,11)*two-r5(10,18)*two-r5(10,19)*four+r4(6,23)+r4(6,24)*two+r4(6,27) &
&                  +r4(15,27)+r4(6,28)*two+r4(15,28)*two+r4(6,39)+r4(6,40)*two-r3(10,24)*two &
&                  -r3(10,25)*four-r3(3,49)*two-r3(3,50)*four+r2(3,14)+r2(3,15)*two+r2(3,18) &
&                  +r2(3,19)*two+r2(1,48)+r2(1,49)*two+r2(1,53)+r2(3,53)+r2(1,54)*two &
&                  +r2(3,54)*two-r1(3,37)*two-r1(3,38)*four+r0(17)+r0(18)*two+r0(22) &
&                  +r0(23)*two)*xz+(+r6(10,11)*two+r6(10,12)-r5(6,19)*four-r5(6,20)*two &
&                  +r4(3,24)*two+r4(3,25)+r4(3,28)*two+r4(10,28)*two+r4(3,29)+r4(10,29) &
&                  -r3(6,25)*four-r3(6,26)*two+r2(5,15)*two+r2(5,16)+r2(5,19)*two+r2(5,20)) &
&                  *zz+(+r5(10,23)*two+r5(10,24)-r4(6,36)*four-r4(6,37)*two+r3(3,37)*two &
&                  +r3(3,38)+r3(3,41)*two+r3(10,41)*two+r3(3,42)+r3(10,42)-r2(3,35)*four &
&                  -r2(3,36)*two+r1(3,15)*two+r1(3,16)+r1(3,19)*two+r1(3,20))*xzz+(+r5(6,24) &
&                  -r4(3,37)*two+r3(1,38)+r3(1,42)+r3(6,42)-r2(5,36)*two+r1(1,16)+r1(1,20)) &
&                  *zzz+rxyz(1)*xzzz
      eri(6,3,1,3)=r013+(+r7(20,3)+r7(20,4)*two-r6(14,7)*two-r6(14,8)*four+r5(9,9) &
&                  +r5(9,10)*two+r5(9,11)+r5(20,11)+r5(9,12)*two+r5(20,12)*two+r5(9,21) &
&                  +r5(9,22)*two-r4(14,12)*two-r4(14,13)*four-r4(5,34)*two-r4(5,35)*four &
&                  +r3(9,7)+r3(9,8)*two+r3(9,9)+r3(9,10)*two+r3(2,35)+r3(2,36)*two+r3(2,39) &
&                  +r3(9,39)+r3(2,40)*two+r3(9,40)*two-r2(6,33)*two-r2(6,34)*four+r1(2,13) &
&                  +r1(2,14)*two+r1(2,17)+r1(2,18)*two)*qz+(+r6(14,11)*two+r6(14,12) &
&                  -r5(9,19)*four-r5(9,20)*two+r4(5,24)*two+r4(5,25)+r4(5,28)*two &
&                  +r4(14,28)*two+r4(5,29)+r4(14,29)-r3(9,25)*four-r3(9,26)*two+r2(6,15)*two &
&                  +r2(6,16)+r2(6,19)*two+r2(6,20))*zz+rxyz(2)*zzz
      eri(1,4,1,3)=r310+(+r7(9,3)*two+r7(9,4)-r6(5,7)*four-r6(5,8)*two+r5(2,9)*two &
&                  +r5(2,10)+r5(2,11)*two+r5(9,11)*two+r5(2,12)+r5(9,12)+r5(9,21)*two &
&                  +r5(9,22)-r4(5,12)*four-r4(5,13)*two-r4(5,34)*four-r4(5,35)*two &
&                  +r3(2,7)*two+r3(2,8)+r3(2,9)*two+r3(2,10)+r3(2,35)*two+r3(2,36) &
&                  +r3(2,39)*two+r3(9,39)*two+r3(2,40)+r3(9,40)-r2(6,33)*four-r2(6,34)*two &
&                  +r1(2,13)*two+r1(2,14)+r1(2,17)*two+r1(2,18))*qx+(+r6(9,10)+r6(9,11)*two &
&                  -r5(5,18)*two-r5(5,19)*four+r4(2,23)+r4(2,24)*two+r4(2,27)+r4(9,27) &
&                  +r4(2,28)*two+r4(9,28)*two-r3(5,24)*two-r3(5,25)*four+r2(4,14) &
&                  +r2(4,15)*two+r2(4,18)+r2(4,19)*two)*xx+rxyz(3)*xxx
      eri(2,4,1,3)=r130+rxyz(8)*qx
      eri(3,4,1,3)=r112+rxyz(19)*qx+(+r7(14,3)*two-r6(9,7)*four+r5(5,9)*two+r5(5,11)*two &
&                  +r5(14,11)*two-r4(9,12)*four+r3(5,7)*two+r3(5,9)*two)*qz+(+r6(14,11)*two &
&                  -r5(9,19)*four+r4(5,24)*two+r4(5,28)*two+r4(14,28)*two-r3(9,25)*four &
&                  +r2(6,15)*two+r2(6,19)*two)*xz+(+r6(9,10)-r5(5,18)*two+r4(2,23)+r4(2,27) &
&                  +r4(9,27)-r3(5,24)*two+r2(4,14)+r2(4,18))*zz+rxyz(3)*xzz
      eri(4,4,1,3)=r220+(+r7(13,3)+r7(13,4)-r6(8,7)*two-r6(8,8)*two+r5(4,9)+r5(4,10) &
&                  +r5(4,11)+r5(13,11)+r5(4,12)+r5(13,12)+r5(6,21)+r5(6,22)-r4(8,12)*two &
&                  -r4(8,13)*two-r4(3,34)*two-r4(3,35)*two+r3(4,7)+r3(4,8)+r3(4,9)+r3(4,10) &
&                  +r3(1,35)+r3(1,36)+r3(1,39)+r3(6,39)+r3(1,40)+r3(6,40)-r2(5,33)*two &
&                  -r2(5,34)*two+r1(1,13)+r1(1,14)+r1(1,17)+r1(1,18))*qx+rxyz(6)*xx
      eri(5,4,1,3)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(9,10)+r6(9,11) &
&                  -r5(5,18)*two-r5(5,19)*two+r4(2,23)+r4(2,24)+r4(2,27)+r4(9,27)+r4(2,28) &
&                  +r4(9,28)-r3(5,24)*two-r3(5,25)*two+r2(4,14)+r2(4,15)+r2(4,18)+r2(4,19)) &
&                  *xz+rxyz(3)*xxz
      eri(6,4,1,3)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,1,3)=r301+(+r7(10,3)*two+r7(10,4)-r6(6,7)*four-r6(6,8)*two+r5(3,9)*two &
&                  +r5(3,10)+r5(3,11)*two+r5(10,11)*two+r5(3,12)+r5(10,12)+r5(10,21)*two &
&                  +r5(10,22)-r4(6,12)*four-r4(6,13)*two-r4(6,34)*four-r4(6,35)*two &
&                  +r3(3,7)*two+r3(3,8)+r3(3,9)*two+r3(3,10)+r3(3,35)*two+r3(3,36) &
&                  +r3(3,39)*two+r3(10,39)*two+r3(3,40)+r3(10,40)-r2(3,33)*four-r2(3,34)*two &
&                  +r1(3,13)*two+r1(3,14)+r1(3,17)*two+r1(3,18))*qx+(+r7(6,4)-r6(3,8)*two &
&                  +r5(1,10)+r5(1,12)+r5(6,12)+r5(6,22)*three-r4(3,13)*two-r4(3,35)*six &
&                  +r3(1,8)+r3(1,10)+r3(1,36)*three+r3(1,40)*three+r3(6,40)*three &
&                  -r2(5,34)*six+r1(1,14)*three+r1(1,18)*three)*qz+(+r6(10,10)+r6(10,11)*two &
&                  -r5(6,18)*two-r5(6,19)*four+r4(3,23)+r4(3,24)*two+r4(3,27)+r4(10,27) &
&                  +r4(3,28)*two+r4(10,28)*two-r3(6,24)*two-r3(6,25)*four+r2(5,14) &
&                  +r2(5,15)*two+r2(5,18)+r2(5,19)*two)*xx+(+r6(6,11)*two+r6(6,12) &
&                  -r5(3,19)*four-r5(3,20)*two+r4(1,24)*two+r4(1,25)+r4(1,28)*two &
&                  +r4(6,28)*two+r4(1,29)+r4(6,29)+r4(6,40)*two+r4(6,41)-r3(3,25)*four &
&                  -r3(3,26)*two-r3(3,50)*four-r3(3,51)*two+r2(1,15)*two+r2(1,16) &
&                  +r2(1,19)*two+r2(1,20)+r2(1,49)*two+r2(1,50)+r2(1,54)*two+r2(3,54)*two &
&                  +r2(1,55)+r2(3,55)-r1(3,38)*four-r1(3,39)*two+r0(18)*two+r0(19)+r0(23)*two &
&                  +r0(24))*xz+(+r5(10,23)-r4(6,36)*two+r3(3,37)+r3(3,41)+r3(10,41) &
&                  -r2(3,35)*two+r1(3,15)+r1(3,19))*xxx+(+r5(6,23)+r5(6,24)*two-r4(3,36)*two &
&                  -r4(3,37)*four+r3(1,37)+r3(1,38)*two+r3(1,41)+r3(6,41)+r3(1,42)*two &
&                  +r3(6,42)*two-r2(5,35)*two-r2(5,36)*four+r1(1,15)+r1(1,16)*two+r1(1,19) &
&                  +r1(1,20)*two)*xxz+rxyz(1)*xxxz
      eri(2,5,1,3)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,1,3)=r103+(+r7(21,4)-r6(15,8)*two+r5(10,10)+r5(10,12)+r5(21,12) &
&                  +r5(10,22)*three-r4(15,13)*two-r4(6,35)*six+r3(10,8)+r3(10,10) &
&                  +r3(3,36)*three+r3(3,40)*three+r3(10,40)*three-r2(3,34)*six+r1(3,14)*three &
&                  +r1(3,18)*three)*qx+(+r7(15,3)*two+r7(15,4)-r6(10,7)*four-r6(10,8)*two &
&                  +r5(6,9)*two+r5(6,10)+r5(6,11)*two+r5(15,11)*two+r5(6,12)+r5(15,12) &
&                  +r5(6,21)*two+r5(6,22)-r4(10,12)*four-r4(10,13)*two-r4(3,34)*four &
&                  -r4(3,35)*two+r3(6,7)*two+r3(6,8)+r3(6,9)*two+r3(6,10)+r3(1,35)*two &
&                  +r3(1,36)+r3(1,39)*two+r3(6,39)*two+r3(1,40)+r3(6,40)-r2(5,33)*four &
&                  -r2(5,34)*two+r1(1,13)*two+r1(1,14)+r1(1,17)*two+r1(1,18))*qz+( &
&                  +r6(15,11)*two+r6(15,12)-r5(10,19)*four-r5(10,20)*two+r4(6,24)*two &
&                  +r4(6,25)+r4(6,28)*two+r4(15,28)*two+r4(6,29)+r4(15,29)+r4(6,40)*two &
&                  +r4(6,41)-r3(10,25)*four-r3(10,26)*two-r3(3,50)*four-r3(3,51)*two &
&                  +r2(3,15)*two+r2(3,16)+r2(3,19)*two+r2(3,20)+r2(1,49)*two+r2(1,50) &
&                  +r2(1,54)*two+r2(3,54)*two+r2(1,55)+r2(3,55)-r1(3,38)*four-r1(3,39)*two &
&                  +r0(18)*two+r0(19)+r0(23)*two+r0(24))*xz+(+r6(10,10)+r6(10,11)*two &
&                  -r5(6,18)*two-r5(6,19)*four+r4(3,23)+r4(3,24)*two+r4(3,27)+r4(10,27) &
&                  +r4(3,28)*two+r4(10,28)*two-r3(6,24)*two-r3(6,25)*four+r2(5,14) &
&                  +r2(5,15)*two+r2(5,18)+r2(5,19)*two)*zz+(+r5(10,23)+r5(10,24)*two &
&                  -r4(6,36)*two-r4(6,37)*four+r3(3,37)+r3(3,38)*two+r3(3,41)+r3(10,41) &
&                  +r3(3,42)*two+r3(10,42)*two-r2(3,35)*two-r2(3,36)*four+r1(3,15) &
&                  +r1(3,16)*two+r1(3,19)+r1(3,20)*two)*xzz+(+r5(6,23)-r4(3,36)*two+r3(1,37) &
&                  +r3(1,41)+r3(6,41)-r2(5,35)*two+r1(1,15)+r1(1,19))*zzz+rxyz(1)*xzzz
      eri(4,5,1,3)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(9,11)+r6(9,12) &
&                  -r5(5,19)*two-r5(5,20)*two+r4(2,24)+r4(2,25)+r4(2,28)+r4(9,28)+r4(2,29) &
&                  +r4(9,29)-r3(5,25)*two-r3(5,26)*two+r2(4,15)+r2(4,16)+r2(4,19)+r2(4,20)) &
&                  *xz+rxyz(2)*xxz
      eri(5,5,1,3)=r202+(+r7(15,3)+r7(15,4)-r6(10,7)*two-r6(10,8)*two+r5(6,9)+r5(6,10) &
&                  +r5(6,11)+r5(15,11)+r5(6,12)+r5(15,12)+r5(6,21)+r5(6,22)-r4(10,12)*two &
&                  -r4(10,13)*two-r4(3,34)*two-r4(3,35)*two+r3(6,7)+r3(6,8)+r3(6,9)+r3(6,10) &
&                  +r3(1,35)+r3(1,36)+r3(1,39)+r3(6,39)+r3(1,40)+r3(6,40)-r2(5,33)*two &
&                  -r2(5,34)*two+r1(1,13)+r1(1,14)+r1(1,17)+r1(1,18))*qx+(+r7(10,3)+r7(10,4) &
&                  -r6(6,7)*two-r6(6,8)*two+r5(3,9)+r5(3,10)+r5(3,11)+r5(10,11)+r5(3,12) &
&                  +r5(10,12)+r5(10,21)+r5(10,22)-r4(6,12)*two-r4(6,13)*two-r4(6,34)*two &
&                  -r4(6,35)*two+r3(3,7)+r3(3,8)+r3(3,9)+r3(3,10)+r3(3,35)+r3(3,36)+r3(3,39) &
&                  +r3(10,39)+r3(3,40)+r3(10,40)-r2(3,33)*two-r2(3,34)*two+r1(3,13)+r1(3,14) &
&                  +r1(3,17)+r1(3,18))*qz+(+r6(15,11)-r5(10,19)*two+r4(6,24)+r4(6,28) &
&                  +r4(15,28)+r4(6,40)-r3(10,25)*two-r3(3,50)*two+r2(3,15)+r2(3,19)+r2(1,49) &
&                  +r2(1,54)+r2(3,54)-r1(3,38)*two+r0(18)+r0(23))*xx+(+r6(10,10) &
&                  +r6(10,11)*two+r6(10,12)-r5(6,18)*two-r5(6,19)*four-r5(6,20)*two+r4(3,23) &
&                  +r4(3,24)*two+r4(3,25)+r4(3,27)+r4(10,27)+r4(3,28)*two+r4(10,28)*two &
&                  +r4(3,29)+r4(10,29)-r3(6,24)*two-r3(6,25)*four-r3(6,26)*two+r2(5,14) &
&                  +r2(5,15)*two+r2(5,16)+r2(5,18)+r2(5,19)*two+r2(5,20))*xz+(+r6(6,11) &
&                  -r5(3,19)*two+r4(1,24)+r4(1,28)+r4(6,28)+r4(6,40)-r3(3,25)*two &
&                  -r3(3,50)*two+r2(1,15)+r2(1,19)+r2(1,49)+r2(1,54)+r2(3,54)-r1(3,38)*two &
&                  +r0(18)+r0(23))*zz+(+r5(10,23)+r5(10,24)-r4(6,36)*two-r4(6,37)*two &
&                  +r3(3,37)+r3(3,38)+r3(3,41)+r3(10,41)+r3(3,42)+r3(10,42)-r2(3,35)*two &
&                  -r2(3,36)*two+r1(3,15)+r1(3,16)+r1(3,19)+r1(3,20))*xxz+(+r5(6,23)+r5(6,24) &
&                  -r4(3,36)*two-r4(3,37)*two+r3(1,37)+r3(1,38)+r3(1,41)+r3(6,41)+r3(1,42) &
&                  +r3(6,42)-r2(5,35)*two-r2(5,36)*two+r1(1,15)+r1(1,16)+r1(1,19)+r1(1,20)) &
&                  *xzz+rxyz(1)*xxzz
      eri(6,5,1,3)=r112+rxyz(19)*qx+(+r7(14,3)+r7(14,4)-r6(9,7)*two-r6(9,8)*two+r5(5,9) &
&                  +r5(5,10)+r5(5,11)+r5(14,11)+r5(5,12)+r5(14,12)-r4(9,12)*two-r4(9,13)*two &
&                  +r3(5,7)+r3(5,8)+r3(5,9)+r3(5,10))*qz+(+r6(14,11)+r6(14,12)-r5(9,19)*two &
&                  -r5(9,20)*two+r4(5,24)+r4(5,25)+r4(5,28)+r4(14,28)+r4(5,29)+r4(14,29) &
&                  -r3(9,25)*two-r3(9,26)*two+r2(6,15)+r2(6,16)+r2(6,19)+r2(6,20))*xz &
&                  +rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,1,3)=r211+(+r7(14,3)*two-r6(9,7)*four+r5(5,9)*two+r5(5,11)*two &
&                  +r5(14,11)*two-r4(9,12)*four+r3(5,7)*two+r3(5,9)*two)*qx+rxyz(16)*qz+( &
&                  +r6(14,10)-r5(9,18)*two+r4(5,23)+r4(5,27)+r4(14,27)-r3(9,24)*two+r2(6,14) &
&                  +r2(6,18))*xx+(+r6(9,11)*two-r5(5,19)*four+r4(2,24)*two+r4(2,28)*two &
&                  +r4(9,28)*two-r3(5,25)*four+r2(4,15)*two+r2(4,19)*two)*xz+rxyz(3)*xxz
      eri(2,6,1,3)=r031+rxyz(8)*qz
      eri(3,6,1,3)=r013+(+r7(20,3)*two+r7(20,4)-r6(14,7)*four-r6(14,8)*two+r5(9,9)*two &
&                  +r5(9,10)+r5(9,11)*two+r5(20,11)*two+r5(9,12)+r5(20,12)+r5(9,21)*two &
&                  +r5(9,22)-r4(14,12)*four-r4(14,13)*two-r4(5,34)*four-r4(5,35)*two &
&                  +r3(9,7)*two+r3(9,8)+r3(9,9)*two+r3(9,10)+r3(2,35)*two+r3(2,36) &
&                  +r3(2,39)*two+r3(9,39)*two+r3(2,40)+r3(9,40)-r2(6,33)*four-r2(6,34)*two &
&                  +r1(2,13)*two+r1(2,14)+r1(2,17)*two+r1(2,18))*qz+(+r6(14,10)+r6(14,11)*two &
&                  -r5(9,18)*two-r5(9,19)*four+r4(5,23)+r4(5,24)*two+r4(5,27)+r4(14,27) &
&                  +r4(5,28)*two+r4(14,28)*two-r3(9,24)*two-r3(9,25)*four+r2(6,14) &
&                  +r2(6,15)*two+r2(6,18)+r2(6,19)*two)*zz+rxyz(3)*zzz
      eri(4,6,1,3)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,1,3)=r112+rxyz(18)*qx+(+r7(14,3)+r7(14,4)-r6(9,7)*two-r6(9,8)*two+r5(5,9) &
&                  +r5(5,10)+r5(5,11)+r5(14,11)+r5(5,12)+r5(14,12)-r4(9,12)*two-r4(9,13)*two &
&                  +r3(5,7)+r3(5,8)+r3(5,9)+r3(5,10))*qz+(+r6(14,10)+r6(14,11)-r5(9,18)*two &
&                  -r5(9,19)*two+r4(5,23)+r4(5,24)+r4(5,27)+r4(14,27)+r4(5,28)+r4(14,28) &
&                  -r3(9,24)*two-r3(9,25)*two+r2(6,14)+r2(6,15)+r2(6,18)+r2(6,19))*xz &
&                  +rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,1,3)=r022+(+r7(19,3)+r7(19,4)-r6(13,7)*two-r6(13,8)*two+r5(8,9)+r5(8,10) &
&                  +r5(8,11)+r5(19,11)+r5(8,12)+r5(19,12)+r5(10,21)+r5(10,22)-r4(13,12)*two &
&                  -r4(13,13)*two-r4(6,34)*two-r4(6,35)*two+r3(8,7)+r3(8,8)+r3(8,9)+r3(8,10) &
&                  +r3(3,35)+r3(3,36)+r3(3,39)+r3(10,39)+r3(3,40)+r3(10,40)-r2(3,33)*two &
&                  -r2(3,34)*two+r1(3,13)+r1(3,14)+r1(3,17)+r1(3,18))*qz+rxyz(6)*zz
!
      r400= r8(13)-r7(8,2)*two+r6(4,3)+r6(4,4)+r6(6,4)+r6(13,9)*six-r5(3,4)*two &
&          -r5(8,17)*p12+r4(1,4)+r4(1,5)+r4(4,22)*six+r4(4,26)*six+r4(6,26)*six &
&          +r4(13,38)*three-r3(3,23)*p12-r3(8,48)*six+r2(1,13)*six+r2(1,17)*six &
&          +r2(2,47)*three+r2(2,52)*three+r2(3,52)*three-r1(3,36)*six+r0(16)*three &
&          +r0(21)*three
      r310= r8(18)-r7(12,2)*two+r6(7,3)+r6(7,4)+r6(9,4)+r6(18,9)*three-r5(5,4)*two &
&          -r5(12,17)*six+r4(2,4)+r4(2,5)+r4(7,22)*three+r4(7,26)*three+r4(9,26)*three &
&          -r3(5,23)*six+r2(4,13)*three+r2(4,17)*three
      r301= r8(19)-r7(13,2)*two+r6(8,3)+r6(8,4)+r6(10,4)+r6(19,9)*three-r5(6,4)*two &
&          -r5(13,17)*six+r4(3,4)+r4(3,5)+r4(8,22)*three+r4(8,26)*three+r4(10,26)*three &
&          -r3(6,23)*six+r2(5,13)*three+r2(5,17)*three
      r220= r8(24)-r7(17,2)*two+r6(11,3)+r6(11,4)+r6(13,4)+r6(13,9)+r6(24,9)-r5(8,4)*two &
&          -r5(8,17)*two-r5(17,17)*two+r4(4,4)+r4(4,5)+r4(4,22)+r4(11,22)+r4(4,26) &
&          +r4(6,26)+r4(11,26)+r4(13,26)+r4(13,38)-r3(3,23)*two-r3(8,23)*two-r3(8,48)*two &
&          +r2(1,13)+r2(2,13)+r2(1,17)+r2(2,17)+r2(2,47)+r2(2,52)+r2(3,52)-r1(3,36)*two &
&          +r0(16)+r0(21)
      r211= r8(25)-r7(18,2)*two+r6(12,3)+r6(12,4)+r6(14,4)+r6(25,9)-r5(9,4)*two &
&          -r5(18,17)*two+r4(5,4)+r4(5,5)+r4(12,22)+r4(12,26)+r4(14,26)-r3(9,23)*two &
&          +r2(6,13)+r2(6,17)
      r202= r8(26)-r7(19,2)*two+r6(13,3)+r6(13,4)+r6(15,4)+r6(13,9)+r6(26,9)-r5(10,4)*two &
&          -r5(8,17)*two-r5(19,17)*two+r4(6,4)+r4(6,5)+r4(4,22)+r4(13,22)+r4(4,26) &
&          +r4(6,26)+r4(13,26)+r4(15,26)+r4(13,38)-r3(3,23)*two-r3(10,23)*two-r3(8,48)*two &
&          +r2(1,13)+r2(3,13)+r2(1,17)+r2(3,17)+r2(2,47)+r2(2,52)+r2(3,52)-r1(3,36)*two &
&          +r0(16)+r0(21)
      r130= r8(31)-r7(23,2)*two+r6(16,3)+r6(16,4)+r6(18,4)+r6(18,9)*three-r5(12,4)*two &
&          -r5(12,17)*six+r4(7,4)+r4(7,5)+r4(7,22)*three+r4(7,26)*three+r4(9,26)*three &
&          -r3(5,23)*six+r2(4,13)*three+r2(4,17)*three
      r121= r8(32)-r7(24,2)*two+r6(17,3)+r6(17,4)+r6(19,4)+r6(19,9)-r5(13,4)*two &
&          -r5(13,17)*two+r4(8,4)+r4(8,5)+r4(8,22)+r4(8,26)+r4(10,26)-r3(6,23)*two &
&          +r2(5,13)+r2(5,17)
      r112= r8(33)-r7(25,2)*two+r6(18,3)+r6(18,4)+r6(20,4)+r6(18,9)-r5(14,4)*two &
&          -r5(12,17)*two+r4(9,4)+r4(9,5)+r4(7,22)+r4(7,26)+r4(9,26)-r3(5,23)*two+r2(4,13) &
&          +r2(4,17)
      r103= r8(34)-r7(26,2)*two+r6(19,3)+r6(19,4)+r6(21,4)+r6(19,9)*three-r5(15,4)*two &
&          -r5(13,17)*six+r4(10,4)+r4(10,5)+r4(8,22)*three+r4(8,26)*three+r4(10,26)*three &
&          -r3(6,23)*six+r2(5,13)*three+r2(5,17)*three
      r040= r8(39)-r7(30,2)*two+r6(22,3)+r6(22,4)+r6(24,4)+r6(24,9)*six-r5(17,4)*two &
&          -r5(17,17)*p12+r4(11,4)+r4(11,5)+r4(11,22)*six+r4(11,26)*six+r4(13,26)*six &
&          +r4(13,38)*three-r3(8,23)*p12-r3(8,48)*six+r2(2,13)*six+r2(2,17)*six &
&          +r2(2,47)*three+r2(2,52)*three+r2(3,52)*three-r1(3,36)*six+r0(16)*three &
&          +r0(21)*three
      r031= r8(40)-r7(31,2)*two+r6(23,3)+r6(23,4)+r6(25,4)+r6(25,9)*three-r5(18,4)*two &
&          -r5(18,17)*six+r4(12,4)+r4(12,5)+r4(12,22)*three+r4(12,26)*three &
&          +r4(14,26)*three-r3(9,23)*six+r2(6,13)*three+r2(6,17)*three
      r022= r8(41)-r7(32,2)*two+r6(24,3)+r6(24,4)+r6(26,4)+r6(24,9)+r6(26,9)-r5(19,4)*two &
&          -r5(17,17)*two-r5(19,17)*two+r4(13,4)+r4(13,5)+r4(11,22)+r4(13,22)+r4(11,26) &
&          +r4(13,26)*two+r4(15,26)+r4(13,38)-r3(8,23)*two-r3(10,23)*two-r3(8,48)*two &
&          +r2(2,13)+r2(3,13)+r2(2,17)+r2(3,17)+r2(2,47)+r2(2,52)+r2(3,52)-r1(3,36)*two &
&          +r0(16)+r0(21)
      r013= r8(42)-r7(33,2)*two+r6(25,3)+r6(25,4)+r6(27,4)+r6(25,9)*three-r5(20,4)*two &
&          -r5(18,17)*six+r4(14,4)+r4(14,5)+r4(12,22)*three+r4(12,26)*three &
&          +r4(14,26)*three-r3(9,23)*six+r2(6,13)*three+r2(6,17)*three
      r004= r8(43)-r7(34,2)*two+r6(26,3)+r6(26,4)+r6(28,4)+r6(26,9)*six-r5(21,4)*two &
&          -r5(19,17)*p12+r4(15,4)+r4(15,5)+r4(13,22)*six+r4(13,26)*six+r4(15,26)*six &
&          +r4(13,38)*three-r3(10,23)*p12-r3(8,48)*six+r2(3,13)*six+r2(3,17)*six &
&          +r2(2,47)*three+r2(2,52)*three+r2(3,52)*three-r1(3,36)*six+r0(16)*three &
&          +r0(21)*three
      rxyz(1)=+r4(13,42)-r3(8,52)*two+r2(2,51)+r2(2,56)+r2(3,56)-r1(3,40)*two+r0(20) &
&             +r0(25)
      rxyz(2)=+r5(18,24)-r4(12,37)*two+r3(7,38)+r3(7,42)+r3(9,42)-r2(6,36)*two+r1(2,16) &
&             +r1(2,20)
      rxyz(3)=+r5(18,23)-r4(12,36)*two+r3(7,37)+r3(7,41)+r3(9,41)-r2(6,35)*two+r1(2,15) &
&             +r1(2,19)
      rxyz(4)=+r6(24,10)-r5(17,18)*two+r4(11,23)+r4(11,27)+r4(13,27)+r4(13,39) &
&             -r3(8,24)*two-r3(8,49)*two+r2(2,14)+r2(2,18)+r2(2,48)+r2(2,53)+r2(3,53) &
&             -r1(3,37)*two+r0(17)+r0(22)
      rxyz(5)=+r6(24,12)-r5(17,20)*two+r4(11,25)+r4(11,29)+r4(13,29)+r4(13,41) &
&             -r3(8,26)*two-r3(8,51)*two+r2(2,16)+r2(2,20)+r2(2,50)+r2(2,55)+r2(3,55) &
&             -r1(3,39)*two+r0(19)+r0(24)
      rxyz(6)=+r6(24,11)-r5(17,19)*two+r4(11,24)+r4(11,28)+r4(13,28)+r4(13,40) &
&             -r3(8,25)*two-r3(8,50)*two+r2(2,15)+r2(2,19)+r2(2,49)+r2(2,54)+r2(3,54) &
&             -r1(3,38)*two+r0(18)+r0(23)
      rxyz(7)=+r7(31,3)-r6(23,7)*two+r5(16,9)+r5(16,11)+r5(18,11)+r5(18,21)*three &
&             -r4(12,12)*two-r4(12,34)*six+r3(7,7)+r3(7,9)+r3(7,35)*three+r3(7,39)*three &
&             +r3(9,39)*three-r2(6,33)*six+r1(2,13)*three+r1(2,17)*three
      rxyz(8)=+r7(31,4)-r6(23,8)*two+r5(16,10)+r5(16,12)+r5(18,12)+r5(18,22)*three &
&             -r4(12,13)*two-r4(12,35)*six+r3(7,8)+r3(7,10)+r3(7,36)*three+r3(7,40)*three &
&             +r3(9,40)*three-r2(6,34)*six+r1(2,14)*three+r1(2,18)*three
      rxyz(9)=+r7(25,3)+r7(25,4)-r6(18,7)*two-r6(18,8)*two+r5(12,9)+r5(12,10)+r5(12,11) &
&             +r5(14,11)+r5(12,12)+r5(14,12)-r4(9,12)*two-r4(9,13)*two+r3(5,7)+r3(5,8) &
&             +r3(5,9)+r3(5,10)
      rxyz(10)=+r6(25,11)-r5(18,19)*two+r4(12,24)+r4(12,28)+r4(14,28)-r3(9,25)*two &
&             +r2(6,15)+r2(6,19)
      rxyz(11)=+r6(18,11)-r5(12,19)*two+r4(7,24)+r4(7,28)+r4(9,28)-r3(5,25)*two+r2(4,15) &
&             +r2(4,19)
      rxyz(12)=+r7(32,3)-r6(24,7)*two+r5(17,9)+r5(17,11)+r5(19,11)+r5(19,21)-r4(13,12)*two &
&             -r4(13,34)*two+r3(8,7)+r3(8,9)+r3(8,35)+r3(8,39)+r3(10,39)-r2(3,33)*two &
&             +r1(3,13)+r1(3,17)
      rxyz(13)=+r7(32,4)-r6(24,8)*two+r5(17,10)+r5(17,12)+r5(19,12)+r5(19,22) &
&             -r4(13,13)*two-r4(13,35)*two+r3(8,8)+r3(8,10)+r3(8,36)+r3(8,40)+r3(10,40) &
&             -r2(3,34)*two+r1(3,14)+r1(3,18)
      rxyz(14)=+r7(24,4)-r6(17,8)*two+r5(11,10)+r5(11,12)+r5(13,12)+r5(13,22)-r4(8,13)*two &
&             -r4(8,35)*two+r3(4,8)+r3(4,10)+r3(4,36)+r3(4,40)+r3(6,40)-r2(5,34)*two &
&             +r1(1,14)+r1(1,18)
      rxyz(15)=+r7(24,3)-r6(17,7)*two+r5(11,9)+r5(11,11)+r5(13,11)+r5(13,21)-r4(8,12)*two &
&             -r4(8,34)*two+r3(4,7)+r3(4,9)+r3(4,35)+r3(4,39)+r3(6,39)-r2(5,33)*two &
&             +r1(1,13)+r1(1,17)
      rxyz(16)=+r7(18,4)-r6(12,8)*two+r5(7,10)+r5(7,12)+r5(9,12)+r5(18,22)-r4(5,13)*two &
&             -r4(12,35)*two+r3(2,8)+r3(2,10)+r3(7,36)+r3(7,40)+r3(9,40)-r2(6,34)*two &
&             +r1(2,14)+r1(2,18)
      rxyz(17)=+r7(18,3)-r6(12,7)*two+r5(7,9)+r5(7,11)+r5(9,11)+r5(18,21)-r4(5,12)*two &
&             -r4(12,34)*two+r3(2,7)+r3(2,9)+r3(7,35)+r3(7,39)+r3(9,39)-r2(6,33)*two &
&             +r1(2,13)+r1(2,17)
      rxyz(18)=+r7(33,3)-r6(25,7)*two+r5(18,9)+r5(18,11)+r5(20,11)+r5(18,21)-r4(14,12)*two &
&             -r4(12,34)*two+r3(9,7)+r3(9,9)+r3(7,35)+r3(7,39)+r3(9,39)-r2(6,33)*two &
&             +r1(2,13)+r1(2,17)
      rxyz(19)=+r7(33,4)-r6(25,8)*two+r5(18,10)+r5(18,12)+r5(20,12)+r5(18,22) &
&             -r4(14,13)*two-r4(12,35)*two+r3(9,8)+r3(9,10)+r3(7,36)+r3(7,40)+r3(9,40) &
&             -r2(6,34)*two+r1(2,14)+r1(2,18)
      rxyz(20)=+r6(19,11)*four-r5(13,19)*eight+r4(8,24)*four+r4(8,28)*four+r4(10,28)*four &
&             -r3(6,25)*eight+r2(5,15)*four+r2(5,19)*four
      eri(1,1,2,3)=r400+(+r7(13,3)*two+r7(13,4)*two-r6(8,7)*four-r6(8,8)*four+r5(4,9)*two &
&                  +r5(4,10)*two+r5(4,11)*two+r5(6,11)*two+r5(4,12)*two+r5(6,12)*two &
&                  +r5(13,21)*six+r5(13,22)*six-r4(3,12)*four-r4(3,13)*four-r4(8,34)*p12 &
&                  -r4(8,35)*p12+r3(1,7)*two+r3(1,8)*two+r3(1,9)*two+r3(1,10)*two &
&                  +r3(4,35)*six+r3(4,36)*six+r3(4,39)*six+r3(6,39)*six+r3(4,40)*six &
&                  +r3(6,40)*six-r2(5,33)*p12-r2(5,34)*p12+r1(1,13)*six+r1(1,14)*six &
&                  +r1(1,17)*six+r1(1,18)*six)*qx+(+r6(13,10)+r6(13,11)*four+r6(13,12) &
&                  -r5(8,18)*two-r5(8,19)*eight-r5(8,20)*two+r4(4,23)+r4(4,24)*four+r4(4,25) &
&                  +r4(4,27)+r4(6,27)+r4(4,28)*four+r4(6,28)*four+r4(4,29)+r4(6,29)+r4(13,39) &
&                  +r4(13,40)*four+r4(13,41)-r3(3,24)*two-r3(3,25)*eight-r3(3,26)*two &
&                  -r3(8,49)*two-r3(8,50)*eight-r3(8,51)*two+r2(1,14)+r2(1,15)*four+r2(1,16) &
&                  +r2(1,18)+r2(1,19)*four+r2(1,20)+r2(2,48)+r2(2,49)*four+r2(2,50)+r2(2,53) &
&                  +r2(3,53)+r2(2,54)*four+r2(3,54)*four+r2(2,55)+r2(3,55)-r1(3,37)*two &
&                  -r1(3,38)*eight-r1(3,39)*two+r0(17)+r0(18)*four+r0(19)+r0(22)+r0(23)*four &
&                  +r0(24))*xx+(+r5(13,23)*two+r5(13,24)*two-r4(8,36)*four-r4(8,37)*four &
&                  +r3(4,37)*two+r3(4,38)*two+r3(4,41)*two+r3(6,41)*two+r3(4,42)*two &
&                  +r3(6,42)*two-r2(5,35)*four-r2(5,36)*four+r1(1,15)*two+r1(1,16)*two &
&                  +r1(1,19)*two+r1(1,20)*two)*xxx+rxyz(1)*xxxx
      eri(2,1,2,3)=r220+(+r7(24,4)*two-r6(17,8)*four+r5(11,10)*two+r5(11,12)*two &
&                  +r5(13,12)*two+r5(13,22)*two-r4(8,13)*four-r4(8,35)*four+r3(4,8)*two &
&                  +r3(4,10)*two+r3(4,36)*two+r3(4,40)*two+r3(6,40)*two-r2(5,34)*four &
&                  +r1(1,14)*two+r1(1,18)*two)*qx+rxyz(5)*xx
      eri(3,1,2,3)=r202+(+r7(26,4)*two-r6(19,8)*four+r5(13,10)*two+r5(13,12)*two &
&                  +r5(15,12)*two+r5(13,22)*two-r4(10,13)*four-r4(8,35)*four+r3(6,8)*two &
&                  +r3(6,10)*two+r3(4,36)*two+r3(4,40)*two+r3(6,40)*two-r2(5,34)*four &
&                  +r1(1,14)*two+r1(1,18)*two)*qx+(+r7(19,3)*two-r6(13,7)*four+r5(8,9)*two &
&                  +r5(8,11)*two+r5(10,11)*two+r5(19,21)*two-r4(6,12)*four-r4(13,34)*four &
&                  +r3(3,7)*two+r3(3,9)*two+r3(8,35)*two+r3(8,39)*two+r3(10,39)*two &
&                  -r2(3,33)*four+r1(3,13)*two+r1(3,17)*two)*qz+(+r6(26,12)-r5(19,20)*two &
&                  +r4(13,25)+r4(13,29)+r4(15,29)+r4(13,41)-r3(10,26)*two-r3(8,51)*two &
&                  +r2(3,16)+r2(3,20)+r2(2,50)+r2(2,55)+r2(3,55)-r1(3,39)*two+r0(19)+r0(24)) &
&                  *xx+rxyz(20)*xz+(+r6(13,10)-r5(8,18)*two+r4(4,23)+r4(4,27)+r4(6,27) &
&                  +r4(13,39)-r3(3,24)*two-r3(8,49)*two+r2(1,14)+r2(1,18)+r2(2,48)+r2(2,53) &
&                  +r2(3,53)-r1(3,37)*two+r0(17)+r0(22))*zz+(+r5(19,24)*two-r4(13,37)*four &
&                  +r3(8,38)*two+r3(8,42)*two+r3(10,42)*two-r2(3,36)*four+r1(3,16)*two &
&                  +r1(3,20)*two)*xxz+(+r5(13,23)*two-r4(8,36)*four+r3(4,37)*two+r3(4,41)*two &
&                  +r3(6,41)*two-r2(5,35)*four+r1(1,15)*two+r1(1,19)*two)*xzz+rxyz(1)*xxzz
      eri(4,1,2,3)=r310+(+r7(18,3)+r7(18,4)*two-r6(12,7)*two-r6(12,8)*four+r5(7,9) &
&                  +r5(7,10)*two+r5(7,11)+r5(9,11)+r5(7,12)*two+r5(9,12)*two+r5(18,21) &
&                  +r5(18,22)*two-r4(5,12)*two-r4(5,13)*four-r4(12,34)*two-r4(12,35)*four &
&                  +r3(2,7)+r3(2,8)*two+r3(2,9)+r3(2,10)*two+r3(7,35)+r3(7,36)*two+r3(7,39) &
&                  +r3(9,39)+r3(7,40)*two+r3(9,40)*two-r2(6,33)*two-r2(6,34)*four+r1(2,13) &
&                  +r1(2,14)*two+r1(2,17)+r1(2,18)*two)*qx+(+r6(18,11)*two+r6(18,12) &
&                  -r5(12,19)*four-r5(12,20)*two+r4(7,24)*two+r4(7,25)+r4(7,28)*two &
&                  +r4(9,28)*two+r4(7,29)+r4(9,29)-r3(5,25)*four-r3(5,26)*two+r2(4,15)*two &
&                  +r2(4,16)+r2(4,19)*two+r2(4,20))*xx+rxyz(2)*xxx
      eri(5,1,2,3)=r301+(+r7(19,3)+r7(19,4)*two-r6(13,7)*two-r6(13,8)*four+r5(8,9) &
&                  +r5(8,10)*two+r5(8,11)+r5(10,11)+r5(8,12)*two+r5(10,12)*two+r5(19,21) &
&                  +r5(19,22)*two-r4(6,12)*two-r4(6,13)*four-r4(13,34)*two-r4(13,35)*four &
&                  +r3(3,7)+r3(3,8)*two+r3(3,9)+r3(3,10)*two+r3(8,35)+r3(8,36)*two+r3(8,39) &
&                  +r3(10,39)+r3(8,40)*two+r3(10,40)*two-r2(3,33)*two-r2(3,34)*four+r1(3,13) &
&                  +r1(3,14)*two+r1(3,17)+r1(3,18)*two)*qx+(+r7(13,3)-r6(8,7)*two+r5(4,9) &
&                  +r5(4,11)+r5(6,11)+r5(13,21)*three-r4(3,12)*two-r4(8,34)*six+r3(1,7) &
&                  +r3(1,9)+r3(4,35)*three+r3(4,39)*three+r3(6,39)*three-r2(5,33)*six &
&                  +r1(1,13)*three+r1(1,17)*three)*qz+(+r6(19,11)*two+r6(19,12) &
&                  -r5(13,19)*four-r5(13,20)*two+r4(8,24)*two+r4(8,25)+r4(8,28)*two &
&                  +r4(10,28)*two+r4(8,29)+r4(10,29)-r3(6,25)*four-r3(6,26)*two+r2(5,15)*two &
&                  +r2(5,16)+r2(5,19)*two+r2(5,20))*xx+(+r6(13,10)+r6(13,11)*two-r5(8,18)*two &
&                  -r5(8,19)*four+r4(4,23)+r4(4,24)*two+r4(4,27)+r4(6,27)+r4(4,28)*two &
&                  +r4(6,28)*two+r4(13,39)+r4(13,40)*two-r3(3,24)*two-r3(3,25)*four &
&                  -r3(8,49)*two-r3(8,50)*four+r2(1,14)+r2(1,15)*two+r2(1,18)+r2(1,19)*two &
&                  +r2(2,48)+r2(2,49)*two+r2(2,53)+r2(3,53)+r2(2,54)*two+r2(3,54)*two &
&                  -r1(3,37)*two-r1(3,38)*four+r0(17)+r0(18)*two+r0(22)+r0(23)*two)*xz+( &
&                  +r5(19,24)-r4(13,37)*two+r3(8,38)+r3(8,42)+r3(10,42)-r2(3,36)*two+r1(3,16) &
&                  +r1(3,20))*xxx+(+r5(13,23)*two+r5(13,24)-r4(8,36)*four-r4(8,37)*two &
&                  +r3(4,37)*two+r3(4,38)+r3(4,41)*two+r3(6,41)*two+r3(4,42)+r3(6,42) &
&                  -r2(5,35)*four-r2(5,36)*two+r1(1,15)*two+r1(1,16)+r1(1,19)*two+r1(1,20)) &
&                  *xxz+rxyz(1)*xxxz
      eri(6,1,2,3)=r211+(+r7(25,4)*two-r6(18,8)*four+r5(12,10)*two+r5(12,12)*two &
&                  +r5(14,12)*two-r4(9,13)*four+r3(5,8)*two+r3(5,10)*two)*qx+rxyz(17)*qz+( &
&                  +r6(25,12)-r5(18,20)*two+r4(12,25)+r4(12,29)+r4(14,29)-r3(9,26)*two &
&                  +r2(6,16)+r2(6,20))*xx+(+r6(18,11)*two-r5(12,19)*four+r4(7,24)*two &
&                  +r4(7,28)*two+r4(9,28)*two-r3(5,25)*four+r2(4,15)*two+r2(4,19)*two)*xz &
&                  +rxyz(2)*xxz
      eri(1,2,2,3)=r220+(+r7(24,3)*two-r6(17,7)*four+r5(11,9)*two+r5(11,11)*two &
&                  +r5(13,11)*two+r5(13,21)*two-r4(8,12)*four-r4(8,34)*four+r3(4,7)*two &
&                  +r3(4,9)*two+r3(4,35)*two+r3(4,39)*two+r3(6,39)*two-r2(5,33)*four &
&                  +r1(1,13)*two+r1(1,17)*two)*qx+rxyz(4)*xx
      eri(2,2,2,3)=r040
      eri(3,2,2,3)=r022+(+r7(32,3)*two-r6(24,7)*four+r5(17,9)*two+r5(17,11)*two &
&                  +r5(19,11)*two+r5(19,21)*two-r4(13,12)*four-r4(13,34)*four+r3(8,7)*two &
&                  +r3(8,9)*two+r3(8,35)*two+r3(8,39)*two+r3(10,39)*two-r2(3,33)*four &
&                  +r1(3,13)*two+r1(3,17)*two)*qz+rxyz(4)*zz
      eri(4,2,2,3)=r130+rxyz(7)*qx
      eri(5,2,2,3)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,2,3)=r031+rxyz(7)*qz
      eri(1,3,2,3)=r202+(+r7(26,3)*two-r6(19,7)*four+r5(13,9)*two+r5(13,11)*two &
&                  +r5(15,11)*two+r5(13,21)*two-r4(10,12)*four-r4(8,34)*four+r3(6,7)*two &
&                  +r3(6,9)*two+r3(4,35)*two+r3(4,39)*two+r3(6,39)*two-r2(5,33)*four &
&                  +r1(1,13)*two+r1(1,17)*two)*qx+(+r7(19,4)*two-r6(13,8)*four+r5(8,10)*two &
&                  +r5(8,12)*two+r5(10,12)*two+r5(19,22)*two-r4(6,13)*four-r4(13,35)*four &
&                  +r3(3,8)*two+r3(3,10)*two+r3(8,36)*two+r3(8,40)*two+r3(10,40)*two &
&                  -r2(3,34)*four+r1(3,14)*two+r1(3,18)*two)*qz+(+r6(26,10)-r5(19,18)*two &
&                  +r4(13,23)+r4(13,27)+r4(15,27)+r4(13,39)-r3(10,24)*two-r3(8,49)*two &
&                  +r2(3,14)+r2(3,18)+r2(2,48)+r2(2,53)+r2(3,53)-r1(3,37)*two+r0(17)+r0(22)) &
&                  *xx+rxyz(20)*xz+(+r6(13,12)-r5(8,20)*two+r4(4,25)+r4(4,29)+r4(6,29) &
&                  +r4(13,41)-r3(3,26)*two-r3(8,51)*two+r2(1,16)+r2(1,20)+r2(2,50)+r2(2,55) &
&                  +r2(3,55)-r1(3,39)*two+r0(19)+r0(24))*zz+(+r5(19,23)*two-r4(13,36)*four &
&                  +r3(8,37)*two+r3(8,41)*two+r3(10,41)*two-r2(3,35)*four+r1(3,15)*two &
&                  +r1(3,19)*two)*xxz+(+r5(13,24)*two-r4(8,37)*four+r3(4,38)*two+r3(4,42)*two &
&                  +r3(6,42)*two-r2(5,36)*four+r1(1,16)*two+r1(1,20)*two)*xzz+rxyz(1)*xxzz
      eri(2,3,2,3)=r022+(+r7(32,4)*two-r6(24,8)*four+r5(17,10)*two+r5(17,12)*two &
&                  +r5(19,12)*two+r5(19,22)*two-r4(13,13)*four-r4(13,35)*four+r3(8,8)*two &
&                  +r3(8,10)*two+r3(8,36)*two+r3(8,40)*two+r3(10,40)*two-r2(3,34)*four &
&                  +r1(3,14)*two+r1(3,18)*two)*qz+rxyz(5)*zz
      eri(3,3,2,3)=r004+(+r7(34,3)*two+r7(34,4)*two-r6(26,7)*four-r6(26,8)*four &
&                  +r5(19,9)*two+r5(19,10)*two+r5(19,11)*two+r5(21,11)*two+r5(19,12)*two &
&                  +r5(21,12)*two+r5(19,21)*six+r5(19,22)*six-r4(15,12)*four-r4(15,13)*four &
&                  -r4(13,34)*p12-r4(13,35)*p12+r3(10,7)*two+r3(10,8)*two+r3(10,9)*two &
&                  +r3(10,10)*two+r3(8,35)*six+r3(8,36)*six+r3(8,39)*six+r3(10,39)*six &
&                  +r3(8,40)*six+r3(10,40)*six-r2(3,33)*p12-r2(3,34)*p12+r1(3,13)*six &
&                  +r1(3,14)*six+r1(3,17)*six+r1(3,18)*six)*qz+(+r6(26,10)+r6(26,11)*four &
&                  +r6(26,12)-r5(19,18)*two-r5(19,19)*eight-r5(19,20)*two+r4(13,23) &
&                  +r4(13,24)*four+r4(13,25)+r4(13,27)+r4(15,27)+r4(13,28)*four &
&                  +r4(15,28)*four+r4(13,29)+r4(15,29)+r4(13,39)+r4(13,40)*four+r4(13,41) &
&                  -r3(10,24)*two-r3(10,25)*eight-r3(10,26)*two-r3(8,49)*two-r3(8,50)*eight &
&                  -r3(8,51)*two+r2(3,14)+r2(3,15)*four+r2(3,16)+r2(3,18)+r2(3,19)*four &
&                  +r2(3,20)+r2(2,48)+r2(2,49)*four+r2(2,50)+r2(2,53)+r2(3,53)+r2(2,54)*four &
&                  +r2(3,54)*four+r2(2,55)+r2(3,55)-r1(3,37)*two-r1(3,38)*eight-r1(3,39)*two &
&                  +r0(17)+r0(18)*four+r0(19)+r0(22)+r0(23)*four+r0(24))*zz+(+r5(19,23)*two &
&                  +r5(19,24)*two-r4(13,36)*four-r4(13,37)*four+r3(8,37)*two+r3(8,38)*two &
&                  +r3(8,41)*two+r3(10,41)*two+r3(8,42)*two+r3(10,42)*two-r2(3,35)*four &
&                  -r2(3,36)*four+r1(3,15)*two+r1(3,16)*two+r1(3,19)*two+r1(3,20)*two)*zzz &
&                  +rxyz(1)*zzzz
      eri(4,3,2,3)=r112+rxyz(18)*qx+(+r7(25,4)*two-r6(18,8)*four+r5(12,10)*two &
&                  +r5(12,12)*two+r5(14,12)*two-r4(9,13)*four+r3(5,8)*two+r3(5,10)*two)*qz+( &
&                  +r6(25,11)*two-r5(18,19)*four+r4(12,24)*two+r4(12,28)*two+r4(14,28)*two &
&                  -r3(9,25)*four+r2(6,15)*two+r2(6,19)*two)*xz+(+r6(18,12)-r5(12,20)*two &
&                  +r4(7,25)+r4(7,29)+r4(9,29)-r3(5,26)*two+r2(4,16)+r2(4,20))*zz+rxyz(2)*xzz
      eri(5,3,2,3)=r103+(+r7(34,3)-r6(26,7)*two+r5(19,9)+r5(19,11)+r5(21,11) &
&                  +r5(19,21)*three-r4(15,12)*two-r4(13,34)*six+r3(10,7)+r3(10,9) &
&                  +r3(8,35)*three+r3(8,39)*three+r3(10,39)*three-r2(3,33)*six+r1(3,13)*three &
&                  +r1(3,17)*three)*qx+(+r7(26,3)+r7(26,4)*two-r6(19,7)*two-r6(19,8)*four &
&                  +r5(13,9)+r5(13,10)*two+r5(13,11)+r5(15,11)+r5(13,12)*two+r5(15,12)*two &
&                  +r5(13,21)+r5(13,22)*two-r4(10,12)*two-r4(10,13)*four-r4(8,34)*two &
&                  -r4(8,35)*four+r3(6,7)+r3(6,8)*two+r3(6,9)+r3(6,10)*two+r3(4,35) &
&                  +r3(4,36)*two+r3(4,39)+r3(6,39)+r3(4,40)*two+r3(6,40)*two-r2(5,33)*two &
&                  -r2(5,34)*four+r1(1,13)+r1(1,14)*two+r1(1,17)+r1(1,18)*two)*qz+(+r6(26,10) &
&                  +r6(26,11)*two-r5(19,18)*two-r5(19,19)*four+r4(13,23)+r4(13,24)*two &
&                  +r4(13,27)+r4(15,27)+r4(13,28)*two+r4(15,28)*two+r4(13,39)+r4(13,40)*two &
&                  -r3(10,24)*two-r3(10,25)*four-r3(8,49)*two-r3(8,50)*four+r2(3,14) &
&                  +r2(3,15)*two+r2(3,18)+r2(3,19)*two+r2(2,48)+r2(2,49)*two+r2(2,53) &
&                  +r2(3,53)+r2(2,54)*two+r2(3,54)*two-r1(3,37)*two-r1(3,38)*four+r0(17) &
&                  +r0(18)*two+r0(22)+r0(23)*two)*xz+(+r6(19,11)*two+r6(19,12)-r5(13,19)*four &
&                  -r5(13,20)*two+r4(8,24)*two+r4(8,25)+r4(8,28)*two+r4(10,28)*two+r4(8,29) &
&                  +r4(10,29)-r3(6,25)*four-r3(6,26)*two+r2(5,15)*two+r2(5,16)+r2(5,19)*two &
&                  +r2(5,20))*zz+(+r5(19,23)*two+r5(19,24)-r4(13,36)*four-r4(13,37)*two &
&                  +r3(8,37)*two+r3(8,38)+r3(8,41)*two+r3(10,41)*two+r3(8,42)+r3(10,42) &
&                  -r2(3,35)*four-r2(3,36)*two+r1(3,15)*two+r1(3,16)+r1(3,19)*two+r1(3,20)) &
&                  *xzz+(+r5(13,24)-r4(8,37)*two+r3(4,38)+r3(4,42)+r3(6,42)-r2(5,36)*two &
&                  +r1(1,16)+r1(1,20))*zzz+rxyz(1)*xzzz
      eri(6,3,2,3)=r013+(+r7(33,3)+r7(33,4)*two-r6(25,7)*two-r6(25,8)*four+r5(18,9) &
&                  +r5(18,10)*two+r5(18,11)+r5(20,11)+r5(18,12)*two+r5(20,12)*two+r5(18,21) &
&                  +r5(18,22)*two-r4(14,12)*two-r4(14,13)*four-r4(12,34)*two-r4(12,35)*four &
&                  +r3(9,7)+r3(9,8)*two+r3(9,9)+r3(9,10)*two+r3(7,35)+r3(7,36)*two+r3(7,39) &
&                  +r3(9,39)+r3(7,40)*two+r3(9,40)*two-r2(6,33)*two-r2(6,34)*four+r1(2,13) &
&                  +r1(2,14)*two+r1(2,17)+r1(2,18)*two)*qz+(+r6(25,11)*two+r6(25,12) &
&                  -r5(18,19)*four-r5(18,20)*two+r4(12,24)*two+r4(12,25)+r4(12,28)*two &
&                  +r4(14,28)*two+r4(12,29)+r4(14,29)-r3(9,25)*four-r3(9,26)*two+r2(6,15)*two &
&                  +r2(6,16)+r2(6,19)*two+r2(6,20))*zz+rxyz(2)*zzz
      eri(1,4,2,3)=r310+(+r7(18,3)*two+r7(18,4)-r6(12,7)*four-r6(12,8)*two+r5(7,9)*two &
&                  +r5(7,10)+r5(7,11)*two+r5(9,11)*two+r5(7,12)+r5(9,12)+r5(18,21)*two &
&                  +r5(18,22)-r4(5,12)*four-r4(5,13)*two-r4(12,34)*four-r4(12,35)*two &
&                  +r3(2,7)*two+r3(2,8)+r3(2,9)*two+r3(2,10)+r3(7,35)*two+r3(7,36) &
&                  +r3(7,39)*two+r3(9,39)*two+r3(7,40)+r3(9,40)-r2(6,33)*four-r2(6,34)*two &
&                  +r1(2,13)*two+r1(2,14)+r1(2,17)*two+r1(2,18))*qx+(+r6(18,10)+r6(18,11)*two &
&                  -r5(12,18)*two-r5(12,19)*four+r4(7,23)+r4(7,24)*two+r4(7,27)+r4(9,27) &
&                  +r4(7,28)*two+r4(9,28)*two-r3(5,24)*two-r3(5,25)*four+r2(4,14) &
&                  +r2(4,15)*two+r2(4,18)+r2(4,19)*two)*xx+rxyz(3)*xxx
      eri(2,4,2,3)=r130+rxyz(8)*qx
      eri(3,4,2,3)=r112+rxyz(19)*qx+(+r7(25,3)*two-r6(18,7)*four+r5(12,9)*two &
&                  +r5(12,11)*two+r5(14,11)*two-r4(9,12)*four+r3(5,7)*two+r3(5,9)*two)*qz+( &
&                  +r6(25,11)*two-r5(18,19)*four+r4(12,24)*two+r4(12,28)*two+r4(14,28)*two &
&                  -r3(9,25)*four+r2(6,15)*two+r2(6,19)*two)*xz+(+r6(18,10)-r5(12,18)*two &
&                  +r4(7,23)+r4(7,27)+r4(9,27)-r3(5,24)*two+r2(4,14)+r2(4,18))*zz+rxyz(3)*xzz
      eri(4,4,2,3)=r220+(+r7(24,3)+r7(24,4)-r6(17,7)*two-r6(17,8)*two+r5(11,9)+r5(11,10) &
&                  +r5(11,11)+r5(13,11)+r5(11,12)+r5(13,12)+r5(13,21)+r5(13,22)-r4(8,12)*two &
&                  -r4(8,13)*two-r4(8,34)*two-r4(8,35)*two+r3(4,7)+r3(4,8)+r3(4,9)+r3(4,10) &
&                  +r3(4,35)+r3(4,36)+r3(4,39)+r3(6,39)+r3(4,40)+r3(6,40)-r2(5,33)*two &
&                  -r2(5,34)*two+r1(1,13)+r1(1,14)+r1(1,17)+r1(1,18))*qx+rxyz(6)*xx
      eri(5,4,2,3)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(18,10)+r6(18,11) &
&                  -r5(12,18)*two-r5(12,19)*two+r4(7,23)+r4(7,24)+r4(7,27)+r4(9,27)+r4(7,28) &
&                  +r4(9,28)-r3(5,24)*two-r3(5,25)*two+r2(4,14)+r2(4,15)+r2(4,18)+r2(4,19)) &
&                  *xz+rxyz(3)*xxz
      eri(6,4,2,3)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,2,3)=r301+(+r7(19,3)*two+r7(19,4)-r6(13,7)*four-r6(13,8)*two+r5(8,9)*two &
&                  +r5(8,10)+r5(8,11)*two+r5(10,11)*two+r5(8,12)+r5(10,12)+r5(19,21)*two &
&                  +r5(19,22)-r4(6,12)*four-r4(6,13)*two-r4(13,34)*four-r4(13,35)*two &
&                  +r3(3,7)*two+r3(3,8)+r3(3,9)*two+r3(3,10)+r3(8,35)*two+r3(8,36) &
&                  +r3(8,39)*two+r3(10,39)*two+r3(8,40)+r3(10,40)-r2(3,33)*four-r2(3,34)*two &
&                  +r1(3,13)*two+r1(3,14)+r1(3,17)*two+r1(3,18))*qx+(+r7(13,4)-r6(8,8)*two &
&                  +r5(4,10)+r5(4,12)+r5(6,12)+r5(13,22)*three-r4(3,13)*two-r4(8,35)*six &
&                  +r3(1,8)+r3(1,10)+r3(4,36)*three+r3(4,40)*three+r3(6,40)*three &
&                  -r2(5,34)*six+r1(1,14)*three+r1(1,18)*three)*qz+(+r6(19,10)+r6(19,11)*two &
&                  -r5(13,18)*two-r5(13,19)*four+r4(8,23)+r4(8,24)*two+r4(8,27)+r4(10,27) &
&                  +r4(8,28)*two+r4(10,28)*two-r3(6,24)*two-r3(6,25)*four+r2(5,14) &
&                  +r2(5,15)*two+r2(5,18)+r2(5,19)*two)*xx+(+r6(13,11)*two+r6(13,12) &
&                  -r5(8,19)*four-r5(8,20)*two+r4(4,24)*two+r4(4,25)+r4(4,28)*two &
&                  +r4(6,28)*two+r4(4,29)+r4(6,29)+r4(13,40)*two+r4(13,41)-r3(3,25)*four &
&                  -r3(3,26)*two-r3(8,50)*four-r3(8,51)*two+r2(1,15)*two+r2(1,16) &
&                  +r2(1,19)*two+r2(1,20)+r2(2,49)*two+r2(2,50)+r2(2,54)*two+r2(3,54)*two &
&                  +r2(2,55)+r2(3,55)-r1(3,38)*four-r1(3,39)*two+r0(18)*two+r0(19)+r0(23)*two &
&                  +r0(24))*xz+(+r5(19,23)-r4(13,36)*two+r3(8,37)+r3(8,41)+r3(10,41) &
&                  -r2(3,35)*two+r1(3,15)+r1(3,19))*xxx+(+r5(13,23)+r5(13,24)*two &
&                  -r4(8,36)*two-r4(8,37)*four+r3(4,37)+r3(4,38)*two+r3(4,41)+r3(6,41) &
&                  +r3(4,42)*two+r3(6,42)*two-r2(5,35)*two-r2(5,36)*four+r1(1,15) &
&                  +r1(1,16)*two+r1(1,19)+r1(1,20)*two)*xxz+rxyz(1)*xxxz
      eri(2,5,2,3)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,2,3)=r103+(+r7(34,4)-r6(26,8)*two+r5(19,10)+r5(19,12)+r5(21,12) &
&                  +r5(19,22)*three-r4(15,13)*two-r4(13,35)*six+r3(10,8)+r3(10,10) &
&                  +r3(8,36)*three+r3(8,40)*three+r3(10,40)*three-r2(3,34)*six+r1(3,14)*three &
&                  +r1(3,18)*three)*qx+(+r7(26,3)*two+r7(26,4)-r6(19,7)*four-r6(19,8)*two &
&                  +r5(13,9)*two+r5(13,10)+r5(13,11)*two+r5(15,11)*two+r5(13,12)+r5(15,12) &
&                  +r5(13,21)*two+r5(13,22)-r4(10,12)*four-r4(10,13)*two-r4(8,34)*four &
&                  -r4(8,35)*two+r3(6,7)*two+r3(6,8)+r3(6,9)*two+r3(6,10)+r3(4,35)*two &
&                  +r3(4,36)+r3(4,39)*two+r3(6,39)*two+r3(4,40)+r3(6,40)-r2(5,33)*four &
&                  -r2(5,34)*two+r1(1,13)*two+r1(1,14)+r1(1,17)*two+r1(1,18))*qz+( &
&                  +r6(26,11)*two+r6(26,12)-r5(19,19)*four-r5(19,20)*two+r4(13,24)*two &
&                  +r4(13,25)+r4(13,28)*two+r4(15,28)*two+r4(13,29)+r4(15,29)+r4(13,40)*two &
&                  +r4(13,41)-r3(10,25)*four-r3(10,26)*two-r3(8,50)*four-r3(8,51)*two &
&                  +r2(3,15)*two+r2(3,16)+r2(3,19)*two+r2(3,20)+r2(2,49)*two+r2(2,50) &
&                  +r2(2,54)*two+r2(3,54)*two+r2(2,55)+r2(3,55)-r1(3,38)*four-r1(3,39)*two &
&                  +r0(18)*two+r0(19)+r0(23)*two+r0(24))*xz+(+r6(19,10)+r6(19,11)*two &
&                  -r5(13,18)*two-r5(13,19)*four+r4(8,23)+r4(8,24)*two+r4(8,27)+r4(10,27) &
&                  +r4(8,28)*two+r4(10,28)*two-r3(6,24)*two-r3(6,25)*four+r2(5,14) &
&                  +r2(5,15)*two+r2(5,18)+r2(5,19)*two)*zz+(+r5(19,23)+r5(19,24)*two &
&                  -r4(13,36)*two-r4(13,37)*four+r3(8,37)+r3(8,38)*two+r3(8,41)+r3(10,41) &
&                  +r3(8,42)*two+r3(10,42)*two-r2(3,35)*two-r2(3,36)*four+r1(3,15) &
&                  +r1(3,16)*two+r1(3,19)+r1(3,20)*two)*xzz+(+r5(13,23)-r4(8,36)*two+r3(4,37) &
&                  +r3(4,41)+r3(6,41)-r2(5,35)*two+r1(1,15)+r1(1,19))*zzz+rxyz(1)*xzzz
      eri(4,5,2,3)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(18,11)+r6(18,12) &
&                  -r5(12,19)*two-r5(12,20)*two+r4(7,24)+r4(7,25)+r4(7,28)+r4(9,28)+r4(7,29) &
&                  +r4(9,29)-r3(5,25)*two-r3(5,26)*two+r2(4,15)+r2(4,16)+r2(4,19)+r2(4,20)) &
&                  *xz+rxyz(2)*xxz
      eri(5,5,2,3)=r202+(+r7(26,3)+r7(26,4)-r6(19,7)*two-r6(19,8)*two+r5(13,9)+r5(13,10) &
&                  +r5(13,11)+r5(15,11)+r5(13,12)+r5(15,12)+r5(13,21)+r5(13,22)-r4(10,12)*two &
&                  -r4(10,13)*two-r4(8,34)*two-r4(8,35)*two+r3(6,7)+r3(6,8)+r3(6,9)+r3(6,10) &
&                  +r3(4,35)+r3(4,36)+r3(4,39)+r3(6,39)+r3(4,40)+r3(6,40)-r2(5,33)*two &
&                  -r2(5,34)*two+r1(1,13)+r1(1,14)+r1(1,17)+r1(1,18))*qx+(+r7(19,3)+r7(19,4) &
&                  -r6(13,7)*two-r6(13,8)*two+r5(8,9)+r5(8,10)+r5(8,11)+r5(10,11)+r5(8,12) &
&                  +r5(10,12)+r5(19,21)+r5(19,22)-r4(6,12)*two-r4(6,13)*two-r4(13,34)*two &
&                  -r4(13,35)*two+r3(3,7)+r3(3,8)+r3(3,9)+r3(3,10)+r3(8,35)+r3(8,36)+r3(8,39) &
&                  +r3(10,39)+r3(8,40)+r3(10,40)-r2(3,33)*two-r2(3,34)*two+r1(3,13)+r1(3,14) &
&                  +r1(3,17)+r1(3,18))*qz+(+r6(26,11)-r5(19,19)*two+r4(13,24)+r4(13,28) &
&                  +r4(15,28)+r4(13,40)-r3(10,25)*two-r3(8,50)*two+r2(3,15)+r2(3,19)+r2(2,49) &
&                  +r2(2,54)+r2(3,54)-r1(3,38)*two+r0(18)+r0(23))*xx+(+r6(19,10) &
&                  +r6(19,11)*two+r6(19,12)-r5(13,18)*two-r5(13,19)*four-r5(13,20)*two &
&                  +r4(8,23)+r4(8,24)*two+r4(8,25)+r4(8,27)+r4(10,27)+r4(8,28)*two &
&                  +r4(10,28)*two+r4(8,29)+r4(10,29)-r3(6,24)*two-r3(6,25)*four-r3(6,26)*two &
&                  +r2(5,14)+r2(5,15)*two+r2(5,16)+r2(5,18)+r2(5,19)*two+r2(5,20))*xz+( &
&                  +r6(13,11)-r5(8,19)*two+r4(4,24)+r4(4,28)+r4(6,28)+r4(13,40)-r3(3,25)*two &
&                  -r3(8,50)*two+r2(1,15)+r2(1,19)+r2(2,49)+r2(2,54)+r2(3,54)-r1(3,38)*two &
&                  +r0(18)+r0(23))*zz+(+r5(19,23)+r5(19,24)-r4(13,36)*two-r4(13,37)*two &
&                  +r3(8,37)+r3(8,38)+r3(8,41)+r3(10,41)+r3(8,42)+r3(10,42)-r2(3,35)*two &
&                  -r2(3,36)*two+r1(3,15)+r1(3,16)+r1(3,19)+r1(3,20))*xxz+(+r5(13,23) &
&                  +r5(13,24)-r4(8,36)*two-r4(8,37)*two+r3(4,37)+r3(4,38)+r3(4,41)+r3(6,41) &
&                  +r3(4,42)+r3(6,42)-r2(5,35)*two-r2(5,36)*two+r1(1,15)+r1(1,16)+r1(1,19) &
&                  +r1(1,20))*xzz+rxyz(1)*xxzz
      eri(6,5,2,3)=r112+rxyz(19)*qx+(+r7(25,3)+r7(25,4)-r6(18,7)*two-r6(18,8)*two &
&                  +r5(12,9)+r5(12,10)+r5(12,11)+r5(14,11)+r5(12,12)+r5(14,12)-r4(9,12)*two &
&                  -r4(9,13)*two+r3(5,7)+r3(5,8)+r3(5,9)+r3(5,10))*qz+(+r6(25,11)+r6(25,12) &
&                  -r5(18,19)*two-r5(18,20)*two+r4(12,24)+r4(12,25)+r4(12,28)+r4(14,28) &
&                  +r4(12,29)+r4(14,29)-r3(9,25)*two-r3(9,26)*two+r2(6,15)+r2(6,16)+r2(6,19) &
&                  +r2(6,20))*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,2,3)=r211+(+r7(25,3)*two-r6(18,7)*four+r5(12,9)*two+r5(12,11)*two &
&                  +r5(14,11)*two-r4(9,12)*four+r3(5,7)*two+r3(5,9)*two)*qx+rxyz(16)*qz+( &
&                  +r6(25,10)-r5(18,18)*two+r4(12,23)+r4(12,27)+r4(14,27)-r3(9,24)*two &
&                  +r2(6,14)+r2(6,18))*xx+(+r6(18,11)*two-r5(12,19)*four+r4(7,24)*two &
&                  +r4(7,28)*two+r4(9,28)*two-r3(5,25)*four+r2(4,15)*two+r2(4,19)*two)*xz &
&                  +rxyz(3)*xxz
      eri(2,6,2,3)=r031+rxyz(8)*qz
      eri(3,6,2,3)=r013+(+r7(33,3)*two+r7(33,4)-r6(25,7)*four-r6(25,8)*two+r5(18,9)*two &
&                  +r5(18,10)+r5(18,11)*two+r5(20,11)*two+r5(18,12)+r5(20,12)+r5(18,21)*two &
&                  +r5(18,22)-r4(14,12)*four-r4(14,13)*two-r4(12,34)*four-r4(12,35)*two &
&                  +r3(9,7)*two+r3(9,8)+r3(9,9)*two+r3(9,10)+r3(7,35)*two+r3(7,36) &
&                  +r3(7,39)*two+r3(9,39)*two+r3(7,40)+r3(9,40)-r2(6,33)*four-r2(6,34)*two &
&                  +r1(2,13)*two+r1(2,14)+r1(2,17)*two+r1(2,18))*qz+(+r6(25,10)+r6(25,11)*two &
&                  -r5(18,18)*two-r5(18,19)*four+r4(12,23)+r4(12,24)*two+r4(12,27)+r4(14,27) &
&                  +r4(12,28)*two+r4(14,28)*two-r3(9,24)*two-r3(9,25)*four+r2(6,14) &
&                  +r2(6,15)*two+r2(6,18)+r2(6,19)*two)*zz+rxyz(3)*zzz
      eri(4,6,2,3)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,2,3)=r112+rxyz(18)*qx+(+r7(25,3)+r7(25,4)-r6(18,7)*two-r6(18,8)*two &
&                  +r5(12,9)+r5(12,10)+r5(12,11)+r5(14,11)+r5(12,12)+r5(14,12)-r4(9,12)*two &
&                  -r4(9,13)*two+r3(5,7)+r3(5,8)+r3(5,9)+r3(5,10))*qz+(+r6(25,10)+r6(25,11) &
&                  -r5(18,18)*two-r5(18,19)*two+r4(12,23)+r4(12,24)+r4(12,27)+r4(14,27) &
&                  +r4(12,28)+r4(14,28)-r3(9,24)*two-r3(9,25)*two+r2(6,14)+r2(6,15)+r2(6,18) &
&                  +r2(6,19))*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,2,3)=r022+(+r7(32,3)+r7(32,4)-r6(24,7)*two-r6(24,8)*two+r5(17,9)+r5(17,10) &
&                  +r5(17,11)+r5(19,11)+r5(17,12)+r5(19,12)+r5(19,21)+r5(19,22)-r4(13,12)*two &
&                  -r4(13,13)*two-r4(13,34)*two-r4(13,35)*two+r3(8,7)+r3(8,8)+r3(8,9) &
&                  +r3(8,10)+r3(8,35)+r3(8,36)+r3(8,39)+r3(10,39)+r3(8,40)+r3(10,40) &
&                  -r2(3,33)*two-r2(3,34)*two+r1(3,13)+r1(3,14)+r1(3,17)+r1(3,18))*qz+rxyz(6) &
&                  *zz
!
      r400= r8(15)-r7(10,1)*two-r7(10,2)*two+r6(6,1)+r6(6,2)*four+r6(6,3)+r6(6,4)*six &
&          +r6(15,9)*six-r5(3,1)*two-r5(3,2)*two-r5(3,3)*six-r5(3,4)*six-r5(10,13)*p12 &
&          -r5(10,17)*p12+r4(1,1)+r4(1,2)+r4(1,3)*four+r4(1,4)+r4(1,5)*three+r4(6,14)*six &
&          +r4(6,18)*p24+r4(6,22)*six+r4(6,26)*p36+r4(15,38)*three-r3(3,11)*p12 &
&          -r3(3,15)*p12-r3(3,19)*p36-r3(3,23)*p36-r3(10,43)*six-r3(10,48)*six+r2(1,1)*six &
&          +r2(1,5)*six+r2(1,9)*p24+r2(1,13)*six+r2(1,17)*p18+r2(3,37)*three+r2(3,42)*p12 &
&          +r2(3,47)*three+r2(3,52)*p18-r1(3,21)*six-r1(3,26)*six-r1(3,31)*p18 &
&          -r1(3,36)*p18+r0(1)*three+r0(6)*three+r0(11)*p12+r0(16)*three+r0(21)*nine
      r310= r8(20)-r7(14,1)*two-r7(14,2)*two+r6(9,1)+r6(9,2)*four+r6(9,3)+r6(9,4)*six &
&          +r6(20,9)*three-r5(5,1)*two-r5(5,2)*two-r5(5,3)*six-r5(5,4)*six-r5(14,13)*six &
&          -r5(14,17)*six+r4(2,1)+r4(2,2)+r4(2,3)*four+r4(2,4)+r4(2,5)*three &
&          +r4(9,14)*three+r4(9,18)*p12+r4(9,22)*three+r4(9,26)*p18-r3(5,11)*six &
&          -r3(5,15)*six-r3(5,19)*p18-r3(5,23)*p18+r2(4,1)*three+r2(4,5)*three+r2(4,9)*p12 &
&          +r2(4,13)*three+r2(4,17)*nine
      r301= r8(21)-r7(15,1)*two-r7(15,2)*two+r6(10,1)+r6(10,2)*four+r6(10,3)+r6(10,4)*six &
&          +r6(21,9)*three-r5(6,1)*two-r5(6,2)*two-r5(6,3)*six-r5(6,4)*six-r5(15,13)*six &
&          -r5(15,17)*six+r4(3,1)+r4(3,2)+r4(3,3)*four+r4(3,4)+r4(3,5)*three &
&          +r4(10,14)*three+r4(10,18)*p12+r4(10,22)*three+r4(10,26)*p18-r3(6,11)*six &
&          -r3(6,15)*six-r3(6,19)*p18-r3(6,23)*p18+r2(5,1)*three+r2(5,5)*three+r2(5,9)*p12 &
&          +r2(5,13)*three+r2(5,17)*nine
      r220= r8(26)-r7(19,1)*two-r7(19,2)*two+r6(13,1)+r6(13,2)*four+r6(13,3)+r6(13,4)*six &
&          +r6(15,9)+r6(26,9)-r5(8,1)*two-r5(8,2)*two-r5(8,3)*six-r5(8,4)*six &
&          -r5(10,13)*two-r5(19,13)*two-r5(10,17)*two-r5(19,17)*two+r4(4,1)+r4(4,2) &
&          +r4(4,3)*four+r4(4,4)+r4(4,5)*three+r4(6,14)+r4(13,14)+r4(6,18)*four &
&          +r4(13,18)*four+r4(6,22)+r4(13,22)+r4(6,26)*six+r4(13,26)*six+r4(15,38) &
&          -r3(3,11)*two-r3(8,11)*two-r3(3,15)*two-r3(8,15)*two-r3(3,19)*six-r3(8,19)*six &
&          -r3(3,23)*six-r3(8,23)*six-r3(10,43)*two-r3(10,48)*two+r2(1,1)+r2(2,1)+r2(1,5) &
&          +r2(2,5)+r2(1,9)*four+r2(2,9)*four+r2(1,13)+r2(2,13)+r2(1,17)*three &
&          +r2(2,17)*three+r2(3,37)+r2(3,42)*four+r2(3,47)+r2(3,52)*six-r1(3,21)*two &
&          -r1(3,26)*two-r1(3,31)*six-r1(3,36)*six+r0(1)+r0(6)+r0(11)*four+r0(16) &
&          +r0(21)*three
      r211= r8(27)-r7(20,1)*two-r7(20,2)*two+r6(14,1)+r6(14,2)*four+r6(14,3)+r6(14,4)*six &
&          +r6(27,9)-r5(9,1)*two-r5(9,2)*two-r5(9,3)*six-r5(9,4)*six-r5(20,13)*two &
&          -r5(20,17)*two+r4(5,1)+r4(5,2)+r4(5,3)*four+r4(5,4)+r4(5,5)*three+r4(14,14) &
&          +r4(14,18)*four+r4(14,22)+r4(14,26)*six-r3(9,11)*two-r3(9,15)*two-r3(9,19)*six &
&          -r3(9,23)*six+r2(6,1)+r2(6,5)+r2(6,9)*four+r2(6,13)+r2(6,17)*three
      r202= r8(28)-r7(21,1)*two-r7(21,2)*two+r6(15,1)+r6(15,2)*four+r6(15,3)+r6(15,4)*six &
&          +r6(15,9)+r6(28,9)-r5(10,1)*two-r5(10,2)*two-r5(10,3)*six-r5(10,4)*six &
&          -r5(10,13)*two-r5(21,13)*two-r5(10,17)*two-r5(21,17)*two+r4(6,1)+r4(6,2) &
&          +r4(6,3)*four+r4(6,4)+r4(6,5)*three+r4(6,14)+r4(15,14)+r4(6,18)*four &
&          +r4(15,18)*four+r4(6,22)+r4(15,22)+r4(6,26)*six+r4(15,26)*six+r4(15,38) &
&          -r3(3,11)*two-r3(10,11)*two-r3(3,15)*two-r3(10,15)*two-r3(3,19)*six &
&          -r3(10,19)*six-r3(3,23)*six-r3(10,23)*six-r3(10,43)*two-r3(10,48)*two+r2(1,1) &
&          +r2(3,1)+r2(1,5)+r2(3,5)+r2(1,9)*four+r2(3,9)*four+r2(1,13)+r2(3,13) &
&          +r2(1,17)*three+r2(3,17)*three+r2(3,37)+r2(3,42)*four+r2(3,47)+r2(3,52)*six &
&          -r1(3,21)*two-r1(3,26)*two-r1(3,31)*six-r1(3,36)*six+r0(1)+r0(6)+r0(11)*four &
&          +r0(16)+r0(21)*three
      r130= r8(33)-r7(25,1)*two-r7(25,2)*two+r6(18,1)+r6(18,2)*four+r6(18,3)+r6(18,4)*six &
&          +r6(20,9)*three-r5(12,1)*two-r5(12,2)*two-r5(12,3)*six-r5(12,4)*six &
&          -r5(14,13)*six-r5(14,17)*six+r4(7,1)+r4(7,2)+r4(7,3)*four+r4(7,4)+r4(7,5)*three &
&          +r4(9,14)*three+r4(9,18)*p12+r4(9,22)*three+r4(9,26)*p18-r3(5,11)*six &
&          -r3(5,15)*six-r3(5,19)*p18-r3(5,23)*p18+r2(4,1)*three+r2(4,5)*three+r2(4,9)*p12 &
&          +r2(4,13)*three+r2(4,17)*nine
      r121= r8(34)-r7(26,1)*two-r7(26,2)*two+r6(19,1)+r6(19,2)*four+r6(19,3)+r6(19,4)*six &
&          +r6(21,9)-r5(13,1)*two-r5(13,2)*two-r5(13,3)*six-r5(13,4)*six-r5(15,13)*two &
&          -r5(15,17)*two+r4(8,1)+r4(8,2)+r4(8,3)*four+r4(8,4)+r4(8,5)*three+r4(10,14) &
&          +r4(10,18)*four+r4(10,22)+r4(10,26)*six-r3(6,11)*two-r3(6,15)*two-r3(6,19)*six &
&          -r3(6,23)*six+r2(5,1)+r2(5,5)+r2(5,9)*four+r2(5,13)+r2(5,17)*three
      r112= r8(35)-r7(27,1)*two-r7(27,2)*two+r6(20,1)+r6(20,2)*four+r6(20,3)+r6(20,4)*six &
&          +r6(20,9)-r5(14,1)*two-r5(14,2)*two-r5(14,3)*six-r5(14,4)*six-r5(14,13)*two &
&          -r5(14,17)*two+r4(9,1)+r4(9,2)+r4(9,3)*four+r4(9,4)+r4(9,5)*three+r4(9,14) &
&          +r4(9,18)*four+r4(9,22)+r4(9,26)*six-r3(5,11)*two-r3(5,15)*two-r3(5,19)*six &
&          -r3(5,23)*six+r2(4,1)+r2(4,5)+r2(4,9)*four+r2(4,13)+r2(4,17)*three
      r103= r8(36)-r7(28,1)*two-r7(28,2)*two+r6(21,1)+r6(21,2)*four+r6(21,3)+r6(21,4)*six &
&          +r6(21,9)*three-r5(15,1)*two-r5(15,2)*two-r5(15,3)*six-r5(15,4)*six &
&          -r5(15,13)*six-r5(15,17)*six+r4(10,1)+r4(10,2)+r4(10,3)*four+r4(10,4) &
&          +r4(10,5)*three+r4(10,14)*three+r4(10,18)*p12+r4(10,22)*three+r4(10,26)*p18 &
&          -r3(6,11)*six-r3(6,15)*six-r3(6,19)*p18-r3(6,23)*p18+r2(5,1)*three &
&          +r2(5,5)*three+r2(5,9)*p12+r2(5,13)*three+r2(5,17)*nine
      r040= r8(41)-r7(32,1)*two-r7(32,2)*two+r6(24,1)+r6(24,2)*four+r6(24,3)+r6(24,4)*six &
&          +r6(26,9)*six-r5(17,1)*two-r5(17,2)*two-r5(17,3)*six-r5(17,4)*six-r5(19,13)*p12 &
&          -r5(19,17)*p12+r4(11,1)+r4(11,2)+r4(11,3)*four+r4(11,4)+r4(11,5)*three &
&          +r4(13,14)*six+r4(13,18)*p24+r4(13,22)*six+r4(13,26)*p36+r4(15,38)*three &
&          -r3(8,11)*p12-r3(8,15)*p12-r3(8,19)*p36-r3(8,23)*p36-r3(10,43)*six &
&          -r3(10,48)*six+r2(2,1)*six+r2(2,5)*six+r2(2,9)*p24+r2(2,13)*six+r2(2,17)*p18 &
&          +r2(3,37)*three+r2(3,42)*p12+r2(3,47)*three+r2(3,52)*p18-r1(3,21)*six &
&          -r1(3,26)*six-r1(3,31)*p18-r1(3,36)*p18+r0(1)*three+r0(6)*three+r0(11)*p12 &
&          +r0(16)*three+r0(21)*nine
      r031= r8(42)-r7(33,1)*two-r7(33,2)*two+r6(25,1)+r6(25,2)*four+r6(25,3)+r6(25,4)*six &
&          +r6(27,9)*three-r5(18,1)*two-r5(18,2)*two-r5(18,3)*six-r5(18,4)*six &
&          -r5(20,13)*six-r5(20,17)*six+r4(12,1)+r4(12,2)+r4(12,3)*four+r4(12,4) &
&          +r4(12,5)*three+r4(14,14)*three+r4(14,18)*p12+r4(14,22)*three+r4(14,26)*p18 &
&          -r3(9,11)*six-r3(9,15)*six-r3(9,19)*p18-r3(9,23)*p18+r2(6,1)*three &
&          +r2(6,5)*three+r2(6,9)*p12+r2(6,13)*three+r2(6,17)*nine
      r022= r8(43)-r7(34,1)*two-r7(34,2)*two+r6(26,1)+r6(26,2)*four+r6(26,3)+r6(26,4)*six &
&          +r6(26,9)+r6(28,9)-r5(19,1)*two-r5(19,2)*two-r5(19,3)*six-r5(19,4)*six &
&          -r5(19,13)*two-r5(21,13)*two-r5(19,17)*two-r5(21,17)*two+r4(13,1)+r4(13,2) &
&          +r4(13,3)*four+r4(13,4)+r4(13,5)*three+r4(13,14)+r4(15,14)+r4(13,18)*four &
&          +r4(15,18)*four+r4(13,22)+r4(15,22)+r4(13,26)*six+r4(15,26)*six+r4(15,38) &
&          -r3(8,11)*two-r3(10,11)*two-r3(8,15)*two-r3(10,15)*two-r3(8,19)*six &
&          -r3(10,19)*six-r3(8,23)*six-r3(10,23)*six-r3(10,43)*two-r3(10,48)*two+r2(2,1) &
&          +r2(3,1)+r2(2,5)+r2(3,5)+r2(2,9)*four+r2(3,9)*four+r2(2,13)+r2(3,13) &
&          +r2(2,17)*three+r2(3,17)*three+r2(3,37)+r2(3,42)*four+r2(3,47)+r2(3,52)*six &
&          -r1(3,21)*two-r1(3,26)*two-r1(3,31)*six-r1(3,36)*six+r0(1)+r0(6)+r0(11)*four &
&          +r0(16)+r0(21)*three
      r013= r8(44)-r7(35,1)*two-r7(35,2)*two+r6(27,1)+r6(27,2)*four+r6(27,3)+r6(27,4)*six &
&          +r6(27,9)*three-r5(20,1)*two-r5(20,2)*two-r5(20,3)*six-r5(20,4)*six &
&          -r5(20,13)*six-r5(20,17)*six+r4(14,1)+r4(14,2)+r4(14,3)*four+r4(14,4) &
&          +r4(14,5)*three+r4(14,14)*three+r4(14,18)*p12+r4(14,22)*three+r4(14,26)*p18 &
&          -r3(9,11)*six-r3(9,15)*six-r3(9,19)*p18-r3(9,23)*p18+r2(6,1)*three &
&          +r2(6,5)*three+r2(6,9)*p12+r2(6,13)*three+r2(6,17)*nine
      r004= r8(45)-r7(36,1)*two-r7(36,2)*two+r6(28,1)+r6(28,2)*four+r6(28,3)+r6(28,4)*six &
&          +r6(28,9)*six-r5(21,1)*two-r5(21,2)*two-r5(21,3)*six-r5(21,4)*six-r5(21,13)*p12 &
&          -r5(21,17)*p12+r4(15,1)+r4(15,2)+r4(15,3)*four+r4(15,4)+r4(15,5)*three &
&          +r4(15,14)*six+r4(15,18)*p24+r4(15,22)*six+r4(15,26)*p36+r4(15,38)*three &
&          -r3(10,11)*p12-r3(10,15)*p12-r3(10,19)*p36-r3(10,23)*p36-r3(10,43)*six &
&          -r3(10,48)*six+r2(3,1)*six+r2(3,5)*six+r2(3,9)*p24+r2(3,13)*six+r2(3,17)*p18 &
&          +r2(3,37)*three+r2(3,42)*p12+r2(3,47)*three+r2(3,52)*p18-r1(3,21)*six &
&          -r1(3,26)*six-r1(3,31)*p18-r1(3,36)*p18+r0(1)*three+r0(6)*three+r0(11)*p12 &
&          +r0(16)*three+r0(21)*nine
      rxyz(1)=+r4(15,42)-r3(10,47)*two-r3(10,52)*two+r2(3,41)+r2(3,46)*four+r2(3,51) &
&             +r2(3,56)*six-r1(3,25)*two-r1(3,30)*two-r1(3,35)*six-r1(3,40)*six+r0(5) &
&             +r0(10)+r0(15)*four+r0(20)+r0(25)*three
      rxyz(2)=+r5(20,24)-r4(14,33)*two-r4(14,37)*two+r3(9,30)+r3(9,34)*four+r3(9,38) &
&             +r3(9,42)*six-r2(6,24)*two-r2(6,28)*two-r2(6,32)*six-r2(6,36)*six+r1(2,4) &
&             +r1(2,8)+r1(2,12)*four+r1(2,16)+r1(2,20)*three
      rxyz(3)=+r5(20,23)-r4(14,32)*two-r4(14,36)*two+r3(9,29)+r3(9,33)*four+r3(9,37) &
&             +r3(9,41)*six-r2(6,23)*two-r2(6,27)*two-r2(6,31)*six-r2(6,35)*six+r1(2,3) &
&             +r1(2,7)+r1(2,11)*four+r1(2,15)+r1(2,19)*three
      rxyz(4)=+r6(26,10)-r5(19,14)*two-r5(19,18)*two+r4(13,15)+r4(13,19)*four+r4(13,23) &
&             +r4(13,27)*six+r4(15,39)-r3(8,12)*two-r3(8,16)*two-r3(8,20)*six-r3(8,24)*six &
&             -r3(10,44)*two-r3(10,49)*two+r2(2,2)+r2(2,6)+r2(2,10)*four+r2(2,14) &
&             +r2(2,18)*three+r2(3,38)+r2(3,43)*four+r2(3,48)+r2(3,53)*six-r1(3,22)*two &
&             -r1(3,27)*two-r1(3,32)*six-r1(3,37)*six+r0(2)+r0(7)+r0(12)*four+r0(17) &
&             +r0(22)*three
      rxyz(5)=+r6(26,12)-r5(19,16)*two-r5(19,20)*two+r4(13,17)+r4(13,21)*four+r4(13,25) &
&             +r4(13,29)*six+r4(15,41)-r3(8,14)*two-r3(8,18)*two-r3(8,22)*six-r3(8,26)*six &
&             -r3(10,46)*two-r3(10,51)*two+r2(2,4)+r2(2,8)+r2(2,12)*four+r2(2,16) &
&             +r2(2,20)*three+r2(3,40)+r2(3,45)*four+r2(3,50)+r2(3,55)*six-r1(3,24)*two &
&             -r1(3,29)*two-r1(3,34)*six-r1(3,39)*six+r0(4)+r0(9)+r0(14)*four+r0(19) &
&             +r0(24)*three
      rxyz(6)=+r6(26,11)-r5(19,15)*two-r5(19,19)*two+r4(13,16)+r4(13,20)*four+r4(13,24) &
&             +r4(13,28)*six+r4(15,40)-r3(8,13)*two-r3(8,17)*two-r3(8,21)*six-r3(8,25)*six &
&             -r3(10,45)*two-r3(10,50)*two+r2(2,3)+r2(2,7)+r2(2,11)*four+r2(2,15) &
&             +r2(2,19)*three+r2(3,39)+r2(3,44)*four+r2(3,49)+r2(3,54)*six-r1(3,23)*two &
&             -r1(3,28)*two-r1(3,33)*six-r1(3,38)*six+r0(3)+r0(8)+r0(13)*four+r0(18) &
&             +r0(23)*three
      rxyz(7)=+r7(33,3)-r6(25,5)*two-r6(25,7)*two+r5(18,5)+r5(18,7)*four+r5(18,9) &
&             +r5(18,11)*six+r5(20,21)*three-r4(12,6)*two-r4(12,8)*two-r4(12,10)*six &
&             -r4(12,12)*six-r4(14,30)*six-r4(14,34)*six+r3(7,1)+r3(7,3)+r3(7,5)*four &
&             +r3(7,7)+r3(7,9)*three+r3(9,27)*three+r3(9,31)*p12+r3(9,35)*three &
&             +r3(9,39)*p18-r2(6,21)*six-r2(6,25)*six-r2(6,29)*p18-r2(6,33)*p18 &
&             +r1(2,1)*three+r1(2,5)*three+r1(2,9)*p12+r1(2,13)*three+r1(2,17)*nine
      rxyz(8)=+r7(33,4)-r6(25,6)*two-r6(25,8)*two+r5(18,6)+r5(18,8)*four+r5(18,10) &
&             +r5(18,12)*six+r5(20,22)*three-r4(12,7)*two-r4(12,9)*two-r4(12,11)*six &
&             -r4(12,13)*six-r4(14,31)*six-r4(14,35)*six+r3(7,2)+r3(7,4)+r3(7,6)*four &
&             +r3(7,8)+r3(7,10)*three+r3(9,28)*three+r3(9,32)*p12+r3(9,36)*three &
&             +r3(9,40)*p18-r2(6,22)*six-r2(6,26)*six-r2(6,30)*p18-r2(6,34)*p18 &
&             +r1(2,2)*three+r1(2,6)*three+r1(2,10)*p12+r1(2,14)*three+r1(2,18)*nine
      rxyz(9)=+r7(27,3)+r7(27,4)-r6(20,5)*two-r6(20,6)*two-r6(20,7)*two-r6(20,8)*two &
&             +r5(14,5)+r5(14,6)+r5(14,7)*four+r5(14,8)*four+r5(14,9)+r5(14,10) &
&             +r5(14,11)*six+r5(14,12)*six-r4(9,6)*two-r4(9,7)*two-r4(9,8)*two-r4(9,9)*two &
&             -r4(9,10)*six-r4(9,11)*six-r4(9,12)*six-r4(9,13)*six+r3(5,1)+r3(5,2)+r3(5,3) &
&             +r3(5,4)+r3(5,5)*four+r3(5,6)*four+r3(5,7)+r3(5,8)+r3(5,9)*three &
&             +r3(5,10)*three
      rxyz(10)=+r6(27,11)-r5(20,15)*two-r5(20,19)*two+r4(14,16)+r4(14,20)*four+r4(14,24) &
&             +r4(14,28)*six-r3(9,13)*two-r3(9,17)*two-r3(9,21)*six-r3(9,25)*six+r2(6,3) &
&             +r2(6,7)+r2(6,11)*four+r2(6,15)+r2(6,19)*three
      rxyz(11)=+r6(20,11)-r5(14,15)*two-r5(14,19)*two+r4(9,16)+r4(9,20)*four+r4(9,24) &
&             +r4(9,28)*six-r3(5,13)*two-r3(5,17)*two-r3(5,21)*six-r3(5,25)*six+r2(4,3) &
&             +r2(4,7)+r2(4,11)*four+r2(4,15)+r2(4,19)*three
      rxyz(12)=+r7(34,3)-r6(26,5)*two-r6(26,7)*two+r5(19,5)+r5(19,7)*four+r5(19,9) &
&             +r5(19,11)*six+r5(21,21)-r4(13,6)*two-r4(13,8)*two-r4(13,10)*six &
&             -r4(13,12)*six-r4(15,30)*two-r4(15,34)*two+r3(8,1)+r3(8,3)+r3(8,5)*four &
&             +r3(8,7)+r3(8,9)*three+r3(10,27)+r3(10,31)*four+r3(10,35)+r3(10,39)*six &
&             -r2(3,21)*two-r2(3,25)*two-r2(3,29)*six-r2(3,33)*six+r1(3,1)+r1(3,5) &
&             +r1(3,9)*four+r1(3,13)+r1(3,17)*three
      rxyz(13)=+r7(34,4)-r6(26,6)*two-r6(26,8)*two+r5(19,6)+r5(19,8)*four+r5(19,10) &
&             +r5(19,12)*six+r5(21,22)-r4(13,7)*two-r4(13,9)*two-r4(13,11)*six &
&             -r4(13,13)*six-r4(15,31)*two-r4(15,35)*two+r3(8,2)+r3(8,4)+r3(8,6)*four &
&             +r3(8,8)+r3(8,10)*three+r3(10,28)+r3(10,32)*four+r3(10,36)+r3(10,40)*six &
&             -r2(3,22)*two-r2(3,26)*two-r2(3,30)*six-r2(3,34)*six+r1(3,2)+r1(3,6) &
&             +r1(3,10)*four+r1(3,14)+r1(3,18)*three
      rxyz(14)=+r7(26,4)-r6(19,6)*two-r6(19,8)*two+r5(13,6)+r5(13,8)*four+r5(13,10) &
&             +r5(13,12)*six+r5(15,22)-r4(8,7)*two-r4(8,9)*two-r4(8,11)*six-r4(8,13)*six &
&             -r4(10,31)*two-r4(10,35)*two+r3(4,2)+r3(4,4)+r3(4,6)*four+r3(4,8) &
&             +r3(4,10)*three+r3(6,28)+r3(6,32)*four+r3(6,36)+r3(6,40)*six-r2(5,22)*two &
&             -r2(5,26)*two-r2(5,30)*six-r2(5,34)*six+r1(1,2)+r1(1,6)+r1(1,10)*four &
&             +r1(1,14)+r1(1,18)*three
      rxyz(15)=+r7(26,3)-r6(19,5)*two-r6(19,7)*two+r5(13,5)+r5(13,7)*four+r5(13,9) &
&             +r5(13,11)*six+r5(15,21)-r4(8,6)*two-r4(8,8)*two-r4(8,10)*six-r4(8,12)*six &
&             -r4(10,30)*two-r4(10,34)*two+r3(4,1)+r3(4,3)+r3(4,5)*four+r3(4,7) &
&             +r3(4,9)*three+r3(6,27)+r3(6,31)*four+r3(6,35)+r3(6,39)*six-r2(5,21)*two &
&             -r2(5,25)*two-r2(5,29)*six-r2(5,33)*six+r1(1,1)+r1(1,5)+r1(1,9)*four &
&             +r1(1,13)+r1(1,17)*three
      rxyz(16)=+r7(20,4)-r6(14,6)*two-r6(14,8)*two+r5(9,6)+r5(9,8)*four+r5(9,10) &
&             +r5(9,12)*six+r5(20,22)-r4(5,7)*two-r4(5,9)*two-r4(5,11)*six-r4(5,13)*six &
&             -r4(14,31)*two-r4(14,35)*two+r3(2,2)+r3(2,4)+r3(2,6)*four+r3(2,8) &
&             +r3(2,10)*three+r3(9,28)+r3(9,32)*four+r3(9,36)+r3(9,40)*six-r2(6,22)*two &
&             -r2(6,26)*two-r2(6,30)*six-r2(6,34)*six+r1(2,2)+r1(2,6)+r1(2,10)*four &
&             +r1(2,14)+r1(2,18)*three
      rxyz(17)=+r7(20,3)-r6(14,5)*two-r6(14,7)*two+r5(9,5)+r5(9,7)*four+r5(9,9) &
&             +r5(9,11)*six+r5(20,21)-r4(5,6)*two-r4(5,8)*two-r4(5,10)*six-r4(5,12)*six &
&             -r4(14,30)*two-r4(14,34)*two+r3(2,1)+r3(2,3)+r3(2,5)*four+r3(2,7) &
&             +r3(2,9)*three+r3(9,27)+r3(9,31)*four+r3(9,35)+r3(9,39)*six-r2(6,21)*two &
&             -r2(6,25)*two-r2(6,29)*six-r2(6,33)*six+r1(2,1)+r1(2,5)+r1(2,9)*four &
&             +r1(2,13)+r1(2,17)*three
      rxyz(18)=+r7(35,3)-r6(27,5)*two-r6(27,7)*two+r5(20,5)+r5(20,7)*four+r5(20,9) &
&             +r5(20,11)*six+r5(20,21)-r4(14,6)*two-r4(14,8)*two-r4(14,10)*six &
&             -r4(14,12)*six-r4(14,30)*two-r4(14,34)*two+r3(9,1)+r3(9,3)+r3(9,5)*four &
&             +r3(9,7)+r3(9,9)*three+r3(9,27)+r3(9,31)*four+r3(9,35)+r3(9,39)*six &
&             -r2(6,21)*two-r2(6,25)*two-r2(6,29)*six-r2(6,33)*six+r1(2,1)+r1(2,5) &
&             +r1(2,9)*four+r1(2,13)+r1(2,17)*three
      rxyz(19)=+r7(35,4)-r6(27,6)*two-r6(27,8)*two+r5(20,6)+r5(20,8)*four+r5(20,10) &
&             +r5(20,12)*six+r5(20,22)-r4(14,7)*two-r4(14,9)*two-r4(14,11)*six &
&             -r4(14,13)*six-r4(14,31)*two-r4(14,35)*two+r3(9,2)+r3(9,4)+r3(9,6)*four &
&             +r3(9,8)+r3(9,10)*three+r3(9,28)+r3(9,32)*four+r3(9,36)+r3(9,40)*six &
&             -r2(6,22)*two-r2(6,26)*two-r2(6,30)*six-r2(6,34)*six+r1(2,2)+r1(2,6) &
&             +r1(2,10)*four+r1(2,14)+r1(2,18)*three
      rxyz(20)=+r6(21,11)*four-r5(15,15)*eight-r5(15,19)*eight+r4(10,16)*four &
&             +r4(10,20)*p16+r4(10,24)*four+r4(10,28)*p24-r3(6,13)*eight-r3(6,17)*eight &
&             -r3(6,21)*p24-r3(6,25)*p24+r2(5,3)*four+r2(5,7)*four+r2(5,11)*p16 &
&             +r2(5,15)*four+r2(5,19)*p12
      eri(1,1,3,3)=r400+(+r7(15,3)*two+r7(15,4)*two-r6(10,5)*four-r6(10,6)*four &
&                  -r6(10,7)*four-r6(10,8)*four+r5(6,5)*two+r5(6,6)*two+r5(6,7)*eight &
&                  +r5(6,8)*eight+r5(6,9)*two+r5(6,10)*two+r5(6,11)*p12+r5(6,12)*p12 &
&                  +r5(15,21)*six+r5(15,22)*six-r4(3,6)*four-r4(3,7)*four-r4(3,8)*four &
&                  -r4(3,9)*four-r4(3,10)*p12-r4(3,11)*p12-r4(3,12)*p12-r4(3,13)*p12 &
&                  -r4(10,30)*p12-r4(10,31)*p12-r4(10,34)*p12-r4(10,35)*p12+r3(1,1)*two &
&                  +r3(1,2)*two+r3(1,3)*two+r3(1,4)*two+r3(1,5)*eight+r3(1,6)*eight &
&                  +r3(1,7)*two+r3(1,8)*two+r3(1,9)*six+r3(1,10)*six+r3(6,27)*six &
&                  +r3(6,28)*six+r3(6,31)*p24+r3(6,32)*p24+r3(6,35)*six+r3(6,36)*six &
&                  +r3(6,39)*p36+r3(6,40)*p36-r2(5,21)*p12-r2(5,22)*p12-r2(5,25)*p12 &
&                  -r2(5,26)*p12-r2(5,29)*p36-r2(5,30)*p36-r2(5,33)*p36-r2(5,34)*p36 &
&                  +r1(1,1)*six+r1(1,2)*six+r1(1,5)*six+r1(1,6)*six+r1(1,9)*p24+r1(1,10)*p24 &
&                  +r1(1,13)*six+r1(1,14)*six+r1(1,17)*p18+r1(1,18)*p18)*qx+(+r6(15,10) &
&                  +r6(15,11)*four+r6(15,12)-r5(10,14)*two-r5(10,15)*eight-r5(10,16)*two &
&                  -r5(10,18)*two-r5(10,19)*eight-r5(10,20)*two+r4(6,15)+r4(6,16)*four &
&                  +r4(6,17)+r4(6,19)*four+r4(6,20)*p16+r4(6,21)*four+r4(6,23)+r4(6,24)*four &
&                  +r4(6,25)+r4(6,27)*six+r4(6,28)*p24+r4(6,29)*six+r4(15,39)+r4(15,40)*four &
&                  +r4(15,41)-r3(3,12)*two-r3(3,13)*eight-r3(3,14)*two-r3(3,16)*two &
&                  -r3(3,17)*eight-r3(3,18)*two-r3(3,20)*six-r3(3,21)*p24-r3(3,22)*six &
&                  -r3(3,24)*six-r3(3,25)*p24-r3(3,26)*six-r3(10,44)*two-r3(10,45)*eight &
&                  -r3(10,46)*two-r3(10,49)*two-r3(10,50)*eight-r3(10,51)*two+r2(1,2) &
&                  +r2(1,3)*four+r2(1,4)+r2(1,6)+r2(1,7)*four+r2(1,8)+r2(1,10)*four &
&                  +r2(1,11)*p16+r2(1,12)*four+r2(1,14)+r2(1,15)*four+r2(1,16)+r2(1,18)*three &
&                  +r2(1,19)*p12+r2(1,20)*three+r2(3,38)+r2(3,39)*four+r2(3,40)+r2(3,43)*four &
&                  +r2(3,44)*p16+r2(3,45)*four+r2(3,48)+r2(3,49)*four+r2(3,50)+r2(3,53)*six &
&                  +r2(3,54)*p24+r2(3,55)*six-r1(3,22)*two-r1(3,23)*eight-r1(3,24)*two &
&                  -r1(3,27)*two-r1(3,28)*eight-r1(3,29)*two-r1(3,32)*six-r1(3,33)*p24 &
&                  -r1(3,34)*six-r1(3,37)*six-r1(3,38)*p24-r1(3,39)*six+r0(2)+r0(3)*four &
&                  +r0(4)+r0(7)+r0(8)*four+r0(9)+r0(12)*four+r0(13)*p16+r0(14)*four+r0(17) &
&                  +r0(18)*four+r0(19)+r0(22)*three+r0(23)*p12+r0(24)*three)*xx+( &
&                  +r5(15,23)*two+r5(15,24)*two-r4(10,32)*four-r4(10,33)*four-r4(10,36)*four &
&                  -r4(10,37)*four+r3(6,29)*two+r3(6,30)*two+r3(6,33)*eight+r3(6,34)*eight &
&                  +r3(6,37)*two+r3(6,38)*two+r3(6,41)*p12+r3(6,42)*p12-r2(5,23)*four &
&                  -r2(5,24)*four-r2(5,27)*four-r2(5,28)*four-r2(5,31)*p12-r2(5,32)*p12 &
&                  -r2(5,35)*p12-r2(5,36)*p12+r1(1,3)*two+r1(1,4)*two+r1(1,7)*two+r1(1,8)*two &
&                  +r1(1,11)*eight+r1(1,12)*eight+r1(1,15)*two+r1(1,16)*two+r1(1,19)*six &
&                  +r1(1,20)*six)*xxx+rxyz(1)*xxxx
      eri(2,1,3,3)=r220+(+r7(26,4)*two-r6(19,6)*four-r6(19,8)*four+r5(13,6)*two &
&                  +r5(13,8)*eight+r5(13,10)*two+r5(13,12)*p12+r5(15,22)*two-r4(8,7)*four &
&                  -r4(8,9)*four-r4(8,11)*p12-r4(8,13)*p12-r4(10,31)*four-r4(10,35)*four &
&                  +r3(4,2)*two+r3(4,4)*two+r3(4,6)*eight+r3(4,8)*two+r3(4,10)*six &
&                  +r3(6,28)*two+r3(6,32)*eight+r3(6,36)*two+r3(6,40)*p12-r2(5,22)*four &
&                  -r2(5,26)*four-r2(5,30)*p12-r2(5,34)*p12+r1(1,2)*two+r1(1,6)*two &
&                  +r1(1,10)*eight+r1(1,14)*two+r1(1,18)*six)*qx+rxyz(5)*xx
      eri(3,1,3,3)=r202+(+r7(28,4)*two-r6(21,6)*four-r6(21,8)*four+r5(15,6)*two &
&                  +r5(15,8)*eight+r5(15,10)*two+r5(15,12)*p12+r5(15,22)*two-r4(10,7)*four &
&                  -r4(10,9)*four-r4(10,11)*p12-r4(10,13)*p12-r4(10,31)*four-r4(10,35)*four &
&                  +r3(6,2)*two+r3(6,4)*two+r3(6,6)*eight+r3(6,8)*two+r3(6,10)*six &
&                  +r3(6,28)*two+r3(6,32)*eight+r3(6,36)*two+r3(6,40)*p12-r2(5,22)*four &
&                  -r2(5,26)*four-r2(5,30)*p12-r2(5,34)*p12+r1(1,2)*two+r1(1,6)*two &
&                  +r1(1,10)*eight+r1(1,14)*two+r1(1,18)*six)*qx+(+r7(21,3)*two-r6(15,5)*four &
&                  -r6(15,7)*four+r5(10,5)*two+r5(10,7)*eight+r5(10,9)*two+r5(10,11)*p12 &
&                  +r5(21,21)*two-r4(6,6)*four-r4(6,8)*four-r4(6,10)*p12-r4(6,12)*p12 &
&                  -r4(15,30)*four-r4(15,34)*four+r3(3,1)*two+r3(3,3)*two+r3(3,5)*eight &
&                  +r3(3,7)*two+r3(3,9)*six+r3(10,27)*two+r3(10,31)*eight+r3(10,35)*two &
&                  +r3(10,39)*p12-r2(3,21)*four-r2(3,25)*four-r2(3,29)*p12-r2(3,33)*p12 &
&                  +r1(3,1)*two+r1(3,5)*two+r1(3,9)*eight+r1(3,13)*two+r1(3,17)*six)*qz+( &
&                  +r6(28,12)-r5(21,16)*two-r5(21,20)*two+r4(15,17)+r4(15,21)*four+r4(15,25) &
&                  +r4(15,29)*six+r4(15,41)-r3(10,14)*two-r3(10,18)*two-r3(10,22)*six &
&                  -r3(10,26)*six-r3(10,46)*two-r3(10,51)*two+r2(3,4)+r2(3,8)+r2(3,12)*four &
&                  +r2(3,16)+r2(3,20)*three+r2(3,40)+r2(3,45)*four+r2(3,50)+r2(3,55)*six &
&                  -r1(3,24)*two-r1(3,29)*two-r1(3,34)*six-r1(3,39)*six+r0(4)+r0(9) &
&                  +r0(14)*four+r0(19)+r0(24)*three)*xx+rxyz(20)*xz+(+r6(15,10)-r5(10,14)*two &
&                  -r5(10,18)*two+r4(6,15)+r4(6,19)*four+r4(6,23)+r4(6,27)*six+r4(15,39) &
&                  -r3(3,12)*two-r3(3,16)*two-r3(3,20)*six-r3(3,24)*six-r3(10,44)*two &
&                  -r3(10,49)*two+r2(1,2)+r2(1,6)+r2(1,10)*four+r2(1,14)+r2(1,18)*three &
&                  +r2(3,38)+r2(3,43)*four+r2(3,48)+r2(3,53)*six-r1(3,22)*two-r1(3,27)*two &
&                  -r1(3,32)*six-r1(3,37)*six+r0(2)+r0(7)+r0(12)*four+r0(17)+r0(22)*three)*zz &
&                  +(+r5(21,24)*two-r4(15,33)*four-r4(15,37)*four+r3(10,30)*two &
&                  +r3(10,34)*eight+r3(10,38)*two+r3(10,42)*p12-r2(3,24)*four-r2(3,28)*four &
&                  -r2(3,32)*p12-r2(3,36)*p12+r1(3,4)*two+r1(3,8)*two+r1(3,12)*eight &
&                  +r1(3,16)*two+r1(3,20)*six)*xxz+(+r5(15,23)*two-r4(10,32)*four &
&                  -r4(10,36)*four+r3(6,29)*two+r3(6,33)*eight+r3(6,37)*two+r3(6,41)*p12 &
&                  -r2(5,23)*four-r2(5,27)*four-r2(5,31)*p12-r2(5,35)*p12+r1(1,3)*two &
&                  +r1(1,7)*two+r1(1,11)*eight+r1(1,15)*two+r1(1,19)*six)*xzz+rxyz(1)*xxzz
      eri(4,1,3,3)=r310+(+r7(20,3)+r7(20,4)*two-r6(14,5)*two-r6(14,6)*four-r6(14,7)*two &
&                  -r6(14,8)*four+r5(9,5)+r5(9,6)*two+r5(9,7)*four+r5(9,8)*eight+r5(9,9) &
&                  +r5(9,10)*two+r5(9,11)*six+r5(9,12)*p12+r5(20,21)+r5(20,22)*two &
&                  -r4(5,6)*two-r4(5,7)*four-r4(5,8)*two-r4(5,9)*four-r4(5,10)*six &
&                  -r4(5,11)*p12-r4(5,12)*six-r4(5,13)*p12-r4(14,30)*two-r4(14,31)*four &
&                  -r4(14,34)*two-r4(14,35)*four+r3(2,1)+r3(2,2)*two+r3(2,3)+r3(2,4)*two &
&                  +r3(2,5)*four+r3(2,6)*eight+r3(2,7)+r3(2,8)*two+r3(2,9)*three+r3(2,10)*six &
&                  +r3(9,27)+r3(9,28)*two+r3(9,31)*four+r3(9,32)*eight+r3(9,35)+r3(9,36)*two &
&                  +r3(9,39)*six+r3(9,40)*p12-r2(6,21)*two-r2(6,22)*four-r2(6,25)*two &
&                  -r2(6,26)*four-r2(6,29)*six-r2(6,30)*p12-r2(6,33)*six-r2(6,34)*p12+r1(2,1) &
&                  +r1(2,2)*two+r1(2,5)+r1(2,6)*two+r1(2,9)*four+r1(2,10)*eight+r1(2,13) &
&                  +r1(2,14)*two+r1(2,17)*three+r1(2,18)*six)*qx+(+r6(20,11)*two+r6(20,12) &
&                  -r5(14,15)*four-r5(14,16)*two-r5(14,19)*four-r5(14,20)*two+r4(9,16)*two &
&                  +r4(9,17)+r4(9,20)*eight+r4(9,21)*four+r4(9,24)*two+r4(9,25)+r4(9,28)*p12 &
&                  +r4(9,29)*six-r3(5,13)*four-r3(5,14)*two-r3(5,17)*four-r3(5,18)*two &
&                  -r3(5,21)*p12-r3(5,22)*six-r3(5,25)*p12-r3(5,26)*six+r2(4,3)*two+r2(4,4) &
&                  +r2(4,7)*two+r2(4,8)+r2(4,11)*eight+r2(4,12)*four+r2(4,15)*two+r2(4,16) &
&                  +r2(4,19)*six+r2(4,20)*three)*xx+rxyz(2)*xxx
      eri(5,1,3,3)=r301+(+r7(21,3)+r7(21,4)*two-r6(15,5)*two-r6(15,6)*four-r6(15,7)*two &
&                  -r6(15,8)*four+r5(10,5)+r5(10,6)*two+r5(10,7)*four+r5(10,8)*eight+r5(10,9) &
&                  +r5(10,10)*two+r5(10,11)*six+r5(10,12)*p12+r5(21,21)+r5(21,22)*two &
&                  -r4(6,6)*two-r4(6,7)*four-r4(6,8)*two-r4(6,9)*four-r4(6,10)*six &
&                  -r4(6,11)*p12-r4(6,12)*six-r4(6,13)*p12-r4(15,30)*two-r4(15,31)*four &
&                  -r4(15,34)*two-r4(15,35)*four+r3(3,1)+r3(3,2)*two+r3(3,3)+r3(3,4)*two &
&                  +r3(3,5)*four+r3(3,6)*eight+r3(3,7)+r3(3,8)*two+r3(3,9)*three+r3(3,10)*six &
&                  +r3(10,27)+r3(10,28)*two+r3(10,31)*four+r3(10,32)*eight+r3(10,35) &
&                  +r3(10,36)*two+r3(10,39)*six+r3(10,40)*p12-r2(3,21)*two-r2(3,22)*four &
&                  -r2(3,25)*two-r2(3,26)*four-r2(3,29)*six-r2(3,30)*p12-r2(3,33)*six &
&                  -r2(3,34)*p12+r1(3,1)+r1(3,2)*two+r1(3,5)+r1(3,6)*two+r1(3,9)*four &
&                  +r1(3,10)*eight+r1(3,13)+r1(3,14)*two+r1(3,17)*three+r1(3,18)*six)*qx+( &
&                  +r7(15,3)-r6(10,5)*two-r6(10,7)*two+r5(6,5)+r5(6,7)*four+r5(6,9) &
&                  +r5(6,11)*six+r5(15,21)*three-r4(3,6)*two-r4(3,8)*two-r4(3,10)*six &
&                  -r4(3,12)*six-r4(10,30)*six-r4(10,34)*six+r3(1,1)+r3(1,3)+r3(1,5)*four &
&                  +r3(1,7)+r3(1,9)*three+r3(6,27)*three+r3(6,31)*p12+r3(6,35)*three &
&                  +r3(6,39)*p18-r2(5,21)*six-r2(5,25)*six-r2(5,29)*p18-r2(5,33)*p18 &
&                  +r1(1,1)*three+r1(1,5)*three+r1(1,9)*p12+r1(1,13)*three+r1(1,17)*nine)*qz &
&                  +(+r6(21,11)*two+r6(21,12)-r5(15,15)*four-r5(15,16)*two-r5(15,19)*four &
&                  -r5(15,20)*two+r4(10,16)*two+r4(10,17)+r4(10,20)*eight+r4(10,21)*four &
&                  +r4(10,24)*two+r4(10,25)+r4(10,28)*p12+r4(10,29)*six-r3(6,13)*four &
&                  -r3(6,14)*two-r3(6,17)*four-r3(6,18)*two-r3(6,21)*p12-r3(6,22)*six &
&                  -r3(6,25)*p12-r3(6,26)*six+r2(5,3)*two+r2(5,4)+r2(5,7)*two+r2(5,8) &
&                  +r2(5,11)*eight+r2(5,12)*four+r2(5,15)*two+r2(5,16)+r2(5,19)*six &
&                  +r2(5,20)*three)*xx+(+r6(15,10)+r6(15,11)*two-r5(10,14)*two-r5(10,15)*four &
&                  -r5(10,18)*two-r5(10,19)*four+r4(6,15)+r4(6,16)*two+r4(6,19)*four &
&                  +r4(6,20)*eight+r4(6,23)+r4(6,24)*two+r4(6,27)*six+r4(6,28)*p12+r4(15,39) &
&                  +r4(15,40)*two-r3(3,12)*two-r3(3,13)*four-r3(3,16)*two-r3(3,17)*four &
&                  -r3(3,20)*six-r3(3,21)*p12-r3(3,24)*six-r3(3,25)*p12-r3(10,44)*two &
&                  -r3(10,45)*four-r3(10,49)*two-r3(10,50)*four+r2(1,2)+r2(1,3)*two+r2(1,6) &
&                  +r2(1,7)*two+r2(1,10)*four+r2(1,11)*eight+r2(1,14)+r2(1,15)*two &
&                  +r2(1,18)*three+r2(1,19)*six+r2(3,38)+r2(3,39)*two+r2(3,43)*four &
&                  +r2(3,44)*eight+r2(3,48)+r2(3,49)*two+r2(3,53)*six+r2(3,54)*p12 &
&                  -r1(3,22)*two-r1(3,23)*four-r1(3,27)*two-r1(3,28)*four-r1(3,32)*six &
&                  -r1(3,33)*p12-r1(3,37)*six-r1(3,38)*p12+r0(2)+r0(3)*two+r0(7)+r0(8)*two &
&                  +r0(12)*four+r0(13)*eight+r0(17)+r0(18)*two+r0(22)*three+r0(23)*six)*xz+( &
&                  +r5(21,24)-r4(15,33)*two-r4(15,37)*two+r3(10,30)+r3(10,34)*four+r3(10,38) &
&                  +r3(10,42)*six-r2(3,24)*two-r2(3,28)*two-r2(3,32)*six-r2(3,36)*six+r1(3,4) &
&                  +r1(3,8)+r1(3,12)*four+r1(3,16)+r1(3,20)*three)*xxx+(+r5(15,23)*two &
&                  +r5(15,24)-r4(10,32)*four-r4(10,33)*two-r4(10,36)*four-r4(10,37)*two &
&                  +r3(6,29)*two+r3(6,30)+r3(6,33)*eight+r3(6,34)*four+r3(6,37)*two+r3(6,38) &
&                  +r3(6,41)*p12+r3(6,42)*six-r2(5,23)*four-r2(5,24)*two-r2(5,27)*four &
&                  -r2(5,28)*two-r2(5,31)*p12-r2(5,32)*six-r2(5,35)*p12-r2(5,36)*six &
&                  +r1(1,3)*two+r1(1,4)+r1(1,7)*two+r1(1,8)+r1(1,11)*eight+r1(1,12)*four &
&                  +r1(1,15)*two+r1(1,16)+r1(1,19)*six+r1(1,20)*three)*xxz+rxyz(1)*xxxz
      eri(6,1,3,3)=r211+(+r7(27,4)*two-r6(20,6)*four-r6(20,8)*four+r5(14,6)*two &
&                  +r5(14,8)*eight+r5(14,10)*two+r5(14,12)*p12-r4(9,7)*four-r4(9,9)*four &
&                  -r4(9,11)*p12-r4(9,13)*p12+r3(5,2)*two+r3(5,4)*two+r3(5,6)*eight &
&                  +r3(5,8)*two+r3(5,10)*six)*qx+rxyz(17)*qz+(+r6(27,12)-r5(20,16)*two &
&                  -r5(20,20)*two+r4(14,17)+r4(14,21)*four+r4(14,25)+r4(14,29)*six &
&                  -r3(9,14)*two-r3(9,18)*two-r3(9,22)*six-r3(9,26)*six+r2(6,4)+r2(6,8) &
&                  +r2(6,12)*four+r2(6,16)+r2(6,20)*three)*xx+(+r6(20,11)*two-r5(14,15)*four &
&                  -r5(14,19)*four+r4(9,16)*two+r4(9,20)*eight+r4(9,24)*two+r4(9,28)*p12 &
&                  -r3(5,13)*four-r3(5,17)*four-r3(5,21)*p12-r3(5,25)*p12+r2(4,3)*two &
&                  +r2(4,7)*two+r2(4,11)*eight+r2(4,15)*two+r2(4,19)*six)*xz+rxyz(2)*xxz
      eri(1,2,3,3)=r220+(+r7(26,3)*two-r6(19,5)*four-r6(19,7)*four+r5(13,5)*two &
&                  +r5(13,7)*eight+r5(13,9)*two+r5(13,11)*p12+r5(15,21)*two-r4(8,6)*four &
&                  -r4(8,8)*four-r4(8,10)*p12-r4(8,12)*p12-r4(10,30)*four-r4(10,34)*four &
&                  +r3(4,1)*two+r3(4,3)*two+r3(4,5)*eight+r3(4,7)*two+r3(4,9)*six &
&                  +r3(6,27)*two+r3(6,31)*eight+r3(6,35)*two+r3(6,39)*p12-r2(5,21)*four &
&                  -r2(5,25)*four-r2(5,29)*p12-r2(5,33)*p12+r1(1,1)*two+r1(1,5)*two &
&                  +r1(1,9)*eight+r1(1,13)*two+r1(1,17)*six)*qx+rxyz(4)*xx
      eri(2,2,3,3)=r040
      eri(3,2,3,3)=r022+(+r7(34,3)*two-r6(26,5)*four-r6(26,7)*four+r5(19,5)*two &
&                  +r5(19,7)*eight+r5(19,9)*two+r5(19,11)*p12+r5(21,21)*two-r4(13,6)*four &
&                  -r4(13,8)*four-r4(13,10)*p12-r4(13,12)*p12-r4(15,30)*four-r4(15,34)*four &
&                  +r3(8,1)*two+r3(8,3)*two+r3(8,5)*eight+r3(8,7)*two+r3(8,9)*six &
&                  +r3(10,27)*two+r3(10,31)*eight+r3(10,35)*two+r3(10,39)*p12-r2(3,21)*four &
&                  -r2(3,25)*four-r2(3,29)*p12-r2(3,33)*p12+r1(3,1)*two+r1(3,5)*two &
&                  +r1(3,9)*eight+r1(3,13)*two+r1(3,17)*six)*qz+rxyz(4)*zz
      eri(4,2,3,3)=r130+rxyz(7)*qx
      eri(5,2,3,3)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,3,3)=r031+rxyz(7)*qz
      eri(1,3,3,3)=r202+(+r7(28,3)*two-r6(21,5)*four-r6(21,7)*four+r5(15,5)*two &
&                  +r5(15,7)*eight+r5(15,9)*two+r5(15,11)*p12+r5(15,21)*two-r4(10,6)*four &
&                  -r4(10,8)*four-r4(10,10)*p12-r4(10,12)*p12-r4(10,30)*four-r4(10,34)*four &
&                  +r3(6,1)*two+r3(6,3)*two+r3(6,5)*eight+r3(6,7)*two+r3(6,9)*six &
&                  +r3(6,27)*two+r3(6,31)*eight+r3(6,35)*two+r3(6,39)*p12-r2(5,21)*four &
&                  -r2(5,25)*four-r2(5,29)*p12-r2(5,33)*p12+r1(1,1)*two+r1(1,5)*two &
&                  +r1(1,9)*eight+r1(1,13)*two+r1(1,17)*six)*qx+(+r7(21,4)*two-r6(15,6)*four &
&                  -r6(15,8)*four+r5(10,6)*two+r5(10,8)*eight+r5(10,10)*two+r5(10,12)*p12 &
&                  +r5(21,22)*two-r4(6,7)*four-r4(6,9)*four-r4(6,11)*p12-r4(6,13)*p12 &
&                  -r4(15,31)*four-r4(15,35)*four+r3(3,2)*two+r3(3,4)*two+r3(3,6)*eight &
&                  +r3(3,8)*two+r3(3,10)*six+r3(10,28)*two+r3(10,32)*eight+r3(10,36)*two &
&                  +r3(10,40)*p12-r2(3,22)*four-r2(3,26)*four-r2(3,30)*p12-r2(3,34)*p12 &
&                  +r1(3,2)*two+r1(3,6)*two+r1(3,10)*eight+r1(3,14)*two+r1(3,18)*six)*qz+( &
&                  +r6(28,10)-r5(21,14)*two-r5(21,18)*two+r4(15,15)+r4(15,19)*four+r4(15,23) &
&                  +r4(15,27)*six+r4(15,39)-r3(10,12)*two-r3(10,16)*two-r3(10,20)*six &
&                  -r3(10,24)*six-r3(10,44)*two-r3(10,49)*two+r2(3,2)+r2(3,6)+r2(3,10)*four &
&                  +r2(3,14)+r2(3,18)*three+r2(3,38)+r2(3,43)*four+r2(3,48)+r2(3,53)*six &
&                  -r1(3,22)*two-r1(3,27)*two-r1(3,32)*six-r1(3,37)*six+r0(2)+r0(7) &
&                  +r0(12)*four+r0(17)+r0(22)*three)*xx+rxyz(20)*xz+(+r6(15,12)-r5(10,16)*two &
&                  -r5(10,20)*two+r4(6,17)+r4(6,21)*four+r4(6,25)+r4(6,29)*six+r4(15,41) &
&                  -r3(3,14)*two-r3(3,18)*two-r3(3,22)*six-r3(3,26)*six-r3(10,46)*two &
&                  -r3(10,51)*two+r2(1,4)+r2(1,8)+r2(1,12)*four+r2(1,16)+r2(1,20)*three &
&                  +r2(3,40)+r2(3,45)*four+r2(3,50)+r2(3,55)*six-r1(3,24)*two-r1(3,29)*two &
&                  -r1(3,34)*six-r1(3,39)*six+r0(4)+r0(9)+r0(14)*four+r0(19)+r0(24)*three)*zz &
&                  +(+r5(21,23)*two-r4(15,32)*four-r4(15,36)*four+r3(10,29)*two &
&                  +r3(10,33)*eight+r3(10,37)*two+r3(10,41)*p12-r2(3,23)*four-r2(3,27)*four &
&                  -r2(3,31)*p12-r2(3,35)*p12+r1(3,3)*two+r1(3,7)*two+r1(3,11)*eight &
&                  +r1(3,15)*two+r1(3,19)*six)*xxz+(+r5(15,24)*two-r4(10,33)*four &
&                  -r4(10,37)*four+r3(6,30)*two+r3(6,34)*eight+r3(6,38)*two+r3(6,42)*p12 &
&                  -r2(5,24)*four-r2(5,28)*four-r2(5,32)*p12-r2(5,36)*p12+r1(1,4)*two &
&                  +r1(1,8)*two+r1(1,12)*eight+r1(1,16)*two+r1(1,20)*six)*xzz+rxyz(1)*xxzz
      eri(2,3,3,3)=r022+(+r7(34,4)*two-r6(26,6)*four-r6(26,8)*four+r5(19,6)*two &
&                  +r5(19,8)*eight+r5(19,10)*two+r5(19,12)*p12+r5(21,22)*two-r4(13,7)*four &
&                  -r4(13,9)*four-r4(13,11)*p12-r4(13,13)*p12-r4(15,31)*four-r4(15,35)*four &
&                  +r3(8,2)*two+r3(8,4)*two+r3(8,6)*eight+r3(8,8)*two+r3(8,10)*six &
&                  +r3(10,28)*two+r3(10,32)*eight+r3(10,36)*two+r3(10,40)*p12-r2(3,22)*four &
&                  -r2(3,26)*four-r2(3,30)*p12-r2(3,34)*p12+r1(3,2)*two+r1(3,6)*two &
&                  +r1(3,10)*eight+r1(3,14)*two+r1(3,18)*six)*qz+rxyz(5)*zz
      eri(3,3,3,3)=r004+(+r7(36,3)*two+r7(36,4)*two-r6(28,5)*four-r6(28,6)*four &
&                  -r6(28,7)*four-r6(28,8)*four+r5(21,5)*two+r5(21,6)*two+r5(21,7)*eight &
&                  +r5(21,8)*eight+r5(21,9)*two+r5(21,10)*two+r5(21,11)*p12+r5(21,12)*p12 &
&                  +r5(21,21)*six+r5(21,22)*six-r4(15,6)*four-r4(15,7)*four-r4(15,8)*four &
&                  -r4(15,9)*four-r4(15,10)*p12-r4(15,11)*p12-r4(15,12)*p12-r4(15,13)*p12 &
&                  -r4(15,30)*p12-r4(15,31)*p12-r4(15,34)*p12-r4(15,35)*p12+r3(10,1)*two &
&                  +r3(10,2)*two+r3(10,3)*two+r3(10,4)*two+r3(10,5)*eight+r3(10,6)*eight &
&                  +r3(10,7)*two+r3(10,8)*two+r3(10,9)*six+r3(10,10)*six+r3(10,27)*six &
&                  +r3(10,28)*six+r3(10,31)*p24+r3(10,32)*p24+r3(10,35)*six+r3(10,36)*six &
&                  +r3(10,39)*p36+r3(10,40)*p36-r2(3,21)*p12-r2(3,22)*p12-r2(3,25)*p12 &
&                  -r2(3,26)*p12-r2(3,29)*p36-r2(3,30)*p36-r2(3,33)*p36-r2(3,34)*p36 &
&                  +r1(3,1)*six+r1(3,2)*six+r1(3,5)*six+r1(3,6)*six+r1(3,9)*p24+r1(3,10)*p24 &
&                  +r1(3,13)*six+r1(3,14)*six+r1(3,17)*p18+r1(3,18)*p18)*qz+(+r6(28,10) &
&                  +r6(28,11)*four+r6(28,12)-r5(21,14)*two-r5(21,15)*eight-r5(21,16)*two &
&                  -r5(21,18)*two-r5(21,19)*eight-r5(21,20)*two+r4(15,15)+r4(15,16)*four &
&                  +r4(15,17)+r4(15,19)*four+r4(15,20)*p16+r4(15,21)*four+r4(15,23) &
&                  +r4(15,24)*four+r4(15,25)+r4(15,27)*six+r4(15,28)*p24+r4(15,29)*six &
&                  +r4(15,39)+r4(15,40)*four+r4(15,41)-r3(10,12)*two-r3(10,13)*eight &
&                  -r3(10,14)*two-r3(10,16)*two-r3(10,17)*eight-r3(10,18)*two-r3(10,20)*six &
&                  -r3(10,21)*p24-r3(10,22)*six-r3(10,24)*six-r3(10,25)*p24-r3(10,26)*six &
&                  -r3(10,44)*two-r3(10,45)*eight-r3(10,46)*two-r3(10,49)*two-r3(10,50)*eight &
&                  -r3(10,51)*two+r2(3,2)+r2(3,3)*four+r2(3,4)+r2(3,6)+r2(3,7)*four+r2(3,8) &
&                  +r2(3,10)*four+r2(3,11)*p16+r2(3,12)*four+r2(3,14)+r2(3,15)*four+r2(3,16) &
&                  +r2(3,18)*three+r2(3,19)*p12+r2(3,20)*three+r2(3,38)+r2(3,39)*four &
&                  +r2(3,40)+r2(3,43)*four+r2(3,44)*p16+r2(3,45)*four+r2(3,48)+r2(3,49)*four &
&                  +r2(3,50)+r2(3,53)*six+r2(3,54)*p24+r2(3,55)*six-r1(3,22)*two &
&                  -r1(3,23)*eight-r1(3,24)*two-r1(3,27)*two-r1(3,28)*eight-r1(3,29)*two &
&                  -r1(3,32)*six-r1(3,33)*p24-r1(3,34)*six-r1(3,37)*six-r1(3,38)*p24 &
&                  -r1(3,39)*six+r0(2)+r0(3)*four+r0(4)+r0(7)+r0(8)*four+r0(9)+r0(12)*four &
&                  +r0(13)*p16+r0(14)*four+r0(17)+r0(18)*four+r0(19)+r0(22)*three+r0(23)*p12 &
&                  +r0(24)*three)*zz+(+r5(21,23)*two+r5(21,24)*two-r4(15,32)*four &
&                  -r4(15,33)*four-r4(15,36)*four-r4(15,37)*four+r3(10,29)*two+r3(10,30)*two &
&                  +r3(10,33)*eight+r3(10,34)*eight+r3(10,37)*two+r3(10,38)*two+r3(10,41)*p12 &
&                  +r3(10,42)*p12-r2(3,23)*four-r2(3,24)*four-r2(3,27)*four-r2(3,28)*four &
&                  -r2(3,31)*p12-r2(3,32)*p12-r2(3,35)*p12-r2(3,36)*p12+r1(3,3)*two &
&                  +r1(3,4)*two+r1(3,7)*two+r1(3,8)*two+r1(3,11)*eight+r1(3,12)*eight &
&                  +r1(3,15)*two+r1(3,16)*two+r1(3,19)*six+r1(3,20)*six)*zzz+rxyz(1)*zzzz
      eri(4,3,3,3)=r112+rxyz(18)*qx+(+r7(27,4)*two-r6(20,6)*four-r6(20,8)*four &
&                  +r5(14,6)*two+r5(14,8)*eight+r5(14,10)*two+r5(14,12)*p12-r4(9,7)*four &
&                  -r4(9,9)*four-r4(9,11)*p12-r4(9,13)*p12+r3(5,2)*two+r3(5,4)*two &
&                  +r3(5,6)*eight+r3(5,8)*two+r3(5,10)*six)*qz+(+r6(27,11)*two-r5(20,15)*four &
&                  -r5(20,19)*four+r4(14,16)*two+r4(14,20)*eight+r4(14,24)*two+r4(14,28)*p12 &
&                  -r3(9,13)*four-r3(9,17)*four-r3(9,21)*p12-r3(9,25)*p12+r2(6,3)*two &
&                  +r2(6,7)*two+r2(6,11)*eight+r2(6,15)*two+r2(6,19)*six)*xz+(+r6(20,12) &
&                  -r5(14,16)*two-r5(14,20)*two+r4(9,17)+r4(9,21)*four+r4(9,25)+r4(9,29)*six &
&                  -r3(5,14)*two-r3(5,18)*two-r3(5,22)*six-r3(5,26)*six+r2(4,4)+r2(4,8) &
&                  +r2(4,12)*four+r2(4,16)+r2(4,20)*three)*zz+rxyz(2)*xzz
      eri(5,3,3,3)=r103+(+r7(36,3)-r6(28,5)*two-r6(28,7)*two+r5(21,5)+r5(21,7)*four &
&                  +r5(21,9)+r5(21,11)*six+r5(21,21)*three-r4(15,6)*two-r4(15,8)*two &
&                  -r4(15,10)*six-r4(15,12)*six-r4(15,30)*six-r4(15,34)*six+r3(10,1)+r3(10,3) &
&                  +r3(10,5)*four+r3(10,7)+r3(10,9)*three+r3(10,27)*three+r3(10,31)*p12 &
&                  +r3(10,35)*three+r3(10,39)*p18-r2(3,21)*six-r2(3,25)*six-r2(3,29)*p18 &
&                  -r2(3,33)*p18+r1(3,1)*three+r1(3,5)*three+r1(3,9)*p12+r1(3,13)*three &
&                  +r1(3,17)*nine)*qx+(+r7(28,3)+r7(28,4)*two-r6(21,5)*two-r6(21,6)*four &
&                  -r6(21,7)*two-r6(21,8)*four+r5(15,5)+r5(15,6)*two+r5(15,7)*four &
&                  +r5(15,8)*eight+r5(15,9)+r5(15,10)*two+r5(15,11)*six+r5(15,12)*p12 &
&                  +r5(15,21)+r5(15,22)*two-r4(10,6)*two-r4(10,7)*four-r4(10,8)*two &
&                  -r4(10,9)*four-r4(10,10)*six-r4(10,11)*p12-r4(10,12)*six-r4(10,13)*p12 &
&                  -r4(10,30)*two-r4(10,31)*four-r4(10,34)*two-r4(10,35)*four+r3(6,1) &
&                  +r3(6,2)*two+r3(6,3)+r3(6,4)*two+r3(6,5)*four+r3(6,6)*eight+r3(6,7) &
&                  +r3(6,8)*two+r3(6,9)*three+r3(6,10)*six+r3(6,27)+r3(6,28)*two &
&                  +r3(6,31)*four+r3(6,32)*eight+r3(6,35)+r3(6,36)*two+r3(6,39)*six &
&                  +r3(6,40)*p12-r2(5,21)*two-r2(5,22)*four-r2(5,25)*two-r2(5,26)*four &
&                  -r2(5,29)*six-r2(5,30)*p12-r2(5,33)*six-r2(5,34)*p12+r1(1,1)+r1(1,2)*two &
&                  +r1(1,5)+r1(1,6)*two+r1(1,9)*four+r1(1,10)*eight+r1(1,13)+r1(1,14)*two &
&                  +r1(1,17)*three+r1(1,18)*six)*qz+(+r6(28,10)+r6(28,11)*two-r5(21,14)*two &
&                  -r5(21,15)*four-r5(21,18)*two-r5(21,19)*four+r4(15,15)+r4(15,16)*two &
&                  +r4(15,19)*four+r4(15,20)*eight+r4(15,23)+r4(15,24)*two+r4(15,27)*six &
&                  +r4(15,28)*p12+r4(15,39)+r4(15,40)*two-r3(10,12)*two-r3(10,13)*four &
&                  -r3(10,16)*two-r3(10,17)*four-r3(10,20)*six-r3(10,21)*p12-r3(10,24)*six &
&                  -r3(10,25)*p12-r3(10,44)*two-r3(10,45)*four-r3(10,49)*two-r3(10,50)*four &
&                  +r2(3,2)+r2(3,3)*two+r2(3,6)+r2(3,7)*two+r2(3,10)*four+r2(3,11)*eight &
&                  +r2(3,14)+r2(3,15)*two+r2(3,18)*three+r2(3,19)*six+r2(3,38)+r2(3,39)*two &
&                  +r2(3,43)*four+r2(3,44)*eight+r2(3,48)+r2(3,49)*two+r2(3,53)*six &
&                  +r2(3,54)*p12-r1(3,22)*two-r1(3,23)*four-r1(3,27)*two-r1(3,28)*four &
&                  -r1(3,32)*six-r1(3,33)*p12-r1(3,37)*six-r1(3,38)*p12+r0(2)+r0(3)*two+r0(7) &
&                  +r0(8)*two+r0(12)*four+r0(13)*eight+r0(17)+r0(18)*two+r0(22)*three &
&                  +r0(23)*six)*xz+(+r6(21,11)*two+r6(21,12)-r5(15,15)*four-r5(15,16)*two &
&                  -r5(15,19)*four-r5(15,20)*two+r4(10,16)*two+r4(10,17)+r4(10,20)*eight &
&                  +r4(10,21)*four+r4(10,24)*two+r4(10,25)+r4(10,28)*p12+r4(10,29)*six &
&                  -r3(6,13)*four-r3(6,14)*two-r3(6,17)*four-r3(6,18)*two-r3(6,21)*p12 &
&                  -r3(6,22)*six-r3(6,25)*p12-r3(6,26)*six+r2(5,3)*two+r2(5,4)+r2(5,7)*two &
&                  +r2(5,8)+r2(5,11)*eight+r2(5,12)*four+r2(5,15)*two+r2(5,16)+r2(5,19)*six &
&                  +r2(5,20)*three)*zz+(+r5(21,23)*two+r5(21,24)-r4(15,32)*four-r4(15,33)*two &
&                  -r4(15,36)*four-r4(15,37)*two+r3(10,29)*two+r3(10,30)+r3(10,33)*eight &
&                  +r3(10,34)*four+r3(10,37)*two+r3(10,38)+r3(10,41)*p12+r3(10,42)*six &
&                  -r2(3,23)*four-r2(3,24)*two-r2(3,27)*four-r2(3,28)*two-r2(3,31)*p12 &
&                  -r2(3,32)*six-r2(3,35)*p12-r2(3,36)*six+r1(3,3)*two+r1(3,4)+r1(3,7)*two &
&                  +r1(3,8)+r1(3,11)*eight+r1(3,12)*four+r1(3,15)*two+r1(3,16)+r1(3,19)*six &
&                  +r1(3,20)*three)*xzz+(+r5(15,24)-r4(10,33)*two-r4(10,37)*two+r3(6,30) &
&                  +r3(6,34)*four+r3(6,38)+r3(6,42)*six-r2(5,24)*two-r2(5,28)*two &
&                  -r2(5,32)*six-r2(5,36)*six+r1(1,4)+r1(1,8)+r1(1,12)*four+r1(1,16) &
&                  +r1(1,20)*three)*zzz+rxyz(1)*xzzz
      eri(6,3,3,3)=r013+(+r7(35,3)+r7(35,4)*two-r6(27,5)*two-r6(27,6)*four-r6(27,7)*two &
&                  -r6(27,8)*four+r5(20,5)+r5(20,6)*two+r5(20,7)*four+r5(20,8)*eight+r5(20,9) &
&                  +r5(20,10)*two+r5(20,11)*six+r5(20,12)*p12+r5(20,21)+r5(20,22)*two &
&                  -r4(14,6)*two-r4(14,7)*four-r4(14,8)*two-r4(14,9)*four-r4(14,10)*six &
&                  -r4(14,11)*p12-r4(14,12)*six-r4(14,13)*p12-r4(14,30)*two-r4(14,31)*four &
&                  -r4(14,34)*two-r4(14,35)*four+r3(9,1)+r3(9,2)*two+r3(9,3)+r3(9,4)*two &
&                  +r3(9,5)*four+r3(9,6)*eight+r3(9,7)+r3(9,8)*two+r3(9,9)*three+r3(9,10)*six &
&                  +r3(9,27)+r3(9,28)*two+r3(9,31)*four+r3(9,32)*eight+r3(9,35)+r3(9,36)*two &
&                  +r3(9,39)*six+r3(9,40)*p12-r2(6,21)*two-r2(6,22)*four-r2(6,25)*two &
&                  -r2(6,26)*four-r2(6,29)*six-r2(6,30)*p12-r2(6,33)*six-r2(6,34)*p12+r1(2,1) &
&                  +r1(2,2)*two+r1(2,5)+r1(2,6)*two+r1(2,9)*four+r1(2,10)*eight+r1(2,13) &
&                  +r1(2,14)*two+r1(2,17)*three+r1(2,18)*six)*qz+(+r6(27,11)*two+r6(27,12) &
&                  -r5(20,15)*four-r5(20,16)*two-r5(20,19)*four-r5(20,20)*two+r4(14,16)*two &
&                  +r4(14,17)+r4(14,20)*eight+r4(14,21)*four+r4(14,24)*two+r4(14,25) &
&                  +r4(14,28)*p12+r4(14,29)*six-r3(9,13)*four-r3(9,14)*two-r3(9,17)*four &
&                  -r3(9,18)*two-r3(9,21)*p12-r3(9,22)*six-r3(9,25)*p12-r3(9,26)*six &
&                  +r2(6,3)*two+r2(6,4)+r2(6,7)*two+r2(6,8)+r2(6,11)*eight+r2(6,12)*four &
&                  +r2(6,15)*two+r2(6,16)+r2(6,19)*six+r2(6,20)*three)*zz+rxyz(2)*zzz
      eri(1,4,3,3)=r310+(+r7(20,3)*two+r7(20,4)-r6(14,5)*four-r6(14,6)*two-r6(14,7)*four &
&                  -r6(14,8)*two+r5(9,5)*two+r5(9,6)+r5(9,7)*eight+r5(9,8)*four+r5(9,9)*two &
&                  +r5(9,10)+r5(9,11)*p12+r5(9,12)*six+r5(20,21)*two+r5(20,22)-r4(5,6)*four &
&                  -r4(5,7)*two-r4(5,8)*four-r4(5,9)*two-r4(5,10)*p12-r4(5,11)*six &
&                  -r4(5,12)*p12-r4(5,13)*six-r4(14,30)*four-r4(14,31)*two-r4(14,34)*four &
&                  -r4(14,35)*two+r3(2,1)*two+r3(2,2)+r3(2,3)*two+r3(2,4)+r3(2,5)*eight &
&                  +r3(2,6)*four+r3(2,7)*two+r3(2,8)+r3(2,9)*six+r3(2,10)*three+r3(9,27)*two &
&                  +r3(9,28)+r3(9,31)*eight+r3(9,32)*four+r3(9,35)*two+r3(9,36)+r3(9,39)*p12 &
&                  +r3(9,40)*six-r2(6,21)*four-r2(6,22)*two-r2(6,25)*four-r2(6,26)*two &
&                  -r2(6,29)*p12-r2(6,30)*six-r2(6,33)*p12-r2(6,34)*six+r1(2,1)*two+r1(2,2) &
&                  +r1(2,5)*two+r1(2,6)+r1(2,9)*eight+r1(2,10)*four+r1(2,13)*two+r1(2,14) &
&                  +r1(2,17)*six+r1(2,18)*three)*qx+(+r6(20,10)+r6(20,11)*two-r5(14,14)*two &
&                  -r5(14,15)*four-r5(14,18)*two-r5(14,19)*four+r4(9,15)+r4(9,16)*two &
&                  +r4(9,19)*four+r4(9,20)*eight+r4(9,23)+r4(9,24)*two+r4(9,27)*six &
&                  +r4(9,28)*p12-r3(5,12)*two-r3(5,13)*four-r3(5,16)*two-r3(5,17)*four &
&                  -r3(5,20)*six-r3(5,21)*p12-r3(5,24)*six-r3(5,25)*p12+r2(4,2)+r2(4,3)*two &
&                  +r2(4,6)+r2(4,7)*two+r2(4,10)*four+r2(4,11)*eight+r2(4,14)+r2(4,15)*two &
&                  +r2(4,18)*three+r2(4,19)*six)*xx+rxyz(3)*xxx
      eri(2,4,3,3)=r130+rxyz(8)*qx
      eri(3,4,3,3)=r112+rxyz(19)*qx+(+r7(27,3)*two-r6(20,5)*four-r6(20,7)*four &
&                  +r5(14,5)*two+r5(14,7)*eight+r5(14,9)*two+r5(14,11)*p12-r4(9,6)*four &
&                  -r4(9,8)*four-r4(9,10)*p12-r4(9,12)*p12+r3(5,1)*two+r3(5,3)*two &
&                  +r3(5,5)*eight+r3(5,7)*two+r3(5,9)*six)*qz+(+r6(27,11)*two-r5(20,15)*four &
&                  -r5(20,19)*four+r4(14,16)*two+r4(14,20)*eight+r4(14,24)*two+r4(14,28)*p12 &
&                  -r3(9,13)*four-r3(9,17)*four-r3(9,21)*p12-r3(9,25)*p12+r2(6,3)*two &
&                  +r2(6,7)*two+r2(6,11)*eight+r2(6,15)*two+r2(6,19)*six)*xz+(+r6(20,10) &
&                  -r5(14,14)*two-r5(14,18)*two+r4(9,15)+r4(9,19)*four+r4(9,23)+r4(9,27)*six &
&                  -r3(5,12)*two-r3(5,16)*two-r3(5,20)*six-r3(5,24)*six+r2(4,2)+r2(4,6) &
&                  +r2(4,10)*four+r2(4,14)+r2(4,18)*three)*zz+rxyz(3)*xzz
      eri(4,4,3,3)=r220+(+r7(26,3)+r7(26,4)-r6(19,5)*two-r6(19,6)*two-r6(19,7)*two &
&                  -r6(19,8)*two+r5(13,5)+r5(13,6)+r5(13,7)*four+r5(13,8)*four+r5(13,9) &
&                  +r5(13,10)+r5(13,11)*six+r5(13,12)*six+r5(15,21)+r5(15,22)-r4(8,6)*two &
&                  -r4(8,7)*two-r4(8,8)*two-r4(8,9)*two-r4(8,10)*six-r4(8,11)*six &
&                  -r4(8,12)*six-r4(8,13)*six-r4(10,30)*two-r4(10,31)*two-r4(10,34)*two &
&                  -r4(10,35)*two+r3(4,1)+r3(4,2)+r3(4,3)+r3(4,4)+r3(4,5)*four+r3(4,6)*four &
&                  +r3(4,7)+r3(4,8)+r3(4,9)*three+r3(4,10)*three+r3(6,27)+r3(6,28) &
&                  +r3(6,31)*four+r3(6,32)*four+r3(6,35)+r3(6,36)+r3(6,39)*six+r3(6,40)*six &
&                  -r2(5,21)*two-r2(5,22)*two-r2(5,25)*two-r2(5,26)*two-r2(5,29)*six &
&                  -r2(5,30)*six-r2(5,33)*six-r2(5,34)*six+r1(1,1)+r1(1,2)+r1(1,5)+r1(1,6) &
&                  +r1(1,9)*four+r1(1,10)*four+r1(1,13)+r1(1,14)+r1(1,17)*three &
&                  +r1(1,18)*three)*qx+rxyz(6)*xx
      eri(5,4,3,3)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(20,10)+r6(20,11) &
&                  -r5(14,14)*two-r5(14,15)*two-r5(14,18)*two-r5(14,19)*two+r4(9,15)+r4(9,16) &
&                  +r4(9,19)*four+r4(9,20)*four+r4(9,23)+r4(9,24)+r4(9,27)*six+r4(9,28)*six &
&                  -r3(5,12)*two-r3(5,13)*two-r3(5,16)*two-r3(5,17)*two-r3(5,20)*six &
&                  -r3(5,21)*six-r3(5,24)*six-r3(5,25)*six+r2(4,2)+r2(4,3)+r2(4,6)+r2(4,7) &
&                  +r2(4,10)*four+r2(4,11)*four+r2(4,14)+r2(4,15)+r2(4,18)*three &
&                  +r2(4,19)*three)*xz+rxyz(3)*xxz
      eri(6,4,3,3)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,3,3)=r301+(+r7(21,3)*two+r7(21,4)-r6(15,5)*four-r6(15,6)*two-r6(15,7)*four &
&                  -r6(15,8)*two+r5(10,5)*two+r5(10,6)+r5(10,7)*eight+r5(10,8)*four &
&                  +r5(10,9)*two+r5(10,10)+r5(10,11)*p12+r5(10,12)*six+r5(21,21)*two &
&                  +r5(21,22)-r4(6,6)*four-r4(6,7)*two-r4(6,8)*four-r4(6,9)*two-r4(6,10)*p12 &
&                  -r4(6,11)*six-r4(6,12)*p12-r4(6,13)*six-r4(15,30)*four-r4(15,31)*two &
&                  -r4(15,34)*four-r4(15,35)*two+r3(3,1)*two+r3(3,2)+r3(3,3)*two+r3(3,4) &
&                  +r3(3,5)*eight+r3(3,6)*four+r3(3,7)*two+r3(3,8)+r3(3,9)*six+r3(3,10)*three &
&                  +r3(10,27)*two+r3(10,28)+r3(10,31)*eight+r3(10,32)*four+r3(10,35)*two &
&                  +r3(10,36)+r3(10,39)*p12+r3(10,40)*six-r2(3,21)*four-r2(3,22)*two &
&                  -r2(3,25)*four-r2(3,26)*two-r2(3,29)*p12-r2(3,30)*six-r2(3,33)*p12 &
&                  -r2(3,34)*six+r1(3,1)*two+r1(3,2)+r1(3,5)*two+r1(3,6)+r1(3,9)*eight &
&                  +r1(3,10)*four+r1(3,13)*two+r1(3,14)+r1(3,17)*six+r1(3,18)*three)*qx+( &
&                  +r7(15,4)-r6(10,6)*two-r6(10,8)*two+r5(6,6)+r5(6,8)*four+r5(6,10) &
&                  +r5(6,12)*six+r5(15,22)*three-r4(3,7)*two-r4(3,9)*two-r4(3,11)*six &
&                  -r4(3,13)*six-r4(10,31)*six-r4(10,35)*six+r3(1,2)+r3(1,4)+r3(1,6)*four &
&                  +r3(1,8)+r3(1,10)*three+r3(6,28)*three+r3(6,32)*p12+r3(6,36)*three &
&                  +r3(6,40)*p18-r2(5,22)*six-r2(5,26)*six-r2(5,30)*p18-r2(5,34)*p18 &
&                  +r1(1,2)*three+r1(1,6)*three+r1(1,10)*p12+r1(1,14)*three+r1(1,18)*nine)*qz &
&                  +(+r6(21,10)+r6(21,11)*two-r5(15,14)*two-r5(15,15)*four-r5(15,18)*two &
&                  -r5(15,19)*four+r4(10,15)+r4(10,16)*two+r4(10,19)*four+r4(10,20)*eight &
&                  +r4(10,23)+r4(10,24)*two+r4(10,27)*six+r4(10,28)*p12-r3(6,12)*two &
&                  -r3(6,13)*four-r3(6,16)*two-r3(6,17)*four-r3(6,20)*six-r3(6,21)*p12 &
&                  -r3(6,24)*six-r3(6,25)*p12+r2(5,2)+r2(5,3)*two+r2(5,6)+r2(5,7)*two &
&                  +r2(5,10)*four+r2(5,11)*eight+r2(5,14)+r2(5,15)*two+r2(5,18)*three &
&                  +r2(5,19)*six)*xx+(+r6(15,11)*two+r6(15,12)-r5(10,15)*four-r5(10,16)*two &
&                  -r5(10,19)*four-r5(10,20)*two+r4(6,16)*two+r4(6,17)+r4(6,20)*eight &
&                  +r4(6,21)*four+r4(6,24)*two+r4(6,25)+r4(6,28)*p12+r4(6,29)*six &
&                  +r4(15,40)*two+r4(15,41)-r3(3,13)*four-r3(3,14)*two-r3(3,17)*four &
&                  -r3(3,18)*two-r3(3,21)*p12-r3(3,22)*six-r3(3,25)*p12-r3(3,26)*six &
&                  -r3(10,45)*four-r3(10,46)*two-r3(10,50)*four-r3(10,51)*two+r2(1,3)*two &
&                  +r2(1,4)+r2(1,7)*two+r2(1,8)+r2(1,11)*eight+r2(1,12)*four+r2(1,15)*two &
&                  +r2(1,16)+r2(1,19)*six+r2(1,20)*three+r2(3,39)*two+r2(3,40)+r2(3,44)*eight &
&                  +r2(3,45)*four+r2(3,49)*two+r2(3,50)+r2(3,54)*p12+r2(3,55)*six &
&                  -r1(3,23)*four-r1(3,24)*two-r1(3,28)*four-r1(3,29)*two-r1(3,33)*p12 &
&                  -r1(3,34)*six-r1(3,38)*p12-r1(3,39)*six+r0(3)*two+r0(4)+r0(8)*two+r0(9) &
&                  +r0(13)*eight+r0(14)*four+r0(18)*two+r0(19)+r0(23)*six+r0(24)*three)*xz+( &
&                  +r5(21,23)-r4(15,32)*two-r4(15,36)*two+r3(10,29)+r3(10,33)*four+r3(10,37) &
&                  +r3(10,41)*six-r2(3,23)*two-r2(3,27)*two-r2(3,31)*six-r2(3,35)*six+r1(3,3) &
&                  +r1(3,7)+r1(3,11)*four+r1(3,15)+r1(3,19)*three)*xxx+(+r5(15,23) &
&                  +r5(15,24)*two-r4(10,32)*two-r4(10,33)*four-r4(10,36)*two-r4(10,37)*four &
&                  +r3(6,29)+r3(6,30)*two+r3(6,33)*four+r3(6,34)*eight+r3(6,37)+r3(6,38)*two &
&                  +r3(6,41)*six+r3(6,42)*p12-r2(5,23)*two-r2(5,24)*four-r2(5,27)*two &
&                  -r2(5,28)*four-r2(5,31)*six-r2(5,32)*p12-r2(5,35)*six-r2(5,36)*p12+r1(1,3) &
&                  +r1(1,4)*two+r1(1,7)+r1(1,8)*two+r1(1,11)*four+r1(1,12)*eight+r1(1,15) &
&                  +r1(1,16)*two+r1(1,19)*three+r1(1,20)*six)*xxz+rxyz(1)*xxxz
      eri(2,5,3,3)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,3,3)=r103+(+r7(36,4)-r6(28,6)*two-r6(28,8)*two+r5(21,6)+r5(21,8)*four &
&                  +r5(21,10)+r5(21,12)*six+r5(21,22)*three-r4(15,7)*two-r4(15,9)*two &
&                  -r4(15,11)*six-r4(15,13)*six-r4(15,31)*six-r4(15,35)*six+r3(10,2)+r3(10,4) &
&                  +r3(10,6)*four+r3(10,8)+r3(10,10)*three+r3(10,28)*three+r3(10,32)*p12 &
&                  +r3(10,36)*three+r3(10,40)*p18-r2(3,22)*six-r2(3,26)*six-r2(3,30)*p18 &
&                  -r2(3,34)*p18+r1(3,2)*three+r1(3,6)*three+r1(3,10)*p12+r1(3,14)*three &
&                  +r1(3,18)*nine)*qx+(+r7(28,3)*two+r7(28,4)-r6(21,5)*four-r6(21,6)*two &
&                  -r6(21,7)*four-r6(21,8)*two+r5(15,5)*two+r5(15,6)+r5(15,7)*eight &
&                  +r5(15,8)*four+r5(15,9)*two+r5(15,10)+r5(15,11)*p12+r5(15,12)*six &
&                  +r5(15,21)*two+r5(15,22)-r4(10,6)*four-r4(10,7)*two-r4(10,8)*four &
&                  -r4(10,9)*two-r4(10,10)*p12-r4(10,11)*six-r4(10,12)*p12-r4(10,13)*six &
&                  -r4(10,30)*four-r4(10,31)*two-r4(10,34)*four-r4(10,35)*two+r3(6,1)*two &
&                  +r3(6,2)+r3(6,3)*two+r3(6,4)+r3(6,5)*eight+r3(6,6)*four+r3(6,7)*two &
&                  +r3(6,8)+r3(6,9)*six+r3(6,10)*three+r3(6,27)*two+r3(6,28)+r3(6,31)*eight &
&                  +r3(6,32)*four+r3(6,35)*two+r3(6,36)+r3(6,39)*p12+r3(6,40)*six &
&                  -r2(5,21)*four-r2(5,22)*two-r2(5,25)*four-r2(5,26)*two-r2(5,29)*p12 &
&                  -r2(5,30)*six-r2(5,33)*p12-r2(5,34)*six+r1(1,1)*two+r1(1,2)+r1(1,5)*two &
&                  +r1(1,6)+r1(1,9)*eight+r1(1,10)*four+r1(1,13)*two+r1(1,14)+r1(1,17)*six &
&                  +r1(1,18)*three)*qz+(+r6(28,11)*two+r6(28,12)-r5(21,15)*four-r5(21,16)*two &
&                  -r5(21,19)*four-r5(21,20)*two+r4(15,16)*two+r4(15,17)+r4(15,20)*eight &
&                  +r4(15,21)*four+r4(15,24)*two+r4(15,25)+r4(15,28)*p12+r4(15,29)*six &
&                  +r4(15,40)*two+r4(15,41)-r3(10,13)*four-r3(10,14)*two-r3(10,17)*four &
&                  -r3(10,18)*two-r3(10,21)*p12-r3(10,22)*six-r3(10,25)*p12-r3(10,26)*six &
&                  -r3(10,45)*four-r3(10,46)*two-r3(10,50)*four-r3(10,51)*two+r2(3,3)*two &
&                  +r2(3,4)+r2(3,7)*two+r2(3,8)+r2(3,11)*eight+r2(3,12)*four+r2(3,15)*two &
&                  +r2(3,16)+r2(3,19)*six+r2(3,20)*three+r2(3,39)*two+r2(3,40)+r2(3,44)*eight &
&                  +r2(3,45)*four+r2(3,49)*two+r2(3,50)+r2(3,54)*p12+r2(3,55)*six &
&                  -r1(3,23)*four-r1(3,24)*two-r1(3,28)*four-r1(3,29)*two-r1(3,33)*p12 &
&                  -r1(3,34)*six-r1(3,38)*p12-r1(3,39)*six+r0(3)*two+r0(4)+r0(8)*two+r0(9) &
&                  +r0(13)*eight+r0(14)*four+r0(18)*two+r0(19)+r0(23)*six+r0(24)*three)*xz+( &
&                  +r6(21,10)+r6(21,11)*two-r5(15,14)*two-r5(15,15)*four-r5(15,18)*two &
&                  -r5(15,19)*four+r4(10,15)+r4(10,16)*two+r4(10,19)*four+r4(10,20)*eight &
&                  +r4(10,23)+r4(10,24)*two+r4(10,27)*six+r4(10,28)*p12-r3(6,12)*two &
&                  -r3(6,13)*four-r3(6,16)*two-r3(6,17)*four-r3(6,20)*six-r3(6,21)*p12 &
&                  -r3(6,24)*six-r3(6,25)*p12+r2(5,2)+r2(5,3)*two+r2(5,6)+r2(5,7)*two &
&                  +r2(5,10)*four+r2(5,11)*eight+r2(5,14)+r2(5,15)*two+r2(5,18)*three &
&                  +r2(5,19)*six)*zz+(+r5(21,23)+r5(21,24)*two-r4(15,32)*two-r4(15,33)*four &
&                  -r4(15,36)*two-r4(15,37)*four+r3(10,29)+r3(10,30)*two+r3(10,33)*four &
&                  +r3(10,34)*eight+r3(10,37)+r3(10,38)*two+r3(10,41)*six+r3(10,42)*p12 &
&                  -r2(3,23)*two-r2(3,24)*four-r2(3,27)*two-r2(3,28)*four-r2(3,31)*six &
&                  -r2(3,32)*p12-r2(3,35)*six-r2(3,36)*p12+r1(3,3)+r1(3,4)*two+r1(3,7) &
&                  +r1(3,8)*two+r1(3,11)*four+r1(3,12)*eight+r1(3,15)+r1(3,16)*two &
&                  +r1(3,19)*three+r1(3,20)*six)*xzz+(+r5(15,23)-r4(10,32)*two-r4(10,36)*two &
&                  +r3(6,29)+r3(6,33)*four+r3(6,37)+r3(6,41)*six-r2(5,23)*two-r2(5,27)*two &
&                  -r2(5,31)*six-r2(5,35)*six+r1(1,3)+r1(1,7)+r1(1,11)*four+r1(1,15) &
&                  +r1(1,19)*three)*zzz+rxyz(1)*xzzz
      eri(4,5,3,3)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(20,11)+r6(20,12) &
&                  -r5(14,15)*two-r5(14,16)*two-r5(14,19)*two-r5(14,20)*two+r4(9,16)+r4(9,17) &
&                  +r4(9,20)*four+r4(9,21)*four+r4(9,24)+r4(9,25)+r4(9,28)*six+r4(9,29)*six &
&                  -r3(5,13)*two-r3(5,14)*two-r3(5,17)*two-r3(5,18)*two-r3(5,21)*six &
&                  -r3(5,22)*six-r3(5,25)*six-r3(5,26)*six+r2(4,3)+r2(4,4)+r2(4,7)+r2(4,8) &
&                  +r2(4,11)*four+r2(4,12)*four+r2(4,15)+r2(4,16)+r2(4,19)*three &
&                  +r2(4,20)*three)*xz+rxyz(2)*xxz
      eri(5,5,3,3)=r202+(+r7(28,3)+r7(28,4)-r6(21,5)*two-r6(21,6)*two-r6(21,7)*two &
&                  -r6(21,8)*two+r5(15,5)+r5(15,6)+r5(15,7)*four+r5(15,8)*four+r5(15,9) &
&                  +r5(15,10)+r5(15,11)*six+r5(15,12)*six+r5(15,21)+r5(15,22)-r4(10,6)*two &
&                  -r4(10,7)*two-r4(10,8)*two-r4(10,9)*two-r4(10,10)*six-r4(10,11)*six &
&                  -r4(10,12)*six-r4(10,13)*six-r4(10,30)*two-r4(10,31)*two-r4(10,34)*two &
&                  -r4(10,35)*two+r3(6,1)+r3(6,2)+r3(6,3)+r3(6,4)+r3(6,5)*four+r3(6,6)*four &
&                  +r3(6,7)+r3(6,8)+r3(6,9)*three+r3(6,10)*three+r3(6,27)+r3(6,28) &
&                  +r3(6,31)*four+r3(6,32)*four+r3(6,35)+r3(6,36)+r3(6,39)*six+r3(6,40)*six &
&                  -r2(5,21)*two-r2(5,22)*two-r2(5,25)*two-r2(5,26)*two-r2(5,29)*six &
&                  -r2(5,30)*six-r2(5,33)*six-r2(5,34)*six+r1(1,1)+r1(1,2)+r1(1,5)+r1(1,6) &
&                  +r1(1,9)*four+r1(1,10)*four+r1(1,13)+r1(1,14)+r1(1,17)*three &
&                  +r1(1,18)*three)*qx+(+r7(21,3)+r7(21,4)-r6(15,5)*two-r6(15,6)*two &
&                  -r6(15,7)*two-r6(15,8)*two+r5(10,5)+r5(10,6)+r5(10,7)*four+r5(10,8)*four &
&                  +r5(10,9)+r5(10,10)+r5(10,11)*six+r5(10,12)*six+r5(21,21)+r5(21,22) &
&                  -r4(6,6)*two-r4(6,7)*two-r4(6,8)*two-r4(6,9)*two-r4(6,10)*six-r4(6,11)*six &
&                  -r4(6,12)*six-r4(6,13)*six-r4(15,30)*two-r4(15,31)*two-r4(15,34)*two &
&                  -r4(15,35)*two+r3(3,1)+r3(3,2)+r3(3,3)+r3(3,4)+r3(3,5)*four+r3(3,6)*four &
&                  +r3(3,7)+r3(3,8)+r3(3,9)*three+r3(3,10)*three+r3(10,27)+r3(10,28) &
&                  +r3(10,31)*four+r3(10,32)*four+r3(10,35)+r3(10,36)+r3(10,39)*six &
&                  +r3(10,40)*six-r2(3,21)*two-r2(3,22)*two-r2(3,25)*two-r2(3,26)*two &
&                  -r2(3,29)*six-r2(3,30)*six-r2(3,33)*six-r2(3,34)*six+r1(3,1)+r1(3,2) &
&                  +r1(3,5)+r1(3,6)+r1(3,9)*four+r1(3,10)*four+r1(3,13)+r1(3,14) &
&                  +r1(3,17)*three+r1(3,18)*three)*qz+(+r6(28,11)-r5(21,15)*two-r5(21,19)*two &
&                  +r4(15,16)+r4(15,20)*four+r4(15,24)+r4(15,28)*six+r4(15,40)-r3(10,13)*two &
&                  -r3(10,17)*two-r3(10,21)*six-r3(10,25)*six-r3(10,45)*two-r3(10,50)*two &
&                  +r2(3,3)+r2(3,7)+r2(3,11)*four+r2(3,15)+r2(3,19)*three+r2(3,39) &
&                  +r2(3,44)*four+r2(3,49)+r2(3,54)*six-r1(3,23)*two-r1(3,28)*two &
&                  -r1(3,33)*six-r1(3,38)*six+r0(3)+r0(8)+r0(13)*four+r0(18)+r0(23)*three)*xx &
&                  +(+r6(21,10)+r6(21,11)*two+r6(21,12)-r5(15,14)*two-r5(15,15)*four &
&                  -r5(15,16)*two-r5(15,18)*two-r5(15,19)*four-r5(15,20)*two+r4(10,15) &
&                  +r4(10,16)*two+r4(10,17)+r4(10,19)*four+r4(10,20)*eight+r4(10,21)*four &
&                  +r4(10,23)+r4(10,24)*two+r4(10,25)+r4(10,27)*six+r4(10,28)*p12 &
&                  +r4(10,29)*six-r3(6,12)*two-r3(6,13)*four-r3(6,14)*two-r3(6,16)*two &
&                  -r3(6,17)*four-r3(6,18)*two-r3(6,20)*six-r3(6,21)*p12-r3(6,22)*six &
&                  -r3(6,24)*six-r3(6,25)*p12-r3(6,26)*six+r2(5,2)+r2(5,3)*two+r2(5,4) &
&                  +r2(5,6)+r2(5,7)*two+r2(5,8)+r2(5,10)*four+r2(5,11)*eight+r2(5,12)*four &
&                  +r2(5,14)+r2(5,15)*two+r2(5,16)+r2(5,18)*three+r2(5,19)*six+r2(5,20)*three &
&                  )*xz+(+r6(15,11)-r5(10,15)*two-r5(10,19)*two+r4(6,16)+r4(6,20)*four &
&                  +r4(6,24)+r4(6,28)*six+r4(15,40)-r3(3,13)*two-r3(3,17)*two-r3(3,21)*six &
&                  -r3(3,25)*six-r3(10,45)*two-r3(10,50)*two+r2(1,3)+r2(1,7)+r2(1,11)*four &
&                  +r2(1,15)+r2(1,19)*three+r2(3,39)+r2(3,44)*four+r2(3,49)+r2(3,54)*six &
&                  -r1(3,23)*two-r1(3,28)*two-r1(3,33)*six-r1(3,38)*six+r0(3)+r0(8) &
&                  +r0(13)*four+r0(18)+r0(23)*three)*zz+(+r5(21,23)+r5(21,24)-r4(15,32)*two &
&                  -r4(15,33)*two-r4(15,36)*two-r4(15,37)*two+r3(10,29)+r3(10,30) &
&                  +r3(10,33)*four+r3(10,34)*four+r3(10,37)+r3(10,38)+r3(10,41)*six &
&                  +r3(10,42)*six-r2(3,23)*two-r2(3,24)*two-r2(3,27)*two-r2(3,28)*two &
&                  -r2(3,31)*six-r2(3,32)*six-r2(3,35)*six-r2(3,36)*six+r1(3,3)+r1(3,4) &
&                  +r1(3,7)+r1(3,8)+r1(3,11)*four+r1(3,12)*four+r1(3,15)+r1(3,16) &
&                  +r1(3,19)*three+r1(3,20)*three)*xxz+(+r5(15,23)+r5(15,24)-r4(10,32)*two &
&                  -r4(10,33)*two-r4(10,36)*two-r4(10,37)*two+r3(6,29)+r3(6,30)+r3(6,33)*four &
&                  +r3(6,34)*four+r3(6,37)+r3(6,38)+r3(6,41)*six+r3(6,42)*six-r2(5,23)*two &
&                  -r2(5,24)*two-r2(5,27)*two-r2(5,28)*two-r2(5,31)*six-r2(5,32)*six &
&                  -r2(5,35)*six-r2(5,36)*six+r1(1,3)+r1(1,4)+r1(1,7)+r1(1,8)+r1(1,11)*four &
&                  +r1(1,12)*four+r1(1,15)+r1(1,16)+r1(1,19)*three+r1(1,20)*three)*xzz &
&                  +rxyz(1)*xxzz
      eri(6,5,3,3)=r112+rxyz(19)*qx+(+r7(27,3)+r7(27,4)-r6(20,5)*two-r6(20,6)*two &
&                  -r6(20,7)*two-r6(20,8)*two+r5(14,5)+r5(14,6)+r5(14,7)*four+r5(14,8)*four &
&                  +r5(14,9)+r5(14,10)+r5(14,11)*six+r5(14,12)*six-r4(9,6)*two-r4(9,7)*two &
&                  -r4(9,8)*two-r4(9,9)*two-r4(9,10)*six-r4(9,11)*six-r4(9,12)*six &
&                  -r4(9,13)*six+r3(5,1)+r3(5,2)+r3(5,3)+r3(5,4)+r3(5,5)*four+r3(5,6)*four &
&                  +r3(5,7)+r3(5,8)+r3(5,9)*three+r3(5,10)*three)*qz+(+r6(27,11)+r6(27,12) &
&                  -r5(20,15)*two-r5(20,16)*two-r5(20,19)*two-r5(20,20)*two+r4(14,16) &
&                  +r4(14,17)+r4(14,20)*four+r4(14,21)*four+r4(14,24)+r4(14,25)+r4(14,28)*six &
&                  +r4(14,29)*six-r3(9,13)*two-r3(9,14)*two-r3(9,17)*two-r3(9,18)*two &
&                  -r3(9,21)*six-r3(9,22)*six-r3(9,25)*six-r3(9,26)*six+r2(6,3)+r2(6,4) &
&                  +r2(6,7)+r2(6,8)+r2(6,11)*four+r2(6,12)*four+r2(6,15)+r2(6,16) &
&                  +r2(6,19)*three+r2(6,20)*three)*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,3,3)=r211+(+r7(27,3)*two-r6(20,5)*four-r6(20,7)*four+r5(14,5)*two &
&                  +r5(14,7)*eight+r5(14,9)*two+r5(14,11)*p12-r4(9,6)*four-r4(9,8)*four &
&                  -r4(9,10)*p12-r4(9,12)*p12+r3(5,1)*two+r3(5,3)*two+r3(5,5)*eight &
&                  +r3(5,7)*two+r3(5,9)*six)*qx+rxyz(16)*qz+(+r6(27,10)-r5(20,14)*two &
&                  -r5(20,18)*two+r4(14,15)+r4(14,19)*four+r4(14,23)+r4(14,27)*six &
&                  -r3(9,12)*two-r3(9,16)*two-r3(9,20)*six-r3(9,24)*six+r2(6,2)+r2(6,6) &
&                  +r2(6,10)*four+r2(6,14)+r2(6,18)*three)*xx+(+r6(20,11)*two-r5(14,15)*four &
&                  -r5(14,19)*four+r4(9,16)*two+r4(9,20)*eight+r4(9,24)*two+r4(9,28)*p12 &
&                  -r3(5,13)*four-r3(5,17)*four-r3(5,21)*p12-r3(5,25)*p12+r2(4,3)*two &
&                  +r2(4,7)*two+r2(4,11)*eight+r2(4,15)*two+r2(4,19)*six)*xz+rxyz(3)*xxz
      eri(2,6,3,3)=r031+rxyz(8)*qz
      eri(3,6,3,3)=r013+(+r7(35,3)*two+r7(35,4)-r6(27,5)*four-r6(27,6)*two-r6(27,7)*four &
&                  -r6(27,8)*two+r5(20,5)*two+r5(20,6)+r5(20,7)*eight+r5(20,8)*four &
&                  +r5(20,9)*two+r5(20,10)+r5(20,11)*p12+r5(20,12)*six+r5(20,21)*two &
&                  +r5(20,22)-r4(14,6)*four-r4(14,7)*two-r4(14,8)*four-r4(14,9)*two &
&                  -r4(14,10)*p12-r4(14,11)*six-r4(14,12)*p12-r4(14,13)*six-r4(14,30)*four &
&                  -r4(14,31)*two-r4(14,34)*four-r4(14,35)*two+r3(9,1)*two+r3(9,2) &
&                  +r3(9,3)*two+r3(9,4)+r3(9,5)*eight+r3(9,6)*four+r3(9,7)*two+r3(9,8) &
&                  +r3(9,9)*six+r3(9,10)*three+r3(9,27)*two+r3(9,28)+r3(9,31)*eight &
&                  +r3(9,32)*four+r3(9,35)*two+r3(9,36)+r3(9,39)*p12+r3(9,40)*six &
&                  -r2(6,21)*four-r2(6,22)*two-r2(6,25)*four-r2(6,26)*two-r2(6,29)*p12 &
&                  -r2(6,30)*six-r2(6,33)*p12-r2(6,34)*six+r1(2,1)*two+r1(2,2)+r1(2,5)*two &
&                  +r1(2,6)+r1(2,9)*eight+r1(2,10)*four+r1(2,13)*two+r1(2,14)+r1(2,17)*six &
&                  +r1(2,18)*three)*qz+(+r6(27,10)+r6(27,11)*two-r5(20,14)*two-r5(20,15)*four &
&                  -r5(20,18)*two-r5(20,19)*four+r4(14,15)+r4(14,16)*two+r4(14,19)*four &
&                  +r4(14,20)*eight+r4(14,23)+r4(14,24)*two+r4(14,27)*six+r4(14,28)*p12 &
&                  -r3(9,12)*two-r3(9,13)*four-r3(9,16)*two-r3(9,17)*four-r3(9,20)*six &
&                  -r3(9,21)*p12-r3(9,24)*six-r3(9,25)*p12+r2(6,2)+r2(6,3)*two+r2(6,6) &
&                  +r2(6,7)*two+r2(6,10)*four+r2(6,11)*eight+r2(6,14)+r2(6,15)*two &
&                  +r2(6,18)*three+r2(6,19)*six)*zz+rxyz(3)*zzz
      eri(4,6,3,3)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,3,3)=r112+rxyz(18)*qx+(+r7(27,3)+r7(27,4)-r6(20,5)*two-r6(20,6)*two &
&                  -r6(20,7)*two-r6(20,8)*two+r5(14,5)+r5(14,6)+r5(14,7)*four+r5(14,8)*four &
&                  +r5(14,9)+r5(14,10)+r5(14,11)*six+r5(14,12)*six-r4(9,6)*two-r4(9,7)*two &
&                  -r4(9,8)*two-r4(9,9)*two-r4(9,10)*six-r4(9,11)*six-r4(9,12)*six &
&                  -r4(9,13)*six+r3(5,1)+r3(5,2)+r3(5,3)+r3(5,4)+r3(5,5)*four+r3(5,6)*four &
&                  +r3(5,7)+r3(5,8)+r3(5,9)*three+r3(5,10)*three)*qz+(+r6(27,10)+r6(27,11) &
&                  -r5(20,14)*two-r5(20,15)*two-r5(20,18)*two-r5(20,19)*two+r4(14,15) &
&                  +r4(14,16)+r4(14,19)*four+r4(14,20)*four+r4(14,23)+r4(14,24)+r4(14,27)*six &
&                  +r4(14,28)*six-r3(9,12)*two-r3(9,13)*two-r3(9,16)*two-r3(9,17)*two &
&                  -r3(9,20)*six-r3(9,21)*six-r3(9,24)*six-r3(9,25)*six+r2(6,2)+r2(6,3) &
&                  +r2(6,6)+r2(6,7)+r2(6,10)*four+r2(6,11)*four+r2(6,14)+r2(6,15) &
&                  +r2(6,18)*three+r2(6,19)*three)*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,3,3)=r022+(+r7(34,3)+r7(34,4)-r6(26,5)*two-r6(26,6)*two-r6(26,7)*two &
&                  -r6(26,8)*two+r5(19,5)+r5(19,6)+r5(19,7)*four+r5(19,8)*four+r5(19,9) &
&                  +r5(19,10)+r5(19,11)*six+r5(19,12)*six+r5(21,21)+r5(21,22)-r4(13,6)*two &
&                  -r4(13,7)*two-r4(13,8)*two-r4(13,9)*two-r4(13,10)*six-r4(13,11)*six &
&                  -r4(13,12)*six-r4(13,13)*six-r4(15,30)*two-r4(15,31)*two-r4(15,34)*two &
&                  -r4(15,35)*two+r3(8,1)+r3(8,2)+r3(8,3)+r3(8,4)+r3(8,5)*four+r3(8,6)*four &
&                  +r3(8,7)+r3(8,8)+r3(8,9)*three+r3(8,10)*three+r3(10,27)+r3(10,28) &
&                  +r3(10,31)*four+r3(10,32)*four+r3(10,35)+r3(10,36)+r3(10,39)*six &
&                  +r3(10,40)*six-r2(3,21)*two-r2(3,22)*two-r2(3,25)*two-r2(3,26)*two &
&                  -r2(3,29)*six-r2(3,30)*six-r2(3,33)*six-r2(3,34)*six+r1(3,1)+r1(3,2) &
&                  +r1(3,5)+r1(3,6)+r1(3,9)*four+r1(3,10)*four+r1(3,13)+r1(3,14) &
&                  +r1(3,17)*three+r1(3,18)*three)*qz+rxyz(6)*zz
!
      r400= r8(9)-r7(5,2)*two+r6(2,3)+r6(2,4)+r6(9,9)*six-r5(5,17)*p12+r4(2,22)*six &
&          +r4(2,26)*six+r4(9,38)*three-r3(5,48)*six+r2(4,47)*three+r2(4,52)*three
      r310= r8(13)-r7(8,2)*two+r6(4,3)+r6(4,4)+r6(13,9)*three-r5(8,17)*six+r4(4,22)*three &
&          +r4(4,26)*three
      r301= r8(14)-r7(9,2)*two+r6(5,3)+r6(5,4)+r6(14,9)*three-r5(9,17)*six+r4(5,22)*three &
&          +r4(5,26)*three
      r220= r8(18)-r7(12,2)*two+r6(7,3)+r6(7,4)+r6(9,9)+r6(18,9)-r5(5,17)*two-r5(12,17)*two &
&          +r4(2,22)+r4(7,22)+r4(2,26)+r4(7,26)+r4(9,38)-r3(5,48)*two+r2(4,47)+r2(4,52)
      r211= r8(19)-r7(13,2)*two+r6(8,3)+r6(8,4)+r6(19,9)-r5(13,17)*two+r4(8,22)+r4(8,26)
      r202= r8(20)-r7(14,2)*two+r6(9,3)+r6(9,4)+r6(9,9)+r6(20,9)-r5(5,17)*two-r5(14,17)*two &
&          +r4(2,22)+r4(9,22)+r4(2,26)+r4(9,26)+r4(9,38)-r3(5,48)*two+r2(4,47)+r2(4,52)
      r130= r8(24)-r7(17,2)*two+r6(11,3)+r6(11,4)+r6(13,9)*three-r5(8,17)*six &
&          +r4(4,22)*three+r4(4,26)*three
      r121= r8(25)-r7(18,2)*two+r6(12,3)+r6(12,4)+r6(14,9)-r5(9,17)*two+r4(5,22)+r4(5,26)
      r112= r8(26)-r7(19,2)*two+r6(13,3)+r6(13,4)+r6(13,9)-r5(8,17)*two+r4(4,22)+r4(4,26)
      r103= r8(27)-r7(20,2)*two+r6(14,3)+r6(14,4)+r6(14,9)*three-r5(9,17)*six &
&          +r4(5,22)*three+r4(5,26)*three
      r040= r8(31)-r7(23,2)*two+r6(16,3)+r6(16,4)+r6(18,9)*six-r5(12,17)*p12+r4(7,22)*six &
&          +r4(7,26)*six+r4(9,38)*three-r3(5,48)*six+r2(4,47)*three+r2(4,52)*three
      r031= r8(32)-r7(24,2)*two+r6(17,3)+r6(17,4)+r6(19,9)*three-r5(13,17)*six &
&          +r4(8,22)*three+r4(8,26)*three
      r022= r8(33)-r7(25,2)*two+r6(18,3)+r6(18,4)+r6(18,9)+r6(20,9)-r5(12,17)*two &
&          -r5(14,17)*two+r4(7,22)+r4(9,22)+r4(7,26)+r4(9,26)+r4(9,38)-r3(5,48)*two &
&          +r2(4,47)+r2(4,52)
      r013= r8(34)-r7(26,2)*two+r6(19,3)+r6(19,4)+r6(19,9)*three-r5(13,17)*six &
&          +r4(8,22)*three+r4(8,26)*three
      r004= r8(35)-r7(27,2)*two+r6(20,3)+r6(20,4)+r6(20,9)*six-r5(14,17)*p12+r4(9,22)*six &
&          +r4(9,26)*six+r4(9,38)*three-r3(5,48)*six+r2(4,47)*three+r2(4,52)*three
      rxyz(1)=+r4(9,42)-r3(5,52)*two+r2(4,51)+r2(4,56)
      rxyz(2)=+r5(13,24)-r4(8,37)*two+r3(4,38)+r3(4,42)
      rxyz(3)=+r5(13,23)-r4(8,36)*two+r3(4,37)+r3(4,41)
      rxyz(4)=+r6(18,10)-r5(12,18)*two+r4(7,23)+r4(7,27)+r4(9,39)-r3(5,49)*two+r2(4,48) &
&             +r2(4,53)
      rxyz(5)=+r6(18,12)-r5(12,20)*two+r4(7,25)+r4(7,29)+r4(9,41)-r3(5,51)*two+r2(4,50) &
&             +r2(4,55)
      rxyz(6)=+r6(18,11)-r5(12,19)*two+r4(7,24)+r4(7,28)+r4(9,40)-r3(5,50)*two+r2(4,49) &
&             +r2(4,54)
      rxyz(7)=+r7(24,3)-r6(17,7)*two+r5(11,9)+r5(11,11)+r5(13,21)*three-r4(8,34)*six &
&             +r3(4,35)*three+r3(4,39)*three
      rxyz(8)=+r7(24,4)-r6(17,8)*two+r5(11,10)+r5(11,12)+r5(13,22)*three-r4(8,35)*six &
&             +r3(4,36)*three+r3(4,40)*three
      rxyz(9)=+r7(19,3)+r7(19,4)-r6(13,7)*two-r6(13,8)*two+r5(8,9)+r5(8,10)+r5(8,11) &
&             +r5(8,12)
      rxyz(10)=+r6(19,11)-r5(13,19)*two+r4(8,24)+r4(8,28)
      rxyz(11)=+r6(13,11)-r5(8,19)*two+r4(4,24)+r4(4,28)
      rxyz(12)=+r7(25,3)-r6(18,7)*two+r5(12,9)+r5(12,11)+r5(14,21)-r4(9,34)*two+r3(5,35) &
&             +r3(5,39)
      rxyz(13)=+r7(25,4)-r6(18,8)*two+r5(12,10)+r5(12,12)+r5(14,22)-r4(9,35)*two+r3(5,36) &
&             +r3(5,40)
      rxyz(14)=+r7(18,4)-r6(12,8)*two+r5(7,10)+r5(7,12)+r5(9,22)-r4(5,35)*two+r3(2,36) &
&             +r3(2,40)
      rxyz(15)=+r7(18,3)-r6(12,7)*two+r5(7,9)+r5(7,11)+r5(9,21)-r4(5,34)*two+r3(2,35) &
&             +r3(2,39)
      rxyz(16)=+r7(13,4)-r6(8,8)*two+r5(4,10)+r5(4,12)+r5(13,22)-r4(8,35)*two+r3(4,36) &
&             +r3(4,40)
      rxyz(17)=+r7(13,3)-r6(8,7)*two+r5(4,9)+r5(4,11)+r5(13,21)-r4(8,34)*two+r3(4,35) &
&             +r3(4,39)
      rxyz(18)=+r7(26,3)-r6(19,7)*two+r5(13,9)+r5(13,11)+r5(13,21)-r4(8,34)*two+r3(4,35) &
&             +r3(4,39)
      rxyz(19)=+r7(26,4)-r6(19,8)*two+r5(13,10)+r5(13,12)+r5(13,22)-r4(8,35)*two+r3(4,36) &
&             +r3(4,40)
      rxyz(20)=+r6(14,11)*four-r5(9,19)*eight+r4(5,24)*four+r4(5,28)*four
      eri(1,1,4,3)=r400+(+r7(9,3)*two+r7(9,4)*two-r6(5,7)*four-r6(5,8)*four+r5(2,9)*two &
&                  +r5(2,10)*two+r5(2,11)*two+r5(2,12)*two+r5(9,21)*six+r5(9,22)*six &
&                  -r4(5,34)*p12-r4(5,35)*p12+r3(2,35)*six+r3(2,36)*six+r3(2,39)*six &
&                  +r3(2,40)*six)*qx+(+r6(9,10)+r6(9,11)*four+r6(9,12)-r5(5,18)*two &
&                  -r5(5,19)*eight-r5(5,20)*two+r4(2,23)+r4(2,24)*four+r4(2,25)+r4(2,27) &
&                  +r4(2,28)*four+r4(2,29)+r4(9,39)+r4(9,40)*four+r4(9,41)-r3(5,49)*two &
&                  -r3(5,50)*eight-r3(5,51)*two+r2(4,48)+r2(4,49)*four+r2(4,50)+r2(4,53) &
&                  +r2(4,54)*four+r2(4,55))*xx+(+r5(9,23)*two+r5(9,24)*two-r4(5,36)*four &
&                  -r4(5,37)*four+r3(2,37)*two+r3(2,38)*two+r3(2,41)*two+r3(2,42)*two)*xxx &
&                  +rxyz(1)*xxxx
      eri(2,1,4,3)=r220+(+r7(18,4)*two-r6(12,8)*four+r5(7,10)*two+r5(7,12)*two &
&                  +r5(9,22)*two-r4(5,35)*four+r3(2,36)*two+r3(2,40)*two)*qx+rxyz(5)*xx
      eri(3,1,4,3)=r202+(+r7(20,4)*two-r6(14,8)*four+r5(9,10)*two+r5(9,12)*two &
&                  +r5(9,22)*two-r4(5,35)*four+r3(2,36)*two+r3(2,40)*two)*qx+(+r7(14,3)*two &
&                  -r6(9,7)*four+r5(5,9)*two+r5(5,11)*two+r5(14,21)*two-r4(9,34)*four &
&                  +r3(5,35)*two+r3(5,39)*two)*qz+(+r6(20,12)-r5(14,20)*two+r4(9,25)+r4(9,29) &
&                  +r4(9,41)-r3(5,51)*two+r2(4,50)+r2(4,55))*xx+rxyz(20)*xz+(+r6(9,10) &
&                  -r5(5,18)*two+r4(2,23)+r4(2,27)+r4(9,39)-r3(5,49)*two+r2(4,48)+r2(4,53)) &
&                  *zz+(+r5(14,24)*two-r4(9,37)*four+r3(5,38)*two+r3(5,42)*two)*xxz+( &
&                  +r5(9,23)*two-r4(5,36)*four+r3(2,37)*two+r3(2,41)*two)*xzz+rxyz(1)*xxzz
      eri(4,1,4,3)=r310+(+r7(13,3)+r7(13,4)*two-r6(8,7)*two-r6(8,8)*four+r5(4,9) &
&                  +r5(4,10)*two+r5(4,11)+r5(4,12)*two+r5(13,21)+r5(13,22)*two-r4(8,34)*two &
&                  -r4(8,35)*four+r3(4,35)+r3(4,36)*two+r3(4,39)+r3(4,40)*two)*qx+( &
&                  +r6(13,11)*two+r6(13,12)-r5(8,19)*four-r5(8,20)*two+r4(4,24)*two+r4(4,25) &
&                  +r4(4,28)*two+r4(4,29))*xx+rxyz(2)*xxx
      eri(5,1,4,3)=r301+(+r7(14,3)+r7(14,4)*two-r6(9,7)*two-r6(9,8)*four+r5(5,9) &
&                  +r5(5,10)*two+r5(5,11)+r5(5,12)*two+r5(14,21)+r5(14,22)*two-r4(9,34)*two &
&                  -r4(9,35)*four+r3(5,35)+r3(5,36)*two+r3(5,39)+r3(5,40)*two)*qx+(+r7(9,3) &
&                  -r6(5,7)*two+r5(2,9)+r5(2,11)+r5(9,21)*three-r4(5,34)*six+r3(2,35)*three &
&                  +r3(2,39)*three)*qz+(+r6(14,11)*two+r6(14,12)-r5(9,19)*four-r5(9,20)*two &
&                  +r4(5,24)*two+r4(5,25)+r4(5,28)*two+r4(5,29))*xx+(+r6(9,10)+r6(9,11)*two &
&                  -r5(5,18)*two-r5(5,19)*four+r4(2,23)+r4(2,24)*two+r4(2,27)+r4(2,28)*two &
&                  +r4(9,39)+r4(9,40)*two-r3(5,49)*two-r3(5,50)*four+r2(4,48)+r2(4,49)*two &
&                  +r2(4,53)+r2(4,54)*two)*xz+(+r5(14,24)-r4(9,37)*two+r3(5,38)+r3(5,42))*xxx &
&                  +(+r5(9,23)*two+r5(9,24)-r4(5,36)*four-r4(5,37)*two+r3(2,37)*two+r3(2,38) &
&                  +r3(2,41)*two+r3(2,42))*xxz+rxyz(1)*xxxz
      eri(6,1,4,3)=r211+(+r7(19,4)*two-r6(13,8)*four+r5(8,10)*two+r5(8,12)*two)*qx &
&                  +rxyz(17)*qz+(+r6(19,12)-r5(13,20)*two+r4(8,25)+r4(8,29))*xx+( &
&                  +r6(13,11)*two-r5(8,19)*four+r4(4,24)*two+r4(4,28)*two)*xz+rxyz(2)*xxz
      eri(1,2,4,3)=r220+(+r7(18,3)*two-r6(12,7)*four+r5(7,9)*two+r5(7,11)*two &
&                  +r5(9,21)*two-r4(5,34)*four+r3(2,35)*two+r3(2,39)*two)*qx+rxyz(4)*xx
      eri(2,2,4,3)=r040
      eri(3,2,4,3)=r022+(+r7(25,3)*two-r6(18,7)*four+r5(12,9)*two+r5(12,11)*two &
&                  +r5(14,21)*two-r4(9,34)*four+r3(5,35)*two+r3(5,39)*two)*qz+rxyz(4)*zz
      eri(4,2,4,3)=r130+rxyz(7)*qx
      eri(5,2,4,3)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,4,3)=r031+rxyz(7)*qz
      eri(1,3,4,3)=r202+(+r7(20,3)*two-r6(14,7)*four+r5(9,9)*two+r5(9,11)*two &
&                  +r5(9,21)*two-r4(5,34)*four+r3(2,35)*two+r3(2,39)*two)*qx+(+r7(14,4)*two &
&                  -r6(9,8)*four+r5(5,10)*two+r5(5,12)*two+r5(14,22)*two-r4(9,35)*four &
&                  +r3(5,36)*two+r3(5,40)*two)*qz+(+r6(20,10)-r5(14,18)*two+r4(9,23)+r4(9,27) &
&                  +r4(9,39)-r3(5,49)*two+r2(4,48)+r2(4,53))*xx+rxyz(20)*xz+(+r6(9,12) &
&                  -r5(5,20)*two+r4(2,25)+r4(2,29)+r4(9,41)-r3(5,51)*two+r2(4,50)+r2(4,55)) &
&                  *zz+(+r5(14,23)*two-r4(9,36)*four+r3(5,37)*two+r3(5,41)*two)*xxz+( &
&                  +r5(9,24)*two-r4(5,37)*four+r3(2,38)*two+r3(2,42)*two)*xzz+rxyz(1)*xxzz
      eri(2,3,4,3)=r022+(+r7(25,4)*two-r6(18,8)*four+r5(12,10)*two+r5(12,12)*two &
&                  +r5(14,22)*two-r4(9,35)*four+r3(5,36)*two+r3(5,40)*two)*qz+rxyz(5)*zz
      eri(3,3,4,3)=r004+(+r7(27,3)*two+r7(27,4)*two-r6(20,7)*four-r6(20,8)*four &
&                  +r5(14,9)*two+r5(14,10)*two+r5(14,11)*two+r5(14,12)*two+r5(14,21)*six &
&                  +r5(14,22)*six-r4(9,34)*p12-r4(9,35)*p12+r3(5,35)*six+r3(5,36)*six &
&                  +r3(5,39)*six+r3(5,40)*six)*qz+(+r6(20,10)+r6(20,11)*four+r6(20,12) &
&                  -r5(14,18)*two-r5(14,19)*eight-r5(14,20)*two+r4(9,23)+r4(9,24)*four &
&                  +r4(9,25)+r4(9,27)+r4(9,28)*four+r4(9,29)+r4(9,39)+r4(9,40)*four+r4(9,41) &
&                  -r3(5,49)*two-r3(5,50)*eight-r3(5,51)*two+r2(4,48)+r2(4,49)*four+r2(4,50) &
&                  +r2(4,53)+r2(4,54)*four+r2(4,55))*zz+(+r5(14,23)*two+r5(14,24)*two &
&                  -r4(9,36)*four-r4(9,37)*four+r3(5,37)*two+r3(5,38)*two+r3(5,41)*two &
&                  +r3(5,42)*two)*zzz+rxyz(1)*zzzz
      eri(4,3,4,3)=r112+rxyz(18)*qx+(+r7(19,4)*two-r6(13,8)*four+r5(8,10)*two &
&                  +r5(8,12)*two)*qz+(+r6(19,11)*two-r5(13,19)*four+r4(8,24)*two+r4(8,28)*two &
&                  )*xz+(+r6(13,12)-r5(8,20)*two+r4(4,25)+r4(4,29))*zz+rxyz(2)*xzz
      eri(5,3,4,3)=r103+(+r7(27,3)-r6(20,7)*two+r5(14,9)+r5(14,11)+r5(14,21)*three &
&                  -r4(9,34)*six+r3(5,35)*three+r3(5,39)*three)*qx+(+r7(20,3)+r7(20,4)*two &
&                  -r6(14,7)*two-r6(14,8)*four+r5(9,9)+r5(9,10)*two+r5(9,11)+r5(9,12)*two &
&                  +r5(9,21)+r5(9,22)*two-r4(5,34)*two-r4(5,35)*four+r3(2,35)+r3(2,36)*two &
&                  +r3(2,39)+r3(2,40)*two)*qz+(+r6(20,10)+r6(20,11)*two-r5(14,18)*two &
&                  -r5(14,19)*four+r4(9,23)+r4(9,24)*two+r4(9,27)+r4(9,28)*two+r4(9,39) &
&                  +r4(9,40)*two-r3(5,49)*two-r3(5,50)*four+r2(4,48)+r2(4,49)*two+r2(4,53) &
&                  +r2(4,54)*two)*xz+(+r6(14,11)*two+r6(14,12)-r5(9,19)*four-r5(9,20)*two &
&                  +r4(5,24)*two+r4(5,25)+r4(5,28)*two+r4(5,29))*zz+(+r5(14,23)*two+r5(14,24) &
&                  -r4(9,36)*four-r4(9,37)*two+r3(5,37)*two+r3(5,38)+r3(5,41)*two+r3(5,42)) &
&                  *xzz+(+r5(9,24)-r4(5,37)*two+r3(2,38)+r3(2,42))*zzz+rxyz(1)*xzzz
      eri(6,3,4,3)=r013+(+r7(26,3)+r7(26,4)*two-r6(19,7)*two-r6(19,8)*four+r5(13,9) &
&                  +r5(13,10)*two+r5(13,11)+r5(13,12)*two+r5(13,21)+r5(13,22)*two &
&                  -r4(8,34)*two-r4(8,35)*four+r3(4,35)+r3(4,36)*two+r3(4,39)+r3(4,40)*two) &
&                  *qz+(+r6(19,11)*two+r6(19,12)-r5(13,19)*four-r5(13,20)*two+r4(8,24)*two &
&                  +r4(8,25)+r4(8,28)*two+r4(8,29))*zz+rxyz(2)*zzz
      eri(1,4,4,3)=r310+(+r7(13,3)*two+r7(13,4)-r6(8,7)*four-r6(8,8)*two+r5(4,9)*two &
&                  +r5(4,10)+r5(4,11)*two+r5(4,12)+r5(13,21)*two+r5(13,22)-r4(8,34)*four &
&                  -r4(8,35)*two+r3(4,35)*two+r3(4,36)+r3(4,39)*two+r3(4,40))*qx+(+r6(13,10) &
&                  +r6(13,11)*two-r5(8,18)*two-r5(8,19)*four+r4(4,23)+r4(4,24)*two+r4(4,27) &
&                  +r4(4,28)*two)*xx+rxyz(3)*xxx
      eri(2,4,4,3)=r130+rxyz(8)*qx
      eri(3,4,4,3)=r112+rxyz(19)*qx+(+r7(19,3)*two-r6(13,7)*four+r5(8,9)*two+r5(8,11)*two &
&                  )*qz+(+r6(19,11)*two-r5(13,19)*four+r4(8,24)*two+r4(8,28)*two)*xz+( &
&                  +r6(13,10)-r5(8,18)*two+r4(4,23)+r4(4,27))*zz+rxyz(3)*xzz
      eri(4,4,4,3)=r220+(+r7(18,3)+r7(18,4)-r6(12,7)*two-r6(12,8)*two+r5(7,9)+r5(7,10) &
&                  +r5(7,11)+r5(7,12)+r5(9,21)+r5(9,22)-r4(5,34)*two-r4(5,35)*two+r3(2,35) &
&                  +r3(2,36)+r3(2,39)+r3(2,40))*qx+rxyz(6)*xx
      eri(5,4,4,3)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(13,10)+r6(13,11) &
&                  -r5(8,18)*two-r5(8,19)*two+r4(4,23)+r4(4,24)+r4(4,27)+r4(4,28))*xz+rxyz(3) &
&                  *xxz
      eri(6,4,4,3)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,4,3)=r301+(+r7(14,3)*two+r7(14,4)-r6(9,7)*four-r6(9,8)*two+r5(5,9)*two &
&                  +r5(5,10)+r5(5,11)*two+r5(5,12)+r5(14,21)*two+r5(14,22)-r4(9,34)*four &
&                  -r4(9,35)*two+r3(5,35)*two+r3(5,36)+r3(5,39)*two+r3(5,40))*qx+(+r7(9,4) &
&                  -r6(5,8)*two+r5(2,10)+r5(2,12)+r5(9,22)*three-r4(5,35)*six+r3(2,36)*three &
&                  +r3(2,40)*three)*qz+(+r6(14,10)+r6(14,11)*two-r5(9,18)*two-r5(9,19)*four &
&                  +r4(5,23)+r4(5,24)*two+r4(5,27)+r4(5,28)*two)*xx+(+r6(9,11)*two+r6(9,12) &
&                  -r5(5,19)*four-r5(5,20)*two+r4(2,24)*two+r4(2,25)+r4(2,28)*two+r4(2,29) &
&                  +r4(9,40)*two+r4(9,41)-r3(5,50)*four-r3(5,51)*two+r2(4,49)*two+r2(4,50) &
&                  +r2(4,54)*two+r2(4,55))*xz+(+r5(14,23)-r4(9,36)*two+r3(5,37)+r3(5,41))*xxx &
&                  +(+r5(9,23)+r5(9,24)*two-r4(5,36)*two-r4(5,37)*four+r3(2,37)+r3(2,38)*two &
&                  +r3(2,41)+r3(2,42)*two)*xxz+rxyz(1)*xxxz
      eri(2,5,4,3)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,4,3)=r103+(+r7(27,4)-r6(20,8)*two+r5(14,10)+r5(14,12)+r5(14,22)*three &
&                  -r4(9,35)*six+r3(5,36)*three+r3(5,40)*three)*qx+(+r7(20,3)*two+r7(20,4) &
&                  -r6(14,7)*four-r6(14,8)*two+r5(9,9)*two+r5(9,10)+r5(9,11)*two+r5(9,12) &
&                  +r5(9,21)*two+r5(9,22)-r4(5,34)*four-r4(5,35)*two+r3(2,35)*two+r3(2,36) &
&                  +r3(2,39)*two+r3(2,40))*qz+(+r6(20,11)*two+r6(20,12)-r5(14,19)*four &
&                  -r5(14,20)*two+r4(9,24)*two+r4(9,25)+r4(9,28)*two+r4(9,29)+r4(9,40)*two &
&                  +r4(9,41)-r3(5,50)*four-r3(5,51)*two+r2(4,49)*two+r2(4,50)+r2(4,54)*two &
&                  +r2(4,55))*xz+(+r6(14,10)+r6(14,11)*two-r5(9,18)*two-r5(9,19)*four &
&                  +r4(5,23)+r4(5,24)*two+r4(5,27)+r4(5,28)*two)*zz+(+r5(14,23)+r5(14,24)*two &
&                  -r4(9,36)*two-r4(9,37)*four+r3(5,37)+r3(5,38)*two+r3(5,41)+r3(5,42)*two) &
&                  *xzz+(+r5(9,23)-r4(5,36)*two+r3(2,37)+r3(2,41))*zzz+rxyz(1)*xzzz
      eri(4,5,4,3)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(13,11)+r6(13,12) &
&                  -r5(8,19)*two-r5(8,20)*two+r4(4,24)+r4(4,25)+r4(4,28)+r4(4,29))*xz+rxyz(2) &
&                  *xxz
      eri(5,5,4,3)=r202+(+r7(20,3)+r7(20,4)-r6(14,7)*two-r6(14,8)*two+r5(9,9)+r5(9,10) &
&                  +r5(9,11)+r5(9,12)+r5(9,21)+r5(9,22)-r4(5,34)*two-r4(5,35)*two+r3(2,35) &
&                  +r3(2,36)+r3(2,39)+r3(2,40))*qx+(+r7(14,3)+r7(14,4)-r6(9,7)*two &
&                  -r6(9,8)*two+r5(5,9)+r5(5,10)+r5(5,11)+r5(5,12)+r5(14,21)+r5(14,22) &
&                  -r4(9,34)*two-r4(9,35)*two+r3(5,35)+r3(5,36)+r3(5,39)+r3(5,40))*qz+( &
&                  +r6(20,11)-r5(14,19)*two+r4(9,24)+r4(9,28)+r4(9,40)-r3(5,50)*two+r2(4,49) &
&                  +r2(4,54))*xx+(+r6(14,10)+r6(14,11)*two+r6(14,12)-r5(9,18)*two &
&                  -r5(9,19)*four-r5(9,20)*two+r4(5,23)+r4(5,24)*two+r4(5,25)+r4(5,27) &
&                  +r4(5,28)*two+r4(5,29))*xz+(+r6(9,11)-r5(5,19)*two+r4(2,24)+r4(2,28) &
&                  +r4(9,40)-r3(5,50)*two+r2(4,49)+r2(4,54))*zz+(+r5(14,23)+r5(14,24) &
&                  -r4(9,36)*two-r4(9,37)*two+r3(5,37)+r3(5,38)+r3(5,41)+r3(5,42))*xxz+( &
&                  +r5(9,23)+r5(9,24)-r4(5,36)*two-r4(5,37)*two+r3(2,37)+r3(2,38)+r3(2,41) &
&                  +r3(2,42))*xzz+rxyz(1)*xxzz
      eri(6,5,4,3)=r112+rxyz(19)*qx+(+r7(19,3)+r7(19,4)-r6(13,7)*two-r6(13,8)*two+r5(8,9) &
&                  +r5(8,10)+r5(8,11)+r5(8,12))*qz+(+r6(19,11)+r6(19,12)-r5(13,19)*two &
&                  -r5(13,20)*two+r4(8,24)+r4(8,25)+r4(8,28)+r4(8,29))*xz+rxyz(11)*zz+rxyz(2) &
&                  *xzz
      eri(1,6,4,3)=r211+(+r7(19,3)*two-r6(13,7)*four+r5(8,9)*two+r5(8,11)*two)*qx &
&                  +rxyz(16)*qz+(+r6(19,10)-r5(13,18)*two+r4(8,23)+r4(8,27))*xx+( &
&                  +r6(13,11)*two-r5(8,19)*four+r4(4,24)*two+r4(4,28)*two)*xz+rxyz(3)*xxz
      eri(2,6,4,3)=r031+rxyz(8)*qz
      eri(3,6,4,3)=r013+(+r7(26,3)*two+r7(26,4)-r6(19,7)*four-r6(19,8)*two+r5(13,9)*two &
&                  +r5(13,10)+r5(13,11)*two+r5(13,12)+r5(13,21)*two+r5(13,22)-r4(8,34)*four &
&                  -r4(8,35)*two+r3(4,35)*two+r3(4,36)+r3(4,39)*two+r3(4,40))*qz+(+r6(19,10) &
&                  +r6(19,11)*two-r5(13,18)*two-r5(13,19)*four+r4(8,23)+r4(8,24)*two+r4(8,27) &
&                  +r4(8,28)*two)*zz+rxyz(3)*zzz
      eri(4,6,4,3)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,4,3)=r112+rxyz(18)*qx+(+r7(19,3)+r7(19,4)-r6(13,7)*two-r6(13,8)*two+r5(8,9) &
&                  +r5(8,10)+r5(8,11)+r5(8,12))*qz+(+r6(19,10)+r6(19,11)-r5(13,18)*two &
&                  -r5(13,19)*two+r4(8,23)+r4(8,24)+r4(8,27)+r4(8,28))*xz+rxyz(11)*zz+rxyz(3) &
&                  *xzz
      eri(6,6,4,3)=r022+(+r7(25,3)+r7(25,4)-r6(18,7)*two-r6(18,8)*two+r5(12,9)+r5(12,10) &
&                  +r5(12,11)+r5(12,12)+r5(14,21)+r5(14,22)-r4(9,34)*two-r4(9,35)*two &
&                  +r3(5,35)+r3(5,36)+r3(5,39)+r3(5,40))*qz+rxyz(6)*zz
!
      r400= r8(10)-r7(6,1)-r7(6,2)*two+r6(3,2)*two+r6(3,3)+r6(3,4)*three+r6(10,9)*six &
&          -r5(1,2)-r5(1,3)-r5(1,4)*two-r5(6,13)*six-r5(6,17)*p12+r4(3,18)*p12 &
&          +r4(3,22)*six+r4(3,26)*p18+r4(10,38)*three-r3(1,15)*six-r3(1,19)*six &
&          -r3(1,23)*p12-r3(6,43)*three-r3(6,48)*six+r2(5,42)*six+r2(5,47)*three &
&          +r2(5,52)*nine-r1(1,26)*three-r1(1,31)*three-r1(1,36)*six
      r310= r8(14)-r7(9,1)-r7(9,2)*two+r6(5,2)*two+r6(5,3)+r6(5,4)*three+r6(14,9)*three &
&          -r5(2,2)-r5(2,3)-r5(2,4)*two-r5(9,13)*three-r5(9,17)*six+r4(5,18)*six &
&          +r4(5,22)*three+r4(5,26)*nine-r3(2,15)*three-r3(2,19)*three-r3(2,23)*six
      r301= r8(15)-r7(10,1)-r7(10,2)*two+r6(6,2)*two+r6(6,3)+r6(6,4)*three+r6(15,9)*three &
&          -r5(3,2)-r5(3,3)-r5(3,4)*two-r5(10,13)*three-r5(10,17)*six+r4(6,18)*six &
&          +r4(6,22)*three+r4(6,26)*nine-r3(3,15)*three-r3(3,19)*three-r3(3,23)*six
      r220= r8(19)-r7(13,1)-r7(13,2)*two+r6(8,2)*two+r6(8,3)+r6(8,4)*three+r6(10,9) &
&          +r6(19,9)-r5(4,2)-r5(4,3)-r5(4,4)*two-r5(6,13)-r5(13,13)-r5(6,17)*two &
&          -r5(13,17)*two+r4(3,18)*two+r4(8,18)*two+r4(3,22)+r4(8,22)+r4(3,26)*three &
&          +r4(8,26)*three+r4(10,38)-r3(1,15)-r3(4,15)-r3(1,19)-r3(4,19)-r3(1,23)*two &
&          -r3(4,23)*two-r3(6,43)-r3(6,48)*two+r2(5,42)*two+r2(5,47)+r2(5,52)*three &
&          -r1(1,26)-r1(1,31)-r1(1,36)*two
      r211= r8(20)-r7(14,1)-r7(14,2)*two+r6(9,2)*two+r6(9,3)+r6(9,4)*three+r6(20,9)-r5(5,2) &
&          -r5(5,3)-r5(5,4)*two-r5(14,13)-r5(14,17)*two+r4(9,18)*two+r4(9,22) &
&          +r4(9,26)*three-r3(5,15)-r3(5,19)-r3(5,23)*two
      r202= r8(21)-r7(15,1)-r7(15,2)*two+r6(10,2)*two+r6(10,3)+r6(10,4)*three+r6(10,9) &
&          +r6(21,9)-r5(6,2)-r5(6,3)-r5(6,4)*two-r5(6,13)-r5(15,13)-r5(6,17)*two &
&          -r5(15,17)*two+r4(3,18)*two+r4(10,18)*two+r4(3,22)+r4(10,22)+r4(3,26)*three &
&          +r4(10,26)*three+r4(10,38)-r3(1,15)-r3(6,15)-r3(1,19)-r3(6,19)-r3(1,23)*two &
&          -r3(6,23)*two-r3(6,43)-r3(6,48)*two+r2(5,42)*two+r2(5,47)+r2(5,52)*three &
&          -r1(1,26)-r1(1,31)-r1(1,36)*two
      r130= r8(25)-r7(18,1)-r7(18,2)*two+r6(12,2)*two+r6(12,3)+r6(12,4)*three &
&          +r6(14,9)*three-r5(7,2)-r5(7,3)-r5(7,4)*two-r5(9,13)*three-r5(9,17)*six &
&          +r4(5,18)*six+r4(5,22)*three+r4(5,26)*nine-r3(2,15)*three-r3(2,19)*three &
&          -r3(2,23)*six
      r121= r8(26)-r7(19,1)-r7(19,2)*two+r6(13,2)*two+r6(13,3)+r6(13,4)*three+r6(15,9) &
&          -r5(8,2)-r5(8,3)-r5(8,4)*two-r5(10,13)-r5(10,17)*two+r4(6,18)*two+r4(6,22) &
&          +r4(6,26)*three-r3(3,15)-r3(3,19)-r3(3,23)*two
      r112= r8(27)-r7(20,1)-r7(20,2)*two+r6(14,2)*two+r6(14,3)+r6(14,4)*three+r6(14,9) &
&          -r5(9,2)-r5(9,3)-r5(9,4)*two-r5(9,13)-r5(9,17)*two+r4(5,18)*two+r4(5,22) &
&          +r4(5,26)*three-r3(2,15)-r3(2,19)-r3(2,23)*two
      r103= r8(28)-r7(21,1)-r7(21,2)*two+r6(15,2)*two+r6(15,3)+r6(15,4)*three &
&          +r6(15,9)*three-r5(10,2)-r5(10,3)-r5(10,4)*two-r5(10,13)*three-r5(10,17)*six &
&          +r4(6,18)*six+r4(6,22)*three+r4(6,26)*nine-r3(3,15)*three-r3(3,19)*three &
&          -r3(3,23)*six
      r040= r8(32)-r7(24,1)-r7(24,2)*two+r6(17,2)*two+r6(17,3)+r6(17,4)*three+r6(19,9)*six &
&          -r5(11,2)-r5(11,3)-r5(11,4)*two-r5(13,13)*six-r5(13,17)*p12+r4(8,18)*p12 &
&          +r4(8,22)*six+r4(8,26)*p18+r4(10,38)*three-r3(4,15)*six-r3(4,19)*six &
&          -r3(4,23)*p12-r3(6,43)*three-r3(6,48)*six+r2(5,42)*six+r2(5,47)*three &
&          +r2(5,52)*nine-r1(1,26)*three-r1(1,31)*three-r1(1,36)*six
      r031= r8(33)-r7(25,1)-r7(25,2)*two+r6(18,2)*two+r6(18,3)+r6(18,4)*three &
&          +r6(20,9)*three-r5(12,2)-r5(12,3)-r5(12,4)*two-r5(14,13)*three-r5(14,17)*six &
&          +r4(9,18)*six+r4(9,22)*three+r4(9,26)*nine-r3(5,15)*three-r3(5,19)*three &
&          -r3(5,23)*six
      r022= r8(34)-r7(26,1)-r7(26,2)*two+r6(19,2)*two+r6(19,3)+r6(19,4)*three+r6(19,9) &
&          +r6(21,9)-r5(13,2)-r5(13,3)-r5(13,4)*two-r5(13,13)-r5(15,13)-r5(13,17)*two &
&          -r5(15,17)*two+r4(8,18)*two+r4(10,18)*two+r4(8,22)+r4(10,22)+r4(8,26)*three &
&          +r4(10,26)*three+r4(10,38)-r3(4,15)-r3(6,15)-r3(4,19)-r3(6,19)-r3(4,23)*two &
&          -r3(6,23)*two-r3(6,43)-r3(6,48)*two+r2(5,42)*two+r2(5,47)+r2(5,52)*three &
&          -r1(1,26)-r1(1,31)-r1(1,36)*two
      r013= r8(35)-r7(27,1)-r7(27,2)*two+r6(20,2)*two+r6(20,3)+r6(20,4)*three &
&          +r6(20,9)*three-r5(14,2)-r5(14,3)-r5(14,4)*two-r5(14,13)*three-r5(14,17)*six &
&          +r4(9,18)*six+r4(9,22)*three+r4(9,26)*nine-r3(5,15)*three-r3(5,19)*three &
&          -r3(5,23)*six
      r004= r8(36)-r7(28,1)-r7(28,2)*two+r6(21,2)*two+r6(21,3)+r6(21,4)*three+r6(21,9)*six &
&          -r5(15,2)-r5(15,3)-r5(15,4)*two-r5(15,13)*six-r5(15,17)*p12+r4(10,18)*p12 &
&          +r4(10,22)*six+r4(10,26)*p18+r4(10,38)*three-r3(6,15)*six-r3(6,19)*six &
&          -r3(6,23)*p12-r3(6,43)*three-r3(6,48)*six+r2(5,42)*six+r2(5,47)*three &
&          +r2(5,52)*nine-r1(1,26)*three-r1(1,31)*three-r1(1,36)*six
      rxyz(1)=+r4(10,42)-r3(6,47)-r3(6,52)*two+r2(5,46)*two+r2(5,51)+r2(5,56)*three &
&             -r1(1,30)-r1(1,35)-r1(1,40)*two
      rxyz(2)=+r5(14,24)-r4(9,33)-r4(9,37)*two+r3(5,34)*two+r3(5,38)+r3(5,42)*three &
&             -r2(4,28)-r2(4,32)-r2(4,36)*two
      rxyz(3)=+r5(14,23)-r4(9,32)-r4(9,36)*two+r3(5,33)*two+r3(5,37)+r3(5,41)*three &
&             -r2(4,27)-r2(4,31)-r2(4,35)*two
      rxyz(4)=+r6(19,10)-r5(13,14)-r5(13,18)*two+r4(8,19)*two+r4(8,23)+r4(8,27)*three &
&             +r4(10,39)-r3(4,16)-r3(4,20)-r3(4,24)*two-r3(6,44)-r3(6,49)*two+r2(5,43)*two &
&             +r2(5,48)+r2(5,53)*three-r1(1,27)-r1(1,32)-r1(1,37)*two
      rxyz(5)=+r6(19,12)-r5(13,16)-r5(13,20)*two+r4(8,21)*two+r4(8,25)+r4(8,29)*three &
&             +r4(10,41)-r3(4,18)-r3(4,22)-r3(4,26)*two-r3(6,46)-r3(6,51)*two+r2(5,45)*two &
&             +r2(5,50)+r2(5,55)*three-r1(1,29)-r1(1,34)-r1(1,39)*two
      rxyz(6)=+r6(19,11)-r5(13,15)-r5(13,19)*two+r4(8,20)*two+r4(8,24)+r4(8,28)*three &
&             +r4(10,40)-r3(4,17)-r3(4,21)-r3(4,25)*two-r3(6,45)-r3(6,50)*two+r2(5,44)*two &
&             +r2(5,49)+r2(5,54)*three-r1(1,28)-r1(1,33)-r1(1,38)*two
      rxyz(7)=+r7(25,3)-r6(18,5)-r6(18,7)*two+r5(12,7)*two+r5(12,9)+r5(12,11)*three &
&             +r5(14,21)*three-r4(7,8)-r4(7,10)-r4(7,12)*two-r4(9,30)*three-r4(9,34)*six &
&             +r3(5,31)*six+r3(5,35)*three+r3(5,39)*nine-r2(4,25)*three-r2(4,29)*three &
&             -r2(4,33)*six
      rxyz(8)=+r7(25,4)-r6(18,6)-r6(18,8)*two+r5(12,8)*two+r5(12,10)+r5(12,12)*three &
&             +r5(14,22)*three-r4(7,9)-r4(7,11)-r4(7,13)*two-r4(9,31)*three-r4(9,35)*six &
&             +r3(5,32)*six+r3(5,36)*three+r3(5,40)*nine-r2(4,26)*three-r2(4,30)*three &
&             -r2(4,34)*six
      rxyz(9)=+r7(20,3)+r7(20,4)-r6(14,5)-r6(14,6)-r6(14,7)*two-r6(14,8)*two+r5(9,7)*two &
&             +r5(9,8)*two+r5(9,9)+r5(9,10)+r5(9,11)*three+r5(9,12)*three-r4(5,8)-r4(5,9) &
&             -r4(5,10)-r4(5,11)-r4(5,12)*two-r4(5,13)*two
      rxyz(10)=+r6(20,11)-r5(14,15)-r5(14,19)*two+r4(9,20)*two+r4(9,24)+r4(9,28)*three &
&             -r3(5,17)-r3(5,21)-r3(5,25)*two
      rxyz(11)=+r6(14,11)-r5(9,15)-r5(9,19)*two+r4(5,20)*two+r4(5,24)+r4(5,28)*three &
&             -r3(2,17)-r3(2,21)-r3(2,25)*two
      rxyz(12)=+r7(26,3)-r6(19,5)-r6(19,7)*two+r5(13,7)*two+r5(13,9)+r5(13,11)*three &
&             +r5(15,21)-r4(8,8)-r4(8,10)-r4(8,12)*two-r4(10,30)-r4(10,34)*two &
&             +r3(6,31)*two+r3(6,35)+r3(6,39)*three-r2(5,25)-r2(5,29)-r2(5,33)*two
      rxyz(13)=+r7(26,4)-r6(19,6)-r6(19,8)*two+r5(13,8)*two+r5(13,10)+r5(13,12)*three &
&             +r5(15,22)-r4(8,9)-r4(8,11)-r4(8,13)*two-r4(10,31)-r4(10,35)*two &
&             +r3(6,32)*two+r3(6,36)+r3(6,40)*three-r2(5,26)-r2(5,30)-r2(5,34)*two
      rxyz(14)=+r7(19,4)-r6(13,6)-r6(13,8)*two+r5(8,8)*two+r5(8,10)+r5(8,12)*three &
&             +r5(10,22)-r4(4,9)-r4(4,11)-r4(4,13)*two-r4(6,31)-r4(6,35)*two+r3(3,32)*two &
&             +r3(3,36)+r3(3,40)*three-r2(1,26)-r2(1,30)-r2(1,34)*two
      rxyz(15)=+r7(19,3)-r6(13,5)-r6(13,7)*two+r5(8,7)*two+r5(8,9)+r5(8,11)*three &
&             +r5(10,21)-r4(4,8)-r4(4,10)-r4(4,12)*two-r4(6,30)-r4(6,34)*two+r3(3,31)*two &
&             +r3(3,35)+r3(3,39)*three-r2(1,25)-r2(1,29)-r2(1,33)*two
      rxyz(16)=+r7(14,4)-r6(9,6)-r6(9,8)*two+r5(5,8)*two+r5(5,10)+r5(5,12)*three+r5(14,22) &
&             -r4(2,9)-r4(2,11)-r4(2,13)*two-r4(9,31)-r4(9,35)*two+r3(5,32)*two+r3(5,36) &
&             +r3(5,40)*three-r2(4,26)-r2(4,30)-r2(4,34)*two
      rxyz(17)=+r7(14,3)-r6(9,5)-r6(9,7)*two+r5(5,7)*two+r5(5,9)+r5(5,11)*three+r5(14,21) &
&             -r4(2,8)-r4(2,10)-r4(2,12)*two-r4(9,30)-r4(9,34)*two+r3(5,31)*two+r3(5,35) &
&             +r3(5,39)*three-r2(4,25)-r2(4,29)-r2(4,33)*two
      rxyz(18)=+r7(27,3)-r6(20,5)-r6(20,7)*two+r5(14,7)*two+r5(14,9)+r5(14,11)*three &
&             +r5(14,21)-r4(9,8)-r4(9,10)-r4(9,12)*two-r4(9,30)-r4(9,34)*two+r3(5,31)*two &
&             +r3(5,35)+r3(5,39)*three-r2(4,25)-r2(4,29)-r2(4,33)*two
      rxyz(19)=+r7(27,4)-r6(20,6)-r6(20,8)*two+r5(14,8)*two+r5(14,10)+r5(14,12)*three &
&             +r5(14,22)-r4(9,9)-r4(9,11)-r4(9,13)*two-r4(9,31)-r4(9,35)*two+r3(5,32)*two &
&             +r3(5,36)+r3(5,40)*three-r2(4,26)-r2(4,30)-r2(4,34)*two
      rxyz(20)=+r6(15,11)*four-r5(10,15)*four-r5(10,19)*eight+r4(6,20)*eight+r4(6,24)*four &
&             +r4(6,28)*p12-r3(3,17)*four-r3(3,21)*four-r3(3,25)*eight
      eri(1,1,5,3)=r400+(+r7(10,3)*two+r7(10,4)*two-r6(6,5)*two-r6(6,6)*two-r6(6,7)*four &
&                  -r6(6,8)*four+r5(3,7)*four+r5(3,8)*four+r5(3,9)*two+r5(3,10)*two &
&                  +r5(3,11)*six+r5(3,12)*six+r5(10,21)*six+r5(10,22)*six-r4(1,8)*two &
&                  -r4(1,9)*two-r4(1,10)*two-r4(1,11)*two-r4(1,12)*four-r4(1,13)*four &
&                  -r4(6,30)*six-r4(6,31)*six-r4(6,34)*p12-r4(6,35)*p12+r3(3,31)*p12 &
&                  +r3(3,32)*p12+r3(3,35)*six+r3(3,36)*six+r3(3,39)*p18+r3(3,40)*p18 &
&                  -r2(1,25)*six-r2(1,26)*six-r2(1,29)*six-r2(1,30)*six-r2(1,33)*p12 &
&                  -r2(1,34)*p12)*qx+(+r6(10,10)+r6(10,11)*four+r6(10,12)-r5(6,14) &
&                  -r5(6,15)*four-r5(6,16)-r5(6,18)*two-r5(6,19)*eight-r5(6,20)*two &
&                  +r4(3,19)*two+r4(3,20)*eight+r4(3,21)*two+r4(3,23)+r4(3,24)*four+r4(3,25) &
&                  +r4(3,27)*three+r4(3,28)*p12+r4(3,29)*three+r4(10,39)+r4(10,40)*four &
&                  +r4(10,41)-r3(1,16)-r3(1,17)*four-r3(1,18)-r3(1,20)-r3(1,21)*four-r3(1,22) &
&                  -r3(1,24)*two-r3(1,25)*eight-r3(1,26)*two-r3(6,44)-r3(6,45)*four-r3(6,46) &
&                  -r3(6,49)*two-r3(6,50)*eight-r3(6,51)*two+r2(5,43)*two+r2(5,44)*eight &
&                  +r2(5,45)*two+r2(5,48)+r2(5,49)*four+r2(5,50)+r2(5,53)*three+r2(5,54)*p12 &
&                  +r2(5,55)*three-r1(1,27)-r1(1,28)*four-r1(1,29)-r1(1,32)-r1(1,33)*four &
&                  -r1(1,34)-r1(1,37)*two-r1(1,38)*eight-r1(1,39)*two)*xx+(+r5(10,23)*two &
&                  +r5(10,24)*two-r4(6,32)*two-r4(6,33)*two-r4(6,36)*four-r4(6,37)*four &
&                  +r3(3,33)*four+r3(3,34)*four+r3(3,37)*two+r3(3,38)*two+r3(3,41)*six &
&                  +r3(3,42)*six-r2(1,27)*two-r2(1,28)*two-r2(1,31)*two-r2(1,32)*two &
&                  -r2(1,35)*four-r2(1,36)*four)*xxx+rxyz(1)*xxxx
      eri(2,1,5,3)=r220+(+r7(19,4)*two-r6(13,6)*two-r6(13,8)*four+r5(8,8)*four &
&                  +r5(8,10)*two+r5(8,12)*six+r5(10,22)*two-r4(4,9)*two-r4(4,11)*two &
&                  -r4(4,13)*four-r4(6,31)*two-r4(6,35)*four+r3(3,32)*four+r3(3,36)*two &
&                  +r3(3,40)*six-r2(1,26)*two-r2(1,30)*two-r2(1,34)*four)*qx+rxyz(5)*xx
      eri(3,1,5,3)=r202+(+r7(21,4)*two-r6(15,6)*two-r6(15,8)*four+r5(10,8)*four &
&                  +r5(10,10)*two+r5(10,12)*six+r5(10,22)*two-r4(6,9)*two-r4(6,11)*two &
&                  -r4(6,13)*four-r4(6,31)*two-r4(6,35)*four+r3(3,32)*four+r3(3,36)*two &
&                  +r3(3,40)*six-r2(1,26)*two-r2(1,30)*two-r2(1,34)*four)*qx+(+r7(15,3)*two &
&                  -r6(10,5)*two-r6(10,7)*four+r5(6,7)*four+r5(6,9)*two+r5(6,11)*six &
&                  +r5(15,21)*two-r4(3,8)*two-r4(3,10)*two-r4(3,12)*four-r4(10,30)*two &
&                  -r4(10,34)*four+r3(6,31)*four+r3(6,35)*two+r3(6,39)*six-r2(5,25)*two &
&                  -r2(5,29)*two-r2(5,33)*four)*qz+(+r6(21,12)-r5(15,16)-r5(15,20)*two &
&                  +r4(10,21)*two+r4(10,25)+r4(10,29)*three+r4(10,41)-r3(6,18)-r3(6,22) &
&                  -r3(6,26)*two-r3(6,46)-r3(6,51)*two+r2(5,45)*two+r2(5,50)+r2(5,55)*three &
&                  -r1(1,29)-r1(1,34)-r1(1,39)*two)*xx+rxyz(20)*xz+(+r6(10,10)-r5(6,14) &
&                  -r5(6,18)*two+r4(3,19)*two+r4(3,23)+r4(3,27)*three+r4(10,39)-r3(1,16) &
&                  -r3(1,20)-r3(1,24)*two-r3(6,44)-r3(6,49)*two+r2(5,43)*two+r2(5,48) &
&                  +r2(5,53)*three-r1(1,27)-r1(1,32)-r1(1,37)*two)*zz+(+r5(15,24)*two &
&                  -r4(10,33)*two-r4(10,37)*four+r3(6,34)*four+r3(6,38)*two+r3(6,42)*six &
&                  -r2(5,28)*two-r2(5,32)*two-r2(5,36)*four)*xxz+(+r5(10,23)*two-r4(6,32)*two &
&                  -r4(6,36)*four+r3(3,33)*four+r3(3,37)*two+r3(3,41)*six-r2(1,27)*two &
&                  -r2(1,31)*two-r2(1,35)*four)*xzz+rxyz(1)*xxzz
      eri(4,1,5,3)=r310+(+r7(14,3)+r7(14,4)*two-r6(9,5)-r6(9,6)*two-r6(9,7)*two &
&                  -r6(9,8)*four+r5(5,7)*two+r5(5,8)*four+r5(5,9)+r5(5,10)*two+r5(5,11)*three &
&                  +r5(5,12)*six+r5(14,21)+r5(14,22)*two-r4(2,8)-r4(2,9)*two-r4(2,10) &
&                  -r4(2,11)*two-r4(2,12)*two-r4(2,13)*four-r4(9,30)-r4(9,31)*two &
&                  -r4(9,34)*two-r4(9,35)*four+r3(5,31)*two+r3(5,32)*four+r3(5,35) &
&                  +r3(5,36)*two+r3(5,39)*three+r3(5,40)*six-r2(4,25)-r2(4,26)*two-r2(4,29) &
&                  -r2(4,30)*two-r2(4,33)*two-r2(4,34)*four)*qx+(+r6(14,11)*two+r6(14,12) &
&                  -r5(9,15)*two-r5(9,16)-r5(9,19)*four-r5(9,20)*two+r4(5,20)*four &
&                  +r4(5,21)*two+r4(5,24)*two+r4(5,25)+r4(5,28)*six+r4(5,29)*three &
&                  -r3(2,17)*two-r3(2,18)-r3(2,21)*two-r3(2,22)-r3(2,25)*four-r3(2,26)*two) &
&                  *xx+rxyz(2)*xxx
      eri(5,1,5,3)=r301+(+r7(15,3)+r7(15,4)*two-r6(10,5)-r6(10,6)*two-r6(10,7)*two &
&                  -r6(10,8)*four+r5(6,7)*two+r5(6,8)*four+r5(6,9)+r5(6,10)*two &
&                  +r5(6,11)*three+r5(6,12)*six+r5(15,21)+r5(15,22)*two-r4(3,8)-r4(3,9)*two &
&                  -r4(3,10)-r4(3,11)*two-r4(3,12)*two-r4(3,13)*four-r4(10,30)-r4(10,31)*two &
&                  -r4(10,34)*two-r4(10,35)*four+r3(6,31)*two+r3(6,32)*four+r3(6,35) &
&                  +r3(6,36)*two+r3(6,39)*three+r3(6,40)*six-r2(5,25)-r2(5,26)*two-r2(5,29) &
&                  -r2(5,30)*two-r2(5,33)*two-r2(5,34)*four)*qx+(+r7(10,3)-r6(6,5) &
&                  -r6(6,7)*two+r5(3,7)*two+r5(3,9)+r5(3,11)*three+r5(10,21)*three-r4(1,8) &
&                  -r4(1,10)-r4(1,12)*two-r4(6,30)*three-r4(6,34)*six+r3(3,31)*six &
&                  +r3(3,35)*three+r3(3,39)*nine-r2(1,25)*three-r2(1,29)*three-r2(1,33)*six) &
&                  *qz+(+r6(15,11)*two+r6(15,12)-r5(10,15)*two-r5(10,16)-r5(10,19)*four &
&                  -r5(10,20)*two+r4(6,20)*four+r4(6,21)*two+r4(6,24)*two+r4(6,25) &
&                  +r4(6,28)*six+r4(6,29)*three-r3(3,17)*two-r3(3,18)-r3(3,21)*two-r3(3,22) &
&                  -r3(3,25)*four-r3(3,26)*two)*xx+(+r6(10,10)+r6(10,11)*two-r5(6,14) &
&                  -r5(6,15)*two-r5(6,18)*two-r5(6,19)*four+r4(3,19)*two+r4(3,20)*four &
&                  +r4(3,23)+r4(3,24)*two+r4(3,27)*three+r4(3,28)*six+r4(10,39)+r4(10,40)*two &
&                  -r3(1,16)-r3(1,17)*two-r3(1,20)-r3(1,21)*two-r3(1,24)*two-r3(1,25)*four &
&                  -r3(6,44)-r3(6,45)*two-r3(6,49)*two-r3(6,50)*four+r2(5,43)*two &
&                  +r2(5,44)*four+r2(5,48)+r2(5,49)*two+r2(5,53)*three+r2(5,54)*six-r1(1,27) &
&                  -r1(1,28)*two-r1(1,32)-r1(1,33)*two-r1(1,37)*two-r1(1,38)*four)*xz+( &
&                  +r5(15,24)-r4(10,33)-r4(10,37)*two+r3(6,34)*two+r3(6,38)+r3(6,42)*three &
&                  -r2(5,28)-r2(5,32)-r2(5,36)*two)*xxx+(+r5(10,23)*two+r5(10,24) &
&                  -r4(6,32)*two-r4(6,33)-r4(6,36)*four-r4(6,37)*two+r3(3,33)*four &
&                  +r3(3,34)*two+r3(3,37)*two+r3(3,38)+r3(3,41)*six+r3(3,42)*three &
&                  -r2(1,27)*two-r2(1,28)-r2(1,31)*two-r2(1,32)-r2(1,35)*four-r2(1,36)*two) &
&                  *xxz+rxyz(1)*xxxz
      eri(6,1,5,3)=r211+(+r7(20,4)*two-r6(14,6)*two-r6(14,8)*four+r5(9,8)*four &
&                  +r5(9,10)*two+r5(9,12)*six-r4(5,9)*two-r4(5,11)*two-r4(5,13)*four)*qx &
&                  +rxyz(17)*qz+(+r6(20,12)-r5(14,16)-r5(14,20)*two+r4(9,21)*two+r4(9,25) &
&                  +r4(9,29)*three-r3(5,18)-r3(5,22)-r3(5,26)*two)*xx+(+r6(14,11)*two &
&                  -r5(9,15)*two-r5(9,19)*four+r4(5,20)*four+r4(5,24)*two+r4(5,28)*six &
&                  -r3(2,17)*two-r3(2,21)*two-r3(2,25)*four)*xz+rxyz(2)*xxz
      eri(1,2,5,3)=r220+(+r7(19,3)*two-r6(13,5)*two-r6(13,7)*four+r5(8,7)*four &
&                  +r5(8,9)*two+r5(8,11)*six+r5(10,21)*two-r4(4,8)*two-r4(4,10)*two &
&                  -r4(4,12)*four-r4(6,30)*two-r4(6,34)*four+r3(3,31)*four+r3(3,35)*two &
&                  +r3(3,39)*six-r2(1,25)*two-r2(1,29)*two-r2(1,33)*four)*qx+rxyz(4)*xx
      eri(2,2,5,3)=r040
      eri(3,2,5,3)=r022+(+r7(26,3)*two-r6(19,5)*two-r6(19,7)*four+r5(13,7)*four &
&                  +r5(13,9)*two+r5(13,11)*six+r5(15,21)*two-r4(8,8)*two-r4(8,10)*two &
&                  -r4(8,12)*four-r4(10,30)*two-r4(10,34)*four+r3(6,31)*four+r3(6,35)*two &
&                  +r3(6,39)*six-r2(5,25)*two-r2(5,29)*two-r2(5,33)*four)*qz+rxyz(4)*zz
      eri(4,2,5,3)=r130+rxyz(7)*qx
      eri(5,2,5,3)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,5,3)=r031+rxyz(7)*qz
      eri(1,3,5,3)=r202+(+r7(21,3)*two-r6(15,5)*two-r6(15,7)*four+r5(10,7)*four &
&                  +r5(10,9)*two+r5(10,11)*six+r5(10,21)*two-r4(6,8)*two-r4(6,10)*two &
&                  -r4(6,12)*four-r4(6,30)*two-r4(6,34)*four+r3(3,31)*four+r3(3,35)*two &
&                  +r3(3,39)*six-r2(1,25)*two-r2(1,29)*two-r2(1,33)*four)*qx+(+r7(15,4)*two &
&                  -r6(10,6)*two-r6(10,8)*four+r5(6,8)*four+r5(6,10)*two+r5(6,12)*six &
&                  +r5(15,22)*two-r4(3,9)*two-r4(3,11)*two-r4(3,13)*four-r4(10,31)*two &
&                  -r4(10,35)*four+r3(6,32)*four+r3(6,36)*two+r3(6,40)*six-r2(5,26)*two &
&                  -r2(5,30)*two-r2(5,34)*four)*qz+(+r6(21,10)-r5(15,14)-r5(15,18)*two &
&                  +r4(10,19)*two+r4(10,23)+r4(10,27)*three+r4(10,39)-r3(6,16)-r3(6,20) &
&                  -r3(6,24)*two-r3(6,44)-r3(6,49)*two+r2(5,43)*two+r2(5,48)+r2(5,53)*three &
&                  -r1(1,27)-r1(1,32)-r1(1,37)*two)*xx+rxyz(20)*xz+(+r6(10,12)-r5(6,16) &
&                  -r5(6,20)*two+r4(3,21)*two+r4(3,25)+r4(3,29)*three+r4(10,41)-r3(1,18) &
&                  -r3(1,22)-r3(1,26)*two-r3(6,46)-r3(6,51)*two+r2(5,45)*two+r2(5,50) &
&                  +r2(5,55)*three-r1(1,29)-r1(1,34)-r1(1,39)*two)*zz+(+r5(15,23)*two &
&                  -r4(10,32)*two-r4(10,36)*four+r3(6,33)*four+r3(6,37)*two+r3(6,41)*six &
&                  -r2(5,27)*two-r2(5,31)*two-r2(5,35)*four)*xxz+(+r5(10,24)*two-r4(6,33)*two &
&                  -r4(6,37)*four+r3(3,34)*four+r3(3,38)*two+r3(3,42)*six-r2(1,28)*two &
&                  -r2(1,32)*two-r2(1,36)*four)*xzz+rxyz(1)*xxzz
      eri(2,3,5,3)=r022+(+r7(26,4)*two-r6(19,6)*two-r6(19,8)*four+r5(13,8)*four &
&                  +r5(13,10)*two+r5(13,12)*six+r5(15,22)*two-r4(8,9)*two-r4(8,11)*two &
&                  -r4(8,13)*four-r4(10,31)*two-r4(10,35)*four+r3(6,32)*four+r3(6,36)*two &
&                  +r3(6,40)*six-r2(5,26)*two-r2(5,30)*two-r2(5,34)*four)*qz+rxyz(5)*zz
      eri(3,3,5,3)=r004+(+r7(28,3)*two+r7(28,4)*two-r6(21,5)*two-r6(21,6)*two &
&                  -r6(21,7)*four-r6(21,8)*four+r5(15,7)*four+r5(15,8)*four+r5(15,9)*two &
&                  +r5(15,10)*two+r5(15,11)*six+r5(15,12)*six+r5(15,21)*six+r5(15,22)*six &
&                  -r4(10,8)*two-r4(10,9)*two-r4(10,10)*two-r4(10,11)*two-r4(10,12)*four &
&                  -r4(10,13)*four-r4(10,30)*six-r4(10,31)*six-r4(10,34)*p12-r4(10,35)*p12 &
&                  +r3(6,31)*p12+r3(6,32)*p12+r3(6,35)*six+r3(6,36)*six+r3(6,39)*p18 &
&                  +r3(6,40)*p18-r2(5,25)*six-r2(5,26)*six-r2(5,29)*six-r2(5,30)*six &
&                  -r2(5,33)*p12-r2(5,34)*p12)*qz+(+r6(21,10)+r6(21,11)*four+r6(21,12) &
&                  -r5(15,14)-r5(15,15)*four-r5(15,16)-r5(15,18)*two-r5(15,19)*eight &
&                  -r5(15,20)*two+r4(10,19)*two+r4(10,20)*eight+r4(10,21)*two+r4(10,23) &
&                  +r4(10,24)*four+r4(10,25)+r4(10,27)*three+r4(10,28)*p12+r4(10,29)*three &
&                  +r4(10,39)+r4(10,40)*four+r4(10,41)-r3(6,16)-r3(6,17)*four-r3(6,18) &
&                  -r3(6,20)-r3(6,21)*four-r3(6,22)-r3(6,24)*two-r3(6,25)*eight-r3(6,26)*two &
&                  -r3(6,44)-r3(6,45)*four-r3(6,46)-r3(6,49)*two-r3(6,50)*eight-r3(6,51)*two &
&                  +r2(5,43)*two+r2(5,44)*eight+r2(5,45)*two+r2(5,48)+r2(5,49)*four+r2(5,50) &
&                  +r2(5,53)*three+r2(5,54)*p12+r2(5,55)*three-r1(1,27)-r1(1,28)*four &
&                  -r1(1,29)-r1(1,32)-r1(1,33)*four-r1(1,34)-r1(1,37)*two-r1(1,38)*eight &
&                  -r1(1,39)*two)*zz+(+r5(15,23)*two+r5(15,24)*two-r4(10,32)*two &
&                  -r4(10,33)*two-r4(10,36)*four-r4(10,37)*four+r3(6,33)*four+r3(6,34)*four &
&                  +r3(6,37)*two+r3(6,38)*two+r3(6,41)*six+r3(6,42)*six-r2(5,27)*two &
&                  -r2(5,28)*two-r2(5,31)*two-r2(5,32)*two-r2(5,35)*four-r2(5,36)*four)*zzz &
&                  +rxyz(1)*zzzz
      eri(4,3,5,3)=r112+rxyz(18)*qx+(+r7(20,4)*two-r6(14,6)*two-r6(14,8)*four &
&                  +r5(9,8)*four+r5(9,10)*two+r5(9,12)*six-r4(5,9)*two-r4(5,11)*two &
&                  -r4(5,13)*four)*qz+(+r6(20,11)*two-r5(14,15)*two-r5(14,19)*four &
&                  +r4(9,20)*four+r4(9,24)*two+r4(9,28)*six-r3(5,17)*two-r3(5,21)*two &
&                  -r3(5,25)*four)*xz+(+r6(14,12)-r5(9,16)-r5(9,20)*two+r4(5,21)*two+r4(5,25) &
&                  +r4(5,29)*three-r3(2,18)-r3(2,22)-r3(2,26)*two)*zz+rxyz(2)*xzz
      eri(5,3,5,3)=r103+(+r7(28,3)-r6(21,5)-r6(21,7)*two+r5(15,7)*two+r5(15,9) &
&                  +r5(15,11)*three+r5(15,21)*three-r4(10,8)-r4(10,10)-r4(10,12)*two &
&                  -r4(10,30)*three-r4(10,34)*six+r3(6,31)*six+r3(6,35)*three+r3(6,39)*nine &
&                  -r2(5,25)*three-r2(5,29)*three-r2(5,33)*six)*qx+(+r7(21,3)+r7(21,4)*two &
&                  -r6(15,5)-r6(15,6)*two-r6(15,7)*two-r6(15,8)*four+r5(10,7)*two &
&                  +r5(10,8)*four+r5(10,9)+r5(10,10)*two+r5(10,11)*three+r5(10,12)*six &
&                  +r5(10,21)+r5(10,22)*two-r4(6,8)-r4(6,9)*two-r4(6,10)-r4(6,11)*two &
&                  -r4(6,12)*two-r4(6,13)*four-r4(6,30)-r4(6,31)*two-r4(6,34)*two &
&                  -r4(6,35)*four+r3(3,31)*two+r3(3,32)*four+r3(3,35)+r3(3,36)*two &
&                  +r3(3,39)*three+r3(3,40)*six-r2(1,25)-r2(1,26)*two-r2(1,29)-r2(1,30)*two &
&                  -r2(1,33)*two-r2(1,34)*four)*qz+(+r6(21,10)+r6(21,11)*two-r5(15,14) &
&                  -r5(15,15)*two-r5(15,18)*two-r5(15,19)*four+r4(10,19)*two+r4(10,20)*four &
&                  +r4(10,23)+r4(10,24)*two+r4(10,27)*three+r4(10,28)*six+r4(10,39) &
&                  +r4(10,40)*two-r3(6,16)-r3(6,17)*two-r3(6,20)-r3(6,21)*two-r3(6,24)*two &
&                  -r3(6,25)*four-r3(6,44)-r3(6,45)*two-r3(6,49)*two-r3(6,50)*four &
&                  +r2(5,43)*two+r2(5,44)*four+r2(5,48)+r2(5,49)*two+r2(5,53)*three &
&                  +r2(5,54)*six-r1(1,27)-r1(1,28)*two-r1(1,32)-r1(1,33)*two-r1(1,37)*two &
&                  -r1(1,38)*four)*xz+(+r6(15,11)*two+r6(15,12)-r5(10,15)*two-r5(10,16) &
&                  -r5(10,19)*four-r5(10,20)*two+r4(6,20)*four+r4(6,21)*two+r4(6,24)*two &
&                  +r4(6,25)+r4(6,28)*six+r4(6,29)*three-r3(3,17)*two-r3(3,18)-r3(3,21)*two &
&                  -r3(3,22)-r3(3,25)*four-r3(3,26)*two)*zz+(+r5(15,23)*two+r5(15,24) &
&                  -r4(10,32)*two-r4(10,33)-r4(10,36)*four-r4(10,37)*two+r3(6,33)*four &
&                  +r3(6,34)*two+r3(6,37)*two+r3(6,38)+r3(6,41)*six+r3(6,42)*three &
&                  -r2(5,27)*two-r2(5,28)-r2(5,31)*two-r2(5,32)-r2(5,35)*four-r2(5,36)*two) &
&                  *xzz+(+r5(10,24)-r4(6,33)-r4(6,37)*two+r3(3,34)*two+r3(3,38) &
&                  +r3(3,42)*three-r2(1,28)-r2(1,32)-r2(1,36)*two)*zzz+rxyz(1)*xzzz
      eri(6,3,5,3)=r013+(+r7(27,3)+r7(27,4)*two-r6(20,5)-r6(20,6)*two-r6(20,7)*two &
&                  -r6(20,8)*four+r5(14,7)*two+r5(14,8)*four+r5(14,9)+r5(14,10)*two &
&                  +r5(14,11)*three+r5(14,12)*six+r5(14,21)+r5(14,22)*two-r4(9,8)-r4(9,9)*two &
&                  -r4(9,10)-r4(9,11)*two-r4(9,12)*two-r4(9,13)*four-r4(9,30)-r4(9,31)*two &
&                  -r4(9,34)*two-r4(9,35)*four+r3(5,31)*two+r3(5,32)*four+r3(5,35) &
&                  +r3(5,36)*two+r3(5,39)*three+r3(5,40)*six-r2(4,25)-r2(4,26)*two-r2(4,29) &
&                  -r2(4,30)*two-r2(4,33)*two-r2(4,34)*four)*qz+(+r6(20,11)*two+r6(20,12) &
&                  -r5(14,15)*two-r5(14,16)-r5(14,19)*four-r5(14,20)*two+r4(9,20)*four &
&                  +r4(9,21)*two+r4(9,24)*two+r4(9,25)+r4(9,28)*six+r4(9,29)*three &
&                  -r3(5,17)*two-r3(5,18)-r3(5,21)*two-r3(5,22)-r3(5,25)*four-r3(5,26)*two) &
&                  *zz+rxyz(2)*zzz
      eri(1,4,5,3)=r310+(+r7(14,3)*two+r7(14,4)-r6(9,5)*two-r6(9,6)-r6(9,7)*four &
&                  -r6(9,8)*two+r5(5,7)*four+r5(5,8)*two+r5(5,9)*two+r5(5,10)+r5(5,11)*six &
&                  +r5(5,12)*three+r5(14,21)*two+r5(14,22)-r4(2,8)*two-r4(2,9)-r4(2,10)*two &
&                  -r4(2,11)-r4(2,12)*four-r4(2,13)*two-r4(9,30)*two-r4(9,31)-r4(9,34)*four &
&                  -r4(9,35)*two+r3(5,31)*four+r3(5,32)*two+r3(5,35)*two+r3(5,36) &
&                  +r3(5,39)*six+r3(5,40)*three-r2(4,25)*two-r2(4,26)-r2(4,29)*two-r2(4,30) &
&                  -r2(4,33)*four-r2(4,34)*two)*qx+(+r6(14,10)+r6(14,11)*two-r5(9,14) &
&                  -r5(9,15)*two-r5(9,18)*two-r5(9,19)*four+r4(5,19)*two+r4(5,20)*four &
&                  +r4(5,23)+r4(5,24)*two+r4(5,27)*three+r4(5,28)*six-r3(2,16)-r3(2,17)*two &
&                  -r3(2,20)-r3(2,21)*two-r3(2,24)*two-r3(2,25)*four)*xx+rxyz(3)*xxx
      eri(2,4,5,3)=r130+rxyz(8)*qx
      eri(3,4,5,3)=r112+rxyz(19)*qx+(+r7(20,3)*two-r6(14,5)*two-r6(14,7)*four &
&                  +r5(9,7)*four+r5(9,9)*two+r5(9,11)*six-r4(5,8)*two-r4(5,10)*two &
&                  -r4(5,12)*four)*qz+(+r6(20,11)*two-r5(14,15)*two-r5(14,19)*four &
&                  +r4(9,20)*four+r4(9,24)*two+r4(9,28)*six-r3(5,17)*two-r3(5,21)*two &
&                  -r3(5,25)*four)*xz+(+r6(14,10)-r5(9,14)-r5(9,18)*two+r4(5,19)*two+r4(5,23) &
&                  +r4(5,27)*three-r3(2,16)-r3(2,20)-r3(2,24)*two)*zz+rxyz(3)*xzz
      eri(4,4,5,3)=r220+(+r7(19,3)+r7(19,4)-r6(13,5)-r6(13,6)-r6(13,7)*two-r6(13,8)*two &
&                  +r5(8,7)*two+r5(8,8)*two+r5(8,9)+r5(8,10)+r5(8,11)*three+r5(8,12)*three &
&                  +r5(10,21)+r5(10,22)-r4(4,8)-r4(4,9)-r4(4,10)-r4(4,11)-r4(4,12)*two &
&                  -r4(4,13)*two-r4(6,30)-r4(6,31)-r4(6,34)*two-r4(6,35)*two+r3(3,31)*two &
&                  +r3(3,32)*two+r3(3,35)+r3(3,36)+r3(3,39)*three+r3(3,40)*three-r2(1,25) &
&                  -r2(1,26)-r2(1,29)-r2(1,30)-r2(1,33)*two-r2(1,34)*two)*qx+rxyz(6)*xx
      eri(5,4,5,3)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(14,10)+r6(14,11)-r5(9,14) &
&                  -r5(9,15)-r5(9,18)*two-r5(9,19)*two+r4(5,19)*two+r4(5,20)*two+r4(5,23) &
&                  +r4(5,24)+r4(5,27)*three+r4(5,28)*three-r3(2,16)-r3(2,17)-r3(2,20) &
&                  -r3(2,21)-r3(2,24)*two-r3(2,25)*two)*xz+rxyz(3)*xxz
      eri(6,4,5,3)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,5,3)=r301+(+r7(15,3)*two+r7(15,4)-r6(10,5)*two-r6(10,6)-r6(10,7)*four &
&                  -r6(10,8)*two+r5(6,7)*four+r5(6,8)*two+r5(6,9)*two+r5(6,10)+r5(6,11)*six &
&                  +r5(6,12)*three+r5(15,21)*two+r5(15,22)-r4(3,8)*two-r4(3,9)-r4(3,10)*two &
&                  -r4(3,11)-r4(3,12)*four-r4(3,13)*two-r4(10,30)*two-r4(10,31) &
&                  -r4(10,34)*four-r4(10,35)*two+r3(6,31)*four+r3(6,32)*two+r3(6,35)*two &
&                  +r3(6,36)+r3(6,39)*six+r3(6,40)*three-r2(5,25)*two-r2(5,26)-r2(5,29)*two &
&                  -r2(5,30)-r2(5,33)*four-r2(5,34)*two)*qx+(+r7(10,4)-r6(6,6)-r6(6,8)*two &
&                  +r5(3,8)*two+r5(3,10)+r5(3,12)*three+r5(10,22)*three-r4(1,9)-r4(1,11) &
&                  -r4(1,13)*two-r4(6,31)*three-r4(6,35)*six+r3(3,32)*six+r3(3,36)*three &
&                  +r3(3,40)*nine-r2(1,26)*three-r2(1,30)*three-r2(1,34)*six)*qz+(+r6(15,10) &
&                  +r6(15,11)*two-r5(10,14)-r5(10,15)*two-r5(10,18)*two-r5(10,19)*four &
&                  +r4(6,19)*two+r4(6,20)*four+r4(6,23)+r4(6,24)*two+r4(6,27)*three &
&                  +r4(6,28)*six-r3(3,16)-r3(3,17)*two-r3(3,20)-r3(3,21)*two-r3(3,24)*two &
&                  -r3(3,25)*four)*xx+(+r6(10,11)*two+r6(10,12)-r5(6,15)*two-r5(6,16) &
&                  -r5(6,19)*four-r5(6,20)*two+r4(3,20)*four+r4(3,21)*two+r4(3,24)*two &
&                  +r4(3,25)+r4(3,28)*six+r4(3,29)*three+r4(10,40)*two+r4(10,41)-r3(1,17)*two &
&                  -r3(1,18)-r3(1,21)*two-r3(1,22)-r3(1,25)*four-r3(1,26)*two-r3(6,45)*two &
&                  -r3(6,46)-r3(6,50)*four-r3(6,51)*two+r2(5,44)*four+r2(5,45)*two &
&                  +r2(5,49)*two+r2(5,50)+r2(5,54)*six+r2(5,55)*three-r1(1,28)*two-r1(1,29) &
&                  -r1(1,33)*two-r1(1,34)-r1(1,38)*four-r1(1,39)*two)*xz+(+r5(15,23) &
&                  -r4(10,32)-r4(10,36)*two+r3(6,33)*two+r3(6,37)+r3(6,41)*three-r2(5,27) &
&                  -r2(5,31)-r2(5,35)*two)*xxx+(+r5(10,23)+r5(10,24)*two-r4(6,32) &
&                  -r4(6,33)*two-r4(6,36)*two-r4(6,37)*four+r3(3,33)*two+r3(3,34)*four &
&                  +r3(3,37)+r3(3,38)*two+r3(3,41)*three+r3(3,42)*six-r2(1,27)-r2(1,28)*two &
&                  -r2(1,31)-r2(1,32)*two-r2(1,35)*two-r2(1,36)*four)*xxz+rxyz(1)*xxxz
      eri(2,5,5,3)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,5,3)=r103+(+r7(28,4)-r6(21,6)-r6(21,8)*two+r5(15,8)*two+r5(15,10) &
&                  +r5(15,12)*three+r5(15,22)*three-r4(10,9)-r4(10,11)-r4(10,13)*two &
&                  -r4(10,31)*three-r4(10,35)*six+r3(6,32)*six+r3(6,36)*three+r3(6,40)*nine &
&                  -r2(5,26)*three-r2(5,30)*three-r2(5,34)*six)*qx+(+r7(21,3)*two+r7(21,4) &
&                  -r6(15,5)*two-r6(15,6)-r6(15,7)*four-r6(15,8)*two+r5(10,7)*four &
&                  +r5(10,8)*two+r5(10,9)*two+r5(10,10)+r5(10,11)*six+r5(10,12)*three &
&                  +r5(10,21)*two+r5(10,22)-r4(6,8)*two-r4(6,9)-r4(6,10)*two-r4(6,11) &
&                  -r4(6,12)*four-r4(6,13)*two-r4(6,30)*two-r4(6,31)-r4(6,34)*four &
&                  -r4(6,35)*two+r3(3,31)*four+r3(3,32)*two+r3(3,35)*two+r3(3,36) &
&                  +r3(3,39)*six+r3(3,40)*three-r2(1,25)*two-r2(1,26)-r2(1,29)*two-r2(1,30) &
&                  -r2(1,33)*four-r2(1,34)*two)*qz+(+r6(21,11)*two+r6(21,12)-r5(15,15)*two &
&                  -r5(15,16)-r5(15,19)*four-r5(15,20)*two+r4(10,20)*four+r4(10,21)*two &
&                  +r4(10,24)*two+r4(10,25)+r4(10,28)*six+r4(10,29)*three+r4(10,40)*two &
&                  +r4(10,41)-r3(6,17)*two-r3(6,18)-r3(6,21)*two-r3(6,22)-r3(6,25)*four &
&                  -r3(6,26)*two-r3(6,45)*two-r3(6,46)-r3(6,50)*four-r3(6,51)*two &
&                  +r2(5,44)*four+r2(5,45)*two+r2(5,49)*two+r2(5,50)+r2(5,54)*six &
&                  +r2(5,55)*three-r1(1,28)*two-r1(1,29)-r1(1,33)*two-r1(1,34)-r1(1,38)*four &
&                  -r1(1,39)*two)*xz+(+r6(15,10)+r6(15,11)*two-r5(10,14)-r5(10,15)*two &
&                  -r5(10,18)*two-r5(10,19)*four+r4(6,19)*two+r4(6,20)*four+r4(6,23) &
&                  +r4(6,24)*two+r4(6,27)*three+r4(6,28)*six-r3(3,16)-r3(3,17)*two-r3(3,20) &
&                  -r3(3,21)*two-r3(3,24)*two-r3(3,25)*four)*zz+(+r5(15,23)+r5(15,24)*two &
&                  -r4(10,32)-r4(10,33)*two-r4(10,36)*two-r4(10,37)*four+r3(6,33)*two &
&                  +r3(6,34)*four+r3(6,37)+r3(6,38)*two+r3(6,41)*three+r3(6,42)*six-r2(5,27) &
&                  -r2(5,28)*two-r2(5,31)-r2(5,32)*two-r2(5,35)*two-r2(5,36)*four)*xzz+( &
&                  +r5(10,23)-r4(6,32)-r4(6,36)*two+r3(3,33)*two+r3(3,37)+r3(3,41)*three &
&                  -r2(1,27)-r2(1,31)-r2(1,35)*two)*zzz+rxyz(1)*xzzz
      eri(4,5,5,3)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(14,11)+r6(14,12)-r5(9,15) &
&                  -r5(9,16)-r5(9,19)*two-r5(9,20)*two+r4(5,20)*two+r4(5,21)*two+r4(5,24) &
&                  +r4(5,25)+r4(5,28)*three+r4(5,29)*three-r3(2,17)-r3(2,18)-r3(2,21) &
&                  -r3(2,22)-r3(2,25)*two-r3(2,26)*two)*xz+rxyz(2)*xxz
      eri(5,5,5,3)=r202+(+r7(21,3)+r7(21,4)-r6(15,5)-r6(15,6)-r6(15,7)*two-r6(15,8)*two &
&                  +r5(10,7)*two+r5(10,8)*two+r5(10,9)+r5(10,10)+r5(10,11)*three &
&                  +r5(10,12)*three+r5(10,21)+r5(10,22)-r4(6,8)-r4(6,9)-r4(6,10)-r4(6,11) &
&                  -r4(6,12)*two-r4(6,13)*two-r4(6,30)-r4(6,31)-r4(6,34)*two-r4(6,35)*two &
&                  +r3(3,31)*two+r3(3,32)*two+r3(3,35)+r3(3,36)+r3(3,39)*three+r3(3,40)*three &
&                  -r2(1,25)-r2(1,26)-r2(1,29)-r2(1,30)-r2(1,33)*two-r2(1,34)*two)*qx+( &
&                  +r7(15,3)+r7(15,4)-r6(10,5)-r6(10,6)-r6(10,7)*two-r6(10,8)*two+r5(6,7)*two &
&                  +r5(6,8)*two+r5(6,9)+r5(6,10)+r5(6,11)*three+r5(6,12)*three+r5(15,21) &
&                  +r5(15,22)-r4(3,8)-r4(3,9)-r4(3,10)-r4(3,11)-r4(3,12)*two-r4(3,13)*two &
&                  -r4(10,30)-r4(10,31)-r4(10,34)*two-r4(10,35)*two+r3(6,31)*two+r3(6,32)*two &
&                  +r3(6,35)+r3(6,36)+r3(6,39)*three+r3(6,40)*three-r2(5,25)-r2(5,26) &
&                  -r2(5,29)-r2(5,30)-r2(5,33)*two-r2(5,34)*two)*qz+(+r6(21,11)-r5(15,15) &
&                  -r5(15,19)*two+r4(10,20)*two+r4(10,24)+r4(10,28)*three+r4(10,40)-r3(6,17) &
&                  -r3(6,21)-r3(6,25)*two-r3(6,45)-r3(6,50)*two+r2(5,44)*two+r2(5,49) &
&                  +r2(5,54)*three-r1(1,28)-r1(1,33)-r1(1,38)*two)*xx+(+r6(15,10) &
&                  +r6(15,11)*two+r6(15,12)-r5(10,14)-r5(10,15)*two-r5(10,16)-r5(10,18)*two &
&                  -r5(10,19)*four-r5(10,20)*two+r4(6,19)*two+r4(6,20)*four+r4(6,21)*two &
&                  +r4(6,23)+r4(6,24)*two+r4(6,25)+r4(6,27)*three+r4(6,28)*six+r4(6,29)*three &
&                  -r3(3,16)-r3(3,17)*two-r3(3,18)-r3(3,20)-r3(3,21)*two-r3(3,22) &
&                  -r3(3,24)*two-r3(3,25)*four-r3(3,26)*two)*xz+(+r6(10,11)-r5(6,15) &
&                  -r5(6,19)*two+r4(3,20)*two+r4(3,24)+r4(3,28)*three+r4(10,40)-r3(1,17) &
&                  -r3(1,21)-r3(1,25)*two-r3(6,45)-r3(6,50)*two+r2(5,44)*two+r2(5,49) &
&                  +r2(5,54)*three-r1(1,28)-r1(1,33)-r1(1,38)*two)*zz+(+r5(15,23)+r5(15,24) &
&                  -r4(10,32)-r4(10,33)-r4(10,36)*two-r4(10,37)*two+r3(6,33)*two+r3(6,34)*two &
&                  +r3(6,37)+r3(6,38)+r3(6,41)*three+r3(6,42)*three-r2(5,27)-r2(5,28) &
&                  -r2(5,31)-r2(5,32)-r2(5,35)*two-r2(5,36)*two)*xxz+(+r5(10,23)+r5(10,24) &
&                  -r4(6,32)-r4(6,33)-r4(6,36)*two-r4(6,37)*two+r3(3,33)*two+r3(3,34)*two &
&                  +r3(3,37)+r3(3,38)+r3(3,41)*three+r3(3,42)*three-r2(1,27)-r2(1,28) &
&                  -r2(1,31)-r2(1,32)-r2(1,35)*two-r2(1,36)*two)*xzz+rxyz(1)*xxzz
      eri(6,5,5,3)=r112+rxyz(19)*qx+(+r7(20,3)+r7(20,4)-r6(14,5)-r6(14,6)-r6(14,7)*two &
&                  -r6(14,8)*two+r5(9,7)*two+r5(9,8)*two+r5(9,9)+r5(9,10)+r5(9,11)*three &
&                  +r5(9,12)*three-r4(5,8)-r4(5,9)-r4(5,10)-r4(5,11)-r4(5,12)*two &
&                  -r4(5,13)*two)*qz+(+r6(20,11)+r6(20,12)-r5(14,15)-r5(14,16)-r5(14,19)*two &
&                  -r5(14,20)*two+r4(9,20)*two+r4(9,21)*two+r4(9,24)+r4(9,25)+r4(9,28)*three &
&                  +r4(9,29)*three-r3(5,17)-r3(5,18)-r3(5,21)-r3(5,22)-r3(5,25)*two &
&                  -r3(5,26)*two)*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,5,3)=r211+(+r7(20,3)*two-r6(14,5)*two-r6(14,7)*four+r5(9,7)*four &
&                  +r5(9,9)*two+r5(9,11)*six-r4(5,8)*two-r4(5,10)*two-r4(5,12)*four)*qx &
&                  +rxyz(16)*qz+(+r6(20,10)-r5(14,14)-r5(14,18)*two+r4(9,19)*two+r4(9,23) &
&                  +r4(9,27)*three-r3(5,16)-r3(5,20)-r3(5,24)*two)*xx+(+r6(14,11)*two &
&                  -r5(9,15)*two-r5(9,19)*four+r4(5,20)*four+r4(5,24)*two+r4(5,28)*six &
&                  -r3(2,17)*two-r3(2,21)*two-r3(2,25)*four)*xz+rxyz(3)*xxz
      eri(2,6,5,3)=r031+rxyz(8)*qz
      eri(3,6,5,3)=r013+(+r7(27,3)*two+r7(27,4)-r6(20,5)*two-r6(20,6)-r6(20,7)*four &
&                  -r6(20,8)*two+r5(14,7)*four+r5(14,8)*two+r5(14,9)*two+r5(14,10) &
&                  +r5(14,11)*six+r5(14,12)*three+r5(14,21)*two+r5(14,22)-r4(9,8)*two-r4(9,9) &
&                  -r4(9,10)*two-r4(9,11)-r4(9,12)*four-r4(9,13)*two-r4(9,30)*two-r4(9,31) &
&                  -r4(9,34)*four-r4(9,35)*two+r3(5,31)*four+r3(5,32)*two+r3(5,35)*two &
&                  +r3(5,36)+r3(5,39)*six+r3(5,40)*three-r2(4,25)*two-r2(4,26)-r2(4,29)*two &
&                  -r2(4,30)-r2(4,33)*four-r2(4,34)*two)*qz+(+r6(20,10)+r6(20,11)*two &
&                  -r5(14,14)-r5(14,15)*two-r5(14,18)*two-r5(14,19)*four+r4(9,19)*two &
&                  +r4(9,20)*four+r4(9,23)+r4(9,24)*two+r4(9,27)*three+r4(9,28)*six-r3(5,16) &
&                  -r3(5,17)*two-r3(5,20)-r3(5,21)*two-r3(5,24)*two-r3(5,25)*four)*zz+rxyz(3) &
&                  *zzz
      eri(4,6,5,3)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,5,3)=r112+rxyz(18)*qx+(+r7(20,3)+r7(20,4)-r6(14,5)-r6(14,6)-r6(14,7)*two &
&                  -r6(14,8)*two+r5(9,7)*two+r5(9,8)*two+r5(9,9)+r5(9,10)+r5(9,11)*three &
&                  +r5(9,12)*three-r4(5,8)-r4(5,9)-r4(5,10)-r4(5,11)-r4(5,12)*two &
&                  -r4(5,13)*two)*qz+(+r6(20,10)+r6(20,11)-r5(14,14)-r5(14,15)-r5(14,18)*two &
&                  -r5(14,19)*two+r4(9,19)*two+r4(9,20)*two+r4(9,23)+r4(9,24)+r4(9,27)*three &
&                  +r4(9,28)*three-r3(5,16)-r3(5,17)-r3(5,20)-r3(5,21)-r3(5,24)*two &
&                  -r3(5,25)*two)*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,5,3)=r022+(+r7(26,3)+r7(26,4)-r6(19,5)-r6(19,6)-r6(19,7)*two-r6(19,8)*two &
&                  +r5(13,7)*two+r5(13,8)*two+r5(13,9)+r5(13,10)+r5(13,11)*three &
&                  +r5(13,12)*three+r5(15,21)+r5(15,22)-r4(8,8)-r4(8,9)-r4(8,10)-r4(8,11) &
&                  -r4(8,12)*two-r4(8,13)*two-r4(10,30)-r4(10,31)-r4(10,34)*two-r4(10,35)*two &
&                  +r3(6,31)*two+r3(6,32)*two+r3(6,35)+r3(6,36)+r3(6,39)*three+r3(6,40)*three &
&                  -r2(5,25)-r2(5,26)-r2(5,29)-r2(5,30)-r2(5,33)*two-r2(5,34)*two)*qz+rxyz(6) &
&                  *zz
!
      r400= r8(14)-r7(9,1)-r7(9,2)*two+r6(5,2)*two+r6(5,3)+r6(5,4)*three+r6(14,9)*six &
&          -r5(2,2)-r5(2,3)-r5(2,4)*two-r5(9,13)*six-r5(9,17)*p12+r4(5,18)*p12 &
&          +r4(5,22)*six+r4(5,26)*p18+r4(14,38)*three-r3(2,15)*six-r3(2,19)*six &
&          -r3(2,23)*p12-r3(9,43)*three-r3(9,48)*six+r2(6,42)*six+r2(6,47)*three &
&          +r2(6,52)*nine-r1(2,26)*three-r1(2,31)*three-r1(2,36)*six
      r310= r8(19)-r7(13,1)-r7(13,2)*two+r6(8,2)*two+r6(8,3)+r6(8,4)*three+r6(19,9)*three &
&          -r5(4,2)-r5(4,3)-r5(4,4)*two-r5(13,13)*three-r5(13,17)*six+r4(8,18)*six &
&          +r4(8,22)*three+r4(8,26)*nine-r3(4,15)*three-r3(4,19)*three-r3(4,23)*six
      r301= r8(20)-r7(14,1)-r7(14,2)*two+r6(9,2)*two+r6(9,3)+r6(9,4)*three+r6(20,9)*three &
&          -r5(5,2)-r5(5,3)-r5(5,4)*two-r5(14,13)*three-r5(14,17)*six+r4(9,18)*six &
&          +r4(9,22)*three+r4(9,26)*nine-r3(5,15)*three-r3(5,19)*three-r3(5,23)*six
      r220= r8(25)-r7(18,1)-r7(18,2)*two+r6(12,2)*two+r6(12,3)+r6(12,4)*three+r6(14,9) &
&          +r6(25,9)-r5(7,2)-r5(7,3)-r5(7,4)*two-r5(9,13)-r5(18,13)-r5(9,17)*two &
&          -r5(18,17)*two+r4(5,18)*two+r4(12,18)*two+r4(5,22)+r4(12,22)+r4(5,26)*three &
&          +r4(12,26)*three+r4(14,38)-r3(2,15)-r3(7,15)-r3(2,19)-r3(7,19)-r3(2,23)*two &
&          -r3(7,23)*two-r3(9,43)-r3(9,48)*two+r2(6,42)*two+r2(6,47)+r2(6,52)*three &
&          -r1(2,26)-r1(2,31)-r1(2,36)*two
      r211= r8(26)-r7(19,1)-r7(19,2)*two+r6(13,2)*two+r6(13,3)+r6(13,4)*three+r6(26,9) &
&          -r5(8,2)-r5(8,3)-r5(8,4)*two-r5(19,13)-r5(19,17)*two+r4(13,18)*two+r4(13,22) &
&          +r4(13,26)*three-r3(8,15)-r3(8,19)-r3(8,23)*two
      r202= r8(27)-r7(20,1)-r7(20,2)*two+r6(14,2)*two+r6(14,3)+r6(14,4)*three+r6(14,9) &
&          +r6(27,9)-r5(9,2)-r5(9,3)-r5(9,4)*two-r5(9,13)-r5(20,13)-r5(9,17)*two &
&          -r5(20,17)*two+r4(5,18)*two+r4(14,18)*two+r4(5,22)+r4(14,22)+r4(5,26)*three &
&          +r4(14,26)*three+r4(14,38)-r3(2,15)-r3(9,15)-r3(2,19)-r3(9,19)-r3(2,23)*two &
&          -r3(9,23)*two-r3(9,43)-r3(9,48)*two+r2(6,42)*two+r2(6,47)+r2(6,52)*three &
&          -r1(2,26)-r1(2,31)-r1(2,36)*two
      r130= r8(32)-r7(24,1)-r7(24,2)*two+r6(17,2)*two+r6(17,3)+r6(17,4)*three &
&          +r6(19,9)*three-r5(11,2)-r5(11,3)-r5(11,4)*two-r5(13,13)*three-r5(13,17)*six &
&          +r4(8,18)*six+r4(8,22)*three+r4(8,26)*nine-r3(4,15)*three-r3(4,19)*three &
&          -r3(4,23)*six
      r121= r8(33)-r7(25,1)-r7(25,2)*two+r6(18,2)*two+r6(18,3)+r6(18,4)*three+r6(20,9) &
&          -r5(12,2)-r5(12,3)-r5(12,4)*two-r5(14,13)-r5(14,17)*two+r4(9,18)*two+r4(9,22) &
&          +r4(9,26)*three-r3(5,15)-r3(5,19)-r3(5,23)*two
      r112= r8(34)-r7(26,1)-r7(26,2)*two+r6(19,2)*two+r6(19,3)+r6(19,4)*three+r6(19,9) &
&          -r5(13,2)-r5(13,3)-r5(13,4)*two-r5(13,13)-r5(13,17)*two+r4(8,18)*two+r4(8,22) &
&          +r4(8,26)*three-r3(4,15)-r3(4,19)-r3(4,23)*two
      r103= r8(35)-r7(27,1)-r7(27,2)*two+r6(20,2)*two+r6(20,3)+r6(20,4)*three &
&          +r6(20,9)*three-r5(14,2)-r5(14,3)-r5(14,4)*two-r5(14,13)*three-r5(14,17)*six &
&          +r4(9,18)*six+r4(9,22)*three+r4(9,26)*nine-r3(5,15)*three-r3(5,19)*three &
&          -r3(5,23)*six
      r040= r8(40)-r7(31,1)-r7(31,2)*two+r6(23,2)*two+r6(23,3)+r6(23,4)*three+r6(25,9)*six &
&          -r5(16,2)-r5(16,3)-r5(16,4)*two-r5(18,13)*six-r5(18,17)*p12+r4(12,18)*p12 &
&          +r4(12,22)*six+r4(12,26)*p18+r4(14,38)*three-r3(7,15)*six-r3(7,19)*six &
&          -r3(7,23)*p12-r3(9,43)*three-r3(9,48)*six+r2(6,42)*six+r2(6,47)*three &
&          +r2(6,52)*nine-r1(2,26)*three-r1(2,31)*three-r1(2,36)*six
      r031= r8(41)-r7(32,1)-r7(32,2)*two+r6(24,2)*two+r6(24,3)+r6(24,4)*three &
&          +r6(26,9)*three-r5(17,2)-r5(17,3)-r5(17,4)*two-r5(19,13)*three-r5(19,17)*six &
&          +r4(13,18)*six+r4(13,22)*three+r4(13,26)*nine-r3(8,15)*three-r3(8,19)*three &
&          -r3(8,23)*six
      r022= r8(42)-r7(33,1)-r7(33,2)*two+r6(25,2)*two+r6(25,3)+r6(25,4)*three+r6(25,9) &
&          +r6(27,9)-r5(18,2)-r5(18,3)-r5(18,4)*two-r5(18,13)-r5(20,13)-r5(18,17)*two &
&          -r5(20,17)*two+r4(12,18)*two+r4(14,18)*two+r4(12,22)+r4(14,22)+r4(12,26)*three &
&          +r4(14,26)*three+r4(14,38)-r3(7,15)-r3(9,15)-r3(7,19)-r3(9,19)-r3(7,23)*two &
&          -r3(9,23)*two-r3(9,43)-r3(9,48)*two+r2(6,42)*two+r2(6,47)+r2(6,52)*three &
&          -r1(2,26)-r1(2,31)-r1(2,36)*two
      r013= r8(43)-r7(34,1)-r7(34,2)*two+r6(26,2)*two+r6(26,3)+r6(26,4)*three &
&          +r6(26,9)*three-r5(19,2)-r5(19,3)-r5(19,4)*two-r5(19,13)*three-r5(19,17)*six &
&          +r4(13,18)*six+r4(13,22)*three+r4(13,26)*nine-r3(8,15)*three-r3(8,19)*three &
&          -r3(8,23)*six
      r004= r8(44)-r7(35,1)-r7(35,2)*two+r6(27,2)*two+r6(27,3)+r6(27,4)*three+r6(27,9)*six &
&          -r5(20,2)-r5(20,3)-r5(20,4)*two-r5(20,13)*six-r5(20,17)*p12+r4(14,18)*p12 &
&          +r4(14,22)*six+r4(14,26)*p18+r4(14,38)*three-r3(9,15)*six-r3(9,19)*six &
&          -r3(9,23)*p12-r3(9,43)*three-r3(9,48)*six+r2(6,42)*six+r2(6,47)*three &
&          +r2(6,52)*nine-r1(2,26)*three-r1(2,31)*three-r1(2,36)*six
      rxyz(1)=+r4(14,42)-r3(9,47)-r3(9,52)*two+r2(6,46)*two+r2(6,51)+r2(6,56)*three &
&             -r1(2,30)-r1(2,35)-r1(2,40)*two
      rxyz(2)=+r5(19,24)-r4(13,33)-r4(13,37)*two+r3(8,34)*two+r3(8,38)+r3(8,42)*three &
&             -r2(2,28)-r2(2,32)-r2(2,36)*two
      rxyz(3)=+r5(19,23)-r4(13,32)-r4(13,36)*two+r3(8,33)*two+r3(8,37)+r3(8,41)*three &
&             -r2(2,27)-r2(2,31)-r2(2,35)*two
      rxyz(4)=+r6(25,10)-r5(18,14)-r5(18,18)*two+r4(12,19)*two+r4(12,23)+r4(12,27)*three &
&             +r4(14,39)-r3(7,16)-r3(7,20)-r3(7,24)*two-r3(9,44)-r3(9,49)*two+r2(6,43)*two &
&             +r2(6,48)+r2(6,53)*three-r1(2,27)-r1(2,32)-r1(2,37)*two
      rxyz(5)=+r6(25,12)-r5(18,16)-r5(18,20)*two+r4(12,21)*two+r4(12,25)+r4(12,29)*three &
&             +r4(14,41)-r3(7,18)-r3(7,22)-r3(7,26)*two-r3(9,46)-r3(9,51)*two+r2(6,45)*two &
&             +r2(6,50)+r2(6,55)*three-r1(2,29)-r1(2,34)-r1(2,39)*two
      rxyz(6)=+r6(25,11)-r5(18,15)-r5(18,19)*two+r4(12,20)*two+r4(12,24)+r4(12,28)*three &
&             +r4(14,40)-r3(7,17)-r3(7,21)-r3(7,25)*two-r3(9,45)-r3(9,50)*two+r2(6,44)*two &
&             +r2(6,49)+r2(6,54)*three-r1(2,28)-r1(2,33)-r1(2,38)*two
      rxyz(7)=+r7(32,3)-r6(24,5)-r6(24,7)*two+r5(17,7)*two+r5(17,9)+r5(17,11)*three &
&             +r5(19,21)*three-r4(11,8)-r4(11,10)-r4(11,12)*two-r4(13,30)*three &
&             -r4(13,34)*six+r3(8,31)*six+r3(8,35)*three+r3(8,39)*nine-r2(2,25)*three &
&             -r2(2,29)*three-r2(2,33)*six
      rxyz(8)=+r7(32,4)-r6(24,6)-r6(24,8)*two+r5(17,8)*two+r5(17,10)+r5(17,12)*three &
&             +r5(19,22)*three-r4(11,9)-r4(11,11)-r4(11,13)*two-r4(13,31)*three &
&             -r4(13,35)*six+r3(8,32)*six+r3(8,36)*three+r3(8,40)*nine-r2(2,26)*three &
&             -r2(2,30)*three-r2(2,34)*six
      rxyz(9)=+r7(26,3)+r7(26,4)-r6(19,5)-r6(19,6)-r6(19,7)*two-r6(19,8)*two+r5(13,7)*two &
&             +r5(13,8)*two+r5(13,9)+r5(13,10)+r5(13,11)*three+r5(13,12)*three-r4(8,8) &
&             -r4(8,9)-r4(8,10)-r4(8,11)-r4(8,12)*two-r4(8,13)*two
      rxyz(10)=+r6(26,11)-r5(19,15)-r5(19,19)*two+r4(13,20)*two+r4(13,24)+r4(13,28)*three &
&             -r3(8,17)-r3(8,21)-r3(8,25)*two
      rxyz(11)=+r6(19,11)-r5(13,15)-r5(13,19)*two+r4(8,20)*two+r4(8,24)+r4(8,28)*three &
&             -r3(4,17)-r3(4,21)-r3(4,25)*two
      rxyz(12)=+r7(33,3)-r6(25,5)-r6(25,7)*two+r5(18,7)*two+r5(18,9)+r5(18,11)*three &
&             +r5(20,21)-r4(12,8)-r4(12,10)-r4(12,12)*two-r4(14,30)-r4(14,34)*two &
&             +r3(9,31)*two+r3(9,35)+r3(9,39)*three-r2(6,25)-r2(6,29)-r2(6,33)*two
      rxyz(13)=+r7(33,4)-r6(25,6)-r6(25,8)*two+r5(18,8)*two+r5(18,10)+r5(18,12)*three &
&             +r5(20,22)-r4(12,9)-r4(12,11)-r4(12,13)*two-r4(14,31)-r4(14,35)*two &
&             +r3(9,32)*two+r3(9,36)+r3(9,40)*three-r2(6,26)-r2(6,30)-r2(6,34)*two
      rxyz(14)=+r7(25,4)-r6(18,6)-r6(18,8)*two+r5(12,8)*two+r5(12,10)+r5(12,12)*three &
&             +r5(14,22)-r4(7,9)-r4(7,11)-r4(7,13)*two-r4(9,31)-r4(9,35)*two+r3(5,32)*two &
&             +r3(5,36)+r3(5,40)*three-r2(4,26)-r2(4,30)-r2(4,34)*two
      rxyz(15)=+r7(25,3)-r6(18,5)-r6(18,7)*two+r5(12,7)*two+r5(12,9)+r5(12,11)*three &
&             +r5(14,21)-r4(7,8)-r4(7,10)-r4(7,12)*two-r4(9,30)-r4(9,34)*two+r3(5,31)*two &
&             +r3(5,35)+r3(5,39)*three-r2(4,25)-r2(4,29)-r2(4,33)*two
      rxyz(16)=+r7(19,4)-r6(13,6)-r6(13,8)*two+r5(8,8)*two+r5(8,10)+r5(8,12)*three &
&             +r5(19,22)-r4(4,9)-r4(4,11)-r4(4,13)*two-r4(13,31)-r4(13,35)*two &
&             +r3(8,32)*two+r3(8,36)+r3(8,40)*three-r2(2,26)-r2(2,30)-r2(2,34)*two
      rxyz(17)=+r7(19,3)-r6(13,5)-r6(13,7)*two+r5(8,7)*two+r5(8,9)+r5(8,11)*three &
&             +r5(19,21)-r4(4,8)-r4(4,10)-r4(4,12)*two-r4(13,30)-r4(13,34)*two &
&             +r3(8,31)*two+r3(8,35)+r3(8,39)*three-r2(2,25)-r2(2,29)-r2(2,33)*two
      rxyz(18)=+r7(34,3)-r6(26,5)-r6(26,7)*two+r5(19,7)*two+r5(19,9)+r5(19,11)*three &
&             +r5(19,21)-r4(13,8)-r4(13,10)-r4(13,12)*two-r4(13,30)-r4(13,34)*two &
&             +r3(8,31)*two+r3(8,35)+r3(8,39)*three-r2(2,25)-r2(2,29)-r2(2,33)*two
      rxyz(19)=+r7(34,4)-r6(26,6)-r6(26,8)*two+r5(19,8)*two+r5(19,10)+r5(19,12)*three &
&             +r5(19,22)-r4(13,9)-r4(13,11)-r4(13,13)*two-r4(13,31)-r4(13,35)*two &
&             +r3(8,32)*two+r3(8,36)+r3(8,40)*three-r2(2,26)-r2(2,30)-r2(2,34)*two
      rxyz(20)=+r6(20,11)*four-r5(14,15)*four-r5(14,19)*eight+r4(9,20)*eight+r4(9,24)*four &
&             +r4(9,28)*p12-r3(5,17)*four-r3(5,21)*four-r3(5,25)*eight
      eri(1,1,6,3)=r400+(+r7(14,3)*two+r7(14,4)*two-r6(9,5)*two-r6(9,6)*two-r6(9,7)*four &
&                  -r6(9,8)*four+r5(5,7)*four+r5(5,8)*four+r5(5,9)*two+r5(5,10)*two &
&                  +r5(5,11)*six+r5(5,12)*six+r5(14,21)*six+r5(14,22)*six-r4(2,8)*two &
&                  -r4(2,9)*two-r4(2,10)*two-r4(2,11)*two-r4(2,12)*four-r4(2,13)*four &
&                  -r4(9,30)*six-r4(9,31)*six-r4(9,34)*p12-r4(9,35)*p12+r3(5,31)*p12 &
&                  +r3(5,32)*p12+r3(5,35)*six+r3(5,36)*six+r3(5,39)*p18+r3(5,40)*p18 &
&                  -r2(4,25)*six-r2(4,26)*six-r2(4,29)*six-r2(4,30)*six-r2(4,33)*p12 &
&                  -r2(4,34)*p12)*qx+(+r6(14,10)+r6(14,11)*four+r6(14,12)-r5(9,14) &
&                  -r5(9,15)*four-r5(9,16)-r5(9,18)*two-r5(9,19)*eight-r5(9,20)*two &
&                  +r4(5,19)*two+r4(5,20)*eight+r4(5,21)*two+r4(5,23)+r4(5,24)*four+r4(5,25) &
&                  +r4(5,27)*three+r4(5,28)*p12+r4(5,29)*three+r4(14,39)+r4(14,40)*four &
&                  +r4(14,41)-r3(2,16)-r3(2,17)*four-r3(2,18)-r3(2,20)-r3(2,21)*four-r3(2,22) &
&                  -r3(2,24)*two-r3(2,25)*eight-r3(2,26)*two-r3(9,44)-r3(9,45)*four-r3(9,46) &
&                  -r3(9,49)*two-r3(9,50)*eight-r3(9,51)*two+r2(6,43)*two+r2(6,44)*eight &
&                  +r2(6,45)*two+r2(6,48)+r2(6,49)*four+r2(6,50)+r2(6,53)*three+r2(6,54)*p12 &
&                  +r2(6,55)*three-r1(2,27)-r1(2,28)*four-r1(2,29)-r1(2,32)-r1(2,33)*four &
&                  -r1(2,34)-r1(2,37)*two-r1(2,38)*eight-r1(2,39)*two)*xx+(+r5(14,23)*two &
&                  +r5(14,24)*two-r4(9,32)*two-r4(9,33)*two-r4(9,36)*four-r4(9,37)*four &
&                  +r3(5,33)*four+r3(5,34)*four+r3(5,37)*two+r3(5,38)*two+r3(5,41)*six &
&                  +r3(5,42)*six-r2(4,27)*two-r2(4,28)*two-r2(4,31)*two-r2(4,32)*two &
&                  -r2(4,35)*four-r2(4,36)*four)*xxx+rxyz(1)*xxxx
      eri(2,1,6,3)=r220+(+r7(25,4)*two-r6(18,6)*two-r6(18,8)*four+r5(12,8)*four &
&                  +r5(12,10)*two+r5(12,12)*six+r5(14,22)*two-r4(7,9)*two-r4(7,11)*two &
&                  -r4(7,13)*four-r4(9,31)*two-r4(9,35)*four+r3(5,32)*four+r3(5,36)*two &
&                  +r3(5,40)*six-r2(4,26)*two-r2(4,30)*two-r2(4,34)*four)*qx+rxyz(5)*xx
      eri(3,1,6,3)=r202+(+r7(27,4)*two-r6(20,6)*two-r6(20,8)*four+r5(14,8)*four &
&                  +r5(14,10)*two+r5(14,12)*six+r5(14,22)*two-r4(9,9)*two-r4(9,11)*two &
&                  -r4(9,13)*four-r4(9,31)*two-r4(9,35)*four+r3(5,32)*four+r3(5,36)*two &
&                  +r3(5,40)*six-r2(4,26)*two-r2(4,30)*two-r2(4,34)*four)*qx+(+r7(20,3)*two &
&                  -r6(14,5)*two-r6(14,7)*four+r5(9,7)*four+r5(9,9)*two+r5(9,11)*six &
&                  +r5(20,21)*two-r4(5,8)*two-r4(5,10)*two-r4(5,12)*four-r4(14,30)*two &
&                  -r4(14,34)*four+r3(9,31)*four+r3(9,35)*two+r3(9,39)*six-r2(6,25)*two &
&                  -r2(6,29)*two-r2(6,33)*four)*qz+(+r6(27,12)-r5(20,16)-r5(20,20)*two &
&                  +r4(14,21)*two+r4(14,25)+r4(14,29)*three+r4(14,41)-r3(9,18)-r3(9,22) &
&                  -r3(9,26)*two-r3(9,46)-r3(9,51)*two+r2(6,45)*two+r2(6,50)+r2(6,55)*three &
&                  -r1(2,29)-r1(2,34)-r1(2,39)*two)*xx+rxyz(20)*xz+(+r6(14,10)-r5(9,14) &
&                  -r5(9,18)*two+r4(5,19)*two+r4(5,23)+r4(5,27)*three+r4(14,39)-r3(2,16) &
&                  -r3(2,20)-r3(2,24)*two-r3(9,44)-r3(9,49)*two+r2(6,43)*two+r2(6,48) &
&                  +r2(6,53)*three-r1(2,27)-r1(2,32)-r1(2,37)*two)*zz+(+r5(20,24)*two &
&                  -r4(14,33)*two-r4(14,37)*four+r3(9,34)*four+r3(9,38)*two+r3(9,42)*six &
&                  -r2(6,28)*two-r2(6,32)*two-r2(6,36)*four)*xxz+(+r5(14,23)*two-r4(9,32)*two &
&                  -r4(9,36)*four+r3(5,33)*four+r3(5,37)*two+r3(5,41)*six-r2(4,27)*two &
&                  -r2(4,31)*two-r2(4,35)*four)*xzz+rxyz(1)*xxzz
      eri(4,1,6,3)=r310+(+r7(19,3)+r7(19,4)*two-r6(13,5)-r6(13,6)*two-r6(13,7)*two &
&                  -r6(13,8)*four+r5(8,7)*two+r5(8,8)*four+r5(8,9)+r5(8,10)*two &
&                  +r5(8,11)*three+r5(8,12)*six+r5(19,21)+r5(19,22)*two-r4(4,8)-r4(4,9)*two &
&                  -r4(4,10)-r4(4,11)*two-r4(4,12)*two-r4(4,13)*four-r4(13,30)-r4(13,31)*two &
&                  -r4(13,34)*two-r4(13,35)*four+r3(8,31)*two+r3(8,32)*four+r3(8,35) &
&                  +r3(8,36)*two+r3(8,39)*three+r3(8,40)*six-r2(2,25)-r2(2,26)*two-r2(2,29) &
&                  -r2(2,30)*two-r2(2,33)*two-r2(2,34)*four)*qx+(+r6(19,11)*two+r6(19,12) &
&                  -r5(13,15)*two-r5(13,16)-r5(13,19)*four-r5(13,20)*two+r4(8,20)*four &
&                  +r4(8,21)*two+r4(8,24)*two+r4(8,25)+r4(8,28)*six+r4(8,29)*three &
&                  -r3(4,17)*two-r3(4,18)-r3(4,21)*two-r3(4,22)-r3(4,25)*four-r3(4,26)*two) &
&                  *xx+rxyz(2)*xxx
      eri(5,1,6,3)=r301+(+r7(20,3)+r7(20,4)*two-r6(14,5)-r6(14,6)*two-r6(14,7)*two &
&                  -r6(14,8)*four+r5(9,7)*two+r5(9,8)*four+r5(9,9)+r5(9,10)*two &
&                  +r5(9,11)*three+r5(9,12)*six+r5(20,21)+r5(20,22)*two-r4(5,8)-r4(5,9)*two &
&                  -r4(5,10)-r4(5,11)*two-r4(5,12)*two-r4(5,13)*four-r4(14,30)-r4(14,31)*two &
&                  -r4(14,34)*two-r4(14,35)*four+r3(9,31)*two+r3(9,32)*four+r3(9,35) &
&                  +r3(9,36)*two+r3(9,39)*three+r3(9,40)*six-r2(6,25)-r2(6,26)*two-r2(6,29) &
&                  -r2(6,30)*two-r2(6,33)*two-r2(6,34)*four)*qx+(+r7(14,3)-r6(9,5) &
&                  -r6(9,7)*two+r5(5,7)*two+r5(5,9)+r5(5,11)*three+r5(14,21)*three-r4(2,8) &
&                  -r4(2,10)-r4(2,12)*two-r4(9,30)*three-r4(9,34)*six+r3(5,31)*six &
&                  +r3(5,35)*three+r3(5,39)*nine-r2(4,25)*three-r2(4,29)*three-r2(4,33)*six) &
&                  *qz+(+r6(20,11)*two+r6(20,12)-r5(14,15)*two-r5(14,16)-r5(14,19)*four &
&                  -r5(14,20)*two+r4(9,20)*four+r4(9,21)*two+r4(9,24)*two+r4(9,25) &
&                  +r4(9,28)*six+r4(9,29)*three-r3(5,17)*two-r3(5,18)-r3(5,21)*two-r3(5,22) &
&                  -r3(5,25)*four-r3(5,26)*two)*xx+(+r6(14,10)+r6(14,11)*two-r5(9,14) &
&                  -r5(9,15)*two-r5(9,18)*two-r5(9,19)*four+r4(5,19)*two+r4(5,20)*four &
&                  +r4(5,23)+r4(5,24)*two+r4(5,27)*three+r4(5,28)*six+r4(14,39)+r4(14,40)*two &
&                  -r3(2,16)-r3(2,17)*two-r3(2,20)-r3(2,21)*two-r3(2,24)*two-r3(2,25)*four &
&                  -r3(9,44)-r3(9,45)*two-r3(9,49)*two-r3(9,50)*four+r2(6,43)*two &
&                  +r2(6,44)*four+r2(6,48)+r2(6,49)*two+r2(6,53)*three+r2(6,54)*six-r1(2,27) &
&                  -r1(2,28)*two-r1(2,32)-r1(2,33)*two-r1(2,37)*two-r1(2,38)*four)*xz+( &
&                  +r5(20,24)-r4(14,33)-r4(14,37)*two+r3(9,34)*two+r3(9,38)+r3(9,42)*three &
&                  -r2(6,28)-r2(6,32)-r2(6,36)*two)*xxx+(+r5(14,23)*two+r5(14,24) &
&                  -r4(9,32)*two-r4(9,33)-r4(9,36)*four-r4(9,37)*two+r3(5,33)*four &
&                  +r3(5,34)*two+r3(5,37)*two+r3(5,38)+r3(5,41)*six+r3(5,42)*three &
&                  -r2(4,27)*two-r2(4,28)-r2(4,31)*two-r2(4,32)-r2(4,35)*four-r2(4,36)*two) &
&                  *xxz+rxyz(1)*xxxz
      eri(6,1,6,3)=r211+(+r7(26,4)*two-r6(19,6)*two-r6(19,8)*four+r5(13,8)*four &
&                  +r5(13,10)*two+r5(13,12)*six-r4(8,9)*two-r4(8,11)*two-r4(8,13)*four)*qx &
&                  +rxyz(17)*qz+(+r6(26,12)-r5(19,16)-r5(19,20)*two+r4(13,21)*two+r4(13,25) &
&                  +r4(13,29)*three-r3(8,18)-r3(8,22)-r3(8,26)*two)*xx+(+r6(19,11)*two &
&                  -r5(13,15)*two-r5(13,19)*four+r4(8,20)*four+r4(8,24)*two+r4(8,28)*six &
&                  -r3(4,17)*two-r3(4,21)*two-r3(4,25)*four)*xz+rxyz(2)*xxz
      eri(1,2,6,3)=r220+(+r7(25,3)*two-r6(18,5)*two-r6(18,7)*four+r5(12,7)*four &
&                  +r5(12,9)*two+r5(12,11)*six+r5(14,21)*two-r4(7,8)*two-r4(7,10)*two &
&                  -r4(7,12)*four-r4(9,30)*two-r4(9,34)*four+r3(5,31)*four+r3(5,35)*two &
&                  +r3(5,39)*six-r2(4,25)*two-r2(4,29)*two-r2(4,33)*four)*qx+rxyz(4)*xx
      eri(2,2,6,3)=r040
      eri(3,2,6,3)=r022+(+r7(33,3)*two-r6(25,5)*two-r6(25,7)*four+r5(18,7)*four &
&                  +r5(18,9)*two+r5(18,11)*six+r5(20,21)*two-r4(12,8)*two-r4(12,10)*two &
&                  -r4(12,12)*four-r4(14,30)*two-r4(14,34)*four+r3(9,31)*four+r3(9,35)*two &
&                  +r3(9,39)*six-r2(6,25)*two-r2(6,29)*two-r2(6,33)*four)*qz+rxyz(4)*zz
      eri(4,2,6,3)=r130+rxyz(7)*qx
      eri(5,2,6,3)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,6,3)=r031+rxyz(7)*qz
      eri(1,3,6,3)=r202+(+r7(27,3)*two-r6(20,5)*two-r6(20,7)*four+r5(14,7)*four &
&                  +r5(14,9)*two+r5(14,11)*six+r5(14,21)*two-r4(9,8)*two-r4(9,10)*two &
&                  -r4(9,12)*four-r4(9,30)*two-r4(9,34)*four+r3(5,31)*four+r3(5,35)*two &
&                  +r3(5,39)*six-r2(4,25)*two-r2(4,29)*two-r2(4,33)*four)*qx+(+r7(20,4)*two &
&                  -r6(14,6)*two-r6(14,8)*four+r5(9,8)*four+r5(9,10)*two+r5(9,12)*six &
&                  +r5(20,22)*two-r4(5,9)*two-r4(5,11)*two-r4(5,13)*four-r4(14,31)*two &
&                  -r4(14,35)*four+r3(9,32)*four+r3(9,36)*two+r3(9,40)*six-r2(6,26)*two &
&                  -r2(6,30)*two-r2(6,34)*four)*qz+(+r6(27,10)-r5(20,14)-r5(20,18)*two &
&                  +r4(14,19)*two+r4(14,23)+r4(14,27)*three+r4(14,39)-r3(9,16)-r3(9,20) &
&                  -r3(9,24)*two-r3(9,44)-r3(9,49)*two+r2(6,43)*two+r2(6,48)+r2(6,53)*three &
&                  -r1(2,27)-r1(2,32)-r1(2,37)*two)*xx+rxyz(20)*xz+(+r6(14,12)-r5(9,16) &
&                  -r5(9,20)*two+r4(5,21)*two+r4(5,25)+r4(5,29)*three+r4(14,41)-r3(2,18) &
&                  -r3(2,22)-r3(2,26)*two-r3(9,46)-r3(9,51)*two+r2(6,45)*two+r2(6,50) &
&                  +r2(6,55)*three-r1(2,29)-r1(2,34)-r1(2,39)*two)*zz+(+r5(20,23)*two &
&                  -r4(14,32)*two-r4(14,36)*four+r3(9,33)*four+r3(9,37)*two+r3(9,41)*six &
&                  -r2(6,27)*two-r2(6,31)*two-r2(6,35)*four)*xxz+(+r5(14,24)*two-r4(9,33)*two &
&                  -r4(9,37)*four+r3(5,34)*four+r3(5,38)*two+r3(5,42)*six-r2(4,28)*two &
&                  -r2(4,32)*two-r2(4,36)*four)*xzz+rxyz(1)*xxzz
      eri(2,3,6,3)=r022+(+r7(33,4)*two-r6(25,6)*two-r6(25,8)*four+r5(18,8)*four &
&                  +r5(18,10)*two+r5(18,12)*six+r5(20,22)*two-r4(12,9)*two-r4(12,11)*two &
&                  -r4(12,13)*four-r4(14,31)*two-r4(14,35)*four+r3(9,32)*four+r3(9,36)*two &
&                  +r3(9,40)*six-r2(6,26)*two-r2(6,30)*two-r2(6,34)*four)*qz+rxyz(5)*zz
      eri(3,3,6,3)=r004+(+r7(35,3)*two+r7(35,4)*two-r6(27,5)*two-r6(27,6)*two &
&                  -r6(27,7)*four-r6(27,8)*four+r5(20,7)*four+r5(20,8)*four+r5(20,9)*two &
&                  +r5(20,10)*two+r5(20,11)*six+r5(20,12)*six+r5(20,21)*six+r5(20,22)*six &
&                  -r4(14,8)*two-r4(14,9)*two-r4(14,10)*two-r4(14,11)*two-r4(14,12)*four &
&                  -r4(14,13)*four-r4(14,30)*six-r4(14,31)*six-r4(14,34)*p12-r4(14,35)*p12 &
&                  +r3(9,31)*p12+r3(9,32)*p12+r3(9,35)*six+r3(9,36)*six+r3(9,39)*p18 &
&                  +r3(9,40)*p18-r2(6,25)*six-r2(6,26)*six-r2(6,29)*six-r2(6,30)*six &
&                  -r2(6,33)*p12-r2(6,34)*p12)*qz+(+r6(27,10)+r6(27,11)*four+r6(27,12) &
&                  -r5(20,14)-r5(20,15)*four-r5(20,16)-r5(20,18)*two-r5(20,19)*eight &
&                  -r5(20,20)*two+r4(14,19)*two+r4(14,20)*eight+r4(14,21)*two+r4(14,23) &
&                  +r4(14,24)*four+r4(14,25)+r4(14,27)*three+r4(14,28)*p12+r4(14,29)*three &
&                  +r4(14,39)+r4(14,40)*four+r4(14,41)-r3(9,16)-r3(9,17)*four-r3(9,18) &
&                  -r3(9,20)-r3(9,21)*four-r3(9,22)-r3(9,24)*two-r3(9,25)*eight-r3(9,26)*two &
&                  -r3(9,44)-r3(9,45)*four-r3(9,46)-r3(9,49)*two-r3(9,50)*eight-r3(9,51)*two &
&                  +r2(6,43)*two+r2(6,44)*eight+r2(6,45)*two+r2(6,48)+r2(6,49)*four+r2(6,50) &
&                  +r2(6,53)*three+r2(6,54)*p12+r2(6,55)*three-r1(2,27)-r1(2,28)*four &
&                  -r1(2,29)-r1(2,32)-r1(2,33)*four-r1(2,34)-r1(2,37)*two-r1(2,38)*eight &
&                  -r1(2,39)*two)*zz+(+r5(20,23)*two+r5(20,24)*two-r4(14,32)*two &
&                  -r4(14,33)*two-r4(14,36)*four-r4(14,37)*four+r3(9,33)*four+r3(9,34)*four &
&                  +r3(9,37)*two+r3(9,38)*two+r3(9,41)*six+r3(9,42)*six-r2(6,27)*two &
&                  -r2(6,28)*two-r2(6,31)*two-r2(6,32)*two-r2(6,35)*four-r2(6,36)*four)*zzz &
&                  +rxyz(1)*zzzz
      eri(4,3,6,3)=r112+rxyz(18)*qx+(+r7(26,4)*two-r6(19,6)*two-r6(19,8)*four &
&                  +r5(13,8)*four+r5(13,10)*two+r5(13,12)*six-r4(8,9)*two-r4(8,11)*two &
&                  -r4(8,13)*four)*qz+(+r6(26,11)*two-r5(19,15)*two-r5(19,19)*four &
&                  +r4(13,20)*four+r4(13,24)*two+r4(13,28)*six-r3(8,17)*two-r3(8,21)*two &
&                  -r3(8,25)*four)*xz+(+r6(19,12)-r5(13,16)-r5(13,20)*two+r4(8,21)*two &
&                  +r4(8,25)+r4(8,29)*three-r3(4,18)-r3(4,22)-r3(4,26)*two)*zz+rxyz(2)*xzz
      eri(5,3,6,3)=r103+(+r7(35,3)-r6(27,5)-r6(27,7)*two+r5(20,7)*two+r5(20,9) &
&                  +r5(20,11)*three+r5(20,21)*three-r4(14,8)-r4(14,10)-r4(14,12)*two &
&                  -r4(14,30)*three-r4(14,34)*six+r3(9,31)*six+r3(9,35)*three+r3(9,39)*nine &
&                  -r2(6,25)*three-r2(6,29)*three-r2(6,33)*six)*qx+(+r7(27,3)+r7(27,4)*two &
&                  -r6(20,5)-r6(20,6)*two-r6(20,7)*two-r6(20,8)*four+r5(14,7)*two &
&                  +r5(14,8)*four+r5(14,9)+r5(14,10)*two+r5(14,11)*three+r5(14,12)*six &
&                  +r5(14,21)+r5(14,22)*two-r4(9,8)-r4(9,9)*two-r4(9,10)-r4(9,11)*two &
&                  -r4(9,12)*two-r4(9,13)*four-r4(9,30)-r4(9,31)*two-r4(9,34)*two &
&                  -r4(9,35)*four+r3(5,31)*two+r3(5,32)*four+r3(5,35)+r3(5,36)*two &
&                  +r3(5,39)*three+r3(5,40)*six-r2(4,25)-r2(4,26)*two-r2(4,29)-r2(4,30)*two &
&                  -r2(4,33)*two-r2(4,34)*four)*qz+(+r6(27,10)+r6(27,11)*two-r5(20,14) &
&                  -r5(20,15)*two-r5(20,18)*two-r5(20,19)*four+r4(14,19)*two+r4(14,20)*four &
&                  +r4(14,23)+r4(14,24)*two+r4(14,27)*three+r4(14,28)*six+r4(14,39) &
&                  +r4(14,40)*two-r3(9,16)-r3(9,17)*two-r3(9,20)-r3(9,21)*two-r3(9,24)*two &
&                  -r3(9,25)*four-r3(9,44)-r3(9,45)*two-r3(9,49)*two-r3(9,50)*four &
&                  +r2(6,43)*two+r2(6,44)*four+r2(6,48)+r2(6,49)*two+r2(6,53)*three &
&                  +r2(6,54)*six-r1(2,27)-r1(2,28)*two-r1(2,32)-r1(2,33)*two-r1(2,37)*two &
&                  -r1(2,38)*four)*xz+(+r6(20,11)*two+r6(20,12)-r5(14,15)*two-r5(14,16) &
&                  -r5(14,19)*four-r5(14,20)*two+r4(9,20)*four+r4(9,21)*two+r4(9,24)*two &
&                  +r4(9,25)+r4(9,28)*six+r4(9,29)*three-r3(5,17)*two-r3(5,18)-r3(5,21)*two &
&                  -r3(5,22)-r3(5,25)*four-r3(5,26)*two)*zz+(+r5(20,23)*two+r5(20,24) &
&                  -r4(14,32)*two-r4(14,33)-r4(14,36)*four-r4(14,37)*two+r3(9,33)*four &
&                  +r3(9,34)*two+r3(9,37)*two+r3(9,38)+r3(9,41)*six+r3(9,42)*three &
&                  -r2(6,27)*two-r2(6,28)-r2(6,31)*two-r2(6,32)-r2(6,35)*four-r2(6,36)*two) &
&                  *xzz+(+r5(14,24)-r4(9,33)-r4(9,37)*two+r3(5,34)*two+r3(5,38) &
&                  +r3(5,42)*three-r2(4,28)-r2(4,32)-r2(4,36)*two)*zzz+rxyz(1)*xzzz
      eri(6,3,6,3)=r013+(+r7(34,3)+r7(34,4)*two-r6(26,5)-r6(26,6)*two-r6(26,7)*two &
&                  -r6(26,8)*four+r5(19,7)*two+r5(19,8)*four+r5(19,9)+r5(19,10)*two &
&                  +r5(19,11)*three+r5(19,12)*six+r5(19,21)+r5(19,22)*two-r4(13,8) &
&                  -r4(13,9)*two-r4(13,10)-r4(13,11)*two-r4(13,12)*two-r4(13,13)*four &
&                  -r4(13,30)-r4(13,31)*two-r4(13,34)*two-r4(13,35)*four+r3(8,31)*two &
&                  +r3(8,32)*four+r3(8,35)+r3(8,36)*two+r3(8,39)*three+r3(8,40)*six-r2(2,25) &
&                  -r2(2,26)*two-r2(2,29)-r2(2,30)*two-r2(2,33)*two-r2(2,34)*four)*qz+( &
&                  +r6(26,11)*two+r6(26,12)-r5(19,15)*two-r5(19,16)-r5(19,19)*four &
&                  -r5(19,20)*two+r4(13,20)*four+r4(13,21)*two+r4(13,24)*two+r4(13,25) &
&                  +r4(13,28)*six+r4(13,29)*three-r3(8,17)*two-r3(8,18)-r3(8,21)*two-r3(8,22) &
&                  -r3(8,25)*four-r3(8,26)*two)*zz+rxyz(2)*zzz
      eri(1,4,6,3)=r310+(+r7(19,3)*two+r7(19,4)-r6(13,5)*two-r6(13,6)-r6(13,7)*four &
&                  -r6(13,8)*two+r5(8,7)*four+r5(8,8)*two+r5(8,9)*two+r5(8,10)+r5(8,11)*six &
&                  +r5(8,12)*three+r5(19,21)*two+r5(19,22)-r4(4,8)*two-r4(4,9)-r4(4,10)*two &
&                  -r4(4,11)-r4(4,12)*four-r4(4,13)*two-r4(13,30)*two-r4(13,31) &
&                  -r4(13,34)*four-r4(13,35)*two+r3(8,31)*four+r3(8,32)*two+r3(8,35)*two &
&                  +r3(8,36)+r3(8,39)*six+r3(8,40)*three-r2(2,25)*two-r2(2,26)-r2(2,29)*two &
&                  -r2(2,30)-r2(2,33)*four-r2(2,34)*two)*qx+(+r6(19,10)+r6(19,11)*two &
&                  -r5(13,14)-r5(13,15)*two-r5(13,18)*two-r5(13,19)*four+r4(8,19)*two &
&                  +r4(8,20)*four+r4(8,23)+r4(8,24)*two+r4(8,27)*three+r4(8,28)*six-r3(4,16) &
&                  -r3(4,17)*two-r3(4,20)-r3(4,21)*two-r3(4,24)*two-r3(4,25)*four)*xx+rxyz(3) &
&                  *xxx
      eri(2,4,6,3)=r130+rxyz(8)*qx
      eri(3,4,6,3)=r112+rxyz(19)*qx+(+r7(26,3)*two-r6(19,5)*two-r6(19,7)*four &
&                  +r5(13,7)*four+r5(13,9)*two+r5(13,11)*six-r4(8,8)*two-r4(8,10)*two &
&                  -r4(8,12)*four)*qz+(+r6(26,11)*two-r5(19,15)*two-r5(19,19)*four &
&                  +r4(13,20)*four+r4(13,24)*two+r4(13,28)*six-r3(8,17)*two-r3(8,21)*two &
&                  -r3(8,25)*four)*xz+(+r6(19,10)-r5(13,14)-r5(13,18)*two+r4(8,19)*two &
&                  +r4(8,23)+r4(8,27)*three-r3(4,16)-r3(4,20)-r3(4,24)*two)*zz+rxyz(3)*xzz
      eri(4,4,6,3)=r220+(+r7(25,3)+r7(25,4)-r6(18,5)-r6(18,6)-r6(18,7)*two-r6(18,8)*two &
&                  +r5(12,7)*two+r5(12,8)*two+r5(12,9)+r5(12,10)+r5(12,11)*three &
&                  +r5(12,12)*three+r5(14,21)+r5(14,22)-r4(7,8)-r4(7,9)-r4(7,10)-r4(7,11) &
&                  -r4(7,12)*two-r4(7,13)*two-r4(9,30)-r4(9,31)-r4(9,34)*two-r4(9,35)*two &
&                  +r3(5,31)*two+r3(5,32)*two+r3(5,35)+r3(5,36)+r3(5,39)*three+r3(5,40)*three &
&                  -r2(4,25)-r2(4,26)-r2(4,29)-r2(4,30)-r2(4,33)*two-r2(4,34)*two)*qx+rxyz(6) &
&                  *xx
      eri(5,4,6,3)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(19,10)+r6(19,11) &
&                  -r5(13,14)-r5(13,15)-r5(13,18)*two-r5(13,19)*two+r4(8,19)*two+r4(8,20)*two &
&                  +r4(8,23)+r4(8,24)+r4(8,27)*three+r4(8,28)*three-r3(4,16)-r3(4,17) &
&                  -r3(4,20)-r3(4,21)-r3(4,24)*two-r3(4,25)*two)*xz+rxyz(3)*xxz
      eri(6,4,6,3)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,6,3)=r301+(+r7(20,3)*two+r7(20,4)-r6(14,5)*two-r6(14,6)-r6(14,7)*four &
&                  -r6(14,8)*two+r5(9,7)*four+r5(9,8)*two+r5(9,9)*two+r5(9,10)+r5(9,11)*six &
&                  +r5(9,12)*three+r5(20,21)*two+r5(20,22)-r4(5,8)*two-r4(5,9)-r4(5,10)*two &
&                  -r4(5,11)-r4(5,12)*four-r4(5,13)*two-r4(14,30)*two-r4(14,31) &
&                  -r4(14,34)*four-r4(14,35)*two+r3(9,31)*four+r3(9,32)*two+r3(9,35)*two &
&                  +r3(9,36)+r3(9,39)*six+r3(9,40)*three-r2(6,25)*two-r2(6,26)-r2(6,29)*two &
&                  -r2(6,30)-r2(6,33)*four-r2(6,34)*two)*qx+(+r7(14,4)-r6(9,6)-r6(9,8)*two &
&                  +r5(5,8)*two+r5(5,10)+r5(5,12)*three+r5(14,22)*three-r4(2,9)-r4(2,11) &
&                  -r4(2,13)*two-r4(9,31)*three-r4(9,35)*six+r3(5,32)*six+r3(5,36)*three &
&                  +r3(5,40)*nine-r2(4,26)*three-r2(4,30)*three-r2(4,34)*six)*qz+(+r6(20,10) &
&                  +r6(20,11)*two-r5(14,14)-r5(14,15)*two-r5(14,18)*two-r5(14,19)*four &
&                  +r4(9,19)*two+r4(9,20)*four+r4(9,23)+r4(9,24)*two+r4(9,27)*three &
&                  +r4(9,28)*six-r3(5,16)-r3(5,17)*two-r3(5,20)-r3(5,21)*two-r3(5,24)*two &
&                  -r3(5,25)*four)*xx+(+r6(14,11)*two+r6(14,12)-r5(9,15)*two-r5(9,16) &
&                  -r5(9,19)*four-r5(9,20)*two+r4(5,20)*four+r4(5,21)*two+r4(5,24)*two &
&                  +r4(5,25)+r4(5,28)*six+r4(5,29)*three+r4(14,40)*two+r4(14,41)-r3(2,17)*two &
&                  -r3(2,18)-r3(2,21)*two-r3(2,22)-r3(2,25)*four-r3(2,26)*two-r3(9,45)*two &
&                  -r3(9,46)-r3(9,50)*four-r3(9,51)*two+r2(6,44)*four+r2(6,45)*two &
&                  +r2(6,49)*two+r2(6,50)+r2(6,54)*six+r2(6,55)*three-r1(2,28)*two-r1(2,29) &
&                  -r1(2,33)*two-r1(2,34)-r1(2,38)*four-r1(2,39)*two)*xz+(+r5(20,23) &
&                  -r4(14,32)-r4(14,36)*two+r3(9,33)*two+r3(9,37)+r3(9,41)*three-r2(6,27) &
&                  -r2(6,31)-r2(6,35)*two)*xxx+(+r5(14,23)+r5(14,24)*two-r4(9,32) &
&                  -r4(9,33)*two-r4(9,36)*two-r4(9,37)*four+r3(5,33)*two+r3(5,34)*four &
&                  +r3(5,37)+r3(5,38)*two+r3(5,41)*three+r3(5,42)*six-r2(4,27)-r2(4,28)*two &
&                  -r2(4,31)-r2(4,32)*two-r2(4,35)*two-r2(4,36)*four)*xxz+rxyz(1)*xxxz
      eri(2,5,6,3)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,6,3)=r103+(+r7(35,4)-r6(27,6)-r6(27,8)*two+r5(20,8)*two+r5(20,10) &
&                  +r5(20,12)*three+r5(20,22)*three-r4(14,9)-r4(14,11)-r4(14,13)*two &
&                  -r4(14,31)*three-r4(14,35)*six+r3(9,32)*six+r3(9,36)*three+r3(9,40)*nine &
&                  -r2(6,26)*three-r2(6,30)*three-r2(6,34)*six)*qx+(+r7(27,3)*two+r7(27,4) &
&                  -r6(20,5)*two-r6(20,6)-r6(20,7)*four-r6(20,8)*two+r5(14,7)*four &
&                  +r5(14,8)*two+r5(14,9)*two+r5(14,10)+r5(14,11)*six+r5(14,12)*three &
&                  +r5(14,21)*two+r5(14,22)-r4(9,8)*two-r4(9,9)-r4(9,10)*two-r4(9,11) &
&                  -r4(9,12)*four-r4(9,13)*two-r4(9,30)*two-r4(9,31)-r4(9,34)*four &
&                  -r4(9,35)*two+r3(5,31)*four+r3(5,32)*two+r3(5,35)*two+r3(5,36) &
&                  +r3(5,39)*six+r3(5,40)*three-r2(4,25)*two-r2(4,26)-r2(4,29)*two-r2(4,30) &
&                  -r2(4,33)*four-r2(4,34)*two)*qz+(+r6(27,11)*two+r6(27,12)-r5(20,15)*two &
&                  -r5(20,16)-r5(20,19)*four-r5(20,20)*two+r4(14,20)*four+r4(14,21)*two &
&                  +r4(14,24)*two+r4(14,25)+r4(14,28)*six+r4(14,29)*three+r4(14,40)*two &
&                  +r4(14,41)-r3(9,17)*two-r3(9,18)-r3(9,21)*two-r3(9,22)-r3(9,25)*four &
&                  -r3(9,26)*two-r3(9,45)*two-r3(9,46)-r3(9,50)*four-r3(9,51)*two &
&                  +r2(6,44)*four+r2(6,45)*two+r2(6,49)*two+r2(6,50)+r2(6,54)*six &
&                  +r2(6,55)*three-r1(2,28)*two-r1(2,29)-r1(2,33)*two-r1(2,34)-r1(2,38)*four &
&                  -r1(2,39)*two)*xz+(+r6(20,10)+r6(20,11)*two-r5(14,14)-r5(14,15)*two &
&                  -r5(14,18)*two-r5(14,19)*four+r4(9,19)*two+r4(9,20)*four+r4(9,23) &
&                  +r4(9,24)*two+r4(9,27)*three+r4(9,28)*six-r3(5,16)-r3(5,17)*two-r3(5,20) &
&                  -r3(5,21)*two-r3(5,24)*two-r3(5,25)*four)*zz+(+r5(20,23)+r5(20,24)*two &
&                  -r4(14,32)-r4(14,33)*two-r4(14,36)*two-r4(14,37)*four+r3(9,33)*two &
&                  +r3(9,34)*four+r3(9,37)+r3(9,38)*two+r3(9,41)*three+r3(9,42)*six-r2(6,27) &
&                  -r2(6,28)*two-r2(6,31)-r2(6,32)*two-r2(6,35)*two-r2(6,36)*four)*xzz+( &
&                  +r5(14,23)-r4(9,32)-r4(9,36)*two+r3(5,33)*two+r3(5,37)+r3(5,41)*three &
&                  -r2(4,27)-r2(4,31)-r2(4,35)*two)*zzz+rxyz(1)*xzzz
      eri(4,5,6,3)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(19,11)+r6(19,12) &
&                  -r5(13,15)-r5(13,16)-r5(13,19)*two-r5(13,20)*two+r4(8,20)*two+r4(8,21)*two &
&                  +r4(8,24)+r4(8,25)+r4(8,28)*three+r4(8,29)*three-r3(4,17)-r3(4,18) &
&                  -r3(4,21)-r3(4,22)-r3(4,25)*two-r3(4,26)*two)*xz+rxyz(2)*xxz
      eri(5,5,6,3)=r202+(+r7(27,3)+r7(27,4)-r6(20,5)-r6(20,6)-r6(20,7)*two-r6(20,8)*two &
&                  +r5(14,7)*two+r5(14,8)*two+r5(14,9)+r5(14,10)+r5(14,11)*three &
&                  +r5(14,12)*three+r5(14,21)+r5(14,22)-r4(9,8)-r4(9,9)-r4(9,10)-r4(9,11) &
&                  -r4(9,12)*two-r4(9,13)*two-r4(9,30)-r4(9,31)-r4(9,34)*two-r4(9,35)*two &
&                  +r3(5,31)*two+r3(5,32)*two+r3(5,35)+r3(5,36)+r3(5,39)*three+r3(5,40)*three &
&                  -r2(4,25)-r2(4,26)-r2(4,29)-r2(4,30)-r2(4,33)*two-r2(4,34)*two)*qx+( &
&                  +r7(20,3)+r7(20,4)-r6(14,5)-r6(14,6)-r6(14,7)*two-r6(14,8)*two+r5(9,7)*two &
&                  +r5(9,8)*two+r5(9,9)+r5(9,10)+r5(9,11)*three+r5(9,12)*three+r5(20,21) &
&                  +r5(20,22)-r4(5,8)-r4(5,9)-r4(5,10)-r4(5,11)-r4(5,12)*two-r4(5,13)*two &
&                  -r4(14,30)-r4(14,31)-r4(14,34)*two-r4(14,35)*two+r3(9,31)*two+r3(9,32)*two &
&                  +r3(9,35)+r3(9,36)+r3(9,39)*three+r3(9,40)*three-r2(6,25)-r2(6,26) &
&                  -r2(6,29)-r2(6,30)-r2(6,33)*two-r2(6,34)*two)*qz+(+r6(27,11)-r5(20,15) &
&                  -r5(20,19)*two+r4(14,20)*two+r4(14,24)+r4(14,28)*three+r4(14,40)-r3(9,17) &
&                  -r3(9,21)-r3(9,25)*two-r3(9,45)-r3(9,50)*two+r2(6,44)*two+r2(6,49) &
&                  +r2(6,54)*three-r1(2,28)-r1(2,33)-r1(2,38)*two)*xx+(+r6(20,10) &
&                  +r6(20,11)*two+r6(20,12)-r5(14,14)-r5(14,15)*two-r5(14,16)-r5(14,18)*two &
&                  -r5(14,19)*four-r5(14,20)*two+r4(9,19)*two+r4(9,20)*four+r4(9,21)*two &
&                  +r4(9,23)+r4(9,24)*two+r4(9,25)+r4(9,27)*three+r4(9,28)*six+r4(9,29)*three &
&                  -r3(5,16)-r3(5,17)*two-r3(5,18)-r3(5,20)-r3(5,21)*two-r3(5,22) &
&                  -r3(5,24)*two-r3(5,25)*four-r3(5,26)*two)*xz+(+r6(14,11)-r5(9,15) &
&                  -r5(9,19)*two+r4(5,20)*two+r4(5,24)+r4(5,28)*three+r4(14,40)-r3(2,17) &
&                  -r3(2,21)-r3(2,25)*two-r3(9,45)-r3(9,50)*two+r2(6,44)*two+r2(6,49) &
&                  +r2(6,54)*three-r1(2,28)-r1(2,33)-r1(2,38)*two)*zz+(+r5(20,23)+r5(20,24) &
&                  -r4(14,32)-r4(14,33)-r4(14,36)*two-r4(14,37)*two+r3(9,33)*two+r3(9,34)*two &
&                  +r3(9,37)+r3(9,38)+r3(9,41)*three+r3(9,42)*three-r2(6,27)-r2(6,28) &
&                  -r2(6,31)-r2(6,32)-r2(6,35)*two-r2(6,36)*two)*xxz+(+r5(14,23)+r5(14,24) &
&                  -r4(9,32)-r4(9,33)-r4(9,36)*two-r4(9,37)*two+r3(5,33)*two+r3(5,34)*two &
&                  +r3(5,37)+r3(5,38)+r3(5,41)*three+r3(5,42)*three-r2(4,27)-r2(4,28) &
&                  -r2(4,31)-r2(4,32)-r2(4,35)*two-r2(4,36)*two)*xzz+rxyz(1)*xxzz
      eri(6,5,6,3)=r112+rxyz(19)*qx+(+r7(26,3)+r7(26,4)-r6(19,5)-r6(19,6)-r6(19,7)*two &
&                  -r6(19,8)*two+r5(13,7)*two+r5(13,8)*two+r5(13,9)+r5(13,10)+r5(13,11)*three &
&                  +r5(13,12)*three-r4(8,8)-r4(8,9)-r4(8,10)-r4(8,11)-r4(8,12)*two &
&                  -r4(8,13)*two)*qz+(+r6(26,11)+r6(26,12)-r5(19,15)-r5(19,16)-r5(19,19)*two &
&                  -r5(19,20)*two+r4(13,20)*two+r4(13,21)*two+r4(13,24)+r4(13,25) &
&                  +r4(13,28)*three+r4(13,29)*three-r3(8,17)-r3(8,18)-r3(8,21)-r3(8,22) &
&                  -r3(8,25)*two-r3(8,26)*two)*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,6,3)=r211+(+r7(26,3)*two-r6(19,5)*two-r6(19,7)*four+r5(13,7)*four &
&                  +r5(13,9)*two+r5(13,11)*six-r4(8,8)*two-r4(8,10)*two-r4(8,12)*four)*qx &
&                  +rxyz(16)*qz+(+r6(26,10)-r5(19,14)-r5(19,18)*two+r4(13,19)*two+r4(13,23) &
&                  +r4(13,27)*three-r3(8,16)-r3(8,20)-r3(8,24)*two)*xx+(+r6(19,11)*two &
&                  -r5(13,15)*two-r5(13,19)*four+r4(8,20)*four+r4(8,24)*two+r4(8,28)*six &
&                  -r3(4,17)*two-r3(4,21)*two-r3(4,25)*four)*xz+rxyz(3)*xxz
      eri(2,6,6,3)=r031+rxyz(8)*qz
      eri(3,6,6,3)=r013+(+r7(34,3)*two+r7(34,4)-r6(26,5)*two-r6(26,6)-r6(26,7)*four &
&                  -r6(26,8)*two+r5(19,7)*four+r5(19,8)*two+r5(19,9)*two+r5(19,10) &
&                  +r5(19,11)*six+r5(19,12)*three+r5(19,21)*two+r5(19,22)-r4(13,8)*two &
&                  -r4(13,9)-r4(13,10)*two-r4(13,11)-r4(13,12)*four-r4(13,13)*two &
&                  -r4(13,30)*two-r4(13,31)-r4(13,34)*four-r4(13,35)*two+r3(8,31)*four &
&                  +r3(8,32)*two+r3(8,35)*two+r3(8,36)+r3(8,39)*six+r3(8,40)*three &
&                  -r2(2,25)*two-r2(2,26)-r2(2,29)*two-r2(2,30)-r2(2,33)*four-r2(2,34)*two) &
&                  *qz+(+r6(26,10)+r6(26,11)*two-r5(19,14)-r5(19,15)*two-r5(19,18)*two &
&                  -r5(19,19)*four+r4(13,19)*two+r4(13,20)*four+r4(13,23)+r4(13,24)*two &
&                  +r4(13,27)*three+r4(13,28)*six-r3(8,16)-r3(8,17)*two-r3(8,20)-r3(8,21)*two &
&                  -r3(8,24)*two-r3(8,25)*four)*zz+rxyz(3)*zzz
      eri(4,6,6,3)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,6,3)=r112+rxyz(18)*qx+(+r7(26,3)+r7(26,4)-r6(19,5)-r6(19,6)-r6(19,7)*two &
&                  -r6(19,8)*two+r5(13,7)*two+r5(13,8)*two+r5(13,9)+r5(13,10)+r5(13,11)*three &
&                  +r5(13,12)*three-r4(8,8)-r4(8,9)-r4(8,10)-r4(8,11)-r4(8,12)*two &
&                  -r4(8,13)*two)*qz+(+r6(26,10)+r6(26,11)-r5(19,14)-r5(19,15)-r5(19,18)*two &
&                  -r5(19,19)*two+r4(13,19)*two+r4(13,20)*two+r4(13,23)+r4(13,24) &
&                  +r4(13,27)*three+r4(13,28)*three-r3(8,16)-r3(8,17)-r3(8,20)-r3(8,21) &
&                  -r3(8,24)*two-r3(8,25)*two)*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,6,3)=r022+(+r7(33,3)+r7(33,4)-r6(25,5)-r6(25,6)-r6(25,7)*two-r6(25,8)*two &
&                  +r5(18,7)*two+r5(18,8)*two+r5(18,9)+r5(18,10)+r5(18,11)*three &
&                  +r5(18,12)*three+r5(20,21)+r5(20,22)-r4(12,8)-r4(12,9)-r4(12,10)-r4(12,11) &
&                  -r4(12,12)*two-r4(12,13)*two-r4(14,30)-r4(14,31)-r4(14,34)*two &
&                  -r4(14,35)*two+r3(9,31)*two+r3(9,32)*two+r3(9,35)+r3(9,36)+r3(9,39)*three &
&                  +r3(9,40)*three-r2(6,25)-r2(6,26)-r2(6,29)-r2(6,30)-r2(6,33)*two &
&                  -r2(6,34)*two)*qz+rxyz(6)*zz
      return
end

!-------------------------------------------------------------
  subroutine int2dddd4(eri,r0,r1,r2,r3,r4,r5,r6,r7,r8,qx,qz)
!-------------------------------------------------------------
!
      implicit none
      integer :: i,j
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00, four=4.0D+00, five=5.0D+00
      real(8),parameter :: six=6.0D+00, seven=7.0D+00, eight=8.0D+00, nine=9.0D+00, ten=1.0D+01
      real(8),parameter :: p11=1.1D+01, p12=1.2D+01, p13=1.3D+01, p15=1.5D+01, p16=1.6D+01
      real(8),parameter :: p18=1.8D+01, p20=2.0D+01, p21=2.1D+01, p24=2.4D+01, p28=2.8D+01
      real(8),parameter :: p30=3.0D+01, p36=3.6D+01, p45=4.5D+01, p105=1.05D+2, p210=2.1D+02
      real(8),parameter :: p420=4.2D+02
      real(8),intent(in) :: r0(25), r1(3,40), r2(6,56), r3(10,52), r4(15,42), r5(21,24)
      real(8),intent(in) :: r6(28,12), r7(36,4), r8(45), qx, qz
      real(8),intent(inout) :: eri(6,6,6,6)
      real(8) :: r400, r310, r301, r220, r211, r202, r130, r121, r112, r103, r040, r031
      real(8) :: r022, r013, r004, rxyz(20), xx, xz, zz, xxx, xxz, xzz, zzz
      real(8) :: xxxx, xxxz, xxzz, xzzz, zzzz
!
      xx= qx*qx
      xz= qx*qz
      zz= qz*qz
      xxx= qx*qx*qx
      xxz= qx*qx*qz
      xzz= qx*qz*qz
      zzz= qz*qz*qz
      xxxx= xx*xx
      xxxz= xx*xz
      xxzz= xx*zz
      xzzz= xz*zz
      zzzz= zz*zz
!
!
      do i= 1,6
        do j= 1,6
          eri(j,i,1,4)= eri(j,i,4,1)
        enddo
      enddo
      do i= 1,6
        do j= 1,6
          eri(j,i,2,4)= eri(j,i,4,2)
        enddo
      enddo
!
      r400= r8(9)-r7(5,1)*two+r6(2,1)+r6(2,4)+r6(9,9)*six-r5(5,13)*p12+r4(2,14)*six &
&          +r4(2,26)*six+r4(9,38)*three-r3(5,43)*six+r2(4,37)*three+r2(4,52)*three
      r310= r8(13)-r7(8,1)*two+r6(4,1)+r6(4,4)+r6(13,9)*three-r5(8,13)*six+r4(4,14)*three &
&          +r4(4,26)*three
      r301= r8(14)-r7(9,1)*two+r6(5,1)+r6(5,4)+r6(14,9)*three-r5(9,13)*six+r4(5,14)*three &
&          +r4(5,26)*three
      r220= r8(18)-r7(12,1)*two+r6(7,1)+r6(7,4)+r6(9,9)+r6(18,9)-r5(5,13)*two-r5(12,13)*two &
&          +r4(2,14)+r4(7,14)+r4(2,26)+r4(7,26)+r4(9,38)-r3(5,43)*two+r2(4,37)+r2(4,52)
      r211= r8(19)-r7(13,1)*two+r6(8,1)+r6(8,4)+r6(19,9)-r5(13,13)*two+r4(8,14)+r4(8,26)
      r202= r8(20)-r7(14,1)*two+r6(9,1)+r6(9,4)+r6(9,9)+r6(20,9)-r5(5,13)*two-r5(14,13)*two &
&          +r4(2,14)+r4(9,14)+r4(2,26)+r4(9,26)+r4(9,38)-r3(5,43)*two+r2(4,37)+r2(4,52)
      r130= r8(24)-r7(17,1)*two+r6(11,1)+r6(11,4)+r6(13,9)*three-r5(8,13)*six &
&          +r4(4,14)*three+r4(4,26)*three
      r121= r8(25)-r7(18,1)*two+r6(12,1)+r6(12,4)+r6(14,9)-r5(9,13)*two+r4(5,14)+r4(5,26)
      r112= r8(26)-r7(19,1)*two+r6(13,1)+r6(13,4)+r6(13,9)-r5(8,13)*two+r4(4,14)+r4(4,26)
      r103= r8(27)-r7(20,1)*two+r6(14,1)+r6(14,4)+r6(14,9)*three-r5(9,13)*six &
&          +r4(5,14)*three+r4(5,26)*three
      r040= r8(31)-r7(23,1)*two+r6(16,1)+r6(16,4)+r6(18,9)*six-r5(12,13)*p12+r4(7,14)*six &
&          +r4(7,26)*six+r4(9,38)*three-r3(5,43)*six+r2(4,37)*three+r2(4,52)*three
      r031= r8(32)-r7(24,1)*two+r6(17,1)+r6(17,4)+r6(19,9)*three-r5(13,13)*six &
&          +r4(8,14)*three+r4(8,26)*three
      r022= r8(33)-r7(25,1)*two+r6(18,1)+r6(18,4)+r6(18,9)+r6(20,9)-r5(12,13)*two &
&          -r5(14,13)*two+r4(7,14)+r4(9,14)+r4(7,26)+r4(9,26)+r4(9,38)-r3(5,43)*two &
&          +r2(4,37)+r2(4,52)
      r013= r8(34)-r7(26,1)*two+r6(19,1)+r6(19,4)+r6(19,9)*three-r5(13,13)*six &
&          +r4(8,14)*three+r4(8,26)*three
      r004= r8(35)-r7(27,1)*two+r6(20,1)+r6(20,4)+r6(20,9)*six-r5(14,13)*p12+r4(9,14)*six &
&          +r4(9,26)*six+r4(9,38)*three-r3(5,43)*six+r2(4,37)*three+r2(4,52)*three
      rxyz(1)=+r4(9,42)-r3(5,47)*two+r2(4,41)+r2(4,56)
      rxyz(2)=+r5(13,24)-r4(8,33)*two+r3(4,30)+r3(4,42)
      rxyz(3)=+r5(13,23)-r4(8,32)*two+r3(4,29)+r3(4,41)
      rxyz(4)=+r6(18,10)-r5(12,14)*two+r4(7,15)+r4(7,27)+r4(9,39)-r3(5,44)*two+r2(4,38) &
&             +r2(4,53)
      rxyz(5)=+r6(18,12)-r5(12,16)*two+r4(7,17)+r4(7,29)+r4(9,41)-r3(5,46)*two+r2(4,40) &
&             +r2(4,55)
      rxyz(6)=+r6(18,11)-r5(12,15)*two+r4(7,16)+r4(7,28)+r4(9,40)-r3(5,45)*two+r2(4,39) &
&             +r2(4,54)
      rxyz(7)=+r7(24,3)-r6(17,5)*two+r5(11,5)+r5(11,11)+r5(13,21)*three-r4(8,30)*six &
&             +r3(4,27)*three+r3(4,39)*three
      rxyz(8)=+r7(24,4)-r6(17,6)*two+r5(11,6)+r5(11,12)+r5(13,22)*three-r4(8,31)*six &
&             +r3(4,28)*three+r3(4,40)*three
      rxyz(9)=+r7(19,3)+r7(19,4)-r6(13,5)*two-r6(13,6)*two+r5(8,5)+r5(8,6)+r5(8,11) &
&             +r5(8,12)
      rxyz(10)=+r6(19,11)-r5(13,15)*two+r4(8,16)+r4(8,28)
      rxyz(11)=+r6(13,11)-r5(8,15)*two+r4(4,16)+r4(4,28)
      rxyz(12)=+r7(25,3)-r6(18,5)*two+r5(12,5)+r5(12,11)+r5(14,21)-r4(9,30)*two+r3(5,27) &
&             +r3(5,39)
      rxyz(13)=+r7(25,4)-r6(18,6)*two+r5(12,6)+r5(12,12)+r5(14,22)-r4(9,31)*two+r3(5,28) &
&             +r3(5,40)
      rxyz(14)=+r7(18,4)-r6(12,6)*two+r5(7,6)+r5(7,12)+r5(9,22)-r4(5,31)*two+r3(2,28) &
&             +r3(2,40)
      rxyz(15)=+r7(18,3)-r6(12,5)*two+r5(7,5)+r5(7,11)+r5(9,21)-r4(5,30)*two+r3(2,27) &
&             +r3(2,39)
      rxyz(16)=+r7(13,4)-r6(8,6)*two+r5(4,6)+r5(4,12)+r5(13,22)-r4(8,31)*two+r3(4,28) &
&             +r3(4,40)
      rxyz(17)=+r7(13,3)-r6(8,5)*two+r5(4,5)+r5(4,11)+r5(13,21)-r4(8,30)*two+r3(4,27) &
&             +r3(4,39)
      rxyz(18)=+r7(26,3)-r6(19,5)*two+r5(13,5)+r5(13,11)+r5(13,21)-r4(8,30)*two+r3(4,27) &
&             +r3(4,39)
      rxyz(19)=+r7(26,4)-r6(19,6)*two+r5(13,6)+r5(13,12)+r5(13,22)-r4(8,31)*two+r3(4,28) &
&             +r3(4,40)
      rxyz(20)=+r6(14,11)*four-r5(9,15)*eight+r4(5,16)*four+r4(5,28)*four
      eri(1,1,3,4)=r400+(+r7(9,3)*two+r7(9,4)*two-r6(5,5)*four-r6(5,6)*four+r5(2,5)*two &
&                  +r5(2,6)*two+r5(2,11)*two+r5(2,12)*two+r5(9,21)*six+r5(9,22)*six &
&                  -r4(5,30)*p12-r4(5,31)*p12+r3(2,27)*six+r3(2,28)*six+r3(2,39)*six &
&                  +r3(2,40)*six)*qx+(+r6(9,10)+r6(9,11)*four+r6(9,12)-r5(5,14)*two &
&                  -r5(5,15)*eight-r5(5,16)*two+r4(2,15)+r4(2,16)*four+r4(2,17)+r4(2,27) &
&                  +r4(2,28)*four+r4(2,29)+r4(9,39)+r4(9,40)*four+r4(9,41)-r3(5,44)*two &
&                  -r3(5,45)*eight-r3(5,46)*two+r2(4,38)+r2(4,39)*four+r2(4,40)+r2(4,53) &
&                  +r2(4,54)*four+r2(4,55))*xx+(+r5(9,23)*two+r5(9,24)*two-r4(5,32)*four &
&                  -r4(5,33)*four+r3(2,29)*two+r3(2,30)*two+r3(2,41)*two+r3(2,42)*two)*xxx &
&                  +rxyz(1)*xxxx
      eri(2,1,3,4)=r220+(+r7(18,4)*two-r6(12,6)*four+r5(7,6)*two+r5(7,12)*two &
&                  +r5(9,22)*two-r4(5,31)*four+r3(2,28)*two+r3(2,40)*two)*qx+rxyz(5)*xx
      eri(3,1,3,4)=r202+(+r7(20,4)*two-r6(14,6)*four+r5(9,6)*two+r5(9,12)*two &
&                  +r5(9,22)*two-r4(5,31)*four+r3(2,28)*two+r3(2,40)*two)*qx+(+r7(14,3)*two &
&                  -r6(9,5)*four+r5(5,5)*two+r5(5,11)*two+r5(14,21)*two-r4(9,30)*four &
&                  +r3(5,27)*two+r3(5,39)*two)*qz+(+r6(20,12)-r5(14,16)*two+r4(9,17)+r4(9,29) &
&                  +r4(9,41)-r3(5,46)*two+r2(4,40)+r2(4,55))*xx+rxyz(20)*xz+(+r6(9,10) &
&                  -r5(5,14)*two+r4(2,15)+r4(2,27)+r4(9,39)-r3(5,44)*two+r2(4,38)+r2(4,53)) &
&                  *zz+(+r5(14,24)*two-r4(9,33)*four+r3(5,30)*two+r3(5,42)*two)*xxz+( &
&                  +r5(9,23)*two-r4(5,32)*four+r3(2,29)*two+r3(2,41)*two)*xzz+rxyz(1)*xxzz
      eri(4,1,3,4)=r310+(+r7(13,3)+r7(13,4)*two-r6(8,5)*two-r6(8,6)*four+r5(4,5) &
&                  +r5(4,6)*two+r5(4,11)+r5(4,12)*two+r5(13,21)+r5(13,22)*two-r4(8,30)*two &
&                  -r4(8,31)*four+r3(4,27)+r3(4,28)*two+r3(4,39)+r3(4,40)*two)*qx+( &
&                  +r6(13,11)*two+r6(13,12)-r5(8,15)*four-r5(8,16)*two+r4(4,16)*two+r4(4,17) &
&                  +r4(4,28)*two+r4(4,29))*xx+rxyz(2)*xxx
      eri(5,1,3,4)=r301+(+r7(14,3)+r7(14,4)*two-r6(9,5)*two-r6(9,6)*four+r5(5,5) &
&                  +r5(5,6)*two+r5(5,11)+r5(5,12)*two+r5(14,21)+r5(14,22)*two-r4(9,30)*two &
&                  -r4(9,31)*four+r3(5,27)+r3(5,28)*two+r3(5,39)+r3(5,40)*two)*qx+(+r7(9,3) &
&                  -r6(5,5)*two+r5(2,5)+r5(2,11)+r5(9,21)*three-r4(5,30)*six+r3(2,27)*three &
&                  +r3(2,39)*three)*qz+(+r6(14,11)*two+r6(14,12)-r5(9,15)*four-r5(9,16)*two &
&                  +r4(5,16)*two+r4(5,17)+r4(5,28)*two+r4(5,29))*xx+(+r6(9,10)+r6(9,11)*two &
&                  -r5(5,14)*two-r5(5,15)*four+r4(2,15)+r4(2,16)*two+r4(2,27)+r4(2,28)*two &
&                  +r4(9,39)+r4(9,40)*two-r3(5,44)*two-r3(5,45)*four+r2(4,38)+r2(4,39)*two &
&                  +r2(4,53)+r2(4,54)*two)*xz+(+r5(14,24)-r4(9,33)*two+r3(5,30)+r3(5,42))*xxx &
&                  +(+r5(9,23)*two+r5(9,24)-r4(5,32)*four-r4(5,33)*two+r3(2,29)*two+r3(2,30) &
&                  +r3(2,41)*two+r3(2,42))*xxz+rxyz(1)*xxxz
      eri(6,1,3,4)=r211+(+r7(19,4)*two-r6(13,6)*four+r5(8,6)*two+r5(8,12)*two)*qx &
&                  +rxyz(17)*qz+(+r6(19,12)-r5(13,16)*two+r4(8,17)+r4(8,29))*xx+( &
&                  +r6(13,11)*two-r5(8,15)*four+r4(4,16)*two+r4(4,28)*two)*xz+rxyz(2)*xxz
      eri(1,2,3,4)=r220+(+r7(18,3)*two-r6(12,5)*four+r5(7,5)*two+r5(7,11)*two &
&                  +r5(9,21)*two-r4(5,30)*four+r3(2,27)*two+r3(2,39)*two)*qx+rxyz(4)*xx
      eri(2,2,3,4)=r040
      eri(3,2,3,4)=r022+(+r7(25,3)*two-r6(18,5)*four+r5(12,5)*two+r5(12,11)*two &
&                  +r5(14,21)*two-r4(9,30)*four+r3(5,27)*two+r3(5,39)*two)*qz+rxyz(4)*zz
      eri(4,2,3,4)=r130+rxyz(7)*qx
      eri(5,2,3,4)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,3,4)=r031+rxyz(7)*qz
      eri(1,3,3,4)=r202+(+r7(20,3)*two-r6(14,5)*four+r5(9,5)*two+r5(9,11)*two &
&                  +r5(9,21)*two-r4(5,30)*four+r3(2,27)*two+r3(2,39)*two)*qx+(+r7(14,4)*two &
&                  -r6(9,6)*four+r5(5,6)*two+r5(5,12)*two+r5(14,22)*two-r4(9,31)*four &
&                  +r3(5,28)*two+r3(5,40)*two)*qz+(+r6(20,10)-r5(14,14)*two+r4(9,15)+r4(9,27) &
&                  +r4(9,39)-r3(5,44)*two+r2(4,38)+r2(4,53))*xx+rxyz(20)*xz+(+r6(9,12) &
&                  -r5(5,16)*two+r4(2,17)+r4(2,29)+r4(9,41)-r3(5,46)*two+r2(4,40)+r2(4,55)) &
&                  *zz+(+r5(14,23)*two-r4(9,32)*four+r3(5,29)*two+r3(5,41)*two)*xxz+( &
&                  +r5(9,24)*two-r4(5,33)*four+r3(2,30)*two+r3(2,42)*two)*xzz+rxyz(1)*xxzz
      eri(2,3,3,4)=r022+(+r7(25,4)*two-r6(18,6)*four+r5(12,6)*two+r5(12,12)*two &
&                  +r5(14,22)*two-r4(9,31)*four+r3(5,28)*two+r3(5,40)*two)*qz+rxyz(5)*zz
      eri(3,3,3,4)=r004+(+r7(27,3)*two+r7(27,4)*two-r6(20,5)*four-r6(20,6)*four &
&                  +r5(14,5)*two+r5(14,6)*two+r5(14,11)*two+r5(14,12)*two+r5(14,21)*six &
&                  +r5(14,22)*six-r4(9,30)*p12-r4(9,31)*p12+r3(5,27)*six+r3(5,28)*six &
&                  +r3(5,39)*six+r3(5,40)*six)*qz+(+r6(20,10)+r6(20,11)*four+r6(20,12) &
&                  -r5(14,14)*two-r5(14,15)*eight-r5(14,16)*two+r4(9,15)+r4(9,16)*four &
&                  +r4(9,17)+r4(9,27)+r4(9,28)*four+r4(9,29)+r4(9,39)+r4(9,40)*four+r4(9,41) &
&                  -r3(5,44)*two-r3(5,45)*eight-r3(5,46)*two+r2(4,38)+r2(4,39)*four+r2(4,40) &
&                  +r2(4,53)+r2(4,54)*four+r2(4,55))*zz+(+r5(14,23)*two+r5(14,24)*two &
&                  -r4(9,32)*four-r4(9,33)*four+r3(5,29)*two+r3(5,30)*two+r3(5,41)*two &
&                  +r3(5,42)*two)*zzz+rxyz(1)*zzzz
      eri(4,3,3,4)=r112+rxyz(18)*qx+(+r7(19,4)*two-r6(13,6)*four+r5(8,6)*two+r5(8,12)*two &
&                  )*qz+(+r6(19,11)*two-r5(13,15)*four+r4(8,16)*two+r4(8,28)*two)*xz+( &
&                  +r6(13,12)-r5(8,16)*two+r4(4,17)+r4(4,29))*zz+rxyz(2)*xzz
      eri(5,3,3,4)=r103+(+r7(27,3)-r6(20,5)*two+r5(14,5)+r5(14,11)+r5(14,21)*three &
&                  -r4(9,30)*six+r3(5,27)*three+r3(5,39)*three)*qx+(+r7(20,3)+r7(20,4)*two &
&                  -r6(14,5)*two-r6(14,6)*four+r5(9,5)+r5(9,6)*two+r5(9,11)+r5(9,12)*two &
&                  +r5(9,21)+r5(9,22)*two-r4(5,30)*two-r4(5,31)*four+r3(2,27)+r3(2,28)*two &
&                  +r3(2,39)+r3(2,40)*two)*qz+(+r6(20,10)+r6(20,11)*two-r5(14,14)*two &
&                  -r5(14,15)*four+r4(9,15)+r4(9,16)*two+r4(9,27)+r4(9,28)*two+r4(9,39) &
&                  +r4(9,40)*two-r3(5,44)*two-r3(5,45)*four+r2(4,38)+r2(4,39)*two+r2(4,53) &
&                  +r2(4,54)*two)*xz+(+r6(14,11)*two+r6(14,12)-r5(9,15)*four-r5(9,16)*two &
&                  +r4(5,16)*two+r4(5,17)+r4(5,28)*two+r4(5,29))*zz+(+r5(14,23)*two+r5(14,24) &
&                  -r4(9,32)*four-r4(9,33)*two+r3(5,29)*two+r3(5,30)+r3(5,41)*two+r3(5,42)) &
&                  *xzz+(+r5(9,24)-r4(5,33)*two+r3(2,30)+r3(2,42))*zzz+rxyz(1)*xzzz
      eri(6,3,3,4)=r013+(+r7(26,3)+r7(26,4)*two-r6(19,5)*two-r6(19,6)*four+r5(13,5) &
&                  +r5(13,6)*two+r5(13,11)+r5(13,12)*two+r5(13,21)+r5(13,22)*two-r4(8,30)*two &
&                  -r4(8,31)*four+r3(4,27)+r3(4,28)*two+r3(4,39)+r3(4,40)*two)*qz+( &
&                  +r6(19,11)*two+r6(19,12)-r5(13,15)*four-r5(13,16)*two+r4(8,16)*two &
&                  +r4(8,17)+r4(8,28)*two+r4(8,29))*zz+rxyz(2)*zzz
      eri(1,4,3,4)=r310+(+r7(13,3)*two+r7(13,4)-r6(8,5)*four-r6(8,6)*two+r5(4,5)*two &
&                  +r5(4,6)+r5(4,11)*two+r5(4,12)+r5(13,21)*two+r5(13,22)-r4(8,30)*four &
&                  -r4(8,31)*two+r3(4,27)*two+r3(4,28)+r3(4,39)*two+r3(4,40))*qx+(+r6(13,10) &
&                  +r6(13,11)*two-r5(8,14)*two-r5(8,15)*four+r4(4,15)+r4(4,16)*two+r4(4,27) &
&                  +r4(4,28)*two)*xx+rxyz(3)*xxx
      eri(2,4,3,4)=r130+rxyz(8)*qx
      eri(3,4,3,4)=r112+rxyz(19)*qx+(+r7(19,3)*two-r6(13,5)*four+r5(8,5)*two+r5(8,11)*two &
&                  )*qz+(+r6(19,11)*two-r5(13,15)*four+r4(8,16)*two+r4(8,28)*two)*xz+( &
&                  +r6(13,10)-r5(8,14)*two+r4(4,15)+r4(4,27))*zz+rxyz(3)*xzz
      eri(4,4,3,4)=r220+(+r7(18,3)+r7(18,4)-r6(12,5)*two-r6(12,6)*two+r5(7,5)+r5(7,6) &
&                  +r5(7,11)+r5(7,12)+r5(9,21)+r5(9,22)-r4(5,30)*two-r4(5,31)*two+r3(2,27) &
&                  +r3(2,28)+r3(2,39)+r3(2,40))*qx+rxyz(6)*xx
      eri(5,4,3,4)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(13,10)+r6(13,11) &
&                  -r5(8,14)*two-r5(8,15)*two+r4(4,15)+r4(4,16)+r4(4,27)+r4(4,28))*xz+rxyz(3) &
&                  *xxz
      eri(6,4,3,4)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,3,4)=r301+(+r7(14,3)*two+r7(14,4)-r6(9,5)*four-r6(9,6)*two+r5(5,5)*two &
&                  +r5(5,6)+r5(5,11)*two+r5(5,12)+r5(14,21)*two+r5(14,22)-r4(9,30)*four &
&                  -r4(9,31)*two+r3(5,27)*two+r3(5,28)+r3(5,39)*two+r3(5,40))*qx+(+r7(9,4) &
&                  -r6(5,6)*two+r5(2,6)+r5(2,12)+r5(9,22)*three-r4(5,31)*six+r3(2,28)*three &
&                  +r3(2,40)*three)*qz+(+r6(14,10)+r6(14,11)*two-r5(9,14)*two-r5(9,15)*four &
&                  +r4(5,15)+r4(5,16)*two+r4(5,27)+r4(5,28)*two)*xx+(+r6(9,11)*two+r6(9,12) &
&                  -r5(5,15)*four-r5(5,16)*two+r4(2,16)*two+r4(2,17)+r4(2,28)*two+r4(2,29) &
&                  +r4(9,40)*two+r4(9,41)-r3(5,45)*four-r3(5,46)*two+r2(4,39)*two+r2(4,40) &
&                  +r2(4,54)*two+r2(4,55))*xz+(+r5(14,23)-r4(9,32)*two+r3(5,29)+r3(5,41))*xxx &
&                  +(+r5(9,23)+r5(9,24)*two-r4(5,32)*two-r4(5,33)*four+r3(2,29)+r3(2,30)*two &
&                  +r3(2,41)+r3(2,42)*two)*xxz+rxyz(1)*xxxz
      eri(2,5,3,4)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,3,4)=r103+(+r7(27,4)-r6(20,6)*two+r5(14,6)+r5(14,12)+r5(14,22)*three &
&                  -r4(9,31)*six+r3(5,28)*three+r3(5,40)*three)*qx+(+r7(20,3)*two+r7(20,4) &
&                  -r6(14,5)*four-r6(14,6)*two+r5(9,5)*two+r5(9,6)+r5(9,11)*two+r5(9,12) &
&                  +r5(9,21)*two+r5(9,22)-r4(5,30)*four-r4(5,31)*two+r3(2,27)*two+r3(2,28) &
&                  +r3(2,39)*two+r3(2,40))*qz+(+r6(20,11)*two+r6(20,12)-r5(14,15)*four &
&                  -r5(14,16)*two+r4(9,16)*two+r4(9,17)+r4(9,28)*two+r4(9,29)+r4(9,40)*two &
&                  +r4(9,41)-r3(5,45)*four-r3(5,46)*two+r2(4,39)*two+r2(4,40)+r2(4,54)*two &
&                  +r2(4,55))*xz+(+r6(14,10)+r6(14,11)*two-r5(9,14)*two-r5(9,15)*four &
&                  +r4(5,15)+r4(5,16)*two+r4(5,27)+r4(5,28)*two)*zz+(+r5(14,23)+r5(14,24)*two &
&                  -r4(9,32)*two-r4(9,33)*four+r3(5,29)+r3(5,30)*two+r3(5,41)+r3(5,42)*two) &
&                  *xzz+(+r5(9,23)-r4(5,32)*two+r3(2,29)+r3(2,41))*zzz+rxyz(1)*xzzz
      eri(4,5,3,4)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(13,11)+r6(13,12) &
&                  -r5(8,15)*two-r5(8,16)*two+r4(4,16)+r4(4,17)+r4(4,28)+r4(4,29))*xz+rxyz(2) &
&                  *xxz
      eri(5,5,3,4)=r202+(+r7(20,3)+r7(20,4)-r6(14,5)*two-r6(14,6)*two+r5(9,5)+r5(9,6) &
&                  +r5(9,11)+r5(9,12)+r5(9,21)+r5(9,22)-r4(5,30)*two-r4(5,31)*two+r3(2,27) &
&                  +r3(2,28)+r3(2,39)+r3(2,40))*qx+(+r7(14,3)+r7(14,4)-r6(9,5)*two &
&                  -r6(9,6)*two+r5(5,5)+r5(5,6)+r5(5,11)+r5(5,12)+r5(14,21)+r5(14,22) &
&                  -r4(9,30)*two-r4(9,31)*two+r3(5,27)+r3(5,28)+r3(5,39)+r3(5,40))*qz+( &
&                  +r6(20,11)-r5(14,15)*two+r4(9,16)+r4(9,28)+r4(9,40)-r3(5,45)*two+r2(4,39) &
&                  +r2(4,54))*xx+(+r6(14,10)+r6(14,11)*two+r6(14,12)-r5(9,14)*two &
&                  -r5(9,15)*four-r5(9,16)*two+r4(5,15)+r4(5,16)*two+r4(5,17)+r4(5,27) &
&                  +r4(5,28)*two+r4(5,29))*xz+(+r6(9,11)-r5(5,15)*two+r4(2,16)+r4(2,28) &
&                  +r4(9,40)-r3(5,45)*two+r2(4,39)+r2(4,54))*zz+(+r5(14,23)+r5(14,24) &
&                  -r4(9,32)*two-r4(9,33)*two+r3(5,29)+r3(5,30)+r3(5,41)+r3(5,42))*xxz+( &
&                  +r5(9,23)+r5(9,24)-r4(5,32)*two-r4(5,33)*two+r3(2,29)+r3(2,30)+r3(2,41) &
&                  +r3(2,42))*xzz+rxyz(1)*xxzz
      eri(6,5,3,4)=r112+rxyz(19)*qx+(+r7(19,3)+r7(19,4)-r6(13,5)*two-r6(13,6)*two+r5(8,5) &
&                  +r5(8,6)+r5(8,11)+r5(8,12))*qz+(+r6(19,11)+r6(19,12)-r5(13,15)*two &
&                  -r5(13,16)*two+r4(8,16)+r4(8,17)+r4(8,28)+r4(8,29))*xz+rxyz(11)*zz+rxyz(2) &
&                  *xzz
      eri(1,6,3,4)=r211+(+r7(19,3)*two-r6(13,5)*four+r5(8,5)*two+r5(8,11)*two)*qx &
&                  +rxyz(16)*qz+(+r6(19,10)-r5(13,14)*two+r4(8,15)+r4(8,27))*xx+( &
&                  +r6(13,11)*two-r5(8,15)*four+r4(4,16)*two+r4(4,28)*two)*xz+rxyz(3)*xxz
      eri(2,6,3,4)=r031+rxyz(8)*qz
      eri(3,6,3,4)=r013+(+r7(26,3)*two+r7(26,4)-r6(19,5)*four-r6(19,6)*two+r5(13,5)*two &
&                  +r5(13,6)+r5(13,11)*two+r5(13,12)+r5(13,21)*two+r5(13,22)-r4(8,30)*four &
&                  -r4(8,31)*two+r3(4,27)*two+r3(4,28)+r3(4,39)*two+r3(4,40))*qz+(+r6(19,10) &
&                  +r6(19,11)*two-r5(13,14)*two-r5(13,15)*four+r4(8,15)+r4(8,16)*two+r4(8,27) &
&                  +r4(8,28)*two)*zz+rxyz(3)*zzz
      eri(4,6,3,4)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,3,4)=r112+rxyz(18)*qx+(+r7(19,3)+r7(19,4)-r6(13,5)*two-r6(13,6)*two+r5(8,5) &
&                  +r5(8,6)+r5(8,11)+r5(8,12))*qz+(+r6(19,10)+r6(19,11)-r5(13,14)*two &
&                  -r5(13,15)*two+r4(8,15)+r4(8,16)+r4(8,27)+r4(8,28))*xz+rxyz(11)*zz+rxyz(3) &
&                  *xzz
      eri(6,6,3,4)=r022+(+r7(25,3)+r7(25,4)-r6(18,5)*two-r6(18,6)*two+r5(12,5)+r5(12,6) &
&                  +r5(12,11)+r5(12,12)+r5(14,21)+r5(14,22)-r4(9,30)*two-r4(9,31)*two &
&                  +r3(5,27)+r3(5,28)+r3(5,39)+r3(5,40))*qz+rxyz(6)*zz
!
      do i= 1,6
        do j= 1,6
          eri(j,i,4,4)= eri(j,i,2,1)
        enddo
      enddo
      do i= 1,6
        do j= 1,6
          eri(j,i,5,4)= eri(j,i,6,1)
        enddo
      enddo
      do i= 1,6
        do j= 1,6
          eri(j,i,6,4)= eri(j,i,5,2)
        enddo
      enddo
!
      r400= r8(3)-r7(1,2)+r6(3,4)*three+r6(3,9)*six-r5(1,4)*three-r5(1,17)*six+r4(3,26)*p18 &
&          +r4(3,38)*three-r3(1,23)*p18-r3(1,48)*three+r2(5,52)*nine-r1(1,36)*nine
      r310= r8(5)-r7(2,2)+r6(5,4)*three+r6(5,9)*three-r5(2,4)*three-r5(2,17)*three &
&          +r4(5,26)*nine-r3(2,23)*nine
      r301= r8(6)-r7(3,2)+r6(6,4)*three+r6(6,9)*three-r5(3,4)*three-r5(3,17)*three &
&          +r4(6,26)*nine-r3(3,23)*nine
      r220= r8(8)-r7(4,2)+r6(8,4)*three+r6(3,9)+r6(8,9)-r5(4,4)*three-r5(1,17)-r5(4,17) &
&          +r4(3,26)*three+r4(8,26)*three+r4(3,38)-r3(1,23)*three-r3(4,23)*three-r3(1,48) &
&          +r2(5,52)*three-r1(1,36)*three
      r211= r8(9)-r7(5,2)+r6(9,4)*three+r6(9,9)-r5(5,4)*three-r5(5,17)+r4(9,26)*three &
&          -r3(5,23)*three
      r202= r8(10)-r7(6,2)+r6(10,4)*three+r6(3,9)+r6(10,9)-r5(6,4)*three-r5(1,17)-r5(6,17) &
&          +r4(3,26)*three+r4(10,26)*three+r4(3,38)-r3(1,23)*three-r3(6,23)*three-r3(1,48) &
&          +r2(5,52)*three-r1(1,36)*three
      r130= r8(12)-r7(7,2)+r6(12,4)*three+r6(5,9)*three-r5(7,4)*three-r5(2,17)*three &
&          +r4(5,26)*nine-r3(2,23)*nine
      r121= r8(13)-r7(8,2)+r6(13,4)*three+r6(6,9)-r5(8,4)*three-r5(3,17)+r4(6,26)*three &
&          -r3(3,23)*three
      r112= r8(14)-r7(9,2)+r6(14,4)*three+r6(5,9)-r5(9,4)*three-r5(2,17)+r4(5,26)*three &
&          -r3(2,23)*three
      r103= r8(15)-r7(10,2)+r6(15,4)*three+r6(6,9)*three-r5(10,4)*three-r5(3,17)*three &
&          +r4(6,26)*nine-r3(3,23)*nine
      r040= r8(17)-r7(11,2)+r6(17,4)*three+r6(8,9)*six-r5(11,4)*three-r5(4,17)*six &
&          +r4(8,26)*p18+r4(3,38)*three-r3(4,23)*p18-r3(1,48)*three+r2(5,52)*nine &
&          -r1(1,36)*nine
      r031= r8(18)-r7(12,2)+r6(18,4)*three+r6(9,9)*three-r5(12,4)*three-r5(5,17)*three &
&          +r4(9,26)*nine-r3(5,23)*nine
      r022= r8(19)-r7(13,2)+r6(19,4)*three+r6(8,9)+r6(10,9)-r5(13,4)*three-r5(4,17) &
&          -r5(6,17)+r4(8,26)*three+r4(10,26)*three+r4(3,38)-r3(4,23)*three-r3(6,23)*three &
&          -r3(1,48)+r2(5,52)*three-r1(1,36)*three
      r013= r8(20)-r7(14,2)+r6(20,4)*three+r6(9,9)*three-r5(14,4)*three-r5(5,17)*three &
&          +r4(9,26)*nine-r3(5,23)*nine
      r004= r8(21)-r7(15,2)+r6(21,4)*three+r6(10,9)*six-r5(15,4)*three-r5(6,17)*six &
&          +r4(10,26)*p18+r4(3,38)*three-r3(6,23)*p18-r3(1,48)*three+r2(5,52)*nine &
&          -r1(1,36)*nine
      rxyz(1)=+r4(3,42)-r3(1,52)+r2(5,56)*three-r1(1,40)*three
      rxyz(2)=+r5(5,24)-r4(2,37)+r3(5,42)*three-r2(4,36)*three
      rxyz(3)=+r5(5,23)-r4(2,36)+r3(5,41)*three-r2(4,35)*three
      rxyz(4)=+r6(8,10)-r5(4,18)+r4(8,27)*three+r4(3,39)-r3(4,24)*three-r3(1,49) &
&             +r2(5,53)*three-r1(1,37)*three
      rxyz(5)=+r6(8,12)-r5(4,20)+r4(8,29)*three+r4(3,41)-r3(4,26)*three-r3(1,51) &
&             +r2(5,55)*three-r1(1,39)*three
      rxyz(6)=+r6(8,11)-r5(4,19)+r4(8,28)*three+r4(3,40)-r3(4,25)*three-r3(1,50) &
&             +r2(5,54)*three-r1(1,38)*three
      rxyz(7)=+r7(12,3)-r6(7,7)+r5(12,11)*three+r5(5,21)*three-r4(7,12)*three &
&             -r4(2,34)*three+r3(5,39)*nine-r2(4,33)*nine
      rxyz(8)=+r7(12,4)-r6(7,8)+r5(12,12)*three+r5(5,22)*three-r4(7,13)*three &
&             -r4(2,35)*three+r3(5,40)*nine-r2(4,34)*nine
      rxyz(9)=+r7(9,3)+r7(9,4)-r6(5,7)-r6(5,8)+r5(9,11)*three+r5(9,12)*three &
&             -r4(5,12)*three-r4(5,13)*three
      rxyz(10)=+r6(9,11)-r5(5,19)+r4(9,28)*three-r3(5,25)*three
      rxyz(11)=+r6(5,11)-r5(2,19)+r4(5,28)*three-r3(2,25)*three
      rxyz(12)=+r7(13,3)-r6(8,7)+r5(13,11)*three+r5(6,21)-r4(8,12)*three-r4(3,34) &
&             +r3(6,39)*three-r2(5,33)*three
      rxyz(13)=+r7(13,4)-r6(8,8)+r5(13,12)*three+r5(6,22)-r4(8,13)*three-r4(3,35) &
&             +r3(6,40)*three-r2(5,34)*three
      rxyz(14)=+r7(8,4)-r6(4,8)+r5(8,12)*three+r5(3,22)-r4(4,13)*three-r4(1,35) &
&             +r3(3,40)*three-r2(1,34)*three
      rxyz(15)=+r7(8,3)-r6(4,7)+r5(8,11)*three+r5(3,21)-r4(4,12)*three-r4(1,34) &
&             +r3(3,39)*three-r2(1,33)*three
      rxyz(16)=+r7(5,4)-r6(2,8)+r5(5,12)*three+r5(5,22)-r4(2,13)*three-r4(2,35) &
&             +r3(5,40)*three-r2(4,34)*three
      rxyz(17)=+r7(5,3)-r6(2,7)+r5(5,11)*three+r5(5,21)-r4(2,12)*three-r4(2,34) &
&             +r3(5,39)*three-r2(4,33)*three
      rxyz(18)=+r7(14,3)-r6(9,7)+r5(14,11)*three+r5(5,21)-r4(9,12)*three-r4(2,34) &
&             +r3(5,39)*three-r2(4,33)*three
      rxyz(19)=+r7(14,4)-r6(9,8)+r5(14,12)*three+r5(5,22)-r4(9,13)*three-r4(2,35) &
&             +r3(5,40)*three-r2(4,34)*three
      rxyz(20)=+r6(6,11)*four-r5(3,19)*four+r4(6,28)*p12-r3(3,25)*p12
      eri(1,1,1,5)=r400+(+r7(3,3)*two+r7(3,4)*two-r6(1,7)*two-r6(1,8)*two+r5(3,11)*six &
&                  +r5(3,12)*six+r5(3,21)*six+r5(3,22)*six-r4(1,12)*six-r4(1,13)*six &
&                  -r4(1,34)*six-r4(1,35)*six+r3(3,39)*p18+r3(3,40)*p18-r2(1,33)*p18 &
&                  -r2(1,34)*p18)*qx+(+r6(3,10)+r6(3,11)*four+r6(3,12)-r5(1,18)-r5(1,19)*four &
&                  -r5(1,20)+r4(3,27)*three+r4(3,28)*p12+r4(3,29)*three+r4(3,39) &
&                  +r4(3,40)*four+r4(3,41)-r3(1,24)*three-r3(1,25)*p12-r3(1,26)*three &
&                  -r3(1,49)-r3(1,50)*four-r3(1,51)+r2(5,53)*three+r2(5,54)*p12 &
&                  +r2(5,55)*three-r1(1,37)*three-r1(1,38)*p12-r1(1,39)*three)*xx+( &
&                  +r5(3,23)*two+r5(3,24)*two-r4(1,36)*two-r4(1,37)*two+r3(3,41)*six &
&                  +r3(3,42)*six-r2(1,35)*six-r2(1,36)*six)*xxx+rxyz(1)*xxxx
      eri(2,1,1,5)=r220+(+r7(8,4)*two-r6(4,8)*two+r5(8,12)*six+r5(3,22)*two-r4(4,13)*six &
&                  -r4(1,35)*two+r3(3,40)*six-r2(1,34)*six)*qx+rxyz(5)*xx
      eri(3,1,1,5)=r202+(+r7(10,4)*two-r6(6,8)*two+r5(10,12)*six+r5(3,22)*two &
&                  -r4(6,13)*six-r4(1,35)*two+r3(3,40)*six-r2(1,34)*six)*qx+(+r7(6,3)*two &
&                  -r6(3,7)*two+r5(6,11)*six+r5(6,21)*two-r4(3,12)*six-r4(3,34)*two &
&                  +r3(6,39)*six-r2(5,33)*six)*qz+(+r6(10,12)-r5(6,20)+r4(10,29)*three &
&                  +r4(3,41)-r3(6,26)*three-r3(1,51)+r2(5,55)*three-r1(1,39)*three)*xx &
&                  +rxyz(20)*xz+(+r6(3,10)-r5(1,18)+r4(3,27)*three+r4(3,39)-r3(1,24)*three &
&                  -r3(1,49)+r2(5,53)*three-r1(1,37)*three)*zz+(+r5(6,24)*two-r4(3,37)*two &
&                  +r3(6,42)*six-r2(5,36)*six)*xxz+(+r5(3,23)*two-r4(1,36)*two+r3(3,41)*six &
&                  -r2(1,35)*six)*xzz+rxyz(1)*xxzz
      eri(4,1,1,5)=r310+(+r7(5,3)+r7(5,4)*two-r6(2,7)-r6(2,8)*two+r5(5,11)*three &
&                  +r5(5,12)*six+r5(5,21)+r5(5,22)*two-r4(2,12)*three-r4(2,13)*six-r4(2,34) &
&                  -r4(2,35)*two+r3(5,39)*three+r3(5,40)*six-r2(4,33)*three-r2(4,34)*six)*qx &
&                  +(+r6(5,11)*two+r6(5,12)-r5(2,19)*two-r5(2,20)+r4(5,28)*six+r4(5,29)*three &
&                  -r3(2,25)*six-r3(2,26)*three)*xx+rxyz(2)*xxx
      eri(5,1,1,5)=r301+(+r7(6,3)+r7(6,4)*two-r6(3,7)-r6(3,8)*two+r5(6,11)*three &
&                  +r5(6,12)*six+r5(6,21)+r5(6,22)*two-r4(3,12)*three-r4(3,13)*six-r4(3,34) &
&                  -r4(3,35)*two+r3(6,39)*three+r3(6,40)*six-r2(5,33)*three-r2(5,34)*six)*qx &
&                  +(+r7(3,3)-r6(1,7)+r5(3,11)*three+r5(3,21)*three-r4(1,12)*three &
&                  -r4(1,34)*three+r3(3,39)*nine-r2(1,33)*nine)*qz+(+r6(6,11)*two+r6(6,12) &
&                  -r5(3,19)*two-r5(3,20)+r4(6,28)*six+r4(6,29)*three-r3(3,25)*six &
&                  -r3(3,26)*three)*xx+(+r6(3,10)+r6(3,11)*two-r5(1,18)-r5(1,19)*two &
&                  +r4(3,27)*three+r4(3,28)*six+r4(3,39)+r4(3,40)*two-r3(1,24)*three &
&                  -r3(1,25)*six-r3(1,49)-r3(1,50)*two+r2(5,53)*three+r2(5,54)*six &
&                  -r1(1,37)*three-r1(1,38)*six)*xz+(+r5(6,24)-r4(3,37)+r3(6,42)*three &
&                  -r2(5,36)*three)*xxx+(+r5(3,23)*two+r5(3,24)-r4(1,36)*two-r4(1,37) &
&                  +r3(3,41)*six+r3(3,42)*three-r2(1,35)*six-r2(1,36)*three)*xxz+rxyz(1)*xxxz
      eri(6,1,1,5)=r211+(+r7(9,4)*two-r6(5,8)*two+r5(9,12)*six-r4(5,13)*six)*qx+rxyz(17) &
&                  *qz+(+r6(9,12)-r5(5,20)+r4(9,29)*three-r3(5,26)*three)*xx+(+r6(5,11)*two &
&                  -r5(2,19)*two+r4(5,28)*six-r3(2,25)*six)*xz+rxyz(2)*xxz
      eri(1,2,1,5)=r220+(+r7(8,3)*two-r6(4,7)*two+r5(8,11)*six+r5(3,21)*two-r4(4,12)*six &
&                  -r4(1,34)*two+r3(3,39)*six-r2(1,33)*six)*qx+rxyz(4)*xx
      eri(2,2,1,5)=r040
      eri(3,2,1,5)=r022+(+r7(13,3)*two-r6(8,7)*two+r5(13,11)*six+r5(6,21)*two &
&                  -r4(8,12)*six-r4(3,34)*two+r3(6,39)*six-r2(5,33)*six)*qz+rxyz(4)*zz
      eri(4,2,1,5)=r130+rxyz(7)*qx
      eri(5,2,1,5)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,1,5)=r031+rxyz(7)*qz
      eri(1,3,1,5)=r202+(+r7(10,3)*two-r6(6,7)*two+r5(10,11)*six+r5(3,21)*two &
&                  -r4(6,12)*six-r4(1,34)*two+r3(3,39)*six-r2(1,33)*six)*qx+(+r7(6,4)*two &
&                  -r6(3,8)*two+r5(6,12)*six+r5(6,22)*two-r4(3,13)*six-r4(3,35)*two &
&                  +r3(6,40)*six-r2(5,34)*six)*qz+(+r6(10,10)-r5(6,18)+r4(10,27)*three &
&                  +r4(3,39)-r3(6,24)*three-r3(1,49)+r2(5,53)*three-r1(1,37)*three)*xx &
&                  +rxyz(20)*xz+(+r6(3,12)-r5(1,20)+r4(3,29)*three+r4(3,41)-r3(1,26)*three &
&                  -r3(1,51)+r2(5,55)*three-r1(1,39)*three)*zz+(+r5(6,23)*two-r4(3,36)*two &
&                  +r3(6,41)*six-r2(5,35)*six)*xxz+(+r5(3,24)*two-r4(1,37)*two+r3(3,42)*six &
&                  -r2(1,36)*six)*xzz+rxyz(1)*xxzz
      eri(2,3,1,5)=r022+(+r7(13,4)*two-r6(8,8)*two+r5(13,12)*six+r5(6,22)*two &
&                  -r4(8,13)*six-r4(3,35)*two+r3(6,40)*six-r2(5,34)*six)*qz+rxyz(5)*zz
      eri(3,3,1,5)=r004+(+r7(15,3)*two+r7(15,4)*two-r6(10,7)*two-r6(10,8)*two &
&                  +r5(15,11)*six+r5(15,12)*six+r5(6,21)*six+r5(6,22)*six-r4(10,12)*six &
&                  -r4(10,13)*six-r4(3,34)*six-r4(3,35)*six+r3(6,39)*p18+r3(6,40)*p18 &
&                  -r2(5,33)*p18-r2(5,34)*p18)*qz+(+r6(10,10)+r6(10,11)*four+r6(10,12) &
&                  -r5(6,18)-r5(6,19)*four-r5(6,20)+r4(10,27)*three+r4(10,28)*p12 &
&                  +r4(10,29)*three+r4(3,39)+r4(3,40)*four+r4(3,41)-r3(6,24)*three &
&                  -r3(6,25)*p12-r3(6,26)*three-r3(1,49)-r3(1,50)*four-r3(1,51) &
&                  +r2(5,53)*three+r2(5,54)*p12+r2(5,55)*three-r1(1,37)*three-r1(1,38)*p12 &
&                  -r1(1,39)*three)*zz+(+r5(6,23)*two+r5(6,24)*two-r4(3,36)*two-r4(3,37)*two &
&                  +r3(6,41)*six+r3(6,42)*six-r2(5,35)*six-r2(5,36)*six)*zzz+rxyz(1)*zzzz
      eri(4,3,1,5)=r112+rxyz(18)*qx+(+r7(9,4)*two-r6(5,8)*two+r5(9,12)*six-r4(5,13)*six) &
&                  *qz+(+r6(9,11)*two-r5(5,19)*two+r4(9,28)*six-r3(5,25)*six)*xz+(+r6(5,12) &
&                  -r5(2,20)+r4(5,29)*three-r3(2,26)*three)*zz+rxyz(2)*xzz
      eri(5,3,1,5)=r103+(+r7(15,3)-r6(10,7)+r5(15,11)*three+r5(6,21)*three &
&                  -r4(10,12)*three-r4(3,34)*three+r3(6,39)*nine-r2(5,33)*nine)*qx+(+r7(10,3) &
&                  +r7(10,4)*two-r6(6,7)-r6(6,8)*two+r5(10,11)*three+r5(10,12)*six+r5(3,21) &
&                  +r5(3,22)*two-r4(6,12)*three-r4(6,13)*six-r4(1,34)-r4(1,35)*two &
&                  +r3(3,39)*three+r3(3,40)*six-r2(1,33)*three-r2(1,34)*six)*qz+(+r6(10,10) &
&                  +r6(10,11)*two-r5(6,18)-r5(6,19)*two+r4(10,27)*three+r4(10,28)*six &
&                  +r4(3,39)+r4(3,40)*two-r3(6,24)*three-r3(6,25)*six-r3(1,49)-r3(1,50)*two &
&                  +r2(5,53)*three+r2(5,54)*six-r1(1,37)*three-r1(1,38)*six)*xz+( &
&                  +r6(6,11)*two+r6(6,12)-r5(3,19)*two-r5(3,20)+r4(6,28)*six+r4(6,29)*three &
&                  -r3(3,25)*six-r3(3,26)*three)*zz+(+r5(6,23)*two+r5(6,24)-r4(3,36)*two &
&                  -r4(3,37)+r3(6,41)*six+r3(6,42)*three-r2(5,35)*six-r2(5,36)*three)*xzz+( &
&                  +r5(3,24)-r4(1,37)+r3(3,42)*three-r2(1,36)*three)*zzz+rxyz(1)*xzzz
      eri(6,3,1,5)=r013+(+r7(14,3)+r7(14,4)*two-r6(9,7)-r6(9,8)*two+r5(14,11)*three &
&                  +r5(14,12)*six+r5(5,21)+r5(5,22)*two-r4(9,12)*three-r4(9,13)*six-r4(2,34) &
&                  -r4(2,35)*two+r3(5,39)*three+r3(5,40)*six-r2(4,33)*three-r2(4,34)*six)*qz &
&                  +(+r6(9,11)*two+r6(9,12)-r5(5,19)*two-r5(5,20)+r4(9,28)*six+r4(9,29)*three &
&                  -r3(5,25)*six-r3(5,26)*three)*zz+rxyz(2)*zzz
      eri(1,4,1,5)=r310+(+r7(5,3)*two+r7(5,4)-r6(2,7)*two-r6(2,8)+r5(5,11)*six &
&                  +r5(5,12)*three+r5(5,21)*two+r5(5,22)-r4(2,12)*six-r4(2,13)*three &
&                  -r4(2,34)*two-r4(2,35)+r3(5,39)*six+r3(5,40)*three-r2(4,33)*six &
&                  -r2(4,34)*three)*qx+(+r6(5,10)+r6(5,11)*two-r5(2,18)-r5(2,19)*two &
&                  +r4(5,27)*three+r4(5,28)*six-r3(2,24)*three-r3(2,25)*six)*xx+rxyz(3)*xxx
      eri(2,4,1,5)=r130+rxyz(8)*qx
      eri(3,4,1,5)=r112+rxyz(19)*qx+(+r7(9,3)*two-r6(5,7)*two+r5(9,11)*six-r4(5,12)*six) &
&                  *qz+(+r6(9,11)*two-r5(5,19)*two+r4(9,28)*six-r3(5,25)*six)*xz+(+r6(5,10) &
&                  -r5(2,18)+r4(5,27)*three-r3(2,24)*three)*zz+rxyz(3)*xzz
      eri(4,4,1,5)=r220+(+r7(8,3)+r7(8,4)-r6(4,7)-r6(4,8)+r5(8,11)*three+r5(8,12)*three &
&                  +r5(3,21)+r5(3,22)-r4(4,12)*three-r4(4,13)*three-r4(1,34)-r4(1,35) &
&                  +r3(3,39)*three+r3(3,40)*three-r2(1,33)*three-r2(1,34)*three)*qx+rxyz(6) &
&                  *xx
      eri(5,4,1,5)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(5,10)+r6(5,11)-r5(2,18) &
&                  -r5(2,19)+r4(5,27)*three+r4(5,28)*three-r3(2,24)*three-r3(2,25)*three)*xz &
&                  +rxyz(3)*xxz
      eri(6,4,1,5)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,1,5)=r301+(+r7(6,3)*two+r7(6,4)-r6(3,7)*two-r6(3,8)+r5(6,11)*six &
&                  +r5(6,12)*three+r5(6,21)*two+r5(6,22)-r4(3,12)*six-r4(3,13)*three &
&                  -r4(3,34)*two-r4(3,35)+r3(6,39)*six+r3(6,40)*three-r2(5,33)*six &
&                  -r2(5,34)*three)*qx+(+r7(3,4)-r6(1,8)+r5(3,12)*three+r5(3,22)*three &
&                  -r4(1,13)*three-r4(1,35)*three+r3(3,40)*nine-r2(1,34)*nine)*qz+(+r6(6,10) &
&                  +r6(6,11)*two-r5(3,18)-r5(3,19)*two+r4(6,27)*three+r4(6,28)*six &
&                  -r3(3,24)*three-r3(3,25)*six)*xx+(+r6(3,11)*two+r6(3,12)-r5(1,19)*two &
&                  -r5(1,20)+r4(3,28)*six+r4(3,29)*three+r4(3,40)*two+r4(3,41)-r3(1,25)*six &
&                  -r3(1,26)*three-r3(1,50)*two-r3(1,51)+r2(5,54)*six+r2(5,55)*three &
&                  -r1(1,38)*six-r1(1,39)*three)*xz+(+r5(6,23)-r4(3,36)+r3(6,41)*three &
&                  -r2(5,35)*three)*xxx+(+r5(3,23)+r5(3,24)*two-r4(1,36)-r4(1,37)*two &
&                  +r3(3,41)*three+r3(3,42)*six-r2(1,35)*three-r2(1,36)*six)*xxz+rxyz(1)*xxxz
      eri(2,5,1,5)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,1,5)=r103+(+r7(15,4)-r6(10,8)+r5(15,12)*three+r5(6,22)*three &
&                  -r4(10,13)*three-r4(3,35)*three+r3(6,40)*nine-r2(5,34)*nine)*qx+( &
&                  +r7(10,3)*two+r7(10,4)-r6(6,7)*two-r6(6,8)+r5(10,11)*six+r5(10,12)*three &
&                  +r5(3,21)*two+r5(3,22)-r4(6,12)*six-r4(6,13)*three-r4(1,34)*two-r4(1,35) &
&                  +r3(3,39)*six+r3(3,40)*three-r2(1,33)*six-r2(1,34)*three)*qz+( &
&                  +r6(10,11)*two+r6(10,12)-r5(6,19)*two-r5(6,20)+r4(10,28)*six &
&                  +r4(10,29)*three+r4(3,40)*two+r4(3,41)-r3(6,25)*six-r3(6,26)*three &
&                  -r3(1,50)*two-r3(1,51)+r2(5,54)*six+r2(5,55)*three-r1(1,38)*six &
&                  -r1(1,39)*three)*xz+(+r6(6,10)+r6(6,11)*two-r5(3,18)-r5(3,19)*two &
&                  +r4(6,27)*three+r4(6,28)*six-r3(3,24)*three-r3(3,25)*six)*zz+(+r5(6,23) &
&                  +r5(6,24)*two-r4(3,36)-r4(3,37)*two+r3(6,41)*three+r3(6,42)*six &
&                  -r2(5,35)*three-r2(5,36)*six)*xzz+(+r5(3,23)-r4(1,36)+r3(3,41)*three &
&                  -r2(1,35)*three)*zzz+rxyz(1)*xzzz
      eri(4,5,1,5)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(5,11)+r6(5,12)-r5(2,19) &
&                  -r5(2,20)+r4(5,28)*three+r4(5,29)*three-r3(2,25)*three-r3(2,26)*three)*xz &
&                  +rxyz(2)*xxz
      eri(5,5,1,5)=r202+(+r7(10,3)+r7(10,4)-r6(6,7)-r6(6,8)+r5(10,11)*three &
&                  +r5(10,12)*three+r5(3,21)+r5(3,22)-r4(6,12)*three-r4(6,13)*three-r4(1,34) &
&                  -r4(1,35)+r3(3,39)*three+r3(3,40)*three-r2(1,33)*three-r2(1,34)*three)*qx &
&                  +(+r7(6,3)+r7(6,4)-r6(3,7)-r6(3,8)+r5(6,11)*three+r5(6,12)*three+r5(6,21) &
&                  +r5(6,22)-r4(3,12)*three-r4(3,13)*three-r4(3,34)-r4(3,35)+r3(6,39)*three &
&                  +r3(6,40)*three-r2(5,33)*three-r2(5,34)*three)*qz+(+r6(10,11)-r5(6,19) &
&                  +r4(10,28)*three+r4(3,40)-r3(6,25)*three-r3(1,50)+r2(5,54)*three &
&                  -r1(1,38)*three)*xx+(+r6(6,10)+r6(6,11)*two+r6(6,12)-r5(3,18)-r5(3,19)*two &
&                  -r5(3,20)+r4(6,27)*three+r4(6,28)*six+r4(6,29)*three-r3(3,24)*three &
&                  -r3(3,25)*six-r3(3,26)*three)*xz+(+r6(3,11)-r5(1,19)+r4(3,28)*three &
&                  +r4(3,40)-r3(1,25)*three-r3(1,50)+r2(5,54)*three-r1(1,38)*three)*zz+( &
&                  +r5(6,23)+r5(6,24)-r4(3,36)-r4(3,37)+r3(6,41)*three+r3(6,42)*three &
&                  -r2(5,35)*three-r2(5,36)*three)*xxz+(+r5(3,23)+r5(3,24)-r4(1,36)-r4(1,37) &
&                  +r3(3,41)*three+r3(3,42)*three-r2(1,35)*three-r2(1,36)*three)*xzz+rxyz(1) &
&                  *xxzz
      eri(6,5,1,5)=r112+rxyz(19)*qx+(+r7(9,3)+r7(9,4)-r6(5,7)-r6(5,8)+r5(9,11)*three &
&                  +r5(9,12)*three-r4(5,12)*three-r4(5,13)*three)*qz+(+r6(9,11)+r6(9,12) &
&                  -r5(5,19)-r5(5,20)+r4(9,28)*three+r4(9,29)*three-r3(5,25)*three &
&                  -r3(5,26)*three)*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,1,5)=r211+(+r7(9,3)*two-r6(5,7)*two+r5(9,11)*six-r4(5,12)*six)*qx+rxyz(16) &
&                  *qz+(+r6(9,10)-r5(5,18)+r4(9,27)*three-r3(5,24)*three)*xx+(+r6(5,11)*two &
&                  -r5(2,19)*two+r4(5,28)*six-r3(2,25)*six)*xz+rxyz(3)*xxz
      eri(2,6,1,5)=r031+rxyz(8)*qz
      eri(3,6,1,5)=r013+(+r7(14,3)*two+r7(14,4)-r6(9,7)*two-r6(9,8)+r5(14,11)*six &
&                  +r5(14,12)*three+r5(5,21)*two+r5(5,22)-r4(9,12)*six-r4(9,13)*three &
&                  -r4(2,34)*two-r4(2,35)+r3(5,39)*six+r3(5,40)*three-r2(4,33)*six &
&                  -r2(4,34)*three)*qz+(+r6(9,10)+r6(9,11)*two-r5(5,18)-r5(5,19)*two &
&                  +r4(9,27)*three+r4(9,28)*six-r3(5,24)*three-r3(5,25)*six)*zz+rxyz(3)*zzz
      eri(4,6,1,5)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,1,5)=r112+rxyz(18)*qx+(+r7(9,3)+r7(9,4)-r6(5,7)-r6(5,8)+r5(9,11)*three &
&                  +r5(9,12)*three-r4(5,12)*three-r4(5,13)*three)*qz+(+r6(9,10)+r6(9,11) &
&                  -r5(5,18)-r5(5,19)+r4(9,27)*three+r4(9,28)*three-r3(5,24)*three &
&                  -r3(5,25)*three)*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,1,5)=r022+(+r7(13,3)+r7(13,4)-r6(8,7)-r6(8,8)+r5(13,11)*three &
&                  +r5(13,12)*three+r5(6,21)+r5(6,22)-r4(8,12)*three-r4(8,13)*three-r4(3,34) &
&                  -r4(3,35)+r3(6,39)*three+r3(6,40)*three-r2(5,33)*three-r2(5,34)*three)*qz &
&                  +rxyz(6)*zz
!
      r400= r8(8)-r7(4,2)+r6(3,4)+r6(8,9)*six-r5(1,4)-r5(4,17)*six+r4(3,26)*six &
&          +r4(8,38)*three-r3(1,23)*six-r3(4,48)*three+r2(5,52)*three-r1(1,36)*three
      r310= r8(12)-r7(7,2)+r6(5,4)+r6(12,9)*three-r5(2,4)-r5(7,17)*three+r4(5,26)*three &
&          -r3(2,23)*three
      r301= r8(13)-r7(8,2)+r6(6,4)+r6(13,9)*three-r5(3,4)-r5(8,17)*three+r4(6,26)*three &
&          -r3(3,23)*three
      r220= r8(17)-r7(11,2)+r6(8,4)+r6(8,9)+r6(17,9)-r5(4,4)-r5(4,17)-r5(11,17)+r4(3,26) &
&          +r4(8,26)+r4(8,38)-r3(1,23)-r3(4,23)-r3(4,48)+r2(5,52)-r1(1,36)
      r211= r8(18)-r7(12,2)+r6(9,4)+r6(18,9)-r5(5,4)-r5(12,17)+r4(9,26)-r3(5,23)
      r202= r8(19)-r7(13,2)+r6(10,4)+r6(8,9)+r6(19,9)-r5(6,4)-r5(4,17)-r5(13,17)+r4(3,26) &
&          +r4(10,26)+r4(8,38)-r3(1,23)-r3(6,23)-r3(4,48)+r2(5,52)-r1(1,36)
      r130= r8(23)-r7(16,2)+r6(12,4)+r6(12,9)*three-r5(7,4)-r5(7,17)*three+r4(5,26)*three &
&          -r3(2,23)*three
      r121= r8(24)-r7(17,2)+r6(13,4)+r6(13,9)-r5(8,4)-r5(8,17)+r4(6,26)-r3(3,23)
      r112= r8(25)-r7(18,2)+r6(14,4)+r6(12,9)-r5(9,4)-r5(7,17)+r4(5,26)-r3(2,23)
      r103= r8(26)-r7(19,2)+r6(15,4)+r6(13,9)*three-r5(10,4)-r5(8,17)*three+r4(6,26)*three &
&          -r3(3,23)*three
      r040= r8(30)-r7(22,2)+r6(17,4)+r6(17,9)*six-r5(11,4)-r5(11,17)*six+r4(8,26)*six &
&          +r4(8,38)*three-r3(4,23)*six-r3(4,48)*three+r2(5,52)*three-r1(1,36)*three
      r031= r8(31)-r7(23,2)+r6(18,4)+r6(18,9)*three-r5(12,4)-r5(12,17)*three+r4(9,26)*three &
&          -r3(5,23)*three
      r022= r8(32)-r7(24,2)+r6(19,4)+r6(17,9)+r6(19,9)-r5(13,4)-r5(11,17)-r5(13,17) &
&          +r4(8,26)+r4(10,26)+r4(8,38)-r3(4,23)-r3(6,23)-r3(4,48)+r2(5,52)-r1(1,36)
      r013= r8(33)-r7(25,2)+r6(20,4)+r6(18,9)*three-r5(14,4)-r5(12,17)*three+r4(9,26)*three &
&          -r3(5,23)*three
      r004= r8(34)-r7(26,2)+r6(21,4)+r6(19,9)*six-r5(15,4)-r5(13,17)*six+r4(10,26)*six &
&          +r4(8,38)*three-r3(6,23)*six-r3(4,48)*three+r2(5,52)*three-r1(1,36)*three
      rxyz(1)=+r4(8,42)-r3(4,52)+r2(5,56)-r1(1,40)
      rxyz(2)=+r5(12,24)-r4(7,37)+r3(5,42)-r2(4,36)
      rxyz(3)=+r5(12,23)-r4(7,36)+r3(5,41)-r2(4,35)
      rxyz(4)=+r6(17,10)-r5(11,18)+r4(8,27)+r4(8,39)-r3(4,24)-r3(4,49)+r2(5,53)-r1(1,37)
      rxyz(5)=+r6(17,12)-r5(11,20)+r4(8,29)+r4(8,41)-r3(4,26)-r3(4,51)+r2(5,55)-r1(1,39)
      rxyz(6)=+r6(17,11)-r5(11,19)+r4(8,28)+r4(8,40)-r3(4,25)-r3(4,50)+r2(5,54)-r1(1,38)
      rxyz(7)=+r7(23,3)-r6(16,7)+r5(12,11)+r5(12,21)*three-r4(7,12)-r4(7,34)*three &
&             +r3(5,39)*three-r2(4,33)*three
      rxyz(8)=+r7(23,4)-r6(16,8)+r5(12,12)+r5(12,22)*three-r4(7,13)-r4(7,35)*three &
&             +r3(5,40)*three-r2(4,34)*three
      rxyz(9)=+r7(18,3)+r7(18,4)-r6(12,7)-r6(12,8)+r5(9,11)+r5(9,12)-r4(5,12)-r4(5,13)
      rxyz(10)=+r6(18,11)-r5(12,19)+r4(9,28)-r3(5,25)
      rxyz(11)=+r6(12,11)-r5(7,19)+r4(5,28)-r3(2,25)
      rxyz(12)=+r7(24,3)-r6(17,7)+r5(13,11)+r5(13,21)-r4(8,12)-r4(8,34)+r3(6,39)-r2(5,33)
      rxyz(13)=+r7(24,4)-r6(17,8)+r5(13,12)+r5(13,22)-r4(8,13)-r4(8,35)+r3(6,40)-r2(5,34)
      rxyz(14)=+r7(17,4)-r6(11,8)+r5(8,12)+r5(8,22)-r4(4,13)-r4(4,35)+r3(3,40)-r2(1,34)
      rxyz(15)=+r7(17,3)-r6(11,7)+r5(8,11)+r5(8,21)-r4(4,12)-r4(4,34)+r3(3,39)-r2(1,33)
      rxyz(16)=+r7(12,4)-r6(7,8)+r5(5,12)+r5(12,22)-r4(2,13)-r4(7,35)+r3(5,40)-r2(4,34)
      rxyz(17)=+r7(12,3)-r6(7,7)+r5(5,11)+r5(12,21)-r4(2,12)-r4(7,34)+r3(5,39)-r2(4,33)
      rxyz(18)=+r7(25,3)-r6(18,7)+r5(14,11)+r5(12,21)-r4(9,12)-r4(7,34)+r3(5,39)-r2(4,33)
      rxyz(19)=+r7(25,4)-r6(18,8)+r5(14,12)+r5(12,22)-r4(9,13)-r4(7,35)+r3(5,40)-r2(4,34)
      rxyz(20)=+r6(13,11)*four-r5(8,19)*four+r4(6,28)*four-r3(3,25)*four
      eri(1,1,2,5)=r400+(+r7(8,3)*two+r7(8,4)*two-r6(4,7)*two-r6(4,8)*two+r5(3,11)*two &
&                  +r5(3,12)*two+r5(8,21)*six+r5(8,22)*six-r4(1,12)*two-r4(1,13)*two &
&                  -r4(4,34)*six-r4(4,35)*six+r3(3,39)*six+r3(3,40)*six-r2(1,33)*six &
&                  -r2(1,34)*six)*qx+(+r6(8,10)+r6(8,11)*four+r6(8,12)-r5(4,18)-r5(4,19)*four &
&                  -r5(4,20)+r4(3,27)+r4(3,28)*four+r4(3,29)+r4(8,39)+r4(8,40)*four+r4(8,41) &
&                  -r3(1,24)-r3(1,25)*four-r3(1,26)-r3(4,49)-r3(4,50)*four-r3(4,51)+r2(5,53) &
&                  +r2(5,54)*four+r2(5,55)-r1(1,37)-r1(1,38)*four-r1(1,39))*xx+(+r5(8,23)*two &
&                  +r5(8,24)*two-r4(4,36)*two-r4(4,37)*two+r3(3,41)*two+r3(3,42)*two &
&                  -r2(1,35)*two-r2(1,36)*two)*xxx+rxyz(1)*xxxx
      eri(2,1,2,5)=r220+(+r7(17,4)*two-r6(11,8)*two+r5(8,12)*two+r5(8,22)*two &
&                  -r4(4,13)*two-r4(4,35)*two+r3(3,40)*two-r2(1,34)*two)*qx+rxyz(5)*xx
      eri(3,1,2,5)=r202+(+r7(19,4)*two-r6(13,8)*two+r5(10,12)*two+r5(8,22)*two &
&                  -r4(6,13)*two-r4(4,35)*two+r3(3,40)*two-r2(1,34)*two)*qx+(+r7(13,3)*two &
&                  -r6(8,7)*two+r5(6,11)*two+r5(13,21)*two-r4(3,12)*two-r4(8,34)*two &
&                  +r3(6,39)*two-r2(5,33)*two)*qz+(+r6(19,12)-r5(13,20)+r4(10,29)+r4(8,41) &
&                  -r3(6,26)-r3(4,51)+r2(5,55)-r1(1,39))*xx+rxyz(20)*xz+(+r6(8,10)-r5(4,18) &
&                  +r4(3,27)+r4(8,39)-r3(1,24)-r3(4,49)+r2(5,53)-r1(1,37))*zz+(+r5(13,24)*two &
&                  -r4(8,37)*two+r3(6,42)*two-r2(5,36)*two)*xxz+(+r5(8,23)*two-r4(4,36)*two &
&                  +r3(3,41)*two-r2(1,35)*two)*xzz+rxyz(1)*xxzz
      eri(4,1,2,5)=r310+(+r7(12,3)+r7(12,4)*two-r6(7,7)-r6(7,8)*two+r5(5,11)+r5(5,12)*two &
&                  +r5(12,21)+r5(12,22)*two-r4(2,12)-r4(2,13)*two-r4(7,34)-r4(7,35)*two &
&                  +r3(5,39)+r3(5,40)*two-r2(4,33)-r2(4,34)*two)*qx+(+r6(12,11)*two+r6(12,12) &
&                  -r5(7,19)*two-r5(7,20)+r4(5,28)*two+r4(5,29)-r3(2,25)*two-r3(2,26))*xx &
&                  +rxyz(2)*xxx
      eri(5,1,2,5)=r301+(+r7(13,3)+r7(13,4)*two-r6(8,7)-r6(8,8)*two+r5(6,11)+r5(6,12)*two &
&                  +r5(13,21)+r5(13,22)*two-r4(3,12)-r4(3,13)*two-r4(8,34)-r4(8,35)*two &
&                  +r3(6,39)+r3(6,40)*two-r2(5,33)-r2(5,34)*two)*qx+(+r7(8,3)-r6(4,7) &
&                  +r5(3,11)+r5(8,21)*three-r4(1,12)-r4(4,34)*three+r3(3,39)*three &
&                  -r2(1,33)*three)*qz+(+r6(13,11)*two+r6(13,12)-r5(8,19)*two-r5(8,20) &
&                  +r4(6,28)*two+r4(6,29)-r3(3,25)*two-r3(3,26))*xx+(+r6(8,10)+r6(8,11)*two &
&                  -r5(4,18)-r5(4,19)*two+r4(3,27)+r4(3,28)*two+r4(8,39)+r4(8,40)*two &
&                  -r3(1,24)-r3(1,25)*two-r3(4,49)-r3(4,50)*two+r2(5,53)+r2(5,54)*two &
&                  -r1(1,37)-r1(1,38)*two)*xz+(+r5(13,24)-r4(8,37)+r3(6,42)-r2(5,36))*xxx+( &
&                  +r5(8,23)*two+r5(8,24)-r4(4,36)*two-r4(4,37)+r3(3,41)*two+r3(3,42) &
&                  -r2(1,35)*two-r2(1,36))*xxz+rxyz(1)*xxxz
      eri(6,1,2,5)=r211+(+r7(18,4)*two-r6(12,8)*two+r5(9,12)*two-r4(5,13)*two)*qx &
&                  +rxyz(17)*qz+(+r6(18,12)-r5(12,20)+r4(9,29)-r3(5,26))*xx+(+r6(12,11)*two &
&                  -r5(7,19)*two+r4(5,28)*two-r3(2,25)*two)*xz+rxyz(2)*xxz
      eri(1,2,2,5)=r220+(+r7(17,3)*two-r6(11,7)*two+r5(8,11)*two+r5(8,21)*two &
&                  -r4(4,12)*two-r4(4,34)*two+r3(3,39)*two-r2(1,33)*two)*qx+rxyz(4)*xx
      eri(2,2,2,5)=r040
      eri(3,2,2,5)=r022+(+r7(24,3)*two-r6(17,7)*two+r5(13,11)*two+r5(13,21)*two &
&                  -r4(8,12)*two-r4(8,34)*two+r3(6,39)*two-r2(5,33)*two)*qz+rxyz(4)*zz
      eri(4,2,2,5)=r130+rxyz(7)*qx
      eri(5,2,2,5)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,2,5)=r031+rxyz(7)*qz
      eri(1,3,2,5)=r202+(+r7(19,3)*two-r6(13,7)*two+r5(10,11)*two+r5(8,21)*two &
&                  -r4(6,12)*two-r4(4,34)*two+r3(3,39)*two-r2(1,33)*two)*qx+(+r7(13,4)*two &
&                  -r6(8,8)*two+r5(6,12)*two+r5(13,22)*two-r4(3,13)*two-r4(8,35)*two &
&                  +r3(6,40)*two-r2(5,34)*two)*qz+(+r6(19,10)-r5(13,18)+r4(10,27)+r4(8,39) &
&                  -r3(6,24)-r3(4,49)+r2(5,53)-r1(1,37))*xx+rxyz(20)*xz+(+r6(8,12)-r5(4,20) &
&                  +r4(3,29)+r4(8,41)-r3(1,26)-r3(4,51)+r2(5,55)-r1(1,39))*zz+(+r5(13,23)*two &
&                  -r4(8,36)*two+r3(6,41)*two-r2(5,35)*two)*xxz+(+r5(8,24)*two-r4(4,37)*two &
&                  +r3(3,42)*two-r2(1,36)*two)*xzz+rxyz(1)*xxzz
      eri(2,3,2,5)=r022+(+r7(24,4)*two-r6(17,8)*two+r5(13,12)*two+r5(13,22)*two &
&                  -r4(8,13)*two-r4(8,35)*two+r3(6,40)*two-r2(5,34)*two)*qz+rxyz(5)*zz
      eri(3,3,2,5)=r004+(+r7(26,3)*two+r7(26,4)*two-r6(19,7)*two-r6(19,8)*two &
&                  +r5(15,11)*two+r5(15,12)*two+r5(13,21)*six+r5(13,22)*six-r4(10,12)*two &
&                  -r4(10,13)*two-r4(8,34)*six-r4(8,35)*six+r3(6,39)*six+r3(6,40)*six &
&                  -r2(5,33)*six-r2(5,34)*six)*qz+(+r6(19,10)+r6(19,11)*four+r6(19,12) &
&                  -r5(13,18)-r5(13,19)*four-r5(13,20)+r4(10,27)+r4(10,28)*four+r4(10,29) &
&                  +r4(8,39)+r4(8,40)*four+r4(8,41)-r3(6,24)-r3(6,25)*four-r3(6,26)-r3(4,49) &
&                  -r3(4,50)*four-r3(4,51)+r2(5,53)+r2(5,54)*four+r2(5,55)-r1(1,37) &
&                  -r1(1,38)*four-r1(1,39))*zz+(+r5(13,23)*two+r5(13,24)*two-r4(8,36)*two &
&                  -r4(8,37)*two+r3(6,41)*two+r3(6,42)*two-r2(5,35)*two-r2(5,36)*two)*zzz &
&                  +rxyz(1)*zzzz
      eri(4,3,2,5)=r112+rxyz(18)*qx+(+r7(18,4)*two-r6(12,8)*two+r5(9,12)*two-r4(5,13)*two &
&                  )*qz+(+r6(18,11)*two-r5(12,19)*two+r4(9,28)*two-r3(5,25)*two)*xz+( &
&                  +r6(12,12)-r5(7,20)+r4(5,29)-r3(2,26))*zz+rxyz(2)*xzz
      eri(5,3,2,5)=r103+(+r7(26,3)-r6(19,7)+r5(15,11)+r5(13,21)*three-r4(10,12) &
&                  -r4(8,34)*three+r3(6,39)*three-r2(5,33)*three)*qx+(+r7(19,3)+r7(19,4)*two &
&                  -r6(13,7)-r6(13,8)*two+r5(10,11)+r5(10,12)*two+r5(8,21)+r5(8,22)*two &
&                  -r4(6,12)-r4(6,13)*two-r4(4,34)-r4(4,35)*two+r3(3,39)+r3(3,40)*two &
&                  -r2(1,33)-r2(1,34)*two)*qz+(+r6(19,10)+r6(19,11)*two-r5(13,18) &
&                  -r5(13,19)*two+r4(10,27)+r4(10,28)*two+r4(8,39)+r4(8,40)*two-r3(6,24) &
&                  -r3(6,25)*two-r3(4,49)-r3(4,50)*two+r2(5,53)+r2(5,54)*two-r1(1,37) &
&                  -r1(1,38)*two)*xz+(+r6(13,11)*two+r6(13,12)-r5(8,19)*two-r5(8,20) &
&                  +r4(6,28)*two+r4(6,29)-r3(3,25)*two-r3(3,26))*zz+(+r5(13,23)*two+r5(13,24) &
&                  -r4(8,36)*two-r4(8,37)+r3(6,41)*two+r3(6,42)-r2(5,35)*two-r2(5,36))*xzz+( &
&                  +r5(8,24)-r4(4,37)+r3(3,42)-r2(1,36))*zzz+rxyz(1)*xzzz
      eri(6,3,2,5)=r013+(+r7(25,3)+r7(25,4)*two-r6(18,7)-r6(18,8)*two+r5(14,11) &
&                  +r5(14,12)*two+r5(12,21)+r5(12,22)*two-r4(9,12)-r4(9,13)*two-r4(7,34) &
&                  -r4(7,35)*two+r3(5,39)+r3(5,40)*two-r2(4,33)-r2(4,34)*two)*qz+( &
&                  +r6(18,11)*two+r6(18,12)-r5(12,19)*two-r5(12,20)+r4(9,28)*two+r4(9,29) &
&                  -r3(5,25)*two-r3(5,26))*zz+rxyz(2)*zzz
      eri(1,4,2,5)=r310+(+r7(12,3)*two+r7(12,4)-r6(7,7)*two-r6(7,8)+r5(5,11)*two+r5(5,12) &
&                  +r5(12,21)*two+r5(12,22)-r4(2,12)*two-r4(2,13)-r4(7,34)*two-r4(7,35) &
&                  +r3(5,39)*two+r3(5,40)-r2(4,33)*two-r2(4,34))*qx+(+r6(12,10)+r6(12,11)*two &
&                  -r5(7,18)-r5(7,19)*two+r4(5,27)+r4(5,28)*two-r3(2,24)-r3(2,25)*two)*xx &
&                  +rxyz(3)*xxx
      eri(2,4,2,5)=r130+rxyz(8)*qx
      eri(3,4,2,5)=r112+rxyz(19)*qx+(+r7(18,3)*two-r6(12,7)*two+r5(9,11)*two-r4(5,12)*two &
&                  )*qz+(+r6(18,11)*two-r5(12,19)*two+r4(9,28)*two-r3(5,25)*two)*xz+( &
&                  +r6(12,10)-r5(7,18)+r4(5,27)-r3(2,24))*zz+rxyz(3)*xzz
      eri(4,4,2,5)=r220+(+r7(17,3)+r7(17,4)-r6(11,7)-r6(11,8)+r5(8,11)+r5(8,12)+r5(8,21) &
&                  +r5(8,22)-r4(4,12)-r4(4,13)-r4(4,34)-r4(4,35)+r3(3,39)+r3(3,40)-r2(1,33) &
&                  -r2(1,34))*qx+rxyz(6)*xx
      eri(5,4,2,5)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(12,10)+r6(12,11)-r5(7,18) &
&                  -r5(7,19)+r4(5,27)+r4(5,28)-r3(2,24)-r3(2,25))*xz+rxyz(3)*xxz
      eri(6,4,2,5)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,2,5)=r301+(+r7(13,3)*two+r7(13,4)-r6(8,7)*two-r6(8,8)+r5(6,11)*two+r5(6,12) &
&                  +r5(13,21)*two+r5(13,22)-r4(3,12)*two-r4(3,13)-r4(8,34)*two-r4(8,35) &
&                  +r3(6,39)*two+r3(6,40)-r2(5,33)*two-r2(5,34))*qx+(+r7(8,4)-r6(4,8) &
&                  +r5(3,12)+r5(8,22)*three-r4(1,13)-r4(4,35)*three+r3(3,40)*three &
&                  -r2(1,34)*three)*qz+(+r6(13,10)+r6(13,11)*two-r5(8,18)-r5(8,19)*two &
&                  +r4(6,27)+r4(6,28)*two-r3(3,24)-r3(3,25)*two)*xx+(+r6(8,11)*two+r6(8,12) &
&                  -r5(4,19)*two-r5(4,20)+r4(3,28)*two+r4(3,29)+r4(8,40)*two+r4(8,41) &
&                  -r3(1,25)*two-r3(1,26)-r3(4,50)*two-r3(4,51)+r2(5,54)*two+r2(5,55) &
&                  -r1(1,38)*two-r1(1,39))*xz+(+r5(13,23)-r4(8,36)+r3(6,41)-r2(5,35))*xxx+( &
&                  +r5(8,23)+r5(8,24)*two-r4(4,36)-r4(4,37)*two+r3(3,41)+r3(3,42)*two &
&                  -r2(1,35)-r2(1,36)*two)*xxz+rxyz(1)*xxxz
      eri(2,5,2,5)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,2,5)=r103+(+r7(26,4)-r6(19,8)+r5(15,12)+r5(13,22)*three-r4(10,13) &
&                  -r4(8,35)*three+r3(6,40)*three-r2(5,34)*three)*qx+(+r7(19,3)*two+r7(19,4) &
&                  -r6(13,7)*two-r6(13,8)+r5(10,11)*two+r5(10,12)+r5(8,21)*two+r5(8,22) &
&                  -r4(6,12)*two-r4(6,13)-r4(4,34)*two-r4(4,35)+r3(3,39)*two+r3(3,40) &
&                  -r2(1,33)*two-r2(1,34))*qz+(+r6(19,11)*two+r6(19,12)-r5(13,19)*two &
&                  -r5(13,20)+r4(10,28)*two+r4(10,29)+r4(8,40)*two+r4(8,41)-r3(6,25)*two &
&                  -r3(6,26)-r3(4,50)*two-r3(4,51)+r2(5,54)*two+r2(5,55)-r1(1,38)*two &
&                  -r1(1,39))*xz+(+r6(13,10)+r6(13,11)*two-r5(8,18)-r5(8,19)*two+r4(6,27) &
&                  +r4(6,28)*two-r3(3,24)-r3(3,25)*two)*zz+(+r5(13,23)+r5(13,24)*two-r4(8,36) &
&                  -r4(8,37)*two+r3(6,41)+r3(6,42)*two-r2(5,35)-r2(5,36)*two)*xzz+(+r5(8,23) &
&                  -r4(4,36)+r3(3,41)-r2(1,35))*zzz+rxyz(1)*xzzz
      eri(4,5,2,5)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(12,11)+r6(12,12)-r5(7,19) &
&                  -r5(7,20)+r4(5,28)+r4(5,29)-r3(2,25)-r3(2,26))*xz+rxyz(2)*xxz
      eri(5,5,2,5)=r202+(+r7(19,3)+r7(19,4)-r6(13,7)-r6(13,8)+r5(10,11)+r5(10,12) &
&                  +r5(8,21)+r5(8,22)-r4(6,12)-r4(6,13)-r4(4,34)-r4(4,35)+r3(3,39)+r3(3,40) &
&                  -r2(1,33)-r2(1,34))*qx+(+r7(13,3)+r7(13,4)-r6(8,7)-r6(8,8)+r5(6,11) &
&                  +r5(6,12)+r5(13,21)+r5(13,22)-r4(3,12)-r4(3,13)-r4(8,34)-r4(8,35)+r3(6,39) &
&                  +r3(6,40)-r2(5,33)-r2(5,34))*qz+(+r6(19,11)-r5(13,19)+r4(10,28)+r4(8,40) &
&                  -r3(6,25)-r3(4,50)+r2(5,54)-r1(1,38))*xx+(+r6(13,10)+r6(13,11)*two &
&                  +r6(13,12)-r5(8,18)-r5(8,19)*two-r5(8,20)+r4(6,27)+r4(6,28)*two+r4(6,29) &
&                  -r3(3,24)-r3(3,25)*two-r3(3,26))*xz+(+r6(8,11)-r5(4,19)+r4(3,28)+r4(8,40) &
&                  -r3(1,25)-r3(4,50)+r2(5,54)-r1(1,38))*zz+(+r5(13,23)+r5(13,24)-r4(8,36) &
&                  -r4(8,37)+r3(6,41)+r3(6,42)-r2(5,35)-r2(5,36))*xxz+(+r5(8,23)+r5(8,24) &
&                  -r4(4,36)-r4(4,37)+r3(3,41)+r3(3,42)-r2(1,35)-r2(1,36))*xzz+rxyz(1)*xxzz
      eri(6,5,2,5)=r112+rxyz(19)*qx+(+r7(18,3)+r7(18,4)-r6(12,7)-r6(12,8)+r5(9,11) &
&                  +r5(9,12)-r4(5,12)-r4(5,13))*qz+(+r6(18,11)+r6(18,12)-r5(12,19)-r5(12,20) &
&                  +r4(9,28)+r4(9,29)-r3(5,25)-r3(5,26))*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,2,5)=r211+(+r7(18,3)*two-r6(12,7)*two+r5(9,11)*two-r4(5,12)*two)*qx &
&                  +rxyz(16)*qz+(+r6(18,10)-r5(12,18)+r4(9,27)-r3(5,24))*xx+(+r6(12,11)*two &
&                  -r5(7,19)*two+r4(5,28)*two-r3(2,25)*two)*xz+rxyz(3)*xxz
      eri(2,6,2,5)=r031+rxyz(8)*qz
      eri(3,6,2,5)=r013+(+r7(25,3)*two+r7(25,4)-r6(18,7)*two-r6(18,8)+r5(14,11)*two &
&                  +r5(14,12)+r5(12,21)*two+r5(12,22)-r4(9,12)*two-r4(9,13)-r4(7,34)*two &
&                  -r4(7,35)+r3(5,39)*two+r3(5,40)-r2(4,33)*two-r2(4,34))*qz+(+r6(18,10) &
&                  +r6(18,11)*two-r5(12,18)-r5(12,19)*two+r4(9,27)+r4(9,28)*two-r3(5,24) &
&                  -r3(5,25)*two)*zz+rxyz(3)*zzz
      eri(4,6,2,5)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,2,5)=r112+rxyz(18)*qx+(+r7(18,3)+r7(18,4)-r6(12,7)-r6(12,8)+r5(9,11) &
&                  +r5(9,12)-r4(5,12)-r4(5,13))*qz+(+r6(18,10)+r6(18,11)-r5(12,18)-r5(12,19) &
&                  +r4(9,27)+r4(9,28)-r3(5,24)-r3(5,25))*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,2,5)=r022+(+r7(24,3)+r7(24,4)-r6(17,7)-r6(17,8)+r5(13,11)+r5(13,12) &
&                  +r5(13,21)+r5(13,22)-r4(8,12)-r4(8,13)-r4(8,34)-r4(8,35)+r3(6,39)+r3(6,40) &
&                  -r2(5,33)-r2(5,34))*qz+rxyz(6)*zz
!
      r400= r8(10)-r7(6,1)*two-r7(6,2)+r6(3,1)+r6(3,2)*two+r6(3,4)*three+r6(10,9)*six &
&          -r5(1,1)-r5(1,3)*two-r5(1,4)-r5(6,13)*p12-r5(6,17)*six+r4(3,14)*six &
&          +r4(3,18)*p12+r4(3,26)*p18+r4(10,38)*three-r3(1,11)*six-r3(1,19)*p12 &
&          -r3(1,23)*six-r3(6,43)*six-r3(6,48)*three+r2(5,37)*three+r2(5,42)*six &
&          +r2(5,52)*nine-r1(1,21)*three-r1(1,31)*six-r1(1,36)*three
      r310= r8(14)-r7(9,1)*two-r7(9,2)+r6(5,1)+r6(5,2)*two+r6(5,4)*three+r6(14,9)*three &
&          -r5(2,1)-r5(2,3)*two-r5(2,4)-r5(9,13)*six-r5(9,17)*three+r4(5,14)*three &
&          +r4(5,18)*six+r4(5,26)*nine-r3(2,11)*three-r3(2,19)*six-r3(2,23)*three
      r301= r8(15)-r7(10,1)*two-r7(10,2)+r6(6,1)+r6(6,2)*two+r6(6,4)*three+r6(15,9)*three &
&          -r5(3,1)-r5(3,3)*two-r5(3,4)-r5(10,13)*six-r5(10,17)*three+r4(6,14)*three &
&          +r4(6,18)*six+r4(6,26)*nine-r3(3,11)*three-r3(3,19)*six-r3(3,23)*three
      r220= r8(19)-r7(13,1)*two-r7(13,2)+r6(8,1)+r6(8,2)*two+r6(8,4)*three+r6(10,9) &
&          +r6(19,9)-r5(4,1)-r5(4,3)*two-r5(4,4)-r5(6,13)*two-r5(13,13)*two-r5(6,17) &
&          -r5(13,17)+r4(3,14)+r4(8,14)+r4(3,18)*two+r4(8,18)*two+r4(3,26)*three &
&          +r4(8,26)*three+r4(10,38)-r3(1,11)-r3(4,11)-r3(1,19)*two-r3(4,19)*two-r3(1,23) &
&          -r3(4,23)-r3(6,43)*two-r3(6,48)+r2(5,37)+r2(5,42)*two+r2(5,52)*three-r1(1,21) &
&          -r1(1,31)*two-r1(1,36)
      r211= r8(20)-r7(14,1)*two-r7(14,2)+r6(9,1)+r6(9,2)*two+r6(9,4)*three+r6(20,9)-r5(5,1) &
&          -r5(5,3)*two-r5(5,4)-r5(14,13)*two-r5(14,17)+r4(9,14)+r4(9,18)*two &
&          +r4(9,26)*three-r3(5,11)-r3(5,19)*two-r3(5,23)
      r202= r8(21)-r7(15,1)*two-r7(15,2)+r6(10,1)+r6(10,2)*two+r6(10,4)*three+r6(10,9) &
&          +r6(21,9)-r5(6,1)-r5(6,3)*two-r5(6,4)-r5(6,13)*two-r5(15,13)*two-r5(6,17) &
&          -r5(15,17)+r4(3,14)+r4(10,14)+r4(3,18)*two+r4(10,18)*two+r4(3,26)*three &
&          +r4(10,26)*three+r4(10,38)-r3(1,11)-r3(6,11)-r3(1,19)*two-r3(6,19)*two-r3(1,23) &
&          -r3(6,23)-r3(6,43)*two-r3(6,48)+r2(5,37)+r2(5,42)*two+r2(5,52)*three-r1(1,21) &
&          -r1(1,31)*two-r1(1,36)
      r130= r8(25)-r7(18,1)*two-r7(18,2)+r6(12,1)+r6(12,2)*two+r6(12,4)*three &
&          +r6(14,9)*three-r5(7,1)-r5(7,3)*two-r5(7,4)-r5(9,13)*six-r5(9,17)*three &
&          +r4(5,14)*three+r4(5,18)*six+r4(5,26)*nine-r3(2,11)*three-r3(2,19)*six &
&          -r3(2,23)*three
      r121= r8(26)-r7(19,1)*two-r7(19,2)+r6(13,1)+r6(13,2)*two+r6(13,4)*three+r6(15,9) &
&          -r5(8,1)-r5(8,3)*two-r5(8,4)-r5(10,13)*two-r5(10,17)+r4(6,14)+r4(6,18)*two &
&          +r4(6,26)*three-r3(3,11)-r3(3,19)*two-r3(3,23)
      r112= r8(27)-r7(20,1)*two-r7(20,2)+r6(14,1)+r6(14,2)*two+r6(14,4)*three+r6(14,9) &
&          -r5(9,1)-r5(9,3)*two-r5(9,4)-r5(9,13)*two-r5(9,17)+r4(5,14)+r4(5,18)*two &
&          +r4(5,26)*three-r3(2,11)-r3(2,19)*two-r3(2,23)
      r103= r8(28)-r7(21,1)*two-r7(21,2)+r6(15,1)+r6(15,2)*two+r6(15,4)*three &
&          +r6(15,9)*three-r5(10,1)-r5(10,3)*two-r5(10,4)-r5(10,13)*six-r5(10,17)*three &
&          +r4(6,14)*three+r4(6,18)*six+r4(6,26)*nine-r3(3,11)*three-r3(3,19)*six &
&          -r3(3,23)*three
      r040= r8(32)-r7(24,1)*two-r7(24,2)+r6(17,1)+r6(17,2)*two+r6(17,4)*three+r6(19,9)*six &
&          -r5(11,1)-r5(11,3)*two-r5(11,4)-r5(13,13)*p12-r5(13,17)*six+r4(8,14)*six &
&          +r4(8,18)*p12+r4(8,26)*p18+r4(10,38)*three-r3(4,11)*six-r3(4,19)*p12 &
&          -r3(4,23)*six-r3(6,43)*six-r3(6,48)*three+r2(5,37)*three+r2(5,42)*six &
&          +r2(5,52)*nine-r1(1,21)*three-r1(1,31)*six-r1(1,36)*three
      r031= r8(33)-r7(25,1)*two-r7(25,2)+r6(18,1)+r6(18,2)*two+r6(18,4)*three &
&          +r6(20,9)*three-r5(12,1)-r5(12,3)*two-r5(12,4)-r5(14,13)*six-r5(14,17)*three &
&          +r4(9,14)*three+r4(9,18)*six+r4(9,26)*nine-r3(5,11)*three-r3(5,19)*six &
&          -r3(5,23)*three
      r022= r8(34)-r7(26,1)*two-r7(26,2)+r6(19,1)+r6(19,2)*two+r6(19,4)*three+r6(19,9) &
&          +r6(21,9)-r5(13,1)-r5(13,3)*two-r5(13,4)-r5(13,13)*two-r5(15,13)*two-r5(13,17) &
&          -r5(15,17)+r4(8,14)+r4(10,14)+r4(8,18)*two+r4(10,18)*two+r4(8,26)*three &
&          +r4(10,26)*three+r4(10,38)-r3(4,11)-r3(6,11)-r3(4,19)*two-r3(6,19)*two-r3(4,23) &
&          -r3(6,23)-r3(6,43)*two-r3(6,48)+r2(5,37)+r2(5,42)*two+r2(5,52)*three-r1(1,21) &
&          -r1(1,31)*two-r1(1,36)
      r013= r8(35)-r7(27,1)*two-r7(27,2)+r6(20,1)+r6(20,2)*two+r6(20,4)*three &
&          +r6(20,9)*three-r5(14,1)-r5(14,3)*two-r5(14,4)-r5(14,13)*six-r5(14,17)*three &
&          +r4(9,14)*three+r4(9,18)*six+r4(9,26)*nine-r3(5,11)*three-r3(5,19)*six &
&          -r3(5,23)*three
      r004= r8(36)-r7(28,1)*two-r7(28,2)+r6(21,1)+r6(21,2)*two+r6(21,4)*three+r6(21,9)*six &
&          -r5(15,1)-r5(15,3)*two-r5(15,4)-r5(15,13)*p12-r5(15,17)*six+r4(10,14)*six &
&          +r4(10,18)*p12+r4(10,26)*p18+r4(10,38)*three-r3(6,11)*six-r3(6,19)*p12 &
&          -r3(6,23)*six-r3(6,43)*six-r3(6,48)*three+r2(5,37)*three+r2(5,42)*six &
&          +r2(5,52)*nine-r1(1,21)*three-r1(1,31)*six-r1(1,36)*three
      rxyz(1)=+r4(10,42)-r3(6,47)*two-r3(6,52)+r2(5,41)+r2(5,46)*two+r2(5,56)*three &
&             -r1(1,25)-r1(1,35)*two-r1(1,40)
      rxyz(2)=+r5(14,24)-r4(9,33)*two-r4(9,37)+r3(5,30)+r3(5,34)*two+r3(5,42)*three &
&             -r2(4,24)-r2(4,32)*two-r2(4,36)
      rxyz(3)=+r5(14,23)-r4(9,32)*two-r4(9,36)+r3(5,29)+r3(5,33)*two+r3(5,41)*three &
&             -r2(4,23)-r2(4,31)*two-r2(4,35)
      rxyz(4)=+r6(19,10)-r5(13,14)*two-r5(13,18)+r4(8,15)+r4(8,19)*two+r4(8,27)*three &
&             +r4(10,39)-r3(4,12)-r3(4,20)*two-r3(4,24)-r3(6,44)*two-r3(6,49)+r2(5,38) &
&             +r2(5,43)*two+r2(5,53)*three-r1(1,22)-r1(1,32)*two-r1(1,37)
      rxyz(5)=+r6(19,12)-r5(13,16)*two-r5(13,20)+r4(8,17)+r4(8,21)*two+r4(8,29)*three &
&             +r4(10,41)-r3(4,14)-r3(4,22)*two-r3(4,26)-r3(6,46)*two-r3(6,51)+r2(5,40) &
&             +r2(5,45)*two+r2(5,55)*three-r1(1,24)-r1(1,34)*two-r1(1,39)
      rxyz(6)=+r6(19,11)-r5(13,15)*two-r5(13,19)+r4(8,16)+r4(8,20)*two+r4(8,28)*three &
&             +r4(10,40)-r3(4,13)-r3(4,21)*two-r3(4,25)-r3(6,45)*two-r3(6,50)+r2(5,39) &
&             +r2(5,44)*two+r2(5,54)*three-r1(1,23)-r1(1,33)*two-r1(1,38)
      rxyz(7)=+r7(25,3)-r6(18,5)*two-r6(18,7)+r5(12,5)+r5(12,7)*two+r5(12,11)*three &
&             +r5(14,21)*three-r4(7,6)-r4(7,10)*two-r4(7,12)-r4(9,30)*six-r4(9,34)*three &
&             +r3(5,27)*three+r3(5,31)*six+r3(5,39)*nine-r2(4,21)*three-r2(4,29)*six &
&             -r2(4,33)*three
      rxyz(8)=+r7(25,4)-r6(18,6)*two-r6(18,8)+r5(12,6)+r5(12,8)*two+r5(12,12)*three &
&             +r5(14,22)*three-r4(7,7)-r4(7,11)*two-r4(7,13)-r4(9,31)*six-r4(9,35)*three &
&             +r3(5,28)*three+r3(5,32)*six+r3(5,40)*nine-r2(4,22)*three-r2(4,30)*six &
&             -r2(4,34)*three
      rxyz(9)=+r7(20,3)+r7(20,4)-r6(14,5)*two-r6(14,6)*two-r6(14,7)-r6(14,8)+r5(9,5) &
&             +r5(9,6)+r5(9,7)*two+r5(9,8)*two+r5(9,11)*three+r5(9,12)*three-r4(5,6) &
&             -r4(5,7)-r4(5,10)*two-r4(5,11)*two-r4(5,12)-r4(5,13)
      rxyz(10)=+r6(20,11)-r5(14,15)*two-r5(14,19)+r4(9,16)+r4(9,20)*two+r4(9,28)*three &
&             -r3(5,13)-r3(5,21)*two-r3(5,25)
      rxyz(11)=+r6(14,11)-r5(9,15)*two-r5(9,19)+r4(5,16)+r4(5,20)*two+r4(5,28)*three &
&             -r3(2,13)-r3(2,21)*two-r3(2,25)
      rxyz(12)=+r7(26,3)-r6(19,5)*two-r6(19,7)+r5(13,5)+r5(13,7)*two+r5(13,11)*three &
&             +r5(15,21)-r4(8,6)-r4(8,10)*two-r4(8,12)-r4(10,30)*two-r4(10,34)+r3(6,27) &
&             +r3(6,31)*two+r3(6,39)*three-r2(5,21)-r2(5,29)*two-r2(5,33)
      rxyz(13)=+r7(26,4)-r6(19,6)*two-r6(19,8)+r5(13,6)+r5(13,8)*two+r5(13,12)*three &
&             +r5(15,22)-r4(8,7)-r4(8,11)*two-r4(8,13)-r4(10,31)*two-r4(10,35)+r3(6,28) &
&             +r3(6,32)*two+r3(6,40)*three-r2(5,22)-r2(5,30)*two-r2(5,34)
      rxyz(14)=+r7(19,4)-r6(13,6)*two-r6(13,8)+r5(8,6)+r5(8,8)*two+r5(8,12)*three &
&             +r5(10,22)-r4(4,7)-r4(4,11)*two-r4(4,13)-r4(6,31)*two-r4(6,35)+r3(3,28) &
&             +r3(3,32)*two+r3(3,40)*three-r2(1,22)-r2(1,30)*two-r2(1,34)
      rxyz(15)=+r7(19,3)-r6(13,5)*two-r6(13,7)+r5(8,5)+r5(8,7)*two+r5(8,11)*three &
&             +r5(10,21)-r4(4,6)-r4(4,10)*two-r4(4,12)-r4(6,30)*two-r4(6,34)+r3(3,27) &
&             +r3(3,31)*two+r3(3,39)*three-r2(1,21)-r2(1,29)*two-r2(1,33)
      rxyz(16)=+r7(14,4)-r6(9,6)*two-r6(9,8)+r5(5,6)+r5(5,8)*two+r5(5,12)*three+r5(14,22) &
&             -r4(2,7)-r4(2,11)*two-r4(2,13)-r4(9,31)*two-r4(9,35)+r3(5,28)+r3(5,32)*two &
&             +r3(5,40)*three-r2(4,22)-r2(4,30)*two-r2(4,34)
      rxyz(17)=+r7(14,3)-r6(9,5)*two-r6(9,7)+r5(5,5)+r5(5,7)*two+r5(5,11)*three+r5(14,21) &
&             -r4(2,6)-r4(2,10)*two-r4(2,12)-r4(9,30)*two-r4(9,34)+r3(5,27)+r3(5,31)*two &
&             +r3(5,39)*three-r2(4,21)-r2(4,29)*two-r2(4,33)
      rxyz(18)=+r7(27,3)-r6(20,5)*two-r6(20,7)+r5(14,5)+r5(14,7)*two+r5(14,11)*three &
&             +r5(14,21)-r4(9,6)-r4(9,10)*two-r4(9,12)-r4(9,30)*two-r4(9,34)+r3(5,27) &
&             +r3(5,31)*two+r3(5,39)*three-r2(4,21)-r2(4,29)*two-r2(4,33)
      rxyz(19)=+r7(27,4)-r6(20,6)*two-r6(20,8)+r5(14,6)+r5(14,8)*two+r5(14,12)*three &
&             +r5(14,22)-r4(9,7)-r4(9,11)*two-r4(9,13)-r4(9,31)*two-r4(9,35)+r3(5,28) &
&             +r3(5,32)*two+r3(5,40)*three-r2(4,22)-r2(4,30)*two-r2(4,34)
      rxyz(20)=+r6(15,11)*four-r5(10,15)*eight-r5(10,19)*four+r4(6,16)*four+r4(6,20)*eight &
&             +r4(6,28)*p12-r3(3,13)*four-r3(3,21)*eight-r3(3,25)*four
      eri(1,1,3,5)=r400+(+r7(10,3)*two+r7(10,4)*two-r6(6,5)*four-r6(6,6)*four-r6(6,7)*two &
&                  -r6(6,8)*two+r5(3,5)*two+r5(3,6)*two+r5(3,7)*four+r5(3,8)*four &
&                  +r5(3,11)*six+r5(3,12)*six+r5(10,21)*six+r5(10,22)*six-r4(1,6)*two &
&                  -r4(1,7)*two-r4(1,10)*four-r4(1,11)*four-r4(1,12)*two-r4(1,13)*two &
&                  -r4(6,30)*p12-r4(6,31)*p12-r4(6,34)*six-r4(6,35)*six+r3(3,27)*six &
&                  +r3(3,28)*six+r3(3,31)*p12+r3(3,32)*p12+r3(3,39)*p18+r3(3,40)*p18 &
&                  -r2(1,21)*six-r2(1,22)*six-r2(1,29)*p12-r2(1,30)*p12-r2(1,33)*six &
&                  -r2(1,34)*six)*qx+(+r6(10,10)+r6(10,11)*four+r6(10,12)-r5(6,14)*two &
&                  -r5(6,15)*eight-r5(6,16)*two-r5(6,18)-r5(6,19)*four-r5(6,20)+r4(3,15) &
&                  +r4(3,16)*four+r4(3,17)+r4(3,19)*two+r4(3,20)*eight+r4(3,21)*two &
&                  +r4(3,27)*three+r4(3,28)*p12+r4(3,29)*three+r4(10,39)+r4(10,40)*four &
&                  +r4(10,41)-r3(1,12)-r3(1,13)*four-r3(1,14)-r3(1,20)*two-r3(1,21)*eight &
&                  -r3(1,22)*two-r3(1,24)-r3(1,25)*four-r3(1,26)-r3(6,44)*two-r3(6,45)*eight &
&                  -r3(6,46)*two-r3(6,49)-r3(6,50)*four-r3(6,51)+r2(5,38)+r2(5,39)*four &
&                  +r2(5,40)+r2(5,43)*two+r2(5,44)*eight+r2(5,45)*two+r2(5,53)*three &
&                  +r2(5,54)*p12+r2(5,55)*three-r1(1,22)-r1(1,23)*four-r1(1,24)-r1(1,32)*two &
&                  -r1(1,33)*eight-r1(1,34)*two-r1(1,37)-r1(1,38)*four-r1(1,39))*xx+( &
&                  +r5(10,23)*two+r5(10,24)*two-r4(6,32)*four-r4(6,33)*four-r4(6,36)*two &
&                  -r4(6,37)*two+r3(3,29)*two+r3(3,30)*two+r3(3,33)*four+r3(3,34)*four &
&                  +r3(3,41)*six+r3(3,42)*six-r2(1,23)*two-r2(1,24)*two-r2(1,31)*four &
&                  -r2(1,32)*four-r2(1,35)*two-r2(1,36)*two)*xxx+rxyz(1)*xxxx
      eri(2,1,3,5)=r220+(+r7(19,4)*two-r6(13,6)*four-r6(13,8)*two+r5(8,6)*two &
&                  +r5(8,8)*four+r5(8,12)*six+r5(10,22)*two-r4(4,7)*two-r4(4,11)*four &
&                  -r4(4,13)*two-r4(6,31)*four-r4(6,35)*two+r3(3,28)*two+r3(3,32)*four &
&                  +r3(3,40)*six-r2(1,22)*two-r2(1,30)*four-r2(1,34)*two)*qx+rxyz(5)*xx
      eri(3,1,3,5)=r202+(+r7(21,4)*two-r6(15,6)*four-r6(15,8)*two+r5(10,6)*two &
&                  +r5(10,8)*four+r5(10,12)*six+r5(10,22)*two-r4(6,7)*two-r4(6,11)*four &
&                  -r4(6,13)*two-r4(6,31)*four-r4(6,35)*two+r3(3,28)*two+r3(3,32)*four &
&                  +r3(3,40)*six-r2(1,22)*two-r2(1,30)*four-r2(1,34)*two)*qx+(+r7(15,3)*two &
&                  -r6(10,5)*four-r6(10,7)*two+r5(6,5)*two+r5(6,7)*four+r5(6,11)*six &
&                  +r5(15,21)*two-r4(3,6)*two-r4(3,10)*four-r4(3,12)*two-r4(10,30)*four &
&                  -r4(10,34)*two+r3(6,27)*two+r3(6,31)*four+r3(6,39)*six-r2(5,21)*two &
&                  -r2(5,29)*four-r2(5,33)*two)*qz+(+r6(21,12)-r5(15,16)*two-r5(15,20) &
&                  +r4(10,17)+r4(10,21)*two+r4(10,29)*three+r4(10,41)-r3(6,14)-r3(6,22)*two &
&                  -r3(6,26)-r3(6,46)*two-r3(6,51)+r2(5,40)+r2(5,45)*two+r2(5,55)*three &
&                  -r1(1,24)-r1(1,34)*two-r1(1,39))*xx+rxyz(20)*xz+(+r6(10,10)-r5(6,14)*two &
&                  -r5(6,18)+r4(3,15)+r4(3,19)*two+r4(3,27)*three+r4(10,39)-r3(1,12) &
&                  -r3(1,20)*two-r3(1,24)-r3(6,44)*two-r3(6,49)+r2(5,38)+r2(5,43)*two &
&                  +r2(5,53)*three-r1(1,22)-r1(1,32)*two-r1(1,37))*zz+(+r5(15,24)*two &
&                  -r4(10,33)*four-r4(10,37)*two+r3(6,30)*two+r3(6,34)*four+r3(6,42)*six &
&                  -r2(5,24)*two-r2(5,32)*four-r2(5,36)*two)*xxz+(+r5(10,23)*two &
&                  -r4(6,32)*four-r4(6,36)*two+r3(3,29)*two+r3(3,33)*four+r3(3,41)*six &
&                  -r2(1,23)*two-r2(1,31)*four-r2(1,35)*two)*xzz+rxyz(1)*xxzz
      eri(4,1,3,5)=r310+(+r7(14,3)+r7(14,4)*two-r6(9,5)*two-r6(9,6)*four-r6(9,7) &
&                  -r6(9,8)*two+r5(5,5)+r5(5,6)*two+r5(5,7)*two+r5(5,8)*four+r5(5,11)*three &
&                  +r5(5,12)*six+r5(14,21)+r5(14,22)*two-r4(2,6)-r4(2,7)*two-r4(2,10)*two &
&                  -r4(2,11)*four-r4(2,12)-r4(2,13)*two-r4(9,30)*two-r4(9,31)*four-r4(9,34) &
&                  -r4(9,35)*two+r3(5,27)+r3(5,28)*two+r3(5,31)*two+r3(5,32)*four &
&                  +r3(5,39)*three+r3(5,40)*six-r2(4,21)-r2(4,22)*two-r2(4,29)*two &
&                  -r2(4,30)*four-r2(4,33)-r2(4,34)*two)*qx+(+r6(14,11)*two+r6(14,12) &
&                  -r5(9,15)*four-r5(9,16)*two-r5(9,19)*two-r5(9,20)+r4(5,16)*two+r4(5,17) &
&                  +r4(5,20)*four+r4(5,21)*two+r4(5,28)*six+r4(5,29)*three-r3(2,13)*two &
&                  -r3(2,14)-r3(2,21)*four-r3(2,22)*two-r3(2,25)*two-r3(2,26))*xx+rxyz(2)*xxx
      eri(5,1,3,5)=r301+(+r7(15,3)+r7(15,4)*two-r6(10,5)*two-r6(10,6)*four-r6(10,7) &
&                  -r6(10,8)*two+r5(6,5)+r5(6,6)*two+r5(6,7)*two+r5(6,8)*four+r5(6,11)*three &
&                  +r5(6,12)*six+r5(15,21)+r5(15,22)*two-r4(3,6)-r4(3,7)*two-r4(3,10)*two &
&                  -r4(3,11)*four-r4(3,12)-r4(3,13)*two-r4(10,30)*two-r4(10,31)*four &
&                  -r4(10,34)-r4(10,35)*two+r3(6,27)+r3(6,28)*two+r3(6,31)*two+r3(6,32)*four &
&                  +r3(6,39)*three+r3(6,40)*six-r2(5,21)-r2(5,22)*two-r2(5,29)*two &
&                  -r2(5,30)*four-r2(5,33)-r2(5,34)*two)*qx+(+r7(10,3)-r6(6,5)*two-r6(6,7) &
&                  +r5(3,5)+r5(3,7)*two+r5(3,11)*three+r5(10,21)*three-r4(1,6)-r4(1,10)*two &
&                  -r4(1,12)-r4(6,30)*six-r4(6,34)*three+r3(3,27)*three+r3(3,31)*six &
&                  +r3(3,39)*nine-r2(1,21)*three-r2(1,29)*six-r2(1,33)*three)*qz+( &
&                  +r6(15,11)*two+r6(15,12)-r5(10,15)*four-r5(10,16)*two-r5(10,19)*two &
&                  -r5(10,20)+r4(6,16)*two+r4(6,17)+r4(6,20)*four+r4(6,21)*two+r4(6,28)*six &
&                  +r4(6,29)*three-r3(3,13)*two-r3(3,14)-r3(3,21)*four-r3(3,22)*two &
&                  -r3(3,25)*two-r3(3,26))*xx+(+r6(10,10)+r6(10,11)*two-r5(6,14)*two &
&                  -r5(6,15)*four-r5(6,18)-r5(6,19)*two+r4(3,15)+r4(3,16)*two+r4(3,19)*two &
&                  +r4(3,20)*four+r4(3,27)*three+r4(3,28)*six+r4(10,39)+r4(10,40)*two &
&                  -r3(1,12)-r3(1,13)*two-r3(1,20)*two-r3(1,21)*four-r3(1,24)-r3(1,25)*two &
&                  -r3(6,44)*two-r3(6,45)*four-r3(6,49)-r3(6,50)*two+r2(5,38)+r2(5,39)*two &
&                  +r2(5,43)*two+r2(5,44)*four+r2(5,53)*three+r2(5,54)*six-r1(1,22) &
&                  -r1(1,23)*two-r1(1,32)*two-r1(1,33)*four-r1(1,37)-r1(1,38)*two)*xz+( &
&                  +r5(15,24)-r4(10,33)*two-r4(10,37)+r3(6,30)+r3(6,34)*two+r3(6,42)*three &
&                  -r2(5,24)-r2(5,32)*two-r2(5,36))*xxx+(+r5(10,23)*two+r5(10,24) &
&                  -r4(6,32)*four-r4(6,33)*two-r4(6,36)*two-r4(6,37)+r3(3,29)*two+r3(3,30) &
&                  +r3(3,33)*four+r3(3,34)*two+r3(3,41)*six+r3(3,42)*three-r2(1,23)*two &
&                  -r2(1,24)-r2(1,31)*four-r2(1,32)*two-r2(1,35)*two-r2(1,36))*xxz+rxyz(1) &
&                  *xxxz
      eri(6,1,3,5)=r211+(+r7(20,4)*two-r6(14,6)*four-r6(14,8)*two+r5(9,6)*two &
&                  +r5(9,8)*four+r5(9,12)*six-r4(5,7)*two-r4(5,11)*four-r4(5,13)*two)*qx &
&                  +rxyz(17)*qz+(+r6(20,12)-r5(14,16)*two-r5(14,20)+r4(9,17)+r4(9,21)*two &
&                  +r4(9,29)*three-r3(5,14)-r3(5,22)*two-r3(5,26))*xx+(+r6(14,11)*two &
&                  -r5(9,15)*four-r5(9,19)*two+r4(5,16)*two+r4(5,20)*four+r4(5,28)*six &
&                  -r3(2,13)*two-r3(2,21)*four-r3(2,25)*two)*xz+rxyz(2)*xxz
      eri(1,2,3,5)=r220+(+r7(19,3)*two-r6(13,5)*four-r6(13,7)*two+r5(8,5)*two &
&                  +r5(8,7)*four+r5(8,11)*six+r5(10,21)*two-r4(4,6)*two-r4(4,10)*four &
&                  -r4(4,12)*two-r4(6,30)*four-r4(6,34)*two+r3(3,27)*two+r3(3,31)*four &
&                  +r3(3,39)*six-r2(1,21)*two-r2(1,29)*four-r2(1,33)*two)*qx+rxyz(4)*xx
      eri(2,2,3,5)=r040
      eri(3,2,3,5)=r022+(+r7(26,3)*two-r6(19,5)*four-r6(19,7)*two+r5(13,5)*two &
&                  +r5(13,7)*four+r5(13,11)*six+r5(15,21)*two-r4(8,6)*two-r4(8,10)*four &
&                  -r4(8,12)*two-r4(10,30)*four-r4(10,34)*two+r3(6,27)*two+r3(6,31)*four &
&                  +r3(6,39)*six-r2(5,21)*two-r2(5,29)*four-r2(5,33)*two)*qz+rxyz(4)*zz
      eri(4,2,3,5)=r130+rxyz(7)*qx
      eri(5,2,3,5)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,3,5)=r031+rxyz(7)*qz
      eri(1,3,3,5)=r202+(+r7(21,3)*two-r6(15,5)*four-r6(15,7)*two+r5(10,5)*two &
&                  +r5(10,7)*four+r5(10,11)*six+r5(10,21)*two-r4(6,6)*two-r4(6,10)*four &
&                  -r4(6,12)*two-r4(6,30)*four-r4(6,34)*two+r3(3,27)*two+r3(3,31)*four &
&                  +r3(3,39)*six-r2(1,21)*two-r2(1,29)*four-r2(1,33)*two)*qx+(+r7(15,4)*two &
&                  -r6(10,6)*four-r6(10,8)*two+r5(6,6)*two+r5(6,8)*four+r5(6,12)*six &
&                  +r5(15,22)*two-r4(3,7)*two-r4(3,11)*four-r4(3,13)*two-r4(10,31)*four &
&                  -r4(10,35)*two+r3(6,28)*two+r3(6,32)*four+r3(6,40)*six-r2(5,22)*two &
&                  -r2(5,30)*four-r2(5,34)*two)*qz+(+r6(21,10)-r5(15,14)*two-r5(15,18) &
&                  +r4(10,15)+r4(10,19)*two+r4(10,27)*three+r4(10,39)-r3(6,12)-r3(6,20)*two &
&                  -r3(6,24)-r3(6,44)*two-r3(6,49)+r2(5,38)+r2(5,43)*two+r2(5,53)*three &
&                  -r1(1,22)-r1(1,32)*two-r1(1,37))*xx+rxyz(20)*xz+(+r6(10,12)-r5(6,16)*two &
&                  -r5(6,20)+r4(3,17)+r4(3,21)*two+r4(3,29)*three+r4(10,41)-r3(1,14) &
&                  -r3(1,22)*two-r3(1,26)-r3(6,46)*two-r3(6,51)+r2(5,40)+r2(5,45)*two &
&                  +r2(5,55)*three-r1(1,24)-r1(1,34)*two-r1(1,39))*zz+(+r5(15,23)*two &
&                  -r4(10,32)*four-r4(10,36)*two+r3(6,29)*two+r3(6,33)*four+r3(6,41)*six &
&                  -r2(5,23)*two-r2(5,31)*four-r2(5,35)*two)*xxz+(+r5(10,24)*two &
&                  -r4(6,33)*four-r4(6,37)*two+r3(3,30)*two+r3(3,34)*four+r3(3,42)*six &
&                  -r2(1,24)*two-r2(1,32)*four-r2(1,36)*two)*xzz+rxyz(1)*xxzz
      eri(2,3,3,5)=r022+(+r7(26,4)*two-r6(19,6)*four-r6(19,8)*two+r5(13,6)*two &
&                  +r5(13,8)*four+r5(13,12)*six+r5(15,22)*two-r4(8,7)*two-r4(8,11)*four &
&                  -r4(8,13)*two-r4(10,31)*four-r4(10,35)*two+r3(6,28)*two+r3(6,32)*four &
&                  +r3(6,40)*six-r2(5,22)*two-r2(5,30)*four-r2(5,34)*two)*qz+rxyz(5)*zz
      eri(3,3,3,5)=r004+(+r7(28,3)*two+r7(28,4)*two-r6(21,5)*four-r6(21,6)*four &
&                  -r6(21,7)*two-r6(21,8)*two+r5(15,5)*two+r5(15,6)*two+r5(15,7)*four &
&                  +r5(15,8)*four+r5(15,11)*six+r5(15,12)*six+r5(15,21)*six+r5(15,22)*six &
&                  -r4(10,6)*two-r4(10,7)*two-r4(10,10)*four-r4(10,11)*four-r4(10,12)*two &
&                  -r4(10,13)*two-r4(10,30)*p12-r4(10,31)*p12-r4(10,34)*six-r4(10,35)*six &
&                  +r3(6,27)*six+r3(6,28)*six+r3(6,31)*p12+r3(6,32)*p12+r3(6,39)*p18 &
&                  +r3(6,40)*p18-r2(5,21)*six-r2(5,22)*six-r2(5,29)*p12-r2(5,30)*p12 &
&                  -r2(5,33)*six-r2(5,34)*six)*qz+(+r6(21,10)+r6(21,11)*four+r6(21,12) &
&                  -r5(15,14)*two-r5(15,15)*eight-r5(15,16)*two-r5(15,18)-r5(15,19)*four &
&                  -r5(15,20)+r4(10,15)+r4(10,16)*four+r4(10,17)+r4(10,19)*two &
&                  +r4(10,20)*eight+r4(10,21)*two+r4(10,27)*three+r4(10,28)*p12 &
&                  +r4(10,29)*three+r4(10,39)+r4(10,40)*four+r4(10,41)-r3(6,12)-r3(6,13)*four &
&                  -r3(6,14)-r3(6,20)*two-r3(6,21)*eight-r3(6,22)*two-r3(6,24)-r3(6,25)*four &
&                  -r3(6,26)-r3(6,44)*two-r3(6,45)*eight-r3(6,46)*two-r3(6,49)-r3(6,50)*four &
&                  -r3(6,51)+r2(5,38)+r2(5,39)*four+r2(5,40)+r2(5,43)*two+r2(5,44)*eight &
&                  +r2(5,45)*two+r2(5,53)*three+r2(5,54)*p12+r2(5,55)*three-r1(1,22) &
&                  -r1(1,23)*four-r1(1,24)-r1(1,32)*two-r1(1,33)*eight-r1(1,34)*two-r1(1,37) &
&                  -r1(1,38)*four-r1(1,39))*zz+(+r5(15,23)*two+r5(15,24)*two-r4(10,32)*four &
&                  -r4(10,33)*four-r4(10,36)*two-r4(10,37)*two+r3(6,29)*two+r3(6,30)*two &
&                  +r3(6,33)*four+r3(6,34)*four+r3(6,41)*six+r3(6,42)*six-r2(5,23)*two &
&                  -r2(5,24)*two-r2(5,31)*four-r2(5,32)*four-r2(5,35)*two-r2(5,36)*two)*zzz &
&                  +rxyz(1)*zzzz
      eri(4,3,3,5)=r112+rxyz(18)*qx+(+r7(20,4)*two-r6(14,6)*four-r6(14,8)*two+r5(9,6)*two &
&                  +r5(9,8)*four+r5(9,12)*six-r4(5,7)*two-r4(5,11)*four-r4(5,13)*two)*qz+( &
&                  +r6(20,11)*two-r5(14,15)*four-r5(14,19)*two+r4(9,16)*two+r4(9,20)*four &
&                  +r4(9,28)*six-r3(5,13)*two-r3(5,21)*four-r3(5,25)*two)*xz+(+r6(14,12) &
&                  -r5(9,16)*two-r5(9,20)+r4(5,17)+r4(5,21)*two+r4(5,29)*three-r3(2,14) &
&                  -r3(2,22)*two-r3(2,26))*zz+rxyz(2)*xzz
      eri(5,3,3,5)=r103+(+r7(28,3)-r6(21,5)*two-r6(21,7)+r5(15,5)+r5(15,7)*two &
&                  +r5(15,11)*three+r5(15,21)*three-r4(10,6)-r4(10,10)*two-r4(10,12) &
&                  -r4(10,30)*six-r4(10,34)*three+r3(6,27)*three+r3(6,31)*six+r3(6,39)*nine &
&                  -r2(5,21)*three-r2(5,29)*six-r2(5,33)*three)*qx+(+r7(21,3)+r7(21,4)*two &
&                  -r6(15,5)*two-r6(15,6)*four-r6(15,7)-r6(15,8)*two+r5(10,5)+r5(10,6)*two &
&                  +r5(10,7)*two+r5(10,8)*four+r5(10,11)*three+r5(10,12)*six+r5(10,21) &
&                  +r5(10,22)*two-r4(6,6)-r4(6,7)*two-r4(6,10)*two-r4(6,11)*four-r4(6,12) &
&                  -r4(6,13)*two-r4(6,30)*two-r4(6,31)*four-r4(6,34)-r4(6,35)*two+r3(3,27) &
&                  +r3(3,28)*two+r3(3,31)*two+r3(3,32)*four+r3(3,39)*three+r3(3,40)*six &
&                  -r2(1,21)-r2(1,22)*two-r2(1,29)*two-r2(1,30)*four-r2(1,33)-r2(1,34)*two) &
&                  *qz+(+r6(21,10)+r6(21,11)*two-r5(15,14)*two-r5(15,15)*four-r5(15,18) &
&                  -r5(15,19)*two+r4(10,15)+r4(10,16)*two+r4(10,19)*two+r4(10,20)*four &
&                  +r4(10,27)*three+r4(10,28)*six+r4(10,39)+r4(10,40)*two-r3(6,12) &
&                  -r3(6,13)*two-r3(6,20)*two-r3(6,21)*four-r3(6,24)-r3(6,25)*two &
&                  -r3(6,44)*two-r3(6,45)*four-r3(6,49)-r3(6,50)*two+r2(5,38)+r2(5,39)*two &
&                  +r2(5,43)*two+r2(5,44)*four+r2(5,53)*three+r2(5,54)*six-r1(1,22) &
&                  -r1(1,23)*two-r1(1,32)*two-r1(1,33)*four-r1(1,37)-r1(1,38)*two)*xz+( &
&                  +r6(15,11)*two+r6(15,12)-r5(10,15)*four-r5(10,16)*two-r5(10,19)*two &
&                  -r5(10,20)+r4(6,16)*two+r4(6,17)+r4(6,20)*four+r4(6,21)*two+r4(6,28)*six &
&                  +r4(6,29)*three-r3(3,13)*two-r3(3,14)-r3(3,21)*four-r3(3,22)*two &
&                  -r3(3,25)*two-r3(3,26))*zz+(+r5(15,23)*two+r5(15,24)-r4(10,32)*four &
&                  -r4(10,33)*two-r4(10,36)*two-r4(10,37)+r3(6,29)*two+r3(6,30)+r3(6,33)*four &
&                  +r3(6,34)*two+r3(6,41)*six+r3(6,42)*three-r2(5,23)*two-r2(5,24) &
&                  -r2(5,31)*four-r2(5,32)*two-r2(5,35)*two-r2(5,36))*xzz+(+r5(10,24) &
&                  -r4(6,33)*two-r4(6,37)+r3(3,30)+r3(3,34)*two+r3(3,42)*three-r2(1,24) &
&                  -r2(1,32)*two-r2(1,36))*zzz+rxyz(1)*xzzz
      eri(6,3,3,5)=r013+(+r7(27,3)+r7(27,4)*two-r6(20,5)*two-r6(20,6)*four-r6(20,7) &
&                  -r6(20,8)*two+r5(14,5)+r5(14,6)*two+r5(14,7)*two+r5(14,8)*four &
&                  +r5(14,11)*three+r5(14,12)*six+r5(14,21)+r5(14,22)*two-r4(9,6)-r4(9,7)*two &
&                  -r4(9,10)*two-r4(9,11)*four-r4(9,12)-r4(9,13)*two-r4(9,30)*two &
&                  -r4(9,31)*four-r4(9,34)-r4(9,35)*two+r3(5,27)+r3(5,28)*two+r3(5,31)*two &
&                  +r3(5,32)*four+r3(5,39)*three+r3(5,40)*six-r2(4,21)-r2(4,22)*two &
&                  -r2(4,29)*two-r2(4,30)*four-r2(4,33)-r2(4,34)*two)*qz+(+r6(20,11)*two &
&                  +r6(20,12)-r5(14,15)*four-r5(14,16)*two-r5(14,19)*two-r5(14,20) &
&                  +r4(9,16)*two+r4(9,17)+r4(9,20)*four+r4(9,21)*two+r4(9,28)*six &
&                  +r4(9,29)*three-r3(5,13)*two-r3(5,14)-r3(5,21)*four-r3(5,22)*two &
&                  -r3(5,25)*two-r3(5,26))*zz+rxyz(2)*zzz
      eri(1,4,3,5)=r310+(+r7(14,3)*two+r7(14,4)-r6(9,5)*four-r6(9,6)*two-r6(9,7)*two &
&                  -r6(9,8)+r5(5,5)*two+r5(5,6)+r5(5,7)*four+r5(5,8)*two+r5(5,11)*six &
&                  +r5(5,12)*three+r5(14,21)*two+r5(14,22)-r4(2,6)*two-r4(2,7)-r4(2,10)*four &
&                  -r4(2,11)*two-r4(2,12)*two-r4(2,13)-r4(9,30)*four-r4(9,31)*two &
&                  -r4(9,34)*two-r4(9,35)+r3(5,27)*two+r3(5,28)+r3(5,31)*four+r3(5,32)*two &
&                  +r3(5,39)*six+r3(5,40)*three-r2(4,21)*two-r2(4,22)-r2(4,29)*four &
&                  -r2(4,30)*two-r2(4,33)*two-r2(4,34))*qx+(+r6(14,10)+r6(14,11)*two &
&                  -r5(9,14)*two-r5(9,15)*four-r5(9,18)-r5(9,19)*two+r4(5,15)+r4(5,16)*two &
&                  +r4(5,19)*two+r4(5,20)*four+r4(5,27)*three+r4(5,28)*six-r3(2,12) &
&                  -r3(2,13)*two-r3(2,20)*two-r3(2,21)*four-r3(2,24)-r3(2,25)*two)*xx+rxyz(3) &
&                  *xxx
      eri(2,4,3,5)=r130+rxyz(8)*qx
      eri(3,4,3,5)=r112+rxyz(19)*qx+(+r7(20,3)*two-r6(14,5)*four-r6(14,7)*two+r5(9,5)*two &
&                  +r5(9,7)*four+r5(9,11)*six-r4(5,6)*two-r4(5,10)*four-r4(5,12)*two)*qz+( &
&                  +r6(20,11)*two-r5(14,15)*four-r5(14,19)*two+r4(9,16)*two+r4(9,20)*four &
&                  +r4(9,28)*six-r3(5,13)*two-r3(5,21)*four-r3(5,25)*two)*xz+(+r6(14,10) &
&                  -r5(9,14)*two-r5(9,18)+r4(5,15)+r4(5,19)*two+r4(5,27)*three-r3(2,12) &
&                  -r3(2,20)*two-r3(2,24))*zz+rxyz(3)*xzz
      eri(4,4,3,5)=r220+(+r7(19,3)+r7(19,4)-r6(13,5)*two-r6(13,6)*two-r6(13,7)-r6(13,8) &
&                  +r5(8,5)+r5(8,6)+r5(8,7)*two+r5(8,8)*two+r5(8,11)*three+r5(8,12)*three &
&                  +r5(10,21)+r5(10,22)-r4(4,6)-r4(4,7)-r4(4,10)*two-r4(4,11)*two-r4(4,12) &
&                  -r4(4,13)-r4(6,30)*two-r4(6,31)*two-r4(6,34)-r4(6,35)+r3(3,27)+r3(3,28) &
&                  +r3(3,31)*two+r3(3,32)*two+r3(3,39)*three+r3(3,40)*three-r2(1,21)-r2(1,22) &
&                  -r2(1,29)*two-r2(1,30)*two-r2(1,33)-r2(1,34))*qx+rxyz(6)*xx
      eri(5,4,3,5)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(14,10)+r6(14,11) &
&                  -r5(9,14)*two-r5(9,15)*two-r5(9,18)-r5(9,19)+r4(5,15)+r4(5,16) &
&                  +r4(5,19)*two+r4(5,20)*two+r4(5,27)*three+r4(5,28)*three-r3(2,12)-r3(2,13) &
&                  -r3(2,20)*two-r3(2,21)*two-r3(2,24)-r3(2,25))*xz+rxyz(3)*xxz
      eri(6,4,3,5)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,3,5)=r301+(+r7(15,3)*two+r7(15,4)-r6(10,5)*four-r6(10,6)*two-r6(10,7)*two &
&                  -r6(10,8)+r5(6,5)*two+r5(6,6)+r5(6,7)*four+r5(6,8)*two+r5(6,11)*six &
&                  +r5(6,12)*three+r5(15,21)*two+r5(15,22)-r4(3,6)*two-r4(3,7)-r4(3,10)*four &
&                  -r4(3,11)*two-r4(3,12)*two-r4(3,13)-r4(10,30)*four-r4(10,31)*two &
&                  -r4(10,34)*two-r4(10,35)+r3(6,27)*two+r3(6,28)+r3(6,31)*four+r3(6,32)*two &
&                  +r3(6,39)*six+r3(6,40)*three-r2(5,21)*two-r2(5,22)-r2(5,29)*four &
&                  -r2(5,30)*two-r2(5,33)*two-r2(5,34))*qx+(+r7(10,4)-r6(6,6)*two-r6(6,8) &
&                  +r5(3,6)+r5(3,8)*two+r5(3,12)*three+r5(10,22)*three-r4(1,7)-r4(1,11)*two &
&                  -r4(1,13)-r4(6,31)*six-r4(6,35)*three+r3(3,28)*three+r3(3,32)*six &
&                  +r3(3,40)*nine-r2(1,22)*three-r2(1,30)*six-r2(1,34)*three)*qz+(+r6(15,10) &
&                  +r6(15,11)*two-r5(10,14)*two-r5(10,15)*four-r5(10,18)-r5(10,19)*two &
&                  +r4(6,15)+r4(6,16)*two+r4(6,19)*two+r4(6,20)*four+r4(6,27)*three &
&                  +r4(6,28)*six-r3(3,12)-r3(3,13)*two-r3(3,20)*two-r3(3,21)*four-r3(3,24) &
&                  -r3(3,25)*two)*xx+(+r6(10,11)*two+r6(10,12)-r5(6,15)*four-r5(6,16)*two &
&                  -r5(6,19)*two-r5(6,20)+r4(3,16)*two+r4(3,17)+r4(3,20)*four+r4(3,21)*two &
&                  +r4(3,28)*six+r4(3,29)*three+r4(10,40)*two+r4(10,41)-r3(1,13)*two-r3(1,14) &
&                  -r3(1,21)*four-r3(1,22)*two-r3(1,25)*two-r3(1,26)-r3(6,45)*four &
&                  -r3(6,46)*two-r3(6,50)*two-r3(6,51)+r2(5,39)*two+r2(5,40)+r2(5,44)*four &
&                  +r2(5,45)*two+r2(5,54)*six+r2(5,55)*three-r1(1,23)*two-r1(1,24) &
&                  -r1(1,33)*four-r1(1,34)*two-r1(1,38)*two-r1(1,39))*xz+(+r5(15,23) &
&                  -r4(10,32)*two-r4(10,36)+r3(6,29)+r3(6,33)*two+r3(6,41)*three-r2(5,23) &
&                  -r2(5,31)*two-r2(5,35))*xxx+(+r5(10,23)+r5(10,24)*two-r4(6,32)*two &
&                  -r4(6,33)*four-r4(6,36)-r4(6,37)*two+r3(3,29)+r3(3,30)*two+r3(3,33)*two &
&                  +r3(3,34)*four+r3(3,41)*three+r3(3,42)*six-r2(1,23)-r2(1,24)*two &
&                  -r2(1,31)*two-r2(1,32)*four-r2(1,35)-r2(1,36)*two)*xxz+rxyz(1)*xxxz
      eri(2,5,3,5)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,3,5)=r103+(+r7(28,4)-r6(21,6)*two-r6(21,8)+r5(15,6)+r5(15,8)*two &
&                  +r5(15,12)*three+r5(15,22)*three-r4(10,7)-r4(10,11)*two-r4(10,13) &
&                  -r4(10,31)*six-r4(10,35)*three+r3(6,28)*three+r3(6,32)*six+r3(6,40)*nine &
&                  -r2(5,22)*three-r2(5,30)*six-r2(5,34)*three)*qx+(+r7(21,3)*two+r7(21,4) &
&                  -r6(15,5)*four-r6(15,6)*two-r6(15,7)*two-r6(15,8)+r5(10,5)*two+r5(10,6) &
&                  +r5(10,7)*four+r5(10,8)*two+r5(10,11)*six+r5(10,12)*three+r5(10,21)*two &
&                  +r5(10,22)-r4(6,6)*two-r4(6,7)-r4(6,10)*four-r4(6,11)*two-r4(6,12)*two &
&                  -r4(6,13)-r4(6,30)*four-r4(6,31)*two-r4(6,34)*two-r4(6,35)+r3(3,27)*two &
&                  +r3(3,28)+r3(3,31)*four+r3(3,32)*two+r3(3,39)*six+r3(3,40)*three &
&                  -r2(1,21)*two-r2(1,22)-r2(1,29)*four-r2(1,30)*two-r2(1,33)*two-r2(1,34)) &
&                  *qz+(+r6(21,11)*two+r6(21,12)-r5(15,15)*four-r5(15,16)*two-r5(15,19)*two &
&                  -r5(15,20)+r4(10,16)*two+r4(10,17)+r4(10,20)*four+r4(10,21)*two &
&                  +r4(10,28)*six+r4(10,29)*three+r4(10,40)*two+r4(10,41)-r3(6,13)*two &
&                  -r3(6,14)-r3(6,21)*four-r3(6,22)*two-r3(6,25)*two-r3(6,26)-r3(6,45)*four &
&                  -r3(6,46)*two-r3(6,50)*two-r3(6,51)+r2(5,39)*two+r2(5,40)+r2(5,44)*four &
&                  +r2(5,45)*two+r2(5,54)*six+r2(5,55)*three-r1(1,23)*two-r1(1,24) &
&                  -r1(1,33)*four-r1(1,34)*two-r1(1,38)*two-r1(1,39))*xz+(+r6(15,10) &
&                  +r6(15,11)*two-r5(10,14)*two-r5(10,15)*four-r5(10,18)-r5(10,19)*two &
&                  +r4(6,15)+r4(6,16)*two+r4(6,19)*two+r4(6,20)*four+r4(6,27)*three &
&                  +r4(6,28)*six-r3(3,12)-r3(3,13)*two-r3(3,20)*two-r3(3,21)*four-r3(3,24) &
&                  -r3(3,25)*two)*zz+(+r5(15,23)+r5(15,24)*two-r4(10,32)*two-r4(10,33)*four &
&                  -r4(10,36)-r4(10,37)*two+r3(6,29)+r3(6,30)*two+r3(6,33)*two+r3(6,34)*four &
&                  +r3(6,41)*three+r3(6,42)*six-r2(5,23)-r2(5,24)*two-r2(5,31)*two &
&                  -r2(5,32)*four-r2(5,35)-r2(5,36)*two)*xzz+(+r5(10,23)-r4(6,32)*two &
&                  -r4(6,36)+r3(3,29)+r3(3,33)*two+r3(3,41)*three-r2(1,23)-r2(1,31)*two &
&                  -r2(1,35))*zzz+rxyz(1)*xzzz
      eri(4,5,3,5)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(14,11)+r6(14,12) &
&                  -r5(9,15)*two-r5(9,16)*two-r5(9,19)-r5(9,20)+r4(5,16)+r4(5,17) &
&                  +r4(5,20)*two+r4(5,21)*two+r4(5,28)*three+r4(5,29)*three-r3(2,13)-r3(2,14) &
&                  -r3(2,21)*two-r3(2,22)*two-r3(2,25)-r3(2,26))*xz+rxyz(2)*xxz
      eri(5,5,3,5)=r202+(+r7(21,3)+r7(21,4)-r6(15,5)*two-r6(15,6)*two-r6(15,7)-r6(15,8) &
&                  +r5(10,5)+r5(10,6)+r5(10,7)*two+r5(10,8)*two+r5(10,11)*three &
&                  +r5(10,12)*three+r5(10,21)+r5(10,22)-r4(6,6)-r4(6,7)-r4(6,10)*two &
&                  -r4(6,11)*two-r4(6,12)-r4(6,13)-r4(6,30)*two-r4(6,31)*two-r4(6,34) &
&                  -r4(6,35)+r3(3,27)+r3(3,28)+r3(3,31)*two+r3(3,32)*two+r3(3,39)*three &
&                  +r3(3,40)*three-r2(1,21)-r2(1,22)-r2(1,29)*two-r2(1,30)*two-r2(1,33) &
&                  -r2(1,34))*qx+(+r7(15,3)+r7(15,4)-r6(10,5)*two-r6(10,6)*two-r6(10,7) &
&                  -r6(10,8)+r5(6,5)+r5(6,6)+r5(6,7)*two+r5(6,8)*two+r5(6,11)*three &
&                  +r5(6,12)*three+r5(15,21)+r5(15,22)-r4(3,6)-r4(3,7)-r4(3,10)*two &
&                  -r4(3,11)*two-r4(3,12)-r4(3,13)-r4(10,30)*two-r4(10,31)*two-r4(10,34) &
&                  -r4(10,35)+r3(6,27)+r3(6,28)+r3(6,31)*two+r3(6,32)*two+r3(6,39)*three &
&                  +r3(6,40)*three-r2(5,21)-r2(5,22)-r2(5,29)*two-r2(5,30)*two-r2(5,33) &
&                  -r2(5,34))*qz+(+r6(21,11)-r5(15,15)*two-r5(15,19)+r4(10,16)+r4(10,20)*two &
&                  +r4(10,28)*three+r4(10,40)-r3(6,13)-r3(6,21)*two-r3(6,25)-r3(6,45)*two &
&                  -r3(6,50)+r2(5,39)+r2(5,44)*two+r2(5,54)*three-r1(1,23)-r1(1,33)*two &
&                  -r1(1,38))*xx+(+r6(15,10)+r6(15,11)*two+r6(15,12)-r5(10,14)*two &
&                  -r5(10,15)*four-r5(10,16)*two-r5(10,18)-r5(10,19)*two-r5(10,20)+r4(6,15) &
&                  +r4(6,16)*two+r4(6,17)+r4(6,19)*two+r4(6,20)*four+r4(6,21)*two &
&                  +r4(6,27)*three+r4(6,28)*six+r4(6,29)*three-r3(3,12)-r3(3,13)*two-r3(3,14) &
&                  -r3(3,20)*two-r3(3,21)*four-r3(3,22)*two-r3(3,24)-r3(3,25)*two-r3(3,26)) &
&                  *xz+(+r6(10,11)-r5(6,15)*two-r5(6,19)+r4(3,16)+r4(3,20)*two+r4(3,28)*three &
&                  +r4(10,40)-r3(1,13)-r3(1,21)*two-r3(1,25)-r3(6,45)*two-r3(6,50)+r2(5,39) &
&                  +r2(5,44)*two+r2(5,54)*three-r1(1,23)-r1(1,33)*two-r1(1,38))*zz+( &
&                  +r5(15,23)+r5(15,24)-r4(10,32)*two-r4(10,33)*two-r4(10,36)-r4(10,37) &
&                  +r3(6,29)+r3(6,30)+r3(6,33)*two+r3(6,34)*two+r3(6,41)*three+r3(6,42)*three &
&                  -r2(5,23)-r2(5,24)-r2(5,31)*two-r2(5,32)*two-r2(5,35)-r2(5,36))*xxz+( &
&                  +r5(10,23)+r5(10,24)-r4(6,32)*two-r4(6,33)*two-r4(6,36)-r4(6,37)+r3(3,29) &
&                  +r3(3,30)+r3(3,33)*two+r3(3,34)*two+r3(3,41)*three+r3(3,42)*three-r2(1,23) &
&                  -r2(1,24)-r2(1,31)*two-r2(1,32)*two-r2(1,35)-r2(1,36))*xzz+rxyz(1)*xxzz
      eri(6,5,3,5)=r112+rxyz(19)*qx+(+r7(20,3)+r7(20,4)-r6(14,5)*two-r6(14,6)*two &
&                  -r6(14,7)-r6(14,8)+r5(9,5)+r5(9,6)+r5(9,7)*two+r5(9,8)*two+r5(9,11)*three &
&                  +r5(9,12)*three-r4(5,6)-r4(5,7)-r4(5,10)*two-r4(5,11)*two-r4(5,12) &
&                  -r4(5,13))*qz+(+r6(20,11)+r6(20,12)-r5(14,15)*two-r5(14,16)*two-r5(14,19) &
&                  -r5(14,20)+r4(9,16)+r4(9,17)+r4(9,20)*two+r4(9,21)*two+r4(9,28)*three &
&                  +r4(9,29)*three-r3(5,13)-r3(5,14)-r3(5,21)*two-r3(5,22)*two-r3(5,25) &
&                  -r3(5,26))*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,3,5)=r211+(+r7(20,3)*two-r6(14,5)*four-r6(14,7)*two+r5(9,5)*two &
&                  +r5(9,7)*four+r5(9,11)*six-r4(5,6)*two-r4(5,10)*four-r4(5,12)*two)*qx &
&                  +rxyz(16)*qz+(+r6(20,10)-r5(14,14)*two-r5(14,18)+r4(9,15)+r4(9,19)*two &
&                  +r4(9,27)*three-r3(5,12)-r3(5,20)*two-r3(5,24))*xx+(+r6(14,11)*two &
&                  -r5(9,15)*four-r5(9,19)*two+r4(5,16)*two+r4(5,20)*four+r4(5,28)*six &
&                  -r3(2,13)*two-r3(2,21)*four-r3(2,25)*two)*xz+rxyz(3)*xxz
      eri(2,6,3,5)=r031+rxyz(8)*qz
      eri(3,6,3,5)=r013+(+r7(27,3)*two+r7(27,4)-r6(20,5)*four-r6(20,6)*two-r6(20,7)*two &
&                  -r6(20,8)+r5(14,5)*two+r5(14,6)+r5(14,7)*four+r5(14,8)*two+r5(14,11)*six &
&                  +r5(14,12)*three+r5(14,21)*two+r5(14,22)-r4(9,6)*two-r4(9,7)-r4(9,10)*four &
&                  -r4(9,11)*two-r4(9,12)*two-r4(9,13)-r4(9,30)*four-r4(9,31)*two &
&                  -r4(9,34)*two-r4(9,35)+r3(5,27)*two+r3(5,28)+r3(5,31)*four+r3(5,32)*two &
&                  +r3(5,39)*six+r3(5,40)*three-r2(4,21)*two-r2(4,22)-r2(4,29)*four &
&                  -r2(4,30)*two-r2(4,33)*two-r2(4,34))*qz+(+r6(20,10)+r6(20,11)*two &
&                  -r5(14,14)*two-r5(14,15)*four-r5(14,18)-r5(14,19)*two+r4(9,15) &
&                  +r4(9,16)*two+r4(9,19)*two+r4(9,20)*four+r4(9,27)*three+r4(9,28)*six &
&                  -r3(5,12)-r3(5,13)*two-r3(5,20)*two-r3(5,21)*four-r3(5,24)-r3(5,25)*two) &
&                  *zz+rxyz(3)*zzz
      eri(4,6,3,5)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,3,5)=r112+rxyz(18)*qx+(+r7(20,3)+r7(20,4)-r6(14,5)*two-r6(14,6)*two &
&                  -r6(14,7)-r6(14,8)+r5(9,5)+r5(9,6)+r5(9,7)*two+r5(9,8)*two+r5(9,11)*three &
&                  +r5(9,12)*three-r4(5,6)-r4(5,7)-r4(5,10)*two-r4(5,11)*two-r4(5,12) &
&                  -r4(5,13))*qz+(+r6(20,10)+r6(20,11)-r5(14,14)*two-r5(14,15)*two-r5(14,18) &
&                  -r5(14,19)+r4(9,15)+r4(9,16)+r4(9,19)*two+r4(9,20)*two+r4(9,27)*three &
&                  +r4(9,28)*three-r3(5,12)-r3(5,13)-r3(5,20)*two-r3(5,21)*two-r3(5,24) &
&                  -r3(5,25))*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,3,5)=r022+(+r7(26,3)+r7(26,4)-r6(19,5)*two-r6(19,6)*two-r6(19,7)-r6(19,8) &
&                  +r5(13,5)+r5(13,6)+r5(13,7)*two+r5(13,8)*two+r5(13,11)*three &
&                  +r5(13,12)*three+r5(15,21)+r5(15,22)-r4(8,6)-r4(8,7)-r4(8,10)*two &
&                  -r4(8,11)*two-r4(8,12)-r4(8,13)-r4(10,30)*two-r4(10,31)*two-r4(10,34) &
&                  -r4(10,35)+r3(6,27)+r3(6,28)+r3(6,31)*two+r3(6,32)*two+r3(6,39)*three &
&                  +r3(6,40)*three-r2(5,21)-r2(5,22)-r2(5,29)*two-r2(5,30)*two-r2(5,33) &
&                  -r2(5,34))*qz+rxyz(6)*zz
!
      r400= r8(5)-r7(2,2)+r6(5,4)+r6(5,9)*six-r5(2,4)-r5(2,17)*six+r4(5,26)*six &
&          +r4(5,38)*three-r3(2,23)*six-r3(2,48)*three+r2(6,52)*three-r1(2,36)*three
      r310= r8(8)-r7(4,2)+r6(8,4)+r6(8,9)*three-r5(4,4)-r5(4,17)*three+r4(8,26)*three &
&          -r3(4,23)*three
      r301= r8(9)-r7(5,2)+r6(9,4)+r6(9,9)*three-r5(5,4)-r5(5,17)*three+r4(9,26)*three &
&          -r3(5,23)*three
      r220= r8(12)-r7(7,2)+r6(12,4)+r6(5,9)+r6(12,9)-r5(7,4)-r5(2,17)-r5(7,17)+r4(5,26) &
&          +r4(12,26)+r4(5,38)-r3(2,23)-r3(7,23)-r3(2,48)+r2(6,52)-r1(2,36)
      r211= r8(13)-r7(8,2)+r6(13,4)+r6(13,9)-r5(8,4)-r5(8,17)+r4(13,26)-r3(8,23)
      r202= r8(14)-r7(9,2)+r6(14,4)+r6(5,9)+r6(14,9)-r5(9,4)-r5(2,17)-r5(9,17)+r4(5,26) &
&          +r4(14,26)+r4(5,38)-r3(2,23)-r3(9,23)-r3(2,48)+r2(6,52)-r1(2,36)
      r130= r8(17)-r7(11,2)+r6(17,4)+r6(8,9)*three-r5(11,4)-r5(4,17)*three+r4(8,26)*three &
&          -r3(4,23)*three
      r121= r8(18)-r7(12,2)+r6(18,4)+r6(9,9)-r5(12,4)-r5(5,17)+r4(9,26)-r3(5,23)
      r112= r8(19)-r7(13,2)+r6(19,4)+r6(8,9)-r5(13,4)-r5(4,17)+r4(8,26)-r3(4,23)
      r103= r8(20)-r7(14,2)+r6(20,4)+r6(9,9)*three-r5(14,4)-r5(5,17)*three+r4(9,26)*three &
&          -r3(5,23)*three
      r040= r8(23)-r7(16,2)+r6(23,4)+r6(12,9)*six-r5(16,4)-r5(7,17)*six+r4(12,26)*six &
&          +r4(5,38)*three-r3(7,23)*six-r3(2,48)*three+r2(6,52)*three-r1(2,36)*three
      r031= r8(24)-r7(17,2)+r6(24,4)+r6(13,9)*three-r5(17,4)-r5(8,17)*three+r4(13,26)*three &
&          -r3(8,23)*three
      r022= r8(25)-r7(18,2)+r6(25,4)+r6(12,9)+r6(14,9)-r5(18,4)-r5(7,17)-r5(9,17)+r4(12,26) &
&          +r4(14,26)+r4(5,38)-r3(7,23)-r3(9,23)-r3(2,48)+r2(6,52)-r1(2,36)
      r013= r8(26)-r7(19,2)+r6(26,4)+r6(13,9)*three-r5(19,4)-r5(8,17)*three+r4(13,26)*three &
&          -r3(8,23)*three
      r004= r8(27)-r7(20,2)+r6(27,4)+r6(14,9)*six-r5(20,4)-r5(9,17)*six+r4(14,26)*six &
&          +r4(5,38)*three-r3(9,23)*six-r3(2,48)*three+r2(6,52)*three-r1(2,36)*three
      rxyz(1)=+r4(5,42)-r3(2,52)+r2(6,56)-r1(2,40)
      rxyz(2)=+r5(8,24)-r4(4,37)+r3(8,42)-r2(2,36)
      rxyz(3)=+r5(8,23)-r4(4,36)+r3(8,41)-r2(2,35)
      rxyz(4)=+r6(12,10)-r5(7,18)+r4(12,27)+r4(5,39)-r3(7,24)-r3(2,49)+r2(6,53)-r1(2,37)
      rxyz(5)=+r6(12,12)-r5(7,20)+r4(12,29)+r4(5,41)-r3(7,26)-r3(2,51)+r2(6,55)-r1(2,39)
      rxyz(6)=+r6(12,11)-r5(7,19)+r4(12,28)+r4(5,40)-r3(7,25)-r3(2,50)+r2(6,54)-r1(2,38)
      rxyz(7)=+r7(17,3)-r6(11,7)+r5(17,11)+r5(8,21)*three-r4(11,12)-r4(4,34)*three &
&             +r3(8,39)*three-r2(2,33)*three
      rxyz(8)=+r7(17,4)-r6(11,8)+r5(17,12)+r5(8,22)*three-r4(11,13)-r4(4,35)*three &
&             +r3(8,40)*three-r2(2,34)*three
      rxyz(9)=+r7(13,3)+r7(13,4)-r6(8,7)-r6(8,8)+r5(13,11)+r5(13,12)-r4(8,12)-r4(8,13)
      rxyz(10)=+r6(13,11)-r5(8,19)+r4(13,28)-r3(8,25)
      rxyz(11)=+r6(8,11)-r5(4,19)+r4(8,28)-r3(4,25)
      rxyz(12)=+r7(18,3)-r6(12,7)+r5(18,11)+r5(9,21)-r4(12,12)-r4(5,34)+r3(9,39)-r2(6,33)
      rxyz(13)=+r7(18,4)-r6(12,8)+r5(18,12)+r5(9,22)-r4(12,13)-r4(5,35)+r3(9,40)-r2(6,34)
      rxyz(14)=+r7(12,4)-r6(7,8)+r5(12,12)+r5(5,22)-r4(7,13)-r4(2,35)+r3(5,40)-r2(4,34)
      rxyz(15)=+r7(12,3)-r6(7,7)+r5(12,11)+r5(5,21)-r4(7,12)-r4(2,34)+r3(5,39)-r2(4,33)
      rxyz(16)=+r7(8,4)-r6(4,8)+r5(8,12)+r5(8,22)-r4(4,13)-r4(4,35)+r3(8,40)-r2(2,34)
      rxyz(17)=+r7(8,3)-r6(4,7)+r5(8,11)+r5(8,21)-r4(4,12)-r4(4,34)+r3(8,39)-r2(2,33)
      rxyz(18)=+r7(19,3)-r6(13,7)+r5(19,11)+r5(8,21)-r4(13,12)-r4(4,34)+r3(8,39)-r2(2,33)
      rxyz(19)=+r7(19,4)-r6(13,8)+r5(19,12)+r5(8,22)-r4(13,13)-r4(4,35)+r3(8,40)-r2(2,34)
      rxyz(20)=+r6(9,11)*four-r5(5,19)*four+r4(9,28)*four-r3(5,25)*four
      eri(1,1,4,5)=r400+(+r7(5,3)*two+r7(5,4)*two-r6(2,7)*two-r6(2,8)*two+r5(5,11)*two &
&                  +r5(5,12)*two+r5(5,21)*six+r5(5,22)*six-r4(2,12)*two-r4(2,13)*two &
&                  -r4(2,34)*six-r4(2,35)*six+r3(5,39)*six+r3(5,40)*six-r2(4,33)*six &
&                  -r2(4,34)*six)*qx+(+r6(5,10)+r6(5,11)*four+r6(5,12)-r5(2,18)-r5(2,19)*four &
&                  -r5(2,20)+r4(5,27)+r4(5,28)*four+r4(5,29)+r4(5,39)+r4(5,40)*four+r4(5,41) &
&                  -r3(2,24)-r3(2,25)*four-r3(2,26)-r3(2,49)-r3(2,50)*four-r3(2,51)+r2(6,53) &
&                  +r2(6,54)*four+r2(6,55)-r1(2,37)-r1(2,38)*four-r1(2,39))*xx+(+r5(5,23)*two &
&                  +r5(5,24)*two-r4(2,36)*two-r4(2,37)*two+r3(5,41)*two+r3(5,42)*two &
&                  -r2(4,35)*two-r2(4,36)*two)*xxx+rxyz(1)*xxxx
      eri(2,1,4,5)=r220+(+r7(12,4)*two-r6(7,8)*two+r5(12,12)*two+r5(5,22)*two &
&                  -r4(7,13)*two-r4(2,35)*two+r3(5,40)*two-r2(4,34)*two)*qx+rxyz(5)*xx
      eri(3,1,4,5)=r202+(+r7(14,4)*two-r6(9,8)*two+r5(14,12)*two+r5(5,22)*two &
&                  -r4(9,13)*two-r4(2,35)*two+r3(5,40)*two-r2(4,34)*two)*qx+(+r7(9,3)*two &
&                  -r6(5,7)*two+r5(9,11)*two+r5(9,21)*two-r4(5,12)*two-r4(5,34)*two &
&                  +r3(9,39)*two-r2(6,33)*two)*qz+(+r6(14,12)-r5(9,20)+r4(14,29)+r4(5,41) &
&                  -r3(9,26)-r3(2,51)+r2(6,55)-r1(2,39))*xx+rxyz(20)*xz+(+r6(5,10)-r5(2,18) &
&                  +r4(5,27)+r4(5,39)-r3(2,24)-r3(2,49)+r2(6,53)-r1(2,37))*zz+(+r5(9,24)*two &
&                  -r4(5,37)*two+r3(9,42)*two-r2(6,36)*two)*xxz+(+r5(5,23)*two-r4(2,36)*two &
&                  +r3(5,41)*two-r2(4,35)*two)*xzz+rxyz(1)*xxzz
      eri(4,1,4,5)=r310+(+r7(8,3)+r7(8,4)*two-r6(4,7)-r6(4,8)*two+r5(8,11)+r5(8,12)*two &
&                  +r5(8,21)+r5(8,22)*two-r4(4,12)-r4(4,13)*two-r4(4,34)-r4(4,35)*two &
&                  +r3(8,39)+r3(8,40)*two-r2(2,33)-r2(2,34)*two)*qx+(+r6(8,11)*two+r6(8,12) &
&                  -r5(4,19)*two-r5(4,20)+r4(8,28)*two+r4(8,29)-r3(4,25)*two-r3(4,26))*xx &
&                  +rxyz(2)*xxx
      eri(5,1,4,5)=r301+(+r7(9,3)+r7(9,4)*two-r6(5,7)-r6(5,8)*two+r5(9,11)+r5(9,12)*two &
&                  +r5(9,21)+r5(9,22)*two-r4(5,12)-r4(5,13)*two-r4(5,34)-r4(5,35)*two &
&                  +r3(9,39)+r3(9,40)*two-r2(6,33)-r2(6,34)*two)*qx+(+r7(5,3)-r6(2,7) &
&                  +r5(5,11)+r5(5,21)*three-r4(2,12)-r4(2,34)*three+r3(5,39)*three &
&                  -r2(4,33)*three)*qz+(+r6(9,11)*two+r6(9,12)-r5(5,19)*two-r5(5,20) &
&                  +r4(9,28)*two+r4(9,29)-r3(5,25)*two-r3(5,26))*xx+(+r6(5,10)+r6(5,11)*two &
&                  -r5(2,18)-r5(2,19)*two+r4(5,27)+r4(5,28)*two+r4(5,39)+r4(5,40)*two &
&                  -r3(2,24)-r3(2,25)*two-r3(2,49)-r3(2,50)*two+r2(6,53)+r2(6,54)*two &
&                  -r1(2,37)-r1(2,38)*two)*xz+(+r5(9,24)-r4(5,37)+r3(9,42)-r2(6,36))*xxx+( &
&                  +r5(5,23)*two+r5(5,24)-r4(2,36)*two-r4(2,37)+r3(5,41)*two+r3(5,42) &
&                  -r2(4,35)*two-r2(4,36))*xxz+rxyz(1)*xxxz
      eri(6,1,4,5)=r211+(+r7(13,4)*two-r6(8,8)*two+r5(13,12)*two-r4(8,13)*two)*qx &
&                  +rxyz(17)*qz+(+r6(13,12)-r5(8,20)+r4(13,29)-r3(8,26))*xx+(+r6(8,11)*two &
&                  -r5(4,19)*two+r4(8,28)*two-r3(4,25)*two)*xz+rxyz(2)*xxz
      eri(1,2,4,5)=r220+(+r7(12,3)*two-r6(7,7)*two+r5(12,11)*two+r5(5,21)*two &
&                  -r4(7,12)*two-r4(2,34)*two+r3(5,39)*two-r2(4,33)*two)*qx+rxyz(4)*xx
      eri(2,2,4,5)=r040
      eri(3,2,4,5)=r022+(+r7(18,3)*two-r6(12,7)*two+r5(18,11)*two+r5(9,21)*two &
&                  -r4(12,12)*two-r4(5,34)*two+r3(9,39)*two-r2(6,33)*two)*qz+rxyz(4)*zz
      eri(4,2,4,5)=r130+rxyz(7)*qx
      eri(5,2,4,5)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,4,5)=r031+rxyz(7)*qz
      eri(1,3,4,5)=r202+(+r7(14,3)*two-r6(9,7)*two+r5(14,11)*two+r5(5,21)*two &
&                  -r4(9,12)*two-r4(2,34)*two+r3(5,39)*two-r2(4,33)*two)*qx+(+r7(9,4)*two &
&                  -r6(5,8)*two+r5(9,12)*two+r5(9,22)*two-r4(5,13)*two-r4(5,35)*two &
&                  +r3(9,40)*two-r2(6,34)*two)*qz+(+r6(14,10)-r5(9,18)+r4(14,27)+r4(5,39) &
&                  -r3(9,24)-r3(2,49)+r2(6,53)-r1(2,37))*xx+rxyz(20)*xz+(+r6(5,12)-r5(2,20) &
&                  +r4(5,29)+r4(5,41)-r3(2,26)-r3(2,51)+r2(6,55)-r1(2,39))*zz+(+r5(9,23)*two &
&                  -r4(5,36)*two+r3(9,41)*two-r2(6,35)*two)*xxz+(+r5(5,24)*two-r4(2,37)*two &
&                  +r3(5,42)*two-r2(4,36)*two)*xzz+rxyz(1)*xxzz
      eri(2,3,4,5)=r022+(+r7(18,4)*two-r6(12,8)*two+r5(18,12)*two+r5(9,22)*two &
&                  -r4(12,13)*two-r4(5,35)*two+r3(9,40)*two-r2(6,34)*two)*qz+rxyz(5)*zz
      eri(3,3,4,5)=r004+(+r7(20,3)*two+r7(20,4)*two-r6(14,7)*two-r6(14,8)*two &
&                  +r5(20,11)*two+r5(20,12)*two+r5(9,21)*six+r5(9,22)*six-r4(14,12)*two &
&                  -r4(14,13)*two-r4(5,34)*six-r4(5,35)*six+r3(9,39)*six+r3(9,40)*six &
&                  -r2(6,33)*six-r2(6,34)*six)*qz+(+r6(14,10)+r6(14,11)*four+r6(14,12) &
&                  -r5(9,18)-r5(9,19)*four-r5(9,20)+r4(14,27)+r4(14,28)*four+r4(14,29) &
&                  +r4(5,39)+r4(5,40)*four+r4(5,41)-r3(9,24)-r3(9,25)*four-r3(9,26)-r3(2,49) &
&                  -r3(2,50)*four-r3(2,51)+r2(6,53)+r2(6,54)*four+r2(6,55)-r1(2,37) &
&                  -r1(2,38)*four-r1(2,39))*zz+(+r5(9,23)*two+r5(9,24)*two-r4(5,36)*two &
&                  -r4(5,37)*two+r3(9,41)*two+r3(9,42)*two-r2(6,35)*two-r2(6,36)*two)*zzz &
&                  +rxyz(1)*zzzz
      eri(4,3,4,5)=r112+rxyz(18)*qx+(+r7(13,4)*two-r6(8,8)*two+r5(13,12)*two-r4(8,13)*two &
&                  )*qz+(+r6(13,11)*two-r5(8,19)*two+r4(13,28)*two-r3(8,25)*two)*xz+( &
&                  +r6(8,12)-r5(4,20)+r4(8,29)-r3(4,26))*zz+rxyz(2)*xzz
      eri(5,3,4,5)=r103+(+r7(20,3)-r6(14,7)+r5(20,11)+r5(9,21)*three-r4(14,12) &
&                  -r4(5,34)*three+r3(9,39)*three-r2(6,33)*three)*qx+(+r7(14,3)+r7(14,4)*two &
&                  -r6(9,7)-r6(9,8)*two+r5(14,11)+r5(14,12)*two+r5(5,21)+r5(5,22)*two &
&                  -r4(9,12)-r4(9,13)*two-r4(2,34)-r4(2,35)*two+r3(5,39)+r3(5,40)*two &
&                  -r2(4,33)-r2(4,34)*two)*qz+(+r6(14,10)+r6(14,11)*two-r5(9,18)-r5(9,19)*two &
&                  +r4(14,27)+r4(14,28)*two+r4(5,39)+r4(5,40)*two-r3(9,24)-r3(9,25)*two &
&                  -r3(2,49)-r3(2,50)*two+r2(6,53)+r2(6,54)*two-r1(2,37)-r1(2,38)*two)*xz+( &
&                  +r6(9,11)*two+r6(9,12)-r5(5,19)*two-r5(5,20)+r4(9,28)*two+r4(9,29) &
&                  -r3(5,25)*two-r3(5,26))*zz+(+r5(9,23)*two+r5(9,24)-r4(5,36)*two-r4(5,37) &
&                  +r3(9,41)*two+r3(9,42)-r2(6,35)*two-r2(6,36))*xzz+(+r5(5,24)-r4(2,37) &
&                  +r3(5,42)-r2(4,36))*zzz+rxyz(1)*xzzz
      eri(6,3,4,5)=r013+(+r7(19,3)+r7(19,4)*two-r6(13,7)-r6(13,8)*two+r5(19,11) &
&                  +r5(19,12)*two+r5(8,21)+r5(8,22)*two-r4(13,12)-r4(13,13)*two-r4(4,34) &
&                  -r4(4,35)*two+r3(8,39)+r3(8,40)*two-r2(2,33)-r2(2,34)*two)*qz+( &
&                  +r6(13,11)*two+r6(13,12)-r5(8,19)*two-r5(8,20)+r4(13,28)*two+r4(13,29) &
&                  -r3(8,25)*two-r3(8,26))*zz+rxyz(2)*zzz
      eri(1,4,4,5)=r310+(+r7(8,3)*two+r7(8,4)-r6(4,7)*two-r6(4,8)+r5(8,11)*two+r5(8,12) &
&                  +r5(8,21)*two+r5(8,22)-r4(4,12)*two-r4(4,13)-r4(4,34)*two-r4(4,35) &
&                  +r3(8,39)*two+r3(8,40)-r2(2,33)*two-r2(2,34))*qx+(+r6(8,10)+r6(8,11)*two &
&                  -r5(4,18)-r5(4,19)*two+r4(8,27)+r4(8,28)*two-r3(4,24)-r3(4,25)*two)*xx &
&                  +rxyz(3)*xxx
      eri(2,4,4,5)=r130+rxyz(8)*qx
      eri(3,4,4,5)=r112+rxyz(19)*qx+(+r7(13,3)*two-r6(8,7)*two+r5(13,11)*two-r4(8,12)*two &
&                  )*qz+(+r6(13,11)*two-r5(8,19)*two+r4(13,28)*two-r3(8,25)*two)*xz+( &
&                  +r6(8,10)-r5(4,18)+r4(8,27)-r3(4,24))*zz+rxyz(3)*xzz
      eri(4,4,4,5)=r220+(+r7(12,3)+r7(12,4)-r6(7,7)-r6(7,8)+r5(12,11)+r5(12,12)+r5(5,21) &
&                  +r5(5,22)-r4(7,12)-r4(7,13)-r4(2,34)-r4(2,35)+r3(5,39)+r3(5,40)-r2(4,33) &
&                  -r2(4,34))*qx+rxyz(6)*xx
      eri(5,4,4,5)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(8,10)+r6(8,11)-r5(4,18) &
&                  -r5(4,19)+r4(8,27)+r4(8,28)-r3(4,24)-r3(4,25))*xz+rxyz(3)*xxz
      eri(6,4,4,5)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,4,5)=r301+(+r7(9,3)*two+r7(9,4)-r6(5,7)*two-r6(5,8)+r5(9,11)*two+r5(9,12) &
&                  +r5(9,21)*two+r5(9,22)-r4(5,12)*two-r4(5,13)-r4(5,34)*two-r4(5,35) &
&                  +r3(9,39)*two+r3(9,40)-r2(6,33)*two-r2(6,34))*qx+(+r7(5,4)-r6(2,8) &
&                  +r5(5,12)+r5(5,22)*three-r4(2,13)-r4(2,35)*three+r3(5,40)*three &
&                  -r2(4,34)*three)*qz+(+r6(9,10)+r6(9,11)*two-r5(5,18)-r5(5,19)*two+r4(9,27) &
&                  +r4(9,28)*two-r3(5,24)-r3(5,25)*two)*xx+(+r6(5,11)*two+r6(5,12) &
&                  -r5(2,19)*two-r5(2,20)+r4(5,28)*two+r4(5,29)+r4(5,40)*two+r4(5,41) &
&                  -r3(2,25)*two-r3(2,26)-r3(2,50)*two-r3(2,51)+r2(6,54)*two+r2(6,55) &
&                  -r1(2,38)*two-r1(2,39))*xz+(+r5(9,23)-r4(5,36)+r3(9,41)-r2(6,35))*xxx+( &
&                  +r5(5,23)+r5(5,24)*two-r4(2,36)-r4(2,37)*two+r3(5,41)+r3(5,42)*two &
&                  -r2(4,35)-r2(4,36)*two)*xxz+rxyz(1)*xxxz
      eri(2,5,4,5)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,4,5)=r103+(+r7(20,4)-r6(14,8)+r5(20,12)+r5(9,22)*three-r4(14,13) &
&                  -r4(5,35)*three+r3(9,40)*three-r2(6,34)*three)*qx+(+r7(14,3)*two+r7(14,4) &
&                  -r6(9,7)*two-r6(9,8)+r5(14,11)*two+r5(14,12)+r5(5,21)*two+r5(5,22) &
&                  -r4(9,12)*two-r4(9,13)-r4(2,34)*two-r4(2,35)+r3(5,39)*two+r3(5,40) &
&                  -r2(4,33)*two-r2(4,34))*qz+(+r6(14,11)*two+r6(14,12)-r5(9,19)*two-r5(9,20) &
&                  +r4(14,28)*two+r4(14,29)+r4(5,40)*two+r4(5,41)-r3(9,25)*two-r3(9,26) &
&                  -r3(2,50)*two-r3(2,51)+r2(6,54)*two+r2(6,55)-r1(2,38)*two-r1(2,39))*xz+( &
&                  +r6(9,10)+r6(9,11)*two-r5(5,18)-r5(5,19)*two+r4(9,27)+r4(9,28)*two &
&                  -r3(5,24)-r3(5,25)*two)*zz+(+r5(9,23)+r5(9,24)*two-r4(5,36)-r4(5,37)*two &
&                  +r3(9,41)+r3(9,42)*two-r2(6,35)-r2(6,36)*two)*xzz+(+r5(5,23)-r4(2,36) &
&                  +r3(5,41)-r2(4,35))*zzz+rxyz(1)*xzzz
      eri(4,5,4,5)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(8,11)+r6(8,12)-r5(4,19) &
&                  -r5(4,20)+r4(8,28)+r4(8,29)-r3(4,25)-r3(4,26))*xz+rxyz(2)*xxz
      eri(5,5,4,5)=r202+(+r7(14,3)+r7(14,4)-r6(9,7)-r6(9,8)+r5(14,11)+r5(14,12)+r5(5,21) &
&                  +r5(5,22)-r4(9,12)-r4(9,13)-r4(2,34)-r4(2,35)+r3(5,39)+r3(5,40)-r2(4,33) &
&                  -r2(4,34))*qx+(+r7(9,3)+r7(9,4)-r6(5,7)-r6(5,8)+r5(9,11)+r5(9,12)+r5(9,21) &
&                  +r5(9,22)-r4(5,12)-r4(5,13)-r4(5,34)-r4(5,35)+r3(9,39)+r3(9,40)-r2(6,33) &
&                  -r2(6,34))*qz+(+r6(14,11)-r5(9,19)+r4(14,28)+r4(5,40)-r3(9,25)-r3(2,50) &
&                  +r2(6,54)-r1(2,38))*xx+(+r6(9,10)+r6(9,11)*two+r6(9,12)-r5(5,18) &
&                  -r5(5,19)*two-r5(5,20)+r4(9,27)+r4(9,28)*two+r4(9,29)-r3(5,24) &
&                  -r3(5,25)*two-r3(5,26))*xz+(+r6(5,11)-r5(2,19)+r4(5,28)+r4(5,40)-r3(2,25) &
&                  -r3(2,50)+r2(6,54)-r1(2,38))*zz+(+r5(9,23)+r5(9,24)-r4(5,36)-r4(5,37) &
&                  +r3(9,41)+r3(9,42)-r2(6,35)-r2(6,36))*xxz+(+r5(5,23)+r5(5,24)-r4(2,36) &
&                  -r4(2,37)+r3(5,41)+r3(5,42)-r2(4,35)-r2(4,36))*xzz+rxyz(1)*xxzz
      eri(6,5,4,5)=r112+rxyz(19)*qx+(+r7(13,3)+r7(13,4)-r6(8,7)-r6(8,8)+r5(13,11) &
&                  +r5(13,12)-r4(8,12)-r4(8,13))*qz+(+r6(13,11)+r6(13,12)-r5(8,19)-r5(8,20) &
&                  +r4(13,28)+r4(13,29)-r3(8,25)-r3(8,26))*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,4,5)=r211+(+r7(13,3)*two-r6(8,7)*two+r5(13,11)*two-r4(8,12)*two)*qx &
&                  +rxyz(16)*qz+(+r6(13,10)-r5(8,18)+r4(13,27)-r3(8,24))*xx+(+r6(8,11)*two &
&                  -r5(4,19)*two+r4(8,28)*two-r3(4,25)*two)*xz+rxyz(3)*xxz
      eri(2,6,4,5)=r031+rxyz(8)*qz
      eri(3,6,4,5)=r013+(+r7(19,3)*two+r7(19,4)-r6(13,7)*two-r6(13,8)+r5(19,11)*two &
&                  +r5(19,12)+r5(8,21)*two+r5(8,22)-r4(13,12)*two-r4(13,13)-r4(4,34)*two &
&                  -r4(4,35)+r3(8,39)*two+r3(8,40)-r2(2,33)*two-r2(2,34))*qz+(+r6(13,10) &
&                  +r6(13,11)*two-r5(8,18)-r5(8,19)*two+r4(13,27)+r4(13,28)*two-r3(8,24) &
&                  -r3(8,25)*two)*zz+rxyz(3)*zzz
      eri(4,6,4,5)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,4,5)=r112+rxyz(18)*qx+(+r7(13,3)+r7(13,4)-r6(8,7)-r6(8,8)+r5(13,11) &
&                  +r5(13,12)-r4(8,12)-r4(8,13))*qz+(+r6(13,10)+r6(13,11)-r5(8,18)-r5(8,19) &
&                  +r4(13,27)+r4(13,28)-r3(8,24)-r3(8,25))*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,4,5)=r022+(+r7(18,3)+r7(18,4)-r6(12,7)-r6(12,8)+r5(18,11)+r5(18,12) &
&                  +r5(9,21)+r5(9,22)-r4(12,12)-r4(12,13)-r4(5,34)-r4(5,35)+r3(9,39)+r3(9,40) &
&                  -r2(6,33)-r2(6,34))*qz+rxyz(6)*zz
!
      r400= r8(6)-r7(3,1)-r7(3,2)+r6(1,2)+r6(1,4)+r6(6,4)+r6(6,9)*six-r5(3,3)-r5(3,4) &
&          -r5(3,13)*six-r5(3,17)*six+r4(1,3)+r4(1,5)+r4(1,18)*six+r4(1,26)*six &
&          +r4(6,26)*six+r4(6,38)*three-r3(3,19)*six-r3(3,23)*six-r3(3,43)*three &
&          -r3(3,48)*three+r2(1,9)*six+r2(1,17)*six+r2(1,42)*three+r2(1,52)*three &
&          +r2(3,52)*three-r1(3,31)*three-r1(3,36)*three+r0(11)*three+r0(21)*three
      r310= r8(9)-r7(5,1)-r7(5,2)+r6(2,2)+r6(2,4)+r6(9,4)+r6(9,9)*three-r5(5,3)-r5(5,4) &
&          -r5(5,13)*three-r5(5,17)*three+r4(2,3)+r4(2,5)+r4(2,18)*three+r4(2,26)*three &
&          +r4(9,26)*three-r3(5,19)*three-r3(5,23)*three+r2(4,9)*three+r2(4,17)*three
      r301= r8(10)-r7(6,1)-r7(6,2)+r6(3,2)+r6(3,4)+r6(10,4)+r6(10,9)*three-r5(6,3)-r5(6,4) &
&          -r5(6,13)*three-r5(6,17)*three+r4(3,3)+r4(3,5)+r4(3,18)*three+r4(3,26)*three &
&          +r4(10,26)*three-r3(6,19)*three-r3(6,23)*three+r2(5,9)*three+r2(5,17)*three
      r220= r8(13)-r7(8,1)-r7(8,2)+r6(4,2)+r6(4,4)+r6(13,4)+r6(6,9)+r6(13,9)-r5(8,3) &
&          -r5(8,4)-r5(3,13)-r5(8,13)-r5(3,17)-r5(8,17)+r4(4,3)+r4(4,5)+r4(1,18)+r4(4,18) &
&          +r4(1,26)+r4(4,26)+r4(6,26)+r4(13,26)+r4(6,38)-r3(3,19)-r3(8,19)-r3(3,23) &
&          -r3(8,23)-r3(3,43)-r3(3,48)+r2(1,9)+r2(2,9)+r2(1,17)+r2(2,17)+r2(1,42)+r2(1,52) &
&          +r2(3,52)-r1(3,31)-r1(3,36)+r0(11)+r0(21)
      r211= r8(14)-r7(9,1)-r7(9,2)+r6(5,2)+r6(5,4)+r6(14,4)+r6(14,9)-r5(9,3)-r5(9,4) &
&          -r5(9,13)-r5(9,17)+r4(5,3)+r4(5,5)+r4(5,18)+r4(5,26)+r4(14,26)-r3(9,19) &
&          -r3(9,23)+r2(6,9)+r2(6,17)
      r202= r8(15)-r7(10,1)-r7(10,2)+r6(6,2)+r6(6,4)+r6(15,4)+r6(6,9)+r6(15,9)-r5(10,3) &
&          -r5(10,4)-r5(3,13)-r5(10,13)-r5(3,17)-r5(10,17)+r4(6,3)+r4(6,5)+r4(1,18) &
&          +r4(6,18)+r4(1,26)+r4(6,26)*two+r4(15,26)+r4(6,38)-r3(3,19)-r3(10,19)-r3(3,23) &
&          -r3(10,23)-r3(3,43)-r3(3,48)+r2(1,9)+r2(3,9)+r2(1,17)+r2(3,17)+r2(1,42) &
&          +r2(1,52)+r2(3,52)-r1(3,31)-r1(3,36)+r0(11)+r0(21)
      r130= r8(18)-r7(12,1)-r7(12,2)+r6(7,2)+r6(7,4)+r6(18,4)+r6(9,9)*three-r5(12,3) &
&          -r5(12,4)-r5(5,13)*three-r5(5,17)*three+r4(7,3)+r4(7,5)+r4(2,18)*three &
&          +r4(2,26)*three+r4(9,26)*three-r3(5,19)*three-r3(5,23)*three+r2(4,9)*three &
&          +r2(4,17)*three
      r121= r8(19)-r7(13,1)-r7(13,2)+r6(8,2)+r6(8,4)+r6(19,4)+r6(10,9)-r5(13,3)-r5(13,4) &
&          -r5(6,13)-r5(6,17)+r4(8,3)+r4(8,5)+r4(3,18)+r4(3,26)+r4(10,26)-r3(6,19) &
&          -r3(6,23)+r2(5,9)+r2(5,17)
      r112= r8(20)-r7(14,1)-r7(14,2)+r6(9,2)+r6(9,4)+r6(20,4)+r6(9,9)-r5(14,3)-r5(14,4) &
&          -r5(5,13)-r5(5,17)+r4(9,3)+r4(9,5)+r4(2,18)+r4(2,26)+r4(9,26)-r3(5,19)-r3(5,23) &
&          +r2(4,9)+r2(4,17)
      r103= r8(21)-r7(15,1)-r7(15,2)+r6(10,2)+r6(10,4)+r6(21,4)+r6(10,9)*three-r5(15,3) &
&          -r5(15,4)-r5(6,13)*three-r5(6,17)*three+r4(10,3)+r4(10,5)+r4(3,18)*three &
&          +r4(3,26)*three+r4(10,26)*three-r3(6,19)*three-r3(6,23)*three+r2(5,9)*three &
&          +r2(5,17)*three
      r040= r8(24)-r7(17,1)-r7(17,2)+r6(11,2)+r6(11,4)+r6(24,4)+r6(13,9)*six-r5(17,3) &
&          -r5(17,4)-r5(8,13)*six-r5(8,17)*six+r4(11,3)+r4(11,5)+r4(4,18)*six+r4(4,26)*six &
&          +r4(13,26)*six+r4(6,38)*three-r3(8,19)*six-r3(8,23)*six-r3(3,43)*three &
&          -r3(3,48)*three+r2(2,9)*six+r2(2,17)*six+r2(1,42)*three+r2(1,52)*three &
&          +r2(3,52)*three-r1(3,31)*three-r1(3,36)*three+r0(11)*three+r0(21)*three
      r031= r8(25)-r7(18,1)-r7(18,2)+r6(12,2)+r6(12,4)+r6(25,4)+r6(14,9)*three-r5(18,3) &
&          -r5(18,4)-r5(9,13)*three-r5(9,17)*three+r4(12,3)+r4(12,5)+r4(5,18)*three &
&          +r4(5,26)*three+r4(14,26)*three-r3(9,19)*three-r3(9,23)*three+r2(6,9)*three &
&          +r2(6,17)*three
      r022= r8(26)-r7(19,1)-r7(19,2)+r6(13,2)+r6(13,4)+r6(26,4)+r6(13,9)+r6(15,9)-r5(19,3) &
&          -r5(19,4)-r5(8,13)-r5(10,13)-r5(8,17)-r5(10,17)+r4(13,3)+r4(13,5)+r4(4,18) &
&          +r4(6,18)+r4(4,26)+r4(6,26)+r4(13,26)+r4(15,26)+r4(6,38)-r3(8,19)-r3(10,19) &
&          -r3(8,23)-r3(10,23)-r3(3,43)-r3(3,48)+r2(2,9)+r2(3,9)+r2(2,17)+r2(3,17) &
&          +r2(1,42)+r2(1,52)+r2(3,52)-r1(3,31)-r1(3,36)+r0(11)+r0(21)
      r013= r8(27)-r7(20,1)-r7(20,2)+r6(14,2)+r6(14,4)+r6(27,4)+r6(14,9)*three-r5(20,3) &
&          -r5(20,4)-r5(9,13)*three-r5(9,17)*three+r4(14,3)+r4(14,5)+r4(5,18)*three &
&          +r4(5,26)*three+r4(14,26)*three-r3(9,19)*three-r3(9,23)*three+r2(6,9)*three &
&          +r2(6,17)*three
      r004= r8(28)-r7(21,1)-r7(21,2)+r6(15,2)+r6(15,4)+r6(28,4)+r6(15,9)*six-r5(21,3) &
&          -r5(21,4)-r5(10,13)*six-r5(10,17)*six+r4(15,3)+r4(15,5)+r4(6,18)*six &
&          +r4(6,26)*six+r4(15,26)*six+r4(6,38)*three-r3(10,19)*six-r3(10,23)*six &
&          -r3(3,43)*three-r3(3,48)*three+r2(3,9)*six+r2(3,17)*six+r2(1,42)*three &
&          +r2(1,52)*three+r2(3,52)*three-r1(3,31)*three-r1(3,36)*three+r0(11)*three &
&          +r0(21)*three
      rxyz(1)=+r4(6,42)-r3(3,47)-r3(3,52)+r2(1,46)+r2(1,56)+r2(3,56)-r1(3,35)-r1(3,40) &
&             +r0(15)+r0(25)
      rxyz(2)=+r5(9,24)-r4(5,33)-r4(5,37)+r3(2,34)+r3(2,42)+r3(9,42)-r2(6,32)-r2(6,36) &
&             +r1(2,12)+r1(2,20)
      rxyz(3)=+r5(9,23)-r4(5,32)-r4(5,36)+r3(2,33)+r3(2,41)+r3(9,41)-r2(6,31)-r2(6,35) &
&             +r1(2,11)+r1(2,19)
      rxyz(4)=+r6(13,10)-r5(8,14)-r5(8,18)+r4(4,19)+r4(4,27)+r4(13,27)+r4(6,39)-r3(8,20) &
&             -r3(8,24)-r3(3,44)-r3(3,49)+r2(2,10)+r2(2,18)+r2(1,43)+r2(1,53)+r2(3,53) &
&             -r1(3,32)-r1(3,37)+r0(12)+r0(22)
      rxyz(5)=+r6(13,12)-r5(8,16)-r5(8,20)+r4(4,21)+r4(4,29)+r4(13,29)+r4(6,41)-r3(8,22) &
&             -r3(8,26)-r3(3,46)-r3(3,51)+r2(2,12)+r2(2,20)+r2(1,45)+r2(1,55)+r2(3,55) &
&             -r1(3,34)-r1(3,39)+r0(14)+r0(24)
      rxyz(6)=+r6(13,11)-r5(8,15)-r5(8,19)+r4(4,20)+r4(4,28)+r4(13,28)+r4(6,40)-r3(8,21) &
&             -r3(8,25)-r3(3,45)-r3(3,50)+r2(2,11)+r2(2,19)+r2(1,44)+r2(1,54)+r2(3,54) &
&             -r1(3,33)-r1(3,38)+r0(13)+r0(23)
      rxyz(7)=+r7(18,3)-r6(12,5)-r6(12,7)+r5(7,7)+r5(7,11)+r5(18,11)+r5(9,21)*three &
&             -r4(12,10)-r4(12,12)-r4(5,30)*three-r4(5,34)*three+r3(7,5)+r3(7,9) &
&             +r3(2,31)*three+r3(2,39)*three+r3(9,39)*three-r2(6,29)*three-r2(6,33)*three &
&             +r1(2,9)*three+r1(2,17)*three
      rxyz(8)=+r7(18,4)-r6(12,6)-r6(12,8)+r5(7,8)+r5(7,12)+r5(18,12)+r5(9,22)*three &
&             -r4(12,11)-r4(12,13)-r4(5,31)*three-r4(5,35)*three+r3(7,6)+r3(7,10) &
&             +r3(2,32)*three+r3(2,40)*three+r3(9,40)*three-r2(6,30)*three-r2(6,34)*three &
&             +r1(2,10)*three+r1(2,18)*three
      rxyz(9)=+r7(14,3)+r7(14,4)-r6(9,5)-r6(9,6)-r6(9,7)-r6(9,8)+r5(5,7)+r5(5,8)+r5(5,11) &
&             +r5(14,11)+r5(5,12)+r5(14,12)-r4(9,10)-r4(9,11)-r4(9,12)-r4(9,13)+r3(5,5) &
&             +r3(5,6)+r3(5,9)+r3(5,10)
      rxyz(10)=+r6(14,11)-r5(9,15)-r5(9,19)+r4(5,20)+r4(5,28)+r4(14,28)-r3(9,21)-r3(9,25) &
&             +r2(6,11)+r2(6,19)
      rxyz(11)=+r6(9,11)-r5(5,15)-r5(5,19)+r4(2,20)+r4(2,28)+r4(9,28)-r3(5,21)-r3(5,25) &
&             +r2(4,11)+r2(4,19)
      rxyz(12)=+r7(19,3)-r6(13,5)-r6(13,7)+r5(8,7)+r5(8,11)+r5(19,11)+r5(10,21)-r4(13,10) &
&             -r4(13,12)-r4(6,30)-r4(6,34)+r3(8,5)+r3(8,9)+r3(3,31)+r3(3,39)+r3(10,39) &
&             -r2(3,29)-r2(3,33)+r1(3,9)+r1(3,17)
      rxyz(13)=+r7(19,4)-r6(13,6)-r6(13,8)+r5(8,8)+r5(8,12)+r5(19,12)+r5(10,22)-r4(13,11) &
&             -r4(13,13)-r4(6,31)-r4(6,35)+r3(8,6)+r3(8,10)+r3(3,32)+r3(3,40)+r3(10,40) &
&             -r2(3,30)-r2(3,34)+r1(3,10)+r1(3,18)
      rxyz(14)=+r7(13,4)-r6(8,6)-r6(8,8)+r5(4,8)+r5(4,12)+r5(13,12)+r5(6,22)-r4(8,11) &
&             -r4(8,13)-r4(3,31)-r4(3,35)+r3(4,6)+r3(4,10)+r3(1,32)+r3(1,40)+r3(6,40) &
&             -r2(5,30)-r2(5,34)+r1(1,10)+r1(1,18)
      rxyz(15)=+r7(13,3)-r6(8,5)-r6(8,7)+r5(4,7)+r5(4,11)+r5(13,11)+r5(6,21)-r4(8,10) &
&             -r4(8,12)-r4(3,30)-r4(3,34)+r3(4,5)+r3(4,9)+r3(1,31)+r3(1,39)+r3(6,39) &
&             -r2(5,29)-r2(5,33)+r1(1,9)+r1(1,17)
      rxyz(16)=+r7(9,4)-r6(5,6)-r6(5,8)+r5(2,8)+r5(2,12)+r5(9,12)+r5(9,22)-r4(5,11) &
&             -r4(5,13)-r4(5,31)-r4(5,35)+r3(2,6)+r3(2,10)+r3(2,32)+r3(2,40)+r3(9,40) &
&             -r2(6,30)-r2(6,34)+r1(2,10)+r1(2,18)
      rxyz(17)=+r7(9,3)-r6(5,5)-r6(5,7)+r5(2,7)+r5(2,11)+r5(9,11)+r5(9,21)-r4(5,10) &
&             -r4(5,12)-r4(5,30)-r4(5,34)+r3(2,5)+r3(2,9)+r3(2,31)+r3(2,39)+r3(9,39) &
&             -r2(6,29)-r2(6,33)+r1(2,9)+r1(2,17)
      rxyz(18)=+r7(20,3)-r6(14,5)-r6(14,7)+r5(9,7)+r5(9,11)+r5(20,11)+r5(9,21)-r4(14,10) &
&             -r4(14,12)-r4(5,30)-r4(5,34)+r3(9,5)+r3(9,9)+r3(2,31)+r3(2,39)+r3(9,39) &
&             -r2(6,29)-r2(6,33)+r1(2,9)+r1(2,17)
      rxyz(19)=+r7(20,4)-r6(14,6)-r6(14,8)+r5(9,8)+r5(9,12)+r5(20,12)+r5(9,22)-r4(14,11) &
&             -r4(14,13)-r4(5,31)-r4(5,35)+r3(9,6)+r3(9,10)+r3(2,32)+r3(2,40)+r3(9,40) &
&             -r2(6,30)-r2(6,34)+r1(2,10)+r1(2,18)
      rxyz(20)=+r6(10,11)*four-r5(6,15)*four-r5(6,19)*four+r4(3,20)*four+r4(3,28)*four &
&             +r4(10,28)*four-r3(6,21)*four-r3(6,25)*four+r2(5,11)*four+r2(5,19)*four
      eri(1,1,5,5)=r400+(+r7(6,3)*two+r7(6,4)*two-r6(3,5)*two-r6(3,6)*two-r6(3,7)*two &
&                  -r6(3,8)*two+r5(1,7)*two+r5(1,8)*two+r5(1,11)*two+r5(6,11)*two &
&                  +r5(1,12)*two+r5(6,12)*two+r5(6,21)*six+r5(6,22)*six-r4(3,10)*two &
&                  -r4(3,11)*two-r4(3,12)*two-r4(3,13)*two-r4(3,30)*six-r4(3,31)*six &
&                  -r4(3,34)*six-r4(3,35)*six+r3(1,5)*two+r3(1,6)*two+r3(1,9)*two &
&                  +r3(1,10)*two+r3(1,31)*six+r3(1,32)*six+r3(1,39)*six+r3(6,39)*six &
&                  +r3(1,40)*six+r3(6,40)*six-r2(5,29)*six-r2(5,30)*six-r2(5,33)*six &
&                  -r2(5,34)*six+r1(1,9)*six+r1(1,10)*six+r1(1,17)*six+r1(1,18)*six)*qx+( &
&                  +r6(6,10)+r6(6,11)*four+r6(6,12)-r5(3,14)-r5(3,15)*four-r5(3,16)-r5(3,18) &
&                  -r5(3,19)*four-r5(3,20)+r4(1,19)+r4(1,20)*four+r4(1,21)+r4(1,27)+r4(6,27) &
&                  +r4(1,28)*four+r4(6,28)*four+r4(1,29)+r4(6,29)+r4(6,39)+r4(6,40)*four &
&                  +r4(6,41)-r3(3,20)-r3(3,21)*four-r3(3,22)-r3(3,24)-r3(3,25)*four-r3(3,26) &
&                  -r3(3,44)-r3(3,45)*four-r3(3,46)-r3(3,49)-r3(3,50)*four-r3(3,51)+r2(1,10) &
&                  +r2(1,11)*four+r2(1,12)+r2(1,18)+r2(1,19)*four+r2(1,20)+r2(1,43) &
&                  +r2(1,44)*four+r2(1,45)+r2(1,53)+r2(3,53)+r2(1,54)*four+r2(3,54)*four &
&                  +r2(1,55)+r2(3,55)-r1(3,32)-r1(3,33)*four-r1(3,34)-r1(3,37)-r1(3,38)*four &
&                  -r1(3,39)+r0(12)+r0(13)*four+r0(14)+r0(22)+r0(23)*four+r0(24))*xx+( &
&                  +r5(6,23)*two+r5(6,24)*two-r4(3,32)*two-r4(3,33)*two-r4(3,36)*two &
&                  -r4(3,37)*two+r3(1,33)*two+r3(1,34)*two+r3(1,41)*two+r3(6,41)*two &
&                  +r3(1,42)*two+r3(6,42)*two-r2(5,31)*two-r2(5,32)*two-r2(5,35)*two &
&                  -r2(5,36)*two+r1(1,11)*two+r1(1,12)*two+r1(1,19)*two+r1(1,20)*two)*xxx &
&                  +rxyz(1)*xxxx
      eri(2,1,5,5)=r220+(+r7(13,4)*two-r6(8,6)*two-r6(8,8)*two+r5(4,8)*two+r5(4,12)*two &
&                  +r5(13,12)*two+r5(6,22)*two-r4(8,11)*two-r4(8,13)*two-r4(3,31)*two &
&                  -r4(3,35)*two+r3(4,6)*two+r3(4,10)*two+r3(1,32)*two+r3(1,40)*two &
&                  +r3(6,40)*two-r2(5,30)*two-r2(5,34)*two+r1(1,10)*two+r1(1,18)*two)*qx &
&                  +rxyz(5)*xx
      eri(3,1,5,5)=r202+(+r7(15,4)*two-r6(10,6)*two-r6(10,8)*two+r5(6,8)*two+r5(6,12)*two &
&                  +r5(15,12)*two+r5(6,22)*two-r4(10,11)*two-r4(10,13)*two-r4(3,31)*two &
&                  -r4(3,35)*two+r3(6,6)*two+r3(6,10)*two+r3(1,32)*two+r3(1,40)*two &
&                  +r3(6,40)*two-r2(5,30)*two-r2(5,34)*two+r1(1,10)*two+r1(1,18)*two)*qx+( &
&                  +r7(10,3)*two-r6(6,5)*two-r6(6,7)*two+r5(3,7)*two+r5(3,11)*two &
&                  +r5(10,11)*two+r5(10,21)*two-r4(6,10)*two-r4(6,12)*two-r4(6,30)*two &
&                  -r4(6,34)*two+r3(3,5)*two+r3(3,9)*two+r3(3,31)*two+r3(3,39)*two &
&                  +r3(10,39)*two-r2(3,29)*two-r2(3,33)*two+r1(3,9)*two+r1(3,17)*two)*qz+( &
&                  +r6(15,12)-r5(10,16)-r5(10,20)+r4(6,21)+r4(6,29)+r4(15,29)+r4(6,41) &
&                  -r3(10,22)-r3(10,26)-r3(3,46)-r3(3,51)+r2(3,12)+r2(3,20)+r2(1,45)+r2(1,55) &
&                  +r2(3,55)-r1(3,34)-r1(3,39)+r0(14)+r0(24))*xx+rxyz(20)*xz+(+r6(6,10) &
&                  -r5(3,14)-r5(3,18)+r4(1,19)+r4(1,27)+r4(6,27)+r4(6,39)-r3(3,20)-r3(3,24) &
&                  -r3(3,44)-r3(3,49)+r2(1,10)+r2(1,18)+r2(1,43)+r2(1,53)+r2(3,53)-r1(3,32) &
&                  -r1(3,37)+r0(12)+r0(22))*zz+(+r5(10,24)*two-r4(6,33)*two-r4(6,37)*two &
&                  +r3(3,34)*two+r3(3,42)*two+r3(10,42)*two-r2(3,32)*two-r2(3,36)*two &
&                  +r1(3,12)*two+r1(3,20)*two)*xxz+(+r5(6,23)*two-r4(3,32)*two-r4(3,36)*two &
&                  +r3(1,33)*two+r3(1,41)*two+r3(6,41)*two-r2(5,31)*two-r2(5,35)*two &
&                  +r1(1,11)*two+r1(1,19)*two)*xzz+rxyz(1)*xxzz
      eri(4,1,5,5)=r310+(+r7(9,3)+r7(9,4)*two-r6(5,5)-r6(5,6)*two-r6(5,7)-r6(5,8)*two &
&                  +r5(2,7)+r5(2,8)*two+r5(2,11)+r5(9,11)+r5(2,12)*two+r5(9,12)*two+r5(9,21) &
&                  +r5(9,22)*two-r4(5,10)-r4(5,11)*two-r4(5,12)-r4(5,13)*two-r4(5,30) &
&                  -r4(5,31)*two-r4(5,34)-r4(5,35)*two+r3(2,5)+r3(2,6)*two+r3(2,9) &
&                  +r3(2,10)*two+r3(2,31)+r3(2,32)*two+r3(2,39)+r3(9,39)+r3(2,40)*two &
&                  +r3(9,40)*two-r2(6,29)-r2(6,30)*two-r2(6,33)-r2(6,34)*two+r1(2,9) &
&                  +r1(2,10)*two+r1(2,17)+r1(2,18)*two)*qx+(+r6(9,11)*two+r6(9,12) &
&                  -r5(5,15)*two-r5(5,16)-r5(5,19)*two-r5(5,20)+r4(2,20)*two+r4(2,21) &
&                  +r4(2,28)*two+r4(9,28)*two+r4(2,29)+r4(9,29)-r3(5,21)*two-r3(5,22) &
&                  -r3(5,25)*two-r3(5,26)+r2(4,11)*two+r2(4,12)+r2(4,19)*two+r2(4,20))*xx &
&                  +rxyz(2)*xxx
      eri(5,1,5,5)=r301+(+r7(10,3)+r7(10,4)*two-r6(6,5)-r6(6,6)*two-r6(6,7)-r6(6,8)*two &
&                  +r5(3,7)+r5(3,8)*two+r5(3,11)+r5(10,11)+r5(3,12)*two+r5(10,12)*two &
&                  +r5(10,21)+r5(10,22)*two-r4(6,10)-r4(6,11)*two-r4(6,12)-r4(6,13)*two &
&                  -r4(6,30)-r4(6,31)*two-r4(6,34)-r4(6,35)*two+r3(3,5)+r3(3,6)*two+r3(3,9) &
&                  +r3(3,10)*two+r3(3,31)+r3(3,32)*two+r3(3,39)+r3(10,39)+r3(3,40)*two &
&                  +r3(10,40)*two-r2(3,29)-r2(3,30)*two-r2(3,33)-r2(3,34)*two+r1(3,9) &
&                  +r1(3,10)*two+r1(3,17)+r1(3,18)*two)*qx+(+r7(6,3)-r6(3,5)-r6(3,7)+r5(1,7) &
&                  +r5(1,11)+r5(6,11)+r5(6,21)*three-r4(3,10)-r4(3,12)-r4(3,30)*three &
&                  -r4(3,34)*three+r3(1,5)+r3(1,9)+r3(1,31)*three+r3(1,39)*three &
&                  +r3(6,39)*three-r2(5,29)*three-r2(5,33)*three+r1(1,9)*three+r1(1,17)*three &
&                  )*qz+(+r6(10,11)*two+r6(10,12)-r5(6,15)*two-r5(6,16)-r5(6,19)*two-r5(6,20) &
&                  +r4(3,20)*two+r4(3,21)+r4(3,28)*two+r4(10,28)*two+r4(3,29)+r4(10,29) &
&                  -r3(6,21)*two-r3(6,22)-r3(6,25)*two-r3(6,26)+r2(5,11)*two+r2(5,12) &
&                  +r2(5,19)*two+r2(5,20))*xx+(+r6(6,10)+r6(6,11)*two-r5(3,14)-r5(3,15)*two &
&                  -r5(3,18)-r5(3,19)*two+r4(1,19)+r4(1,20)*two+r4(1,27)+r4(6,27) &
&                  +r4(1,28)*two+r4(6,28)*two+r4(6,39)+r4(6,40)*two-r3(3,20)-r3(3,21)*two &
&                  -r3(3,24)-r3(3,25)*two-r3(3,44)-r3(3,45)*two-r3(3,49)-r3(3,50)*two &
&                  +r2(1,10)+r2(1,11)*two+r2(1,18)+r2(1,19)*two+r2(1,43)+r2(1,44)*two &
&                  +r2(1,53)+r2(3,53)+r2(1,54)*two+r2(3,54)*two-r1(3,32)-r1(3,33)*two &
&                  -r1(3,37)-r1(3,38)*two+r0(12)+r0(13)*two+r0(22)+r0(23)*two)*xz+(+r5(10,24) &
&                  -r4(6,33)-r4(6,37)+r3(3,34)+r3(3,42)+r3(10,42)-r2(3,32)-r2(3,36)+r1(3,12) &
&                  +r1(3,20))*xxx+(+r5(6,23)*two+r5(6,24)-r4(3,32)*two-r4(3,33)-r4(3,36)*two &
&                  -r4(3,37)+r3(1,33)*two+r3(1,34)+r3(1,41)*two+r3(6,41)*two+r3(1,42) &
&                  +r3(6,42)-r2(5,31)*two-r2(5,32)-r2(5,35)*two-r2(5,36)+r1(1,11)*two &
&                  +r1(1,12)+r1(1,19)*two+r1(1,20))*xxz+rxyz(1)*xxxz
      eri(6,1,5,5)=r211+(+r7(14,4)*two-r6(9,6)*two-r6(9,8)*two+r5(5,8)*two+r5(5,12)*two &
&                  +r5(14,12)*two-r4(9,11)*two-r4(9,13)*two+r3(5,6)*two+r3(5,10)*two)*qx &
&                  +rxyz(17)*qz+(+r6(14,12)-r5(9,16)-r5(9,20)+r4(5,21)+r4(5,29)+r4(14,29) &
&                  -r3(9,22)-r3(9,26)+r2(6,12)+r2(6,20))*xx+(+r6(9,11)*two-r5(5,15)*two &
&                  -r5(5,19)*two+r4(2,20)*two+r4(2,28)*two+r4(9,28)*two-r3(5,21)*two &
&                  -r3(5,25)*two+r2(4,11)*two+r2(4,19)*two)*xz+rxyz(2)*xxz
      eri(1,2,5,5)=r220+(+r7(13,3)*two-r6(8,5)*two-r6(8,7)*two+r5(4,7)*two+r5(4,11)*two &
&                  +r5(13,11)*two+r5(6,21)*two-r4(8,10)*two-r4(8,12)*two-r4(3,30)*two &
&                  -r4(3,34)*two+r3(4,5)*two+r3(4,9)*two+r3(1,31)*two+r3(1,39)*two &
&                  +r3(6,39)*two-r2(5,29)*two-r2(5,33)*two+r1(1,9)*two+r1(1,17)*two)*qx &
&                  +rxyz(4)*xx
      eri(2,2,5,5)=r040
      eri(3,2,5,5)=r022+(+r7(19,3)*two-r6(13,5)*two-r6(13,7)*two+r5(8,7)*two+r5(8,11)*two &
&                  +r5(19,11)*two+r5(10,21)*two-r4(13,10)*two-r4(13,12)*two-r4(6,30)*two &
&                  -r4(6,34)*two+r3(8,5)*two+r3(8,9)*two+r3(3,31)*two+r3(3,39)*two &
&                  +r3(10,39)*two-r2(3,29)*two-r2(3,33)*two+r1(3,9)*two+r1(3,17)*two)*qz &
&                  +rxyz(4)*zz
      eri(4,2,5,5)=r130+rxyz(7)*qx
      eri(5,2,5,5)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,5,5)=r031+rxyz(7)*qz
      eri(1,3,5,5)=r202+(+r7(15,3)*two-r6(10,5)*two-r6(10,7)*two+r5(6,7)*two+r5(6,11)*two &
&                  +r5(15,11)*two+r5(6,21)*two-r4(10,10)*two-r4(10,12)*two-r4(3,30)*two &
&                  -r4(3,34)*two+r3(6,5)*two+r3(6,9)*two+r3(1,31)*two+r3(1,39)*two &
&                  +r3(6,39)*two-r2(5,29)*two-r2(5,33)*two+r1(1,9)*two+r1(1,17)*two)*qx+( &
&                  +r7(10,4)*two-r6(6,6)*two-r6(6,8)*two+r5(3,8)*two+r5(3,12)*two &
&                  +r5(10,12)*two+r5(10,22)*two-r4(6,11)*two-r4(6,13)*two-r4(6,31)*two &
&                  -r4(6,35)*two+r3(3,6)*two+r3(3,10)*two+r3(3,32)*two+r3(3,40)*two &
&                  +r3(10,40)*two-r2(3,30)*two-r2(3,34)*two+r1(3,10)*two+r1(3,18)*two)*qz+( &
&                  +r6(15,10)-r5(10,14)-r5(10,18)+r4(6,19)+r4(6,27)+r4(15,27)+r4(6,39) &
&                  -r3(10,20)-r3(10,24)-r3(3,44)-r3(3,49)+r2(3,10)+r2(3,18)+r2(1,43)+r2(1,53) &
&                  +r2(3,53)-r1(3,32)-r1(3,37)+r0(12)+r0(22))*xx+rxyz(20)*xz+(+r6(6,12) &
&                  -r5(3,16)-r5(3,20)+r4(1,21)+r4(1,29)+r4(6,29)+r4(6,41)-r3(3,22)-r3(3,26) &
&                  -r3(3,46)-r3(3,51)+r2(1,12)+r2(1,20)+r2(1,45)+r2(1,55)+r2(3,55)-r1(3,34) &
&                  -r1(3,39)+r0(14)+r0(24))*zz+(+r5(10,23)*two-r4(6,32)*two-r4(6,36)*two &
&                  +r3(3,33)*two+r3(3,41)*two+r3(10,41)*two-r2(3,31)*two-r2(3,35)*two &
&                  +r1(3,11)*two+r1(3,19)*two)*xxz+(+r5(6,24)*two-r4(3,33)*two-r4(3,37)*two &
&                  +r3(1,34)*two+r3(1,42)*two+r3(6,42)*two-r2(5,32)*two-r2(5,36)*two &
&                  +r1(1,12)*two+r1(1,20)*two)*xzz+rxyz(1)*xxzz
      eri(2,3,5,5)=r022+(+r7(19,4)*two-r6(13,6)*two-r6(13,8)*two+r5(8,8)*two+r5(8,12)*two &
&                  +r5(19,12)*two+r5(10,22)*two-r4(13,11)*two-r4(13,13)*two-r4(6,31)*two &
&                  -r4(6,35)*two+r3(8,6)*two+r3(8,10)*two+r3(3,32)*two+r3(3,40)*two &
&                  +r3(10,40)*two-r2(3,30)*two-r2(3,34)*two+r1(3,10)*two+r1(3,18)*two)*qz &
&                  +rxyz(5)*zz
      eri(3,3,5,5)=r004+(+r7(21,3)*two+r7(21,4)*two-r6(15,5)*two-r6(15,6)*two &
&                  -r6(15,7)*two-r6(15,8)*two+r5(10,7)*two+r5(10,8)*two+r5(10,11)*two &
&                  +r5(21,11)*two+r5(10,12)*two+r5(21,12)*two+r5(10,21)*six+r5(10,22)*six &
&                  -r4(15,10)*two-r4(15,11)*two-r4(15,12)*two-r4(15,13)*two-r4(6,30)*six &
&                  -r4(6,31)*six-r4(6,34)*six-r4(6,35)*six+r3(10,5)*two+r3(10,6)*two &
&                  +r3(10,9)*two+r3(10,10)*two+r3(3,31)*six+r3(3,32)*six+r3(3,39)*six &
&                  +r3(10,39)*six+r3(3,40)*six+r3(10,40)*six-r2(3,29)*six-r2(3,30)*six &
&                  -r2(3,33)*six-r2(3,34)*six+r1(3,9)*six+r1(3,10)*six+r1(3,17)*six &
&                  +r1(3,18)*six)*qz+(+r6(15,10)+r6(15,11)*four+r6(15,12)-r5(10,14) &
&                  -r5(10,15)*four-r5(10,16)-r5(10,18)-r5(10,19)*four-r5(10,20)+r4(6,19) &
&                  +r4(6,20)*four+r4(6,21)+r4(6,27)+r4(15,27)+r4(6,28)*four+r4(15,28)*four &
&                  +r4(6,29)+r4(15,29)+r4(6,39)+r4(6,40)*four+r4(6,41)-r3(10,20) &
&                  -r3(10,21)*four-r3(10,22)-r3(10,24)-r3(10,25)*four-r3(10,26)-r3(3,44) &
&                  -r3(3,45)*four-r3(3,46)-r3(3,49)-r3(3,50)*four-r3(3,51)+r2(3,10) &
&                  +r2(3,11)*four+r2(3,12)+r2(3,18)+r2(3,19)*four+r2(3,20)+r2(1,43) &
&                  +r2(1,44)*four+r2(1,45)+r2(1,53)+r2(3,53)+r2(1,54)*four+r2(3,54)*four &
&                  +r2(1,55)+r2(3,55)-r1(3,32)-r1(3,33)*four-r1(3,34)-r1(3,37)-r1(3,38)*four &
&                  -r1(3,39)+r0(12)+r0(13)*four+r0(14)+r0(22)+r0(23)*four+r0(24))*zz+( &
&                  +r5(10,23)*two+r5(10,24)*two-r4(6,32)*two-r4(6,33)*two-r4(6,36)*two &
&                  -r4(6,37)*two+r3(3,33)*two+r3(3,34)*two+r3(3,41)*two+r3(10,41)*two &
&                  +r3(3,42)*two+r3(10,42)*two-r2(3,31)*two-r2(3,32)*two-r2(3,35)*two &
&                  -r2(3,36)*two+r1(3,11)*two+r1(3,12)*two+r1(3,19)*two+r1(3,20)*two)*zzz &
&                  +rxyz(1)*zzzz
      eri(4,3,5,5)=r112+rxyz(18)*qx+(+r7(14,4)*two-r6(9,6)*two-r6(9,8)*two+r5(5,8)*two &
&                  +r5(5,12)*two+r5(14,12)*two-r4(9,11)*two-r4(9,13)*two+r3(5,6)*two &
&                  +r3(5,10)*two)*qz+(+r6(14,11)*two-r5(9,15)*two-r5(9,19)*two+r4(5,20)*two &
&                  +r4(5,28)*two+r4(14,28)*two-r3(9,21)*two-r3(9,25)*two+r2(6,11)*two &
&                  +r2(6,19)*two)*xz+(+r6(9,12)-r5(5,16)-r5(5,20)+r4(2,21)+r4(2,29)+r4(9,29) &
&                  -r3(5,22)-r3(5,26)+r2(4,12)+r2(4,20))*zz+rxyz(2)*xzz
      eri(5,3,5,5)=r103+(+r7(21,3)-r6(15,5)-r6(15,7)+r5(10,7)+r5(10,11)+r5(21,11) &
&                  +r5(10,21)*three-r4(15,10)-r4(15,12)-r4(6,30)*three-r4(6,34)*three &
&                  +r3(10,5)+r3(10,9)+r3(3,31)*three+r3(3,39)*three+r3(10,39)*three &
&                  -r2(3,29)*three-r2(3,33)*three+r1(3,9)*three+r1(3,17)*three)*qx+(+r7(15,3) &
&                  +r7(15,4)*two-r6(10,5)-r6(10,6)*two-r6(10,7)-r6(10,8)*two+r5(6,7) &
&                  +r5(6,8)*two+r5(6,11)+r5(15,11)+r5(6,12)*two+r5(15,12)*two+r5(6,21) &
&                  +r5(6,22)*two-r4(10,10)-r4(10,11)*two-r4(10,12)-r4(10,13)*two-r4(3,30) &
&                  -r4(3,31)*two-r4(3,34)-r4(3,35)*two+r3(6,5)+r3(6,6)*two+r3(6,9) &
&                  +r3(6,10)*two+r3(1,31)+r3(1,32)*two+r3(1,39)+r3(6,39)+r3(1,40)*two &
&                  +r3(6,40)*two-r2(5,29)-r2(5,30)*two-r2(5,33)-r2(5,34)*two+r1(1,9) &
&                  +r1(1,10)*two+r1(1,17)+r1(1,18)*two)*qz+(+r6(15,10)+r6(15,11)*two &
&                  -r5(10,14)-r5(10,15)*two-r5(10,18)-r5(10,19)*two+r4(6,19)+r4(6,20)*two &
&                  +r4(6,27)+r4(15,27)+r4(6,28)*two+r4(15,28)*two+r4(6,39)+r4(6,40)*two &
&                  -r3(10,20)-r3(10,21)*two-r3(10,24)-r3(10,25)*two-r3(3,44)-r3(3,45)*two &
&                  -r3(3,49)-r3(3,50)*two+r2(3,10)+r2(3,11)*two+r2(3,18)+r2(3,19)*two &
&                  +r2(1,43)+r2(1,44)*two+r2(1,53)+r2(3,53)+r2(1,54)*two+r2(3,54)*two &
&                  -r1(3,32)-r1(3,33)*two-r1(3,37)-r1(3,38)*two+r0(12)+r0(13)*two+r0(22) &
&                  +r0(23)*two)*xz+(+r6(10,11)*two+r6(10,12)-r5(6,15)*two-r5(6,16) &
&                  -r5(6,19)*two-r5(6,20)+r4(3,20)*two+r4(3,21)+r4(3,28)*two+r4(10,28)*two &
&                  +r4(3,29)+r4(10,29)-r3(6,21)*two-r3(6,22)-r3(6,25)*two-r3(6,26) &
&                  +r2(5,11)*two+r2(5,12)+r2(5,19)*two+r2(5,20))*zz+(+r5(10,23)*two+r5(10,24) &
&                  -r4(6,32)*two-r4(6,33)-r4(6,36)*two-r4(6,37)+r3(3,33)*two+r3(3,34) &
&                  +r3(3,41)*two+r3(10,41)*two+r3(3,42)+r3(10,42)-r2(3,31)*two-r2(3,32) &
&                  -r2(3,35)*two-r2(3,36)+r1(3,11)*two+r1(3,12)+r1(3,19)*two+r1(3,20))*xzz+( &
&                  +r5(6,24)-r4(3,33)-r4(3,37)+r3(1,34)+r3(1,42)+r3(6,42)-r2(5,32)-r2(5,36) &
&                  +r1(1,12)+r1(1,20))*zzz+rxyz(1)*xzzz
      eri(6,3,5,5)=r013+(+r7(20,3)+r7(20,4)*two-r6(14,5)-r6(14,6)*two-r6(14,7) &
&                  -r6(14,8)*two+r5(9,7)+r5(9,8)*two+r5(9,11)+r5(20,11)+r5(9,12)*two &
&                  +r5(20,12)*two+r5(9,21)+r5(9,22)*two-r4(14,10)-r4(14,11)*two-r4(14,12) &
&                  -r4(14,13)*two-r4(5,30)-r4(5,31)*two-r4(5,34)-r4(5,35)*two+r3(9,5) &
&                  +r3(9,6)*two+r3(9,9)+r3(9,10)*two+r3(2,31)+r3(2,32)*two+r3(2,39)+r3(9,39) &
&                  +r3(2,40)*two+r3(9,40)*two-r2(6,29)-r2(6,30)*two-r2(6,33)-r2(6,34)*two &
&                  +r1(2,9)+r1(2,10)*two+r1(2,17)+r1(2,18)*two)*qz+(+r6(14,11)*two+r6(14,12) &
&                  -r5(9,15)*two-r5(9,16)-r5(9,19)*two-r5(9,20)+r4(5,20)*two+r4(5,21) &
&                  +r4(5,28)*two+r4(14,28)*two+r4(5,29)+r4(14,29)-r3(9,21)*two-r3(9,22) &
&                  -r3(9,25)*two-r3(9,26)+r2(6,11)*two+r2(6,12)+r2(6,19)*two+r2(6,20))*zz &
&                  +rxyz(2)*zzz
      eri(1,4,5,5)=r310+(+r7(9,3)*two+r7(9,4)-r6(5,5)*two-r6(5,6)-r6(5,7)*two-r6(5,8) &
&                  +r5(2,7)*two+r5(2,8)+r5(2,11)*two+r5(9,11)*two+r5(2,12)+r5(9,12) &
&                  +r5(9,21)*two+r5(9,22)-r4(5,10)*two-r4(5,11)-r4(5,12)*two-r4(5,13) &
&                  -r4(5,30)*two-r4(5,31)-r4(5,34)*two-r4(5,35)+r3(2,5)*two+r3(2,6) &
&                  +r3(2,9)*two+r3(2,10)+r3(2,31)*two+r3(2,32)+r3(2,39)*two+r3(9,39)*two &
&                  +r3(2,40)+r3(9,40)-r2(6,29)*two-r2(6,30)-r2(6,33)*two-r2(6,34)+r1(2,9)*two &
&                  +r1(2,10)+r1(2,17)*two+r1(2,18))*qx+(+r6(9,10)+r6(9,11)*two-r5(5,14) &
&                  -r5(5,15)*two-r5(5,18)-r5(5,19)*two+r4(2,19)+r4(2,20)*two+r4(2,27) &
&                  +r4(9,27)+r4(2,28)*two+r4(9,28)*two-r3(5,20)-r3(5,21)*two-r3(5,24) &
&                  -r3(5,25)*two+r2(4,10)+r2(4,11)*two+r2(4,18)+r2(4,19)*two)*xx+rxyz(3)*xxx
      eri(2,4,5,5)=r130+rxyz(8)*qx
      eri(3,4,5,5)=r112+rxyz(19)*qx+(+r7(14,3)*two-r6(9,5)*two-r6(9,7)*two+r5(5,7)*two &
&                  +r5(5,11)*two+r5(14,11)*two-r4(9,10)*two-r4(9,12)*two+r3(5,5)*two &
&                  +r3(5,9)*two)*qz+(+r6(14,11)*two-r5(9,15)*two-r5(9,19)*two+r4(5,20)*two &
&                  +r4(5,28)*two+r4(14,28)*two-r3(9,21)*two-r3(9,25)*two+r2(6,11)*two &
&                  +r2(6,19)*two)*xz+(+r6(9,10)-r5(5,14)-r5(5,18)+r4(2,19)+r4(2,27)+r4(9,27) &
&                  -r3(5,20)-r3(5,24)+r2(4,10)+r2(4,18))*zz+rxyz(3)*xzz
      eri(4,4,5,5)=r220+(+r7(13,3)+r7(13,4)-r6(8,5)-r6(8,6)-r6(8,7)-r6(8,8)+r5(4,7) &
&                  +r5(4,8)+r5(4,11)+r5(13,11)+r5(4,12)+r5(13,12)+r5(6,21)+r5(6,22)-r4(8,10) &
&                  -r4(8,11)-r4(8,12)-r4(8,13)-r4(3,30)-r4(3,31)-r4(3,34)-r4(3,35)+r3(4,5) &
&                  +r3(4,6)+r3(4,9)+r3(4,10)+r3(1,31)+r3(1,32)+r3(1,39)+r3(6,39)+r3(1,40) &
&                  +r3(6,40)-r2(5,29)-r2(5,30)-r2(5,33)-r2(5,34)+r1(1,9)+r1(1,10)+r1(1,17) &
&                  +r1(1,18))*qx+rxyz(6)*xx
      eri(5,4,5,5)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(9,10)+r6(9,11)-r5(5,14) &
&                  -r5(5,15)-r5(5,18)-r5(5,19)+r4(2,19)+r4(2,20)+r4(2,27)+r4(9,27)+r4(2,28) &
&                  +r4(9,28)-r3(5,20)-r3(5,21)-r3(5,24)-r3(5,25)+r2(4,10)+r2(4,11)+r2(4,18) &
&                  +r2(4,19))*xz+rxyz(3)*xxz
      eri(6,4,5,5)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,5,5)=r301+(+r7(10,3)*two+r7(10,4)-r6(6,5)*two-r6(6,6)-r6(6,7)*two-r6(6,8) &
&                  +r5(3,7)*two+r5(3,8)+r5(3,11)*two+r5(10,11)*two+r5(3,12)+r5(10,12) &
&                  +r5(10,21)*two+r5(10,22)-r4(6,10)*two-r4(6,11)-r4(6,12)*two-r4(6,13) &
&                  -r4(6,30)*two-r4(6,31)-r4(6,34)*two-r4(6,35)+r3(3,5)*two+r3(3,6) &
&                  +r3(3,9)*two+r3(3,10)+r3(3,31)*two+r3(3,32)+r3(3,39)*two+r3(10,39)*two &
&                  +r3(3,40)+r3(10,40)-r2(3,29)*two-r2(3,30)-r2(3,33)*two-r2(3,34) &
&                  +r1(3,9)*two+r1(3,10)+r1(3,17)*two+r1(3,18))*qx+(+r7(6,4)-r6(3,6)-r6(3,8) &
&                  +r5(1,8)+r5(1,12)+r5(6,12)+r5(6,22)*three-r4(3,11)-r4(3,13)-r4(3,31)*three &
&                  -r4(3,35)*three+r3(1,6)+r3(1,10)+r3(1,32)*three+r3(1,40)*three &
&                  +r3(6,40)*three-r2(5,30)*three-r2(5,34)*three+r1(1,10)*three &
&                  +r1(1,18)*three)*qz+(+r6(10,10)+r6(10,11)*two-r5(6,14)-r5(6,15)*two &
&                  -r5(6,18)-r5(6,19)*two+r4(3,19)+r4(3,20)*two+r4(3,27)+r4(10,27) &
&                  +r4(3,28)*two+r4(10,28)*two-r3(6,20)-r3(6,21)*two-r3(6,24)-r3(6,25)*two &
&                  +r2(5,10)+r2(5,11)*two+r2(5,18)+r2(5,19)*two)*xx+(+r6(6,11)*two+r6(6,12) &
&                  -r5(3,15)*two-r5(3,16)-r5(3,19)*two-r5(3,20)+r4(1,20)*two+r4(1,21) &
&                  +r4(1,28)*two+r4(6,28)*two+r4(1,29)+r4(6,29)+r4(6,40)*two+r4(6,41) &
&                  -r3(3,21)*two-r3(3,22)-r3(3,25)*two-r3(3,26)-r3(3,45)*two-r3(3,46) &
&                  -r3(3,50)*two-r3(3,51)+r2(1,11)*two+r2(1,12)+r2(1,19)*two+r2(1,20) &
&                  +r2(1,44)*two+r2(1,45)+r2(1,54)*two+r2(3,54)*two+r2(1,55)+r2(3,55) &
&                  -r1(3,33)*two-r1(3,34)-r1(3,38)*two-r1(3,39)+r0(13)*two+r0(14)+r0(23)*two &
&                  +r0(24))*xz+(+r5(10,23)-r4(6,32)-r4(6,36)+r3(3,33)+r3(3,41)+r3(10,41) &
&                  -r2(3,31)-r2(3,35)+r1(3,11)+r1(3,19))*xxx+(+r5(6,23)+r5(6,24)*two-r4(3,32) &
&                  -r4(3,33)*two-r4(3,36)-r4(3,37)*two+r3(1,33)+r3(1,34)*two+r3(1,41) &
&                  +r3(6,41)+r3(1,42)*two+r3(6,42)*two-r2(5,31)-r2(5,32)*two-r2(5,35) &
&                  -r2(5,36)*two+r1(1,11)+r1(1,12)*two+r1(1,19)+r1(1,20)*two)*xxz+rxyz(1) &
&                  *xxxz
      eri(2,5,5,5)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,5,5)=r103+(+r7(21,4)-r6(15,6)-r6(15,8)+r5(10,8)+r5(10,12)+r5(21,12) &
&                  +r5(10,22)*three-r4(15,11)-r4(15,13)-r4(6,31)*three-r4(6,35)*three &
&                  +r3(10,6)+r3(10,10)+r3(3,32)*three+r3(3,40)*three+r3(10,40)*three &
&                  -r2(3,30)*three-r2(3,34)*three+r1(3,10)*three+r1(3,18)*three)*qx+( &
&                  +r7(15,3)*two+r7(15,4)-r6(10,5)*two-r6(10,6)-r6(10,7)*two-r6(10,8) &
&                  +r5(6,7)*two+r5(6,8)+r5(6,11)*two+r5(15,11)*two+r5(6,12)+r5(15,12) &
&                  +r5(6,21)*two+r5(6,22)-r4(10,10)*two-r4(10,11)-r4(10,12)*two-r4(10,13) &
&                  -r4(3,30)*two-r4(3,31)-r4(3,34)*two-r4(3,35)+r3(6,5)*two+r3(6,6) &
&                  +r3(6,9)*two+r3(6,10)+r3(1,31)*two+r3(1,32)+r3(1,39)*two+r3(6,39)*two &
&                  +r3(1,40)+r3(6,40)-r2(5,29)*two-r2(5,30)-r2(5,33)*two-r2(5,34)+r1(1,9)*two &
&                  +r1(1,10)+r1(1,17)*two+r1(1,18))*qz+(+r6(15,11)*two+r6(15,12) &
&                  -r5(10,15)*two-r5(10,16)-r5(10,19)*two-r5(10,20)+r4(6,20)*two+r4(6,21) &
&                  +r4(6,28)*two+r4(15,28)*two+r4(6,29)+r4(15,29)+r4(6,40)*two+r4(6,41) &
&                  -r3(10,21)*two-r3(10,22)-r3(10,25)*two-r3(10,26)-r3(3,45)*two-r3(3,46) &
&                  -r3(3,50)*two-r3(3,51)+r2(3,11)*two+r2(3,12)+r2(3,19)*two+r2(3,20) &
&                  +r2(1,44)*two+r2(1,45)+r2(1,54)*two+r2(3,54)*two+r2(1,55)+r2(3,55) &
&                  -r1(3,33)*two-r1(3,34)-r1(3,38)*two-r1(3,39)+r0(13)*two+r0(14)+r0(23)*two &
&                  +r0(24))*xz+(+r6(10,10)+r6(10,11)*two-r5(6,14)-r5(6,15)*two-r5(6,18) &
&                  -r5(6,19)*two+r4(3,19)+r4(3,20)*two+r4(3,27)+r4(10,27)+r4(3,28)*two &
&                  +r4(10,28)*two-r3(6,20)-r3(6,21)*two-r3(6,24)-r3(6,25)*two+r2(5,10) &
&                  +r2(5,11)*two+r2(5,18)+r2(5,19)*two)*zz+(+r5(10,23)+r5(10,24)*two-r4(6,32) &
&                  -r4(6,33)*two-r4(6,36)-r4(6,37)*two+r3(3,33)+r3(3,34)*two+r3(3,41) &
&                  +r3(10,41)+r3(3,42)*two+r3(10,42)*two-r2(3,31)-r2(3,32)*two-r2(3,35) &
&                  -r2(3,36)*two+r1(3,11)+r1(3,12)*two+r1(3,19)+r1(3,20)*two)*xzz+(+r5(6,23) &
&                  -r4(3,32)-r4(3,36)+r3(1,33)+r3(1,41)+r3(6,41)-r2(5,31)-r2(5,35)+r1(1,11) &
&                  +r1(1,19))*zzz+rxyz(1)*xzzz
      eri(4,5,5,5)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(9,11)+r6(9,12)-r5(5,15) &
&                  -r5(5,16)-r5(5,19)-r5(5,20)+r4(2,20)+r4(2,21)+r4(2,28)+r4(9,28)+r4(2,29) &
&                  +r4(9,29)-r3(5,21)-r3(5,22)-r3(5,25)-r3(5,26)+r2(4,11)+r2(4,12)+r2(4,19) &
&                  +r2(4,20))*xz+rxyz(2)*xxz
      eri(5,5,5,5)=r202+(+r7(15,3)+r7(15,4)-r6(10,5)-r6(10,6)-r6(10,7)-r6(10,8)+r5(6,7) &
&                  +r5(6,8)+r5(6,11)+r5(15,11)+r5(6,12)+r5(15,12)+r5(6,21)+r5(6,22)-r4(10,10) &
&                  -r4(10,11)-r4(10,12)-r4(10,13)-r4(3,30)-r4(3,31)-r4(3,34)-r4(3,35)+r3(6,5) &
&                  +r3(6,6)+r3(6,9)+r3(6,10)+r3(1,31)+r3(1,32)+r3(1,39)+r3(6,39)+r3(1,40) &
&                  +r3(6,40)-r2(5,29)-r2(5,30)-r2(5,33)-r2(5,34)+r1(1,9)+r1(1,10)+r1(1,17) &
&                  +r1(1,18))*qx+(+r7(10,3)+r7(10,4)-r6(6,5)-r6(6,6)-r6(6,7)-r6(6,8)+r5(3,7) &
&                  +r5(3,8)+r5(3,11)+r5(10,11)+r5(3,12)+r5(10,12)+r5(10,21)+r5(10,22) &
&                  -r4(6,10)-r4(6,11)-r4(6,12)-r4(6,13)-r4(6,30)-r4(6,31)-r4(6,34)-r4(6,35) &
&                  +r3(3,5)+r3(3,6)+r3(3,9)+r3(3,10)+r3(3,31)+r3(3,32)+r3(3,39)+r3(10,39) &
&                  +r3(3,40)+r3(10,40)-r2(3,29)-r2(3,30)-r2(3,33)-r2(3,34)+r1(3,9)+r1(3,10) &
&                  +r1(3,17)+r1(3,18))*qz+(+r6(15,11)-r5(10,15)-r5(10,19)+r4(6,20)+r4(6,28) &
&                  +r4(15,28)+r4(6,40)-r3(10,21)-r3(10,25)-r3(3,45)-r3(3,50)+r2(3,11) &
&                  +r2(3,19)+r2(1,44)+r2(1,54)+r2(3,54)-r1(3,33)-r1(3,38)+r0(13)+r0(23))*xx+( &
&                  +r6(10,10)+r6(10,11)*two+r6(10,12)-r5(6,14)-r5(6,15)*two-r5(6,16)-r5(6,18) &
&                  -r5(6,19)*two-r5(6,20)+r4(3,19)+r4(3,20)*two+r4(3,21)+r4(3,27)+r4(10,27) &
&                  +r4(3,28)*two+r4(10,28)*two+r4(3,29)+r4(10,29)-r3(6,20)-r3(6,21)*two &
&                  -r3(6,22)-r3(6,24)-r3(6,25)*two-r3(6,26)+r2(5,10)+r2(5,11)*two+r2(5,12) &
&                  +r2(5,18)+r2(5,19)*two+r2(5,20))*xz+(+r6(6,11)-r5(3,15)-r5(3,19)+r4(1,20) &
&                  +r4(1,28)+r4(6,28)+r4(6,40)-r3(3,21)-r3(3,25)-r3(3,45)-r3(3,50)+r2(1,11) &
&                  +r2(1,19)+r2(1,44)+r2(1,54)+r2(3,54)-r1(3,33)-r1(3,38)+r0(13)+r0(23))*zz+( &
&                  +r5(10,23)+r5(10,24)-r4(6,32)-r4(6,33)-r4(6,36)-r4(6,37)+r3(3,33)+r3(3,34) &
&                  +r3(3,41)+r3(10,41)+r3(3,42)+r3(10,42)-r2(3,31)-r2(3,32)-r2(3,35)-r2(3,36) &
&                  +r1(3,11)+r1(3,12)+r1(3,19)+r1(3,20))*xxz+(+r5(6,23)+r5(6,24)-r4(3,32) &
&                  -r4(3,33)-r4(3,36)-r4(3,37)+r3(1,33)+r3(1,34)+r3(1,41)+r3(6,41)+r3(1,42) &
&                  +r3(6,42)-r2(5,31)-r2(5,32)-r2(5,35)-r2(5,36)+r1(1,11)+r1(1,12)+r1(1,19) &
&                  +r1(1,20))*xzz+rxyz(1)*xxzz
      eri(6,5,5,5)=r112+rxyz(19)*qx+(+r7(14,3)+r7(14,4)-r6(9,5)-r6(9,6)-r6(9,7)-r6(9,8) &
&                  +r5(5,7)+r5(5,8)+r5(5,11)+r5(14,11)+r5(5,12)+r5(14,12)-r4(9,10)-r4(9,11) &
&                  -r4(9,12)-r4(9,13)+r3(5,5)+r3(5,6)+r3(5,9)+r3(5,10))*qz+(+r6(14,11) &
&                  +r6(14,12)-r5(9,15)-r5(9,16)-r5(9,19)-r5(9,20)+r4(5,20)+r4(5,21)+r4(5,28) &
&                  +r4(14,28)+r4(5,29)+r4(14,29)-r3(9,21)-r3(9,22)-r3(9,25)-r3(9,26)+r2(6,11) &
&                  +r2(6,12)+r2(6,19)+r2(6,20))*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,5,5)=r211+(+r7(14,3)*two-r6(9,5)*two-r6(9,7)*two+r5(5,7)*two+r5(5,11)*two &
&                  +r5(14,11)*two-r4(9,10)*two-r4(9,12)*two+r3(5,5)*two+r3(5,9)*two)*qx &
&                  +rxyz(16)*qz+(+r6(14,10)-r5(9,14)-r5(9,18)+r4(5,19)+r4(5,27)+r4(14,27) &
&                  -r3(9,20)-r3(9,24)+r2(6,10)+r2(6,18))*xx+(+r6(9,11)*two-r5(5,15)*two &
&                  -r5(5,19)*two+r4(2,20)*two+r4(2,28)*two+r4(9,28)*two-r3(5,21)*two &
&                  -r3(5,25)*two+r2(4,11)*two+r2(4,19)*two)*xz+rxyz(3)*xxz
      eri(2,6,5,5)=r031+rxyz(8)*qz
      eri(3,6,5,5)=r013+(+r7(20,3)*two+r7(20,4)-r6(14,5)*two-r6(14,6)-r6(14,7)*two &
&                  -r6(14,8)+r5(9,7)*two+r5(9,8)+r5(9,11)*two+r5(20,11)*two+r5(9,12) &
&                  +r5(20,12)+r5(9,21)*two+r5(9,22)-r4(14,10)*two-r4(14,11)-r4(14,12)*two &
&                  -r4(14,13)-r4(5,30)*two-r4(5,31)-r4(5,34)*two-r4(5,35)+r3(9,5)*two+r3(9,6) &
&                  +r3(9,9)*two+r3(9,10)+r3(2,31)*two+r3(2,32)+r3(2,39)*two+r3(9,39)*two &
&                  +r3(2,40)+r3(9,40)-r2(6,29)*two-r2(6,30)-r2(6,33)*two-r2(6,34)+r1(2,9)*two &
&                  +r1(2,10)+r1(2,17)*two+r1(2,18))*qz+(+r6(14,10)+r6(14,11)*two-r5(9,14) &
&                  -r5(9,15)*two-r5(9,18)-r5(9,19)*two+r4(5,19)+r4(5,20)*two+r4(5,27) &
&                  +r4(14,27)+r4(5,28)*two+r4(14,28)*two-r3(9,20)-r3(9,21)*two-r3(9,24) &
&                  -r3(9,25)*two+r2(6,10)+r2(6,11)*two+r2(6,18)+r2(6,19)*two)*zz+rxyz(3)*zzz
      eri(4,6,5,5)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,5,5)=r112+rxyz(18)*qx+(+r7(14,3)+r7(14,4)-r6(9,5)-r6(9,6)-r6(9,7)-r6(9,8) &
&                  +r5(5,7)+r5(5,8)+r5(5,11)+r5(14,11)+r5(5,12)+r5(14,12)-r4(9,10)-r4(9,11) &
&                  -r4(9,12)-r4(9,13)+r3(5,5)+r3(5,6)+r3(5,9)+r3(5,10))*qz+(+r6(14,10) &
&                  +r6(14,11)-r5(9,14)-r5(9,15)-r5(9,18)-r5(9,19)+r4(5,19)+r4(5,20)+r4(5,27) &
&                  +r4(14,27)+r4(5,28)+r4(14,28)-r3(9,20)-r3(9,21)-r3(9,24)-r3(9,25)+r2(6,10) &
&                  +r2(6,11)+r2(6,18)+r2(6,19))*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,5,5)=r022+(+r7(19,3)+r7(19,4)-r6(13,5)-r6(13,6)-r6(13,7)-r6(13,8)+r5(8,7) &
&                  +r5(8,8)+r5(8,11)+r5(19,11)+r5(8,12)+r5(19,12)+r5(10,21)+r5(10,22) &
&                  -r4(13,10)-r4(13,11)-r4(13,12)-r4(13,13)-r4(6,30)-r4(6,31)-r4(6,34) &
&                  -r4(6,35)+r3(8,5)+r3(8,6)+r3(8,9)+r3(8,10)+r3(3,31)+r3(3,32)+r3(3,39) &
&                  +r3(10,39)+r3(3,40)+r3(10,40)-r2(3,29)-r2(3,30)-r2(3,33)-r2(3,34)+r1(3,9) &
&                  +r1(3,10)+r1(3,17)+r1(3,18))*qz+rxyz(6)*zz
!
      r400= r8(9)-r7(5,1)-r7(5,2)+r6(2,2)+r6(2,4)+r6(9,9)*six-r5(5,13)*six-r5(5,17)*six &
&          +r4(2,18)*six+r4(2,26)*six+r4(9,38)*three-r3(5,43)*three-r3(5,48)*three &
&          +r2(4,42)*three+r2(4,52)*three
      r310= r8(13)-r7(8,1)-r7(8,2)+r6(4,2)+r6(4,4)+r6(13,9)*three-r5(8,13)*three &
&          -r5(8,17)*three+r4(4,18)*three+r4(4,26)*three
      r301= r8(14)-r7(9,1)-r7(9,2)+r6(5,2)+r6(5,4)+r6(14,9)*three-r5(9,13)*three &
&          -r5(9,17)*three+r4(5,18)*three+r4(5,26)*three
      r220= r8(18)-r7(12,1)-r7(12,2)+r6(7,2)+r6(7,4)+r6(9,9)+r6(18,9)-r5(5,13)-r5(12,13) &
&          -r5(5,17)-r5(12,17)+r4(2,18)+r4(7,18)+r4(2,26)+r4(7,26)+r4(9,38)-r3(5,43) &
&          -r3(5,48)+r2(4,42)+r2(4,52)
      r211= r8(19)-r7(13,1)-r7(13,2)+r6(8,2)+r6(8,4)+r6(19,9)-r5(13,13)-r5(13,17)+r4(8,18) &
&          +r4(8,26)
      r202= r8(20)-r7(14,1)-r7(14,2)+r6(9,2)+r6(9,4)+r6(9,9)+r6(20,9)-r5(5,13)-r5(14,13) &
&          -r5(5,17)-r5(14,17)+r4(2,18)+r4(9,18)+r4(2,26)+r4(9,26)+r4(9,38)-r3(5,43) &
&          -r3(5,48)+r2(4,42)+r2(4,52)
      r130= r8(24)-r7(17,1)-r7(17,2)+r6(11,2)+r6(11,4)+r6(13,9)*three-r5(8,13)*three &
&          -r5(8,17)*three+r4(4,18)*three+r4(4,26)*three
      r121= r8(25)-r7(18,1)-r7(18,2)+r6(12,2)+r6(12,4)+r6(14,9)-r5(9,13)-r5(9,17)+r4(5,18) &
&          +r4(5,26)
      r112= r8(26)-r7(19,1)-r7(19,2)+r6(13,2)+r6(13,4)+r6(13,9)-r5(8,13)-r5(8,17)+r4(4,18) &
&          +r4(4,26)
      r103= r8(27)-r7(20,1)-r7(20,2)+r6(14,2)+r6(14,4)+r6(14,9)*three-r5(9,13)*three &
&          -r5(9,17)*three+r4(5,18)*three+r4(5,26)*three
      r040= r8(31)-r7(23,1)-r7(23,2)+r6(16,2)+r6(16,4)+r6(18,9)*six-r5(12,13)*six &
&          -r5(12,17)*six+r4(7,18)*six+r4(7,26)*six+r4(9,38)*three-r3(5,43)*three &
&          -r3(5,48)*three+r2(4,42)*three+r2(4,52)*three
      r031= r8(32)-r7(24,1)-r7(24,2)+r6(17,2)+r6(17,4)+r6(19,9)*three-r5(13,13)*three &
&          -r5(13,17)*three+r4(8,18)*three+r4(8,26)*three
      r022= r8(33)-r7(25,1)-r7(25,2)+r6(18,2)+r6(18,4)+r6(18,9)+r6(20,9)-r5(12,13) &
&          -r5(14,13)-r5(12,17)-r5(14,17)+r4(7,18)+r4(9,18)+r4(7,26)+r4(9,26)+r4(9,38) &
&          -r3(5,43)-r3(5,48)+r2(4,42)+r2(4,52)
      r013= r8(34)-r7(26,1)-r7(26,2)+r6(19,2)+r6(19,4)+r6(19,9)*three-r5(13,13)*three &
&          -r5(13,17)*three+r4(8,18)*three+r4(8,26)*three
      r004= r8(35)-r7(27,1)-r7(27,2)+r6(20,2)+r6(20,4)+r6(20,9)*six-r5(14,13)*six &
&          -r5(14,17)*six+r4(9,18)*six+r4(9,26)*six+r4(9,38)*three-r3(5,43)*three &
&          -r3(5,48)*three+r2(4,42)*three+r2(4,52)*three
      rxyz(1)=+r4(9,42)-r3(5,47)-r3(5,52)+r2(4,46)+r2(4,56)
      rxyz(2)=+r5(13,24)-r4(8,33)-r4(8,37)+r3(4,34)+r3(4,42)
      rxyz(3)=+r5(13,23)-r4(8,32)-r4(8,36)+r3(4,33)+r3(4,41)
      rxyz(4)=+r6(18,10)-r5(12,14)-r5(12,18)+r4(7,19)+r4(7,27)+r4(9,39)-r3(5,44)-r3(5,49) &
&             +r2(4,43)+r2(4,53)
      rxyz(5)=+r6(18,12)-r5(12,16)-r5(12,20)+r4(7,21)+r4(7,29)+r4(9,41)-r3(5,46)-r3(5,51) &
&             +r2(4,45)+r2(4,55)
      rxyz(6)=+r6(18,11)-r5(12,15)-r5(12,19)+r4(7,20)+r4(7,28)+r4(9,40)-r3(5,45)-r3(5,50) &
&             +r2(4,44)+r2(4,54)
      rxyz(7)=+r7(24,3)-r6(17,5)-r6(17,7)+r5(11,7)+r5(11,11)+r5(13,21)*three &
&             -r4(8,30)*three-r4(8,34)*three+r3(4,31)*three+r3(4,39)*three
      rxyz(8)=+r7(24,4)-r6(17,6)-r6(17,8)+r5(11,8)+r5(11,12)+r5(13,22)*three &
&             -r4(8,31)*three-r4(8,35)*three+r3(4,32)*three+r3(4,40)*three
      rxyz(9)=+r7(19,3)+r7(19,4)-r6(13,5)-r6(13,6)-r6(13,7)-r6(13,8)+r5(8,7)+r5(8,8) &
&             +r5(8,11)+r5(8,12)
      rxyz(10)=+r6(19,11)-r5(13,15)-r5(13,19)+r4(8,20)+r4(8,28)
      rxyz(11)=+r6(13,11)-r5(8,15)-r5(8,19)+r4(4,20)+r4(4,28)
      rxyz(12)=+r7(25,3)-r6(18,5)-r6(18,7)+r5(12,7)+r5(12,11)+r5(14,21)-r4(9,30)-r4(9,34) &
&             +r3(5,31)+r3(5,39)
      rxyz(13)=+r7(25,4)-r6(18,6)-r6(18,8)+r5(12,8)+r5(12,12)+r5(14,22)-r4(9,31)-r4(9,35) &
&             +r3(5,32)+r3(5,40)
      rxyz(14)=+r7(18,4)-r6(12,6)-r6(12,8)+r5(7,8)+r5(7,12)+r5(9,22)-r4(5,31)-r4(5,35) &
&             +r3(2,32)+r3(2,40)
      rxyz(15)=+r7(18,3)-r6(12,5)-r6(12,7)+r5(7,7)+r5(7,11)+r5(9,21)-r4(5,30)-r4(5,34) &
&             +r3(2,31)+r3(2,39)
      rxyz(16)=+r7(13,4)-r6(8,6)-r6(8,8)+r5(4,8)+r5(4,12)+r5(13,22)-r4(8,31)-r4(8,35) &
&             +r3(4,32)+r3(4,40)
      rxyz(17)=+r7(13,3)-r6(8,5)-r6(8,7)+r5(4,7)+r5(4,11)+r5(13,21)-r4(8,30)-r4(8,34) &
&             +r3(4,31)+r3(4,39)
      rxyz(18)=+r7(26,3)-r6(19,5)-r6(19,7)+r5(13,7)+r5(13,11)+r5(13,21)-r4(8,30)-r4(8,34) &
&             +r3(4,31)+r3(4,39)
      rxyz(19)=+r7(26,4)-r6(19,6)-r6(19,8)+r5(13,8)+r5(13,12)+r5(13,22)-r4(8,31)-r4(8,35) &
&             +r3(4,32)+r3(4,40)
      rxyz(20)=+r6(14,11)*four-r5(9,15)*four-r5(9,19)*four+r4(5,20)*four+r4(5,28)*four
      eri(1,1,6,5)=r400+(+r7(9,3)*two+r7(9,4)*two-r6(5,5)*two-r6(5,6)*two-r6(5,7)*two &
&                  -r6(5,8)*two+r5(2,7)*two+r5(2,8)*two+r5(2,11)*two+r5(2,12)*two &
&                  +r5(9,21)*six+r5(9,22)*six-r4(5,30)*six-r4(5,31)*six-r4(5,34)*six &
&                  -r4(5,35)*six+r3(2,31)*six+r3(2,32)*six+r3(2,39)*six+r3(2,40)*six)*qx+( &
&                  +r6(9,10)+r6(9,11)*four+r6(9,12)-r5(5,14)-r5(5,15)*four-r5(5,16)-r5(5,18) &
&                  -r5(5,19)*four-r5(5,20)+r4(2,19)+r4(2,20)*four+r4(2,21)+r4(2,27) &
&                  +r4(2,28)*four+r4(2,29)+r4(9,39)+r4(9,40)*four+r4(9,41)-r3(5,44) &
&                  -r3(5,45)*four-r3(5,46)-r3(5,49)-r3(5,50)*four-r3(5,51)+r2(4,43) &
&                  +r2(4,44)*four+r2(4,45)+r2(4,53)+r2(4,54)*four+r2(4,55))*xx+(+r5(9,23)*two &
&                  +r5(9,24)*two-r4(5,32)*two-r4(5,33)*two-r4(5,36)*two-r4(5,37)*two &
&                  +r3(2,33)*two+r3(2,34)*two+r3(2,41)*two+r3(2,42)*two)*xxx+rxyz(1)*xxxx
      eri(2,1,6,5)=r220+(+r7(18,4)*two-r6(12,6)*two-r6(12,8)*two+r5(7,8)*two+r5(7,12)*two &
&                  +r5(9,22)*two-r4(5,31)*two-r4(5,35)*two+r3(2,32)*two+r3(2,40)*two)*qx &
&                  +rxyz(5)*xx
      eri(3,1,6,5)=r202+(+r7(20,4)*two-r6(14,6)*two-r6(14,8)*two+r5(9,8)*two+r5(9,12)*two &
&                  +r5(9,22)*two-r4(5,31)*two-r4(5,35)*two+r3(2,32)*two+r3(2,40)*two)*qx+( &
&                  +r7(14,3)*two-r6(9,5)*two-r6(9,7)*two+r5(5,7)*two+r5(5,11)*two &
&                  +r5(14,21)*two-r4(9,30)*two-r4(9,34)*two+r3(5,31)*two+r3(5,39)*two)*qz+( &
&                  +r6(20,12)-r5(14,16)-r5(14,20)+r4(9,21)+r4(9,29)+r4(9,41)-r3(5,46) &
&                  -r3(5,51)+r2(4,45)+r2(4,55))*xx+rxyz(20)*xz+(+r6(9,10)-r5(5,14)-r5(5,18) &
&                  +r4(2,19)+r4(2,27)+r4(9,39)-r3(5,44)-r3(5,49)+r2(4,43)+r2(4,53))*zz+( &
&                  +r5(14,24)*two-r4(9,33)*two-r4(9,37)*two+r3(5,34)*two+r3(5,42)*two)*xxz+( &
&                  +r5(9,23)*two-r4(5,32)*two-r4(5,36)*two+r3(2,33)*two+r3(2,41)*two)*xzz &
&                  +rxyz(1)*xxzz
      eri(4,1,6,5)=r310+(+r7(13,3)+r7(13,4)*two-r6(8,5)-r6(8,6)*two-r6(8,7)-r6(8,8)*two &
&                  +r5(4,7)+r5(4,8)*two+r5(4,11)+r5(4,12)*two+r5(13,21)+r5(13,22)*two &
&                  -r4(8,30)-r4(8,31)*two-r4(8,34)-r4(8,35)*two+r3(4,31)+r3(4,32)*two &
&                  +r3(4,39)+r3(4,40)*two)*qx+(+r6(13,11)*two+r6(13,12)-r5(8,15)*two-r5(8,16) &
&                  -r5(8,19)*two-r5(8,20)+r4(4,20)*two+r4(4,21)+r4(4,28)*two+r4(4,29))*xx &
&                  +rxyz(2)*xxx
      eri(5,1,6,5)=r301+(+r7(14,3)+r7(14,4)*two-r6(9,5)-r6(9,6)*two-r6(9,7)-r6(9,8)*two &
&                  +r5(5,7)+r5(5,8)*two+r5(5,11)+r5(5,12)*two+r5(14,21)+r5(14,22)*two &
&                  -r4(9,30)-r4(9,31)*two-r4(9,34)-r4(9,35)*two+r3(5,31)+r3(5,32)*two &
&                  +r3(5,39)+r3(5,40)*two)*qx+(+r7(9,3)-r6(5,5)-r6(5,7)+r5(2,7)+r5(2,11) &
&                  +r5(9,21)*three-r4(5,30)*three-r4(5,34)*three+r3(2,31)*three &
&                  +r3(2,39)*three)*qz+(+r6(14,11)*two+r6(14,12)-r5(9,15)*two-r5(9,16) &
&                  -r5(9,19)*two-r5(9,20)+r4(5,20)*two+r4(5,21)+r4(5,28)*two+r4(5,29))*xx+( &
&                  +r6(9,10)+r6(9,11)*two-r5(5,14)-r5(5,15)*two-r5(5,18)-r5(5,19)*two &
&                  +r4(2,19)+r4(2,20)*two+r4(2,27)+r4(2,28)*two+r4(9,39)+r4(9,40)*two &
&                  -r3(5,44)-r3(5,45)*two-r3(5,49)-r3(5,50)*two+r2(4,43)+r2(4,44)*two &
&                  +r2(4,53)+r2(4,54)*two)*xz+(+r5(14,24)-r4(9,33)-r4(9,37)+r3(5,34)+r3(5,42) &
&                  )*xxx+(+r5(9,23)*two+r5(9,24)-r4(5,32)*two-r4(5,33)-r4(5,36)*two-r4(5,37) &
&                  +r3(2,33)*two+r3(2,34)+r3(2,41)*two+r3(2,42))*xxz+rxyz(1)*xxxz
      eri(6,1,6,5)=r211+(+r7(19,4)*two-r6(13,6)*two-r6(13,8)*two+r5(8,8)*two+r5(8,12)*two &
&                  )*qx+rxyz(17)*qz+(+r6(19,12)-r5(13,16)-r5(13,20)+r4(8,21)+r4(8,29))*xx+( &
&                  +r6(13,11)*two-r5(8,15)*two-r5(8,19)*two+r4(4,20)*two+r4(4,28)*two)*xz &
&                  +rxyz(2)*xxz
      eri(1,2,6,5)=r220+(+r7(18,3)*two-r6(12,5)*two-r6(12,7)*two+r5(7,7)*two+r5(7,11)*two &
&                  +r5(9,21)*two-r4(5,30)*two-r4(5,34)*two+r3(2,31)*two+r3(2,39)*two)*qx &
&                  +rxyz(4)*xx
      eri(2,2,6,5)=r040
      eri(3,2,6,5)=r022+(+r7(25,3)*two-r6(18,5)*two-r6(18,7)*two+r5(12,7)*two &
&                  +r5(12,11)*two+r5(14,21)*two-r4(9,30)*two-r4(9,34)*two+r3(5,31)*two &
&                  +r3(5,39)*two)*qz+rxyz(4)*zz
      eri(4,2,6,5)=r130+rxyz(7)*qx
      eri(5,2,6,5)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,6,5)=r031+rxyz(7)*qz
      eri(1,3,6,5)=r202+(+r7(20,3)*two-r6(14,5)*two-r6(14,7)*two+r5(9,7)*two+r5(9,11)*two &
&                  +r5(9,21)*two-r4(5,30)*two-r4(5,34)*two+r3(2,31)*two+r3(2,39)*two)*qx+( &
&                  +r7(14,4)*two-r6(9,6)*two-r6(9,8)*two+r5(5,8)*two+r5(5,12)*two &
&                  +r5(14,22)*two-r4(9,31)*two-r4(9,35)*two+r3(5,32)*two+r3(5,40)*two)*qz+( &
&                  +r6(20,10)-r5(14,14)-r5(14,18)+r4(9,19)+r4(9,27)+r4(9,39)-r3(5,44) &
&                  -r3(5,49)+r2(4,43)+r2(4,53))*xx+rxyz(20)*xz+(+r6(9,12)-r5(5,16)-r5(5,20) &
&                  +r4(2,21)+r4(2,29)+r4(9,41)-r3(5,46)-r3(5,51)+r2(4,45)+r2(4,55))*zz+( &
&                  +r5(14,23)*two-r4(9,32)*two-r4(9,36)*two+r3(5,33)*two+r3(5,41)*two)*xxz+( &
&                  +r5(9,24)*two-r4(5,33)*two-r4(5,37)*two+r3(2,34)*two+r3(2,42)*two)*xzz &
&                  +rxyz(1)*xxzz
      eri(2,3,6,5)=r022+(+r7(25,4)*two-r6(18,6)*two-r6(18,8)*two+r5(12,8)*two &
&                  +r5(12,12)*two+r5(14,22)*two-r4(9,31)*two-r4(9,35)*two+r3(5,32)*two &
&                  +r3(5,40)*two)*qz+rxyz(5)*zz
      eri(3,3,6,5)=r004+(+r7(27,3)*two+r7(27,4)*two-r6(20,5)*two-r6(20,6)*two &
&                  -r6(20,7)*two-r6(20,8)*two+r5(14,7)*two+r5(14,8)*two+r5(14,11)*two &
&                  +r5(14,12)*two+r5(14,21)*six+r5(14,22)*six-r4(9,30)*six-r4(9,31)*six &
&                  -r4(9,34)*six-r4(9,35)*six+r3(5,31)*six+r3(5,32)*six+r3(5,39)*six &
&                  +r3(5,40)*six)*qz+(+r6(20,10)+r6(20,11)*four+r6(20,12)-r5(14,14) &
&                  -r5(14,15)*four-r5(14,16)-r5(14,18)-r5(14,19)*four-r5(14,20)+r4(9,19) &
&                  +r4(9,20)*four+r4(9,21)+r4(9,27)+r4(9,28)*four+r4(9,29)+r4(9,39) &
&                  +r4(9,40)*four+r4(9,41)-r3(5,44)-r3(5,45)*four-r3(5,46)-r3(5,49) &
&                  -r3(5,50)*four-r3(5,51)+r2(4,43)+r2(4,44)*four+r2(4,45)+r2(4,53) &
&                  +r2(4,54)*four+r2(4,55))*zz+(+r5(14,23)*two+r5(14,24)*two-r4(9,32)*two &
&                  -r4(9,33)*two-r4(9,36)*two-r4(9,37)*two+r3(5,33)*two+r3(5,34)*two &
&                  +r3(5,41)*two+r3(5,42)*two)*zzz+rxyz(1)*zzzz
      eri(4,3,6,5)=r112+rxyz(18)*qx+(+r7(19,4)*two-r6(13,6)*two-r6(13,8)*two+r5(8,8)*two &
&                  +r5(8,12)*two)*qz+(+r6(19,11)*two-r5(13,15)*two-r5(13,19)*two+r4(8,20)*two &
&                  +r4(8,28)*two)*xz+(+r6(13,12)-r5(8,16)-r5(8,20)+r4(4,21)+r4(4,29))*zz &
&                  +rxyz(2)*xzz
      eri(5,3,6,5)=r103+(+r7(27,3)-r6(20,5)-r6(20,7)+r5(14,7)+r5(14,11)+r5(14,21)*three &
&                  -r4(9,30)*three-r4(9,34)*three+r3(5,31)*three+r3(5,39)*three)*qx+( &
&                  +r7(20,3)+r7(20,4)*two-r6(14,5)-r6(14,6)*two-r6(14,7)-r6(14,8)*two+r5(9,7) &
&                  +r5(9,8)*two+r5(9,11)+r5(9,12)*two+r5(9,21)+r5(9,22)*two-r4(5,30) &
&                  -r4(5,31)*two-r4(5,34)-r4(5,35)*two+r3(2,31)+r3(2,32)*two+r3(2,39) &
&                  +r3(2,40)*two)*qz+(+r6(20,10)+r6(20,11)*two-r5(14,14)-r5(14,15)*two &
&                  -r5(14,18)-r5(14,19)*two+r4(9,19)+r4(9,20)*two+r4(9,27)+r4(9,28)*two &
&                  +r4(9,39)+r4(9,40)*two-r3(5,44)-r3(5,45)*two-r3(5,49)-r3(5,50)*two &
&                  +r2(4,43)+r2(4,44)*two+r2(4,53)+r2(4,54)*two)*xz+(+r6(14,11)*two+r6(14,12) &
&                  -r5(9,15)*two-r5(9,16)-r5(9,19)*two-r5(9,20)+r4(5,20)*two+r4(5,21) &
&                  +r4(5,28)*two+r4(5,29))*zz+(+r5(14,23)*two+r5(14,24)-r4(9,32)*two-r4(9,33) &
&                  -r4(9,36)*two-r4(9,37)+r3(5,33)*two+r3(5,34)+r3(5,41)*two+r3(5,42))*xzz+( &
&                  +r5(9,24)-r4(5,33)-r4(5,37)+r3(2,34)+r3(2,42))*zzz+rxyz(1)*xzzz
      eri(6,3,6,5)=r013+(+r7(26,3)+r7(26,4)*two-r6(19,5)-r6(19,6)*two-r6(19,7) &
&                  -r6(19,8)*two+r5(13,7)+r5(13,8)*two+r5(13,11)+r5(13,12)*two+r5(13,21) &
&                  +r5(13,22)*two-r4(8,30)-r4(8,31)*two-r4(8,34)-r4(8,35)*two+r3(4,31) &
&                  +r3(4,32)*two+r3(4,39)+r3(4,40)*two)*qz+(+r6(19,11)*two+r6(19,12) &
&                  -r5(13,15)*two-r5(13,16)-r5(13,19)*two-r5(13,20)+r4(8,20)*two+r4(8,21) &
&                  +r4(8,28)*two+r4(8,29))*zz+rxyz(2)*zzz
      eri(1,4,6,5)=r310+(+r7(13,3)*two+r7(13,4)-r6(8,5)*two-r6(8,6)-r6(8,7)*two-r6(8,8) &
&                  +r5(4,7)*two+r5(4,8)+r5(4,11)*two+r5(4,12)+r5(13,21)*two+r5(13,22) &
&                  -r4(8,30)*two-r4(8,31)-r4(8,34)*two-r4(8,35)+r3(4,31)*two+r3(4,32) &
&                  +r3(4,39)*two+r3(4,40))*qx+(+r6(13,10)+r6(13,11)*two-r5(8,14)-r5(8,15)*two &
&                  -r5(8,18)-r5(8,19)*two+r4(4,19)+r4(4,20)*two+r4(4,27)+r4(4,28)*two)*xx &
&                  +rxyz(3)*xxx
      eri(2,4,6,5)=r130+rxyz(8)*qx
      eri(3,4,6,5)=r112+rxyz(19)*qx+(+r7(19,3)*two-r6(13,5)*two-r6(13,7)*two+r5(8,7)*two &
&                  +r5(8,11)*two)*qz+(+r6(19,11)*two-r5(13,15)*two-r5(13,19)*two+r4(8,20)*two &
&                  +r4(8,28)*two)*xz+(+r6(13,10)-r5(8,14)-r5(8,18)+r4(4,19)+r4(4,27))*zz &
&                  +rxyz(3)*xzz
      eri(4,4,6,5)=r220+(+r7(18,3)+r7(18,4)-r6(12,5)-r6(12,6)-r6(12,7)-r6(12,8)+r5(7,7) &
&                  +r5(7,8)+r5(7,11)+r5(7,12)+r5(9,21)+r5(9,22)-r4(5,30)-r4(5,31)-r4(5,34) &
&                  -r4(5,35)+r3(2,31)+r3(2,32)+r3(2,39)+r3(2,40))*qx+rxyz(6)*xx
      eri(5,4,6,5)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(13,10)+r6(13,11)-r5(8,14) &
&                  -r5(8,15)-r5(8,18)-r5(8,19)+r4(4,19)+r4(4,20)+r4(4,27)+r4(4,28))*xz &
&                  +rxyz(3)*xxz
      eri(6,4,6,5)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,6,5)=r301+(+r7(14,3)*two+r7(14,4)-r6(9,5)*two-r6(9,6)-r6(9,7)*two-r6(9,8) &
&                  +r5(5,7)*two+r5(5,8)+r5(5,11)*two+r5(5,12)+r5(14,21)*two+r5(14,22) &
&                  -r4(9,30)*two-r4(9,31)-r4(9,34)*two-r4(9,35)+r3(5,31)*two+r3(5,32) &
&                  +r3(5,39)*two+r3(5,40))*qx+(+r7(9,4)-r6(5,6)-r6(5,8)+r5(2,8)+r5(2,12) &
&                  +r5(9,22)*three-r4(5,31)*three-r4(5,35)*three+r3(2,32)*three &
&                  +r3(2,40)*three)*qz+(+r6(14,10)+r6(14,11)*two-r5(9,14)-r5(9,15)*two &
&                  -r5(9,18)-r5(9,19)*two+r4(5,19)+r4(5,20)*two+r4(5,27)+r4(5,28)*two)*xx+( &
&                  +r6(9,11)*two+r6(9,12)-r5(5,15)*two-r5(5,16)-r5(5,19)*two-r5(5,20) &
&                  +r4(2,20)*two+r4(2,21)+r4(2,28)*two+r4(2,29)+r4(9,40)*two+r4(9,41) &
&                  -r3(5,45)*two-r3(5,46)-r3(5,50)*two-r3(5,51)+r2(4,44)*two+r2(4,45) &
&                  +r2(4,54)*two+r2(4,55))*xz+(+r5(14,23)-r4(9,32)-r4(9,36)+r3(5,33)+r3(5,41) &
&                  )*xxx+(+r5(9,23)+r5(9,24)*two-r4(5,32)-r4(5,33)*two-r4(5,36)-r4(5,37)*two &
&                  +r3(2,33)+r3(2,34)*two+r3(2,41)+r3(2,42)*two)*xxz+rxyz(1)*xxxz
      eri(2,5,6,5)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,6,5)=r103+(+r7(27,4)-r6(20,6)-r6(20,8)+r5(14,8)+r5(14,12)+r5(14,22)*three &
&                  -r4(9,31)*three-r4(9,35)*three+r3(5,32)*three+r3(5,40)*three)*qx+( &
&                  +r7(20,3)*two+r7(20,4)-r6(14,5)*two-r6(14,6)-r6(14,7)*two-r6(14,8) &
&                  +r5(9,7)*two+r5(9,8)+r5(9,11)*two+r5(9,12)+r5(9,21)*two+r5(9,22) &
&                  -r4(5,30)*two-r4(5,31)-r4(5,34)*two-r4(5,35)+r3(2,31)*two+r3(2,32) &
&                  +r3(2,39)*two+r3(2,40))*qz+(+r6(20,11)*two+r6(20,12)-r5(14,15)*two &
&                  -r5(14,16)-r5(14,19)*two-r5(14,20)+r4(9,20)*two+r4(9,21)+r4(9,28)*two &
&                  +r4(9,29)+r4(9,40)*two+r4(9,41)-r3(5,45)*two-r3(5,46)-r3(5,50)*two &
&                  -r3(5,51)+r2(4,44)*two+r2(4,45)+r2(4,54)*two+r2(4,55))*xz+(+r6(14,10) &
&                  +r6(14,11)*two-r5(9,14)-r5(9,15)*two-r5(9,18)-r5(9,19)*two+r4(5,19) &
&                  +r4(5,20)*two+r4(5,27)+r4(5,28)*two)*zz+(+r5(14,23)+r5(14,24)*two-r4(9,32) &
&                  -r4(9,33)*two-r4(9,36)-r4(9,37)*two+r3(5,33)+r3(5,34)*two+r3(5,41) &
&                  +r3(5,42)*two)*xzz+(+r5(9,23)-r4(5,32)-r4(5,36)+r3(2,33)+r3(2,41))*zzz &
&                  +rxyz(1)*xzzz
      eri(4,5,6,5)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(13,11)+r6(13,12)-r5(8,15) &
&                  -r5(8,16)-r5(8,19)-r5(8,20)+r4(4,20)+r4(4,21)+r4(4,28)+r4(4,29))*xz &
&                  +rxyz(2)*xxz
      eri(5,5,6,5)=r202+(+r7(20,3)+r7(20,4)-r6(14,5)-r6(14,6)-r6(14,7)-r6(14,8)+r5(9,7) &
&                  +r5(9,8)+r5(9,11)+r5(9,12)+r5(9,21)+r5(9,22)-r4(5,30)-r4(5,31)-r4(5,34) &
&                  -r4(5,35)+r3(2,31)+r3(2,32)+r3(2,39)+r3(2,40))*qx+(+r7(14,3)+r7(14,4) &
&                  -r6(9,5)-r6(9,6)-r6(9,7)-r6(9,8)+r5(5,7)+r5(5,8)+r5(5,11)+r5(5,12) &
&                  +r5(14,21)+r5(14,22)-r4(9,30)-r4(9,31)-r4(9,34)-r4(9,35)+r3(5,31)+r3(5,32) &
&                  +r3(5,39)+r3(5,40))*qz+(+r6(20,11)-r5(14,15)-r5(14,19)+r4(9,20)+r4(9,28) &
&                  +r4(9,40)-r3(5,45)-r3(5,50)+r2(4,44)+r2(4,54))*xx+(+r6(14,10) &
&                  +r6(14,11)*two+r6(14,12)-r5(9,14)-r5(9,15)*two-r5(9,16)-r5(9,18) &
&                  -r5(9,19)*two-r5(9,20)+r4(5,19)+r4(5,20)*two+r4(5,21)+r4(5,27) &
&                  +r4(5,28)*two+r4(5,29))*xz+(+r6(9,11)-r5(5,15)-r5(5,19)+r4(2,20)+r4(2,28) &
&                  +r4(9,40)-r3(5,45)-r3(5,50)+r2(4,44)+r2(4,54))*zz+(+r5(14,23)+r5(14,24) &
&                  -r4(9,32)-r4(9,33)-r4(9,36)-r4(9,37)+r3(5,33)+r3(5,34)+r3(5,41)+r3(5,42)) &
&                  *xxz+(+r5(9,23)+r5(9,24)-r4(5,32)-r4(5,33)-r4(5,36)-r4(5,37)+r3(2,33) &
&                  +r3(2,34)+r3(2,41)+r3(2,42))*xzz+rxyz(1)*xxzz
      eri(6,5,6,5)=r112+rxyz(19)*qx+(+r7(19,3)+r7(19,4)-r6(13,5)-r6(13,6)-r6(13,7) &
&                  -r6(13,8)+r5(8,7)+r5(8,8)+r5(8,11)+r5(8,12))*qz+(+r6(19,11)+r6(19,12) &
&                  -r5(13,15)-r5(13,16)-r5(13,19)-r5(13,20)+r4(8,20)+r4(8,21)+r4(8,28) &
&                  +r4(8,29))*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,6,5)=r211+(+r7(19,3)*two-r6(13,5)*two-r6(13,7)*two+r5(8,7)*two+r5(8,11)*two &
&                  )*qx+rxyz(16)*qz+(+r6(19,10)-r5(13,14)-r5(13,18)+r4(8,19)+r4(8,27))*xx+( &
&                  +r6(13,11)*two-r5(8,15)*two-r5(8,19)*two+r4(4,20)*two+r4(4,28)*two)*xz &
&                  +rxyz(3)*xxz
      eri(2,6,6,5)=r031+rxyz(8)*qz
      eri(3,6,6,5)=r013+(+r7(26,3)*two+r7(26,4)-r6(19,5)*two-r6(19,6)-r6(19,7)*two &
&                  -r6(19,8)+r5(13,7)*two+r5(13,8)+r5(13,11)*two+r5(13,12)+r5(13,21)*two &
&                  +r5(13,22)-r4(8,30)*two-r4(8,31)-r4(8,34)*two-r4(8,35)+r3(4,31)*two &
&                  +r3(4,32)+r3(4,39)*two+r3(4,40))*qz+(+r6(19,10)+r6(19,11)*two-r5(13,14) &
&                  -r5(13,15)*two-r5(13,18)-r5(13,19)*two+r4(8,19)+r4(8,20)*two+r4(8,27) &
&                  +r4(8,28)*two)*zz+rxyz(3)*zzz
      eri(4,6,6,5)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,6,5)=r112+rxyz(18)*qx+(+r7(19,3)+r7(19,4)-r6(13,5)-r6(13,6)-r6(13,7) &
&                  -r6(13,8)+r5(8,7)+r5(8,8)+r5(8,11)+r5(8,12))*qz+(+r6(19,10)+r6(19,11) &
&                  -r5(13,14)-r5(13,15)-r5(13,18)-r5(13,19)+r4(8,19)+r4(8,20)+r4(8,27) &
&                  +r4(8,28))*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,6,5)=r022+(+r7(25,3)+r7(25,4)-r6(18,5)-r6(18,6)-r6(18,7)-r6(18,8)+r5(12,7) &
&                  +r5(12,8)+r5(12,11)+r5(12,12)+r5(14,21)+r5(14,22)-r4(9,30)-r4(9,31) &
&                  -r4(9,34)-r4(9,35)+r3(5,31)+r3(5,32)+r3(5,39)+r3(5,40))*qz+rxyz(6)*zz
      return
end

!-------------------------------------------------------------
  subroutine int2dddd5(eri,r0,r1,r2,r3,r4,r5,r6,r7,r8,qx,qz)
!-------------------------------------------------------------
!
      implicit none
      integer :: i,j
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00, four=4.0D+00, five=5.0D+00
      real(8),parameter :: six=6.0D+00, seven=7.0D+00, eight=8.0D+00, nine=9.0D+00, ten=1.0D+01
      real(8),parameter :: p11=1.1D+01, p12=1.2D+01, p13=1.3D+01, p15=1.5D+01, p16=1.6D+01
      real(8),parameter :: p18=1.8D+01, p20=2.0D+01, p21=2.1D+01, p24=2.4D+01, p28=2.8D+01
      real(8),parameter :: p30=3.0D+01, p36=3.6D+01, p45=4.5D+01, p105=1.05D+2, p210=2.1D+02
      real(8),parameter :: p420=4.2D+02
      real(8),intent(in) :: r0(25), r1(3,40), r2(6,56), r3(10,52), r4(15,42), r5(21,24)
      real(8),intent(in) :: r6(28,12), r7(36,4), r8(45), qx, qz
      real(8),intent(inout) :: eri(6,6,6,6)
      real(8) :: r400, r310, r301, r220, r211, r202, r130, r121, r112, r103, r040, r031
      real(8) :: r022, r013, r004, rxyz(20), xx, xz, zz, xxx, xxz, xzz, zzz
      real(8) :: xxxx, xxxz, xxzz, xzzz, zzzz
!
      xx= qx*qx
      xz= qx*qz
      zz= qz*qz
      xxx= qx*qx*qx
      xxz= qx*qx*qz
      xzz= qx*qz*qz
      zzz= qz*qz*qz
      xxxx= xx*xx
      xxxz= xx*xz
      xxzz= xx*zz
      xzzz= xz*zz
      zzzz= zz*zz
!
!
      do i= 1,6
        do j= 1,6
          eri(j,i,1,6)= eri(j,i,4,5)
        enddo
      enddo
!
      r400= r8(12)-r7(7,2)+r6(5,4)*three+r6(12,9)*six-r5(2,4)*three-r5(7,17)*six &
&          +r4(5,26)*p18+r4(12,38)*three-r3(2,23)*p18-r3(7,48)*three+r2(6,52)*nine &
&          -r1(2,36)*nine
      r310= r8(17)-r7(11,2)+r6(8,4)*three+r6(17,9)*three-r5(4,4)*three-r5(11,17)*three &
&          +r4(8,26)*nine-r3(4,23)*nine
      r301= r8(18)-r7(12,2)+r6(9,4)*three+r6(18,9)*three-r5(5,4)*three-r5(12,17)*three &
&          +r4(9,26)*nine-r3(5,23)*nine
      r220= r8(23)-r7(16,2)+r6(12,4)*three+r6(12,9)+r6(23,9)-r5(7,4)*three-r5(7,17) &
&          -r5(16,17)+r4(5,26)*three+r4(12,26)*three+r4(12,38)-r3(2,23)*three &
&          -r3(7,23)*three-r3(7,48)+r2(6,52)*three-r1(2,36)*three
      r211= r8(24)-r7(17,2)+r6(13,4)*three+r6(24,9)-r5(8,4)*three-r5(17,17)+r4(13,26)*three &
&          -r3(8,23)*three
      r202= r8(25)-r7(18,2)+r6(14,4)*three+r6(12,9)+r6(25,9)-r5(9,4)*three-r5(7,17) &
&          -r5(18,17)+r4(5,26)*three+r4(14,26)*three+r4(12,38)-r3(2,23)*three &
&          -r3(9,23)*three-r3(7,48)+r2(6,52)*three-r1(2,36)*three
      r130= r8(30)-r7(22,2)+r6(17,4)*three+r6(17,9)*three-r5(11,4)*three-r5(11,17)*three &
&          +r4(8,26)*nine-r3(4,23)*nine
      r121= r8(31)-r7(23,2)+r6(18,4)*three+r6(18,9)-r5(12,4)*three-r5(12,17)+r4(9,26)*three &
&          -r3(5,23)*three
      r112= r8(32)-r7(24,2)+r6(19,4)*three+r6(17,9)-r5(13,4)*three-r5(11,17)+r4(8,26)*three &
&          -r3(4,23)*three
      r103= r8(33)-r7(25,2)+r6(20,4)*three+r6(18,9)*three-r5(14,4)*three-r5(12,17)*three &
&          +r4(9,26)*nine-r3(5,23)*nine
      r040= r8(38)-r7(29,2)+r6(23,4)*three+r6(23,9)*six-r5(16,4)*three-r5(16,17)*six &
&          +r4(12,26)*p18+r4(12,38)*three-r3(7,23)*p18-r3(7,48)*three+r2(6,52)*nine &
&          -r1(2,36)*nine
      r031= r8(39)-r7(30,2)+r6(24,4)*three+r6(24,9)*three-r5(17,4)*three-r5(17,17)*three &
&          +r4(13,26)*nine-r3(8,23)*nine
      r022= r8(40)-r7(31,2)+r6(25,4)*three+r6(23,9)+r6(25,9)-r5(18,4)*three-r5(16,17) &
&          -r5(18,17)+r4(12,26)*three+r4(14,26)*three+r4(12,38)-r3(7,23)*three &
&          -r3(9,23)*three-r3(7,48)+r2(6,52)*three-r1(2,36)*three
      r013= r8(41)-r7(32,2)+r6(26,4)*three+r6(24,9)*three-r5(19,4)*three-r5(17,17)*three &
&          +r4(13,26)*nine-r3(8,23)*nine
      r004= r8(42)-r7(33,2)+r6(27,4)*three+r6(25,9)*six-r5(20,4)*three-r5(18,17)*six &
&          +r4(14,26)*p18+r4(12,38)*three-r3(9,23)*p18-r3(7,48)*three+r2(6,52)*nine &
&          -r1(2,36)*nine
      rxyz(1)=+r4(12,42)-r3(7,52)+r2(6,56)*three-r1(2,40)*three
      rxyz(2)=+r5(17,24)-r4(11,37)+r3(8,42)*three-r2(2,36)*three
      rxyz(3)=+r5(17,23)-r4(11,36)+r3(8,41)*three-r2(2,35)*three
      rxyz(4)=+r6(23,10)-r5(16,18)+r4(12,27)*three+r4(12,39)-r3(7,24)*three-r3(7,49) &
&             +r2(6,53)*three-r1(2,37)*three
      rxyz(5)=+r6(23,12)-r5(16,20)+r4(12,29)*three+r4(12,41)-r3(7,26)*three-r3(7,51) &
&             +r2(6,55)*three-r1(2,39)*three
      rxyz(6)=+r6(23,11)-r5(16,19)+r4(12,28)*three+r4(12,40)-r3(7,25)*three-r3(7,50) &
&             +r2(6,54)*three-r1(2,38)*three
      rxyz(7)=+r7(30,3)-r6(22,7)+r5(17,11)*three+r5(17,21)*three-r4(11,12)*three &
&             -r4(11,34)*three+r3(8,39)*nine-r2(2,33)*nine
      rxyz(8)=+r7(30,4)-r6(22,8)+r5(17,12)*three+r5(17,22)*three-r4(11,13)*three &
&             -r4(11,35)*three+r3(8,40)*nine-r2(2,34)*nine
      rxyz(9)=+r7(24,3)+r7(24,4)-r6(17,7)-r6(17,8)+r5(13,11)*three+r5(13,12)*three &
&             -r4(8,12)*three-r4(8,13)*three
      rxyz(10)=+r6(24,11)-r5(17,19)+r4(13,28)*three-r3(8,25)*three
      rxyz(11)=+r6(17,11)-r5(11,19)+r4(8,28)*three-r3(4,25)*three
      rxyz(12)=+r7(31,3)-r6(23,7)+r5(18,11)*three+r5(18,21)-r4(12,12)*three-r4(12,34) &
&             +r3(9,39)*three-r2(6,33)*three
      rxyz(13)=+r7(31,4)-r6(23,8)+r5(18,12)*three+r5(18,22)-r4(12,13)*three-r4(12,35) &
&             +r3(9,40)*three-r2(6,34)*three
      rxyz(14)=+r7(23,4)-r6(16,8)+r5(12,12)*three+r5(12,22)-r4(7,13)*three-r4(7,35) &
&             +r3(5,40)*three-r2(4,34)*three
      rxyz(15)=+r7(23,3)-r6(16,7)+r5(12,11)*three+r5(12,21)-r4(7,12)*three-r4(7,34) &
&             +r3(5,39)*three-r2(4,33)*three
      rxyz(16)=+r7(17,4)-r6(11,8)+r5(8,12)*three+r5(17,22)-r4(4,13)*three-r4(11,35) &
&             +r3(8,40)*three-r2(2,34)*three
      rxyz(17)=+r7(17,3)-r6(11,7)+r5(8,11)*three+r5(17,21)-r4(4,12)*three-r4(11,34) &
&             +r3(8,39)*three-r2(2,33)*three
      rxyz(18)=+r7(32,3)-r6(24,7)+r5(19,11)*three+r5(17,21)-r4(13,12)*three-r4(11,34) &
&             +r3(8,39)*three-r2(2,33)*three
      rxyz(19)=+r7(32,4)-r6(24,8)+r5(19,12)*three+r5(17,22)-r4(13,13)*three-r4(11,35) &
&             +r3(8,40)*three-r2(2,34)*three
      rxyz(20)=+r6(18,11)*four-r5(12,19)*four+r4(9,28)*p12-r3(5,25)*p12
      eri(1,1,2,6)=r400+(+r7(12,3)*two+r7(12,4)*two-r6(7,7)*two-r6(7,8)*two+r5(5,11)*six &
&                  +r5(5,12)*six+r5(12,21)*six+r5(12,22)*six-r4(2,12)*six-r4(2,13)*six &
&                  -r4(7,34)*six-r4(7,35)*six+r3(5,39)*p18+r3(5,40)*p18-r2(4,33)*p18 &
&                  -r2(4,34)*p18)*qx+(+r6(12,10)+r6(12,11)*four+r6(12,12)-r5(7,18) &
&                  -r5(7,19)*four-r5(7,20)+r4(5,27)*three+r4(5,28)*p12+r4(5,29)*three &
&                  +r4(12,39)+r4(12,40)*four+r4(12,41)-r3(2,24)*three-r3(2,25)*p12 &
&                  -r3(2,26)*three-r3(7,49)-r3(7,50)*four-r3(7,51)+r2(6,53)*three &
&                  +r2(6,54)*p12+r2(6,55)*three-r1(2,37)*three-r1(2,38)*p12-r1(2,39)*three) &
&                  *xx+(+r5(12,23)*two+r5(12,24)*two-r4(7,36)*two-r4(7,37)*two+r3(5,41)*six &
&                  +r3(5,42)*six-r2(4,35)*six-r2(4,36)*six)*xxx+rxyz(1)*xxxx
      eri(2,1,2,6)=r220+(+r7(23,4)*two-r6(16,8)*two+r5(12,12)*six+r5(12,22)*two &
&                  -r4(7,13)*six-r4(7,35)*two+r3(5,40)*six-r2(4,34)*six)*qx+rxyz(5)*xx
      eri(3,1,2,6)=r202+(+r7(25,4)*two-r6(18,8)*two+r5(14,12)*six+r5(12,22)*two &
&                  -r4(9,13)*six-r4(7,35)*two+r3(5,40)*six-r2(4,34)*six)*qx+(+r7(18,3)*two &
&                  -r6(12,7)*two+r5(9,11)*six+r5(18,21)*two-r4(5,12)*six-r4(12,34)*two &
&                  +r3(9,39)*six-r2(6,33)*six)*qz+(+r6(25,12)-r5(18,20)+r4(14,29)*three &
&                  +r4(12,41)-r3(9,26)*three-r3(7,51)+r2(6,55)*three-r1(2,39)*three)*xx &
&                  +rxyz(20)*xz+(+r6(12,10)-r5(7,18)+r4(5,27)*three+r4(12,39)-r3(2,24)*three &
&                  -r3(7,49)+r2(6,53)*three-r1(2,37)*three)*zz+(+r5(18,24)*two-r4(12,37)*two &
&                  +r3(9,42)*six-r2(6,36)*six)*xxz+(+r5(12,23)*two-r4(7,36)*two+r3(5,41)*six &
&                  -r2(4,35)*six)*xzz+rxyz(1)*xxzz
      eri(4,1,2,6)=r310+(+r7(17,3)+r7(17,4)*two-r6(11,7)-r6(11,8)*two+r5(8,11)*three &
&                  +r5(8,12)*six+r5(17,21)+r5(17,22)*two-r4(4,12)*three-r4(4,13)*six &
&                  -r4(11,34)-r4(11,35)*two+r3(8,39)*three+r3(8,40)*six-r2(2,33)*three &
&                  -r2(2,34)*six)*qx+(+r6(17,11)*two+r6(17,12)-r5(11,19)*two-r5(11,20) &
&                  +r4(8,28)*six+r4(8,29)*three-r3(4,25)*six-r3(4,26)*three)*xx+rxyz(2)*xxx
      eri(5,1,2,6)=r301+(+r7(18,3)+r7(18,4)*two-r6(12,7)-r6(12,8)*two+r5(9,11)*three &
&                  +r5(9,12)*six+r5(18,21)+r5(18,22)*two-r4(5,12)*three-r4(5,13)*six &
&                  -r4(12,34)-r4(12,35)*two+r3(9,39)*three+r3(9,40)*six-r2(6,33)*three &
&                  -r2(6,34)*six)*qx+(+r7(12,3)-r6(7,7)+r5(5,11)*three+r5(12,21)*three &
&                  -r4(2,12)*three-r4(7,34)*three+r3(5,39)*nine-r2(4,33)*nine)*qz+( &
&                  +r6(18,11)*two+r6(18,12)-r5(12,19)*two-r5(12,20)+r4(9,28)*six &
&                  +r4(9,29)*three-r3(5,25)*six-r3(5,26)*three)*xx+(+r6(12,10)+r6(12,11)*two &
&                  -r5(7,18)-r5(7,19)*two+r4(5,27)*three+r4(5,28)*six+r4(12,39)+r4(12,40)*two &
&                  -r3(2,24)*three-r3(2,25)*six-r3(7,49)-r3(7,50)*two+r2(6,53)*three &
&                  +r2(6,54)*six-r1(2,37)*three-r1(2,38)*six)*xz+(+r5(18,24)-r4(12,37) &
&                  +r3(9,42)*three-r2(6,36)*three)*xxx+(+r5(12,23)*two+r5(12,24)-r4(7,36)*two &
&                  -r4(7,37)+r3(5,41)*six+r3(5,42)*three-r2(4,35)*six-r2(4,36)*three)*xxz &
&                  +rxyz(1)*xxxz
      eri(6,1,2,6)=r211+(+r7(24,4)*two-r6(17,8)*two+r5(13,12)*six-r4(8,13)*six)*qx &
&                  +rxyz(17)*qz+(+r6(24,12)-r5(17,20)+r4(13,29)*three-r3(8,26)*three)*xx+( &
&                  +r6(17,11)*two-r5(11,19)*two+r4(8,28)*six-r3(4,25)*six)*xz+rxyz(2)*xxz
      eri(1,2,2,6)=r220+(+r7(23,3)*two-r6(16,7)*two+r5(12,11)*six+r5(12,21)*two &
&                  -r4(7,12)*six-r4(7,34)*two+r3(5,39)*six-r2(4,33)*six)*qx+rxyz(4)*xx
      eri(2,2,2,6)=r040
      eri(3,2,2,6)=r022+(+r7(31,3)*two-r6(23,7)*two+r5(18,11)*six+r5(18,21)*two &
&                  -r4(12,12)*six-r4(12,34)*two+r3(9,39)*six-r2(6,33)*six)*qz+rxyz(4)*zz
      eri(4,2,2,6)=r130+rxyz(7)*qx
      eri(5,2,2,6)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,2,6)=r031+rxyz(7)*qz
      eri(1,3,2,6)=r202+(+r7(25,3)*two-r6(18,7)*two+r5(14,11)*six+r5(12,21)*two &
&                  -r4(9,12)*six-r4(7,34)*two+r3(5,39)*six-r2(4,33)*six)*qx+(+r7(18,4)*two &
&                  -r6(12,8)*two+r5(9,12)*six+r5(18,22)*two-r4(5,13)*six-r4(12,35)*two &
&                  +r3(9,40)*six-r2(6,34)*six)*qz+(+r6(25,10)-r5(18,18)+r4(14,27)*three &
&                  +r4(12,39)-r3(9,24)*three-r3(7,49)+r2(6,53)*three-r1(2,37)*three)*xx &
&                  +rxyz(20)*xz+(+r6(12,12)-r5(7,20)+r4(5,29)*three+r4(12,41)-r3(2,26)*three &
&                  -r3(7,51)+r2(6,55)*three-r1(2,39)*three)*zz+(+r5(18,23)*two-r4(12,36)*two &
&                  +r3(9,41)*six-r2(6,35)*six)*xxz+(+r5(12,24)*two-r4(7,37)*two+r3(5,42)*six &
&                  -r2(4,36)*six)*xzz+rxyz(1)*xxzz
      eri(2,3,2,6)=r022+(+r7(31,4)*two-r6(23,8)*two+r5(18,12)*six+r5(18,22)*two &
&                  -r4(12,13)*six-r4(12,35)*two+r3(9,40)*six-r2(6,34)*six)*qz+rxyz(5)*zz
      eri(3,3,2,6)=r004+(+r7(33,3)*two+r7(33,4)*two-r6(25,7)*two-r6(25,8)*two &
&                  +r5(20,11)*six+r5(20,12)*six+r5(18,21)*six+r5(18,22)*six-r4(14,12)*six &
&                  -r4(14,13)*six-r4(12,34)*six-r4(12,35)*six+r3(9,39)*p18+r3(9,40)*p18 &
&                  -r2(6,33)*p18-r2(6,34)*p18)*qz+(+r6(25,10)+r6(25,11)*four+r6(25,12) &
&                  -r5(18,18)-r5(18,19)*four-r5(18,20)+r4(14,27)*three+r4(14,28)*p12 &
&                  +r4(14,29)*three+r4(12,39)+r4(12,40)*four+r4(12,41)-r3(9,24)*three &
&                  -r3(9,25)*p12-r3(9,26)*three-r3(7,49)-r3(7,50)*four-r3(7,51) &
&                  +r2(6,53)*three+r2(6,54)*p12+r2(6,55)*three-r1(2,37)*three-r1(2,38)*p12 &
&                  -r1(2,39)*three)*zz+(+r5(18,23)*two+r5(18,24)*two-r4(12,36)*two &
&                  -r4(12,37)*two+r3(9,41)*six+r3(9,42)*six-r2(6,35)*six-r2(6,36)*six)*zzz &
&                  +rxyz(1)*zzzz
      eri(4,3,2,6)=r112+rxyz(18)*qx+(+r7(24,4)*two-r6(17,8)*two+r5(13,12)*six &
&                  -r4(8,13)*six)*qz+(+r6(24,11)*two-r5(17,19)*two+r4(13,28)*six-r3(8,25)*six &
&                  )*xz+(+r6(17,12)-r5(11,20)+r4(8,29)*three-r3(4,26)*three)*zz+rxyz(2)*xzz
      eri(5,3,2,6)=r103+(+r7(33,3)-r6(25,7)+r5(20,11)*three+r5(18,21)*three &
&                  -r4(14,12)*three-r4(12,34)*three+r3(9,39)*nine-r2(6,33)*nine)*qx+( &
&                  +r7(25,3)+r7(25,4)*two-r6(18,7)-r6(18,8)*two+r5(14,11)*three+r5(14,12)*six &
&                  +r5(12,21)+r5(12,22)*two-r4(9,12)*three-r4(9,13)*six-r4(7,34)-r4(7,35)*two &
&                  +r3(5,39)*three+r3(5,40)*six-r2(4,33)*three-r2(4,34)*six)*qz+(+r6(25,10) &
&                  +r6(25,11)*two-r5(18,18)-r5(18,19)*two+r4(14,27)*three+r4(14,28)*six &
&                  +r4(12,39)+r4(12,40)*two-r3(9,24)*three-r3(9,25)*six-r3(7,49)-r3(7,50)*two &
&                  +r2(6,53)*three+r2(6,54)*six-r1(2,37)*three-r1(2,38)*six)*xz+( &
&                  +r6(18,11)*two+r6(18,12)-r5(12,19)*two-r5(12,20)+r4(9,28)*six &
&                  +r4(9,29)*three-r3(5,25)*six-r3(5,26)*three)*zz+(+r5(18,23)*two+r5(18,24) &
&                  -r4(12,36)*two-r4(12,37)+r3(9,41)*six+r3(9,42)*three-r2(6,35)*six &
&                  -r2(6,36)*three)*xzz+(+r5(12,24)-r4(7,37)+r3(5,42)*three-r2(4,36)*three) &
&                  *zzz+rxyz(1)*xzzz
      eri(6,3,2,6)=r013+(+r7(32,3)+r7(32,4)*two-r6(24,7)-r6(24,8)*two+r5(19,11)*three &
&                  +r5(19,12)*six+r5(17,21)+r5(17,22)*two-r4(13,12)*three-r4(13,13)*six &
&                  -r4(11,34)-r4(11,35)*two+r3(8,39)*three+r3(8,40)*six-r2(2,33)*three &
&                  -r2(2,34)*six)*qz+(+r6(24,11)*two+r6(24,12)-r5(17,19)*two-r5(17,20) &
&                  +r4(13,28)*six+r4(13,29)*three-r3(8,25)*six-r3(8,26)*three)*zz+rxyz(2)*zzz
      eri(1,4,2,6)=r310+(+r7(17,3)*two+r7(17,4)-r6(11,7)*two-r6(11,8)+r5(8,11)*six &
&                  +r5(8,12)*three+r5(17,21)*two+r5(17,22)-r4(4,12)*six-r4(4,13)*three &
&                  -r4(11,34)*two-r4(11,35)+r3(8,39)*six+r3(8,40)*three-r2(2,33)*six &
&                  -r2(2,34)*three)*qx+(+r6(17,10)+r6(17,11)*two-r5(11,18)-r5(11,19)*two &
&                  +r4(8,27)*three+r4(8,28)*six-r3(4,24)*three-r3(4,25)*six)*xx+rxyz(3)*xxx
      eri(2,4,2,6)=r130+rxyz(8)*qx
      eri(3,4,2,6)=r112+rxyz(19)*qx+(+r7(24,3)*two-r6(17,7)*two+r5(13,11)*six &
&                  -r4(8,12)*six)*qz+(+r6(24,11)*two-r5(17,19)*two+r4(13,28)*six-r3(8,25)*six &
&                  )*xz+(+r6(17,10)-r5(11,18)+r4(8,27)*three-r3(4,24)*three)*zz+rxyz(3)*xzz
      eri(4,4,2,6)=r220+(+r7(23,3)+r7(23,4)-r6(16,7)-r6(16,8)+r5(12,11)*three &
&                  +r5(12,12)*three+r5(12,21)+r5(12,22)-r4(7,12)*three-r4(7,13)*three &
&                  -r4(7,34)-r4(7,35)+r3(5,39)*three+r3(5,40)*three-r2(4,33)*three &
&                  -r2(4,34)*three)*qx+rxyz(6)*xx
      eri(5,4,2,6)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(17,10)+r6(17,11) &
&                  -r5(11,18)-r5(11,19)+r4(8,27)*three+r4(8,28)*three-r3(4,24)*three &
&                  -r3(4,25)*three)*xz+rxyz(3)*xxz
      eri(6,4,2,6)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,2,6)=r301+(+r7(18,3)*two+r7(18,4)-r6(12,7)*two-r6(12,8)+r5(9,11)*six &
&                  +r5(9,12)*three+r5(18,21)*two+r5(18,22)-r4(5,12)*six-r4(5,13)*three &
&                  -r4(12,34)*two-r4(12,35)+r3(9,39)*six+r3(9,40)*three-r2(6,33)*six &
&                  -r2(6,34)*three)*qx+(+r7(12,4)-r6(7,8)+r5(5,12)*three+r5(12,22)*three &
&                  -r4(2,13)*three-r4(7,35)*three+r3(5,40)*nine-r2(4,34)*nine)*qz+(+r6(18,10) &
&                  +r6(18,11)*two-r5(12,18)-r5(12,19)*two+r4(9,27)*three+r4(9,28)*six &
&                  -r3(5,24)*three-r3(5,25)*six)*xx+(+r6(12,11)*two+r6(12,12)-r5(7,19)*two &
&                  -r5(7,20)+r4(5,28)*six+r4(5,29)*three+r4(12,40)*two+r4(12,41)-r3(2,25)*six &
&                  -r3(2,26)*three-r3(7,50)*two-r3(7,51)+r2(6,54)*six+r2(6,55)*three &
&                  -r1(2,38)*six-r1(2,39)*three)*xz+(+r5(18,23)-r4(12,36)+r3(9,41)*three &
&                  -r2(6,35)*three)*xxx+(+r5(12,23)+r5(12,24)*two-r4(7,36)-r4(7,37)*two &
&                  +r3(5,41)*three+r3(5,42)*six-r2(4,35)*three-r2(4,36)*six)*xxz+rxyz(1)*xxxz
      eri(2,5,2,6)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,2,6)=r103+(+r7(33,4)-r6(25,8)+r5(20,12)*three+r5(18,22)*three &
&                  -r4(14,13)*three-r4(12,35)*three+r3(9,40)*nine-r2(6,34)*nine)*qx+( &
&                  +r7(25,3)*two+r7(25,4)-r6(18,7)*two-r6(18,8)+r5(14,11)*six+r5(14,12)*three &
&                  +r5(12,21)*two+r5(12,22)-r4(9,12)*six-r4(9,13)*three-r4(7,34)*two-r4(7,35) &
&                  +r3(5,39)*six+r3(5,40)*three-r2(4,33)*six-r2(4,34)*three)*qz+( &
&                  +r6(25,11)*two+r6(25,12)-r5(18,19)*two-r5(18,20)+r4(14,28)*six &
&                  +r4(14,29)*three+r4(12,40)*two+r4(12,41)-r3(9,25)*six-r3(9,26)*three &
&                  -r3(7,50)*two-r3(7,51)+r2(6,54)*six+r2(6,55)*three-r1(2,38)*six &
&                  -r1(2,39)*three)*xz+(+r6(18,10)+r6(18,11)*two-r5(12,18)-r5(12,19)*two &
&                  +r4(9,27)*three+r4(9,28)*six-r3(5,24)*three-r3(5,25)*six)*zz+(+r5(18,23) &
&                  +r5(18,24)*two-r4(12,36)-r4(12,37)*two+r3(9,41)*three+r3(9,42)*six &
&                  -r2(6,35)*three-r2(6,36)*six)*xzz+(+r5(12,23)-r4(7,36)+r3(5,41)*three &
&                  -r2(4,35)*three)*zzz+rxyz(1)*xzzz
      eri(4,5,2,6)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(17,11)+r6(17,12) &
&                  -r5(11,19)-r5(11,20)+r4(8,28)*three+r4(8,29)*three-r3(4,25)*three &
&                  -r3(4,26)*three)*xz+rxyz(2)*xxz
      eri(5,5,2,6)=r202+(+r7(25,3)+r7(25,4)-r6(18,7)-r6(18,8)+r5(14,11)*three &
&                  +r5(14,12)*three+r5(12,21)+r5(12,22)-r4(9,12)*three-r4(9,13)*three &
&                  -r4(7,34)-r4(7,35)+r3(5,39)*three+r3(5,40)*three-r2(4,33)*three &
&                  -r2(4,34)*three)*qx+(+r7(18,3)+r7(18,4)-r6(12,7)-r6(12,8)+r5(9,11)*three &
&                  +r5(9,12)*three+r5(18,21)+r5(18,22)-r4(5,12)*three-r4(5,13)*three &
&                  -r4(12,34)-r4(12,35)+r3(9,39)*three+r3(9,40)*three-r2(6,33)*three &
&                  -r2(6,34)*three)*qz+(+r6(25,11)-r5(18,19)+r4(14,28)*three+r4(12,40) &
&                  -r3(9,25)*three-r3(7,50)+r2(6,54)*three-r1(2,38)*three)*xx+(+r6(18,10) &
&                  +r6(18,11)*two+r6(18,12)-r5(12,18)-r5(12,19)*two-r5(12,20)+r4(9,27)*three &
&                  +r4(9,28)*six+r4(9,29)*three-r3(5,24)*three-r3(5,25)*six-r3(5,26)*three) &
&                  *xz+(+r6(12,11)-r5(7,19)+r4(5,28)*three+r4(12,40)-r3(2,25)*three-r3(7,50) &
&                  +r2(6,54)*three-r1(2,38)*three)*zz+(+r5(18,23)+r5(18,24)-r4(12,36) &
&                  -r4(12,37)+r3(9,41)*three+r3(9,42)*three-r2(6,35)*three-r2(6,36)*three) &
&                  *xxz+(+r5(12,23)+r5(12,24)-r4(7,36)-r4(7,37)+r3(5,41)*three+r3(5,42)*three &
&                  -r2(4,35)*three-r2(4,36)*three)*xzz+rxyz(1)*xxzz
      eri(6,5,2,6)=r112+rxyz(19)*qx+(+r7(24,3)+r7(24,4)-r6(17,7)-r6(17,8)+r5(13,11)*three &
&                  +r5(13,12)*three-r4(8,12)*three-r4(8,13)*three)*qz+(+r6(24,11)+r6(24,12) &
&                  -r5(17,19)-r5(17,20)+r4(13,28)*three+r4(13,29)*three-r3(8,25)*three &
&                  -r3(8,26)*three)*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,2,6)=r211+(+r7(24,3)*two-r6(17,7)*two+r5(13,11)*six-r4(8,12)*six)*qx &
&                  +rxyz(16)*qz+(+r6(24,10)-r5(17,18)+r4(13,27)*three-r3(8,24)*three)*xx+( &
&                  +r6(17,11)*two-r5(11,19)*two+r4(8,28)*six-r3(4,25)*six)*xz+rxyz(3)*xxz
      eri(2,6,2,6)=r031+rxyz(8)*qz
      eri(3,6,2,6)=r013+(+r7(32,3)*two+r7(32,4)-r6(24,7)*two-r6(24,8)+r5(19,11)*six &
&                  +r5(19,12)*three+r5(17,21)*two+r5(17,22)-r4(13,12)*six-r4(13,13)*three &
&                  -r4(11,34)*two-r4(11,35)+r3(8,39)*six+r3(8,40)*three-r2(2,33)*six &
&                  -r2(2,34)*three)*qz+(+r6(24,10)+r6(24,11)*two-r5(17,18)-r5(17,19)*two &
&                  +r4(13,27)*three+r4(13,28)*six-r3(8,24)*three-r3(8,25)*six)*zz+rxyz(3)*zzz
      eri(4,6,2,6)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,2,6)=r112+rxyz(18)*qx+(+r7(24,3)+r7(24,4)-r6(17,7)-r6(17,8)+r5(13,11)*three &
&                  +r5(13,12)*three-r4(8,12)*three-r4(8,13)*three)*qz+(+r6(24,10)+r6(24,11) &
&                  -r5(17,18)-r5(17,19)+r4(13,27)*three+r4(13,28)*three-r3(8,24)*three &
&                  -r3(8,25)*three)*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,2,6)=r022+(+r7(31,3)+r7(31,4)-r6(23,7)-r6(23,8)+r5(18,11)*three &
&                  +r5(18,12)*three+r5(18,21)+r5(18,22)-r4(12,12)*three-r4(12,13)*three &
&                  -r4(12,34)-r4(12,35)+r3(9,39)*three+r3(9,40)*three-r2(6,33)*three &
&                  -r2(6,34)*three)*qz+rxyz(6)*zz
!
      r400= r8(14)-r7(9,1)*two-r7(9,2)+r6(5,1)+r6(5,2)*two+r6(5,4)*three+r6(14,9)*six &
&          -r5(2,1)-r5(2,3)*two-r5(2,4)-r5(9,13)*p12-r5(9,17)*six+r4(5,14)*six &
&          +r4(5,18)*p12+r4(5,26)*p18+r4(14,38)*three-r3(2,11)*six-r3(2,19)*p12 &
&          -r3(2,23)*six-r3(9,43)*six-r3(9,48)*three+r2(6,37)*three+r2(6,42)*six &
&          +r2(6,52)*nine-r1(2,21)*three-r1(2,31)*six-r1(2,36)*three
      r310= r8(19)-r7(13,1)*two-r7(13,2)+r6(8,1)+r6(8,2)*two+r6(8,4)*three+r6(19,9)*three &
&          -r5(4,1)-r5(4,3)*two-r5(4,4)-r5(13,13)*six-r5(13,17)*three+r4(8,14)*three &
&          +r4(8,18)*six+r4(8,26)*nine-r3(4,11)*three-r3(4,19)*six-r3(4,23)*three
      r301= r8(20)-r7(14,1)*two-r7(14,2)+r6(9,1)+r6(9,2)*two+r6(9,4)*three+r6(20,9)*three &
&          -r5(5,1)-r5(5,3)*two-r5(5,4)-r5(14,13)*six-r5(14,17)*three+r4(9,14)*three &
&          +r4(9,18)*six+r4(9,26)*nine-r3(5,11)*three-r3(5,19)*six-r3(5,23)*three
      r220= r8(25)-r7(18,1)*two-r7(18,2)+r6(12,1)+r6(12,2)*two+r6(12,4)*three+r6(14,9) &
&          +r6(25,9)-r5(7,1)-r5(7,3)*two-r5(7,4)-r5(9,13)*two-r5(18,13)*two-r5(9,17) &
&          -r5(18,17)+r4(5,14)+r4(12,14)+r4(5,18)*two+r4(12,18)*two+r4(5,26)*three &
&          +r4(12,26)*three+r4(14,38)-r3(2,11)-r3(7,11)-r3(2,19)*two-r3(7,19)*two-r3(2,23) &
&          -r3(7,23)-r3(9,43)*two-r3(9,48)+r2(6,37)+r2(6,42)*two+r2(6,52)*three-r1(2,21) &
&          -r1(2,31)*two-r1(2,36)
      r211= r8(26)-r7(19,1)*two-r7(19,2)+r6(13,1)+r6(13,2)*two+r6(13,4)*three+r6(26,9) &
&          -r5(8,1)-r5(8,3)*two-r5(8,4)-r5(19,13)*two-r5(19,17)+r4(13,14)+r4(13,18)*two &
&          +r4(13,26)*three-r3(8,11)-r3(8,19)*two-r3(8,23)
      r202= r8(27)-r7(20,1)*two-r7(20,2)+r6(14,1)+r6(14,2)*two+r6(14,4)*three+r6(14,9) &
&          +r6(27,9)-r5(9,1)-r5(9,3)*two-r5(9,4)-r5(9,13)*two-r5(20,13)*two-r5(9,17) &
&          -r5(20,17)+r4(5,14)+r4(14,14)+r4(5,18)*two+r4(14,18)*two+r4(5,26)*three &
&          +r4(14,26)*three+r4(14,38)-r3(2,11)-r3(9,11)-r3(2,19)*two-r3(9,19)*two-r3(2,23) &
&          -r3(9,23)-r3(9,43)*two-r3(9,48)+r2(6,37)+r2(6,42)*two+r2(6,52)*three-r1(2,21) &
&          -r1(2,31)*two-r1(2,36)
      r130= r8(32)-r7(24,1)*two-r7(24,2)+r6(17,1)+r6(17,2)*two+r6(17,4)*three &
&          +r6(19,9)*three-r5(11,1)-r5(11,3)*two-r5(11,4)-r5(13,13)*six-r5(13,17)*three &
&          +r4(8,14)*three+r4(8,18)*six+r4(8,26)*nine-r3(4,11)*three-r3(4,19)*six &
&          -r3(4,23)*three
      r121= r8(33)-r7(25,1)*two-r7(25,2)+r6(18,1)+r6(18,2)*two+r6(18,4)*three+r6(20,9) &
&          -r5(12,1)-r5(12,3)*two-r5(12,4)-r5(14,13)*two-r5(14,17)+r4(9,14)+r4(9,18)*two &
&          +r4(9,26)*three-r3(5,11)-r3(5,19)*two-r3(5,23)
      r112= r8(34)-r7(26,1)*two-r7(26,2)+r6(19,1)+r6(19,2)*two+r6(19,4)*three+r6(19,9) &
&          -r5(13,1)-r5(13,3)*two-r5(13,4)-r5(13,13)*two-r5(13,17)+r4(8,14)+r4(8,18)*two &
&          +r4(8,26)*three-r3(4,11)-r3(4,19)*two-r3(4,23)
      r103= r8(35)-r7(27,1)*two-r7(27,2)+r6(20,1)+r6(20,2)*two+r6(20,4)*three &
&          +r6(20,9)*three-r5(14,1)-r5(14,3)*two-r5(14,4)-r5(14,13)*six-r5(14,17)*three &
&          +r4(9,14)*three+r4(9,18)*six+r4(9,26)*nine-r3(5,11)*three-r3(5,19)*six &
&          -r3(5,23)*three
      r040= r8(40)-r7(31,1)*two-r7(31,2)+r6(23,1)+r6(23,2)*two+r6(23,4)*three+r6(25,9)*six &
&          -r5(16,1)-r5(16,3)*two-r5(16,4)-r5(18,13)*p12-r5(18,17)*six+r4(12,14)*six &
&          +r4(12,18)*p12+r4(12,26)*p18+r4(14,38)*three-r3(7,11)*six-r3(7,19)*p12 &
&          -r3(7,23)*six-r3(9,43)*six-r3(9,48)*three+r2(6,37)*three+r2(6,42)*six &
&          +r2(6,52)*nine-r1(2,21)*three-r1(2,31)*six-r1(2,36)*three
      r031= r8(41)-r7(32,1)*two-r7(32,2)+r6(24,1)+r6(24,2)*two+r6(24,4)*three &
&          +r6(26,9)*three-r5(17,1)-r5(17,3)*two-r5(17,4)-r5(19,13)*six-r5(19,17)*three &
&          +r4(13,14)*three+r4(13,18)*six+r4(13,26)*nine-r3(8,11)*three-r3(8,19)*six &
&          -r3(8,23)*three
      r022= r8(42)-r7(33,1)*two-r7(33,2)+r6(25,1)+r6(25,2)*two+r6(25,4)*three+r6(25,9) &
&          +r6(27,9)-r5(18,1)-r5(18,3)*two-r5(18,4)-r5(18,13)*two-r5(20,13)*two-r5(18,17) &
&          -r5(20,17)+r4(12,14)+r4(14,14)+r4(12,18)*two+r4(14,18)*two+r4(12,26)*three &
&          +r4(14,26)*three+r4(14,38)-r3(7,11)-r3(9,11)-r3(7,19)*two-r3(9,19)*two-r3(7,23) &
&          -r3(9,23)-r3(9,43)*two-r3(9,48)+r2(6,37)+r2(6,42)*two+r2(6,52)*three-r1(2,21) &
&          -r1(2,31)*two-r1(2,36)
      r013= r8(43)-r7(34,1)*two-r7(34,2)+r6(26,1)+r6(26,2)*two+r6(26,4)*three &
&          +r6(26,9)*three-r5(19,1)-r5(19,3)*two-r5(19,4)-r5(19,13)*six-r5(19,17)*three &
&          +r4(13,14)*three+r4(13,18)*six+r4(13,26)*nine-r3(8,11)*three-r3(8,19)*six &
&          -r3(8,23)*three
      r004= r8(44)-r7(35,1)*two-r7(35,2)+r6(27,1)+r6(27,2)*two+r6(27,4)*three+r6(27,9)*six &
&          -r5(20,1)-r5(20,3)*two-r5(20,4)-r5(20,13)*p12-r5(20,17)*six+r4(14,14)*six &
&          +r4(14,18)*p12+r4(14,26)*p18+r4(14,38)*three-r3(9,11)*six-r3(9,19)*p12 &
&          -r3(9,23)*six-r3(9,43)*six-r3(9,48)*three+r2(6,37)*three+r2(6,42)*six &
&          +r2(6,52)*nine-r1(2,21)*three-r1(2,31)*six-r1(2,36)*three
      rxyz(1)=+r4(14,42)-r3(9,47)*two-r3(9,52)+r2(6,41)+r2(6,46)*two+r2(6,56)*three &
&             -r1(2,25)-r1(2,35)*two-r1(2,40)
      rxyz(2)=+r5(19,24)-r4(13,33)*two-r4(13,37)+r3(8,30)+r3(8,34)*two+r3(8,42)*three &
&             -r2(2,24)-r2(2,32)*two-r2(2,36)
      rxyz(3)=+r5(19,23)-r4(13,32)*two-r4(13,36)+r3(8,29)+r3(8,33)*two+r3(8,41)*three &
&             -r2(2,23)-r2(2,31)*two-r2(2,35)
      rxyz(4)=+r6(25,10)-r5(18,14)*two-r5(18,18)+r4(12,15)+r4(12,19)*two+r4(12,27)*three &
&             +r4(14,39)-r3(7,12)-r3(7,20)*two-r3(7,24)-r3(9,44)*two-r3(9,49)+r2(6,38) &
&             +r2(6,43)*two+r2(6,53)*three-r1(2,22)-r1(2,32)*two-r1(2,37)
      rxyz(5)=+r6(25,12)-r5(18,16)*two-r5(18,20)+r4(12,17)+r4(12,21)*two+r4(12,29)*three &
&             +r4(14,41)-r3(7,14)-r3(7,22)*two-r3(7,26)-r3(9,46)*two-r3(9,51)+r2(6,40) &
&             +r2(6,45)*two+r2(6,55)*three-r1(2,24)-r1(2,34)*two-r1(2,39)
      rxyz(6)=+r6(25,11)-r5(18,15)*two-r5(18,19)+r4(12,16)+r4(12,20)*two+r4(12,28)*three &
&             +r4(14,40)-r3(7,13)-r3(7,21)*two-r3(7,25)-r3(9,45)*two-r3(9,50)+r2(6,39) &
&             +r2(6,44)*two+r2(6,54)*three-r1(2,23)-r1(2,33)*two-r1(2,38)
      rxyz(7)=+r7(32,3)-r6(24,5)*two-r6(24,7)+r5(17,5)+r5(17,7)*two+r5(17,11)*three &
&             +r5(19,21)*three-r4(11,6)-r4(11,10)*two-r4(11,12)-r4(13,30)*six &
&             -r4(13,34)*three+r3(8,27)*three+r3(8,31)*six+r3(8,39)*nine-r2(2,21)*three &
&             -r2(2,29)*six-r2(2,33)*three
      rxyz(8)=+r7(32,4)-r6(24,6)*two-r6(24,8)+r5(17,6)+r5(17,8)*two+r5(17,12)*three &
&             +r5(19,22)*three-r4(11,7)-r4(11,11)*two-r4(11,13)-r4(13,31)*six &
&             -r4(13,35)*three+r3(8,28)*three+r3(8,32)*six+r3(8,40)*nine-r2(2,22)*three &
&             -r2(2,30)*six-r2(2,34)*three
      rxyz(9)=+r7(26,3)+r7(26,4)-r6(19,5)*two-r6(19,6)*two-r6(19,7)-r6(19,8)+r5(13,5) &
&             +r5(13,6)+r5(13,7)*two+r5(13,8)*two+r5(13,11)*three+r5(13,12)*three-r4(8,6) &
&             -r4(8,7)-r4(8,10)*two-r4(8,11)*two-r4(8,12)-r4(8,13)
      rxyz(10)=+r6(26,11)-r5(19,15)*two-r5(19,19)+r4(13,16)+r4(13,20)*two+r4(13,28)*three &
&             -r3(8,13)-r3(8,21)*two-r3(8,25)
      rxyz(11)=+r6(19,11)-r5(13,15)*two-r5(13,19)+r4(8,16)+r4(8,20)*two+r4(8,28)*three &
&             -r3(4,13)-r3(4,21)*two-r3(4,25)
      rxyz(12)=+r7(33,3)-r6(25,5)*two-r6(25,7)+r5(18,5)+r5(18,7)*two+r5(18,11)*three &
&             +r5(20,21)-r4(12,6)-r4(12,10)*two-r4(12,12)-r4(14,30)*two-r4(14,34)+r3(9,27) &
&             +r3(9,31)*two+r3(9,39)*three-r2(6,21)-r2(6,29)*two-r2(6,33)
      rxyz(13)=+r7(33,4)-r6(25,6)*two-r6(25,8)+r5(18,6)+r5(18,8)*two+r5(18,12)*three &
&             +r5(20,22)-r4(12,7)-r4(12,11)*two-r4(12,13)-r4(14,31)*two-r4(14,35)+r3(9,28) &
&             +r3(9,32)*two+r3(9,40)*three-r2(6,22)-r2(6,30)*two-r2(6,34)
      rxyz(14)=+r7(25,4)-r6(18,6)*two-r6(18,8)+r5(12,6)+r5(12,8)*two+r5(12,12)*three &
&             +r5(14,22)-r4(7,7)-r4(7,11)*two-r4(7,13)-r4(9,31)*two-r4(9,35)+r3(5,28) &
&             +r3(5,32)*two+r3(5,40)*three-r2(4,22)-r2(4,30)*two-r2(4,34)
      rxyz(15)=+r7(25,3)-r6(18,5)*two-r6(18,7)+r5(12,5)+r5(12,7)*two+r5(12,11)*three &
&             +r5(14,21)-r4(7,6)-r4(7,10)*two-r4(7,12)-r4(9,30)*two-r4(9,34)+r3(5,27) &
&             +r3(5,31)*two+r3(5,39)*three-r2(4,21)-r2(4,29)*two-r2(4,33)
      rxyz(16)=+r7(19,4)-r6(13,6)*two-r6(13,8)+r5(8,6)+r5(8,8)*two+r5(8,12)*three &
&             +r5(19,22)-r4(4,7)-r4(4,11)*two-r4(4,13)-r4(13,31)*two-r4(13,35)+r3(8,28) &
&             +r3(8,32)*two+r3(8,40)*three-r2(2,22)-r2(2,30)*two-r2(2,34)
      rxyz(17)=+r7(19,3)-r6(13,5)*two-r6(13,7)+r5(8,5)+r5(8,7)*two+r5(8,11)*three &
&             +r5(19,21)-r4(4,6)-r4(4,10)*two-r4(4,12)-r4(13,30)*two-r4(13,34)+r3(8,27) &
&             +r3(8,31)*two+r3(8,39)*three-r2(2,21)-r2(2,29)*two-r2(2,33)
      rxyz(18)=+r7(34,3)-r6(26,5)*two-r6(26,7)+r5(19,5)+r5(19,7)*two+r5(19,11)*three &
&             +r5(19,21)-r4(13,6)-r4(13,10)*two-r4(13,12)-r4(13,30)*two-r4(13,34)+r3(8,27) &
&             +r3(8,31)*two+r3(8,39)*three-r2(2,21)-r2(2,29)*two-r2(2,33)
      rxyz(19)=+r7(34,4)-r6(26,6)*two-r6(26,8)+r5(19,6)+r5(19,8)*two+r5(19,12)*three &
&             +r5(19,22)-r4(13,7)-r4(13,11)*two-r4(13,13)-r4(13,31)*two-r4(13,35)+r3(8,28) &
&             +r3(8,32)*two+r3(8,40)*three-r2(2,22)-r2(2,30)*two-r2(2,34)
      rxyz(20)=+r6(20,11)*four-r5(14,15)*eight-r5(14,19)*four+r4(9,16)*four+r4(9,20)*eight &
&             +r4(9,28)*p12-r3(5,13)*four-r3(5,21)*eight-r3(5,25)*four
      eri(1,1,3,6)=r400+(+r7(14,3)*two+r7(14,4)*two-r6(9,5)*four-r6(9,6)*four-r6(9,7)*two &
&                  -r6(9,8)*two+r5(5,5)*two+r5(5,6)*two+r5(5,7)*four+r5(5,8)*four &
&                  +r5(5,11)*six+r5(5,12)*six+r5(14,21)*six+r5(14,22)*six-r4(2,6)*two &
&                  -r4(2,7)*two-r4(2,10)*four-r4(2,11)*four-r4(2,12)*two-r4(2,13)*two &
&                  -r4(9,30)*p12-r4(9,31)*p12-r4(9,34)*six-r4(9,35)*six+r3(5,27)*six &
&                  +r3(5,28)*six+r3(5,31)*p12+r3(5,32)*p12+r3(5,39)*p18+r3(5,40)*p18 &
&                  -r2(4,21)*six-r2(4,22)*six-r2(4,29)*p12-r2(4,30)*p12-r2(4,33)*six &
&                  -r2(4,34)*six)*qx+(+r6(14,10)+r6(14,11)*four+r6(14,12)-r5(9,14)*two &
&                  -r5(9,15)*eight-r5(9,16)*two-r5(9,18)-r5(9,19)*four-r5(9,20)+r4(5,15) &
&                  +r4(5,16)*four+r4(5,17)+r4(5,19)*two+r4(5,20)*eight+r4(5,21)*two &
&                  +r4(5,27)*three+r4(5,28)*p12+r4(5,29)*three+r4(14,39)+r4(14,40)*four &
&                  +r4(14,41)-r3(2,12)-r3(2,13)*four-r3(2,14)-r3(2,20)*two-r3(2,21)*eight &
&                  -r3(2,22)*two-r3(2,24)-r3(2,25)*four-r3(2,26)-r3(9,44)*two-r3(9,45)*eight &
&                  -r3(9,46)*two-r3(9,49)-r3(9,50)*four-r3(9,51)+r2(6,38)+r2(6,39)*four &
&                  +r2(6,40)+r2(6,43)*two+r2(6,44)*eight+r2(6,45)*two+r2(6,53)*three &
&                  +r2(6,54)*p12+r2(6,55)*three-r1(2,22)-r1(2,23)*four-r1(2,24)-r1(2,32)*two &
&                  -r1(2,33)*eight-r1(2,34)*two-r1(2,37)-r1(2,38)*four-r1(2,39))*xx+( &
&                  +r5(14,23)*two+r5(14,24)*two-r4(9,32)*four-r4(9,33)*four-r4(9,36)*two &
&                  -r4(9,37)*two+r3(5,29)*two+r3(5,30)*two+r3(5,33)*four+r3(5,34)*four &
&                  +r3(5,41)*six+r3(5,42)*six-r2(4,23)*two-r2(4,24)*two-r2(4,31)*four &
&                  -r2(4,32)*four-r2(4,35)*two-r2(4,36)*two)*xxx+rxyz(1)*xxxx
      eri(2,1,3,6)=r220+(+r7(25,4)*two-r6(18,6)*four-r6(18,8)*two+r5(12,6)*two &
&                  +r5(12,8)*four+r5(12,12)*six+r5(14,22)*two-r4(7,7)*two-r4(7,11)*four &
&                  -r4(7,13)*two-r4(9,31)*four-r4(9,35)*two+r3(5,28)*two+r3(5,32)*four &
&                  +r3(5,40)*six-r2(4,22)*two-r2(4,30)*four-r2(4,34)*two)*qx+rxyz(5)*xx
      eri(3,1,3,6)=r202+(+r7(27,4)*two-r6(20,6)*four-r6(20,8)*two+r5(14,6)*two &
&                  +r5(14,8)*four+r5(14,12)*six+r5(14,22)*two-r4(9,7)*two-r4(9,11)*four &
&                  -r4(9,13)*two-r4(9,31)*four-r4(9,35)*two+r3(5,28)*two+r3(5,32)*four &
&                  +r3(5,40)*six-r2(4,22)*two-r2(4,30)*four-r2(4,34)*two)*qx+(+r7(20,3)*two &
&                  -r6(14,5)*four-r6(14,7)*two+r5(9,5)*two+r5(9,7)*four+r5(9,11)*six &
&                  +r5(20,21)*two-r4(5,6)*two-r4(5,10)*four-r4(5,12)*two-r4(14,30)*four &
&                  -r4(14,34)*two+r3(9,27)*two+r3(9,31)*four+r3(9,39)*six-r2(6,21)*two &
&                  -r2(6,29)*four-r2(6,33)*two)*qz+(+r6(27,12)-r5(20,16)*two-r5(20,20) &
&                  +r4(14,17)+r4(14,21)*two+r4(14,29)*three+r4(14,41)-r3(9,14)-r3(9,22)*two &
&                  -r3(9,26)-r3(9,46)*two-r3(9,51)+r2(6,40)+r2(6,45)*two+r2(6,55)*three &
&                  -r1(2,24)-r1(2,34)*two-r1(2,39))*xx+rxyz(20)*xz+(+r6(14,10)-r5(9,14)*two &
&                  -r5(9,18)+r4(5,15)+r4(5,19)*two+r4(5,27)*three+r4(14,39)-r3(2,12) &
&                  -r3(2,20)*two-r3(2,24)-r3(9,44)*two-r3(9,49)+r2(6,38)+r2(6,43)*two &
&                  +r2(6,53)*three-r1(2,22)-r1(2,32)*two-r1(2,37))*zz+(+r5(20,24)*two &
&                  -r4(14,33)*four-r4(14,37)*two+r3(9,30)*two+r3(9,34)*four+r3(9,42)*six &
&                  -r2(6,24)*two-r2(6,32)*four-r2(6,36)*two)*xxz+(+r5(14,23)*two &
&                  -r4(9,32)*four-r4(9,36)*two+r3(5,29)*two+r3(5,33)*four+r3(5,41)*six &
&                  -r2(4,23)*two-r2(4,31)*four-r2(4,35)*two)*xzz+rxyz(1)*xxzz
      eri(4,1,3,6)=r310+(+r7(19,3)+r7(19,4)*two-r6(13,5)*two-r6(13,6)*four-r6(13,7) &
&                  -r6(13,8)*two+r5(8,5)+r5(8,6)*two+r5(8,7)*two+r5(8,8)*four+r5(8,11)*three &
&                  +r5(8,12)*six+r5(19,21)+r5(19,22)*two-r4(4,6)-r4(4,7)*two-r4(4,10)*two &
&                  -r4(4,11)*four-r4(4,12)-r4(4,13)*two-r4(13,30)*two-r4(13,31)*four &
&                  -r4(13,34)-r4(13,35)*two+r3(8,27)+r3(8,28)*two+r3(8,31)*two+r3(8,32)*four &
&                  +r3(8,39)*three+r3(8,40)*six-r2(2,21)-r2(2,22)*two-r2(2,29)*two &
&                  -r2(2,30)*four-r2(2,33)-r2(2,34)*two)*qx+(+r6(19,11)*two+r6(19,12) &
&                  -r5(13,15)*four-r5(13,16)*two-r5(13,19)*two-r5(13,20)+r4(8,16)*two &
&                  +r4(8,17)+r4(8,20)*four+r4(8,21)*two+r4(8,28)*six+r4(8,29)*three &
&                  -r3(4,13)*two-r3(4,14)-r3(4,21)*four-r3(4,22)*two-r3(4,25)*two-r3(4,26)) &
&                  *xx+rxyz(2)*xxx
      eri(5,1,3,6)=r301+(+r7(20,3)+r7(20,4)*two-r6(14,5)*two-r6(14,6)*four-r6(14,7) &
&                  -r6(14,8)*two+r5(9,5)+r5(9,6)*two+r5(9,7)*two+r5(9,8)*four+r5(9,11)*three &
&                  +r5(9,12)*six+r5(20,21)+r5(20,22)*two-r4(5,6)-r4(5,7)*two-r4(5,10)*two &
&                  -r4(5,11)*four-r4(5,12)-r4(5,13)*two-r4(14,30)*two-r4(14,31)*four &
&                  -r4(14,34)-r4(14,35)*two+r3(9,27)+r3(9,28)*two+r3(9,31)*two+r3(9,32)*four &
&                  +r3(9,39)*three+r3(9,40)*six-r2(6,21)-r2(6,22)*two-r2(6,29)*two &
&                  -r2(6,30)*four-r2(6,33)-r2(6,34)*two)*qx+(+r7(14,3)-r6(9,5)*two-r6(9,7) &
&                  +r5(5,5)+r5(5,7)*two+r5(5,11)*three+r5(14,21)*three-r4(2,6)-r4(2,10)*two &
&                  -r4(2,12)-r4(9,30)*six-r4(9,34)*three+r3(5,27)*three+r3(5,31)*six &
&                  +r3(5,39)*nine-r2(4,21)*three-r2(4,29)*six-r2(4,33)*three)*qz+( &
&                  +r6(20,11)*two+r6(20,12)-r5(14,15)*four-r5(14,16)*two-r5(14,19)*two &
&                  -r5(14,20)+r4(9,16)*two+r4(9,17)+r4(9,20)*four+r4(9,21)*two+r4(9,28)*six &
&                  +r4(9,29)*three-r3(5,13)*two-r3(5,14)-r3(5,21)*four-r3(5,22)*two &
&                  -r3(5,25)*two-r3(5,26))*xx+(+r6(14,10)+r6(14,11)*two-r5(9,14)*two &
&                  -r5(9,15)*four-r5(9,18)-r5(9,19)*two+r4(5,15)+r4(5,16)*two+r4(5,19)*two &
&                  +r4(5,20)*four+r4(5,27)*three+r4(5,28)*six+r4(14,39)+r4(14,40)*two &
&                  -r3(2,12)-r3(2,13)*two-r3(2,20)*two-r3(2,21)*four-r3(2,24)-r3(2,25)*two &
&                  -r3(9,44)*two-r3(9,45)*four-r3(9,49)-r3(9,50)*two+r2(6,38)+r2(6,39)*two &
&                  +r2(6,43)*two+r2(6,44)*four+r2(6,53)*three+r2(6,54)*six-r1(2,22) &
&                  -r1(2,23)*two-r1(2,32)*two-r1(2,33)*four-r1(2,37)-r1(2,38)*two)*xz+( &
&                  +r5(20,24)-r4(14,33)*two-r4(14,37)+r3(9,30)+r3(9,34)*two+r3(9,42)*three &
&                  -r2(6,24)-r2(6,32)*two-r2(6,36))*xxx+(+r5(14,23)*two+r5(14,24) &
&                  -r4(9,32)*four-r4(9,33)*two-r4(9,36)*two-r4(9,37)+r3(5,29)*two+r3(5,30) &
&                  +r3(5,33)*four+r3(5,34)*two+r3(5,41)*six+r3(5,42)*three-r2(4,23)*two &
&                  -r2(4,24)-r2(4,31)*four-r2(4,32)*two-r2(4,35)*two-r2(4,36))*xxz+rxyz(1) &
&                  *xxxz
      eri(6,1,3,6)=r211+(+r7(26,4)*two-r6(19,6)*four-r6(19,8)*two+r5(13,6)*two &
&                  +r5(13,8)*four+r5(13,12)*six-r4(8,7)*two-r4(8,11)*four-r4(8,13)*two)*qx &
&                  +rxyz(17)*qz+(+r6(26,12)-r5(19,16)*two-r5(19,20)+r4(13,17)+r4(13,21)*two &
&                  +r4(13,29)*three-r3(8,14)-r3(8,22)*two-r3(8,26))*xx+(+r6(19,11)*two &
&                  -r5(13,15)*four-r5(13,19)*two+r4(8,16)*two+r4(8,20)*four+r4(8,28)*six &
&                  -r3(4,13)*two-r3(4,21)*four-r3(4,25)*two)*xz+rxyz(2)*xxz
      eri(1,2,3,6)=r220+(+r7(25,3)*two-r6(18,5)*four-r6(18,7)*two+r5(12,5)*two &
&                  +r5(12,7)*four+r5(12,11)*six+r5(14,21)*two-r4(7,6)*two-r4(7,10)*four &
&                  -r4(7,12)*two-r4(9,30)*four-r4(9,34)*two+r3(5,27)*two+r3(5,31)*four &
&                  +r3(5,39)*six-r2(4,21)*two-r2(4,29)*four-r2(4,33)*two)*qx+rxyz(4)*xx
      eri(2,2,3,6)=r040
      eri(3,2,3,6)=r022+(+r7(33,3)*two-r6(25,5)*four-r6(25,7)*two+r5(18,5)*two &
&                  +r5(18,7)*four+r5(18,11)*six+r5(20,21)*two-r4(12,6)*two-r4(12,10)*four &
&                  -r4(12,12)*two-r4(14,30)*four-r4(14,34)*two+r3(9,27)*two+r3(9,31)*four &
&                  +r3(9,39)*six-r2(6,21)*two-r2(6,29)*four-r2(6,33)*two)*qz+rxyz(4)*zz
      eri(4,2,3,6)=r130+rxyz(7)*qx
      eri(5,2,3,6)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,3,6)=r031+rxyz(7)*qz
      eri(1,3,3,6)=r202+(+r7(27,3)*two-r6(20,5)*four-r6(20,7)*two+r5(14,5)*two &
&                  +r5(14,7)*four+r5(14,11)*six+r5(14,21)*two-r4(9,6)*two-r4(9,10)*four &
&                  -r4(9,12)*two-r4(9,30)*four-r4(9,34)*two+r3(5,27)*two+r3(5,31)*four &
&                  +r3(5,39)*six-r2(4,21)*two-r2(4,29)*four-r2(4,33)*two)*qx+(+r7(20,4)*two &
&                  -r6(14,6)*four-r6(14,8)*two+r5(9,6)*two+r5(9,8)*four+r5(9,12)*six &
&                  +r5(20,22)*two-r4(5,7)*two-r4(5,11)*four-r4(5,13)*two-r4(14,31)*four &
&                  -r4(14,35)*two+r3(9,28)*two+r3(9,32)*four+r3(9,40)*six-r2(6,22)*two &
&                  -r2(6,30)*four-r2(6,34)*two)*qz+(+r6(27,10)-r5(20,14)*two-r5(20,18) &
&                  +r4(14,15)+r4(14,19)*two+r4(14,27)*three+r4(14,39)-r3(9,12)-r3(9,20)*two &
&                  -r3(9,24)-r3(9,44)*two-r3(9,49)+r2(6,38)+r2(6,43)*two+r2(6,53)*three &
&                  -r1(2,22)-r1(2,32)*two-r1(2,37))*xx+rxyz(20)*xz+(+r6(14,12)-r5(9,16)*two &
&                  -r5(9,20)+r4(5,17)+r4(5,21)*two+r4(5,29)*three+r4(14,41)-r3(2,14) &
&                  -r3(2,22)*two-r3(2,26)-r3(9,46)*two-r3(9,51)+r2(6,40)+r2(6,45)*two &
&                  +r2(6,55)*three-r1(2,24)-r1(2,34)*two-r1(2,39))*zz+(+r5(20,23)*two &
&                  -r4(14,32)*four-r4(14,36)*two+r3(9,29)*two+r3(9,33)*four+r3(9,41)*six &
&                  -r2(6,23)*two-r2(6,31)*four-r2(6,35)*two)*xxz+(+r5(14,24)*two &
&                  -r4(9,33)*four-r4(9,37)*two+r3(5,30)*two+r3(5,34)*four+r3(5,42)*six &
&                  -r2(4,24)*two-r2(4,32)*four-r2(4,36)*two)*xzz+rxyz(1)*xxzz
      eri(2,3,3,6)=r022+(+r7(33,4)*two-r6(25,6)*four-r6(25,8)*two+r5(18,6)*two &
&                  +r5(18,8)*four+r5(18,12)*six+r5(20,22)*two-r4(12,7)*two-r4(12,11)*four &
&                  -r4(12,13)*two-r4(14,31)*four-r4(14,35)*two+r3(9,28)*two+r3(9,32)*four &
&                  +r3(9,40)*six-r2(6,22)*two-r2(6,30)*four-r2(6,34)*two)*qz+rxyz(5)*zz
      eri(3,3,3,6)=r004+(+r7(35,3)*two+r7(35,4)*two-r6(27,5)*four-r6(27,6)*four &
&                  -r6(27,7)*two-r6(27,8)*two+r5(20,5)*two+r5(20,6)*two+r5(20,7)*four &
&                  +r5(20,8)*four+r5(20,11)*six+r5(20,12)*six+r5(20,21)*six+r5(20,22)*six &
&                  -r4(14,6)*two-r4(14,7)*two-r4(14,10)*four-r4(14,11)*four-r4(14,12)*two &
&                  -r4(14,13)*two-r4(14,30)*p12-r4(14,31)*p12-r4(14,34)*six-r4(14,35)*six &
&                  +r3(9,27)*six+r3(9,28)*six+r3(9,31)*p12+r3(9,32)*p12+r3(9,39)*p18 &
&                  +r3(9,40)*p18-r2(6,21)*six-r2(6,22)*six-r2(6,29)*p12-r2(6,30)*p12 &
&                  -r2(6,33)*six-r2(6,34)*six)*qz+(+r6(27,10)+r6(27,11)*four+r6(27,12) &
&                  -r5(20,14)*two-r5(20,15)*eight-r5(20,16)*two-r5(20,18)-r5(20,19)*four &
&                  -r5(20,20)+r4(14,15)+r4(14,16)*four+r4(14,17)+r4(14,19)*two &
&                  +r4(14,20)*eight+r4(14,21)*two+r4(14,27)*three+r4(14,28)*p12 &
&                  +r4(14,29)*three+r4(14,39)+r4(14,40)*four+r4(14,41)-r3(9,12)-r3(9,13)*four &
&                  -r3(9,14)-r3(9,20)*two-r3(9,21)*eight-r3(9,22)*two-r3(9,24)-r3(9,25)*four &
&                  -r3(9,26)-r3(9,44)*two-r3(9,45)*eight-r3(9,46)*two-r3(9,49)-r3(9,50)*four &
&                  -r3(9,51)+r2(6,38)+r2(6,39)*four+r2(6,40)+r2(6,43)*two+r2(6,44)*eight &
&                  +r2(6,45)*two+r2(6,53)*three+r2(6,54)*p12+r2(6,55)*three-r1(2,22) &
&                  -r1(2,23)*four-r1(2,24)-r1(2,32)*two-r1(2,33)*eight-r1(2,34)*two-r1(2,37) &
&                  -r1(2,38)*four-r1(2,39))*zz+(+r5(20,23)*two+r5(20,24)*two-r4(14,32)*four &
&                  -r4(14,33)*four-r4(14,36)*two-r4(14,37)*two+r3(9,29)*two+r3(9,30)*two &
&                  +r3(9,33)*four+r3(9,34)*four+r3(9,41)*six+r3(9,42)*six-r2(6,23)*two &
&                  -r2(6,24)*two-r2(6,31)*four-r2(6,32)*four-r2(6,35)*two-r2(6,36)*two)*zzz &
&                  +rxyz(1)*zzzz
      eri(4,3,3,6)=r112+rxyz(18)*qx+(+r7(26,4)*two-r6(19,6)*four-r6(19,8)*two &
&                  +r5(13,6)*two+r5(13,8)*four+r5(13,12)*six-r4(8,7)*two-r4(8,11)*four &
&                  -r4(8,13)*two)*qz+(+r6(26,11)*two-r5(19,15)*four-r5(19,19)*two &
&                  +r4(13,16)*two+r4(13,20)*four+r4(13,28)*six-r3(8,13)*two-r3(8,21)*four &
&                  -r3(8,25)*two)*xz+(+r6(19,12)-r5(13,16)*two-r5(13,20)+r4(8,17) &
&                  +r4(8,21)*two+r4(8,29)*three-r3(4,14)-r3(4,22)*two-r3(4,26))*zz+rxyz(2) &
&                  *xzz
      eri(5,3,3,6)=r103+(+r7(35,3)-r6(27,5)*two-r6(27,7)+r5(20,5)+r5(20,7)*two &
&                  +r5(20,11)*three+r5(20,21)*three-r4(14,6)-r4(14,10)*two-r4(14,12) &
&                  -r4(14,30)*six-r4(14,34)*three+r3(9,27)*three+r3(9,31)*six+r3(9,39)*nine &
&                  -r2(6,21)*three-r2(6,29)*six-r2(6,33)*three)*qx+(+r7(27,3)+r7(27,4)*two &
&                  -r6(20,5)*two-r6(20,6)*four-r6(20,7)-r6(20,8)*two+r5(14,5)+r5(14,6)*two &
&                  +r5(14,7)*two+r5(14,8)*four+r5(14,11)*three+r5(14,12)*six+r5(14,21) &
&                  +r5(14,22)*two-r4(9,6)-r4(9,7)*two-r4(9,10)*two-r4(9,11)*four-r4(9,12) &
&                  -r4(9,13)*two-r4(9,30)*two-r4(9,31)*four-r4(9,34)-r4(9,35)*two+r3(5,27) &
&                  +r3(5,28)*two+r3(5,31)*two+r3(5,32)*four+r3(5,39)*three+r3(5,40)*six &
&                  -r2(4,21)-r2(4,22)*two-r2(4,29)*two-r2(4,30)*four-r2(4,33)-r2(4,34)*two) &
&                  *qz+(+r6(27,10)+r6(27,11)*two-r5(20,14)*two-r5(20,15)*four-r5(20,18) &
&                  -r5(20,19)*two+r4(14,15)+r4(14,16)*two+r4(14,19)*two+r4(14,20)*four &
&                  +r4(14,27)*three+r4(14,28)*six+r4(14,39)+r4(14,40)*two-r3(9,12) &
&                  -r3(9,13)*two-r3(9,20)*two-r3(9,21)*four-r3(9,24)-r3(9,25)*two &
&                  -r3(9,44)*two-r3(9,45)*four-r3(9,49)-r3(9,50)*two+r2(6,38)+r2(6,39)*two &
&                  +r2(6,43)*two+r2(6,44)*four+r2(6,53)*three+r2(6,54)*six-r1(2,22) &
&                  -r1(2,23)*two-r1(2,32)*two-r1(2,33)*four-r1(2,37)-r1(2,38)*two)*xz+( &
&                  +r6(20,11)*two+r6(20,12)-r5(14,15)*four-r5(14,16)*two-r5(14,19)*two &
&                  -r5(14,20)+r4(9,16)*two+r4(9,17)+r4(9,20)*four+r4(9,21)*two+r4(9,28)*six &
&                  +r4(9,29)*three-r3(5,13)*two-r3(5,14)-r3(5,21)*four-r3(5,22)*two &
&                  -r3(5,25)*two-r3(5,26))*zz+(+r5(20,23)*two+r5(20,24)-r4(14,32)*four &
&                  -r4(14,33)*two-r4(14,36)*two-r4(14,37)+r3(9,29)*two+r3(9,30)+r3(9,33)*four &
&                  +r3(9,34)*two+r3(9,41)*six+r3(9,42)*three-r2(6,23)*two-r2(6,24) &
&                  -r2(6,31)*four-r2(6,32)*two-r2(6,35)*two-r2(6,36))*xzz+(+r5(14,24) &
&                  -r4(9,33)*two-r4(9,37)+r3(5,30)+r3(5,34)*two+r3(5,42)*three-r2(4,24) &
&                  -r2(4,32)*two-r2(4,36))*zzz+rxyz(1)*xzzz
      eri(6,3,3,6)=r013+(+r7(34,3)+r7(34,4)*two-r6(26,5)*two-r6(26,6)*four-r6(26,7) &
&                  -r6(26,8)*two+r5(19,5)+r5(19,6)*two+r5(19,7)*two+r5(19,8)*four &
&                  +r5(19,11)*three+r5(19,12)*six+r5(19,21)+r5(19,22)*two-r4(13,6) &
&                  -r4(13,7)*two-r4(13,10)*two-r4(13,11)*four-r4(13,12)-r4(13,13)*two &
&                  -r4(13,30)*two-r4(13,31)*four-r4(13,34)-r4(13,35)*two+r3(8,27) &
&                  +r3(8,28)*two+r3(8,31)*two+r3(8,32)*four+r3(8,39)*three+r3(8,40)*six &
&                  -r2(2,21)-r2(2,22)*two-r2(2,29)*two-r2(2,30)*four-r2(2,33)-r2(2,34)*two) &
&                  *qz+(+r6(26,11)*two+r6(26,12)-r5(19,15)*four-r5(19,16)*two-r5(19,19)*two &
&                  -r5(19,20)+r4(13,16)*two+r4(13,17)+r4(13,20)*four+r4(13,21)*two &
&                  +r4(13,28)*six+r4(13,29)*three-r3(8,13)*two-r3(8,14)-r3(8,21)*four &
&                  -r3(8,22)*two-r3(8,25)*two-r3(8,26))*zz+rxyz(2)*zzz
      eri(1,4,3,6)=r310+(+r7(19,3)*two+r7(19,4)-r6(13,5)*four-r6(13,6)*two-r6(13,7)*two &
&                  -r6(13,8)+r5(8,5)*two+r5(8,6)+r5(8,7)*four+r5(8,8)*two+r5(8,11)*six &
&                  +r5(8,12)*three+r5(19,21)*two+r5(19,22)-r4(4,6)*two-r4(4,7)-r4(4,10)*four &
&                  -r4(4,11)*two-r4(4,12)*two-r4(4,13)-r4(13,30)*four-r4(13,31)*two &
&                  -r4(13,34)*two-r4(13,35)+r3(8,27)*two+r3(8,28)+r3(8,31)*four+r3(8,32)*two &
&                  +r3(8,39)*six+r3(8,40)*three-r2(2,21)*two-r2(2,22)-r2(2,29)*four &
&                  -r2(2,30)*two-r2(2,33)*two-r2(2,34))*qx+(+r6(19,10)+r6(19,11)*two &
&                  -r5(13,14)*two-r5(13,15)*four-r5(13,18)-r5(13,19)*two+r4(8,15) &
&                  +r4(8,16)*two+r4(8,19)*two+r4(8,20)*four+r4(8,27)*three+r4(8,28)*six &
&                  -r3(4,12)-r3(4,13)*two-r3(4,20)*two-r3(4,21)*four-r3(4,24)-r3(4,25)*two) &
&                  *xx+rxyz(3)*xxx
      eri(2,4,3,6)=r130+rxyz(8)*qx
      eri(3,4,3,6)=r112+rxyz(19)*qx+(+r7(26,3)*two-r6(19,5)*four-r6(19,7)*two &
&                  +r5(13,5)*two+r5(13,7)*four+r5(13,11)*six-r4(8,6)*two-r4(8,10)*four &
&                  -r4(8,12)*two)*qz+(+r6(26,11)*two-r5(19,15)*four-r5(19,19)*two &
&                  +r4(13,16)*two+r4(13,20)*four+r4(13,28)*six-r3(8,13)*two-r3(8,21)*four &
&                  -r3(8,25)*two)*xz+(+r6(19,10)-r5(13,14)*two-r5(13,18)+r4(8,15) &
&                  +r4(8,19)*two+r4(8,27)*three-r3(4,12)-r3(4,20)*two-r3(4,24))*zz+rxyz(3) &
&                  *xzz
      eri(4,4,3,6)=r220+(+r7(25,3)+r7(25,4)-r6(18,5)*two-r6(18,6)*two-r6(18,7)-r6(18,8) &
&                  +r5(12,5)+r5(12,6)+r5(12,7)*two+r5(12,8)*two+r5(12,11)*three &
&                  +r5(12,12)*three+r5(14,21)+r5(14,22)-r4(7,6)-r4(7,7)-r4(7,10)*two &
&                  -r4(7,11)*two-r4(7,12)-r4(7,13)-r4(9,30)*two-r4(9,31)*two-r4(9,34) &
&                  -r4(9,35)+r3(5,27)+r3(5,28)+r3(5,31)*two+r3(5,32)*two+r3(5,39)*three &
&                  +r3(5,40)*three-r2(4,21)-r2(4,22)-r2(4,29)*two-r2(4,30)*two-r2(4,33) &
&                  -r2(4,34))*qx+rxyz(6)*xx
      eri(5,4,3,6)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(19,10)+r6(19,11) &
&                  -r5(13,14)*two-r5(13,15)*two-r5(13,18)-r5(13,19)+r4(8,15)+r4(8,16) &
&                  +r4(8,19)*two+r4(8,20)*two+r4(8,27)*three+r4(8,28)*three-r3(4,12)-r3(4,13) &
&                  -r3(4,20)*two-r3(4,21)*two-r3(4,24)-r3(4,25))*xz+rxyz(3)*xxz
      eri(6,4,3,6)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,3,6)=r301+(+r7(20,3)*two+r7(20,4)-r6(14,5)*four-r6(14,6)*two-r6(14,7)*two &
&                  -r6(14,8)+r5(9,5)*two+r5(9,6)+r5(9,7)*four+r5(9,8)*two+r5(9,11)*six &
&                  +r5(9,12)*three+r5(20,21)*two+r5(20,22)-r4(5,6)*two-r4(5,7)-r4(5,10)*four &
&                  -r4(5,11)*two-r4(5,12)*two-r4(5,13)-r4(14,30)*four-r4(14,31)*two &
&                  -r4(14,34)*two-r4(14,35)+r3(9,27)*two+r3(9,28)+r3(9,31)*four+r3(9,32)*two &
&                  +r3(9,39)*six+r3(9,40)*three-r2(6,21)*two-r2(6,22)-r2(6,29)*four &
&                  -r2(6,30)*two-r2(6,33)*two-r2(6,34))*qx+(+r7(14,4)-r6(9,6)*two-r6(9,8) &
&                  +r5(5,6)+r5(5,8)*two+r5(5,12)*three+r5(14,22)*three-r4(2,7)-r4(2,11)*two &
&                  -r4(2,13)-r4(9,31)*six-r4(9,35)*three+r3(5,28)*three+r3(5,32)*six &
&                  +r3(5,40)*nine-r2(4,22)*three-r2(4,30)*six-r2(4,34)*three)*qz+(+r6(20,10) &
&                  +r6(20,11)*two-r5(14,14)*two-r5(14,15)*four-r5(14,18)-r5(14,19)*two &
&                  +r4(9,15)+r4(9,16)*two+r4(9,19)*two+r4(9,20)*four+r4(9,27)*three &
&                  +r4(9,28)*six-r3(5,12)-r3(5,13)*two-r3(5,20)*two-r3(5,21)*four-r3(5,24) &
&                  -r3(5,25)*two)*xx+(+r6(14,11)*two+r6(14,12)-r5(9,15)*four-r5(9,16)*two &
&                  -r5(9,19)*two-r5(9,20)+r4(5,16)*two+r4(5,17)+r4(5,20)*four+r4(5,21)*two &
&                  +r4(5,28)*six+r4(5,29)*three+r4(14,40)*two+r4(14,41)-r3(2,13)*two-r3(2,14) &
&                  -r3(2,21)*four-r3(2,22)*two-r3(2,25)*two-r3(2,26)-r3(9,45)*four &
&                  -r3(9,46)*two-r3(9,50)*two-r3(9,51)+r2(6,39)*two+r2(6,40)+r2(6,44)*four &
&                  +r2(6,45)*two+r2(6,54)*six+r2(6,55)*three-r1(2,23)*two-r1(2,24) &
&                  -r1(2,33)*four-r1(2,34)*two-r1(2,38)*two-r1(2,39))*xz+(+r5(20,23) &
&                  -r4(14,32)*two-r4(14,36)+r3(9,29)+r3(9,33)*two+r3(9,41)*three-r2(6,23) &
&                  -r2(6,31)*two-r2(6,35))*xxx+(+r5(14,23)+r5(14,24)*two-r4(9,32)*two &
&                  -r4(9,33)*four-r4(9,36)-r4(9,37)*two+r3(5,29)+r3(5,30)*two+r3(5,33)*two &
&                  +r3(5,34)*four+r3(5,41)*three+r3(5,42)*six-r2(4,23)-r2(4,24)*two &
&                  -r2(4,31)*two-r2(4,32)*four-r2(4,35)-r2(4,36)*two)*xxz+rxyz(1)*xxxz
      eri(2,5,3,6)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,3,6)=r103+(+r7(35,4)-r6(27,6)*two-r6(27,8)+r5(20,6)+r5(20,8)*two &
&                  +r5(20,12)*three+r5(20,22)*three-r4(14,7)-r4(14,11)*two-r4(14,13) &
&                  -r4(14,31)*six-r4(14,35)*three+r3(9,28)*three+r3(9,32)*six+r3(9,40)*nine &
&                  -r2(6,22)*three-r2(6,30)*six-r2(6,34)*three)*qx+(+r7(27,3)*two+r7(27,4) &
&                  -r6(20,5)*four-r6(20,6)*two-r6(20,7)*two-r6(20,8)+r5(14,5)*two+r5(14,6) &
&                  +r5(14,7)*four+r5(14,8)*two+r5(14,11)*six+r5(14,12)*three+r5(14,21)*two &
&                  +r5(14,22)-r4(9,6)*two-r4(9,7)-r4(9,10)*four-r4(9,11)*two-r4(9,12)*two &
&                  -r4(9,13)-r4(9,30)*four-r4(9,31)*two-r4(9,34)*two-r4(9,35)+r3(5,27)*two &
&                  +r3(5,28)+r3(5,31)*four+r3(5,32)*two+r3(5,39)*six+r3(5,40)*three &
&                  -r2(4,21)*two-r2(4,22)-r2(4,29)*four-r2(4,30)*two-r2(4,33)*two-r2(4,34)) &
&                  *qz+(+r6(27,11)*two+r6(27,12)-r5(20,15)*four-r5(20,16)*two-r5(20,19)*two &
&                  -r5(20,20)+r4(14,16)*two+r4(14,17)+r4(14,20)*four+r4(14,21)*two &
&                  +r4(14,28)*six+r4(14,29)*three+r4(14,40)*two+r4(14,41)-r3(9,13)*two &
&                  -r3(9,14)-r3(9,21)*four-r3(9,22)*two-r3(9,25)*two-r3(9,26)-r3(9,45)*four &
&                  -r3(9,46)*two-r3(9,50)*two-r3(9,51)+r2(6,39)*two+r2(6,40)+r2(6,44)*four &
&                  +r2(6,45)*two+r2(6,54)*six+r2(6,55)*three-r1(2,23)*two-r1(2,24) &
&                  -r1(2,33)*four-r1(2,34)*two-r1(2,38)*two-r1(2,39))*xz+(+r6(20,10) &
&                  +r6(20,11)*two-r5(14,14)*two-r5(14,15)*four-r5(14,18)-r5(14,19)*two &
&                  +r4(9,15)+r4(9,16)*two+r4(9,19)*two+r4(9,20)*four+r4(9,27)*three &
&                  +r4(9,28)*six-r3(5,12)-r3(5,13)*two-r3(5,20)*two-r3(5,21)*four-r3(5,24) &
&                  -r3(5,25)*two)*zz+(+r5(20,23)+r5(20,24)*two-r4(14,32)*two-r4(14,33)*four &
&                  -r4(14,36)-r4(14,37)*two+r3(9,29)+r3(9,30)*two+r3(9,33)*two+r3(9,34)*four &
&                  +r3(9,41)*three+r3(9,42)*six-r2(6,23)-r2(6,24)*two-r2(6,31)*two &
&                  -r2(6,32)*four-r2(6,35)-r2(6,36)*two)*xzz+(+r5(14,23)-r4(9,32)*two &
&                  -r4(9,36)+r3(5,29)+r3(5,33)*two+r3(5,41)*three-r2(4,23)-r2(4,31)*two &
&                  -r2(4,35))*zzz+rxyz(1)*xzzz
      eri(4,5,3,6)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(19,11)+r6(19,12) &
&                  -r5(13,15)*two-r5(13,16)*two-r5(13,19)-r5(13,20)+r4(8,16)+r4(8,17) &
&                  +r4(8,20)*two+r4(8,21)*two+r4(8,28)*three+r4(8,29)*three-r3(4,13)-r3(4,14) &
&                  -r3(4,21)*two-r3(4,22)*two-r3(4,25)-r3(4,26))*xz+rxyz(2)*xxz
      eri(5,5,3,6)=r202+(+r7(27,3)+r7(27,4)-r6(20,5)*two-r6(20,6)*two-r6(20,7)-r6(20,8) &
&                  +r5(14,5)+r5(14,6)+r5(14,7)*two+r5(14,8)*two+r5(14,11)*three &
&                  +r5(14,12)*three+r5(14,21)+r5(14,22)-r4(9,6)-r4(9,7)-r4(9,10)*two &
&                  -r4(9,11)*two-r4(9,12)-r4(9,13)-r4(9,30)*two-r4(9,31)*two-r4(9,34) &
&                  -r4(9,35)+r3(5,27)+r3(5,28)+r3(5,31)*two+r3(5,32)*two+r3(5,39)*three &
&                  +r3(5,40)*three-r2(4,21)-r2(4,22)-r2(4,29)*two-r2(4,30)*two-r2(4,33) &
&                  -r2(4,34))*qx+(+r7(20,3)+r7(20,4)-r6(14,5)*two-r6(14,6)*two-r6(14,7) &
&                  -r6(14,8)+r5(9,5)+r5(9,6)+r5(9,7)*two+r5(9,8)*two+r5(9,11)*three &
&                  +r5(9,12)*three+r5(20,21)+r5(20,22)-r4(5,6)-r4(5,7)-r4(5,10)*two &
&                  -r4(5,11)*two-r4(5,12)-r4(5,13)-r4(14,30)*two-r4(14,31)*two-r4(14,34) &
&                  -r4(14,35)+r3(9,27)+r3(9,28)+r3(9,31)*two+r3(9,32)*two+r3(9,39)*three &
&                  +r3(9,40)*three-r2(6,21)-r2(6,22)-r2(6,29)*two-r2(6,30)*two-r2(6,33) &
&                  -r2(6,34))*qz+(+r6(27,11)-r5(20,15)*two-r5(20,19)+r4(14,16)+r4(14,20)*two &
&                  +r4(14,28)*three+r4(14,40)-r3(9,13)-r3(9,21)*two-r3(9,25)-r3(9,45)*two &
&                  -r3(9,50)+r2(6,39)+r2(6,44)*two+r2(6,54)*three-r1(2,23)-r1(2,33)*two &
&                  -r1(2,38))*xx+(+r6(20,10)+r6(20,11)*two+r6(20,12)-r5(14,14)*two &
&                  -r5(14,15)*four-r5(14,16)*two-r5(14,18)-r5(14,19)*two-r5(14,20)+r4(9,15) &
&                  +r4(9,16)*two+r4(9,17)+r4(9,19)*two+r4(9,20)*four+r4(9,21)*two &
&                  +r4(9,27)*three+r4(9,28)*six+r4(9,29)*three-r3(5,12)-r3(5,13)*two-r3(5,14) &
&                  -r3(5,20)*two-r3(5,21)*four-r3(5,22)*two-r3(5,24)-r3(5,25)*two-r3(5,26)) &
&                  *xz+(+r6(14,11)-r5(9,15)*two-r5(9,19)+r4(5,16)+r4(5,20)*two+r4(5,28)*three &
&                  +r4(14,40)-r3(2,13)-r3(2,21)*two-r3(2,25)-r3(9,45)*two-r3(9,50)+r2(6,39) &
&                  +r2(6,44)*two+r2(6,54)*three-r1(2,23)-r1(2,33)*two-r1(2,38))*zz+( &
&                  +r5(20,23)+r5(20,24)-r4(14,32)*two-r4(14,33)*two-r4(14,36)-r4(14,37) &
&                  +r3(9,29)+r3(9,30)+r3(9,33)*two+r3(9,34)*two+r3(9,41)*three+r3(9,42)*three &
&                  -r2(6,23)-r2(6,24)-r2(6,31)*two-r2(6,32)*two-r2(6,35)-r2(6,36))*xxz+( &
&                  +r5(14,23)+r5(14,24)-r4(9,32)*two-r4(9,33)*two-r4(9,36)-r4(9,37)+r3(5,29) &
&                  +r3(5,30)+r3(5,33)*two+r3(5,34)*two+r3(5,41)*three+r3(5,42)*three-r2(4,23) &
&                  -r2(4,24)-r2(4,31)*two-r2(4,32)*two-r2(4,35)-r2(4,36))*xzz+rxyz(1)*xxzz
      eri(6,5,3,6)=r112+rxyz(19)*qx+(+r7(26,3)+r7(26,4)-r6(19,5)*two-r6(19,6)*two &
&                  -r6(19,7)-r6(19,8)+r5(13,5)+r5(13,6)+r5(13,7)*two+r5(13,8)*two &
&                  +r5(13,11)*three+r5(13,12)*three-r4(8,6)-r4(8,7)-r4(8,10)*two-r4(8,11)*two &
&                  -r4(8,12)-r4(8,13))*qz+(+r6(26,11)+r6(26,12)-r5(19,15)*two-r5(19,16)*two &
&                  -r5(19,19)-r5(19,20)+r4(13,16)+r4(13,17)+r4(13,20)*two+r4(13,21)*two &
&                  +r4(13,28)*three+r4(13,29)*three-r3(8,13)-r3(8,14)-r3(8,21)*two &
&                  -r3(8,22)*two-r3(8,25)-r3(8,26))*xz+rxyz(11)*zz+rxyz(2)*xzz
      eri(1,6,3,6)=r211+(+r7(26,3)*two-r6(19,5)*four-r6(19,7)*two+r5(13,5)*two &
&                  +r5(13,7)*four+r5(13,11)*six-r4(8,6)*two-r4(8,10)*four-r4(8,12)*two)*qx &
&                  +rxyz(16)*qz+(+r6(26,10)-r5(19,14)*two-r5(19,18)+r4(13,15)+r4(13,19)*two &
&                  +r4(13,27)*three-r3(8,12)-r3(8,20)*two-r3(8,24))*xx+(+r6(19,11)*two &
&                  -r5(13,15)*four-r5(13,19)*two+r4(8,16)*two+r4(8,20)*four+r4(8,28)*six &
&                  -r3(4,13)*two-r3(4,21)*four-r3(4,25)*two)*xz+rxyz(3)*xxz
      eri(2,6,3,6)=r031+rxyz(8)*qz
      eri(3,6,3,6)=r013+(+r7(34,3)*two+r7(34,4)-r6(26,5)*four-r6(26,6)*two-r6(26,7)*two &
&                  -r6(26,8)+r5(19,5)*two+r5(19,6)+r5(19,7)*four+r5(19,8)*two+r5(19,11)*six &
&                  +r5(19,12)*three+r5(19,21)*two+r5(19,22)-r4(13,6)*two-r4(13,7) &
&                  -r4(13,10)*four-r4(13,11)*two-r4(13,12)*two-r4(13,13)-r4(13,30)*four &
&                  -r4(13,31)*two-r4(13,34)*two-r4(13,35)+r3(8,27)*two+r3(8,28)+r3(8,31)*four &
&                  +r3(8,32)*two+r3(8,39)*six+r3(8,40)*three-r2(2,21)*two-r2(2,22) &
&                  -r2(2,29)*four-r2(2,30)*two-r2(2,33)*two-r2(2,34))*qz+(+r6(26,10) &
&                  +r6(26,11)*two-r5(19,14)*two-r5(19,15)*four-r5(19,18)-r5(19,19)*two &
&                  +r4(13,15)+r4(13,16)*two+r4(13,19)*two+r4(13,20)*four+r4(13,27)*three &
&                  +r4(13,28)*six-r3(8,12)-r3(8,13)*two-r3(8,20)*two-r3(8,21)*four-r3(8,24) &
&                  -r3(8,25)*two)*zz+rxyz(3)*zzz
      eri(4,6,3,6)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,3,6)=r112+rxyz(18)*qx+(+r7(26,3)+r7(26,4)-r6(19,5)*two-r6(19,6)*two &
&                  -r6(19,7)-r6(19,8)+r5(13,5)+r5(13,6)+r5(13,7)*two+r5(13,8)*two &
&                  +r5(13,11)*three+r5(13,12)*three-r4(8,6)-r4(8,7)-r4(8,10)*two-r4(8,11)*two &
&                  -r4(8,12)-r4(8,13))*qz+(+r6(26,10)+r6(26,11)-r5(19,14)*two-r5(19,15)*two &
&                  -r5(19,18)-r5(19,19)+r4(13,15)+r4(13,16)+r4(13,19)*two+r4(13,20)*two &
&                  +r4(13,27)*three+r4(13,28)*three-r3(8,12)-r3(8,13)-r3(8,20)*two &
&                  -r3(8,21)*two-r3(8,24)-r3(8,25))*xz+rxyz(11)*zz+rxyz(3)*xzz
      eri(6,6,3,6)=r022+(+r7(33,3)+r7(33,4)-r6(25,5)*two-r6(25,6)*two-r6(25,7)-r6(25,8) &
&                  +r5(18,5)+r5(18,6)+r5(18,7)*two+r5(18,8)*two+r5(18,11)*three &
&                  +r5(18,12)*three+r5(20,21)+r5(20,22)-r4(12,6)-r4(12,7)-r4(12,10)*two &
&                  -r4(12,11)*two-r4(12,12)-r4(12,13)-r4(14,30)*two-r4(14,31)*two-r4(14,34) &
&                  -r4(14,35)+r3(9,27)+r3(9,28)+r3(9,31)*two+r3(9,32)*two+r3(9,39)*three &
&                  +r3(9,40)*three-r2(6,21)-r2(6,22)-r2(6,29)*two-r2(6,30)*two-r2(6,33) &
&                  -r2(6,34))*qz+rxyz(6)*zz
!
      do i= 1,6
        do j= 1,6
          eri(j,i,4,6)= eri(j,i,2,5)
        enddo
      enddo
      do i= 1,6
        do j= 1,6
          eri(j,i,5,6)= eri(j,i,6,5)
        enddo
      enddo
!
      r400= r8(13)-r7(8,1)-r7(8,2)+r6(4,2)+r6(4,4)+r6(6,4)+r6(13,9)*six-r5(3,3)-r5(3,4) &
&          -r5(8,13)*six-r5(8,17)*six+r4(1,3)+r4(1,5)+r4(4,18)*six+r4(4,26)*six &
&          +r4(6,26)*six+r4(13,38)*three-r3(3,19)*six-r3(3,23)*six-r3(8,43)*three &
&          -r3(8,48)*three+r2(1,9)*six+r2(1,17)*six+r2(2,42)*three+r2(2,52)*three &
&          +r2(3,52)*three-r1(3,31)*three-r1(3,36)*three+r0(11)*three+r0(21)*three
      r310= r8(18)-r7(12,1)-r7(12,2)+r6(7,2)+r6(7,4)+r6(9,4)+r6(18,9)*three-r5(5,3)-r5(5,4) &
&          -r5(12,13)*three-r5(12,17)*three+r4(2,3)+r4(2,5)+r4(7,18)*three+r4(7,26)*three &
&          +r4(9,26)*three-r3(5,19)*three-r3(5,23)*three+r2(4,9)*three+r2(4,17)*three
      r301= r8(19)-r7(13,1)-r7(13,2)+r6(8,2)+r6(8,4)+r6(10,4)+r6(19,9)*three-r5(6,3) &
&          -r5(6,4)-r5(13,13)*three-r5(13,17)*three+r4(3,3)+r4(3,5)+r4(8,18)*three &
&          +r4(8,26)*three+r4(10,26)*three-r3(6,19)*three-r3(6,23)*three+r2(5,9)*three &
&          +r2(5,17)*three
      r220= r8(24)-r7(17,1)-r7(17,2)+r6(11,2)+r6(11,4)+r6(13,4)+r6(13,9)+r6(24,9)-r5(8,3) &
&          -r5(8,4)-r5(8,13)-r5(17,13)-r5(8,17)-r5(17,17)+r4(4,3)+r4(4,5)+r4(4,18) &
&          +r4(11,18)+r4(4,26)+r4(6,26)+r4(11,26)+r4(13,26)+r4(13,38)-r3(3,19)-r3(8,19) &
&          -r3(3,23)-r3(8,23)-r3(8,43)-r3(8,48)+r2(1,9)+r2(2,9)+r2(1,17)+r2(2,17)+r2(2,42) &
&          +r2(2,52)+r2(3,52)-r1(3,31)-r1(3,36)+r0(11)+r0(21)
      r211= r8(25)-r7(18,1)-r7(18,2)+r6(12,2)+r6(12,4)+r6(14,4)+r6(25,9)-r5(9,3)-r5(9,4) &
&          -r5(18,13)-r5(18,17)+r4(5,3)+r4(5,5)+r4(12,18)+r4(12,26)+r4(14,26)-r3(9,19) &
&          -r3(9,23)+r2(6,9)+r2(6,17)
      r202= r8(26)-r7(19,1)-r7(19,2)+r6(13,2)+r6(13,4)+r6(15,4)+r6(13,9)+r6(26,9)-r5(10,3) &
&          -r5(10,4)-r5(8,13)-r5(19,13)-r5(8,17)-r5(19,17)+r4(6,3)+r4(6,5)+r4(4,18) &
&          +r4(13,18)+r4(4,26)+r4(6,26)+r4(13,26)+r4(15,26)+r4(13,38)-r3(3,19)-r3(10,19) &
&          -r3(3,23)-r3(10,23)-r3(8,43)-r3(8,48)+r2(1,9)+r2(3,9)+r2(1,17)+r2(3,17) &
&          +r2(2,42)+r2(2,52)+r2(3,52)-r1(3,31)-r1(3,36)+r0(11)+r0(21)
      r130= r8(31)-r7(23,1)-r7(23,2)+r6(16,2)+r6(16,4)+r6(18,4)+r6(18,9)*three-r5(12,3) &
&          -r5(12,4)-r5(12,13)*three-r5(12,17)*three+r4(7,3)+r4(7,5)+r4(7,18)*three &
&          +r4(7,26)*three+r4(9,26)*three-r3(5,19)*three-r3(5,23)*three+r2(4,9)*three &
&          +r2(4,17)*three
      r121= r8(32)-r7(24,1)-r7(24,2)+r6(17,2)+r6(17,4)+r6(19,4)+r6(19,9)-r5(13,3)-r5(13,4) &
&          -r5(13,13)-r5(13,17)+r4(8,3)+r4(8,5)+r4(8,18)+r4(8,26)+r4(10,26)-r3(6,19) &
&          -r3(6,23)+r2(5,9)+r2(5,17)
      r112= r8(33)-r7(25,1)-r7(25,2)+r6(18,2)+r6(18,4)+r6(20,4)+r6(18,9)-r5(14,3)-r5(14,4) &
&          -r5(12,13)-r5(12,17)+r4(9,3)+r4(9,5)+r4(7,18)+r4(7,26)+r4(9,26)-r3(5,19) &
&          -r3(5,23)+r2(4,9)+r2(4,17)
      r103= r8(34)-r7(26,1)-r7(26,2)+r6(19,2)+r6(19,4)+r6(21,4)+r6(19,9)*three-r5(15,3) &
&          -r5(15,4)-r5(13,13)*three-r5(13,17)*three+r4(10,3)+r4(10,5)+r4(8,18)*three &
&          +r4(8,26)*three+r4(10,26)*three-r3(6,19)*three-r3(6,23)*three+r2(5,9)*three &
&          +r2(5,17)*three
      r040= r8(39)-r7(30,1)-r7(30,2)+r6(22,2)+r6(22,4)+r6(24,4)+r6(24,9)*six-r5(17,3) &
&          -r5(17,4)-r5(17,13)*six-r5(17,17)*six+r4(11,3)+r4(11,5)+r4(11,18)*six &
&          +r4(11,26)*six+r4(13,26)*six+r4(13,38)*three-r3(8,19)*six-r3(8,23)*six &
&          -r3(8,43)*three-r3(8,48)*three+r2(2,9)*six+r2(2,17)*six+r2(2,42)*three &
&          +r2(2,52)*three+r2(3,52)*three-r1(3,31)*three-r1(3,36)*three+r0(11)*three &
&          +r0(21)*three
      r031= r8(40)-r7(31,1)-r7(31,2)+r6(23,2)+r6(23,4)+r6(25,4)+r6(25,9)*three-r5(18,3) &
&          -r5(18,4)-r5(18,13)*three-r5(18,17)*three+r4(12,3)+r4(12,5)+r4(12,18)*three &
&          +r4(12,26)*three+r4(14,26)*three-r3(9,19)*three-r3(9,23)*three+r2(6,9)*three &
&          +r2(6,17)*three
      r022= r8(41)-r7(32,1)-r7(32,2)+r6(24,2)+r6(24,4)+r6(26,4)+r6(24,9)+r6(26,9)-r5(19,3) &
&          -r5(19,4)-r5(17,13)-r5(19,13)-r5(17,17)-r5(19,17)+r4(13,3)+r4(13,5)+r4(11,18) &
&          +r4(13,18)+r4(11,26)+r4(13,26)*two+r4(15,26)+r4(13,38)-r3(8,19)-r3(10,19) &
&          -r3(8,23)-r3(10,23)-r3(8,43)-r3(8,48)+r2(2,9)+r2(3,9)+r2(2,17)+r2(3,17) &
&          +r2(2,42)+r2(2,52)+r2(3,52)-r1(3,31)-r1(3,36)+r0(11)+r0(21)
      r013= r8(42)-r7(33,1)-r7(33,2)+r6(25,2)+r6(25,4)+r6(27,4)+r6(25,9)*three-r5(20,3) &
&          -r5(20,4)-r5(18,13)*three-r5(18,17)*three+r4(14,3)+r4(14,5)+r4(12,18)*three &
&          +r4(12,26)*three+r4(14,26)*three-r3(9,19)*three-r3(9,23)*three+r2(6,9)*three &
&          +r2(6,17)*three
      r004= r8(43)-r7(34,1)-r7(34,2)+r6(26,2)+r6(26,4)+r6(28,4)+r6(26,9)*six-r5(21,3) &
&          -r5(21,4)-r5(19,13)*six-r5(19,17)*six+r4(15,3)+r4(15,5)+r4(13,18)*six &
&          +r4(13,26)*six+r4(15,26)*six+r4(13,38)*three-r3(10,19)*six-r3(10,23)*six &
&          -r3(8,43)*three-r3(8,48)*three+r2(3,9)*six+r2(3,17)*six+r2(2,42)*three &
&          +r2(2,52)*three+r2(3,52)*three-r1(3,31)*three-r1(3,36)*three+r0(11)*three &
&          +r0(21)*three
      rxyz(1)=+r4(13,42)-r3(8,47)-r3(8,52)+r2(2,46)+r2(2,56)+r2(3,56)-r1(3,35)-r1(3,40) &
&             +r0(15)+r0(25)
      rxyz(2)=+r5(18,24)-r4(12,33)-r4(12,37)+r3(7,34)+r3(7,42)+r3(9,42)-r2(6,32)-r2(6,36) &
&             +r1(2,12)+r1(2,20)
      rxyz(3)=+r5(18,23)-r4(12,32)-r4(12,36)+r3(7,33)+r3(7,41)+r3(9,41)-r2(6,31)-r2(6,35) &
&             +r1(2,11)+r1(2,19)
      rxyz(4)=+r6(24,10)-r5(17,14)-r5(17,18)+r4(11,19)+r4(11,27)+r4(13,27)+r4(13,39) &
&             -r3(8,20)-r3(8,24)-r3(8,44)-r3(8,49)+r2(2,10)+r2(2,18)+r2(2,43)+r2(2,53) &
&             +r2(3,53)-r1(3,32)-r1(3,37)+r0(12)+r0(22)
      rxyz(5)=+r6(24,12)-r5(17,16)-r5(17,20)+r4(11,21)+r4(11,29)+r4(13,29)+r4(13,41) &
&             -r3(8,22)-r3(8,26)-r3(8,46)-r3(8,51)+r2(2,12)+r2(2,20)+r2(2,45)+r2(2,55) &
&             +r2(3,55)-r1(3,34)-r1(3,39)+r0(14)+r0(24)
      rxyz(6)=+r6(24,11)-r5(17,15)-r5(17,19)+r4(11,20)+r4(11,28)+r4(13,28)+r4(13,40) &
&             -r3(8,21)-r3(8,25)-r3(8,45)-r3(8,50)+r2(2,11)+r2(2,19)+r2(2,44)+r2(2,54) &
&             +r2(3,54)-r1(3,33)-r1(3,38)+r0(13)+r0(23)
      rxyz(7)=+r7(31,3)-r6(23,5)-r6(23,7)+r5(16,7)+r5(16,11)+r5(18,11)+r5(18,21)*three &
&             -r4(12,10)-r4(12,12)-r4(12,30)*three-r4(12,34)*three+r3(7,5)+r3(7,9) &
&             +r3(7,31)*three+r3(7,39)*three+r3(9,39)*three-r2(6,29)*three-r2(6,33)*three &
&             +r1(2,9)*three+r1(2,17)*three
      rxyz(8)=+r7(31,4)-r6(23,6)-r6(23,8)+r5(16,8)+r5(16,12)+r5(18,12)+r5(18,22)*three &
&             -r4(12,11)-r4(12,13)-r4(12,31)*three-r4(12,35)*three+r3(7,6)+r3(7,10) &
&             +r3(7,32)*three+r3(7,40)*three+r3(9,40)*three-r2(6,30)*three-r2(6,34)*three &
&             +r1(2,10)*three+r1(2,18)*three
      rxyz(9)=+r7(25,3)+r7(25,4)-r6(18,5)-r6(18,6)-r6(18,7)-r6(18,8)+r5(12,7)+r5(12,8) &
&             +r5(12,11)+r5(14,11)+r5(12,12)+r5(14,12)-r4(9,10)-r4(9,11)-r4(9,12)-r4(9,13) &
&             +r3(5,5)+r3(5,6)+r3(5,9)+r3(5,10)
      rxyz(10)=+r6(25,11)-r5(18,15)-r5(18,19)+r4(12,20)+r4(12,28)+r4(14,28)-r3(9,21) &
&             -r3(9,25)+r2(6,11)+r2(6,19)
      rxyz(11)=+r6(18,11)-r5(12,15)-r5(12,19)+r4(7,20)+r4(7,28)+r4(9,28)-r3(5,21)-r3(5,25) &
&             +r2(4,11)+r2(4,19)
      rxyz(12)=+r7(32,3)-r6(24,5)-r6(24,7)+r5(17,7)+r5(17,11)+r5(19,11)+r5(19,21) &
&             -r4(13,10)-r4(13,12)-r4(13,30)-r4(13,34)+r3(8,5)+r3(8,9)+r3(8,31)+r3(8,39) &
&             +r3(10,39)-r2(3,29)-r2(3,33)+r1(3,9)+r1(3,17)
      rxyz(13)=+r7(32,4)-r6(24,6)-r6(24,8)+r5(17,8)+r5(17,12)+r5(19,12)+r5(19,22) &
&             -r4(13,11)-r4(13,13)-r4(13,31)-r4(13,35)+r3(8,6)+r3(8,10)+r3(8,32)+r3(8,40) &
&             +r3(10,40)-r2(3,30)-r2(3,34)+r1(3,10)+r1(3,18)
      rxyz(14)=+r7(24,4)-r6(17,6)-r6(17,8)+r5(11,8)+r5(11,12)+r5(13,12)+r5(13,22)-r4(8,11) &
&             -r4(8,13)-r4(8,31)-r4(8,35)+r3(4,6)+r3(4,10)+r3(4,32)+r3(4,40)+r3(6,40) &
&             -r2(5,30)-r2(5,34)+r1(1,10)+r1(1,18)
      rxyz(15)=+r7(24,3)-r6(17,5)-r6(17,7)+r5(11,7)+r5(11,11)+r5(13,11)+r5(13,21)-r4(8,10) &
&             -r4(8,12)-r4(8,30)-r4(8,34)+r3(4,5)+r3(4,9)+r3(4,31)+r3(4,39)+r3(6,39) &
&             -r2(5,29)-r2(5,33)+r1(1,9)+r1(1,17)
      rxyz(16)=+r7(18,4)-r6(12,6)-r6(12,8)+r5(7,8)+r5(7,12)+r5(9,12)+r5(18,22)-r4(5,11) &
&             -r4(5,13)-r4(12,31)-r4(12,35)+r3(2,6)+r3(2,10)+r3(7,32)+r3(7,40)+r3(9,40) &
&             -r2(6,30)-r2(6,34)+r1(2,10)+r1(2,18)
      rxyz(17)=+r7(18,3)-r6(12,5)-r6(12,7)+r5(7,7)+r5(7,11)+r5(9,11)+r5(18,21)-r4(5,10) &
&             -r4(5,12)-r4(12,30)-r4(12,34)+r3(2,5)+r3(2,9)+r3(7,31)+r3(7,39)+r3(9,39) &
&             -r2(6,29)-r2(6,33)+r1(2,9)+r1(2,17)
      rxyz(18)=+r7(33,3)-r6(25,5)-r6(25,7)+r5(18,7)+r5(18,11)+r5(20,11)+r5(18,21) &
&             -r4(14,10)-r4(14,12)-r4(12,30)-r4(12,34)+r3(9,5)+r3(9,9)+r3(7,31)+r3(7,39) &
&             +r3(9,39)-r2(6,29)-r2(6,33)+r1(2,9)+r1(2,17)
      rxyz(19)=+r7(33,4)-r6(25,6)-r6(25,8)+r5(18,8)+r5(18,12)+r5(20,12)+r5(18,22) &
&             -r4(14,11)-r4(14,13)-r4(12,31)-r4(12,35)+r3(9,6)+r3(9,10)+r3(7,32)+r3(7,40) &
&             +r3(9,40)-r2(6,30)-r2(6,34)+r1(2,10)+r1(2,18)
      rxyz(20)=+r6(19,11)*four-r5(13,15)*four-r5(13,19)*four+r4(8,20)*four+r4(8,28)*four &
&             +r4(10,28)*four-r3(6,21)*four-r3(6,25)*four+r2(5,11)*four+r2(5,19)*four
      eri(1,1,6,6)=r400+(+r7(13,3)*two+r7(13,4)*two-r6(8,5)*two-r6(8,6)*two-r6(8,7)*two &
&                  -r6(8,8)*two+r5(4,7)*two+r5(4,8)*two+r5(4,11)*two+r5(6,11)*two &
&                  +r5(4,12)*two+r5(6,12)*two+r5(13,21)*six+r5(13,22)*six-r4(3,10)*two &
&                  -r4(3,11)*two-r4(3,12)*two-r4(3,13)*two-r4(8,30)*six-r4(8,31)*six &
&                  -r4(8,34)*six-r4(8,35)*six+r3(1,5)*two+r3(1,6)*two+r3(1,9)*two &
&                  +r3(1,10)*two+r3(4,31)*six+r3(4,32)*six+r3(4,39)*six+r3(6,39)*six &
&                  +r3(4,40)*six+r3(6,40)*six-r2(5,29)*six-r2(5,30)*six-r2(5,33)*six &
&                  -r2(5,34)*six+r1(1,9)*six+r1(1,10)*six+r1(1,17)*six+r1(1,18)*six)*qx+( &
&                  +r6(13,10)+r6(13,11)*four+r6(13,12)-r5(8,14)-r5(8,15)*four-r5(8,16) &
&                  -r5(8,18)-r5(8,19)*four-r5(8,20)+r4(4,19)+r4(4,20)*four+r4(4,21)+r4(4,27) &
&                  +r4(6,27)+r4(4,28)*four+r4(6,28)*four+r4(4,29)+r4(6,29)+r4(13,39) &
&                  +r4(13,40)*four+r4(13,41)-r3(3,20)-r3(3,21)*four-r3(3,22)-r3(3,24) &
&                  -r3(3,25)*four-r3(3,26)-r3(8,44)-r3(8,45)*four-r3(8,46)-r3(8,49) &
&                  -r3(8,50)*four-r3(8,51)+r2(1,10)+r2(1,11)*four+r2(1,12)+r2(1,18) &
&                  +r2(1,19)*four+r2(1,20)+r2(2,43)+r2(2,44)*four+r2(2,45)+r2(2,53)+r2(3,53) &
&                  +r2(2,54)*four+r2(3,54)*four+r2(2,55)+r2(3,55)-r1(3,32)-r1(3,33)*four &
&                  -r1(3,34)-r1(3,37)-r1(3,38)*four-r1(3,39)+r0(12)+r0(13)*four+r0(14)+r0(22) &
&                  +r0(23)*four+r0(24))*xx+(+r5(13,23)*two+r5(13,24)*two-r4(8,32)*two &
&                  -r4(8,33)*two-r4(8,36)*two-r4(8,37)*two+r3(4,33)*two+r3(4,34)*two &
&                  +r3(4,41)*two+r3(6,41)*two+r3(4,42)*two+r3(6,42)*two-r2(5,31)*two &
&                  -r2(5,32)*two-r2(5,35)*two-r2(5,36)*two+r1(1,11)*two+r1(1,12)*two &
&                  +r1(1,19)*two+r1(1,20)*two)*xxx+rxyz(1)*xxxx
      eri(2,1,6,6)=r220+(+r7(24,4)*two-r6(17,6)*two-r6(17,8)*two+r5(11,8)*two &
&                  +r5(11,12)*two+r5(13,12)*two+r5(13,22)*two-r4(8,11)*two-r4(8,13)*two &
&                  -r4(8,31)*two-r4(8,35)*two+r3(4,6)*two+r3(4,10)*two+r3(4,32)*two &
&                  +r3(4,40)*two+r3(6,40)*two-r2(5,30)*two-r2(5,34)*two+r1(1,10)*two &
&                  +r1(1,18)*two)*qx+rxyz(5)*xx
      eri(3,1,6,6)=r202+(+r7(26,4)*two-r6(19,6)*two-r6(19,8)*two+r5(13,8)*two &
&                  +r5(13,12)*two+r5(15,12)*two+r5(13,22)*two-r4(10,11)*two-r4(10,13)*two &
&                  -r4(8,31)*two-r4(8,35)*two+r3(6,6)*two+r3(6,10)*two+r3(4,32)*two &
&                  +r3(4,40)*two+r3(6,40)*two-r2(5,30)*two-r2(5,34)*two+r1(1,10)*two &
&                  +r1(1,18)*two)*qx+(+r7(19,3)*two-r6(13,5)*two-r6(13,7)*two+r5(8,7)*two &
&                  +r5(8,11)*two+r5(10,11)*two+r5(19,21)*two-r4(6,10)*two-r4(6,12)*two &
&                  -r4(13,30)*two-r4(13,34)*two+r3(3,5)*two+r3(3,9)*two+r3(8,31)*two &
&                  +r3(8,39)*two+r3(10,39)*two-r2(3,29)*two-r2(3,33)*two+r1(3,9)*two &
&                  +r1(3,17)*two)*qz+(+r6(26,12)-r5(19,16)-r5(19,20)+r4(13,21)+r4(13,29) &
&                  +r4(15,29)+r4(13,41)-r3(10,22)-r3(10,26)-r3(8,46)-r3(8,51)+r2(3,12) &
&                  +r2(3,20)+r2(2,45)+r2(2,55)+r2(3,55)-r1(3,34)-r1(3,39)+r0(14)+r0(24))*xx &
&                  +rxyz(20)*xz+(+r6(13,10)-r5(8,14)-r5(8,18)+r4(4,19)+r4(4,27)+r4(6,27) &
&                  +r4(13,39)-r3(3,20)-r3(3,24)-r3(8,44)-r3(8,49)+r2(1,10)+r2(1,18)+r2(2,43) &
&                  +r2(2,53)+r2(3,53)-r1(3,32)-r1(3,37)+r0(12)+r0(22))*zz+(+r5(19,24)*two &
&                  -r4(13,33)*two-r4(13,37)*two+r3(8,34)*two+r3(8,42)*two+r3(10,42)*two &
&                  -r2(3,32)*two-r2(3,36)*two+r1(3,12)*two+r1(3,20)*two)*xxz+(+r5(13,23)*two &
&                  -r4(8,32)*two-r4(8,36)*two+r3(4,33)*two+r3(4,41)*two+r3(6,41)*two &
&                  -r2(5,31)*two-r2(5,35)*two+r1(1,11)*two+r1(1,19)*two)*xzz+rxyz(1)*xxzz
      eri(4,1,6,6)=r310+(+r7(18,3)+r7(18,4)*two-r6(12,5)-r6(12,6)*two-r6(12,7) &
&                  -r6(12,8)*two+r5(7,7)+r5(7,8)*two+r5(7,11)+r5(9,11)+r5(7,12)*two &
&                  +r5(9,12)*two+r5(18,21)+r5(18,22)*two-r4(5,10)-r4(5,11)*two-r4(5,12) &
&                  -r4(5,13)*two-r4(12,30)-r4(12,31)*two-r4(12,34)-r4(12,35)*two+r3(2,5) &
&                  +r3(2,6)*two+r3(2,9)+r3(2,10)*two+r3(7,31)+r3(7,32)*two+r3(7,39)+r3(9,39) &
&                  +r3(7,40)*two+r3(9,40)*two-r2(6,29)-r2(6,30)*two-r2(6,33)-r2(6,34)*two &
&                  +r1(2,9)+r1(2,10)*two+r1(2,17)+r1(2,18)*two)*qx+(+r6(18,11)*two+r6(18,12) &
&                  -r5(12,15)*two-r5(12,16)-r5(12,19)*two-r5(12,20)+r4(7,20)*two+r4(7,21) &
&                  +r4(7,28)*two+r4(9,28)*two+r4(7,29)+r4(9,29)-r3(5,21)*two-r3(5,22) &
&                  -r3(5,25)*two-r3(5,26)+r2(4,11)*two+r2(4,12)+r2(4,19)*two+r2(4,20))*xx &
&                  +rxyz(2)*xxx
      eri(5,1,6,6)=r301+(+r7(19,3)+r7(19,4)*two-r6(13,5)-r6(13,6)*two-r6(13,7) &
&                  -r6(13,8)*two+r5(8,7)+r5(8,8)*two+r5(8,11)+r5(10,11)+r5(8,12)*two &
&                  +r5(10,12)*two+r5(19,21)+r5(19,22)*two-r4(6,10)-r4(6,11)*two-r4(6,12) &
&                  -r4(6,13)*two-r4(13,30)-r4(13,31)*two-r4(13,34)-r4(13,35)*two+r3(3,5) &
&                  +r3(3,6)*two+r3(3,9)+r3(3,10)*two+r3(8,31)+r3(8,32)*two+r3(8,39)+r3(10,39) &
&                  +r3(8,40)*two+r3(10,40)*two-r2(3,29)-r2(3,30)*two-r2(3,33)-r2(3,34)*two &
&                  +r1(3,9)+r1(3,10)*two+r1(3,17)+r1(3,18)*two)*qx+(+r7(13,3)-r6(8,5)-r6(8,7) &
&                  +r5(4,7)+r5(4,11)+r5(6,11)+r5(13,21)*three-r4(3,10)-r4(3,12) &
&                  -r4(8,30)*three-r4(8,34)*three+r3(1,5)+r3(1,9)+r3(4,31)*three &
&                  +r3(4,39)*three+r3(6,39)*three-r2(5,29)*three-r2(5,33)*three+r1(1,9)*three &
&                  +r1(1,17)*three)*qz+(+r6(19,11)*two+r6(19,12)-r5(13,15)*two-r5(13,16) &
&                  -r5(13,19)*two-r5(13,20)+r4(8,20)*two+r4(8,21)+r4(8,28)*two+r4(10,28)*two &
&                  +r4(8,29)+r4(10,29)-r3(6,21)*two-r3(6,22)-r3(6,25)*two-r3(6,26) &
&                  +r2(5,11)*two+r2(5,12)+r2(5,19)*two+r2(5,20))*xx+(+r6(13,10)+r6(13,11)*two &
&                  -r5(8,14)-r5(8,15)*two-r5(8,18)-r5(8,19)*two+r4(4,19)+r4(4,20)*two &
&                  +r4(4,27)+r4(6,27)+r4(4,28)*two+r4(6,28)*two+r4(13,39)+r4(13,40)*two &
&                  -r3(3,20)-r3(3,21)*two-r3(3,24)-r3(3,25)*two-r3(8,44)-r3(8,45)*two &
&                  -r3(8,49)-r3(8,50)*two+r2(1,10)+r2(1,11)*two+r2(1,18)+r2(1,19)*two &
&                  +r2(2,43)+r2(2,44)*two+r2(2,53)+r2(3,53)+r2(2,54)*two+r2(3,54)*two &
&                  -r1(3,32)-r1(3,33)*two-r1(3,37)-r1(3,38)*two+r0(12)+r0(13)*two+r0(22) &
&                  +r0(23)*two)*xz+(+r5(19,24)-r4(13,33)-r4(13,37)+r3(8,34)+r3(8,42) &
&                  +r3(10,42)-r2(3,32)-r2(3,36)+r1(3,12)+r1(3,20))*xxx+(+r5(13,23)*two &
&                  +r5(13,24)-r4(8,32)*two-r4(8,33)-r4(8,36)*two-r4(8,37)+r3(4,33)*two &
&                  +r3(4,34)+r3(4,41)*two+r3(6,41)*two+r3(4,42)+r3(6,42)-r2(5,31)*two &
&                  -r2(5,32)-r2(5,35)*two-r2(5,36)+r1(1,11)*two+r1(1,12)+r1(1,19)*two &
&                  +r1(1,20))*xxz+rxyz(1)*xxxz
      eri(6,1,6,6)=r211+(+r7(25,4)*two-r6(18,6)*two-r6(18,8)*two+r5(12,8)*two &
&                  +r5(12,12)*two+r5(14,12)*two-r4(9,11)*two-r4(9,13)*two+r3(5,6)*two &
&                  +r3(5,10)*two)*qx+rxyz(17)*qz+(+r6(25,12)-r5(18,16)-r5(18,20)+r4(12,21) &
&                  +r4(12,29)+r4(14,29)-r3(9,22)-r3(9,26)+r2(6,12)+r2(6,20))*xx+( &
&                  +r6(18,11)*two-r5(12,15)*two-r5(12,19)*two+r4(7,20)*two+r4(7,28)*two &
&                  +r4(9,28)*two-r3(5,21)*two-r3(5,25)*two+r2(4,11)*two+r2(4,19)*two)*xz &
&                  +rxyz(2)*xxz
      eri(1,2,6,6)=r220+(+r7(24,3)*two-r6(17,5)*two-r6(17,7)*two+r5(11,7)*two &
&                  +r5(11,11)*two+r5(13,11)*two+r5(13,21)*two-r4(8,10)*two-r4(8,12)*two &
&                  -r4(8,30)*two-r4(8,34)*two+r3(4,5)*two+r3(4,9)*two+r3(4,31)*two &
&                  +r3(4,39)*two+r3(6,39)*two-r2(5,29)*two-r2(5,33)*two+r1(1,9)*two &
&                  +r1(1,17)*two)*qx+rxyz(4)*xx
      eri(2,2,6,6)=r040
      eri(3,2,6,6)=r022+(+r7(32,3)*two-r6(24,5)*two-r6(24,7)*two+r5(17,7)*two &
&                  +r5(17,11)*two+r5(19,11)*two+r5(19,21)*two-r4(13,10)*two-r4(13,12)*two &
&                  -r4(13,30)*two-r4(13,34)*two+r3(8,5)*two+r3(8,9)*two+r3(8,31)*two &
&                  +r3(8,39)*two+r3(10,39)*two-r2(3,29)*two-r2(3,33)*two+r1(3,9)*two &
&                  +r1(3,17)*two)*qz+rxyz(4)*zz
      eri(4,2,6,6)=r130+rxyz(7)*qx
      eri(5,2,6,6)=r121+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,6,6)=r031+rxyz(7)*qz
      eri(1,3,6,6)=r202+(+r7(26,3)*two-r6(19,5)*two-r6(19,7)*two+r5(13,7)*two &
&                  +r5(13,11)*two+r5(15,11)*two+r5(13,21)*two-r4(10,10)*two-r4(10,12)*two &
&                  -r4(8,30)*two-r4(8,34)*two+r3(6,5)*two+r3(6,9)*two+r3(4,31)*two &
&                  +r3(4,39)*two+r3(6,39)*two-r2(5,29)*two-r2(5,33)*two+r1(1,9)*two &
&                  +r1(1,17)*two)*qx+(+r7(19,4)*two-r6(13,6)*two-r6(13,8)*two+r5(8,8)*two &
&                  +r5(8,12)*two+r5(10,12)*two+r5(19,22)*two-r4(6,11)*two-r4(6,13)*two &
&                  -r4(13,31)*two-r4(13,35)*two+r3(3,6)*two+r3(3,10)*two+r3(8,32)*two &
&                  +r3(8,40)*two+r3(10,40)*two-r2(3,30)*two-r2(3,34)*two+r1(3,10)*two &
&                  +r1(3,18)*two)*qz+(+r6(26,10)-r5(19,14)-r5(19,18)+r4(13,19)+r4(13,27) &
&                  +r4(15,27)+r4(13,39)-r3(10,20)-r3(10,24)-r3(8,44)-r3(8,49)+r2(3,10) &
&                  +r2(3,18)+r2(2,43)+r2(2,53)+r2(3,53)-r1(3,32)-r1(3,37)+r0(12)+r0(22))*xx &
&                  +rxyz(20)*xz+(+r6(13,12)-r5(8,16)-r5(8,20)+r4(4,21)+r4(4,29)+r4(6,29) &
&                  +r4(13,41)-r3(3,22)-r3(3,26)-r3(8,46)-r3(8,51)+r2(1,12)+r2(1,20)+r2(2,45) &
&                  +r2(2,55)+r2(3,55)-r1(3,34)-r1(3,39)+r0(14)+r0(24))*zz+(+r5(19,23)*two &
&                  -r4(13,32)*two-r4(13,36)*two+r3(8,33)*two+r3(8,41)*two+r3(10,41)*two &
&                  -r2(3,31)*two-r2(3,35)*two+r1(3,11)*two+r1(3,19)*two)*xxz+(+r5(13,24)*two &
&                  -r4(8,33)*two-r4(8,37)*two+r3(4,34)*two+r3(4,42)*two+r3(6,42)*two &
&                  -r2(5,32)*two-r2(5,36)*two+r1(1,12)*two+r1(1,20)*two)*xzz+rxyz(1)*xxzz
      eri(2,3,6,6)=r022+(+r7(32,4)*two-r6(24,6)*two-r6(24,8)*two+r5(17,8)*two &
&                  +r5(17,12)*two+r5(19,12)*two+r5(19,22)*two-r4(13,11)*two-r4(13,13)*two &
&                  -r4(13,31)*two-r4(13,35)*two+r3(8,6)*two+r3(8,10)*two+r3(8,32)*two &
&                  +r3(8,40)*two+r3(10,40)*two-r2(3,30)*two-r2(3,34)*two+r1(3,10)*two &
&                  +r1(3,18)*two)*qz+rxyz(5)*zz
      eri(3,3,6,6)=r004+(+r7(34,3)*two+r7(34,4)*two-r6(26,5)*two-r6(26,6)*two &
&                  -r6(26,7)*two-r6(26,8)*two+r5(19,7)*two+r5(19,8)*two+r5(19,11)*two &
&                  +r5(21,11)*two+r5(19,12)*two+r5(21,12)*two+r5(19,21)*six+r5(19,22)*six &
&                  -r4(15,10)*two-r4(15,11)*two-r4(15,12)*two-r4(15,13)*two-r4(13,30)*six &
&                  -r4(13,31)*six-r4(13,34)*six-r4(13,35)*six+r3(10,5)*two+r3(10,6)*two &
&                  +r3(10,9)*two+r3(10,10)*two+r3(8,31)*six+r3(8,32)*six+r3(8,39)*six &
&                  +r3(10,39)*six+r3(8,40)*six+r3(10,40)*six-r2(3,29)*six-r2(3,30)*six &
&                  -r2(3,33)*six-r2(3,34)*six+r1(3,9)*six+r1(3,10)*six+r1(3,17)*six &
&                  +r1(3,18)*six)*qz+(+r6(26,10)+r6(26,11)*four+r6(26,12)-r5(19,14) &
&                  -r5(19,15)*four-r5(19,16)-r5(19,18)-r5(19,19)*four-r5(19,20)+r4(13,19) &
&                  +r4(13,20)*four+r4(13,21)+r4(13,27)+r4(15,27)+r4(13,28)*four &
&                  +r4(15,28)*four+r4(13,29)+r4(15,29)+r4(13,39)+r4(13,40)*four+r4(13,41) &
&                  -r3(10,20)-r3(10,21)*four-r3(10,22)-r3(10,24)-r3(10,25)*four-r3(10,26) &
&                  -r3(8,44)-r3(8,45)*four-r3(8,46)-r3(8,49)-r3(8,50)*four-r3(8,51)+r2(3,10) &
&                  +r2(3,11)*four+r2(3,12)+r2(3,18)+r2(3,19)*four+r2(3,20)+r2(2,43) &
&                  +r2(2,44)*four+r2(2,45)+r2(2,53)+r2(3,53)+r2(2,54)*four+r2(3,54)*four &
&                  +r2(2,55)+r2(3,55)-r1(3,32)-r1(3,33)*four-r1(3,34)-r1(3,37)-r1(3,38)*four &
&                  -r1(3,39)+r0(12)+r0(13)*four+r0(14)+r0(22)+r0(23)*four+r0(24))*zz+( &
&                  +r5(19,23)*two+r5(19,24)*two-r4(13,32)*two-r4(13,33)*two-r4(13,36)*two &
&                  -r4(13,37)*two+r3(8,33)*two+r3(8,34)*two+r3(8,41)*two+r3(10,41)*two &
&                  +r3(8,42)*two+r3(10,42)*two-r2(3,31)*two-r2(3,32)*two-r2(3,35)*two &
&                  -r2(3,36)*two+r1(3,11)*two+r1(3,12)*two+r1(3,19)*two+r1(3,20)*two)*zzz &
&                  +rxyz(1)*zzzz
      eri(4,3,6,6)=r112+rxyz(18)*qx+(+r7(25,4)*two-r6(18,6)*two-r6(18,8)*two+r5(12,8)*two &
&                  +r5(12,12)*two+r5(14,12)*two-r4(9,11)*two-r4(9,13)*two+r3(5,6)*two &
&                  +r3(5,10)*two)*qz+(+r6(25,11)*two-r5(18,15)*two-r5(18,19)*two &
&                  +r4(12,20)*two+r4(12,28)*two+r4(14,28)*two-r3(9,21)*two-r3(9,25)*two &
&                  +r2(6,11)*two+r2(6,19)*two)*xz+(+r6(18,12)-r5(12,16)-r5(12,20)+r4(7,21) &
&                  +r4(7,29)+r4(9,29)-r3(5,22)-r3(5,26)+r2(4,12)+r2(4,20))*zz+rxyz(2)*xzz
      eri(5,3,6,6)=r103+(+r7(34,3)-r6(26,5)-r6(26,7)+r5(19,7)+r5(19,11)+r5(21,11) &
&                  +r5(19,21)*three-r4(15,10)-r4(15,12)-r4(13,30)*three-r4(13,34)*three &
&                  +r3(10,5)+r3(10,9)+r3(8,31)*three+r3(8,39)*three+r3(10,39)*three &
&                  -r2(3,29)*three-r2(3,33)*three+r1(3,9)*three+r1(3,17)*three)*qx+(+r7(26,3) &
&                  +r7(26,4)*two-r6(19,5)-r6(19,6)*two-r6(19,7)-r6(19,8)*two+r5(13,7) &
&                  +r5(13,8)*two+r5(13,11)+r5(15,11)+r5(13,12)*two+r5(15,12)*two+r5(13,21) &
&                  +r5(13,22)*two-r4(10,10)-r4(10,11)*two-r4(10,12)-r4(10,13)*two-r4(8,30) &
&                  -r4(8,31)*two-r4(8,34)-r4(8,35)*two+r3(6,5)+r3(6,6)*two+r3(6,9) &
&                  +r3(6,10)*two+r3(4,31)+r3(4,32)*two+r3(4,39)+r3(6,39)+r3(4,40)*two &
&                  +r3(6,40)*two-r2(5,29)-r2(5,30)*two-r2(5,33)-r2(5,34)*two+r1(1,9) &
&                  +r1(1,10)*two+r1(1,17)+r1(1,18)*two)*qz+(+r6(26,10)+r6(26,11)*two &
&                  -r5(19,14)-r5(19,15)*two-r5(19,18)-r5(19,19)*two+r4(13,19)+r4(13,20)*two &
&                  +r4(13,27)+r4(15,27)+r4(13,28)*two+r4(15,28)*two+r4(13,39)+r4(13,40)*two &
&                  -r3(10,20)-r3(10,21)*two-r3(10,24)-r3(10,25)*two-r3(8,44)-r3(8,45)*two &
&                  -r3(8,49)-r3(8,50)*two+r2(3,10)+r2(3,11)*two+r2(3,18)+r2(3,19)*two &
&                  +r2(2,43)+r2(2,44)*two+r2(2,53)+r2(3,53)+r2(2,54)*two+r2(3,54)*two &
&                  -r1(3,32)-r1(3,33)*two-r1(3,37)-r1(3,38)*two+r0(12)+r0(13)*two+r0(22) &
&                  +r0(23)*two)*xz+(+r6(19,11)*two+r6(19,12)-r5(13,15)*two-r5(13,16) &
&                  -r5(13,19)*two-r5(13,20)+r4(8,20)*two+r4(8,21)+r4(8,28)*two+r4(10,28)*two &
&                  +r4(8,29)+r4(10,29)-r3(6,21)*two-r3(6,22)-r3(6,25)*two-r3(6,26) &
&                  +r2(5,11)*two+r2(5,12)+r2(5,19)*two+r2(5,20))*zz+(+r5(19,23)*two+r5(19,24) &
&                  -r4(13,32)*two-r4(13,33)-r4(13,36)*two-r4(13,37)+r3(8,33)*two+r3(8,34) &
&                  +r3(8,41)*two+r3(10,41)*two+r3(8,42)+r3(10,42)-r2(3,31)*two-r2(3,32) &
&                  -r2(3,35)*two-r2(3,36)+r1(3,11)*two+r1(3,12)+r1(3,19)*two+r1(3,20))*xzz+( &
&                  +r5(13,24)-r4(8,33)-r4(8,37)+r3(4,34)+r3(4,42)+r3(6,42)-r2(5,32)-r2(5,36) &
&                  +r1(1,12)+r1(1,20))*zzz+rxyz(1)*xzzz
      eri(6,3,6,6)=r013+(+r7(33,3)+r7(33,4)*two-r6(25,5)-r6(25,6)*two-r6(25,7) &
&                  -r6(25,8)*two+r5(18,7)+r5(18,8)*two+r5(18,11)+r5(20,11)+r5(18,12)*two &
&                  +r5(20,12)*two+r5(18,21)+r5(18,22)*two-r4(14,10)-r4(14,11)*two-r4(14,12) &
&                  -r4(14,13)*two-r4(12,30)-r4(12,31)*two-r4(12,34)-r4(12,35)*two+r3(9,5) &
&                  +r3(9,6)*two+r3(9,9)+r3(9,10)*two+r3(7,31)+r3(7,32)*two+r3(7,39)+r3(9,39) &
&                  +r3(7,40)*two+r3(9,40)*two-r2(6,29)-r2(6,30)*two-r2(6,33)-r2(6,34)*two &
&                  +r1(2,9)+r1(2,10)*two+r1(2,17)+r1(2,18)*two)*qz+(+r6(25,11)*two+r6(25,12) &
&                  -r5(18,15)*two-r5(18,16)-r5(18,19)*two-r5(18,20)+r4(12,20)*two+r4(12,21) &
&                  +r4(12,28)*two+r4(14,28)*two+r4(12,29)+r4(14,29)-r3(9,21)*two-r3(9,22) &
&                  -r3(9,25)*two-r3(9,26)+r2(6,11)*two+r2(6,12)+r2(6,19)*two+r2(6,20))*zz &
&                  +rxyz(2)*zzz
      eri(1,4,6,6)=r310+(+r7(18,3)*two+r7(18,4)-r6(12,5)*two-r6(12,6)-r6(12,7)*two &
&                  -r6(12,8)+r5(7,7)*two+r5(7,8)+r5(7,11)*two+r5(9,11)*two+r5(7,12)+r5(9,12) &
&                  +r5(18,21)*two+r5(18,22)-r4(5,10)*two-r4(5,11)-r4(5,12)*two-r4(5,13) &
&                  -r4(12,30)*two-r4(12,31)-r4(12,34)*two-r4(12,35)+r3(2,5)*two+r3(2,6) &
&                  +r3(2,9)*two+r3(2,10)+r3(7,31)*two+r3(7,32)+r3(7,39)*two+r3(9,39)*two &
&                  +r3(7,40)+r3(9,40)-r2(6,29)*two-r2(6,30)-r2(6,33)*two-r2(6,34)+r1(2,9)*two &
&                  +r1(2,10)+r1(2,17)*two+r1(2,18))*qx+(+r6(18,10)+r6(18,11)*two-r5(12,14) &
&                  -r5(12,15)*two-r5(12,18)-r5(12,19)*two+r4(7,19)+r4(7,20)*two+r4(7,27) &
&                  +r4(9,27)+r4(7,28)*two+r4(9,28)*two-r3(5,20)-r3(5,21)*two-r3(5,24) &
&                  -r3(5,25)*two+r2(4,10)+r2(4,11)*two+r2(4,18)+r2(4,19)*two)*xx+rxyz(3)*xxx
      eri(2,4,6,6)=r130+rxyz(8)*qx
      eri(3,4,6,6)=r112+rxyz(19)*qx+(+r7(25,3)*two-r6(18,5)*two-r6(18,7)*two+r5(12,7)*two &
&                  +r5(12,11)*two+r5(14,11)*two-r4(9,10)*two-r4(9,12)*two+r3(5,5)*two &
&                  +r3(5,9)*two)*qz+(+r6(25,11)*two-r5(18,15)*two-r5(18,19)*two+r4(12,20)*two &
&                  +r4(12,28)*two+r4(14,28)*two-r3(9,21)*two-r3(9,25)*two+r2(6,11)*two &
&                  +r2(6,19)*two)*xz+(+r6(18,10)-r5(12,14)-r5(12,18)+r4(7,19)+r4(7,27) &
&                  +r4(9,27)-r3(5,20)-r3(5,24)+r2(4,10)+r2(4,18))*zz+rxyz(3)*xzz
      eri(4,4,6,6)=r220+(+r7(24,3)+r7(24,4)-r6(17,5)-r6(17,6)-r6(17,7)-r6(17,8)+r5(11,7) &
&                  +r5(11,8)+r5(11,11)+r5(13,11)+r5(11,12)+r5(13,12)+r5(13,21)+r5(13,22) &
&                  -r4(8,10)-r4(8,11)-r4(8,12)-r4(8,13)-r4(8,30)-r4(8,31)-r4(8,34)-r4(8,35) &
&                  +r3(4,5)+r3(4,6)+r3(4,9)+r3(4,10)+r3(4,31)+r3(4,32)+r3(4,39)+r3(6,39) &
&                  +r3(4,40)+r3(6,40)-r2(5,29)-r2(5,30)-r2(5,33)-r2(5,34)+r1(1,9)+r1(1,10) &
&                  +r1(1,17)+r1(1,18))*qx+rxyz(6)*xx
      eri(5,4,6,6)=r211+rxyz(9)*qx+rxyz(17)*qz+rxyz(10)*xx+(+r6(18,10)+r6(18,11) &
&                  -r5(12,14)-r5(12,15)-r5(12,18)-r5(12,19)+r4(7,19)+r4(7,20)+r4(7,27) &
&                  +r4(9,27)+r4(7,28)+r4(9,28)-r3(5,20)-r3(5,21)-r3(5,24)-r3(5,25)+r2(4,10) &
&                  +r2(4,11)+r2(4,18)+r2(4,19))*xz+rxyz(3)*xxz
      eri(6,4,6,6)=r121+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,6,6)=r301+(+r7(19,3)*two+r7(19,4)-r6(13,5)*two-r6(13,6)-r6(13,7)*two &
&                  -r6(13,8)+r5(8,7)*two+r5(8,8)+r5(8,11)*two+r5(10,11)*two+r5(8,12) &
&                  +r5(10,12)+r5(19,21)*two+r5(19,22)-r4(6,10)*two-r4(6,11)-r4(6,12)*two &
&                  -r4(6,13)-r4(13,30)*two-r4(13,31)-r4(13,34)*two-r4(13,35)+r3(3,5)*two &
&                  +r3(3,6)+r3(3,9)*two+r3(3,10)+r3(8,31)*two+r3(8,32)+r3(8,39)*two &
&                  +r3(10,39)*two+r3(8,40)+r3(10,40)-r2(3,29)*two-r2(3,30)-r2(3,33)*two &
&                  -r2(3,34)+r1(3,9)*two+r1(3,10)+r1(3,17)*two+r1(3,18))*qx+(+r7(13,4) &
&                  -r6(8,6)-r6(8,8)+r5(4,8)+r5(4,12)+r5(6,12)+r5(13,22)*three-r4(3,11) &
&                  -r4(3,13)-r4(8,31)*three-r4(8,35)*three+r3(1,6)+r3(1,10)+r3(4,32)*three &
&                  +r3(4,40)*three+r3(6,40)*three-r2(5,30)*three-r2(5,34)*three &
&                  +r1(1,10)*three+r1(1,18)*three)*qz+(+r6(19,10)+r6(19,11)*two-r5(13,14) &
&                  -r5(13,15)*two-r5(13,18)-r5(13,19)*two+r4(8,19)+r4(8,20)*two+r4(8,27) &
&                  +r4(10,27)+r4(8,28)*two+r4(10,28)*two-r3(6,20)-r3(6,21)*two-r3(6,24) &
&                  -r3(6,25)*two+r2(5,10)+r2(5,11)*two+r2(5,18)+r2(5,19)*two)*xx+( &
&                  +r6(13,11)*two+r6(13,12)-r5(8,15)*two-r5(8,16)-r5(8,19)*two-r5(8,20) &
&                  +r4(4,20)*two+r4(4,21)+r4(4,28)*two+r4(6,28)*two+r4(4,29)+r4(6,29) &
&                  +r4(13,40)*two+r4(13,41)-r3(3,21)*two-r3(3,22)-r3(3,25)*two-r3(3,26) &
&                  -r3(8,45)*two-r3(8,46)-r3(8,50)*two-r3(8,51)+r2(1,11)*two+r2(1,12) &
&                  +r2(1,19)*two+r2(1,20)+r2(2,44)*two+r2(2,45)+r2(2,54)*two+r2(3,54)*two &
&                  +r2(2,55)+r2(3,55)-r1(3,33)*two-r1(3,34)-r1(3,38)*two-r1(3,39)+r0(13)*two &
&                  +r0(14)+r0(23)*two+r0(24))*xz+(+r5(19,23)-r4(13,32)-r4(13,36)+r3(8,33) &
&                  +r3(8,41)+r3(10,41)-r2(3,31)-r2(3,35)+r1(3,11)+r1(3,19))*xxx+(+r5(13,23) &
&                  +r5(13,24)*two-r4(8,32)-r4(8,33)*two-r4(8,36)-r4(8,37)*two+r3(4,33) &
&                  +r3(4,34)*two+r3(4,41)+r3(6,41)+r3(4,42)*two+r3(6,42)*two-r2(5,31) &
&                  -r2(5,32)*two-r2(5,35)-r2(5,36)*two+r1(1,11)+r1(1,12)*two+r1(1,19) &
&                  +r1(1,20)*two)*xxz+rxyz(1)*xxxz
      eri(2,5,6,6)=r121+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,6,6)=r103+(+r7(34,4)-r6(26,6)-r6(26,8)+r5(19,8)+r5(19,12)+r5(21,12) &
&                  +r5(19,22)*three-r4(15,11)-r4(15,13)-r4(13,31)*three-r4(13,35)*three &
&                  +r3(10,6)+r3(10,10)+r3(8,32)*three+r3(8,40)*three+r3(10,40)*three &
&                  -r2(3,30)*three-r2(3,34)*three+r1(3,10)*three+r1(3,18)*three)*qx+( &
&                  +r7(26,3)*two+r7(26,4)-r6(19,5)*two-r6(19,6)-r6(19,7)*two-r6(19,8) &
&                  +r5(13,7)*two+r5(13,8)+r5(13,11)*two+r5(15,11)*two+r5(13,12)+r5(15,12) &
&                  +r5(13,21)*two+r5(13,22)-r4(10,10)*two-r4(10,11)-r4(10,12)*two-r4(10,13) &
&                  -r4(8,30)*two-r4(8,31)-r4(8,34)*two-r4(8,35)+r3(6,5)*two+r3(6,6) &
&                  +r3(6,9)*two+r3(6,10)+r3(4,31)*two+r3(4,32)+r3(4,39)*two+r3(6,39)*two &
&                  +r3(4,40)+r3(6,40)-r2(5,29)*two-r2(5,30)-r2(5,33)*two-r2(5,34)+r1(1,9)*two &
&                  +r1(1,10)+r1(1,17)*two+r1(1,18))*qz+(+r6(26,11)*two+r6(26,12) &
&                  -r5(19,15)*two-r5(19,16)-r5(19,19)*two-r5(19,20)+r4(13,20)*two+r4(13,21) &
&                  +r4(13,28)*two+r4(15,28)*two+r4(13,29)+r4(15,29)+r4(13,40)*two+r4(13,41) &
&                  -r3(10,21)*two-r3(10,22)-r3(10,25)*two-r3(10,26)-r3(8,45)*two-r3(8,46) &
&                  -r3(8,50)*two-r3(8,51)+r2(3,11)*two+r2(3,12)+r2(3,19)*two+r2(3,20) &
&                  +r2(2,44)*two+r2(2,45)+r2(2,54)*two+r2(3,54)*two+r2(2,55)+r2(3,55) &
&                  -r1(3,33)*two-r1(3,34)-r1(3,38)*two-r1(3,39)+r0(13)*two+r0(14)+r0(23)*two &
&                  +r0(24))*xz+(+r6(19,10)+r6(19,11)*two-r5(13,14)-r5(13,15)*two-r5(13,18) &
&                  -r5(13,19)*two+r4(8,19)+r4(8,20)*two+r4(8,27)+r4(10,27)+r4(8,28)*two &
&                  +r4(10,28)*two-r3(6,20)-r3(6,21)*two-r3(6,24)-r3(6,25)*two+r2(5,10) &
&                  +r2(5,11)*two+r2(5,18)+r2(5,19)*two)*zz+(+r5(19,23)+r5(19,24)*two &
&                  -r4(13,32)-r4(13,33)*two-r4(13,36)-r4(13,37)*two+r3(8,33)+r3(8,34)*two &
&                  +r3(8,41)+r3(10,41)+r3(8,42)*two+r3(10,42)*two-r2(3,31)-r2(3,32)*two &
&                  -r2(3,35)-r2(3,36)*two+r1(3,11)+r1(3,12)*two+r1(3,19)+r1(3,20)*two)*xzz+( &
&                  +r5(13,23)-r4(8,32)-r4(8,36)+r3(4,33)+r3(4,41)+r3(6,41)-r2(5,31)-r2(5,35) &
&                  +r1(1,11)+r1(1,19))*zzz+rxyz(1)*xzzz
      eri(4,5,6,6)=r211+rxyz(9)*qx+rxyz(16)*qz+rxyz(10)*xx+(+r6(18,11)+r6(18,12) &
&                  -r5(12,15)-r5(12,16)-r5(12,19)-r5(12,20)+r4(7,20)+r4(7,21)+r4(7,28) &
&                  +r4(9,28)+r4(7,29)+r4(9,29)-r3(5,21)-r3(5,22)-r3(5,25)-r3(5,26)+r2(4,11) &
&                  +r2(4,12)+r2(4,19)+r2(4,20))*xz+rxyz(2)*xxz
      eri(5,5,6,6)=r202+(+r7(26,3)+r7(26,4)-r6(19,5)-r6(19,6)-r6(19,7)-r6(19,8)+r5(13,7) &
&                  +r5(13,8)+r5(13,11)+r5(15,11)+r5(13,12)+r5(15,12)+r5(13,21)+r5(13,22) &
&                  -r4(10,10)-r4(10,11)-r4(10,12)-r4(10,13)-r4(8,30)-r4(8,31)-r4(8,34) &
&                  -r4(8,35)+r3(6,5)+r3(6,6)+r3(6,9)+r3(6,10)+r3(4,31)+r3(4,32)+r3(4,39) &
&                  +r3(6,39)+r3(4,40)+r3(6,40)-r2(5,29)-r2(5,30)-r2(5,33)-r2(5,34)+r1(1,9) &
&                  +r1(1,10)+r1(1,17)+r1(1,18))*qx+(+r7(19,3)+r7(19,4)-r6(13,5)-r6(13,6) &
&                  -r6(13,7)-r6(13,8)+r5(8,7)+r5(8,8)+r5(8,11)+r5(10,11)+r5(8,12)+r5(10,12) &
&                  +r5(19,21)+r5(19,22)-r4(6,10)-r4(6,11)-r4(6,12)-r4(6,13)-r4(13,30) &
&                  -r4(13,31)-r4(13,34)-r4(13,35)+r3(3,5)+r3(3,6)+r3(3,9)+r3(3,10)+r3(8,31) &
&                  +r3(8,32)+r3(8,39)+r3(10,39)+r3(8,40)+r3(10,40)-r2(3,29)-r2(3,30)-r2(3,33) &
&                  -r2(3,34)+r1(3,9)+r1(3,10)+r1(3,17)+r1(3,18))*qz+(+r6(26,11)-r5(19,15) &
&                  -r5(19,19)+r4(13,20)+r4(13,28)+r4(15,28)+r4(13,40)-r3(10,21)-r3(10,25) &
&                  -r3(8,45)-r3(8,50)+r2(3,11)+r2(3,19)+r2(2,44)+r2(2,54)+r2(3,54)-r1(3,33) &
&                  -r1(3,38)+r0(13)+r0(23))*xx+(+r6(19,10)+r6(19,11)*two+r6(19,12)-r5(13,14) &
&                  -r5(13,15)*two-r5(13,16)-r5(13,18)-r5(13,19)*two-r5(13,20)+r4(8,19) &
&                  +r4(8,20)*two+r4(8,21)+r4(8,27)+r4(10,27)+r4(8,28)*two+r4(10,28)*two &
&                  +r4(8,29)+r4(10,29)-r3(6,20)-r3(6,21)*two-r3(6,22)-r3(6,24)-r3(6,25)*two &
&                  -r3(6,26)+r2(5,10)+r2(5,11)*two+r2(5,12)+r2(5,18)+r2(5,19)*two+r2(5,20)) &
&                  *xz+(+r6(13,11)-r5(8,15)-r5(8,19)+r4(4,20)+r4(4,28)+r4(6,28)+r4(13,40) &
&                  -r3(3,21)-r3(3,25)-r3(8,45)-r3(8,50)+r2(1,11)+r2(1,19)+r2(2,44)+r2(2,54) &
&                  +r2(3,54)-r1(3,33)-r1(3,38)+r0(13)+r0(23))*zz+(+r5(19,23)+r5(19,24) &
&                  -r4(13,32)-r4(13,33)-r4(13,36)-r4(13,37)+r3(8,33)+r3(8,34)+r3(8,41) &
&                  +r3(10,41)+r3(8,42)+r3(10,42)-r2(3,31)-r2(3,32)-r2(3,35)-r2(3,36)+r1(3,11) &
&                  +r1(3,12)+r1(3,19)+r1(3,20))*xxz+(+r5(13,23)+r5(13,24)-r4(8,32)-r4(8,33) &
&                  -r4(8,36)-r4(8,37)+r3(4,33)+r3(4,34)+r3(4,41)+r3(6,41)+r3(4,42)+r3(6,42) &
&                  -r2(5,31)-r2(5,32)-r2(5,35)-r2(5,36)+r1(1,11)+r1(1,12)+r1(1,19)+r1(1,20)) &
&                  *xzz+rxyz(1)*xxzz
      eri(6,5,6,6)=r112+rxyz(19)*qx+(+r7(25,3)+r7(25,4)-r6(18,5)-r6(18,6)-r6(18,7) &
&                  -r6(18,8)+r5(12,7)+r5(12,8)+r5(12,11)+r5(14,11)+r5(12,12)+r5(14,12) &
&                  -r4(9,10)-r4(9,11)-r4(9,12)-r4(9,13)+r3(5,5)+r3(5,6)+r3(5,9)+r3(5,10))*qz &
&                  +(+r6(25,11)+r6(25,12)-r5(18,15)-r5(18,16)-r5(18,19)-r5(18,20)+r4(12,20) &
&                  +r4(12,21)+r4(12,28)+r4(14,28)+r4(12,29)+r4(14,29)-r3(9,21)-r3(9,22) &
&                  -r3(9,25)-r3(9,26)+r2(6,11)+r2(6,12)+r2(6,19)+r2(6,20))*xz+rxyz(11)*zz &
&                  +rxyz(2)*xzz
      eri(1,6,6,6)=r211+(+r7(25,3)*two-r6(18,5)*two-r6(18,7)*two+r5(12,7)*two &
&                  +r5(12,11)*two+r5(14,11)*two-r4(9,10)*two-r4(9,12)*two+r3(5,5)*two &
&                  +r3(5,9)*two)*qx+rxyz(16)*qz+(+r6(25,10)-r5(18,14)-r5(18,18)+r4(12,19) &
&                  +r4(12,27)+r4(14,27)-r3(9,20)-r3(9,24)+r2(6,10)+r2(6,18))*xx+( &
&                  +r6(18,11)*two-r5(12,15)*two-r5(12,19)*two+r4(7,20)*two+r4(7,28)*two &
&                  +r4(9,28)*two-r3(5,21)*two-r3(5,25)*two+r2(4,11)*two+r2(4,19)*two)*xz &
&                  +rxyz(3)*xxz
      eri(2,6,6,6)=r031+rxyz(8)*qz
      eri(3,6,6,6)=r013+(+r7(33,3)*two+r7(33,4)-r6(25,5)*two-r6(25,6)-r6(25,7)*two &
&                  -r6(25,8)+r5(18,7)*two+r5(18,8)+r5(18,11)*two+r5(20,11)*two+r5(18,12) &
&                  +r5(20,12)+r5(18,21)*two+r5(18,22)-r4(14,10)*two-r4(14,11)-r4(14,12)*two &
&                  -r4(14,13)-r4(12,30)*two-r4(12,31)-r4(12,34)*two-r4(12,35)+r3(9,5)*two &
&                  +r3(9,6)+r3(9,9)*two+r3(9,10)+r3(7,31)*two+r3(7,32)+r3(7,39)*two &
&                  +r3(9,39)*two+r3(7,40)+r3(9,40)-r2(6,29)*two-r2(6,30)-r2(6,33)*two &
&                  -r2(6,34)+r1(2,9)*two+r1(2,10)+r1(2,17)*two+r1(2,18))*qz+(+r6(25,10) &
&                  +r6(25,11)*two-r5(18,14)-r5(18,15)*two-r5(18,18)-r5(18,19)*two+r4(12,19) &
&                  +r4(12,20)*two+r4(12,27)+r4(14,27)+r4(12,28)*two+r4(14,28)*two-r3(9,20) &
&                  -r3(9,21)*two-r3(9,24)-r3(9,25)*two+r2(6,10)+r2(6,11)*two+r2(6,18) &
&                  +r2(6,19)*two)*zz+rxyz(3)*zzz
      eri(4,6,6,6)=r121+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,6,6)=r112+rxyz(18)*qx+(+r7(25,3)+r7(25,4)-r6(18,5)-r6(18,6)-r6(18,7) &
&                  -r6(18,8)+r5(12,7)+r5(12,8)+r5(12,11)+r5(14,11)+r5(12,12)+r5(14,12) &
&                  -r4(9,10)-r4(9,11)-r4(9,12)-r4(9,13)+r3(5,5)+r3(5,6)+r3(5,9)+r3(5,10))*qz &
&                  +(+r6(25,10)+r6(25,11)-r5(18,14)-r5(18,15)-r5(18,18)-r5(18,19)+r4(12,19) &
&                  +r4(12,20)+r4(12,27)+r4(14,27)+r4(12,28)+r4(14,28)-r3(9,20)-r3(9,21) &
&                  -r3(9,24)-r3(9,25)+r2(6,10)+r2(6,11)+r2(6,18)+r2(6,19))*xz+rxyz(11)*zz &
&                  +rxyz(3)*xzz
      eri(6,6,6,6)=r022+(+r7(32,3)+r7(32,4)-r6(24,5)-r6(24,6)-r6(24,7)-r6(24,8)+r5(17,7) &
&                  +r5(17,8)+r5(17,11)+r5(19,11)+r5(17,12)+r5(19,12)+r5(19,21)+r5(19,22) &
&                  -r4(13,10)-r4(13,11)-r4(13,12)-r4(13,13)-r4(13,30)-r4(13,31)-r4(13,34) &
&                  -r4(13,35)+r3(8,5)+r3(8,6)+r3(8,9)+r3(8,10)+r3(8,31)+r3(8,32)+r3(8,39) &
&                  +r3(10,39)+r3(8,40)+r3(10,40)-r2(3,29)-r2(3,30)-r2(3,33)-r2(3,34)+r1(3,9) &
&                  +r1(3,10)+r1(3,17)+r1(3,18))*qz+rxyz(6)*zz
      return
end
