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
  subroutine int2dddp(phmdint,exfac1,exfac2,xyziq,xzkl,rot,nijkl,nbfijkl)
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
      real(8),parameter :: p11=1.1D+01, p12=1.2D+01, p13=1.3D+01, p15=1.5D+01, p18=1.8D+01
      real(8),parameter :: p21=2.1D+01, p45=4.5D+01, p105=1.05D+2
      real(8),parameter :: sqrt3=1.73205080756888D+00, sqrt3h=8.660254037844386D-01
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*), xzkl(2), rot(3,3)
      real(8),intent(out) :: phmdint(6,6,6,6)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:7)
      real(8) :: f0(3), f1(2,8), f2(3,9), f3(4,9), f4(5,9), f5(6,6), f6(7,3), f7(8), ftw(7,9)
      real(8) :: r0(15), r1(3,27), r2(6,34), r3(10,31), r4(15,21), r5(21,11), r6(28,4), r7(36)
      real(8) :: ex12, ex34, ex43, ex14, ex41, expq, expq2, expq3, expq4, ex3q, ex4q, c12, c34
      real(8) :: zip, ziq, xiq, yiq, xiq2, yiq2, xyiq, xiq4, yiq4, xyiq2, xiq6, yiq6, x4y2
      real(8) :: x2y4, xypq2, zpq, zpq2, zpq3, zpq4, zpq5, zpq6, fac, ex33q, ex34q, ex44q, zjp
      real(8) :: pmd, pmd2, qmd, qmd2, qmd3, qmd4, qmd4x, qmd4y, qmd4xy, qx, qz
      real(8) :: eri(6,6,6,3), work(15), rot2(6,6), rot3(6,5)
      real(8) :: f1w(3,6), f2w(6,8), f3w(10,9), f4w(15,9), f5w(21,6), f6w(28,3)
!
! Zero-clear
!
      r0(1:15)     = zero
      r1(1:3 ,1:27)= zero
      r2(1:6 ,1:34)= zero
      r3(1:10,1:31)= zero
      r4(1:15,1:21)= zero
      r5(1:21,1:11)= zero
      r6(1:28,1:4) = zero
      r7(1:36)     = zero
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
        x4y2 = xiq4*yiq2
        x2y4 = xiq2*yiq4
        xypq2= xiq2+yiq2
        f0(1:3)    = zero
        f1(1:2,1:8)= zero
        f2(1:3,1:9)= zero
        f3(1:4,1:9)= zero
        f4(1:5,1:9)= zero
        f5(1:6,1:6)= zero
        f6(1:7,1:3)= zero
        f7(1:8)    = zero
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
            do ii= 0,7
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
          endif
          zpq3= zpq2*zpq
          zpq4= zpq2*zpq2
          zpq5= zpq2*zpq3
          zpq6= zpq3*zpq3
          pmd2= pmd*pmd
          ftw(1,1)= zjp*zjp*zip
          ftw(1,2)= pmd*zjp
          ftw(1,3)= pmd*zip
          ftw(1,4)= pmd*zjp*zjp
          ftw(1,5)= pmd*zjp*zip
          ftw(1,6)= pmd2
          ftw(1,7)= pmd2*zjp
          ftw(1,8)= pmd2*zip
          ftw(1,9)= pmd2*pmd
          do i= 1,9
            ftw(2,i)= ftw(1,i)*zpq
            ftw(3,i)= ftw(1,i)*zpq2
            ftw(4,i)= ftw(1,i)*zpq3
            ftw(5,i)= ftw(1,i)*zpq4
            ftw(6,i)= ftw(1,i)*zpq5
            ftw(7,i)= ftw(1,i)*zpq6
          enddo
          f0(1)= f0(1)+ft(0)*ftw(1,1)
          f0(2)= f0(2)+ft(0)*ftw(1,2)
          f0(3)= f0(3)+ft(0)*ftw(1,3)
          do i= 1,6
            f1(1,i)= f1(1,i)-ft(1)*ftw(1,i)
            f1(2,i)= f1(2,i)-ft(1)*ftw(2,i)
          enddo
          f1(1,7)= f1(1,7)-ft(1)*ftw(1,7)
          f1(1,8)= f1(1,8)-ft(1)*ftw(1,8)
          do i= 1,9
            f2(1,i)= f2(1,i)+ft(2)*ftw(1,i)
            f2(2,i)= f2(2,i)+ft(2)*ftw(2,i)
            f2(3,i)= f2(3,i)+ft(2)*ftw(3,i)
          enddo
          do i= 1,9
            do j= 1,4
              f3(j,i)= f3(j,i)-ft(3)*ftw(j,i)
            enddo
          enddo
          do i= 1,9
            do j= 1,5
              f4(j,i)= f4(j,i)+ft(4)*ftw(j,i)
            enddo
          enddo
          do i= 1,6
            do j= 1,6
              f5(j,i)= f5(j,i)-ft(5)*ftw(j,i+3)
            enddo
          enddo
          do i= 1,3
            do j= 1,7
              f6(j,i)= f6(j,i)+ft(6)*ftw(j,i+6)
            enddo
          enddo
          f7(1)= f7(1)-ft(7)*ftw(1,9)
          f7(2)= f7(2)-ft(7)*ftw(2,9)
          f7(3)= f7(3)-ft(7)*ftw(3,9)
          f7(4)= f7(4)-ft(7)*ftw(4,9)
          f7(5)= f7(5)-ft(7)*ftw(5,9)
          f7(6)= f7(6)-ft(7)*ftw(6,9)
          f7(7)= f7(7)-ft(7)*ftw(7,9)
          f7(8)= f7(8)-ft(7)*ftw(7,9)*zpq
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
        enddo
!
        do i= 1,6
          f1w(1,i)= f1(1,i)*xiq
          f1w(2,i)= f1(1,i)*yiq
          f1w(3,i)= f1(2,i)
        enddo
        do i= 1,4
          do j= 1,3
            r1(j,i  )= r1(j,i  )+f1w(j,1)*work(i+5)
            r1(j,i+4)= r1(j,i+4)+f1w(j,2)*work(i+5)
            r1(j,i+8)= r1(j,i+8)+f1w(j,3)*work(i+5)
          enddo
        enddo
        do i= 1,5
          do j= 1,3
            r1(j,i+12)= r1(j,i+12)+f1w(j,4)*work(i)
            r1(j,i+17)= r1(j,i+17)+f1w(j,5)*work(i)
            r1(j,i+22)= r1(j,i+22)+f1w(j,6)*work(i)
          enddo
        enddo
!
        do i= 1,8
          f2w(1,i)= f2(1,i)*xiq2+f1(1,i)
          f2w(2,i)= f2(1,i)*yiq2+f1(1,i)
          f2w(3,i)= f2(3,i)     +f1(1,i)
          f2w(4,i)= f2(1,i)*xyiq
          f2w(5,i)= f2(2,i)*xiq
          f2w(6,i)= f2(2,i)*yiq
        enddo
        do i= 1,4
          do j= 1,6
            r2(j,i  )= r2(j,i  )+f2w(j,1)*work(i+9)
            r2(j,i+4)= r2(j,i+4)+f2w(j,2)*work(i+9)
            r2(j,i+8)= r2(j,i+8)+f2w(j,3)*work(i+9)
          enddo
        enddo
        do i= 1,4
          do j= 1,6
            r2(j,i+12)= r2(j,i+12)+f2w(j,4)*work(i+5)
            r2(j,i+16)= r2(j,i+16)+f2w(j,5)*work(i+5)
            r2(j,i+20)= r2(j,i+20)+f2w(j,6)*work(i+5)
          enddo
        enddo
        do i= 1,5
          do j= 1,6
            r2(j,i+24)= r2(j,i+24)+f2w(j,7)*work(i)
            r2(j,i+29)= r2(j,i+29)+f2w(j,8)*work(i)
          enddo
        enddo
! 
        do i= 1,9
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
            r3(j,i  )= r3(j,i  )+f3w(j,1)*work(i+13)
            r3(j,i+2)= r3(j,i+2)+f3w(j,2)*work(i+13)
            r3(j,i+4)= r3(j,i+4)+f3w(j,3)*work(i+13)
          enddo
        enddo
        do i= 1,4
          do j= 1,10
            r3(j,i+ 6)= r3(j,i+ 6)+f3w(j,4)*work(i+9)
            r3(j,i+10)= r3(j,i+10)+f3w(j,5)*work(i+9)
            r3(j,i+14)= r3(j,i+14)+f3w(j,6)*work(i+9)
          enddo
        enddo
        do i= 1,4
          do j= 1,10
            r3(j,i+18)= r3(j,i+18)+f3w(j,7)*work(i+5)
            r3(j,i+22)= r3(j,i+22)+f3w(j,8)*work(i+5)
          enddo
        enddo
        do i= 27,31
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,9)*work(i-26)
          enddo
        enddo
!
        do i= 1,9
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
        do j= 1,15
          r4(j,1)= r4(j,1)+f4w(j,1)*qmd4
          r4(j,2)= r4(j,2)+f4w(j,2)*qmd4
          r4(j,3)= r4(j,3)+f4w(j,3)*qmd4
        enddo
        do i= 1,2
          do j= 1,15
            r4(j,i+3)= r4(j,i+3)+f4w(j,4)*work(i+13)
            r4(j,i+5)= r4(j,i+5)+f4w(j,5)*work(i+13)
            r4(j,i+7)= r4(j,i+7)+f4w(j,6)*work(i+13)
          enddo
        enddo
        do i= 1,4
          do j= 1,15
            r4(j,i+ 9)= r4(j,i+ 9)+f4w(j,7)*work(i+9)
            r4(j,i+13)= r4(j,i+13)+f4w(j,8)*work(i+9)
          enddo
        enddo
        do i= 18,21
          do j= 1,15
            r4(j,i)= r4(j,i)+f4w(j,9)*work(i-12)
          enddo
        enddo
!
        do i= 1,6
          f5w( 1,i)=(f5(1,i)*xiq4 +f4(1,i+3)*xiq2*ten                 +f3(1,i+3)*p15  )*xiq
          f5w( 2,i)=(f5(1,i)*xiq4 +f4(1,i+3)*xiq2*six                 +f3(1,i+3)*three)*yiq
          f5w( 3,i)=(f5(2,i)*xiq4 +f4(2,i+3)*xiq2*six                 +f3(2,i+3)*three)
          f5w( 4,i)=(f5(1,i)*xyiq2+f4(1,i+3)*xiq2+f4(1,i+3)*yiq2*three+f3(1,i+3)*three)*xiq
          f5w( 5,i)=(f5(2,i)*xiq2 +f4(2,i+3)*three                                    )*xyiq
          f5w( 6,i)=(f5(3,i)*xiq2 +f4(1,i+3)*xiq2+f4(3,i+3)*three     +f3(1,i+3)*three)*xiq
          f5w( 7,i)=(f5(1,i)*xyiq2+f4(1,i+3)*xiq2*three+f4(1,i+3)*yiq2+f3(1,i+3)*three)*yiq
          f5w( 8,i)=(f5(2,i)*xyiq2+f4(2,i+3)*xiq2+f4(2,i+3)*yiq2      +f3(2,i+3)      )
          f5w( 9,i)=(f5(3,i)*xiq2 +f4(1,i+3)*xiq2+f4(3,i+3)           +f3(1,i+3)      )*yiq
          f5w(10,i)=(f5(4,i)*xiq2 +f4(2,i+3)*xiq2*three+f4(4,i+3)     +f3(2,i+3)*three)
          f5w(11,i)=(f5(1,i)*yiq4 +f4(1,i+3)*yiq2*six                 +f3(1,i+3)*three)*xiq
          f5w(12,i)=(f5(2,i)*yiq2 +f4(2,i+3)*three                                    )*xyiq
          f5w(13,i)=(f5(3,i)*yiq2 +f4(1,i+3)*yiq2+f4(3,i+3)           +f3(1,i+3)      )*xiq
          f5w(14,i)=(f5(4,i)      +f4(2,i+3)*three                                    )*xyiq
          f5w(15,i)=(f5(5,i)      +f4(3,i+3)*six                      +f3(1,i+3)*three)*xiq
          f5w(16,i)=(f5(1,i)*yiq4 +f4(1,i+3)*yiq2*ten                 +f3(1,i+3)*p15  )*yiq
          f5w(17,i)=(f5(2,i)*yiq4 +f4(2,i+3)*yiq2*six                 +f3(2,i+3)*three)
          f5w(18,i)=(f5(3,i)*yiq2 +f4(1,i+3)*yiq2+f4(3,i+3)*three     +f3(1,i+3)*three)*yiq
          f5w(19,i)=(f5(4,i)*yiq2 +f4(2,i+3)*yiq2*three+f4(4,i+3)     +f3(2,i+3)*three)
          f5w(20,i)=(f5(5,i)      +f4(3,i+3)*six                      +f3(1,i+3)*three)*yiq
          f5w(21,i)=(f5(6,i)      +f4(4,i+3)*ten                      +f3(2,i+3)*p15  )
        enddo
        do j= 1,21
          r5(j,1)= r5(j,1)+f5w(j,1)*qmd4
          r5(j,2)= r5(j,2)+f5w(j,2)*qmd4
          r5(j,3)= r5(j,3)+f5w(j,3)*qmd4
        enddo
        do i= 1,2
          do j= 1,21
            r5(j,i+3)= r5(j,i+3)+f5w(j,4)*work(i+13)
            r5(j,i+5)= r5(j,i+5)+f5w(j,5)*work(i+13)
          enddo
        enddo
        do i= 8,11
          do j= 1,21
            r5(j,i)= r5(j,i)+f5w(j,6)*work(i+2)
          enddo
        enddo
!
        do i= 1,3
          f6w( 1,i)=(f6(1,i)*xiq6+f5(1,i+3)*xiq4*p15+f4(1,i+6)*xiq2*p45 &
&                   +f3(1,i+6)*p15)
          f6w( 2,i)=(f6(1,i)*xiq4+f5(1,i+3)*xiq2*ten+f4(1,i+6)*p15)*xyiq
          f6w( 3,i)=(f6(2,i)*xiq4+f5(2,i+3)*xiq2*ten+f4(2,i+6)*p15)*xiq
          f6w( 4,i)=(f6(1,i)*x4y2+f5(1,i+3)*xiq4+f5(1,i+3)*xyiq2*six &
&                   +f4(1,i+6)*xiq2*six+f4(1,i+6)*yiq2*three+f3(1,i+6)*three)
          f6w( 5,i)=(f6(2,i)*xiq4+f5(2,i+3)*xiq2*six+f4(2,i+6)*three)*yiq
          f6w( 6,i)=(f6(3,i)*xiq4+f5(1,i+3)*xiq4+f5(3,i+3)*xiq2*six &
&                   +f4(1,i+6)*xiq2*six+f4(3,i+6)*three+f3(1,i+6)*three)
          f6w( 7,i)=(f6(1,i)*xyiq2+f5(1,i+3)*xiq2*three+f5(1,i+3)*yiq2*three &
&                   +f4(1,i+6)*nine)*xyiq
          f6w( 8,i)=(f6(2,i)*xyiq2+f5(2,i+3)*xiq2+f5(2,i+3)*yiq2*three &
&                   +f4(2,i+6)*three)*xiq
          f6w( 9,i)=(f6(3,i)*xiq2+f5(1,i+3)*xiq2+f5(3,i+3)*three &
&                   +f4(1,i+6)*three)*xyiq
          f6w(10,i)=(f6(4,i)*xiq2+f5(2,i+3)*xiq2*three+f5(4,i+3)*three &
&                   +f4(2,i+6)*nine)*xiq
          f6w(11,i)=(f6(1,i)*x2y4+f5(1,i+3)*xyiq2*six+f5(1,i+3)*yiq4 &
&                   +f4(1,i+6)*xiq2*three+f4(1,i+6)*yiq2*six+f3(1,i+6)*three)
          f6w(12,i)=(f6(2,i)*xyiq2+f5(2,i+3)*xiq2*three+f5(2,i+3)*yiq2 &
&                   +f4(2,i+6)*three)*yiq
          f6w(13,i)=(f6(3,i)*xyiq2+f5(1,i+3)*xyiq2+f5(3,i+3)*xiq2+f5(3,i+3)*yiq2 &
&                   +f4(1,i+6)*xiq2+f4(1,i+6)*yiq2+f4(3,i+6)+f3(1,i+6))
          f6w(14,i)=(f6(4,i)*xiq2+f5(2,i+3)*xiq2*three+f5(4,i+3) &
&                   +f4(2,i+6)*three)*yiq
          f6w(15,i)=(f6(5,i)*xiq2+f5(3,i+3)*xiq2*six+f5(5,i+3) &
&                   +f4(1,i+6)*xiq2*three+f4(3,i+6)*six+f3(1,i+6)*three)
          f6w(16,i)=(f6(1,i)*yiq4+f5(1,i+3)*yiq2*ten+f4(1,i+6)*p15)*xyiq
          f6w(17,i)=(f6(2,i)*yiq4+f5(2,i+3)*yiq2*six+f4(2,i+6)*three)*xiq
          f6w(18,i)=(f6(3,i)*yiq2+f5(1,i+3)*yiq2+f5(3,i+3)*three &
&                   +f4(1,i+6)*three)*xyiq
          f6w(19,i)=(f6(4,i)*yiq2+f5(2,i+3)*yiq2*three+f5(4,i+3) &
&                   +f4(2,i+6)*three)*xiq
          f6w(20,i)=(f6(5,i)+f5(3,i+3)*six+f4(1,i+6)*three)*xyiq
          f6w(21,i)=(f6(6,i)+f5(4,i+3)*ten+f4(2,i+6)*p15)*xiq
          f6w(22,i)=(f6(1,i)*yiq6+f5(1,i+3)*yiq4*p15+f4(1,i+6)*yiq2*p45 &
&                   +f3(1,i+6)*p15)
          f6w(23,i)=(f6(2,i)*yiq4+f5(2,i+3)*yiq2*ten+f4(2,i+6)*p15)*yiq
          f6w(24,i)=(f6(3,i)*yiq4+f5(1,i+3)*yiq4+f5(3,i+3)*yiq2*six &
&                   +f4(1,i+6)*yiq2*six+f4(3,i+6)*three+f3(1,i+6)*three)
          f6w(25,i)=(f6(4,i)*yiq2+f5(2,i+3)*yiq2*three+f5(4,i+3)*three &
&                   +f4(2,i+6)*nine)*yiq
          f6w(26,i)=(f6(5,i)*yiq2+f5(3,i+3)*yiq2*six+f5(5,i+3) &
&                   +f4(1,i+6)*yiq2*three+f4(3,i+6)*six+f3(1,i+6)*three)
          f6w(27,i)=(f6(6,i)+f5(4,i+3)*ten+f4(2,i+6)*p15)*yiq
          f6w(28,i)=(f6(7,i)+f5(5,i+3)*p15+f4(3,i+6)*p45+f3(1,i+6)*p15)
        enddo
        do j= 1,28
          r6(j,1)= r6(j,1)+f6w(j,1)*qmd4
          r6(j,2)= r6(j,2)+f6w(j,2)*qmd4
          r6(j,3)= r6(j,3)+f6w(j,3)*work(14)
          r6(j,4)= r6(j,4)+f6w(j,3)*work(15)
        enddo
!
        r7( 1)= r7( 1)+(f7(1)*xiq6+f6(1,3)*xiq4*p21+f5(1,6)*xiq2*p105+f4(1,9)*p105)*qmd4x
        r7( 2)= r7( 2)+(f7(1)*xiq6+f6(1,3)*xiq4*p15+f5(1,6)*xiq2*p45+f4(1,9)*p15)*qmd4y
        r7( 3)= r7( 3)+(f7(2)*xiq6+f6(2,3)*xiq4*p15+f5(2,6)*xiq2*p45+f4(2,9)*p15)*qmd4
        r7( 4)= r7( 4)+(f7(1)*x4y2+f6(1,3)*xiq4+f6(1,3)*xyiq2*ten+f5(1,6)*xiq2*ten &
&                      +f5(1,6)*yiq2*p15+f4(1,9)*p15)*qmd4x
        r7( 5)= r7( 5)+(f7(2)*xiq4+f6(2,3)*xiq2*ten+f5(2,6)*p15)*qmd4xy
        r7( 6)= r7( 6)+(f7(3)*xiq4+f6(1,3)*xiq4+f6(3,3)*xiq2*ten+f5(1,6)*xiq2*ten &
&                      +f5(3,6)*p15+f4(1,9)*p15)*qmd4x
        r7( 7)= r7( 7)+(f7(1)*x4y2+f6(1,3)*xiq4*three+f6(1,3)*xyiq2*six+f5(1,6)*xiq2*p18 &
&                      +f5(1,6)*yiq2*three+f4(1,9)*nine)*qmd4y
        r7( 8)= r7( 8)+(f7(2)*x4y2+f6(2,3)*xiq4+f6(2,3)*xyiq2*six+f5(2,6)*xiq2*six &
&                      +f5(2,6)*yiq2*three+f4(2,9)*three)*qmd4
        r7( 9)= r7( 9)+(f7(3)*xiq4+f6(1,3)*xiq4+f6(3,3)*xiq2*six+f5(1,6)*xiq2*six &
&                      +f5(3,6)*three+f4(1,9)*three)*qmd4y
        r7(10)= r7(10)+(f7(4)*xiq4+f6(2,3)*xiq4*three+f6(4,3)*xiq2*six+f5(2,6)*xiq2*p18 &
&                      +f5(4,6)*three+f4(2,9)*nine)*qmd4
        r7(11)= r7(11)+(f7(1)*x2y4+f6(1,3)*xyiq2*six+f6(1,3)*yiq4*three+f5(1,6)*xiq2*three &
&                      +f5(1,6)*yiq2*p18+f4(1,9)*nine)*qmd4x
        r7(12)= r7(12)+(f7(2)*xyiq2+f6(2,3)*xiq2*three+f6(2,3)*yiq2*three+f5(2,6)*nine)*qmd4xy
        r7(13)= r7(13)+(f7(3)*xyiq2+f6(1,3)*xyiq2+f6(3,3)*xiq2+f6(3,3)*yiq2*three &
&                      +f5(1,6)*xiq2+f5(1,6)*yiq2*three+f5(3,6)*three+f4(1,9)*three)*qmd4x
        r7(14)= r7(14)+(f7(4)*xiq2+f6(2,3)*xiq2*three+f6(4,3)*three+f5(2,6)*nine)*qmd4xy
        r7(15)= r7(15)+(f7(5)*xiq2+f6(3,3)*xiq2*six+f6(5,3)*three+f5(1,6)*xiq2*three &
&                      +f5(3,6)*p18+f4(1,9)*nine)*qmd4x
        r7(16)= r7(16)+(f7(1)*x2y4+f6(1,3)*xyiq2*ten+f6(1,3)*yiq4+f5(1,6)*xiq2*p15 &
&                      +f5(1,6)*yiq2*ten+f4(1,9)*p15)*qmd4y
        r7(17)= r7(17)+(f7(2)*x2y4+f6(2,3)*xyiq2*six+f6(2,3)*yiq4+f5(2,6)*xiq2*three &
&                      +f5(2,6)*yiq2*six+f4(2,9)*three)*qmd4
        r7(18)= r7(18)+(f7(3)*xyiq2+f6(1,3)*xyiq2+f6(3,3)*xiq2*three+f6(3,3)*yiq2 &
&                      +f5(1,6)*xiq2*three+f5(1,6)*yiq2+f5(3,6)*three+f4(1,9)*three)*qmd4y
        r7(19)= r7(19)+(f7(4)*xyiq2+f6(2,3)*xyiq2*three+f6(4,3)*xiq2+f6(4,3)*yiq2 &
&                      +f5(2,6)*xiq2*three+f5(2,6)*yiq2*three+f5(4,6)+f4(2,9)*three)*qmd4
        r7(20)= r7(20)+(f7(5)*xiq2+f6(3,3)*xiq2*six+f6(5,3)+f5(1,6)*xiq2*three+f5(3,6)*six &
&                      +f4(1,9)*three)*qmd4y
        r7(21)= r7(21)+(f7(6)*xiq2+f6(4,3)*xiq2*ten+f6(6,3)+f5(2,6)*xiq2*p15+f5(4,6)*ten &
&                      +f4(2,9)*p15)*qmd4
        r7(22)= r7(22)+(f7(1)*yiq6+f6(1,3)*yiq4*p15+f5(1,6)*yiq2*p45+f4(1,9)*p15)*qmd4x
        r7(23)= r7(23)+(f7(2)*yiq4+f6(2,3)*yiq2*ten+f5(2,6)*p15)*qmd4xy
        r7(24)= r7(24)+(f7(3)*yiq4+f6(1,3)*yiq4+f6(3,3)*yiq2*six+f5(1,6)*yiq2*six &
&                      +f5(3,6)*three+f4(1,9)*three)*qmd4x
        r7(25)= r7(25)+(f7(4)*yiq2+f6(2,3)*yiq2*three+f6(4,3)*three+f5(2,6)*nine)*qmd4xy
        r7(26)= r7(26)+(f7(5)*yiq2+f6(3,3)*yiq2*six+f6(5,3)+f5(1,6)*yiq2*three+f5(3,6)*six &
&                      +f4(1,9)*three)*qmd4x
        r7(27)= r7(27)+(f7(6)+f6(4,3)*ten+f5(2,6)*p15)*qmd4xy
        r7(28)= r7(28)+(f7(7)+f6(5,3)*p15+f5(3,6)*p45+f4(1,9)*p15)*qmd4x
        r7(29)= r7(29)+(f7(1)*yiq6+f6(1,3)*yiq4*p21+f5(1,6)*yiq2*p105+f4(1,9)*p105)*qmd4y
        r7(30)= r7(30)+(f7(2)*yiq6+f6(2,3)*yiq4*p15+f5(2,6)*yiq2*p45+f4(2,9)*p15)*qmd4
        r7(31)= r7(31)+(f7(3)*yiq4+f6(1,3)*yiq4+f6(3,3)*yiq2*ten+f5(1,6)*yiq2*ten &
&                      +f5(3,6)*p15+f4(1,9)*p15)*qmd4y
        r7(32)= r7(32)+(f7(4)*yiq4+f6(2,3)*yiq4*three+f6(4,3)*yiq2*six+f5(2,6)*yiq2*p18 &
&                      +f5(4,6)*three+f4(2,9)*nine)*qmd4
        r7(33)= r7(33)+(f7(5)*yiq2+f6(3,3)*yiq2*six+f6(5,3)*three+f5(1,6)*yiq2*three &
&                      +f5(3,6)*p18+f4(1,9)*nine)*qmd4y
        r7(34)= r7(34)+(f7(6)*yiq2+f6(4,3)*yiq2*ten+f6(6,3)+f5(2,6)*yiq2*p15+f5(4,6)*ten &
&                      +f4(2,9)*p15)*qmd4
        r7(35)= r7(35)+(f7(7)+f6(5,3)*p15+f5(3,6)*p45+f4(1,9)*p15)*qmd4y
        r7(36)= r7(36)+(f7(8)+f6(6,3)*p21+f5(4,6)*p105+f4(2,9)*p105)*qmd4
      enddo
!
      qx= xzkl(1)
      qz= xzkl(2)
      call int2dddp1(eri,r0,r1,r2,r3,r4,r5,r6,r7,qx,qz)
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
      do j= 1,6
        do k= 1,6
          do l= 1,6
            work(1)= eri(l,k,j,1)
            work(2)= eri(l,k,j,2)
            work(3)= eri(l,k,j,3)
            eri(l,k,j,1)= work(1)*rot(1,1)+work(2)*rot(2,1)+work(3)*rot(3,1)
            eri(l,k,j,2)= work(1)*rot(1,2)+work(2)*rot(2,2)+work(3)*rot(3,2)
            eri(l,k,j,3)= work(1)*rot(1,3)+work(2)*rot(2,3)+work(3)*rot(3,3)
          enddo
        enddo
      enddo
!
      if(nbfijkl(2) == 6)then
        do i= 1,3
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
        do i= 1,3
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
        do i= 1,3
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
        do i= 1,3
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
        do i= 1,3
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
        do i= 1,3
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


!----------------------------------------------------------
  subroutine int2dddp1(eri,r0,r1,r2,r3,r4,r5,r6,r7,qx,qz)
!----------------------------------------------------------
!
      implicit none
      integer :: i, j
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00, four=4.0D+00, five=5.0D+00
      real(8),parameter :: six=6.0D+00, seven=7.0D+00, eight=8.0D+00, nine=9.0D+00, ten=1.0D+01
      real(8),parameter :: p11=1.1D+01, p12=1.2D+01, p13=1.3D+01, p15=1.5D+01, p18=1.8D+01
      real(8),parameter :: p21=2.1D+01, p45=4.5D+01, p105=1.05D+2
      real(8),parameter :: sqrt3=1.73205080756888D+00
      real(8),intent(in) :: r0(15), r1(3,27), r2(6,34), r3(10,31), r4(15,21), r5(21,11)
      real(8),intent(in) :: r6(28,4), r7(36), qx, qz
      real(8),intent(out) :: eri(6,6,6,3)
      real(8) :: rrr(15), rxyz(45), xx, xz, zz, xxx, xxz, xzz, zzz
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
      rrr( 1)=-r7(1)-r5(1,3)*three-r5(1,8)*six-r3(1,15)*p18-r3(1,27)*three-r1(1,23)*nine
      rrr( 2)=-r7(2)-r5(2,3)*three-r5(2,8)*three-r3(2,15)*nine
      rrr( 3)=-r7(3)-r5(3,3)*three-r5(3,8)*three-r3(3,15)*nine
      rrr( 4)=-r7(4)-r5(4,3)*three-r5(1,8)-r5(4,8)-r3(1,15)*three-r3(4,15)*three-r3(1,27) &
&          -r1(1,23)*three
      rrr( 5)=-r7(5)-r5(5,3)*three-r5(5,8)-r3(5,15)*three
      rrr( 6)=-r7(6)-r5(6,3)*three-r5(1,8)-r5(6,8)-r3(1,15)*three-r3(6,15)*three-r3(1,27) &
&          -r1(1,23)*three
      rrr( 7)=-r7(7)-r5(7,3)*three-r5(2,8)*three-r3(2,15)*nine
      rrr( 8)=-r7(8)-r5(8,3)*three-r5(3,8)-r3(3,15)*three
      rrr( 9)=-r7(9)-r5(9,3)*three-r5(2,8)-r3(2,15)*three
      rrr(10)=-r7(10)-r5(10,3)*three-r5(3,8)*three-r3(3,15)*nine
      rrr(11)=-r7(11)-r5(11,3)*three-r5(4,8)*six-r3(4,15)*p18-r3(1,27)*three-r1(1,23)*nine
      rrr(12)=-r7(12)-r5(12,3)*three-r5(5,8)*three-r3(5,15)*nine
      rrr(13)=-r7(13)-r5(13,3)*three-r5(4,8)-r5(6,8)-r3(4,15)*three-r3(6,15)*three-r3(1,27) &
&          -r1(1,23)*three
      rrr(14)=-r7(14)-r5(14,3)*three-r5(5,8)*three-r3(5,15)*nine
      rrr(15)=-r7(15)-r5(15,3)*three-r5(6,8)*six-r3(6,15)*p18-r3(1,27)*three-r1(1,23)*nine
      rxyz(1)=-r3(1,31)-r1(1,27)*three
      rxyz(2)=-r4(2,21)-r2(4,24)*three
      rxyz(3)=-r4(2,20)-r2(4,23)*three
      rxyz(4)=-r5(4,9)-r3(4,16)*three-r3(1,28)-r1(1,24)*three
      rxyz(5)=-r5(4,11)-r3(4,18)*three-r3(1,30)-r1(1,26)*three
      rxyz(6)=-r5(4,10)-r3(4,17)*three-r3(1,29)-r1(1,25)*three
      rxyz(7)=-r6(7,3)-r4(7,8)*three-r4(2,18)*three-r2(4,21)*nine
      rxyz(8)=-r6(7,4)-r4(7,9)*three-r4(2,19)*three-r2(4,22)*nine
      rxyz(9)=-r5(5,9)-r3(5,16)*three
      rxyz(10)=-r5(5,10)-r3(5,17)*three
      rxyz(11)=-r5(5,11)-r3(5,18)*three
      rxyz(12)=-r6(8,3)-r4(8,8)*three-r4(3,18)-r2(5,21)*three
      rxyz(13)=-r6(8,4)-r4(8,9)*three-r4(3,19)-r2(5,22)*three
      rxyz(14)=-r6(4,4)-r4(4,9)*three-r4(1,19)-r2(1,22)*three
      rxyz(15)=-r6(4,3)-r4(4,8)*three-r4(1,18)-r2(1,21)*three
      rxyz(16)=-r6(2,4)-r4(2,9)*three-r4(2,19)-r2(4,22)*three
      rxyz(17)=-r6(2,3)-r4(2,8)*three-r4(2,18)-r2(4,21)*three
      rxyz(18)=-r6(9,3)-r4(9,8)*three-r4(2,18)-r2(4,21)*three
      rxyz(19)=-r6(9,4)-r4(9,9)*three-r4(2,19)-r2(4,22)*three
      rxyz(20)=-r5(3,9)-r3(3,16)*three
      rxyz(21)=-r5(3,10)*two-r3(3,17)*six
      rxyz(22)=-r5(3,11)-r3(3,18)*three
      rxyz(23)=-r5(1,9)-r3(1,16)*three-r3(1,28)-r1(1,24)*three
      rxyz(24)=-r5(1,10)-r3(1,17)*three-r3(1,29)-r1(1,25)*three
      rxyz(25)=-r5(1,11)-r3(1,18)*three-r3(1,30)-r1(1,26)*three
      rxyz(26)=-r4(1,20)-r2(1,23)*three
      rxyz(27)=-r4(1,21)-r2(1,24)*three
      rxyz(28)=-r6(6,3)-r4(6,8)*three-r4(1,18)-r2(1,21)*three
      rxyz(29)=-r6(6,4)-r4(6,9)*three-r4(1,19)-r2(1,22)*three
      rxyz(30)=-r6(3,3)-r4(3,8)*three-r4(3,18)-r2(5,21)*three
      rxyz(31)=-r6(3,4)-r4(3,9)*three-r4(3,19)-r2(5,22)*three
      rxyz(32)=-r5(6,9)-r3(6,16)*three-r3(1,28)-r1(1,24)*three
      rxyz(33)=-r5(6,10)-r3(6,17)*three-r3(1,29)-r1(1,25)*three
      rxyz(34)=-r5(6,11)-r3(6,18)*three-r3(1,30)-r1(1,26)*three
      rxyz(35)=-r4(3,20)-r2(5,23)*three
      rxyz(36)=-r4(3,21)-r2(5,24)*three
      rxyz(37)=-r6(1,3)-r4(1,8)*three-r4(1,18)*three-r2(1,21)*nine
      rxyz(38)=-r6(1,4)-r4(1,9)*three-r4(1,19)*three-r2(1,22)*nine
      rxyz(39)=-r6(5,3)-r4(5,8)*three
      rxyz(40)=-r6(5,4)-r4(5,9)*three
      rxyz(41)=-r5(2,9)-r3(2,16)*three
      rxyz(42)=-r5(2,10)-r3(2,17)*three
      rxyz(43)=-r5(2,11)-r3(2,18)*three
      rxyz(44)=-r6(10,3)-r4(10,8)*three-r4(3,18)*three-r2(5,21)*nine
      rxyz(45)=-r6(10,4)-r4(10,9)*three-r4(3,19)*three-r2(5,22)*nine
      call int2dddp2(eri,rrr,rxyz,qx,qz,1,1)
!
      rrr( 1)=-r7(4)-r5(1,3)-r5(4,8)*six-r3(1,15)*six-r3(4,27)*three-r1(1,23)*three
      rrr( 2)=-r7(7)-r5(2,3)-r5(7,8)*three-r3(2,15)*three
      rrr( 3)=-r7(8)-r5(3,3)-r5(8,8)*three-r3(3,15)*three
      rrr( 4)=-r7(11)-r5(4,3)-r5(4,8)-r5(11,8)-r3(1,15)-r3(4,15)-r3(4,27)-r1(1,23)
      rrr( 5)=-r7(12)-r5(5,3)-r5(12,8)-r3(5,15)
      rrr( 6)=-r7(13)-r5(6,3)-r5(4,8)-r5(13,8)-r3(1,15)-r3(6,15)-r3(4,27)-r1(1,23)
      rrr( 7)=-r7(16)-r5(7,3)-r5(7,8)*three-r3(2,15)*three
      rrr( 8)=-r7(17)-r5(8,3)-r5(8,8)-r3(3,15)
      rrr( 9)=-r7(18)-r5(9,3)-r5(7,8)-r3(2,15)
      rrr(10)=-r7(19)-r5(10,3)-r5(8,8)*three-r3(3,15)*three
      rrr(11)=-r7(22)-r5(11,3)-r5(11,8)*six-r3(4,15)*six-r3(4,27)*three-r1(1,23)*three
      rrr(12)=-r7(23)-r5(12,3)-r5(12,8)*three-r3(5,15)*three
      rrr(13)=-r7(24)-r5(13,3)-r5(11,8)-r5(13,8)-r3(4,15)-r3(6,15)-r3(4,27)-r1(1,23)
      rrr(14)=-r7(25)-r5(14,3)-r5(12,8)*three-r3(5,15)*three
      rrr(15)=-r7(26)-r5(15,3)-r5(13,8)*six-r3(6,15)*six-r3(4,27)*three-r1(1,23)*three
      rxyz(1)=-r3(4,31)-r1(1,27)
      rxyz(2)=-r4(7,21)-r2(4,24)
      rxyz(3)=-r4(7,20)-r2(4,23)
      rxyz(4)=-r5(11,9)-r3(4,16)-r3(4,28)-r1(1,24)
      rxyz(5)=-r5(11,11)-r3(4,18)-r3(4,30)-r1(1,26)
      rxyz(6)=-r5(11,10)-r3(4,17)-r3(4,29)-r1(1,25)
      rxyz(7)=-r6(16,3)-r4(7,8)-r4(7,18)*three-r2(4,21)*three
      rxyz(8)=-r6(16,4)-r4(7,9)-r4(7,19)*three-r2(4,22)*three
      rxyz(9)=-r5(12,9)-r3(5,16)
      rxyz(10)=-r5(12,10)-r3(5,17)
      rxyz(11)=-r5(12,11)-r3(5,18)
      rxyz(12)=-r6(17,3)-r4(8,8)-r4(8,18)-r2(5,21)
      rxyz(13)=-r6(17,4)-r4(8,9)-r4(8,19)-r2(5,22)
      rxyz(14)=-r6(11,4)-r4(4,9)-r4(4,19)-r2(1,22)
      rxyz(15)=-r6(11,3)-r4(4,8)-r4(4,18)-r2(1,21)
      rxyz(16)=-r6(7,4)-r4(2,9)-r4(7,19)-r2(4,22)
      rxyz(17)=-r6(7,3)-r4(2,8)-r4(7,18)-r2(4,21)
      rxyz(18)=-r6(18,3)-r4(9,8)-r4(7,18)-r2(4,21)
      rxyz(19)=-r6(18,4)-r4(9,9)-r4(7,19)-r2(4,22)
      rxyz(20)=-r5(8,9)-r3(3,16)
      rxyz(21)=-r5(8,10)*two-r3(3,17)*two
      rxyz(22)=-r5(8,11)-r3(3,18)
      rxyz(23)=-r5(4,9)-r3(1,16)-r3(4,28)-r1(1,24)
      rxyz(24)=-r5(4,10)-r3(1,17)-r3(4,29)-r1(1,25)
      rxyz(25)=-r5(4,11)-r3(1,18)-r3(4,30)-r1(1,26)
      rxyz(26)=-r4(4,20)-r2(1,23)
      rxyz(27)=-r4(4,21)-r2(1,24)
      rxyz(28)=-r6(13,3)-r4(6,8)-r4(4,18)-r2(1,21)
      rxyz(29)=-r6(13,4)-r4(6,9)-r4(4,19)-r2(1,22)
      rxyz(30)=-r6(8,3)-r4(3,8)-r4(8,18)-r2(5,21)
      rxyz(31)=-r6(8,4)-r4(3,9)-r4(8,19)-r2(5,22)
      rxyz(32)=-r5(13,9)-r3(6,16)-r3(4,28)-r1(1,24)
      rxyz(33)=-r5(13,10)-r3(6,17)-r3(4,29)-r1(1,25)
      rxyz(34)=-r5(13,11)-r3(6,18)-r3(4,30)-r1(1,26)
      rxyz(35)=-r4(8,20)-r2(5,23)
      rxyz(36)=-r4(8,21)-r2(5,24)
      rxyz(37)=-r6(4,3)-r4(1,8)-r4(4,18)*three-r2(1,21)*three
      rxyz(38)=-r6(4,4)-r4(1,9)-r4(4,19)*three-r2(1,22)*three
      rxyz(39)=-r6(12,3)-r4(5,8)
      rxyz(40)=-r6(12,4)-r4(5,9)
      rxyz(41)=-r5(7,9)-r3(2,16)
      rxyz(42)=-r5(7,10)-r3(2,17)
      rxyz(43)=-r5(7,11)-r3(2,18)
      rxyz(44)=-r6(19,3)-r4(10,8)-r4(8,18)*three-r2(5,21)*three
      rxyz(45)=-r6(19,4)-r4(10,9)-r4(8,19)*three-r2(5,22)*three
      call int2dddp2(eri,rrr,rxyz,qx,qz,2,1)
!
      rrr( 1)=-r7(6)+r6(3,1)*two-r5(1,1)-r5(1,3)-r5(6,8)*six+r4(3,10)*p12-r3(1,7)*six &
&          -r3(1,15)*six-r3(6,27)*three+r2(5,25)*six-r1(1,13)*three-r1(1,23)*three
      rrr( 2)=-r7(9)+r6(5,1)*two-r5(2,1)-r5(2,3)-r5(9,8)*three+r4(5,10)*six-r3(2,7)*three &
&          -r3(2,15)*three
      rrr( 3)=-r7(10)+r6(6,1)*two-r5(3,1)-r5(3,3)-r5(10,8)*three+r4(6,10)*six-r3(3,7)*three &
&          -r3(3,15)*three
      rrr( 4)=-r7(13)+r6(8,1)*two-r5(4,1)-r5(4,3)-r5(6,8)-r5(13,8)+r4(3,10)*two+r4(8,10)*two &
&          -r3(1,7)-r3(4,7)-r3(1,15)-r3(4,15)-r3(6,27)+r2(5,25)*two-r1(1,13)-r1(1,23)
      rrr( 5)=-r7(14)+r6(9,1)*two-r5(5,1)-r5(5,3)-r5(14,8)+r4(9,10)*two-r3(5,7)-r3(5,15)
      rrr( 6)=-r7(15)+r6(10,1)*two-r5(6,1)-r5(6,3)-r5(6,8)-r5(15,8)+r4(3,10)*two &
&          +r4(10,10)*two-r3(1,7)-r3(6,7)-r3(1,15)-r3(6,15)-r3(6,27)+r2(5,25)*two-r1(1,13) &
&          -r1(1,23)
      rrr( 7)=-r7(18)+r6(12,1)*two-r5(7,1)-r5(7,3)-r5(9,8)*three+r4(5,10)*six-r3(2,7)*three &
&          -r3(2,15)*three
      rrr( 8)=-r7(19)+r6(13,1)*two-r5(8,1)-r5(8,3)-r5(10,8)+r4(6,10)*two-r3(3,7)-r3(3,15)
      rrr( 9)=-r7(20)+r6(14,1)*two-r5(9,1)-r5(9,3)-r5(9,8)+r4(5,10)*two-r3(2,7)-r3(2,15)
      rrr(10)=-r7(21)+r6(15,1)*two-r5(10,1)-r5(10,3)-r5(10,8)*three+r4(6,10)*six &
&          -r3(3,7)*three-r3(3,15)*three
      rrr(11)=-r7(24)+r6(17,1)*two-r5(11,1)-r5(11,3)-r5(13,8)*six+r4(8,10)*p12-r3(4,7)*six &
&          -r3(4,15)*six-r3(6,27)*three+r2(5,25)*six-r1(1,13)*three-r1(1,23)*three
      rrr(12)=-r7(25)+r6(18,1)*two-r5(12,1)-r5(12,3)-r5(14,8)*three+r4(9,10)*six &
&          -r3(5,7)*three-r3(5,15)*three
      rrr(13)=-r7(26)+r6(19,1)*two-r5(13,1)-r5(13,3)-r5(13,8)-r5(15,8)+r4(8,10)*two &
&          +r4(10,10)*two-r3(4,7)-r3(6,7)-r3(4,15)-r3(6,15)-r3(6,27)+r2(5,25)*two-r1(1,13) &
&          -r1(1,23)
      rrr(14)=-r7(27)+r6(20,1)*two-r5(14,1)-r5(14,3)-r5(14,8)*three+r4(9,10)*six &
&          -r3(5,7)*three-r3(5,15)*three
      rrr(15)=-r7(28)+r6(21,1)*two-r5(15,1)-r5(15,3)-r5(15,8)*six+r4(10,10)*p12-r3(6,7)*six &
&          -r3(6,15)*six-r3(6,27)*three+r2(5,25)*six-r1(1,13)*three-r1(1,23)*three
      rxyz(1)=-r3(6,31)+r2(5,29)*two-r1(1,17)-r1(1,27)
      rxyz(2)=-r4(9,21)+r3(5,22)*two-r2(4,16)-r2(4,24)
      rxyz(3)=-r4(9,20)+r3(5,21)*two-r2(4,15)-r2(4,23)
      rxyz(4)=-r5(13,9)+r4(8,11)*two-r3(4,8)-r3(4,16)-r3(6,28)+r2(5,26)*two-r1(1,14) &
&             -r1(1,24)
      rxyz(5)=-r5(13,11)+r4(8,13)*two-r3(4,10)-r3(4,18)-r3(6,30)+r2(5,28)*two-r1(1,16) &
&             -r1(1,26)
      rxyz(6)=-r5(13,10)+r4(8,12)*two-r3(4,9)-r3(4,17)-r3(6,29)+r2(5,27)*two-r1(1,15) &
&             -r1(1,25)
      rxyz(7)=-r6(18,3)+r5(12,4)*two-r4(7,4)-r4(7,8)-r4(9,18)*three+r3(5,19)*six &
&             -r2(4,13)*three-r2(4,21)*three
      rxyz(8)=-r6(18,4)+r5(12,5)*two-r4(7,5)-r4(7,9)-r4(9,19)*three+r3(5,20)*six &
&             -r2(4,14)*three-r2(4,22)*three
      rxyz(9)=-r5(14,9)+r4(9,11)*two-r3(5,8)-r3(5,16)
      rxyz(10)=-r5(14,10)+r4(9,12)*two-r3(5,9)-r3(5,17)
      rxyz(11)=-r5(14,11)+r4(9,13)*two-r3(5,10)-r3(5,18)
      rxyz(12)=-r6(19,3)+r5(13,4)*two-r4(8,4)-r4(8,8)-r4(10,18)+r3(6,19)*two-r2(5,13) &
&             -r2(5,21)
      rxyz(13)=-r6(19,4)+r5(13,5)*two-r4(8,5)-r4(8,9)-r4(10,19)+r3(6,20)*two-r2(5,14) &
&             -r2(5,22)
      rxyz(14)=-r6(13,4)+r5(8,5)*two-r4(4,5)-r4(4,9)-r4(6,19)+r3(3,20)*two-r2(1,14) &
&             -r2(1,22)
      rxyz(15)=-r6(13,3)+r5(8,4)*two-r4(4,4)-r4(4,8)-r4(6,18)+r3(3,19)*two-r2(1,13) &
&             -r2(1,21)
      rxyz(16)=-r6(9,4)+r5(5,5)*two-r4(2,5)-r4(2,9)-r4(9,19)+r3(5,20)*two-r2(4,14) &
&             -r2(4,22)
      rxyz(17)=-r6(9,3)+r5(5,4)*two-r4(2,4)-r4(2,8)-r4(9,18)+r3(5,19)*two-r2(4,13) &
&             -r2(4,21)
      rxyz(18)=-r6(20,3)+r5(14,4)*two-r4(9,4)-r4(9,8)-r4(9,18)+r3(5,19)*two-r2(4,13) &
&             -r2(4,21)
      rxyz(19)=-r6(20,4)+r5(14,5)*two-r4(9,5)-r4(9,9)-r4(9,19)+r3(5,20)*two-r2(4,14) &
&             -r2(4,22)
      rxyz(20)=-r5(10,9)+r4(6,11)*two-r3(3,8)-r3(3,16)
      rxyz(21)=-r5(10,10)*two+r4(6,12)*four-r3(3,9)*two-r3(3,17)*two
      rxyz(22)=-r5(10,11)+r4(6,13)*two-r3(3,10)-r3(3,18)
      rxyz(23)=-r5(6,9)+r4(3,11)*two-r3(1,8)-r3(1,16)-r3(6,28)+r2(5,26)*two-r1(1,14)-r1(1,24)
      rxyz(24)=-r5(6,10)+r4(3,12)*two-r3(1,9)-r3(1,17)-r3(6,29)+r2(5,27)*two-r1(1,15)-r1(1,25)
      rxyz(25)=-r5(6,11)+r4(3,13)*two-r3(1,10)-r3(1,18)-r3(6,30)+r2(5,28)*two-r1(1,16)-r1(1,26)
      rxyz(26)=-r4(6,20)+r3(3,21)*two-r2(1,15)-r2(1,23)
      rxyz(27)=-r4(6,21)+r3(3,22)*two-r2(1,16)-r2(1,24)
      rxyz(28)=-r6(15,3)+r5(10,4)*two-r4(6,4)-r4(6,8)-r4(6,18)+r3(3,19)*two-r2(1,13)-r2(1,21)
      rxyz(29)=-r6(15,4)+r5(10,5)*two-r4(6,5)-r4(6,9)-r4(6,19)+r3(3,20)*two-r2(1,14)-r2(1,22)
      rxyz(30)=-r6(10,3)+r5(6,4)*two-r4(3,4)-r4(3,8)-r4(10,18)+r3(6,19)*two-r2(5,13)-r2(5,21)
      rxyz(31)=-r6(10,4)+r5(6,5)*two-r4(3,5)-r4(3,9)-r4(10,19)+r3(6,20)*two-r2(5,14)-r2(5,22)
      rxyz(32)=-r5(15,9)+r4(10,11)*two-r3(6,8)-r3(6,16)-r3(6,28)+r2(5,26)*two-r1(1,14)-r1(1,24)
      rxyz(33)=-r5(15,10)+r4(10,12)*two-r3(6,9)-r3(6,17)-r3(6,29)+r2(5,27)*two-r1(1,15)-r1(1,25)
      rxyz(34)=-r5(15,11)+r4(10,13)*two-r3(6,10)-r3(6,18)-r3(6,30)+r2(5,28)*two-r1(1,16)-r1(1,26)
      rxyz(35)=-r4(10,20)+r3(6,21)*two-r2(5,15)-r2(5,23)
      rxyz(36)=-r4(10,21)+r3(6,22)*two-r2(5,16)-r2(5,24)
      rxyz(37)=-r6(6,3)+r5(3,4)*two-r4(1,4)-r4(1,8)-r4(6,18)*three+r3(3,19)*six-r2(1,13)*three &
&              -r2(1,21)*three
      rxyz(38)=-r6(6,4)+r5(3,5)*two-r4(1,5)-r4(1,9)-r4(6,19)*three+r3(3,20)*six-r2(1,14)*three &
&              -r2(1,22)*three
      rxyz(39)=-r6(14,3)+r5(9,4)*two-r4(5,4)-r4(5,8)
      rxyz(40)=-r6(14,4)+r5(9,5)*two-r4(5,5)-r4(5,9)
      rxyz(41)=-r5(9,9)+r4(5,11)*two-r3(2,8)-r3(2,16)
      rxyz(42)=-r5(9,10)+r4(5,12)*two-r3(2,9)-r3(2,17)
      rxyz(43)=-r5(9,11)+r4(5,13)*two-r3(2,10)-r3(2,18)
      rxyz(44)=-r6(21,3)+r5(15,4)*two-r4(10,4)-r4(10,8)-r4(10,18)*three+r3(6,19)*six &
&              -r2(5,13)*three-r2(5,21)*three
      rxyz(45)=-r6(21,4)+r5(15,5)*two-r4(10,5)-r4(10,9)-r4(10,19)*three+r3(6,20)*six &
&              -r2(5,14)*three-r2(5,22)*three
      call int2dddp2(eri,rrr,rxyz,qx,qz,3,1)
!
      rrr( 1)=-r7(2)-r5(2,3)-r5(2,8)*six-r3(2,15)*six-r3(2,27)*three-r1(2,23)*three
      rrr( 2)=-r7(4)-r5(4,3)-r5(4,8)*three-r3(4,15)*three
      rrr( 3)=-r7(5)-r5(5,3)-r5(5,8)*three-r3(5,15)*three
      rrr( 4)=-r7(7)-r5(7,3)-r5(2,8)-r5(7,8)-r3(2,15)-r3(7,15)-r3(2,27)-r1(2,23)
      rrr( 5)=-r7(8)-r5(8,3)-r5(8,8)-r3(8,15)
      rrr( 6)=-r7(9)-r5(9,3)-r5(2,8)-r5(9,8)-r3(2,15)-r3(9,15)-r3(2,27)-r1(2,23)
      rrr( 7)=-r7(11)-r5(11,3)-r5(4,8)*three-r3(4,15)*three
      rrr( 8)=-r7(12)-r5(12,3)-r5(5,8)-r3(5,15)
      rrr( 9)=-r7(13)-r5(13,3)-r5(4,8)-r3(4,15)
      rrr(10)=-r7(14)-r5(14,3)-r5(5,8)*three-r3(5,15)*three
      rrr(11)=-r7(16)-r5(16,3)-r5(7,8)*six-r3(7,15)*six-r3(2,27)*three-r1(2,23)*three
      rrr(12)=-r7(17)-r5(17,3)-r5(8,8)*three-r3(8,15)*three
      rrr(13)=-r7(18)-r5(18,3)-r5(7,8)-r5(9,8)-r3(7,15)-r3(9,15)-r3(2,27)-r1(2,23)
      rrr(14)=-r7(19)-r5(19,3)-r5(8,8)*three-r3(8,15)*three
      rrr(15)=-r7(20)-r5(20,3)-r5(9,8)*six-r3(9,15)*six-r3(2,27)*three-r1(2,23)*three
      rxyz(1)=-r3(2,31)-r1(2,27)
      rxyz(2)=-r4(4,21)-r2(2,24)
      rxyz(3)=-r4(4,20)-r2(2,23)
      rxyz(4)=-r5(7,9)-r3(7,16)-r3(2,28)-r1(2,24)
      rxyz(5)=-r5(7,11)-r3(7,18)-r3(2,30)-r1(2,26)
      rxyz(6)=-r5(7,10)-r3(7,17)-r3(2,29)-r1(2,25)
      rxyz(7)=-r6(11,3)-r4(11,8)-r4(4,18)*three-r2(2,21)*three
      rxyz(8)=-r6(11,4)-r4(11,9)-r4(4,19)*three-r2(2,22)*three
      rxyz(9)=-r5(8,9)-r3(8,16)
      rxyz(10)=-r5(8,10)-r3(8,17)
      rxyz(11)=-r5(8,11)-r3(8,18)
      rxyz(12)=-r6(12,3)-r4(12,8)-r4(5,18)-r2(6,21)
      rxyz(13)=-r6(12,4)-r4(12,9)-r4(5,19)-r2(6,22)
      rxyz(14)=-r6(7,4)-r4(7,9)-r4(2,19)-r2(4,22)
      rxyz(15)=-r6(7,3)-r4(7,8)-r4(2,18)-r2(4,21)
      rxyz(16)=-r6(4,4)-r4(4,9)-r4(4,19)-r2(2,22)
      rxyz(17)=-r6(4,3)-r4(4,8)-r4(4,18)-r2(2,21)
      rxyz(18)=-r6(13,3)-r4(13,8)-r4(4,18)-r2(2,21)
      rxyz(19)=-r6(13,4)-r4(13,9)-r4(4,19)-r2(2,22)
      rxyz(20)=-r5(5,9)-r3(5,16)
      rxyz(21)=-r5(5,10)*two-r3(5,17)*two
      rxyz(22)=-r5(5,11)-r3(5,18)
      rxyz(23)=-r5(2,9)-r3(2,16)-r3(2,28)-r1(2,24)
      rxyz(24)=-r5(2,10)-r3(2,17)-r3(2,29)-r1(2,25)
      rxyz(25)=-r5(2,11)-r3(2,18)-r3(2,30)-r1(2,26)
      rxyz(26)=-r4(2,20)-r2(4,23)
      rxyz(27)=-r4(2,21)-r2(4,24)
      rxyz(28)=-r6(9,3)-r4(9,8)-r4(2,18)-r2(4,21)
      rxyz(29)=-r6(9,4)-r4(9,9)-r4(2,19)-r2(4,22)
      rxyz(30)=-r6(5,3)-r4(5,8)-r4(5,18)-r2(6,21)
      rxyz(31)=-r6(5,4)-r4(5,9)-r4(5,19)-r2(6,22)
      rxyz(32)=-r5(9,9)-r3(9,16)-r3(2,28)-r1(2,24)
      rxyz(33)=-r5(9,10)-r3(9,17)-r3(2,29)-r1(2,25)
      rxyz(34)=-r5(9,11)-r3(9,18)-r3(2,30)-r1(2,26)
      rxyz(35)=-r4(5,20)-r2(6,23)
      rxyz(36)=-r4(5,21)-r2(6,24)
      rxyz(37)=-r6(2,3)-r4(2,8)-r4(2,18)*three-r2(4,21)*three
      rxyz(38)=-r6(2,4)-r4(2,9)-r4(2,19)*three-r2(4,22)*three
      rxyz(39)=-r6(8,3)-r4(8,8)
      rxyz(40)=-r6(8,4)-r4(8,9)
      rxyz(41)=-r5(4,9)-r3(4,16)
      rxyz(42)=-r5(4,10)-r3(4,17)
      rxyz(43)=-r5(4,11)-r3(4,18)
      rxyz(44)=-r6(14,3)-r4(14,8)-r4(5,18)*three-r2(6,21)*three
      rxyz(45)=-r6(14,4)-r4(14,9)-r4(5,19)*three-r2(6,22)*three
      call int2dddp2(eri,rrr,rxyz,qx,qz,4,1)
!
      rrr( 1)=-r7(3)+r6(1,1)-r5(3,3)-r5(3,8)*six+r4(1,2)+r4(1,10)*six-r3(3,15)*six &
&          -r3(3,27)*three+r2(1,5)*six+r2(1,25)*three-r1(3,23)*three+r0(6)*three
      rrr( 2)=-r7(5)+r6(2,1)-r5(5,3)-r5(5,8)*three+r4(2,2)+r4(2,10)*three-r3(5,15)*three &
&          +r2(4,5)*three
      rrr( 3)=-r7(6)+r6(3,1)-r5(6,3)-r5(6,8)*three+r4(3,2)+r4(3,10)*three-r3(6,15)*three &
&          +r2(5,5)*three
      rrr( 4)=-r7(8)+r6(4,1)-r5(8,3)-r5(3,8)-r5(8,8)+r4(4,2)+r4(1,10)+r4(4,10)-r3(3,15) &
&          -r3(8,15)-r3(3,27)+r2(1,5)+r2(2,5)+r2(1,25)-r1(3,23)+r0(6)
      rrr( 5)=-r7(9)+r6(5,1)-r5(9,3)-r5(9,8)+r4(5,2)+r4(5,10)-r3(9,15)+r2(6,5)
      rrr( 6)=-r7(10)+r6(6,1)-r5(10,3)-r5(3,8)-r5(10,8)+r4(6,2)+r4(1,10)+r4(6,10)-r3(3,15) &
&          -r3(10,15)-r3(3,27)+r2(1,5)+r2(3,5)+r2(1,25)-r1(3,23)+r0(6)
      rrr( 7)=-r7(12)+r6(7,1)-r5(12,3)-r5(5,8)*three+r4(7,2)+r4(2,10)*three-r3(5,15)*three &
&          +r2(4,5)*three
      rrr( 8)=-r7(13)+r6(8,1)-r5(13,3)-r5(6,8)+r4(8,2)+r4(3,10)-r3(6,15)+r2(5,5)
      rrr( 9)=-r7(14)+r6(9,1)-r5(14,3)-r5(5,8)+r4(9,2)+r4(2,10)-r3(5,15)+r2(4,5)
      rrr(10)=-r7(15)+r6(10,1)-r5(15,3)-r5(6,8)*three+r4(10,2)+r4(3,10)*three-r3(6,15)*three &
&          +r2(5,5)*three
      rrr(11)=-r7(17)+r6(11,1)-r5(17,3)-r5(8,8)*six+r4(11,2)+r4(4,10)*six-r3(8,15)*six &
&          -r3(3,27)*three+r2(2,5)*six+r2(1,25)*three-r1(3,23)*three+r0(6)*three
      rrr(12)=-r7(18)+r6(12,1)-r5(18,3)-r5(9,8)*three+r4(12,2)+r4(5,10)*three-r3(9,15)*three &
&          +r2(6,5)*three
      rrr(13)=-r7(19)+r6(13,1)-r5(19,3)-r5(8,8)-r5(10,8)+r4(13,2)+r4(4,10)+r4(6,10)-r3(8,15) &
&          -r3(10,15)-r3(3,27)+r2(2,5)+r2(3,5)+r2(1,25)-r1(3,23)+r0(6)
      rrr(14)=-r7(20)+r6(14,1)-r5(20,3)-r5(9,8)*three+r4(14,2)+r4(5,10)*three-r3(9,15)*three &
&          +r2(6,5)*three
      rrr(15)=-r7(21)+r6(15,1)-r5(21,3)-r5(10,8)*six+r4(15,2)+r4(6,10)*six-r3(10,15)*six &
&          -r3(3,27)*three+r2(3,5)*six+r2(1,25)*three-r1(3,23)*three+r0(6)*three
      rxyz(1)=-r3(3,31)+r2(1,29)-r1(3,27)+r0(10)
      rxyz(2)=-r4(5,21)+r3(2,22)-r2(6,24)+r1(2,8)
      rxyz(3)=-r4(5,20)+r3(2,21)-r2(6,23)+r1(2,7)
      rxyz(4)=-r5(8,9)+r4(4,11)-r3(8,16)-r3(3,28)+r2(2,6)+r2(1,26)-r1(3,24)+r0(7)
      rxyz(5)=-r5(8,11)+r4(4,13)-r3(8,18)-r3(3,30)+r2(2,8)+r2(1,28)-r1(3,26)+r0(9)
      rxyz(6)=-r5(8,10)+r4(4,12)-r3(8,17)-r3(3,29)+r2(2,7)+r2(1,27)-r1(3,25)+r0(8)
      rxyz(7)=-r6(12,3)+r5(7,4)-r4(12,8)-r4(5,18)*three+r3(7,3)+r3(2,19)*three &
&             -r2(6,21)*three+r1(2,5)*three
      rxyz(8)=-r6(12,4)+r5(7,5)-r4(12,9)-r4(5,19)*three+r3(7,4)+r3(2,20)*three &
&             -r2(6,22)*three+r1(2,6)*three
      rxyz(9)=-r5(9,9)+r4(5,11)-r3(9,16)+r2(6,6)
      rxyz(10)=-r5(9,10)+r4(5,12)-r3(9,17)+r2(6,7)
      rxyz(11)=-r5(9,11)+r4(5,13)-r3(9,18)+r2(6,8)
      rxyz(12)=-r6(13,3)+r5(8,4)-r4(13,8)-r4(6,18)+r3(8,3)+r3(3,19)-r2(3,21)+r1(3,5)
      rxyz(13)=-r6(13,4)+r5(8,5)-r4(13,9)-r4(6,19)+r3(8,4)+r3(3,20)-r2(3,22)+r1(3,6)
      rxyz(14)=-r6(8,4)+r5(4,5)-r4(8,9)-r4(3,19)+r3(4,4)+r3(1,20)-r2(5,22)+r1(1,6)
      rxyz(15)=-r6(8,3)+r5(4,4)-r4(8,8)-r4(3,18)+r3(4,3)+r3(1,19)-r2(5,21)+r1(1,5)
      rxyz(16)=-r6(5,4)+r5(2,5)-r4(5,9)-r4(5,19)+r3(2,4)+r3(2,20)-r2(6,22)+r1(2,6)
      rxyz(17)=-r6(5,3)+r5(2,4)-r4(5,8)-r4(5,18)+r3(2,3)+r3(2,19)-r2(6,21)+r1(2,5)
      rxyz(18)=-r6(14,3)+r5(9,4)-r4(14,8)-r4(5,18)+r3(9,3)+r3(2,19)-r2(6,21)+r1(2,5)
      rxyz(19)=-r6(14,4)+r5(9,5)-r4(14,9)-r4(5,19)+r3(9,4)+r3(2,20)-r2(6,22)+r1(2,6)
      rxyz(20)=-r5(6,9)+r4(3,11)-r3(6,16)+r2(5,6)
      rxyz(21)=-r5(6,10)*two+r4(3,12)*two-r3(6,17)*two+r2(5,7)*two
      rxyz(22)=-r5(6,11)+r4(3,13)-r3(6,18)+r2(5,8)
      rxyz(23)=-r5(3,9)+r4(1,11)-r3(3,16)-r3(3,28)+r2(1,6)+r2(1,26)-r1(3,24)+r0(7)
      rxyz(24)=-r5(3,10)+r4(1,12)-r3(3,17)-r3(3,29)+r2(1,7)+r2(1,27)-r1(3,25)+r0(8)
      rxyz(25)=-r5(3,11)+r4(1,13)-r3(3,18)-r3(3,30)+r2(1,8)+r2(1,28)-r1(3,26)+r0(9)
      rxyz(26)=-r4(3,20)+r3(1,21)-r2(5,23)+r1(1,7)
      rxyz(27)=-r4(3,21)+r3(1,22)-r2(5,24)+r1(1,8)
      rxyz(28)=-r6(10,3)+r5(6,4)-r4(10,8)-r4(3,18)+r3(6,3)+r3(1,19)-r2(5,21)+r1(1,5)
      rxyz(29)=-r6(10,4)+r5(6,5)-r4(10,9)-r4(3,19)+r3(6,4)+r3(1,20)-r2(5,22)+r1(1,6)
      rxyz(30)=-r6(6,3)+r5(3,4)-r4(6,8)-r4(6,18)+r3(3,3)+r3(3,19)-r2(3,21)+r1(3,5)
      rxyz(31)=-r6(6,4)+r5(3,5)-r4(6,9)-r4(6,19)+r3(3,4)+r3(3,20)-r2(3,22)+r1(3,6)
      rxyz(32)=-r5(10,9)+r4(6,11)-r3(10,16)-r3(3,28)+r2(3,6)+r2(1,26)-r1(3,24)+r0(7)
      rxyz(33)=-r5(10,10)+r4(6,12)-r3(10,17)-r3(3,29)+r2(3,7)+r2(1,27)-r1(3,25)+r0(8)
      rxyz(34)=-r5(10,11)+r4(6,13)-r3(10,18)-r3(3,30)+r2(3,8)+r2(1,28)-r1(3,26)+r0(9)
      rxyz(35)=-r4(6,20)+r3(3,21)-r2(3,23)+r1(3,7)
      rxyz(36)=-r4(6,21)+r3(3,22)-r2(3,24)+r1(3,8)
      rxyz(37)=-r6(3,3)+r5(1,4)-r4(3,8)-r4(3,18)*three+r3(1,3)+r3(1,19)*three-r2(5,21)*three &
&              +r1(1,5)*three
      rxyz(38)=-r6(3,4)+r5(1,5)-r4(3,9)-r4(3,19)*three+r3(1,4)+r3(1,20)*three-r2(5,22)*three &
&              +r1(1,6)*three
      rxyz(39)=-r6(9,3)+r5(5,4)-r4(9,8)+r3(5,3)
      rxyz(40)=-r6(9,4)+r5(5,5)-r4(9,9)+r3(5,4)
      rxyz(41)=-r5(5,9)+r4(2,11)-r3(5,16)+r2(4,6)
      rxyz(42)=-r5(5,10)+r4(2,12)-r3(5,17)+r2(4,7)
      rxyz(43)=-r5(5,11)+r4(2,13)-r3(5,18)+r2(4,8)
      rxyz(44)=-r6(15,3)+r5(10,4)-r4(15,8)-r4(6,18)*three+r3(10,3)+r3(3,19)*three-r2(3,21)*three &
&              +r1(3,5)*three
      rxyz(45)=-r6(15,4)+r5(10,5)-r4(15,9)-r4(6,19)*three+r3(10,4)+r3(3,20)*three-r2(3,22)*three &
&              +r1(3,6)*three
      call int2dddp2(eri,rrr,rxyz,qx,qz,5,1)
!
      rrr( 1)=-r7(5)+r6(2,1)-r5(5,8)*six+r4(2,10)*six-r3(5,27)*three+r2(4,25)*three
      rrr( 2)=-r7(8)+r6(4,1)-r5(8,8)*three+r4(4,10)*three
      rrr( 3)=-r7(9)+r6(5,1)-r5(9,8)*three+r4(5,10)*three
      rrr( 4)=-r7(12)+r6(7,1)-r5(5,8)-r5(12,8)+r4(2,10)+r4(7,10)-r3(5,27)+r2(4,25)
      rrr( 5)=-r7(13)+r6(8,1)-r5(13,8)+r4(8,10)
      rrr( 6)=-r7(14)+r6(9,1)-r5(5,8)-r5(14,8)+r4(2,10)+r4(9,10)-r3(5,27)+r2(4,25)
      rrr( 7)=-r7(17)+r6(11,1)-r5(8,8)*three+r4(4,10)*three
      rrr( 8)=-r7(18)+r6(12,1)-r5(9,8)+r4(5,10)
      rrr( 9)=-r7(19)+r6(13,1)-r5(8,8)+r4(4,10)
      rrr(10)=-r7(20)+r6(14,1)-r5(9,8)*three+r4(5,10)*three
      rrr(11)=-r7(23)+r6(16,1)-r5(12,8)*six+r4(7,10)*six-r3(5,27)*three+r2(4,25)*three
      rrr(12)=-r7(24)+r6(17,1)-r5(13,8)*three+r4(8,10)*three
      rrr(13)=-r7(25)+r6(18,1)-r5(12,8)-r5(14,8)+r4(7,10)+r4(9,10)-r3(5,27)+r2(4,25)
      rrr(14)=-r7(26)+r6(19,1)-r5(13,8)*three+r4(8,10)*three
      rrr(15)=-r7(27)+r6(20,1)-r5(14,8)*six+r4(9,10)*six-r3(5,27)*three+r2(4,25)*three
      rxyz(1)=-r3(5,31)+r2(4,29)
      rxyz(2)=-r4(8,21)+r3(4,22)
      rxyz(3)=-r4(8,20)+r3(4,21)
      rxyz(4)=-r5(12,9)+r4(7,11)-r3(5,28)+r2(4,26)
      rxyz(5)=-r5(12,11)+r4(7,13)-r3(5,30)+r2(4,28)
      rxyz(6)=-r5(12,10)+r4(7,12)-r3(5,29)+r2(4,27)
      rxyz(7)=-r6(17,3)+r5(11,4)-r4(8,18)*three+r3(4,19)*three
      rxyz(8)=-r6(17,4)+r5(11,5)-r4(8,19)*three+r3(4,20)*three
      rxyz(9)=-r5(13,9)+r4(8,11)
      rxyz(10)=-r5(13,10)+r4(8,12)
      rxyz(11)=-r5(13,11)+r4(8,13)
      rxyz(12)=-r6(18,3)+r5(12,4)-r4(9,18)+r3(5,19)
      rxyz(13)=-r6(18,4)+r5(12,5)-r4(9,19)+r3(5,20)
      rxyz(14)=-r6(12,4)+r5(7,5)-r4(5,19)+r3(2,20)
      rxyz(15)=-r6(12,3)+r5(7,4)-r4(5,18)+r3(2,19)
      rxyz(16)=-r6(8,4)+r5(4,5)-r4(8,19)+r3(4,20)
      rxyz(17)=-r6(8,3)+r5(4,4)-r4(8,18)+r3(4,19)
      rxyz(18)=-r6(19,3)+r5(13,4)-r4(8,18)+r3(4,19)
      rxyz(19)=-r6(19,4)+r5(13,5)-r4(8,19)+r3(4,20)
      rxyz(20)=-r5(9,9)+r4(5,11)
      rxyz(21)=-r5(9,10)*two+r4(5,12)*two
      rxyz(22)=-r5(9,11)+r4(5,13)
      rxyz(23)=-r5(5,9)+r4(2,11)-r3(5,28)+r2(4,26)
      rxyz(24)=-r5(5,10)+r4(2,12)-r3(5,29)+r2(4,27)
      rxyz(25)=-r5(5,11)+r4(2,13)-r3(5,30)+r2(4,28)
      rxyz(26)=-r4(5,20)+r3(2,21)
      rxyz(27)=-r4(5,21)+r3(2,22)
      rxyz(28)=-r6(14,3)+r5(9,4)-r4(5,18)+r3(2,19)
      rxyz(29)=-r6(14,4)+r5(9,5)-r4(5,19)+r3(2,20)
      rxyz(30)=-r6(9,3)+r5(5,4)-r4(9,18)+r3(5,19)
      rxyz(31)=-r6(9,4)+r5(5,5)-r4(9,19)+r3(5,20)
      rxyz(32)=-r5(14,9)+r4(9,11)-r3(5,28)+r2(4,26)
      rxyz(33)=-r5(14,10)+r4(9,12)-r3(5,29)+r2(4,27)
      rxyz(34)=-r5(14,11)+r4(9,13)-r3(5,30)+r2(4,28)
      rxyz(35)=-r4(9,20)+r3(5,21)
      rxyz(36)=-r4(9,21)+r3(5,22)
      rxyz(37)=-r6(5,3)+r5(2,4)-r4(5,18)*three+r3(2,19)*three
      rxyz(38)=-r6(5,4)+r5(2,5)-r4(5,19)*three+r3(2,20)*three
      rxyz(39)=-r6(13,3)+r5(8,4)
      rxyz(40)=-r6(13,4)+r5(8,5)
      rxyz(41)=-r5(8,9)+r4(4,11)
      rxyz(42)=-r5(8,10)+r4(4,12)
      rxyz(43)=-r5(8,11)+r4(4,13)
      rxyz(44)=-r6(20,3)+r5(14,4)-r4(9,18)*three+r3(5,19)*three
      rxyz(45)=-r6(20,4)+r5(14,5)-r4(9,19)*three+r3(5,20)*three
      call int2dddp2(eri,rrr,rxyz,qx,qz,6,1)
!
      do i= 1,6
        do j= 1,6
          eri(j,i,1,2)= eri(j,i,4,1)
        enddo
      enddo
!
      rrr( 1)=-r7(7)-r5(2,3)*three-r5(7,8)*six-r3(2,15)*p18-r3(7,27)*three-r1(2,23)*nine
      rrr( 2)=-r7(11)-r5(4,3)*three-r5(11,8)*three-r3(4,15)*nine
      rrr( 3)=-r7(12)-r5(5,3)*three-r5(12,8)*three-r3(5,15)*nine
      rrr( 4)=-r7(16)-r5(7,3)*three-r5(7,8)-r5(16,8)-r3(2,15)*three-r3(7,15)*three-r3(7,27) &
&          -r1(2,23)*three
      rrr( 5)=-r7(17)-r5(8,3)*three-r5(17,8)-r3(8,15)*three
      rrr( 6)=-r7(18)-r5(9,3)*three-r5(7,8)-r5(18,8)-r3(2,15)*three-r3(9,15)*three-r3(7,27) &
&          -r1(2,23)*three
      rrr( 7)=-r7(22)-r5(11,3)*three-r5(11,8)*three-r3(4,15)*nine
      rrr( 8)=-r7(23)-r5(12,3)*three-r5(12,8)-r3(5,15)*three
      rrr( 9)=-r7(24)-r5(13,3)*three-r5(11,8)-r3(4,15)*three
      rrr(10)=-r7(25)-r5(14,3)*three-r5(12,8)*three-r3(5,15)*nine
      rrr(11)=-r7(29)-r5(16,3)*three-r5(16,8)*six-r3(7,15)*p18-r3(7,27)*three-r1(2,23)*nine
      rrr(12)=-r7(30)-r5(17,3)*three-r5(17,8)*three-r3(8,15)*nine
      rrr(13)=-r7(31)-r5(18,3)*three-r5(16,8)-r5(18,8)-r3(7,15)*three-r3(9,15)*three-r3(7,27) &
&          -r1(2,23)*three
      rrr(14)=-r7(32)-r5(19,3)*three-r5(17,8)*three-r3(8,15)*nine
      rrr(15)=-r7(33)-r5(20,3)*three-r5(18,8)*six-r3(9,15)*p18-r3(7,27)*three-r1(2,23)*nine
      rxyz(1)=-r3(7,31)-r1(2,27)*three
      rxyz(2)=-r4(11,21)-r2(2,24)*three
      rxyz(3)=-r4(11,20)-r2(2,23)*three
      rxyz(4)=-r5(16,9)-r3(7,16)*three-r3(7,28)-r1(2,24)*three
      rxyz(5)=-r5(16,11)-r3(7,18)*three-r3(7,30)-r1(2,26)*three
      rxyz(6)=-r5(16,10)-r3(7,17)*three-r3(7,29)-r1(2,25)*three
      rxyz(7)=-r6(22,3)-r4(11,8)*three-r4(11,18)*three-r2(2,21)*nine
      rxyz(8)=-r6(22,4)-r4(11,9)*three-r4(11,19)*three-r2(2,22)*nine
      rxyz(9)=-r5(17,9)-r3(8,16)*three
      rxyz(10)=-r5(17,10)-r3(8,17)*three
      rxyz(11)=-r5(17,11)-r3(8,18)*three
      rxyz(12)=-r6(23,3)-r4(12,8)*three-r4(12,18)-r2(6,21)*three
      rxyz(13)=-r6(23,4)-r4(12,9)*three-r4(12,19)-r2(6,22)*three
      rxyz(14)=-r6(16,4)-r4(7,9)*three-r4(7,19)-r2(4,22)*three
      rxyz(15)=-r6(16,3)-r4(7,8)*three-r4(7,18)-r2(4,21)*three
      rxyz(16)=-r6(11,4)-r4(4,9)*three-r4(11,19)-r2(2,22)*three
      rxyz(17)=-r6(11,3)-r4(4,8)*three-r4(11,18)-r2(2,21)*three
      rxyz(18)=-r6(24,3)-r4(13,8)*three-r4(11,18)-r2(2,21)*three
      rxyz(19)=-r6(24,4)-r4(13,9)*three-r4(11,19)-r2(2,22)*three
      rxyz(20)=-r5(12,9)-r3(5,16)*three
      rxyz(21)=-r5(12,10)*two-r3(5,17)*six
      rxyz(22)=-r5(12,11)-r3(5,18)*three
      rxyz(23)=-r5(7,9)-r3(2,16)*three-r3(7,28)-r1(2,24)*three
      rxyz(24)=-r5(7,10)-r3(2,17)*three-r3(7,29)-r1(2,25)*three
      rxyz(25)=-r5(7,11)-r3(2,18)*three-r3(7,30)-r1(2,26)*three
      rxyz(26)=-r4(7,20)-r2(4,23)*three
      rxyz(27)=-r4(7,21)-r2(4,24)*three
      rxyz(28)=-r6(18,3)-r4(9,8)*three-r4(7,18)-r2(4,21)*three
      rxyz(29)=-r6(18,4)-r4(9,9)*three-r4(7,19)-r2(4,22)*three
      rxyz(30)=-r6(12,3)-r4(5,8)*three-r4(12,18)-r2(6,21)*three
      rxyz(31)=-r6(12,4)-r4(5,9)*three-r4(12,19)-r2(6,22)*three
      rxyz(32)=-r5(18,9)-r3(9,16)*three-r3(7,28)-r1(2,24)*three
      rxyz(33)=-r5(18,10)-r3(9,17)*three-r3(7,29)-r1(2,25)*three
      rxyz(34)=-r5(18,11)-r3(9,18)*three-r3(7,30)-r1(2,26)*three
      rxyz(35)=-r4(12,20)-r2(6,23)*three
      rxyz(36)=-r4(12,21)-r2(6,24)*three
      rxyz(37)=-r6(7,3)-r4(2,8)*three-r4(7,18)*three-r2(4,21)*nine
      rxyz(38)=-r6(7,4)-r4(2,9)*three-r4(7,19)*three-r2(4,22)*nine
      rxyz(39)=-r6(17,3)-r4(8,8)*three
      rxyz(40)=-r6(17,4)-r4(8,9)*three
      rxyz(41)=-r5(11,9)-r3(4,16)*three
      rxyz(42)=-r5(11,10)-r3(4,17)*three
      rxyz(43)=-r5(11,11)-r3(4,18)*three
      rxyz(44)=-r6(25,3)-r4(14,8)*three-r4(12,18)*three-r2(6,21)*nine
      rxyz(45)=-r6(25,4)-r4(14,9)*three-r4(12,19)*three-r2(6,22)*nine
      call int2dddp2(eri,rrr,rxyz,qx,qz,2,2)
!
      rrr( 1)=-r7(9)+r6(5,1)*two-r5(2,1)-r5(2,3)-r5(9,8)*six+r4(5,10)*p12-r3(2,7)*six &
&          -r3(2,15)*six-r3(9,27)*three+r2(6,25)*six-r1(2,13)*three-r1(2,23)*three
      rrr( 2)=-r7(13)+r6(8,1)*two-r5(4,1)-r5(4,3)-r5(13,8)*three+r4(8,10)*six-r3(4,7)*three &
&          -r3(4,15)*three
      rrr( 3)=-r7(14)+r6(9,1)*two-r5(5,1)-r5(5,3)-r5(14,8)*three+r4(9,10)*six-r3(5,7)*three &
&          -r3(5,15)*three
      rrr( 4)=-r7(18)+r6(12,1)*two-r5(7,1)-r5(7,3)-r5(9,8)-r5(18,8)+r4(5,10)*two &
&          +r4(12,10)*two-r3(2,7)-r3(7,7)-r3(2,15)-r3(7,15)-r3(9,27)+r2(6,25)*two-r1(2,13) &
&          -r1(2,23)
      rrr( 5)=-r7(19)+r6(13,1)*two-r5(8,1)-r5(8,3)-r5(19,8)+r4(13,10)*two-r3(8,7)-r3(8,15)
      rrr( 6)=-r7(20)+r6(14,1)*two-r5(9,1)-r5(9,3)-r5(9,8)-r5(20,8)+r4(5,10)*two &
&          +r4(14,10)*two-r3(2,7)-r3(9,7)-r3(2,15)-r3(9,15)-r3(9,27)+r2(6,25)*two-r1(2,13) &
&          -r1(2,23)
      rrr( 7)=-r7(24)+r6(17,1)*two-r5(11,1)-r5(11,3)-r5(13,8)*three+r4(8,10)*six &
&          -r3(4,7)*three-r3(4,15)*three
      rrr( 8)=-r7(25)+r6(18,1)*two-r5(12,1)-r5(12,3)-r5(14,8)+r4(9,10)*two-r3(5,7)-r3(5,15)
      rrr( 9)=-r7(26)+r6(19,1)*two-r5(13,1)-r5(13,3)-r5(13,8)+r4(8,10)*two-r3(4,7)-r3(4,15)
      rrr(10)=-r7(27)+r6(20,1)*two-r5(14,1)-r5(14,3)-r5(14,8)*three+r4(9,10)*six &
&          -r3(5,7)*three-r3(5,15)*three
      rrr(11)=-r7(31)+r6(23,1)*two-r5(16,1)-r5(16,3)-r5(18,8)*six+r4(12,10)*p12-r3(7,7)*six &
&          -r3(7,15)*six-r3(9,27)*three+r2(6,25)*six-r1(2,13)*three-r1(2,23)*three
      rrr(12)=-r7(32)+r6(24,1)*two-r5(17,1)-r5(17,3)-r5(19,8)*three+r4(13,10)*six &
&          -r3(8,7)*three-r3(8,15)*three
      rrr(13)=-r7(33)+r6(25,1)*two-r5(18,1)-r5(18,3)-r5(18,8)-r5(20,8)+r4(12,10)*two &
&          +r4(14,10)*two-r3(7,7)-r3(9,7)-r3(7,15)-r3(9,15)-r3(9,27)+r2(6,25)*two-r1(2,13) &
&          -r1(2,23)
      rrr(14)=-r7(34)+r6(26,1)*two-r5(19,1)-r5(19,3)-r5(19,8)*three+r4(13,10)*six &
&          -r3(8,7)*three-r3(8,15)*three
      rrr(15)=-r7(35)+r6(27,1)*two-r5(20,1)-r5(20,3)-r5(20,8)*six+r4(14,10)*p12-r3(9,7)*six &
&          -r3(9,15)*six-r3(9,27)*three+r2(6,25)*six-r1(2,13)*three-r1(2,23)*three
      rxyz(1)=-r3(9,31)+r2(6,29)*two-r1(2,17)-r1(2,27)
      rxyz(2)=-r4(13,21)+r3(8,22)*two-r2(2,16)-r2(2,24)
      rxyz(3)=-r4(13,20)+r3(8,21)*two-r2(2,15)-r2(2,23)
      rxyz(4)=-r5(18,9)+r4(12,11)*two-r3(7,8)-r3(7,16)-r3(9,28)+r2(6,26)*two-r1(2,14) &
&             -r1(2,24)
      rxyz(5)=-r5(18,11)+r4(12,13)*two-r3(7,10)-r3(7,18)-r3(9,30)+r2(6,28)*two-r1(2,16) &
&             -r1(2,26)
      rxyz(6)=-r5(18,10)+r4(12,12)*two-r3(7,9)-r3(7,17)-r3(9,29)+r2(6,27)*two-r1(2,15) &
&             -r1(2,25)
      rxyz(7)=-r6(24,3)+r5(17,4)*two-r4(11,4)-r4(11,8)-r4(13,18)*three+r3(8,19)*six &
&             -r2(2,13)*three-r2(2,21)*three
      rxyz(8)=-r6(24,4)+r5(17,5)*two-r4(11,5)-r4(11,9)-r4(13,19)*three+r3(8,20)*six &
&             -r2(2,14)*three-r2(2,22)*three
      rxyz(9)=-r5(19,9)+r4(13,11)*two-r3(8,8)-r3(8,16)
      rxyz(10)=-r5(19,10)+r4(13,12)*two-r3(8,9)-r3(8,17)
      rxyz(11)=-r5(19,11)+r4(13,13)*two-r3(8,10)-r3(8,18)
      rxyz(12)=-r6(25,3)+r5(18,4)*two-r4(12,4)-r4(12,8)-r4(14,18)+r3(9,19)*two-r2(6,13) &
&             -r2(6,21)
      rxyz(13)=-r6(25,4)+r5(18,5)*two-r4(12,5)-r4(12,9)-r4(14,19)+r3(9,20)*two-r2(6,14) &
&             -r2(6,22)
      rxyz(14)=-r6(18,4)+r5(12,5)*two-r4(7,5)-r4(7,9)-r4(9,19)+r3(5,20)*two-r2(4,14) &
&             -r2(4,22)
      rxyz(15)=-r6(18,3)+r5(12,4)*two-r4(7,4)-r4(7,8)-r4(9,18)+r3(5,19)*two-r2(4,13) &
&             -r2(4,21)
      rxyz(16)=-r6(13,4)+r5(8,5)*two-r4(4,5)-r4(4,9)-r4(13,19)+r3(8,20)*two-r2(2,14) &
&             -r2(2,22)
      rxyz(17)=-r6(13,3)+r5(8,4)*two-r4(4,4)-r4(4,8)-r4(13,18)+r3(8,19)*two-r2(2,13) &
&             -r2(2,21)
      rxyz(18)=-r6(26,3)+r5(19,4)*two-r4(13,4)-r4(13,8)-r4(13,18)+r3(8,19)*two-r2(2,13) &
&             -r2(2,21)
      rxyz(19)=-r6(26,4)+r5(19,5)*two-r4(13,5)-r4(13,9)-r4(13,19)+r3(8,20)*two-r2(2,14) &
&             -r2(2,22)
      rxyz(20)=-r5(14,9)+r4(9,11)*two-r3(5,8)-r3(5,16)
      rxyz(21)=-r5(14,10)*two+r4(9,12)*four-r3(5,9)*two-r3(5,17)*two
      rxyz(22)=-r5(14,11)+r4(9,13)*two-r3(5,10)-r3(5,18)
      rxyz(23)=-r5(9,9)+r4(5,11)*two-r3(2,8)-r3(2,16)-r3(9,28)+r2(6,26)*two-r1(2,14)-r1(2,24)
      rxyz(24)=-r5(9,10)+r4(5,12)*two-r3(2,9)-r3(2,17)-r3(9,29)+r2(6,27)*two-r1(2,15)-r1(2,25)
      rxyz(25)=-r5(9,11)+r4(5,13)*two-r3(2,10)-r3(2,18)-r3(9,30)+r2(6,28)*two-r1(2,16)-r1(2,26)
      rxyz(26)=-r4(9,20)+r3(5,21)*two-r2(4,15)-r2(4,23)
      rxyz(27)=-r4(9,21)+r3(5,22)*two-r2(4,16)-r2(4,24)
      rxyz(28)=-r6(20,3)+r5(14,4)*two-r4(9,4)-r4(9,8)-r4(9,18)+r3(5,19)*two-r2(4,13)-r2(4,21)
      rxyz(29)=-r6(20,4)+r5(14,5)*two-r4(9,5)-r4(9,9)-r4(9,19)+r3(5,20)*two-r2(4,14)-r2(4,22)
      rxyz(30)=-r6(14,3)+r5(9,4)*two-r4(5,4)-r4(5,8)-r4(14,18)+r3(9,19)*two-r2(6,13)-r2(6,21)
      rxyz(31)=-r6(14,4)+r5(9,5)*two-r4(5,5)-r4(5,9)-r4(14,19)+r3(9,20)*two-r2(6,14)-r2(6,22)
      rxyz(32)=-r5(20,9)+r4(14,11)*two-r3(9,8)-r3(9,16)-r3(9,28)+r2(6,26)*two-r1(2,14)-r1(2,24)
      rxyz(33)=-r5(20,10)+r4(14,12)*two-r3(9,9)-r3(9,17)-r3(9,29)+r2(6,27)*two-r1(2,15)-r1(2,25)
      rxyz(34)=-r5(20,11)+r4(14,13)*two-r3(9,10)-r3(9,18)-r3(9,30)+r2(6,28)*two-r1(2,16)-r1(2,26)
      rxyz(35)=-r4(14,20)+r3(9,21)*two-r2(6,15)-r2(6,23)
      rxyz(36)=-r4(14,21)+r3(9,22)*two-r2(6,16)-r2(6,24)
      rxyz(37)=-r6(9,3)+r5(5,4)*two-r4(2,4)-r4(2,8)-r4(9,18)*three+r3(5,19)*six-r2(4,13)*three &
&              -r2(4,21)*three
      rxyz(38)=-r6(9,4)+r5(5,5)*two-r4(2,5)-r4(2,9)-r4(9,19)*three+r3(5,20)*six-r2(4,14)*three &
&              -r2(4,22)*three
      rxyz(39)=-r6(19,3)+r5(13,4)*two-r4(8,4)-r4(8,8)
      rxyz(40)=-r6(19,4)+r5(13,5)*two-r4(8,5)-r4(8,9)
      rxyz(41)=-r5(13,9)+r4(8,11)*two-r3(4,8)-r3(4,16)
      rxyz(42)=-r5(13,10)+r4(8,12)*two-r3(4,9)-r3(4,17)
      rxyz(43)=-r5(13,11)+r4(8,13)*two-r3(4,10)-r3(4,18)
      rxyz(44)=-r6(27,3)+r5(20,4)*two-r4(14,4)-r4(14,8)-r4(14,18)*three+r3(9,19)*six &
&              -r2(6,13)*three-r2(6,21)*three
      rxyz(45)=-r6(27,4)+r5(20,5)*two-r4(14,5)-r4(14,9)-r4(14,19)*three+r3(9,20)*six &
&                  -r2(6,14)*three-r2(6,22)*three
      call int2dddp2(eri,rrr,rxyz,qx,qz,3,2)
!
      do i= 1,6
        do j= 1,6
          eri(j,i,4,2)= eri(j,i,2,1)
        enddo
      enddo
      do i= 1,6
        do j= 1,6
          eri(j,i,5,2)= eri(j,i,6,1)
        enddo
      enddo
!
      rrr( 1)=-r7(8)+r6(4,1)-r5(3,3)-r5(8,8)*six+r4(1,2)+r4(4,10)*six-r3(3,15)*six &
&          -r3(8,27)*three+r2(1,5)*six+r2(2,25)*three-r1(3,23)*three+r0(6)*three
      rrr( 2)=-r7(12)+r6(7,1)-r5(5,3)-r5(12,8)*three+r4(2,2)+r4(7,10)*three-r3(5,15)*three &
&          +r2(4,5)*three
      rrr( 3)=-r7(13)+r6(8,1)-r5(6,3)-r5(13,8)*three+r4(3,2)+r4(8,10)*three-r3(6,15)*three &
&          +r2(5,5)*three
      rrr( 4)=-r7(17)+r6(11,1)-r5(8,3)-r5(8,8)-r5(17,8)+r4(4,2)+r4(4,10)+r4(11,10)-r3(3,15) &
&          -r3(8,15)-r3(8,27)+r2(1,5)+r2(2,5)+r2(2,25)-r1(3,23)+r0(6)
      rrr( 5)=-r7(18)+r6(12,1)-r5(9,3)-r5(18,8)+r4(5,2)+r4(12,10)-r3(9,15)+r2(6,5)
      rrr( 6)=-r7(19)+r6(13,1)-r5(10,3)-r5(8,8)-r5(19,8)+r4(6,2)+r4(4,10)+r4(13,10)-r3(3,15) &
&          -r3(10,15)-r3(8,27)+r2(1,5)+r2(3,5)+r2(2,25)-r1(3,23)+r0(6)
      rrr( 7)=-r7(23)+r6(16,1)-r5(12,3)-r5(12,8)*three+r4(7,2)+r4(7,10)*three-r3(5,15)*three &
&          +r2(4,5)*three
      rrr( 8)=-r7(24)+r6(17,1)-r5(13,3)-r5(13,8)+r4(8,2)+r4(8,10)-r3(6,15)+r2(5,5)
      rrr( 9)=-r7(25)+r6(18,1)-r5(14,3)-r5(12,8)+r4(9,2)+r4(7,10)-r3(5,15)+r2(4,5)
      rrr(10)=-r7(26)+r6(19,1)-r5(15,3)-r5(13,8)*three+r4(10,2)+r4(8,10)*three-r3(6,15)*three &
&          +r2(5,5)*three
      rrr(11)=-r7(30)+r6(22,1)-r5(17,3)-r5(17,8)*six+r4(11,2)+r4(11,10)*six-r3(8,15)*six &
&          -r3(8,27)*three+r2(2,5)*six+r2(2,25)*three-r1(3,23)*three+r0(6)*three
      rrr(12)=-r7(31)+r6(23,1)-r5(18,3)-r5(18,8)*three+r4(12,2)+r4(12,10)*three &
&          -r3(9,15)*three+r2(6,5)*three
      rrr(13)=-r7(32)+r6(24,1)-r5(19,3)-r5(17,8)-r5(19,8)+r4(13,2)+r4(11,10)+r4(13,10) &
&          -r3(8,15)-r3(10,15)-r3(8,27)+r2(2,5)+r2(3,5)+r2(2,25)-r1(3,23)+r0(6)
      rrr(14)=-r7(33)+r6(25,1)-r5(20,3)-r5(18,8)*three+r4(14,2)+r4(12,10)*three &
&          -r3(9,15)*three+r2(6,5)*three
      rrr(15)=-r7(34)+r6(26,1)-r5(21,3)-r5(19,8)*six+r4(15,2)+r4(13,10)*six-r3(10,15)*six &
&          -r3(8,27)*three+r2(3,5)*six+r2(2,25)*three-r1(3,23)*three+r0(6)*three
      rxyz(1)=-r3(8,31)+r2(2,29)-r1(3,27)+r0(10)
      rxyz(2)=-r4(12,21)+r3(7,22)-r2(6,24)+r1(2,8)
      rxyz(3)=-r4(12,20)+r3(7,21)-r2(6,23)+r1(2,7)
      rxyz(4)=-r5(17,9)+r4(11,11)-r3(8,16)-r3(8,28)+r2(2,6)+r2(2,26)-r1(3,24)+r0(7)
      rxyz(5)=-r5(17,11)+r4(11,13)-r3(8,18)-r3(8,30)+r2(2,8)+r2(2,28)-r1(3,26)+r0(9)
      rxyz(6)=-r5(17,10)+r4(11,12)-r3(8,17)-r3(8,29)+r2(2,7)+r2(2,27)-r1(3,25)+r0(8)
      rxyz(7)=-r6(23,3)+r5(16,4)-r4(12,8)-r4(12,18)*three+r3(7,3)+r3(7,19)*three &
&             -r2(6,21)*three+r1(2,5)*three
      rxyz(8)=-r6(23,4)+r5(16,5)-r4(12,9)-r4(12,19)*three+r3(7,4)+r3(7,20)*three &
&             -r2(6,22)*three+r1(2,6)*three
      rxyz(9)=-r5(18,9)+r4(12,11)-r3(9,16)+r2(6,6)
      rxyz(10)=-r5(18,10)+r4(12,12)-r3(9,17)+r2(6,7)
      rxyz(11)=-r5(18,11)+r4(12,13)-r3(9,18)+r2(6,8)
      rxyz(12)=-r6(24,3)+r5(17,4)-r4(13,8)-r4(13,18)+r3(8,3)+r3(8,19)-r2(3,21)+r1(3,5)
      rxyz(13)=-r6(24,4)+r5(17,5)-r4(13,9)-r4(13,19)+r3(8,4)+r3(8,20)-r2(3,22)+r1(3,6)
      rxyz(14)=-r6(17,4)+r5(11,5)-r4(8,9)-r4(8,19)+r3(4,4)+r3(4,20)-r2(5,22)+r1(1,6)
      rxyz(15)=-r6(17,3)+r5(11,4)-r4(8,8)-r4(8,18)+r3(4,3)+r3(4,19)-r2(5,21)+r1(1,5)
      rxyz(16)=-r6(12,4)+r5(7,5)-r4(5,9)-r4(12,19)+r3(2,4)+r3(7,20)-r2(6,22)+r1(2,6)
      rxyz(17)=-r6(12,3)+r5(7,4)-r4(5,8)-r4(12,18)+r3(2,3)+r3(7,19)-r2(6,21)+r1(2,5)
      rxyz(18)=-r6(25,3)+r5(18,4)-r4(14,8)-r4(12,18)+r3(9,3)+r3(7,19)-r2(6,21)+r1(2,5)
      rxyz(19)=-r6(25,4)+r5(18,5)-r4(14,9)-r4(12,19)+r3(9,4)+r3(7,20)-r2(6,22)+r1(2,6)
      rxyz(20)=-r5(13,9)+r4(8,11)-r3(6,16)+r2(5,6)
      rxyz(21)=-r5(13,10)*two+r4(8,12)*two-r3(6,17)*two+r2(5,7)*two
      rxyz(22)=-r5(13,11)+r4(8,13)-r3(6,18)+r2(5,8)
      rxyz(23)=-r5(8,9)+r4(4,11)-r3(3,16)-r3(8,28)+r2(1,6)+r2(2,26)-r1(3,24)+r0(7)
      rxyz(24)=-r5(8,10)+r4(4,12)-r3(3,17)-r3(8,29)+r2(1,7)+r2(2,27)-r1(3,25)+r0(8)
      rxyz(25)=-r5(8,11)+r4(4,13)-r3(3,18)-r3(8,30)+r2(1,8)+r2(2,28)-r1(3,26)+r0(9)
      rxyz(26)=-r4(8,20)+r3(4,21)-r2(5,23)+r1(1,7)
      rxyz(27)=-r4(8,21)+r3(4,22)-r2(5,24)+r1(1,8)
      rxyz(28)=-r6(19,3)+r5(13,4)-r4(10,8)-r4(8,18)+r3(6,3)+r3(4,19)-r2(5,21)+r1(1,5)
      rxyz(29)=-r6(19,4)+r5(13,5)-r4(10,9)-r4(8,19)+r3(6,4)+r3(4,20)-r2(5,22)+r1(1,6)
      rxyz(30)=-r6(13,3)+r5(8,4)-r4(6,8)-r4(13,18)+r3(3,3)+r3(8,19)-r2(3,21)+r1(3,5)
      rxyz(31)=-r6(13,4)+r5(8,5)-r4(6,9)-r4(13,19)+r3(3,4)+r3(8,20)-r2(3,22)+r1(3,6)
      rxyz(32)=-r5(19,9)+r4(13,11)-r3(10,16)-r3(8,28)+r2(3,6)+r2(2,26)-r1(3,24)+r0(7)
      rxyz(33)=-r5(19,10)+r4(13,12)-r3(10,17)-r3(8,29)+r2(3,7)+r2(2,27)-r1(3,25)+r0(8)
      rxyz(34)=-r5(19,11)+r4(13,13)-r3(10,18)-r3(8,30)+r2(3,8)+r2(2,28)-r1(3,26)+r0(9)
      rxyz(35)=-r4(13,20)+r3(8,21)-r2(3,23)+r1(3,7)
      rxyz(36)=-r4(13,21)+r3(8,22)-r2(3,24)+r1(3,8)
      rxyz(37)=-r6(8,3)+r5(4,4)-r4(3,8)-r4(8,18)*three+r3(1,3)+r3(4,19)*three-r2(5,21)*three &
&              +r1(1,5)*three
      rxyz(38)=-r6(8,4)+r5(4,5)-r4(3,9)-r4(8,19)*three+r3(1,4)+r3(4,20)*three-r2(5,22)*three &
&              +r1(1,6)*three
      rxyz(39)=-r6(18,3)+r5(12,4)-r4(9,8)+r3(5,3)
      rxyz(40)=-r6(18,4)+r5(12,5)-r4(9,9)+r3(5,4)
      rxyz(41)=-r5(12,9)+r4(7,11)-r3(5,16)+r2(4,6)
      rxyz(42)=-r5(12,10)+r4(7,12)-r3(5,17)+r2(4,7)
      rxyz(43)=-r5(12,11)+r4(7,13)-r3(5,18)+r2(4,8)
      rxyz(44)=-r6(26,3)+r5(19,4)-r4(15,8)-r4(13,18)*three+r3(10,3)+r3(8,19)*three &
&              -r2(3,21)*three+r1(3,5)*three
      rxyz(45)=-r6(26,4)+r5(19,5)-r4(15,9)-r4(13,19)*three+r3(10,4)+r3(8,20)*three &
&                  -r2(3,22)*three+r1(3,6)*three
      call int2dddp2(eri,rrr,rxyz,qx,qz,6,2)
!
      rrr( 1)=-r7(3)+r6(1,2)-r5(3,3)-r5(3,8)*six+r4(1,3)+r4(1,14)*six-r3(3,15)*six &
&          -r3(3,27)*three+r2(1,9)*six+r2(1,30)*three-r1(3,23)*three+r0(11)*three
      rrr( 2)=-r7(5)+r6(2,2)-r5(5,3)-r5(5,8)*three+r4(2,3)+r4(2,14)*three-r3(5,15)*three &
&          +r2(4,9)*three
      rrr( 3)=-r7(6)+r6(3,2)-r5(6,3)-r5(6,8)*three+r4(3,3)+r4(3,14)*three-r3(6,15)*three &
&          +r2(5,9)*three
      rrr( 4)=-r7(8)+r6(4,2)-r5(8,3)-r5(3,8)-r5(8,8)+r4(4,3)+r4(1,14)+r4(4,14)-r3(3,15) &
&          -r3(8,15)-r3(3,27)+r2(1,9)+r2(2,9)+r2(1,30)-r1(3,23)+r0(11)
      rrr( 5)=-r7(9)+r6(5,2)-r5(9,3)-r5(9,8)+r4(5,3)+r4(5,14)-r3(9,15)+r2(6,9)
      rrr( 6)=-r7(10)+r6(6,2)-r5(10,3)-r5(3,8)-r5(10,8)+r4(6,3)+r4(1,14)+r4(6,14)-r3(3,15) &
&          -r3(10,15)-r3(3,27)+r2(1,9)+r2(3,9)+r2(1,30)-r1(3,23)+r0(11)
      rrr( 7)=-r7(12)+r6(7,2)-r5(12,3)-r5(5,8)*three+r4(7,3)+r4(2,14)*three-r3(5,15)*three &
&          +r2(4,9)*three
      rrr( 8)=-r7(13)+r6(8,2)-r5(13,3)-r5(6,8)+r4(8,3)+r4(3,14)-r3(6,15)+r2(5,9)
      rrr( 9)=-r7(14)+r6(9,2)-r5(14,3)-r5(5,8)+r4(9,3)+r4(2,14)-r3(5,15)+r2(4,9)
      rrr(10)=-r7(15)+r6(10,2)-r5(15,3)-r5(6,8)*three+r4(10,3)+r4(3,14)*three-r3(6,15)*three &
&          +r2(5,9)*three
      rrr(11)=-r7(17)+r6(11,2)-r5(17,3)-r5(8,8)*six+r4(11,3)+r4(4,14)*six-r3(8,15)*six &
&          -r3(3,27)*three+r2(2,9)*six+r2(1,30)*three-r1(3,23)*three+r0(11)*three
      rrr(12)=-r7(18)+r6(12,2)-r5(18,3)-r5(9,8)*three+r4(12,3)+r4(5,14)*three-r3(9,15)*three &
&          +r2(6,9)*three
      rrr(13)=-r7(19)+r6(13,2)-r5(19,3)-r5(8,8)-r5(10,8)+r4(13,3)+r4(4,14)+r4(6,14)-r3(8,15) &
&          -r3(10,15)-r3(3,27)+r2(2,9)+r2(3,9)+r2(1,30)-r1(3,23)+r0(11)
      rrr(14)=-r7(20)+r6(14,2)-r5(20,3)-r5(9,8)*three+r4(14,3)+r4(5,14)*three-r3(9,15)*three &
&          +r2(6,9)*three
      rrr(15)=-r7(21)+r6(15,2)-r5(21,3)-r5(10,8)*six+r4(15,3)+r4(6,14)*six-r3(10,15)*six &
&          -r3(3,27)*three+r2(3,9)*six+r2(1,30)*three-r1(3,23)*three+r0(11)*three
      rxyz(1)=-r3(3,31)+r2(1,34)-r1(3,27)+r0(15)
      rxyz(2)=-r4(5,21)+r3(2,26)-r2(6,24)+r1(2,12)
      rxyz(3)=-r4(5,20)+r3(2,25)-r2(6,23)+r1(2,11)
      rxyz(4)=-r5(8,9)+r4(4,15)-r3(8,16)-r3(3,28)+r2(2,10)+r2(1,31)-r1(3,24)+r0(12)
      rxyz(5)=-r5(8,11)+r4(4,17)-r3(8,18)-r3(3,30)+r2(2,12)+r2(1,33)-r1(3,26)+r0(14)
      rxyz(6)=-r5(8,10)+r4(4,16)-r3(8,17)-r3(3,29)+r2(2,11)+r2(1,32)-r1(3,25)+r0(13)
      rxyz(7)=-r6(12,3)+r5(7,6)-r4(12,8)-r4(5,18)*three+r3(7,5)+r3(2,23)*three &
&             -r2(6,21)*three+r1(2,9)*three
      rxyz(8)=-r6(12,4)+r5(7,7)-r4(12,9)-r4(5,19)*three+r3(7,6)+r3(2,24)*three &
&             -r2(6,22)*three+r1(2,10)*three
      rxyz(9)=-r5(9,9)+r4(5,15)-r3(9,16)+r2(6,10)
      rxyz(10)=-r5(9,10)+r4(5,16)-r3(9,17)+r2(6,11)
      rxyz(11)=-r5(9,11)+r4(5,17)-r3(9,18)+r2(6,12)
      rxyz(12)=-r6(13,3)+r5(8,6)-r4(13,8)-r4(6,18)+r3(8,5)+r3(3,23)-r2(3,21)+r1(3,9)
      rxyz(13)=-r6(13,4)+r5(8,7)-r4(13,9)-r4(6,19)+r3(8,6)+r3(3,24)-r2(3,22)+r1(3,10)
      rxyz(14)=-r6(8,4)+r5(4,7)-r4(8,9)-r4(3,19)+r3(4,6)+r3(1,24)-r2(5,22)+r1(1,10)
      rxyz(15)=-r6(8,3)+r5(4,6)-r4(8,8)-r4(3,18)+r3(4,5)+r3(1,23)-r2(5,21)+r1(1,9)
      rxyz(16)=-r6(5,4)+r5(2,7)-r4(5,9)-r4(5,19)+r3(2,6)+r3(2,24)-r2(6,22)+r1(2,10)
      rxyz(17)=-r6(5,3)+r5(2,6)-r4(5,8)-r4(5,18)+r3(2,5)+r3(2,23)-r2(6,21)+r1(2,9)
      rxyz(18)=-r6(14,3)+r5(9,6)-r4(14,8)-r4(5,18)+r3(9,5)+r3(2,23)-r2(6,21)+r1(2,9)
      rxyz(19)=-r6(14,4)+r5(9,7)-r4(14,9)-r4(5,19)+r3(9,6)+r3(2,24)-r2(6,22)+r1(2,10)
      rxyz(20)=-r5(6,9)+r4(3,15)-r3(6,16)+r2(5,10)
      rxyz(21)=-r5(6,10)*two+r4(3,16)*two-r3(6,17)*two+r2(5,11)*two
      rxyz(22)=-r5(6,11)+r4(3,17)-r3(6,18)+r2(5,12)
      rxyz(23)=-r5(3,9)+r4(1,15)-r3(3,16)-r3(3,28)+r2(1,10)+r2(1,31)-r1(3,24)+r0(12)
      rxyz(24)=-r5(3,10)+r4(1,16)-r3(3,17)-r3(3,29)+r2(1,11)+r2(1,32)-r1(3,25)+r0(13)
      rxyz(25)=-r5(3,11)+r4(1,17)-r3(3,18)-r3(3,30)+r2(1,12)+r2(1,33)-r1(3,26)+r0(14)
      rxyz(26)=-r4(3,20)+r3(1,25)-r2(5,23)+r1(1,11)
      rxyz(27)=-r4(3,21)+r3(1,26)-r2(5,24)+r1(1,12)
      rxyz(28)=-r6(10,3)+r5(6,6)-r4(10,8)-r4(3,18)+r3(6,5)+r3(1,23)-r2(5,21)+r1(1,9)
      rxyz(29)=-r6(10,4)+r5(6,7)-r4(10,9)-r4(3,19)+r3(6,6)+r3(1,24)-r2(5,22)+r1(1,10)
      rxyz(30)=-r6(6,3)+r5(3,6)-r4(6,8)-r4(6,18)+r3(3,5)+r3(3,23)-r2(3,21)+r1(3,9)
      rxyz(31)=-r6(6,4)+r5(3,7)-r4(6,9)-r4(6,19)+r3(3,6)+r3(3,24)-r2(3,22)+r1(3,10)
      rxyz(32)=-r5(10,9)+r4(6,15)-r3(10,16)-r3(3,28)+r2(3,10)+r2(1,31)-r1(3,24)+r0(12)
      rxyz(33)=-r5(10,10)+r4(6,16)-r3(10,17)-r3(3,29)+r2(3,11)+r2(1,32)-r1(3,25)+r0(13)
      rxyz(34)=-r5(10,11)+r4(6,17)-r3(10,18)-r3(3,30)+r2(3,12)+r2(1,33) -r1(3,26)+r0(14)
      rxyz(35)=-r4(6,20)+r3(3,25)-r2(3,23)+r1(3,11)
      rxyz(36)=-r4(6,21)+r3(3,26)-r2(3,24)+r1(3,12)
      rxyz(37)=-r6(3,3)+r5(1,6)-r4(3,8)-r4(3,18)*three+r3(1,5)+r3(1,23)*three-r2(5,21)*three &
&              +r1(1,9)*three
      rxyz(38)=-r6(3,4)+r5(1,7)-r4(3,9)-r4(3,19)*three+r3(1,6)+r3(1,24)*three-r2(5,22)*three &
&              +r1(1,10)*three
      rxyz(39)=-r6(9,3)+r5(5,6)-r4(9,8)+r3(5,5)
      rxyz(40)=-r6(9,4)+r5(5,7)-r4(9,9)+r3(5,6)
      rxyz(41)=-r5(5,9)+r4(2,15)-r3(5,16)+r2(4,10)
      rxyz(42)=-r5(5,10)+r4(2,16)-r3(5,17)+r2(4,11)
      rxyz(43)=-r5(5,11)+r4(2,17)-r3(5,18)+r2(4,12)
      rxyz(44)=-r6(15,3)+r5(10,6)-r4(15,8)-r4(6,18)*three+r3(10,5)+r3(3,23)*three-r2(3,21)*three &
&              +r1(3,9)*three
      rxyz(45)=-r6(15,4)+r5(10,7)-r4(15,9)-r4(6,19)*three+r3(10,6)+r3(3,24)*three-r2(3,22)*three &
&              +r1(3,10)*three
      call int2dddp2(eri,rrr,rxyz,qx,qz,1,3)
!
      rrr( 1)=-r7(8)+r6(4,2)-r5(3,3)-r5(8,8)*six+r4(1,3)+r4(4,14)*six-r3(3,15)*six &
&          -r3(8,27)*three+r2(1,9)*six+r2(2,30)*three-r1(3,23)*three+r0(11)*three
      rrr( 2)=-r7(12)+r6(7,2)-r5(5,3)-r5(12,8)*three+r4(2,3)+r4(7,14)*three-r3(5,15)*three &
&          +r2(4,9)*three
      rrr( 3)=-r7(13)+r6(8,2)-r5(6,3)-r5(13,8)*three+r4(3,3)+r4(8,14)*three-r3(6,15)*three &
&          +r2(5,9)*three
      rrr( 4)=-r7(17)+r6(11,2)-r5(8,3)-r5(8,8)-r5(17,8)+r4(4,3)+r4(4,14)+r4(11,14)-r3(3,15) &
&          -r3(8,15)-r3(8,27)+r2(1,9)+r2(2,9)+r2(2,30)-r1(3,23)+r0(11)
      rrr( 5)=-r7(18)+r6(12,2)-r5(9,3)-r5(18,8)+r4(5,3)+r4(12,14)-r3(9,15)+r2(6,9)
      rrr( 6)=-r7(19)+r6(13,2)-r5(10,3)-r5(8,8)-r5(19,8)+r4(6,3)+r4(4,14)+r4(13,14)-r3(3,15) &
&          -r3(10,15)-r3(8,27)+r2(1,9)+r2(3,9)+r2(2,30)-r1(3,23)+r0(11)
      rrr( 7)=-r7(23)+r6(16,2)-r5(12,3)-r5(12,8)*three+r4(7,3)+r4(7,14)*three-r3(5,15)*three &
&          +r2(4,9)*three
      rrr( 8)=-r7(24)+r6(17,2)-r5(13,3)-r5(13,8)+r4(8,3)+r4(8,14)-r3(6,15)+r2(5,9)
      rrr( 9)=-r7(25)+r6(18,2)-r5(14,3)-r5(12,8)+r4(9,3)+r4(7,14)-r3(5,15)+r2(4,9)
      rrr(10)=-r7(26)+r6(19,2)-r5(15,3)-r5(13,8)*three+r4(10,3)+r4(8,14)*three-r3(6,15)*three &
&          +r2(5,9)*three
      rrr(11)=-r7(30)+r6(22,2)-r5(17,3)-r5(17,8)*six+r4(11,3)+r4(11,14)*six-r3(8,15)*six &
&          -r3(8,27)*three+r2(2,9)*six+r2(2,30)*three-r1(3,23)*three+r0(11)*three
      rrr(12)=-r7(31)+r6(23,2)-r5(18,3)-r5(18,8)*three+r4(12,3)+r4(12,14)*three &
&          -r3(9,15)*three+r2(6,9)*three
      rrr(13)=-r7(32)+r6(24,2)-r5(19,3)-r5(17,8)-r5(19,8)+r4(13,3)+r4(11,14)+r4(13,14) &
&          -r3(8,15)-r3(10,15)-r3(8,27)+r2(2,9)+r2(3,9)+r2(2,30)-r1(3,23)+r0(11)
      rrr(14)=-r7(33)+r6(25,2)-r5(20,3)-r5(18,8)*three+r4(14,3)+r4(12,14)*three &
&          -r3(9,15)*three+r2(6,9)*three
      rrr(15)=-r7(34)+r6(26,2)-r5(21,3)-r5(19,8)*six+r4(15,3)+r4(13,14)*six-r3(10,15)*six &
&          -r3(8,27)*three+r2(3,9)*six+r2(2,30)*three-r1(3,23)*three+r0(11)*three
      rxyz(1)=-r3(8,31)+r2(2,34)-r1(3,27)+r0(15)
      rxyz(2)=-r4(12,21)+r3(7,26)-r2(6,24)+r1(2,12)
      rxyz(3)=-r4(12,20)+r3(7,25)-r2(6,23)+r1(2,11)
      rxyz(4)=-r5(17,9)+r4(11,15)-r3(8,16)-r3(8,28)+r2(2,10)+r2(2,31)-r1(3,24)+r0(12)
      rxyz(5)=-r5(17,11)+r4(11,17)-r3(8,18)-r3(8,30)+r2(2,12)+r2(2,33)-r1(3,26)+r0(14)
      rxyz(6)=-r5(17,10)+r4(11,16)-r3(8,17)-r3(8,29)+r2(2,11)+r2(2,32)-r1(3,25)+r0(13)
      rxyz(7)=-r6(23,3)+r5(16,6)-r4(12,8)-r4(12,18)*three+r3(7,5)+r3(7,23)*three &
&             -r2(6,21)*three+r1(2,9)*three
      rxyz(8)=-r6(23,4)+r5(16,7)-r4(12,9)-r4(12,19)*three+r3(7,6)+r3(7,24)*three &
&             -r2(6,22)*three+r1(2,10)*three
      rxyz(9)=-r5(18,9)+r4(12,15)-r3(9,16)+r2(6,10)
      rxyz(10)=-r5(18,10)+r4(12,16)-r3(9,17)+r2(6,11)
      rxyz(11)=-r5(18,11)+r4(12,17)-r3(9,18)+r2(6,12)
      rxyz(12)=-r6(24,3)+r5(17,6)-r4(13,8)-r4(13,18)+r3(8,5)+r3(8,23)-r2(3,21)+r1(3,9)
      rxyz(13)=-r6(24,4)+r5(17,7)-r4(13,9)-r4(13,19)+r3(8,6)+r3(8,24)-r2(3,22)+r1(3,10)
      rxyz(14)=-r6(17,4)+r5(11,7)-r4(8,9)-r4(8,19)+r3(4,6)+r3(4,24)-r2(5,22)+r1(1,10)
      rxyz(15)=-r6(17,3)+r5(11,6)-r4(8,8)-r4(8,18)+r3(4,5)+r3(4,23)-r2(5,21)+r1(1,9)
      rxyz(16)=-r6(12,4)+r5(7,7)-r4(5,9)-r4(12,19)+r3(2,6)+r3(7,24)-r2(6,22)+r1(2,10)
      rxyz(17)=-r6(12,3)+r5(7,6)-r4(5,8)-r4(12,18)+r3(2,5)+r3(7,23)-r2(6,21)+r1(2,9)
      rxyz(18)=-r6(25,3)+r5(18,6)-r4(14,8)-r4(12,18)+r3(9,5)+r3(7,23)-r2(6,21)+r1(2,9)
      rxyz(19)=-r6(25,4)+r5(18,7)-r4(14,9)-r4(12,19)+r3(9,6)+r3(7,24)-r2(6,22)+r1(2,10)
      rxyz(20)=-r5(13,9)+r4(8,15)-r3(6,16)+r2(5,10)
      rxyz(21)=-r5(13,10)*two+r4(8,16)*two-r3(6,17)*two+r2(5,11)*two
      rxyz(22)=-r5(13,11)+r4(8,17)-r3(6,18)+r2(5,12)
      rxyz(23)=-r5(8,9)+r4(4,15)-r3(3,16)-r3(8,28)+r2(1,10)+r2(2,31)-r1(3,24)+r0(12)
      rxyz(24)=-r5(8,10)+r4(4,16)-r3(3,17)-r3(8,29)+r2(1,11)+r2(2,32)-r1(3,25)+r0(13)
      rxyz(25)=-r5(8,11)+r4(4,17)-r3(3,18)-r3(8,30)+r2(1,12)+r2(2,33)-r1(3,26)+r0(14)
      rxyz(26)=-r4(8,20)+r3(4,25)-r2(5,23)+r1(1,11)
      rxyz(27)=-r4(8,21)+r3(4,26)-r2(5,24)+r1(1,12)
      rxyz(28)=-r6(19,3)+r5(13,6)-r4(10,8)-r4(8,18)+r3(6,5)+r3(4,23)-r2(5,21)+r1(1,9)
      rxyz(29)=-r6(19,4)+r5(13,7)-r4(10,9)-r4(8,19)+r3(6,6)+r3(4,24)-r2(5,22)+r1(1,10)
      rxyz(30)=-r6(13,3)+r5(8,6)-r4(6,8)-r4(13,18)+r3(3,5)+r3(8,23)-r2(3,21)+r1(3,9)
      rxyz(31)=-r6(13,4)+r5(8,7)-r4(6,9)-r4(13,19)+r3(3,6)+r3(8,24)-r2(3,22)+r1(3,10)
      rxyz(32)=-r5(19,9)+r4(13,15)-r3(10,16)-r3(8,28)+r2(3,10)+r2(2,31)-r1(3,24)+r0(12)
      rxyz(33)=-r5(19,10)+r4(13,16)-r3(10,17)-r3(8,29)+r2(3,11)+r2(2,32)-r1(3,25)+r0(13)
      rxyz(34)=-r5(19,11)+r4(13,17)-r3(10,18)-r3(8,30)+r2(3,12)+r2(2,33)-r1(3,26)+r0(14)
      rxyz(35)=-r4(13,20)+r3(8,25)-r2(3,23)+r1(3,11)
      rxyz(36)=-r4(13,21)+r3(8,26)-r2(3,24)+r1(3,12)
      rxyz(37)=-r6(8,3)+r5(4,6)-r4(3,8)-r4(8,18)*three+r3(1,5)+r3(4,23)*three-r2(5,21)*three &
&              +r1(1,9)*three
      rxyz(38)=-r6(8,4)+r5(4,7)-r4(3,9)-r4(8,19)*three+r3(1,6)+r3(4,24)*three-r2(5,22)*three &
&              +r1(1,10)*three
      rxyz(39)=-r6(18,3)+r5(12,6)-r4(9,8)+r3(5,5)
      rxyz(40)=-r6(18,4)+r5(12,7)-r4(9,9)+r3(5,6)
      rxyz(41)=-r5(12,9)+r4(7,15)-r3(5,16)+r2(4,10)
      rxyz(42)=-r5(12,10)+r4(7,16)-r3(5,17)+r2(4,11)
      rxyz(43)=-r5(12,11)+r4(7,17)-r3(5,18)+r2(4,12)
      rxyz(44)=-r6(26,3)+r5(19,6)-r4(15,8)-r4(13,18)*three+r3(10,5)+r3(8,23)*three &
&              -r2(3,21)*three+r1(3,9)*three
      rxyz(45)=-r6(26,4)+r5(19,7)-r4(15,9)-r4(13,19)*three+r3(10,6)+r3(8,24)*three &
&              -r2(3,22)*three+r1(3,10)*three
      call int2dddp2(eri,rrr,rxyz,qx,qz,2,3)
!
      rrr( 1)=-r7(10)+r6(6,1)*two+r6(6,2)-r5(3,1)-r5(3,2)*two-r5(3,3)*three-r5(10,8)*six &
&          +r4(1,1)+r4(1,2)*two+r4(1,3)+r4(6,10)*p12+r4(6,14)*six-r3(3,7)*six-r3(3,11)*p12 &
&          -r3(3,15)*p18-r3(10,27)*three+r2(1,1)*six+r2(1,5)*p12+r2(1,9)*six+r2(3,25)*six &
&          +r2(3,30)*three-r1(3,13)*three-r1(3,18)*six-r1(3,23)*nine+r0(1)*three+r0(6)*six &
&          +r0(11)*three
      rrr( 2)=-r7(14)+r6(9,1)*two+r6(9,2)-r5(5,1)-r5(5,2)*two-r5(5,3)*three-r5(14,8)*three &
&          +r4(2,1)+r4(2,2)*two+r4(2,3)+r4(9,10)*six+r4(9,14)*three-r3(5,7)*three &
&          -r3(5,11)*six-r3(5,15)*nine+r2(4,1)*three+r2(4,5)*six+r2(4,9)*three
      rrr( 3)=-r7(15)+r6(10,1)*two+r6(10,2)-r5(6,1)-r5(6,2)*two-r5(6,3)*three-r5(15,8)*three &
&          +r4(3,1)+r4(3,2)*two+r4(3,3)+r4(10,10)*six+r4(10,14)*three-r3(6,7)*three &
&          -r3(6,11)*six-r3(6,15)*nine+r2(5,1)*three+r2(5,5)*six+r2(5,9)*three
      rrr( 4)=-r7(19)+r6(13,1)*two+r6(13,2)-r5(8,1)-r5(8,2)*two-r5(8,3)*three-r5(10,8) &
&          -r5(19,8)+r4(4,1)+r4(4,2)*two+r4(4,3)+r4(6,10)*two+r4(13,10)*two+r4(6,14) &
&          +r4(13,14)-r3(3,7)-r3(8,7)-r3(3,11)*two-r3(8,11)*two-r3(3,15)*three &
&          -r3(8,15)*three-r3(10,27)+r2(1,1)+r2(2,1)+r2(1,5)*two+r2(2,5)*two+r2(1,9) &
&          +r2(2,9)+r2(3,25)*two+r2(3,30)-r1(3,13)-r1(3,18)*two-r1(3,23)*three+r0(1) &
&          +r0(6)*two+r0(11)
      rrr( 5)=-r7(20)+r6(14,1)*two+r6(14,2)-r5(9,1)-r5(9,2)*two-r5(9,3)*three-r5(20,8) &
&          +r4(5,1)+r4(5,2)*two+r4(5,3)+r4(14,10)*two+r4(14,14)-r3(9,7)-r3(9,11)*two &
&          -r3(9,15)*three+r2(6,1)+r2(6,5)*two+r2(6,9)
      rrr( 6)=-r7(21)+r6(15,1)*two+r6(15,2)-r5(10,1)-r5(10,2)*two-r5(10,3)*three-r5(10,8) &
&          -r5(21,8)+r4(6,1)+r4(6,2)*two+r4(6,3)+r4(6,10)*two+r4(15,10)*two+r4(6,14) &
&          +r4(15,14)-r3(3,7)-r3(10,7)-r3(3,11)*two-r3(10,11)*two-r3(3,15)*three &
&          -r3(10,15)*three-r3(10,27)+r2(1,1)+r2(3,1)+r2(1,5)*two+r2(3,5)*two+r2(1,9) &
&          +r2(3,9)+r2(3,25)*two+r2(3,30)-r1(3,13)-r1(3,18)*two-r1(3,23)*three+r0(1) &
&          +r0(6)*two+r0(11)
      rrr( 7)=-r7(25)+r6(18,1)*two+r6(18,2)-r5(12,1)-r5(12,2)*two-r5(12,3)*three &
&          -r5(14,8)*three+r4(7,1)+r4(7,2)*two+r4(7,3)+r4(9,10)*six+r4(9,14)*three &
&          -r3(5,7)*three-r3(5,11)*six-r3(5,15)*nine+r2(4,1)*three+r2(4,5)*six &
&          +r2(4,9)*three
      rrr( 8)=-r7(26)+r6(19,1)*two+r6(19,2)-r5(13,1)-r5(13,2)*two-r5(13,3)*three-r5(15,8) &
&          +r4(8,1)+r4(8,2)*two+r4(8,3)+r4(10,10)*two+r4(10,14)-r3(6,7)-r3(6,11)*two &
&          -r3(6,15)*three+r2(5,1)+r2(5,5)*two+r2(5,9)
      rrr( 9)=-r7(27)+r6(20,1)*two+r6(20,2)-r5(14,1)-r5(14,2)*two-r5(14,3)*three-r5(14,8) &
&          +r4(9,1)+r4(9,2)*two+r4(9,3)+r4(9,10)*two+r4(9,14)-r3(5,7)-r3(5,11)*two &
&          -r3(5,15)*three+r2(4,1)+r2(4,5)*two+r2(4,9)
      rrr(10)=-r7(28)+r6(21,1)*two+r6(21,2)-r5(15,1)-r5(15,2)*two-r5(15,3)*three &
&          -r5(15,8)*three+r4(10,1)+r4(10,2)*two+r4(10,3)+r4(10,10)*six+r4(10,14)*three &
&          -r3(6,7)*three-r3(6,11)*six-r3(6,15)*nine+r2(5,1)*three+r2(5,5)*six &
&          +r2(5,9)*three
      rrr(11)=-r7(32)+r6(24,1)*two+r6(24,2)-r5(17,1)-r5(17,2)*two-r5(17,3)*three-r5(19,8)*six &
&          +r4(11,1)+r4(11,2)*two+r4(11,3)+r4(13,10)*p12+r4(13,14)*six-r3(8,7)*six &
&          -r3(8,11)*p12-r3(8,15)*p18-r3(10,27)*three+r2(2,1)*six+r2(2,5)*p12+r2(2,9)*six &
&          +r2(3,25)*six+r2(3,30)*three-r1(3,13)*three-r1(3,18)*six-r1(3,23)*nine &
&          +r0(1)*three+r0(6)*six+r0(11)*three
      rrr(12)=-r7(33)+r6(25,1)*two+r6(25,2)-r5(18,1)-r5(18,2)*two-r5(18,3)*three &
&          -r5(20,8)*three+r4(12,1)+r4(12,2)*two+r4(12,3)+r4(14,10)*six+r4(14,14)*three &
&          -r3(9,7)*three-r3(9,11)*six-r3(9,15)*nine+r2(6,1)*three+r2(6,5)*six &
&          +r2(6,9)*three
      rrr(13)=-r7(34)+r6(26,1)*two+r6(26,2)-r5(19,1)-r5(19,2)*two-r5(19,3)*three-r5(19,8) &
&          -r5(21,8)+r4(13,1)+r4(13,2)*two+r4(13,3)+r4(13,10)*two+r4(15,10)*two+r4(13,14) &
&          +r4(15,14)-r3(8,7)-r3(10,7)-r3(8,11)*two-r3(10,11)*two-r3(8,15)*three &
&          -r3(10,15)*three-r3(10,27)+r2(2,1)+r2(3,1)+r2(2,5)*two+r2(3,5)*two+r2(2,9) &
&          +r2(3,9)+r2(3,25)*two+r2(3,30)-r1(3,13)-r1(3,18)*two-r1(3,23)*three+r0(1) &
&          +r0(6)*two+r0(11)
      rrr(14)=-r7(35)+r6(27,1)*two+r6(27,2)-r5(20,1)-r5(20,2)*two-r5(20,3)*three &
&          -r5(20,8)*three+r4(14,1)+r4(14,2)*two+r4(14,3)+r4(14,10)*six+r4(14,14)*three &
&          -r3(9,7)*three-r3(9,11)*six-r3(9,15)*nine+r2(6,1)*three+r2(6,5)*six &
&          +r2(6,9)*three
      rrr(15)=-r7(36)+r6(28,1)*two+r6(28,2)-r5(21,1)-r5(21,2)*two-r5(21,3)*three-r5(21,8)*six &
&          +r4(15,1)+r4(15,2)*two+r4(15,3)+r4(15,10)*p12+r4(15,14)*six-r3(10,7)*six &
&          -r3(10,11)*p12-r3(10,15)*p18-r3(10,27)*three+r2(3,1)*six+r2(3,5)*p12 &
&          +r2(3,9)*six+r2(3,25)*six+r2(3,30)*three-r1(3,13)*three-r1(3,18)*six &
&          -r1(3,23)*nine+r0(1)*three+r0(6)*six+r0(11)*three
      rxyz(1)=-r3(10,31)+r2(3,29)*two+r2(3,34)-r1(3,17)-r1(3,22)*two-r1(3,27)*three+r0(5) &
&             +r0(10)*two+r0(15)
      rxyz(2)=-r4(14,21)+r3(9,22)*two+r3(9,26)-r2(6,16)-r2(6,20)*two-r2(6,24)*three &
&             +r1(2,4)+r1(2,8)*two+r1(2,12)
      rxyz(3)=-r4(14,20)+r3(9,21)*two+r3(9,25)-r2(6,15)-r2(6,19)*two-r2(6,23)*three &
&             +r1(2,3)+r1(2,7)*two+r1(2,11)
      rxyz(4)=-r5(19,9)+r4(13,11)*two+r4(13,15)-r3(8,8)-r3(8,12)*two-r3(8,16)*three &
&             -r3(10,28)+r2(2,2)+r2(2,6)*two+r2(2,10)+r2(3,26)*two+r2(3,31)-r1(3,14) &
&             -r1(3,19)*two-r1(3,24)*three+r0(2)+r0(7)*two+r0(12)
      rxyz(5)=-r5(19,11)+r4(13,13)*two+r4(13,17)-r3(8,10)-r3(8,14)*two-r3(8,18)*three &
&             -r3(10,30)+r2(2,4)+r2(2,8)*two+r2(2,12)+r2(3,28)*two+r2(3,33)-r1(3,16) &
&             -r1(3,21)*two-r1(3,26)*three+r0(4)+r0(9)*two+r0(14)
      rxyz(6)=-r5(19,10)+r4(13,12)*two+r4(13,16)-r3(8,9)-r3(8,13)*two-r3(8,17)*three &
&             -r3(10,29)+r2(2,3)+r2(2,7)*two+r2(2,11)+r2(3,27)*two+r2(3,32)-r1(3,15) &
&             -r1(3,20)*two-r1(3,25)*three+r0(3)+r0(8)*two+r0(13)
      rxyz(7)=-r6(25,3)+r5(18,4)*two+r5(18,6)-r4(12,4)-r4(12,6)*two-r4(12,8)*three &
&             -r4(14,18)*three+r3(7,1)+r3(7,3)*two+r3(7,5)+r3(9,19)*six+r3(9,23)*three &
&             -r2(6,13)*three-r2(6,17)*six-r2(6,21)*nine+r1(2,1)*three+r1(2,5)*six &
&             +r1(2,9)*three
      rxyz(8)=-r6(25,4)+r5(18,5)*two+r5(18,7)-r4(12,5)-r4(12,7)*two-r4(12,9)*three &
&             -r4(14,19)*three+r3(7,2)+r3(7,4)*two+r3(7,6)+r3(9,20)*six+r3(9,24)*three &
&             -r2(6,14)*three-r2(6,18)*six-r2(6,22)*nine+r1(2,2)*three+r1(2,6)*six &
&             +r1(2,10)*three
      rxyz(9)=-r5(20,9)+r4(14,11)*two+r4(14,15)-r3(9,8)-r3(9,12)*two-r3(9,16)*three+r2(6,2) &
&             +r2(6,6)*two+r2(6,10)
      rxyz(10)=-r5(20,10)+r4(14,12)*two+r4(14,16)-r3(9,9)-r3(9,13)*two-r3(9,17)*three &
&              +r2(6,3)+r2(6,7)*two+r2(6,11)
      rxyz(11)=-r5(20,11)+r4(14,13)*two+r4(14,17)-r3(9,10)-r3(9,14)*two-r3(9,18)*three &
&              +r2(6,4)+r2(6,8)*two+r2(6,12)
      rxyz(12)=-r6(26,3)+r5(19,4)*two+r5(19,6)-r4(13,4)-r4(13,6)*two-r4(13,8)*three &
&             -r4(15,18)+r3(8,1)+r3(8,3)*two+r3(8,5)+r3(10,19)*two+r3(10,23)-r2(3,13) &
&             -r2(3,17)*two-r2(3,21)*three+r1(3,1)+r1(3,5)*two+r1(3,9)
      rxyz(13)=-r6(26,4)+r5(19,5)*two+r5(19,7)-r4(13,5)-r4(13,7)*two-r4(13,9)*three &
&             -r4(15,19)+r3(8,2)+r3(8,4)*two+r3(8,6)+r3(10,20)*two+r3(10,24)-r2(3,14) &
&             -r2(3,18)*two-r2(3,22)*three+r1(3,2)+r1(3,6)*two+r1(3,10)
      rxyz(14)=-r6(19,4)+r5(13,5)*two+r5(13,7)-r4(8,5)-r4(8,7)*two-r4(8,9)*three-r4(10,19) &
&             +r3(4,2)+r3(4,4)*two+r3(4,6)+r3(6,20)*two+r3(6,24)-r2(5,14)-r2(5,18)*two &
&             -r2(5,22)*three+r1(1,2)+r1(1,6)*two+r1(1,10)
      rxyz(15)=-r6(19,3)+r5(13,4)*two+r5(13,6)-r4(8,4)-r4(8,6)*two-r4(8,8)*three-r4(10,18) &
&             +r3(4,1)+r3(4,3)*two+r3(4,5)+r3(6,19)*two+r3(6,23)-r2(5,13)-r2(5,17)*two &
&             -r2(5,21)*three+r1(1,1)+r1(1,5)*two+r1(1,9)
      rxyz(16)=-r6(14,4)+r5(9,5)*two+r5(9,7)-r4(5,5)-r4(5,7)*two-r4(5,9)*three-r4(14,19) &
&             +r3(2,2)+r3(2,4)*two+r3(2,6)+r3(9,20)*two+r3(9,24)-r2(6,14)-r2(6,18)*two &
&             -r2(6,22)*three+r1(2,2)+r1(2,6)*two+r1(2,10)
      rxyz(17)=-r6(14,3)+r5(9,4)*two+r5(9,6)-r4(5,4)-r4(5,6)*two-r4(5,8)*three-r4(14,18) &
&             +r3(2,1)+r3(2,3)*two+r3(2,5)+r3(9,19)*two+r3(9,23)-r2(6,13)-r2(6,17)*two &
&             -r2(6,21)*three+r1(2,1)+r1(2,5)*two+r1(2,9)
      rxyz(18)=-r6(27,3)+r5(20,4)*two+r5(20,6)-r4(14,4)-r4(14,6)*two-r4(14,8)*three &
&             -r4(14,18)+r3(9,1)+r3(9,3)*two+r3(9,5)+r3(9,19)*two+r3(9,23)-r2(6,13) &
&             -r2(6,17)*two-r2(6,21)*three+r1(2,1)+r1(2,5)*two+r1(2,9)
      rxyz(19)=-r6(27,4)+r5(20,5)*two+r5(20,7)-r4(14,5)-r4(14,7)*two-r4(14,9)*three &
&             -r4(14,19)+r3(9,2)+r3(9,4)*two+r3(9,6)+r3(9,20)*two+r3(9,24)-r2(6,14) &
&             -r2(6,18)*two-r2(6,22)*three+r1(2,2)+r1(2,6)*two+r1(2,10)
      rxyz(20)=-r5(15,9)+r4(10,11)*two+r4(10,15)-r3(6,8)-r3(6,12)*two-r3(6,16)*three+r2(5,2) &
&              +r2(5,6)*two+r2(5,10)
      rxyz(21)=-r5(15,10)*two+r4(10,12)*four+r4(10,16)*two-r3(6,9)*two-r3(6,13)*four &
&              -r3(6,17)*six+r2(5,3)*two+r2(5,7)*four+r2(5,11)*two
      rxyz(22)=-r5(15,11)+r4(10,13)*two+r4(10,17)-r3(6,10)-r3(6,14)*two-r3(6,18)*three+r2(5,4) &
&              +r2(5,8)*two+r2(5,12)
      rxyz(23)=-r5(10,9)+r4(6,11)*two+r4(6,15)-r3(3,8)-r3(3,12)*two-r3(3,16)*three-r3(10,28) &
&              +r2(1,2)+r2(1,6)*two+r2(1,10)+r2(3,26)*two+r2(3,31)-r1(3,14)-r1(3,19)*two &
&              -r1(3,24)*three+r0(2)+r0(7)*two+r0(12)
      rxyz(24)=-r5(10,10)+r4(6,12)*two+r4(6,16)-r3(3,9)-r3(3,13)*two-r3(3,17)*three-r3(10,29) &
&              +r2(1,3)+r2(1,7)*two+r2(1,11)+r2(3,27)*two+r2(3,32)-r1(3,15)-r1(3,20)*two &
&              -r1(3,25)*three+r0(3)+r0(8)*two+r0(13)
      rxyz(25)=-r5(10,11)+r4(6,13)*two+r4(6,17)-r3(3,10)-r3(3,14)*two-r3(3,18)*three-r3(10,30) &
&              +r2(1,4)+r2(1,8)*two+r2(1,12)+r2(3,28)*two+r2(3,33)-r1(3,16)-r1(3,21)*two &
&              -r1(3,26)*three+r0(4)+r0(9)*two+r0(14)
      rxyz(26)=-r4(10,20)+r3(6,21)*two+r3(6,25)-r2(5,15)-r2(5,19)*two-r2(5,23)*three+r1(1,3) &
&               +r1(1,7)*two+r1(1,11)
      rxyz(27)=-r4(10,21)+r3(6,22)*two+r3(6,26)-r2(5,16)-r2(5,20)*two-r2(5,24)*three+r1(1,4) &
&               +r1(1,8)*two+r1(1,12)
      rxyz(28)=-r6(21,3)+r5(15,4)*two+r5(15,6)-r4(10,4)-r4(10,6)*two-r4(10,8)*three-r4(10,18) &
&               +r3(6,1)+r3(6,3)*two+r3(6,5)+r3(6,19)*two+r3(6,23)-r2(5,13)-r2(5,17)*two &
&               -r2(5,21)*three+r1(1,1)+r1(1,5)*two+r1(1,9)
      rxyz(29)=-r6(21,4)+r5(15,5)*two+r5(15,7)-r4(10,5)-r4(10,7)*two-r4(10,9)*three-r4(10,19) &
&               +r3(6,2)+r3(6,4)*two+r3(6,6)+r3(6,20)*two+r3(6,24)-r2(5,14)-r2(5,18)*two &
&               -r2(5,22)*three+r1(1,2)+r1(1,6)*two+r1(1,10)
      rxyz(30)=-r6(15,3)+r5(10,4)*two+r5(10,6)-r4(6,4)-r4(6,6)*two-r4(6,8)*three-r4(15,18) &
&              +r3(3,1)+r3(3,3)*two+r3(3,5)+r3(10,19)*two+r3(10,23)-r2(3,13)-r2(3,17)*two &
&              -r2(3,21)*three+r1(3,1)+r1(3,5)*two+r1(3,9)
      rxyz(31)=-r6(15,4)+r5(10,5)*two+r5(10,7)-r4(6,5)-r4(6,7)*two-r4(6,9)*three-r4(15,19) &
&              +r3(3,2)+r3(3,4)*two+r3(3,6)+r3(10,20)*two+r3(10,24)-r2(3,14)-r2(3,18)*two &
&              -r2(3,22)*three+r1(3,2)+r1(3,6)*two+r1(3,10)
      rxyz(32)=-r5(21,9)+r4(15,11)*two+r4(15,15)-r3(10,8)-r3(10,12)*two-r3(10,16)*three &
&              -r3(10,28)+r2(3,2)+r2(3,6)*two+r2(3,10)+r2(3,26)*two+r2(3,31)-r1(3,14) &
&              -r1(3,19)*two-r1(3,24)*three+r0(2)+r0(7)*two+r0(12)
      rxyz(33)=-r5(21,10)+r4(15,12)*two+r4(15,16)-r3(10,9)-r3(10,13)*two-r3(10,17)*three &
&              -r3(10,29)+r2(3,3)+r2(3,7)*two+r2(3,11)+r2(3,27)*two+r2(3,32)-r1(3,15) &
&              -r1(3,20)*two-r1(3,25)*three+r0(3)+r0(8)*two+r0(13)
      rxyz(34)=-r5(21,11)+r4(15,13)*two+r4(15,17)-r3(10,10)-r3(10,14)*two-r3(10,18)*three &
&              -r3(10,30)+r2(3,4)+r2(3,8)*two+r2(3,12)+r2(3,28)*two+r2(3,33)-r1(3,16) &
&              -r1(3,21)*two-r1(3,26)*three+r0(4)+r0(9)*two+r0(14)
      rxyz(35)=-r4(15,20)+r3(10,21)*two+r3(10,25)-r2(3,15)-r2(3,19)*two-r2(3,23)*three+r1(3,3) &
&              +r1(3,7)*two+r1(3,11)
      rxyz(36)=-r4(15,21)+r3(10,22)*two+r3(10,26)-r2(3,16)-r2(3,20)*two-r2(3,24)*three+r1(3,4) &
&              +r1(3,8)*two+r1(3,12)
      rxyz(37)=-r6(10,3)+r5(6,4)*two+r5(6,6)-r4(3,4)-r4(3,6)*two-r4(3,8)*three-r4(10,18)*three &
&              +r3(1,1)+r3(1,3)*two+r3(1,5)+r3(6,19)*six+r3(6,23)*three-r2(5,13)*three &
&              -r2(5,17)*six-r2(5,21)*nine+r1(1,1)*three+r1(1,5)*six+r1(1,9)*three
      rxyz(38)=-r6(10,4)+r5(6,5)*two+r5(6,7)-r4(3,5)-r4(3,7)*two-r4(3,9)*three-r4(10,19)*three &
&              +r3(1,2)+r3(1,4)*two+r3(1,6)+r3(6,20)*six+r3(6,24)*three-r2(5,14)*three &
&              -r2(5,18)*six-r2(5,22)*nine+r1(1,2)*three+r1(1,6)*six+r1(1,10)*three
      rxyz(39)=-r6(20,3)+r5(14,4)*two+r5(14,6)-r4(9,4)-r4(9,6)*two-r4(9,8)*three+r3(5,1) &
&              +r3(5,3)*two+r3(5,5)
      rxyz(40)=-r6(20,4)+r5(14,5)*two+r5(14,7)-r4(9,5)-r4(9,7)*two-r4(9,9)*three+r3(5,2) &
&              +r3(5,4)*two+r3(5,6)
      rxyz(41)=-r5(14,9)+r4(9,11)*two+r4(9,15)-r3(5,8)-r3(5,12)*two-r3(5,16)*three+r2(4,2) &
&              +r2(4,6)*two+r2(4,10)
      rxyz(42)=-r5(14,10)+r4(9,12)*two+r4(9,16)-r3(5,9)-r3(5,13)*two-r3(5,17)*three &
&              +r2(4,3)+r2(4,7)*two+r2(4,11)
      rxyz(43)=-r5(14,11)+r4(9,13)*two+r4(9,17)-r3(5,10)-r3(5,14)*two-r3(5,18)*three+r2(4,4) &
&              +r2(4,8)*two+r2(4,12)
      rxyz(44)=-r6(28,3)+r5(21,4)*two+r5(21,6)-r4(15,4)-r4(15,6)*two-r4(15,8)*three &
&              -r4(15,18)*three+r3(10,1)+r3(10,3)*two+r3(10,5)+r3(10,19)*six+r3(10,23)*three &
&              -r2(3,13)*three-r2(3,17)*six-r2(3,21)*nine+r1(3,1)*three+r1(3,5)*six+r1(3,9)*three
      rxyz(45)=-r6(28,4)+r5(21,5)*two+r5(21,7)-r4(15,5)-r4(15,7)*two-r4(15,9)*three &
&              -r4(15,19)*three+r3(10,2)+r3(10,4)*two+r3(10,6)+r3(10,20)*six+r3(10,24)*three &
&              -r2(3,14)*three-r2(3,18)*six-r2(3,22)*nine+r1(3,2)*three+r1(3,6)*six+r1(3,10)*three
      call int2dddp2(eri,rrr,rxyz,qx,qz,3,3)
!
      rrr( 1)=-r7(5)+r6(2,2)-r5(5,8)*six+r4(2,14)*six-r3(5,27)*three+r2(4,30)*three
      rrr( 2)=-r7(8)+r6(4,2)-r5(8,8)*three+r4(4,14)*three
      rrr( 3)=-r7(9)+r6(5,2)-r5(9,8)*three+r4(5,14)*three
      rrr( 4)=-r7(12)+r6(7,2)-r5(5,8)-r5(12,8)+r4(2,14)+r4(7,14)-r3(5,27)+r2(4,30)
      rrr( 5)=-r7(13)+r6(8,2)-r5(13,8)+r4(8,14)
      rrr( 6)=-r7(14)+r6(9,2)-r5(5,8)-r5(14,8)+r4(2,14)+r4(9,14)-r3(5,27)+r2(4,30)
      rrr( 7)=-r7(17)+r6(11,2)-r5(8,8)*three+r4(4,14)*three
      rrr( 8)=-r7(18)+r6(12,2)-r5(9,8)+r4(5,14)
      rrr( 9)=-r7(19)+r6(13,2)-r5(8,8)+r4(4,14)
      rrr(10)=-r7(20)+r6(14,2)-r5(9,8)*three+r4(5,14)*three
      rrr(11)=-r7(23)+r6(16,2)-r5(12,8)*six+r4(7,14)*six-r3(5,27)*three+r2(4,30)*three
      rrr(12)=-r7(24)+r6(17,2)-r5(13,8)*three+r4(8,14)*three
      rrr(13)=-r7(25)+r6(18,2)-r5(12,8)-r5(14,8)+r4(7,14)+r4(9,14)-r3(5,27)+r2(4,30)
      rrr(14)=-r7(26)+r6(19,2)-r5(13,8)*three+r4(8,14)*three
      rrr(15)=-r7(27)+r6(20,2)-r5(14,8)*six+r4(9,14)*six-r3(5,27)*three+r2(4,30)*three
      rxyz(1)=-r3(5,31)+r2(4,34)
      rxyz(2)=-r4(8,21)+r3(4,26)
      rxyz(3)=-r4(8,20)+r3(4,25)
      rxyz(4)=-r5(12,9)+r4(7,15)-r3(5,28)+r2(4,31)
      rxyz(5)=-r5(12,11)+r4(7,17)-r3(5,30)+r2(4,33)
      rxyz(6)=-r5(12,10)+r4(7,16)-r3(5,29)+r2(4,32)
      rxyz(7)=-r6(17,3)+r5(11,6)-r4(8,18)*three+r3(4,23)*three
      rxyz(8)=-r6(17,4)+r5(11,7)-r4(8,19)*three+r3(4,24)*three
      rxyz(9)=-r5(13,9)+r4(8,15)
      rxyz(10)=-r5(13,10)+r4(8,16)
      rxyz(11)=-r5(13,11)+r4(8,17)
      rxyz(12)=-r6(18,3)+r5(12,6)-r4(9,18)+r3(5,23)
      rxyz(13)=-r6(18,4)+r5(12,7)-r4(9,19)+r3(5,24)
      rxyz(14)=-r6(12,4)+r5(7,7)-r4(5,19)+r3(2,24)
      rxyz(15)=-r6(12,3)+r5(7,6)-r4(5,18)+r3(2,23)
      rxyz(16)=-r6(8,4)+r5(4,7)-r4(8,19)+r3(4,24)
      rxyz(17)=-r6(8,3)+r5(4,6)-r4(8,18)+r3(4,23)
      rxyz(18)=-r6(19,3)+r5(13,6)-r4(8,18)+r3(4,23)
      rxyz(19)=-r6(19,4)+r5(13,7)-r4(8,19)+r3(4,24)
      rxyz(20)=-r5(9,9)+r4(5,15)
      rxyz(21)=-r5(9,10)*two+r4(5,16)*two
      rxyz(22)=-r5(9,11)+r4(5,17)
      rxyz(23)=-r5(5,9)+r4(2,15)-r3(5,28)+r2(4,31)
      rxyz(24)=-r5(5,10)+r4(2,16)-r3(5,29)+r2(4,32)
      rxyz(25)=-r5(5,11)+r4(2,17)-r3(5,30)+r2(4,33)
      rxyz(26)=-r4(5,20)+r3(2,25)
      rxyz(27)=-r4(5,21)+r3(2,26)
      rxyz(28)=-r6(14,3)+r5(9,6)-r4(5,18)+r3(2,23)
      rxyz(29)=-r6(14,4)+r5(9,7)-r4(5,19)+r3(2,24)
      rxyz(30)=-r6(9,3)+r5(5,6)-r4(9,18)+r3(5,23)
      rxyz(31)=-r6(9,4)+r5(5,7)-r4(9,19)+r3(5,24)
      rxyz(32)=-r5(14,9)+r4(9,15)-r3(5,28)+r2(4,31)
      rxyz(33)=-r5(14,10)+r4(9,16)-r3(5,29)+r2(4,32)
      rxyz(34)=-r5(14,11)+r4(9,17)-r3(5,30)+r2(4,33)
      rxyz(35)=-r4(9,20)+r3(5,25)
      rxyz(36)=-r4(9,21)+r3(5,26)
      rxyz(37)=-r6(5,3)+r5(2,6)-r4(5,18)*three+r3(2,23)*three
      rxyz(38)=-r6(5,4)+r5(2,7)-r4(5,19)*three+r3(2,24)*three
      rxyz(39)=-r6(13,3)+r5(8,6)
      rxyz(40)=-r6(13,4)+r5(8,7)
      rxyz(41)=-r5(8,9)+r4(4,15)
      rxyz(42)=-r5(8,10)+r4(4,16)
      rxyz(43)=-r5(8,11)+r4(4,17)
      rxyz(44)=-r6(20,3)+r5(14,6)-r4(9,18)*three+r3(5,23)*three
      rxyz(45)=-r6(20,4)+r5(14,7)-r4(9,19)*three+r3(5,24)*three
      call int2dddp2(eri,rrr,rxyz,qx,qz,4,3)
!
      rrr( 1)=-r7(6)+r6(3,1)+r6(3,2)-r5(1,2)-r5(1,3)-r5(6,8)*six+r4(3,10)*six+r4(3,14)*six &
&          -r3(1,11)*six-r3(1,15)*six-r3(6,27)*three+r2(5,25)*three+r2(5,30)*three &
&          -r1(1,18)*three-r1(1,23)*three
      rrr( 2)=-r7(9)+r6(5,1)+r6(5,2)-r5(2,2)-r5(2,3)-r5(9,8)*three+r4(5,10)*three &
&          +r4(5,14)*three-r3(2,11)*three-r3(2,15)*three
      rrr( 3)=-r7(10)+r6(6,1)+r6(6,2)-r5(3,2)-r5(3,3)-r5(10,8)*three+r4(6,10)*three &
&          +r4(6,14)*three-r3(3,11)*three-r3(3,15)*three
      rrr( 4)=-r7(13)+r6(8,1)+r6(8,2)-r5(4,2)-r5(4,3)-r5(6,8)-r5(13,8)+r4(3,10)+r4(8,10) &
&          +r4(3,14)+r4(8,14)-r3(1,11)-r3(4,11)-r3(1,15)-r3(4,15)-r3(6,27)+r2(5,25) &
&          +r2(5,30)-r1(1,18)-r1(1,23)
      rrr( 5)=-r7(14)+r6(9,1)+r6(9,2)-r5(5,2)-r5(5,3)-r5(14,8)+r4(9,10)+r4(9,14)-r3(5,11) &
&          -r3(5,15)
      rrr( 6)=-r7(15)+r6(10,1)+r6(10,2)-r5(6,2)-r5(6,3)-r5(6,8)-r5(15,8)+r4(3,10)+r4(10,10) &
&          +r4(3,14)+r4(10,14)-r3(1,11)-r3(6,11)-r3(1,15)-r3(6,15)-r3(6,27)+r2(5,25) &
&          +r2(5,30)-r1(1,18)-r1(1,23)
      rrr( 7)=-r7(18)+r6(12,1)+r6(12,2)-r5(7,2)-r5(7,3)-r5(9,8)*three+r4(5,10)*three &
&          +r4(5,14)*three-r3(2,11)*three-r3(2,15)*three
      rrr( 8)=-r7(19)+r6(13,1)+r6(13,2)-r5(8,2)-r5(8,3)-r5(10,8)+r4(6,10)+r4(6,14)-r3(3,11) &
&          -r3(3,15)
      rrr( 9)=-r7(20)+r6(14,1)+r6(14,2)-r5(9,2)-r5(9,3)-r5(9,8)+r4(5,10)+r4(5,14)-r3(2,11) &
&          -r3(2,15)
      rrr(10)=-r7(21)+r6(15,1)+r6(15,2)-r5(10,2)-r5(10,3)-r5(10,8)*three+r4(6,10)*three &
&          +r4(6,14)*three-r3(3,11)*three-r3(3,15)*three
      rrr(11)=-r7(24)+r6(17,1)+r6(17,2)-r5(11,2)-r5(11,3)-r5(13,8)*six+r4(8,10)*six &
&          +r4(8,14)*six-r3(4,11)*six-r3(4,15)*six-r3(6,27)*three+r2(5,25)*three &
&          +r2(5,30)*three-r1(1,18)*three-r1(1,23)*three
      rrr(12)=-r7(25)+r6(18,1)+r6(18,2)-r5(12,2)-r5(12,3)-r5(14,8)*three+r4(9,10)*three &
&          +r4(9,14)*three-r3(5,11)*three-r3(5,15)*three
      rrr(13)=-r7(26)+r6(19,1)+r6(19,2)-r5(13,2)-r5(13,3)-r5(13,8)-r5(15,8)+r4(8,10) &
&          +r4(10,10)+r4(8,14)+r4(10,14)-r3(4,11)-r3(6,11)-r3(4,15)-r3(6,15)-r3(6,27) &
&          +r2(5,25)+r2(5,30)-r1(1,18)-r1(1,23)
      rrr(14)=-r7(27)+r6(20,1)+r6(20,2)-r5(14,2)-r5(14,3)-r5(14,8)*three+r4(9,10)*three &
&          +r4(9,14)*three-r3(5,11)*three-r3(5,15)*three
      rrr(15)=-r7(28)+r6(21,1)+r6(21,2)-r5(15,2)-r5(15,3)-r5(15,8)*six+r4(10,10)*six &
&          +r4(10,14)*six-r3(6,11)*six-r3(6,15)*six-r3(6,27)*three+r2(5,25)*three &
&          +r2(5,30)*three-r1(1,18)*three-r1(1,23)*three
      rxyz(1)=-r3(6,31)+r2(5,29)+r2(5,34)-r1(1,22)-r1(1,27)
      rxyz(2)=-r4(9,21)+r3(5,22)+r3(5,26)-r2(4,20)-r2(4,24)
      rxyz(3)=-r4(9,20)+r3(5,21)+r3(5,25)-r2(4,19)-r2(4,23)
      rxyz(4)=-r5(13,9)+r4(8,11)+r4(8,15)-r3(4,12)-r3(4,16)-r3(6,28)+r2(5,26)+r2(5,31) &
&             -r1(1,19)-r1(1,24)
      rxyz(5)=-r5(13,11)+r4(8,13)+r4(8,17)-r3(4,14)-r3(4,18)-r3(6,30)+r2(5,28)+r2(5,33) &
&             -r1(1,21)-r1(1,26)
      rxyz(6)=-r5(13,10)+r4(8,12)+r4(8,16)-r3(4,13)-r3(4,17)-r3(6,29)+r2(5,27)+r2(5,32) &
&             -r1(1,20)-r1(1,25)
      rxyz(7)=-r6(18,3)+r5(12,4)+r5(12,6)-r4(7,6)-r4(7,8)-r4(9,18)*three+r3(5,19)*three &
&             +r3(5,23)*three-r2(4,17)*three-r2(4,21)*three
      rxyz(8)=-r6(18,4)+r5(12,5)+r5(12,7)-r4(7,7)-r4(7,9)-r4(9,19)*three+r3(5,20)*three &
&             +r3(5,24)*three-r2(4,18)*three-r2(4,22)*three
      rxyz(9)=-r5(14,9)+r4(9,11)+r4(9,15)-r3(5,12)-r3(5,16)
      rxyz(10)=-r5(14,10)+r4(9,12)+r4(9,16)-r3(5,13)-r3(5,17)
      rxyz(11)=-r5(14,11)+r4(9,13)+r4(9,17)-r3(5,14)-r3(5,18)
      rxyz(12)=-r6(19,3)+r5(13,4)+r5(13,6)-r4(8,6)-r4(8,8)-r4(10,18)+r3(6,19)+r3(6,23) &
&             -r2(5,17)-r2(5,21)
      rxyz(13)=-r6(19,4)+r5(13,5)+r5(13,7)-r4(8,7)-r4(8,9)-r4(10,19)+r3(6,20)+r3(6,24) &
&             -r2(5,18)-r2(5,22)
      rxyz(14)=-r6(13,4)+r5(8,5)+r5(8,7)-r4(4,7)-r4(4,9)-r4(6,19)+r3(3,20)+r3(3,24) &
&             -r2(1,18)-r2(1,22)
      rxyz(15)=-r6(13,3)+r5(8,4)+r5(8,6)-r4(4,6)-r4(4,8)-r4(6,18)+r3(3,19)+r3(3,23) &
&             -r2(1,17)-r2(1,21)
      rxyz(16)=-r6(9,4)+r5(5,5)+r5(5,7)-r4(2,7)-r4(2,9)-r4(9,19)+r3(5,20)+r3(5,24) &
&             -r2(4,18)-r2(4,22)
      rxyz(17)=-r6(9,3)+r5(5,4)+r5(5,6)-r4(2,6)-r4(2,8)-r4(9,18)+r3(5,19)+r3(5,23) &
&             -r2(4,17)-r2(4,21)
      rxyz(18)=-r6(20,3)+r5(14,4)+r5(14,6)-r4(9,6)-r4(9,8)-r4(9,18)+r3(5,19)+r3(5,23) &
&             -r2(4,17)-r2(4,21)
      rxyz(19)=-r6(20,4)+r5(14,5)+r5(14,7)-r4(9,7)-r4(9,9)-r4(9,19)+r3(5,20)+r3(5,24) &
&             -r2(4,18)-r2(4,22)
      rxyz(20)=-r5(10,9)+r4(6,11)+r4(6,15)-r3(3,12)-r3(3,16)
      rxyz(21)=-r5(10,10)*two+r4(6,12)*two+r4(6,16)*two-r3(3,13)*two-r3(3,17)*two
      rxyz(22)=-r5(10,11)+r4(6,13)+r4(6,17)-r3(3,14)-r3(3,18)
      rxyz(23)=-r5(6,9)+r4(3,11)+r4(3,15)-r3(1,12)-r3(1,16)-r3(6,28)+r2(5,26)+r2(5,31) &
&              -r1(1,19)-r1(1,24)
      rxyz(24)=-r5(6,10)+r4(3,12)+r4(3,16)-r3(1,13)-r3(1,17)-r3(6,29)+r2(5,27)+r2(5,32) &
&              -r1(1,20)-r1(1,25)
      rxyz(25)=-r5(6,11)+r4(3,13)+r4(3,17)-r3(1,14)-r3(1,18)-r3(6,30)+r2(5,28)+r2(5,33) &
&              -r1(1,21)-r1(1,26)
      rxyz(26)=-r4(6,20)+r3(3,21)+r3(3,25)-r2(1,19)-r2(1,23)
      rxyz(27)=-r4(6,21)+r3(3,22)+r3(3,26)-r2(1,20)-r2(1,24)
      rxyz(28)=-r6(15,3)+r5(10,4)+r5(10,6)-r4(6,6)-r4(6,8)-r4(6,18)+r3(3,19)+r3(3,23)-r2(1,17) &
&              -r2(1,21)
      rxyz(29)=-r6(15,4)+r5(10,5)+r5(10,7)-r4(6,7)-r4(6,9)-r4(6,19)+r3(3,20)+r3(3,24)-r2(1,18) &
&              -r2(1,22)
      rxyz(30)=-r6(10,3)+r5(6,4)+r5(6,6)-r4(3,6)-r4(3,8)-r4(10,18)+r3(6,19)+r3(6,23)-r2(5,17) &
&              -r2(5,21)
      rxyz(31)=-r6(10,4)+r5(6,5)+r5(6,7)-r4(3,7)-r4(3,9)-r4(10,19)+r3(6,20)+r3(6,24)-r2(5,18) &
&              -r2(5,22)
      rxyz(32)=-r5(15,9)+r4(10,11)+r4(10,15)-r3(6,12)-r3(6,16)-r3(6,28)+r2(5,26)+r2(5,31) &
&              -r1(1,19)-r1(1,24)
      rxyz(33)=-r5(15,10)+r4(10,12)+r4(10,16)-r3(6,13)-r3(6,17)-r3(6,29)+r2(5,27)+r2(5,32) &
&              -r1(1,20)-r1(1,25)
      rxyz(34)=-r5(15,11)+r4(10,13)+r4(10,17)-r3(6,14)-r3(6,18)-r3(6,30)+r2(5,28)+r2(5,33) &
&              -r1(1,21)-r1(1,26)
      rxyz(35)=-r4(10,20)+r3(6,21)+r3(6,25)-r2(5,19)-r2(5,23)
      rxyz(36)=-r4(10,21)+r3(6,22)+r3(6,26)-r2(5,20)-r2(5,24)
      rxyz(37)=-r6(6,3)+r5(3,4)+r5(3,6)-r4(1,6)-r4(1,8)-r4(6,18)*three+r3(3,19)*three &
&              +r3(3,23)*three-r2(1,17)*three-r2(1,21)*three
      rxyz(38)=-r6(6,4)+r5(3,5)+r5(3,7)-r4(1,7)-r4(1,9)-r4(6,19)*three+r3(3,20)*three &
&              +r3(3,24)*three-r2(1,18)*three-r2(1,22)*three
      rxyz(39)=-r6(14,3)+r5(9,4)+r5(9,6)-r4(5,6)-r4(5,8)
      rxyz(40)=-r6(14,4)+r5(9,5)+r5(9,7)-r4(5,7)-r4(5,9)
      rxyz(41)=-r5(9,9)+r4(5,11)+r4(5,15)-r3(2,12)-r3(2,16)
      rxyz(42)=-r5(9,10)+r4(5,12)+r4(5,16)-r3(2,13)-r3(2,17)
      rxyz(43)=-r5(9,11)+r4(5,13)+r4(5,17)-r3(2,14)-r3(2,18)
      rxyz(44)=-r6(21,3)+r5(15,4)+r5(15,6)-r4(10,6)-r4(10,8)-r4(10,18)*three+r3(6,19)*three &
&              +r3(6,23)*three-r2(5,17)*three-r2(5,21)*three
      rxyz(45)=-r6(21,4)+r5(15,5)+r5(15,7)-r4(10,7)-r4(10,9)-r4(10,19)*three+r3(6,20)*three &
&              +r3(6,24)*three-r2(5,18)*three-r2(5,22)*three
      call int2dddp2(eri,rrr,rxyz,qx,qz,5,3)
!
      rrr( 1)=-r7(9)+r6(5,1)+r6(5,2)-r5(2,2)-r5(2,3)-r5(9,8)*six+r4(5,10)*six+r4(5,14)*six &
&          -r3(2,11)*six-r3(2,15)*six-r3(9,27)*three+r2(6,25)*three+r2(6,30)*three &
&          -r1(2,18)*three-r1(2,23)*three
      rrr( 2)=-r7(13)+r6(8,1)+r6(8,2)-r5(4,2)-r5(4,3)-r5(13,8)*three+r4(8,10)*three &
&          +r4(8,14)*three-r3(4,11)*three-r3(4,15)*three
      rrr( 3)=-r7(14)+r6(9,1)+r6(9,2)-r5(5,2)-r5(5,3)-r5(14,8)*three+r4(9,10)*three &
&          +r4(9,14)*three-r3(5,11)*three-r3(5,15)*three
      rrr( 4)=-r7(18)+r6(12,1)+r6(12,2)-r5(7,2)-r5(7,3)-r5(9,8)-r5(18,8)+r4(5,10)+r4(12,10) &
&          +r4(5,14)+r4(12,14)-r3(2,11)-r3(7,11)-r3(2,15)-r3(7,15)-r3(9,27)+r2(6,25) &
&          +r2(6,30)-r1(2,18)-r1(2,23)
      rrr( 5)=-r7(19)+r6(13,1)+r6(13,2)-r5(8,2)-r5(8,3)-r5(19,8)+r4(13,10)+r4(13,14)-r3(8,11) &
&          -r3(8,15)
      rrr( 6)=-r7(20)+r6(14,1)+r6(14,2)-r5(9,2)-r5(9,3)-r5(9,8)-r5(20,8)+r4(5,10)+r4(14,10) &
&          +r4(5,14)+r4(14,14)-r3(2,11)-r3(9,11)-r3(2,15)-r3(9,15)-r3(9,27)+r2(6,25) &
&          +r2(6,30)-r1(2,18)-r1(2,23)
      rrr( 7)=-r7(24)+r6(17,1)+r6(17,2)-r5(11,2)-r5(11,3)-r5(13,8)*three+r4(8,10)*three &
&          +r4(8,14)*three-r3(4,11)*three-r3(4,15)*three
      rrr( 8)=-r7(25)+r6(18,1)+r6(18,2)-r5(12,2)-r5(12,3)-r5(14,8)+r4(9,10)+r4(9,14)-r3(5,11) &
&          -r3(5,15)
      rrr( 9)=-r7(26)+r6(19,1)+r6(19,2)-r5(13,2)-r5(13,3)-r5(13,8)+r4(8,10)+r4(8,14)-r3(4,11) &
&          -r3(4,15)
      rrr(10)=-r7(27)+r6(20,1)+r6(20,2)-r5(14,2)-r5(14,3)-r5(14,8)*three+r4(9,10)*three &
&          +r4(9,14)*three-r3(5,11)*three-r3(5,15)*three
      rrr(11)=-r7(31)+r6(23,1)+r6(23,2)-r5(16,2)-r5(16,3)-r5(18,8)*six+r4(12,10)*six &
&          +r4(12,14)*six-r3(7,11)*six-r3(7,15)*six-r3(9,27)*three+r2(6,25)*three &
&          +r2(6,30)*three-r1(2,18)*three-r1(2,23)*three
      rrr(12)=-r7(32)+r6(24,1)+r6(24,2)-r5(17,2)-r5(17,3)-r5(19,8)*three+r4(13,10)*three &
&          +r4(13,14)*three-r3(8,11)*three-r3(8,15)*three
      rrr(13)=-r7(33)+r6(25,1)+r6(25,2)-r5(18,2)-r5(18,3)-r5(18,8)-r5(20,8)+r4(12,10) &
&          +r4(14,10)+r4(12,14)+r4(14,14)-r3(7,11)-r3(9,11)-r3(7,15)-r3(9,15)-r3(9,27) &
&          +r2(6,25)+r2(6,30)-r1(2,18)-r1(2,23)
      rrr(14)=-r7(34)+r6(26,1)+r6(26,2)-r5(19,2)-r5(19,3)-r5(19,8)*three+r4(13,10)*three &
&          +r4(13,14)*three-r3(8,11)*three-r3(8,15)*three
      rrr(15)=-r7(35)+r6(27,1)+r6(27,2)-r5(20,2)-r5(20,3)-r5(20,8)*six+r4(14,10)*six &
&          +r4(14,14)*six-r3(9,11)*six-r3(9,15)*six-r3(9,27)*three+r2(6,25)*three &
&          +r2(6,30)*three-r1(2,18)*three-r1(2,23)*three
      rxyz(1)=-r3(9,31)+r2(6,29)+r2(6,34)-r1(2,22)-r1(2,27)
      rxyz(2)=-r4(13,21)+r3(8,22)+r3(8,26)-r2(2,20)-r2(2,24)
      rxyz(3)=-r4(13,20)+r3(8,21)+r3(8,25)-r2(2,19)-r2(2,23)
      rxyz(4)=-r5(18,9)+r4(12,11)+r4(12,15)-r3(7,12)-r3(7,16)-r3(9,28)+r2(6,26)+r2(6,31) &
&             -r1(2,19)-r1(2,24)
      rxyz(5)=-r5(18,11)+r4(12,13)+r4(12,17)-r3(7,14)-r3(7,18)-r3(9,30)+r2(6,28)+r2(6,33) &
&             -r1(2,21)-r1(2,26)
      rxyz(6)=-r5(18,10)+r4(12,12)+r4(12,16)-r3(7,13)-r3(7,17)-r3(9,29)+r2(6,27)+r2(6,32) &
&             -r1(2,20)-r1(2,25)
      rxyz(7)=-r6(24,3)+r5(17,4)+r5(17,6)-r4(11,6)-r4(11,8)-r4(13,18)*three+r3(8,19)*three &
&             +r3(8,23)*three-r2(2,17)*three-r2(2,21)*three
      rxyz(8)=-r6(24,4)+r5(17,5)+r5(17,7)-r4(11,7)-r4(11,9)-r4(13,19)*three+r3(8,20)*three &
&             +r3(8,24)*three-r2(2,18)*three-r2(2,22)*three
      rxyz(9)=-r5(19,9)+r4(13,11)+r4(13,15)-r3(8,12)-r3(8,16)
      rxyz(10)=-r5(19,10)+r4(13,12)+r4(13,16)-r3(8,13)-r3(8,17)
      rxyz(11)=-r5(19,11)+r4(13,13)+r4(13,17)-r3(8,14)-r3(8,18)
      rxyz(12)=-r6(25,3)+r5(18,4)+r5(18,6)-r4(12,6)-r4(12,8)-r4(14,18)+r3(9,19)+r3(9,23) &
&             -r2(6,17)-r2(6,21)
      rxyz(13)=-r6(25,4)+r5(18,5)+r5(18,7)-r4(12,7)-r4(12,9)-r4(14,19)+r3(9,20)+r3(9,24) &
&             -r2(6,18)-r2(6,22)
      rxyz(14)=-r6(18,4)+r5(12,5)+r5(12,7)-r4(7,7)-r4(7,9)-r4(9,19)+r3(5,20)+r3(5,24) &
&             -r2(4,18)-r2(4,22)
      rxyz(15)=-r6(18,3)+r5(12,4)+r5(12,6)-r4(7,6)-r4(7,8)-r4(9,18)+r3(5,19)+r3(5,23) &
&             -r2(4,17)-r2(4,21)
      rxyz(16)=-r6(13,4)+r5(8,5)+r5(8,7)-r4(4,7)-r4(4,9)-r4(13,19)+r3(8,20)+r3(8,24) &
&             -r2(2,18)-r2(2,22)
      rxyz(17)=-r6(13,3)+r5(8,4)+r5(8,6)-r4(4,6)-r4(4,8)-r4(13,18)+r3(8,19)+r3(8,23) &
&             -r2(2,17)-r2(2,21)
      rxyz(18)=-r6(26,3)+r5(19,4)+r5(19,6)-r4(13,6)-r4(13,8)-r4(13,18)+r3(8,19)+r3(8,23) &
&             -r2(2,17)-r2(2,21)
      rxyz(19)=-r6(26,4)+r5(19,5)+r5(19,7)-r4(13,7)-r4(13,9)-r4(13,19)+r3(8,20)+r3(8,24) &
&             -r2(2,18)-r2(2,22)
      rxyz(20)=-r5(14,9)+r4(9,11)+r4(9,15)-r3(5,12)-r3(5,16)
      rxyz(21)=-r5(14,10)*two+r4(9,12)*two+r4(9,16)*two-r3(5,13)*two-r3(5,17)*two
      rxyz(22)=-r5(14,11)+r4(9,13)+r4(9,17)-r3(5,14)-r3(5,18)
      rxyz(23)=-r5(9,9)+r4(5,11)+r4(5,15)-r3(2,12)-r3(2,16)-r3(9,28)+r2(6,26)+r2(6,31) &
&              -r1(2,19)-r1(2,24)
      rxyz(24)=-r5(9,10)+r4(5,12)+r4(5,16)-r3(2,13)-r3(2,17)-r3(9,29)+r2(6,27)+r2(6,32) &
&              -r1(2,20)-r1(2,25)
      rxyz(25)=-r5(9,11)+r4(5,13)+r4(5,17)-r3(2,14)-r3(2,18)-r3(9,30)+r2(6,28)+r2(6,33) &
&              -r1(2,21)-r1(2,26)
      rxyz(26)=-r4(9,20)+r3(5,21)+r3(5,25)-r2(4,19)-r2(4,23)
      rxyz(27)=-r4(9,21)+r3(5,22)+r3(5,26)-r2(4,20)-r2(4,24)
      rxyz(28)=-r6(20,3)+r5(14,4)+r5(14,6)-r4(9,6)-r4(9,8)-r4(9,18)+r3(5,19)+r3(5,23)-r2(4,17) &
&              -r2(4,21)
      rxyz(29)=-r6(20,4)+r5(14,5)+r5(14,7)-r4(9,7)-r4(9,9)-r4(9,19)+r3(5,20)+r3(5,24)-r2(4,18) &
&              -r2(4,22)
      rxyz(30)=-r6(14,3)+r5(9,4)+r5(9,6)-r4(5,6)-r4(5,8)-r4(14,18)+r3(9,19)+r3(9,23)-r2(6,17) &
&              -r2(6,21)
      rxyz(31)=-r6(14,4)+r5(9,5)+r5(9,7)-r4(5,7)-r4(5,9)-r4(14,19)+r3(9,20)+r3(9,24)-r2(6,18) &
&              -r2(6,22)
      rxyz(32)=-r5(20,9)+r4(14,11)+r4(14,15)-r3(9,12)-r3(9,16)-r3(9,28)+r2(6,26)+r2(6,31) &
&              -r1(2,19)-r1(2,24)
      rxyz(33)=-r5(20,10)+r4(14,12)+r4(14,16)-r3(9,13)-r3(9,17)-r3(9,29)+r2(6,27)+r2(6,32) &
&              -r1(2,20)-r1(2,25)
      rxyz(34)=-r5(20,11)+r4(14,13)+r4(14,17)-r3(9,14)-r3(9,18)-r3(9,30)+r2(6,28)+r2(6,33) &
&              -r1(2,21)-r1(2,26)
      rxyz(35)=-r4(14,20)+r3(9,21)+r3(9,25)-r2(6,19)-r2(6,23)
      rxyz(36)=-r4(14,21)+r3(9,22)+r3(9,26)-r2(6,20)-r2(6,24)
      rxyz(37)=-r6(9,3)+r5(5,4)+r5(5,6)-r4(2,6)-r4(2,8)-r4(9,18)*three+r3(5,19)*three &
&              +r3(5,23)*three-r2(4,17)*three-r2(4,21)*three
      rxyz(38)=-r6(9,4)+r5(5,5)+r5(5,7)-r4(2,7)-r4(2,9)-r4(9,19)*three+r3(5,20)*three &
&              +r3(5,24)*three-r2(4,18)*three-r2(4,22)*three
      rxyz(39)=-r6(19,3)+r5(13,4)+r5(13,6)-r4(8,6)-r4(8,8)
      rxyz(40)=-r6(19,4)+r5(13,5)+r5(13,7)-r4(8,7)-r4(8,9)
      rxyz(41)=-r5(13,9)+r4(8,11)+r4(8,15)-r3(4,12)-r3(4,16)
      rxyz(42)=-r5(13,10)+r4(8,12)+r4(8,16)-r3(4,13)-r3(4,17)
      rxyz(43)=-r5(13,11)+r4(8,13)+r4(8,17)-r3(4,14)-r3(4,18)
      rxyz(44)=-r6(27,3)+r5(20,4)+r5(20,6)-r4(14,6)-r4(14,8)-r4(14,18)*three+r3(9,19)*three &
&              +r3(9,23)*three-r2(6,17)*three-r2(6,21)*three
      rxyz(45)=-r6(27,4)+r5(20,5)+r5(20,7)-r4(14,7)-r4(14,9)-r4(14,19)*three+r3(9,20)*three &
&               +r3(9,24)*three-r2(6,18)*three-r2(6,22)*three
      call int2dddp2(eri,rrr,rxyz,qx,qz,6,3)
      return
end


!-----------------------------------------------
  subroutine int2dddp2(eri,rrr,rxyz,qx,qz,i,j)
!-----------------------------------------------
!
      implicit none
      integer,intent(in) :: i, j
      real(8),parameter :: two=2.0D+00, four=4.0D+00
      real(8),intent(in) ::rrr(15), rxyz(45), qx, qz
      real(8),intent(inout) :: eri(6,6,6,3)
      real(8) :: xx, xz, zz, xxx, xxz, xzz, zzz, xxxx, xxxz, xxzz, xzzz, zzzz
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
      eri(1,1,i,j)=rrr( 1)+(rxyz(37)+rxyz(38))*two*qx+(rxyz(23)+rxyz(24)*four+rxyz(25))*xx &
&                  +(rxyz(26)+rxyz(27))*two*xxx+rxyz(1)*xxxx
      eri(2,1,i,j)=rrr( 4)+rxyz(14)*two*qx+rxyz(5)*xx
      eri(3,1,i,j)=rrr( 6)+rxyz(29)*two*qx+rxyz(30)*two*qz+rxyz(34)*xx+rxyz(21)*two*xz &
&                  +rxyz(23)*zz+rxyz(36)*two*xxz+rxyz(26)*two*xzz+rxyz(1)*xxzz
      eri(4,1,i,j)=rrr( 2)+(rxyz(16)*two+rxyz(17))*qx+(rxyz(42)*two+rxyz(43))*xx+rxyz(2)*xxx
      eri(5,1,i,j)=rrr( 3)+(rxyz(30)+rxyz(31)*two)*qx+rxyz(37)*qz+(rxyz(21)+rxyz(22))*xx &
&                  +(rxyz(23)+rxyz(24)*two)*xz+rxyz(36)*xxx+(rxyz(26)*two+rxyz(27))*xxz &
&                  +rxyz(1)*xxxz
      eri(6,1,i,j)=rrr( 5)+rxyz(40)*two*qx+rxyz(17)*qz+rxyz(11)*xx+rxyz(42)*two*xz+rxyz(2)*xxz
      eri(1,2,i,j)=rrr( 4)+rxyz(15)*two*qx+rxyz(4)*xx
      eri(2,2,i,j)=rrr(11)
      eri(3,2,i,j)=rrr(13)+rxyz(12)*two*qz+rxyz(4)*zz
      eri(4,2,i,j)=rrr( 7)+rxyz(7)*qx
      eri(5,2,i,j)=rrr( 8)+rxyz(12)*qx+rxyz(15)*qz+rxyz(4)*xz
      eri(6,2,i,j)=rrr(12)+rxyz(7)*qz
      eri(1,3,i,j)=rrr( 6)+rxyz(28)*two*qx+rxyz(31)*two*qz+rxyz(32)*xx+rxyz(21)*two*xz &
&                  +rxyz(25)*zz+rxyz(35)*two*xxz+rxyz(27)*two*xzz+rxyz(1)*xxzz
      eri(2,3,i,j)=rrr(13)+rxyz(13)*two*qz+rxyz(5)*zz
      eri(3,3,i,j)=rrr(15)+(rxyz(44)+rxyz(45))*two*qz+(rxyz(32)+rxyz(33)*four+rxyz(34))*zz &
&                  +(rxyz(35)+rxyz(36))*two*zzz+rxyz(1)*zzzz
      eri(4,3,i,j)=rrr( 9)+rxyz(18)*qx+rxyz(40)*two*qz+rxyz(10)*two*xz+rxyz(43)*zz+rxyz(2)*xzz
      eri(5,3,i,j)=rrr(10)+rxyz(44)*qx+(rxyz(28)+rxyz(29)*two)*qz+(rxyz(32)+rxyz(33)*two)*xz &
&                  +(rxyz(21)+rxyz(22))*zz+(rxyz(35)*two+rxyz(36))*xzz+rxyz(27)*zzz+rxyz(1)*xzzz
      eri(6,3,i,j)=rrr(14)+(rxyz(18)+rxyz(19)*two)*qz+(rxyz(10)*two+rxyz(11))*zz+rxyz(2)*zzz
      eri(1,4,i,j)=rrr( 2)+(rxyz(16)+rxyz(17)*two)*qx+(rxyz(41)+rxyz(42)*two)*xx+rxyz(3)*xxx
      eri(2,4,i,j)=rrr( 7)+rxyz(8)*qx
      eri(3,4,i,j)=rrr( 9)+rxyz(19)*qx+rxyz(39)*two*qz+rxyz(10)*two*xz+rxyz(41)*zz+rxyz(3)*xzz
      eri(4,4,i,j)=rrr( 4)+(rxyz(14)+rxyz(15))*qx+rxyz(6)*xx
      eri(5,4,i,j)=rrr( 5)+(rxyz(39)+rxyz(40))*qx+rxyz(17)*qz+rxyz(10)*xx+(rxyz(41)+rxyz(42))*xz &
&                  +rxyz(3)*xxz
      eri(6,4,i,j)=rrr( 8)+rxyz(13)*qx+rxyz(15)*qz+rxyz(6)*xz
      eri(1,5,i,j)=rrr( 3)+(rxyz(30)*two+rxyz(31))*qx+rxyz(38)*qz+(rxyz(20)+rxyz(21))*xx &
&                  +(rxyz(24)*two+rxyz(25))*xz+rxyz(35)*xxx+(rxyz(26)+rxyz(27)*two)*xxz &
&                  +rxyz(1)*xxxz
      eri(2,5,i,j)=rrr( 8)+rxyz(13)*qx+rxyz(14)*qz+rxyz(5)*xz
      eri(3,5,i,j)=rrr(10)+rxyz(45)*qx+(rxyz(28)*two+rxyz(29))*qz+(rxyz(33)*two+rxyz(34))*xz &
&                  +(rxyz(20)+rxyz(21))*zz+(rxyz(35)+rxyz(36)*two)*xzz+rxyz(26)*zzz+rxyz(1)*xzzz
      eri(4,5,i,j)=rrr( 5)+(rxyz(39)+rxyz(40))*qx+rxyz(16)*qz+rxyz(10)*xx+(rxyz(42)+rxyz(43))*xz &
&                  +rxyz(2)*xxz
      eri(5,5,i,j)=rrr( 6)+(rxyz(28)+rxyz(29))*qx+(rxyz(30)+rxyz(31))*qz+rxyz(33)*xx &
&                  +(rxyz(20)+rxyz(21)+rxyz(22))*xz+rxyz(24)*zz+(rxyz(35)+rxyz(36))*xxz &
&                  +(rxyz(26)+rxyz(27))*xzz+rxyz(1)*xxzz
      eri(6,5,i,j)=rrr( 9)+rxyz(19)*qx+(rxyz(39)+rxyz(40))*qz+(rxyz(10)+rxyz(11))*xz+rxyz(42)*zz &
&                  +rxyz(2)*xzz
      eri(1,6,i,j)=rrr( 5)+rxyz(39)*two*qx+rxyz(16)*qz+rxyz(9)*xx+rxyz(42)*two*xz+rxyz(3)*xxz
      eri(2,6,i,j)=rrr(12)+rxyz(8)*qz
      eri(3,6,i,j)=rrr(14)+(rxyz(18)*two+rxyz(19))*qz+(rxyz(9)+rxyz(10)*two)*zz+rxyz(3)*zzz
      eri(4,6,i,j)=rrr( 8)+rxyz(12)*qx+rxyz(14)*qz+rxyz(6)*xz
      eri(5,6,i,j)=rrr( 9)+rxyz(18)*qx+(rxyz(39)+rxyz(40))*qz+(rxyz(9)+rxyz(10))*xz+rxyz(42)*zz &
&                  +rxyz(3)*xzz
      eri(6,6,i,j)=rrr(13)+(rxyz(12)+rxyz(13))*qz+rxyz(6)*zz
!
      return
end


