!--------------------------------------------------------------------------
  subroutine int2ddps(phmdint,exfac1,exfac2,xyziq,xzkl,rot,nijkl,nbfijkl)
!--------------------------------------------------------------------------
!
! Calculate (dd|ps) integrals
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
      real(8),parameter :: six=6.0D+00, seven=7.0D+00, nine=9.0D+00, ten=1.0D+01, p15=1.5D+01
      real(8),parameter :: sqrt3=1.73205080756888D+00, sqrt3h=8.660254037844386D-01
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*), xzkl(2), rot(3,3)
      real(8),intent(out) :: phmdint(6,6,6,6)
      real(8) :: tval, tval1, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:5)
      real(8) :: f0, f1(2,2), f2(3,2), f3(4,2), f4(5,2), f5(6), ftw(5,2)
      real(8) :: r0(5), r1(3,9), r2(6,8), r3(10,6), r4(15,3), r5(21)
      real(8) :: ex12, ex34, ex43, ex41, expq, expq2, expq4, ex3q, ex4q, c12, c34, zip
      real(8) :: xiq, yiq, ziq, xiq2, yiq2, xyiq, xiq4, yiq4, xyiq2, xypq2, zpq, zpq2, fac
      real(8) :: ex33q, ex34q, ex44q, zjp, pmd, qmd, qmd2, qmd3, qmd4, qmd4x, qmd4y, qmd4xy
      real(8) :: qx, qz, xx, xz, zz, xxx, xxz, xzz, zzz, xxxx, xxxz, xxzz, xzzz, zzzz
      real(8) :: eri(6,6,3), work(15), f1w(3,2), f2w(6,2), f3w(10,2), f4w(15,2), rot2(6,6)
      real(8) :: rot3(6,5)
!
! Zero-clear
!
      r0(1:5)     = zero
      r1(1:3,1:9) = zero
      r2(1:6,1:8) = zero
      r3(1:10,1:6)= zero
      r4(1:15,1:3)= zero
      r5(1:21)    = zero
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
        xypq2= xiq2+yiq2
        f0         = zero
        f1(1:2,1:2)= zero
        f2(1:3,1:2)= zero
        f3(1:4,1:2)= zero
        f4(1:5,1:2)= zero
        f5(1:6)    = zero
        do ij= 1,nijkl(1)
          ex12= exfac1(1,ij)
          pmd = exfac1(2,ij)
          zjp = exfac1(3,ij)
          zip = exfac1(4,ij)
          c12 = exfac1(5,ij)
          ex41= one/(ex12+ex34)
          zpq = ziq-zip
          zpq2= zpq*zpq
          expq= ex12*ex34*ex41
          tval=(xypq2+zpq2)*expq
!
! Calculate Fm(T)
!
          if(tval >= threshtval) then
            tinv= one/tval
            ft(0)= c12*sqrtpi4*sqrt(ex41*tval)*tinv
            ft(1)= ft(0)*tinv*expq
            ft(2)= ft(1)*tinv*expq*three
            ft(3)= ft(2)*tinv*expq*five
            ft(4)= ft(3)*tinv*expq*seven
            ft(5)= ft(4)*tinv*expq*nine
          else
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
            do ii= 0,5
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
            fac= c12*sqrt(ex41)
            expq= expq*two
            expq2= expq*expq
            expq4= expq2*expq2
            ft(0)= ft(0)*fac
            ft(1)= ft(1)*fac*expq
            ft(2)= ft(2)*fac*expq2
            ft(3)= ft(3)*fac*expq2*expq
            ft(4)= ft(4)*fac*expq4
            ft(5)= ft(5)*fac*expq4*expq
          endif
          ftw(1,1)= zjp
          ftw(1,2)= pmd
          do i= 1,2
            ftw(2,i)= ftw(1,i)*zpq
            ftw(3,i)= ftw(1,i)*zpq2
            ftw(4,i)= ftw(1,i)*zpq2*zpq
            ftw(5,i)= ftw(1,i)*zpq2*zpq2
          enddo
          f0= f0+ft(0)*ftw(1,1)
          do i= 1,2
            f1(1,i)= f1(1,i)-ft(1)*ftw(1,i)
            f1(2,i)= f1(2,i)-ft(1)*ftw(2,i)
          enddo
          do i= 1,2
            f2(1,i)= f2(1,i)+ft(2)*ftw(1,i)
            f2(2,i)= f2(2,i)+ft(2)*ftw(2,i)
            f2(3,i)= f2(3,i)+ft(2)*ftw(3,i)
          enddo
          do i= 1,2
            do j= 1,4
              f3(j,i)= f3(j,i)-ft(3)*ftw(j,i)
            enddo
          enddo
          do i= 1,2
            do j= 1,5
              f4(j,i)= f4(j,i)+ft(4)*ftw(j,i)
            enddo
          enddo
          f5(1)= f5(1)-ft(5)*ftw(1,2)
          f5(2)= f5(2)-ft(5)*ftw(2,2)
          f5(3)= f5(3)-ft(5)*ftw(3,2)
          f5(4)= f5(4)-ft(5)*ftw(4,2)
          f5(5)= f5(5)-ft(5)*ftw(5,2)
          f5(6)= f5(6)-ft(5)*ftw(5,2)*zpq
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
          r0(i)= r0(i)+f0*work(i)
        enddo
!
        do i= 1,2
          f1w(1,i)= f1(1,i)*xiq
          f1w(2,i)= f1(1,i)*yiq
          f1w(3,i)= f1(2,i)
        enddo
        do i= 1,4
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,1)*work(i+5)
          enddo
        enddo
        do i= 5,9
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,2)*work(i-4)
          enddo
        enddo
!
        do i= 1,2
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
            r2(j,i)= r2(j,i)+f2w(j,2)*work(i+1)
          enddo
        enddo
! 
        do i= 1,2
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
        do i= 3,6
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,2)*work(i+7)
          enddo
        enddo
!
        do i= 1,2
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
          r4(j,2)= r4(j,2)+f4w(j,2)*work(14)
          r4(j,3)= r4(j,3)+f4w(j,2)*work(15)
        enddo
!
        r5( 1)= r5( 1)+(f5(1)*xiq4 +f4(1,2)*xiq2*ten               +f3(1,2)*p15  )*qmd4x
        r5( 2)= r5( 2)+(f5(1)*xiq4 +f4(1,2)*xiq2*six               +f3(1,2)*three)*qmd4y
        r5( 3)= r5( 3)+(f5(2)*xiq4 +f4(2,2)*xiq2*six               +f3(2,2)*three)*qmd4
        r5( 4)= r5( 4)+(f5(1)*xyiq2+f4(1,2)*xiq2+f4(1,2)*yiq2*three+f3(1,2)*three)*qmd4x
        r5( 5)= r5( 5)+(f5(2)*xiq2 +f4(2,2)*three                                )*qmd4xy
        r5( 6)= r5( 6)+(f5(3)*xiq2 +f4(1,2)*xiq2+f4(3,2)*three     +f3(1,2)*three)*qmd4x
        r5( 7)= r5( 7)+(f5(1)*xyiq2+f4(1,2)*xiq2*three+f4(1,2)*yiq2+f3(1,2)*three)*qmd4y
        r5( 8)= r5( 8)+(f5(2)*xyiq2+f4(2,2)*xiq2+f4(2,2)*yiq2      +f3(2,2)      )*qmd4
        r5( 9)= r5( 9)+(f5(3)*xiq2 +f4(1,2)*xiq2+f4(3,2)           +f3(1,2)      )*qmd4y
        r5(10)= r5(10)+(f5(4)*xiq2 +f4(2,2)*xiq2*three+f4(4,2)     +f3(2,2)*three)*qmd4
        r5(11)= r5(11)+(f5(1)*yiq4 +f4(1,2)*yiq2*six               +f3(1,2)*three)*qmd4x
        r5(12)= r5(12)+(f5(2)*yiq2 +f4(2,2)*three                                )*qmd4xy
        r5(13)= r5(13)+(f5(3)*yiq2 +f4(1,2)*yiq2+f4(3,2)           +f3(1,2)      )*qmd4x
        r5(14)= r5(14)+(f5(4)      +f4(2,2)*three                                )*qmd4xy
        r5(15)= r5(15)+(f5(5)      +f4(3,2)*six                    +f3(1,2)*three)*qmd4x
        r5(16)= r5(16)+(f5(1)*yiq4 +f4(1,2)*yiq2*ten               +f3(1,2)*p15  )*qmd4y
        r5(17)= r5(17)+(f5(2)*yiq4 +f4(2,2)*yiq2*six               +f3(2,2)*three)*qmd4
        r5(18)= r5(18)+(f5(3)*yiq2 +f4(1,2)*yiq2+f4(3,2)*three     +f3(1,2)*three)*qmd4y
        r5(19)= r5(19)+(f5(4)*yiq2 +f4(2,2)*yiq2*three+f4(4,2)     +f3(2,2)*three)*qmd4
        r5(20)= r5(20)+(f5(5)      +f4(3,2)*six                    +f3(1,2)*three)*qmd4y
        r5(21)= r5(21)+(f5(6)      +f4(4,2)*ten                    +f3(2,2)*p15  )*qmd4
      enddo
!
      qx= xzkl(1)
      qz= xzkl(2)
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
      eri(1,1,1)=-r5(1)-r3(1,3)*six-r1(1,5)*three+(-r4(1,2)*two-r4(1,3)*two-r2(1,5)*six &
&                -r2(1,6)*six)*qx+(-r3(1,4)-r3(1,5)*four-r3(1,6)-r1(1,6)-r1(1,7)*four-r1(1,8) &
&                )*xx+(-r2(1,7)*two-r2(1,8)*two)*xxx-r1(1,9)*xxxx
      eri(2,1,1)=-r5(4)-r3(1,3)-r3(4,3)-r1(1,5)+(-r4(4,3)*two-r2(1,6)*two)*qx+(-r3(4,6) &
&                -r1(1,8))*xx
      eri(3,1,1)=-r5(6)-r3(1,3)-r3(6,3)-r1(1,5)+(-r4(6,3)*two-r2(1,6)*two)*qx+( &
&                -r4(3,2)*two-r2(5,5)*two)*qz+(-r3(6,6)-r1(1,8))*xx-r3(3,5)*four*xz+(-r3(1,4) &
&                -r1(1,6))*zz-r2(5,8)*two*xxz-r2(1,7)*two*xzz-r1(1,9)*xxzz
      eri(4,1,1)=-r5(2)-r3(2,3)*three+(-r4(2,2)-r4(2,3)*two-r2(4,5)-r2(4,6)*two)*qx+( &
&                -r3(2,5)*two-r3(2,6))*xx-r2(4,8)*xxx
      eri(5,1,1)=-r5(3)-r3(3,3)*three+(-r4(3,2)-r4(3,3)*two-r2(5,5)-r2(5,6)*two)*qx+( &
&                -r4(1,2)-r2(1,5)*three)*qz+(-r3(3,5)*two-r3(3,6))*xx+(-r3(1,4)-r3(1,5)*two &
&                -r1(1,6)-r1(1,7)*two)*xz-r2(5,8)*xxx+(-r2(1,7)*two-r2(1,8))*xxz-r1(1,9)*xxxz
      eri(6,1,1)=-r5(5)-r3(5,3)-r4(5,3)*two*qx+(-r4(2,2)-r2(4,5))*qz-r3(5,6)*xx &
&                -r3(2,5)*two*xz-r2(4,8)*xxz
      eri(1,2,1)=-r5(4)-r3(1,3)-r3(4,3)-r1(1,5)+(-r4(4,2)*two-r2(1,5)*two)*qx+(-r3(4,4) &
&                -r1(1,6))*xx
      eri(2,2,1)=-r5(11)-r3(4,3)*six-r1(1,5)*three
      eri(3,2,1)=-r5(13)-r3(4,3)-r3(6,3)-r1(1,5)+(-r4(8,2)*two-r2(5,5)*two)*qz+(-r3(4,4) &
&                -r1(1,6))*zz
      eri(4,2,1)=-r5(7)-r3(2,3)*three+(-r4(7,2)-r2(4,5)*three)*qx
      eri(5,2,1)=-r5(8)-r3(3,3)+(-r4(8,2)-r2(5,5))*qx+(-r4(4,2)-r2(1,5))*qz+(-r3(4,4) &
&                -r1(1,6))*xz
      eri(6,2,1)=-r5(12)-r3(5,3)*three+(-r4(7,2)-r2(4,5)*three)*qz
      eri(1,3,1)=-r5(6)-r3(1,3)-r3(6,3)-r1(1,5)+(-r4(6,2)*two-r2(1,5)*two)*qx+( &
&                -r4(3,3)*two-r2(5,6)*two)*qz+(-r3(6,4)-r1(1,6))*xx-r3(3,5)*four*xz+(-r3(1,6) &
&                -r1(1,8))*zz-r2(5,7)*two*xxz-r2(1,8)*two*xzz-r1(1,9)*xxzz
      eri(2,3,1)=-r5(13)-r3(4,3)-r3(6,3)-r1(1,5)+(-r4(8,3)*two-r2(5,6)*two)*qz+(-r3(4,6) &
&                -r1(1,8))*zz
      eri(3,3,1)=-r5(15)-r3(6,3)*six-r1(1,5)*three+(-r4(10,2)*two-r4(10,3)*two-r2(5,5)*six &
&                -r2(5,6)*six)*qz+(-r3(6,4)-r3(6,5)*four-r3(6,6)-r1(1,6)-r1(1,7)*four-r1(1,8) &
&                )*zz+(-r2(5,7)*two-r2(5,8)*two)*zzz-r1(1,9)*zzzz
      eri(4,3,1)=-r5(9)-r3(2,3)+(-r4(9,2)-r2(4,5))*qx-r4(5,3)*two*qz-r3(5,5)*two*xz &
&                -r3(2,6)*zz-r2(4,8)*xzz
      eri(5,3,1)=-r5(10)-r3(3,3)*three+(-r4(10,2)-r2(5,5)*three)*qx+(-r4(6,2)-r4(6,3)*two &
&                -r2(1,5)-r2(1,6)*two)*qz+(-r3(6,4)-r3(6,5)*two-r1(1,6)-r1(1,7)*two)*xz+( &
&                -r3(3,5)*two-r3(3,6))*zz+(-r2(5,7)*two-r2(5,8))*xzz-r2(1,8)*zzz-r1(1,9)*xzzz
      eri(6,3,1)=-r5(14)-r3(5,3)*three+(-r4(9,2)-r4(9,3)*two-r2(4,5)-r2(4,6)*two)*qz+( &
&                -r3(5,5)*two-r3(5,6))*zz-r2(4,8)*zzz
      eri(1,4,1)=-r5(2)-r3(2,3)*three+(-r4(2,2)*two-r4(2,3)-r2(4,5)*two-r2(4,6))*qx+( &
&                -r3(2,4)-r3(2,5)*two)*xx-r2(4,7)*xxx
      eri(2,4,1)=-r5(7)-r3(2,3)*three+(-r4(7,3)-r2(4,6)*three)*qx
      eri(3,4,1)=-r5(9)-r3(2,3)+(-r4(9,3)-r2(4,6))*qx-r4(5,2)*two*qz-r3(5,5)*two*xz &
&                -r3(2,4)*zz-r2(4,7)*xzz
      eri(4,4,1)=-r5(4)-r3(1,3)-r3(4,3)-r1(1,5)+(-r4(4,2)-r4(4,3)-r2(1,5)-r2(1,6))*qx+( &
&                -r3(4,5)-r1(1,7))*xx
      eri(5,4,1)=-r5(5)-r3(5,3)+(-r4(5,2)-r4(5,3))*qx+(-r4(2,2)-r2(4,5))*qz-r3(5,5)*xx+( &
&                -r3(2,4)-r3(2,5))*xz-r2(4,7)*xxz
      eri(6,4,1)=-r5(8)-r3(3,3)+(-r4(8,3)-r2(5,6))*qx+(-r4(4,2)-r2(1,5))*qz+(-r3(4,5) &
&                -r1(1,7))*xz
      eri(1,5,1)=-r5(3)-r3(3,3)*three+(-r4(3,2)*two-r4(3,3)-r2(5,5)*two-r2(5,6))*qx+( &
&                -r4(1,3)-r2(1,6)*three)*qz+(-r3(3,4)-r3(3,5)*two)*xx+(-r3(1,5)*two-r3(1,6) &
&                -r1(1,7)*two-r1(1,8))*xz-r2(5,7)*xxx+(-r2(1,7)-r2(1,8)*two)*xxz-r1(1,9)*xxxz
      eri(2,5,1)=-r5(8)-r3(3,3)+(-r4(8,3)-r2(5,6))*qx+(-r4(4,3)-r2(1,6))*qz+(-r3(4,6) &
&                -r1(1,8))*xz
      eri(3,5,1)=-r5(10)-r3(3,3)*three+(-r4(10,3)-r2(5,6)*three)*qx+(-r4(6,2)*two-r4(6,3) &
&                -r2(1,5)*two-r2(1,6))*qz+(-r3(6,5)*two-r3(6,6)-r1(1,7)*two-r1(1,8))*xz+( &
&                -r3(3,4)-r3(3,5)*two)*zz+(-r2(5,7)-r2(5,8)*two)*xzz-r2(1,7)*zzz-r1(1,9)*xzzz
      eri(4,5,1)=-r5(5)-r3(5,3)+(-r4(5,2)-r4(5,3))*qx+(-r4(2,3)-r2(4,6))*qz-r3(5,5)*xx+( &
&                -r3(2,5)-r3(2,6))*xz-r2(4,8)*xxz
      eri(5,5,1)=-r5(6)-r3(1,3)-r3(6,3)-r1(1,5)+(-r4(6,2)-r4(6,3)-r2(1,5)-r2(1,6))*qx+( &
&                -r4(3,2)-r4(3,3)-r2(5,5)-r2(5,6))*qz+(-r3(6,5)-r1(1,7))*xx+(-r3(3,4) &
&                -r3(3,5)*two-r3(3,6))*xz+(-r3(1,5)-r1(1,7))*zz+(-r2(5,7)-r2(5,8))*xxz+( &
&                -r2(1,7)-r2(1,8))*xzz-r1(1,9)*xxzz
      eri(6,5,1)=-r5(9)-r3(2,3)+(-r4(9,3)-r2(4,6))*qx+(-r4(5,2)-r4(5,3))*qz+(-r3(5,5) &
&                -r3(5,6))*xz-r3(2,5)*zz-r2(4,8)*xzz
      eri(1,6,1)=-r5(5)-r3(5,3)-r4(5,2)*two*qx+(-r4(2,3)-r2(4,6))*qz-r3(5,4)*xx &
&                -r3(2,5)*two*xz-r2(4,7)*xxz
      eri(2,6,1)=-r5(12)-r3(5,3)*three+(-r4(7,3)-r2(4,6)*three)*qz
      eri(3,6,1)=-r5(14)-r3(5,3)*three+(-r4(9,2)*two-r4(9,3)-r2(4,5)*two-r2(4,6))*qz+( &
&                -r3(5,4)-r3(5,5)*two)*zz-r2(4,7)*zzz
      eri(4,6,1)=-r5(8)-r3(3,3)+(-r4(8,2)-r2(5,5))*qx+(-r4(4,3)-r2(1,6))*qz+(-r3(4,5) &
&                -r1(1,7))*xz
      eri(5,6,1)=-r5(9)-r3(2,3)+(-r4(9,2)-r2(4,5))*qx+(-r4(5,2)-r4(5,3))*qz+(-r3(5,4) &
&                -r3(5,5))*xz-r3(2,5)*zz-r2(4,7)*xzz
      eri(6,6,1)=-r5(13)-r3(4,3)-r3(6,3)-r1(1,5)+(-r4(8,2)-r4(8,3)-r2(5,5)-r2(5,6))*qz+( &
&                -r3(4,5)-r1(1,7))*zz
      eri(1,1,2)=-r5(2)-r3(2,3)*six-r1(2,5)*three+(-r4(2,2)*two-r4(2,3)*two-r2(4,5)*six &
&                -r2(4,6)*six)*qx+(-r3(2,4)-r3(2,5)*four-r3(2,6)-r1(2,6)-r1(2,7)*four-r1(2,8) &
&                )*xx+(-r2(4,7)*two-r2(4,8)*two)*xxx-r1(2,9)*xxxx
      eri(2,1,2)=-r5(7)-r3(2,3)-r3(7,3)-r1(2,5)+(-r4(7,3)*two-r2(4,6)*two)*qx+(-r3(7,6) &
&                -r1(2,8))*xx
      eri(3,1,2)=-r5(9)-r3(2,3)-r3(9,3)-r1(2,5)+(-r4(9,3)*two-r2(4,6)*two)*qx+( &
&                -r4(5,2)*two-r2(6,5)*two)*qz+(-r3(9,6)-r1(2,8))*xx-r3(5,5)*four*xz+(-r3(2,4) &
&                -r1(2,6))*zz-r2(6,8)*two*xxz-r2(4,7)*two*xzz-r1(2,9)*xxzz
      eri(4,1,2)=-r5(4)-r3(4,3)*three+(-r4(4,2)-r4(4,3)*two-r2(2,5)-r2(2,6)*two)*qx+( &
&                -r3(4,5)*two-r3(4,6))*xx-r2(2,8)*xxx
      eri(5,1,2)=-r5(5)-r3(5,3)*three+(-r4(5,2)-r4(5,3)*two-r2(6,5)-r2(6,6)*two)*qx+( &
&                -r4(2,2)-r2(4,5)*three)*qz+(-r3(5,5)*two-r3(5,6))*xx+(-r3(2,4)-r3(2,5)*two &
&                -r1(2,6)-r1(2,7)*two)*xz-r2(6,8)*xxx+(-r2(4,7)*two-r2(4,8))*xxz-r1(2,9)*xxxz
      eri(6,1,2)=-r5(8)-r3(8,3)-r4(8,3)*two*qx+(-r4(4,2)-r2(2,5))*qz-r3(8,6)*xx &
&                -r3(4,5)*two*xz-r2(2,8)*xxz
      eri(1,2,2)=-r5(7)-r3(2,3)-r3(7,3)-r1(2,5)+(-r4(7,2)*two-r2(4,5)*two)*qx+(-r3(7,4) &
&                -r1(2,6))*xx
      eri(2,2,2)=-r5(16)-r3(7,3)*six-r1(2,5)*three
      eri(3,2,2)=-r5(18)-r3(7,3)-r3(9,3)-r1(2,5)+(-r4(12,2)*two-r2(6,5)*two)*qz+(-r3(7,4) &
&                -r1(2,6))*zz
      eri(4,2,2)=-r5(11)-r3(4,3)*three+(-r4(11,2)-r2(2,5)*three)*qx
      eri(5,2,2)=-r5(12)-r3(5,3)+(-r4(12,2)-r2(6,5))*qx+(-r4(7,2)-r2(4,5))*qz+(-r3(7,4) &
&                -r1(2,6))*xz
      eri(6,2,2)=-r5(17)-r3(8,3)*three+(-r4(11,2)-r2(2,5)*three)*qz
      eri(1,3,2)=-r5(9)-r3(2,3)-r3(9,3)-r1(2,5)+(-r4(9,2)*two-r2(4,5)*two)*qx+( &
&                -r4(5,3)*two-r2(6,6)*two)*qz+(-r3(9,4)-r1(2,6))*xx-r3(5,5)*four*xz+(-r3(2,6) &
&                -r1(2,8))*zz-r2(6,7)*two*xxz-r2(4,8)*two*xzz-r1(2,9)*xxzz
      eri(2,3,2)=-r5(18)-r3(7,3)-r3(9,3)-r1(2,5)+(-r4(12,3)*two-r2(6,6)*two)*qz+(-r3(7,6) &
&                -r1(2,8))*zz
      eri(3,3,2)=-r5(20)-r3(9,3)*six-r1(2,5)*three+(-r4(14,2)*two-r4(14,3)*two-r2(6,5)*six &
&                -r2(6,6)*six)*qz+(-r3(9,4)-r3(9,5)*four-r3(9,6)-r1(2,6)-r1(2,7)*four-r1(2,8) &
&                )*zz+(-r2(6,7)*two-r2(6,8)*two)*zzz-r1(2,9)*zzzz
      eri(4,3,2)=-r5(13)-r3(4,3)+(-r4(13,2)-r2(2,5))*qx-r4(8,3)*two*qz-r3(8,5)*two*xz &
&                -r3(4,6)*zz-r2(2,8)*xzz
      eri(5,3,2)=-r5(14)-r3(5,3)*three+(-r4(14,2)-r2(6,5)*three)*qx+(-r4(9,2)-r4(9,3)*two &
&                -r2(4,5)-r2(4,6)*two)*qz+(-r3(9,4)-r3(9,5)*two-r1(2,6)-r1(2,7)*two)*xz+( &
&                -r3(5,5)*two-r3(5,6))*zz+(-r2(6,7)*two-r2(6,8))*xzz-r2(4,8)*zzz-r1(2,9)*xzzz
      eri(6,3,2)=-r5(19)-r3(8,3)*three+(-r4(13,2)-r4(13,3)*two-r2(2,5)-r2(2,6)*two)*qz+( &
&                -r3(8,5)*two-r3(8,6))*zz-r2(2,8)*zzz
      eri(1,4,2)=-r5(4)-r3(4,3)*three+(-r4(4,2)*two-r4(4,3)-r2(2,5)*two-r2(2,6))*qx+( &
&                -r3(4,4)-r3(4,5)*two)*xx-r2(2,7)*xxx
      eri(2,4,2)=-r5(11)-r3(4,3)*three+(-r4(11,3)-r2(2,6)*three)*qx
      eri(3,4,2)=-r5(13)-r3(4,3)+(-r4(13,3)-r2(2,6))*qx-r4(8,2)*two*qz-r3(8,5)*two*xz &
&                -r3(4,4)*zz-r2(2,7)*xzz
      eri(4,4,2)=-r5(7)-r3(2,3)-r3(7,3)-r1(2,5)+(-r4(7,2)-r4(7,3)-r2(4,5)-r2(4,6))*qx+( &
&                -r3(7,5)-r1(2,7))*xx
      eri(5,4,2)=-r5(8)-r3(8,3)+(-r4(8,2)-r4(8,3))*qx+(-r4(4,2)-r2(2,5))*qz-r3(8,5)*xx+( &
&                -r3(4,4)-r3(4,5))*xz-r2(2,7)*xxz
      eri(6,4,2)=-r5(12)-r3(5,3)+(-r4(12,3)-r2(6,6))*qx+(-r4(7,2)-r2(4,5))*qz+(-r3(7,5) &
&                -r1(2,7))*xz
      eri(1,5,2)=-r5(5)-r3(5,3)*three+(-r4(5,2)*two-r4(5,3)-r2(6,5)*two-r2(6,6))*qx+( &
&                -r4(2,3)-r2(4,6)*three)*qz+(-r3(5,4)-r3(5,5)*two)*xx+(-r3(2,5)*two-r3(2,6) &
&                -r1(2,7)*two-r1(2,8))*xz-r2(6,7)*xxx+(-r2(4,7)-r2(4,8)*two)*xxz-r1(2,9)*xxxz
      eri(2,5,2)=-r5(12)-r3(5,3)+(-r4(12,3)-r2(6,6))*qx+(-r4(7,3)-r2(4,6))*qz+(-r3(7,6) &
&                -r1(2,8))*xz
      eri(3,5,2)=-r5(14)-r3(5,3)*three+(-r4(14,3)-r2(6,6)*three)*qx+(-r4(9,2)*two-r4(9,3) &
&                -r2(4,5)*two-r2(4,6))*qz+(-r3(9,5)*two-r3(9,6)-r1(2,7)*two-r1(2,8))*xz+( &
&                -r3(5,4)-r3(5,5)*two)*zz+(-r2(6,7)-r2(6,8)*two)*xzz-r2(4,7)*zzz-r1(2,9)*xzzz
      eri(4,5,2)=-r5(8)-r3(8,3)+(-r4(8,2)-r4(8,3))*qx+(-r4(4,3)-r2(2,6))*qz-r3(8,5)*xx+( &
&                -r3(4,5)-r3(4,6))*xz-r2(2,8)*xxz
      eri(5,5,2)=-r5(9)-r3(2,3)-r3(9,3)-r1(2,5)+(-r4(9,2)-r4(9,3)-r2(4,5)-r2(4,6))*qx+( &
&                -r4(5,2)-r4(5,3)-r2(6,5)-r2(6,6))*qz+(-r3(9,5)-r1(2,7))*xx+(-r3(5,4) &
&                -r3(5,5)*two-r3(5,6))*xz+(-r3(2,5)-r1(2,7))*zz+(-r2(6,7)-r2(6,8))*xxz+( &
&                -r2(4,7)-r2(4,8))*xzz-r1(2,9)*xxzz
      eri(6,5,2)=-r5(13)-r3(4,3)+(-r4(13,3)-r2(2,6))*qx+(-r4(8,2)-r4(8,3))*qz+(-r3(8,5) &
&                -r3(8,6))*xz-r3(4,5)*zz-r2(2,8)*xzz
      eri(1,6,2)=-r5(8)-r3(8,3)-r4(8,2)*two*qx+(-r4(4,3)-r2(2,6))*qz-r3(8,4)*xx &
&                -r3(4,5)*two*xz-r2(2,7)*xxz
      eri(2,6,2)=-r5(17)-r3(8,3)*three+(-r4(11,3)-r2(2,6)*three)*qz
      eri(3,6,2)=-r5(19)-r3(8,3)*three+(-r4(13,2)*two-r4(13,3)-r2(2,5)*two-r2(2,6))*qz+( &
&                -r3(8,4)-r3(8,5)*two)*zz-r2(2,7)*zzz
      eri(4,6,2)=-r5(12)-r3(5,3)+(-r4(12,2)-r2(6,5))*qx+(-r4(7,3)-r2(4,6))*qz+(-r3(7,5) &
&                -r1(2,7))*xz
      eri(5,6,2)=-r5(13)-r3(4,3)+(-r4(13,2)-r2(2,5))*qx+(-r4(8,2)-r4(8,3))*qz+(-r3(8,4) &
&                -r3(8,5))*xz-r3(4,5)*zz-r2(2,7)*xzz
      eri(6,6,2)=-r5(18)-r3(7,3)-r3(9,3)-r1(2,5)+(-r4(12,2)-r4(12,3)-r2(6,5)-r2(6,6))*qz+( &
&                -r3(7,5)-r1(2,7))*zz
      eri(1,1,3)=-r5(3)+r4(1,1)-r3(3,3)*six+r2(1,1)*six-r1(3,5)*three+r0(1)*three+( &
&                -r4(3,2)*two-r4(3,3)*two+r3(1,1)*two+r3(1,2)*two-r2(5,5)*six-r2(5,6)*six &
&                +r1(1,1)*six+r1(1,2)*six)*qx+(-r3(3,4)-r3(3,5)*four-r3(3,6)+r2(1,2) &
&                +r2(1,3)*four+r2(1,4)-r1(3,6)-r1(3,7)*four-r1(3,8)+r0(2)+r0(3)*four+r0(4)) &
&                *xx+(-r2(5,7)*two-r2(5,8)*two+r1(1,3)*two+r1(1,4)*two)*xxx+(-r1(3,9)+r0(5)) &
&                *xxxx
      eri(2,1,3)=-r5(8)+r4(4,1)-r3(3,3)-r3(8,3)+r2(1,1)+r2(2,1)-r1(3,5)+r0(1)+( &
&                -r4(8,3)*two+r3(4,2)*two-r2(5,6)*two+r1(1,2)*two)*qx+(-r3(8,6)+r2(2,4) &
&                -r1(3,8)+r0(4))*xx
      eri(3,1,3)=-r5(10)+r4(6,1)-r3(3,3)-r3(10,3)+r2(1,1)+r2(3,1)-r1(3,5)+r0(1)+( &
&                -r4(10,3)*two+r3(6,2)*two-r2(5,6)*two+r1(1,2)*two)*qx+(-r4(6,2)*two &
&                +r3(3,1)*two-r2(3,5)*two+r1(3,1)*two)*qz+(-r3(10,6)+r2(3,4)-r1(3,8)+r0(4)) &
&                *xx+(-r3(6,5)*four+r2(5,3)*four)*xz+(-r3(3,4)+r2(1,2)-r1(3,6)+r0(2))*zz+( &
&                -r2(3,8)*two+r1(3,4)*two)*xxz+(-r2(5,7)*two+r1(1,3)*two)*xzz+(-r1(3,9)+r0(5) &
&                )*xxzz
      eri(4,1,3)=-r5(5)+r4(2,1)-r3(5,3)*three+r2(4,1)*three+(-r4(5,2)-r4(5,3)*two+r3(2,1) &
&                +r3(2,2)*two-r2(6,5)-r2(6,6)*two+r1(2,1)+r1(2,2)*two)*qx+(-r3(5,5)*two &
&                -r3(5,6)+r2(4,3)*two+r2(4,4))*xx+(-r2(6,8)+r1(2,4))*xxx
      eri(5,1,3)=-r5(6)+r4(3,1)-r3(6,3)*three+r2(5,1)*three+(-r4(6,2)-r4(6,3)*two+r3(3,1) &
&                +r3(3,2)*two-r2(3,5)-r2(3,6)*two+r1(3,1)+r1(3,2)*two)*qx+(-r4(3,2)+r3(1,1) &
&                -r2(5,5)*three+r1(1,1)*three)*qz+(-r3(6,5)*two-r3(6,6)+r2(5,3)*two+r2(5,4)) &
&                *xx+(-r3(3,4)-r3(3,5)*two+r2(1,2)+r2(1,3)*two-r1(3,6)-r1(3,7)*two+r0(2) &
&                +r0(3)*two)*xz+(-r2(3,8)+r1(3,4))*xxx+(-r2(5,7)*two-r2(5,8)+r1(1,3)*two &
&                +r1(1,4))*xxz+(-r1(3,9)+r0(5))*xxxz
      eri(6,1,3)=-r5(9)+r4(5,1)-r3(9,3)+r2(6,1)+(-r4(9,3)*two+r3(5,2)*two)*qx+(-r4(5,2) &
&                +r3(2,1)-r2(6,5)+r1(2,1))*qz+(-r3(9,6)+r2(6,4))*xx+(-r3(5,5)*two+r2(4,3)*two &
&                )*xz+(-r2(6,8)+r1(2,4))*xxz
      eri(1,2,3)=-r5(8)+r4(4,1)-r3(3,3)-r3(8,3)+r2(1,1)+r2(2,1)-r1(3,5)+r0(1)+( &
&                -r4(8,2)*two+r3(4,1)*two-r2(5,5)*two+r1(1,1)*two)*qx+(-r3(8,4)+r2(2,2) &
&                -r1(3,6)+r0(2))*xx
      eri(2,2,3)=-r5(17)+r4(11,1)-r3(8,3)*six+r2(2,1)*six-r1(3,5)*three+r0(1)*three
      eri(3,2,3)=-r5(19)+r4(13,1)-r3(8,3)-r3(10,3)+r2(2,1)+r2(3,1)-r1(3,5)+r0(1)+( &
&                -r4(13,2)*two+r3(8,1)*two-r2(3,5)*two+r1(3,1)*two)*qz+(-r3(8,4)+r2(2,2) &
&                -r1(3,6)+r0(2))*zz
      eri(4,2,3)=-r5(12)+r4(7,1)-r3(5,3)*three+r2(4,1)*three+(-r4(12,2)+r3(7,1) &
&                -r2(6,5)*three+r1(2,1)*three)*qx
      eri(5,2,3)=-r5(13)+r4(8,1)-r3(6,3)+r2(5,1)+(-r4(13,2)+r3(8,1)-r2(3,5)+r1(3,1))*qx+( &
&                -r4(8,2)+r3(4,1)-r2(5,5)+r1(1,1))*qz+(-r3(8,4)+r2(2,2)-r1(3,6)+r0(2))*xz
      eri(6,2,3)=-r5(18)+r4(12,1)-r3(9,3)*three+r2(6,1)*three+(-r4(12,2)+r3(7,1) &
&                -r2(6,5)*three+r1(2,1)*three)*qz
      eri(1,3,3)=-r5(10)+r4(6,1)-r3(3,3)-r3(10,3)+r2(1,1)+r2(3,1)-r1(3,5)+r0(1)+( &
&                -r4(10,2)*two+r3(6,1)*two-r2(5,5)*two+r1(1,1)*two)*qx+(-r4(6,3)*two &
&                +r3(3,2)*two-r2(3,6)*two+r1(3,2)*two)*qz+(-r3(10,4)+r2(3,2)-r1(3,6)+r0(2)) &
&                *xx+(-r3(6,5)*four+r2(5,3)*four)*xz+(-r3(3,6)+r2(1,4)-r1(3,8)+r0(4))*zz+( &
&                -r2(3,7)*two+r1(3,3)*two)*xxz+(-r2(5,8)*two+r1(1,4)*two)*xzz+(-r1(3,9)+r0(5) &
&                )*xxzz
      eri(2,3,3)=-r5(19)+r4(13,1)-r3(8,3)-r3(10,3)+r2(2,1)+r2(3,1)-r1(3,5)+r0(1)+( &
&                -r4(13,3)*two+r3(8,2)*two-r2(3,6)*two+r1(3,2)*two)*qz+(-r3(8,6)+r2(2,4) &
&                -r1(3,8)+r0(4))*zz
      eri(3,3,3)=-r5(21)+r4(15,1)-r3(10,3)*six+r2(3,1)*six-r1(3,5)*three+r0(1)*three+( &
&                -r4(15,2)*two-r4(15,3)*two+r3(10,1)*two+r3(10,2)*two-r2(3,5)*six-r2(3,6)*six &
&                +r1(3,1)*six+r1(3,2)*six)*qz+(-r3(10,4)-r3(10,5)*four-r3(10,6)+r2(3,2) &
&                +r2(3,3)*four+r2(3,4)-r1(3,6)-r1(3,7)*four-r1(3,8)+r0(2)+r0(3)*four+r0(4)) &
&                *zz+(-r2(3,7)*two-r2(3,8)*two+r1(3,3)*two+r1(3,4)*two)*zzz+(-r1(3,9)+r0(5)) &
&                *zzzz
      eri(4,3,3)=-r5(14)+r4(9,1)-r3(5,3)+r2(4,1)+(-r4(14,2)+r3(9,1)-r2(6,5)+r1(2,1))*qx+( &
&                -r4(9,3)*two+r3(5,2)*two)*qz+(-r3(9,5)*two+r2(6,3)*two)*xz+(-r3(5,6)+r2(4,4) &
&                )*zz+(-r2(6,8)+r1(2,4))*xzz
      eri(5,3,3)=-r5(15)+r4(10,1)-r3(6,3)*three+r2(5,1)*three+(-r4(15,2)+r3(10,1) &
&                -r2(3,5)*three+r1(3,1)*three)*qx+(-r4(10,2)-r4(10,3)*two+r3(6,1)+r3(6,2)*two &
&                -r2(5,5)-r2(5,6)*two+r1(1,1)+r1(1,2)*two)*qz+(-r3(10,4)-r3(10,5)*two+r2(3,2) &
&                +r2(3,3)*two-r1(3,6)-r1(3,7)*two+r0(2)+r0(3)*two)*xz+(-r3(6,5)*two-r3(6,6) &
&                +r2(5,3)*two+r2(5,4))*zz+(-r2(3,7)*two-r2(3,8)+r1(3,3)*two+r1(3,4))*xzz+( &
&                -r2(5,8)+r1(1,4))*zzz+(-r1(3,9)+r0(5))*xzzz
      eri(6,3,3)=-r5(20)+r4(14,1)-r3(9,3)*three+r2(6,1)*three+(-r4(14,2)-r4(14,3)*two &
&                +r3(9,1)+r3(9,2)*two-r2(6,5)-r2(6,6)*two+r1(2,1)+r1(2,2)*two)*qz+( &
&                -r3(9,5)*two-r3(9,6)+r2(6,3)*two+r2(6,4))*zz+(-r2(6,8)+r1(2,4))*zzz
      eri(1,4,3)=-r5(5)+r4(2,1)-r3(5,3)*three+r2(4,1)*three+(-r4(5,2)*two-r4(5,3) &
&                +r3(2,1)*two+r3(2,2)-r2(6,5)*two-r2(6,6)+r1(2,1)*two+r1(2,2))*qx+(-r3(5,4) &
&                -r3(5,5)*two+r2(4,2)+r2(4,3)*two)*xx+(-r2(6,7)+r1(2,3))*xxx
      eri(2,4,3)=-r5(12)+r4(7,1)-r3(5,3)*three+r2(4,1)*three+(-r4(12,3)+r3(7,2) &
&                -r2(6,6)*three+r1(2,2)*three)*qx
      eri(3,4,3)=-r5(14)+r4(9,1)-r3(5,3)+r2(4,1)+(-r4(14,3)+r3(9,2)-r2(6,6)+r1(2,2))*qx+( &
&                -r4(9,2)*two+r3(5,1)*two)*qz+(-r3(9,5)*two+r2(6,3)*two)*xz+(-r3(5,4)+r2(4,2) &
&                )*zz+(-r2(6,7)+r1(2,3))*xzz
      eri(4,4,3)=-r5(8)+r4(4,1)-r3(3,3)-r3(8,3)+r2(1,1)+r2(2,1)-r1(3,5)+r0(1)+(-r4(8,2) &
&                -r4(8,3)+r3(4,1)+r3(4,2)-r2(5,5)-r2(5,6)+r1(1,1)+r1(1,2))*qx+(-r3(8,5) &
&                +r2(2,3)-r1(3,7)+r0(3))*xx
      eri(5,4,3)=-r5(9)+r4(5,1)-r3(9,3)+r2(6,1)+(-r4(9,2)-r4(9,3)+r3(5,1)+r3(5,2))*qx+( &
&                -r4(5,2)+r3(2,1)-r2(6,5)+r1(2,1))*qz+(-r3(9,5)+r2(6,3))*xx+(-r3(5,4)-r3(5,5) &
&                +r2(4,2)+r2(4,3))*xz+(-r2(6,7)+r1(2,3))*xxz
      eri(6,4,3)=-r5(13)+r4(8,1)-r3(6,3)+r2(5,1)+(-r4(13,3)+r3(8,2)-r2(3,6)+r1(3,2))*qx+( &
&                -r4(8,2)+r3(4,1)-r2(5,5)+r1(1,1))*qz+(-r3(8,5)+r2(2,3)-r1(3,7)+r0(3))*xz
      eri(1,5,3)=-r5(6)+r4(3,1)-r3(6,3)*three+r2(5,1)*three+(-r4(6,2)*two-r4(6,3) &
&                +r3(3,1)*two+r3(3,2)-r2(3,5)*two-r2(3,6)+r1(3,1)*two+r1(3,2))*qx+(-r4(3,3) &
&                +r3(1,2)-r2(5,6)*three+r1(1,2)*three)*qz+(-r3(6,4)-r3(6,5)*two+r2(5,2) &
&                +r2(5,3)*two)*xx+(-r3(3,5)*two-r3(3,6)+r2(1,3)*two+r2(1,4)-r1(3,7)*two &
&                -r1(3,8)+r0(3)*two+r0(4))*xz+(-r2(3,7)+r1(3,3))*xxx+(-r2(5,7)-r2(5,8)*two &
&                +r1(1,3)+r1(1,4)*two)*xxz+(-r1(3,9)+r0(5))*xxxz
      eri(2,5,3)=-r5(13)+r4(8,1)-r3(6,3)+r2(5,1)+(-r4(13,3)+r3(8,2)-r2(3,6)+r1(3,2))*qx+( &
&                -r4(8,3)+r3(4,2)-r2(5,6)+r1(1,2))*qz+(-r3(8,6)+r2(2,4)-r1(3,8)+r0(4))*xz
      eri(3,5,3)=-r5(15)+r4(10,1)-r3(6,3)*three+r2(5,1)*three+(-r4(15,3)+r3(10,2) &
&                -r2(3,6)*three+r1(3,2)*three)*qx+(-r4(10,2)*two-r4(10,3)+r3(6,1)*two+r3(6,2) &
&                -r2(5,5)*two-r2(5,6)+r1(1,1)*two+r1(1,2))*qz+(-r3(10,5)*two-r3(10,6) &
&                +r2(3,3)*two+r2(3,4)-r1(3,7)*two-r1(3,8)+r0(3)*two+r0(4))*xz+(-r3(6,4) &
&                -r3(6,5)*two+r2(5,2)+r2(5,3)*two)*zz+(-r2(3,7)-r2(3,8)*two+r1(3,3) &
&                +r1(3,4)*two)*xzz+(-r2(5,7)+r1(1,3))*zzz+(-r1(3,9)+r0(5))*xzzz
      eri(4,5,3)=-r5(9)+r4(5,1)-r3(9,3)+r2(6,1)+(-r4(9,2)-r4(9,3)+r3(5,1)+r3(5,2))*qx+( &
&                -r4(5,3)+r3(2,2)-r2(6,6)+r1(2,2))*qz+(-r3(9,5)+r2(6,3))*xx+(-r3(5,5)-r3(5,6) &
&                +r2(4,3)+r2(4,4))*xz+(-r2(6,8)+r1(2,4))*xxz
      eri(5,5,3)=-r5(10)+r4(6,1)-r3(3,3)-r3(10,3)+r2(1,1)+r2(3,1)-r1(3,5)+r0(1)+(-r4(10,2) &
&                -r4(10,3)+r3(6,1)+r3(6,2)-r2(5,5)-r2(5,6)+r1(1,1)+r1(1,2))*qx+(-r4(6,2) &
&                -r4(6,3)+r3(3,1)+r3(3,2)-r2(3,5)-r2(3,6)+r1(3,1)+r1(3,2))*qz+(-r3(10,5) &
&                +r2(3,3)-r1(3,7)+r0(3))*xx+(-r3(6,4)-r3(6,5)*two-r3(6,6)+r2(5,2)+r2(5,3)*two &
&                +r2(5,4))*xz+(-r3(3,5)+r2(1,3)-r1(3,7)+r0(3))*zz+(-r2(3,7)-r2(3,8)+r1(3,3) &
&                +r1(3,4))*xxz+(-r2(5,7)-r2(5,8)+r1(1,3)+r1(1,4))*xzz+(-r1(3,9)+r0(5))*xxzz
      eri(6,5,3)=-r5(14)+r4(9,1)-r3(5,3)+r2(4,1)+(-r4(14,3)+r3(9,2)-r2(6,6)+r1(2,2))*qx+( &
&                -r4(9,2)-r4(9,3)+r3(5,1)+r3(5,2))*qz+(-r3(9,5)-r3(9,6)+r2(6,3)+r2(6,4))*xz+( &
&                -r3(5,5)+r2(4,3))*zz+(-r2(6,8)+r1(2,4))*xzz
      eri(1,6,3)=-r5(9)+r4(5,1)-r3(9,3)+r2(6,1)+(-r4(9,2)*two+r3(5,1)*two)*qx+(-r4(5,3) &
&                +r3(2,2)-r2(6,6)+r1(2,2))*qz+(-r3(9,4)+r2(6,2))*xx+(-r3(5,5)*two+r2(4,3)*two &
&                )*xz+(-r2(6,7)+r1(2,3))*xxz
      eri(2,6,3)=-r5(18)+r4(12,1)-r3(9,3)*three+r2(6,1)*three+(-r4(12,3)+r3(7,2) &
&                -r2(6,6)*three+r1(2,2)*three)*qz
      eri(3,6,3)=-r5(20)+r4(14,1)-r3(9,3)*three+r2(6,1)*three+(-r4(14,2)*two-r4(14,3) &
&                +r3(9,1)*two+r3(9,2)-r2(6,5)*two-r2(6,6)+r1(2,1)*two+r1(2,2))*qz+(-r3(9,4) &
&                -r3(9,5)*two+r2(6,2)+r2(6,3)*two)*zz+(-r2(6,7)+r1(2,3))*zzz
      eri(4,6,3)=-r5(13)+r4(8,1)-r3(6,3)+r2(5,1)+(-r4(13,2)+r3(8,1)-r2(3,5)+r1(3,1))*qx+( &
&                -r4(8,3)+r3(4,2)-r2(5,6)+r1(1,2))*qz+(-r3(8,5)+r2(2,3)-r1(3,7)+r0(3))*xz
      eri(5,6,3)=-r5(14)+r4(9,1)-r3(5,3)+r2(4,1)+(-r4(14,2)+r3(9,1)-r2(6,5)+r1(2,1))*qx+( &
&                -r4(9,2)-r4(9,3)+r3(5,1)+r3(5,2))*qz+(-r3(9,4)-r3(9,5)+r2(6,2)+r2(6,3))*xz+( &
&                -r3(5,5)+r2(4,3))*zz+(-r2(6,7)+r1(2,3))*xzz
      eri(6,6,3)=-r5(19)+r4(13,1)-r3(8,3)-r3(10,3)+r2(2,1)+r2(3,1)-r1(3,5)+r0(1)+( &
&                -r4(13,2)-r4(13,3)+r3(8,1)+r3(8,2)-r2(3,5)-r2(3,6)+r1(3,1)+r1(3,2))*qz+( &
&                -r3(8,5)+r2(2,3)-r1(3,7)+r0(3))*zz
!
      do l= 1,3
        rot2(1,l)= rot(1,l)*rot(1,l)
        rot2(2,l)= rot(2,l)*rot(2,l)
        rot2(3,l)= rot(3,l)*rot(3,l)
        rot2(4,l)= rot(1,l)*rot(2,l)*two
        rot2(5,l)= rot(1,l)*rot(3,l)*two
        rot2(6,l)= rot(2,l)*rot(3,l)*two
      enddo
      do l= 4,5
        rot2(1,l)= rot(1,1)*rot(1,l-2)
        rot2(2,l)= rot(2,1)*rot(2,l-2)
        rot2(3,l)= rot(3,1)*rot(3,l-2)
        rot2(4,l)= rot(1,1)*rot(2,l-2)+rot(2,1)*rot(1,l-2)
        rot2(5,l)= rot(1,1)*rot(3,l-2)+rot(3,1)*rot(1,l-2)
        rot2(6,l)= rot(2,1)*rot(3,l-2)+rot(3,1)*rot(2,l-2)
      enddo
      rot2(1,6)= rot(1,2)*rot(1,3)
      rot2(2,6)= rot(2,2)*rot(2,3)
      rot2(3,6)= rot(3,2)*rot(3,3)
      rot2(4,6)= rot(1,2)*rot(2,3)+rot(2,2)*rot(1,3)
      rot2(5,6)= rot(1,2)*rot(3,3)+rot(3,2)*rot(1,3)
      rot2(6,6)= rot(2,2)*rot(3,3)+rot(3,2)*rot(2,3)
      do l= 4,6
        do k= 1,6
          rot2(k,l)= rot2(k,l)*sqrt3
        enddo
      enddo
      do l= 1,6
        rot3(l,1)= rot2(l,3)-(rot2(l,1)+rot2(l,2))*half
        rot3(l,2)= rot2(l,5)
        rot3(l,3)= rot2(l,6)
        rot3(l,4)=(rot2(l,1)-rot2(l,2))*sqrt3h
        rot3(l,5)= rot2(l,4)
      enddo
!
      do k= 1,6
        do l= 1,6
          work(1)= eri(l,k,1)
          work(2)= eri(l,k,2)
          work(3)= eri(l,k,3)
          eri(l,k,1)= work(1)*rot(1,1)+work(2)*rot(2,1)+work(3)*rot(3,1)
          eri(l,k,2)= work(1)*rot(1,2)+work(2)*rot(2,2)+work(3)*rot(3,2)
          eri(l,k,3)= work(1)*rot(1,3)+work(2)*rot(2,3)+work(3)*rot(3,3)
        enddo
      enddo
!
      if(nbfijkl(3) == 6)then
        do j= 1,3
          do l= 1,6
            do k= 1,6
              work(k)= eri(l,k,j)
            enddo
            do k= 1,6
              eri(l,k,j)= work(1)*rot2(1,k)+work(2)*rot2(2,k)+work(3)*rot2(3,k) &
&                        +work(4)*rot2(4,k)+work(5)*rot2(5,k)+work(6)*rot2(6,k) 
            enddo
          enddo
        enddo
      else
        do j= 1,3
          do l= 1,6
            do k= 1,6
              work(k)= eri(l,k,j)
            enddo
            do k= 1,5
              eri(l,k,j)= work(1)*rot3(1,k)+work(2)*rot3(2,k)+work(3)*rot3(3,k) &
&                        +work(4)*rot3(4,k)+work(5)*rot3(5,k)+work(6)*rot3(6,k)
            enddo
          enddo
        enddo
      endif
!
      if(nbfijkl(4) == 6)then
        do j= 1,3
          do k= 1,nbfijkl(3)
            do l= 1,6
              phmdint(l,k,j,1)= eri(1,k,j)*rot2(1,l)+eri(2,k,j)*rot2(2,l)+eri(3,k,j)*rot2(3,l) &
&                              +eri(4,k,j)*rot2(4,l)+eri(5,k,j)*rot2(5,l)+eri(6,k,j)*rot2(6,l) 
            enddo
          enddo
        enddo
      else
        do j= 1,3
          do k= 1,nbfijkl(3)
            do l= 1,5
              phmdint(l,k,j,1)= eri(1,k,j)*rot3(1,l)+eri(2,k,j)*rot3(2,l)+eri(3,k,j)*rot3(3,l) &
&                              +eri(4,k,j)*rot3(4,l)+eri(5,k,j)*rot3(5,l)+eri(6,k,j)*rot3(6,l) 
            enddo
          enddo
        enddo
      endif
!
      return
end
