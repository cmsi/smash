!--------------------------------------------------------------------------
  subroutine int2dpds(phmdint,exfac1,exfac2,xyziq,xzkl,rot,nijkl,nbfijkl)
!--------------------------------------------------------------------------
!
! Calculate (dp|ds) integrals
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
      real(8) :: f0(2), f1(2,4), f2(3,4), f3(4,4), f4(5,2), f5(6), ftw(5,4)
      real(8) :: r0(6), r1(3,9), r2(6,10), r3(10,7), r4(15,3), r5(21)
      real(8) :: ex12, ex34, ex43, ex41, expq, expq2, expq4, ex3q, ex4q, c12, c34, zip
      real(8) :: xiq, yiq, ziq, xiq2, yiq2, xyiq, xiq4, yiq4, xyiq2, xypq2, zpq, zpq2, fac
      real(8) :: ex33q, ex34q, zjp, pmd, qmd, qmd2, qmd3, qmd3x, qmd3y, qmd3xy
      real(8) :: qx, qz, xx, xz, zz, xxx, xxz, xzz, zzz, eri(6,3,6), work(8)
      real(8) :: f1w(3,4), f2w(6,4), f3w(10,4), f4w(15,2), rot2(6,6), rot3(6,5)
!
! Zero-clear
!
      r0(1:6)     = zero
      r1(1:3,1:9) = zero
      r2(1:6,1:10)= zero
      r3(1:10,1:7)= zero
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
        f0(1:2)    = zero
        f1(1:2,1:4)= zero
        f2(1:3,1:4)= zero
        f3(1:4,1:4)= zero
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
          ftw(1,1)= zjp*zjp
          ftw(1,2)= pmd
          ftw(1,3)= pmd*zjp
          ftw(1,4)= pmd*pmd
          do i= 1,4
            ftw(2,i)= ftw(1,i)*zpq
            ftw(3,i)= ftw(1,i)*zpq2
            ftw(4,i)= ftw(1,i)*zpq2*zpq
            ftw(5,i)= ftw(1,i)*zpq2*zpq2
          enddo
          f0(1)= f0(1)+ft(0)*ftw(1,1)
          f0(2)= f0(2)+ft(0)*ftw(1,2)
          do i= 1,3
            f1(1,i)= f1(1,i)-ft(1)*ftw(1,i)
            f1(2,i)= f1(2,i)-ft(1)*ftw(2,i)
          enddo
          f1(1,4)= f1(1,4)-ft(1)*ftw(1,4)
          do i= 1,4
            f2(1,i)= f2(1,i)+ft(2)*ftw(1,i)
            f2(2,i)= f2(2,i)+ft(2)*ftw(2,i)
            f2(3,i)= f2(3,i)+ft(2)*ftw(3,i)
          enddo
          do i= 1,4
            do j= 1,4
              f3(j,i)= f3(j,i)-ft(3)*ftw(j,i)
            enddo
          enddo
          do i= 1,2
            do j= 1,5
              f4(j,i)= f4(j,i)+ft(4)*ftw(j,i+2)
            enddo
          enddo
          f5(1)= f5(1)-ft(5)*ftw(1,4)
          f5(2)= f5(2)-ft(5)*ftw(2,4)
          f5(3)= f5(3)-ft(5)*ftw(3,4)
          f5(4)= f5(4)-ft(5)*ftw(4,4)
          f5(5)= f5(5)-ft(5)*ftw(5,4)
          f5(6)= f5(6)-ft(5)*ftw(5,4)*zpq
        enddo
!
        qmd = ex43*c34
        qmd2= qmd*ex43
        qmd3= qmd*ex43*ex43
        ex33q= ex3q*ex3q
        ex34q= ex3q*ex4q
        qmd3x= qmd3*xiq
        qmd3y= qmd3*yiq
        qmd3xy=qmd3*xiq*yiq
!
        work( 1)= qmd*ex3q
        work( 2)= qmd*ex4q
        work( 3)= ex33q*ex4q*c34
        work( 4)= qmd2
        work( 5)= qmd*ex33q
        work( 6)= qmd*ex34q
        work( 7)= qmd2*ex3q
        work( 8)= qmd2*ex4q
!
        do i= 1,3
          r0(i  )= r0(i  )+f0(1)*work(i)
          r0(i+3)= r0(i+3)+f0(2)*work(i)
        enddo
!
        do i= 1,3
          f1w(1,i)= f1(1,i)*xiq
          f1w(2,i)= f1(1,i)*yiq
          f1w(3,i)= f1(2,i)
        enddo
        do i= 1,3
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,1)*work(i+3)
          enddo
        enddo
        do i= 4,6
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,2)*work(i)
          enddo
        enddo
        do i= 7,9
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,3)*work(i-6)
          enddo
        enddo
!
        do i= 1,4
          f2w(1,i)= f2(1,i)*xiq2+f1(1,i)
          f2w(2,i)= f2(1,i)*yiq2+f1(1,i)
          f2w(3,i)= f2(3,i)     +f1(1,i)
          f2w(4,i)= f2(1,i)*xyiq
          f2w(5,i)= f2(2,i)*xiq
          f2w(6,i)= f2(2,i)*yiq
        enddo
        do i= 1,2
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,1)*work(i+6)
          enddo
        enddo
        do i= 3,4
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,2)*work(i+4)
          enddo
        enddo
        do i= 5,7
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,3)*work(i-1)
          enddo
        enddo
        do i= 8,10
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,4)*work(i-7)
          enddo
        enddo
! 
        do i= 1,4
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
        do j= 1,10
          r3(j,1)= r3(j,1)+f3w(j,1)*qmd3
          r3(j,2)= r3(j,2)+f3w(j,2)*qmd3
        enddo
        do i= 3,4
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,3)*work(i+4)
          enddo
        enddo
        do i= 5,7
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,4)*work(i-1)
          enddo
        enddo
!
        do i= 1,2
          f4w( 1,i)=(f4(1,i)*xiq4 +f3(1,i+2)*xiq2*six           +f2(1,i+2)*three)
          f4w( 2,i)=(f4(1,i)*xiq2 +f3(1,i+2)*three                              )*xyiq
          f4w( 3,i)=(f4(2,i)*xiq2 +f3(2,i+2)*three                              )*xiq
          f4w( 4,i)=(f4(1,i)*xyiq2+f3(1,i+2)*xiq2+f3(1,i+2)*yiq2+f2(1,i+2)      )
          f4w( 5,i)=(f4(2,i)*xiq2 +f3(2,i+2)                                    )*yiq
          f4w( 6,i)=(f4(3,i)*xiq2 +f3(1,i+2)*xiq2+f3(3,i+2)     +f2(1,i+2)      )
          f4w( 7,i)=(f4(1,i)*yiq2 +f3(1,i+2)*three                              )*xyiq
          f4w( 8,i)=(f4(2,i)*yiq2 +f3(2,i+2)                                    )*xiq
          f4w( 9,i)=(f4(3,i)      +f3(1,i+2)                                    )*xyiq
          f4w(10,i)=(f4(4,i)      +f3(2,i+2)*three                              )*xiq
          f4w(11,i)=(f4(1,i)*yiq4 +f3(1,i+2)*yiq2*six           +f2(1,i+2)*three)
          f4w(12,i)=(f4(2,i)*yiq2 +f3(2,i+2)*three                              )*yiq
          f4w(13,i)=(f4(3,i)*yiq2 +f3(1,i+2)*yiq2+f3(3,i+2)     +f2(1,i+2)      )
          f4w(14,i)=(f4(4,i)      +f3(2,i+2)*three                              )*yiq
          f4w(15,i)=(f4(5,i)      +f3(3,i+2)*six                +f2(1,i+2)*three)
        enddo
        do j= 1,15
          r4(j,1)= r4(j,1)+f4w(j,1)*qmd3
          r4(j,2)= r4(j,2)+f4w(j,2)*work(7)
          r4(j,3)= r4(j,3)+f4w(j,2)*work(8)
        enddo
!
        r5( 1)= r5( 1)+(f5(1)*xiq4 +f4(1,2)*xiq2*ten               +f3(1,4)*p15  )*qmd3x
        r5( 2)= r5( 2)+(f5(1)*xiq4 +f4(1,2)*xiq2*six               +f3(1,4)*three)*qmd3y
        r5( 3)= r5( 3)+(f5(2)*xiq4 +f4(2,2)*xiq2*six               +f3(2,4)*three)*qmd3
        r5( 4)= r5( 4)+(f5(1)*xyiq2+f4(1,2)*xiq2+f4(1,2)*yiq2*three+f3(1,4)*three)*qmd3x
        r5( 5)= r5( 5)+(f5(2)*xiq2 +f4(2,2)*three                                )*qmd3xy
        r5( 6)= r5( 6)+(f5(3)*xiq2 +f4(1,2)*xiq2+f4(3,2)*three     +f3(1,4)*three)*qmd3x
        r5( 7)= r5( 7)+(f5(1)*xyiq2+f4(1,2)*xiq2*three+f4(1,2)*yiq2+f3(1,4)*three)*qmd3y
        r5( 8)= r5( 8)+(f5(2)*xyiq2+f4(2,2)*xiq2+f4(2,2)*yiq2      +f3(2,4)      )*qmd3
        r5( 9)= r5( 9)+(f5(3)*xiq2 +f4(1,2)*xiq2+f4(3,2)           +f3(1,4)      )*qmd3y
        r5(10)= r5(10)+(f5(4)*xiq2 +f4(2,2)*xiq2*three+f4(4,2)     +f3(2,4)*three)*qmd3
        r5(11)= r5(11)+(f5(1)*yiq4 +f4(1,2)*yiq2*six               +f3(1,4)*three)*qmd3x
        r5(12)= r5(12)+(f5(2)*yiq2 +f4(2,2)*three                                )*qmd3xy
        r5(13)= r5(13)+(f5(3)*yiq2 +f4(1,2)*yiq2+f4(3,2)           +f3(1,4)      )*qmd3x
        r5(14)= r5(14)+(f5(4)      +f4(2,2)*three                                )*qmd3xy
        r5(15)= r5(15)+(f5(5)      +f4(3,2)*six                    +f3(1,4)*three)*qmd3x
        r5(16)= r5(16)+(f5(1)*yiq4 +f4(1,2)*yiq2*ten               +f3(1,4)*p15  )*qmd3y
        r5(17)= r5(17)+(f5(2)*yiq4 +f4(2,2)*yiq2*six               +f3(2,4)*three)*qmd3
        r5(18)= r5(18)+(f5(3)*yiq2 +f4(1,2)*yiq2+f4(3,2)*three     +f3(1,4)*three)*qmd3y
        r5(19)= r5(19)+(f5(4)*yiq2 +f4(2,2)*yiq2*three+f4(4,2)     +f3(2,4)*three)*qmd3
        r5(20)= r5(20)+(f5(5)      +f4(3,2)*six                    +f3(1,4)*three)*qmd3y
        r5(21)= r5(21)+(f5(6)      +f4(4,2)*ten                    +f3(2,4)*p15  )*qmd3
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
      eri(1,1,1)= r5(1)+r3(1,2)+r3(1,5)*three+r1(1,4)*three+(+r4(1,2)*two+r4(1,3) &
&                +r2(1,3)*two+r2(1,4)+r2(1,8)*two+r2(1,9)+r0(4)*two+r0(5))*qx+(+r3(1,6) &
&                +r3(1,7)*two+r1(1,5)+r1(1,6)*two)*xx+(+r2(1,10)+r0(6))*xxx
      eri(2,1,1)= r5(4)+r3(4,2)+r3(1,5)+r1(1,4)+(+r4(4,3)+r2(2,4)+r2(1,9)+r0(5))*qx
      eri(3,1,1)= r5(6)+r3(6,2)+r3(1,5)+r1(1,4)+(+r4(6,3)+r2(3,4)+r2(1,9)+r0(5))*qx+( &
&                +r4(3,2)*two+r2(5,3)*two)*qz+(+r3(3,7)*two+r1(3,6)*two)*xz+(+r3(1,6)+r1(1,5) &
&                )*zz+(+r2(1,10)+r0(6))*xzz
      eri(4,1,1)= r5(2)+r3(2,2)+r3(2,5)+r1(2,4)+(+r4(2,2)+r4(2,3)+r2(4,3)+r2(4,4))*qx+( &
&                +r3(2,7)+r1(2,6))*xx
      eri(5,1,1)= r5(3)+r3(3,2)+r3(3,5)+r1(3,4)+(+r4(3,2)+r4(3,3)+r2(5,3)+r2(5,4))*qx+( &
&                +r4(1,2)+r2(1,3)+r2(1,8)+r0(4))*qz+(+r3(3,7)+r1(3,6))*xx+(+r3(1,6)+r3(1,7) &
&                +r1(1,5)+r1(1,6))*xz+(+r2(1,10)+r0(6))*xxz
      eri(6,1,1)= r5(5)+r3(5,2)+(+r4(5,3)+r2(6,4))*qx+(+r4(2,2)+r2(4,3))*qz+(+r3(2,7) &
&                +r1(2,6))*xz
      eri(1,2,1)= r5(2)+r3(2,2)+r3(2,5)+r1(2,4)+(+r4(2,2)*two+r2(4,3)*two)*qx+(+r3(2,6) &
&                +r1(2,5))*xx
      eri(2,2,1)= r5(7)+r3(7,2)+r3(2,5)*three+r1(2,4)*three
      eri(3,2,1)= r5(9)+r3(9,2)+r3(2,5)+r1(2,4)+(+r4(5,2)*two+r2(6,3)*two)*qz+(+r3(2,6) &
&                +r1(2,5))*zz
      eri(4,2,1)= r5(4)+r3(4,2)+r3(1,5)+r1(1,4)+(+r4(4,2)+r2(2,3)+r2(1,8)+r0(4))*qx
      eri(5,2,1)= r5(5)+r3(5,2)+(+r4(5,2)+r2(6,3))*qx+(+r4(2,2)+r2(4,3))*qz+(+r3(2,6) &
&                +r1(2,5))*xz
      eri(6,2,1)= r5(8)+r3(8,2)+r3(3,5)+r1(3,4)+(+r4(4,2)+r2(2,3)+r2(1,8)+r0(4))*qz
      eri(1,3,1)= r5(3)+r3(3,2)+r3(3,5)+r1(3,4)+(+r4(3,2)*two+r2(5,3)*two)*qx+(+r4(1,3) &
&                +r2(1,4)+r2(1,9)+r0(5))*qz+(+r3(3,6)+r1(3,5))*xx+(+r3(1,7)*two+r1(1,6)*two) &
&                *xz+(+r2(1,10)+r0(6))*xxz
      eri(2,3,1)= r5(8)+r3(8,2)+r3(3,5)+r1(3,4)+(+r4(4,3)+r2(2,4)+r2(1,9)+r0(5))*qz
      eri(3,3,1)= r5(10)+r3(10,2)+r3(3,5)*three+r1(3,4)*three+(+r4(6,2)*two+r4(6,3) &
&                +r2(3,3)*two+r2(3,4)+r2(1,8)*two+r2(1,9)+r0(4)*two+r0(5))*qz+(+r3(3,6) &
&                +r3(3,7)*two+r1(3,5)+r1(3,6)*two)*zz+(+r2(1,10)+r0(6))*zzz
      eri(4,3,1)= r5(5)+r3(5,2)+(+r4(5,2)+r2(6,3))*qx+(+r4(2,3)+r2(4,4))*qz+(+r3(2,7) &
&                +r1(2,6))*xz
      eri(5,3,1)= r5(6)+r3(6,2)+r3(1,5)+r1(1,4)+(+r4(6,2)+r2(3,3)+r2(1,8)+r0(4))*qx+( &
&                +r4(3,2)+r4(3,3)+r2(5,3)+r2(5,4))*qz+(+r3(3,6)+r3(3,7)+r1(3,5)+r1(3,6))*xz+( &
&                +r3(1,7)+r1(1,6))*zz+(+r2(1,10)+r0(6))*xzz
      eri(6,3,1)= r5(9)+r3(9,2)+r3(2,5)+r1(2,4)+(+r4(5,2)+r4(5,3)+r2(6,3)+r2(6,4))*qz+( &
&                +r3(2,7)+r1(2,6))*zz
      eri(1,1,2)= r5(4)+r3(1,2)+r3(4,5)*three+r1(1,4)*three+(+r4(4,2)*two+r4(4,3) &
&                +r2(1,3)*two+r2(1,4)+r2(2,8)*two+r2(2,9)+r0(4)*two+r0(5))*qx+(+r3(4,6) &
&                +r3(4,7)*two+r1(1,5)+r1(1,6)*two)*xx+(+r2(2,10)+r0(6))*xxx
      eri(2,1,2)= r5(11)+r3(4,2)+r3(4,5)+r1(1,4)+(+r4(11,3)+r2(2,4)+r2(2,9)+r0(5))*qx
      eri(3,1,2)= r5(13)+r3(6,2)+r3(4,5)+r1(1,4)+(+r4(13,3)+r2(3,4)+r2(2,9)+r0(5))*qx+( &
&                +r4(8,2)*two+r2(5,3)*two)*qz+(+r3(8,7)*two+r1(3,6)*two)*xz+(+r3(4,6)+r1(1,5) &
&                )*zz+(+r2(2,10)+r0(6))*xzz
      eri(4,1,2)= r5(7)+r3(2,2)+r3(7,5)+r1(2,4)+(+r4(7,2)+r4(7,3)+r2(4,3)+r2(4,4))*qx+( &
&                +r3(7,7)+r1(2,6))*xx
      eri(5,1,2)= r5(8)+r3(3,2)+r3(8,5)+r1(3,4)+(+r4(8,2)+r4(8,3)+r2(5,3)+r2(5,4))*qx+( &
&                +r4(4,2)+r2(1,3)+r2(2,8)+r0(4))*qz+(+r3(8,7)+r1(3,6))*xx+(+r3(4,6)+r3(4,7) &
&                +r1(1,5)+r1(1,6))*xz+(+r2(2,10)+r0(6))*xxz
      eri(6,1,2)= r5(12)+r3(5,2)+(+r4(12,3)+r2(6,4))*qx+(+r4(7,2)+r2(4,3))*qz+(+r3(7,7) &
&                +r1(2,6))*xz
      eri(1,2,2)= r5(7)+r3(2,2)+r3(7,5)+r1(2,4)+(+r4(7,2)*two+r2(4,3)*two)*qx+(+r3(7,6) &
&                +r1(2,5))*xx
      eri(2,2,2)= r5(16)+r3(7,2)+r3(7,5)*three+r1(2,4)*three
      eri(3,2,2)= r5(18)+r3(9,2)+r3(7,5)+r1(2,4)+(+r4(12,2)*two+r2(6,3)*two)*qz+(+r3(7,6) &
&                +r1(2,5))*zz
      eri(4,2,2)= r5(11)+r3(4,2)+r3(4,5)+r1(1,4)+(+r4(11,2)+r2(2,3)+r2(2,8)+r0(4))*qx
      eri(5,2,2)= r5(12)+r3(5,2)+(+r4(12,2)+r2(6,3))*qx+(+r4(7,2)+r2(4,3))*qz+(+r3(7,6) &
&                +r1(2,5))*xz
      eri(6,2,2)= r5(17)+r3(8,2)+r3(8,5)+r1(3,4)+(+r4(11,2)+r2(2,3)+r2(2,8)+r0(4))*qz
      eri(1,3,2)= r5(8)+r3(3,2)+r3(8,5)+r1(3,4)+(+r4(8,2)*two+r2(5,3)*two)*qx+(+r4(4,3) &
&                +r2(1,4)+r2(2,9)+r0(5))*qz+(+r3(8,6)+r1(3,5))*xx+(+r3(4,7)*two+r1(1,6)*two) &
&                *xz+(+r2(2,10)+r0(6))*xxz
      eri(2,3,2)= r5(17)+r3(8,2)+r3(8,5)+r1(3,4)+(+r4(11,3)+r2(2,4)+r2(2,9)+r0(5))*qz
      eri(3,3,2)= r5(19)+r3(10,2)+r3(8,5)*three+r1(3,4)*three+(+r4(13,2)*two+r4(13,3) &
&                +r2(3,3)*two+r2(3,4)+r2(2,8)*two+r2(2,9)+r0(4)*two+r0(5))*qz+(+r3(8,6) &
&                +r3(8,7)*two+r1(3,5)+r1(3,6)*two)*zz+(+r2(2,10)+r0(6))*zzz
      eri(4,3,2)= r5(12)+r3(5,2)+(+r4(12,2)+r2(6,3))*qx+(+r4(7,3)+r2(4,4))*qz+(+r3(7,7) &
&                +r1(2,6))*xz
      eri(5,3,2)= r5(13)+r3(6,2)+r3(4,5)+r1(1,4)+(+r4(13,2)+r2(3,3)+r2(2,8)+r0(4))*qx+( &
&                +r4(8,2)+r4(8,3)+r2(5,3)+r2(5,4))*qz+(+r3(8,6)+r3(8,7)+r1(3,5)+r1(3,6))*xz+( &
&                +r3(4,7)+r1(1,6))*zz+(+r2(2,10)+r0(6))*xzz
      eri(6,3,2)= r5(18)+r3(9,2)+r3(7,5)+r1(2,4)+(+r4(12,2)+r4(12,3)+r2(6,3)+r2(6,4))*qz+( &
&                +r3(7,7)+r1(2,6))*zz
      eri(1,1,3)= r5(6)-r4(3,1)*two+r3(1,1)+r3(1,2)+r3(6,5)*three-r2(5,5)*six &
&                +r1(1,1)*three+r1(1,4)*three+(+r4(6,2)*two+r4(6,3)-r3(3,3)*four-r3(3,4)*two &
&                +r2(1,1)*two+r2(1,2)+r2(1,3)*two+r2(1,4)+r2(3,8)*two+r2(3,9)-r1(3,7)*four &
&                -r1(3,8)*two+r0(1)*two+r0(2)+r0(4)*two+r0(5))*qx+(+r3(6,6)+r3(6,7)*two &
&                -r2(5,6)*two-r2(5,7)*four+r1(1,2)+r1(1,3)*two+r1(1,5)+r1(1,6)*two)*xx+( &
&                +r2(3,10)-r1(3,9)*two+r0(3)+r0(6))*xxx
      eri(2,1,3)= r5(13)-r4(8,1)*two+r3(4,1)+r3(4,2)+r3(6,5)-r2(5,5)*two+r1(1,1)+r1(1,4)+( &
&                +r4(13,3)-r3(8,4)*two+r2(2,2)+r2(2,4)+r2(3,9)-r1(3,8)*two+r0(2)+r0(5))*qx
      eri(3,1,3)= r5(15)-r4(10,1)*two+r3(6,1)+r3(6,2)+r3(6,5)-r2(5,5)*two+r1(1,1)+r1(1,4) &
&                +(+r4(15,3)-r3(10,4)*two+r2(3,2)+r2(3,4)+r2(3,9)-r1(3,8)*two+r0(2)+r0(5))*qx &
&                +(+r4(10,2)*two-r3(6,3)*four+r2(5,1)*two+r2(5,3)*two)*qz+(+r3(10,7)*two &
&                -r2(3,7)*four+r1(3,3)*two+r1(3,6)*two)*xz+(+r3(6,6)-r2(5,6)*two+r1(1,2) &
&                +r1(1,5))*zz+(+r2(3,10)-r1(3,9)*two+r0(3)+r0(6))*xzz
      eri(4,1,3)= r5(9)-r4(5,1)*two+r3(2,1)+r3(2,2)+r3(9,5)-r2(6,5)*two+r1(2,1)+r1(2,4)+( &
&                +r4(9,2)+r4(9,3)-r3(5,3)*two-r3(5,4)*two+r2(4,1)+r2(4,2)+r2(4,3)+r2(4,4))*qx &
&                +(+r3(9,7)-r2(6,7)*two+r1(2,3)+r1(2,6))*xx
      eri(5,1,3)= r5(10)-r4(6,1)*two+r3(3,1)+r3(3,2)+r3(10,5)-r2(3,5)*two+r1(3,1)+r1(3,4) &
&                +(+r4(10,2)+r4(10,3)-r3(6,3)*two-r3(6,4)*two+r2(5,1)+r2(5,2)+r2(5,3)+r2(5,4) &
&                )*qx+(+r4(6,2)-r3(3,3)*two+r2(1,1)+r2(1,3)+r2(3,8)-r1(3,7)*two+r0(1)+r0(4)) &
&                *qz+(+r3(10,7)-r2(3,7)*two+r1(3,3)+r1(3,6))*xx+(+r3(6,6)+r3(6,7)-r2(5,6)*two &
&                -r2(5,7)*two+r1(1,2)+r1(1,3)+r1(1,5)+r1(1,6))*xz+(+r2(3,10)-r1(3,9)*two &
&                +r0(3)+r0(6))*xxz
      eri(6,1,3)= r5(14)-r4(9,1)*two+r3(5,1)+r3(5,2)+(+r4(14,3)-r3(9,4)*two+r2(6,2) &
&                +r2(6,4))*qx+(+r4(9,2)-r3(5,3)*two+r2(4,1)+r2(4,3))*qz+(+r3(9,7)-r2(6,7)*two &
&                +r1(2,3)+r1(2,6))*xz
      eri(1,2,3)= r5(9)-r4(5,1)*two+r3(2,1)+r3(2,2)+r3(9,5)-r2(6,5)*two+r1(2,1)+r1(2,4)+( &
&                +r4(9,2)*two-r3(5,3)*four+r2(4,1)*two+r2(4,3)*two)*qx+(+r3(9,6)-r2(6,6)*two &
&                +r1(2,2)+r1(2,5))*xx
      eri(2,2,3)= r5(18)-r4(12,1)*two+r3(7,1)+r3(7,2)+r3(9,5)*three-r2(6,5)*six &
&                +r1(2,1)*three+r1(2,4)*three
      eri(3,2,3)= r5(20)-r4(14,1)*two+r3(9,1)+r3(9,2)+r3(9,5)-r2(6,5)*two+r1(2,1)+r1(2,4) &
&                +(+r4(14,2)*two-r3(9,3)*four+r2(6,1)*two+r2(6,3)*two)*qz+(+r3(9,6) &
&                -r2(6,6)*two+r1(2,2)+r1(2,5))*zz
      eri(4,2,3)= r5(13)-r4(8,1)*two+r3(4,1)+r3(4,2)+r3(6,5)-r2(5,5)*two+r1(1,1)+r1(1,4)+( &
&                +r4(13,2)-r3(8,3)*two+r2(2,1)+r2(2,3)+r2(3,8)-r1(3,7)*two+r0(1)+r0(4))*qx
      eri(5,2,3)= r5(14)-r4(9,1)*two+r3(5,1)+r3(5,2)+(+r4(14,2)-r3(9,3)*two+r2(6,1) &
&                +r2(6,3))*qx+(+r4(9,2)-r3(5,3)*two+r2(4,1)+r2(4,3))*qz+(+r3(9,6)-r2(6,6)*two &
&                +r1(2,2)+r1(2,5))*xz
      eri(6,2,3)= r5(19)-r4(13,1)*two+r3(8,1)+r3(8,2)+r3(10,5)-r2(3,5)*two+r1(3,1)+r1(3,4) &
&                +(+r4(13,2)-r3(8,3)*two+r2(2,1)+r2(2,3)+r2(3,8)-r1(3,7)*two+r0(1)+r0(4))*qz
      eri(1,3,3)= r5(10)-r4(6,1)*two+r3(3,1)+r3(3,2)+r3(10,5)-r2(3,5)*two+r1(3,1)+r1(3,4) &
&                +(+r4(10,2)*two-r3(6,3)*four+r2(5,1)*two+r2(5,3)*two)*qx+(+r4(6,3) &
&                -r3(3,4)*two+r2(1,2)+r2(1,4)+r2(3,9)-r1(3,8)*two+r0(2)+r0(5))*qz+(+r3(10,6) &
&                -r2(3,6)*two+r1(3,2)+r1(3,5))*xx+(+r3(6,7)*two-r2(5,7)*four+r1(1,3)*two &
&                +r1(1,6)*two)*xz+(+r2(3,10)-r1(3,9)*two+r0(3)+r0(6))*xxz
      eri(2,3,3)= r5(19)-r4(13,1)*two+r3(8,1)+r3(8,2)+r3(10,5)-r2(3,5)*two+r1(3,1)+r1(3,4) &
&                +(+r4(13,3)-r3(8,4)*two+r2(2,2)+r2(2,4)+r2(3,9)-r1(3,8)*two+r0(2)+r0(5))*qz
      eri(3,3,3)= r5(21)-r4(15,1)*two+r3(10,1)+r3(10,2)+r3(10,5)*three-r2(3,5)*six &
&                +r1(3,1)*three+r1(3,4)*three+(+r4(15,2)*two+r4(15,3)-r3(10,3)*four &
&                -r3(10,4)*two+r2(3,1)*two+r2(3,2)+r2(3,3)*two+r2(3,4)+r2(3,8)*two+r2(3,9) &
&                -r1(3,7)*four-r1(3,8)*two+r0(1)*two+r0(2)+r0(4)*two+r0(5))*qz+(+r3(10,6) &
&                +r3(10,7)*two-r2(3,6)*two-r2(3,7)*four+r1(3,2)+r1(3,3)*two+r1(3,5) &
&                +r1(3,6)*two)*zz+(+r2(3,10)-r1(3,9)*two+r0(3)+r0(6))*zzz
      eri(4,3,3)= r5(14)-r4(9,1)*two+r3(5,1)+r3(5,2)+(+r4(14,2)-r3(9,3)*two+r2(6,1) &
&                +r2(6,3))*qx+(+r4(9,3)-r3(5,4)*two+r2(4,2)+r2(4,4))*qz+(+r3(9,7)-r2(6,7)*two &
&                +r1(2,3)+r1(2,6))*xz
      eri(5,3,3)= r5(15)-r4(10,1)*two+r3(6,1)+r3(6,2)+r3(6,5)-r2(5,5)*two+r1(1,1)+r1(1,4) &
&                +(+r4(15,2)-r3(10,3)*two+r2(3,1)+r2(3,3)+r2(3,8)-r1(3,7)*two+r0(1)+r0(4))*qx &
&                +(+r4(10,2)+r4(10,3)-r3(6,3)*two-r3(6,4)*two+r2(5,1)+r2(5,2)+r2(5,3)+r2(5,4) &
&                )*qz+(+r3(10,6)+r3(10,7)-r2(3,6)*two-r2(3,7)*two+r1(3,2)+r1(3,3)+r1(3,5) &
&                +r1(3,6))*xz+(+r3(6,7)-r2(5,7)*two+r1(1,3)+r1(1,6))*zz+(+r2(3,10) &
&                -r1(3,9)*two+r0(3)+r0(6))*xzz
      eri(6,3,3)= r5(20)-r4(14,1)*two+r3(9,1)+r3(9,2)+r3(9,5)-r2(6,5)*two+r1(2,1)+r1(2,4) &
&                +(+r4(14,2)+r4(14,3)-r3(9,3)*two-r3(9,4)*two+r2(6,1)+r2(6,2)+r2(6,3)+r2(6,4) &
&                )*qz+(+r3(9,7)-r2(6,7)*two+r1(2,3)+r1(2,6))*zz
      eri(1,1,4)= r5(2)+r3(2,5)*three+(+r4(2,2)*two+r4(2,3)+r2(4,8)*two+r2(4,9))*qx+( &
&                +r3(2,6)+r3(2,7)*two)*xx+r2(4,10)*xxx
      eri(2,1,4)= r5(7)+r3(2,5)+(+r4(7,3)+r2(4,9))*qx
      eri(3,1,4)= r5(9)+r3(2,5)+(+r4(9,3)+r2(4,9))*qx+r4(5,2)*two*qz+r3(5,7)*two*xz &
&                +r3(2,6)*zz+r2(4,10)*xzz
      eri(4,1,4)= r5(4)+r3(4,5)+(+r4(4,2)+r4(4,3))*qx+r3(4,7)*xx
      eri(5,1,4)= r5(5)+r3(5,5)+(+r4(5,2)+r4(5,3))*qx+(+r4(2,2)+r2(4,8))*qz+r3(5,7)*xx+( &
&                +r3(2,6)+r3(2,7))*xz+r2(4,10)*xxz
      eri(6,1,4)= r5(8)+r4(8,3)*qx+r4(4,2)*qz+r3(4,7)*xz
      eri(1,2,4)= r5(4)+r3(4,5)+r4(4,2)*two*qx+r3(4,6)*xx
      eri(2,2,4)= r5(11)+r3(4,5)*three
      eri(3,2,4)= r5(13)+r3(4,5)+r4(8,2)*two*qz+r3(4,6)*zz
      eri(4,2,4)= r5(7)+r3(2,5)+(+r4(7,2)+r2(4,8))*qx
      eri(5,2,4)= r5(8)+r4(8,2)*qx+r4(4,2)*qz+r3(4,6)*xz
      eri(6,2,4)= r5(12)+r3(5,5)+(+r4(7,2)+r2(4,8))*qz
      eri(1,3,4)= r5(5)+r3(5,5)+r4(5,2)*two*qx+(+r4(2,3)+r2(4,9))*qz+r3(5,6)*xx &
&                +r3(2,7)*two*xz+r2(4,10)*xxz
      eri(2,3,4)= r5(12)+r3(5,5)+(+r4(7,3)+r2(4,9))*qz
      eri(3,3,4)= r5(14)+r3(5,5)*three+(+r4(9,2)*two+r4(9,3)+r2(4,8)*two+r2(4,9))*qz+( &
&                +r3(5,6)+r3(5,7)*two)*zz+r2(4,10)*zzz
      eri(4,3,4)= r5(8)+r4(8,2)*qx+r4(4,3)*qz+r3(4,7)*xz
      eri(5,3,4)= r5(9)+r3(2,5)+(+r4(9,2)+r2(4,8))*qx+(+r4(5,2)+r4(5,3))*qz+(+r3(5,6) &
&                +r3(5,7))*xz+r3(2,7)*zz+r2(4,10)*xzz
      eri(6,3,4)= r5(13)+r3(4,5)+(+r4(8,2)+r4(8,3))*qz+r3(4,7)*zz
      eri(1,1,5)= r5(3)-r4(1,1)+r3(3,5)*three-r2(1,5)*three+(+r4(3,2)*two+r4(3,3) &
&                -r3(1,3)*two-r3(1,4)+r2(5,8)*two+r2(5,9)-r1(1,7)*two-r1(1,8))*qx+(+r3(3,6) &
&                +r3(3,7)*two-r2(1,6)-r2(1,7)*two)*xx+(+r2(5,10)-r1(1,9))*xxx
      eri(2,1,5)= r5(8)-r4(4,1)+r3(3,5)-r2(1,5)+(+r4(8,3)-r3(4,4)+r2(5,9)-r1(1,8))*qx
      eri(3,1,5)= r5(10)-r4(6,1)+r3(3,5)-r2(1,5)+(+r4(10,3)-r3(6,4)+r2(5,9)-r1(1,8))*qx+( &
&                +r4(6,2)*two-r3(3,3)*two)*qz+(+r3(6,7)*two-r2(5,7)*two)*xz+(+r3(3,6)-r2(1,6) &
&                )*zz+(+r2(5,10)-r1(1,9))*xzz
      eri(4,1,5)= r5(5)-r4(2,1)+r3(5,5)-r2(4,5)+(+r4(5,2)+r4(5,3)-r3(2,3)-r3(2,4))*qx+( &
&                +r3(5,7)-r2(4,7))*xx
      eri(5,1,5)= r5(6)-r4(3,1)+r3(6,5)-r2(5,5)+(+r4(6,2)+r4(6,3)-r3(3,3)-r3(3,4))*qx+( &
&                +r4(3,2)-r3(1,3)+r2(5,8)-r1(1,7))*qz+(+r3(6,7)-r2(5,7))*xx+(+r3(3,6)+r3(3,7) &
&                -r2(1,6)-r2(1,7))*xz+(+r2(5,10)-r1(1,9))*xxz
      eri(6,1,5)= r5(9)-r4(5,1)+(+r4(9,3)-r3(5,4))*qx+(+r4(5,2)-r3(2,3))*qz+(+r3(5,7) &
&                -r2(4,7))*xz
      eri(1,2,5)= r5(5)-r4(2,1)+r3(5,5)-r2(4,5)+(+r4(5,2)*two-r3(2,3)*two)*qx+(+r3(5,6) &
&                -r2(4,6))*xx
      eri(2,2,5)= r5(12)-r4(7,1)+r3(5,5)*three-r2(4,5)*three
      eri(3,2,5)= r5(14)-r4(9,1)+r3(5,5)-r2(4,5)+(+r4(9,2)*two-r3(5,3)*two)*qz+(+r3(5,6) &
&                -r2(4,6))*zz
      eri(4,2,5)= r5(8)-r4(4,1)+r3(3,5)-r2(1,5)+(+r4(8,2)-r3(4,3)+r2(5,8)-r1(1,7))*qx
      eri(5,2,5)= r5(9)-r4(5,1)+(+r4(9,2)-r3(5,3))*qx+(+r4(5,2)-r3(2,3))*qz+(+r3(5,6) &
&                -r2(4,6))*xz
      eri(6,2,5)= r5(13)-r4(8,1)+r3(6,5)-r2(5,5)+(+r4(8,2)-r3(4,3)+r2(5,8)-r1(1,7))*qz
      eri(1,3,5)= r5(6)-r4(3,1)+r3(6,5)-r2(5,5)+(+r4(6,2)*two-r3(3,3)*two)*qx+(+r4(3,3) &
&                -r3(1,4)+r2(5,9)-r1(1,8))*qz+(+r3(6,6)-r2(5,6))*xx+(+r3(3,7)*two-r2(1,7)*two &
&                )*xz+(+r2(5,10)-r1(1,9))*xxz
      eri(2,3,5)= r5(13)-r4(8,1)+r3(6,5)-r2(5,5)+(+r4(8,3)-r3(4,4)+r2(5,9)-r1(1,8))*qz
      eri(3,3,5)= r5(15)-r4(10,1)+r3(6,5)*three-r2(5,5)*three+(+r4(10,2)*two+r4(10,3) &
&                -r3(6,3)*two-r3(6,4)+r2(5,8)*two+r2(5,9)-r1(1,7)*two-r1(1,8))*qz+(+r3(6,6) &
&                +r3(6,7)*two-r2(5,6)-r2(5,7)*two)*zz+(+r2(5,10)-r1(1,9))*zzz
      eri(4,3,5)= r5(9)-r4(5,1)+(+r4(9,2)-r3(5,3))*qx+(+r4(5,3)-r3(2,4))*qz+(+r3(5,7) &
&                -r2(4,7))*xz
      eri(5,3,5)= r5(10)-r4(6,1)+r3(3,5)-r2(1,5)+(+r4(10,2)-r3(6,3)+r2(5,8)-r1(1,7))*qx+( &
&                +r4(6,2)+r4(6,3)-r3(3,3)-r3(3,4))*qz+(+r3(6,6)+r3(6,7)-r2(5,6)-r2(5,7))*xz+( &
&                +r3(3,7)-r2(1,7))*zz+(+r2(5,10)-r1(1,9))*xzz
      eri(6,3,5)= r5(14)-r4(9,1)+r3(5,5)-r2(4,5)+(+r4(9,2)+r4(9,3)-r3(5,3)-r3(5,4))*qz+( &
&                +r3(5,7)-r2(4,7))*zz
      eri(1,1,6)= r5(5)-r4(2,1)+r3(5,5)*three-r2(4,5)*three+(+r4(5,2)*two+r4(5,3) &
&                -r3(2,3)*two-r3(2,4)+r2(6,8)*two+r2(6,9)-r1(2,7)*two-r1(2,8))*qx+(+r3(5,6) &
&                +r3(5,7)*two-r2(4,6)-r2(4,7)*two)*xx+(+r2(6,10)-r1(2,9))*xxx
      eri(2,1,6)= r5(12)-r4(7,1)+r3(5,5)-r2(4,5)+(+r4(12,3)-r3(7,4)+r2(6,9)-r1(2,8))*qx
      eri(3,1,6)= r5(14)-r4(9,1)+r3(5,5)-r2(4,5)+(+r4(14,3)-r3(9,4)+r2(6,9)-r1(2,8))*qx+( &
&                +r4(9,2)*two-r3(5,3)*two)*qz+(+r3(9,7)*two-r2(6,7)*two)*xz+(+r3(5,6)-r2(4,6) &
&                )*zz+(+r2(6,10)-r1(2,9))*xzz
      eri(4,1,6)= r5(8)-r4(4,1)+r3(8,5)-r2(2,5)+(+r4(8,2)+r4(8,3)-r3(4,3)-r3(4,4))*qx+( &
&                +r3(8,7)-r2(2,7))*xx
      eri(5,1,6)= r5(9)-r4(5,1)+r3(9,5)-r2(6,5)+(+r4(9,2)+r4(9,3)-r3(5,3)-r3(5,4))*qx+( &
&                +r4(5,2)-r3(2,3)+r2(6,8)-r1(2,7))*qz+(+r3(9,7)-r2(6,7))*xx+(+r3(5,6)+r3(5,7) &
&                -r2(4,6)-r2(4,7))*xz+(+r2(6,10)-r1(2,9))*xxz
      eri(6,1,6)= r5(13)-r4(8,1)+(+r4(13,3)-r3(8,4))*qx+(+r4(8,2)-r3(4,3))*qz+(+r3(8,7) &
&                -r2(2,7))*xz
      eri(1,2,6)= r5(8)-r4(4,1)+r3(8,5)-r2(2,5)+(+r4(8,2)*two-r3(4,3)*two)*qx+(+r3(8,6) &
&                -r2(2,6))*xx
      eri(2,2,6)= r5(17)-r4(11,1)+r3(8,5)*three-r2(2,5)*three
      eri(3,2,6)= r5(19)-r4(13,1)+r3(8,5)-r2(2,5)+(+r4(13,2)*two-r3(8,3)*two)*qz+(+r3(8,6) &
&                -r2(2,6))*zz
      eri(4,2,6)= r5(12)-r4(7,1)+r3(5,5)-r2(4,5)+(+r4(12,2)-r3(7,3)+r2(6,8)-r1(2,7))*qx
      eri(5,2,6)= r5(13)-r4(8,1)+(+r4(13,2)-r3(8,3))*qx+(+r4(8,2)-r3(4,3))*qz+(+r3(8,6) &
&                -r2(2,6))*xz
      eri(6,2,6)= r5(18)-r4(12,1)+r3(9,5)-r2(6,5)+(+r4(12,2)-r3(7,3)+r2(6,8)-r1(2,7))*qz
      eri(1,3,6)= r5(9)-r4(5,1)+r3(9,5)-r2(6,5)+(+r4(9,2)*two-r3(5,3)*two)*qx+(+r4(5,3) &
&                -r3(2,4)+r2(6,9)-r1(2,8))*qz+(+r3(9,6)-r2(6,6))*xx+(+r3(5,7)*two-r2(4,7)*two &
&                )*xz+(+r2(6,10)-r1(2,9))*xxz
      eri(2,3,6)= r5(18)-r4(12,1)+r3(9,5)-r2(6,5)+(+r4(12,3)-r3(7,4)+r2(6,9)-r1(2,8))*qz
      eri(3,3,6)= r5(20)-r4(14,1)+r3(9,5)*three-r2(6,5)*three+(+r4(14,2)*two+r4(14,3) &
&                -r3(9,3)*two-r3(9,4)+r2(6,8)*two+r2(6,9)-r1(2,7)*two-r1(2,8))*qz+(+r3(9,6) &
&                +r3(9,7)*two-r2(6,6)-r2(6,7)*two)*zz+(+r2(6,10)-r1(2,9))*zzz
      eri(4,3,6)= r5(13)-r4(8,1)+(+r4(13,2)-r3(8,3))*qx+(+r4(8,3)-r3(4,4))*qz+(+r3(8,7) &
&                -r2(2,7))*xz
      eri(5,3,6)= r5(14)-r4(9,1)+r3(5,5)-r2(4,5)+(+r4(14,2)-r3(9,3)+r2(6,8)-r1(2,7))*qx+( &
&                +r4(9,2)+r4(9,3)-r3(5,3)-r3(5,4))*qz+(+r3(9,6)+r3(9,7)-r2(6,6)-r2(6,7))*xz+( &
&                +r3(5,7)-r2(4,7))*zz+(+r2(6,10)-r1(2,9))*xzz
      eri(6,3,6)= r5(19)-r4(13,1)+r3(8,5)-r2(2,5)+(+r4(13,2)+r4(13,3)-r3(8,3)-r3(8,4))*qz &
&                +(+r3(8,7)-r2(2,7))*zz
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
      if(nbfijkl(2) == 6)then
        do k= 1,3
          do l= 1,6
            do j= 1,6
              work(j)= eri(l,k,j)
            enddo
            do j= 1,6
              eri(l,k,j)= work(1)*rot2(1,j)+work(2)*rot2(2,j)+work(3)*rot2(3,j) &
&                        +work(4)*rot2(4,j)+work(5)*rot2(5,j)+work(6)*rot2(6,j) 
            enddo
          enddo
        enddo
      else
        do k= 1,3
          do l= 1,6
            do j= 1,6
              work(j)= eri(l,k,j)
            enddo
            do j= 1,5
              eri(l,k,j)= work(1)*rot3(1,j)+work(2)*rot3(2,j)+work(3)*rot3(3,j) &
&                        +work(4)*rot3(4,j)+work(5)*rot3(5,j)+work(6)*rot3(6,j)
            enddo
          enddo
        enddo
      endif
!
      do j= 1,nbfijkl(2)
        do l= 1,6
          work(1)= eri(l,1,j)
          work(2)= eri(l,2,j)
          work(3)= eri(l,3,j)
          eri(l,1,j)= work(1)*rot(1,1)+work(2)*rot(2,1)+work(3)*rot(3,1)
          eri(l,2,j)= work(1)*rot(1,2)+work(2)*rot(2,2)+work(3)*rot(3,2)
          eri(l,3,j)= work(1)*rot(1,3)+work(2)*rot(2,3)+work(3)*rot(3,3)
        enddo
      enddo
!
      if(nbfijkl(4) == 6)then
        do j= 1,nbfijkl(2)
          do k= 1,3
            do l= 1,6
              phmdint(l,k,j,1)= eri(1,k,j)*rot2(1,l)+eri(2,k,j)*rot2(2,l)+eri(3,k,j)*rot2(3,l) &
&                              +eri(4,k,j)*rot2(4,l)+eri(5,k,j)*rot2(5,l)+eri(6,k,j)*rot2(6,l) 
            enddo
          enddo
        enddo
      else
        do j= 1,nbfijkl(2)
          do k= 1,3
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
