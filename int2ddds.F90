!--------------------------------------------------------------------------
  subroutine int2ddds(phmdint,exfac1,exfac2,xyziq,xzkl,rot,nijkl,nbfijkl)
!--------------------------------------------------------------------------
!
! Calculate (dd|ds) integrals
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
      real(8),parameter :: p11=1.1D+01, p12=1.2D+01, p15=1.5D+01, p45=4.5D+01
      real(8),parameter :: sqrt3=1.73205080756888D+00, sqrt3h=8.660254037844386D-01
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*), xzkl(2), rot(3,3)
      real(8),intent(out) :: phmdint(6,6,6,6)
      real(8) :: tval, tval1, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:6)
      real(8) :: f0(2), f1(2,4), f2(3,4), f3(4,4), f4(5,4), f5(6,2), f6(7), ftw(6,4)
      real(8) :: r0(10), r1(3,13), r2(6,17), r3(10,12), r4(15,8), r5(21,3), r6(28)
      real(8) :: ex12, ex34, ex43, ex41, expq, expq2, expq4, ex3q, ex4q, c12, c34, zip
      real(8) :: xiq, yiq, ziq, xiq2, yiq2, xyiq, xiq4, yiq4, xyiq2, xiq6, yiq6, x4y2, x2y4
      real(8) :: xypq2, zpq, zpq2, zpq3, zpq4, zpq5, fac, ex33q, ex34q, ex44q, zjp
      real(8) :: pmd, qmd, qmd2, qmd3, qmd4, qmd4x, qmd4y, qmd4xy
      real(8) :: qx, qz, xx, xz, zz, xxx, xxz, xzz, zzz, xxxx, xxxz, xxzz, xzzz, zzzz
      real(8) :: eri(6,6,6), work(15), f1w(3,3), f2w(6,4), f3w(10,4), f4w(15,4), f5w(21,2)
      real(8) :: rot2(6,6), rot3(6,5)
!
! Zero-clear
!
      r0(1:10)     = zero
      r1(1:3 ,1:13)= zero
      r2(1:6 ,1:17)= zero
      r3(1:10,1:12)= zero
      r4(1:15,1:8) = zero
      r5(1:21,1:3) = zero
      r6(1:28)     = zero
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
        f0(1:2)    = zero
        f1(1:2,1:4)= zero
        f2(1:3,1:4)= zero
        f3(1:4,1:4)= zero
        f4(1:5,1:4)= zero
        f5(1:6,1:2)= zero
        f6(1:7)    = zero
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
            ft(6)= ft(5)*tinv*expq*p11
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
            do ii= 0,6
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
            ft(6)= ft(6)*fac*expq4*expq2
          endif
          zpq3= zpq2*zpq
          zpq4= zpq2*zpq2
          zpq5= zpq2*zpq2*zpq
          ftw(1,1)= zjp*zjp
          ftw(1,2)= pmd
          ftw(1,3)= pmd*zjp
          ftw(1,4)= pmd*pmd
          do i= 1,4
            ftw(2,i)= ftw(1,i)*zpq
            ftw(3,i)= ftw(1,i)*zpq2
            ftw(4,i)= ftw(1,i)*zpq3
            ftw(5,i)= ftw(1,i)*zpq4
            ftw(6,i)= ftw(1,i)*zpq5
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
          do i= 1,4
            do j= 1,5
              f4(j,i)= f4(j,i)+ft(4)*ftw(j,i)
            enddo
          enddo
          do i= 1,2
            do j= 1,6
              f5(j,i)= f5(j,i)-ft(5)*ftw(j,i+2)
            enddo
          enddo
          f6(1)= f6(1)+ft(6)*ftw(1,4)
          f6(2)= f6(2)+ft(6)*ftw(2,4)
          f6(3)= f6(3)+ft(6)*ftw(3,4)
          f6(4)= f6(4)+ft(6)*ftw(4,4)
          f6(5)= f6(5)+ft(6)*ftw(5,4)
          f6(6)= f6(6)+ft(6)*ftw(6,4)
          f6(7)= f6(7)+ft(6)*ftw(6,4)*zpq
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
          r0(i  )= r0(i  )+f0(1)*work(i)
          r0(i+5)= r0(i+5)+f0(2)*work(i)
        enddo
!
        do i= 1,3
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
        do i= 9,13
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,3)*work(i-8)
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
            r2(j,i)= r2(j,i)+f2w(j,3)*work(i-3)
          enddo
        enddo
        do i= 13,17
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,4)*work(i-12)
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
        do i= 5,8
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,3)*work(i+5)
          enddo
        enddo
        do i= 9,12
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,4)*work(i-3)
          enddo
        enddo
!
        do i= 1,4
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
          r4(j,3)= r4(j,3)+f4w(j,3)*work(14)
          r4(j,4)= r4(j,4)+f4w(j,3)*work(15)
        enddo
        do i= 5,8
          do j= 1,15
            r4(j,i)= r4(j,i)+f4w(j,4)*work(i+5)
          enddo
        enddo
!
        do i= 1,2
          f5w( 1,i)=(f5(1,i)*xiq4 +f4(1,i+2)*xiq2*ten                 +f3(1,i+2)*p15  )*xiq
          f5w( 2,i)=(f5(1,i)*xiq4 +f4(1,i+2)*xiq2*six                 +f3(1,i+2)*three)*yiq
          f5w( 3,i)=(f5(2,i)*xiq4 +f4(2,i+2)*xiq2*six                 +f3(2,i+2)*three)
          f5w( 4,i)=(f5(1,i)*xyiq2+f4(1,i+2)*xiq2+f4(1,i+2)*yiq2*three+f3(1,i+2)*three)*xiq
          f5w( 5,i)=(f5(2,i)*xiq2 +f4(2,i+2)*three                                    )*xyiq
          f5w( 6,i)=(f5(3,i)*xiq2 +f4(1,i+2)*xiq2+f4(3,i+2)*three     +f3(1,i+2)*three)*xiq
          f5w( 7,i)=(f5(1,i)*xyiq2+f4(1,i+2)*xiq2*three+f4(1,i+2)*yiq2+f3(1,i+2)*three)*yiq
          f5w( 8,i)=(f5(2,i)*xyiq2+f4(2,i+2)*xiq2+f4(2,i+2)*yiq2      +f3(2,i+2)      )
          f5w( 9,i)=(f5(3,i)*xiq2 +f4(1,i+2)*xiq2+f4(3,i+2)           +f3(1,i+2)      )*yiq
          f5w(10,i)=(f5(4,i)*xiq2 +f4(2,i+2)*xiq2*three+f4(4,i+2)     +f3(2,i+2)*three)
          f5w(11,i)=(f5(1,i)*yiq4 +f4(1,i+2)*yiq2*six                 +f3(1,i+2)*three)*xiq
          f5w(12,i)=(f5(2,i)*yiq2 +f4(2,i+2)*three                                    )*xyiq
          f5w(13,i)=(f5(3,i)*yiq2 +f4(1,i+2)*yiq2+f4(3,i+2)           +f3(1,i+2)      )*xiq
          f5w(14,i)=(f5(4,i)      +f4(2,i+2)*three                                    )*xyiq
          f5w(15,i)=(f5(5,i)      +f4(3,i+2)*six                      +f3(1,i+2)*three)*xiq
          f5w(16,i)=(f5(1,i)*yiq4 +f4(1,i+2)*yiq2*ten                 +f3(1,i+2)*p15  )*yiq
          f5w(17,i)=(f5(2,i)*yiq4 +f4(2,i+2)*yiq2*six                 +f3(2,i+2)*three)
          f5w(18,i)=(f5(3,i)*yiq2 +f4(1,i+2)*yiq2+f4(3,i+2)*three     +f3(1,i+2)*three)*yiq
          f5w(19,i)=(f5(4,i)*yiq2 +f4(2,i+2)*yiq2*three+f4(4,i+2)     +f3(2,i+2)*three)
          f5w(20,i)=(f5(5,i)      +f4(3,i+2)*six                      +f3(1,i+2)*three)*yiq
          f5w(21,i)=(f5(6,i)      +f4(4,i+2)*ten                      +f3(2,i+2)*p15  )
        enddo
        do j= 1,21
          r5(j,1)= r5(j,1)+f5w(j,1)*qmd4
          r5(j,2)= r5(j,2)+f5w(j,2)*work(14)
          r5(j,3)= r5(j,3)+f5w(j,2)*work(15)
        enddo
!
        r6( 1)= r6( 1)+(f6(1)*xiq6+f5(1,2)*xiq4*p15+f4(1,4)*xiq2*p45+f3(1,4)*p15)*qmd4
        r6( 2)= r6( 2)+(f6(1)*xiq4+f5(1,2)*xiq2*ten+f4(1,4)*p15)*qmd4xy
        r6( 3)= r6( 3)+(f6(2)*xiq4+f5(2,2)*xiq2*ten+f4(2,4)*p15)*qmd4x
        r6( 4)= r6( 4)+(f6(1)*x4y2+f5(1,2)*xiq4+f5(1,2)*xyiq2*six+f4(1,4)*xiq2*six &
&                      +f4(1,4)*yiq2*three+f3(1,4)*three)*qmd4
        r6( 5)= r6( 5)+(f6(2)*xiq4+f5(2,2)*xiq2*six+f4(2,4)*three)*qmd4y
        r6( 6)= r6( 6)+(f6(3)*xiq4+f5(1,2)*xiq4+f5(3,2)*xiq2*six+f4(1,4)*xiq2*six &
&                      +f4(3,4)*three+f3(1,4)*three)*qmd4
        r6( 7)= r6( 7)+(f6(1)*xyiq2+f5(1,2)*xiq2*three+f5(1,2)*yiq2*three+f4(1,4)*nine)*qmd4xy
        r6( 8)= r6( 8)+(f6(2)*xyiq2+f5(2,2)*xiq2+f5(2,2)*yiq2*three+f4(2,4)*three)*qmd4x
        r6( 9)= r6( 9)+(f6(3)*xiq2+f5(1,2)*xiq2+f5(3,2)*three+f4(1,4)*three)*qmd4xy
        r6(10)= r6(10)+(f6(4)*xiq2+f5(2,2)*xiq2*three+f5(4,2)*three+f4(2,4)*nine)*qmd4x
        r6(11)= r6(11)+(f6(1)*x2y4+f5(1,2)*xyiq2*six+f5(1,2)*yiq4+f4(1,4)*xiq2*three &
&                      +f4(1,4)*yiq2*six+f3(1,4)*three)*qmd4
        r6(12)= r6(12)+(f6(2)*xyiq2+f5(2,2)*xiq2*three+f5(2,2)*yiq2+f4(2,4)*three)*qmd4y
        r6(13)= r6(13)+(f6(3)*xyiq2+f5(1,2)*xyiq2+f5(3,2)*xiq2+f5(3,2)*yiq2+f4(1,4)*xiq2 &
&                      +f4(1,4)*yiq2+f4(3,4)+f3(1,4))*qmd4
        r6(14)= r6(14)+(f6(4)*xiq2+f5(2,2)*xiq2*three+f5(4,2)+f4(2,4)*three)*qmd4y
        r6(15)= r6(15)+(f6(5)*xiq2+f5(3,2)*xiq2*six+f5(5,2)+f4(1,4)*xiq2*three+f4(3,4)*six &
&                      +f3(1,4)*three)*qmd4
        r6(16)= r6(16)+(f6(1)*yiq4+f5(1,2)*yiq2*ten+f4(1,4)*p15)*qmd4xy
        r6(17)= r6(17)+(f6(2)*yiq4+f5(2,2)*yiq2*six+f4(2,4)*three)*qmd4x
        r6(18)= r6(18)+(f6(3)*yiq2+f5(1,2)*yiq2+f5(3,2)*three+f4(1,4)*three)*qmd4xy
        r6(19)= r6(19)+(f6(4)*yiq2+f5(2,2)*yiq2*three+f5(4,2)+f4(2,4)*three)*qmd4x
        r6(20)= r6(20)+(f6(5)+f5(3,2)*six+f4(1,4)*three)*qmd4xy
        r6(21)= r6(21)+(f6(6)+f5(4,2)*ten+f4(2,4)*p15)*qmd4x
        r6(22)= r6(22)+(f6(1)*yiq6+f5(1,2)*yiq4*p15+f4(1,4)*yiq2*p45+f3(1,4)*p15)*qmd4
        r6(23)= r6(23)+(f6(2)*yiq4+f5(2,2)*yiq2*ten+f4(2,4)*p15)*qmd4y
        r6(24)= r6(24)+(f6(3)*yiq4+f5(1,2)*yiq4+f5(3,2)*yiq2*six+f4(1,4)*yiq2*six &
&                      +f4(3,4)*three+f3(1,4)*three)*qmd4
        r6(25)= r6(25)+(f6(4)*yiq2+f5(2,2)*yiq2*three+f5(4,2)*three+f4(2,4)*nine)*qmd4y
        r6(26)= r6(26)+(f6(5)*yiq2+f5(3,2)*yiq2*six+f5(5,2)+f4(1,4)*yiq2*three+f4(3,4)*six &
&                      +f3(1,4)*three)*qmd4
        r6(27)= r6(27)+(f6(6)+f5(4,2)*ten+f4(2,4)*p15)*qmd4y
        r6(28)= r6(28)+(f6(7)+f5(5,2)*p15+f4(3,4)*p45+f3(1,4)*p15)*qmd4
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
      eri(1,1,1)= r6(1)+r4(1,2)+r4(1,5)*six+r2(1,5)*six+r2(1,13)*three+r0(6)*three+( &
&                +r5(1,2)*two+r5(1,3)*two+r3(1,3)*two+r3(1,4)*two+r3(1,9)*six+r3(1,10)*six &
&                +r1(1,5)*six+r1(1,6)*six)*qx+(+r4(1,6)+r4(1,7)*four+r4(1,8)+r2(1,6) &
&                +r2(1,7)*four+r2(1,8)+r2(1,14)+r2(1,15)*four+r2(1,16)+r0(7)+r0(8)*four+r0(9) &
&                )*xx+(+r3(1,11)*two+r3(1,12)*two+r1(1,7)*two+r1(1,8)*two)*xxx+(+r2(1,17) &
&                +r0(10))*xxxx
      eri(2,1,1)= r6(4)+r4(4,2)+r4(1,5)+r4(4,5)+r2(1,5)+r2(2,5)+r2(1,13)+r0(6)+( &
&                +r5(4,3)*two+r3(4,4)*two+r3(1,10)*two+r1(1,6)*two)*qx+(+r4(4,8)+r2(2,8) &
&                +r2(1,16)+r0(9))*xx
      eri(3,1,1)= r6(6)+r4(6,2)+r4(1,5)+r4(6,5)+r2(1,5)+r2(3,5)+r2(1,13)+r0(6)+( &
&                +r5(6,3)*two+r3(6,4)*two+r3(1,10)*two+r1(1,6)*two)*qx+(+r5(3,2)*two &
&                +r3(3,3)*two+r3(3,9)*two+r1(3,5)*two)*qz+(+r4(6,8)+r2(3,8)+r2(1,16)+r0(9)) &
&                *xx+(+r4(3,7)*four+r2(5,7)*four)*xz+(+r4(1,6)+r2(1,6)+r2(1,14)+r0(7))*zz+( &
&                +r3(3,12)*two+r1(3,8)*two)*xxz+(+r3(1,11)*two+r1(1,7)*two)*xzz+(+r2(1,17) &
&                +r0(10))*xxzz
      eri(4,1,1)= r6(2)+r4(2,2)+r4(2,5)*three+r2(4,5)*three+(+r5(2,2)+r5(2,3)*two+r3(2,3) &
&                +r3(2,4)*two+r3(2,9)+r3(2,10)*two+r1(2,5)+r1(2,6)*two)*qx+(+r4(2,7)*two &
&                +r4(2,8)+r2(4,7)*two+r2(4,8))*xx+(+r3(2,12)+r1(2,8))*xxx
      eri(5,1,1)= r6(3)+r4(3,2)+r4(3,5)*three+r2(5,5)*three+(+r5(3,2)+r5(3,3)*two+r3(3,3) &
&                +r3(3,4)*two+r3(3,9)+r3(3,10)*two+r1(3,5)+r1(3,6)*two)*qx+(+r5(1,2)+r3(1,3) &
&                +r3(1,9)*three+r1(1,5)*three)*qz+(+r4(3,7)*two+r4(3,8)+r2(5,7)*two+r2(5,8)) &
&                *xx+(+r4(1,6)+r4(1,7)*two+r2(1,6)+r2(1,7)*two+r2(1,14)+r2(1,15)*two+r0(7) &
&                +r0(8)*two)*xz+(+r3(3,12)+r1(3,8))*xxx+(+r3(1,11)*two+r3(1,12)+r1(1,7)*two &
&                +r1(1,8))*xxz+(+r2(1,17)+r0(10))*xxxz
      eri(6,1,1)= r6(5)+r4(5,2)+r4(5,5)+r2(6,5)+(+r5(5,3)*two+r3(5,4)*two)*qx+(+r5(2,2) &
&                +r3(2,3)+r3(2,9)+r1(2,5))*qz+(+r4(5,8)+r2(6,8))*xx+(+r4(2,7)*two+r2(4,7)*two &
&                )*xz+(+r3(2,12)+r1(2,8))*xxz
      eri(1,2,1)= r6(4)+r4(4,2)+r4(1,5)+r4(4,5)+r2(1,5)+r2(2,5)+r2(1,13)+r0(6)+( &
&                +r5(4,2)*two+r3(4,3)*two+r3(1,9)*two+r1(1,5)*two)*qx+(+r4(4,6)+r2(2,6) &
&                +r2(1,14)+r0(7))*xx
      eri(2,2,1)= r6(11)+r4(11,2)+r4(4,5)*six+r2(2,5)*six+r2(1,13)*three+r0(6)*three
      eri(3,2,1)= r6(13)+r4(13,2)+r4(4,5)+r4(6,5)+r2(2,5)+r2(3,5)+r2(1,13)+r0(6)+( &
&                +r5(8,2)*two+r3(8,3)*two+r3(3,9)*two+r1(3,5)*two)*qz+(+r4(4,6)+r2(2,6) &
&                +r2(1,14)+r0(7))*zz
      eri(4,2,1)= r6(7)+r4(7,2)+r4(2,5)*three+r2(4,5)*three+(+r5(7,2)+r3(7,3)+r3(2,9)*three &
&                +r1(2,5)*three)*qx
      eri(5,2,1)= r6(8)+r4(8,2)+r4(3,5)+r2(5,5)+(+r5(8,2)+r3(8,3)+r3(3,9)+r1(3,5))*qx+( &
&                +r5(4,2)+r3(4,3)+r3(1,9)+r1(1,5))*qz+(+r4(4,6)+r2(2,6)+r2(1,14)+r0(7))*xz
      eri(6,2,1)= r6(12)+r4(12,2)+r4(5,5)*three+r2(6,5)*three+(+r5(7,2)+r3(7,3) &
&                +r3(2,9)*three+r1(2,5)*three)*qz
      eri(1,3,1)= r6(6)+r4(6,2)+r4(1,5)+r4(6,5)+r2(1,5)+r2(3,5)+r2(1,13)+r0(6)+( &
&                +r5(6,2)*two+r3(6,3)*two+r3(1,9)*two+r1(1,5)*two)*qx+(+r5(3,3)*two &
&                +r3(3,4)*two+r3(3,10)*two+r1(3,6)*two)*qz+(+r4(6,6)+r2(3,6)+r2(1,14)+r0(7)) &
&                *xx+(+r4(3,7)*four+r2(5,7)*four)*xz+(+r4(1,8)+r2(1,8)+r2(1,16)+r0(9))*zz+( &
&                +r3(3,11)*two+r1(3,7)*two)*xxz+(+r3(1,12)*two+r1(1,8)*two)*xzz+(+r2(1,17) &
&                +r0(10))*xxzz
      eri(2,3,1)= r6(13)+r4(13,2)+r4(4,5)+r4(6,5)+r2(2,5)+r2(3,5)+r2(1,13)+r0(6)+( &
&                +r5(8,3)*two+r3(8,4)*two+r3(3,10)*two+r1(3,6)*two)*qz+(+r4(4,8)+r2(2,8) &
&                +r2(1,16)+r0(9))*zz
      eri(3,3,1)= r6(15)+r4(15,2)+r4(6,5)*six+r2(3,5)*six+r2(1,13)*three+r0(6)*three+( &
&                +r5(10,2)*two+r5(10,3)*two+r3(10,3)*two+r3(10,4)*two+r3(3,9)*six &
&                +r3(3,10)*six+r1(3,5)*six+r1(3,6)*six)*qz+(+r4(6,6)+r4(6,7)*four+r4(6,8) &
&                +r2(3,6)+r2(3,7)*four+r2(3,8)+r2(1,14)+r2(1,15)*four+r2(1,16)+r0(7) &
&                +r0(8)*four+r0(9))*zz+(+r3(3,11)*two+r3(3,12)*two+r1(3,7)*two+r1(3,8)*two) &
&                *zzz+(+r2(1,17)+r0(10))*zzzz
      eri(4,3,1)= r6(9)+r4(9,2)+r4(2,5)+r2(4,5)+(+r5(9,2)+r3(9,3)+r3(2,9)+r1(2,5))*qx+( &
&                +r5(5,3)*two+r3(5,4)*two)*qz+(+r4(5,7)*two+r2(6,7)*two)*xz+(+r4(2,8)+r2(4,8) &
&                )*zz+(+r3(2,12)+r1(2,8))*xzz
      eri(5,3,1)= r6(10)+r4(10,2)+r4(3,5)*three+r2(5,5)*three+(+r5(10,2)+r3(10,3) &
&                +r3(3,9)*three+r1(3,5)*three)*qx+(+r5(6,2)+r5(6,3)*two+r3(6,3)+r3(6,4)*two &
&                +r3(1,9)+r3(1,10)*two+r1(1,5)+r1(1,6)*two)*qz+(+r4(6,6)+r4(6,7)*two+r2(3,6) &
&                +r2(3,7)*two+r2(1,14)+r2(1,15)*two+r0(7)+r0(8)*two)*xz+(+r4(3,7)*two+r4(3,8) &
&                +r2(5,7)*two+r2(5,8))*zz+(+r3(3,11)*two+r3(3,12)+r1(3,7)*two+r1(3,8))*xzz+( &
&                +r3(1,12)+r1(1,8))*zzz+(+r2(1,17)+r0(10))*xzzz
      eri(6,3,1)= r6(14)+r4(14,2)+r4(5,5)*three+r2(6,5)*three+(+r5(9,2)+r5(9,3)*two+r3(9,3) &
&                +r3(9,4)*two+r3(2,9)+r3(2,10)*two+r1(2,5)+r1(2,6)*two)*qz+(+r4(5,7)*two &
&                +r4(5,8)+r2(6,7)*two+r2(6,8))*zz+(+r3(2,12)+r1(2,8))*zzz
      eri(1,4,1)= r6(2)+r4(2,2)+r4(2,5)*three+r2(4,5)*three+(+r5(2,2)*two+r5(2,3) &
&                +r3(2,3)*two+r3(2,4)+r3(2,9)*two+r3(2,10)+r1(2,5)*two+r1(2,6))*qx+(+r4(2,6) &
&                +r4(2,7)*two+r2(4,6)+r2(4,7)*two)*xx+(+r3(2,11)+r1(2,7))*xxx
      eri(2,4,1)= r6(7)+r4(7,2)+r4(2,5)*three+r2(4,5)*three+(+r5(7,3)+r3(7,4) &
&                +r3(2,10)*three+r1(2,6)*three)*qx
      eri(3,4,1)= r6(9)+r4(9,2)+r4(2,5)+r2(4,5)+(+r5(9,3)+r3(9,4)+r3(2,10)+r1(2,6))*qx+( &
&                +r5(5,2)*two+r3(5,3)*two)*qz+(+r4(5,7)*two+r2(6,7)*two)*xz+(+r4(2,6)+r2(4,6) &
&                )*zz+(+r3(2,11)+r1(2,7))*xzz
      eri(4,4,1)= r6(4)+r4(4,2)+r4(1,5)+r4(4,5)+r2(1,5)+r2(2,5)+r2(1,13)+r0(6)+(+r5(4,2) &
&                +r5(4,3)+r3(4,3)+r3(4,4)+r3(1,9)+r3(1,10)+r1(1,5)+r1(1,6))*qx+(+r4(4,7) &
&                +r2(2,7)+r2(1,15)+r0(8))*xx
      eri(5,4,1)= r6(5)+r4(5,2)+r4(5,5)+r2(6,5)+(+r5(5,2)+r5(5,3)+r3(5,3)+r3(5,4))*qx+( &
&                +r5(2,2)+r3(2,3)+r3(2,9)+r1(2,5))*qz+(+r4(5,7)+r2(6,7))*xx+(+r4(2,6)+r4(2,7) &
&                +r2(4,6)+r2(4,7))*xz+(+r3(2,11)+r1(2,7))*xxz
      eri(6,4,1)= r6(8)+r4(8,2)+r4(3,5)+r2(5,5)+(+r5(8,3)+r3(8,4)+r3(3,10)+r1(3,6))*qx+( &
&                +r5(4,2)+r3(4,3)+r3(1,9)+r1(1,5))*qz+(+r4(4,7)+r2(2,7)+r2(1,15)+r0(8))*xz
      eri(1,5,1)= r6(3)+r4(3,2)+r4(3,5)*three+r2(5,5)*three+(+r5(3,2)*two+r5(3,3) &
&                +r3(3,3)*two+r3(3,4)+r3(3,9)*two+r3(3,10)+r1(3,5)*two+r1(3,6))*qx+(+r5(1,3) &
&                +r3(1,4)+r3(1,10)*three+r1(1,6)*three)*qz+(+r4(3,6)+r4(3,7)*two+r2(5,6) &
&                +r2(5,7)*two)*xx+(+r4(1,7)*two+r4(1,8)+r2(1,7)*two+r2(1,8)+r2(1,15)*two &
&                +r2(1,16)+r0(8)*two+r0(9))*xz+(+r3(3,11)+r1(3,7))*xxx+(+r3(1,11) &
&                +r3(1,12)*two+r1(1,7)+r1(1,8)*two)*xxz+(+r2(1,17)+r0(10))*xxxz
      eri(2,5,1)= r6(8)+r4(8,2)+r4(3,5)+r2(5,5)+(+r5(8,3)+r3(8,4)+r3(3,10)+r1(3,6))*qx+( &
&                +r5(4,3)+r3(4,4)+r3(1,10)+r1(1,6))*qz+(+r4(4,8)+r2(2,8)+r2(1,16)+r0(9))*xz
      eri(3,5,1)= r6(10)+r4(10,2)+r4(3,5)*three+r2(5,5)*three+(+r5(10,3)+r3(10,4) &
&                +r3(3,10)*three+r1(3,6)*three)*qx+(+r5(6,2)*two+r5(6,3)+r3(6,3)*two+r3(6,4) &
&                +r3(1,9)*two+r3(1,10)+r1(1,5)*two+r1(1,6))*qz+(+r4(6,7)*two+r4(6,8) &
&                +r2(3,7)*two+r2(3,8)+r2(1,15)*two+r2(1,16)+r0(8)*two+r0(9))*xz+(+r4(3,6) &
&                +r4(3,7)*two+r2(5,6)+r2(5,7)*two)*zz+(+r3(3,11)+r3(3,12)*two+r1(3,7) &
&                +r1(3,8)*two)*xzz+(+r3(1,11)+r1(1,7))*zzz+(+r2(1,17)+r0(10))*xzzz
      eri(4,5,1)= r6(5)+r4(5,2)+r4(5,5)+r2(6,5)+(+r5(5,2)+r5(5,3)+r3(5,3)+r3(5,4))*qx+( &
&                +r5(2,3)+r3(2,4)+r3(2,10)+r1(2,6))*qz+(+r4(5,7)+r2(6,7))*xx+(+r4(2,7) &
&                +r4(2,8)+r2(4,7)+r2(4,8))*xz+(+r3(2,12)+r1(2,8))*xxz
      eri(5,5,1)= r6(6)+r4(6,2)+r4(1,5)+r4(6,5)+r2(1,5)+r2(3,5)+r2(1,13)+r0(6)+(+r5(6,2) &
&                +r5(6,3)+r3(6,3)+r3(6,4)+r3(1,9)+r3(1,10)+r1(1,5)+r1(1,6))*qx+(+r5(3,2) &
&                +r5(3,3)+r3(3,3)+r3(3,4)+r3(3,9)+r3(3,10)+r1(3,5)+r1(3,6))*qz+(+r4(6,7) &
&                +r2(3,7)+r2(1,15)+r0(8))*xx+(+r4(3,6)+r4(3,7)*two+r4(3,8)+r2(5,6) &
&                +r2(5,7)*two+r2(5,8))*xz+(+r4(1,7)+r2(1,7)+r2(1,15)+r0(8))*zz+(+r3(3,11) &
&                +r3(3,12)+r1(3,7)+r1(3,8))*xxz+(+r3(1,11)+r3(1,12)+r1(1,7)+r1(1,8))*xzz+( &
&                +r2(1,17)+r0(10))*xxzz
      eri(6,5,1)= r6(9)+r4(9,2)+r4(2,5)+r2(4,5)+(+r5(9,3)+r3(9,4)+r3(2,10)+r1(2,6))*qx+( &
&                +r5(5,2)+r5(5,3)+r3(5,3)+r3(5,4))*qz+(+r4(5,7)+r4(5,8)+r2(6,7)+r2(6,8))*xz+( &
&                +r4(2,7)+r2(4,7))*zz+(+r3(2,12)+r1(2,8))*xzz
      eri(1,6,1)= r6(5)+r4(5,2)+r4(5,5)+r2(6,5)+(+r5(5,2)*two+r3(5,3)*two)*qx+(+r5(2,3) &
&                +r3(2,4)+r3(2,10)+r1(2,6))*qz+(+r4(5,6)+r2(6,6))*xx+(+r4(2,7)*two &
&                +r2(4,7)*two)*xz+(+r3(2,11)+r1(2,7))*xxz
      eri(2,6,1)= r6(12)+r4(12,2)+r4(5,5)*three+r2(6,5)*three+(+r5(7,3)+r3(7,4) &
&                +r3(2,10)*three+r1(2,6)*three)*qz
      eri(3,6,1)= r6(14)+r4(14,2)+r4(5,5)*three+r2(6,5)*three+(+r5(9,2)*two+r5(9,3) &
&                +r3(9,3)*two+r3(9,4)+r3(2,9)*two+r3(2,10)+r1(2,5)*two+r1(2,6))*qz+(+r4(5,6) &
&                +r4(5,7)*two+r2(6,6)+r2(6,7)*two)*zz+(+r3(2,11)+r1(2,7))*zzz
      eri(4,6,1)= r6(8)+r4(8,2)+r4(3,5)+r2(5,5)+(+r5(8,2)+r3(8,3)+r3(3,9)+r1(3,5))*qx+( &
&                +r5(4,3)+r3(4,4)+r3(1,10)+r1(1,6))*qz+(+r4(4,7)+r2(2,7)+r2(1,15)+r0(8))*xz
      eri(5,6,1)= r6(9)+r4(9,2)+r4(2,5)+r2(4,5)+(+r5(9,2)+r3(9,3)+r3(2,9)+r1(2,5))*qx+( &
&                +r5(5,2)+r5(5,3)+r3(5,3)+r3(5,4))*qz+(+r4(5,6)+r4(5,7)+r2(6,6)+r2(6,7))*xz+( &
&                +r4(2,7)+r2(4,7))*zz+(+r3(2,11)+r1(2,7))*xzz
      eri(6,6,1)= r6(13)+r4(13,2)+r4(4,5)+r4(6,5)+r2(2,5)+r2(3,5)+r2(1,13)+r0(6)+(+r5(8,2) &
&                +r5(8,3)+r3(8,3)+r3(8,4)+r3(3,9)+r3(3,10)+r1(3,5)+r1(3,6))*qz+(+r4(4,7) &
&                +r2(2,7)+r2(1,15)+r0(8))*zz
      eri(1,1,2)= r6(4)+r4(1,2)+r4(4,5)*six+r2(1,5)*six+r2(2,13)*three+r0(6)*three+( &
&                +r5(4,2)*two+r5(4,3)*two+r3(1,3)*two+r3(1,4)*two+r3(4,9)*six+r3(4,10)*six &
&                +r1(1,5)*six+r1(1,6)*six)*qx+(+r4(4,6)+r4(4,7)*four+r4(4,8)+r2(1,6) &
&                +r2(1,7)*four+r2(1,8)+r2(2,14)+r2(2,15)*four+r2(2,16)+r0(7)+r0(8)*four+r0(9) &
&                )*xx+(+r3(4,11)*two+r3(4,12)*two+r1(1,7)*two+r1(1,8)*two)*xxx+(+r2(2,17) &
&                +r0(10))*xxxx
      eri(2,1,2)= r6(11)+r4(4,2)+r4(4,5)+r4(11,5)+r2(1,5)+r2(2,5)+r2(2,13)+r0(6)+( &
&                +r5(11,3)*two+r3(4,4)*two+r3(4,10)*two+r1(1,6)*two)*qx+(+r4(11,8)+r2(2,8) &
&                +r2(2,16)+r0(9))*xx
      eri(3,1,2)= r6(13)+r4(6,2)+r4(4,5)+r4(13,5)+r2(1,5)+r2(3,5)+r2(2,13)+r0(6)+( &
&                +r5(13,3)*two+r3(6,4)*two+r3(4,10)*two+r1(1,6)*two)*qx+(+r5(8,2)*two &
&                +r3(3,3)*two+r3(8,9)*two+r1(3,5)*two)*qz+(+r4(13,8)+r2(3,8)+r2(2,16)+r0(9)) &
&                *xx+(+r4(8,7)*four+r2(5,7)*four)*xz+(+r4(4,6)+r2(1,6)+r2(2,14)+r0(7))*zz+( &
&                +r3(8,12)*two+r1(3,8)*two)*xxz+(+r3(4,11)*two+r1(1,7)*two)*xzz+(+r2(2,17) &
&                +r0(10))*xxzz
      eri(4,1,2)= r6(7)+r4(2,2)+r4(7,5)*three+r2(4,5)*three+(+r5(7,2)+r5(7,3)*two+r3(2,3) &
&                +r3(2,4)*two+r3(7,9)+r3(7,10)*two+r1(2,5)+r1(2,6)*two)*qx+(+r4(7,7)*two &
&                +r4(7,8)+r2(4,7)*two+r2(4,8))*xx+(+r3(7,12)+r1(2,8))*xxx
      eri(5,1,2)= r6(8)+r4(3,2)+r4(8,5)*three+r2(5,5)*three+(+r5(8,2)+r5(8,3)*two+r3(3,3) &
&                +r3(3,4)*two+r3(8,9)+r3(8,10)*two+r1(3,5)+r1(3,6)*two)*qx+(+r5(4,2)+r3(1,3) &
&                +r3(4,9)*three+r1(1,5)*three)*qz+(+r4(8,7)*two+r4(8,8)+r2(5,7)*two+r2(5,8)) &
&                *xx+(+r4(4,6)+r4(4,7)*two+r2(1,6)+r2(1,7)*two+r2(2,14)+r2(2,15)*two+r0(7) &
&                +r0(8)*two)*xz+(+r3(8,12)+r1(3,8))*xxx+(+r3(4,11)*two+r3(4,12)+r1(1,7)*two &
&                +r1(1,8))*xxz+(+r2(2,17)+r0(10))*xxxz
      eri(6,1,2)= r6(12)+r4(5,2)+r4(12,5)+r2(6,5)+(+r5(12,3)*two+r3(5,4)*two)*qx+(+r5(7,2) &
&                +r3(2,3)+r3(7,9)+r1(2,5))*qz+(+r4(12,8)+r2(6,8))*xx+(+r4(7,7)*two &
&                +r2(4,7)*two)*xz+(+r3(7,12)+r1(2,8))*xxz
      eri(1,2,2)= r6(11)+r4(4,2)+r4(4,5)+r4(11,5)+r2(1,5)+r2(2,5)+r2(2,13)+r0(6)+( &
&                +r5(11,2)*two+r3(4,3)*two+r3(4,9)*two+r1(1,5)*two)*qx+(+r4(11,6)+r2(2,6) &
&                +r2(2,14)+r0(7))*xx
      eri(2,2,2)= r6(22)+r4(11,2)+r4(11,5)*six+r2(2,5)*six+r2(2,13)*three+r0(6)*three
      eri(3,2,2)= r6(24)+r4(13,2)+r4(11,5)+r4(13,5)+r2(2,5)+r2(3,5)+r2(2,13)+r0(6)+( &
&                +r5(17,2)*two+r3(8,3)*two+r3(8,9)*two+r1(3,5)*two)*qz+(+r4(11,6)+r2(2,6) &
&                +r2(2,14)+r0(7))*zz
      eri(4,2,2)= r6(16)+r4(7,2)+r4(7,5)*three+r2(4,5)*three+(+r5(16,2)+r3(7,3) &
&                +r3(7,9)*three+r1(2,5)*three)*qx
      eri(5,2,2)= r6(17)+r4(8,2)+r4(8,5)+r2(5,5)+(+r5(17,2)+r3(8,3)+r3(8,9)+r1(3,5))*qx+( &
&                +r5(11,2)+r3(4,3)+r3(4,9)+r1(1,5))*qz+(+r4(11,6)+r2(2,6)+r2(2,14)+r0(7))*xz
      eri(6,2,2)= r6(23)+r4(12,2)+r4(12,5)*three+r2(6,5)*three+(+r5(16,2)+r3(7,3) &
&                +r3(7,9)*three+r1(2,5)*three)*qz
      eri(1,3,2)= r6(13)+r4(6,2)+r4(4,5)+r4(13,5)+r2(1,5)+r2(3,5)+r2(2,13)+r0(6)+( &
&                +r5(13,2)*two+r3(6,3)*two+r3(4,9)*two+r1(1,5)*two)*qx+(+r5(8,3)*two &
&                +r3(3,4)*two+r3(8,10)*two+r1(3,6)*two)*qz+(+r4(13,6)+r2(3,6)+r2(2,14)+r0(7)) &
&                *xx+(+r4(8,7)*four+r2(5,7)*four)*xz+(+r4(4,8)+r2(1,8)+r2(2,16)+r0(9))*zz+( &
&                +r3(8,11)*two+r1(3,7)*two)*xxz+(+r3(4,12)*two+r1(1,8)*two)*xzz+(+r2(2,17) &
&                +r0(10))*xxzz
      eri(2,3,2)= r6(24)+r4(13,2)+r4(11,5)+r4(13,5)+r2(2,5)+r2(3,5)+r2(2,13)+r0(6)+( &
&                +r5(17,3)*two+r3(8,4)*two+r3(8,10)*two+r1(3,6)*two)*qz+(+r4(11,8)+r2(2,8) &
&                +r2(2,16)+r0(9))*zz
      eri(3,3,2)= r6(26)+r4(15,2)+r4(13,5)*six+r2(3,5)*six+r2(2,13)*three+r0(6)*three+( &
&                +r5(19,2)*two+r5(19,3)*two+r3(10,3)*two+r3(10,4)*two+r3(8,9)*six &
&                +r3(8,10)*six+r1(3,5)*six+r1(3,6)*six)*qz+(+r4(13,6)+r4(13,7)*four+r4(13,8) &
&                +r2(3,6)+r2(3,7)*four+r2(3,8)+r2(2,14)+r2(2,15)*four+r2(2,16)+r0(7) &
&                +r0(8)*four+r0(9))*zz+(+r3(8,11)*two+r3(8,12)*two+r1(3,7)*two+r1(3,8)*two) &
&                *zzz+(+r2(2,17)+r0(10))*zzzz
      eri(4,3,2)= r6(18)+r4(9,2)+r4(7,5)+r2(4,5)+(+r5(18,2)+r3(9,3)+r3(7,9)+r1(2,5))*qx+( &
&                +r5(12,3)*two+r3(5,4)*two)*qz+(+r4(12,7)*two+r2(6,7)*two)*xz+(+r4(7,8) &
&                +r2(4,8))*zz+(+r3(7,12)+r1(2,8))*xzz
      eri(5,3,2)= r6(19)+r4(10,2)+r4(8,5)*three+r2(5,5)*three+(+r5(19,2)+r3(10,3) &
&                +r3(8,9)*three+r1(3,5)*three)*qx+(+r5(13,2)+r5(13,3)*two+r3(6,3)+r3(6,4)*two &
&                +r3(4,9)+r3(4,10)*two+r1(1,5)+r1(1,6)*two)*qz+(+r4(13,6)+r4(13,7)*two &
&                +r2(3,6)+r2(3,7)*two+r2(2,14)+r2(2,15)*two+r0(7)+r0(8)*two)*xz+(+r4(8,7)*two &
&                +r4(8,8)+r2(5,7)*two+r2(5,8))*zz+(+r3(8,11)*two+r3(8,12)+r1(3,7)*two+r1(3,8) &
&                )*xzz+(+r3(4,12)+r1(1,8))*zzz+(+r2(2,17)+r0(10))*xzzz
      eri(6,3,2)= r6(25)+r4(14,2)+r4(12,5)*three+r2(6,5)*three+(+r5(18,2)+r5(18,3)*two &
&                +r3(9,3)+r3(9,4)*two+r3(7,9)+r3(7,10)*two+r1(2,5)+r1(2,6)*two)*qz+( &
&                +r4(12,7)*two+r4(12,8)+r2(6,7)*two+r2(6,8))*zz+(+r3(7,12)+r1(2,8))*zzz
      eri(1,4,2)= r6(7)+r4(2,2)+r4(7,5)*three+r2(4,5)*three+(+r5(7,2)*two+r5(7,3) &
&                +r3(2,3)*two+r3(2,4)+r3(7,9)*two+r3(7,10)+r1(2,5)*two+r1(2,6))*qx+(+r4(7,6) &
&                +r4(7,7)*two+r2(4,6)+r2(4,7)*two)*xx+(+r3(7,11)+r1(2,7))*xxx
      eri(2,4,2)= r6(16)+r4(7,2)+r4(7,5)*three+r2(4,5)*three+(+r5(16,3)+r3(7,4) &
&                +r3(7,10)*three+r1(2,6)*three)*qx
      eri(3,4,2)= r6(18)+r4(9,2)+r4(7,5)+r2(4,5)+(+r5(18,3)+r3(9,4)+r3(7,10)+r1(2,6))*qx+( &
&                +r5(12,2)*two+r3(5,3)*two)*qz+(+r4(12,7)*two+r2(6,7)*two)*xz+(+r4(7,6) &
&                +r2(4,6))*zz+(+r3(7,11)+r1(2,7))*xzz
      eri(4,4,2)= r6(11)+r4(4,2)+r4(4,5)+r4(11,5)+r2(1,5)+r2(2,5)+r2(2,13)+r0(6)+(+r5(11,2) &
&                +r5(11,3)+r3(4,3)+r3(4,4)+r3(4,9)+r3(4,10)+r1(1,5)+r1(1,6))*qx+(+r4(11,7) &
&                +r2(2,7)+r2(2,15)+r0(8))*xx
      eri(5,4,2)= r6(12)+r4(5,2)+r4(12,5)+r2(6,5)+(+r5(12,2)+r5(12,3)+r3(5,3)+r3(5,4))*qx+( &
&                +r5(7,2)+r3(2,3)+r3(7,9)+r1(2,5))*qz+(+r4(12,7)+r2(6,7))*xx+(+r4(7,6) &
&                +r4(7,7)+r2(4,6)+r2(4,7))*xz+(+r3(7,11)+r1(2,7))*xxz
      eri(6,4,2)= r6(17)+r4(8,2)+r4(8,5)+r2(5,5)+(+r5(17,3)+r3(8,4)+r3(8,10)+r1(3,6))*qx+( &
&                +r5(11,2)+r3(4,3)+r3(4,9)+r1(1,5))*qz+(+r4(11,7)+r2(2,7)+r2(2,15)+r0(8))*xz
      eri(1,5,2)= r6(8)+r4(3,2)+r4(8,5)*three+r2(5,5)*three+(+r5(8,2)*two+r5(8,3) &
&                +r3(3,3)*two+r3(3,4)+r3(8,9)*two+r3(8,10)+r1(3,5)*two+r1(3,6))*qx+(+r5(4,3) &
&                +r3(1,4)+r3(4,10)*three+r1(1,6)*three)*qz+(+r4(8,6)+r4(8,7)*two+r2(5,6) &
&                +r2(5,7)*two)*xx+(+r4(4,7)*two+r4(4,8)+r2(1,7)*two+r2(1,8)+r2(2,15)*two &
&                +r2(2,16)+r0(8)*two+r0(9))*xz+(+r3(8,11)+r1(3,7))*xxx+(+r3(4,11) &
&                +r3(4,12)*two+r1(1,7)+r1(1,8)*two)*xxz+(+r2(2,17)+r0(10))*xxxz
      eri(2,5,2)= r6(17)+r4(8,2)+r4(8,5)+r2(5,5)+(+r5(17,3)+r3(8,4)+r3(8,10)+r1(3,6))*qx+( &
&                +r5(11,3)+r3(4,4)+r3(4,10)+r1(1,6))*qz+(+r4(11,8)+r2(2,8)+r2(2,16)+r0(9))*xz
      eri(3,5,2)= r6(19)+r4(10,2)+r4(8,5)*three+r2(5,5)*three+(+r5(19,3)+r3(10,4) &
&                +r3(8,10)*three+r1(3,6)*three)*qx+(+r5(13,2)*two+r5(13,3)+r3(6,3)*two &
&                +r3(6,4)+r3(4,9)*two+r3(4,10)+r1(1,5)*two+r1(1,6))*qz+(+r4(13,7)*two &
&                +r4(13,8)+r2(3,7)*two+r2(3,8)+r2(2,15)*two+r2(2,16)+r0(8)*two+r0(9))*xz+( &
&                +r4(8,6)+r4(8,7)*two+r2(5,6)+r2(5,7)*two)*zz+(+r3(8,11)+r3(8,12)*two+r1(3,7) &
&                +r1(3,8)*two)*xzz+(+r3(4,11)+r1(1,7))*zzz+(+r2(2,17)+r0(10))*xzzz
      eri(4,5,2)= r6(12)+r4(5,2)+r4(12,5)+r2(6,5)+(+r5(12,2)+r5(12,3)+r3(5,3)+r3(5,4))*qx+( &
&                +r5(7,3)+r3(2,4)+r3(7,10)+r1(2,6))*qz+(+r4(12,7)+r2(6,7))*xx+(+r4(7,7) &
&                +r4(7,8)+r2(4,7)+r2(4,8))*xz+(+r3(7,12)+r1(2,8))*xxz
      eri(5,5,2)= r6(13)+r4(6,2)+r4(4,5)+r4(13,5)+r2(1,5)+r2(3,5)+r2(2,13)+r0(6)+(+r5(13,2) &
&                +r5(13,3)+r3(6,3)+r3(6,4)+r3(4,9)+r3(4,10)+r1(1,5)+r1(1,6))*qx+(+r5(8,2) &
&                +r5(8,3)+r3(3,3)+r3(3,4)+r3(8,9)+r3(8,10)+r1(3,5)+r1(3,6))*qz+(+r4(13,7) &
&                +r2(3,7)+r2(2,15)+r0(8))*xx+(+r4(8,6)+r4(8,7)*two+r4(8,8)+r2(5,6) &
&                +r2(5,7)*two+r2(5,8))*xz+(+r4(4,7)+r2(1,7)+r2(2,15)+r0(8))*zz+(+r3(8,11) &
&                +r3(8,12)+r1(3,7)+r1(3,8))*xxz+(+r3(4,11)+r3(4,12)+r1(1,7)+r1(1,8))*xzz+( &
&                +r2(2,17)+r0(10))*xxzz
      eri(6,5,2)= r6(18)+r4(9,2)+r4(7,5)+r2(4,5)+(+r5(18,3)+r3(9,4)+r3(7,10)+r1(2,6))*qx+( &
&                +r5(12,2)+r5(12,3)+r3(5,3)+r3(5,4))*qz+(+r4(12,7)+r4(12,8)+r2(6,7)+r2(6,8)) &
&                *xz+(+r4(7,7)+r2(4,7))*zz+(+r3(7,12)+r1(2,8))*xzz
      eri(1,6,2)= r6(12)+r4(5,2)+r4(12,5)+r2(6,5)+(+r5(12,2)*two+r3(5,3)*two)*qx+(+r5(7,3) &
&                +r3(2,4)+r3(7,10)+r1(2,6))*qz+(+r4(12,6)+r2(6,6))*xx+(+r4(7,7)*two &
&                +r2(4,7)*two)*xz+(+r3(7,11)+r1(2,7))*xxz
      eri(2,6,2)= r6(23)+r4(12,2)+r4(12,5)*three+r2(6,5)*three+(+r5(16,3)+r3(7,4) &
&                +r3(7,10)*three+r1(2,6)*three)*qz
      eri(3,6,2)= r6(25)+r4(14,2)+r4(12,5)*three+r2(6,5)*three+(+r5(18,2)*two+r5(18,3) &
&                +r3(9,3)*two+r3(9,4)+r3(7,9)*two+r3(7,10)+r1(2,5)*two+r1(2,6))*qz+(+r4(12,6) &
&                +r4(12,7)*two+r2(6,6)+r2(6,7)*two)*zz+(+r3(7,11)+r1(2,7))*zzz
      eri(4,6,2)= r6(17)+r4(8,2)+r4(8,5)+r2(5,5)+(+r5(17,2)+r3(8,3)+r3(8,9)+r1(3,5))*qx+( &
&                +r5(11,3)+r3(4,4)+r3(4,10)+r1(1,6))*qz+(+r4(11,7)+r2(2,7)+r2(2,15)+r0(8))*xz
      eri(5,6,2)= r6(18)+r4(9,2)+r4(7,5)+r2(4,5)+(+r5(18,2)+r3(9,3)+r3(7,9)+r1(2,5))*qx+( &
&                +r5(12,2)+r5(12,3)+r3(5,3)+r3(5,4))*qz+(+r4(12,6)+r4(12,7)+r2(6,6)+r2(6,7)) &
&                *xz+(+r4(7,7)+r2(4,7))*zz+(+r3(7,11)+r1(2,7))*xzz
      eri(6,6,2)= r6(24)+r4(13,2)+r4(11,5)+r4(13,5)+r2(2,5)+r2(3,5)+r2(2,13)+r0(6)+( &
&                +r5(17,2)+r5(17,3)+r3(8,3)+r3(8,4)+r3(8,9)+r3(8,10)+r1(3,5)+r1(3,6))*qz+( &
&                +r4(11,7)+r2(2,7)+r2(2,15)+r0(8))*zz
      eri(1,1,3)= r6(6)-r5(3,1)*two+r4(1,1)+r4(1,2)+r4(6,5)*six-r3(3,5)*p12+r2(1,1)*six &
&                +r2(1,5)*six+r2(3,13)*three-r1(3,9)*six+r0(1)*three+r0(6)*three+( &
&                +r5(6,2)*two+r5(6,3)*two-r4(3,3)*four-r4(3,4)*four+r3(1,1)*two+r3(1,2)*two &
&                +r3(1,3)*two+r3(1,4)*two+r3(6,9)*six+r3(6,10)*six-r2(5,9)*p12-r2(5,10)*p12 &
&                +r1(1,1)*six+r1(1,2)*six+r1(1,5)*six+r1(1,6)*six)*qx+(+r4(6,6)+r4(6,7)*four &
&                +r4(6,8)-r3(3,6)*two-r3(3,7)*eight-r3(3,8)*two+r2(1,2)+r2(1,3)*four+r2(1,4) &
&                +r2(1,6)+r2(1,7)*four+r2(1,8)+r2(3,14)+r2(3,15)*four+r2(3,16)-r1(3,10)*two &
&                -r1(3,11)*eight-r1(3,12)*two+r0(2)+r0(3)*four+r0(4)+r0(7)+r0(8)*four+r0(9)) &
&                *xx+(+r3(6,11)*two+r3(6,12)*two-r2(5,11)*four-r2(5,12)*four+r1(1,3)*two &
&                +r1(1,4)*two+r1(1,7)*two+r1(1,8)*two)*xxx+(+r2(3,17)-r1(3,13)*two+r0(5) &
&                +r0(10))*xxxx
      eri(2,1,3)= r6(13)-r5(8,1)*two+r4(4,1)+r4(4,2)+r4(6,5)+r4(13,5)-r3(3,5)*two &
&                -r3(8,5)*two+r2(1,1)+r2(2,1)+r2(1,5)+r2(2,5)+r2(3,13)-r1(3,9)*two+r0(1) &
&                +r0(6)+(+r5(13,3)*two-r4(8,4)*four+r3(4,2)*two+r3(4,4)*two+r3(6,10)*two &
&                -r2(5,10)*four+r1(1,2)*two+r1(1,6)*two)*qx+(+r4(13,8)-r3(8,8)*two+r2(2,4) &
&                +r2(2,8)+r2(3,16)-r1(3,12)*two+r0(4)+r0(9))*xx
      eri(3,1,3)= r6(15)-r5(10,1)*two+r4(6,1)+r4(6,2)+r4(6,5)+r4(15,5)-r3(3,5)*two &
&                -r3(10,5)*two+r2(1,1)+r2(3,1)+r2(1,5)+r2(3,5)+r2(3,13)-r1(3,9)*two+r0(1) &
&                +r0(6)+(+r5(15,3)*two-r4(10,4)*four+r3(6,2)*two+r3(6,4)*two+r3(6,10)*two &
&                -r2(5,10)*four+r1(1,2)*two+r1(1,6)*two)*qx+(+r5(10,2)*two-r4(6,3)*four &
&                +r3(3,1)*two+r3(3,3)*two+r3(10,9)*two-r2(3,9)*four+r1(3,1)*two+r1(3,5)*two) &
&                *qz+(+r4(15,8)-r3(10,8)*two+r2(3,4)+r2(3,8)+r2(3,16)-r1(3,12)*two+r0(4) &
&                +r0(9))*xx+(+r4(10,7)*four-r3(6,7)*eight+r2(5,3)*four+r2(5,7)*four)*xz+( &
&                +r4(6,6)-r3(3,6)*two+r2(1,2)+r2(1,6)+r2(3,14)-r1(3,10)*two+r0(2)+r0(7))*zz+( &
&                +r3(10,12)*two-r2(3,12)*four+r1(3,4)*two+r1(3,8)*two)*xxz+(+r3(6,11)*two &
&                -r2(5,11)*four+r1(1,3)*two+r1(1,7)*two)*xzz+(+r2(3,17)-r1(3,13)*two+r0(5) &
&                +r0(10))*xxzz
      eri(4,1,3)= r6(9)-r5(5,1)*two+r4(2,1)+r4(2,2)+r4(9,5)*three-r3(5,5)*six+r2(4,1)*three &
&                +r2(4,5)*three+(+r5(9,2)+r5(9,3)*two-r4(5,3)*two-r4(5,4)*four+r3(2,1) &
&                +r3(2,2)*two+r3(2,3)+r3(2,4)*two+r3(9,9)+r3(9,10)*two-r2(6,9)*two &
&                -r2(6,10)*four+r1(2,1)+r1(2,2)*two+r1(2,5)+r1(2,6)*two)*qx+(+r4(9,7)*two &
&                +r4(9,8)-r3(5,7)*four-r3(5,8)*two+r2(4,3)*two+r2(4,4)+r2(4,7)*two+r2(4,8)) &
&                *xx+(+r3(9,12)-r2(6,12)*two+r1(2,4)+r1(2,8))*xxx
      eri(5,1,3)= r6(10)-r5(6,1)*two+r4(3,1)+r4(3,2)+r4(10,5)*three-r3(6,5)*six &
&                +r2(5,1)*three+r2(5,5)*three+(+r5(10,2)+r5(10,3)*two-r4(6,3)*two &
&                -r4(6,4)*four+r3(3,1)+r3(3,2)*two+r3(3,3)+r3(3,4)*two+r3(10,9)+r3(10,10)*two &
&                -r2(3,9)*two-r2(3,10)*four+r1(3,1)+r1(3,2)*two+r1(3,5)+r1(3,6)*two)*qx+( &
&                +r5(6,2)-r4(3,3)*two+r3(1,1)+r3(1,3)+r3(6,9)*three-r2(5,9)*six+r1(1,1)*three &
&                +r1(1,5)*three)*qz+(+r4(10,7)*two+r4(10,8)-r3(6,7)*four-r3(6,8)*two &
&                +r2(5,3)*two+r2(5,4)+r2(5,7)*two+r2(5,8))*xx+(+r4(6,6)+r4(6,7)*two &
&                -r3(3,6)*two-r3(3,7)*four+r2(1,2)+r2(1,3)*two+r2(1,6)+r2(1,7)*two+r2(3,14) &
&                +r2(3,15)*two-r1(3,10)*two-r1(3,11)*four+r0(2)+r0(3)*two+r0(7)+r0(8)*two)*xz &
&                +(+r3(10,12)-r2(3,12)*two+r1(3,4)+r1(3,8))*xxx+(+r3(6,11)*two+r3(6,12) &
&                -r2(5,11)*four-r2(5,12)*two+r1(1,3)*two+r1(1,4)+r1(1,7)*two+r1(1,8))*xxz+( &
&                +r2(3,17)-r1(3,13)*two+r0(5)+r0(10))*xxxz
      eri(6,1,3)= r6(14)-r5(9,1)*two+r4(5,1)+r4(5,2)+r4(14,5)-r3(9,5)*two+r2(6,1)+r2(6,5)+( &
&                +r5(14,3)*two-r4(9,4)*four+r3(5,2)*two+r3(5,4)*two)*qx+(+r5(9,2)-r4(5,3)*two &
&                +r3(2,1)+r3(2,3)+r3(9,9)-r2(6,9)*two+r1(2,1)+r1(2,5))*qz+(+r4(14,8) &
&                -r3(9,8)*two+r2(6,4)+r2(6,8))*xx+(+r4(9,7)*two-r3(5,7)*four+r2(4,3)*two &
&                +r2(4,7)*two)*xz+(+r3(9,12)-r2(6,12)*two+r1(2,4)+r1(2,8))*xxz
      eri(1,2,3)= r6(13)-r5(8,1)*two+r4(4,1)+r4(4,2)+r4(6,5)+r4(13,5)-r3(3,5)*two &
&                -r3(8,5)*two+r2(1,1)+r2(2,1)+r2(1,5)+r2(2,5)+r2(3,13)-r1(3,9)*two+r0(1) &
&                +r0(6)+(+r5(13,2)*two-r4(8,3)*four+r3(4,1)*two+r3(4,3)*two+r3(6,9)*two &
&                -r2(5,9)*four+r1(1,1)*two+r1(1,5)*two)*qx+(+r4(13,6)-r3(8,6)*two+r2(2,2) &
&                +r2(2,6)+r2(3,14)-r1(3,10)*two+r0(2)+r0(7))*xx
      eri(2,2,3)= r6(24)-r5(17,1)*two+r4(11,1)+r4(11,2)+r4(13,5)*six-r3(8,5)*p12 &
&                +r2(2,1)*six+r2(2,5)*six+r2(3,13)*three-r1(3,9)*six+r0(1)*three+r0(6)*three
      eri(3,2,3)= r6(26)-r5(19,1)*two+r4(13,1)+r4(13,2)+r4(13,5)+r4(15,5)-r3(8,5)*two &
&                -r3(10,5)*two+r2(2,1)+r2(3,1)+r2(2,5)+r2(3,5)+r2(3,13)-r1(3,9)*two+r0(1) &
&                +r0(6)+(+r5(19,2)*two-r4(13,3)*four+r3(8,1)*two+r3(8,3)*two+r3(10,9)*two &
&                -r2(3,9)*four+r1(3,1)*two+r1(3,5)*two)*qz+(+r4(13,6)-r3(8,6)*two+r2(2,2) &
&                +r2(2,6)+r2(3,14)-r1(3,10)*two+r0(2)+r0(7))*zz
      eri(4,2,3)= r6(18)-r5(12,1)*two+r4(7,1)+r4(7,2)+r4(9,5)*three-r3(5,5)*six &
&                +r2(4,1)*three+r2(4,5)*three+(+r5(18,2)-r4(12,3)*two+r3(7,1)+r3(7,3) &
&                +r3(9,9)*three-r2(6,9)*six+r1(2,1)*three+r1(2,5)*three)*qx
      eri(5,2,3)= r6(19)-r5(13,1)*two+r4(8,1)+r4(8,2)+r4(10,5)-r3(6,5)*two+r2(5,1)+r2(5,5) &
&                +(+r5(19,2)-r4(13,3)*two+r3(8,1)+r3(8,3)+r3(10,9)-r2(3,9)*two+r1(3,1) &
&                +r1(3,5))*qx+(+r5(13,2)-r4(8,3)*two+r3(4,1)+r3(4,3)+r3(6,9)-r2(5,9)*two &
&                +r1(1,1)+r1(1,5))*qz+(+r4(13,6)-r3(8,6)*two+r2(2,2)+r2(2,6)+r2(3,14) &
&                -r1(3,10)*two+r0(2)+r0(7))*xz
      eri(6,2,3)= r6(25)-r5(18,1)*two+r4(12,1)+r4(12,2)+r4(14,5)*three-r3(9,5)*six &
&                +r2(6,1)*three+r2(6,5)*three+(+r5(18,2)-r4(12,3)*two+r3(7,1)+r3(7,3) &
&                +r3(9,9)*three-r2(6,9)*six+r1(2,1)*three+r1(2,5)*three)*qz
      eri(1,3,3)= r6(15)-r5(10,1)*two+r4(6,1)+r4(6,2)+r4(6,5)+r4(15,5)-r3(3,5)*two &
&                -r3(10,5)*two+r2(1,1)+r2(3,1)+r2(1,5)+r2(3,5)+r2(3,13)-r1(3,9)*two+r0(1) &
&                +r0(6)+(+r5(15,2)*two-r4(10,3)*four+r3(6,1)*two+r3(6,3)*two+r3(6,9)*two &
&                -r2(5,9)*four+r1(1,1)*two+r1(1,5)*two)*qx+(+r5(10,3)*two-r4(6,4)*four &
&                +r3(3,2)*two+r3(3,4)*two+r3(10,10)*two-r2(3,10)*four+r1(3,2)*two+r1(3,6)*two &
&                )*qz+(+r4(15,6)-r3(10,6)*two+r2(3,2)+r2(3,6)+r2(3,14)-r1(3,10)*two+r0(2) &
&                +r0(7))*xx+(+r4(10,7)*four-r3(6,7)*eight+r2(5,3)*four+r2(5,7)*four)*xz+( &
&                +r4(6,8)-r3(3,8)*two+r2(1,4)+r2(1,8)+r2(3,16)-r1(3,12)*two+r0(4)+r0(9))*zz+( &
&                +r3(10,11)*two-r2(3,11)*four+r1(3,3)*two+r1(3,7)*two)*xxz+(+r3(6,12)*two &
&                -r2(5,12)*four+r1(1,4)*two+r1(1,8)*two)*xzz+(+r2(3,17)-r1(3,13)*two+r0(5) &
&                +r0(10))*xxzz
      eri(2,3,3)= r6(26)-r5(19,1)*two+r4(13,1)+r4(13,2)+r4(13,5)+r4(15,5)-r3(8,5)*two &
&                -r3(10,5)*two+r2(2,1)+r2(3,1)+r2(2,5)+r2(3,5)+r2(3,13)-r1(3,9)*two+r0(1) &
&                +r0(6)+(+r5(19,3)*two-r4(13,4)*four+r3(8,2)*two+r3(8,4)*two+r3(10,10)*two &
&                -r2(3,10)*four+r1(3,2)*two+r1(3,6)*two)*qz+(+r4(13,8)-r3(8,8)*two+r2(2,4) &
&                +r2(2,8)+r2(3,16)-r1(3,12)*two+r0(4)+r0(9))*zz
      eri(3,3,3)= r6(28)-r5(21,1)*two+r4(15,1)+r4(15,2)+r4(15,5)*six-r3(10,5)*p12 &
&                +r2(3,1)*six+r2(3,5)*six+r2(3,13)*three-r1(3,9)*six+r0(1)*three+r0(6)*three &
&                +(+r5(21,2)*two+r5(21,3)*two-r4(15,3)*four-r4(15,4)*four+r3(10,1)*two &
&                +r3(10,2)*two+r3(10,3)*two+r3(10,4)*two+r3(10,9)*six+r3(10,10)*six &
&                -r2(3,9)*p12-r2(3,10)*p12+r1(3,1)*six+r1(3,2)*six+r1(3,5)*six+r1(3,6)*six) &
&                *qz+(+r4(15,6)+r4(15,7)*four+r4(15,8)-r3(10,6)*two-r3(10,7)*eight &
&                -r3(10,8)*two+r2(3,2)+r2(3,3)*four+r2(3,4)+r2(3,6)+r2(3,7)*four+r2(3,8) &
&                +r2(3,14)+r2(3,15)*four+r2(3,16)-r1(3,10)*two-r1(3,11)*eight-r1(3,12)*two &
&                +r0(2)+r0(3)*four+r0(4)+r0(7)+r0(8)*four+r0(9))*zz+(+r3(10,11)*two &
&                +r3(10,12)*two-r2(3,11)*four-r2(3,12)*four+r1(3,3)*two+r1(3,4)*two &
&                +r1(3,7)*two+r1(3,8)*two)*zzz+(+r2(3,17)-r1(3,13)*two+r0(5)+r0(10))*zzzz
      eri(4,3,3)= r6(20)-r5(14,1)*two+r4(9,1)+r4(9,2)+r4(9,5)-r3(5,5)*two+r2(4,1)+r2(4,5)+( &
&                +r5(20,2)-r4(14,3)*two+r3(9,1)+r3(9,3)+r3(9,9)-r2(6,9)*two+r1(2,1)+r1(2,5)) &
&                *qx+(+r5(14,3)*two-r4(9,4)*four+r3(5,2)*two+r3(5,4)*two)*qz+(+r4(14,7)*two &
&                -r3(9,7)*four+r2(6,3)*two+r2(6,7)*two)*xz+(+r4(9,8)-r3(5,8)*two+r2(4,4) &
&                +r2(4,8))*zz+(+r3(9,12)-r2(6,12)*two+r1(2,4)+r1(2,8))*xzz
      eri(5,3,3)= r6(21)-r5(15,1)*two+r4(10,1)+r4(10,2)+r4(10,5)*three-r3(6,5)*six &
&                +r2(5,1)*three+r2(5,5)*three+(+r5(21,2)-r4(15,3)*two+r3(10,1)+r3(10,3) &
&                +r3(10,9)*three-r2(3,9)*six+r1(3,1)*three+r1(3,5)*three)*qx+(+r5(15,2) &
&                +r5(15,3)*two-r4(10,3)*two-r4(10,4)*four+r3(6,1)+r3(6,2)*two+r3(6,3) &
&                +r3(6,4)*two+r3(6,9)+r3(6,10)*two-r2(5,9)*two-r2(5,10)*four+r1(1,1) &
&                +r1(1,2)*two+r1(1,5)+r1(1,6)*two)*qz+(+r4(15,6)+r4(15,7)*two-r3(10,6)*two &
&                -r3(10,7)*four+r2(3,2)+r2(3,3)*two+r2(3,6)+r2(3,7)*two+r2(3,14)+r2(3,15)*two &
&                -r1(3,10)*two-r1(3,11)*four+r0(2)+r0(3)*two+r0(7)+r0(8)*two)*xz+( &
&                +r4(10,7)*two+r4(10,8)-r3(6,7)*four-r3(6,8)*two+r2(5,3)*two+r2(5,4) &
&                +r2(5,7)*two+r2(5,8))*zz+(+r3(10,11)*two+r3(10,12)-r2(3,11)*four &
&                -r2(3,12)*two+r1(3,3)*two+r1(3,4)+r1(3,7)*two+r1(3,8))*xzz+(+r3(6,12) &
&                -r2(5,12)*two+r1(1,4)+r1(1,8))*zzz+(+r2(3,17)-r1(3,13)*two+r0(5)+r0(10)) &
&                *xzzz
      eri(6,3,3)= r6(27)-r5(20,1)*two+r4(14,1)+r4(14,2)+r4(14,5)*three-r3(9,5)*six &
&                +r2(6,1)*three+r2(6,5)*three+(+r5(20,2)+r5(20,3)*two-r4(14,3)*two &
&                -r4(14,4)*four+r3(9,1)+r3(9,2)*two+r3(9,3)+r3(9,4)*two+r3(9,9)+r3(9,10)*two &
&                -r2(6,9)*two-r2(6,10)*four+r1(2,1)+r1(2,2)*two+r1(2,5)+r1(2,6)*two)*qz+( &
&                +r4(14,7)*two+r4(14,8)-r3(9,7)*four-r3(9,8)*two+r2(6,3)*two+r2(6,4) &
&                +r2(6,7)*two+r2(6,8))*zz+(+r3(9,12)-r2(6,12)*two+r1(2,4)+r1(2,8))*zzz
      eri(1,4,3)= r6(9)-r5(5,1)*two+r4(2,1)+r4(2,2)+r4(9,5)*three-r3(5,5)*six+r2(4,1)*three &
&                +r2(4,5)*three+(+r5(9,2)*two+r5(9,3)-r4(5,3)*four-r4(5,4)*two+r3(2,1)*two &
&                +r3(2,2)+r3(2,3)*two+r3(2,4)+r3(9,9)*two+r3(9,10)-r2(6,9)*four-r2(6,10)*two &
&                +r1(2,1)*two+r1(2,2)+r1(2,5)*two+r1(2,6))*qx+(+r4(9,6)+r4(9,7)*two &
&                -r3(5,6)*two-r3(5,7)*four+r2(4,2)+r2(4,3)*two+r2(4,6)+r2(4,7)*two)*xx+( &
&                +r3(9,11)-r2(6,11)*two+r1(2,3)+r1(2,7))*xxx
      eri(2,4,3)= r6(18)-r5(12,1)*two+r4(7,1)+r4(7,2)+r4(9,5)*three-r3(5,5)*six &
&                +r2(4,1)*three+r2(4,5)*three+(+r5(18,3)-r4(12,4)*two+r3(7,2)+r3(7,4) &
&                +r3(9,10)*three-r2(6,10)*six+r1(2,2)*three+r1(2,6)*three)*qx
      eri(3,4,3)= r6(20)-r5(14,1)*two+r4(9,1)+r4(9,2)+r4(9,5)-r3(5,5)*two+r2(4,1)+r2(4,5)+( &
&                +r5(20,3)-r4(14,4)*two+r3(9,2)+r3(9,4)+r3(9,10)-r2(6,10)*two+r1(2,2)+r1(2,6) &
&                )*qx+(+r5(14,2)*two-r4(9,3)*four+r3(5,1)*two+r3(5,3)*two)*qz+(+r4(14,7)*two &
&                -r3(9,7)*four+r2(6,3)*two+r2(6,7)*two)*xz+(+r4(9,6)-r3(5,6)*two+r2(4,2) &
&                +r2(4,6))*zz+(+r3(9,11)-r2(6,11)*two+r1(2,3)+r1(2,7))*xzz
      eri(4,4,3)= r6(13)-r5(8,1)*two+r4(4,1)+r4(4,2)+r4(6,5)+r4(13,5)-r3(3,5)*two &
&                -r3(8,5)*two+r2(1,1)+r2(2,1)+r2(1,5)+r2(2,5)+r2(3,13)-r1(3,9)*two+r0(1) &
&                +r0(6)+(+r5(13,2)+r5(13,3)-r4(8,3)*two-r4(8,4)*two+r3(4,1)+r3(4,2)+r3(4,3) &
&                +r3(4,4)+r3(6,9)+r3(6,10)-r2(5,9)*two-r2(5,10)*two+r1(1,1)+r1(1,2)+r1(1,5) &
&                +r1(1,6))*qx+(+r4(13,7)-r3(8,7)*two+r2(2,3)+r2(2,7)+r2(3,15)-r1(3,11)*two &
&                +r0(3)+r0(8))*xx
      eri(5,4,3)= r6(14)-r5(9,1)*two+r4(5,1)+r4(5,2)+r4(14,5)-r3(9,5)*two+r2(6,1)+r2(6,5)+( &
&                +r5(14,2)+r5(14,3)-r4(9,3)*two-r4(9,4)*two+r3(5,1)+r3(5,2)+r3(5,3)+r3(5,4)) &
&                *qx+(+r5(9,2)-r4(5,3)*two+r3(2,1)+r3(2,3)+r3(9,9)-r2(6,9)*two+r1(2,1) &
&                +r1(2,5))*qz+(+r4(14,7)-r3(9,7)*two+r2(6,3)+r2(6,7))*xx+(+r4(9,6)+r4(9,7) &
&                -r3(5,6)*two-r3(5,7)*two+r2(4,2)+r2(4,3)+r2(4,6)+r2(4,7))*xz+(+r3(9,11) &
&                -r2(6,11)*two+r1(2,3)+r1(2,7))*xxz
      eri(6,4,3)= r6(19)-r5(13,1)*two+r4(8,1)+r4(8,2)+r4(10,5)-r3(6,5)*two+r2(5,1)+r2(5,5) &
&                +(+r5(19,3)-r4(13,4)*two+r3(8,2)+r3(8,4)+r3(10,10)-r2(3,10)*two+r1(3,2) &
&                +r1(3,6))*qx+(+r5(13,2)-r4(8,3)*two+r3(4,1)+r3(4,3)+r3(6,9)-r2(5,9)*two &
&                +r1(1,1)+r1(1,5))*qz+(+r4(13,7)-r3(8,7)*two+r2(2,3)+r2(2,7)+r2(3,15) &
&                -r1(3,11)*two+r0(3)+r0(8))*xz
      eri(1,5,3)= r6(10)-r5(6,1)*two+r4(3,1)+r4(3,2)+r4(10,5)*three-r3(6,5)*six &
&                +r2(5,1)*three+r2(5,5)*three+(+r5(10,2)*two+r5(10,3)-r4(6,3)*four &
&                -r4(6,4)*two+r3(3,1)*two+r3(3,2)+r3(3,3)*two+r3(3,4)+r3(10,9)*two+r3(10,10) &
&                -r2(3,9)*four-r2(3,10)*two+r1(3,1)*two+r1(3,2)+r1(3,5)*two+r1(3,6))*qx+( &
&                +r5(6,3)-r4(3,4)*two+r3(1,2)+r3(1,4)+r3(6,10)*three-r2(5,10)*six &
&                +r1(1,2)*three+r1(1,6)*three)*qz+(+r4(10,6)+r4(10,7)*two-r3(6,6)*two &
&                -r3(6,7)*four+r2(5,2)+r2(5,3)*two+r2(5,6)+r2(5,7)*two)*xx+(+r4(6,7)*two &
&                +r4(6,8)-r3(3,7)*four-r3(3,8)*two+r2(1,3)*two+r2(1,4)+r2(1,7)*two+r2(1,8) &
&                +r2(3,15)*two+r2(3,16)-r1(3,11)*four-r1(3,12)*two+r0(3)*two+r0(4)+r0(8)*two &
&                +r0(9))*xz+(+r3(10,11)-r2(3,11)*two+r1(3,3)+r1(3,7))*xxx+(+r3(6,11) &
&                +r3(6,12)*two-r2(5,11)*two-r2(5,12)*four+r1(1,3)+r1(1,4)*two+r1(1,7) &
&                +r1(1,8)*two)*xxz+(+r2(3,17)-r1(3,13)*two+r0(5)+r0(10))*xxxz
      eri(2,5,3)= r6(19)-r5(13,1)*two+r4(8,1)+r4(8,2)+r4(10,5)-r3(6,5)*two+r2(5,1)+r2(5,5) &
&                +(+r5(19,3)-r4(13,4)*two+r3(8,2)+r3(8,4)+r3(10,10)-r2(3,10)*two+r1(3,2) &
&                +r1(3,6))*qx+(+r5(13,3)-r4(8,4)*two+r3(4,2)+r3(4,4)+r3(6,10)-r2(5,10)*two &
&                +r1(1,2)+r1(1,6))*qz+(+r4(13,8)-r3(8,8)*two+r2(2,4)+r2(2,8)+r2(3,16) &
&                -r1(3,12)*two+r0(4)+r0(9))*xz
      eri(3,5,3)= r6(21)-r5(15,1)*two+r4(10,1)+r4(10,2)+r4(10,5)*three-r3(6,5)*six &
&                +r2(5,1)*three+r2(5,5)*three+(+r5(21,3)-r4(15,4)*two+r3(10,2)+r3(10,4) &
&                +r3(10,10)*three-r2(3,10)*six+r1(3,2)*three+r1(3,6)*three)*qx+(+r5(15,2)*two &
&                +r5(15,3)-r4(10,3)*four-r4(10,4)*two+r3(6,1)*two+r3(6,2)+r3(6,3)*two+r3(6,4) &
&                +r3(6,9)*two+r3(6,10)-r2(5,9)*four-r2(5,10)*two+r1(1,1)*two+r1(1,2) &
&                +r1(1,5)*two+r1(1,6))*qz+(+r4(15,7)*two+r4(15,8)-r3(10,7)*four-r3(10,8)*two &
&                +r2(3,3)*two+r2(3,4)+r2(3,7)*two+r2(3,8)+r2(3,15)*two+r2(3,16)-r1(3,11)*four &
&                -r1(3,12)*two+r0(3)*two+r0(4)+r0(8)*two+r0(9))*xz+(+r4(10,6)+r4(10,7)*two &
&                -r3(6,6)*two-r3(6,7)*four+r2(5,2)+r2(5,3)*two+r2(5,6)+r2(5,7)*two)*zz+( &
&                +r3(10,11)+r3(10,12)*two-r2(3,11)*two-r2(3,12)*four+r1(3,3)+r1(3,4)*two &
&                +r1(3,7)+r1(3,8)*two)*xzz+(+r3(6,11)-r2(5,11)*two+r1(1,3)+r1(1,7))*zzz+( &
&                +r2(3,17)-r1(3,13)*two+r0(5)+r0(10))*xzzz
      eri(4,5,3)= r6(14)-r5(9,1)*two+r4(5,1)+r4(5,2)+r4(14,5)-r3(9,5)*two+r2(6,1)+r2(6,5)+( &
&                +r5(14,2)+r5(14,3)-r4(9,3)*two-r4(9,4)*two+r3(5,1)+r3(5,2)+r3(5,3)+r3(5,4)) &
&                *qx+(+r5(9,3)-r4(5,4)*two+r3(2,2)+r3(2,4)+r3(9,10)-r2(6,10)*two+r1(2,2) &
&                +r1(2,6))*qz+(+r4(14,7)-r3(9,7)*two+r2(6,3)+r2(6,7))*xx+(+r4(9,7)+r4(9,8) &
&                -r3(5,7)*two-r3(5,8)*two+r2(4,3)+r2(4,4)+r2(4,7)+r2(4,8))*xz+(+r3(9,12) &
&                -r2(6,12)*two+r1(2,4)+r1(2,8))*xxz
      eri(5,5,3)= r6(15)-r5(10,1)*two+r4(6,1)+r4(6,2)+r4(6,5)+r4(15,5)-r3(3,5)*two &
&                -r3(10,5)*two+r2(1,1)+r2(3,1)+r2(1,5)+r2(3,5)+r2(3,13)-r1(3,9)*two+r0(1) &
&                +r0(6)+(+r5(15,2)+r5(15,3)-r4(10,3)*two-r4(10,4)*two+r3(6,1)+r3(6,2)+r3(6,3) &
&                +r3(6,4)+r3(6,9)+r3(6,10)-r2(5,9)*two-r2(5,10)*two+r1(1,1)+r1(1,2)+r1(1,5) &
&                +r1(1,6))*qx+(+r5(10,2)+r5(10,3)-r4(6,3)*two-r4(6,4)*two+r3(3,1)+r3(3,2) &
&                +r3(3,3)+r3(3,4)+r3(10,9)+r3(10,10)-r2(3,9)*two-r2(3,10)*two+r1(3,1)+r1(3,2) &
&                +r1(3,5)+r1(3,6))*qz+(+r4(15,7)-r3(10,7)*two+r2(3,3)+r2(3,7)+r2(3,15) &
&                -r1(3,11)*two+r0(3)+r0(8))*xx+(+r4(10,6)+r4(10,7)*two+r4(10,8)-r3(6,6)*two &
&                -r3(6,7)*four-r3(6,8)*two+r2(5,2)+r2(5,3)*two+r2(5,4)+r2(5,6)+r2(5,7)*two &
&                +r2(5,8))*xz+(+r4(6,7)-r3(3,7)*two+r2(1,3)+r2(1,7)+r2(3,15)-r1(3,11)*two &
&                +r0(3)+r0(8))*zz+(+r3(10,11)+r3(10,12)-r2(3,11)*two-r2(3,12)*two+r1(3,3) &
&                +r1(3,4)+r1(3,7)+r1(3,8))*xxz+(+r3(6,11)+r3(6,12)-r2(5,11)*two-r2(5,12)*two &
&                +r1(1,3)+r1(1,4)+r1(1,7)+r1(1,8))*xzz+(+r2(3,17)-r1(3,13)*two+r0(5)+r0(10)) &
&                *xxzz
      eri(6,5,3)= r6(20)-r5(14,1)*two+r4(9,1)+r4(9,2)+r4(9,5)-r3(5,5)*two+r2(4,1)+r2(4,5)+( &
&                +r5(20,3)-r4(14,4)*two+r3(9,2)+r3(9,4)+r3(9,10)-r2(6,10)*two+r1(2,2)+r1(2,6) &
&                )*qx+(+r5(14,2)+r5(14,3)-r4(9,3)*two-r4(9,4)*two+r3(5,1)+r3(5,2)+r3(5,3) &
&                +r3(5,4))*qz+(+r4(14,7)+r4(14,8)-r3(9,7)*two-r3(9,8)*two+r2(6,3)+r2(6,4) &
&                +r2(6,7)+r2(6,8))*xz+(+r4(9,7)-r3(5,7)*two+r2(4,3)+r2(4,7))*zz+(+r3(9,12) &
&                -r2(6,12)*two+r1(2,4)+r1(2,8))*xzz
      eri(1,6,3)= r6(14)-r5(9,1)*two+r4(5,1)+r4(5,2)+r4(14,5)-r3(9,5)*two+r2(6,1)+r2(6,5)+( &
&                +r5(14,2)*two-r4(9,3)*four+r3(5,1)*two+r3(5,3)*two)*qx+(+r5(9,3)-r4(5,4)*two &
&                +r3(2,2)+r3(2,4)+r3(9,10)-r2(6,10)*two+r1(2,2)+r1(2,6))*qz+(+r4(14,6) &
&                -r3(9,6)*two+r2(6,2)+r2(6,6))*xx+(+r4(9,7)*two-r3(5,7)*four+r2(4,3)*two &
&                +r2(4,7)*two)*xz+(+r3(9,11)-r2(6,11)*two+r1(2,3)+r1(2,7))*xxz
      eri(2,6,3)= r6(25)-r5(18,1)*two+r4(12,1)+r4(12,2)+r4(14,5)*three-r3(9,5)*six &
&                +r2(6,1)*three+r2(6,5)*three+(+r5(18,3)-r4(12,4)*two+r3(7,2)+r3(7,4) &
&                +r3(9,10)*three-r2(6,10)*six+r1(2,2)*three+r1(2,6)*three)*qz
      eri(3,6,3)= r6(27)-r5(20,1)*two+r4(14,1)+r4(14,2)+r4(14,5)*three-r3(9,5)*six &
&                +r2(6,1)*three+r2(6,5)*three+(+r5(20,2)*two+r5(20,3)-r4(14,3)*four &
&                -r4(14,4)*two+r3(9,1)*two+r3(9,2)+r3(9,3)*two+r3(9,4)+r3(9,9)*two+r3(9,10) &
&                -r2(6,9)*four-r2(6,10)*two+r1(2,1)*two+r1(2,2)+r1(2,5)*two+r1(2,6))*qz+( &
&                +r4(14,6)+r4(14,7)*two-r3(9,6)*two-r3(9,7)*four+r2(6,2)+r2(6,3)*two+r2(6,6) &
&                +r2(6,7)*two)*zz+(+r3(9,11)-r2(6,11)*two+r1(2,3)+r1(2,7))*zzz
      eri(4,6,3)= r6(19)-r5(13,1)*two+r4(8,1)+r4(8,2)+r4(10,5)-r3(6,5)*two+r2(5,1)+r2(5,5) &
&                +(+r5(19,2)-r4(13,3)*two+r3(8,1)+r3(8,3)+r3(10,9)-r2(3,9)*two+r1(3,1) &
&                +r1(3,5))*qx+(+r5(13,3)-r4(8,4)*two+r3(4,2)+r3(4,4)+r3(6,10)-r2(5,10)*two &
&                +r1(1,2)+r1(1,6))*qz+(+r4(13,7)-r3(8,7)*two+r2(2,3)+r2(2,7)+r2(3,15) &
&                -r1(3,11)*two+r0(3)+r0(8))*xz
      eri(5,6,3)= r6(20)-r5(14,1)*two+r4(9,1)+r4(9,2)+r4(9,5)-r3(5,5)*two+r2(4,1)+r2(4,5)+( &
&                +r5(20,2)-r4(14,3)*two+r3(9,1)+r3(9,3)+r3(9,9)-r2(6,9)*two+r1(2,1)+r1(2,5)) &
&                *qx+(+r5(14,2)+r5(14,3)-r4(9,3)*two-r4(9,4)*two+r3(5,1)+r3(5,2)+r3(5,3) &
&                +r3(5,4))*qz+(+r4(14,6)+r4(14,7)-r3(9,6)*two-r3(9,7)*two+r2(6,2)+r2(6,3) &
&                +r2(6,6)+r2(6,7))*xz+(+r4(9,7)-r3(5,7)*two+r2(4,3)+r2(4,7))*zz+(+r3(9,11) &
&                -r2(6,11)*two+r1(2,3)+r1(2,7))*xzz
      eri(6,6,3)= r6(26)-r5(19,1)*two+r4(13,1)+r4(13,2)+r4(13,5)+r4(15,5)-r3(8,5)*two &
&                -r3(10,5)*two+r2(2,1)+r2(3,1)+r2(2,5)+r2(3,5)+r2(3,13)-r1(3,9)*two+r0(1) &
&                +r0(6)+(+r5(19,2)+r5(19,3)-r4(13,3)*two-r4(13,4)*two+r3(8,1)+r3(8,2)+r3(8,3) &
&                +r3(8,4)+r3(10,9)+r3(10,10)-r2(3,9)*two-r2(3,10)*two+r1(3,1)+r1(3,2)+r1(3,5) &
&                +r1(3,6))*qz+(+r4(13,7)-r3(8,7)*two+r2(2,3)+r2(2,7)+r2(3,15)-r1(3,11)*two &
&                +r0(3)+r0(8))*zz
      eri(1,1,4)= r6(2)+r4(2,5)*six+r2(4,13)*three+(+r5(2,2)*two+r5(2,3)*two+r3(2,9)*six &
&                +r3(2,10)*six)*qx+(+r4(2,6)+r4(2,7)*four+r4(2,8)+r2(4,14)+r2(4,15)*four &
&                +r2(4,16))*xx+(+r3(2,11)*two+r3(2,12)*two)*xxx+r2(4,17)*xxxx
      eri(2,1,4)= r6(7)+r4(2,5)+r4(7,5)+r2(4,13)+(+r5(7,3)*two+r3(2,10)*two)*qx+(+r4(7,8) &
&                +r2(4,16))*xx
      eri(3,1,4)= r6(9)+r4(2,5)+r4(9,5)+r2(4,13)+(+r5(9,3)*two+r3(2,10)*two)*qx+( &
&                +r5(5,2)*two+r3(5,9)*two)*qz+(+r4(9,8)+r2(4,16))*xx+r4(5,7)*four*xz+( &
&                +r4(2,6)+r2(4,14))*zz+r3(5,12)*two*xxz+r3(2,11)*two*xzz+r2(4,17)*xxzz
      eri(4,1,4)= r6(4)+r4(4,5)*three+(+r5(4,2)+r5(4,3)*two+r3(4,9)+r3(4,10)*two)*qx+( &
&                +r4(4,7)*two+r4(4,8))*xx+r3(4,12)*xxx
      eri(5,1,4)= r6(5)+r4(5,5)*three+(+r5(5,2)+r5(5,3)*two+r3(5,9)+r3(5,10)*two)*qx+( &
&                +r5(2,2)+r3(2,9)*three)*qz+(+r4(5,7)*two+r4(5,8))*xx+(+r4(2,6)+r4(2,7)*two &
&                +r2(4,14)+r2(4,15)*two)*xz+r3(5,12)*xxx+(+r3(2,11)*two+r3(2,12))*xxz &
&                +r2(4,17)*xxxz
      eri(6,1,4)= r6(8)+r4(8,5)+r5(8,3)*two*qx+(+r5(4,2)+r3(4,9))*qz+r4(8,8)*xx+r4(4,7)*two &
&                *xz+r3(4,12)*xxz
      eri(1,2,4)= r6(7)+r4(2,5)+r4(7,5)+r2(4,13)+(+r5(7,2)*two+r3(2,9)*two)*qx+(+r4(7,6) &
&                +r2(4,14))*xx
      eri(2,2,4)= r6(16)+r4(7,5)*six+r2(4,13)*three
      eri(3,2,4)= r6(18)+r4(7,5)+r4(9,5)+r2(4,13)+(+r5(12,2)*two+r3(5,9)*two)*qz+(+r4(7,6) &
&                +r2(4,14))*zz
      eri(4,2,4)= r6(11)+r4(4,5)*three+(+r5(11,2)+r3(4,9)*three)*qx
      eri(5,2,4)= r6(12)+r4(5,5)+(+r5(12,2)+r3(5,9))*qx+(+r5(7,2)+r3(2,9))*qz+(+r4(7,6) &
&                +r2(4,14))*xz
      eri(6,2,4)= r6(17)+r4(8,5)*three+(+r5(11,2)+r3(4,9)*three)*qz
      eri(1,3,4)= r6(9)+r4(2,5)+r4(9,5)+r2(4,13)+(+r5(9,2)*two+r3(2,9)*two)*qx+( &
&                +r5(5,3)*two+r3(5,10)*two)*qz+(+r4(9,6)+r2(4,14))*xx+r4(5,7)*four*xz+( &
&                +r4(2,8)+r2(4,16))*zz+r3(5,11)*two*xxz+r3(2,12)*two*xzz+r2(4,17)*xxzz
      eri(2,3,4)= r6(18)+r4(7,5)+r4(9,5)+r2(4,13)+(+r5(12,3)*two+r3(5,10)*two)*qz+(+r4(7,8) &
&                +r2(4,16))*zz
      eri(3,3,4)= r6(20)+r4(9,5)*six+r2(4,13)*three+(+r5(14,2)*two+r5(14,3)*two+r3(5,9)*six &
&                +r3(5,10)*six)*qz+(+r4(9,6)+r4(9,7)*four+r4(9,8)+r2(4,14)+r2(4,15)*four &
&                +r2(4,16))*zz+(+r3(5,11)*two+r3(5,12)*two)*zzz+r2(4,17)*zzzz
      eri(4,3,4)= r6(13)+r4(4,5)+(+r5(13,2)+r3(4,9))*qx+r5(8,3)*two*qz+r4(8,7)*two*xz &
&                +r4(4,8)*zz+r3(4,12)*xzz
      eri(5,3,4)= r6(14)+r4(5,5)*three+(+r5(14,2)+r3(5,9)*three)*qx+(+r5(9,2)+r5(9,3)*two &
&                +r3(2,9)+r3(2,10)*two)*qz+(+r4(9,6)+r4(9,7)*two+r2(4,14)+r2(4,15)*two)*xz+( &
&                +r4(5,7)*two+r4(5,8))*zz+(+r3(5,11)*two+r3(5,12))*xzz+r3(2,12)*zzz+r2(4,17) &
&                *xzzz
      eri(6,3,4)= r6(19)+r4(8,5)*three+(+r5(13,2)+r5(13,3)*two+r3(4,9)+r3(4,10)*two)*qz+( &
&                +r4(8,7)*two+r4(8,8))*zz+r3(4,12)*zzz
      eri(1,4,4)= r6(4)+r4(4,5)*three+(+r5(4,2)*two+r5(4,3)+r3(4,9)*two+r3(4,10))*qx+( &
&                +r4(4,6)+r4(4,7)*two)*xx+r3(4,11)*xxx
      eri(2,4,4)= r6(11)+r4(4,5)*three+(+r5(11,3)+r3(4,10)*three)*qx
      eri(3,4,4)= r6(13)+r4(4,5)+(+r5(13,3)+r3(4,10))*qx+r5(8,2)*two*qz+r4(8,7)*two*xz &
&                +r4(4,6)*zz+r3(4,11)*xzz
      eri(4,4,4)= r6(7)+r4(2,5)+r4(7,5)+r2(4,13)+(+r5(7,2)+r5(7,3)+r3(2,9)+r3(2,10))*qx+( &
&                +r4(7,7)+r2(4,15))*xx
      eri(5,4,4)= r6(8)+r4(8,5)+(+r5(8,2)+r5(8,3))*qx+(+r5(4,2)+r3(4,9))*qz+r4(8,7)*xx+( &
&                +r4(4,6)+r4(4,7))*xz+r3(4,11)*xxz
      eri(6,4,4)= r6(12)+r4(5,5)+(+r5(12,3)+r3(5,10))*qx+(+r5(7,2)+r3(2,9))*qz+(+r4(7,7) &
&                +r2(4,15))*xz
      eri(1,5,4)= r6(5)+r4(5,5)*three+(+r5(5,2)*two+r5(5,3)+r3(5,9)*two+r3(5,10))*qx+( &
&                +r5(2,3)+r3(2,10)*three)*qz+(+r4(5,6)+r4(5,7)*two)*xx+(+r4(2,7)*two+r4(2,8) &
&                +r2(4,15)*two+r2(4,16))*xz+r3(5,11)*xxx+(+r3(2,11)+r3(2,12)*two)*xxz &
&                +r2(4,17)*xxxz
      eri(2,5,4)= r6(12)+r4(5,5)+(+r5(12,3)+r3(5,10))*qx+(+r5(7,3)+r3(2,10))*qz+(+r4(7,8) &
&                +r2(4,16))*xz
      eri(3,5,4)= r6(14)+r4(5,5)*three+(+r5(14,3)+r3(5,10)*three)*qx+(+r5(9,2)*two+r5(9,3) &
&                +r3(2,9)*two+r3(2,10))*qz+(+r4(9,7)*two+r4(9,8)+r2(4,15)*two+r2(4,16))*xz+( &
&                +r4(5,6)+r4(5,7)*two)*zz+(+r3(5,11)+r3(5,12)*two)*xzz+r3(2,11)*zzz+r2(4,17) &
&                *xzzz
      eri(4,5,4)= r6(8)+r4(8,5)+(+r5(8,2)+r5(8,3))*qx+(+r5(4,3)+r3(4,10))*qz+r4(8,7)*xx+( &
&                +r4(4,7)+r4(4,8))*xz+r3(4,12)*xxz
      eri(5,5,4)= r6(9)+r4(2,5)+r4(9,5)+r2(4,13)+(+r5(9,2)+r5(9,3)+r3(2,9)+r3(2,10))*qx+( &
&                +r5(5,2)+r5(5,3)+r3(5,9)+r3(5,10))*qz+(+r4(9,7)+r2(4,15))*xx+(+r4(5,6) &
&                +r4(5,7)*two+r4(5,8))*xz+(+r4(2,7)+r2(4,15))*zz+(+r3(5,11)+r3(5,12))*xxz+( &
&                +r3(2,11)+r3(2,12))*xzz+r2(4,17)*xxzz
      eri(6,5,4)= r6(13)+r4(4,5)+(+r5(13,3)+r3(4,10))*qx+(+r5(8,2)+r5(8,3))*qz+(+r4(8,7) &
&                +r4(8,8))*xz+r4(4,7)*zz+r3(4,12)*xzz
      eri(1,6,4)= r6(8)+r4(8,5)+r5(8,2)*two*qx+(+r5(4,3)+r3(4,10))*qz+r4(8,6)*xx &
&                +r4(4,7)*two*xz+r3(4,11)*xxz
      eri(2,6,4)= r6(17)+r4(8,5)*three+(+r5(11,3)+r3(4,10)*three)*qz
      eri(3,6,4)= r6(19)+r4(8,5)*three+(+r5(13,2)*two+r5(13,3)+r3(4,9)*two+r3(4,10))*qz+( &
&                +r4(8,6)+r4(8,7)*two)*zz+r3(4,11)*zzz
      eri(4,6,4)= r6(12)+r4(5,5)+(+r5(12,2)+r3(5,9))*qx+(+r5(7,3)+r3(2,10))*qz+(+r4(7,7) &
&                +r2(4,15))*xz
      eri(5,6,4)= r6(13)+r4(4,5)+(+r5(13,2)+r3(4,9))*qx+(+r5(8,2)+r5(8,3))*qz+(+r4(8,6) &
&                +r4(8,7))*xz+r4(4,7)*zz+r3(4,11)*xzz
      eri(6,6,4)= r6(18)+r4(7,5)+r4(9,5)+r2(4,13)+(+r5(12,2)+r5(12,3)+r3(5,9)+r3(5,10))*qz &
&                +(+r4(7,7)+r2(4,15))*zz
      eri(1,1,5)= r6(3)-r5(1,1)+r4(3,5)*six-r3(1,5)*six+r2(5,13)*three-r1(1,9)*three+( &
&                +r5(3,2)*two+r5(3,3)*two-r4(1,3)*two-r4(1,4)*two+r3(3,9)*six+r3(3,10)*six &
&                -r2(1,9)*six-r2(1,10)*six)*qx+(+r4(3,6)+r4(3,7)*four+r4(3,8)-r3(1,6) &
&                -r3(1,7)*four-r3(1,8)+r2(5,14)+r2(5,15)*four+r2(5,16)-r1(1,10)-r1(1,11)*four &
&                -r1(1,12))*xx+(+r3(3,11)*two+r3(3,12)*two-r2(1,11)*two-r2(1,12)*two)*xxx+( &
&                +r2(5,17)-r1(1,13))*xxxx
      eri(2,1,5)= r6(8)-r5(4,1)+r4(3,5)+r4(8,5)-r3(1,5)-r3(4,5)+r2(5,13)-r1(1,9)+( &
&                +r5(8,3)*two-r4(4,4)*two+r3(3,10)*two-r2(1,10)*two)*qx+(+r4(8,8)-r3(4,8) &
&                +r2(5,16)-r1(1,12))*xx
      eri(3,1,5)= r6(10)-r5(6,1)+r4(3,5)+r4(10,5)-r3(1,5)-r3(6,5)+r2(5,13)-r1(1,9)+( &
&                +r5(10,3)*two-r4(6,4)*two+r3(3,10)*two-r2(1,10)*two)*qx+(+r5(6,2)*two &
&                -r4(3,3)*two+r3(6,9)*two-r2(5,9)*two)*qz+(+r4(10,8)-r3(6,8)+r2(5,16) &
&                -r1(1,12))*xx+(+r4(6,7)*four-r3(3,7)*four)*xz+(+r4(3,6)-r3(1,6)+r2(5,14) &
&                -r1(1,10))*zz+(+r3(6,12)*two-r2(5,12)*two)*xxz+(+r3(3,11)*two-r2(1,11)*two) &
&                *xzz+(+r2(5,17)-r1(1,13))*xxzz
      eri(4,1,5)= r6(5)-r5(2,1)+r4(5,5)*three-r3(2,5)*three+(+r5(5,2)+r5(5,3)*two-r4(2,3) &
&                -r4(2,4)*two+r3(5,9)+r3(5,10)*two-r2(4,9)-r2(4,10)*two)*qx+(+r4(5,7)*two &
&                +r4(5,8)-r3(2,7)*two-r3(2,8))*xx+(+r3(5,12)-r2(4,12))*xxx
      eri(5,1,5)= r6(6)-r5(3,1)+r4(6,5)*three-r3(3,5)*three+(+r5(6,2)+r5(6,3)*two-r4(3,3) &
&                -r4(3,4)*two+r3(6,9)+r3(6,10)*two-r2(5,9)-r2(5,10)*two)*qx+(+r5(3,2)-r4(1,3) &
&                +r3(3,9)*three-r2(1,9)*three)*qz+(+r4(6,7)*two+r4(6,8)-r3(3,7)*two-r3(3,8)) &
&                *xx+(+r4(3,6)+r4(3,7)*two-r3(1,6)-r3(1,7)*two+r2(5,14)+r2(5,15)*two-r1(1,10) &
&                -r1(1,11)*two)*xz+(+r3(6,12)-r2(5,12))*xxx+(+r3(3,11)*two+r3(3,12) &
&                -r2(1,11)*two-r2(1,12))*xxz+(+r2(5,17)-r1(1,13))*xxxz
      eri(6,1,5)= r6(9)-r5(5,1)+r4(9,5)-r3(5,5)+(+r5(9,3)*two-r4(5,4)*two)*qx+(+r5(5,2) &
&                -r4(2,3)+r3(5,9)-r2(4,9))*qz+(+r4(9,8)-r3(5,8))*xx+(+r4(5,7)*two-r3(2,7)*two &
&                )*xz+(+r3(5,12)-r2(4,12))*xxz
      eri(1,2,5)= r6(8)-r5(4,1)+r4(3,5)+r4(8,5)-r3(1,5)-r3(4,5)+r2(5,13)-r1(1,9)+( &
&                +r5(8,2)*two-r4(4,3)*two+r3(3,9)*two-r2(1,9)*two)*qx+(+r4(8,6)-r3(4,6) &
&                +r2(5,14)-r1(1,10))*xx
      eri(2,2,5)= r6(17)-r5(11,1)+r4(8,5)*six-r3(4,5)*six+r2(5,13)*three-r1(1,9)*three
      eri(3,2,5)= r6(19)-r5(13,1)+r4(8,5)+r4(10,5)-r3(4,5)-r3(6,5)+r2(5,13)-r1(1,9)+( &
&                +r5(13,2)*two-r4(8,3)*two+r3(6,9)*two-r2(5,9)*two)*qz+(+r4(8,6)-r3(4,6) &
&                +r2(5,14)-r1(1,10))*zz
      eri(4,2,5)= r6(12)-r5(7,1)+r4(5,5)*three-r3(2,5)*three+(+r5(12,2)-r4(7,3) &
&                +r3(5,9)*three-r2(4,9)*three)*qx
      eri(5,2,5)= r6(13)-r5(8,1)+r4(6,5)-r3(3,5)+(+r5(13,2)-r4(8,3)+r3(6,9)-r2(5,9))*qx+( &
&                +r5(8,2)-r4(4,3)+r3(3,9)-r2(1,9))*qz+(+r4(8,6)-r3(4,6)+r2(5,14)-r1(1,10))*xz
      eri(6,2,5)= r6(18)-r5(12,1)+r4(9,5)*three-r3(5,5)*three+(+r5(12,2)-r4(7,3) &
&                +r3(5,9)*three-r2(4,9)*three)*qz
      eri(1,3,5)= r6(10)-r5(6,1)+r4(3,5)+r4(10,5)-r3(1,5)-r3(6,5)+r2(5,13)-r1(1,9)+( &
&                +r5(10,2)*two-r4(6,3)*two+r3(3,9)*two-r2(1,9)*two)*qx+(+r5(6,3)*two &
&                -r4(3,4)*two+r3(6,10)*two-r2(5,10)*two)*qz+(+r4(10,6)-r3(6,6)+r2(5,14) &
&                -r1(1,10))*xx+(+r4(6,7)*four-r3(3,7)*four)*xz+(+r4(3,8)-r3(1,8)+r2(5,16) &
&                -r1(1,12))*zz+(+r3(6,11)*two-r2(5,11)*two)*xxz+(+r3(3,12)*two-r2(1,12)*two) &
&                *xzz+(+r2(5,17)-r1(1,13))*xxzz
      eri(2,3,5)= r6(19)-r5(13,1)+r4(8,5)+r4(10,5)-r3(4,5)-r3(6,5)+r2(5,13)-r1(1,9)+( &
&                +r5(13,3)*two-r4(8,4)*two+r3(6,10)*two-r2(5,10)*two)*qz+(+r4(8,8)-r3(4,8) &
&                +r2(5,16)-r1(1,12))*zz
      eri(3,3,5)= r6(21)-r5(15,1)+r4(10,5)*six-r3(6,5)*six+r2(5,13)*three-r1(1,9)*three+( &
&                +r5(15,2)*two+r5(15,3)*two-r4(10,3)*two-r4(10,4)*two+r3(6,9)*six &
&                +r3(6,10)*six-r2(5,9)*six-r2(5,10)*six)*qz+(+r4(10,6)+r4(10,7)*four+r4(10,8) &
&                -r3(6,6)-r3(6,7)*four-r3(6,8)+r2(5,14)+r2(5,15)*four+r2(5,16)-r1(1,10) &
&                -r1(1,11)*four-r1(1,12))*zz+(+r3(6,11)*two+r3(6,12)*two-r2(5,11)*two &
&                -r2(5,12)*two)*zzz+(+r2(5,17)-r1(1,13))*zzzz
      eri(4,3,5)= r6(14)-r5(9,1)+r4(5,5)-r3(2,5)+(+r5(14,2)-r4(9,3)+r3(5,9)-r2(4,9))*qx+( &
&                +r5(9,3)*two-r4(5,4)*two)*qz+(+r4(9,7)*two-r3(5,7)*two)*xz+(+r4(5,8)-r3(2,8) &
&                )*zz+(+r3(5,12)-r2(4,12))*xzz
      eri(5,3,5)= r6(15)-r5(10,1)+r4(6,5)*three-r3(3,5)*three+(+r5(15,2)-r4(10,3) &
&                +r3(6,9)*three-r2(5,9)*three)*qx+(+r5(10,2)+r5(10,3)*two-r4(6,3)-r4(6,4)*two &
&                +r3(3,9)+r3(3,10)*two-r2(1,9)-r2(1,10)*two)*qz+(+r4(10,6)+r4(10,7)*two &
&                -r3(6,6)-r3(6,7)*two+r2(5,14)+r2(5,15)*two-r1(1,10)-r1(1,11)*two)*xz+( &
&                +r4(6,7)*two+r4(6,8)-r3(3,7)*two-r3(3,8))*zz+(+r3(6,11)*two+r3(6,12) &
&                -r2(5,11)*two-r2(5,12))*xzz+(+r3(3,12)-r2(1,12))*zzz+(+r2(5,17)-r1(1,13)) &
&                *xzzz
      eri(6,3,5)= r6(20)-r5(14,1)+r4(9,5)*three-r3(5,5)*three+(+r5(14,2)+r5(14,3)*two &
&                -r4(9,3)-r4(9,4)*two+r3(5,9)+r3(5,10)*two-r2(4,9)-r2(4,10)*two)*qz+( &
&                +r4(9,7)*two+r4(9,8)-r3(5,7)*two-r3(5,8))*zz+(+r3(5,12)-r2(4,12))*zzz
      eri(1,4,5)= r6(5)-r5(2,1)+r4(5,5)*three-r3(2,5)*three+(+r5(5,2)*two+r5(5,3) &
&                -r4(2,3)*two-r4(2,4)+r3(5,9)*two+r3(5,10)-r2(4,9)*two-r2(4,10))*qx+(+r4(5,6) &
&                +r4(5,7)*two-r3(2,6)-r3(2,7)*two)*xx+(+r3(5,11)-r2(4,11))*xxx
      eri(2,4,5)= r6(12)-r5(7,1)+r4(5,5)*three-r3(2,5)*three+(+r5(12,3)-r4(7,4) &
&                +r3(5,10)*three-r2(4,10)*three)*qx
      eri(3,4,5)= r6(14)-r5(9,1)+r4(5,5)-r3(2,5)+(+r5(14,3)-r4(9,4)+r3(5,10)-r2(4,10))*qx+( &
&                +r5(9,2)*two-r4(5,3)*two)*qz+(+r4(9,7)*two-r3(5,7)*two)*xz+(+r4(5,6)-r3(2,6) &
&                )*zz+(+r3(5,11)-r2(4,11))*xzz
      eri(4,4,5)= r6(8)-r5(4,1)+r4(3,5)+r4(8,5)-r3(1,5)-r3(4,5)+r2(5,13)-r1(1,9)+(+r5(8,2) &
&                +r5(8,3)-r4(4,3)-r4(4,4)+r3(3,9)+r3(3,10)-r2(1,9)-r2(1,10))*qx+(+r4(8,7) &
&                -r3(4,7)+r2(5,15)-r1(1,11))*xx
      eri(5,4,5)= r6(9)-r5(5,1)+r4(9,5)-r3(5,5)+(+r5(9,2)+r5(9,3)-r4(5,3)-r4(5,4))*qx+( &
&                +r5(5,2)-r4(2,3)+r3(5,9)-r2(4,9))*qz+(+r4(9,7)-r3(5,7))*xx+(+r4(5,6)+r4(5,7) &
&                -r3(2,6)-r3(2,7))*xz+(+r3(5,11)-r2(4,11))*xxz
      eri(6,4,5)= r6(13)-r5(8,1)+r4(6,5)-r3(3,5)+(+r5(13,3)-r4(8,4)+r3(6,10)-r2(5,10))*qx+( &
&                +r5(8,2)-r4(4,3)+r3(3,9)-r2(1,9))*qz+(+r4(8,7)-r3(4,7)+r2(5,15)-r1(1,11))*xz
      eri(1,5,5)= r6(6)-r5(3,1)+r4(6,5)*three-r3(3,5)*three+(+r5(6,2)*two+r5(6,3) &
&                -r4(3,3)*two-r4(3,4)+r3(6,9)*two+r3(6,10)-r2(5,9)*two-r2(5,10))*qx+(+r5(3,3) &
&                -r4(1,4)+r3(3,10)*three-r2(1,10)*three)*qz+(+r4(6,6)+r4(6,7)*two-r3(3,6) &
&                -r3(3,7)*two)*xx+(+r4(3,7)*two+r4(3,8)-r3(1,7)*two-r3(1,8)+r2(5,15)*two &
&                +r2(5,16)-r1(1,11)*two-r1(1,12))*xz+(+r3(6,11)-r2(5,11))*xxx+(+r3(3,11) &
&                +r3(3,12)*two-r2(1,11)-r2(1,12)*two)*xxz+(+r2(5,17)-r1(1,13))*xxxz
      eri(2,5,5)= r6(13)-r5(8,1)+r4(6,5)-r3(3,5)+(+r5(13,3)-r4(8,4)+r3(6,10)-r2(5,10))*qx+( &
&                +r5(8,3)-r4(4,4)+r3(3,10)-r2(1,10))*qz+(+r4(8,8)-r3(4,8)+r2(5,16)-r1(1,12)) &
&                *xz
      eri(3,5,5)= r6(15)-r5(10,1)+r4(6,5)*three-r3(3,5)*three+(+r5(15,3)-r4(10,4) &
&                +r3(6,10)*three-r2(5,10)*three)*qx+(+r5(10,2)*two+r5(10,3)-r4(6,3)*two &
&                -r4(6,4)+r3(3,9)*two+r3(3,10)-r2(1,9)*two-r2(1,10))*qz+(+r4(10,7)*two &
&                +r4(10,8)-r3(6,7)*two-r3(6,8)+r2(5,15)*two+r2(5,16)-r1(1,11)*two-r1(1,12)) &
&                *xz+(+r4(6,6)+r4(6,7)*two-r3(3,6)-r3(3,7)*two)*zz+(+r3(6,11)+r3(6,12)*two &
&                -r2(5,11)-r2(5,12)*two)*xzz+(+r3(3,11)-r2(1,11))*zzz+(+r2(5,17)-r1(1,13)) &
&                *xzzz
      eri(4,5,5)= r6(9)-r5(5,1)+r4(9,5)-r3(5,5)+(+r5(9,2)+r5(9,3)-r4(5,3)-r4(5,4))*qx+( &
&                +r5(5,3)-r4(2,4)+r3(5,10)-r2(4,10))*qz+(+r4(9,7)-r3(5,7))*xx+(+r4(5,7) &
&                +r4(5,8)-r3(2,7)-r3(2,8))*xz+(+r3(5,12)-r2(4,12))*xxz
      eri(5,5,5)= r6(10)-r5(6,1)+r4(3,5)+r4(10,5)-r3(1,5)-r3(6,5)+r2(5,13)-r1(1,9)+( &
&                +r5(10,2)+r5(10,3)-r4(6,3)-r4(6,4)+r3(3,9)+r3(3,10)-r2(1,9)-r2(1,10))*qx+( &
&                +r5(6,2)+r5(6,3)-r4(3,3)-r4(3,4)+r3(6,9)+r3(6,10)-r2(5,9)-r2(5,10))*qz+( &
&                +r4(10,7)-r3(6,7)+r2(5,15)-r1(1,11))*xx+(+r4(6,6)+r4(6,7)*two+r4(6,8) &
&                -r3(3,6)-r3(3,7)*two-r3(3,8))*xz+(+r4(3,7)-r3(1,7)+r2(5,15)-r1(1,11))*zz+( &
&                +r3(6,11)+r3(6,12)-r2(5,11)-r2(5,12))*xxz+(+r3(3,11)+r3(3,12)-r2(1,11) &
&                -r2(1,12))*xzz+(+r2(5,17)-r1(1,13))*xxzz
      eri(6,5,5)= r6(14)-r5(9,1)+r4(5,5)-r3(2,5)+(+r5(14,3)-r4(9,4)+r3(5,10)-r2(4,10))*qx+( &
&                +r5(9,2)+r5(9,3)-r4(5,3)-r4(5,4))*qz+(+r4(9,7)+r4(9,8)-r3(5,7)-r3(5,8))*xz+( &
&                +r4(5,7)-r3(2,7))*zz+(+r3(5,12)-r2(4,12))*xzz
      eri(1,6,5)= r6(9)-r5(5,1)+r4(9,5)-r3(5,5)+(+r5(9,2)*two-r4(5,3)*two)*qx+(+r5(5,3) &
&                -r4(2,4)+r3(5,10)-r2(4,10))*qz+(+r4(9,6)-r3(5,6))*xx+(+r4(5,7)*two &
&                -r3(2,7)*two)*xz+(+r3(5,11)-r2(4,11))*xxz
      eri(2,6,5)= r6(18)-r5(12,1)+r4(9,5)*three-r3(5,5)*three+(+r5(12,3)-r4(7,4) &
&                +r3(5,10)*three-r2(4,10)*three)*qz
      eri(3,6,5)= r6(20)-r5(14,1)+r4(9,5)*three-r3(5,5)*three+(+r5(14,2)*two+r5(14,3) &
&                -r4(9,3)*two-r4(9,4)+r3(5,9)*two+r3(5,10)-r2(4,9)*two-r2(4,10))*qz+(+r4(9,6) &
&                +r4(9,7)*two-r3(5,6)-r3(5,7)*two)*zz+(+r3(5,11)-r2(4,11))*zzz
      eri(4,6,5)= r6(13)-r5(8,1)+r4(6,5)-r3(3,5)+(+r5(13,2)-r4(8,3)+r3(6,9)-r2(5,9))*qx+( &
&                +r5(8,3)-r4(4,4)+r3(3,10)-r2(1,10))*qz+(+r4(8,7)-r3(4,7)+r2(5,15)-r1(1,11)) &
&                *xz
      eri(5,6,5)= r6(14)-r5(9,1)+r4(5,5)-r3(2,5)+(+r5(14,2)-r4(9,3)+r3(5,9)-r2(4,9))*qx+( &
&                +r5(9,2)+r5(9,3)-r4(5,3)-r4(5,4))*qz+(+r4(9,6)+r4(9,7)-r3(5,6)-r3(5,7))*xz+( &
&                +r4(5,7)-r3(2,7))*zz+(+r3(5,11)-r2(4,11))*xzz
      eri(6,6,5)= r6(19)-r5(13,1)+r4(8,5)+r4(10,5)-r3(4,5)-r3(6,5)+r2(5,13)-r1(1,9)+( &
&                +r5(13,2)+r5(13,3)-r4(8,3)-r4(8,4)+r3(6,9)+r3(6,10)-r2(5,9)-r2(5,10))*qz+( &
&                +r4(8,7)-r3(4,7)+r2(5,15)-r1(1,11))*zz
      eri(1,1,6)= r6(5)-r5(2,1)+r4(5,5)*six-r3(2,5)*six+r2(6,13)*three-r1(2,9)*three+( &
&                +r5(5,2)*two+r5(5,3)*two-r4(2,3)*two-r4(2,4)*two+r3(5,9)*six+r3(5,10)*six &
&                -r2(4,9)*six-r2(4,10)*six)*qx+(+r4(5,6)+r4(5,7)*four+r4(5,8)-r3(2,6) &
&                -r3(2,7)*four-r3(2,8)+r2(6,14)+r2(6,15)*four+r2(6,16)-r1(2,10)-r1(2,11)*four &
&                -r1(2,12))*xx+(+r3(5,11)*two+r3(5,12)*two-r2(4,11)*two-r2(4,12)*two)*xxx+( &
&                +r2(6,17)-r1(2,13))*xxxx
      eri(2,1,6)= r6(12)-r5(7,1)+r4(5,5)+r4(12,5)-r3(2,5)-r3(7,5)+r2(6,13)-r1(2,9)+( &
&                +r5(12,3)*two-r4(7,4)*two+r3(5,10)*two-r2(4,10)*two)*qx+(+r4(12,8)-r3(7,8) &
&                +r2(6,16)-r1(2,12))*xx
      eri(3,1,6)= r6(14)-r5(9,1)+r4(5,5)+r4(14,5)-r3(2,5)-r3(9,5)+r2(6,13)-r1(2,9)+( &
&                +r5(14,3)*two-r4(9,4)*two+r3(5,10)*two-r2(4,10)*two)*qx+(+r5(9,2)*two &
&                -r4(5,3)*two+r3(9,9)*two-r2(6,9)*two)*qz+(+r4(14,8)-r3(9,8)+r2(6,16) &
&                -r1(2,12))*xx+(+r4(9,7)*four-r3(5,7)*four)*xz+(+r4(5,6)-r3(2,6)+r2(6,14) &
&                -r1(2,10))*zz+(+r3(9,12)*two-r2(6,12)*two)*xxz+(+r3(5,11)*two-r2(4,11)*two) &
&                *xzz+(+r2(6,17)-r1(2,13))*xxzz
      eri(4,1,6)= r6(8)-r5(4,1)+r4(8,5)*three-r3(4,5)*three+(+r5(8,2)+r5(8,3)*two-r4(4,3) &
&                -r4(4,4)*two+r3(8,9)+r3(8,10)*two-r2(2,9)-r2(2,10)*two)*qx+(+r4(8,7)*two &
&                +r4(8,8)-r3(4,7)*two-r3(4,8))*xx+(+r3(8,12)-r2(2,12))*xxx
      eri(5,1,6)= r6(9)-r5(5,1)+r4(9,5)*three-r3(5,5)*three+(+r5(9,2)+r5(9,3)*two-r4(5,3) &
&                -r4(5,4)*two+r3(9,9)+r3(9,10)*two-r2(6,9)-r2(6,10)*two)*qx+(+r5(5,2)-r4(2,3) &
&                +r3(5,9)*three-r2(4,9)*three)*qz+(+r4(9,7)*two+r4(9,8)-r3(5,7)*two-r3(5,8)) &
&                *xx+(+r4(5,6)+r4(5,7)*two-r3(2,6)-r3(2,7)*two+r2(6,14)+r2(6,15)*two-r1(2,10) &
&                -r1(2,11)*two)*xz+(+r3(9,12)-r2(6,12))*xxx+(+r3(5,11)*two+r3(5,12) &
&                -r2(4,11)*two-r2(4,12))*xxz+(+r2(6,17)-r1(2,13))*xxxz
      eri(6,1,6)= r6(13)-r5(8,1)+r4(13,5)-r3(8,5)+(+r5(13,3)*two-r4(8,4)*two)*qx+(+r5(8,2) &
&                -r4(4,3)+r3(8,9)-r2(2,9))*qz+(+r4(13,8)-r3(8,8))*xx+(+r4(8,7)*two &
&                -r3(4,7)*two)*xz+(+r3(8,12)-r2(2,12))*xxz
      eri(1,2,6)= r6(12)-r5(7,1)+r4(5,5)+r4(12,5)-r3(2,5)-r3(7,5)+r2(6,13)-r1(2,9)+( &
&                +r5(12,2)*two-r4(7,3)*two+r3(5,9)*two-r2(4,9)*two)*qx+(+r4(12,6)-r3(7,6) &
&                +r2(6,14)-r1(2,10))*xx
      eri(2,2,6)= r6(23)-r5(16,1)+r4(12,5)*six-r3(7,5)*six+r2(6,13)*three-r1(2,9)*three
      eri(3,2,6)= r6(25)-r5(18,1)+r4(12,5)+r4(14,5)-r3(7,5)-r3(9,5)+r2(6,13)-r1(2,9)+( &
&                +r5(18,2)*two-r4(12,3)*two+r3(9,9)*two-r2(6,9)*two)*qz+(+r4(12,6)-r3(7,6) &
&                +r2(6,14)-r1(2,10))*zz
      eri(4,2,6)= r6(17)-r5(11,1)+r4(8,5)*three-r3(4,5)*three+(+r5(17,2)-r4(11,3) &
&                +r3(8,9)*three-r2(2,9)*three)*qx
      eri(5,2,6)= r6(18)-r5(12,1)+r4(9,5)-r3(5,5)+(+r5(18,2)-r4(12,3)+r3(9,9)-r2(6,9))*qx+( &
&                +r5(12,2)-r4(7,3)+r3(5,9)-r2(4,9))*qz+(+r4(12,6)-r3(7,6)+r2(6,14)-r1(2,10)) &
&                *xz
      eri(6,2,6)= r6(24)-r5(17,1)+r4(13,5)*three-r3(8,5)*three+(+r5(17,2)-r4(11,3) &
&                +r3(8,9)*three-r2(2,9)*three)*qz
      eri(1,3,6)= r6(14)-r5(9,1)+r4(5,5)+r4(14,5)-r3(2,5)-r3(9,5)+r2(6,13)-r1(2,9)+( &
&                +r5(14,2)*two-r4(9,3)*two+r3(5,9)*two-r2(4,9)*two)*qx+(+r5(9,3)*two &
&                -r4(5,4)*two+r3(9,10)*two-r2(6,10)*two)*qz+(+r4(14,6)-r3(9,6)+r2(6,14) &
&                -r1(2,10))*xx+(+r4(9,7)*four-r3(5,7)*four)*xz+(+r4(5,8)-r3(2,8)+r2(6,16) &
&                -r1(2,12))*zz+(+r3(9,11)*two-r2(6,11)*two)*xxz+(+r3(5,12)*two-r2(4,12)*two) &
&                *xzz+(+r2(6,17)-r1(2,13))*xxzz
      eri(2,3,6)= r6(25)-r5(18,1)+r4(12,5)+r4(14,5)-r3(7,5)-r3(9,5)+r2(6,13)-r1(2,9)+( &
&                +r5(18,3)*two-r4(12,4)*two+r3(9,10)*two-r2(6,10)*two)*qz+(+r4(12,8)-r3(7,8) &
&                +r2(6,16)-r1(2,12))*zz
      eri(3,3,6)= r6(27)-r5(20,1)+r4(14,5)*six-r3(9,5)*six+r2(6,13)*three-r1(2,9)*three+( &
&                +r5(20,2)*two+r5(20,3)*two-r4(14,3)*two-r4(14,4)*two+r3(9,9)*six &
&                +r3(9,10)*six-r2(6,9)*six-r2(6,10)*six)*qz+(+r4(14,6)+r4(14,7)*four+r4(14,8) &
&                -r3(9,6)-r3(9,7)*four-r3(9,8)+r2(6,14)+r2(6,15)*four+r2(6,16)-r1(2,10) &
&                -r1(2,11)*four-r1(2,12))*zz+(+r3(9,11)*two+r3(9,12)*two-r2(6,11)*two &
&                -r2(6,12)*two)*zzz+(+r2(6,17)-r1(2,13))*zzzz
      eri(4,3,6)= r6(19)-r5(13,1)+r4(8,5)-r3(4,5)+(+r5(19,2)-r4(13,3)+r3(8,9)-r2(2,9))*qx+( &
&                +r5(13,3)*two-r4(8,4)*two)*qz+(+r4(13,7)*two-r3(8,7)*two)*xz+(+r4(8,8) &
&                -r3(4,8))*zz+(+r3(8,12)-r2(2,12))*xzz
      eri(5,3,6)= r6(20)-r5(14,1)+r4(9,5)*three-r3(5,5)*three+(+r5(20,2)-r4(14,3) &
&                +r3(9,9)*three-r2(6,9)*three)*qx+(+r5(14,2)+r5(14,3)*two-r4(9,3)-r4(9,4)*two &
&                +r3(5,9)+r3(5,10)*two-r2(4,9)-r2(4,10)*two)*qz+(+r4(14,6)+r4(14,7)*two &
&                -r3(9,6)-r3(9,7)*two+r2(6,14)+r2(6,15)*two-r1(2,10)-r1(2,11)*two)*xz+( &
&                +r4(9,7)*two+r4(9,8)-r3(5,7)*two-r3(5,8))*zz+(+r3(9,11)*two+r3(9,12) &
&                -r2(6,11)*two-r2(6,12))*xzz+(+r3(5,12)-r2(4,12))*zzz+(+r2(6,17)-r1(2,13)) &
&                *xzzz
      eri(6,3,6)= r6(26)-r5(19,1)+r4(13,5)*three-r3(8,5)*three+(+r5(19,2)+r5(19,3)*two &
&                -r4(13,3)-r4(13,4)*two+r3(8,9)+r3(8,10)*two-r2(2,9)-r2(2,10)*two)*qz+( &
&                +r4(13,7)*two+r4(13,8)-r3(8,7)*two-r3(8,8))*zz+(+r3(8,12)-r2(2,12))*zzz
      eri(1,4,6)= r6(8)-r5(4,1)+r4(8,5)*three-r3(4,5)*three+(+r5(8,2)*two+r5(8,3) &
&                -r4(4,3)*two-r4(4,4)+r3(8,9)*two+r3(8,10)-r2(2,9)*two-r2(2,10))*qx+(+r4(8,6) &
&                +r4(8,7)*two-r3(4,6)-r3(4,7)*two)*xx+(+r3(8,11)-r2(2,11))*xxx
      eri(2,4,6)= r6(17)-r5(11,1)+r4(8,5)*three-r3(4,5)*three+(+r5(17,3)-r4(11,4) &
&                +r3(8,10)*three-r2(2,10)*three)*qx
      eri(3,4,6)= r6(19)-r5(13,1)+r4(8,5)-r3(4,5)+(+r5(19,3)-r4(13,4)+r3(8,10)-r2(2,10))*qx &
&                +(+r5(13,2)*two-r4(8,3)*two)*qz+(+r4(13,7)*two-r3(8,7)*two)*xz+(+r4(8,6) &
&                -r3(4,6))*zz+(+r3(8,11)-r2(2,11))*xzz
      eri(4,4,6)= r6(12)-r5(7,1)+r4(5,5)+r4(12,5)-r3(2,5)-r3(7,5)+r2(6,13)-r1(2,9)+( &
&                +r5(12,2)+r5(12,3)-r4(7,3)-r4(7,4)+r3(5,9)+r3(5,10)-r2(4,9)-r2(4,10))*qx+( &
&                +r4(12,7)-r3(7,7)+r2(6,15)-r1(2,11))*xx
      eri(5,4,6)= r6(13)-r5(8,1)+r4(13,5)-r3(8,5)+(+r5(13,2)+r5(13,3)-r4(8,3)-r4(8,4))*qx+( &
&                +r5(8,2)-r4(4,3)+r3(8,9)-r2(2,9))*qz+(+r4(13,7)-r3(8,7))*xx+(+r4(8,6) &
&                +r4(8,7)-r3(4,6)-r3(4,7))*xz+(+r3(8,11)-r2(2,11))*xxz
      eri(6,4,6)= r6(18)-r5(12,1)+r4(9,5)-r3(5,5)+(+r5(18,3)-r4(12,4)+r3(9,10)-r2(6,10))*qx &
&                +(+r5(12,2)-r4(7,3)+r3(5,9)-r2(4,9))*qz+(+r4(12,7)-r3(7,7)+r2(6,15)-r1(2,11) &
&                )*xz
      eri(1,5,6)= r6(9)-r5(5,1)+r4(9,5)*three-r3(5,5)*three+(+r5(9,2)*two+r5(9,3) &
&                -r4(5,3)*two-r4(5,4)+r3(9,9)*two+r3(9,10)-r2(6,9)*two-r2(6,10))*qx+(+r5(5,3) &
&                -r4(2,4)+r3(5,10)*three-r2(4,10)*three)*qz+(+r4(9,6)+r4(9,7)*two-r3(5,6) &
&                -r3(5,7)*two)*xx+(+r4(5,7)*two+r4(5,8)-r3(2,7)*two-r3(2,8)+r2(6,15)*two &
&                +r2(6,16)-r1(2,11)*two-r1(2,12))*xz+(+r3(9,11)-r2(6,11))*xxx+(+r3(5,11) &
&                +r3(5,12)*two-r2(4,11)-r2(4,12)*two)*xxz+(+r2(6,17)-r1(2,13))*xxxz
      eri(2,5,6)= r6(18)-r5(12,1)+r4(9,5)-r3(5,5)+(+r5(18,3)-r4(12,4)+r3(9,10)-r2(6,10))*qx &
&                +(+r5(12,3)-r4(7,4)+r3(5,10)-r2(4,10))*qz+(+r4(12,8)-r3(7,8)+r2(6,16) &
&                -r1(2,12))*xz
      eri(3,5,6)= r6(20)-r5(14,1)+r4(9,5)*three-r3(5,5)*three+(+r5(20,3)-r4(14,4) &
&                +r3(9,10)*three-r2(6,10)*three)*qx+(+r5(14,2)*two+r5(14,3)-r4(9,3)*two &
&                -r4(9,4)+r3(5,9)*two+r3(5,10)-r2(4,9)*two-r2(4,10))*qz+(+r4(14,7)*two &
&                +r4(14,8)-r3(9,7)*two-r3(9,8)+r2(6,15)*two+r2(6,16)-r1(2,11)*two-r1(2,12)) &
&                *xz+(+r4(9,6)+r4(9,7)*two-r3(5,6)-r3(5,7)*two)*zz+(+r3(9,11)+r3(9,12)*two &
&                -r2(6,11)-r2(6,12)*two)*xzz+(+r3(5,11)-r2(4,11))*zzz+(+r2(6,17)-r1(2,13)) &
&                *xzzz
      eri(4,5,6)= r6(13)-r5(8,1)+r4(13,5)-r3(8,5)+(+r5(13,2)+r5(13,3)-r4(8,3)-r4(8,4))*qx+( &
&                +r5(8,3)-r4(4,4)+r3(8,10)-r2(2,10))*qz+(+r4(13,7)-r3(8,7))*xx+(+r4(8,7) &
&                +r4(8,8)-r3(4,7)-r3(4,8))*xz+(+r3(8,12)-r2(2,12))*xxz
      eri(5,5,6)= r6(14)-r5(9,1)+r4(5,5)+r4(14,5)-r3(2,5)-r3(9,5)+r2(6,13)-r1(2,9)+( &
&                +r5(14,2)+r5(14,3)-r4(9,3)-r4(9,4)+r3(5,9)+r3(5,10)-r2(4,9)-r2(4,10))*qx+( &
&                +r5(9,2)+r5(9,3)-r4(5,3)-r4(5,4)+r3(9,9)+r3(9,10)-r2(6,9)-r2(6,10))*qz+( &
&                +r4(14,7)-r3(9,7)+r2(6,15)-r1(2,11))*xx+(+r4(9,6)+r4(9,7)*two+r4(9,8) &
&                -r3(5,6)-r3(5,7)*two-r3(5,8))*xz+(+r4(5,7)-r3(2,7)+r2(6,15)-r1(2,11))*zz+( &
&                +r3(9,11)+r3(9,12)-r2(6,11)-r2(6,12))*xxz+(+r3(5,11)+r3(5,12)-r2(4,11) &
&                -r2(4,12))*xzz+(+r2(6,17)-r1(2,13))*xxzz
      eri(6,5,6)= r6(19)-r5(13,1)+r4(8,5)-r3(4,5)+(+r5(19,3)-r4(13,4)+r3(8,10)-r2(2,10))*qx &
&                +(+r5(13,2)+r5(13,3)-r4(8,3)-r4(8,4))*qz+(+r4(13,7)+r4(13,8)-r3(8,7)-r3(8,8) &
&                )*xz+(+r4(8,7)-r3(4,7))*zz+(+r3(8,12)-r2(2,12))*xzz
      eri(1,6,6)= r6(13)-r5(8,1)+r4(13,5)-r3(8,5)+(+r5(13,2)*two-r4(8,3)*two)*qx+(+r5(8,3) &
&                -r4(4,4)+r3(8,10)-r2(2,10))*qz+(+r4(13,6)-r3(8,6))*xx+(+r4(8,7)*two &
&                -r3(4,7)*two)*xz+(+r3(8,11)-r2(2,11))*xxz
      eri(2,6,6)= r6(24)-r5(17,1)+r4(13,5)*three-r3(8,5)*three+(+r5(17,3)-r4(11,4) &
&                +r3(8,10)*three-r2(2,10)*three)*qz
      eri(3,6,6)= r6(26)-r5(19,1)+r4(13,5)*three-r3(8,5)*three+(+r5(19,2)*two+r5(19,3) &
&                -r4(13,3)*two-r4(13,4)+r3(8,9)*two+r3(8,10)-r2(2,9)*two-r2(2,10))*qz+( &
&                +r4(13,6)+r4(13,7)*two-r3(8,6)-r3(8,7)*two)*zz+(+r3(8,11)-r2(2,11))*zzz
      eri(4,6,6)= r6(18)-r5(12,1)+r4(9,5)-r3(5,5)+(+r5(18,2)-r4(12,3)+r3(9,9)-r2(6,9))*qx+( &
&                +r5(12,3)-r4(7,4)+r3(5,10)-r2(4,10))*qz+(+r4(12,7)-r3(7,7)+r2(6,15)-r1(2,11) &
&                )*xz
      eri(5,6,6)= r6(19)-r5(13,1)+r4(8,5)-r3(4,5)+(+r5(19,2)-r4(13,3)+r3(8,9)-r2(2,9))*qx+( &
&                +r5(13,2)+r5(13,3)-r4(8,3)-r4(8,4))*qz+(+r4(13,6)+r4(13,7)-r3(8,6)-r3(8,7)) &
&                *xz+(+r4(8,7)-r3(4,7))*zz+(+r3(8,11)-r2(2,11))*xzz
      eri(6,6,6)= r6(25)-r5(18,1)+r4(12,5)+r4(14,5)-r3(7,5)-r3(9,5)+r2(6,13)-r1(2,9)+( &
&                +r5(18,2)+r5(18,3)-r4(12,3)-r4(12,4)+r3(9,9)+r3(9,10)-r2(6,9)-r2(6,10))*qz+( &
&                +r4(12,7)-r3(7,7)+r2(6,15)-r1(2,11))*zz
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
        do k= 1,6
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
        do k= 1,6
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
      if(nbfijkl(3) == 6)then
        do j= 1,nbfijkl(2)
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
        do j= 1,nbfijkl(2)
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
        do j= 1,nbfijkl(2)
          do k= 1,nbfijkl(3)
            do l= 1,6
              phmdint(l,k,j,1)= eri(1,k,j)*rot2(1,l)+eri(2,k,j)*rot2(2,l)+eri(3,k,j)*rot2(3,l) &
&                              +eri(4,k,j)*rot2(4,l)+eri(5,k,j)*rot2(5,l)+eri(6,k,j)*rot2(6,l) 
            enddo
          enddo
        enddo
      else
        do j= 1,nbfijkl(2)
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
