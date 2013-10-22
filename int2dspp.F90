!--------------------------------------------------------------------------
  subroutine int2dspp(phmdint,exfac1,exfac2,xyziq,xzkl,rot,nijkl,nbfijkl)
!--------------------------------------------------------------------------
!
! Calculate (ds|pp) integrals
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
      real(8),parameter :: six=6.0D+00, seven=7.0D+00, sqrt3=1.73205080756888D+00
      real(8),parameter :: sqrt3h=8.660254037844386D-01
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*), xzkl(2), rot(3,3)
      real(8),intent(out) :: phmdint(6,6,6,6)
      real(8) :: tval, tval1, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:4), ftw(4,5)
      real(8) :: f0(2), f1(2,5), f2(3,5), f3(4,3), f4(5), r0(4), r1(3,6), r2(6,6), r3(10,3), r4(15)
      real(8) :: ex12, ex34, ex43, ex41, expq, expq2, ex3q, ex4q, c12, c34, zip, xiq, yiq, ziq
      real(8) :: xiq2, yiq2, xyiq, xiq4, yiq4, xyiq2, xypq2, zpq, zpq2, zpq3, fac, ex33q, ex3qmd
      real(8) :: zjp, pmd, qmd, qmd2, qmd2x, qmd2y, qmd2xy, qx, qz, xx, xz, zz
      real(8) :: eri(6,3,3), work(6), f1w(3,4), f2w(6,5), f3w(10,3), rot2(6,6), rot3(6,5)
!
! Zero-clear
!
      r0(1:4)     = zero
      r1(1:3,1:6) = zero
      r2(1:6,1:6) = zero
      r3(1:10,1:3)= zero
      r4(1:15)    = zero
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
        f1(1:2,1:5)= zero
        f2(1:3,1:5)= zero
        f3(1:4,1:3)= zero
        f4(1:5)= zero
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
            do ii= 0,4
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
            fac= c12*sqrt(ex41)
            expq= expq*two
            expq2= expq*expq
            ft(0)= ft(0)*fac
            ft(1)= ft(1)*fac*expq
            ft(2)= ft(2)*fac*expq2
            ft(3)= ft(3)*fac*expq2*expq
            ft(4)= ft(4)*fac*expq2*expq2
          endif
          zpq3= zpq2*zpq
          ftw(1,1)= zjp*zip
          ftw(1,2)= pmd
          ftw(1,3)= pmd*zjp
          ftw(1,4)= pmd*zip
          ftw(1,5)= pmd*pmd
          do i= 1,5
            ftw(2,i)= ftw(1,i)*zpq
            ftw(3,i)= ftw(1,i)*zpq2
            ftw(4,i)= ftw(1,i)*zpq3
          enddo
          f0(1)= f0(1)+ft(0)*ftw(1,1)
          f0(2)= f0(2)+ft(0)*ftw(1,2)
          do i= 1,4
            f1(1,i)= f1(1,i)-ft(1)*ftw(1,i)
            f1(2,i)= f1(2,i)-ft(1)*ftw(2,i)
          enddo
          f1(1,5)= f1(1,5)-ft(1)*ftw(1,5)
          do i= 1,5
            f2(1,i)= f2(1,i)+ft(2)*ftw(1,i)
            f2(2,i)= f2(2,i)+ft(2)*ftw(2,i)
            f2(3,i)= f2(3,i)+ft(2)*ftw(3,i)
          enddo
          do i= 1,3
            f3(1,i)= f3(1,i)-ft(3)*ftw(1,i+2)
            f3(2,i)= f3(2,i)-ft(3)*ftw(2,i+2)
            f3(3,i)= f3(3,i)-ft(3)*ftw(3,i+2)
            f3(4,i)= f3(4,i)-ft(3)*ftw(4,i+2)
          enddo
          f4(1)= f4(1)+ft(4)*ftw(1,5)
          f4(2)= f4(2)+ft(4)*ftw(2,5)
          f4(3)= f4(3)+ft(4)*ftw(3,5)
          f4(4)= f4(4)+ft(4)*ftw(4,5)
          f4(5)= f4(5)+ft(4)*ftw(4,5)*zpq
        enddo
!
        qmd = ex43*c34
        qmd2= qmd*ex43
        ex33q= ex3q*ex3q*c34
        qmd2x= qmd2*xiq
        qmd2y= qmd2*yiq
        qmd2xy=qmd2*xiq*yiq
        ex3qmd=ex3q*qmd
!
        r0(1)= r0(1)+f0(1)*qmd
        r0(2)= r0(2)+f0(1)*ex33q
        r0(3)= r0(3)+f0(2)*qmd
        r0(4)= r0(4)+f0(2)*ex33q
!
        do i= 1,4
          f1w(1,i)= f1(1,i)*xiq
          f1w(2,i)= f1(1,i)*yiq
          f1w(3,i)= f1(2,i)
        enddo
        do i= 1,3
          r1(i,1)= r1(i,1)+f1w(i,1)*ex3qmd
          r1(i,2)= r1(i,2)+f1w(i,2)*ex3qmd
          r1(i,3)= r1(i,3)+f1w(i,3)*qmd
          r1(i,4)= r1(i,4)+f1w(i,3)*ex33q
          r1(i,5)= r1(i,5)+f1w(i,4)*qmd
          r1(i,6)= r1(i,6)+f1w(i,4)*ex33q
        enddo
!
        do i= 1,5
          f2w(1,i)= f2(1,i)*xiq2+f1(1,i)
          f2w(2,i)= f2(1,i)*yiq2+f1(1,i)
          f2w(3,i)= f2(3,i)     +f1(1,i)
          f2w(4,i)= f2(1,i)*xyiq
          f2w(5,i)= f2(2,i)*xiq
          f2w(6,i)= f2(2,i)*yiq
        enddo
        do i= 1,6
          r2(i,1)= r2(i,1)+f2w(i,1)*qmd2
          r2(i,2)= r2(i,2)+f2w(i,2)*qmd2
          r2(i,3)= r2(i,3)+f2w(i,3)*ex3qmd
          r2(i,4)= r2(i,4)+f2w(i,4)*ex3qmd
          r2(i,5)= r2(i,5)+f2w(i,5)*qmd
          r2(i,6)= r2(i,6)+f2w(i,5)*ex33q
        enddo
! 
        do i= 1,3
          f3w( 1,i)=(f3(1,i)*xiq2+f2(1,i+2)*three)*xiq
          f3w( 2,i)=(f3(1,i)*xiq2+f2(1,i+2)      )*yiq
          f3w( 3,i)=(f3(2,i)*xiq2+f2(2,i+2)      )
          f3w( 4,i)=(f3(1,i)*yiq2+f2(1,i+2)      )*xiq
          f3w( 5,i)=(f3(2,i)*xyiq                )
          f3w( 6,i)=(f3(3,i)     +f2(1,i+2)      )*xiq
          f3w( 7,i)=(f3(1,i)*yiq2+f2(1,i+2)*three)*yiq
          f3w( 8,i)=(f3(2,i)*yiq2+f2(2,i+2)      )
          f3w( 9,i)=(f3(3,i)     +f2(1,i+2)      )*yiq
          f3w(10,i)=(f3(4,i)     +f2(2,i+2)*three)
        enddo
        do i= 1,10
          r3(i,1)= r3(i,1)+f3w(i,1)*qmd2
          r3(i,2)= r3(i,2)+f3w(i,2)*qmd2
          r3(i,3)= r3(i,3)+f3w(i,3)*ex3qmd
        enddo
!
        r4( 1)= r4( 1)+(f4(1)*xiq4 +f3(1,3)*xiq2*six         +f2(1,5)*three)*qmd2
        r4( 2)= r4( 2)+(f4(1)*xiq2 +f3(1,3)*three                          )*qmd2xy
        r4( 3)= r4( 3)+(f4(2)*xiq2 +f3(2,3)*three                          )*qmd2x
        r4( 4)= r4( 4)+(f4(1)*xyiq2+f3(1,3)*xiq2+f3(1,3)*yiq2+f2(1,5)      )*qmd2
        r4( 5)= r4( 5)+(f4(2)*xiq2 +f3(2,3)                                )*qmd2y
        r4( 6)= r4( 6)+(f4(3)*xiq2 +f3(1,3)*xiq2+f3(3,3)     +f2(1,5)      )*qmd2
        r4( 7)= r4( 7)+(f4(1)*yiq2 +f3(1,3)*three                          )*qmd2xy
        r4( 8)= r4( 8)+(f4(2)*yiq2 +f3(2,3)                                )*qmd2x
        r4( 9)= r4( 9)+(f4(3)      +f3(1,3)                                )*qmd2xy
        r4(10)= r4(10)+(f4(4)      +f3(2,3)*three                          )*qmd2x
        r4(11)= r4(11)+(f4(1)*yiq4 +f3(1,3)*yiq2*six         +f2(1,5)*three)*qmd2
        r4(12)= r4(12)+(f4(2)*yiq2 +f3(2,3)*three                          )*qmd2y
        r4(13)= r4(13)+(f4(3)*yiq2 +f3(1,3)*yiq2+f3(3,3)     +f2(1,5)      )*qmd2
        r4(14)= r4(14)+(f4(4)      +f3(2,3)*three                          )*qmd2y
        r4(15)= r4(15)+(f4(5)      +f3(3,3)*six              +f2(1,5)*three)*qmd2
      enddo
!
      qx= xzkl(1)
      qz= xzkl(2)
      xx= qx*qx
      xz= qx*qz
      zz= qz*qz
      eri(1,1,1)= r4(1)+r2(1,2)+r2(1,5)+r0(3)+(+r3(1,3)*two+r1(1,2)*two)*qx+(+r2(1,6) &
&                +r0(4))*xx
      eri(2,1,1)= r4(4)+r2(2,2)+r2(1,5)+r0(3)
      eri(3,1,1)= r4(6)+r2(3,2)+r2(1,5)+r0(3)+(+r3(3,3)*two+r1(3,2)*two)*qz+(+r2(1,6) &
&                +r0(4))*zz
      eri(4,1,1)= r4(2)+r2(4,2)+(+r3(2,3)+r1(2,2))*qx
      eri(5,1,1)= r4(3)+r2(5,2)+(+r3(3,3)+r1(3,2))*qx+(+r3(1,3)+r1(1,2))*qz+(+r2(1,6) &
&                +r0(4))*xz
      eri(6,1,1)= r4(5)+r2(6,2)+(+r3(2,3)+r1(2,2))*qz
      eri(1,2,1)= r4(2)+r2(4,5)+r3(2,3)*two*qx+r2(4,6)*xx
      eri(2,2,1)= r4(7)+r2(4,5)
      eri(3,2,1)= r4(9)+r2(4,5)+r3(5,3)*two*qz+r2(4,6)*zz
      eri(4,2,1)= r4(4)+r3(4,3)*qx
      eri(5,2,1)= r4(5)+r3(5,3)*qx+r3(2,3)*qz+r2(4,6)*xz
      eri(6,2,1)= r4(8)+r3(4,3)*qz
      eri(1,3,1)= r4(3)-r3(1,1)+r2(5,5)-r1(1,3)+(+r3(3,3)*two-r2(1,3)*two)*qx+(+r2(5,6) &
&                -r1(1,4))*xx
      eri(2,3,1)= r4(8)-r3(4,1)+r2(5,5)-r1(1,3)
      eri(3,3,1)= r4(10)-r3(6,1)+r2(5,5)-r1(1,3)+(+r3(6,3)*two-r2(5,3)*two)*qz+(+r2(5,6) &
&                -r1(1,4))*zz
      eri(4,3,1)= r4(5)-r3(2,1)+(+r3(5,3)-r2(4,3))*qx
      eri(5,3,1)= r4(6)-r3(3,1)+(+r3(6,3)-r2(5,3))*qx+(+r3(3,3)-r2(1,3))*qz+(+r2(5,6) &
&                -r1(1,4))*xz
      eri(6,3,1)= r4(9)-r3(5,1)+(+r3(5,3)-r2(4,3))*qz
      eri(1,1,2)= eri(1,2,1)
      eri(2,1,2)= eri(2,2,1)
      eri(3,1,2)= eri(3,2,1)
      eri(4,1,2)= eri(4,2,1)
      eri(5,1,2)= eri(5,2,1)
      eri(6,1,2)= eri(6,2,1)
      eri(1,2,2)= r4(4)+r2(1,2)+r2(2,5)+r0(3)+(+r3(4,3)*two+r1(1,2)*two)*qx+(+r2(2,6) &
&                +r0(4))*xx
      eri(2,2,2)= r4(11)+r2(2,2)+r2(2,5)+r0(3)
      eri(3,2,2)= r4(13)+r2(3,2)+r2(2,5)+r0(3)+(+r3(8,3)*two+r1(3,2)*two)*qz+(+r2(2,6) &
&                +r0(4))*zz
      eri(4,2,2)= r4(7)+r2(4,2)+(+r3(7,3)+r1(2,2))*qx
      eri(5,2,2)= r4(8)+r2(5,2)+(+r3(8,3)+r1(3,2))*qx+(+r3(4,3)+r1(1,2))*qz+(+r2(2,6) &
&                +r0(4))*xz
      eri(6,2,2)= r4(12)+r2(6,2)+(+r3(7,3)+r1(2,2))*qz
      eri(1,3,2)= r4(5)-r3(2,1)+r2(6,5)-r1(2,3)+(+r3(5,3)*two-r2(4,3)*two)*qx+(+r2(6,6) &
&                -r1(2,4))*xx
      eri(2,3,2)= r4(12)-r3(7,1)+r2(6,5)-r1(2,3)
      eri(3,3,2)= r4(14)-r3(9,1)+r2(6,5)-r1(2,3)+(+r3(9,3)*two-r2(6,3)*two)*qz+(+r2(6,6) &
&                -r1(2,4))*zz
      eri(4,3,2)= r4(8)-r3(4,1)+(+r3(8,3)-r2(2,3))*qx
      eri(5,3,2)= r4(9)-r3(5,1)+(+r3(9,3)-r2(6,3))*qx+(+r3(5,3)-r2(4,3))*qz+(+r2(6,6) &
&                -r1(2,4))*xz
      eri(6,3,2)= r4(13)-r3(8,1)+(+r3(8,3)-r2(2,3))*qz
      eri(1,1,3)= r4(3)-r3(1,2)+r2(5,5)-r1(1,5)+(+r3(3,3)*two-r2(1,4)*two)*qx+(+r2(5,6) &
&                -r1(1,6))*xx
      eri(2,1,3)= r4(8)-r3(4,2)+r2(5,5)-r1(1,5)
      eri(3,1,3)= r4(10)-r3(6,2)+r2(5,5)-r1(1,5)+(+r3(6,3)*two-r2(5,4)*two)*qz+(+r2(5,6) &
&                -r1(1,6))*zz
      eri(4,1,3)= r4(5)-r3(2,2)+(+r3(5,3)-r2(4,4))*qx
      eri(5,1,3)= r4(6)-r3(3,2)+(+r3(6,3)-r2(5,4))*qx+(+r3(3,3)-r2(1,4))*qz+(+r2(5,6) &
&                -r1(1,6))*xz
      eri(6,1,3)= r4(9)-r3(5,2)+(+r3(5,3)-r2(4,4))*qz
      eri(1,2,3)= r4(5)-r3(2,2)+r2(6,5)-r1(2,5)+(+r3(5,3)*two-r2(4,4)*two)*qx+(+r2(6,6) &
&                -r1(2,6))*xx
      eri(2,2,3)= r4(12)-r3(7,2)+r2(6,5)-r1(2,5)
      eri(3,2,3)= r4(14)-r3(9,2)+r2(6,5)-r1(2,5)+(+r3(9,3)*two-r2(6,4)*two)*qz+(+r2(6,6) &
&                -r1(2,6))*zz
      eri(4,2,3)= r4(8)-r3(4,2)+(+r3(8,3)-r2(2,4))*qx
      eri(5,2,3)= r4(9)-r3(5,2)+(+r3(9,3)-r2(6,4))*qx+(+r3(5,3)-r2(4,4))*qz+(+r2(6,6) &
&                -r1(2,6))*xz
      eri(6,2,3)= r4(13)-r3(8,2)+(+r3(8,3)-r2(2,4))*qz
      eri(1,3,3)= r4(6)-r3(3,1)-r3(3,2)+r2(1,1)+r2(1,2)+r2(3,5)-r1(3,3)-r1(3,5)+r0(1) &
&                +r0(3)+(+r3(6,3)*two-r2(5,3)*two-r2(5,4)*two+r1(1,1)*two+r1(1,2)*two)*qx+( &
&                +r2(3,6)-r1(3,4)-r1(3,6)+r0(2)+r0(4))*xx
      eri(2,3,3)= r4(13)-r3(8,1)-r3(8,2)+r2(2,1)+r2(2,2)+r2(3,5)-r1(3,3)-r1(3,5)+r0(1)+r0(3)
      eri(3,3,3)= r4(15)-r3(10,1)-r3(10,2)+r2(3,1)+r2(3,2)+r2(3,5)-r1(3,3)-r1(3,5)+r0(1) &
&                +r0(3)+(+r3(10,3)*two-r2(3,3)*two-r2(3,4)*two+r1(3,1)*two+r1(3,2)*two)*qz+( &
&                +r2(3,6)-r1(3,4)-r1(3,6)+r0(2)+r0(4))*zz
      eri(4,3,3)= r4(9)-r3(5,1)-r3(5,2)+r2(4,1)+r2(4,2)+(+r3(9,3)-r2(6,3)-r2(6,4)+r1(2,1) &
&                +r1(2,2))*qx
      eri(5,3,3)= r4(10)-r3(6,1)-r3(6,2)+r2(5,1)+r2(5,2)+(+r3(10,3)-r2(3,3)-r2(3,4) &
&                +r1(3,1)+r1(3,2))*qx+(+r3(6,3)-r2(5,3)-r2(5,4)+r1(1,1)+r1(1,2))*qz+(+r2(3,6) &
&                -r1(3,4)-r1(3,6)+r0(2)+r0(4))*xz
      eri(6,3,3)= r4(14)-r3(9,1)-r3(9,2)+r2(6,1)+r2(6,2)+(+r3(9,3)-r2(6,3)-r2(6,4)+r1(2,1) &
&                +r1(2,2))*qz
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
!
      do j= 1,3
        do l= 1,6
          work(1)= eri(l,j,1)
          work(2)= eri(l,j,2)
          work(3)= eri(l,j,3)
          eri(l,j,1)= work(1)*rot(1,1)+work(2)*rot(2,1)+work(3)*rot(3,1)
          eri(l,j,2)= work(1)*rot(1,2)+work(2)*rot(2,2)+work(3)*rot(3,2)
          eri(l,j,3)= work(1)*rot(1,3)+work(2)*rot(2,3)+work(3)*rot(3,3)
        enddo
      enddo
      do i= 1,3
        do l= 1,6
          work(1)= eri(l,1,i)
          work(2)= eri(l,2,i)
          work(3)= eri(l,3,i)
          eri(l,1,i)= work(1)*rot(1,1)+work(2)*rot(2,1)+work(3)*rot(3,1)
          eri(l,2,i)= work(1)*rot(1,2)+work(2)*rot(2,2)+work(3)*rot(3,2)
          eri(l,3,i)= work(1)*rot(1,3)+work(2)*rot(2,3)+work(3)*rot(3,3)
        enddo
      enddo
      if(nbfijkl(4) == 6)then
        do i= 1,3
          do j= 1,3
            do l= 1,6
              phmdint(l,1,j,i)= eri(1,j,i)*rot2(1,l)+eri(2,j,i)*rot2(2,l)+eri(3,j,i)*rot2(3,l) &
&                              +eri(4,j,i)*rot2(4,l)+eri(5,j,i)*rot2(5,l)+eri(6,j,i)*rot2(6,l) 
            enddo
          enddo
        enddo
      else
        do l= 1,6
          rot3(l,1)= rot2(l,3)-(rot2(l,1)+rot2(l,2))*half
          rot3(l,2)= rot2(l,5)
          rot3(l,3)= rot2(l,6)
          rot3(l,4)=(rot2(l,1)-rot2(l,2))*sqrt3h
          rot3(l,5)= rot2(l,4)
        enddo
        do i= 1,3
          do j= 1,3
            do l= 1,5
              phmdint(l,1,j,i)= eri(1,j,i)*rot3(1,l)+eri(2,j,i)*rot3(2,l)+eri(3,j,i)*rot3(3,l) &
&                              +eri(4,j,i)*rot3(4,l)+eri(5,j,i)*rot3(5,l)+eri(6,j,i)*rot3(6,l)
            enddo
          enddo
        enddo
      endif
!
      return
end
