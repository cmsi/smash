!--------------------------------------------------------------------------
  subroutine int2ddss(phmdint,exfac1,exfac2,xyziq,xzkl,rot,nijkl,nbfijkl)
!--------------------------------------------------------------------------
!
! Calculate (dd|ss) integrals
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
      real(8) :: tinv, ft(0:4)
      real(8) :: f0, f1(2), f2(3), f3(4), f4(5), r0(5), r1(3,4), r2(6,4), r3(10,2), r4(15)
      real(8) :: ex12, ex34, ex43, ex41, expq, expq2, ex3q, ex4q, c12, c34, zip, xiq, yiq, ziq
      real(8) :: xiq2, yiq2, xyiq, xiq4, yiq4, xyiq2, xypq2, zpq, zpq2, fac, ex33q, ex34q, ex44q
      real(8) :: qmd, qmd2, qmd3, qmd4, qmd4x, qmd4y, qmd4xy,qx, qz, xx, xz, zz, xxx, xxz, xzz
      real(8) :: zzz, xxxx, xxxz, xxzz, xzzz, zzzz, eri(6,6), work(6), f1w(3), f2w(6), f3w(10)
      real(8) :: rot2(6,6), rot3(6,5)
!
! Zero-clear
!
      r0(1:5)     = zero
      r1(1:3,1:4) = zero
      r2(1:6,1:4) = zero
      r3(1:10,1:2)= zero
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
        f0   = zero
        f1(1:2)= zero
        f2(1:3)= zero
        f3(1:4)= zero
        f4(1:5)= zero
        do ij= 1,nijkl(1)
          ex12= exfac1(1,ij)
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
          f0= f0+ft(0)
          f1(1)= f1(1)-ft(1)
          f1(2)= f1(2)-ft(1)*zpq
          f2(1)= f2(1)+ft(2)
          f2(2)= f2(2)+ft(2)*zpq
          f2(3)= f2(3)+ft(2)*zpq2
          f3(1)= f3(1)-ft(3)
          f3(2)= f3(2)-ft(3)*zpq
          f3(3)= f3(3)-ft(3)*zpq2
          f3(4)= f3(4)-ft(3)*zpq2*zpq
          f4(1)= f4(1)+ft(4)
          f4(2)= f4(2)+ft(4)*zpq
          f4(3)= f4(3)+ft(4)*zpq2
          f4(4)= f4(4)+ft(4)*zpq2*zpq
          f4(5)= f4(5)+ft(4)*zpq2*zpq2
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
        r0(1)= r0(1)+f0*qmd2
        r0(2)= r0(2)+f0*qmd*ex33q
        r0(3)= r0(3)+f0*qmd*ex34q
        r0(4)= r0(4)+f0*qmd*ex44q
        r0(5)= r0(5)+f0*ex33q*ex44q*c34
!
        work(1)= qmd2*ex3q
        work(2)= qmd2*ex4q
        work(3)= qmd*ex34q*ex3q
        work(4)= qmd*ex34q*ex4q
        f1w(1)= f1(1)*xiq
        f1w(2)= f1(1)*yiq
        f1w(3)= f1(2)
        do j= 1,4
          do i= 1,3
            r1(i,j)= r1(i,j)+f1w(i)*work(j)
          enddo
        enddo
!
        work(1)= qmd3
        work(2)= qmd2*ex33q
        work(3)= qmd2*ex34q
        work(4)= qmd2*ex44q
        f2w(1)= f2(1)*xiq2+f1(1)
        f2w(2)= f2(1)*yiq2+f1(1)
        f2w(3)= f2(3)     +f1(1)
        f2w(4)= f2(1)*xyiq
        f2w(5)= f2(2)*xiq
        f2w(6)= f2(2)*yiq
        do j= 1,4
          do i= 1,6
            r2(i,j)= r2(i,j)+f2w(i)*work(j)
          enddo
        enddo
! 
        work(1)= qmd3*ex3q
        work(2)= qmd3*ex4q
        f3w( 1)= (f3(1)*xiq2+f2(1)*three)*xiq
        f3w( 2)= (f3(1)*xiq2+f2(1)      )*yiq
        f3w( 3)= (f3(2)*xiq2+f2(2)      )
        f3w( 4)= (f3(1)*yiq2+f2(1)      )*xiq
        f3w( 5)= (f3(2)*xyiq            )
        f3w( 6)= (f3(3)     +f2(1)      )*xiq
        f3w( 7)= (f3(1)*yiq2+f2(1)*three)*yiq
        f3w( 8)= (f3(2)*yiq2+f2(2)      )
        f3w( 9)= (f3(3)     +f2(1)      )*yiq
        f3w(10)= (f3(4)     +f2(2)*three)
        do j= 1,2
          do i= 1,10
            r3(i,j)= r3(i,j)+f3w(i)*work(j)
          enddo
        enddo
!
        r4( 1)= r4( 1)+(f4(1)*xiq4 +f3(1)*xiq2*six       +f2(1)*three)*qmd4
        r4( 2)= r4( 2)+(f4(1)*xiq2 +f3(1)*three                      )*qmd4xy
        r4( 3)= r4( 3)+(f4(2)*xiq2 +f3(2)*three                      )*qmd4x
        r4( 4)= r4( 4)+(f4(1)*xyiq2+f3(1)*xiq2+f3(1)*yiq2+f2(1)      )*qmd4
        r4( 5)= r4( 5)+(f4(2)*xiq2 +f3(2)                            )*qmd4y
        r4( 6)= r4( 6)+(f4(3)*xiq2 +f3(1)*xiq2+f3(3)     +f2(1)      )*qmd4
        r4( 7)= r4( 7)+(f4(1)*yiq2 +f3(1)*three                      )*qmd4xy
        r4( 8)= r4( 8)+(f4(2)*yiq2 +f3(2)                            )*qmd4x
        r4( 9)= r4( 9)+(f4(3)      +f3(1)                            )*qmd4xy
        r4(10)= r4(10)+(f4(4)      +f3(2)*three                      )*qmd4x
        r4(11)= r4(11)+(f4(1)*yiq4 +f3(1)*yiq2*six       +f2(1)*three)*qmd4
        r4(12)= r4(12)+(f4(2)*yiq2 +f3(2)*three                      )*qmd4y
        r4(13)= r4(13)+(f4(3)*yiq2 +f3(1)*yiq2+f3(3)     +f2(1)      )*qmd4
        r4(14)= r4(14)+(f4(4)      +f3(2)*three                      )*qmd4y
        r4(15)= r4(15)+(f4(5)      +f3(3)*six            +f2(1)*three)*qmd4
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
      eri(1,1)= r4(1)+r2(1,1)*six+r0(1)*three+(+r3(1,1)*two+r3(1,2)*two+r1(1,1)*six &
&              +r1(1,2)*six)*qx+(+r2(1,2)+r2(1,3)*four+r2(1,4)+r0(2)+r0(3)*four+r0(4))*xx+( &
&              +r1(1,3)*two+r1(1,4)*two)*xxx+r0(5)*xxxx
      eri(2,1)= r4(4)+r2(1,1)+r2(2,1)+r0(1)+(+r3(4,2)*two+r1(1,2)*two)*qx+(+r2(2,4)+r0(4)) &
&              *xx
      eri(3,1)= r4(6)+r2(1,1)+r2(3,1)+r0(1)+(+r3(6,2)*two+r1(1,2)*two)*qx+(+r3(3,1)*two &
&              +r1(3,1)*two)*qz+(+r2(3,4)+r0(4))*xx+r2(5,3)*four*xz+(+r2(1,2)+r0(2))*zz &
&              +r1(3,4)*two*xxz+r1(1,3)*two*xzz+r0(5)*xxzz
      eri(4,1)= r4(2)+r2(4,1)*three+(+r3(2,1)+r3(2,2)*two+r1(2,1)+r1(2,2)*two)*qx+( &
&              +r2(4,3)*two+r2(4,4))*xx+r1(2,4)*xxx
      eri(5,1)= r4(3)+r2(5,1)*three+(+r3(3,1)+r3(3,2)*two+r1(3,1)+r1(3,2)*two)*qx+( &
&              +r3(1,1)+r1(1,1)*three)*qz+(+r2(5,3)*two+r2(5,4))*xx+(+r2(1,2)+r2(1,3)*two &
&              +r0(2)+r0(3)*two)*xz+r1(3,4)*xxx+(+r1(1,3)*two+r1(1,4))*xxz+r0(5)*xxxz
      eri(6,1)= r4(5)+r2(6,1)+r3(5,2)*two*qx+(+r3(2,1)+r1(2,1))*qz+r2(6,4)*xx+r2(4,3)*two &
&              *xz+r1(2,4)*xxz
      eri(1,2)= r4(4)+r2(1,1)+r2(2,1)+r0(1)+(+r3(4,1)*two+r1(1,1)*two)*qx+(+r2(2,2)+r0(2)) &
&              *xx
      eri(2,2)= r4(11)+r2(2,1)*six+r0(1)*three
      eri(3,2)= r4(13)+r2(2,1)+r2(3,1)+r0(1)+(+r3(8,1)*two+r1(3,1)*two)*qz+(+r2(2,2)+r0(2) &
&              )*zz
      eri(4,2)= r4(7)+r2(4,1)*three+(+r3(7,1)+r1(2,1)*three)*qx
      eri(5,2)= r4(8)+r2(5,1)+(+r3(8,1)+r1(3,1))*qx+(+r3(4,1)+r1(1,1))*qz+(+r2(2,2)+r0(2)) &
&              *xz
      eri(6,2)= r4(12)+r2(6,1)*three+(+r3(7,1)+r1(2,1)*three)*qz
      eri(1,3)= r4(6)+r2(1,1)+r2(3,1)+r0(1)+(+r3(6,1)*two+r1(1,1)*two)*qx+(+r3(3,2)*two &
&              +r1(3,2)*two)*qz+(+r2(3,2)+r0(2))*xx+r2(5,3)*four*xz+(+r2(1,4)+r0(4))*zz &
&              +r1(3,3)*two*xxz+r1(1,4)*two*xzz+r0(5)*xxzz
      eri(2,3)= r4(13)+r2(2,1)+r2(3,1)+r0(1)+(+r3(8,2)*two+r1(3,2)*two)*qz+(+r2(2,4)+r0(4) &
&              )*zz
      eri(3,3)= r4(15)+r2(3,1)*six+r0(1)*three+(+r3(10,1)*two+r3(10,2)*two+r1(3,1)*six &
&              +r1(3,2)*six)*qz+(+r2(3,2)+r2(3,3)*four+r2(3,4)+r0(2)+r0(3)*four+r0(4))*zz+( &
&              +r1(3,3)*two+r1(3,4)*two)*zzz+r0(5)*zzzz
      eri(4,3)= r4(9)+r2(4,1)+(+r3(9,1)+r1(2,1))*qx+r3(5,2)*two*qz+r2(6,3)*two*xz+r2(4,4) &
&              *zz+r1(2,4)*xzz
      eri(5,3)= r4(10)+r2(5,1)*three+(+r3(10,1)+r1(3,1)*three)*qx+(+r3(6,1)+r3(6,2)*two &
&              +r1(1,1)+r1(1,2)*two)*qz+(+r2(3,2)+r2(3,3)*two+r0(2)+r0(3)*two)*xz+( &
&              +r2(5,3)*two+r2(5,4))*zz+(+r1(3,3)*two+r1(3,4))*xzz+r1(1,4)*zzz+r0(5)*xzzz
      eri(6,3)= r4(14)+r2(6,1)*three+(+r3(9,1)+r3(9,2)*two+r1(2,1)+r1(2,2)*two)*qz+( &
&              +r2(6,3)*two+r2(6,4))*zz+r1(2,4)*zzz
      eri(1,4)= r4(2)+r2(4,1)*three+(+r3(2,1)*two+r3(2,2)+r1(2,1)*two+r1(2,2))*qx+( &
&              +r2(4,2)+r2(4,3)*two)*xx+r1(2,3)*xxx
      eri(2,4)= r4(7)+r2(4,1)*three+(+r3(7,2)+r1(2,2)*three)*qx
      eri(3,4)= r4(9)+r2(4,1)+(+r3(9,2)+r1(2,2))*qx+r3(5,1)*two*qz+r2(6,3)*two*xz+r2(4,2) &
&              *zz+r1(2,3)*xzz
      eri(4,4)= r4(4)+r2(1,1)+r2(2,1)+r0(1)+(+r3(4,1)+r3(4,2)+r1(1,1)+r1(1,2))*qx+( &
&              +r2(2,3)+r0(3))*xx
      eri(5,4)= r4(5)+r2(6,1)+(+r3(5,1)+r3(5,2))*qx+(+r3(2,1)+r1(2,1))*qz+r2(6,3)*xx+( &
&              +r2(4,2)+r2(4,3))*xz+r1(2,3)*xxz
      eri(6,4)= r4(8)+r2(5,1)+(+r3(8,2)+r1(3,2))*qx+(+r3(4,1)+r1(1,1))*qz+(+r2(2,3)+r0(3)) &
&              *xz
      eri(1,5)= r4(3)+r2(5,1)*three+(+r3(3,1)*two+r3(3,2)+r1(3,1)*two+r1(3,2))*qx+( &
&              +r3(1,2)+r1(1,2)*three)*qz+(+r2(5,2)+r2(5,3)*two)*xx+(+r2(1,3)*two+r2(1,4) &
&              +r0(3)*two+r0(4))*xz+r1(3,3)*xxx+(+r1(1,3)+r1(1,4)*two)*xxz+r0(5)*xxxz
      eri(2,5)= r4(8)+r2(5,1)+(+r3(8,2)+r1(3,2))*qx+(+r3(4,2)+r1(1,2))*qz+(+r2(2,4)+r0(4)) &
&              *xz
      eri(3,5)= r4(10)+r2(5,1)*three+(+r3(10,2)+r1(3,2)*three)*qx+(+r3(6,1)*two+r3(6,2) &
&              +r1(1,1)*two+r1(1,2))*qz+(+r2(3,3)*two+r2(3,4)+r0(3)*two+r0(4))*xz+(+r2(5,2) &
&              +r2(5,3)*two)*zz+(+r1(3,3)+r1(3,4)*two)*xzz+r1(1,3)*zzz+r0(5)*xzzz
      eri(4,5)= r4(5)+r2(6,1)+(+r3(5,1)+r3(5,2))*qx+(+r3(2,2)+r1(2,2))*qz+r2(6,3)*xx+( &
&              +r2(4,3)+r2(4,4))*xz+r1(2,4)*xxz
      eri(5,5)= r4(6)+r2(1,1)+r2(3,1)+r0(1)+(+r3(6,1)+r3(6,2)+r1(1,1)+r1(1,2))*qx+( &
&              +r3(3,1)+r3(3,2)+r1(3,1)+r1(3,2))*qz+(+r2(3,3)+r0(3))*xx+(+r2(5,2) &
&              +r2(5,3)*two+r2(5,4))*xz+(+r2(1,3)+r0(3))*zz+(+r1(3,3)+r1(3,4))*xxz+( &
&              +r1(1,3)+r1(1,4))*xzz+r0(5)*xxzz
      eri(6,5)= r4(9)+r2(4,1)+(+r3(9,2)+r1(2,2))*qx+(+r3(5,1)+r3(5,2))*qz+(+r2(6,3) &
&              +r2(6,4))*xz+r2(4,3)*zz+r1(2,4)*xzz
      eri(1,6)= r4(5)+r2(6,1)+r3(5,1)*two*qx+(+r3(2,2)+r1(2,2))*qz+r2(6,2)*xx+r2(4,3)*two &
&              *xz+r1(2,3)*xxz
      eri(2,6)= r4(12)+r2(6,1)*three+(+r3(7,2)+r1(2,2)*three)*qz
      eri(3,6)= r4(14)+r2(6,1)*three+(+r3(9,1)*two+r3(9,2)+r1(2,1)*two+r1(2,2))*qz+( &
&              +r2(6,2)+r2(6,3)*two)*zz+r1(2,3)*zzz
      eri(4,6)= r4(8)+r2(5,1)+(+r3(8,1)+r1(3,1))*qx+(+r3(4,2)+r1(1,2))*qz+(+r2(2,3)+r0(3)) &
&              *xz
      eri(5,6)= r4(9)+r2(4,1)+(+r3(9,1)+r1(2,1))*qx+(+r3(5,1)+r3(5,2))*qz+(+r2(6,2) &
&              +r2(6,3))*xz+r2(4,3)*zz+r1(2,3)*xzz
      eri(6,6)= r4(13)+r2(2,1)+r2(3,1)+r0(1)+(+r3(8,1)+r3(8,2)+r1(3,1)+r1(3,2))*qz+( &
&              +r2(2,3)+r0(3))*zz
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
      if(nbfijkl(3) == 6)then
        do l= 1,6
          do k= 1,6
            work(k)= eri(l,k)
          enddo
          do k= 1,6
            eri(l,k)= work(1)*rot2(1,k)+work(2)*rot2(2,k)+work(3)*rot2(3,k) &
&                    +work(4)*rot2(4,k)+work(5)*rot2(5,k)+work(6)*rot2(6,k) 
          enddo
        enddo
      else
        do l= 1,6
          do k= 1,6
            work(k)= eri(l,k)
          enddo
          do k= 1,5
            eri(l,k)= work(1)*rot3(1,k)+work(2)*rot3(2,k)+work(3)*rot3(3,k) &
&                    +work(4)*rot3(4,k)+work(5)*rot3(5,k)+work(6)*rot3(6,k)
          enddo
        enddo
      endif
!
      if(nbfijkl(4) == 6)then
        do k= 1,nbfijkl(3)
          do l= 1,6
            phmdint(l,k,1,1)= eri(1,k)*rot2(1,l)+eri(2,k)*rot2(2,l)+eri(3,k)*rot2(3,l) &
&                            +eri(4,k)*rot2(4,l)+eri(5,k)*rot2(5,l)+eri(6,k)*rot2(6,l) 
          enddo
        enddo
      else
        do k= 1,nbfijkl(3)
          do l= 1,5
            phmdint(l,k,1,1)= eri(1,k)*rot3(1,l)+eri(2,k)*rot3(2,l)+eri(3,k)*rot3(3,l) &
&                            +eri(4,k)*rot3(4,l)+eri(5,k)*rot3(5,l)+eri(6,k)*rot3(6,l) 
          enddo
        enddo
      endif
!
      return
end
