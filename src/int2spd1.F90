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
  subroutine int2dsss(phmdint,exfac1,exfac2,xyziq,xzkl,rot,nijkl,nbfijkl)
!--------------------------------------------------------------------------
!
! Calculate (ds|ss) integrals
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
      integer :: ij, kl, igrid, l, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, sqrtpi4=0.8862269254527580D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00
      real(8),parameter :: sqrt3=1.73205080756888D+00, sqrt3h=8.660254037844386D-01
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*), xzkl(2), rot(3,3)
      real(8),intent(out) :: phmdint(6,6,6,6)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:2)
      real(8) :: f0, f1(2), f2(3), r0(2), r1(3), r2(6)
      real(8) :: ex12, ex34, ex43, ex14, ex41, expq, ex3q, ex4q, c12, c34, zip, xiq, yiq, ziq
      real(8) :: xiq2, yiq2, xypq2, zpq, zpq2, fac, qmd, qmd2, qx, qz, xx, xz, zz
      real(8) :: eri(6), work(6), rot2(6,6), rot3(6,5)
!
! Zero-clear
!
      r0(1:2)= zero
      r1(1:3)= zero
      r2(1:6)= zero
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
        xypq2= xiq2+yiq2
        f0   = zero
        f1(1:2)= zero
        f2(1:3)= zero
        do ij= 1,nijkl(1)
          ex12= exfac1(1,ij)
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
            do ii= 0,2
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
            fac= c12*ex41
            expq= expq*two*ex41*ex41
            ft(0)= ft(0)*fac
            ft(1)= ft(1)*fac*expq
            ft(2)= ft(2)*fac*expq*expq
          endif
          f0= f0+ft(0)
          f1(1)= f1(1)-ft(1)
          f1(2)= f1(2)-ft(1)*zpq
          f2(1)= f2(1)+ft(2)
          f2(2)= f2(2)+ft(2)*zpq
          f2(3)= f2(3)+ft(2)*zpq2
        enddo
!
        qmd = ex43*c34
        qmd2= qmd*ex43
        work(1)= qmd*ex3q
        r0(1)= r0(1)+f0*qmd
        r0(2)= r0(2)+f0*ex3q*ex3q*c34
        r1(1)= r1(1)+f1(1)*work(1)*xiq
        r1(2)= r1(2)+f1(1)*work(1)*yiq
        r1(3)= r1(3)+f1(2)*work(1)
        r2(1)= r2(1)+(f2(1)*xiq2+f1(1))*qmd2
        r2(2)= r2(2)+(f2(1)*yiq2+f1(1))*qmd2
        r2(3)= r2(3)+(f2(3)     +f1(1))*qmd2
        r2(4)= r2(4)+(f2(1)*xiq*yiq   )*qmd2
        r2(5)= r2(5)+(f2(2)*xiq       )*qmd2
        r2(6)= r2(6)+(f2(2)*yiq       )*qmd2
      enddo
!
      qx= xzkl(1)
      qz= xzkl(2)
      xx= qx*qx
      xz= qx*qz
      zz= qz*qz
      eri(1)= r2(1)+r1(1)*qx*two+r0(1)+r0(2)*xx
      eri(2)= r2(2)             +r0(1)
      eri(3)= r2(3)+r1(3)*qz*two+r0(1)+r0(2)*zz
      eri(4)= r2(4)+r1(2)*qx
      eri(5)= r2(5)+r1(1)*qz+r1(3)*qx +r0(2)*xz
      eri(6)= r2(6)+r1(2)*qz
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
      if(nbfijkl(4) == 6)then
        do l= 1,6
          phmdint(l,1,1,1)= eri(1)*rot2(1,l)+eri(2)*rot2(2,l)+eri(3)*rot2(3,l) &
&                          +eri(4)*rot2(4,l)+eri(5)*rot2(5,l)+eri(6)*rot2(6,l) 
        enddo
      else
        do l= 1,6
          rot3(l,1)= rot2(l,2)
          rot3(l,2)= rot2(l,5)
          rot3(l,3)= rot2(l,6)-(rot2(l,1)+rot2(l,4))*half
          rot3(l,4)= rot2(l,3)
          rot3(l,5)=(rot2(l,1)-rot2(l,4))*sqrt3h
        enddo
        do l= 1,5
          phmdint(l,1,1,1)= eri(1)*rot3(1,l)+eri(2)*rot3(2,l)+eri(3)*rot3(3,l) &
&                          +eri(4)*rot3(4,l)+eri(5)*rot3(5,l)+eri(6)*rot3(6,l) 
        enddo
      endif
!
      return
end


!--------------------------------------------------------------------------
  subroutine int2dpss(phmdint,exfac1,exfac2,xyziq,xzkl,rot,nijkl,nbfijkl)
!--------------------------------------------------------------------------
!
! Calculate (dp|ss) integrals
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
      integer :: ij, kl, igrid, i, k, l, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, sqrtpi4=0.8862269254527580D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00, five=5.0D+00
      real(8),parameter :: sqrt3=1.73205080756888D+00, sqrt3h=8.660254037844386D-01
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*), xzkl(2), rot(3,3)
      real(8),intent(out) :: phmdint(6,6,6,6)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:3)
      real(8) :: f0, f1(2), f2(3), f3(4), r0(3), r1(3,3), r2(6,2), r3(10)
      real(8) :: ex12, ex34, ex43, ex14, ex41, expq, expq2, ex3q, ex4q, c12, c34, zip, xiq, yiq
      real(8) :: ziq, xiq2, yiq2, xyiq, xypq2, zpq, zpq2, fac, qmd, qmd2, qmd3, qmd3x, qmd3y
      real(8) :: qx, qz, xx, xz, zz, xxx, xxz, xzz, zzz, eri(6,3), work(6), rot2(6,6), rot3(6,5)
!
! Zero-clear
!
      r0(1:3)=     zero
      r1(1:3,1:3)= zero
      r2(1:6,1:2)= zero
      r3(1:10)=    zero
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
        xypq2= xiq2+yiq2
        f0   = zero
        f1(1:2)= zero
        f2(1:3)= zero
        f3(1:4)= zero
        do ij= 1,nijkl(1)
          ex12= exfac1(1,ij)
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
            do ii= 0,3
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
            fac= c12*ex41
            expq= expq*two*ex41*ex41
            expq2= expq*expq
            ft(0)= ft(0)*fac
            ft(1)= ft(1)*fac*expq
            ft(2)= ft(2)*fac*expq2
            ft(3)= ft(3)*fac*expq2*expq
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
        enddo
!
        qmd = ex43*c34
        qmd2= qmd*ex43
        qmd3= qmd2*ex43
        qmd3x=qmd3*xiq
        qmd3y=qmd3*yiq
        work(1)= qmd2
        work(2)= qmd*ex3q*ex3q
        work(3)= qmd*ex3q*ex4q
        work(4)= qmd2*ex3q
        work(5)= qmd2*ex4q
!
        r0(1)= r0(1)+f0*qmd*ex3q
        r0(2)= r0(2)+f0*qmd*ex4q
        r0(3)= r0(3)+f0*ex3q*ex3q*ex4q*c34
!
        do i= 1,3
          r1(1,i)= r1(1,i)+f1(1)*work(i)*xiq
          r1(2,i)= r1(2,i)+f1(1)*work(i)*yiq
          r1(3,i)= r1(3,i)+f1(2)*work(i)
        enddo
!
        do i= 1,2
          r2(1,i)= r2(1,i)+(f2(1)*xiq2+f1(1))*work(i+3)
          r2(2,i)= r2(2,i)+(f2(1)*yiq2+f1(1))*work(i+3)
          r2(3,i)= r2(3,i)+(f2(3)     +f1(1))*work(i+3)
          r2(4,i)= r2(4,i)+(f2(1)*xiq*yiq   )*work(i+3)
          r2(5,i)= r2(5,i)+(f2(2)*xiq       )*work(i+3)
          r2(6,i)= r2(6,i)+(f2(2)*yiq       )*work(i+3)
        enddo
!
        r3( 1)= r3( 1)+(f3(1)*xiq2+f2(1)*three)*qmd3x
        r3( 2)= r3( 2)+(f3(1)*xiq2+f2(1)      )*qmd3y
        r3( 3)= r3( 3)+(f3(2)*xiq2+f2(2)      )*qmd3
        r3( 4)= r3( 4)+(f3(1)*yiq2+f2(1)      )*qmd3x
        r3( 5)= r3( 5)+(f3(2)*xyiq            )*qmd3
        r3( 6)= r3( 6)+(f3(3)     +f2(1)      )*qmd3x
        r3( 7)= r3( 7)+(f3(1)*yiq2+f2(1)*three)*qmd3y
        r3( 8)= r3( 8)+(f3(2)*yiq2+f2(2)      )*qmd3
        r3( 9)= r3( 9)+(f3(3)     +f2(1)      )*qmd3y
        r3(10)= r3(10)+(f3(4)     +f2(2)*three)*qmd3
      enddo
!
      qx= xzkl(1)
      qz= xzkl(2)
      xx= qx*qx
      xz= qx*qz
      zz= qz*qz
      xxx=qx*qx*qx
      xxz=qx*qx*qz
      xzz=qx*qz*qz
      zzz=qz*qz*qz
      eri(1,1)= r3(1)+r1(1,1)*three+(+r2(1,1)*two+r2(1,2)+r0(1)*two+r0(2))*qx+(+r1(1,2) &
&              +r1(1,3)*two)*xx+r0(3)*xxx
      eri(2,1)= r3(4)+r1(1,1)+(+r2(2,2)+r0(2))*qx
      eri(3,1)= r3(6)+r1(1,1)+(+r2(3,2)+r0(2))*qx+r2(5,1)*two*qz+r1(3,3)*two*xz+r1(1,2)*zz &
&              +r0(3)*xzz
      eri(4,1)= r3(2)+r1(2,1)+(+r2(4,1)+r2(4,2))*qx+r1(2,3)*xx
      eri(5,1)= r3(3)+r1(3,1)+(+r2(5,1)+r2(5,2))*qx+(+r2(1,1)+r0(1))*qz+r1(3,3)*xx+( &
&              +r1(1,2)+r1(1,3))*xz+r0(3)*xxz
      eri(6,1)= r3(5)+r2(6,2)*qx+r2(4,1)*qz+r1(2,3)*xz
      eri(1,2)= r3(2)+r1(2,1)+r2(4,1)*two*qx+r1(2,2)*xx
      eri(2,2)= r3(7)+r1(2,1)*three
      eri(3,2)= r3(9)+r1(2,1)+r2(6,1)*two*qz+r1(2,2)*zz
      eri(4,2)= r3(4)+r1(1,1)+(+r2(2,1)+r0(1))*qx
      eri(5,2)= r3(5)+r2(6,1)*qx+r2(4,1)*qz+r1(2,2)*xz
      eri(6,2)= r3(8)+r1(3,1)+(+r2(2,1)+r0(1))*qz
      eri(1,3)= r3(3)+r1(3,1)+r2(5,1)*two*qx+(+r2(1,2)+r0(2))*qz+r1(3,2)*xx+r1(1,3)*two*xz &
&              +r0(3)*xxz
      eri(2,3)= r3(8)+r1(3,1)+(+r2(2,2)+r0(2))*qz
      eri(3,3)= r3(10)+r1(3,1)*three+(+r2(3,1)*two+r2(3,2)+r0(1)*two+r0(2))*qz+(+r1(3,2) &
&              +r1(3,3)*two)*zz+r0(3)*zzz
      eri(4,3)= r3(5)+r2(6,1)*qx+r2(4,2)*qz+r1(2,3)*xz
      eri(5,3)= r3(6)+r1(1,1)+(+r2(3,1)+r0(1))*qx+(+r2(5,1)+r2(5,2))*qz+(+r1(3,2)+r1(3,3)) &
&              *xz+r1(1,3)*zz+r0(3)*xzz
      eri(6,3)= r3(9)+r1(2,1)+(+r2(6,1)+r2(6,2))*qz+r1(2,3)*zz
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
      do l=1,6
        work(1)= eri(l,1)
        work(2)= eri(l,2)
        work(3)= eri(l,3)
        eri(l,1)= work(1)*rot(1,1)+work(2)*rot(2,1)+work(3)*rot(3,1)
        eri(l,2)= work(1)*rot(1,2)+work(2)*rot(2,2)+work(3)*rot(3,2)
        eri(l,3)= work(1)*rot(1,3)+work(2)*rot(2,3)+work(3)*rot(3,3)
      enddo
!
      if(nbfijkl(4) == 6)then
        do k= 1,3
          do l= 1,6
            phmdint(l,k,1,1)= eri(1,k)*rot2(1,l)+eri(2,k)*rot2(2,l)+eri(3,k)*rot2(3,l) &
&                            +eri(4,k)*rot2(4,l)+eri(5,k)*rot2(5,l)+eri(6,k)*rot2(6,l) 
          enddo
        enddo
      else
        do l= 1,6
          rot3(l,1)= rot2(l,2)
          rot3(l,2)= rot2(l,5)
          rot3(l,3)= rot2(l,6)-(rot2(l,1)+rot2(l,4))*half
          rot3(l,4)= rot2(l,3)
          rot3(l,5)=(rot2(l,1)-rot2(l,4))*sqrt3h
        enddo
        do k= 1,3
          do l= 1,5
            phmdint(l,k,1,1)= eri(1,k)*rot3(1,l)+eri(2,k)*rot3(2,l)+eri(3,k)*rot3(3,l) &
&                            +eri(4,k)*rot3(4,l)+eri(5,k)*rot3(5,l)+eri(6,k)*rot3(6,l) 
          enddo
        enddo
      endif
!
      return
end


!--------------------------------------------------------------------------
  subroutine int2dsps(phmdint,exfac1,exfac2,xyziq,xzkl,rot,nijkl,nbfijkl)
!--------------------------------------------------------------------------
!
! Calculate (dp|ss) integrals
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
      integer :: ij, kl, igrid, i, j, l, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, sqrtpi4=0.8862269254527580D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00, five=5.0D+00
      real(8),parameter :: sqrt3=1.73205080756888D+00, sqrt3h=8.660254037844386D-01
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*), xzkl(2), rot(3,3)
      real(8),intent(out) :: phmdint(6,6,6,6)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:3)
      real(8) :: f0, f1(2,2), f2(3,2), f3(4), r0(2), r1(3,3), r2(6,2), r3(10)
      real(8) :: ex12, ex34, ex43, ex14, ex41, expq, expq2, ex3q, ex4q, c12, c34, zip, xiq, yiq
      real(8) :: ziq, xiq2, yiq2, xyiq, xypq2, zpq, zpq2, fac, qmd, qmd2, qmd2x, qmd2y, pmd, zjp
      real(8) :: qx, qz, xx, xz, zz, eri(6,3), work(6), rot2(6,6), rot3(6,5)
!
! Zero-clear
!
      r0(1:2)=     zero
      r1(1:3,1:3)= zero
      r2(1:6,1:2)= zero
      r3(1:10)=    zero
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
        xypq2= xiq2+yiq2
        f0= zero
        f1(1:2,1:2)= zero
        f2(1:3,1:2)= zero
        f3(1:4)= zero
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
            do ii= 0,3
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
            fac= c12*ex41
            expq= expq*two*ex41*ex41
            expq2= expq*expq
            ft(0)= ft(0)*fac
            ft(1)= ft(1)*fac*expq
            ft(2)= ft(2)*fac*expq2
            ft(3)= ft(3)*fac*expq2*expq
          endif
          f0     = f0     +ft(0)*zjp
          f1(1,1)= f1(1,1)-ft(1)*zjp
          f1(2,1)= f1(2,1)-ft(1)*zjp*zpq
          f1(1,2)= f1(1,2)-ft(1)*pmd
          f1(2,2)= f1(2,2)-ft(1)*pmd*zpq
          f2(1,1)= f2(1,1)+ft(2)*zjp
          f2(2,1)= f2(2,1)+ft(2)*zjp*zpq
          f2(3,1)= f2(3,1)+ft(2)*zjp*zpq2
          f2(1,2)= f2(1,2)+ft(2)*pmd
          f2(2,2)= f2(2,2)+ft(2)*pmd*zpq
          f2(3,2)= f2(3,2)+ft(2)*pmd*zpq2
          f3(1)  = f3(1)  -ft(3)*pmd
          f3(2)  = f3(2)  -ft(3)*pmd*zpq
          f3(3)  = f3(3)  -ft(3)*pmd*zpq2
          f3(4)  = f3(4)  -ft(3)*pmd*zpq2*zpq
        enddo
!
        qmd = ex43*c34
        qmd2= qmd*ex43
        qmd2x=qmd2*xiq
        qmd2y=qmd2*yiq
        work(1)= qmd
        work(2)= ex3q*ex3q*c34
        work(3)= qmd2
        work(4)= qmd*ex3q
!
        r0(1)= r0(1)+f0*work(1)
        r0(2)= r0(2)+f0*work(2)
!
        r1(1,1)= r1(1,1)+f1(1,1)*work(4)*xiq
        r1(2,1)= r1(2,1)+f1(1,1)*work(4)*yiq
        r1(3,1)= r1(3,1)+f1(2,1)*work(4)
        r1(1,2)= r1(1,2)+f1(1,2)*work(1)*xiq
        r1(2,2)= r1(2,2)+f1(1,2)*work(1)*yiq
        r1(3,2)= r1(3,2)+f1(2,2)*work(1)
        r1(1,3)= r1(1,3)+f1(1,2)*work(2)*xiq
        r1(2,3)= r1(2,3)+f1(1,2)*work(2)*yiq
        r1(3,3)= r1(3,3)+f1(2,2)*work(2)
!
        do i= 1,2
          r2(1,i)= r2(1,i)+(f2(1,i)*xiq2+f1(1,i))*work(i+2)
          r2(2,i)= r2(2,i)+(f2(1,i)*yiq2+f1(1,i))*work(i+2)
          r2(3,i)= r2(3,i)+(f2(3,i)     +f1(1,i))*work(i+2)
          r2(4,i)= r2(4,i)+(f2(1,i)*xiq*yiq     )*work(i+2)
          r2(5,i)= r2(5,i)+(f2(2,i)*xiq         )*work(i+2)
          r2(6,i)= r2(6,i)+(f2(2,i)*yiq         )*work(i+2)
        enddo
!
        r3( 1)= r3( 1)+(f3(1)*xiq2+f2(1,2)*three)*qmd2x
        r3( 2)= r3( 2)+(f3(1)*xiq2+f2(1,2)      )*qmd2y
        r3( 3)= r3( 3)+(f3(2)*xiq2+f2(2,2)      )*qmd2
        r3( 4)= r3( 4)+(f3(1)*yiq2+f2(1,2)      )*qmd2x
        r3( 5)= r3( 5)+(f3(2)*xyiq              )*qmd2
        r3( 6)= r3( 6)+(f3(3)     +f2(1,2)      )*qmd2x
        r3( 7)= r3( 7)+(f3(1)*yiq2+f2(1,2)*three)*qmd2y
        r3( 8)= r3( 8)+(f3(2)*yiq2+f2(2,2)      )*qmd2
        r3( 9)= r3( 9)+(f3(3)     +f2(1,2)      )*qmd2y
        r3(10)= r3(10)+(f3(4)     +f2(2,2)*three)*qmd2
      enddo
!
      qx= xzkl(1)
      qz= xzkl(2)
      xx= qx*qx
      xz= qx*qz
      zz= qz*qz
      eri(1,1)= -r3(1)-r1(1,2)-r2(1,2)*two*qx-r1(1,3)*xx
      eri(2,1)= -r3(4)-r1(1,2)
      eri(3,1)= -r3(6)-r1(1,2)-r2(5,2)*two*qz-r1(1,3)*zz
      eri(4,1)= -r3(2)-r2(4,2)*qx
      eri(5,1)= -r3(3)-r2(5,2)*qx-r2(1,2)*qz-r1(1,3)*xz
      eri(6,1)= -r3(5)-r2(4,2)*qz
      eri(1,2)= -r3(2)-r1(2,2)-r2(4,2)*two*qx-r1(2,3)*xx
      eri(2,2)= -r3(7)-r1(2,2)
      eri(3,2)= -r3(9)-r1(2,2)-r2(6,2)*two*qz-r1(2,3)*zz
      eri(4,2)= -r3(4)-r2(2,2)*qx
      eri(5,2)= -r3(5)-r2(6,2)*qx-r2(4,2)*qz-r1(2,3)*xz
      eri(6,2)= -r3(8)-r2(2,2)*qz
      eri(1,3)= -r3(3)+r2(1,1)-r1(3,2)+r0(1)+(-r2(5,2)*two+r1(1,1)*two)*qx+(-r1(3,3)+r0(2))*xx
      eri(2,3)= -r3(8)+r2(2,1)-r1(3,2)+r0(1)
      eri(3,3)= -r3(10)+r2(3,1)-r1(3,2)+r0(1)+(-r2(3,2)*two+r1(3,1)*two)*qz+(-r1(3,3)+r0(2))*zz
      eri(4,3)= -r3(5)+r2(4,1)+(-r2(6,2)+r1(2,1))*qx
      eri(5,3)= -r3(6)+r2(5,1)+(-r2(3,2)+r1(3,1))*qx+(-r2(5,2)+r1(1,1))*qz+(-r1(3,3)+r0(2))*xz
      eri(6,3)= -r3(9)+r2(6,1)+(-r2(6,2)+r1(2,1))*qz
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
      do l=1,6
        work(1)= eri(l,1)
        work(2)= eri(l,2)
        work(3)= eri(l,3)
        eri(l,1)= work(1)*rot(1,1)+work(2)*rot(2,1)+work(3)*rot(3,1)
        eri(l,2)= work(1)*rot(1,2)+work(2)*rot(2,2)+work(3)*rot(3,2)
        eri(l,3)= work(1)*rot(1,3)+work(2)*rot(2,3)+work(3)*rot(3,3)
      enddo
!
      if(nbfijkl(4) == 6) then
        do j= 1,3
          do l= 1,6
            phmdint(l,1,j,1)= eri(1,j)*rot2(1,l)+eri(2,j)*rot2(2,l)+eri(3,j)*rot2(3,l) &
&                            +eri(4,j)*rot2(4,l)+eri(5,j)*rot2(5,l)+eri(6,j)*rot2(6,l) 
          enddo
        enddo
      else
        do l= 1,6
          rot3(l,1)= rot2(l,2)
          rot3(l,2)= rot2(l,5)
          rot3(l,3)= rot2(l,6)-(rot2(l,1)+rot2(l,4))*half
          rot3(l,4)= rot2(l,3)
          rot3(l,5)=(rot2(l,1)-rot2(l,4))*sqrt3h
        enddo
        do j= 1,3
          do l= 1,5
            phmdint(l,1,j,1)= eri(1,j)*rot3(1,l)+eri(2,j)*rot3(2,l)+eri(3,j)*rot3(3,l) &
&                            +eri(4,j)*rot3(4,l)+eri(5,j)*rot3(5,l)+eri(6,j)*rot3(6,l) 
          enddo
        enddo
      endif
!
      return
end


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
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:4)
      real(8) :: f0, f1(2), f2(3), f3(4), f4(5), r0(5), r1(3,4), r2(6,4), r3(10,2), r4(15)
      real(8) :: ex12, ex34, ex43, ex14, ex41, expq, expq2, ex3q, ex4q, c12, c34, zip
      real(8) :: xiq, yiq, ziq, xiq2, yiq2, xyiq, xiq4, yiq4, xyiq2, xypq2, zpq, zpq2, fac
      real(8) :: ex33q, ex34q, ex44q, qmd, qmd2, qmd3, qmd4, qmd4x, qmd4y, qmd4xy
      real(8) :: qx, qz, xx, xz, zz, xxx, xxz, xzz, zzz, xxxx, xxxz, xxzz, xzzz, zzzz
      real(8) :: eri(6,6), work(6), f1w(3), f2w(6), f3w(10), rot2(6,6), rot3(6,5)
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
            do ii= 0,4
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
            fac= c12*ex41
            expq= expq*two*ex41*ex41
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
      if(nbfijkl(3) == 6) then
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


!--------------------------------------------------------------------------
  subroutine int2dpps(phmdint,exfac1,exfac2,xyziq,xzkl,rot,nijkl,nbfijkl)
!--------------------------------------------------------------------------
!
! Calculate (dp|ps) integrals
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
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:4)
      real(8) :: f0, f1(2,2), f2(3,2), f3(4,2), f4(5), r0(3), r1(3,6), r2(6,5), r3(10,3), r4(15)
      real(8) :: ex12, ex34, ex43, ex14, ex41, expq, expq2, ex3q, ex4q, c12, c34, zip
      real(8) :: xiq, yiq, ziq, xiq2, yiq2, xyiq, xiq4, yiq4, xyiq2, xypq2, zpq, zpq2, fac
      real(8) :: ex33q, ex34q, zjp, pmd, qmd, qmd2, qmd3, qmd3x, qmd3y, qmd3xy, qx, qz
      real(8) :: xx, xz, zz, xxx, xxz, xzz, zzz, eri(6,3,3), work(9), f2w(6,2), f3w(10,2)
      real(8) :: rot2(6,6), rot3(6,5)
!
! Zero-clear
!
      r0(1:3)     = zero
      r1(1:3,1:6) = zero
      r2(1:6,1:5) = zero
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
        f0   = zero
        f1(1:2,1:2)= zero
        f2(1:3,1:2)= zero
        f3(1:4,1:2)= zero
        f4(1:5)= zero
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
            do ii= 0,4
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
            fac= c12*ex41
            expq= expq*two*ex41*ex41
            expq2= expq*expq
            ft(0)= ft(0)*fac
            ft(1)= ft(1)*fac*expq
            ft(2)= ft(2)*fac*expq2
            ft(3)= ft(3)*fac*expq2*expq
            ft(4)= ft(4)*fac*expq2*expq2
          endif
          ft(4)= ft(4)*pmd
          f0     = f0     +ft(0)*zjp
          f1(1,1)= f1(1,1)-ft(1)*zjp
          f1(2,1)= f1(2,1)-ft(1)*zjp*zpq
          f1(1,2)= f1(1,2)-ft(1)*pmd
          f1(2,2)= f1(2,2)-ft(1)*pmd*zpq
          f2(1,1)= f2(1,1)+ft(2)*zjp
          f2(2,1)= f2(2,1)+ft(2)*zjp*zpq
          f2(3,1)= f2(3,1)+ft(2)*zjp*zpq2
          f2(1,2)= f2(1,2)+ft(2)*pmd
          f2(2,2)= f2(2,2)+ft(2)*pmd*zpq
          f2(3,2)= f2(3,2)+ft(2)*pmd*zpq2
          f3(1,1)= f3(1,1)-ft(3)*zjp
          f3(2,1)= f3(2,1)-ft(3)*zjp*zpq
          f3(3,1)= f3(3,1)-ft(3)*zjp*zpq2
          f3(4,1)= f3(4,1)-ft(3)*zjp*zpq2*zpq
          f3(1,2)= f3(1,2)-ft(3)*pmd
          f3(2,2)= f3(2,2)-ft(3)*pmd*zpq
          f3(3,2)= f3(3,2)-ft(3)*pmd*zpq2
          f3(4,2)= f3(4,2)-ft(3)*pmd*zpq2*zpq
          f4(1)  = f4(1)  +ft(4)
          f4(2)  = f4(2)  +ft(4)*zpq
          f4(3)  = f4(3)  +ft(4)*zpq2
          f4(4)  = f4(4)  +ft(4)*zpq2*zpq
          f4(5)  = f4(5)  +ft(4)*zpq2*zpq2
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
        work(1)= qmd3
        work(2)= qmd2*ex3q
        work(3)= qmd2*ex4q
        work(4)= qmd2
        work(5)= qmd*ex33q
        work(6)= qmd*ex34q
        work(7)= qmd*ex3q
        work(8)= qmd*ex4q
        work(9)= ex33q*ex4q*c34
!
        r0(1)= r0(1)+f0*work(7)
        r0(2)= r0(2)+f0*work(8)
        r0(3)= r0(3)+f0*work(9)
!
        do j= 1,3
          r1(1,j)= r1(1,j)+f1(1,1)*xiq*work(j+3)
          r1(2,j)= r1(2,j)+f1(1,1)*yiq*work(j+3)
          r1(3,j)= r1(3,j)+f1(2,1)    *work(j+3)
        enddo
        do j= 4,6
          r1(1,j)= r1(1,j)+f1(1,2)*xiq*work(j+3)
          r1(2,j)= r1(2,j)+f1(1,2)*yiq*work(j+3)
          r1(3,j)= r1(3,j)+f1(2,2)    *work(j+3)
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
        do j= 1,2
          do i= 1,6
            r2(i,j)= r2(i,j)+f2w(i,1)*work(j+1)
          enddo
        enddo
        do j= 3,5
          do i= 1,6
            r2(i,j)= r2(i,j)+f2w(i,2)*work(j+1)
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
        do i= 1,10
          r3(i,1)= r3(i,1)+f3w(i,1)*work(1)
        enddo
        do j= 2,3
          do i= 1,10
            r3(i,j)= r3(i,j)+f3w(i,2)*work(j)
          enddo
        enddo
!
        r4( 1)= r4( 1)+(f4(1)*xiq4 +f3(1,2)*xiq2*six         +f2(1,2)*three)*qmd3
        r4( 2)= r4( 2)+(f4(1)*xiq2 +f3(1,2)*three                          )*qmd3xy
        r4( 3)= r4( 3)+(f4(2)*xiq2 +f3(2,2)*three                          )*qmd3x
        r4( 4)= r4( 4)+(f4(1)*xyiq2+f3(1,2)*xiq2+f3(1,2)*yiq2+f2(1,2)      )*qmd3
        r4( 5)= r4( 5)+(f4(2)*xiq2 +f3(2,2)                                )*qmd3y
        r4( 6)= r4( 6)+(f4(3)*xiq2 +f3(1,2)*xiq2+f3(3,2)     +f2(1,2)      )*qmd3
        r4( 7)= r4( 7)+(f4(1)*yiq2 +f3(1,2)*three                          )*qmd3xy
        r4( 8)= r4( 8)+(f4(2)*yiq2 +f3(2,2)                                )*qmd3x
        r4( 9)= r4( 9)+(f4(3)      +f3(1,2)                                )*qmd3xy
        r4(10)= r4(10)+(f4(4)      +f3(2,2)*three                          )*qmd3x
        r4(11)= r4(11)+(f4(1)*yiq4 +f3(1,2)*yiq2*six         +f2(1,2)*three)*qmd3
        r4(12)= r4(12)+(f4(2)*yiq2 +f3(2,2)*three                          )*qmd3y
        r4(13)= r4(13)+(f4(3)*yiq2 +f3(1,2)*yiq2+f3(3,2)     +f2(1,2)      )*qmd3
        r4(14)= r4(14)+(f4(4)      +f3(2,2)*three                          )*qmd3y
        r4(15)= r4(15)+(f4(5)      +f3(3,2)*six              +f2(1,2)*three)*qmd3
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
      eri(1,1,1)=-r4(1)-r2(1,3)*three+(-r3(1,2)*two-r3(1,3)-r1(1,4)*two-r1(1,5))*qx+( &
&                -r2(1,4)-r2(1,5)*two)*xx-r1(1,6)*xxx
      eri(2,1,1)=-r4(4)-r2(1,3)+(-r3(4,3)-r1(1,5))*qx
      eri(3,1,1)=-r4(6)-r2(1,3)+(-r3(6,3)-r1(1,5))*qx-r3(3,2)*two*qz-r2(5,5)*two*xz &
&                -r2(1,4)*zz-r1(1,6)*xzz
      eri(4,1,1)=-r4(2)-r2(4,3)+(-r3(2,2)-r3(2,3))*qx-r2(4,5)*xx
      eri(5,1,1)=-r4(3)-r2(5,3)+(-r3(3,2)-r3(3,3))*qx+(-r3(1,2)-r1(1,4))*qz-r2(5,5)*xx+( &
&                -r2(1,4)-r2(1,5))*xz-r1(1,6)*xxz
      eri(6,1,1)=-r4(5)-r3(5,3)*qx-r3(2,2)*qz-r2(4,5)*xz
      eri(1,2,1)=-r4(2)-r2(4,3)-r3(2,2)*two*qx-r2(4,4)*xx
      eri(2,2,1)=-r4(7)-r2(4,3)*three
      eri(3,2,1)=-r4(9)-r2(4,3)-r3(5,2)*two*qz-r2(4,4)*zz
      eri(4,2,1)=-r4(4)-r2(1,3)+(-r3(4,2)-r1(1,4))*qx
      eri(5,2,1)=-r4(5)-r3(5,2)*qx-r3(2,2)*qz-r2(4,4)*xz
      eri(6,2,1)=-r4(8)-r2(5,3)+(-r3(4,2)-r1(1,4))*qz
      eri(1,3,1)=-r4(3)-r2(5,3)-r3(3,2)*two*qx+(-r3(1,3)-r1(1,5))*qz-r2(5,4)*xx &
&                -r2(1,5)*two*xz-r1(1,6)*xxz
      eri(2,3,1)=-r4(8)-r2(5,3)+(-r3(4,3)-r1(1,5))*qz
      eri(3,3,1)=-r4(10)-r2(5,3)*three+(-r3(6,2)*two-r3(6,3)-r1(1,4)*two-r1(1,5))*qz+( &
&                -r2(5,4)-r2(5,5)*two)*zz-r1(1,6)*zzz
      eri(4,3,1)=-r4(5)-r3(5,2)*qx-r3(2,3)*qz-r2(4,5)*xz
      eri(5,3,1)=-r4(6)-r2(1,3)+(-r3(6,2)-r1(1,4))*qx+(-r3(3,2)-r3(3,3))*qz+(-r2(5,4) &
&                -r2(5,5))*xz-r2(1,5)*zz-r1(1,6)*xzz
      eri(6,3,1)=-r4(9)-r2(4,3)+(-r3(5,2)-r3(5,3))*qz-r2(4,5)*zz
      eri(1,1,2)=-r4(2)-r2(4,3)*three+(-r3(2,2)*two-r3(2,3)-r1(2,4)*two-r1(2,5))*qx+( &
&                -r2(4,4)-r2(4,5)*two)*xx-r1(2,6)*xxx
      eri(2,1,2)=-r4(7)-r2(4,3)+(-r3(7,3)-r1(2,5))*qx
      eri(3,1,2)=-r4(9)-r2(4,3)+(-r3(9,3)-r1(2,5))*qx-r3(5,2)*two*qz-r2(6,5)*two*xz &
&                -r2(4,4)*zz-r1(2,6)*xzz
      eri(4,1,2)=-r4(4)-r2(2,3)+(-r3(4,2)-r3(4,3))*qx-r2(2,5)*xx
      eri(5,1,2)=-r4(5)-r2(6,3)+(-r3(5,2)-r3(5,3))*qx+(-r3(2,2)-r1(2,4))*qz-r2(6,5)*xx+( &
&                -r2(4,4)-r2(4,5))*xz-r1(2,6)*xxz
      eri(6,1,2)=-r4(8)-r3(8,3)*qx-r3(4,2)*qz-r2(2,5)*xz
      eri(1,2,2)=-r4(4)-r2(2,3)-r3(4,2)*two*qx-r2(2,4)*xx
      eri(2,2,2)=-r4(11)-r2(2,3)*three
      eri(3,2,2)=-r4(13)-r2(2,3)-r3(8,2)*two*qz-r2(2,4)*zz
      eri(4,2,2)=-r4(7)-r2(4,3)+(-r3(7,2)-r1(2,4))*qx
      eri(5,2,2)=-r4(8)-r3(8,2)*qx-r3(4,2)*qz-r2(2,4)*xz
      eri(6,2,2)=-r4(12)-r2(6,3)+(-r3(7,2)-r1(2,4))*qz
      eri(1,3,2)=-r4(5)-r2(6,3)-r3(5,2)*two*qx+(-r3(2,3)-r1(2,5))*qz-r2(6,4)*xx &
&                -r2(4,5)*two*xz-r1(2,6)*xxz
      eri(2,3,2)=-r4(12)-r2(6,3)+(-r3(7,3)-r1(2,5))*qz
      eri(3,3,2)=-r4(14)-r2(6,3)*three+(-r3(9,2)*two-r3(9,3)-r1(2,4)*two-r1(2,5))*qz+( &
&                -r2(6,4)-r2(6,5)*two)*zz-r1(2,6)*zzz
      eri(4,3,2)=-r4(8)-r3(8,2)*qx-r3(4,3)*qz-r2(2,5)*xz
      eri(5,3,2)=-r4(9)-r2(4,3)+(-r3(9,2)-r1(2,4))*qx+(-r3(5,2)-r3(5,3))*qz+(-r2(6,4) &
&                -r2(6,5))*xz-r2(4,5)*zz-r1(2,6)*xzz
      eri(6,3,2)=-r4(13)-r2(2,3)+(-r3(8,2)-r3(8,3))*qz-r2(2,5)*zz
      eri(1,1,3)=-r4(3)+r3(1,1)-r2(5,3)*three+r1(1,1)*three+(-r3(3,2)*two-r3(3,3) &
&                +r2(1,1)*two+r2(1,2)-r1(3,4)*two-r1(3,5)+r0(1)*two+r0(2))*qx+(-r2(5,4) &
&                -r2(5,5)*two+r1(1,2)+r1(1,3)*two)*xx+(-r1(3,6)+r0(3))*xxx
      eri(2,1,3)=-r4(8)+r3(4,1)-r2(5,3)+r1(1,1)+(-r3(8,3)+r2(2,2)-r1(3,5)+r0(2))*qx
      eri(3,1,3)=-r4(10)+r3(6,1)-r2(5,3)+r1(1,1)+(-r3(10,3)+r2(3,2)-r1(3,5)+r0(2))*qx+( &
&                -r3(6,2)*two+r2(5,1)*two)*qz+(-r2(3,5)*two+r1(3,3)*two)*xz+(-r2(5,4)+r1(1,2) &
&                )*zz+(-r1(3,6)+r0(3))*xzz
      eri(4,1,3)=-r4(5)+r3(2,1)-r2(6,3)+r1(2,1)+(-r3(5,2)-r3(5,3)+r2(4,1)+r2(4,2))*qx+( &
&                -r2(6,5)+r1(2,3))*xx
      eri(5,1,3)=-r4(6)+r3(3,1)-r2(3,3)+r1(3,1)+(-r3(6,2)-r3(6,3)+r2(5,1)+r2(5,2))*qx+( &
&                -r3(3,2)+r2(1,1)-r1(3,4)+r0(1))*qz+(-r2(3,5)+r1(3,3))*xx+(-r2(5,4)-r2(5,5) &
&                +r1(1,2)+r1(1,3))*xz+(-r1(3,6)+r0(3))*xxz
      eri(6,1,3)=-r4(9)+r3(5,1)+(-r3(9,3)+r2(6,2))*qx+(-r3(5,2)+r2(4,1))*qz+(-r2(6,5) &
&                +r1(2,3))*xz
      eri(1,2,3)=-r4(5)+r3(2,1)-r2(6,3)+r1(2,1)+(-r3(5,2)*two+r2(4,1)*two)*qx+(-r2(6,4) &
&                +r1(2,2))*xx
      eri(2,2,3)=-r4(12)+r3(7,1)-r2(6,3)*three+r1(2,1)*three
      eri(3,2,3)=-r4(14)+r3(9,1)-r2(6,3)+r1(2,1)+(-r3(9,2)*two+r2(6,1)*two)*qz+(-r2(6,4) &
&                +r1(2,2))*zz
      eri(4,2,3)=-r4(8)+r3(4,1)-r2(5,3)+r1(1,1)+(-r3(8,2)+r2(2,1)-r1(3,4)+r0(1))*qx
      eri(5,2,3)=-r4(9)+r3(5,1)+(-r3(9,2)+r2(6,1))*qx+(-r3(5,2)+r2(4,1))*qz+(-r2(6,4) &
&                +r1(2,2))*xz
      eri(6,2,3)=-r4(13)+r3(8,1)-r2(3,3)+r1(3,1)+(-r3(8,2)+r2(2,1)-r1(3,4)+r0(1))*qz
      eri(1,3,3)=-r4(6)+r3(3,1)-r2(3,3)+r1(3,1)+(-r3(6,2)*two+r2(5,1)*two)*qx+(-r3(3,3) &
&                +r2(1,2)-r1(3,5)+r0(2))*qz+(-r2(3,4)+r1(3,2))*xx+(-r2(5,5)*two+r1(1,3)*two) &
&                *xz+(-r1(3,6)+r0(3))*xxz
      eri(2,3,3)=-r4(13)+r3(8,1)-r2(3,3)+r1(3,1)+(-r3(8,3)+r2(2,2)-r1(3,5)+r0(2))*qz
      eri(3,3,3)=-r4(15)+r3(10,1)-r2(3,3)*three+r1(3,1)*three+(-r3(10,2)*two-r3(10,3) &
&                +r2(3,1)*two+r2(3,2)-r1(3,4)*two-r1(3,5)+r0(1)*two+r0(2))*qz+(-r2(3,4) &
&                -r2(3,5)*two+r1(3,2)+r1(3,3)*two)*zz+(-r1(3,6)+r0(3))*zzz
      eri(4,3,3)=-r4(9)+r3(5,1)+(-r3(9,2)+r2(6,1))*qx+(-r3(5,3)+r2(4,2))*qz+(-r2(6,5) &
&                +r1(2,3))*xz
      eri(5,3,3)=-r4(10)+r3(6,1)-r2(5,3)+r1(1,1)+(-r3(10,2)+r2(3,1)-r1(3,4)+r0(1))*qx+( &
&                -r3(6,2)-r3(6,3)+r2(5,1)+r2(5,2))*qz+(-r2(3,4)-r2(3,5)+r1(3,2)+r1(3,3))*xz+( &
&                -r2(5,5)+r1(1,3))*zz+(-r1(3,6)+r0(3))*xzz
      eri(6,3,3)=-r4(14)+r3(9,1)-r2(6,3)+r1(2,1)+(-r3(9,2)-r3(9,3)+r2(6,1)+r2(6,2))*qz+( &
&                -r2(6,5)+r1(2,3))*zz
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
      do k= 1,3
        do l= 1,6
          work(1)= eri(l,k,1)
          work(2)= eri(l,k,2)
          work(3)= eri(l,k,3)
          eri(l,k,1)= work(1)*rot(1,1)+work(2)*rot(2,1)+work(3)*rot(3,1)
          eri(l,k,2)= work(1)*rot(1,2)+work(2)*rot(2,2)+work(3)*rot(3,2)
          eri(l,k,3)= work(1)*rot(1,3)+work(2)*rot(2,3)+work(3)*rot(3,3)
        enddo
      enddo
      do j= 1,3
        do l= 1,6
          work(1)= eri(l,1,j)
          work(2)= eri(l,2,j)
          work(3)= eri(l,3,j)
          eri(l,1,j)= work(1)*rot(1,1)+work(2)*rot(2,1)+work(3)*rot(3,1)
          eri(l,2,j)= work(1)*rot(1,2)+work(2)*rot(2,2)+work(3)*rot(3,2)
          eri(l,3,j)= work(1)*rot(1,3)+work(2)*rot(2,3)+work(3)*rot(3,3)
        enddo
      enddo
      if(nbfijkl(4) == 6)then
        do j= 1,3
          do k= 1,3
            do l= 1,6
              phmdint(l,k,j,1)= eri(1,k,j)*rot2(1,l)+eri(2,k,j)*rot2(2,l)+eri(3,k,j)*rot2(3,l) &
&                              +eri(4,k,j)*rot2(4,l)+eri(5,k,j)*rot2(5,l)+eri(6,k,j)*rot2(6,l) 
            enddo
          enddo
        enddo
      else
        do l= 1,6
          rot3(l,1)= rot2(l,2)
          rot3(l,2)= rot2(l,5)
          rot3(l,3)= rot2(l,6)-(rot2(l,1)+rot2(l,4))*half
          rot3(l,4)= rot2(l,3)
          rot3(l,5)=(rot2(l,1)-rot2(l,4))*sqrt3h
        enddo
        do j= 1,3
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


!--------------------------------------------------------------------------
  subroutine int2dsds(phmdint,exfac1,exfac2,xyziq,xzkl,rot,nijkl,nbfijkl)
!--------------------------------------------------------------------------
!
! Calculate (ds|ds) integrals
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
      integer :: ij, kl, igrid, i, j, l, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, sqrtpi4=0.8862269254527580D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00, four=4.0D+00, five=5.0D+00
      real(8),parameter :: six=6.0D+00, seven=7.0D+00, sqrt3=1.73205080756888D+00
      real(8),parameter :: sqrt3h=8.660254037844386D-01
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*), xzkl(2), rot(3,3)
      real(8),intent(out) :: phmdint(6,6,6,6)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:4), f0(2), f1(2,4), f2(3,4), f3(4,2), f4(5)
      real(8) :: r0(4), r1(3,4), r2(6,5), r3(10,2), r4(15)
      real(8) :: ex12, ex34, ex43, ex14, ex41, expq, expq2, ex3q, ex4q, c12, c34, zip
      real(8) :: xiq, yiq, ziq, xiq2, yiq2, xyiq, xiq4, yiq4, xyiq2, xypq2, zpq, zpq2, zpq3
      real(8) :: fac, ex33q, ex3qmd, zjp, pmd, qmd, qmd2, qmd2x, qmd2y, qmd2xy, qx, qz
      real(8) :: xx, xz, zz, eri(6,6), work(6), f2w(6,4), f3w(10,2), rot2(6,6), rot3(6,5)
!
! Zero-clear
!
      r0(1:4)     = zero
      r1(1:3,1:4) = zero
      r2(1:6,1:5) = zero
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
        f0(1:2)    = zero
        f1(1:2,1:4)= zero
        f2(1:3,1:4)= zero
        f3(1:4,1:2)= zero
        f4(1:5)= zero
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
            do ii= 0,4
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
            fac= c12*ex41
            expq= expq*two*ex41*ex41
            expq2= expq*expq
            ft(0)= ft(0)*fac
            ft(1)= ft(1)*fac*expq
            ft(2)= ft(2)*fac*expq2
            ft(3)= ft(3)*fac*expq2*expq
            ft(4)= ft(4)*fac*expq2*expq2
          endif
          zpq3= zpq2*zpq
          work(1)= zjp*zjp
          work(2)= pmd
          work(3)= pmd*zjp
          work(4)= pmd*pmd
          ft(4)= ft(4)*work(4)
          f0(1)  = f0(1)  +ft(0)*work(1)
          f0(2)  = f0(2)  +ft(0)*work(2)
          do i= 1,4
            f1(1,i)= f1(1,i)-ft(1)*work(i)
            f1(2,i)= f1(2,i)-ft(1)*work(i)*zpq
          enddo
          do i= 1,4
            f2(1,i)= f2(1,i)+ft(2)*work(i)
            f2(2,i)= f2(2,i)+ft(2)*work(i)*zpq
            f2(3,i)= f2(3,i)+ft(2)*work(i)*zpq2
          enddo
          do i= 1,2
            f3(1,i)= f3(1,i)-ft(3)*work(i+2)
            f3(2,i)= f3(2,i)-ft(3)*work(i+2)*zpq
            f3(3,i)= f3(3,i)-ft(3)*work(i+2)*zpq2
            f3(4,i)= f3(4,i)-ft(3)*work(i+2)*zpq3
          enddo
          f4(1)= f4(1)+ft(4)
          f4(2)= f4(2)+ft(4)*zpq
          f4(3)= f4(3)+ft(4)*zpq2
          f4(4)= f4(4)+ft(4)*zpq3
          f4(5)= f4(5)+ft(4)*zpq3*zpq
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
        r1(1,1)= r1(1,1)+f1(1,1)*xiq*ex3qmd
        r1(2,1)= r1(2,1)+f1(1,1)*yiq*ex3qmd
        r1(3,1)= r1(3,1)+f1(2,1)    *ex3qmd
        r1(1,2)= r1(1,2)+f1(1,2)*xiq*ex3qmd
        r1(2,2)= r1(2,2)+f1(1,2)*yiq*ex3qmd
        r1(3,2)= r1(3,2)+f1(2,2)    *ex3qmd
        r1(1,3)= r1(1,3)+f1(1,3)*xiq*qmd
        r1(2,3)= r1(2,3)+f1(1,3)*yiq*qmd
        r1(3,3)= r1(3,3)+f1(2,3)    *qmd
        r1(1,4)= r1(1,4)+f1(1,3)*xiq*ex33q
        r1(2,4)= r1(2,4)+f1(1,3)*yiq*ex33q
        r1(3,4)= r1(3,4)+f1(2,3)    *ex33q
!
        do i= 1,4
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
          r2(i,4)= r2(i,4)+f2w(i,4)*qmd
          r2(i,5)= r2(i,5)+f2w(i,4)*ex33q
        enddo
! 
        do i= 1,2
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
          r3(i,2)= r3(i,2)+f3w(i,2)*ex3qmd
        enddo
!
        r4( 1)= r4( 1)+(f4(1)*xiq4 +f3(1,2)*xiq2*six         +f2(1,4)*three)*qmd2
        r4( 2)= r4( 2)+(f4(1)*xiq2 +f3(1,2)*three                          )*qmd2xy
        r4( 3)= r4( 3)+(f4(2)*xiq2 +f3(2,2)*three                          )*qmd2x
        r4( 4)= r4( 4)+(f4(1)*xyiq2+f3(1,2)*xiq2+f3(1,2)*yiq2+f2(1,4)      )*qmd2
        r4( 5)= r4( 5)+(f4(2)*xiq2 +f3(2,2)                                )*qmd2y
        r4( 6)= r4( 6)+(f4(3)*xiq2 +f3(1,2)*xiq2+f3(3,2)     +f2(1,4)      )*qmd2
        r4( 7)= r4( 7)+(f4(1)*yiq2 +f3(1,2)*three                          )*qmd2xy
        r4( 8)= r4( 8)+(f4(2)*yiq2 +f3(2,2)                                )*qmd2x
        r4( 9)= r4( 9)+(f4(3)      +f3(1,2)                                )*qmd2xy
        r4(10)= r4(10)+(f4(4)      +f3(2,2)*three                          )*qmd2x
        r4(11)= r4(11)+(f4(1)*yiq4 +f3(1,2)*yiq2*six         +f2(1,4)*three)*qmd2
        r4(12)= r4(12)+(f4(2)*yiq2 +f3(2,2)*three                          )*qmd2y
        r4(13)= r4(13)+(f4(3)*yiq2 +f3(1,2)*yiq2+f3(3,2)     +f2(1,4)      )*qmd2
        r4(14)= r4(14)+(f4(4)      +f3(2,2)*three                          )*qmd2y
        r4(15)= r4(15)+(f4(5)      +f3(3,2)*six              +f2(1,4)*three)*qmd2
      enddo
!
      qx= xzkl(1)
      qz= xzkl(2)
      xx= qx*qx
      xz= qx*qz
      zz= qz*qz
      eri(1,1)= r4(1)+r2(1,2)+r2(1,4)+r0(3)+(+r3(1,2)*two+r1(1,2)*two)*qx+(+r2(1,5)+r0(4))*xx
      eri(2,1)= r4(4)+r2(2,2)+r2(1,4)+r0(3)
      eri(3,1)= r4(6)+r2(3,2)+r2(1,4)+r0(3)+(+r3(3,2)*two+r1(3,2)*two)*qz+(+r2(1,5)+r0(4))*zz
      eri(4,1)= r4(2)+r2(4,2)+(+r3(2,2)+r1(2,2))*qx
      eri(5,1)= r4(3)+r2(5,2)+(+r3(3,2)+r1(3,2))*qx+(+r3(1,2)+r1(1,2))*qz+(+r2(1,5)+r0(4))*xz
      eri(6,1)= r4(5)+r2(6,2)+(+r3(2,2)+r1(2,2))*qz
      eri(1,2)= r4(4)+r2(1,2)+r2(2,4)+r0(3)+(+r3(4,2)*two+r1(1,2)*two)*qx+(+r2(2,5)+r0(4))*xx
      eri(2,2)= r4(11)+r2(2,2)+r2(2,4)+r0(3)
      eri(3,2)= r4(13)+r2(3,2)+r2(2,4)+r0(3)+(+r3(8,2)*two+r1(3,2)*two)*qz+(+r2(2,5)+r0(4))*zz
      eri(4,2)= r4(7)+r2(4,2)+(+r3(7,2)+r1(2,2))*qx
      eri(5,2)= r4(8)+r2(5,2)+(+r3(8,2)+r1(3,2))*qx+(+r3(4,2)+r1(1,2))*qz+(+r2(2,5)+r0(4))*xz
      eri(6,2)= r4(12)+r2(6,2)+(+r3(7,2)+r1(2,2))*qz
      eri(1,3)= r4(6)-r3(3,1)*two+r2(1,1)+r2(1,2)+r2(3,4)-r1(3,3)*two+r0(1)+r0(3)+( &
&              +r3(6,2)*two-r2(5,3)*four+r1(1,1)*two+r1(1,2)*two)*qx+(+r2(3,5)-r1(3,4)*two &
&              +r0(2)+r0(4))*xx
      eri(2,3)= r4(13)-r3(8,1)*two+r2(2,1)+r2(2,2)+r2(3,4)-r1(3,3)*two+r0(1)+r0(3)
      eri(3,3)= r4(15)-r3(10,1)*two+r2(3,1)+r2(3,2)+r2(3,4)-r1(3,3)*two+r0(1)+r0(3)+( &
&              +r3(10,2)*two-r2(3,3)*four+r1(3,1)*two+r1(3,2)*two)*qz+(+r2(3,5)-r1(3,4)*two &
&              +r0(2)+r0(4))*zz
      eri(4,3)= r4(9)-r3(5,1)*two+r2(4,1)+r2(4,2)+(+r3(9,2)-r2(6,3)*two+r1(2,1)+r1(2,2))*qx
      eri(5,3)= r4(10)-r3(6,1)*two+r2(5,1)+r2(5,2)+(+r3(10,2)-r2(3,3)*two+r1(3,1)+r1(3,2)) &
&              *qx+(+r3(6,2)-r2(5,3)*two+r1(1,1)+r1(1,2))*qz+(+r2(3,5)-r1(3,4)*two+r0(2) &
&              +r0(4))*xz
      eri(6,3)= r4(14)-r3(9,1)*two+r2(6,1)+r2(6,2)+(+r3(9,2)-r2(6,3)*two+r1(2,1)+r1(2,2))*qz
      eri(1,4)= r4(2)+r2(4,4)+r3(2,2)*two*qx+r2(4,5)*xx
      eri(2,4)= r4(7)+r2(4,4)
      eri(3,4)= r4(9)+r2(4,4)+r3(5,2)*two*qz+r2(4,5)*zz
      eri(4,4)= r4(4)+r3(4,2)*qx
      eri(5,4)= r4(5)+r3(5,2)*qx+r3(2,2)*qz+r2(4,5)*xz
      eri(6,4)= r4(8)+r3(4,2)*qz
      eri(1,5)= r4(3)-r3(1,1)+r2(5,4)-r1(1,3)+(+r3(3,2)*two-r2(1,3)*two)*qx+(+r2(5,5) &
&              -r1(1,4))*xx
      eri(2,5)= r4(8)-r3(4,1)+r2(5,4)-r1(1,3)
      eri(3,5)= r4(10)-r3(6,1)+r2(5,4)-r1(1,3)+(+r3(6,2)*two-r2(5,3)*two)*qz+(+r2(5,5) &
&              -r1(1,4))*zz
      eri(4,5)= r4(5)-r3(2,1)+(+r3(5,2)-r2(4,3))*qx
      eri(5,5)= r4(6)-r3(3,1)+(+r3(6,2)-r2(5,3))*qx+(+r3(3,2)-r2(1,3))*qz+(+r2(5,5)-r1(1,4))*xz
      eri(6,5)= r4(9)-r3(5,1)+(+r3(5,2)-r2(4,3))*qz
      eri(1,6)= r4(5)-r3(2,1)+r2(6,4)-r1(2,3)+(+r3(5,2)*two-r2(4,3)*two)*qx+(+r2(6,5) &
&              -r1(2,4))*xx
      eri(2,6)= r4(12)-r3(7,1)+r2(6,4)-r1(2,3)
      eri(3,6)= r4(14)-r3(9,1)+r2(6,4)-r1(2,3)+(+r3(9,2)*two-r2(6,3)*two)*qz+(+r2(6,5) &
&              -r1(2,4))*zz
      eri(4,6)= r4(8)-r3(4,1)+(+r3(8,2)-r2(2,3))*qx
      eri(5,6)= r4(9)-r3(5,1)+(+r3(9,2)-r2(6,3))*qx+(+r3(5,2)-r2(4,3))*qz+(+r2(6,5)-r1(2,4))*xz
      eri(6,6)= r4(13)-r3(8,1)+(+r3(8,2)-r2(2,3))*qz
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
      if(nbfijkl(2) == 6)then
        do l= 1,6
          do j= 1,6
            work(j)= eri(l,j)
          enddo
          do j= 1,6
            eri(l,j)= work(1)*rot2(1,j)+work(2)*rot2(2,j)+work(3)*rot2(3,j) &
&                    +work(4)*rot2(4,j)+work(5)*rot2(5,j)+work(6)*rot2(6,j)
          enddo
        enddo
      else
        do l= 1,6
          do j= 1,6
            work(j)= eri(l,j)
          enddo
          do j= 1,5
            eri(l,j)= work(1)*rot3(1,j)+work(2)*rot3(2,j)+work(3)*rot3(3,j) &
&                    +work(4)*rot3(4,j)+work(5)*rot3(5,j)+work(6)*rot3(6,j)
          enddo
        enddo
      endif
!
      if(nbfijkl(4) == 6)then
        do j= 1,6
          do l= 1,6
            phmdint(l,1,j,1)= eri(1,j)*rot2(1,l)+eri(2,j)*rot2(2,l)+eri(3,j)*rot2(3,l) &
&                            +eri(4,j)*rot2(4,l)+eri(5,j)*rot2(5,l)+eri(6,j)*rot2(6,l) 
          enddo
        enddo
      else
        do j= 1,6
          do l= 1,5
            phmdint(l,1,j,1)= eri(1,j)*rot3(1,l)+eri(2,j)*rot3(2,l)+eri(3,j)*rot3(3,l) &
&                            +eri(4,j)*rot3(4,l)+eri(5,j)*rot3(5,l)+eri(6,j)*rot3(6,l) 
          enddo
        enddo
      endif
!
      return
end


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
      integer :: ij, kl, igrid, i, j, l, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, sqrtpi4=0.8862269254527580D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00, four=4.0D+00, five=5.0D+00
      real(8),parameter :: six=6.0D+00, seven=7.0D+00, sqrt3=1.73205080756888D+00
      real(8),parameter :: sqrt3h=8.660254037844386D-01
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*), xzkl(2), rot(3,3)
      real(8),intent(out) :: phmdint(6,6,6,6)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:4), ftw(4,5), f0(2), f1(2,5), f2(3,5), f3(4,3), f4(5)
      real(8) :: r0(4), r1(3,6), r2(6,6), r3(10,3), r4(15) 
      real(8) :: ex12, ex34, ex43, ex14, ex41, expq, expq2, ex3q, ex4q, c12, c34, zip
      real(8) :: xiq, yiq, ziq, xiq2, yiq2, xyiq, xiq4, yiq4, xyiq2, xypq2, zpq, zpq2, zpq3, fac
      real(8) :: ex33q, ex3qmd, zjp, pmd, qmd, qmd2, qmd2x, qmd2y, qmd2xy, qx, qz, xx, xz, zz
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
            do ii= 0,4
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
            fac= c12*ex41
            expq= expq*two*ex41*ex41
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
          rot3(l,1)= rot2(l,2)
          rot3(l,2)= rot2(l,5)
          rot3(l,3)= rot2(l,6)-(rot2(l,1)+rot2(l,4))*half
          rot3(l,4)= rot2(l,3)
          rot3(l,5)=(rot2(l,1)-rot2(l,4))*sqrt3h
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
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:5)
      real(8) :: f0, f1(2,2), f2(3,2), f3(4,2), f4(5,2), f5(6), ftw(5,2)
      real(8) :: r0(5), r1(3,9), r2(6,8), r3(10,6), r4(15,3), r5(21)
      real(8) :: ex12, ex34, ex43, ex14, ex41, expq, expq2, expq4, ex3q, ex4q, c12, c34, zip
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
            do ii= 0,5
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
            fac= c12*ex41
            expq= expq*two*ex41*ex41
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
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:5)
      real(8) :: f0(2), f1(2,4), f2(3,4), f3(4,4), f4(5,2), f5(6), ftw(5,4)
      real(8) :: r0(6), r1(3,9), r2(6,10), r3(10,7), r4(15,3), r5(21)
      real(8) :: ex12, ex34, ex43, ex14, ex41, expq, expq2, expq4, ex3q, ex4q, c12, c34, zip
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
            do ii= 0,5
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
            fac= c12*ex41
            expq= expq*two*ex41*ex41
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


!--------------------------------------------------------------------------
  subroutine int2dppp(phmdint,exfac1,exfac2,xyziq,xzkl,rot,nijkl,nbfijkl)
!--------------------------------------------------------------------------
!
! Calculate (dp|pp) integrals
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
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:5)
      real(8) :: f0(2), f1(2,5), f2(3,5), f3(4,5), f4(5,3), f5(6), ftw(5,5)
      real(8) :: r0(6), r1(3,12), r2(6,13), r3(10,9), r4(15,4), r5(21)
      real(8) :: ex12, ex34, ex43, ex14, ex41, expq, expq2, expq4, ex3q, ex4q, c12, c34, zip
      real(8) :: xiq, yiq, ziq, xiq2, yiq2, xyiq, xiq4, yiq4, xyiq2, xypq2, zpq, zpq2, fac
      real(8) :: ex33q, ex34q, zjp, pmd, qmd, qmd2, qmd3, qmd3x, qmd3y, qmd3xy
      real(8) :: qx, qz, xx, xz, zz, xxx, xxz, xzz, zzz, eri(6,3,3,3), work(8)
      real(8) :: f1w(3,4), f2w(6,5), f3w(10,5), f4w(15,3), rot2(6,6), rot3(6,5)
!
! Zero-clear
!
      r0(1:6)     = zero
      r1(1:3,1:12)= zero
      r2(1:6,1:13)= zero
      r3(1:10,1:9)= zero
      r4(1:15,1:4)= zero
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
        f1(1:2,1:5)= zero
        f2(1:3,1:5)= zero
        f3(1:4,1:5)= zero
        f4(1:5,1:3)= zero
        f5(1:6)    = zero
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
            do ii= 0,5
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
            fac= c12*ex41
            expq= expq*two*ex41*ex41
            expq2= expq*expq
            expq4= expq2*expq2
            ft(0)= ft(0)*fac
            ft(1)= ft(1)*fac*expq
            ft(2)= ft(2)*fac*expq2
            ft(3)= ft(3)*fac*expq2*expq
            ft(4)= ft(4)*fac*expq4
            ft(5)= ft(5)*fac*expq4*expq
          endif
          ftw(1,1)= zjp*zip
          ftw(1,2)= pmd
          ftw(1,3)= pmd*zjp
          ftw(1,4)= pmd*zip
          ftw(1,5)= pmd*pmd
          do i= 1,5
            ftw(2,i)= ftw(1,i)*zpq
            ftw(3,i)= ftw(1,i)*zpq2
            ftw(4,i)= ftw(1,i)*zpq2*zpq
            ftw(5,i)= ftw(1,i)*zpq2*zpq2
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
          do i= 1,5
            do j= 1,4
              f3(j,i)= f3(j,i)-ft(3)*ftw(j,i)
            enddo
          enddo
          do i= 1,3
            do j= 1,5
              f4(j,i)= f4(j,i)+ft(4)*ftw(j,i+2)
            enddo
          enddo
          f5(1)= f5(1)-ft(5)*ftw(1,5)
          f5(2)= f5(2)-ft(5)*ftw(2,5)
          f5(3)= f5(3)-ft(5)*ftw(3,5)
          f5(4)= f5(4)-ft(5)*ftw(4,5)
          f5(5)= f5(5)-ft(5)*ftw(5,5)
          f5(6)= f5(6)-ft(5)*ftw(5,5)*zpq
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
        do i= 1,4
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
        do i= 10,12
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,4)*work(i-9)
          enddo
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
            r2(j,i)= r2(j,i)+f2w(j,4)*work(i-4)
          enddo
        enddo
        do i= 11,13
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,5)*work(i-10)
          enddo
        enddo
! 
        do i= 1,5
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
        do i= 5,6
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,4)*work(i+2)
          enddo
        enddo
        do i= 7,9
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,5)*work(i-3)
          enddo
        enddo
!
        do i= 1,3
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
          r4(j,2)= r4(j,2)+f4w(j,2)*qmd3
          r4(j,3)= r4(j,3)+f4w(j,3)*work(7)
          r4(j,4)= r4(j,4)+f4w(j,3)*work(8)
        enddo
!
        r5( 1)= r5( 1)+(f5(1)*xiq4 +f4(1,3)*xiq2*ten               +f3(1,5)*p15  )*qmd3x
        r5( 2)= r5( 2)+(f5(1)*xiq4 +f4(1,3)*xiq2*six               +f3(1,5)*three)*qmd3y
        r5( 3)= r5( 3)+(f5(2)*xiq4 +f4(2,3)*xiq2*six               +f3(2,5)*three)*qmd3
        r5( 4)= r5( 4)+(f5(1)*xyiq2+f4(1,3)*xiq2+f4(1,3)*yiq2*three+f3(1,5)*three)*qmd3x
        r5( 5)= r5( 5)+(f5(2)*xiq2 +f4(2,3)*three                                )*qmd3xy
        r5( 6)= r5( 6)+(f5(3)*xiq2 +f4(1,3)*xiq2+f4(3,3)*three     +f3(1,5)*three)*qmd3x
        r5( 7)= r5( 7)+(f5(1)*xyiq2+f4(1,3)*xiq2*three+f4(1,3)*yiq2+f3(1,5)*three)*qmd3y
        r5( 8)= r5( 8)+(f5(2)*xyiq2+f4(2,3)*xiq2+f4(2,3)*yiq2      +f3(2,5)      )*qmd3
        r5( 9)= r5( 9)+(f5(3)*xiq2 +f4(1,3)*xiq2+f4(3,3)           +f3(1,5)      )*qmd3y
        r5(10)= r5(10)+(f5(4)*xiq2 +f4(2,3)*xiq2*three+f4(4,3)     +f3(2,5)*three)*qmd3
        r5(11)= r5(11)+(f5(1)*yiq4 +f4(1,3)*yiq2*six               +f3(1,5)*three)*qmd3x
        r5(12)= r5(12)+(f5(2)*yiq2 +f4(2,3)*three                                )*qmd3xy
        r5(13)= r5(13)+(f5(3)*yiq2 +f4(1,3)*yiq2+f4(3,3)           +f3(1,5)      )*qmd3x
        r5(14)= r5(14)+(f5(4)      +f4(2,3)*three                                )*qmd3xy
        r5(15)= r5(15)+(f5(5)      +f4(3,3)*six                    +f3(1,5)*three)*qmd3x
        r5(16)= r5(16)+(f5(1)*yiq4 +f4(1,3)*yiq2*ten               +f3(1,5)*p15  )*qmd3y
        r5(17)= r5(17)+(f5(2)*yiq4 +f4(2,3)*yiq2*six               +f3(2,5)*three)*qmd3
        r5(18)= r5(18)+(f5(3)*yiq2 +f4(1,3)*yiq2+f4(3,3)*three     +f3(1,5)*three)*qmd3y
        r5(19)= r5(19)+(f5(4)*yiq2 +f4(2,3)*yiq2*three+f4(4,3)     +f3(2,5)*three)*qmd3
        r5(20)= r5(20)+(f5(5)      +f4(3,3)*six                    +f3(1,5)*three)*qmd3y
        r5(21)= r5(21)+(f5(6)      +f4(4,3)*ten                    +f3(2,5)*p15  )*qmd3
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
      eri(1,1,1,1)= r5(1)+r3(1,2)+r3(1,7)*three+r1(1,4)*three+(+r4(1,3)*two+r4(1,4) &
&                  +r2(1,3)*two+r2(1,4)+r2(1,11)*two+r2(1,12)+r0(4)*two+r0(5))*qx+(+r3(1,8) &
&                  +r3(1,9)*two+r1(1,5)+r1(1,6)*two)*xx+(+r2(1,13)+r0(6))*xxx
      eri(2,1,1,1)= r5(4)+r3(4,2)+r3(1,7)+r1(1,4)+(+r4(4,4)+r2(2,4)+r2(1,12)+r0(5))*qx
      eri(3,1,1,1)= r5(6)+r3(6,2)+r3(1,7)+r1(1,4)+(+r4(6,4)+r2(3,4)+r2(1,12)+r0(5))*qx+( &
&                  +r4(3,3)*two+r2(5,3)*two)*qz+(+r3(3,9)*two+r1(3,6)*two)*xz+(+r3(1,8) &
&                  +r1(1,5))*zz+(+r2(1,13)+r0(6))*xzz
      eri(4,1,1,1)= r5(2)+r3(2,2)+r3(2,7)+r1(2,4)+(+r4(2,3)+r4(2,4)+r2(4,3)+r2(4,4))*qx+( &
&                  +r3(2,9)+r1(2,6))*xx
      eri(5,1,1,1)= r5(3)+r3(3,2)+r3(3,7)+r1(3,4)+(+r4(3,3)+r4(3,4)+r2(5,3)+r2(5,4))*qx+( &
&                  +r4(1,3)+r2(1,3)+r2(1,11)+r0(4))*qz+(+r3(3,9)+r1(3,6))*xx+(+r3(1,8) &
&                  +r3(1,9)+r1(1,5)+r1(1,6))*xz+(+r2(1,13)+r0(6))*xxz
      eri(6,1,1,1)= r5(5)+r3(5,2)+(+r4(5,4)+r2(6,4))*qx+(+r4(2,3)+r2(4,3))*qz+(+r3(2,9) &
&                  +r1(2,6))*xz
      eri(1,2,1,1)= r5(2)+r3(2,2)+r3(2,7)+r1(2,4)+(+r4(2,3)*two+r2(4,3)*two)*qx+(+r3(2,8) &
&                  +r1(2,5))*xx
      eri(2,2,1,1)= r5(7)+r3(7,2)+r3(2,7)*three+r1(2,4)*three
      eri(3,2,1,1)= r5(9)+r3(9,2)+r3(2,7)+r1(2,4)+(+r4(5,3)*two+r2(6,3)*two)*qz+(+r3(2,8) &
&                  +r1(2,5))*zz
      eri(4,2,1,1)= r5(4)+r3(4,2)+r3(1,7)+r1(1,4)+(+r4(4,3)+r2(2,3)+r2(1,11)+r0(4))*qx
      eri(5,2,1,1)= r5(5)+r3(5,2)+(+r4(5,3)+r2(6,3))*qx+(+r4(2,3)+r2(4,3))*qz+(+r3(2,8) &
&                  +r1(2,5))*xz
      eri(6,2,1,1)= r5(8)+r3(8,2)+r3(3,7)+r1(3,4)+(+r4(4,3)+r2(2,3)+r2(1,11)+r0(4))*qz
      eri(1,3,1,1)= r5(3)+r3(3,2)+r3(3,7)+r1(3,4)+(+r4(3,3)*two+r2(5,3)*two)*qx+(+r4(1,4) &
&                  +r2(1,4)+r2(1,12)+r0(5))*qz+(+r3(3,8)+r1(3,5))*xx+(+r3(1,9)*two &
&                  +r1(1,6)*two)*xz+(+r2(1,13)+r0(6))*xxz
      eri(2,3,1,1)= r5(8)+r3(8,2)+r3(3,7)+r1(3,4)+(+r4(4,4)+r2(2,4)+r2(1,12)+r0(5))*qz
      eri(3,3,1,1)= r5(10)+r3(10,2)+r3(3,7)*three+r1(3,4)*three+(+r4(6,3)*two+r4(6,4) &
&                  +r2(3,3)*two+r2(3,4)+r2(1,11)*two+r2(1,12)+r0(4)*two+r0(5))*qz+(+r3(3,8) &
&                  +r3(3,9)*two+r1(3,5)+r1(3,6)*two)*zz+(+r2(1,13)+r0(6))*zzz
      eri(4,3,1,1)= r5(5)+r3(5,2)+(+r4(5,3)+r2(6,3))*qx+(+r4(2,4)+r2(4,4))*qz+(+r3(2,9) &
&                  +r1(2,6))*xz
      eri(5,3,1,1)= r5(6)+r3(6,2)+r3(1,7)+r1(1,4)+(+r4(6,3)+r2(3,3)+r2(1,11)+r0(4))*qx+( &
&                  +r4(3,3)+r4(3,4)+r2(5,3)+r2(5,4))*qz+(+r3(3,8)+r3(3,9)+r1(3,5)+r1(3,6))*xz &
&                  +(+r3(1,9)+r1(1,6))*zz+(+r2(1,13)+r0(6))*xzz
      eri(6,3,1,1)= r5(9)+r3(9,2)+r3(2,7)+r1(2,4)+(+r4(5,3)+r4(5,4)+r2(6,3)+r2(6,4))*qz+( &
&                  +r3(2,9)+r1(2,6))*zz
      eri(1,1,2,1)= r5(2)+r3(2,7)*three+(+r4(2,3)*two+r4(2,4)+r2(4,11)*two+r2(4,12))*qx+( &
&                  +r3(2,8)+r3(2,9)*two)*xx+r2(4,13)*xxx
      eri(2,1,2,1)= r5(7)+r3(2,7)+(+r4(7,4)+r2(4,12))*qx
      eri(3,1,2,1)= r5(9)+r3(2,7)+(+r4(9,4)+r2(4,12))*qx+r4(5,3)*two*qz+r3(5,9)*two*xz &
&                  +r3(2,8)*zz+r2(4,13)*xzz
      eri(4,1,2,1)= r5(4)+r3(4,7)+(+r4(4,3)+r4(4,4))*qx+r3(4,9)*xx
      eri(5,1,2,1)= r5(5)+r3(5,7)+(+r4(5,3)+r4(5,4))*qx+(+r4(2,3)+r2(4,11))*qz+r3(5,9)*xx &
&                  +(+r3(2,8)+r3(2,9))*xz+r2(4,13)*xxz
      eri(6,1,2,1)= r5(8)+r4(8,4)*qx+r4(4,3)*qz+r3(4,9)*xz
      eri(1,2,2,1)= r5(4)+r3(4,7)+r4(4,3)*two*qx+r3(4,8)*xx
      eri(2,2,2,1)= r5(11)+r3(4,7)*three
      eri(3,2,2,1)= r5(13)+r3(4,7)+r4(8,3)*two*qz+r3(4,8)*zz
      eri(4,2,2,1)= r5(7)+r3(2,7)+(+r4(7,3)+r2(4,11))*qx
      eri(5,2,2,1)= r5(8)+r4(8,3)*qx+r4(4,3)*qz+r3(4,8)*xz
      eri(6,2,2,1)= r5(12)+r3(5,7)+(+r4(7,3)+r2(4,11))*qz
      eri(1,3,2,1)= r5(5)+r3(5,7)+r4(5,3)*two*qx+(+r4(2,4)+r2(4,12))*qz+r3(5,8)*xx &
&                  +r3(2,9)*two*xz+r2(4,13)*xxz
      eri(2,3,2,1)= r5(12)+r3(5,7)+(+r4(7,4)+r2(4,12))*qz
      eri(3,3,2,1)= r5(14)+r3(5,7)*three+(+r4(9,3)*two+r4(9,4)+r2(4,11)*two+r2(4,12))*qz+( &
&                  +r3(5,8)+r3(5,9)*two)*zz+r2(4,13)*zzz
      eri(4,3,2,1)= r5(8)+r4(8,3)*qx+r4(4,4)*qz+r3(4,9)*xz
      eri(5,3,2,1)= r5(9)+r3(2,7)+(+r4(9,3)+r2(4,11))*qx+(+r4(5,3)+r4(5,4))*qz+(+r3(5,8) &
&                  +r3(5,9))*xz+r3(2,9)*zz+r2(4,13)*xzz
      eri(6,3,2,1)= r5(13)+r3(4,7)+(+r4(8,3)+r4(8,4))*qz+r3(4,9)*zz
      eri(1,1,3,1)= r5(3)-r4(1,1)+r3(3,7)*three-r2(1,5)*three+(+r4(3,3)*two+r4(3,4) &
&                  -r3(1,3)*two-r3(1,4)+r2(5,11)*two+r2(5,12)-r1(1,7)*two-r1(1,8))*qx+( &
&                  +r3(3,8)+r3(3,9)*two-r2(1,6)-r2(1,7)*two)*xx+(+r2(5,13)-r1(1,9))*xxx
      eri(2,1,3,1)= r5(8)-r4(4,1)+r3(3,7)-r2(1,5)+(+r4(8,4)-r3(4,4)+r2(5,12)-r1(1,8))*qx
      eri(3,1,3,1)= r5(10)-r4(6,1)+r3(3,7)-r2(1,5)+(+r4(10,4)-r3(6,4)+r2(5,12)-r1(1,8))*qx &
&                  +(+r4(6,3)*two-r3(3,3)*two)*qz+(+r3(6,9)*two-r2(5,7)*two)*xz+(+r3(3,8) &
&                  -r2(1,6))*zz+(+r2(5,13)-r1(1,9))*xzz
      eri(4,1,3,1)= r5(5)-r4(2,1)+r3(5,7)-r2(4,5)+(+r4(5,3)+r4(5,4)-r3(2,3)-r3(2,4))*qx+( &
&                  +r3(5,9)-r2(4,7))*xx
      eri(5,1,3,1)= r5(6)-r4(3,1)+r3(6,7)-r2(5,5)+(+r4(6,3)+r4(6,4)-r3(3,3)-r3(3,4))*qx+( &
&                  +r4(3,3)-r3(1,3)+r2(5,11)-r1(1,7))*qz+(+r3(6,9)-r2(5,7))*xx+(+r3(3,8) &
&                  +r3(3,9)-r2(1,6)-r2(1,7))*xz+(+r2(5,13)-r1(1,9))*xxz
      eri(6,1,3,1)= r5(9)-r4(5,1)+(+r4(9,4)-r3(5,4))*qx+(+r4(5,3)-r3(2,3))*qz+(+r3(5,9) &
&                  -r2(4,7))*xz
      eri(1,2,3,1)= r5(5)-r4(2,1)+r3(5,7)-r2(4,5)+(+r4(5,3)*two-r3(2,3)*two)*qx+(+r3(5,8) &
&                  -r2(4,6))*xx
      eri(2,2,3,1)= r5(12)-r4(7,1)+r3(5,7)*three-r2(4,5)*three
      eri(3,2,3,1)= r5(14)-r4(9,1)+r3(5,7)-r2(4,5)+(+r4(9,3)*two-r3(5,3)*two)*qz+(+r3(5,8) &
&                  -r2(4,6))*zz
      eri(4,2,3,1)= r5(8)-r4(4,1)+r3(3,7)-r2(1,5)+(+r4(8,3)-r3(4,3)+r2(5,11)-r1(1,7))*qx
      eri(5,2,3,1)= r5(9)-r4(5,1)+(+r4(9,3)-r3(5,3))*qx+(+r4(5,3)-r3(2,3))*qz+(+r3(5,8) &
&                  -r2(4,6))*xz
      eri(6,2,3,1)= r5(13)-r4(8,1)+r3(6,7)-r2(5,5)+(+r4(8,3)-r3(4,3)+r2(5,11)-r1(1,7))*qz
      eri(1,3,3,1)= r5(6)-r4(3,1)+r3(6,7)-r2(5,5)+(+r4(6,3)*two-r3(3,3)*two)*qx+(+r4(3,4) &
&                  -r3(1,4)+r2(5,12)-r1(1,8))*qz+(+r3(6,8)-r2(5,6))*xx+(+r3(3,9)*two &
&                  -r2(1,7)*two)*xz+(+r2(5,13)-r1(1,9))*xxz
      eri(2,3,3,1)= r5(13)-r4(8,1)+r3(6,7)-r2(5,5)+(+r4(8,4)-r3(4,4)+r2(5,12)-r1(1,8))*qz
      eri(3,3,3,1)= r5(15)-r4(10,1)+r3(6,7)*three-r2(5,5)*three+(+r4(10,3)*two+r4(10,4) &
&                  -r3(6,3)*two-r3(6,4)+r2(5,11)*two+r2(5,12)-r1(1,7)*two-r1(1,8))*qz+( &
&                  +r3(6,8)+r3(6,9)*two-r2(5,6)-r2(5,7)*two)*zz+(+r2(5,13)-r1(1,9))*zzz
      eri(4,3,3,1)= r5(9)-r4(5,1)+(+r4(9,3)-r3(5,3))*qx+(+r4(5,4)-r3(2,4))*qz+(+r3(5,9) &
&                  -r2(4,7))*xz
      eri(5,3,3,1)= r5(10)-r4(6,1)+r3(3,7)-r2(1,5)+(+r4(10,3)-r3(6,3)+r2(5,11)-r1(1,7))*qx &
&                  +(+r4(6,3)+r4(6,4)-r3(3,3)-r3(3,4))*qz+(+r3(6,8)+r3(6,9)-r2(5,6)-r2(5,7)) &
&                  *xz+(+r3(3,9)-r2(1,7))*zz+(+r2(5,13)-r1(1,9))*xzz
      eri(6,3,3,1)= r5(14)-r4(9,1)+r3(5,7)-r2(4,5)+(+r4(9,3)+r4(9,4)-r3(5,3)-r3(5,4))*qz+( &
&                  +r3(5,9)-r2(4,7))*zz
      eri(1,1,1,2)= eri(1,1,2,1)
      eri(2,1,1,2)= eri(2,1,2,1)
      eri(3,1,1,2)= eri(3,1,2,1)
      eri(4,1,1,2)= eri(4,1,2,1)
      eri(5,1,1,2)= eri(5,1,2,1)
      eri(6,1,1,2)= eri(6,1,2,1)
      eri(1,2,1,2)= eri(1,2,2,1)
      eri(2,2,1,2)= eri(2,2,2,1)
      eri(3,2,1,2)= eri(3,2,2,1)
      eri(4,2,1,2)= eri(4,2,2,1)
      eri(5,2,1,2)= eri(5,2,2,1)
      eri(6,2,1,2)= eri(6,2,2,1)
      eri(1,3,1,2)= eri(1,3,2,1)
      eri(2,3,1,2)= eri(2,3,2,1)
      eri(3,3,1,2)= eri(3,3,2,1)
      eri(4,3,1,2)= eri(4,3,2,1)
      eri(5,3,1,2)= eri(5,3,2,1)
      eri(6,3,1,2)= eri(6,3,2,1)
      eri(1,1,2,2)= r5(4)+r3(1,2)+r3(4,7)*three+r1(1,4)*three+(+r4(4,3)*two+r4(4,4) &
&                  +r2(1,3)*two+r2(1,4)+r2(2,11)*two+r2(2,12)+r0(4)*two+r0(5))*qx+(+r3(4,8) &
&                  +r3(4,9)*two+r1(1,5)+r1(1,6)*two)*xx+(+r2(2,13)+r0(6))*xxx
      eri(2,1,2,2)= r5(11)+r3(4,2)+r3(4,7)+r1(1,4)+(+r4(11,4)+r2(2,4)+r2(2,12)+r0(5))*qx
      eri(3,1,2,2)= r5(13)+r3(6,2)+r3(4,7)+r1(1,4)+(+r4(13,4)+r2(3,4)+r2(2,12)+r0(5))*qx+( &
&                  +r4(8,3)*two+r2(5,3)*two)*qz+(+r3(8,9)*two+r1(3,6)*two)*xz+(+r3(4,8) &
&                  +r1(1,5))*zz+(+r2(2,13)+r0(6))*xzz
      eri(4,1,2,2)= r5(7)+r3(2,2)+r3(7,7)+r1(2,4)+(+r4(7,3)+r4(7,4)+r2(4,3)+r2(4,4))*qx+( &
&                  +r3(7,9)+r1(2,6))*xx
      eri(5,1,2,2)= r5(8)+r3(3,2)+r3(8,7)+r1(3,4)+(+r4(8,3)+r4(8,4)+r2(5,3)+r2(5,4))*qx+( &
&                  +r4(4,3)+r2(1,3)+r2(2,11)+r0(4))*qz+(+r3(8,9)+r1(3,6))*xx+(+r3(4,8) &
&                  +r3(4,9)+r1(1,5)+r1(1,6))*xz+(+r2(2,13)+r0(6))*xxz
      eri(6,1,2,2)= r5(12)+r3(5,2)+(+r4(12,4)+r2(6,4))*qx+(+r4(7,3)+r2(4,3))*qz+(+r3(7,9) &
&                  +r1(2,6))*xz
      eri(1,2,2,2)= r5(7)+r3(2,2)+r3(7,7)+r1(2,4)+(+r4(7,3)*two+r2(4,3)*two)*qx+(+r3(7,8) &
&                  +r1(2,5))*xx
      eri(2,2,2,2)= r5(16)+r3(7,2)+r3(7,7)*three+r1(2,4)*three
      eri(3,2,2,2)= r5(18)+r3(9,2)+r3(7,7)+r1(2,4)+(+r4(12,3)*two+r2(6,3)*two)*qz+( &
&                  +r3(7,8)+r1(2,5))*zz
      eri(4,2,2,2)= r5(11)+r3(4,2)+r3(4,7)+r1(1,4)+(+r4(11,3)+r2(2,3)+r2(2,11)+r0(4))*qx
      eri(5,2,2,2)= r5(12)+r3(5,2)+(+r4(12,3)+r2(6,3))*qx+(+r4(7,3)+r2(4,3))*qz+(+r3(7,8) &
&                  +r1(2,5))*xz
      eri(6,2,2,2)= r5(17)+r3(8,2)+r3(8,7)+r1(3,4)+(+r4(11,3)+r2(2,3)+r2(2,11)+r0(4))*qz
      eri(1,3,2,2)= r5(8)+r3(3,2)+r3(8,7)+r1(3,4)+(+r4(8,3)*two+r2(5,3)*two)*qx+(+r4(4,4) &
&                  +r2(1,4)+r2(2,12)+r0(5))*qz+(+r3(8,8)+r1(3,5))*xx+(+r3(4,9)*two &
&                  +r1(1,6)*two)*xz+(+r2(2,13)+r0(6))*xxz
      eri(2,3,2,2)= r5(17)+r3(8,2)+r3(8,7)+r1(3,4)+(+r4(11,4)+r2(2,4)+r2(2,12)+r0(5))*qz
      eri(3,3,2,2)= r5(19)+r3(10,2)+r3(8,7)*three+r1(3,4)*three+(+r4(13,3)*two+r4(13,4) &
&                  +r2(3,3)*two+r2(3,4)+r2(2,11)*two+r2(2,12)+r0(4)*two+r0(5))*qz+(+r3(8,8) &
&                  +r3(8,9)*two+r1(3,5)+r1(3,6)*two)*zz+(+r2(2,13)+r0(6))*zzz
      eri(4,3,2,2)= r5(12)+r3(5,2)+(+r4(12,3)+r2(6,3))*qx+(+r4(7,4)+r2(4,4))*qz+(+r3(7,9) &
&                  +r1(2,6))*xz
      eri(5,3,2,2)= r5(13)+r3(6,2)+r3(4,7)+r1(1,4)+(+r4(13,3)+r2(3,3)+r2(2,11)+r0(4))*qx+( &
&                  +r4(8,3)+r4(8,4)+r2(5,3)+r2(5,4))*qz+(+r3(8,8)+r3(8,9)+r1(3,5)+r1(3,6))*xz &
&                  +(+r3(4,9)+r1(1,6))*zz+(+r2(2,13)+r0(6))*xzz
      eri(6,3,2,2)= r5(18)+r3(9,2)+r3(7,7)+r1(2,4)+(+r4(12,3)+r4(12,4)+r2(6,3)+r2(6,4))*qz &
&                  +(+r3(7,9)+r1(2,6))*zz
      eri(1,1,3,2)= r5(5)-r4(2,1)+r3(5,7)*three-r2(4,5)*three+(+r4(5,3)*two+r4(5,4) &
&                  -r3(2,3)*two-r3(2,4)+r2(6,11)*two+r2(6,12)-r1(2,7)*two-r1(2,8))*qx+( &
&                  +r3(5,8)+r3(5,9)*two-r2(4,6)-r2(4,7)*two)*xx+(+r2(6,13)-r1(2,9))*xxx
      eri(2,1,3,2)= r5(12)-r4(7,1)+r3(5,7)-r2(4,5)+(+r4(12,4)-r3(7,4)+r2(6,12)-r1(2,8))*qx
      eri(3,1,3,2)= r5(14)-r4(9,1)+r3(5,7)-r2(4,5)+(+r4(14,4)-r3(9,4)+r2(6,12)-r1(2,8))*qx &
&                  +(+r4(9,3)*two-r3(5,3)*two)*qz+(+r3(9,9)*two-r2(6,7)*two)*xz+(+r3(5,8) &
&                  -r2(4,6))*zz+(+r2(6,13)-r1(2,9))*xzz
      eri(4,1,3,2)= r5(8)-r4(4,1)+r3(8,7)-r2(2,5)+(+r4(8,3)+r4(8,4)-r3(4,3)-r3(4,4))*qx+( &
&                  +r3(8,9)-r2(2,7))*xx
      eri(5,1,3,2)= r5(9)-r4(5,1)+r3(9,7)-r2(6,5)+(+r4(9,3)+r4(9,4)-r3(5,3)-r3(5,4))*qx+( &
&                  +r4(5,3)-r3(2,3)+r2(6,11)-r1(2,7))*qz+(+r3(9,9)-r2(6,7))*xx+(+r3(5,8) &
&                  +r3(5,9)-r2(4,6)-r2(4,7))*xz+(+r2(6,13)-r1(2,9))*xxz
      eri(6,1,3,2)= r5(13)-r4(8,1)+(+r4(13,4)-r3(8,4))*qx+(+r4(8,3)-r3(4,3))*qz+(+r3(8,9) &
&                  -r2(2,7))*xz
      eri(1,2,3,2)= r5(8)-r4(4,1)+r3(8,7)-r2(2,5)+(+r4(8,3)*two-r3(4,3)*two)*qx+(+r3(8,8) &
&                  -r2(2,6))*xx
      eri(2,2,3,2)= r5(17)-r4(11,1)+r3(8,7)*three-r2(2,5)*three
      eri(3,2,3,2)= r5(19)-r4(13,1)+r3(8,7)-r2(2,5)+(+r4(13,3)*two-r3(8,3)*two)*qz+( &
&                  +r3(8,8)-r2(2,6))*zz
      eri(4,2,3,2)= r5(12)-r4(7,1)+r3(5,7)-r2(4,5)+(+r4(12,3)-r3(7,3)+r2(6,11)-r1(2,7))*qx
      eri(5,2,3,2)= r5(13)-r4(8,1)+(+r4(13,3)-r3(8,3))*qx+(+r4(8,3)-r3(4,3))*qz+(+r3(8,8) &
&                  -r2(2,6))*xz
      eri(6,2,3,2)= r5(18)-r4(12,1)+r3(9,7)-r2(6,5)+(+r4(12,3)-r3(7,3)+r2(6,11)-r1(2,7)) &
&                  *qz
      eri(1,3,3,2)= r5(9)-r4(5,1)+r3(9,7)-r2(6,5)+(+r4(9,3)*two-r3(5,3)*two)*qx+(+r4(5,4) &
&                  -r3(2,4)+r2(6,12)-r1(2,8))*qz+(+r3(9,8)-r2(6,6))*xx+(+r3(5,9)*two &
&                  -r2(4,7)*two)*xz+(+r2(6,13)-r1(2,9))*xxz
      eri(2,3,3,2)= r5(18)-r4(12,1)+r3(9,7)-r2(6,5)+(+r4(12,4)-r3(7,4)+r2(6,12)-r1(2,8)) &
&                  *qz
      eri(3,3,3,2)= r5(20)-r4(14,1)+r3(9,7)*three-r2(6,5)*three+(+r4(14,3)*two+r4(14,4) &
&                  -r3(9,3)*two-r3(9,4)+r2(6,11)*two+r2(6,12)-r1(2,7)*two-r1(2,8))*qz+( &
&                  +r3(9,8)+r3(9,9)*two-r2(6,6)-r2(6,7)*two)*zz+(+r2(6,13)-r1(2,9))*zzz
      eri(4,3,3,2)= r5(13)-r4(8,1)+(+r4(13,3)-r3(8,3))*qx+(+r4(8,4)-r3(4,4))*qz+(+r3(8,9) &
&                  -r2(2,7))*xz
      eri(5,3,3,2)= r5(14)-r4(9,1)+r3(5,7)-r2(4,5)+(+r4(14,3)-r3(9,3)+r2(6,11)-r1(2,7))*qx &
&                  +(+r4(9,3)+r4(9,4)-r3(5,3)-r3(5,4))*qz+(+r3(9,8)+r3(9,9)-r2(6,6)-r2(6,7)) &
&                  *xz+(+r3(5,9)-r2(4,7))*zz+(+r2(6,13)-r1(2,9))*xzz
      eri(6,3,3,2)= r5(19)-r4(13,1)+r3(8,7)-r2(2,5)+(+r4(13,3)+r4(13,4)-r3(8,3)-r3(8,4)) &
&                  *qz+(+r3(8,9)-r2(2,7))*zz
      eri(1,1,1,3)= r5(3)-r4(1,2)+r3(3,7)*three-r2(1,8)*three+(+r4(3,3)*two+r4(3,4) &
&                  -r3(1,5)*two-r3(1,6)+r2(5,11)*two+r2(5,12)-r1(1,10)*two-r1(1,11))*qx+( &
&                  +r3(3,8)+r3(3,9)*two-r2(1,9)-r2(1,10)*two)*xx+(+r2(5,13)-r1(1,12))*xxx
      eri(2,1,1,3)= r5(8)-r4(4,2)+r3(3,7)-r2(1,8)+(+r4(8,4)-r3(4,6)+r2(5,12)-r1(1,11))*qx
      eri(3,1,1,3)= r5(10)-r4(6,2)+r3(3,7)-r2(1,8)+(+r4(10,4)-r3(6,6)+r2(5,12)-r1(1,11)) &
&                  *qx+(+r4(6,3)*two-r3(3,5)*two)*qz+(+r3(6,9)*two-r2(5,10)*two)*xz+(+r3(3,8) &
&                  -r2(1,9))*zz+(+r2(5,13)-r1(1,12))*xzz
      eri(4,1,1,3)= r5(5)-r4(2,2)+r3(5,7)-r2(4,8)+(+r4(5,3)+r4(5,4)-r3(2,5)-r3(2,6))*qx+( &
&                  +r3(5,9)-r2(4,10))*xx
      eri(5,1,1,3)= r5(6)-r4(3,2)+r3(6,7)-r2(5,8)+(+r4(6,3)+r4(6,4)-r3(3,5)-r3(3,6))*qx+( &
&                  +r4(3,3)-r3(1,5)+r2(5,11)-r1(1,10))*qz+(+r3(6,9)-r2(5,10))*xx+(+r3(3,8) &
&                  +r3(3,9)-r2(1,9)-r2(1,10))*xz+(+r2(5,13)-r1(1,12))*xxz
      eri(6,1,1,3)= r5(9)-r4(5,2)+(+r4(9,4)-r3(5,6))*qx+(+r4(5,3)-r3(2,5))*qz+(+r3(5,9) &
&                  -r2(4,10))*xz
      eri(1,2,1,3)= r5(5)-r4(2,2)+r3(5,7)-r2(4,8)+(+r4(5,3)*two-r3(2,5)*two)*qx+(+r3(5,8) &
&                  -r2(4,9))*xx
      eri(2,2,1,3)= r5(12)-r4(7,2)+r3(5,7)*three-r2(4,8)*three
      eri(3,2,1,3)= r5(14)-r4(9,2)+r3(5,7)-r2(4,8)+(+r4(9,3)*two-r3(5,5)*two)*qz+(+r3(5,8) &
&                  -r2(4,9))*zz
      eri(4,2,1,3)= r5(8)-r4(4,2)+r3(3,7)-r2(1,8)+(+r4(8,3)-r3(4,5)+r2(5,11)-r1(1,10))*qx
      eri(5,2,1,3)= r5(9)-r4(5,2)+(+r4(9,3)-r3(5,5))*qx+(+r4(5,3)-r3(2,5))*qz+(+r3(5,8) &
&                  -r2(4,9))*xz
      eri(6,2,1,3)= r5(13)-r4(8,2)+r3(6,7)-r2(5,8)+(+r4(8,3)-r3(4,5)+r2(5,11)-r1(1,10))*qz
      eri(1,3,1,3)= r5(6)-r4(3,2)+r3(6,7)-r2(5,8)+(+r4(6,3)*two-r3(3,5)*two)*qx+(+r4(3,4) &
&                  -r3(1,6)+r2(5,12)-r1(1,11))*qz+(+r3(6,8)-r2(5,9))*xx+(+r3(3,9)*two &
&                  -r2(1,10)*two)*xz+(+r2(5,13)-r1(1,12))*xxz
      eri(2,3,1,3)= r5(13)-r4(8,2)+r3(6,7)-r2(5,8)+(+r4(8,4)-r3(4,6)+r2(5,12)-r1(1,11))*qz
      eri(3,3,1,3)= r5(15)-r4(10,2)+r3(6,7)*three-r2(5,8)*three+(+r4(10,3)*two+r4(10,4) &
&                  -r3(6,5)*two-r3(6,6)+r2(5,11)*two+r2(5,12)-r1(1,10)*two-r1(1,11))*qz+( &
&                  +r3(6,8)+r3(6,9)*two-r2(5,9)-r2(5,10)*two)*zz+(+r2(5,13)-r1(1,12))*zzz
      eri(4,3,1,3)= r5(9)-r4(5,2)+(+r4(9,3)-r3(5,5))*qx+(+r4(5,4)-r3(2,6))*qz+(+r3(5,9) &
&                  -r2(4,10))*xz
      eri(5,3,1,3)= r5(10)-r4(6,2)+r3(3,7)-r2(1,8)+(+r4(10,3)-r3(6,5)+r2(5,11)-r1(1,10)) &
&                  *qx+(+r4(6,3)+r4(6,4)-r3(3,5)-r3(3,6))*qz+(+r3(6,8)+r3(6,9)-r2(5,9) &
&                  -r2(5,10))*xz+(+r3(3,9)-r2(1,10))*zz+(+r2(5,13)-r1(1,12))*xzz
      eri(6,3,1,3)= r5(14)-r4(9,2)+r3(5,7)-r2(4,8)+(+r4(9,3)+r4(9,4)-r3(5,5)-r3(5,6))*qz+( &
&                  +r3(5,9)-r2(4,10))*zz
      eri(1,1,2,3)= r5(5)-r4(2,2)+r3(5,7)*three-r2(4,8)*three+(+r4(5,3)*two+r4(5,4) &
&                  -r3(2,5)*two-r3(2,6)+r2(6,11)*two+r2(6,12)-r1(2,10)*two-r1(2,11))*qx+( &
&                  +r3(5,8)+r3(5,9)*two-r2(4,9)-r2(4,10)*two)*xx+(+r2(6,13)-r1(2,12))*xxx
      eri(2,1,2,3)= r5(12)-r4(7,2)+r3(5,7)-r2(4,8)+(+r4(12,4)-r3(7,6)+r2(6,12)-r1(2,11)) &
&                  *qx
      eri(3,1,2,3)= r5(14)-r4(9,2)+r3(5,7)-r2(4,8)+(+r4(14,4)-r3(9,6)+r2(6,12)-r1(2,11)) &
&                  *qx+(+r4(9,3)*two-r3(5,5)*two)*qz+(+r3(9,9)*two-r2(6,10)*two)*xz+(+r3(5,8) &
&                  -r2(4,9))*zz+(+r2(6,13)-r1(2,12))*xzz
      eri(4,1,2,3)= r5(8)-r4(4,2)+r3(8,7)-r2(2,8)+(+r4(8,3)+r4(8,4)-r3(4,5)-r3(4,6))*qx+( &
&                  +r3(8,9)-r2(2,10))*xx
      eri(5,1,2,3)= r5(9)-r4(5,2)+r3(9,7)-r2(6,8)+(+r4(9,3)+r4(9,4)-r3(5,5)-r3(5,6))*qx+( &
&                  +r4(5,3)-r3(2,5)+r2(6,11)-r1(2,10))*qz+(+r3(9,9)-r2(6,10))*xx+(+r3(5,8) &
&                  +r3(5,9)-r2(4,9)-r2(4,10))*xz+(+r2(6,13)-r1(2,12))*xxz
      eri(6,1,2,3)= r5(13)-r4(8,2)+(+r4(13,4)-r3(8,6))*qx+(+r4(8,3)-r3(4,5))*qz+(+r3(8,9) &
&                  -r2(2,10))*xz
      eri(1,2,2,3)= r5(8)-r4(4,2)+r3(8,7)-r2(2,8)+(+r4(8,3)*two-r3(4,5)*two)*qx+(+r3(8,8) &
&                  -r2(2,9))*xx
      eri(2,2,2,3)= r5(17)-r4(11,2)+r3(8,7)*three-r2(2,8)*three
      eri(3,2,2,3)= r5(19)-r4(13,2)+r3(8,7)-r2(2,8)+(+r4(13,3)*two-r3(8,5)*two)*qz+( &
&                  +r3(8,8)-r2(2,9))*zz
      eri(4,2,2,3)= r5(12)-r4(7,2)+r3(5,7)-r2(4,8)+(+r4(12,3)-r3(7,5)+r2(6,11)-r1(2,10)) &
&                  *qx
      eri(5,2,2,3)= r5(13)-r4(8,2)+(+r4(13,3)-r3(8,5))*qx+(+r4(8,3)-r3(4,5))*qz+(+r3(8,8) &
&                  -r2(2,9))*xz
      eri(6,2,2,3)= r5(18)-r4(12,2)+r3(9,7)-r2(6,8)+(+r4(12,3)-r3(7,5)+r2(6,11)-r1(2,10)) &
&                  *qz
      eri(1,3,2,3)= r5(9)-r4(5,2)+r3(9,7)-r2(6,8)+(+r4(9,3)*two-r3(5,5)*two)*qx+(+r4(5,4) &
&                  -r3(2,6)+r2(6,12)-r1(2,11))*qz+(+r3(9,8)-r2(6,9))*xx+(+r3(5,9)*two &
&                  -r2(4,10)*two)*xz+(+r2(6,13)-r1(2,12))*xxz
      eri(2,3,2,3)= r5(18)-r4(12,2)+r3(9,7)-r2(6,8)+(+r4(12,4)-r3(7,6)+r2(6,12)-r1(2,11)) &
&                  *qz
      eri(3,3,2,3)= r5(20)-r4(14,2)+r3(9,7)*three-r2(6,8)*three+(+r4(14,3)*two+r4(14,4) &
&                  -r3(9,5)*two-r3(9,6)+r2(6,11)*two+r2(6,12)-r1(2,10)*two-r1(2,11))*qz+( &
&                  +r3(9,8)+r3(9,9)*two-r2(6,9)-r2(6,10)*two)*zz+(+r2(6,13)-r1(2,12))*zzz
      eri(4,3,2,3)= r5(13)-r4(8,2)+(+r4(13,3)-r3(8,5))*qx+(+r4(8,4)-r3(4,6))*qz+(+r3(8,9) &
&                  -r2(2,10))*xz
      eri(5,3,2,3)= r5(14)-r4(9,2)+r3(5,7)-r2(4,8)+(+r4(14,3)-r3(9,5)+r2(6,11)-r1(2,10)) &
&                  *qx+(+r4(9,3)+r4(9,4)-r3(5,5)-r3(5,6))*qz+(+r3(9,8)+r3(9,9)-r2(6,9) &
&                  -r2(6,10))*xz+(+r3(5,9)-r2(4,10))*zz+(+r2(6,13)-r1(2,12))*xzz
      eri(6,3,2,3)= r5(19)-r4(13,2)+r3(8,7)-r2(2,8)+(+r4(13,3)+r4(13,4)-r3(8,5)-r3(8,6)) &
&                  *qz+(+r3(8,9)-r2(2,10))*zz
      eri(1,1,3,3)= r5(6)-r4(3,1)-r4(3,2)+r3(1,1)+r3(1,2)+r3(6,7)*three-r2(5,5)*three &
&                  -r2(5,8)*three+r1(1,1)*three+r1(1,4)*three+(+r4(6,3)*two+r4(6,4) &
&                  -r3(3,3)*two-r3(3,4)-r3(3,5)*two-r3(3,6)+r2(1,1)*two+r2(1,2)+r2(1,3)*two &
&                  +r2(1,4)+r2(3,11)*two+r2(3,12)-r1(3,7)*two-r1(3,8)-r1(3,10)*two-r1(3,11) &
&                  +r0(1)*two+r0(2)+r0(4)*two+r0(5))*qx+(+r3(6,8)+r3(6,9)*two-r2(5,6) &
&                  -r2(5,7)*two-r2(5,9)-r2(5,10)*two+r1(1,2)+r1(1,3)*two+r1(1,5)+r1(1,6)*two) &
&                  *xx+(+r2(3,13)-r1(3,9)-r1(3,12)+r0(3)+r0(6))*xxx
      eri(2,1,3,3)= r5(13)-r4(8,1)-r4(8,2)+r3(4,1)+r3(4,2)+r3(6,7)-r2(5,5)-r2(5,8)+r1(1,1) &
&                  +r1(1,4)+(+r4(13,4)-r3(8,4)-r3(8,6)+r2(2,2)+r2(2,4)+r2(3,12)-r1(3,8) &
&                  -r1(3,11)+r0(2)+r0(5))*qx
      eri(3,1,3,3)= r5(15)-r4(10,1)-r4(10,2)+r3(6,1)+r3(6,2)+r3(6,7)-r2(5,5)-r2(5,8) &
&                  +r1(1,1)+r1(1,4)+(+r4(15,4)-r3(10,4)-r3(10,6)+r2(3,2)+r2(3,4)+r2(3,12) &
&                  -r1(3,8)-r1(3,11)+r0(2)+r0(5))*qx+(+r4(10,3)*two-r3(6,3)*two-r3(6,5)*two &
&                  +r2(5,1)*two+r2(5,3)*two)*qz+(+r3(10,9)*two-r2(3,7)*two-r2(3,10)*two &
&                  +r1(3,3)*two+r1(3,6)*two)*xz+(+r3(6,8)-r2(5,6)-r2(5,9)+r1(1,2)+r1(1,5))*zz &
&                  +(+r2(3,13)-r1(3,9)-r1(3,12)+r0(3)+r0(6))*xzz
      eri(4,1,3,3)= r5(9)-r4(5,1)-r4(5,2)+r3(2,1)+r3(2,2)+r3(9,7)-r2(6,5)-r2(6,8)+r1(2,1) &
&                  +r1(2,4)+(+r4(9,3)+r4(9,4)-r3(5,3)-r3(5,4)-r3(5,5)-r3(5,6)+r2(4,1)+r2(4,2) &
&                  +r2(4,3)+r2(4,4))*qx+(+r3(9,9)-r2(6,7)-r2(6,10)+r1(2,3)+r1(2,6))*xx
      eri(5,1,3,3)= r5(10)-r4(6,1)-r4(6,2)+r3(3,1)+r3(3,2)+r3(10,7)-r2(3,5)-r2(3,8) &
&                  +r1(3,1)+r1(3,4)+(+r4(10,3)+r4(10,4)-r3(6,3)-r3(6,4)-r3(6,5)-r3(6,6) &
&                  +r2(5,1)+r2(5,2)+r2(5,3)+r2(5,4))*qx+(+r4(6,3)-r3(3,3)-r3(3,5)+r2(1,1) &
&                  +r2(1,3)+r2(3,11)-r1(3,7)-r1(3,10)+r0(1)+r0(4))*qz+(+r3(10,9)-r2(3,7) &
&                  -r2(3,10)+r1(3,3)+r1(3,6))*xx+(+r3(6,8)+r3(6,9)-r2(5,6)-r2(5,7)-r2(5,9) &
&                  -r2(5,10)+r1(1,2)+r1(1,3)+r1(1,5)+r1(1,6))*xz+(+r2(3,13)-r1(3,9)-r1(3,12) &
&                  +r0(3)+r0(6))*xxz
      eri(6,1,3,3)= r5(14)-r4(9,1)-r4(9,2)+r3(5,1)+r3(5,2)+(+r4(14,4)-r3(9,4)-r3(9,6) &
&                  +r2(6,2)+r2(6,4))*qx+(+r4(9,3)-r3(5,3)-r3(5,5)+r2(4,1)+r2(4,3))*qz+( &
&                  +r3(9,9)-r2(6,7)-r2(6,10)+r1(2,3)+r1(2,6))*xz
      eri(1,2,3,3)= r5(9)-r4(5,1)-r4(5,2)+r3(2,1)+r3(2,2)+r3(9,7)-r2(6,5)-r2(6,8)+r1(2,1) &
&                  +r1(2,4)+(+r4(9,3)*two-r3(5,3)*two-r3(5,5)*two+r2(4,1)*two+r2(4,3)*two)*qx &
&                  +(+r3(9,8)-r2(6,6)-r2(6,9)+r1(2,2)+r1(2,5))*xx
      eri(2,2,3,3)= r5(18)-r4(12,1)-r4(12,2)+r3(7,1)+r3(7,2)+r3(9,7)*three-r2(6,5)*three &
&                  -r2(6,8)*three+r1(2,1)*three+r1(2,4)*three
      eri(3,2,3,3)= r5(20)-r4(14,1)-r4(14,2)+r3(9,1)+r3(9,2)+r3(9,7)-r2(6,5)-r2(6,8) &
&                  +r1(2,1)+r1(2,4)+(+r4(14,3)*two-r3(9,3)*two-r3(9,5)*two+r2(6,1)*two &
&                  +r2(6,3)*two)*qz+(+r3(9,8)-r2(6,6)-r2(6,9)+r1(2,2)+r1(2,5))*zz
      eri(4,2,3,3)= r5(13)-r4(8,1)-r4(8,2)+r3(4,1)+r3(4,2)+r3(6,7)-r2(5,5)-r2(5,8)+r1(1,1) &
&                  +r1(1,4)+(+r4(13,3)-r3(8,3)-r3(8,5)+r2(2,1)+r2(2,3)+r2(3,11)-r1(3,7) &
&                  -r1(3,10)+r0(1)+r0(4))*qx
      eri(5,2,3,3)= r5(14)-r4(9,1)-r4(9,2)+r3(5,1)+r3(5,2)+(+r4(14,3)-r3(9,3)-r3(9,5) &
&                  +r2(6,1)+r2(6,3))*qx+(+r4(9,3)-r3(5,3)-r3(5,5)+r2(4,1)+r2(4,3))*qz+( &
&                  +r3(9,8)-r2(6,6)-r2(6,9)+r1(2,2)+r1(2,5))*xz
      eri(6,2,3,3)= r5(19)-r4(13,1)-r4(13,2)+r3(8,1)+r3(8,2)+r3(10,7)-r2(3,5)-r2(3,8) &
&                  +r1(3,1)+r1(3,4)+(+r4(13,3)-r3(8,3)-r3(8,5)+r2(2,1)+r2(2,3)+r2(3,11) &
&                  -r1(3,7)-r1(3,10)+r0(1)+r0(4))*qz
      eri(1,3,3,3)= r5(10)-r4(6,1)-r4(6,2)+r3(3,1)+r3(3,2)+r3(10,7)-r2(3,5)-r2(3,8) &
&                  +r1(3,1)+r1(3,4)+(+r4(10,3)*two-r3(6,3)*two-r3(6,5)*two+r2(5,1)*two &
&                  +r2(5,3)*two)*qx+(+r4(6,4)-r3(3,4)-r3(3,6)+r2(1,2)+r2(1,4)+r2(3,12) &
&                  -r1(3,8)-r1(3,11)+r0(2)+r0(5))*qz+(+r3(10,8)-r2(3,6)-r2(3,9)+r1(3,2) &
&                  +r1(3,5))*xx+(+r3(6,9)*two-r2(5,7)*two-r2(5,10)*two+r1(1,3)*two &
&                  +r1(1,6)*two)*xz+(+r2(3,13)-r1(3,9)-r1(3,12)+r0(3)+r0(6))*xxz
      eri(2,3,3,3)= r5(19)-r4(13,1)-r4(13,2)+r3(8,1)+r3(8,2)+r3(10,7)-r2(3,5)-r2(3,8) &
&                  +r1(3,1)+r1(3,4)+(+r4(13,4)-r3(8,4)-r3(8,6)+r2(2,2)+r2(2,4)+r2(3,12) &
&                  -r1(3,8)-r1(3,11)+r0(2)+r0(5))*qz
      eri(3,3,3,3)= r5(21)-r4(15,1)-r4(15,2)+r3(10,1)+r3(10,2)+r3(10,7)*three &
&                  -r2(3,5)*three-r2(3,8)*three+r1(3,1)*three+r1(3,4)*three+(+r4(15,3)*two &
&                  +r4(15,4)-r3(10,3)*two-r3(10,4)-r3(10,5)*two-r3(10,6)+r2(3,1)*two+r2(3,2) &
&                  +r2(3,3)*two+r2(3,4)+r2(3,11)*two+r2(3,12)-r1(3,7)*two-r1(3,8) &
&                  -r1(3,10)*two-r1(3,11)+r0(1)*two+r0(2)+r0(4)*two+r0(5))*qz+(+r3(10,8) &
&                  +r3(10,9)*two-r2(3,6)-r2(3,7)*two-r2(3,9)-r2(3,10)*two+r1(3,2)+r1(3,3)*two &
&                  +r1(3,5)+r1(3,6)*two)*zz+(+r2(3,13)-r1(3,9)-r1(3,12)+r0(3)+r0(6))*zzz
      eri(4,3,3,3)= r5(14)-r4(9,1)-r4(9,2)+r3(5,1)+r3(5,2)+(+r4(14,3)-r3(9,3)-r3(9,5) &
&                  +r2(6,1)+r2(6,3))*qx+(+r4(9,4)-r3(5,4)-r3(5,6)+r2(4,2)+r2(4,4))*qz+( &
&                  +r3(9,9)-r2(6,7)-r2(6,10)+r1(2,3)+r1(2,6))*xz
      eri(5,3,3,3)= r5(15)-r4(10,1)-r4(10,2)+r3(6,1)+r3(6,2)+r3(6,7)-r2(5,5)-r2(5,8) &
&                  +r1(1,1)+r1(1,4)+(+r4(15,3)-r3(10,3)-r3(10,5)+r2(3,1)+r2(3,3)+r2(3,11) &
&                  -r1(3,7)-r1(3,10)+r0(1)+r0(4))*qx+(+r4(10,3)+r4(10,4)-r3(6,3)-r3(6,4) &
&                  -r3(6,5)-r3(6,6)+r2(5,1)+r2(5,2)+r2(5,3)+r2(5,4))*qz+(+r3(10,8)+r3(10,9) &
&                  -r2(3,6)-r2(3,7)-r2(3,9)-r2(3,10)+r1(3,2)+r1(3,3)+r1(3,5)+r1(3,6))*xz+( &
&                  +r3(6,9)-r2(5,7)-r2(5,10)+r1(1,3)+r1(1,6))*zz+(+r2(3,13)-r1(3,9)-r1(3,12) &
&                  +r0(3)+r0(6))*xzz
      eri(6,3,3,3)= r5(20)-r4(14,1)-r4(14,2)+r3(9,1)+r3(9,2)+r3(9,7)-r2(6,5)-r2(6,8) &
&                  +r1(2,1)+r1(2,4)+(+r4(14,3)+r4(14,4)-r3(9,3)-r3(9,4)-r3(9,5)-r3(9,6) &
&                  +r2(6,1)+r2(6,2)+r2(6,3)+r2(6,4))*qz+(+r3(9,9)-r2(6,7)-r2(6,10)+r1(2,3) &
&                  +r1(2,6))*zz
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
      do j= 1,3
        do k= 1,3
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
      do i= 1,3
        do k= 1,3
          do l= 1,6
            work(1)= eri(l,k,1,i)
            work(2)= eri(l,k,2,i)
            work(3)= eri(l,k,3,i)
            eri(l,k,1,i)= work(1)*rot(1,1)+work(2)*rot(2,1)+work(3)*rot(3,1)
            eri(l,k,2,i)= work(1)*rot(1,2)+work(2)*rot(2,2)+work(3)*rot(3,2)
            eri(l,k,3,i)= work(1)*rot(1,3)+work(2)*rot(2,3)+work(3)*rot(3,3)
          enddo
        enddo
      enddo
      do i= 1,3
        do j= 1,3
          do l= 1,6
            work(1)= eri(l,1,j,i)
            work(2)= eri(l,2,j,i)
            work(3)= eri(l,3,j,i)
            eri(l,1,j,i)= work(1)*rot(1,1)+work(2)*rot(2,1)+work(3)*rot(3,1)
            eri(l,2,j,i)= work(1)*rot(1,2)+work(2)*rot(2,2)+work(3)*rot(3,2)
            eri(l,3,j,i)= work(1)*rot(1,3)+work(2)*rot(2,3)+work(3)*rot(3,3)
          enddo
        enddo
      enddo
!
      if(nbfijkl(4) == 6)then
        do i= 1,3
          do j= 1,3
            do k= 1,3
              do l= 1,6
                phmdint(l,k,j,i)= &
&                   eri(1,k,j,i)*rot2(1,l)+eri(2,k,j,i)*rot2(2,l)+eri(3,k,j,i)*rot2(3,l) &
&                  +eri(4,k,j,i)*rot2(4,l)+eri(5,k,j,i)*rot2(5,l)+eri(6,k,j,i)*rot2(6,l) 
              enddo
            enddo
          enddo
        enddo
      else
        do l= 1,6
          rot3(l,1)= rot2(l,2)
          rot3(l,2)= rot2(l,5)
          rot3(l,3)= rot2(l,6)-(rot2(l,1)+rot2(l,4))*half
          rot3(l,4)= rot2(l,3)
          rot3(l,5)=(rot2(l,1)-rot2(l,4))*sqrt3h
        enddo
        do i= 1,3
          do j= 1,3
            do k= 1,3
              do l= 1,5
                phmdint(l,k,j,i)= &
                    eri(1,k,j,i)*rot3(1,l)+eri(2,k,j,i)*rot3(2,l)+eri(3,k,j,i)*rot3(3,l) &
&                  +eri(4,k,j,i)*rot3(4,l)+eri(5,k,j,i)*rot3(5,l)+eri(6,k,j,i)*rot3(6,l)
              enddo
            enddo
          enddo
        enddo
      endif
!
      return
end


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
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:6)
      real(8) :: f0(2), f1(2,4), f2(3,4), f3(4,4), f4(5,4), f5(6,2), f6(7), ftw(6,4)
      real(8) :: r0(10), r1(3,13), r2(6,17), r3(10,12), r4(15,8), r5(21,3), r6(28)
      real(8) :: ex12, ex34, ex43, ex14, ex41, expq, expq2, expq4, ex3q, ex4q, c12, c34, zip
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
            do ii= 0,6
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
            fac= c12*ex41
            expq= expq*two*ex41*ex41
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
