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
!--------------------------------------------------------------
  subroutine int2ssss(phmdint,exfac1,exfac2,xyziq,work,nijkl)
!--------------------------------------------------------------
!
! Calculate (ss|ss) integral
!
! In  : exfac12, exfac34 (exponents and coefficients of primitive pair functions)
!       xyziq    (x,y,z length of i and q)
!       nij, nkl (number of primitive pair functions)
! Out : phmdint  (two-electron integral)
!       work     (work space)
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,parameter :: mxprsh2=30
      integer,intent(in) :: nijkl(2)
      integer :: ij, kl, igrid, i1
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, sqrtpi4=0.8862269254527580D+00
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*)
      real(8),intent(out) :: phmdint(6,6,6,6), work(2,*)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, f0, ft0, r0
      real(8) :: ex12, ex34, ex14, ex41, expq, c12, c34, zip, xiq, yiq, ziq, xypq2, zpq, zpq2
!
! Zero-clear
!
      r0= zero
!
      do kl= 1,nijkl(2)
        ex34= exfac2(1,kl)
        c34 = exfac2(5,kl)
        xiq = xyziq(1,kl)
        yiq = xyziq(2,kl)
        ziq = xyziq(3,kl)
        xypq2= xiq*xiq+yiq*yiq
        f0= zero
        i1= 0
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
            i1= i1+1
            work(1,i1)= tval
            work(2,i1)= c12
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
            ft0= fgrid(0,0,igrid)      +fgrid( 1,0,igrid)*tval  +fgrid( 2,0,igrid)*tval2 &
&               +fgrid(3,0,igrid)*tval3+fgrid( 4,0,igrid)*tval4 +fgrid( 5,0,igrid)*tval5 &
&               +fgrid(6,0,igrid)*tval6+fgrid( 7,0,igrid)*tval7 +fgrid( 8,0,igrid)*tval8 &
&               +fgrid(9,0,igrid)*tval9+fgrid(10,0,igrid)*tval10
            ft0= ft0*c12*ex41
            f0= f0+ft0
          endif
        enddo
        do ij= 1,i1
          tinv= work(2,ij)/sqrt(work(1,ij))
          f0= f0+sqrtpi4*tinv
        enddo
        r0= r0+f0*c34
      enddo
!
      phmdint(1,1,1,1)= r0
      return
end


!-----------------------------------------------------------------------
  subroutine int2psss(phmdint,exfac1,exfac2,xyziq,xzkl,rot,work,nijkl)
!-----------------------------------------------------------------------
!
! Calculate (ps|ss) integrals
!
! In  : exfac12, exfac34 (exponents and coefficients of primitive pair functions)
!       xyziq    (x,y,z elements of i and q)
!       xzkl     (x,z elements of k and l)
!       nij, nkl (number of primitive pair functions)
! Out : phmdint  (two-electron integral)
!       work     (work space)
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,parameter :: mxprsh2=30
      integer,intent(in) :: nijkl(2)
      integer :: ij, kl, igrid, i1
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, sqrtpi4=0.8862269254527580D+00
      real(8),parameter :: two=2.0D+00, half=0.5D+00
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*), xzkl(2), rot(3,3)
      real(8),intent(out) :: phmdint(6,6,6,6), work(4,*)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:1)
      real(8) :: f0, f1(2), r0, r1(3)
      real(8) :: ex12, ex34, ex43, ex41,ex14, expq, ex3q, c12, c34, zip, xiq, yiq, ziq
      real(8) :: xypq2, zpq, zpq2, fac, exq, qx, qz, eri(3)
!
! Zero-clear
!
      r0= zero
      r1(1)= zero
      r1(2)= zero
      r1(3)= zero
!
      do kl= 1,nijkl(2)
        ex34= exfac2(1,kl)
        ex43= exfac2(2,kl)
        ex3q= exfac2(3,kl)
        c34 = exfac2(5,kl)
        xiq = xyziq(1,kl)
        yiq = xyziq(2,kl)
        ziq = xyziq(3,kl)
        xypq2= xiq*xiq+yiq*yiq
        f0= zero
        f1(1)= zero
        f1(2)= zero
        i1=0
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
            i1= i1+1
            work(1,i1)= tval
            work(2,i1)= c12
            work(3,i1)= expq
            work(4,i1)= zpq
          else
            ex41=one/sqrt(ex14)
            tval=tval*ex41*ex41
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
            ft(0)= fgrid(0,0,igrid)      +fgrid( 1,0,igrid)*tval  +fgrid( 2,0,igrid)*tval2 &
&                 +fgrid(3,0,igrid)*tval3+fgrid( 4,0,igrid)*tval4 +fgrid( 5,0,igrid)*tval5 &
&                 +fgrid(6,0,igrid)*tval6+fgrid( 7,0,igrid)*tval7 +fgrid( 8,0,igrid)*tval8 &
&                 +fgrid(9,0,igrid)*tval9+fgrid(10,0,igrid)*tval10
            ft(1)= fgrid(0,1,igrid)      +fgrid( 1,1,igrid)*tval  +fgrid( 2,1,igrid)*tval2 &
&                 +fgrid(3,1,igrid)*tval3+fgrid( 4,1,igrid)*tval4 +fgrid( 5,1,igrid)*tval5 &
&                 +fgrid(6,1,igrid)*tval6+fgrid( 7,1,igrid)*tval7 +fgrid( 8,1,igrid)*tval8 &
&                 +fgrid(9,1,igrid)*tval9+fgrid(10,1,igrid)*tval10
            fac= c12*ex41
            ft(1)= ft(1)*fac*expq*ex41*ex41*two
            f0= f0+ft(0)*fac
            f1(1)= f1(1)-ft(1)
            f1(2)= f1(2)-ft(1)*zpq
          endif
        enddo
        do ij= 1,i1
          tinv= one/sqrt(work(1,ij))
          ft(0)= work(2,ij)*sqrtpi4*tinv
          ft(1)= ft(0)*tinv*tinv*work(3,ij)
          f0= f0+ft(0)
          f1(1)= f1(1)-ft(1)
          f1(2)= f1(2)-ft(1)*work(4,ij)
        enddo
        exq= ex43*c34
        r0= r0+f0*ex3q*c34
        r1(1)= r1(1)+f1(1)*exq*xiq
        r1(2)= r1(2)+f1(1)*exq*yiq
        r1(3)= r1(3)+f1(2)*exq
      enddo
!
      qx= xzkl(1)
      qz= xzkl(2)
      eri(1)= r1(1)+r0*qx
      eri(2)= r1(2)
      eri(3)= r1(3)+r0*qz
!
! Rotate back to the original axis
!
      phmdint(1,1,1,1)= eri(1)*rot(1,1)+eri(2)*rot(2,1)+eri(3)*rot(3,1)
      phmdint(2,1,1,1)= eri(1)*rot(1,2)+eri(2)*rot(2,2)+eri(3)*rot(3,2)
      phmdint(3,1,1,1)= eri(1)*rot(1,3)+eri(2)*rot(2,3)+eri(3)*rot(3,3)
      return
end


!------------------------------------------------------------------------
  subroutine int2ppss(phmdint,exfac1,exfac2,xyziq,xzkl,rot,work2,nijkl)
!------------------------------------------------------------------------
!
! Calculate (pp|ss) integrals
!
! In  : exfac12, exfac34 (exponents and coefficients of primitive pair functions)
!       xyziq    (x,y,z elements of i and q)
!       xzkl     (x,z elements of k and l)
!       nij, nkl (number of primitive pair functions)
! Out : phmdint  (two-electron integral)
!       work2    (work space)
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,parameter :: mxprsh2=30
      integer,intent(in) :: nijkl(2)
      integer :: ij, kl, igrid, ii, i, k, l, i1
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, sqrtpi4=0.8862269254527580D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*), xzkl(2), rot(3,3)
      real(8),intent(out) :: phmdint(6,6,6,6), work2(4,*)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:2)
      real(8) :: f0, f1(2), f2(3), r0(2), r1(3,2), r2(6)
      real(8) :: ex12, ex34, ex43, ex14, ex41, expq, ex3q, ex4q, c12, c34, zip, xiq, yiq, ziq
      real(8) :: xiq2, yiq2, xypq2, zpq, zpq2, fac, qmd, qmd2, qx, qz, eri(3,3), work(3)
!
! Zero-clear
!
      r0(1:2)= zero
      r1(1:3,1:2)= zero
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
        i1= 0
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
            i1= i1+1
            work2(1,i1)= tval
            work2(2,i1)= c12
            work2(3,i1)= expq
            work2(4,i1)= zpq
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
            f0= f0+ft(0)
            f1(1)= f1(1)-ft(1)
            f1(2)= f1(2)-ft(1)*zpq
            f2(1)= f2(1)+ft(2)
            f2(2)= f2(2)+ft(2)*zpq
            f2(3)= f2(3)+ft(2)*zpq2
          endif
        enddo
        do ij= 1,i1
          tinv= one/sqrt(work2(1,ij))
          ft(0)= work2(2,ij)*sqrtpi4*tinv
          expq= work2(3,ij)*tinv*tinv
          ft(1)= ft(0)*expq
          ft(2)= ft(1)*expq*three
          f0= f0+ft(0)
          f1(1)= f1(1)-ft(1)
          f1(2)= f1(2)-ft(1)*work2(4,ij)
          f2(1)= f2(1)+ft(2)
          f2(2)= f2(2)+ft(2)*work2(4,ij)
          f2(3)= f2(3)+ft(2)*work2(4,ij)*work2(4,ij)
        enddo
        qmd = ex43*c34
        qmd2= qmd*ex43
        r0(1)= r0(1)+f0*qmd
        r0(2)= r0(2)+f0*ex3q*ex4q*c34
        work(1)= ex3q*qmd
        work(2)= ex4q*qmd
        do i= 1,2
          r1(1,i)= r1(1,i)+f1(1)*work(i)*xiq
          r1(2,i)= r1(2,i)+f1(1)*work(i)*yiq
          r1(3,i)= r1(3,i)+f1(2)*work(i)
        enddo
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
      eri(1,1)= r2(1)+r1(1,1)*qx+r1(1,2)*qx+r0(1)+r0(2)*qx*qx
      eri(2,1)= r2(4)           +r1(2,2)*qx
      eri(3,1)= r2(5)+r1(1,1)*qz+r1(3,2)*qx      +r0(2)*qx*qz
      eri(1,2)= r2(4)+r1(2,1)*qx
      eri(2,2)= r2(2)                      +r0(1)
      eri(3,2)= r2(6)+r1(2,1)*qz
      eri(1,3)= r2(5)+r1(3,1)*qx+r1(1,2)*qz      +r0(2)*qx*qz
      eri(2,3)= r2(6)           +r1(2,2)*qz
      eri(3,3)= r2(3)+r1(3,1)*qz+r1(3,2)*qz+r0(1)+r0(2)*qz*qz
!
      do l= 1,3
        work(1)= eri(l,1)
        work(2)= eri(l,2)
        work(3)= eri(l,3)
        eri(l,1)= work(1)*rot(1,1)+work(2)*rot(2,1)+work(3)*rot(3,1)
        eri(l,2)= work(1)*rot(1,2)+work(2)*rot(2,2)+work(3)*rot(3,2)
        eri(l,3)= work(1)*rot(1,3)+work(2)*rot(2,3)+work(3)*rot(3,3)
      enddo
      do k= 1,3
        work(1)= eri(1,k)
        work(2)= eri(2,k)
        work(3)= eri(3,k)
        eri(1,k)= work(1)*rot(1,1)+work(2)*rot(2,1)+work(3)*rot(3,1)
        eri(2,k)= work(1)*rot(1,2)+work(2)*rot(2,2)+work(3)*rot(3,2)
        eri(3,k)= work(1)*rot(1,3)+work(2)*rot(2,3)+work(3)*rot(3,3)
      enddo
!
      do k= 1,3
        do l= 1,3
          phmdint(l,k,1,1)= eri(l,k)
        enddo
      enddo
      return
end


!-----------------------------------------------------------------------
  subroutine int2psps(phmdint,exfac1,exfac2,xyziq,xzkl,rot,work,nijkl)
!-----------------------------------------------------------------------
!
! Calculate (ps|ps) integrals
!
! In  : exfac12, exfac34 (exponents and coefficients of primitive pair functions)
!       xyziq    (x,y,z elements of i and q)
!       xzkl     (x,z elements of k and l)
!       nij, nkl (number of primitive pair functions)
! Out : phmdint  (two-electron integral)
!       work     (work space)
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,parameter :: mxprsh2=30
      integer,intent(in) :: nijkl(2)
      integer :: ij, kl, igrid, ii, j, l, i1
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, sqrtpi4=0.8862269254527580D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*), xzkl(2), rot(3,3)
      real(8),intent(out) :: phmdint(6,6,6,6), work(6,*)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:2)
      real(8) :: f0, f1(4), f2(3), r0, r1(3,2), r2(6)
      real(8) :: ex12, ex34, ex43, ex14, ex41, expq, ex3q, c12, c34, zip, zjp, xiq, yiq, ziq
      real(8) :: xiq2, yiq2, xypq2, zpq, zpq2, fac, pmd, qmd, qx, qz, eri(3,3), eri2(3,3)
!
! Zero-clear
!
      r0= zero
      r1(1:3,1:2)= zero
      r2(1:6)= zero
!
      do kl= 1,nijkl(2)
        ex34= exfac2(1,kl)
        ex43= exfac2(2,kl)
        ex3q= exfac2(3,kl)
        c34 = exfac2(5,kl)
        xiq = xyziq(1,kl)
        yiq = xyziq(2,kl)
        ziq = xyziq(3,kl)
        xiq2= xiq*xiq
        yiq2= yiq*yiq
        xypq2= xiq2+yiq2
        f0= zero
        f1(1:4)= zero
        f2(1:3)= zero
        i1= 0
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
            i1= i1+1
            work(1,i1)= tval
            work(2,i1)= c12
            work(3,i1)= expq
            work(4,i1)= zjp
            work(5,i1)= zpq
            work(6,i1)= pmd
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
            f0   = f0   +ft(0)*zjp
            f1(1)= f1(1)-ft(1)*zjp
            f1(2)= f1(2)-ft(1)*zjp*zpq
            f1(3)= f1(3)-ft(1)*pmd
            f1(4)= f1(4)-ft(1)*pmd*zpq
            f2(1)= f2(1)+ft(2)*pmd
            f2(2)= f2(2)+ft(2)*pmd*zpq
            f2(3)= f2(3)+ft(2)*pmd*zpq2
          endif
        enddo
        do ij= 1,i1
          tinv= one/sqrt(work(1,ij))
          ft(0)= work(2,ij)*sqrtpi4*tinv
          expq= work(3,ij)*tinv*tinv
          ft(1)= ft(0)*expq
          ft(2)= ft(1)*expq*three
          f0   = f0   +ft(0)*work(4,ij)
          f1(1)= f1(1)-ft(1)*work(4,ij)
          f1(2)= f1(2)-ft(1)*work(4,ij)*work(5,ij)
          f1(3)= f1(3)-ft(1)*work(6,ij)
          f1(4)= f1(4)-ft(1)*work(6,ij)*work(5,ij)
          f2(1)= f2(1)+ft(2)*work(6,ij)
          f2(2)= f2(2)+ft(2)*work(6,ij)*work(5,ij)
          f2(3)= f2(3)+ft(2)*work(6,ij)*work(5,ij)*work(5,ij)
        enddo
        qmd = ex43*c34
        ex3q= ex3q*c34
        r0     = r0+f0*ex3q
        r1(1,1)= r1(1,1)+f1(1)*qmd*xiq
        r1(2,1)= r1(2,1)+f1(1)*qmd*yiq
        r1(3,1)= r1(3,1)+f1(2)*qmd
        r1(1,2)= r1(1,2)+f1(3)*ex3q*xiq
        r1(2,2)= r1(2,2)+f1(3)*ex3q*yiq
        r1(3,2)= r1(3,2)+f1(4)*ex3q
        r2(1)  = r2(1)+(f2(1)*xiq2+f1(3))*qmd
        r2(2)  = r2(2)+(f2(1)*yiq2+f1(3))*qmd
        r2(3)  = r2(3)+(f2(3)     +f1(3))*qmd
        r2(4)  = r2(4)+(f2(1)*xiq*yiq   )*qmd
        r2(5)  = r2(5)+(f2(2)*xiq       )*qmd
        r2(6)  = r2(6)+(f2(2)*yiq       )*qmd
      enddo
!
      qx= xzkl(1)
      qz= xzkl(2)
      eri(1,1)=-r2(1)        -r1(1,2)*qx
      eri(2,1)=-r2(4)
      eri(3,1)=-r2(5)        -r1(1,2)*qz
      eri(1,2)=-r2(4)        -r1(2,2)*qx
      eri(2,2)=-r2(2)
      eri(3,2)=-r2(6)        -r1(2,2)*qz
      eri(1,3)=-r2(5)+r1(1,1)-r1(3,2)*qx+r0*qx
      eri(2,3)=-r2(6)+r1(2,1)
      eri(3,3)=-r2(3)+r1(3,1)-r1(3,2)*qz+r0*qz
!
      do j=1,3
        do l= 1,3
          eri2(l,j)= eri(l,1)*rot(1,j)+eri(l,2)*rot(2,j)+eri(l,3)*rot(3,j)
        enddo
      enddo
      do j= 1,3
        do l=1,3
         eri(l,j)= eri2(1,j)*rot(1,l)+eri2(2,j)*rot(2,l)+eri2(3,j)*rot(3,l)
        enddo
      enddo
!
      do j= 1,3
        do l= 1,3
          phmdint(l,1,j,1)= eri(l,j)
        enddo
      enddo
      return
end


!------------------------------------------------------------------
  subroutine int2ppps(phmdint,exfac1,exfac2,xyziq,xzkl,rot,nijkl)
!------------------------------------------------------------------
!
! Calculate (pp|ps) integrals
!
! In  : exfac12, exfac34 (exponents and coefficients of primitive pair functions)
!       xyziq    (x,y,z elements of i and q)
!       xzkl     (x,z elements of k and l)
!       nij, nkl (number of primitive pair functions)
! Out : phmdint  (two-electron integral)
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nijkl(2)
      integer :: ij, kl, igrid, ii, i, j, k, l
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, sqrtpi4=0.8862269254527580D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00, five=5.0D+00
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*), xzkl(2), rot(3,3)
      real(8),intent(out) :: phmdint(6,6,6,6)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:3)
      real(8) :: f0, f1(2,2), f2(3,2), f3(4), r0(2), r1(3,4), r2(6,3), r3(10)
      real(8) :: ex12, ex34, ex43, ex14, ex41, expq, expq2, ex3q, ex4q, c12, c34, zip, zjp
      real(8) :: xiq, yiq, ziq, xiq2, yiq2, xyiq, xypq2, zpq, zpq2, fac, pmd
      real(8) :: qmd, qmd2, qmd2x, qmd2y, qx, qz, xx, xz, zz, eri(3,3,3), work(4), fw(6)
      real(8) :: eri2(3,3,3)
!
! Zero-clear
!
      r0(1:2)= zero
      r1(1:3,1:4)= zero
      r2(1:6,1:3)= zero
      r3(1:10)= zero
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
            expq = expq*two*ex41*ex41
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
        qmd  = ex43*c34
        qmd2 = ex43*qmd
        qmd2x= qmd2*xiq
        qmd2y= qmd2*yiq
        work(1)= ex3q*ex4q*c34
        work(2)= qmd
        work(3)= qmd*ex3q
        work(4)= qmd*ex4q
!
        r0(1)= r0(1)+f0*work(1)
        r0(2)= r0(2)+f0*work(2)
!
        do i= 1,2
           r1(1,i)= r1(1,i)+f1(1,1)*work(i+2)*xiq
           r1(2,i)= r1(2,i)+f1(1,1)*work(i+2)*yiq
           r1(3,i)= r1(3,i)+f1(2,1)*work(i+2)
        enddo
        do i= 3,4
           r1(1,i)= r1(1,i)+f1(1,2)*work(i-2)*xiq
           r1(2,i)= r1(2,i)+f1(1,2)*work(i-2)*yiq
           r1(3,i)= r1(3,i)+f1(2,2)*work(i-2)
        enddo
!
        r2(1,1)= r2(1,1)+(f2(1,1)*xiq2+f1(1,1))*qmd2
        r2(2,1)= r2(2,1)+(f2(1,1)*yiq2+f1(1,1))*qmd2
        r2(3,1)= r2(3,1)+(f2(3,1)     +f1(1,1))*qmd2
        r2(4,1)= r2(4,1)+(f2(1,1)*xyiq        )*qmd2
        r2(5,1)= r2(5,1)+(f2(2,1)*xiq         )*qmd2
        r2(6,1)= r2(6,1)+(f2(2,1)*yiq         )*qmd2
        fw(1)= f2(1,2)*xiq2+f1(1,2)
        fw(2)= f2(1,2)*yiq2+f1(1,2)
        fw(3)= f2(3,2)     +f1(1,2)
        fw(4)= f2(1,2)*xyiq
        fw(5)= f2(2,2)*xiq
        fw(6)= f2(2,2)*yiq
        do i= 2,3
          do j= 1,6
            r2(j,i)= r2(j,i)+fw(j)*work(i+1)
          enddo
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
      eri(1,1,1)=-r3(1)-r2(1,2)*qx-r2(1,3)*qx-r1(1,3)*xx-r1(1,4)
      eri(2,1,1)=-r3(2)           -r2(4,3)*qx
      eri(3,1,1)=-r3(3)-r2(1,2)*qz-r2(5,3)*qx-r1(1,3)*xz
      eri(1,2,1)=-r3(2)-r2(4,2)*qx
      eri(2,2,1)=-r3(4)                                 -r1(1,4)
      eri(3,2,1)=-r3(5)-r2(4,2)*qz
      eri(1,3,1)=-r3(3)-r2(5,2)*qx-r2(1,3)*qz-r1(1,3)*xz
      eri(2,3,1)=-r3(5)           -r2(4,3)*qz
      eri(3,3,1)=-r3(6)-r2(5,2)*qz-r2(5,3)*qz-r1(1,3)*zz-r1(1,4)
      eri(1,1,2)=-r3(2)-r2(4,2)*qx-r2(4,3)*qx-r1(2,3)*xx-r1(2,4)
      eri(2,1,2)=-r3(4)           -r2(2,3)*qx
      eri(3,1,2)=-r3(5)-r2(4,2)*qz-r2(6,3)*qx-r1(2,3)*xz
      eri(1,2,2)=-r3(4)-r2(2,2)*qx
      eri(2,2,2)=-r3(7)                                 -r1(2,4)
      eri(3,2,2)=-r3(8)-r2(2,2)*qz
      eri(1,3,2)=-r3(5)-r2(6,2)*qx-r2(4,3)*qz-r1(2,3)*xz
      eri(2,3,2)=-r3(8)           -r2(2,3)*qz
      eri(3,3,2)=-r3(9)-r2(6,2)*qz-r2(6,3)*qz-r1(2,3)*zz-r1(2,4)
      eri(1,1,3)=-r3(3)+r2(1,1)-r2(5,2)*qx-r2(5,3)*qx+r1(1,1)*qx+r1(1,2)*qx &
&                -r1(3,3)*xx-r1(3,4)+r0(1)*xx+r0(2)
      eri(2,1,3)=-r3(5)+r2(4,1)-r2(6,3)*qx+r1(2,2)*qx
      eri(3,1,3)=-r3(6)+r2(5,1)-r2(5,2)*qz-r2(3,3)*qx+r1(1,1)*qz+r1(3,2)*qx &
&                -r1(3,3)*xz+r0(1)*xz
      eri(1,2,3)=-r3(5)+r2(4,1)-r2(6,2)*qx+r1(2,1)*qx
      eri(2,2,3)=-r3(8)+r2(2,1)-r1(3,4)+r0(2)
      eri(3,2,3)=-r3(9)+r2(6,1)-r2(6,2)*qz+r1(2,1)*qz
      eri(1,3,3)=-r3(6)+r2(5,1)-r2(3,2)*qx-r2(5,3)*qz+r1(3,1)*qx+r1(1,2)*qz &
&                -r1(3,3)*xz+r0(1)*xz
      eri(2,3,3)=-r3(9)+r2(6,1)-r2(6,3)*qz+r1(2,2)*qz
      eri(3,3,3)=-r3(10)+r2(3,1)-r2(3,2)*qz-r2(3,3)*qz+r1(3,1)*qz+r1(3,2)*qz &
&                -r1(3,3)*zz-r1(3,4)+r0(1)*zz+r0(2)
!
      do j=1,3
        do k= 1,3
          do l= 1,3
            eri2(l,k,j)= eri(l,k,1)*rot(1,j)+eri(l,k,2)*rot(2,j)+eri(l,k,3)*rot(3,j)
          enddo
        enddo
      enddo
      do j= 1,3
        do k=1,3
          do l= 1,3
            eri(l,k,j)= eri2(l,1,j)*rot(1,k)+eri2(l,2,j)*rot(2,k)+eri2(l,3,j)*rot(3,k)
          enddo
        enddo
      enddo
      do j= 1,3
        do k= 1,3
          do l= 1,3
            eri2(l,k,j)= eri(1,k,j)*rot(1,l)+eri(2,k,j)*rot(2,l)+eri(3,k,j)*rot(3,l)
          enddo
        enddo
      enddo
!
      do j= 1,3
        do k= 1,3
          do l= 1,3
            phmdint(l,k,j,1)= eri2(l,k,j)
          enddo
        enddo
      enddo
!
      return
end


!------------------------------------------------------------------
  subroutine int2pppp(phmdint,exfac1,exfac2,xyziq,xzkl,rot,nijkl)
!------------------------------------------------------------------
!
! Calculate (pp|pp) integrals
!
! In  : exfac12, exfac34 (exponents and coefficients of primitive pair functions)
!       xyziq    (x,y,z elements of i and q)
!       xzkl     (x,z elements of k and l)
!       nij, nkl (number of primitive pair functions)
! Out : phmdint  (two-electron integral)
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nijkl(2)
      integer :: ij, kl, igrid, ii, i, j, k, l
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, sqrtpi4=0.8862269254527580D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00, five=5.0D+00
      real(8),parameter :: six=6.0D+00, seven=7.0D+00
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*), xzkl(2), rot(3,3)
      real(8),intent(out) :: phmdint(6,6,6,6)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:4)
      real(8) :: f0(2), f1(2,5), f2(3,5), f3(4,3), f4(5), r0(4), r1(3,8), r2(6,8)
      real(8) :: r3(10,4), r4(15), f1w(3,4), f2w(6,5), f3w(10,3), ftw(4,5)
      real(8) :: ex12, ex34, ex43, ex14, ex41, expq, expq2, ex3q, ex4q, c12, c34, zip, zjp
      real(8) :: xiq, yiq, ziq, xiq2, yiq2, xiq4, yiq4, xyiq, xyiq2, xypq2, zpq, zpq2, zpq3
      real(8) :: fac, pmd, qmd, qmd2, qmd2x, qmd2y, qmd2xy, qx, qz, xx, xz, zz, eri(3,3,3,3)
      real(8) :: work(5)
!
! Zero-clear
!
      r0(1:4)= zero
      r1(1:3,1:8)= zero
      r2(1:6,1:8)= zero
      r3(1:10,1:4)= zero
      r4(1:15)= zero
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
        xiq4= xiq2*xiq2
        yiq4= yiq2*yiq2
        xyiq= xiq*yiq
        xyiq2= xyiq*xyiq
        xypq2= xiq2+yiq2
        f0(1:2)= zero
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
          zpq3= zpq*zpq2
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
            tval3= tval *tval *tval
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
            expq = expq*two*ex41*ex41
            expq2= expq*expq
            ft(0)= ft(0)*fac
            ft(1)= ft(1)*fac*expq
            ft(2)= ft(2)*fac*expq2
            ft(3)= ft(3)*fac*expq2*expq
            ft(4)= ft(4)*fac*expq2*expq2
          endif
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
        qmd   = ex43*c34
        qmd2  = ex43*qmd
        qmd2x = qmd2*xiq
        qmd2y = qmd2*yiq
        qmd2xy= qmd2*xyiq
        work(1)= ex3q*ex4q*c34
        work(2)= qmd
        work(3)= qmd*ex3q
        work(4)= qmd*ex4q
        work(5)= qmd2
!
        r0(1)= r0(1)+f0(1)*work(1)
        r0(2)= r0(2)+f0(1)*work(2)
        r0(3)= r0(3)+f0(2)*work(1)
        r0(4)= r0(4)+f0(2)*work(2)
!
        do i= 1,4
          f1w(1,i)= f1(1,i)*xiq
          f1w(2,i)= f1(1,i)*yiq
          f1w(3,i)= f1(2,i)
        enddo
        do i= 1,2
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,1)*work(i+2)
          enddo
        enddo
        do i= 3,4
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,2)*work(i)
          enddo
        enddo
        do i= 5,6
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,3)*work(i-4)
          enddo
        enddo
        do i= 7,8
          do j= 1,3
            r1(j,i)= r1(j,i)+f1w(j,4)*work(i-6)
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
            r2(j,i)= r2(j,i)+f2w(j,i)*work(5)
          enddo
        enddo
        do i= 3,4
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,3)*work(i)
          enddo
        enddo
        do i= 5,6
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,4)*work(i-2)
          enddo
        enddo
        do i= 7,8
          do j= 1,6
            r2(j,i)= r2(j,i)+f2w(j,5)*work(i-6)
          enddo
        enddo
!
        do i= 1,3
          f3w( 1,i)=(f3(1,i)*xiq2+f2(1,i+2)*three)*xiq
          f3w( 2,i)=(f3(1,i)*xiq2+f2(1,i+2)      )*yiq
          f3w( 3,i)= f3(2,i)*xiq2+f2(2,i+2)
          f3w( 4,i)=(f3(1,i)*yiq2+f2(1,i+2)      )*xiq
          f3w( 5,i)= f3(2,i)*xyiq
          f3w( 6,i)=(f3(3,i)     +f2(1,i+2)      )*xiq
          f3w( 7,i)=(f3(1,i)*yiq2+f2(1,i+2)*three)*yiq
          f3w( 8,i)= f3(2,i)*yiq2+f2(2,i+2)
          f3w( 9,i)=(f3(3,i)     +f2(1,i+2)      )*yiq
          f3w(10,i)= f3(4,i)     +f2(2,i+2)*three
        enddo
        do j= 1,10
          r3(j,1)= r3(j,1)+f3w(j,1)*work(5)
        enddo
        do j= 1,10
          r3(j,2)= r3(j,2)+f3w(j,2)*work(5)
        enddo
        do i= 3,4
          do j= 1,10
            r3(j,i)= r3(j,i)+f3w(j,3)*work(i)
          enddo
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
      eri(1,1,1,1)= r4( 1)+r3(1,3)*qx+r3(1,4)*qx+r2(1,2)+r2(1,7)*xx+r2(1,8)+r1(1,3)*qx &
&                  +r1(1,4)*qx+r0(3)*xx+r0(4)
      eri(2,1,1,1)= r4( 2)+r3(2,4)*qx+r2(4,2)+r1(2,4)*qx
      eri(3,1,1,1)= r4( 3)+r3(1,3)*qz+r3(3,4)*qx+r2(5,2)+r2(1,7)*xz+r1(1,3)*qz+r1(3,4)*qx+r0(3)*xz
      eri(1,2,1,1)= r4( 2)+r3(2,3)*qx+r2(4,2)+r1(2,3)*qx
      eri(2,2,1,1)= r4( 4)+r2(2,2)+r2(1,8)+r0(4)
      eri(3,2,1,1)= r4( 5)+r3(2,3)*qz+r2(6,2)+r1(2,3)*qz
      eri(1,3,1,1)= r4( 3)+r3(3,3)*qx+r3(1,4)*qz+r2(5,2)+r2(1,7)*xz+r1(3,3)*qx+r1(1,4)*qz+r0(3)*xz
      eri(2,3,1,1)= r4( 5)+r3(2,4)*qz+r2(6,2)+r1(2,4)*qz
      eri(3,3,1,1)= r4( 6)+r3(3,3)*qz+r3(3,4)*qz+r2(3,2)+r2(1,7)*zz+r2(1,8)+r1(3,3)*qz &
&                  +r1(3,4)*qz+r0(3)*zz+r0(4)
      eri(1,1,2,1)= r4( 2)+r3(2,3)*qx+r3(2,4)*qx+r2(4,7)*xx+r2(4,8)
      eri(2,1,2,1)= r4( 4)+r3(4,4)*qx
      eri(3,1,2,1)= r4( 5)+r3(2,3)*qz+r3(5,4)*qx+r2(4,7)*xz
      eri(1,2,2,1)= r4( 4)+r3(4,3)*qx
      eri(2,2,2,1)= r4( 7)+r2(4,8)
      eri(3,2,2,1)= r4( 8)+r3(4,3)*qz
      eri(1,3,2,1)= r4( 5)+r3(5,3)*qx+r3(2,4)*qz+r2(4,7)*xz
      eri(2,3,2,1)= r4( 8)+r3(4,4)*qz
      eri(3,3,2,1)= r4( 9)+r3(5,3)*qz+r3(5,4)*qz+r2(4,7)*zz+r2(4,8)
      eri(1,1,3,1)= r4( 3)-r3(1,1)+r3(3,3)*qx+r3(3,4)*qx-r2(1,3)*qx-r2(1,4)*qx+r2(5,7)*xx &
&                  +r2(5,8)-r1(1,5)*xx-r1(1,6)
      eri(2,1,3,1)= r4( 5)-r3(2,1)+r3(5,4)*qx-r2(4,4)*qx
      eri(3,1,3,1)= r4( 6)-r3(3,1)+r3(3,3)*qz+r3(6,4)*qx-r2(1,3)*qz-r2(5,4)*qx+r2(5,7)*xz &
&                  -r1(1,5)*xz
      eri(1,2,3,1)= r4( 5)-r3(2,1)+r3(5,3)*qx-r2(4,3)*qx
      eri(2,2,3,1)= r4( 8)-r3(4,1)+r2(5,8)-r1(1,6)
      eri(3,2,3,1)= r4( 9)-r3(5,1)+r3(5,3)*qz-r2(4,3)*qz
      eri(1,3,3,1)= r4( 6)-r3(3,1)+r3(6,3)*qx+r3(3,4)*qz-r2(5,3)*qx-r2(1,4)*qz+r2(5,7)*xz &
&                  -r1(1,5)*xz
      eri(2,3,3,1)= r4( 9)-r3(5,1)+r3(5,4)*qz-r2(4,4)*qz
      eri(3,3,3,1)= r4(10)-r3(6,1)+r3(6,3)*qz+r3(6,4)*qz-r2(5,3)*qz-r2(5,4)*qz+r2(5,7)*zz &
&                  +r2(5,8)-r1(1,5)*zz-r1(1,6)
      eri(1,1,1,2)= eri(1,1,2,1)
      eri(2,1,1,2)= eri(2,1,2,1)
      eri(3,1,1,2)= eri(3,1,2,1)
      eri(1,2,1,2)= eri(1,2,2,1)
      eri(2,2,1,2)= eri(2,2,2,1)
      eri(3,2,1,2)= eri(3,2,2,1)
      eri(1,3,1,2)= eri(1,3,2,1)
      eri(2,3,1,2)= eri(2,3,2,1)
      eri(3,3,1,2)= eri(3,3,2,1)
      eri(1,1,2,2)= r4( 4)+r3(4,3)*qx+r3(4,4)*qx+r2(1,2)+r2(2,7)*xx+r2(2,8)+r1(1,3)*qx &
&                  +r1(1,4)*qx+r0(3)*xx+r0(4)
      eri(2,1,2,2)= r4( 7)+r3(7,4)*qx+r2(4,2)+r1(2,4)*qx
      eri(3,1,2,2)= r4( 8)+r3(4,3)*qz+r3(8,4)*qx+r2(5,2)+r2(2,7)*xz+r1(1,3)*qz+r1(3,4)*qx &
&                  +r0(3)*xz
      eri(1,2,2,2)= r4( 7)+r3(7,3)*qx+r2(4,2)+r1(2,3)*qx
      eri(2,2,2,2)= r4(11)+r2(2,2)+r2(2,8)+r0(4)
      eri(3,2,2,2)= r4(12)+r3(7,3)*qz+r2(6,2)+r1(2,3)*qz
      eri(1,3,2,2)= r4( 8)+r3(8,3)*qx+r3(4,4)*qz+r2(5,2)+r2(2,7)*xz+r1(3,3)*qx+r1(1,4)*qz &
&                  +r0(3)*xz
      eri(2,3,2,2)= r4(12)+r3(7,4)*qz+r2(6,2)+r1(2,4)*qz
      eri(3,3,2,2)= r4(13)+r3(8,3)*qz+r3(8,4)*qz+r2(3,2)+r2(2,7)*zz+r2(2,8)+r1(3,3)*qz &
&                  +r1(3,4)*qz+r0(3)*zz+r0(4)
      eri(1,1,3,2)= r4( 5)-r3(2,1)+r3(5,3)*qx+r3(5,4)*qx-r2(4,3)*qx-r2(4,4)*qx+r2(6,7)*xx &
&                  +r2(6,8)-r1(2,5)*xx-r1(2,6)
      eri(2,1,3,2)= r4( 8)-r3(4,1)+r3(8,4)*qx-r2(2,4)*qx
      eri(3,1,3,2)= r4( 9)-r3(5,1)+r3(5,3)*qz+r3(9,4)*qx-r2(4,3)*qz-r2(6,4)*qx+r2(6,7)*xz &
&                  -r1(2,5)*xz
      eri(1,2,3,2)= r4( 8)-r3(4,1)+r3(8,3)*qx-r2(2,3)*qx
      eri(2,2,3,2)= r4(12)-r3(7,1)+r2(6,8)-r1(2,6)
      eri(3,2,3,2)= r4(13)-r3(8,1)+r3(8,3)*qz-r2(2,3)*qz
      eri(1,3,3,2)= r4( 9)-r3(5,1)+r3(9,3)*qx+r3(5,4)*qz-r2(6,3)*qx-r2(4,4)*qz+r2(6,7)*xz &
&                  -r1(2,5)*xz
      eri(2,3,3,2)= r4(13)-r3(8,1)+r3(8,4)*qz-r2(2,4)*qz
      eri(3,3,3,2)= r4(14)-r3(9,1)+r3(9,3)*qz+r3(9,4)*qz-r2(6,3)*qz-r2(6,4)*qz+r2(6,7)*zz &
&                  +r2(6,8)-r1(2,5)*zz-r1(2,6)
      eri(1,1,1,3)= r4( 3)-r3(1,2)+r3(3,3)*qx+r3(3,4)*qx-r2(1,5)*qx-r2(1,6)*qx+r2(5,7)*xx &
&                  +r2(5,8)-r1(1,7)*xx-r1(1,8)
      eri(2,1,1,3)= r4( 5)-r3(2,2)+r3(5,4)*qx-r2(4,6)*qx
      eri(3,1,1,3)= r4( 6)-r3(3,2)+r3(3,3)*qz+r3(6,4)*qx-r2(1,5)*qz-r2(5,6)*qx+r2(5,7)*xz &
&                  -r1(1,7)*xz
      eri(1,2,1,3)= r4( 5)-r3(2,2)+r3(5,3)*qx-r2(4,5)*qx
      eri(2,2,1,3)= r4( 8)-r3(4,2)+r2(5,8)-r1(1,8)
      eri(3,2,1,3)= r4( 9)-r3(5,2)+r3(5,3)*qz-r2(4,5)*qz
      eri(1,3,1,3)= r4( 6)-r3(3,2)+r3(6,3)*qx+r3(3,4)*qz-r2(5,5)*qx-r2(1,6)*qz+r2(5,7)*xz &
&                  -r1(1,7)*xz
      eri(2,3,1,3)= r4( 9)-r3(5,2)+r3(5,4)*qz-r2(4,6)*qz
      eri(3,3,1,3)= r4(10)-r3(6,2)+r3(6,3)*qz+r3(6,4)*qz-r2(5,5)*qz-r2(5,6)*qz+r2(5,7)*zz &
&                  +r2(5,8)-r1(1,7)*zz-r1(1,8)
      eri(1,1,2,3)= r4( 5)-r3(2,2)+r3(5,3)*qx+r3(5,4)*qx-r2(4,5)*qx-r2(4,6)*qx+r2(6,7)*xx &
&                  +r2(6,8)-r1(2,7)*xx-r1(2,8)
      eri(2,1,2,3)= r4( 8)-r3(4,2)+r3(8,4)*qx-r2(2,6)*qx
      eri(3,1,2,3)= r4( 9)-r3(5,2)+r3(5,3)*qz+r3(9,4)*qx-r2(4,5)*qz-r2(6,6)*qx+r2(6,7)*xz &
&                  -r1(2,7)*xz
      eri(1,2,2,3)= r4( 8)-r3(4,2)+r3(8,3)*qx-r2(2,5)*qx
      eri(2,2,2,3)= r4(12)-r3(7,2)+r2(6,8)-r1(2,8)
      eri(3,2,2,3)= r4(13)-r3(8,2)+r3(8,3)*qz-r2(2,5)*qz
      eri(1,3,2,3)= r4( 9)-r3(5,2)+r3(9,3)*qx+r3(5,4)*qz-r2(6,5)*qx-r2(4,6)*qz+r2(6,7)*xz &
&                  -r1(2,7)*xz
      eri(2,3,2,3)= r4(13)-r3(8,2)+r3(8,4)*qz-r2(2,6)*qz
      eri(3,3,2,3)= r4(14)-r3(9,2)+r3(9,3)*qz+r3(9,4)*qz-r2(6,5)*qz-r2(6,6)*qz+r2(6,7)*zz &
&                  +r2(6,8)-r1(2,7)*zz-r1(2,8)
      eri(1,1,3,3)= r4( 6)-r3(3,1)-r3(3,2)+r3(6,3)*qx+r3(6,4)*qx+r2(1,1)+r2(1,2) &
&                  -r2(5,3)*qx-r2(5,4)*qx-r2(5,5)*qx-r2(5,6)*qx+r2(3,7)*xx+r2(3,8) &
&                  +r1(1,1)*qx+r1(1,2)*qx+r1(1,3)*qx+r1(1,4)*qx-r1(3,5)*xx-r1(3,6) &
&                  -r1(3,7)*xx-r1(3,8)+r0(1)*xx+r0(2)+r0(3)*xx+r0(4)
      eri(2,1,3,3)= r4( 9)-r3(5,1)-r3(5,2)+r3(9,4)*qx+r2(4,1)+r2(4,2)-r2(6,4)*qx &
&                  -r2(6,6)*qx+r1(2,2)*qx+r1(2,4)*qx
      eri(3,1,3,3)= r4(10)-r3(6,1)-r3(6,2)+r3(6,3)*qz+r3(10,4)*qx+r2(5,1)+r2(5,2) &
&                  -r2(5,3)*qz-r2(3,4)*qx-r2(5,5)*qz-r2(3,6)*qx+r2(3,7)*xz+r1(1,1)*qz &
&                  +r1(3,2)*qx+r1(1,3)*qz+r1(3,4)*qx-r1(3,5)*xz-r1(3,7)*xz+r0(1)*xz+r0(3)*xz
      eri(1,2,3,3)= r4( 9)-r3(5,1)-r3(5,2)+r3(9,3)*qx+r2(4,1)+r2(4,2)-r2(6,3)*qx &
&                  -r2(6,5)*qx+r1(2,1)*qx+r1(2,3)*qx
      eri(2,2,3,3)= r4(13)-r3(8,1)-r3(8,2)+r2(2,1)+r2(2,2)+r2(3,8)-r1(3,6)-r1(3,8)+r0(2)+r0(4)
      eri(3,2,3,3)= r4(14)-r3(9,1)-r3(9,2)+r3(9,3)*qz+r2(6,1)+r2(6,2)-r2(6,3)*qz &
&                  -r2(6,5)*qz+r1(2,1)*qz+r1(2,3)*qz
      eri(1,3,3,3)= r4(10)-r3(6,1)-r3(6,2)+r3(10,3)*qx+r3(6,4)*qz+r2(5,1)+r2(5,2) &
&                  -r2(3,3)*qx-r2(5,4)*qz-r2(3,5)*qx-r2(5,6)*qz+r2(3,7)*xz+r1(3,1)*qx &
&                  +r1(1,2)*qz+r1(3,3)*qx+r1(1,4)*qz-r1(3,5)*xz-r1(3,7)*xz+r0(1)*xz+r0(3)*xz
      eri(2,3,3,3)= r4(14)-r3(9,1)-r3(9,2)+r3(9,4)*qz+r2(6,1)+r2(6,2)-r2(6,4)*qz &
&                  -r2(6,6)*qz+r1(2,2)*qz+r1(2,4)*qz
      eri(3,3,3,3)= r4(15)-r3(10,2)-r3(10,1)+r3(10,3)*qz+r3(10,4)*qz+r2(3,1)+r2(3,2) &
&                  -r2(3,3)*qz-r2(3,4)*qz-r2(3,5)*qz-r2(3,6)*qz+r2(3,7)*zz+r2(3,8) &
&                  +r1(3,1)*qz+r1(3,2)*qz+r1(3,3)*qz+r1(3,4)*qz-r1(3,5)*zz-r1(3,6) &
&                  -r1(3,7)*zz-r1(3,8)+r0(1)*zz+r0(2)+r0(3)*zz+r0(4)
!
      do j= 1,3
        do k= 1,3
          do l= 1,3
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
          do l= 1,3
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
          do l= 1,3
            work(1)= eri(l,1,j,i)
            work(2)= eri(l,2,j,i)
            work(3)= eri(l,3,j,i)
            eri(l,1,j,i)= work(1)*rot(1,1)+work(2)*rot(2,1)+work(3)*rot(3,1)
            eri(l,2,j,i)= work(1)*rot(1,2)+work(2)*rot(2,2)+work(3)*rot(3,2)
            eri(l,3,j,i)= work(1)*rot(1,3)+work(2)*rot(2,3)+work(3)*rot(3,3)
          enddo
        enddo
      enddo
      do i= 1,3
        do j= 1,3
          do k= 1,3
            work(1)= eri(1,k,j,i)
            work(2)= eri(2,k,j,i)
            work(3)= eri(3,k,j,i)
            eri(1,k,j,i)= work(1)*rot(1,1)+work(2)*rot(2,1)+work(3)*rot(3,1)
            eri(2,k,j,i)= work(1)*rot(1,2)+work(2)*rot(2,2)+work(3)*rot(3,2)
            eri(3,k,j,i)= work(1)*rot(1,3)+work(2)*rot(2,3)+work(3)*rot(3,3)
          enddo
        enddo
      enddo
!
      do i= 1,3
        do j= 1,3
          do k= 1,3
            do l= 1,3
              phmdint(l,k,j,i)= eri(l,k,j,i)
            enddo
          enddo
        enddo
      enddo
      return
end
