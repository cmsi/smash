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
      integer :: ij, kl, igrid, i, j, k, l, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, sqrtpi4=0.8862269254527580D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00, five=5.0D+00
      real(8),parameter :: sqrt3=1.73205080756888D+00, sqrt3h=8.660254037844386D-01
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*), xzkl(2), rot(3,3)
      real(8),intent(out) :: phmdint(6,6,6,6)
      real(8) :: tval, tval1, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:3)
      real(8) :: f0, f1(2,2), f2(3,2), f3(4), r0(2), r1(3,3), r2(6,2), r3(10)
      real(8) :: ex12, ex34, ex43, ex41, expq, expq2, ex3q, ex4q, c12, c34, zip, xiq, yiq, ziq
      real(8) :: xiq2, yiq2, xyiq, xypq2, zpq, zpq2, fac, qmd, qmd2, qmd2x, qmd2y, pmd, zjp
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
            do ii= 0,3
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
        do j= 1,3
          do l= 1,6
            phmdint(l,1,j,1)= eri(1,j)*rot2(1,l)+eri(2,j)*rot2(2,l)+eri(3,j)*rot2(3,l) &
&                            +eri(4,j)*rot2(4,l)+eri(5,j)*rot2(5,l)+eri(6,j)*rot2(6,l) 
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
        do j= 1,3
          do l= 1,5
            phmdint(l,1,j,1)= eri(1,j)*rot3(1,l)+eri(2,j)*rot3(2,l)+eri(3,j)*rot3(3,l) &
&                            +eri(4,j)*rot3(4,l)+eri(5,j)*rot3(5,l)+eri(6,j)*rot3(6,l) 
          enddo
        enddo
      endif
      return
end
