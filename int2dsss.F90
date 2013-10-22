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
      integer :: ij, kl, igrid, k, l, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, sqrtpi4=0.8862269254527580D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00
      real(8),parameter :: sqrt3=1.73205080756888D+00, sqrt3h=8.660254037844386D-01
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*), xzkl(2), rot(3,3)
      real(8),intent(out) :: phmdint(6,6,6,6)
      real(8) :: tval, tval1, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:2)
      real(8) :: f0, f1(2), f2(3), r0(2), r1(3), r2(6)
      real(8) :: ex12, ex34, ex43, ex41, expq, ex3q, ex4q, c12, c34, zip, xiq, yiq, ziq
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
            do ii= 0,2
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
            fac= c12*sqrt(ex41)
            expq= expq*two
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
      if(nbfijkl(4) == 6)then
        do l= 1,6
          phmdint(l,1,1,1)= eri(1)*rot2(1,l)+eri(2)*rot2(2,l)+eri(3)*rot2(3,l) &
&                          +eri(4)*rot2(4,l)+eri(5)*rot2(5,l)+eri(6)*rot2(6,l) 
        enddo
      else
        do l= 1,6
          rot3(l,1)= rot2(l,3)-(rot2(l,1)+rot2(l,2))*half
          rot3(l,2)= rot2(l,5)
          rot3(l,3)= rot2(l,6)
          rot3(l,4)=(rot2(l,1)-rot2(l,2))*sqrt3h
          rot3(l,5)= rot2(l,4)
        enddo
        do l= 1,5
          phmdint(l,1,1,1)= eri(1)*rot3(1,l)+eri(2)*rot3(2,l)+eri(3)*rot3(3,l) &
&                          +eri(4)*rot3(4,l)+eri(5)*rot3(5,l)+eri(6)*rot3(6,l) 
        enddo
      endif
!
      return
end
