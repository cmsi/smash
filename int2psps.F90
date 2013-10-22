!------------------------------------------------------------------
  subroutine int2psps(phmdint,exfac1,exfac2,xyziq,xzkl,rot,nijkl)
!------------------------------------------------------------------
!
! Calculate (ps|ps) integrals
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
      integer :: ij, kl, igrid, ii, j, l
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, sqrtpi4=0.8862269254527580D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*), xzkl(2), rot(3,3)
      real(8),intent(out) :: phmdint(6,6,6,6)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, ft(0:2)
      real(8) :: f0, f1(4), f2(3), r0, r1(3,2), r2(6)
      real(8) :: ex12, ex34, ex43, ex41, expq, ex3q, c12, c34, zip, zjp, xiq, yiq, ziq
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
          f0   = f0   +ft(0)*zjp
          f1(1)= f1(1)-ft(1)*zjp
          f1(2)= f1(2)-ft(1)*zjp*zpq
          f1(3)= f1(3)-ft(1)*pmd
          f1(4)= f1(4)-ft(1)*pmd*zpq
          f2(1)= f2(1)+ft(2)*pmd
          f2(2)= f2(2)+ft(2)*pmd*zpq
          f2(3)= f2(3)+ft(2)*pmd*zpq2
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
