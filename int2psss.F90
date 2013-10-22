!------------------------------------------------------------------
  subroutine int2psss(phmdint,exfac1,exfac2,xyziq,xzkl,rot,nijkl)
!------------------------------------------------------------------
!
! Calculate (ps|ss) integrals
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
      integer :: ij, kl, igrid
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, sqrtpi4=0.8862269254527580D+00
      real(8),parameter :: two=2.0D+00, half=0.5D+00
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*), xzkl(2), rot(3,3)
      real(8),intent(out) :: phmdint(6,6,6,6)
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
        do ij= 1,nijkl(1)
          ex12= exfac1(1,ij)
          zip = exfac1(4,ij)
          c12 = exfac1(5,ij)
          ex14= ex12+ex34
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
            ft(0)= fgrid(0,0,igrid)      +fgrid( 1,0,igrid)*tval  +fgrid( 2,0,igrid)*tval2 &
&                 +fgrid(3,0,igrid)*tval3+fgrid( 4,0,igrid)*tval4 +fgrid( 5,0,igrid)*tval5 &
&                 +fgrid(6,0,igrid)*tval6+fgrid( 7,0,igrid)*tval7 +fgrid( 8,0,igrid)*tval8 &
&                 +fgrid(9,0,igrid)*tval9+fgrid(10,0,igrid)*tval10
            ft(1)= fgrid(0,1,igrid)      +fgrid( 1,1,igrid)*tval  +fgrid( 2,1,igrid)*tval2 &
&                 +fgrid(3,1,igrid)*tval3+fgrid( 4,1,igrid)*tval4 +fgrid( 5,1,igrid)*tval5 &
&                 +fgrid(6,1,igrid)*tval6+fgrid( 7,1,igrid)*tval7 +fgrid( 8,1,igrid)*tval8 &
&                 +fgrid(9,1,igrid)*tval9+fgrid(10,1,igrid)*tval10
            fac= c12*sqrt(ex41)
            ft(0)= ft(0)*fac
            ft(1)= ft(1)*fac*expq*two
          endif
          f0= f0+ft(0)
          f1(1)= f1(1)-ft(1)
          f1(2)= f1(2)-ft(1)*zpq
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
