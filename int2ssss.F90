!---------------------------------------------------------
  subroutine int2ssss(phmdint,exfac1,exfac2,xyziq,nijkl)
!---------------------------------------------------------
!
! Calculate (ss|ss) integral
!
! In  : exfac12, exfac34 (exponents and coefficients of primitive pair functions)
!       xyziq    (x,y,z length of i and q)
!       nij, nkl (number of primitive pair functions)
! Out : phmdint  (two-electron integral)
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nijkl(2)
      integer :: ij, kl, igrid
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, sqrtpi4=0.8862269254527580D+00
      real(8),intent(in) :: exfac1(5,*), exfac2(5,*)
      real(8),intent(in) :: xyziq(3,*)
      real(8),intent(out) :: phmdint(6,6,6,6)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: tinv, f0, ft0, r0
      real(8) :: ex12, ex34, ex41, expq, c12, c34, zip, xiq, yiq, ziq, xypq2, zpq, zpq2
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
            ft0= c12*sqrtpi4*sqrt(ex41*tval)*tinv
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
            ft0= fgrid(0,0,igrid)      +fgrid( 1,0,igrid)*tval  +fgrid( 2,0,igrid)*tval2 &
&               +fgrid(3,0,igrid)*tval3+fgrid( 4,0,igrid)*tval4 +fgrid( 5,0,igrid)*tval5 &
&               +fgrid(6,0,igrid)*tval6+fgrid( 7,0,igrid)*tval7 +fgrid( 8,0,igrid)*tval8 &
&               +fgrid(9,0,igrid)*tval9+fgrid(10,0,igrid)*tval10
            ft0= ft0*c12*sqrt(ex41)
          endif
          f0= f0+ft0
        enddo
        r0= r0+f0*c34
      enddo
!
      phmdint(1,1,1,1)= r0
      return
end
