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
!----------------------------------------------------------------------------------------
  subroutine int2elec(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                     mxprsh,threshex)
!----------------------------------------------------------------------------------------
!
! Wrapper of two-electron integral routines
!
! In  : exijkl    (Exponents of basis functions)
!       coijkl    (Coefficients of basis functions)
!       xyzijkl   (x,y,z coordinates)
!       nprimijkl (Numbers of primitive functions)
!       nangijkl  (Degrees of angular momentum)
!       nbfijkl   (Numbers of basis functions)
!       maxdim    (Dimension of two-electron integral array)
!       mxprsh    (Size of primitive function array)
!       threshex  (Threshold of exponential calculation)
! Out : twoeri    (Two-electron integrals
!                  The values are stored in order of (lsh,ksh,jsh,ish).)
!
      implicit none
      integer,intent(in) :: nprimijkl(4), nangijkl(4), nbfijkl(4), maxdim, mxprsh
      real(8),intent(in) :: exijkl(mxprsh,4), coijkl(mxprsh,4), xyzijkl(3,4), threshex
      real(8),intent(out) :: twoeri(maxdim,maxdim,maxdim,maxdim)
!
      if((nangijkl(1)<=2).and.(nangijkl(2)<=2).and.(nangijkl(3)<=2).and.(nangijkl(4)<=2)) then
!
! Pople-Hehre and McMurchie-Davidson scheme
!
        call int2phmd(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                     mxprsh,threshex)
      elseif((nangijkl(1)<=6).and.(nangijkl(2)<=6).and.(nangijkl(3)<=6).and.(nangijkl(4)<=6)) then
!
! Rys quadrature scheme
!
        call int2rys(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                    mxprsh,threshex)
      else
        write(*,'(" Error! Subroutine int2elec supports up to g function.")')
      endif
!
      return
end


!----------------------------------------------------------------------------------------
  subroutine int2phmd(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                     mxprsh,threshex)
!----------------------------------------------------------------------------------------
!
! Driver of two-electron integrals from (ss|ss) to (dd|dd)
! using McMurchie-Davidson and Pople-Hehre mathods
!
! In  : exijkl    (Exponents of basis functions)
!       coijkl    (Coefficients of basis functions)
!       xyzijkl   (x,y,z coordinates)
!       nprimijkl (Numbers of primitive functions)
!       nangijkl  (Degrees of angular momentum)
!       nbfijkl   (Numbers of basis functions)
!       maxdim    (Dimension of two-electron integral array)
!       mxprsh    (Size of primitive function array)
!       threshex  (Threshold of exponential calculation)
! Out : twoeri    (Two-electron integrals
!                  The values are stored in order of (lsh,ksh,jsh,ish).)
!
      implicit none
      integer,intent(in) :: nprimijkl(4), nangijkl(4), nbfijkl(4), maxdim, mxprsh
      integer,parameter :: mxprsh2=30
      integer :: nangtotal, inttype, intorder
      integer :: iprim, jprim, kprim, lprim, i, j, k, l, nijkl(2)
      integer :: intijkl(4,8), inttable(81,2), ii, jj, kk, ll, nbfijkl2(4)
      real(8),parameter :: zero=0.0d+00, one=1.0D+00, half=0.5D+00
      real(8),parameter :: pi52=0.3498683665524973D+02 ! 2.0*pi**2.5
      real(8),intent(in) :: exijkl(mxprsh,4), coijkl(mxprsh,4), xyzijkl(3,4), threshex
      real(8),intent(out) :: twoeri(maxdim,maxdim,maxdim,maxdim)
      real(8) :: phmdint(6,6,6,6)
      real(8) :: r12, r34, rij, rkl, tmp, rot(3,3), veckl(3), cosph, sinph, abscos
      real(8) :: xyzij(3), xyzkl(3), xyzik(3), xzkl(2)
      real(8) :: exi, exj, ci, ex12, ex21, r12exi, rijexi, rijexj, cij, r12ex
      real(8) :: exk, exl, ck, ex34, r34exk, ckl, r34ex
      real(8) :: exfac(5,mxprsh2*mxprsh2*2), xyziq(3,mxprsh2*mxprsh2)
      real(8) :: work(mxprsh2*mxprsh2*6)
      data intijkl/1,2,3,4, 2,1,3,4, 1,2,4,3, 3,4,1,2, 2,1,4,3, 3,4,2,1, 4,3,1,2, 4,3,2,1/
      data inttable/ 1, 2, 7, 2, 3, 8, 7, 8,10, 2, 4, 9, 4, 5,11, 9,11,14, 7, 9,&
&                   12, 9,13,15,12,15,17, 2, 4, 9, 4, 5,11, 9,11,14, 3, 5,13, 5,&
&                    6,16,13,16,18, 8,11,15,11,16,19,15,19,20, 7, 9,12, 9,13,15,&
&                   12,15,17, 8,11,15,11,16,19,15,19,20,10,14,17,14,18,20,17,20,21,&
&                   1,1,1,3,1,1,3,3,1,4,1,1,3,1,1,3,3,1,4,4,&
&                   1,7,4,1,3,3,1,6,2,2,5,2,2,5,5,2,4,4,1,7,&
&                   1,1,3,3,1,4,4,4,7,4,1,7,3,1,6,6,2,8,6,2,&
&                   5,5,2,6,6,6,8,6,2,8,5,2,4,4,4,7,4,4,7,7,1/
!
      if(mxprsh > mxprsh2) then
        write(*,'(" Error! Parameter mxprsh2 in int2phmd is small!")')
        call exit
      endif
!
      nangtotal= nangijkl(1)*27+nangijkl(2)*9+nangijkl(3)*3+nangijkl(4)+1
      inttype  = inttable(nangtotal,1)
      intorder = inttable(nangtotal,2)
!
      ii= intijkl(1,intorder)
      jj= intijkl(2,intorder)
      kk= intijkl(3,intorder)
      ll= intijkl(4,intorder)
!
      do i= 1,3
        xyzij(i)= xyzijkl(i,jj)-xyzijkl(i,ii)
        xyzkl(i)= xyzijkl(i,ll)-xyzijkl(i,kk)
        xyzik(i)= xyzijkl(i,kk)-xyzijkl(i,ii)
      enddo
      r12= xyzij(1)*xyzij(1)+xyzij(2)*xyzij(2)+xyzij(3)*xyzij(3)
      r34= xyzkl(1)*xyzkl(1)+xyzkl(2)*xyzkl(2)+xyzkl(3)*xyzkl(3)
!
! Calculate z-axis of Pople-Hehre algorithm
!
      if(r12.ne.zero) then
        rij= sqrt(r12)
        tmp= one/rij
        rot(1,3)= xyzij(1)*tmp
        rot(2,3)= xyzij(2)*tmp
        rot(3,3)= xyzij(3)*tmp
      else
        rij= zero
        rot(1,3)= zero
        rot(2,3)= zero
        rot(3,3)= one
      endif
!
! Calculate KL vector
!
      if(r34.ne.zero) then
        rkl= sqrt(r34)
        tmp= one/rkl
        veckl(1)= xyzkl(1)*tmp
        veckl(2)= xyzkl(2)*tmp
        veckl(3)= xyzkl(3)*tmp
      else
        rkl= zero
        veckl(1)= zero
        veckl(2)= zero
        veckl(3)= one
      endif
!
! Calculate y-axis of Pople-Hehre algorithm
!
      rot(1,2)= veckl(3)*rot(2,3)-veckl(2)*rot(3,3)
      rot(2,2)= veckl(1)*rot(3,3)-veckl(3)*rot(1,3)
      rot(3,2)= veckl(2)*rot(1,3)-veckl(1)*rot(2,3)
!
      cosph= veckl(1)*rot(1,3)+veckl(2)*rot(2,3)+veckl(3)*rot(3,3)
      if(cosph.gt.one)then
        cosph=one
      elseif(cosph.lt.-one)then
        cosph=-one
      endif
      abscos=abs(cosph)
      if(abscos.gt.0.9D+00) then
        sinph= sqrt(rot(1,2)*rot(1,2)+rot(2,2)*rot(2,2)+rot(3,2)*rot(3,2))
      else
        sinph= sqrt(one-cosph*cosph)
      endif
      if((abscos.le.0.9D+00).or.(sinph.ge.1.0D-12)) then
        tmp= one/sinph
        rot(1,2)= rot(1,2)*tmp
        rot(2,2)= rot(2,2)*tmp
        rot(3,2)= rot(3,2)*tmp
      else
        i= 3
        if(abs(rot(1,3)).le.0.7D+00) i= 1
        tmp= rot(i,3)*rot(i,3)
        if(tmp.ge.one)then
          tmp=zero
        else
          tmp= one/sqrt(one-tmp)
        endif
        if(abs(rot(1,3)).le.0.7D+00) then
          rot(1,2)= zero
          rot(2,2)= rot(3,3)*tmp
          rot(3,2)=-rot(2,3)*tmp
        else
          rot(1,2)= rot(2,3)*tmp
          rot(2,2)=-rot(1,3)*tmp
          rot(3,2)= zero
        endif
      endif
!
! Calculate x-axis of Pople-Hehre algorithm
!
      rot(1,1)= rot(2,2)*rot(3,3)-rot(3,2)*rot(2,3)
      rot(2,1)= rot(3,2)*rot(1,3)-rot(1,2)*rot(3,3)
      rot(3,1)= rot(1,2)*rot(2,3)-rot(2,2)*rot(1,3)
!
      veckl(1)= xyzik(1)*rot(1,1)+xyzik(2)*rot(2,1)+xyzik(3)*rot(3,1)
      veckl(2)= xyzik(1)*rot(1,2)+xyzik(2)*rot(2,2)+xyzik(3)*rot(3,2)
      veckl(3)= xyzik(1)*rot(1,3)+xyzik(2)*rot(2,3)+xyzik(3)*rot(3,3)
      xyzik(1)= veckl(1)
      xyzik(2)= veckl(2)
      xyzik(3)= veckl(3)
!
      xzkl(1)= rkl*sinph
      xzkl(2)= rkl*cosph
!
      nijkl(:)= 0
      do iprim= 1,nprimijkl(ii)
        exi= exijkl(iprim,ii)
        ci = coijkl(iprim,ii)*pi52
        r12exi= exi*r12
        rijexi=-exi*rij
        do jprim= 1,nprimijkl(jj)
          exj = exijkl(jprim,jj)
          cij = ci*coijkl(jprim,jj)
          ex12= exi+exj
          rijexj= rij*exj
          r12ex= r12exi*exj
          if(r12ex.gt.threshex*ex12)cycle
          nijkl(1)= nijkl(1)+1
          exfac(1,nijkl(1))= ex12
          exfac(3,nijkl(1))= rijexi
          exfac(4,nijkl(1))= rijexj
          exfac(5,nijkl(1))= cij
          work(nijkl(1))=-r12ex
        enddo
      enddo
!
      do kprim= 1,nprimijkl(kk)
        exk= exijkl(kprim,kk)
        ck = coijkl(kprim,kk)
        r34exk= r34*exk
        do lprim= 1,nprimijkl(ll)
          exl = exijkl(lprim,ll)
          ckl = ck*coijkl(lprim,ll)
          ex34= exk+exl
          r34ex= r34exk*exl
          if(r34ex.gt.threshex*ex34)cycle
          nijkl(2)= nijkl(2)+1
          exfac(1,nijkl(1)+nijkl(2))= ex34
          exfac(3,nijkl(1)+nijkl(2))=-exk
          exfac(4,nijkl(1)+nijkl(2))= exl
          exfac(5,nijkl(1)+nijkl(2))= ckl
          work(nijkl(1)+nijkl(2))=-r34ex
        enddo
      enddo
!
      do i=1,nijkl(1)+nijkl(2)
        ex21= one/exfac(1,i)
        exfac(2,i)= ex21*half
        exfac(3,i)= exfac(3,i)*ex21
        exfac(4,i)= exfac(4,i)*ex21
        exfac(5,i)= exfac(5,i)*ex21*exp(work(i)*ex21)
      enddo
      do i=1,nijkl(2)
        xyziq(1,i)= xyzik(1)+xzkl(1)*exfac(4,nijkl(1)+i)
        xyziq(2,i)= xyzik(2)
        xyziq(3,i)= xyzik(3)+xzkl(2)*exfac(4,nijkl(1)+i)
      enddo
!
! Transpose rotational matrix
!
      do j= 1,2
        do i= j+1,3
          tmp= rot(i,j)
          rot(i,j)= rot(j,i)
          rot(j,i)= tmp
        enddo
      enddo
!
      nbfijkl2(1)= nbfijkl(ii)
      nbfijkl2(2)= nbfijkl(jj)
      nbfijkl2(3)= nbfijkl(kk)
      nbfijkl2(4)= nbfijkl(ll)
!
      select case(inttype)
        case (1)
          call int2ssss(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,work,nijkl)
        case (2)
          call int2psss(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,work,nijkl)
        case (3)
          call int2ppss(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,work,nijkl)
        case (4)
          call int2psps(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,work,nijkl)
        case (5)
          call int2ppps(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl)
        case (6)
          call int2pppp(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl)
        case (7)
          call int2dsss(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl,nbfijkl2)
        case (8)
          call int2dpss(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl,nbfijkl2)
        case (9)
          call int2dsps(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl,nbfijkl2)
        case (10)
          call int2ddss(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl,nbfijkl2)
        case (11)
          call int2dpps(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl,nbfijkl2)
        case (12)
          call int2dsds(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl,nbfijkl2)
        case (13)
          call int2dspp(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl,nbfijkl2)
        case (14)
          call int2ddps(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl,nbfijkl2)
        case (15)
          call int2dpds(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl,nbfijkl2)
        case (16)
          call int2dppp(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl,nbfijkl2)
        case (17)
          call int2ddds(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl,nbfijkl2)
        case (18)
          call int2ddpp(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl,nbfijkl2)
        case (19)
          call int2dpdp(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl,nbfijkl2)
        case (20)
          call int2dddp(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl,nbfijkl2)
        case (21)
          call int2dddd(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl,nbfijkl2)
      end select
!
      select case(intorder)
        case(1)
          do i= 1,nbfijkl2(1)
            do j= 1,nbfijkl2(2)
              do k= 1,nbfijkl2(3)
                do l= 1,nbfijkl2(4)
                  twoeri(l,k,j,i)= phmdint(l,k,j,i)
                enddo
              enddo
            enddo
          enddo
        case(2)
          do j= 1,nbfijkl2(1)
            do i= 1,nbfijkl2(2)
              do k= 1,nbfijkl2(3)
                do l= 1,nbfijkl2(4)
                  twoeri(l,k,j,i)= phmdint(l,k,i,j)
                enddo
              enddo
            enddo
          enddo
        case(3)
        do i= 1,nbfijkl2(1)
          do j= 1,nbfijkl2(2)
            do l= 1,nbfijkl2(3)
              do k= 1,nbfijkl2(4)
                  twoeri(l,k,j,i)= phmdint(k,l,j,i)
                enddo
              enddo
            enddo
          enddo
        case(4)
          do k= 1,nbfijkl2(1)
            do l= 1,nbfijkl2(2)
              do i= 1,nbfijkl2(3)
                do j= 1,nbfijkl2(4)
                  twoeri(l,k,j,i)= phmdint(j,i,l,k)
                enddo
              enddo
            enddo
          enddo
        case(5)
          do j= 1,nbfijkl2(1)
            do i= 1,nbfijkl2(2)
              do l= 1,nbfijkl2(3)
                do k= 1,nbfijkl2(4)
                  twoeri(l,k,j,i)= phmdint(k,l,i,j)
                enddo
              enddo
            enddo
          enddo
        case(6)
          do k= 1,nbfijkl2(1)
            do l= 1,nbfijkl2(2)
              do j= 1,nbfijkl2(3)
                do i= 1,nbfijkl2(4)
                  twoeri(l,k,j,i)= phmdint(i,j,l,k)
                enddo
              enddo
            enddo
          enddo
        case(7)
          do l= 1,nbfijkl2(1)
            do k= 1,nbfijkl2(2)
              do i= 1,nbfijkl2(3)
                do j= 1,nbfijkl2(4)
                  twoeri(l,k,j,i)= phmdint(j,i,k,l)
                enddo
              enddo
            enddo
          enddo
        case(8)
          do l= 1,nbfijkl2(1)
            do k= 1,nbfijkl2(2)
              do j= 1,nbfijkl2(3)
                do i= 1,nbfijkl2(4)
                  twoeri(l,k,j,i)= phmdint(i,j,k,l)
                enddo
              enddo
            enddo
          enddo
      end select
!
      return
end


!---------------------------------------------------------------------------------------
  subroutine int2rys(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                    mxprsh,threshex)
!---------------------------------------------------------------------------------------
!
! Calculate two-electron integrals using Rys quadrature
!
! In  : exijkl    (Exponents of basis functions)
!       coijkl    (Coefficients of basis functions)
!       xyzijkl   (x,y,z coordinates)
!       nprimijkl (Numbers of primitive functions)
!       nangijkl  (Degrees of angular momentum)
!       nbfijkl   (Numbers of basis functions)
!       maxdim    (Dimension of two-electron integral array)
!       mxprsh    (Size of primitive function array)
!       threshex  (Threshold of exponential calculation)
! Out : twoeri    (Two-electron integrals
!                  The values are stored in order of (lsh,ksh,jsh,ish).)
!
      implicit none
      integer,intent(in) :: nprimijkl(4), nangijkl(4), nbfijkl(4), maxdim, mxprsh
      integer,parameter :: ncart(0:6)=(/1,3,6,10,15,21,28/), maxdim2=28, mxprsh2=30
      integer :: ncartijkl(4)
      real(8),intent(in) :: exijkl(mxprsh,4), coijkl(mxprsh,4), xyzijkl(3,4), threshex
      real(8),intent(out) :: twoeri(maxdim,maxdim,maxdim,maxdim)
      real(8) :: eritmp(maxdim2*maxdim2*maxdim2*maxdim2)
      real(8) :: ex12(3,mxprsh2*mxprsh2), ex34(3,mxprsh2*mxprsh2)
      real(8) :: xyza(3,mxprsh2*mxprsh2), xyzb(3,mxprsh2*mxprsh2)
      real(8) :: xyzaj(3,mxprsh2*mxprsh2), xyzbl(3,mxprsh2*mxprsh2)
!
      if(mxprsh > mxprsh2) then
        write(*,'(" Error! Parameter mxprsh2 in int2rys is small!")')
        call exit
      elseif(maxdim > maxdim2) then
        write(*,'(" Error! Parameter maxdim2 in int2rys is small!")')
        call exit
      endif
!
      ncartijkl(1)= ncart(nangijkl(1))
      ncartijkl(2)= ncart(nangijkl(2))
      ncartijkl(3)= ncart(nangijkl(3))
      ncartijkl(4)= ncart(nangijkl(4))
!
      call int2rys2(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                   mxprsh,threshex,eritmp,ex12,ex34,xyza,xyzb,xyzaj,xyzbl,ncartijkl)
      return
end


!------------------------------------------------------------------------------------------
  subroutine int2rys2(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                     mxprsh,threshex,eritmp,ex12,ex34,xyza,xyzb,xyzaj,xyzbl,ncartijkl)
!------------------------------------------------------------------------------------------
!
! Calculate two-electron integrals using Rys quadrature
!
      implicit none
      integer,intent(in) :: nprimijkl(4), nangijkl(4), nbfijkl(4), maxdim, mxprsh
      integer,intent(in) :: ncartijkl(4)
      integer :: nroots, mmax, nmax
      integer :: i, j, k, l, nij, nkl, ij, kl, iprim, jprim, kprim, lprim, iroot
      integer :: ii, jj, kk, ll, ix, iy, iz, jx, jy, jz, kx, ky, kz, lx, ly, lz, icount, nn
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, two=2.0D+00
      real(8),parameter :: three=3.0D+00, four=4.0D+00, five=5.0D+00, six=6.0D+00
      real(8),parameter :: eight=8.0D+00, p9=9.0D+00, ten=10.0D+00, twelve=12.0D+00
      real(8),parameter :: p15=15.0D+00, p16=16.0D+00, p18=18.0D+00, p20=20.0D+00 
      real(8),parameter :: p24=24.0D+00, p30=30.0D+00, p32=32.0D+00, p40=40.0D+00
      real(8),parameter :: p60=60.0D+00, p90=90.0D+00, p120=1.2D+02, p180=1.8D+02
      real(8),parameter :: eighth=0.125D+00, sixteenth=6.25D-02
      real(8),parameter :: pi52=0.3498683665524973D+02 ! 2.0*pi**2.5
      real(8),parameter :: sqrt3=1.732050807568877D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrt5=2.236067977499790D+00, sqrt15=3.872983346207417D+00
      real(8),parameter :: sqrt7=2.645751311064590D+00, sqrt35=5.916079783099616D+00
      real(8),parameter :: sqrt35third=3.415650255319866D+00
      real(8),parameter :: sqrt21=4.582575694955840D+00, sqrt63=7.937253933193772D+00
      real(8),parameter :: sqrt105=1.024695076595960D+01, sqrt11=3.316624790355400D+00
      real(8),parameter :: sqrt33=5.744562646538029D+00, sqrt99=9.949874371066200D+00
      real(8),parameter :: sqrt231=1.519868415357066D+01, sqrt231fifth=6.797058187186571D+00
      real(8),parameter :: sqrt385=1.962141687034858D+01
      real(8),parameter :: facf1=0.79056941504209483D+00 ! sqrt(5/2)/2
      real(8),parameter :: facf2=3.87298334620741688D+00 ! sqrt(15)
      real(8),parameter :: facf3=0.61237243569579452D+00 ! sqrt(3/2)/2
      real(8),parameter :: facf4=1.93649167310370844D+00 ! sqrt(15)/2
      real(8),parameter :: facg1=2.95803989154980802D+00 ! sqrt(35)/2
      real(8),parameter :: facg2=2.09165006633518887D+00 ! sqrt(35/2)/2
      real(8),parameter :: facg3=1.11803398874989484D+00 ! sqrt(5)/2
      real(8),parameter :: facg4=0.79056941504209483D+00 ! sqrt(5/2)/2
      real(8),parameter :: facg5=0.55901699437494742D+00 ! sqrt(5)/4
      real(8),parameter :: facg6=0.73950997288745200D+00 ! sqrt(35)/8
      real(8),parameter :: fach1=0.70156076002011400D+00 ! sqrt(63/2)/8
      real(8),parameter :: fach2=2.21852991866235601D+00 ! sqrt(315)/8
      real(8),parameter :: fach3=0.52291251658379721D+00 ! sqrt(35/2)/8
      real(8),parameter :: fach4=2.56173769148989959D+00 ! sqrt(105)/4
      real(8),parameter :: fach5=0.48412291827592711D+00 ! sqrt(15)/8
      real(8),parameter :: faci1=0.67169328938139615D+00 ! sqrt(231/2)/16
      real(8),parameter :: faci2=2.32681380862328561D+00 ! sqrt(693/2)/8
      real(8),parameter :: faci3=0.49607837082461073D+00 ! sqrt(63)/16
      real(8),parameter :: faci4=0.90571104663683991D+00 ! sqrt(105/2)/8
      real(8),parameter :: faci5=0.45285552331841995D+00 ! sqrt(105/2)/16
      real(8),parameter :: faci6=0.57282196186948000D+00 ! sqrt(21)/8
      real(8),intent(in) :: exijkl(mxprsh,4), coijkl(mxprsh,4), xyzijkl(3,4), threshex
      real(8),intent(out) :: twoeri(maxdim,maxdim,maxdim,maxdim)
      real(8),intent(out) :: eritmp(ncartijkl(4),ncartijkl(3),ncartijkl(2),ncartijkl(1))
      real(8),intent(out) :: ex12(3,*), xyza(3,*), xyzaj(3,*)
      real(8),intent(out) :: ex34(3,*), xyzb(3,*), xyzbl(3,*)
      real(8) :: xyzij(3), xyzkl(3), r12, r34, exi, exj, exk, exl, ci, cj, ck, cl
      real(8) :: xi, yi, zi, xk, yk, zk, exij, exji, exkl, exlk, r12ex, r34ex, dijkl
      real(8) :: ex21h, xyzab(3)
      real(8) :: ex41, ex41h, ex43h, b10t, bp01t, cx, cy, cz, cpx, cpy, cpz, tval
      real(8) :: trys(13), wrys(13), b00, b10, bp01, c00(3),cp00(3)
      real(8) :: xyzint(3,0:12,0:6,0:12), rysint(3,8,0:6,0:6,0:12,0:6,13), work(28)
!
      nroots=(nangijkl(1)+nangijkl(2)+nangijkl(3)+nangijkl(4))/2+1
      mmax= nangijkl(1)+nangijkl(2)
      nmax= nangijkl(3)+nangijkl(4)
      do i= 1,3
        xyzij(i)= xyzijkl(i,1)-xyzijkl(i,2)
        xyzkl(i)= xyzijkl(i,3)-xyzijkl(i,4)
      enddo
      r12= xyzij(1)*xyzij(1)+xyzij(2)*xyzij(2)+xyzij(3)*xyzij(3)
      r34= xyzkl(1)*xyzkl(1)+xyzkl(2)*xyzkl(2)+xyzkl(3)*xyzkl(3)
!
      do i= 1,ncartijkl(1)
        do j= 1,ncartijkl(2)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              eritmp(l,k,j,i)= zero
            enddo
          enddo
        enddo
      enddo
!
      nkl= 0
      do kprim= 1,nprimijkl(3)
        exk= exijkl(kprim,3)
        ck = coijkl(kprim,3)*pi52
        xk = exk*xyzijkl(1,3)
        yk = exk*xyzijkl(2,3)
        zk = exk*xyzijkl(3,3)
        do lprim= 1,nprimijkl(4)
          exl= exijkl(lprim,4)
          cl = coijkl(lprim,4)
          exkl= exk+exl
          exlk=one/exkl
          r34ex=exk*exl*exlk*r34
          if(r34ex >= threshex)cycle
          nkl= nkl+1
          ex34(1,nkl)= exkl
          ex34(2,nkl)= exlk*half
          ex34(3,nkl)= ck*cl*exp(-r34ex)*exlk
          xyzb(1,nkl)=(xk+exl*xyzijkl(1,4))*exlk
          xyzb(2,nkl)=(yk+exl*xyzijkl(2,4))*exlk
          xyzb(3,nkl)=(zk+exl*xyzijkl(3,4))*exlk
          xyzbl(1,nkl)= xyzb(1,nkl)-xyzijkl(1,4)
          xyzbl(2,nkl)= xyzb(2,nkl)-xyzijkl(2,4)
          xyzbl(3,nkl)= xyzb(3,nkl)-xyzijkl(3,4)
        enddo
      enddo
!
      nij= 0
      do iprim= 1,nprimijkl(1)
        exi= exijkl(iprim,1)
        ci = coijkl(iprim,1)
        xi = exi*xyzijkl(1,1)
        yi = exi*xyzijkl(2,1)
        zi = exi*xyzijkl(3,1)
        do jprim= 1,nprimijkl(2)
          exj= exijkl(jprim,2)
          cj = coijkl(jprim,2)
          exij= exi+exj
          exji= one/exij
          r12ex=exi*exj*exji*r12
          if(r12ex >= threshex) cycle
          nij= nij+1
          ex12(1,nij)= exij
          ex12(2,nij)= exji*half
          ex12(3,nij)= ci*cj*exp(-r12ex)*exji
          xyza(1,nij)=(xi+exj*xyzijkl(1,2))*exji
          xyza(2,nij)=(yi+exj*xyzijkl(2,2))*exji
          xyza(3,nij)=(zi+exj*xyzijkl(3,2))*exji
          xyzaj(1,nij)= xyza(1,nij)-xyzijkl(1,2)
          xyzaj(2,nij)= xyza(2,nij)-xyzijkl(2,2)
          xyzaj(3,nij)= xyza(3,nij)-xyzijkl(3,2)
        enddo
      enddo
!
      icount= 0
      do ij= 1,nij
        ex21h= ex12(2,ij)
        do kl= 1,nkl
          icount= icount+1
          xyzab(1)= xyza(1,ij)-xyzb(1,kl)
          xyzab(2)= xyza(2,ij)-xyzb(2,kl)
          xyzab(3)= xyza(3,ij)-xyzb(3,kl)
          ex41=one/(ex12(1,ij)+ex34(1,kl))
          ex41h=half*ex41
          ex43h=ex34(2,kl)
          b10t= ex12(1,ij)*ex43h*ex41
          bp01t=ex34(1,kl)*ex21h*ex41
          cx = ex12(1,ij)*xyzab(1)*ex41
          cy = ex12(1,ij)*xyzab(2)*ex41
          cz = ex12(1,ij)*xyzab(3)*ex41
          cpx= ex34(1,kl)*xyzab(1)*ex41
          cpy= ex34(1,kl)*xyzab(2)*ex41
          cpz= ex34(1,kl)*xyzab(3)*ex41
          dijkl= ex12(3,ij)*ex34(3,kl)*sqrt(ex41)
          tval=ex12(1,ij)*ex34(1,kl)*ex41*(xyzab(1)*xyzab(1)+xyzab(2)*xyzab(2)+xyzab(3)*xyzab(3))
!
          call rysquad(tval,trys,wrys,nroots)
          do iroot= 1,nroots
            b00 = ex41h*trys(iroot)
            b10 = ex43h-b10t*trys(iroot)
            bp01= ex21h-bp01t*trys(iroot)
            c00(1)= xyzbl(1,kl)+cx*trys(iroot)
            c00(2)= xyzbl(2,kl)+cy*trys(iroot)
            c00(3)= xyzbl(3,kl)+cz*trys(iroot)
            cp00(1) = xyzaj(1,ij)-cpx*trys(iroot)
            cp00(2) = xyzaj(2,ij)-cpy*trys(iroot)
            cp00(3) = xyzaj(3,ij)-cpz*trys(iroot)
! I(0,0)
            xyzint(1,0,0,0)= one
            xyzint(2,0,0,0)= one
            xyzint(3,0,0,0)= wrys(iroot)
!
            if(nmax >= 1) then
! I(1,0)
              xyzint(1,1,0,0)= c00(1)
              xyzint(2,1,0,0)= c00(2)
              xyzint(3,1,0,0)= c00(3)*wrys(iroot)
! I(n,0)
              do l= 2,nmax
                 xyzint(1,l,0,0)=(l-1)*b10*xyzint(1,l-2,0,0)+c00(1)*xyzint(1,l-1,0,0)
                 xyzint(2,l,0,0)=(l-1)*b10*xyzint(2,l-2,0,0)+c00(2)*xyzint(2,l-1,0,0)
                 xyzint(3,l,0,0)=(l-1)*b10*xyzint(3,l-2,0,0)+c00(3)*xyzint(3,l-1,0,0)
              enddo
            endif
!
            if(mmax >= 1) then
! I(0,1)
              xyzint(1,0,0,1)= cp00(1)
              xyzint(2,0,0,1)= cp00(2)
              xyzint(3,0,0,1)= cp00(3)*wrys(iroot)
! I(0,m)
              do j= 2,mmax
                xyzint(1,0,0,j)=(j-1)*bp01*xyzint(1,0,0,j-2)+cp00(1)*xyzint(1,0,0,j-1)
                xyzint(2,0,0,j)=(j-1)*bp01*xyzint(2,0,0,j-2)+cp00(2)*xyzint(2,0,0,j-1)
                xyzint(3,0,0,j)=(j-1)*bp01*xyzint(3,0,0,j-2)+cp00(3)*xyzint(3,0,0,j-1)
              enddo
              if(nmax >= 1) then
! I(1,1)
                xyzint(1,1,0,1)= b00+c00(1)*cp00(1)
                xyzint(2,1,0,1)= b00+c00(2)*cp00(2)
                xyzint(3,1,0,1)=(b00+c00(3)*cp00(3))*wrys(iroot)
! I(1,m)
                do j= 2,mmax
                  xyzint(1,1,0,j)=(j-1)*bp01*xyzint(1,1,0,j-2)+b00*xyzint(1,0,0,j-1) &
&                                   +cp00(1)*xyzint(1,1,0,j-1)
                  xyzint(2,1,0,j)=(j-1)*bp01*xyzint(2,1,0,j-2)+b00*xyzint(2,0,0,j-1) &
&                                   +cp00(2)*xyzint(2,1,0,j-1)
                  xyzint(3,1,0,j)=(j-1)*bp01*xyzint(3,1,0,j-2)+b00*xyzint(3,0,0,j-1) &
&                                   +cp00(3)*xyzint(3,1,0,j-1)
                enddo
              endif
            endif
! I(n,m)
            do j= 1,mmax
              do l= 2,nmax
                xyzint(1,l,0,j)=(l-1)*b10*xyzint(1,l-2,0,j)+j*b00*xyzint(1,l-1,0,j-1) &
&                                 +c00(1)*xyzint(1,l-1,0,j)
                xyzint(2,l,0,j)=(l-1)*b10*xyzint(2,l-2,0,j)+j*b00*xyzint(2,l-1,0,j-1) &
&                                 +c00(2)*xyzint(2,l-1,0,j)
                xyzint(3,l,0,j)=(l-1)*b10*xyzint(3,l-2,0,j)+j*b00*xyzint(3,l-1,0,j-1) &
&                                 +c00(3)*xyzint(3,l-1,0,j)
              enddo
            enddo
!
! I(l,k,m)
            do j= 0,mmax
              do k=1,nangijkl(3)
                do l=0,nmax-k
                  xyzint(1,l,k,j)= xyzint(1,l+1,k-1,j)-xyzkl(1)*xyzint(1,l,k-1,j)
                  xyzint(2,l,k,j)= xyzint(2,l+1,k-1,j)-xyzkl(2)*xyzint(2,l,k-1,j)
                  xyzint(3,l,k,j)= xyzint(3,l+1,k-1,j)-xyzkl(3)*xyzint(3,l,k-1,j)
                enddo
              enddo
            enddo
!
! I(l,k,j,i)
            do j= 0,mmax
              do k=0,nangijkl(3)
                do l=0,nangijkl(4)
                  rysint(1,icount,l,k,j,0,iroot)= xyzint(1,l,k,j)*dijkl
                  rysint(2,icount,l,k,j,0,iroot)= xyzint(2,l,k,j)
                  rysint(3,icount,l,k,j,0,iroot)= xyzint(3,l,k,j)
                enddo
              enddo
            enddo
            do i= 1,nangijkl(1)
              do j= 0,mmax-i
                do k=0,nangijkl(3)
                  do l=0,nangijkl(4)
                    rysint(1,icount,l,k,j,i,iroot)= rysint(1,icount,l,k,j+1,i-1,iroot) &
&                                                  -rysint(1,icount,l,k,j  ,i-1,iroot)*xyzij(1)
                    rysint(2,icount,l,k,j,i,iroot)= rysint(2,icount,l,k,j+1,i-1,iroot) &
&                                                  -rysint(2,icount,l,k,j  ,i-1,iroot)*xyzij(2) 
                    rysint(3,icount,l,k,j,i,iroot)= rysint(3,icount,l,k,j+1,i-1,iroot) &
&                                                  -rysint(3,icount,l,k,j  ,i-1,iroot)*xyzij(3) 
                  enddo
                enddo
              enddo
            enddo
          enddo
!
          if(icount == 8) then
            ii= 0
            do ix= nangijkl(1),0,-1
              do iy= nangijkl(1)-ix,0,-1
                iz= nangijkl(1)-ix-iy
                ii= ii+1
                jj= 0
                do jx= nangijkl(2),0,-1
                  do jy= nangijkl(2)-jx,0,-1
                    jz= nangijkl(2)-jx-jy
                    jj= jj+1
                    kk= 0
                    do kx= nangijkl(3),0,-1
                      do ky= nangijkl(3)-kx,0,-1
                        kz= nangijkl(3)-kx-ky
                        kk= kk+1
                        ll= 0
                        do lx= nangijkl(4),0,-1
                          do ly= nangijkl(4)-lx,0,-1
                            lz= nangijkl(4)-lx-ly
                            ll= ll+1
                            do iroot= 1,nroots
                              do nn= 1,8
                                eritmp(ll,kk,jj,ii)= eritmp(ll,kk,jj,ii) &
&                                                  +(rysint(1,nn,lx,kx,jx,ix,iroot) &
&                                                   *rysint(2,nn,ly,ky,jy,iy,iroot) &
&                                                   *rysint(3,nn,lz,kz,jz,iz,iroot))
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
            icount= 0
          endif
!
        enddo
      enddo
!
      if(icount /= 0) then
        ii= 0
        do ix= nangijkl(1),0,-1
          do iy= nangijkl(1)-ix,0,-1
            iz= nangijkl(1)-ix-iy
            ii= ii+1
            jj= 0
            do jx= nangijkl(2),0,-1
              do jy= nangijkl(2)-jx,0,-1
                jz= nangijkl(2)-jx-jy
                jj= jj+1
                kk= 0
                do kx= nangijkl(3),0,-1
                  do ky= nangijkl(3)-kx,0,-1
                    kz= nangijkl(3)-kx-ky
                    kk= kk+1
                    ll= 0
                    do lx= nangijkl(4),0,-1
                      do ly= nangijkl(4)-lx,0,-1
                        lz= nangijkl(4)-lx-ly
                        ll= ll+1
                        do iroot= 1,nroots
                          do nn= 1,icount
                            eritmp(ll,kk,jj,ii)= eritmp(ll,kk,jj,ii) &
&                                              +(rysint(1,nn,lx,kx,jx,ix,iroot) &
&                                               *rysint(2,nn,ly,ky,jy,iy,iroot) &
&                                               *rysint(3,nn,lz,kz,jz,iz,iroot))
                          enddo
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      endif
!
! Adjust coefficients for over d funcions.
!
      if(nbfijkl(1) == 5)then
        do j= 1,ncartijkl(2)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              do i= 1,6
                work(i)= eritmp(l,k,j,i)
              enddo
              eritmp(l,k,j,1)= work(2)*sqrt3
              eritmp(l,k,j,2)= work(5)*sqrt3
              eritmp(l,k,j,3)=(work(6)*two-work(1)-work(4))*half
              eritmp(l,k,j,4)= work(3)*sqrt3
              eritmp(l,k,j,5)=(work(1)-work(4))*sqrt3h
            enddo
          enddo
        enddo
      elseif(nbfijkl(1) == 6) then
        do j= 1,ncartijkl(2)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              eritmp(l,k,j,2)= eritmp(l,k,j,2)*sqrt3
              eritmp(l,k,j,3)= eritmp(l,k,j,3)*sqrt3
              eritmp(l,k,j,5)= eritmp(l,k,j,5)*sqrt3
            enddo
          enddo
        enddo
!
      elseif(nbfijkl(1) == 7) then
        do j= 1,ncartijkl(2)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              do i= 1,10
                work(i)= eritmp(l,k,j,i)
              enddo
              eritmp(l,k,j,1)=(-work(7)+three*work(2)                   )*facf1
              eritmp(l,k,j,2)=  work(5)                                  *facf2
              eritmp(l,k,j,3)=(-work(7)-work(2)+four*work(9)            )*facf3
              eritmp(l,k,j,4)=( two*work(10)-three*work(3)-three*work(8))*half
              eritmp(l,k,j,5)=(-work(1)-work(4)+four*work(6)            )*facf3
              eritmp(l,k,j,6)=( work(3)-work(8)                         )*facf4
              eritmp(l,k,j,7)=( work(1)-three*work(4)                   )*facf1
            enddo
          enddo
        enddo
      elseif(nbfijkl(1) == 10) then
        do j= 1,ncartijkl(2)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              eritmp(l,k,j, 2)= eritmp(l,k,j, 2)*sqrt5
              eritmp(l,k,j, 3)= eritmp(l,k,j, 3)*sqrt5
              eritmp(l,k,j, 4)= eritmp(l,k,j, 4)*sqrt5
              eritmp(l,k,j, 5)= eritmp(l,k,j, 5)*sqrt15
              eritmp(l,k,j, 6)= eritmp(l,k,j, 6)*sqrt5
              eritmp(l,k,j, 8)= eritmp(l,k,j, 8)*sqrt5
              eritmp(l,k,j, 9)= eritmp(l,k,j, 9)*sqrt5
            enddo
          enddo
        enddo
!
      elseif(nbfijkl(1) == 9) then
        do j= 1,ncartijkl(2)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              do i= 1,15
                work(i)= eritmp(l,k,j,i)
              enddo
              eritmp(l,k,j,1)=(work(2)-work(7))*facg1
              eritmp(l,k,j,2)=(-work(12)+work(5)*three)*facg2
              eritmp(l,k,j,3)=(-work(2)-work(7)+work(9)*six)*facg3
              eritmp(l,k,j,4)=(-work(12)*three+work(14)*four-work(5)*three)*facg4
              eritmp(l,k,j,5)=(work(1)*three+work(11)*three+work(15)*eight+work(4)*six &
&                              -work(6)*p24-work(13)*p24)*eighth
              eritmp(l,k,j,6)=(-work(3)*three+work(10)*four-work(8)*three)*facg4
              eritmp(l,k,j,7)=(-work(1)+work(11)+work(6)*six-work(13)*six)*facg5
              eritmp(l,k,j,8)=(work(3)-work(8)*three)*facg2
              eritmp(l,k,j,9)=(work(1)+work(11)-work(4)*six)*facg6
            enddo
          enddo
        enddo
      elseif(nbfijkl(1) == 15) then
        do j= 1,ncartijkl(2)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              eritmp(l,k,j, 2)= eritmp(l,k,j, 2)*sqrt7
              eritmp(l,k,j, 3)= eritmp(l,k,j, 3)*sqrt7
              eritmp(l,k,j, 4)= eritmp(l,k,j, 4)*sqrt35third
              eritmp(l,k,j, 5)= eritmp(l,k,j, 5)*sqrt35
              eritmp(l,k,j, 6)= eritmp(l,k,j, 6)*sqrt35third
              eritmp(l,k,j, 7)= eritmp(l,k,j, 7)*sqrt7
              eritmp(l,k,j, 8)= eritmp(l,k,j, 8)*sqrt35
              eritmp(l,k,j, 9)= eritmp(l,k,j, 9)*sqrt35
              eritmp(l,k,j,10)= eritmp(l,k,j,10)*sqrt7
              eritmp(l,k,j,12)= eritmp(l,k,j,12)*sqrt7
              eritmp(l,k,j,13)= eritmp(l,k,j,13)*sqrt35third
              eritmp(l,k,j,14)= eritmp(l,k,j,14)*sqrt7
            enddo
          enddo
        enddo
!
      elseif(nbfijkl(1) == 11) then
        do j= 1,ncartijkl(2)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              do i= 1,21
                work(i)= eritmp(l,k,j,i)
              enddo
              eritmp(l,k,j, 1)=(work(2)*five-work(7)*ten+work(16))*fach1
              eritmp(l,k,j, 2)=(work(5)*four-work(12)*four)*fach2
              eritmp(l,k,j, 3)=(-work(2)*three-work(7)*two+work(9)*p24+work(16) &
&                              -work(18)*eight)*fach3
              eritmp(l,k,j, 4)=(-work(5)*two-work(12)*two+work(14)*four)*fach4
              eritmp(l,k,j, 5)=(work(2)+work(7)*two-work(9)*twelve+work(16)-work(18)*twelve &
&                              +work(20)*eight)*fach5
              eritmp(l,k,j, 6)=(work(3)*p15+work(8)*p30-work(10)*p40+work(17)*p15-work(19)*p40 &
&                              +work(21)*eight)*eighth
              eritmp(l,k,j, 7)=(work(1)+work(4)*two-work(6)*twelve+work(11)-work(13)*twelve &
&                              +work(15)*eight)*fach5
              eritmp(l,k,j, 8)=(-work(3)+work(10)*two+work(17)-work(19)*two)*fach4
              eritmp(l,k,j, 9)=(-work(1)+work(4)*two+work(6)*eight+work(11)*three &
&                              -work(13)*p24)*fach3
              eritmp(l,k,j,10)=(work(3)-work(8)*six+work(17))*fach2
              eritmp(l,k,j,11)=(work(1)-work(4)*ten+work(11)*five)*fach1
            enddo
          enddo
        enddo
      elseif(nbfijkl(1) == 21) then
        do j= 1,ncartijkl(2)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              eritmp(l,k,j, 2)= eritmp(l,k,j, 2)*three
              eritmp(l,k,j, 3)= eritmp(l,k,j, 3)*three
              eritmp(l,k,j, 4)= eritmp(l,k,j, 4)*sqrt21
              eritmp(l,k,j, 5)= eritmp(l,k,j, 5)*sqrt63
              eritmp(l,k,j, 6)= eritmp(l,k,j, 6)*sqrt21
              eritmp(l,k,j, 7)= eritmp(l,k,j, 7)*sqrt21
              eritmp(l,k,j, 8)= eritmp(l,k,j, 8)*sqrt105
              eritmp(l,k,j, 9)= eritmp(l,k,j, 9)*sqrt105
              eritmp(l,k,j,10)= eritmp(l,k,j,10)*sqrt21
              eritmp(l,k,j,11)= eritmp(l,k,j,11)*three
              eritmp(l,k,j,12)= eritmp(l,k,j,12)*sqrt63
              eritmp(l,k,j,13)= eritmp(l,k,j,13)*sqrt105
              eritmp(l,k,j,14)= eritmp(l,k,j,14)*sqrt63
              eritmp(l,k,j,15)= eritmp(l,k,j,15)*three
              eritmp(l,k,j,17)= eritmp(l,k,j,17)*three
              eritmp(l,k,j,18)= eritmp(l,k,j,18)*sqrt21
              eritmp(l,k,j,19)= eritmp(l,k,j,19)*sqrt21
              eritmp(l,k,j,20)= eritmp(l,k,j,20)*three
            enddo
          enddo
        enddo
!
      elseif(nbfijkl(1) == 13) then
        do j= 1,ncartijkl(2)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              do i= 1,28
                work(i)= eritmp(l,k,j,i)
              enddo
              eritmp(l,k,j, 1)=(work(2)*six-work(7)*p20+work(16)*six)*faci1
              eritmp(l,k,j, 2)=(work(5)*five-work(12)*ten+work(23))*faci2
              eritmp(l,k,j, 3)=(-work(2)*four+work(9)*p40+work(16)*four-work(18)*p40)*faci3
              eritmp(l,k,j, 4)=(-work(5)*p9-work(12)*six+work(14)*p24+work(23)*three &
&                              -work(25)*eight)*faci4
              eritmp(l,k,j, 5)=(work(2)*two+work(7)*four-work(9)*p32+work(16)*two &
&                              -work(18)*p32+work(20)*p32)*faci5
              eritmp(l,k,j, 6)=(work(5)*five+work(12)*ten-work(14)*p20+work(23)*five &
&                              -work(25)*p20+work(27)*eight)*faci6
              eritmp(l,k,j, 7)=(-work(1)*five-work(4)*p15+work(6)*p90-work(11)*p15 &
&                              +work(13)*p180-work(15)*p120-work(22)*five+work(24)*p90 &
&                              -work(26)*p120+work(28)*p16)*sixteenth
              eritmp(l,k,j, 8)=(work(3)*five+work(8)*ten-work(10)*p20+work(17)*five &
&                              -work(19)*p20+work(21)*eight)*faci6
              eritmp(l,k,j, 9)=(work(1)+work(4)-work(6)*p16-work(11)+work(15)*p16 &
&                              -work(22)+work(24)*p16-work(26)*p16)*faci5
              eritmp(l,k,j,10)=(-work(3)*three+work(8)*six+work(10)*eight+work(17)*p9 &
&                              -work(19)*p24)*faci4
              eritmp(l,k,j,11)=(-work(1)+work(4)*five+work(6)*ten+work(11)*five &
&                              -work(13)*p60-work(22)+work(24)*ten)*faci3
              eritmp(l,k,j,12)=(work(3)-work(8)*ten+work(17)*five)*faci2
              eritmp(l,k,j,13)=(work(1)-work(4)*p15+work(11)*p15-work(22))*faci1
            enddo
          enddo
        enddo
      elseif(nbfijkl(1) == 28) then
        do j= 1,ncartijkl(2)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              eritmp(l,k,j, 2)= eritmp(l,k,j, 2)*sqrt11
              eritmp(l,k,j, 3)= eritmp(l,k,j, 3)*sqrt11
              eritmp(l,k,j, 4)= eritmp(l,k,j, 4)*sqrt33
              eritmp(l,k,j, 5)= eritmp(l,k,j, 5)*sqrt99
              eritmp(l,k,j, 6)= eritmp(l,k,j, 6)*sqrt33
              eritmp(l,k,j, 7)= eritmp(l,k,j, 7)*sqrt231fifth
              eritmp(l,k,j, 8)= eritmp(l,k,j, 8)*sqrt231
              eritmp(l,k,j, 9)= eritmp(l,k,j, 9)*sqrt231
              eritmp(l,k,j,10)= eritmp(l,k,j,10)*sqrt231fifth
              eritmp(l,k,j,11)= eritmp(l,k,j,11)*sqrt33
              eritmp(l,k,j,12)= eritmp(l,k,j,12)*sqrt231
              eritmp(l,k,j,13)= eritmp(l,k,j,13)*sqrt385
              eritmp(l,k,j,14)= eritmp(l,k,j,14)*sqrt231
              eritmp(l,k,j,15)= eritmp(l,k,j,15)*sqrt33
              eritmp(l,k,j,16)= eritmp(l,k,j,16)*sqrt11
              eritmp(l,k,j,17)= eritmp(l,k,j,17)*sqrt99
              eritmp(l,k,j,18)= eritmp(l,k,j,18)*sqrt231
              eritmp(l,k,j,19)= eritmp(l,k,j,19)*sqrt231
              eritmp(l,k,j,20)= eritmp(l,k,j,20)*sqrt99
              eritmp(l,k,j,21)= eritmp(l,k,j,21)*sqrt11
              eritmp(l,k,j,23)= eritmp(l,k,j,23)*sqrt11
              eritmp(l,k,j,24)= eritmp(l,k,j,24)*sqrt33
              eritmp(l,k,j,25)= eritmp(l,k,j,25)*sqrt231fifth
              eritmp(l,k,j,26)= eritmp(l,k,j,26)*sqrt33
              eritmp(l,k,j,27)= eritmp(l,k,j,27)*sqrt11
            enddo
          enddo
        enddo
      endif
!
      if(nbfijkl(2) == 5) then
        do i= 1,nbfijkl(1)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              do j= 1,6
                work(j)= eritmp(l,k,j,i)
              enddo
              eritmp(l,k,1,i)= work(2)*sqrt3
              eritmp(l,k,2,i)= work(5)*sqrt3
              eritmp(l,k,3,i)=(work(6)*two-work(1)-work(4))*half
              eritmp(l,k,4,i)= work(3)*sqrt3
              eritmp(l,k,5,i)=(work(1)-work(4))*sqrt3h
            enddo
          enddo
        enddo
      elseif(nbfijkl(2) == 6) then
        do i= 1,nbfijkl(1)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              eritmp(l,k,2,i)= eritmp(l,k,2,i)*sqrt3
              eritmp(l,k,3,i)= eritmp(l,k,3,i)*sqrt3
              eritmp(l,k,5,i)= eritmp(l,k,5,i)*sqrt3
            enddo
          enddo
        enddo
!
      elseif(nbfijkl(2) == 7) then
        do i= 1,nbfijkl(1)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              do j= 1,10
                work(j)= eritmp(l,k,j,i)
              enddo
              eritmp(l,k,1,i)=(-work(7)+three*work(2)                   )*facf1
              eritmp(l,k,2,i)=  work(5)                                  *facf2
              eritmp(l,k,3,i)=(-work(7)-work(2)+four*work(9)            )*facf3
              eritmp(l,k,4,i)=( two*work(10)-three*work(3)-three*work(8))*half
              eritmp(l,k,5,i)=(-work(1)-work(4)+four*work(6)            )*facf3
              eritmp(l,k,6,i)=( work(3)-work(8)                         )*facf4
              eritmp(l,k,7,i)=( work(1)-three*work(4)                   )*facf1
            enddo
          enddo
        enddo
      elseif(nbfijkl(2) == 10) then
        do i= 1,nbfijkl(1)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              eritmp(l,k, 2,i)= eritmp(l,k, 2,i)*sqrt5
              eritmp(l,k, 3,i)= eritmp(l,k, 3,i)*sqrt5
              eritmp(l,k, 4,i)= eritmp(l,k, 4,i)*sqrt5
              eritmp(l,k, 5,i)= eritmp(l,k, 5,i)*sqrt15
              eritmp(l,k, 6,i)= eritmp(l,k, 6,i)*sqrt5
              eritmp(l,k, 8,i)= eritmp(l,k, 8,i)*sqrt5
              eritmp(l,k, 9,i)= eritmp(l,k, 9,i)*sqrt5
            enddo
          enddo
        enddo
!
      elseif(nbfijkl(2) == 9) then
        do i= 1,nbfijkl(1)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              do j= 1,15
                work(j)= eritmp(l,k,j,i)
              enddo
              eritmp(l,k,1,i)=(work(2)-work(7))*facg1
              eritmp(l,k,2,i)=(-work(12)+work(5)*three)*facg2
              eritmp(l,k,3,i)=(-work(2)-work(7)+work(9)*six)*facg3
              eritmp(l,k,4,i)=(-work(12)*three+work(14)*four-work(5)*three)*facg4
              eritmp(l,k,5,i)=(work(1)*three+work(11)*three+work(15)*eight+work(4)*six &
&                              -work(6)*p24-work(13)*p24)*eighth
              eritmp(l,k,6,i)=(-work(3)*three+work(10)*four-work(8)*three)*facg4
              eritmp(l,k,7,i)=(-work(1)+work(11)+work(6)*six-work(13)*six)*facg5
              eritmp(l,k,8,i)=(work(3)-work(8)*three)*facg2
              eritmp(l,k,9,i)=(work(1)+work(11)-work(4)*six)*facg6
            enddo
          enddo
        enddo
      elseif(nbfijkl(2) == 15) then
        do i= 1,nbfijkl(1)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              eritmp(l,k, 2,i)= eritmp(l,k, 2,i)*sqrt7
              eritmp(l,k, 3,i)= eritmp(l,k, 3,i)*sqrt7
              eritmp(l,k, 4,i)= eritmp(l,k, 4,i)*sqrt35third
              eritmp(l,k, 5,i)= eritmp(l,k, 5,i)*sqrt35
              eritmp(l,k, 6,i)= eritmp(l,k, 6,i)*sqrt35third
              eritmp(l,k, 7,i)= eritmp(l,k, 7,i)*sqrt7
              eritmp(l,k, 8,i)= eritmp(l,k, 8,i)*sqrt35
              eritmp(l,k, 9,i)= eritmp(l,k, 9,i)*sqrt35
              eritmp(l,k,10,i)= eritmp(l,k,10,i)*sqrt7
              eritmp(l,k,12,i)= eritmp(l,k,12,i)*sqrt7
              eritmp(l,k,13,i)= eritmp(l,k,13,i)*sqrt35third
              eritmp(l,k,14,i)= eritmp(l,k,14,i)*sqrt7
            enddo
          enddo
        enddo
!
      elseif(nbfijkl(2) == 11) then
        do i= 1,nbfijkl(1)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              do j= 1,21
                work(j)= eritmp(l,k,j,i)
              enddo
              eritmp(l,k, 1,i)=(work(2)*five-work(7)*ten+work(16))*fach1
              eritmp(l,k, 2,i)=(work(5)*four-work(12)*four)*fach2
              eritmp(l,k, 3,i)=(-work(2)*three-work(7)*two+work(9)*p24+work(16) &
&                              -work(18)*eight)*fach3
              eritmp(l,k, 4,i)=(-work(5)*two-work(12)*two+work(14)*four)*fach4
              eritmp(l,k, 5,i)=(work(2)+work(7)*two-work(9)*twelve+work(16)-work(18)*twelve &
&                              +work(20)*eight)*fach5
              eritmp(l,k, 6,i)=(work(3)*p15+work(8)*p30-work(10)*p40+work(17)*p15-work(19)*p40 &
&                              +work(21)*eight)*eighth
              eritmp(l,k, 7,i)=(work(1)+work(4)*two-work(6)*twelve+work(11)-work(13)*twelve &
&                              +work(15)*eight)*fach5
              eritmp(l,k, 8,i)=(-work(3)+work(10)*two+work(17)-work(19)*two)*fach4
              eritmp(l,k, 9,i)=(-work(1)+work(4)*two+work(6)*eight+work(11)*three &
&                              -work(13)*p24)*fach3
              eritmp(l,k,10,i)=(work(3)-work(8)*six+work(17))*fach2
              eritmp(l,k,11,i)=(work(1)-work(4)*ten+work(11)*five)*fach1
            enddo
          enddo
        enddo
      elseif(nbfijkl(2) == 21) then
        do i= 1,nbfijkl(1)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              eritmp(l,k, 2,i)= eritmp(l,k, 2,i)*three
              eritmp(l,k, 3,i)= eritmp(l,k, 3,i)*three
              eritmp(l,k, 4,i)= eritmp(l,k, 4,i)*sqrt21
              eritmp(l,k, 5,i)= eritmp(l,k, 5,i)*sqrt63
              eritmp(l,k, 6,i)= eritmp(l,k, 6,i)*sqrt21
              eritmp(l,k, 7,i)= eritmp(l,k, 7,i)*sqrt21
              eritmp(l,k, 8,i)= eritmp(l,k, 8,i)*sqrt105
              eritmp(l,k, 9,i)= eritmp(l,k, 9,i)*sqrt105
              eritmp(l,k,10,i)= eritmp(l,k,10,i)*sqrt21
              eritmp(l,k,11,i)= eritmp(l,k,11,i)*three
              eritmp(l,k,12,i)= eritmp(l,k,12,i)*sqrt63
              eritmp(l,k,13,i)= eritmp(l,k,13,i)*sqrt105
              eritmp(l,k,14,i)= eritmp(l,k,14,i)*sqrt63
              eritmp(l,k,15,i)= eritmp(l,k,15,i)*three
              eritmp(l,k,17,i)= eritmp(l,k,17,i)*three
              eritmp(l,k,18,i)= eritmp(l,k,18,i)*sqrt21
              eritmp(l,k,19,i)= eritmp(l,k,19,i)*sqrt21
              eritmp(l,k,20,i)= eritmp(l,k,20,i)*three
            enddo
          enddo
        enddo
!
      elseif(nbfijkl(2) == 13) then
        do i= 1,nbfijkl(1)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              do j= 1,28
                work(j)= eritmp(l,k,j,i)
              enddo
              eritmp(l,k, 1,i)=(work(2)*six-work(7)*p20+work(16)*six)*faci1
              eritmp(l,k, 2,i)=(work(5)*five-work(12)*ten+work(23))*faci2
              eritmp(l,k, 3,i)=(-work(2)*four+work(9)*p40+work(16)*four-work(18)*p40)*faci3
              eritmp(l,k, 4,i)=(-work(5)*p9-work(12)*six+work(14)*p24+work(23)*three &
&                              -work(25)*eight)*faci4
              eritmp(l,k, 5,i)=(work(2)*two+work(7)*four-work(9)*p32+work(16)*two &
&                              -work(18)*p32+work(20)*p32)*faci5
              eritmp(l,k, 6,i)=(work(5)*five+work(12)*ten-work(14)*p20+work(23)*five &
&                              -work(25)*p20+work(27)*eight)*faci6
              eritmp(l,k, 7,i)=(-work(1)*five-work(4)*p15+work(6)*p90-work(11)*p15 &
&                              +work(13)*p180-work(15)*p120-work(22)*five+work(24)*p90 &
&                              -work(26)*p120+work(28)*p16)*sixteenth
              eritmp(l,k, 8,i)=(work(3)*five+work(8)*ten-work(10)*p20+work(17)*five &
&                              -work(19)*p20+work(21)*eight)*faci6
              eritmp(l,k, 9,i)=(work(1)+work(4)-work(6)*p16-work(11)+work(15)*p16 &
&                              -work(22)+work(24)*p16-work(26)*p16)*faci5
              eritmp(l,k,10,i)=(-work(3)*three+work(8)*six+work(10)*eight+work(17)*p9 &
&                              -work(19)*p24)*faci4
              eritmp(l,k,11,i)=(-work(1)+work(4)*five+work(6)*ten+work(11)*five &
&                              -work(13)*p60-work(22)+work(24)*ten)*faci3
              eritmp(l,k,12,i)=(work(3)-work(8)*ten+work(17)*five)*faci2
              eritmp(l,k,13,i)=(work(1)-work(4)*p15+work(11)*p15-work(22))*faci1
            enddo
          enddo
        enddo
      elseif(nbfijkl(2) == 28) then
        do i= 1,nbfijkl(1)
          do k= 1,ncartijkl(3)
            do l= 1,ncartijkl(4)
              eritmp(l,k, 2,i)= eritmp(l,k, 2,i)*sqrt11
              eritmp(l,k, 3,i)= eritmp(l,k, 3,i)*sqrt11
              eritmp(l,k, 4,i)= eritmp(l,k, 4,i)*sqrt33
              eritmp(l,k, 5,i)= eritmp(l,k, 5,i)*sqrt99
              eritmp(l,k, 6,i)= eritmp(l,k, 6,i)*sqrt33
              eritmp(l,k, 7,i)= eritmp(l,k, 7,i)*sqrt231fifth
              eritmp(l,k, 8,i)= eritmp(l,k, 8,i)*sqrt231
              eritmp(l,k, 9,i)= eritmp(l,k, 9,i)*sqrt231
              eritmp(l,k,10,i)= eritmp(l,k,10,i)*sqrt231fifth
              eritmp(l,k,11,i)= eritmp(l,k,11,i)*sqrt33
              eritmp(l,k,12,i)= eritmp(l,k,12,i)*sqrt231
              eritmp(l,k,13,i)= eritmp(l,k,13,i)*sqrt385
              eritmp(l,k,14,i)= eritmp(l,k,14,i)*sqrt231
              eritmp(l,k,15,i)= eritmp(l,k,15,i)*sqrt33
              eritmp(l,k,16,i)= eritmp(l,k,16,i)*sqrt11
              eritmp(l,k,17,i)= eritmp(l,k,17,i)*sqrt99
              eritmp(l,k,18,i)= eritmp(l,k,18,i)*sqrt231
              eritmp(l,k,19,i)= eritmp(l,k,19,i)*sqrt231
              eritmp(l,k,20,i)= eritmp(l,k,20,i)*sqrt99
              eritmp(l,k,21,i)= eritmp(l,k,21,i)*sqrt11
              eritmp(l,k,23,i)= eritmp(l,k,23,i)*sqrt11
              eritmp(l,k,24,i)= eritmp(l,k,24,i)*sqrt33
              eritmp(l,k,25,i)= eritmp(l,k,25,i)*sqrt231fifth
              eritmp(l,k,26,i)= eritmp(l,k,26,i)*sqrt33
              eritmp(l,k,27,i)= eritmp(l,k,27,i)*sqrt11
            enddo
          enddo
        enddo
      endif
!
      if(nbfijkl(3) == 5)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do l= 1,ncartijkl(4)
              do k= 1,6
                work(k)= eritmp(l,k,j,i)
              enddo
              eritmp(l,1,j,i)= work(2)*sqrt3
              eritmp(l,2,j,i)= work(5)*sqrt3
              eritmp(l,3,j,i)=(work(6)*two-work(1)-work(4))*half
              eritmp(l,4,j,i)= work(3)*sqrt3
              eritmp(l,5,j,i)=(work(1)-work(4))*sqrt3h
            enddo
          enddo
        enddo
      elseif(nbfijkl(3) == 6)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do l= 1,ncartijkl(4)
              eritmp(l,2,j,i)= eritmp(l,2,j,i)*sqrt3
              eritmp(l,3,j,i)= eritmp(l,3,j,i)*sqrt3
              eritmp(l,5,j,i)= eritmp(l,5,j,i)*sqrt3
            enddo
          enddo
        enddo
!
      elseif(nbfijkl(3) == 7)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do l= 1,ncartijkl(4)
              do k= 1,10
                work(k)= eritmp(l,k,j,i)
              enddo
              eritmp(l,1,j,i)=(-work(7)+three*work(2)                   )*facf1
              eritmp(l,2,j,i)=  work(5)                                  *facf2
              eritmp(l,3,j,i)=(-work(7)-work(2)+four*work(9)            )*facf3
              eritmp(l,4,j,i)=( two*work(10)-three*work(3)-three*work(8))*half
              eritmp(l,5,j,i)=(-work(1)-work(4)+four*work(6)            )*facf3
              eritmp(l,6,j,i)=( work(3)-work(8)                         )*facf4
              eritmp(l,7,j,i)=( work(1)-three*work(4)                   )*facf1
            enddo
          enddo
        enddo
      elseif(nbfijkl(3) == 10)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do l= 1,ncartijkl(4)
              eritmp(l, 2,j,i)= eritmp(l, 2,j,i)*sqrt5
              eritmp(l, 3,j,i)= eritmp(l, 3,j,i)*sqrt5
              eritmp(l, 4,j,i)= eritmp(l, 4,j,i)*sqrt5
              eritmp(l, 5,j,i)= eritmp(l, 5,j,i)*sqrt15
              eritmp(l, 6,j,i)= eritmp(l, 6,j,i)*sqrt5
              eritmp(l, 8,j,i)= eritmp(l, 8,j,i)*sqrt5
              eritmp(l, 9,j,i)= eritmp(l, 9,j,i)*sqrt5
            enddo
          enddo
        enddo
!
      elseif(nbfijkl(3) == 9)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do l= 1,ncartijkl(4)
              do k= 1,15
                work(k)= eritmp(l,k,j,i)
              enddo
              eritmp(l,1,j,i)=(work(2)-work(7))*facg1
              eritmp(l,2,j,i)=(-work(12)+work(5)*three)*facg2
              eritmp(l,3,j,i)=(-work(2)-work(7)+work(9)*six)*facg3
              eritmp(l,4,j,i)=(-work(12)*three+work(14)*four-work(5)*three)*facg4
              eritmp(l,5,j,i)=(work(1)*three+work(11)*three+work(15)*eight+work(4)*six &
&                              -work(6)*p24-work(13)*p24)*eighth
              eritmp(l,6,j,i)=(-work(3)*three+work(10)*four-work(8)*three)*facg4
              eritmp(l,7,j,i)=(-work(1)+work(11)+work(6)*six-work(13)*six)*facg5
              eritmp(l,8,j,i)=(work(3)-work(8)*three)*facg2
              eritmp(l,9,j,i)=(work(1)+work(11)-work(4)*six)*facg6
            enddo
          enddo
        enddo
      elseif(nbfijkl(3) == 15)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do l= 1,ncartijkl(4)
              eritmp(l, 2,j,i)= eritmp(l, 2,j,i)*sqrt7
              eritmp(l, 3,j,i)= eritmp(l, 3,j,i)*sqrt7
              eritmp(l, 4,j,i)= eritmp(l, 4,j,i)*sqrt35third
              eritmp(l, 5,j,i)= eritmp(l, 5,j,i)*sqrt35
              eritmp(l, 6,j,i)= eritmp(l, 6,j,i)*sqrt35third
              eritmp(l, 7,j,i)= eritmp(l, 7,j,i)*sqrt7
              eritmp(l, 8,j,i)= eritmp(l, 8,j,i)*sqrt35
              eritmp(l, 9,j,i)= eritmp(l, 9,j,i)*sqrt35
              eritmp(l,10,j,i)= eritmp(l,10,j,i)*sqrt7
              eritmp(l,12,j,i)= eritmp(l,12,j,i)*sqrt7
              eritmp(l,13,j,i)= eritmp(l,13,j,i)*sqrt35third
              eritmp(l,14,j,i)= eritmp(l,14,j,i)*sqrt7
            enddo
          enddo
        enddo
!
      elseif(nbfijkl(3) == 11)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do l= 1,ncartijkl(4)
              do k= 1,21
                work(k)= eritmp(l,k,j,i)
              enddo
              eritmp(l, 1,j,i)=(work(2)*five-work(7)*ten+work(16))*fach1
              eritmp(l, 2,j,i)=(work(5)*four-work(12)*four)*fach2
              eritmp(l, 3,j,i)=(-work(2)*three-work(7)*two+work(9)*p24+work(16) &
&                              -work(18)*eight)*fach3
              eritmp(l, 4,j,i)=(-work(5)*two-work(12)*two+work(14)*four)*fach4
              eritmp(l, 5,j,i)=(work(2)+work(7)*two-work(9)*twelve+work(16)-work(18)*twelve &
&                              +work(20)*eight)*fach5
              eritmp(l, 6,j,i)=(work(3)*p15+work(8)*p30-work(10)*p40+work(17)*p15-work(19)*p40 &
&                              +work(21)*eight)*eighth
              eritmp(l, 7,j,i)=(work(1)+work(4)*two-work(6)*twelve+work(11)-work(13)*twelve &
&                              +work(15)*eight)*fach5
              eritmp(l, 8,j,i)=(-work(3)+work(10)*two+work(17)-work(19)*two)*fach4
              eritmp(l, 9,j,i)=(-work(1)+work(4)*two+work(6)*eight+work(11)*three &
&                              -work(13)*p24)*fach3
              eritmp(l,10,j,i)=(work(3)-work(8)*six+work(17))*fach2
              eritmp(l,11,j,i)=(work(1)-work(4)*ten+work(11)*five)*fach1
            enddo
          enddo
        enddo
      elseif(nbfijkl(3) == 21)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do l= 1,ncartijkl(4)
              eritmp(l, 2,j,i)= eritmp(l, 2,j,i)*three
              eritmp(l, 3,j,i)= eritmp(l, 3,j,i)*three
              eritmp(l, 4,j,i)= eritmp(l, 4,j,i)*sqrt21
              eritmp(l, 5,j,i)= eritmp(l, 5,j,i)*sqrt63
              eritmp(l, 6,j,i)= eritmp(l, 6,j,i)*sqrt21
              eritmp(l, 7,j,i)= eritmp(l, 7,j,i)*sqrt21
              eritmp(l, 8,j,i)= eritmp(l, 8,j,i)*sqrt105
              eritmp(l, 9,j,i)= eritmp(l, 9,j,i)*sqrt105
              eritmp(l,10,j,i)= eritmp(l,10,j,i)*sqrt21
              eritmp(l,11,j,i)= eritmp(l,11,j,i)*three
              eritmp(l,12,j,i)= eritmp(l,12,j,i)*sqrt63
              eritmp(l,13,j,i)= eritmp(l,13,j,i)*sqrt105
              eritmp(l,14,j,i)= eritmp(l,14,j,i)*sqrt63
              eritmp(l,15,j,i)= eritmp(l,15,j,i)*three
              eritmp(l,17,j,i)= eritmp(l,17,j,i)*three
              eritmp(l,18,j,i)= eritmp(l,18,j,i)*sqrt21
              eritmp(l,19,j,i)= eritmp(l,19,j,i)*sqrt21
              eritmp(l,20,j,i)= eritmp(l,20,j,i)*three
            enddo
          enddo
        enddo
!
      elseif(nbfijkl(3) == 13)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do l= 1,ncartijkl(4)
              do k= 1,28
                work(k)= eritmp(l,k,j,i)
              enddo
              eritmp(l, 1,j,i)=(work(2)*six-work(7)*p20+work(16)*six)*faci1
              eritmp(l, 2,j,i)=(work(5)*five-work(12)*ten+work(23))*faci2
              eritmp(l, 3,j,i)=(-work(2)*four+work(9)*p40+work(16)*four-work(18)*p40)*faci3
              eritmp(l, 4,j,i)=(-work(5)*p9-work(12)*six+work(14)*p24+work(23)*three &
&                              -work(25)*eight)*faci4
              eritmp(l, 5,j,i)=(work(2)*two+work(7)*four-work(9)*p32+work(16)*two &
&                              -work(18)*p32+work(20)*p32)*faci5
              eritmp(l, 6,j,i)=(work(5)*five+work(12)*ten-work(14)*p20+work(23)*five &
&                              -work(25)*p20+work(27)*eight)*faci6
              eritmp(l, 7,j,i)=(-work(1)*five-work(4)*p15+work(6)*p90-work(11)*p15 &
&                              +work(13)*p180-work(15)*p120-work(22)*five+work(24)*p90 &
&                              -work(26)*p120+work(28)*p16)*sixteenth
              eritmp(l, 8,j,i)=(work(3)*five+work(8)*ten-work(10)*p20+work(17)*five &
&                              -work(19)*p20+work(21)*eight)*faci6
              eritmp(l, 9,j,i)=(work(1)+work(4)-work(6)*p16-work(11)+work(15)*p16 &
&                              -work(22)+work(24)*p16-work(26)*p16)*faci5
              eritmp(l,10,j,i)=(-work(3)*three+work(8)*six+work(10)*eight+work(17)*p9 &
&                              -work(19)*p24)*faci4
              eritmp(l,11,j,i)=(-work(1)+work(4)*five+work(6)*ten+work(11)*five &
&                              -work(13)*p60-work(22)+work(24)*ten)*faci3
              eritmp(l,12,j,i)=(work(3)-work(8)*ten+work(17)*five)*faci2
              eritmp(l,13,j,i)=(work(1)-work(4)*p15+work(11)*p15-work(22))*faci1
            enddo
          enddo
        enddo
      elseif(nbfijkl(3) == 28)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do l= 1,ncartijkl(4)
              eritmp(l, 2,j,i)= eritmp(l, 2,j,i)*sqrt11
              eritmp(l, 3,j,i)= eritmp(l, 3,j,i)*sqrt11
              eritmp(l, 4,j,i)= eritmp(l, 4,j,i)*sqrt33
              eritmp(l, 5,j,i)= eritmp(l, 5,j,i)*sqrt99
              eritmp(l, 6,j,i)= eritmp(l, 6,j,i)*sqrt33
              eritmp(l, 7,j,i)= eritmp(l, 7,j,i)*sqrt231fifth
              eritmp(l, 8,j,i)= eritmp(l, 8,j,i)*sqrt231
              eritmp(l, 9,j,i)= eritmp(l, 9,j,i)*sqrt231
              eritmp(l,10,j,i)= eritmp(l,10,j,i)*sqrt231fifth
              eritmp(l,11,j,i)= eritmp(l,11,j,i)*sqrt33
              eritmp(l,12,j,i)= eritmp(l,12,j,i)*sqrt231
              eritmp(l,13,j,i)= eritmp(l,13,j,i)*sqrt385
              eritmp(l,14,j,i)= eritmp(l,14,j,i)*sqrt231
              eritmp(l,15,j,i)= eritmp(l,15,j,i)*sqrt33
              eritmp(l,16,j,i)= eritmp(l,16,j,i)*sqrt11
              eritmp(l,17,j,i)= eritmp(l,17,j,i)*sqrt99
              eritmp(l,18,j,i)= eritmp(l,18,j,i)*sqrt231
              eritmp(l,19,j,i)= eritmp(l,19,j,i)*sqrt231
              eritmp(l,20,j,i)= eritmp(l,20,j,i)*sqrt99
              eritmp(l,21,j,i)= eritmp(l,21,j,i)*sqrt11
              eritmp(l,23,j,i)= eritmp(l,23,j,i)*sqrt11
              eritmp(l,24,j,i)= eritmp(l,24,j,i)*sqrt33
              eritmp(l,25,j,i)= eritmp(l,25,j,i)*sqrt231fifth
              eritmp(l,26,j,i)= eritmp(l,26,j,i)*sqrt33
              eritmp(l,27,j,i)= eritmp(l,27,j,i)*sqrt11
            enddo
          enddo
        enddo
      endif
!
      if(nbfijkl(4) <= 3)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do k= 1,nbfijkl(3)
              do l= 1,nbfijkl(4)
                twoeri(l,k,j,i)= eritmp(l,k,j,i)
              enddo
            enddo
          enddo
        enddo
      elseif(nbfijkl(4) == 5)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do k= 1,nbfijkl(3)
              do l= 1,6
                work(l)= eritmp(l,k,j,i)
              enddo
              twoeri(1,k,j,i)= work(2)*sqrt3
              twoeri(2,k,j,i)= work(5)*sqrt3
              twoeri(3,k,j,i)=(work(6)*two-work(1)-work(4))*half
              twoeri(4,k,j,i)= work(3)*sqrt3
              twoeri(5,k,j,i)=(work(1)-work(4))*sqrt3h
            enddo
          enddo
        enddo
      elseif(nbfijkl(4) == 6)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do k= 1,nbfijkl(3)
              twoeri(1,k,j,i)= eritmp(1,k,j,i)
              twoeri(2,k,j,i)= eritmp(2,k,j,i)*sqrt3
              twoeri(3,k,j,i)= eritmp(3,k,j,i)*sqrt3
              twoeri(4,k,j,i)= eritmp(4,k,j,i)
              twoeri(5,k,j,i)= eritmp(5,k,j,i)*sqrt3
              twoeri(6,k,j,i)= eritmp(6,k,j,i)
            enddo
          enddo
        enddo
!
      elseif(nbfijkl(4) == 7)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do k= 1,nbfijkl(3)
              do l= 1,10
                work(l)= eritmp(l,k,j,i)
              enddo
              twoeri(1,k,j,i)=(-work(7)+three*work(2)                   )*facf1
              twoeri(2,k,j,i)=  work(5)                                  *facf2
              twoeri(3,k,j,i)=(-work(7)-work(2)+four*work(9)            )*facf3
              twoeri(4,k,j,i)=( two*work(10)-three*work(3)-three*work(8))*half
              twoeri(5,k,j,i)=(-work(1)-work(4)+four*work(6)            )*facf3
              twoeri(6,k,j,i)=( work(3)-work(8)                         )*facf4
              twoeri(7,k,j,i)=( work(1)-three*work(4)                   )*facf1
            enddo
          enddo
        enddo
      elseif(nbfijkl(4) == 10)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do k= 1,nbfijkl(3)
              twoeri( 1,k,j,i)= eritmp( 1,k,j,i)
              twoeri( 2,k,j,i)= eritmp( 2,k,j,i)*sqrt5
              twoeri( 3,k,j,i)= eritmp( 3,k,j,i)*sqrt5
              twoeri( 4,k,j,i)= eritmp( 4,k,j,i)*sqrt5
              twoeri( 5,k,j,i)= eritmp( 5,k,j,i)*sqrt15
              twoeri( 6,k,j,i)= eritmp( 6,k,j,i)*sqrt5
              twoeri( 7,k,j,i)= eritmp( 7,k,j,i)
              twoeri( 8,k,j,i)= eritmp( 8,k,j,i)*sqrt5
              twoeri( 9,k,j,i)= eritmp( 9,k,j,i)*sqrt5
              twoeri(10,k,j,i)= eritmp(10,k,j,i)
            enddo
          enddo
        enddo
!
      elseif(nbfijkl(4) == 9)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do k= 1,nbfijkl(3)
              do l= 1,15
                work(l)= eritmp(l,k,j,i)
              enddo
              twoeri(1,k,j,i)=(work(2)-work(7))*facg1
              twoeri(2,k,j,i)=(-work(12)+work(5)*three)*facg2
              twoeri(3,k,j,i)=(-work(2)-work(7)+work(9)*six)*facg3
              twoeri(4,k,j,i)=(-work(12)*three+work(14)*four-work(5)*three)*facg4
              twoeri(5,k,j,i)=(work(1)*three+work(11)*three+work(15)*eight+work(4)*six &
&                              -work(6)*p24-work(13)*p24)*eighth
              twoeri(6,k,j,i)=(-work(3)*three+work(10)*four-work(8)*three)*facg4
              twoeri(7,k,j,i)=(-work(1)+work(11)+work(6)*six-work(13)*six)*facg5
              twoeri(8,k,j,i)=(work(3)-work(8)*three)*facg2
              twoeri(9,k,j,i)=(work(1)+work(11)-work(4)*six)*facg6
            enddo
          enddo
        enddo
      elseif(nbfijkl(4) == 15)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do k= 1,nbfijkl(3)
              twoeri( 1,k,j,i)= eritmp( 1,k,j,i)
              twoeri( 2,k,j,i)= eritmp( 2,k,j,i)*sqrt7
              twoeri( 3,k,j,i)= eritmp( 3,k,j,i)*sqrt7
              twoeri( 4,k,j,i)= eritmp( 4,k,j,i)*sqrt35third
              twoeri( 5,k,j,i)= eritmp( 5,k,j,i)*sqrt35
              twoeri( 6,k,j,i)= eritmp( 6,k,j,i)*sqrt35third
              twoeri( 7,k,j,i)= eritmp( 7,k,j,i)*sqrt7
              twoeri( 8,k,j,i)= eritmp( 8,k,j,i)*sqrt35
              twoeri( 9,k,j,i)= eritmp( 9,k,j,i)*sqrt35
              twoeri(10,k,j,i)= eritmp(10,k,j,i)*sqrt7
              twoeri(11,k,j,i)= eritmp(11,k,j,i)
              twoeri(12,k,j,i)= eritmp(12,k,j,i)*sqrt7
              twoeri(13,k,j,i)= eritmp(13,k,j,i)*sqrt35third
              twoeri(14,k,j,i)= eritmp(14,k,j,i)*sqrt7
              twoeri(15,k,j,i)= eritmp(15,k,j,i)
            enddo
          enddo
        enddo
!
      elseif(nbfijkl(4) == 11)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do k= 1,nbfijkl(3)
              do l= 1,21
                work(l)= eritmp(l,k,j,i)
              enddo
              twoeri( 1,k,j,i)=(work(2)*five-work(7)*ten+work(16))*fach1
              twoeri( 2,k,j,i)=(work(5)*four-work(12)*four)*fach2
              twoeri( 3,k,j,i)=(-work(2)*three-work(7)*two+work(9)*p24+work(16) &
&                              -work(18)*eight)*fach3
              twoeri( 4,k,j,i)=(-work(5)*two-work(12)*two+work(14)*four)*fach4
              twoeri( 5,k,j,i)=(work(2)+work(7)*two-work(9)*twelve+work(16)-work(18)*twelve &
&                              +work(20)*eight)*fach5
              twoeri( 6,k,j,i)=(work(3)*p15+work(8)*p30-work(10)*p40+work(17)*p15-work(19)*p40 &
&                              +work(21)*eight)*eighth
              twoeri( 7,k,j,i)=(work(1)+work(4)*two-work(6)*twelve+work(11)-work(13)*twelve &
&                              +work(15)*eight)*fach5
              twoeri( 8,k,j,i)=(-work(3)+work(10)*two+work(17)-work(19)*two)*fach4
              twoeri( 9,k,j,i)=(-work(1)+work(4)*two+work(6)*eight+work(11)*three &
&                              -work(13)*p24)*fach3
              twoeri(10,k,j,i)=(work(3)-work(8)*six+work(17))*fach2
              twoeri(11,k,j,i)=(work(1)-work(4)*ten+work(11)*five)*fach1
            enddo
          enddo
        enddo
      elseif(nbfijkl(4) == 21)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do k= 1,nbfijkl(3)
              twoeri( 1,k,j,i)= eritmp( 1,k,j,i)
              twoeri( 2,k,j,i)= eritmp( 2,k,j,i)*three
              twoeri( 3,k,j,i)= eritmp( 3,k,j,i)*three
              twoeri( 4,k,j,i)= eritmp( 4,k,j,i)*sqrt21
              twoeri( 5,k,j,i)= eritmp( 5,k,j,i)*sqrt63
              twoeri( 6,k,j,i)= eritmp( 6,k,j,i)*sqrt21
              twoeri( 7,k,j,i)= eritmp( 7,k,j,i)*sqrt21
              twoeri( 8,k,j,i)= eritmp( 8,k,j,i)*sqrt105
              twoeri( 9,k,j,i)= eritmp( 9,k,j,i)*sqrt105
              twoeri(10,k,j,i)= eritmp(10,k,j,i)*sqrt21
              twoeri(11,k,j,i)= eritmp(11,k,j,i)*three
              twoeri(12,k,j,i)= eritmp(12,k,j,i)*sqrt63
              twoeri(13,k,j,i)= eritmp(13,k,j,i)*sqrt105
              twoeri(14,k,j,i)= eritmp(14,k,j,i)*sqrt63
              twoeri(15,k,j,i)= eritmp(15,k,j,i)*three
              twoeri(16,k,j,i)= eritmp(16,k,j,i)
              twoeri(17,k,j,i)= eritmp(17,k,j,i)*three
              twoeri(18,k,j,i)= eritmp(18,k,j,i)*sqrt21
              twoeri(19,k,j,i)= eritmp(19,k,j,i)*sqrt21
              twoeri(20,k,j,i)= eritmp(20,k,j,i)*three
              twoeri(21,k,j,i)= eritmp(21,k,j,i)
            enddo
          enddo
        enddo
!
      elseif(nbfijkl(4) == 13)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do k= 1,nbfijkl(3)
              do l= 1,28
                work(l)= eritmp(l,k,j,i)
              enddo
              twoeri( 1,k,j,i)=(work(2)*six-work(7)*p20+work(16)*six)*faci1
              twoeri( 2,k,j,i)=(work(5)*five-work(12)*ten+work(23))*faci2
              twoeri( 3,k,j,i)=(-work(2)*four+work(9)*p40+work(16)*four-work(18)*p40)*faci3
              twoeri( 4,k,j,i)=(-work(5)*p9-work(12)*six+work(14)*p24+work(23)*three &
&                              -work(25)*eight)*faci4
              twoeri( 5,k,j,i)=(work(2)*two+work(7)*four-work(9)*p32+work(16)*two &
&                              -work(18)*p32+work(20)*p32)*faci5
              twoeri( 6,k,j,i)=(work(5)*five+work(12)*ten-work(14)*p20+work(23)*five &
&                              -work(25)*p20+work(27)*eight)*faci6
              twoeri( 7,k,j,i)=(-work(1)*five-work(4)*p15+work(6)*p90-work(11)*p15 &
&                              +work(13)*p180-work(15)*p120-work(22)*five+work(24)*p90 &
&                              -work(26)*p120+work(28)*p16)*sixteenth
              twoeri( 8,k,j,i)=(work(3)*five+work(8)*ten-work(10)*p20+work(17)*five &
&                              -work(19)*p20+work(21)*eight)*faci6
              twoeri( 9,k,j,i)=(work(1)+work(4)-work(6)*p16-work(11)+work(15)*p16 &
&                              -work(22)+work(24)*p16-work(26)*p16)*faci5
              twoeri(10,k,j,i)=(-work(3)*three+work(8)*six+work(10)*eight+work(17)*p9 &
&                              -work(19)*p24)*faci4
              twoeri(11,k,j,i)=(-work(1)+work(4)*five+work(6)*ten+work(11)*five &
&                              -work(13)*p60-work(22)+work(24)*ten)*faci3
              twoeri(12,k,j,i)=(work(3)-work(8)*ten+work(17)*five)*faci2
              twoeri(13,k,j,i)=(work(1)-work(4)*p15+work(11)*p15-work(22))*faci1
            enddo
          enddo
        enddo
      elseif(nbfijkl(4) == 28)then
        do i= 1,nbfijkl(1)
          do j= 1,nbfijkl(2)
            do k= 1,nbfijkl(3)
              twoeri( 1,k,j,i)= eritmp( 1,k,j,i)
              twoeri( 2,k,j,i)= eritmp( 2,k,j,i)*sqrt11
              twoeri( 3,k,j,i)= eritmp( 3,k,j,i)*sqrt11
              twoeri( 4,k,j,i)= eritmp( 4,k,j,i)*sqrt33
              twoeri( 5,k,j,i)= eritmp( 5,k,j,i)*sqrt99
              twoeri( 6,k,j,i)= eritmp( 6,k,j,i)*sqrt33
              twoeri( 7,k,j,i)= eritmp( 7,k,j,i)*sqrt231fifth
              twoeri( 8,k,j,i)= eritmp( 8,k,j,i)*sqrt231
              twoeri( 9,k,j,i)= eritmp( 9,k,j,i)*sqrt231
              twoeri(10,k,j,i)= eritmp(10,k,j,i)*sqrt231fifth
              twoeri(11,k,j,i)= eritmp(11,k,j,i)*sqrt33
              twoeri(12,k,j,i)= eritmp(12,k,j,i)*sqrt231
              twoeri(13,k,j,i)= eritmp(13,k,j,i)*sqrt385
              twoeri(14,k,j,i)= eritmp(14,k,j,i)*sqrt231
              twoeri(15,k,j,i)= eritmp(15,k,j,i)*sqrt33
              twoeri(16,k,j,i)= eritmp(16,k,j,i)*sqrt11
              twoeri(17,k,j,i)= eritmp(17,k,j,i)*sqrt99
              twoeri(18,k,j,i)= eritmp(18,k,j,i)*sqrt231
              twoeri(19,k,j,i)= eritmp(19,k,j,i)*sqrt231
              twoeri(20,k,j,i)= eritmp(20,k,j,i)*sqrt99
              twoeri(21,k,j,i)= eritmp(21,k,j,i)*sqrt11
              twoeri(22,k,j,i)= eritmp(22,k,j,i)
              twoeri(23,k,j,i)= eritmp(23,k,j,i)*sqrt11
              twoeri(24,k,j,i)= eritmp(24,k,j,i)*sqrt33
              twoeri(25,k,j,i)= eritmp(25,k,j,i)*sqrt231fifth
              twoeri(26,k,j,i)= eritmp(26,k,j,i)*sqrt33
              twoeri(27,k,j,i)= eritmp(27,k,j,i)*sqrt11
              twoeri(28,k,j,i)= eritmp(28,k,j,i)
            enddo
          enddo
        enddo
      endif
!
      return
end
