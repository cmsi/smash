! Copyright 2014  Kazuya Ishimura
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
!----------------------------------------------------------------
  subroutine calc2eri(twoeri,ish,jsh,ksh,lsh,maxdim,lrint,emu2)
!----------------------------------------------------------------
!
! Driver of two-electron repulsion integrals 
! 
! In  : ish, jsh, ksh, lsh (Shell indices)
!       lrint  (Flag of long range correction)
!       emu2   (mu*mu value for long range correction)
! Out : twoeri (Two-electron repulsion integrals)
!       idxeri (Indices of two-electron repulsion integrals)
!       numeri (Number of two-electron repulsion integrals)
!
      use modparam, only : mxprsh
      use modmolecule, only : coord
      use modbasis, only : locatom, locprim, mprim, mbf, mtype, ex, coeff
      use modthresh, only : threshex
      implicit none
      integer,intent(in) :: ish, jsh, ksh, lsh, maxdim
      integer :: nangijkl(4), nprimijkl(4), nbfijkl(4)
      integer :: iloc, jloc, kloc, lloc, i, j, k, l
      real(8),intent(in) :: emu2
      real(8),intent(out) :: twoeri(maxdim,maxdim,maxdim,maxdim)
      real(8) :: xyzijkl(3,4), exijkl(mxprsh,4), coijkl(mxprsh,4)
      logical,intent(in) :: lrint
!
      if(lrint) then
        write(*,'(" Error! Long range correction is not supported now.")')
        write(*,'(" emu2=",d10.3)')emu2
      endif
!
      nangijkl(1)= mtype(ish)
      nangijkl(2)= mtype(jsh)
      nangijkl(3)= mtype(ksh)
      nangijkl(4)= mtype(lsh)
!
      iloc  = locprim(ish)
      jloc  = locprim(jsh)
      kloc  = locprim(ksh)
      lloc  = locprim(lsh)
      nbfijkl(1)= mbf(ish)
      nbfijkl(2)= mbf(jsh)
      nbfijkl(3)= mbf(ksh)
      nbfijkl(4)= mbf(lsh)
      nprimijkl(1)= mprim(ish)
      nprimijkl(2)= mprim(jsh)
      nprimijkl(3)= mprim(ksh)
      nprimijkl(4)= mprim(lsh)
!
      do i= 1,3
        xyzijkl(i,1)= coord(i,locatom(ish))
        xyzijkl(i,2)= coord(i,locatom(jsh))
        xyzijkl(i,3)= coord(i,locatom(ksh))
        xyzijkl(i,4)= coord(i,locatom(lsh))
      enddo
!
      do i= 1,nprimijkl(1)
        exijkl(i,1)=ex(iloc+i)
        coijkl(i,1)=coeff(iloc+i)
      enddo
      do j= 1,nprimijkl(2)
        exijkl(j,2)=ex(jloc+j)
        coijkl(j,2)=coeff(jloc+j)
      enddo
      do k= 1,nprimijkl(3)
        exijkl(k,3)=ex(kloc+k)
        coijkl(k,3)=coeff(kloc+k)
      enddo
      do l= 1,nprimijkl(4)
        exijkl(l,4)=ex(lloc+l)
        coijkl(l,4)=coeff(lloc+l)
      enddo
!
      if((nangijkl(1)<=2).and.(nangijkl(2)<=2).and.(nangijkl(3)<=2).and.(nangijkl(4)<=2)) then
!
! Pople-Hehre and McMurchie-Davidson scheme
!
        call int2phmd(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                     mxprsh,threshex,lrint,emu2)
      else
!
! Rys quadrature scheme
!
        call int2rys(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                    mxprsh,threshex,lrint,emu2)
      endif
      return
end


!----------------------------------------------------------------------------------------
  subroutine int2phmd(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                     mxprsh,threshex,lrint,emu2)
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
!       lrint     (Flag of long range correction)
!       emu2      (mu*mu value for long range correction)
! Out : twoeri    (Two-electron integrals
!                  The values are stored in order of (lsh,ksh,jsh,ish).)
!
      implicit none
      integer,intent(in) :: nprimijkl(4), nangijkl(4), nbfijkl(4), maxdim, mxprsh
      integer,parameter :: mxprsh2=40
      integer :: nangtotal, inttype, intorder
      integer :: iprim, jprim, kprim, lprim, i, j, k, l, nijkl(2)
      integer :: intijkl(4,8), inttable(81,2), ii, jj, kk, ll, nbfijkl2(4)
      real(8),parameter :: zero=0.0d+00, one=1.0D+00, half=0.5D+00, pt7=0.7D+00
      real(8),parameter :: pi52=0.3498683665524973D+02 ! 2.0*pi**2.5
      real(8),intent(in) :: exijkl(mxprsh,4), coijkl(mxprsh,4), xyzijkl(3,4), threshex
      real(8),intent(in) :: emu2
      real(8),intent(out) :: twoeri(maxdim,maxdim,maxdim,maxdim)
      real(8) :: phmdint(6,6,6,6)
      real(8) :: r12, r34, rij, rkl, tmp, rot(3,3), veckl(3), cosph, sinph, abscos
      real(8) :: xyzij(3), xyzkl(3), xyzik(3), xzkl(2)
      real(8) :: exi, exj, ci, cj, ex12, ex21, ex1p, ex2p, r12ex
      real(8) :: exk, exl, ck, cl, ex34, ex43, ex3q, ex4q, r34ex
      real(8) :: exfac(5,mxprsh2*mxprsh2*2), xyziq(3,mxprsh2*mxprsh2)
      real(8) :: tmpexp(mxprsh2*mxprsh2)
      logical,intent(in) :: lrint
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
      if(lrint) then
        write(*,'(" Error! Long range correction is not supported now.")')
        write(*,'(" emu2=",d10.3)')emu2
      endif
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
        if(abs(rot(1,3)).le.pt7) then
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
        do jprim= 1,nprimijkl(jj)
          exj = exijkl(jprim,jj)
          cj  = coijkl(jprim,jj)
          ex12= exi+exj
          ex21= one/ex12
          ex1p= exi*ex21
          ex2p= exj*ex21
          r12ex= r12*exj*ex1p
          if(r12ex.gt.threshex)cycle
          nijkl(1)= nijkl(1)+1
          exfac(1,nijkl(1))= ex12
          exfac(2,nijkl(1))= ex21*half
          exfac(3,nijkl(1))=-ex1p*rij
          exfac(4,nijkl(1))= ex2p*rij
!         exfac(5,nijkl(1))= ex21*ci*cj*exp(-r12ex)
          exfac(5,nijkl(1))= ex21*ci*cj
          tmpexp(nijkl(1))=-r12ex
        enddo
      enddo
!
!     do i=1,nijkl(1)
!       exfac1(5,i)=exfac1(5,i)*exp(tmpexp(i))
!     enddo
!
!     nijkl(2)= 0
      do kprim= 1,nprimijkl(kk)
        exk= exijkl(kprim,kk)
        ck = coijkl(kprim,kk)
        do lprim= 1,nprimijkl(ll)
          exl = exijkl(lprim,ll)
          cl  = coijkl(lprim,ll)
          ex34= exk+exl
          ex43= one/ex34
          ex3q= exk*ex43
          ex4q= exl*ex43
          r34ex= r34*exl*ex3q
          if(r34ex.gt.threshex)cycle
          nijkl(2)= nijkl(2)+1
          exfac(1,nijkl(1)+nijkl(2))= ex34
          exfac(2,nijkl(1)+nijkl(2))= ex43*half
          exfac(3,nijkl(1)+nijkl(2))=-ex3q
          exfac(4,nijkl(1)+nijkl(2))= ex4q
!         exfac(5,nijkl(1)+nijkl(2))= ex43*ck*cl*exp(-r34ex)
          exfac(5,nijkl(1)+nijkl(2))= ex43*ck*cl
          tmpexp(nijkl(1)+nijkl(2))=-r34ex
          xyziq(1,nijkl(2))= xyzik(1)+xzkl(1)*ex4q
          xyziq(2,nijkl(2))= xyzik(2)
          xyziq(3,nijkl(2))= xyzik(3)+xzkl(2)*ex4q
        enddo
      enddo
!
      do i=1,nijkl(1)+nijkl(2)
        exfac(5,i)=exfac(5,i)*exp(tmpexp(i))
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
          call int2ssss(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,nijkl)
        case (2)
          call int2psss(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl)
        case (3)
          call int2ppss(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl)
        case (4)
          call int2psps(phmdint,exfac,exfac(1,nijkl(1)+1),xyziq,xzkl,rot,nijkl)
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


!-----------------------------------------------------------------------
  subroutine calcschwarzeri(xint,xinttmp,maxdim,nproc,myrank,mpi_comm)
!-----------------------------------------------------------------------
!
! Driver of (ij|ij) integral calculation
!
! Out : xint    ((ij|ij) integrals for Schwarz screening)
!       xinttmp (Work array)
!
      use modbasis, only : nshell, mtype, mbf
      use modthresh, only : cutint2, threshex
      implicit none
      integer,intent(in) :: maxdim, nproc, myrank
      integer(4),intent(in) :: mpi_comm
      integer :: ish, jsh, nbfi, nbfj, i, j, ii, ij
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, two=2.0D+00
      real(8),intent(out) :: xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: xinttmp(nshell*(nshell+1)/2)
      real(8) :: twoeri(maxdim,maxdim,maxdim,maxdim), xintmax, cutsave, val
!
      cutsave= cutint2
      cutint2= cutint2*cutint2
      threshex=threshex*two
!
      xinttmp(:)= zero
!
!$OMP parallel private(jsh,xintmax,twoeri,nbfi,nbfj,i,j,val,ij,ii)
      do ish= nshell-myrank,1,-nproc
        nbfi= mbf(ish)
        ii  = ish*(ish-1)/2
!$OMP do
        do jsh = 1,ish
          nbfj = mbf(jsh)
          call calc2eri(twoeri,ish,jsh,ish,jsh,maxdim,.false.,zero)
          xintmax= zero
          do i= 1,nbfi
            do j= 1,nbfj
              val= twoeri(j,i,j,i)
              if(val.gt.xintmax) xintmax= val
            enddo
          enddo
          ij= ii+jsh
          xinttmp(ij)=sqrt(xintmax)
        enddo
!$OMP end do
      enddo
!$OMP end parallel
!
      call para_allreducer(xinttmp,xint,nshell*(nshell+1)/2,mpi_comm)
!
      cutint2= cutsave
      threshex=threshex*half
      return
end

