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
!----------------------------------------------------------------------------
  subroutine grad2eri(egrad,egrad2,fulldmtrx1,fulldmtrx2,xint,hfexchange, &
&                     maxdim,maxgraddim,nproc,myrank,itype)
!----------------------------------------------------------------------------
!
! Main driver of derivatives for two-electron integrals
!
! In    : fulldmtrx1(Full alpha density matrix (itype=1)
!                    Full alpha+beta (itype=2)
!                    Full HF alpha density matrix (itype=3))
!         fulldmtrx2(Full alpha density matrix (itype=1)
!                    Full alpha-beta (itype=2)
!                    Full MP2 alpha density matrix (itype=3))
!         maxdim    (Maximum dimension of pdmtrx and dtwoeri)
!         maxgraddim(Maximum dimension of twoeri)
!         hfexchange(Hartree-Fock exchange scaling factor)
!         itype     (1:RHF, 2:UHF)
! Inout : egrad2    (Energy gradient values)
!
      use modbasis, only : nshell, nao
      use modthresh, only : cutint2
      use modmolecule, only : natom
      implicit none
      integer,intent(in) :: maxdim, maxgraddim, nproc, myrank, itype
      integer :: ijsh, ish, jsh, ksh, lsh, ij, kl
      integer :: ii, kk, kstart, ishcheck
      integer(8) :: ncount, icount(nshell)
      real(8),parameter :: zero=0.0D+00, four=4.0D+00
      real(8),intent(in) :: fulldmtrx1(nao,nao), fulldmtrx2(nao,nao)
      real(8),intent(in) :: xint(nshell*(nshell+1)/2), hfexchange
      real(8),intent(out) :: egrad2(3*natom)
      real(8),intent(inout) :: egrad(3*natom)
      real(8) :: xijkl, twoeri(maxgraddim**4), dtwoeri(maxdim**4,3), pdmtrx(maxdim**4)
      real(8) :: pdmax
!
      egrad2(:)= zero
!
      ncount= 0
      ncount= ncount+(2*nshell**3+3*nshell**2+nshell)/6+myrank
      do ish= 1,nshell
        icount(ish)= ncount-(2*ish*ish*ish-3*ish*ish+ish)/6
      enddo
!
      ish= nshell
      ii= ish*(ish-1)/2
!
!$OMP parallel do schedule(dynamic,1) &
!$OMP private(ijsh,jsh,ksh,lsh,ij,kl,xijkl,twoeri,dtwoeri,pdmtrx,pdmax,kk,kstart) &
!$OMP firstprivate(ish,ii) reduction(+:egrad2)
      do ijsh= nshell*(nshell+1)/2,1,-1
        do ishcheck=1,nshell
          if(ijsh > ii) then
            jsh= ijsh-ii
            exit
          else
            ish= ish-1
            ii= ish*(ish-1)/2
          endif
        enddo
!
        ij= ii+jsh
        kstart=mod(icount(ish)-ish*(jsh-1),nproc)+1
 kloop: do ksh= kstart,ish,nproc
          kk= ksh*(ksh-1)/2
          do lsh= 1,ksh
            kl= kk+lsh
            if(kl.gt.ij) exit kloop
            xijkl=xint(ij)*xint(kl)
            if(xijkl.lt.cutint2) cycle
            call calcpdmtrx(fulldmtrx1,fulldmtrx2,pdmtrx,pdmax,hfexchange, &
&                           ish,jsh,ksh,lsh,maxdim,itype)
            if((xijkl*pdmax).lt.cutint2) cycle
            call calcd2eri(egrad2,pdmtrx,twoeri,dtwoeri,ish,jsh,ksh,lsh,maxdim,maxgraddim)
          enddo
        enddo kloop
      enddo
!$OMP end parallel do
!
      do ii= 1,3*natom
        egrad(ii)= egrad(ii)+egrad2(ii)*four
      enddo
      return
end


!---------------------------------------------------------------------------------------
  subroutine calcd2eri(egrad2,pdmtrx,twoeri,dtwoeri,ish,jsh,ksh,lsh,maxdim,maxgraddim)
!---------------------------------------------------------------------------------------
!
! Driver of derivatives for two-electron integrals
!
! In    : pdmtrx    (Products of density matrix)
!         maxdim    (Maximum dimension of pdmtrx and dtwoeri)
!         maxgraddim(Maximum dimension of twoeri)
!         ish,jsh,ksh,lsh (Shell indices)
! Out   : twoeri    (Derivatives for two-electron repulsion integrals)
!         dtwoeri   (Derivatives for two-electron repulsion integrals)
! Inout : egrad2    (Energy gradient values)
!
      use modparam, only : mxprsh
      use modmolecule, only : coord, natom
      use modbasis, only : locatom, locprim, mprim, mbf, mtype, ex, coeff
      use modthresh, only : threshex
      implicit none
      integer,parameter :: ncart(0:6)=(/1,3,6,10,15,21,28/)
      integer,intent(in) :: ish, jsh, ksh, lsh, maxdim, maxgraddim
      integer :: i, j, k, l, iatom, jatom, katom, latom, iloc, jloc, kloc, lloc, ider
      integer :: nangijkl(4), nbfijkl(4), nprimijkl(4)
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, two=2.0D+00, three=3.0D+00
      real(8),parameter :: four=4.0D+00, five=5.0D+00, six=6.0D+00, eight=8.0D+00, ten=10.0D+00
      real(8),parameter :: twelve=12.0D+00, p15=15.0D+00, p24=24.0D+00, p30=30.0D+00
      real(8),parameter :: p40=40.0D+00
      real(8),parameter :: third=3.333333333333333D-01, eighth=0.125D+00
      real(8),parameter :: sqrt3=1.732050807568877D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrtthird=0.5773502691896258D+00, sqrtfifth=0.4472135954999579D+00
      real(8),parameter :: sqrt3fifth=0.7745966692414834D+00, sqrtseventh=0.3779644730092272D+00 
      real(8),parameter :: sqrtinv35=0.1690308509457033D+00, sqrt3inv35=0.2927700218845599D+00
      real(8),parameter :: sqrt5=2.236067977499790D+00, sqrt15=3.872983346207417D+00
      real(8),parameter :: sqrt7=2.645751311064591D+00, sqrt35=5.916079783099616D+00
      real(8),parameter :: sqrt35third=3.415650255319866D+00
      real(8),parameter :: sqrtinv15=2.581988897471611D-01, sqrtinv21=2.182178902359924D-01
      real(8),parameter :: sqrtinv63=1.259881576697424D-01, sqrtinv105=9.759000729485332D-02
      real(8),parameter :: sqrtinv11=3.015113445777636D-01, sqrtinv33=1.740776559556978D-01
      real(8),parameter :: sqrtinv99=1.005037815259212D-01, sqrtinv231=6.579516949597690D-02
      real(8),parameter :: sqrtinv385=5.096471914376255D-02, sqrt5inv231=1.471224715841249D-01
      real(8),parameter :: sqrt21=4.582575694955840D+00, sqrt63=7.937253933193772D+00
      real(8),parameter :: sqrt105=1.024695076595960D+01
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
      real(8),intent(in) :: pdmtrx(maxdim,maxdim,maxdim,maxdim)
      real(8),intent(out) :: twoeri(maxgraddim,maxgraddim,maxgraddim,maxgraddim)
      real(8),intent(out) :: dtwoeri(maxdim,maxdim,maxdim,maxdim,3)
      real(8),intent(inout) :: egrad2(3,natom)
      real(8) :: gradtwo(3,4), xyzijkl(3,4), exijkl(mxprsh,4), coijkl(mxprsh,4), work(21)
!
      iatom= locatom(ish)
      jatom= locatom(jsh)
      katom= locatom(ksh)
      latom= locatom(lsh)
!
      if((iatom == jatom).and.(iatom == katom).and.(iatom == latom)) return
      gradtwo(:,:)= zero
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
        xyzijkl(i,1)= coord(i,iatom)
        xyzijkl(i,2)= coord(i,jatom)
        xyzijkl(i,3)= coord(i,katom)
        xyzijkl(i,4)= coord(i,latom)
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
! Lsh derivative 
!
      nangijkl(4)= mtype(lsh)+1
      nbfijkl(4) = ncart(nangijkl(4))
      do l= 1,nprimijkl(4)
        coijkl(l,4)= two*ex(lloc+l)*coeff(lloc+l)
      enddo
!
! Two-electron integral calculation
!
      call int2elec(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxgraddim, &
&                   mxprsh,threshex)
!
      select case(nangijkl(4))
        case(1)
          do i= 1,nbfijkl(1)
            do j= 1,nbfijkl(2)
              do k= 1,nbfijkl(3)
                dtwoeri(1,k,j,i,1)= twoeri(1,k,j,i)
                dtwoeri(1,k,j,i,2)= twoeri(2,k,j,i)
                dtwoeri(1,k,j,i,3)= twoeri(3,k,j,i)
              enddo
            enddo
          enddo
        case(2)
          do i= 1,nbfijkl(1)
            do j= 1,nbfijkl(2)
              do k= 1,nbfijkl(3)
                dtwoeri(1,k,j,i,1)= twoeri(1,k,j,i)
                dtwoeri(2,k,j,i,1)= twoeri(2,k,j,i)*sqrtthird
                dtwoeri(3,k,j,i,1)= twoeri(3,k,j,i)*sqrtthird
                dtwoeri(1,k,j,i,2)= twoeri(2,k,j,i)*sqrtthird
                dtwoeri(2,k,j,i,2)= twoeri(4,k,j,i)
                dtwoeri(3,k,j,i,2)= twoeri(5,k,j,i)*sqrtthird
                dtwoeri(1,k,j,i,3)= twoeri(3,k,j,i)*sqrtthird
                dtwoeri(2,k,j,i,3)= twoeri(5,k,j,i)*sqrtthird
                dtwoeri(3,k,j,i,3)= twoeri(6,k,j,i)
              enddo
            enddo
          enddo
        case(3)
          do i= 1,nbfijkl(1)
            do j= 1,nbfijkl(2)
              do k= 1,nbfijkl(3)
                dtwoeri(1,k,j,i,1)= twoeri( 1,k,j,i)
                dtwoeri(2,k,j,i,1)= twoeri( 2,k,j,i)*sqrt3fifth
                dtwoeri(3,k,j,i,1)= twoeri( 3,k,j,i)*sqrt3fifth
                dtwoeri(4,k,j,i,1)= twoeri( 4,k,j,i)*sqrtfifth
                dtwoeri(5,k,j,i,1)= twoeri( 5,k,j,i)*sqrtfifth
                dtwoeri(6,k,j,i,1)= twoeri( 6,k,j,i)*sqrtfifth
                dtwoeri(1,k,j,i,2)= twoeri( 2,k,j,i)*sqrtfifth
                dtwoeri(2,k,j,i,2)= twoeri( 4,k,j,i)*sqrt3fifth
                dtwoeri(3,k,j,i,2)= twoeri( 5,k,j,i)*sqrtfifth
                dtwoeri(4,k,j,i,2)= twoeri( 7,k,j,i)
                dtwoeri(5,k,j,i,2)= twoeri( 8,k,j,i)*sqrt3fifth
                dtwoeri(6,k,j,i,2)= twoeri( 9,k,j,i)*sqrtfifth
                dtwoeri(1,k,j,i,3)= twoeri( 3,k,j,i)*sqrtfifth
                dtwoeri(2,k,j,i,3)= twoeri( 5,k,j,i)*sqrtfifth
                dtwoeri(3,k,j,i,3)= twoeri( 6,k,j,i)*sqrt3fifth
                dtwoeri(4,k,j,i,3)= twoeri( 8,k,j,i)*sqrtfifth
                dtwoeri(5,k,j,i,3)= twoeri( 9,k,j,i)*sqrt3fifth
                dtwoeri(6,k,j,i,3)= twoeri(10,k,j,i)
              enddo
            enddo
          enddo
        case(4)
          do i= 1,nbfijkl(1)
            do j= 1,nbfijkl(2)
              do k= 1,nbfijkl(3)
                dtwoeri( 1,k,j,i,1)= twoeri( 1,k,j,i)
                dtwoeri( 2,k,j,i,1)= twoeri( 2,k,j,i)*sqrtseventh
                dtwoeri( 3,k,j,i,1)= twoeri( 3,k,j,i)*sqrtseventh
                dtwoeri( 4,k,j,i,1)= twoeri( 4,k,j,i)*sqrt3inv35
                dtwoeri( 5,k,j,i,1)= twoeri( 5,k,j,i)*sqrtinv35
                dtwoeri( 6,k,j,i,1)= twoeri( 6,k,j,i)*sqrt3inv35
                dtwoeri( 7,k,j,i,1)= twoeri( 7,k,j,i)*sqrtseventh
                dtwoeri( 8,k,j,i,1)= twoeri( 8,k,j,i)*sqrtinv35
                dtwoeri( 9,k,j,i,1)= twoeri( 9,k,j,i)*sqrtinv35
                dtwoeri(10,k,j,i,1)= twoeri(10,k,j,i)*sqrtseventh
                dtwoeri( 1,k,j,i,2)= twoeri( 2,k,j,i)*sqrtseventh
                dtwoeri( 2,k,j,i,2)= twoeri( 4,k,j,i)*sqrt3inv35
                dtwoeri( 3,k,j,i,2)= twoeri( 5,k,j,i)*sqrtinv35
                dtwoeri( 4,k,j,i,2)= twoeri( 7,k,j,i)*sqrtseventh
                dtwoeri( 5,k,j,i,2)= twoeri( 8,k,j,i)*sqrtinv35
                dtwoeri( 6,k,j,i,2)= twoeri( 9,k,j,i)*sqrtinv35
                dtwoeri( 7,k,j,i,2)= twoeri(11,k,j,i)
                dtwoeri( 8,k,j,i,2)= twoeri(12,k,j,i)*sqrtseventh
                dtwoeri( 9,k,j,i,2)= twoeri(13,k,j,i)*sqrt3inv35
                dtwoeri(10,k,j,i,2)= twoeri(14,k,j,i)*sqrtseventh
                dtwoeri( 1,k,j,i,3)= twoeri( 3,k,j,i)*sqrtseventh
                dtwoeri( 2,k,j,i,3)= twoeri( 5,k,j,i)*sqrtinv35
                dtwoeri( 3,k,j,i,3)= twoeri( 6,k,j,i)*sqrt3inv35
                dtwoeri( 4,k,j,i,3)= twoeri( 8,k,j,i)*sqrtinv35
                dtwoeri( 5,k,j,i,3)= twoeri( 9,k,j,i)*sqrtinv35
                dtwoeri( 6,k,j,i,3)= twoeri(10,k,j,i)*sqrtseventh
                dtwoeri( 7,k,j,i,3)= twoeri(12,k,j,i)*sqrtseventh
                dtwoeri( 8,k,j,i,3)= twoeri(13,k,j,i)*sqrt3inv35
                dtwoeri( 9,k,j,i,3)= twoeri(14,k,j,i)*sqrtseventh
                dtwoeri(10,k,j,i,3)= twoeri(15,k,j,i)
              enddo
            enddo
          enddo
        case(5)
          do i= 1,nbfijkl(1)
            do j= 1,nbfijkl(2)
              do k= 1,nbfijkl(3)
                dtwoeri( 1,k,j,i,1)= twoeri( 1,k,j,i)
                dtwoeri( 2,k,j,i,1)= twoeri( 2,k,j,i)*third
                dtwoeri( 3,k,j,i,1)= twoeri( 3,k,j,i)*third
                dtwoeri( 4,k,j,i,1)= twoeri( 4,k,j,i)*sqrtinv21
                dtwoeri( 5,k,j,i,1)= twoeri( 5,k,j,i)*sqrtinv63
                dtwoeri( 6,k,j,i,1)= twoeri( 6,k,j,i)*sqrtinv21
                dtwoeri( 7,k,j,i,1)= twoeri( 7,k,j,i)*sqrtinv21
                dtwoeri( 8,k,j,i,1)= twoeri( 8,k,j,i)*sqrtinv105
                dtwoeri( 9,k,j,i,1)= twoeri( 9,k,j,i)*sqrtinv105
                dtwoeri(10,k,j,i,1)= twoeri(10,k,j,i)*sqrtinv21
                dtwoeri(11,k,j,i,1)= twoeri(11,k,j,i)*third
                dtwoeri(12,k,j,i,1)= twoeri(12,k,j,i)*sqrtinv63
                dtwoeri(13,k,j,i,1)= twoeri(13,k,j,i)*sqrtinv105
                dtwoeri(14,k,j,i,1)= twoeri(14,k,j,i)*sqrtinv63
                dtwoeri(15,k,j,i,1)= twoeri(15,k,j,i)*third
                dtwoeri( 1,k,j,i,2)= twoeri( 2,k,j,i)*third
                dtwoeri( 2,k,j,i,2)= twoeri( 4,k,j,i)*sqrtinv21
                dtwoeri( 3,k,j,i,2)= twoeri( 5,k,j,i)*sqrtinv63
                dtwoeri( 4,k,j,i,2)= twoeri( 7,k,j,i)*sqrtinv21
                dtwoeri( 5,k,j,i,2)= twoeri( 8,k,j,i)*sqrtinv105
                dtwoeri( 6,k,j,i,2)= twoeri( 9,k,j,i)*sqrtinv105
                dtwoeri( 7,k,j,i,2)= twoeri(11,k,j,i)*third
                dtwoeri( 8,k,j,i,2)= twoeri(12,k,j,i)*sqrtinv63
                dtwoeri( 9,k,j,i,2)= twoeri(13,k,j,i)*sqrtinv105
                dtwoeri(10,k,j,i,2)= twoeri(14,k,j,i)*sqrtinv63
                dtwoeri(11,k,j,i,2)= twoeri(16,k,j,i)
                dtwoeri(12,k,j,i,2)= twoeri(17,k,j,i)*third
                dtwoeri(13,k,j,i,2)= twoeri(18,k,j,i)*sqrtinv21
                dtwoeri(14,k,j,i,2)= twoeri(19,k,j,i)*sqrtinv21
                dtwoeri(15,k,j,i,2)= twoeri(20,k,j,i)*third
                dtwoeri( 1,k,j,i,3)= twoeri( 3,k,j,i)*third
                dtwoeri( 2,k,j,i,3)= twoeri( 5,k,j,i)*sqrtinv63
                dtwoeri( 3,k,j,i,3)= twoeri( 6,k,j,i)*sqrtinv21
                dtwoeri( 4,k,j,i,3)= twoeri( 8,k,j,i)*sqrtinv105
                dtwoeri( 5,k,j,i,3)= twoeri( 9,k,j,i)*sqrtinv105
                dtwoeri( 6,k,j,i,3)= twoeri(10,k,j,i)*sqrtinv21
                dtwoeri( 7,k,j,i,3)= twoeri(12,k,j,i)*sqrtinv63
                dtwoeri( 8,k,j,i,3)= twoeri(13,k,j,i)*sqrtinv105
                dtwoeri( 9,k,j,i,3)= twoeri(14,k,j,i)*sqrtinv63
                dtwoeri(10,k,j,i,3)= twoeri(15,k,j,i)*third
                dtwoeri(11,k,j,i,3)= twoeri(17,k,j,i)*third
                dtwoeri(12,k,j,i,3)= twoeri(18,k,j,i)*sqrtinv21
                dtwoeri(13,k,j,i,3)= twoeri(19,k,j,i)*sqrtinv21
                dtwoeri(14,k,j,i,3)= twoeri(20,k,j,i)*third
                dtwoeri(15,k,j,i,3)= twoeri(21,k,j,i)
              enddo
            enddo
          enddo
        case(6)
          do i= 1,nbfijkl(1)
            do j= 1,nbfijkl(2)
              do k= 1,nbfijkl(3)
                dtwoeri( 1,k,j,i,1)= twoeri( 1,k,j,i)
                dtwoeri( 2,k,j,i,1)= twoeri( 2,k,j,i)*sqrtinv11
                dtwoeri( 3,k,j,i,1)= twoeri( 3,k,j,i)*sqrtinv11
                dtwoeri( 4,k,j,i,1)= twoeri( 4,k,j,i)*sqrtinv33
                dtwoeri( 5,k,j,i,1)= twoeri( 5,k,j,i)*sqrtinv99
                dtwoeri( 6,k,j,i,1)= twoeri( 6,k,j,i)*sqrtinv33
                dtwoeri( 7,k,j,i,1)= twoeri( 7,k,j,i)*sqrt5inv231
                dtwoeri( 8,k,j,i,1)= twoeri( 8,k,j,i)*sqrtinv231
                dtwoeri( 9,k,j,i,1)= twoeri( 9,k,j,i)*sqrtinv231
                dtwoeri(10,k,j,i,1)= twoeri(10,k,j,i)*sqrt5inv231
                dtwoeri(11,k,j,i,1)= twoeri(11,k,j,i)*sqrtinv33
                dtwoeri(12,k,j,i,1)= twoeri(12,k,j,i)*sqrtinv231
                dtwoeri(13,k,j,i,1)= twoeri(13,k,j,i)*sqrtinv385
                dtwoeri(14,k,j,i,1)= twoeri(14,k,j,i)*sqrtinv231
                dtwoeri(15,k,j,i,1)= twoeri(15,k,j,i)*sqrtinv33
                dtwoeri(16,k,j,i,1)= twoeri(16,k,j,i)*sqrtinv11
                dtwoeri(17,k,j,i,1)= twoeri(17,k,j,i)*sqrtinv99
                dtwoeri(18,k,j,i,1)= twoeri(18,k,j,i)*sqrtinv231
                dtwoeri(19,k,j,i,1)= twoeri(19,k,j,i)*sqrtinv231
                dtwoeri(20,k,j,i,1)= twoeri(20,k,j,i)*sqrtinv99
                dtwoeri(21,k,j,i,1)= twoeri(21,k,j,i)*sqrtinv11
                dtwoeri( 1,k,j,i,2)= twoeri( 2,k,j,i)*sqrtinv11
                dtwoeri( 2,k,j,i,2)= twoeri( 4,k,j,i)*sqrtinv33
                dtwoeri( 3,k,j,i,2)= twoeri( 5,k,j,i)*sqrtinv99
                dtwoeri( 4,k,j,i,2)= twoeri( 7,k,j,i)*sqrt5inv231
                dtwoeri( 5,k,j,i,2)= twoeri( 8,k,j,i)*sqrtinv231
                dtwoeri( 6,k,j,i,2)= twoeri( 9,k,j,i)*sqrtinv231
                dtwoeri( 7,k,j,i,2)= twoeri(11,k,j,i)*sqrtinv33
                dtwoeri( 8,k,j,i,2)= twoeri(12,k,j,i)*sqrtinv231
                dtwoeri( 9,k,j,i,2)= twoeri(13,k,j,i)*sqrtinv385
                dtwoeri(10,k,j,i,2)= twoeri(14,k,j,i)*sqrtinv231
                dtwoeri(11,k,j,i,2)= twoeri(16,k,j,i)*sqrtinv11 
                dtwoeri(12,k,j,i,2)= twoeri(17,k,j,i)*sqrtinv99 
                dtwoeri(13,k,j,i,2)= twoeri(18,k,j,i)*sqrtinv231
                dtwoeri(14,k,j,i,2)= twoeri(19,k,j,i)*sqrtinv231
                dtwoeri(15,k,j,i,2)= twoeri(20,k,j,i)*sqrtinv99
                dtwoeri(16,k,j,i,2)= twoeri(22,k,j,i)
                dtwoeri(17,k,j,i,2)= twoeri(23,k,j,i)*sqrtinv11
                dtwoeri(18,k,j,i,2)= twoeri(24,k,j,i)*sqrtinv33
                dtwoeri(19,k,j,i,2)= twoeri(25,k,j,i)*sqrt5inv231
                dtwoeri(20,k,j,i,2)= twoeri(26,k,j,i)*sqrtinv33
                dtwoeri(21,k,j,i,2)= twoeri(27,k,j,i)*sqrtinv11
                dtwoeri( 1,k,j,i,3)= twoeri( 3,k,j,i)*sqrtinv11
                dtwoeri( 2,k,j,i,3)= twoeri( 5,k,j,i)*sqrtinv99
                dtwoeri( 3,k,j,i,3)= twoeri( 6,k,j,i)*sqrtinv33
                dtwoeri( 4,k,j,i,3)= twoeri( 8,k,j,i)*sqrtinv231
                dtwoeri( 5,k,j,i,3)= twoeri( 9,k,j,i)*sqrtinv231
                dtwoeri( 6,k,j,i,3)= twoeri(10,k,j,i)*sqrt5inv231
                dtwoeri( 7,k,j,i,3)= twoeri(12,k,j,i)*sqrtinv231
                dtwoeri( 8,k,j,i,3)= twoeri(13,k,j,i)*sqrtinv385
                dtwoeri( 9,k,j,i,3)= twoeri(14,k,j,i)*sqrtinv231
                dtwoeri(10,k,j,i,3)= twoeri(15,k,j,i)*sqrtinv33
                dtwoeri(11,k,j,i,3)= twoeri(17,k,j,i)*sqrtinv99
                dtwoeri(12,k,j,i,3)= twoeri(18,k,j,i)*sqrtinv231
                dtwoeri(13,k,j,i,3)= twoeri(19,k,j,i)*sqrtinv231
                dtwoeri(14,k,j,i,3)= twoeri(20,k,j,i)*sqrtinv99
                dtwoeri(15,k,j,i,3)= twoeri(21,k,j,i)*sqrtinv11
                dtwoeri(16,k,j,i,3)= twoeri(23,k,j,i)*sqrtinv11
                dtwoeri(17,k,j,i,3)= twoeri(24,k,j,i)*sqrtinv33
                dtwoeri(18,k,j,i,3)= twoeri(25,k,j,i)*sqrt5inv231
                dtwoeri(19,k,j,i,3)= twoeri(26,k,j,i)*sqrtinv33
                dtwoeri(20,k,j,i,3)= twoeri(27,k,j,i)*sqrtinv11
                dtwoeri(21,k,j,i,3)= twoeri(28,k,j,i)
              enddo
            enddo
          enddo
        case default
          write(*,'(" Error! This program supports up to h function in calcd2eri")')
          call iabort
      end select
!
      if(mtype(lsh) >= 1) then
        nangijkl(4)= mtype(lsh)-1
        nbfijkl(4) = ncart(nangijkl(4))
        do l= 1,nprimijkl(4)
          coijkl(l,4)= coeff(lloc+l)
        enddo
!
! Two-electron integral calculation
!
        call int2elec(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxgraddim, &
&                     mxprsh,threshex)
!
        select case(nangijkl(4))
          case(0)
            do i= 1,nbfijkl(1)
              do j= 1,nbfijkl(2)
                do k= 1,nbfijkl(3)
                  dtwoeri(1,k,j,i,1)= dtwoeri(1,k,j,i,1)-twoeri(1,k,j,i)
                  dtwoeri(2,k,j,i,2)= dtwoeri(2,k,j,i,2)-twoeri(1,k,j,i)
                  dtwoeri(3,k,j,i,3)= dtwoeri(3,k,j,i,3)-twoeri(1,k,j,i)
                enddo
              enddo
            enddo
          case(1)
            do i= 1,nbfijkl(1)
              do j= 1,nbfijkl(2)
                do k= 1,nbfijkl(3)
                  dtwoeri(1,k,j,i,1)= dtwoeri(1,k,j,i,1)-twoeri(1,k,j,i)*two
                  dtwoeri(2,k,j,i,1)= dtwoeri(2,k,j,i,1)-twoeri(2,k,j,i)*sqrt3
                  dtwoeri(3,k,j,i,1)= dtwoeri(3,k,j,i,1)-twoeri(3,k,j,i)*sqrt3
                  dtwoeri(2,k,j,i,2)= dtwoeri(2,k,j,i,2)-twoeri(1,k,j,i)*sqrt3
                  dtwoeri(4,k,j,i,2)= dtwoeri(4,k,j,i,2)-twoeri(2,k,j,i)*two
                  dtwoeri(5,k,j,i,2)= dtwoeri(5,k,j,i,2)-twoeri(3,k,j,i)*sqrt3
                  dtwoeri(3,k,j,i,3)= dtwoeri(3,k,j,i,3)-twoeri(1,k,j,i)*sqrt3
                  dtwoeri(5,k,j,i,3)= dtwoeri(5,k,j,i,3)-twoeri(2,k,j,i)*sqrt3
                  dtwoeri(6,k,j,i,3)= dtwoeri(6,k,j,i,3)-twoeri(3,k,j,i)*two
                enddo
              enddo
            enddo
            if(mbf(lsh) == 5) then
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do j= 1,nbfijkl(2)
                    do k= 1,nbfijkl(3)
                      do l= 1,6
                        work(l)= dtwoeri(l,k,j,i,ider)
                      enddo
                      dtwoeri(1,k,j,i,ider)= work(2)
                      dtwoeri(2,k,j,i,ider)= work(5)
                      dtwoeri(3,k,j,i,ider)=(work(6)*two-work(1)-work(4))*half
                      dtwoeri(4,k,j,i,ider)= work(3)
                      dtwoeri(5,k,j,i,ider)=(work(1)-work(4))*sqrt3h
                    enddo
                  enddo
                enddo
              enddo
            endif
          case(2)
            do i= 1,nbfijkl(1)
              do j= 1,nbfijkl(2)
                do k= 1,nbfijkl(3)
                  dtwoeri( 1,k,j,i,1)= dtwoeri( 1,k,j,i,1)-twoeri(1,k,j,i)*three
                  dtwoeri( 2,k,j,i,1)= dtwoeri( 2,k,j,i,1)-twoeri(2,k,j,i)*two*sqrtthird
                  dtwoeri( 3,k,j,i,1)= dtwoeri( 3,k,j,i,1)-twoeri(3,k,j,i)*two*sqrtthird
                  dtwoeri( 4,k,j,i,1)= dtwoeri( 4,k,j,i,1)-twoeri(4,k,j,i)
                  dtwoeri( 5,k,j,i,1)= dtwoeri( 5,k,j,i,1)-twoeri(5,k,j,i)*sqrtthird
                  dtwoeri( 6,k,j,i,1)= dtwoeri( 6,k,j,i,1)-twoeri(6,k,j,i)
                  dtwoeri( 2,k,j,i,2)= dtwoeri( 2,k,j,i,2)-twoeri(1,k,j,i)
                  dtwoeri( 4,k,j,i,2)= dtwoeri( 4,k,j,i,2)-twoeri(2,k,j,i)*two*sqrtthird
                  dtwoeri( 5,k,j,i,2)= dtwoeri( 5,k,j,i,2)-twoeri(3,k,j,i)*sqrtthird
                  dtwoeri( 7,k,j,i,2)= dtwoeri( 7,k,j,i,2)-twoeri(4,k,j,i)*three
                  dtwoeri( 8,k,j,i,2)= dtwoeri( 8,k,j,i,2)-twoeri(5,k,j,i)*two*sqrtthird
                  dtwoeri( 9,k,j,i,2)= dtwoeri( 9,k,j,i,2)-twoeri(6,k,j,i)
                  dtwoeri( 3,k,j,i,3)= dtwoeri( 3,k,j,i,3)-twoeri(1,k,j,i)
                  dtwoeri( 5,k,j,i,3)= dtwoeri( 5,k,j,i,3)-twoeri(2,k,j,i)*sqrtthird
                  dtwoeri( 6,k,j,i,3)= dtwoeri( 6,k,j,i,3)-twoeri(3,k,j,i)*two*sqrtthird
                  dtwoeri( 8,k,j,i,3)= dtwoeri( 8,k,j,i,3)-twoeri(4,k,j,i)
                  dtwoeri( 9,k,j,i,3)= dtwoeri( 9,k,j,i,3)-twoeri(5,k,j,i)*two*sqrtthird
                  dtwoeri(10,k,j,i,3)= dtwoeri(10,k,j,i,3)-twoeri(6,k,j,i)*three
                enddo
              enddo
            enddo
            if(mbf(lsh) == 7) then
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do j= 1,nbfijkl(2)
                    do k= 1,nbfijkl(3)
                      do l= 1,10
                        work(l)= dtwoeri(l,k,j,i,ider)
                      enddo
                      dtwoeri(1,k,j,i,ider)=(-work(7)+work(2)*three                   )*facf1
                      dtwoeri(2,k,j,i,ider)=  work(5)                                  *facf2
                      dtwoeri(3,k,j,i,ider)=(-work(7)-work(2)+work(9)*four            )*facf3
                      dtwoeri(4,k,j,i,ider)=( work(10)*two-work(3)*three-work(8)*three)*half
                      dtwoeri(5,k,j,i,ider)=(-work(1)-work(4)+work(6)*four            )*facf3
                      dtwoeri(6,k,j,i,ider)=( work(3)-work(8)                         )*facf4
                      dtwoeri(7,k,j,i,ider)=( work(1)-work(4)*three                   )*facf1
                    enddo
                  enddo
                enddo
              enddo
            else
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do j= 1,nbfijkl(2)
                    do k= 1,nbfijkl(3)
                      dtwoeri(2,k,j,i,ider)= dtwoeri(2,k,j,i,ider)*sqrt5
                      dtwoeri(3,k,j,i,ider)= dtwoeri(3,k,j,i,ider)*sqrt5
                      dtwoeri(4,k,j,i,ider)= dtwoeri(4,k,j,i,ider)*sqrt5
                      dtwoeri(5,k,j,i,ider)= dtwoeri(5,k,j,i,ider)*sqrt15
                      dtwoeri(6,k,j,i,ider)= dtwoeri(6,k,j,i,ider)*sqrt5
                      dtwoeri(8,k,j,i,ider)= dtwoeri(8,k,j,i,ider)*sqrt5
                      dtwoeri(9,k,j,i,ider)= dtwoeri(9,k,j,i,ider)*sqrt5
                    enddo
                  enddo
                enddo
              enddo
            endif
          case(3)
            do i= 1,nbfijkl(1)
              do j= 1,nbfijkl(2)
                do k= 1,nbfijkl(3)
                  dtwoeri( 1,k,j,i,1)= dtwoeri( 1,k,j,i,1)-twoeri( 1,k,j,i)*four
                  dtwoeri( 2,k,j,i,1)= dtwoeri( 2,k,j,i,1)-twoeri( 2,k,j,i)*three*sqrtfifth
                  dtwoeri( 3,k,j,i,1)= dtwoeri( 3,k,j,i,1)-twoeri( 3,k,j,i)*three*sqrtfifth
                  dtwoeri( 4,k,j,i,1)= dtwoeri( 4,k,j,i,1)-twoeri( 4,k,j,i)*two*sqrtfifth
                  dtwoeri( 5,k,j,i,1)= dtwoeri( 5,k,j,i,1)-twoeri( 5,k,j,i)*two*sqrtinv15
                  dtwoeri( 6,k,j,i,1)= dtwoeri( 6,k,j,i,1)-twoeri( 6,k,j,i)*two*sqrtfifth
                  dtwoeri( 7,k,j,i,1)= dtwoeri( 7,k,j,i,1)-twoeri( 7,k,j,i)
                  dtwoeri( 8,k,j,i,1)= dtwoeri( 8,k,j,i,1)-twoeri( 8,k,j,i)*sqrtfifth
                  dtwoeri( 9,k,j,i,1)= dtwoeri( 9,k,j,i,1)-twoeri( 9,k,j,i)*sqrtfifth
                  dtwoeri(10,k,j,i,1)= dtwoeri(10,k,j,i,1)-twoeri(10,k,j,i)
                  dtwoeri( 2,k,j,i,2)= dtwoeri( 2,k,j,i,2)-twoeri( 1,k,j,i)
                  dtwoeri( 4,k,j,i,2)= dtwoeri( 4,k,j,i,2)-twoeri( 2,k,j,i)*two*sqrtfifth
                  dtwoeri( 5,k,j,i,2)= dtwoeri( 5,k,j,i,2)-twoeri( 3,k,j,i)*sqrtfifth
                  dtwoeri( 7,k,j,i,2)= dtwoeri( 7,k,j,i,2)-twoeri( 4,k,j,i)*three*sqrtfifth
                  dtwoeri( 8,k,j,i,2)= dtwoeri( 8,k,j,i,2)-twoeri( 5,k,j,i)*two*sqrtinv15
                  dtwoeri( 9,k,j,i,2)= dtwoeri( 9,k,j,i,2)-twoeri( 6,k,j,i)*sqrtfifth
                  dtwoeri(11,k,j,i,2)= dtwoeri(11,k,j,i,2)-twoeri( 7,k,j,i)*four
                  dtwoeri(12,k,j,i,2)= dtwoeri(12,k,j,i,2)-twoeri( 8,k,j,i)*three*sqrtfifth
                  dtwoeri(13,k,j,i,2)= dtwoeri(13,k,j,i,2)-twoeri( 9,k,j,i)*two*sqrtfifth
                  dtwoeri(14,k,j,i,2)= dtwoeri(14,k,j,i,2)-twoeri(10,k,j,i)
                  dtwoeri( 3,k,j,i,3)= dtwoeri( 3,k,j,i,3)-twoeri( 1,k,j,i)
                  dtwoeri( 5,k,j,i,3)= dtwoeri( 5,k,j,i,3)-twoeri( 2,k,j,i)*sqrtfifth
                  dtwoeri( 6,k,j,i,3)= dtwoeri( 6,k,j,i,3)-twoeri( 3,k,j,i)*two*sqrtfifth
                  dtwoeri( 8,k,j,i,3)= dtwoeri( 8,k,j,i,3)-twoeri( 4,k,j,i)*sqrtfifth
                  dtwoeri( 9,k,j,i,3)= dtwoeri( 9,k,j,i,3)-twoeri( 5,k,j,i)*two*sqrtinv15
                  dtwoeri(10,k,j,i,3)= dtwoeri(10,k,j,i,3)-twoeri( 6,k,j,i)*three*sqrtfifth
                  dtwoeri(12,k,j,i,3)= dtwoeri(12,k,j,i,3)-twoeri( 7,k,j,i)
                  dtwoeri(13,k,j,i,3)= dtwoeri(13,k,j,i,3)-twoeri( 8,k,j,i)*two*sqrtfifth
                  dtwoeri(14,k,j,i,3)= dtwoeri(14,k,j,i,3)-twoeri( 9,k,j,i)*three*sqrtfifth
                  dtwoeri(15,k,j,i,3)= dtwoeri(15,k,j,i,3)-twoeri(10,k,j,i)*four
                enddo
              enddo
            enddo
            if(mbf(lsh) == 9) then
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do j= 1,nbfijkl(2)
                    do k= 1,nbfijkl(3)
                      do l= 1,15
                        work(l)= dtwoeri(l,k,j,i,ider)
                      enddo
                      dtwoeri(1,k,j,i,ider)=(work(2)-work(7))*facg1
                      dtwoeri(2,k,j,i,ider)=(-work(12)+work(5)*three)*facg2
                      dtwoeri(3,k,j,i,ider)=(-work(2)-work(7)+work(9)*six)*facg3
                      dtwoeri(4,k,j,i,ider)=(-work(12)*three+work(14)*four-work(5)*three)*facg4
                      dtwoeri(5,k,j,i,ider)=(work(1)*three+work(11)*three+work(15)*eight &
&                                           +work(4)*six-work(6)*p24-work(13)*p24)*eighth
                      dtwoeri(6,k,j,i,ider)=(-work(3)*three+work(10)*four-work(8)*three)*facg4
                      dtwoeri(7,k,j,i,ider)=(-work(1)+work(11)+work(6)*six-work(13)*six)*facg5
                      dtwoeri(8,k,j,i,ider)=(work(3)-work(8)*three)*facg2
                      dtwoeri(9,k,j,i,ider)=(work(1)+work(11)-work(4)*six)*facg6
                    enddo
                  enddo
                enddo
              enddo
            else
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do j= 1,nbfijkl(2)
                    do k= 1,nbfijkl(3)
                      dtwoeri( 2,k,j,i,ider)= dtwoeri( 2,k,j,i,ider)*sqrt7
                      dtwoeri( 3,k,j,i,ider)= dtwoeri( 3,k,j,i,ider)*sqrt7
                      dtwoeri( 4,k,j,i,ider)= dtwoeri( 4,k,j,i,ider)*sqrt35third
                      dtwoeri( 5,k,j,i,ider)= dtwoeri( 5,k,j,i,ider)*sqrt35
                      dtwoeri( 6,k,j,i,ider)= dtwoeri( 6,k,j,i,ider)*sqrt35third
                      dtwoeri( 7,k,j,i,ider)= dtwoeri( 7,k,j,i,ider)*sqrt7
                      dtwoeri( 8,k,j,i,ider)= dtwoeri( 8,k,j,i,ider)*sqrt35
                      dtwoeri( 9,k,j,i,ider)= dtwoeri( 9,k,j,i,ider)*sqrt35
                      dtwoeri(10,k,j,i,ider)= dtwoeri(10,k,j,i,ider)*sqrt7
                      dtwoeri(12,k,j,i,ider)= dtwoeri(12,k,j,i,ider)*sqrt7
                      dtwoeri(13,k,j,i,ider)= dtwoeri(13,k,j,i,ider)*sqrt35third
                      dtwoeri(14,k,j,i,ider)= dtwoeri(14,k,j,i,ider)*sqrt7
                    enddo
                  enddo
                enddo
              enddo
            endif
          case(4)
            do i= 1,nbfijkl(1)
              do j= 1,nbfijkl(2)
                do k= 1,nbfijkl(3)
                  dtwoeri( 1,k,j,i,1)= dtwoeri( 1,k,j,i,1)-twoeri( 1,k,j,i)*five
                  dtwoeri( 2,k,j,i,1)= dtwoeri( 2,k,j,i,1)-twoeri( 2,k,j,i)*four*sqrtseventh
                  dtwoeri( 3,k,j,i,1)= dtwoeri( 3,k,j,i,1)-twoeri( 3,k,j,i)*four*sqrtseventh
                  dtwoeri( 4,k,j,i,1)= dtwoeri( 4,k,j,i,1)-twoeri( 4,k,j,i)*three*sqrt3inv35
                  dtwoeri( 5,k,j,i,1)= dtwoeri( 5,k,j,i,1)-twoeri( 5,k,j,i)*three*sqrtinv35
                  dtwoeri( 6,k,j,i,1)= dtwoeri( 6,k,j,i,1)-twoeri( 6,k,j,i)*three*sqrt3inv35
                  dtwoeri( 7,k,j,i,1)= dtwoeri( 7,k,j,i,1)-twoeri( 7,k,j,i)*two*sqrtseventh
                  dtwoeri( 8,k,j,i,1)= dtwoeri( 8,k,j,i,1)-twoeri( 8,k,j,i)*two*sqrtinv35
                  dtwoeri( 9,k,j,i,1)= dtwoeri( 9,k,j,i,1)-twoeri( 9,k,j,i)*two*sqrtinv35
                  dtwoeri(10,k,j,i,1)= dtwoeri(10,k,j,i,1)-twoeri(10,k,j,i)*two*sqrtseventh
                  dtwoeri(11,k,j,i,1)= dtwoeri(11,k,j,i,1)-twoeri(11,k,j,i)
                  dtwoeri(12,k,j,i,1)= dtwoeri(12,k,j,i,1)-twoeri(12,k,j,i)*sqrtseventh
                  dtwoeri(13,k,j,i,1)= dtwoeri(13,k,j,i,1)-twoeri(13,k,j,i)*sqrt3inv35
                  dtwoeri(14,k,j,i,1)= dtwoeri(14,k,j,i,1)-twoeri(14,k,j,i)*sqrtseventh
                  dtwoeri(15,k,j,i,1)= dtwoeri(15,k,j,i,1)-twoeri(15,k,j,i)
                  dtwoeri( 2,k,j,i,2)= dtwoeri( 2,k,j,i,2)-twoeri( 1,k,j,i)
                  dtwoeri( 4,k,j,i,2)= dtwoeri( 4,k,j,i,2)-twoeri( 2,k,j,i)*two*sqrtseventh
                  dtwoeri( 5,k,j,i,2)= dtwoeri( 5,k,j,i,2)-twoeri( 3,k,j,i)*sqrtseventh
                  dtwoeri( 7,k,j,i,2)= dtwoeri( 7,k,j,i,2)-twoeri( 4,k,j,i)*three*sqrt3inv35
                  dtwoeri( 8,k,j,i,2)= dtwoeri( 8,k,j,i,2)-twoeri( 5,k,j,i)*two*sqrtinv35
                  dtwoeri( 9,k,j,i,2)= dtwoeri( 9,k,j,i,2)-twoeri( 6,k,j,i)*sqrt3inv35
                  dtwoeri(11,k,j,i,2)= dtwoeri(11,k,j,i,2)-twoeri( 7,k,j,i)*four*sqrtseventh
                  dtwoeri(12,k,j,i,2)= dtwoeri(12,k,j,i,2)-twoeri( 8,k,j,i)*three*sqrtinv35
                  dtwoeri(13,k,j,i,2)= dtwoeri(13,k,j,i,2)-twoeri( 9,k,j,i)*two*sqrtinv35
                  dtwoeri(14,k,j,i,2)= dtwoeri(14,k,j,i,2)-twoeri(10,k,j,i)*sqrtseventh
                  dtwoeri(16,k,j,i,2)= dtwoeri(16,k,j,i,2)-twoeri(11,k,j,i)*five
                  dtwoeri(17,k,j,i,2)= dtwoeri(17,k,j,i,2)-twoeri(12,k,j,i)*four*sqrtseventh
                  dtwoeri(18,k,j,i,2)= dtwoeri(18,k,j,i,2)-twoeri(13,k,j,i)*three*sqrt3inv35
                  dtwoeri(19,k,j,i,2)= dtwoeri(19,k,j,i,2)-twoeri(14,k,j,i)*two*sqrtseventh
                  dtwoeri(20,k,j,i,2)= dtwoeri(20,k,j,i,2)-twoeri(15,k,j,i)
                  dtwoeri( 3,k,j,i,3)= dtwoeri( 3,k,j,i,3)-twoeri( 1,k,j,i)
                  dtwoeri( 5,k,j,i,3)= dtwoeri( 5,k,j,i,3)-twoeri( 2,k,j,i)*sqrtseventh
                  dtwoeri( 6,k,j,i,3)= dtwoeri( 6,k,j,i,3)-twoeri( 3,k,j,i)*two*sqrtseventh
                  dtwoeri( 8,k,j,i,3)= dtwoeri( 8,k,j,i,3)-twoeri( 4,k,j,i)*sqrt3inv35
                  dtwoeri( 9,k,j,i,3)= dtwoeri( 9,k,j,i,3)-twoeri( 5,k,j,i)*two*sqrtinv35
                  dtwoeri(10,k,j,i,3)= dtwoeri(10,k,j,i,3)-twoeri( 6,k,j,i)*three*sqrt3inv35
                  dtwoeri(12,k,j,i,3)= dtwoeri(12,k,j,i,3)-twoeri( 7,k,j,i)*sqrtseventh
                  dtwoeri(13,k,j,i,3)= dtwoeri(13,k,j,i,3)-twoeri( 8,k,j,i)*two*sqrtinv35
                  dtwoeri(14,k,j,i,3)= dtwoeri(14,k,j,i,3)-twoeri( 9,k,j,i)*three*sqrtinv35
                  dtwoeri(15,k,j,i,3)= dtwoeri(15,k,j,i,3)-twoeri(10,k,j,i)*four*sqrtseventh
                  dtwoeri(17,k,j,i,3)= dtwoeri(17,k,j,i,3)-twoeri(11,k,j,i)
                  dtwoeri(18,k,j,i,3)= dtwoeri(18,k,j,i,3)-twoeri(12,k,j,i)*two*sqrtseventh
                  dtwoeri(19,k,j,i,3)= dtwoeri(19,k,j,i,3)-twoeri(13,k,j,i)*three*sqrt3inv35
                  dtwoeri(20,k,j,i,3)= dtwoeri(20,k,j,i,3)-twoeri(14,k,j,i)*four*sqrtseventh
                  dtwoeri(21,k,j,i,3)= dtwoeri(21,k,j,i,3)-twoeri(15,k,j,i)*five
                enddo
              enddo
            enddo
            if(mbf(lsh) == 11) then
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do j= 1,nbfijkl(2)
                    do k= 1,nbfijkl(3)
                      do l= 1,21
                        work(l)= dtwoeri(l,k,j,i,ider)
                      enddo
                      dtwoeri( 1,k,j,i,ider)=(work(2)*five-work(7)*ten+work(16))*fach1
                      dtwoeri( 2,k,j,i,ider)=(work(5)*four-work(12)*four)*fach2
                      dtwoeri( 3,k,j,i,ider)=(-work(2)*three-work(7)*two+work(9)*p24 &
&                                            +work(16)-work(18)*eight)*fach3
                      dtwoeri( 4,k,j,i,ider)=(-work(5)*two-work(12)*two+work(14)*four)*fach4
                      dtwoeri( 5,k,j,i,ider)=(work(2)+work(7)*two-work(9)*twelve+work(16) &
&                                            -work(18)*twelve+work(20)*eight)*fach5
                      dtwoeri( 6,k,j,i,ider)=(work(3)*p15+work(8)*p30-work(10)*p40 &
&                                            +work(17)*p15-work(19)*p40+work(21)*eight)*eighth
                      dtwoeri( 7,k,j,i,ider)=(work(1)+work(4)*two-work(6)*twelve+work(11) &
&                                            -work(13)*twelve+work(15)*eight)*fach5
                      dtwoeri( 8,k,j,i,ider)=(-work(3)+work(10)*two+work(17)-work(19)*two)*fach4
                      dtwoeri( 9,k,j,i,ider)=(-work(1)+work(4)*two+work(6)*eight &
&                                            +work(11)*three-work(13)*p24)*fach3
                      dtwoeri(10,k,j,i,ider)=(work(3)-work(8)*six+work(17))*fach2
                      dtwoeri(11,k,j,i,ider)=(work(1)-work(4)*ten+work(11)*five)*fach1
                    enddo
                  enddo
                enddo
              enddo
            else
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do j= 1,nbfijkl(2)
                    do k= 1,nbfijkl(3)
                      dtwoeri( 2,k,j,i,ider)= dtwoeri( 2,k,j,i,ider)*three
                      dtwoeri( 3,k,j,i,ider)= dtwoeri( 3,k,j,i,ider)*three
                      dtwoeri( 4,k,j,i,ider)= dtwoeri( 4,k,j,i,ider)*sqrt21
                      dtwoeri( 5,k,j,i,ider)= dtwoeri( 5,k,j,i,ider)*sqrt63
                      dtwoeri( 6,k,j,i,ider)= dtwoeri( 6,k,j,i,ider)*sqrt21
                      dtwoeri( 7,k,j,i,ider)= dtwoeri( 7,k,j,i,ider)*sqrt21
                      dtwoeri( 8,k,j,i,ider)= dtwoeri( 8,k,j,i,ider)*sqrt105
                      dtwoeri( 9,k,j,i,ider)= dtwoeri( 9,k,j,i,ider)*sqrt105
                      dtwoeri(10,k,j,i,ider)= dtwoeri(10,k,j,i,ider)*sqrt21
                      dtwoeri(11,k,j,i,ider)= dtwoeri(11,k,j,i,ider)*three
                      dtwoeri(12,k,j,i,ider)= dtwoeri(12,k,j,i,ider)*sqrt63
                      dtwoeri(13,k,j,i,ider)= dtwoeri(13,k,j,i,ider)*sqrt105
                      dtwoeri(14,k,j,i,ider)= dtwoeri(14,k,j,i,ider)*sqrt63
                      dtwoeri(15,k,j,i,ider)= dtwoeri(15,k,j,i,ider)*three
                      dtwoeri(17,k,j,i,ider)= dtwoeri(17,k,j,i,ider)*three
                      dtwoeri(18,k,j,i,ider)= dtwoeri(18,k,j,i,ider)*sqrt21
                      dtwoeri(19,k,j,i,ider)= dtwoeri(19,k,j,i,ider)*sqrt21
                      dtwoeri(20,k,j,i,ider)= dtwoeri(20,k,j,i,ider)*three
                    enddo
                  enddo
                enddo
              enddo
            endif
        end select
      endif
!
      nangijkl(4)= mtype(lsh)
      nbfijkl(4) = mbf(lsh)
      do l= 1,nprimijkl(4)
        coijkl(l,4)= coeff(lloc+l)
      enddo
      do i= 1,nbfijkl(1)
        do j= 1,nbfijkl(2)
          do k= 1,nbfijkl(3)
            do l= 1,nbfijkl(4)
              gradtwo(1,4)= gradtwo(1,4)+pdmtrx(l,k,j,i)*dtwoeri(l,k,j,i,1)
              gradtwo(2,4)= gradtwo(2,4)+pdmtrx(l,k,j,i)*dtwoeri(l,k,j,i,2)
              gradtwo(3,4)= gradtwo(3,4)+pdmtrx(l,k,j,i)*dtwoeri(l,k,j,i,3)
            enddo
          enddo
        enddo
      enddo
!
! Ksh derivative
!
      nangijkl(3)= mtype(ksh)+1
      nbfijkl(3) = ncart(nangijkl(3))
      do k= 1,nprimijkl(3)
        coijkl(k,3)= two*ex(kloc+k)*coeff(kloc+k)
      enddo
!
! Two-electron integral calculation
!
      call int2elec(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxgraddim, &
&                   mxprsh,threshex)
!
      select case(nangijkl(3))
        case(1)
          do i= 1,nbfijkl(1)
            do j= 1,nbfijkl(2)
              do l= 1,nbfijkl(4)
                dtwoeri(l,1,j,i,1)= twoeri(l,1,j,i)
                dtwoeri(l,1,j,i,2)= twoeri(l,2,j,i)
                dtwoeri(l,1,j,i,3)= twoeri(l,3,j,i)
              enddo
            enddo
          enddo
        case(2)
          do i= 1,nbfijkl(1)
            do j= 1,nbfijkl(2)
              do l= 1,nbfijkl(4)
                dtwoeri(l,1,j,i,1)= twoeri(l,1,j,i)
                dtwoeri(l,2,j,i,1)= twoeri(l,2,j,i)*sqrtthird
                dtwoeri(l,3,j,i,1)= twoeri(l,3,j,i)*sqrtthird
                dtwoeri(l,1,j,i,2)= twoeri(l,2,j,i)*sqrtthird
                dtwoeri(l,2,j,i,2)= twoeri(l,4,j,i)
                dtwoeri(l,3,j,i,2)= twoeri(l,5,j,i)*sqrtthird
                dtwoeri(l,1,j,i,3)= twoeri(l,3,j,i)*sqrtthird
                dtwoeri(l,2,j,i,3)= twoeri(l,5,j,i)*sqrtthird
                dtwoeri(l,3,j,i,3)= twoeri(l,6,j,i)
              enddo
            enddo
          enddo
        case(3)
          do i= 1,nbfijkl(1)
            do j= 1,nbfijkl(2)
              do l= 1,nbfijkl(4)
                dtwoeri(l,1,j,i,1)= twoeri(l, 1,j,i)
                dtwoeri(l,2,j,i,1)= twoeri(l, 2,j,i)*sqrt3fifth
                dtwoeri(l,3,j,i,1)= twoeri(l, 3,j,i)*sqrt3fifth
                dtwoeri(l,4,j,i,1)= twoeri(l, 4,j,i)*sqrtfifth
                dtwoeri(l,5,j,i,1)= twoeri(l, 5,j,i)*sqrtfifth
                dtwoeri(l,6,j,i,1)= twoeri(l, 6,j,i)*sqrtfifth
                dtwoeri(l,1,j,i,2)= twoeri(l, 2,j,i)*sqrtfifth
                dtwoeri(l,2,j,i,2)= twoeri(l, 4,j,i)*sqrt3fifth
                dtwoeri(l,3,j,i,2)= twoeri(l, 5,j,i)*sqrtfifth
                dtwoeri(l,4,j,i,2)= twoeri(l, 7,j,i)
                dtwoeri(l,5,j,i,2)= twoeri(l, 8,j,i)*sqrt3fifth
                dtwoeri(l,6,j,i,2)= twoeri(l, 9,j,i)*sqrtfifth
                dtwoeri(l,1,j,i,3)= twoeri(l, 3,j,i)*sqrtfifth
                dtwoeri(l,2,j,i,3)= twoeri(l, 5,j,i)*sqrtfifth
                dtwoeri(l,3,j,i,3)= twoeri(l, 6,j,i)*sqrt3fifth
                dtwoeri(l,4,j,i,3)= twoeri(l, 8,j,i)*sqrtfifth
                dtwoeri(l,5,j,i,3)= twoeri(l, 9,j,i)*sqrt3fifth
                dtwoeri(l,6,j,i,3)= twoeri(l,10,j,i)
              enddo
            enddo
          enddo
        case(4)
          do i= 1,nbfijkl(1)
            do j= 1,nbfijkl(2)
              do l= 1,nbfijkl(4)
                dtwoeri(l, 1,j,i,1)= twoeri(l, 1,j,i)
                dtwoeri(l, 2,j,i,1)= twoeri(l, 2,j,i)*sqrtseventh
                dtwoeri(l, 3,j,i,1)= twoeri(l, 3,j,i)*sqrtseventh
                dtwoeri(l, 4,j,i,1)= twoeri(l, 4,j,i)*sqrt3inv35
                dtwoeri(l, 5,j,i,1)= twoeri(l, 5,j,i)*sqrtinv35
                dtwoeri(l, 6,j,i,1)= twoeri(l, 6,j,i)*sqrt3inv35
                dtwoeri(l, 7,j,i,1)= twoeri(l, 7,j,i)*sqrtseventh
                dtwoeri(l, 8,j,i,1)= twoeri(l, 8,j,i)*sqrtinv35
                dtwoeri(l, 9,j,i,1)= twoeri(l, 9,j,i)*sqrtinv35
                dtwoeri(l,10,j,i,1)= twoeri(l,10,j,i)*sqrtseventh
                dtwoeri(l, 1,j,i,2)= twoeri(l, 2,j,i)*sqrtseventh
                dtwoeri(l, 2,j,i,2)= twoeri(l, 4,j,i)*sqrt3inv35
                dtwoeri(l, 3,j,i,2)= twoeri(l, 5,j,i)*sqrtinv35
                dtwoeri(l, 4,j,i,2)= twoeri(l, 7,j,i)*sqrtseventh
                dtwoeri(l, 5,j,i,2)= twoeri(l, 8,j,i)*sqrtinv35
                dtwoeri(l, 6,j,i,2)= twoeri(l, 9,j,i)*sqrtinv35
                dtwoeri(l, 7,j,i,2)= twoeri(l,11,j,i)
                dtwoeri(l, 8,j,i,2)= twoeri(l,12,j,i)*sqrtseventh
                dtwoeri(l, 9,j,i,2)= twoeri(l,13,j,i)*sqrt3inv35
                dtwoeri(l,10,j,i,2)= twoeri(l,14,j,i)*sqrtseventh
                dtwoeri(l, 1,j,i,3)= twoeri(l, 3,j,i)*sqrtseventh
                dtwoeri(l, 2,j,i,3)= twoeri(l, 5,j,i)*sqrtinv35
                dtwoeri(l, 3,j,i,3)= twoeri(l, 6,j,i)*sqrt3inv35
                dtwoeri(l, 4,j,i,3)= twoeri(l, 8,j,i)*sqrtinv35
                dtwoeri(l, 5,j,i,3)= twoeri(l, 9,j,i)*sqrtinv35
                dtwoeri(l, 6,j,i,3)= twoeri(l,10,j,i)*sqrtseventh
                dtwoeri(l, 7,j,i,3)= twoeri(l,12,j,i)*sqrtseventh
                dtwoeri(l, 8,j,i,3)= twoeri(l,13,j,i)*sqrt3inv35
                dtwoeri(l, 9,j,i,3)= twoeri(l,14,j,i)*sqrtseventh
                dtwoeri(l,10,j,i,3)= twoeri(l,15,j,i)
              enddo
            enddo
          enddo
        case(5)
          do i= 1,nbfijkl(1)
            do j= 1,nbfijkl(2)
              do l= 1,nbfijkl(4)
                dtwoeri(l, 1,j,i,1)= twoeri(l, 1,j,i)
                dtwoeri(l, 2,j,i,1)= twoeri(l, 2,j,i)*third
                dtwoeri(l, 3,j,i,1)= twoeri(l, 3,j,i)*third
                dtwoeri(l, 4,j,i,1)= twoeri(l, 4,j,i)*sqrtinv21
                dtwoeri(l, 5,j,i,1)= twoeri(l, 5,j,i)*sqrtinv63
                dtwoeri(l, 6,j,i,1)= twoeri(l, 6,j,i)*sqrtinv21
                dtwoeri(l, 7,j,i,1)= twoeri(l, 7,j,i)*sqrtinv21
                dtwoeri(l, 8,j,i,1)= twoeri(l, 8,j,i)*sqrtinv105
                dtwoeri(l, 9,j,i,1)= twoeri(l, 9,j,i)*sqrtinv105
                dtwoeri(l,10,j,i,1)= twoeri(l,10,j,i)*sqrtinv21
                dtwoeri(l,11,j,i,1)= twoeri(l,11,j,i)*third
                dtwoeri(l,12,j,i,1)= twoeri(l,12,j,i)*sqrtinv63
                dtwoeri(l,13,j,i,1)= twoeri(l,13,j,i)*sqrtinv105
                dtwoeri(l,14,j,i,1)= twoeri(l,14,j,i)*sqrtinv63
                dtwoeri(l,15,j,i,1)= twoeri(l,15,j,i)*third
                dtwoeri(l, 1,j,i,2)= twoeri(l, 2,j,i)*third
                dtwoeri(l, 2,j,i,2)= twoeri(l, 4,j,i)*sqrtinv21
                dtwoeri(l, 3,j,i,2)= twoeri(l, 5,j,i)*sqrtinv63
                dtwoeri(l, 4,j,i,2)= twoeri(l, 7,j,i)*sqrtinv21
                dtwoeri(l, 5,j,i,2)= twoeri(l, 8,j,i)*sqrtinv105
                dtwoeri(l, 6,j,i,2)= twoeri(l, 9,j,i)*sqrtinv105
                dtwoeri(l, 7,j,i,2)= twoeri(l,11,j,i)*third
                dtwoeri(l, 8,j,i,2)= twoeri(l,12,j,i)*sqrtinv63
                dtwoeri(l, 9,j,i,2)= twoeri(l,13,j,i)*sqrtinv105
                dtwoeri(l,10,j,i,2)= twoeri(l,14,j,i)*sqrtinv63
                dtwoeri(l,11,j,i,2)= twoeri(l,16,j,i)
                dtwoeri(l,12,j,i,2)= twoeri(l,17,j,i)*third
                dtwoeri(l,13,j,i,2)= twoeri(l,18,j,i)*sqrtinv21
                dtwoeri(l,14,j,i,2)= twoeri(l,19,j,i)*sqrtinv21
                dtwoeri(l,15,j,i,2)= twoeri(l,20,j,i)*third
                dtwoeri(l, 1,j,i,3)= twoeri(l, 3,j,i)*third
                dtwoeri(l, 2,j,i,3)= twoeri(l, 5,j,i)*sqrtinv63
                dtwoeri(l, 3,j,i,3)= twoeri(l, 6,j,i)*sqrtinv21
                dtwoeri(l, 4,j,i,3)= twoeri(l, 8,j,i)*sqrtinv105
                dtwoeri(l, 5,j,i,3)= twoeri(l, 9,j,i)*sqrtinv105
                dtwoeri(l, 6,j,i,3)= twoeri(l,10,j,i)*sqrtinv21
                dtwoeri(l, 7,j,i,3)= twoeri(l,12,j,i)*sqrtinv63
                dtwoeri(l, 8,j,i,3)= twoeri(l,13,j,i)*sqrtinv105
                dtwoeri(l, 9,j,i,3)= twoeri(l,14,j,i)*sqrtinv63
                dtwoeri(l,10,j,i,3)= twoeri(l,15,j,i)*third
                dtwoeri(l,11,j,i,3)= twoeri(l,17,j,i)*third
                dtwoeri(l,12,j,i,3)= twoeri(l,18,j,i)*sqrtinv21
                dtwoeri(l,13,j,i,3)= twoeri(l,19,j,i)*sqrtinv21
                dtwoeri(l,14,j,i,3)= twoeri(l,20,j,i)*third
                dtwoeri(l,15,j,i,3)= twoeri(l,21,j,i)
              enddo
            enddo
          enddo
        case(6)
          do i= 1,nbfijkl(1)
            do j= 1,nbfijkl(2)
              do l= 1,nbfijkl(4)
                dtwoeri(l, 1,j,i,1)= twoeri(l, 1,j,i)
                dtwoeri(l, 2,j,i,1)= twoeri(l, 2,j,i)*sqrtinv11
                dtwoeri(l, 3,j,i,1)= twoeri(l, 3,j,i)*sqrtinv11
                dtwoeri(l, 4,j,i,1)= twoeri(l, 4,j,i)*sqrtinv33
                dtwoeri(l, 5,j,i,1)= twoeri(l, 5,j,i)*sqrtinv99
                dtwoeri(l, 6,j,i,1)= twoeri(l, 6,j,i)*sqrtinv33
                dtwoeri(l, 7,j,i,1)= twoeri(l, 7,j,i)*sqrt5inv231
                dtwoeri(l, 8,j,i,1)= twoeri(l, 8,j,i)*sqrtinv231
                dtwoeri(l, 9,j,i,1)= twoeri(l, 9,j,i)*sqrtinv231
                dtwoeri(l,10,j,i,1)= twoeri(l,10,j,i)*sqrt5inv231
                dtwoeri(l,11,j,i,1)= twoeri(l,11,j,i)*sqrtinv33
                dtwoeri(l,12,j,i,1)= twoeri(l,12,j,i)*sqrtinv231
                dtwoeri(l,13,j,i,1)= twoeri(l,13,j,i)*sqrtinv385
                dtwoeri(l,14,j,i,1)= twoeri(l,14,j,i)*sqrtinv231
                dtwoeri(l,15,j,i,1)= twoeri(l,15,j,i)*sqrtinv33
                dtwoeri(l,16,j,i,1)= twoeri(l,16,j,i)*sqrtinv11
                dtwoeri(l,17,j,i,1)= twoeri(l,17,j,i)*sqrtinv99
                dtwoeri(l,18,j,i,1)= twoeri(l,18,j,i)*sqrtinv231
                dtwoeri(l,19,j,i,1)= twoeri(l,19,j,i)*sqrtinv231
                dtwoeri(l,20,j,i,1)= twoeri(l,20,j,i)*sqrtinv99
                dtwoeri(l,21,j,i,1)= twoeri(l,21,j,i)*sqrtinv11
                dtwoeri(l, 1,j,i,2)= twoeri(l, 2,j,i)*sqrtinv11
                dtwoeri(l, 2,j,i,2)= twoeri(l, 4,j,i)*sqrtinv33
                dtwoeri(l, 3,j,i,2)= twoeri(l, 5,j,i)*sqrtinv99
                dtwoeri(l, 4,j,i,2)= twoeri(l, 7,j,i)*sqrt5inv231
                dtwoeri(l, 5,j,i,2)= twoeri(l, 8,j,i)*sqrtinv231
                dtwoeri(l, 6,j,i,2)= twoeri(l, 9,j,i)*sqrtinv231
                dtwoeri(l, 7,j,i,2)= twoeri(l,11,j,i)*sqrtinv33
                dtwoeri(l, 8,j,i,2)= twoeri(l,12,j,i)*sqrtinv231
                dtwoeri(l, 9,j,i,2)= twoeri(l,13,j,i)*sqrtinv385
                dtwoeri(l,10,j,i,2)= twoeri(l,14,j,i)*sqrtinv231
                dtwoeri(l,11,j,i,2)= twoeri(l,16,j,i)*sqrtinv11
                dtwoeri(l,12,j,i,2)= twoeri(l,17,j,i)*sqrtinv99
                dtwoeri(l,13,j,i,2)= twoeri(l,18,j,i)*sqrtinv231
                dtwoeri(l,14,j,i,2)= twoeri(l,19,j,i)*sqrtinv231
                dtwoeri(l,15,j,i,2)= twoeri(l,20,j,i)*sqrtinv99
                dtwoeri(l,16,j,i,2)= twoeri(l,22,j,i)
                dtwoeri(l,17,j,i,2)= twoeri(l,23,j,i)*sqrtinv11
                dtwoeri(l,18,j,i,2)= twoeri(l,24,j,i)*sqrtinv33
                dtwoeri(l,19,j,i,2)= twoeri(l,25,j,i)*sqrt5inv231
                dtwoeri(l,20,j,i,2)= twoeri(l,26,j,i)*sqrtinv33
                dtwoeri(l,21,j,i,2)= twoeri(l,27,j,i)*sqrtinv11
                dtwoeri(l, 1,j,i,3)= twoeri(l, 3,j,i)*sqrtinv11
                dtwoeri(l, 2,j,i,3)= twoeri(l, 5,j,i)*sqrtinv99
                dtwoeri(l, 3,j,i,3)= twoeri(l, 6,j,i)*sqrtinv33
                dtwoeri(l, 4,j,i,3)= twoeri(l, 8,j,i)*sqrtinv231
                dtwoeri(l, 5,j,i,3)= twoeri(l, 9,j,i)*sqrtinv231
                dtwoeri(l, 6,j,i,3)= twoeri(l,10,j,i)*sqrt5inv231
                dtwoeri(l, 7,j,i,3)= twoeri(l,12,j,i)*sqrtinv231
                dtwoeri(l, 8,j,i,3)= twoeri(l,13,j,i)*sqrtinv385
                dtwoeri(l, 9,j,i,3)= twoeri(l,14,j,i)*sqrtinv231
                dtwoeri(l,10,j,i,3)= twoeri(l,15,j,i)*sqrtinv33
                dtwoeri(l,11,j,i,3)= twoeri(l,17,j,i)*sqrtinv99
                dtwoeri(l,12,j,i,3)= twoeri(l,18,j,i)*sqrtinv231
                dtwoeri(l,13,j,i,3)= twoeri(l,19,j,i)*sqrtinv231
                dtwoeri(l,14,j,i,3)= twoeri(l,20,j,i)*sqrtinv99
                dtwoeri(l,15,j,i,3)= twoeri(l,21,j,i)*sqrtinv11
                dtwoeri(l,16,j,i,3)= twoeri(l,23,j,i)*sqrtinv11
                dtwoeri(l,17,j,i,3)= twoeri(l,24,j,i)*sqrtinv33
                dtwoeri(l,18,j,i,3)= twoeri(l,25,j,i)*sqrt5inv231
                dtwoeri(l,19,j,i,3)= twoeri(l,26,j,i)*sqrtinv33
                dtwoeri(l,20,j,i,3)= twoeri(l,27,j,i)*sqrtinv11
                dtwoeri(l,21,j,i,3)= twoeri(l,28,j,i)
              enddo
            enddo
          enddo
        case default
          write(*,'(" Error! This program supports up to h function in calcd2eri")')
          call iabort
      end select
!
      if(mtype(ksh) >= 1) then
        nangijkl(3)= mtype(ksh)-1
        nbfijkl(3) = ncart(nangijkl(3))
        do k= 1,nprimijkl(3)
          coijkl(k,3)= coeff(kloc+k)
        enddo
!
! Two-electron integral calculation
!
        call int2elec(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxgraddim, &
&                     mxprsh,threshex)
!
        select case(nangijkl(3))
          case(0)
            do i= 1,nbfijkl(1)
              do j= 1,nbfijkl(2)
                do l= 1,nbfijkl(4)
                  dtwoeri(l,1,j,i,1)= dtwoeri(l,1,j,i,1)-twoeri(l,1,j,i)
                  dtwoeri(l,2,j,i,2)= dtwoeri(l,2,j,i,2)-twoeri(l,1,j,i)
                  dtwoeri(l,3,j,i,3)= dtwoeri(l,3,j,i,3)-twoeri(l,1,j,i)
                enddo
              enddo
            enddo
          case(1)
            do i= 1,nbfijkl(1)
              do j= 1,nbfijkl(2)
                do l= 1,nbfijkl(4)
                  dtwoeri(l,1,j,i,1)= dtwoeri(l,1,j,i,1)-twoeri(l,1,j,i)*two
                  dtwoeri(l,2,j,i,1)= dtwoeri(l,2,j,i,1)-twoeri(l,2,j,i)*sqrt3
                  dtwoeri(l,3,j,i,1)= dtwoeri(l,3,j,i,1)-twoeri(l,3,j,i)*sqrt3
                  dtwoeri(l,2,j,i,2)= dtwoeri(l,2,j,i,2)-twoeri(l,1,j,i)*sqrt3
                  dtwoeri(l,4,j,i,2)= dtwoeri(l,4,j,i,2)-twoeri(l,2,j,i)*two
                  dtwoeri(l,5,j,i,2)= dtwoeri(l,5,j,i,2)-twoeri(l,3,j,i)*sqrt3
                  dtwoeri(l,3,j,i,3)= dtwoeri(l,3,j,i,3)-twoeri(l,1,j,i)*sqrt3
                  dtwoeri(l,5,j,i,3)= dtwoeri(l,5,j,i,3)-twoeri(l,2,j,i)*sqrt3
                  dtwoeri(l,6,j,i,3)= dtwoeri(l,6,j,i,3)-twoeri(l,3,j,i)*two
                enddo
              enddo
            enddo
            if(mbf(ksh) == 5) then
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do j= 1,nbfijkl(2)
                    do l= 1,nbfijkl(4)
                      do k= 1,6
                        work(k)= dtwoeri(l,k,j,i,ider)
                      enddo
                      dtwoeri(l,1,j,i,ider)= work(2)
                      dtwoeri(l,2,j,i,ider)= work(5)
                      dtwoeri(l,3,j,i,ider)=(work(6)*two-work(1)-work(4))*half
                      dtwoeri(l,4,j,i,ider)= work(3)
                      dtwoeri(l,5,j,i,ider)=(work(1)-work(4))*sqrt3h
                    enddo
                  enddo
                enddo
              enddo
            endif
          case(2)
            do i= 1,nbfijkl(1)
              do j= 1,nbfijkl(2)
                do l= 1,nbfijkl(4)
                  dtwoeri(l, 1,j,i,1)= dtwoeri(l, 1,j,i,1)-twoeri(l,1,j,i)*three
                  dtwoeri(l, 2,j,i,1)= dtwoeri(l, 2,j,i,1)-twoeri(l,2,j,i)*two*sqrtthird
                  dtwoeri(l, 3,j,i,1)= dtwoeri(l, 3,j,i,1)-twoeri(l,3,j,i)*two*sqrtthird
                  dtwoeri(l, 4,j,i,1)= dtwoeri(l, 4,j,i,1)-twoeri(l,4,j,i)
                  dtwoeri(l, 5,j,i,1)= dtwoeri(l, 5,j,i,1)-twoeri(l,5,j,i)*sqrtthird
                  dtwoeri(l, 6,j,i,1)= dtwoeri(l, 6,j,i,1)-twoeri(l,6,j,i)
                  dtwoeri(l, 2,j,i,2)= dtwoeri(l, 2,j,i,2)-twoeri(l,1,j,i)
                  dtwoeri(l, 4,j,i,2)= dtwoeri(l, 4,j,i,2)-twoeri(l,2,j,i)*two*sqrtthird
                  dtwoeri(l, 5,j,i,2)= dtwoeri(l, 5,j,i,2)-twoeri(l,3,j,i)*sqrtthird
                  dtwoeri(l, 7,j,i,2)= dtwoeri(l, 7,j,i,2)-twoeri(l,4,j,i)*three
                  dtwoeri(l, 8,j,i,2)= dtwoeri(l, 8,j,i,2)-twoeri(l,5,j,i)*two*sqrtthird
                  dtwoeri(l, 9,j,i,2)= dtwoeri(l, 9,j,i,2)-twoeri(l,6,j,i)
                  dtwoeri(l, 3,j,i,3)= dtwoeri(l, 3,j,i,3)-twoeri(l,1,j,i)
                  dtwoeri(l, 5,j,i,3)= dtwoeri(l, 5,j,i,3)-twoeri(l,2,j,i)*sqrtthird
                  dtwoeri(l, 6,j,i,3)= dtwoeri(l, 6,j,i,3)-twoeri(l,3,j,i)*two*sqrtthird
                  dtwoeri(l, 8,j,i,3)= dtwoeri(l, 8,j,i,3)-twoeri(l,4,j,i)
                  dtwoeri(l, 9,j,i,3)= dtwoeri(l, 9,j,i,3)-twoeri(l,5,j,i)*two*sqrtthird
                  dtwoeri(l,10,j,i,3)= dtwoeri(l,10,j,i,3)-twoeri(l,6,j,i)*three
                enddo
              enddo
            enddo
            if(mbf(ksh) == 7) then
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do j= 1,nbfijkl(2)
                    do l= 1,nbfijkl(4)
                      do k= 1,10
                        work(k)= dtwoeri(l,k,j,i,ider)
                      enddo
                      dtwoeri(l,1,j,i,ider)=(-work(7)+work(2)*three                   )*facf1
                      dtwoeri(l,2,j,i,ider)=  work(5)                                  *facf2
                      dtwoeri(l,3,j,i,ider)=(-work(7)-work(2)+work(9)*four            )*facf3
                      dtwoeri(l,4,j,i,ider)=( work(10)*two-work(3)*three-work(8)*three)*half
                      dtwoeri(l,5,j,i,ider)=(-work(1)-work(4)+work(6)*four            )*facf3
                      dtwoeri(l,6,j,i,ider)=( work(3)-work(8)                         )*facf4
                      dtwoeri(l,7,j,i,ider)=( work(1)-work(4)*three                   )*facf1
                    enddo
                  enddo
                enddo
              enddo
            else
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do j= 1,nbfijkl(2)
                    do l= 1,nbfijkl(4)
                      dtwoeri(l,2,j,i,ider)= dtwoeri(l,2,j,i,ider)*sqrt5
                      dtwoeri(l,3,j,i,ider)= dtwoeri(l,3,j,i,ider)*sqrt5
                      dtwoeri(l,4,j,i,ider)= dtwoeri(l,4,j,i,ider)*sqrt5
                      dtwoeri(l,5,j,i,ider)= dtwoeri(l,5,j,i,ider)*sqrt15
                      dtwoeri(l,6,j,i,ider)= dtwoeri(l,6,j,i,ider)*sqrt5
                      dtwoeri(l,8,j,i,ider)= dtwoeri(l,8,j,i,ider)*sqrt5
                      dtwoeri(l,9,j,i,ider)= dtwoeri(l,9,j,i,ider)*sqrt5
                    enddo
                  enddo
                enddo
              enddo
            endif
          case(3)
            do i= 1,nbfijkl(1)
              do j= 1,nbfijkl(2)
                do l= 1,nbfijkl(4)
                  dtwoeri(l, 1,j,i,1)= dtwoeri(l, 1,j,i,1)-twoeri(l, 1,j,i)*four
                  dtwoeri(l, 2,j,i,1)= dtwoeri(l, 2,j,i,1)-twoeri(l, 2,j,i)*three*sqrtfifth
                  dtwoeri(l, 3,j,i,1)= dtwoeri(l, 3,j,i,1)-twoeri(l, 3,j,i)*three*sqrtfifth
                  dtwoeri(l, 4,j,i,1)= dtwoeri(l, 4,j,i,1)-twoeri(l, 4,j,i)*two*sqrtfifth
                  dtwoeri(l, 5,j,i,1)= dtwoeri(l, 5,j,i,1)-twoeri(l, 5,j,i)*two*sqrtinv15
                  dtwoeri(l, 6,j,i,1)= dtwoeri(l, 6,j,i,1)-twoeri(l, 6,j,i)*two*sqrtfifth
                  dtwoeri(l, 7,j,i,1)= dtwoeri(l, 7,j,i,1)-twoeri(l, 7,j,i)
                  dtwoeri(l, 8,j,i,1)= dtwoeri(l, 8,j,i,1)-twoeri(l, 8,j,i)*sqrtfifth
                  dtwoeri(l, 9,j,i,1)= dtwoeri(l, 9,j,i,1)-twoeri(l, 9,j,i)*sqrtfifth
                  dtwoeri(l,10,j,i,1)= dtwoeri(l,10,j,i,1)-twoeri(l,10,j,i)
                  dtwoeri(l, 2,j,i,2)= dtwoeri(l, 2,j,i,2)-twoeri(l, 1,j,i)
                  dtwoeri(l, 4,j,i,2)= dtwoeri(l, 4,j,i,2)-twoeri(l, 2,j,i)*two*sqrtfifth
                  dtwoeri(l, 5,j,i,2)= dtwoeri(l, 5,j,i,2)-twoeri(l, 3,j,i)*sqrtfifth
                  dtwoeri(l, 7,j,i,2)= dtwoeri(l, 7,j,i,2)-twoeri(l, 4,j,i)*three*sqrtfifth
                  dtwoeri(l, 8,j,i,2)= dtwoeri(l, 8,j,i,2)-twoeri(l, 5,j,i)*two*sqrtinv15
                  dtwoeri(l, 9,j,i,2)= dtwoeri(l, 9,j,i,2)-twoeri(l, 6,j,i)*sqrtfifth
                  dtwoeri(l,11,j,i,2)= dtwoeri(l,11,j,i,2)-twoeri(l, 7,j,i)*four
                  dtwoeri(l,12,j,i,2)= dtwoeri(l,12,j,i,2)-twoeri(l, 8,j,i)*three*sqrtfifth
                  dtwoeri(l,13,j,i,2)= dtwoeri(l,13,j,i,2)-twoeri(l, 9,j,i)*two*sqrtfifth
                  dtwoeri(l,14,j,i,2)= dtwoeri(l,14,j,i,2)-twoeri(l,10,j,i)
                  dtwoeri(l, 3,j,i,3)= dtwoeri(l, 3,j,i,3)-twoeri(l, 1,j,i)
                  dtwoeri(l, 5,j,i,3)= dtwoeri(l, 5,j,i,3)-twoeri(l, 2,j,i)*sqrtfifth
                  dtwoeri(l, 6,j,i,3)= dtwoeri(l, 6,j,i,3)-twoeri(l, 3,j,i)*two*sqrtfifth
                  dtwoeri(l, 8,j,i,3)= dtwoeri(l, 8,j,i,3)-twoeri(l, 4,j,i)*sqrtfifth
                  dtwoeri(l, 9,j,i,3)= dtwoeri(l, 9,j,i,3)-twoeri(l, 5,j,i)*two*sqrtinv15
                  dtwoeri(l,10,j,i,3)= dtwoeri(l,10,j,i,3)-twoeri(l, 6,j,i)*three*sqrtfifth
                  dtwoeri(l,12,j,i,3)= dtwoeri(l,12,j,i,3)-twoeri(l, 7,j,i)
                  dtwoeri(l,13,j,i,3)= dtwoeri(l,13,j,i,3)-twoeri(l, 8,j,i)*two*sqrtfifth
                  dtwoeri(l,14,j,i,3)= dtwoeri(l,14,j,i,3)-twoeri(l, 9,j,i)*three*sqrtfifth
                  dtwoeri(l,15,j,i,3)= dtwoeri(l,15,j,i,3)-twoeri(l,10,j,i)*four
                enddo
              enddo
            enddo
            if(mbf(ksh) == 9) then
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do j= 1,nbfijkl(2)
                    do l= 1,nbfijkl(4)
                      do k= 1,15
                        work(k)= dtwoeri(l,k,j,i,ider)
                      enddo
                      dtwoeri(l,1,j,i,ider)=(work(2)-work(7))*facg1
                      dtwoeri(l,2,j,i,ider)=(-work(12)+work(5)*three)*facg2
                      dtwoeri(l,3,j,i,ider)=(-work(2)-work(7)+work(9)*six)*facg3
                      dtwoeri(l,4,j,i,ider)=(-work(12)*three+work(14)*four-work(5)*three)*facg4
                      dtwoeri(l,5,j,i,ider)=(work(1)*three+work(11)*three+work(15)*eight &
&                                           +work(4)*six-work(6)*p24-work(13)*p24)*eighth
                      dtwoeri(l,6,j,i,ider)=(-work(3)*three+work(10)*four-work(8)*three)*facg4
                      dtwoeri(l,7,j,i,ider)=(-work(1)+work(11)+work(6)*six-work(13)*six)*facg5
                      dtwoeri(l,8,j,i,ider)=(work(3)-work(8)*three)*facg2
                      dtwoeri(l,9,j,i,ider)=(work(1)+work(11)-work(4)*six)*facg6
                    enddo
                  enddo
                enddo
              enddo
            else
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do j= 1,nbfijkl(2)
                    do l= 1,nbfijkl(4)
                      dtwoeri(l, 2,j,i,ider)= dtwoeri(l, 2,j,i,ider)*sqrt7
                      dtwoeri(l, 3,j,i,ider)= dtwoeri(l, 3,j,i,ider)*sqrt7
                      dtwoeri(l, 4,j,i,ider)= dtwoeri(l, 4,j,i,ider)*sqrt35third
                      dtwoeri(l, 5,j,i,ider)= dtwoeri(l, 5,j,i,ider)*sqrt35
                      dtwoeri(l, 6,j,i,ider)= dtwoeri(l, 6,j,i,ider)*sqrt35third
                      dtwoeri(l, 7,j,i,ider)= dtwoeri(l, 7,j,i,ider)*sqrt7
                      dtwoeri(l, 8,j,i,ider)= dtwoeri(l, 8,j,i,ider)*sqrt35
                      dtwoeri(l, 9,j,i,ider)= dtwoeri(l, 9,j,i,ider)*sqrt35
                      dtwoeri(l,10,j,i,ider)= dtwoeri(l,10,j,i,ider)*sqrt7
                      dtwoeri(l,12,j,i,ider)= dtwoeri(l,12,j,i,ider)*sqrt7
                      dtwoeri(l,13,j,i,ider)= dtwoeri(l,13,j,i,ider)*sqrt35third
                      dtwoeri(l,14,j,i,ider)= dtwoeri(l,14,j,i,ider)*sqrt7
                    enddo
                  enddo
                enddo
              enddo
            endif
          case(4)
            do i= 1,nbfijkl(1)
              do j= 1,nbfijkl(2)
                do l= 1,nbfijkl(4)
                  dtwoeri(l, 1,j,i,1)= dtwoeri(l, 1,j,i,1)-twoeri(l, 1,j,i)*five
                  dtwoeri(l, 2,j,i,1)= dtwoeri(l, 2,j,i,1)-twoeri(l, 2,j,i)*four*sqrtseventh
                  dtwoeri(l, 3,j,i,1)= dtwoeri(l, 3,j,i,1)-twoeri(l, 3,j,i)*four*sqrtseventh
                  dtwoeri(l, 4,j,i,1)= dtwoeri(l, 4,j,i,1)-twoeri(l, 4,j,i)*three*sqrt3inv35
                  dtwoeri(l, 5,j,i,1)= dtwoeri(l, 5,j,i,1)-twoeri(l, 5,j,i)*three*sqrtinv35
                  dtwoeri(l, 6,j,i,1)= dtwoeri(l, 6,j,i,1)-twoeri(l, 6,j,i)*three*sqrt3inv35
                  dtwoeri(l, 7,j,i,1)= dtwoeri(l, 7,j,i,1)-twoeri(l, 7,j,i)*two*sqrtseventh
                  dtwoeri(l, 8,j,i,1)= dtwoeri(l, 8,j,i,1)-twoeri(l, 8,j,i)*two*sqrtinv35
                  dtwoeri(l, 9,j,i,1)= dtwoeri(l, 9,j,i,1)-twoeri(l, 9,j,i)*two*sqrtinv35
                  dtwoeri(l,10,j,i,1)= dtwoeri(l,10,j,i,1)-twoeri(l,10,j,i)*two*sqrtseventh
                  dtwoeri(l,11,j,i,1)= dtwoeri(l,11,j,i,1)-twoeri(l,11,j,i)
                  dtwoeri(l,12,j,i,1)= dtwoeri(l,12,j,i,1)-twoeri(l,12,j,i)*sqrtseventh
                  dtwoeri(l,13,j,i,1)= dtwoeri(l,13,j,i,1)-twoeri(l,13,j,i)*sqrt3inv35
                  dtwoeri(l,14,j,i,1)= dtwoeri(l,14,j,i,1)-twoeri(l,14,j,i)*sqrtseventh
                  dtwoeri(l,15,j,i,1)= dtwoeri(l,15,j,i,1)-twoeri(l,15,j,i)
                  dtwoeri(l, 2,j,i,2)= dtwoeri(l, 2,j,i,2)-twoeri(l, 1,j,i)
                  dtwoeri(l, 4,j,i,2)= dtwoeri(l, 4,j,i,2)-twoeri(l, 2,j,i)*two*sqrtseventh
                  dtwoeri(l, 5,j,i,2)= dtwoeri(l, 5,j,i,2)-twoeri(l, 3,j,i)*sqrtseventh
                  dtwoeri(l, 7,j,i,2)= dtwoeri(l, 7,j,i,2)-twoeri(l, 4,j,i)*three*sqrt3inv35
                  dtwoeri(l, 8,j,i,2)= dtwoeri(l, 8,j,i,2)-twoeri(l, 5,j,i)*two*sqrtinv35
                  dtwoeri(l, 9,j,i,2)= dtwoeri(l, 9,j,i,2)-twoeri(l, 6,j,i)*sqrt3inv35
                  dtwoeri(l,11,j,i,2)= dtwoeri(l,11,j,i,2)-twoeri(l, 7,j,i)*four*sqrtseventh
                  dtwoeri(l,12,j,i,2)= dtwoeri(l,12,j,i,2)-twoeri(l, 8,j,i)*three*sqrtinv35
                  dtwoeri(l,13,j,i,2)= dtwoeri(l,13,j,i,2)-twoeri(l, 9,j,i)*two*sqrtinv35
                  dtwoeri(l,14,j,i,2)= dtwoeri(l,14,j,i,2)-twoeri(l,10,j,i)*sqrtseventh
                  dtwoeri(l,16,j,i,2)= dtwoeri(l,16,j,i,2)-twoeri(l,11,j,i)*five
                  dtwoeri(l,17,j,i,2)= dtwoeri(l,17,j,i,2)-twoeri(l,12,j,i)*four*sqrtseventh
                  dtwoeri(l,18,j,i,2)= dtwoeri(l,18,j,i,2)-twoeri(l,13,j,i)*three*sqrt3inv35
                  dtwoeri(l,19,j,i,2)= dtwoeri(l,19,j,i,2)-twoeri(l,14,j,i)*two*sqrtseventh
                  dtwoeri(l,20,j,i,2)= dtwoeri(l,20,j,i,2)-twoeri(l,15,j,i)
                  dtwoeri(l, 3,j,i,3)= dtwoeri(l, 3,j,i,3)-twoeri(l, 1,j,i)
                  dtwoeri(l, 5,j,i,3)= dtwoeri(l, 5,j,i,3)-twoeri(l, 2,j,i)*sqrtseventh
                  dtwoeri(l, 6,j,i,3)= dtwoeri(l, 6,j,i,3)-twoeri(l, 3,j,i)*two*sqrtseventh
                  dtwoeri(l, 8,j,i,3)= dtwoeri(l, 8,j,i,3)-twoeri(l, 4,j,i)*sqrt3inv35
                  dtwoeri(l, 9,j,i,3)= dtwoeri(l, 9,j,i,3)-twoeri(l, 5,j,i)*two*sqrtinv35
                  dtwoeri(l,10,j,i,3)= dtwoeri(l,10,j,i,3)-twoeri(l, 6,j,i)*three*sqrt3inv35
                  dtwoeri(l,12,j,i,3)= dtwoeri(l,12,j,i,3)-twoeri(l, 7,j,i)*sqrtseventh
                  dtwoeri(l,13,j,i,3)= dtwoeri(l,13,j,i,3)-twoeri(l, 8,j,i)*two*sqrtinv35
                  dtwoeri(l,14,j,i,3)= dtwoeri(l,14,j,i,3)-twoeri(l, 9,j,i)*three*sqrtinv35
                  dtwoeri(l,15,j,i,3)= dtwoeri(l,15,j,i,3)-twoeri(l,10,j,i)*four*sqrtseventh
                  dtwoeri(l,17,j,i,3)= dtwoeri(l,17,j,i,3)-twoeri(l,11,j,i)
                  dtwoeri(l,18,j,i,3)= dtwoeri(l,18,j,i,3)-twoeri(l,12,j,i)*two*sqrtseventh
                  dtwoeri(l,19,j,i,3)= dtwoeri(l,19,j,i,3)-twoeri(l,13,j,i)*three*sqrt3inv35
                  dtwoeri(l,20,j,i,3)= dtwoeri(l,20,j,i,3)-twoeri(l,14,j,i)*four*sqrtseventh
                  dtwoeri(l,21,j,i,3)= dtwoeri(l,21,j,i,3)-twoeri(l,15,j,i)*five
                enddo
              enddo
            enddo
            if(mbf(ksh) == 11) then
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do j= 1,nbfijkl(2)
                    do l= 1,nbfijkl(4)
                      do k= 1,21
                        work(k)= dtwoeri(l,k,j,i,ider)
                      enddo
                      dtwoeri(l, 1,j,i,ider)=(work(2)*five-work(7)*ten+work(16))*fach1
                      dtwoeri(l, 2,j,i,ider)=(work(5)*four-work(12)*four)*fach2
                      dtwoeri(l, 3,j,i,ider)=(-work(2)*three-work(7)*two+work(9)*p24 &
&                                            +work(16)-work(18)*eight)*fach3
                      dtwoeri(l, 4,j,i,ider)=(-work(5)*two-work(12)*two+work(14)*four)*fach4
                      dtwoeri(l, 5,j,i,ider)=(work(2)+work(7)*two-work(9)*twelve+work(16) &
&                                            -work(18)*twelve+work(20)*eight)*fach5
                      dtwoeri(l, 6,j,i,ider)=(work(3)*p15+work(8)*p30-work(10)*p40 &
&                                            +work(17)*p15-work(19)*p40+work(21)*eight)*eighth
                      dtwoeri(l, 7,j,i,ider)=(work(1)+work(4)*two-work(6)*twelve+work(11) &
&                                            -work(13)*twelve+work(15)*eight)*fach5
                      dtwoeri(l, 8,j,i,ider)=(-work(3)+work(10)*two+work(17)-work(19)*two)*fach4
                      dtwoeri(l, 9,j,i,ider)=(-work(1)+work(4)*two+work(6)*eight &
&                                            +work(11)*three-work(13)*p24)*fach3
                      dtwoeri(l,10,j,i,ider)=(work(3)-work(8)*six+work(17))*fach2
                      dtwoeri(l,11,j,i,ider)=(work(1)-work(4)*ten+work(11)*five)*fach1
                    enddo
                  enddo
                enddo
              enddo
            else
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do j= 1,nbfijkl(2)
                    do l= 1,nbfijkl(4)
                      dtwoeri(l, 2,j,i,ider)= dtwoeri(l, 2,j,i,ider)*three
                      dtwoeri(l, 3,j,i,ider)= dtwoeri(l, 3,j,i,ider)*three
                      dtwoeri(l, 4,j,i,ider)= dtwoeri(l, 4,j,i,ider)*sqrt21
                      dtwoeri(l, 5,j,i,ider)= dtwoeri(l, 5,j,i,ider)*sqrt63
                      dtwoeri(l, 6,j,i,ider)= dtwoeri(l, 6,j,i,ider)*sqrt21
                      dtwoeri(l, 7,j,i,ider)= dtwoeri(l, 7,j,i,ider)*sqrt21
                      dtwoeri(l, 8,j,i,ider)= dtwoeri(l, 8,j,i,ider)*sqrt105
                      dtwoeri(l, 9,j,i,ider)= dtwoeri(l, 9,j,i,ider)*sqrt105
                      dtwoeri(l,10,j,i,ider)= dtwoeri(l,10,j,i,ider)*sqrt21
                      dtwoeri(l,11,j,i,ider)= dtwoeri(l,11,j,i,ider)*three
                      dtwoeri(l,12,j,i,ider)= dtwoeri(l,12,j,i,ider)*sqrt63
                      dtwoeri(l,13,j,i,ider)= dtwoeri(l,13,j,i,ider)*sqrt105
                      dtwoeri(l,14,j,i,ider)= dtwoeri(l,14,j,i,ider)*sqrt63
                      dtwoeri(l,15,j,i,ider)= dtwoeri(l,15,j,i,ider)*three
                      dtwoeri(l,17,j,i,ider)= dtwoeri(l,17,j,i,ider)*three
                      dtwoeri(l,18,j,i,ider)= dtwoeri(l,18,j,i,ider)*sqrt21
                      dtwoeri(l,19,j,i,ider)= dtwoeri(l,19,j,i,ider)*sqrt21
                      dtwoeri(l,20,j,i,ider)= dtwoeri(l,20,j,i,ider)*three
                    enddo
                  enddo
                enddo
              enddo
            endif
        end select
      endif
!
      nangijkl(3)= mtype(ksh)
      nbfijkl(3) = mbf(ksh)
      do k= 1,nprimijkl(3)
        coijkl(k,3)= coeff(kloc+k)
      enddo
      do i= 1,nbfijkl(1)
        do j= 1,nbfijkl(2)
          do k= 1,nbfijkl(3)
            do l= 1,nbfijkl(4)
              gradtwo(1,3)= gradtwo(1,3)+pdmtrx(l,k,j,i)*dtwoeri(l,k,j,i,1)
              gradtwo(2,3)= gradtwo(2,3)+pdmtrx(l,k,j,i)*dtwoeri(l,k,j,i,2)
              gradtwo(3,3)= gradtwo(3,3)+pdmtrx(l,k,j,i)*dtwoeri(l,k,j,i,3)
            enddo
          enddo
        enddo
      enddo
!
! Jsh derivative
!
      nangijkl(2)= mtype(jsh)+1
      nbfijkl(2) = ncart(nangijkl(2))
      do j= 1,nprimijkl(2)
        coijkl(j,2)= two*ex(jloc+j)*coeff(jloc+j)
      enddo
!
! Two-electron integral calculation
!
      call int2elec(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxgraddim, &
&                   mxprsh,threshex)
!
      select case(nangijkl(2))
        case(1)
          do i= 1,nbfijkl(1)
            do k= 1,nbfijkl(3)
              do l= 1,nbfijkl(4)
                dtwoeri(l,k,1,i,1)= twoeri(l,k,1,i)
                dtwoeri(l,k,1,i,2)= twoeri(l,k,2,i)
                dtwoeri(l,k,1,i,3)= twoeri(l,k,3,i)
              enddo
            enddo
          enddo
        case(2)
          do i= 1,nbfijkl(1)
            do k= 1,nbfijkl(3)
              do l= 1,nbfijkl(4)
                dtwoeri(l,k,1,i,1)= twoeri(l,k,1,i)
                dtwoeri(l,k,2,i,1)= twoeri(l,k,2,i)*sqrtthird
                dtwoeri(l,k,3,i,1)= twoeri(l,k,3,i)*sqrtthird
                dtwoeri(l,k,1,i,2)= twoeri(l,k,2,i)*sqrtthird
                dtwoeri(l,k,2,i,2)= twoeri(l,k,4,i)
                dtwoeri(l,k,3,i,2)= twoeri(l,k,5,i)*sqrtthird
                dtwoeri(l,k,1,i,3)= twoeri(l,k,3,i)*sqrtthird
                dtwoeri(l,k,2,i,3)= twoeri(l,k,5,i)*sqrtthird
                dtwoeri(l,k,3,i,3)= twoeri(l,k,6,i)
              enddo
            enddo
          enddo
        case(3)
          do i= 1,nbfijkl(1)
            do k= 1,nbfijkl(3)
              do l= 1,nbfijkl(4)
                dtwoeri(l,k,1,i,1)= twoeri(l,k, 1,i)
                dtwoeri(l,k,2,i,1)= twoeri(l,k, 2,i)*sqrt3fifth
                dtwoeri(l,k,3,i,1)= twoeri(l,k, 3,i)*sqrt3fifth
                dtwoeri(l,k,4,i,1)= twoeri(l,k, 4,i)*sqrtfifth
                dtwoeri(l,k,5,i,1)= twoeri(l,k, 5,i)*sqrtfifth
                dtwoeri(l,k,6,i,1)= twoeri(l,k, 6,i)*sqrtfifth
                dtwoeri(l,k,1,i,2)= twoeri(l,k, 2,i)*sqrtfifth
                dtwoeri(l,k,2,i,2)= twoeri(l,k, 4,i)*sqrt3fifth
                dtwoeri(l,k,3,i,2)= twoeri(l,k, 5,i)*sqrtfifth
                dtwoeri(l,k,4,i,2)= twoeri(l,k, 7,i)
                dtwoeri(l,k,5,i,2)= twoeri(l,k, 8,i)*sqrt3fifth
                dtwoeri(l,k,6,i,2)= twoeri(l,k, 9,i)*sqrtfifth
                dtwoeri(l,k,1,i,3)= twoeri(l,k, 3,i)*sqrtfifth
                dtwoeri(l,k,2,i,3)= twoeri(l,k, 5,i)*sqrtfifth
                dtwoeri(l,k,3,i,3)= twoeri(l,k, 6,i)*sqrt3fifth
                dtwoeri(l,k,4,i,3)= twoeri(l,k, 8,i)*sqrtfifth
                dtwoeri(l,k,5,i,3)= twoeri(l,k, 9,i)*sqrt3fifth
                dtwoeri(l,k,6,i,3)= twoeri(l,k,10,i)
              enddo
            enddo
          enddo
        case(4)
          do i= 1,nbfijkl(1)
            do k= 1,nbfijkl(3)
              do l= 1,nbfijkl(4)
                dtwoeri(l,k, 1,i,1)= twoeri(l,k, 1,i)
                dtwoeri(l,k, 2,i,1)= twoeri(l,k, 2,i)*sqrtseventh
                dtwoeri(l,k, 3,i,1)= twoeri(l,k, 3,i)*sqrtseventh
                dtwoeri(l,k, 4,i,1)= twoeri(l,k, 4,i)*sqrt3inv35
                dtwoeri(l,k, 5,i,1)= twoeri(l,k, 5,i)*sqrtinv35
                dtwoeri(l,k, 6,i,1)= twoeri(l,k, 6,i)*sqrt3inv35
                dtwoeri(l,k, 7,i,1)= twoeri(l,k, 7,i)*sqrtseventh
                dtwoeri(l,k, 8,i,1)= twoeri(l,k, 8,i)*sqrtinv35
                dtwoeri(l,k, 9,i,1)= twoeri(l,k, 9,i)*sqrtinv35
                dtwoeri(l,k,10,i,1)= twoeri(l,k,10,i)*sqrtseventh
                dtwoeri(l,k, 1,i,2)= twoeri(l,k, 2,i)*sqrtseventh
                dtwoeri(l,k, 2,i,2)= twoeri(l,k, 4,i)*sqrt3inv35
                dtwoeri(l,k, 3,i,2)= twoeri(l,k, 5,i)*sqrtinv35
                dtwoeri(l,k, 4,i,2)= twoeri(l,k, 7,i)*sqrtseventh
                dtwoeri(l,k, 5,i,2)= twoeri(l,k, 8,i)*sqrtinv35
                dtwoeri(l,k, 6,i,2)= twoeri(l,k, 9,i)*sqrtinv35
                dtwoeri(l,k, 7,i,2)= twoeri(l,k,11,i)
                dtwoeri(l,k, 8,i,2)= twoeri(l,k,12,i)*sqrtseventh
                dtwoeri(l,k, 9,i,2)= twoeri(l,k,13,i)*sqrt3inv35
                dtwoeri(l,k,10,i,2)= twoeri(l,k,14,i)*sqrtseventh
                dtwoeri(l,k, 1,i,3)= twoeri(l,k, 3,i)*sqrtseventh
                dtwoeri(l,k, 2,i,3)= twoeri(l,k, 5,i)*sqrtinv35
                dtwoeri(l,k, 3,i,3)= twoeri(l,k, 6,i)*sqrt3inv35
                dtwoeri(l,k, 4,i,3)= twoeri(l,k, 8,i)*sqrtinv35
                dtwoeri(l,k, 5,i,3)= twoeri(l,k, 9,i)*sqrtinv35
                dtwoeri(l,k, 6,i,3)= twoeri(l,k,10,i)*sqrtseventh
                dtwoeri(l,k, 7,i,3)= twoeri(l,k,12,i)*sqrtseventh
                dtwoeri(l,k, 8,i,3)= twoeri(l,k,13,i)*sqrt3inv35
                dtwoeri(l,k, 9,i,3)= twoeri(l,k,14,i)*sqrtseventh
                dtwoeri(l,k,10,i,3)= twoeri(l,k,15,i)
              enddo
            enddo
          enddo
        case(5)
          do i= 1,nbfijkl(1)
            do k= 1,nbfijkl(3)
              do l= 1,nbfijkl(4)
                dtwoeri(l,k, 1,i,1)= twoeri(l,k, 1,i)
                dtwoeri(l,k, 2,i,1)= twoeri(l,k, 2,i)*third
                dtwoeri(l,k, 3,i,1)= twoeri(l,k, 3,i)*third
                dtwoeri(l,k, 4,i,1)= twoeri(l,k, 4,i)*sqrtinv21
                dtwoeri(l,k, 5,i,1)= twoeri(l,k, 5,i)*sqrtinv63
                dtwoeri(l,k, 6,i,1)= twoeri(l,k, 6,i)*sqrtinv21
                dtwoeri(l,k, 7,i,1)= twoeri(l,k, 7,i)*sqrtinv21
                dtwoeri(l,k, 8,i,1)= twoeri(l,k, 8,i)*sqrtinv105
                dtwoeri(l,k, 9,i,1)= twoeri(l,k, 9,i)*sqrtinv105
                dtwoeri(l,k,10,i,1)= twoeri(l,k,10,i)*sqrtinv21
                dtwoeri(l,k,11,i,1)= twoeri(l,k,11,i)*third
                dtwoeri(l,k,12,i,1)= twoeri(l,k,12,i)*sqrtinv63
                dtwoeri(l,k,13,i,1)= twoeri(l,k,13,i)*sqrtinv105
                dtwoeri(l,k,14,i,1)= twoeri(l,k,14,i)*sqrtinv63
                dtwoeri(l,k,15,i,1)= twoeri(l,k,15,i)*third
                dtwoeri(l,k, 1,i,2)= twoeri(l,k, 2,i)*third
                dtwoeri(l,k, 2,i,2)= twoeri(l,k, 4,i)*sqrtinv21
                dtwoeri(l,k, 3,i,2)= twoeri(l,k, 5,i)*sqrtinv63
                dtwoeri(l,k, 4,i,2)= twoeri(l,k, 7,i)*sqrtinv21
                dtwoeri(l,k, 5,i,2)= twoeri(l,k, 8,i)*sqrtinv105
                dtwoeri(l,k, 6,i,2)= twoeri(l,k, 9,i)*sqrtinv105
                dtwoeri(l,k, 7,i,2)= twoeri(l,k,11,i)*third
                dtwoeri(l,k, 8,i,2)= twoeri(l,k,12,i)*sqrtinv63
                dtwoeri(l,k, 9,i,2)= twoeri(l,k,13,i)*sqrtinv105
                dtwoeri(l,k,10,i,2)= twoeri(l,k,14,i)*sqrtinv63
                dtwoeri(l,k,11,i,2)= twoeri(l,k,16,i)
                dtwoeri(l,k,12,i,2)= twoeri(l,k,17,i)*third
                dtwoeri(l,k,13,i,2)= twoeri(l,k,18,i)*sqrtinv21
                dtwoeri(l,k,14,i,2)= twoeri(l,k,19,i)*sqrtinv21
                dtwoeri(l,k,15,i,2)= twoeri(l,k,20,i)*third
                dtwoeri(l,k, 1,i,3)= twoeri(l,k, 3,i)*third
                dtwoeri(l,k, 2,i,3)= twoeri(l,k, 5,i)*sqrtinv63
                dtwoeri(l,k, 3,i,3)= twoeri(l,k, 6,i)*sqrtinv21
                dtwoeri(l,k, 4,i,3)= twoeri(l,k, 8,i)*sqrtinv105
                dtwoeri(l,k, 5,i,3)= twoeri(l,k, 9,i)*sqrtinv105
                dtwoeri(l,k, 6,i,3)= twoeri(l,k,10,i)*sqrtinv21
                dtwoeri(l,k, 7,i,3)= twoeri(l,k,12,i)*sqrtinv63
                dtwoeri(l,k, 8,i,3)= twoeri(l,k,13,i)*sqrtinv105
                dtwoeri(l,k, 9,i,3)= twoeri(l,k,14,i)*sqrtinv63
                dtwoeri(l,k,10,i,3)= twoeri(l,k,15,i)*third
                dtwoeri(l,k,11,i,3)= twoeri(l,k,17,i)*third
                dtwoeri(l,k,12,i,3)= twoeri(l,k,18,i)*sqrtinv21
                dtwoeri(l,k,13,i,3)= twoeri(l,k,19,i)*sqrtinv21
                dtwoeri(l,k,14,i,3)= twoeri(l,k,20,i)*third
                dtwoeri(l,k,15,i,3)= twoeri(l,k,21,i)
              enddo
            enddo
          enddo
        case(6)
          do i= 1,nbfijkl(1)
            do k= 1,nbfijkl(3)
              do l= 1,nbfijkl(4)
                dtwoeri(l,k, 1,i,1)= twoeri(l,k, 1,i)
                dtwoeri(l,k, 2,i,1)= twoeri(l,k, 2,i)*sqrtinv11
                dtwoeri(l,k, 3,i,1)= twoeri(l,k, 3,i)*sqrtinv11
                dtwoeri(l,k, 4,i,1)= twoeri(l,k, 4,i)*sqrtinv33
                dtwoeri(l,k, 5,i,1)= twoeri(l,k, 5,i)*sqrtinv99
                dtwoeri(l,k, 6,i,1)= twoeri(l,k, 6,i)*sqrtinv33
                dtwoeri(l,k, 7,i,1)= twoeri(l,k, 7,i)*sqrt5inv231
                dtwoeri(l,k, 8,i,1)= twoeri(l,k, 8,i)*sqrtinv231
                dtwoeri(l,k, 9,i,1)= twoeri(l,k, 9,i)*sqrtinv231
                dtwoeri(l,k,10,i,1)= twoeri(l,k,10,i)*sqrt5inv231
                dtwoeri(l,k,11,i,1)= twoeri(l,k,11,i)*sqrtinv33
                dtwoeri(l,k,12,i,1)= twoeri(l,k,12,i)*sqrtinv231
                dtwoeri(l,k,13,i,1)= twoeri(l,k,13,i)*sqrtinv385
                dtwoeri(l,k,14,i,1)= twoeri(l,k,14,i)*sqrtinv231
                dtwoeri(l,k,15,i,1)= twoeri(l,k,15,i)*sqrtinv33
                dtwoeri(l,k,16,i,1)= twoeri(l,k,16,i)*sqrtinv11
                dtwoeri(l,k,17,i,1)= twoeri(l,k,17,i)*sqrtinv99
                dtwoeri(l,k,18,i,1)= twoeri(l,k,18,i)*sqrtinv231
                dtwoeri(l,k,19,i,1)= twoeri(l,k,19,i)*sqrtinv231
                dtwoeri(l,k,20,i,1)= twoeri(l,k,20,i)*sqrtinv99
                dtwoeri(l,k,21,i,1)= twoeri(l,k,21,i)*sqrtinv11
                dtwoeri(l,k, 1,i,2)= twoeri(l,k, 2,i)*sqrtinv11
                dtwoeri(l,k, 2,i,2)= twoeri(l,k, 4,i)*sqrtinv33
                dtwoeri(l,k, 3,i,2)= twoeri(l,k, 5,i)*sqrtinv99
                dtwoeri(l,k, 4,i,2)= twoeri(l,k, 7,i)*sqrt5inv231
                dtwoeri(l,k, 5,i,2)= twoeri(l,k, 8,i)*sqrtinv231
                dtwoeri(l,k, 6,i,2)= twoeri(l,k, 9,i)*sqrtinv231
                dtwoeri(l,k, 7,i,2)= twoeri(l,k,11,i)*sqrtinv33
                dtwoeri(l,k, 8,i,2)= twoeri(l,k,12,i)*sqrtinv231
                dtwoeri(l,k, 9,i,2)= twoeri(l,k,13,i)*sqrtinv385
                dtwoeri(l,k,10,i,2)= twoeri(l,k,14,i)*sqrtinv231
                dtwoeri(l,k,11,i,2)= twoeri(l,k,16,i)*sqrtinv11
                dtwoeri(l,k,12,i,2)= twoeri(l,k,17,i)*sqrtinv99
                dtwoeri(l,k,13,i,2)= twoeri(l,k,18,i)*sqrtinv231
                dtwoeri(l,k,14,i,2)= twoeri(l,k,19,i)*sqrtinv231
                dtwoeri(l,k,15,i,2)= twoeri(l,k,20,i)*sqrtinv99
                dtwoeri(l,k,16,i,2)= twoeri(l,k,22,i)
                dtwoeri(l,k,17,i,2)= twoeri(l,k,23,i)*sqrtinv11
                dtwoeri(l,k,18,i,2)= twoeri(l,k,24,i)*sqrtinv33
                dtwoeri(l,k,19,i,2)= twoeri(l,k,25,i)*sqrt5inv231
                dtwoeri(l,k,20,i,2)= twoeri(l,k,26,i)*sqrtinv33
                dtwoeri(l,k,21,i,2)= twoeri(l,k,27,i)*sqrtinv11
                dtwoeri(l,k, 1,i,3)= twoeri(l,k, 3,i)*sqrtinv11
                dtwoeri(l,k, 2,i,3)= twoeri(l,k, 5,i)*sqrtinv99
                dtwoeri(l,k, 3,i,3)= twoeri(l,k, 6,i)*sqrtinv33
                dtwoeri(l,k, 4,i,3)= twoeri(l,k, 8,i)*sqrtinv231
                dtwoeri(l,k, 5,i,3)= twoeri(l,k, 9,i)*sqrtinv231
                dtwoeri(l,k, 6,i,3)= twoeri(l,k,10,i)*sqrt5inv231
                dtwoeri(l,k, 7,i,3)= twoeri(l,k,12,i)*sqrtinv231
                dtwoeri(l,k, 8,i,3)= twoeri(l,k,13,i)*sqrtinv385
                dtwoeri(l,k, 9,i,3)= twoeri(l,k,14,i)*sqrtinv231
                dtwoeri(l,k,10,i,3)= twoeri(l,k,15,i)*sqrtinv33
                dtwoeri(l,k,11,i,3)= twoeri(l,k,17,i)*sqrtinv99
                dtwoeri(l,k,12,i,3)= twoeri(l,k,18,i)*sqrtinv231
                dtwoeri(l,k,13,i,3)= twoeri(l,k,19,i)*sqrtinv231
                dtwoeri(l,k,14,i,3)= twoeri(l,k,20,i)*sqrtinv99
                dtwoeri(l,k,15,i,3)= twoeri(l,k,21,i)*sqrtinv11
                dtwoeri(l,k,16,i,3)= twoeri(l,k,23,i)*sqrtinv11
                dtwoeri(l,k,17,i,3)= twoeri(l,k,24,i)*sqrtinv33
                dtwoeri(l,k,18,i,3)= twoeri(l,k,25,i)*sqrt5inv231
                dtwoeri(l,k,19,i,3)= twoeri(l,k,26,i)*sqrtinv33
                dtwoeri(l,k,20,i,3)= twoeri(l,k,27,i)*sqrtinv11
                dtwoeri(l,k,21,i,3)= twoeri(l,k,28,i)
              enddo
            enddo
          enddo
        case default
          write(*,'(" Error! This program supports up to h function in calcd2eri")')
          call iabort
      end select
!
      if(mtype(jsh) >= 1) then
        nangijkl(2)= mtype(jsh)-1
        nbfijkl(2) = ncart(nangijkl(2))
        do j= 1,nprimijkl(2)
          coijkl(j,2)= coeff(jloc+j)
        enddo
!
! Two-electron integral calculation
!
        call int2elec(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxgraddim, &
&                     mxprsh,threshex)
!
        select case(nangijkl(2))
          case(0)
            do i= 1,nbfijkl(1)
              do k= 1,nbfijkl(3)
                do l= 1,nbfijkl(4)
                  dtwoeri(l,k,1,i,1)= dtwoeri(l,k,1,i,1)-twoeri(l,k,1,i)
                  dtwoeri(l,k,2,i,2)= dtwoeri(l,k,2,i,2)-twoeri(l,k,1,i)
                  dtwoeri(l,k,3,i,3)= dtwoeri(l,k,3,i,3)-twoeri(l,k,1,i)
                enddo
              enddo
            enddo
          case(1)
            do i= 1,nbfijkl(1)
              do k= 1,nbfijkl(3)
                do l= 1,nbfijkl(4)
                  dtwoeri(l,k,1,i,1)= dtwoeri(l,k,1,i,1)-twoeri(l,k,1,i)*two
                  dtwoeri(l,k,2,i,1)= dtwoeri(l,k,2,i,1)-twoeri(l,k,2,i)*sqrt3
                  dtwoeri(l,k,3,i,1)= dtwoeri(l,k,3,i,1)-twoeri(l,k,3,i)*sqrt3
                  dtwoeri(l,k,2,i,2)= dtwoeri(l,k,2,i,2)-twoeri(l,k,1,i)*sqrt3
                  dtwoeri(l,k,4,i,2)= dtwoeri(l,k,4,i,2)-twoeri(l,k,2,i)*two
                  dtwoeri(l,k,5,i,2)= dtwoeri(l,k,5,i,2)-twoeri(l,k,3,i)*sqrt3
                  dtwoeri(l,k,3,i,3)= dtwoeri(l,k,3,i,3)-twoeri(l,k,1,i)*sqrt3
                  dtwoeri(l,k,5,i,3)= dtwoeri(l,k,5,i,3)-twoeri(l,k,2,i)*sqrt3
                  dtwoeri(l,k,6,i,3)= dtwoeri(l,k,6,i,3)-twoeri(l,k,3,i)*two
                enddo
              enddo
            enddo
            if(mbf(jsh) == 5) then
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do k= 1,nbfijkl(3)
                    do l= 1,nbfijkl(4)
                      do j= 1,6
                        work(j)= dtwoeri(l,k,j,i,ider)
                      enddo
                      dtwoeri(l,k,1,i,ider)= work(2)
                      dtwoeri(l,k,2,i,ider)= work(5)
                      dtwoeri(l,k,3,i,ider)=(work(6)*two-work(1)-work(4))*half
                      dtwoeri(l,k,4,i,ider)= work(3)
                      dtwoeri(l,k,5,i,ider)=(work(1)-work(4))*sqrt3h
                    enddo
                  enddo
                enddo
              enddo
            endif
          case(2)
            do i= 1,nbfijkl(1)
              do k= 1,nbfijkl(3)
                do l= 1,nbfijkl(4)
                  dtwoeri(l,k, 1,i,1)= dtwoeri(l,k, 1,i,1)-twoeri(l,k,1,i)*three
                  dtwoeri(l,k, 2,i,1)= dtwoeri(l,k, 2,i,1)-twoeri(l,k,2,i)*two*sqrtthird
                  dtwoeri(l,k, 3,i,1)= dtwoeri(l,k, 3,i,1)-twoeri(l,k,3,i)*two*sqrtthird
                  dtwoeri(l,k, 4,i,1)= dtwoeri(l,k, 4,i,1)-twoeri(l,k,4,i)
                  dtwoeri(l,k, 5,i,1)= dtwoeri(l,k, 5,i,1)-twoeri(l,k,5,i)*sqrtthird
                  dtwoeri(l,k, 6,i,1)= dtwoeri(l,k, 6,i,1)-twoeri(l,k,6,i)
                  dtwoeri(l,k, 2,i,2)= dtwoeri(l,k, 2,i,2)-twoeri(l,k,1,i)
                  dtwoeri(l,k, 4,i,2)= dtwoeri(l,k, 4,i,2)-twoeri(l,k,2,i)*two*sqrtthird
                  dtwoeri(l,k, 5,i,2)= dtwoeri(l,k, 5,i,2)-twoeri(l,k,3,i)*sqrtthird
                  dtwoeri(l,k, 7,i,2)= dtwoeri(l,k, 7,i,2)-twoeri(l,k,4,i)*three
                  dtwoeri(l,k, 8,i,2)= dtwoeri(l,k, 8,i,2)-twoeri(l,k,5,i)*two*sqrtthird
                  dtwoeri(l,k, 9,i,2)= dtwoeri(l,k, 9,i,2)-twoeri(l,k,6,i)
                  dtwoeri(l,k, 3,i,3)= dtwoeri(l,k, 3,i,3)-twoeri(l,k,1,i)
                  dtwoeri(l,k, 5,i,3)= dtwoeri(l,k, 5,i,3)-twoeri(l,k,2,i)*sqrtthird
                  dtwoeri(l,k, 6,i,3)= dtwoeri(l,k, 6,i,3)-twoeri(l,k,3,i)*two*sqrtthird
                  dtwoeri(l,k, 8,i,3)= dtwoeri(l,k, 8,i,3)-twoeri(l,k,4,i)
                  dtwoeri(l,k, 9,i,3)= dtwoeri(l,k, 9,i,3)-twoeri(l,k,5,i)*two*sqrtthird
                  dtwoeri(l,k,10,i,3)= dtwoeri(l,k,10,i,3)-twoeri(l,k,6,i)*three
                enddo
              enddo
            enddo
            if(mbf(jsh) == 7) then
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do k= 1,nbfijkl(3)
                    do l= 1,nbfijkl(4)
                      do j= 1,10
                        work(j)= dtwoeri(l,k,j,i,ider)
                      enddo
                      dtwoeri(l,k,1,i,ider)=(-work(7)+work(2)*three                   )*facf1
                      dtwoeri(l,k,2,i,ider)=  work(5)                                  *facf2
                      dtwoeri(l,k,3,i,ider)=(-work(7)-work(2)+work(9)*four            )*facf3
                      dtwoeri(l,k,4,i,ider)=( work(10)*two-work(3)*three-work(8)*three)*half
                      dtwoeri(l,k,5,i,ider)=(-work(1)-work(4)+work(6)*four            )*facf3
                      dtwoeri(l,k,6,i,ider)=( work(3)-work(8)                         )*facf4
                      dtwoeri(l,k,7,i,ider)=( work(1)-work(4)*three                   )*facf1
                    enddo
                  enddo
                enddo
              enddo
            else
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do k= 1,nbfijkl(3)
                    do l= 1,nbfijkl(4)
                      dtwoeri(l,k,2,i,ider)= dtwoeri(l,k,2,i,ider)*sqrt5
                      dtwoeri(l,k,3,i,ider)= dtwoeri(l,k,3,i,ider)*sqrt5
                      dtwoeri(l,k,4,i,ider)= dtwoeri(l,k,4,i,ider)*sqrt5
                      dtwoeri(l,k,5,i,ider)= dtwoeri(l,k,5,i,ider)*sqrt15
                      dtwoeri(l,k,6,i,ider)= dtwoeri(l,k,6,i,ider)*sqrt5
                      dtwoeri(l,k,8,i,ider)= dtwoeri(l,k,8,i,ider)*sqrt5
                      dtwoeri(l,k,9,i,ider)= dtwoeri(l,k,9,i,ider)*sqrt5
                    enddo
                  enddo
                enddo
              enddo
            endif
          case(3)
            do i= 1,nbfijkl(1)
              do k= 1,nbfijkl(3)
                do l= 1,nbfijkl(4)
                  dtwoeri(l,k, 1,i,1)= dtwoeri(l,k, 1,i,1)-twoeri(l,k, 1,i)*four
                  dtwoeri(l,k, 2,i,1)= dtwoeri(l,k, 2,i,1)-twoeri(l,k, 2,i)*three*sqrtfifth
                  dtwoeri(l,k, 3,i,1)= dtwoeri(l,k, 3,i,1)-twoeri(l,k, 3,i)*three*sqrtfifth
                  dtwoeri(l,k, 4,i,1)= dtwoeri(l,k, 4,i,1)-twoeri(l,k, 4,i)*two*sqrtfifth
                  dtwoeri(l,k, 5,i,1)= dtwoeri(l,k, 5,i,1)-twoeri(l,k, 5,i)*two*sqrtinv15
                  dtwoeri(l,k, 6,i,1)= dtwoeri(l,k, 6,i,1)-twoeri(l,k, 6,i)*two*sqrtfifth
                  dtwoeri(l,k, 7,i,1)= dtwoeri(l,k, 7,i,1)-twoeri(l,k, 7,i)
                  dtwoeri(l,k, 8,i,1)= dtwoeri(l,k, 8,i,1)-twoeri(l,k, 8,i)*sqrtfifth
                  dtwoeri(l,k, 9,i,1)= dtwoeri(l,k, 9,i,1)-twoeri(l,k, 9,i)*sqrtfifth
                  dtwoeri(l,k,10,i,1)= dtwoeri(l,k,10,i,1)-twoeri(l,k,10,i)
                  dtwoeri(l,k, 2,i,2)= dtwoeri(l,k, 2,i,2)-twoeri(l,k, 1,i)
                  dtwoeri(l,k, 4,i,2)= dtwoeri(l,k, 4,i,2)-twoeri(l,k, 2,i)*two*sqrtfifth
                  dtwoeri(l,k, 5,i,2)= dtwoeri(l,k, 5,i,2)-twoeri(l,k, 3,i)*sqrtfifth
                  dtwoeri(l,k, 7,i,2)= dtwoeri(l,k, 7,i,2)-twoeri(l,k, 4,i)*three*sqrtfifth
                  dtwoeri(l,k, 8,i,2)= dtwoeri(l,k, 8,i,2)-twoeri(l,k, 5,i)*two*sqrtinv15
                  dtwoeri(l,k, 9,i,2)= dtwoeri(l,k, 9,i,2)-twoeri(l,k, 6,i)*sqrtfifth
                  dtwoeri(l,k,11,i,2)= dtwoeri(l,k,11,i,2)-twoeri(l,k, 7,i)*four
                  dtwoeri(l,k,12,i,2)= dtwoeri(l,k,12,i,2)-twoeri(l,k, 8,i)*three*sqrtfifth
                  dtwoeri(l,k,13,i,2)= dtwoeri(l,k,13,i,2)-twoeri(l,k, 9,i)*two*sqrtfifth
                  dtwoeri(l,k,14,i,2)= dtwoeri(l,k,14,i,2)-twoeri(l,k,10,i)
                  dtwoeri(l,k, 3,i,3)= dtwoeri(l,k, 3,i,3)-twoeri(l,k, 1,i)
                  dtwoeri(l,k, 5,i,3)= dtwoeri(l,k, 5,i,3)-twoeri(l,k, 2,i)*sqrtfifth
                  dtwoeri(l,k, 6,i,3)= dtwoeri(l,k, 6,i,3)-twoeri(l,k, 3,i)*two*sqrtfifth
                  dtwoeri(l,k, 8,i,3)= dtwoeri(l,k, 8,i,3)-twoeri(l,k, 4,i)*sqrtfifth
                  dtwoeri(l,k, 9,i,3)= dtwoeri(l,k, 9,i,3)-twoeri(l,k, 5,i)*two*sqrtinv15
                  dtwoeri(l,k,10,i,3)= dtwoeri(l,k,10,i,3)-twoeri(l,k, 6,i)*three*sqrtfifth
                  dtwoeri(l,k,12,i,3)= dtwoeri(l,k,12,i,3)-twoeri(l,k, 7,i)
                  dtwoeri(l,k,13,i,3)= dtwoeri(l,k,13,i,3)-twoeri(l,k, 8,i)*two*sqrtfifth
                  dtwoeri(l,k,14,i,3)= dtwoeri(l,k,14,i,3)-twoeri(l,k, 9,i)*three*sqrtfifth
                  dtwoeri(l,k,15,i,3)= dtwoeri(l,k,15,i,3)-twoeri(l,k,10,i)*four
                enddo
              enddo
            enddo
            if(mbf(jsh) == 9) then
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do k= 1,nbfijkl(3)
                    do l= 1,nbfijkl(4)
                      do j= 1,15
                        work(j)= dtwoeri(l,k,j,i,ider)
                      enddo
                      dtwoeri(l,k,1,i,ider)=(work(2)-work(7))*facg1
                      dtwoeri(l,k,2,i,ider)=(-work(12)+work(5)*three)*facg2
                      dtwoeri(l,k,3,i,ider)=(-work(2)-work(7)+work(9)*six)*facg3
                      dtwoeri(l,k,4,i,ider)=(-work(12)*three+work(14)*four-work(5)*three)*facg4
                      dtwoeri(l,k,5,i,ider)=(work(1)*three+work(11)*three+work(15)*eight &
&                                           +work(4)*six-work(6)*p24-work(13)*p24)*eighth
                      dtwoeri(l,k,6,i,ider)=(-work(3)*three+work(10)*four-work(8)*three)*facg4
                      dtwoeri(l,k,7,i,ider)=(-work(1)+work(11)+work(6)*six-work(13)*six)*facg5
                      dtwoeri(l,k,8,i,ider)=(work(3)-work(8)*three)*facg2
                      dtwoeri(l,k,9,i,ider)=(work(1)+work(11)-work(4)*six)*facg6
                    enddo
                  enddo
                enddo
              enddo
            else
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do k= 1,nbfijkl(3)
                    do l= 1,nbfijkl(4)
                      dtwoeri(l,k, 2,i,ider)= dtwoeri(l,k, 2,i,ider)*sqrt7
                      dtwoeri(l,k, 3,i,ider)= dtwoeri(l,k, 3,i,ider)*sqrt7
                      dtwoeri(l,k, 4,i,ider)= dtwoeri(l,k, 4,i,ider)*sqrt35third
                      dtwoeri(l,k, 5,i,ider)= dtwoeri(l,k, 5,i,ider)*sqrt35
                      dtwoeri(l,k, 6,i,ider)= dtwoeri(l,k, 6,i,ider)*sqrt35third
                      dtwoeri(l,k, 7,i,ider)= dtwoeri(l,k, 7,i,ider)*sqrt7
                      dtwoeri(l,k, 8,i,ider)= dtwoeri(l,k, 8,i,ider)*sqrt35
                      dtwoeri(l,k, 9,i,ider)= dtwoeri(l,k, 9,i,ider)*sqrt35
                      dtwoeri(l,k,10,i,ider)= dtwoeri(l,k,10,i,ider)*sqrt7
                      dtwoeri(l,k,12,i,ider)= dtwoeri(l,k,12,i,ider)*sqrt7
                      dtwoeri(l,k,13,i,ider)= dtwoeri(l,k,13,i,ider)*sqrt35third
                      dtwoeri(l,k,14,i,ider)= dtwoeri(l,k,14,i,ider)*sqrt7
                    enddo
                  enddo
                enddo
              enddo
            endif
          case(4)
            do i= 1,nbfijkl(1)
              do k= 1,nbfijkl(3)
                do l= 1,nbfijkl(4)
                  dtwoeri(l,k, 1,i,1)= dtwoeri(l,k, 1,i,1)-twoeri(l,k, 1,i)*five
                  dtwoeri(l,k, 2,i,1)= dtwoeri(l,k, 2,i,1)-twoeri(l,k, 2,i)*four*sqrtseventh
                  dtwoeri(l,k, 3,i,1)= dtwoeri(l,k, 3,i,1)-twoeri(l,k, 3,i)*four*sqrtseventh
                  dtwoeri(l,k, 4,i,1)= dtwoeri(l,k, 4,i,1)-twoeri(l,k, 4,i)*three*sqrt3inv35
                  dtwoeri(l,k, 5,i,1)= dtwoeri(l,k, 5,i,1)-twoeri(l,k, 5,i)*three*sqrtinv35
                  dtwoeri(l,k, 6,i,1)= dtwoeri(l,k, 6,i,1)-twoeri(l,k, 6,i)*three*sqrt3inv35
                  dtwoeri(l,k, 7,i,1)= dtwoeri(l,k, 7,i,1)-twoeri(l,k, 7,i)*two*sqrtseventh
                  dtwoeri(l,k, 8,i,1)= dtwoeri(l,k, 8,i,1)-twoeri(l,k, 8,i)*two*sqrtinv35
                  dtwoeri(l,k, 9,i,1)= dtwoeri(l,k, 9,i,1)-twoeri(l,k, 9,i)*two*sqrtinv35
                  dtwoeri(l,k,10,i,1)= dtwoeri(l,k,10,i,1)-twoeri(l,k,10,i)*two*sqrtseventh
                  dtwoeri(l,k,11,i,1)= dtwoeri(l,k,11,i,1)-twoeri(l,k,11,i)
                  dtwoeri(l,k,12,i,1)= dtwoeri(l,k,12,i,1)-twoeri(l,k,12,i)*sqrtseventh
                  dtwoeri(l,k,13,i,1)= dtwoeri(l,k,13,i,1)-twoeri(l,k,13,i)*sqrt3inv35
                  dtwoeri(l,k,14,i,1)= dtwoeri(l,k,14,i,1)-twoeri(l,k,14,i)*sqrtseventh
                  dtwoeri(l,k,15,i,1)= dtwoeri(l,k,15,i,1)-twoeri(l,k,15,i)
                  dtwoeri(l,k, 2,i,2)= dtwoeri(l,k, 2,i,2)-twoeri(l,k, 1,i)
                  dtwoeri(l,k, 4,i,2)= dtwoeri(l,k, 4,i,2)-twoeri(l,k, 2,i)*two*sqrtseventh
                  dtwoeri(l,k, 5,i,2)= dtwoeri(l,k, 5,i,2)-twoeri(l,k, 3,i)*sqrtseventh
                  dtwoeri(l,k, 7,i,2)= dtwoeri(l,k, 7,i,2)-twoeri(l,k, 4,i)*three*sqrt3inv35
                  dtwoeri(l,k, 8,i,2)= dtwoeri(l,k, 8,i,2)-twoeri(l,k, 5,i)*two*sqrtinv35
                  dtwoeri(l,k, 9,i,2)= dtwoeri(l,k, 9,i,2)-twoeri(l,k, 6,i)*sqrt3inv35
                  dtwoeri(l,k,11,i,2)= dtwoeri(l,k,11,i,2)-twoeri(l,k, 7,i)*four*sqrtseventh
                  dtwoeri(l,k,12,i,2)= dtwoeri(l,k,12,i,2)-twoeri(l,k, 8,i)*three*sqrtinv35
                  dtwoeri(l,k,13,i,2)= dtwoeri(l,k,13,i,2)-twoeri(l,k, 9,i)*two*sqrtinv35
                  dtwoeri(l,k,14,i,2)= dtwoeri(l,k,14,i,2)-twoeri(l,k,10,i)*sqrtseventh
                  dtwoeri(l,k,16,i,2)= dtwoeri(l,k,16,i,2)-twoeri(l,k,11,i)*five
                  dtwoeri(l,k,17,i,2)= dtwoeri(l,k,17,i,2)-twoeri(l,k,12,i)*four*sqrtseventh
                  dtwoeri(l,k,18,i,2)= dtwoeri(l,k,18,i,2)-twoeri(l,k,13,i)*three*sqrt3inv35
                  dtwoeri(l,k,19,i,2)= dtwoeri(l,k,19,i,2)-twoeri(l,k,14,i)*two*sqrtseventh
                  dtwoeri(l,k,20,i,2)= dtwoeri(l,k,20,i,2)-twoeri(l,k,15,i)
                  dtwoeri(l,k, 3,i,3)= dtwoeri(l,k, 3,i,3)-twoeri(l,k, 1,i)
                  dtwoeri(l,k, 5,i,3)= dtwoeri(l,k, 5,i,3)-twoeri(l,k, 2,i)*sqrtseventh
                  dtwoeri(l,k, 6,i,3)= dtwoeri(l,k, 6,i,3)-twoeri(l,k, 3,i)*two*sqrtseventh
                  dtwoeri(l,k, 8,i,3)= dtwoeri(l,k, 8,i,3)-twoeri(l,k, 4,i)*sqrt3inv35
                  dtwoeri(l,k, 9,i,3)= dtwoeri(l,k, 9,i,3)-twoeri(l,k, 5,i)*two*sqrtinv35
                  dtwoeri(l,k,10,i,3)= dtwoeri(l,k,10,i,3)-twoeri(l,k, 6,i)*three*sqrt3inv35
                  dtwoeri(l,k,12,i,3)= dtwoeri(l,k,12,i,3)-twoeri(l,k, 7,i)*sqrtseventh
                  dtwoeri(l,k,13,i,3)= dtwoeri(l,k,13,i,3)-twoeri(l,k, 8,i)*two*sqrtinv35
                  dtwoeri(l,k,14,i,3)= dtwoeri(l,k,14,i,3)-twoeri(l,k, 9,i)*three*sqrtinv35
                  dtwoeri(l,k,15,i,3)= dtwoeri(l,k,15,i,3)-twoeri(l,k,10,i)*four*sqrtseventh
                  dtwoeri(l,k,17,i,3)= dtwoeri(l,k,17,i,3)-twoeri(l,k,11,i)
                  dtwoeri(l,k,18,i,3)= dtwoeri(l,k,18,i,3)-twoeri(l,k,12,i)*two*sqrtseventh
                  dtwoeri(l,k,19,i,3)= dtwoeri(l,k,19,i,3)-twoeri(l,k,13,i)*three*sqrt3inv35
                  dtwoeri(l,k,20,i,3)= dtwoeri(l,k,20,i,3)-twoeri(l,k,14,i)*four*sqrtseventh
                  dtwoeri(l,k,21,i,3)= dtwoeri(l,k,21,i,3)-twoeri(l,k,15,i)*five
                enddo
              enddo
            enddo
            if(mbf(jsh) == 11) then
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do k= 1,nbfijkl(3)
                    do l= 1,nbfijkl(4)
                      do j= 1,21
                        work(j)= dtwoeri(l,k,j,i,ider)
                      enddo
                      dtwoeri(l,k, 1,i,ider)=(work(2)*five-work(7)*ten+work(16))*fach1
                      dtwoeri(l,k, 2,i,ider)=(work(5)*four-work(12)*four)*fach2
                      dtwoeri(l,k, 3,i,ider)=(-work(2)*three-work(7)*two+work(9)*p24 &
&                                            +work(16)-work(18)*eight)*fach3
                      dtwoeri(l,k, 4,i,ider)=(-work(5)*two-work(12)*two+work(14)*four)*fach4
                      dtwoeri(l,k, 5,i,ider)=(work(2)+work(7)*two-work(9)*twelve+work(16) &
&                                            -work(18)*twelve+work(20)*eight)*fach5
                      dtwoeri(l,k, 6,i,ider)=(work(3)*p15+work(8)*p30-work(10)*p40 &
&                                            +work(17)*p15-work(19)*p40+work(21)*eight)*eighth
                      dtwoeri(l,k, 7,i,ider)=(work(1)+work(4)*two-work(6)*twelve+work(11) &
&                                            -work(13)*twelve+work(15)*eight)*fach5
                      dtwoeri(l,k, 8,i,ider)=(-work(3)+work(10)*two+work(17)-work(19)*two)*fach4
                      dtwoeri(l,k, 9,i,ider)=(-work(1)+work(4)*two+work(6)*eight &
&                                            +work(11)*three-work(13)*p24)*fach3
                      dtwoeri(l,k,10,i,ider)=(work(3)-work(8)*six+work(17))*fach2
                      dtwoeri(l,k,11,i,ider)=(work(1)-work(4)*ten+work(11)*five)*fach1
                    enddo
                  enddo
                enddo
              enddo
            else
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do k= 1,nbfijkl(3)
                    do l= 1,nbfijkl(4)
                      dtwoeri(l,k, 2,i,ider)= dtwoeri(l,k, 2,i,ider)*three
                      dtwoeri(l,k, 3,i,ider)= dtwoeri(l,k, 3,i,ider)*three
                      dtwoeri(l,k, 4,i,ider)= dtwoeri(l,k, 4,i,ider)*sqrt21
                      dtwoeri(l,k, 5,i,ider)= dtwoeri(l,k, 5,i,ider)*sqrt63
                      dtwoeri(l,k, 6,i,ider)= dtwoeri(l,k, 6,i,ider)*sqrt21
                      dtwoeri(l,k, 7,i,ider)= dtwoeri(l,k, 7,i,ider)*sqrt21
                      dtwoeri(l,k, 8,i,ider)= dtwoeri(l,k, 8,i,ider)*sqrt105
                      dtwoeri(l,k, 9,i,ider)= dtwoeri(l,k, 9,i,ider)*sqrt105
                      dtwoeri(l,k,10,i,ider)= dtwoeri(l,k,10,i,ider)*sqrt21
                      dtwoeri(l,k,11,i,ider)= dtwoeri(l,k,11,i,ider)*three
                      dtwoeri(l,k,12,i,ider)= dtwoeri(l,k,12,i,ider)*sqrt63
                      dtwoeri(l,k,13,i,ider)= dtwoeri(l,k,13,i,ider)*sqrt105
                      dtwoeri(l,k,14,i,ider)= dtwoeri(l,k,14,i,ider)*sqrt63
                      dtwoeri(l,k,15,i,ider)= dtwoeri(l,k,15,i,ider)*three
                      dtwoeri(l,k,17,i,ider)= dtwoeri(l,k,17,i,ider)*three
                      dtwoeri(l,k,18,i,ider)= dtwoeri(l,k,18,i,ider)*sqrt21
                      dtwoeri(l,k,19,i,ider)= dtwoeri(l,k,19,i,ider)*sqrt21
                      dtwoeri(l,k,20,i,ider)= dtwoeri(l,k,20,i,ider)*three
                    enddo
                  enddo
                enddo
              enddo
            endif
        end select
      endif
!
      nbfijkl(2) = mbf(jsh)
      do i= 1,nbfijkl(1)
        do j= 1,nbfijkl(2)
          do k= 1,nbfijkl(3)
            do l= 1,nbfijkl(4)
              gradtwo(1,2)= gradtwo(1,2)+pdmtrx(l,k,j,i)*dtwoeri(l,k,j,i,1)
              gradtwo(2,2)= gradtwo(2,2)+pdmtrx(l,k,j,i)*dtwoeri(l,k,j,i,2)
              gradtwo(3,2)= gradtwo(3,2)+pdmtrx(l,k,j,i)*dtwoeri(l,k,j,i,3)
            enddo
          enddo
        enddo
      enddo
!
! Ish derivative
!
      gradtwo(1,1)=-gradtwo(1,2)-gradtwo(1,3)-gradtwo(1,4)
      gradtwo(2,1)=-gradtwo(2,2)-gradtwo(2,3)-gradtwo(2,4)
      gradtwo(3,1)=-gradtwo(3,2)-gradtwo(3,3)-gradtwo(3,4)
!
! Add derivative terms of 2-ERI
!
      do i= 1,3
        egrad2(i,iatom)= egrad2(i,iatom)+gradtwo(i,1)
        egrad2(i,jatom)= egrad2(i,jatom)+gradtwo(i,2)
        egrad2(i,katom)= egrad2(i,katom)+gradtwo(i,3)
        egrad2(i,latom)= egrad2(i,latom)+gradtwo(i,4)
      enddo
!
      return
end


!-------------------------------------------------------------------------
  subroutine calcpdmtrx(fulldmtrx1,fulldmtrx2,pdmtrx,pdmax,hfexchange, &
&                       ish,jsh,ksh,lsh,maxdim,itype)
!-------------------------------------------------------------------------
!
! Calculate 4*(4*Dij*Dkl-Dil*Djk-Dik*Djl)
!
      use modbasis, only : nao, mbf, locbf
      implicit none
      integer,intent(in) :: ish, jsh, ksh, lsh, maxdim, itype
      integer :: i, j, k, l, ilocbf, jlocbf, klocbf, llocbf, ii, jj, kk, ll
      integer :: nbfijkl(4)
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, four=4.0D+00
      real(8),intent(in) :: fulldmtrx1(nao,nao), fulldmtrx2(nao,nao), hfexchange
      real(8),intent(out) :: pdmtrx(maxdim,maxdim,maxdim,maxdim), pdmax
      real(8) :: factor
!
      pdmax= zero
      nbfijkl(1)= mbf(ish)
      nbfijkl(2)= mbf(jsh)
      nbfijkl(3)= mbf(ksh)
      nbfijkl(4)= mbf(lsh)
      ilocbf= locbf(ish)
      jlocbf= locbf(jsh)
      klocbf= locbf(ksh)
      llocbf= locbf(lsh)

!
! Check ish, jsh, ksh, lsh
!
      factor= four
      factor= one
      factor= one/four
      if(ish == jsh) factor= factor*half
      if(ksh == lsh) factor= factor*half
      if((ish == ksh).and.(jsh == lsh)) factor= factor*half
!
! HF closed-shell
!
      select case(itype)
        case(1)
!
! 4*(4*Dij*Dkl-Dil*Djk-Dik*Djl)
!
          do i= 1,nbfijkl(1)
            ii= ilocbf+i
            do j= 1,nbfijkl(2)
              jj= jlocbf+j
              do k= 1,nbfijkl(3)
                kk= klocbf+k
                do l= 1,nbfijkl(4)
                  ll= llocbf+l
                  pdmtrx(l,k,j,i)= factor*(four*fulldmtrx1(jj,ii)*fulldmtrx1(ll,kk) &
&                                   -hfexchange*fulldmtrx1(kk,ii)*fulldmtrx1(ll,jj) &
&                                   -hfexchange*fulldmtrx1(ll,ii)*fulldmtrx1(kk,jj))
                  if(pdmax < abs(pdmtrx(l,k,j,i))) pdmax= abs(pdmtrx(l,k,j,i))
                enddo
              enddo
            enddo
          enddo
!
! HF open-shell
!
        case(2)
!
! 4*Dija*Dkla+4*Dijb*Dklb+8*Dija*Dklb-2*Dika*Djla-2*Dila*Djka-2*Dikb*Djlb-2*Dilb*Djkb
!
          do i= 1,nbfijkl(1)
            ii= ilocbf+i
            do j= 1,nbfijkl(2)
              jj= jlocbf+j
              do k= 1,nbfijkl(3)
                kk= klocbf+k
                do l= 1,nbfijkl(4)
                  ll= llocbf+l
                  pdmtrx(l,k,j,i)= factor*(four*fulldmtrx1(jj,ii)*fulldmtrx1(ll,kk) &
&                                   -hfexchange*fulldmtrx1(kk,ii)*fulldmtrx1(ll,jj) &
&                                   -hfexchange*fulldmtrx1(ll,ii)*fulldmtrx1(kk,jj) &
&                                   -hfexchange*fulldmtrx2(kk,ii)*fulldmtrx2(ll,jj) &
&                                   -hfexchange*fulldmtrx2(ll,ii)*fulldmtrx2(kk,jj))
                  if(pdmax < abs(pdmtrx(l,k,j,i))) pdmax= abs(pdmtrx(l,k,j,i))
                enddo
              enddo
            enddo
          enddo
!
! MP2 closed-shell
!
        case(3)
          do i= 1,nbfijkl(1)
            ii= ilocbf+i
            do j= 1,nbfijkl(2)
              jj= jlocbf+j
              do k= 1,nbfijkl(3)
                kk= klocbf+k
                do l= 1,nbfijkl(4)
                  ll= llocbf+l
                  pdmtrx(l,k,j,i)= factor &
&                                 *(four*((fulldmtrx1(jj,ii)+fulldmtrx2(jj,ii))*fulldmtrx1(ll,kk) &
&                                         +fulldmtrx1(jj,ii)*fulldmtrx2(ll,kk)) &
&                                       -((fulldmtrx1(kk,ii)+fulldmtrx2(kk,ii))*fulldmtrx1(ll,jj) &
&                                        +(fulldmtrx1(ll,ii)+fulldmtrx2(ll,ii))*fulldmtrx1(kk,jj) &
&                                         +fulldmtrx1(kk,ii)*fulldmtrx2(ll,jj) &
&                                         +fulldmtrx1(ll,ii)*fulldmtrx2(kk,jj)))
                  if(pdmax < abs(pdmtrx(l,k,j,i))) pdmax= abs(pdmtrx(l,k,j,i))
                enddo
              enddo
            enddo
          enddo
      end select
!
      return
end


