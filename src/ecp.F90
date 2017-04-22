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
!--------------------------------------------
  subroutine oneeiecp(hstmat2,nproc,myrank)
!--------------------------------------------
!
! Calculate ECP integrals and add them into Hamiltonian matrix
!
! Out  : hstmat2 (one electron Hamiltonian matrix)
!
      use modecp, only : nterm1, nterm2, maxangecp
      use modbasis, only : nao, mtype, nshell
      use modmolecule, only : natom
      implicit none
      integer,intent(in) :: nproc, myrank
      integer :: maxfunc(0:6)=(/1,3,6,10,15,21,28/)
      integer :: maxbasis, numtbasis, maxdim, llmax, maxecpdim, nsizecp1, nsizecp2, ish, jsh
      integer,allocatable :: label1ecp(:), num1ecp(:), label2ecp(:), num2ecp(:)
      real(8),intent(inout) :: hstmat2((nao*(nao+1))/2)
      real(8),allocatable :: term1ecp(:), term2ecp(:), term0ecp(:), xyzintecp(:)
!
      maxbasis= maxval(mtype(1:nshell))
      maxdim= maxfunc(maxbasis)
      llmax= maxval(maxangecp(1:natom))
      if(maxbasis >= 6) then
        write(*,'(" Error! This program supports up to h function in ecp calculation.")')
        call exit
      endif
      if(llmax >= 5) then
        write(*,'(" Error! This program supports up to SPDFG core potentials.")')
        call exit
      endif
      maxecpdim= max(maxbasis,llmax-1)
      numtbasis=(maxecpdim+1)*(maxecpdim+2)*(maxecpdim+3)/6
      nsizecp1=(numtbasis*numtbasis+numtbasis)/2+1
      nsizecp2= llmax*llmax*numtbasis+1
!
      call memset(nterm1*10+nterm2*7+nsizecp1+nsizecp2*2+25*25*25)
      allocate(term1ecp(nterm1),term2ecp(nterm2),term0ecp(nsizecp2),xyzintecp(25*25*25), &
&              label1ecp(nterm1*9),label2ecp(nterm2*6),num1ecp(nsizecp1),num2ecp(nsizecp2))
!
! Prepare intermediate integrals for ECP
!
      call ecpinit(term1ecp,term2ecp,term0ecp,xyzintecp,label1ecp,label2ecp, &
&                  num1ecp,num2ecp,maxecpdim,llmax,numtbasis)
!
!$OMP parallel
      do ish= nshell-myrank,1,-nproc
!$OMP do
        do jsh= 1,ish
          call calcintecp(hstmat2,term1ecp,term2ecp,term0ecp,xyzintecp,label1ecp,label2ecp, &
&                         num1ecp,num2ecp,numtbasis,ish,jsh,maxdim)
        enddo
!$OMP enddo
      enddo
!$OMP end parallel
!
      deallocate(term1ecp,term2ecp,term0ecp,xyzintecp, &
&                label1ecp,label2ecp,num1ecp,num2ecp)
      call memunset(nterm1*10+nterm2*7+nsizecp1+nsizecp2*2+25*25*25)
!
      return
end


!---------------------------------------------------------------------------------
  subroutine ecpinit(term1ecp,term2ecp,term0ecp,xyzintecp,label1ecp,label2ecp, &
&                    num1ecp,num2ecp,maxecpdim,llmax,numtbasis)
!---------------------------------------------------------------------------------
!
! Prepare values for ECP calculation
!
      use modecp, only : nterm1, nterm2
      implicit none
      integer,parameter :: lenecp(0:6)=(/1,4,10,20,35,56,84/)
      integer,intent(in) :: maxecpdim, llmax, numtbasis
      integer,intent(out) :: label1ecp(9*nterm1), label2ecp(6*nterm2)
      integer,intent(out) :: num1ecp(*), num2ecp(*)
      real(8),intent(out) :: term0ecp(*)
      real(8),intent(out) :: term1ecp(nterm1), term2ecp(nterm2), xyzintecp(25*25*25)
!
      call ecpzlm
      call ecpxyz(xyzintecp,maxecpdim)
      call ecpangint0(term0ecp,xyzintecp,maxecpdim,llmax,numtbasis)
      call ecpangint1(term1ecp,label1ecp,num1ecp,xyzintecp,maxecpdim)
      call ecpangint2(term2ecp,label2ecp,num2ecp,xyzintecp,maxecpdim,llmax,numtbasis)
!
      return
end


!-----------------------------------------------------------------------------------------
  subroutine calcintecp(hmat,term1ecp,term2ecp,term0ecp,xyzintecp,label1ecp,label2ecp, &
&                       num1ecp,num2ecp,numtbasis,ish,jsh,len1)
!-----------------------------------------------------------------------------------------
!
! Driver of effective core potential integrals
!   (j|Uecp|i)
!
      use modparam, only : mxprsh
      use modmolecule, only : natom, coord
      use modbasis, only : locatom, locprim, locbf, mprim, mbf, mtype, ex, coeff, nao
      use modecp, only : maxangecp, mtypeecp, locecp, mprimecp, execp, coeffecp, nterm1, nterm2
      implicit none
      integer,intent(in) :: label1ecp(nterm1*9), label2ecp(nterm2*6), num1ecp(*), num2ecp(*)
      integer,intent(in) :: numtbasis, ish, jsh, len1
      integer :: nangij(2), nprimij(2), nbfij(2), iatom, jatom, katom
      integer :: iloc, jloc, ilocbf, jlocbf, iprim, jprim, i, j, ii, ij, maxj
      integer :: ijkatom(3), nangkecp(mxprsh,6), nprimkecp(6), lmaxecp
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: term1ecp(nterm1), term2ecp(nterm2), term0ecp(*), xyzintecp(25*25*25)
      real(8),intent(inout) :: hmat((nao*(nao+1))/2)
      real(8) :: exij(mxprsh,2), coij(mxprsh,2), coordijk(3,3), ecpint(len1,len1)
      real(8) :: exkecp(mxprsh,6), cokecp(mxprsh,6), ecpintsum(len1,len1)
      logical :: iandj
!
      iandj=(ish == jsh)
      nangij(1)= mtype(ish)
      nangij(2)= mtype(jsh)
      nprimij(1)= mprim(ish)
      nprimij(2)= mprim(jsh)
      nbfij(1)  = mbf(ish)
      nbfij(2)  = mbf(jsh)
      iatom = locatom(ish)
      iloc  = locprim(ish)
      ilocbf= locbf(ish)
      jatom = locatom(jsh)
      jloc  = locprim(jsh)
      jlocbf= locbf(jsh)
!
      if((nangij(1) > 6).or.(nangij(2) > 6))then
        write(*,'(" Error! This program supports up to i function in calcintecp")')
        call exit
      endif
!
      do i= 1,3
        coordijk(i,1)= coord(i,iatom)
        coordijk(i,2)= coord(i,jatom)
      enddo
      do iprim= 1,nprimij(1)
        exij(iprim,1)= ex(iloc+iprim)
        coij(iprim,1)= coeff(iloc+iprim)
      enddo
      do jprim= 1,nprimij(2)
        exij(jprim,2)= ex(jloc+jprim)
        coij(jprim,2)= coeff(jloc+jprim)
      enddo
      ijkatom(1)= iatom
      ijkatom(2)= jatom
!
      ecpintsum(1:nbfij(2),1:nbfij(1))= zero
      do katom= 1,natom
        if(maxangecp(katom) == -1) cycle
        do i= 1,3
          coordijk(i,3)= coord(i,katom)
        enddo
        ijkatom(3)= katom
        lmaxecp= maxangecp(katom)
        do i= 1,maxangecp(katom)+1
          nprimkecp(i)= mprimecp(i-1,katom)
          do j= 1,nprimkecp(i)
            exkecp(j,i)= execp(locecp(i-1,katom)+j)
            cokecp(j,i)= coeffecp(locecp(i-1,katom)+j)
            nangkecp(j,i)= mtypeecp(locecp(i-1,katom)+j)
          enddo
        enddo
        call intecp(ecpint,exij,coij,coordijk,nprimij,nangij,nbfij, &
&                   exkecp,cokecp,nprimkecp,nangkecp,lmaxecp, &
&                   term1ecp,term2ecp,term0ecp,xyzintecp,label1ecp,label2ecp, &
&                   num1ecp,num2ecp,numtbasis,ijkatom,len1,mxprsh)
        do i= 1,nbfij(1)
          do j= 1,nbfij(2)
            ecpintsum(j,i)= ecpintsum(j,i)+ecpint(j,i)
          enddo
        enddo
      enddo
!
      maxj= nbfij(2)
      do i= 1,nbfij(1)
        if(iandj) maxj= i
        ii= ilocbf+i
        ij= ii*(ii-1)/2+jlocbf
        do j= 1,maxj
          hmat(ij+j)= hmat(ij+j)+ecpintsum(j,i)
        enddo
      enddo
!
      return
end


!--------------------------------------------------------------------------------
  subroutine intecp(ecpint,exij,coij,coordijk,nprimij,nangij,nbfij, &
&                   exkecp,cokecp,nprimkecp,nangkecp,lmaxecp, &
&                   term1ecp,term2ecp,term0ecp,xyzintecp,label1ecp,label2ecp, &
&                   num1ecp,num2ecp,numtbasis,ijkatom,len1,mxprsh)
!---------------------------------------------------------------------------------
!
! Wrapper of ECP integral routine 
!
      implicit none
      integer,intent(in) :: mxprsh, nprimij(2), nangij(2), nbfij(2), nprimkecp(6)
      integer,intent(in) :: nangkecp(mxprsh,6), lmaxecp, label1ecp(*), label2ecp(*)
      integer,intent(in) :: num1ecp(*), num2ecp(*), numtbasis, ijkatom(3), len1
      integer :: icab
      real(8),intent(in) :: exij(mxprsh,2), coij(mxprsh,2), coordijk(3,3)
      real(8),intent(in) :: exkecp(mxprsh,6), cokecp(mxprsh,6), term1ecp(*), term2ecp(*)
      real(8),intent(in) :: term0ecp(*), xyzintecp(25*25*25)
      real(8),intent(inout) :: ecpint(len1,len1)
!
      if(ijkatom(1) == ijkatom(3)) then
        if(ijkatom(2) == ijkatom(3)) then
          call ecpaaa(ecpint,exij,coij,nprimij,nangij,nbfij, &
&                     exkecp,cokecp,nprimkecp,nangkecp,lmaxecp, &
&                     term0ecp,xyzintecp,numtbasis,len1,mxprsh)
        else
          icab= 2
          call ecpcaa(ecpint,exij,coij,coordijk,nprimij,nangij,nbfij, &
&                     exkecp,cokecp,nprimkecp,nangkecp,lmaxecp, &
&                     term1ecp,term2ecp,term0ecp,label1ecp,label2ecp, &
&                     num1ecp,num2ecp,numtbasis,len1,mxprsh,icab)
        endif
      else
        if(ijkatom(2) == ijkatom(3)) then
          icab= 3
          call ecpcaa(ecpint,exij,coij,coordijk,nprimij,nangij,nbfij, &
&                     exkecp,cokecp,nprimkecp,nangkecp,lmaxecp, &
&                     term1ecp,term2ecp,term0ecp,label1ecp,label2ecp, &
&                     num1ecp,num2ecp,numtbasis,len1,mxprsh,icab)
        else
          call ecpcab(ecpint,exij,coij,coordijk,nprimij,nangij,nbfij, &
&                     exkecp,cokecp,nprimkecp,nangkecp,lmaxecp, &
&                     term1ecp,term2ecp,label1ecp,label2ecp, &
&                     num1ecp,num2ecp,numtbasis,len1,mxprsh)
        endif
      endif
!
      return
end


!----------------------------------------------------------------
  subroutine ecpaaa(ecpint,exij,coij,nprimij,nangij,nbfij, &
&                   exkecp,cokecp,nprimkecp,nangkecp,lmaxecp, &
&                   term0ecp,xyzintecp,numtbasis,len1,mxprsh)
!----------------------------------------------------------------
!
! Calculate  <A|A|A>-type ECP integral
!
      use modecp, only : nx, ny, nz
      implicit none
      integer :: locxyzecp(0:7)=(/0,1,4,10,20,35,56,84/)
      integer :: ncart(0:6)=(/1,3,6,10,15,21,28/)
      integer,intent(in) :: len1, mxprsh, nprimij(2), nangij(2), nbfij(2), nprimkecp(6)
      integer,intent(in) :: nangkecp(mxprsh,6), lmaxecp, numtbasis
      integer :: nang12, ncarti, ncartj, iprim, jprim, kk, iindx, jindx, ll, lambda, nlm1
      integer :: ii, jj, i2, j2, mx, my, mz, mu, nlm
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00, four=4.0D+00
      real(8),parameter :: five=5.0D+00, six=6.0D+00, eight=8.0D+00, p9=9.0D+00, ten=10.0D+00
      real(8),parameter :: twelve=12.0D+00, p15=15.0D+00, p16=16.0D+00
      real(8),parameter :: p20=20.0D+00, p24=24.0D+00, p30=30.0D+00, p32=32.0D+00
      real(8),parameter :: p40=40.0D+00, p60=60.0D+00, p90=90.0D+00, p120=1.2D+02
      real(8),parameter :: p180=1.8D+02, eighth=0.125D+00, sixteenth=6.25D-02
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
      real(8),intent(in) :: exij(mxprsh,2), coij(mxprsh,2),exkecp(mxprsh,6)
      real(8),intent(in) :: cokecp(mxprsh,6), term0ecp(*), xyzintecp(0:24,0:24,0:24)
      real(8),intent(out) :: ecpint(len1,len1)
      real(8) :: ecptmp(28,28), work(28), rad1, ex12, co12, tmp1, rad2, tmp2, type1da0, ang2
!
      nang12= nangij(1)+nangij(2)
      if(mod(nang12,2) == 1) then
        ecpint(1:nbfij(2),1:nbfij(1))= zero
        return
      endif
!
      ncarti= ncart(nangij(1))
      ncartj= ncart(nangij(2))
      ecptmp(1:ncartj,1:ncarti)= zero
!
! Calculate radial integral of type1
!
      rad1= zero
      do iprim= 1,nprimij(1)
        do jprim=1,nprimij(2)
          ex12= exij(iprim,1)+exij(jprim,2)
          co12= coij(iprim,1)*coij(jprim,2)
          tmp1= zero
          do kk=1,nprimkecp(1)
            tmp1= tmp1+type1da0(nang12+nangkecp(kk,1),ex12+exkecp(kk,1))*cokecp(kk,1)
          enddo
          rad1= rad1+tmp1*co12
        enddo
      enddo
!
! Calculate type1 term
!
      iindx= locxyzecp(nangij(1))
      jindx= locxyzecp(nangij(2))
      do ii= 1,ncarti
        do jj= 1,ncartj
          mx= nx(iindx+ii)+nx(jindx+jj)
          my= ny(iindx+ii)+ny(jindx+jj)
          mz= nz(iindx+ii)+nz(jindx+jj)
          ecptmp(jj,ii)= ecptmp(jj,ii)+rad1*xyzintecp(mz,my,mx)
        enddo
      enddo
!
      do ll= 2,lmaxecp+1
        lambda= ll-2
        nlm1=(lambda*(lambda+1))
!
! Calculate radial integral of type2
!
        rad2= zero
        do iprim= 1,nprimij(1)
          do jprim= 1,nprimij(2)
            ex12= exij(iprim,1)+exij(jprim,2)
            co12= coij(iprim,1)*coij(jprim,2)
            tmp2= zero
            do kk= 1,nprimkecp(ll)
             tmp2= tmp2+type1da0(nang12+nangkecp(kk,ll),ex12+exkecp(kk,ll))*cokecp(kk,ll)
            enddo
            rad2 = rad2+tmp2*co12
          enddo
        enddo
!
! Calculate type2 term
!
        do ii= 1,ncarti
          i2= iindx+ii
          do jj= 1,ncartj
            j2= jindx+jj
            ang2= zero
            do mu=-lambda,lambda
              nlm=(nlm1-mu)*numtbasis
              ang2= ang2+term0ecp(nlm+i2)*term0ecp(nlm+j2)
            enddo
            ecptmp(jj,ii)= ecptmp(jj,ii)+rad2*ang2
          enddo
        enddo
      enddo
!
! Normalize ECP integrals
!
! Bra part
!
      select case(nbfij(2))
! D function
        case(5)
          do ii= 1,ncarti
            do jj= 1,6
              work(jj)= ecptmp(jj,ii)
            enddo
            ecptmp(1,ii)= work(4)*sqrt3
            ecptmp(2,ii)= work(6)*sqrt3
            ecptmp(3,ii)=(work(3)*two-work(1)-work(2))*half
            ecptmp(4,ii)= work(5)*sqrt3
            ecptmp(5,ii)=(work(1)-work(2))*sqrt3h
          enddo
        case(6)
          do ii= 1,ncarti
            work(1:6)= ecptmp(1:6,ii)
            ecptmp(1,ii)= work(1)
            ecptmp(2,ii)= work(4)*sqrt3
            ecptmp(3,ii)= work(5)*sqrt3
            ecptmp(4,ii)= work(2)
            ecptmp(5,ii)= work(6)*sqrt3
            ecptmp(6,ii)= work(3)
          enddo
! F function
        case(7)
          do ii= 1,ncarti
            do jj= 1,10
              work(jj)= ecptmp(jj,ii)
            enddo
            ecptmp(1,ii)=(-work(2)+three*work(4)                  )*facf1
            ecptmp(2,ii)=  work(10)                                *facf2
            ecptmp(3,ii)=(-work(2)-work(4)+four*work(9)           )*facf3
            ecptmp(4,ii)=( two*work(3)-three*work(5)-three*work(7))*half
            ecptmp(5,ii)=(-work(1)-work(6)+four*work(8)           )*facf3
            ecptmp(6,ii)=( work(5)-work(7)                        )*facf4
            ecptmp(7,ii)=( work(1)-three*work(6)                  )*facf1
          enddo
        case(10)
          do ii= 1,ncarti
            work(1:10)= ecptmp(1:10,ii)
            ecptmp( 1,ii)= work( 1)
            ecptmp( 2,ii)= work( 4)*sqrt5
            ecptmp( 3,ii)= work( 5)*sqrt5
            ecptmp( 4,ii)= work( 6)*sqrt5
            ecptmp( 5,ii)= work(10)*sqrt15
            ecptmp( 6,ii)= work( 8)*sqrt5
            ecptmp( 7,ii)= work( 2)
            ecptmp( 8,ii)= work( 7)*sqrt5
            ecptmp( 9,ii)= work( 9)*sqrt5
            ecptmp(10,ii)= work( 3)
          enddo
! G function
        case(9)
          do ii= 1,ncarti
            do jj= 1,15
              work(jj)= ecptmp(jj,ii)
            enddo
            ecptmp(1,ii)=(work(4)-work(6))*facg1
            ecptmp(2,ii)=(-work(7)+work(13)*three)*facg2
            ecptmp(3,ii)=(-work(4)-work(6)+work(15)*six)*facg3
            ecptmp(4,ii)=(-work(7)*three+work(9)*four-work(13)*three)*facg4
            ecptmp(5,ii)=(work(1)*three+work(2)*three+work(3)*eight+work(10)*six &
&                        -work(11)*p24-work(12)*p24)*eighth
            ecptmp(6,ii)=(-work(5)*three+work(8)*four-work(14)*three)*facg4
            ecptmp(7,ii)=(-work(1)+work(2)+work(11)*six-work(12)*six)*facg5
            ecptmp(8,ii)=(work(5)-work(14)*three)*facg2
            ecptmp(9,ii)=(work(1)+work(2)-work(10)*six)*facg6
          enddo
        case(15)
          do ii= 1,ncarti
            work(1:15)= ecptmp(1:15,ii)
            ecptmp( 1,ii)= work( 1)
            ecptmp( 2,ii)= work( 4)*sqrt7
            ecptmp( 3,ii)= work( 5)*sqrt7
            ecptmp( 4,ii)= work(10)*sqrt35third
            ecptmp( 5,ii)= work(13)*sqrt35
            ecptmp( 6,ii)= work(11)*sqrt35third
            ecptmp( 7,ii)= work( 6)*sqrt7
            ecptmp( 8,ii)= work(14)*sqrt35
            ecptmp( 9,ii)= work(15)*sqrt35
            ecptmp(10,ii)= work( 8)*sqrt7
            ecptmp(11,ii)= work( 2)
            ecptmp(12,ii)= work( 7)*sqrt7
            ecptmp(13,ii)= work(12)*sqrt35third
            ecptmp(14,ii)= work( 9)*sqrt7
            ecptmp(15,ii)= work( 3)
          enddo
! H function
        case(11)
          do ii= 1,ncarti
            work(1:21)= ecptmp(1:21,ii)
            ecptmp( 1,ii)=(work(4)*five-work(12)*ten+work(2))*fach1
            ecptmp( 2,ii)=(work(16)*four-work(17)*four)*fach2
            ecptmp( 3,ii)=(-work(4)*three-work(12)*two+work(20)*p24+work(2)-work(13)*eight)*fach3
            ecptmp( 4,ii)=(-work(16)*two-work(17)*two+work(18)*four)*fach4
            ecptmp( 5,ii)=(work(4)+work(12)*two-work(20)*twelve+work(2)-work(13)*twelve &
&                         +work(9)*eight)*fach5
            ecptmp( 6,ii)=(work(5)*p15+work(19)*p30-work(14)*p40+work(7)*p15-work(15)*p40 &
&                         +work(3)*eight)*eighth
            ecptmp( 7,ii)=(work(1)+work(10)*two-work(11)*twelve+work(6)-work(21)*twelve &
&                         +work(8)*eight)*fach5
            ecptmp( 8,ii)=(-work(5)+work(14)*two+work(7)-work(15)*two)*fach4
            ecptmp( 9,ii)=(-work(1)+work(10)*two+work(11)*eight+work(6)*three-work(21)*p24)*fach3
            ecptmp(10,ii)=(work(5)-work(19)*six+work(7))*fach2
            ecptmp(11,ii)=(work(1)-work(10)*ten+work(6)*five)*fach1
          enddo
        case(21)
          do ii= 1,ncarti
            work(1:21)= ecptmp(1:21,ii)
            ecptmp( 1,ii)= work( 1)
            ecptmp( 2,ii)= work( 4)*three
            ecptmp( 3,ii)= work( 5)*three
            ecptmp( 4,ii)= work(10)*sqrt21
            ecptmp( 5,ii)= work(16)*sqrt63
            ecptmp( 6,ii)= work(11)*sqrt21
            ecptmp( 7,ii)= work(12)*sqrt21
            ecptmp( 8,ii)= work(19)*sqrt105
            ecptmp( 9,ii)= work(20)*sqrt105
            ecptmp(10,ii)= work(14)*sqrt21
            ecptmp(11,ii)= work( 6)*three
            ecptmp(12,ii)= work(17)*sqrt63
            ecptmp(13,ii)= work(21)*sqrt105
            ecptmp(14,ii)= work(18)*sqrt63
            ecptmp(15,ii)= work( 8)*three
            ecptmp(16,ii)= work( 2)
            ecptmp(17,ii)= work( 7)*three
            ecptmp(18,ii)= work(13)*sqrt21
            ecptmp(19,ii)= work(15)*sqrt21
            ecptmp(20,ii)= work( 9)*three
            ecptmp(21,ii)= work( 3)
          enddo
! I function
        case(13)
          do ii= 1,ncarti
            work(1:28)= ecptmp(1:28,ii)
            ecptmp( 1,ii)=(work(4)*six-work(19)*p20+work(6)*six)*faci1
            ecptmp( 2,ii)=(work(16)*five-work(24)*ten+work(7))*faci2
            ecptmp( 3,ii)=(-work(4)*four+work(23)*p40+work(6)*four-work(25)*p40)*faci3
            ecptmp( 4,ii)=(-work(16)*p9-work(24)*six+work(26)*p24+work(7)*three &
&                         -work(21)*eight)*faci4
            ecptmp( 5,ii)=(work(4)*two+work(19)*four-work(23)*p32+work(6)*two &
&                         -work(25)*p32+work(18)*p32)*faci5
            ecptmp( 6,ii)=(work(16)*five+work(24)*ten-work(26)*p20+work(7)*five &
&                         -work(21)*p20+work(9)*eight)*faci6
            ecptmp( 7,ii)=(-work(1)*five-work(10)*p15+work(11)*p90-work(12)*p15 &
&                         +work(28)*p180-work(14)*p120-work(2)*five+work(13)*p90 &
&                         -work(15)*p120+work(3)*p16)*sixteenth
            ecptmp( 8,ii)=(work(5)*five+work(22)*ten-work(20)*p20+work(17)*five &
&                         -work(27)*p20+work(8)*eight)*faci6
            ecptmp( 9,ii)=(work(1)+work(10)-work(11)*p16-work(12)+work(14)*p16 &
&                         -work(2)+work(13)*p16-work(15)*p16)*faci5
            ecptmp(10,ii)=(-work(5)*three+work(22)*six+work(20)*eight+work(17)*p9 &
&                         -work(27)*p24)*faci4
            ecptmp(11,ii)=(-work(1)+work(10)*five+work(11)*ten+work(12)*five &
&                         -work(28)*p60-work(2)+work(13)*ten)*faci3
            ecptmp(12,ii)=(work(5)-work(22)*ten+work(17)*five)*faci2
            ecptmp(13,ii)=(work(1)-work(10)*p15+work(12)*p15-work(2))*faci1
          enddo
        case(28)
          do ii= 1,ncarti
            work(1:28)= ecptmp(1:28,ii)
            ecptmp( 1,ii)= work( 1)
            ecptmp( 2,ii)= work( 4)*sqrt11
            ecptmp( 3,ii)= work( 5)*sqrt11
            ecptmp( 4,ii)= work(10)*sqrt33
            ecptmp( 5,ii)= work(16)*sqrt99
            ecptmp( 6,ii)= work(11)*sqrt33
            ecptmp( 7,ii)= work(19)*sqrt231fifth
            ecptmp( 8,ii)= work(22)*sqrt231
            ecptmp( 9,ii)= work(23)*sqrt231
            ecptmp(10,ii)= work(20)*sqrt231fifth
            ecptmp(11,ii)= work(12)*sqrt33
            ecptmp(12,ii)= work(24)*sqrt231
            ecptmp(13,ii)= work(28)*sqrt385
            ecptmp(14,ii)= work(26)*sqrt231
            ecptmp(15,ii)= work(14)*sqrt33
            ecptmp(16,ii)= work( 6)*sqrt11
            ecptmp(17,ii)= work(17)*sqrt99
            ecptmp(18,ii)= work(25)*sqrt231
            ecptmp(19,ii)= work(27)*sqrt231
            ecptmp(20,ii)= work(18)*sqrt99
            ecptmp(21,ii)= work( 8)*sqrt11
            ecptmp(22,ii)= work( 2)
            ecptmp(23,ii)= work( 7)*sqrt11
            ecptmp(24,ii)= work(13)*sqrt33
            ecptmp(25,ii)= work(21)*sqrt231fifth
            ecptmp(26,ii)= work(15)*sqrt33
            ecptmp(27,ii)= work( 9)*sqrt11
            ecptmp(28,ii)= work( 3)
          enddo
      end select
!
! Ket part
!
      select case(nbfij(1))
! D function
        case(5)
          do jj= 1,nbfij(2)
            do ii= 1,6
              work(ii)= ecptmp(jj,ii)
            enddo
            ecptmp(jj,1)= work(4)*sqrt3
            ecptmp(jj,2)= work(6)*sqrt3
            ecptmp(jj,3)=(work(3)*two-work(1)-work(2))*half
            ecptmp(jj,4)= work(5)*sqrt3
            ecptmp(jj,5)=(work(1)-work(2))*sqrt3h
          enddo
        case(6)
          do jj= 1,nbfij(2)
            work(1:6)= ecptmp(jj,1:6)
            ecptmp(jj,1)= work(1)
            ecptmp(jj,2)= work(4)*sqrt3
            ecptmp(jj,3)= work(5)*sqrt3
            ecptmp(jj,4)= work(2)
            ecptmp(jj,5)= work(6)*sqrt3
            ecptmp(jj,6)= work(3)
          enddo
! F function
        case(7)
          do jj= 1,nbfij(2)
            do ii= 1,10
              work(ii)= ecptmp(jj,ii)
            enddo
            ecptmp(jj,1)=(-work(2)+three*work(4)                  )*facf1
            ecptmp(jj,2)=  work(10)                                *facf2
            ecptmp(jj,3)=(-work(2)-work(4)+four*work(9)           )*facf3
            ecptmp(jj,4)=( two*work(3)-three*work(5)-three*work(7))*half
            ecptmp(jj,5)=(-work(1)-work(6)+four*work(8)           )*facf3
            ecptmp(jj,6)=( work(5)-work(7)                        )*facf4
            ecptmp(jj,7)=( work(1)-three*work(6)                  )*facf1
          enddo
        case(10)
          do jj= 1,nbfij(2)
            work(1:10)= ecptmp(jj,1:10)
            ecptmp(jj, 1)= work( 1)
            ecptmp(jj, 2)= work( 4)*sqrt5
            ecptmp(jj, 3)= work( 5)*sqrt5
            ecptmp(jj, 4)= work( 6)*sqrt5
            ecptmp(jj, 5)= work(10)*sqrt15
            ecptmp(jj, 6)= work( 8)*sqrt5
            ecptmp(jj, 7)= work( 2)
            ecptmp(jj, 8)= work( 7)*sqrt5
            ecptmp(jj, 9)= work( 9)*sqrt5
            ecptmp(jj,10)= work( 3)
          enddo
! G function
        case(9)
          do jj= 1,nbfij(2)
            do ii= 1,15
              work(ii)= ecptmp(jj,ii)
            enddo
            ecptmp(jj,1)=(work(4)-work(6))*facg1
            ecptmp(jj,2)=(-work(7)+work(13)*three)*facg2
            ecptmp(jj,3)=(-work(4)-work(6)+work(15)*six)*facg3
            ecptmp(jj,4)=(-work(7)*three+work(9)*four-work(13)*three)*facg4
            ecptmp(jj,5)=(work(1)*three+work(2)*three+work(3)*eight+work(10)*six &
&                       -work(11)*p24-work(12)*p24)*eighth
            ecptmp(jj,6)=(-work(5)*three+work(8)*four-work(14)*three)*facg4
            ecptmp(jj,7)=(-work(1)+work(2)+work(11)*six-work(12)*six)*facg5
            ecptmp(jj,8)=(work(5)-work(14)*three)*facg2
            ecptmp(jj,9)=(work(1)+work(2)-work(10)*six)*facg6
          enddo
        case(15)
          do jj= 1,nbfij(2)
            work(1:15)= ecptmp(jj,1:15)
            ecptmp(jj, 1)= work( 1)
            ecptmp(jj, 2)= work( 4)*sqrt7
            ecptmp(jj, 3)= work( 5)*sqrt7
            ecptmp(jj, 4)= work(10)*sqrt35third
            ecptmp(jj, 5)= work(13)*sqrt35
            ecptmp(jj, 6)= work(11)*sqrt35third
            ecptmp(jj, 7)= work( 6)*sqrt7
            ecptmp(jj, 8)= work(14)*sqrt35
            ecptmp(jj, 9)= work(15)*sqrt35
            ecptmp(jj,10)= work( 8)*sqrt7
            ecptmp(jj,11)= work( 2)
            ecptmp(jj,12)= work( 7)*sqrt7
            ecptmp(jj,13)= work(12)*sqrt35third
            ecptmp(jj,14)= work( 9)*sqrt7
            ecptmp(jj,15)= work( 3)
          enddo
! H function
        case(11)
          do jj= 1,nbfij(2)
            work(1:21)= ecptmp(jj,1:21)
            ecptmp(jj, 1)=(work(4)*five-work(12)*ten+work(2))*fach1
            ecptmp(jj, 2)=(work(16)*four-work(17)*four)*fach2
            ecptmp(jj, 3)=(-work(4)*three-work(12)*two+work(20)*p24+work(2)-work(13)*eight)*fach3
            ecptmp(jj, 4)=(-work(16)*two-work(17)*two+work(18)*four)*fach4
            ecptmp(jj, 5)=(work(4)+work(12)*two-work(20)*twelve+work(2)-work(13)*twelve &
&                         +work(9)*eight)*fach5
            ecptmp(jj, 6)=(work(5)*p15+work(19)*p30-work(14)*p40+work(7)*p15-work(15)*p40 &
&                         +work(3)*eight)*eighth
            ecptmp(jj, 7)=(work(1)+work(10)*two-work(11)*twelve+work(6)-work(21)*twelve &
&                         +work(8)*eight)*fach5
            ecptmp(jj, 8)=(-work(5)+work(14)*two+work(7)-work(15)*two)*fach4
            ecptmp(jj, 9)=(-work(1)+work(10)*two+work(11)*eight+work(6)*three-work(21)*p24)*fach3
            ecptmp(jj,10)=(work(5)-work(19)*six+work(7))*fach2
            ecptmp(jj,11)=(work(1)-work(10)*ten+work(6)*five)*fach1
          enddo
        case(21)
          do jj= 1,nbfij(2)
            work(1:21)= ecptmp(jj,1:21)
            ecptmp(jj, 1)= work( 1)
            ecptmp(jj, 2)= work( 4)*three
            ecptmp(jj, 3)= work( 5)*three
            ecptmp(jj, 4)= work(10)*sqrt21
            ecptmp(jj, 5)= work(16)*sqrt63
            ecptmp(jj, 6)= work(11)*sqrt21
            ecptmp(jj, 7)= work(12)*sqrt21
            ecptmp(jj, 8)= work(19)*sqrt105
            ecptmp(jj, 9)= work(20)*sqrt105
            ecptmp(jj,10)= work(14)*sqrt21
            ecptmp(jj,11)= work( 6)*three
            ecptmp(jj,12)= work(17)*sqrt63
            ecptmp(jj,13)= work(21)*sqrt105
            ecptmp(jj,14)= work(18)*sqrt63
            ecptmp(jj,15)= work( 8)*three
            ecptmp(jj,16)= work( 2)
            ecptmp(jj,17)= work( 7)*three
            ecptmp(jj,18)= work(13)*sqrt21
            ecptmp(jj,19)= work(15)*sqrt21
            ecptmp(jj,20)= work( 9)*three
            ecptmp(jj,21)= work( 3)
          enddo
! I function
        case(13)
          do jj= 1,nbfij(2)
            work(1:28)= ecptmp(jj,1:28)
            ecptmp(jj, 1)=(work(4)*six-work(19)*p20+work(6)*six)*faci1
            ecptmp(jj, 2)=(work(16)*five-work(24)*ten+work(7))*faci2
            ecptmp(jj, 3)=(-work(4)*four+work(23)*p40+work(6)*four-work(25)*p40)*faci3
            ecptmp(jj, 4)=(-work(16)*p9-work(24)*six+work(26)*p24+work(7)*three &
&                         -work(21)*eight)*faci4
            ecptmp(jj, 5)=(work(4)*two+work(19)*four-work(23)*p32+work(6)*two &
&                         -work(25)*p32+work(18)*p32)*faci5
            ecptmp(jj, 6)=(work(16)*five+work(24)*ten-work(26)*p20+work(7)*five &
&                         -work(21)*p20+work(9)*eight)*faci6
            ecptmp(jj, 7)=(-work(1)*five-work(10)*p15+work(11)*p90-work(12)*p15 &
&                         +work(28)*p180-work(14)*p120-work(2)*five+work(13)*p90 &
&                         -work(15)*p120+work(3)*p16)*sixteenth
            ecptmp(jj, 8)=(work(5)*five+work(22)*ten-work(20)*p20+work(17)*five &
&                         -work(27)*p20+work(8)*eight)*faci6
            ecptmp(jj, 9)=(work(1)+work(10)-work(11)*p16-work(12)+work(14)*p16 &
&                         -work(2)+work(13)*p16-work(15)*p16)*faci5
            ecptmp(jj,10)=(-work(5)*three+work(22)*six+work(20)*eight+work(17)*p9 &
&                         -work(27)*p24)*faci4
            ecptmp(jj,11)=(-work(1)+work(10)*five+work(11)*ten+work(12)*five &
&                         -work(28)*p60-work(2)+work(13)*ten)*faci3
            ecptmp(jj,12)=(work(5)-work(22)*ten+work(17)*five)*faci2
            ecptmp(jj,13)=(work(1)-work(10)*p15+work(12)*p15-work(2))*faci1
          enddo
        case(28)
          do jj= 1,nbfij(2)
            work(1:28)= ecptmp(jj,1:28)
            ecptmp(jj, 1)= work( 1)
            ecptmp(jj, 2)= work( 4)*sqrt11
            ecptmp(jj, 3)= work( 5)*sqrt11
            ecptmp(jj, 4)= work(10)*sqrt33
            ecptmp(jj, 5)= work(16)*sqrt99
            ecptmp(jj, 6)= work(11)*sqrt33
            ecptmp(jj, 7)= work(19)*sqrt231fifth
            ecptmp(jj, 8)= work(22)*sqrt231
            ecptmp(jj, 9)= work(23)*sqrt231
            ecptmp(jj,10)= work(20)*sqrt231fifth
            ecptmp(jj,11)= work(12)*sqrt33
            ecptmp(jj,12)= work(24)*sqrt231
            ecptmp(jj,13)= work(28)*sqrt385
            ecptmp(jj,14)= work(26)*sqrt231
            ecptmp(jj,15)= work(14)*sqrt33
            ecptmp(jj,16)= work( 6)*sqrt11
            ecptmp(jj,17)= work(17)*sqrt99
            ecptmp(jj,18)= work(25)*sqrt231
            ecptmp(jj,19)= work(27)*sqrt231
            ecptmp(jj,20)= work(18)*sqrt99
            ecptmp(jj,21)= work( 8)*sqrt11
            ecptmp(jj,22)= work( 2)
            ecptmp(jj,23)= work( 7)*sqrt11
            ecptmp(jj,24)= work(13)*sqrt33
            ecptmp(jj,25)= work(21)*sqrt231fifth
            ecptmp(jj,26)= work(15)*sqrt33
            ecptmp(jj,27)= work( 9)*sqrt11
            ecptmp(jj,28)= work( 3)
          enddo
      end select
!
      do ii= 1,nbfij(1)
        do jj= 1,nbfij(2)
          ecpint(jj,ii)= ecptmp(jj,ii)
        enddo
      enddo
!
      return
end


!----------------------------------------------------------------------
  subroutine ecpcaa(ecpint,exij,coij,coordijk,nprimij,nangij,nbfij, &
&                   exkecp,cokecp,nprimkecp,nangkecp,lmaxecp, &
&                   term1ecp,term2ecp,term0ecp,label1ecp,label2ecp, &
&                   num1ecp,num2ecp,numtbasis,len1,mxprsh,icab)
!----------------------------------------------------------------------
!
! Calculate  <C|A|A> or <A|A|B>-type ECP integral
!
      use modecp, only : nx, ny, nz
      implicit none
      integer,intent(in) :: len1, mxprsh, nprimij(2), nangij(2), nbfij(2), nprimkecp(6)
      integer,intent(in) :: nangkecp(mxprsh,6), lmaxecp, numtbasis
      integer,intent(in) :: label1ecp(9,*), label2ecp(6,*), num1ecp(*), num2ecp(*), icab
      integer :: ncart(0:6)=(/1,3,6,10,15,21,28/), locxyzecp(0:7)=(/0,1,4,10,20,35,56,84/)
      integer :: nangbasis, ncarti, ncartj, ii, jj, nangbasis2, iprim, jprim, nangk, nn
      integer :: iindex, jindex, i1, j1, indx, label1f, label1l, lambda, ijang, nxyz
      integer :: kk, mm, icount, ka, la, mu, kx, ky, kz, kxp, kyp, kzp, lmindx, ll, ll2
      integer :: ksum, nangll, mx, my, mz, nlm, label2f, label2l
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, pi4=1.256637061435917D+01
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00, four=4.0D+00
      real(8),parameter :: five=5.0D+00, six=6.0D+00, eight=8.0D+00, p9=9.0D+00, ten=10.0D+00
      real(8),parameter :: twelve=12.0D+00, p15=15.0D+00, p16=16.0D+00
      real(8),parameter :: p20=20.0D+00, p24=24.0D+00, p30=30.0D+00, p32=32.0D+00
      real(8),parameter :: p40=40.0D+00, p60=60.0D+00, p90=90.0D+00, p120=1.2D+02
      real(8),parameter :: p180=1.8D+02, eighth=0.125D+00, sixteenth=6.25D-02
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
      real(8),intent(in) :: exij(mxprsh,2), coij(mxprsh,2), exkecp(mxprsh,6)
      real(8),intent(in) :: cokecp(mxprsh,6), term0ecp(*), term1ecp(*), term2ecp(*)
      real(8),intent(in) :: coordijk(3,3)
      real(8),intent(out) :: ecpint(len1,len1)
      real(8) :: ecptmp(28,28), work(28), unitvec(3), bax, bay, baz, ba, cax, cay, caz, ca
      real(8) :: abx(0:6), aby(0:6), abz(0:6), acx(0:6), acy(0:6), acz(0:6)
      real(8) :: type12rad(78), ex01, ex02, ex12, co12, dazeta, da2zeta, besselarginv
      real(8) :: xp0, sqrtinvzeta, alpha, xp1, co123, radint(78), zzlm(121)
      real(8) :: type12ang(11,11), type1xyz, type1sum, xyz0, type2xyz, type2sum
!
      nangbasis= nangij(1)+nangij(2)
      ncarti= ncart(nangij(1))
      ncartj= ncart(nangij(2))
      ecptmp(1:ncartj,1:ncarti)= zero
!
      if(icab == 2) then
        bax= coordijk(1,2)-coordijk(1,3)
        bay= coordijk(2,2)-coordijk(2,3)
        baz= coordijk(3,2)-coordijk(3,3)
        ba = sqrt(bax*bax+bay*bay+baz*baz)
        unitvec(1)= bax/ba
        unitvec(2)= bay/ba
        unitvec(3)= baz/ba
        abx(0)= one
        aby(0)= one
        abz(0)= one
        do jj= 1,nangij(2)
          abx(jj)=-bax*abx(jj-1)
          aby(jj)=-bay*aby(jj-1)
          abz(jj)=-baz*abz(jj-1)
        enddo
      else
        cax= coordijk(1,1)-coordijk(1,3)
        cay= coordijk(2,1)-coordijk(2,3)
        caz= coordijk(3,1)-coordijk(3,3)
        ca = sqrt(cax*cax+cay*cay+caz*caz)
        unitvec(1)= cax/ca
        unitvec(2)= cay/ca
        unitvec(3)= caz/ca
        acx(0)= one
        acy(0)= one
        acz(0)= one
        do ii= 1,nangij(1)
          acx(ii)=-cax*acx(ii-1)
          acy(ii)=-cay*acy(ii-1)
          acz(ii)=-caz*acz(ii-1)
        enddo
      endif
!
! Type 1
!
      nangbasis2=((nangbasis+1)*(nangbasis+2))/2
      type12rad(1:nangbasis2)= zero
!
! Radial integral of Type 1
!
      do iprim= 1,nprimij(1)
        ex01= exij(iprim,1)
        do jprim= 1,nprimij(2)
          ex02= exij(jprim,2)
          ex12= ex01+ex02
          co12= coij(iprim,1)*coij(jprim,2)
          if(icab == 2) then
            dazeta= ex02*ba
            da2zeta= dazeta*ba
          elseif(icab == 3) then
            dazeta= ex01*ca
            da2zeta= dazeta*ca
          endif
          besselarginv= one/(dazeta+dazeta)
          xp0 = exp(-da2zeta)
!
          do kk= 1,nprimkecp(1)
            sqrtinvzeta = one/sqrt(ex12+exkecp(kk,1))
            alpha= dazeta*sqrtinvzeta
            xp1= exp(-da2zeta+alpha*alpha)
            nangk= nangkecp(kk,1)
            co123= cokecp(kk,1)*co12
            call calctype1rad(radint,besselarginv,alpha,sqrtinvzeta, &
&                             xp0,xp1,nangk,nangbasis)
            do nn= 1,nangbasis2
              type12rad(nn)= type12rad(nn)+radint(nn)*co123
            enddo
          enddo
        enddo
      enddo
!
! Angular integral of Type 1
!
      call calczspher(zzlm,unitvec,nangbasis)
      iindex= locxyzecp(nangij(1))
      jindex= locxyzecp(nangij(2))
      do i1= 1,ncarti
        ii= iindex+i1
        do j1= 1,ncartj
          jj= jindex+j1
          if(ii >= jj) then
            indx= jj+(ii*(ii-1))/2
          else
            indx= ii+(jj*(jj-1))/2
          endif
          label1f= num1ecp(indx)
          label1l= num1ecp(indx+1)-1
          do lambda= 1,nangbasis+1
            do ijang= lambda,nangbasis+1
              type12ang(ijang,lambda)= zero
            enddo
          enddo
          nxyz = nx(ii)+ny(ii)+nz(ii)+nx(jj)+ny(jj)+nz(jj)
          kk= ii
          mm= 3
          if(icab == 2) then
             kk= jj
             mm= 6
          endif
          if(ii >= jj) mm= 9-mm
          do icount= label1f,label1l
            ka = label1ecp(1,icount)+1
            la = label1ecp(2,icount)+1
            mu = label1ecp(3,icount)
            kx = label1ecp(mm+1,icount)+nx(kk)
            ky = label1ecp(mm+2,icount)+ny(kk)
            kz = label1ecp(mm+3,icount)+nz(kk)
            if((kx+ky+kz) /= nxyz) cycle
            kxp= label1ecp(4,icount)+label1ecp(7,icount)
            kyp= label1ecp(5,icount)+label1ecp(8,icount)
            kzp= label1ecp(6,icount)+label1ecp(9,icount)
            lmindx = la*(la-1)-mu+1
            if(icab == 2) then
              type1xyz= abx(kx-kxp)*aby(ky-kyp)*abz(kz-kzp)*zzlm(lmindx)
            else
              type1xyz= acx(kx-kxp)*acy(ky-kyp)*acz(kz-kzp)*zzlm(lmindx)
            endif
            type12ang(ka,la)= type12ang(ka,la)+term1ecp(icount)*type1xyz
          enddo
          type1sum= zero
          nn= 1
          do ijang= 1,nangbasis+1
            do lambda= 1,ijang
              type1sum= type1sum+type12rad(nn)*type12ang(ijang,lambda)
              nn= nn+1
            enddo
          enddo
          ecptmp(j1,i1)= ecptmp(j1,i1)+pi4*type1sum
        enddo
      enddo
!
! Type 2
!
      do ll= 2,lmaxecp+1
!
! Radial integral of Type 2
!
        type12rad(1:nangbasis2)= zero
        do iprim= 1,nprimij(1)
          ex01= exij(iprim,1)
          do jprim= 1,nprimij(2)
            ex02= exij(jprim,2)
            ex12= ex01+ex02
            co12= coij(iprim,1)*coij(jprim,2)
            if(icab == 2) then
              dazeta= ex02*ba
              da2zeta= dazeta*ba
            elseif(icab == 3) then
              dazeta= ex01*ca
              da2zeta= dazeta*ca
            endif
            besselarginv= one/(dazeta*two)
            xp0= exp(-da2zeta)
!
            do kk= 1,nprimkecp(ll)
              sqrtinvzeta= one/sqrt(ex12+exkecp(kk,ll))
              alpha= dazeta*sqrtinvzeta
              xp1= exp(-da2zeta+alpha*alpha)
              nangk= nangkecp(kk,ll)
              co123= cokecp(kk,ll)*co12
              call calctype1rad(radint,besselarginv,alpha,sqrtinvzeta, &
&                               xp0,xp1,nangk,nangbasis)
              do nn= 1,nangbasis2
                type12rad(nn)= type12rad(nn)+radint(nn)*co123
              enddo
            enddo
          enddo
        enddo
!
! Angular integral of Type 2
!
        ll2=(ll-2)*(ll-1)
        call calczspher(zzlm,unitvec,nangbasis)
        do i1= 1,ncarti
          ii= iindex+i1
          do j1= 1,ncartj
            jj= jindex+j1
            do lambda= 1,nangbasis+1
              do ijang= lambda,nangbasis+1
                type12ang(ijang,lambda)= zero
              enddo
            enddo
            mx= nx(ii)+nx(jj)
            my= ny(ii)+ny(jj)
            mz= nz(ii)+nz(jj)
            if(icab == 2) then
              kk= ii
            else
              kk= jj
            endif
            ksum= nx(kk)+ny(kk)+nz(kk)+1
            nangll= ll-2
            do mm= nangll,-nangll,-1
              nlm=(ll2-mm)*numtbasis
              xyz0= term0ecp(nlm+kk)
              nlm= nlm+ii+jj
              label2f = num2ecp(nlm-kk)
              label2l = num2ecp(nlm-kk+1)-1
              do icount= label2f,label2l
                la= label2ecp(1,icount)+1
                ka= label2ecp(2,icount)+ksum
                mu= label2ecp(3,icount)
                kx= label2ecp(4,icount)+nx(kk)
                ky= label2ecp(5,icount)+ny(kk)
                kz= label2ecp(6,icount)+nz(kk)
                lmindx= la*(la-1)-mu+1
                if(icab == 2) then
                  type2xyz= xyz0*abx(mx-kx)*aby(my-ky)*abz(mz-kz)*zzlm(lmindx)
                else
                  type2xyz= xyz0*acx(mx-kx)*acy(my-ky)*acz(mz-kz)*zzlm(lmindx)
                endif
                type12ang(ka,la)= type12ang(ka,la)+term2ecp(icount)*type2xyz
              enddo
            enddo
            type2sum= zero
            nn= 1
            do ijang= 1,nangbasis+1
              do lambda= 1,ijang
                type2sum= type2sum+type12rad(nn)*type12ang(ijang,lambda)
                nn= nn+1
              enddo
            enddo
            ecptmp(j1,i1)= ecptmp(j1,i1)+pi4*type2sum
          enddo
        enddo
      enddo
!
! Normalize ECP integrals
!
! Bra part
!
      select case(nbfij(2))
! D function
        case(5)
          do ii= 1,ncarti
            do jj= 1,6
              work(jj)= ecptmp(jj,ii)
            enddo
            ecptmp(1,ii)= work(4)*sqrt3
            ecptmp(2,ii)= work(6)*sqrt3
            ecptmp(3,ii)=(work(3)*two-work(1)-work(2))*half
            ecptmp(4,ii)= work(5)*sqrt3
            ecptmp(5,ii)=(work(1)-work(2))*sqrt3h
          enddo
        case(6)
          do ii= 1,ncarti
            work(1:6)= ecptmp(1:6,ii)
            ecptmp(1,ii)= work(1)
            ecptmp(2,ii)= work(4)*sqrt3
            ecptmp(3,ii)= work(5)*sqrt3
            ecptmp(4,ii)= work(2)
            ecptmp(5,ii)= work(6)*sqrt3
            ecptmp(6,ii)= work(3)
          enddo
! F function
        case(7)
          do ii= 1,ncarti
            do jj= 1,10
              work(jj)= ecptmp(jj,ii)
            enddo
            ecptmp(1,ii)=(-work(2)+three*work(4)                  )*facf1
            ecptmp(2,ii)=  work(10)                                *facf2
            ecptmp(3,ii)=(-work(2)-work(4)+four*work(9)           )*facf3
            ecptmp(4,ii)=( two*work(3)-three*work(5)-three*work(7))*half
            ecptmp(5,ii)=(-work(1)-work(6)+four*work(8)           )*facf3
            ecptmp(6,ii)=( work(5)-work(7)                        )*facf4
            ecptmp(7,ii)=( work(1)-three*work(6)                  )*facf1
          enddo
        case(10)
          do ii= 1,ncarti
            work(1:10)= ecptmp(1:10,ii)
            ecptmp( 1,ii)= work( 1)
            ecptmp( 2,ii)= work( 4)*sqrt5
            ecptmp( 3,ii)= work( 5)*sqrt5
            ecptmp( 4,ii)= work( 6)*sqrt5
            ecptmp( 5,ii)= work(10)*sqrt15
            ecptmp( 6,ii)= work( 8)*sqrt5
            ecptmp( 7,ii)= work( 2)
            ecptmp( 8,ii)= work( 7)*sqrt5
            ecptmp( 9,ii)= work( 9)*sqrt5
            ecptmp(10,ii)= work( 3)
          enddo
! G function
        case(9)
          do ii= 1,ncarti
            do jj= 1,15
              work(jj)= ecptmp(jj,ii)
            enddo
            ecptmp(1,ii)=(work(4)-work(6))*facg1
            ecptmp(2,ii)=(-work(7)+work(13)*three)*facg2
            ecptmp(3,ii)=(-work(4)-work(6)+work(15)*six)*facg3
            ecptmp(4,ii)=(-work(7)*three+work(9)*four-work(13)*three)*facg4
            ecptmp(5,ii)=(work(1)*three+work(2)*three+work(3)*eight+work(10)*six &
&                        -work(11)*p24-work(12)*p24)*eighth
            ecptmp(6,ii)=(-work(5)*three+work(8)*four-work(14)*three)*facg4
            ecptmp(7,ii)=(-work(1)+work(2)+work(11)*six-work(12)*six)*facg5
            ecptmp(8,ii)=(work(5)-work(14)*three)*facg2
            ecptmp(9,ii)=(work(1)+work(2)-work(10)*six)*facg6
          enddo
        case(15)
          do ii= 1,ncarti
            work(1:15)= ecptmp(1:15,ii)
            ecptmp( 1,ii)= work( 1)
            ecptmp( 2,ii)= work( 4)*sqrt7
            ecptmp( 3,ii)= work( 5)*sqrt7
            ecptmp( 4,ii)= work(10)*sqrt35third
            ecptmp( 5,ii)= work(13)*sqrt35
            ecptmp( 6,ii)= work(11)*sqrt35third
            ecptmp( 7,ii)= work( 6)*sqrt7
            ecptmp( 8,ii)= work(14)*sqrt35
            ecptmp( 9,ii)= work(15)*sqrt35
            ecptmp(10,ii)= work( 8)*sqrt7
            ecptmp(11,ii)= work( 2)
            ecptmp(12,ii)= work( 7)*sqrt7
            ecptmp(13,ii)= work(12)*sqrt35third
            ecptmp(14,ii)= work( 9)*sqrt7
            ecptmp(15,ii)= work( 3)
          enddo
! H function
        case(11)
          do ii= 1,ncarti
            work(1:21)= ecptmp(1:21,ii)
            ecptmp( 1,ii)=(work(4)*five-work(12)*ten+work(2))*fach1
            ecptmp( 2,ii)=(work(16)*four-work(17)*four)*fach2
            ecptmp( 3,ii)=(-work(4)*three-work(12)*two+work(20)*p24+work(2)-work(13)*eight)*fach3
            ecptmp( 4,ii)=(-work(16)*two-work(17)*two+work(18)*four)*fach4
            ecptmp( 5,ii)=(work(4)+work(12)*two-work(20)*twelve+work(2)-work(13)*twelve &
&                         +work(9)*eight)*fach5
            ecptmp( 6,ii)=(work(5)*p15+work(19)*p30-work(14)*p40+work(7)*p15-work(15)*p40 &
&                         +work(3)*eight)*eighth
            ecptmp( 7,ii)=(work(1)+work(10)*two-work(11)*twelve+work(6)-work(21)*twelve &
&                         +work(8)*eight)*fach5
            ecptmp( 8,ii)=(-work(5)+work(14)*two+work(7)-work(15)*two)*fach4
            ecptmp( 9,ii)=(-work(1)+work(10)*two+work(11)*eight+work(6)*three-work(21)*p24)*fach3
            ecptmp(10,ii)=(work(5)-work(19)*six+work(7))*fach2
            ecptmp(11,ii)=(work(1)-work(10)*ten+work(6)*five)*fach1
          enddo
        case(21)
          do ii= 1,ncarti
            work(1:21)= ecptmp(1:21,ii)
            ecptmp( 1,ii)= work( 1)
            ecptmp( 2,ii)= work( 4)*three
            ecptmp( 3,ii)= work( 5)*three
            ecptmp( 4,ii)= work(10)*sqrt21
            ecptmp( 5,ii)= work(16)*sqrt63
            ecptmp( 6,ii)= work(11)*sqrt21
            ecptmp( 7,ii)= work(12)*sqrt21
            ecptmp( 8,ii)= work(19)*sqrt105
            ecptmp( 9,ii)= work(20)*sqrt105
            ecptmp(10,ii)= work(14)*sqrt21
            ecptmp(11,ii)= work( 6)*three
            ecptmp(12,ii)= work(17)*sqrt63
            ecptmp(13,ii)= work(21)*sqrt105
            ecptmp(14,ii)= work(18)*sqrt63
            ecptmp(15,ii)= work( 8)*three
            ecptmp(16,ii)= work( 2)
            ecptmp(17,ii)= work( 7)*three
            ecptmp(18,ii)= work(13)*sqrt21
            ecptmp(19,ii)= work(15)*sqrt21
            ecptmp(20,ii)= work( 9)*three
            ecptmp(21,ii)= work( 3)
          enddo
! I function
        case(13)
          do ii= 1,ncarti
            work(1:28)= ecptmp(1:28,ii)
            ecptmp( 1,ii)=(work(4)*six-work(19)*p20+work(6)*six)*faci1
            ecptmp( 2,ii)=(work(16)*five-work(24)*ten+work(7))*faci2
            ecptmp( 3,ii)=(-work(4)*four+work(23)*p40+work(6)*four-work(25)*p40)*faci3
            ecptmp( 4,ii)=(-work(16)*p9-work(24)*six+work(26)*p24+work(7)*three &
&                         -work(21)*eight)*faci4
            ecptmp( 5,ii)=(work(4)*two+work(19)*four-work(23)*p32+work(6)*two &
&                         -work(25)*p32+work(18)*p32)*faci5
            ecptmp( 6,ii)=(work(16)*five+work(24)*ten-work(26)*p20+work(7)*five &
&                         -work(21)*p20+work(9)*eight)*faci6
            ecptmp( 7,ii)=(-work(1)*five-work(10)*p15+work(11)*p90-work(12)*p15 &
&                         +work(28)*p180-work(14)*p120-work(2)*five+work(13)*p90 &
&                         -work(15)*p120+work(3)*p16)*sixteenth
            ecptmp( 8,ii)=(work(5)*five+work(22)*ten-work(20)*p20+work(17)*five &
&                         -work(27)*p20+work(8)*eight)*faci6
            ecptmp( 9,ii)=(work(1)+work(10)-work(11)*p16-work(12)+work(14)*p16 &
&                         -work(2)+work(13)*p16-work(15)*p16)*faci5
            ecptmp(10,ii)=(-work(5)*three+work(22)*six+work(20)*eight+work(17)*p9 &
&                         -work(27)*p24)*faci4
            ecptmp(11,ii)=(-work(1)+work(10)*five+work(11)*ten+work(12)*five &
&                         -work(28)*p60-work(2)+work(13)*ten)*faci3
            ecptmp(12,ii)=(work(5)-work(22)*ten+work(17)*five)*faci2
            ecptmp(13,ii)=(work(1)-work(10)*p15+work(12)*p15-work(2))*faci1
          enddo
        case(28)
          do ii= 1,ncarti
            work(1:28)= ecptmp(1:28,ii)
            ecptmp( 1,ii)= work( 1)
            ecptmp( 2,ii)= work( 4)*sqrt11
            ecptmp( 3,ii)= work( 5)*sqrt11
            ecptmp( 4,ii)= work(10)*sqrt33
            ecptmp( 5,ii)= work(16)*sqrt99
            ecptmp( 6,ii)= work(11)*sqrt33
            ecptmp( 7,ii)= work(19)*sqrt231fifth
            ecptmp( 8,ii)= work(22)*sqrt231
            ecptmp( 9,ii)= work(23)*sqrt231
            ecptmp(10,ii)= work(20)*sqrt231fifth
            ecptmp(11,ii)= work(12)*sqrt33
            ecptmp(12,ii)= work(24)*sqrt231
            ecptmp(13,ii)= work(28)*sqrt385
            ecptmp(14,ii)= work(26)*sqrt231
            ecptmp(15,ii)= work(14)*sqrt33
            ecptmp(16,ii)= work( 6)*sqrt11
            ecptmp(17,ii)= work(17)*sqrt99
            ecptmp(18,ii)= work(25)*sqrt231
            ecptmp(19,ii)= work(27)*sqrt231
            ecptmp(20,ii)= work(18)*sqrt99
            ecptmp(21,ii)= work( 8)*sqrt11
            ecptmp(22,ii)= work( 2)
            ecptmp(23,ii)= work( 7)*sqrt11
            ecptmp(24,ii)= work(13)*sqrt33
            ecptmp(25,ii)= work(21)*sqrt231fifth
            ecptmp(26,ii)= work(15)*sqrt33
            ecptmp(27,ii)= work( 9)*sqrt11
            ecptmp(28,ii)= work( 3)
          enddo
      end select
!
! Ket part
!
      select case(nbfij(1))
! D function
        case(5)
          do jj= 1,nbfij(2)
            do ii= 1,6
              work(ii)= ecptmp(jj,ii)
            enddo
            ecptmp(jj,1)= work(4)*sqrt3
            ecptmp(jj,2)= work(6)*sqrt3
            ecptmp(jj,3)=(work(3)*two-work(1)-work(2))*half
            ecptmp(jj,4)= work(5)*sqrt3
            ecptmp(jj,5)=(work(1)-work(2))*sqrt3h
          enddo
        case(6)
          do jj= 1,nbfij(2)
            work(1:6)= ecptmp(jj,1:6)
            ecptmp(jj,1)= work(1)
            ecptmp(jj,2)= work(4)*sqrt3
            ecptmp(jj,3)= work(5)*sqrt3
            ecptmp(jj,4)= work(2)
            ecptmp(jj,5)= work(6)*sqrt3
            ecptmp(jj,6)= work(3)
          enddo
! F function
        case(7)
          do jj= 1,nbfij(2)
            do ii= 1,10
              work(ii)= ecptmp(jj,ii)
            enddo
            ecptmp(jj,1)=(-work(2)+three*work(4)                  )*facf1
            ecptmp(jj,2)=  work(10)                                *facf2
            ecptmp(jj,3)=(-work(2)-work(4)+four*work(9)           )*facf3
            ecptmp(jj,4)=( two*work(3)-three*work(5)-three*work(7))*half
            ecptmp(jj,5)=(-work(1)-work(6)+four*work(8)           )*facf3
            ecptmp(jj,6)=( work(5)-work(7)                        )*facf4
            ecptmp(jj,7)=( work(1)-three*work(6)                  )*facf1
          enddo
        case(10)
          do jj= 1,nbfij(2)
            work(1:10)= ecptmp(jj,1:10)
            ecptmp(jj, 1)= work( 1)
            ecptmp(jj, 2)= work( 4)*sqrt5
            ecptmp(jj, 3)= work( 5)*sqrt5
            ecptmp(jj, 4)= work( 6)*sqrt5
            ecptmp(jj, 5)= work(10)*sqrt15
            ecptmp(jj, 6)= work( 8)*sqrt5
            ecptmp(jj, 7)= work( 2)
            ecptmp(jj, 8)= work( 7)*sqrt5
            ecptmp(jj, 9)= work( 9)*sqrt5
            ecptmp(jj,10)= work( 3)
          enddo
! G function
        case(9)
          do jj= 1,nbfij(2)
            do ii= 1,15
              work(ii)= ecptmp(jj,ii)
            enddo
            ecptmp(jj,1)=(work(4)-work(6))*facg1
            ecptmp(jj,2)=(-work(7)+work(13)*three)*facg2
            ecptmp(jj,3)=(-work(4)-work(6)+work(15)*six)*facg3
            ecptmp(jj,4)=(-work(7)*three+work(9)*four-work(13)*three)*facg4
            ecptmp(jj,5)=(work(1)*three+work(2)*three+work(3)*eight+work(10)*six &
&                       -work(11)*p24-work(12)*p24)*eighth
            ecptmp(jj,6)=(-work(5)*three+work(8)*four-work(14)*three)*facg4
            ecptmp(jj,7)=(-work(1)+work(2)+work(11)*six-work(12)*six)*facg5
            ecptmp(jj,8)=(work(5)-work(14)*three)*facg2
            ecptmp(jj,9)=(work(1)+work(2)-work(10)*six)*facg6
          enddo
        case(15)
          do jj= 1,nbfij(2)
            work(1:15)= ecptmp(jj,1:15)
            ecptmp(jj, 1)= work( 1)
            ecptmp(jj, 2)= work( 4)*sqrt7
            ecptmp(jj, 3)= work( 5)*sqrt7
            ecptmp(jj, 4)= work(10)*sqrt35third
            ecptmp(jj, 5)= work(13)*sqrt35
            ecptmp(jj, 6)= work(11)*sqrt35third
            ecptmp(jj, 7)= work( 6)*sqrt7
            ecptmp(jj, 8)= work(14)*sqrt35
            ecptmp(jj, 9)= work(15)*sqrt35
            ecptmp(jj,10)= work( 8)*sqrt7
            ecptmp(jj,11)= work( 2)
            ecptmp(jj,12)= work( 7)*sqrt7
            ecptmp(jj,13)= work(12)*sqrt35third
            ecptmp(jj,14)= work( 9)*sqrt7
            ecptmp(jj,15)= work( 3)
          enddo
! H function
        case(11)
          do jj= 1,nbfij(2)
            work(1:21)= ecptmp(jj,1:21)
            ecptmp(jj, 1)=(work(4)*five-work(12)*ten+work(2))*fach1
            ecptmp(jj, 2)=(work(16)*four-work(17)*four)*fach2
            ecptmp(jj, 3)=(-work(4)*three-work(12)*two+work(20)*p24+work(2)-work(13)*eight)*fach3
            ecptmp(jj, 4)=(-work(16)*two-work(17)*two+work(18)*four)*fach4
            ecptmp(jj, 5)=(work(4)+work(12)*two-work(20)*twelve+work(2)-work(13)*twelve &
&                         +work(9)*eight)*fach5
            ecptmp(jj, 6)=(work(5)*p15+work(19)*p30-work(14)*p40+work(7)*p15-work(15)*p40 &
&                         +work(3)*eight)*eighth
            ecptmp(jj, 7)=(work(1)+work(10)*two-work(11)*twelve+work(6)-work(21)*twelve &
&                         +work(8)*eight)*fach5
            ecptmp(jj, 8)=(-work(5)+work(14)*two+work(7)-work(15)*two)*fach4
            ecptmp(jj, 9)=(-work(1)+work(10)*two+work(11)*eight+work(6)*three-work(21)*p24)*fach3
            ecptmp(jj,10)=(work(5)-work(19)*six+work(7))*fach2
            ecptmp(jj,11)=(work(1)-work(10)*ten+work(6)*five)*fach1
          enddo
        case(21)
          do jj= 1,nbfij(2)
            work(1:21)= ecptmp(jj,1:21)
            ecptmp(jj, 1)= work( 1)
            ecptmp(jj, 2)= work( 4)*three
            ecptmp(jj, 3)= work( 5)*three
            ecptmp(jj, 4)= work(10)*sqrt21
            ecptmp(jj, 5)= work(16)*sqrt63
            ecptmp(jj, 6)= work(11)*sqrt21
            ecptmp(jj, 7)= work(12)*sqrt21
            ecptmp(jj, 8)= work(19)*sqrt105
            ecptmp(jj, 9)= work(20)*sqrt105
            ecptmp(jj,10)= work(14)*sqrt21
            ecptmp(jj,11)= work( 6)*three
            ecptmp(jj,12)= work(17)*sqrt63
            ecptmp(jj,13)= work(21)*sqrt105
            ecptmp(jj,14)= work(18)*sqrt63
            ecptmp(jj,15)= work( 8)*three
            ecptmp(jj,16)= work( 2)
            ecptmp(jj,17)= work( 7)*three
            ecptmp(jj,18)= work(13)*sqrt21
            ecptmp(jj,19)= work(15)*sqrt21
            ecptmp(jj,20)= work( 9)*three
            ecptmp(jj,21)= work( 3)
          enddo
! I function
        case(13)
          do jj= 1,nbfij(2)
            work(1:28)= ecptmp(jj,1:28)
            ecptmp(jj, 1)=(work(4)*six-work(19)*p20+work(6)*six)*faci1
            ecptmp(jj, 2)=(work(16)*five-work(24)*ten+work(7))*faci2
            ecptmp(jj, 3)=(-work(4)*four+work(23)*p40+work(6)*four-work(25)*p40)*faci3
            ecptmp(jj, 4)=(-work(16)*p9-work(24)*six+work(26)*p24+work(7)*three &
&                         -work(21)*eight)*faci4
            ecptmp(jj, 5)=(work(4)*two+work(19)*four-work(23)*p32+work(6)*two &
&                         -work(25)*p32+work(18)*p32)*faci5
            ecptmp(jj, 6)=(work(16)*five+work(24)*ten-work(26)*p20+work(7)*five &
&                         -work(21)*p20+work(9)*eight)*faci6
            ecptmp(jj, 7)=(-work(1)*five-work(10)*p15+work(11)*p90-work(12)*p15 &
&                         +work(28)*p180-work(14)*p120-work(2)*five+work(13)*p90 &
&                         -work(15)*p120+work(3)*p16)*sixteenth
            ecptmp(jj, 8)=(work(5)*five+work(22)*ten-work(20)*p20+work(17)*five &
&                         -work(27)*p20+work(8)*eight)*faci6
            ecptmp(jj, 9)=(work(1)+work(10)-work(11)*p16-work(12)+work(14)*p16 &
&                         -work(2)+work(13)*p16-work(15)*p16)*faci5
            ecptmp(jj,10)=(-work(5)*three+work(22)*six+work(20)*eight+work(17)*p9 &
&                         -work(27)*p24)*faci4
            ecptmp(jj,11)=(-work(1)+work(10)*five+work(11)*ten+work(12)*five &
&                         -work(28)*p60-work(2)+work(13)*ten)*faci3
            ecptmp(jj,12)=(work(5)-work(22)*ten+work(17)*five)*faci2
            ecptmp(jj,13)=(work(1)-work(10)*p15+work(12)*p15-work(2))*faci1
          enddo
        case(28)
          do jj= 1,nbfij(2)
            work(1:28)= ecptmp(jj,1:28)
            ecptmp(jj, 1)= work( 1)
            ecptmp(jj, 2)= work( 4)*sqrt11
            ecptmp(jj, 3)= work( 5)*sqrt11
            ecptmp(jj, 4)= work(10)*sqrt33
            ecptmp(jj, 5)= work(16)*sqrt99
            ecptmp(jj, 6)= work(11)*sqrt33
            ecptmp(jj, 7)= work(19)*sqrt231fifth
            ecptmp(jj, 8)= work(22)*sqrt231
            ecptmp(jj, 9)= work(23)*sqrt231
            ecptmp(jj,10)= work(20)*sqrt231fifth
            ecptmp(jj,11)= work(12)*sqrt33
            ecptmp(jj,12)= work(24)*sqrt231
            ecptmp(jj,13)= work(28)*sqrt385
            ecptmp(jj,14)= work(26)*sqrt231
            ecptmp(jj,15)= work(14)*sqrt33
            ecptmp(jj,16)= work( 6)*sqrt11
            ecptmp(jj,17)= work(17)*sqrt99
            ecptmp(jj,18)= work(25)*sqrt231
            ecptmp(jj,19)= work(27)*sqrt231
            ecptmp(jj,20)= work(18)*sqrt99
            ecptmp(jj,21)= work( 8)*sqrt11
            ecptmp(jj,22)= work( 2)
            ecptmp(jj,23)= work( 7)*sqrt11
            ecptmp(jj,24)= work(13)*sqrt33
            ecptmp(jj,25)= work(21)*sqrt231fifth
            ecptmp(jj,26)= work(15)*sqrt33
            ecptmp(jj,27)= work( 9)*sqrt11
            ecptmp(jj,28)= work( 3)
          enddo
      end select
!
      do ii= 1,nbfij(1)
        do jj= 1,nbfij(2)
          ecpint(jj,ii)= ecptmp(jj,ii)
        enddo
      enddo
!
      return
end


!----------------------------------------------------------------------
  subroutine ecpcab(ecpint,exij,coij,coordijk,nprimij,nangij,nbfij, &
&                   exkecp,cokecp,nprimkecp,nangkecp,lmaxecp, &
&                   term1ecp,term2ecp,label1ecp,label2ecp, &
&                   num1ecp,num2ecp,numtbasis,len1,mxprsh)
!----------------------------------------------------------------------
!
! Calculate  <C|A|B>-type ECP integral
!
!
      use modecp, only : nx, ny, nz
      implicit none
      integer,intent(in) :: len1, mxprsh, nprimij(2), nangij(2), nbfij(2), nprimkecp(6)
      integer,intent(in) :: nangkecp(mxprsh,6), lmaxecp, numtbasis
      integer,intent(in) :: label1ecp(9,*), label2ecp(6,*), num1ecp(*), num2ecp(*)
      integer :: ncart(0:6)=(/1,3,6,10,15,21,28/), locxyzecp(0:7)=(/0,1,4,10,20,35,56,84/)
      integer :: nangbasis, ncarti, ncartj, nangbasis3, iprim, jprim, kk, nangk, nlm
      integer :: nn, iangij, ll, ll2, mu, nk, iindex, jindex, i1, j1, ii, jj, indx
      integer :: label1f, label1l, nxc, nyc, nzc, nxb, nyb, nzb, ka, la, kx, ky, kz
      integer :: kxp, kyp, kzp, nangll, nangmaxij, nangsum, nangca, nangba, nmax
      integer :: ltmax, lemax, mm, m1, m2, m3, m3f, m4, nangmaxll1, label2f, label2l
      integer :: icount, lmindx
      real(8),parameter :: zero=0.0D+00, quarter=0.25D+00, one=1.0D+00
      real(8),parameter :: pi4=1.256637061435917D+01
      real(8),parameter :: sqrtinvpi4=0.2820947917738781D+00
      real(8),parameter :: abthresh=1.0D-01, dathresh=1.0D-07
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00, four=4.0D+00
      real(8),parameter :: five=5.0D+00, six=6.0D+00, eight=8.0D+00, p9=9.0D+00, ten=10.0D+00
      real(8),parameter :: twelve=12.0D+00, p15=15.0D+00, p16=16.0D+00
      real(8),parameter :: p20=20.0D+00, p24=24.0D+00, p30=30.0D+00, p32=32.0D+00
      real(8),parameter :: p40=40.0D+00, p60=60.0D+00, p90=90.0D+00, p120=1.2D+02
      real(8),parameter :: p180=1.8D+02, eighth=0.125D+00, sixteenth=6.25D-02
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
      real(8),intent(in) :: exij(mxprsh,2), coij(mxprsh,2), exkecp(mxprsh,6)
      real(8),intent(in) :: cokecp(mxprsh,6), term1ecp(*), term2ecp(*), coordijk(3,3)
      real(8),intent(out) :: ecpint(len1,len1)
      real(8) :: ecptmp(28,28), work(28), r12, unitvec(3,3), bax, bay, baz, ba, cax, cay
      real(8) :: caz, ca, abx(0:6), aby(0:6), abz(0:6), acx(0:6), acy(0:6), acz(0:6)
      real(8) :: rad12int(2*11*11*11), ex01, ex02, ex12, co12, ex21, y01, y02, y12
      real(8) :: dax, day, daz, da, dazeta, da2zeta, besselarginv, zzlm(121), xp0
      real(8) :: sqrtinvzeta, alpha, xp1, co123, rad1int(78), radtmp, zeta, type1da0
      real(8) :: type12ang(23,12,12), type1sum, cazeta, ca2zeta, alfbessel, bazeta, ba2zeta
      real(8) :: betbessel, albe, beta, xka, xkb, a1, a2, xp1p, xp1m, rad2int(11,11,11)
      real(8) :: zzlmc(121), zzlmb(121), type2angc(11,11), type2angb(11,11), type2sum, angcc
!
      nangbasis= nangij(1)+nangij(2)
      ncarti= ncart(nangij(1))
      ncartj= ncart(nangij(2))
      ecptmp(1:ncartj,1:ncarti)= zero
!
      r12=(coordijk(1,1)-coordijk(1,2))**2+(coordijk(2,1)-coordijk(2,2))**2 &
&        +(coordijk(3,1)-coordijk(3,2))**2
      cax= coordijk(1,1)-coordijk(1,3)
      cay= coordijk(2,1)-coordijk(2,3)
      caz= coordijk(3,1)-coordijk(3,3)
      ca = sqrt(cax*cax+cay*cay+caz*caz)
      unitvec(1,2)= cax/ca
      unitvec(2,2)= cay/ca
      unitvec(3,2)= caz/ca
      bax= coordijk(1,2)-coordijk(1,3)
      bay= coordijk(2,2)-coordijk(2,3)
      baz= coordijk(3,2)-coordijk(3,3)
      ba = sqrt(bax*bax+bay*bay+baz*baz)
      unitvec(1,3)= bax/ba
      unitvec(2,3)= bay/ba
      unitvec(3,3)= baz/ba
!
      acx(0)= one
      acy(0)= one
      acz(0)= one
      do ii= 1,nangij(1)
        acx(ii)=-cax*acx(ii-1)
        acy(ii)=-cay*acy(ii-1)
        acz(ii)=-caz*acz(ii-1)
      enddo
      abx(0)= one
      aby(0)= one
      abz(0)= one
      do jj= 1,nangij(2)
        abx(jj)=-bax*abx(jj-1)
        aby(jj)=-bay*aby(jj-1)
        abz(jj)=-baz*abz(jj-1)
      enddo
!
! Type 1
!
      nangbasis3=((nangbasis+1)*(nangbasis+2)*(2*nangbasis+3))/6
      rad12int(1:nangbasis3)= zero
!
! Radial integral of Type 1
!
      do iprim= 1,nprimij(1)
        ex01= exij(iprim,1)
        do jprim= 1,nprimij(2)
          ex02= exij(jprim,2)
          ex12= ex01+ex02
          co12= coij(iprim,1)*coij(jprim,2)
          ex21= one/ex12
          y01= ex01*ex21
          y02= one-y01
          y12= y01*ex02
          dax= coordijk(1,1)+(coordijk(1,2)-coordijk(1,1))*y02-coordijk(1,3)
          day= coordijk(2,1)+(coordijk(2,2)-coordijk(2,1))*y02-coordijk(2,3)
          daz= coordijk(3,1)+(coordijk(3,2)-coordijk(3,1))*y02-coordijk(3,3)
          da = sqrt(dax*dax+day*day+daz*daz)
          co12= co12*exp(-r12*y12)
          if(da >= dathresh) then
            dazeta= ex12*da
            da2zeta= dazeta*da
            besselarginv= one/(dazeta+dazeta)
            unitvec(1,1)= dax/da
            unitvec(2,1)= day/da
            unitvec(3,1)= daz/da
            call calczspher(zzlm,unitvec,nangbasis)
            xp0 = exp(-da2zeta)
            do kk= 1,nprimkecp(1)
              sqrtinvzeta= one/sqrt(ex12+exkecp(kk,1))
              alpha= dazeta*sqrtinvzeta
              xp1= exp(-da2zeta+alpha*alpha)
              nangk= nangkecp(kk,1)
              co123= cokecp(kk,1)*co12
              call calctype1rad(rad1int,besselarginv,alpha,sqrtinvzeta, &
&                               xp0,xp1,nangk,nangbasis)
              nlm= 1
              nn= 1
              do iangij= 1,nangbasis+1
                do ll= 1,iangij
                  radtmp= rad1int(nn)*co123
                  ll2= ll*(ll-1)+1+ll
                  do mu= 1,2*ll-1
                    rad12int(nlm)= rad12int(nlm)+radtmp*zzlm(ll2-mu)
                    nlm= nlm+1
                  enddo
                  nn= nn+1
                enddo
              enddo
            enddo
          else
            do kk= 1,nprimkecp(1)
              zeta= ex12+exkecp(kk,1)
              nangk= nangkecp(kk,1)
              co123= cokecp(kk,1)*sqrtinvpi4*co12
              do nn= 1,nangbasis+1
                nk=(nn*(nn-1)*(2*nn-1))/6
                rad12int(nk)= rad12int(nk)+type1da0(nn+nangk,zeta)*co123
              enddo
            enddo
          endif
        enddo
      enddo
!
! Angular integral of Type 1
!
      iindex= locxyzecp(nangij(1))
      jindex= locxyzecp(nangij(2))
      do i1= 1,ncarti
        ii= iindex+i1
        do j1= 1,ncartj
          jj= jindex+j1
          if(ii >= jj) then
            indx= jj+(ii*(ii-1))/2
          else
            indx= ii+(jj*(jj-1))/2
          endif
          label1f= num1ecp(indx)
          label1l= num1ecp(indx+1)-1
          nxc= nx(ii)
          nyc= ny(ii)
          nzc= nz(ii)
          nxb= nx(jj)
          nyb= ny(jj)
          nzb= nz(jj)
          do ka= 1,nangbasis+2
            do la= 1,ka
              do mu= 1,2*la-1
                type12ang(mu,la,ka)= zero
              enddo
            enddo
          enddo
!
          do kk= label1f,label1l
            ka= label1ecp(1,kk)+1
            la= label1ecp(2,kk)+1
            mu= label1ecp(3,kk)
            if(ii >= jj) then
              kx = label1ecp(4,kk)
              ky = label1ecp(5,kk)
              kz = label1ecp(6,kk)
              kxp= label1ecp(7,kk)
              kyp= label1ecp(8,kk)
              kzp= label1ecp(9,kk)
            else
              kxp= label1ecp(4,kk)
              kyp= label1ecp(5,kk)
              kzp= label1ecp(6,kk)
              kx = label1ecp(7,kk)
              ky = label1ecp(8,kk)
              kz = label1ecp(9,kk)
            endif
            type12ang(la+mu,la,ka)= type12ang(la+mu,la,ka)+term1ecp(kk) &
&                                  *acx(nxc-kx) *acy(nyc-ky) *acz(nzc-kz) &
&                                  *abx(nxb-kxp)*aby(nyb-kyp)*abz(nzb-kzp)
          enddo
!
          type1sum= zero
          nn= 1
          do ka=1,nangbasis+1
            do la=1,ka
              do mu=1,2*la-1
                type1sum= type1sum+rad12int(nn)*type12ang(mu,la,ka)
                nn= nn+1
              enddo
            enddo
          enddo
          ecptmp(j1,i1)= ecptmp(j1,i1)+pi4*type1sum
        enddo
      enddo
!
! Type 2
!
      do ll= 2,lmaxecp+1
        nangll= ll-2
        nangmaxij  = max(nangij(1),nangij(2))+1
        nangsum= nangmaxij+nangll
        nangca= nangij(1)+nangll+1
        nangba= nangij(2)+nangll+1
        nmax = nangbasis*nangca*(nangba/2+1)+6
        ltmax= max(nangbasis+1,nangsum)
        lemax = max(nangll,ltmax/2)
!
! Radial integral of Type 2
!
        rad12int(1:nmax)= zero
        do iprim= 1,nprimij(1)
          ex01= exij(iprim,1)
          cazeta= ex01*ca
          ca2zeta= cazeta*ca
          alfbessel= one/(cazeta+cazeta)
          do jprim= 1,nprimij(2)
            ex02= exij(jprim,2)
            ex12= ex01+ex02
            bazeta= ex02*ba
            ba2zeta= bazeta*ba
            betbessel= one/(bazeta+bazeta)
            albe= ca2zeta+ba2zeta
            xp0 = exp(-albe)
            co12= coij(iprim,1)*coij(jprim,2)
            do kk= 1,nprimkecp(ll)
              sqrtinvzeta= one/sqrt(ex12+exkecp(kk,ll))
              alpha= cazeta*sqrtinvzeta
              beta= bazeta*sqrtinvzeta
              if(alpha*beta > abthresh) then
                xka= cazeta*two
                xkb= bazeta*two
                a1 = alpha+beta
                a2 = alpha-beta
                zeta= ex12+exkecp(kk,ll)
                xp1p= exp(-albe+(alpha+beta)**2)
                xp1m= exp(-albe+(alpha-beta)**2)
              else
                xp1p= exp(-albe+alpha*alpha)
                xp1m= exp(-albe+beta*beta)
              endif
              nangk= nangkecp(kk,ll)
              co123= cokecp(kk,ll)*co12
              call calctype2rad(rad2int,alfbessel,betbessel,alpha,beta,sqrtinvzeta, &
&                               xp1p,xp1m,xp0,xka,xkb,a1,a2,zeta,nangk,nangbasis+1,ltmax,lemax)
              do mm= 1,lemax+1
                rad12int(mm)= rad12int(mm)+rad2int(mm,mm,1)*co123
              enddo
              nn= 7
              do m1= 2,nangbasis+1
                m3f= mod(m1,2)
                do m2= 1,nangca
                  m3f= 1-m3f
                  do m3= m3f+1,nangba,2
                    rad12int(nn)= rad12int(nn)+rad2int(m3,m2,m1)*co123
                    nn= nn+1
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
!
! Angular integral of Type 2
!
        nangmaxll1= max(1,nangll+1)
        ll2=(ll-1)*(ll-2)
        call calczspher(zzlmc,unitvec(1,2),nangca-1)
        call calczspher(zzlmb,unitvec(1,3),nangba-1)
        do i1= 1,ncarti
          ii= iindex+i1
          do j1= 1,ncartj
            jj= jindex+j1
            nxc= nx(ii)
            nyc= ny(ii)
            nzc= nz(ii)
            nxb= nx(jj)
            nyb= ny(jj)
            nzb= nz(jj)
            type12ang(1:nangba,1:nangca,1:(2*nangmaxij-1))= zero
            nlm=(ll2-nangll)*numtbasis
            do mm= nangll,-nangll,-1
              type2angc(1:nangsum,1:nangmaxij)= zero
              type2angb(1:nangsum,1:nangmaxij)= zero
              label2f= num2ecp(nlm+ii)
              label2l= num2ecp(nlm+ii+1)-1
              do icount= label2f,label2l
                la= label2ecp(1,icount)+1
                ka= label2ecp(2,icount)+1
                mu= label2ecp(3,icount)
                kx= label2ecp(4,icount)
                ky= label2ecp(5,icount)
                kz= label2ecp(6,icount)
                lmindx= la*(la-1)-mu+1
                type2angc(la,ka)= type2angc(la,ka)+term2ecp(icount)*zzlmc(lmindx) &
&                                                 *acx(nxc-kx)*acy(nyc-ky)*acz(nzc-kz)
              enddo
              label2f= num2ecp(nlm+jj)
              label2l= num2ecp(nlm+jj+1)-1
              do icount= label2f,label2l
                la= label2ecp(1,icount)+1
                ka= label2ecp(2,icount)+1
                mu= label2ecp(3,icount)
                kx= label2ecp(4,icount)
                ky= label2ecp(5,icount)
                kz= label2ecp(6,icount)
                lmindx= la*(la-1)-mu+1
                type2angb(la,ka)= type2angb(la,ka)+term2ecp(icount)*zzlmb(lmindx) &
&                                                 *abx(nxb-kx)*aby(nyb-ky)*abz(nzb-kz)
              enddo
              do m1=1,nangmaxij
                do m2=1,nangca
                  angcc= type2angc(m2,m1)
                  do m3=m1,m1+nangmaxij-1
                    do m4=1,nangba
                      type12ang(m4,m2,m3)= type12ang(m4,m2,m3)+angcc*type2angb(m4,m3+1-m1)
                    enddo
                  enddo
                enddo
              enddo
              nlm= nlm+numtbasis
            enddo
!
            type2sum= zero
            do mm= 1,nangmaxll1
              type2sum= type2sum+rad12int(mm)*type12ang(mm,mm,1)
            enddo
            nn= 7
            do m1=2,nangbasis+1
              m3f= mod(m1,2)
              do m2=1,nangca
                m3f= 1-m3f
                do m3= m3f+1,nangba,2
                  type2sum= type2sum+rad12int(nn)*type12ang(m3,m2,m1)
                  nn= nn+1
                enddo
              enddo
            enddo
            ecptmp(j1,i1)= ecptmp(j1,i1)+pi4*pi4*type2sum
          enddo
        enddo
      enddo
!
! Normalize ECP integrals
!
! Bra part
!
      select case(nbfij(2))
! D function
        case(5)
          do ii= 1,ncarti
            do jj= 1,6
              work(jj)= ecptmp(jj,ii)
            enddo
            ecptmp(1,ii)= work(4)*sqrt3
            ecptmp(2,ii)= work(6)*sqrt3
            ecptmp(3,ii)=(work(3)*two-work(1)-work(2))*half
            ecptmp(4,ii)= work(5)*sqrt3
            ecptmp(5,ii)=(work(1)-work(2))*sqrt3h
          enddo
        case(6)
          do ii= 1,ncarti
            work(1:6)= ecptmp(1:6,ii)
            ecptmp(1,ii)= work(1)
            ecptmp(2,ii)= work(4)*sqrt3
            ecptmp(3,ii)= work(5)*sqrt3
            ecptmp(4,ii)= work(2)
            ecptmp(5,ii)= work(6)*sqrt3
            ecptmp(6,ii)= work(3)
          enddo
! F function
        case(7)
          do ii= 1,ncarti
            do jj= 1,10
              work(jj)= ecptmp(jj,ii)
            enddo
            ecptmp(1,ii)=(-work(2)+three*work(4)                  )*facf1
            ecptmp(2,ii)=  work(10)                                *facf2
            ecptmp(3,ii)=(-work(2)-work(4)+four*work(9)           )*facf3
            ecptmp(4,ii)=( two*work(3)-three*work(5)-three*work(7))*half
            ecptmp(5,ii)=(-work(1)-work(6)+four*work(8)           )*facf3
            ecptmp(6,ii)=( work(5)-work(7)                        )*facf4
            ecptmp(7,ii)=( work(1)-three*work(6)                  )*facf1
          enddo
        case(10)
          do ii= 1,ncarti
            work(1:10)= ecptmp(1:10,ii)
            ecptmp( 1,ii)= work( 1)
            ecptmp( 2,ii)= work( 4)*sqrt5
            ecptmp( 3,ii)= work( 5)*sqrt5
            ecptmp( 4,ii)= work( 6)*sqrt5
            ecptmp( 5,ii)= work(10)*sqrt15
            ecptmp( 6,ii)= work( 8)*sqrt5
            ecptmp( 7,ii)= work( 2)
            ecptmp( 8,ii)= work( 7)*sqrt5
            ecptmp( 9,ii)= work( 9)*sqrt5
            ecptmp(10,ii)= work( 3)
          enddo
! G function
        case(9)
          do ii= 1,ncarti
            do jj= 1,15
              work(jj)= ecptmp(jj,ii)
            enddo
            ecptmp(1,ii)=(work(4)-work(6))*facg1
            ecptmp(2,ii)=(-work(7)+work(13)*three)*facg2
            ecptmp(3,ii)=(-work(4)-work(6)+work(15)*six)*facg3
            ecptmp(4,ii)=(-work(7)*three+work(9)*four-work(13)*three)*facg4
            ecptmp(5,ii)=(work(1)*three+work(2)*three+work(3)*eight+work(10)*six &
&                        -work(11)*p24-work(12)*p24)*eighth
            ecptmp(6,ii)=(-work(5)*three+work(8)*four-work(14)*three)*facg4
            ecptmp(7,ii)=(-work(1)+work(2)+work(11)*six-work(12)*six)*facg5
            ecptmp(8,ii)=(work(5)-work(14)*three)*facg2
            ecptmp(9,ii)=(work(1)+work(2)-work(10)*six)*facg6
          enddo
        case(15)
          do ii= 1,ncarti
            work(1:15)= ecptmp(1:15,ii)
            ecptmp( 1,ii)= work( 1)
            ecptmp( 2,ii)= work( 4)*sqrt7
            ecptmp( 3,ii)= work( 5)*sqrt7
            ecptmp( 4,ii)= work(10)*sqrt35third
            ecptmp( 5,ii)= work(13)*sqrt35
            ecptmp( 6,ii)= work(11)*sqrt35third
            ecptmp( 7,ii)= work( 6)*sqrt7
            ecptmp( 8,ii)= work(14)*sqrt35
            ecptmp( 9,ii)= work(15)*sqrt35
            ecptmp(10,ii)= work( 8)*sqrt7
            ecptmp(11,ii)= work( 2)
            ecptmp(12,ii)= work( 7)*sqrt7
            ecptmp(13,ii)= work(12)*sqrt35third
            ecptmp(14,ii)= work( 9)*sqrt7
            ecptmp(15,ii)= work( 3)
          enddo
! H function
        case(11)
          do ii= 1,ncarti
            work(1:21)= ecptmp(1:21,ii)
            ecptmp( 1,ii)=(work(4)*five-work(12)*ten+work(2))*fach1
            ecptmp( 2,ii)=(work(16)*four-work(17)*four)*fach2
            ecptmp( 3,ii)=(-work(4)*three-work(12)*two+work(20)*p24+work(2)-work(13)*eight)*fach3
            ecptmp( 4,ii)=(-work(16)*two-work(17)*two+work(18)*four)*fach4
            ecptmp( 5,ii)=(work(4)+work(12)*two-work(20)*twelve+work(2)-work(13)*twelve &
&                         +work(9)*eight)*fach5
            ecptmp( 6,ii)=(work(5)*p15+work(19)*p30-work(14)*p40+work(7)*p15-work(15)*p40 &
&                         +work(3)*eight)*eighth
            ecptmp( 7,ii)=(work(1)+work(10)*two-work(11)*twelve+work(6)-work(21)*twelve &
&                         +work(8)*eight)*fach5
            ecptmp( 8,ii)=(-work(5)+work(14)*two+work(7)-work(15)*two)*fach4
            ecptmp( 9,ii)=(-work(1)+work(10)*two+work(11)*eight+work(6)*three-work(21)*p24)*fach3
            ecptmp(10,ii)=(work(5)-work(19)*six+work(7))*fach2
            ecptmp(11,ii)=(work(1)-work(10)*ten+work(6)*five)*fach1
          enddo
        case(21)
          do ii= 1,ncarti
            work(1:21)= ecptmp(1:21,ii)
            ecptmp( 1,ii)= work( 1)
            ecptmp( 2,ii)= work( 4)*three
            ecptmp( 3,ii)= work( 5)*three
            ecptmp( 4,ii)= work(10)*sqrt21
            ecptmp( 5,ii)= work(16)*sqrt63
            ecptmp( 6,ii)= work(11)*sqrt21
            ecptmp( 7,ii)= work(12)*sqrt21
            ecptmp( 8,ii)= work(19)*sqrt105
            ecptmp( 9,ii)= work(20)*sqrt105
            ecptmp(10,ii)= work(14)*sqrt21
            ecptmp(11,ii)= work( 6)*three
            ecptmp(12,ii)= work(17)*sqrt63
            ecptmp(13,ii)= work(21)*sqrt105
            ecptmp(14,ii)= work(18)*sqrt63
            ecptmp(15,ii)= work( 8)*three
            ecptmp(16,ii)= work( 2)
            ecptmp(17,ii)= work( 7)*three
            ecptmp(18,ii)= work(13)*sqrt21
            ecptmp(19,ii)= work(15)*sqrt21
            ecptmp(20,ii)= work( 9)*three
            ecptmp(21,ii)= work( 3)
          enddo
! I function
        case(13)
          do ii= 1,ncarti
            work(1:28)= ecptmp(1:28,ii)
            ecptmp( 1,ii)=(work(4)*six-work(19)*p20+work(6)*six)*faci1
            ecptmp( 2,ii)=(work(16)*five-work(24)*ten+work(7))*faci2
            ecptmp( 3,ii)=(-work(4)*four+work(23)*p40+work(6)*four-work(25)*p40)*faci3
            ecptmp( 4,ii)=(-work(16)*p9-work(24)*six+work(26)*p24+work(7)*three &
&                         -work(21)*eight)*faci4
            ecptmp( 5,ii)=(work(4)*two+work(19)*four-work(23)*p32+work(6)*two &
&                         -work(25)*p32+work(18)*p32)*faci5
            ecptmp( 6,ii)=(work(16)*five+work(24)*ten-work(26)*p20+work(7)*five &
&                         -work(21)*p20+work(9)*eight)*faci6
            ecptmp( 7,ii)=(-work(1)*five-work(10)*p15+work(11)*p90-work(12)*p15 &
&                         +work(28)*p180-work(14)*p120-work(2)*five+work(13)*p90 &
&                         -work(15)*p120+work(3)*p16)*sixteenth
            ecptmp( 8,ii)=(work(5)*five+work(22)*ten-work(20)*p20+work(17)*five &
&                         -work(27)*p20+work(8)*eight)*faci6
            ecptmp( 9,ii)=(work(1)+work(10)-work(11)*p16-work(12)+work(14)*p16 &
&                         -work(2)+work(13)*p16-work(15)*p16)*faci5
            ecptmp(10,ii)=(-work(5)*three+work(22)*six+work(20)*eight+work(17)*p9 &
&                         -work(27)*p24)*faci4
            ecptmp(11,ii)=(-work(1)+work(10)*five+work(11)*ten+work(12)*five &
&                         -work(28)*p60-work(2)+work(13)*ten)*faci3
            ecptmp(12,ii)=(work(5)-work(22)*ten+work(17)*five)*faci2
            ecptmp(13,ii)=(work(1)-work(10)*p15+work(12)*p15-work(2))*faci1
          enddo
        case(28)
          do ii= 1,ncarti
            work(1:28)= ecptmp(1:28,ii)
            ecptmp( 1,ii)= work( 1)
            ecptmp( 2,ii)= work( 4)*sqrt11
            ecptmp( 3,ii)= work( 5)*sqrt11
            ecptmp( 4,ii)= work(10)*sqrt33
            ecptmp( 5,ii)= work(16)*sqrt99
            ecptmp( 6,ii)= work(11)*sqrt33
            ecptmp( 7,ii)= work(19)*sqrt231fifth
            ecptmp( 8,ii)= work(22)*sqrt231
            ecptmp( 9,ii)= work(23)*sqrt231
            ecptmp(10,ii)= work(20)*sqrt231fifth
            ecptmp(11,ii)= work(12)*sqrt33
            ecptmp(12,ii)= work(24)*sqrt231
            ecptmp(13,ii)= work(28)*sqrt385
            ecptmp(14,ii)= work(26)*sqrt231
            ecptmp(15,ii)= work(14)*sqrt33
            ecptmp(16,ii)= work( 6)*sqrt11
            ecptmp(17,ii)= work(17)*sqrt99
            ecptmp(18,ii)= work(25)*sqrt231
            ecptmp(19,ii)= work(27)*sqrt231
            ecptmp(20,ii)= work(18)*sqrt99
            ecptmp(21,ii)= work( 8)*sqrt11
            ecptmp(22,ii)= work( 2)
            ecptmp(23,ii)= work( 7)*sqrt11
            ecptmp(24,ii)= work(13)*sqrt33
            ecptmp(25,ii)= work(21)*sqrt231fifth
            ecptmp(26,ii)= work(15)*sqrt33
            ecptmp(27,ii)= work( 9)*sqrt11
            ecptmp(28,ii)= work( 3)
          enddo
      end select
!
! Ket part
!
      select case(nbfij(1))
! D function
        case(5)
          do jj= 1,nbfij(2)
            do ii= 1,6
              work(ii)= ecptmp(jj,ii)
            enddo
            ecptmp(jj,1)= work(4)*sqrt3
            ecptmp(jj,2)= work(6)*sqrt3
            ecptmp(jj,3)=(work(3)*two-work(1)-work(2))*half
            ecptmp(jj,4)= work(5)*sqrt3
            ecptmp(jj,5)=(work(1)-work(2))*sqrt3h
          enddo
        case(6)
          do jj= 1,nbfij(2)
            work(1:6)= ecptmp(jj,1:6)
            ecptmp(jj,1)= work(1)
            ecptmp(jj,2)= work(4)*sqrt3
            ecptmp(jj,3)= work(5)*sqrt3
            ecptmp(jj,4)= work(2)
            ecptmp(jj,5)= work(6)*sqrt3
            ecptmp(jj,6)= work(3)
          enddo
! F function
        case(7)
          do jj= 1,nbfij(2)
            do ii= 1,10
              work(ii)= ecptmp(jj,ii)
            enddo
            ecptmp(jj,1)=(-work(2)+three*work(4)                  )*facf1
            ecptmp(jj,2)=  work(10)                                *facf2
            ecptmp(jj,3)=(-work(2)-work(4)+four*work(9)           )*facf3
            ecptmp(jj,4)=( two*work(3)-three*work(5)-three*work(7))*half
            ecptmp(jj,5)=(-work(1)-work(6)+four*work(8)           )*facf3
            ecptmp(jj,6)=( work(5)-work(7)                        )*facf4
            ecptmp(jj,7)=( work(1)-three*work(6)                  )*facf1
          enddo
        case(10)
          do jj= 1,nbfij(2)
            work(1:10)= ecptmp(jj,1:10)
            ecptmp(jj, 1)= work( 1)
            ecptmp(jj, 2)= work( 4)*sqrt5
            ecptmp(jj, 3)= work( 5)*sqrt5
            ecptmp(jj, 4)= work( 6)*sqrt5
            ecptmp(jj, 5)= work(10)*sqrt15
            ecptmp(jj, 6)= work( 8)*sqrt5
            ecptmp(jj, 7)= work( 2)
            ecptmp(jj, 8)= work( 7)*sqrt5
            ecptmp(jj, 9)= work( 9)*sqrt5
            ecptmp(jj,10)= work( 3)
          enddo
! G function
        case(9)
          do jj= 1,nbfij(2)
            do ii= 1,15
              work(ii)= ecptmp(jj,ii)
            enddo
            ecptmp(jj,1)=(work(4)-work(6))*facg1
            ecptmp(jj,2)=(-work(7)+work(13)*three)*facg2
            ecptmp(jj,3)=(-work(4)-work(6)+work(15)*six)*facg3
            ecptmp(jj,4)=(-work(7)*three+work(9)*four-work(13)*three)*facg4
            ecptmp(jj,5)=(work(1)*three+work(2)*three+work(3)*eight+work(10)*six &
&                       -work(11)*p24-work(12)*p24)*eighth
            ecptmp(jj,6)=(-work(5)*three+work(8)*four-work(14)*three)*facg4
            ecptmp(jj,7)=(-work(1)+work(2)+work(11)*six-work(12)*six)*facg5
            ecptmp(jj,8)=(work(5)-work(14)*three)*facg2
            ecptmp(jj,9)=(work(1)+work(2)-work(10)*six)*facg6
          enddo
        case(15)
          do jj= 1,nbfij(2)
            work(1:15)= ecptmp(jj,1:15)
            ecptmp(jj, 1)= work( 1)
            ecptmp(jj, 2)= work( 4)*sqrt7
            ecptmp(jj, 3)= work( 5)*sqrt7
            ecptmp(jj, 4)= work(10)*sqrt35third
            ecptmp(jj, 5)= work(13)*sqrt35
            ecptmp(jj, 6)= work(11)*sqrt35third
            ecptmp(jj, 7)= work( 6)*sqrt7
            ecptmp(jj, 8)= work(14)*sqrt35
            ecptmp(jj, 9)= work(15)*sqrt35
            ecptmp(jj,10)= work( 8)*sqrt7
            ecptmp(jj,11)= work( 2)
            ecptmp(jj,12)= work( 7)*sqrt7
            ecptmp(jj,13)= work(12)*sqrt35third
            ecptmp(jj,14)= work( 9)*sqrt7
            ecptmp(jj,15)= work( 3)
          enddo
! H function
        case(11)
          do jj= 1,nbfij(2)
            work(1:21)= ecptmp(jj,1:21)
            ecptmp(jj, 1)=(work(4)*five-work(12)*ten+work(2))*fach1
            ecptmp(jj, 2)=(work(16)*four-work(17)*four)*fach2
            ecptmp(jj, 3)=(-work(4)*three-work(12)*two+work(20)*p24+work(2)-work(13)*eight)*fach3
            ecptmp(jj, 4)=(-work(16)*two-work(17)*two+work(18)*four)*fach4
            ecptmp(jj, 5)=(work(4)+work(12)*two-work(20)*twelve+work(2)-work(13)*twelve &
&                         +work(9)*eight)*fach5
            ecptmp(jj, 6)=(work(5)*p15+work(19)*p30-work(14)*p40+work(7)*p15-work(15)*p40 &
&                         +work(3)*eight)*eighth
            ecptmp(jj, 7)=(work(1)+work(10)*two-work(11)*twelve+work(6)-work(21)*twelve &
&                         +work(8)*eight)*fach5
            ecptmp(jj, 8)=(-work(5)+work(14)*two+work(7)-work(15)*two)*fach4
            ecptmp(jj, 9)=(-work(1)+work(10)*two+work(11)*eight+work(6)*three-work(21)*p24)*fach3
            ecptmp(jj,10)=(work(5)-work(19)*six+work(7))*fach2
            ecptmp(jj,11)=(work(1)-work(10)*ten+work(6)*five)*fach1
          enddo
        case(21)
          do jj= 1,nbfij(2)
            work(1:21)= ecptmp(jj,1:21)
            ecptmp(jj, 1)= work( 1)
            ecptmp(jj, 2)= work( 4)*three
            ecptmp(jj, 3)= work( 5)*three
            ecptmp(jj, 4)= work(10)*sqrt21
            ecptmp(jj, 5)= work(16)*sqrt63
            ecptmp(jj, 6)= work(11)*sqrt21
            ecptmp(jj, 7)= work(12)*sqrt21
            ecptmp(jj, 8)= work(19)*sqrt105
            ecptmp(jj, 9)= work(20)*sqrt105
            ecptmp(jj,10)= work(14)*sqrt21
            ecptmp(jj,11)= work( 6)*three
            ecptmp(jj,12)= work(17)*sqrt63
            ecptmp(jj,13)= work(21)*sqrt105
            ecptmp(jj,14)= work(18)*sqrt63
            ecptmp(jj,15)= work( 8)*three
            ecptmp(jj,16)= work( 2)
            ecptmp(jj,17)= work( 7)*three
            ecptmp(jj,18)= work(13)*sqrt21
            ecptmp(jj,19)= work(15)*sqrt21
            ecptmp(jj,20)= work( 9)*three
            ecptmp(jj,21)= work( 3)
          enddo
! I function
        case(13)
          do jj= 1,nbfij(2)
            work(1:28)= ecptmp(jj,1:28)
            ecptmp(jj, 1)=(work(4)*six-work(19)*p20+work(6)*six)*faci1
            ecptmp(jj, 2)=(work(16)*five-work(24)*ten+work(7))*faci2
            ecptmp(jj, 3)=(-work(4)*four+work(23)*p40+work(6)*four-work(25)*p40)*faci3
            ecptmp(jj, 4)=(-work(16)*p9-work(24)*six+work(26)*p24+work(7)*three &
&                         -work(21)*eight)*faci4
            ecptmp(jj, 5)=(work(4)*two+work(19)*four-work(23)*p32+work(6)*two &
&                         -work(25)*p32+work(18)*p32)*faci5
            ecptmp(jj, 6)=(work(16)*five+work(24)*ten-work(26)*p20+work(7)*five &
&                         -work(21)*p20+work(9)*eight)*faci6
            ecptmp(jj, 7)=(-work(1)*five-work(10)*p15+work(11)*p90-work(12)*p15 &
&                         +work(28)*p180-work(14)*p120-work(2)*five+work(13)*p90 &
&                         -work(15)*p120+work(3)*p16)*sixteenth
            ecptmp(jj, 8)=(work(5)*five+work(22)*ten-work(20)*p20+work(17)*five &
&                         -work(27)*p20+work(8)*eight)*faci6
            ecptmp(jj, 9)=(work(1)+work(10)-work(11)*p16-work(12)+work(14)*p16 &
&                         -work(2)+work(13)*p16-work(15)*p16)*faci5
            ecptmp(jj,10)=(-work(5)*three+work(22)*six+work(20)*eight+work(17)*p9 &
&                         -work(27)*p24)*faci4
            ecptmp(jj,11)=(-work(1)+work(10)*five+work(11)*ten+work(12)*five &
&                         -work(28)*p60-work(2)+work(13)*ten)*faci3
            ecptmp(jj,12)=(work(5)-work(22)*ten+work(17)*five)*faci2
            ecptmp(jj,13)=(work(1)-work(10)*p15+work(12)*p15-work(2))*faci1
          enddo
        case(28)
          do jj= 1,nbfij(2)
            work(1:28)= ecptmp(jj,1:28)
            ecptmp(jj, 1)= work( 1)
            ecptmp(jj, 2)= work( 4)*sqrt11
            ecptmp(jj, 3)= work( 5)*sqrt11
            ecptmp(jj, 4)= work(10)*sqrt33
            ecptmp(jj, 5)= work(16)*sqrt99
            ecptmp(jj, 6)= work(11)*sqrt33
            ecptmp(jj, 7)= work(19)*sqrt231fifth
            ecptmp(jj, 8)= work(22)*sqrt231
            ecptmp(jj, 9)= work(23)*sqrt231
            ecptmp(jj,10)= work(20)*sqrt231fifth
            ecptmp(jj,11)= work(12)*sqrt33
            ecptmp(jj,12)= work(24)*sqrt231
            ecptmp(jj,13)= work(28)*sqrt385
            ecptmp(jj,14)= work(26)*sqrt231
            ecptmp(jj,15)= work(14)*sqrt33
            ecptmp(jj,16)= work( 6)*sqrt11
            ecptmp(jj,17)= work(17)*sqrt99
            ecptmp(jj,18)= work(25)*sqrt231
            ecptmp(jj,19)= work(27)*sqrt231
            ecptmp(jj,20)= work(18)*sqrt99
            ecptmp(jj,21)= work( 8)*sqrt11
            ecptmp(jj,22)= work( 2)
            ecptmp(jj,23)= work( 7)*sqrt11
            ecptmp(jj,24)= work(13)*sqrt33
            ecptmp(jj,25)= work(21)*sqrt231fifth
            ecptmp(jj,26)= work(15)*sqrt33
            ecptmp(jj,27)= work( 9)*sqrt11
            ecptmp(jj,28)= work( 3)
          enddo
      end select
!
      do ii= 1,nbfij(1)
        do jj= 1,nbfij(2)
          ecpint(jj,ii)= ecptmp(jj,ii)
        enddo
      enddo
!
      return
end


!--------------------
  subroutine ecpzlm
!--------------------
!
! Set up the coefficients for the real spherical harmonics
!
      use modecp, only : zlm
      implicit none
      integer :: i
      real(8),parameter :: one=1.0D+00, two=2.0D+00, three=3.0D+00, four=4.0D+00
      real(8),parameter :: five=5.0D+00, six=6.0D+00, seven=7.0D+00, eight=8.0D+00
      real(8),parameter :: nine=9.0D+00, ten=10.0D+00, eleven=11.0D+00, p13=13.0D+00
      real(8),parameter :: p14=14.0D+00, p15=15.0D+00, p17=17.0D+00, p19=19.0D+00
      real(8),parameter :: p20=20.0D+00, half=0.5D+00, pi4=1.256637061435917D+01
      real(8),parameter :: sr3=1.732050807568877D+00, sr5=2.236067977499790D+00
      real(8),parameter :: sr7=2.645751311064591D+00, sr11=3.316624790355400D+00
      real(8),parameter :: sr13=3.605551275463989D+00, sr15=3.872983346207417D+00
      real(8),parameter :: sr17=4.123105625617661D+00, sr19=4.358898943540674D+00
      real(8),parameter :: sr21=4.582575694955840D+00, sr35=5.916079783099616D+00
      real(8),parameter :: sr55=7.416198487095663D+00, sr77=8.774964387392122D+00
      real(8),parameter :: sr91=9.539392014169456D+00, sr119=1.090871211463571D+01
      real(8),parameter :: sr221=1.486606874731851D+01, sr247=1.571623364550171D+01
      real(8),parameter :: sr323=1.797220075561143D+01, sr385=1.962141687034858D+01
      real(8),parameter :: sr1001=3.163858403911275D+01, sr1309=3.618010503025109D+01
      real(8),parameter :: sr1365=3.694590640382233D+01, sr2431=4.930517214248420D+01
      real(8),parameter :: sr3003=5.479963503528103D+01, sr5005=7.074602462329597D+01
      real(8),parameter :: s13585=1.1655470818461174D+02
      real(8) :: fap(0:17), fat(0:17), tmp
!
      fap(0)= sqrt(one/pi4)
      fap(1)= fap(0)*sqrt(half)
      do i= 2,16,2
         fap(i  )= fap(i-2)*half
         fap(i+1)= fap(i-1)*half
      enddo
      do i= 0,17
        fat(i)= fap(i)*three
      enddo
!
! L= 0, ML= 0
!
      zlm(  1)= fap(0)
!
! L= 1, ML=+1,0,-1
!
      zlm(  2)= sr3*fap(0)
      zlm(  3)= zlm(  2)
      zlm(  4)= zlm(  2)
!
! L= 2, ML=+2...-2
!
      zlm(  5)= sr15*fap(2)
      zlm(  6)=-zlm(  5)
      zlm(  7)= zlm(  5)*two
      tmp     = sr5*fap(2)
      zlm(  8)= tmp*three
      zlm(  9)=-tmp
      zlm( 10)= zlm(  7)
      zlm( 11)= zlm(  7)
!
! L= 3, ML=+3...-3
!
      zlm( 12)= sr35*fap(3)
      zlm( 13)=-zlm( 12)*three
      zlm( 14)= sr21*sr5*fap(2)
      zlm( 15)=-zlm( 14)
      tmp     = sr21*fap(3)
      zlm( 16)= tmp*five
      zlm( 17)=-tmp
      tmp     = sr7*fap(2)
      zlm( 18)= tmp*five
      zlm( 19)=-tmp*three
      zlm( 20)= zlm( 16)
      zlm( 21)= zlm( 17)
      zlm( 22)= zlm( 14)*two
      zlm( 23)=-zlm( 13)
      zlm( 24)=-zlm( 12)
!
! L= 4, ML=+4...-4
!
      zlm( 25)= sr35*fat(6)
      zlm( 26)=-zlm( 25)*six
      zlm( 27)= zlm( 25)
      zlm( 28)= sr35*fat(3)
      zlm( 29)=-zlm( 28)*three
      tmp     = sr5*fat(4)
      zlm( 30)= tmp*seven
      zlm( 31)=-zlm( 30)
      zlm( 32)= tmp
      zlm( 33)=-tmp
      tmp     = sr5*fat(3)
      zlm( 34)= tmp*seven
      zlm( 35)=-tmp*three
      tmp     = fat(6)
      zlm( 36)= tmp*five*seven
      zlm( 37)=-tmp*five*six
      zlm( 38)= tmp*three
      zlm( 39)= zlm( 34)
      zlm( 40)= zlm( 35)
      tmp     = sr5*fat(2)
      zlm( 41)= tmp*seven
      zlm( 42)=-tmp
      zlm( 43)=-zlm( 29)
      zlm( 44)=-zlm( 28)
      zlm( 45)= sr35*fat(2)
      zlm( 46)=-zlm( 45)
!
! L= 5, ML=+5...-5
!
      zlm( 47)= sr77*fat(7)
      zlm( 48)=-zlm( 47)*ten
      zlm( 49)= zlm( 47)*five
      zlm( 50)= sr385*fat(6)
      zlm( 51)=-zlm( 50)*six
      zlm( 52)= zlm( 50)
      tmp     = sr385*fap(7)
      zlm( 53)= tmp*three*three
      zlm( 54)=-tmp*three*nine
      zlm( 55)= tmp*three
      zlm( 56)=-tmp
      tmp     = sr385*sr3*fap(4)
      zlm( 57)= tmp*three
      zlm( 58)=-zlm( 57)
      zlm( 59)=-tmp
      zlm( 60)= tmp
      tmp     = sr15*sr11*fap(6)
      zlm( 61)= tmp*2.1D+01
      zlm( 62)=-tmp*p14
      zlm( 63)= tmp
      tmp     = sr11*fap(6)
      zlm( 64)= tmp*6.3D+01
      zlm( 65)=-tmp*7.0D+01
      zlm( 66)= tmp*p15
      zlm( 67)= zlm( 61)
      zlm( 68)= zlm( 62)
      zlm( 69)= zlm( 63)
      tmp     = sr385*sr3*fap(2)
      zlm( 70)= tmp*three
      zlm( 71)=-tmp
      zlm( 72)=-zlm( 54)
      zlm( 73)=-zlm( 55)
      zlm( 74)=-zlm( 53)
      zlm( 75)=-zlm( 56)
      zlm( 76)= sr385*fat(2)
      zlm( 77)=-zlm( 76)
      zlm( 78)= zlm( 49)
      zlm( 79)= zlm( 48)
      zlm( 80)= zlm( 47)
!
! L= 6, ML=+6...-6
!
      zlm( 81)= sr3003*fap(9)
      zlm( 82)=-zlm( 81)*p15
      zlm( 83)=-zlm( 82)
      zlm( 84)=-zlm( 81)
      zlm( 85)= sr1001*fat(7)
      zlm( 86)=-zlm( 85)*ten
      zlm( 87)= zlm( 85)*five
      tmp     = sr91*fat(8)
      zlm( 88)= tmp*eleven
      zlm( 89)=-tmp*eleven*six
      zlm( 90)= zlm( 88)
      zlm( 91)=-tmp
      zlm( 92)= tmp*six
      zlm( 93)=-tmp
      tmp     = sr1365*fap(7)
      zlm( 94)= tmp*eleven
      zlm( 95)=-tmp*three*eleven
      zlm( 96)= tmp*three*three
      zlm( 97)=-tmp*three
      tmp     = sr1365*fap(9)
      zlm( 98)= tmp*three*eleven
      zlm( 99)=-zlm( 98)
      zlm(100)=-tmp*three*six
      zlm(101)= tmp*three*six
      zlm(102)= tmp
      zlm(103)=-tmp
      tmp     = sr21*sr13*fap(6)
      zlm(104)= tmp*three*eleven
      zlm(105)=-tmp*three*ten
      zlm(106)= tmp*five
      tmp     = sr13*fap(8)
      zlm(107)= tmp*2.31D+02
      zlm(108)=-tmp*3.15D+02
      zlm(109)= tmp*1.05D+02
      zlm(110)=-tmp*five
      zlm(111)= zlm(104)
      zlm(112)= zlm(105)
      zlm(113)= zlm(106)
      tmp     = sr1365*fap(7)
      zlm(114)= tmp*three*eleven
      zlm(115)=-tmp*three*six
      zlm(116)= tmp
      zlm(117)=-zlm( 95)
      zlm(118)=-zlm( 94)
      zlm(119)=-zlm( 96)
      zlm(120)=-zlm( 97)
      tmp     = sr91*fat(4)
      zlm(121)= tmp*eleven
      zlm(122)=-zlm(121)
      zlm(123)=-tmp
      zlm(124)= tmp
      zlm(125)= zlm( 87)
      zlm(126)= zlm( 86)
      zlm(127)= zlm( 85)
      tmp     = sr3003*fap(9)
      zlm(128)= tmp*six
      zlm(129)=-tmp*p20
      zlm(130)= zlm(128)
!
! L= 7, ML=-7
      zlm(131)= sr55*sr13*fat(10)
      zlm(132)= zlm(131)*seven*five
      zlm(133)=-zlm(131)*seven*three
      zlm(134)=-zlm(131)*seven
! L= 7, ML=-6
      zlm(135)= sr5005*fat(9)
      zlm(136)=-zlm(135)*p15
      zlm(137)=-zlm(135)
      zlm(138)=-zlm(136)
! L= 7, ML=-5
      tmp     = sr385*fat(10)
      zlm(139)= tmp*p13
      zlm(140)=-zlm(139)*ten
      zlm(141)= zlm(139)*five
      zlm(142)=-tmp
      zlm(143)=-zlm(142)*ten
      zlm(144)= zlm(142)*five
! L= 7, ML=-4
      tmp     = sr385*fat(8)
      zlm(145)= tmp*p13
      zlm(146)=-tmp*p13*six
      zlm(147)= zlm(145)
      zlm(148)=-tmp*three
      zlm(149)= tmp*three*six
      zlm(150)= zlm(148)
! L= 7, ML=-3
      tmp     = sr35*fat(10)
      zlm(151)= tmp*eleven*p13
      zlm(152)=-zlm(151)*three
      zlm(153)=-tmp*eleven*six
      zlm(154)=-zlm(153)*three
      zlm(155)= tmp*three
      zlm(156)=-tmp*nine
! L= 7, ML=-2
      tmp     = sr35*fat(9)
      zlm(157)= tmp*eleven*p13
      zlm(158)=-zlm(157)
      zlm(159)=-tmp*eleven*ten
      zlm(160)=-zlm(159)
      zlm(161)= tmp*p15
      zlm(162)=-zlm(161)
! L= 7, ML=-1
      tmp     = sr7*sr15*fap(10)
      zlm(163)= tmp*4.29D+02
      zlm(164)=-tmp*4.95D+02
      zlm(165)= tmp*1.35D+02
      zlm(166)=-tmp*five
! L= 7, ML= 0
      tmp     = sr15*fap(8)
      zlm(167)= tmp*4.29D+02
      zlm(168)=-tmp*6.93D+02
      zlm(169)= tmp*3.15D+02
      zlm(170)=-tmp*3.5D+01
! L= 7, ML= 1
      zlm(171)= zlm(163)
      zlm(172)= zlm(164)
      zlm(173)= zlm(165)
      zlm(174)= zlm(166)
! L= 7, ML= 2
      zlm(175)= zlm(157)*two
      zlm(176)= zlm(159)*two
      zlm(177)= zlm(161)*two
! L= 7, ML= 3
      zlm(178)=-zlm(152)
      zlm(179)=-zlm(151)
      zlm(180)=-zlm(154)
      zlm(181)=-zlm(153)
      zlm(182)=-zlm(156)
      zlm(183)=-zlm(155)
! L= 7, ML= 4
      zlm(184)= zlm(145)*four
      zlm(185)=-zlm(184)
      zlm(186)= zlm(148)*four
      zlm(187)=-zlm(186)
! L= 7, ML= 5
      zlm(188)= zlm(141)
      zlm(189)= zlm(140)
      zlm(190)= zlm(139)
      zlm(191)= zlm(144)
      zlm(192)= zlm(143)
      zlm(193)= zlm(142)
! L= 7, ML= 6
      tmp     = sr5005*fat(7)
      zlm(194)= tmp*three
      zlm(195)=-tmp*ten
      zlm(196)= zlm(194)
! L= 7, ML= 7
      zlm(197)=-zlm(131)
      zlm(198)=-zlm(132)
      zlm(199)=-zlm(133)
      zlm(200)=-zlm(134)
!
! L= 8, ML=-8
      zlm(201)= sr2431*sr5*fat(14)
      zlm(202)=-zlm(201)*2.8D+01
      zlm(203)= zlm(201)*7.0D+01
      zlm(204)= zlm(202)
      zlm(205)= zlm(201)
! L= 8, ML=-7
      zlm(206)= sr2431*sr5*fat(10)
      zlm(207)=-zlm(206)*seven*three
      zlm(208)= zlm(206)*seven*five
      zlm(209)=-zlm(206)*seven
! L= 8, ML=-6
      tmp     = sr2431*sr3*fap(11)
      zlm(210)= tmp*p15
      zlm(211)=-zlm(210)*p15
      zlm(212)=-zlm(211)
      zlm(213)=-zlm(210)
      zlm(214)=-tmp
      zlm(215)= zlm(210)
      zlm(216)=-zlm(215)
      zlm(217)= tmp
! L= 8, ML=-5
      tmp     = sr1001*sr17*fat(10)
      zlm(218)= tmp*five
      zlm(219)=-zlm(218)*ten
      zlm(220)= zlm(218)*five
      zlm(221)=-tmp
      zlm(222)= tmp*ten
      zlm(223)= zlm(221)*five
! L= 8, ML=-4
      tmp     = sr1309*fat(12)
      zlm(224)= tmp*p13*five
      zlm(225)= zlm(224)
      zlm(226)=-zlm(224)*six
      zlm(227)=-tmp*p13*two
      zlm(228)= zlm(227)
      zlm(229)=-zlm(227)*six
      zlm(230)= tmp
      zlm(231)= tmp
      zlm(232)=-tmp*six
! L= 8, ML=-3
      tmp     = sr1309*sr15*fap(10)
      zlm(233)= tmp*p13*three
      zlm(234)=-zlm(233)*three
      zlm(235)=-tmp*p13*two
      zlm(236)=-zlm(235)*three
      zlm(237)= tmp*three
      zlm(238)=-zlm(237)*three
! L= 8, ML=-2
      tmp     = sr119*sr5*fat(11)
      zlm(239)= tmp*eleven*p13
      zlm(240)=-zlm(239)
      zlm(241)= zlm(240)
      zlm(242)= zlm(239)
      zlm(243)= tmp*eleven*three
      zlm(244)=-zlm(243)
      zlm(245)=-tmp
      zlm(246)= tmp
! L= 8, ML=-1
      tmp     = sr17*fat(10)
      zlm(247)= tmp*7.15D+02
      zlm(248)=-tmp*1.001D+03
      zlm(249)= tmp*3.85D+02
      zlm(250)=-tmp*3.5D+01
! L= 8, ML= 0
      tmp     = sr17*fap(14)
      zlm(251)= tmp*6.435D+03
      zlm(252)=-tmp*1.2012D+04
      zlm(253)= tmp*6.930D+03
      zlm(254)=-tmp*1.260D+03
      zlm(255)= tmp*3.5D+01
! L= 8, ML= 1
      zlm(256)= zlm(247)
      zlm(257)= zlm(248)
      zlm(258)= zlm(249)
      zlm(259)= zlm(250)
! L= 8, ML= 2
      tmp     = sr119*sr5*fat(9)
      zlm(260)= tmp*eleven*p13
      zlm(261)=-zlm(260)
      zlm(262)= tmp*eleven*three
      zlm(263)=-tmp
! L= 8, ML= 3
      zlm(264)=-zlm(234)
      zlm(265)=-zlm(233)
      zlm(266)=-zlm(236)
      zlm(267)=-zlm(235)
      zlm(268)=-zlm(238)
      zlm(269)=-zlm(237)
! L= 8, ML= 4
      tmp     = sr1309*fat(8)
      zlm(270)= tmp*p13*five
      zlm(271)=-zlm(270)
      zlm(272)=-tmp*p13*two
      zlm(273)=-zlm(272)
      zlm(274)= tmp
      zlm(275)=-tmp
! L= 8, ML= 5
      tmp     = sr1001*sr17*fat(10)
      zlm(276)= tmp*five*five
      zlm(277)=-tmp*five*ten
      zlm(278)= tmp*five
      zlm(279)=-tmp*five
      zlm(280)= tmp*ten
      zlm(281)=-tmp
! L= 8, ML= 6
      tmp     = sr2431*sr3*fap(9)
      zlm(282)= tmp*p15*three
      zlm(283)=-tmp*p15*ten
      zlm(284)= zlm(282)
      zlm(285)=-tmp*three
      zlm(286)= tmp*ten
      zlm(287)= zlm(285)
! L= 8, ML= 7
      zlm(288)=-zlm(209)
      zlm(289)=-zlm(208)
      zlm(290)=-zlm(207)
      zlm(291)=-zlm(206)
! L= 8, ML= 8
      zlm(292)= zlm(201)*eight
      zlm(293)=-zlm(292)*seven
      zlm(294)=-zlm(293)
      zlm(295)=-zlm(292)
!
! L= 9, ML=-9
      zlm(296)= s13585*sr17*fap(15)
      zlm(297)=-zlm(296)*3.6D+01
      zlm(298)= zlm(296)*1.26D+02
      zlm(299)=-zlm(296)*8.4D+01
      zlm(300)= zlm(296)*nine
! L= 9, ML=-8
      zlm(301)= s13585*sr17*fat(14)
      zlm(302)=-zlm(301)*seven*four
      zlm(303)= zlm(301)*seven*ten
      zlm(304)= zlm(302)
      zlm(305)= zlm(301)
! L= 9, ML=-7
      tmp     = s13585*fat(15)
      zlm(306)= tmp*p17
      zlm(307)=-zlm(306)*seven*three
      zlm(308)= zlm(306)*seven*five
      zlm(309)=-zlm(306)*seven
      zlm(310)=-tmp
      zlm(311)= tmp*seven*three
      zlm(312)=-tmp*seven*five
      zlm(313)= tmp*seven
! L= 9, ML=-6
      tmp     = sr247*sr11*sr15*fap(11)
      zlm(314)= tmp*p17
      zlm(315)=-zlm(314)*p15
      zlm(316)=-zlm(315)
      zlm(317)=-zlm(314)
      zlm(318)=-tmp*three
      zlm(319)=-zlm(318)*p15
      zlm(320)=-zlm(319)
      zlm(321)=-zlm(318)
! L= 9, ML=-5
      tmp     = sr247*sr11*fat(13)
      zlm(322)= tmp*five*p17
      zlm(323)=-zlm(322)*ten
      zlm(324)= zlm(322)*five
      zlm(325)=-tmp*five*six
      zlm(326)=-zlm(325)*ten
      zlm(327)= zlm(325)*five
      zlm(328)= tmp
      zlm(329)=-tmp*ten
      zlm(330)= tmp*five
! L= 9, ML=-4
      tmp     = sr5005*sr19*fat(12)
      zlm(331)= tmp*p17
      zlm(332)= zlm(331)
      zlm(333)=-zlm(331)*six
      zlm(334)=-tmp*ten
      zlm(335)= zlm(334)
      zlm(336)=-zlm(334)*six
      zlm(337)= tmp
      zlm(338)= tmp
      zlm(339)=-tmp*six
! L= 9, ML=-3
      tmp     = sr15*sr77*sr19*fap(13)
      zlm(340)= tmp*2.21D+02
      zlm(341)=-zlm(340)*three
      zlm(342)=-tmp*p13*three*five
      zlm(343)=-zlm(342)*three
      zlm(344)= tmp*p13*three
      zlm(345)=-zlm(344)*three
      zlm(346)=-tmp
      zlm(347)= tmp*three
! L= 9, ML=-2
      tmp     = sr55*sr19*fat(11)
      zlm(348)= tmp*2.21D+02
      zlm(349)=-zlm(348)
      zlm(350)=-tmp*seven*p13*three
      zlm(351)=-zlm(350)
      zlm(352)= tmp*seven*p13
      zlm(353)=-zlm(352)
      zlm(354)=-tmp*seven
      zlm(355)=-zlm(354)
! L= 9, ML=-1
      tmp     = sr5*sr19*fat(14)
      zlm(356)= tmp*2.431D+03
      zlm(357)=-tmp*4.004D+03
      zlm(358)= tmp*2.002D+03
      zlm(359)=-tmp*3.08D+02
      zlm(360)= tmp*seven
! L= 9, ML= 0
      tmp     = sr19*fap(14)
      zlm(361)= tmp*1.2155D+04
      zlm(362)=-tmp*2.5740D+04
      zlm(363)= tmp*1.8018D+04
      zlm(364)=-tmp*4.620D+03
      zlm(365)= tmp*3.15D+02
! L= 9, ML= 1
      zlm(366)= zlm(356)
      zlm(367)= zlm(357)
      zlm(368)= zlm(358)
      zlm(369)= zlm(359)
      zlm(370)= zlm(360)
! L= 9, ML= 2
      tmp     = sr55*sr19*fat(9)
      zlm(371)= tmp*2.21D+02
      zlm(372)=-tmp*2.73D+02
      zlm(373)= tmp*9.1D+01
      zlm(374)=-tmp*seven
! L= 9, ML= 3
      zlm(375)=-zlm(341)
      zlm(376)=-zlm(340)
      zlm(377)=-zlm(343)
      zlm(378)=-zlm(342)
      zlm(379)=-zlm(345)
      zlm(380)=-zlm(344)
      zlm(381)=-zlm(347)
      zlm(382)=-zlm(346)
! L= 9, ML= 4
      tmp     = sr5005*sr19*fat(8)
      zlm(383)= tmp*p17
      zlm(384)=-zlm(383)
      zlm(385)=-tmp*ten
      zlm(386)=-zlm(385)
      zlm(387)= tmp
      zlm(388)=-zlm(387)
! L= 9, ML= 5
      zlm(389)= zlm(324)
      zlm(390)= zlm(323)
      zlm(391)= zlm(322)
      zlm(392)= zlm(327)
      zlm(393)= zlm(326)
      zlm(394)= zlm(325)
      zlm(395)= zlm(330)
      zlm(396)= zlm(329)
      zlm(397)= zlm(328)
! L= 9, ML= 6
      tmp     = sr247*sr15*sr11*fap(9)
      zlm(398)= tmp*p17*three
      zlm(399)=-tmp*p17*ten
      zlm(400)= zlm(398)
      zlm(401)=-tmp*three*three
      zlm(402)= tmp*three*ten
      zlm(403)= zlm(401)
! L= 9, ML= 7
      zlm(404)=-zlm(309)
      zlm(405)=-zlm(308)
      zlm(406)=-zlm(307)
      zlm(407)=-zlm(306)
      zlm(408)=-zlm(313)
      zlm(409)=-zlm(312)
      zlm(410)=-zlm(311)
      zlm(411)=-zlm(310)
! L= 9, ML= 8
      zlm(412)= zlm(301)*eight
      zlm(413)=-zlm(412)*seven
      zlm(414)=-zlm(413)
      zlm(415)=-zlm(412)
! L= 9, ML= 9
      zlm(416)= zlm(300)
      zlm(417)= zlm(299)
      zlm(418)= zlm(298)
      zlm(419)= zlm(297)
      zlm(420)= zlm(296)
!
! L=10, ML=-10
      zlm(421)= sr323*sr3003*fap(17)
      zlm(422)=-zlm(421)*4.5D+01
      zlm(423)= zlm(421)*2.1D+02
      zlm(424)=-zlm(423)
      zlm(425)=-zlm(422)
      zlm(426)=-zlm(421)
! L=10, ML=-9
      zlm(427)= sr323*sr3003*sr5*fap(15)
      zlm(428)=-zlm(427)*3.6D+01
      zlm(429)= zlm(427)*1.26D+02
      zlm(430)=-zlm(427)*8.4D+01
      zlm(431)= zlm(427)*nine
! L=10, ML=-8
      tmp     = sr5005*sr17*sr3*fap(16)
      zlm(432)= tmp*p19
      zlm(433)=-zlm(432)*seven*four
      zlm(434)= zlm(432)*seven*ten
      zlm(435)= zlm(433)
      zlm(436)= zlm(432)
      zlm(437)=-tmp
      zlm(438)=-zlm(437)*seven*four
      zlm(439)= zlm(437)*seven*ten
      zlm(440)= zlm(438)
      zlm(441)= zlm(437)
! L=10, ML=-7
      tmp     = sr5005*sr17*fat(15)
      zlm(442)= tmp*p19
      zlm(443)=-zlm(442)*seven*three
      zlm(444)= zlm(442)*seven*five
      zlm(445)=-zlm(442)*seven
      zlm(446)=-tmp*three
      zlm(447)=-zlm(446)*seven*three
      zlm(448)= zlm(446)*seven*five
      zlm(449)=-zlm(446)*seven
! L=10, ML=-6
      tmp     = sr5005*fat(17)
      zlm(450)= tmp*3.23D+02
      zlm(451)=-zlm(450)*p15
      zlm(452)=-zlm(451)
      zlm(453)=-zlm(450)
      zlm(454)=-tmp*1.02D+02
      zlm(455)=-zlm(454)*p15
      zlm(456)=-zlm(455)
      zlm(457)=-zlm(454)
      zlm(458)= tmp*three
      zlm(459)=-zlm(458)*p15
      zlm(460)=-zlm(459)
      zlm(461)=-zlm(458)
! L=10, ML=-5
      tmp     = sr1001*fat(13)
      zlm(462)= tmp*3.23D+02
      zlm(463)=-zlm(462)*ten
      zlm(464)= zlm(462)*five
      zlm(465)=-tmp*1.70D+02
      zlm(466)=-zlm(465)*ten
      zlm(467)= zlm(465)*five
      zlm(468)= tmp*p15
      zlm(469)=-zlm(468)*ten
      zlm(470)= zlm(468)*five
! L=10, ML=-4
      tmp     = sr5005*fat(14)
      zlm(471)= tmp*3.23D+02
      zlm(472)= zlm(471)
      zlm(473)=-zlm(471)*six
      zlm(474)=-tmp*2.55D+02
      zlm(475)= zlm(474)
      zlm(476)=-zlm(474)*six
      zlm(477)= tmp*4.5D+01
      zlm(478)= zlm(477)
      zlm(479)=-zlm(477)*six
      zlm(480)=-tmp
      zlm(481)= zlm(480)
      zlm(482)=-zlm(480)*six
! L=10, ML=-3
      tmp     = sr5005*fat(13)
      zlm(483)= tmp*3.23D+02
      zlm(484)=-zlm(483)*three
      zlm(485)=-tmp*3.57D+02
      zlm(486)=-zlm(485)*three
      zlm(487)= tmp*1.05D+02
      zlm(488)=-zlm(487)*three
      zlm(489)=-tmp*seven
      zlm(490)=-zlm(489)*three
! L=10, ML=-2
      tmp     = sr385*fat(16)
      zlm(491)= tmp*4.199D+03
      zlm(492)=-zlm(491)
      zlm(493)=-tmp*6.188D+03
      zlm(494)=-zlm(493)
      zlm(495)= tmp*2.730D+03
      zlm(496)=-zlm(495)
      zlm(497)=-tmp*3.64D+02
      zlm(498)=-zlm(497)
      zlm(499)= tmp*seven
      zlm(500)=-zlm(499)
! L=10, ML=-1
      tmp     = sr385*sr3*fap(14)
      zlm(501)= tmp*4.199D+03
      zlm(502)=-tmp*7.956D+03
      zlm(503)= tmp*4.914D+03
      zlm(504)=-tmp*1.092D+03
      zlm(505)= tmp*6.3D+01
! L=10, ML= 0
      tmp     = sr21*fap(16)
      zlm(506)= tmp*4.6189D+04
      zlm(507)=-tmp*1.09395D+05
      zlm(508)= tmp*9.0090D+04
      zlm(509)=-tmp*3.0030D+04
      zlm(510)= tmp*3.465D+03
      zlm(511)=-tmp*6.3D+01
! L=10, ML= 1
      zlm(512)= zlm(501)
      zlm(513)= zlm(502)
      zlm(514)= zlm(503)
      zlm(515)= zlm(504)
      zlm(516)= zlm(505)
! L=10, ML= 2
      zlm(517)= zlm(491)*two
      zlm(518)= zlm(493)*two
      zlm(519)= zlm(495)*two
      zlm(520)= zlm(497)*two
      zlm(521)= zlm(499)*two
! L=10, ML= 3
      zlm(522)=-zlm(484)
      zlm(523)=-zlm(483)
      zlm(524)=-zlm(486)
      zlm(525)=-zlm(485)
      zlm(526)=-zlm(488)
      zlm(527)=-zlm(487)
      zlm(528)=-zlm(490)
      zlm(529)=-zlm(489)
! L=10, ML= 4
      zlm(530)= zlm(471)*four
      zlm(531)=-zlm(530)
      zlm(532)= zlm(474)*four
      zlm(533)=-zlm(532)
      zlm(534)= zlm(477)*four
      zlm(535)=-zlm(534)
      zlm(536)= zlm(480)*four
      zlm(537)=-zlm(536)
! L=10, ML= 5
      zlm(538)= zlm(462)
      zlm(539)= zlm(463)
      zlm(540)= zlm(464)
      zlm(541)= zlm(465)
      zlm(542)= zlm(466)
      zlm(543)= zlm(467)
      zlm(544)= zlm(468)
      zlm(545)= zlm(469)
      zlm(546)= zlm(470)
! L=10, ML= 6
      zlm(547)= zlm(450)*six
      zlm(548)=-zlm(450)*p20
      zlm(549)= zlm(547)
      zlm(550)= zlm(454)*six
      zlm(551)=-zlm(454)*p20
      zlm(552)= zlm(550)
      zlm(553)= zlm(458)*six
      zlm(554)=-zlm(458)*p20
      zlm(555)= zlm(553)
! L=10, ML= 7
      zlm(556)=-zlm(445)
      zlm(557)=-zlm(444)
      zlm(558)=-zlm(443)
      zlm(559)=-zlm(442)
      zlm(560)=-zlm(449)
      zlm(561)=-zlm(448)
      zlm(562)=-zlm(447)
      zlm(563)=-zlm(446)
! L=10, ML= 8
      zlm(564)= zlm(432)*eight
      zlm(565)=-zlm(432)*5.6D+01
      zlm(566)=-zlm(565)
      zlm(567)=-zlm(564)
      zlm(568)= zlm(437)*eight
      zlm(569)=-zlm(437)*5.6D+01
      zlm(570)=-zlm(569)
      zlm(571)=-zlm(568)
! L=10, ML= 9
      zlm(572)= zlm(431)
      zlm(573)= zlm(430)
      zlm(574)= zlm(429)
      zlm(575)= zlm(428)
      zlm(576)= zlm(427)
! L=10, ML= 10
      zlm(577)= zlm(421)*ten
      zlm(578)=-zlm(421)*1.20D+02
      zlm(579)= zlm(421)*2.52D+02
      zlm(580)= zlm(578)
      zlm(581)= zlm(577)
!
      return
end


!-----------------------------------------
  subroutine ecpxyz(xyzintecp,maxecpdim)
!-----------------------------------------
!
! Set up angular parts for ECP calculations
!
      implicit none
      integer,intent(in) :: maxecpdim
      integer :: maxdim, i, j, k
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00, three=3.0D+00
      real(8),parameter :: pi4=1.256637061435917D+01
      real(8),intent(out) :: xyzintecp(0:24,0:24,0:24)
      real(8) :: factor(0:72), denom, xx, yy, zz
!
      maxdim= maxecpdim*4
!
      if(maxdim > 24) then
        write(*,'(" Error! This program supports up to i function in Subroutine ecpxyz!")')
        call exit
      endif
!
      denom= three
      do i= 0,3*maxdim-2,2
        factor(i)=one/denom
        denom= denom+two
      enddo
      xyzintecp(:,:,:)= zero
      xyzintecp(0,0,0)= pi4

      xx=-one
      do i= 0,maxdim,2
        yy=-one
        do j= 0,maxdim,2
          zz=-one
          do k= 0,maxdim,2
            if(i > 0) then
              xyzintecp(k,j,i)= xx*xyzintecp(k,j,i-2)*factor(k+j+i-2)
            elseif(j > 0) then
              xyzintecp(k,j,i)= yy*xyzintecp(k,j-2,i)*factor(k+j-2+i)
            elseif(k > 0) then
              xyzintecp(k,j,i)= zz*xyzintecp(k-2,j,i)*factor(k-2+j+i)
            endif
            zz= zz+two
          enddo
          yy= yy+two
        enddo
        xx= xx+two
      enddo
!
      return
end


!------------------------------------------------------------------------
  subroutine ecpangint1(term1ecp,label1ecp,num1ecp,xyzintecp,maxecpdim)
!------------------------------------------------------------------------
!
! Generate angular parts of term1 for ECP calculations
!
      use modecp, only : nx, ny, nz, lmf, lmx, lmy, lmz, zlm, nterm1
      implicit none
      integer,parameter :: itri(0:6)=(/0,1,3,6,10,15,21/), mbasis(0:7)=(/1,2,5,11,21,36,57,85/)
      integer,intent(in) :: maxecpdim
      integer,intent(out) :: label1ecp(9,nterm1), num1ecp(*)
      integer :: icount, ibasis, jbasis, mx1, my1, mz1, mx2, my2, mz2, lmax, lambda, mu, lmindex
      integer :: kx, ky, kz, kxp, kyp, kzp, ksum, lm, ijx, ijy, ijz, ii, iit, jj, ijindex
      real(8),parameter :: combi(0:27)= &
&     (/1.0D+00, 1.0D+00,1.0D+00,  1.0D+00,2.0D+00,1.0D+00, 1.0D+00,3.0D+00,3.0D+00,1.0D+00, &
&       1.0D+00,4.0D+00,6.0D+00,4.0D+00,1.0D+00, 1.0D+00,5.0D+00,10.0D+00,10.0D+00,5.0D+00, &
&       1.0D+00, 1.0D+00,6.0D+00,15.0D+00,20.0D+00,15.0D+00,6.0D+00,1.0D+00/)
      real(8),parameter :: zero=0.0D+00, tol=1.0D-15
      real(8),intent(in) :: xyzintecp(0:24,0:24,0:24)
      real(8),intent(out) :: term1ecp(nterm1)
      real(8) :: angint1
!
      icount= 0
      num1ecp(1)= 1
!
      do ii= 0,maxecpdim
        do ibasis= mbasis(ii),mbasis(ii+1)-1
          mx1= itri(nx(ibasis))
          my1= itri(ny(ibasis))
          mz1= itri(nz(ibasis))
          iit= ibasis*(ibasis-1)/2
          do jj= 0,ii
            do jbasis= mbasis(jj),min(mbasis(jj+1)-1,ibasis)
              mx2= itri(nx(jbasis))
              my2= itri(ny(jbasis))
              mz2= itri(nz(jbasis))
              lmax= ii+jj
              ijindex= iit+jbasis
              do lambda= 0,lmax
                do mu=-lambda,lambda
                  lmindex= lambda*(lambda+1)-mu+1
                  do kx= 0,nx(ibasis)
                    do ky= 0,ny(ibasis)
                      do kz= 0,nz(ibasis)
                        do kxp= 0,nx(jbasis)
                          do kyp= 0,ny(jbasis)
                            do kzp= 0,nz(jbasis)
                              ksum=kx+ky+kz+kxp+kyp+kzp
                              if((mod(lambda+ksum,2) /= 1).and.(lambda <= ksum)) then
                                angint1= zero
                                do lm= lmf(lmindex),lmf(lmindex+1)-1
                                  ijx= kx+kxp+lmx(lm)
                                  ijy= ky+kyp+lmy(lm)
                                  ijz= kz+kzp+lmz(lm)
                                  angint1= angint1+zlm(lm)*xyzintecp(ijz,ijy,ijx)
                                enddo
                                if(abs(angint1) > tol)then
                                  icount= icount+1
                                  label1ecp(1,icount)= ksum
                                  label1ecp(2,icount)= lambda
                                  label1ecp(3,icount)= mu
                                  label1ecp(4,icount)= kx
                                  label1ecp(5,icount)= ky
                                  label1ecp(6,icount)= kz
                                  label1ecp(7,icount)= kxp
                                  label1ecp(8,icount)= kyp
                                  label1ecp(9,icount)= kzp
                                  term1ecp(icount)= combi(mx1+kx) *combi(my1+ky) *combi(mz1+kz) * &
&                                                   combi(mx2+kxp)*combi(my2+kyp)*combi(mz2+kzp)* &
&                                                   angint1
                                endif
                              endif
                            enddo
                          enddo
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              enddo
              num1ecp(ijindex+1)= icount+1
            enddo
          enddo
        enddo
      enddo
      if(icount > nterm1) then
        write(*,'(" Error! Icount exceeds",i8," in Subroutine ecpangint1!")')nterm1
        call exit
      endif
!
      return
end


!----------------------------------------------------------------------------------------
  subroutine ecpangint2(term2ecp,label2ecp,num2ecp,xyzintecp,maxecpdim,llmax,numtbasis)
!----------------------------------------------------------------------------------------
!
! Generate angular parts of term2 for ECP calculations
!
      use modecp, only : nx, ny, nz, lmf, lmx, lmy, lmz, zlm, nterm2
      implicit none
      integer,parameter :: itri(0:6)=(/0,1,3,6,10,15,21/), mbasis(0:7)=(/1,2,5,11,21,36,57,85/)
      integer,intent(in) :: maxecpdim, llmax, numtbasis
      integer,intent(out) :: label2ecp(6,nterm2), num2ecp(*)
      integer :: icount, lambda, mu, lmindex1, lmlabel, ii, ltotal, ibasis
      integer :: mx, my, mz, ll, mm, lmindex2, kx, ky, kz, ksum, lm1, ijx, ijy, ijz, lm2
      real(8),parameter :: combi(0:27)= &
&     (/1.0D+00, 1.0D+00,1.0D+00,  1.0D+00,2.0D+00,1.0D+00, 1.0D+00,3.0D+00,3.0D+00,1.0D+00, &
&       1.0D+00,4.0D+00,6.0D+00,4.0D+00,1.0D+00, 1.0D+00,5.0D+00,10.0D+00,10.0D+00,5.0D+00, &
&       1.0D+00, 1.0D+00,6.0D+00,15.0D+00,20.0D+00,15.0D+00,6.0D+00,1.0D+00/)
      real(8),parameter :: zero=0.0D+00, tol=1.0D-15
      real(8),intent(in) :: xyzintecp(0:24,0:24,0:24)
      real(8),intent(out) :: term2ecp(nterm2)
      real(8) :: angint2, tmp
!
      icount= 0
      num2ecp(1)= 1
!
      do lambda= 0,llmax-1
        do mu=-lambda,lambda
          lmindex1= lambda*(lambda+1)+mu+1
          lmlabel=(lambda*(lambda+1)+mu)*numtbasis
          do ii= 0,maxecpdim
            ltotal= lambda+ii
            do ibasis= mbasis(ii),mbasis(ii+1)-1
              mx= itri(nx(ibasis))
              my= itri(ny(ibasis))
              mz= itri(nz(ibasis))
              do ll= 0,ltotal
                do mm=-ll,ll
                  lmindex2= ll*(ll+1)-mm+1
                  do kx= 0,nx(ibasis)
                    do ky= 0,ny(ibasis)
                      do kz= 0,nz(ibasis)
                      ksum= kx+ky+kz
                      angint2= zero
                      do lm1= lmf(lmindex1),lmf(lmindex1+1)-1
                        ijx= kx+lmx(lm1)
                        ijy= ky+lmy(lm1)
                        ijz= kz+lmz(lm1)
                        tmp= zero
                        do lm2= lmf(lmindex2),lmf(lmindex2+1)-1
                          tmp= tmp+zlm(lm2)*xyzintecp(ijz+lmz(lm2),ijy+lmy(lm2),ijx+lmx(lm2))
                        enddo
                        angint2= angint2+zlm(lm1)*tmp
                      enddo
                      if(abs(angint2) > tol) then
                        icount= icount+1
                        label2ecp(1,icount)= ll
                        label2ecp(2,icount)= ksum
                        label2ecp(3,icount)= mm
                        label2ecp(4,icount)= kx
                        label2ecp(5,icount)= ky
                        label2ecp(6,icount)= kz
                        term2ecp(icount)= combi(mx+kx)*combi(my+ky)*combi(mz+kz)*angint2
                      endif
                      enddo
                    enddo
                  enddo
                enddo
              enddo
              num2ecp(lmlabel+ibasis+1)= icount+1
            enddo
          enddo
        enddo
      enddo
      if(icount > nterm2) then
        write(*,'(" Error! Icount exceeds",i8," in Subroutine ecpangint2!")')nterm2
        call exit
      endif
!
      return
end


!----------------------------------------------------------------------
  subroutine ecpangint0(term0ecp,xyzintecp,maxecpdim,llmax,numtbasis)
!----------------------------------------------------------------------
!
! Generate angular parts of term0 for ECP calculations with same center
!
      use modecp, only : nx, ny, nz, lmf, lmx, lmy, lmz, zlm
      implicit none
      integer,parameter :: mbasis(0:7)=(/1,2,5,11,21,36,57,85/)
      integer,intent(in) :: maxecpdim, llmax, numtbasis
      integer :: lambda, mu, lmindex, lmlabel, i, ibasis, lm, ijx, ijy, ijz
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: xyzintecp(0:24,0:24,0:24)
      real(8),intent(out) :: term0ecp(*)
      real(8) :: angint0
!
      do lambda= 0,llmax-1
        do mu=-lambda,lambda
          lmindex= lambda*(lambda+1)+mu+1
          lmlabel=(lambda*(lambda+1)+mu)*numtbasis
          do i= 0,maxecpdim
            do ibasis= mbasis(i),mbasis(i+1)-1
              angint0= zero
              do lm= lmf(lmindex),lmf(lmindex+1)-1
                ijx= nx(ibasis)+lmx(lm)
                ijy= ny(ibasis)+lmy(lm)
                ijz= nz(ibasis)+lmz(lm)
                angint0= angint0+zlm(lm)*xyzintecp(ijz,ijy,ijx)
              enddo
              term0ecp(lmlabel+ibasis)= angint0
            enddo
          enddo
        enddo
      enddo
!
      return
end


!-----------------------
  function dawson(val)
!-----------------------
!
! Evaluate the Dawson function
!
      implicit none
      integer :: igrid, last, i
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, ten=1.0D+01
      real(8),parameter :: griddaw=5.0D+00
      real(8),intent(in) :: val
      real(8) :: dawson, val1, val2, val3, val4, val5, val6, val7, val8, val9
      real(8) :: val2inv, factor, total
      real(8) :: dawgrid(0:9,0:49)
      data dawgrid / &
&    6.051722509208561D-16,   9.999999999994005D-01,   9.777221521513968D-11, &
&   -6.666666728319638D-01,   1.965028896818846D-07,   2.666630955078337D-01, &
&    3.902825044190826D-05,  -7.644842301092071D-02,   9.777389039626225D-04, &
&    1.519711057206845D-02,   2.953748577202614D-08,   9.999989604795872D-01, &
&    1.632904172907173D-05,  -6.668173339622047D-01,   9.033269735325170D-04, &
&    2.629969013292320D-01,   1.018881107249708D-02,  -9.510163851001064D-02, &
&    2.186581860556165D-02,   4.206328735375168D-03,   5.022734342717940D-06, &
&    9.998964691434666D-01,   9.591681774814235D-04,  -6.719264553715743D-01, &
&    1.890290652354323D-02,   2.201777674225539D-01,   7.906931394260616D-02, &
&   -1.674426202438445D-01,   6.691761616978956D-02,  -8.476321032168930D-03, &
&    3.669932078830077D-05,   9.994086664526127D-01,   4.313802985030914D-03, &
&   -6.854492931066708D-01,   5.411705735859217D-02,   1.587484871726615D-01, &
&    1.508505021326560D-01,  -2.216150267576290D-01,   9.087352189345086D-02, &
&   -1.320491480520023D-02,  -9.240529554626878D-04,   1.009683221114372D+00, &
&   -4.461635662652316D-02,  -5.492383023037898D-01,  -1.901858785183352D-01, &
&    4.515560484687740D-01,  -8.369896234661265D-02,  -1.005130134305425D-01, &
&    5.429770655133835D-02,  -8.280903040450580D-03,  -1.029775943821545D-02, &
&    1.092399269442080D+00,  -3.697768628624183D-01,   1.981878620126198D-01, &
&   -1.297394353718835D+00,   1.547787220212697D+00,  -8.091518439150200D-01, &
&    2.089255992186689D-01,  -2.290165554070681D-02,   3.020315018984405D-04, &
&   -3.635700884765439D-02,   1.287744381657310D+00,  -1.021778458451531D+00, &
&    1.469941481761009D+00,  -2.894980361086176D+00,   2.888159254985652D+00, &
&   -1.560225529053189D+00,   4.799668635758824D-01,  -8.005976523151092D-02, &
&    5.668640259844556D-03,  -2.430672682332134D-02,   1.218214029581446D+00, &
&   -8.454866418703851D-01,   1.212966595174907D+00,  -2.658748167631678D+00, &
&    2.747139007450606D+00,  -1.506200145273536D+00,   4.674338456679980D-01, &
&   -7.853618452992746D-02,   5.604485954876450D-03,   2.086326626721770D-01, &
&   -7.659095371377669D-02,   2.356228131334686D+00,  -3.409613779527687D+00, &
&    1.635769789964874D+00,   8.474024580105899D-02,  -4.047464050401199D-01, &
&    1.742063022481715D-01,  -3.295397297470270D-02,   2.452076490858092D-03, &
&    7.819875729224917D-01,  -2.937321717501164D+00,   8.705241547587803D+00, &
&   -1.163605938718322D+01,   8.493710605657423D+00,  -3.729830062303937D+00, &
&    1.010943062539432D+00,  -1.638322120981495D-01,   1.416996484227287D-02, &
&   -4.700124848768664D-04,   1.304119931791932D+00,  -5.308704826333093D+00, &
&    1.349480033480972D+01,  -1.728228980794192D+01,   1.277511030971419D+01, &
&   -5.895380025040220D+00,   1.741581652446646D+00,  -3.223913333975149D-01, &
&    3.425315766198099D-02,  -1.601172584406407D-03,   8.453967279834436D-01, &
&   -3.473857303563860D+00,   1.023282373069887D+01,  -1.389938194430536D+01, &
&    1.051971338675129D+01,  -4.892911707033860D+00,   1.444531320532356D+00, &
&   -2.658060541288357D-01,   2.796549043503804D-02,  -1.290655862378755D-03, &
&   -1.168868023013557D+00,   4.051680384856491D+00,  -2.268990395019065D+00, &
&   -1.778797427179115D+00,   2.962044850454449D+00,  -1.749791207253024D+00, &
&    5.726770875040680D-01,  -1.102658329522576D-01,   1.177126126343049D-02, &
&   -5.409385300534173D-04,  -3.983678090018918D+00,   1.380643270044792D+01, &
&   -1.729954790908351D+01,   1.173647405757583D+01,  -4.853561872641671D+00, &
&    1.264482812164978D+00,  -2.026489282503666D-01,   1.798868732943707D-02, &
&   -6.095550168767980D-04,  -9.547785170957528D-06,  -5.814425354970969D+00, &
&    1.973235274865224D+01,  -2.582728197811491D+01,   1.889725061452376D+01, &
&   -8.720183404777921D+00,   2.656809882256313D+00,  -5.369883591179374D-01, &
&    6.961546456260244D-02,  -5.261176970280754D-03,   1.767788777026804D-04, &
&   -5.329689633295526D+00,   1.831861036236389D+01,  -2.399576707118949D+01, &
&    1.751397831213208D+01,  -8.049000441611372D+00,   2.439847689603503D+00, &
&   -4.902670710796377D-01,   6.315275734302311D-02,  -4.740151937244599D-03, &
&    1.581271836076713D-04,  -2.742874229129735D+00,   1.105990371185430D+01, &
&   -1.494089841574973D+01,   1.092324271927255D+01,  -4.964303598486863D+00, &
&    1.477099885778468D+00,  -2.898963332247264D-01,   3.633739182918682D-02, &
&   -2.646216460201094D-03,   8.543732483697165D-05,   5.334358295001977D-01, &
&    2.380065505895000D+00,  -4.718367967869862D+00,   3.898594616817992D+00, &
&   -1.860403980125197D+00,   5.625617963018695D-01,  -1.102134285396532D-01, &
&    1.363727777246505D-02,  -9.729401144281072D-04,   3.060634340168288D-05, &
&    3.080620069423487D+00,  -4.005453782311422D+00,   2.397694500496021D+00, &
&   -7.283035437468496D-01,   7.398635616053244D-02,   2.330434398922989D-02, &
&   -9.972367942442934D-03,   1.656144541521961D-03,  -1.374241748923812D-04, &
&    4.705353090582574D-06,   4.300353044663805D+00,  -6.909569738185462D+00, &
&    5.471347664874574D+00,  -2.626260134921639D+00,   8.275246738550115D-01, &
&   -1.761785998740552D-01,   2.523933144701089D-02,  -2.340119146346428D-03, &
&    1.271874623248163D-04,  -3.083104605409867D-06,   4.396530707222811D+00, &
&   -7.134196626551821D+00,   5.704257216451345D+00,  -2.766989300405958D+00, &
&    8.821357230714314D-01,  -1.902941758333474D-01,   2.766965259236612D-02, &
&   -2.608903835934439D-03,   1.445153102345324D-04,  -3.579244353753892D-06, &
&    3.896858160212054D+00,  -6.065851873388711D+00,   4.688902043479053D+00, &
&   -2.203994088320108D+00,   6.814255243986297D-01,  -1.425844799207229D-01, &
&    2.010800406935943D-02,  -1.838346719741425D-03,   9.870384880496008D-05, &
&   -2.368575908529228D-06,   3.246922643039968D+00,  -4.735897812566311D+00, &
&    3.479189852464805D+00,  -1.562037274177388D+00,   4.623945424588869D-01, &
&   -9.275624721235535D-02,   1.254985334572399D-02,  -1.101238261520193D-03, &
&    5.676415575057387D-05,  -1.307863413499134D-06,   2.675504458023278D+00, &
&   -3.616587250649823D+00,   2.504601603004815D+00,  -1.066969118221491D+00, &
&    3.007058498254446D-01,  -5.754679238803369D-02,   7.437700565744401D-03, &
&   -6.240180589474892D-04,   3.077414469741080D-05,  -6.786944578519707D-07, &
&    2.243099420626957D+00,  -2.804674807408287D+00,   1.826964650222504D+00, &
&   -7.370149564337481D-01,   1.974115492186115D-01,  -3.598626016422626D-02, &
&    4.437136966000534D-03,  -3.555380013897935D-04,   1.675926942372285D-05, &
&   -3.535061057294843D-07,   1.932029874176854D+00,  -2.243964191783463D+00, &
&    1.377717335005858D+00,  -5.270258066140214D-01,   1.343055779563454D-01, &
&   -2.334179536965425D-02,   2.747914965941448D-03,  -2.104488646833919D-04, &
&    9.489056293487565D-06,  -1.915773367624041D-07,   1.706708849350314D+00, &
&   -1.853496313508908D+00,   1.076949899974014D+00,  -3.918691452174971D-01, &
&    9.525721411066638D-02,  -1.582000633173914D-02,   1.781881255989836D-03, &
&   -1.306821199166410D-04,   5.646571822689784D-06,  -1.093031483485892D-07, &
&    1.537669196618148D+00,  -1.571461484279969D+00,   8.677913416794271D-01, &
&   -3.013779063895062D-01,   7.008664461576448D-02,  -1.115202231884632D-02, &
&    1.204695757535073D-03,  -8.479858763893169D-05,   3.518648920327665D-06, &
&   -6.543865306252688D-08,   1.405549318697144D+00,  -1.358932771279227D+00, &
&    7.158335591164556D-01,  -2.379933527355135D-01,   5.308868832732922D-02, &
&   -8.112838465920787D-03,   8.423992317003719D-04,  -5.703190720497197D-05, &
&    2.277172150589343D-06,  -4.076640715826195D-08,   1.298576035768416D+00, &
&   -1.192809905599154D+00,   6.011669931292463D-01,  -1.918194294121542D-01, &
&    4.113486133628721D-02,  -6.049542193356675D-03,   6.049555543509496D-04, &
&   -3.946441816132111D-05,   1.518924456429258D-06,  -2.621967464720331D-08, &
&    1.209577830183879D+00,  -1.059221212897871D+00,   5.120400290656999D-01, &
&   -1.571298707954724D-01,   3.245452949501371D-02,  -4.601384001133819D-03, &
&    4.438769850238281D-04,  -2.794561679474710D-05,   1.038387435358890D-06, &
&   -1.730927703431548D-08,   1.133997350075610D+00,  -9.494408737898824D-01, &
&    4.411658021914067D-01,  -1.304366957767376D-01,   2.599117585807237D-02, &
&   -3.557972489308219D-03,   3.315734790826391D-04,  -2.017462082284952D-05, &
&    7.246921372578375D-07,  -1.168083844809918D-08,   1.068780177926012D+00, &
&   -8.576790733538303D-01,   3.837794250426504D-01,  -1.095003008494338D-01, &
&    2.108052719188805D-02,  -2.790056105619201D-03,   2.515113598636368D-04, &
&   -1.480820888493274D-05,   5.148527578287036D-07,  -8.033830684925624D-09, &
&    1.011779626492335D+00,  -7.799125024116445D-01,   3.366217910576227D-01, &
&   -9.281801625720301D-02,   1.728649366216113D-02,  -2.214770949664605D-03, &
&    1.933544248618440D-04,  -1.102846916299806D-05,   3.715460199775798D-07, &
&   -5.618835277160696D-09,   9.614096479074705D-01,  -7.132165323271270D-01, &
&    2.973689581221718D-01,  -7.934126227493824D-02,   1.431181444060736D-02, &
&   -1.777017303985521D-03,   1.504053144171016D-04,  -8.319406490500501D-06, &
&    2.718617342711667D-07,  -3.988496802653806D-09,   9.165526309123725D-01, &
&   -6.555196224552808D-01,   2.643839356220485D-01,  -6.834051008923832D-02, &
&    1.195314538712647D-02,  -1.439850416169431D-03,   1.182719583461593D-04, &
&   -6.350585515470334D-06,   2.014901001305661D-07,  -2.870524681794359D-09, &
&    8.762449238894968D-01,  -6.051156577245481D-01,   2.363694886681162D-01, &
&   -5.925729795070593D-02,   1.005977937101519D-02,  -1.176725686645202D-03, &
&    9.389268779118348D-05,  -4.898415677199682D-06,   1.510293060240204D-07, &
&   -2.091176767916272D-09,   8.397973374736750D-01,  -5.607717546712178D-01, &
&    2.123900746196760D-01,  -5.169276657599643D-02,   8.525650051352768D-03, &
&   -9.692951718191823D-04,   7.519387477589534D-05,  -3.814758613123273D-06, &
&    1.143932740559919D-07,  -1.540667804508020D-09,   8.066453081016311D-01, &
&   -5.214997953854644D-01,   1.917127957710592D-01,  -4.534177013409804D-02, &
&    7.271568748050543D-03,  -8.041981986548870D-04,   6.070341686089757D-05, &
&   -2.997121855987458D-06,   8.747939440781839D-08,  -1.146908820242620D-09, &
&    7.763328700748650D-01,  -4.865131195188174D-01,   1.737645584369129D-01, &
&   -3.997049326403647D-02,   6.238170884938409D-03,  -6.716462489699137D-04, &
&    4.936811186568123D-05,  -2.373940361615341D-06,   6.749309265777579D-08, &
&   -8.620120298597762D-10,   7.484872168611639D-01,  -4.551777208971689D-01, &
&    1.580916582586361D-01,  -3.539751489571131D-02,   5.380378405930295D-03, &
&   -5.643726412581538D-04,   4.042412621397579D-05,  -1.894535342555393D-06, &
&    5.250281623647732D-08,  -6.536814794056658D-10,   7.228002002617220D-01, &
&   -4.269770225009569D-01,   1.443309279425442D-01,  -3.148047647059620D-02, &
&    4.663564236690544D-03,  -4.769181755117127D-04,   3.331061449736308D-05, &
&   -1.522556290172317D-06,   4.115564531122071D-08,  -4.998334049864280D-10, &
&    6.990146269097030D-01,  -4.014859925399112D-01,   1.321887743925377D-01, &
&   -2.810653190187287D-02,   4.060849619196801D-03,  -4.051368747339830D-04, &
&    2.761111481851752D-05,  -1.231622827617666D-06,   3.249232324267357D-08, &
&   -3.851739792444945D-10,   6.769139442201537D-01,  -3.783518223820482D-01, &
&    1.214257096016929D-01,  -2.518540056502214D-02,   3.551170329362119D-03, &
&   -3.458487212194843D-04,   2.301316857550481D-05,  -1.002382719404352D-06, &
&    2.582501734241042D-08,  -2.989865561703131D-10,   6.563143660892078D-01, &
&   -3.572792987091244D-01,   1.118447759141316D-01,  -2.264424363481021D-02, &
&    3.117873830445191D-03,  -2.965922633845162D-04,   1.928010426958492D-05, &
&   -8.204973417272264D-07,   2.065535195607614D-08,  -2.336795023974850D-10, &
&    6.370587839081798D-01,  -3.380195885781189D-01,   1.032827760925571D-01, &
&   -2.042383776876847D-02,   2.747688779923439D-03,  -2.554461321436570D-04, &
&    1.623106686867063D-05,  -6.752440992628909D-07,   1.661870682317888D-08, &
&   -1.838200346615760D-10,   6.190120033975020D-01,  -3.203615420730855D-01, &
&    9.560355253336551D-02,  -1.847568500303102D-02,   2.429959362224784D-03, &
&   -2.208987753263688D-04,   1.372671585753052D-05,  -5.585351186952662D-07, &
&    1.344589462985757D-08,  -1.454832089595699D-10,   6.020569792050626D-01, &
&   -3.041248766244393D-01,   8.869278306835553D-02,  -1.675980697559348D-02, &
&    2.156070393944811D-03,  -1.917523666799123D-04,   1.165887107556909D-05, &
&   -4.642205187486950D-07,   1.093649208022802D-08,  -1.158080645769749D-10, &
&    5.860918095088875D-01,  -2.891547845930740D-01,   8.245391253436112D-02, &
&   -1.524304526077679D-02,   1.919011474773674D-03,  -1.670512229576893D-04, &
&    9.942940784081151D-06,  -3.875885889643397D-07,   8.940082061173694D-09, &
&   -9.269166157280035D-11,   5.710273159428749D-01,  -2.753176292173221D-01, &
&    7.680494417517188D-02,  -1.389774006501874D-02,   1.713043581841007D-03, &
&   -1.460279888984735D-04,   8.512332244579305D-06,  -3.250037532305761D-07, &
&    7.342928218127729D-09,  -7.457599172953888D-11/
!
      val1= abs(val)
      if(val1 < ten) then
        igrid= aint(val1*griddaw)
        val2= val1*val1
        val3= val1*val1*val1
        val4= val2*val2
        val5= val2*val3
        val6= val3*val3
        val7= val3*val4
        val8= val4*val4
        val9= val4*val4*val1
        dawson= dawgrid(0,igrid)     +dawgrid(1,igrid)*val1+dawgrid(2,igrid)*val2 &
&              +dawgrid(3,igrid)*val3+dawgrid(4,igrid)*val4+dawgrid(5,igrid)*val5 &
&              +dawgrid(6,igrid)*val6+dawgrid(7,igrid)*val7+dawgrid(8,igrid)*val8 &
&              +dawgrid(9,igrid)*val9
        if(val < zero) dawson=-dawson
      else
        val2inv= one/(val1*val1)
        last= 9
        factor=(last+half)*val2inv
        total= factor
        do i= last,1,-1
          factor= factor-val2inv
          total = factor*(one+total)
        enddo
        dawson= val*factor*(one+total)
      endif
!
      return
end


!-----------------------
  function dawerf(val)
!-----------------------
!
! Evaluate the Dawson-Error function
!
      implicit none
      integer :: igrid, last, i
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, ten=1.0D+01
      real(8),parameter :: griddaw=5.0D+00
      real(8),intent(in) :: val
      real(8) :: dawerf, val1, val2, val3, val4, val5, val6, val7, val8, val9
      real(8) :: val2inv, factor, total
      real(8) :: dawerfgrid(0:9,0:49)
      data dawerfgrid / &
&   -1.464050227112128D-15,   1.467214183101818D-12,   5.641895833050387D-01, &
&    1.558305800769015D-08,  -3.761268974805824D-01,   9.539122392875117D-06, &
&    1.440726531464679D-01,   7.747133582156654D-04,  -4.274238191602858D-02, &
&    8.037433984828787D-03,  -9.951223039681468D-09,   3.665620582051934D-07, &
&    5.641835068621823D-01,   5.980782869308470D-05,  -3.765143811500408D-01, &
&    1.739802654227826D-03,   1.386911329865945D-01,   1.214475306151317D-02, &
&   -5.768347019627200D-02,   1.737017358803486D-02,   3.505730303251949D-06, &
&   -6.903611739405604D-05,   5.647934766198526D-01,  -3.073493000942437D-03, &
&   -3.661454290548012D-01,  -2.118740661715355D-02,   1.725643524039445D-01, &
&   -2.009535893707317D-02,  -3.975174841729590D-02,   1.293250194958627D-02, &
&    1.753521010416664D-04,  -2.509136076865635D-03,   5.802654485366705D-01, &
&   -6.059720330548457D-02,  -2.278842721571183D-01,  -2.440739289353691D-01, &
&    4.136483694069031D-01,  -1.888713464726069D-01,   2.966074254247206D-02, &
&    1.521467396299521D-04,   1.493399476360794D-03,  -1.710408420590437D-02, &
&    6.523669837322738D-01,  -2.691935415762747D-01,   1.616337380445307D-01, &
&   -7.309645117937897D-01,   8.210629072070301D-01,  -4.089382992885937D-01, &
&    9.928847617461442D-02,  -9.678922904382312D-03,   7.022349491783161D-04, &
&   -1.140029281594338D-02,   6.351255670624923D-01,  -2.418066058845062D-01, &
&    1.396558590557061D-01,  -7.280085092378150D-01,   8.317740595851261D-01, &
&   -4.189875032602776D-01,   1.031772638668831D-01,  -1.026592845497949D-02, &
&   -3.507828429386371D-02,   2.514441281612858D-01,  -2.243057271902791D-01, &
&    1.399931801042554D+00,  -1.879615077893551D+00,   9.304166624505790D-01, &
&   -7.777252579099405D-02,  -9.776713912647873D-02,   3.688751816961262D-02, &
&   -4.175226370687775D-03,  -1.635192995852709D-01,   1.072530607692978D+00, &
&   -2.560291425726306D+00,   5.281908470059739D+00,  -6.032401344817990D+00, &
&    3.896121221224385D+00,  -1.491676614240306D+00,   3.361654360188453D-01, &
&   -4.090475833958504D-02,   2.031481573663880D-03,  -2.969114851113720D-01, &
&    1.833603926731921D+00,  -4.491789807629176D+00,   8.143636663140860D+00, &
&   -8.760254633771366D+00,   5.630966936454882D+00,  -2.227785302457666D+00, &
&    5.371041861773432D-01,  -7.292479631053337D-02,   4.300877192829566D-03, &
&   -1.170572792872378D-03,   3.848954456804926D-01,  -1.336433123369666D+00, &
&    4.132965957640049D+00,  -5.481656626887961D+00,   3.843409548629058D+00, &
&   -1.577751169995361D+00,   3.850753397958487D-01,  -5.217397376886755D-02, &
&    3.041464194679571D-03,   1.174905657055873D+00,  -4.883509768159756D+00, &
&    9.159562518680080D+00,  -8.072949374706283D+00,   3.649353826728024D+00, &
&   -7.134423639453488D-01,  -6.066654191095896D-02,   6.016752490052443D-02, &
&   -1.155594302949185D-02,   7.831369932113179D-04,   2.807721647232346D+00, &
&   -1.157888576964923D+01,   2.136828102040118D+01,  -2.106630531361268D+01, &
&    1.254393798005916D+01,  -4.774847360791843D+00,   1.176340252893272D+00, &
&   -1.821681364947661D-01,   1.615253878701492D-02,  -6.256883756474562D-04, &
&    3.348333345889982D+00,  -1.366040981604734D+01,   2.492960224723117D+01, &
&   -2.462000868400332D+01,   1.482320891370940D+01,  -5.749289456829261D+00, &
&    1.454034752069000D+00,  -2.330356890766367D-01,   2.158734562461293D-02, &
&   -8.837375538626612D-04,   1.419529100885629D+00,  -7.035253702007297D+00, &
&    1.481226734600577D+01,  -1.560431657184057D+01,   9.656736117340066D+00, &
&   -3.774835259282677D+00,   9.508132419060235D-01,  -1.505577918487977D-01, &
&    1.369900313160904D-02,  -5.483047252040330D-04,  -2.544030075609960D+00, &
&    5.695392362927285D+00,  -3.367263672606511D+00,  -4.553991907628051D-01, &
&    1.538810907283111D+00,  -8.737041337394308D-01,   2.593827681045920D-01, &
&   -4.458505602238055D-02,   4.221228573272723D-03,  -1.714398979549848D-04, &
&   -6.348319107215338D+00,   1.714154675262704D+01,  -1.867786552351831D+01, &
&    1.149466787193640D+01,  -4.458975430819412D+00,   1.133760418179475D+00, &
&   -1.886815669411439D-01,   1.972420037565071D-02,  -1.164553455654676D-03, &
&    2.908398028937653D-05,  -7.854162583305704D+00,   2.142212197657538D+01, &
&   -2.408684341726629D+01,   1.548232235238074D+01,  -6.349167351893407D+00, &
&    1.731172160849376D+00,  -3.145799647601360D-01,   3.678306401855272D-02, &
&   -2.513096665959857D-03,   7.647153503181144D-05,  -6.587158759431202D+00, &
&    1.809711106668978D+01,  -2.020811667944246D+01,   1.284255197793339D+01, &
&   -5.194060318853190D+00,   1.394154820937142D+00,  -2.490172759971632D-01, &
&    2.858253201482559D-02,  -1.914670041028719D-03,   5.705978923133756D-05, &
&   -3.695796883540154D+00,   1.087367011914940D+01,  -1.218591014395560D+01, &
&    7.644375518853737D+00,  -3.028284844164388D+00,   7.924612082670810D-01, &
&   -1.375526789103503D-01,   1.530541438488127D-02,  -9.919323992952312D-04, &
&    2.855217532838284D-05,  -7.482906093884771D-01,   3.883066123961691D+00, &
&   -4.815793457346782D+00,   3.110883017076421D+00,  -1.235253709919397D+00, &
&    3.196005001244074D-01,  -5.440141038991109D-02,   5.903846745328878D-03, &
&   -3.717356719999930D-04,   1.036535141718575D-05,   1.294527932794106D+00, &
&   -7.257393636901485D-01,  -1.937132305125310D-01,   4.064505512011646D-01, &
&   -2.178302178553769D-01,   6.438426042200402D-02,  -1.171442818237512D-02, &
&    1.313255790973573D-03,  -8.371099516452187D-05,   2.332327334871007D-06, &
&    2.281922905825749D+00,  -2.850457119940149D+00,   1.838606956259343D+00, &
&   -7.276786094946810D-01,   1.890904180882148D-01,  -3.296386797751393D-02, &
&    3.813581584120174D-03,  -2.792448147132996D-04,   1.157321481202174D-05, &
&   -2.018677537492912D-07,   2.522613418072983D+00,  -3.347167236826221D+00, &
&    2.294223284946532D+00,  -9.714859435212977D-01,   2.729674006790809D-01, &
&   -5.220281515163948D-02,   6.755712070554076D-03,  -5.685059381042721D-04, &
&    2.816402860507764D-05,  -6.248236956132250D-07,   2.392463694043553D+00, &
&   -3.093955395832736D+00,   2.075255779445555D+00,  -8.610196537332251D-01, &
&    2.371386057955971D-01,  -4.445496071742221D-02,   5.638648897267226D-03, &
&   -4.649614499635758D-04,   2.256475450295823D-05,  -4.902398768025890D-07, &
&    2.145937647444013D+00,  -2.631801405379542D+00,   1.690150488548531D+00, &
&   -6.738047236875862D-01,   1.786235824007151D-01,  -3.226072843616258D-02, &
&    3.944304415481221D-03,  -3.136001516209016D-04,   1.467622326425910D-05, &
&   -3.074945512748113D-07,   1.902545939559848D+00,  -2.193377544986805D+00, &
&    1.339117421330826D+00,  -5.098338957524816D-01,   1.293802434109594D-01, &
&   -2.240054608830456D-02,   2.627931671017149D-03,  -2.006114668182249D-04, &
&    9.018321589031551D-06,  -1.815611881036022D-07,   1.698760592612476D+00, &
&   -1.840334358360287D+00,   1.067257847340245D+00,  -3.877038263953083D-01, &
&    9.410586515745797D-02,  -1.560774149743725D-02,   1.755780460356919D-03, &
&   -1.286180118934224D-04,   5.551312083652253D-06,  -1.073484549117828D-07, &
&    1.535758007637357D+00,  -1.568404091098504D+00,   8.656166142754264D-01, &
&   -3.004751789766771D-01,   6.984565529807552D-02,  -1.110911649012001D-02, &
&    1.199601203814210D-03,  -8.440957173131123D-05,   3.501314999817940D-06, &
&   -6.509526124657302D-08,   1.405138044663853D+00,  -1.358296585353194D+00, &
&    7.153960262808818D-01,  -2.378177606181386D-01,   5.304337147871756D-02, &
&   -8.105039014296582D-03,   8.415040453986343D-04,  -5.696583648760865D-05, &
&    2.274326714246108D-06,  -4.071192810443949D-08,   1.298496606911207D+00, &
&   -1.192690990251207D+00,   6.010878440391463D-01,  -1.917886898653514D-01, &
&    4.112718443148547D-02,  -6.048263688424006D-03,   6.048135693374317D-04, &
&   -3.945427884611453D-05,   1.518501982745022D-06,  -2.621184913044049D-08, &
&    1.209564029814801D+00,  -1.059201198185802D+00,   5.120271247188780D-01, &
&   -1.571250162707902D-01,   3.245335520362024D-02,  -4.601194585987844D-03, &
&    4.438566117780785D-04,  -2.794420778070469D-05,   1.038330578930851D-06, &
&   -1.730825714974198D-08,   1.133995188487756D+00,  -9.494378342324285D-01, &
&    4.411639021663465D-01,  -1.304360028035252D-01,   2.599101334929775D-02, &
&   -3.557947077695818D-03,   3.315708294828669D-04,  -2.017444318953867D-05, &
&    7.246851891792878D-07,  -1.168071763897105D-08,   1.068779872142251D+00, &
&   -8.576786561036843D-01,   3.837791719517366D-01,  -1.095002112816680D-01, &
&    2.108050681125562D-02,  -2.790053013426468D-03,   2.515110470422076D-04, &
&   -1.480818853741970D-05,   5.148519856610921D-07,  -8.033817659372651D-09, &
&    1.011779587362307D+00,  -7.799124505569585D-01,   3.366217605117409D-01, &
&   -9.281800575932100D-02,   1.728649134246527D-02,  -2.214770607895047D-03, &
&    1.933543912875351D-04,  -1.102846704240979D-05,   3.715459418360281D-07, &
&   -5.618833997239144D-09,   9.613901115608500D-01,  -7.131910271662479D-01, &
&    2.973541593751507D-01,  -7.933625347946173D-02,   1.431072463057332D-02, &
&   -1.776859225880105D-03,   1.503900283396192D-04,  -8.318456258752728D-06, &
&    2.718272775078469D-07,  -3.987941497667665D-09,   9.165526309123725D-01, &
&   -6.555196224552808D-01,   2.643839356220485D-01,  -6.834051008923832D-02, &
&    1.195314538712647D-02,  -1.439850416169431D-03,   1.182719583461593D-04, &
&   -6.350585515470334D-06,   2.014901001305661D-07,  -2.870524681794359D-09, &
&    8.762449238894968D-01,  -6.051156577245481D-01,   2.363694886681162D-01, &
&   -5.925729795070593D-02,   1.005977937101519D-02,  -1.176725686645202D-03, &
&    9.389268779118348D-05,  -4.898415677199682D-06,   1.510293060240204D-07, &
&   -2.091176767916272D-09,   8.397973374736750D-01,  -5.607717546712178D-01, &
&    2.123900746196760D-01,  -5.169276657599643D-02,   8.525650051352768D-03, &
&   -9.692951718191823D-04,   7.519387477589534D-05,  -3.814758613123273D-06, &
&    1.143932740559919D-07,  -1.540667804508020D-09,   8.066453081016311D-01, &
&   -5.214997953854644D-01,   1.917127957710592D-01,  -4.534177013409804D-02, &
&    7.271568748050543D-03,  -8.041981986548870D-04,   6.070341686089757D-05, &
&   -2.997121855987458D-06,   8.747939440781839D-08,  -1.146908820242620D-09, &
&    7.763328700748650D-01,  -4.865131195188174D-01,   1.737645584369129D-01, &
&   -3.997049326403647D-02,   6.238170884938409D-03,  -6.716462489699137D-04, &
&    4.936811186568123D-05,  -2.373940361615341D-06,   6.749309265777579D-08, &
&   -8.620120298597762D-10,   7.484872168611639D-01,  -4.551777208971689D-01, &
&    1.580916582586361D-01,  -3.539751489571131D-02,   5.380378405930295D-03, &
&   -5.643726412581538D-04,   4.042412621397579D-05,  -1.894535342555393D-06, &
&    5.250281623647732D-08,  -6.536814794056658D-10,   7.228002002617220D-01, &
&   -4.269770225009569D-01,   1.443309279425442D-01,  -3.148047647059620D-02, &
&    4.663564236690544D-03,  -4.769181755117127D-04,   3.331061449736308D-05, &
&   -1.522556290172317D-06,   4.115564531122071D-08,  -4.998334049864280D-10, &
&    6.990146269097030D-01,  -4.014859925399112D-01,   1.321887743925377D-01, &
&   -2.810653190187287D-02,   4.060849619196801D-03,  -4.051368747339830D-04, &
&    2.761111481851752D-05,  -1.231622827617666D-06,   3.249232324267357D-08, &
&   -3.851739792444945D-10,   6.769139442201537D-01,  -3.783518223820482D-01, &
&    1.214257096016929D-01,  -2.518540056502214D-02,   3.551170329362119D-03, &
&   -3.458487212194843D-04,   2.301316857550481D-05,  -1.002382719404352D-06, &
&    2.582501734241042D-08,  -2.989865561703131D-10,   6.563143660892078D-01, &
&   -3.572792987091244D-01,   1.118447759141316D-01,  -2.264424363481021D-02, &
&    3.117873830445191D-03,  -2.965922633845162D-04,   1.928010426958492D-05, &
&   -8.204973417272264D-07,   2.065535195607614D-08,  -2.336795023974850D-10, &
&    6.370587839081798D-01,  -3.380195885781189D-01,   1.032827760925571D-01, &
&   -2.042383776876847D-02,   2.747688779923439D-03,  -2.554461321436570D-04, &
&    1.623106686867063D-05,  -6.752440992628909D-07,   1.661870682317888D-08, &
&   -1.838200346615760D-10,   6.190120033975020D-01,  -3.203615420730855D-01, &
&    9.560355253336551D-02,  -1.847568500303102D-02,   2.429959362224784D-03, &
&   -2.208987753263688D-04,   1.372671585753052D-05,  -5.585351186952662D-07, &
&    1.344589462985757D-08,  -1.454832089595699D-10,   6.020569792050626D-01, &
&   -3.041248766244393D-01,   8.869278306835553D-02,  -1.675980697559348D-02, &
&    2.156070393944811D-03,  -1.917523666799123D-04,   1.165887107556909D-05, &
&   -4.642205187486950D-07,   1.093649208022802D-08,  -1.158080645769749D-10, &
&    5.860918095088875D-01,  -2.891547845930740D-01,   8.245391253436112D-02, &
&   -1.524304526077679D-02,   1.919011474773674D-03,  -1.670512229576893D-04, &
&    9.942940784081151D-06,  -3.875885889643397D-07,   8.940082061173694D-09, &
&   -9.269166157280035D-11,   5.710273159428749D-01,  -2.753176292173221D-01, &
&    7.680494417517188D-02,  -1.389774006501874D-02,   1.713043581841007D-03, &
&   -1.460279888984735D-04,   8.512332244579305D-06,  -3.250037532305761D-07, &
&    7.342928218127729D-09,  -7.457599172953888D-11/
!
      val1= abs(val)
      if(val1 < ten) then
        igrid= aint(val1*griddaw)
        val2= val1*val1
        val3= val1*val1*val1
        val4= val2*val2
        val5= val2*val3
        val6= val3*val3
        val7= val3*val4
        val8= val4*val4
        val9= val4*val4*val1
        dawerf= dawerfgrid(0,igrid)     +dawerfgrid(1,igrid)*val1+dawerfgrid(2,igrid)*val2 &
&              +dawerfgrid(3,igrid)*val3+dawerfgrid(4,igrid)*val4+dawerfgrid(5,igrid)*val5 &
&              +dawerfgrid(6,igrid)*val6+dawerfgrid(7,igrid)*val7+dawerfgrid(8,igrid)*val8 &
&              +dawerfgrid(9,igrid)*val9
      else
        val2inv= one/(val1*val1)
        last= 9
        factor=(last+half)*val2inv
        total= factor
        do i= last,1,-1
          factor= factor-val2inv
          total = factor*(one+total)
        enddo
        dawerf= val1*factor*(one+total)
      endif
!
      return
end


!---------------------------------
  function type1da0(norder,zeta)
!---------------------------------
!
! Calculate type 1 radial integral for the case of DA=zero
!
      implicit none
      integer,intent(in) :: norder
      integer :: nhalf
      real(8),parameter :: sqrtpi=1.772453850905516D+00
      real(8),intent(in) :: zeta
      real(8) :: type1da0
      real(8) :: gammao(0:8), gammae(0:8), tmp
      data gammao/0.5D+00, 0.5D+00, 1.0D+00, 3.0D+00, 12.0D+00, 60.0D+00, &
&                  360.0D+00, 2520.0D+00, 20160.0D+00/
      data gammae/0.5D+00, 0.25D+00, 0.375D+00, 0.9375D+00, 3.28125D+00, 14.765625D+00, &
&                 81.2109375D+00, 527.87109375D+00, 3959.033203125D+00/
!
      nhalf= norder/2
!
! N : odd
!
      if(mod(norder,2) /= 0) then
        tmp= gammao(nhalf)
!
! N : even
!
      else
        tmp= gammae(nhalf)*sqrtpi*sqrt(zeta)
      endif
      type1da0= tmp/(zeta**(nhalf+1))
!
      return
end


!-------------------------------------------
  subroutine calczspher(zzlm,unitvec,lmax)
!-------------------------------------------
!
! Calculate the real spherical harmonics given the value of the three 
! cartesian components of a unit vector.
!
      use modecp, only : zlm, lmf, lmx, lmy, lmz
      implicit none
      integer,intent(in) :: lmax
      integer :: lambda, mu, imn, imx, ii, id
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: unitvec(3)
      real(8),intent(out) :: zzlm(121)
      real(8) :: sumz, tmp
!
      if(lmax <= 10) then
        do lambda= 0,lmax
          id= lambda*(lambda+1)+lambda+1
          do mu=-lambda,lambda
            imn= lmf(id)
            imx= lmf(id+1)-1
            sumz= zero
            do ii=imn,imx
              tmp = zlm(ii)
              if(lmx(ii) > 0) tmp= tmp*(unitvec(1)**lmx(ii))
              if(lmy(ii) > 0) tmp= tmp*(unitvec(2)**lmy(ii))
              if(lmz(ii) > 0) tmp= tmp*(unitvec(3)**lmz(ii))
              sumz= sumz+tmp
            enddo
            zzlm(id)= sumz
            id= id-1
          enddo
        enddo
      else
        write(*,'(" Error! Lmax is too big in calczspher,",i3,".")') lmax
        call exit
      endif
!
      return
end


!--------------------------------------------
  function type1radint(nn,ll,alpha,xp0,xp1)
!--------------------------------------------
!
! Calculate core terms for integrals with exp and modified spherical Bessel function
!
      implicit none
      integer,intent(in) :: nn, ll
      integer :: np1
      real(8),parameter :: athresh=1.0D-01, sqrtpi=1.772453850905516D+00, half=0.5D+00
      real(8),parameter :: one=1.0D+00, p15=1.5D+00, two=2.0D+00, p25=2.5D+00, three=3.0D+00
      real(8),parameter :: five=5.0D+00, six=6.0D+00, p75=7.5D+00, fifteen=15.0D+00
      real(8),intent(in) :: alpha, xp0, xp1
      real(8) :: type1radint, radtmp, dawrad, dawson, type1radintps, alpha2, xalpha
      real(8) :: errrad, alpha3, alpha4, erf
!
      if(ll == 0)then
        if(alpha >= athresh) then
          np1= nn+1
          if(np1 == 2) then
            radtmp= sqrtpi*xp1*erf(alpha)/alpha
          else
            dawrad=dawson(alpha)
            radtmp= sqrtpi*xp1*dawrad/alpha
          endif
        else
          radtmp= type1radintps(nn,0,alpha,xp0)
        endif
      elseif(ll == 1)then
        if(alpha >= athresh) then
          np1= nn+1
          alpha2= alpha*alpha
          xalpha= alpha2+alpha2
          dawrad= dawson(alpha)
          errrad= sqrtpi*xp1*erf(alpha)
          if(np1 < 2) then
            radtmp=(half*errrad-xp0*alpha)/alpha2
          elseif(np1 == 2) then
            radtmp= sqrtpi*xp1*(alpha-dawrad)/alpha2
          else
            radtmp=(two*xp0*alpha+(xalpha-one)*errrad)/alpha2
          endif
        else
          radtmp= type1radintps(nn,1,alpha,xp0)
        endif
      elseif(ll == 2)then
        if(alpha >= athresh) then
          np1= nn+1
          alpha2= alpha*alpha
          alpha3= alpha*alpha2
          xalpha= alpha2+alpha2
          dawrad= dawson(alpha)
          errrad= sqrtpi*xp1*erf(alpha)
          if(np1 == 1)then
            radtmp= half*sqrtpi*xp1*(p15*alpha-(alpha2+p15)*dawrad)/alpha3
          elseif(np1 == 2) then
            radtmp=(three*xp0*alpha+(alpha2-p15)*errrad)/alpha3
          elseif(np1 == 3) then
            radtmp= sqrtpi*xp1*(alpha*(xalpha-three)+three*dawrad)/alpha3
          elseif(np1 == 4) then
            radtmp=(two*xp0*alpha*(xalpha-three)+(xalpha*(xalpha-two)+three)*errrad)/alpha3
          endif
        else
          radtmp= type1radintps(nn,2,alpha,xp0)
        endif
      elseif(ll == 3) then
        if(alpha >= athresh) then
          np1= nn+1
          alpha2= alpha*alpha
          alpha3= alpha*alpha2
          alpha4= alpha*alpha3
          xalpha= alpha2+alpha2
          dawrad= dawson(alpha)
          errrad= sqrtpi*xp1*erf(alpha)
          if(np1 == 1) then
            radtmp=(two*xp0*alpha*(xalpha+p75)/six+half*(alpha2-p25)*errrad)/alpha4
          elseif(np1 == 2) then
            radtmp= half*sqrtpi*xp1*(alpha*(xalpha-p75)+p15*(xalpha+five)*dawrad)/alpha4
          elseif(np1 == 3) then
            radtmp=(two*xp0*alpha*(alpha2-p75)+(xalpha*(alpha2-three)+p75)*errrad)/alpha4
          elseif(np1 == 4) then
            radtmp= sqrtpi*xp1*(alpha*(xalpha*(xalpha-five)+fifteen)-fifteen*dawrad)/alpha4
          elseif(np1 == 5) then
            radtmp=(two*xp0*alpha*(xalpha*(xalpha-one-three)+fifteen)+ &
&                  (xalpha*(xalpha*(xalpha-three)+six+three)-fifteen)*errrad)/alpha4
          endif
        else
          radtmp= type1radintps(nn,3,alpha,xp0)
        endif
      endif
      type1radint= radtmp
!
      return
end


!------------------------------------------
  function type1radintps(nn,ll,alpha,xp0)
!------------------------------------------
!
! Calculate core terms for integrals with exp and modified spherical Bessel function
! by a power series
!
      implicit none
      integer,intent(in) :: nn, ll
      integer,parameter :: maxii=10
      integer :: nl, lambda, l2, ii
      real(8),parameter :: one=1.0D+00, two=2.0D+00, sqrtpi=1.772453850905516D+00
      real(8) :: fctrle(10),fctrlo(10)
      real(8),intent(in) :: alpha, xp0
      real(8) :: rgamma, fac, alpha2, tmp, type1radintps
      data fctrle/1.0D+00,1.0D+00,2.0D+00,6.0D+00,24.0D+00,120.0D+00, &
&                 720.0D+00,5040.0D+00,40320.0D+00,362880.0D+00/
      data fctrlo/1.0D+00,1.0D+00,3.0D+00,15.0D+00,105.0D+00,945.0D+00, &
&                 10395.0D+00,135135.0D+00,2027025.0D+00,34459425.0D+00/
!
      nl= nn+ll
      lambda= nl/2
!
! Odd integral
!
      if(mod(nl,2) /= 0) then
        rgamma  =  fctrle(lambda+1)/fctrlo(ll+2)
        fac= two**nl
      else
!
! Even integral
!
        rgamma  = fctrlo(lambda+1)/fctrlo(ll+2)
        fac=(two**lambda)*sqrtpi
      endif
      l2= ll+ll+1
      alpha2 = alpha*alpha
      tmp= rgamma
      do ii=1,maxii
        rgamma =(rgamma*alpha2*(ii+ii+nl-1))/(ii*(ii+ii+l2))
        tmp= tmp+rgamma
      enddo
      type1radintps = fac*xp0*(alpha**ll)*tmp
!
      return
end


!--------------------------------------------------------------------------------------------
  subroutine calctype1rad(rad1int,besselarginv,alpha,sqrtinvzeta,xp0,xp1,nangecp,nangbasis)
!--------------------------------------------------------------------------------------------
!
! Calculate terms for type1 integrals with exp and modified spherical Bessel function
!
      integer,intent(in) :: nangecp, nangbasis
      integer :: nangtotal, nangeven, nangodd, nk1, nk2, nn, ll
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00
      real(8),intent(in) :: besselarginv, alpha, sqrtinvzeta, xp0, xp1
      real(8),intent(out) :: rad1int(78)
      real(8) :: sizeta2, a1, a2, radtmp(19,2), type1radint, t1, t2, tmp
!
      nangtotal= nangecp+nangbasis
      if(mod(nangtotal,2) /= 0) then
        nangeven= nangtotal+1
        nangodd = nangtotal
      else
        nangeven= nangtotal
        nangodd = nangtotal+1
      endif
      if(nangtotal == 0) nangodd=-1
!
      sizeta2= sqrtinvzeta*half
      a1= sqrtinvzeta*alpha
      a2= sqrtinvzeta*sizeta2
!
! Generate table of elements of even lattice
!
      if(nangeven >=0) then
        radtmp(1,1)= type1radint(0,0,alpha,xp0,xp1)*sizeta2
        if(nangeven >= 1) then
          radtmp(2,2)= type1radint(1,1,alpha,xp0,xp1)*sizeta2*sizeta2
          t1= a2
          do nn= 2,nangeven,2
            t2= t1-a2
            radtmp(nn+1,1)= a1*radtmp(nn  ,2)+t1*radtmp(nn-1,1)
            radtmp(nn+2,2)= a1*radtmp(nn+1,1)+t2*radtmp(nn  ,2)
            t1= t1+a2*two
          enddo
        endif
      endif
!
! Generate table of elelments of odd lattice
!
      if(nangodd >= 0) then
        radtmp(1,2)= type1radint(0,1,alpha,xp0,xp1)*sizeta2
        if(nangodd >= 1) then
          radtmp(2,1)= type1radint(1,0,alpha,xp0,xp1)*sizeta2*sizeta2
          t1= a2*two
          do nn= 2,nangodd,2
            t2= t1-a2*three
            radtmp(nn+1,2)= a1*radtmp(nn  ,1)+t2*radtmp(nn-1,2)
            radtmp(nn+2,1)= a1*radtmp(nn+1,2)+t1*radtmp(nn  ,1)
            t1= t1+a2*two
          enddo
        endif
      endif
!
! Retrieve required integrals from tables
!
      rad1int(1)= radtmp(nangecp+1,1)
      rad1int(2)= radtmp(nangecp+2,1)
      rad1int(3)= radtmp(nangecp+2,2)
      nk1= 3
      nk2=1
      do nn= 3,nangbasis+1
        rad1int(nk1+1)= radtmp(nangecp+nn,1)
        rad1int(nk1+2)= radtmp(nangecp+nn,2)
        tmp= besselarginv*three
        do ll= 3,nn
          rad1int(nk1+ll)= rad1int(nk1+ll-2)-tmp*rad1int(nk2+ll-1)
          tmp= tmp+besselarginv*two
        enddo
        nk2= nk1
        nk1= nk1+nn
      enddo
!
      return
end


!---------------------------------------------------
  subroutine rad2recur(rad2core,nmin,nmax,lmax,xx)
!---------------------------------------------------
!
! Generate table for type2 radial integral
!
      implicit none
      integer,intent(in) :: nmin, nmax, lmax
      integer :: nn
      real(8),parameter :: two=2.0D+00
      real(8),intent(in) :: xx
      real(8),intent(inout) :: rad2core(30,7)
      real(8) :: xx2, t01, t02
!
      xx2= xx*two
      t01=(nmin-two)*two
!
      select case(lmax)
        case (0)
          do nn= nmin,nmax,2
            t02= t01+two
            rad2core(nn+1,1)= t02*rad2core(nn-1,1)+xx2*rad2core(nn  ,2)
            rad2core(nn+2,2)= t01*rad2core(nn  ,2)+xx2*rad2core(nn+1,1)
            t01= t02+two
          enddo
        case (1)
          do nn= nmin,nmax,2
            t02= t01+two
            rad2core(nn+1,1)= t02*rad2core(nn-1,1)+xx2*rad2core(nn  ,2)
            rad2core(nn+2,2)= t01*rad2core(nn  ,2)+xx2*rad2core(nn+1,1)
            rad2core(nn+3,3)= t01*rad2core(nn+1,3)+xx2*rad2core(nn+2,2)
            t01= t02+two
          enddo
        case (2)
          do nn= nmin,nmax,2
            t02= t01+two
            rad2core(nn+1,1)= t02*rad2core(nn-1,1)+xx2*rad2core(nn  ,2)
            rad2core(nn+2,2)= t01*rad2core(nn  ,2)+xx2*rad2core(nn+1,1)
            rad2core(nn+3,3)= t01*rad2core(nn+1,3)+xx2*rad2core(nn+2,2)
            rad2core(nn+4,4)= t01*rad2core(nn+2,4)+xx2*rad2core(nn+3,3)
            t01= t02+two
          enddo
        case (3)
          do nn= nmin,nmax,2
            t02= t01+two
            rad2core(nn+1,1)= t02*rad2core(nn-1,1)+xx2*rad2core(nn  ,2)
            rad2core(nn+2,2)= t01*rad2core(nn  ,2)+xx2*rad2core(nn+1,1)
            rad2core(nn+3,3)= t01*rad2core(nn+1,3)+xx2*rad2core(nn+2,2)
            rad2core(nn+4,4)= t01*rad2core(nn+2,4)+xx2*rad2core(nn+3,3)
            rad2core(nn+5,5)= t01*rad2core(nn+3,5)+xx2*rad2core(nn+4,4)
            t01= t02+two
          enddo
        case (4)
          do nn= nmin,nmax,2
            t02= t01+two
            rad2core(nn+1,1)= t02*rad2core(nn-1,1)+xx2*rad2core(nn  ,2)
            rad2core(nn+2,2)= t01*rad2core(nn  ,2)+xx2*rad2core(nn+1,1)
            rad2core(nn+3,3)= t01*rad2core(nn+1,3)+xx2*rad2core(nn+2,2)
            rad2core(nn+4,4)= t01*rad2core(nn+2,4)+xx2*rad2core(nn+3,3)
            rad2core(nn+5,5)= t01*rad2core(nn+3,5)+xx2*rad2core(nn+4,4)
            rad2core(nn+6,6)= t01*rad2core(nn+4,6)+xx2*rad2core(nn+5,5)
            t01= t02+two
          enddo
        case (5)
          do nn= nmin,nmax,2
            t02= t01+two
            rad2core(nn+1,1)= t02*rad2core(nn-1,1)+xx2*rad2core(nn  ,2)
            rad2core(nn+2,2)= t01*rad2core(nn  ,2)+xx2*rad2core(nn+1,1)
            rad2core(nn+3,3)= t01*rad2core(nn+1,3)+xx2*rad2core(nn+2,2)
            rad2core(nn+4,4)= t01*rad2core(nn+2,4)+xx2*rad2core(nn+3,3)
            rad2core(nn+5,5)= t01*rad2core(nn+3,5)+xx2*rad2core(nn+4,4)
            rad2core(nn+6,6)= t01*rad2core(nn+4,6)+xx2*rad2core(nn+5,5)
            rad2core(nn+7,7)= t01*rad2core(nn+5,7)+xx2*rad2core(nn+6,6)
            t01= t02+two
          enddo
        case default
          write(*,*)"lmax=",lmax
          call exit
      end select
!
      return
end


!----------------------------------------------------------------------
  subroutine rad2table(rad2core,alpha,beta,xp1p,xp1m,xp0,lemax,lomax)
!----------------------------------------------------------------------
!
! Prepare table for type2 radial integral
!
      implicit none
      integer,intent(in) :: lemax, lomax
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: alpha, beta, xp1p, xp1m, xp0
      real(8),intent(out) :: rad2core(30,7)
      real(8) :: alf, xp1, type1radint, type1radintps
!
      if((alpha-beta) <= zero) then
        alf= beta
        xp1 = xp1m
      else
        alf= alpha
        xp1 = xp1p
      endif
!
! Generate even lattice
!
      rad2core(1,1)= type1radint(0,0,alf,xp0,xp1)
      rad2core(2,2)= type1radint(1,1,alf,xp0,xp1)
      if(lemax >= 2) rad2core(3,3)= type1radint(2,2,alf,xp0,xp1)
      if(lemax >= 3) rad2core(4,4)= type1radint(3,3,alf,xp0,xp1)
      if(lemax >= 4) rad2core(5,5)= type1radintps(4,4,alf,xp0)
      if(lemax >= 5) rad2core(6,6)= type1radintps(5,5,alf,xp0)
      if(lemax >= 6) rad2core(7,7)= type1radintps(6,6,alf,xp0)
      call rad2recur(rad2core,2,22,lemax,alf)
!
! Generate odd lattice
!
      rad2core(1,2)= type1radint(0,1,alf,xp0,xp1)
      rad2core(2,1)= type1radint(1,0,alf,xp0,xp1)
      rad2core(3,2)= type1radint(2,1,alf,xp0,xp1)
      if(lomax >= 2) rad2core(4,3)= type1radint(3,2,alf,xp0,xp1)
      if(lomax >= 3) rad2core(5,4)= type1radint(4,3,alf,xp0,xp1)
      if(lomax >= 4) rad2core(6,5)= type1radintps(4,5,alf,xp0)
      if(lomax >= 5) rad2core(7,6)= type1radintps(5,6,alf,xp0)
      if(lomax >= 6) rad2core(8,7)= type1radintps(6,7,alf,xp0)
      call rad2recur(rad2core,3,21,lomax,alf)
!
      return
end


!--------------------------------------------------------------------------
  function type2radintps(rad2core,nn,lalpha,lbeta,alpha,beta,sqrtinvzeta)
!--------------------------------------------------------------------------
!
! Calculate terms for type2 radial integrals
!
      implicit none
      integer,parameter :: maxii=10
      integer,intent(in) :: nn, lalpha, lbeta
      integer :: ll, l1, l2, l3, ii
      real(8),parameter :: zero=0.0D+00, half=0.5D+00
      real(8),intent(in) :: rad2core(30,7), alpha, beta, sqrtinvzeta
      real(8) :: type2radintps, xx, tmp, sumrad2
      real(8) :: fctrlo(7)=(/1.0D+00,3.0D+00,15.0D+00,105.0D+00,945.0D+00,1.0395D+04,1.35135D+05/)
!
      if((alpha-beta) > zero) then
        ll= lbeta+1
        xx= beta
      else
        ll= lalpha+1
        xx= alpha
      endif
      tmp = xx**(ll-1)/fctrlo(ll)
      xx = xx*xx
      l1= lalpha+lbeta+2-ll
      l2= ll+nn
      l3= ll+ll-1
!
      sumrad2 = tmp*rad2core(l2,l1)
      do ii=1,maxii
        tmp= tmp*xx/(2*ii*(ii+ii+l3))
        sumrad2= sumrad2+tmp*rad2core(ii+ii+l2,l1)
      enddo
      type2radintps= sumrad2*(half*sqrtinvzeta)**(nn+1)
!
      return
end


!-----------------------------------------------------------------------
  subroutine rad2hyperbolic(rho,sga,sgb,tau,sqrtinvzeta,xp0,xka,xkb, &
&                           gamma1,gamma2,a1,a2,zeta,nmax)
!-----------------------------------------------------------------------
!
! Calculate rho, sigma, sigma bar, tau for type2 radial integral
! See J.Chem.Phys. 111, 8778-8784 (1999).
!
      implicit none
      integer,intent(in) :: nmax
      integer :: ii
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00
      real(8),parameter :: sqrtpi=1.772453850905516D+00
      real(8),parameter :: sqrtpi2=3.544907701811032D+00
      real(8),intent(in) :: sqrtinvzeta, xp0, xka, xkb, gamma1, gamma2, a1, a2, zeta
      real(8),intent(out) :: rho(15), sga(15), sgb(15), tau(15)
      real(8) :: factor, factbc(13), b1(15), b2(15), c1(15), c2(15), dawerf, dawson, erf
      real(8) :: zeta2, xkabp, xkabm, xp0p, xp0m, facp, denom, facm, facb, facc, facbc
      data factbc/1.0D+00,4.0D+00,1.8D+01,9.6D+01,6.0D+02,4.32D+03, &
&                 3.528D+04,3.2256D+05,3.26592D+06,3.6288D+07, &
&                 4.390848D+08,5.7480192D+09,8.09512704D+10/
!
      factor= sqrtpi*sqrtinvzeta
      b1(1) = factor*gamma1
      c1(1) = factor*gamma1*erf(a1)
      b2(1) = factor*gamma2
      c2(1) = factor*gamma2*erf(a2)
      b1(2) = sqrtpi2*gamma1*dawerf(a1)
      c1(2) = sqrtpi2*gamma1*dawson(a1)
      b2(2) = sqrtpi2*gamma2*dawerf(a2)
      c2(2) = sqrtpi2*gamma2*dawson(a2)
!
      zeta2= zeta+zeta
      xkabp= xka+xkb
      xkabm= xka-xkb
      xp0p= xp0
      xp0m= xp0
      facp= zero
      denom= one
      do ii= 3,nmax
        xp0p= xp0p*xkabp
        xp0m= xp0m*xkabm
        facm= half-facp
        facb= facp/factbc(ii-2)
        facc= facm/factbc(ii-2)
        facbc= one/denom
        b1(ii)= xp0p*facb-(zeta2*b1(ii-2)-xkabp*c1(ii-1))*facbc
        c1(ii)= xp0p*facc-(zeta2*c1(ii-2)-xkabp*b1(ii-1))*facbc
        b2(ii)= xp0m*facb-(zeta2*b2(ii-2)-xkabm*c2(ii-1))*facbc
        c2(ii)= xp0m*facc-(zeta2*c2(ii-2)-xkabm*b2(ii-1))*facbc
        facp= facm
        denom= denom+one
      enddo
!
      do ii= 1,nmax
        rho(ii)= b1(ii)-b2(ii)
        sga(ii)= c1(ii)+c2(ii)
        sgb(ii)= c1(ii)-c2(ii)
        tau(ii)= b1(ii)+b2(ii)
      enddo
!
      return
      end


!-----------------------------------------------------------------------------------------
  subroutine calctype2rad(rad2int,alfbessel,betbessel,alpha,beta,sqrtinvzeta, &
&                         xp1p,xp1m,xp0,xka,xkb,a1,a2,zeta,nangecp,nangbasis,lmax,lemax)
!-----------------------------------------------------------------------------------------
!
! Calculate type2 radial integral
!
      implicit none
      integer,intent(in) :: nangecp, nangbasis, lmax, lemax
      integer :: nangtotal, nemax, nomax, lomax, ii, jj, kk, nn
      real(8),parameter :: zero=0.0D+00, p01=0.1D+00, quarter=0.25D+00, half=0.5D+00
      real(8),parameter :: one=1.0D+00, two=2.0D+00, three=3.0D+00, five=5.0D+00
      real(8),parameter :: six=6.0D+00, ten=1.0D+01, p35=3.5D+00, abthresh=1.0D-01
      real(8),intent(in) :: alfbessel, betbessel, alpha, beta, sqrtinvzeta
      real(8),intent(in) :: xp1p, xp1m, xp0, xka, xkb, a1, a2, zeta
      real(8),intent(out) :: rad2int(11,11,11)
      real(8) :: rho(15), sga(15), sgb(15), tau(15), radtmp(19,2,2), rad2core(30,7)
      real(8) :: xkakb, denomab, denomab2, gamma1, gamma2, zzalpha, zzbeta, zetainv
      real(8) :: t1, t2, t3, xka2, xkb2, xka4, xkb4, xka6, xkb6, xkakb2
      real(8) :: t00, t01, t02, t03, t04, t05, t06, qa2, qb2, qa4, qb4, qa6, qb6, qab
      real(8) :: ra0, ra1, ra2, ra3, ra4, type2radintps
!
      nangtotal= nangecp+nangbasis-1
      nemax= nangtotal
      nomax= nangtotal
      if(nangtotal == 0) nomax=-1
      lomax= lemax
!
      if(alpha*beta > abthresh) then
        xkakb = xka*xkb
        denomab= one/xkakb
        denomab2= denomab*denomab
        if(lemax > 6) then
          write(*,'(" Error! Lemax is too big, lemax=",i3)')lemax
          call exit
        endif
        do ii= 1,15
          rho(ii)= zero
          sga(ii)= zero
          sgb(ii)= zero
          tau(ii)= zero
        enddo
        gamma1= xp1p*quarter
        gamma2= xp1m*quarter
        call rad2hyperbolic(rho,sga,sgb,tau,sqrtinvzeta,xp0,xka,xkb, &
&                           gamma1,gamma2,a1,a2,zeta,2*lemax+3)
        radtmp(1,1,1)= rho(3)*denomab
        radtmp(2,1,1)= rho(2)*denomab
        radtmp(1,2,1)=(sgb(3)-rho(4)*denomab*xkb)*denomab
        radtmp(2,2,1)=(sgb(2)-rho(3)*denomab*xkb)*denomab
        radtmp(1,1,2)=(sga(3)-rho(4)*denomab*xka)*denomab
        radtmp(2,1,2)=(sga(2)-rho(3)*denomab*xka)*denomab
        radtmp(1,2,2)=(rho(5)-xka*sgb(4)-xkb*sga(4)+xkakb*tau(3))*denomab2
        radtmp(2,2,2)=(rho(4)-xka*sgb(3)-xkb*sga(3)+xkakb*tau(2))*denomab2
      else
        call rad2table(rad2core,alpha,beta,xp1p,xp1m,xp0,lemax,lomax)
        radtmp(1,1,1)=type2radintps(rad2core,0,0,0,alpha,beta,sqrtinvzeta)
        radtmp(2,1,1)=type2radintps(rad2core,1,0,0,alpha,beta,sqrtinvzeta)
        radtmp(1,2,1)=type2radintps(rad2core,0,1,0,alpha,beta,sqrtinvzeta)
        radtmp(2,2,1)=type2radintps(rad2core,1,1,0,alpha,beta,sqrtinvzeta)
        radtmp(1,1,2)=type2radintps(rad2core,0,0,1,alpha,beta,sqrtinvzeta)
        radtmp(2,1,2)=type2radintps(rad2core,1,0,1,alpha,beta,sqrtinvzeta)
        radtmp(1,2,2)=type2radintps(rad2core,0,1,1,alpha,beta,sqrtinvzeta)
        radtmp(2,2,2)=type2radintps(rad2core,1,1,1,alpha,beta,sqrtinvzeta)
      endif
      zzalpha= sqrtinvzeta*alpha
      zzbeta= sqrtinvzeta*beta
      zetainv= sqrtinvzeta*sqrtinvzeta*half
!
! Generate table of elements of even lattice
!
      t2= zero
      do nn= 2,nemax,2
        t1= t2+zetainv
        t3= t2-zetainv*three
        radtmp(nn+1,1,1)= zzalpha*radtmp(nn  ,2,1)+zzbeta*radtmp(nn  ,1,2)+t1*radtmp(nn-1,1,1)
        radtmp(nn+1,2,2)= zzalpha*radtmp(nn  ,1,2)+zzbeta*radtmp(nn  ,2,1)+t3*radtmp(nn-1,2,2)
        radtmp(nn+2,2,1)= zzalpha*radtmp(nn+1,1,1)+zzbeta*radtmp(nn+1,2,2)+t2*radtmp(nn  ,2,1)
        radtmp(nn+2,1,2)= zzalpha*radtmp(nn+1,2,2)+zzbeta*radtmp(nn+1,1,1)+t2*radtmp(nn  ,1,2)
        t2= t2+zetainv*two
      enddo
!
! Generate table of elements of odd lattice
!
      t2=-zetainv
      do nn= 2,nomax,2
        t1= t2-zetainv
        t3= t2+zetainv*three
        radtmp(nn+1,2,1)= zzalpha*radtmp(nn  ,1,1)+zzbeta*radtmp(nn  ,2,2)+t2*radtmp(nn-1,2,1)
        radtmp(nn+1,1,2)= zzalpha*radtmp(nn  ,2,2)+zzbeta*radtmp(nn  ,1,1)+t2*radtmp(nn-1,1,2)
        radtmp(nn+2,1,1)= zzalpha*radtmp(nn+1,2,1)+zzbeta*radtmp(nn+1,1,2)+t3*radtmp(nn  ,1,1)
        radtmp(nn+2,2,2)= zzalpha*radtmp(nn+1,1,2)+zzbeta*radtmp(nn+1,2,1)+t1*radtmp(nn  ,2,2)
        t2= t2+zetainv*two
      enddo
!
! Retrieve required type2 radial integrals from tables
!
      if((lemax >= 2).and.(lemax <= 6)) then
        if(alpha*beta > abthresh) then
          if(nangecp > 2) then
            write(*,'(" Error! Nangecp is too big,",i3,", lemax=",i3,".")')nangecp,lemax
            call exit
          endif
        endif
      elseif(lemax > 6) then
        call exit
      endif
!
      nn= nangecp
      do ii= 1,lemax+1
        if(ii <  3) then
          rad2int(ii,ii,1)= radtmp(nn+1,ii,ii)
          cycle
        elseif(alpha*beta <= abthresh) then
          rad2int(ii,ii,1)= type2radintps(rad2core,nn,ii-1,ii-1,alpha,beta,sqrtinvzeta)
          cycle
        endif
        if(ii == 3) then
          xka2 = xka*xka
          xkb2 = xkb*xkb
          xkakb2= xkakb*xkakb
          ra4= denomab2
          qa2= 3.0D+00*xka2
          qb2= 3.0D+00*xkb2
          t02= qa2+qb2
          t00= 9.0D+00
          ra0= t00*rho( 7-nn)+t02*rho( 5-nn)+xkakb2*rho( 3-nn)
          ra1= t00*sga( 6-nn)+qa2*sga( 4-nn)
          ra2= t00*sgb( 6-nn)+qb2*sgb( 4-nn)
          ra3= t00*tau( 5-nn)
        elseif(ii == 4) then
          qa2= 1.5D+01*xka2
          qb2= 1.5D+01*xkb2
          t02= qa2+qb2
          qab= six*xkakb2
          t00= 2.25D+02
          ra0= t00*rho( 9-nn)+(t02*rho( 7-nn)+rho( 5-nn)*qab)*six
          ra1= t00*sga( 8-nn)+(t02+five* qa2)*sga( 6-nn)+sga( 4-nn)*qab
          ra2= t00*sgb( 8-nn)+(t02+five* qb2)*sgb( 6-nn)+sgb( 4-nn)*qab
          ra3= t00*tau( 7-nn)+ t02*tau( 5-nn)+tau( 3-nn)*xkakb2
        elseif(ii == 5) then
          xka4 = xka2*xka2
          xkb4 = xkb2*xkb2
          qa2= 4.5D+01*xka2
          qb2= 4.5D+01*xkb2
          t02= qa2+qb2
          t03= 1.05D+03*(xka2+xkb2)
          qa4= 1.05D+02*xka4
          qb4= 1.05D+02*xkb4
          t04= qa4+qb4
          qab= ten*xkakb2
          t00= 1.1025D+04
          ra0 = t00*rho(11-nn)+1.05D+02*t02*rho( 9-nn)+t04*rho( 7-nn)+ &
&              (2.025D+03*rho( 7-nn)+t02*rho( 5-nn)+rho( 3-nn)*xkakb2)*xkakb2
          ra1 = t00*sga(10-nn)+(t03+3.675D+03*xka2)*sga( 8-nn)+ &
&               qa4*sga( 6-nn)+(4.5D+01*sga( 6-nn)+xka2*sga( 4-nn))*qab
          ra2 = t00*sgb(10-nn)+(t03+3.675D+03*xkb2)*sgb( 8-nn)+ &
&               qb4*sgb( 6-nn)+(4.5D+01*sgb( 6-nn)+xkb2*sgb( 4-nn))*qab
          ra3 = t00*tau( 9-nn)+t03*tau( 7-nn)+ten*tau( 5-nn)*qab
        elseif(ii == 6) then
          qa2= 4.2D+02*xka2
          qb2= 4.2D+02*xkb2
          t02= qa2+qb2
          t03= 9.9225D+04*(xka2+xkb2)
          qa4= 9.45D+02*xka4
          qb4= 9.45D+02*xkb4
          t04= qa4+qb4
          qab= 1.5D+01*xkakb2
          t00= 8.93025D+05
          ra0= t00*rho(13-nn)+4.0D+00*t03*rho(11-nn)+1.5D+01*t04*rho( 9-nn)+ &
&             (1.176D+04*rho( 9-nn)+t02*rho( 7-nn)+rho( 5-nn)*qab)*qab
          ra1= t00*sga(12-nn)+(t03+7.0875D+02*qa2)*sga(10-nn)+ &
&             (t04+1.4D+01*qa4)*sga( 8-nn)+(4.41D+04*sga( 8-nn)+ &
&             (t02+2.75D+00*qa2)*sga( 6-nn)+sga( 4-nn)*qab)*xkakb2
          ra2= t00*sgb(12-nn)+(t03+7.0875D+02*qb2)*sgb(10-nn)+ &
              (t04+1.4D+01*qb4)*sgb( 8-nn)+(4.41D+04*sgb( 8-nn)+ &
&             (t02+2.75D+00*qb2)*sgb( 6-nn)+sgb( 4-nn)*qab)*xkakb2
          ra3= t00*tau(11-nn)+t03*tau( 9-nn)+t04*tau( 7-nn)+ &
&              (1.1025D+04*tau( 7-nn)+2.5D-01*t02*tau( 5-nn)+tau( 3-nn)*xkakb2)*xkakb2
        elseif(ii == 7) then
          xka6 = xka4*xka2
          xkb6 = xkb4*xkb2
          qa2= p35*xka2
          qb2= p35*xkb2
          t02= qa2+qb2
          t03= 6.496875D+01*t02
          qa4= p35*xka4
          qb4= p35*xkb4
          t04= 3.75D-01*(qa4+qb4)
          t05= 7.7D+00*t04
          qa6= 4.8125D-02*xka6
          qb6= 4.8125D-02*xkb6
          t06= qa6+qb6
          qab= p01*xkakb2/six
          t00= 5.00259375D+02
          ra0= t00*rho(15-nn)+t03*rho(13-nn)+t05*rho(11-nn)+t06*rho( 9-nn)+ &
&             (6.2015625D+03*rho(11-nn)+7.875D+01*t02*rho( 9-nn)+t04*rho( 7-nn)+ &
              (7.35D+02*rho( 7-nn)+t02*rho( 5-nn)+rho( 3-nn)*qab)*qab)*qab
          ra1= t00*sga(14-nn)+(t03-1.66753125D+02*xkb2)*sga(12-nn)+ &
&             (t05-2.59875D+00*qb4)*sga(10-nn)+qa6*sga( 8-nn)+ &
&             (1.65375D+04*sga(10-nn)+2.1D+02*(t02-6.25D-01*qb2)*sga( 8-nn)+ &
&              qa4*sga( 6-nn)+(7.35D+02*sga( 6-nn)+Qa2*sga( 4-nn))*qab)*qab*p01
          ra2= t00*sgb(14-nn)+(t03-1.66753125D+02*xka2)*sgb(12-nn)+ &
&             (t05-2.59875D+00*qa4)*sgb(10-nn)+qb6*sgb( 8-nn)+ &
&             (1.65375D+04*sgb(10-nn)+2.1D+02*(t02-6.25D-01*qa2)*sgb( 8-nn)+ &
&              qb4*sgb( 6-nn)+(7.35D+02*sgb( 6-nn)+qb2*sgb( 4-nn))*qab)*qab*p01
          ra3= t00*tau(13-nn)+1.7325D+01*t02*tau(11-nn)+p01*t05 *tau( 9-nn)+ &
&             (2.1D+02*tau( 9-nn)+t02*tau( 7-nn)+p35*tau( 5-nn)*qab)*qab*2.1D+00
          ra4= ra4*2.16D+05
         endif
         ra4= ra4*denomab
         rad2int(ii,ii,1)=(ra0-ra1*xkb-ra2*xka+ra3*xkakb)*ra4
      enddo
!
      do ii= 2,nangbasis
        do jj= 1,2
          do kk= 1,2
            rad2int(kk,jj,ii)= radtmp(nangecp+ii,jj,kk)
          enddo
        enddo
        t01= alfbessel*three
        do jj= 3,lmax
          rad2int( 1,jj,ii)= rad2int(1,jj-2,ii)-t01*rad2int(1,jj-1,ii-1)
          rad2int( 2,jj,ii)= rad2int(2,jj-2,ii)-t01*rad2int(2,jj-1,ii-1)
          t01= t01+alfbessel*two
        enddo
        t02= betbessel*three
        do kk= 3,lmax
          rad2int(kk, 1,ii)= rad2int(kk-2,1,ii)-t02*rad2int(kk-1,1,ii-1)
          rad2int(kk, 2,ii)= rad2int(kk-2,2,ii)-t02*rad2int(kk-1,2,ii-1)
          t02= t02+betbessel*two
        enddo
        t01= alfbessel*three
        do jj= 3,lmax
          t02= betbessel*three
          do kk= 3,jj-1
            rad2int(kk,jj,ii)= rad2int(kk,jj-2,ii)-t01*rad2int(kk,jj-1,ii-1)
            t02= t02+betbessel*two
          enddo
          do kk= jj,lmax
            rad2int(kk,jj,ii)= rad2int(kk-2,jj,ii)-t02*rad2int(kk-1,jj,ii-1)
            t02= t02+betbessel*two
          enddo
          t01= t01+alfbessel*two
        enddo
      enddo
!
      return
end
