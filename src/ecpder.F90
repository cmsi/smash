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
!---------------------------------------------------------
  subroutine gradoneeiecp(egrad1,fulldmtrx,nproc,myrank)
!---------------------------------------------------------
!
! Calculate ECP derivative terms and add them into energy gradient matrix
!
! Inout  : egrad1 (one electron gradient matrix)
!
      use modecp, only : nterm1, nterm2, maxangecp
      use modbasis, only : nao, mtype, nshell
      use modmolecule, only : natom
      implicit none
      integer,intent(in) :: nproc, myrank
      integer :: maxfunc(0:6)=(/1,3,6,10,15,21,28/)
      integer :: maxbasis, numtbasis, maxdim, llmax, maxecpdim, nsizecp1, nsizecp2, ish, jsh
      integer,allocatable :: label1ecp(:), num1ecp(:), label2ecp(:), num2ecp(:)
      real(8),intent(in) :: fulldmtrx(nao*nao)
      real(8),intent(inout) :: egrad1(natom*3)
      real(8),allocatable :: term1ecp(:), term2ecp(:), term0ecp(:), xyzintecp(:)
!
      maxbasis= maxval(mtype(1:nshell))+1
      maxdim= maxfunc(maxbasis)
      llmax= maxval(maxangecp(1:natom))
      if(maxbasis >= 6) then
        write(*,'(" Error! This program supports up to g function in ecp derivative calculation.")')
        call iabort
      endif
      if(llmax >= 5) then
        write(*,'(" Error! This program supports up to SPDFG core potentials in ecp calculation.")')
        call iabort
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
!$OMP parallel reduction(+:egrad1)
      do ish= nshell-myrank,1,-nproc
!$OMP do
        do jsh= 1,nshell
          call calcdintecp(egrad1,fulldmtrx,term1ecp,term2ecp,term0ecp,xyzintecp, &
&                          label1ecp,label2ecp,num1ecp,num2ecp,numtbasis,ish,jsh,maxdim)
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


!-------------------------------------------------------------------------------------
  subroutine calcdintecp(egrad1,fulldmtrx,term1ecp,term2ecp,term0ecp,xyzintecp, &
&                        label1ecp,label2ecp,num1ecp,num2ecp,numtbasis,ish,jsh,len1)
!-------------------------------------------------------------------------------------
!
! Driver of ECP derivative terms
!
! In : fulldmtrx (density matrix)
!    : ish, jsh  (shell indices)
! Inout : egrad1 (energy gradient value)
!
      use modparam, only : mxprsh
      use modbasis, only : locatom, locprim, locbf, mprim, mbf, mtype, ex, coeff, nao
      use modmolecule, only : natom, coord
      use modecp, only : maxangecp, mtypeecp, locecp, mprimecp, execp, coeffecp, nterm1, nterm2
      implicit none
      integer,intent(in) :: label1ecp(nterm1*9), label2ecp(nterm2*6), num1ecp(*), num2ecp(*)
      integer,intent(in) :: numtbasis, ish, jsh, len1
      integer :: nangij(2), nprimij(2), nbfij(2), iatom, jatom, katom, iloc, jloc, ilocbf, jlocbf
      integer :: iprim, jprim, i, j, k, ii, ncart(0:6)
      integer :: ijkatom(3), nangkecp(mxprsh,6), nprimkecp(6), lmaxecp
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00, half=0.5D+00, three=3.0D+00
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
      real(8),intent(in) :: fulldmtrx(nao,nao), term1ecp(nterm1), term2ecp(nterm2)
      real(8),intent(in) :: term0ecp(*), xyzintecp(25*25*25)
      real(8),intent(inout) :: egrad1(3,natom)
      real(8) :: exij(mxprsh,2), coij(mxprsh,2), coordijk(3,3), ecpint1(len1,len1)
      real(8) :: ecpint2(len1,len1), decpint(28,28,3), work(28)
      real(8) :: exkecp(mxprsh,6), cokecp(mxprsh,6)
      data ncart/1,3,6,10,15,21,28/
!
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
        write(*,'(" Error! This program supports up to h function in calcdintecp")')
        call iabort
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
        enddo
!
      ijkatom(1)= iatom
      ijkatom(2)= jatom
!
      do katom= 1,natom
        if(maxangecp(katom) == -1) cycle
        if(katom == jatom) cycle
        do i= 1,3
          coordijk(i,3)= coord(i,katom)
        enddo
        nangij(2)= mtype(jsh)+1
        nbfij(2)= ncart(nangij(2))
        do jprim= 1,nprimij(2)
          coij(jprim,2)= coeff(jloc+jprim)*two*ex(jloc+jprim)
        enddo
        ijkatom(3)= katom
        lmaxecp= maxangecp(katom)
        do i= 1,maxangecp(katom)+1
          nprimkecp(i)= mprimecp(i-1,katom)
          do jprim= 1,nprimkecp(i)
            exkecp(jprim,i)= execp(locecp(i-1,katom)+jprim)
            cokecp(jprim,i)= coeffecp(locecp(i-1,katom)+jprim)
            nangkecp(jprim,i)= mtypeecp(locecp(i-1,katom)+jprim)
          enddo
        enddo
        call intecp(ecpint1,exij,coij,coordijk,nprimij,nangij,nbfij, &
&                   exkecp,cokecp,nprimkecp,nangkecp,lmaxecp, &
&                   term1ecp,term2ecp,term0ecp,xyzintecp,label1ecp,label2ecp, &
&                   num1ecp,num2ecp,numtbasis,ijkatom,len1,mxprsh)
!
        if(mtype(jsh) >= 1) then
          nangij(2)= mtype(jsh)-1
          nbfij(2)= ncart(nangij(2))
          do jprim= 1,nprimij(2)
            coij(jprim,2) = coeff(jloc+jprim)
          enddo
          call intecp(ecpint2,exij,coij,coordijk,nprimij,nangij,nbfij, &
&                     exkecp,cokecp,nprimkecp,nangkecp,lmaxecp, &
&                     term1ecp,term2ecp,term0ecp,xyzintecp,label1ecp,label2ecp, &
&                     num1ecp,num2ecp,numtbasis,ijkatom,len1,mxprsh)
        else
          ecpint2(1:nbfij(2),1:nbfij(1))= zero
        endif
!
        select case(mtype(jsh))
          case (0)
            do i= 1,nbfij(1)
              decpint(1,i,1)= ecpint1(1,i)
              decpint(1,i,2)= ecpint1(2,i)
              decpint(1,i,3)= ecpint1(3,i)
            enddo
          case (1)
            do i= 1,nbfij(1)
              decpint(1,i,1)= ecpint1(1,i)          -ecpint2(1,i)
              decpint(2,i,1)= ecpint1(2,i)*sqrtthird
              decpint(3,i,1)= ecpint1(3,i)*sqrtthird
              decpint(1,i,2)= ecpint1(2,i)*sqrtthird
              decpint(2,i,2)= ecpint1(4,i)          -ecpint2(1,i)
              decpint(3,i,2)= ecpint1(5,i)*sqrtthird
              decpint(1,i,3)= ecpint1(3,i)*sqrtthird
              decpint(2,i,3)= ecpint1(5,i)*sqrtthird
              decpint(3,i,3)= ecpint1(6,i)          -ecpint2(1,i)
            enddo
          case (2)
            do i= 1,nbfij(1)
              decpint(1,i,1)= ecpint1( 1,i)           -ecpint2(1,i)*two
              decpint(2,i,1)= ecpint1( 2,i)*sqrt3fifth-ecpint2(2,i)*sqrt3
              decpint(3,i,1)= ecpint1( 3,i)*sqrt3fifth-ecpint2(3,i)*sqrt3
              decpint(4,i,1)= ecpint1( 4,i)*sqrtfifth
              decpint(5,i,1)= ecpint1( 5,i)*sqrtfifth
              decpint(6,i,1)= ecpint1( 6,i)*sqrtfifth
              decpint(1,i,2)= ecpint1( 2,i)*sqrtfifth
              decpint(2,i,2)= ecpint1( 4,i)*sqrt3fifth-ecpint2(1,i)*sqrt3
              decpint(3,i,2)= ecpint1( 5,i)*sqrtfifth
              decpint(4,i,2)= ecpint1( 7,i)           -ecpint2(2,i)*two
              decpint(5,i,2)= ecpint1( 8,i)*sqrt3fifth-ecpint2(3,i)*sqrt3
              decpint(6,i,2)= ecpint1( 9,i)*sqrtfifth
              decpint(1,i,3)= ecpint1( 3,i)*sqrtfifth
              decpint(2,i,3)= ecpint1( 5,i)*sqrtfifth
              decpint(3,i,3)= ecpint1( 6,i)*sqrt3fifth-ecpint2(1,i)*sqrt3
              decpint(4,i,3)= ecpint1( 8,i)*sqrtfifth
              decpint(5,i,3)= ecpint1( 9,i)*sqrt3fifth-ecpint2(2,i)*sqrt3
              decpint(6,i,3)= ecpint1(10,i)           -ecpint2(3,i)*two
            enddo
          case (3)
            do i= 1,nbfij(1)
              decpint( 1,i,1)= ecpint1( 1,i)            -ecpint2(1,i)*three
              decpint( 2,i,1)= ecpint1( 2,i)*sqrtseventh-ecpint2(2,i)*two*sqrtthird
              decpint( 3,i,1)= ecpint1( 3,i)*sqrtseventh-ecpint2(3,i)*two*sqrtthird
              decpint( 4,i,1)= ecpint1( 4,i)*sqrt3inv35 -ecpint2(4,i)
              decpint( 5,i,1)= ecpint1( 5,i)*sqrtinv35  -ecpint2(5,i)*sqrtthird
              decpint( 6,i,1)= ecpint1( 6,i)*sqrt3inv35 -ecpint2(6,i)
              decpint( 7,i,1)= ecpint1( 7,i)*sqrtseventh
              decpint( 8,i,1)= ecpint1( 8,i)*sqrtinv35
              decpint( 9,i,1)= ecpint1( 9,i)*sqrtinv35
              decpint(10,i,1)= ecpint1(10,i)*sqrtseventh
              decpint( 1,i,2)= ecpint1( 2,i)*sqrtseventh
              decpint( 2,i,2)= ecpint1( 4,i)*sqrt3inv35 -ecpint2(1,i)
              decpint( 3,i,2)= ecpint1( 5,i)*sqrtinv35
              decpint( 4,i,2)= ecpint1( 7,i)*sqrtseventh-ecpint2(2,i)*two*sqrtthird
              decpint( 5,i,2)= ecpint1( 8,i)*sqrtinv35  -ecpint2(3,i)*sqrtthird
              decpint( 6,i,2)= ecpint1( 9,i)*sqrtinv35
              decpint( 7,i,2)= ecpint1(11,i)            -ecpint2(4,i)*three
              decpint( 8,i,2)= ecpint1(12,i)*sqrtseventh-ecpint2(5,i)*two*sqrtthird
              decpint( 9,i,2)= ecpint1(13,i)*sqrt3inv35 -ecpint2(6,i)
              decpint(10,i,2)= ecpint1(14,i)*sqrtseventh
              decpint( 1,i,3)= ecpint1( 3,i)*sqrtseventh
              decpint( 2,i,3)= ecpint1( 5,i)*sqrtinv35
              decpint( 3,i,3)= ecpint1( 6,i)*sqrt3inv35 -ecpint2(1,i)
              decpint( 4,i,3)= ecpint1( 8,i)*sqrtinv35
              decpint( 5,i,3)= ecpint1( 9,i)*sqrtinv35  -ecpint2(2,i)*sqrtthird
              decpint( 6,i,3)= ecpint1(10,i)*sqrtseventh-ecpint2(3,i)*two*sqrtthird
              decpint( 7,i,3)= ecpint1(12,i)*sqrtseventh
              decpint( 8,i,3)= ecpint1(13,i)*sqrt3inv35 -ecpint2(4,i)
              decpint( 9,i,3)= ecpint1(14,i)*sqrtseventh-ecpint2(5,i)*two*sqrtthird
              decpint(10,i,3)= ecpint1(15,i)            -ecpint2(6,i)*three
            enddo
        case (4)
          do i= 1,nbfij(1)
            decpint( 1,i,1)= ecpint1( 1,i)           -ecpint2( 1,i)*four
            decpint( 2,i,1)= ecpint1( 2,i)*third     -ecpint2( 2,i)*three*sqrtfifth
            decpint( 3,i,1)= ecpint1( 3,i)*third     -ecpint2( 3,i)*three*sqrtfifth
            decpint( 4,i,1)= ecpint1( 4,i)*sqrtinv21 -ecpint2( 4,i)*two*sqrtfifth
            decpint( 5,i,1)= ecpint1( 5,i)*sqrtinv63 -ecpint2( 5,i)*two*sqrtinv15
            decpint( 6,i,1)= ecpint1( 6,i)*sqrtinv21 -ecpint2( 6,i)*two*sqrtfifth
            decpint( 7,i,1)= ecpint1( 7,i)*sqrtinv21 -ecpint2( 7,i)
            decpint( 8,i,1)= ecpint1( 8,i)*sqrtinv105-ecpint2( 8,i)*sqrtfifth
            decpint( 9,i,1)= ecpint1( 9,i)*sqrtinv105-ecpint2( 9,i)*sqrtfifth
            decpint(10,i,1)= ecpint1(10,i)*sqrtinv21 -ecpint2(10,i)
            decpint(11,i,1)= ecpint1(11,i)*third
            decpint(12,i,1)= ecpint1(12,i)*sqrtinv63
            decpint(13,i,1)= ecpint1(13,i)*sqrtinv105
            decpint(14,i,1)= ecpint1(14,i)*sqrtinv63
            decpint(15,i,1)= ecpint1(15,i)*third
            decpint( 1,i,2)= ecpint1( 2,i)*third
            decpint( 2,i,2)= ecpint1( 4,i)*sqrtinv21 -ecpint2( 1,i)
            decpint( 3,i,2)= ecpint1( 5,i)*sqrtinv63
            decpint( 4,i,2)= ecpint1( 7,i)*sqrtinv21 -ecpint2( 2,i)*two*sqrtfifth
            decpint( 5,i,2)= ecpint1( 8,i)*sqrtinv105-ecpint2( 3,i)*sqrtfifth
            decpint( 6,i,2)= ecpint1( 9,i)*sqrtinv105
            decpint( 7,i,2)= ecpint1(11,i)*third     -ecpint2( 4,i)*three*sqrtfifth
            decpint( 8,i,2)= ecpint1(12,i)*sqrtinv63 -ecpint2( 5,i)*two*sqrtinv15
            decpint( 9,i,2)= ecpint1(13,i)*sqrtinv105-ecpint2( 6,i)*sqrtfifth
            decpint(10,i,2)= ecpint1(14,i)*sqrtinv63
            decpint(11,i,2)= ecpint1(16,i)           -ecpint2( 7,i)*four
            decpint(12,i,2)= ecpint1(17,i)*third     -ecpint2( 8,i)*three*sqrtfifth
            decpint(13,i,2)= ecpint1(18,i)*sqrtinv21 -ecpint2( 9,i)*two*sqrtfifth
            decpint(14,i,2)= ecpint1(19,i)*sqrtinv21 -ecpint2(10,i)
            decpint(15,i,2)= ecpint1(20,i)*third
            decpint( 1,i,3)= ecpint1( 3,i)*third
            decpint( 2,i,3)= ecpint1( 5,i)*sqrtinv63
            decpint( 3,i,3)= ecpint1( 6,i)*sqrtinv21 -ecpint2( 1,i)
            decpint( 4,i,3)= ecpint1( 8,i)*sqrtinv105
            decpint( 5,i,3)= ecpint1( 9,i)*sqrtinv105-ecpint2( 2,i)*sqrtfifth
            decpint( 6,i,3)= ecpint1(10,i)*sqrtinv21 -ecpint2( 3,i)*two*sqrtfifth
            decpint( 7,i,3)= ecpint1(12,i)*sqrtinv63
            decpint( 8,i,3)= ecpint1(13,i)*sqrtinv105-ecpint2( 4,i)*sqrtfifth
            decpint( 9,i,3)= ecpint1(14,i)*sqrtinv63 -ecpint2( 5,i)*two*sqrtinv15
            decpint(10,i,3)= ecpint1(15,i)*third     -ecpint2( 6,i)*three*sqrtfifth
            decpint(11,i,3)= ecpint1(17,i)*third
            decpint(12,i,3)= ecpint1(18,i)*sqrtinv21 -ecpint2( 7,i)
            decpint(13,i,3)= ecpint1(19,i)*sqrtinv21 -ecpint2( 8,i)*two*sqrtfifth
            decpint(14,i,3)= ecpint1(20,i)*third     -ecpint2( 9,i)*three*sqrtfifth
            decpint(15,i,3)= ecpint1(21,i)           -ecpint2(10,i)*four
          enddo
        case (5)
          do i= 1,nbfij(1)
            decpint( 1,i,1)= ecpint1( 1,i)            -ecpint2( 1,i)*five
            decpint( 2,i,1)= ecpint1( 2,i)*sqrtinv11  -ecpint2( 2,i)*four*sqrtseventh
            decpint( 3,i,1)= ecpint1( 3,i)*sqrtinv11  -ecpint2( 3,i)*four*sqrtseventh
            decpint( 4,i,1)= ecpint1( 4,i)*sqrtinv33  -ecpint2( 4,i)*three*sqrt3inv35
            decpint( 5,i,1)= ecpint1( 5,i)*sqrtinv99  -ecpint2( 5,i)*three*sqrtinv35
            decpint( 6,i,1)= ecpint1( 6,i)*sqrtinv33  -ecpint2( 6,i)*three*sqrt3inv35
            decpint( 7,i,1)= ecpint1( 7,i)*sqrt5inv231-ecpint2( 7,i)*two*sqrtseventh
            decpint( 8,i,1)= ecpint1( 8,i)*sqrtinv231 -ecpint2( 8,i)*two*sqrtinv35
            decpint( 9,i,1)= ecpint1( 9,i)*sqrtinv231 -ecpint2( 9,i)*two*sqrtinv35
            decpint(10,i,1)= ecpint1(10,i)*sqrt5inv231-ecpint2(10,i)*two*sqrtseventh
            decpint(11,i,1)= ecpint1(11,i)*sqrtinv33  -ecpint2(11,i)
            decpint(12,i,1)= ecpint1(12,i)*sqrtinv231 -ecpint2(12,i)*sqrtseventh
            decpint(13,i,1)= ecpint1(13,i)*sqrtinv385 -ecpint2(13,i)*sqrt3inv35
            decpint(14,i,1)= ecpint1(14,i)*sqrtinv231 -ecpint2(14,i)*sqrtseventh
            decpint(15,i,1)= ecpint1(15,i)*sqrtinv33  -ecpint2(15,i)
            decpint(16,i,1)= ecpint1(16,i)*sqrtinv11
            decpint(17,i,1)= ecpint1(17,i)*sqrtinv99
            decpint(18,i,1)= ecpint1(18,i)*sqrtinv231
            decpint(19,i,1)= ecpint1(19,i)*sqrtinv231
            decpint(20,i,1)= ecpint1(20,i)*sqrtinv99
            decpint(21,i,1)= ecpint1(21,i)*sqrtinv11
            decpint( 1,i,2)= ecpint1( 2,i)*sqrtinv11
            decpint( 2,i,2)= ecpint1( 4,i)*sqrtinv33  -ecpint2( 1,i)
            decpint( 3,i,2)= ecpint1( 5,i)*sqrtinv99
            decpint( 4,i,2)= ecpint1( 7,i)*sqrt5inv231-ecpint2( 2,i)*two*sqrtseventh
            decpint( 5,i,2)= ecpint1( 8,i)*sqrtinv231 -ecpint2( 3,i)*sqrtseventh
            decpint( 6,i,2)= ecpint1( 9,i)*sqrtinv231
            decpint( 7,i,2)= ecpint1(11,i)*sqrtinv33  -ecpint2( 4,i)*three*sqrt3inv35
            decpint( 8,i,2)= ecpint1(12,i)*sqrtinv231 -ecpint2( 5,i)*two*sqrtinv35
            decpint( 9,i,2)= ecpint1(13,i)*sqrtinv385 -ecpint2( 6,i)*sqrt3inv35
            decpint(10,i,2)= ecpint1(14,i)*sqrtinv231
            decpint(11,i,2)= ecpint1(16,i)*sqrtinv11  -ecpint2( 7,i)*four*sqrtseventh
            decpint(12,i,2)= ecpint1(17,i)*sqrtinv99  -ecpint2( 8,i)*three*sqrtinv35
            decpint(13,i,2)= ecpint1(18,i)*sqrtinv231 -ecpint2( 9,i)*two*sqrtinv35
            decpint(14,i,2)= ecpint1(19,i)*sqrtinv231 -ecpint2(10,i)*sqrtseventh
            decpint(15,i,2)= ecpint1(20,i)*sqrtinv99
            decpint(16,i,2)= ecpint1(22,i)            -ecpint2(11,i)*five
            decpint(17,i,2)= ecpint1(23,i)*sqrtinv11  -ecpint2(12,i)*four*sqrtseventh
            decpint(18,i,2)= ecpint1(24,i)*sqrtinv33  -ecpint2(13,i)*three*sqrt3inv35
            decpint(19,i,2)= ecpint1(25,i)*sqrt5inv231-ecpint2(14,i)*two*sqrtseventh
            decpint(20,i,2)= ecpint1(26,i)*sqrtinv33  -ecpint2(15,i)
            decpint(21,i,2)= ecpint1(27,i)*sqrtinv11
            decpint( 1,i,3)= ecpint1( 3,i)*sqrtinv11
            decpint( 2,i,3)= ecpint1( 5,i)*sqrtinv99
            decpint( 3,i,3)= ecpint1( 6,i)*sqrtinv33  -ecpint2( 1,i)
            decpint( 4,i,3)= ecpint1( 8,i)*sqrtinv231
            decpint( 5,i,3)= ecpint1( 9,i)*sqrtinv231 -ecpint2( 2,i)*sqrtseventh
            decpint( 6,i,3)= ecpint1(10,i)*sqrt5inv231-ecpint2( 3,i)*two*sqrtseventh
            decpint( 7,i,3)= ecpint1(12,i)*sqrtinv231
            decpint( 8,i,3)= ecpint1(13,i)*sqrtinv385 -ecpint2( 4,i)*sqrt3inv35
            decpint( 9,i,3)= ecpint1(14,i)*sqrtinv231 -ecpint2( 5,i)*two*sqrtinv35
            decpint(10,i,3)= ecpint1(15,i)*sqrtinv33  -ecpint2( 6,i)*three*sqrt3inv35
            decpint(11,i,3)= ecpint1(17,i)*sqrtinv99
            decpint(12,i,3)= ecpint1(18,i)*sqrtinv231 -ecpint2( 7,i)*sqrtseventh
            decpint(13,i,3)= ecpint1(19,i)*sqrtinv231 -ecpint2( 8,i)*two*sqrtinv35
            decpint(14,i,3)= ecpint1(20,i)*sqrtinv99  -ecpint2( 9,i)*three*sqrtinv35
            decpint(15,i,3)= ecpint1(21,i)*sqrtinv11  -ecpint2(10,i)*four*sqrtseventh
            decpint(16,i,3)= ecpint1(23,i)*sqrtinv11
            decpint(17,i,3)= ecpint1(24,i)*sqrtinv33  -ecpint2(11,i)
            decpint(18,i,3)= ecpint1(25,i)*sqrt5inv231-ecpint2(12,i)*two*sqrtseventh
            decpint(19,i,3)= ecpint1(26,i)*sqrtinv33  -ecpint2(13,i)*three*sqrt3inv35
            decpint(20,i,3)= ecpint1(27,i)*sqrtinv11  -ecpint2(14,i)*four*sqrtseventh
            decpint(21,i,3)= ecpint1(28,i)            -ecpint2(15,i)*five
          enddo
        end select
!
        nbfij(2)  = mbf(jsh)
        select case(nbfij(2))
          case(5)
            do k= 1,3
              do i= 1,nbfij(1)
                do j= 1,6
                  work(j)= decpint(j,i,k)
                enddo
                decpint(1,i,k)= work(2)
                decpint(2,i,k)= work(5)
                decpint(3,i,k)=(work(6)*two-work(1)-work(4))*half
                decpint(4,i,k)= work(3)
                decpint(5,i,k)=(work(1)-work(4))*sqrt3h
              enddo
            enddo
          case(7)
            do k= 1,3
              do i= 1,nbfij(1)
                do j= 1,10
                  work(j)= decpint(j,i,k)
                enddo
                decpint(1,i,k)=(-work(7)+work(2)*three                   )*facf1
                decpint(2,i,k)=  work(5)                                  *facf2
                decpint(3,i,k)=(-work(7)-work(2)+work(9)*four            )*facf3
                decpint(4,i,k)=( work(10)*two-work(3)*three-work(8)*three)*half
                decpint(5,i,k)=(-work(1)-work(4)+work(6)*four            )*facf3
                decpint(6,i,k)=( work(3)-work(8)                         )*facf4
                decpint(7,i,k)=( work(1)-work(4)*three                   )*facf1
              enddo
            enddo
          case(10)
            do k= 1,3
              do i= 1,nbfij(1)
                decpint(2,i,k)= decpint(2,i,k)*sqrt5
                decpint(3,i,k)= decpint(3,i,k)*sqrt5
                decpint(4,i,k)= decpint(4,i,k)*sqrt5
                decpint(5,i,k)= decpint(5,i,k)*sqrt15
                decpint(6,i,k)= decpint(6,i,k)*sqrt5
                decpint(8,i,k)= decpint(8,i,k)*sqrt5
                decpint(9,i,k)= decpint(9,i,k)*sqrt5
              enddo
            enddo
          case(9)
            do k= 1,3
              do i= 1,nbfij(1)
                do j= 1,15
                  work(j)= decpint(j,i,k)
                enddo
                decpint(1,i,k)=(work(2)-work(7))*facg1
                decpint(2,i,k)=(-work(12)+work(5)*three)*facg2
                decpint(3,i,k)=(-work(2)-work(7)+work(9)*six)*facg3
                decpint(4,i,k)=(-work(12)*three+work(14)*four-work(5)*three)*facg4
                decpint(5,i,k)=(work(1)*three+work(11)*three+work(15)*eight+work(4)*six &
&                            -work(6)*p24-work(13)*p24)*eighth
                decpint(6,i,k)=(-work(3)*three+work(10)*four-work(8)*three)*facg4
                decpint(7,i,k)=(-work(1)+work(11)+work(6)*six-work(13)*six)*facg5
                decpint(8,i,k)=(work(3)-work(8)*three)*facg2
                decpint(9,i,k)=(work(1)+work(11)-work(4)*six)*facg6
              enddo
            enddo
          case(15)
            do k= 1,3
              do i= 1,nbfij(1)
                decpint( 2,i,k)= decpint( 2,i,k)*sqrt7
                decpint( 3,i,k)= decpint( 3,i,k)*sqrt7
                decpint( 4,i,k)= decpint( 4,i,k)*sqrt35third
                decpint( 5,i,k)= decpint( 5,i,k)*sqrt35
                decpint( 6,i,k)= decpint( 6,i,k)*sqrt35third
                decpint( 7,i,k)= decpint( 7,i,k)*sqrt7
                decpint( 8,i,k)= decpint( 8,i,k)*sqrt35
                decpint( 9,i,k)= decpint( 9,i,k)*sqrt35
                decpint(10,i,k)= decpint(10,i,k)*sqrt7
                decpint(12,i,k)= decpint(12,i,k)*sqrt7
                decpint(13,i,k)= decpint(13,i,k)*sqrt35third
                decpint(14,i,k)= decpint(14,i,k)*sqrt7
              enddo
            enddo
          case(11)
            do k= 1,3
              do i= 1,nbfij(1)
                do j= 1,21
                  work(j)= decpint(j,i,k)
                enddo
                decpint( 1,i,k)=(work(2)*five-work(7)*ten+work(16))*fach1
                decpint( 2,i,k)=(work(5)*four-work(12)*four)*fach2
                decpint( 3,i,k)=(-work(2)*three-work(7)*two+work(9)*p24+work(16) &
&                               -work(18)*eight)*fach3
                decpint( 4,i,k)=(-work(5)*two-work(12)*two+work(14)*four)*fach4
                decpint( 5,i,k)=(work(2)+work(7)*two-work(9)*twelve+work(16)-work(18)*twelve &
&                               +work(20)*eight)*fach5
                decpint( 6,i,k)=(work(3)*p15+work(8)*p30-work(10)*p40+work(17)*p15 &
&                               -work(19)*p40+work(21)*eight)*eighth
                decpint( 7,i,k)=(work(1)+work(4)*two-work(6)*twelve+work(11)-work(13)*twelve &
&                               +work(15)*eight)*fach5
                decpint( 8,i,k)=(-work(3)+work(10)*two+work(17)-work(19)*two)*fach4
                decpint( 9,i,k)=(-work(1)+work(4)*two+work(6)*eight+work(11)*three &
&                               -work(13)*p24)*fach3
                decpint(10,i,k)=(work(3)-work(8)*six+work(17))*fach2
                decpint(11,i,k)=(work(1)-work(4)*ten+work(11)*five)*fach1
              enddo
            enddo
          case(21)
            do k= 1,3
              do i= 1,nbfij(1)
                decpint( 2,i,k)= decpint( 2,i,k)*three
                decpint( 3,i,k)= decpint( 3,i,k)*three
                decpint( 4,i,k)= decpint( 4,i,k)*sqrt21
                decpint( 5,i,k)= decpint( 5,i,k)*sqrt63
                decpint( 6,i,k)= decpint( 6,i,k)*sqrt21
                decpint( 7,i,k)= decpint( 7,i,k)*sqrt21
                decpint( 8,i,k)= decpint( 8,i,k)*sqrt105
                decpint( 9,i,k)= decpint( 9,i,k)*sqrt105
                decpint(10,i,k)= decpint(10,i,k)*sqrt21
                decpint(11,i,k)= decpint(11,i,k)*three
                decpint(12,i,k)= decpint(12,i,k)*sqrt63
                decpint(13,i,k)= decpint(13,i,k)*sqrt105
                decpint(14,i,k)= decpint(14,i,k)*sqrt63
                decpint(15,i,k)= decpint(15,i,k)*three
                decpint(17,i,k)= decpint(17,i,k)*three
                decpint(18,i,k)= decpint(18,i,k)*sqrt21
                decpint(19,i,k)= decpint(19,i,k)*sqrt21
                decpint(20,i,k)= decpint(20,i,k)*three
              enddo
            enddo
        end select
!
        do i= 1,nbfij(1)
          ii= ilocbf+i
          do j= 1,nbfij(2)
            egrad1(1,jatom)= egrad1(1,jatom)+fulldmtrx(jlocbf+j,ii)*decpint(j,i,1)
            egrad1(2,jatom)= egrad1(2,jatom)+fulldmtrx(jlocbf+j,ii)*decpint(j,i,2)
            egrad1(3,jatom)= egrad1(3,jatom)+fulldmtrx(jlocbf+j,ii)*decpint(j,i,3)
            egrad1(1,katom)= egrad1(1,katom)-fulldmtrx(jlocbf+j,ii)*decpint(j,i,1)
            egrad1(2,katom)= egrad1(2,katom)-fulldmtrx(jlocbf+j,ii)*decpint(j,i,2)
            egrad1(3,katom)= egrad1(3,katom)-fulldmtrx(jlocbf+j,ii)*decpint(j,i,3)
          enddo
        enddo
!
      enddo
      return
end
