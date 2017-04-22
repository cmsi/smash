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
!---------------------------------------------------------------
  subroutine calcrmp2(cmo,energymo,xint,nproc,myrank,mpi_comm)
!---------------------------------------------------------------
!
! Driver of restricted MP2 calculation
!
! In  : cmo     (MO coefficient matrix)
!       energymo(MO energies)
!       xint    (Exchange integral matrix)
!       
      use modparallel, only : master
      use modbasis, only : nshell, nao, mbf
      use modmolecule, only : neleca, nmo
      use modenergy, only : escf, emp2, escsmp2
      use modmp2, only : ncore, nvfz
      implicit none
      integer,intent(in) :: nproc, myrank, mpi_comm
      integer :: maxdim, noac, nvac, noac3, iproc, icount
      integer :: idis(0:nproc-1,4), maxsize, ish, jsh, msize, memneed
      integer :: numocc3, npass
      real(8),parameter :: zero=0.0D+00, three=3.0D+00, p12=1.2D+00
      real(8),intent(in) :: cmo(nao,nao), energymo(nmo), xint(nshell*(nshell+1)/2)
      real(8) :: emp2st(2), emp2stsum(2)
!
      emp2= zero
      emp2st(:)= zero
      maxdim= maxval(mbf(1:nshell))
      noac= neleca-ncore
      nvac= nmo-neleca-nvfz
      noac3= noac*(noac+1)/2
!
      if(master) then
        write(*,'(" ----------------------------------------------")')
        write(*,'("   MP2 calculation ")')
        write(*,'(" ----------------------------------------------")')
        write(*,'("     Ncore=",i4,",   Nvfz=",i4)')ncore,nvfz
        write(*,'(" ----------------------------------------------")')
        write(*,'("     Number of basis functions         =",i5)')nao
        write(*,'("     Number of basis shells            =",i5)')nshell
        write(*,'("     Number of correlated occupied MOs =",i5)')noac
        write(*,'("     Number of active virtual MOs      =",i5)')nvac
        write(*,'(" ----------------------------------------------")')
      endif
!
      icount= 0
      idis(:,:)= 0
      ish= 1
      jsh= 1
!
! Calculate distributed AO indices
!
      do iproc= 0,nproc-1
        idis(iproc,1)= ish
        idis(iproc,2)= jsh
        jsh= jsh+1
        if(jsh > nshell) then
          ish= ish+1
          jsh= 1
        endif
      enddo
!
      iproc= 0
      if(nproc /= 1) then
        do ish= 1,nshell
          do jsh= 1,nshell
            idis(iproc,3)= idis(iproc,3)+mbf(ish)*mbf(jsh)
            idis(iproc,4)= idis(iproc,4)+1
            icount= icount+1
            if(mod(icount,nproc) == 0) then
              iproc= iproc+2
            else
              iproc= iproc+1
            endif
            if(iproc >= nproc) iproc= iproc-nproc
          enddo
        enddo
      else
        idis(iproc,3)= nao*nao
        idis(iproc,4)= nshell*nshell
      endif
      maxsize= maxval(idis(0:nproc-1,3))
!
! Check available memory size and judge the number of passes
!
      call memrest(msize)
      memneed= max(nao*noac+noac*maxdim**3+noac*nao*maxdim**2,nao*nao*2)
      numocc3=(msize-memneed)/maxsize
!
      if(numocc3 <= 0) then
        if(master) then
          call memset(maxsize+memneed)
          call iabort
        endif
!
! Single pass
!
      elseif(numocc3 >= noac3) then
        if(master) write(*,'(" == Single pass calculation ==")')
        call mp2single(emp2st,cmo,energymo,xint,noac,nvac,ncore,maxsize,maxdim,idis, &
&                      nproc,myrank,mpi_comm)
!
! Multiple pass
!
      else
        npass=(noac3-1)/numocc3+1
        numocc3=(noac3-1)/npass+1
        if(master) then
          write(*,'(" == Multiple pass calculation ==")')
          write(*,'("    Number of passes :",i5)')npass
        endif
        call mp2multi(emp2st,cmo,energymo,xint,noac,nvac,ncore,maxsize,maxdim,idis, &
&                     npass,numocc3,nproc,myrank,mpi_comm)
      endif
!
      call para_allreducer(emp2st,emp2stsum,2,mpi_comm)
      emp2= emp2stsum(1)+emp2stsum(2)
      escsmp2= emp2stsum(1)*p12+emp2stsum(2)/three
!
      if(master) then
        write(*,'(" -------------------------------------------------")')
        write(*,'("   HF Energy                  =",f17.9)')escf
        write(*,'("   MP2 Correlation Energy     =",f17.9)')emp2
        write(*,'("   HF + MP2 Energy            =",f17.9)')escf+emp2
        write(*,'("   SCS-MP2 Correlation Energy =",f17.9)')escsmp2
        write(*,'("   HF + SCS-MP2 Energy        =",f17.9)')escf+escsmp2
        write(*,'(" -------------------------------------------------")')
      endif
      return
end


!-----------------------
  function ncorecalc()
!-----------------------
!
! Calculate the number of core MOs
!
      use modmolecule, only : natom, numatomic
      use modecp, only : izcore
      implicit none
      integer :: ncorecalc, ncore, iatom, numcore(-9:137)
      data numcore/12*0, 8*1, 8*5, 18*9, 18*18, 32*27, 32*43, 19*59/
!
      ncore= 0
      do iatom= 1,natom
        if(izcore(iatom) == 0) then
          ncore= ncore+numcore(numatomic(iatom))
        endif
      enddo
      ncorecalc= ncore
      return
end


!---------------------------------------------------------------------------------------
  subroutine mp2single(emp2st,cmo,energymo,xint,noac,nvac,ncore,maxsize,maxdim,idis, &
&                      nproc,myrank,mpi_comm)
!---------------------------------------------------------------------------------------
!
! Driver of single pass MP2 energy calculation
!
! In  : cmo     (MO coefficient matrix)
!       energymo(MO energies)
!       xint    (Exchange integral matrix)
!       noac    (Number of active occupied MOs)
!       nvac    (Number of active virtual MOs)
!       ncore   (Number of core occupied MOs)
!       maxsize (Maximum size of basis function pairs)
!       maxdim  (Maximum size of basis functions in a shell)
!       idis    (Information for parallelization)
! Out : emp2st  (Partial MP2 energies)
!
      use modparallel, only : master
      use modbasis, only : nshell, nao
      use modmolecule, only : neleca, nmo
      implicit none
      integer,intent(in) :: noac, nvac, ncore, maxsize, maxdim, nproc, myrank, mpi_comm
      integer,intent(in) :: idis(0:nproc-1,4)
      integer :: noac3, nao2, msize, mlsize
      real(8),intent(in) :: cmo(nao,nao), energymo(nmo), xint(nshell*(nshell+1)/2)
      real(8),intent(inout) :: emp2st(2)
      real(8),allocatable :: trint2(:), cmowrk(:), trint1a(:), trint1b(:)
      real(8),allocatable :: trint3(:), trint4(:)
!
      noac3= noac*(noac+1)/2
      nao2 = nao*nao
!
      call memset(maxsize*noac3)
      allocate(trint2(maxsize*noac3))
      call memrest(msize)
      mlsize=(msize-nao*noac-noac*maxdim**3)/(noac*nao)
      if(mlsize > maxsize) then
        mlsize= maxsize
      elseif(mlsize < maxdim*maxdim) then
        mlsize= maxdim*maxdim
      endif
      call memset(nao*noac+noac*maxdim**3+mlsize*noac*nao)
      allocate(cmowrk(nao*noac),trint1a(noac*maxdim**3),trint1b(mlsize*noac*nao))
!
! AO intengral generation and first and second integral transformations
!
      if(master) write(*,'("    Start first and second integral transformations")')
      call mp2trans12(cmo(1,ncore+1),cmowrk,trint1a,trint1b,trint2,xint, &
&                     mlsize,noac,maxdim,idis,nproc,myrank)
!
      deallocate(cmowrk,trint1a,trint1b)
      call memunset(nao*noac+noac*maxdim**3+mlsize*noac*nao)
!
      call memset(2*nao2)
      allocate(trint3(nao2),trint4(nao2))
!
! Third and fourth integral transformations and MP2 energy calculation
!
      if(master) write(*,'("    Start third and fourth integral transformations")')
      call mp2trans34(cmo(1,neleca+1),energymo,trint2,trint3,trint4,emp2st,noac,nvac,ncore, &
&                     idis,nproc,myrank,mpi_comm)
!
      deallocate(trint3,trint4)
      call memunset(2*nao2)
      deallocate(trint2)
      call memunset(maxsize*noac3)
!
      return
end


!---------------------------------------------------------------------------------------
  subroutine mp2multi(emp2st,cmo,energymo,xint,noac,nvac,ncore,maxsize,maxdim,idis, &
&                     npass,numocc3,nproc,myrank,mpi_comm)
!---------------------------------------------------------------------------------------
!
! Driver of single pass MP2 energy calculation
!
! In  : cmo     (MO coefficient matrix)
!       energymo(MO energies)
!       xint    (Exchange integral matrix)
!       noac    (Number of active occupied MOs)
!       nvac    (Number of active virtual MOs)
!       ncore   (Number of core occupied MOs)
!       maxsize (Maximum size of basis function pairs)
!       maxdim  (Maximum size of basis functions in a shell)
!       idis    (Information for parallelization)
!       npass   (Number of passes)
!       numocc3 (Number of active occupied MO pairs per pass)
! Out : emp2st  (Partial MP2 energies)
!
      use modparallel, only : master
      use modbasis, only : nshell, nao
      use modmolecule, only : neleca, nmo
      implicit none
      integer,intent(in) :: noac, nvac, ncore, maxsize, maxdim, nproc, myrank, mpi_comm
      integer,intent(in) :: idis(0:nproc-1,4), npass, numocc3
      integer :: noac3, nao2, msize, mlsize, ijindex(4,npass), icount, ipass, moi, moj, numij
      real(8),intent(in) :: cmo(nao,nao), energymo(nmo), xint(nshell*(nshell+1)/2)
      real(8),intent(inout) :: emp2st(2)
      real(8),allocatable :: trint2(:), cmowrk(:), trint1a(:), trint1b(:)
      real(8),allocatable :: trint3(:), trint4(:)
!
      noac3= noac*(noac+1)/2
      nao2 = nao*nao
      icount= 0
      ipass= 1
!
! ijindex(1,*): first moi
! ijindex(2,*): first moj
! ijindex(3,*): last moi
! ijindex(4,*): last moj
!
      do moi= 1,noac
        do moj= 1,moi
          icount=icount+1
          if(icount == numocc3) then
            ijindex(3,ipass)= moi
            ijindex(4,ipass)= moj
            icount= 0
            ipass= ipass+1
          endif
        enddo
      enddo
      ijindex(3,npass)=noac
      ijindex(4,npass)=noac
      ijindex(1,1)=1
      ijindex(2,1)=1
      do ipass=2,npass
        moi= ijindex(3,ipass-1)
        moj= ijindex(4,ipass-1)+1
        if(moj > moi) then
          moi= moi+1
          moj= 1
        endif
        ijindex(1,ipass)=moi
        ijindex(2,ipass)=moj
      enddo
!
      call memset(maxsize*numocc3)
      allocate(trint2(maxsize*numocc3))
      call memrest(msize)
      mlsize=(msize-nao*noac-noac*maxdim**3)/(noac*nao)
      if(mlsize > maxsize) then
        mlsize= maxsize
      elseif(mlsize < maxdim*maxdim) then
        mlsize= maxdim*maxdim
      endif
!
!
! Start multiple pass
!
      numij= numocc3
      do ipass= 1,npass
        if(ipass == npass) numij= noac3-(npass-1)*numocc3
        call memset(nao*noac+noac*maxdim**3+mlsize*noac*nao)
        allocate(cmowrk(nao*noac),trint1a(noac*maxdim**3),trint1b(mlsize*noac*nao))
!
! AO intengral generation and first and second integral transformations
!
        if(master) &
&         write(*,'("    Start first and second integral transformations of Pass",i5)')ipass
        call mp2trans12m(cmo(1,ncore+1),cmowrk,trint1a,trint1b,trint2,xint, &
&                        mlsize,noac,maxdim,idis,numij,ijindex(1,ipass),nproc,myrank)
!
        deallocate(cmowrk,trint1a,trint1b)
        call memunset(nao*noac+noac*maxdim**3+mlsize*noac*nao)
!
        call memset(2*nao2)
        allocate(trint3(nao2),trint4(nao2))
!
! Third and fourth integral transformations and MP2 energy calculation
!
        if(master) &
&         write(*,'("    Start third and fourth integral transformations of Pass",i5)')ipass
        call mp2trans34m(cmo(1,neleca+1),energymo,trint2,trint3,trint4,emp2st,nvac,ncore, &
&                        idis,numij,ijindex(1,ipass),nproc,myrank,mpi_comm)
!
        deallocate(trint3,trint4)
        call memunset(2*nao2)
      enddo
      deallocate(trint2)
      call memunset(maxsize*numocc3)
!
      return
end


!---------------------------------------------------------------------
  subroutine mp2trans12(cmoocc,cmowrk,trint1a,trint1b,trint2,xint, &
&                       mlsize,noac,maxdim,idis,nproc,myrank)
!---------------------------------------------------------------------
!
! Driver of AO intengral generation and first and second integral transformations
! for single pass
!
! In  : cmoocc  (MO coefficient matrix)
!       xint    (Exchange integral matrix)
!       mlsize  (Block size of first transformed integrals)
!       noac    (Number of active occupied MOs)
!       maxdim  (Maximum size of basis functions in a shell)
!       idis    (Information for parallelization)
! Out : trint2  (Second transformed integrals)
! Work: cmowrk  (Transposed MO coefficient matrix)
!       trint1a (First transformed integrals)
!       trint1b (First transformed integrals)
!       
      use modbasis, only : nshell, nao, mbf
      implicit none
      integer,intent(in) :: maxdim, mlsize, noac, nproc, myrank, idis(0:nproc-1,4)
      integer :: ish, ksh, ish1, ksh1, mlcount, mlstart, mlshell, numshell, ii, jcount
      integer :: jcount1
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: cmoocc(nao,noac), xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: cmowrk(noac,nao), trint1a(noac,maxdim,maxdim**2)
      real(8),intent(out) :: trint1b(mlsize*noac*nao)
      real(8),intent(out) :: trint2(idis(myrank,3)*noac*(noac+1)/2)
!
      cmowrk=transpose(cmoocc)
!
      ish= idis(myrank,1)
      ksh= idis(myrank,2)
      ish1= ish
      ksh1= ksh
      jcount= 0
      jcount1= 0
      mlcount= 0
      mlstart= 1
      mlshell=0
      do numshell= 1,idis(myrank,4)
        mlshell= mlshell+1
        mlcount= mlcount+mbf(ish)*mbf(ksh)
        if(numshell == idis(myrank,4)) then
!
! AO intengral generation and first integral transformation
!
          call transmoint1(trint1a,trint1b,cmowrk,xint,ish1,ksh1,maxdim,noac,jcount1, &
&                          mlshell,mlsize,nproc,myrank)
!
! Second integral transformation
!
          call transmoint2(trint2,trint1b,cmoocc,noac,mlcount,mlstart,mlsize,idis,nproc,myrank)
        else
          if(jcount == myrank) then
            ksh= ksh+2*nproc-1
          else
            ksh= ksh+nproc-1
          endif
          jcount= jcount+1
          if(jcount == nproc) jcount= 0
!
          if(ksh > nshell) then
            do ii= 1,nshell
              ksh= ksh-nshell
              ish= ish+1
              if(ksh <= nshell) exit
            enddo
          endif
          if(mlcount+mbf(ish)*mbf(ksh) > mlsize) then
            call transmoint1(trint1a,trint1b,cmowrk,xint,ish1,ksh1,maxdim,noac,jcount1, &
&                            mlshell,mlsize,nproc,myrank)
            call transmoint2(trint2,trint1b,cmoocc,noac,mlcount,mlstart,mlsize,idis,nproc,myrank)
            mlstart= mlstart+mlcount
            mlshell= 0
            mlcount= 0
            jcount1= jcount
            ish1= ish
            ksh1= ksh
          endif
        endif
      enddo
!
!     if(iprint >= 3) then
!       write(*,'("Time for Tr1=",f9.3,"  Time for Tr2=",f9.3,i3)')tr1,tr2,myrank
!     endif
      return
end


!---------------------------------------------------------------------------------------
  subroutine mp2trans34(cmovir,energymo,trint2,trint3,trint4,emp2st,noac,nvac,ncore, &
&                       idis,nproc,myrank,mpi_comm)
!---------------------------------------------------------------------------------------
!
! Driver of third and fourth integral transformations and MP2 energy calculation
!
! In  : cmovir  (MO coefficient matrix)
!       energymo(MO energies)
!       trint2  (Second transformed integrals)
!       noac    (Number of active occupied MOs)
!       nvac    (Number of virtual MOs)
!       idis    (Information for parallelization)
! Work: trint3  (Third transformed integrals)
!       trint4  (Fourth transformed integrals)
!
      use modbasis, only : nao
      use modmolecule, only : nmo
      implicit none
      integer,intent(in) :: noac, nvac, ncore, nproc, myrank, mpi_comm, idis(0:nproc-1,4)
      integer :: numrecv, iproc, irecv(0:nproc-1), noac3, ncycle, icycle, myij
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: cmovir(nao*nvac), energymo(nmo)
      real(8),intent(in) :: trint2(idis(myrank,3),noac*(noac+1)/2)
      real(8),intent(out) :: trint3(nao*nao), trint4(nao*nao), emp2st(2)
!
      numrecv= 1
      do iproc= 0,nproc-1
        irecv(iproc)= numrecv
        numrecv= numrecv+idis(iproc,3)
      enddo
      noac3= noac*(noac+1)/2
      ncycle=(noac3-1)/nproc+1
!
      do icycle= 1,ncycle    
        call mp2int_sendrecv(trint2,trint3,trint4,icycle,irecv,noac3, &
&                            idis,nproc,myrank,mpi_comm)
!
        myij=(icycle-1)*nproc+1+myrank
        if(myij > noac3) cycle
!
! Third integral transformation
!   (mi|lj) trint4[l,m] -> (ai|lj) trint3[l,a]
!
        call dgemm('N','N',nao,nvac,nao,one,trint4,nao,cmovir,nao,zero,trint3,nao)
!
! Fourth integral transformation
!   (ai|lj) trint3[l,a] -> (ai|bj) trint4[b,a]
!
        call dgemm('T','N',nvac,nvac,nao,one,cmovir,nao,trint3,nao,zero,trint4,nvac)
!
! MP2 energy calculation
!
        call calcrmp2energy(trint4,energymo,emp2st,noac,nvac,ncore,icycle,nproc,myrank)
      enddo
      return
end


!-----------------------------------------------------------------------------------
  subroutine transmoint1(trint1a,trint1b,cmowrk,xint,ish,ksh,maxdim,numi,jcount, &
&                        mlshell,mlsize,nproc,myrank)
!-----------------------------------------------------------------------------------
!
! AO intengral generation and first-quarter integral transformation
!    (mn|ls) -> (mi|ls)
!
! In  : cmowrk  (Transposed MO coefficient matrix)
!       xint    (Exchange integral matrix)
!       maxdim  (Maximum size of basis functions in a shell)
!       numi    (Number of active occupied MOs)
!       jcount  (Counter of basis shell pair)
!       mlshell (Number of shell pairs)
!       mlsize  (Size of trint1b)
! Out : trint1b (First transformed integrals, [s,i,ml])
! Work: trint1a (First transformed integrals, [i,s,ml])
! Inout : ish,ksh (Basis shell indices)
!
      use modbasis, only : nao, nshell, mbf, locbf
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: maxdim, numi, mlshell, mlsize, nproc, myrank
      integer,intent(inout) :: ish, ksh, jcount
      integer :: mlcount, ml, mlindex(3,mlshell), jsh, lsh, nbfi, nbfj, nbfk, nbfl
      integer :: nbfik, locbfj, locbfl, moi, i, j, k, l, ik, kl, ij, ii, jloc, lloc
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: cmowrk(numi,nao), xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: trint1a(numi,maxdim,maxdim*maxdim), trint1b(nao,numi,mlsize)
      real(8) :: twoeri(maxdim,maxdim,maxdim,maxdim)
!
      mlcount= 0
      do ml= 1,mlshell
        mlindex(1,ml)= ish
        mlindex(2,ml)= ksh
        mlindex(3,ml)= mlcount
        mlcount= mlcount+mbf(ish)*mbf(ksh)
        if(jcount == myrank) then
          ksh= ksh+2*nproc-1
        else
          ksh= ksh+nproc-1
        endif
        jcount= jcount+1
        if(jcount == nproc) jcount= 0
!
        if(ksh > nshell) then
          do ii= 1,nshell
            ksh= ksh-nshell
            ish= ish+1
            if(ksh <= nshell) exit
          enddo
        endif
      enddo
!
!$OMP parallel do schedule(dynamic,1) collapse(2) private(ish,ksh,mlcount,nbfi,nbfk,nbfl, &
!$OMP nbfik,locbfl,trint1a,kl,nbfj,locbfj,ij,twoeri,ii,ik,jloc,lloc)
      do ml= 1,mlshell
        do lsh= 1,nshell
          ish= mlindex(1,ml)
          ksh= mlindex(2,ml)
          mlcount= mlindex(3,ml)
          nbfi= mbf(ish)
          nbfk= mbf(ksh)
          nbfl= mbf(lsh)
          nbfik= nbfi*nbfk
          locbfl= locbf(lsh)
          do ik= 1,nbfik
            do l= 1,nbfl
              do moi= 1,numi
                trint1a(moi,l,ik)= zero
              enddo
            enddo
          enddo
          if(ksh >= lsh) then
            kl= ksh*(ksh-1)/2+lsh
          else
            kl= lsh*(lsh-1)/2+ksh
          endif
!
          do jsh= 1,nshell
            nbfj  = mbf(jsh)
            locbfj= locbf(jsh)
            if(ish >= jsh) then
              ij= ish*(ish-1)/2+jsh
            else
              ij= jsh*(jsh-1)/2+ish
            endif
!
! AO integral calculation
!
            if(xint(ij)*xint(kl) < cutint2) cycle
            call calc2eri(twoeri,ksh,lsh,ish,jsh,maxdim)
!
! First integral transformation
!
            do i= 1,nbfi
              ii=(i-1)*nbfk
              do k= 1,nbfk
                ik= ii+k
                do l= 1,nbfl
                  do j= 1,nbfj
                    jloc= locbfj+j
                    if(abs(twoeri(j,i,l,k)) < cutint2) cycle
                    do moi= 1,numi
                      trint1a(moi,l,ik)= trint1a(moi,l,ik)+twoeri(j,i,l,k)*cmowrk(moi,jloc)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
!
! Reorder of first transformed integrals
!
          do ik= 1,nbfik
            do l= 1,nbfl
              lloc= locbfl+l
              do moi= 1,numi
                trint1b(lloc,moi,mlcount+ik)= trint1a(moi,l,ik)
              enddo
            enddo
          enddo
!
        enddo
      enddo
!$OMP end parallel do
!
      return
end
!


!----------------------------------------------------------------------------------------------
  subroutine transmoint2(trint2,trint1b,cmoocc,noac,mlcount,mlstart,mlsize,idis,nproc,myrank)
!----------------------------------------------------------------------------------------------
!
! Second-quarter integral transformation
!    (mi|ls) -> (mi|lj)
!
! In  : trint1b (First transformed integrals, [s,i,ml])
!       cmoocc  (MO coefficient matrix)
!       noac    (Number of active occupied MOs)
!       mlcount (Number of transformed AOs)
!       mlstart (First index of trint2)
!       mlsize  (Size of trint1b)
!       idis    (Information for parallelization)
! Out : trint2  (Second transformed integrals, [ml,ij])
!
      use modbasis, only : nao
      implicit none
      integer,intent(in) :: noac, mlcount, mlstart, mlsize, nproc, myrank, idis(0:nproc-1,4)
      integer :: moi, moij
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: trint1b(nao,noac,mlsize), cmoocc(nao,noac)
      real(8),intent(inout) :: trint2(idis(myrank,3),noac*(noac+1)/2)
!
      do moi= 1,noac
        moij=moi*(moi-1)/2+1
        call dgemm('T','N',mlcount,moi,nao,one,trint1b(1,moi,1),nao*noac,cmoocc,nao,zero,&
&                  trint2(mlstart,moij),idis(myrank,3))
      enddo
      return
end


!------------------------------------------------------------------------
  subroutine mp2int_sendrecv(trint2,trint3,trint4,icycle,irecv,numij, &
&                            idis,nproc,myrank,mpi_comm)
!------------------------------------------------------------------------
!
! Send and Receive second transformed integrals (mi|lj)
!
! In  : trint2  (Second transformed integrals)
!       icycle  (Mp2trans2 cycle number)
!       irecv   (Numbers of receiving data)
!       numij   (Number of active occupied MOs)
!       idis    (Information for parallelization)
! Out : trint4  (Second transformed integrals, [l,m])
! Work: trint3  (Receiving data)
!
      use modbasis, only : nao, nshell, mbf, locbf
      implicit none
      integer,intent(in) :: icycle, nproc, myrank, mpi_comm, irecv(0:nproc-1)
      integer,intent(in) :: numij, idis(0:nproc-1,4)
      integer :: iproc, jproc, ijstart, myij, ij, nsend, nrecv, ish, ksh, nbfi, nbfk
      integer :: locbfi, locbfk, i, k, ik, num, ii, jcount
      real(8),intent(in) :: trint2(idis(myrank,3),numij)
      real(8),intent(out) :: trint3(nao*nao), trint4(nao,nao)
!
      jproc= myrank
      ijstart=(icycle-1)*nproc+1
      myij= ijstart+myrank
!
! Send and receive second transformed integrals
!
      do iproc= myrank+1,nproc-1
        jproc= jproc-1
        if(jproc < 0) jproc= nproc-1
        nsend= idis(myrank,3)
        nrecv= idis(jproc,3)
        ij= ijstart+iproc
        if(ij > numij) then
          nsend= 0
          ij   = 1
        endif
        if(myij > numij) nrecv=0
        call para_sendrecvr(trint2(1,ij),nsend,iproc,myrank, &
&                           trint3(irecv(jproc)),nrecv,jproc,jproc,mpi_comm)
      enddo
!
      do iproc= 0, myrank-1
        jproc= jproc-1
        if(jproc < 0) jproc= nproc-1
        nsend= idis(myrank,3)
        nrecv= idis(jproc,3)
        ij= ijstart+iproc
        if(ij > numij) then
          nsend= 0
          ij   = 1
        endif
        if(myij > numij) nrecv=0
        call para_sendrecvr(trint2(1,ij),nsend,iproc,myrank, &
&                           trint3(irecv(jproc)),nrecv,jproc,jproc,mpi_comm)
      enddo
!
      nsend= idis(myrank,3)
      if(myij <= numij) call dcopy(nsend,trint2(1,myij),1,trint3(irecv(myrank)),1)
!
! Reorder of received data
!
      if(myij <= numij) then
        ik= 0
        do iproc= 0,nproc-1
          ish= idis(iproc,1)
          ksh= idis(iproc,2)
          jcount= 0
          do num= 1,idis(iproc,4)
            nbfi= mbf(ish)
            nbfk= mbf(ksh)
            locbfi= locbf(ish)
            locbfk= locbf(ksh)
            do i= 1,nbfi
              do k= 1,nbfk
                ik= ik+1
                trint4(locbfk+k,locbfi+i)= trint3(ik)
              enddo
            enddo
!
            if(jcount == iproc) then
              ksh= ksh+2*nproc-1
            else
              ksh= ksh+nproc-1
            endif
            if(ksh > nshell) then
              do ii= 1,nshell
                ksh= ksh-nshell
                ish= ish+1
                if(ksh <= nshell) exit
              enddo
            endif
            jcount= jcount+1
            if(jcount == nproc) jcount= 0
          enddo
        enddo
      endif
      return
end


!----------------------------------------------------------------------------------------
  subroutine calcrmp2energy(trint4,energymo,emp2st,noac,nvac,ncore,icycle,nproc,myrank)
!----------------------------------------------------------------------------------------
!
! Calculate MP2 energy
!
! In  : trint4   (MO integrals, [b,a])
!       energymo (MO energies)
!       noac     (Number of active occupied Mos)
!       nvac     (Number of virtual MOs)
!       icycle   (Mp2trans2 cycle number)
!
      use modmolecule, only : neleca, nmo
      implicit none
      integer,intent(in) :: noac, nvac, ncore, icycle, nproc, myrank
      integer :: moi, moj, myij, ii, moa, mob
      real(8),parameter:: zero=0.0D+00, one=1.0D+00, two=2.0D+00
      real(8),intent(in) :: trint4(nvac,nvac), energymo(nmo)
      real(8),intent(inout) :: emp2st(2)
      real(8) :: eij, eija, eijab, etmp(2)
!
! Calculate occupied MOs, i and j
!
      myij=(icycle-1)*nproc+1+myrank
      do ii= 1,noac
        if(myij <=(ii*(ii+1)/2)) then
          moi= ii
          moj= myij-ii*(ii-1)/2
          exit
        endif
      enddo
!
! Calculate MP2 energy
!
      eij= energymo(moi+ncore)+energymo(moj+ncore)
      etmp(:)= zero
!$OMP parallel do private(eija,eijab) reduction(+:etmp)
      do moa= 1,nvac
        eija= eij-energymo(moa+neleca)
        do mob= 1,nvac
          eijab= one/(eija-energymo(mob+neleca))
          etmp(1)= etmp(1)+trint4(mob,moa)*trint4(mob,moa)*eijab
          etmp(2)= etmp(2)+trint4(mob,moa)*(trint4(mob,moa)-trint4(moa,mob))*eijab
        enddo
      enddo
!$OMP end parallel do
!
      if(moi /= moj) then
        emp2st(1)= emp2st(1)+etmp(1)*two
        emp2st(2)= emp2st(2)+etmp(2)*two
      else
        emp2st(1)= emp2st(1)+etmp(1)
        emp2st(2)= emp2st(2)+etmp(2)
      endif
      return
end


!-----------------------------------------------------------------------------
  subroutine mp2trans12m(cmoocc,cmowrk,trint1a,trint1b,trint2,xint, &
&                        mlsize,noac,maxdim,idis,numij,ijindex,nproc,myrank)
!-----------------------------------------------------------------------------
!
! Driver of AO intengral generation and first and second integral transformations
! for multiple pass
!
! In  : cmoocc  (MO coefficient matrix)
!       xint    (Exchange integral matrix)
!       mlsize  (Block size of first transformed integrals)
!       noac    (Number of active occupied MOs)
!       maxdim  (Maximum size of basis functions in a shell)
!       idis    (Information for parallelization)
!       numij   (Number of active occupied MO pairs)
!       ijindex (First and last indices of active occupied MO pairs)
! Out : trint2  (Second transformed integrals)
! Work: cmowrk  (Transposed MO coefficient matrix)
!       trint1a (First transformed integrals)
!       trint1b (First transformed integrals)
!       
      use modbasis, only : nshell, nao, mbf
      implicit none
      integer,intent(in) :: maxdim, mlsize, noac, nproc, myrank, idis(0:nproc-1,4)
      integer,intent(in) :: numij, ijindex(4)
      integer :: ish, ksh, ish1, ksh1, mlcount, mlstart, mlshell, numshell, ii, jj, jcount
      integer :: jcount1, numi
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: cmoocc(nao,noac), xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: cmowrk(ijindex(3)-ijindex(1)+1,nao)
      real(8),intent(out) :: trint1a(noac,maxdim,maxdim**2)
      real(8),intent(out) :: trint1b(mlsize*noac*nao)
      real(8),intent(out) :: trint2(idis(myrank,3)*numij)
!
      numi= ijindex(3)-ijindex(1)+1
!$OMP parallel do
      do ii= 1,numi
        do jj= 1,nao
          cmowrk(ii,jj)= cmoocc(jj,ijindex(1)+ii-1)
        enddo
      enddo
!$OMP end parallel do
!
      ish= idis(myrank,1)
      ksh= idis(myrank,2)
      ish1= ish
      ksh1= ksh
      jcount= 0
      jcount1= 0
      mlcount= 0
      mlstart= 1
      mlshell=0
      do numshell= 1,idis(myrank,4)
        mlshell= mlshell+1
        mlcount= mlcount+mbf(ish)*mbf(ksh)
        if(numshell == idis(myrank,4)) then
!
! AO intengral generation and first integral transformation
!
          call transmoint1(trint1a,trint1b,cmowrk,xint,ish1,ksh1,maxdim,numi,jcount1, &
&                          mlshell,mlsize,nproc,myrank)
!
! Second integral transformation
!
          call transmoint2m(trint2,trint1b,cmoocc,noac,mlcount,mlstart,mlsize,idis, &
&                           numij,ijindex,nproc,myrank)
        else
          if(jcount == myrank) then
            ksh= ksh+2*nproc-1
          else
            ksh= ksh+nproc-1
          endif
          jcount= jcount+1
          if(jcount == nproc) jcount= 0
!
          if(ksh > nshell) then
            do ii= 1,nshell
              ksh= ksh-nshell
              ish= ish+1
              if(ksh <= nshell) exit
            enddo
          endif
          if(mlcount+mbf(ish)*mbf(ksh) > mlsize) then
            call transmoint1(trint1a,trint1b,cmowrk,xint,ish1,ksh1,maxdim,numi,jcount1, &
&                            mlshell,mlsize,nproc,myrank)
            call transmoint2m(trint2,trint1b,cmoocc,noac,mlcount,mlstart,mlsize,idis, &
&                             numij,ijindex,nproc,myrank)
            mlstart= mlstart+mlcount
            mlshell= 0
            mlcount= 0
            jcount1= jcount
            ish1= ish
            ksh1= ksh
          endif
        endif
      enddo
!
!     if(iprint >= 3) then
!       write(*,'("Time for Tr1=",f9.3,"  Time for Tr2=",f9.3,i3)')tr1,tr2,myrank
!     endif
      return
end


!---------------------------------------------------------------------------------------
  subroutine mp2trans34m(cmovir,energymo,trint2,trint3,trint4,emp2st,nvac,ncore, &
&                        idis,numij,ijindex,nproc,myrank,mpi_comm)
!---------------------------------------------------------------------------------------
!
! Driver of third and fourth integral transformations and MP2 energy calculation
!
! In  : cmovir  (MO coefficient matrix)
!       energymo(MO energies)
!       trint2  (Second transformed integrals)
!       nvac    (Number of virtual MOs)
!       idis    (Information for parallelization)
!       numij   (Number of active occupied MO pairs)
!       ijindex (First and last indices of active occupied MO pairs)
! Work: trint3  (Third transformed integrals)
!       trint4  (Fourth transformed integrals)
!
      use modbasis, only : nao
      use modmolecule, only : nmo
      implicit none
      integer,intent(in) :: nvac, ncore, nproc, myrank, mpi_comm, idis(0:nproc-1,4)
      integer,intent(in) :: numij, ijindex(4)
      integer :: numrecv, iproc, irecv(0:nproc-1), ncycle, icycle, myij
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: cmovir(nao*nvac), energymo(nmo)
      real(8),intent(in) :: trint2(idis(myrank,3),numij)
      real(8),intent(out) :: trint3(nao*nao), trint4(nao*nao), emp2st(2)
!
      numrecv= 1
      do iproc= 0,nproc-1
        irecv(iproc)= numrecv
        numrecv= numrecv+idis(iproc,3)
      enddo
      ncycle=(numij-1)/nproc+1
!
      do icycle= 1,ncycle    
        call mp2int_sendrecv(trint2,trint3,trint4,icycle,irecv,numij, &
&                            idis,nproc,myrank,mpi_comm)
!
        myij=(icycle-1)*nproc+1+myrank
        if(myij > numij) cycle
!
! Third integral transformation
!   (mi|lj) trint4[l,m] -> (ai|lj) trint3[l,a]
!
        call dgemm('N','N',nao,nvac,nao,one,trint4,nao,cmovir,nao,zero,trint3,nao)
!
! Fourth integral transformation
!   (ai|lj) trint3[l,a] -> (ai|bj) trint4[b,a]
!
        call dgemm('T','N',nvac,nvac,nao,one,cmovir,nao,trint3,nao,zero,trint4,nvac)
!
! MP2 energy calculation
!
        call calcrmp2energym(trint4,energymo,emp2st,nvac,ncore,icycle,ijindex,nproc,myrank)
      enddo
      return
end


!------------------------------------------------------------------------------------
  subroutine transmoint2m(trint2,trint1b,cmoocc,noac,mlcount,mlstart,mlsize,idis, &
&                         numij,ijindex,nproc,myrank)
!------------------------------------------------------------------------------------
!
! Second-quarter integral transformation for multiple pass
!    (mi|ls) -> (mi|lj)
!
! In  : trint1b (First transformed integrals, [s,i,ml])
!       cmoocc  (MO coefficient matrix)
!       noac    (Number of active occupied MOs)
!       mlcount (Number of transformed AOs)
!       mlstart (First index of trint2)
!       mlsize  (Size of trint1b)
!       idis    (Information for parallelization)
!       numij   (Number of active occupied MO pairs)
!       ijindex (First and last indices of active occupied MO pairs)
! Out : trint2  (Second transformed integrals, [ml,ij])
!
      use modbasis, only : nao
      implicit none
      integer,intent(in) :: noac, mlcount, mlstart, mlsize, nproc, myrank, idis(0:nproc-1,4)
      integer,intent(in) :: numij, ijindex(4)
      integer :: numi, moij, moi, mojf, numj
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: trint1b(nao,ijindex(3)-ijindex(1)+1,mlsize), cmoocc(nao,noac)
      real(8),intent(inout) :: trint2(idis(myrank,3),numij)
!
      numi= ijindex(3)-ijindex(1)+1
      moij= 1
      do moi= ijindex(1),ijindex(3)
        mojf= 1
        numj= moi
        if(moi == ijindex(1)) then
          mojf= ijindex(2)
          numj= moi-mojf+1
        endif
        if(moi == ijindex(3)) then
          numj= ijindex(4)-mojf+1
        endif
        call dgemm('T','N',mlcount,numj,nao,one,trint1b(1,moi-ijindex(1)+1,1),nao*numi, &
&                  cmoocc(1,mojf),nao,zero,trint2(mlstart,moij),idis(myrank,3))
        moij= moij+numj
      enddo
!
      return
end


!--------------------------------------------------------------------------------------------
  subroutine calcrmp2energym(trint4,energymo,emp2st,nvac,ncore,icycle,ijindex,nproc,myrank)
!--------------------------------------------------------------------------------------------
!
! Calculate MP2 energy
!
! In  : trint4   (MO integrals, [b,a])
!       energymo (MO energies)
!       nvac     (Number of virtual MOs)
!       icycle   (Mp2trans2 cycle number)
!       ijindex  (First and last indices of active occupied MO pairs)
! Out : emp2st   (Partial MP2 energy)
!
      use modmolecule, only : neleca, nmo
      implicit none
      integer,intent(in) :: nvac, ncore, icycle, ijindex(4), nproc, myrank
      integer :: moi, moj, ii, moa, mob
      real(8),parameter:: zero=0.0D+00, one=1.0D+00, two=2.0D+00
      real(8),intent(in) :: trint4(nvac,nvac), energymo(nmo)
      real(8),intent(inout) :: emp2st(2)
      real(8) :: eij, eija, eijab, etmp(2)
!
! Calculate occupied MOs i and j
!
      moi= ijindex(1)
      moj=(icycle-1)*nproc+1+myrank+ijindex(2)-1
      do ii= ijindex(1),ijindex(3)
        if(moi < moj) then
          moj= moj-moi
          moi= moi+1
        endif
      enddo
!
! Calculate MP2 energy
!
      eij= energymo(moi+ncore)+energymo(moj+ncore)
      etmp(:)= zero
!$OMP parallel do private(eija,eijab) reduction(+:etmp)
      do moa= 1,nvac
        eija= eij-energymo(moa+neleca)
        do mob= 1,nvac
          eijab= one/(eija-energymo(mob+neleca))
          etmp(1)= etmp(1)+trint4(mob,moa)*trint4(mob,moa)*eijab
          etmp(2)= etmp(2)+trint4(mob,moa)*(trint4(mob,moa)-trint4(moa,mob))*eijab
        enddo
      enddo
!$OMP end parallel do
!
      if(moi /= moj) then
        emp2st(1)= emp2st(1)+etmp(1)*two
        emp2st(2)= emp2st(2)+etmp(2)*two
      else
        emp2st(1)= emp2st(1)+etmp(1)
        emp2st(2)= emp2st(2)+etmp(2)
      endif
      return
end


!------------------------------------------------------------------
  subroutine mp2int_isendrecv(trint2,recvint,icycle,irecv,numij, &
&                             idis,ireq,nproc,myrank,mpi_comm)
!------------------------------------------------------------------
!
! Non-blocking communications of second transformed integrals (mi|lj)
!
! In  : trint2  (Second transformed integrals)
!       icycle  (Mp2trans2 cycle number)
!       irecv   (Numbers of receiving data)
!       numij   (Number of active occupied MOs)
!       idis    (Information for parallelization)
! Out : recvint (Receiving data)
!       ireq    (Request for MPI_Wait)
      use modbasis, only : nao
      implicit none
      integer,intent(in) :: icycle, nproc, myrank, mpi_comm, irecv(0:nproc-1)
      integer,intent(in) :: numij, idis(0:nproc-1,4)
      integer,intent(out) :: ireq(0:nproc-1,2)
      integer :: iproc, jproc, ijstart, myij, ij, nsend, nrecv
      real(8),intent(in) :: trint2(idis(myrank,3),numij)
      real(8),intent(out) :: recvint(nao*nao)
!
      jproc= myrank
      ijstart=(icycle-1)*nproc+1
      myij= ijstart+myrank
!
! Send and receive second transformed integrals
!
      if(nproc > 1) then
        do iproc= myrank+1,nproc-1
          jproc= jproc-1
          if(jproc < 0) jproc= nproc-1
          nsend= idis(myrank,3)
          nrecv= idis(jproc,3)
          ij= ijstart+iproc
          if(ij > numij) then
            nsend= 0
            ij   = 1
          endif
          if(myij > numij) nrecv=0
          call para_isendr(trint2(1,ij),nsend,iproc,myrank,mpi_comm,ireq(iproc,1))
          call para_irecvr(recvint(irecv(jproc)),nrecv,jproc,jproc,mpi_comm,ireq(iproc,2))
        enddo
!
        do iproc= 0, myrank-1
          jproc= jproc-1
          if(jproc < 0) jproc= nproc-1
          nsend= idis(myrank,3)
          nrecv= idis(jproc,3)
          ij= ijstart+iproc
          if(ij > numij) then
            nsend= 0
            ij   = 1
          endif
          if(myij > numij) nrecv=0
          call para_isendr(trint2(1,ij),nsend,iproc,myrank,mpi_comm,ireq(iproc,1))
          call para_irecvr(recvint(irecv(jproc)),nrecv,jproc,jproc,mpi_comm,ireq(iproc,2))
        enddo
!  
        nsend= idis(myrank,3)
        if(myij > numij) nsend= 0
        call para_isendr(trint2(1,myij),nsend,myrank,myrank,mpi_comm,ireq(myrank,1))
        call para_irecvr(recvint(irecv(myrank)),nsend,myrank,myrank,mpi_comm,ireq(myrank,2))
      else
        nsend= idis(myrank,3)
        if(myij <= numij) call dcopy(nsend,trint2(1,myij),1,recvint(irecv(myrank)),1)
      endif
!
      return
end


!--------------------------------------------------------------------------------------
  subroutine mp2int_sort(recvint,trint4,icycle,numij,idis,ireq,nproc,myrank,mpi_comm)
!--------------------------------------------------------------------------------------
!
! MPI_Waitall and sorting of receieved second-transformed integrals
!
      use modbasis, only : nao, nshell, mbf, locbf
      implicit none
      integer,intent(in) :: icycle, nproc, myrank, mpi_comm, numij
      integer,intent(in) :: idis(0:nproc-1,4), ireq(2*nproc)
      integer :: iproc, ijstart, myij, ish, ksh, nbfi, nbfk
      integer :: locbfi, locbfk, ii, kk, ik, num, jcount
      real(8),intent(in) :: recvint(nao*nao)
      real(8),intent(out) :: trint4(nao,nao)
!
      ijstart=(icycle-1)*nproc+1
      myij= ijstart+myrank
!
      if(nproc > 1) call para_waitall(2*nproc,ireq)
!
! Reorder of received data
!
      if(myij <= numij) then
        ik= 0
        do iproc= 0,nproc-1
          ish= idis(iproc,1)
          ksh= idis(iproc,2)
          jcount= 0
          do num= 1,idis(iproc,4)
            nbfi= mbf(ish)
            nbfk= mbf(ksh)
            locbfi= locbf(ish)
            locbfk= locbf(ksh)
            do ii= 1,nbfi
              do kk= 1,nbfk
                ik= ik+1
                trint4(locbfk+kk,locbfi+ii)= recvint(ik)
              enddo
            enddo
!
            if(jcount == iproc) then
              ksh= ksh+2*nproc-1
            else
              ksh= ksh+nproc-1
            endif
            if(ksh > nshell) then
              do ii= 1,nshell
                ksh= ksh-nshell
                ish= ish+1
                if(ksh <= nshell) exit
              enddo
            endif
            jcount= jcount+1
            if(jcount == nproc) jcount= 0
          enddo
        enddo
      endif
!
      return
end

