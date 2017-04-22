! Copyright 2016-2017  Kazuya Ishimura
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
!-------------------------------------------------------------------------
  subroutine calcgradrmp2(cmo,energymo,xint,egrad,nproc,myrank,mpi_comm)
!-------------------------------------------------------------------------
!
! Main driver of closed-shell MP2 energy gradient calculation
!
! In  : cmo     (MO coefficient matrix)
!       energymo(MO energies)
!       xint    (Exchange integral matrix)
! Out : egrad   (MP2 energy gradients)
!
      use modparallel, only : master
      use modbasis, only : nao, nshell, mbf, mtype
      use modmolecule, only : nmo, natom, neleca, numatomic
      use modenergy, only : escf, emp2, escsmp2
      use modmp2, only : ncore, nvfz, maxmp2diis
      implicit none
      integer,intent(in) :: nproc, myrank, mpi_comm
      integer :: maxfunc(0:7), maxdim, maxgraddim, nocc, nvir, noac, nvac, icount, ish, jsh
      integer :: idis(0:nproc-1,8), iproc, nvac2, numab, numov, maxsize, numirecv
      integer :: nocc2, nocc3, nvir2, nao2, nao3, mem1, mem2, mem2i, mem3, mem3i, mem4, mem5
      integer :: mem5i, mem6, memmin345, memmin16, msize, memmin3, memmin4, memmin5
      integer :: numi3, numi4, numi5, numi, npass, ii, jj
      real(8),parameter :: zero=0.0D+00, three=3.0D+00, p12=1.2D+00
      real(8),intent(in) :: cmo(nao,nao), energymo(nao), xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: egrad(3,natom)
      real(8) :: emp2st(2), emp2stsum(2), egradtmp(3*natom)
      character(len=3) :: table(-9:112)= &
&     (/'Bq9','Bq8','Bq7','Bq6','Bq5','Bq4','Bq3','Bq2','Bq ','X  ',&
&       'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ','Al ','Si ','P  ',&
&       'S  ','Cl ','Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',&
&       'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ',&
&       'Pd ','Ag ','Cd ','In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',&
&       'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ','Ta ','W  ','Re ',&
&       'Os ','Ir ','Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',&
&       'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ','Md ','No ','Lr ','Rf ','Db ',&
&       'Sg ','Bh ','Hs ','Mt ','Uun','Uuu','Uub'/)
      data maxfunc/1,3,6,10,15,21,28,36/
!
      emp2= zero
      emp2st(:)= zero
      egrad(:,:)= zero
      egradtmp(:)= zero
      maxdim= maxval(mtype(1:nshell))
      maxgraddim= maxfunc(maxdim+1)
      maxdim= maxfunc(maxdim)
      nocc= neleca
      nvir= nmo-neleca
      noac= nocc-ncore
      nvac= nvir-nvfz
!
      if(master) then
        write(*,'(" ----------------------------------------------")')
        write(*,'("   MP2 energy gradient calculation ")')
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
!
      nvac2= nvac*nvac
      numab=(nvac2-1)/nproc+1
      numov=(nocc*nvir-1)/nproc+1
      do iproc= 0,nproc-1
        if(numab*(iproc+1) <= nvac2) then
          idis(iproc,5)= numab
          idis(iproc,6)= numab*iproc
        elseif(numab*iproc < nvac2) then
          idis(iproc,5)= nvac2-numab*iproc
          idis(iproc,6)= numab*iproc
        else
          idis(iproc,5)= 0
          idis(iproc,6)= 0
        endif
        if(numov*(iproc+1) <= nocc*nvir) then
          idis(iproc,7)= numov
          idis(iproc,8)= numov*iproc
        elseif(numov*iproc < nocc*nvir) then
          idis(iproc,7)= nocc*nvir-numov*iproc
          idis(iproc,8)= numov*iproc
        else
          idis(iproc,7)= 0
          idis(iproc,8)= 0
        endif
      enddo
!
! Check available memory size and decide the number of passes
!
      maxsize= maxval(idis(0:nproc-1,3))
      numirecv= nproc/noac
      if(mod(nproc,noac) /= 0) numirecv= numirecv+2
!
      nocc2= nocc*nocc
      nocc3= nocc*(nocc+1)/2
      nvir2= nvir*nvir
      nao2= nao*nao
      nao3= nao*(nao+1)/2
      mem1= nocc2+nvir2+nocc*nvac+nocc*nvir+nao2
      mem2= nocc3+nvir2
      mem2i=maxsize*nocc
      mem3= nao*noac
      mem3i= maxdim**3+nao*maxdim**2
      mem4= 2*nproc*maxsize+nao2+nocc*nvac+2*numab*noac*numirecv+2*numab
      mem5= 2*nao2+nao*maxdim**2
      mem5i= nao*maxdim**2+nao
      mem6= 2*nocc*nvir+2*nao3+nshell*(nshell+1)/2+2*numov*maxmp2diis &
&          +maxmp2diis*(maxmp2diis+1)/2+2*nao2
!
      memmin345=max(mem3+mem3i,mem4,mem5+mem5i)
      memmin16= max(mem1+mem2+mem2i+memmin345,mem1+mem6)
!
      call memrest(msize)
      if((msize < memmin16).and.master) then
        write(*,'(" Error! Available memory size for MP2 energy gradient is small!")')
        call memset(memmin16)
        call iabort
      endif
!
      memmin3= msize-(mem1+mem2+mem3)
      memmin4= msize-(mem1+mem2+mem4)
      memmin5= msize-(mem1+mem2+mem5)
      numi3= memmin3/(mem2i+mem3i)
      numi4= memmin4/(mem2i)
      numi5= memmin5/(mem2i+mem5i)
      numi= min(numi3,numi4,numi5)
      if(numi > noac) numi= noac
!
      if(numi <= 0) then
        if(master) then
          write(*,'(" Error! Available memory size for MP2 energy gradient is small!")')
          call iabort
        endif
      else
        npass=(noac-1)/numi+1
        if(master) then
          if(npass == 1) then
            write(*,'(" == Single pass calculation ==")')
          else
            write(*,'(" == Multiple pass calculation ==")')
            write(*,'("    Number of passes :",i5)')npass
          endif
        endif
        call mp2gradmulti(emp2st,egradtmp,egrad,cmo,energymo,xint,nocc,noac,nvir,nvac, &
&                         ncore,nvfz,maxsize,maxdim,maxgraddim,idis,npass,numi,numab,numirecv, &
&                         nproc,myrank,mpi_comm)
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
!
      call para_allreducer(egradtmp,egrad,3*natom,mpi_comm)
!
      if(master) then
        write(*,'(" ----------------------------------------------------")')
        write(*,'("          Gradient (Hartree/Bohr)")')
        write(*,'("  Atom            X             Y             Z")')
        write(*,'(" ----------------------------------------------------")')
        do ii= 1,natom
          write(*,'(3x,a3,3x,3f14.7)')table(numatomic(ii)),(egrad(jj,ii),jj=1,3)
        enddo
        write(*,'(" ----------------------------------------------------")')
      endif
!
      return
end


!-------------------------------------------------------------------------------------------------
  subroutine mp2gradmulti(emp2st,egrad,egradtmp,cmo,energymo,xint,nocc,noac,nvir,nvac, &
&                         ncore,nvfz,maxsize,maxdim,maxgraddim,idis,npass,numi,numab,numirecv, &
&                         nproc,myrank,mpi_comm)
!-------------------------------------------------------------------------------------------------
!
! Driver of single pass MP2 energy calculation
!
! In  : cmo     (MO coefficient matrix)
!       energymo(MO energies)
!       xint    (Exchange integral matrix)
!       nocc    (Number of occupied MOs)
!       noac    (Number of active occupied MOs)
!       nvir    (Number of virtual MOs)
!       nvac    (Number of active virtual MOs)
!       ncore   (Number of core occupied MOs)
!       nvfz    (Number of frozen vitrual MOs)
!       maxsize (Maximum size of basis function pairs)
!       maxdim  (Maximum size of basis functions in a shell)
!       maxgraddim(Maximum size of basis functions in a shell for derivative calculation)
!       idis    (Information for parallelization)
!       npass   (Number of passes)
!       numi    (Number of active occupied MOs per pass)
!       numab   (Number of active virtual MO pairs per process)
!       numirecv(Number of active occupied MOs for receiving tijab+(ia|jb))
! Out : emp2st  (Partial MP2 energies)
!       egrad   (Partial MP2 energy gradients)
! Work: egradtmp(Work space for MP2 energy gradients)
!
      use modparallel, only : master
      use modbasis, only : nshell, nao
      use modmolecule, only : nmo, natom
      use modmp2, only : maxmp2diis
      implicit none
      integer,intent(in) :: nocc, noac, nvir, nvac, ncore, nvfz, maxsize, maxdim, maxgraddim
      integer,intent(in) :: nproc, myrank, mpi_comm, idis(0:nproc-1,8), npass
      integer,intent(in) :: numi, numab, numirecv
      integer :: nao2, nao3, nocc2, nocc3, nvir2, numitrans, ipass, msize, mlsize, istart
      integer :: mlsize2
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: cmo(nao,nao), energymo(nmo), xint(nshell*(nshell+1)/2)
      real(8),intent(inout) :: emp2st(2), egrad(3*natom), egradtmp(3*natom)
      real(8),allocatable :: wij(:), wab(:), wai(:), xlai(:), pmn(:)
      real(8),allocatable :: trint2(:), trint2core(:), pij(:), pab(:)
      real(8),allocatable :: cmowrk(:), trint1a(:), trint1b(:)
      real(8),allocatable :: xaibj(:), xaikj(:), tijab(:), recvint(:), sendint(:)
      real(8),allocatable :: recvt(:), storet(:)
      real(8),allocatable :: tisml(:), xlmi(:), xlmn(:), work1(:)
      real(8),allocatable :: pai(:), paiprev(:), pls(:), pmax(:), paifock(:)
      real(8),allocatable :: errdiis(:), paidiis(:), diismtrx(:), work2(:)
!
      nao2= nao*nao
      nao3= nao*(nao+1)/2
      nocc2= nocc*nocc
      nocc3= nocc*(nocc+1)/2
      nvir2= nvir*nvir
      numitrans= numi
! mem1
      call memset(nocc2+nvir2+nocc*nvac+nocc*nvir+nao2)
      allocate(wij(nocc2),wab(nvir2),wai(nocc*nvac),xlai(nocc*nvir),pmn(nao2))
! mem2
      call memset(maxsize*nocc*numi+nocc3+nvir2)
      allocate(trint2(maxsize*noac*numi),trint2core(maxsize*ncore*numi),pij(nocc3),pab(nvir2))
!
      call memrest(msize)
      mlsize=(msize-nao*noac-numi*maxdim**3)/(nao*numi)
      if(mlsize > maxsize) then
        mlsize= maxsize
      elseif(mlsize < maxdim*maxdim) then
        mlsize= maxdim*maxdim
      endif
!
      pij(:)= zero
      pab(:)= zero
      wij(:)= zero
      wab(:)= zero
      wai(:)= zero
      xlai(:)= zero
!
! Start multiple pass
!
      do ipass= 1,npass
        if(master) write(*,'("    Start Pass",i5)')ipass
        if(ipass == npass) numitrans= noac-(npass-1)*numi
        istart=(ipass-1)*numi
!
! AO intengral generation and first and second integral transformations
!
! mem3
        call memset(nao*noac+numi*maxdim**3+mlsize*nao*numi)
        allocate(cmowrk(nao*noac),trint1a(numi*maxdim**3),trint1b(mlsize*nao*numi))
!
        call mp2gradtrans12(cmo,cmowrk,trint1a,trint1b,trint2,trint2core,xint,istart,mlsize, &
&                           noac,ncore,maxdim,numitrans,idis,nproc,myrank)
!
        deallocate(cmowrk,trint1a,trint1b)
        call memunset(nao*noac+numi*maxdim**3+mlsize*nao*numi)
!ishimura
!call tstamp(1)
!
! Third and fourth integral transformations, and
! MP2 energy, Pij, Pab, Wij[I], Wab[I], Wai[I](=Lai3), Tijab half back-transformation calculations
!
! mem4
        call memset(2*nproc*maxsize+nao2+nocc*nvac+2*numab*noac*numirecv+2*numab)
        allocate(xaibj(nproc*maxsize),xaikj(nocc*nvac),tijab(nao2),recvint(2*numab*noac*numirecv), &
&                sendint(2*numab),recvt(nproc*maxsize),storet(noac*maxsize))
!
        call mp2gradtrans34(emp2st,cmo,energymo,trint2,trint2core,xaibj,xaikj, &
&                           tijab,pij,pab,wij,wab,wai,recvint,sendint,recvt,storet, &
&                           nocc,noac,nvir,nvac,ncore,nvfz,istart,numitrans,numirecv, &
&                           numab,maxsize,idis,nproc,myrank,mpi_comm)
!
        deallocate(xaibj,xaikj,tijab,recvint, &
&                  sendint,recvt,storet)
        call memunset(2*nproc*maxsize+nao2+nocc*nvac+2*numab*noac*numirecv+2*numab)
!ishimura
!call tstamp(1)
!
! Second-half back-transformation (tnsml), and
! Lai1,2,4, tnsml*(mn|ls)', PiJ, Wij[II], and Wab[II] calculations
!
! mem5
        call memrest(msize)
        mlsize2=(msize-nao*numi-2*nao2)/(nao*(numi+1))
        if(mlsize2 > maxsize) mlsize2= maxsize
        call memset(numi*mlsize2*nao+nao*numi+2*nao2+nao*mlsize2)
        allocate(tisml(numi*mlsize2*nao),xlmi(nao*numi),xlmn(nao2),work1(nao2),work2(nao*mlsize2))
!
        if(ipass /= npass) then
          call mp2gradbacktrans1(egrad,egradtmp,xlai,tisml,xlmi,cmo,xint,trint2,work1,work2, &
&                                nocc,noac,nvir,ncore,numitrans,istart,maxdim,maxgraddim, &
&                                mlsize2,idis,nproc,myrank)
        else
          call mp2gradbacktrans2(egrad,egradtmp,wij,wab,wai,xlai,tisml,xlmi,xlmn,cmo,xint,trint2, &
&                                pij,pab,pmn,energymo,work1,work2,nocc,noac,nvir,ncore,numitrans, &
&                                istart,maxdim,maxgraddim,mlsize2,idis,nproc,myrank,mpi_comm)
        endif
!
        deallocate(tisml,xlmi,xlmn,work1,work2)
        call memunset(numi*mlsize2*nao+nao*numi+2*nao2+nao*mlsize2)
!ishimura
!call tstamp(1)
      enddo
!
      deallocate(trint2,trint2core,pij,pab)
      call memunset(maxsize*nocc*numi+nocc3+nvir2)
! mem6
      call memset(2*nocc*nvir+2*nao3+nshell*(nshell+1)/2+2*idis(myrank,7)*maxmp2diis &
&                 +maxmp2diis*(maxmp2diis+1)/2+2*nao2)
      allocate(pai(nocc*nvir),paiprev(nocc*nvir),pls(nao3),pmax(nshell*(nshell+1)/2), &
&              paifock(nao3),errdiis(idis(myrank,7)*maxmp2diis), &
&              paidiis(idis(myrank,7)*maxmp2diis),diismtrx(maxmp2diis*(maxmp2diis+1)/2), &
&              work1(nao2),work2(nao2))
!
! Solve CPHF equation and obtain Pai
!
      call mp2gradcphf(pai,xlai,cmo,xint,energymo,paiprev,pls,pmax,paifock,errdiis,paidiis, &
&                      diismtrx,work1,work2,nocc,nvir,maxdim,idis,nproc,myrank,mpi_comm)
!ishimura
!call tstamp(1)
!
! Calculate Wij[III] and Wai[II]
!
      call mp2gradwij3(wij,wai,pmn,pai,paifock,cmo,energymo,xint,pls,pmax,work1,work2, &
&                      nocc,nvir,nvac,maxdim,nproc,myrank,mpi_comm)
!ishimura
!call tstamp(1)
!
! Calculate integral derivatives and their MP2 energy gradient
!
      call mp2graddint(egrad,egradtmp,pmn,wij,wab,wai,xint,cmo,energymo,pls,work1,work2, &
&                      nocc,nvir,maxdim,maxgraddim,nproc,myrank,mpi_comm)
!ishimura
!call tstamp(1)
!
      deallocate(pai,paiprev,pls,pmax, &
&                paifock,errdiis, &
&                paidiis,diismtrx, &
&                work1,work2)
      call memunset(2*nocc*nvir+2*nao3+nshell*(nshell+1)/2+2*idis(myrank,7)*maxmp2diis &
&                   +maxmp2diis*(maxmp2diis+1)/2+2*nao2)
!
      deallocate(wij,wab,wai,xlai,pmn)
      call memunset(nocc2+nvir2+nocc*nvac+nocc*nvir+nao2)
!
      return
end


!-----------------------------------------------------------------------------------------------
  subroutine mp2gradtrans12(cmo,cmowrk,trint1a,trint1b,trint2,trint2core,xint,istart,mlsize, &
&                           noac,ncore,maxdim,numitrans,idis,nproc,myrank)
!-----------------------------------------------------------------------------------------------
!
! Driver of AO intengral generation and first and second integral transformations
! for multiple pass of MP2 energy gradient calculation
!
! In  : cmo      (MO coefficient matrix)
!       xint     (Exchange integral matrix)
!       istart   (First index of transformed occupied MOs)
!       mlsize   (Block size of first transformed integrals)
!       noac     (Number of active occupied MOs)
!       ncore    (Number of core occupied MOs)
!       maxdim   (Maximum size of basis functions in a shell)
!       numitrans(Number of transformed active occupied MOs)
!       idis     (Information for parallelization)
! Out : trint2   (Second transformed integrals of active occupied MOs)
!       trint2core(Second transformed integrals of core occupied MOs)
! Work: cmowrk   (Transposed MO coefficient matrix)
!       trint1a  (First transformed integrals)
!       trint1b  (First transformed integrals)
!
      use modbasis, only : nshell, nao, mbf
      implicit none
      integer,intent(in) :: istart, mlsize, noac, ncore, maxdim, numitrans
      integer,intent(in) :: nproc, myrank, idis(0:nproc-1,8)
      integer :: ii, jj, ish, ksh, ish1, ksh1, jcount, jcount1, mlcount, mlstart, mlshell
      integer :: numshell
      real(8),intent(in) :: cmo(nao,nao), xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: cmowrk(numitrans,nao)
      real(8),intent(out) :: trint1a(numitrans,maxdim,maxdim**2)
      real(8),intent(out) :: trint1b(mlsize*numitrans*nao)
      real(8),intent(out) :: trint2(idis(myrank,3)*noac*numitrans)
      real(8),intent(out) :: trint2core(idis(myrank,3)*ncore*numitrans)
!
!$OMP parallel do
      do jj= 1,numitrans
        do ii= 1,nao
          cmowrk(jj,ii)= cmo(ii,jj+istart+ncore)
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
! AO integral generation and first integral transformation
!
          call transmoint1(trint1a,trint1b,cmowrk,xint,ish1,ksh1,maxdim,numitrans,jcount1, &
&                          mlshell,mlsize,nproc,myrank)
!
! Second integral transformation
!
          call transmointgrad2(trint2,trint2core,trint1b,cmo,noac,ncore,numitrans,mlcount, &
&                              mlstart,mlsize,idis,nproc,myrank)
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
            call transmoint1(trint1a,trint1b,cmowrk,xint,ish1,ksh1,maxdim,numitrans,jcount1, &
&                            mlshell,mlsize,nproc,myrank)
            call transmointgrad2(trint2,trint2core,trint1b,cmo,noac,ncore,numitrans,mlcount, &
&                                mlstart,mlsize,idis,nproc,myrank)
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
      return
end


!-------------------------------------------------------------------------------------------
  subroutine transmointgrad2(trint2,trint2core,trint1b,cmo,noac,ncore,numitrans,mlcount, &
&                            mlstart,mlsize,idis,nproc,myrank)
!-------------------------------------------------------------------------------------------
!
! Second-quarter integral transformation for MP2 energy gradient calculation
!    (mi|ls) -> (mi|lj)
!
! In  : trint1b  (First transformed integrals, [s,i,ml])
!       cmo      (MO coefficient matrix)
!       numitrans(Number of transformed occupied MOs)
!       mlcount  (Number of transformed AOs)
!       mlstart  (First index of trint2)
!       mlsize   (Size of trint1b)
!       ncore    (Number of core occupied MOs)
!       idis     (Information for parallelization)
! Out : trint2   (Second transformed integrals of active MOs, [ml,ij])
!       trint2core(Second transformed integrals of core MOs, [ml,ij])
!
      use modbasis, only : nao
      implicit none
      integer,intent(in) :: noac, ncore, numitrans, mlcount, mlstart, mlsize, nproc, myrank
      integer,intent(in) :: idis(0:nproc-1,8)
      integer :: moi
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: trint1b(nao,numitrans,mlsize), cmo(nao,nao)
      real(8),intent(inout) :: trint2(idis(myrank,3),noac,numitrans)
      real(8),intent(inout) :: trint2core(idis(myrank,3),ncore,numitrans)
!
      do moi= 1,numitrans
        call dgemm('T','N',mlcount,noac,nao,one,trint1b(1,moi,1),nao*numitrans,cmo(1,ncore+1), &
&                  nao,zero,trint2(mlstart,1,moi),idis(myrank,3))
      enddo
!
      if(ncore /= 0) then
        do moi= 1,numitrans
          call dgemm('T','N',mlcount,ncore,nao,one,trint1b(1,moi,1),nao*numitrans,cmo, &
&                    nao,zero,trint2core(mlstart,1,moi),idis(myrank,3))
        enddo
      endif
!
      return
end

!--------------------------------------------------------------------------------------------
  subroutine mp2gradtrans34(emp2st,cmo,energymo,trint2,trint2core,xaibj,xaikj, &
&                           tijab,pij,pab,wij,wab,wai,recvint,sendint,recvt,storet, &
&                           nocc,noac,nvir,nvac,ncore,nvfz,istart,numitrans,numirecv, &
&                           numab,maxsize,idis,nproc,myrank,mpi_comm)
!--------------------------------------------------------------------------------------------
!
! Driver of third and fourth integral transformations, and
! MP2 energy, Pij, Pab, Wij[I], Wab[I], Wai[I](=Lai3), and Tijab back-transformation calculations
!
! In  : cmo      (MO coefficient matrix)
!       energymo (MO energies)
!       trint2core(Second transformed integrals of core MOs)
!       nocc     (Number of occupied MOs)
!       noac     (Number of active occupied MOs)
!       nvir     (Number of virtual MOs)
!       nvac     (Number of active virtual MOs)
!       ncore    (Number of frozen core MOs)
!       nvfz     (Number of frozen virtual MOs)
!       istart   (First index of transformed occupied MOs)
!       numitrans(Number of transformed active occupied MOs)
!       numirecv (Buffer size of receiving xaibj and trint2)
!       numab    (Number of virtual MO pairs of receiving xaibj and trint2)
!       maxsize  (Maximum size of basis function pairs)
!       idis     (Information for parallelization)
! Inout:emp2st   (Partial MP2 energies)
!       trint2   (In: Second transformed integrals)
!                (Work: half back-transformed MP2 amplitudes)
!       pij      (Pij)
!       pab      (Pab)
!       wij      (Wij)
!       wab      (Wab)
!       wai      (Wai)
! Work: xaibj    ((mi|lj), (ai|bj), and (ai|bj)/Dijab)
!       xaikj    ((ai|kj))
!       tijab    ((ai|lj), tijab and tijml)
!       recvint  (Receiving buffer for tijab and (ai|bj)/Dijab)
!       sendint  (Sending buffer for tijab and (ai|bj)/Dijab)
!       recvt    (Receiving buffer for tijml)
!       storet   (Storing buffer for tijml)
!
      use modbasis, only : nao
      use modmolecule, only : nmo
      implicit none
      integer,intent(in) :: nocc, noac, nvir, nvac, ncore, nvfz, istart, numitrans, numirecv
      integer,intent(in) :: numab, maxsize, nproc, myrank, mpi_comm, idis(0:nproc-1,8)
      integer :: numij, numrecv, iproc, irecv(0:nproc-1), ncycle, icycle, myij, moi, moj
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: cmo(nao,nao), energymo(nmo)
      real(8),intent(out) :: xaibj(nproc*maxsize), xaikj(nocc*nvac), tijab(nao*nao)
      real(8),intent(out) :: recvint(2*numab*noac*numirecv), sendint(2*numab)
      real(8),intent(out) :: recvt(nproc*maxsize), storet(noac*maxsize)
      real(8),intent(inout) :: emp2st(2), trint2(idis(myrank,3),noac,numitrans)
      real(8),intent(inout) :: trint2core(idis(myrank,3),ncore,numitrans)
      real(8),intent(inout) :: pij(nocc*(nocc+1)/2), pab(nvir*nvir), wij(nocc,nocc)
      real(8),intent(inout) :: wab(nvir*nvir), wai(nocc*nvac)
!
      numij= noac*numitrans
      numrecv= 1
      do iproc= 0,nproc-1
        irecv(iproc)= numrecv
        numrecv= numrecv+idis(iproc,3)
      enddo
      ncycle=(numij-1)/nproc+1
!
      do icycle= 1,ncycle
        call mp2int_sendrecv(trint2,tijab,xaibj,icycle,irecv,numij, &
&                            idis,nproc,myrank,mpi_comm)
!
        myij=(icycle-1)*nproc+1+myrank
        if(myij <= numij) then
          moi=(myij-1)/noac+1+istart
          moj= mod(myij-1,noac)+1
!
! Third integral transformation
!   (mi|lj) xaibj[l,m] -> (ai|lj) tijab[l,a]
!
          call dgemm('N','N',nao,nvac,nao,one,xaibj,nao,cmo(1,nocc+1),nao,zero,tijab,nao)
!
! Fourth integral transformation
!   (ai|lj) tijab[l,a] -> (ai|bj) xaibj[a,b]
!   (ai|lj) tijab[l,a] -> (ai|kj) xaikj[a,k]
!
          call dgemm('T','N',nvac,nvir,nao,one,tijab,nao,cmo(1,nocc+1),nao,zero,xaibj,nvac)
          call dgemm('T','N',nvac,nocc,nao,one,tijab,nao,cmo,nao,zero,xaikj,nvac)
!
! Calculate MP2 amplitute and energy, Pab[I], Wab[I], and Wai[I](=Lai3) terms
!
          call mp2gradamp(tijab,pab,wab,wai,xaibj,xaikj,emp2st,energymo, &
&                         nocc,ncore,nvir,nvac,nvfz,moi,moj,nproc,myrank)
        endif
!
! Send and receive tijab and (ai|bj)/Dijab, and Calculate Pij
!
        call mp2gradpij(pij,tijab,xaibj,recvint,sendint,nocc,noac,ncore,nvac, &
&                       numirecv,numab,icycle,numij,idis,nproc,myrank,mpi_comm)
!
! First half back-transformation
!   tijab[a,b] -> tijml[l,m]
!
        if(myij <= numij) then
          call dgemm('N','T',nvac,nao,nvac,one,tijab,nvac,cmo(1,nocc+1),nao,zero,xaibj,nvac)
          call dgemm('T','T',nao,nao,nvac,one,xaibj,nvac,cmo(1,nocc+1),nao,zero,tijab,nao)
        endif
!
! Send and receive tijml, and calculate Wij[I]
!
        call mp2gradwij1(wij,tijab,trint2,trint2core,xaibj,recvt,storet,nocc, &
&                        noac,ncore,numitrans,maxsize,icycle,numij,idis,nproc,myrank,mpi_comm)
      enddo
!
      return
end


!-------------------------------------------------------------------------
  subroutine mp2gradamp(tijab,pab,wab,wai,xaibj,xaikj,emp2st,energymo, &
&                       nocc,ncore,nvir,nvac,nvfz,moi,moj,nproc,myrank)
!-------------------------------------------------------------------------
!
! Calculate MP2 amplitute and energy, Pab[I], Wab[I], and Wai[I](=Lai3) terms
!
! In:   xaikj    ((ai|kj))
! Inout:pab      (Pab)
!       wab      (Wab)
!       wai      (Wai)
!       xaibj    (In:(ai|bj), Out:(ai|bj)/Dijab)
!       emp2st   (Partial MP2 energies)
! Out:  tijab    (tijab)
!
      use modmolecule, only : nmo
      implicit none
      integer,intent(in) :: nocc, ncore, nvir, nvac, nvfz, moi, moj, nproc, myrank
      integer :: moa, mob, moc
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00, four=4.0D+00
      real(8),intent(in) :: xaikj(nvac,nocc), energymo(nmo)
      real(8),intent(inout) :: pab(nvir,nvir), wab(nvir*nvir), wai(nocc*nvac)
      real(8),intent(inout) :: xaibj(nvac,nvir), emp2st(2)
      real(8),intent(out) :: tijab(nvac,nvac)
      real(8) :: eij, eijb, eijab, etmp(2), eab
!
! Calculate MP2 amplitute and energy
!
      eij= energymo(moi+ncore)+energymo(moj+ncore)
      etmp(:)= zero
!$OMP parallel do private(eijb,eijab) reduction(+:etmp)
      do mob= 1,nvac
        eijb= eij-energymo(mob+nocc)
        do moa= 1,nvac
          eijab= one/(eijb-energymo(moa+nocc))
          tijab(moa,mob)=(two*xaibj(moa,mob)-xaibj(mob,moa))*eijab
          etmp(1)= etmp(1)+xaibj(moa,mob)*xaibj(moa,mob)*eijab
          etmp(2)= etmp(2)+xaibj(moa,mob)*(xaibj(moa,mob)-xaibj(mob,moa))*eijab
        enddo
      enddo
!$OMP end parallel do
!
      emp2st(1)= emp2st(1)+etmp(1)
      emp2st(2)= emp2st(2)+etmp(2)
!
! Calculate Wab[I]
!
      call dgemm('N','T',nvac,nvac,nvac,-two,tijab,nvac,xaibj,nvac,one,wab,nvir)
!
! Calculate Wai[I] (=Lai3)
!
      call dgemm('T','N',nocc,nvac,nvac,-four,xaikj,nvac,tijab,nvac,one,wai,nocc)
!
! Calculate (ai|bj)/Dijab
!
!$OMP parallel do private(eijb)
      do mob= 1,nvac
        eijb= eij-energymo(mob+nocc)
        do moa= 1,nvac
          xaibj(moa,mob)= xaibj(moa,mob)/(eijb-energymo(moa+nocc))
        enddo
      enddo
!$OMP end parallel do
!
! Calculate Pab
!
      call dgemm('T','N',nvac,nvac,nvac,two,tijab,nvac,xaibj,nvac,one,pab,nvir)
!
      if(nvfz /= 0 ) then
!$OMP parallel do private(eab)
        do mob= nvac+myrank+1,nvir,nproc
          do moa= 1,nvac
            eab= two/(energymo(moa+nocc)-energymo(mob+nocc))
            do moc= 1,nvac
              pab(moa,mob)= pab(moa,mob)+eab*tijab(moc,moa)*xaibj(moc,mob)
            enddo
          enddo
        enddo
!$OMP end parallel do
      endif
!
      return
end


!--------------------------------------------------------------------------------
  subroutine mp2gradpij(pij,tijab,xaibj,recvint,sendint,nocc,noac,ncore,nvac, &
&                       numirecv,numab,icycle,numij,idis,nproc,myrank,mpi_comm)
!--------------------------------------------------------------------------------
!
! Calculate Pij term
!
! In:   tijab    (tijab)
!       xaibj    ((ai|bj))
! Inout:pij      (Pij)
! Work: recvint  (Buffer for receiving tijab and (ai|bj))
!       sendint  (Buffer for sending tijab and (ai|bj))
!
      implicit none
      integer,intent(in) :: nocc, noac, ncore, nvac, numirecv, numab, icycle, numij
      integer,intent(in) :: nproc, myrank, mpi_comm, idis(0:nproc-1,8)
      integer :: iproc, jproc, moab, moijrecv, moirecv, mojrecv
      integer :: moijfirst, moijlast, moifirst, moilast, mojlast, kfirst, klast, jlast
      integer :: moi, moj, moij, mok
      real(8),parameter :: two=2.0D+00
      real(8),intent(in) :: tijab(nvac*nvac), xaibj(nvac*nvac)
      real(8),intent(out) :: recvint(numab,2,noac,numirecv), sendint(numab,2)
      real(8),intent(inout) :: pij(nocc*(nocc+1)/2)
!
      jproc= myrank
!
! Send and receive tijab and xaibj
!
      do iproc= myrank+1,nproc-1
        jproc= jproc-1
        if(jproc < 0) jproc= nproc-1
        do moab= 1,idis(iproc,5)
          sendint(moab,1)= tijab(idis(iproc,6)+moab)
          sendint(moab,2)= xaibj(idis(iproc,6)+moab)
        enddo 

        moijrecv=(icycle-1)*nproc+1+jproc
        moirecv= mod((moijrecv-1)/noac,numirecv)+1
        mojrecv= mod(moijrecv-1,noac)+1
        call para_sendrecvr(sendint,2*numab,iproc,myrank, &
&                           recvint(1,1,mojrecv,moirecv),2*numab,jproc,jproc,mpi_comm)
      enddo
!
      do iproc= 0, myrank-1
        jproc= jproc-1
        if(jproc < 0) jproc= nproc-1
        do moab= 1,idis(iproc,5)
          sendint(moab,1)= tijab(idis(iproc,6)+moab)
          sendint(moab,2)= xaibj(idis(iproc,6)+moab)
        enddo

        moijrecv=(icycle-1)*nproc+1+jproc
        moirecv= mod((moijrecv-1)/noac,numirecv)+1
        mojrecv= mod((moijrecv-1),noac)+1
        call para_sendrecvr(sendint,2*numab,iproc,myrank, &
&                           recvint(1,1,mojrecv,moirecv),2*numab,jproc,jproc,mpi_comm)
      enddo
!
      moijrecv=(icycle-1)*nproc+1+myrank
      moirecv= mod((moijrecv-1)/noac,numirecv)+1
      mojrecv= mod(moijrecv-1,noac)+1
      if(moijrecv <= numij) then
        do moab= 1,idis(myrank,5)
          recvint(moab,1,mojrecv,moirecv)= tijab(idis(myrank,6)+moab)
          recvint(moab,2,mojrecv,moirecv)= xaibj(idis(myrank,6)+moab)
        enddo
      endif
!
! Calculate Pij
!
      moijfirst=(icycle-1)*nproc+1
      moijlast = icycle*nproc
      if(moijlast > numij) moijlast= numij
      moifirst=(moijfirst-1)/noac+1
      moilast =(moijlast -1)/noac+1
      mojlast = mod((moijlast-1),noac)+1
      kfirst= mod((moifirst-1),numirecv)+1
      klast = mod((moilast-1),numirecv)+1
      jlast = mod((mojlast-1),noac)+1
!
      if(kfirst == klast) then
        if(jlast == noac) then
!$OMP parallel do private(moij)
          do moi= 1,noac
            do moj= 1,moi
              moij=(moi+ncore)*(moi+ncore-1)/2+moj+ncore
              do moab= 1,idis(myrank,5)
                pij(moij)= pij(moij)-two*recvint(moab,1,moj,kfirst)*recvint(moab,2,moi,kfirst)
              enddo
            enddo
          enddo
!$OMP end parallel do
        endif
      elseif(kfirst < klast) then
        if(jlast /= noac) klast= klast-1
!!$OMP parallel private(moij)
        do mok= kfirst,klast
!!$OMP do
          do moi= 1,noac
            do moj= 1,moi
              moij=(moi+ncore)*(moi+ncore-1)/2+moj+ncore
              do moab= 1,idis(myrank,5)
                pij(moij)= pij(moij)-two*recvint(moab,1,moj,mok)*recvint(moab,2,moi,mok)
              enddo
            enddo
          enddo
!!$OMP end do
        enddo
!!$OMP end parallel
      else
        if(jlast /= noac) klast= klast-1
!!$OMP parallel private(moij)
        do mok= kfirst,numirecv
!!$OMP do
          do moi= 1,noac
            do moj= 1,moi
              moij=(moi+ncore)*(moi+ncore-1)/2+moj+ncore
              do moab= 1,idis(myrank,5)
                pij(moij)= pij(moij)-two*recvint(moab,1,moj,mok)*recvint(moab,2,moi,mok)
              enddo
            enddo
          enddo
!!$OMP end do
        enddo
        do mok= 1,klast
!!$OMP do
          do moi= 1,noac
            do moj= 1,moi
              moij=(moi+ncore)*(moi+ncore-1)/2+moj+ncore
              do moab= 1,idis(myrank,5)
                pij(moij)= pij(moij)-two*recvint(moab,1,moj,mok)*recvint(moab,2,moi,mok)
              enddo
            enddo
          enddo
!!$OMP end do
        enddo
!!$OMP end parallel
      endif
!
      return
end


!-----------------------------------------------------------------------------------------------
  subroutine mp2gradwij1(wij,tijml,trint2,trint2core,sendt,recvt,storet,nocc, &
&                        noac,ncore,numitrans,maxsize,icycle,numij,idis,nproc,myrank,mpi_comm)
!-----------------------------------------------------------------------------------------------
!
! Send and receive tijml, and calculate Wij[I]
!
! In:   tijml     (tijml)
!       trint2core(Second transformed integrals of core MOs)
! Out:  sendt     (Sending buffer for tijml)
!       recvt     (Receiving buffer for tijml)
!       storet    (Storing buffer for tijml)
! Inout:wij       (Wij)
!       trint2    (In:Second transformed integrals, Out:tijml)
!
      use modbasis, only : nao, mbf, locbf, nshell
      implicit none
      integer,intent(in) :: nocc, noac, ncore, numitrans, maxsize, icycle, numij
      integer,intent(in) :: nproc, myrank, mpi_comm, idis(0:nproc-1,8)
      integer :: myij, iproc, jproc, ish, ksh, jcount, num, nbfi, nbfk, locbfi, locbfk
      integer :: ik, ii, kk, ijstart, nrecv, nsend
      integer :: moikfirst, moiklast, mokfirst, moklast, moifirst, moilast
      integer :: moi, mok, mlao, moirecvt, numml
      real(8),parameter :: one=1.0D+00, two=2.0D+00
      real(8),intent(in) :: tijml(nao,nao), trint2core(idis(myrank,3),ncore,numitrans)
      real(8),intent(out) :: sendt(maxsize,0:nproc-1), recvt(maxsize,0:nproc-1)
      real(8),intent(out) :: storet(maxsize,noac)
      real(8),intent(inout) :: wij(nocc,nocc)
      real(8),intent(inout) :: trint2(idis(myrank,3),noac,numitrans)
!
! Copy tijml to sending buffer
!

      myij=(icycle-1)*nproc+1+myrank
      if(myij <= numij) then
        do iproc= 0,nproc-1
          ik= 0
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
                sendt(ik,iproc)= tijml(locbfk+kk,locbfi+ii)
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
! Send and receive tijml
!
      jproc= myrank
      ijstart=(icycle-1)*nproc+1
      nrecv= idis(myrank,3)
!
      do iproc= myrank+1,nproc-1
        jproc= jproc-1
        if(jproc < 0) jproc= nproc-1
        nsend= idis(iproc,3)
        call para_sendrecvr(sendt(1,iproc),nsend,iproc,myrank, &
&                           recvt(1,jproc),nrecv,jproc,jproc,mpi_comm)
      enddo
!
      do iproc= 0,myrank-1
        jproc= jproc-1
        if(jproc < 0) jproc= nproc-1
        nsend= idis(iproc,3)
!ishimura
!       if((ijstart+iproc) > numij) nsend= 0
        call para_sendrecvr(sendt(1,iproc),nsend,iproc,myrank, &
&                           recvt(1,jproc),nrecv,jproc,jproc,mpi_comm)
      enddo
!
      if(myij <= numij) call dcopy(nrecv,sendt(1,myrank),1,recvt(1,myrank),1)
!
! Calculate Wij[I] and copy tijml to trint2
!
      moikfirst=(icycle-1)*nproc+1
      moiklast = icycle*nproc
      if(moiklast > numij) moiklast= numij
      mokfirst= mod(((moikfirst-1)/noac),numitrans)+1
      moklast = mod(((moiklast -1)/noac),numitrans)+1
      moifirst= mod((moikfirst-1),noac)+1
      moilast = mod((moiklast -1),noac)+1
!
      if((mokfirst /= moklast).or.(moilast == noac)) then
        numml= idis(myrank,3)
        if((numml /= 0).and.(moifirst /= 1)) then
          call dgemm('T','N',noac,moifirst-1,numml,-two,trint2(1,1,mokfirst),numml, &
&                    storet(1,1),maxsize,one,wij(ncore+1,ncore+1),nocc)
          call dgemm('T','N',ncore,moifirst-1,numml,-two,trint2core(1,1,mokfirst),numml, &
&                    storet(1,1),maxsize,one,wij(1,ncore+1),nocc)
        endif
        if(numml /= 0) then
          call dgemm('T','N',noac,noac-moifirst+1,numml,-two,trint2(1,1,mokfirst),numml, &
&                    recvt(1,0),maxsize,one,wij(ncore+1,ncore+moifirst),nocc)
          call dgemm('T','N',ncore,noac-moifirst+1,numml,-two,trint2core(1,1,mokfirst),numml, &
&                    recvt(1,0),maxsize,one,wij(1,ncore+moifirst),nocc)
        endif
!
        do moi= 1,moifirst-1
          do mlao= 1,idis(myrank,3)
            trint2(mlao,moi,mokfirst)= storet(mlao,moi)
          enddo
        enddo
        do moi= moifirst,noac
          do mlao= 1,idis(myrank,3)
            trint2(mlao,moi,mokfirst)= recvt(mlao,moi-moifirst)
          enddo
        enddo
      endif
!
!ishimura
      do mok= mokfirst+1,moklast
        if((mok == moklast).and.(moilast /= noac)) cycle
        if(numml /= 0) then
          moirecvt= noac*(mok-mokfirst-1)+(noac-moifirst)+1
          call dgemm('T','N',noac,noac,numml,-two,trint2(1,1,mok),numml, &
&                    recvt(1,moirecvt),maxsize,one,wij(ncore+1,ncore+1),nocc)
          call dgemm('T','N',ncore,noac,numml,-two,trint2core(1,1,mok),numml, &
&                    recvt(1,moirecvt),maxsize,one,wij(1,ncore+1),nocc)
        endif
        do moi= 1,noac
          moirecvt= noac*(mok-mokfirst-1)+(noac-moifirst)+moi
          do mlao= 1,numml
            trint2(mlao,moi,mok)= recvt(mlao,moirecvt)
          enddo
        enddo
      enddo
!
      if(moilast /= noac) then
        if(mokfirst == moklast) then
          do moi= moifirst,moilast
            do mlao= 1,idis(myrank,3)
              storet(mlao,moi)= recvt(mlao,moi-moifirst)
            enddo
          enddo
        else
          do moi= 1,moilast
            moirecvt= moi+noac*(moklast-mokfirst-1)+(noac-moifirst)
            do mlao= 1,idis(myrank,3)
              storet(mlao,moi)= recvt(mlao,moirecvt)
            enddo
          enddo
        endif
      endif
!
      return
end


!-------------------------------------------------------------------------------------------
  subroutine mp2gradbacktrans1(egrad,egradtmp,xlai,tisml,xlmi,cmo,xint,tijml,work,work2, &
&                              nocc,noac,nvir,ncore,numitrans,istart,maxdim,maxgraddim, &
&                              mlsize2,idis,nproc,myrank)
!-------------------------------------------------------------------------------------------
!
! Driver of second-half back-transformation (tnsml), and Lai4 and tnsml*(mn|ls)' terms
!
! In:   cmo       (MO coefficient matrix)
!       xint      (Exchange integral matrix)
!       tijml     (tijml)
!       istart    (First index of transformed occupied MOs)
!       maxdim    (Maximum size of basis functions in a shell)
!       maxgraddim(Maximum size of basis functions in a shell for derivative calculation)
! Inout:egrad     (MP2 energy gradients)
!       xlai      (Lai)
! Work: egradtmp  (Workspace for MP2 energy gradients)
!       tisml     (tisml)
!       xlmi      (Lmi)
!
      use modbasis, only : nao, nshell, mbf
      use modmolecule, only : natom
      implicit none
      integer,intent(in) :: nocc, noac, nvir, ncore, numitrans, istart, maxdim, maxgraddim
      integer,intent(in) :: mlsize2, nproc, myrank, idis(0:nproc-1,8)
      integer :: mlcount, ish, ksh, ml, mlindex(4,idis(myrank,4)), jcount, ii, numshell
      integer :: numml, mlstart, mlend
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00, four=4.0D+00
      real(8),intent(in) :: cmo(nao,nao), xint(nshell*(nshell+1)/2)
      real(8),intent(in) :: tijml(idis(myrank,3),noac,numitrans)
      real(8),intent(out) :: egradtmp(3*natom), tisml(nao*mlsize2*numitrans)
      real(8),intent(out) :: xlmi(numitrans,nao), work(nao*nao), work2(nao*mlsize2)
      real(8),intent(inout) :: egrad(3*natom), xlai(nocc,nvir) 
!
!
      mlcount= 0
      ish= idis(myrank,1)
      ksh= idis(myrank,2)
      jcount= 0
      do ml= 1,idis(myrank,4)
        mlindex(1,ml)= ish
        mlindex(2,ml)= ksh
        mlindex(3,ml)= mlcount
        mlindex(4,ml)= mbf(ish)*mbf(ksh)
        mlcount= mlcount+mlindex(4,ml)
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
      numml= 0
      mlstart= 1
!
      do numshell= 1,idis(myrank,4)
        numml= numml+mlindex(4,numshell)
        if(numshell == idis(myrank,4)) then
          mlend= numshell
          call mp2gradbt1(xlai,egrad,tijml,cmo,xint,tisml,xlmi,egradtmp,work,work2, &
&                         mlindex,numml,mlstart,mlend,numitrans,nocc,noac,ncore,nvir,maxdim, &
&                         maxgraddim,istart,idis,nproc,myrank)
          exit
        endif
!
        if(numml+mlindex(4,numshell+1) > mlsize2) then
          mlend= numshell
          call mp2gradbt1(xlai,egrad,tijml,cmo,xint,tisml,xlmi,egradtmp,work,work2, &
&                         mlindex,numml,mlstart,mlend,numitrans,nocc,noac,ncore,nvir,maxdim, &
&                         maxgraddim,istart,idis,nproc,myrank)
          mlstart= numshell+1
          numml= 0
        endif
      enddo
!
      return
end


!-------------------------------------------------------------------------------------------------
  subroutine mp2gradbacktrans2(egrad,egradtmp,wij,wab,wai,xlai,tisml,xlmi,xlmn,cmo,xint,tijml, &
&                              pij,pab,pmn,energymo,work,work2,nocc,noac,nvir,ncore,numitrans, &
&                              istart,maxdim,maxgraddim,mlsize2,idis,nproc,myrank,mpi_comm)
!-------------------------------------------------------------------------------------------------
!
! Driver of second-half back-transformation (tnsml), and PiJ, Lai1,2,4 and tnsml*(mn|ls)' terms,
! and accumulation of Pij and Pab
!
! In:   wai       (Wai)
!       cmo       (MO coefficient matrix)
!       xint      (Exchange integral matrix)
!       tijml     (tijml)
!       pab       (Pab)
!       energymo  (MO energies)
!       istart    (First index of transformed occupied MOs)
!       maxdim    (Maximum size of basis functions in a shell)
!       maxgraddim(Maximum size of basis functions in a shell for derivative calculation)
!       mlsize2   (Maximum size of basis function pairs for second-half back-transformation)
! Out:  pmn       (Pmn)
! Inout:egrad     (MP2 energy gradients)
!       wij       (Wij)
!       wab       (Wab)
!       xlai      (Lai)
!       pij       (Pij)
! Work: egradtmp  (Workspace for MP2 energy gradients)
!       tisml     (tisml)
!       xlmi      (Lmi)
!       xlmn      (Lmn)
!       work      (Work space)
!       work2     (Work space)
!
      use modbasis, only : nao, nshell, mbf
      use modmolecule, only : natom
      implicit none
      integer,intent(in) :: nocc, noac, nvir, ncore, numitrans, istart, maxdim, maxgraddim
      integer,intent(in) :: mlsize2, nproc, myrank, mpi_comm, idis(0:nproc-1,8)
      integer :: mlcount, ish, ksh, ml, mlindex(4,idis(myrank,4)), jcount, ii, numshell
      integer :: numml, moi, moj
      integer :: moij, ij, jj, mlstart, mlend
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, two=2.0D+00, four=4.0D+00
      real(8),intent(in) :: wai(nocc*nvir), cmo(nao,nao), xint(nshell*(nshell+1)/2)
      real(8),intent(in) :: tijml(idis(myrank,3),noac,numitrans)
      real(8),intent(in) :: pab(nvir*nvir), energymo(nao)
      real(8),intent(out) :: egradtmp(3*natom), tisml(nao*mlsize2*numitrans)
      real(8),intent(out) :: xlmi(nao,numitrans), xlmn(nao,nao), pmn(nao,nao), work(nao*nao)
      real(8),intent(out) :: work2(nao*mlsize2)
      real(8),intent(inout) :: egrad(3*natom), wij(nocc,nocc), wab(nvir,nvir)
      real(8),intent(inout) :: xlai(nocc,nvir), pij(nocc*(nocc+1)/2)
      real(8) :: eij
!
      mlcount= 0
      ish= idis(myrank,1)
      ksh= idis(myrank,2)
      jcount= 0
      do ml= 1,idis(myrank,4)
        mlindex(1,ml)= ish
        mlindex(2,ml)= ksh
        mlindex(3,ml)= mlcount
        mlindex(4,ml)= mbf(ish)*mbf(ksh)
        mlcount= mlcount+mlindex(4,ml)
        if(jcount == myrank) then
          ksh= ksh+2*nproc-1
        else
          ksh= ksh+nproc-1
        endif
        jcount= jcount+1
        if(jcount == nproc) jcount= 0
        if(ksh > nshell) then
          do ii= 1,nshell
            ksh= ksh-nshell
            ish= ish+1
            if(ksh <= nshell) exit
          enddo
        endif
      enddo

!
! Calculate PiJ
!
      do moi= ncore+1,nocc
        do moj= 1,ncore
          moij= moi*(moi-1)/2+moj
          eij=-one/(energymo(moi)-energymo(moj))
          pij(moij)= pij(moij)+eij*wij(moj,moi)
        enddo
      enddo
!
! Accumulate Pij, and calculate Wij[II] and Pmn
!
      call para_allreducer(pij,work,nocc*(nocc+1)/2,mpi_comm)
      do ii= 1,nocc
        ij= ii*(ii-1)/2
        do jj= 1,ii
          xlmn(jj,ii)= work(ij+jj)
        enddo
      enddo
      do ii= myrank+1,nocc,nproc
        do jj= 1,ii-1
          wij(ii,jj)= wij(ii,jj)-half*xlmn(jj,ii)*(energymo(ii)+energymo(jj))
          wij(jj,ii)= wij(jj,ii)-half*xlmn(jj,ii)*(energymo(ii)+energymo(jj))
        enddo
        wij(ii,ii)= wij(ii,ii)-half*xlmn(ii,ii)*(energymo(ii)+energymo(ii))
      enddo
      call dsymm('R','U',nao,nocc,one,xlmn,nao,cmo,nao,zero,work,nao)
      call dgemm('N','T',nao,nao,nocc,one,work,nao,cmo,nao,zero,pmn,nao)
!
! Accumulate Pab, and calculate Wab[II] and Pmn
!
      call para_allreducer(pab,work,nvir*nvir,mpi_comm)
      do ii= 1,nvir
        ij= nvir*(ii-1)
        do jj= 1,ii
          xlmn(jj,ii)= work(ij+jj)
        enddo
      enddo
      do ii= myrank+1,nvir,nproc
        ij= nvir*(ii-1)
        do jj= 1,ii-1
          wab(jj,ii)= wab(jj,ii)-half*xlmn(jj,ii)*(energymo(ii+nocc)+energymo(jj+nocc))
          wab(ii,jj)= wab(ii,jj)-half*xlmn(jj,ii)*(energymo(ii+nocc)+energymo(jj+nocc))
        enddo
        wab(ii,ii)= wab(ii,ii)-half*xlmn(ii,ii)*(energymo(ii+nocc)+energymo(ii+nocc))
      enddo
      call dsymm('R','U',nao,nvir,one,xlmn,nao,cmo(1,nocc+1),nao,zero,work,nao)
      call dgemm('N','T',nao,nao,nvir,one,work,nao,cmo(1,nocc+1),nao,one,pmn,nao)
!
!
      numml= 0
      mlstart= 1
!
      do numshell= 1,idis(myrank,4)
        numml= numml+mlindex(4,numshell)
        if(numshell == idis(myrank,4)) then
          mlend= numshell
          call mp2gradbt2(xlai,egrad,tijml,cmo,pmn,xint,tisml,xlmi,xlmn,egradtmp,work,work2, &
&                         mlindex,numml,mlstart,mlend,numitrans,nocc,noac,ncore,nvir,maxdim, &
&                         maxgraddim,istart,idis,nproc,myrank)
          exit
        endif
!
        if(numml+mlindex(4,numshell+1) > mlsize2) then
          mlend= numshell
          call mp2gradbt2(xlai,egrad,tijml,cmo,pmn,xint,tisml,xlmi,xlmn,egradtmp,work,work2, &
&                         mlindex,numml,mlstart,mlend,numitrans,nocc,noac,ncore,nvir,maxdim, &
&                         maxgraddim,istart,idis,nproc,myrank)
          mlstart= numshell+1
          numml= 0
        endif
      enddo
!
!
! Add Wai[I] as Lai4
!
      call daxpy(nocc*nvir,one,wai,1,xlai,1)
      if(nproc /= 1) then
        call para_allreducer(xlai,work,nocc*nvir,mpi_comm)
        call dcopy(nocc*nvir,work,1,xlai,1)
      endif
!
      return
end

!------------------------------------------------------------------------------------------------
  subroutine mp2gradcphf(pai,xlai,cmo,xint,energymo,paiprev,pls,pmax,paifock,errdiis,paidiis, &
&                        diismtrx,work1,work2,nocc,nvir,maxdim,idis,nproc,myrank,mpi_comm)
!------------------------------------------------------------------------------------------------
!
! Solve CPHF equation for MP2 energy gradient
!
! In:   xlai      (Lai)
!       cmo       (MO coefficient matrix)
!       xint      (Exchange integral matrix)
!       energymo  (MO energies)
! Out:  pai       (Pai)
! Work: pairev    (Previous Pai)
!       pls       (Pls (AO basis))
!       pmax      (Maximum Pls in a shell pair)
!       paifock   (New Pls)
!       errdiis   (Error matrix for DIIS)
!       paidiis   (Old Pai for DIIS)
!       diismtrx  (DIIS matrix)
!
      use modparallel, only : master
      use modbasis, only : nao, nshell
      use modmp2, only : maxmp2diis, maxmp2iter
      use modthresh, only : threshmp2cphf
      use modparallel, only : master
      implicit none
      integer,intent(in) :: nocc, nvir, maxdim, nproc, myrank, mpi_comm, idis(0:nproc-1,8)
      integer :: moa, moi, itdiis, iter, ii, ij, jj, ia
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: xlai(nocc,nvir), cmo(nao,nao), xint(nshell*(nshell+1)/2)
      real(8),intent(in) :: energymo(nao)
      real(8),intent(out) :: pai(nocc,nvir), paiprev(nocc,nvir), pls(nao*(nao+1)/2)
      real(8),intent(out) :: pmax(nshell*(nshell+1)/2), paifock(nao*(nao+1)/2)
      real(8),intent(out) :: errdiis(idis(myrank,7)*maxmp2diis)
      real(8),intent(out) :: paidiis(idis(myrank,7)*maxmp2diis)
      real(8),intent(out) :: diismtrx(maxmp2diis*(maxmp2diis+1)/2)
      real(8),intent(out) :: work1(nao,nao), work2(nao*nao)
      real(8) :: deltapai
!
! Set initial Pai
!
      do moa= 1,nvir
        do moi= 1,nocc
          paiprev(moi,moa)=-xlai(moi,moa)/(energymo(moa+nocc)-energymo(moi))
        enddo
      enddo
!
! Start CPHF iteration
!
      itdiis= 0
      do iter= 1,maxmp2iter
!
! Calculate Pls 
!
        call dgemm('N','T',nocc,nao,nvir,one,paiprev,nocc,cmo(1,nocc+1),nao,zero,work2,nocc)
        call dgemm('N','N',nao,nao,nocc,one,cmo,nao,work2,nocc,zero,work1,nao)
!
!$OMP parallel do private(ij)
        do ii= 1,nao
          ij= ii*(ii-1)/2
          do jj= 1,ii
            pls(ij+jj)= work1(ii,jj)+work1(jj,ii)
          enddo
        enddo
!$OMP end parallel do
!
! Calculate maximum Pls of each shell
!
        call calcrdmax(pls,pmax,work1,nproc,myrank,mpi_comm)
!
! Calculate two-electron integrals and Fock-like matrix
!
        call formrfock(paifock,work1,pls,pmax,xint,maxdim,nproc,myrank,mpi_comm)
!
! Transform Fock-like matrix to MO basis
!
        do ii= 1,nao
          ij= ii*(ii-1)/2
          do jj= 1,ii
            work1(jj,ii)= paifock(ij+jj)
          enddo
        enddo
!            
        call dsymm('L','U',nao,nocc,one,work1,nao,cmo,nao,zero,work2,nao)
        call dgemm('T','N',nocc,nvir,nao,one,work2,nao,cmo(1,nocc+1),nao,zero,pai,nocc)
!
! Calculate new Pai and error matrix
!
!$OMP parallel do private(ia)
        do moa= 1,nvir
          ia=(moa-1)*nocc
          do moi= 1,nocc
            pai(moi,moa)=-(xlai(moi,moa)+pai(moi,moa))/(energymo(moa+nocc)-energymo(moi))
            work2(ia+moi)= pai(moi,moa)-paiprev(moi,moa)
          enddo
        enddo
!$OMP end parallel do
!
        if(iter >= 2) then
          itdiis= itdiis+1
          call mp2graddiis(pai,work2,errdiis,paidiis,diismtrx,work1,nocc,nvir,maxmp2diis, &
&                          itdiis,idis,nproc,myrank,mpi_comm)
        endif
!
! Check convergence
!
        deltapai= zero
!$OMP parallel do reduction(+:deltapai)
        do moa= 1,nvir
          do moi= 1,nocc
            deltapai= deltapai+(pai(moi,moa)-paiprev(moi,moa))**2
          enddo
        enddo
!$OMP end parallel do
!
        deltapai= sqrt(deltapai/(nocc*nvir))
!ishimura
        if(master) write(*,'(6x,"Cycle",i3,3x,"Z-Vector error=",1p,d11.3)')iter,deltapai
        if(deltapai < threshmp2cphf) exit
!
        if(itdiis == maxmp2diis) itdiis= 0
        if(iter == maxmp2iter) then
          if(master) write(*,'(" Error! MP2 CPHF calculation did not converge.")')
          call iabort
        endif
!
        do moa= 1,nvir
          do moi= 1,nocc
            paiprev(moi,moa)= pai(moi,moa)
          enddo
        enddo
      enddo
!
      return
end

!---------------------------------------------------------------------------------------
  subroutine mp2graddiis(pai,err,errdiis,paidiis,diismtrx,work,nocc,nvir,maxmp2diis, &
&                        itdiis,idis,nproc,myrank,mpi_comm)
!---------------------------------------------------------------------------------------
!
! DIIS interpolation to solve MP2 CPHF equation
!
! In   :err       (Error matrix for DIIS)
! Inout:pai       (Pai)
!       errdiis   (Storage of error matrix for DIIS)
!       paidiis   (Storage of Pai for DIIS)
!       diismtrx  (DIIS matrix)
! Work :work
!
      use modbasis, only : nao
      implicit none
      integer,intent(in) :: nocc, nvir, maxmp2diis, itdiis, nproc, myrank, mpi_comm
      integer,intent(in) :: idis(0:nproc-1,8)
      integer :: num, istart, ii, jj, ij, ipiv(maxmp2diis+1), info
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(inout) :: err(nocc*nvir)
      real(8),intent(out) :: work(nao*nao)
      real(8),intent(inout) :: pai(nocc*nvir), errdiis(idis(myrank,7),maxmp2diis)
      real(8),intent(inout) :: paidiis(idis(myrank,7),maxmp2diis)
      real(8),intent(inout) :: diismtrx(maxmp2diis*(maxmp2diis+1)/2)
      real(8) :: ddot, diisb(maxmp2diis+1,maxmp2diis+1), diiscoeff(maxmp2diis+1)
!
      num= idis(myrank,7)
      istart= idis(myrank,8)+1
      work(1:itdiis)= zero
      if(num > 0) then
        call dcopy(num,err(istart),1,errdiis(1,itdiis),1)
        do ii= 1,itdiis
          work(ii)= ddot(num,errdiis(1,ii),1,errdiis(1,itdiis),1)
        enddo
      endif
      call para_allreducer(work,diismtrx(itdiis*(itdiis-1)/2+1),itdiis,mpi_comm)
!
      do ii= 1,itdiis
        do jj= 1,ii
          ij= ii*(ii-1)/2+jj
          diisb(jj,ii)= diismtrx(ij)
          diisb(ii,jj)= diismtrx(ij)
        enddo
      enddo
      do jj= 1,itdiis
        diisb(itdiis+1,jj)=-one
        diisb(jj,itdiis+1)=-one
        diiscoeff(jj)= zero
      enddo
      diisb(itdiis+1,itdiis+1)= zero
      diiscoeff(itdiis+1)=-one
!
      call dgesv(itdiis+1,1,diisb,maxmp2diis+1,ipiv,diiscoeff,maxmp2diis+1,info)
!
      if(num > 0) then
        call dcopy(num,pai(istart),1,paidiis(1,itdiis),1)
        work(1:num)= zero
        do ii= 1,itdiis
          call daxpy(num,diiscoeff(ii),paidiis(1,ii),1,work,1)
        enddo
      endif
      call para_allgathervr(work,num,pai,idis(0,7),idis(0,8),nproc,mpi_comm)
!
      return
end


!--------------------------------------------------------------------------------------------
  subroutine mp2gradwij3(wij,wai,pmn,pai,pmnfock,cmo,energymo,xint,pmn2,pmax,work1,work2, &
&                        nocc,nvir,nvac,maxdim,nproc,myrank,mpi_comm)
!--------------------------------------------------------------------------------------------
!
! Calculate Wij[III] and Wai[II]
!
! In   :pai       (Pai)
!       cmo       (MO coefficient matrix)
!       energymo  (MO energies)
!       xint      (Exchange integral matrix)
! Inout:wij       (Wij)
!       wai       (Wai)
!       pmn       (Pmn (AO bais)
! Work :pmnfock,pmn2,pmax,work1,work2
!
      use modbasis, only : nao, nshell
      implicit none
      integer,intent(in) :: nocc, nvir, nvac, maxdim, nproc, myrank, mpi_comm
      integer :: moa, moi, moj, ii, ij, jj
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00
      real(8),intent(in) :: pai(nocc,nvir), cmo(nao,nao), energymo(nao), xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: pmnfock(nao*(nao+1)/2), pmn2(nao*(nao+1)/2)
      real(8),intent(out) :: pmax(nshell*(nshell+1)/2), work1(nao,nao), work2(nao*nao)
      real(8),intent(inout) :: wij(nocc,nocc), wai(nocc,nvac), pmn(nao,nao)
!
! Calculate Wai[II]
!
      do moa= myrank+1,nvac,nproc
        do moi= 1,nocc
          wai(moi,moa)= wai(moi,moa)-pai(moi,moa)*energymo(moi)
        enddo
      enddo
!
! Add Pai to Pmn
!
      call dgemm('N','T',nocc,nao,nvir,one,pai,nocc,cmo(1,nocc+1),nao,zero,work1,nocc)
      call dgemm('N','N',nao,nao,nocc,half,cmo,nao,work1,nocc,one,pmn,nao)
      call dgemm('T','T',nao,nao,nocc,half,work1,nocc,cmo,nao,one,pmn,nao)
!
! Calculate maximum Pmn of each shell
!
!$OMP parallel do private(ij)
      do ii= 1,nao
        ij= ii*(ii-1)/2
        do jj= 1,ii
          pmn2(ij+jj)= pmn(jj,ii)
        enddo
      enddo
!$OMP end parallel do
!
      call calcrdmax(pmn2,pmax,work1,nproc,myrank,mpi_comm)
!
! Calculate two-electron integrals and Fock-like matrix
!
      call formrfock(pmnfock,work1,pmn2,pmax,xint,maxdim,nproc,myrank,mpi_comm)
!
! Calculate Wij[III]
!
!$OMP parallel do private(ij)
      do ii= 1,nao
        ij= ii*(ii-1)/2
        do jj= 1,ii
          work1(jj,ii)= pmnfock(ij+jj)
        enddo
      enddo
!$OMP end parallel do
!
      call dsymm('L','U',nao,nocc,one,work1,nao,cmo,nao,zero,work2,nao)
      call dgemm('T','N',nocc,nocc,nao,-one,work2,nao,cmo,nao,zero,work1,nao)
!
      do moi= myrank+1,nocc,nproc
        do moj= 1,moi-1
          wij(moj,moi)= wij(moj,moi)+work1(moj,moi)
          wij(moi,moj)= wij(moi,moj)+work1(moj,moi)
        enddo
        wij(moi,moi)= wij(moi,moi)+work1(moi,moi)
      enddo
!
      return
end


!--------------------------------------------------------------------------------------------
  subroutine mp2graddint(egrad,egradtmp,pmn,wij,wab,wai,xint,cmo,energymo,wmn,pmnhf,work, &
&                        nocc,nvir,maxdim,maxgraddim,nproc,myrank,mpi_comm)
!--------------------------------------------------------------------------------------------
!
! Calculate integral derivatives and MP2 energy gradient
!
! In   :wij       (Wij)
!       wab       (Wab)
!       wai       (Wai)
!       xint      (Exchange integral matrix)
!       cmo       (MO coefficient matrix)
!       energymo  (MO energies)
! Inout:egrad     (MP2 energy gradients)
!       pmn       (Pmn)
! Work :egradtmp,wmn,work
!
      use modbasis, only : nao, nshell
      use modmolecule, only : natom
      implicit none
      integer,intent(in) :: nocc, nvir, maxdim, maxgraddim, nproc, myrank, mpi_comm
      integer :: ii, jj, ij, kk
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, two=2.0D+00
      real(8),intent(in) :: wij(nocc*nocc), wab(nvir*nvir), wai(nocc*nvir)
      real(8),intent(in) :: xint(nshell*(nshell+1)/2), cmo(nao,nao), energymo(nao)
      real(8),intent(out) :: egradtmp(3*natom), wmn(nao*(nao+1)/2), pmnhf(nao,nao), work(nao,nao)
      real(8),intent(inout) :: egrad(3*natom), pmn(nao,nao)
      real(8) :: egradtmp2(3*natom), ewtmp
!
      egradtmp(:)= zero
!
! Calculate energy gradient of nuclear repulsion
!
      call nucgradient(egradtmp,nproc,myrank)
!
! Calculate HF density matrix
!
      call dgemm('N','T',nao,nao,nocc,two,cmo,nao,cmo,nao,zero,pmnhf,nao)
!
! Calculate derivatives for two-electron integrals
!
      call grad2eri(egradtmp,egradtmp2,pmnhf,pmn,xint,one, &
&                   maxdim,maxgraddim,nproc,myrank,3)
!
! Calculate HF+MP2 density matrix
!
      call daxpy(nao*nao,one,pmnhf,1,pmn,1)
!
! Calculate Wmn
!
      call dgemm('N','T',nocc,nao,nocc,one,wij,nocc,cmo,nao,zero,pmnhf,nocc)
      call dgemm('N','T',nocc,nao,nvir,one,wai,nocc,cmo(1,nocc+1),nao,one,pmnhf,nocc)
      call dgemm('N','N',nao,nao,nocc,one,cmo,nao,pmnhf,nocc,zero,work,nao)
      call dgemm('N','T',nvir,nao,nvir,one,wab,nvir,cmo(1,nocc+1),nao,zero,pmnhf,nvir)
      call dgemm('N','N',nao,nao,nvir,one,cmo(1,nocc+1),nao,pmnhf,nvir,one,work,nao)
      call para_allreducer(work,pmnhf,nao*nao,mpi_comm)
!
      do ii= 1,nocc
        do jj= 1,nao
          work(ii,jj)= cmo(jj,ii)
        enddo
      enddo
      do ii= nao,1,-1
        ij= ii*(ii-1)/2
        do jj= 1,ii
          ewtmp= zero
          do kk= 1,nocc
            ewtmp= ewtmp-work(kk,ii)*work(kk,jj)*energymo(kk)
          enddo
          wmn(ij+jj)= ewtmp*two+(pmnhf(ii,jj)+pmnhf(jj,ii))*half
        enddo
      enddo
!
! Calculate derivatives of one-electron integrals
!
      call gradoneei(egradtmp,egradtmp2,pmn,wmn,nproc,myrank)
!
      do ii= 1,3*natom
        egrad(ii)= egrad(ii)+egradtmp(ii)
      enddo
!
      return
end


!---------------------------------------------------------------------------------------------
  subroutine mp2gradbt1(xlai,egrad,tijml,cmo,xint,tisml,xlmi,egradtmp,cmotrans,work2, &
&                       mlindex,numml,mlstart,mlend,numitrans,nocc,noac,ncore,nvir,maxdim, &
&                       maxgraddim,istart,idis,nproc,myrank)
!---------------------------------------------------------------------------------------------
!
! Calculate second-half back-transformation (tnsml), Lai4 and tnsml*(mn|ls)' terms
!
! In:   tijml     (tijml)
!       cmo       (MO coefficient matrix)
!       xint      (Exchange integral matrix)
!       istart    (First index of transformed occupied MOs)
!       maxdim    (Maximum size of basis functions in a shell)
!       maxgraddim(Maximum size of basis functions in a shell for derivative calculation)
! Inout:xlai      (Lai)
!       egrad     (MP2 energy gradients)
! Work: egradtmp  (Workspace for MP2 energy gradients)
!       tisml     (tisml)
!       xlmi      (Lmi)
!       cmotrans  (transposed matrix of cmo)
!       work2     (work space)
!
      use modbasis, only : nao, nshell, mbf, locbf
      use modmolecule, only : natom
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: numml, mlstart, mlend, numitrans, nocc, noac, ncore, nvir, maxdim
      integer,intent(in) :: maxgraddim, istart, nproc, myrank, idis(0:nproc-1,8)
      integer,intent(in) :: mlindex(4,idis(myrank,4))
      integer :: moi, mlcount, numshell, ish, jsh, ksh, lsh, nbfi, nbfj, nbfk, nbfl
      integer :: locbfi, locbfj, locbfk, locbfl, ij, kl, ii, jj, kk, ll, i2, ik, ml
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00, four=4.0D+00
      real(8),intent(in) :: tijml(idis(myrank,3),noac,numitrans), cmo(nao,nao)
      real(8),intent(in) :: xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: tisml(numitrans,nao,numml), xlmi(numitrans,nao)
      real(8),intent(out) :: egradtmp(3*natom), cmotrans(noac,*), work2(nao,numml)
      real(8),intent(inout) :: xlai(nocc,nvir), egrad(3*natom)
      real(8) :: twoeri(maxgraddim**4), twork(maxdim,maxdim,maxdim,maxdim)
      real(8) :: dtwoeri(maxdim,maxdim,maxdim,maxdim,3), tmax, tmp, xijkl
!
      do moi= 1,numitrans
        call dgemm('N','T',nao,numml,noac,one,cmo(1,ncore+1),nao, &
&                  tijml(mlindex(3,mlstart)+1,1,moi),idis(myrank,3),zero,work2,nao)
!$OMP parallel do
        do ml= 1,numml
          do ll= 1,nao
            tisml(moi,ll,ml)=work2(ll,ml)
          enddo
        enddo
!$OMP end parallel do
      enddo
!
      do ii= 1,noac
        do jj= 1,nao
          cmotrans(ii,jj)= cmo(jj,ii+ncore)
        enddo
      enddo
!
      xlmi(:,:)= zero
      egradtmp(:)= zero
!$OMP parallel do schedule(dynamic,1) collapse(2) private(ish,ksh,mlcount,nbfi,nbfj,nbfk,nbfl, &
!$OMP locbfi,locbfj,locbfk,locbfl,ij,kl,i2,ik,xijkl,tmax,tmp,twoeri,twork,dtwoeri) &
!$OMP reduction(+:xlmi,egradtmp)
      do numshell= mlstart,mlend
        do lsh= 1,nshell
          ish= mlindex(1,numshell)
          ksh= mlindex(2,numshell)
          mlcount= mlindex(3,numshell)-mlindex(3,mlstart)
          nbfi= mbf(ish)
          nbfk= mbf(ksh)
          nbfl= mbf(lsh)
          locbfi= locbf(ish)
          locbfk= locbf(ksh)
          locbfl= locbf(lsh)
          if(ksh >= lsh) then
            kl= ksh*(ksh-1)/2+lsh
          else
            kl= lsh*(lsh-1)/2+ksh
          endif
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
            xijkl= xint(ij)*xint(kl)
            if(xijkl < cutint2) cycle
            call calc2eri(dtwoeri,ish,jsh,ksh,lsh,maxdim)
!
            tmax= zero
            do ii= 1,nbfi
              i2=(ii-1)*nbfk+mlcount
              do kk= 1,nbfk
                ik= i2+kk
                do jj= 1,nbfj
                  do ll= 1,nbfl
                    tmp= zero
                    if(abs(dtwoeri(ll,kk,jj,ii,1)) > cutint2) then
                      do moi= 1,numitrans
                        xlmi(moi,locbfj+jj)= xlmi(moi,locbfj+jj) &
&                                           +tisml(moi,locbfl+ll,ik)*dtwoeri(ll,kk,jj,ii,1)
                        tmp= tmp+tisml(moi,locbfl+ll,ik)*cmotrans(istart+moi,locbfj+jj)
                      enddo
                    endif
                    twork(ll,kk,jj,ii)= tmp*two
                    if(abs(tmp) > tmax) tmax= abs(tmp)
                  enddo
                enddo
              enddo
            enddo
!
            if(xijkl*tmax < cutint2) cycle
            call calcd2eri(egradtmp,twork,twoeri,dtwoeri,ish,jsh,ksh,lsh,maxdim,maxgraddim)
          enddo
        enddo
      enddo
!$OMP end parallel do
!
      call daxpy(3*natom,one,egradtmp,1,egrad,1)
      call dgemm('N','N',numitrans,nvir,nao,four,xlmi,numitrans,cmo(1,nocc+1),nao,one, &
&                xlai(istart+ncore+1,1),nocc)
!
      return
end


!-------------------------------------------------------------------------------------------------
  subroutine mp2gradbt2(xlai,egrad,tijml,cmo,pmn,xint,tisml,xlmi,xlmn,egradtmp,cmotrans,work2, &
&                       mlindex,numml,mlstart,mlend,numitrans,nocc,noac,ncore,nvir,maxdim, &
&                       maxgraddim,istart,idis,nproc,myrank)
!-------------------------------------------------------------------------------------------------
!
! Calculate second-half back-transformation (tnsml), Lai1,2,4 and tnsml*(mn|ls)' terms
!
! In:   cmo       (MO coefficient matrix)
!       xint      (Exchange integral matrix)
!       tijml     (tijml)
!       pmn       (Pmn)
!       istart    (First index of transformed occupied MOs)
!       maxdim    (Maximum size of basis functions in a shell)
!       maxgraddim(Maximum size of basis functions in a shell for derivative calculation)
! Inout:egrad     (MP2 energy gradients)
!       xlai      (Lai)
! Work: egradtmp  (Workspace for MP2 energy gradients)
!       tisml     (tisml)
!       xlmi      (Lmi)
!       xlmn      (Lmn)
!       cmotrans  (transposed matrix of cmo)
!       work2     (Work space)
!
      use modbasis, only : nao, nshell, mbf, locbf
      use modmolecule, only : natom
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: numml, mlstart, mlend, numitrans, nocc, noac, ncore, nvir, maxdim
      integer,intent(in) :: maxgraddim, istart, nproc, myrank, idis(0:nproc-1,8)
      integer,intent(in) :: mlindex(4,idis(myrank,4))
      integer :: moi, mlcount, numshell, ish, jsh, ksh, lsh, nbfi, nbfj, nbfk, nbfl
      integer :: locbfi, locbfj, locbfk, locbfl, ij, kl, ii, jj, kk, ll, i2, ik, ml
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00, four=4.0D+00
      real(8),intent(in) :: tijml(idis(myrank,3),noac,numitrans), cmo(nao,nao), pmn(nao,nao)
      real(8),intent(in) :: xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: tisml(numitrans,nao,numml), xlmi(numitrans,nao), xlmn(nao,nao)
      real(8),intent(out) :: egradtmp(3*natom), cmotrans(noac,*), work2(nao,numml)
      real(8),intent(inout) :: xlai(nocc,nvir), egrad(3*natom)
      real(8) :: twoeri(maxgraddim**4), twork(maxdim,maxdim,maxdim,maxdim)
      real(8) :: dtwoeri(maxdim,maxdim,maxdim,maxdim,3), tmax, tmp, xijkl
!
      do moi= 1,numitrans
        call dgemm('N','T',nao,numml,noac,one,cmo(1,ncore+1),nao, &
&                  tijml(mlindex(3,mlstart)+1,1,moi),idis(myrank,3),zero,work2,nao)
!$OMP parallel do
        do ml= 1,numml
          do ll= 1,nao
            tisml(moi,ll,ml)=work2(ll,ml)
          enddo
        enddo
!$OMP end parallel do
      enddo
!
      do ii= 1,noac
        do jj= 1,nao
          cmotrans(ii,jj)= cmo(jj,ii+ncore)
        enddo
      enddo
!
      xlmi(:,:)= zero
      xlmn(:,:)= zero
      egradtmp(:)= zero
!$OMP parallel do schedule(dynamic,1) collapse(2) private(ish,ksh,mlcount,nbfi,nbfj,nbfk,nbfl, &
!$OMP locbfi,locbfj,locbfk,locbfl,ij,kl,i2,ik,xijkl,tmax,tmp,twoeri,twork,dtwoeri) &
!$OMP reduction(+:xlmi,xlmn,egradtmp)
      do numshell= mlstart,mlend
        do lsh= 1,nshell
          ish= mlindex(1,numshell)
          ksh= mlindex(2,numshell)
          mlcount= mlindex(3,numshell)-mlindex(3,mlstart)
          nbfi= mbf(ish)
          nbfk= mbf(ksh)
          nbfl= mbf(lsh)
          locbfi= locbf(ish)
          locbfk= locbf(ksh)
          locbfl= locbf(lsh)
          if(ksh >= lsh) then
            kl= ksh*(ksh-1)/2+lsh
          else
            kl= lsh*(lsh-1)/2+ksh
          endif
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
            xijkl= xint(ij)*xint(kl)
            if(xijkl < cutint2) cycle
            call calc2eri(dtwoeri,ish,jsh,ksh,lsh,maxdim)
!
            tmax= zero
            do ii= 1,nbfi
              i2=(ii-1)*nbfk+mlcount
              do kk= 1,nbfk
                ik= i2+kk
                do jj= 1,nbfj
                  do ll= 1,nbfl
                    tmp= zero
                    if(abs(dtwoeri(ll,kk,jj,ii,1)) > cutint2) then
                      do moi= 1,numitrans
                        xlmi(moi,locbfj+jj)= xlmi(moi,locbfj+jj) &
&                                           +tisml(moi,locbfl+ll,ik)*dtwoeri(ll,kk,jj,ii,1)
                        tmp= tmp+tisml(moi,locbfl+ll,ik)*cmotrans(istart+moi,locbfj+jj)
                      enddo
                    endif
                    twork(ll,kk,jj,ii)= tmp*two
                    if(abs(tmp) > tmax) tmax= abs(tmp)
                  enddo
                enddo
              enddo
            enddo
!
            do ii= 1,nbfi
              do kk= 1,nbfk
                do jj= 1,nbfj
                  do ll= 1,nbfl
                    xlmn(locbfl+ll,locbfk+kk)= xlmn(locbfl+ll,locbfk+kk)+pmn(locbfj+jj,locbfi+ii) &
&                                              *(four*dtwoeri(ll,kk,jj,ii,1))
                    xlmn(locbfl+ll,locbfi+ii)= xlmn(locbfl+ll,locbfi+ii)+pmn(locbfj+jj,locbfk+kk) &
&                                              *(-two*dtwoeri(ll,kk,jj,ii,1))
                  enddo
                enddo
              enddo
            enddo
!
            if(xijkl*tmax < cutint2) cycle
            call calcd2eri(egradtmp,twork,twoeri,dtwoeri,ish,jsh,ksh,lsh,maxdim,maxgraddim)
          enddo
        enddo
      enddo
!$OMP end parallel do
!
      call daxpy(3*natom,one,egradtmp,1,egrad,1)
      call dgemm('N','N',numitrans,nvir,nao,four,xlmi,numitrans,cmo(1,nocc+1),nao,one, &
&                xlai(istart+ncore+1,1),nocc)
      call dgemm('T','N',nocc,nao,nao,one,cmo,nao,xlmn,nao,zero,cmotrans,nocc)
      call dgemm('N','N',nocc,nvir,nao,one,cmotrans,nocc,cmo(1,nocc+1),nao,one,xlai,nocc)
!
      return
end
