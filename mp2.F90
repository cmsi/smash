!------------------------------------------
  subroutine calcrmp2(cmo,energymo,xint)
!------------------------------------------
!
! Driver of restricted MP2 calculation
!
! In  : cmo     (MO coefficient matrix)
!       energymo(MO energies)
!       xint    (Exchange integral matrix)
!       
      use modparallel
      use modbasis, only : nshell, nao, mbf
      use modmolecule, only : neleca, nmo
      use modenergy, only : enuc, escf, escfe, emp2, escsmp2
      implicit none
      integer :: ncore, ncorecalc, maxdim, nocc, nvir, nao2, ndis, iproc, icount, icountao
      integer :: icountsh, idis(8,0:nproc-1), maxsize, ish, jsh, ishs, jshs, msize, mlsize
      real(8),parameter :: zero=0.0D+00, three=3.0D+00, p12=1.2D+00
      real(8),intent(in) :: cmo(nao,nao), energymo(nmo), xint(nshell*(nshell+1)/2)
      real(8),allocatable :: trint2(:), cmowrk(:), trint1a(:), trint1b(:)
      real(8),allocatable :: trint3(:), trint4(:)
      real(8) :: emp2st(2), emp2stsum(2)
!ishimura
      integer :: jproc
!
      ncore= ncorecalc()
      emp2= zero
      emp2st(:)= zero
      maxdim=maxval(mbf(1:nshell))
      nocc= neleca-ncore
      nvir= nmo-neleca
      nao2= nao*nao
!
      if(master) then
        write(*,'(" ----------------------------------------------")')
        write(*,'("   MP2 calculation ")')
        write(*,'("     Number of basis functions         =",i5)')nao
        write(*,'("     Number of basis shells            =",i5)')nshell
        write(*,'("     Number of correlated occupied MOs =",i5)')nocc
        write(*,'("     Number of virtual MOs             =",i5)')nvir
        write(*,'(" ----------------------------------------------",/)')
      endif
!
      ndis= nao2/nproc
      icount= 0
      icountao= 0
      icountsh= 0
      idis= 0
      ishs= 1
      jshs= 1
!
      do iproc= 0,nproc-1
        idis(1,iproc)= ishs
        idis(2,iproc)= jshs
        jshs= jshs+1
        if(jshs > ishs) then
          ishs= ishs+1
          jshs= 1
        endif
      enddo
!
      iproc= 0
      do ish= 1,nshell
        do jsh= 1,ish
          idis(3,iproc)= idis(3,iproc)+mbf(ish)*mbf(jsh)
          idis(4,iproc)= idis(4,iproc)+1
          iproc= iproc+1
          if(iproc == nproc) iproc= 0
        enddo
      enddo
!ishimura
      ishs=1
      jshs=2
      do jproc=iproc,nproc-1
        idis(5,jproc)= ishs
        idis(6,jproc)= jshs
        jshs= jshs+1
        if(jshs > nshell) then
          ishs= ishs+1
          jshs= ishs+1
        endif
      enddo
      do jproc=0,iproc-1
        idis(5,jproc)= ishs
        idis(6,jproc)= jshs
        jshs= jshs+1
        if(jshs > nshell) then
          ishs= ishs+1
          jshs= ishs+1
        endif
      enddo

      do ish= 1,nshell
        do jsh= ish+1,nshell
          idis(7,iproc)= idis(7,iproc)+mbf(ish)*mbf(jsh)
          idis(8,iproc)= idis(8,iproc)+1
          iproc= iproc+1
          if(iproc == nproc) iproc= 0
        enddo
      enddo
!ishimura
!  if(master)write(*,'(8i6)')idis
!
      ishs=0
      do iproc=0,nproc-1
        ishs=max(idis(3,iproc)+idis(7,iproc),ishs)
      enddo
      maxsize= ishs
      call memset(maxsize*nocc*(nocc+1)/2)
      allocate(trint2(maxsize*nocc*(nocc+1)/2))
      call memset(nao*nocc+nocc*maxdim**3)
      allocate(cmowrk(nao*nocc),trint1a(nocc*maxdim**3))
!
      call memrest(msize)
      mlsize= msize/(nocc*nao)
      if(mlsize > maxsize) then
         mlsize= maxsize
      elseif(mlsize < maxdim*maxdim) then
         mlsize= maxdim*maxdim
      endif
      call memset(mlsize*nocc*nao)
      allocate(trint1b(mlsize*nocc*nao))
!
      call mp2trans1(cmo(1,ncore+1),cmowrk,trint1a,trint1b,trint2,xint, &
&                    mlsize,nocc,maxdim,idis)
!
      deallocate(trint1b)
      call memunset(mlsize*nocc*nao)
      deallocate(cmowrk,trint1a)
      call memunset(nao*nocc+nocc*maxdim**3)
!
      call memset(2*nao2)
      allocate(trint3(nao2),trint4(nao2))
!
      call mp2trans2(cmo(1,neleca+1),energymo,trint2,trint3,trint4,emp2st,nocc,nvir,ncore,idis)
!
      deallocate(trint3,trint4)
      call memunset(2*nao2)
      deallocate(trint2)
      call memunset(maxsize*nocc*(nocc+1)/2)
!
      call para_allreducer(emp2st,emp2stsum,2,MPI_COMM_WORLD)
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
      integer :: ncorecalc, ncore, iatom, numcore(137)
      data numcore/2*0, 8*1, 8*5, 18*9, 18*18, 32*27, 32*43, 19*59/
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


!----------------------------------------------------------------------
  subroutine mp2trans1(cmoocc,cmowrk,trint1a,trint1b,trint2,xint, &
&                      mlsize,nocc,maxdim,idis)
!----------------------------------------------------------------------
!
! Driver of AO intengral generation and first and second integral transformations
!
! In  : cmoocc  (MO coefficient matrix)
!       xint    (Exchange integral matrix)
!       mlsize  (Block size of first-transformed integrals)
!       nocc    (Number of active occupied MOs)
!       maxdim  (Maximum size of basis functions in a shell)
!       idis    (Information for parallelization)
! Out : trint2  (Second-transformed integrals)
! Work: cmowrk  (Transposed MO coefficient matrix)
!       trint1a (First-transformed integrals)
!       trint1b (First-transformed integrals)
!       
      use modparallel
      use modbasis, only : nshell, nao, mbf, locbf
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: maxdim, mlsize, nocc, idis(8,0:nproc-1)
      integer :: ish, ksh, mlcount, mlstart, numshell, ii
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: cmoocc(nao,nocc), xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: cmowrk(nocc,nao), trint1a(nocc,maxdim,maxdim**2)
      real(8),intent(out) :: trint1b(mlsize*nocc*nao), trint2((idis(3,myrank)+idis(7,myrank))*nocc*(nocc+1)/2)
!ishimura
    real(8) :: t1,t2,t3,tr1=0.0D0,tr2=0.0D0
!
      cmowrk=transpose(cmoocc)
!
      ish= idis(1,myrank)
      ksh= idis(2,myrank)
      mlcount= 0
      mlstart= 1
      do numshell= 1,idis(4,myrank)
!ishimura
    call cpu_time(t1)
!
! AO intengral generation and first integral transformation
!
        call transmoint1(trint1a,trint1b,cmowrk,xint,ish,ksh,maxdim,nocc,mlcount,mlsize)
!ishimura
    call cpu_time(t2)
        mlcount= mlcount+mbf(ish)*mbf(ksh)
        if(numshell == idis(4,myrank)) then
          call transmoint2(trint2,trint1b,cmoocc,nocc,mlcount,mlstart,mlsize,idis)
            mlstart= mlstart+mlcount
            mlcount= 0
        else
          ksh= ksh+nproc
          if(ksh > ish) then
            do ii= 1,nshell
              ksh= ksh-ish
              ish= ish+1
              if(ksh <= ish) exit
            enddo
          endif
!         ksh= ksh+1
!         if(ksh == nshell+1) then
!           ish= ish+1
!           ksh= 1
!         endif
          if(mlcount+mbf(ish)*mbf(ksh) > mlsize) then
!
! Second integral transformation
!
            call transmoint2(trint2,trint1b,cmoocc,nocc,mlcount,mlstart,mlsize,idis)
            mlstart= mlstart+mlcount
            mlcount= 0
          endif
        endif
!ishimura
    call cpu_time(t3)
    tr1=tr1+t2-t1
    tr2=tr2+t3-t2
      enddo
!ishimura
      ish= idis(5,myrank)
      ksh= idis(6,myrank)
      do numshell= 1,idis(8,myrank)
!ishimura
    call cpu_time(t1)
!
! AO intengral generation and first integral transformation
!
        call transmoint1(trint1a,trint1b,cmowrk,xint,ish,ksh,maxdim,nocc,mlcount,mlsize)
!ishimura
    call cpu_time(t2)
        mlcount= mlcount+mbf(ish)*mbf(ksh)
        if(numshell == idis(8,myrank)) then
          call transmoint2(trint2,trint1b,cmoocc,nocc,mlcount,mlstart,mlsize,idis)
        else
          ksh= ksh+nproc
          if(ksh > nshell) then
            do ii= 1,nshell
              ish= ish+1
              ksh= ksh-nshell+ish
              if(ksh <= nshell) exit
            enddo
          endif
!         ksh= ksh+1
!         if(ksh == nshell+1) then
!           ish= ish+1
!           ksh= 1
!         endif
          if(mlcount+mbf(ish)*mbf(ksh) > mlsize) then
!
! Second integral transformation
!
            call transmoint2(trint2,trint1b,cmoocc,nocc,mlcount,mlstart,mlsize,idis)
            mlstart= mlstart+mlcount
            mlcount= 0
          endif
        endif
!ishimura
    call cpu_time(t3)
    tr1=tr1+t2-t1
    tr2=tr2+t3-t2
      enddo

    write(*,'("Time for Tr1=",f9.3,"  Time for Tr2=",f9.3,i3)')tr1,tr2,myrank
      return
end


!-----------------------------------------------------------------------------------------
  subroutine mp2trans2(cmovir,energymo,trint2,trint3,trint4,emp2st,nocc,nvir,ncore,idis)
!-----------------------------------------------------------------------------------------
!
! Driver of third and fourth integral transformations and MP2 energy calculation
!
! In  : cmovir  (MO coefficient matrix)
!       energymo(MO energies)
!       trint2  (Second-transformed integrals)
!       nocc    (Number of active occupied MOs)
!       nvir    (Number of virtual MOs)
!       idis    (Information for parallelization)
! Work: trint3  (Third-transformed integrals)
!       trint4  (Fourth-transformed integrals)
!
      use modparallel
      use modbasis, only : nao
      use modmolecule, only : nmo
      implicit none
      integer,intent(in) :: nocc, nvir, ncore, idis(8,0:nproc-1)
      integer :: numrecv, iproc, irecv(0:nproc-1), nocc2, ncycle, icycle, myij
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: cmovir(nao*nvir), energymo(nmo)
      real(8),intent(in) :: trint2((idis(3,myrank)+idis(7,myrank)),nocc*(nocc+1)/2)
      real(8),intent(out) :: trint3(nao*nao), trint4(nao*nao), emp2st(2)
!
      numrecv= 1
      do iproc= 0,nproc-1
        irecv(iproc)= numrecv
        numrecv= numrecv+idis(3,iproc)+idis(7,iproc)
      enddo
      nocc2= nocc*(nocc+1)/2
      ncycle=(nocc2-1)/nproc+1
!
      do icycle= 1,ncycle    
        call mp2int_sendrecv(trint2,trint3,trint4,icycle,irecv,nocc2,idis)
!
        myij=(icycle-1)*nproc+1+myrank
        if(myij > nocc2) cycle
!
! Third integral transformation
!   (mi|lj) trint4[l,m] -> (ai|lj) trint3[l,a]
!
        call dgemm('N','N',nao,nvir,nao,one,trint4,nao,cmovir,nao,zero,trint3,nao)
!
! Fourth integral transformation
!   (ai|lj) trint3[l,a] -> (ai|bj) trint4[b,a]
!
        call dgemm('T','N',nvir,nvir,nao,one,cmovir,nao,trint3,nao,zero,trint4,nvir)
!
! MP2 energy calculation
!
        call calcmp2energy(trint4,energymo,emp2st,nocc,nvir,ncore,icycle)
      enddo
      return
end


!------------------------------------------------------------------------------------------
  subroutine transmoint1(trint1a,trint1b,cmowrk,xint,ish,ksh,maxdim,nocc,mlcount,mlsize)
!------------------------------------------------------------------------------------------
!
! AO intengral generation and first-quarter integral transformation
!    (mn|ls) -> (mi|ls)
!
! In  : cmowrk  (Transposed MO coefficient matrix)
!       xint    (Exchange integral matrix)
!       ish,lsh (Basis shell indices)
!       maxdim  (Maximum size of basis functions in a shell)
!       nocc    (Number of active occupied MOs)
!       mlcount (Start index of transformed AOs)
!       mlsize  (Size of trint1b)
! Out : trint1b (First-transformed integrals, [s,i,ml])
! Work: trint1a (First-transformed integrals, [i,s,ml])
!
      use modparallel
      use modbasis, only : nao, nshell, mbf, locbf
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: ish, ksh, maxdim, nocc, mlcount, mlsize
      integer :: nbfi, nbfj, nbfk, nbfl, nbfik, locbfj, locbfl, jsh, lsh, ik, i, j, k, l
      integer :: imo, kl, ij, ii, jloc, lloc
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: cmowrk(nocc,nao), xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: trint1a(nocc,maxdim,maxdim*maxdim), trint1b(nao,nocc,mlsize)
      real(8) :: twoeri(maxdim,maxdim,maxdim,maxdim)
!
      nbfi= mbf(ish)
      nbfk= mbf(ksh)
      nbfik= nbfi*nbfk
!
!$OMP parallel do schedule(dynamic,1) private(twoeri,trint1a,ij,kl,ik,ii,nbfj,nbfl,locbfj,locbfl,jloc,lloc)
      do lsh= 1,nshell
        nbfl  = mbf(lsh)
        locbfl= locbf(lsh)
        do ik= 1,nbfik
          do l= 1,nbfl
            do imo= 1,nocc
              trint1a(imo,l,ik)= zero
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
                  do imo= 1,nocc
                    trint1a(imo,l,ik)= trint1a(imo,l,ik)+twoeri(j,i,l,k)*cmowrk(imo,jloc)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
!
! Reorder of first-transformed integrals
!
        do ik= 1,nbfik
          do l= 1,nbfl
            lloc= locbfl+l
            do imo= 1,nocc
              trint1b(lloc,imo,mlcount+ik)= trint1a(imo,l,ik)
            enddo
          enddo
        enddo
      enddo
!$OMP end parallel do
      return
end
!


!----------------------------------------------------------------------------------
  subroutine transmoint2(trint2,trint1b,cmoocc,nocc,mlcount,mlstart,mlsize,idis)
!----------------------------------------------------------------------------------
!
! Second-quarter integral transformation
!    (mi|ls) -> (mi|lj)
!
! In  : trint1b (First-transformed integrals, [s,i,ml])
!       cmoocc  (MO coefficient matrix)
!       nocc    (Number of active occupied MOs)
!       mlcount (Number of transformed AOs)
!       mlstart (First index of trint2)
!       mlsize  (Size of trint1b)
!       idis    (Information for parallelization)
! Out : trint2  (Second-transformed integrals, [ml,ij])
!
      use modparallel
      use modbasis, only : nao
      implicit none
      integer,intent(in) :: nocc, mlcount, mlstart, mlsize, idis(8,0:nproc-1)
      integer :: imo, ijmo
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: trint1b(nao,nocc,mlsize), cmoocc(nao,nocc)
      real(8),intent(inout) :: trint2((idis(3,myrank)+idis(7,myrank)),nocc*(nocc+1)/2)
!
      do imo= 1,nocc
        ijmo=imo*(imo-1)/2+1
        call dgemm('T','N',mlcount,imo,nao,one,trint1b(1,imo,1),nao*nocc,cmoocc,nao,zero,&
&                  trint2(mlstart,ijmo),idis(3,myrank)+idis(7,myrank))
      enddo
      return
end


!---------------------------------------------------------------------------
  subroutine mp2int_sendrecv(trint2,trint3,trint4,icycle,irecv,nocc2,idis)
!---------------------------------------------------------------------------
!
! Send and Receive second-transformed integrals (mi|lj)
!
! In  : trint2  (Second-transformed integrals)
!       icycle  (Mp2trans2 cycle number)
!       irecv   (Numbers of receiving data)
!       nocc2   (nocc*(nocc+1)/2, nocc:Number of active occupied MOs)
!       idis    (Information for parallelization)
! Out : trint4  (Second-transformed integrals, [l,m])
! Work: trint3  (Receiving data)
!
      use modparallel
      use modbasis, only : nao, nshell, mbf, locbf
      implicit none
      integer,intent(in) :: icycle, irecv(0:nproc-1), nocc2, idis(8,0:nproc-1)
      integer :: iproc, jproc, ijstart, myij, ij, nsend, nrecv, ish, ksh, nbfi, nbfk
      integer :: locbfi, locbfk, i, k, ik, num, ii
      real(8),intent(in) :: trint2(idis(3,myrank)+idis(7,myrank),nocc2)
      real(8),intent(out) :: trint3(nao*nao), trint4(nao,nao)
!
      jproc= myrank
      ijstart=(icycle-1)*nproc+1
      myij= ijstart+myrank
!
! Send and receive second-transformed integrals
!
      do iproc= myrank+1,nproc-1
        jproc= jproc-1
        if(jproc < 0) jproc= nproc-1
        nsend= idis(3,myrank)+idis(7,myrank)
        nrecv= idis(3,jproc)+idis(7,jproc)
        ij= ijstart+iproc
        if(ij > nocc2) then
          nsend= 0
          ij   = 1
        endif
        if(myij > nocc2) nrecv=0
        call para_sendrecv(trint2(1,ij),nsend,iproc,myrank, &
&                          trint3(irecv(jproc)),nrecv,jproc,jproc,MPI_COMM_WORLD)
      enddo
!
      do iproc= 0, myrank-1
        jproc= jproc-1
        if(jproc < 0) jproc= nproc-1
        nsend= idis(3,myrank)+idis(7,myrank)
        nrecv= idis(3,jproc)+idis(7,jproc)
        ij= ijstart+iproc
        if(ij > nocc2) then
          nsend= 0
          ij   = 1
        endif
        if(myij > nocc2) nrecv=0
        call para_sendrecv(trint2(1,ij),nsend,iproc,myrank, &
&                          trint3(irecv(jproc)),nrecv,jproc,jproc,MPI_COMM_WORLD)
      enddo
!
      nsend= idis(3,myrank)+idis(7,myrank)
      if(myij <= nocc2) call dcopy(nsend,trint2(1,myij),1,trint3(irecv(myrank)),1)
!
! Reorder of received data
!
      if(myij <= nocc2) then
        ik= 0
        do iproc= 0,nproc-1
          ish= idis(1,iproc)
          ksh= idis(2,iproc)
          do num= 1,idis(4,iproc)
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
            ksh= ksh+nproc
            if(ksh > ish) then
              do ii= 1,nshell
                ksh= ksh-ish
                ish= ish+1
                if(ksh <= ish) exit
              enddo
            endif
          enddo

          ish= idis(5,iproc)
          ksh= idis(6,iproc)
          do num= 1,idis(8,iproc)
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
            ksh= ksh+nproc
            if(ksh > nshell) then
              do ii= 1,nshell
                ish= ish+1
                ksh= ksh-nshell+ish
                if(ksh <= nshell) exit
              enddo
            endif
          enddo
        enddo

!!$OMP parallel do private(nbfi,nbfk,locbfi,locbfk,ik)
!        do ish= 1,nshell
!          nbfi= mbf(ish)
!          locbfi= locbf(ish)
!          ik= locbfi*nao
!          do ksh= 1,nshell
!            nbfk= mbf(ksh)
!            locbfk= locbf(ksh)
!            do i= 1,nbfi
!              do k= 1,nbfk
!                ik= ik+1
!                trint4(locbfk+k,locbfi+i)= trint3(ik)
!              enddo
!            enddo
!          enddo
!        enddo
!!$OMP end parallel do
      endif
      return
end


!--------------------------------------------------------------------------
  subroutine calcmp2energy(trint4,energymo,emp2st,nocc,nvir,ncore,icycle)
!--------------------------------------------------------------------------
!
! Calculate MP2 energy
!
! In  : trint4   (MO integrals, [b,a])
!       energymo (MO energies)
!       nocc     (Number of active occupied Mos)
!       nvir     (Number of virtual MOs)
!       icycle   (Mp2trans2 cycle number)
!
      use modparallel
      use modmolecule, only : neleca, nmo
      implicit none
      integer,intent(in) :: nocc, nvir, ncore, icycle
      integer :: moi, moj, myij, ii, moa, mob
      real(8),parameter:: zero=0.0D+00, one=1.0D+00, two=2.0D+00
      real(8),intent(in) :: trint4(nvir,nvir), energymo(nmo)
      real(8),intent(inout) :: emp2st(2)
      real(8) :: eij, eija, eijab, etmp(2)
!
! Calculate occupied MOs, i and j
!
      myij=(icycle-1)*nproc+1+myrank
      do ii= 1,nocc
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
      do moa= 1,nvir
        eija= eij-energymo(moa+neleca)
        do mob= 1,nvir
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



