!----------------------------------------------------------
  subroutine calcrhf(h1mtrx,cmo,ortho,overlap,xint,eigen)
!----------------------------------------------------------
!
! Driver of restricted Hartree-Fock calculation
!
! In  : h1mtrx  (1-electron integral matrix)
!       ortho   (Orthogonalization matrix)
!       overlap (Overlap integral matrix)
! Out : xint    ((ij|ij) integral matrix)
!       eigen   (MO energy)
! Inout : cmo   (MO coefficient matrix)
!
      use modparallel
      use modbasis, only : nshell, nao, mtype
      use modmolecule, only : neleca, nmo
      use modscf, only : maxiter, fdiff, dconv, maxdiis, maxsoscf, diis, extrap
      use modenergy, only : enuc, escf, escfe
      use modthresh, only : threshsoscf, cutint2, threshex, threshover, thresherr
      use modprint, only : iprint
      implicit none
      integer :: nao2, nao3, nshell3, maxdim, maxfunc(0:6), iter, i, itsub, itdiis
      integer :: itextra, itsoscf, nocc, nvir
      integer :: idis(nproc,14), isize1, isize2, isize3
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, two=2.0D+00
      real(8),parameter :: small=1.0D-10
      real(8),intent(in) :: h1mtrx(nao*(nao+1)/2), ortho(nao*nao), overlap(nao*(nao+1)/2)
      real(8),intent(out) :: xint(nshell*(nshell+1)/2), eigen(nao)
      real(8),intent(inout) :: cmo(nao*nao)
      real(8),allocatable :: fock(:), fockprev(:), dmtrx(:), dmtrxprev(:), dmax(:), work(:)
      real(8),allocatable :: fockdiis(:), errdiis(:), diismtrx(:), work2(:)
      real(8),allocatable :: hstart(:), sograd(:,:), sodisp(:), sovecy(:)
      real(8) :: escfprev, diffmax, tridot, deltae, errmax, sogradmax, sodispmax
      real(8) :: time1, time2, time3, time4
      logical :: convsoscf
      data maxfunc/1,3,6,10,15,21,28/
!
      nao2= nao*nao
      nao3= nao*(nao+1)/2
      nshell3= nshell*(nshell+1)/2
      nocc= neleca
      nvir= nmo-neleca
      maxdim=maxfunc(maxval(mtype(1:nshell)))
!
! Distribute fock and error arrays for DIIS
!
      call distarray(idis,nmo,nao,nao3,nocc,nvir,nocc,nvir,nproc)
!
! Set arrays
!
      call memset(nao2+nao3*4+nshell3)
      allocate(fock(nao3),fockprev(nao3),dmtrx(nao3),dmtrxprev(nao3),dmax(nshell3),work(nao2))
!
      isize1= max(idis(myrank+1,3),idis(myrank+1,7)*nao,maxdiis)
      isize2= idis(myrank+1,3)
      isize3= idis(myrank+1,5)
      call memset(isize3*maxdiis+isize1+isize2*maxdiis+maxdiis*(maxdiis+1)/2)
      allocate(fockdiis(isize3*maxdiis), errdiis(isize2*maxdiis), &
&              diismtrx(maxdiis*(maxdiis+1)/2), work2(isize1))
      if(.not.diis) then
       call memset(nocc*nvir*3*maxsoscf)
       allocate(hstart(nocc*nvir),sograd(nocc*nvir,maxsoscf),sodisp(nocc*nvir*maxsoscf), &
&               sovecy(nocc*nvir*(maxsoscf-1)))
      endif
!
      escfprev= zero
      itdiis =0
      itextra=0
      itsoscf=0
      convsoscf=.false.
!
! Calculate initial density matrix
!
        call calcdmtrx(cmo,dmtrx,work,nao,neleca)
!
! Set 1-electron Hamiltonian
!
      call dcopy(nao3,h1mtrx,1,fockprev,1)
      call dcopy(nao3,dmtrx,1,dmtrxprev,1)
!
! Calculate (ij|ij) integrals
!
      call calcschwarzeri(xint,work,maxdim)
!
      if(master) then
        write(*,'(1x,74("-"))')
        write(*,'("   Hartree-Fock calculation")')
        write(*,'(1x,74("-"))')
        write(*,'("   DIIS = ",l1,",   SOSCF = ",l1)')diis,.not.diis
        write(*,'("   Dconv      =",1p,d9.2,",  MaxIter    = ",i9  ,",  ThreshSOSCF=",d9.2)') &
&                    dconv, maxiter, threshsoscf
        write(*,'("   Cutint2    =",1p,d9.2,",  ThreshEx   = ",d9.2,",  ThreshOver =",d9.2)') &
&                    cutint2, threshex, threshover
        write(*,'(1x,74("-"))')
        write(*,'(" ====================")')
        write(*,'("    SCF Iteration")')
        write(*,'(" ====================")')
        if(diis) then
          write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                      "Delta Density     DIIS Error")')
        else
          write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                      "Delta Density    Orbital Grad")')
        endif
      endif
!
! Start SCF iteration
!
      do iter= 1,maxiter
!
! Calculate maximum density matrix elements
!
        call calcdmax(dmtrx,dmax,work)
!
! Calculate two-electron integrals and Fock matrix
!
        call cpu_time(time1)
        call formfock(fock,work,dmtrx,dmax,xint,maxdim)
        call dscal(nao3,half,fock,1)
        do i= 1,nao
          fock(i*(i+1)/2)= fock(i*(i+1)/2)*two
        enddo
        call cpu_time(time2)
!
! Form full Fock matrix
!
        call focksum(fock,fockprev,nao3)
!
! Calculate SCF ENERGY
!
        escfe=(tridot(dmtrxprev,fock,nao)+tridot(dmtrxprev,h1mtrx,nao))*half
        escf = escfe+enuc
        deltae = escf-escfprev
        escfprev= escf
!
        if(diis) then
!
! DIIS interpolation
!
          call calcdiiserr(fock,dmtrxprev,overlap,ortho,cmo,work,work2,errmax,idis,nao,nmo)
          if(((itdiis /= 0).or.(errmax <= thresherr)).and.(errmax > small))then
            itdiis= itdiis+1
            call calcdiis(fock,errdiis,fockdiis,diismtrx,cmo,work2,idis,itdiis,nao,maxdiis)
          endif
!
! Extrapolate Fock matrix
!
          if(extrap.and.itdiis == 0) &
&           call fockextrap(fock,fockdiis,work,cmo,dmtrx,idis,itextra,nao,maxdiis)
!
! Diagonalize Fock matrix
!
          call diagfock(fock,work,ortho,cmo,work2,eigen,idis)
        else
!
! Approximated Second-order SCF method
!
          if((itsoscf == 0).or.(convsoscf)) then
            call diagfock(fock,work,ortho,cmo,work2,eigen,idis)
            sogradmax= zero
            itsoscf= itsoscf+1
          else
            call expand(fock,work,nao)
            call soscfgrad(work,work2,sograd(1,itsoscf),cmo,nocc,nvir,sogradmax,idis,nao,1)
            if(sogradmax <= threshsoscf) then
              if(itsoscf == 1) call soscfinith(hstart,eigen,nocc,nvir,nao)
              call soscfnewh(hstart,sograd,sodisp,sovecy,nocc,nvir,itsoscf,maxsoscf,sodispmax)
              call soscfupdate(cmo,sodisp,work,work2,nocc,nvir,itsoscf,maxsoscf,idis, &
&                              nao,nmo,sodispmax)
              itsoscf= itsoscf+1
            else
              call diagfock(fock,work,ortho,cmo,work2,eigen,idis)
            endif
          endif
        endif
        call cpu_time(time3)
!
! Copy previous density matrix and calculate new density matrix
!
        call calcdmtrx(cmo,dmtrx,work,nao,neleca)
        call ddiff(dmtrx,dmtrxprev,work,nao3,diffmax)
        if(diis) then
          if(extrap.and.(itdiis==0)) then
            itsub= itextra
          else
            itsub= itdiis
          endif
        else
          itsub= itsoscf
        endif
        if(master) then
          if(diis) then
            write(*,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,errmax
          else
            write(*,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,sogradmax
          endif
        endif
!
! Check SCF convergence
!
        if(diis) then
          if(diffmax.lt.dconv) exit
          if(itdiis >= maxdiis) itdiis= 0
        else
          if((diffmax.lt.dconv).and.(convsoscf)) exit
          if((diffmax.lt.dconv).and.(itsoscf==1)) exit
          if((diffmax.lt.dconv).and.(.not.convsoscf)) convsoscf=.true.
          if(itsoscf >= maxsoscf) itsoscf= 0
        endif
        if(iter.eq.maxiter) then
          if(master) then
            write(*,'(" Not Converged.")')
            call iabort
          endif
        endif
        call dcopy(nao3,dmtrx,1,dmtrxprev,1)
        call dcopy(nao3,work,1,dmtrx,1)
        call cpu_time(time4)
        if(master.and.(iprint >= 2)) write(*,'(10x,6f8.3)')time2-time1,time3-time2,time4-time3
      enddo
!
      if(master) then
        write(*,'(" -----------------------------------------")')
        write(*,'("    SCF Converged.")')
        write(*,'("    RHF Energy = ",f17.9," a.u.")')escf
        write(*,'(" -----------------------------------------"/)')
      endif
!
! Unset arrays
!
      if(.not.diis) then
       call memunset(nocc*nvir*3*maxsoscf)
       deallocate(hstart,sograd,sodisp,sovecy)
      endif
      deallocate(fockdiis,errdiis,diismtrx,work2)
      call memunset(isize3*maxdiis+isize1+isize2*maxdiis+maxdiis*(maxdiis+1)/2)
      deallocate(fock,fockprev,dmtrx,dmtrxprev,dmax,work)
      call memunset(nao2+nao3*4+nshell3)
      return
end


!------------------------------------------------------------
  subroutine diagfock(fock,work,ortho,cmo,work2,eigen,idis)
!------------------------------------------------------------
!
! Driver of Fock matrix diagonalization
!
      use modparallel
      use modbasis, only : nao
      use modmolecule, only : nmo
      implicit none
      integer,intent(in) :: idis(nproc,14)
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: fock(nao*(nao+1)/2), ortho(nao,nao)
      real(8),intent(out) :: work(nao,nao), cmo(nao*nao), work2(*), eigen(nao)
!
! Canonicalize Fock matrix
!
        call expand(fock,work,nao)
        call canonicalizep(work,ortho,cmo,work2,nao,nmo,idis)
!
! Diagonalize canonicalized matrix
!
        call diag('V','U',nmo,work,nao,eigen)
!
! Backtransform to AO basis
!
      if(idis(myrank+1,1) /= 0) &
&       call dgemm('N','N',nao,idis(myrank+1,1),nmo,one,ortho,nao, &
&                  work(1,idis(myrank+1,2)+1),nao,zero,work2,nao)
      call para_allgatherv(work2,idis(myrank+1,3),cmo,idis(1,3),idis(1,4),MPI_COMM_WORLD)
!
      return
end


!-------------------------------------------------------------
  subroutine formfock(focktotal,fock,dmtrx,dmax,xint,maxdim)
!-------------------------------------------------------------
!
! Driver of Fock matrix formation from two-electron intgrals
!
      use modparallel
      use modbasis, only : nshell, nao
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: maxdim
      integer :: ish, jsh, ksh, lsh, ij, kl, ik, il, jk, jl
      integer :: ii, jj, kk, kstart
      integer(8) :: ncount, icount
      real(8),parameter :: zero=0.0D+00, four=4.0D+00
      real(8),intent(in) :: dmtrx(nao*(nao+1)/2), dmax(nshell*(nshell+1)/2)
      real(8),intent(in) :: xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: focktotal(nao*(nao+1)/2), fock(nao*(nao+1)/2)
      real(8) :: xijkl, denmax, twoeri(maxdim**4), denmax1
      integer :: last, ltmp(nshell),lnum,ll
!
      fock(:)= zero
!
      ncount= 0
      ncount= ncount+(2*nshell**3+3*nshell**2+nshell)/6+myrank
!
!$OMP parallel do schedule(dynamic,1) &
!$OMP private(jsh,ksh,lsh,ij,kl,ik,il,jk,jl,xijkl,denmax1,denmax,twoeri,ii,jj,kk, &
!$OMP icount,kstart,last,ltmp,lnum,ll) reduction(+:fock)
      do ish= nshell,1,-1
        ii= ish*(ish-1)/2
        icount=ncount-(2*ish*ish*ish-3*ish*ish+ish)/6
        do jsh= 1,ish
          ij= ii+jsh
          jj= jsh*(jsh-1)/2
          kstart=mod(icount-ish*(jsh-1),nproc)+1
          do ksh= kstart,ish,nproc
            kk= ksh*(ksh-1)/2
            ik= ii+ksh
            jk= jj+ksh
            if(jsh.lt.ksh) jk= kk+jsh
            denmax1=max(four*dmax(ij),dmax(ik),dmax(jk))
            last= ksh
            if(ish == ksh) last= jsh
            ll=min(jsh,ksh)
            lnum=0
!           do lsh= 1,ksh
            do lsh= 1,ll
              kl= kk+lsh
              il= ii+lsh
              jl= jj+lsh
              xijkl=xint(ij)*xint(kl)
              denmax=max(denmax1,four*dmax(kl),dmax(il),dmax(jl))
              if(xijkl*denmax.ge.cutint2) then
                lnum=lnum+1
                ltmp(lnum)=lsh
              endif
            enddo
            do lsh= ll+1,last
              kl= kk+lsh
              il= ii+lsh
              jl= lsh*(lsh-1)/2+jsh
              xijkl=xint(ij)*xint(kl)
              denmax=max(denmax1,four*dmax(kl),dmax(il),dmax(jl))
              if(xijkl*denmax.ge.cutint2) then
                lnum=lnum+1
                ltmp(lnum)=lsh
              endif
            enddo
            do lsh= 1,lnum
              call calc2eri(twoeri,ish,jsh,ksh,ltmp(lsh),maxdim)
              call fockeri(fock,dmtrx,twoeri,ish,jsh,ksh,ltmp(lsh),maxdim)
            enddo
          enddo
        enddo
      enddo
!$OMP end parallel do
!
      call para_allreduce(fock,focktotal,nao*(nao+1)/2,MPI_SUM,MPI_COMM_WORLD)
      return
end


!---------------------------------------------------------------
  subroutine fockeri(fock,dmtrx,twoeri,ish,jsh,ksh,lsh,maxdim)
!---------------------------------------------------------------
!
! Form Fock matrix from two-electron intgrals
!
      use modbasis, only : nshell, nao, mbf, locbf
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: ish, jsh, ksh, lsh, maxdim
      integer :: nbfi, nbfj, nbfk, nbfl
      integer :: locbfi, locbfj, locbfk, locbfl, jmax, lmax, i, j, k, l, ij, kl
      integer :: nij, nkl, nik, nil, njk, njl
      integer :: iloc, jloc, kloc, lloc, iloc2, jloc2, kloc2, lloc2
      real(8),parameter :: half=0.5D+00, four=4.0D+00
      real(8),intent(in) :: dmtrx(nao*(nao+1)/2), twoeri(maxdim,maxdim,maxdim,maxdim)
      real(8),intent(inout) :: fock(nao*(nao+1)/2)
      real(8) :: val, val4
      logical :: ieqj, keql, ieqk, jeql, ikandjl, ijorkl
!
      nbfi = mbf(ish)
      nbfj = mbf(jsh)
      nbfk = mbf(ksh)
      nbfl = mbf(lsh)
      locbfi= locbf(ish)
      locbfj= locbf(jsh)
      locbfk= locbf(ksh)
      locbfl= locbf(lsh)
!
      ieqj= ish.eq.jsh
      keql= ksh.eq.lsh
      ieqk= ish.eq.ksh
      jeql= jsh.eq.lsh
      ikandjl= ieqk.and.jeql
      ijorkl= ieqj.or.keql
      jmax= nbfj
      lmax= nbfl
      ij= 0
!
      do i= 1,nbfi
        iloc= locbfi+i
        iloc2= iloc*(iloc-1)/2
        if(ieqj) jmax= i
        do j= 1,jmax
          jloc= locbfj+j
          jloc2= jloc*(jloc-1)/2
          ij= ij+1
          kl= 0
  kloop:  do k= 1,nbfk
            kloc= locbfk+k
            kloc2= kloc*(kloc-1)/2
            if(keql) lmax= k
            do l= 1,lmax
              kl= kl+1
              if(ikandjl.and.kl.gt.ij) exit kloop
              val= twoeri(l,k,j,i)
              if(abs(val).lt.cutint2) cycle
              lloc= locbfl+l
              nij= iloc2+jloc
              nkl= kloc2+lloc
              nik= iloc2+kloc
              nil= iloc2+lloc
              njk= jloc2+kloc
              njl= jloc2+lloc
              if(jloc.lt.kloc) njk= kloc2+jloc
              if(jloc.lt.lloc) njl= lloc*(lloc-1)/2+jloc
              if(ieqk) then
                if(iloc.lt.kloc) then
                  lloc2= lloc*(lloc-1)/2
                  nij= kloc2+lloc
                  nkl= iloc2+jloc
                  nik= kloc2+iloc
                  nil= kloc2+jloc
                  njk= lloc2+iloc
                  njl= lloc2+jloc
                  if(lloc.lt.iloc) njk= iloc2+lloc
                  if(lloc.lt.jloc) njl= jloc2+lloc
                elseif(iloc.eq.kloc.and.jloc.eq.lloc) then
                  val= val*half
                endif
              endif
              if(ijorkl) then
                if(iloc.eq.jloc) val= val*half
                if(kloc.eq.lloc) val= val*half
              endif
              val4= val*four
!
              fock(nij)= fock(nij)+val4*dmtrx(nkl)
              fock(nkl)= fock(nkl)+val4*dmtrx(nij)
              fock(nik)= fock(nik)-val *dmtrx(njl)
              fock(nil)= fock(nil)-val *dmtrx(njk)
              fock(njk)= fock(njk)-val *dmtrx(nil)
              fock(njl)= fock(njl)-val *dmtrx(nik)
            enddo
          enddo kloop
        enddo
      enddo
      return
end


!----------------------------------------------------------------------------
  subroutine formrdftfock(focktotal,fock,dmtrx,dmax,xint,maxdim,hfexchange)
!----------------------------------------------------------------------------
!
! Driver of DFT Fock matrix formation from two-electron intgrals
!
      use modparallel
      use modbasis, only : nshell, nao
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: maxdim
      integer :: ish, jsh, ksh, lsh, ij, kl, ik, il, jk, jl
      integer :: ii, jj, kk, kstart
      integer(8) :: ncount, icount
      real(8),parameter :: zero=0.0D+00, four=4.0D+00
      real(8),intent(in) :: dmtrx(nao*(nao+1)/2), dmax(nshell*(nshell+1)/2)
      real(8),intent(in) :: xint(nshell*(nshell+1)/2), hfexchange
      real(8),intent(out) :: focktotal(nao*(nao+1)/2), fock(nao*(nao+1)/2)
      real(8) :: xijkl, denmax, twoeri(maxdim**4),denmax1
      integer :: last, ltmp(nshell), lnum, ll
!
      fock(:)= zero
!
      ncount=(2*nshell**3+3*nshell**2+nshell)/6+myrank
!
!$OMP parallel do schedule(dynamic,1) &
!$OMP private(jsh,ksh,lsh,ij,kl,ik,il,jk,jl,xijkl,denmax1,denmax,twoeri,ii,jj,kk, &
!$OMP icount,kstart,last,ltmp,lnum,ll) reduction(+:fock)
      do ish= nshell,1,-1
        ii= ish*(ish-1)/2
        icount=ncount-(2*ish*ish*ish-3*ish*ish+ish)/6
        do jsh= 1,ish
          ij= ii+jsh
          jj= jsh*(jsh-1)/2
          kstart=mod(icount-ish*(jsh-1),nproc)+1
          do ksh= kstart,ish,nproc
            kk= ksh*(ksh-1)/2
            ik= ii+ksh
            jk= jj+ksh
            if(jsh.lt.ksh) jk= kk+jsh
            denmax1=max(four*dmax(ij),dmax(ik),dmax(jk))
            last= ksh
            if(ish == ksh) last= jsh
            ll=min(jsh,ksh)
            lnum=0
!           do lsh= 1,ksh
            do lsh= 1,ll
              kl= kk+lsh
              il= ii+lsh
              jl= jj+lsh
              xijkl=xint(ij)*xint(kl)
              denmax=max(denmax1,four*dmax(kl),dmax(il),dmax(jl))
              if(xijkl*denmax.ge.cutint2) then
                lnum=lnum+1
                ltmp(lnum)=lsh
              endif
            enddo
            do lsh= ll+1,last
              kl= kk+lsh
              il= ii+lsh
              jl= lsh*(lsh-1)/2+jsh
              xijkl=xint(ij)*xint(kl)
              denmax=max(denmax1,four*dmax(kl),dmax(il),dmax(jl))
              if(xijkl*denmax.ge.cutint2) then
                lnum=lnum+1
                ltmp(lnum)=lsh
              endif
            enddo
            do lsh= 1,lnum
              call calc2eri(twoeri,ish,jsh,ksh,ltmp(lsh),maxdim)
              call dftfockeri(fock,dmtrx,twoeri,ish,jsh,ksh,ltmp(lsh),maxdim,hfexchange)
            enddo
          enddo
        enddo
      enddo
!$OMP end parallel do
!
      call para_allreduce(fock,focktotal,nao*(nao+1)/2,MPI_SUM,MPI_COMM_WORLD)
      return
end


!-----------------------------------------------------------------------------
  subroutine dftfockeri(fock,dmtrx,twoeri,ish,jsh,ksh,lsh,maxdim,hfexchange)
!-----------------------------------------------------------------------------
!
! Form DFT Fock matrix from two-electron intgrals
!
      use modbasis, only : nshell, nao, mbf, locbf
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: ish, jsh, ksh, lsh, maxdim
      integer :: nbfi, nbfj, nbfk, nbfl
      integer :: locbfi, locbfj, locbfk, locbfl, jmax, lmax, i, j, k, l, ij, kl
      integer :: ii, jj, kk, ll, ii2, jj2, kk2, n, nij, nkl, nik, nil, njk, njl
      integer :: iloc, jloc, kloc, iloc2, jloc2, kloc2
      real(8),parameter :: half=0.5D+00, four=4.0D+00
      real(8),intent(in) :: dmtrx(nao*(nao+1)/2), twoeri(maxdim,maxdim,maxdim,maxdim)
      real(8),intent(in) :: hfexchange
      real(8),intent(inout) :: fock(nao*(nao+1)/2)
      real(8) :: val, val4
      logical :: ieqj, keql, ieqk, jeql, ikandjl, ijorkl
!
      nbfi = mbf(ish)
      nbfj = mbf(jsh)
      nbfk = mbf(ksh)
      nbfl = mbf(lsh)
      locbfi= locbf(ish)
      locbfj= locbf(jsh)
      locbfk= locbf(ksh)
      locbfl= locbf(lsh)
!
      ieqj= ish.eq.jsh
      keql= ksh.eq.lsh
      ieqk= ish.eq.ksh
      jeql= jsh.eq.lsh
      ikandjl= ieqk.and.jeql
      ijorkl= ieqj.or.keql
      jmax= nbfj
      lmax= nbfl
      ij= 0
!
      do i= 1,nbfi
        iloc=locbfi+i
        iloc2= iloc*(iloc-1)/2
        if(ieqj) jmax= i
        do j= 1,jmax
          jloc=locbfj+j
          jloc2= jloc*(jloc-1)/2
          ij= ij+1
          kl= 0
  kloop:  do k= 1,nbfk
            kloc=locbfk+k
            kloc2= kloc*(kloc-1)/2
            if(keql) lmax= k
            do l= 1,lmax
              kl= kl+1
              if(ikandjl.and.kl.gt.ij) exit kloop
              val= twoeri(l,k,j,i)
              if(abs(val).lt.cutint2) cycle
              ii= iloc
              jj= jloc
              kk= kloc
              ll= locbfl+l
              ii2= iloc2
              jj2= jloc2
              kk2= kloc2
              if(ieqk) then
                if(ii.lt.kk) then
                  n = ii
                  ii= kk
                  kk= n
                  n = jj
                  jj= ll
                  ll= n
                  n  = ii2
                  ii2= kk2
                  kk2= n
                  jj2= ll*(ll-1)/2
                elseif(ii.eq.kk.and.jj.eq.ll) then
                  val= val*half
                endif
              endif
              if(ijorkl) then
                if(ii.eq.jj) val= val*half
                if(kk.eq.ll) val= val*half
              endif
              nij= ii2+jj
              nkl= kk2+ll
              nik= ii2+kk
              nil= ii2+ll
              njk= jj2+kk
              njl= jj2+ll
              if(jj.lt.kk) njk= kk2+jj
              if(jj.lt.ll) njl= ll*(ll-1)/2+jj
              val4= val*four
!
              fock(nij)= fock(nij)+val4*dmtrx(nkl)
              fock(nkl)= fock(nkl)+val4*dmtrx(nij)
              fock(nik)= fock(nik)-val *dmtrx(njl)*hfexchange
              fock(nil)= fock(nil)-val *dmtrx(njk)*hfexchange
              fock(njk)= fock(njk)-val *dmtrx(nil)*hfexchange
              fock(njl)= fock(njl)-val *dmtrx(nik)*hfexchange
            enddo
          enddo kloop
        enddo
      enddo
      return
end



!-----------------------------------------------------------
  subroutine calcrdft(h1mtrx,cmo,ortho,overlap,xint,eigen)
!-----------------------------------------------------------
!
! Driver of restricted DFT calculation
!
! In  : h1mtrx  (1-electron integral matrix)
!       ortho   (Orthogonalization matrix)
!       overlap (Overlap integral matrix)
! Out : xint    ((ij|ij) integral matrix)
!       eigen   (MO energy)
! Inout : cmo   (MO coefficient matrix)
!
      use modparallel
      use moddft, only : idft, nrad, nleb, hfexchange
      use modatom, only : atomrad
      use modbasis, only : nshell, nao, mtype
      use modmolecule, only : neleca, nmo, natom, numatomic
      use modscf, only : maxiter, fdiff, dconv, maxdiis, maxsoscf, diis, extrap
      use modenergy, only : enuc, escf, escfe
      use modthresh, only : threshsoscf, cutint2, threshex, threshover, thresherr, &
&                           threshrho, threshdfock
      use modprint, only : iprint
      use modunit, only : tobohr
      implicit none
      integer :: nao2, nao3, nshell3, maxdim, maxfunc(0:6), iter, i, itsub, itdiis
      integer :: itextra, itsoscf, nocc, nvir
      integer :: idis(nproc,14), isize1, isize2, isize3, iatom
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, two=2.0D+00
      real(8),parameter :: small=1.0D-10
      real(8),intent(in) :: h1mtrx(nao*(nao+1)/2), ortho(nao*nao), overlap(nao*(nao+1)/2)
      real(8),intent(out) :: xint(nshell*(nshell+1)/2), eigen(nao)
      real(8),intent(inout) :: cmo(nao*nao)
      real(8),allocatable :: fock(:), fockprev(:), dmtrx(:), dmtrxprev(:), dmax(:), work(:)
      real(8),allocatable :: fockdiis(:), errdiis(:), diismtrx(:), work2(:)
      real(8),allocatable :: hstart(:), sograd(:,:), sodisp(:), sovecy(:)
      real(8),allocatable :: fockd(:), rad(:), atomvec(:), surface(:),  radpt(:), angpt(:)
      real(8),allocatable :: ptweight(:), xyzpt(:), rsqrd(:), vao(:), vmo(:)
      real(8) :: escfprev, diffmax, tridot, deltae, errmax, sogradmax, sodispmax
      real(8) :: edft, totalelec
      real(8) :: time1, time2, time3, time4
      logical :: convsoscf
      data maxfunc/1,3,6,10,15,21,28/
!
      nao2= nao*nao
      nao3= nao*(nao+1)/2
      nshell3= nshell*(nshell+1)/2
      nocc= neleca
      nvir= nmo-neleca
      maxdim=maxfunc(maxval(mtype(1:nshell)))
!
! Distribute fock and error arrays for DIIS
!
      call distarray(idis,nmo,nao,nao3,nocc,nvir,nocc,nvir,nproc)
!
! Set arrays
!
      call memset(nao2+nao3*4+nshell3)
      allocate(fock(nao3),fockprev(nao3),dmtrx(nao3),dmtrxprev(nao3),dmax(nshell3),work(nao2))
      call memset(nao3+natom*5+6*natom*natom+2*nrad+4*nleb+natom*nrad*nleb+4*nao+4*nocc)
      allocate(fockd(nao3),rad(natom),atomvec(5*natom*natom),surface(natom*natom), &
&              radpt(2*nrad),angpt(4*nleb),ptweight(natom*nrad*nleb),xyzpt(3*natom), &
&              rsqrd(natom),vao(4*nao),vmo(4*nocc))
!
      isize1= max(idis(myrank+1,3),idis(myrank+1,7)*nao,3*natom,nao,maxdiis)
      isize2= idis(myrank+1,3)
      isize3=idis(myrank+1,5)
      call memset(isize3*maxdiis+isize1+isize2*maxdiis+maxdiis*(maxdiis+1)/2)
      allocate(fockdiis(isize3*maxdiis), errdiis(isize2*maxdiis), &
&              diismtrx(maxdiis*(maxdiis+1)/2), work2(isize1))
      if(.not.diis) then
       call memset(nocc*nvir*3*maxsoscf)
       allocate(hstart(nocc*nvir),sograd(nocc*nvir,maxsoscf),sodisp(nocc*nvir*maxsoscf), &
&               sovecy(nocc*nvir*(maxsoscf-1)))
      endif
!
      escfprev= zero
      itdiis =0
      itextra=0
      itsoscf=0
      convsoscf=.false.
!
! Calculate DFT information
!
      call calcatomvec(atomvec,surface)
      call calcradpt(radpt,nrad)
      call calclebpt(angpt,nleb)
      do iatom= 1,natom
        rad(iatom)= atomrad(numatomic(iatom))*tobohr
      enddo
      call calcgridweight(ptweight,rad,radpt,angpt,atomvec,surface,xyzpt,work2)
!
! Calculate initial density matrix
!
        call calcdmtrx(cmo,dmtrx,work,nao,neleca)
!
! Set 1-electron Hamiltonian
!
      call dcopy(nao3,h1mtrx,1,fockprev,1)
      call dcopy(nao3,dmtrx,1,dmtrxprev,1)
!
! Calculate (ij|ij) integrals
!
      call calcschwarzeri(xint,work,maxdim)
!
      if(master) then
        write(*,'(1x,74("-"))')
        write(*,'("   DFT calculation")')
        write(*,'(1x,74("-"))')
        write(*,'("   DIIS = ",l1,",   SOSCF = ",l1)')diis,.not.diis
        write(*,'("   Dconv      =",1p,d9.2,",  MaxIter    = ",i9  ,",  ThreshSOSCF=",d9.2)') &
&                    dconv, maxiter, threshsoscf
        write(*,'("   Cutint2    =",1p,d9.2,",  ThreshEx   = ",d9.2,",  ThreshOver =",d9.2)') &
&                    cutint2, threshex, threshover
        write(*,'("   Nrad       =",1p,i9  ,",  Nleb       = ",i9  ,",  ThreshRho  =",d9.2)') &
&                    nrad, nleb, threshrho
        write(*,'("   ThreshDfock=",1p,d9.2)') threshdfock
        write(*,'(1x,74("-"))')
        write(*,'(" ====================")')
        write(*,'("    SCF Iteration")')
        write(*,'(" ====================")')
        if(diis) then
          write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                      "Delta Density     DIIS Error")')
        else
          write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                      "Delta Density    Orbital Grad")')
        endif
      endif
!
! Start SCF iteration
!
      do iter= 1,maxiter
!
! Calculate maximum density matrix elements
!
        call calcdmax(dmtrx,dmax,work)
        call cpu_time(time1)
!
! Calculate exchange-correlation terms
!
        call formrfockexcor(fockd,fock,edft,totalelec,cmo,atomvec,radpt,angpt, &
&                           rad,ptweight,vao,vmo,xyzpt,rsqrd,work,work2,idft)
!
! Calculate two-electron integrals
!
        call formrdftfock(fock,work,dmtrx,dmax,xint,maxdim,hfexchange)
        call dscal(nao3,half,fock,1)
        do i= 1,nao
          fock(i*(i+1)/2)= fock(i*(i+1)/2)*two
        enddo
        call cpu_time(time2)
!
! Form Fock matrix
!
        call focksum(fock,fockprev,nao3)
!
! Calculate DFT ENERGY
!
        escfe=(tridot(dmtrxprev,fock,nao)+tridot(dmtrxprev,h1mtrx,nao))*half+edft
        escf = escfe+enuc
        deltae = escf-escfprev
        escfprev= escf
!
! Form full Fock matrix
!
        call daxpy(nao3,one,fockd,1,fock,1)
!
        if(diis) then
!
! DIIS interpolation
!
          call calcdiiserr(fock,dmtrxprev,overlap,ortho,cmo,work,work2,errmax,idis,nao,nmo)
          if(((itdiis /= 0).or.(errmax <= thresherr)).and.(errmax > small))then
            itdiis= itdiis+1
            call calcdiis(fock,errdiis,fockdiis,diismtrx,cmo,work2,idis,itdiis,nao,maxdiis)
          endif
!
! Extrapolate Fock matrix
!
          if(extrap.and.itdiis == 0) &
&           call fockextrap(fock,fockdiis,work,cmo,dmtrx,idis,itextra,nao,maxdiis)
!
! Diagonalize Fock matrix
!
          call diagfock(fock,work,ortho,cmo,work2,eigen,idis)
        else
!
! Approximated Second-order SCF method
!
          if((itsoscf == 0).or.(convsoscf)) then
            call diagfock(fock,work,ortho,cmo,work2,eigen,idis)
            sogradmax= zero
            itsoscf= itsoscf+1
          else
            call expand(fock,work,nao)
            call soscfgrad(work,work2,sograd(1,itsoscf),cmo,nocc,nvir,sogradmax,idis,nao,1)
            if(sogradmax <= threshsoscf) then
              if(itsoscf == 1) call soscfinith(hstart,eigen,nocc,nvir,nao)
              call soscfnewh(hstart,sograd,sodisp,sovecy,nocc,nvir,itsoscf,maxsoscf,sodispmax)
              call soscfupdate(cmo,sodisp,work,work2,nocc,nvir,itsoscf,maxsoscf,idis, &
&                              nao,nmo,sodispmax)
              itsoscf= itsoscf+1
            else
              call diagfock(fock,work,ortho,cmo,work2,eigen,idis)
            endif
          endif
        endif
        call cpu_time(time3)
!
! Copy previous density matrix and calculate new density matrix
!
        call calcdmtrx(cmo,dmtrx,work,nao,neleca)
        call ddiff(dmtrx,dmtrxprev,work,nao3,diffmax)
        if(diis) then
          if(extrap.and.(itdiis==0)) then
            itsub= itextra
          else
            itsub= itdiis
          endif
        else
          itsub= itsoscf
        endif
        if(master) then
          if(diis) then
            write(*,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,errmax
          else
            write(*,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,sogradmax
          endif
        endif
!
! Check SCF convergence
!
        if(diis) then
          if(diffmax.lt.dconv) exit
          if(itdiis >= maxdiis) itdiis= 0
        else
          if((diffmax.lt.dconv).and.(convsoscf)) exit
          if((diffmax.lt.dconv).and.(itsoscf==1)) exit
          if((diffmax.lt.dconv).and.(.not.convsoscf)) convsoscf=.true.
          if(itsoscf >= maxsoscf) itsoscf= 0
        endif
        if(iter.eq.maxiter) then
          if(master) then
            write(*,'(" Not Converged.")')
            call iabort
          endif
        endif
        call dcopy(nao3,dmtrx,1,dmtrxprev,1)
        call dcopy(nao3,work,1,dmtrx,1)
        call cpu_time(time4)
        if(master.and.(iprint >= 2)) write(*,'(10x,6f8.3)')time2-time1,time3-time2,time4-time3
      enddo
!
      if(master) then
        write(*,'(" -----------------------------------------------------------")')
        write(*,'("    SCF Converged.")')
        write(*,'("    DFT Energy = ",f17.9," a.u.")')escf
        write(*,'("    Exchange + Correlation energy = ",f17.9," a.u.")')edft
        write(*,'("    Number of electrons           = ",f17.9)')totalelec
        write(*,'(" -----------------------------------------------------------"/)')
      endif
!
! Unset arrays
!
      if(.not.diis) then
       call memunset(nocc*nvir*3*maxsoscf)
       deallocate(hstart,sograd,sodisp,sovecy)
      endif
      deallocate(fockd,rad,atomvec,surface,radpt,angpt,ptweight,xyzpt,rsqrd,vao,vmo)
      call memunset(nao3+natom*5+6*natom*natom+2*nrad+4*nleb+natom*nrad*nleb+4*nao+4*nocc)
      deallocate(fockdiis,errdiis,diismtrx,work2)
      call memunset(isize3*maxdiis+isize1+isize2*maxdiis+maxdiis*(maxdiis+1)/2)
      deallocate(fock,fockprev,dmtrx,dmtrxprev,dmax,work)
      call memunset(nao2+nao3*4+nshell3)
      return
end


!------------------------------------------------------------------------
  subroutine calcuhf(h1mtrx,cmoa,cmob,ortho,overlap,xint,eigena,eigenb)
!------------------------------------------------------------------------
!
! Driver of unrestricted Hartree-Fock calculation
!
! In  : h1mtrx  (1-electron integral matrix)
!       ortho   (Orthogonalization matrix)
!       overlap (Overlap integral matrix)
! Out : xint    ((ij|ij) integral matrix)
!       eigena  (Alpha MO energy)
!       eigenb  (Beta MO energy)
! Inout : cmoa  (Alpha MO coefficient matrix)
!         cmob  (Beta MO coefficient matrix)
!
      use modparallel
      use modbasis, only : nshell, nao, mtype
      use modmolecule, only : neleca, nelecb, nmo
      use modscf, only : maxiter, fdiff, dconv, maxdiis, maxsoscf, diis, extrap
      use modenergy, only : enuc, escf, escfe
      use modthresh, only : threshsoscf, cutint2, threshex, threshover, thresherr
      use modprint, only : iprint
      implicit none
      integer :: nao3, nshell3, maxdim, maxfunc(0:6), iter, i, itsub, itdiis
      integer :: itextra, itsoscf, nocca, nvira, noccb, nvirb
      integer :: idis(nproc,14), isize1, isize2, isize3
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, two=2.0D+00
      real(8),parameter :: small=1.0D-10
      real(8),intent(in) :: h1mtrx(nao*(nao+1)/2), ortho(nao*nao), overlap(nao*(nao+1)/2)
      real(8),intent(out) :: xint(nshell*(nshell+1)/2), eigena(nao), eigenb(nao)
      real(8),intent(inout) :: cmoa(nao*nao), cmob(nao*nao)
      real(8),allocatable :: focka(:), fockb(:), fockpreva(:), fockprevb(:), dmtrxa(:)
      real(8),allocatable :: dmtrxb(:), dmtrxpreva(:), dmtrxprevb(:), dmax(:), work(:)
      real(8),allocatable :: fockdiisa(:), fockdiisb(:), errdiisa(:), errdiisb(:)
      real(8),allocatable :: diismtrx(:), work2(:), work3(:)
      real(8),allocatable :: hstarta(:), hstartb(:), sograda(:,:), sogradb(:,:)
      real(8),allocatable :: sodispa(:), sodispb(:), sovecya(:), sovecyb(:)
      real(8) :: escfprev, diffmax, diffmaxa, diffmaxb, tridot, deltae
      real(8) :: errmax, errmaxa, errmaxb, sogradmax, sogradmaxa, sogradmaxb, sodispmax
      real(8) :: s2, sz
      real(8) :: time1, time2, time3, time4
      logical :: convsoscf
      data maxfunc/1,3,6,10,15,21,28/
!
      nao3= nao*(nao+1)/2
      nshell3= nshell*(nshell+1)/2
      nocca= neleca
      nvira= nmo-neleca
      noccb= nelecb
      nvirb= nmo-nelecb
      maxdim=maxfunc(maxval(mtype(1:nshell)))
!
! Distribute fock and error arrays for DIIS
!
      call distarray(idis,nmo,nao,nao3,nocca,nvira,noccb,nvirb,nproc)
!
! Set arrays
!
      call memset(nao3*10+nshell3)
      allocate(focka(nao3),fockb(nao3),fockpreva(nao3),fockprevb(nao3),dmtrxa(nao3), &
&              dmtrxb(nao3),dmtrxpreva(nao3),dmtrxprevb(nao3),dmax(nshell3),work(nao3*2))
!
      isize1= max(idis(myrank+1,3),idis(myrank+1,7)*nao,idis(myrank+1,11)*nao,maxdiis)
      isize2=idis(myrank+1,3)
      isize3=idis(myrank+1,5)
      call memset(isize3*maxdiis*2+isize1+isize2*maxdiis*2+maxdiis*(maxdiis+1)/2)
      allocate(fockdiisa(isize3*maxdiis),fockdiisb(isize3*maxdiis),errdiisa(isize2*maxdiis),&
&              errdiisb(isize2*maxdiis),diismtrx(maxdiis*(maxdiis+1)/2),work2(isize1))
      if(.not.diis) then
       call memset(nocca*nvira*3*maxsoscf+noccb*nvirb*3*maxsoscf)
       allocate(hstarta(nocca*nvira),hstartb(noccb*nvirb), &
&               sograda(nocca*nvira,maxsoscf),sogradb(noccb*nvirb,maxsoscf), &
&               sodispa(nocca*nvira*maxsoscf),sodispb(noccb*nvirb*maxsoscf), &
&               sovecya(nocca*nvira*(maxsoscf-1)),sovecyb(noccb*nvirb*(maxsoscf-1)))
      endif
!
      escfprev= zero
      itdiis =0
      itextra=0
      itsoscf=0
      convsoscf=.false.
!
! Calculate initial density matrix
!
        call calcudmtrx(cmoa,cmob,dmtrxa,dmtrxb,work,nao,neleca,nelecb)
!
! Set 1-electron Hamiltonian
!
      call dcopy(nao3,h1mtrx,1,fockpreva,1)
      call dcopy(nao3,h1mtrx,1,fockprevb,1)
      call dcopy(nao3,dmtrxa,1,dmtrxpreva,1)
      call dcopy(nao3,dmtrxb,1,dmtrxprevb,1)
!
! Calculate (ij|ij) integrals
!
      call calcschwarzeri(xint,work,maxdim)
!
      if(master) then
        write(*,'(1x,74("-"))')
        write(*,'("   Unrestricted Hartree-Fock calculation")')
        write(*,'(1x,74("-"))')
        write(*,'("   DIIS = ",l1,",   SOSCF = ",l1)')diis,.not.diis
        write(*,'("   Dconv      =",1p,d9.2,",  MaxIter    = ",i9  ,",  ThreshSOSCF=",d9.2)') &
&                    dconv, maxiter, threshsoscf
        write(*,'("   Cutint2    =",1p,d9.2,",  ThreshEx   = ",d9.2,",  ThreshOver =",d9.2)') &
&                    cutint2, threshex, threshover
        write(*,'(1x,74("-"))')
        write(*,'(" ====================")')
        write(*,'("    SCF Iteration")')
        write(*,'(" ====================")')
        if(diis) then
          write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                      "Delta Density     DIIS Error")')
        else
          write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                      "Delta Density    Orbital Grad")')
        endif
      endif
!
! Start SCF iteration
!
      do iter= 1,maxiter
!
! Calculate maximum density matrix elements
!
        call calcudmax(dmtrxa,dmtrxb,dmax,work)
!
! Calculate two-electron integrals and Fock matrix
!
        call cpu_time(time1)
        call formufock(focka,fockb,work,dmtrxa,dmtrxb,dmax,xint,maxdim)
        call dscal(nao3,half,focka,1)
        call dscal(nao3,half,fockb,1)
        do i= 1,nao
          focka(i*(i+1)/2)= focka(i*(i+1)/2)*two
          fockb(i*(i+1)/2)= fockb(i*(i+1)/2)*two
        enddo
        call cpu_time(time2)
!
! Form full Fock matrix
!
        call focksum(focka,fockpreva,nao3)
        call focksum(fockb,fockprevb,nao3)
!
! Calculate SCF ENERGY
!
        escfe=(tridot(dmtrxpreva,focka,nao)+tridot(dmtrxpreva,h1mtrx,nao))*half &
&            +(tridot(dmtrxprevb,fockb,nao)+tridot(dmtrxprevb,h1mtrx,nao))*half
        escf = escfe+enuc
        deltae = escf-escfprev
        escfprev= escf
!
        if(diis) then
!
! DIIS interpolation
!
          call calcdiiserr(focka,dmtrxpreva,overlap,ortho,cmoa,work,work2,errmaxa,idis,nao,nmo)
          call calcdiiserr(fockb,dmtrxprevb,overlap,ortho,cmob,work,work2,errmaxb,idis,nao,nmo)
          errmax= max(errmaxa,errmaxb)
          if(((itdiis /= 0).or.(errmax <= thresherr)).and.(errmax > small))then
            itdiis= itdiis+1
            call calcudiis(focka,fockb,errdiisa,errdiisb,fockdiisa,fockdiisb, &
&                          diismtrx,cmoa,cmob,work2,idis,itdiis,nao,maxdiis)
          endif
!
! Extrapolate Fock matrix
!
          if(extrap.and.itdiis == 0) then
            call fockextrap(focka,fockdiisa,work,cmoa,dmtrxa,idis,itextra,nao,maxdiis)
            call fockextrap(fockb,fockdiisb,work,cmob,dmtrxb,idis,itextra,nao,maxdiis)
          endif
!
! Diagonalize Fock matrix
!
          call diagfock(focka,work,ortho,cmoa,work2,eigena,idis)
          call diagfock(fockb,work,ortho,cmob,work2,eigenb,idis)
        else
!
! Approximated Second-order SCF method
!
          if((itsoscf == 0).or.(convsoscf)) then
            call diagfock(focka,work,ortho,cmoa,work2,eigena,idis)
            call diagfock(fockb,work,ortho,cmob,work2,eigenb,idis)
            sogradmax= zero
            itsoscf= itsoscf+1
          else
            call expand(focka,work,nao)
            call soscfgrad(work,work2,sograda(1,itsoscf),cmoa,nocca,nvira,sogradmaxa,idis,nao,1)
            if(noccb /= 0) then
              call expand(fockb,work,nao)
              call soscfgrad(work,work2,sogradb(1,itsoscf),cmob,noccb,nvirb,sogradmaxb,idis,nao,2)
            else
              sogradmaxb= zero
            endif
            sogradmax= max(sogradmaxa,sogradmaxb)
            if(sogradmax <= threshsoscf) then
              if(itsoscf == 1) then
                call soscfinith(hstarta,eigena,nocca,nvira,nao)
                if(noccb /= 0) call soscfinith(hstartb,eigenb,noccb,nvirb,nao)
              endif
              call soscfunewh(hstarta,hstartb,sograda,sogradb,sodispa,sodispb,sovecya,sovecyb, &
&                             nocca,noccb,nvira,nvirb,itsoscf,maxsoscf,sodispmax)
              call soscfupdate(cmoa,sodispa,work,work2,nocca,nvira,itsoscf,maxsoscf,idis, &
&                              nao,nmo,sodispmax)
              if(noccb /= 0 ) &
&               call soscfupdate(cmob,sodispb,work,work2,noccb,nvirb,itsoscf,maxsoscf,idis, &
&                                nao,nmo,sodispmax)
              itsoscf= itsoscf+1
            else
              call diagfock(focka,work,ortho,cmoa,work2,eigena,idis)
              call diagfock(fockb,work,ortho,cmob,work2,eigenb,idis)
            endif
          endif
        endif
        call cpu_time(time3)
!
! Copy previous density matrix and calculate new density matrix
!
        call calcudmtrx(cmoa,cmob,dmtrxa,dmtrxb,work,nao,neleca,nelecb)
        call ddiff(dmtrxa,dmtrxpreva,work(1),nao3,diffmaxa)
        call ddiff(dmtrxb,dmtrxprevb,work(nao3+1),nao3,diffmaxb)
        diffmax= diffmaxa+diffmaxb
        if(diis) then
          if(extrap.and.(itdiis==0)) then
            itsub= itextra
          else
            itsub= itdiis
          endif
        else
          itsub= itsoscf
        endif
        if(master) then
          if(diis) then
            write(*,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,errmax
          else
            write(*,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,sogradmax
          endif
        endif
!
! Check SCF convergence
!
        if(diis) then
          if(diffmax.lt.dconv) exit
          if(itdiis >= maxdiis) itdiis= 0
        else
          if((diffmax.lt.dconv).and.(convsoscf)) exit
          if((diffmax.lt.dconv).and.(itsoscf==1)) exit
          if((diffmax.lt.dconv).and.(.not.convsoscf)) convsoscf=.true.
          if(itsoscf >= maxsoscf) itsoscf= 0
        endif
        if(iter.eq.maxiter) then
          if(master) then
            write(*,'(" SCF Not Converged.")')
            call iabort
          endif
        endif
        call dcopy(nao3,dmtrxa,1,dmtrxpreva,1)
        call dcopy(nao3,dmtrxb,1,dmtrxprevb,1)
        call dcopy(nao3,work(1),1,dmtrxa,1)
        call dcopy(nao3,work(nao3+1),1,dmtrxb,1)
        call cpu_time(time4)
        if(master.and.(iprint >= 2)) write(*,'(10x,6f8.3)')time2-time1,time3-time2,time4-time3
      enddo
!
      if(master) then
        write(*,'(" -----------------------------------------")')
        write(*,'("    SCF Converged.")')
        write(*,'("    UHF Energy = ",f17.9," a.u.")')escf
        write(*,'(" -----------------------------------------"/)')
      endif
!
! Unset arrays
!
      if(.not.diis) then
       call memunset(nocca*nvira*3*maxsoscf+noccb*nvirb*3*maxsoscf)
       deallocate(hstarta,hstartb, &
&                 sograda,sogradb, &
&                 sodispa,sodispb, &
&                 sovecya,sovecyb)
      endif
      deallocate(fockdiisa,fockdiisb,errdiisa, &
&                errdiisb,diismtrx,work2)
      call memunset(isize3*maxdiis*2+isize1+isize2*maxdiis*2+maxdiis*(maxdiis+1)/2)
!
! Set arrays
!
      call memset(nao*nao+idis(myrank+1,3))
      allocate(work2(nao*nao),work3(idis(myrank+1,3)))
!
! Calculate spin expectation values
!
      call calcspin(sz,s2,dmtrxa,dmtrxb,overlap,work,work2,work3,neleca,nelecb,nao,idis)
!
      if(master) then
        write(*,'(" -------------------------------")')
        write(*,'("    Sz =",f7.3)')sz
        write(*,'("    S-squared =",f7.3)')s2
        write(*,'(" -------------------------------"/)')
      endif
!
! Unset arrays
!
      deallocate(work2,work3)
      call memunset(nao*nao+idis(myrank+1,3))
!
      deallocate(focka,fockb,fockpreva,fockprevb,dmtrxa, &
&                dmtrxb,dmtrxpreva,dmtrxprevb,dmax,work)
      call memunset(nao3*10+nshell3)
      return
end


!-------------------------------------------------------------------------
  subroutine formufock(fock1,fock2,fock3,dmtrxa,dmtrxb,dmax,xint,maxdim)
!-------------------------------------------------------------------------
!
! Driver of Fock matrix formation from two-electron intgrals
!
! In  : dmtrxa (Alpha density matrix)
!       dmtrxb (Beta density matrix)
! Out : fock1  (Alpha Fock matrix)
!       fock2  (Beta Fock matrix)
!       fock3  (Work space)
!
      use modparallel
      use modbasis, only : nshell, nao
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: maxdim
      integer :: ish, jsh, ksh, lsh, ij, kl, ik, il, jk, jl
      integer :: ii, jj, kk, kstart
      integer(8) :: ncount, icount
      real(8),parameter :: zero=0.0D+00, four=4.0D+00
      real(8),intent(in) :: dmtrxa(nao*(nao+1)/2), dmtrxb(nao*(nao+1)/2)
      real(8),intent(in) :: dmax(nshell*(nshell+1)/2), xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: fock1(nao*(nao+1)/2), fock2(nao*(nao+1)/2), fock3(nao*(nao+1)/2)
      real(8) :: xijkl, denmax, twoeri(maxdim**4), denmax1
      integer :: last, ltmp(nshell),lnum,ll
!
      fock2(:)= zero
      fock3(:)= zero
      ncount= 0
      ncount= ncount+(2*nshell**3+3*nshell**2+nshell)/6+myrank
!
!$OMP parallel do schedule(dynamic,1) &
!$OMP private(jsh,ksh,lsh,ij,kl,ik,il,jk,jl,xijkl,denmax1,denmax,twoeri,ii,jj,kk, &
!$OMP icount,kstart,last,ltmp,lnum,ll) reduction(+:fock2,fock3)
      do ish= nshell,1,-1
        ii= ish*(ish-1)/2
        icount=ncount-(2*ish*ish*ish-3*ish*ish+ish)/6
        do jsh= 1,ish
          ij= ii+jsh
          jj= jsh*(jsh-1)/2
          kstart=mod(icount-ish*(jsh-1),nproc)+1
          do ksh= kstart,ish,nproc
            kk= ksh*(ksh-1)/2
            ik= ii+ksh
            jk= jj+ksh
            if(jsh.lt.ksh) jk= kk+jsh
            denmax1=max(four*dmax(ij),dmax(ik),dmax(jk))
            last= ksh
            if(ish == ksh) last= jsh
            ll=min(jsh,ksh)
            lnum=0
!           do lsh= 1,ksh
            do lsh= 1,ll
              kl= kk+lsh
              il= ii+lsh
              jl= jj+lsh
              xijkl=xint(ij)*xint(kl)
              denmax=max(denmax1,four*dmax(kl),dmax(il),dmax(jl))
              if(xijkl*denmax.ge.cutint2) then
                lnum=lnum+1
                ltmp(lnum)=lsh
              endif
            enddo
            do lsh= ll+1,last
              kl= kk+lsh
              il= ii+lsh
              jl= lsh*(lsh-1)/2+jsh
              xijkl=xint(ij)*xint(kl)
              denmax=max(denmax1,four*dmax(kl),dmax(il),dmax(jl))
              if(xijkl*denmax.ge.cutint2) then
                lnum=lnum+1
                ltmp(lnum)=lsh
              endif
            enddo
            do lsh= 1,lnum
              call calc2eri(twoeri,ish,jsh,ksh,ltmp(lsh),maxdim)
              call ufockeri(fock2,fock3,dmtrxa,dmtrxb,twoeri,ish,jsh,ksh,ltmp(lsh),maxdim)
            enddo
          enddo
        enddo
      enddo
!$OMP end parallel do
!
      call para_allreduce(fock2,fock1,nao*(nao+1)/2,MPI_SUM,MPI_COMM_WORLD)
      call para_allreduce(fock3,fock2,nao*(nao+1)/2,MPI_SUM,MPI_COMM_WORLD)
      return
end


!--------------------------------------------------------------------------------
  subroutine ufockeri(focka,fockb,dmtrxa,dmtrxb,twoeri,ish,jsh,ksh,lsh,maxdim)
!--------------------------------------------------------------------------------
!
! Form unrestricted Fock matrix from two-electron intgrals
!
      use modbasis, only : nshell, nao, mbf, locbf
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: ish, jsh, ksh, lsh, maxdim
      integer :: nbfi, nbfj, nbfk, nbfl
      integer :: locbfi, locbfj, locbfk, locbfl, jmax, lmax, i, j, k, l, ij, kl
      integer :: nij, nkl, nik, nil, njk, njl
      integer :: iloc, jloc, kloc, lloc, iloc2, jloc2, kloc2, lloc2
      real(8),parameter :: half=0.5D+00, two=2.0D+00, four=4.0D+00
      real(8),intent(in) :: dmtrxa(nao*(nao+1)/2), dmtrxb(nao*(nao+1)/2)
      real(8),intent(in) :: twoeri(maxdim,maxdim,maxdim,maxdim)
      real(8),intent(inout) :: focka(nao*(nao+1)/2), fockb(nao*(nao+1)/2)
      real(8) :: val, val2, val4
      logical :: ieqj, keql, ieqk, jeql, ikandjl, ijorkl
!
      nbfi = mbf(ish)
      nbfj = mbf(jsh)
      nbfk = mbf(ksh)
      nbfl = mbf(lsh)
      locbfi= locbf(ish)
      locbfj= locbf(jsh)
      locbfk= locbf(ksh)
      locbfl= locbf(lsh)
!
      ieqj= ish.eq.jsh
      keql= ksh.eq.lsh
      ieqk= ish.eq.ksh
      jeql= jsh.eq.lsh
      ikandjl= ieqk.and.jeql
      ijorkl= ieqj.or.keql
      jmax= nbfj
      lmax= nbfl
      ij= 0
!
      do i= 1,nbfi
        iloc= locbfi+i
        iloc2= iloc*(iloc-1)/2
        if(ieqj) jmax= i
        do j= 1,jmax
          jloc= locbfj+j
          jloc2= jloc*(jloc-1)/2
          ij= ij+1
          kl= 0
  kloop:  do k= 1,nbfk
            kloc= locbfk+k
            kloc2= kloc*(kloc-1)/2
            if(keql) lmax= k
            do l= 1,lmax
              kl= kl+1
              if(ikandjl.and.kl.gt.ij) exit kloop
              val= twoeri(l,k,j,i)
              if(abs(val).lt.cutint2) cycle
              lloc= locbfl+l
              nij= iloc2+jloc
              nkl= kloc2+lloc
              nik= iloc2+kloc
              nil= iloc2+lloc
              njk= jloc2+kloc
              njl= jloc2+lloc
              if(jloc.lt.kloc) njk= kloc2+jloc
              if(jloc.lt.lloc) njl= lloc*(lloc-1)/2+jloc
              if(ieqk) then
                if(iloc.lt.kloc) then
                  lloc2= lloc*(lloc-1)/2
                  nij= kloc2+lloc
                  nkl= iloc2+jloc
                  nik= kloc2+iloc
                  nil= kloc2+jloc
                  njk= lloc2+iloc
                  njl= lloc2+jloc
                  if(lloc.lt.iloc) njk= iloc2+lloc
                  if(lloc.lt.jloc) njl= jloc2+lloc
                elseif(iloc.eq.kloc.and.jloc.eq.lloc) then
                  val= val*half
                endif
              endif
              if(ijorkl) then
                if(iloc.eq.jloc) val= val*half
                if(kloc.eq.lloc) val= val*half
              endif
              val2= val*two
              val4= val*four
!
              focka(nij)= focka(nij)+val4*(dmtrxa(nkl)+dmtrxb(nkl))
              fockb(nij)= fockb(nij)+val4*(dmtrxa(nkl)+dmtrxb(nkl))
              focka(nkl)= focka(nkl)+val4*(dmtrxa(nij)+dmtrxb(nij))
              fockb(nkl)= fockb(nkl)+val4*(dmtrxa(nij)+dmtrxb(nij))
!
              focka(nik)= focka(nik)-val2*dmtrxa(njl)
              focka(nil)= focka(nil)-val2*dmtrxa(njk)
              focka(njk)= focka(njk)-val2*dmtrxa(nil)
              focka(njl)= focka(njl)-val2*dmtrxa(nik)
!
              fockb(nik)= fockb(nik)-val2*dmtrxb(njl)
              fockb(nil)= fockb(nil)-val2*dmtrxb(njk)
              fockb(njk)= fockb(njk)-val2*dmtrxb(nil)
              fockb(njl)= fockb(njl)-val2*dmtrxb(nik)
            enddo
          enddo kloop
        enddo
      enddo
      return
end


!-------------------------------------------------------------------------
  subroutine calcudft(h1mtrx,cmoa,cmob,ortho,overlap,xint,eigena,eigenb)
!-------------------------------------------------------------------------
!
! Driver of unrestricted DFT calculation
!
! In  : h1mtrx  (1-electron integral matrix)
!       ortho   (Orthogonalization matrix)
!       overlap (Overlap integral matrix)
! Out : xint    ((ij|ij) integral matrix)
!       eigena  (Alpha MO energy)
!       eigenb  (Beta MO energy)
! Inout : cmoa  (Alpha MO coefficient matrix)
!         cmob  (Beta MO coefficient matrix)
!
      use modparallel
      use moddft, only : idft, nrad, nleb, hfexchange
      use modatom, only : atomrad
      use modbasis, only : nshell, nao, mtype
      use modmolecule, only : neleca, nelecb, nmo, natom, numatomic
      use modscf, only : maxiter, fdiff, dconv, maxdiis, maxsoscf, diis, extrap
      use modenergy, only : enuc, escf, escfe
      use modthresh, only : threshsoscf, cutint2, threshex, threshover, thresherr, &
&                           threshrho, threshdfock
      use modprint, only : iprint
      use modunit, only : tobohr
      implicit none
      integer :: nao3, nshell3, maxdim, maxfunc(0:6), iter, i, itsub, itdiis
      integer :: itextra, itsoscf, nocca, nvira, noccb, nvirb
      integer :: idis(nproc,14), isize1, isize2, isize3, iatom
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, two=2.0D+00
      real(8),parameter :: small=1.0D-10
      real(8),intent(in) :: h1mtrx(nao*(nao+1)/2), ortho(nao*nao), overlap(nao*(nao+1)/2)
      real(8),intent(out) :: xint(nshell*(nshell+1)/2), eigena(nao), eigenb(nao)
      real(8),intent(inout) :: cmoa(nao*nao), cmob(nao*nao)
      real(8),allocatable :: focka(:), fockb(:), fockpreva(:), fockprevb(:), dmtrxa(:)
      real(8),allocatable :: dmtrxb(:), dmtrxpreva(:), dmtrxprevb(:), dmax(:), work(:)
      real(8),allocatable :: fockdiisa(:), fockdiisb(:), errdiisa(:), errdiisb(:)
      real(8),allocatable :: diismtrx(:), work2(:), work3(:)
      real(8),allocatable :: hstarta(:), hstartb(:), sograda(:,:), sogradb(:,:)
      real(8),allocatable :: sodispa(:), sodispb(:), sovecya(:), sovecyb(:)
      real(8),allocatable :: fockda(:), fockdb(:), rad(:), atomvec(:), surface(:)
      real(8),allocatable :: radpt(:), angpt(:), ptweight(:), xyzpt(:), rsqrd(:)
      real(8),allocatable :: vao(:), vmoa(:), vmob(:)
      real(8) :: escfprev, diffmax, diffmaxa, diffmaxb, tridot, deltae
      real(8) :: errmax, errmaxa, errmaxb, sogradmax, sogradmaxa, sogradmaxb, sodispmax
      real(8) :: s2, sz, edft, totalelec
      real(8) :: time1, time2, time3, time4
      logical :: convsoscf
      data maxfunc/1,3,6,10,15,21,28/
!
      nao3= nao*(nao+1)/2
      nshell3= nshell*(nshell+1)/2
      nocca= neleca
      nvira= nmo-neleca
      noccb= nelecb
      nvirb= nmo-nelecb
      maxdim=maxfunc(maxval(mtype(1:nshell)))
!
! Distribute fock and error arrays for DIIS
!
      call distarray(idis,nmo,nao,nao3,nocca,nvira,noccb,nvirb,nproc)
!
! Set arrays
!
      call memset(nao3*10+nshell3)
      allocate(focka(nao3),fockb(nao3),fockpreva(nao3),fockprevb(nao3),dmtrxa(nao3), &
&              dmtrxb(nao3),dmtrxpreva(nao3),dmtrxprevb(nao3),dmax(nshell3),work(nao3*2))
      call memset(nao3*2+natom*5+6*natom*natom+2*nrad+4*nleb+natom*nrad*nleb+4*nao &
&                 +4*nocca+4*noccb)
      allocate(fockda(nao3),fockdb(nao3),rad(natom),atomvec(5*natom*natom), &
&              surface(natom*natom),radpt(2*nrad),angpt(4*nleb),ptweight(natom*nrad*nleb), &
&              xyzpt(3*natom),rsqrd(natom),vao(4*nao),vmoa(4*nocca),vmob(4*noccb))
!
      isize1= max(idis(myrank+1,3),idis(myrank+1,7)*nao,idis(myrank+1,11)*nao,3*natom,2*nao,maxdiis)
      isize2=idis(myrank+1,3)
      isize3=idis(myrank+1,5)
      call memset(isize3*maxdiis*2+isize1+isize2*maxdiis*2+maxdiis*(maxdiis+1)/2)
      allocate(fockdiisa(isize3*maxdiis),fockdiisb(isize3*maxdiis),errdiisa(isize2*maxdiis),&
&              errdiisb(isize2*maxdiis),diismtrx(maxdiis*(maxdiis+1)/2),work2(isize1))
      if(.not.diis) then
       call memset(nocca*nvira*3*maxsoscf+noccb*nvirb*3*maxsoscf)
       allocate(hstarta(nocca*nvira),hstartb(noccb*nvirb), &
&               sograda(nocca*nvira,maxsoscf),sogradb(noccb*nvirb,maxsoscf), &
&               sodispa(nocca*nvira*maxsoscf),sodispb(noccb*nvirb*maxsoscf), &
&               sovecya(nocca*nvira*(maxsoscf-1)),sovecyb(noccb*nvirb*(maxsoscf-1)))
      endif
!
      escfprev= zero
      itdiis =0
      itextra=0
      itsoscf=0
      convsoscf=.false.
!
! Calculate DFT information
!
      call calcatomvec(atomvec,surface)
      call calcradpt(radpt,nrad)
      call calclebpt(angpt,nleb)
      do iatom= 1,natom
        rad(iatom)= atomrad(numatomic(iatom))*tobohr
      enddo
      call calcgridweight(ptweight,rad,radpt,angpt,atomvec,surface,xyzpt,work2)
!
! Calculate initial density matrix
!
        call calcudmtrx(cmoa,cmob,dmtrxa,dmtrxb,work,nao,neleca,nelecb)
!
! Set 1-electron Hamiltonian
!
      call dcopy(nao3,h1mtrx,1,fockpreva,1)
      call dcopy(nao3,h1mtrx,1,fockprevb,1)
      call dcopy(nao3,dmtrxa,1,dmtrxpreva,1)
      call dcopy(nao3,dmtrxb,1,dmtrxprevb,1)
!
! Calculate (ij|ij) integrals
!
      call calcschwarzeri(xint,work,maxdim)
!
      if(master) then
        write(*,'(1x,74("-"))')
        write(*,'("   Unrestricted DFT calculation")')
        write(*,'(1x,74("-"))')
        write(*,'("   DIIS = ",l1,",   SOSCF = ",l1)')diis,.not.diis
        write(*,'("   Dconv      =",1p,d9.2,",  MaxIter    = ",i9  ,",  ThreshSOSCF=",d9.2)') &
&                    dconv, maxiter, threshsoscf
        write(*,'("   Cutint2    =",1p,d9.2,",  ThreshEx   = ",d9.2,",  ThreshOver =",d9.2)') &
&                    cutint2, threshex, threshover
        write(*,'("   Nrad       =",1p,i9  ,",  Nleb       = ",i9  ,",  ThreshRho  =",d9.2)') &
&                    nrad, nleb, threshrho
        write(*,'("   ThreshDfock=",1p,d9.2)') threshdfock
        write(*,'(1x,74("-"))')
        write(*,'(" ====================")')
        write(*,'("    SCF Iteration")')
        write(*,'(" ====================")')
        if(diis) then
          write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                      "Delta Density     DIIS Error")')
        else
          write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                      "Delta Density    Orbital Grad")')
        endif
      endif
!
! Start SCF iteration
!
      do iter= 1,maxiter
!
! Calculate maximum density matrix elements
!
        call calcudmax(dmtrxa,dmtrxb,dmax,work)
        call cpu_time(time1)
!
! Calculate exchange-correlation terms
!
        call formufockexcor(fockda,fockdb,focka,edft,totalelec,cmoa,cmob,atomvec,&
&                           radpt,angpt,rad,ptweight,vao,vmoa,vmob,xyzpt,rsqrd, &
&                           work,work(neleca*nao+1),work2,idft)
!
! Calculate two-electron integrals
!
        call formudftfock(focka,fockb,work,dmtrxa,dmtrxb,dmax,xint,maxdim,hfexchange)
        call dscal(nao3,half,focka,1)
        call dscal(nao3,half,fockb,1)
        do i= 1,nao
          focka(i*(i+1)/2)= focka(i*(i+1)/2)*two
          fockb(i*(i+1)/2)= fockb(i*(i+1)/2)*two
        enddo
        call cpu_time(time2)
!
! Form Fock matrix
!
        call focksum(focka,fockpreva,nao3)
        call focksum(fockb,fockprevb,nao3)
!
! Calculate DFT ENERGY
!
        escfe=(tridot(dmtrxpreva,focka,nao)+tridot(dmtrxpreva,h1mtrx,nao))*half &
&            +(tridot(dmtrxprevb,fockb,nao)+tridot(dmtrxprevb,h1mtrx,nao))*half+ edft
        escf = escfe+enuc
        deltae = escf-escfprev
        escfprev= escf
!
! Form full Fock matrix
!
        call daxpy(nao3,one,fockda,1,focka,1)
        call daxpy(nao3,one,fockdb,1,fockb,1)
!
        if(diis) then
!
! DIIS interpolation
!
          call calcdiiserr(focka,dmtrxpreva,overlap,ortho,cmoa,work,work2,errmaxa,idis,nao,nmo)
          call calcdiiserr(fockb,dmtrxprevb,overlap,ortho,cmob,work,work2,errmaxb,idis,nao,nmo)
          errmax= max(errmaxa,errmaxb)
          if(((itdiis /= 0).or.(errmax <= thresherr)).and.(errmax > small))then
            itdiis= itdiis+1
            call calcudiis(focka,fockb,errdiisa,errdiisb,fockdiisa,fockdiisb, &
&                          diismtrx,cmoa,cmob,work2,idis,itdiis,nao,maxdiis)
          endif
!
! Extrapolate Fock matrix
!
          if(extrap.and.itdiis == 0) then
            call fockextrap(focka,fockdiisa,work,cmoa,dmtrxa,idis,itextra,nao,maxdiis)
            call fockextrap(fockb,fockdiisb,work,cmob,dmtrxb,idis,itextra,nao,maxdiis)
          endif
!
! Diagonalize Fock matrix
!
          call diagfock(focka,work,ortho,cmoa,work2,eigena,idis)
          call diagfock(fockb,work,ortho,cmob,work2,eigenb,idis)
        else
!
! Approximated Second-order SCF method
!
          if((itsoscf == 0).or.(convsoscf)) then
            call diagfock(focka,work,ortho,cmoa,work2,eigena,idis)
            call diagfock(fockb,work,ortho,cmob,work2,eigenb,idis)
            sogradmax= zero
            itsoscf= itsoscf+1
          else
            call expand(focka,work,nao)
            call soscfgrad(work,work2,sograda(1,itsoscf),cmoa,nocca,nvira,sogradmaxa,idis,nao,1)
            if(noccb /= 0) then
              call expand(fockb,work,nao)
              call soscfgrad(work,work2,sogradb(1,itsoscf),cmob,noccb,nvirb,sogradmaxb,idis,nao,2)
            else
              sogradmaxb= zero
            endif
            sogradmax= max(sogradmaxa,sogradmaxb)
            if(sogradmax <= threshsoscf) then
              if(itsoscf == 1) then
                call soscfinith(hstarta,eigena,nocca,nvira,nao)
                if(noccb /= 0) call soscfinith(hstartb,eigenb,noccb,nvirb,nao)
              endif
              call soscfunewh(hstarta,hstartb,sograda,sogradb,sodispa,sodispb,sovecya,sovecyb, &
&                             nocca,noccb,nvira,nvirb,itsoscf,maxsoscf,sodispmax)
              call soscfupdate(cmoa,sodispa,work,work2,nocca,nvira,itsoscf,maxsoscf,idis, &
&                              nao,nmo,sodispmax)
              if(noccb /= 0 ) &
&               call soscfupdate(cmob,sodispb,work,work2,noccb,nvirb,itsoscf,maxsoscf,idis, &
&                                nao,nmo,sodispmax)
              itsoscf= itsoscf+1
            else
              call diagfock(focka,work,ortho,cmoa,work2,eigena,idis)
              call diagfock(fockb,work,ortho,cmob,work2,eigenb,idis)
            endif
          endif
        endif
        call cpu_time(time3)
!
! Copy previous density matrix and calculate new density matrix
!
        call calcudmtrx(cmoa,cmob,dmtrxa,dmtrxb,work,nao,neleca,nelecb)
        call ddiff(dmtrxa,dmtrxpreva,work(1),nao3,diffmaxa)
        call ddiff(dmtrxb,dmtrxprevb,work(nao3+1),nao3,diffmaxb)
        diffmax= diffmaxa+diffmaxb
        if(diis) then
          if(extrap.and.(itdiis==0)) then
            itsub= itextra
          else
            itsub= itdiis
          endif
        else
          itsub= itsoscf
        endif
        if(master) then
          if(diis) then
            write(*,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,errmax
          else
            write(*,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,sogradmax
          endif
        endif
!
! Check SCF convergence
!
        if(diis) then
          if(diffmax.lt.dconv) exit
          if(itdiis >= maxdiis) itdiis= 0
        else
          if((diffmax.lt.dconv).and.(convsoscf)) exit
          if((diffmax.lt.dconv).and.(itsoscf==1)) exit
          if((diffmax.lt.dconv).and.(.not.convsoscf)) convsoscf=.true.
          if(itsoscf >= maxsoscf) itsoscf= 0
        endif
        if(iter.eq.maxiter) then
          if(master) then
            write(*,'(" SCF Not Converged.")')
            call iabort
          endif
        endif
        call dcopy(nao3,dmtrxa,1,dmtrxpreva,1)
        call dcopy(nao3,dmtrxb,1,dmtrxprevb,1)
        call dcopy(nao3,work(1),1,dmtrxa,1)
        call dcopy(nao3,work(nao3+1),1,dmtrxb,1)
        call cpu_time(time4)
        if(master.and.(iprint >= 2)) write(*,'(10x,6f8.3)')time2-time1,time3-time2,time4-time3
      enddo
!
      if(master) then
        write(*,'(" -----------------------------------------")')
        write(*,'("    SCF Converged.")')
        write(*,'("    DFT Energy = ",f17.9," a.u.")')escf
        write(*,'("    Exchange + Correlation energy = ",f17.9," a.u.")')edft
        write(*,'("    Number of electrons           = ",f17.9)')totalelec
        write(*,'(" -----------------------------------------"/)')
      endif
!
! Unset arrays
!
      if(.not.diis) then
       call memunset(nocca*nvira*3*maxsoscf+noccb*nvirb*3*maxsoscf)
       deallocate(hstarta,hstartb, &
&                 sograda,sogradb, &
&                 sodispa,sodispb, &
&                 sovecya,sovecyb)
      endif
      deallocate(fockdiisa,fockdiisb,errdiisa, &
&                errdiisb,diismtrx,work2)
      call memunset(isize3*maxdiis*2+isize1+isize2*maxdiis*2+maxdiis*(maxdiis+1)/2)
      deallocate(fockda,fockdb,rad,atomvec, &
&                surface,radpt,angpt,ptweight, &
&                xyzpt,rsqrd,vao,vmoa,vmob)
      call memunset(nao3*2+natom*5+6*natom*natom+2*nrad+4*nleb+natom*nrad*nleb+4*nao &
&                   +4*nocca+4*noccb)
!
! Set arrays
!
      call memset(nao*nao+idis(myrank+1,3))
      allocate(work2(nao*nao),work3(idis(myrank+1,3)))
!
! Calculate spin expectation values
!
      call calcspin(sz,s2,dmtrxa,dmtrxb,overlap,work,work2,work3,neleca,nelecb,nao,idis)
!
      if(master) then
        write(*,'(" -------------------------------")')
        write(*,'("    Sz =",f7.3)')sz
        write(*,'("    S-squared =",f7.3)')s2
        write(*,'(" -------------------------------"/)')
      endif
!
! Unset arrays
!
      deallocate(work2,work3)
      call memunset(nao*nao+idis(myrank+1,3))
!
      deallocate(focka,fockb,fockpreva,fockprevb,dmtrxa, &
&                dmtrxb,dmtrxpreva,dmtrxprevb,dmax,work)
      call memunset(nao3*10+nshell3)
      return
end


!---------------------------------------------------------------------------------------
  subroutine formudftfock(fock1,fock2,fock3,dmtrxa,dmtrxb,dmax,xint,maxdim,hfexchange)
!---------------------------------------------------------------------------------------
!
! Driver of DFT Fock matrix formation from two-electron intgrals
!
! In  : dmtrxa (Alpha density matrix)
!       dmtrxb (Beta density matrix)
! Out : fock1  (Alpha Fock matrix)
!       fock2  (Beta Fock matrix)
!       fock3  (Work space)
!
      use modparallel
      use modbasis, only : nshell, nao
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: maxdim
      integer :: ish, jsh, ksh, lsh, ij, kl, ik, il, jk, jl
      integer :: ii, jj, kk, kstart
      integer(8) :: ncount, icount
      real(8),parameter :: zero=0.0D+00, four=4.0D+00
      real(8),intent(in) :: dmtrxa(nao*(nao+1)/2), dmtrxb(nao*(nao+1)/2)
      real(8),intent(in) :: dmax(nshell*(nshell+1)/2), xint(nshell*(nshell+1)/2), hfexchange
      real(8),intent(out) :: fock1(nao*(nao+1)/2), fock2(nao*(nao+1)/2), fock3(nao*(nao+1)/2)
      real(8) :: xijkl, denmax, twoeri(maxdim**4),denmax1
      integer :: last, ltmp(nshell), lnum, ll
!
      fock2(:)= zero
      fock3(:)= zero
!
      ncount=(2*nshell**3+3*nshell**2+nshell)/6+myrank
!
!$OMP parallel do schedule(dynamic,1) &
!$OMP private(jsh,ksh,lsh,ij,kl,ik,il,jk,jl,xijkl,denmax1,denmax,twoeri,ii,jj,kk, &
!$OMP icount,kstart,last,ltmp,lnum,ll) reduction(+:fock2,fock3)
      do ish= nshell,1,-1
        ii= ish*(ish-1)/2
        icount=ncount-(2*ish*ish*ish-3*ish*ish+ish)/6
        do jsh= 1,ish
          ij= ii+jsh
          jj= jsh*(jsh-1)/2
          kstart=mod(icount-ish*(jsh-1),nproc)+1
          do ksh= kstart,ish,nproc
            kk= ksh*(ksh-1)/2
            ik= ii+ksh
            jk= jj+ksh
            if(jsh.lt.ksh) jk= kk+jsh
            denmax1=max(four*dmax(ij),dmax(ik),dmax(jk))
            last= ksh
            if(ish == ksh) last= jsh
            ll=min(jsh,ksh)
            lnum=0
!           do lsh= 1,ksh
            do lsh= 1,ll
              kl= kk+lsh
              il= ii+lsh
              jl= jj+lsh
              xijkl=xint(ij)*xint(kl)
              denmax=max(denmax1,four*dmax(kl),dmax(il),dmax(jl))
              if(xijkl*denmax.ge.cutint2) then
                lnum=lnum+1
                ltmp(lnum)=lsh
              endif
            enddo
            do lsh= ll+1,last
              kl= kk+lsh
              il= ii+lsh
              jl= lsh*(lsh-1)/2+jsh
              xijkl=xint(ij)*xint(kl)
              denmax=max(denmax1,four*dmax(kl),dmax(il),dmax(jl))
              if(xijkl*denmax.ge.cutint2) then
                lnum=lnum+1
                ltmp(lnum)=lsh
              endif
            enddo
            do lsh= 1,lnum
              call calc2eri(twoeri,ish,jsh,ksh,ltmp(lsh),maxdim)
              call udftfockeri(fock2,fock3,dmtrxa,dmtrxb,twoeri,ish,jsh,ksh,ltmp(lsh),maxdim, &
&                              hfexchange)
            enddo
          enddo
        enddo
      enddo
!$OMP end parallel do
!
      call para_allreduce(fock2,fock1,nao*(nao+1)/2,MPI_SUM,MPI_COMM_WORLD)
      call para_allreduce(fock3,fock2,nao*(nao+1)/2,MPI_SUM,MPI_COMM_WORLD)
      return
end


!------------------------------------------------------------------------------------
  subroutine udftfockeri(focka,fockb,dmtrxa,dmtrxb,twoeri,ish,jsh,ksh,lsh,maxdim, &
&                        hfexchange)
!------------------------------------------------------------------------------------
!
! Form unrestricted DFT Fock matrix from two-electron intgrals
!
      use modbasis, only : nshell, nao, mbf, locbf
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: ish, jsh, ksh, lsh, maxdim
      integer :: nbfi, nbfj, nbfk, nbfl
      integer :: locbfi, locbfj, locbfk, locbfl, jmax, lmax, i, j, k, l, ij, kl
      integer :: nij, nkl, nik, nil, njk, njl
      integer :: iloc, jloc, kloc, lloc, iloc2, jloc2, kloc2, lloc2
      real(8),parameter :: half=0.5D+00, two=2.0D+00, four=4.0D+00
      real(8),intent(in) :: dmtrxa(nao*(nao+1)/2), dmtrxb(nao*(nao+1)/2)
      real(8),intent(in) :: twoeri(maxdim,maxdim,maxdim,maxdim), hfexchange
      real(8),intent(inout) :: focka(nao*(nao+1)/2), fockb(nao*(nao+1)/2)
      real(8) :: val, val2, val4
      logical :: ieqj, keql, ieqk, jeql, ikandjl, ijorkl
!
      nbfi = mbf(ish)
      nbfj = mbf(jsh)
      nbfk = mbf(ksh)
      nbfl = mbf(lsh)
      locbfi= locbf(ish)
      locbfj= locbf(jsh)
      locbfk= locbf(ksh)
      locbfl= locbf(lsh)
!
      ieqj= ish.eq.jsh
      keql= ksh.eq.lsh
      ieqk= ish.eq.ksh
      jeql= jsh.eq.lsh
      ikandjl= ieqk.and.jeql
      ijorkl= ieqj.or.keql
      jmax= nbfj
      lmax= nbfl
      ij= 0
!
      do i= 1,nbfi
        iloc= locbfi+i
        iloc2= iloc*(iloc-1)/2
        if(ieqj) jmax= i
        do j= 1,jmax
          jloc= locbfj+j
          jloc2= jloc*(jloc-1)/2
          ij= ij+1
          kl= 0
  kloop:  do k= 1,nbfk
            kloc= locbfk+k
            kloc2= kloc*(kloc-1)/2
            if(keql) lmax= k
            do l= 1,lmax
              kl= kl+1
              if(ikandjl.and.kl.gt.ij) exit kloop
              val= twoeri(l,k,j,i)
              if(abs(val).lt.cutint2) cycle
              lloc= locbfl+l
              nij= iloc2+jloc
              nkl= kloc2+lloc
              nik= iloc2+kloc
              nil= iloc2+lloc
              njk= jloc2+kloc
              njl= jloc2+lloc
              if(jloc.lt.kloc) njk= kloc2+jloc
              if(jloc.lt.lloc) njl= lloc*(lloc-1)/2+jloc
              if(ieqk) then
                if(iloc.lt.kloc) then
                  lloc2= lloc*(lloc-1)/2
                  nij= kloc2+lloc
                  nkl= iloc2+jloc
                  nik= kloc2+iloc
                  nil= kloc2+jloc
                  njk= lloc2+iloc
                  njl= lloc2+jloc
                  if(lloc.lt.iloc) njk= iloc2+lloc
                  if(lloc.lt.jloc) njl= jloc2+lloc
                elseif(iloc.eq.kloc.and.jloc.eq.lloc) then
                  val= val*half
                endif
              endif
              if(ijorkl) then
                if(iloc.eq.jloc) val= val*half
                if(kloc.eq.lloc) val= val*half
              endif
              val2= val*two*hfexchange
              val4= val*four
!
              focka(nij)= focka(nij)+val4*(dmtrxa(nkl)+dmtrxb(nkl))
              fockb(nij)= fockb(nij)+val4*(dmtrxa(nkl)+dmtrxb(nkl))
              focka(nkl)= focka(nkl)+val4*(dmtrxa(nij)+dmtrxb(nij))
              fockb(nkl)= fockb(nkl)+val4*(dmtrxa(nij)+dmtrxb(nij))
!
              focka(nik)= focka(nik)-val2*dmtrxa(njl)
              focka(nil)= focka(nil)-val2*dmtrxa(njk)
              focka(njk)= focka(njk)-val2*dmtrxa(nil)
              focka(njl)= focka(njl)-val2*dmtrxa(nik)
!
              fockb(nik)= fockb(nik)-val2*dmtrxb(njl)
              fockb(nil)= fockb(nil)-val2*dmtrxb(njk)
              fockb(njk)= fockb(njk)-val2*dmtrxb(nil)
              fockb(njl)= fockb(njl)-val2*dmtrxb(nik)
            enddo
          enddo kloop
        enddo
      enddo
      return
end






