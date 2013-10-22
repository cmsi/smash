!-----------------------------------------------------------------
  subroutine calcrhf(h1mtrx,cmoa,ortho,overlap,dmtrx,xint,eigen)
!-----------------------------------------------------------------
!
! Driver of restricted Hartree-Fock calculation
!
! In  : h1mtrx  (1-electron integral matrix)
!       cmoa    (MO coefficient matrix)
!       ortho   (orthogonalization matrix)
!       overlap (overlap integral matrix)
!       dmtrx   (density matrix)
!       xint    (exchange integral matrix)
!       eigen   (MO energy)
!
!     use procpar, only : master, nproc, myrank
      use procpar
      use basis, only : nshell, nao, mtype
      use molecule, only : neleca, nmo
      use scf, only : maxiter, fdiff, dconv, maxdiis, maxsoscf, dodiis
      use energy, only : enuc, escf, escfe
      use iofile, only : iout
      use thresh, only : threshsoscf, cutint2, threshex, threshover
      implicit none
      integer :: nao2, nao3, nshell3, maxdim, maxfunc(0:6), iter, i, itsub, itdiis
      integer :: itextra, itsoscf, nocc, nvir
!ishimura
      integer :: idis(nproc,10), isize1, isize2
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, two=2.0D+00
      real(8),intent(in) :: h1mtrx(nao*(nao+1)/2), ortho(nao*nao), overlap(nao*(nao+1)/2)
      real(8),intent(out) :: xint(nshell*(nshell+1)/2), eigen(nao)
      real(8),intent(inout) :: cmoa(nao*nao), dmtrx(nao*(nao+1)/2)
      real(8),allocatable :: fock(:), fockprev(:), dmtrxprev(:), dmax(:), work(:,:)
      real(8),allocatable :: fockdiis(:), errdiis(:), diismtrx(:), work2(:)
      real(8),allocatable :: hstart(:), sograd(:,:), sodisp(:), sovecy(:)
      real(8) :: escfprev, diffmax, tridot, deltae, errmax, sogradmax
      logical :: convsoscf
      data maxfunc/1,3,6,10,15,21,28/
!ishimura
    real(8) :: t1,t2,t3,t4
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
      call distarray(idis,nmo,nao,nao3,nocc,nvir,nproc)
!
! Set arrays
!
      call memset(nao2+nao3*3+nshell3)
      allocate(fock(nao3),fockprev(nao3),dmtrxprev(nao3),dmax(nshell3),work(nao,nao))
!
      isize1= max(idis(myrank+1,3),idis(myrank+1,7)*nao)
      isize2=idis(myrank+1,5)
      call memset(isize2*maxdiis+isize1*(maxdiis+1)+maxdiis*(maxdiis+1)/2)
      allocate(fockdiis(isize2*maxdiis), errdiis(isize1*maxdiis), &
&              diismtrx(maxdiis*(maxdiis+1)/2), work2(isize1))
      if(.not.dodiis) then
       call memset(nocc*nvir*3*maxsoscf)
       allocate(hstart(nocc*nvir),sograd(nocc*nvir,maxsoscf),sodisp(nocc*nvir*maxsoscf), &
&               sovecy(nocc*nvir*(maxsoscf-1)))
      endif
!
      escfprev= zero
      itdiis =1
      itextra=1
      itsoscf=0
      convsoscf=.false.
!
! Set 1-electron Hamiltonian
!
      call dcopy(nao3,h1mtrx,1,fockprev,1)
      call dcopy(nao3,dmtrx,1,dmtrxprev,1)
      call zeroclr(diismtrx,maxdiis*(maxdiis+1)/2)
!
! Calculate exchange integrals
!
      call calcexchange(xint,maxdim)
!
      if(master) then
        write(iout,'(1x,74(1H-))')
        write(iout,'("   Hartree-Fock calculation")')
        write(iout,'(1x,74(1H-))')
        write(iout,'("   DIIS = ",l1,",   SOSCF = ",l1)')dodiis,.not.dodiis
        write(iout,'("   Dconv      =",1p,d9.2,",  MaxIter    = ",i9  ,",  ThreshSOSCF=",d9.2)') &
&                    dconv, maxiter, threshsoscf
        write(iout,'("   Cutint2    =",1p,d9.2,",  ThreshEx   = ",d9.2,",  ThreshOver =",d9.2)') &
&                    cutint2, threshex, threshover
        write(iout,'(1x,74(1H-))')
        write(iout,'(" ====================")')
        write(iout,'("    SCF Iteration")')
        write(iout,'(" ====================")')
        if(dodiis) then
          write(iout,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                      "Delta Density     DIIS Error")')
        else
          write(iout,'(" Iter SubIt   Total Energy      Delta Energy      ", &
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
        call calcdmax(dmtrx,dmax)
!
! Calculate two-electron integrals and Fock matrix
!
        call cpu_time(t1)
        call zeroclr(fock,nao3)
        call formfock(fock,dmtrx,dmax,xint,maxdim)
        call dscal(nao3,half,fock,1)
        do i= 1,nao
          fock(i*(i+1)/2)= fock(i*(i+1)/2)*two
        enddo
        call cpu_time(t2)
!
! Form full Fock matrix
!
        call daxpy(nao3,one,fockprev,1,fock,1)
        call dcopy(nao3,fock,1,fockprev,1)
!
! Calculate SCF ENERGY
!
        escfe=(tridot(dmtrxprev,fock,nao)+tridot(dmtrxprev,h1mtrx,nao))*half
        escf = escfe+enuc
        deltae = escf-escfprev
        escfprev= escf
!
        if(dodiis) then
!
! DIIS interpolation
!
          call calcdiis(fock,dmtrxprev,overlap,ortho,errdiis,fockdiis,diismtrx, &
&                       cmoa,work,work2,errmax,idis,itdiis)
!
! Extrapolate Fock matrix
!
          if(itdiis == 1) call extrap(fock,fockdiis,work,cmoa,dmtrx,idis,itextra)
!
! Diagonalize Fock matrix
!
          call diagfock(fock,work,ortho,cmoa,work2,eigen,idis)
        else
!
! Approximated Second-order SCF method
!
          if((itsoscf == 0).or.(convsoscf)) then
            call diagfock(fock,work,ortho,cmoa,work2,eigen,idis)
            sogradmax= zero
            itsoscf= itsoscf+1
          else
            call expand(fock,work,nao)
            call soscfgrad(work,work2,sograd(1,itsoscf),cmoa,nocc,nvir,sogradmax,idis)
            if(sogradmax <= threshsoscf) then
              if(itsoscf == 1) call soscfhess(hstart,eigen,nocc,nvir)
              call soscfupdate(cmoa,hstart,sograd,sodisp,sovecy,work,work2, &
&                              nocc,nvir,itsoscf,maxsoscf,idis)
              itsoscf= itsoscf+1
            else
              call diagfock(fock,work,ortho,cmoa,work2,eigen,idis)
            endif
          endif
        endif
        call cpu_time(t3)
!
! Copy previous density matrix and calculate new density matrix
!
        call calcdmtrx(cmoa,dmtrx,work,nao,neleca)
        call ddiff(dmtrx,dmtrxprev,work,nao3,diffmax)
        if(dodiis) then
          if(itdiis > 1) then
            itsub= itdiis-1
          else
            itsub= itextra-1
          endif
        else
          itsub= itsoscf
        endif
        if(master) then
          if(dodiis) then
            write(iout,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,errmax
          else
            write(iout,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,sogradmax
          endif
        endif
!
! Check SCF convergence
!
        if(dodiis) then
          if(diffmax.lt.dconv) exit
          if(itdiis >= maxdiis) itdiis= 1
        else
          if((diffmax.lt.dconv).and.(convsoscf)) exit
          if((diffmax.lt.dconv).and.(.not.convsoscf)) convsoscf=.true.
          if(itsoscf >= maxsoscf) itsoscf= 0
        endif
        if(iter.eq.maxiter) then
          if(master) then
            write(iout,'(" Not Converged.")')
            call iabort
          endif
        endif
        call dcopy(nao3,dmtrx,1,dmtrxprev,1)
        call dcopy(nao3,work,1,dmtrx,1)
        call cpu_time(t4)
        if(master) write(*,'(10x,6f8.3)')t2-t1,t3-t2,t4-t3
      enddo
!
      if(master) then
        write(iout,'(" -----------------------------------------")')
        write(iout,'("    SCF Converged.")')
        write(iout,'("    RHF Energy = ",f17.9," a.u.")')escf
        write(iout,'(" -----------------------------------------"/)')
      endif
!
! Unset arrays
!
      if(.not.dodiis) then
       call memunset(nocc*nvir*3*maxsoscf)
       deallocate(hstart,sograd,sodisp,sovecy)
      endif
      deallocate(fockdiis,errdiis,diismtrx,work2)
      call memunset(isize2*maxdiis+isize1*(maxdiis+1)+maxdiis*(maxdiis+1)/2)
      deallocate(fock,fockprev,dmtrxprev,dmax,work)
      call memunset(nao2+nao3*3+nshell3)
      return
end


!-------------------------------------------------------------
  subroutine diagfock(fock,work,ortho,cmoa,work2,eigen,idis)
!-------------------------------------------------------------
!
! Driver of Fock matrix diagonalization
!
      use procpar
      use basis, only : nao
      use molecule, only : nmo
      implicit none
      integer,intent(in) :: idis(nproc,10)
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: fock(nao*(nao+1)/2), ortho(nao,nao)
      real(8),intent(out) :: work(nao,nao), cmoa(nao*nao), work2(*), eigen(nao)
!
! Canonicalize Fock matrix
!
        call expand(fock,work,nao)
        call canonicalizep(work,ortho,cmoa,work2,nao,nmo,idis)
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
      call para_allgatherv(work2,idis(myrank+1,3),'D',cmoa,idis(1,3),idis(1,4),MPI_COMM_WORLD)
!
      return
end


!---------------------------------------------------
  subroutine formfock(fock,dmtrx,dmax,xint,maxdim)
!---------------------------------------------------
!
! Driver of Fock matrix formation from two-electron intgrals
!
      use procpar, only : nproc, myrank, MPI_SUM, MPI_COMM_WORLD
      use basis, only : nshell, nao
      use thresh, only : cutint2
      implicit none
      integer,intent(in) :: maxdim
      integer :: ish, jsh, ksh, lsh, ij, kl, ik, il, jk, jl
      integer :: ii, jj, kk, kstart
      integer(8) :: ncount, icount
      real(8),parameter :: four=4.0D+00
      real(8),intent(in) :: dmtrx(nao*(nao+1)/2), dmax(nshell*(nshell+1)/2)
      real(8),intent(in) :: xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: fock(nao*(nao+1)/2)
      real(8) :: xijkl, denmax, twoeri(maxdim**4), denmax1
      integer :: last, ltmp(nshell),lnum,ll
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
      call para_allreduce(fock,xijkl,nao*(nao+1)/2,'D',MPI_SUM,MPI_COMM_WORLD,1)
      return
end


!---------------------------------------------------------------
  subroutine fockeri(fock,dmtrx,twoeri,ish,jsh,ksh,lsh,maxdim)
!---------------------------------------------------------------
!
! Form Fock matrix from two-electron intgrals
!
      use basis, only : nshell, nao, mbf, locbf
      use thresh, only : cutint2
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


!-----------------------------------------------------------------
  subroutine formdftfock(fock,dmtrx,dmax,xint,maxdim,hfexchange)
!-----------------------------------------------------------------
!
! Driver of DFT Fock matrix formation from two-electron intgrals
!
      use procpar, only : nproc, myrank, MPI_SUM, MPI_COMM_WORLD
      use basis, only : nshell, nao
      use thresh, only : cutint2
      implicit none
      integer,intent(in) :: maxdim
      integer :: ish, jsh, ksh, lsh, ij, kl, ik, il, jk, jl
      integer :: ii, jj, kk, kstart
      integer(8) :: ncount, icount
      real(8),parameter :: four=4.0D+00
      real(8),intent(in) :: dmtrx(nao*(nao+1)/2), dmax(nshell*(nshell+1)/2)
      real(8),intent(in) :: xint(nshell*(nshell+1)/2), hfexchange
      real(8),intent(out) :: fock(nao*(nao+1)/2)
      real(8) :: xijkl, denmax, twoeri(maxdim**4),denmax1
!
      ncount=(2*nshell**3+3*nshell**2+nshell)/6+myrank
!
!$OMP parallel do schedule(dynamic,1) &
!$OMP private(jsh,ksh,lsh,ij,kl,ik,il,jk,jl,xijkl,denmax1,denmax,twoeri,ii,jj,kk,icount,kstart) &
!$OMP reduction(+:fock)
      do ish= nshell,1,-1
        ii= ish*(ish-1)/2
        icount=ncount-(2*ish*ish*ish-3*ish*ish+ish)/6
        do jsh= 1,ish
          ij= ii+jsh
          jj= jsh*(jsh-1)/2
          kstart=mod(icount-ish*(jsh-1),nproc)+1
 kloop:   do ksh= kstart,ish,nproc
            kk= ksh*(ksh-1)/2
            ik= ii+ksh
            jk= jj+ksh
            if(jsh.lt.ksh) jk= kk+jsh
            denmax1=max(four*dmax(ij),dmax(ik),dmax(jk))
            do lsh= 1,ksh
              kl= kk+lsh
              if(kl.gt.ij) exit kloop
              xijkl=xint(ij)*xint(kl)
              il= ii+lsh
              jl= jj+lsh
              if(jsh.lt.lsh) jl= lsh*(lsh-1)/2+jsh
              denmax=max(denmax1,four*dmax(kl),dmax(il),dmax(jl))
              if(xijkl*denmax.lt.cutint2) cycle
              call calc2eri(twoeri,ish,jsh,ksh,lsh,maxdim)
              call fockdfteri(fock,dmtrx,twoeri,ish,jsh,ksh,lsh,maxdim,hfexchange)
            enddo
          enddo kloop
        enddo
      enddo
!$OMP end parallel do
!
      call para_allreduce(fock,xijkl,nao*(nao+1)/2,'D',MPI_SUM,MPI_COMM_WORLD,1)
      return
end


!-----------------------------------------------------------------------------
  subroutine fockdfteri(fock,dmtrx,twoeri,ish,jsh,ksh,lsh,maxdim,hfexchange)
!-----------------------------------------------------------------------------
!
! Form DFT Fock matrix from two-electron intgrals
!
      use basis, only : nshell, nao, mbf, locbf
      use thresh, only : cutint2
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
!ishimura
              if(ikandjl.and.kl.gt.ij) exit kloop
!             if(ikandjl.and.kl.gt.ij) cycle
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



!------------------------------------------------------------------
  subroutine calcrdft(h1mtrx,cmoa,ortho,overlap,dmtrx,xint,eigen)
!------------------------------------------------------------------
!
! Driver of restricted DFT calculation
!
! In  : h1mtrx  (1-electron integral matrix)
!       cmoa    (MO coefficient matrix)
!       ortho   (orthogonalization matrix)
!       overlap (overlap integral matrix)
!       dmtrx   (density matrix)
!       xint    (exchange integral matrix)
!       eigen   (MO energy)
!
!     use procpar, only : master, nproc, myrank
      use procpar
      use dft, only : nrad, nleb, hfexchange
      use atominfo, only : atomrad
      use basis, only : nshell, nao, mtype
      use molecule, only : neleca, nmo, natom, numatomic
      use scf, only : maxiter, fdiff, dconv, maxdiis, maxsoscf, dodiis
      use energy, only : enuc, escf, escfe
      use iofile, only : iout
      use thresh, only : threshsoscf, cutint2, threshex, threshover, threshrho, threshdfock
      implicit none
      integer :: nao2, nao3, nshell3, maxdim, maxfunc(0:6), iter, i, itsub, itdiis
      integer :: itextra, itsoscf, nocc, nvir
      integer :: idis(nproc,10), isize1, isize2, iatom
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, two=2.0D+00
      real(8),intent(in) :: h1mtrx(nao*(nao+1)/2), ortho(nao*nao), overlap(nao*(nao+1)/2)
      real(8),intent(out) :: xint(nshell*(nshell+1)/2), eigen(nao)
      real(8),intent(inout) :: cmoa(nao*nao), dmtrx(nao*(nao+1)/2)
      real(8),allocatable :: fock(:), fockprev(:), dmtrxprev(:), dmax(:), work(:,:)
      real(8),allocatable :: fockdiis(:), errdiis(:), diismtrx(:), work2(:)
      real(8),allocatable :: hstart(:), sograd(:,:), sodisp(:), sovecy(:)
      real(8),allocatable :: fockd(:), rad(:), atomvec(:), surface(:)
      real(8),allocatable :: radpt(:), angpt(:), ptweight(:)
      real(8) :: escfprev, diffmax, tridot, deltae, errmax, sogradmax, edft, totalelec
      logical :: convsoscf
      data maxfunc/1,3,6,10,15,21,28/
!ishimura
    real(8) :: t1,t2,t3,t4
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
      call distarray(idis,nmo,nao,nao3,nocc,nvir,nproc)
!
! Set arrays
!
      call memset(nao2+nao3*3+nshell3)
      allocate(fock(nao3),fockprev(nao3),dmtrxprev(nao3),dmax(nshell3),work(nao,nao))
      call memset(nao3+natom+6*natom*natom+2*nrad+4*nleb+natom*nrad*nleb)
      allocate(fockd(nao3),rad(natom),atomvec(5*natom*natom), &
&              surface(natom*natom),radpt(2*nrad),angpt(4*nleb),ptweight(natom*nrad*nleb))
!
      isize1= max(idis(myrank+1,3),idis(myrank+1,7)*nao)
      isize2=idis(myrank+1,5)
      call memset(isize2*maxdiis+isize1*(maxdiis+1)+maxdiis*(maxdiis+1)/2)
      allocate(fockdiis(isize2*maxdiis), errdiis(isize1*maxdiis), &
&              diismtrx(maxdiis*(maxdiis+1)/2), work2(isize1))
      if(.not.dodiis) then
       call memset(nocc*nvir*3*maxsoscf)
       allocate(hstart(nocc*nvir),sograd(nocc*nvir,maxsoscf),sodisp(nocc*nvir*maxsoscf), &
&               sovecy(nocc*nvir*(maxsoscf-1)))
      endif
!
      escfprev= zero
      itdiis =1
      itextra=1
      itsoscf=0
      convsoscf=.false.
!
! Calculate DFT information
!
      hfexchange=0.2D+00
      call calcatomvec(atomvec,surface)
      call calcradpt(radpt,nrad)
      call calclebpt(angpt,nleb)
      do iatom= 1,natom
        rad(iatom)= atomrad(numatomic(iatom))
      enddo
      call calcgridweight(ptweight,rad,radpt,angpt,atomvec,surface)
!
! Set 1-electron Hamiltonian
!
      call dcopy(nao3,h1mtrx,1,fockprev,1)
      call dcopy(nao3,dmtrx,1,dmtrxprev,1)
      call zeroclr(diismtrx,maxdiis*(maxdiis+1)/2)
!
! Calculate exchange integrals
!
      call calcexchange(xint,maxdim)
!
      if(master) then
        write(iout,'(1x,74(1H-))')
        write(iout,'("   DFT calculation")')
        write(iout,'(1x,74(1H-))')
        write(iout,'("   DIIS = ",l1,",   SOSCF = ",l1)')dodiis,.not.dodiis
        write(iout,'("   Dconv      =",1p,d9.2,",  MaxIter    = ",i9  ,",  ThreshSOSCF=",d9.2)') &
&                    dconv, maxiter, threshsoscf
        write(iout,'("   Cutint2    =",1p,d9.2,",  ThreshEx   = ",d9.2,",  ThreshOver =",d9.2)') &
&                    cutint2, threshex, threshover
        write(iout,'("   Nrad       =",1p,i9  ,",  Nleb       = ",i9  ,",  ThreshRho  =",d9.2)') &
&                    nrad, nleb, threshrho
        write(iout,'("   ThreshDfock=",1p,d9.2)') threshdfock
        write(iout,'(1x,74(1H-))')
        write(iout,'(" ====================")')
        write(iout,'("    SCF Iteration")')
        write(iout,'(" ====================")')
        if(dodiis) then
          write(iout,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                      "Delta Density     DIIS Error")')
        else
          write(iout,'(" Iter SubIt   Total Energy      Delta Energy      ", &
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
        call calcdmax(dmtrx,dmax)
!
! Calculate two-electron integrals and Fock matrix
!
        call cpu_time(t1)
        call zeroclr(fock,nao3)
        call zeroclr(fockd,nao3)
        call formdftfock(fock,dmtrx,dmax,xint,maxdim,hfexchange)
!       call formfockexcor(fockd,edft,cmoa,atomvec,radpt,angpt,rad,ptweight)
        call formfockexcor(fockd,edft,totalelec,cmoa,atomvec,radpt,angpt,rad,ptweight)
        call dscal(nao3,half,fock,1)
        do i= 1,nao
          fock(i*(i+1)/2)= fock(i*(i+1)/2)*two
        enddo
        call cpu_time(t2)
!
! Form Fock matrix
!
        call daxpy(nao3,one,fockprev,1,fock,1)
        call dcopy(nao3,fock,1,fockprev,1)
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
        if(dodiis) then
!
! DIIS interpolation
!
          call calcdiis(fock,dmtrxprev,overlap,ortho,errdiis,fockdiis,diismtrx, &
&                       cmoa,work,work2,errmax,idis,itdiis)
!
! Extrapolate Fock matrix
!
          if(itdiis == 1) call extrap(fock,fockdiis,work,cmoa,dmtrx,idis,itextra)
!
! Diagonalize Fock matrix
!
          call diagfock(fock,work,ortho,cmoa,work2,eigen,idis)
        else
!
! Approximated Second-order SCF method
!
          if((itsoscf == 0).or.(convsoscf)) then
            call diagfock(fock,work,ortho,cmoa,work2,eigen,idis)
            sogradmax= zero
            itsoscf= itsoscf+1
          else
            call expand(fock,work,nao)
            call soscfgrad(work,work2,sograd(1,itsoscf),cmoa,nocc,nvir,sogradmax,idis)
            if(sogradmax <= threshsoscf) then
              if(itsoscf == 1) call soscfhess(hstart,eigen,nocc,nvir)
              call soscfupdate(cmoa,hstart,sograd,sodisp,sovecy,work,work2, &
&                              nocc,nvir,itsoscf,maxsoscf,idis)
              itsoscf= itsoscf+1
            else
              call diagfock(fock,work,ortho,cmoa,work2,eigen,idis)
            endif
          endif
        endif
        call cpu_time(t3)
!
! Copy previous density matrix and calculate new density matrix
!
        call calcdmtrx(cmoa,dmtrx,work,nao,neleca)
        call ddiff(dmtrx,dmtrxprev,work,nao3,diffmax)
        if(dodiis) then
          if(itdiis > 1) then
            itsub= itdiis-1
          else
            itsub= itextra-1
          endif
        else
          itsub= itsoscf
        endif
        if(master) then
          if(dodiis) then
            write(iout,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,errmax
          else
            write(iout,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,sogradmax
          endif
        endif
!
! Check SCF convergence
!
        if(dodiis) then
          if(diffmax.lt.dconv) exit
          if(itdiis >= maxdiis) itdiis= 1
        else
          if((diffmax.lt.dconv).and.(convsoscf)) exit
          if((diffmax.lt.dconv).and.(.not.convsoscf)) convsoscf=.true.
          if(itsoscf >= maxsoscf) itsoscf= 0
        endif
        if(iter.eq.maxiter) then
          if(master) then
            write(iout,'(" Not Converged.")')
            call iabort
          endif
        endif
        call dcopy(nao3,dmtrx,1,dmtrxprev,1)
        call dcopy(nao3,work,1,dmtrx,1)
        call cpu_time(t4)
        if(master) write(*,'(10x,6f8.3)')t2-t1,t3-t2,t4-t3
      enddo
!
      if(master) then
        write(iout,'(" -----------------------------------------------------------")')
        write(iout,'("    SCF Converged.")')
        write(iout,'("    DFT Energy = ",f17.9," a.u.")')escf
        write(iout,'("    Exchange + Correlation energy = ",f17.9," a.u.")')edft
        write(iout,'("    Number of electrons           = ",f17.9)')totalelec
        write(iout,'(" -----------------------------------------------------------"/)')
      endif
!
! Unset arrays
!
      if(.not.dodiis) then
       call memunset(nocc*nvir*3*maxsoscf)
       deallocate(hstart,sograd,sodisp,sovecy)
      endif
      deallocate(fockd,rad,atomvec,surface,radpt,angpt,ptweight)
      call memunset(nao3+natom+6*natom*natom+2*nrad+4*nleb+natom*nrad*nleb)
      deallocate(fockdiis,errdiis,diismtrx,work2)
      call memunset(isize2*maxdiis+isize1*(maxdiis+1)+maxdiis*(maxdiis+1)/2)
      deallocate(fock,fockprev,dmtrxprev,dmax,work)
      call memunset(nao2+nao3*3+nshell3)
      return
end













