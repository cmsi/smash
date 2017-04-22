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
!------------------------------------------------------------------------
  subroutine calcrhf(h1mtrx,cmo,ortho,overlap,dmtrx,xint,eigen, &
&                    nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!------------------------------------------------------------------------
!
! Driver of restricted Hartree-Fock calculation
!
! In  : h1mtrx  (1-electron integral matrix)
!       ortho   (Orthogonalization matrix)
!       overlap (Overlap integral matrix)
! Out : dmtrx   (Density matrix)
!       xint    ((ij|ij) integral matrix)
!       eigen   (MO energy)
! Inout : cmo   (MO coefficient matrix)
!
      use modparallel, only : master
      use modbasis, only : nshell, nao, mtype
      use modmolecule, only : neleca, nmo
      use modscf, only : maxiter, fdiff, dconv, maxdiis, maxsoscf, maxqc, maxqcdiag, &
&                        maxqcdiagsub, scfconv, extrap
      use modenergy, only : enuc, escf, escfe
      use modthresh, only : threshsoscf, threshqc, cutint2, threshex, threshover, threshdiis
      use modprint, only : iprint
      use modwarn, only : nwarn
      implicit none
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2, mpi_comm1, mpi_comm2
      integer :: nao2, nao3, nshell3, maxdim, maxfunc(0:6), iter, itsub, itdiis
      integer :: itextra, itsoscf, itqc, nocc, nvir
      integer :: idis(nproc2,14), isize1, isize2, isize3
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00
      real(8),parameter :: small=1.0D-10
      real(8),intent(in) :: h1mtrx(nao*(nao+1)/2), ortho(nao*nao), overlap(nao*(nao+1)/2)
      real(8),intent(out) :: dmtrx(nao*(nao+1)/2), xint(nshell*(nshell+1)/2), eigen(nao)
      real(8),intent(inout) :: cmo(nao*nao)
      real(8),allocatable :: fock(:), fockprev(:), dmtrxprev(:), dmax(:), work(:)
      real(8),allocatable :: fockdiis(:), errdiis(:), diismtrx(:), work2(:)
      real(8),allocatable :: hstart(:), sograd(:,:), sodisp(:), sovecy(:)
      real(8),allocatable :: qcvec(:), qcmat(:), qcmatsave(:), qceigen(:), qcgmn(:)
      real(8) :: escfprev, diffmax, tridot, deltae, errmax, sogradmax, sodispmax
      real(8) :: time1, time2, time3, time4
      logical :: convsoscf, convqc
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
      call distarray(idis,nmo,nao,nao3,nocc,nvir,nocc,nvir,nproc2)
!
! Set arrays
!
      isize1= max(idis(myrank2+1,3),idis(myrank2+1,7)*nao,maxdiis)
      isize2= idis(myrank2+1,3)
      isize3= idis(myrank2+1,5)
      select case(scfconv)
        case('DIIS')
          call memset(nao2+nao3*3+nshell3+isize3*maxdiis+isize1+isize2*maxdiis &
&                    +maxdiis*(maxdiis+1)/2)
          allocate(fock(nao3),fockprev(nao3),dmtrxprev(nao3),dmax(nshell3),work(nao2), &
&                  fockdiis(isize3*maxdiis), errdiis(isize2*maxdiis), &
&                  diismtrx(maxdiis*(maxdiis+1)/2),work2(isize1))
        case('SOSCF')
          call memset(nao2+nao3*3+nshell3+nocc*nvir*3*maxsoscf+isize1)
          allocate(fock(nao3),fockprev(nao3),dmtrxprev(nao3),dmax(nshell3),work(nao2), &
                   hstart(nocc*nvir),sograd(nocc*nvir,maxsoscf),sodisp(nocc*nvir*maxsoscf), &
&                  sovecy(nocc*nvir*(maxsoscf-1)),work2(isize1))
        case('QC')
          call memset(nao2*2+nao3*2+nshell3+(nocc*nvir+1)*(maxqcdiagsub+1)*2+maxqcdiagsub**2 &
&                    +maxqcdiagsub*(maxqcdiagsub+1)/2+maxqcdiagsub)
          allocate(fock(nao2),fockprev(nao3),dmtrxprev(nao3),dmax(nshell3),work(nao2), &
                   qcvec((nocc*nvir+1)*(maxqcdiagsub+1)*2),qcmat(maxqcdiagsub**2), &
&                  qcmatsave(maxqcdiagsub*(maxqcdiagsub+1)/2),qceigen(maxqcdiagsub), &
&                  qcgmn(nao3),work2(nao2))
        case default
          if(master) then
            write(*,'(" SCFConv=",a12,"is not supported.")') 
            call iabort
          endif
      end select
!
      escfprev= zero
      itdiis =0
      itextra=0
      itsoscf=0
      itqc   =0
      convsoscf=.false.
      convqc=.false.
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
      call calcschwarzeri(xint,work,maxdim,nproc2,myrank2,mpi_comm2)
!
      if(master) then
        write(*,'(1x,74("-"))')
        write(*,'("   Restricted Hartree-Fock calculation")')
        write(*,'(1x,74("-"))')
        write(*,'("   SCFConv    = ",a8,",  Dconv      =",1p,d9.2,",  MaxIter    =",i9)') &
&                    scfconv, dconv, maxiter
        write(*,'("   Cutint2    =",1p,d9.2,",  ThreshEx   =",d9.2,",  ThreshOver =",d9.2)') &
&                    cutint2, threshex, threshover
        select case(scfconv)
          case('DIIS')
            write(*,'("   MaxDIIS    =",i9,",  ThreshDIIS =",1p,d9.2)') &
&                    maxdiis, threshdiis
          case('SOSCF')
            write(*,'("   MaxSOSCF   =",i9,",  ThreshSOSCF=",1p,d9.2)') &
&                    maxdiis, threshdiis
          case('QC')
            write(*,'("   MaxQC      =",i9,",  ThreshQC   =",1p,d9.2,",  MaxQCDiag  =",i9)') &
&                    maxqc, threshqc, maxqcdiag
            write(*,'("   MaxQCDiagSub=",i8)') &
&                    maxqcdiagsub
        end select
        write(*,'(1x,74("-"))')
!
        if((nao /= nmo).and.master) then
          write(*,'(/," Warning! Number of MOs is reduced from",i5," to",i5,/)') nao, nmo
          nwarn= nwarn+1
        endif
!
        write(*,'(" ====================")')
        write(*,'("    SCF Iteration")')
        write(*,'(" ====================")')
        select case(scfconv)
          case('DIIS')
            write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                     "Delta Density     DIIS Error")')
          case('SOSCF')
            write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                     "Delta Density    Orbital Grad")')
          case('QC')
            write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                     "Delta Density")')
        end select
      endif
!
! Start SCF iteration
!
      do iter= 1,maxiter
!
! Calculate maximum density matrix elements
!
        call calcrdmax(dmtrx,dmax,work,nproc2,myrank2,mpi_comm2)
!
! Calculate two-electron integrals and Fock matrix
!
        call cpu_time(time1)
        call formrfock(fock,work,dmtrx,dmax,xint,maxdim,nproc1,myrank1,mpi_comm1)
        call dscal(nao3,half,fock,1)
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
        select case(scfconv)
!
! DIIS interpolation
!
          case('DIIS')
            call calcdiiserr(fock,dmtrxprev,overlap,ortho,cmo,work,work2,errmax,nao,nmo, &
&                            idis,nproc2,myrank2,mpi_comm2)
            if(((itdiis /= 0).or.(errmax <= threshdiis)).and.(errmax > small))then
              itdiis= itdiis+1
              call calcrdiis(fock,errdiis,fockdiis,diismtrx,cmo,work2,itdiis,nao,maxdiis, &
&                            idis,nproc2,myrank2,mpi_comm2)
            endif
!
! Extrapolate Fock matrix
!
            if(extrap.and.itdiis == 0) &
&             call fockextrap(fock,fockdiis,work,cmo,dmtrx,itextra,nao,maxdiis, &
&                             idis,nproc2,myrank2,mpi_comm2)
!
! Diagonalize Fock matrix
!
            call diagfock(fock,work,ortho,cmo,work2,eigen,idis,nproc2,myrank2,mpi_comm2)
!
! Approximated Second-order SCF method
!
          case('SOSCF')
            if((itsoscf == 0).or.(convsoscf)) then
              call diagfock(fock,work,ortho,cmo,work2,eigen,idis,nproc2,myrank2,mpi_comm2)
              sogradmax= zero
              itsoscf= itsoscf+1
            else
              call expand(fock,work,nao)
              call soscfgrad(work,work2,sograd(1,itsoscf),cmo,nocc,nvir,sogradmax,nao, &
&                            idis,nproc2,myrank2,mpi_comm2,1)
              if(sogradmax <= threshsoscf) then
                if(itsoscf == 1) call soscfinith(hstart,eigen,nocc,nvir,nao)
                call soscfnewh(hstart,sograd,sodisp,sovecy,nocc,nvir,itsoscf,maxsoscf,sodispmax)
                call soscfupdate(cmo,sodisp,work,work2,nocc,nvir,itsoscf,maxsoscf, &
&                                nao,nmo,sodispmax,idis,nproc2,myrank2,mpi_comm2)
                itsoscf= itsoscf+1
              else
                call diagfock(fock,work,ortho,cmo,work2,eigen,idis,nproc2,myrank2,mpi_comm2)
              endif
            endif
!
! Quadratically convergent SCF method
!
          case('QC')
            if((itqc == 0).or.(convqc)) then
              call diagfock(fock,work,ortho,cmo,work2,eigen,idis,nproc2,myrank2,mpi_comm2)
              itqc= itqc+1
            elseif(itqc == maxqc) then
              call diagfock(fock,work,ortho,cmo,work2,eigen,idis,nproc2,myrank2,mpi_comm2)
              itqc= 1
            else
              call rhfqc(fock,cmo,dmax,qcgmn,qcvec,qcmat,qcmatsave,qceigen,overlap,xint, &
&                        work,work2,one,nao,nmo,nocc,nvir,nshell,maxdim,maxqcdiag, &
&                        maxqcdiagsub,threshqc,nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
              itqc= itqc+1
            endif
        end select
        call cpu_time(time3)
!
! Copy previous density matrix and calculate new density matrix
!
        call calcdmtrx(cmo,dmtrx,work,nao,neleca)
        call ddiff(dmtrx,dmtrxprev,work,nao3,diffmax)
!
        select case(scfconv)
          case('DIIS')
            if(extrap.and.(itdiis==0)) then
              itsub= itextra
            else
              itsub= itdiis
            endif
            if(master) &
&             write(*,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,errmax
          case('SOSCF')
            itsub= itsoscf
            if(master) &
&             write(*,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,sogradmax
          case('QC')
            itsub= itqc
            if(master) &
&             write(*,'(1x,i3,2x,i3,2(1x,f17.9),1f17.9)')iter,itsub,escf,deltae,diffmax
        end select
!
! Check SCF convergence
!
        select case(scfconv)
          case('DIIS')
            if(diffmax < dconv) exit
            if(itdiis >= maxdiis) itdiis= 0
          case('SOSCF')
            if((diffmax < dconv).and.(convsoscf)) exit
            if((diffmax < dconv).and.(itsoscf == 1)) exit
            if((diffmax < dconv).and.(.not.convsoscf)) convsoscf=.true.
            if(itsoscf >= maxsoscf) itsoscf= 0
          case('QC')
            if((diffmax < dconv).and.(convqc)) exit
            if((diffmax < dconv).and.(itqc == 1)) exit
            if((diffmax < dconv).and.(.not.convqc)) convqc=.true.
            if((diffmax >= dconv).and.(convqc)) then
              itqc= 0
              convqc=.false.
            endif
        end select
!
        if(iter == maxiter) then
          if(master) then
            write(*,'(" SCF did not converge.")')
            call iabort
          endif
        endif
        call dcopy(nao3,dmtrx,1,dmtrxprev,1)
        call dcopy(nao3,work,1,dmtrx,1)
        call cpu_time(time4)
        if(master.and.(iprint >= 3)) write(*,'(10x,6f8.3)')time2-time1,time3-time2,time4-time3
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
      select case(scfconv)
        case('DIIS')
          call memunset(nao2+nao3*3+nshell3+isize3*maxdiis+isize1+isize2*maxdiis+ &
&                       maxdiis*(maxdiis+1)/2)
          deallocate(fock,fockprev,dmtrxprev,dmax,work, &
&                    fockdiis,errdiis, &
&                    diismtrx,work2)
        case('SOSCF')
          call memunset(nao2+nao3*3+nshell3+nocc*nvir*3*maxsoscf+isize1)
          deallocate(fock,fockprev,dmtrxprev,dmax,work, &
                     hstart,sograd,sodisp, &
&                    sovecy,work2)
        case('QC')
          call memunset(nao2*2+nao3*2+nshell3+(nocc*nvir+1)*(maxqcdiagsub+1)*2+maxqcdiagsub**2 &
&                      +maxqcdiagsub*(maxqcdiagsub+1)/2+maxqcdiagsub)
          deallocate(fock,fockprev,dmtrxprev,dmax,work, &
&                    qcvec,qcmat, &
&                    qcmatsave,qceigen, &
&                    qcgmn,work2)
      end select
!
      return
end


!----------------------------------------------------------------------------------
  subroutine diagfock(fock,work,ortho,cmo,work2,eigen,idis,nproc,myrank,mpi_comm)
!----------------------------------------------------------------------------------
!
! Driver of Fock matrix diagonalization
!
! In  : fock  (Fock matrix)
!       ortho (Orthogonalization matrix)
! Out : cmo   (Canonical MO matrx)
!       work, work2, eigen (work space)
!
      use modbasis, only : nao
      use modmolecule, only : nmo
      implicit none
      integer,intent(in) :: nproc, myrank, mpi_comm, idis(nproc,14)
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: fock(nao*(nao+1)/2), ortho(nao,nao)
      real(8),intent(out) :: work(nao,nao), cmo(nao*nao), work2(*), eigen(nao)
!
! Canonicalize Fock matrix
!
        call expand(fock,work,nao)
        call canonicalizep(work,ortho,cmo,work2,nao,nmo,idis,nproc,myrank,mpi_comm)
!
! Diagonalize canonicalized matrix
!
        call diag('V','U',nmo,work,nao,eigen,nproc,myrank,mpi_comm)
!
! Backtransform to AO basis
!
      if(idis(myrank+1,1) /= 0) &
&       call dgemm('N','N',nao,idis(myrank+1,1),nmo,one,ortho,nao, &
&                  work(1,idis(myrank+1,2)+1),nao,zero,work2,nao)
      call para_allgathervr(work2,idis(myrank+1,3),cmo,idis(1,3),idis(1,4),nproc,mpi_comm)
!
      return
end


!------------------------------------------------------------------------------------
  subroutine formrfock(focktotal,fock,dmtrx,dmax,xint,maxdim,nproc,myrank,mpi_comm)
!------------------------------------------------------------------------------------
!
! Driver of Fock matrix formation from two-electron intgrals
!
      use modbasis, only : nshell, nao
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: maxdim, nproc, myrank, mpi_comm
      integer :: ijsh, ish, jsh, ksh, lsh, ij, kl, ik, il, jk, jl
      integer :: ii, jj, kk, kstart, ishcheck
      integer(8) :: ncount, icount(nshell)
      real(8),parameter :: zero=0.0D+00, two=2.0D+00, four=4.0D+00
      real(8),intent(in) :: dmtrx(nao*(nao+1)/2), dmax(nshell*(nshell+1)/2)
      real(8),intent(in) :: xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: focktotal(nao*(nao+1)/2), fock(nao*(nao+1)/2)
      real(8) :: xijkl, denmax, twoeri(maxdim**4), denmax1
      integer :: last, ltmp(nshell), lnum, ll
!
      fock(:)= zero
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
!$OMP private(ijsh,jsh,ksh,lsh,ij,kl,ik,il,jk,jl,xijkl,denmax1,denmax,twoeri,jj,kk, &
!$OMP kstart,last,ltmp,lnum,ll) firstprivate(ish,ii) reduction(+:fock)
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
        jj= jsh*(jsh-1)/2
        kstart=mod(icount(ish)-ish*(jsh-1),nproc)+1
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
!         do lsh= 1,ksh
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
            call rfockeri(fock,dmtrx,twoeri,ish,jsh,ksh,ltmp(lsh),maxdim)
          enddo
        enddo
      enddo
!$OMP end parallel do
!
      do ii= 1,nao
        fock(ii*(ii+1)/2)= fock(ii*(ii+1)/2)*two
      enddo
!
      call para_allreducer(fock,focktotal,nao*(nao+1)/2,mpi_comm)
      return
end


!----------------------------------------------------------------
  subroutine rfockeri(fock,dmtrx,twoeri,ish,jsh,ksh,lsh,maxdim)
!----------------------------------------------------------------
!
! Form Fock matrix from two-electron intgrals
!
      use modbasis, only : nao, mbf, locbf
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: ish, jsh, ksh, lsh, maxdim
      integer :: nbfi, nbfj, nbfk, nbfl
      integer :: locbfi, locbfj, locbfk, locbfl, jmax, lmax, i, j, k, l, ij, kl
      integer :: nij, nkl, nik, nil, njk, njl
      integer :: iloc, jloc, kloc, lloc, iloc2, jloc2, kloc2, lloc2, kloc0, jloc0
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
      jloc0= locbfj*(locbfj-1)/2
      kloc0= locbfk*(locbfk-1)/2
!
      if(.not.ieqk)then
        do i= 1,nbfi
          iloc= locbfi+i
          iloc2= iloc*(iloc-1)/2
          jloc2= jloc0
          if(ieqj) jmax= i
          do j= 1,jmax
            jloc= locbfj+j
            jloc2= jloc2+jloc-1
            kloc2= kloc0
            nij= iloc2+jloc
            do k= 1,nbfk
              kloc= locbfk+k
              kloc2= kloc2+kloc-1
              nik= iloc2+kloc
              njk= jloc2+kloc
              if(jloc.lt.kloc) njk= kloc2+jloc
              if(keql) lmax= k
              do l= 1,lmax
                val= twoeri(l,k,j,i)
                if(abs(val).lt.cutint2) cycle
                lloc= locbfl+l
                nkl= kloc2+lloc
                nil= iloc2+lloc
                njl= jloc2+lloc
                if(jloc.lt.lloc) njl= lloc*(lloc-1)/2+jloc
                if(ijorkl) then
                  if(iloc.eq.jloc) val= val*half
                  if(kloc.eq.lloc) val= val*half
                endif
                val4= val*four
                fock(nij)= fock(nij)+val4*dmtrx(nkl)
                fock(nkl)= fock(nkl)+val4*dmtrx(nij)
                fock(nik)= fock(nik)-val *dmtrx(njl)
                fock(nil)= fock(nil)-val *dmtrx(njk)
                fock(njk)= fock(njk)-val *dmtrx(nil)
                fock(njl)= fock(njl)-val *dmtrx(nik)
              enddo
            enddo
          enddo
        enddo
      else
        do i= 1,nbfi
          iloc= locbfi+i
          iloc2= iloc*(iloc-1)/2
          jloc2= jloc0
          if(ieqj) jmax= i
          do j= 1,jmax
            jloc= locbfj+j
            jloc2= jloc2+jloc-1
            kloc2= kloc0
            ij= ij+1
            kl= 0
     kloop: do k= 1,nbfk
              kloc= locbfk+k
              kloc2= kloc2+kloc-1
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
      endif
!
      return
end


!------------------------------------------------------------------------------
  subroutine formrdftfock(focktotal,fock,dmtrx,dmax,xint,maxdim,hfexchange, &
&                         nproc,myrank,mpi_comm)
!------------------------------------------------------------------------------
!
! Driver of DFT Fock matrix formation from two-electron intgrals
!
      use modbasis, only : nshell, nao
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: maxdim, nproc, myrank, mpi_comm
      integer :: ijsh, ish, jsh, ksh, lsh, ij, kl, ik, il, jk, jl
      integer :: ii, jj, kk, kstart, ishcheck
      integer(8) :: ncount, icount(nshell)
      real(8),parameter :: zero=0.0D+00, two=2.0D+00, four=4.0D+00
      real(8),intent(in) :: dmtrx(nao*(nao+1)/2), dmax(nshell*(nshell+1)/2)
      real(8),intent(in) :: xint(nshell*(nshell+1)/2), hfexchange
      real(8),intent(out) :: focktotal(nao*(nao+1)/2), fock(nao*(nao+1)/2)
      real(8) :: xijkl, denmax, twoeri(maxdim**4), denmax1
      integer :: last, ltmp(nshell), lnum, ll
!
      fock(:)= zero
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
!$OMP private(ijsh,jsh,ksh,lsh,ij,kl,ik,il,jk,jl,xijkl,denmax1,denmax,twoeri,jj,kk, &
!$OMP kstart,last,ltmp,lnum,ll) firstprivate(ish,ii) reduction(+:fock)
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
        jj= jsh*(jsh-1)/2
        kstart=mod(icount(ish)-ish*(jsh-1),nproc)+1
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
!         do lsh= 1,ksh
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
            call rdftfockeri(fock,dmtrx,twoeri,ish,jsh,ksh,ltmp(lsh),maxdim,hfexchange)
          enddo
        enddo
      enddo
!$OMP end parallel do
!
      do ii= 1,nao
        fock(ii*(ii+1)/2)= fock(ii*(ii+1)/2)*two
      enddo
!
      call para_allreducer(fock,focktotal,nao*(nao+1)/2,mpi_comm)
      return
end


!------------------------------------------------------------------------------
  subroutine rdftfockeri(fock,dmtrx,twoeri,ish,jsh,ksh,lsh,maxdim,hfexchange)
!------------------------------------------------------------------------------
!
! Form DFT Fock matrix from two-electron intgrals
!
      use modbasis, only : nao, mbf, locbf
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: ish, jsh, ksh, lsh, maxdim
      integer :: nbfi, nbfj, nbfk, nbfl
      integer :: locbfi, locbfj, locbfk, locbfl, jmax, lmax, i, j, k, l, ij, kl
      integer :: nij, nkl, nik, nil, njk, njl
      integer :: iloc, jloc, kloc, lloc, iloc2, jloc2, kloc2, lloc2, jloc0, kloc0
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
      jloc0= locbfj*(locbfj-1)/2
      kloc0= locbfk*(locbfk-1)/2
!
      if(.not.ieqk)then
        do i= 1,nbfi
          iloc= locbfi+i
          iloc2= iloc*(iloc-1)/2
          jloc2= jloc0
          if(ieqj) jmax= i
          do j= 1,jmax
            jloc= locbfj+j
            jloc2= jloc2+jloc-1
            kloc2= kloc0
            nij= iloc2+jloc
            do k= 1,nbfk
              kloc= locbfk+k
              kloc2= kloc2+kloc-1
              nik= iloc2+kloc
              njk= jloc2+kloc
              if(jloc.lt.kloc) njk= kloc2+jloc
              if(keql) lmax= k
              do l= 1,lmax
                val= twoeri(l,k,j,i)
                if(abs(val).lt.cutint2) cycle
                lloc= locbfl+l
                nkl= kloc2+lloc
                nil= iloc2+lloc
                njl= jloc2+lloc
                if(jloc.lt.lloc) njl= lloc*(lloc-1)/2+jloc
                if(ijorkl) then
                  if(iloc.eq.jloc) val= val*half
                  if(kloc.eq.lloc) val= val*half
                endif
                val4= val*four
                val = val*hfexchange
                fock(nij)= fock(nij)+val4*dmtrx(nkl)
                fock(nkl)= fock(nkl)+val4*dmtrx(nij)
                fock(nik)= fock(nik)-val *dmtrx(njl)
                fock(nil)= fock(nil)-val *dmtrx(njk)
                fock(njk)= fock(njk)-val *dmtrx(nil)
                fock(njl)= fock(njl)-val *dmtrx(nik)
              enddo
            enddo
          enddo
        enddo
      else
        do i= 1,nbfi
          iloc= locbfi+i
          iloc2= iloc*(iloc-1)/2
          jloc2= jloc0
          if(ieqj) jmax= i
          do j= 1,jmax
            jloc= locbfj+j
            jloc2= jloc2+jloc-1
            kloc2= kloc0
            ij= ij+1
            kl= 0
     kloop: do k= 1,nbfk
              kloc= locbfk+k
              kloc2= kloc2+kloc-1
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
                if(ijorkl) then
                  if(iloc.eq.jloc) val= val*half
                  if(kloc.eq.lloc) val= val*half
                endif
                val4= val*four
                val = val*hfexchange
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
      endif
!
      return
end


!-------------------------------------------------------------------------
  subroutine calcrdft(h1mtrx,cmo,ortho,overlap,dmtrx,xint,eigen, &
&                     nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!-------------------------------------------------------------------------
!
! Driver of restricted DFT calculation
!
! In  : h1mtrx  (1-electron integral matrix)
!       ortho   (Orthogonalization matrix)
!       overlap (Overlap integral matrix)
! Out : dmtrx   (Density matrix)
!       xint    ((ij|ij) integral matrix)
!       eigen   (MO energy)
! Inout : cmo   (MO coefficient matrix)
!
      use modparallel, only : master
      use moddft, only : idftex, idftcor, nrad, nleb, hfexchange
      use modatom, only : atomrad
      use modbasis, only : nshell, nao, mtype
      use modmolecule, only : neleca, nmo, natom, numatomic
      use modscf, only : maxiter, fdiff, dconv, maxdiis, maxsoscf, maxqc, maxqcdiag, &
                         maxqcdiagsub, scfconv, extrap
      use modenergy, only : enuc, escf, escfe
      use modthresh, only : threshsoscf, threshqc, cutint2, threshex, threshover, threshdiis, &
&                           threshweight, threshrho, threshdfock, threshdftao
      use modprint, only : iprint
      use modunit, only : tobohr
      use modwarn, only : nwarn
      implicit none
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2, mpi_comm1, mpi_comm2
      integer :: nao2, nao3, nshell3, maxdim, maxfunc(0:6), iter, itsub, itdiis
      integer :: itextra, itsoscf, itqc, nocc, nvir
      integer :: idis(nproc2,14), isize1, isize2, isize3, iatom
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00
      real(8),parameter :: small=1.0D-10
      real(8),intent(in) :: h1mtrx(nao*(nao+1)/2), ortho(nao*nao), overlap(nao*(nao+1)/2)
      real(8),intent(out) :: dmtrx(nao*(nao+1)/2), xint(nshell*(nshell+1)/2), eigen(nao)
      real(8),intent(inout) :: cmo(nao*nao)
      real(8),allocatable :: fock(:), fockprev(:), dmtrxprev(:), dmax(:), work(:)
      real(8),allocatable :: fockdiis(:), errdiis(:), diismtrx(:), work2(:)
      real(8),allocatable :: hstart(:), sograd(:,:), sodisp(:), sovecy(:)
      real(8),allocatable :: fockd(:), rad(:), atomvec(:), surface(:),  radpt(:), angpt(:)
      real(8),allocatable :: ptweight(:), xyzpt(:), rsqrd(:), vao(:), vmo(:)
      real(8),allocatable :: qcvec(:), qcmat(:), qcmatsave(:), qceigen(:), qcgmn(:)
      real(8) :: escfprev, diffmax, tridot, deltae, errmax, sogradmax, sodispmax
      real(8) :: edft, totalelec
      real(8) :: time1, time2, time3, time4
      logical :: convsoscf, convqc
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
      call distarray(idis,nmo,nao,nao3,nocc,nvir,nocc,nvir,nproc2)
!
! Set arrays
!
      isize1= max(idis(myrank2+1,3),idis(myrank2+1,7)*nao,3*natom,nao,maxdiis)
      isize2= idis(myrank2+1,3)
      isize3=idis(myrank2+1,5)
      select case(scfconv)
        case('DIIS')
          call memset(nao2+nao3*4+nshell3+natom*5+natom*natom*6+nrad*2+nleb*4+natom*nrad*nleb &
&                    +nao*4+nocc*4+isize3*maxdiis+isize1+isize2*maxdiis+maxdiis*(maxdiis+1)/2)
          allocate(fock(nao3),fockprev(nao3),dmtrxprev(nao3),dmax(nshell3),work(nao2), &
&                  fockd(nao3),rad(natom),atomvec(5*natom*natom),surface(natom*natom), &
&                  radpt(2*nrad),angpt(4*nleb),ptweight(natom*nrad*nleb),xyzpt(3*natom), &
&                  rsqrd(natom),vao(4*nao),vmo(4*nocc),fockdiis(isize3*maxdiis), &
&                  errdiis(isize2*maxdiis),diismtrx(maxdiis*(maxdiis+1)/2),work2(isize1))
        case('SOSCF')
          call memset(nao2+nao3*4+nshell3+natom*5+natom*natom*6+nrad*2+nleb*4+natom*nrad*nleb &
&                    +nao*4+nocc*4+isize1+nocc*nvir*3*maxsoscf)
          allocate(fock(nao3),fockprev(nao3),dmtrxprev(nao3),dmax(nshell3),work(nao2), &
&                  fockd(nao3),rad(natom),atomvec(5*natom*natom),surface(natom*natom), &
&                  radpt(2*nrad),angpt(4*nleb),ptweight(natom*nrad*nleb),xyzpt(3*natom), &
&                  rsqrd(natom),vao(4*nao),vmo(4*nocc),work2(isize1), &
&                  hstart(nocc*nvir),sograd(nocc*nvir,maxsoscf),sodisp(nocc*nvir*maxsoscf), &
&                  sovecy(nocc*nvir*(maxsoscf-1)))
        case ('QC')
          isize1= max(isize1,nao2)
          call memset(nao2*2+nao3*4+nshell3+natom*5+natom*natom*6+nrad*2+nleb*4+natom*nrad*nleb &
&                    +nao*4+nocc*4+isize1+(nocc*nvir+1)*(maxqcdiagsub+1)*2 &
&                    +maxqcdiagsub*(maxqcdiagsub*3+1)/2+maxqcdiagsub)
          allocate(fock(nao2),fockprev(nao3),dmtrxprev(nao3),dmax(nshell3),work(nao2), &
&                  fockd(nao3),rad(natom),atomvec(5*natom*natom),surface(natom*natom), &
&                  radpt(2*nrad),angpt(4*nleb),ptweight(natom*nrad*nleb),xyzpt(3*natom), &
&                  rsqrd(natom),vao(4*nao),vmo(4*nocc),work2(isize1), &
&                  qcvec((nocc*nvir+1)*(maxqcdiagsub+1)*2),qcmat(maxqcdiagsub*maxqcdiagsub), &
&                  qcmatsave(maxqcdiagsub*(maxqcdiagsub+1)/2),qceigen(maxqcdiagsub),qcgmn(nao3))
        case default
          if(master) then
            write(*,'(" SCFConv=",a12,"is not supported.")')
            call iabort
          endif
      end select
!
      escfprev= zero
      itdiis =0
      itextra=0
      itsoscf=0
      itqc   =0
      convsoscf=.false.
      convqc=.false.
!
! Calculate DFT information
!
      call calcatomvec(atomvec,surface)
      call calcradpt(radpt,nrad)
      call calclebpt(angpt,nleb)
      do iatom= 1,natom
        rad(iatom)= atomrad(numatomic(iatom))*tobohr
      enddo
      call calcgridweight(ptweight,rad,radpt,angpt,atomvec,surface,xyzpt,work2,nproc1,myrank1)
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
      call calcschwarzeri(xint,work,maxdim,nproc2,myrank2,mpi_comm2)
!
      if(master) then
        write(*,'(1x,74("-"))')
        write(*,'("   Restricted DFT calculation")')
        write(*,'(1x,74("-"))')
        write(*,'("   SCFConv    = ",a8,",  Dconv      =",1p,d9.2,",  MaxIter    =",i9)') &
&                    scfconv, dconv, maxiter
        write(*,'("   Cutint2    =",1p,d9.2,",  ThreshEx   =",d9.2,",  ThreshOver =",d9.2)') &
&                    cutint2, threshex, threshover
        write(*,'("   Nrad       =",i9  ,",  Nleb       =",i9  ,",  ThreshRho  =",1p,d9.2)') &
&                    nrad, nleb, threshrho
        write(*,'("   ThreshDfock=",1p,d9.2,",  Threshdftao=",d9.2,",  ThreshWeight=",d8.2)') &
&                    threshdfock, threshdftao, threshweight
        select case(scfconv)
          case('DIIS')
            write(*,'("   MaxDIIS    =",i9,",  ThreshDIIS =",1p,d9.2)') &
&                    maxdiis, threshdiis
          case('SOSCF')
            write(*,'("   MaxSOSCF   =",i9,",  ThreshSOSCF=",1p,d9.2)') &
&                    maxdiis, threshdiis
          case('QC')
            write(*,'("   MaxQC      =",i9,",  ThreshQC   =",1p,d9.2,",  MaxQCDiag  =",i9)') &
&                    maxqc, threshqc, maxqcdiag
            write(*,'("   MaxQCDiagSub=",i8)') &
&                    maxqcdiagsub
        end select
        write(*,'(1x,74("-"))')
!
        if((nao /= nmo).and.master) then
          write(*,'(/," Warning! Number of MOs is reduced from",i5," to",i5,/)') nao, nmo
          nwarn= nwarn+1
        endif
!
        write(*,'(" ====================")')
        write(*,'("    SCF Iteration")')
        write(*,'(" ====================")')
        select case(scfconv)
          case('DIIS')
            write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                     "Delta Density     DIIS Error")')
          case('SOSCF')
            write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                     "Delta Density    Orbital Grad")')
          case('QC')
            write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                     "Delta Density")')
        end select
      endif
!
! Start SCF iteration
!
      do iter= 1,maxiter
!
! Calculate maximum density matrix elements
!
        call calcrdmax(dmtrx,dmax,work,nproc2,myrank2,mpi_comm2)
        call cpu_time(time1)
!
! Calculate exchange-correlation terms
!
        call formrfockexcor(fockd,fock,edft,totalelec,cmo,atomvec,radpt,angpt, &
&                           rad,ptweight,vao,vmo,xyzpt,rsqrd,work,work2,idftex,idftcor, &
&                           nproc1,myrank1,mpi_comm1)
!
! Calculate two-electron integrals
!
        call formrdftfock(fock,work,dmtrx,dmax,xint,maxdim,hfexchange, &
&                         nproc1,myrank1,mpi_comm1)
        call dscal(nao3,half,fock,1)
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
        select case(scfconv)
!
! DIIS interpolation
!
          case('DIIS')
            call calcdiiserr(fock,dmtrxprev,overlap,ortho,cmo,work,work2,errmax,nao,nmo, &
&                            idis,nproc2,myrank2,mpi_comm2)
            if(((itdiis /= 0).or.(errmax <= threshdiis)).and.(errmax > small))then
              itdiis= itdiis+1
              call calcrdiis(fock,errdiis,fockdiis,diismtrx,cmo,work2,itdiis,nao,maxdiis, &
&                            idis,nproc2,myrank2,mpi_comm2)
            endif
!
! Extrapolate Fock matrix
!
            if(extrap.and.itdiis == 0) &
&             call fockextrap(fock,fockdiis,work,cmo,dmtrx,itextra,nao,maxdiis, &
&                             idis,nproc2,myrank2,mpi_comm2)
!
! Diagonalize Fock matrix
!
            call diagfock(fock,work,ortho,cmo,work2,eigen,idis,nproc2,myrank2,mpi_comm2)
!
! Approximated Second-order SCF method
!
          case('SOSCF')
            if((itsoscf == 0).or.(convsoscf)) then
              call diagfock(fock,work,ortho,cmo,work2,eigen,idis,nproc2,myrank2,mpi_comm2)
              sogradmax= zero
              itsoscf= itsoscf+1
            else
              call expand(fock,work,nao)
              call soscfgrad(work,work2,sograd(1,itsoscf),cmo,nocc,nvir,sogradmax,nao, &
&                            idis,nproc2,myrank2,mpi_comm2,1)
              if(sogradmax <= threshsoscf) then
                if(itsoscf == 1) call soscfinith(hstart,eigen,nocc,nvir,nao)
                call soscfnewh(hstart,sograd,sodisp,sovecy,nocc,nvir,itsoscf,maxsoscf,sodispmax)
                call soscfupdate(cmo,sodisp,work,work2,nocc,nvir,itsoscf,maxsoscf, &
&                                nao,nmo,sodispmax,idis,nproc2,myrank2,mpi_comm2)
                itsoscf= itsoscf+1
              else
                call diagfock(fock,work,ortho,cmo,work2,eigen,idis,nproc2,myrank2,mpi_comm2)
              endif
            endif
!
! Quadratically convergent SCF method
!
          case('QC')
            if((itqc == 0).or.(convqc)) then
              call diagfock(fock,work,ortho,cmo,work2,eigen,idis,nproc2,myrank2,mpi_comm2)
              itqc= itqc+1
            elseif((itqc == maxqc)) then
              call diagfock(fock,work,ortho,cmo,work2,eigen,idis,nproc2,myrank2,mpi_comm2)
              itqc= 1
            else
              call rhfqc(fock,cmo,dmax,qcgmn,qcvec,qcmat,qcmatsave,qceigen,overlap,xint, &
&                        work,work2,hfexchange,nao,nmo,nocc,nvir,nshell,maxdim,maxqcdiag, &
&                        maxqcdiagsub,threshqc,nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
              itqc= itqc+1
            endif
        end select
        call cpu_time(time3)
!
! Copy previous density matrix and calculate new density matrix
!
        call calcdmtrx(cmo,dmtrx,work,nao,neleca)
        call ddiff(dmtrx,dmtrxprev,work,nao3,diffmax)
!
        select case(scfconv)
          case('DIIS')
            if(extrap.and.(itdiis==0)) then
              itsub= itextra
            else
              itsub= itdiis
            endif
            if(master) &
&             write(*,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,errmax
          case('SOSCF')
            itsub= itsoscf
            if(master) &
&             write(*,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,sogradmax
          case('QC')
            itsub= itqc
            if(master) &
&             write(*,'(1x,i3,2x,i3,2(1x,f17.9),1f17.9)')iter,itsub,escf,deltae,diffmax
        end select
!
! Check SCF convergence
!
        select case(scfconv)
          case('DIIS')
            if(diffmax < dconv) exit
            if(itdiis >= maxdiis) itdiis= 0
          case('SOSCF')
            if((diffmax < dconv).and.(convsoscf)) exit
            if((diffmax < dconv).and.(itsoscf == 1)) exit
            if((diffmax < dconv).and.(.not.convsoscf)) convsoscf=.true.
            if(itsoscf >= maxsoscf) itsoscf= 0
          case('QC')
            if((diffmax < dconv).and.(convqc)) exit
            if((diffmax < dconv).and.(itqc == 1)) exit
            if((diffmax < dconv).and.(.not.convqc)) convqc=.true.
            if((diffmax >= dconv).and.(convqc)) then
              itqc= 0
              convqc=.false.
            endif
        end select
!
        if(iter == maxiter) then
          if(master) then
            write(*,'(" SCF did not converge.")')
            call iabort
          endif
        endif
        call dcopy(nao3,dmtrx,1,dmtrxprev,1)
        call dcopy(nao3,work,1,dmtrx,1)
        call cpu_time(time4)
        if(master.and.(iprint >= 3)) write(*,'(10x,6f8.3)')time2-time1,time3-time2,time4-time3
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
      select case(scfconv)
        case('DIIS')
          deallocate(fock,fockprev,dmtrxprev,dmax,work, &
&                    fockd,rad,atomvec,surface, &
&                    radpt,angpt,ptweight,xyzpt, &
&                    rsqrd,vao,vmo,fockdiis, &
&                    errdiis,diismtrx,work2)
          call memunset(nao2+nao3*4+nshell3+natom*5+natom*natom*6+nrad*2+nleb*4+natom*nrad*nleb &
&                      +nao*4+nocc*4+isize3*maxdiis+isize1+isize2*maxdiis+maxdiis*(maxdiis+1)/2)
        case('SOSCF')
          deallocate(fock,fockprev,dmtrxprev,dmax,work, &
&                    fockd,rad,atomvec,surface, &
&                    radpt,angpt,ptweight,xyzpt, &
&                    rsqrd,vao,vmo,work2, &
&                    hstart,sograd,sodisp, &
&                    sovecy)
          call memunset(nao2+nao3*4+nshell3+natom*5+natom*natom*6+nrad*2+nleb*4+natom*nrad*nleb &
&                      +nao*4+nocc*4+isize1+nocc*nvir*3*maxsoscf)
        case('QC')
          deallocate(fock,fockprev,dmtrxprev,dmax,work, &
&                    fockd,rad,atomvec,surface, &
&                    radpt,angpt,ptweight,xyzpt, &
&                    rsqrd,vao,vmo,work2, &
&                    qcvec,qcmat, &
&                    qcmatsave,qceigen,qcgmn)
          call memunset(nao2*2+nao3*4+nshell3+natom*5+natom*natom*6+nrad*2+nleb*4+natom*nrad*nleb &
&                      +nao*4+nocc*4+isize1+(nocc*nvir+1)*(maxqcdiagsub+1)*2 &
&                      +maxqcdiagsub*(maxqcdiagsub*3+1)/2+maxqcdiagsub)
      end select
!
      return
end


!----------------------------------------------------------------------------------------
  subroutine calcuhf(h1mtrx,cmoa,cmob,ortho,overlap,dmtrxa,dmtrxb,xint,eigena,eigenb, &
&                    nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!----------------------------------------------------------------------------------------
!
! Driver of unrestricted Hartree-Fock calculation
!
! In  : h1mtrx  (1-electron integral matrix)
!       ortho   (Orthogonalization matrix)
!       overlap (Overlap integral matrix)
! Out : dmtrxa  (Alpha density matrix)
!       dmtrxb  (Beta density matrix)
!       xint    ((ij|ij) integral matrix)
!       eigena  (Alpha MO energy)
!       eigenb  (Beta MO energy)
! Inout : cmoa  (Alpha MO coefficient matrix)
!         cmob  (Beta MO coefficient matrix)
!
      use modparallel, only : master
      use modbasis, only : nshell, nao, mtype
      use modmolecule, only : neleca, nelecb, nmo
      use modscf, only : maxiter, fdiff, dconv, maxdiis, maxsoscf, maxqc, maxqcdiag, &
&                        maxqcdiagsub, scfconv, extrap
      use modenergy, only : enuc, escf, escfe
      use modthresh, only : threshsoscf, threshqc, cutint2, threshex, threshover, threshdiis
      use modprint, only : iprint
      use modwarn, only : nwarn
      implicit none
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2, mpi_comm1, mpi_comm2
      integer :: nao2, nao3, nshell3, maxdim, maxfunc(0:6), iter, itsub, itdiis
      integer :: itextra, itsoscf, itqc, nocca, nvira, noccb, nvirb
      integer :: idis(nproc2,14), isize1, isize2, isize3
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00
      real(8),parameter :: small=1.0D-10
      real(8),intent(in) :: h1mtrx(nao*(nao+1)/2), ortho(nao*nao), overlap(nao*(nao+1)/2)
      real(8),intent(out) :: dmtrxa(nao*(nao+1)/2), dmtrxb(nao*(nao+1)/2)
      real(8),intent(out) :: xint(nshell*(nshell+1)/2), eigena(nao), eigenb(nao)
      real(8),intent(inout) :: cmoa(nao*nao), cmob(nao*nao)
      real(8),allocatable :: focka(:), fockb(:), fockpreva(:), fockprevb(:) 
      real(8),allocatable :: dmtrxpreva(:), dmtrxprevb(:), dmax(:), work(:)
      real(8),allocatable :: fockdiisa(:), fockdiisb(:), errdiisa(:), errdiisb(:)
      real(8),allocatable :: diismtrx(:), work2(:), work3(:)
      real(8),allocatable :: hstarta(:), hstartb(:), sograda(:,:), sogradb(:,:)
      real(8),allocatable :: sodispa(:), sodispb(:), sovecya(:), sovecyb(:)
      real(8),allocatable :: qcvec(:), qcmat(:), qcmatsave(:), qceigen(:)
      real(8),allocatable :: qcgmna(:), qcgmnb(:)
      real(8) :: escfprev, diffmax, diffmaxa, diffmaxb, tridot, deltae
      real(8) :: errmax, errmaxa, errmaxb, sogradmax, sogradmaxa, sogradmaxb, sodispmax
      real(8) :: s2, sz
      real(8) :: time1, time2, time3, time4
      logical :: convsoscf, convqc
      data maxfunc/1,3,6,10,15,21,28/
!
      nao2= nao*nao
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
      call distarray(idis,nmo,nao,nao3,nocca,nvira,noccb,nvirb,nproc2)
!
! Set arrays
!
      isize1= max(nao2,maxdiis)
      isize2=idis(myrank2+1,3)
      isize3=idis(myrank2+1,5)
      select case(scfconv)
        case('DIIS')
          call memset(nao3*8+nshell3+isize3*maxdiis*2+isize1+isize2*maxdiis*2 &
&                    +maxdiis*(maxdiis+1)/2+isize2)
          allocate(focka(nao3),fockb(nao3),fockpreva(nao3),fockprevb(nao3), &
&                  dmtrxpreva(nao3),dmtrxprevb(nao3),dmax(nshell3), &
&                  fockdiisa(isize3*maxdiis),fockdiisb(isize3*maxdiis),errdiisa(isize2*maxdiis),&
&                  errdiisb(isize2*maxdiis),diismtrx(maxdiis*(maxdiis+1)/2), &
&                  work(nao3*2),work2(isize1),work3(isize2))
        case('SOSCF')
          call memset(nao3*8+nshell3+nocca*nvira*3*maxsoscf+noccb*nvirb*3*maxsoscf &
&                    +isize1+isize2)
          allocate(focka(nao3),fockb(nao3),fockpreva(nao3),fockprevb(nao3), &
&                  dmtrxpreva(nao3),dmtrxprevb(nao3),dmax(nshell3), &
&                  hstarta(nocca*nvira),hstartb(noccb*nvirb), &
&                  sograda(nocca*nvira,maxsoscf),sogradb(noccb*nvirb,maxsoscf), &
&                  sodispa(nocca*nvira*maxsoscf),sodispb(noccb*nvirb*maxsoscf), &
&                  sovecya(nocca*nvira*(maxsoscf-1)),sovecyb(noccb*nvirb*(maxsoscf-1)), &
&                  work(nao3*2),work2(isize1),work3(isize2))
        case('QC')
          call memset(nao2*4+nao3*8+(nocca*nvira+noccb*nvirb+1)*(maxqcdiagsub+1)*2 &
&                    +maxqcdiagsub**2+maxqcdiagsub*(maxqcdiagsub+1)/2+maxqcdiagsub)
          allocate(focka(nao2),fockb(nao2),fockpreva(nao3),fockprevb(nao3), &
&                  dmtrxpreva(nao3),dmtrxprevb(nao3),dmax(nshell3), &
&                  qcvec((nocca*nvira+noccb*nvirb+1)*(maxqcdiagsub+1)*2), &
&                  qcmat(maxqcdiagsub**2),qcmatsave(maxqcdiagsub*(maxqcdiagsub+1)/2), &
&                  qceigen(maxqcdiagsub),qcgmna(nao3),qcgmnb(nao3), &
&                  work(nao3*2),work2(nao2),work3(nao2))
        case default
          if(master) then
            write(*,'(" SCFConv=",a12,"is not supported.")')
            call iabort
          endif
      end select
!
      escfprev= zero
      itdiis =0
      itextra=0
      itsoscf=0
      itqc   =0
      convsoscf=.false.
      convqc=.false.
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
      call calcschwarzeri(xint,work,maxdim,nproc2,myrank2,mpi_comm2)
!
      if(master) then
        write(*,'(1x,74("-"))')
        write(*,'("   Unrestricted Hartree-Fock calculation")')
        write(*,'(1x,74("-"))')
        write(*,'("   SCFConv    = ",a8,",  Dconv      =",1p,d9.2,",  MaxIter    =",i9)') &
&                    scfconv, dconv, maxiter
        write(*,'("   Cutint2    =",1p,d9.2,",  ThreshEx   =",d9.2,",  ThreshOver =",d9.2)') &
&                    cutint2, threshex, threshover
        select case(scfconv)
          case('DIIS')
            write(*,'("   MaxDIIS    =",i9,",  ThreshDIIS =",1p,d9.2)') &
&                    maxdiis, threshdiis
          case('SOSCF')
            write(*,'("   MaxSOSCF   =",i9,",  ThreshSOSCF=",1p,d9.2)') &
&                    maxdiis, threshdiis
          case('QC')
            write(*,'("   MaxQC      =",i9,",  ThreshQC   =",1p,d9.2,",  MaxQCDiag  =",i9)') &
&                    maxqc, threshqc, maxqcdiag
            write(*,'("   MaxQCDiagSub=",i8)') &
&                    maxqcdiagsub
        end select
        write(*,'(1x,74("-"))')
!
        if((nao /= nmo).and.master) then
          write(*,'(/," Warning! Number of MOs is reduced from",i5," to",i5,/)') nao, nmo
          nwarn= nwarn+1
        endif
!
        write(*,'(" ====================")')
        write(*,'("    SCF Iteration")')
        write(*,'(" ====================")')
        select case(scfconv)
          case('DIIS')
            write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                     "Delta Density     DIIS Error")')
          case('SOSCF')
            write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                     "Delta Density    Orbital Grad")')
          case('QC')
            write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                     "Delta Density")')
        end select
      endif
!
! Start SCF iteration
!
      do iter= 1,maxiter
!
! Calculate maximum density matrix elements
!
        call calcudmax(dmtrxa,dmtrxb,dmax,work,nproc2,myrank2,mpi_comm2)
!
! Calculate two-electron integrals and Fock matrix
!
        call cpu_time(time1)
        call formufock(focka,fockb,work,dmtrxa,dmtrxb,dmax,xint,maxdim,nproc1,myrank1,mpi_comm1)
        call dscal(nao3,half,focka,1)
        call dscal(nao3,half,fockb,1)
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
        select case(scfconv)
!
! DIIS interpolation
!
          case('DIIS')
            call calcdiiserr(focka,dmtrxpreva,overlap,ortho,cmoa,work,work2,errmaxa,nao,nmo, &
&                            idis,nproc2,myrank2,mpi_comm2)
            call calcdiiserr(fockb,dmtrxprevb,overlap,ortho,cmob,work,work2,errmaxb,nao,nmo, &
&                            idis,nproc2,myrank2,mpi_comm2)
            errmax= max(errmaxa,errmaxb)
            if(((itdiis /= 0).or.(errmax <= threshdiis)).and.(errmax > small))then
              itdiis= itdiis+1
              call calcudiis(focka,fockb,errdiisa,errdiisb,fockdiisa,fockdiisb, &
&                            diismtrx,cmoa,cmob,work2,itdiis,nao,maxdiis, &
&                            idis,nproc2,myrank2,mpi_comm2)
            endif
!
! Extrapolate Fock matrix
!
            if(extrap.and.itdiis == 0) then
              call fockextrap(focka,fockdiisa,work,cmoa,dmtrxa,itextra,nao,maxdiis, &
&                             idis,nproc2,myrank2,mpi_comm2)
              call fockextrap(fockb,fockdiisb,work,cmob,dmtrxb,itextra,nao,maxdiis, &
&                             idis,nproc2,myrank2,mpi_comm2)
            endif
!
! Diagonalize Fock matrix
!
            call diagfock(focka,work,ortho,cmoa,work2,eigena,idis,nproc2,myrank2,mpi_comm2)
            call diagfock(fockb,work,ortho,cmob,work2,eigenb,idis,nproc2,myrank2,mpi_comm2)
!
! Approximated Second-order SCF method
!
          case('SOSCF')
            if((itsoscf == 0).or.(convsoscf)) then
              call diagfock(focka,work,ortho,cmoa,work2,eigena,idis,nproc2,myrank2,mpi_comm2)
              call diagfock(fockb,work,ortho,cmob,work2,eigenb,idis,nproc2,myrank2,mpi_comm2)
              sogradmax= zero
              itsoscf= itsoscf+1
            else
              call expand(focka,work,nao)
              call soscfgrad(work,work2,sograda(1,itsoscf),cmoa,nocca,nvira,sogradmaxa,nao, &
&                            idis,nproc2,myrank2,mpi_comm2,1)
              if(noccb /= 0) then
                call expand(fockb,work,nao)
                call soscfgrad(work,work2,sogradb(1,itsoscf),cmob,noccb,nvirb,sogradmaxb,nao, &
&                              idis,nproc2,myrank2,mpi_comm2,2)
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
&                               nocca,noccb,nvira,nvirb,itsoscf,maxsoscf,sodispmax)
                call soscfupdate(cmoa,sodispa,work,work2,nocca,nvira,itsoscf,maxsoscf, &
&                                nao,nmo,sodispmax,idis,nproc2,myrank2,mpi_comm2)
                if(noccb /= 0 ) &
&                 call soscfupdate(cmob,sodispb,work,work2,noccb,nvirb,itsoscf,maxsoscf, &
&                                  nao,nmo,sodispmax,idis,nproc2,myrank2,mpi_comm2)
                itsoscf= itsoscf+1
              else
                call diagfock(focka,work,ortho,cmoa,work2,eigena,idis,nproc2,myrank2,mpi_comm2)
                call diagfock(fockb,work,ortho,cmob,work2,eigenb,idis,nproc2,myrank2,mpi_comm2)
              endif
            endif
!
! Quadratically convergent SCF method
!
          case('QC')
            if((itqc == 0).or.(convqc)) then
              call diagfock(focka,work,ortho,cmoa,work2,eigena,idis,nproc2,myrank2,mpi_comm2)
              call diagfock(fockb,work,ortho,cmob,work2,eigenb,idis,nproc2,myrank2,mpi_comm2)
              itqc= itqc+1
            elseif((itqc == maxqc)) then
              call diagfock(focka,work,ortho,cmoa,work2,eigena,idis,nproc2,myrank2,mpi_comm2)
              call diagfock(fockb,work,ortho,cmob,work2,eigenb,idis,nproc2,myrank2,mpi_comm2)
              itqc= 1
            else
              call uhfqc(focka,fockb,cmoa,cmob,dmax,qcgmna,qcgmnb,qcvec, &
&                        qcmat,qcmatsave,qceigen,overlap,xint, &
&                        work,work2,work3,one,nao,nmo,nocca,noccb,nvira,nvirb,nshell, &
&                        maxdim,maxqcdiag,maxqcdiagsub,threshqc, &
&                        nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
              itqc= itqc+1
            endif
        end select
        call cpu_time(time3)
!
! Copy previous density matrix and calculate new density matrix
!
        call calcudmtrx(cmoa,cmob,dmtrxa,dmtrxb,work,nao,neleca,nelecb)
        call ddiff(dmtrxa,dmtrxpreva,work(1),nao3,diffmaxa)
        call ddiff(dmtrxb,dmtrxprevb,work(nao3+1),nao3,diffmaxb)
        diffmax= diffmaxa+diffmaxb
!
        select case(scfconv)
          case('DIIS')
            if(extrap.and.(itdiis==0)) then
              itsub= itextra
            else
              itsub= itdiis
            endif
            if(master) &
&             write(*,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,errmax
          case('SOSCF')
            itsub= itsoscf
            if(master) &
&             write(*,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,sogradmax
          case('QC')
            itsub= itqc
            if(master) &
&             write(*,'(1x,i3,2x,i3,2(1x,f17.9),1f17.9)')iter,itsub,escf,deltae,diffmax
        end select
!
! Check SCF convergence
!
        select case(scfconv)
          case('DIIS')
            if(diffmax.lt.dconv) exit
            if(itdiis >= maxdiis) itdiis= 0
          case('SOSCF')
            if((diffmax.lt.dconv).and.(convsoscf)) exit
            if((diffmax.lt.dconv).and.(itsoscf==1)) exit
            if((diffmax.lt.dconv).and.(.not.convsoscf)) convsoscf=.true.
            if(itsoscf >= maxsoscf) itsoscf= 0
          case('QC')
            if((diffmax < dconv).and.(convqc)) exit
            if((diffmax < dconv).and.(itqc == 1)) exit
            if((diffmax < dconv).and.(.not.convqc)) convqc=.true.
            if((diffmax >= dconv).and.(convqc)) then
              itqc= 0
              convqc=.false.
            endif
        end select
!
        if(iter == maxiter) then
          if(master) then
            write(*,'(" SCF did not converge.")')
            call iabort
          endif
        endif
        call dcopy(nao3,dmtrxa,1,dmtrxpreva,1)
        call dcopy(nao3,dmtrxb,1,dmtrxprevb,1)
        call dcopy(nao3,work(1),1,dmtrxa,1)
        call dcopy(nao3,work(nao3+1),1,dmtrxb,1)
        call cpu_time(time4)
        if(master.and.(iprint >= 3)) write(*,'(10x,6f8.3)')time2-time1,time3-time2,time4-time3
      enddo
!
      if(master) then
        write(*,'(" -----------------------------------------")')
        write(*,'("    SCF Converged.")')
        write(*,'("    UHF Energy = ",f17.9," a.u.")')escf
        write(*,'(" -----------------------------------------"/)')
      endif
!
! Calculate spin expectation values
!
      call calcspin(sz,s2,dmtrxa,dmtrxb,overlap,work,work2,work3,neleca,nelecb,nao, &
&                   idis,nproc2,myrank2,mpi_comm2)
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
      select case(scfconv)
        case('DIIS')
          deallocate(focka,fockb,fockpreva,fockprevb, &
&                    dmtrxpreva,dmtrxprevb,dmax, &
&                    fockdiisa,fockdiisb,errdiisa,&
&                    errdiisb,diismtrx, &
&                    work,work2,work3)
          call memunset(nao3*8+nshell3+isize3*maxdiis*2+isize1+isize2*maxdiis*2 &
&                      +maxdiis*(maxdiis+1)/2+isize2)
        case('SOSCF')
          deallocate(focka,fockb,fockpreva,fockprevb, &
&                    dmtrxpreva,dmtrxprevb,dmax, &
&                    hstarta,hstartb, &
&                    sograda,sogradb, &
&                    sodispa,sodispb, &
&                    sovecya,sovecyb, &
&                    work,work2,work3)
          call memunset(nao3*8+nshell3+nocca*nvira*3*maxsoscf+noccb*nvirb*3*maxsoscf &
&                      +isize1+isize2)
        case('QC')
          deallocate(focka,fockb,fockpreva,fockprevb, &
&                    dmtrxpreva,dmtrxprevb,dmax, &
&                    qcvec, &
&                    qcmat,qcmatsave, &
&                    qceigen,qcgmna,qcgmnb, &
&                    work,work2,work3)
          call memunset(nao2*4+nao3*8+(nocca*nvira+noccb*nvirb+1)*(maxqcdiagsub+1)*2 &
&                      +maxqcdiagsub**2+maxqcdiagsub*(maxqcdiagsub+1)/2+maxqcdiagsub)
      end select
!
      return
end


!-----------------------------------------------------------------------------------------------
  subroutine formufock(fock1,fock2,fock3,dmtrxa,dmtrxb,dmax,xint,maxdim,nproc,myrank,mpi_comm)
!-----------------------------------------------------------------------------------------------
!
! Driver of Fock matrix formation from two-electron intgrals
!
! In  : dmtrxa (Alpha density matrix)
!       dmtrxb (Beta density matrix)
! Out : fock1  (Alpha Fock matrix)
!       fock2  (Beta Fock matrix)
!       fock3  (Work space)
!
      use modbasis, only : nshell, nao
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: maxdim, nproc, myrank, mpi_comm
      integer :: ijsh, ish, jsh, ksh, lsh, ij, kl, ik, il, jk, jl
      integer :: ii, jj, kk, kstart, ishcheck
      integer(8) :: ncount, icount(nshell)
      real(8),parameter :: zero=0.0D+00, two=2.0D+00, four=4.0D+00
      real(8),intent(in) :: dmtrxa(nao*(nao+1)/2), dmtrxb(nao*(nao+1)/2)
      real(8),intent(in) :: dmax(nshell*(nshell+1)/2), xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: fock1(nao*(nao+1)/2), fock2(nao*(nao+1)/2), fock3(nao*(nao+1)/2)
      real(8) :: xijkl, denmax, twoeri(maxdim**4), denmax1
      integer :: last, ltmp(nshell), lnum, ll
!
      fock2(:)= zero
      fock3(:)= zero
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
!$OMP private(ijsh,jsh,ksh,lsh,ij,kl,ik,il,jk,jl,xijkl,denmax1,denmax,twoeri,jj,kk, &
!$OMP kstart,last,ltmp,lnum,ll) firstprivate(ish,ii) reduction(+:fock2,fock3)
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
        jj= jsh*(jsh-1)/2
        kstart=mod(icount(ish)-ish*(jsh-1),nproc)+1
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
!         do lsh= 1,ksh
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
!$OMP end parallel do
!
      do ii= 1,nao
        fock2(ii*(ii+1)/2)= fock2(ii*(ii+1)/2)*two
        fock3(ii*(ii+1)/2)= fock3(ii*(ii+1)/2)*two
      enddo
!
      call para_allreducer(fock2,fock1,nao*(nao+1)/2,mpi_comm)
      call para_allreducer(fock3,fock2,nao*(nao+1)/2,mpi_comm)
      return
end


!-------------------------------------------------------------------------------
  subroutine ufockeri(focka,fockb,dmtrxa,dmtrxb,twoeri,ish,jsh,ksh,lsh,maxdim)
!-------------------------------------------------------------------------------
!
! Form unrestricted Fock matrix from two-electron intgrals
!
      use modbasis, only : nao, mbf, locbf
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: ish, jsh, ksh, lsh, maxdim
      integer :: nbfi, nbfj, nbfk, nbfl
      integer :: locbfi, locbfj, locbfk, locbfl, jmax, lmax, i, j, k, l, ij, kl
      integer :: nij, nkl, nik, nil, njk, njl
      integer :: iloc, jloc, kloc, lloc, iloc2, jloc2, kloc2, lloc2, jloc0, kloc0
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
      jloc0= locbfj*(locbfj-1)/2
      kloc0= locbfk*(locbfk-1)/2
!
      if(.not.ieqk)then
        do i= 1,nbfi
          iloc= locbfi+i
          iloc2= iloc*(iloc-1)/2
          jloc2= jloc0
          if(ieqj) jmax= i
          do j= 1,jmax
            jloc= locbfj+j
            jloc2= jloc2+jloc-1
            kloc2= kloc0
            nij= iloc2+jloc
            do k= 1,nbfk
              kloc= locbfk+k
              kloc2= kloc2+kloc-1
              nik= iloc2+kloc
              njk= jloc2+kloc
              if(jloc.lt.kloc) njk= kloc2+jloc
              if(keql) lmax= k
              do l= 1,lmax
                val= twoeri(l,k,j,i)
                if(abs(val).lt.cutint2) cycle
                lloc= locbfl+l
                nkl= kloc2+lloc
                nil= iloc2+lloc
                njl= jloc2+lloc
                if(jloc.lt.lloc) njl= lloc*(lloc-1)/2+jloc
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
            enddo
          enddo
        enddo
      else
        do i= 1,nbfi
          iloc= locbfi+i
          iloc2= iloc*(iloc-1)/2
          jloc2= jloc0
          if(ieqj) jmax= i
          do j= 1,jmax
            jloc= locbfj+j
            jloc2= jloc2+jloc-1
            kloc2= kloc0
            ij= ij+1
            kl= 0
     kloop: do k= 1,nbfk
              kloc= locbfk+k
              kloc2= kloc2+kloc-1
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
      endif
!
      return
end


!-----------------------------------------------------------------------------------------
  subroutine calcudft(h1mtrx,cmoa,cmob,ortho,overlap,dmtrxa,dmtrxb,xint,eigena,eigenb, &
&                     nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!-----------------------------------------------------------------------------------------
!
! Driver of unrestricted DFT calculation
!
! In  : h1mtrx  (1-electron integral matrix)
!       ortho   (Orthogonalization matrix)
!       overlap (Overlap integral matrix)
! Out : dmtrxa  (Alpha density matrix)
!       dmtrxb  (Beta density matrix)
!       xint    ((ij|ij) integral matrix)
!       eigena  (Alpha MO energy)
!       eigenb  (Beta MO energy)
! Inout : cmoa  (Alpha MO coefficient matrix)
!         cmob  (Beta MO coefficient matrix)
!
      use modparallel, only : master
      use moddft, only : idftex, idftcor, nrad, nleb, hfexchange
      use modatom, only : atomrad
      use modbasis, only : nshell, nao, mtype
      use modmolecule, only : neleca, nelecb, nmo, natom, numatomic
      use modscf, only : maxiter, fdiff, dconv, maxdiis, maxsoscf, maxqc, maxqcdiag, &
&                        maxqcdiagsub, scfconv, extrap
      use modenergy, only : enuc, escf, escfe
      use modthresh, only : threshsoscf, threshqc, cutint2, threshex, threshover, threshdiis, &
&                           threshweight, threshrho, threshdfock, threshdftao
      use modprint, only : iprint
      use modunit, only : tobohr
      use modwarn, only : nwarn
      implicit none
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2, mpi_comm1, mpi_comm2
      integer :: nao2, nao3, nshell3, maxdim, maxfunc(0:6), numwork, iter, itsub, itdiis
      integer :: itextra, itsoscf, itqc, nocca, nvira, noccb, nvirb
      integer :: idis(nproc2,14), isize1, isize2, isize3, iatom
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00
      real(8),parameter :: small=1.0D-10
      real(8),intent(in) :: h1mtrx(nao*(nao+1)/2), ortho(nao*nao), overlap(nao*(nao+1)/2)
      real(8),intent(out) :: dmtrxa(nao*(nao+1)/2), dmtrxb(nao*(nao+1)/2)
      real(8),intent(out) :: xint(nshell*(nshell+1)/2), eigena(nao), eigenb(nao)
      real(8),intent(inout) :: cmoa(nao*nao), cmob(nao*nao)
      real(8),allocatable :: focka(:), fockb(:), fockpreva(:), fockprevb(:)
      real(8),allocatable :: dmtrxpreva(:), dmtrxprevb(:), dmax(:), work(:)
      real(8),allocatable :: fockdiisa(:), fockdiisb(:), errdiisa(:), errdiisb(:)
      real(8),allocatable :: diismtrx(:), work2(:), work3(:)
      real(8),allocatable :: hstarta(:), hstartb(:), sograda(:,:), sogradb(:,:)
      real(8),allocatable :: sodispa(:), sodispb(:), sovecya(:), sovecyb(:)
      real(8),allocatable :: fockda(:), fockdb(:), rad(:), atomvec(:), surface(:)
      real(8),allocatable :: radpt(:), angpt(:), ptweight(:), xyzpt(:), rsqrd(:)
      real(8),allocatable :: vao(:), vmoa(:), vmob(:)
      real(8),allocatable :: qcvec(:), qcmat(:), qcmatsave(:), qceigen(:)
      real(8),allocatable :: qcgmna(:), qcgmnb(:)
      real(8) :: escfprev, diffmax, diffmaxa, diffmaxb, tridot, deltae
      real(8) :: errmax, errmaxa, errmaxb, sogradmax, sogradmaxa, sogradmaxb, sodispmax
      real(8) :: s2, sz, edft, totalelec
      real(8) :: time1, time2, time3, time4
      logical :: convsoscf, convqc
      data maxfunc/1,3,6,10,15,21,28/
!
      nao2= nao*nao
      nao3= nao*(nao+1)/2
      nshell3= nshell*(nshell+1)/2
      nocca= neleca
      nvira= nmo-neleca
      noccb= nelecb
      nvirb= nmo-nelecb
      maxdim=maxfunc(maxval(mtype(1:nshell)))
      numwork= max(nao3*2,(neleca+nelecb)*nao)
!
! Distribute fock and error arrays for DIIS
!
      call distarray(idis,nmo,nao,nao3,nocca,nvira,noccb,nvirb,nproc2)
!
! Set arrays
!
      isize1= max(nao2,idis(myrank2+1,3),idis(myrank2+1,7)*nao,idis(myrank2+1,11)*nao, &
&                 natom*3,nao*2,maxdiis)
      isize2=idis(myrank2+1,3)
      isize3=idis(myrank2+1,5)
      select case(scfconv)
        case('DIIS')
          call memset(nao3*8+nshell3+numwork+natom*5+natom*natom*6+nrad*2+nleb*4 &
&                    +natom*nrad*nleb+nao*4+nocca*4+noccb*4+isize3*maxdiis*2+isize1 &
&                    +isize2*maxdiis*2+maxdiis*(maxdiis+1)/2+idis(myrank2+1,3))
          allocate(focka(nao3),fockb(nao3),fockpreva(nao3),fockprevb(nao3), &
&                  dmtrxpreva(nao3),dmtrxprevb(nao3),dmax(nshell3),work(numwork), &
&                  fockda(nao3),fockdb(nao3),rad(natom),atomvec(5*natom*natom), &
&                  surface(natom*natom),radpt(2*nrad),angpt(4*nleb),ptweight(natom*nrad*nleb), &
&                  xyzpt(3*natom),rsqrd(natom),vao(4*nao),vmoa(4*nocca),vmob(4*noccb), &
&                  fockdiisa(isize3*maxdiis),fockdiisb(isize3*maxdiis),errdiisa(isize2*maxdiis), &
&                  errdiisb(isize2*maxdiis),diismtrx(maxdiis*(maxdiis+1)/2),work2(isize1), &
&                  work3(idis(myrank2+1,3)))
        case('SOSCF')
          call memset(nao3*8+nshell3+numwork+natom*5+natom*natom*6+nrad*2+nleb*4 &
&                    +natom*nrad*nleb+nao*4+nocca*4+noccb*4+isize1+idis(myrank2+1,3) &
&                    +nocca*nvira*3*maxsoscf+noccb*nvirb*3*maxsoscf)
          allocate(focka(nao3),fockb(nao3),fockpreva(nao3),fockprevb(nao3), &
&                  dmtrxpreva(nao3),dmtrxprevb(nao3),dmax(nshell3),work(numwork), &
&                  fockda(nao3),fockdb(nao3),rad(natom),atomvec(5*natom*natom), &
&                  surface(natom*natom),radpt(2*nrad),angpt(4*nleb),ptweight(natom*nrad*nleb), &
&                  xyzpt(3*natom),rsqrd(natom),vao(4*nao),vmoa(4*nocca),vmob(4*noccb), &
&                  work2(isize1),work3(idis(myrank2+1,3)), &
&                  hstarta(nocca*nvira),hstartb(noccb*nvirb), &
&                  sograda(nocca*nvira,maxsoscf),sogradb(noccb*nvirb,maxsoscf), &
&                  sodispa(nocca*nvira*maxsoscf),sodispb(noccb*nvirb*maxsoscf), &
&                  sovecya(nocca*nvira*(maxsoscf-1)),sovecyb(noccb*nvirb*(maxsoscf-1)))
        case('QC')
          isize1= max(isize1,nao2)
          call memset(nao2*3+nao3*8+nshell3+numwork+natom*5+natom*natom*6+nrad*2+nleb*4 &
&                    +natom*nrad*nleb+nao*4+nocca*4+noccb*4+isize1 &
&                    +(nocca*nvira+noccb*nvirb+1)*(maxqcdiagsub+1)*2 &
&                    +maxqcdiagsub*(maxqcdiagsub*3+1)/2+maxqcdiagsub)
          allocate(focka(nao2),fockb(nao2),fockpreva(nao3),fockprevb(nao3), &
&                  dmtrxpreva(nao3),dmtrxprevb(nao3),dmax(nshell3),work(numwork), &
&                  fockda(nao3),fockdb(nao3),rad(natom),atomvec(5*natom*natom), &
&                  surface(natom*natom),radpt(2*nrad),angpt(4*nleb),ptweight(natom*nrad*nleb), &
&                  xyzpt(3*natom),rsqrd(natom),vao(4*nao),vmoa(4*nocca),vmob(4*noccb), &
&                  qcvec((nocca*nvira+noccb*nvirb+1)*(maxqcdiagsub+1)*2), &
&                  qcmat(maxqcdiagsub**2),qcmatsave(maxqcdiagsub*(maxqcdiagsub+1)/2), &
&                  qceigen(maxqcdiagsub),qcgmna(nao3),qcgmnb(nao3),work2(isize1),work3(nao2))
        case default
          if(master) then
            write(*,'(" SCFConv=",a12,"is not supported.")')
            call iabort
          endif
      end select
!
      escfprev= zero
      itdiis =0
      itextra=0
      itsoscf=0
      itqc   =0
      convsoscf=.false.
      convqc   =.false.
!
! Calculate DFT information
!
      call calcatomvec(atomvec,surface)
      call calcradpt(radpt,nrad)
      call calclebpt(angpt,nleb)
      do iatom= 1,natom
        rad(iatom)= atomrad(numatomic(iatom))*tobohr
      enddo
      call calcgridweight(ptweight,rad,radpt,angpt,atomvec,surface,xyzpt,work2,nproc1,myrank1)
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
      call calcschwarzeri(xint,work,maxdim,nproc2,myrank2,mpi_comm2)
!
      if(master) then
        write(*,'(1x,74("-"))')
        write(*,'("   Unrestricted DFT calculation")')
        write(*,'(1x,74("-"))')
        write(*,'("   SCFConv    = ",a8,",  Dconv      =",1p,d9.2,",  MaxIter    =",i9)') &
&                    scfconv, dconv, maxiter
        write(*,'("   Cutint2    =",1p,d9.2,",  ThreshEx   =",d9.2,",  ThreshOver =",d9.2)') &
&                    cutint2, threshex, threshover
        write(*,'("   Nrad       =",i9  ,",  Nleb       =",i9  ,",  ThreshRho  =",1p,d9.2)') &
&                    nrad, nleb, threshrho
        write(*,'("   ThreshDfock=",1p,d9.2,",  Threshdftao=",d9.2,",  ThreshWeight=",d8.2)') &
&                    threshdfock, threshdftao, threshweight
        select case(scfconv)
          case('DIIS')
            write(*,'("   MaxDIIS    =",i9,",  ThreshDIIS =",1p,d9.2)') &
&                    maxdiis, threshdiis
          case('SOSCF')
            write(*,'("   MaxSOSCF   =",i9,",  ThreshSOSCF=",1p,d9.2)') &
&                    maxdiis, threshdiis
          case('QC')
            write(*,'("   MaxQC      =",i9,",  ThreshQC   =",1p,d9.2,",  MaxQCDiag  =",i9)') &
&                    maxqc, threshqc, maxqcdiag
            write(*,'("   MaxQCDiagSub=",i8)') &
&                    maxqcdiagsub
        end select
        write(*,'(1x,74("-"))')
!
        if((nao /= nmo).and.master) then
          write(*,'(/," Warning! Number of MOs is reduced from",i5," to",i5,/)') nao, nmo
          nwarn= nwarn+1
        endif
!
        write(*,'(" ====================")')
        write(*,'("    SCF Iteration")')
        write(*,'(" ====================")')
        select case(scfconv)
          case('DIIS')
            write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                     "Delta Density     DIIS Error")')
          case('SOSCF')
            write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                     "Delta Density    Orbital Grad")')
          case('QC')
            write(*,'(" Iter SubIt   Total Energy      Delta Energy      ", &
&                     "Delta Density")')
        end select
      endif
!
! Start SCF iteration
!
      do iter= 1,maxiter
!
! Calculate maximum density matrix elements
!
        call calcudmax(dmtrxa,dmtrxb,dmax,work,nproc2,myrank2,mpi_comm2)
        call cpu_time(time1)
!
! Calculate exchange-correlation terms
!
        call formufockexcor(fockda,fockdb,focka,edft,totalelec,cmoa,cmob,atomvec,&
&                           radpt,angpt,rad,ptweight,vao,vmoa,vmob,xyzpt,rsqrd, &
&                           work,work(neleca*nao+1),work2,idftex,idftcor,nproc1,myrank1,mpi_comm1)
!
! Calculate two-electron integrals
!
        call formudftfock(focka,fockb,work,dmtrxa,dmtrxb,dmax,xint,maxdim,hfexchange, &
&                         nproc1,myrank1,mpi_comm1)
        call dscal(nao3,half,focka,1)
        call dscal(nao3,half,fockb,1)
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
        select case(scfconv)
!
! DIIS interpolation
!
          case('DIIS')
            call calcdiiserr(focka,dmtrxpreva,overlap,ortho,cmoa,work,work2,errmaxa,nao,nmo, &
&                            idis,nproc2,myrank2,mpi_comm2)
            call calcdiiserr(fockb,dmtrxprevb,overlap,ortho,cmob,work,work2,errmaxb,nao,nmo, &
&                            idis,nproc2,myrank2,mpi_comm2)
            errmax= max(errmaxa,errmaxb)
            if(((itdiis /= 0).or.(errmax <= threshdiis)).and.(errmax > small))then
              itdiis= itdiis+1
              call calcudiis(focka,fockb,errdiisa,errdiisb,fockdiisa,fockdiisb, &
&                            diismtrx,cmoa,cmob,work2,itdiis,nao,maxdiis, &
&                            idis,nproc2,myrank2,mpi_comm2)
            endif
!
! Extrapolate Fock matrix
!
            if(extrap.and.itdiis == 0) then
              call fockextrap(focka,fockdiisa,work,cmoa,dmtrxa,itextra,nao,maxdiis, &
&                             idis,nproc2,myrank2,mpi_comm2)
              call fockextrap(fockb,fockdiisb,work,cmob,dmtrxb,itextra,nao,maxdiis, &
&                             idis,nproc2,myrank2,mpi_comm2)
            endif
!
! Diagonalize Fock matrix
!
            call diagfock(focka,work,ortho,cmoa,work2,eigena,idis,nproc2,myrank2,mpi_comm2)
            call diagfock(fockb,work,ortho,cmob,work2,eigenb,idis,nproc2,myrank2,mpi_comm2)
!
! Approximated Second-order SCF method
!
          case('SOSCF')
            if((itsoscf == 0).or.(convsoscf)) then
              call diagfock(focka,work,ortho,cmoa,work2,eigena,idis,nproc2,myrank2,mpi_comm2)
              call diagfock(fockb,work,ortho,cmob,work2,eigenb,idis,nproc2,myrank2,mpi_comm2)
              sogradmax= zero
              itsoscf= itsoscf+1
            else
              call expand(focka,work,nao)
              call soscfgrad(work,work2,sograda(1,itsoscf),cmoa,nocca,nvira,sogradmaxa,nao, &
&                            idis,nproc2,myrank2,mpi_comm2,1)
              if(noccb /= 0) then
                call expand(fockb,work,nao)
                call soscfgrad(work,work2,sogradb(1,itsoscf),cmob,noccb,nvirb,sogradmaxb,nao, &
&                              idis,nproc2,myrank2,mpi_comm2,2)
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
&                               nocca,noccb,nvira,nvirb,itsoscf,maxsoscf,sodispmax)
                call soscfupdate(cmoa,sodispa,work,work2,nocca,nvira,itsoscf,maxsoscf, &
&                                nao,nmo,sodispmax,idis,nproc2,myrank2,mpi_comm2)
                if(noccb /= 0 ) &
&                 call soscfupdate(cmob,sodispb,work,work2,noccb,nvirb,itsoscf,maxsoscf, &
&                                  nao,nmo,sodispmax,idis,nproc2,myrank2,mpi_comm2)
                itsoscf= itsoscf+1
              else
                call diagfock(focka,work,ortho,cmoa,work2,eigena,idis,nproc2,myrank2,mpi_comm2)
                call diagfock(fockb,work,ortho,cmob,work2,eigenb,idis,nproc2,myrank2,mpi_comm2)
              endif
            endif
!
! Quadratically convergent SCF method
!
          case('QC')
            if((itqc == 0).or.(convqc)) then
              call diagfock(focka,work,ortho,cmoa,work2,eigena,idis,nproc2,myrank2,mpi_comm2)
              call diagfock(fockb,work,ortho,cmob,work2,eigenb,idis,nproc2,myrank2,mpi_comm2)
              itqc= itqc+1
            elseif((itqc == maxqc)) then
              call diagfock(focka,work,ortho,cmoa,work2,eigena,idis,nproc2,myrank2,mpi_comm2)
              call diagfock(fockb,work,ortho,cmob,work2,eigenb,idis,nproc2,myrank2,mpi_comm2)
              itqc= 1
            else
              call uhfqc(focka,fockb,cmoa,cmob,dmax,qcgmna,qcgmnb,qcvec, &
&                        qcmat,qcmatsave,qceigen,overlap,xint, &
&                        work,work2,work3,hfexchange,nao,nmo,nocca,noccb,nvira,nvirb,nshell, &
&                        maxdim,maxqcdiag,maxqcdiagsub,threshqc, &
&                        nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
              itqc= itqc+1
            endif
        end select
        call cpu_time(time3)
!
! Copy previous density matrix and calculate new density matrix
!
        call calcudmtrx(cmoa,cmob,dmtrxa,dmtrxb,work,nao,neleca,nelecb)
        call ddiff(dmtrxa,dmtrxpreva,work(1),nao3,diffmaxa)
        call ddiff(dmtrxb,dmtrxprevb,work(nao3+1),nao3,diffmaxb)
        diffmax= diffmaxa+diffmaxb
!
        select case(scfconv)
          case('DIIS')
            if(extrap.and.(itdiis==0)) then
              itsub= itextra
            else
              itsub= itdiis
            endif
            if(master) &
&             write(*,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,errmax
          case('SOSCF')
            itsub= itsoscf
            if(master) &
&             write(*,'(1x,i3,2x,i3,2(1x,f17.9),2f17.9)')iter,itsub,escf,deltae,diffmax,sogradmax
          case('QC')
            itsub= itqc
            if(master) &
&             write(*,'(1x,i3,2x,i3,2(1x,f17.9),1f17.9)')iter,itsub,escf,deltae,diffmax
        end select
!
! Check SCF convergence
!
        select case(scfconv)
          case('DIIS')
            if(diffmax < dconv) exit
            if(itdiis >= maxdiis) itdiis= 0
          case('SOSCF')
            if((diffmax < dconv).and.(convsoscf)) exit
            if((diffmax < dconv).and.(itsoscf == 1)) exit
            if((diffmax < dconv).and.(.not.convsoscf)) convsoscf=.true.
            if(itsoscf >= maxsoscf) itsoscf= 0
          case('QC')
            if((diffmax < dconv).and.(convqc)) exit
            if((diffmax < dconv).and.(itqc == 1)) exit
            if((diffmax < dconv).and.(.not.convqc)) convqc=.true.
            if((diffmax >= dconv).and.(convqc)) then
              itqc= 0
              convqc=.false.
            endif
        end select
!
        if(iter == maxiter) then
          if(master) then
            write(*,'(" SCF did not converge.")')
            call iabort
          endif
        endif
        call dcopy(nao3,dmtrxa,1,dmtrxpreva,1)
        call dcopy(nao3,dmtrxb,1,dmtrxprevb,1)
        call dcopy(nao3,work(1),1,dmtrxa,1)
        call dcopy(nao3,work(nao3+1),1,dmtrxb,1)
        call cpu_time(time4)
        if(master.and.(iprint >= 3)) write(*,'(10x,6f8.3)')time2-time1,time3-time2,time4-time3
      enddo
!
      if(master) then
        write(*,'(" -----------------------------------------------------------")')
        write(*,'("    SCF Converged.")')
        write(*,'("    DFT Energy = ",f17.9," a.u.")')escf
        write(*,'("    Exchange + Correlation energy = ",f17.9," a.u.")')edft
        write(*,'("    Number of electrons           = ",f17.9)')totalelec
        write(*,'(" -----------------------------------------------------------")')
      endif
!
! Calculate spin expectation values
!
      call calcspin(sz,s2,dmtrxa,dmtrxb,overlap,work,work2,work3,neleca,nelecb,nao, &
&                   idis,nproc2,myrank2,mpi_comm2)
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
      select case(scfconv)
        case('DIIS')
          deallocate(focka,fockb,fockpreva,fockprevb, &
&                    dmtrxpreva,dmtrxprevb,dmax,work, &
&                    fockda,fockdb,rad,atomvec, &
&                    surface,radpt,angpt,ptweight, &
&                    xyzpt,rsqrd,vao,vmoa,vmob, &
&                    fockdiisa,fockdiisb,errdiisa, &
&                    errdiisb,diismtrx,work2, &
&                    work3)
          call memunset(nao3*8+nshell3+numwork+natom*5+natom*natom*6+nrad*2+nleb*4 &
&                      +natom*nrad*nleb+nao*4+nocca*4+noccb*4+isize3*maxdiis*2+isize1 &
&                      +isize2*maxdiis*2+maxdiis*(maxdiis+1)/2+idis(myrank2+1,3))
        case('SOSCF')
          deallocate(focka,fockb,fockpreva,fockprevb, &
&                    dmtrxpreva,dmtrxprevb,dmax,work, &
&                    fockda,fockdb,rad,atomvec, &
&                    surface,radpt,angpt,ptweight, &
&                    xyzpt,rsqrd,vao,vmoa,vmob, &
&                    work2,work3, &
&                    hstarta,hstartb, &
&                    sograda,sogradb, &
&                    sodispa,sodispb, &
&                    sovecya,sovecyb)
          call memunset(nao3*8+nshell3+numwork+natom*5+natom*natom*6+nrad*2+nleb*4 &
&                      +natom*nrad*nleb+nao*4+nocca*4+noccb*4+isize1+idis(myrank2+1,3) &
&                      +nocca*nvira*3*maxsoscf+noccb*nvirb*3*maxsoscf)
        case('QC')
          deallocate(focka,fockb,fockpreva,fockprevb, &
&                    dmtrxpreva,dmtrxprevb,dmax,work, &
&                    fockda,fockdb,rad,atomvec, &
&                    surface,radpt,angpt,ptweight, &
&                    xyzpt,rsqrd,vao,vmoa,vmob, &
&                    qcvec, &
&                    qcmat,qcmatsave, &
&                    qceigen,qcgmna,qcgmnb,work2,work3)
          call memunset(nao2*3+nao3*8+nshell3+numwork+natom*5+natom*natom*6+nrad*2+nleb*4 &
&                      +natom*nrad*nleb+nao*4+nocca*4+noccb*4+isize1 &
&                      +(nocca*nvira+noccb*nvirb+1)*(maxqcdiagsub+1)*2 &
&                      +maxqcdiagsub*(maxqcdiagsub*3+1)/2+maxqcdiagsub)
      end select
      return
end


!-----------------------------------------------------------------------------------------
  subroutine formudftfock(fock1,fock2,fock3,dmtrxa,dmtrxb,dmax,xint,maxdim,hfexchange, &
&                         nproc,myrank,mpi_comm)
!-----------------------------------------------------------------------------------------
!
! Driver of DFT Fock matrix formation from two-electron intgrals
!
! In  : dmtrxa (Alpha density matrix)
!       dmtrxb (Beta density matrix)
! Out : fock1  (Alpha Fock matrix)
!       fock2  (Beta Fock matrix)
!       fock3  (Work space)
!
      use modbasis, only : nshell, nao
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: maxdim, nproc, myrank, mpi_comm
      integer :: ijsh, ish, jsh, ksh, lsh, ij, kl, ik, il, jk, jl
      integer :: ii, jj, kk, kstart, ishcheck
      integer(8) :: ncount, icount(nshell)
      real(8),parameter :: zero=0.0D+00, two=2.0D+00, four=4.0D+00
      real(8),intent(in) :: dmtrxa(nao*(nao+1)/2), dmtrxb(nao*(nao+1)/2)
      real(8),intent(in) :: dmax(nshell*(nshell+1)/2), xint(nshell*(nshell+1)/2), hfexchange
      real(8),intent(out) :: fock1(nao*(nao+1)/2), fock2(nao*(nao+1)/2), fock3(nao*(nao+1)/2)
      real(8) :: xijkl, denmax, twoeri(maxdim**4), denmax1
      integer :: last, ltmp(nshell), lnum, ll
!
      fock2(:)= zero
      fock3(:)= zero
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
!$OMP private(ijsh,jsh,ksh,lsh,ij,kl,ik,il,jk,jl,xijkl,denmax1,denmax,twoeri,jj,kk, &
!$OMP kstart,last,ltmp,lnum,ll) firstprivate(ish,ii) reduction(+:fock2,fock3)
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
        jj= jsh*(jsh-1)/2
        kstart=mod(icount(ish)-ish*(jsh-1),nproc)+1
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
!         do lsh= 1,ksh
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
&                            hfexchange)
          enddo
        enddo
      enddo
!$OMP end parallel do
!
      do ii= 1,nao
        fock2(ii*(ii+1)/2)= fock2(ii*(ii+1)/2)*two
        fock3(ii*(ii+1)/2)= fock3(ii*(ii+1)/2)*two
      enddo
!
      call para_allreducer(fock2,fock1,nao*(nao+1)/2,mpi_comm)
      call para_allreducer(fock3,fock2,nao*(nao+1)/2,mpi_comm)
      return
end


!------------------------------------------------------------------------------------
  subroutine udftfockeri(focka,fockb,dmtrxa,dmtrxb,twoeri,ish,jsh,ksh,lsh,maxdim, &
&                        hfexchange)
!------------------------------------------------------------------------------------
!
! Form unrestricted DFT Fock matrix from two-electron intgrals
!
      use modbasis, only : nao, mbf, locbf
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: ish, jsh, ksh, lsh, maxdim
      integer :: nbfi, nbfj, nbfk, nbfl
      integer :: locbfi, locbfj, locbfk, locbfl, jmax, lmax, i, j, k, l, ij, kl
      integer :: nij, nkl, nik, nil, njk, njl
      integer :: iloc, jloc, kloc, lloc, iloc2, jloc2, kloc2, lloc2, jloc0, kloc0
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
      jloc0= locbfj*(locbfj-1)/2
      kloc0= locbfk*(locbfk-1)/2
!
      if(.not.ieqk)then
        do i= 1,nbfi
          iloc= locbfi+i
          iloc2= iloc*(iloc-1)/2
          jloc2= jloc0
          if(ieqj) jmax= i
          do j= 1,jmax
            jloc= locbfj+j
            jloc2= jloc2+jloc-1
            kloc2= kloc0
            nij= iloc2+jloc
            do k= 1,nbfk
              kloc= locbfk+k
              kloc2= kloc2+kloc-1
              nik= iloc2+kloc
              njk= jloc2+kloc
              if(jloc.lt.kloc) njk= kloc2+jloc
              if(keql) lmax= k
              do l= 1,lmax
                val= twoeri(l,k,j,i)
                if(abs(val).lt.cutint2) cycle
                lloc= locbfl+l
                nkl= kloc2+lloc
                nil= iloc2+lloc
                njl= jloc2+lloc
                if(jloc.lt.lloc) njl= lloc*(lloc-1)/2+jloc
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
            enddo
          enddo
        enddo
      else
        do i= 1,nbfi
          iloc= locbfi+i
          iloc2= iloc*(iloc-1)/2
          jloc2= jloc0
          if(ieqj) jmax= i
          do j= 1,jmax
            jloc= locbfj+j
            jloc2= jloc2+jloc-1
            kloc2= kloc0
            ij= ij+1
            kl= 0
     kloop: do k= 1,nbfk
              kloc= locbfk+k
              kloc2= kloc2+kloc-1
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
      endif
!
      return
end






