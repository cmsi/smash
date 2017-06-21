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
!---------------------------------------------------
  subroutine calcdmtrx(cmo,dmtrx,work,ndim,neleca)
!---------------------------------------------------
!
! Calculate density matrix for closed shell
!
! In  : cmo   (MO coefficients)
!       ndim  (Dimension of basis functions)
!       neleca(Number of alpha electrons)
! Out : dmtrx (Density matrix)
!
      implicit none
      integer,intent(in) :: ndim, neleca
      integer :: i, j, ij
      real(8),parameter :: zero=0.0D+00, two=2.0D+00
      real(8),intent(in) :: cmo(ndim*ndim)
      real(8),intent(out) :: dmtrx(ndim*(ndim+1)/2), work(ndim,ndim)
!
      call dgemm('N','T',ndim,ndim,neleca,two,cmo,ndim,cmo,ndim,zero,work,ndim)
!
!$OMP parallel do schedule(static,1) private(ij)
      do i= 1,ndim
        ij= i*(i-1)/2
        do j= 1,i
          dmtrx(ij+j)= work(j,i)
        enddo
      enddo
!$OMP end parallel do
      return
end


!-------------------------------------------------------------------------
  subroutine calcudmtrx(cmoa,cmob,dmtrxa,dmtrxb,work,ndim,neleca,nelecb)
!-------------------------------------------------------------------------
!
! Calculate density matrix for closed shell
!
! In  : cmoa   (Alpha MO coefficients)
!       cmob   (Beta MO coefficients)
!       ndim   (Dimension of basis functions)
!       neleca (Number of alpha electrons)
!       nelecb (Number of beta electrons)
! Out : dmtrxa (Alpha density matrix)
!       dmtrxb (Beta density matrix)
!
      implicit none
      integer,intent(in) :: ndim, neleca, nelecb
      integer :: i, j, ij
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: cmoa(ndim*ndim), cmob(ndim*ndim)
      real(8),intent(out) :: dmtrxa(ndim*(ndim+1)/2), dmtrxb(ndim*(ndim+1)/2)
      real(8),intent(out) :: work(ndim,ndim)
!
! Alpha electron
!
      call dgemm('N','T',ndim,ndim,neleca,one,cmoa,ndim,cmoa,ndim,zero,work,ndim)
!
!$OMP parallel do schedule(static,1) private(ij)
      do i= 1,ndim
        ij= i*(i-1)/2
        do j= 1,i
          dmtrxa(ij+j)= work(j,i)
        enddo
      enddo
!$OMP end parallel do
!
! Beta electron
!
      if(nelecb /= 0)then
        call dgemm('N','T',ndim,ndim,nelecb,one,cmob,ndim,cmob,ndim,zero,work,ndim)
!
!$OMP parallel do schedule(static,1) private(ij)
        do i= 1,ndim
          ij= i*(i-1)/2
          do j= 1,i
            dmtrxb(ij+j)= work(j,i)
          enddo
        enddo
!$OMP end parallel do
      else
        dmtrxb(:)= zero
      endif
!
      return
end


!-----------------------------------------------------------------
  subroutine calcrdmax(dmtrx,dmax,dmaxtmp,nproc,myrank,mpi_comm)
!-----------------------------------------------------------------
!
! Calculate maximum density matrix element for each shell
!
! In  : dmtrx   (Density matrix elements)
! Out : dmax    (Maximum density matrix elements)
!       dmaxtmp (Work array)
!
      use modbasis, only : nshell, nao, mbf, locbf
      implicit none
      integer,intent(in) :: nproc, myrank, mpi_comm
      integer :: ish, jsh, ijsh, locbfi, locbfj, nbfi, nbfj
      integer :: jnbf, i, j, ii, ij
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: dmtrx(nao*(nao+1)/2)
      real(8),intent(out) :: dmax(nshell*(nshell+1)/2), dmaxtmp(nshell*(nshell+1)/2)
      real(8) :: dtmp
!
      dmaxtmp(:)= zero
!
!$OMP parallel private(locbfi,locbfj,nbfi,nbfj,jsh,ijsh,dtmp,i,j,ii,ij,jnbf)
      do ish= nshell-myrank,1,-nproc
        locbfi= locbf(ish)
        nbfi  = mbf(ish)
!$OMP do
        do jsh= 1,ish
          locbfj= locbf(jsh)
          nbfj  = mbf(jsh)
          ijsh= ish*(ish-1)/2+jsh
          dtmp= zero
          do i= 1,nbfi
            jnbf= nbfj
            if(ish.eq.jsh) jnbf=i
            ii= locbfi+i
            ij= ii*(ii-1)/2+locbfj
            do j= 1,jnbf
              if(abs(dmtrx(ij+j)).GT.dtmp) dtmp= abs(dmtrx(ij+j))
            enddo
          enddo
          dmaxtmp(ijsh)= dtmp
        enddo
!$OMP enddo
      enddo
!$OMP end parallel
!
      call para_allreducer(dmaxtmp,dmax,nshell*(nshell+1)/2,mpi_comm)
      return
end


!-------------------------------------------------------------------------
  subroutine calcudmax(dmtrxa,dmtrxb,dmax,dmaxtmp,nproc,myrank,mpi_comm)
!-------------------------------------------------------------------------
!
! Calculate maximum unrestricted density matrix element for each shell
!
! In  : dmtrxa  (Alpha density matrix elements)
!       dmtrxb  (Beta density matrix elements)
! Out : dmax    (Maximum density matrix elements)
!       dmaxtmp (Work array)
!
      use modbasis, only : nshell, nao, mbf, locbf
      implicit none
      integer,intent(in) :: nproc, myrank, mpi_comm
      integer :: ish, jsh, ijsh, locbfi, locbfj, nbfi, nbfj
      integer :: jnbf, i, j, ii, ij
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: dmtrxa(nao*(nao+1)/2), dmtrxb(nao*(nao+1)/2)
      real(8),intent(out) :: dmax(nshell*(nshell+1)/2), dmaxtmp(nshell*(nshell+1)/2)
      real(8) :: dtmp
!
      dmaxtmp(:)= zero
!
!$OMP parallel private(locbfi,locbfj,nbfi,nbfj,jsh,ijsh,dtmp,i,j,ii,ij,jnbf)
      do ish= nshell-myrank,1,-nproc
        locbfi= locbf(ish)
        nbfi  = mbf(ish)
!$OMP do
        do jsh= 1,ish
          locbfj= locbf(jsh)
          nbfj  = mbf(jsh)
          ijsh= ish*(ish-1)/2+jsh
          dtmp= zero
          do i= 1,nbfi
            jnbf= nbfj
            if(ish.eq.jsh) jnbf=i
            ii= locbfi+i
            ij= ii*(ii-1)/2+locbfj
            do j= 1,jnbf
              if(abs(dmtrxa(ij+j)).GT.dtmp) dtmp= abs(dmtrxa(ij+j))
              if(abs(dmtrxb(ij+j)).GT.dtmp) dtmp= abs(dmtrxb(ij+j))
            enddo
          enddo
          dmaxtmp(ijsh)= dtmp
        enddo
!$OMP enddo
      enddo
!$OMP end parallel
!
      call para_allreducer(dmaxtmp,dmax,nshell*(nshell+1)/2,mpi_comm)
      return
end


!------------------------------------------------------
  subroutine ddiff(dmtrx,dmtrxprev,work,ndim,diffmax)
!------------------------------------------------------
!
! Calculate largest absolute change in density matrix
!
! In  : dmtrx    (Density matrix elements)
!       dmtrxprev(Previous density matrix elements)
!       ndim     (Dimensions of dmtrx and dmtrxprev)
! Out : work     (Change in density matrix)
!       diffmax  (Largest absolute change in density matrix)
!
      implicit none
      integer,intent(in) :: ndim
      integer :: ii
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: dmtrx(ndim), dmtrxprev(ndim)
      real(8),intent(out) :: diffmax, work(ndim)
!
      diffmax= zero
!$OMP parallel do reduction(max:diffmax)
      do ii= 1,ndim
        work(ii)= dmtrx(ii)-dmtrxprev(ii)
        if(abs(work(ii)) > diffmax) diffmax= abs(work(ii))
      enddo
!$OMP end parallel do
      return
end


!-----------------------------------------
  subroutine focksum(fock,fockprev,ndim)
!-----------------------------------------
!
! Sum Fock and previous Fock matrices
!
! In  : fock     (Fock matrix with density matrix difference)
!       fockprev (previous Fock matrix)
! Out : fock     (incremental Fock matrix)
!       fockprev (incremental Fock matrix)
!
      implicit none
      integer,intent(in) :: ndim
      integer :: ii
      real(8),intent(inout) :: fock(ndim), fockprev(ndim)
!
!$OMP parallel do
      do ii= 1,ndim
        fock(ii)= fock(ii)+fockprev(ii)
        fockprev(ii)= fock(ii)
      enddo
!$OMP end parallel do
      return
end


!-------------------------------------------------------------------------------
  subroutine fockextrap(fock,fockwork,work1,work2,work3,itextra,nao,maxdiis, &
&                       idis,nproc,myrank,mpi_comm)
!-------------------------------------------------------------------------------
!
! Extrapolate Fock matrix
!
! In    : itextra (Counter for extrapolation)
! InOut : fock    (Fock matrix)
!         fockwork(Previous Fock and work matrix)
!
      use modparallel, only : master
      implicit none
      integer,intent(in) :: nao, maxdiis, nproc, myrank, mpi_comm, idis(nproc,14) 
      integer,intent(inout) :: itextra
      integer :: num, istart, i, iskip, nao3
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00, four=4.0D+00
      real(8),parameter :: pm10=1.0D-10, pm17=1.0D+17, pm7=1.0D-7, pm19=1.9D+00, pm99=0.99D+00
      real(8),intent(inout) :: fock(nao*(nao+1)/2)
      real(8),intent(inout) :: fockwork(idis(myrank+1,5),maxdiis)
      real(8),intent(out) :: work1(*), work2(*), work3(*)
      real(8) :: tridot
      real(8) :: sp11, sp12, sp13, sp22, sp23, sp33, dp1, dp2, dp3, cosphi, cospsi
      real(8) :: xyz(3), xy1, xy2, xy3, xy4
!
      if(maxdiis < 6) then
        if(master) then
          write(*,'(" Set Maxdiis for more than 6.")')
          call iabort
        endif
      endif
!
      num= idis(myrank+1,5)
      istart= idis(myrank+1,6)+1
      if(itextra.le.2) then
        itextra= itextra+1
        if(num /= 0) call dcopy(num,fock(istart),1,fockwork(1,itextra),1)
        return
      endif
!
      iskip= 0
      nao3= nao*(nao+1)/2 
!
      do i= 1,num
        fockwork(i,4)= fock(istart+i-1)-fockwork(i,3)
        fockwork(i,5)= fockwork(i,3)   -fockwork(i,2)
        fockwork(i,6)= fockwork(i,2)   -fockwork(i,1)
      enddo
      call para_allgathervr(fockwork(1,4),num,work1,idis(1,5),idis(1,6),nproc,mpi_comm)
      call para_allgathervr(fockwork(1,5),num,work2,idis(1,5),idis(1,6),nproc,mpi_comm)
      call para_allgathervr(fockwork(1,6),num,work3,idis(1,5),idis(1,6),nproc,mpi_comm)
      sp11= tridot(work1,work1,nao)
      sp12= tridot(work2,work1,nao)
      sp13= tridot(work3,work1,nao)
      sp22= tridot(work2,work2,nao)
      sp23= tridot(work3,work2,nao)
      sp33= tridot(work3,work3,nao)
      dp1 = sqrt(sp11)
      dp2 = sqrt(sp22)
      dp3 = sqrt(sp33)
      if(dp1.lt.pm10) iskip= 1
      if(dp2.lt.pm10) iskip= 1
      if(dp3.lt.pm10) iskip= 1
!
      cosphi= sp12/(dp1*dp2)
!
      xyz(3)= one/(sp11*sp22-sp12*sp12)
      if(abs(xyz(3)).gt.pm17) iskip= 1
      xyz(1)=(sp13*sp22-sp12*sp23)*xyz(3)
      xyz(2)=(sp23*sp11-sp12*sp13)*xyz(3)
      cospsi= sqrt(xyz(1)*xyz(1)*sp11+xyz(2)*xyz(2)*sp22+two*xyz(1)*xyz(2)*sp12)/dp3
!
      if(cospsi.le.pm7) iskip= 1
!
      xyz(1)= one/xyz(1)
      xyz(2)=-xyz(2)*xyz(1)
!
      xy1= xyz(2)*xyz(2)+four*xyz(1)
      if(xy1.lt.zero) then
        iskip= 1
        xy2= zero
      else
        xy2= abs(xyz(2))+sqrt(xy1)
      endif
      if(xy2.gt.pm19.and.abs(cosphi).le.pm99) iskip= 1
!
      if(iskip.eq.0) then
        if(xy2.gt.pm19) then
          xyz(1)= dp1/(dp2*cosphi-dp1)
          do i= 1,nao3
            fock(i)= fock(i)+xyz(1)*work1(i)
          enddo
        else
          xy3= xyz(1)/(one-xyz(1)-xyz(2))
          xy4=(xyz(1)+xyz(2))/(one-xyz(1)-xyz(2))       
          do i= 1,nao3
            fock(i)= fock(i)+xy3*work2(i)+xy4*work1(i)
          enddo
        endif
        itextra=0
      else
        if(num /= 0) then
          call dcopy(num,fockwork(1,2),1,fockwork(1,1),1)
          call dcopy(num,fockwork(1,3),1,fockwork(1,2),1)
          call dcopy(num,fock(istart) ,1,fockwork(1,3),1)
        endif
        itextra= itextra+1
      endif
      return
end


!--------------------------------------------------------------------------------------
  subroutine calcdiiserr(fock,dmtrx,overlap,ortho,work1,work2,work3,errmax,nao,nmo, &
&                        idis,nproc,myrank,mpi_comm)
!--------------------------------------------------------------------------------------
!
! Calculate error matrix (=FDS-SDF) for Direct Inversion in the Iterative Subspace
! (DIIS) interporation
!
! In    : fock    (Interporated Fock matrix (out))
!         dmtrx   (Density matrix)
!         overlap (Overlap matrix)
!         ortho   (Orthogonalization matrix)
! Out   : work1   (Error matrix)
!
      implicit none
      integer,intent(in) :: nao, nmo, nproc, myrank, mpi_comm, idis(nproc,14)
      integer :: num, istart, i, j
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: fock(nao*(nao+1)/2), dmtrx(nao*(nao+1)/2)
      real(8),intent(in) :: overlap(nao*(nao+1)/2), ortho(nao,nao)
      real(8),intent(out) :: work1(nao,*), work2(*), work3(*), errmax
      real(8) :: diff
!
! Calculate error matrix (=FDS-SDF)
!
      num=idis(myrank+1,1)
      istart=idis(myrank+1,2)+1
      if(num > 0)then
        call expand(overlap,work1,nao)
        call dsymm('L','U',nao,num,one,work1,nao,ortho(1,istart),nao,zero,work2,nao)
        call expand(dmtrx,work1,nao)
        call dsymm('L','U',nao,num,one,work1,nao,work2,nao,zero,work3,nao)
        call expand(fock,work1,nao)
        call dsymm('L','U',nao,num,one,work1,nao,work3,nao,zero,work2,nao)
        call dgemm('T','N',nmo,num,nao,one,ortho,nao,work2,nao,zero,work3,nao)
      endif
      call para_allgathervr(work3,idis(myrank+1,3),work1,idis(1,3),idis(1,4),nproc,mpi_comm)
!
! Next calculation should be parallelized.
!
      errmax= zero
!$OMP parallel do private(diff) reduction(max:errmax)
      do i= 1,nmo
        do j= 1,i
          diff= work1(j,i)-work1(i,j)
          work1(j,i)= diff
          work1(i,j)=-diff
          if(abs(diff) > errmax) errmax= abs(diff)
        enddo
      enddo
!$OMP end parallel do
!
      if(nmo /= nao) then
        do i= 1,nmo
          do j= nmo+1,nao
            work1(j,i)= zero
          enddo
        enddo
      endif
!
      return
end


!----------------------------------------------------------------------------------------
  subroutine calcrdiis(fock,errdiis,fockdiis,diismtrx,work1,work2,itdiis,nao,maxdiis, &
&                      idis,nproc,myrank,mpi_comm)
!----------------------------------------------------------------------------------------
!
! Direct Inversion in the Iterative Subspace (DIIS) interporation for closed-shell
!
! In    : fock    (Fock matrix)
! Out   : fock    (Interporated Fock matrix)
! Inout : fockdiis(Previous Fock matrix)
!         errdiis (DIIS error matrix)
!
      implicit none
      integer,intent(in) :: itdiis, nao, maxdiis, nproc, myrank, mpi_comm, idis(nproc,14)
      integer :: num, istart, i, j, ij, ipiv(maxdiis+1), info
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(inout) :: fock(nao*(nao+1)/2), errdiis(idis(myrank+1,3),maxdiis)
      real(8),intent(inout) :: fockdiis(idis(myrank+1,5),maxdiis)
      real(8),intent(inout) :: diismtrx(maxdiis*(maxdiis+1)/2), work1(nao*nao)
      real(8),intent(out) :: work2(*)
      real(8) :: diisb(maxdiis+1,maxdiis+1), diiscoeff(maxdiis+1), ddot
!
      num=idis(myrank+1,1)
      istart=idis(myrank+1,4)+1
      work2(1:itdiis)= zero
      if(num > 0) then
        call dcopy(num*nao,work1(istart),1,errdiis(1,itdiis),1)
        do i= 1,itdiis
          work2(i)= ddot(num*nao,errdiis(1,i),1,errdiis(1,itdiis),1)
        enddo
      endif
      call para_allreducer(work2,diismtrx(itdiis*(itdiis-1)/2+1),itdiis,mpi_comm)
!
      do i= 1,itdiis
        do j= 1,i
          ij= i*(i-1)/2+j
          diisb(j,i)= diismtrx(ij)
          diisb(i,j)= diismtrx(ij)
        enddo
      enddo
      do j= 1,itdiis
        diisb(itdiis+1,j)=-one
        diisb(j,itdiis+1)=-one
        diiscoeff(j)= zero
      enddo
      diisb(itdiis+1,itdiis+1)= zero
      diiscoeff(itdiis+1)=-one
!
      call dgesv(itdiis+1,1,diisb,maxdiis+1,ipiv,diiscoeff,maxdiis+1,info)
!
      num=idis(myrank+1,5)
      istart=idis(myrank+1,6)+1
      if(num > 0) then
        call dcopy(num,fock(istart),1,fockdiis(1,itdiis),1)
        work1(1:num)= zero
        do i= 1,itdiis
          call daxpy(num,diiscoeff(i),fockdiis(1,i),1,work1,1)
        enddo
      endif
      call para_allgathervr(work1,num,fock,idis(1,5),idis(1,6),nproc,mpi_comm)
!
      return
end


!----------------------------------------------------------------------------
  subroutine calcudiis(focka,fockb,errdiisa,errdiisb,fockdiisa,fockdiisb, &
&                      diismtrx,worka,workb,work2,itdiis,nao,maxdiis, &
&                      idis,nproc,myrank,mpi_comm)
!----------------------------------------------------------------------------
!
! Direct Inversion in the Iterative Subspace (DIIS) interporation for open-shell
!
! In    : focka    (Alpha Fock matrix)
!         fockb    (Beta Fock matrix)
!         worka    (Alpha DIIS error matrix)
!         workb    (Beta DIIS error matrix)
! Out   : focka    (Interporated alpha Fock matrix)
!         fockb    (Interporated beta Fock matrix)
! Inout : fockdiisa(Previous alpha Fock matrix)
!         fockdiisb(Previous beta Fock matrix)
!         errdiisa (History of Alpha DIIS error matrix)
!         errdiisb (History of Beta DIIS error matrix)
!
      implicit none
      integer,intent(in) :: itdiis, nao, maxdiis, nproc, myrank, mpi_comm, idis(nproc,14)
      integer :: num, istart, i, j, ij, ipiv(maxdiis+1), info
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00
      real(8),intent(inout) :: focka(nao*(nao+1)/2), fockb(nao*(nao+1)/2)
      real(8),intent(inout) :: errdiisa(idis(myrank+1,3),maxdiis)
      real(8),intent(inout) :: errdiisb(idis(myrank+1,3),maxdiis)
      real(8),intent(inout) :: fockdiisa(idis(myrank+1,5),maxdiis)
      real(8),intent(inout) :: fockdiisb(idis(myrank+1,5),maxdiis)
      real(8),intent(inout) :: diismtrx(maxdiis*(maxdiis+1)/2)
      real(8),intent(inout) :: worka(nao*nao), workb(nao*nao)
      real(8),intent(out) :: work2(*)
      real(8) :: diisb(maxdiis+1,maxdiis+1), diiscoeff(maxdiis+1), ddot
!
      num=idis(myrank+1,1)
      istart=idis(myrank+1,4)+1
      work2(1:itdiis)= zero
      if(num > 0) then
        call dcopy(num*nao,worka(istart),1,errdiisa(1,itdiis),1)
        do i= 1,itdiis
          work2(i)= ddot(num*nao,errdiisa(1,i),1,errdiisa(1,itdiis),1)
        enddo
        call dcopy(num*nao,workb(istart),1,errdiisb(1,itdiis),1)
        do i= 1,itdiis
          work2(i)= work2(i)+ddot(num*nao,errdiisb(1,i),1,errdiisb(1,itdiis),1)
        enddo
      endif
      call para_allreducer(work2,diismtrx(itdiis*(itdiis-1)/2+1),itdiis,mpi_comm)
!
      do i= 1,itdiis
        do j= 1,i
          ij= i*(i-1)/2+j
          diisb(j,i)= diismtrx(ij)*half
          diisb(i,j)= diismtrx(ij)*half
        enddo
      enddo
      do j= 1,itdiis
        diisb(itdiis+1,j)=-one
        diisb(j,itdiis+1)=-one
        diiscoeff(j)= zero
      enddo
      diisb(itdiis+1,itdiis+1)= zero
      diiscoeff(itdiis+1)=-one
!
      call dgesv(itdiis+1,1,diisb,maxdiis+1,ipiv,diiscoeff,maxdiis+1,info)
!
      num=idis(myrank+1,5)
      istart=idis(myrank+1,6)+1
!
! Interporate alpha Fock matrix
!
      if(num > 0) then
        call dcopy(num,focka(istart),1,fockdiisa(1,itdiis),1)
        worka(1:num)= zero
        do i= 1,itdiis
          call daxpy(num,diiscoeff(i),fockdiisa(1,i),1,worka,1)
        enddo
      endif
      call para_allgathervr(worka,num,focka,idis(1,5),idis(1,6),nproc,mpi_comm)
!
! Interporate beta Fock matrix
!
      if(num > 0) then
        call dcopy(num,fockb(istart),1,fockdiisb(1,itdiis),1)
        workb(1:num)= zero
        do i= 1,itdiis
          call daxpy(num,diiscoeff(i),fockdiisb(1,i),1,workb,1)
        enddo
      endif
      call para_allgathervr(workb,num,fockb,idis(1,5),idis(1,6),nproc,mpi_comm)
!
      return
end


!------------------------------------------------------------------------
  subroutine soscfgrad(work,work2,sograd,cmo,nocc,nvir,sogradmax,nao, &
&                      idis,nproc,myrank,mpi_comm,itype)
!------------------------------------------------------------------------
!
! Calculate orbital gradient <occ|F|vir>
!
! In  : work      (Fock matrix (filled in upper triangle))
!       cmo       (MO coeffient matrix)
!       itype     (1:Alpha electron, 2:Beta electron)
! Out : sograd    (SOSCF gradient matrix)
!       sogradmax (Maximum SOSCF gradient
!
      implicit none
      integer,intent(in) :: nocc, nvir, nao, itype, nproc, myrank, mpi_comm, idis(nproc,14)
      integer :: numwork, num, istart, isomax, idamax
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: cmo(nao,nao)
      real(8),intent(out) :: sograd(nocc*nvir), sogradmax
      real(8),intent(inout) :: work(nao*nao), work2(*)
!
! Alpha electron gradient
!
      if(itype == 1) then
        num= idis(myrank+1,7)
        istart= idis(myrank+1,8)+1+nocc
        if(num > 0) then
          call dsymm('L','U',nao,num,one,work,nao,cmo(1,istart),nao,zero,work2,nao)
          call dgemm('T','N',nocc,num,nao,one,cmo,nao,work2,nao,zero,work,nocc)
        endif
        numwork= idis(myrank+1,9)
        call para_allgathervr(work,numwork,sograd,idis(1,9),idis(1,10),nproc,mpi_comm)
!
! Beta electron gradient
!
      else
        num= idis(myrank+1,11)
        istart= idis(myrank+1,12)+1+nocc
        if(num > 0) then
          call dsymm('L','U',nao,num,one,work,nao,cmo(1,istart),nao,zero,work2,nao)
          call dgemm('T','N',nocc,num,nao,one,cmo,nao,work2,nao,zero,work,nocc)
        endif
!
        numwork= idis(myrank+1,13)
        call para_allgathervr(work,numwork,sograd,idis(1,13),idis(1,14),nproc,mpi_comm)
      endif
!
! Obtain maximum sograd value
!
      isomax= idamax(nocc*nvir,sograd,1)
      sogradmax= abs(sograd(isomax))
      return
end


!----------------------------------------------------
  subroutine soscfinith(hstart,eigen,nocc,nvir,nao)
!----------------------------------------------------
!
! Set initial inverse Hessian matrix for SOSCF
!
      implicit none
      integer,intent(in) :: nocc, nvir, nao
      integer :: imo, jmo
      real(8),parameter :: one=1.0D+00, small=5.0D-03, p200=2.0D+02
      real(8),intent(in) :: eigen(nao)
      real(8),intent(out) :: hstart(nocc,nvir)
      real(8) :: dele
!
!$OMP parallel do private(dele)
      do imo= 1,nvir
        do jmo= 1,nocc
          dele= eigen(imo+nocc)-eigen(jmo)
          if(abs(dele) > small) then
            hstart(jmo,imo)= one/dele
          else
            hstart(jmo,imo)= p200
          endif
        enddo
      enddo
!$OMP end parallel do
      return
end


!-----------------------------------------------------------------------------------------
  subroutine soscfnewh(hstart,sograd,sodisp,sovecy,nocc,nvir,itsoscf,maxsoscf,sodispmax)
!-----------------------------------------------------------------------------------------
!
! Update inverse Hessian matrix for approximated SOSCF method
!
      implicit none
      integer,intent(in) :: nocc, nvir, itsoscf, maxsoscf
      integer :: ijmo, it, idamax
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: hstart(nocc*nvir)
      real(8),intent(inout) :: sograd(nocc*nvir,maxsoscf), sodisp(nocc*nvir,maxsoscf)
      real(8),intent(inout) :: sovecy(nocc*nvir,maxsoscf-1)
      real(8),intent(out) :: sodispmax
      real(8) :: s1, s2, s3, s4, s5, s6, s1s2, t1, t2, t3, t4
!
! Initialize displacement vector
!
!$OMP parallel do
      do ijmo= 1,nocc*nvir
        sodisp(ijmo,itsoscf)= hstart(ijmo)*sograd(ijmo,itsoscf)
      enddo
!$OMP end parallel do
!
      if(itsoscf >= 2) then
!$OMP parallel do
        do ijmo= 1,nocc*nvir
          sograd(ijmo,itsoscf-1)= sograd(ijmo,itsoscf)-sograd(ijmo,itsoscf-1)
          sovecy(ijmo,itsoscf-1)=-hstart(ijmo)*sograd(ijmo,itsoscf-1)
        enddo
!$OMP end parallel do
!
        do it= 1,itsoscf-2
          s1= zero
          s2= zero
          s3= zero
          s4= zero
          s5= zero
          s6= zero
!$OMP parallel do reduction(+:s1,s2,s3,s4,s5,s6)
          do ijmo= 1,nocc*nvir
            s1= s1+sodisp(ijmo,it)*sograd(ijmo,it)
            s2= s2+sograd(ijmo,it)*sovecy(ijmo,it)
            s3= s3+sodisp(ijmo,it)*sograd(ijmo,itsoscf)
            s4= s4+sovecy(ijmo,it)*sograd(ijmo,itsoscf)
            s5= s5+sodisp(ijmo,it)*sograd(ijmo,itsoscf-1)
            s6= s6+sovecy(ijmo,it)*sograd(ijmo,itsoscf-1)
          enddo
!$OMP end parallel do
          s1= one/s1
          s2= one/s2
          s1s2= one+s1/s2
          t2= s1*s3
          t4= s1*s5
          t1= s1s2*t2-s1*s4
          t3= s1s2*t4-s1*s6
!$OMP parallel do
          do ijmo= 1,nocc*nvir
            sodisp(ijmo,itsoscf)= sodisp(ijmo,itsoscf)-t1*sodisp(ijmo,it) &
&                                                     +t2*sovecy(ijmo,it)
            sovecy(ijmo,itsoscf-1)= sovecy(ijmo,itsoscf-1)+t3*sodisp(ijmo,it) &
&                                                         -t4*sovecy(ijmo,it)
          enddo
!$OMP end parallel do
        enddo
!
        s1= zero
        s2= zero
        s3= zero
        s4= zero
!$OMP parallel do reduction(+:s1,s2,s3,s4)
        do ijmo= 1,nocc*nvir
          s1= s1+sodisp(ijmo,itsoscf-1)*sograd(ijmo,itsoscf-1)
          s2= s2+sograd(ijmo,itsoscf-1)*sovecy(ijmo,itsoscf-1)
          s3= s3+sodisp(ijmo,itsoscf-1)*sograd(ijmo,itsoscf)
          s4= s4+sovecy(ijmo,itsoscf-1)*sograd(ijmo,itsoscf)
        enddo
!$OMP end parallel do
        s1= one/s1
        s2= one/s2
        s1s2= one+s1/s2
        t1= s1s2*s1*s3-s1*s4
        t2= s1*s3
!$OMP parallel do
        do ijmo= 1,nocc*nvir
          sodisp(ijmo,itsoscf)= sodisp(ijmo,itsoscf)-t1*sodisp(ijmo,itsoscf-1) &
&                                                   +t2*sovecy(ijmo,itsoscf-1)
        enddo
!$OMP end parallel do
      endif
!
      ijmo= idamax(nocc*nvir,sodisp(1,itsoscf),1)
      sodispmax= sodisp(ijmo,itsoscf)
!
      return
end


!-------------------------------------------------------------------------------------------
  subroutine soscfunewh(hstarta,hstartb,sograda,sogradb,sodispa,sodispb,sovecya,sovecyb, &
&                       nocca,noccb,nvira,nvirb,itsoscf,maxsoscf,sodispmax)
!-------------------------------------------------------------------------------------------
!
! Update inverse Hessian matrix for approximated SOSCF method (unrestricted version)
!
      implicit none
      integer,intent(in) :: nocca, noccb, nvira, nvirb, itsoscf, maxsoscf
      integer :: ijmo, it, idamax
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: hstarta(nocca*nvira), hstartb(noccb*nvirb)
      real(8),intent(inout) :: sograda(nocca*nvira,maxsoscf), sogradb(noccb*nvirb,maxsoscf)
      real(8),intent(inout) :: sodispa(nocca*nvira,maxsoscf), sodispb(noccb*nvirb,maxsoscf)
      real(8),intent(inout) :: sovecya(nocca*nvira,maxsoscf-1), sovecyb(noccb*nvirb,maxsoscf-1)
      real(8),intent(out) :: sodispmax
      real(8) :: s1, s2, s3, s4, s5, s6, s1s2, t1, t2, t3, t4
!
! Initialize displacement vector
!
!$OMP parallel do
      do ijmo= 1,nocca*nvira
        sodispa(ijmo,itsoscf)= hstarta(ijmo)*sograda(ijmo,itsoscf)
      enddo
!$OMP end parallel do
!$OMP parallel do
      do ijmo= 1,noccb*nvirb
          sodispb(ijmo,itsoscf)= hstartb(ijmo)*sogradb(ijmo,itsoscf)
      enddo
!$OMP end parallel do
!
      if(itsoscf >= 2) then
!$OMP parallel do
        do ijmo= 1,nocca*nvira
          sograda(ijmo,itsoscf-1)= sograda(ijmo,itsoscf)-sograda(ijmo,itsoscf-1)
          sovecya(ijmo,itsoscf-1)=-hstarta(ijmo)*sograda(ijmo,itsoscf-1)
        enddo
!$OMP end parallel do
!$OMP parallel do
        do ijmo= 1,noccb*nvirb
          sogradb(ijmo,itsoscf-1)= sogradb(ijmo,itsoscf)-sogradb(ijmo,itsoscf-1)
          sovecyb(ijmo,itsoscf-1)=-hstartb(ijmo)*sogradb(ijmo,itsoscf-1)
        enddo
!$OMP end parallel do
!
        do it= 1,itsoscf-2
          s1= zero
          s2= zero
          s3= zero
          s4= zero
          s5= zero
          s6= zero
!$OMP parallel reduction(+:s1,s2,s3,s4,s5,s6) 
!$OMP do
          do ijmo= 1,nocca*nvira
            s1= s1+sodispa(ijmo,it)*sograda(ijmo,it)
            s2= s2+sograda(ijmo,it)*sovecya(ijmo,it)
            s3= s3+sodispa(ijmo,it)*sograda(ijmo,itsoscf)
            s4= s4+sovecya(ijmo,it)*sograda(ijmo,itsoscf)
            s5= s5+sodispa(ijmo,it)*sograda(ijmo,itsoscf-1)
            s6= s6+sovecya(ijmo,it)*sograda(ijmo,itsoscf-1)
          enddo
!$OMP end do
!$OMP do 
          do ijmo= 1,noccb*nvirb
            s1= s1+sodispb(ijmo,it)*sogradb(ijmo,it)
            s2= s2+sogradb(ijmo,it)*sovecyb(ijmo,it)
            s3= s3+sodispb(ijmo,it)*sogradb(ijmo,itsoscf)
            s4= s4+sovecyb(ijmo,it)*sogradb(ijmo,itsoscf)
            s5= s5+sodispb(ijmo,it)*sogradb(ijmo,itsoscf-1)
            s6= s6+sovecyb(ijmo,it)*sogradb(ijmo,itsoscf-1)
          enddo
!$OMP end do
!$OMP end parallel
          s1= one/s1
          s2= one/s2
          s1s2= one+s1/s2
          t2= s1*s3
          t4= s1*s5
          t1= s1s2*t2-s1*s4
          t3= s1s2*t4-s1*s6
!$OMP parallel do
          do ijmo= 1,nocca*nvira
            sodispa(ijmo,itsoscf)= sodispa(ijmo,itsoscf)-t1*sodispa(ijmo,it) &
&                                                       +t2*sovecya(ijmo,it)
            sovecya(ijmo,itsoscf-1)= sovecya(ijmo,itsoscf-1)+t3*sodispa(ijmo,it) &
&                                                           -t4*sovecya(ijmo,it)
          enddo
!$OMP end parallel do
!$OMP parallel do
          do ijmo= 1,noccb*nvirb
            sodispb(ijmo,itsoscf)= sodispb(ijmo,itsoscf)-t1*sodispb(ijmo,it) &
&                                                       +t2*sovecyb(ijmo,it)
            sovecyb(ijmo,itsoscf-1)= sovecyb(ijmo,itsoscf-1)+t3*sodispb(ijmo,it) &
&                                                           -t4*sovecyb(ijmo,it)
          enddo
!$OMP end parallel do
        enddo
!
        s1= zero
        s2= zero
        s3= zero
        s4= zero
!$OMP parallel reduction(+:s1,s2,s3,s4)
!$OMP do
        do ijmo= 1,nocca*nvira
          s1= s1+sodispa(ijmo,itsoscf-1)*sograda(ijmo,itsoscf-1)
          s2= s2+sograda(ijmo,itsoscf-1)*sovecya(ijmo,itsoscf-1)
          s3= s3+sodispa(ijmo,itsoscf-1)*sograda(ijmo,itsoscf)
          s4= s4+sovecya(ijmo,itsoscf-1)*sograda(ijmo,itsoscf)
        enddo
!$OMP end do
!$OMP do
        do ijmo= 1,noccb*nvirb
          s1= s1+sodispb(ijmo,itsoscf-1)*sogradb(ijmo,itsoscf-1)
          s2= s2+sogradb(ijmo,itsoscf-1)*sovecyb(ijmo,itsoscf-1)
          s3= s3+sodispb(ijmo,itsoscf-1)*sogradb(ijmo,itsoscf)
          s4= s4+sovecyb(ijmo,itsoscf-1)*sogradb(ijmo,itsoscf)
        enddo
!$OMP end do
!$OMP end parallel
        s1= one/s1
        s2= one/s2
        s1s2= one+s1/s2
        t1= s1s2*s1*s3-s1*s4
        t2= s1*s3
!$OMP parallel do
        do ijmo= 1,nocca*nvira
          sodispa(ijmo,itsoscf)= sodispa(ijmo,itsoscf)-t1*sodispa(ijmo,itsoscf-1) &
&                                                     +t2*sovecya(ijmo,itsoscf-1)
        enddo
!$OMP end parallel do
!$OMP parallel do
        do ijmo= 1,noccb*nvirb
          sodispb(ijmo,itsoscf)= sodispb(ijmo,itsoscf)-t1*sodispb(ijmo,itsoscf-1) &
&                                                     +t2*sovecyb(ijmo,itsoscf-1)
        enddo
!$OMP end parallel do
!
!ishimura
      t1=dot_product(sodispa(1:nocca*nvira,itsoscf),sodispa(1:nocca*nvira,itsoscf)) &
&       +dot_product(sodispb(1:noccb*nvirb,itsoscf),sodispb(1:noccb*nvirb,itsoscf))
      t1=sqrt(t1)
      ijmo= idamax(nocca*nvira,sodispa(1,itsoscf),1)
      it= idamax(noccb*nvirb,sodispb(1,itsoscf),1)
      t4=max(sodispa(ijmo,itsoscf),sodispb(it,itsoscf))
     
!  if(t1.gt.1.0)then
!    call dscal(nocca*nvira,1.0/t1,sodispa(1,itsoscf),1)
!    call dscal(noccb*nvirb,1.0/t1,sodispb(1,itsoscf),1)
!  endif
!ishimura
      endif
!
      ijmo= idamax(nocca*nvira,sodispa(1,itsoscf),1)
      sodispmax= abs(sodispa(ijmo,itsoscf))
      ijmo= idamax(noccb*nvirb,sodispb(1,itsoscf),1)
      sodispmax= max(sodispmax,abs(sodispb(ijmo,itsoscf)))
!
      return
end


!-----------------------------------------------------------------------------
  subroutine soscfupdate(cmo,sodisp,work,work2,nocc,nvir,itsoscf,maxsoscf, &
&                        nao,nmo,sodispmax,idis,nproc,myrank,mpi_comm)
!-----------------------------------------------------------------------------
!
! Update molecular orbitals using approximated SOSCF method
!
      implicit none
      integer,intent(in) :: nocc, nvir, itsoscf, maxsoscf, nao, nmo, nproc, myrank, mpi_comm
      integer,intent(in) :: idis(nproc,14)
      integer :: imo, jmo, num, istart
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, p01= 1.0D-01
      real(8),intent(in) :: sodisp(nocc,nvir,maxsoscf), sodispmax
      real(8),intent(inout) :: cmo(nao,nao)
      real(8),intent(out) :: work(nmo,nmo), work2(*)
      real(8) :: factor
!
! U = exp(A) = I + A
!
!ishimura
      work(:,:)= zero
!ishimura
!      call dgemm('N','T',nocc,nocc,nvir,-half,sodisp(1,1,itsoscf),nocc,sodisp(1,1,itsoscf), &
!&                nocc,zero,work,nmo)
!      call dgemm('T','N',nvir,nvir,nocc,-half,sodisp(1,1,itsoscf),nocc,sodisp(1,1,itsoscf), &
!&                nmo,zero,work(nocc+1,nocc+1),nmo)
      if(sodispmax < p01) then
!$OMP parallel do
        do imo= 1,nvir
          do jmo= 1,nocc
            work(jmo,imo+nocc)= sodisp(jmo,imo,itsoscf)
            work(imo+nocc,jmo)=-sodisp(jmo,imo,itsoscf)
          enddo
        enddo
!$OMP end parallel do
      else
        factor= p01/sodispmax
!$OMP parallel do
        do imo= 1,nvir
          do jmo= 1,nocc
            work(jmo,imo+nocc)= sodisp(jmo,imo,itsoscf)*factor
            work(imo+nocc,jmo)=-sodisp(jmo,imo,itsoscf)*factor
          enddo
        enddo
!$OMP end parallel do
      endif
!
      do imo= 1,nmo
        work(imo,imo)= work(imo,imo)+one
      enddo
!
      call orthonorm(work,nmo,nproc,myrank,mpi_comm)
!
      num= idis(myrank+1,1)
      istart= idis(myrank+1,2)+1
      if(num > 0) call dgemm('N','N',nao,num,nmo,one,cmo,nao,work(1,istart),nmo,zero, &
&                            work2,nao)
      call para_allgathervr(work2,idis(myrank+1,3),cmo,idis(1,3),idis(1,4),nproc,mpi_comm)
      return
end


!----------------------------------------------------------------------------------------
  subroutine calcspin(sz,s2,dmtrxa,dmtrxb,overlap,work,work2,work3,neleca,nelecb,nao, &
&                     idis,nproc,myrank,mpi_comm)
!----------------------------------------------------------------------------------------
!
! Calculate spin expectation values, sz and S^2
!
      implicit none
      integer,intent(in) :: neleca, nelecb, nao, nproc, myrank, mpi_comm, idis(nproc,14)
      integer :: num, istart, ij, i, j
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, two=2.0D+00
      real(8),intent(in) :: dmtrxa(nao*(nao+1)/2), dmtrxb(nao*(nao+1)/2), overlap(nao*(nao+1)/2)
      real(8),intent(out) :: sz, s2, work(nao,nao), work2(nao,nao), work3(*)
!
      call expand(dmtrxa,work,nao)
      call expand2(overlap,work2,nao)
!
! Calcualte S*Da*S
!
      num=idis(myrank+1,1)
      istart=idis(myrank+1,2)+1
      if(num > 0) then
        call dsymm('L','U',nao,num,one,work,nao,work2(1,istart),nao,zero,work3,nao)
        call dgemm('N','N',nao,num,nao,one,work2,nao,work3,nao,zero,work,nao)
      endif
      call para_allgathervr(work,num*nao,work2,idis(1,3),idis(1,4),nproc,mpi_comm)
!
! Trace(Db*S*Da*S)
!
      s2= zero
!$OMP parallel do private(ij) reduction(+:s2)
      do i= 1,nao
        ij= i*(i-1)/2
        do j= 1,i-1
          ij= ij+1
           s2= s2+dmtrxb(ij)*work2(j,i)
        enddo
        ij= ij+1
        s2= s2+dmtrxb(ij)*work2(i,i)*half
      enddo
!$OMP end parallel do
      s2= s2*two
!
      sz= half*(neleca-nelecb)
      s2= sz*sz+half*(neleca+nelecb)-s2
!
      return
end


!------------------------------------------
  subroutine calcrmulliken(dmtrx,overlap)
!------------------------------------------
!
! Execute Mulliken population analysis for closed-shell
!
      use modparallel, only : master
      use modmolecule, only : natom, znuc, numatomic
      use modbasis, only : locatom, nao, nshell, locbf, mbf
      implicit none
      integer :: ii, jj, ij, ish, iatom, locbfi
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: dmtrx(nao*(nao+1)/2), overlap(nao*(nao+1)/2)
      real(8) :: grossorb(nao), grossatom(natom), totalgross
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
!
      totalgross= zero
      grossorb(:)= zero
      grossatom(:)= zero
!
! Calculate Gross orbital population
!
!$OMP parallel do schedule(static,1) private(ij) reduction(+:grossorb)
      do ii= 1,nao
        ij= ii*(ii-1)/2
        do jj= 1,ii-1
          ij= ij+1
          grossorb(ii)= grossorb(ii)+dmtrx(ij)*overlap(ij)
          grossorb(jj)= grossorb(jj)+dmtrx(ij)*overlap(ij)
        enddo
        ij= ij+1
        grossorb(ii)= grossorb(ii)+dmtrx(ij)*overlap(ij)
      enddo
!$OMP end parallel do
!
! Calculate Gross atom population
!
      do ish= 1,nshell
        iatom= locatom(ish)
        locbfi= locbf(ish)
        do ii= 1,mbf(ish)
          grossatom(iatom)= grossatom(iatom)+grossorb(locbfi+ii)
        enddo
      enddo
      do iatom= 1,natom
        totalgross= totalgross+znuc(iatom)-grossatom(iatom)
      enddo
!
      if(master) then
        write(*,'(" -------------------------------------")')
        write(*,'("      Mulliken Population Analysis")')
        write(*,'("     Atom     Population     Charge")')
        write(*,'(" -------------------------------------")')
        do iatom= 1,natom
          write(*,'(1x,i4,2x,a3,2f13.6)')iatom,table(numatomic(iatom)), &
&                                        grossatom(iatom),znuc(iatom)-grossatom(iatom)
        enddo
        write(*,'(" -------------------------------------")')
        write(*,'("     Total",13x,f13.6)')totalgross
        write(*,'(" -------------------------------------")')
        write(*,*)
      endif
!
      return
end


!--------------------------------------------------
  subroutine calcumulliken(dmtrxa,dmtrxb,overlap)
!--------------------------------------------------
!
! Execute Mulliken population analysis for open-shell
!
      use modparallel, only : master
      use modmolecule, only : natom, znuc, numatomic
      use modbasis, only : locatom, nao, nshell, locbf, mbf
      implicit none
      integer :: ii, jj, ij, ish, iatom, locbfi
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: dmtrxa(nao*(nao+1)/2), dmtrxb(nao*(nao+1)/2), overlap(nao*(nao+1)/2)
      real(8) :: grossorb(nao), grossatom(natom), totalgross
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
!
      totalgross= zero
      grossorb(:)= zero
      grossatom(:)= zero
!
! Calculate Gross orbital population
!
!$OMP parallel do schedule(static,1) private(ij) reduction(+:grossorb)
      do ii= 1,nao
        ij= ii*(ii-1)/2
        do jj= 1,ii-1
          ij= ij+1
          grossorb(ii)= grossorb(ii)+(dmtrxa(ij)+dmtrxb(ij))*overlap(ij)
          grossorb(jj)= grossorb(jj)+(dmtrxa(ij)+dmtrxb(ij))*overlap(ij)
        enddo
        ij= ij+1
        grossorb(ii)= grossorb(ii)+(dmtrxa(ij)+dmtrxb(ij))*overlap(ij)
      enddo
!$OMP end parallel do
!
! Calculate Gross atom population
!
      do ish= 1,nshell
        iatom= locatom(ish)
        locbfi= locbf(ish)
        do ii= 1,mbf(ish)
          grossatom(iatom)= grossatom(iatom)+grossorb(locbfi+ii)
        enddo
      enddo
      do iatom= 1,natom
        totalgross= totalgross+znuc(iatom)-grossatom(iatom)
      enddo
!
      if(master) then
        write(*,'(" -------------------------------------")')
        write(*,'("      Mulliken Population Analysis")')
        write(*,'("     Atom     Population     Charge")')
        write(*,'(" -------------------------------------")')
        do iatom= 1,natom
          write(*,'(1x,i4,2x,a3,2f13.6)')iatom,table(numatomic(iatom)), &
&                                        grossatom(iatom),znuc(iatom)-grossatom(iatom)
        enddo
        write(*,'(" -------------------------------------")')
        write(*,'("     Total",13x,f13.6)')totalgross
        write(*,'(" -------------------------------------")')
        write(*,*)
      endif
!
      return
end


!------------------------------------------------------------------
  subroutine calcrdipole(dipmat,work,dmtrx,nproc,myrank,mpi_comm)
!------------------------------------------------------------------
!
! Driver of dipole moment calculation for closed-shell
!
      use modparallel, only : master
      use modbasis, only : nao
      use modunit, only : todebye
      use modmolecule, only : natom, coord, znuc
      implicit none
      integer,intent(in) :: nproc, myrank, mpi_comm
      integer :: iatom
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: dmtrx((nao*(nao+1))/2)
      real(8),intent(out) :: dipmat((nao*(nao+1))/2,3), work((nao*(nao+1))/2,3)
      real(8) :: dipcenter(3), xdip, ydip, zdip, totaldip, tridot
      real(8) :: xdipplus, ydipplus, zdipplus, xdipminus, ydipminus, zdipminus
!
! Nuclear part
!
      xdipplus= zero
      ydipplus= zero
      zdipplus= zero
!
      do iatom= 1,natom
        xdipplus= xdipplus+coord(1,iatom)*znuc(iatom)
        ydipplus= ydipplus+coord(2,iatom)*znuc(iatom)
        zdipplus= zdipplus+coord(3,iatom)*znuc(iatom)
      enddo
!
! Electron part
!
      dipcenter(:)= zero
!
      call calcmatdipole(dipmat,work,dipcenter,nproc,myrank,mpi_comm)
!
      xdipminus=-tridot(dmtrx,dipmat(1,1),nao)
      ydipminus=-tridot(dmtrx,dipmat(1,2),nao)
      zdipminus=-tridot(dmtrx,dipmat(1,3),nao)
!
! Sum Nuclear and Electron parts
!
      xdip=(xdipplus+xdipminus)*todebye
      ydip=(ydipplus+ydipminus)*todebye
      zdip=(zdipplus+zdipminus)*todebye
      totaldip= sqrt(xdip*xdip+ydip*ydip+zdip*zdip)
!
      if(master) then
        write(*,'("  ----------------------------------------------")')
        write(*,'("                Dipole Momemt (Debye)")')
        write(*,'("         X          Y          Z       Total")')
        write(*,'("  ----------------------------------------------")')
        write(*,'(2x,4f11.4)')xdip, ydip, zdip, totaldip
        write(*,'("  ----------------------------------------------",/)')
      endif
!
      return
end


!-----------------------------------------------------------------
  subroutine calcroctupole(dipmat,quadpmat,octpmat,work,dmtrx, &
&                          nproc,myrank,mpi_comm)
!-----------------------------------------------------------------
!
! Driver of dipole, quadrupole, and octupole moment calculation for closed-shell
!
      use modparallel, only : master
      use modbasis, only : nao
      use modunit, only : todebye, toang
      use modmolecule, only : natom, coord, znuc
      implicit none
      integer,intent(in) :: nproc, myrank, mpi_comm
      integer :: iatom, ii
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, two=2.0D+00, three=3.0D+00
      real(8),parameter :: four=4.0D+00, five=5.0D+00
      real(8),intent(in) :: dmtrx((nao*(nao+1))/2)
      real(8),intent(out) :: dipmat((nao*(nao+1))/2,3), quadpmat((nao*(nao+1))/2,6)
      real(8),intent(out) :: octpmat((nao*(nao+1))/2,10), work((nao*(nao+1))/2*10)
      real(8) :: dipcenter(3), xdip, ydip, zdip, totaldip, tridot
      real(8) :: xdipplus, ydipplus, zdipplus, xdipminus, ydipminus, zdipminus
      real(8) :: xxquadpplus, xyquadpplus, xzquadpplus, yyquadpplus, yzquadpplus, zzquadpplus
      real(8) :: xxquadpminus, xyquadpminus, xzquadpminus, yyquadpminus, yzquadpminus, zzquadpminus
      real(8) :: xxquadp, xyquadp, xzquadp, yyquadp, yzquadp, zzquadp, quadp(6), facp
      real(8) :: xxxoctpplus, xxyoctpplus, xxzoctpplus, xyyoctpplus, xyzoctpplus, xzzoctpplus
      real(8) :: yyyoctpplus, yyzoctpplus, yzzoctpplus, zzzoctpplus
      real(8) :: xxxoctpminus, xxyoctpminus, xxzoctpminus, xyyoctpminus, xyzoctpminus, xzzoctpminus
      real(8) :: yyyoctpminus, yyzoctpminus, yzzoctpminus, zzzoctpminus
      real(8) :: xxxoctp, xxyoctp, xxzoctp, xyyoctp, xyzoctp, xzzoctp
      real(8) :: yyyoctp, yyzoctp, yzzoctp, zzzoctp, octp(10), faco
!
! Nuclear part
!
      xdipplus= zero
      ydipplus= zero
      zdipplus= zero
      xxquadpplus= zero
      xyquadpplus= zero
      xzquadpplus= zero
      yyquadpplus= zero
      yzquadpplus= zero
      zzquadpplus= zero
      xxxoctpplus= zero
      xxyoctpplus= zero
      xxzoctpplus= zero
      xyyoctpplus= zero
      xyzoctpplus= zero
      xzzoctpplus= zero
      yyyoctpplus= zero
      yyzoctpplus= zero
      yzzoctpplus= zero
      zzzoctpplus= zero
!
      do iatom= 1,natom
        xdipplus= xdipplus+coord(1,iatom)*znuc(iatom)
        ydipplus= ydipplus+coord(2,iatom)*znuc(iatom)
        zdipplus= zdipplus+coord(3,iatom)*znuc(iatom)
        xxquadpplus= xxquadpplus+coord(1,iatom)*coord(1,iatom)*znuc(iatom)
        xyquadpplus= xyquadpplus+coord(1,iatom)*coord(2,iatom)*znuc(iatom)
        xzquadpplus= xzquadpplus+coord(1,iatom)*coord(3,iatom)*znuc(iatom)
        yyquadpplus= yyquadpplus+coord(2,iatom)*coord(2,iatom)*znuc(iatom)
        yzquadpplus= yzquadpplus+coord(2,iatom)*coord(3,iatom)*znuc(iatom)
        zzquadpplus= zzquadpplus+coord(3,iatom)*coord(3,iatom)*znuc(iatom)
        xxxoctpplus= xxxoctpplus+coord(1,iatom)*coord(1,iatom)*coord(1,iatom)*znuc(iatom)
        xxyoctpplus= xxyoctpplus+coord(1,iatom)*coord(1,iatom)*coord(2,iatom)*znuc(iatom)
        xxzoctpplus= xxzoctpplus+coord(1,iatom)*coord(1,iatom)*coord(3,iatom)*znuc(iatom)
        xyyoctpplus= xyyoctpplus+coord(1,iatom)*coord(2,iatom)*coord(2,iatom)*znuc(iatom)
        xyzoctpplus= xyzoctpplus+coord(1,iatom)*coord(2,iatom)*coord(3,iatom)*znuc(iatom)
        xzzoctpplus= xzzoctpplus+coord(1,iatom)*coord(3,iatom)*coord(3,iatom)*znuc(iatom)
        yyyoctpplus= yyyoctpplus+coord(2,iatom)*coord(2,iatom)*coord(2,iatom)*znuc(iatom)
        yyzoctpplus= yyzoctpplus+coord(2,iatom)*coord(2,iatom)*coord(3,iatom)*znuc(iatom)
        yzzoctpplus= yzzoctpplus+coord(2,iatom)*coord(3,iatom)*coord(3,iatom)*znuc(iatom)
        zzzoctpplus= zzzoctpplus+coord(3,iatom)*coord(3,iatom)*coord(3,iatom)*znuc(iatom)
      enddo
!
! Electron part
!
      dipcenter(:)= zero
!
      call calcmatoctupole(dipmat,quadpmat,octpmat,work,dipcenter,nproc,myrank,mpi_comm)
!
      xdipminus=-tridot(dmtrx,dipmat(1,1),nao)
      ydipminus=-tridot(dmtrx,dipmat(1,2),nao)
      zdipminus=-tridot(dmtrx,dipmat(1,3),nao)
      xxquadpminus=-tridot(dmtrx,quadpmat(1,1),nao)
      xyquadpminus=-tridot(dmtrx,quadpmat(1,2),nao)
      xzquadpminus=-tridot(dmtrx,quadpmat(1,3),nao)
      yyquadpminus=-tridot(dmtrx,quadpmat(1,4),nao)
      yzquadpminus=-tridot(dmtrx,quadpmat(1,5),nao)
      zzquadpminus=-tridot(dmtrx,quadpmat(1,6),nao)
      xxxoctpminus=-tridot(dmtrx,octpmat(1, 1),nao)
      xxyoctpminus=-tridot(dmtrx,octpmat(1, 2),nao)
      xxzoctpminus=-tridot(dmtrx,octpmat(1, 3),nao)
      xyyoctpminus=-tridot(dmtrx,octpmat(1, 4),nao)
      xyzoctpminus=-tridot(dmtrx,octpmat(1, 5),nao)
      xzzoctpminus=-tridot(dmtrx,octpmat(1, 6),nao)
      yyyoctpminus=-tridot(dmtrx,octpmat(1, 7),nao)
      yyzoctpminus=-tridot(dmtrx,octpmat(1, 8),nao)
      yzzoctpminus=-tridot(dmtrx,octpmat(1, 9),nao)
      zzzoctpminus=-tridot(dmtrx,octpmat(1,10),nao)
!
! Sum Nuclear and Electron parts
!
      facp= todebye*toang*half
      faco= todebye*toang*toang*half
      xdip=(xdipplus+xdipminus)*todebye
      ydip=(ydipplus+ydipminus)*todebye
      zdip=(zdipplus+zdipminus)*todebye
      xxquadp=(xxquadpplus+xxquadpminus)*facp
      xyquadp=(xyquadpplus+xyquadpminus)*facp
      xzquadp=(xzquadpplus+xzquadpminus)*facp
      yyquadp=(yyquadpplus+yyquadpminus)*facp
      yzquadp=(yzquadpplus+yzquadpminus)*facp
      zzquadp=(zzquadpplus+zzquadpminus)*facp
      xxxoctp=(xxxoctpplus+xxxoctpminus)*faco
      xxyoctp=(xxyoctpplus+xxyoctpminus)*faco
      xxzoctp=(xxzoctpplus+xxzoctpminus)*faco
      xyyoctp=(xyyoctpplus+xyyoctpminus)*faco
      xyzoctp=(xyzoctpplus+xyzoctpminus)*faco
      xzzoctp=(xzzoctpplus+xzzoctpminus)*faco
      yyyoctp=(yyyoctpplus+yyyoctpminus)*faco
      yyzoctp=(yyzoctpplus+yyzoctpminus)*faco
      yzzoctp=(yzzoctpplus+yzzoctpminus)*faco
      zzzoctp=(zzzoctpplus+zzzoctpminus)*faco
      totaldip= sqrt(xdip*xdip+ydip*ydip+zdip*zdip)
!
      quadp(1)= xxquadp*two-yyquadp    -zzquadp
      quadp(2)= xyquadp*three
      quadp(3)= xzquadp*three
      quadp(4)=-xxquadp    +yyquadp*two-zzquadp
      quadp(5)= yzquadp*three
      quadp(6)=-xxquadp    -yyquadp    +zzquadp*two
      octp( 1)= xxxoctp*two  -xyyoctp*three-xzzoctp*three
      octp( 2)= xxyoctp*four -yyyoctp      -yzzoctp
      octp( 3)= xxzoctp*four -yyzoctp      -zzzoctp
      octp( 4)=-xxxoctp      +xyyoctp*four -xzzoctp
      octp( 5)= xyzoctp*five
      octp( 6)=-xxxoctp      -xyyoctp      +xzzoctp*four
      octp( 7)=-xxyoctp*three+yyyoctp*two  -yzzoctp*three
      octp( 8)=-xxzoctp      +yyzoctp*four -zzzoctp
      octp( 9)=-xxyoctp      -yyyoctp      +yzzoctp*four
      octp(10)=-xxzoctp*three-yyzoctp*three+zzzoctp*two
!
      if(master) then
        write(*,'("  ----------------------------------------------")')
        write(*,'("                Dipole Momemt (Debye)")')
        write(*,'("         X          Y          Z       Total")')
        write(*,'("  ----------------------------------------------")')
        write(*,'(2x,4f11.4)')xdip, ydip, zdip, totaldip
        write(*,'("  ----------------------------------------------",/)')
        write(*,'("  --------------------------------------------------------------------")')
        write(*,'("                Quadrupole Momemt (Debye*Angstrom)")')
        write(*,'("         XX         XY         XZ         YY         YZ         ZZ")')
        write(*,'("  --------------------------------------------------------------------")')
        write(*,'(2x,6f11.4)')(quadp(ii),ii=1,6)
        write(*,'("  --------------------------------------------------------------------",/)')
        write(*,'("  ---------------------------------------------------------")')
        write(*,'("                Octupole Momemt (Debye*Angstrom^2)")')
        write(*,'("        XXX        XXY        XXZ       XYY        XYZ")')
        write(*,'("  ---------------------------------------------------------")')
        write(*,'(2x,5f11.4)')(octp(ii),ii=1,5)
        write(*,'("  ---------------------------------------------------------")')
        write(*,'("        XZZ        YYY        YYZ       YZZ        ZZZ")')
        write(*,'("  ---------------------------------------------------------")')
        write(*,'(2x,5f11.4)')(octp(ii),ii=6,10)
        write(*,'("  ---------------------------------------------------------",/)')
      endif
!
      return
end


!--------------------------------------------------------------------------
  subroutine calcudipole(dipmat,work,dmtrxa,dmtrxb,nproc,myrank,mpi_comm)
!--------------------------------------------------------------------------
!
! Driver of dipole moment calculation for open-shell
!
      use modparallel, only : master
      use modbasis, only : nao
      use modunit, only : todebye
      use modmolecule, only : natom, coord, znuc
      implicit none
      integer,intent(in) :: nproc, myrank, mpi_comm
      integer :: iatom
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: dmtrxa((nao*(nao+1))/2), dmtrxb((nao*(nao+1))/2)
      real(8),intent(out) :: dipmat((nao*(nao+1))/2,3), work((nao*(nao+1))/2,3)
      real(8) :: dipcenter(3), xdip, ydip, zdip, totaldip, tridot
      real(8) :: xdipplus, ydipplus, zdipplus, xdipminus, ydipminus, zdipminus
!
! Nuclear part
!
      xdipplus= zero
      ydipplus= zero
      zdipplus= zero
!
      do iatom= 1,natom
        xdipplus= xdipplus+coord(1,iatom)*znuc(iatom)
        ydipplus= ydipplus+coord(2,iatom)*znuc(iatom)
        zdipplus= zdipplus+coord(3,iatom)*znuc(iatom)
      enddo
!
! Electron part
!
      dipcenter(:)= zero
!
      call calcmatdipole(dipmat,work,dipcenter,nproc,myrank,mpi_comm)
!
      xdipminus=-tridot(dmtrxa,dipmat(1,1),nao)-tridot(dmtrxb,dipmat(1,1),nao)
      ydipminus=-tridot(dmtrxa,dipmat(1,2),nao)-tridot(dmtrxb,dipmat(1,2),nao)
      zdipminus=-tridot(dmtrxa,dipmat(1,3),nao)-tridot(dmtrxb,dipmat(1,3),nao)
!
! Sum Nuclear and Electron parts
!
      xdip=(xdipplus+xdipminus)*todebye
      ydip=(ydipplus+ydipminus)*todebye
      zdip=(zdipplus+zdipminus)*todebye
      totaldip= sqrt(xdip*xdip+ydip*ydip+zdip*zdip)
!
      if(master) then
        write(*,'("  ----------------------------------------------")')
        write(*,'("                Dipole Momemt (Debye)")')
        write(*,'("         X          Y          Z       Total")')
        write(*,'("  ----------------------------------------------")')
        write(*,'(2x,4f11.4)')xdip, ydip, zdip, totaldip
        write(*,'("  ----------------------------------------------",/)')
      endif
!
      return
end


!-------------------------------------------------------------------------
  subroutine calcuoctupole(dipmat,quadpmat,octpmat,work,dmtrxa,dmtrxb, &
&                          nproc,myrank,mpi_comm)
!-------------------------------------------------------------------------
!
! Driver of dipole, quadrupole, and octupole moment calculation for open-shell
!
      use modparallel, only : master
      use modbasis, only : nao
      use modunit, only : todebye, toang
      use modmolecule, only : natom, coord, znuc
      implicit none
      integer,intent(in) :: nproc, myrank, mpi_comm
      integer :: iatom, ii
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, two=2.0D+00, three=3.0D+00
      real(8),parameter :: four=4.0D+00, five=5.0D+00
      real(8),intent(in) :: dmtrxa((nao*(nao+1))/2), dmtrxb((nao*(nao+1))/2)
      real(8),intent(out) :: dipmat((nao*(nao+1))/2,3), quadpmat((nao*(nao+1))/2,6)
      real(8),intent(out) :: octpmat((nao*(nao+1))/2,10), work((nao*(nao+1))/2*10)
      real(8) :: dipcenter(3), xdip, ydip, zdip, totaldip, tridot
      real(8) :: xdipplus, ydipplus, zdipplus, xdipminus, ydipminus, zdipminus
      real(8) :: xxquadpplus, xyquadpplus, xzquadpplus, yyquadpplus, yzquadpplus, zzquadpplus
      real(8) :: xxquadpminus, xyquadpminus, xzquadpminus, yyquadpminus, yzquadpminus, zzquadpminus
      real(8) :: xxquadp, xyquadp, xzquadp, yyquadp, yzquadp, zzquadp, quadp(6), facp
      real(8) :: xxxoctpplus, xxyoctpplus, xxzoctpplus, xyyoctpplus, xyzoctpplus, xzzoctpplus
      real(8) :: yyyoctpplus, yyzoctpplus, yzzoctpplus, zzzoctpplus
      real(8) :: xxxoctpminus, xxyoctpminus, xxzoctpminus, xyyoctpminus, xyzoctpminus, xzzoctpminus
      real(8) :: yyyoctpminus, yyzoctpminus, yzzoctpminus, zzzoctpminus
      real(8) :: xxxoctp, xxyoctp, xxzoctp, xyyoctp, xyzoctp, xzzoctp
      real(8) :: yyyoctp, yyzoctp, yzzoctp, zzzoctp, octp(10), faco
!
! Nuclear part
!
      xdipplus= zero
      ydipplus= zero
      zdipplus= zero
      xxquadpplus= zero
      xyquadpplus= zero
      xzquadpplus= zero
      yyquadpplus= zero
      yzquadpplus= zero
      zzquadpplus= zero
      xxxoctpplus= zero
      xxyoctpplus= zero
      xxzoctpplus= zero
      xyyoctpplus= zero
      xyzoctpplus= zero
      xzzoctpplus= zero
      yyyoctpplus= zero
      yyzoctpplus= zero
      yzzoctpplus= zero
      zzzoctpplus= zero
!
      do iatom= 1,natom
        xdipplus= xdipplus+coord(1,iatom)*znuc(iatom)
        ydipplus= ydipplus+coord(2,iatom)*znuc(iatom)
        zdipplus= zdipplus+coord(3,iatom)*znuc(iatom)
        xxquadpplus= xxquadpplus+coord(1,iatom)*coord(1,iatom)*znuc(iatom)
        xyquadpplus= xyquadpplus+coord(1,iatom)*coord(2,iatom)*znuc(iatom)
        xzquadpplus= xzquadpplus+coord(1,iatom)*coord(3,iatom)*znuc(iatom)
        yyquadpplus= yyquadpplus+coord(2,iatom)*coord(2,iatom)*znuc(iatom)
        yzquadpplus= yzquadpplus+coord(2,iatom)*coord(3,iatom)*znuc(iatom)
        zzquadpplus= zzquadpplus+coord(3,iatom)*coord(3,iatom)*znuc(iatom)
        xxxoctpplus= xxxoctpplus+coord(1,iatom)*coord(1,iatom)*coord(1,iatom)*znuc(iatom)
        xxyoctpplus= xxyoctpplus+coord(1,iatom)*coord(1,iatom)*coord(2,iatom)*znuc(iatom)
        xxzoctpplus= xxzoctpplus+coord(1,iatom)*coord(1,iatom)*coord(3,iatom)*znuc(iatom)
        xyyoctpplus= xyyoctpplus+coord(1,iatom)*coord(2,iatom)*coord(2,iatom)*znuc(iatom)
        xyzoctpplus= xyzoctpplus+coord(1,iatom)*coord(2,iatom)*coord(3,iatom)*znuc(iatom)
        xzzoctpplus= xzzoctpplus+coord(1,iatom)*coord(3,iatom)*coord(3,iatom)*znuc(iatom)
        yyyoctpplus= yyyoctpplus+coord(2,iatom)*coord(2,iatom)*coord(2,iatom)*znuc(iatom)
        yyzoctpplus= yyzoctpplus+coord(2,iatom)*coord(2,iatom)*coord(3,iatom)*znuc(iatom)
        yzzoctpplus= yzzoctpplus+coord(2,iatom)*coord(3,iatom)*coord(3,iatom)*znuc(iatom)
        zzzoctpplus= zzzoctpplus+coord(3,iatom)*coord(3,iatom)*coord(3,iatom)*znuc(iatom)
      enddo
!
! Electron part
!
      dipcenter(:)= zero
!
      call calcmatoctupole(dipmat,quadpmat,octpmat,work,dipcenter,nproc,myrank,mpi_comm)
!
      do ii= 1,nao*(nao+1)/2
        work(ii)= dmtrxa(ii)+dmtrxb(ii)
      enddo
!
      xdipminus=-tridot(work,dipmat(1,1),nao)
      ydipminus=-tridot(work,dipmat(1,2),nao)
      zdipminus=-tridot(work,dipmat(1,3),nao)
      xxquadpminus=-tridot(work,quadpmat(1,1),nao)
      xyquadpminus=-tridot(work,quadpmat(1,2),nao)
      xzquadpminus=-tridot(work,quadpmat(1,3),nao)
      yyquadpminus=-tridot(work,quadpmat(1,4),nao)
      yzquadpminus=-tridot(work,quadpmat(1,5),nao)
      zzquadpminus=-tridot(work,quadpmat(1,6),nao)
      xxxoctpminus=-tridot(work,octpmat(1, 1),nao)
      xxyoctpminus=-tridot(work,octpmat(1, 2),nao)
      xxzoctpminus=-tridot(work,octpmat(1, 3),nao)
      xyyoctpminus=-tridot(work,octpmat(1, 4),nao)
      xyzoctpminus=-tridot(work,octpmat(1, 5),nao)
      xzzoctpminus=-tridot(work,octpmat(1, 6),nao)
      yyyoctpminus=-tridot(work,octpmat(1, 7),nao)
      yyzoctpminus=-tridot(work,octpmat(1, 8),nao)
      yzzoctpminus=-tridot(work,octpmat(1, 9),nao)
      zzzoctpminus=-tridot(work,octpmat(1,10),nao)
!
! Sum Nuclear and Electron parts
!
      facp= todebye*toang*half
      faco= todebye*toang*toang*half
      xdip=(xdipplus+xdipminus)*todebye
      ydip=(ydipplus+ydipminus)*todebye
      zdip=(zdipplus+zdipminus)*todebye
      xxquadp=(xxquadpplus+xxquadpminus)*facp
      xyquadp=(xyquadpplus+xyquadpminus)*facp
      xzquadp=(xzquadpplus+xzquadpminus)*facp
      yyquadp=(yyquadpplus+yyquadpminus)*facp
      yzquadp=(yzquadpplus+yzquadpminus)*facp
      zzquadp=(zzquadpplus+zzquadpminus)*facp
      xxxoctp=(xxxoctpplus+xxxoctpminus)*faco
      xxyoctp=(xxyoctpplus+xxyoctpminus)*faco
      xxzoctp=(xxzoctpplus+xxzoctpminus)*faco
      xyyoctp=(xyyoctpplus+xyyoctpminus)*faco
      xyzoctp=(xyzoctpplus+xyzoctpminus)*faco
      xzzoctp=(xzzoctpplus+xzzoctpminus)*faco
      yyyoctp=(yyyoctpplus+yyyoctpminus)*faco
      yyzoctp=(yyzoctpplus+yyzoctpminus)*faco
      yzzoctp=(yzzoctpplus+yzzoctpminus)*faco
      zzzoctp=(zzzoctpplus+zzzoctpminus)*faco
      totaldip= sqrt(xdip*xdip+ydip*ydip+zdip*zdip)
!
      quadp(1)= xxquadp*two-yyquadp    -zzquadp
      quadp(2)= xyquadp*three
      quadp(3)= xzquadp*three
      quadp(4)=-xxquadp    +yyquadp*two-zzquadp
      quadp(5)= yzquadp*three
      quadp(6)=-xxquadp    -yyquadp    +zzquadp*two
      octp( 1)= xxxoctp*two  -xyyoctp*three-xzzoctp*three
      octp( 2)= xxyoctp*four -yyyoctp      -yzzoctp
      octp( 3)= xxzoctp*four -yyzoctp      -zzzoctp
      octp( 4)=-xxxoctp      +xyyoctp*four -xzzoctp
      octp( 5)= xyzoctp*five
      octp( 6)=-xxxoctp      -xyyoctp      +xzzoctp*four
      octp( 7)=-xxyoctp*three+yyyoctp*two  -yzzoctp*three
      octp( 8)=-xxzoctp      +yyzoctp*four -zzzoctp
      octp( 9)=-xxyoctp      -yyyoctp      +yzzoctp*four
      octp(10)=-xxzoctp*three-yyzoctp*three+zzzoctp*two
!
      if(master) then
        write(*,'("  ----------------------------------------------")')
        write(*,'("                Dipole Momemt (Debye)")')
        write(*,'("         X          Y          Z       Total")')
        write(*,'("  ----------------------------------------------")')
        write(*,'(2x,4f11.4)')xdip, ydip, zdip, totaldip
        write(*,'("  ----------------------------------------------",/)')
        write(*,'("  --------------------------------------------------------------------")')
        write(*,'("                Quadrupole Momemt (Debye*Angstrom)")')
        write(*,'("         XX         XY         XZ         YY         YZ         ZZ")')
        write(*,'("  --------------------------------------------------------------------")')
        write(*,'(2x,6f11.4)')(quadp(ii),ii=1,6)
        write(*,'("  --------------------------------------------------------------------",/)')
        write(*,'("  ---------------------------------------------------------")')
        write(*,'("                Octupole Momemt (Debye*Angstrom^2)")')
        write(*,'("        XXX        XXY        XXZ       XYY        XYZ")')
        write(*,'("  ---------------------------------------------------------")')
        write(*,'(2x,5f11.4)')(octp(ii),ii=1,5)
        write(*,'("  ---------------------------------------------------------")')
        write(*,'("        XZZ        YYY        YYZ       YZZ        ZZZ")')
        write(*,'("  ---------------------------------------------------------")')
        write(*,'(2x,5f11.4)')(octp(ii),ii=6,10)
        write(*,'("  ---------------------------------------------------------",/)')
      endif
!
      return
end


!------------------------------------------------------------------------------
  subroutine calcqcrmn(qcrmn,qcvec,cmo,work,nao,nocc,nvir,itdav,maxqcdiagsub)
!------------------------------------------------------------------------------
!
! Calculate Rmn for quadratically convergent method
!
      implicit none
      integer,intent(in) :: nao, nocc, nvir, itdav, maxqcdiagsub
      integer :: ii, jj, ij
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: qcvec(nocc*nvir+1,maxqcdiagsub+1,2), cmo(nao,nao)
      real(8),intent(out) :: qcrmn(nao*nao), work(nao,nao)
!
      call dgemm('N','N',nao,nvir,nocc,one,cmo,nao,qcvec(2,itdav,1),nocc,zero,qcrmn,nao)
      call dgemm('N','T',nao,nao,nvir,one,qcrmn,nao,cmo(1,nocc+1),nao,zero,work,nao)
!
      ij= 0
!$OMP parallel do private(ij)
      do ii= 1,nao
        ij= ii*(ii-1)/2
        do jj= 1,ii
          ij= ij+1
          qcrmn(ij)= work(jj,ii)+work(ii,jj)
        enddo
      enddo
!$OMP end parallel do
!
      return
end


!------------------------------------------------------------------------------------------
  subroutine calcqcurmn(qcrmna,qcrmnb,qcvec,cmoa,cmob,work,nao,nocca,noccb,nvira,nvirb, &
&                       itdav,maxqcdiagsub)
!------------------------------------------------------------------------------------------
!
! Calculate Rmn for quadratically convergent method of UHF
!
      implicit none
      integer,intent(in) :: nao, nocca, noccb, nvira, nvirb, itdav, maxqcdiagsub
      integer :: ii, jj, ij
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: qcvec(nocca*nvira+noccb*nvirb+1,maxqcdiagsub+1,2)
      real(8),intent(in) :: cmoa(nao,nao), cmob(nao,nao)
      real(8),intent(out) :: qcrmna(nao*nao), qcrmnb(nao*nao), work(nao,nao)
!
      call dgemm('N','N',nao,nvira,nocca,one,cmoa,nao,qcvec(2,itdav,1),nocca,zero,qcrmna,nao)
      call dgemm('N','T',nao,nao,nvira,one,qcrmna,nao,cmoa(1,nocca+1),nao,zero,work,nao)
!
      ij= 0
!$OMP parallel do private(ij)
      do ii= 1,nao
        ij= ii*(ii-1)/2
        do jj= 1,ii
          qcrmna(ij+jj)= work(jj,ii)+work(ii,jj)
        enddo
      enddo
!$OMP end parallel do
!
      call dgemm('N','N',nao,nvirb,noccb,one,cmob,nao,qcvec(nocca*nvira+2,itdav,1),noccb,zero, &
&                qcrmnb,nao)
      call dgemm('N','T',nao,nao,nvirb,one,qcrmnb,nao,cmob(1,noccb+1),nao,zero,work,nao)
!
      ij= 0
!$OMP parallel do private(ij)
      do ii= 1,nao
        ij= ii*(ii-1)/2
        do jj= 1,ii
          qcrmnb(ij+jj)= work(jj,ii)+work(ii,jj)
        enddo
      enddo
!$OMP end parallel do
!
      return
end


!--------------------------------------------------------------------------------------------
  subroutine rhfqc(fock,cmo,qcrmax,qcgmn,qcvec,qcmat,qcmatsave,qceigen,overlap,xint, &
&                  qcwork,work,hfexchange,nao,nmo,nocc,nvir,nshell,maxdim,maxqcdiag, &
&                  maxqcdiagsub,threshqc,nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!--------------------------------------------------------------------------------------------
!
! Driver of Davidson diagonalization for quadratically convergent of RHF
!
      use modparallel, only : master
      implicit none
      integer,intent(in) :: nao, nmo, nocc, nvir, nshell, maxdim, maxqcdiag, maxqcdiagsub
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2, mpi_comm1, mpi_comm2
      integer :: itdav, itqcdiag, ii, ij, jj, kk, ia, ib, istart, icount
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: overlap(nao*(nao+1)/2), xint(nshell*(nshell+1)/2)
      real(8),intent(in) :: hfexchange, threshqc
      real(8),intent(out) :: qcvec(nocc*nvir+1,maxqcdiagsub+1,2), qcrmax(nshell*(nshell+1)/2)
      real(8),intent(out) :: qcwork(nao,nao), qcmat(maxqcdiagsub,maxqcdiagsub)
      real(8),intent(out) :: qcmatsave(maxqcdiagsub*(maxqcdiagsub+1)/2), qceigen(maxqcdiagsub)
      real(8),intent(out) :: work(nao,nao)
      real(8),intent(inout) :: fock(nao,nao), cmo(nao,nao)
      real(8),intent(inout) :: qcgmn(nao*(nao+1)/2)
      real(8) :: tmp, qcnorm, ddot, rotqc
!
      qcmatsave(:)= zero
!
! Calculate Fock matrix in MO basis
!
      call expand(fock,qcwork,nao)
      call dsymm('L','U',nao,nocc+nvir,one,qcwork,nao,cmo,nao,zero,work,nao)
      call dgemm('T','N',nocc+nvir,nocc+nvir,nao,one,cmo,nao,work,nao,zero,fock,nao)
!
! Calculate initial (first and second) qc vector
!
      qcvec(1,1,1)= one
      qcvec(2:nocc*nvir+1,1,1)= zero
      qcvec(1,1,2)= zero
      qcvec(1,2,1)= zero
!
      qcnorm= zero
!$OMP parallel do private(ij) reduction(+:qcnorm)
      do ii= 1,nvir
        ij=(ii-1)*nocc+1
        do jj= 1,nocc
          qcvec(ij+jj,1,2)= fock(jj,ii+nocc)
          qcvec(ij+jj,2,1)= fock(jj,ii+nocc)
          qcnorm= qcnorm+fock(jj,ii+nocc)*fock(jj,ii+nocc)
        enddo
      enddo
!$OMP end parallel do
      qcnorm= one/sqrt(qcnorm)
!$OMP parallel do
      do ii= 2,nocc*nvir+1
        qcvec(ii,2,1)= qcvec(ii,2,1)*qcnorm
      enddo
!$OMP end parallel do
!
! Start Davidson diagonalization
!
      itdav= 2
      do itqcdiag= 2,maxqcdiag
!
! Calculate Gmn
!
        call calcqcrmn(qcwork,qcvec,cmo,work,nao,nocc,nvir,itdav,maxqcdiagsub)
        call calcrdmax(qcwork,qcrmax,work,nproc2,myrank2,mpi_comm2)
        call formrdftfock(qcgmn,work,qcwork,qcrmax,xint,maxdim,hfexchange, &
&                         nproc1,myrank1,mpi_comm1)
!
! Add two-electron integral contribution
!
        call expand(qcgmn,qcwork,nao)
        call dsymm('L','U',nao,nvir,one,qcwork,nao,cmo(1,nocc+1),nao,zero,work,nao)
        call dgemm('T','N',nocc,nvir,nao,one,cmo,nao,work,nao,zero,qcvec(2,itdav,2),nocc)
!
! Add Fock matrix element contribution
!
        tmp= zero
!$OMP parallel do collapse(2) private(kk) reduction(+:tmp)
        do ia= 1,nvir
          do ii= 1,nocc
            kk= (ia-1)*nocc+ii+1
            tmp= tmp+fock(ii,ia+nocc)*qcvec(kk,itdav,1)
            qcvec(kk,itdav,2)= qcvec(kk,itdav,2)+fock(ii,ia+nocc)*qcvec(1,itdav,1)
            do ib= 1,nvir
              qcvec(kk,itdav,2)= qcvec(kk,itdav,2)+fock(ib+nocc,ia+nocc) &
&                                                 *qcvec((ib-1)*nocc+ii+1,itdav,1)
            enddo
            do ij= 1,nocc
              qcvec(kk,itdav,2)= qcvec(kk,itdav,2)-fock(ij,ii)*qcvec((ia-1)*nocc+ij+1,itdav,1)
            enddo
          enddo
        enddo
!$OMP end parallel do
        qcvec(1,itdav,2)= tmp
!
! Calculate small matrix <b|A|b>
!
        qceigen(:)= zero
!$OMP parallel reduction(+:qceigen)
        do ii= 1,itdav
!$OMP do
          do jj= 1,nocc*nvir+1
            qceigen(ii)= qceigen(ii)+qcvec(jj,ii,1)*qcvec(jj,itdav,2)
          enddo
!$OMP end do
        enddo
!$OMP end parallel
!
        istart= itdav*(itdav-1)/2
        do ii= 1,itdav
          qcmatsave(istart+ii)= qceigen(ii)
        enddo
        icount= 0
        do ii= 1,itdav
          do jj= 1,ii
            icount= icount+1
            qcmat(jj,ii)= qcmatsave(icount)
          enddo
        enddo
!
! Diagonalize small matrix
!
        call diag('V','U',itdav,qcmat,maxqcdiagsub,qceigen,nproc2,myrank2,mpi_comm2)
!
! Form correction vector
!
        qcnorm= zero
!$OMP parallel do reduction(+:qcnorm)
        do jj= 1,nocc*nvir+1
          qcvec(jj,itdav+1,1)= zero
          do ii= 1,itdav
            qcvec(jj,itdav+1,1)= qcvec(jj,itdav+1,1) &
&                               +qcmat(ii,1)*(qcvec(jj,ii,2)-qceigen(1)*qcvec(jj,ii,1))
          enddo
          qcnorm= qcnorm+qcvec(jj,itdav+1,1)*qcvec(jj,itdav+1,1)
        enddo
!$OMP end parallel do
        qcnorm= sqrt(qcnorm)
!
        if(master) write(*,'(5x,"QC Cycle ",i2," : Norm = ",1p,d10.3)') itdav-1, qcnorm
!
! Check convergence
!
        if(qcnorm < threshqc) exit
!
        qcvec(1,itdav+1,1)=-qcvec(1,itdav+1,1)/qceigen(1)
!$OMP parallel do private(ij)
        do ii= 1,nvir
          ij=(ii-1)*nocc+1
          do jj= 1,nocc
            qcvec(ij+jj,itdav+1,1)= qcvec(ij+jj,itdav+1,1)/ &
&                                  (fock(ii+nocc,ii+nocc)-fock(jj,jj)-qceigen(1))
          enddo
        enddo
!$OMP end parallel do
!
! Schmidt orthogonalization
!
        do ii= 1,itdav
          tmp= zero
          do jj= 1,nocc*nvir+1
            tmp= tmp-qcvec(jj,ii,1)*qcvec(jj,itdav+1,1)
          enddo
          do jj= 1,nocc*nvir+1
            qcvec(jj,itdav+1,1)= qcvec(jj,itdav+1,1)+tmp*qcvec(jj,ii,1)
          enddo
        enddo
!
! Renormalization of new vector
!
        qcnorm= sqrt(ddot(nocc*nvir+1,qcvec(1,itdav+1,1),1,qcvec(1,itdav+1,1),1))
!
! Check convergence
!
        if(qcnorm < threshqc) exit
!
! Reset Davidson diagonalization
!
        if(itdav == maxqcdiagsub) then
          do jj= 1,nocc*nvir+1
            qcvec(jj,1,1)= qcvec(jj,1,1)*qcmat(1,1)
          enddo
          do ii= 2,itdav
            do jj= 1,nocc*nvir+1
              qcvec(jj,1,1)= qcvec(jj,1,1)+qcvec(jj,ii,1)*qcmat(ii,1)
            enddo
          enddo
!    
          itdav=1
          cycle
        endif
!
        if(itqcdiag ==(maxqcdiag)) then
          if(master) then
            write(*,'(" Error! Number of iteration for Quadratically convergent ",&
&                     "method exceeds maxqcdiag=",i3,".")') maxqcdiag
            write(*,'(" Set larger value for maxqcdiag in scf section.")')
          endif
          call iabort
        endif
!
        qcnorm= one/qcnorm
!$OMP parallel do
        do ii= 1,nocc*nvir+1
          qcvec(ii,itdav+1,1)= qcvec(ii,itdav+1,1)*qcnorm
        enddo
!$OMP end parallel do
        itdav= itdav+1
      enddo
!
! End of Davidson diagonalization
!
! Calculate orbital rotation matrix
!
!$OMP parallel private(rotqc,tmp)
!$OMP do
      do jj= 1,nocc*nvir+1
        qcvec(jj,1,1)= qcvec(jj,1,1)*qcmat(1,1)
      enddo
!$OMP end do
      do ii= 2,itdav
!$OMP do
        do jj= 1,nocc*nvir+1
          qcvec(jj,1,1)= qcvec(jj,1,1)+qcvec(jj,ii,1)*qcmat(ii,1)
        enddo
!$OMP end do
      enddo
!
      tmp= one/qcvec(1,1,1)
!$OMP do
      do jj= 2,nocc*nvir+1
        qcvec(jj,1,1)= qcvec(jj,1,1)*tmp
      enddo
!$OMP end do
!
! Rotate occupied MOs
!
!$OMP do
      do ii= 1,nocc
        do jj= 1,nvir
          rotqc= qcvec((jj-1)*nocc+ii+1,1,1)
          do kk= 1,nao
            cmo(kk,ii)= cmo(kk,ii)+rotqc*cmo(kk,jj+nocc)
          enddo
        enddo
      enddo
!$OMP end do
!$OMP end parallel
!
! Schmidt orthonormalization of new occupied MOs
!
!$OMP parallel private(ij,tmp)
!$OMP do
      do ii= 1,nao
        ij= ii*(ii-1)/2
        do jj= 1,ii
          work(jj,ii)= overlap(ij+jj)
          work(ii,jj)= overlap(ij+jj)
        enddo
      enddo
!$OMP end do
!
      do ii= 1,nmo
!$OMP do
        do jj= 1,nao
          qcwork(jj,1)= zero
          do kk= 1,nao
            qcwork(jj,1)= qcwork(jj,1)+work(kk,jj)*cmo(kk,ii)
          enddo
        enddo
!$OMP end do
        tmp= zero
        do kk= 1,nao
          tmp= tmp+qcwork(kk,1)*cmo(kk,ii)
        enddo
        tmp= one/sqrt(tmp)
!$OMP barrier
!$OMP do
        do kk= 1,nao
          cmo(kk,ii)= cmo(kk,ii)*tmp
          qcwork(kk,1)= qcwork(kk,1)*tmp
        enddo
!$OMP end do
        if(ii == nmo) cycle
!
!$OMP do
        do jj= ii+1,nmo
          tmp= zero
          do kk= 1,nao
            tmp= tmp-qcwork(kk,1)*cmo(kk,jj)
          enddo
          do kk= 1,nao
            cmo(kk,jj)= cmo(kk,jj)+tmp*cmo(kk,ii)
          enddo
        enddo
!$OMP end do
      enddo
!$OMP end parallel
!
      return
end


!---------------------------------------------------------------------------------------------
  subroutine uhfqc(focka,fockb,cmoa,cmob,qcrmax,qcgmna,qcgmnb,qcvec, &
&                  qcmat,qcmatsave,qceigen,overlap,xint, &
&                  qcworka,qcworkb,work,hfexchange,nao,nmo,nocca,noccb,nvira,nvirb,nshell, &
&                  maxdim,maxqcdiag,maxqcdiagsub,threshqc, &
&                  nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!---------------------------------------------------------------------------------------------
!
! Driver of Davidson diagonalization for quadratically convergent of UHF
!
      use modparallel, only : master
      implicit none
      integer,intent(in) :: nao, nmo, nocca, noccb, nvira, nvirb, nshell, maxdim, maxqcdiag
      integer,intent(in) :: maxqcdiagsub, nproc1, nproc2, myrank1, myrank2, mpi_comm1, mpi_comm2
      integer :: itdav, itqcdiag, ii, ij, jj, kk, ia, ib, istart, icount
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: overlap(nao*(nao+1)/2), xint(nshell*(nshell+1)/2)
      real(8),intent(in) :: hfexchange, threshqc
      real(8),intent(inout) :: qcrmax(nshell*(nshell+1)/2)
      real(8),intent(inout) :: qcgmna(nao*(nao+1)/2), qcgmnb(nao*(nao+1)/2)
      real(8),intent(inout) :: qcvec(nocca*nvira+noccb*nvirb+1,maxqcdiagsub+1,2)
      real(8),intent(inout) :: qcmat(maxqcdiagsub,maxqcdiagsub)
      real(8),intent(inout) :: qcmatsave(maxqcdiagsub*(maxqcdiagsub+1)/2)
      real(8),intent(inout) :: qceigen(maxqcdiagsub)
      real(8),intent(inout) :: qcworka(nao,nao), qcworkb(nao,nao), work(nao,nao)
      real(8),intent(inout) :: focka(nao,nao),fockb(nao,nao), cmoa(nao,nao), cmob(nao,nao)
      real(8) :: tmp, tmpa, tmpb, qcnorm, ddot, rotqc
!
      qcmatsave(:)= zero
!
! Calculate Fock matrix in MO basis
!
      call expand(focka,qcworka,nao)
      call dsymm('L','U',nao,nmo,one,qcworka,nao,cmoa,nao,zero,work,nao)
      call dgemm('T','N',nmo,nmo,nao,one,cmoa,nao,work,nao,zero,focka,nao)
      call expand(fockb,qcworka,nao)
      call dsymm('L','U',nao,nmo,one,qcworka,nao,cmob,nao,zero,work,nao)
      call dgemm('T','N',nmo,nmo,nao,one,cmob,nao,work,nao,zero,fockb,nao)
!
! Calculate initial (first and second) qc vector
!
      qcvec(1,1,1)= one
      qcvec(2:nocca*nvira+noccb*nvirb+1,1,1)= zero
      qcvec(1,1,2)= zero
      qcvec(1,2,1)= zero
!
      qcnorm= zero
!$OMP parallel private(ij) reduction(+:qcnorm)
!$OMP do
      do ii= 1,nvira
        ij=(ii-1)*nocca+1
        do jj= 1,nocca
          qcvec(ij+jj,1,2)= focka(jj,ii+nocca)
          qcvec(ij+jj,2,1)= focka(jj,ii+nocca)
          qcnorm= qcnorm+focka(jj,ii+nocca)*focka(jj,ii+nocca)
        enddo
      enddo
!$OMP end do
!$OMP do
      do ii= 1,nvirb
        ij=(ii-1)*noccb+nocca*nvira+1
        do jj= 1,noccb
          qcvec(ij+jj,1,2)= fockb(jj,ii+noccb)
          qcvec(ij+jj,2,1)= fockb(jj,ii+noccb)
          qcnorm= qcnorm+fockb(jj,ii+noccb)*fockb(jj,ii+noccb)
        enddo
      enddo
!$OMP end do
!$OMP end parallel
      qcnorm= one/sqrt(qcnorm)
!$OMP parallel do
      do ii= 2,nocca*nvira+noccb*nvirb+1
        qcvec(ii,2,1)= qcvec(ii,2,1)*qcnorm
      enddo
!$OMP end parallel do

!
! Start Davidson diagonalization
!
      itdav= 2
      do itqcdiag= 2,maxqcdiag
!
! Calculate Gmn
!
        call calcqcurmn(qcworka,qcworkb,qcvec,cmoa,cmob,work,nao,nocca,noccb,nvira,nvirb, &
&                       itdav,maxqcdiagsub)
        call calcudmax(qcworka,qcworkb,qcrmax,work,nproc2,myrank2,mpi_comm2)
        call calcqcugmn(qcgmna,qcgmnb,work,qcworka,qcworkb,qcrmax,xint,hfexchange,maxdim, &
&                       nao,nshell,nproc1,myrank1,mpi_comm1)
!
! Add two-electron integral contribution
!
        call expand(qcgmna,qcworka,nao)
        call dsymm('L','U',nao,nvira,one,qcworka,nao,cmoa(1,nocca+1),nao,zero,work,nao)
        call dgemm('T','N',nocca,nvira,nao,one,cmoa,nao,work,nao,zero,qcvec(2,itdav,2),nocca)
        call expand(qcgmnb,qcworka,nao)
        call dsymm('L','U',nao,nvirb,one,qcworka,nao,cmob(1,noccb+1),nao,zero,work,nao)
        call dgemm('T','N',noccb,nvirb,nao,one,cmob,nao,work,nao,zero, &
&                  qcvec(nocca*nvira+2,itdav,2),noccb)
!
! Add Fock matrix element contribution
!
        tmp= zero
!$OMP parallel private(kk) reduction(+:tmp)
!$OMP do collapse(2)
        do ia= 1,nvira
          do ii= 1,nocca
            kk= (ia-1)*nocca+ii+1
            tmp= tmp+focka(ii,ia+nocca)*qcvec(kk,itdav,1)
            qcvec(kk,itdav,2)= qcvec(kk,itdav,2)+focka(ii,ia+nocca)*qcvec(1,itdav,1)
            do ib= 1,nvira
              qcvec(kk,itdav,2)= qcvec(kk,itdav,2)+focka(ib+nocca,ia+nocca) &
&                                                 *qcvec((ib-1)*nocca+ii+1,itdav,1)
            enddo
            do ij= 1,nocca
              qcvec(kk,itdav,2)= qcvec(kk,itdav,2)-focka(ij,ii)*qcvec((ia-1)*nocca+ij+1,itdav,1)
            enddo
          enddo
        enddo
!$OMP end do
!$OMP do collapse(2)
        do ia= 1,nvirb
          do ii= 1,noccb
            kk= (ia-1)*noccb+ii+nocca*nvira+1
            tmp= tmp+fockb(ii,ia+noccb)*qcvec(kk,itdav,1)
            qcvec(kk,itdav,2)= qcvec(kk,itdav,2)+fockb(ii,ia+noccb)*qcvec(1,itdav,1)
            do ib= 1,nvirb
              qcvec(kk,itdav,2)= qcvec(kk,itdav,2)+fockb(ib+noccb,ia+noccb) &
&                                                 *qcvec((ib-1)*noccb+ii+nocca*nvira+1,itdav,1)
            enddo
            do ij= 1,noccb
              qcvec(kk,itdav,2)= qcvec(kk,itdav,2)-fockb(ij,ii) &
&                                                 *qcvec((ia-1)*noccb+ij+nocca*nvira+1,itdav,1)
            enddo
          enddo
        enddo
!$OMP end do
!$OMP end parallel
        qcvec(1,itdav,2)= tmp
!
! Calculate small matrix <b|A|b>
!
        qceigen(:)= zero
!$OMP parallel reduction(+:qceigen)
        do ii= 1,itdav
!$OMP do
          do jj= 1,nocca*nvira+noccb*nvirb+1
            qceigen(ii)= qceigen(ii)+qcvec(jj,ii,1)*qcvec(jj,itdav,2)
          enddo
!$OMP end do
        enddo
!$OMP end parallel
!
        istart= itdav*(itdav-1)/2
        do ii= 1,itdav
          qcmatsave(istart+ii)= qceigen(ii)
        enddo
        icount= 0
        do ii= 1,itdav
          do jj= 1,ii
            icount= icount+1
            qcmat(jj,ii)= qcmatsave(icount)
          enddo
        enddo
!
! Diagonalize small matrix
!
        call diag('V','U',itdav,qcmat,maxqcdiagsub,qceigen,nproc2,myrank2,mpi_comm2)
!
! Form correction vector
!
        qcnorm= zero
!$OMP parallel do reduction(+:qcnorm)
        do jj= 1,nocca*nvira+noccb*nvirb+1
          qcvec(jj,itdav+1,1)= zero
          do ii= 1,itdav
            qcvec(jj,itdav+1,1)= qcvec(jj,itdav+1,1) &
&                               +qcmat(ii,1)*(qcvec(jj,ii,2)-qceigen(1)*qcvec(jj,ii,1))
          enddo
          qcnorm= qcnorm+qcvec(jj,itdav+1,1)*qcvec(jj,itdav+1,1)
        enddo
!$OMP end parallel do
        qcnorm= sqrt(qcnorm)
!
        if(master) write(*,'(5x,"QC Cycle ",i2," : Norm = ",1p,d10.3)') itdav-1, qcnorm
!
! Check convergence
!
        if(qcnorm < threshqc) exit
!
        qcvec(1,itdav+1,1)=-qcvec(1,itdav+1,1)/qceigen(1)
!$OMP parallel private(ij)
!$OMP do
        do ii= 1,nvira
          ij=(ii-1)*nocca+1
          do jj= 1,nocca
            qcvec(ij+jj,itdav+1,1)= qcvec(ij+jj,itdav+1,1)/ &
&                                  (focka(ii+nocca,ii+nocca)-focka(jj,jj)-qceigen(1))
          enddo
        enddo
!$OMP end do
!$OMP do
        do ii= 1,nvirb
          ij=(ii-1)*noccb+nocca*nvira+1
          do jj= 1,noccb
            qcvec(ij+jj,itdav+1,1)= qcvec(ij+jj,itdav+1,1)/ &
&                                  (fockb(ii+noccb,ii+noccb)-fockb(jj,jj)-qceigen(1))
          enddo
        enddo
!$OMP end do
!$OMP end parallel
!
! Schmidt orthogonalization
!
        do ii= 1,itdav
          tmp= zero
          do jj= 1,nocca*nvira+noccb*nvirb+1
            tmp= tmp-qcvec(jj,ii,1)*qcvec(jj,itdav+1,1)
          enddo
          do jj= 1,nocca*nvira+noccb*nvirb+1
            qcvec(jj,itdav+1,1)= qcvec(jj,itdav+1,1)+tmp*qcvec(jj,ii,1)
          enddo
        enddo
!
! Renormalization of new vector
!
        qcnorm= sqrt(ddot(nocca*nvira+noccb*nvirb+1,qcvec(1,itdav+1,1),1,qcvec(1,itdav+1,1),1))
!
! Check convergence
!
        if(qcnorm < threshqc) exit
!
        if(itqcdiag ==(maxqcdiag)) then
          if(master) then
            write(*,'(" Error! Number of iteration for Quadratically convergent ",&
&                     "method exceeds maxqcdiag=",i3,".")') maxqcdiag
            write(*,'(" Set larger value for maxqcdiag in scf section.")')
          endif
          call iabort
        endif
!
! Reset Davidson diagonalization
!
        if(itdav == maxqcdiagsub) then
          do jj= 1,nocca*nvira+noccb*nvirb+1
            qcvec(jj,1,1)= qcvec(jj,1,1)*qcmat(1,1)
          enddo
          do ii= 2,itdav
            do jj= 1,nocca*nvira+noccb*nvirb+1
              qcvec(jj,1,1)= qcvec(jj,1,1)+qcvec(jj,ii,1)*qcmat(ii,1)
            enddo
          enddo
!
          itdav=1
          cycle
        endif
!
        qcnorm= one/qcnorm
!$OMP parallel do
        do ii= 1,nocca*nvira+noccb*nvirb+1
          qcvec(ii,itdav+1,1)= qcvec(ii,itdav+1,1)*qcnorm
        enddo
!$OMP end parallel do
        itdav= itdav+1
      enddo
!
! End of Davidson diagonalization
!
! Calculate orbital rotation matrix
!
!$OMP parallel private(rotqc,tmp)
!$OMP do
      do jj= 1,nocca*nvira+noccb*nvirb+1
        qcvec(jj,1,1)= qcvec(jj,1,1)*qcmat(1,1)
      enddo
!$OMP end do
      do ii= 2,itdav
!$OMP do
        do jj= 1,nocca*nvira+noccb*nvirb+1
          qcvec(jj,1,1)= qcvec(jj,1,1)+qcvec(jj,ii,1)*qcmat(ii,1)
        enddo
!$OMP end do
      enddo
!
      tmp= one/qcvec(1,1,1)
!$OMP do
      do jj= 2,nocca*nvira+noccb*nvirb+1
        qcvec(jj,1,1)= qcvec(jj,1,1)*tmp
      enddo
!$OMP end do
!
! Rotate occupied MOs
!
!$OMP do
      do ii= 1,nocca
        do jj= 1,nvira
          rotqc= qcvec((jj-1)*nocca+ii+1,1,1)
          do kk= 1,nao
            cmoa(kk,ii)= cmoa(kk,ii)+rotqc*cmoa(kk,jj+nocca)
          enddo
        enddo
      enddo
!$OMP end do
!$OMP do
      do ii= 1,noccb
        do jj= 1,nvirb
          rotqc= qcvec((jj-1)*noccb+ii+nocca*nvira+1,1,1)
          do kk= 1,nao
            cmob(kk,ii)= cmob(kk,ii)+rotqc*cmob(kk,jj+noccb)
          enddo
        enddo
      enddo
!$OMP end do
!$OMP end parallel
!
! Schmidt orthonormalization of new occupied MOs
!
!$OMP parallel private(ij,tmpa,tmpb)
!$OMP do
      do ii= 1,nao
        ij= ii*(ii-1)/2
        do jj= 1,ii
          work(jj,ii)= overlap(ij+jj)
          work(ii,jj)= overlap(ij+jj)
        enddo
      enddo
!$OMP end do
!
      do ii= 1,nmo
!$OMP do
        do jj= 1,nao
          qcworka(jj,1)= zero
          qcworkb(jj,1)= zero
          do kk= 1,nao
            qcworka(jj,1)= qcworka(jj,1)+work(kk,jj)*cmoa(kk,ii)
            qcworkb(jj,1)= qcworkb(jj,1)+work(kk,jj)*cmob(kk,ii)
          enddo
        enddo
!$OMP end do
        tmpa= zero
        tmpb= zero
        do kk= 1,nao
          tmpa= tmpa+qcworka(kk,1)*cmoa(kk,ii)
          tmpb= tmpb+qcworkb(kk,1)*cmob(kk,ii)
        enddo
        tmpa= one/sqrt(tmpa)
        tmpb= one/sqrt(tmpb)
!$OMP barrier
!$OMP do
        do kk= 1,nao
          cmoa(kk,ii)= cmoa(kk,ii)*tmpa
          cmob(kk,ii)= cmob(kk,ii)*tmpb
          qcworka(kk,1)= qcworka(kk,1)*tmpa
          qcworkb(kk,1)= qcworkb(kk,1)*tmpb
        enddo
!$OMP end do
        if(ii == nmo) cycle
!
!$OMP do
        do jj= ii+1,nmo
          tmpa= zero
          tmpb= zero
          do kk= 1,nao
            tmpa= tmpa-qcworka(kk,1)*cmoa(kk,jj)
            tmpb= tmpb-qcworkb(kk,1)*cmob(kk,jj)
          enddo
          do kk= 1,nao
            cmoa(kk,jj)= cmoa(kk,jj)+tmpa*cmoa(kk,ii)
            cmob(kk,jj)= cmob(kk,jj)+tmpb*cmob(kk,ii)
          enddo
        enddo
!$OMP end do
      enddo
!$OMP end parallel
!
      return
end


!---------------------------------------------------------------------------------
  subroutine calcqcugmn(gmn1,gmn2,gmn3,rmna,rmnb,rmtrx,xint,hfexchange,maxdim, &
&                       nao,nshell,nproc,myrank,mpi_comm)
!---------------------------------------------------------------------------------
!
! Driver of Gmn matrix formation from two-electron intgrals
!
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: nshell, nao, maxdim, nproc, myrank, mpi_comm
      integer :: ijsh, ish, jsh, ksh, lsh, ij, kl, ik, il, jk, jl
      integer :: ii, jj, kk, kstart, ishcheck
      integer(8) :: ncount, icount(nshell)
      real(8),parameter :: zero=0.0D+00, two=2.0D+00
      real(8),intent(in) :: rmna(nao*(nao+1)/2), rmnb(nao*(nao+1)/2)
      real(8),intent(in) :: rmtrx(nshell*(nshell+1)/2), xint(nshell*(nshell+1)/2)
      real(8),intent(in) :: hfexchange
      real(8),intent(out) :: gmn1(nao*(nao+1)/2), gmn2(nao*(nao+1)/2), gmn3(nao*(nao+1)/2)
      real(8) :: xijkl, rmax, twoeri(maxdim**4), rmax1
      integer :: last, ltmp(nshell), lnum, ll
!
      gmn2(:)= zero
      gmn3(:)= zero
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
!$OMP private(ijsh,jsh,ksh,lsh,ij,kl,ik,il,jk,jl,xijkl,rmax1,rmax,twoeri,jj,kk, &
!$OMP kstart,last,ltmp,lnum,ll) firstprivate(ish,ii) reduction(+:gmn2,gmn3)
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
          rmax1=max(two*rmtrx(ij),rmtrx(ik),rmtrx(jk))
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
            rmax=max(rmax1,two*rmtrx(kl),rmtrx(il),rmtrx(jl))
            if(xijkl*rmax.ge.cutint2) then
              lnum=lnum+1
              ltmp(lnum)=lsh
            endif
          enddo
          do lsh= ll+1,last
            kl= kk+lsh
            il= ii+lsh
            jl= lsh*(lsh-1)/2+jsh
            xijkl=xint(ij)*xint(kl)
            rmax=max(rmax1,two*rmtrx(kl),rmtrx(il),rmtrx(jl))
            if(xijkl*rmax.ge.cutint2) then
              lnum=lnum+1
              ltmp(lnum)=lsh
            endif
          enddo
          do lsh= 1,lnum
            call calc2eri(twoeri,ish,jsh,ksh,ltmp(lsh),maxdim)
            call ugmneri(gmn2,gmn3,rmna,rmnb,twoeri,ish,jsh,ksh,ltmp(lsh),maxdim,hfexchange)
          enddo
        enddo
      enddo
!$OMP end parallel do
!
      do ii= 1,nao
        gmn2(ii*(ii+1)/2)= gmn2(ii*(ii+1)/2)*two
        gmn3(ii*(ii+1)/2)= gmn3(ii*(ii+1)/2)*two
      enddo
!
      call para_allreducer(gmn2,gmn1,nao*(nao+1)/2,mpi_comm)
      call para_allreducer(gmn3,gmn2,nao*(nao+1)/2,mpi_comm)
      return
end


!-----------------------------------------------------------------------------------
  subroutine ugmneri(gmna,gmnb,rmna,rmnb,twoeri,ish,jsh,ksh,lsh,maxdim,hfexchange)
!-----------------------------------------------------------------------------------
!
! Form Gmn matrix from two-electron intgrals
!
      use modbasis, only : nao, mbf, locbf
      use modthresh, only : cutint2
      implicit none
      integer,intent(in) :: ish, jsh, ksh, lsh, maxdim
      integer :: nbfi, nbfj, nbfk, nbfl
      integer :: locbfi, locbfj, locbfk, locbfl, jmax, lmax, i, j, k, l, ij, kl
      integer :: nij, nkl, nik, nil, njk, njl
      integer :: iloc, jloc, kloc, lloc, iloc2, jloc2, kloc2, lloc2, kloc0, jloc0
      real(8),parameter :: half=0.5D+00, two=2.0D+00
      real(8),intent(in) :: rmna(nao*(nao+1)/2), rmnb(nao*(nao+1)/2)
      real(8),intent(in) :: twoeri(maxdim,maxdim,maxdim,maxdim), hfexchange
      real(8),intent(inout) :: gmna(nao*(nao+1)/2), gmnb(nao*(nao+1)/2)
      real(8) :: val, val2
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
                val = val*hfexchange
                gmna(nij)= gmna(nij)+val2*rmna(nkl)+val2*rmnb(nkl)
                gmna(nkl)= gmna(nkl)+val2*rmna(nij)+val2*rmnb(nij)
                gmna(nik)= gmna(nik)-val *rmna(njl)
                gmna(nil)= gmna(nil)-val *rmna(njk)
                gmna(njk)= gmna(njk)-val *rmna(nil)
                gmna(njl)= gmna(njl)-val *rmna(nik)
                gmnb(nij)= gmnb(nij)+val2*rmnb(nkl)+val2*rmna(nkl)
                gmnb(nkl)= gmnb(nkl)+val2*rmnb(nij)+val2*rmna(nij)
                gmnb(nik)= gmnb(nik)-val *rmnb(njl)
                gmnb(nil)= gmnb(nil)-val *rmnb(njk)
                gmnb(njk)= gmnb(njk)-val *rmnb(nil)
                gmnb(njl)= gmnb(njl)-val *rmnb(nik)
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
                val = val*hfexchange
                gmna(nij)= gmna(nij)+val2*rmna(nkl)+val2*rmnb(nkl)
                gmna(nkl)= gmna(nkl)+val2*rmna(nij)+val2*rmnb(nij)
                gmna(nik)= gmna(nik)-val *rmna(njl)
                gmna(nil)= gmna(nil)-val *rmna(njk)
                gmna(njk)= gmna(njk)-val *rmna(nil)
                gmna(njl)= gmna(njl)-val *rmna(nik)
                gmnb(nij)= gmnb(nij)+val2*rmnb(nkl)+val2*rmna(nkl)
                gmnb(nkl)= gmnb(nkl)+val2*rmnb(nij)+val2*rmna(nij)
                gmnb(nik)= gmnb(nik)-val *rmnb(njl)
                gmnb(nil)= gmnb(nil)-val *rmnb(njk)
                gmnb(njk)= gmnb(njk)-val *rmnb(nil)
                gmnb(njl)= gmnb(njl)-val *rmnb(nik)
              enddo
            enddo kloop
          enddo
        enddo
      endif
!
      return
end

