!----------------------------------------------------
  subroutine calcdmtrx(cmo,dmtrx,work,ndim,neleca)
!----------------------------------------------------
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
!$OMP parallel do private(ij)
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
!$OMP parallel do private(ij)
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
!$OMP parallel do private(ij)
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
      integer,intent(in) :: nproc, myrank
      integer(4),intent(in) :: mpi_comm
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
      integer,intent(in) :: nproc, myrank
      integer(4),intent(in) :: mpi_comm
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


!----------------------------------------------------
  subroutine ddiff(dmtrx,dmtrxprev,work,ndim,diffmax)
!----------------------------------------------------
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
      integer :: i, imax, idamax
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: dmtrx(ndim), dmtrxprev(ndim)
      real(8),intent(out) :: diffmax, work(ndim)
!
      do i= 1,ndim
        work(i)= dmtrx(i)-dmtrxprev(i)
      enddo
      imax= idamax(ndim,work,1)
      diffmax=abs(work(imax))
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
      integer,intent(in) :: nao, maxdiis, nproc, myrank, idis(nproc,14) 
      integer(4),intent(in) :: mpi_comm
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
! In    : fock    (Interporated Fock matrx (out))
!         dmtrx   (Density matrix)
!         overlap (Overlap matrix)
!         ortho   (Orthogonalization matrix)
! Out   : work1   (Error matrix)
!
      implicit none
      integer,intent(in) :: nao, nmo, nproc, myrank, idis(nproc,14)
      integer(4),intent(in) :: mpi_comm
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
      do i= 1,nmo
        do j= 1,i
          diff= work1(j,i)-work1(i,j)
          work1(j,i)= diff
          work1(i,j)=-diff
          if(abs(diff) > errmax) errmax= abs(diff)
        enddo
      enddo
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
! In    : fock    (Fock matrx)
! Out   : fock    (Interporated Fock matrx)
! Inout : fockdiis(Previous Fock matrix)
!         errdiis (DIIS error matrix)
!
      implicit none
      integer,intent(in) :: itdiis, nao, maxdiis, nproc, myrank, idis(nproc,14)
      integer(4),intent(in) :: mpi_comm
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


!---------------------------------------------------------------------------
  subroutine calcudiis(focka,fockb,errdiisa,errdiisb,fockdiisa,fockdiisb, &
&                      diismtrx,worka,workb,work2,itdiis,nao,maxdiis, &
&                      idis,nproc,myrank,mpi_comm)
!---------------------------------------------------------------------------
!
! Direct Inversion in the Iterative Subspace (DIIS) interporation for open-shell
!
! In    : focka    (Alpha Fock matrx)
!         fockb    (Beta Fock matrx)
!         worka    (Alpha DIIS error matrix)
!         workb    (Beta DIIS error matrix)
! Out   : focka    (Interporated alpha Fock matrx)
!         fockb    (Interporated beta Fock matrx)
! Inout : fockdiisa(Previous alpha Fock matrix)
!         fockdiisb(Previous beta Fock matrix)
!         errdiisa (History of Alpha DIIS error matrix)
!         errdiisb (History of Beta DIIS error matrix)
!
      implicit none
      integer,intent(in) :: itdiis, nao, maxdiis, nproc, myrank, idis(nproc,14)
      integer(4),intent(in) :: mpi_comm
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
      integer,intent(in) :: nocc, nvir, nao, itype, nproc, myrank, idis(nproc,14)
      integer(4),intent(in) :: mpi_comm
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
common/ishimura/t4
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
&    +dot_product(sodispb(1:noccb*nvirb,itsoscf),sodispb(1:noccb*nvirb,itsoscf))
   t1=sqrt(t1)
      ijmo= idamax(nocca*nvira,sodispa(1,itsoscf),1)
      it= idamax(noccb*nvirb,sodispb(1,itsoscf),1)
      t4=max(sodispa(ijmo,itsoscf),sodispb(it,itsoscf))
     
!  if(t1.gt.1.0)then
!rite(*,*)"ishi",t1
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
      integer,intent(in) :: nocc, nvir, itsoscf, maxsoscf, nao, nmo, nproc, myrank, idis(nproc,14)
      integer(4),intent(in) :: mpi_comm
      integer :: imo, jmo, num, istart
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, p01= 1.0D-01
      real(8),intent(in) :: sodisp(nocc,nvir,maxsoscf), sodispmax
      real(8),intent(inout) :: cmo(nao,nao)
      real(8),intent(out) :: work(nmo,nmo), work2(*)
      real(8) :: factor
!
! U = exp(A) = I + A + A*A/2
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
      integer,intent(in) :: neleca, nelecb, nao, nproc, myrank, idis(nproc,14)
      integer(4),intent(in) :: mpi_comm
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
      ij= 0
      s2= zero
      do i= 1,nao
        do j= 1,i-1
          ij= ij+1
           s2= s2+dmtrxb(ij)*work2(j,i)
        enddo
        ij= ij+1
        s2= s2+dmtrxb(ij)*work2(i,i)*half
      enddo
      s2= s2*two
!
      sz= half*(neleca-nelecb)
      s2= sz*sz+half*(neleca+nelecb)-s2
!
      return
end


!--------------------------------------------------
  subroutine writeeigenvalue(eigena,eigenb,itype)
!--------------------------------------------------
!
! Write eigenvalues
!
      use modparallel, only : master
      use modmolecule, only : nmo, neleca, nelecb
      use modprint, only : iprint
      integer,intent(in) :: itype
      integer :: imo
      real(8),intent(in) :: eigena(nmo), eigenb(nmo)
!
! Closed-shell
!
      if(master.and.(iprint >= 1)) then
        if(itype == 1) then
          write(*,'(1x,80("-"))')
          write(*,'("   Eigenvalues (Hartree)")')
          write(*,'(1x,80("-"))')
          write(*,'("   Alpha Occupied: ",5f12.5)')(eigena(imo),imo=1,neleca)
          write(*,'("   Alpha Virtual : ",5f12.5)')(eigena(imo),imo=neleca+1,nmo)
          write(*,'(1x,80("-"))')
!
! Open-shell
!
        elseif(itype == 2) then
          write(*,'(1x,80("-"))')
          write(*,'("   Eigenvalues (Hartree)")')
          write(*,'(1x,80("-"))')
          write(*,'("   Alpha Occupied: ",5f12.5)')(eigena(imo),imo=1,neleca)
          write(*,'("   Alpha Virtual : ",5f12.5)')(eigena(imo),imo=neleca+1,nmo)
          write(*,'("   Beta  Occupied: ",5f12.5)')(eigenb(imo),imo=1,nelecb)
          write(*,'("   Beta  Virtual : ",5f12.5)')(eigenb(imo),imo=nelecb+1,nmo)
          write(*,'(1x,80("-"))')
        endif
      endif
!
      return
end


!-----------------------------------------
  subroutine writeeigenvector(cmo,eigen)
!-----------------------------------------
!
! Write eigenvalues
!
      use modparallel, only : master
      use modmolecule, only : nmo, neleca, numatomic
      use modbasis, only : nao, nshell, mtype, spher, locatom
      use modprint, only : iprint
      use modparam, only : mxao
      implicit none
      integer :: maxmo, imin, imax, ii, jj, kk, iao, iatom
      real(8),intent(in) :: cmo(nao,nao), eigen(nmo)
      character(len=8) :: atomlabel(mxao)
      character(len=5) :: bflabel(mxao)
      character(len=3) :: table(112)= &
&     (/'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ','Al ','Si ','P  ',&
&       'S  ','Cl ','Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',&
&       'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ',&
&       'Pd ','Ag ','Cd ','In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',&
&       'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ','Ta ','W  ','Re ',&
&       'Os ','Ir ','Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',&
&       'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ','Md ','No ','Lr ','Rf ','Db ',&
&       'Sg ','Bh ','Hs ','Mt ','Uun','Uuu','Uub'/)
!
      character(len=5) :: anglabel(56)= &
&     (/'S    ','Px   ','Py   ','Pz   ','Dxx  ','Dyy  ','Dzz  ','Dxy  ','Dxz  ','Dyz  ', &
&       'D 0  ','D+1  ','D-1  ','D+2  ','D-2  ','Fxxx ','Fyyy ','Fzzz ','Fxxy ','Fxxz ', &
&       'Fxyy ','Fyyz ','Fxzz ','Fyzz ','Fxyz ','F 0  ','F+1  ','F-1  ','F+2  ','F-2  ', &
&       'F+3  ','F-3  ','Gxxxx','Gyyyy','Gzzzz','Gxxxy','Gxxxz','Gxyyy','Gyyyz','Gxzzz', &
&       'Gyzzz','Gxxyy','Gxxzz','Gyyzz','Gxxyz','Gxyyz','Gxyzz','G 0  ','G+1  ','G-1  ', &
&       'G+2  ','G-2  ','G+3  ','G-3  ','G+4  ','G-4  '/)
!
      if(maxval(mtype(1:nshell)) >= 5) then
        if(master) write(*,'(" Sorry! This program can not display MOs of h functions now.")')
        return
      endif
!
      atomlabel(1:nao)= ''
      iao= 1
      iatom= 0
      do ii= 1,nshell
        select case(mtype(ii))
          case(0)
            bflabel(iao)= anglabel(1)
            if(locatom(ii) /= iatom) then
              iatom= locatom(ii)
              write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
            endif
            iao= iao+1
          case(1)
            bflabel(iao:iao+2)= anglabel(2:4)
            if(locatom(ii) /= iatom) then
              iatom= locatom(ii)
              write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
            endif
            iao= iao+3
          case(2)
            if(spher) then
              bflabel(iao:iao+4)= anglabel(11:15)
              if(locatom(ii) /= iatom) then
                iatom= locatom(ii)
                write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
              endif
              iao= iao+5
            else
              bflabel(iao:iao+5)= anglabel(5:10)
              if(locatom(ii) /= iatom) then
                iatom= locatom(ii)
                write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
              endif
              iao= iao+6
            endif
          case(3)
            if(spher) then
              bflabel(iao:iao+6)= anglabel(26:32)
              if(locatom(ii) /= iatom) then
                iatom= locatom(ii)
                write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
              endif
              iao= iao+7
            else
              bflabel(iao:iao+9)= anglabel(16:25)
              if(locatom(ii) /= iatom) then
                iatom= locatom(ii)
                write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
              endif
              iao= iao+10
            endif
          case(4)
            if(spher) then
              bflabel(iao:iao+8)= anglabel(48:56)
              if(locatom(ii) /= iatom) then
                iatom= locatom(ii)
                write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
              endif
              iao= iao+9
            else
              bflabel(iao:iao+14)= anglabel(33:47)
              if(locatom(ii) /= iatom) then
                iatom= locatom(ii)
                write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
              endif
              iao= iao+15
            endif
        end select
      enddo
!
      maxmo=min(nmo,neleca+20)
      imin= 1
      imax= 5
      if(master) then
        do ii= 1,(maxmo-1)/5+1
          if(imax > maxmo) imax= maxmo
          write(*,*)
          write(*,'(18x,5(5x,i4,2x))')(jj,jj=imin,imax)
          write(*,'(18x,5f11.4)')(eigen(jj),jj=imin,imax)
          do kk= 1,nao
            write(*,'(i5,a8,a5,5f11.4)')kk,atomlabel(kk),bflabel(kk),(cmo(kk,jj),jj=imin,imax)
          enddo
          imin= imin+5
          imax= imax+5
        enddo
        write(*,*)
      endif
!
      return
end








