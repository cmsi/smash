!----------------------------------------------------
  subroutine calcdmtrx(cmoa,dmtrx,work,ndim,neleca)
!----------------------------------------------------
!
! Calculate density matrix for closed shell
!
! In  : cmoa  (MO coefficients)
!       ndim  (dimension of basis functions)
!       neleca(number of alpha electrons)
! Out : dmtrx (density matrix)
!
      implicit none
      integer,intent(in) :: ndim, neleca
      integer :: i, j, ij
      real(8),parameter :: zero=0.0D+00, two=2.0D+00
      real(8),intent(in) :: cmoa(ndim*ndim)
      real(8),intent(out) :: dmtrx(ndim*(ndim+1)/2), work(ndim,ndim)
!
      call dgemm('N','T',ndim,ndim,neleca,two,cmoa,ndim,cmoa,ndim,zero,work,ndim)
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


!----------------------------------
  subroutine calcdmax(dmtrx,dmax)
!----------------------------------
!
! Calculate maximum density matrix element for each shell
!
! In  : dmtrx (density matrix elements)
! Out : dmax  (maximum density matrix elements)
!
      use procpar, only : nproc, myrank, MPI_SUM, MPI_COMM_WORLD
      use basis, only : nshell, nao, mbf, locbf
      implicit none
      integer :: ish, jsh, ijsh, locbfi, locbfj, nbfi, nbfj
      integer :: jnbf, i, j, ii, ij
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: dmtrx(nao*(nao+1)/2)
      real(8),intent(out) :: dmax(nshell*(nshell+1)/2)
      real(8) :: dtmp
!
      call zeroclr(dmax,nshell*(nshell+1)/2)
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
          dmax(ijsh)= dtmp
        enddo
!$OMP enddo
      enddo
!$OMP end parallel
!
      call para_allreduce(dmax,dtmp,nshell*(nshell+1)/2,"D",MPI_SUM,MPI_COMM_WORLD,1)
      return
end


!----------------------------------------------------
  subroutine ddiff(dmtrx,dmtrxprev,work,ndim,diffmax)
!----------------------------------------------------
!
! Calculate largest absolute change in density matrix
!
! In  : dmtrx    (density matrix elements)
!       dmtrxprev(previous density matrix elements)
!       ndim     (dimensions of dmtrx and dmtrxprev)
! Out : work     (change in density matrix)
!       diffmax  (largest absolute change in density matrix)
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


!------------------------------------------------------------------
  subroutine extrap(fock,fockwork,work1,work2,work3,idis,itextra)
!------------------------------------------------------------------
!
! Extrapolate Fock matrix
!
! In    : itextra (counter for extrapolation)
! InOut : fock    (Fock matrix)
!         fockwork(previous Fock and work matrix)
!
      use procpar
      use basis, only : nao
      use scf, only : maxdiis
      use iofile, only : iout
      implicit none
      integer,intent(in) :: idis(nproc,10)
      integer,intent(inout) :: itextra
      integer :: num, istart, i, iskip, nao3
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00, four=4.0D+00
      real(8),parameter :: pm10=1.0D-10, pm17=1.0D+17, pm7=1.0D-7, pm19=1.9D+00, pm99=0.99D+00
      real(8),intent(inout) :: fock(nao*(nao+1)/2)
      real(8),intent(inout) :: fockwork(idis(myrank+1,5),maxdiis)
      real(8),intent(out) :: work1(*), work2(*), work3(*)
      real(8) :: tridot
      real(8) :: sp11, sp12, sp13, sp22, sp23, sp33, dp1, dp2, dp3, cosphi, cospsi
      real(8) :: x, y, z, xy, xy2, xxx, yyy
!
      if(maxdiis < 6) then
        if(master) then
          write(iout,'(" Set Maxdiis for more than 6.")')
          call iabort
        endif
      endif
!
      num= idis(myrank+1,5)
      istart= idis(myrank+1,6)+1
      if(itextra.le.3) then
        if(num /= 0) call dcopy(num,fock(istart),1,fockwork(1,itextra),1)
        itextra= itextra+1
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
      call para_allgatherv(fockwork(1,4),num,'D',work1,idis(1,5),idis(1,6),MPI_COMM_WORLD)
      call para_allgatherv(fockwork(1,5),num,'D',work2,idis(1,5),idis(1,6),MPI_COMM_WORLD)
      call para_allgatherv(fockwork(1,6),num,'D',work3,idis(1,5),idis(1,6),MPI_COMM_WORLD)
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
!     ----- FIND COSINE OF ANGLE BETWEEN SUCCESSIVE DISPLACEMENTS -----
!
      cosphi= sp12/(dp1*dp2)
!
      z= one/(sp11*sp22-sp12*sp12)
      if(abs(z).gt.pm17) iskip= 1
      x=(sp13*sp22-sp12*sp23)*z
      y=(sp23*sp11-sp12*sp13)*z
      cospsi= sqrt(x*x*sp11+y*y*sp22+two*x*y*sp12)/dp3
!
!     ----- DO NOT EXTRAPOLATE UNLESS -4- CONSECUTIVE POINTS ARE
!           NEARLY COPLANAR -----
!
      if(cospsi.le.pm7) iskip= 1
!
!     ----- EXPRESS -DP(1)- AS X*DP(3)(PROJECTED)+Y*DP(2) -----
!
      x= one/x
      y=-y*x
!
!     ----- TEST IF 2*2 MATRIX HAS REAL EIGENVALUES
!           BETWEEN -TOL/2 AND +TOL/2 -----
!
      xy= y*y+four*x
      if(xy.lt.zero) then
        iskip= 1
        xy2= zero
      else
        xy2= abs(y)+sqrt(xy)
      endif
      if(xy2.gt.pm19.and.abs(cosphi).le.pm99) iskip= 1
!
!     ----- IF -4- POINT EXTRAPOLATION IS NOT POSSIBLE,
!           TRY -3- POINT
!
      if(iskip.eq.0) then
        if(xy2.gt.pm19) then
          x= dp1/(dp2*cosphi-dp1)
          do i= 1,nao3
            fock(i)= fock(i)+x*work1(i)
          enddo
        else
          xxx= x/(one-x-y)
          yyy=(x+y)/(one-x-y)       
          do i= 1,nao3
            fock(i)= fock(i)+xxx*work2(i)+yyy*work1(i)
          enddo
        endif
        itextra=1
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


!----------------------------------------------------------------------------
  subroutine calcdiis(fock,dmtrx,overlap,ortho,errdiis,fockdiis,diismtrx, &
&                     work1,work2,work3,errmax,idis,itdiis) 
!----------------------------------------------------------------------------
!
! Direct Inversion in the Iterative Subspace (DIIS) interporation
!
      use procpar
      use basis, only : nao
      use scf, only : maxdiis
      use molecule, only : nmo
      use iofile, only : iout
      use thresh, only : thresherr
      implicit none
      integer,intent(in) :: idis(nproc,10)
      integer,intent(inout) :: itdiis
      integer :: num, istart, i, j, imax, ij, ipiv(maxdiis+1), info, idamax
      real(8),intent(in) :: dmtrx(nao*(nao+1)/2), overlap(nao*(nao+1)/2), ortho(nao,nao)
      real(8),intent(inout) :: fock(nao*(nao+1)/2), errdiis(idis(myrank+1,3),maxdiis)
      real(8),intent(inout) :: fockdiis(idis(myrank+1,5),maxdiis), diismtrx(maxdiis*(maxdiis+1)/2)
      real(8),intent(out) :: work1(nao,*), work2(*), work3(*), errmax
      real(8) :: diisb(maxdiis+1,maxdiis+1), diiscoeff(maxdiis+1), diff, ddot
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
!
      if(itdiis > maxdiis) then
        itdiis=1
        call zeroclr(diismtrx,maxdiis*(maxdiis+1)/2)
      endif
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
      call para_allgatherv(work3,idis(myrank+1,3),'D',work1,idis(1,3),idis(1,4),MPI_COMM_WORLD)
      do i=1,nmo
        do j=1,i
          diff= work1(j,i)-work1(i,j)
          work1(j,i)= diff
          work1(i,j)=-diff
        enddo
      enddo
!
      imax= idamax(nmo*nmo,work1,1)
      errmax=abs(work1(mod(imax-1,nmo)+1,(imax-1)/nmo+1))
!
!ishimura
      if((itdiis == 1).and.(errmax >= thresherr)) return
!
      if(num > 0) then
        call dcopy(num*nmo,work1(1,istart),1,errdiis(1,itdiis),1)
        do i= 1,itdiis
          diismtrx(itdiis*(itdiis-1)/2+i)= ddot(num*nmo,errdiis(1,i),1,errdiis(1,itdiis),1)
        enddo
      endif
      call para_allreduce(diismtrx(itdiis*(itdiis-1)/2+1),diff,itdiis,'D',MPI_SUM,MPI_COMM_WORLD,1)
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
        call zeroclr(work1,num)
        do i= 1,itdiis
          call daxpy(num,diiscoeff(i),fockdiis(1,i),1,work1,1)
        enddo
      endif
      call para_allgatherv(work1,num,'D',fock,idis(1,5),idis(1,6),MPI_COMM_WORLD)
!
      itdiis=itdiis+1
!
      return
end


!------------------------------------------------------------------------
  subroutine soscfgrad(work,work2,sograd,cmoa,nocc,nvir,sogradmax,idis)
!------------------------------------------------------------------------
!
! Calculate orbital gradient <occ|F|vir>
!
      use procpar
      use basis, only : nao
      implicit none
      integer,intent(in) :: nocc, nvir, idis(nproc,10)
      integer :: numwork, num, istart, isomax, idamax
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: cmoa(nao,nao)
      real(8),intent(out) :: sograd(nocc*nvir), sogradmax
      real(8),intent(inout) :: work(nao*nao), work2(*)
!
      num= idis(myrank+1,7)
      istart= idis(myrank+1,8)+1+nocc
      if(num > 0) then
        call dsymm('L','U',nao,num,one,work,nao,cmoa(1,istart),nao,zero,work2,nao)
        call dgemm('T','N',nocc,num,nao,one,cmoa,nao,work2,nao,zero,work,nocc)
      endif
!
      numwork= idis(myrank+1,9)
      call para_allgatherv(work,numwork,'D',sograd,idis(1,9),idis(1,10),MPI_COMM_WORLD)
!
! Obtain maximum sograd value
!
      isomax= idamax(nocc*nvir,sograd,1)
      sogradmax= abs(sograd(isomax))
      return
end


!-----------------------------------------------
  subroutine soscfhess(hstart,eigen,nocc,nvir)
!-----------------------------------------------
!
! Set initial inverse Hessian matrix for SOSCF
!
      use basis, only : nao
      implicit none
      integer,intent(in) :: nocc, nvir
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


!------------------------------------------------------------------------
  subroutine soscfupdate(cmoa,hstart,sograd,sodisp,sovecy,work,work2, &
&                        nocc,nvir,itsoscf,maxsoscf,idis)
!------------------------------------------------------------------------
!
! Update molecular orbitals using approximated SOSCF method
!
      use procpar
      use basis, only : nao
      use molecule, only : nmo
      implicit none
      integer,intent(in) :: nocc, nvir, itsoscf, maxsoscf, idis(nproc,10)
      integer :: imo, jmo, it, num, istart
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00
      real(8),intent(in) :: hstart(nocc,nvir)
      real(8),intent(inout) :: cmoa(nao,nao), sograd(nocc,nvir,maxsoscf)
      real(8),intent(inout) :: sodisp(nocc,nvir,maxsoscf), sovecy(nocc,nvir,maxsoscf-1)
      real(8),intent(inout) :: work(nmo,nmo), work2(*)
      real(8) :: s1, s2, s3, s4, s5, s6, s1s2, t1, t2, t3, t4
!
! Initialize displacement vector
!
!$OMP parallel do
      do imo= 1,nvir
        do jmo= 1,nocc
          sodisp(jmo,imo,itsoscf)= hstart(jmo,imo)*sograd(jmo,imo,itsoscf)
        enddo
      enddo
!$OMP end parallel do
!
      if(itsoscf >= 2) then
!$OMP parallel do
        do imo= 1,nvir
          do jmo= 1,nocc
            sograd(jmo,imo,itsoscf-1)= sograd(jmo,imo,itsoscf)-sograd(jmo,imo,itsoscf-1)
            sovecy(jmo,imo,itsoscf-1)=-hstart(jmo,imo)*sograd(jmo,imo,itsoscf-1)
          enddo
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
          do imo= 1,nvir
            do jmo= 1,nocc
              s1= s1+sodisp(jmo,imo,it)*sograd(jmo,imo,it)
              s2= s2+sograd(jmo,imo,it)*sovecy(jmo,imo,it)
              s3= s3+sodisp(jmo,imo,it)*sograd(jmo,imo,itsoscf)
              s4= s4+sovecy(jmo,imo,it)*sograd(jmo,imo,itsoscf)
              s5= s5+sodisp(jmo,imo,it)*sograd(jmo,imo,itsoscf-1)
              s6= s6+sovecy(jmo,imo,it)*sograd(jmo,imo,itsoscf-1)
            enddo
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
          do imo= 1,nvir
            do jmo= 1,nocc
              sodisp(jmo,imo,itsoscf)= sodisp(jmo,imo,itsoscf)-t1*sodisp(jmo,imo,it) &
&                                                             +t2*sovecy(jmo,imo,it)
              sovecy(jmo,imo,itsoscf-1)= sovecy(jmo,imo,itsoscf-1)+t3*sodisp(jmo,imo,it) &
&                                                                 -t4*sovecy(jmo,imo,it)
            enddo
          enddo
!$OMP end parallel do
        enddo
!
        s1= zero
        s2= zero
        s3= zero
        s4= zero
!$OMP parallel do reduction(+:s1,s2,s3,s4)
        do imo= 1,nvir
          do jmo= 1,nocc
            s1= s1+sodisp(jmo,imo,itsoscf-1)*sograd(jmo,imo,itsoscf-1)
            s2= s2+sograd(jmo,imo,itsoscf-1)*sovecy(jmo,imo,itsoscf-1)
            s3= s3+sodisp(jmo,imo,itsoscf-1)*sograd(jmo,imo,itsoscf)
            s4= s4+sovecy(jmo,imo,itsoscf-1)*sograd(jmo,imo,itsoscf)
          enddo
        enddo
!$OMP end parallel do
        s1= one/s1
        s2= one/s2
        s1s2= one+s1/s2
        t1= s1s2*s1*s3-s1*s4
        t2= s1*s3
!$OMP parallel do
        do imo= 1,nvir
          do jmo= 1,nocc
            sodisp(jmo,imo,itsoscf)= sodisp(jmo,imo,itsoscf)-t1*sodisp(jmo,imo,itsoscf-1) &
&                                                           +t2*sovecy(jmo,imo,itsoscf-1)
          enddo
        enddo
!$OMP end parallel do
      endif
!
! U = exp(A) = I + A + A*A/2
!
!ishimura
      work= zero
!      call dgemm('N','T',nocc,nocc,nvir,-half,sodisp(1,1,itsoscf),nocc,sodisp(1,1,itsoscf), &
!&                nocc,zero,work,nmo)
!      call dgemm('T','N',nvir,nvir,nocc,-half,sodisp(1,1,itsoscf),nocc,sodisp(1,1,itsoscf), &
!&                nmo,zero,work(nocc+1,nocc+1),nmo)
!$OMP parallel do
      do imo= 1,nvir
        do jmo= 1,nocc
          work(jmo,imo+nocc)= sodisp(jmo,imo,itsoscf)
          work(imo+nocc,jmo)=-sodisp(jmo,imo,itsoscf)
        enddo
      enddo
!$OMP end parallel do
      do imo= 1,nmo
        work(imo,imo)= work(imo,imo)+one
      enddo
!
      call orthonorm(work,nmo)
!
      num= idis(myrank+1,1)
      istart= idis(myrank+1,2)+1
      if(num > 0) call dgemm('N','N',nao,num,nmo,one,cmoa,nao,work(1,istart),nmo,zero, &
&                            work2,nao)
      call para_allgatherv(work2,idis(myrank+1,3),'D',cmoa,idis(1,3),idis(1,4),MPI_COMM_WORLD)
      return
end



















