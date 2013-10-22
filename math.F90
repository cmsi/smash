!----------------------------------
  subroutine zeroclr(array,isize)
!----------------------------------
!
! Zero clear of array
!
      implicit none
      integer :: isize, i
      real(8),parameter :: zero=0.0D+0
      real(8) :: array(isize)
!
!$OMP parallel do
      do i= 1,isize
        array(i)= zero
      enddo
!$OMP end parallel do
      return
end


!-------------------------------------------------------
  subroutine ghquad(xyzint, expgh, xyzpij, iang, jang)
!-------------------------------------------------------
!
! Calculate Gauss-Hermite quadrature
!
      use hermite, only : hnode, hweight, minh, maxh
      implicit none
      integer,intent(in) :: iang, jang
      integer :: nroot, ij, i, j
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: expgh, xyzpij(3,2)
      real(8),intent(out) :: xyzint(3)
      real(8) :: ghxyz(3), exnode, pxyz(3)
!
      do i= 1,3
        xyzint(i)= zero
      enddo
!
      nroot=(iang+jang)/2+1
      do ij= minh(nroot),maxh(nroot)
        do i= 1,3
          ghxyz(i)= hweight(ij)
        enddo
        exnode= hnode(ij)*expgh
        if(iang >= 1) then
          do i= 1,3
            pxyz(i)= exnode+xyzpij(i,1)
          enddo
          do i= 1,iang
            do j= 1,3
              ghxyz(j)= ghxyz(j)*pxyz(j)
            enddo
          enddo
        endif
        if(jang >= 1) then
          do i= 1,3
            pxyz(i)= exnode+xyzpij(i,2)
          enddo
          do i= 1,jang
            do j= 1,3
              ghxyz(j)= ghxyz(j)*pxyz(j)
            enddo
          enddo
        endif
        do i= 1,3
        xyzint(i)= xyzint(i)+ghxyz(i)
        enddo
      enddo
      return
end


!-------------------------------------------
 subroutine fullmtrx(trimat,fullmat,ndim)
!-------------------------------------------
!
! Copy triangler size matrix to full size matrix
!
      use guess, only : nao_g, nao_v
      implicit none
      integer,intent(in) :: ndim
      integer :: i, j, ii
      real(8),intent(in) :: trimat((ndim*(ndim+1))/2)
      real(8),intent(out) :: fullmat(ndim,ndim)
!
!$OMP parallel do private(ii)
      do i= 1,ndim
        ii= i*(i-1)/2
        do j= 1,i
          fullmat(j,i)= trimat(ii+j)
        enddo
      enddo
!$OMP end parallel do
      return
end


!-------------------------------------------------------
  subroutine mtrxcanon(ortho,overlap,eigen,ndim,newdim)
!-------------------------------------------------------
!
! Calculate canonicalization matrix
!
! The ortho matrix satisfiles (Ortho)-daggar * S * (Ortho) = I
! where S is the overlap matrix.
!
      use iofile,only : iout
      use procpar, only : master
      use thresh, only : threshover
      implicit none
      integer,intent(in) :: ndim
      integer,intent(out) :: newdim
      integer :: i, j, icount
      real(8),parameter :: one=1.0D+00
      real(8),intent(inout) :: overlap(ndim,ndim)
      real(8),intent(out) :: ortho(ndim,ndim), eigen(ndim)
      real(8) :: ecanon
!
! Diagonalize ortho matrix
!
      call diag('V','U',ndim,overlap,ndim,eigen)
!
! Eliminate eigenvectors with small eigenvalues
!
      icount= 0
      do i= 1,ndim
        if(eigen(i) < threshover) icount= icount+1
      enddo
      newdim= ndim-icount
!
! Form canonical orthonormal matrix
!
!$OMP parallel do private(ecanon)
      do i= 1,newdim
        ecanon= one/sqrt(eigen(i+icount))
        do j= 1,ndim
          ortho(j,i)=overlap(j,i+icount)*ecanon
        enddo
      enddo
!$OMP end parallel do
      return
end


!---------------------------------------------------
  subroutine canonicalize(hmat,ortho,work,nao,nmo)
!---------------------------------------------------
!
! Canonicalize symmetric Hamiltonian matrix
!   In : hmat(nao*nao), elements are in upper triangle.
!   Out: hmat(nmo*nmo)
!
      implicit none
      integer,intent(in) :: nao, nmo
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: ortho(nao*nmo)
      real(8),intent(inout) :: hmat(nao*nao), work(nao*nao)
!
      call dsymm('L','U',nao,nmo,one,hmat,nao,ortho,nao,zero,work,nao)
      call dgemm('T','N',nmo,nmo,nao,one,ortho,nao,work,nao,zero,hmat,nao)
      return
end


!---------------------------------------------------------------
  subroutine canonicalizep(hmat,ortho,work,work2,nao,nmo,idis)
!---------------------------------------------------------------
!
! Canonicalize symmetric Hamiltonian matrix
!   In : hmat(nao*nao), elements are in upper triangle.
!   Out: hmat(nmo*nmo), matrix size is nao*nao
!
      use procpar
      implicit none
      integer,intent(in) :: nao, nmo, idis(nproc,4)
      integer :: num, istart
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: ortho(nao,nmo)
      real(8),intent(inout) :: hmat(nao*nao), work(nao*nao), work2(*)
!
      num=idis(myrank+1,1)
      istart=idis(myrank+1,2)+1
      if(num > 0)then
        call dsymm('L','U',nao,num,one,hmat,nao,ortho(1,istart),nao,zero,work,nao)
        call dgemm('T','N',nmo,num,nao,one,ortho,nao,work,nao,zero,work2,nao)
      endif
      call para_allgatherv(work2,num*nao,'D',hmat,idis(1,3),idis(1,4),MPI_COMM_WORLD)
!     call dsymm('L','U',nao,nmo,one,hmat,nao,ortho,nao,zero,work,nao)
!     call dgemm('T','N',nmo,nmo,nao,one,ortho,nao,work,nao,zero,hmat,nao)
      return
end


!----------------------------------
  subroutine orthonorm(hmat,ndim)
!----------------------------------
!
! Gram-Schmidt orthonormalization
!
      use procpar
      implicit none
      integer,intent(in) :: ndim
      integer :: i, j,k
      real(8),intent(inout) :: hmat(ndim,ndim)
      real(8) :: dum, ddot
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
!
      do i= 1,ndim
        dum= zero
        do k= 1,ndim
          dum= dum+hmat(k,i)*hmat(k,i)
        enddo
        dum= one/sqrt(dum)
        do k= 1,ndim
          hmat(k,i)= hmat(k,i)*dum
        enddo
        if(i == ndim) cycle
!$OMP parallel do private(dum)
        do j= i+1,ndim
          if(mod(j,nproc).ne.myrank)cycle
          dum=-ddot(ndim,hmat(1,i),1,hmat(1,j),1)
          call daxpy(ndim,dum,hmat(1,i),1,hmat(1,j),1)
        enddo
        call para_bcast(hmat(1,i+1),ndim,"D",mod(i+1,nproc),MPI_COMM_WORLD)
      enddo
      return
end


!--------------------------------------------------------
  subroutine diag(jobz,uplo,ndim,vector,lda,eigen)
!--------------------------------------------------------
!
! Diagonalize matrix
!
      use iofile,only : iout
!      use procpar, only : master
      use iso_c_binding
      use procpar
      use rokko

      implicit none
      character(1),intent(in) :: jobz, uplo
      integer,intent(in) :: ndim, lda
      real(8),intent(out) :: eigen(lda)
!      real(8),intent(inout) :: vector(*)
      real(8),intent(inout) :: vector(lda, lda)
      real(8)  :: array(lda, lda)

      integer :: i, j, ierr, iproc
!
      type(distributed_matrix) :: mat, Z
      type(grid) :: g
      type(solver) :: solver_
!
      do i = 1, ndim
         do j = 1, i
            vector(i,j) = vector(j,i)
         end do
      end do
!
      print*, "ROKKO ndim=", ndim
      call set_solver(solver_, 'scalapack') ! scalapack') !solver_name)
      print*, "after_solver"
      call set_grid(g, MPI_COMM_WORLD)
      print*, "after_grid"      
      call set_distributed_matrix(mat, ndim, ndim, g, solver_)
      call set_distributed_matrix(Z, ndim, ndim, g, solver_)
      print*, "after_distributed_matrix"
      
      call generate_array_distributed_matrix(vector, mat, &
     &ndim, ndim, lda)
      print*, "after_generate"
!      call print_distributed_matrix(mat)

      call mpi_barrier(mpi_comm_world, ierr)
      call diagonalize(solver_, mat, eigen, Z)
      print*, "after_diagonalize"
     call all_gather(Z, array)
     call mpi_barrier(mpi_comm_world, ierr)

!      array = matmul(transpose(vector), vector)
!  print*, "array=", array                                                                                                                                                    
      do iproc = 0, nproc
          if (iproc == myrank) then
          do i = 1, ndim
           write(*,'(10f8.4)') (array(i, j), j=1, ndim)
          end do
          print*
       endif
     call mpi_barrier(mpi_comm_world, ierr)
     call sleep(1)
     end do

      call del_distributed_matrix(mat)
      call del_distributed_matrix(Z)
      print*, "after_del_distributed"
      call del_solver(solver_)

      return
end


!-------------------------------------------------------------
  subroutine mtrxcanoninv(ortho,overinv,overlap,ndim,newdim)
!-------------------------------------------------------------
!
! Calculate canonicalization matrix and inverse matrix of overlap
!
! The ortho matrix satisfiles (Ortho)-daggar * S * (Ortho) = I
! where S is the overlap matrix.
!
      use procpar, only : master
      use thresh, only : threshover
      implicit none
      integer,intent(in) :: ndim
      integer,intent(out) :: newdim
      integer :: i, j, icount
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(out) :: ortho(ndim,ndim), overinv(ndim,ndim)
      real(8),intent(inout) :: overlap(ndim,ndim)
      real(8),allocatable :: eigen(:), ecanon(:)
!
! Set arrays
!
      call memset(ndim*2)
      allocate(eigen(ndim),ecanon(ndim))
!
! Diagonalize ortho matrix
!
      call diag('V','U',ndim,overlap,ndim,eigen)
!
! Eliminate eigenvectors with small eigenvalues
!
      icount= 0
      do i= 1,ndim
        if(eigen(i) < threshover) icount= icount+1
      enddo
      newdim= ndim-icount
!
! Form canonicalization and inverse matrix
!
!$OMP parallel
!$OMP do
      do i= 1,ndim
        eigen(i)  = one/eigen(i)
      enddo
!$OMP enddo
!$OMP do
      do i= 1,newdim
        ecanon(i)= sqrt(eigen(i+icount))
      enddo
!$OMP enddo
!$OMP do
      do i= 1,ndim
        do j= 1,ndim
          ortho(j,i) =overlap(j,i)*eigen(i)
        enddo
      enddo
!$OMP enddo
!$OMP end parallel
      call dgemm('N','T',ndim,ndim,ndim,one,overlap,ndim,ortho,ndim,zero,overinv,ndim)
!$OMP parallel
!$OMP do
      do i= 1,newdim
        do j= 1,ndim
          ortho(j,i)=overlap(j,i+icount)*ecanon(i)
        enddo
      enddo
!$OMP enddo
!$OMP end parallel
!
! Unset arrays
!
      call memunset(ndim*2)
      deallocate(eigen,ecanon)
      return
end


!--------------------------------------
  function tridot(array1,array2,ndim)
!--------------------------------------
!
! Calculate scf energy from Fock and density matrices.
!
      implicit none
      integer,intent(in) :: ndim
      integer :: i, j, ndim2
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: array1(ndim*(ndim+1)/2), array2(ndim*(ndim+1)/2) 
      real(8) :: tridot, tmp
!
      ndim2= ndim*(ndim+1)/2
      tmp= zero
      j= 0
!
!$OMP parallel do reduction(+:tmp)
      do i= 1,ndim2
        tmp= tmp+array1(i)*array2(i)
      enddo
!$OMP end parallel do
      tmp= tmp+tmp
      do i= 1, ndim
        j= j+i
        tmp= tmp-array1(j)*array2(j)
      enddo
      tridot= tmp
      return
end


!------------------------------------
  subroutine dadd(ndim,array1,array2)
!------------------------------------
!
! array1 = array1 + array2
!
      implicit none
      integer,intent(in) :: ndim
      integer :: i
      real(8),intent(in) :: array2(ndim)
      real(8),intent(inout) :: array1(ndim)
!
      do i= 1,ndim
        array1(i)= array1(i)+array2(i)
      enddo
      return
end


!------------------------------------
  subroutine dsub(ndim,array1,array2)
!------------------------------------
!
! array1 = array1 - array2
!
      implicit none
      integer,intent(in) :: ndim
      integer :: i
      real(8),intent(in) :: array2(ndim)
      real(8),intent(inout) :: array1(ndim)
!
      do i= 1,ndim
        array1(i)= array1(i)-array2(i)
      enddo
      return
end


!----------------------------------------
  subroutine expand(array1,array2,ndim)
!----------------------------------------
!
! Expand upper-triangle matrix to full matrix
!
      implicit none
      integer,intent(in) :: ndim
      integer :: ij, i, j
      real(8),intent(in) :: array1(ndim*(ndim+1)/2)
      real(8),intent(out) :: array2(ndim,ndim)
!
      ij= 0
      do i= 1,ndim
        do j= 1,i
          ij= ij+1
          array2(j,i)= array1(ij)
        enddo
      enddo
      return
end


!----------------------------------------
  subroutine expand2(array1,array2,ndim)
!----------------------------------------
!
! Expand upper-triangle matrix to full matrix
!
      implicit none
      integer,intent(in) :: ndim
      integer :: ij, i, j
      real(8),intent(in) :: array1(ndim*(ndim+1)/2)
      real(8),intent(out) :: array2(ndim,ndim)
!
      ij= 0
      do i= 1,ndim
        do j= 1,i
          ij= ij+1
          array2(j,i)= array1(ij)
          array2(i,j)= array1(ij)
        enddo
      enddo
      return
end


!----------------------------------------
  subroutine icopy(num,ix,incx,iy,incy)
!----------------------------------------
!
! Copy integer array
!
      implicit none
      integer,intent(in) :: num, ix(*), incx, incy
      integer,intent(out) :: iy(*)
      integer :: j
!
      if((incx == 1).and.(incy == 1)) then
        do j= 1,num
          iy(j)= ix(j)
        enddo
      else
        do j= 1,num
          iy(incy*(j-1)+1)=ix(incx*(j-1)+1)
        enddo
      endif
      return
end



!----------------------------------------------------------
  subroutine distarray(idis,nmo,nao,nao3,nocc,nvir,nproc)
!----------------------------------------------------------
!
! Distribute arrays
!
      implicit none
      integer,intent(in) :: nmo, nao, nao3, nocc, nvir, nproc
      integer,intent(out) :: idis(nproc,10)
      integer :: isize1, isize2, isize3, i, istart, iend
!
      isize1= (nmo -1)/nproc+1
      isize2= (nao3-1)/nproc+1
      isize3= (nvir-1)/nproc+1
      do i=1,nproc
        istart=isize1*(i-1)
        if(istart >= nmo) then
          idis(i,1)= 0
          idis(i,2)= 1
          idis(i,3)= 0
          idis(i,4)= 1
        else
          iend  =isize1*i
          if(iend > nmo) iend=nmo
          idis(i,1)= iend-istart
          idis(i,2)= istart
          idis(i,3)=(iend-istart)*nao
          idis(i,4)= istart*nao
        endif
        istart=isize2*(i-1)
        if(istart >= nao3) then
          idis(i,5)= 0
          idis(i,6)= 1
        else
          iend  =isize2*i
          if(iend > nao3) iend=nao3
          idis(i,5)= iend-istart
          idis(i,6)= istart
        endif
        istart=isize3*(i-1)
        if(istart >= nvir) then
          idis(i, 7)= 0
          idis(i, 8)= 1
          idis(i, 9)= 0
          idis(i,10)= 1
        else
          iend  =isize3*i
          if(iend > nvir) iend=nvir
          idis(i, 7)= iend-istart
          idis(i, 8)= istart
          idis(i, 9)=(iend-istart)*nocc
          idis(i,10)= istart*nocc
        endif
      enddo
!
      return
end


!----------------------------------------------------------------------
  subroutine hessianbfgs(ehess,coord,coordold,egrad,egradold,vec,ndim)
!----------------------------------------------------------------------
!
! Update hessian matrix using BFGS method
!
      implicit none
      integer,intent(in) :: ndim
      integer :: i, j, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: coord(ndim), coordold(ndim), egrad(ndim), egradold(ndim)
      real(8),intent(inout) :: ehess(ndim*(ndim+1)/2), vec(ndim,3)
      real(8) :: denom, factor, ddot
!
      do i= 1,ndim
        vec(i,1)= coord(i)-coordold(i)
        vec(i,2)= egrad(i)-egradold(i)
        vec(i,3)= zero
      enddo
      denom= ddot(ndim,vec(1,1),1,vec(1,2),1)
      denom= one/denom
!
      do i= 1,ndim
        ii= i*(i-1)/2
        do j= 1,i-1
          vec(i,3)= vec(i,3)+ehess(ii+j)*vec(j,1)
          vec(j,3)= vec(j,3)+ehess(ii+j)*vec(i,1)
        enddo
        vec(i,3)= vec(i,3)+ehess(ii+i)*vec(i,1)
      enddo
      factor= ddot(ndim,vec(1,1),1,vec(1,3),1)
      factor= one/factor
      do i= 1,ndim
        ii= i*(i-1)/2
        do j= 1,i
          ehess(ii+j)= ehess(ii+j)-factor*vec(j,3)*vec(i,3) &
&                                 +denom*vec(j,2)*vec(i,2)
        enddo
      enddo

!      do i= 1,ndim
!        ii= i*(i-1)/2
!        do j= 1,i-1
!          vec(i,3)= vec(i,3)+ehess(ii+j)*vec(j,2)
!          vec(j,3)= vec(j,3)+ehess(ii+j)*vec(i,2)
!        enddo
!        vec(i,3)= vec(i,3)+ehess(ii+i)*vec(i,2)
!      enddo
!      factor= ddot(ndim,vec(1,2),1,vec(1,3),1)
!      factor= denom+factor*denom*denom
!!
!!$OMP parallel do private(ii)
!      do i= 1,ndim
!        ii= i*(i-1)/2
!        do j= 1,i
!          ehess(ii+j)= ehess(ii+j)+factor*vec(j,1)*vec(i,1) &
!&                                 -denom*(vec(j,1)*vec(i,3)+vec(i,1)*vec(j,3))
!        enddo
!      enddo
!!$OMP end parallel do
      return
end



