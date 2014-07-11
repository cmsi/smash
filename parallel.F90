!-----------------------
  subroutine para_init
!-----------------------
#ifdef MPI
      use modparallel
      implicit none
      integer(4) :: ierr
!
      call mpi_init(ierr)
#endif
      return
end


!-----------------------------------
  subroutine para_comm_size(np)
!-----------------------------------
#ifdef MPI
      use modparallel
      implicit none
      integer :: np
      integer(4) :: np4, ierr
!
      call mpi_comm_size(MPI_COMM_WORLD,np4,ierr)
      np= np4
#endif
      return
end


!-------------------------------
subroutine para_comm_rank(nid)
!-------------------------------
#ifdef MPI
      use modparallel
      implicit none
      integer :: nid
      integer(4) :: nid4, ierr
!
      call mpi_comm_rank(MPI_COMM_WORLD,nid4,ierr)
      nid= nid4
#endif
      return
end


!---------------------------
  subroutine para_finalize
!---------------------------
#ifdef MPI
      use modparallel
      implicit none
      integer(4) :: ierr
!
      call mpi_finalize(ierr)
#endif
      return
end


!------------------------
  subroutine para_abort
!------------------------
#ifdef MPI
      use modparallel
      implicit none
      integer(4) :: icode, ierr
!
      icode=9
      call mpi_abort(MPI_COMM_WORLD,icode,ierr)
#endif
      return
end


!-----------------------------------------------
  subroutine para_bcastd1(buff,num,irank,comm)
!-----------------------------------------------
#ifdef MPI
      use modparallel
      implicit none
      integer :: num, irank
      integer(4) :: num4, irank4, comm, ierr
      real(8) :: buff(*)
!
      num4= num
      irank4= irank
!
      call mpi_bcast(buff,num4,mpi_real8,irank4,comm,ierr)
#endif
      return
end


!-----------------------------------------------
  subroutine para_bcastd2(buff,num,irank,comm)
!-----------------------------------------------
#ifdef MPI
      use modparallel
      implicit none
      integer :: num, irank
      integer(4) :: num4, irank4, comm, ierr
      real(8) :: buff(*)
!
      num4= num
      irank4= irank
!
      call mpi_bcast(buff,num4,mpi_real8,irank4,comm,ierr)
#endif
      return
end


!------------------------------------------------
  subroutine para_bcasti41(buff,num,irank,comm)
!------------------------------------------------
#ifdef MPI
      use modparallel
      implicit none
      integer :: num, irank
      integer(4) :: num4, irank4, comm, ierr
      integer(4) :: buff(*)
!
      num4= num
      irank4= irank
!
      call mpi_bcast(buff,num4,mpi_integer4,irank4,comm,ierr)
#endif
      return
end


!------------------------------------------------
  subroutine para_bcasti42(buff,num,irank,comm)
!------------------------------------------------
#ifdef MPI
      use modparallel
      implicit none
      integer :: num, irank
      integer(4) :: num4, irank4, comm, ierr
      integer(4) :: buff(*)
!
      num4= num
      irank4= irank
!
      call mpi_bcast(buff,num4,mpi_integer4,irank4,comm,ierr)
#endif
      return
end

!------------------------------------------------
  subroutine para_bcasti81(buff,num,irank,comm)
!------------------------------------------------
#ifdef MPI
      use modparallel
      implicit none
      integer :: num, irank
      integer(4) :: num4, irank4, comm, ierr
      integer(8) :: buff(*)
!
      num4= num
      irank4= irank
!
      call mpi_bcast(buff,num4,mpi_integer8,irank4,comm,ierr)
#endif
      return
end


!------------------------------------------------
  subroutine para_bcasti82(buff,num,irank,comm)
!------------------------------------------------
#ifdef MPI
      use modparallel
      implicit none
      integer :: num, irank
      integer(4) :: num4, irank4, comm, ierr
      integer(8) :: buff(*)
!
      num4= num
      irank4= irank
!
      call mpi_bcast(buff,num4,mpi_integer8,irank4,comm,ierr)
#endif
      return
end


!-----------------------------------------------
  subroutine para_bcastc1(buff,num,irank,comm)
!-----------------------------------------------
#ifdef MPI
      use modparallel
      implicit none
      integer :: num, irank
      integer(4) :: num4, irank4, comm, ierr
      character(*) :: buff(*)
!
      num4= num
      irank4= irank
!
      call mpi_bcast(buff,num4,mpi_character,irank4,comm,ierr)
#endif
      return
end


!-----------------------------------------------
  subroutine para_bcastc2(buff,num,irank,comm)
!-----------------------------------------------
#ifdef MPI
      use modparallel
      implicit none
      integer :: num, irank
      integer(4) :: num4, irank4, comm, ierr
      character(*) :: buff(*)
!
      num4= num
      irank4= irank
!
      call mpi_bcast(buff,num4,mpi_character,irank4,comm,ierr)
#endif
      return
end


!-----------------------------------------------
  subroutine para_bcastl1(buff,num,irank,comm)
!-----------------------------------------------
#ifdef MPI
      use modparallel
      implicit none
      integer :: num, irank, ii
      integer(4) :: num4, irank4, comm, ierr
      logical :: buff(*)
      integer(4) :: itmp(num)
!
      num4= num
      irank4= irank
!
      if(master) then
        do ii= 1,num
          if(buff(ii)) then
            itmp(ii)= 1
          else
            itmp(ii)= 0
          endif
        enddo
      endif
!
      call mpi_bcast(itmp,num4,mpi_integer4,irank4,comm,ierr)
!
      do ii= 1,num
        if(itmp(ii) == 1) then
          buff(ii)= .true.
        else
          buff(ii)= .false.
        endif
      enddo
#endif
      return
end


!-----------------------------------------------
  subroutine para_bcastl2(buff,num,irank,comm)
!-----------------------------------------------
#ifdef MPI
      use modparallel
      implicit none
      integer :: num, irank, ii
      integer(4) :: num4, irank4, comm, ierr
      logical :: buff(*)
      integer(4) :: itmp(num)
!
      num4= num
      irank4= irank
!
      if(master) then
        do ii= 1,num
          if(buff(ii)) then
            itmp(ii)= 1
          else
            itmp(ii)= 0
          endif
        enddo
      endif
!
      call mpi_bcast(itmp,num4,mpi_integer4,irank4,comm,ierr)
!
      do ii= 1,num
        if(itmp(ii) == 1) then
          buff(ii)= .true.
        else
          buff(ii)= .false.
        endif
      enddo
#endif
      return
end


!-------------------------------------------------------
  subroutine para_allreduced1(sbuff,rbuff,num,op,comm)
!-------------------------------------------------------
      use modparallel
      implicit none
      integer,intent(in) :: num
      integer(4),intent(in) :: op, comm
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out) :: rbuff(*)
#ifdef MPI
      integer(4) :: num4, ierr
!
      num4= num
      call mpi_allreduce(sbuff,rbuff,num4,mpi_real8,op,comm,ierr)
#else
      call dcopy(num,sbuff,1,rbuff,1)
#endif
      return
end


!-------------------------------------------------------
  subroutine para_allreduced2(sbuff,rbuff,num,op,comm)
!-------------------------------------------------------
      use modparallel
      implicit none
      integer,intent(in) :: num
      integer(4),intent(in) :: op, comm
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out) :: rbuff(*)
#ifdef MPI
      integer(4) :: num4, ierr
!
      num4= num
      call mpi_allreduce(sbuff,rbuff,num4,mpi_real8,op,comm,ierr)
#else
      call dcopy(num,sbuff,1,rbuff,1)
#endif
      return
end


!-------------------------------------------------------
  subroutine para_allreducei1(sbuff,rbuff,num,op,comm)
!-------------------------------------------------------
      use modparallel
      implicit none
      integer,intent(in) :: num
      integer(4),intent(in) :: op, comm
      integer,intent(in) :: sbuff(*)
      integer,intent(out) :: rbuff(*)
#ifdef MPI
      integer(4) :: num4, ierr
!
      num4= num
      call mpi_allreduce(sbuff,rbuff,num4,mpi_integer8,op,comm,ierr)
#else
      integer ii
!
      do ii= 1,num
        rbuff(ii)= sbuff(ii)
      enddo
#endif
      return
end


!-------------------------------------------------------
  subroutine para_allreducei2(sbuff,rbuff,num,op,comm)
!-------------------------------------------------------
      use modparallel
      implicit none
      integer,intent(in) :: num
      integer(4),intent(in) :: op, comm
      integer,intent(in) :: sbuff(*)
      integer,intent(out) :: rbuff(*)
#ifdef MPI
      integer(4) :: num4, ierr
!
      num4= num
      call mpi_allreduce(sbuff,rbuff,num4,mpi_integer8,op,comm,ierr)
#else
      integer ii
!
      do ii= 1,num
        rbuff(ii)= sbuff(ii)
      enddo
#endif
      return
end


!---------------------------------------------------------------
  subroutine para_allgatherv(sbuff,num,rbuff,idisa,idisb,comm)
!---------------------------------------------------------------
      use modparallel
      implicit none
      integer,intent(in) :: num, idisa(nproc), idisb(nproc)
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out):: rbuff(*)
      integer(4),intent(in) :: comm
#ifdef MPI
      integer :: i
      integer(4) :: num4, idisa4(nproc), idisb4(nproc), ierr
!
      num4= num
      do i= 1,nproc
        idisa4(i)= idisa(i)
        idisb4(i)= idisb(i)
      enddo
      call mpi_allgatherv(sbuff,num4,mpi_real8,rbuff,idisa4,idisb4,mpi_real8,comm,ierr)
#else
      call dcopy(num,sbuff,1,rbuff,1)
#endif
!
      return
end


!---------------------------------------------------------------------------------
  subroutine para_sendrecvd1(sbuff,nums,dest,ntags,rbuff,numr,source,ntagr,comm)
!---------------------------------------------------------------------------------
      use modparallel
      implicit none
      integer,intent(in) :: nums, dest, ntags, numr, source, ntagr
      integer(4),intent(in) :: comm
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out) :: rbuff(*)
#ifdef MPI
      integer(4) :: nums4, dest4, ntags4, numr4, source4, ntagr4, ierr
      integer(4) :: STATUS(MPI_STATUS_SIZE)
!
      nums4= nums
      dest4= dest
      ntags4= ntags
      numr4= numr
      source4= source
      ntagr4 = ntagr
      call mpi_sendrecv(sbuff,nums4,MPI_DOUBLE_PRECISION,dest4,ntags4, &
&                           rbuff,numr4,MPI_DOUBLE_PRECISION,source4,ntagr4,comm,STATUS,ierr)
#else
      call dcopy(nums,sbuff,1,rbuff,1)
#endif
      return
end


!---------------------------------------------------------------------------------
  subroutine para_sendrecvd2(sbuff,nums,dest,ntags,rbuff,numr,source,ntagr,comm)
!---------------------------------------------------------------------------------
      use modparallel
      implicit none
      integer,intent(in) :: nums, dest, ntags, numr, source, ntagr
      integer(4),intent(in) :: comm
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out) :: rbuff(*)
#ifdef MPI
      integer(4) :: nums4, dest4, ntags4, numr4, source4, ntagr4, ierr
      integer(4) :: STATUS(MPI_STATUS_SIZE)
!
      nums4= nums
      dest4= dest
      ntags4= ntags
      numr4= numr
      source4= source
      ntagr4 = ntagr
      call mpi_sendrecv(sbuff,nums4,MPI_DOUBLE_PRECISION,dest4,ntags4, &
&                       rbuff,numr4,MPI_DOUBLE_PRECISION,source4,ntagr4,comm,STATUS,ierr)
#else
      call dcopy(nums,sbuff,1,rbuff,1)
#endif
      return
end



