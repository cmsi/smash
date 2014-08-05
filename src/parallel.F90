! Copyright 2014  Kazuya Ishimura
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
!----------------------------------
  subroutine para_init(mpi_comm1)
!----------------------------------
#ifdef MPI
      use mpi
      implicit none
      integer(4),intent(out) :: mpi_comm1
      integer(4) :: ierr
!
      call mpi_init(ierr)
      mpi_comm1= MPI_COMM_WORLD
#endif
      return
end


!--------------------------------------------
  subroutine para_comm_size(nproc,mpi_comm)
!--------------------------------------------
#ifdef MPI
      use mpi
      implicit none
      integer(4),intent(in) :: mpi_comm
      integer,intent(out) :: nproc
      integer(4) :: nproc4, ierr
!
      call mpi_comm_size(mpi_comm,nproc4,ierr)
      nproc= nproc4
#endif
      return
end


!---------------------------------------------
  subroutine para_comm_rank(myrank,mpi_comm)
!---------------------------------------------
#ifdef MPI
      use mpi
      implicit none
      integer(4),intent(in) :: mpi_comm
      integer,intent(out) :: myrank
      integer(4) :: myrank4, ierr
!
      call mpi_comm_rank(mpi_comm,myrank4,ierr)
      myrank= myrank4
#endif
      return
end


!---------------------------
  subroutine para_finalize
!---------------------------
#ifdef MPI
      use mpi
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
      use mpi
      implicit none
      integer(4) :: icode, ierr
!
      icode=9
      call mpi_abort(MPI_COMM_WORLD,icode,ierr)
#endif
      return
end


!----------------------------------
  subroutine checkintsize4(isize)
!----------------------------------
      implicit none
      integer(4),intent(out) :: isize
!
      isize= 4
end


!----------------------------------
  subroutine checkintsize8(isize)
!----------------------------------
      implicit none
      integer(8),intent(out) :: isize
!
      isize= 8
end


!---------------------------------------------------
  subroutine para_bcastr(buff,num,myrank,mpi_comm)
!---------------------------------------------------
#ifdef MPI
      use mpi
      implicit none
      integer,intent(in) :: num, myrank
      integer(4),intent(in) :: mpi_comm
      integer(4) :: num4, myrank4, ierr
      real(8),intent(inout) :: buff(*)
!
      num4= num
      myrank4= myrank
!
      call mpi_bcast(buff,num4,mpi_real8,myrank4,mpi_comm,ierr)
#endif
      return
end


!----------------------------------------------------
  subroutine para_bcasti(ibuff,num,myrank,mpi_comm)
!----------------------------------------------------
#ifdef MPI
      use mpi
      use modparallel, only : checkintsize
      implicit none
      integer,intent(in) :: num, myrank
      integer(4),intent(in) :: mpi_comm
      integer(4) :: num4, myrank4, ierr
      integer,intent(inout) :: ibuff(*)
      integer :: isize
!
      num4= num
      myrank4= myrank
!
      call checkintsize(isize)
      if(isize == 4) then
        call mpi_bcast(ibuff,num4,mpi_integer4,myrank4,mpi_comm,ierr)
      elseif(isize == 8) then
        call mpi_bcast(ibuff,num4,mpi_integer8,myrank4,mpi_comm,ierr)
      endif
#endif
      return
end


!---------------------------------------------------
  subroutine para_bcastc(buff,num,myrank,mpi_comm)
!---------------------------------------------------
#ifdef MPI
      use mpi
      implicit none
      integer,intent(in) :: num, myrank
      integer(4),intent(in) :: mpi_comm
      integer(4) :: num4, myrank4, ierr
      character(*),intent(inout) :: buff(*)
!
      num4= num
      myrank4= myrank
!
      call mpi_bcast(buff,num4,mpi_character,myrank4,mpi_comm,ierr)
#endif
      return
end


!---------------------------------------------------
  subroutine para_bcastl(buff,num,myrank,mpi_comm)
!---------------------------------------------------
#ifdef MPI
      use mpi
      implicit none
      integer,intent(in) :: num, myrank
      integer(4),intent(in) :: mpi_comm
      integer :: ii
      integer(4) :: num4, myrank4, ierr, itmp(num)
      logical,intent(inout) :: buff(*)
!
      num4= num
      myrank4= myrank
!
      if(myrank == 0) then
        do ii= 1,num
          if(buff(ii)) then
            itmp(ii)= 1
          else
            itmp(ii)= 0
          endif
        enddo
      endif
!
      call mpi_bcast(itmp,num4,mpi_integer4,myrank4,mpi_comm,ierr)
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
  subroutine para_allreducer(sbuff,rbuff,num,mpi_comm)
!-------------------------------------------------------
#ifdef MPI
      use mpi
#endif
      implicit none
      integer,intent(in) :: num
      integer(4),intent(in) :: mpi_comm
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out) :: rbuff(*)
#ifdef MPI
      integer(4) :: num4, ierr
!
      num4= num
      call mpi_allreduce(sbuff,rbuff,num4,mpi_real8,MPI_SUM,mpi_comm,ierr)
#else
      call dcopy(num,sbuff,1,rbuff,1)
#endif
      return
end


!-------------------------------------------------------
  subroutine para_allreducei(sbuff,rbuff,num,mpi_comm)
!-------------------------------------------------------
#ifdef MPI
      use mpi
      use modparallel, only : checkintsize
#endif
      implicit none
      integer,intent(in) :: num
      integer(4),intent(in) :: mpi_comm
      integer,intent(in) :: sbuff(*)
      integer,intent(out) :: rbuff(*)
      integer :: isize
#ifdef MPI
      integer(4) :: num4, ierr
!
      num4= num
!
      call checkintsize(isize)
      if(isize == 4) then
        call mpi_allreduce(sbuff,rbuff,num4,mpi_integer4,MPI_SUM,mpi_comm,ierr)
      elseif(isize == 8) then
        call mpi_allreduce(sbuff,rbuff,num4,mpi_integer8,MPI_SUM,mpi_comm,ierr)
      endif
#else
      integer ii
!
      do ii= 1,num
        rbuff(ii)= sbuff(ii)
      enddo
#endif
      return
end


!--------------------------------------------------------------------
  subroutine para_allgathervr(sbuff,num,rbuff,idisa,idisb,nproc,mpi_comm)
!--------------------------------------------------------------------
#ifdef MPI
      use mpi
#endif
      implicit none
      integer,intent(in) :: num, nproc, idisa(nproc), idisb(nproc)
      integer(4),intent(in) :: mpi_comm
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out):: rbuff(*)
#ifdef MPI
      integer :: i
      integer(4) :: num4, idisa4(nproc), idisb4(nproc), ierr
!
      num4= num
      do i= 1,nproc
        idisa4(i)= idisa(i)
        idisb4(i)= idisb(i)
      enddo
      call mpi_allgatherv(sbuff,num4,mpi_real8,rbuff,idisa4,idisb4,mpi_real8,mpi_comm,ierr)
#else
      call dcopy(num,sbuff,1,rbuff,1)
#endif
!
      return
end


!------------------------------------------------------------------------------------
  subroutine para_sendrecvr(sbuff,nums,dest,ntags,rbuff,numr,source,ntagr,mpi_comm)
!------------------------------------------------------------------------------------
#ifdef MPI
      use mpi
#endif
      implicit none
      integer,intent(in) :: nums, dest, ntags, numr, source, ntagr
      integer(4),intent(in) :: mpi_comm
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
&                       rbuff,numr4,MPI_DOUBLE_PRECISION,source4,ntagr4,mpi_comm,STATUS,ierr)
#else
      call dcopy(nums,sbuff,1,rbuff,1)
#endif
      return
end




