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
!
! Start MPI and set mpi_comm1=MPI_COMM_WORLD
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer(selected_int_kind(9)) :: ierr
#else
      implicit none
      include "mpif.h"
      integer(selected_int_kind(18)) :: ierr
#endif
      integer,intent(out) :: mpi_comm1
!
      call mpi_init(ierr)
      mpi_comm1= MPI_COMM_WORLD
#endif
      return
end


!--------------------------------------------
  subroutine para_comm_size(nproc,mpi_comm)
!--------------------------------------------
!
! Return the number of processes in mpi_comm
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer(selected_int_kind(9)) :: mpi_comm4, nproc4, ierr
#else
      implicit none
      integer(selected_int_kind(18)) :: mpi_comm4, nproc4, ierr
#endif
      integer,intent(in) :: mpi_comm
      integer,intent(out) :: nproc
!
      mpi_comm4= mpi_comm
      call mpi_comm_size(mpi_comm4,nproc4,ierr)
      nproc= nproc4
#endif
      return
end


!---------------------------------------------
  subroutine para_comm_rank(myrank,mpi_comm)
!---------------------------------------------
!
! Return the MPI rank in mpi_comm
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer(selected_int_kind(9)) :: mpi_comm4, myrank4, ierr
#else
      implicit none
      integer(selected_int_kind(18)) :: mpi_comm4, myrank4, ierr
#endif
      integer,intent(in) :: mpi_comm
      integer,intent(out) :: myrank
!
      mpi_comm4= mpi_comm
      call mpi_comm_rank(mpi_comm4,myrank4,ierr)
      myrank= myrank4
#endif
      return
end


!---------------------------
  subroutine para_finalize
!---------------------------
!
! Finalize MPI
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer(selected_int_kind(9)) :: ierr
#else
      implicit none
      integer(selected_int_kind(18)) :: ierr
#endif
!
      call mpi_finalize(ierr)
#endif
      return
end


!------------------------
  subroutine para_abort
!------------------------
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer(selected_int_kind(9)) :: icode, ierr
#else
      implicit none
      include "mpif.h"
      integer(selected_int_kind(18)) :: icode, ierr
#endif
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
      integer(selected_int_kind(9)),intent(out) :: isize
!
      isize= 4
end


!----------------------------------
  subroutine checkintsize8(isize)
!----------------------------------
      implicit none
      integer(selected_int_kind(18)),intent(out) :: isize
!
      isize= 8
end


!--------------------------------------------------
  subroutine para_bcastr(buff,num,irank,mpi_comm)
!--------------------------------------------------
!
! Broadcast real(8) data
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer(selected_int_kind(9)) :: num4, irank4, mpi_comm4, ierr
#else
      implicit none
      include "mpif.h"
      integer(selected_int_kind(18)) :: num4, irank4, mpi_comm4, ierr
#endif
      integer,intent(in) :: num, irank, mpi_comm
      real(8),intent(inout) :: buff(*)
!
      num4= num
      irank4= irank
      mpi_comm4= mpi_comm
!
      call mpi_bcast(buff,num4,mpi_real8,irank4,mpi_comm4,ierr)
#endif
      return
end


!---------------------------------------------------
  subroutine para_bcasti(ibuff,num,irank,mpi_comm)
!---------------------------------------------------
!
! Broadcast integer data
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      use modparallel, only : checkintsize
      implicit none
      integer(selected_int_kind(9)) :: num4, irank4, mpi_comm4, ierr
#else
      use modparallel, only : checkintsize
      implicit none
      include "mpif.h"
      integer(selected_int_kind(18)) :: num4, irank4, mpi_comm4, ierr
#endif
      integer,intent(in) :: num, irank, mpi_comm
      integer,intent(inout) :: ibuff(*)
      integer :: isize
!
      num4= num
      irank4= irank
      mpi_comm4= mpi_comm
!
      call checkintsize(isize)
      if(isize == 4) then
        call mpi_bcast(ibuff,num4,mpi_integer4,irank4,mpi_comm4,ierr)
      elseif(isize == 8) then
        call mpi_bcast(ibuff,num4,mpi_integer8,irank4,mpi_comm4,ierr)
      endif
#endif
      return
end


!--------------------------------------------------
  subroutine para_bcastc(buff,num,irank,mpi_comm)
!--------------------------------------------------
!
! Broadcast character data
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer(selected_int_kind(9)) :: num4, irank4, mpi_comm4, ierr
#else
      implicit none
      include "mpif.h"
      integer(selected_int_kind(18)) :: num4, irank4, mpi_comm4, ierr
#endif
      integer,intent(in) :: num, irank, mpi_comm
      character(*),intent(inout) :: buff(*)
!
      num4= num
      irank4= irank
      mpi_comm4= mpi_comm
!
      call mpi_bcast(buff,num4,mpi_character,irank4,mpi_comm4,ierr)
#endif
      return
end


!--------------------------------------------------
  subroutine para_bcastl(buff,num,irank,mpi_comm)
!--------------------------------------------------
!
! Broadcast logical data
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer(selected_int_kind(9)) :: num4, irank4, mpi_comm4, ierr, itmp(num)
#else
      implicit none
      include "mpif.h"
      integer(selected_int_kind(18)) :: num4, irank4, mpi_comm4, ierr
      integer(selected_int_kind(9)) :: itmp(num)
#endif
      integer,intent(in) :: num, irank, mpi_comm
      integer :: ii, myrank
      logical,intent(inout) :: buff(*)
!
      call para_comm_rank(myrank,mpi_comm)
!
      if(irank == myrank) then
        do ii= 1,num
          if(buff(ii)) then
            itmp(ii)= 1
          else
            itmp(ii)= 0
          endif
        enddo
      endif
!
      num4= num
      irank4= irank
      mpi_comm4= mpi_comm
!
      call mpi_bcast(itmp,num4,mpi_integer4,irank4,mpi_comm4,ierr)
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
!
! Combine real(8) values from all processes and 
! distributes the result back to all processes in mpi_comm
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer(selected_int_kind(9)) :: num4, mpi_comm4, ierr
#else
      implicit none
      include "mpif.h"
      integer(selected_int_kind(18)) :: num4, mpi_comm4, ierr
#endif
      integer,intent(in) :: num, mpi_comm
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out) :: rbuff(*)
!
      num4= num
      mpi_comm4= mpi_comm
      call mpi_allreduce(sbuff,rbuff,num4,mpi_real8,MPI_SUM,mpi_comm4,ierr)
#else
      implicit none
      integer,intent(in) :: num, mpi_comm
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out) :: rbuff(*)
!
      call dcopy(num,sbuff,1,rbuff,1)
#endif
      return
end


!-------------------------------------------------------
  subroutine para_allreducei(sbuff,rbuff,num,mpi_comm)
!-------------------------------------------------------
!
! Combine integer values from all processes and 
! distributes the result back to all processes in mpi_comm
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      use modparallel, only : checkintsize
      implicit none
      integer(selected_int_kind(9)) :: num4, mpi_comm4, ierr
#else
      use modparallel, only : checkintsize
      implicit none
      include "mpif.h"
      integer(selected_int_kind(18)) :: num4, mpi_comm4, ierr
#endif
      integer,intent(in) :: num, mpi_comm
      integer,intent(in) :: sbuff(*)
      integer,intent(out) :: rbuff(*)
      integer :: isize
!
      num4= num
      mpi_comm4= mpi_comm
!
      call checkintsize(isize)
      if(isize == 4) then
        call mpi_allreduce(sbuff,rbuff,num4,mpi_integer4,MPI_SUM,mpi_comm4,ierr)
      elseif(isize == 8) then
        call mpi_allreduce(sbuff,rbuff,num4,mpi_integer8,MPI_SUM,mpi_comm4,ierr)
      endif
#else
      implicit none
      integer,intent(in) :: num, mpi_comm
      integer,intent(in) :: sbuff(*)
      integer,intent(out) :: rbuff(*)
      integer ii
!
      do ii= 1,num
        rbuff(ii)= sbuff(ii)
      enddo
#endif
      return
end


!--------------------------------------------------------------------------
  subroutine para_allgathervr(sbuff,num,rbuff,idisa,idisb,nproc,mpi_comm)
!--------------------------------------------------------------------------
!
! Gather data from all tasks and deliver the combined data to all tasks in mpi_comm
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer(selected_int_kind(9)) :: num4, idisa4(nproc), idisb4(nproc), mpi_comm4, ierr
#else
      implicit none
      include "mpif.h"
      integer(selected_int_kind(18)) :: num4, idisa4(nproc), idisb4(nproc), mpi_comm4, ierr
#endif
      integer,intent(in) :: num, nproc, idisa(nproc), idisb(nproc), mpi_comm
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out):: rbuff(*)
      integer :: ii
!
      num4= num
      mpi_comm4= mpi_comm
      do ii= 1,nproc
        idisa4(ii)= idisa(ii)
        idisb4(ii)= idisb(ii)
      enddo
      call mpi_allgatherv(sbuff,num4,mpi_real8,rbuff,idisa4,idisb4,mpi_real8,mpi_comm4,ierr)
#else
      implicit none
      integer,intent(in) :: num, nproc, idisa(nproc), idisb(nproc), mpi_comm
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out):: rbuff(*)
!
      call dcopy(num,sbuff,1,rbuff,1)
#endif
!
      return
end


!------------------------------------------------------------------------------------
  subroutine para_sendrecvr(sbuff,nums,idest,ntags,rbuff,numr,isource,ntagr,mpi_comm)
!------------------------------------------------------------------------------------
!
! Send and receive real(8) data
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer(selected_int_kind(9)) :: nums4, idest4, ntags4, numr4, isource4, ntagr4
      integer(selected_int_kind(9)) :: mpi_comm4, ierr, STATUS(MPI_STATUS_SIZE)
#else
      implicit none
      include "mpif.h"
      integer(selected_int_kind(18)) :: nums4, idest4, ntags4, numr4, isource4, ntagr4
      integer(selected_int_kind(18)) :: mpi_comm4, ierr, STATUS(MPI_STATUS_SIZE)
#endif
      integer,intent(in) :: nums, idest, ntags, numr, isource, ntagr, mpi_comm
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out) :: rbuff(*)
!
      nums4= nums
      idest4= idest
      ntags4= ntags
      numr4= numr
      isource4= isource
      ntagr4 = ntagr
      mpi_comm4= mpi_comm
      call mpi_sendrecv(sbuff,nums4,MPI_DOUBLE_PRECISION,idest4,ntags4, &
&                       rbuff,numr4,MPI_DOUBLE_PRECISION,isource4,ntagr4,mpi_comm4,STATUS,ierr)
#else
      implicit none
      integer,intent(in) :: nums, idest, ntags, numr, isource, ntagr, mpi_comm
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out) :: rbuff(*)
!
      call dcopy(nums,sbuff,1,rbuff,1)
#endif
      return
end
