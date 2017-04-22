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


!-----------------------------------------------
  subroutine para_comm_size(nproc,mpi_commin)
!-----------------------------------------------
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
      integer,intent(in) :: mpi_commin
      integer,intent(out) :: nproc
!
      mpi_comm4= mpi_commin
      call mpi_comm_size(mpi_comm4,nproc4,ierr)
      nproc= nproc4
#else
      integer,intent(in) :: mpi_commin
      integer,intent(out) :: nproc
!
      nproc= 1
#endif
      return
end


!-----------------------------------------------
  subroutine para_comm_rank(myrank,mpi_commin)
!-----------------------------------------------
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
      integer,intent(in) :: mpi_commin
      integer,intent(out) :: myrank
!
      mpi_comm4= mpi_commin
      call mpi_comm_rank(mpi_comm4,myrank4,ierr)
      myrank= myrank4
#else
      integer,intent(in) :: mpi_commin
      integer,intent(out) :: myrank
!
      myrank= 0
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


!----------------------------------------------------
  subroutine para_bcastr(buff,num,irank,mpi_commin)
!----------------------------------------------------
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
      integer,intent(in) :: num, irank, mpi_commin
      real(8),intent(inout) :: buff(*)
!
      num4= num
      irank4= irank
      mpi_comm4= mpi_commin
!
      call mpi_bcast(buff,num4,mpi_real8,irank4,mpi_comm4,ierr)
#endif
      return
end


!-----------------------------------------------------
  subroutine para_bcasti(ibuff,num,irank,mpi_commin)
!-----------------------------------------------------
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
      integer,intent(in) :: num, irank, mpi_commin
      integer,intent(inout) :: ibuff(*)
      integer :: isize
!
      num4= num
      irank4= irank
      mpi_comm4= mpi_commin
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


!----------------------------------------------------
  subroutine para_bcastc(buff,num,irank,mpi_commin)
!----------------------------------------------------
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
      integer,intent(in) :: num, irank, mpi_commin
      character(*),intent(inout) :: buff(*)
!
      num4= num
      irank4= irank
      mpi_comm4= mpi_commin
!
      call mpi_bcast(buff,num4,mpi_character,irank4,mpi_comm4,ierr)
#endif
      return
end


!----------------------------------------------------
  subroutine para_bcastl(buff,num,irank,mpi_commin)
!----------------------------------------------------
!
! Broadcast logical data
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer,intent(in) :: num, irank, mpi_commin
      integer(selected_int_kind(9)) :: num4, irank4, mpi_comm4, ierr, itmp(num)
#else
      implicit none
      include "mpif.h"
      integer,intent(in) :: num, irank, mpi_commin
      integer(selected_int_kind(18)) :: num4, irank4, mpi_comm4, ierr
      integer(selected_int_kind(9)) :: itmp(num)
#endif
      integer :: ii, myrank
      logical,intent(inout) :: buff(*)
!
      call para_comm_rank(myrank,mpi_commin)
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
      mpi_comm4= mpi_commin
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


!---------------------------------------------------------
  subroutine para_allreducer(sbuff,rbuff,num,mpi_commin)
!---------------------------------------------------------
!
! Accumulate real(8) values from all processes and 
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
      integer,intent(in) :: num, mpi_commin
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out) :: rbuff(*)
!
      num4= num
      mpi_comm4= mpi_commin
      call mpi_allreduce(sbuff,rbuff,num4,mpi_real8,MPI_SUM,mpi_comm4,ierr)
#else
      implicit none
      integer,intent(in) :: num, mpi_commin
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out) :: rbuff(*)
!
      call dcopy(num,sbuff,1,rbuff,1)
#endif
      return
end


!---------------------------------------------------------
  subroutine para_allreducei(sbuff,rbuff,num,mpi_commin)
!---------------------------------------------------------
!
! Accumulate integer values from all processes and 
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
      integer,intent(in) :: num, mpi_commin
      integer,intent(in) :: sbuff(*)
      integer,intent(out) :: rbuff(*)
      integer :: isize
!
      num4= num
      mpi_comm4= mpi_commin
!
      call checkintsize(isize)
      if(isize == 4) then
        call mpi_allreduce(sbuff,rbuff,num4,mpi_integer4,MPI_SUM,mpi_comm4,ierr)
      elseif(isize == 8) then
        call mpi_allreduce(sbuff,rbuff,num4,mpi_integer8,MPI_SUM,mpi_comm4,ierr)
      endif
#else
      implicit none
      integer,intent(in) :: num, mpi_commin
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


!----------------------------------------------------------------------------
  subroutine para_allgathervr(sbuff,num,rbuff,idisa,idisb,nproc,mpi_commin)
!----------------------------------------------------------------------------
!
! Gather data from all tasks and deliver the combined data to all tasks in mpi_comm
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer,intent(in) :: num, nproc, idisa(nproc), idisb(nproc), mpi_commin
      integer(selected_int_kind(9)) :: num4, idisa4(nproc), idisb4(nproc), mpi_comm4, ierr
#else
      implicit none
      include "mpif.h"
      integer,intent(in) :: num, nproc, idisa(nproc), idisb(nproc), mpi_commin
      integer(selected_int_kind(18)) :: num4, idisa4(nproc), idisb4(nproc), mpi_comm4, ierr
#endif
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out):: rbuff(*)
      integer :: ii
!
      num4= num
      mpi_comm4= mpi_commin
      do ii= 1,nproc
        idisa4(ii)= idisa(ii)
        idisb4(ii)= idisb(ii)
      enddo
      call mpi_allgatherv(sbuff,num4,mpi_real8,rbuff,idisa4,idisb4,mpi_real8,mpi_comm4,ierr)
#else
      implicit none
      integer,intent(in) :: num, nproc, idisa(nproc), idisb(nproc), mpi_commin
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out):: rbuff(*)
!
      call dcopy(num,sbuff,1,rbuff,1)
#endif
!
      return
end


!----------------------------------------------------------------------------------------
  subroutine para_sendrecvr(sbuff,nums,idest,ntags,rbuff,numr,isource,ntagr,mpi_commin)
!----------------------------------------------------------------------------------------
!
! Send and receive real(8) data
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer,intent(in) :: nums, idest, ntags, numr, isource, ntagr, mpi_commin
      integer(selected_int_kind(9)) :: nums4, idest4, ntags4, numr4, isource4, ntagr4
      integer(selected_int_kind(9)) :: mpi_comm4, ierr, STATUS(MPI_STATUS_SIZE)
#else
      implicit none
      include "mpif.h"
      integer,intent(in) :: nums, idest, ntags, numr, isource, ntagr, mpi_commin
      integer(selected_int_kind(18)) :: nums4, idest4, ntags4, numr4, isource4, ntagr4
      integer(selected_int_kind(18)) :: mpi_comm4, ierr, STATUS(MPI_STATUS_SIZE)
#endif
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out) :: rbuff(*)
!
      nums4= nums
      idest4= idest
      ntags4= ntags
      numr4= numr
      isource4= isource
      ntagr4 = ntagr
      mpi_comm4= mpi_commin
      call mpi_sendrecv(sbuff,nums4,MPI_DOUBLE_PRECISION,idest4,ntags4, &
&                       rbuff,numr4,MPI_DOUBLE_PRECISION,isource4,ntagr4,mpi_comm4,STATUS,ierr)
#else
      implicit none
      integer,intent(in) :: nums, idest, ntags, numr, isource, ntagr, mpi_commin
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out) :: rbuff(*)
!
      call dcopy(nums,sbuff,1,rbuff,1)
#endif
      return
end


!--------------------------------------------------------------
  subroutine para_isendr(buff,num,idest,ntag,mpi_commin,ireq)
!--------------------------------------------------------------
!
! Non-blocking MPI_Isend of real(8) data
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer,intent(in) :: num, idest, ntag, mpi_commin
      integer,intent(out) :: ireq
      integer(selected_int_kind(9)) :: num4, idest4, ntag4, ireq4, mpi_comm4, ierr
#else
      implicit none
      include "mpif.h"
      integer,intent(in) :: num, idest, ntag, mpi_commin
      integer,intent(out) :: ireq
      integer(selected_int_kind(18)) :: num4, idest4, ntag4, ireq4, mpi_comm4, ierr
#endif
      real(8),intent(in) :: buff(*)
!
      num4= num
      idest4= idest
      ntag4= ntag
      mpi_comm4= mpi_commin
!
      call mpi_isend(buff,num4,MPI_DOUBLE_PRECISION,idest4,ntag4,mpi_comm4,ireq4,ierr)
!
      ireq= ireq4
#endif
      return
end


!----------------------------------------------------------------
  subroutine para_irecvr(buff,num,isource,ntag,mpi_commin,ireq)
!----------------------------------------------------------------
!
! Non-blocking MPI_Irecv of real(8) data
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer,intent(in) :: num, isource, ntag, mpi_commin
      integer,intent(out) :: ireq
      integer(selected_int_kind(9)) :: num4, isource4, ntag4, ireq4, mpi_comm4, ierr
#else
      implicit none
      include "mpif.h"
      integer,intent(in) :: num, isource, ntag, mpi_commin
      integer,intent(out) :: ireq
      integer(selected_int_kind(18)) :: num4, isource4, ntag4, ireq4, mpi_comm4, ierr
#endif
      real(8),intent(out) :: buff(*)
!
      num4= num
      isource4= isource
      ntag4= ntag
      mpi_comm4= mpi_commin
!
      call mpi_irecv(buff,num4,MPI_DOUBLE_PRECISION,isource4,ntag4,mpi_comm4,ireq4,ierr)
!
      ireq= ireq4
#endif
      return
end


!-------------------------------------
  subroutine para_waitall(nump,ireq)
!-------------------------------------
!
!  MPI_Waitall
!
#ifndef noMPI
#ifndef ILP64
      use mpi
      implicit none
      integer,intent(in) :: nump, ireq(nump)
      integer(selected_int_kind(9)) :: nump4, ireq4(nump), STATUS(MPI_STATUS_SIZE,nump), ierr
#else
      implicit none
      include "mpif.h"
      integer,intent(in) :: nump, ireq(nump)
      integer(selected_int_kind(18)) :: nump4, ireq4(nump), STATUS(MPI_STATUS_SIZE,nump), ierr
#endif
      integer :: ii
!
      nump4= nump
      do ii= 1,nump
        ireq4(ii)= ireq(ii)
      enddo
!
      call mpi_waitall(nump4,ireq4,STATUS,ierr)
#endif
      return
end
