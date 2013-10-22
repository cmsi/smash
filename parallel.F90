!-----------------------
  subroutine para_init
!-----------------------
#ifdef MPI
      use procpar
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
      use procpar
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
      use procpar
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
      use procpar
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
      use procpar
      implicit none
      integer(4) :: icode, ierr
!
      icode=9
      call mpi_abort(MPI_COMM_WORLD,icode,ierr)
#endif
      return
end


!---------------------------------------------------
  subroutine para_bcast(buff,num,vtype,irank,comm)
!---------------------------------------------------
#ifdef MPI
      use procpar
      use iofile, only : iout
      implicit none
      integer :: num, irank
      integer(4) :: num4, irank4, comm, ierr
      real(8) :: buff(*)
      character(len=1) :: vtype
!
      num4= num
      irank4= irank
!
      if(vtype.eq."D" .or. vtype.eq."d") then
         call mpi_bcast(buff,num4,mpi_real8,irank4,comm,ierr)
      elseif(vtype.eq."I" .or. vtype.eq."i") then
         call mpi_bcast(buff,num4,mpi_integer8,irank4,comm,ierr)
      elseif(vtype.eq."C" .or. vtype.eq."c") then
         call mpi_bcast(buff,num4,mpi_character,irank4,comm,ierr)
      else
         if(master) write(iout,'("The 3rd argument of para_bcast is not supported.")')
         call iabort
      endif
#endif
      return
end


!-----------------------------------------------------------------
  subroutine para_allreduce(sbuff,rbuff,num,vtype,op,comm,itype)
!-----------------------------------------------------------------
      use procpar
      use iofile, only : iout
      implicit none
      integer,intent(in) :: num, itype
      integer(4),intent(in) :: op, comm
      real(8),intent(inout) :: sbuff(*)
      real(8),intent(out) :: rbuff(*)
      character(len=1),intent(in) :: vtype
#ifdef MPI
      integer(4) :: num4, mpityp, ierr
!
      if(op.ne.MPI_SUM) then
        if(master) write(iout,'("The 5th argument of para_allreduce is not supported.")')
        call iabort
      endif
!
      if(itype == 1) then
        if(vtype == "D" .or. vtype == "d") then
           call para_allreduced(sbuff,num,op,comm)
        elseif(vtype == "I" .or. vtype =="i") then
           call para_allreducei(sbuff,num,op,comm)
        else
          if(master) write(iout,'("The 4th argument of para_allreduce is not supported.")')
          call iabort
        endif
!
      elseif(itype == 2) then
        num4= num
        if(vtype == "D" .or. vtype == "d") then
          mpityp= mpi_real8
        elseif(vtype == "I" .or. vtype =="i") then
          mpityp= mpi_integer
        else
          if(master) write(iout,'("The 4th argument of para_allreduce is not supported.")')
          call iabort
        endif
        call mpi_allreduce(sbuff,rbuff,num4,mpityp,op,comm,ierr)
!
      else
        if(master) write(iout,'("The 7th argument of para_allreduce is not supported.")')
        call iabort
      endif
!
#else
      if(itype == 2) then
        if(vtype == "D" .or. vtype == "d") then
          call dcopy(num,sbuff,1,rbuff,1)
        elseif(vtype == "I" .or. vtype =="i") then
          call icopy(num,sbuff,1,rbuff,1)
        else
          if(master) write(iout,'("The 4th argument of para_bcast is not supported.")')
          call iabort
        endif
      endif
#endif
      return
end


!------------------------------------------------
  subroutine para_allreduced(sbuff,num,op,comm)
!------------------------------------------------
#ifdef MPI
      use procpar
      implicit none
      integer :: num, icycle, ncycle, nlast, istart, ncopy
      integer(4) :: num4, op, comm, ierr
      real(8) :: sbuff(*), buff(nbuf)
!
      num4= num
      ncycle=(num-1)/nbuf+1
      nlast=mod(num,nbuf)
!
      do icycle=1,ncycle
        istart=1+nbuf*(icycle-1)
        num4=nbuf
        if(icycle.eq.ncycle) num4=nlast
        call mpi_allreduce(sbuff(istart),buff,num4,mpi_real8,op,comm,ierr)
        ncopy=num4
        call dcopy(ncopy,buff,1,sbuff(istart),1)
      enddo
#endif
      return
end


!------------------------------------------------
  subroutine para_allreducei(sbuff,num,op,comm)
!------------------------------------------------
#ifdef MPI
      use procpar
      implicit none
      integer :: num, icycle, ncycle, nlast, istart, ncopy
      integer(4) :: num4, op, comm, ierr
      integer :: sbuff(*), buff(nbuf)
!
      num4= num
      ncycle=(num-1)/nbuf+1
      nlast=mod(num,nbuf)
!
      do icycle=1,ncycle
        istart=1+nbuf*(icycle-1)
        num4=nbuf
        if(icycle.eq.ncycle) num4=nlast
        call mpi_allreduce(sbuff(istart),buff,num4,mpi_integer,op,comm,ierr)
        ncopy=num4
        call icopy(ncopy,buff,1,sbuff(istart),1)
      enddo
#endif
      return
end


!---------------------------------------------------------------------
  subroutine para_allgatherv(sbuff,num,vtype,rbuff,idisa,idisb,comm)
!---------------------------------------------------------------------
      use procpar
!     use procpar, only : nproc, mpi_real8, master
      use iofile, only : iout
      implicit none
      integer,intent(in) :: num, idisa(nproc), idisb(nproc)
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out):: rbuff(*)
      character(len=1),intent(in) :: vtype
      integer(4),intent(in) :: comm
#ifdef MPI
      integer :: i
      integer(4) :: num4, idisa4(nproc), idisb4(nproc), ierr
#endif
!
      if((vtype == 'D').or.(vtype == 'd'))then
#ifdef MPI
        num4= num
        do i= 1,nproc
          idisa4(i)= idisa(i)
          idisb4(i)= idisb(i)
        enddo
        call mpi_allgatherv(sbuff,num4,mpi_real8,rbuff,idisa4,idisb4,mpi_real8,comm,ierr)
#else
        call dcopy(num,sbuff,1,rbuff,1)
#endif
      else
        if(master) write(iout,'(" Error! The 3rd argument of para_allgatherv is not supported.")')
        call iabort
      endif
      return
end


!-------------------------------------------------------------------------------------
  subroutine para_sendrecv(sbuff,nums,vtype,dest,ntags,rbuff,numr,source,ntagr,comm)
!-------------------------------------------------------------------------------------
#ifdef MPI
      use procpar
      use iofile, only : iout
      implicit none
      integer,intent(in) :: nums, dest, ntags, numr, source, ntagr
      integer(4),intent(in) :: comm
      real(8),intent(in) :: sbuff(*)
      real(8),intent(out) :: rbuff(*)
      character(len=1),intent(in) :: vtype
      integer(4) :: nums4, dest4, ntags4, numr4, source4, ntagr4, ierr
      integer(4) :: STATUS(MPI_STATUS_SIZE)
!
      nums4= nums
      dest4= dest
      ntags4= ntags
      numr4= numr
      source4= source
      ntagr4 = ntagr
      if(vtype == "D" .or. vtype == "d") then
        call mpi_sendrecv(sbuff,nums4,MPI_DOUBLE_PRECISION,dest4,ntags4, &
&                         rbuff,numr4,MPI_DOUBLE_PRECISION,source4,ntagr4,comm,STATUS,ierr)
      else
        if(master) write(iout,'(" Error! Only DOUBLE_PRECISION is supported in para_sendrecv.")')
        call abort
      endif
#endif
      return
end
