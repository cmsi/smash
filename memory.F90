!---------------------------
  subroutine memset(msize)
!---------------------------
!
! Allocate requested memory size, "msize".
!
      use memory, only : memmax, memused, memusedmax
      use iofile, only : iout
      use procpar, only : master
      implicit none
      integer,intent(in) :: msize
!
      if(msize < 0) then
          write(iout,'(" Required memory size is negative!",i6,"MB")')msize/125000
          call iabort
      endif
      memused= memused+msize
      if(master) then
        if(memused > memmax) then
          write(iout,'(" Error! Required memory size exceeds.")')
          write(iout,'(" Required:",i6,"MB,  Available:",i6,"MB")')memused/125000, memmax/125000
          call iabort
        endif
      endif
      memusedmax=max(memusedmax,memused)
      return
end


!-----------------------------
  subroutine memunset(msize)
!-----------------------------
!
! Deallocate requested memory size, "msize".
!
      use memory, only : memmax, memused
      use procpar, only : master
      use warn, only : nwarn
      use iofile, only : iout
      implicit none
      integer,intent(in) :: msize
!
      memused= memused-msize
!     if(master) then
        if(memused < 0) then
          nwarn= nwarn+1
          write(iout,'(" Warning! Msize in memunset is less than 0.")')
        endif
!     endif
      return
end


!----------------------
  subroutine memcheck
!----------------------
!
! Check memory deallocation
!
      use memory, only : memused
      use procpar, only : master
      use warn, only : nwarn
      use iofile, only : iout
!
      if(master) then
        if(memused /= 0) then
          nwarn= nwarn+1
          write(iout,'(" Warning! Memory deallocation is not completed.")')
        endif
      endif
      return
end


!----------------------------
  subroutine memrest(msize)
!----------------------------
!
! Check available memory size
!
      use memory, only : memmax, memused
      use procpar, only : master
      use warn, only : nwarn
      use iofile, only : iout
      implicit none
      integer,intent(out) :: msize
!
      msize= memmax-memused
      return
end
