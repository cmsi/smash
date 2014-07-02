!-----------------------
  subroutine maxmemset
!-----------------------
!
! Set maximum memory size
!
      use modmemory, only : memmax, memory, mem
      implicit none
      integer :: lock, locm, locg, loct, locb
!
      if(len_trim(memory) /= 0) then
        lock=scan(memory,'K')
        locm=scan(memory,'M')
        locg=scan(memory,'G')
        loct=scan(memory,'T')
        locb=scan(memory,'B')
        if(lock /= 0) then
          read(memory(1:lock-1),*)memmax
          memmax=memmax*125
        elseif(locm /= 0) then
          read(memory(1:locm-1),*)memmax
          memmax=memmax*125000
        elseif(locg /= 0) then
          read(memory(1:locg-1),*)memmax
          memmax=memmax*125000000
        elseif(loct /= 0) then
          read(memory(1:loct-1),*)memmax
          memmax=memmax*125000000*1000
        elseif(locb /= 0) then
          read(memory(1:locb-1),*)memmax
          memmax=memmax/8
        else
          read(memory,*)memmax
          memmax=memmax/8
        endif
      endif
      return
end


!---------------------------
  subroutine memset(msize)
!---------------------------
!
! Allocate requested memory size, "msize".
!
      use modmemory, only : memmax, memused, memusedmax
      use modparallel, only : master
      implicit none
      integer,intent(in) :: msize
!
      if(msize < 0) then
          write(*,'(" Required memory size is negative!",i6,"MB")')msize/125000
          call iabort
      endif
      memused= memused+msize
      if(master) then
        if(memused > memmax) then
          write(*,'(" Error! Required memory size exceeds.")')
          write(*,'(" Required:",i6,"MB,  Available:",i6,"MB")')memused/125000, memmax/125000
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
      use modmemory, only : memmax, memused
      use modparallel, only : master
      use modwarn, only : nwarn
      implicit none
      integer,intent(in) :: msize
!
      memused= memused-msize
!     if(master) then
        if(memused < 0) then
          nwarn= nwarn+1
          write(*,'(" Warning! Msize in memunset is less than 0.")')
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
      use modmemory, only : memused
      use modparallel, only : master
      use modwarn, only : nwarn
!
      if(master) then
        if(memused /= 0) then
          nwarn= nwarn+1
          write(*,'(" Warning! Memory deallocation is not completed.")')
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
      use modmemory, only : memmax, memused
      use modparallel, only : master
      use modwarn, only : nwarn
      implicit none
      integer,intent(out) :: msize
!
      msize= memmax-memused
      return
end
