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
!-----------------------
  subroutine maxmemset
!-----------------------
!
! Set maximum memory size
!
      use modmemory, only : memmax, memory
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
      use modparallel, only : master
      use modmemory, only : memmax, memused, memusedmax
      implicit none
      integer,intent(in) :: msize
!
      if(msize < 0) then
        if(master) write(*,'(" Required memory size is negative!",i6,"MB")')msize/125000
        call iabort
      endif
      memused= memused+msize
      if(memused > memmax) then
        if(master) then
          write(*,'(" Error! Required memory size exceeds.")')
          write(*,'(" Required:",i6,"MB,  Available:",i6,"MB")')memused/125000, memmax/125000
        endif
        call iabort
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
      use modparallel, only : master
      use modmemory, only : memused
      use modwarn, only : nwarn
      implicit none
      integer,intent(in) :: msize
!
      memused= memused-msize
      if(memused < 0) then
        nwarn= nwarn+1
        if(master) write(*,'(" Warning! Msize in memunset is less than 0.")')
      endif
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
      if(memused /= 0) then
        nwarn= nwarn+1
        if(master) write(*,'(" Warning! Memory deallocation is not completed.")')
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
      implicit none
      integer,intent(out) :: msize
!
      msize= memmax-memused
      return
end

