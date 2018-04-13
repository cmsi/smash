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
  subroutine maxmemset(datacomp)
!-----------------------
!
! Set maximum memory size
!
      use modjob, only : memory
      use modtype, only : typecomp
      implicit none
      type(typecomp),intent(inout) :: datacomp
      integer :: lock, locm, locg, loct, locb
!
      if(len_trim(memory) /= 0) then
        lock=scan(memory,'K')
        locm=scan(memory,'M')
        locg=scan(memory,'G')
        loct=scan(memory,'T')
        locb=scan(memory,'B')
        if(lock /= 0) then
          read(memory(1:lock-1),*)datacomp%memmax
          datacomp%memmax=datacomp%memmax*125
        elseif(locm /= 0) then
          read(memory(1:locm-1),*)datacomp%memmax
          datacomp%memmax=datacomp%memmax*125000
        elseif(locg /= 0) then
          read(memory(1:locg-1),*)datacomp%memmax
          datacomp%memmax=datacomp%memmax*125000000
        elseif(loct /= 0) then
          read(memory(1:loct-1),*)datacomp%memmax
          datacomp%memmax=datacomp%memmax*125000000*1000
        elseif(locb /= 0) then
          read(memory(1:locb-1),*)datacomp%memmax
          datacomp%memmax=datacomp%memmax/8
        else
          read(memory,*)datacomp%memmax
          datacomp%memmax=datacomp%memmax/8
        endif
      endif
      return
end


!---------------------------
  subroutine memset(msize,datacomp)
!---------------------------
!
! Allocate requested memory size, "msize".
!
      use modparallel, only : master
      use modtype, only : typecomp
      implicit none
      type(typecomp),intent(inout) :: datacomp
      integer,intent(in) :: msize
!
      if(msize < 0) then
        if(master) write(*,'(" Required memory size is negative!",i6,"MB")')msize/125000
        call iabort
      endif
      datacomp%memused= datacomp%memused+msize
      if(datacomp%memused > datacomp%memmax) then
        if(master) then
          write(*,'(" Error! Required memory size exceeds.")')
          write(*,'(" Required:",i6,"MB,  Available:",i6,"MB")') &
&               datacomp%memused/125000, datacomp%memmax/125000
        endif
        call iabort
      endif
      datacomp%memusedmax=max(datacomp%memusedmax,datacomp%memused)
      return
end


!-----------------------------
  subroutine memunset(msize,datacomp)
!-----------------------------
!
! Deallocate requested memory size, "msize".
!
      use modparallel, only : master
      use modtype, only : typecomp
      implicit none
      type(typecomp),intent(inout) :: datacomp
      integer,intent(in) :: msize
!
      datacomp%memused= datacomp%memused-msize
      if(datacomp%memused < 0) then
        datacomp%nwarn= datacomp%nwarn+1
        if(master) write(*,'(" Warning! Msize in memunset is less than 0.")')
      endif
      return
end


!----------------------
  subroutine memcheck(datacomp)
!----------------------
!
! Check memory deallocation
!
      use modparallel, only : master
      use modtype, only : typecomp
      implicit none
      type(typecomp),intent(inout) :: datacomp
!
      if(datacomp%memused /= 0) then
        datacomp%nwarn= datacomp%nwarn+1
        if(master) write(*,'(" Warning! Memory deallocation is not completed.")')
      endif
      return
end


!----------------------------
  subroutine memrest(msize,datacomp)
!----------------------------
!
! Check available memory size
!
      use modtype, only : typecomp
      implicit none
      type(typecomp),intent(inout) :: datacomp
      integer,intent(out) :: msize
!
      msize= datacomp%memmax-datacomp%memused
      return
end

