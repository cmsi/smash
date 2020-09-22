! Copyright 2014-2020  Kazuya Ishimura
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
!-----------------------------------------
  subroutine maxmemset(datajob,datacomp)
!-----------------------------------------
!
! Set maximum memory size
!
      use modtype, only : typejob, typecomp
      implicit none
      type(typejob),intent(in) :: datajob
      type(typecomp),intent(inout) :: datacomp
      integer :: lock, locm, locg, loct, locb
      integer(selected_int_kind(18)) :: mem64
!
      if(len_trim(datajob%memory) /= 0) then
        lock=scan(datajob%memory,'K')
        locm=scan(datajob%memory,'M')
        locg=scan(datajob%memory,'G')
        loct=scan(datajob%memory,'T')
        locb=scan(datajob%memory,'B')
        if(lock /= 0) then
          read(datajob%memory(1:lock-1),*) datacomp%memmax
          mem64= datacomp%memmax
          mem64= mem64*125
          datacomp%memmax= datacomp%memmax*125
        elseif(locm /= 0) then
          read(datajob%memory(1:locm-1),*) datacomp%memmax
          mem64= datacomp%memmax
          mem64= mem64*125000
          datacomp%memmax= datacomp%memmax*125000
        elseif(locg /= 0) then
          read(datajob%memory(1:locg-1),*) datacomp%memmax
          mem64= datacomp%memmax
          mem64= mem64*125000000
          datacomp%memmax= datacomp%memmax*125000000
        elseif(loct /= 0) then
          read(datajob%memory(1:loct-1),*) datacomp%memmax
          mem64= datacomp%memmax
          mem64= mem64*125000000*1000
          datacomp%memmax= datacomp%memmax*125000000*1000
        elseif(locb /= 0) then
          read(datajob%memory(1:locb-1),*) datacomp%memmax
          mem64= datacomp%memmax
          mem64= mem64/8
          datacomp%memmax= datacomp%memmax/8
        else
          read(datajob%memory,*) datacomp%memmax
          mem64= datacomp%memmax
          mem64= mem64/8
          datacomp%memmax= datacomp%memmax/8
        endif
        if(datacomp%memmax /= mem64) then
          if(datacomp%master) then
            write(*,'(" Error! Compilation with 32-bit integer supports up to 17.1GB memory.",/, &
&                     " Reduce memory size, or compile with 64-bit integer.",/)')
            call iabort
          endif
        endif
      endif
      return
end


!------------------------------------
  subroutine memset(msize,datacomp)
!------------------------------------
!
! Allocate requested memory size, "msize".
!
      use modtype, only : typecomp
      implicit none
      type(typecomp),intent(inout) :: datacomp
      integer,intent(in) :: msize
!
      if(msize < 0) then
        if(datacomp%master) write(*,'(" Required memory size is negative!",i6,"MB")')msize/125000
        call iabort
      endif
      datacomp%memused= datacomp%memused+msize
      if(datacomp%memused > datacomp%memmax) then
        if(datacomp%master) then
          write(*,'(" Error! Required memory size exceeds.")')
          write(*,'(" Required:",i6,"MB,  Available:",i6,"MB")') &
&               datacomp%memused/125000, datacomp%memmax/125000
        endif
        call iabort
      endif
      datacomp%memusedmax=max(datacomp%memusedmax,datacomp%memused)
      return
end


!--------------------------------------
  subroutine memunset(msize,datacomp)
!--------------------------------------
!
! Deallocate requested memory size, "msize".
!
      use modtype, only : typecomp
      implicit none
      type(typecomp),intent(inout) :: datacomp
      integer,intent(in) :: msize
!
      datacomp%memused= datacomp%memused-msize
      if(datacomp%memused < 0) then
        datacomp%nwarn= datacomp%nwarn+1
        if(datacomp%master) write(datacomp%iout,'(" Warning! Msize in memunset is less than 0.")')
      endif
      return
end


!--------------------------------
  subroutine memcheck(datacomp)
!--------------------------------
!
! Check memory deallocation
!
      use modtype, only : typecomp
      implicit none
      type(typecomp),intent(inout) :: datacomp
!
      if(datacomp%memused /= 0) then
        datacomp%nwarn= datacomp%nwarn+1
        if(datacomp%master) write(datacomp%iout,'(" Warning! Memory deallocation is not completed.")')
      endif
      return
end


!-------------------------------------
  subroutine memrest(msize,datacomp)
!-------------------------------------
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

