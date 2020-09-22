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
!-------------------------------------
  subroutine tstamp(indext,datacomp)
!-------------------------------------
!
! Print CPU and elapsed times
! indext = 0 : at the beginning of the program
! indext = 1 : at the end of a step
! indext = 2 : at the end of the program
!
      use modtype, only : typecomp
      implicit none
      type(typecomp),intent(inout) :: datacomp
      integer :: indext, iwall2, iwrate, iwmax, iday, ihour, imin
      real(8) :: cpu2, wall0, wall1, sec
      character(len=24) :: tdate
!
      if(.not.datacomp%master) return
      if(indext == 0) then
        call cpu_time(datacomp%cpu0)
        call system_clock(datacomp%iwall0,iwrate,iwmax)
        call fdate(tdate)
        write(datacomp%iout,'(" The job started at ",a)')tdate
        datacomp%cpu1= datacomp%cpu0
        datacomp%iwall1= datacomp%iwall0
!
      elseif(indext == 1) then
        call cpu_time(cpu2)
        call system_clock(iwall2,iwrate,iwmax)
        if(iwall2 < datacomp%iwall1) then
          datacomp%iwall0= datacomp%iwall0-iwmax
          datacomp%iwall1= datacomp%iwall1-iwmax
        endif
        wall0= dble(iwall2-datacomp%iwall0)/dble(iwrate)
        wall1= dble(iwall2-datacomp%iwall1)/dble(iwrate)
        call fdate(tdate)
        write(datacomp%iout,'(1x,"Step CPU :",f10.1,", Total CPU :",f10.1,&
&             " of Master node")') cpu2-datacomp%cpu1, cpu2-datacomp%cpu0
        write(datacomp%iout,'(1x,"Step Wall :",f9.1", Total Wall :",f9.1,&
&             " at ",a24,/)') wall1,wall0,tdate
        datacomp%cpu1 = cpu2
        datacomp%iwall1= iwall2
!
      elseif(indext == 2) then
        call cpu_time(cpu2)
        call system_clock(iwall2,iwrate,iwmax)
        if(iwall2 < datacomp%iwall1) then
          datacomp%iwall0= datacomp%iwall0-iwmax
          datacomp%iwall1= datacomp%iwall1-iwmax
        endif
        wall0= dble(iwall2-datacomp%iwall0)/dble(iwrate)
        call fdate(tdate)
        iday =(iwall2-datacomp%iwall0)/iwrate/86400
        ihour= mod((iwall2-datacomp%iwall0)/iwrate,86400)/3600
        imin = mod((iwall2-datacomp%iwall0)/iwrate,3600)/60
        sec  = dble(iwall2-datacomp%iwall0)/dble(iwrate)-dble(86400*iday+3600*ihour+60*imin)
        write(datacomp%iout,'(1x,"Total CPU time :",f11.1," seconds")') cpu2
        write(datacomp%iout,'(1x,"Total Wall time:",f11.1," seconds")') wall0
        write(datacomp%iout,'(17x,"(",i2," days",i3," hours",i3," minutes",f5.1," seconds)")')&
&                   iday, ihour, imin, sec
        write(datacomp%iout,'(" The job finished at ",a)')tdate
      endif
      return
end


!------------------------------------
  subroutine parallelinfo(datacomp)
!------------------------------------
!
! Write hostname of master node and numbers of processes and threads
!
      use modtype, only : typecomp
!$    use omp_lib
      implicit none
      type(typecomp),intent(inout) :: datacomp
      integer :: nthread, istat, hostnm, llen, len_trim
      character(len=128) :: hostname
!
      nthread=1
!$OMP parallel
!$OMP master
!$    nthread= omp_get_num_threads()
!$OMP end master
!$OMP end parallel
!
      if(datacomp%master) then
        istat= hostnm(hostname)
        llen= len_trim(hostname)
        write(datacomp%iout,'(" Master node is ",a)')hostname(1:llen)
!
        write(datacomp%iout,'(" Number of processes =",i6  )')datacomp%nproc1
        write(datacomp%iout,'(" Number of threads   =",i6,/)')nthread
      endif
      return
end
 

!-----------------------------------
  subroutine opendatfile(datacomp)
!-----------------------------------
!
! Open temporary file
!
      use modtype, only : typecomp
      implicit none
      type(typecomp),intent(inout) :: datacomp
      integer(selected_int_kind(9)) :: getpid, iprocess
      character(len=128) :: filename
!
      if(datacomp%master) then
        iprocess= getpid()
        write(filename,*) iprocess
        filename= "input.dat"//adjustl(filename)
        open(unit=datacomp%inpcopy,file=filename,status='replace')
      endif
!
      return
end


!------------------------------------
  subroutine closedatfile(datacomp)
!------------------------------------
!
! Close temporary file
!
      use modtype, only : typecomp
      implicit none
      type(typecomp),intent(inout) :: datacomp
!
      if(datacomp%master) then
        close(datacomp%inpcopy,status='delete')
      endif
!
      return
end


!---------------------------------------------
  subroutine opencheckfile(datajob,datacomp)
!---------------------------------------------
!
! Open checkpoint file
!
      use modtype, only : typejob, typecomp
      implicit none
      type(typejob),intent(in) :: datajob
      type(typecomp),intent(in) :: datacomp
      integer :: llen
!
      if(datacomp%master) then
        llen= len_trim(datajob%check)
        open(unit=datacomp%icheck,file=datajob%check(1:llen),form='unformatted',status='unknown')
      endif
!
      return
end


!--------------------------------------------
  subroutine closecheckfile(datacomp)
!--------------------------------------------
!
! Close checkpoint file
!
      use modtype, only : typecomp
      implicit none
      type(typecomp),intent(inout) :: datacomp
!
      if(datacomp%master) then
        close(datacomp%icheck)
      endif
!
      return
end


!-------------------------------------------
  subroutine openxyzfile(datajob,datacomp)
!-------------------------------------------
!
! Open xyz file
!
      use modtype, only : typejob, typecomp
      implicit none
      type(typejob),intent(in) :: datajob
      type(typecomp),intent(in) :: datacomp
      integer :: llen
!
      if(datacomp%master) then
        llen= len_trim(datajob%xyz)
        open(unit=datacomp%ixyz,file=datajob%xyz(1:llen),status='replace')
      endif
!
      return
end


!------------------------------------
  subroutine closexyzfile(datacomp)
!------------------------------------
!
! Close xyz file
!
      use modtype, only : typecomp
      implicit none
      type(typecomp),intent(inout) :: datacomp
!
      if(datacomp%master) then
        close(datacomp%ixyz)
      endif
!
      return
end


!-------------------------------------
  subroutine openinputfile(datacomp)
!-------------------------------------
!
! Open argument input file
!
      use modtype, only : typecomp
      implicit none
      type(typecomp),intent(inout) :: datacomp
      integer :: llen
      character(len=256) :: filename
!
      if(datacomp%master) then
        call get_command_argument(1,filename)
        llen= len_trim(filename)
        datacomp%inpstd= 9
        open(unit=datacomp%inpstd,file=filename(1:llen),status='old',err=9999)
      endif
!
      return
9999 write(*,'(" Error! File ",a," is not found.")') filename(1:llen)
     call iabort
end


!--------------------------------------
  subroutine closeinputfile(datacomp)
!--------------------------------------
!
! Close argument input file
!
      use modtype, only : typecomp
      implicit none
      type(typecomp),intent(inout) :: datacomp
!
      if(datacomp%master) then
        close(datacomp%inpstd)
      endif
!
      return
end


!--------------------
  subroutine iabort
!--------------------
!
! Abort the calculation
!
      implicit none
      character(len=24) :: tdate
!
      call fdate(tdate)
      write(*,'(" The job finished at ",a)')tdate
      write(*,'(" ===========================================")')
      write(*,'("   Error! Calculation finished abnormally.")')
      write(*,'(" ===========================================")')
      call para_abort
      call exit
end
