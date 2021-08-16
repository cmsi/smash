! Copyright 2014-2021  Kazuya Ishimura
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
!---------------------------------
  subroutine smashmain(datacomp)
!---------------------------------
!
! This is the main driver of Scalable Molecular Analysis Solver 
! for High performance computing systems (SMASH).
!
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typecomp),intent(inout) :: datacomp
      type(typejob) :: datajob
      type(typemol) :: datamol
      type(typebasis) :: databasis
      integer :: nao, nao2
      real(8),allocatable :: egrad(:), cmo(:), cmob(:), energymo(:), energymob(:)
!
! Write SMASH version, starting time, and parallel information
!
      if(datacomp%master) then
        write(datacomp%iout,&
&           '(" *******************************************",/,&
&             "    Scalable Molecular Analysis Solver for",/,&
&             "      High performance computing systems",/,&
&             "            SMASH Version ",a10/,&
&             "           written by K. ISHIMURA",/,&
&             " *******************************************",/)') datajob%version
      endif
      call tstamp(0,datacomp)
      call parallelinfo(datacomp)
!
! Open input.dat file and read input data
!
      if(datacomp%master) call opendatfile(datacomp)
      call readinput(datajob,datamol,databasis,datacomp)
!
! Open checkpoint file if necessary
!
      if((datacomp%master).and.(datajob%check /= '')) call opencheckfile(datajob,datacomp)
!
! Set parameters and write computaional conditions
!
      call setdetails(datajob,datamol,databasis,datacomp)
!
! Start calculations
!
      nao = databasis%nao
      nao2= nao*nao
      select case(datajob%runtype)
        case('ENERGY')
          select case(datajob%scftype)
            case('RHF')
              call memset(nao2+nao,datacomp)
              allocate(cmo(nao2),energymo(nao))
              call calcrenergy(cmo,energymo,datajob,datamol,databasis,datacomp)
              deallocate(cmo,energymo)
              call memunset(nao2+nao,datacomp)
            case('UHF')
              call memset(nao2*2+nao*2,datacomp)
              allocate(cmo(nao2),cmob(nao2),energymo(nao),energymob(nao))
              call calcuenergy(cmo,cmob,energymo,energymob,datajob,datamol,databasis,datacomp)
              deallocate(cmo,cmob,energymo,energymob)
              call memunset(nao2*2+nao*2,datacomp)
          end select
        case('GRADIENT')
          select case(datajob%scftype)
            case('RHF')
              call memset(datamol%natom*3+nao2+nao,datacomp)
              allocate(egrad(datamol%natom*3),cmo(nao2),energymo(nao))
              call calcrgradient(egrad,cmo,energymo,datajob,datamol,databasis,datacomp)
              deallocate(egrad,cmo,energymo)
              call memunset(datamol%natom*3+nao2+nao,datacomp)
            case('UHF')
              call memset(datamol%natom*3+nao2*2+nao*2,datacomp)
              allocate(egrad(datamol%natom*3),cmo(nao2),cmob(nao2),energymo(nao),energymob(nao))
              call calcugradient(egrad,cmo,cmob,energymo,energymob, &
&                                datajob,datamol,databasis,datacomp)
              deallocate(egrad,cmo,cmob,energymo,energymob)
              call memunset(datamol%natom*3+nao2*2+nao*2,datacomp)
          end select
        case('OPT')
          select case(datajob%scftype)
            case('RHF')
              call memset(datamol%natom*3+nao2+nao,datacomp)
              allocate(egrad(datamol%natom*3),cmo(nao2),energymo(nao))
              call calcrgeometry(egrad,cmo,energymo,datajob,datamol,databasis,datacomp)
              deallocate(egrad,cmo,energymo)
              call memunset(datamol%natom*3+nao2+nao,datacomp)
            case('UHF')
              call memset(datamol%natom*3+nao2*2+nao*2,datacomp)
              allocate(egrad(datamol%natom*3),cmo(nao2),cmob(nao2),energymo(nao),energymob(nao))
              call calcugeometry(egrad,cmo,cmob,energymo,energymob, &
&                                datajob,datamol,databasis,datacomp)
              deallocate(egrad,cmo,cmob,energymo,energymob)
              call memunset(datamol%natom*3+nao2*2+nao*2,datacomp)
          end select
        case default
          if(datacomp%master) then
            write(datacomp%iout,'(" Error! This program does not support runtype= ",a16,".")') datajob%runtype
            call iabort(datacomp)
          endif
      end select
!
! Open, write, and close xyz file if necessary
!
      if(datacomp%master.and.(datajob%xyz /= '')) then
        call openxyzfile(datajob,datacomp)
        call writexyzfile(datamol,datacomp)
        call closexyzfile(datacomp)
      endif
!
! Close input.dat and checkpoint files
!
      if(datacomp%master) call closedatfile(datacomp)
      if(datacomp%master.and.(datajob%check /= '')) call closecheckfile(datacomp)
!
      call memcheck(datacomp)
      call tstamp(2,datacomp)
!
! Write final message
!
      if(datacomp%master) then
        write(datacomp%iout,'(" Used memory :",i8," MB")') datacomp%memusedmax/125000
        if(.not.datacomp%convergedscf) then
          write(datacomp%iout, &
&           '(" =============================================",/, &
&             "    Error! SCF calculation did not converge!",/, &
&             " =============================================")')
          return
        elseif(.not.datacomp%convergedgeom) then
          write(datacomp%iout, &
&           '(" ================================================================",/, &
&             "    Error! Geometry optimization calculation did not converge!",/, &
&             " ================================================================")')
          return
        else
          if(datacomp%nwarn <= 1) then
            write(datacomp%iout, &
&             '(" Your calculation finished successfully with",i3," warning.")') datacomp%nwarn
          else
            write(datacomp%iout, &
&             '(" Your calculation finished successfully with",i3," warnings.")') datacomp%nwarn
          endif
        endif
      endif
      return
end


!-----------------------------------
  subroutine setparallel(datacomp)
!-----------------------------------
!
!  Initialize MPI execution environment
!
      use modtype, only : typecomp
      implicit none
      type(typecomp),intent(inout) :: datacomp
!
! Start MPI parallelization and set mpi_comm1=MPI_COMM_WORLD
!
      call para_init(datacomp%mpi_comm1)
      call para_comm_size(datacomp%nproc1,datacomp%mpi_comm1)
      call para_comm_rank(datacomp%myrank1,datacomp%mpi_comm1)
!
      datacomp%nproc2   = datacomp%nproc1
      datacomp%myrank2  = datacomp%myrank1
      datacomp%mpi_comm2= datacomp%mpi_comm1
!
      if(datacomp%nproc1.gt.1) datacomp%master =(datacomp%myrank1 == 0)
!
      return
end


!------------------------------------------------------------
  subroutine setdetails(datajob,datamol,databasis,datacomp)
!------------------------------------------------------------
!
! Set variables
!
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typemol),intent(inout) :: datamol
      type(typebasis),intent(inout) :: databasis
      type(typecomp),intent(inout) :: datacomp
!
! Read atomic data
!
      call readatom(datajob,datamol,datacomp)
!
! Set basis and ECP functions
!
      call setbasis(datajob,datamol,databasis,datacomp)
      if(datajob%flagecp) call setecp(datamol,databasis,datacomp)
!
! Set maximum memory size
!
      call maxmemset(datajob,datacomp)
!
! Set number of electrons
!
      call setelectron(datajob,datamol,databasis,datacomp)
!
! Reset defaults after reading input file
!
      call setdefault(datajob,datacomp)
!
! Set functional information and adjust the number of DFT grids
!
      call setdft(datajob,datamol,databasis,datacomp)
!
! Set MP2 calculation condition
!
      call setmp2(datajob,datamol,databasis)
!
! Write input data
!
      call writecondition(datajob,datamol,databasis,datacomp)
      call writegeom(datamol,datacomp)
      call writebasis(datajob,datamol,databasis,datacomp)
      if(datajob%flagecp) call writeecp(datajob,datamol,databasis,datacomp)
!
! Set atom charge including dummy atom
!
      call setcharge(datamol,datacomp)
!
! Check keywords
!
      if(datacomp%master) then
        if((datajob%scftype /= 'RHF').and.(datajob%scftype /= 'UHF')) then
          write(datacomp%iout,'(" Error! SCFtype=",a16," is not supported.")') datajob%scftype
          call iabort(datacomp)
        endif
!
        if((datajob%multipole == 'OCTUPOLE').and.(datajob%method == 'MP2')) then
          write(datacomp%iout,'(" Error! This program does not suupport MP2 octupole calculation.")')
          call iabort(datacomp)
        endif
!
        select case(datajob%pop)
          case('MULLIKEN','NONE','')
          case('NPA','NBO')
            if(datajob%method == 'MP2') then
              write(datacomp%iout,'(" Error! This program does not support MP2 ",a3," calculation.")') &
&                   datajob%pop
              call iabort(datacomp)
            endif
          case default
            write(datacomp%iout,'(" Error! This program does not support pop= ", a16,".")') datajob%pop
            call iabort(datacomp)
        end select
!
        select case(datajob%multipole)
          case('DIPOLE','OCTUPOLE','NONE','')
          case default
            write(datacomp%iout,'(" Error! This program does not support multipole= ", a16,".")') &
&                 datajob%multipole
            call iabort(datacomp)
        end select
      endif
!
      return
end


!------------------------------------------
  subroutine setdefault(datajob,datacomp)
!------------------------------------------
!
! Reset defaults after reading input file
!
      use modtype, only : typejob, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typecomp),intent(inout) :: datacomp
      real(8),parameter :: zero= 0.0D+00
!
      select case(datajob%precision)
        case('HIGH')
          if(datajob%cutint2 < zero) datajob%cutint2= 1.0D-12
          if(datajob%dconv   < zero) datajob%dconv  = 5.0D-06
          if(datajob%optconv < zero) datajob%optconv= 1.0D-05
          if(datajob%threshweight < zero) datajob%threshweight=1.0D-08
          if(datajob%threshrho    < zero) datajob%threshrho   =1.0D-06
          if(datajob%threshdfock  < zero) datajob%threshdfock =1.0D-05
          if(datajob%threshdftao  < zero) datajob%threshdftao =1.0D-04
          if(datajob%nrad == 0) datajob%nrad= 150
          if(datajob%nleb == 0) datajob%nleb= 590
        case('MEDIUM')
          if(datajob%cutint2 < zero) datajob%cutint2= 1.0D-11
          if(datajob%dconv   < zero) datajob%dconv  = 5.0D-06
          if(datajob%optconv < zero) datajob%optconv= 1.0D-04
          if(datajob%threshweight < zero) datajob%threshweight=1.0D-08
          if(datajob%threshrho    < zero) datajob%threshrho   =1.0D-05
          if(datajob%threshdfock  < zero) datajob%threshdfock =1.0D-04
          if(datajob%threshdftao  < zero) datajob%threshdftao =1.0D-03
          if(datajob%nrad == 0) datajob%nrad= 96
          if(datajob%nleb == 0) datajob%nleb= 302
        case('LOW')
          if(datajob%cutint2 < zero) datajob%cutint2= 1.0D-10
          if(datajob%dconv   < zero) datajob%dconv  = 1.0D-05
          if(datajob%optconv < zero) datajob%optconv= 1.0D-03
          if(datajob%threshweight < zero) datajob%threshweight=1.0D-08
          if(datajob%threshrho    < zero) datajob%threshrho   =1.0D-04
          if(datajob%threshdfock  < zero) datajob%threshdfock =1.0D-04
          if(datajob%threshdftao  < zero) datajob%threshdftao =1.0D-02
          if(datajob%nrad == 0) datajob%nrad= 72
          if(datajob%nleb == 0) datajob%nleb= 302
        case default
          if(datacomp%master) write(datacomp%iout,'(" Error! This program does not support precision= ", &
&                                                   a16,".")') datajob%precision
          call iabort(datacomp)
      end select
!
      if(datajob%iprint == 0) then
        select case(datajob%print)
          case('MINIMUM')
            datajob%iprint= 1
          case('COMPACT')
            datajob%iprint= 2
          case('NORMAL')
            datajob%iprint= 3
          case('FULLMO')
            datajob%iprint= 4
          case('VERBOSE')
            datajob%iprint= 5
          case('')
            datajob%iprint= 3
          case default
            if(datacomp%master) write(datacomp%iout,'(" Error! This program does not support print= ", &
&                                                      a16,".")') datajob%print
            call iabort(datacomp)
        end select
      endif
!
      return
end


!-------------------------------------------------------------
  subroutine setelectron(datajob,datamol,databasis,datacomp)
!-------------------------------------------------------------
!
! Set number of electrons
!
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typemol),intent(inout) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer :: nume, ii
!
! Calculate total number of electrons
!
      nume= -nint(datamol%charge)
      do ii= 1,datamol%natom
        if(datamol%numatomic(ii) > 0) nume= nume+datamol%numatomic(ii)
      enddo
!
! Subtract electrons of core potentials
!
      if(datajob%flagecp) then
        do ii= 1,datamol%natom
          nume= nume-databasis%izcore(ii)
        enddo
      endif

!
! Calculate numbers of alpha and beta electrons
!
      if((datajob%scftype == 'RHF').and.(datamol%multi /= 1)) then
        if(datacomp%master) write(datacomp%iout,'(" Warning! SCFtype changes from RHF to UHF.")')
        datajob%scftype = 'UHF'
        datacomp%nwarn= datacomp%nwarn+1
      endif
!
      datamol%neleca=(nume+datamol%multi-1)/2
      datamol%nelecb=(nume-datamol%multi+1)/2
      if((datamol%neleca+datamol%nelecb)/= nume) then
        if(datacomp%master) write(datacomp%iout,'(" Error! Spin multiplicity is ",i2, &
&                                 ", but number of elctrons is ",i5,".")') datamol%multi, nume
        call iabort(datacomp)
      endif
!
      return
end


!--------------------------------------------------------------------------
  subroutine calcrenergy(cmo,energymo,datajob,datamol,databasis,datacomp)
!--------------------------------------------------------------------------
!
! Driver of closed-shell energy calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
!$    use omp_lib
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typemol),intent(inout) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer :: nao, nao2, nao3, nshell3, nthread
      real(8),intent(out) :: cmo(databasis%nao**2), energymo(databasis%nao)
      real(8), allocatable :: h1mtrx(:), fock(:), smtrx(:), tmtrx(:), ortho(:), dmtrx(:)
      real(8), allocatable :: xint(:), overinv(:), work(:)
      real(8) :: savedconv, savecutint2, dummy
!
      nao    = databasis%nao
      nao2   = nao*nao
      nao3   =(nao*(nao+1))/2
      nshell3=(databasis%nshell*(databasis%nshell+1))/2
      nthread= 1
!$omp parallel
!$    nthread= omp_get_num_threads()
!$omp end parallel 
!
! Set arrays 1 (actually, fock(nao3*nthread))
!
      call memset(nao3*(4+nthread)+nao2+nshell3,datacomp)
      allocate(h1mtrx(nao3),fock(nao3),smtrx(nao3),tmtrx(nao3),ortho(nao2),dmtrx(nao3), &
&              xint(nshell3))
!
! Calculate nuclear repulsion energy
!
      call nucenergy(datajob%threshatom,datamol,datacomp)
      if(datacomp%master) then
        write(datacomp%iout,'(" Nuclear repulsion energy =",f18.9," Hartree",/)') datamol%enuc
      endif
!
! Set arrays 2
!
      call memset(nao2*2,datacomp)
      allocate(overinv(nao2),work(nao2))
!
! Calculate overlap and 1-electron integrals
!
      call oneei(h1mtrx,smtrx,tmtrx,work,datacomp%nproc2,datacomp%myrank2,datacomp%mpi_comm2, &
&                datajob,datamol,databasis)
!
! Calculate canonical transformation matrix and inverse matrix of overlap
!
      call fullmtrx(smtrx,work,nao)
      call mtrxcanoninv(ortho,overinv,work,nao,datamol%nmo,datajob%threshover,datacomp)
!
! Calculate initial MOs
!
      call guessmo(cmo,cmo,overinv,h1mtrx,ortho,datajob,datamol,databasis,datacomp)
!
! Unset arrays 2
!
      deallocate(overinv,work)
      call memunset(nao2*2,datacomp)
      call tstamp(1,datacomp)
!
! Start SCF
!
      if(datajob%method == 'HF') then
        call calcrhf(h1mtrx,cmo,fock,ortho,smtrx,dmtrx,xint,energymo, &
&                    datajob,datamol,databasis,datacomp)
        if(.not.datacomp%convergedscf) go to 9999
        call writeeigenvalue(energymo,energymo,1,datamol,datacomp)
        if(datacomp%master.and.(mod(datajob%iprint,10) >= 2)) &
&         call writeeigenvector(cmo,energymo,1,datajob,datamol,databasis,datacomp)
        call tstamp(1,datacomp)
!
      elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
        if(datajob%guess == 'HUCKEL') then
          savedconv= datajob%dconv
          savecutint2= datajob%cutint2
          datajob%dconv= min(datajob%dconv*1.0D+04,1.0D-1)
          datajob%cutint2= min(datajob%cutint2*5.0D+02,1.0D-8)
          call calcrhf(h1mtrx,cmo,fock,ortho,smtrx,dmtrx,xint,energymo, &
&                      datajob,datamol,databasis,datacomp)
          if(.not.datacomp%convergedscf) go to 9999
          datajob%dconv= savedconv
          datajob%cutint2= savecutint2
          call tstamp(1,datacomp)
        endif
        call calcrdft(h1mtrx,cmo,fock,ortho,smtrx,dmtrx,xint,energymo, &
&                     datajob,datamol,databasis,datacomp)
        if(.not.datacomp%convergedscf) go to 9999
        call writeeigenvalue(energymo,energymo,1,datamol,datacomp)
        if(datacomp%master.and.(mod(datajob%iprint,10) >= 2)) &
&         call writeeigenvector(cmo,energymo,1,datajob,datamol,databasis,datacomp)
        call tstamp(1,datacomp)
!
      elseif(datajob%method == 'MP2') then
        call calcrhf(h1mtrx,cmo,fock,ortho,smtrx,dmtrx,xint,energymo, &
&                    datajob,datamol,databasis,datacomp)
        if(.not.datacomp%convergedscf) go to 9999
        call writeeigenvalue(energymo,energymo,1,datamol,datacomp)
        if(datacomp%master.and.(mod(datajob%iprint,10) >= 2)) &
&         call writeeigenvector(cmo,energymo,1,datajob,datamol,databasis,datacomp)
        call tstamp(1,datacomp)
        call calcrmp2(cmo,energymo,xint,datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1, &
&                     datajob,datamol,databasis,datacomp)
        call tstamp(1,datacomp)
      else
        if(datacomp%master) then
          write(datacomp%iout,'(" Error! This program does not support method= ",a16,".")') datajob%method
          call iabort(datacomp)
        endif
      endif
!
      if(datajob%method /= 'MP2') then
        if(datacomp%master) then
          write(datacomp%iout,'(" ========================")')
          write(datacomp%iout,'("   Property calculation")')
          write(datacomp%iout,'(" ========================")')
        endif
!
! Calculate Mulliken and/or NBO charge
!
        select case(datajob%pop)
          case('MULLIKEN')
            call calcrmulliken(dmtrx,smtrx,datamol,databasis,datacomp)
          case('NPA','NBO')
            call calcrmulliken(dmtrx,smtrx,datamol,databasis,datacomp)
            call calcrnpa(dmtrx,fock,smtrx,datamol,databasis,datacomp)
        end select
!
! Calculate dipole, quadrupole, and octupole moments
!
        select case(datajob%multipole)
          case('DIPOLE')
            call memset(nao3*6,datacomp)
            allocate(work(nao3*6))
            call calcrdipole(work,work(nao3*3+1),dmtrx,datacomp%nproc1,datacomp%myrank1, &
&                            datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
            deallocate(work)
            call memunset(nao3*6,datacomp)
          case('OCTUPOLE')
            call memset(nao3*29,datacomp)
            allocate(work(nao3*29))
            call calcroctupole(work,work(nao3*3+1),work(nao3*9+1),work(nao3*19+1),dmtrx, &
&                              datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1, &
&                              datajob,datamol,databasis,datacomp)
            deallocate(work)
            call memunset(nao3*29,datacomp)
        end select
      else
        if(datacomp%master) then
          write(datacomp%iout,'(" MP2 property calculations are skipped.",/)')
        endif
      endif
!
9999 continue
!
! Write checkpoint file
!
      if(datacomp%master.and.(datajob%check /= '')) &
&       call writecheck(cmo,dummy,dmtrx,dummy,energymo,dummy,fock,dummy,smtrx,h1mtrx, &
&                       dummy,datajob,datamol,databasis,datacomp)
!
! Unset arrays 1
!
      deallocate(h1mtrx,fock,smtrx,tmtrx,ortho,dmtrx, &
&                xint)
      call memunset(nao3*(4+nthread)+nao2+nshell3,datacomp)
      call tstamp(1,datacomp)
      return
end


!-------------------------------------------------------------------------------------------
  subroutine calcuenergy(cmoa,cmob,energymoa,energymob,datajob,datamol,databasis,datacomp)
!-------------------------------------------------------------------------------------------
!
! Driver of open-shell energy calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
!$    use omp_lib
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typemol),intent(inout) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer :: nao, nao2, nao3, nshell3, nthread
      real(8), intent(out) :: cmoa(databasis%nao**2), cmob(databasis%nao**2)
      real(8), intent(out) :: energymoa(databasis%nao), energymob(databasis%nao)
      real(8), allocatable :: h1mtrx(:), focka(:), fockb(:), smtrx(:), tmtrx(:), ortho(:)
      real(8), allocatable :: dmtrxa(:), dmtrxb(:), xint(:), overinv(:), work(:)
      real(8) :: savedconv, savecutint2, dummy
!
      nao    = databasis%nao
      nao2   = nao*nao
      nao3   =(nao*(nao+1))/2
      nshell3=(databasis%nshell*(databasis%nshell+1))/2
      nthread= 1
!$omp parallel
!$    nthread= omp_get_num_threads()
!$omp end parallel
!
! Set arrays 1 (actually, focka(nao3*nthread), fockb(nao3*nthread))
!
      call memset(nao3*(5+nthread*2)+nao2+nshell3,datacomp)
      allocate(h1mtrx(nao3),focka(nao3),fockb(nao3),smtrx(nao3),tmtrx(nao3), &
&              ortho(nao2),dmtrxa(nao3),dmtrxb(nao3),xint(nshell3))
!
! Calculate nuclear repulsion energy
!
      call nucenergy(datajob%threshatom,datamol,datacomp)
      if(datacomp%master) then
        write(datacomp%iout,'(" Nuclear repulsion energy =",f18.9," Hartree",/)') datamol%enuc
      endif
!
! Set arrays 2
!
      call memset(nao2*2,datacomp)
      allocate(overinv(nao2),work(nao2))
!
! Calculate overlap and 1-electron integrals
!
      call oneei(h1mtrx,smtrx,tmtrx,work,datacomp%nproc2,datacomp%myrank2,datacomp%mpi_comm2, &
&                datajob,datamol,databasis)
!
! Calculate canonical transformation matrix and inverse matrix of overlap
!
      call fullmtrx(smtrx,work,nao)
      call mtrxcanoninv(ortho,overinv,work,nao,datamol%nmo,datajob%threshover,datacomp)
!
! Calculate initial MOs
!
      call guessmo(cmoa,cmob,overinv,h1mtrx,ortho,datajob,datamol,databasis,datacomp)
!
! Unset arrays 2
!
      deallocate(overinv,work)
      call memunset(nao2*2,datacomp)
      call tstamp(1,datacomp)
!
! Start SCF
!
      if(datajob%method == 'HF') then
        call calcuhf(h1mtrx,cmoa,cmob,focka,fockb,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa, &
&                    energymob,datajob,datamol,databasis,datacomp)
        if(.not.datacomp%convergedscf) go to 9999
        call writeeigenvalue(energymoa,energymob,2,datamol,datacomp)
        if(datacomp%master.and.(mod(datajob%iprint,10) >= 2)) then
          call writeeigenvector(cmoa,energymoa,1,datajob,datamol,databasis,datacomp)
          call writeeigenvector(cmob,energymob,2,datajob,datamol,databasis,datacomp)
        endif
        call tstamp(1,datacomp)
      elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
        if(datajob%guess == 'HUCKEL') then
          savedconv= datajob%dconv
          savecutint2= datajob%cutint2
          datajob%dconv= min(datajob%dconv*1.0D+04,1.0D-1)
          datajob%cutint2= min(datajob%cutint2*5.0D+02,1.0D-8)
          call calcuhf(h1mtrx,cmoa,cmob,focka,fockb,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa, &
&                      energymob,datajob,datamol,databasis,datacomp)
          if(.not.datacomp%convergedscf) go to 9999
          datajob%dconv= savedconv
          datajob%cutint2= savecutint2
          call tstamp(1,datacomp)
        endif
        call calcudft(h1mtrx,cmoa,cmob,focka,fockb,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa, &
&                     energymob,datajob,datamol,databasis,datacomp)
        if(.not.datacomp%convergedscf) go to 9999
        call writeeigenvalue(energymoa,energymob,2,datamol,datacomp)
        if(datacomp%master.and.(mod(datajob%iprint,10) >= 2)) then
          call writeeigenvector(cmoa,energymoa,1,datajob,datamol,databasis,datacomp)
          call writeeigenvector(cmob,energymob,2,datajob,datamol,databasis,datacomp)
        endif
        call tstamp(1,datacomp)
!     elseif(method == 'MP2') then
      else
        if(datacomp%master) then
          write(datacomp%iout,'(" Error! This program does not support method= ",a16,".")')datajob%method
          call iabort(datacomp)
        endif
      endif
!
      if(datajob%method /= 'MP2') then
        if(datacomp%master) then
          write(datacomp%iout,'(" ========================")')
          write(datacomp%iout,'("   Property calculation")')
          write(datacomp%iout,'(" ========================")')
        endif
!
! Calculate Mulliken and/or NBO charge
!
        select case(datajob%pop)
          case('MULLIKEN')
            call calcumulliken(dmtrxa,dmtrxb,smtrx,datamol,databasis,datacomp)
          case('NPA','NBO')
            call calcumulliken(dmtrxa,dmtrxb,smtrx,datamol,databasis,datacomp)
            call calcunpa(dmtrxa,dmtrxb,focka,fockb,smtrx,datamol,databasis,datacomp)
        end select
!
! Calculate dipole, quadrupole, and octupole moments
!
        select case(datajob%multipole)
          case('DIPOLE')
            call memset(nao3*6,datacomp)
            allocate(work(nao3*6))
            call calcudipole(work,work(nao3*3+1),dmtrxa,dmtrxb,datacomp%nproc1,datacomp%myrank1, &
&                            datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
            deallocate(work)
            call memunset(nao3*6,datacomp)
          case('OCTUPOLE')
            call memset(nao3*29,datacomp)
            allocate(work(nao3*29))
            call calcuoctupole(work,work(nao3*3+1),work(nao3*9+1),work(nao3*19+1),dmtrxa,dmtrxb, &
&                              datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1, &
&                              datajob,datamol,databasis,datacomp)
            deallocate(work)
            call memunset(nao3*29,datacomp)
        end select
      endif
!
9999 continue
!
! Write checkpoint file
!
      if(datacomp%master.and.(datajob%check /= '')) &
&       call writecheck(cmoa,cmob,dmtrxa,dmtrxb,energymoa,energymob,focka,fockb,smtrx,h1mtrx, &
&                       dummy,datajob,datamol,databasis,datacomp)
!
! Unset arrays 1
!
      deallocate(h1mtrx,focka,fockb,smtrx,tmtrx, &
&                ortho,dmtrxa,dmtrxb,xint)
      call memunset(nao3*(5+nthread*2)+nao2+nshell3,datacomp)
      return
end


!--------------------------------------------------------------------------------
  subroutine calcrgradient(egrad,cmo,energymo,datajob,datamol,databasis,datacomp)
!--------------------------------------------------------------------------------
!
! Driver of energy gradient calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
!$    use omp_lib
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typemol),intent(inout) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer :: nao, nao2, nao3, nshell3, nthread
      real(8), intent(out) :: egrad(datamol%natom*3), cmo(databasis%nao**2), energymo(databasis%nao)
      real(8), allocatable :: h1mtrx(:), fock(:), smtrx(:), tmtrx(:), ortho(:), dmtrx(:)
      real(8), allocatable :: xint(:)
      real(8), allocatable :: overinv(:), work(:)
      real(8) :: egradmax, egradrms
      real(8) :: savedconv, savecutint2, dummy
!
      nao    = databasis%nao
      nao2   = nao*nao
      nao3   =(nao*(nao+1))/2
      nshell3=(databasis%nshell*(databasis%nshell+1))/2
      nthread= 1
!$omp parallel
!$    nthread= omp_get_num_threads()
!$omp end parallel
!
! Set arrays 1 (actually, fock(nao3*nthread))
!
      call memset(nao3*(4+nthread)+nao2+nshell3,datacomp)
      allocate(h1mtrx(nao3),fock(nao3),smtrx(nao3),tmtrx(nao3),ortho(nao2),dmtrx(nao3), &
&              xint(nshell3))
!
! Calculate nuclear repulsion energy
!
      call nucenergy(datajob%threshatom,datamol,datacomp)
      if(datacomp%master) then
        write(datacomp%iout,'(" Nuclear repulsion energy =",f18.9," Hartree",/)') datamol%enuc
      endif
!
! Set arrays 2
!
      call memset(nao2*2,datacomp)
      allocate(overinv(nao2),work(nao2))
!
! Calculate overlap and 1-electron integrals
!
      call oneei(h1mtrx,smtrx,tmtrx,work,datacomp%nproc2,datacomp%myrank2,datacomp%mpi_comm2, &
&                datajob,datamol,databasis)
!
! Calculate canonical transformation matrix and inverse matrix of overlap
!
      call fullmtrx(smtrx,work,nao)
      call mtrxcanoninv(ortho,overinv,work,nao,datamol%nmo,datajob%threshover,datacomp)
!
! Calculate initial MOs
!
      call guessmo(cmo,cmo,overinv,h1mtrx,ortho,datajob,datamol,databasis,datacomp)
!
! Unset arrays 2
!
      deallocate(overinv,work)
      call memunset(nao2*2,datacomp)
      call tstamp(1,datacomp)
!
! Start SCF
!
      if((datajob%method == 'HF').or.(datajob%method == 'MP2')) then
        call calcrhf(h1mtrx,cmo,fock,ortho,smtrx,dmtrx,xint,energymo, &
&                    datajob,datamol,databasis,datacomp)
        if(.not.datacomp%convergedscf) go to 9999
      elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
        if(datajob%guess == 'HUCKEL') then
          savedconv= datajob%dconv
          savecutint2= datajob%cutint2
          datajob%dconv= min(datajob%dconv*1.0D+04,1.0D-1)
          datajob%cutint2= min(datajob%cutint2*5.0D+02,1.0D-8)
          call calcrhf(h1mtrx,cmo,fock,ortho,smtrx,dmtrx,xint,energymo, &
&                      datajob,datamol,databasis,datacomp)
          if(.not.datacomp%convergedscf) go to 9999
          datajob%dconv= savedconv
          datajob%cutint2= savecutint2
          call tstamp(1,datacomp)
        endif
        call calcrdft(h1mtrx,cmo,fock,ortho,smtrx,dmtrx,xint,energymo, &
&                     datajob,datamol,databasis,datacomp)
        if(.not.datacomp%convergedscf) go to 9999
      else
        if(datacomp%master) then
          write(datacomp%iout,'(" Error! This program does not support method= ",a16,".")') datajob%method
          call iabort(datacomp)
        endif
      endif
!
      call writeeigenvalue(energymo,energymo,1,datamol,datacomp)
      if(mod(datajob%iprint,10) >= 2) then
        call writeeigenvector(cmo,energymo,1,datajob,datamol,databasis,datacomp)
      endif
      call tstamp(1,datacomp)
!
! Calculate energy gradient
!
      if(datajob%method == 'HF') then
        call calcgradrhf(cmo,energymo,xint,egrad,datacomp%nproc1,datacomp%myrank1, &
&                        datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
      elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
        call calcgradrdft(cmo,energymo,xint,egrad,datacomp%nproc1,datacomp%myrank1, &
&                         datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
      elseif(datajob%method == 'MP2') then
        call calcgradrmp2(cmo,energymo,xint,egrad,datacomp%nproc1,datacomp%myrank1, &
&                         datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
      else
        if(datacomp%master) then
          write(datacomp%iout,'(" Error! This program does not support method= ",a16,".")') datajob%method
          call iabort(datacomp)
        endif
      endif
!
! Calculate maximum and root mean square gradient values
!
      call calcmaxgrad(egradmax,egradrms,egrad,datamol%natom*3)
      if(datacomp%master) &
&       write(datacomp%iout,'("   Maximum gradient =",f13.8,/,"   RMS gradient     =",f13.8,/, &
&                 " ----------------------------------------------------",/)') egradmax,egradrms
      call tstamp(1,datacomp)
!
      if(datajob%method /= 'MP2') then
        if(datacomp%master) then
          write(datacomp%iout,'(" ========================")')
          write(datacomp%iout,'("   Property calculation")')
          write(datacomp%iout,'(" ========================")')
        endif
!
! Calculate Mulliken and/or NBO charge
!
        select case(datajob%pop)
          case('MULLIKEN')
            call calcrmulliken(dmtrx,smtrx,datamol,databasis,datacomp)
          case('NPA','NBO')
            call calcrmulliken(dmtrx,smtrx,datamol,databasis,datacomp)
            call calcrnpa(dmtrx,fock,smtrx,datamol,databasis,datacomp)
        end select
!
! Calculate dipole, quadrupole, and octupole moments
!
        select case(datajob%multipole)
          case('DIPOLE')
            call memset(nao3*6,datacomp)
            allocate(work(nao3*6))
            call calcrdipole(work,work(nao3*3+1),dmtrx,datacomp%nproc1,datacomp%myrank1, &
&                            datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
            deallocate(work)
            call memunset(nao3*6,datacomp)
          case('OCTUPOLE')
            call memset(nao3*29,datacomp)
            allocate(work(nao3*29))
            call calcroctupole(work,work(nao3*3+1),work(nao3*9+1),work(nao3*19+1),dmtrx, &
&                              datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1, &
&                              datajob,datamol,databasis,datacomp)
            deallocate(work)
            call memunset(nao3*29,datacomp)
        end select
      else
        if(datacomp%master) then
          write(datacomp%iout,'(" MP2 property calculations are skipped.",/)')
        endif
      endif
!
9999 continue
!
! Write checkpoint file
!
      if(datacomp%master.and.(datajob%check /= '')) &
&       call writecheck(cmo,dummy,dmtrx,dummy,energymo,dummy,fock,dummy,smtrx,h1mtrx, &
&                       egrad,datajob,datamol,databasis,datacomp)
!
! Unset arrays 1
!
      deallocate(h1mtrx,fock,smtrx,tmtrx,ortho,dmtrx, &
&                xint)
      call memunset(nao3*(4+nthread)+nao2+nshell3,datacomp)
      call tstamp(1,datacomp)
      return
end


!------------------------------------------------------------------
  subroutine calcugradient(egrad,cmoa,cmob,energymoa,energymob, &
&                          datajob,datamol,databasis,datacomp)
!------------------------------------------------------------------
!
! Driver of open-shell energy gradient calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
!$    use omp_lib
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typemol),intent(inout) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer :: nao, nao2, nao3, nshell3, nthread
      real(8), intent(out) :: egrad(datamol%natom*3), cmoa(databasis%nao**2), cmob(databasis%nao**2)
      real(8), intent(out) :: energymoa(databasis%nao), energymob(databasis%nao)
      real(8), allocatable :: h1mtrx(:), focka(:), fockb(:), smtrx(:), tmtrx(:), ortho(:)
      real(8), allocatable :: dmtrxa(:), dmtrxb(:), xint(:)
      real(8), allocatable :: overinv(:), work(:)
      real(8) :: egradmax, egradrms
      real(8) :: savedconv, savecutint2
!
      nao    = databasis%nao
      nao2   = nao*nao
      nao3   =(nao*(nao+1))/2
      nshell3=(databasis%nshell*(databasis%nshell+1))/2
      nthread= 1
!$omp parallel
!$    nthread= omp_get_num_threads()
!$omp end parallel
!
! Set arrays 1 (actually, focka(nao3*nthread), fockb(nao3*nthread))
!
      call memset(nao3*(5+nthread*2)+nao2+nshell3,datacomp)
      allocate(h1mtrx(nao3),focka(nao3),fockb(nao3),smtrx(nao3),tmtrx(nao3), &
&              ortho(nao2),dmtrxa(nao3),dmtrxb(nao3),xint(nshell3))
!
! Calculate nuclear repulsion energy
!
      call nucenergy(datajob%threshatom,datamol,datacomp)
      if(datacomp%master) then
        write(datacomp%iout,'(" Nuclear repulsion energy =",f18.9," Hartree",/)') datamol%enuc
      endif
!
! Set arrays 2
!
      call memset(nao2*2,datacomp)
      allocate(overinv(nao2),work(nao2))
!
! Calculate overlap and 1-electron integrals
!
      call oneei(h1mtrx,smtrx,tmtrx,work,datacomp%nproc2,datacomp%myrank2,datacomp%mpi_comm2, &
&                datajob,datamol,databasis)
!
! Calculate canonical transformation matrix and inverse matrix of overlap
!
      call fullmtrx(smtrx,work,nao)
      call mtrxcanoninv(ortho,overinv,work,nao,datamol%nmo,datajob%threshover,datacomp)
!
! Calculate initial MOs
!
      call guessmo(cmoa,cmob,overinv,h1mtrx,ortho,datajob,datamol,databasis,datacomp)
!
! Unset arrays 2
!
      deallocate(overinv,work)
      call memunset(nao2*2,datacomp)
      call tstamp(1,datacomp)
!
! Start SCF
!
      if(datajob%method == 'HF') then
        call calcuhf(h1mtrx,cmoa,cmob,focka,fockb,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa, &
&                    energymob,datajob,datamol,databasis,datacomp)
        if(.not.datacomp%convergedscf) go to 9999
      elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
        if(datajob%guess == 'HUCKEL') then
          savedconv= datajob%dconv
          savecutint2= datajob%cutint2
          datajob%dconv= min(datajob%dconv*1.0D+04,1.0D-1)
          datajob%cutint2= min(datajob%cutint2*5.0D+02,1.0D-8)
          call calcuhf(h1mtrx,cmoa,cmob,focka,fockb,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa, &
&                      energymob,datajob,datamol,databasis,datacomp)
          if(.not.datacomp%convergedscf) go to 9999
          datajob%dconv= savedconv
          datajob%cutint2= savecutint2
          call tstamp(1,datacomp)
        endif
        call calcudft(h1mtrx,cmoa,cmob,focka,fockb,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa, &
&                     energymob,datajob,datamol,databasis,datacomp)
        if(.not.datacomp%convergedscf) go to 9999
      else
        if(datacomp%master) then
          write(datacomp%iout,'(" Error! This program does not support method= ",a16,".")') datajob%method
          call iabort(datacomp)
        endif
      endif
!
      call writeeigenvalue(energymoa,energymob,2,datamol,datacomp)
      if(mod(datajob%iprint,10) >= 2) then
        call writeeigenvector(cmoa,energymoa,1,datajob,datamol,databasis,datacomp)
        call writeeigenvector(cmob,energymob,2,datajob,datamol,databasis,datacomp)
      endif
      call tstamp(1,datacomp)
!
! Calculate energy gradient
!
      if(datajob%method == 'HF') then
        call calcgraduhf(cmoa,cmob,energymoa,energymob,xint,egrad,datacomp%nproc1, &
&                        datacomp%myrank1,datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
      elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
        call calcgradudft(cmoa,cmob,energymoa,energymob,xint,egrad,datacomp%nproc1, &
&                         datacomp%myrank1,datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
      else
        if(datacomp%master) then
          write(datacomp%iout,'(" Error! This program does not support method= ",a16, &
&                               " in energy gradient.")') datajob%method
          call iabort(datacomp)
        endif
      endif
!
! Calculate maximum and root mean square gradient values
!
      call calcmaxgrad(egradmax,egradrms,egrad,datamol%natom*3)
      if(datacomp%master) &
&       write(datacomp%iout,'("   Maximum gradient =",f13.8,/,"   RMS gradient     =",f13.8,/, &
&                 " ----------------------------------------------------",/)') egradmax,egradrms
      call tstamp(1,datacomp)
!
      if(datajob%method /= 'MP2') then
        if(datacomp%master) then
          write(datacomp%iout,'(" ========================")')
          write(datacomp%iout,'("   Property calculation")')
          write(datacomp%iout,'(" ========================")')
        endif
!
! Calculate Mulliken and/or NBO charge
!
        select case(datajob%pop)
          case('MULLIKEN')
            call calcumulliken(dmtrxa,dmtrxb,smtrx,datamol,databasis,datacomp)
          case('NPA','NBO')
            call calcumulliken(dmtrxa,dmtrxb,smtrx,datamol,databasis,datacomp)
            call calcunpa(dmtrxa,dmtrxb,focka,fockb,smtrx,datamol,databasis,datacomp)
        end select
!
! Calculate dipole, quadrupole, and octupole moments
!
        select case(datajob%multipole)
          case('DIPOLE')
            call memset(nao3*6,datacomp)
            allocate(work(nao3*6))
            call calcudipole(work,work(nao3*3+1),dmtrxa,dmtrxb,datacomp%nproc1,datacomp%myrank1, &
&                            datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
            deallocate(work)
            call memunset(nao3*6,datacomp)
          case('OCTUPOLE')
            call memset(nao3*29,datacomp)
            allocate(work(nao3*29))
            call calcuoctupole(work,work(nao3*3+1),work(nao3*9+1),work(nao3*19+1),dmtrxa,dmtrxb, &
&                              datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1, &
&                              datajob,datamol,databasis,datacomp)
            deallocate(work)
            call memunset(nao3*29,datacomp)
        end select
      endif
!
9999 continue
!
! Write checkpoint file
!
      if(datacomp%master.and.(datajob%check /= '')) &
&       call writecheck(cmoa,cmob,dmtrxa,dmtrxb,energymoa,energymob,focka,fockb,smtrx,h1mtrx, &
&                       egrad,datajob,datamol,databasis,datacomp)
!
! Unset arrays 1
!
      deallocate(h1mtrx,focka,fockb,smtrx,tmtrx, &
&                ortho,dmtrxa,dmtrxb,xint)
      call memunset(nao3*(5+nthread*2)+nao2+nshell3,datacomp)
      call tstamp(1,datacomp)
      return
end


!----------------------------------------------------------------------------------
  subroutine calcrgeometry(egrad,cmo,energymo,datajob,datamol,databasis,datacomp)
!----------------------------------------------------------------------------------
!
! Driver of geometry optimization calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
!$    use omp_lib
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typemol),intent(inout) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer,allocatable :: iredun(:)
      integer :: nao, nao2, nao3, nshell3, natom3, nthread, ii, iopt
      integer :: isizered, numbond, numangle, numtorsion, numredun, maxredun
      real(8), parameter :: third=0.3333333333333333D+00
      real(8), intent(out) :: egrad(datamol%natom*3), cmo(databasis%nao**2), energymo(databasis%nao)
      real(8), allocatable :: h1mtrx(:), fock(:), smtrx(:), tmtrx(:), ortho(:), dmtrx(:)
      real(8), allocatable :: xint(:)
      real(8), allocatable :: egradold(:), ehess(:)
      real(8), allocatable :: overinv(:), work(:,:)
      real(8), allocatable :: workv(:), coordredun(:), egradredun(:)
      real(8) :: egradmax, egradrms
      real(8) :: savedconv, savecutint2, dummy
      logical :: exceed
!
      nao    = databasis%nao
      nao2   = nao*nao
      nao3   =(nao*(nao+1))/2
      nshell3=(databasis%nshell*(databasis%nshell+1))/2
      natom3 = datamol%natom*3
      nthread= 1
!$omp parallel
!$    nthread= omp_get_num_threads()
!$omp end parallel
!
      datacomp%convergedgeom=.false.
!
! Calculate redundant coordinate
!
      if(.not.datajob%cartesian) then
        isizered= datamol%natom*4*10
        call memset(isizered,datacomp)
        allocate(iredun(isizered))
        do ii= 1,10
          call setredundantcoord(iredun,isizered,numbond,numangle,numtorsion,exceed, &
&                                datajob,datamol,datacomp)
          if(.not.exceed) exit
          call memunset(isizered,datacomp)
          deallocate(iredun)
          isizered= isizered*2
          call memset(isizered,datacomp)
          allocate(iredun(isizered))
          if(ii == 10) then
            write(datacomp%iout,'(" Error! The array size for redundant coordinate is too large.")')
            call iabort(datacomp)
          endif
        enddo
        numredun= numbond+numangle+numtorsion
        maxredun= max(numredun,natom3)
      endif
!
! Set arrays for energy (actually, fock(nao3*nthread))
!
      call memset(nao3*(4+nthread)+nao2+nshell3,datacomp)
      allocate(h1mtrx(nao3),fock(nao3),smtrx(nao3),tmtrx(nao3),ortho(nao2),dmtrx(nao3), &
&              xint(nshell3))
!
! Set arrays for energy gradient and geometry optimization
!
      if(datajob%cartesian) then
        call memset(natom3+natom3*(natom3+1)/2,datacomp)
        allocate(egradold(natom3),ehess(natom3*(natom3+1)/2))
      else
        call memset(numredun*4+numredun*(numredun+1)/2,datacomp)
        allocate(coordredun(numredun*2),egradredun(numredun*2),ehess(numredun*(numredun+1)/2))
      endif
!
! Start geometry optimization cycle
!
      do iopt= 1,datajob%nopt
!
! Calculate nuclear repulsion energy
!
        call nucenergy(datajob%threshatom,datamol,datacomp)
        if(datacomp%master) then
          write(datacomp%iout,'(" Nuclear repulsion energy =",f18.9," Hartree",/)') datamol%enuc
        endif
!
! Set work arrays 1
!
        call memset(nao2*2,datacomp)
        allocate(overinv(nao2),work(nao,nao))
!
! Calculate overlap and 1-electron integrals
!
        call oneei(h1mtrx,smtrx,tmtrx,work,datacomp%nproc2,datacomp%myrank2,datacomp%mpi_comm2, &
&                  datajob,datamol,databasis)
!
! Calculate canonical transformation matrix and inverse matrix of overlap
!
        call fullmtrx(smtrx,work,nao)
        call mtrxcanoninv(ortho,overinv,work,nao,datamol%nmo,datajob%threshover,datacomp)
!
! Calculate initial MOs
!
        if(iopt == 1) then
          call guessmo(cmo,cmo,overinv,h1mtrx,ortho,datajob,datamol,databasis,datacomp)
        endif
!
! Unset work arrays 1
!
        deallocate(overinv,work)
        call memunset(nao2*2,datacomp)
        call tstamp(1,datacomp)
!
! Calculate energy
!
        if((datajob%method == 'HF').or.(datajob%method == 'MP2')) then
          call calcrhf(h1mtrx,cmo,fock,ortho,smtrx,dmtrx,xint,energymo, &
&                      datajob,datamol,databasis,datacomp)
          if(.not.datacomp%convergedscf) go to 9999
        elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
          if((iopt == 1).and.(datajob%guess == 'HUCKEL')) then
            savedconv= datajob%dconv
            savecutint2= datajob%cutint2
            datajob%dconv= min(datajob%dconv*1.0D+04,1.0D-1)
            datajob%cutint2= min(datajob%cutint2*5.0D+02,1.0D-8)
            call calcrhf(h1mtrx,cmo,fock,ortho,smtrx,dmtrx,xint,energymo, &
&                        datajob,datamol,databasis,datacomp)
            if(.not.datacomp%convergedscf) go to 9999
            datajob%dconv= savedconv
            datajob%cutint2= savecutint2
            call tstamp(1,datacomp)
          endif
          call calcrdft(h1mtrx,cmo,fock,ortho,smtrx,dmtrx,xint,energymo, &
&                       datajob,datamol,databasis,datacomp)
          if(.not.datacomp%convergedscf) go to 9999
        else
          if(datacomp%master) then
            write(datacomp%iout,'(" Error! This program does not support method= ",a16,".")') datajob%method
            call iabort(datacomp)
          endif
        endif
!
        if(((iopt == 1).and.(mod(datajob%iprint,10) >= 2)).or. &
&          (mod(datajob%iprint,10) >= 3)) then
          call writeeigenvalue(energymo,energymo,1,datamol,datacomp)
        endif
        if(((iopt == 1).and.(mod(datajob%iprint,10) >= 2)).or. &
&          (mod(datajob%iprint,10) >= 5)) then
          call writeeigenvector(cmo,energymo,1,datajob,datamol,databasis,datacomp)
        endif
        call tstamp(1,datacomp)
!
! Calculate energy gradient
!
        if(datajob%method == 'HF') then
          call calcgradrhf(cmo,energymo,xint,egrad,datacomp%nproc1,datacomp%myrank1, &
&                          datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
        elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
          call calcgradrdft(cmo,energymo,xint,egrad,datacomp%nproc1,datacomp%myrank1, &
&                           datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
        elseif(datajob%method == 'MP2') then
          call calcgradrmp2(cmo,energymo,xint,egrad,datacomp%nproc1,datacomp%myrank1, &
&                           datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
        else
          if(datacomp%master) then
            write(datacomp%iout,'(" Error! This program does not support method= ",a16,".")') datajob%method
            call iabort(datacomp)
          endif
        endif
!
! Calculate maximum and root mean square gradient values
!
        call calcmaxgrad(egradmax,egradrms,egrad,natom3)
        if(datacomp%master) &
&         write(datacomp%iout, &
&                 '(" -----------------------------------",/, &
&                   "   Geometry convergence check",/, &
&                   "     Optimization cycle:",i4," /",i4,/, &
&                   "     Optconv           =",f10.6,/, &
&                   " -----------------------------------",/, &
&                   "     Maximum gradient  =",f10.6,/, &
&                   "     RMS gradient      =",f10.6,/, &
&                   " -----------------------------------")') &
&                   iopt,datajob%nopt,datajob%optconv,egradmax,egradrms
!
! Check convergence
!
        if((egradmax <= datajob%optconv).and.(egradrms <= datajob%optconv*third)) then
          if(datacomp%master) &
&           write(datacomp%iout,'("   ============================",/, &
&                                 "   ==== Geometry converged ====",/, &
&                                 "   ============================",/)')
          datacomp%convergedgeom=.true.
          call tstamp(1,datacomp)
          exit
        endif
        call tstamp(1,datacomp)
!
! Write checkpoint file
!
        if(datacomp%master.and.(datajob%check /= '')) &
&         call writecheck(cmo,dummy,dmtrx,dummy,energymo,dummy,fock,dummy,smtrx,h1mtrx, &
&                         egrad,datajob,datamol,databasis,datacomp)
!
! Set work arrays 2
!
        if(datajob%cartesian) then
          call memset(natom3*3,datacomp)
          allocate(workv(natom3*3))
        else
          call memset(maxredun*maxredun*4+maxredun*3,datacomp)
          allocate(work(maxredun*maxredun,4),workv(maxredun*3))
        endif
!
! Calculate new coordinate
!
        if(datajob%cartesian) then
          call calcnewcoord(datamol%coord,datamol%coordold,egrad,egradold,ehess,workv,natom3,iopt, &
&                           datajob,datamol,datacomp)
        else
          call calcnewcoordred(datamol%coord,datamol%coordold,coordredun,egrad,egradredun,ehess, &
&                              work(1,1),work(1,2),work(1,3),work(1,4),workv,iopt,iredun,isizered, &
&                              maxredun,numbond,numangle,numtorsion,numredun, &
&                              datajob,datamol,datacomp)
        endif
!
! Unset work arrays 2
!
        if(datajob%cartesian) then
          deallocate(workv)
          call memunset(natom3*3,datacomp)
        else
          deallocate(work,workv)
          call memunset(maxredun*maxredun*4+maxredun*3,datacomp)
        endif
!!
!! Set guess MO calculation flag from Huckel to projection
!!
!       call setnextopt(coordold,datamol%natom,iopt,datajob)
!
        if(iopt == datajob%nopt) then
          if(datacomp%master) then
            write(datacomp%iout,'(" ==============================")')
            write(datacomp%iout,'("   Geometry did not converge.")')
            write(datacomp%iout,'(" =============================="/)')
            call tstamp(1,datacomp)
          endif
          go to 9999
        endif
!
! Print geometry
!
         call writegeom(datamol,datacomp)
        call tstamp(1,datacomp)
      enddo
!
! End of optimization cycle 
!
!
! Print MO energies and coefficients if necessary
!
      if((iopt /= 1).and.(mod(datajob%iprint,10) >= 2)) then
        call writeeigenvalue(energymo,energymo,1,datamol,datacomp)
      endif
      if(((iopt /= 1).and.(mod(datajob%iprint,10) >= 2)).and.(mod(datajob%iprint,10) <= 4)) then
        call writeeigenvector(cmo,energymo,1,datajob,datamol,databasis,datacomp)
      endif
!
      if(datajob%method /= 'MP2') then
        if(datacomp%master) then
          write(datacomp%iout,'(" ========================")')
          write(datacomp%iout,'("   Property calculation")')
          write(datacomp%iout,'(" ========================")')
        endif
!
! Calculate Mulliken and/or NBO charge
!
        select case(datajob%pop)
          case('MULLIKEN')
            call calcrmulliken(dmtrx,smtrx,datamol,databasis,datacomp)
          case('NPA','NBO')
            call calcrmulliken(dmtrx,smtrx,datamol,databasis,datacomp)
            call calcrnpa(dmtrx,fock,smtrx,datamol,databasis,datacomp)
        end select
!
! Calculate dipole, quadrupole, and octupole moments
!
        select case(datajob%multipole)
          case('DIPOLE')
            call memset(nao3*6,datacomp)
            allocate(work(nao3,6))
            call calcrdipole(work,work(1,4),dmtrx,datacomp%nproc1,datacomp%myrank1, &
&                            datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
            deallocate(work)
            call memunset(nao3*6,datacomp)
          case('OCTUPOLE')
            call memset(nao3*29,datacomp)
            allocate(work(nao3,29))
            call calcroctupole(work,work(1,4),work(1,10),work(1,20),dmtrx, &
&                              datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1, &
&                              datajob,datamol,databasis,datacomp)
            deallocate(work)
            call memunset(nao3*29,datacomp)
        end select
      else
        if(datacomp%master) then
          write(datacomp%iout,'(" MP2 property calculations are skipped.",/)')
        endif
      endif
!
! Write optimized geometry
!
      if(datacomp%master) then
        write(datacomp%iout,'(" ===========")')
        write(datacomp%iout,'("   Summary")')
        write(datacomp%iout,'(" ===========")')
        write(datacomp%iout,'(" ----------------------------------------------------")')
        if(datajob%method == 'HF') then
          write(datacomp%iout,'("   Final HF Energy =",f18.9," Hartree")') datamol%escf
        elseif(datajob%method == 'MP2') then
          write(datacomp%iout,'("   Final MP2 Energy =",f18.9," Hartree")')datamol%escf+datamol%emp2
        elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
          write(datacomp%iout,'("   Final DFT Energy =",f18.9," Hartree")')datamol%escf
        endif
        call writeoptgeom(datamol,datacomp)
      endif
!
9999 continue
!
! Write checkpoint file
!
      if(datacomp%master.and.(datajob%check /= '')) &
&       call writecheck(cmo,dummy,dmtrx,dummy,energymo,dummy,fock,dummy,smtrx,h1mtrx, &
&                       egrad,datajob,datamol,databasis,datacomp)
!
! Unset arrays for energy gradient and geometry optimization
!
      if(datajob%cartesian) then
        deallocate(egradold,ehess)
        call memunset(natom3+natom3*(natom3+1)/2,datacomp)
      else
        deallocate(coordredun,egradredun,ehess)
        call memunset(numredun*4+numredun*(numredun+1)/2,datacomp)
      endif
!
! Unset arrays for energy
!
      deallocate(h1mtrx,fock,smtrx,tmtrx,ortho,dmtrx, &
&                xint)
      call memunset(nao3*(4+nthread)+nao2+nshell3,datacomp)
!
! Unset array for redundant coordinate
!
      if(.not.datajob%cartesian) then
        deallocate(iredun)
        call memunset(isizered,datacomp)
      endif
!
      call tstamp(1,datacomp)
      return
end


!------------------------------------------------------------------
  subroutine calcugeometry(egrad,cmoa,cmob,energymoa,energymob, &
&                          datajob,datamol,databasis,datacomp)
!------------------------------------------------------------------
!
! Driver of open-shell geometry optimization calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
!$    use omp_lib
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typemol),intent(inout) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer,allocatable :: iredun(:)
      integer :: nao, nao2, nao3, nshell3, natom3, nthread, ii, iopt
      integer :: isizered, numbond, numangle, numtorsion, numredun, maxredun
      real(8), parameter :: third=0.3333333333333333D+00
      real(8), intent(out) :: egrad(datamol%natom*3), cmoa(databasis%nao**2), cmob(databasis%nao**2)
      real(8), intent(out) :: energymoa(databasis%nao), energymob(databasis%nao)
      real(8), allocatable :: h1mtrx(:), focka(:), fockb(:), smtrx(:), tmtrx(:), ortho(:)
      real(8), allocatable :: dmtrxa(:), dmtrxb(:), xint(:)
      real(8), allocatable :: egradold(:), ehess(:)
      real(8), allocatable :: overinv(:,:), work(:,:)
      real(8), allocatable :: workv(:), coordredun(:), egradredun(:)
      real(8) :: egradmax, egradrms
      real(8) :: savedconv, savecutint2
      logical :: exceed
!
      nao    = databasis%nao
      nao2   = nao*nao
      nao3   =(nao*(nao+1))/2
      nshell3=(databasis%nshell*(databasis%nshell+1))/2
      natom3 = datamol%natom*3
      nthread= 1
!$omp parallel
!$    nthread= omp_get_num_threads()
!$omp end parallel
!
      datacomp%convergedgeom=.false.
!
! Calculate redundant coordinate
!
      if(.not.datajob%cartesian) then
        isizered= datamol%natom*4*10
        call memset(isizered,datacomp)
        allocate(iredun(isizered))
        do ii= 1,10
          call setredundantcoord(iredun,isizered,numbond,numangle,numtorsion,exceed, &
&                                datajob,datamol,datacomp)
          if(.not.exceed) exit
          call memunset(isizered,datacomp)
          deallocate(iredun)
          isizered= isizered*2
          call memset(isizered,datacomp)
          allocate(iredun(isizered))
          if(ii == 10) then
            write(datacomp%iout,'(" Error! The array size for redundant coordinate is too large.")')
            call iabort(datacomp)
          endif
        enddo
        numredun= numbond+numangle+numtorsion
        maxredun= max(numredun,natom3)
      endif
!
! Set arrays for energy (actually, focka(nao3*nthread), fockb(nao3*nthread))
!
      call memset(nao3*(5+nthread*2)+nao2+nshell3,datacomp)
      allocate(h1mtrx(nao3),focka(nao3),fockb(nao3),smtrx(nao3),tmtrx(nao3), &
&              ortho(nao2),dmtrxa(nao3),dmtrxb(nao3),xint(nshell3))
!
! Set arrays for energy gradient and geometry optimization
!
      if(datajob%cartesian) then
        call memset(natom3+natom3*(natom3+1)/2,datacomp)
        allocate(egradold(natom3),ehess(natom3*(natom3+1)/2))
      else
        call memset(numredun*4+numredun*(numredun+1)/2,datacomp)
        allocate(coordredun(numredun*2),egradredun(numredun*2),ehess(numredun*(numredun+1)/2))
      endif
!
! Start geometry optimization cycle
!
      do iopt= 1,datajob%nopt
!
! Print geometry
!
        if(iopt >= 2) call writegeom(datamol,datacomp)
!
! Calculate nuclear repulsion energy
!
        call nucenergy(datajob%threshatom,datamol,datacomp)
        if(datacomp%master) then
          write(datacomp%iout,'(" Nuclear repulsion energy =",f18.9," Hartree",/)') datamol%enuc
        endif
!
! Set arrays 1
!
        call memset(nao2*2,datacomp)
        allocate(overinv(nao,nao),work(nao,nao))
!
! Calculate overlap and 1-electron integrals
!
        call oneei(h1mtrx,smtrx,tmtrx,work,datacomp%nproc2,datacomp%myrank2,datacomp%mpi_comm2, &
&                  datajob,datamol,databasis)
!
! Calculate canonical transformation matrix and inverse matrix of overlap
!
        call fullmtrx(smtrx,work,nao)
        call mtrxcanoninv(ortho,overinv,work,nao,datamol%nmo,datajob%threshover,datacomp)
!
! Calculate initial MOs
!
        if(iopt == 1) then
          call guessmo(cmoa,cmob,overinv,h1mtrx,ortho,datajob,datamol,databasis,datacomp)
        endif
!
! Unset arrays 1
!
        deallocate(overinv,work)
        call memunset(nao2*2,datacomp)
        call tstamp(1,datacomp)
!
! Calculate energy
!
        if(datajob%method == 'HF') then
          call calcuhf(h1mtrx,cmoa,cmob,focka,fockb,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa, &
&                      energymob,datajob,datamol,databasis,datacomp)
          if(.not.datacomp%convergedscf) go to 9999
        elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
          if((iopt == 1).and.(datajob%guess == 'HUCKEL')) then
            savedconv= datajob%dconv
            savecutint2= datajob%cutint2
            datajob%dconv= min(datajob%dconv*1.0D+04,1.0D-1)
            datajob%cutint2= min(datajob%cutint2*5.0D+02,1.0D-8)
            call calcuhf(h1mtrx,cmoa,cmob,focka,fockb,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa, &
&                        energymob,datajob,datamol,databasis,datacomp)
            if(.not.datacomp%convergedscf) go to 9999
            datajob%dconv= savedconv
            datajob%cutint2= savecutint2
            call tstamp(1,datacomp)
          endif
          call calcudft(h1mtrx,cmoa,cmob,focka,fockb,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa, &
&                       energymob,datajob,datamol,databasis,datacomp)
          if(.not.datacomp%convergedscf) go to 9999
        else
          if(datacomp%master) then
            write(datacomp%iout,'(" Error! This program does not support method= ",a16,".")') datajob%method
            call iabort(datacomp)
          endif
        endif
!
        if(((iopt == 1).and.(mod(datajob%iprint,10) >= 2)).or. &
&          (mod(datajob%iprint,10) >= 3)) then
          call writeeigenvalue(energymoa,energymob,2,datamol,datacomp)
        endif
        if(((iopt == 1).and.(mod(datajob%iprint,10) >= 2)).or. &
&          (mod(datajob%iprint,10) >= 5)) then
          call writeeigenvector(cmoa,energymoa,1,datajob,datamol,databasis,datacomp)
          call writeeigenvector(cmob,energymob,2,datajob,datamol,databasis,datacomp)
        endif
        call tstamp(1,datacomp)
!
! Calculate energy gradient
!
        if(datajob%method == 'HF') then
          call calcgraduhf(cmoa,cmob,energymoa,energymob,xint,egrad,datacomp%nproc1, &
&                          datacomp%myrank1,datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
        elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
          call calcgradudft(cmoa,cmob,energymoa,energymob,xint,egrad,datacomp%nproc1, &
&                           datacomp%myrank1,datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
        else
          if(datacomp%master) then
            write(datacomp%iout,'(" Error! This program does not support method= ",a16, &
&                                 " in energy gradient.")') datajob%method
            call iabort(datacomp)
          endif
        endif
!
! Calculate maximum and root mean square gradient values
!
        call calcmaxgrad(egradmax,egradrms,egrad,natom3)
        if(datacomp%master) &
&         write(datacomp%iout, &
&                 '(" -----------------------------------",/, &
&                   "   Geometry convergence check",/, &
&                   "     Optimization cycle:",i4," /",i4,/, &
&                   "     Optconv           =",f10.6,/, &
&                   " -----------------------------------",/, &
&                   "     Maximum gradient  =",f10.6,/, &
&                   "     RMS gradient      =",f10.6,/, &
&                   " -----------------------------------")') &
&                   iopt,datajob%nopt,datajob%optconv,egradmax,egradrms
!
! Check convergence
!
        if((egradmax <= datajob%optconv).and.(egradrms <= datajob%optconv*third)) then
          if(datacomp%master) &
&           write(datacomp%iout,'("   ============================",/, &
&                                 "   ==== Geometry converged ====",/, &
&                                 "   ============================",/)')
          datacomp%convergedgeom=.true.
          call tstamp(1,datacomp)
          exit
        endif
        call tstamp(1,datacomp)
!
! Write checkpoint file
!
        if(datacomp%master.and.(datajob%check /= '')) &
&         call writecheck(cmoa,cmob,dmtrxa,dmtrxb,energymoa,energymob,focka,fockb,smtrx,h1mtrx, &
&                         egrad,datajob,datamol,databasis,datacomp)
!
! Set work arrays 2
!
        if(datajob%cartesian) then
          call memset(natom3*3,datacomp)
          allocate(workv(natom3*3))
        else
          call memset(maxredun*maxredun*4+maxredun*3,datacomp)
          allocate(work(maxredun*maxredun,4),workv(maxredun*3))
        endif
!
! Calculate new coordinate
!
        if(datajob%cartesian) then
          call calcnewcoord(datamol%coord,datamol%coordold,egrad,egradold,ehess,workv,natom3,iopt, &
&                           datajob,datamol,datacomp)
        else
          call calcnewcoordred(datamol%coord,datamol%coordold,coordredun,egrad,egradredun,ehess, &
&                              work(1,1),work(1,2),work(1,3),work(1,4),workv,iopt,iredun,isizered, &
&                              maxredun,numbond,numangle,numtorsion,numredun, &
&                              datajob,datamol,datacomp)
        endif
!
! Unset work arrays 2
!
        if(datajob%cartesian) then
          deallocate(workv)
          call memunset(natom3*3,datacomp)
        else
          deallocate(work,workv)
          call memunset(maxredun*maxredun*4+maxredun*3,datacomp)
        endif
!
! Set guess MO calculation flag from Huckel to projection
!
!       call setnextopt(coordold,datamol%natom,iopt,datajob)
!
        if(iopt == datajob%nopt) then
          if(datacomp%master) then
            write(datacomp%iout,'(" ==============================")')
            write(datacomp%iout,'("   Geometry did not converge.")')
            write(datacomp%iout,'(" =============================="/)')
            call tstamp(1,datacomp)
          endif
          go to 9999
        endif
        call tstamp(1,datacomp)
      enddo
!
! End of optimization cycle 
!
!
! Print MO energies and coefficients if necessary
!
      if((iopt /= 1).and.(mod(datajob%iprint,10) >= 2)) then
        call writeeigenvalue(energymoa,energymob,2,datamol,datacomp)
      endif
      if(((iopt /= 1).and.(mod(datajob%iprint,10) >= 2)).and.(mod(datajob%iprint,10) <= 4)) then
        call writeeigenvector(cmoa,energymoa,1,datajob,datamol,databasis,datacomp)
        call writeeigenvector(cmob,energymob,2,datajob,datamol,databasis,datacomp)
      endif
!
      if(datajob%method /= 'MP2') then
        if(datacomp%master) then
          write(datacomp%iout,'(" ========================")')
          write(datacomp%iout,'("   Property calculation")')
          write(datacomp%iout,'(" ========================")')
        endif
!
! Calculate Mulliken and/or NBO charge
!
        select case(datajob%pop)
          case('MULLIKEN')
            call calcumulliken(dmtrxa,dmtrxb,smtrx,datamol,databasis,datacomp)
          case('NPA','NBO')
            call calcumulliken(dmtrxa,dmtrxb,smtrx,datamol,databasis,datacomp)
            call calcunpa(dmtrxa,dmtrxb,focka,fockb,smtrx,datamol,databasis,datacomp)
        end select
!
! Calculate dipole, quadrupole, and octupole moments
!
        select case(datajob%multipole)
          case('DIPOLE')
            call memset(nao3*6,datacomp)
            allocate(work(nao3,6))
            call calcudipole(work,work(1,4),dmtrxa,dmtrxb,datacomp%nproc1,datacomp%myrank1, &
&                            datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
            deallocate(work)
            call memunset(nao3*6,datacomp)
          case('OCTUPOLE')
            call memset(nao3*29,datacomp)
            allocate(work(nao3,29))
            call calcuoctupole(work,work(1,4),work(1,10),work(1,20),dmtrxa,dmtrxb, &
&                              datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1, &
&                              datajob,datamol,databasis,datacomp)
            deallocate(work)
            call memunset(nao3*29,datacomp)
        end select
      endif
!
! Write optimized geometry
!
      if(datacomp%master) then
        write(datacomp%iout,'(" ===========")')
        write(datacomp%iout,'("   Summary")')
        write(datacomp%iout,'(" ===========")')
        write(datacomp%iout,'(" ----------------------------------------------------")')
        if(datajob%method == 'HF') then
          write(datacomp%iout,'("   Final HF Energy =",f18.9," Hartree")') datamol%escf
        elseif(datajob%method == 'MP2') then
          write(datacomp%iout,'("   Final MP2 Energy =",f18.9," Hartree")')datamol%escf+datamol%emp2
        elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
          write(datacomp%iout,'("   Final DFT Energy =",f18.9," Hartree")') datamol%escf
        endif
        call writeoptgeom(datamol,datacomp)
      endif
!
9999 continue
!
! Write checkpoint file
!
      if(datacomp%master.and.(datajob%check /= '')) &
&       call writecheck(cmoa,cmob,dmtrxa,dmtrxb,energymoa,energymob,focka,fockb,smtrx,h1mtrx, &
&                       egrad,datajob,datamol,databasis,datacomp)
!
! Unset arrays for energy gradient and geometry optimization
!
      if(datajob%cartesian) then
        deallocate(egradold,ehess)
        call memunset(natom3+natom3*(natom3+1)/2,datacomp)
      else
        deallocate(coordredun,egradredun,ehess)
        call memunset(numredun*4+numredun*(numredun+1)/2,datacomp)
      endif
!
! Unset arrays for energy
!
      deallocate(h1mtrx,focka,fockb,smtrx,tmtrx, &
&                ortho,dmtrxa,dmtrxb,xint)
      call memunset(nao3*(5+nthread*2)+nao2+nshell3,datacomp)
      call tstamp(1,datacomp)
!
! Unset array for redundant coordinate
!
      if(.not.datajob%cartesian) then
        deallocate(iredun)
        call memunset(isizered,datacomp)
      endif
!
      return
end


!!---------------------------------------------
!  subroutine setnextopt(coordold,natom,iopt,datajob)
!!---------------------------------------------
!!
!! Set parameters for next optimization step
!!
!      use modguess, only : coord_g, nmo_g
!      use modmolecule, only : nmo
!      use modtype, only : typejob
!      implicit none
!      type(typejob),intent(inout) :: datajob
!      integer,intent(in) :: natom, iopt
!      integer :: iatom
!      real(8),intent(in) :: coordold(3,natom)
!!
!! Set MO projection as initial MO calculation
!!
!      datajob%guess= 'UPDATE'
!!
!! Copy coordinate and energy gradient
!!
!      do iatom= 1,natom
!        coord_g(1,iatom)= coordold(1,iatom)
!        coord_g(2,iatom)= coordold(2,iatom)
!        coord_g(3,iatom)= coordold(3,iatom)
!      enddo
!!
!      nmo_g= nmo
!!
!      return
!end


!--------------------------------------------------------
  subroutine setdft(datajob,datamol,databasis,datacomp)
!--------------------------------------------------------
!
! Set functional information
! Adjust the numbe of DFT grids when heavy elements are included
!
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typemol),intent(inout) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer :: ii, maxelem
!
      do ii= 1,9
        datamol%atomrad(-ii)= datajob%bqrad(ii)
      enddo
!
      select case(datajob%method)
        case('B3LYP')
          datajob%idftex = 1
          datajob%idftcor= 1
          datajob%hfexchange= 0.2D+00
        case('B3LYP5')
          datajob%idftex = 1
          datajob%idftcor= 2
          datajob%hfexchange= 0.2D+00
        case('HF','MP2')
        case default
          if(datacomp%master) then
            write(datacomp%iout,'(" Error! This program does not support method= ",a16,".")') datajob%method
            call iabort(datacomp)
          endif
      endselect
!
      if((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
        maxelem= maxval(datamol%numatomic(1:datamol%natom))
        if(((maxelem >= 55).or.(databasis%nao >= 2000)).and. &
&          ((datajob%nrad == 96).and.(datajob%nleb == 302))) then
          datacomp%nwarn= datacomp%nwarn+1
          if(datacomp%master) &
&           write(datacomp%iout,'(" Warning! The number of DFT grids may not be enough.")')
        endif
      endif
!
!     if(((datajob%idftex == 0).and.(datajob%idftcor == 0)).and.(datajob%guess == 'HF')) then
!       datajob%guess= 'HUCKEL'
!       if(datacomp%master) write(datacomp%iout,'(" Warning! Guess changes from HF to HUCKEL.")')
!       datacomp%nwarn= datacomp%nwarn+1
!     endif
!
      return
end


!-----------------------------------------------
  subroutine setmp2(datajob,datamol,databasis)
!-----------------------------------------------
!
! Set MP2 data
!
      use modtype, only : typejob, typemol, typebasis
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasis
      integer :: ncorecalc
!
      if(datajob%ncore == -1) datajob%ncore= ncorecalc(datamol,databasis)
!
      return
end
