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
!-----------------------------------------------------------
  subroutine smashmain(datajob,datamol,databasis,datacomp)
!-----------------------------------------------------------
!
! This is the main driver of Scalable Molecular Analysis Solver 
! for High performance computing systems (SMASH).
!
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typemol),intent(inout) :: datamol
      type(typebasis),intent(inout) :: databasis
      type(typecomp),intent(inout) :: datacomp
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
! Read input data and open argument input file if necessary
!
      if(datacomp%master) then
        call opendatfile(datacomp)
        if(command_argument_count() >= 1) call openinputfile(datacomp)
      endif
      call readinput(datajob,datamol,databasis,datacomp)
!
! Open checkpoint file and close argument input file if necessary
!
      if(datacomp%master) then
        if(datajob%check /= '') call opencheckfile(datajob,datacomp)
        if(command_argument_count() >= 1) call closeinputfile(datacomp)
      endif
!
! Set parameters and write computaional conditions
!
      call setdetails(datajob,datamol,databasis,datacomp)
!
! Start calculations
!
      select case(datajob%scftype)
        case('RHF')
          select case(datajob%runtype)
            case('ENERGY')
              call calcrenergy(datajob,datamol,databasis,datacomp)
            case('GRADIENT')
              call calcrgradient(datajob,datamol,databasis,datacomp)
            case('OPT')
              call calcrgeometry(datajob,datamol,databasis,datacomp)
            case default
              if(datacomp%master) then
                write(*,'(" Error! This program does not support runtype= ",a16,".")') &
&                       datajob%runtype
                call iabort
              endif
          end select
        case('UHF')
          select case(datajob%runtype)
            case('ENERGY')
              call calcuenergy(datajob,datamol,databasis,datacomp)
            case('GRADIENT')
              call calcugradient(datajob,datamol,databasis,datacomp)
            case('OPT')
              call calcugeometry(datajob,datamol,databasis,datacomp)
            case default
              if(datacomp%master) then
                write(*,'(" Error! This program does not support runtype= ",a16,".")') &
&                       datajob%runtype
                call iabort
              endif
          end select
        case default
          if(datacomp%master) &
&           write(*,'(" Error! SCFtype=",a16," is not supported.")') datajob%scftype
          call iabort
      end select
!
! Open, write, and close xyz file if set
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
! Check convergence of SCF and geometry optimization
!
      if(datacomp%master) then
        if(.not.datacomp%convergedscf) then
          write(datacomp%iout, &
&           '(" Used memory :",i8," MB",/, &
&             " =============================================",/, &
&             "    Error! SCF calculation did not converge!",/, &
&             " =============================================")') datacomp%memusedmax/125000
          return
!
        elseif(.not.datacomp%convergedgeom) then
          write(datacomp%iout, &
&           '(" Used memory :",i8," MB",/, &
&             " ================================================================",/, &
&             "    Error! Geometry optimization calculation did not converge!",/, &
&             " ================================================================")') &
&             datacomp%memusedmax/125000
          return
!
        else
          write(datacomp%iout, &
&           '(" Used memory :",i8," MB",/, &
&             " Your calculation finished successfully with",i3," warning(s).")') &
&             datacomp%memusedmax/125000, datacomp%nwarn
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
! Initialize variable for parallelization
!
      datacomp%master = .true.
!
! Start MPI parallelization and set mpi_comm1=MPI_COMM_WORLD
!
      call para_init(datacomp%mpi_comm1)
      call para_comm_size(datacomp%nproc1,datacomp%mpi_comm1)
      call para_comm_rank(datacomp%myrank1,datacomp%mpi_comm1)
!
      datacomp%nproc2= datacomp%nproc1
      datacomp%myrank2= datacomp%myrank1
      datacomp%mpi_comm2= datacomp%mpi_comm1
!
      if(datacomp%nproc1.gt.1) then
        datacomp%master =(datacomp%myrank1 == 0)
      endif
!
      return
end


!------------------------------------------------------------
  subroutine setdetails(datajob,datamol,databasis,datacomp)
!------------------------------------------------------------
!
! Read input file and set variables
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
! Set basis functions
!
      call setbasis(datajob,datamol,databasis,datacomp)
!
! Set ECP functions
!
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
          if(datajob%threshweight < zero) datajob%threshweight=1.0D-08
          if(datajob%threshrho    < zero) datajob%threshrho   =1.0D-06
          if(datajob%threshdfock  < zero) datajob%threshdfock =1.0D-05
          if(datajob%threshdftao  < zero) datajob%threshdftao =1.0D-04
          if(datajob%nrad == 0) datajob%nrad= 150
          if(datajob%nleb == 0) datajob%nleb= 590
        case('MEDIUM')
          if(datajob%cutint2 < zero) datajob%cutint2= 1.0D-11
          if(datajob%dconv   < zero) datajob%dconv  = 5.0D-06
          if(datajob%threshweight < zero) datajob%threshweight=1.0D-08
          if(datajob%threshrho    < zero) datajob%threshrho   =1.0D-05
          if(datajob%threshdfock  < zero) datajob%threshdfock =1.0D-04
          if(datajob%threshdftao  < zero) datajob%threshdftao =1.0D-03
          if(datajob%nrad == 0) datajob%nrad= 96
          if(datajob%nleb == 0) datajob%nleb= 302
        case('LOW')
          if(datajob%cutint2 < zero) datajob%cutint2= 1.0D-10
          if(datajob%dconv   < zero) datajob%dconv  = 1.0D-05
          if(datajob%threshweight < zero) datajob%threshweight=1.0D-08
          if(datajob%threshrho    < zero) datajob%threshrho   =1.0D-04
          if(datajob%threshdfock  < zero) datajob%threshdfock =1.0D-04
          if(datajob%threshdftao  < zero) datajob%threshdftao =1.0D-02
          if(datajob%nrad == 0) datajob%nrad= 72
          if(datajob%nleb == 0) datajob%nleb= 302
        case default
          if(datacomp%master) write(*,'(" Error! This program does not support precision= ", &
&                                       a16,".")') datajob%precision
          call iabort
      end select
!
      select case(datajob%output)
        case('MINIMUM')
          datajob%iprint= 1
        case('COMPACT')
          datajob%iprint= 2
        case('STANDARD')
          datajob%iprint= 3
        case('ALL')
          datajob%iprint= 4
        case('VERBOSE')
          datajob%iprint= 5
        case default
          if(datacomp%master) write(*,'(" Error! This program does not support output= ", &
&                                       a16,".")') datajob%output
          call iabort
      end select
!
      if((datajob%method == 'HF').and.(datajob%guess == 'HF')) then
        datajob%guess= 'HUCKEL'
        if(datacomp%master) write(datacomp%iout,'(" Warning! Guess changes from HF to HUCKEL.")')
        datacomp%nwarn= datacomp%nwarn+1
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
        if(datacomp%master) write(*,'(" Error! Spin multiplicity is ",i2, &
&                                     ", but number of elctrons is ",i5,".")') datamol%multi, nume
        call iabort
      endif
!
      return
end


!-------------------------------------------------------------
  subroutine calcrenergy(datajob,datamol,databasis,datacomp)
!-------------------------------------------------------------
!
! Driver of closed-shell energy calculation
!
! Parallel information
!
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typemol),intent(inout) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer :: nao, nao2, nao3, nshell3
      real(8), allocatable :: h1mtrx(:), fock(:), smtrx(:), tmtrx(:), cmo(:), ortho(:), dmtrx(:)
      real(8), allocatable :: xint(:), energymo(:)
      real(8), allocatable :: overinv(:), work(:)
      real(8) :: savedconv, savecutint2
!
      nao= databasis%nao
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(databasis%nshell*(databasis%nshell+1))/2
!
! Set arrays 1
!
      call memset(nao3*5+nao2*2+nshell3+nao,datacomp)
      allocate(h1mtrx(nao3),fock(nao3),smtrx(nao3),tmtrx(nao3),cmo(nao2),ortho(nao2),dmtrx(nao3),&
&              xint(nshell3),energymo(nao))
!
! Calculate nuclear repulsion energy
!
      call nucenergy(datajob%threshatom,datamol,datacomp)
      if(datacomp%master) then
        write(datacomp%iout,'(" Nuclear repulsion energy =",f15.8," Hartree",/)') datamol%enuc
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
! Calculate canonicalization and inverse overlap matrices
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
        if(.not.datacomp%convergedscf) then
          call tstamp(1,datacomp)
          return
        endif
        call writeeigenvalue(energymo,energymo,1,datajob,datamol,datacomp)
        call tstamp(1,datacomp)
      elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
        if(datajob%guess == 'HF') then
          savedconv= datajob%dconv
          savecutint2= datajob%cutint2
          datajob%dconv= max(datajob%dconv,1.0D-2)
          datajob%cutint2= max(datajob%cutint2,1.0D-9)
          call calcrhf(h1mtrx,cmo,fock,ortho,smtrx,dmtrx,xint,energymo, &
&                      datajob,datamol,databasis,datacomp)
          if(.not.datacomp%convergedscf) then
            call tstamp(1,datacomp)
            return
          endif
          datajob%dconv= savedconv
          datajob%cutint2= savecutint2
          call tstamp(1,datacomp)
        endif
        call calcrdft(h1mtrx,cmo,fock,ortho,smtrx,dmtrx,xint,energymo, &
&                     datajob,datamol,databasis,datacomp)
        if(.not.datacomp%convergedscf) then
          call tstamp(1,datacomp)
          return
        endif
        call writeeigenvalue(energymo,energymo,1,datajob,datamol,datacomp)
        call tstamp(1,datacomp)
      elseif(datajob%method == 'MP2') then
        call calcrhf(h1mtrx,cmo,fock,ortho,smtrx,dmtrx,xint,energymo, &
&                    datajob,datamol,databasis,datacomp)
        if(.not.datacomp%convergedscf) then
          call tstamp(1,datacomp)
          return
        endif
        call writeeigenvalue(energymo,energymo,1,datajob,datamol,datacomp)
        call tstamp(1,datacomp)
        call calcrmp2(cmo,energymo,xint,datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1, &
&                     datajob,datamol,databasis,datacomp)
        call tstamp(1,datacomp)
      else
        if(datacomp%master) then
          write(*,'(" Error! This program does not support method= ",a16,".")') datajob%method
          call iabort
        endif
      endif
!
! Print MOs
!
      if(datacomp%master.and.(mod(datajob%iprint,10) >= 2)) then
        write(datacomp%iout,'("  -------------------")')
        write(datacomp%iout,'("    MO coefficients")')
        write(datacomp%iout,'("  -------------------")')
        call writeeigenvector(cmo,energymo,datajob,datamol,databasis,datacomp)
      endif
!
! Calculate Mulliken charge
!
      call calcrmulliken(dmtrx,smtrx,datamol,databasis,datacomp)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(datajob%octupole) then
        call memset(nao3*29,datacomp)
        allocate(work(nao3*29))
        call calcroctupole(work,work(nao3*3+1),work(nao3*9+1),work(nao3*19+1),dmtrx, &
&                          datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1, &
&                          datajob,datamol,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*29,datacomp)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6,datacomp)
        allocate(work(nao3*6))
        call calcrdipole(work,work(nao3*3+1),dmtrx,datacomp%nproc1,datacomp%myrank1, &
&                        datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*6,datacomp)
      endif
!
! Write checkpoint file
!
      if(datacomp%master.and.(datajob%check /= '')) &
&       call writecheck(cmo,cmo,dmtrx,dmtrx,energymo,energymo,datajob,datamol,databasis,datacomp)
!
! Unset arrays 1
!
      deallocate(h1mtrx,fock,smtrx,tmtrx,cmo,ortho,dmtrx, &
&                xint,energymo)
      call memunset(nao3*5+nao2*2+nshell3+nao,datacomp)
      call tstamp(1,datacomp)
      return
end


!-------------------------------------------------------------
  subroutine calcuenergy(datajob,datamol,databasis,datacomp)
!-------------------------------------------------------------
!
! Driver of open-shell energy calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typemol),intent(inout) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer :: nao, nao2, nao3, nshell3
      real(8), allocatable :: h1mtrx(:), focka(:), fockb(:), smtrx(:), tmtrx(:), cmoa(:), cmob(:), ortho(:)
      real(8), allocatable :: dmtrxa(:), dmtrxb(:), xint(:), energymoa(:), energymob(:)
      real(8), allocatable :: overinv(:), work(:)
      real(8) :: savedconv, savecutint2
!
      nao= databasis%nao
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(databasis%nshell*(databasis%nshell+1))/2
!
! Set arrays 1
!
      call memset(nao3*7+nao2*3+nshell3+nao*2,datacomp)
      allocate(h1mtrx(nao3),focka(nao3),fockb(nao3),smtrx(nao3),tmtrx(nao3),cmoa(nao2),cmob(nao2),ortho(nao2),&
&              dmtrxa(nao3),dmtrxb(nao3),xint(nshell3),energymoa(nao),energymob(nao))
!
! Calculate nuclear repulsion energy
!
      call nucenergy(datajob%threshatom,datamol,datacomp)
      if(datacomp%master) then
        write(datacomp%iout,'(" Nuclear repulsion energy =",f15.8," Hartree",/)') datamol%enuc
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
! Calculate canonicalization and inverse overlap matrices
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
        call calcuhf(h1mtrx,cmoa,cmob,focka,fockb,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                    datajob,datamol,databasis,datacomp)
        if(.not.datacomp%convergedscf) then
          call tstamp(1,datacomp)
          return
        endif
        call writeeigenvalue(energymoa,energymob,2,datajob,datamol,datacomp)
        call tstamp(1,datacomp)
      elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
        if(datajob%guess == 'HF') then
          savedconv= datajob%dconv
          savecutint2= datajob%cutint2
          datajob%dconv= max(datajob%dconv,1.0D-2)
          datajob%cutint2= max(datajob%cutint2,1.0D-9)
          call calcuhf(h1mtrx,cmoa,cmob,focka,fockb,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                      datajob,datamol,databasis,datacomp)
          if(.not.datacomp%convergedscf) then
            call tstamp(1,datacomp)
            return
          endif
          datajob%dconv= savedconv
          datajob%cutint2= savecutint2
          call tstamp(1,datacomp)
        endif
        call calcudft(h1mtrx,cmoa,cmob,focka,fockb,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                     datajob,datamol,databasis,datacomp)
        if(.not.datacomp%convergedscf) then
          call tstamp(1,datacomp)
          return
        endif
        call writeeigenvalue(energymoa,energymob,2,datajob,datamol,datacomp)
        call tstamp(1,datacomp)
!     elseif(method == 'MP2') then
!       call calcuhf(h1mtrx,cmoa,ortho,smtrx,xint,energymoa)
!       call tstamp(1)
!       call calcump2(cmoa,energymoa,xint)
!       call tstamp(1)
      else
        if(datacomp%master) then
          write(*,'(" Error! This program does not support method= ",a16,".")')datajob%method
          call iabort
        endif
      endif
!
! Print MOs
!
      if(datacomp%master.and.(mod(datajob%iprint,10) >= 2)) then
        write(datacomp%iout,'("  -------------------------")')
        write(datacomp%iout,'("    Alpha MO coefficients")')
        write(datacomp%iout,'("  -------------------------")')
        call writeeigenvector(cmoa,energymoa,datajob,datamol,databasis,datacomp)
        write(datacomp%iout,'("  ------------------------")')
        write(datacomp%iout,'("    Beta MO coefficients")')
        write(datacomp%iout,'("  ------------------------")')
        call writeeigenvector(cmob,energymob,datajob,datamol,databasis,datacomp)
      endif
!
! Calculate Mulliken charge
!
      call calcumulliken(dmtrxa,dmtrxb,smtrx,datamol,databasis,datacomp)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(datajob%octupole) then
        call memset(nao3*29,datacomp)
        allocate(work(nao3*29))
        call calcuoctupole(work,work(nao3*3+1),work(nao3*9+1),work(nao3*19+1),dmtrxa,dmtrxb, &
&                          datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1, &
&                          datajob,datamol,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*29,datacomp)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6,datacomp)
        allocate(work(nao3*6))
        call calcudipole(work,work(nao3*3+1),dmtrxa,dmtrxb,datacomp%nproc1,datacomp%myrank1, &
&                        datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*6,datacomp)
      endif
!
! Write checkpoint file
!
      if(datacomp%master.and.(datajob%check /= '')) &
&        call writecheck(cmoa,cmob,dmtrxa,dmtrxb,energymoa,energymob,datajob,datamol,databasis,datacomp)
!
! Unset arrays 1
!
      deallocate(h1mtrx,focka,fockb,smtrx,tmtrx,cmoa,cmob,ortho, &
&                dmtrxa,dmtrxb,xint,energymoa,energymob)
      call memunset(nao3*7+nao2*3+nshell3+nao*2,datacomp)
      return
end


!---------------------------------------------------------------
  subroutine calcrgradient(datajob,datamol,databasis,datacomp)
!---------------------------------------------------------------
!
! Driver of energy gradient calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typemol),intent(inout) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer :: nao, nao2, nao3, nshell3
      real(8), allocatable :: h1mtrx(:), fock(:), smtrx(:), tmtrx(:), cmo(:), ortho(:), dmtrx(:)
      real(8), allocatable :: xint(:), energymo(:)
      real(8), allocatable :: overinv(:), work(:)
      real(8), allocatable :: egrad(:)
      real(8) :: egradmax, egradrms
      real(8) :: savedconv, savecutint2
!
      nao= databasis%nao
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(databasis%nshell*(databasis%nshell+1))/2
!
! Set arrays 1
!
      call memset(nao3*5+nao2*2+nshell3+nao,datacomp)
      allocate(h1mtrx(nao3),fock(nao3),smtrx(nao3),tmtrx(nao3),cmo(nao2),ortho(nao2),dmtrx(nao3),&
&              xint(nshell3),energymo(nao))
!
! Calculate nuclear repulsion energy
!
      call nucenergy(datajob%threshatom,datamol,datacomp)
      if(datacomp%master) then
        write(datacomp%iout,'(" Nuclear repulsion energy =",f15.8," Hartree",/)') datamol%enuc
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
! Calculate canonicalization and inverse overlap matrices
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
        if(.not.datacomp%convergedscf) then
          call tstamp(1,datacomp)
          return
        endif
        call writeeigenvalue(energymo,energymo,1,datajob,datamol,datacomp)
        call tstamp(1,datacomp)
      elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
        if(datajob%guess == 'HF') then
          savedconv= datajob%dconv
          savecutint2= datajob%cutint2
          datajob%dconv= max(datajob%dconv,1.0D-2)
          datajob%cutint2= max(datajob%cutint2,1.0D-9)
          call calcrhf(h1mtrx,cmo,fock,ortho,smtrx,dmtrx,xint,energymo, &
&                      datajob,datamol,databasis,datacomp)
          if(.not.datacomp%convergedscf) then
            call tstamp(1,datacomp)
            return
          endif
          datajob%dconv= savedconv
          datajob%cutint2= savecutint2
          call tstamp(1,datacomp)
        endif
        call calcrdft(h1mtrx,cmo,fock,ortho,smtrx,dmtrx,xint,energymo, &
&                     datajob,datamol,databasis,datacomp)
        if(.not.datacomp%convergedscf) then
          call tstamp(1,datacomp)
          return
        endif
        call writeeigenvalue(energymo,energymo,1,datajob,datamol,datacomp)
        call tstamp(1,datacomp)
      else
        if(datacomp%master) then
          write(*,'(" Error! This program does not support method= ",a16,".")') datajob%method
          call iabort
        endif
      endif
!
! Set arrays 3
!
      call memset(datamol%natom*3,datacomp)
      allocate(egrad(datamol%natom*3))
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
          write(*,'(" Error! This program does not support method= ",a16,".")') datajob%method
          call iabort
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
! Unset arrays 3
!
      deallocate(egrad)
      call memunset(datamol%natom*3,datacomp)
!
! Print MOs
!
      if(datacomp%master.and.(mod(datajob%iprint,10) >= 2)) then
        write(datacomp%iout,'("  -------------------")')
        write(datacomp%iout,'("    MO coefficients")')
        write(datacomp%iout,'("  -------------------")')
        call writeeigenvector(cmo,energymo,datajob,datamol,databasis,datacomp)
      endif
!
! Calculate Mulliken charge
!
      call calcrmulliken(dmtrx,smtrx,datamol,databasis,datacomp)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(datajob%octupole) then
        call memset(nao3*29,datacomp)
        allocate(work(nao3*29))
        call calcroctupole(work,work(nao3*3+1),work(nao3*9+1),work(nao3*19+1),dmtrx, &
&                          datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1, &
&                          datajob,datamol,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*29,datacomp)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6,datacomp)
        allocate(work(nao3*6))
        call calcrdipole(work,work(nao3*3+1),dmtrx,datacomp%nproc1,datacomp%myrank1, &
&                        datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*6,datacomp)
      endif
!
! Write checkpoint file
!
      if(datacomp%master.and.(datajob%check /= '')) &
&        call writecheck(cmo,cmo,dmtrx,dmtrx,energymo,energymo,datajob,datamol,databasis,datacomp)
!
! Unset arrays 1
!
      deallocate(h1mtrx,fock,smtrx,tmtrx,cmo,ortho,dmtrx, &
&                xint,energymo)
      call memunset(nao3*5+nao2*2+nshell3+nao,datacomp)
      call tstamp(1,datacomp)
      return
end


!---------------------------------------------------------------
  subroutine calcugradient(datajob,datamol,databasis,datacomp)
!---------------------------------------------------------------
!
! Driver of open-shell energy gradient calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typemol),intent(inout) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer :: nao, nao2, nao3, nshell3
      real(8), allocatable :: h1mtrx(:), focka(:), fockb(:), smtrx(:), tmtrx(:), cmoa(:), cmob(:), ortho(:)
      real(8), allocatable :: dmtrxa(:), dmtrxb(:), xint(:), energymoa(:), energymob(:)
      real(8), allocatable :: overinv(:), work(:)
      real(8), allocatable :: egrad(:)
      real(8) :: egradmax, egradrms
      real(8) :: savedconv, savecutint2
!
      nao= databasis%nao
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(databasis%nshell*(databasis%nshell+1))/2
!
! Set arrays 1
!
      call memset(nao3*7+nao2*3+nshell3+nao*2,datacomp)
      allocate(h1mtrx(nao3),focka(nao3),fockb(nao3),smtrx(nao3),tmtrx(nao3),cmoa(nao2),cmob(nao2),ortho(nao2),&
&              dmtrxa(nao3),dmtrxb(nao3),xint(nshell3),energymoa(nao),energymob(nao))
!
! Calculate nuclear repulsion energy
!
      call nucenergy(datajob%threshatom,datamol,datacomp)
      if(datacomp%master) then
        write(datacomp%iout,'(" Nuclear repulsion energy =",f15.8," Hartree",/)') datamol%enuc
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
! Calculate canonicalization and inverse overlap matrices
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
        call calcuhf(h1mtrx,cmoa,cmob,focka,fockb,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                    datajob,datamol,databasis,datacomp)
        if(.not.datacomp%convergedscf) then
          call tstamp(1,datacomp)
          return
        endif
        call writeeigenvalue(energymoa,energymob,2,datajob,datamol,datacomp)
        call tstamp(1,datacomp)
      elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
        if(datajob%guess == 'HF') then
          savedconv= datajob%dconv
          savecutint2= datajob%cutint2
          datajob%dconv= max(datajob%dconv,1.0D-2)
          datajob%cutint2= max(datajob%cutint2,1.0D-9)
          call calcuhf(h1mtrx,cmoa,cmob,focka,fockb,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                      datajob,datamol,databasis,datacomp)
          if(.not.datacomp%convergedscf) then
            call tstamp(1,datacomp)
            return
          endif
          datajob%dconv= savedconv
          datajob%cutint2= savecutint2
          call tstamp(1,datacomp)
        endif
        call calcudft(h1mtrx,cmoa,cmob,focka,fockb,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                     datajob,datamol,databasis,datacomp)
        if(.not.datacomp%convergedscf) then
          call tstamp(1,datacomp)
          return
        endif
        call writeeigenvalue(energymoa,energymob,2,datajob,datamol,datacomp)
        call tstamp(1,datacomp)
      else
        if(datacomp%master) then
          write(*,'(" Error! This program does not support method= ",a16,".")') datajob%method
          call iabort
        endif
      endif
!
! Set arrays 3
!
      call memset(datamol%natom*3,datacomp)
      allocate(egrad(datamol%natom*3))
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
          write(*,'(" Error! This program does not support method= ",a16," in energy gradient.")') &
&               datajob%method
          call iabort
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
! Unset arrays 3
!
      deallocate(egrad)
      call memunset(datamol%natom*3,datacomp)
!
! Print MOs
!
      if(datacomp%master.and.(mod(datajob%iprint,10) >= 2)) then
        write(datacomp%iout,'("  -------------------------")')
        write(datacomp%iout,'("    Alpha MO coefficients")')
        write(datacomp%iout,'("  -------------------------")')
        call writeeigenvector(cmoa,energymoa,datajob,datamol,databasis,datacomp)
        write(datacomp%iout,'("  ------------------------")')
        write(datacomp%iout,'("    Beta MO coefficients")')
        write(datacomp%iout,'("  ------------------------")')
        call writeeigenvector(cmob,energymob,datajob,datamol,databasis,datacomp)
      endif
!
! Calculate Mulliken charge
!
      call calcumulliken(dmtrxa,dmtrxb,smtrx,datamol,databasis,datacomp)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(datajob%octupole) then
        call memset(nao3*29,datacomp)
        allocate(work(nao3*29))
        call calcuoctupole(work,work(nao3*3+1),work(nao3*9+1),work(nao3*19+1),dmtrxa,dmtrxb, &
&                          datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1, &
&                          datajob,datamol,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*29,datacomp)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6,datacomp)
        allocate(work(nao3*6))
        call calcudipole(work,work(nao3*3+1),dmtrxa,dmtrxb,datacomp%nproc1,datacomp%myrank1, &
&                        datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*6,datacomp)
      endif
!
! Write checkpoint file
!
      if(datacomp%master.and.(datajob%check /= '')) &
&        call writecheck(cmoa,cmob,dmtrxa,dmtrxb,energymoa,energymob,datajob,datamol,databasis,datacomp)
!
! Unset arrays 1
!
      deallocate(h1mtrx,focka,fockb,smtrx,tmtrx,cmoa,cmob,ortho, &
&                dmtrxa,dmtrxb,xint,energymoa,energymob)
      call memunset(nao3*7+nao2*3+nshell3+nao*2,datacomp)
      call tstamp(1,datacomp)
      return
end


!-------------------------------------------------------------------------
  subroutine calcrgeometry(datajob,datamol,databasis,datacomp)
!-------------------------------------------------------------------------
!
! Driver of geometry optimization calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typemol),intent(inout) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer,allocatable :: iredun(:)
      integer :: nao, nao2, nao3, nshell3, natom3, ii, iopt
      integer :: isizered, numbond, numangle, numtorsion, numredun, maxredun
      real(8), parameter :: third=0.3333333333333333D+00
      real(8), allocatable :: h1mtrx(:), fock(:), smtrx(:), tmtrx(:), cmo(:), ortho(:), dmtrx(:)
      real(8), allocatable :: xint(:), energymo(:)
      real(8), allocatable :: egrad(:), egradold(:), ehess(:)
      real(8), allocatable :: overinv(:), work(:,:)
      real(8), allocatable :: workv(:), coordredun(:), egradredun(:)
      real(8) :: egradmax, egradrms
      real(8) :: savedconv, savecutint2
      logical :: exceed
!
      nao= databasis%nao
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(databasis%nshell*(databasis%nshell+1))/2
      natom3= datamol%natom*3
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
            write(*,'(" Error! The array size for redundant coordinate is too large.")')
            call iabort
          endif
        enddo
        numredun= numbond+numangle+numtorsion
        maxredun= max(numredun,natom3)
      endif
!
! Set arrays for energy
!
      call memset(nao3*5+nao2*2+nshell3+nao,datacomp)
      allocate(h1mtrx(nao3),fock(nao3),smtrx(nao3),tmtrx(nao3),cmo(nao2),ortho(nao2),dmtrx(nao3), &
&              xint(nshell3),energymo(nao))
!
! Set arrays for energy gradient and geometry optimization
!
      if(datajob%cartesian) then
        call memset(natom3*2+natom3*(natom3+1)/2,datacomp)
        allocate(egrad(natom3),egradold(natom3),ehess(natom3*(natom3+1)/2))
      else
        call memset(natom3+numredun*4+numredun*(numredun+1)/2,datacomp)
        allocate(egrad(natom3),coordredun(numredun*2),egradredun(numredun*2), &
&                ehess(numredun*(numredun+1)/2))
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
          write(datacomp%iout,'(" Nuclear repulsion energy =",f15.8," Hartree",/)') datamol%enuc
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
! Calculate canonicalization and inverse overlap matrices
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
          if(.not.datacomp%convergedscf) then
            call tstamp(1,datacomp)
            return
          endif
          if(iopt == 1) call writeeigenvalue(energymo,energymo,1,datajob,datamol,datacomp)
          call tstamp(1,datacomp)
        elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
          if((iopt == 1).and.(datajob%guess == 'HF')) then
            savedconv= datajob%dconv
            savecutint2= datajob%cutint2
            datajob%dconv= max(datajob%dconv,1.0D-2)
            datajob%cutint2= max(datajob%cutint2,1.0D-9)
            call calcrhf(h1mtrx,cmo,fock,ortho,smtrx,dmtrx,xint,energymo, &
&                        datajob,datamol,databasis,datacomp)
            if(.not.datacomp%convergedscf) then
              call tstamp(1,datacomp)
              return
            endif
            datajob%dconv= savedconv
            datajob%cutint2= savecutint2
            call tstamp(1,datacomp)
          endif
          call calcrdft(h1mtrx,cmo,fock,ortho,smtrx,dmtrx,xint,energymo, &
&                       datajob,datamol,databasis,datacomp)
          if(.not.datacomp%convergedscf) then
            call tstamp(1,datacomp)
            return
          endif
          if(iopt == 1) call writeeigenvalue(energymo,energymo,1,datajob,datamol,datacomp)
          call tstamp(1,datacomp)
        else
          if(datacomp%master) then
            write(*,'(" Error! This program does not support method= ",a16,".")') datajob%method
            call iabort
          endif
        endif
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
            write(*,'(" Error! This program does not support method= ",a16,".")') datajob%method
            call iabort
          endif
        endif
!
! Calculate maximum and root mean square gradient values
!
        call calcmaxgrad(egradmax,egradrms,egrad,natom3)
        if(datacomp%master) &
&         write(datacomp%iout,'(" -----------------------------------",/, &
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
          if(datacomp%master) write(datacomp%iout,'("   ==== Geometry converged ====",/)')
          datacomp%convergedgeom=.true.
          call tstamp(1,datacomp)
          exit
        endif
        call tstamp(1,datacomp)
!
! Write checkpoint file
!
        if(datacomp%master.and.(datajob%check /= '')) &
&         call writecheck(cmo,cmo,dmtrx,dmtrx,energymo,energymo,datajob,datamol,databasis,datacomp)
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
            write(datacomp%iout,'(" ------------------------------")')
            write(datacomp%iout,'("   Geometry did not converge.")')
            write(datacomp%iout,'(" ------------------------------"/)')
            call tstamp(1,datacomp)
          endif
          return
        endif
        call tstamp(1,datacomp)
      enddo
!
! End of optimization cycle 
!
      call writeeigenvalue(energymo,energymo,1,datajob,datamol,datacomp)
!
! Print MOs
!
      if(datacomp%master.and.(mod(datajob%iprint,10) >= 2)) then
        write(datacomp%iout,'("  -------------------")')
        write(datacomp%iout,'("    MO coefficients")')
        write(datacomp%iout,'("  -------------------")')
        call writeeigenvector(cmo,energymo,datajob,datamol,databasis,datacomp)
      endif
!
! Calculate Mulliken charge
!
      call calcrmulliken(dmtrx,smtrx,datamol,databasis,datacomp)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(datajob%octupole) then
        call memset(nao3*29,datacomp)
        allocate(work(nao3,29))
        call calcroctupole(work,work(1,4),work(1,10),work(1,20),dmtrx, &
&                          datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1, &
&                          datajob,datamol,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*29,datacomp)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6,datacomp)
        allocate(work(nao3,6))
        call calcrdipole(work,work(1,4),dmtrx,datacomp%nproc1,datacomp%myrank1, &
&                        datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*6,datacomp)
      endif
!
! Write optimized geometry
!
      if(datacomp%master) then
        write(datacomp%iout,'(" ==========================")')
        write(datacomp%iout,'("     Optimized Geometry")')
        write(datacomp%iout,'(" ==========================")')
        call writegeom(datamol,datacomp)
      endif
!
! Write checkpoint file
!
      if(datacomp%master.and.(datajob%check /= '')) &
&       call writecheck(cmo,cmo,dmtrx,dmtrx,energymo,energymo,datajob,datamol,databasis,datacomp)
!
! Unset arrays for energy gradient and geometry optimization
!
      if(datajob%cartesian) then
        deallocate(egrad,egradold,ehess)
        call memunset(natom3*2+natom3*(natom3+1)/2,datacomp)
      else
        deallocate(egrad,coordredun,egradredun, &
&                  ehess)
        call memunset(natom3+numredun*4+numredun*(numredun+1)/2,datacomp)
      endif
!
! Unset arrays for energy
!
      deallocate(h1mtrx,fock,smtrx,tmtrx,cmo,ortho,dmtrx, &
&                xint,energymo)
      call memunset(nao3*5+nao2*2+nshell3+nao,datacomp)
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


!----------------------------------------------------------------------------------------
  subroutine calcugeometry(datajob,datamol,databasis,datacomp)
!----------------------------------------------------------------------------------------
!
! Driver of open-shell geometry optimization calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typemol),intent(inout) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer,allocatable :: iredun(:)
      integer :: nao, nao2, nao3, nshell3, natom3, ii, iopt
      integer :: isizered, numbond, numangle, numtorsion, numredun, maxredun
      real(8), parameter :: third=0.3333333333333333D+00
      real(8), allocatable :: h1mtrx(:), focka(:), fockb(:), smtrx(:), tmtrx(:), cmoa(:), cmob(:), ortho(:)
      real(8), allocatable :: dmtrxa(:), dmtrxb(:), xint(:), energymoa(:), energymob(:)
      real(8), allocatable :: egrad(:), egradold(:), ehess(:)
      real(8), allocatable :: overinv(:,:), work(:,:)
      real(8), allocatable :: workv(:), coordredun(:), egradredun(:)
      real(8) :: egradmax, egradrms
      real(8) :: savedconv, savecutint2
      logical :: exceed
!
      nao= databasis%nao
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(databasis%nshell*(databasis%nshell+1))/2
      natom3= datamol%natom*3
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
            write(*,'(" Error! The array size for redundant coordinate is too large.")')
            call iabort
          endif
        enddo
        numredun= numbond+numangle+numtorsion
        maxredun= max(numredun,natom3)
      endif
!
! Set arrays for energy
!
      call memset(nao3*7+nao2*3+nshell3+nao*2,datacomp)
      allocate(h1mtrx(nao3),focka(nao3),fockb(nao3),smtrx(nao3),tmtrx(nao3),cmoa(nao2),cmob(nao2),ortho(nao2), &
&              dmtrxa(nao3),dmtrxb(nao3),xint(nshell3),energymoa(nao),energymob(nao))
!
! Set arrays for energy gradient and geometry optimization
!
      if(datajob%cartesian) then
        call memset(natom3*2+natom3*(natom3+1)/2,datacomp)
        allocate(egrad(natom3),egradold(natom3),ehess(natom3*(natom3+1)/2))
      else
        call memset(natom3+numredun*4+numredun*(numredun+1)/2,datacomp)
        allocate(egrad(natom3),coordredun(numredun*2),egradredun(numredun*2), &
&                ehess(numredun*(numredun+1)/2))
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
          write(datacomp%iout,'(" Nuclear repulsion energy =",f15.8," Hartree",/)') datamol%enuc
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
! Calculate canonicalization and inverse overlap matrices
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
          call calcuhf(h1mtrx,cmoa,cmob,focka,fockb,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                      datajob,datamol,databasis,datacomp)
          if(.not.datacomp%convergedscf) then
            call tstamp(1,datacomp)
            return
          endif
          if(iopt == 1) call writeeigenvalue(energymoa,energymob,2,datajob,datamol,datacomp)
          call tstamp(1,datacomp)
        elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
          if((iopt == 1).and.(datajob%guess == 'HF')) then
            savedconv= datajob%dconv
            savecutint2= datajob%cutint2
            datajob%dconv= max(datajob%dconv,1.0D-2)
            datajob%cutint2= max(datajob%cutint2,1.0D-9)
            call calcuhf(h1mtrx,cmoa,cmob,focka,fockb,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                        datajob,datamol,databasis,datacomp)
            if(.not.datacomp%convergedscf) then
              call tstamp(1,datacomp)
              return
            endif
            datajob%dconv= savedconv
            datajob%cutint2= savecutint2
            call tstamp(1,datacomp)
          endif
          call calcudft(h1mtrx,cmoa,cmob,focka,fockb,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                       datajob,datamol,databasis,datacomp)
          if(.not.datacomp%convergedscf) then
            call tstamp(1,datacomp)
            return
          endif
          if(iopt == 1) call writeeigenvalue(energymoa,energymob,2,datajob,datamol,datacomp)
          call tstamp(1,datacomp)
        else
          if(datacomp%master) then
            write(*,'(" Error! This program does not support method= ",a16,".")') datajob%method
            call iabort
          endif
        endif
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
            write(*,'(" Error! This program does not support method= ",a16," in energy gradient.")') &
&                 datajob%method
            call iabort
          endif
        endif
!
! Calculate maximum and root mean square gradient values
!
        call calcmaxgrad(egradmax,egradrms,egrad,natom3)
        if(datacomp%master) &
&         write(datacomp%iout,'(" -----------------------------------",/, &
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
          if(datacomp%master) write(datacomp%iout,'(" Geometry converged.",/)')
          datacomp%convergedgeom=.true.
          exit
          call tstamp(1,datacomp)
        endif
        call tstamp(1,datacomp)
!
! Write checkpoint file
!
        if(datacomp%master.and.(datajob%check /= '')) &
&          call writecheck(cmoa,cmob,dmtrxa,dmtrxb,energymoa,energymob,datajob,datamol,databasis,datacomp)
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
            write(datacomp%iout,'(" ------------------------------")')
            write(datacomp%iout,'("   Geometry did not converge.")')
            write(datacomp%iout,'(" ------------------------------"/)')
            call tstamp(1,datacomp)
          endif
          return
        endif
        call tstamp(1,datacomp)
      enddo
!
! End of optimization cycle 
!
      call writeeigenvalue(energymoa,energymob,2,datajob,datamol,datacomp)
!
! Print MOs
!
      if(datacomp%master.and.(mod(datajob%iprint,10) >= 2)) then
        write(datacomp%iout,'("  -------------------------")')
        write(datacomp%iout,'("    Alpha MO coefficients")')
        write(datacomp%iout,'("  -------------------------")')
        call writeeigenvector(cmoa,energymoa,datajob,datamol,databasis,datacomp)
        write(datacomp%iout,'("  ------------------------")')
        write(datacomp%iout,'("    Beta MO coefficients")')
        write(datacomp%iout,'("  ------------------------")')
        call writeeigenvector(cmob,energymob,datajob,datamol,databasis,datacomp)
      endif
!
! Calculate Mulliken charge
!
      call calcumulliken(dmtrxa,dmtrxb,smtrx,datamol,databasis,datacomp)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(datajob%octupole) then
        call memset(nao3*29,datacomp)
        allocate(work(nao3,29))
        call calcuoctupole(work,work(1,4),work(1,10),work(1,20),dmtrxa,dmtrxb, &
&                          datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1, &
&                          datajob,datamol,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*29,datacomp)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6,datacomp)
        allocate(work(nao3,6))
        call calcudipole(work,work(1,4),dmtrxa,dmtrxb,datacomp%nproc1,datacomp%myrank1, &
&                        datacomp%mpi_comm1,datajob,datamol,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*6,datacomp)
      endif
!
! Write optimized geometry
!
      if(datacomp%master) then
        write(datacomp%iout,'(" ==========================")')
        write(datacomp%iout,'("     Optimized Geometry")')
        write(datacomp%iout,'(" ==========================")')
        call writegeom(datamol,datacomp)
      endif
!
! Write checkpoint file
!
      if(datacomp%master.and.(datajob%check /= '')) &
&        call writecheck(cmoa,cmob,dmtrxa,dmtrxb,energymoa,energymob,datajob,datamol,databasis,datacomp)
!
! Unset arrays for energy gradient and geometry optimization
!
      if(datajob%cartesian) then
        deallocate(egrad,egradold,ehess)
        call memunset(natom3*2+natom3*(natom3+1)/2,datacomp)
      else
        deallocate(egrad,coordredun,egradredun, &
&                ehess)
        call memunset(natom3+numredun*4+numredun*(numredun+1)/2,datacomp)
      endif
!
! Unset arrays for energy
!
      deallocate(h1mtrx,focka,fockb,smtrx,tmtrx,cmoa,cmob,ortho, &
&                dmtrxa,dmtrxb,xint,energymoa,energymob)
      call memunset(nao3*7+nao2*3+nshell3+nao*2,datacomp)
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
            write(*,'(" Error! This program does not support method= ",a16,".")') datajob%method
            call iabort
          endif
      endselect
!
      if((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
        maxelem= maxval(datamol%numatomic(1:datamol%natom))
        if(((maxelem >= 55).or.(databasis%nao >= 2000)).and. &
&          ((datajob%nrad == 96).and.(datajob%nleb == 302))) then
          datacomp%nwarn= datacomp%nwarn+1
          if(datacomp%master) write(datacomp%iout,'(" Warning! The number of DFT grids may not be enough.")')
        endif
      endif
!
      return
end


!-----------------------------------------------
  subroutine setmp2(datajob,datamol,databasis)
!-----------------------------------------------
!
! Set MP2 information
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
