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
!---------------
  program main
!---------------
!
! This is the main driver of Scalable Molecular Analysis Solver 
! for High performance computing systems (SMASH).
!
      use modparam, only : input, icheck
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob) :: datajob
      type(typemol) :: datamol
      type(typebasis) :: databasis
      type(typecomp) :: datacomp
      logical :: converged
!
      call setparallel(datacomp)
!
      if(datacomp%master) then
        write(*,&
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
! Read input file and set details
!
      call setdetails(datajob,datamol,databasis,datacomp)
!
! Start calculations
!
      select case(datajob%scftype)
        case('RHF')
          select case(datajob%runtype)
            case('ENERGY')
              call calcrenergy(datajob,databasis,datacomp)
            case('GRADIENT')
              call calcrgradient(datajob,databasis,datacomp)
            case('OPTIMIZE')
              call calcrgeometry(converged,datajob,databasis,datacomp)
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
              call calcuenergy(datajob,databasis,datacomp)
            case('GRADIENT')
              call calcugradient(datajob,databasis,datacomp)
            case('OPTIMIZE')
              call calcugeometry(converged,datajob,databasis,datacomp)
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
! Close input.dat and checkpoint files
!
      if(datacomp%master) close(unit=input,status='DELETE')
      if(datacomp%master.and.(datajob%check /= '')) close(unit=icheck)
!
      call para_finalize
      call memcheck(datacomp)
      call tstamp(2,datacomp)
!
      if(datacomp%master) then
        write(*,'(" Used memory :",1x,i6," MB")') datacomp%memusedmax/125000
        if((datajob%runtype =='OPTIMIZE').and.(.not.converged)) then
          write(*,'(/," ============================================================")')
          write(*,'("  Geometry optimization did not finish with",i3," warning(s)!")') &
&                 datacomp%nwarn
          write(*,'(" ============================================================")')
        else
          write(*,'(" Your calculation finished with",i3," warning(s).")') datacomp%nwarn
        endif
      endif
end program main


!-------------------------
  subroutine setparallel(datacomp)
!-------------------------
!
!  Initialize MPI execution environment
!
      use modtype, only : typecomp
      implicit none
      type(typecomp),intent(inout) :: datacomp
!
! Initialize variables for parallelization
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


!----------------------------------
  subroutine setdetails(datajob,datamol,databasis,datacomp)
!----------------------------------
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
! Set defaults before reading input file
!
      call setdefault1(datacomp)
!
! Read input data and open checkpoint file if necessary
!
      if(datacomp%master) call opendatfile(datacomp)
      call readinput(datajob,datamol,databasis,datacomp)
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
      call setelectron(datajob,databasis,datacomp)
!
! Reset defaults after reading input file
!
      call setdefault2(datajob,datacomp)
!
! Set functional information and adjust the number of DFT grids
!
      call setdft(datajob,databasis,datacomp)
!
! Set functional information and adjust the number of DFT grids
!
      call setmp2(datajob,databasis)
!
! Write input data
!
      call writecondition(datajob,databasis,datacomp)
      call writegeom(datacomp)
      call writebasis(datajob,databasis,datacomp)
      if(datajob%flagecp) call writeecp(datajob,databasis,datacomp)
!
! Set atom charge including dummy atom
!
      call setcharge(datacomp)
!
      return
end


!-------------------------
  subroutine setdefault1(datacomp)
!-------------------------
!
! Set defaults before reading input file
!
      use modmolecule, only : multi, charge
      use modtype, only : typecomp
      implicit none
      type(typecomp),intent(inout) :: datacomp
!
!     datacomp%nwarn= 0
!     datacomp%memmax = 1000000000
!     datacomp%memused= 0
!     datacomp%memusedmax= 0
!     memory = ''
!     maxiter= 150
!     maxdiis= 20
!     maxsoscf= 20
!     maxqc   = 15
!     maxqcdiag= 100
!     maxqcdiagsub= 10
!     fdiff  =.true.
!     scfconv='DIIS'
!     extrap =.false.
!     threshsoscf= 0.25D+00
!     threshqc   = 1.0D-05
!     threshover = 1.0D-06
!     threshatom = 2.0D-01
!     threshdiis = 6.0D-01
!     threshmp2cphf=1.0D-10
!     idftex = 0
!     idftcor= 0
!     iprint = 2
!     bohr   =.false.
!     spher  =.true.
!     spher_g=.true.
!     nopt   = 100
!     optconv= 1.0D-04
!     cartesian=.false.
      multi  = 1
      charge = 0.0D+00
!     bqrad(:)=1.0D+00
!     nvfz= 0
!     maxmp2diis= 20
!     maxmp2iter= 100
!
!     cutint2=-1.0d+00
!     nrad = 0
!     nleb = 0
!     ncore= -1
!     dconv=-1.0D+00
!     threshweight=-1.0D+00
!     threshrho=-1.0D+00
!     threshdfock=-1.0D+00
!     threshdftao=-1.0D+00
!
!     precision='MEDIUM'
!     flagecp= .false.
!     scftype='RHF'
!     method='HARTREE-FOCK'
!     runtype='ENERGY'
!     basis='STO-3G'
!     guess='HUCKEL'
!     ecp=''
!     check=''
!     octupole=.false.
!
      return
end


!-------------------------
  subroutine setdefault2(datajob,datacomp)
!-------------------------
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
      return
end


!-------------------------
  subroutine setelectron(datajob,databasis,datacomp)
!-------------------------
!
! Set number of electrons
!
      use modmolecule, only : numatomic, neleca, nelecb, natom, multi, charge
      use modtype, only : typejob, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer :: nume, ii
!
! Calculate total number of electrons
!
      nume= -nint(charge)
      do ii= 1,natom
        if(numatomic(ii) > 0) nume= nume+numatomic(ii)
      enddo
!
! Subtract electrons of core potentials
!
      if(datajob%flagecp) then
        do ii= 1,natom
          nume= nume-databasis%izcore(ii)
        enddo
      endif

!
! Calculate numbers of alpha and beta electrons
!
      if((datajob%scftype == 'RHF').and.(multi /= 1)) then
        if(datacomp%master) write(*,'(" Warning! SCFtype changes from RHF to UHF.")')
        datajob%scftype = 'UHF'
        datacomp%nwarn= datacomp%nwarn+1
      endif
!
      neleca=(nume+multi-1)/2
      nelecb=(nume-multi+1)/2
      if((neleca+nelecb)/= nume) then
        if(datacomp%master) write(*,'(" Error! Spin multiplicity is ",i2, &
&                                     ", but number of elctrons is ",i5,".")')multi, nume
        call iabort
      endif
!
      return
end


!----------------------------------------------------------------------------
  subroutine calcrenergy(datajob,databasis,datacomp)
!----------------------------------------------------------------------------
!
! Driver of closed-shell energy calculation
!
! Parallel information
!
      use modmolecule, only : nmo, enuc
      use modtype, only : typejob, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer :: nao, nao2, nao3, nshell3
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmo(:), ortho(:), dmtrx(:)
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
      call memset(nao3*4+nao2*2+nshell3+nao,datacomp)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmo(nao2),ortho(nao2),dmtrx(nao3),&
&              xint(nshell3),energymo(nao))
!
! Calculate nuclear repulsion energy
!
      call nucenergy(datajob%threshatom,datacomp)
      if(datacomp%master) then
        write(*,'(" Nuclear repulsion energy =",f15.8," a.u.",/)') enuc
      endif
!
! Set arrays 2
!
      call memset(nao2*2,datacomp)
      allocate(overinv(nao2),work(nao2))
!
! Calculate overlap and 1-electron integrals
!
      call oneei(h1mtrx,smtrx,tmtrx,work,datacomp%nproc2,datacomp%myrank2,datacomp%mpi_comm2,datajob,databasis)
!
! Calculate canonicalization and inverse overlap matrices
!
      call fullmtrx(smtrx,work,nao)
      call mtrxcanoninv(ortho,overinv,work,nao,nmo,datajob%threshover,datacomp)
!
! Calculate initial MOs
!
      call guessmo(cmo,cmo,overinv,h1mtrx,ortho,datajob,databasis,datacomp)
!
! Unset arrays 2
!
      deallocate(overinv,work)
      call memunset(nao2*2,datacomp)
      call tstamp(1,datacomp)
!
! Start SCF
!
      if(datajob%method == 'HARTREE-FOCK') then
        call calcrhf(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo,datajob,databasis,datacomp)
        call writeeigenvalue(energymo,energymo,1,datajob,datacomp)
        call tstamp(1,datacomp)
      elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
        if(datajob%guess == 'HUCKEL') then
          savedconv= datajob%dconv
          savecutint2= datajob%cutint2
          datajob%dconv= max(datajob%dconv,1.0D-2)
          datajob%cutint2= max(datajob%cutint2,1.0D-9)
          call calcrhf(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo,datajob,databasis,datacomp)
          datajob%dconv= savedconv
          datajob%cutint2= savecutint2
          call tstamp(1,datacomp)
        endif
        call calcrdft(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo,datajob,databasis,datacomp)
        call writeeigenvalue(energymo,energymo,1,datajob,datacomp)
        call tstamp(1,datacomp)
      elseif(datajob%method == 'MP2') then
        call calcrhf(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo,datajob,databasis,datacomp)
        call writeeigenvalue(energymo,energymo,1,datajob,datacomp)
        call tstamp(1,datacomp)
        call calcrmp2(cmo,energymo,xint,datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
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
      if(datacomp%master.and.(datajob%iprint >= 2)) then
        write(*,'("  -------------------")')
        write(*,'("    MO coefficients")')
        write(*,'("  -------------------")')
        call writeeigenvector(cmo,energymo,databasis,datacomp)
      endif
!
! Calculate Mulliken charge
!
      call calcrmulliken(dmtrx,smtrx,databasis,datacomp)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(datajob%octupole) then
        call memset(nao3*29,datacomp)
        allocate(work(nao3*29))
        call calcroctupole(work,work(nao3*3+1),work(nao3*9+1),work(nao3*19+1),dmtrx, &
&                          datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*29,datacomp)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6,datacomp)
        allocate(work(nao3*6))
        call calcrdipole(work,work(nao3*3+1),dmtrx,datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*6,datacomp)
      endif
!
! Write checkpoint file
!
      if(datacomp%master.and.(datajob%check /= '')) call writecheck(cmo,cmo,dmtrx,dmtrx,energymo,energymo,datajob,databasis)
!
! Unset arrays 1
!
      deallocate(h1mtrx,smtrx,tmtrx,cmo,ortho,dmtrx, &
&                xint,energymo)
      call memunset(nao3*4+nao2*2+nshell3+nao,datacomp)
      call tstamp(1,datacomp)
      return
end


!----------------------------------------------------------------------------
  subroutine calcuenergy(datajob,databasis,datacomp)
!----------------------------------------------------------------------------
!
! Driver of open-shell energy calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
      use modmolecule, only : nmo, enuc
      use modtype, only : typejob, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer :: nao, nao2, nao3, nshell3
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmoa(:), cmob(:), ortho(:)
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
      call memset(nao3*5+nao2*3+nshell3+nao*2,datacomp)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmoa(nao2),cmob(nao2),ortho(nao2),&
&              dmtrxa(nao3),dmtrxb(nao3),xint(nshell3),energymoa(nao),energymob(nao))
!
! Calculate nuclear repulsion energy
!
      call nucenergy(datajob%threshatom,datacomp)
      if(datacomp%master) then
        write(*,'(" Nuclear repulsion energy =",f15.8," a.u.",/)') enuc
      endif
!
! Set arrays 2
!
      call memset(nao2*2,datacomp)
      allocate(overinv(nao2),work(nao2))
!
! Calculate overlap and 1-electron integrals
!
      call oneei(h1mtrx,smtrx,tmtrx,work,datacomp%nproc2,datacomp%myrank2,datacomp%mpi_comm2,datajob,databasis)
!
! Calculate canonicalization and inverse overlap matrices
!
      call fullmtrx(smtrx,work,nao)
      call mtrxcanoninv(ortho,overinv,work,nao,nmo,datajob%threshover,datacomp)
!
! Calculate initial MOs
!
      call guessmo(cmoa,cmob,overinv,h1mtrx,ortho,datajob,databasis,datacomp)
!
! Unset arrays 2
!
      deallocate(overinv,work)
      call memunset(nao2*2,datacomp)
      call tstamp(1,datacomp)
!
! Start SCF
!
      if(datajob%method == 'HARTREE-FOCK') then
        call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob,datajob,databasis,datacomp)
        call writeeigenvalue(energymoa,energymob,2,datajob,datacomp)
        call tstamp(1,datacomp)
      elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
        if(datajob%guess == 'HUCKEL') then
          savedconv= datajob%dconv
          savecutint2= datajob%cutint2
          datajob%dconv= max(datajob%dconv,1.0D-2)
          datajob%cutint2= max(datajob%cutint2,1.0D-9)
          call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob,datajob,databasis,datacomp)
          datajob%dconv= savedconv
          datajob%cutint2= savecutint2
          call tstamp(1,datacomp)
        endif
        call calcudft(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob,datajob,databasis,datacomp)
        call writeeigenvalue(energymoa,energymob,2,datajob,datacomp)
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
      if(datacomp%master.and.(datajob%iprint >= 2)) then
        write(*,'("  -------------------------")')
        write(*,'("    Alpha MO coefficients")')
        write(*,'("  -------------------------")')
        call writeeigenvector(cmoa,energymoa,databasis,datacomp)
        write(*,'("  ------------------------")')
        write(*,'("    Beta MO coefficients")')
        write(*,'("  ------------------------")')
        call writeeigenvector(cmob,energymob,databasis,datacomp)
      endif
!
! Calculate Mulliken charge
!
      call calcumulliken(dmtrxa,dmtrxb,smtrx,databasis,datacomp)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(datajob%octupole) then
        call memset(nao3*29,datacomp)
        allocate(work(nao3*29))
        call calcuoctupole(work,work(nao3*3+1),work(nao3*9+1),work(nao3*19+1),dmtrxa,dmtrxb, &
&                          datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*29,datacomp)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6,datacomp)
        allocate(work(nao3*6))
        call calcudipole(work,work(nao3*3+1),dmtrxa,dmtrxb,datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*6,datacomp)
      endif
!
! Write checkpoint file
!
      if(datacomp%master.and.(datajob%check /= '')) call writecheck(cmoa,cmob,dmtrxa,dmtrxb,energymoa,energymob,datajob,databasis)
!
! Unset arrays 1
!
      deallocate(h1mtrx,smtrx,tmtrx,cmoa,cmob,ortho, &
&                dmtrxa,dmtrxb,xint,energymoa,energymob)
      call memunset(nao3*5+nao2*3+nshell3+nao*2,datacomp)
      return
end


!------------------------------------------------------------------------------
  subroutine calcrgradient(datajob,databasis,datacomp)
!------------------------------------------------------------------------------
!
! Driver of energy gradient calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
      use modmolecule, only : nmo, natom, enuc
      use modtype, only : typejob, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer :: nao, nao2, nao3, nshell3
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmo(:), ortho(:), dmtrx(:)
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
      call memset(nao3*4+nao2*2+nshell3+nao,datacomp)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmo(nao2),ortho(nao2),dmtrx(nao3),&
&              xint(nshell3),energymo(nao))
!
! Calculate nuclear repulsion energy
!
      call nucenergy(datajob%threshatom,datacomp)
      if(datacomp%master) then
        write(*,'(" Nuclear repulsion energy =",f15.8," a.u.",/)') enuc
      endif
!
! Set arrays 2
!
      call memset(nao2*2,datacomp)
      allocate(overinv(nao2),work(nao2))
!
! Calculate overlap and 1-electron integrals
!
      call oneei(h1mtrx,smtrx,tmtrx,work,datacomp%nproc2,datacomp%myrank2,datacomp%mpi_comm2,datajob,databasis)
!
! Calculate canonicalization and inverse overlap matrices
!
      call fullmtrx(smtrx,work,nao)
      call mtrxcanoninv(ortho,overinv,work,nao,nmo,datajob%threshover,datacomp)
!
! Calculate initial MOs
!
      call guessmo(cmo,cmo,overinv,h1mtrx,ortho,datajob,databasis,datacomp)
!
! Unset arrays 2
!
      deallocate(overinv,work)
      call memunset(nao2*2,datacomp)
      call tstamp(1,datacomp)
!
! Start SCF
!
      if((datajob%method == 'HARTREE-FOCK').or.(datajob%method == 'MP2')) then
        call calcrhf(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo,datajob,databasis,datacomp)
        call writeeigenvalue(energymo,energymo,1,datajob,datacomp)
        call tstamp(1,datacomp)
      elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
        if(datajob%guess == 'HUCKEL') then
          savedconv= datajob%dconv
          savecutint2= datajob%cutint2
          datajob%dconv= max(datajob%dconv,1.0D-2)
          datajob%cutint2= max(datajob%cutint2,1.0D-9)
          call calcrhf(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo,datajob,databasis,datacomp)
          datajob%dconv= savedconv
          datajob%cutint2= savecutint2
          call tstamp(1,datacomp)
        endif
        call calcrdft(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo,datajob,databasis,datacomp)
        call writeeigenvalue(energymo,energymo,1,datajob,datacomp)
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
      call memset(natom*3,datacomp)
      allocate(egrad(natom*3))
!
! Calculate energy gradient
!
      if(datajob%method == 'HARTREE-FOCK') then
        call calcgradrhf(cmo,energymo,xint,egrad,datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
        call tstamp(1,datacomp)
      elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
        call calcgradrdft(cmo,energymo,xint,egrad,datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
        call tstamp(1,datacomp)
      elseif(datajob%method == 'MP2') then
        call calcgradrmp2(cmo,energymo,xint,egrad,datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
        call tstamp(1,datacomp)
      else
        if(datacomp%master) then
          write(*,'(" Error! This program does not support method= ",a16,".")') datajob%method
          call iabort
        endif
      endif
!
! Calculate maximum and root mean square gradient values
!
      call calcmaxgrad(egradmax,egradrms,egrad,natom*3)
      if(datacomp%master) write(*,'(" Maximum gradient =",f13.8,"  RMS gradient =",f13.8,/)') &
&                      egradmax,egradrms
!
! Unset arrays 3
!
      deallocate(egrad)
      call memunset(natom*3,datacomp)
!
! Print MOs
!
      if(datacomp%master.and.(datajob%iprint >= 2)) then
        write(*,'("  -------------------")')
        write(*,'("    MO coefficients")')
        write(*,'("  -------------------")')
        call writeeigenvector(cmo,energymo,databasis,datacomp)
      endif
!
! Calculate Mulliken charge
!
      call calcrmulliken(dmtrx,smtrx,databasis,datacomp)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(datajob%octupole) then
        call memset(nao3*29,datacomp)
        allocate(work(nao3*29))
        call calcroctupole(work,work(nao3*3+1),work(nao3*9+1),work(nao3*19+1),dmtrx, &
&                          datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*29,datacomp)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6,datacomp)
        allocate(work(nao3*6))
        call calcrdipole(work,work(nao3*3+1),dmtrx,datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*6,datacomp)
      endif
!
! Write checkpoint file
!
      if(datacomp%master.and.(datajob%check /= '')) call writecheck(cmo,cmo,dmtrx,dmtrx,energymo,energymo,datajob,databasis)
!
! Unset arrays 1
!
      deallocate(h1mtrx,smtrx,tmtrx,cmo,ortho,dmtrx, &
&                xint,energymo)
      call memunset(nao3*4+nao2*2+nshell3+nao,datacomp)
      call tstamp(1,datacomp)
      return
end


!------------------------------------------------------------------------------
  subroutine calcugradient(datajob,databasis,datacomp)
!------------------------------------------------------------------------------
!
! Driver of open-shell energy gradient calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
      use modmolecule, only : nmo, natom, enuc
      use modtype, only : typejob, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer :: nao, nao2, nao3, nshell3
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmoa(:), cmob(:), ortho(:)
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
      call memset(nao3*5+nao2*3+nshell3+nao*2,datacomp)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmoa(nao2),cmob(nao2),ortho(nao2),&
&              dmtrxa(nao3),dmtrxb(nao3),xint(nshell3),energymoa(nao),energymob(nao))
!
! Calculate nuclear repulsion energy
!
      call nucenergy(datajob%threshatom,datacomp)
      if(datacomp%master) then
        write(*,'(" Nuclear repulsion energy =",f15.8," a.u.",/)') enuc
      endif
!
! Set arrays 2
!
      call memset(nao2*2,datacomp)
      allocate(overinv(nao2),work(nao2))
!
! Calculate overlap and 1-electron integrals
!
      call oneei(h1mtrx,smtrx,tmtrx,work,datacomp%nproc2,datacomp%myrank2,datacomp%mpi_comm2,datajob,databasis)
!
! Calculate canonicalization and inverse overlap matrices
!
      call fullmtrx(smtrx,work,nao)
      call mtrxcanoninv(ortho,overinv,work,nao,nmo,datajob%threshover,datacomp)
!
! Calculate initial MOs
!
      call guessmo(cmoa,cmob,overinv,h1mtrx,ortho,datajob,databasis,datacomp)
!
! Unset arrays 2
!
      deallocate(overinv,work)
      call memunset(nao2*2,datacomp)
      call tstamp(1,datacomp)
!
! Start SCF
!
      if(datajob%method == 'HARTREE-FOCK') then
        call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob,datajob,databasis,datacomp)
        call writeeigenvalue(energymoa,energymob,2,datajob,datacomp)
        call tstamp(1,datacomp)
      elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
        if(datajob%guess == 'HUCKEL') then
          savedconv= datajob%dconv
          savecutint2= datajob%cutint2
          datajob%dconv= max(datajob%dconv,1.0D-2)
          datajob%cutint2= max(datajob%cutint2,1.0D-9)
          call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob,datajob,databasis,datacomp)
          datajob%dconv= savedconv
          datajob%cutint2= savecutint2
          call tstamp(1,datacomp)
        endif
        call calcudft(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob,datajob,databasis,datacomp)
        call writeeigenvalue(energymoa,energymob,2,datajob,datacomp)
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
      call memset(natom*3,datacomp)
      allocate(egrad(natom*3))
!
! Calculate energy gradient
!
      if(datajob%method == 'HARTREE-FOCK') then
        call calcgraduhf(cmoa,cmob,energymoa,energymob,xint,egrad,datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
      elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
        call calcgradudft(cmoa,cmob,energymoa,energymob,xint,egrad,datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
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
      call calcmaxgrad(egradmax,egradrms,egrad,natom*3)
      if(datacomp%master) write(*,'(" Maximum gradient =",f13.8,"  RMS gradient =",f13.8,/)') &
&                      egradmax,egradrms
!
! Unset arrays 3
!
      deallocate(egrad)
      call memunset(natom*3,datacomp)
!
! Print MOs
!
      if(datacomp%master.and.(datajob%iprint >= 2)) then
        write(*,'("  -------------------------")')
        write(*,'("    Alpha MO coefficients")')
        write(*,'("  -------------------------")')
        call writeeigenvector(cmoa,energymoa,databasis,datacomp)
        write(*,'("  ------------------------")')
        write(*,'("    Beta MO coefficients")')
        write(*,'("  ------------------------")')
        call writeeigenvector(cmob,energymob,databasis,datacomp)
      endif
!
! Calculate Mulliken charge
!
      call calcumulliken(dmtrxa,dmtrxb,smtrx,databasis,datacomp)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(datajob%octupole) then
        call memset(nao3*29,datacomp)
        allocate(work(nao3*29))
        call calcuoctupole(work,work(nao3*3+1),work(nao3*9+1),work(nao3*19+1),dmtrxa,dmtrxb, &
&                          datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*29,datacomp)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6,datacomp)
        allocate(work(nao3*6))
        call calcudipole(work,work(nao3*3+1),dmtrxa,dmtrxb,datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*6,datacomp)
      endif
!
! Write checkpoint file
!
      if(datacomp%master.and.(datajob%check /= '')) call writecheck(cmoa,cmob,dmtrxa,dmtrxb,energymoa,energymob,datajob,databasis)
!
! Unset arrays 1
!
      deallocate(h1mtrx,smtrx,tmtrx,cmoa,cmob,ortho, &
&                dmtrxa,dmtrxb,xint,energymoa,energymob)
      call memunset(nao3*5+nao2*3+nshell3+nao*2,datacomp)
      call tstamp(1,datacomp)
      return
end


!----------------------------------------------------------------------------------------
  subroutine calcrgeometry(converged,datajob,databasis,datacomp)
!----------------------------------------------------------------------------------------
!
! Driver of geometry optimization calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
      use modmolecule, only : nmo, natom, coord, coordold, enuc
      use modtype, only : typejob, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer,allocatable :: iredun(:)
      integer :: nao, nao2, nao3, nshell3, natom3, ii, iopt
      integer :: isizered, numbond, numangle, numtorsion, numredun, maxredun
      real(8), parameter :: third=0.3333333333333333D+00
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmo(:), ortho(:), dmtrx(:)
      real(8), allocatable :: xint(:), energymo(:)
      real(8), allocatable :: egrad(:), egradold(:), ehess(:)
      real(8), allocatable :: overinv(:), work(:,:)
      real(8), allocatable :: workv(:), coordredun(:), egradredun(:)
      real(8) :: egradmax, egradrms
      real(8) :: savedconv, savecutint2
      logical,intent(out) :: converged
      logical :: exceed
!
      nao= databasis%nao
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(databasis%nshell*(databasis%nshell+1))/2
      natom3= natom*3
      converged=.false.
!
! Calculate redundant coordinate
!
      if(.not.datajob%cartesian) then
        isizered= natom*4*10
        call memset(isizered,datacomp)
        allocate(iredun(isizered))
        do ii= 1,10
          call setredundantcoord(iredun,isizered,numbond,numangle,numtorsion,exceed,datacomp)
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
      call memset(nao3*4+nao2*2+nshell3+nao,datacomp)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmo(nao2),ortho(nao2),dmtrx(nao3), &
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
        if(iopt >= 2) call writegeom(datacomp)
!
! Calculate nuclear repulsion energy
!
        call nucenergy(datajob%threshatom,datacomp)
        if(datacomp%master) then
          write(*,'(" Nuclear repulsion energy =",f15.8," a.u.",/)') enuc
        endif
!
! Set work arrays 1
!
        call memset(nao2*2,datacomp)
        allocate(overinv(nao2),work(nao,nao))
!
! Calculate overlap and 1-electron integrals
!
        call oneei(h1mtrx,smtrx,tmtrx,work,datacomp%nproc2,datacomp%myrank2,datacomp%mpi_comm2,datajob,databasis)
!
! Calculate canonicalization and inverse overlap matrices
!
        call fullmtrx(smtrx,work,nao)
        call mtrxcanoninv(ortho,overinv,work,nao,nmo,datajob%threshover,datacomp)
!
! Calculate initial MOs
!
        if(iopt == 1) then
          call guessmo(cmo,cmo,overinv,h1mtrx,ortho,datajob,databasis,datacomp)
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
        if((datajob%method == 'HARTREE-FOCK').or.(datajob%method == 'MP2')) then
          call calcrhf(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo,datajob,databasis,datacomp)
          if(iopt == 1) call writeeigenvalue(energymo,energymo,1,datajob,datacomp)
          call tstamp(1,datacomp)
        elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
          if((iopt == 1).and.(datajob%guess == 'HUCKEL')) then
            savedconv= datajob%dconv
            savecutint2= datajob%cutint2
            datajob%dconv= max(datajob%dconv,1.0D-2)
            datajob%cutint2= max(datajob%cutint2,1.0D-9)
            call calcrhf(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo,datajob,databasis,datacomp)
            datajob%dconv= savedconv
            datajob%cutint2= savecutint2
            call tstamp(1,datacomp)
          endif
          call calcrdft(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo,datajob,databasis,datacomp)
          if(iopt == 1) call writeeigenvalue(energymo,energymo,1,datajob,datacomp)
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
        if(datajob%method == 'HARTREE-FOCK') then
          call calcgradrhf(cmo,energymo,xint,egrad,datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
          call tstamp(1,datacomp)
        elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
          call calcgradrdft(cmo,energymo,xint,egrad,datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
          call tstamp(1,datacomp)
        elseif(datajob%method == 'MP2') then
          call calcgradrmp2(cmo,energymo,xint,egrad,datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
          call tstamp(1,datacomp)
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
        if(datacomp%master) write(*,'(" Optimization Cycle",i4,4x,"Maximum gradient =",f11.6,4x, &
&                            "RMS gradient =",f11.6,/)') iopt,egradmax,egradrms
!
! Check convergence
!
        if((egradmax <= datajob%optconv).and.(egradrms <= datajob%optconv*third)) then
          if(datacomp%master) write(*,'(" Geometry converged.",/)')
          converged=.true.
          exit
        endif
!
! Write checkpoint file
!
        if(datacomp%master.and.(datajob%check /= '')) call writecheck(cmo,cmo,dmtrx,dmtrx,energymo,energymo,datajob,databasis)
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
          call calcnewcoord(coord,coordold,egrad,egradold,ehess,workv,natom3,iopt,datajob,datacomp)
        else
          call calcnewcoordred(coord,coordold,coordredun,egrad,egradredun,ehess,work(1,1), &
&                              work(1,2),work(1,3),work(1,4),workv,iopt,iredun,isizered, &
&                              maxredun,numbond,numangle,numtorsion,numredun,datajob,datacomp)
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
        call setnextopt(coordold,natom,iopt,datajob)
!
        if((iopt == datajob%nopt).and.datacomp%master) then
          write(*,'("Warning! Geometry did not converge.")')
          datacomp%nwarn= datacomp%nwarn+1
          exit
        endif
        call tstamp(1,datacomp)
      enddo
!
! End of optimization cycle 
!
      call writeeigenvalue(energymo,energymo,1,datajob,datacomp)
!
! Print MOs
!
      if(datacomp%master.and.(datajob%iprint >= 2)) then
        write(*,'("  -------------------")')
        write(*,'("    MO coefficients")')
        write(*,'("  -------------------")')
        call writeeigenvector(cmo,energymo,databasis,datacomp)
      endif
!
! Calculate Mulliken charge
!
      call calcrmulliken(dmtrx,smtrx,databasis,datacomp)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(datajob%octupole) then
        call memset(nao3*29,datacomp)
        allocate(work(nao3,29))
        call calcroctupole(work,work(1,4),work(1,10),work(1,20),dmtrx, &
&                          datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*29,datacomp)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6,datacomp)
        allocate(work(nao3,6))
        call calcrdipole(work,work(1,4),dmtrx,datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*6,datacomp)
      endif
!
! Write optimized geometry
!
      if(datacomp%master.and.converged) then
        write(*,'(" ==========================")')
        write(*,'("     Optimized Geometry")')
        write(*,'(" ==========================")')
        call writegeom(datacomp)
      endif
!
! Write checkpoint file
!
      if(datacomp%master.and.(datajob%check /= '')) call writecheck(cmo,cmo,dmtrx,dmtrx,energymo,energymo,datajob,databasis)
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
      deallocate(h1mtrx,smtrx,tmtrx,cmo,ortho,dmtrx, &
&                xint,energymo)
      call memunset(nao3*4+nao2*2+nshell3+nao,datacomp)
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
  subroutine calcugeometry(converged,datajob,databasis,datacomp)
!----------------------------------------------------------------------------------------
!
! Driver of open-shell geometry optimization calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
      use modmolecule, only : nmo, natom, coord, coordold, enuc
      use modtype, only : typejob, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer,allocatable :: iredun(:)
      integer :: nao, nao2, nao3, nshell3, natom3, ii, iopt
      integer :: isizered, numbond, numangle, numtorsion, numredun, maxredun
      real(8), parameter :: third=0.3333333333333333D+00
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmoa(:), cmob(:), ortho(:)
      real(8), allocatable :: dmtrxa(:), dmtrxb(:), xint(:), energymoa(:), energymob(:)
      real(8), allocatable :: egrad(:), egradold(:), ehess(:)
      real(8), allocatable :: overinv(:,:), work(:,:)
      real(8), allocatable :: workv(:), coordredun(:), egradredun(:)
      real(8) :: egradmax, egradrms
      real(8) :: savedconv, savecutint2
      logical,intent(out) :: converged
      logical :: exceed
!
      nao= databasis%nao
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(databasis%nshell*(databasis%nshell+1))/2
      natom3= natom*3
      converged=.false.
!
! Calculate redundant coordinate
!
      if(.not.datajob%cartesian) then
        isizered= natom*4*10
        call memset(isizered,datacomp)
        allocate(iredun(isizered))
        do ii= 1,10
          call setredundantcoord(iredun,isizered,numbond,numangle,numtorsion,exceed,datacomp)
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
      call memset(nao3*5+nao2*3+nshell3+nao*2,datacomp)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmoa(nao2),cmob(nao2),ortho(nao2), &
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
        if(iopt >= 2) call writegeom(datacomp)
!
! Calculate nuclear repulsion energy
!
        call nucenergy(datajob%threshatom,datacomp)
        if(datacomp%master) then
          write(*,'(" Nuclear repulsion energy =",f15.8," a.u.",/)') enuc
        endif
!
! Set arrays 1
!
        call memset(nao2*2,datacomp)
        allocate(overinv(nao,nao),work(nao,nao))
!
! Calculate overlap and 1-electron integrals
!
        call oneei(h1mtrx,smtrx,tmtrx,work,datacomp%nproc2,datacomp%myrank2,datacomp%mpi_comm2,datajob,databasis)
!
! Calculate canonicalization and inverse overlap matrices
!
        call fullmtrx(smtrx,work,nao)
        call mtrxcanoninv(ortho,overinv,work,nao,nmo,datajob%threshover,datacomp)
!
! Calculate initial MOs
!
        if(iopt == 1) then
          call guessmo(cmoa,cmob,overinv,h1mtrx,ortho,datajob,databasis,datacomp)
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
        if(datajob%method == 'HARTREE-FOCK') then
          call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob,datajob,databasis,datacomp)
          if(iopt == 1) call writeeigenvalue(energymoa,energymob,2,datajob,datacomp)
          call tstamp(1,datacomp)
        elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
          if((iopt == 1).and.(datajob%guess == 'HUCKEL')) then
            savedconv= datajob%dconv
            savecutint2= datajob%cutint2
            datajob%dconv= max(datajob%dconv,1.0D-2)
            datajob%cutint2= max(datajob%cutint2,1.0D-9)
            call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob,datajob,databasis,datacomp)
            datajob%dconv= savedconv
            datajob%cutint2= savecutint2
            call tstamp(1,datacomp)
          endif
          call calcudft(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob,datajob,databasis,datacomp)
          if(iopt == 1) call writeeigenvalue(energymoa,energymob,2,datajob,datacomp)
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
        if(datajob%method == 'HARTREE-FOCK') then
          call calcgraduhf(cmoa,cmob,energymoa,energymob,xint,egrad,datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
        elseif((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
          call calcgradudft(cmoa,cmob,energymoa,energymob,xint,egrad,datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
        else
          if(datacomp%master) then
            write(*,'(" Error! This program does not support method= ",a16," in energy gradient.")') &
&                 datajob%method
            call iabort
          endif
        endif
        call tstamp(1,datacomp)
!
! Calculate maximum and root mean square gradient values
!
        call calcmaxgrad(egradmax,egradrms,egrad,natom3)
        if(datacomp%master) write(*,'(" Optimization Cycle",i4,4x,"Maximum gradient =",f11.6,4x, &
&                            "RMS gradient =",f11.6,/)') iopt,egradmax,egradrms
!
! Write checkpoint file
!
        if(datacomp%master.and.(datajob%check /= '')) call writecheck(cmoa,cmob,dmtrxa,dmtrxb,energymoa,energymob,datajob,databasis)
!
! Check convergence
!
        if((egradmax <= datajob%optconv).and.(egradrms <= datajob%optconv*third)) then
          if(datacomp%master) write(*,'(" Geometry converged.",/)')
          converged=.true.
          exit
        endif
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
          call calcnewcoord(coord,coordold,egrad,egradold,ehess,workv,natom3,iopt,datajob,datacomp)
        else
          call calcnewcoordred(coord,coordold,coordredun,egrad,egradredun,ehess,work(1,1), &
&                              work(1,2),work(1,3),work(1,4),workv,iopt,iredun,isizered, &
&                              maxredun,numbond,numangle,numtorsion,numredun,datajob,datacomp)
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
        call setnextopt(coordold,natom,iopt,datajob)
!
        if((iopt == datajob%nopt).and.datacomp%master) then
          write(*,'("Warning! Geometry did not converge.")')
          datacomp%nwarn= datacomp%nwarn+1
          exit
        endif
        call tstamp(1,datacomp)
      enddo
!
! End of optimization cycle 
!
      call writeeigenvalue(energymoa,energymob,2,datajob,datacomp)
!
! Print MOs
!
      if(datacomp%master.and.(datajob%iprint >= 2)) then
        write(*,'("  -------------------------")')
        write(*,'("    Alpha MO coefficients")')
        write(*,'("  -------------------------")')
        call writeeigenvector(cmoa,energymoa,databasis,datacomp)
        write(*,'("  ------------------------")')
        write(*,'("    Beta MO coefficients")')
        write(*,'("  ------------------------")')
        call writeeigenvector(cmob,energymob,databasis,datacomp)
      endif
!
! Calculate Mulliken charge
!
      call calcumulliken(dmtrxa,dmtrxb,smtrx,databasis,datacomp)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(datajob%octupole) then
        call memset(nao3*29,datacomp)
        allocate(work(nao3,29))
        call calcuoctupole(work,work(1,4),work(1,10),work(1,20),dmtrxa,dmtrxb, &
&                          datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*29,datacomp)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6,datacomp)
        allocate(work(nao3,6))
        call calcudipole(work,work(1,4),dmtrxa,dmtrxb,datacomp%nproc1,datacomp%myrank1,datacomp%mpi_comm1,datajob,databasis,datacomp)
        deallocate(work)
        call memunset(nao3*6,datacomp)
      endif
!
! Write optimized geometry
!
      if(datacomp%master.and.converged) then
        write(*,'(" ==========================")')
        write(*,'("     Optimized Geometry")')
        write(*,'(" ==========================")')
        call writegeom(datacomp)
      endif
!
! Write checkpoint file
!
      if(datacomp%master.and.(datajob%check /= '')) call writecheck(cmoa,cmob,dmtrxa,dmtrxb,energymoa,energymob,datajob,databasis)
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
      deallocate(h1mtrx,smtrx,tmtrx,cmoa,cmob,ortho, &
&                dmtrxa,dmtrxb,xint,energymoa,energymob)
      call memunset(nao3*5+nao2*3+nshell3+nao*2,datacomp)
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


!---------------------------------------------
  subroutine setnextopt(coordold,natom,iopt,datajob)
!---------------------------------------------
!
! Set parameters for next optimization step
!
      use modguess, only : coord_g, nmo_g
      use modmolecule, only : nmo
      use modtype, only : typejob
      implicit none
      type(typejob),intent(inout) :: datajob
      integer,intent(in) :: natom, iopt
      integer :: iatom
      real(8),intent(in) :: coordold(3,natom)
!
! Set MO projection as initial MO calculation
!
      datajob%guess= 'UPDATE'
!
! Copy coordinate and energy gradient
!
      do iatom= 1,natom
        coord_g(1,iatom)= coordold(1,iatom)
        coord_g(2,iatom)= coordold(2,iatom)
        coord_g(3,iatom)= coordold(3,iatom)
      enddo
!
      nmo_g= nmo
!
      return
end


!--------------------
  subroutine setdft(datajob,databasis,datacomp)
!--------------------
!
! Set functional information
! Adjust the numbe of DFT grids when heavy elements are included
!
      use modmolecule, only : natom, numatomic, atomrad
      use modtype, only : typejob, typebasis, typecomp
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer :: ii, maxelem
!
      do ii= 1,9
        atomrad(-ii)= datajob%bqrad(ii)
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
        case('HARTREE-FOCK','MP2')
        case default
          if(datacomp%master) then
            write(*,'(" Error! This program does not support method= ",a16,".")') datajob%method
            call iabort
          endif
      endselect
!
      if((datajob%idftex >= 1).or.(datajob%idftcor >= 1)) then
        maxelem= maxval(numatomic(1:natom))
        if(((maxelem >= 55).or.(databasis%nao >= 2000)).and.((datajob%nrad == 96).and.(datajob%nleb == 302))) then
          datacomp%nwarn= datacomp%nwarn+1
          if(datacomp%master) write(*,'(" Warning! The number of DFT grids may not be enough.")')
        endif
      endif
!
      return
end


!--------------------
  subroutine setmp2(datajob,databasis)
!--------------------
!
! Set MP2 information
!
      use modtype, only : typejob, typebasis
      implicit none
      type(typejob),intent(inout) :: datajob
      type(typebasis),intent(in) :: databasis
      integer :: ncorecalc
!
      if(datajob%ncore == -1) datajob%ncore= ncorecalc(databasis)
!
      return
end
