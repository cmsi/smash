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
      use modparallel, only : master, nproc1, nproc2, myrank1, myrank2, mpi_comm1, mpi_comm2
      use modwarn, only : nwarn
      use modmemory, only : memusedmax
      use modjob, only : runtype, scftype
      use modiofile, only : input, icheck, check, version
      implicit none
      logical :: converged
!
      call setparallel
      version='2.2.0'
!
      if(master) then
        write(*,&
&           '(" *******************************************",/,&
&             "    Scalable Molecular Analysis Solver for",/,&
&             "      High performance computing systems",/,&
&             "            SMASH Version ",a10/,&
&             "           written by K. ISHIMURA",/,&
&             " *******************************************",/)') version
      endif
      call tstamp(0)
      call parallelinfo
!
! Read input file and set details
!
      call setdetails(mpi_comm1)
!
! Start calculations
!
      if(scftype == 'RHF') then
        if(runtype == 'ENERGY') then
          call calcrenergy(nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        elseif(runtype == 'GRADIENT') then
          call calcrgradient(nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        elseif(runtype == 'OPTIMIZE') then
          call calcrgeometry(converged,nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        else
          if(master) then
            write(*,'(" Error! This program does not support runtype= ",a16,".")')runtype
            call iabort
          endif
        endif
      elseif(scftype == 'UHF') then
        if(runtype == 'ENERGY') then
          call calcuenergy(nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        elseif(runtype == 'GRADIENT') then
          call calcugradient(nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        elseif(runtype == 'OPTIMIZE') then
          call calcugeometry(converged,nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        endif
      else
        if(master) write(*,'(" Error! SCFtype=",a16," is not supported.")')scftype
        call iabort
      endif
!
! Close input.dat and checkpoint files
!
      if(master) close(unit=input,status='DELETE')
      if(master.and.(check /= '')) close(unit=icheck)
!
      call para_finalize
      call memcheck
      call tstamp(2)
      if(master) then
        write(*,'(" Used memory :",1x,i6," MB")')memusedmax/125000
        if((runtype =='OPTIMIZE').and.(.not.converged))then
          write(*,'(/," ============================================================")')
          write(*,'("  Geometry optimization did not finish with",i3," warning(s)!")')nwarn
          write(*,'(" ============================================================")')
        else
          write(*,'(" Your calculation finished with",i3," warning(s).")')nwarn
        endif
      endif
end program main


!-------------------------
  subroutine setparallel
!-------------------------
!
!  Initialize MPI execution environment
!
      use modparallel, only : master, nproc1, nproc2, myrank1, myrank2, &
&                             mpi_comm1, mpi_comm2
      implicit none
!
! Initialize variables for parallelization
!
      master = .true.
!
! Start MPI parallelization and set mpi_comm1=MPI_COMM_WORLD
!
      call para_init(mpi_comm1)
      call para_comm_size(nproc1,mpi_comm1)
      call para_comm_rank(myrank1,mpi_comm1)
!
      nproc2= nproc1
      myrank2= myrank1
      mpi_comm2= mpi_comm1
!
      if(nproc1.gt.1) then
        master =(myrank1 == 0)
      endif
!
      return
end


!----------------------------------
  subroutine setdetails(mpi_comm)
!----------------------------------
!
! Read input file and set variables
!
      use modparallel, only : master
      use modecp, only : flagecp
      implicit none
      integer,intent(in) :: mpi_comm
!
! Set defaults before reading input file
!
      call setdefault1
!
! Read input data and open checkpoint file if necessary
!
      if(master) call opendatfile
      call readinput(mpi_comm)
!
! Set basis functions
!
      call setbasis(mpi_comm)
!
! Set ECP functions
!
      if(flagecp) call setecp(mpi_comm)
!
! Set maximum memory size
!
      call maxmemset
!
! Set number of electrons
!
      call setelectron
!
! Reset defaults after reading input file
!
      call setdefault2
!
! Set functional information and adjust the number of DFT grids
!
      call setdft
!
! Set functional information and adjust the number of DFT grids
!
      call setmp2
!
! Write input data
!
      call writecondition
      call writegeom
      call writebasis
      if(flagecp) call writeecp
!
! Set atom charge including dummy atom
!
      call setcharge(mpi_comm)
!
      return
end


!-------------------------
  subroutine setdefault1
!-------------------------
!
! Set defaults before reading input file
!
      use modiofile, only : check
      use modwarn, only : nwarn
      use modguess, only : spher_g, guess
      use modmemory, only : memmax, memused, memusedmax, memory
      use modprint, only : iprint
      use modunit, only : bohr
      use modbasis, only : spher, basis
      use modscf, only : maxiter, dconv, fdiff, scfconv, maxdiis, maxsoscf, maxqc, &
&                        maxqcdiag, maxqcdiagsub, extrap
      use modthresh, only : precision, cutint2, threshsoscf, threshqc, threshover, threshatom, &
&                           threshdiis, threshweight, threshrho, threshdfock, threshdftao, &
&                           threshmp2cphf
      use moddft, only : idftex, idftcor, nrad, nleb, bqrad
      use modopt, only : nopt, optconv, cartesian
      use modecp, only : ecp, flagecp
      use modjob, only : scftype, runtype, method
      use modmolecule, only : multi, charge
      use modmp2, only : ncore, nvfz, maxmp2diis, maxmp2iter
      use modprop, only : octupole
      implicit none
!
      nwarn  = 0
      memmax = 1000000000
      memused= 0
      memusedmax= 0
      memory = ''
      maxiter= 150
      maxdiis= 20
      maxsoscf= 20
      maxqc   = 15
      maxqcdiag= 100
      maxqcdiagsub= 10
      fdiff  =.true.
      scfconv='DIIS'
      extrap =.false.
      threshsoscf= 0.25D+00
      threshqc   = 1.0D-05
      threshover = 1.0D-06
      threshatom = 2.0D-01
      threshdiis = 6.0D-01
      threshmp2cphf=1.0D-10
      idftex = 0
      idftcor= 0
      iprint = 2
      bohr   =.false.
      spher  =.true.
      spher_g=.true.
      nopt   = 100
      optconv= 1.0D-04
      cartesian=.false.
      multi  = 1
      charge = 0.0D+00
      bqrad(:)=1.0D+00
      nvfz= 0
      maxmp2diis= 20
      maxmp2iter= 100
!
      cutint2=-1.0d+00
      nrad = 0
      nleb = 0
      ncore= -1
      dconv=-1.0D+00
      threshweight=-1.0D+00
      threshrho=-1.0D+00
      threshdfock=-1.0D+00
      threshdftao=-1.0D+00
!
      precision='MEDIUM'
      flagecp= .false.
      scftype='RHF'
      method='HARTREE-FOCK'
      runtype='ENERGY'
      basis='STO-3G'
      guess='HUCKEL'
      ecp=''
      check=''
      octupole=.false.
!
      return
end


!-------------------------
  subroutine setdefault2
!-------------------------
!
! Reset defaults after reading input file
!
      use modparallel, only : master
      use modthresh, only : precision, cutint2, threshweight, threshrho, threshdfock, threshdftao
      use moddft, only : nrad, nleb
      use modscf, only : dconv
      implicit none
      real(8),parameter :: zero= 0.0D+00

      select case(precision)
        case('HIGH')
          if(cutint2 < zero) cutint2= 1.0D-12
          if(dconv   < zero) dconv  = 5.0D-06
          if(threshweight < zero) threshweight=1.0D-08
          if(threshrho    < zero) threshrho   =1.0D-06
          if(threshdfock  < zero) threshdfock =1.0D-05
          if(threshdftao  < zero) threshdftao =1.0D-04
          if(nrad == 0) nrad= 150
          if(nleb == 0) nleb= 590
        case('MEDIUM')
          if(cutint2 < zero) cutint2= 1.0D-11
          if(dconv   < zero) dconv  = 5.0D-06
          if(threshweight < zero) threshweight=1.0D-08
          if(threshrho    < zero) threshrho   =1.0D-05
          if(threshdfock  < zero) threshdfock =1.0D-04
          if(threshdftao  < zero) threshdftao =1.0D-03
          if(nrad == 0) nrad= 96
          if(nleb == 0) nleb= 302
        case('LOW')
          if(cutint2 < zero) cutint2= 1.0D-10
          if(dconv   < zero) dconv  = 1.0D-05
          if(threshweight < zero) threshweight=1.0D-08
          if(threshrho    < zero) threshrho   =1.0D-04
          if(threshdfock  < zero) threshdfock =1.0D-04
          if(threshdftao  < zero) threshdftao =1.0D-02
          if(nrad == 0) nrad= 72
          if(nleb == 0) nleb= 302
        case default
          if(master) write(*,'(" Error! This program does not support precision= ", &
&                          a16,".")') precision
          call iabort
      end select
!
      return
end


!-------------------------
  subroutine setelectron
!-------------------------
!
! Set number of electrons
!
      use modparallel, only : master
      use modmolecule, only : numatomic, neleca, nelecb, natom, multi, charge
      use modecp, only : flagecp, izcore
      use modjob, only : scftype
      use modwarn, only : nwarn
      implicit none
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
      if(flagecp) then
        do ii= 1,natom
          nume= nume-izcore(ii)
        enddo
      endif

!
! Calculate numbers of alpha and beta electrons
!
      if((scftype == 'RHF').and.(multi /= 1)) then
        if(master) write(*,'(" Warning! SCFtype changes from RHF to UHF.")')
        scftype = 'UHF'
        nwarn= nwarn+1
      endif
!
      neleca=(nume+multi-1)/2
      nelecb=(nume-multi+1)/2
      if((neleca+nelecb)/= nume) then
        if(master) write(*,'(" Error! Spin multiplicity is ",i2, &
&                               ", but number of elctrons is ",i5,".")')multi, nume
        call iabort
      endif
!
      return
end


!----------------------------------------------------------------------------
  subroutine calcrenergy(nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!----------------------------------------------------------------------------
!
! Driver of closed-shell energy calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
      use modparallel, only : master
      use modiofile, only : check
      use modbasis, only : nao, nshell
      use modenergy, only : enuc
      use modmolecule, only : nmo
      use modjob, only : method
      use moddft, only : idftex, idftcor
      use modguess, only : guess
      use modprint, only : iprint
      use modscf, only : dconv
      use modthresh, only : cutint2, threshover
      use modprop, only : octupole
      implicit none
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2, mpi_comm1, mpi_comm2
      integer :: nao2, nao3, nshell3
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmo(:), ortho(:), dmtrx(:)
      real(8), allocatable :: xint(:), energymo(:)
      real(8), allocatable :: overinv(:), work(:)
      real(8) :: savedconv, savecutint2
!
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(nshell*(nshell+1))/2
!
! Set arrays 1
!
      call memset(nao3*4+nao2*2+nshell3+nao)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmo(nao2),ortho(nao2),dmtrx(nao3),&
&              xint(nshell3),energymo(nao))
!
! Calculate nuclear repulsion energy
!
      call nucenergy
      if(master) then
        write(*,'(" Nuclear repulsion energy =",f15.8," a.u.",/)') enuc
      endif
!
! Set arrays 2
!
      call memset(nao2*2)
      allocate(overinv(nao2),work(nao2))
!
! Calculate overlap and 1-electron integrals
!
      call oneei(h1mtrx,smtrx,tmtrx,work,nproc2,myrank2,mpi_comm2)
!
! Calculate canonicalization and inverse overlap matrices
!
      call fullmtrx(smtrx,work,nao)
      call mtrxcanoninv(ortho,overinv,work,nao,nmo,threshover,nproc2,myrank2,mpi_comm2)
!
! Calculate initial MOs
!
      call guessmo(cmo,cmo,overinv,h1mtrx,ortho, &
&                  nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!
! Unset arrays 2
!
      deallocate(overinv,work)
      call memunset(nao2*2)
      call tstamp(1)
!
! Start SCF
!
      if(method == 'HARTREE-FOCK') then
        call calcrhf(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo, &
&                    nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymo,energymo,1)
        call tstamp(1)
      elseif((idftex >= 1).or.(idftcor >= 1)) then
        if(guess == 'HUCKEL') then
          savedconv= dconv
          savecutint2= cutint2
          dconv= max(dconv,1.0D-2)
          cutint2= max(cutint2,1.0D-9)
          call calcrhf(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo, &
&                      nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
          dconv= savedconv
          cutint2= savecutint2
          call tstamp(1)
        endif
        call calcrdft(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo, &
&                     nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymo,energymo,1)
        call tstamp(1)
      elseif(method == 'MP2') then
        call calcrhf(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo, &
&                    nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymo,energymo,1)
        call tstamp(1)
        call calcrmp2(cmo,energymo,xint,nproc1,myrank1,mpi_comm1)
        call tstamp(1)
      else
        if(master) then
          write(*,'(" Error! This program does not support method= ",a16,".")')method
          call iabort
        endif
      endif
!
! Print MOs
!
      if(master.and.(iprint >= 2)) then
        write(*,'("  -------------------")')
        write(*,'("    MO coefficients")')
        write(*,'("  -------------------")')
        call writeeigenvector(cmo,energymo)
      endif
!
! Calculate Mulliken charge
!
      call calcrmulliken(dmtrx,smtrx)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(octupole) then
        call memset(nao3*29)
        allocate(work(nao3*29))
        call calcroctupole(work,work(nao3*3+1),work(nao3*9+1),work(nao3*19+1),dmtrx, &
&                          nproc1,myrank1,mpi_comm1)
        deallocate(work)
        call memunset(nao3*29)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6)
        allocate(work(nao3*6))
        call calcrdipole(work,work(nao3*3+1),dmtrx,nproc1,myrank1,mpi_comm1)
        deallocate(work)
        call memunset(nao3*6)
      endif
!
! Write checkpoint file
!
      if(master.and.(check /= '')) call writecheck(cmo,cmo,dmtrx,dmtrx,energymo,energymo)
!
! Unset arrays 1
!
      deallocate(h1mtrx,smtrx,tmtrx,cmo,ortho,dmtrx, &
&                xint,energymo)
      call memunset(nao3*4+nao2*2+nshell3+nao)
      call tstamp(1)
      return
end


!----------------------------------------------------------------------------
  subroutine calcuenergy(nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!----------------------------------------------------------------------------
!
! Driver of open-shell energy calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
      use modparallel, only : master
      use modiofile, only : check
      use modbasis, only : nao, nshell
      use modenergy, only : enuc
      use modmolecule, only : nmo
      use modjob, only : method
      use moddft, only : idftex, idftcor
      use modguess, only : guess
      use modprint, only : iprint
      use modscf, only : dconv
      use modthresh, only : cutint2, threshover
      use modprop, only : octupole
      implicit none
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2, mpi_comm1, mpi_comm2
      integer :: nao2, nao3, nshell3
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmoa(:), cmob(:), ortho(:)
      real(8), allocatable :: dmtrxa(:), dmtrxb(:), xint(:), energymoa(:), energymob(:)
      real(8), allocatable :: overinv(:), work(:)
      real(8) :: savedconv, savecutint2
!
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(nshell*(nshell+1))/2
!
! Set arrays 1
!
      call memset(nao3*5+nao2*3+nshell3+nao*2)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmoa(nao2),cmob(nao2),ortho(nao2),&
&              dmtrxa(nao3),dmtrxb(nao3),xint(nshell3),energymoa(nao),energymob(nao))
!
! Calculate nuclear repulsion energy
!
      call nucenergy
      if(master) then
        write(*,'(" Nuclear repulsion energy =",f15.8," a.u.",/)') enuc
      endif
!
! Set arrays 2
!
      call memset(nao2*2)
      allocate(overinv(nao2),work(nao2))
!
! Calculate overlap and 1-electron integrals
!
      call oneei(h1mtrx,smtrx,tmtrx,work,nproc2,myrank2,mpi_comm2)
!
! Calculate canonicalization and inverse overlap matrices
!
      call fullmtrx(smtrx,work,nao)
      call mtrxcanoninv(ortho,overinv,work,nao,nmo,threshover,nproc2,myrank2,mpi_comm2)
!
! Calculate initial MOs
!
      call guessmo(cmoa,cmob,overinv,h1mtrx,ortho, &
&                  nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!
! Unset arrays 2
!
      deallocate(overinv,work)
      call memunset(nao2*2)
      call tstamp(1)
!
! Start SCF
!
      if(method == 'HARTREE-FOCK') then
        call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                    nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymoa,energymob,2)
        call tstamp(1)
      elseif((idftex >= 1).or.(idftcor >= 1)) then
        if(guess == 'HUCKEL') then
          savedconv= dconv
          savecutint2= cutint2
          dconv= max(dconv,1.0D-2)
          cutint2= max(cutint2,1.0D-9)
          call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                      nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
          dconv= savedconv
          cutint2= savecutint2
          call tstamp(1)
        endif
        call calcudft(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                     nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymoa,energymob,2)
        call tstamp(1)
!     elseif(method == 'MP2') then
!       call calcuhf(h1mtrx,cmoa,ortho,smtrx,xint,energymoa)
!       call tstamp(1)
!       call calcump2(cmoa,energymoa,xint)
!       call tstamp(1)
      else
        if(master) then
          write(*,'(" Error! This program does not support method= ",a16,".")')method
          call iabort
        endif
      endif
!
! Print MOs
!
      if(master.and.(iprint >= 2)) then
        write(*,'("  -------------------------")')
        write(*,'("    Alpha MO coefficients")')
        write(*,'("  -------------------------")')
        call writeeigenvector(cmoa,energymoa)
        write(*,'("  ------------------------")')
        write(*,'("    Beta MO coefficients")')
        write(*,'("  ------------------------")')
        call writeeigenvector(cmob,energymob)
      endif
!
! Calculate Mulliken charge
!
      call calcumulliken(dmtrxa,dmtrxb,smtrx)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(octupole) then
        call memset(nao3*29)
        allocate(work(nao3*29))
        call calcuoctupole(work,work(nao3*3+1),work(nao3*9+1),work(nao3*19+1),dmtrxa,dmtrxb, &
&                          nproc1,myrank1,mpi_comm1)
        deallocate(work)
        call memunset(nao3*29)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6)
        allocate(work(nao3*6))
        call calcudipole(work,work(nao3*3+1),dmtrxa,dmtrxb,nproc1,myrank1,mpi_comm1)
        deallocate(work)
        call memunset(nao3*6)
      endif
!
! Write checkpoint file
!
      if(master.and.(check /= '')) call writecheck(cmoa,cmob,dmtrxa,dmtrxb,energymoa,energymob)
!
! Unset arrays 1
!
      deallocate(h1mtrx,smtrx,tmtrx,cmoa,cmob,ortho, &
&                dmtrxa,dmtrxb,xint,energymoa,energymob)
      call memunset(nao3*5+nao2*3+nshell3+nao*2)
      return
end


!------------------------------------------------------------------------------
  subroutine calcrgradient(nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!------------------------------------------------------------------------------
!
! Driver of energy gradient calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
      use modparallel, only : master
      use modiofile, only : check
      use modbasis, only : nao, nshell
      use modenergy, only : enuc
      use modmolecule, only : nmo, natom
      use modjob, only : method
      use moddft, only : idftex, idftcor
      use modguess, only : guess
      use modprint, only : iprint
      use modscf, only : dconv
      use modthresh, only : cutint2, threshover
      use modprop, only : octupole
      implicit none
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2, mpi_comm1, mpi_comm2
      integer :: nao2, nao3, nshell3
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmo(:), ortho(:), dmtrx(:)
      real(8), allocatable :: xint(:), energymo(:)
      real(8), allocatable :: overinv(:), work(:)
      real(8), allocatable :: egrad(:)
      real(8) :: egradmax, egradrms
      real(8) :: savedconv, savecutint2
!
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(nshell*(nshell+1))/2
!
! Set arrays 1
!
      call memset(nao3*4+nao2*2+nshell3+nao)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmo(nao2),ortho(nao2),dmtrx(nao3),&
&              xint(nshell3),energymo(nao))
!
! Calculate nuclear repulsion energy
!
      call nucenergy
      if(master) then
        write(*,'(" Nuclear repulsion energy =",f15.8," a.u.",/)') enuc
      endif
!
! Set arrays 2
!
      call memset(nao2*2)
      allocate(overinv(nao2),work(nao2))
!
! Calculate overlap and 1-electron integrals
!
      call oneei(h1mtrx,smtrx,tmtrx,work,nproc2,myrank2,mpi_comm2)
!
! Calculate canonicalization and inverse overlap matrices
!
      call fullmtrx(smtrx,work,nao)
      call mtrxcanoninv(ortho,overinv,work,nao,nmo,threshover,nproc2,myrank2,mpi_comm2)
!
! Calculate initial MOs
!
      call guessmo(cmo,cmo,overinv,h1mtrx,ortho, &
&                  nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!
! Unset arrays 2
!
      deallocate(overinv,work)
      call memunset(nao2*2)
      call tstamp(1)
!
! Start SCF
!
      if((method == 'HARTREE-FOCK').or.(method == 'MP2')) then
        call calcrhf(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo, &
&                    nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymo,energymo,1)
        call tstamp(1)
      elseif((idftex >= 1).or.(idftcor >= 1)) then
        if(guess == 'HUCKEL') then
          savedconv= dconv
          savecutint2= cutint2
          dconv= max(dconv,1.0D-2)
          cutint2= max(cutint2,1.0D-9)
          call calcrhf(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo, &
&                      nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
          dconv= savedconv
          cutint2= savecutint2
          call tstamp(1)
        endif
        call calcrdft(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo, &
&                     nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymo,energymo,1)
        call tstamp(1)
      else
        if(master) then
          write(*,'(" Error! This program does not support method= ",a16,".")')method
          call iabort
        endif
      endif
!
! Set arrays 3
!
      call memset(natom*3)
      allocate(egrad(natom*3))
!
! Calculate energy gradient
!
      if(method == 'HARTREE-FOCK') then
        call calcgradrhf(cmo,energymo,xint,egrad,nproc1,myrank1,mpi_comm1)
        call tstamp(1)
      elseif((idftex >= 1).or.(idftcor >= 1)) then
        call calcgradrdft(cmo,energymo,xint,egrad,nproc1,myrank1,mpi_comm1)
        call tstamp(1)
      elseif(method == 'MP2') then
        call calcgradrmp2(cmo,energymo,xint,egrad,nproc1,myrank1,mpi_comm1)
        call tstamp(1)
      else
        if(master) then
          write(*,'(" Error! This program does not support method= ",a16,".")')method
          call iabort
        endif
      endif
!
! Calculate maximum and root mean square gradient values
!
      call calcmaxgrad(egradmax,egradrms,egrad,natom*3)
      if(master) write(*,'(" Maximum gradient =",f13.8,"  RMS gradient =",f13.8,/)') &
&                      egradmax,egradrms
!
! Unset arrays 3
!
      deallocate(egrad)
      call memunset(natom*3)
!
! Print MOs
!
      if(master.and.(iprint >= 2)) then
        write(*,'("  -------------------")')
        write(*,'("    MO coefficients")')
        write(*,'("  -------------------")')
        call writeeigenvector(cmo,energymo)
      endif
!
! Calculate Mulliken charge
!
      call calcrmulliken(dmtrx,smtrx)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(octupole) then
        call memset(nao3*29)
        allocate(work(nao3*29))
        call calcroctupole(work,work(nao3*3+1),work(nao3*9+1),work(nao3*19+1),dmtrx, &
&                          nproc1,myrank1,mpi_comm1)
        deallocate(work)
        call memunset(nao3*29)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6)
        allocate(work(nao3*6))
        call calcrdipole(work,work(nao3*3+1),dmtrx,nproc1,myrank1,mpi_comm1)
        deallocate(work)
        call memunset(nao3*6)
      endif
!
! Write checkpoint file
!
      if(master.and.(check /= '')) call writecheck(cmo,cmo,dmtrx,dmtrx,energymo,energymo)
!
! Unset arrays 1
!
      deallocate(h1mtrx,smtrx,tmtrx,cmo,ortho,dmtrx, &
&                xint,energymo)
      call memunset(nao3*4+nao2*2+nshell3+nao)
      call tstamp(1)
      return
end


!------------------------------------------------------------------------------
  subroutine calcugradient(nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!------------------------------------------------------------------------------
!
! Driver of open-shell energy gradient calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
      use modparallel, only : master
      use modiofile, only : check
      use modbasis, only : nao, nshell
      use modenergy, only : enuc
      use modmolecule, only : nmo, natom
      use modjob, only : method
      use moddft, only : idftex, idftcor
      use modguess, only : guess
      use modprint, only : iprint
      use modscf, only : dconv
      use modthresh, only : cutint2, threshover
      use modprop, only : octupole
      implicit none
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2, mpi_comm1, mpi_comm2
      integer :: nao2, nao3, nshell3
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmoa(:), cmob(:), ortho(:)
      real(8), allocatable :: dmtrxa(:), dmtrxb(:), xint(:), energymoa(:), energymob(:)
      real(8), allocatable :: overinv(:), work(:)
      real(8), allocatable :: egrad(:)
      real(8) :: egradmax, egradrms
      real(8) :: savedconv, savecutint2
!
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(nshell*(nshell+1))/2
!
! Set arrays 1
!
      call memset(nao3*5+nao2*3+nshell3+nao*2)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmoa(nao2),cmob(nao2),ortho(nao2),&
&              dmtrxa(nao3),dmtrxb(nao3),xint(nshell3),energymoa(nao),energymob(nao))
!
! Calculate nuclear repulsion energy
!
      call nucenergy
      if(master) then
        write(*,'(" Nuclear repulsion energy =",f15.8," a.u.",/)') enuc
      endif
!
! Set arrays 2
!
      call memset(nao2*2)
      allocate(overinv(nao2),work(nao2))
!
! Calculate overlap and 1-electron integrals
!
      call oneei(h1mtrx,smtrx,tmtrx,work,nproc2,myrank2,mpi_comm2)
!
! Calculate canonicalization and inverse overlap matrices
!
      call fullmtrx(smtrx,work,nao)
      call mtrxcanoninv(ortho,overinv,work,nao,nmo,threshover,nproc2,myrank2,mpi_comm2)
!
! Calculate initial MOs
!
      call guessmo(cmoa,cmob,overinv,h1mtrx,ortho, &
&                  nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!
! Unset arrays 2
!
      deallocate(overinv,work)
      call memunset(nao2*2)
      call tstamp(1)
!
! Start SCF
!
      if(method == 'HARTREE-FOCK') then
        call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                    nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymoa,energymob,2)
        call tstamp(1)
      elseif((idftex >= 1).or.(idftcor >= 1)) then
        if(guess == 'HUCKEL') then
          savedconv= dconv
          savecutint2= cutint2
          dconv= max(dconv,1.0D-2)
          cutint2= max(cutint2,1.0D-9)
          call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                      nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
          dconv= savedconv
          cutint2= savecutint2
          call tstamp(1)
        endif
        call calcudft(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                     nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymoa,energymob,2)
        call tstamp(1)
!     elseif(method == 'MP2') then
!       call calcuhf(h1mtrx,cmoa,ortho,smtrx,xint,energymoa)
!       call tstamp(1)
!       call calcump2(cmoa,energymoa,xint)
!       call tstamp(1)
      else
        if(master) then
          write(*,'(" Error! This program does not support method= ",a16,".")')method
          call iabort
        endif
      endif
!
! Set arrays 3
!
      call memset(natom*3)
      allocate(egrad(natom*3))
!
! Calculate energy gradient
!
      if(method == 'HARTREE-FOCK') then
        call calcgraduhf(cmoa,cmob,energymoa,energymob,xint,egrad,nproc1,myrank1,mpi_comm1)
      elseif((idftex >= 1).or.(idftcor >= 1)) then
        call calcgradudft(cmoa,cmob,energymoa,energymob,xint,egrad,nproc1,myrank1,mpi_comm1)
      else
        if(master) then
          write(*,'(" Error! This program does not support method= ",a16," in energy gradient.")') &
&               method
          call iabort
        endif
      endif
!
! Calculate maximum and root mean square gradient values
!
      call calcmaxgrad(egradmax,egradrms,egrad,natom*3)
      if(master) write(*,'(" Maximum gradient =",f13.8,"  RMS gradient =",f13.8,/)') &
&                      egradmax,egradrms
!
! Unset arrays 3
!
      deallocate(egrad)
      call memunset(natom*3)
!
! Print MOs
!
      if(master.and.(iprint >= 2)) then
        write(*,'("  -------------------------")')
        write(*,'("    Alpha MO coefficients")')
        write(*,'("  -------------------------")')
        call writeeigenvector(cmoa,energymoa)
        write(*,'("  ------------------------")')
        write(*,'("    Beta MO coefficients")')
        write(*,'("  ------------------------")')
        call writeeigenvector(cmob,energymob)
      endif
!
! Calculate Mulliken charge
!
      call calcumulliken(dmtrxa,dmtrxb,smtrx)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(octupole) then
        call memset(nao3*29)
        allocate(work(nao3*29))
        call calcuoctupole(work,work(nao3*3+1),work(nao3*9+1),work(nao3*19+1),dmtrxa,dmtrxb, &
&                          nproc1,myrank1,mpi_comm1)
        deallocate(work)
        call memunset(nao3*29)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6)
        allocate(work(nao3*6))
        call calcudipole(work,work(nao3*3+1),dmtrxa,dmtrxb,nproc1,myrank1,mpi_comm1)
        deallocate(work)
        call memunset(nao3*6)
      endif
!
! Write checkpoint file
!
      if(master.and.(check /= '')) call writecheck(cmoa,cmob,dmtrxa,dmtrxb,energymoa,energymob)
!
! Unset arrays 1
!
      deallocate(h1mtrx,smtrx,tmtrx,cmoa,cmob,ortho, &
&                dmtrxa,dmtrxb,xint,energymoa,energymob)
      call memunset(nao3*5+nao2*3+nshell3+nao*2)
      call tstamp(1)
      return
end


!----------------------------------------------------------------------------------------
  subroutine calcrgeometry(converged,nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!----------------------------------------------------------------------------------------
!
! Driver of geometry optimization calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
      use modparallel, only : master
      use modiofile, only : check
      use modbasis, only : nao, nshell
      use modenergy, only : enuc
      use modmolecule, only : nmo, natom, coord, coordold
      use modopt, only : nopt, optconv, cartesian
      use modwarn, only : nwarn
      use modjob, only : method
      use moddft, only : idftex, idftcor
      use modguess, only : guess
      use modprint, only : iprint
      use modscf, only : dconv
      use modthresh, only : cutint2, threshover
      use modprop, only : octupole
      implicit none
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2, mpi_comm1, mpi_comm2
      integer,allocatable :: iredun(:)
      integer :: nao2, nao3, nshell3, natom3, ii, iopt
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
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(nshell*(nshell+1))/2
      natom3= natom*3
      converged=.false.
!
! Calculate redundant coordinate
!
      if(.not.cartesian) then
        isizered= natom*4*10
        call memset(isizered)
        allocate(iredun(isizered))
        do ii= 1,10
          call setredundantcoord(iredun,isizered,numbond,numangle,numtorsion,exceed)
          if(.not.exceed) exit
          call memunset(isizered)
          deallocate(iredun)
          isizered= isizered*2
          call memset(isizered)
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
      call memset(nao3*4+nao2*2+nshell3+nao)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmo(nao2),ortho(nao2),dmtrx(nao3), &
&              xint(nshell3),energymo(nao))
!
! Set arrays for energy gradient and geometry optimization
!
      if(cartesian) then
        call memset(natom3*2+natom3*(natom3+1)/2)
        allocate(egrad(natom3),egradold(natom3),ehess(natom3*(natom3+1)/2))
      else
        call memset(natom3+numredun*4+numredun*(numredun+1)/2)
        allocate(egrad(natom3),coordredun(numredun*2),egradredun(numredun*2), &
&                ehess(numredun*(numredun+1)/2))
      endif
!
! Start geometry optimization cycle
!
      do iopt= 1,nopt
!
! Print geometry
!
        if(iopt >= 2) call writegeom
!
! Calculate nuclear repulsion energy
!
        call nucenergy
        if(master) then
          write(*,'(" Nuclear repulsion energy =",f15.8," a.u.",/)') enuc
        endif
!
! Set work arrays 1
!
        call memset(nao2*2)
        allocate(overinv(nao2),work(nao,nao))
!
! Calculate overlap and 1-electron integrals
!
        call oneei(h1mtrx,smtrx,tmtrx,work,nproc2,myrank2,mpi_comm2)
!
! Calculate canonicalization and inverse overlap matrices
!
        call fullmtrx(smtrx,work,nao)
        call mtrxcanoninv(ortho,overinv,work,nao,nmo,threshover,nproc2,myrank2,mpi_comm2)
!
! Calculate initial MOs
!
        if(iopt == 1) then
          call guessmo(cmo,cmo,overinv,h1mtrx,ortho, &
&                      nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        endif
!
! Unset work arrays 1
!
        deallocate(overinv,work)
        call memunset(nao2*2)
        call tstamp(1)
!
! Calculate energy
!
        if((method == 'HARTREE-FOCK').or.(method == 'MP2')) then
          call calcrhf(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo, &
&                      nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
          if(iopt == 1) call writeeigenvalue(energymo,energymo,1)
          call tstamp(1)
        elseif((idftex >= 1).or.(idftcor >= 1)) then
          if((iopt == 1).and.(guess == 'HUCKEL')) then
            savedconv= dconv
            savecutint2= cutint2
            dconv= max(dconv,1.0D-2)
            cutint2= max(cutint2,1.0D-9)
            call calcrhf(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo, &
&                        nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
            dconv= savedconv
            cutint2= savecutint2
            call tstamp(1)
          endif
          call calcrdft(h1mtrx,cmo,ortho,smtrx,dmtrx,xint,energymo, &
&                       nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
          if(iopt == 1) call writeeigenvalue(energymo,energymo,1)
          call tstamp(1)
        else
          if(master) then
            write(*,'(" Error! This program does not support method= ",a16,".")')method
            call iabort
          endif
        endif
!
! Calculate energy gradient
!
        if(method == 'HARTREE-FOCK') then
          call calcgradrhf(cmo,energymo,xint,egrad,nproc1,myrank1,mpi_comm1)
          call tstamp(1)
        elseif((idftex >= 1).or.(idftcor >= 1)) then
          call calcgradrdft(cmo,energymo,xint,egrad,nproc1,myrank1,mpi_comm1)
          call tstamp(1)
        elseif(method == 'MP2') then
          call calcgradrmp2(cmo,energymo,xint,egrad,nproc1,myrank1,mpi_comm1)
          call tstamp(1)
        else
          if(master) then
            write(*,'(" Error! This program does not support method= ",a16,".")')method
            call iabort
          endif
        endif
!
! Calculate maximum and root mean square gradient values
!
        call calcmaxgrad(egradmax,egradrms,egrad,natom3)
        if(master) write(*,'(" Optimization Cycle",i4,4x,"Maximum gradient =",f11.6,4x, &
&                            "RMS gradient =",f11.6,/)') iopt,egradmax,egradrms
!
! Check convergence
!
        if((egradmax <= optconv).and.(egradrms <= optconv*third)) then
          if(master) write(*,'(" Geometry converged.",/)')
          converged=.true.
          exit
        endif
!
! Write checkpoint file
!
        if(master.and.(check /= '')) call writecheck(cmo,cmo,dmtrx,dmtrx,energymo,energymo)
!
! Set work arrays 2
!
        if(cartesian) then
          call memset(natom3*3)
          allocate(workv(natom3*3))
        else
          call memset(maxredun*maxredun*4+maxredun*3)
          allocate(work(maxredun*maxredun,4),workv(maxredun*3))
        endif
!
! Calculate new coordinate
!
        if(cartesian) then
          call calcnewcoord(coord,coordold,egrad,egradold,ehess,workv,natom3,iopt, &
&                           nproc2,myrank2,mpi_comm2)
        else
          call calcnewcoordred(coord,coordold,coordredun,egrad,egradredun,ehess,work(1,1), &
&                              work(1,2),work(1,3),work(1,4),workv,iopt,iredun,isizered, &
&                              maxredun,numbond,numangle,numtorsion,numredun, &
&                              nproc2,myrank2,mpi_comm2)
        endif
!
! Unset work arrays 2
!
        if(cartesian) then
          deallocate(workv)
          call memunset(natom3*3)
        else
          deallocate(work,workv)
          call memunset(maxredun*maxredun*4+maxredun*3)
        endif
!
! Set guess MO calculation flag from Huckel to projection
!
        call setnextopt(coordold,natom,iopt)
!
        if((iopt == nopt).and.master) then
          write(*,'("Warning! Geometry did not converge.")')
          nwarn= nwarn+1
          exit
        endif
        call tstamp(1)
      enddo
!
! End of optimization cycle 
!
      call writeeigenvalue(energymo,energymo,1)
!
! Print MOs
!
      if(master.and.(iprint >= 2)) then
        write(*,'("  -------------------")')
        write(*,'("    MO coefficients")')
        write(*,'("  -------------------")')
        call writeeigenvector(cmo,energymo)
      endif
!
! Calculate Mulliken charge
!
      call calcrmulliken(dmtrx,smtrx)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(octupole) then
        call memset(nao3*29)
        allocate(work(nao3,29))
        call calcroctupole(work,work(1,4),work(1,10),work(1,20),dmtrx, &
&                          nproc1,myrank1,mpi_comm1)
        deallocate(work)
        call memunset(nao3*29)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6)
        allocate(work(nao3,6))
        call calcrdipole(work,work(1,4),dmtrx,nproc1,myrank1,mpi_comm1)
        deallocate(work)
        call memunset(nao3*6)
      endif
!
! Write optimized geometry
!
      if(master.and.converged) then
        write(*,'(" ==========================")')
        write(*,'("     Optimized Geometry")')
        write(*,'(" ==========================")')
        call writegeom
      endif
!
! Write checkpoint file
!
      if(master.and.(check /= '')) call writecheck(cmo,cmo,dmtrx,dmtrx,energymo,energymo)
!
! Unset arrays for energy gradient and geometry optimization
!
      if(cartesian) then
        deallocate(egrad,egradold,ehess)
        call memunset(natom3*2+natom3*(natom3+1)/2)
      else
        deallocate(egrad,coordredun,egradredun, &
&                  ehess)
        call memunset(natom3+numredun*4+numredun*(numredun+1)/2)
      endif
!
! Unset arrays for energy
!
      deallocate(h1mtrx,smtrx,tmtrx,cmo,ortho,dmtrx, &
&                xint,energymo)
      call memunset(nao3*4+nao2*2+nshell3+nao)
!
! Unset array for redundant coordinate
!
      if(.not.cartesian) then
        deallocate(iredun)
        call memunset(isizered)
      endif
!
      call tstamp(1)
      return
end


!----------------------------------------------------------------------------------------
  subroutine calcugeometry(converged,nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!----------------------------------------------------------------------------------------
!
! Driver of open-shell geometry optimization calculation
!
! Parallel information
!   nproc1, myrank1, mpi_comm1 : MPI_COMM_WORLD (all nodes)
!   nproc2, myrank2, mpi_comm2 : new communicator for matrix operations
!                               (default: MPI_COMM_WORLD)
!
      use modparallel, only : master
      use modiofile, only : check
      use modbasis, only : nao, nshell
      use modenergy, only : enuc
      use modmolecule, only : nmo, natom, coord, coordold
      use modopt, only : nopt, optconv, cartesian
      use modwarn, only : nwarn
      use modjob, only : method
      use moddft, only : idftex, idftcor
      use modguess, only : guess
      use modprint, only : iprint
      use modscf, only : dconv
      use modthresh, only : cutint2, threshover
      use modprop, only : octupole
      implicit none
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2, mpi_comm1, mpi_comm2
      integer,allocatable :: iredun(:)
      integer :: nao2, nao3, nshell3, natom3, ii, iopt
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
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(nshell*(nshell+1))/2
      natom3= natom*3
      converged=.false.
!
! Calculate redundant coordinate
!
      if(.not.cartesian) then
        isizered= natom*4*10
        call memset(isizered)
        allocate(iredun(isizered))
        do ii= 1,10
          call setredundantcoord(iredun,isizered,numbond,numangle,numtorsion,exceed)
          if(.not.exceed) exit
          call memunset(isizered)
          deallocate(iredun)
          isizered= isizered*2
          call memset(isizered)
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
      call memset(nao3*5+nao2*3+nshell3+nao*2)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmoa(nao2),cmob(nao2),ortho(nao2), &
&              dmtrxa(nao3),dmtrxb(nao3),xint(nshell3),energymoa(nao),energymob(nao))
!
! Set arrays for energy gradient and geometry optimization
!
      if(cartesian) then
        call memset(natom3*2+natom3*(natom3+1)/2)
        allocate(egrad(natom3),egradold(natom3),ehess(natom3*(natom3+1)/2))
      else
        call memset(natom3+numredun*4+numredun*(numredun+1)/2)
        allocate(egrad(natom3),coordredun(numredun*2),egradredun(numredun*2), &
&                ehess(numredun*(numredun+1)/2))
      endif
!
! Start geometry optimization cycle
!
      do iopt= 1,nopt
!
! Print geometry
!
        if(iopt >= 2) call writegeom
!
! Calculate nuclear repulsion energy
!
        call nucenergy
        if(master) then
          write(*,'(" Nuclear repulsion energy =",f15.8," a.u.",/)') enuc
        endif
!
! Set arrays 1
!
        call memset(nao2*2)
        allocate(overinv(nao,nao),work(nao,nao))
!
! Calculate overlap and 1-electron integrals
!
        call oneei(h1mtrx,smtrx,tmtrx,work,nproc2,myrank2,mpi_comm2)
!
! Calculate canonicalization and inverse overlap matrices
!
        call fullmtrx(smtrx,work,nao)
        call mtrxcanoninv(ortho,overinv,work,nao,nmo,threshover,nproc2,myrank2,mpi_comm2)
!
! Calculate initial MOs
!
        if(iopt == 1) then
          call guessmo(cmoa,cmob,overinv,h1mtrx,ortho, &
&                      nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        endif
!
! Unset arrays 1
!
        deallocate(overinv,work)
        call memunset(nao2*2)
        call tstamp(1)
!
! Calculate energy
!
        if(method == 'HARTREE-FOCK') then
          call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                      nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
          if(iopt == 1) call writeeigenvalue(energymoa,energymob,2)
          call tstamp(1)
        elseif((idftex >= 1).or.(idftcor >= 1)) then
          if((iopt == 1).and.(guess == 'HUCKEL')) then
            savedconv= dconv
            savecutint2= cutint2
            dconv= max(dconv,1.0D-2)
            cutint2= max(cutint2,1.0D-9)
            call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                        nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
            dconv= savedconv
            cutint2= savecutint2
            call tstamp(1)
          endif
          call calcudft(h1mtrx,cmoa,cmob,ortho,smtrx,dmtrxa,dmtrxb,xint,energymoa,energymob, &
&                       nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
          if(iopt == 1) call writeeigenvalue(energymoa,energymob,2)
          call tstamp(1)
!       elseif(method == 'MP2') then
!         call calcuhf(h1mtrx,cmoa,ortho,smtrx,xint,energymoa)
!         call tstamp(1)
!         call calcump2(cmoa,energymoa,xint)
!         call tstamp(1)
        else
          if(master) then
            write(*,'(" Error! This program does not support method= ",a16,".")')method
            call iabort
          endif
        endif
!
! Calculate energy gradient
!
        if(method == 'HARTREE-FOCK') then
          call calcgraduhf(cmoa,cmob,energymoa,energymob,xint,egrad,nproc1,myrank1,mpi_comm1)
        elseif((idftex >= 1).or.(idftcor >= 1)) then
          call calcgradudft(cmoa,cmob,energymoa,energymob,xint,egrad,nproc1,myrank1,mpi_comm1)
        else
          if(master) then
            write(*,'(" Error! This program does not support method= ",a16," in energy gradient.")') &
&                 method
            call iabort
          endif
        endif
        call tstamp(1)
!
! Calculate maximum and root mean square gradient values
!
        call calcmaxgrad(egradmax,egradrms,egrad,natom3)
        if(master) write(*,'(" Optimization Cycle",i4,4x,"Maximum gradient =",f11.6,4x, &
&                            "RMS gradient =",f11.6,/)') iopt,egradmax,egradrms
!
! Write checkpoint file
!
        if(master.and.(check /= '')) call writecheck(cmoa,cmob,dmtrxa,dmtrxb,energymoa,energymob)
!
! Check convergence
!
        if((egradmax <= optconv).and.(egradrms <= optconv*third)) then
          if(master) write(*,'(" Geometry converged.",/)')
          converged=.true.
          exit
        endif
!
! Set work arrays 2
!
        if(cartesian) then
          call memset(natom3*3)
          allocate(workv(natom3*3))
        else
          call memset(maxredun*maxredun*4+maxredun*3)
          allocate(work(maxredun*maxredun,4),workv(maxredun*3))
        endif
!
! Calculate new coordinate
!
        if(cartesian) then
          call calcnewcoord(coord,coordold,egrad,egradold,ehess,workv,natom3,iopt, &
&                           nproc2,myrank2,mpi_comm2)
        else
          call calcnewcoordred(coord,coordold,coordredun,egrad,egradredun,ehess,work(1,1), &
&                              work(1,2),work(1,3),work(1,4),workv,iopt,iredun,isizered, &
&                              maxredun,numbond,numangle,numtorsion,numredun, &
&                              nproc2,myrank2,mpi_comm2)
        endif
!
! Unset work arrays 2
!
        if(cartesian) then
          deallocate(workv)
          call memunset(natom3*3)
        else
          deallocate(work,workv)
          call memunset(maxredun*maxredun*4+maxredun*3)
        endif
!
! Set guess MO calculation flag from Huckel to projection
!
        call setnextopt(coordold,natom,iopt)
!
        if((iopt == nopt).and.master) then
          write(*,'("Warning! Geometry did not converge.")')
          nwarn= nwarn+1
          exit
        endif
        call tstamp(1)
      enddo
!
! End of optimization cycle 
!
      call writeeigenvalue(energymoa,energymob,2)
!
! Print MOs
!
      if(master.and.(iprint >= 2)) then
        write(*,'("  -------------------------")')
        write(*,'("    Alpha MO coefficients")')
        write(*,'("  -------------------------")')
        call writeeigenvector(cmoa,energymoa)
        write(*,'("  ------------------------")')
        write(*,'("    Beta MO coefficients")')
        write(*,'("  ------------------------")')
        call writeeigenvector(cmob,energymob)
      endif
!
! Calculate Mulliken charge
!
      call calcumulliken(dmtrxa,dmtrxb,smtrx)
!
! Calculate dipole, quadrupole, and octupole moments
!
      if(octupole) then
        call memset(nao3*29)
        allocate(work(nao3,29))
        call calcuoctupole(work,work(1,4),work(1,10),work(1,20),dmtrxa,dmtrxb, &
&                          nproc1,myrank1,mpi_comm1)
        deallocate(work)
        call memunset(nao3*29)
      else
!
! Calculate dipole moment
!
        call memset(nao3*6)
        allocate(work(nao3,6))
        call calcudipole(work,work(1,4),dmtrxa,dmtrxb,nproc1,myrank1,mpi_comm1)
        deallocate(work)
        call memunset(nao3*6)
      endif
!
! Write optimized geometry
!
      if(master.and.converged) then
        write(*,'(" ==========================")')
        write(*,'("     Optimized Geometry")')
        write(*,'(" ==========================")')
        call writegeom
      endif
!
! Write checkpoint file
!
      if(master.and.(check /= '')) call writecheck(cmoa,cmob,dmtrxa,dmtrxb,energymoa,energymob)
!
! Unset arrays for energy gradient and geometry optimization
!
      if(cartesian) then
        deallocate(egrad,egradold,ehess)
        call memunset(natom3*2+natom3*(natom3+1)/2)
      else
        deallocate(egrad,coordredun,egradredun, &
&                ehess)
        call memunset(natom3+numredun*4+numredun*(numredun+1)/2)
      endif
!
! Unset arrays for energy
!
      deallocate(h1mtrx,smtrx,tmtrx,cmoa,cmob,ortho, &
&                dmtrxa,dmtrxb,xint,energymoa,energymob)
      call memunset(nao3*5+nao2*3+nshell3+nao*2)
      call tstamp(1)
!
! Unset array for redundant coordinate
!
      if(.not.cartesian) then
        deallocate(iredun)
        call memunset(isizered)
      endif
!
      return
end


!---------------------------------------------
  subroutine setnextopt(coordold,natom,iopt)
!---------------------------------------------
!
! Set parameters for next optimization step
!
      use modbasis, only : ex, coeff, nshell, nao, nprim, locprim, locbf, &
&                          locatom, mprim, mbf, mtype, spher
      use modguess, only : ex_g, coeff_g, nshell_g, nao_g, nmo_g, nprim_g, locprim_g, locbf_g, &
&                          locatom_g, mprim_g, mbf_g, mtype_g, spher_g, coord_g, guess
      use modmolecule, only : nmo
      implicit none
      integer,intent(in) :: natom, iopt
      integer :: iprim, ishell, iatom
      real(8),intent(in) :: coordold(3,natom)
!
! Set MO projection as initial MO calculation
!
      guess= 'UPDATE'
!
! Copy coordinate and energy gradient
!
      do iatom= 1,natom
        coord_g(1,iatom)= coordold(1,iatom)
        coord_g(2,iatom)= coordold(2,iatom)
        coord_g(3,iatom)= coordold(3,iatom)
      enddo
!
! Copy basis set information
!
      nmo_g= nmo
!
      if(iopt == 1) then
        nao_g= nao
        nprim_g= nprim
        nshell_g= nshell
!
        do iprim= 1,nprim
          ex_g(iprim)= ex(iprim)
          coeff_g(iprim)= coeff(iprim)
        enddo
!
        do ishell= 1,nshell
          locprim_g(ishell)= locprim(ishell)
          locbf_g(ishell)  = locbf(ishell)
          locatom_g(ishell)= locatom(ishell)
          mprim_g(ishell)  = mprim(ishell)
          mbf_g(ishell)    = mbf(ishell)
          mtype_g(ishell)  = mtype(ishell)
        enddo
!
        spher_g= spher
      endif
!
      return
end


!--------------------
  subroutine setdft
!--------------------
!
! Set functional information
! Adjust the numbe of DFT grids when heavy elements are included
!
      use modparallel, only : master
      use moddft, only : idftex, idftcor, nrad, nleb, hfexchange, bqrad
      use modatom, only : atomrad
      use modmolecule, only : natom, numatomic
      use modjob, only : method
      use modwarn, only : nwarn
      use modbasis, only : nao
      implicit none
      integer :: ii, maxelem
!
      do ii= 1,9
        atomrad(-ii)= bqrad(ii)
      enddo
!
      select case(method)
        case('B3LYP')
          idftex = 1
          idftcor= 1
          hfexchange= 0.2D+00
        case('B3LYP5')
          idftex = 1
          idftcor= 2
          hfexchange= 0.2D+00
        case('HARTREE-FOCK','MP2')
        case default
          if(master) then
            write(*,'(" Error! This program does not support method= ",a16,".")') method
            call iabort
          endif
      endselect
!
      if((idftex >= 1).or.(idftcor >= 1)) then
        maxelem= maxval(numatomic(1:natom))
        if(((maxelem >= 55).or.(nao >= 2000)).and.((nrad == 96).and.(nleb == 302))) then
          nwarn= nwarn+1
          if(master) write(*,'(" Warning! The number of DFT grids may not be enough.")')
        endif
      endif
!
      return
end


!--------------------
  subroutine setmp2
!--------------------
!
! Set MP2 information
!
      use modmp2, only : ncore
      implicit none
      integer :: ncorecalc
!
      if(ncore == -1) ncore= ncorecalc()
!
      return
end
