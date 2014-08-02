!---------------
  program main
!---------------
!
! This is the main driver of Scalable Molecular Analysis Solver 
! for High performance computing (SMASH).
!
      use modparallel
      use modiofile, only : in
      use modwarn, only : nwarn
      use modmemory, only : memusedmax
      use modjob, only : runtype, method, scftype
      use modecp, only : flagecp
      implicit none
      logical :: converged
!
!
      call start
!
      if(master) then
        write(*,&
&           '(" ****************************************",/,&
&             "    Scalable Molecular Analysis Solver",/,&
&             "      for High performance computing",/,&
&             "            SMASH Version 0.1",/,&
&             "          written by K. ISHIMURA",/,&
&             " ****************************************",/)')
      endif
      call tstamp(0)
      call gethostnm
      call parallelinfo
!
! Read input data
!
      if(master) open(unit=in,file='input.dat',status='replace')
      call readinput(mpi_comm1)
!
! Set maximum memory size
!
      call maxmemset
!
! Set basis functions
!
      call setbasis(mpi_comm1)
!
! Set number of electrons
!
      call setelectron
!
! Set ECP functions
!
      if(flagecp) call setecp(mpi_comm1)
!
! Set functional information and adjust the number of DFT grids
!
      call setdft
!
! Write input data
!
      call writecondition
      call writegeom
      call writebasis
      if(flagecp) call writeecp
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
            write(*,'(" Error! This program does not support ",a16,".")')runtype
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
      if(master) close(in)
      call para_finalize
      call memcheck
      call tstamp(2)
      if(master) then
        write(*,'(" Used memory :",1x,i6," MB")')memusedmax/125000
        if((runtype =='OPTIMIZE').and.(.not.converged))then
          write(*,'(" Geometry optimization did not finish with",i3," warning(s).")')nwarn
        else
          write(*,'(" Your calculation finished with",i3," warning(s).")')nwarn
        endif
      endif
end program main


!-------------------
  subroutine start
!-------------------
!
! Set computational data and machine information
!
      use modparallel
      use modwarn, only : nwarn
      use modguess, only : iguess, spher_g, guess
      use modmemory, only : memmax, memused, memusedmax, memory
      use modprint, only : iprint
      use modunit, only : bohr
      use modbasis, only : spher, basis
      use modscf, only : maxiter, dconv, fdiff, diis, maxdiis, maxsoscf, extrap
      use modthresh, only : cutint2, threshsoscf
      use moddft, only : idft, nrad, nleb
      use modopt, only : nopt, optconv, cartesian
      use modecp, only : ecp, flagecp
      use modjob, only : scftype, runtype, method
      use modmolecule, only : multi, charge
      implicit none
!
! Initialize valuables for parallelization
!
      nproc  = 1
      myrank = 0
      master = .true.
      parallel = .false.
!
      mpi_comm1= MPI_COMM_WORLD
      call para_init
      call para_comm_size(nproc,mpi_comm1)
      call para_comm_rank(myrank,mpi_comm1)
!
      mpi_comm2= mpi_comm2
      nproc1= nproc
      nproc2= nproc
      myrank1= myrank
      myrank2= myrank
!
      if(nproc.gt.1) then
        master =(myrank == 0)
        parallel = .true.
      endif
!
      nwarn  = 0
      iguess = 1
      memmax = 250000000
      memused= 0
      memusedmax= 0
      memory = ''
      maxiter= 100
      maxdiis= 20
      maxsoscf= 20
      fdiff  =.true.
      diis =.true.
      extrap =.false.
      cutint2= 1.0d-12
      threshsoscf= 0.25d0
      dconv  = 5.0d-6
      idft   = 0
      nrad   = 96
      nleb   = 302
      iprint = 1
      bohr   =.false.
      spher  =.true.
      spher_g=.true.
      nopt   = 50
      optconv= 1.0D-04
      cartesian=.false.
      multi  = 1
      charge = 0.0D+00
!
      flagecp= .false.
      scftype='RHF'
      method='HF'
      runtype='ENERGY'
      basis='STO-3G'
      guess='HUCKEL'
      ecp=''
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
      use modjob, only : scftype
      use modbasis, only : nao
      use modwarn, only : nwarn
      implicit none
      integer :: nume, ii
!
! Calculate total number of electrons
!
      nume= -charge
      do ii= 1,natom
        nume= nume+numatomic(ii)
      enddo
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
!      if(neleca > nao) then
!        if(master) write(*,'(" Error! The number of electrons is larger than the number ",&
!&                            "of basis functions.")')
!        call iabort
!      endif
      return
end


!----------------------------------------------------------------------------
  subroutine calcrenergy(nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!----------------------------------------------------------------------------
!
! Driver of closed-shell energy calculation
!
    
      use modbasis, only : nao, nshell
      use modenergy, only : enuc
      use modparallel, only : master
      use modmolecule, only : nmo, neleca
      use modjob, only : method
      use moddft, only : idft
      implicit none
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2
      integer(4),intent(in) :: mpi_comm1, mpi_comm2
      integer :: nao2, nao3, nshell3
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmo(:), ortho(:)
      real(8), allocatable :: xint(:), energymo(:)
      real(8), allocatable :: overinv(:), work(:)
!
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(nshell*(nshell+1))/2
!
! Set arrays 1
!
      call memset(nao3*3+nao2*2+nshell3+nao)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmo(nao2),ortho(nao2),&
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
      call mtrxcanoninv(ortho,overinv,work,nao,nmo,nproc2,myrank2,mpi_comm2)
!
! Calculate initial MOs
!
      call guessmo(cmo,overinv,nproc2,myrank2,mpi_comm2)
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
        call calcrhf(h1mtrx,cmo,ortho,smtrx,xint,energymo, &
&                    nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymo,energymo,1)
        call tstamp(1)
      elseif(idft >= 1) then
        call calcrhf(h1mtrx,cmo,ortho,smtrx,xint,energymo, &
&                    nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call tstamp(1)
        call calcrdft(h1mtrx,cmo,ortho,smtrx,xint,energymo, &
&                     nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymo,energymo,1)
        call tstamp(1)
      elseif(method == 'MP2') then
        call calcrhf(h1mtrx,cmo,ortho,smtrx,xint,energymo, &
&                    nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymo,energymo,1)
        call tstamp(1)
        call calcrmp2(cmo,energymo,xint,nproc1,myrank1,mpi_comm1)
        call tstamp(1)
      else
        if(master) then
          write(*,'(" Error! This program does not support ",a16,".")')method
          call iabort
        endif
      endif
!
! Print MOs
!
      if(master) then
        write(*,'("-----------------")')
        write(*,'(" MO coefficients")')
        write(*,'("-----------------")')
      endif
      call writeeigenvector(cmo,energymo)
!
! Unset arrays 1
!
      deallocate(h1mtrx,smtrx,tmtrx,cmo,ortho, &
&                xint,energymo)
      call memunset(nao3*3+nao2*2+nshell3+nao)
      return
end


!----------------------------------------------------------------------------
  subroutine calcuenergy(nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!----------------------------------------------------------------------------
!
! Driver of open-shell energy calculation
!
      use modparallel, only : master
      use modbasis, only : nao, nshell
      use modenergy, only : enuc
      use modmolecule, only : nmo
      use modjob, only : method
      use moddft, only : idft
      implicit none
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2
      integer(4),intent(in) :: mpi_comm1, mpi_comm2
      integer :: nao2, nao3, nshell3
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmoa(:), cmob(:), ortho(:)
      real(8), allocatable :: xint(:), energymoa(:), energymob(:)
      real(8), allocatable :: overinv(:), work(:)
!
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(nshell*(nshell+1))/2
!
! Set arrays 1
!
      call memset(nao3*3+nao2*3+nshell3+nao*2)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmoa(nao2),cmob(nao2),ortho(nao2),&
&              xint(nshell3),energymoa(nao),energymob(nao))
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
      call mtrxcanoninv(ortho,overinv,work,nao,nmo,nproc2,myrank2,mpi_comm2)
!
! Calculate initial MOs
!
      call guessmo(cmoa,overinv,nproc2,myrank2,mpi_comm2)
      cmob(:)= cmoa(:)
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
        call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,xint,energymoa,energymob, &
&                    nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymoa,energymob,2)
        call tstamp(1)
      elseif(idft >= 1) then
        call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,xint,energymoa,energymob, &
&                    nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call tstamp(1)
        call calcudft(h1mtrx,cmoa,cmob,ortho,smtrx,xint,energymoa,energymob, &
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
          write(*,'(" Error! This program does not support ",a16,".")')method
          call iabort
        endif
      endif
!
! Print MOs
!
      if(master) then
        write(*,'("-----------------------")')
        write(*,'(" Alpha MO coefficients")')
        write(*,'("-----------------------")')
      endif
      call writeeigenvector(cmoa,energymoa)
      if(master) then
        write(*,'("-----------------------")')
        write(*,'(" Beta MO coefficients")')
        write(*,'("-----------------------")')
      endif
      call writeeigenvector(cmob,energymob)
!
! Unset arrays 1
!
      deallocate(h1mtrx,smtrx,tmtrx,cmoa,cmob,ortho, &
&                xint,energymoa,energymob)
      call memunset(nao3*3+nao2*3+nshell3+nao*2)
      return
end


!------------------------------------------------------------------------------
  subroutine calcrgradient(nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!------------------------------------------------------------------------------
!
! Driver of energy gradient calculation
!
      use modparallel, only : master
      use modbasis, only : nao, nshell
      use modenergy, only : enuc
      use modmolecule, only : nmo, natom
      use modjob, only : method
      use moddft, only : idft
      implicit none
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2
      integer(4),intent(in) :: mpi_comm1, mpi_comm2
      integer :: nao2, nao3, nshell3
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmo(:), ortho(:)
      real(8), allocatable :: xint(:), energymo(:)
      real(8), allocatable :: overinv(:), work(:)
      real(8), allocatable :: egrad(:)
      real(8) :: egradmax, egradrms
!
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(nshell*(nshell+1))/2
!
! Set arrays 1
!
      call memset(nao3*3+nao2*2+nshell3+nao)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmo(nao2),ortho(nao2),&
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
      call mtrxcanoninv(ortho,overinv,work,nao,nmo,nproc2,myrank2,mpi_comm2)
!
! Calculate initial MOs
!
      call guessmo(cmo,overinv,nproc2,myrank2,mpi_comm2)
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
        call calcrhf(h1mtrx,cmo,ortho,smtrx,xint,energymo, &
&                    nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymo,energymo,1)
        call tstamp(1)
      elseif(idft >= 1) then
        call calcrhf(h1mtrx,cmo,ortho,smtrx,xint,energymo, &
&                    nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call tstamp(1)
        call calcrdft(h1mtrx,cmo,ortho,smtrx,xint,energymo, &
&                     nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymo,energymo,1)
        call tstamp(1)
      else
        if(master) then
          write(*,'(" Error! This program does not support ",a16,".")')method
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
      elseif(idft >= 1) then
        call calcgradrdft(cmo,energymo,xint,egrad,nproc1,myrank1,mpi_comm1)
      else
        if(master) then
          write(*,'(" Error! This program does not support ",a16,".")')method
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
      if(master) then
        write(*,'("-----------------")')
        write(*,'(" MO coefficients")')
        write(*,'("-----------------")')
      endif
      call writeeigenvector(cmo,energymo)
!
! Unset arrays 1
!
      deallocate(h1mtrx,smtrx,tmtrx,cmo,ortho, &
&                xint,energymo)
      call memunset(nao3*3+nao2*2+nshell3+nao)
      call tstamp(1)
      return
end


!------------------------------------------------------------------------------
  subroutine calcugradient(nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!------------------------------------------------------------------------------
!
! Driver of open-shell energy gradient calculation
!
      use modparallel, only : master
      use modbasis, only : nao, nshell
      use modenergy, only : enuc
      use modmolecule, only : nmo, natom
      use modjob, only : method
      use moddft, only : idft
      implicit none
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2
      integer(4),intent(in) :: mpi_comm1, mpi_comm2
      integer :: nao2, nao3, nshell3
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmoa(:), cmob(:), ortho(:)
      real(8), allocatable :: xint(:), energymoa(:), energymob(:)
      real(8), allocatable :: overinv(:), work(:)
      real(8), allocatable :: egrad(:)
      real(8) :: egradmax, egradrms
!
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(nshell*(nshell+1))/2
!
! Set arrays 1
!
      call memset(nao3*3+nao2*3+nshell3+nao*2)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmoa(nao2),cmob(nao2),ortho(nao2),&
&              xint(nshell3),energymoa(nao),energymob(nao))
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
      call mtrxcanoninv(ortho,overinv,work,nao,nmo,nproc2,myrank2,,mpi_comm2)
!
! Calculate initial MOs
!
      call guessmo(cmoa,overinv,nproc2,myrank2,mpi_comm2)
      cmob(:)= cmoa(:)
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
        call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,xint,energymoa,energymob, &
&                    nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call writeeigenvalue(energymoa,energymob,2)
        call tstamp(1)
      elseif(idft >= 1) then
        call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,xint,energymoa,energymob, &
&                    nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
        call tstamp(1)
        call calcudft(h1mtrx,cmoa,cmob,ortho,smtrx,xint,energymoa,energymob, &
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
          write(*,'(" Error! This program does not support ",a16,".")')method
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
      elseif(idft >= 1) then
        call calcgradudft(cmoa,cmob,energymoa,energymob,xint,egrad,nproc1,myrank1,mpi_comm1)
      else
        if(master) then
          write(*,'(" Error! This program does not support ",a16," in energy gradient.")') &
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
      if(master) then
        write(*,'("-----------------------")')
        write(*,'(" Alpha MO coefficients")')
        write(*,'("-----------------------")')
      endif
      call writeeigenvector(cmoa,energymoa)
      if(master) then
        write(*,'("-----------------------")')
        write(*,'(" Beta MO coefficients")')
        write(*,'("-----------------------")')
      endif
      call writeeigenvector(cmob,energymob)
!
! Unset arrays 1
!
      deallocate(h1mtrx,smtrx,tmtrx,cmoa,cmob,ortho, &
&                xint,energymoa,energymob)
      call memunset(nao3*3+nao2*3+nshell3+nao*2)
      call tstamp(1)
      return
end



!----------------------------------------------------------------------------------------
  subroutine calcrgeometry(converged,nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
!----------------------------------------------------------------------------------------
!
! Driver of geometry optimization calculation
!
      use modparallel, only : master
      use modbasis, only : nao, nshell
      use modenergy, only : enuc
      use modmolecule, only : nmo, natom, coord, coordold
      use modopt, only : nopt, optconv, cartesian
      use modwarn, only : nwarn
      use modjob, only : method
      use moddft, only : idft
      implicit none
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2
      integer(4),intent(in) :: mpi_comm1, mpi_comm2
      integer,allocatable :: iredun(:)
      integer :: nao2, nao3, nshell3, natom3, ii, iopt
      integer :: isizered, numbond, numangle, numtorsion, numredun, maxredun
      real(8), parameter :: third=0.3333333333333333D+00
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmo(:), ortho(:)
      real(8), allocatable :: xint(:), energymo(:)
      real(8), allocatable :: egrad(:), egradold(:), ehess(:)
      real(8), allocatable :: overinv(:), work(:,:)
      real(8), allocatable :: workv(:), coordredun(:), egradredun(:)
      real(8) :: egradmax, egradrms
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
      call memset(nao3*3+nao2*2+nshell3+nao)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmo(nao2),ortho(nao2), &
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
        call mtrxcanoninv(ortho,overinv,work,nao,nmo,nproc2,myrank2,mpi_comm2)
!
! Calculate initial MOs
!
!       call guessmo(cmo,overinv,nproc2,myrank2,mpi_comm2)
        if(iopt == 1) call guessmo(cmo,overinv,nproc2,myrank2,mpi_comm2)
!
! Unset work arrays 1
!
        deallocate(overinv,work)
        call memunset(nao2*2)
        call tstamp(1)
!
! Calculate energy
!
        if(method == 'HARTREE-FOCK') then
          call calcrhf(h1mtrx,cmo,ortho,smtrx,xint,energymo, &
&                      nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
          if(iopt == 1) call writeeigenvalue(energymo,energymo,1)
          call tstamp(1)
        elseif(idft >= 1) then
          if(iopt == 1) then
            call calcrhf(h1mtrx,cmo,ortho,smtrx,xint,energymo, &
&                        nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
            call tstamp(1)
          endif
          call calcrdft(h1mtrx,cmo,ortho,smtrx,xint,energymo, &
&                       nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
          if(iopt == 1) call writeeigenvalue(energymo,energymo,1)
          call tstamp(1)
        else
          if(master) then
            write(*,'(" Error! This program does not support ",a16,".")')method
            call iabort
          endif
        endif
!
! Calculate energy gradient
!
        if(method == 'HARTREE-FOCK') then
          call calcgradrhf(cmo,energymo,xint,egrad,nproc1,myrank1,mpi_comm1)
          call tstamp(1)
        elseif(idft >= 1) then
          call calcgradrdft(cmo,energymo,xint,egrad,nproc1,myrank1,mpi_comm1)
          call tstamp(1)
        else
          if(master) then
            write(*,'(" Error! This program does not support ",a16,".")')method
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
          if(master) write(*,'(" Geometry is converged.",/)')
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
          write(*,'("Warning! Geometry is not converged.")')
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
      if(master) then
        write(*,'("-----------------")')
        write(*,'(" MO coefficients")')
        write(*,'("-----------------")')
      endif
      call writeeigenvector(cmo,energymo)
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
      deallocate(h1mtrx,smtrx,tmtrx,cmo,ortho, &
&                xint,energymo)
      call memunset(nao3*3+nao2*2+nshell3+nao)
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
      use modparallel, only : master
      use modbasis, only : nao, nshell
      use modenergy, only : enuc
      use modmolecule, only : nmo, natom, coord, coordold
      use modopt, only : nopt, optconv, cartesian
      use modwarn, only : nwarn
      use modguess, only : iguess
      use modjob, only : method
      use moddft, only : idft
      implicit none
      integer,intent(in) :: nproc1, nproc2, myrank1, myrank2
      integer(4),intent(in) :: mpi_comm1, mpi_comm2
      integer,allocatable :: iredun(:)
      integer :: nao2, nao3, nshell3, natom3, ii, iopt
      integer :: isizered, numbond, numangle, numtorsion, numredun, maxredun
      real(8), parameter :: third=0.3333333333333333D+00
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmoa(:), cmob(:), ortho(:)
      real(8), allocatable :: xint(:), energymoa(:), energymob(:)
      real(8), allocatable :: egrad(:), egradold(:), ehess(:)
      real(8), allocatable :: overinv(:,:), work(:,:)
      real(8), allocatable :: workv(:), coordredun(:), egradredun(:)
      real(8) :: egradmax, egradrms
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
      call memset(nao3*3+nao2*3+nshell3+nao*2)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmoa(nao2),cmob(nao2),ortho(nao2), &
&              xint(nshell3),energymoa(nao),energymob(nao))
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
        call mtrxcanoninv(ortho,overinv,work,nao,nmo,nproc2,myrank2,mpi_comm2)
!
! Calculate initial MOs
!
        if(iguess == 1) then
          call guessmo(cmoa,overinv,nproc2,myrank2,mpi_comm2)
          cmob(:)= cmoa(:)
!       elseif(iguess == 2) then
!         work(:,:)= overinv(:,:)
!         call guessmo(cmoa,overinv,nproc2,myrank2,mpi_comm2)
!         call guessmo(cmob,work,nproc2,myrank2,mpi_comm2)
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
          call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,xint,energymoa,energymob, &
&                      nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
          if(iopt == 1) call writeeigenvalue(energymoa,energymob,2)
          call tstamp(1)
        elseif(idft >= 1) then
          if(iopt == 1) then
            call calcuhf(h1mtrx,cmoa,cmob,ortho,smtrx,xint,energymoa,energymob, &
&                        nproc1,nproc2,myrank1,myrank2,mpi_comm1,mpi_comm2)
            call tstamp(1)
          endif
          call calcudft(h1mtrx,cmoa,cmob,ortho,smtrx,xint,energymoa,energymob, &
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
            write(*,'(" Error! This program does not support ",a16,".")')method
            call iabort
          endif
        endif
!
! Calculate energy gradient
!
        if(method == 'HARTREE-FOCK') then
          call calcgraduhf(cmoa,cmob,energymoa,energymob,xint,egrad,nproc1,myrank1,mpi_comm1)
        elseif(idft >= 1) then
          call calcgradudft(cmoa,cmob,energymoa,energymob,xint,egrad,nproc1,myrank1,mpi_comm1)
        else
          if(master) then
            write(*,'(" Error! This program does not support ",a16," in energy gradient.")') &
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
! Check convergence
!
        if((egradmax <= optconv).and.(egradrms <= optconv*third)) then
          if(master) write(*,'(" Geometry is converged.",/)')
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
          write(*,'("Warning! Geometry is not converged.")')
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
      if(master) then
        write(*,'("-----------------------")')
        write(*,'(" Alpha MO coefficients")')
        write(*,'("-----------------------")')
      endif
      call writeeigenvector(cmoa,energymoa)
      if(master) then
        write(*,'("-----------------------")')
        write(*,'(" Beta MO coefficients")')
        write(*,'("-----------------------")')
      endif
      call writeeigenvector(cmob,energymob)
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
&                xint,energymoa,energymob)
      call memunset(nao3*3+nao2*3+nshell3+nao*2)
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
&                          locatom_g, mprim_g, mbf_g, mtype_g, spher_g, coord_g, iguess
      use modmolecule, only : nmo
      implicit none
      integer,intent(in) :: natom, iopt
      integer :: iprim, ishell, iatom
      real(8),intent(in) :: coordold(3,natom)
!
! Set MO projection as initial MO calculation
!
      iguess= 2
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
      use moddft, only : idft, nrad, nleb, hfexchange
      use modmolecule, only : natom, numatomic
      use modjob, only : method
      use modwarn, only : nwarn
      implicit none
      integer :: maxelem
!
      select case(method)
        case('B3LYP')
          idft= 1
          hfexchange= 0.2D+00
        case('B3PW91')
          idft= 2
          hfexchange= 0.2D+00
      endselect
!
      if(idft >= 1) then
        maxelem= maxval(numatomic(1:natom))
        if((maxelem >= 55).and.(nrad == 96).and.(nleb == 302)) then
          nrad= 150
          nleb= 590
          nwarn= nwarn+1
          if(master) then
            write(*,'(" Warning! The number of DFT grids is changed.")')
          endif
        endif
      endif
!
      return
end

