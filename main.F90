!---------------
  program main
!---------------
!
! This is a highly parallelized quantum chemistry calculation
! program package.
!
      use procpar, only : master
      use iofile, only : in, iout
      use warn, only : nwarn
      use memory, only : memusedmax
      use control, only : runtype, method
      implicit none
!
!
      call start
!
      if(master) then
        write(iout,&
&           '(" ***************************************************",/,&
&             "    Quantum chemistry calculation program package",   /,&
&             "                   Revision 0.1",                     /,&
&             "               written by K. ISHIMURA",               /,&
&             " ***************************************************")')
      endif
      call tstamp(0)
      call gethostnm
      call parallelinfo
!
! Read input data
!
      if(master) open(unit=in,file='input.dat',status='replace')
      call readinput
!
! Set basis functions
!
      call setbasis
!
! Write input data
!
      call writegeom
      call writebasis
      call writecondition
!
! Start energy calculation
!
      if(runtype == "ENERGY") then
        call calcenergy
      elseif(runtype == "GRADIENT") then
        call calcgradient
      elseif(runtype == "OPTIMIZE") then
        if(method == "HF") then
          call calcgeometry
        else
          if(master) write(iout,'(" Error! Optimization supports only HF method.")')
          call iabort
        endif
      endif
!
!ISHIMURA
!     if(master) close(in,status='delete')
      if(master) close(in)
      call para_finalize
      call memcheck
      call tstamp(2)
      if(master) then
        write(iout,'(" Used memory :",1x,i6," MB")')memusedmax/125000
        write(iout,'(" Your calculation finished with",i3," warning(s).")')nwarn
      endif
end program main


!-------------------
  subroutine start
!-------------------
!
! Set computational data and machine information
!
      use procpar, only : master, parallel, nproc, myrank
      use warn, only : nwarn
      use guess, only : iguess, spher_g
      use memory, only : memmax, memused, memusedmax
      use print, only : iprint
      use units, only : bohr
      use basis, only : spher
      use scf, only : maxiter, dconv, fdiff, dodiis, maxdiis, maxsoscf
      use thresh, only : cutint2, threshsoscf
      use dft, only : nrad, nleb
      use opt, only : nopt, optconv
      use ecp, only : flagecp
      implicit none
!
! Initialize valuables for parallelization
!
      nproc  = 1
      myrank = 0
      master = .true.
      parallel = .false.
!
      call para_init
      call para_comm_size(nproc)
      call para_comm_rank(myrank)
!
      if(nproc.gt.1) then
        master =(myrank == 0)
        parallel = .true.
      endif
!
      nwarn  = 0
      iguess = 1
      memmax = 2000000000
      memused= 0
      memusedmax= 0
      maxiter= 100
      maxdiis= 20
      maxsoscf= 20
      fdiff  =.true.
      dodiis =.true.
      dodiis =.false.
      cutint2= 5.0d-11
      threshsoscf= 0.25d0
!     dconv  = 1.0d-5    ! for energy calculation
      nrad   = 96
      nleb   = 302
!     nrad   = 24
!     nleb   = 110
      dconv  = 5.0d-6    ! for energy gradient calculation
!ishimura
      iprint = 2
      bohr   =.false.
      spher  =.false.
!     spher  =.true.
!     spher_g=.false.
      spher_g=.true.
      nopt   = 50
      optconv= 1.0D-04
!
      flagecp= .false.
      return
end


!------------------------
  subroutine calcenergy
!------------------------
!
! Driver of energy calculation
!
      use basis, only : nao, nshell
      use energy, only : enuc
      use iofile, only : iout
      use procpar, only : master
      use molecule, only : nmo, neleca
      use control, only : method
      implicit none
      integer :: nao2, nao3, nshell3
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmoa(:), ortho(:)
      real(8), allocatable :: dmtrx(:), xint(:), energymo(:)
      real(8), allocatable :: overinv(:), work(:)
!
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(nshell*(nshell+1))/2
!
! Set arrays 1
!
      call memset(nao3*4+nao2*2+nshell3+nao)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmoa(nao2),ortho(nao2),&
&              dmtrx(nao3),xint(nshell3),energymo(nao))
!
! Calculate nuclear repulsion energy
!
      call nucenergy
      if(master) then
        write(iout,'(" Nuclear repulsion energy =",f15.8," a.u.",/)') enuc
      endif
!
! Calculate overlap and 1-electron integrals
!
      call oneei(h1mtrx,smtrx,tmtrx)
!
! Set arrays 2
!
      call memset(nao2*2)
      allocate(overinv(nao2),work(nao2))
!
! Calculate canonicalization and inverse overlap matrices
!
      call fullmtrx(smtrx,work,nao)
      call mtrxcanoninv(ortho,overinv,work,nao,nmo)
!
! Calculate initial MOs and density matrix
!
      call guessmo(cmoa,overinv)
      call calcdmtrx(cmoa,dmtrx,work,nao,neleca)
!
! Unset arrays 2
!
      deallocate(overinv,work)
      call memunset(nao2*2)
!
! Start SCF
!
      if(method == "HF") then
        call calcrhf(h1mtrx,cmoa,ortho,smtrx,dmtrx,xint,energymo)
      elseif(method == "B3LYP") then
        call calcrhf(h1mtrx,cmoa,ortho,smtrx,dmtrx,xint,energymo)
        call tstamp(1)
        call calcrdft(h1mtrx,cmoa,ortho,smtrx,dmtrx,xint,energymo)
        call tstamp(1)
      elseif(method == "MP2") then
        call calcrhf(h1mtrx,cmoa,ortho,smtrx,dmtrx,xint,energymo)
        call tstamp(1)
        call calcrmp2(cmoa,energymo,xint)
        call tstamp(1)
      endif
!
! Unset arrays 1
!
      deallocate(h1mtrx,smtrx,tmtrx,cmoa,ortho, &
&                dmtrx,xint,energymo)
      call memunset(nao3*4+nao2*2+nshell3+nao)
      return
end


!--------------------------
  subroutine calcgradient
!--------------------------
!
! Driver of energy gradient calculation
!
      use basis, only : nao, nshell
      use energy, only : enuc
      use iofile, only : iout
      use procpar, only : master
      use molecule, only : nmo, neleca, natom
      implicit none
      integer :: nao2, nao3, nshell3
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmoa(:), ortho(:)
      real(8), allocatable :: dmtrx(:), xint(:), energymo(:)
      real(8), allocatable :: overinv(:), work(:)
      real(8), allocatable :: egrad(:), fulldmtrx(:), ewdmtrx(:)
!
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(nshell*(nshell+1))/2
!
! Set arrays 1
!
      call memset(nao3*4+nao2*2+nshell3+nao)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmoa(nao2),ortho(nao2),&
&              dmtrx(nao3),xint(nshell3),energymo(nao))
!
! Calculate nuclear repulsion energy
!
      call nucenergy
      if(master) then
        write(iout,'(" Nuclear repulsion energy =",f15.8," a.u.",/)') enuc
      endif
!
! Calculate overlap and 1-electron integrals
!
      call oneei(h1mtrx,smtrx,tmtrx)
!
! Set arrays 2
!
      call memset(nao2*2)
      allocate(overinv(nao2),work(nao2))
!
! Calculate canonicalization and inverse overlap matrices
!
      call fullmtrx(smtrx,work,nao)
      call mtrxcanoninv(ortho,overinv,work,nao,nmo)
!
! Calculate initial MOs and density matrix
!
      call guessmo(cmoa,overinv)
      call calcdmtrx(cmoa,dmtrx,work,nao,neleca)
!
! Unset arrays 2
!
      deallocate(overinv,work)
      call memunset(nao2*2)
!
! Start SCF
!
      call calcrhf(h1mtrx,cmoa,ortho,smtrx,dmtrx,xint,energymo)
      call tstamp(1)
!
! Set arrays 3
!
      call memset(natom*3+nao2+nao3)
      allocate(egrad(natom*3),fulldmtrx(nao2),ewdmtrx(nao3))
!
! Calculate energy gradient
!
      call calcgradrhf(cmoa,energymo,xint,egrad,fulldmtrx,ewdmtrx)
      call tstamp(1)
!
! Unset arrays 3
!
      deallocate(egrad,fulldmtrx,ewdmtrx)
      call memunset(natom*3+nao2+nao3)
!
! Unset arrays 1
!
      deallocate(h1mtrx,smtrx,tmtrx,cmoa,ortho, &
&                dmtrx,xint,energymo)
      call memunset(nao3*4+nao2*2+nshell3+nao)
      return
end


!--------------------------
  subroutine calcgeometry
!--------------------------
!
! Driver of geometry optimization calculation
!
      use basis, only : nao, nshell
      use energy, only : enuc
      use iofile, only : iout
      use procpar, only : master
      use molecule, only : nmo, neleca, natom, coord
      use opt, only : nopt, optconv
      use warn, only : nwarn
      implicit none
      integer :: nao2, nao3, nshell3, natom3, i, iopt
      real(8), allocatable :: h1mtrx(:), smtrx(:), tmtrx(:), cmoa(:), ortho(:)
      real(8), allocatable :: dmtrx(:), xint(:), energymo(:)
      real(8), allocatable :: egrad(:), egradold(:), displc(:), ehess(:), coordold(:)
      real(8), allocatable :: overinv(:), work(:)
      real(8), allocatable :: fulldmtrx(:), ewdmtrx(:)
      real(8), parameter :: zero=0.0D+00, fifth=0.2D+00, third=0.3333333333333333D+00
      real(8) :: egradmax, egradrms
!
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      nshell3=(nshell*(nshell+1))/2
      natom3= natom*3
!
! Set arrays 1
!
      call memset(nao3*4+nao2*2+nshell3+nao+natom3*6+natom3*(natom3+1)/2)
      allocate(h1mtrx(nao3),smtrx(nao3),tmtrx(nao3),cmoa(nao2),ortho(nao2), &
&              dmtrx(nao3),xint(nshell3),energymo(nao),egrad(natom3),egradold(natom3), &
&              displc(natom3*3),ehess(natom3*(natom3+1)/2),coordold(natom3))
!
! Set initial hessian
!
      ehess(:)= zero
      do i= 1,natom3
        ehess(i*(i+1)/2)= 1.5d0
      enddo
!ishimura
!     call kazuya(ehess,natom,natom3)
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
          write(iout,'(" Nuclear repulsion energy =",f15.8," a.u.",/)') enuc
        endif
!
! Calculate overlap and 1-electron integrals
!
        call oneei(h1mtrx,smtrx,tmtrx)
!
! Set arrays 2
!
        call memset(nao2*2)
        allocate(overinv(nao2),work(nao2))
!
! Calculate canonicalization and inverse overlap matrices
!
        call fullmtrx(smtrx,work,nao)
        call mtrxcanoninv(ortho,overinv,work,nao,nmo)
!
! Calculate initial MOs and density matrix
!
        call guessmo(cmoa,overinv)
        call calcdmtrx(cmoa,dmtrx,work,nao,neleca)
!
! Unset arrays 2
!
        deallocate(overinv,work)
        call memunset(nao2*2)
!
! Start SCF
!
        call calcrhf(h1mtrx,cmoa,ortho,smtrx,dmtrx,xint,energymo)
        call tstamp(1)
!
! Set arrays 3
!
        call memset(nao2+nao3)
        allocate(fulldmtrx(nao2),ewdmtrx(nao3))
!
! Calculate energy gradient
!
        call calcgradrhf(cmoa,energymo,xint,egrad,fulldmtrx,ewdmtrx)
        call tstamp(1)
!
! Unset arrays 3
!
        deallocate(fulldmtrx,ewdmtrx)
        call memunset(nao2+nao3)
!
! Calculate maximum and root mean square gradient values
!
        call calcmaxgrad(egradmax,egradrms,egrad,natom3)
        if(master) write(iout,'(" Cycle ",i3,"    Maximum gradient =",f13.8, &
&                               "  RMS gradient =",f13.8,/)')iopt,egradmax,egradrms
!
! Check convergence
!
        if((egradmax <= optconv).and.(egradrms <= optconv*third)) then
          if(master) write(iout,'(" Geometry optimization is converged.",/)')
          exit
        endif
!
! Updata hessian matrix
!
        if(iopt >= 2) call hessianbfgs(ehess,coord,coordold,egrad,egradold,displc,natom3)
!ishimura
!   call ishi(ehess,natom3)
!
! Calculate new coordinate
!
        call calcnewcoord(coord,coordold,egrad,ehess,natom3)
!ishimura
!       call ccmethod(coord,coordold,egrad,egradold,displc,natom3,iopt)
!
! Set guess MO calculation flag from Huckel to projection
! Copy energy gradient values
!
        call setnextopt(egrad,egradold,coordold,natom)
!
        if((iopt == nopt).and.master) then
          nwarn= nwarn+1
          write(iout,'("Warning! Geometry optimization is not converged.")')
        endif
      enddo
!
! End of optimization cycle 
!
!
! Unset arrays 1
!
      deallocate(h1mtrx,smtrx,tmtrx,cmoa,ortho, &
&                dmtrx,xint,energymo,egrad,egradold, &
&                displc,ehess,coordold)
      call memunset(nao3*4+nao2*2+nshell3+nao+natom3*6+natom3*(natom3+1)/2)
      return
end

subroutine ccmethod(coord,coordold,egrad,egradold,displc,natom3,iopt)
real(8):: coord(natom3),coordold(natom3),egrad(natom3),egradold(natom3),displc(natom3)
real(8) :: ddot,tmp1,tmp2
coordold=coord
if(iopt==1)then
  do i=1,natom3
     displc(i)=-egrad(i)
     coord(i)=coord(i)+displc(i)
  enddo
else
  tmp1=0.0D0
  do i=1,natom3
    tmp1=tmp1+egrad(i)*(egrad(i)-egradold(i))
  enddo
! tmp1=ddot(natom3,egrad,1,egrad,1)
  tmp2=ddot(natom3,egradold,1,egradold,1)
  tmp1=tmp1/tmp2
  do i=1,natom3
    displc(i)=-egrad(i)+tmp1*displc(i)
    coord(i)=coord(i)+displc(i)
  enddo
endif
end


subroutine ishi(ehess,ndim)
real(8) :: ehess(ndim*(ndim+1)/2),work(ndim,ndim),eigen(ndim)

ij=0
do i=1,ndim
  do j=1,i
    ij=ij+1
    work(j,i)=ehess(ij)
  enddo
enddo
      call diag('V','U',ndim,work,ndim,eigen)
    write(*,*)"EHESS"
    write(*,'(5f17.8)')ehess
    write(*,*)"Eigen"
    write(*,'(5f17.8)')eigen
end

subroutine kazuya(ehess,natom,natom3)
      use molecule, only :  coord
real(8) :: ehess(natom3*(natom3+1)/2)
real(8) :: work1(natom3,natom3),eigen(natom3),work2(natom3,natom3),work3(natom3,natom3)
real(8),parameter :: zero=0.0D0,one=1.0D0

       open(unit=20,file='helfey')
      do i=1,natom3
       read(20,*) (work1(j,i),j=1,natom3)
      enddo

!     call diag('V','U',natom3,work1,natom3,eigen)
! write(*,*)eigen
!     do i= 1,natom3
!       eigeninv= one/eigen(i)
!      if(abs(eigen(i)).lt.3.d-3)eigeninv=0.0
!       do j= 1,natom3
!         work2(j,i)= work1(j,i)*eigeninv
!       enddo
!     enddo
!     call dgemm('N','T',natom3,natom3,natom3,one,work1,natom3,work2,natom3,zero,work3,natom3)

      ij=0
      do i= 1,natom3
        do j= 1,i
          ij=ij+1
          ehess(ij)= work1(j,i)
        enddo
      enddo
      close(20)
end


!--------------------------------------------------------
  subroutine setnextopt(egrad,egradold,coordold,natom)
!--------------------------------------------------------
!
! Set parameters for next optimization step
!
      use basis, only : ex, coeff, nshell, nao, nprim, locprim, locbf, &
&                       locatom, mprim, mbf, mtype, spher
      use guess, only : ex_g, coeff_g, nshell_g, nao_g, nprim_g, locprim_g, locbf_g, &
&                       locatom_g, mprim_g, mbf_g, mtype_g, spher_g, coord_g, iguess
      implicit none
      integer,intent(in) :: natom
      integer :: iprim, ishell, iatom
      real(8),intent(in) :: egrad(natom*3), coordold(3,natom)
      real(8),intent(out) :: egradold(natom*3)
!
! Set MO projection as initial MO calculation
!
      iguess= 2
!
! Copy basis set, coordinate, sherical information
!
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
      do iatom= 1,natom
        coord_g(1,iatom)= coordold(1,iatom)
        coord_g(2,iatom)= coordold(2,iatom)
        coord_g(3,iatom)= coordold(3,iatom)
      enddo
!
      spher_g= spher
!
! Copy energy gradient
!
      egradold(:)= egrad(:)
!
      return
end







