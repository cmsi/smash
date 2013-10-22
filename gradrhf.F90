!---------------------------------------------------------------------
  subroutine calcgradrhf(cmoa,energymo,xint,egrad,fulldmtrx,ewdmtrx)
!---------------------------------------------------------------------
!
! Driver of RHF energy gradient calculation
!
      use procpar
      use basis, only : nshell, nao, mtype
      use molecule, only : natom, neleca, numatomic
      use iofile, only : iout
      implicit none
      integer :: maxdim, maxfunc(0:7), i, j
      real(8),intent(in) :: cmoa(nao*nao), energymo(nao), xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: egrad(3,natom), fulldmtrx(nao*nao), ewdmtrx(nao*(nao+1)/2)
      real(8) :: dummy
      character(3) :: table(112)= &
&     (/'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ','Al ','Si ','P  ',&
&       'S  ','Cl ','Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',&
&       'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ',&
&       'Pd ','Ag ','Cd ','In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',&
&       'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ','Ta ','W  ','Re ',&
&       'Os ','Ir ','Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',&
&       'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ','Md ','No ','Lr ','Rf ','Db ',&
&       'Sg ','Bh ','Hs ','Mt ','Uun','Uuu','Uub'/)
      data maxfunc/1,3,6,10,15,21,28,36/
!
      call zeroclr(egrad,natom*3)
!
! Calculate energy gradient of nuclear repulsion 
!
      call nucgradient(egrad)

!ishimura
!     if(master) then
!       write(iout,'(" ------------------------------------------------")')
!       write(iout,'("          Force   nucgrad             ")')
!       write(iout,'("  Atom           X            Y            Z   ")')
!       write(iout,'(" ------------------------------------------------")')
!       do i= 1,natom
!         write(iout,'(3x,a3,3x,3f15.8)')table(numatomic(i)),(egrad(j,i),j=1,3)
!       enddo
!       write(iout,'(" ------------------------------------------------",/)')
!     endif
!
! Calculate energy-weighted and full density matrix
!
      call calcewdmtrx(cmoa,energymo,fulldmtrx,ewdmtrx,nao,neleca)
!
! Calculate derivatives for one-electron integrals
!
      call gradoneei(egrad,fulldmtrx,ewdmtrx)

!ishimura
!     if(master) then
!       write(iout,'(" ------------------------------------------------")')
!       write(iout,'("          Force    onegrad            ")')
!       write(iout,'("  Atom           X            Y            Z   ")')
!       write(iout,'(" ------------------------------------------------")')
!       do i= 1,natom
!         write(iout,'(3x,a3,3x,3f15.8)')table(numatomic(i)),(egrad(j,i),j=1,3)
!       enddo
!       write(iout,'(" ------------------------------------------------",/)')
!     endif

!
! Calculate derivatives for two-electron integrals
!
      maxdim=maxfunc(maxval(mtype(1:nshell))+1)
      call grad2eri(egrad,fulldmtrx,xint,maxdim)
!
      call para_allreduce(egrad,dummy,3*natom,'D',MPI_SUM,MPI_COMM_WORLD,1)
!
      if(master) then
        write(iout,'(" ------------------------------------------------")')
        write(iout,'("          Force                       ")')
        write(iout,'("  Atom           X            Y            Z   ")')
        write(iout,'(" ------------------------------------------------")')
        do i= 1,natom
          write(iout,'(3x,a3,3x,3f15.8)')table(numatomic(i)),(egrad(j,i),j=1,3)
        enddo
        write(iout,'(" ------------------------------------------------",/)')
      endif
      return
end


!---------------------------------------------------------
  subroutine calcmaxgrad(egradmax,egradrms,egrad,natom3)
!---------------------------------------------------------
!
! Calculate maximum and root mean square gradient values
!
      implicit none
      integer,intent(in) :: natom3
      integer :: ipos, idamax
      real(8),intent(in) :: egrad(natom3)
      real(8),intent(out) :: egradmax, egradrms
      real(8) :: ddot
!
      ipos= idamax(natom3,egrad,1)
      egradmax= abs(egrad(ipos))
      egradrms= ddot(natom3,egrad,1,egrad,1)
      egradrms= sqrt(egradrms/natom3)
!
      return
end


!-------------------------------------------------------------
  subroutine calcnewcoord(coord,coordold,egrad,ehess,natom3)
!-------------------------------------------------------------
!
! Calculate new coordinate using gradient and hessian
!
      use procpar
      implicit none
      integer,intent(in) :: natom3
      integer :: i, j, ii
!     real(8),intent(in) :: egrad(natom3), ehess(natom3*(natom3+1)/2)
      real(8),intent(inout) :: egrad(natom3), ehess(natom3*(natom3+1)/2)
      real(8),intent(inout) :: coord(natom3), coordold(natom3)
!
      real(8) :: work(natom3,natom3),eigen(natom3),work2(natom3,natom3),work3(natom3,natom3)
 
      do i=1,natom3
        ii= i*(i-1)/2
        do j=1,i
          work(j,i)=ehess(ii+j)
        enddo
      enddo
      call diag('V','U',natom3,work,natom3,eigen)

!ishimura
!     if(master) then
!       write(*,*)"eigen"
!       write(*,'(5f15.9)')eigen 
!     endif
!     do i=1,6
!       eigen(i)=0.0D0
!     enddo
!     do i=7,natom3
!       eigen(i)=1.0D0/eigen(i)
!     enddo
      do i=1,natom3
        eigen(i)=1.0D0/eigen(i)
      enddo
      do i=1,natom3
        do j=1,natom3
          work2(j,i)=work(j,i)*eigen(i)
        enddo
      enddo
      call dgemm('N','T',natom3,natom3,natom3,1.0D0,work,natom3,work2,natom3,0.0D0,work3,natom3)
!
! Copy old coordinate
!
      coordold(:)= coord(:)
!
      do i=1,natom3
        do j=1,natom3
          coord(i)=coord(i)-work3(i,j)*egrad(j)
        enddo
      enddo
!
! Print delta xyz
      if(master) then
        write(*,'("Delta x")')
        write(*,'(3f15.9)')(coord(i)-coordold(i),i=1,natom3)
      endif
        
!
! Copy old coordinate
!
!     coordold(:)= coord(:)
!
! Update coordinate
!
!     do i= 1,natom3
!       ii= i*(i-1)/2
!       do j= 1,i-1
!         coord(i)= coord(i)-egrad(j)*ehess(ii+j)
!         coord(j)= coord(j)-egrad(i)*ehess(ii+j)
!       enddo
!       coord(i)= coord(i)-egrad(i)*ehess(ii+i)
!     enddo
!
      return
end

