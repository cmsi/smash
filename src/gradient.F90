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
!---------------------------------------------------------------------------
  subroutine calcgradrhf(cmo,energymo,xint,egrad,nproc1,myrank1,mpi_comm1)
!---------------------------------------------------------------------------
!
! Driver of RHF energy gradient calculation
!
      use modparallel, only : master
      use modbasis, only : nshell, nao, mtype
      use modmolecule, only : natom, neleca, numatomic
      implicit none
      integer,intent(in) :: nproc1, myrank1, mpi_comm1
      integer :: nao2, nao3, maxdim, maxgraddim, maxfunc(0:7), i, j
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: cmo(nao*nao), energymo(nao), xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: egrad(3,natom)
      real(8),allocatable :: fulldmtrx(:), ewdmtrx(:), egradtmp(:)
      character(len=3) :: table(-9:112)= &
&     (/'Bq9','Bq8','Bq7','Bq6','Bq5','Bq4','Bq3','Bq2','Bq ','X  ',&
&       'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ','Al ','Si ','P  ',&
&       'S  ','Cl ','Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',&
&       'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ',&
&       'Pd ','Ag ','Cd ','In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',&
&       'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ','Ta ','W  ','Re ',&
&       'Os ','Ir ','Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',&
&       'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ','Md ','No ','Lr ','Rf ','Db ',&
&       'Sg ','Bh ','Hs ','Mt ','Uun','Uuu','Uub'/)
      data maxfunc/1,3,6,10,15,21,28,36/
!
      if(master) then
        write(*,'(" --------------------------------------------")')
        write(*,'("   Hartree-Fock energy gradient calculation")')
        write(*,'(" --------------------------------------------")')
      endif
!
! Set arrays
!
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      call memset(nao2+nao3+natom*3)
      allocate(fulldmtrx(nao2),ewdmtrx(nao3),egradtmp(natom*3))
!
      egradtmp(:)= zero
!
! Calculate energy gradient of nuclear repulsion 
!
      call nucgradient(egradtmp,nproc1,myrank1)
!
! Calculate energy-weighted and full density matrix
!
      call calcewdmtrx(cmo,energymo,fulldmtrx,ewdmtrx,nao,neleca)
!
! Calculate derivatives of one-electron integrals
!
      call gradoneei(egradtmp,egrad,fulldmtrx,ewdmtrx,nproc1,myrank1)
!
! Calculate derivatives for two-electron integrals
!
      maxdim= maxval(mtype(1:nshell))
      maxgraddim= maxfunc(maxdim+1)
      maxdim= maxfunc(maxdim)
      call grad2eri(egradtmp,egrad,fulldmtrx,fulldmtrx,xint,one, &
&                   maxdim,maxgraddim,nproc1,myrank1,1)
!
      call para_allreducer(egradtmp(1),egrad(1,1),3*natom,mpi_comm1)
!
      if(master) then
        write(*,'(" ----------------------------------------------------")')
        write(*,'("          Gradient (Hartree/Bohr)")')
        write(*,'("  Atom            X             Y             Z")')
        write(*,'(" ----------------------------------------------------")')
        do i= 1,natom
          write(*,'(3x,a3,3x,3f14.7)')table(numatomic(i)),(egrad(j,i),j=1,3)
        enddo
        write(*,'(" ----------------------------------------------------")')
      endif
!
! Unset arrays
!
      deallocate(fulldmtrx,ewdmtrx,egradtmp)
      call memunset(nao2+nao3+natom*3)
!
      return
end


!--------------------------------------------------------------------------------------------
  subroutine calcgraduhf(cmoa,cmob,energymoa,energymob,xint,egrad,nproc1,myrank1,mpi_comm1)
!--------------------------------------------------------------------------------------------
!
! Driver of UHF energy gradient calculation
!
      use modparallel, only : master
      use modbasis, only : nshell, nao, mtype
      use modmolecule, only : natom, neleca, nelecb, numatomic
      implicit none
      integer,intent(in) :: nproc1, myrank1, mpi_comm1
      integer :: nao2, nao3, maxdim, maxgraddim, maxfunc(0:7), i, j
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: cmoa(nao*nao), cmob(nao*nao), energymoa(nao), energymob(nao)
      real(8),intent(in) :: xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: egrad(3,natom)
      real(8),allocatable :: fulldmtrx1(:), fulldmtrx2(:), ewdmtrx(:), egradtmp(:)
      character(len=3) :: table(-9:112)= &
&     (/'Bq9','Bq8','Bq7','Bq6','Bq5','Bq4','Bq3','Bq2','Bq ','X  ',&
&       'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ','Al ','Si ','P  ',&
&       'S  ','Cl ','Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',&
&       'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ',&
&       'Pd ','Ag ','Cd ','In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',&
&       'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ','Ta ','W  ','Re ',&
&       'Os ','Ir ','Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',&
&       'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ','Md ','No ','Lr ','Rf ','Db ',&
&       'Sg ','Bh ','Hs ','Mt ','Uun','Uuu','Uub'/)
      data maxfunc/1,3,6,10,15,21,28,36/
!
      if(master) then
        write(*,'(" ---------------------------------------------------------")')
        write(*,'("   Unrestricted Hartree-Fock energy gradient calculation")')
        write(*,'(" ---------------------------------------------------------")')
      endif
!
! Set arrays
!
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      call memset(nao2*2+nao3+natom*3)
      allocate(fulldmtrx1(nao2),fulldmtrx2(nao2),ewdmtrx(nao3),egradtmp(natom*3))
!
      egradtmp(:)= zero
!
! Calculate energy gradient of nuclear repulsion 
!
      call nucgradient(egradtmp,nproc1,myrank1)
!
! Calculate energy-weighted and full density matrix
!
      call calcuewdmtrx(cmoa,cmob,energymoa,energymob,fulldmtrx1,fulldmtrx2,ewdmtrx, &
&                       nao,neleca,nelecb)
!
! Calculate derivatives for one-electron integrals
!
      call gradoneei(egradtmp,egrad,fulldmtrx1,ewdmtrx,nproc1,myrank1)
!
! Calculate derivatives for two-electron integrals
!
      maxdim= maxval(mtype(1:nshell))
      maxgraddim= maxfunc(maxdim+1)
      maxdim= maxfunc(maxdim)
      call grad2eri(egradtmp,egrad,fulldmtrx1,fulldmtrx2,xint,one, &
&                   maxdim,maxgraddim,nproc1,myrank1,2)
!
      call para_allreducer(egradtmp(1),egrad(1,1),3*natom,mpi_comm1)
!
      if(master) then
        write(*,'(" ----------------------------------------------------")')
        write(*,'("          Gradient (Hartree/Bohr)")')
        write(*,'("  Atom            X             Y             Z")')
        write(*,'(" ----------------------------------------------------")')
        do i= 1,natom
          write(*,'(3x,a3,3x,3f14.7)')table(numatomic(i)),(egrad(j,i),j=1,3)
        enddo
        write(*,'(" ----------------------------------------------------")')
      endif
!
! Unset arrays
!
      deallocate(fulldmtrx1,fulldmtrx2,ewdmtrx,egradtmp)
      call memunset(nao2*2+nao3+natom*3)
!
      return
end


!----------------------------------------------------------------------------
  subroutine calcgradrdft(cmo,energymo,xint,egrad,nproc1,myrank1,mpi_comm1)
!----------------------------------------------------------------------------
!
! Driver of closed-shell DFT energy gradient calculation
!
      use modparallel, only : master
      use modbasis, only : nshell, nao, mtype
      use modmolecule, only : natom, neleca, numatomic
      use moddft, only : nrad, nleb
      use modatom, only : atomrad
      use modunit, only : tobohr
      use moddft, only : idftex, idftcor, hfexchange
      implicit none
      integer,intent(in) :: nproc1, myrank1, mpi_comm1
      integer :: nao2, nao3, maxdim, maxgraddim, maxfunc(0:7), i, j, iatom
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: cmo(nao*nao), energymo(nao), xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: egrad(3,natom)
      real(8),allocatable :: fulldmtrx(:), ewdmtrx(:), egradtmp(:,:)
      real(8),allocatable :: atomvec(:), surface(:), radpt(:), angpt(:), rad(:), ptweight(:)
      real(8),allocatable :: xyzpt(:), rsqrd(:), rr(:), uvec(:), vao(:), vmo(:)
      real(8),allocatable :: dweight(:), dpa(:), pa(:), work(:)
      character(len=3) :: table(-9:112)= &
&     (/'Bq9','Bq8','Bq7','Bq6','Bq5','Bq4','Bq3','Bq2','Bq ','X  ',&
&       'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ','Al ','Si ','P  ',&
&       'S  ','Cl ','Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',&
&       'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ',&
&       'Pd ','Ag ','Cd ','In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',&
&       'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ','Ta ','W  ','Re ',&
&       'Os ','Ir ','Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',&
&       'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ','Md ','No ','Lr ','Rf ','Db ',&
&       'Sg ','Bh ','Hs ','Mt ','Uun','Uuu','Uub'/)
      data maxfunc/1,3,6,10,15,21,28,36/
!
      if(master) then
        write(*,'(" -----------------------------------")')
        write(*,'("   DFT energy gradient calculation")')
        write(*,'(" -----------------------------------")')
      endif
!
! Set arrays
!
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      call memset(nao2+nao3+natom*3)
      allocate(fulldmtrx(nao2),ewdmtrx(nao3),egradtmp(3,natom))
      call memset(natom*natom*9+natom*13+nrad*2+nleb*4+nleb*nrad*natom+nao*10+neleca*(nao+4))
      allocate(atomvec(5*natom*natom),surface(natom*natom),radpt(2*nrad),angpt(4*nleb), &
&              rad(natom),ptweight(nleb*nrad*natom),xyzpt(3*natom),rsqrd(natom),rr(natom), &
&              uvec(3*natom),vao(10*nao),vmo(4*neleca),dweight(3*natom), &
&              dpa(3*natom*natom),pa(natom),work(neleca*nao))
!
      egradtmp(:,:)= zero
!
! Calculate energy gradient of nuclear repulsion 
!
      call nucgradient(egradtmp,nproc1,myrank1)
!
! Calculate energy-weighted and full density matrix
!
      call calcewdmtrx(cmo,energymo,fulldmtrx,ewdmtrx,nao,neleca)
!
! Calculate derivatives of one-electron integrals
!
      call gradoneei(egradtmp,egrad,fulldmtrx,ewdmtrx,nproc1,myrank1)
!
! Calculate DFT information
!
      call calcatomvec(atomvec,surface)
      call calcradpt(radpt,nrad)
      call calclebpt(angpt,nleb)
      do iatom= 1,natom
        rad(iatom)= atomrad(numatomic(iatom))*tobohr
      enddo
      call calcgridweight(ptweight,rad,radpt,angpt,atomvec,surface,xyzpt,dweight,nproc1,myrank1)
!
! Calculate derivatives of two-electron integrals
!
      maxdim= maxval(mtype(1:nshell))
      maxgraddim= maxfunc(maxdim+1)
      maxdim= maxfunc(maxdim)
      call grad2eri(egradtmp,egrad,fulldmtrx,fulldmtrx,xint,hfexchange, &
&                   maxdim,maxgraddim,nproc1,myrank1,1)
!
! Calculate derivatives of exchange-correlation terms 
!
      call gradrexcor(egradtmp,egrad,cmo,fulldmtrx,atomvec,surface,radpt,angpt,rad,ptweight, &
&                     xyzpt,rsqrd,rr,uvec,vao,vmo,dweight,dpa,pa,work,idftex,idftcor, &
&                     nproc1,myrank1)
!
      call para_allreducer(egradtmp(1,1),egrad(1,1),3*natom,mpi_comm1)
!
      if(master) then
        write(*,'(" ----------------------------------------------------")')
        write(*,'("          Gradient (Hartree/Bohr)")')
        write(*,'("  Atom            X             Y             Z")')
        write(*,'(" ----------------------------------------------------")')
        do i= 1,natom
          write(*,'(3x,a3,3x,3f14.7)')table(numatomic(i)),(egrad(j,i),j=1,3)
        enddo
        write(*,'(" ----------------------------------------------------")')
      endif
!
! Unset arrays
!
      deallocate(atomvec,surface,radpt,angpt, &
&                rad,ptweight,xyzpt,rsqrd,rr, &
&                uvec,vao,vmo,dweight, &
&                dpa,pa,work)
      call memunset(natom*natom*9+natom*13+nrad*2+nleb*4+nleb*nrad*natom+nao*10+neleca*(nao+4))
      deallocate(fulldmtrx,ewdmtrx,egradtmp)
      call memunset(nao2+nao3+natom*3)
!
      return
end


!---------------------------------------------------------------------------------------------
  subroutine calcgradudft(cmoa,cmob,energymoa,energymob,xint,egrad,nproc1,myrank1,mpi_comm1)
!---------------------------------------------------------------------------------------------
!
! Driver of open-shell DFT energy gradient calculation
!
      use modparallel, only : master
      use modbasis, only : nshell, nao, mtype
      use modmolecule, only : natom, neleca, nelecb, numatomic
      use moddft, only : nrad, nleb
      use modatom, only : atomrad
      use modunit, only : tobohr
      use moddft, only : idftex, idftcor, hfexchange
      implicit none
      integer :: nao2, nao3, maxdim, maxgraddim, maxfunc(0:7), i, j, iatom
      integer,intent(in) :: nproc1, myrank1, mpi_comm1
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: cmoa(nao*nao), cmob(nao*nao), energymoa(nao), energymob(nao)
      real(8),intent(in) :: xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: egrad(3,natom)
      real(8),allocatable :: fulldmtrx1(:), fulldmtrx2(:), ewdmtrx(:), egradtmp(:,:)
      real(8),allocatable :: atomvec(:), surface(:), radpt(:), angpt(:), rad(:), ptweight(:)
      real(8),allocatable :: xyzpt(:), rsqrd(:), rr(:), uvec(:), vao(:), vmoa(:)
      real(8),allocatable :: vmob(:), dweight(:), dpa(:), pa(:), work(:)
      character(len=3) :: table(-9:112)= &
&     (/'Bq9','Bq8','Bq7','Bq6','Bq5','Bq4','Bq3','Bq2','Bq ','X  ',&
&       'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ','Al ','Si ','P  ',&
&       'S  ','Cl ','Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',&
&       'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ',&
&       'Pd ','Ag ','Cd ','In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',&
&       'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ','Ta ','W  ','Re ',&
&       'Os ','Ir ','Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',&
&       'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ','Md ','No ','Lr ','Rf ','Db ',&
&       'Sg ','Bh ','Hs ','Mt ','Uun','Uuu','Uub'/)
      data maxfunc/1,3,6,10,15,21,28,36/
!
      if(master) then
        write(*,'(" ------------------------------------------------")')
        write(*,'("   Unrestricted DFT energy gradient calculation")')
        write(*,'(" ------------------------------------------------")')
      endif
!
! Set arrays
!
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      call memset(nao2*2+nao3+natom*3)
      allocate(fulldmtrx1(nao2),fulldmtrx2(nao2),ewdmtrx(nao3),egradtmp(3,natom))
      call memset(natom*natom*9+natom*13+nrad*2+nleb*4+nleb*nrad*natom+nao*10 &
&                +neleca*4+nelecb*4+(neleca+nelecb)*nao)
      allocate(atomvec(5*natom*natom),surface(natom*natom),radpt(2*nrad),angpt(4*nleb), &
&              rad(natom),ptweight(nleb*nrad*natom),xyzpt(3*natom),rsqrd(natom),rr(natom), &
&              uvec(3*natom),vao(10*nao),vmoa(4*neleca),vmob(4*nelecb), &
&              dweight(3*natom),dpa(3*natom*natom),pa(natom),work((neleca+nelecb)*nao))
!
      egradtmp(:,:)= zero
!
! Calculate energy gradient of nuclear repulsion 
!
      call nucgradient(egradtmp,nproc1,myrank1)
!
! Calculate energy-weighted and full density matrix
!
      call calcuewdmtrx(cmoa,cmob,energymoa,energymob,fulldmtrx1,fulldmtrx2,ewdmtrx, &
&                       nao,neleca,nelecb)
!
! Calculate derivatives of one-electron integrals
!
      call gradoneei(egradtmp,egrad,fulldmtrx1,ewdmtrx,nproc1,myrank1)
!
! Calculate DFT information
!
      call calcatomvec(atomvec,surface)
      call calcradpt(radpt,nrad)
      call calclebpt(angpt,nleb)
      do iatom= 1,natom
        rad(iatom)= atomrad(numatomic(iatom))*tobohr
      enddo
      call calcgridweight(ptweight,rad,radpt,angpt,atomvec,surface,xyzpt,dweight,nproc1,myrank1)
!
! Calculate derivatives of two-electron integrals
!
      maxdim= maxval(mtype(1:nshell))
      maxgraddim= maxfunc(maxdim+1)
      maxdim= maxfunc(maxdim)
      call grad2eri(egradtmp,egrad,fulldmtrx1,fulldmtrx2,xint,hfexchange, &
&                   maxdim,maxgraddim,nproc1,myrank1,2)
!
! Calculate derivatives of exchange-correlation terms 
!
      call graduexcor(egradtmp,egrad,cmoa,cmob,fulldmtrx1,fulldmtrx2,atomvec,surface,radpt, &
&                     angpt,rad,ptweight,xyzpt,rsqrd,rr,uvec,vao,vmoa,vmob,dweight, &
&                     dpa,pa,work,work(neleca*nao+1),idftex,idftcor,nproc1,myrank1)
!
      call para_allreducer(egradtmp(1,1),egrad(1,1),3*natom,mpi_comm1)
!
      if(master) then
        write(*,'(" ----------------------------------------------------")')
        write(*,'("          Gradient (Hartree/Bohr)")')
        write(*,'("  Atom            X             Y             Z")')
        write(*,'(" ----------------------------------------------------")')
        do i= 1,natom
          write(*,'(3x,a3,3x,3f14.7)')table(numatomic(i)),(egrad(j,i),j=1,3)
        enddo
        write(*,'(" ----------------------------------------------------")')
      endif
!
! Unset arrays
!
      deallocate(atomvec,surface,radpt,angpt, &
&                rad,ptweight,xyzpt,rsqrd,rr, &
&                uvec,vao,vmoa,vmob, &
&                dweight,dpa,pa,work)
      call memunset(natom*natom*9+natom*13+nrad*2+nleb*4+nleb*nrad*natom+nao*10 &
&                  +neleca*4+nelecb*4+(neleca+nelecb)*nao)
      deallocate(fulldmtrx1,fulldmtrx2,ewdmtrx,egradtmp)
      call memunset(nao2*2+nao3+natom*3)
!
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


