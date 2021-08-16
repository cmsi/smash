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
!-----------------------------------------------------------------------
  subroutine calcgradrhf(cmo,energymo,xint,egrad,nproc1,myrank1, &
&                        mpi_comm1,datajob,datamol,databasis,datacomp)
!-----------------------------------------------------------------------
!
! Driver of RHF energy gradient calculation
!
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(in) :: datajob
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer,intent(in) :: nproc1, myrank1, mpi_comm1
      integer :: nao, natom, nao2, nao3, maxdim, maxgraddim, maxfunc(0:7), i, j
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: cmo(databasis%nao*databasis%nao), energymo(databasis%nao)
      real(8),intent(in) :: xint(databasis%nshell*(databasis%nshell+1)/2)
      real(8),intent(out) :: egrad(3,datamol%natom)
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
&       'Sg ','Bh ','Hs ','Mt ','Ds ','Rg ','Cn '/)
      data maxfunc/1,3,6,10,15,21,28,36/
!
      if(datacomp%master) then
        write(datacomp%iout,'(" ============================================")')
        write(datacomp%iout,'("   Hartree-Fock energy gradient calculation")')
        write(datacomp%iout,'(" ============================================")')
      endif
!
! Set arrays
!
      nao= databasis%nao
      natom= datamol%natom
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      call memset(nao2+nao3+natom*3,datacomp)
      allocate(fulldmtrx(nao2),ewdmtrx(nao3),egradtmp(natom*3))
!
      egradtmp(:)= zero
!
! Calculate energy gradient of nuclear repulsion 
!
      call nucgradient(egradtmp,nproc1,myrank1,datamol)
!
! Calculate energy-weighted and full density matrix
!
      call calcewdmtrx(cmo,energymo,fulldmtrx,ewdmtrx,nao,datamol%neleca)
!
! Calculate derivatives of one-electron integrals
!
      call gradoneei(egradtmp,egrad,fulldmtrx,ewdmtrx,nproc1,myrank1, &
&                    datajob,datamol,databasis,datacomp)
!
! Calculate derivatives for two-electron integrals
!
      maxdim= maxval(databasis%mtype(1:databasis%nshell))
      maxgraddim= maxfunc(maxdim+1)
      maxdim= maxfunc(maxdim)
      call grad2eri(egradtmp,egrad,fulldmtrx,fulldmtrx,xint,one, &
&                   maxdim,maxgraddim,nproc1,myrank1,1,datajob,datamol,databasis,datacomp)
!
      call para_allreducer(egradtmp(1),egrad(1,1),3*natom,mpi_comm1)
!
      if(datacomp%master) then
        write(datacomp%iout,'(" ----------------------------------------------------")')
        write(datacomp%iout,'("          Gradient (Hartree/Bohr)")')
        write(datacomp%iout,'("  Atom            X             Y             Z")')
        write(datacomp%iout,'(" ----------------------------------------------------")')
        do i= 1,natom
          write(datacomp%iout,'(3x,a3,3x,3f14.7)') table(datamol%numatomic(i)),(egrad(j,i),j=1,3)
        enddo
        write(datacomp%iout,'(" ----------------------------------------------------")')
      endif
!
! Unset arrays
!
      deallocate(fulldmtrx,ewdmtrx,egradtmp)
      call memunset(nao2+nao3+natom*3,datacomp)
!
      return
end


!-------------------------------------------------------------------------------
  subroutine calcgraduhf(cmoa,cmob,energymoa,energymob,xint,egrad,nproc1, &
&                        myrank1,mpi_comm1,datajob,datamol,databasis,datacomp)
!-------------------------------------------------------------------------------
!
! Driver of UHF energy gradient calculation
!
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(in) :: datajob
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer,intent(in) :: nproc1, myrank1, mpi_comm1
      integer :: nao, natom, nao2, nao3, maxdim, maxgraddim, maxfunc(0:7), i, j
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: cmoa(databasis%nao*databasis%nao), cmob(databasis%nao*databasis%nao)
      real(8),intent(in) :: energymoa(databasis%nao), energymob(databasis%nao)
      real(8),intent(in) :: xint(databasis%nshell*(databasis%nshell+1)/2)
      real(8),intent(out) :: egrad(3,datamol%natom)
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
&       'Sg ','Bh ','Hs ','Mt ','Ds ','Rg ','Cn '/)
      data maxfunc/1,3,6,10,15,21,28,36/
!
      if(datacomp%master) then
        write(datacomp%iout,'(" =========================================================")')
        write(datacomp%iout,'("   Unrestricted Hartree-Fock energy gradient calculation")')
        write(datacomp%iout,'(" =========================================================")')
      endif
!
! Set arrays
!
      nao= databasis%nao
      natom= datamol%natom
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      call memset(nao2*2+nao3+natom*3,datacomp)
      allocate(fulldmtrx1(nao2),fulldmtrx2(nao2),ewdmtrx(nao3),egradtmp(natom*3))
!
      egradtmp(:)= zero
!
! Calculate energy gradient of nuclear repulsion 
!
      call nucgradient(egradtmp,nproc1,myrank1,datamol)
!
! Calculate energy-weighted and full density matrix
!
      call calcuewdmtrx(cmoa,cmob,energymoa,energymob,fulldmtrx1,fulldmtrx2,ewdmtrx, &
&                       nao,datamol%neleca,datamol%nelecb)
!
! Calculate derivatives for one-electron integrals
!
      call gradoneei(egradtmp,egrad,fulldmtrx1,ewdmtrx,nproc1,myrank1, &
&                    datajob,datamol,databasis,datacomp)
!
! Calculate derivatives for two-electron integrals
!
      maxdim= maxval(databasis%mtype(1:databasis%nshell))
      maxgraddim= maxfunc(maxdim+1)
      maxdim= maxfunc(maxdim)
      call grad2eri(egradtmp,egrad,fulldmtrx1,fulldmtrx2,xint,one, &
&                   maxdim,maxgraddim,nproc1,myrank1,2,datajob,datamol,databasis,datacomp)
!
      call para_allreducer(egradtmp(1),egrad(1,1),3*natom,mpi_comm1)
!
      if(datacomp%master) then
        write(datacomp%iout,'(" ----------------------------------------------------")')
        write(datacomp%iout,'("          Gradient (Hartree/Bohr)")')
        write(datacomp%iout,'("  Atom            X             Y             Z")')
        write(datacomp%iout,'(" ----------------------------------------------------")')
        do i= 1,natom
          write(datacomp%iout,'(3x,a3,3x,3f14.7)') table(datamol%numatomic(i)),(egrad(j,i),j=1,3)
        enddo
        write(datacomp%iout,'(" ----------------------------------------------------")')
      endif
!
! Unset arrays
!
      deallocate(fulldmtrx1,fulldmtrx2,ewdmtrx,egradtmp)
      call memunset(nao2*2+nao3+natom*3,datacomp)
!
      return
end


!------------------------------------------------------------------------
  subroutine calcgradrdft(cmo,energymo,xint,egrad,nproc1,myrank1, &
&                         mpi_comm1,datajob,datamol,databasis,datacomp)
!------------------------------------------------------------------------
!
! Driver of closed-shell DFT energy gradient calculation
!
      use modparam, only : tobohr
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(in) :: datajob
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer,intent(in) :: nproc1, myrank1, mpi_comm1
      integer :: nrad, nleb, nao, ndftatom, nao2, nao3, maxdim, maxgraddim, maxfunc(0:7), i, j, iatom
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: cmo(databasis%nao*databasis%nao), energymo(databasis%nao)
      real(8),intent(in) :: xint(databasis%nshell*(databasis%nshell+1)/2)
      real(8),intent(out) :: egrad(3,datamol%natom)
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
&       'Sg ','Bh ','Hs ','Mt ','Ds ','Rg ','Cn '/)
      data maxfunc/1,3,6,10,15,21,28,36/
!
      if(datacomp%master) then
        write(datacomp%iout,'(" ===================================")')
        write(datacomp%iout,'("   DFT energy gradient calculation")')
        write(datacomp%iout,'(" ===================================")')
      endif
!
! Set arrays
!
      nrad= datajob%nrad
      nleb= datajob%nleb
      nao= databasis%nao
      ndftatom= datamol%natom-datamol%ndummyatom
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      call memset(nao2+nao3+datamol%natom*3,datacomp)
      allocate(fulldmtrx(nao2),ewdmtrx(nao3),egradtmp(3,datamol%natom))
      call memset(ndftatom*ndftatom*9+ndftatom*13+nrad*2+nleb*4+nleb*nrad*ndftatom+nao*10 &
&                +datamol%neleca*(nao+4),datacomp)
      allocate(atomvec(5*ndftatom*ndftatom),surface(ndftatom*ndftatom),radpt(2*nrad), &
&              angpt(4*nleb),rad(ndftatom),ptweight(nleb*nrad*ndftatom),xyzpt(3*ndftatom), &
&              rsqrd(ndftatom),rr(ndftatom),uvec(3*ndftatom),vao(10*nao),vmo(4*datamol%neleca), &
&              dweight(3*ndftatom),dpa(3*ndftatom*ndftatom),pa(ndftatom),work(datamol%neleca*nao))
!
      egradtmp(:,:)= zero
!
! Calculate energy gradient of nuclear repulsion 
!
      call nucgradient(egradtmp,nproc1,myrank1,datamol)
!
! Calculate energy-weighted and full density matrix
!
      call calcewdmtrx(cmo,energymo,fulldmtrx,ewdmtrx,nao,datamol%neleca)
!
! Calculate derivatives of one-electron integrals
!
      call gradoneei(egradtmp,egrad,fulldmtrx,ewdmtrx,nproc1,myrank1, &
&                    datajob,datamol,databasis,datacomp)
!
! Calculate DFT information
!
      call calcatomvec(atomvec,surface,ndftatom,datamol)
      call calcradpt(radpt,nrad)
      call calclebpt(angpt,nleb,datacomp)
      do iatom= 1,ndftatom
        rad(iatom)= datamol%atomrad(datamol%numatomic(iatom))*tobohr
      enddo
      call calcgridweight(ptweight,rad,radpt,angpt,atomvec,surface,xyzpt,dweight,nrad,nleb, &
&                         ndftatom,nproc1,myrank1)
!
! Calculate derivatives of two-electron integrals
!
      maxdim= maxval(databasis%mtype(1:databasis%nshell))
      maxgraddim= maxfunc(maxdim+1)
      maxdim= maxfunc(maxdim)
      call grad2eri(egradtmp,egrad,fulldmtrx,fulldmtrx,xint,datajob%hfexchange, &
&                   maxdim,maxgraddim,nproc1,myrank1,1,datajob,datamol,databasis,datacomp)
!
! Calculate derivatives of exchange-correlation terms 
!
      call gradrexcor(egradtmp,egrad,cmo,fulldmtrx,atomvec,surface,radpt,angpt,rad,ptweight, &
&                     xyzpt,rsqrd,rr,uvec,vao,vmo,dweight,dpa,pa,work,ndftatom,datajob%idftex, &
&                     datajob%idftcor,nproc1,myrank1,datajob,datamol,databasis,datacomp)
!
      call para_allreducer(egradtmp(1,1),egrad(1,1),3*datamol%natom,mpi_comm1)
!
      if(datacomp%master) then
        write(datacomp%iout,'(" ----------------------------------------------------")')
        write(datacomp%iout,'("          Gradient (Hartree/Bohr)")')
        write(datacomp%iout,'("  Atom            X             Y             Z")')
        write(datacomp%iout,'(" ----------------------------------------------------")')
        do i= 1,datamol%natom
          write(datacomp%iout,'(3x,a3,3x,3f14.7)') table(datamol%numatomic(i)),(egrad(j,i),j=1,3)
        enddo
        write(datacomp%iout,'(" ----------------------------------------------------")')
      endif
!
! Unset arrays
!
      deallocate(atomvec,surface,radpt, &
&                angpt,rad,ptweight,xyzpt, &
&                rsqrd,rr,uvec,vao,vmo, &
&                dweight,dpa,pa,work)
      call memunset(ndftatom*ndftatom*9+ndftatom*13+nrad*2+nleb*4+nleb*nrad*ndftatom+nao*10 &
&                  +datamol%neleca*(nao+4),datacomp)
      deallocate(fulldmtrx,ewdmtrx,egradtmp)
      call memunset(nao2+nao3+datamol%natom*3,datacomp)
!
      return
end


!--------------------------------------------------------------------------------
  subroutine calcgradudft(cmoa,cmob,energymoa,energymob,xint,egrad,nproc1, &
&                         myrank1,mpi_comm1,datajob,datamol,databasis,datacomp)
!--------------------------------------------------------------------------------
!
! Driver of open-shell DFT energy gradient calculation
!
      use modparam, only : tobohr
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(in) :: datajob
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer,intent(in) :: nproc1, myrank1, mpi_comm1
      integer :: nrad, nleb, nao, ndftatom, nao2, nao3, maxdim, maxgraddim, maxfunc(0:7)
      integer :: i, j, iatom
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: cmoa(databasis%nao*databasis%nao), cmob(databasis%nao*databasis%nao)
      real(8),intent(in) :: energymoa(databasis%nao), energymob(databasis%nao)
      real(8),intent(in) :: xint(databasis%nshell*(databasis%nshell+1)/2)
      real(8),intent(out) :: egrad(3,datamol%natom)
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
&       'Sg ','Bh ','Hs ','Mt ','Ds ','Rg ','Cn '/)
      data maxfunc/1,3,6,10,15,21,28,36/
!
      if(datacomp%master) then
        write(datacomp%iout,'(" ================================================")')
        write(datacomp%iout,'("   Unrestricted DFT energy gradient calculation")')
        write(datacomp%iout,'(" ================================================")')
      endif
!
! Set arrays
!
      nrad= datajob%nrad
      nleb= datajob%nleb
      nao= databasis%nao
      ndftatom= datamol%natom-datamol%ndummyatom
      nao2= nao*nao
      nao3=(nao*(nao+1))/2
      call memset(nao2*2+nao3+datamol%natom*3,datacomp)
      allocate(fulldmtrx1(nao2),fulldmtrx2(nao2),ewdmtrx(nao3),egradtmp(3,datamol%natom))
      call memset(ndftatom*ndftatom*9+ndftatom*13+nrad*2+nleb*4+nleb*nrad*ndftatom+nao*10 &
&                +datamol%neleca*4+datamol%nelecb*4+(datamol%neleca+datamol%nelecb)*nao,datacomp)
      allocate(atomvec(5*ndftatom*ndftatom),surface(ndftatom*ndftatom),radpt(2*nrad), &
&              angpt(4*nleb),rad(ndftatom),ptweight(nleb*nrad*ndftatom),xyzpt(3*ndftatom), &
&              rsqrd(ndftatom),rr(ndftatom),uvec(3*ndftatom),vao(10*nao),vmoa(4*datamol%neleca), &
&              vmob(4*datamol%nelecb),dweight(3*ndftatom),dpa(3*ndftatom*ndftatom),pa(ndftatom), &
&              work((datamol%neleca+datamol%nelecb)*nao))
!
      egradtmp(:,:)= zero
!
! Calculate energy gradient of nuclear repulsion 
!
      call nucgradient(egradtmp,nproc1,myrank1,datamol)
!
! Calculate energy-weighted and full density matrix
!
      call calcuewdmtrx(cmoa,cmob,energymoa,energymob,fulldmtrx1,fulldmtrx2,ewdmtrx, &
&                       nao,datamol%neleca,datamol%nelecb)
!
! Calculate derivatives of one-electron integrals
!
      call gradoneei(egradtmp,egrad,fulldmtrx1,ewdmtrx,nproc1,myrank1, &
&                    datajob,datamol,databasis,datacomp)
!
! Calculate DFT information
!
      call calcatomvec(atomvec,surface,ndftatom,datamol)
      call calcradpt(radpt,nrad)
      call calclebpt(angpt,nleb,datacomp)
      do iatom= 1,ndftatom
        rad(iatom)= datamol%atomrad(datamol%numatomic(iatom))*tobohr
      enddo
      call calcgridweight(ptweight,rad,radpt,angpt,atomvec,surface,xyzpt,dweight,nrad,nleb, &
&                         ndftatom,nproc1,myrank1)
!
! Calculate derivatives of two-electron integrals
!
      maxdim= maxval(databasis%mtype(1:databasis%nshell))
      maxgraddim= maxfunc(maxdim+1)
      maxdim= maxfunc(maxdim)
      call grad2eri(egradtmp,egrad,fulldmtrx1,fulldmtrx2,xint,datajob%hfexchange, &
&                   maxdim,maxgraddim,nproc1,myrank1,2,datajob,datamol,databasis,datacomp)
!
! Calculate derivatives of exchange-correlation terms 
!
      call graduexcor(egradtmp,egrad,cmoa,cmob,fulldmtrx1,fulldmtrx2,atomvec,surface,radpt, &
&                     angpt,rad,ptweight,xyzpt,rsqrd,rr,uvec,vao,vmoa,vmob,dweight,dpa,pa, &
&                     work,work(datamol%neleca*nao+1),ndftatom,datajob%idftex,datajob%idftcor, &
&                     nproc1,myrank1,datajob,datamol,databasis,datacomp)
!
      call para_allreducer(egradtmp(1,1),egrad(1,1),3*datamol%natom,mpi_comm1)
!
      if(datacomp%master) then
        write(datacomp%iout,'(" ----------------------------------------------------")')
        write(datacomp%iout,'("          Gradient (Hartree/Bohr)")')
        write(datacomp%iout,'("  Atom            X             Y             Z")')
        write(datacomp%iout,'(" ----------------------------------------------------")')
        do i= 1,datamol%natom
          write(datacomp%iout,'(3x,a3,3x,3f14.7)') table(datamol%numatomic(i)),(egrad(j,i),j=1,3)
        enddo
        write(datacomp%iout,'(" ----------------------------------------------------")')
      endif
!
! Unset arrays
!
      deallocate(atomvec,surface,radpt, &
&                angpt,rad,ptweight,xyzpt, &
&                rsqrd,rr,uvec,vao,vmoa, &
&                vmob,dweight,dpa,pa, &
&                work)
      call memunset(ndftatom*ndftatom*9+ndftatom*13+nrad*2+nleb*4+nleb*nrad*ndftatom+nao*10 &
&                  +datamol%neleca*4+datamol%nelecb*4+(datamol%neleca+datamol%nelecb)*nao,datacomp)
      deallocate(fulldmtrx1,fulldmtrx2,ewdmtrx,egradtmp)
      call memunset(nao2*2+nao3+datamol%natom*3,datacomp)
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


