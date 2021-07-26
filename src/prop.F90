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
!---------------------------------------------------------------------
  subroutine calcrmulliken(dmtrx,overlap,datamol,databasis,datacomp)
!---------------------------------------------------------------------
!
! Execute Mulliken population analysis for closed-shell
!
      use modtype, only : typemol, typebasis, typecomp
      implicit none
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(in) :: datacomp
      integer :: ii, jj, ij, ish, iatom, locbfi
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: dmtrx(databasis%nao*(databasis%nao+1)/2)
      real(8),intent(in) :: overlap(databasis%nao*(databasis%nao+1)/2)
      real(8) :: grossorb(databasis%nao), grossatom(datamol%natom), totalgross
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
!
      totalgross= zero
      grossorb(:)= zero
      grossatom(:)= zero
!
! Calculate Gross orbital population
!
!$OMP parallel do schedule(static,1) private(ij) reduction(+:grossorb)
      do ii= 1,databasis%nao
        ij= ii*(ii-1)/2
        do jj= 1,ii-1
          ij= ij+1
          grossorb(ii)= grossorb(ii)+dmtrx(ij)*overlap(ij)
          grossorb(jj)= grossorb(jj)+dmtrx(ij)*overlap(ij)
        enddo
        ij= ij+1
        grossorb(ii)= grossorb(ii)+dmtrx(ij)*overlap(ij)
      enddo
!$OMP end parallel do
!
! Calculate Gross atom population
!
      do ish= 1,databasis%nshell
        iatom= databasis%locatom(ish)
        locbfi= databasis%locbf(ish)
        do ii= 1,databasis%mbf(ish)
          grossatom(iatom)= grossatom(iatom)+grossorb(locbfi+ii)
        enddo
      enddo
      do iatom= 1,datamol%natom
        totalgross= totalgross+datamol%znuc(iatom)-grossatom(iatom)
      enddo
!
      if(datacomp%master) then
        write(datacomp%iout,'(" -------------------------------------")')
        write(datacomp%iout,'("      Mulliken Population Analysis")')
        write(datacomp%iout,'("     Atom       Charge     Population")')
        write(datacomp%iout,'(" -------------------------------------")')
        do iatom= 1,datamol%natom
          write(datacomp%iout,'(1x,i4,2x,a3,2f13.6)')iatom,table(datamol%numatomic(iatom)), &
&                                        datamol%znuc(iatom)-grossatom(iatom),grossatom(iatom)
        enddo
        write(datacomp%iout,'(" -------------------------------------")')
        write(datacomp%iout,'("     Total",f13.6)')totalgross
        write(datacomp%iout,'(" -------------------------------------")')
        write(datacomp%iout,*)
      endif
!
      return
end


!-----------------------------------------------------------------------------
  subroutine calcumulliken(dmtrxa,dmtrxb,overlap,datamol,databasis,datacomp)
!-----------------------------------------------------------------------------
!
! Execute Mulliken population analysis for open-shell
!
      use modtype, only : typemol, typebasis, typecomp
      implicit none
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(in) :: datacomp
      integer :: ii, jj, ij, ish, iatom, locbfi
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: dmtrxa(databasis%nao*(databasis%nao+1)/2)
      real(8),intent(in) :: dmtrxb(databasis%nao*(databasis%nao+1)/2)
      real(8),intent(in) :: overlap(databasis%nao*(databasis%nao+1)/2)
      real(8) :: grossorb(databasis%nao), grossatom(datamol%natom), totalgross
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
!
      totalgross= zero
      grossorb(:)= zero
      grossatom(:)= zero
!
! Calculate Gross orbital population
!
!$OMP parallel do schedule(static,1) private(ij) reduction(+:grossorb)
      do ii= 1,databasis%nao
        ij= ii*(ii-1)/2
        do jj= 1,ii-1
          ij= ij+1
          grossorb(ii)= grossorb(ii)+(dmtrxa(ij)+dmtrxb(ij))*overlap(ij)
          grossorb(jj)= grossorb(jj)+(dmtrxa(ij)+dmtrxb(ij))*overlap(ij)
        enddo
        ij= ij+1
        grossorb(ii)= grossorb(ii)+(dmtrxa(ij)+dmtrxb(ij))*overlap(ij)
      enddo
!$OMP end parallel do
!
! Calculate Gross atom population
!
      do ish= 1,databasis%nshell
        iatom= databasis%locatom(ish)
        locbfi= databasis%locbf(ish)
        do ii= 1,databasis%mbf(ish)
          grossatom(iatom)= grossatom(iatom)+grossorb(locbfi+ii)
        enddo
      enddo
      do iatom= 1,datamol%natom
        totalgross= totalgross+datamol%znuc(iatom)-grossatom(iatom)
      enddo
!
      if(datacomp%master) then
        write(datacomp%iout,'(" -------------------------------------")')
        write(datacomp%iout,'("      Mulliken Population Analysis")')
        write(datacomp%iout,'("     Atom       Charge     Population")')
        write(datacomp%iout,'(" -------------------------------------")')
        do iatom= 1,datamol%natom
          write(datacomp%iout,'(1x,i4,2x,a3,2f13.6)')iatom,table(datamol%numatomic(iatom)), &
&                                        datamol%znuc(iatom)-grossatom(iatom),grossatom(iatom)
        enddo
        write(datacomp%iout,'(" -------------------------------------")')
        write(datacomp%iout,'("     Total",f13.6)')totalgross
        write(datacomp%iout,'(" -------------------------------------")')
        write(datacomp%iout,*)
      endif
!
      return
end


!----------------------------------------------------------------------
  subroutine calcrdipole(dipmat,work,dmtrx,nproc,myrank, &
&                        mpi_comm,datajob,datamol,databasis,datacomp)
!----------------------------------------------------------------------
!
! Driver of dipole moment calculation for closed-shell
!
      use modparam, only : todebye
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(in) :: datajob
      type(typemol),intent(inout) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(inout) :: datacomp
      integer,intent(in) :: nproc, myrank, mpi_comm
      integer :: iatom, ii
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: dmtrx((databasis%nao*(databasis%nao+1))/2)
      real(8),intent(out) :: dipmat((databasis%nao*(databasis%nao+1))/2,3)
      real(8),intent(out) :: work((databasis%nao*(databasis%nao+1))/2,3)
      real(8) :: dipcenter(3), tridot
      real(8) :: xdipplus, ydipplus, zdipplus, xdipminus, ydipminus, zdipminus
!
! Nuclear part
!
      xdipplus= zero
      ydipplus= zero
      zdipplus= zero
!
      do iatom= 1,datamol%natom
        xdipplus= xdipplus+datamol%coord(1,iatom)*datamol%znuc(iatom)
        ydipplus= ydipplus+datamol%coord(2,iatom)*datamol%znuc(iatom)
        zdipplus= zdipplus+datamol%coord(3,iatom)*datamol%znuc(iatom)
      enddo
!
! Electron part
!
      dipcenter(:)= zero
!
      call calcmatdipole(dipmat,work,dipcenter,nproc,myrank,mpi_comm,datajob,datamol,databasis)
!
      xdipminus=-tridot(dmtrx,dipmat(1,1),databasis%nao)
      ydipminus=-tridot(dmtrx,dipmat(1,2),databasis%nao)
      zdipminus=-tridot(dmtrx,dipmat(1,3),databasis%nao)
!
! Sum Nuclear and Electron parts
!
      datamol%dipole(1)=(xdipplus+xdipminus)
      datamol%dipole(2)=(ydipplus+ydipminus)
      datamol%dipole(3)=(zdipplus+zdipminus)
      datamol%dipole(4)= sqrt(datamol%dipole(1)**2+datamol%dipole(2)**2+datamol%dipole(3)**2)
!
      if(datacomp%master) then
        write(datacomp%iout,'(" ----------------------------------------------")')
        write(datacomp%iout,'("               Dipole Momemt (Debye)")')
        write(datacomp%iout,'("        X          Y          Z       Total")')
        write(datacomp%iout,'(" ----------------------------------------------")')
        write(datacomp%iout,'(1x,4f11.4)') (datamol%dipole(ii)*todebye, ii=1,4)
        write(datacomp%iout,'(" ----------------------------------------------",/)')
      endif
!
      return
end


!-----------------------------------------------------------------
  subroutine calcroctupole(dipmat,quadpmat,octpmat,work,dmtrx, &
&                          nproc,myrank,mpi_comm, &
&                          datajob,datamol,databasis,datacomp)
!-----------------------------------------------------------------
!
! Driver of dipole, quadrupole, and octupole moment calculation for closed-shell
!
      use modparam, only : todebye, toang
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(in) :: datajob
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(in) :: datacomp
      integer,intent(in) :: nproc, myrank, mpi_comm
      integer :: iatom, ii
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, two=2.0D+00, three=3.0D+00
      real(8),parameter :: four=4.0D+00, five=5.0D+00
      real(8),intent(in) :: dmtrx((databasis%nao*(databasis%nao+1))/2)
      real(8),intent(out) :: dipmat((databasis%nao*(databasis%nao+1))/2,3)
      real(8),intent(out) :: quadpmat((databasis%nao*(databasis%nao+1))/2,6)
      real(8),intent(out) :: octpmat((databasis%nao*(databasis%nao+1))/2,10)
      real(8),intent(out) :: work((databasis%nao*(databasis%nao+1))/2*10)
      real(8) :: dipcenter(3), xdip, ydip, zdip, totaldip, tridot
      real(8) :: xdipplus, ydipplus, zdipplus, xdipminus, ydipminus, zdipminus
      real(8) :: xxquadpplus, xyquadpplus, xzquadpplus, yyquadpplus, yzquadpplus, zzquadpplus
      real(8) :: xxquadpminus, xyquadpminus, xzquadpminus, yyquadpminus, yzquadpminus, zzquadpminus
      real(8) :: xxquadp, xyquadp, xzquadp, yyquadp, yzquadp, zzquadp, quadp(6), facp
      real(8) :: xxxoctpplus, xxyoctpplus, xxzoctpplus, xyyoctpplus, xyzoctpplus, xzzoctpplus
      real(8) :: yyyoctpplus, yyzoctpplus, yzzoctpplus, zzzoctpplus
      real(8) :: xxxoctpminus, xxyoctpminus, xxzoctpminus, xyyoctpminus, xyzoctpminus, xzzoctpminus
      real(8) :: yyyoctpminus, yyzoctpminus, yzzoctpminus, zzzoctpminus
      real(8) :: xxxoctp, xxyoctp, xxzoctp, xyyoctp, xyzoctp, xzzoctp
      real(8) :: yyyoctp, yyzoctp, yzzoctp, zzzoctp, octp(10), faco
!
! Nuclear part
!
      xdipplus= zero
      ydipplus= zero
      zdipplus= zero
      xxquadpplus= zero
      xyquadpplus= zero
      xzquadpplus= zero
      yyquadpplus= zero
      yzquadpplus= zero
      zzquadpplus= zero
      xxxoctpplus= zero
      xxyoctpplus= zero
      xxzoctpplus= zero
      xyyoctpplus= zero
      xyzoctpplus= zero
      xzzoctpplus= zero
      yyyoctpplus= zero
      yyzoctpplus= zero
      yzzoctpplus= zero
      zzzoctpplus= zero
!
      do iatom= 1,datamol%natom
        xdipplus= xdipplus+datamol%coord(1,iatom)*datamol%znuc(iatom)
        ydipplus= ydipplus+datamol%coord(2,iatom)*datamol%znuc(iatom)
        zdipplus= zdipplus+datamol%coord(3,iatom)*datamol%znuc(iatom)
        xxquadpplus= xxquadpplus+datamol%coord(1,iatom)*datamol%coord(1,iatom)*datamol%znuc(iatom)
        xyquadpplus= xyquadpplus+datamol%coord(1,iatom)*datamol%coord(2,iatom)*datamol%znuc(iatom)
        xzquadpplus= xzquadpplus+datamol%coord(1,iatom)*datamol%coord(3,iatom)*datamol%znuc(iatom)
        yyquadpplus= yyquadpplus+datamol%coord(2,iatom)*datamol%coord(2,iatom)*datamol%znuc(iatom)
        yzquadpplus= yzquadpplus+datamol%coord(2,iatom)*datamol%coord(3,iatom)*datamol%znuc(iatom)
        zzquadpplus= zzquadpplus+datamol%coord(3,iatom)*datamol%coord(3,iatom)*datamol%znuc(iatom)
        xxxoctpplus= xxxoctpplus+datamol%coord(1,iatom)*datamol%coord(1,iatom) &
&                   *datamol%coord(1,iatom)*datamol%znuc(iatom)
        xxyoctpplus= xxyoctpplus+datamol%coord(1,iatom)*datamol%coord(1,iatom) &
&                   *datamol%coord(2,iatom)*datamol%znuc(iatom)
        xxzoctpplus= xxzoctpplus+datamol%coord(1,iatom)*datamol%coord(1,iatom) &
&                   *datamol%coord(3,iatom)*datamol%znuc(iatom)
        xyyoctpplus= xyyoctpplus+datamol%coord(1,iatom)*datamol%coord(2,iatom) &
&                   *datamol%coord(2,iatom)*datamol%znuc(iatom)
        xyzoctpplus= xyzoctpplus+datamol%coord(1,iatom)*datamol%coord(2,iatom) &
&                   *datamol%coord(3,iatom)*datamol%znuc(iatom)
        xzzoctpplus= xzzoctpplus+datamol%coord(1,iatom)*datamol%coord(3,iatom) &
&                   *datamol%coord(3,iatom)*datamol%znuc(iatom)
        yyyoctpplus= yyyoctpplus+datamol%coord(2,iatom)*datamol%coord(2,iatom) &
&                   *datamol%coord(2,iatom)*datamol%znuc(iatom)
        yyzoctpplus= yyzoctpplus+datamol%coord(2,iatom)*datamol%coord(2,iatom) &
&                   *datamol%coord(3,iatom)*datamol%znuc(iatom)
        yzzoctpplus= yzzoctpplus+datamol%coord(2,iatom)*datamol%coord(3,iatom) &
&                   *datamol%coord(3,iatom)*datamol%znuc(iatom)
        zzzoctpplus= zzzoctpplus+datamol%coord(3,iatom)*datamol%coord(3,iatom) &
&                   *datamol%coord(3,iatom)*datamol%znuc(iatom)
      enddo
!
! Electron part
!
      dipcenter(:)= zero
!
      call calcmatoctupole(dipmat,quadpmat,octpmat,work,dipcenter,nproc,myrank,mpi_comm, &
&                          datajob,datamol,databasis)
!
      xdipminus=-tridot(dmtrx,dipmat(1,1),databasis%nao)
      ydipminus=-tridot(dmtrx,dipmat(1,2),databasis%nao)
      zdipminus=-tridot(dmtrx,dipmat(1,3),databasis%nao)
      xxquadpminus=-tridot(dmtrx,quadpmat(1,1),databasis%nao)
      xyquadpminus=-tridot(dmtrx,quadpmat(1,2),databasis%nao)
      xzquadpminus=-tridot(dmtrx,quadpmat(1,3),databasis%nao)
      yyquadpminus=-tridot(dmtrx,quadpmat(1,4),databasis%nao)
      yzquadpminus=-tridot(dmtrx,quadpmat(1,5),databasis%nao)
      zzquadpminus=-tridot(dmtrx,quadpmat(1,6),databasis%nao)
      xxxoctpminus=-tridot(dmtrx,octpmat(1, 1),databasis%nao)
      xxyoctpminus=-tridot(dmtrx,octpmat(1, 2),databasis%nao)
      xxzoctpminus=-tridot(dmtrx,octpmat(1, 3),databasis%nao)
      xyyoctpminus=-tridot(dmtrx,octpmat(1, 4),databasis%nao)
      xyzoctpminus=-tridot(dmtrx,octpmat(1, 5),databasis%nao)
      xzzoctpminus=-tridot(dmtrx,octpmat(1, 6),databasis%nao)
      yyyoctpminus=-tridot(dmtrx,octpmat(1, 7),databasis%nao)
      yyzoctpminus=-tridot(dmtrx,octpmat(1, 8),databasis%nao)
      yzzoctpminus=-tridot(dmtrx,octpmat(1, 9),databasis%nao)
      zzzoctpminus=-tridot(dmtrx,octpmat(1,10),databasis%nao)
!
! Sum Nuclear and Electron parts
!
      facp= todebye*toang*half
      faco= todebye*toang*toang*half
      xdip=(xdipplus+xdipminus)*todebye
      ydip=(ydipplus+ydipminus)*todebye
      zdip=(zdipplus+zdipminus)*todebye
      xxquadp=(xxquadpplus+xxquadpminus)*facp
      xyquadp=(xyquadpplus+xyquadpminus)*facp
      xzquadp=(xzquadpplus+xzquadpminus)*facp
      yyquadp=(yyquadpplus+yyquadpminus)*facp
      yzquadp=(yzquadpplus+yzquadpminus)*facp
      zzquadp=(zzquadpplus+zzquadpminus)*facp
      xxxoctp=(xxxoctpplus+xxxoctpminus)*faco
      xxyoctp=(xxyoctpplus+xxyoctpminus)*faco
      xxzoctp=(xxzoctpplus+xxzoctpminus)*faco
      xyyoctp=(xyyoctpplus+xyyoctpminus)*faco
      xyzoctp=(xyzoctpplus+xyzoctpminus)*faco
      xzzoctp=(xzzoctpplus+xzzoctpminus)*faco
      yyyoctp=(yyyoctpplus+yyyoctpminus)*faco
      yyzoctp=(yyzoctpplus+yyzoctpminus)*faco
      yzzoctp=(yzzoctpplus+yzzoctpminus)*faco
      zzzoctp=(zzzoctpplus+zzzoctpminus)*faco
      totaldip= sqrt(xdip*xdip+ydip*ydip+zdip*zdip)
!
      quadp(1)= xxquadp*two-yyquadp    -zzquadp
      quadp(2)= xyquadp*three
      quadp(3)= xzquadp*three
      quadp(4)=-xxquadp    +yyquadp*two-zzquadp
      quadp(5)= yzquadp*three
      quadp(6)=-xxquadp    -yyquadp    +zzquadp*two
      octp( 1)= xxxoctp*two  -xyyoctp*three-xzzoctp*three
      octp( 2)= xxyoctp*four -yyyoctp      -yzzoctp
      octp( 3)= xxzoctp*four -yyzoctp      -zzzoctp
      octp( 4)=-xxxoctp      +xyyoctp*four -xzzoctp
      octp( 5)= xyzoctp*five
      octp( 6)=-xxxoctp      -xyyoctp      +xzzoctp*four
      octp( 7)=-xxyoctp*three+yyyoctp*two  -yzzoctp*three
      octp( 8)=-xxzoctp      +yyzoctp*four -zzzoctp
      octp( 9)=-xxyoctp      -yyyoctp      +yzzoctp*four
      octp(10)=-xxzoctp*three-yyzoctp*three+zzzoctp*two
!
      if(datacomp%master) then
        write(datacomp%iout,'(" ----------------------------------------------")')
        write(datacomp%iout,'("               Dipole Momemt (Debye)")')
        write(datacomp%iout,'("        X          Y          Z       Total")')
        write(datacomp%iout,'(" ----------------------------------------------")')
        write(datacomp%iout,'(1x,4f11.4)')xdip, ydip, zdip, totaldip
        write(datacomp%iout,'(" ----------------------------------------------",/)')
        write(datacomp%iout,'(1x,68("-"))')
        write(datacomp%iout,'("               Quadrupole Momemt (Debye*Angstrom)")')
        write(datacomp%iout,'("        XX         XY         XZ         YY         YZ         ZZ")')
        write(datacomp%iout,'(1x,68("-"))')
        write(datacomp%iout,'(1x,6f11.4)')(quadp(ii),ii=1,6)
        write(datacomp%iout,'(1x,68("-"),/)')
        write(datacomp%iout,'(" ---------------------------------------------------------")')
        write(datacomp%iout,'("               Octupole Momemt (Debye*Angstrom^2)")')
        write(datacomp%iout,'("       XXX        XXY        XXZ       XYY        XYZ")')
        write(datacomp%iout,'(" ---------------------------------------------------------")')
        write(datacomp%iout,'(1x,5f11.4)')(octp(ii),ii=1,5)
        write(datacomp%iout,'(" ---------------------------------------------------------")')
        write(datacomp%iout,'("       XZZ        YYY        YYZ       YZZ        ZZZ")')
        write(datacomp%iout,'(" ---------------------------------------------------------")')
        write(datacomp%iout,'(1x,5f11.4)')(octp(ii),ii=6,10)
        write(datacomp%iout,'(" ---------------------------------------------------------",/)')
      endif
!
      return
end


!----------------------------------------------------------------------
  subroutine calcudipole(dipmat,work,dmtrxa,dmtrxb,nproc,myrank, &
&                        mpi_comm,datajob,datamol,databasis,datacomp)
!----------------------------------------------------------------------
!
! Driver of dipole moment calculation for open-shell
!
      use modparam, only : todebye
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(in) :: datajob
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(in) :: datacomp
      integer,intent(in) :: nproc, myrank, mpi_comm
      integer :: iatom
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: dmtrxa((databasis%nao*(databasis%nao+1))/2)
      real(8),intent(in) :: dmtrxb((databasis%nao*(databasis%nao+1))/2)
      real(8),intent(out) :: dipmat((databasis%nao*(databasis%nao+1))/2,3)
      real(8),intent(out) :: work((databasis%nao*(databasis%nao+1))/2,3)
      real(8) :: dipcenter(3), xdip, ydip, zdip, totaldip, tridot
      real(8) :: xdipplus, ydipplus, zdipplus, xdipminus, ydipminus, zdipminus
!
! Nuclear part
!
      xdipplus= zero
      ydipplus= zero
      zdipplus= zero
!
      do iatom= 1,datamol%natom
        xdipplus= xdipplus+datamol%coord(1,iatom)*datamol%znuc(iatom)
        ydipplus= ydipplus+datamol%coord(2,iatom)*datamol%znuc(iatom)
        zdipplus= zdipplus+datamol%coord(3,iatom)*datamol%znuc(iatom)
      enddo
!
! Electron part
!
      dipcenter(:)= zero
!
      call calcmatdipole(dipmat,work,dipcenter,nproc,myrank,mpi_comm,datajob,datamol,databasis)
!
      xdipminus=-tridot(dmtrxa,dipmat(1,1),databasis%nao)-tridot(dmtrxb,dipmat(1,1),databasis%nao)
      ydipminus=-tridot(dmtrxa,dipmat(1,2),databasis%nao)-tridot(dmtrxb,dipmat(1,2),databasis%nao)
      zdipminus=-tridot(dmtrxa,dipmat(1,3),databasis%nao)-tridot(dmtrxb,dipmat(1,3),databasis%nao)
!
! Sum Nuclear and Electron parts
!
      xdip=(xdipplus+xdipminus)*todebye
      ydip=(ydipplus+ydipminus)*todebye
      zdip=(zdipplus+zdipminus)*todebye
      totaldip= sqrt(xdip*xdip+ydip*ydip+zdip*zdip)
!
      if(datacomp%master) then
        write(datacomp%iout,'(" ----------------------------------------------")')
        write(datacomp%iout,'("               Dipole Momemt (Debye)")')
        write(datacomp%iout,'("        X          Y          Z       Total")')
        write(datacomp%iout,'(" ----------------------------------------------")')
        write(datacomp%iout,'(1x,4f11.4)')xdip, ydip, zdip, totaldip
        write(datacomp%iout,'(" ----------------------------------------------",/)')
      endif
!
      return
end


!-------------------------------------------------------------------------
  subroutine calcuoctupole(dipmat,quadpmat,octpmat,work,dmtrxa,dmtrxb, &
&                          nproc,myrank,mpi_comm, &
&                          datajob,datamol,databasis,datacomp)
!-------------------------------------------------------------------------
!
! Driver of dipole, quadrupole, and octupole moment calculation for open-shell
!
      use modparam, only : todebye, toang
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob),intent(in) :: datajob
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(in) :: datacomp
      integer,intent(in) :: nproc, myrank, mpi_comm
      integer :: iatom, ii
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, two=2.0D+00, three=3.0D+00
      real(8),parameter :: four=4.0D+00, five=5.0D+00
      real(8),intent(in) :: dmtrxa((databasis%nao*(databasis%nao+1))/2)
      real(8),intent(in) :: dmtrxb((databasis%nao*(databasis%nao+1))/2)
      real(8),intent(out) :: dipmat((databasis%nao*(databasis%nao+1))/2,3)
      real(8),intent(out) :: quadpmat((databasis%nao*(databasis%nao+1))/2,6)
      real(8),intent(out) :: octpmat((databasis%nao*(databasis%nao+1))/2,10)
      real(8),intent(out) :: work((databasis%nao*(databasis%nao+1))/2*10)
      real(8) :: dipcenter(3), xdip, ydip, zdip, totaldip, tridot
      real(8) :: xdipplus, ydipplus, zdipplus, xdipminus, ydipminus, zdipminus
      real(8) :: xxquadpplus, xyquadpplus, xzquadpplus, yyquadpplus, yzquadpplus, zzquadpplus
      real(8) :: xxquadpminus, xyquadpminus, xzquadpminus, yyquadpminus, yzquadpminus, zzquadpminus
      real(8) :: xxquadp, xyquadp, xzquadp, yyquadp, yzquadp, zzquadp, quadp(6), facp
      real(8) :: xxxoctpplus, xxyoctpplus, xxzoctpplus, xyyoctpplus, xyzoctpplus, xzzoctpplus
      real(8) :: yyyoctpplus, yyzoctpplus, yzzoctpplus, zzzoctpplus
      real(8) :: xxxoctpminus, xxyoctpminus, xxzoctpminus, xyyoctpminus, xyzoctpminus, xzzoctpminus
      real(8) :: yyyoctpminus, yyzoctpminus, yzzoctpminus, zzzoctpminus
      real(8) :: xxxoctp, xxyoctp, xxzoctp, xyyoctp, xyzoctp, xzzoctp
      real(8) :: yyyoctp, yyzoctp, yzzoctp, zzzoctp, octp(10), faco
!
! Nuclear part
!
      xdipplus= zero
      ydipplus= zero
      zdipplus= zero
      xxquadpplus= zero
      xyquadpplus= zero
      xzquadpplus= zero
      yyquadpplus= zero
      yzquadpplus= zero
      zzquadpplus= zero
      xxxoctpplus= zero
      xxyoctpplus= zero
      xxzoctpplus= zero
      xyyoctpplus= zero
      xyzoctpplus= zero
      xzzoctpplus= zero
      yyyoctpplus= zero
      yyzoctpplus= zero
      yzzoctpplus= zero
      zzzoctpplus= zero
!
      do iatom= 1,datamol%natom
        xdipplus= xdipplus+datamol%coord(1,iatom)*datamol%znuc(iatom)
        ydipplus= ydipplus+datamol%coord(2,iatom)*datamol%znuc(iatom)
        zdipplus= zdipplus+datamol%coord(3,iatom)*datamol%znuc(iatom)
        xxquadpplus= xxquadpplus+datamol%coord(1,iatom)*datamol%coord(1,iatom)*datamol%znuc(iatom)
        xyquadpplus= xyquadpplus+datamol%coord(1,iatom)*datamol%coord(2,iatom)*datamol%znuc(iatom)
        xzquadpplus= xzquadpplus+datamol%coord(1,iatom)*datamol%coord(3,iatom)*datamol%znuc(iatom)
        yyquadpplus= yyquadpplus+datamol%coord(2,iatom)*datamol%coord(2,iatom)*datamol%znuc(iatom)
        yzquadpplus= yzquadpplus+datamol%coord(2,iatom)*datamol%coord(3,iatom)*datamol%znuc(iatom)
        zzquadpplus= zzquadpplus+datamol%coord(3,iatom)*datamol%coord(3,iatom)*datamol%znuc(iatom)
        xxxoctpplus= xxxoctpplus+datamol%coord(1,iatom)*datamol%coord(1,iatom) &
&                   *datamol%coord(1,iatom)*datamol%znuc(iatom)
        xxyoctpplus= xxyoctpplus+datamol%coord(1,iatom)*datamol%coord(1,iatom) &
&                   *datamol%coord(2,iatom)*datamol%znuc(iatom)
        xxzoctpplus= xxzoctpplus+datamol%coord(1,iatom)*datamol%coord(1,iatom) &
&                   *datamol%coord(3,iatom)*datamol%znuc(iatom)
        xyyoctpplus= xyyoctpplus+datamol%coord(1,iatom)*datamol%coord(2,iatom) &
&                   *datamol%coord(2,iatom)*datamol%znuc(iatom)
        xyzoctpplus= xyzoctpplus+datamol%coord(1,iatom)*datamol%coord(2,iatom) &
&                   *datamol%coord(3,iatom)*datamol%znuc(iatom)
        xzzoctpplus= xzzoctpplus+datamol%coord(1,iatom)*datamol%coord(3,iatom) &
&                   *datamol%coord(3,iatom)*datamol%znuc(iatom)
        yyyoctpplus= yyyoctpplus+datamol%coord(2,iatom)*datamol%coord(2,iatom) &
&                   *datamol%coord(2,iatom)*datamol%znuc(iatom)
        yyzoctpplus= yyzoctpplus+datamol%coord(2,iatom)*datamol%coord(2,iatom) &
&                   *datamol%coord(3,iatom)*datamol%znuc(iatom)
        yzzoctpplus= yzzoctpplus+datamol%coord(2,iatom)*datamol%coord(3,iatom) &
&                   *datamol%coord(3,iatom)*datamol%znuc(iatom)
        zzzoctpplus= zzzoctpplus+datamol%coord(3,iatom)*datamol%coord(3,iatom) &
&                   *datamol%coord(3,iatom)*datamol%znuc(iatom)
      enddo
!
! Electron part
!
      dipcenter(:)= zero
!
      call calcmatoctupole(dipmat,quadpmat,octpmat,work,dipcenter,nproc,myrank,mpi_comm, &
&                          datajob,datamol,databasis)
!
      do ii= 1,databasis%nao*(databasis%nao+1)/2
        work(ii)= dmtrxa(ii)+dmtrxb(ii)
      enddo
!
      xdipminus=-tridot(work,dipmat(1,1),databasis%nao)
      ydipminus=-tridot(work,dipmat(1,2),databasis%nao)
      zdipminus=-tridot(work,dipmat(1,3),databasis%nao)
      xxquadpminus=-tridot(work,quadpmat(1,1),databasis%nao)
      xyquadpminus=-tridot(work,quadpmat(1,2),databasis%nao)
      xzquadpminus=-tridot(work,quadpmat(1,3),databasis%nao)
      yyquadpminus=-tridot(work,quadpmat(1,4),databasis%nao)
      yzquadpminus=-tridot(work,quadpmat(1,5),databasis%nao)
      zzquadpminus=-tridot(work,quadpmat(1,6),databasis%nao)
      xxxoctpminus=-tridot(work,octpmat(1, 1),databasis%nao)
      xxyoctpminus=-tridot(work,octpmat(1, 2),databasis%nao)
      xxzoctpminus=-tridot(work,octpmat(1, 3),databasis%nao)
      xyyoctpminus=-tridot(work,octpmat(1, 4),databasis%nao)
      xyzoctpminus=-tridot(work,octpmat(1, 5),databasis%nao)
      xzzoctpminus=-tridot(work,octpmat(1, 6),databasis%nao)
      yyyoctpminus=-tridot(work,octpmat(1, 7),databasis%nao)
      yyzoctpminus=-tridot(work,octpmat(1, 8),databasis%nao)
      yzzoctpminus=-tridot(work,octpmat(1, 9),databasis%nao)
      zzzoctpminus=-tridot(work,octpmat(1,10),databasis%nao)
!
! Sum Nuclear and Electron parts
!
      facp= todebye*toang*half
      faco= todebye*toang*toang*half
      xdip=(xdipplus+xdipminus)*todebye
      ydip=(ydipplus+ydipminus)*todebye
      zdip=(zdipplus+zdipminus)*todebye
      xxquadp=(xxquadpplus+xxquadpminus)*facp
      xyquadp=(xyquadpplus+xyquadpminus)*facp
      xzquadp=(xzquadpplus+xzquadpminus)*facp
      yyquadp=(yyquadpplus+yyquadpminus)*facp
      yzquadp=(yzquadpplus+yzquadpminus)*facp
      zzquadp=(zzquadpplus+zzquadpminus)*facp
      xxxoctp=(xxxoctpplus+xxxoctpminus)*faco
      xxyoctp=(xxyoctpplus+xxyoctpminus)*faco
      xxzoctp=(xxzoctpplus+xxzoctpminus)*faco
      xyyoctp=(xyyoctpplus+xyyoctpminus)*faco
      xyzoctp=(xyzoctpplus+xyzoctpminus)*faco
      xzzoctp=(xzzoctpplus+xzzoctpminus)*faco
      yyyoctp=(yyyoctpplus+yyyoctpminus)*faco
      yyzoctp=(yyzoctpplus+yyzoctpminus)*faco
      yzzoctp=(yzzoctpplus+yzzoctpminus)*faco
      zzzoctp=(zzzoctpplus+zzzoctpminus)*faco
      totaldip= sqrt(xdip*xdip+ydip*ydip+zdip*zdip)
!
      quadp(1)= xxquadp*two-yyquadp    -zzquadp
      quadp(2)= xyquadp*three
      quadp(3)= xzquadp*three
      quadp(4)=-xxquadp    +yyquadp*two-zzquadp
      quadp(5)= yzquadp*three
      quadp(6)=-xxquadp    -yyquadp    +zzquadp*two
      octp( 1)= xxxoctp*two  -xyyoctp*three-xzzoctp*three
      octp( 2)= xxyoctp*four -yyyoctp      -yzzoctp
      octp( 3)= xxzoctp*four -yyzoctp      -zzzoctp
      octp( 4)=-xxxoctp      +xyyoctp*four -xzzoctp
      octp( 5)= xyzoctp*five
      octp( 6)=-xxxoctp      -xyyoctp      +xzzoctp*four
      octp( 7)=-xxyoctp*three+yyyoctp*two  -yzzoctp*three
      octp( 8)=-xxzoctp      +yyzoctp*four -zzzoctp
      octp( 9)=-xxyoctp      -yyyoctp      +yzzoctp*four
      octp(10)=-xxzoctp*three-yyzoctp*three+zzzoctp*two
!
      if(datacomp%master) then
        write(datacomp%iout,'(" ----------------------------------------------")')
        write(datacomp%iout,'("               Dipole Momemt (Debye)")')
        write(datacomp%iout,'("        X          Y          Z       Total")')
        write(datacomp%iout,'(" ----------------------------------------------")')
        write(datacomp%iout,'(1x,4f11.4)')xdip, ydip, zdip, totaldip
        write(datacomp%iout,'(" ----------------------------------------------",/)')
        write(datacomp%iout,'(1x,68("-"))')
        write(datacomp%iout,'("               Quadrupole Momemt (Debye*Angstrom)")')
        write(datacomp%iout,'("        XX         XY         XZ         YY         YZ         ZZ")')
        write(datacomp%iout,'(1x,68("-"))')
        write(datacomp%iout,'(1x,6f11.4)')(quadp(ii),ii=1,6)
        write(datacomp%iout,'(1x,68("-"),/)')
        write(datacomp%iout,'(" ---------------------------------------------------------")')
        write(datacomp%iout,'("               Octupole Momemt (Debye*Angstrom^2)")')
        write(datacomp%iout,'("       XXX        XXY        XXZ       XYY        XYZ")')
        write(datacomp%iout,'(" ---------------------------------------------------------")')
        write(datacomp%iout,'(1x,5f11.4)')(octp(ii),ii=1,5)
        write(datacomp%iout,'(" ---------------------------------------------------------")')
        write(datacomp%iout,'("       XZZ        YYY        YYZ       YZZ        ZZZ")')
        write(datacomp%iout,'(" ---------------------------------------------------------")')
        write(datacomp%iout,'(1x,5f11.4)')(octp(ii),ii=6,10)
        write(datacomp%iout,'(" ---------------------------------------------------------",/)')
      endif
!
      return
end


!---------------------------------------------------------------------
  subroutine calcrnpa(dmtrx,fock,overlap,datamol,databasis,datacomp)
!---------------------------------------------------------------------
!
! Driver of closed-shell Natural Population Analysis calculation
!
      use modtype, only : typemol, typebasis, typecomp
      implicit none
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(in) :: datacomp
      integer :: maxang, nao, nao2, maxsize
      integer,allocatable :: infobasis(:,:,:), infonmb(:,:,:), list1(:), list2(:)
      real(8),intent(in) :: dmtrx(databasis%nao*(databasis%nao+1)/2)
      real(8),intent(in) :: fock(databasis%nao*(databasis%nao+1)/2)
      real(8),intent(in) :: overlap(databasis%nao*(databasis%nao+1)/2)
      real(8),allocatable :: pnao(:), snao(:), trans(:), work1(:), work2(:), work3(:)
      real(8),allocatable :: wnao(:)
!
      maxang= maxval(databasis%mtype(1:databasis%nshell))
      nao   = databasis%nao
      nao2  = nao*nao
!
! Allocate arrays
!
      call memset(nao2*6+nao*3+3*(maxang+8)*datamol%natom,datacomp)
      allocate(pnao(nao2),snao(nao2),trans(nao2),work1(nao2),work2(nao2), &
&              work3(nao2),wnao(nao),infobasis(3,maxang+1,datamol%natom), &
&              infonmb(3,7,datamol%natom),list1(nao),list2(nao))
!
! Calculate Natural Population
!
      call calcnpa(dmtrx,overlap,pnao,snao,trans,work1,work2,work3,wnao, &
&                  infobasis,infonmb,list1,list2,maxang,maxsize,datamol,databasis,datacomp)
!
! Transform fock matrix and print result
!
      call printrnpa(pnao,trans,fock,work1,work2,infonmb,infobasis, &
&                    maxang,maxsize,datamol,databasis,datacomp)
!
! Deallocate arrays
!
      call memunset(nao2*6+nao*3+3*(maxang+8)*datamol%natom,datacomp)
      deallocate(pnao,snao,trans,work1,work2, &
&                work3,wnao,infobasis, &
&                infonmb,list1,list2)
!
      return
end


!------------------------------------------------------------------------------------
  subroutine calcunpa(dmtrxa,dmtrxb,focka,fockb,overlap,datamol,databasis,datacomp)
!------------------------------------------------------------------------------------
!
! Driver of open-shell Natural Population Analysis calculation
!
      use modtype, only : typemol, typebasis, typecomp
      implicit none
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(in) :: datacomp
      integer :: maxang, nao, nao2, nao3, maxsize
      integer,allocatable :: infobasis(:,:,:), infonmb(:,:,:), list1(:), list2(:)
      real(8),intent(in) :: dmtrxa(databasis%nao*(databasis%nao+1)/2)
      real(8),intent(in) :: dmtrxb(databasis%nao*(databasis%nao+1)/2)
      real(8),intent(in) :: focka(databasis%nao*(databasis%nao+1)/2)
      real(8),intent(in) :: fockb(databasis%nao*(databasis%nao+1)/2)
      real(8),intent(in) :: overlap(databasis%nao*(databasis%nao+1)/2)
      real(8),allocatable :: pnao(:), snao(:), trans(:), work1(:), work2(:), work3(:), dmtrxab(:)
      real(8),allocatable :: wnao(:)
!
      maxang= maxval(databasis%mtype(1:databasis%nshell))
      nao   = databasis%nao
      nao2  = nao*nao
      nao3  = nao*(nao+1)/2
!
! Allocate arrays
!
      call memset(nao2*6+nao3+nao*3+3*(maxang+8)*datamol%natom,datacomp)
      allocate(pnao(nao2),snao(nao2),trans(nao2),work1(nao2),work2(nao2), &
&              work3(nao2),dmtrxab(nao3),wnao(nao),infobasis(3,maxang+1,datamol%natom), &
&              infonmb(3,7,datamol%natom),list1(nao),list2(nao))
!
! Calculate Natural Population
!
      dmtrxab(1:nao3)= dmtrxa(1:nao3)+dmtrxb(1:nao3)
      call calcnpa(dmtrxab,overlap,pnao,snao,trans,work1,work2,work3,wnao, &
&                  infobasis,infonmb,list1,list2,maxang,maxsize,datamol,databasis,datacomp)
!
! Transform fock matrix and print result
!
      call printunpa(focka,fockb,dmtrxa,dmtrxb,overlap,trans,pnao,work1,work2,work3,infonmb,infobasis, &
&                    maxang,maxsize,datamol,databasis,datacomp)
!
! Deallocate arrays
!
      call memunset(nao2*6+nao3+nao*3+3*(maxang+8)*datamol%natom,datacomp)
      deallocate(pnao,snao,trans,work1,work2, &
&                work3,dmtrxab,wnao,infobasis, &
&                infonmb,list1,list2)
!
      return
end


!----------------------------------------------------------------------------------------------
  subroutine calcnpa(dmtrx,overlap,pnao,snao,trans,work1,work2,work3,wnao, &
&                    infobasis,infonmb,list1,list2,maxang,maxsize,datamol,databasis,datacomp)
!----------------------------------------------------------------------------------------------
!
! Main driver of Natural Population Analysis calculation
!
      use modtype, only : typemol, typebasis, typecomp
      implicit none
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(in) :: datacomp
      type(typebasis) :: databasisnpa0, databasisnpa1, databasisnpa2
      integer :: nao, nao2, numnmb, numnrb, numnmbshell, iao, numlnrb, numsnrb
      integer,intent(in) :: maxang
      integer,intent(out) :: infobasis(3,maxang+1,datamol%natom), infonmb(3,7,datamol%natom), maxsize
      integer,intent(out) :: list1(databasis%nao), list2(databasis%nao)
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: dmtrx(databasis%nao*(databasis%nao+1)/2)
      real(8),intent(in) :: overlap(databasis%nao*(databasis%nao+1)/2)
      real(8),intent(out) :: pnao(databasis%nao**2), snao(databasis%nao**2)
      real(8),intent(out) :: trans(databasis%nao,databasis%nao), work1(databasis%nao**2)
      real(8),intent(out) :: work2(databasis%nao**2), work3(databasis%nao**2), wnao(databasis%nao)
!
      nao   = databasis%nao
      nao2  = nao*nao
!
! Set arrays
!
      infobasis(:,:,:)= 0
      infonmb(:,:,:)= 0
!
! Set full matrices of density, overlap, and NPA transformation
!
      call expand(dmtrx,pnao,nao)
      call expand2(overlap,snao,nao)
      call dsymm('L','U',nao,nao,one,pnao,nao,snao,nao,zero,work1,nao)
      call dsymm('L','U',nao,nao,one,snao,nao,work1,nao,zero,pnao,nao)
      trans(1:nao,1:nao)= zero
      do iao= 1,nao
        trans(iao,iao)= one
      enddo
!
! Transform basis functions from Cartesian to Spherial Harmonics if necessary
!
      if(.not.databasis%spher) then
        if(maxang >= 5) then
          write(*,'(" Error! This program supports up to g Cartesian functions", &
&                   " in Natural Population Analysis.",/, &
&                   " Set spher=.true. in Control section.")')
          call iabort
        endif
        call setnpaspher(pnao,snao,trans,work1,databasis,databasisnpa0,datacomp)
      else
        databasisnpa0= databasis
      endif
!
! Set numbers of core and valance Natural Minimal Basis (NMB) sets
!
      call setnmb(infonmb,numnmb,numnrb,numnmbshell,datamol,databasisnpa0,datacomp)
!
! Sort basis by atoms and angular momenta
!
      call sortbasisnpa1(pnao,snao,trans,work1,maxang,infobasis,datamol,databasisnpa0, &
&                        databasisnpa1)
      maxsize= maxval(infobasis(2,1:maxang+1,1:datamol%natom))
!
! Calculate pre-NAOs and seperate them into NMB and Natural Rydberg Basis (NRB) sets
!
      call calcprenao1(pnao,snao,trans,work1,work2,wnao,maxang,maxsize,infonmb,infobasis, &
&                      list1,numnmb,numnmbshell,datamol,databasisnpa1,databasisnpa2,datacomp)
!
! Orthogonalize NMB sets
!
      call orthonmb(pnao,snao,trans,wnao,work1,work2,work3,numnmb,databasisnpa2,datacomp)
!
! Orthogonalize NRB sets to NMB sets
!
      call orthonrb1(pnao,snao,trans,work1,work2,numnmb,numnrb,databasisnpa2)
!
! Calculate pre-NAOs of NRB sets
!
       call calcprenao2(pnao,snao,trans,work1,work2,wnao,maxang,maxsize,infobasis,list1, &
&                       list2,numnmb,numnmbshell,numlnrb,numsnrb,datamol,databasisnpa2,datacomp)
!
! Orthogonalize NRB sets with large w
!
      call orthonrb2(pnao,snao,trans,wnao,work1,work2,work3,numnmb,numlnrb, &
&                    databasisnpa2,datacomp)
!
! Orthogonalize NRB sets with small w
!
      call orthonrb3(pnao,snao,trans,work1,work2,work3,wnao,numnmb,numnrb,numlnrb,numsnrb, &
&                    databasisnpa2,datacomp)
!
! Sort basis by atoms and angular momenta
!
      call sortbasisnpa2(pnao,snao,trans,work1,list2,databasisnpa2)
!
! Calculate pre-NAOs
!
      call calcprenao3(pnao,snao,trans,work1,work2,maxang,maxsize,infobasis, &
&                      datamol,databasisnpa1,datacomp)
!
      return
end


!--------------------------------------------------------------------------------
  subroutine setnpaspher(pnao,snao,trans,work,databasis,databasisnpa0,datacomp)
!--------------------------------------------------------------------------------
!
! Transform basis functions from Cartesian to Spherical Harminics
!
      use modtype, only : typebasis, typecomp
      implicit none
      type(typebasis),intent(in) :: databasis
      type(typebasis),intent(out) :: databasisnpa0
      type(typecomp),intent(in) :: datacomp
      integer :: newshell, ishell, loctrans, nao, ii, jj
      integer :: intftrans(10,10), intgtrans(15,15)
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, halfm=-0.5D+00, one=1.0D+00
      real(8),parameter :: sqrtfifth=4.472135954999579D-01, sqrtfifteen=2.581988897471611D-01
      real(8),parameter :: sqrt3h=8.660254037844386D-01, sqrt3hm=-8.660254037844386D-01
      real(8),parameter :: sqrtseventh=3.779644730092272D-01, sqrtinv35=1.690308509457033D-01
      real(8),parameter :: sqrt3inv35=2.927700218845599D-01
      real(8),intent(inout) :: pnao(databasis%nao,databasis%nao)
      real(8),intent(inout) :: snao(databasis%nao,databasis%nao)
      real(8),intent(inout) :: trans(databasis%nao,databasis%nao)
      real(8),intent(out) :: work(databasis%nao,databasis%nao)
      real(8) :: dtrans(6,6), ffac1(10), ffac2(10), gfac1(15), gfac2(15)
      data dtrans / sqrtfifth, zero, zero, sqrtfifth, zero, sqrtfifth, &  ! s
&                   zero     , one , zero, zero     , zero, zero     , &  ! d-2
&                   halfm    , zero, zero, halfm    , zero, one      , &  ! d-1
&                   zero     , zero, zero, zero     , one , zero     , &  ! d0
&                   zero     , zero, one , zero     , zero, zero     , &  ! d+1
&                   sqrt3h   , zero, zero, sqrt3hm  , zero, zero     /    ! d+2
      data ffac1 / 6.546536707079771D-01, 6.546536707079771D-01, & ! sqrt(3/7) , sqrt(3/7)
&                  6.546536707079771D-01, 7.905694150420948D-01, & ! sqrt(3/7) , sqrt(5/2)/2 
&                  3.872983346207416D+00, 6.123724356957945D-01, & ! sqrt(15)  , sqrt(3/2)/2 
&                  5.000000000000000D-01, 6.123724356957945D-01, & ! (1/2)     , sqrt(3/2)/2 
&                  1.936491673103708D+00, 7.905694150420948D-01/   ! sqrt(15)/2, sqrt(5/2)/2
      data ffac2 / one, sqrtfifth, sqrtfifth, sqrtfifth, sqrtfifteen, &
&                  sqrtfifth, one, sqrtfifth, sqrtfifth, one/
      data gfac1 / 3.333333333333333D-01, 1.290994448735805D+00, & ! (1/3)     , sqrt(5/3) 
&                  1.290994448735805D+00, 3.726779962499649D-01, & ! sqrt(5/3) , sqrt(5)/6
&                  1.290994448735805D+00, 6.454972243679028D-01, & ! sqrt(5/3) , sqrt(5/12)
&                  2.958039891549808D+00, 2.091650066335188D+00, & ! sqrt(35)/2, sqrt(35/2)/2
&                  1.118033988749894D+00, 7.905694150420948D-01, & ! sqrt(5)/2 , sqrt(5/2)/2
&                  1.250000000000000D-01, 7.905694150420948D-01, & ! (1/8)     , sqrt(5/2)/2
&                  5.590169943749474D-01, 2.091650066335188D+00, & ! sqrt(5)/4 , sqrt(35/2)/2 
&                  7.395099728874520D-01/                          ! sqrt(35)/80
      data gfac2 / one, sqrtseventh, sqrtseventh, sqrt3inv35, sqrtinv35, &
&                  sqrt3inv35, sqrtseventh, sqrtinv35, sqrtinv35, sqrtseventh,&
&                  one, sqrtseventh, sqrt3inv35, sqrtseventh, one/
      data intftrans / 1,  0,  0,  1,  0,    1,  0,  0,  0,  0, & ! p-1
&                      0,  1,  0,  0,  0,    0,  1,  0,  1,  0, & ! p0
&                      0,  0,  1,  0,  0,    0,  0,  1,  0,  1, & ! p+1
&                      0,  3,  0,  0,  0,    0, -1,  0,  0,  0, & ! f-3
&                      0,  0,  0,  0,  1,    0,  0,  0,  0,  0, & ! f-2
&                      0, -1,  0,  0,  0,    0, -1,  0,  4,  0, & ! f-1
&                      0,  0, -3,  0,  0,    0,  0, -3,  0,  2, & ! f0
&                     -1,  0,  0, -1,  0,    4,  0,  0,  0,  0, & ! f+1
&                      0,  0,  1,  0,  0,    0,  0, -1,  0,  0, & ! f+2
&                      1,  0,  0, -3,  0,    0,  0,  0,  0,  0/   ! f+3
      data intgtrans / 1,  0,  0,  2,  0,    2,  0,  0,  0,  0,    1,  0,  2,  0,  1, & ! s
&                      0,  1,  0,  0,  0,    0,  1,  0,  1,  0,    0,  0,  0,  0,  0, & ! d-2
&                      0,  0,  0,  0,  1,    0,  0,  0,  0,  0,    0,  1,  0,  1,  0, & ! d-1
&                     -1,  0,  0, -2,  0,    1,  0,  0,  0,  0,   -1,  0,  1,  0,  2, & ! d0
&                      0,  0,  1,  0,  0,    0,  0,  1,  0,  1,    0,  0,  0,  0,  0, & ! d+1
&                      1,  0,  0,  0,  0,    1,  0,  0,  0,  0,   -1,  0, -1,  0,  0, & ! d+2
&                      0,  1,  0,  0,  0,    0, -1,  0,  0,  0,    0,  0,  0,  0,  0, & ! g-4
&                      0,  0,  0,  0,  3,    0,  0,  0,  0,  0,    0, -1,  0,  0,  0, & ! g-3
&                      0, -1,  0,  0,  0,    0, -1,  0,  6,  0,    0,  0,  0,  0,  0, & ! g-2
&                      0,  0,  0,  0, -3,    0,  0,  0,  0,  0,    0, -3,  0,  4,  0, & ! g-1
&                      3,  0,  0,  6,  0,  -24,  0,  0,  0,  0,    3,  0,-24,  0,  8, & ! g0
&                      0,  0, -3,  0,  0,    0,  0, -3,  0,  4,    0,  0,  0,  0,  0, & ! g+1
&                     -1,  0,  0,  0,  0,    6,  0,  0,  0,  0,    1,  0, -6,  0,  0, & ! g+2
&                      0,  0,  1,  0,  0,    0,  0, -3,  0,  0,    0,  0,  0,  0,  0, & ! g+3
&                      1,  0,  0, -6,  0,    0,  0,  0,  0,  0,    1,  0,  0,  0,  0/   ! g+4
!
      newshell= 0
      do ishell= 1,databasis%nshell
        select case(databasis%mtype(ishell))
          case(0:1) ! s or p shell
            newshell= newshell+1
            databasisnpa0%mtype(newshell)  = databasis%mtype(ishell)
            databasisnpa0%mbf(newshell)    = databasis%mbf(ishell)
            databasisnpa0%locatom(newshell)= databasis%locatom(ishell)
            databasisnpa0%locbf(newshell)  = databasis%locbf(ishell)
          case(2)  ! d shell
            newshell= newshell+1
            databasisnpa0%mtype(newshell)  = 0
            databasisnpa0%mbf(newshell)    = 1
            databasisnpa0%locatom(newshell)= databasis%locatom(ishell)
            databasisnpa0%locbf(newshell)  = databasis%locbf(ishell)
            newshell= newshell+1
            databasisnpa0%mtype(newshell)  = 2
            databasisnpa0%mbf(newshell)    = 5
            databasisnpa0%locatom(newshell)= databasis%locatom(ishell)
            databasisnpa0%locbf(newshell)  = databasis%locbf(ishell)+1
            loctrans= databasis%locbf(ishell)
            trans(loctrans+1:loctrans+6,loctrans+1:loctrans+6)= dtrans(1:6,1:6)
          case(3)  ! f shell
            newshell= newshell+1
            databasisnpa0%mtype(newshell)  = 1
            databasisnpa0%mbf(newshell)    = 3
            databasisnpa0%locatom(newshell)= databasis%locatom(ishell)
            databasisnpa0%locbf(newshell)  = databasis%locbf(ishell)
            newshell= newshell+1
            databasisnpa0%mtype(newshell)  = 3
            databasisnpa0%mbf(newshell)    = 7
            databasisnpa0%locatom(newshell)= databasis%locatom(ishell)
            databasisnpa0%locbf(newshell)  = databasis%locbf(ishell)+3
            loctrans= databasis%locbf(ishell)
            do ii= 1,10
              do jj= 1,10
                trans(loctrans+jj,loctrans+ii)= dble(intftrans(jj,ii))*ffac1(ii)*ffac2(jj)
              enddo
            enddo
          case(4)  ! g shell
            newshell= newshell+1
            databasisnpa0%mtype(newshell)  = 0
            databasisnpa0%mbf(newshell)    = 1
            databasisnpa0%locatom(newshell)= databasis%locatom(ishell)
            databasisnpa0%locbf(newshell)  = databasis%locbf(ishell)
            newshell= newshell+1
            databasisnpa0%mtype(newshell)  = 2
            databasisnpa0%mbf(newshell)    = 5
            databasisnpa0%locatom(newshell)= databasis%locatom(ishell)
            databasisnpa0%locbf(newshell)  = databasis%locbf(ishell)+1
            newshell= newshell+1
            databasisnpa0%mtype(newshell)  = 4
            databasisnpa0%mbf(newshell)    = 9
            databasisnpa0%locatom(newshell)= databasis%locatom(ishell)
            databasisnpa0%locbf(newshell)  = databasis%locbf(ishell)+6
            loctrans= databasis%locbf(ishell)
            do ii= 1,15
              do jj= 1,15
                trans(loctrans+jj,loctrans+ii)= dble(intgtrans(jj,ii))*gfac1(ii)*gfac2(jj)
              enddo
            enddo
          case(5:)  ! h or higher shell
            if(datacomp%master) then
              write(*,'(" Error! This program supports up to g Cartesian functions", &
&                       " in Natural Population Analysis.")')
              call iabort
            endif
        end select
      enddo
      databasisnpa0%nao   = databasis%nao
      databasisnpa0%nshell= newshell
!
      nao= databasis%nao
      call dsymm('L','U',nao,nao,one,pnao,nao,trans,nao,zero,work,nao)
      call dgemm('T','N',nao,nao,nao,one,trans,nao,work,nao,zero,pnao,nao)
      call dsymm('L','U',nao,nao,one,snao,nao,trans,nao,zero,work,nao)
      call dgemm('T','N',nao,nao,nao,one,trans,nao,work,nao,zero,snao,nao)
!
      return
end


!------------------------------------------------------------------------------------------
  subroutine sortbasisnpa1(pnao,snao,trans,work,maxang,infobasis,datamol,databasisnpa0, &
&                          databasisnpa1)
!------------------------------------------------------------------------------------------
!
! Sort basis functions by atoms and angular momenta
!
      use modtype, only : typemol, typebasis
      implicit none
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasisnpa0
      type(typebasis),intent(out) :: databasisnpa1
      integer,intent(in) :: maxang
      integer,intent(out) :: infobasis(3,0:maxang,datamol%natom)
      integer :: jao, jloc, jshell, iatom, iang, ishell, iao, list(databasisnpa0%nao)
      real(8),intent(inout) :: pnao(databasisnpa0%nao,databasisnpa0%nao)
      real(8),intent(inout) :: snao(databasisnpa0%nao,databasisnpa0%nao)
      real(8),intent(inout) :: trans(databasisnpa0%nao,databasisnpa0%nao)
      real(8),intent(out) :: work(databasisnpa0%nao,databasisnpa0%nao)
!
      databasisnpa1%nshell= databasisnpa0%nshell
      databasisnpa1%nao= databasisnpa0%nao
!
! Make reordering list
!
      jao= 0
      jloc= 0
      jshell= 0
      do iatom= 1,datamol%natom
        do iang= 0,maxang
          do ishell= 1,databasisnpa0%nshell
            if((databasisnpa0%locatom(ishell) == iatom).and. &
&              (databasisnpa0%mtype(ishell) == iang)) then
              infobasis(1,iang,iatom)= infobasis(1,iang,iatom)+1
              infobasis(2,iang,iatom)= infobasis(2,iang,iatom)+2*iang+1
              jshell= jshell+1
!ishimura
!             databasisnpa1%locatom(jshell)= iatom
!             databasisnpa1%mtype(jshell)= iang
!             databasisnpa1%mbf(jshell)= databasisnpa0%mbf(ishell)
              databasisnpa1%locbf(jshell)= jloc
              jloc= jloc+2*iang+1
              do iao= 1,2*iang+1
                jao= jao+1
                list(jao)= databasisnpa0%locbf(ishell)+iao
              enddo
            endif
          enddo
        enddo
      enddo
!
! Reorder pnao, snao, and trans matrices
!
!$OMP parallel
!$OMP do
      do iao= 1,databasisnpa0%nao
        do jao= 1,databasisnpa0%nao
          work(iao,jao)= pnao(jao,list(iao))
        enddo
      enddo
!$OMP enddo
!$OMP do
      do iao= 1,databasisnpa0%nao
        do jao= 1,databasisnpa0%nao
          pnao(jao,iao)= work(jao,list(iao))
        enddo
      enddo
!$OMP enddo
!
!$OMP do
      do iao= 1,databasisnpa0%nao
        do jao= 1,databasisnpa0%nao
          work(iao,jao)= snao(jao,list(iao))
        enddo
      enddo
!$OMP enddo
!$OMP do
      do iao= 1,databasisnpa0%nao
        do jao= 1,databasisnpa0%nao
          snao(jao,iao)= work(jao,list(iao))
        enddo
      enddo
!$OMP enddo
!
!$OMP do
      do iao= 1,databasisnpa0%nao
        do jao= 1,databasisnpa0%nao
          work(jao,iao)= trans(jao,list(iao))
        enddo
      enddo
!$OMP enddo
!$OMP do
      do iao= 1,databasisnpa0%nao
        do jao= 1,databasisnpa0%nao
          trans(jao,iao)= work(jao,iao)
        enddo
      enddo
!$OMP enddo
!$OMP end parallel
!
      return
end


!-------------------------------------------------------------------------------------------------
  subroutine calcprenao1(pnao,snao,trans,worktrans,work,wnao,maxang,maxsize,infonmb,infobasis, &
&                        list1,numnmb,numnmbshell,datamol,databasisnpa1,databasisnpa2,datacomp)
!-------------------------------------------------------------------------------------------------
!
! Calculate pre-NAOs and separate them into Natural Minimal Basis (NMB) and 
! Natural Rydberg Basis (NRB) sets
!
      use modtype, only : typemol, typebasis, typecomp
      implicit none
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasisnpa1
      type(typebasis),intent(out) :: databasisnpa2
      type(typecomp),intent(in) :: datacomp
      integer,intent(in) :: maxang, maxsize, infonmb(3,0:6,datamol%natom), numnmb, numnmbshell
      integer,intent(inout) :: infobasis(3,0:maxang,datamol%natom)
      integer,intent(out) :: list1(databasisnpa1%nao)
      integer :: locshell, iatom, iang, ishell, jshell, iao, jao, numshell, locnmb, locnrb, nao
      integer :: locnmbshell, locnrbshell
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(inout) :: pnao(databasisnpa1%nao,databasisnpa1%nao)
      real(8),intent(inout) :: snao(databasisnpa1%nao,databasisnpa1%nao)
      real(8),intent(inout) :: trans(databasisnpa1%nao,databasisnpa1%nao)
      real(8),intent(out) :: worktrans(databasisnpa1%nao,databasisnpa1%nao)
      real(8),intent(out) :: work(databasisnpa1%nao,databasisnpa1%nao)
      real(8),intent(out) :: wnao(databasisnpa1%nao)
      real(8) :: pblock(maxsize,maxsize), sblock(maxsize,maxsize), wblock(maxsize)
      real(8) :: pelem
!
      nao= databasisnpa1%nao
      worktrans(1:nao,1:nao)= zero
      databasisnpa2= databasisnpa1
!
! Calculate pre-NAO transformation matrix
!
      locshell= 0
      locnmb= 0
      locnrb= numnmb
      locnmbshell= 0
      locnrbshell= numnmbshell
      do iatom= 1,datamol%natom
        do iang= 0,maxang
          numshell= infobasis(1,iang,iatom)
          if(numshell == 0) cycle
!
!   Construct and diagonalize P-Al and S-Al subblocks
!
          do ishell= 1,numshell
            do jshell= 1,ishell
              pelem= zero
              do iao= 1,2*iang+1
                pelem= pelem+pnao(databasisnpa1%locbf(locshell+jshell)+iao, &
&                                 databasisnpa1%locbf(locshell+ishell)+iao)
              enddo
              pblock(jshell,ishell)= pelem/dble(2*iang+1)
              sblock(jshell,ishell)= snao(databasisnpa1%locbf(locshell+jshell)+1, &
&                                         databasisnpa1%locbf(locshell+ishell)+1)
            enddo
          enddo
          call gendiag(1,'V','U',numshell,pblock,maxsize,sblock,maxsize, &
&                      wblock,datacomp)
!
!   Set NMB transformation
!
          do ishell= 1,infonmb(1,iang,iatom)
            do iao= 1,2*iang+1
              wnao(locnmb+iao)= wblock(numshell-ishell+1)
              list1(locnmb+iao)= databasisnpa1%locbf(locshell+ishell)+iao
            enddo
            do jshell= 1,numshell
              do jao= 1,2*iang+1
                worktrans(databasisnpa1%locbf(locshell+jshell)+jao, &
&                         locnmb+jao)= pblock(jshell,numshell-ishell+1)
              enddo
            enddo
!ishimura
!           databasisnpa2%locatom(locnmbshell+1)= iatom
!           databasisnpa2%mtype(locnmbshell+1)  = iang
!           databasisnpa2%mbf(locnmbshell+1)    = 2*iang+1
            databasisnpa2%locbf(locnmbshell+1)  = locnmb
            locnmb= locnmb+2*iang+1
            locnmbshell= locnmbshell+1
          enddo
!
!   Set NRB transformation
!
          do ishell= infonmb(1,iang,iatom)+1,numshell
            do iao= 1,2*iang+1
              wnao(locnrb+iao)= wblock(numshell-ishell+1)
              list1(locnrb+iao)= databasisnpa1%locbf(locshell+ishell)+iao
            enddo
            do jshell= 1,numshell
              do jao= 1,2*iang+1
                worktrans(databasisnpa1%locbf(locshell+jshell)+jao, &
&                         locnrb+jao)= pblock(jshell,numshell-ishell+1)
              enddo
            enddo
!ishimura
!           databasisnpa2%locatom(locnrbshell+1)= iatom
!           databasisnpa2%mtype(locnrbshell+1)  = iang
!           databasisnpa2%mbf(locnrbshell+1)    = 2*iang+1
            databasisnpa2%locbf(locnrbshell+1)  = locnrb
            locnrb= locnrb+2*iang+1
            locnrbshell= locnrbshell+1
            infobasis(3,iang,iatom)= infobasis(3,iang,iatom)+1
          enddo
          locshell= locshell+numshell
        enddo
      enddo
!
! Transform into pre-NAOs and separate into NMB and NRB sets
!
      call dsymm('L','U',nao,nao,one,pnao,nao,worktrans,nao,zero,work,nao)
      call dgemm('T','N',nao,nao,nao,one,worktrans,nao,work,nao,zero,pnao,nao)
      call dsymm('L','U',nao,nao,one,snao,nao,worktrans,nao,zero,work,nao)
      call dgemm('T','N',nao,nao,nao,one,worktrans,nao,work,nao,zero,snao,nao)
      call dgemm('N','N',nao,nao,nao,one,trans,nao,worktrans,nao,zero,work,nao)
      trans(1:nao,1:nao)= work(1:nao,1:nao)
!
      return
end


!--------------------------------------------------------------------------------------------------
  subroutine calcprenao2(pnao,snao,trans,worktrans,work,wnao,maxang,maxsize,infobasis,list1, &
&                        list2,numnmb,numnmbshell,numlnrb,numsnrb,datamol,databasisnpa2,datacomp)
!--------------------------------------------------------------------------------------------------
!
! Calculate pre-NAOs of NRB sets and separate them into 2 groups which have large or small Ws.
!
      use modtype, only : typemol, typebasis, typecomp
      implicit none
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasisnpa2
      type(typecomp),intent(in) :: datacomp
      integer,intent(in) :: maxang, maxsize, numnmb, numnmbshell
      integer,intent(in) :: infobasis(3,0:maxang,datamol%natom)
      integer,intent(in) :: list1(databasisnpa2%nao)
      integer,intent(out) :: list2(databasisnpa2%nao), numlnrb, numsnrb
      integer :: iatom, iang, ishell, jshell, iao, jao, numshell, nao
      integer :: locnrbshell, listtmp(databasisnpa2%nao), loclnrb, locsnrb
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, wthresh=1.0D-03
      real(8),intent(inout) :: pnao(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(inout) :: snao(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(inout) :: trans(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(out) :: worktrans(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(out) :: work(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(out) :: wnao(databasisnpa2%nao)
      real(8) :: pblock(maxsize,maxsize), sblock(maxsize,maxsize), wblock(maxsize)
      real(8) :: pelem, selem
!
      nao= databasisnpa2%nao
      worktrans(1:nao,1:nao)= zero
      work(1:nao,1:(nao-numnmb))= zero
      do iao= 1,numnmb
        worktrans(iao,iao)= one
        list2(iao)= list1(iao)
      enddo
!
! Calculate pre-NAO transformation matrix
!
      loclnrb= numnmb
      locsnrb= 0
      locnrbshell= numnmbshell
      numlnrb= 0
      numsnrb= 0
      do iatom= 1,datamol%natom
        do iang= 0,maxang
          numshell= infobasis(3,iang,iatom)
          if(numshell == 0) cycle
!
!   Construct and diagonalize P-Al and S-Al subblocks
!
          do ishell= 1,numshell
            do jshell= 1,ishell
              pelem= zero
              selem= zero
              do iao= 1,2*iang+1
                pelem= pelem+pnao(databasisnpa2%locbf(locnrbshell+jshell)+iao, &
&                                 databasisnpa2%locbf(locnrbshell+ishell)+iao)
                selem= selem+snao(databasisnpa2%locbf(locnrbshell+jshell)+iao, &
&                                 databasisnpa2%locbf(locnrbshell+ishell)+iao)
              enddo
              pblock(jshell,ishell)= pelem/dble(2*iang+1)
              sblock(jshell,ishell)= selem/dble(2*iang+1)
            enddo
          enddo
          call gendiag(1,'V','U',numshell,pblock,maxsize,sblock,maxsize, &
&                      wblock,datacomp)
!
!   Set NRB transformation
!
          do ishell= 1,numshell
            if(wblock(ishell) >= wthresh) then
              do iao= 1,2*iang+1
                wnao(loclnrb+iao)= wblock(ishell)
                list2(loclnrb+iao)= list1(databasisnpa2%locbf(locnrbshell+ishell)+iao)
              enddo
              do jshell= 1,numshell
                do jao= 1,2*iang+1
                  worktrans(databasisnpa2%locbf(locnrbshell+jshell)+jao, &
&                           loclnrb+jao)= pblock(jshell,ishell)
                enddo
              enddo
              loclnrb= loclnrb+2*iang+1
              numlnrb= numlnrb+2*iang+1
            else
              do iao= 1,2*iang+1
                listtmp(locsnrb+iao)= list1(databasisnpa2%locbf(locnrbshell+ishell)+iao)
              enddo
              do jshell= 1,numshell
                do jao= 1,2*iang+1
                  work(databasisnpa2%locbf(locnrbshell+jshell)+jao, &
&                       locsnrb+jao)= pblock(jshell,ishell)
                enddo
              enddo
              locsnrb= locsnrb+2*iang+1
              numsnrb= numsnrb+2*iang+1
            endif
          enddo
          locnrbshell= locnrbshell+numshell
        enddo
      enddo
!
      worktrans(numnmb+1:nao,numnmb+numlnrb+1:nao)= work(numnmb+1:nao,1:numsnrb)
      list2(numnmb+numlnrb+1:nao)= listtmp(1:numsnrb)
!
! Transform all NRB sets into pre-NAOs and separate them into large-w and small-w NRB groups
!
      call dsymm('L','U',nao,nao,one,pnao,nao,worktrans,nao,zero,work,nao)
      call dgemm('T','N',nao,nao,nao,one,worktrans,nao,work,nao,zero,pnao,nao)
      call dsymm('L','U',nao,nao,one,snao,nao,worktrans,nao,zero,work,nao)
      call dgemm('T','N',nao,nao,nao,one,worktrans,nao,work,nao,zero,snao,nao)
      call dgemm('N','N',nao,nao,nao,one,trans,nao,worktrans,nao,zero,work,nao)
      trans(1:nao,1:nao)= work(1:nao,1:nao)
!
      return
end


!------------------------------------------------------------------------------------
  subroutine calcprenao3(pnao,snao,trans,worktrans,work,maxang,maxsize,infobasis, &
&                        datamol,databasisnpa1,datacomp)
!------------------------------------------------------------------------------------
!
! Calculate pre-NAOs
!
      use modtype, only : typemol, typebasis, typecomp
      implicit none
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasisnpa1
      type(typecomp),intent(in) :: datacomp
      integer,intent(in) :: maxang, maxsize, infobasis(3,0:maxang,datamol%natom)
      integer :: locshell, iatom, iang, ishell, jshell, iao, jao, numshell, nao
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(inout) :: pnao(databasisnpa1%nao,databasisnpa1%nao)
      real(8),intent(inout) :: snao(databasisnpa1%nao,databasisnpa1%nao)
      real(8),intent(inout) :: trans(databasisnpa1%nao,databasisnpa1%nao)
      real(8),intent(out) :: worktrans(databasisnpa1%nao,databasisnpa1%nao)
      real(8),intent(out) :: work(databasisnpa1%nao,databasisnpa1%nao)
      real(8) :: pblock(maxsize,maxsize), sblock(maxsize,maxsize), wblock(maxsize)
      real(8) :: pelem
!
      nao= databasisnpa1%nao
      worktrans(1:nao,1:nao)= zero
!
! Calculate pre-NAO transformation matrix
!
      locshell= 0
      do iatom= 1,datamol%natom
        do iang= 0,maxang
          numshell= infobasis(1,iang,iatom)
          if(numshell == 0) cycle
!
!   Construct and diagonalize P-Al and S-Al subblocks
!
          do ishell= 1,numshell
            do jshell= 1,ishell
              pelem= zero
              do iao= 1,2*iang+1
                pelem= pelem+pnao(databasisnpa1%locbf(locshell+jshell)+iao, &
&                                 databasisnpa1%locbf(locshell+ishell)+iao)
              enddo
              pblock(jshell,ishell)= pelem/dble(2*iang+1)
!ishimura
!!!!!!!!!!not average ok?
              sblock(jshell,ishell)= snao(databasisnpa1%locbf(locshell+jshell)+1, &
&                                         databasisnpa1%locbf(locshell+ishell)+1)
            enddo
          enddo
          call gendiag(1,'V','U',numshell,pblock,maxsize,sblock,maxsize, &
&                      wblock,datacomp)
!
!   Set final transformation
!
          do ishell= 1,numshell
            do jshell= 1,numshell
              do jao= 1,2*iang+1
                worktrans(databasisnpa1%locbf(locshell+jshell)+jao, &
&                         databasisnpa1%locbf(locshell+ishell)+jao)= pblock(jshell,ishell)
              enddo
            enddo
          enddo
          locshell= locshell+numshell
        enddo
      enddo
!
! Transform into pre-NAOs and separate into NMB and NRB sets
!
      call dsymm('L','U',nao,nao,one,pnao,nao,worktrans,nao,zero,work,nao)
      call dgemm('T','N',nao,nao,nao,one,worktrans,nao,work,nao,zero,pnao,nao)
      call dgemm('N','N',nao,nao,nao,one,trans,nao,worktrans,nao,zero,work,nao)
      trans(1:nao,1:nao)= work(1:nao,1:nao)
!
      return
end


!--------------------------------------------------------------------------------------
  subroutine setnmb(infonmb,numnmb,numnrb,numnmbshell,datamol,databasisnpa0,datacomp)
!--------------------------------------------------------------------------------------
!
! Set numbers of core and valance Natural Minimal Basis (NMB) sets
!
!   infonmb(1,*,*) : number of NMBs
!   infonmb(2,*,*) : number of cores
!   infonmb(3,*,*) : number of frozen cores
!
      use modtype, only : typemol, typebasis, typecomp
      implicit none
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasisnpa0
      type(typecomp),intent(in) :: datacomp
      integer,intent(out) :: infonmb(3,0:6,datamol%natom), numnmb, numnrb, numnmbshell
      integer :: iatom, iang
!
      do iatom= 1,datamol%natom
        select case(datamol%numatomic(iatom))
! H - He
          case(1:2)
            infonmb(1,0,iatom)= 1
! Li - Be
          case(3:4)
            infonmb(1,0,iatom)= 2
            infonmb(2,0,iatom)= 1
! B - Ne
          case(5:10)
            infonmb(1,0,iatom)= 2
            infonmb(2,0,iatom)= 1
            infonmb(1,1,iatom)= 1
! Na - Mg
          case(11:12)
            infonmb(1,0,iatom)= 3
            infonmb(2,0,iatom)= 2
            infonmb(1,1,iatom)= 1
            infonmb(2,1,iatom)= 1
! Al - Ar
          case(13:18)
            infonmb(1,0,iatom)= 3
            infonmb(2,0,iatom)= 2
            infonmb(1,1,iatom)= 2
            infonmb(2,1,iatom)= 1
! K - Ca
          case(19:20)
            infonmb(1,0,iatom)= 4
            infonmb(2,0,iatom)= 3
            infonmb(1,1,iatom)= 2
            infonmb(2,1,iatom)= 2
! Sc - Zn
          case(21:30)
            infonmb(1,0,iatom)= 4
            infonmb(2,0,iatom)= 3
            infonmb(1,1,iatom)= 3
            infonmb(2,1,iatom)= 2
            infonmb(1,2,iatom)= 1
! Ga - Kr
          case(31:36)
            infonmb(1,0,iatom)= 4
            infonmb(2,0,iatom)= 3
            infonmb(1,1,iatom)= 3
            infonmb(2,1,iatom)= 2
            infonmb(1,2,iatom)= 1
            infonmb(2,2,iatom)= 1
! Rb - Sr
          case(37:38)
            infonmb(1,0,iatom)= 5
            infonmb(2,0,iatom)= 4
            infonmb(1,1,iatom)= 3
            infonmb(2,1,iatom)= 3
            infonmb(1,2,iatom)= 1
            infonmb(2,2,iatom)= 1
! Y - Cd
          case(39:48)
            infonmb(1,0,iatom)= 5
            infonmb(2,0,iatom)= 4
            infonmb(1,1,iatom)= 4
            infonmb(2,1,iatom)= 3
            infonmb(1,2,iatom)= 2
            infonmb(2,2,iatom)= 1
! In - Xe
          case(49:54)
            infonmb(1,0,iatom)= 5
            infonmb(2,0,iatom)= 4
            infonmb(1,1,iatom)= 4
            infonmb(2,1,iatom)= 3
            infonmb(1,2,iatom)= 2
            infonmb(2,2,iatom)= 2
! Cs - Ba
          case(55:56)
            infonmb(1,0,iatom)= 6
            infonmb(2,0,iatom)= 5
            infonmb(1,1,iatom)= 4
            infonmb(2,1,iatom)= 4
            infonmb(1,2,iatom)= 2
            infonmb(2,2,iatom)= 2
! La - Yb
          case(57:70)
            infonmb(1,0,iatom)= 6
            infonmb(2,0,iatom)= 5
            infonmb(1,1,iatom)= 5
            infonmb(2,1,iatom)= 4
            infonmb(1,2,iatom)= 3
            infonmb(2,2,iatom)= 2
            infonmb(1,3,iatom)= 1
! Lu - Hg
          case(71:80)
            infonmb(1,0,iatom)= 6
            infonmb(2,0,iatom)= 5
            infonmb(1,1,iatom)= 5
            infonmb(2,1,iatom)= 4
            infonmb(1,2,iatom)= 3
            infonmb(2,2,iatom)= 2
            infonmb(1,3,iatom)= 1
            infonmb(1,3,iatom)= 1
! Tl - Rn
          case(81:86)
            infonmb(1,0,iatom)= 6
            infonmb(2,0,iatom)= 5
            infonmb(1,1,iatom)= 5
            infonmb(2,1,iatom)= 4
            infonmb(1,2,iatom)= 3
            infonmb(2,2,iatom)= 3
            infonmb(1,3,iatom)= 1
            infonmb(1,3,iatom)= 1
          case(87:)
            if(datacomp%master) then
              write(*,'(" Error! Supportted atom in Natural Population Analysis is up to Rn.")')
              call iabort
            endif
        end select
        if(databasisnpa0%izcore(iatom) /= 0) then
          select case(databasisnpa0%izcore(iatom))
            case(2)
              infonmb(1,0,iatom)= infonmb(1,0,iatom)-1
              infonmb(2,0,iatom)= infonmb(2,0,iatom)-1
              infonmb(3,0,iatom)= 1
            case(10)
              infonmb(1,0,iatom)= infonmb(1,0,iatom)-2
              infonmb(2,0,iatom)= infonmb(2,0,iatom)-2
              infonmb(3,0,iatom)= 2
              infonmb(1,1,iatom)= infonmb(1,1,iatom)-1
              infonmb(2,1,iatom)= infonmb(2,1,iatom)-1
              infonmb(3,1,iatom)= 1
            case(18)
              infonmb(1,0,iatom)= infonmb(1,0,iatom)-3
              infonmb(2,0,iatom)= infonmb(2,0,iatom)-3
              infonmb(3,0,iatom)= 3
              infonmb(1,1,iatom)= infonmb(1,1,iatom)-2
              infonmb(2,1,iatom)= infonmb(2,1,iatom)-2
              infonmb(3,1,iatom)= 2
            case(28)
              infonmb(1,0,iatom)= infonmb(1,0,iatom)-3
              infonmb(2,0,iatom)= infonmb(2,0,iatom)-3
              infonmb(3,0,iatom)= 3
              infonmb(1,1,iatom)= infonmb(1,1,iatom)-2
              infonmb(2,1,iatom)= infonmb(2,1,iatom)-2
              infonmb(3,1,iatom)= 2
              infonmb(1,2,iatom)= infonmb(1,2,iatom)-1
              infonmb(2,2,iatom)= infonmb(2,2,iatom)-1
              infonmb(3,2,iatom)= 1
            case(36)
              infonmb(1,0,iatom)= infonmb(1,0,iatom)-4
              infonmb(2,0,iatom)= infonmb(2,0,iatom)-4
              infonmb(3,0,iatom)= 4
              infonmb(1,1,iatom)= infonmb(1,1,iatom)-3
              infonmb(2,1,iatom)= infonmb(2,1,iatom)-3
              infonmb(3,1,iatom)= 3
              infonmb(1,2,iatom)= infonmb(1,2,iatom)-1
              infonmb(2,2,iatom)= infonmb(2,2,iatom)-1
              infonmb(3,2,iatom)= 1
            case(46)
              infonmb(1,0,iatom)= infonmb(1,0,iatom)-4
              infonmb(2,0,iatom)= infonmb(2,0,iatom)-4
              infonmb(3,0,iatom)= 4
              infonmb(1,1,iatom)= infonmb(1,1,iatom)-3
              infonmb(2,1,iatom)= infonmb(2,1,iatom)-3
              infonmb(3,1,iatom)= 3
              infonmb(1,2,iatom)= infonmb(1,2,iatom)-2
              infonmb(2,2,iatom)= infonmb(2,2,iatom)-2
              infonmb(3,2,iatom)= 2
            case(60)
              infonmb(1,0,iatom)= infonmb(1,0,iatom)-4
              infonmb(2,0,iatom)= infonmb(2,0,iatom)-4
              infonmb(3,0,iatom)= 4
              infonmb(1,1,iatom)= infonmb(1,1,iatom)-3
              infonmb(2,1,iatom)= infonmb(2,1,iatom)-3
              infonmb(3,1,iatom)= 3
              infonmb(1,2,iatom)= infonmb(1,2,iatom)-2
              infonmb(2,2,iatom)= infonmb(2,2,iatom)-2
              infonmb(3,2,iatom)= 2
              infonmb(1,3,iatom)= infonmb(1,3,iatom)-1
              infonmb(2,3,iatom)= infonmb(2,3,iatom)-1
              infonmb(3,3,iatom)= 1
            case(68)
              infonmb(1,0,iatom)= infonmb(1,0,iatom)-5
              infonmb(2,0,iatom)= infonmb(2,0,iatom)-5
              infonmb(3,0,iatom)= 5
              infonmb(1,1,iatom)= infonmb(1,1,iatom)-4
              infonmb(2,1,iatom)= infonmb(2,1,iatom)-4
              infonmb(3,1,iatom)= 4
              infonmb(1,2,iatom)= infonmb(1,2,iatom)-2
              infonmb(2,2,iatom)= infonmb(2,2,iatom)-2
              infonmb(3,2,iatom)= 2
              infonmb(1,3,iatom)= infonmb(1,3,iatom)-1
              infonmb(2,3,iatom)= infonmb(2,3,iatom)-1
              infonmb(3,3,iatom)= 1
            case(78)
              infonmb(1,0,iatom)= infonmb(1,0,iatom)-5
              infonmb(2,0,iatom)= infonmb(2,0,iatom)-5
              infonmb(3,0,iatom)= 5
              infonmb(1,1,iatom)= infonmb(1,1,iatom)-4
              infonmb(2,1,iatom)= infonmb(2,1,iatom)-4
              infonmb(3,1,iatom)= 4
              infonmb(1,2,iatom)= infonmb(1,2,iatom)-3
              infonmb(2,2,iatom)= infonmb(2,2,iatom)-3
              infonmb(3,2,iatom)= 3
              infonmb(1,3,iatom)= infonmb(1,3,iatom)-1
              infonmb(2,3,iatom)= infonmb(2,3,iatom)-1
              infonmb(3,3,iatom)= 1
            case default 
              if(datacomp%master) then
                write(*,'(" Error! This ECP is not supported in Natural Population Analysis.")')
                call iabort
              endif
          end select
        endif
      enddo
!
      numnmb= 0
      numnmbshell= 0
      do iatom= 1,datamol%natom
        do iang= 0,3
          numnmb     = numnmb+infonmb(1,iang,iatom)*(2*iang+1)
          numnmbshell= numnmbshell+infonmb(1,iang,iatom)
        enddo
      enddo
      numnrb= databasisnpa0%nao-numnmb
!
      return
end


!--------------------------------------------------------------------------------------------
  subroutine orthonmb(pnao,snao,trans,wnao,work1,work2,work3,numnmb,databasisnpa2,datacomp)
!--------------------------------------------------------------------------------------------
!
! Orthogonalize NMB sets
!
      use modtype, only : typebasis, typecomp
      implicit none
      type(typebasis),intent(in) :: databasisnpa2
      type(typecomp),intent(in) :: datacomp
      integer,intent(in) :: numnmb
      integer :: nao, iao, jao
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: wnao(databasisnpa2%nao)
      real(8),intent(out) :: work1(databasisnpa2%nao,databasisnpa2%nao)  !  work1(nunmnb,numnmb)
      real(8),intent(out) :: work2(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(out) :: work3(databasisnpa2%nao,databasisnpa2%nao)  !  work3(nunmnb,numnmb)
      real(8),intent(inout) :: pnao(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(inout) :: snao(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(inout) :: trans(databasisnpa2%nao,databasisnpa2%nao)
!
! Calculate WSW matrix
!
      nao= databasisnpa2%nao
      do iao= 1,numnmb
        do jao= 1,numnmb
          work1(jao,iao)= snao(jao,iao)*wnao(iao)*wnao(jao)
        enddo
      enddo
!
! Diagonalize WSW matrix and calculate W(WSW)^-1/2
!
      call diag('V','U',numnmb,work1,nao,work3,datacomp)
      do iao= 1,numnmb
        do jao= 1,numnmb
          work2(jao,iao)= work1(jao,iao)/(sqrt(work3(iao,1)))
        enddo
      enddo
      call dgemm('N','T',numnmb,numnmb,numnmb,one,work2,nao,work1,nao,zero,work3,nao)
      do iao= 1,numnmb
        do jao= 1,numnmb
          work3(jao,iao)= work3(jao,iao)*wnao(jao)
        enddo
      enddo
!
! Occupancy-weighted symmetric orthonalization
!
      call dgemm('N','N',nao,numnmb,numnmb,one,pnao,nao,work3,nao,zero,work2,nao)
      pnao(1:nao,1:numnmb)= work2(1:nao,1:numnmb)
      call dgemm('T','N',numnmb,nao,numnmb,one,work3,nao,pnao,nao,zero,work2,nao)
      pnao(1:numnmb,1:nao)= work2(1:numnmb,1:nao)
      call dgemm('N','N',nao,numnmb,numnmb,one,snao,nao,work3,nao,zero,work2,nao)
      snao(1:nao,1:numnmb)= work2(1:nao,1:numnmb)
      call dgemm('T','N',numnmb,nao,numnmb,one,work3,nao,snao,nao,zero,work2,nao)
      snao(1:numnmb,1:nao)= work2(1:numnmb,1:nao)
      call dgemm('N','N',nao,numnmb,numnmb,one,trans,nao,work3,nao,zero,work2,nao)
      trans(1:nao,1:numnmb)= work2(1:nao,1:numnmb)
!
      return
end


!-----------------------------------------------------------------------------------
  subroutine orthonrb1(pnao,snao,trans,worktrans,work,numnmb,numnrb,databasisnpa2)
!-----------------------------------------------------------------------------------
!
! Orthogonalize NRB sets to NMB sets by Gram-Schmidt process
!
      use modtype, only : typebasis
      implicit none
      type(typebasis),intent(in) :: databasisnpa2
      integer,intent(in) :: numnmb, numnrb
      integer :: nao, iao, jao
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(out) :: worktrans(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(out) :: work(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(inout) :: pnao(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(inout) :: snao(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(inout) :: trans(databasisnpa2%nao,databasisnpa2%nao)
!
      if(numnrb == 0) return
!
      nao= databasisnpa2%nao
      worktrans(1:nao,1:nao)= zero
!
      do iao= 1,nao
        worktrans(iao,iao)= one
      enddo 
      do iao= numnmb+1,nao
        do jao= 1,numnmb
          worktrans(iao,iao)= worktrans(iao,iao)-snao(jao,iao)*snao(jao,iao)
        enddo 
        worktrans(iao,iao)= one/sqrt(worktrans(iao,iao))
        do jao= 1,numnmb
          worktrans(jao,iao)=-snao(jao,iao)*worktrans(iao,iao)
        enddo 
      enddo 
!
      call dgemm('N','N',nao,nao,nao,one,pnao,nao,worktrans,nao,zero,work,nao)
      call dgemm('T','N',nao,nao,nao,one,worktrans,nao,work,nao,zero,pnao,nao)
      call dgemm('N','N',nao,nao,nao,one,snao,nao,worktrans,nao,zero,work,nao)
      call dgemm('T','N',nao,nao,nao,one,worktrans,nao,work,nao,zero,snao,nao)
      call dgemm('N','N',nao,nao,nao,one,trans,nao,worktrans,nao,zero,work,nao)
      trans(1:nao,1:nao)= work(1:nao,1:nao)
!
      return
end


!--------------------------------------------------------------------------------
  subroutine orthonrb2(pnao,snao,trans,wnao,work1,work2,work3,numnmb,numlnrb, &
&                      databasisnpa2,datacomp)
!--------------------------------------------------------------------------------
!
! Orthogonalize NRB sets with large w
!
      use modtype, only : typebasis, typecomp
      implicit none
      type(typebasis),intent(in) :: databasisnpa2
      type(typecomp),intent(in) :: datacomp
      integer,intent(in) :: numnmb, numlnrb
      integer :: nao, iao, jao
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: wnao(databasisnpa2%nao)
      real(8),intent(out) :: work1(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(out) :: work2(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(out) :: work3(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(inout) :: pnao(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(inout) :: snao(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(inout) :: trans(databasisnpa2%nao,databasisnpa2%nao)
!
      if(numlnrb == 0) return
!
! Calculate WSW matrix
!
      nao= databasisnpa2%nao
      do iao= 1,numlnrb
        do jao= 1,numlnrb
          work1(jao,iao)= snao(numnmb+jao,numnmb+iao)*wnao(numnmb+iao)*wnao(numnmb+jao)
        enddo
      enddo
!
! Diagonalize WSW matrix and calculate W(WSW)^-1/2
!
      call diag('V','U',numlnrb,work1,nao,work3,datacomp)
      do iao= 1,numlnrb
        do jao= 1,numlnrb
          work2(jao,iao)= work1(jao,iao)/(sqrt(work3(iao,1)))
        enddo
      enddo
      call dgemm('N','T',numlnrb,numlnrb,numlnrb,one,work2,nao,work1,nao,zero,work3,nao)
      do iao= 1,numlnrb
        do jao= 1,numlnrb
          work3(jao,iao)= work3(jao,iao)*wnao(numnmb+jao)
        enddo
      enddo
!
! Occupancy-weighted symmetric orthonalization
!
      call dgemm('N','N',nao,numlnrb,numlnrb,one,pnao(1,numnmb+1),nao,work3,nao,zero,work2,nao)
      pnao(1:nao,numnmb+1:numnmb+numlnrb)= work2(1:nao,1:numlnrb)
      call dgemm('T','N',numlnrb,nao,numlnrb,one,work3,nao,pnao(numnmb+1,1),nao,zero,work2,nao)
      pnao(numnmb+1:numnmb+numlnrb,1:nao)= work2(1:numlnrb,1:nao)
      call dgemm('N','N',nao,numlnrb,numlnrb,one,snao(1,numnmb+1),nao,work3,nao,zero,work2,nao)
      snao(1:nao,numnmb+1:numnmb+numlnrb)= work2(1:nao,1:numlnrb)
      call dgemm('T','N',numlnrb,nao,numlnrb,one,work3,nao,snao(numnmb+1,1),nao,zero,work2,nao)
      snao(numnmb+1:numnmb+numlnrb,1:nao)= work2(1:numlnrb,1:nao)
      call dgemm('N','N',nao,numlnrb,numlnrb,one,trans(1,numnmb+1),nao,work3,nao,zero,work2,nao)
      trans(1:nao,numnmb+1:numnmb+numlnrb)= work2(1:nao,1:numlnrb)
!
      return
end


!----------------------------------------------------------------------------------------------------
  subroutine orthonrb3(pnao,snao,trans,worktrans,work2,work3,work4,numnmb,numnrb,numlnrb,numsnrb, &
&                      databasisnpa2,datacomp)
!----------------------------------------------------------------------------------------------------
!
! Orthogonalize NRB sets with small w
!
      use modtype, only : typebasis, typecomp
      implicit none
      type(typebasis),intent(in) :: databasisnpa2
      type(typecomp),intent(in) :: datacomp
      integer,intent(in) :: numnmb, numnrb, numlnrb, numsnrb
      integer :: nao, iao, jao
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(out) :: worktrans(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(out) :: work2(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(out) :: work3(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(out) :: work4(databasisnpa2%nao)
      real(8),intent(inout) :: pnao(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(inout) :: snao(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(inout) :: trans(databasisnpa2%nao,databasisnpa2%nao)
!
      if(numsnrb == 0) return
!
! Orthogonalize small-w NRB sets to large-w NRB sets by Gram-Schmidt process
!
      nao= databasisnpa2%nao
      worktrans(1:numnrb,1:numsnrb)= zero
!
      do iao= 1,numsnrb
        worktrans(numlnrb+iao,iao)= one
      enddo
      do iao= 1,numsnrb !numnmb+1,nao
        do jao= 1,numlnrb
          worktrans(numlnrb+iao,iao)= worktrans(numlnrb+iao,iao) &
&                                    -snao(numnmb+jao,numnmb+numlnrb+iao)**2
        enddo
        worktrans(numlnrb+iao,iao)= one/sqrt(worktrans(numlnrb+iao,iao))
        do jao= 1,numlnrb
          worktrans(jao,iao)=-snao(numnmb+jao,numnmb+numlnrb+iao)*worktrans(numlnrb+iao,iao)
        enddo
      enddo
!
      call dgemm('N','N',nao,numsnrb,numnrb,one,pnao(1,numnmb+1),nao,worktrans,nao,zero,work2,nao)
      pnao(1:nao,numnmb+numlnrb+1:nao)= work2(1:nao,1:numsnrb)
      call dgemm('T','N',numsnrb,nao,numnrb,one,worktrans,nao,pnao(numnmb+1,1),nao,zero,work2,nao)
      pnao(numnmb+numlnrb+1:nao,1:nao)= work2(1:numsnrb,1:nao)
      call dgemm('N','N',nao,numsnrb,numnrb,one,snao(1,numnmb+1),nao,worktrans,nao,zero,work2,nao)
      snao(1:nao,numnmb+numlnrb+1:nao)= work2(1:nao,1:numsnrb)
      call dgemm('T','N',numsnrb,nao,numnrb,one,worktrans,nao,snao(numnmb+1,1),nao,zero,work2,nao)
      snao(numnmb+numlnrb+1:nao,1:nao)= work2(1:numsnrb,1:nao)
      call dgemm('N','N',nao,numsnrb,numnrb,one,trans(1,numnmb+1),nao,worktrans,nao,zero,work2,nao)
      trans(1:nao,numnmb+numlnrb+1:nao)= work2(1:nao,1:numsnrb)
!
! Orthogonalize within small-w NRB sets
!
      do iao= 1,numsnrb
        do jao= 1,numsnrb
          work2(jao,iao)= snao(numnmb+numlnrb+jao,numnmb+numlnrb+iao)
        enddo
      enddo
!
!   Diagonalize S matrix and calculate S^-1/2
!
      call diag('V','U',numsnrb,work2,nao,work4,datacomp)
      do iao= 1,numsnrb
        do jao= 1,numsnrb
          work3(jao,iao)= work2(jao,iao)/(sqrt(work4(iao)))
        enddo
      enddo
      call dgemm('N','T',numsnrb,numsnrb,numsnrb,one,work3,nao,work2,nao,zero,worktrans,nao)
!
! Occupancy-weighted symmetric orthonalization
!
      call dgemm('N','N',nao,numsnrb,numsnrb,one,pnao(1,numnmb+numlnrb+1),nao,worktrans,nao, &
&                zero,work2,nao)
      pnao(1:nao,numnmb+numlnrb+1:nao)= work2(1:nao,1:numsnrb)
      call dgemm('T','N',numsnrb,nao,numsnrb,one,worktrans,nao,pnao(numnmb+numlnrb+1,1),nao &
&                ,zero,work2,nao)
      pnao(numnmb+numlnrb+1:nao,1:nao)= work2(1:numsnrb,1:nao)
      call dgemm('N','N',nao,numsnrb,numsnrb,one,snao(1,numnmb+numlnrb+1),nao,worktrans,nao &
&                ,zero,work2,nao)
      snao(1:nao,numnmb+numlnrb+1:nao)= work2(1:nao,1:numsnrb)
      call dgemm('T','N',numsnrb,nao,numsnrb,one,worktrans,nao,snao(numnmb+numlnrb+1,1),nao, &
&                zero,work2,nao)
      snao(numnmb+numlnrb+1:nao,1:nao)= work2(1:numsnrb,1:nao)
      call dgemm('N','N',nao,numsnrb,numsnrb,one,trans(1,numnmb+numlnrb+1),nao,worktrans,nao, &
&                zero,work2,nao)
      trans(1:nao,numnmb+numlnrb+1:nao)= work2(1:nao,1:numsnrb)
!
      return
end


!---------------------------------------------------------------------
  subroutine sortbasisnpa2(pnao,snao,trans,work,list2,databasisnpa2)
!---------------------------------------------------------------------
!
! Resort NMB + NRB sets by atoms and angular momenta
!
      use modtype, only : typebasis
      implicit none
      type(typebasis),intent(in) :: databasisnpa2
      integer,intent(in) :: list2(databasisnpa2%nao)
      integer :: iao, jao
      real(8),intent(out) :: work(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(inout) :: pnao(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(inout) :: snao(databasisnpa2%nao,databasisnpa2%nao)
      real(8),intent(inout) :: trans(databasisnpa2%nao,databasisnpa2%nao)
!

      do iao= 1,databasisnpa2%nao
        do jao= 1,databasisnpa2%nao
          work(jao,list2(iao))= pnao(jao,iao)
        enddo
      enddo
      do iao= 1,databasisnpa2%nao
        do jao= 1,databasisnpa2%nao
          pnao(list2(jao),iao)= work(jao,iao)
        enddo
      enddo
      do iao= 1,databasisnpa2%nao
        do jao= 1,databasisnpa2%nao
          work(jao,list2(iao))= snao(jao,iao)
        enddo
      enddo
      do iao= 1,databasisnpa2%nao
        do jao= 1,databasisnpa2%nao
          snao(list2(jao),iao)= work(jao,iao)
        enddo
      enddo
      do iao= 1,databasisnpa2%nao
        do jao= 1,databasisnpa2%nao
          work(jao,list2(iao))= trans(jao,iao)
        enddo
      enddo
      do iao= 1,databasisnpa2%nao
        do jao= 1,databasisnpa2%nao
          trans(jao,iao)= work(jao,iao)
        enddo
      enddo
!
      return
end


!----------------------------------------------------------------------
  subroutine printrnpa(pnao,trans,fock,fnao,work,infonmb,infobasis, &
&                      maxang,maxsize,datamol,databasis,datacomp)
!----------------------------------------------------------------------
!
! Print closed-shell Natural Population Analysis Result
!
      use modtype, only : typemol, typebasis, typecomp
      implicit none
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(in) :: datacomp
      integer,intent(in) :: maxang, maxsize, infonmb(3,0:6,datamol%natom)
      integer,intent(in) :: infobasis(3,0:maxang,datamol%natom)
      integer :: nao, iatom, iang, ishell, iao, jao, numshell, locao, itmp, list(maxsize)
      integer :: ii, jj, imo
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: pnao(databasis%nao,databasis%nao)
      real(8),intent(in) :: trans(databasis%nao,databasis%nao)
      real(8),intent(in) :: fock(databasis%nao*(databasis%nao+1)/2)
      real(8),intent(out) :: fnao(databasis%nao,databasis%nao)
      real(8),intent(out) :: work(databasis%nao,databasis%nao)
      real(8) :: faverage, fdiag(maxsize), pdiag(maxsize), fshell(maxsize), tmp
      real(8) :: chargenpa(3,0:maxang,datamol%natom), atomnpa(0:3,0:datamol%natom), totalcharge
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
      character(len=3) :: anglabel(91)= &
&     (/'s  ','   ','   ','   ','   ','   ','   ','   ','   ','   ','   ','   ','   ', &
&       'px ','py ','pz ','   ','   ','   ','   ','   ','   ','   ','   ','   ','   ', &
&       'd-2','d-1','d0 ','d+1','d+2','   ','   ','   ','   ','   ','   ','   ','   ', &
&       'f-3','f-2','f-1','f0 ','f+1','f+2','f+3','   ','   ','   ','   ','   ','   ', &
&       'g-4','g-3','g-2','g-1','g0 ','g+1','g+2','g+3','g+4','   ','   ','   ','   ', &
&       'h-5','h-4','h-3','h-2','h-1','h0 ','h+1','h+2','h+3','h+4','h+5','   ','   ', &
&       'i-6','i-5','i-4','i-3','i-2','i-1','i0 ','i+1','i+2','i+3','i+4','i+5','i+6'/)
      character(len=3) :: motype(3)=(/'Cor','Val','Ryd'/)
!
      chargenpa(1:3,0:maxang,1:datamol%natom)= zero
      atomnpa(0:3,0:datamol%natom)= zero
      totalcharge= zero
!
      if(datacomp%master) then
        write(datacomp%iout, &
&         '(" -----------------------------------------------------------",/ &
&           "      Natural Population Analysis (Orbitals)",/ &
&           "    NAO   Atom   nlm    Type    Occupancy        Energy",/ &
&           " -----------------------------------------------------------")')
      endif
!
! Transform Fock matrix to Natural Atomic Orbital basis
!
      nao= databasis%nao
      call expand(fock,fnao,nao)
      call dsymm('L','U',nao,nao,one,fnao,nao,trans,nao,zero,work,nao)
      call dgemm('T','N',nao,nao,nao,one,trans,nao,work,nao,zero,fnao,nao)
!
! Summarize Natural Atomic Orbital Energies
!
      locao= 0
      jao  = 0
      do iatom= 1,datamol%natom
        do iang= 0,maxang
          numshell= infobasis(1,iang,iatom)
!
          do ishell= 1,numshell
            faverage= zero
            do iao= 1,2*iang+1
              fdiag((ishell-1)*(2*iang+1)+iao)= fnao(locao+iao,locao+iao)
              pdiag((ishell-1)*(2*iang+1)+iao)= pnao(locao+iao,locao+iao)
              faverage= faverage+fnao(locao+iao,locao+iao)
            enddo
            fshell(ishell)= faverage/dble(2*iang+1)
            list(ishell)= ishell-1
            locao= locao+2*iang+1
          enddo
!
          do ii= 1,numshell-1
            do jj= ii+1, numshell
              if(fshell(ii) > fshell(jj)) then
                tmp= fshell(ii)
                fshell(ii)= fshell(jj)
                fshell(jj)= tmp
                itmp= list(ii)
                list(ii)= list(jj)
                list(jj)= itmp
              endif
            enddo
          enddo
!
          do ishell= 1,numshell
            if(ishell <= infonmb(1,iang,iatom)) then
              if(ishell <= infonmb(2,iang,iatom)) then
                imo= 1
              else
                imo= 2
              endif
            else
              imo= 3
            endif
            do iao= 1,2*iang+1
              jao= jao+1
              write(datacomp%iout,'(i6,i5,1x,a3,1x,i2,a3,4x,a3,f13.6,f16.6)') &
&                   jao, iatom, table(datamol%numatomic(iatom)), &
&                   ishell+infonmb(3,iang,iatom)+iang, anglabel(13*iang+iao), motype(imo),  &
&                   pdiag(list(ishell)*(2*iang+1)+iao), fdiag(list(ishell)*(2*iang+1)+iao)
              chargenpa(imo,iang,iatom)= chargenpa(imo,iang,iatom) &
&                                       +pdiag(list(ishell)*(2*iang+1)+iao)
            enddo
          enddo
        enddo
      enddo
!
      write(datacomp%iout, &
&       '(" -----------------------------------------------------------",/, &
&         " ----------------------------------------------------------------------------",/, &
&         "      Natural Population Analysis (Atoms)",/, &
&         "     Atom       Charge     Pop(Total)     Core        Valence      Rydberg",/, &
&         " ----------------------------------------------------------------------------")')
      do iatom= 1,datamol%natom
        do imo= 1,3
          do iang= 0,maxang
            atomnpa(imo,iatom)= atomnpa(imo,iatom)+chargenpa(imo,iang,iatom)
          enddo
        enddo
        atomnpa(0,iatom)= atomnpa(1,iatom)+atomnpa(2,iatom)+atomnpa(3,iatom)
        atomnpa(0,0)= atomnpa(0,0)+atomnpa(0,iatom)
        atomnpa(1,0)= atomnpa(1,0)+atomnpa(1,iatom)
        atomnpa(2,0)= atomnpa(2,0)+atomnpa(2,iatom)
        atomnpa(3,0)= atomnpa(3,0)+atomnpa(3,iatom)
        totalcharge= totalcharge+datamol%znuc(iatom)-atomnpa(0,iatom)
        write(datacomp%iout,'(i6,1x,a3,5f13.6)') &
&             iatom, table(datamol%numatomic(iatom)), datamol%znuc(iatom)-atomnpa(0,iatom), &
&             atomnpa(0:3,iatom)
      enddo
      write(datacomp%iout, &
&       '(" ----------------------------------------------------------------------------",/, &
&         "    Total ",5f13.6,/, &
&         " ----------------------------------------------------------------------------",/)') &
&         totalcharge, atomnpa(0:3,0)
!
      return
end


!----------------------------------------------------------------------------------------------------------
  subroutine printunpa(focka,fockb,dmtrxa,dmtrxb,overlap,trans,pnao,fnao,worka,workb,infonmb,infobasis, &
&                      maxang,maxsize,datamol,databasis,datacomp)
!----------------------------------------------------------------------------------------------------------
!
! Print open-shell Natural Population Analysis Result
!
      use modtype, only : typemol, typebasis, typecomp
      implicit none
      type(typemol),intent(in) :: datamol
      type(typebasis),intent(in) :: databasis
      type(typecomp),intent(in) :: datacomp
      integer,intent(in) :: maxang, maxsize, infonmb(3,0:6,datamol%natom)
      integer,intent(in) :: infobasis(3,0:maxang,datamol%natom)
      integer :: nao, iatom, iang, ishell, iao, jao, numshell, locao, itmp, list(maxsize)
      integer :: ii, jj, imo
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00
      real(8),intent(in) :: focka(databasis%nao*(databasis%nao+1)/2)
      real(8),intent(in) :: fockb(databasis%nao*(databasis%nao+1)/2)
      real(8),intent(in) :: dmtrxa(databasis%nao*(databasis%nao+1)/2)
      real(8),intent(in) :: dmtrxb(databasis%nao*(databasis%nao+1)/2)
      real(8),intent(in) :: overlap(databasis%nao*(databasis%nao+1)/2)
      real(8),intent(in) :: trans(databasis%nao,databasis%nao)
      real(8),intent(out) :: pnao(databasis%nao,databasis%nao)
      real(8),intent(out) :: fnao(databasis%nao,databasis%nao)
      real(8),intent(out) :: worka(databasis%nao*databasis%nao)
      real(8),intent(out) :: workb(databasis%nao*databasis%nao)
      real(8) :: faverage, fdiag(maxsize), pdiag(maxsize), fshell(maxsize), tmp
      real(8) :: chargenpaa(3,0:maxang,datamol%natom), chargenpab(3,0:maxang,datamol%natom)
      real(8) :: atomnpaa(0:3,0:datamol%natom), atomnpab(0:3,0:datamol%natom), totalchargea, totalchargeb
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
      character(len=3) :: anglabel(91)= &
&     (/'s  ','   ','   ','   ','   ','   ','   ','   ','   ','   ','   ','   ','   ', &
&       'px ','py ','pz ','   ','   ','   ','   ','   ','   ','   ','   ','   ','   ', &
&       'd-2','d-1','d0 ','d+1','d+2','   ','   ','   ','   ','   ','   ','   ','   ', &
&       'f-3','f-2','f-1','f0 ','f+1','f+2','f+3','   ','   ','   ','   ','   ','   ', &
&       'g-4','g-3','g-2','g-1','g0 ','g+1','g+2','g+3','g+4','   ','   ','   ','   ', &
&       'h-5','h-4','h-3','h-2','h-1','h0 ','h+1','h+2','h+3','h+4','h+5','   ','   ', &
&       'i-6','i-5','i-4','i-3','i-2','i-1','i0 ','i+1','i+2','i+3','i+4','i+5','i+6'/)
      character(len=3) :: motype(3)=(/'Cor','Val','Ryd'/)
!
      chargenpaa(1:3,0:maxang,1:datamol%natom)= zero
      chargenpab(1:3,0:maxang,1:datamol%natom)= zero
      atomnpaa(0:3,0:datamol%natom)= zero
      atomnpab(0:3,0:datamol%natom)= zero
      totalchargea= zero
      totalchargeb= zero
!
! Alpha electrons
!
      if(datacomp%master) then
        write(datacomp%iout, &
&         '(" -----------------------------------------------------------",/ &
&           "      Natural Population Analysis (Orbitals, Alpha)",/ &
&           "    NAO   Atom   nlm    Type    Occupancy        Energy",/ &
&           " -----------------------------------------------------------")')
      endif
!
!   Transform alpha Fock matrix to Natural Atomic Orbital basis
!
      nao= databasis%nao
      call expand(dmtrxa,pnao,nao)
      call expand2(overlap,fnao,nao)
      call dsymm('L','U',nao,nao,one,pnao,nao,fnao,nao,zero,worka,nao)
      call dsymm('L','U',nao,nao,one,fnao,nao,worka,nao,zero,pnao,nao)
      call dsymm('L','U',nao,nao,one,pnao,nao,trans,nao,zero,worka,nao)
      call dgemm('T','N',nao,nao,nao,one,trans,nao,worka,nao,zero,pnao,nao)
      call expand(focka,fnao,nao)
      call dsymm('L','U',nao,nao,one,fnao,nao,trans,nao,zero,worka,nao)
      call dgemm('T','N',nao,nao,nao,one,trans,nao,worka,nao,zero,fnao,nao)
!
!   Summarize alpha Natural Atomic Orbital Energies
!
      locao= 0
      jao  = 0
      do iatom= 1,datamol%natom
        do iang= 0,maxang
          numshell= infobasis(1,iang,iatom)
!
          do ishell= 1,numshell
            faverage= zero
            do iao= 1,2*iang+1
              fdiag((ishell-1)*(2*iang+1)+iao)= fnao(locao+iao,locao+iao)
              pdiag((ishell-1)*(2*iang+1)+iao)= pnao(locao+iao,locao+iao)
              faverage= faverage+fnao(locao+iao,locao+iao)
            enddo
            fshell(ishell)= faverage/dble(2*iang+1)
            list(ishell)= ishell-1
            locao= locao+2*iang+1
          enddo
!
          do ii= 1,numshell-1
            do jj= ii+1, numshell
              if(fshell(ii) > fshell(jj)) then
                tmp= fshell(ii)
                fshell(ii)= fshell(jj)
                fshell(jj)= tmp
                itmp= list(ii)
                list(ii)= list(jj)
                list(jj)= itmp
              endif
            enddo
          enddo
!
          do ishell= 1,numshell
            if(ishell <= infonmb(1,iang,iatom)) then
              if(ishell <= infonmb(2,iang,iatom)) then
                imo= 1
              else
                imo= 2
              endif
            else
              imo= 3
            endif
            do iao= 1,2*iang+1
              jao= jao+1
              write(datacomp%iout,'(i6,i5,1x,a3,1x,i2,a3,4x,a3,f13.6,f16.6)') &
&                   jao, iatom, table(datamol%numatomic(iatom)), &
&                   ishell+infonmb(3,iang,iatom)+iang, anglabel(13*iang+iao), motype(imo),  &
&                   pdiag(list(ishell)*(2*iang+1)+iao), fdiag(list(ishell)*(2*iang+1)+iao)
              chargenpaa(imo,iang,iatom)= chargenpaa(imo,iang,iatom) &
&                                       +pdiag(list(ishell)*(2*iang+1)+iao)
              worka(jao)= pdiag(list(ishell)*(2*iang+1)+iao)
            enddo
          enddo
        enddo
      enddo
!
      write(datacomp%iout, &
&       '(" -----------------------------------------------------------",/, &
&         " ----------------------------------------------------------------------------",/, &
&         "      Natural Population Analysis (Atoms, Alpha)",/, &
&         "     Atom       Charge     Pop(Total)     Core        Valence      Rydberg",/, &
&         " ----------------------------------------------------------------------------")')
      do iatom= 1,datamol%natom
        do imo= 1,3
          do iang= 0,maxang
            atomnpaa(imo,iatom)= atomnpaa(imo,iatom)+chargenpaa(imo,iang,iatom)
          enddo
        enddo
        atomnpaa(0,iatom)= atomnpaa(1,iatom)+atomnpaa(2,iatom)+atomnpaa(3,iatom)
        atomnpaa(0,0)= atomnpaa(0,0)+atomnpaa(0,iatom)
        atomnpaa(1,0)= atomnpaa(1,0)+atomnpaa(1,iatom)
        atomnpaa(2,0)= atomnpaa(2,0)+atomnpaa(2,iatom)
        atomnpaa(3,0)= atomnpaa(3,0)+atomnpaa(3,iatom)
        totalchargea= totalchargea+datamol%znuc(iatom)*half-atomnpaa(0,iatom)
        write(datacomp%iout,'(i6,1x,a3,5f13.6)') &
&             iatom, table(datamol%numatomic(iatom)), datamol%znuc(iatom)*half-atomnpaa(0,iatom), &
&             atomnpaa(0:3,iatom)
      enddo
      write(datacomp%iout, &
&       '(" ----------------------------------------------------------------------------",/, &
&         "    Total ",5f13.6,/, &
&         " ----------------------------------------------------------------------------",/)') &
&         totalchargea, atomnpaa(0:3,0)
!
! Beta electrons
!
      if(datacomp%master) then
        write(datacomp%iout, &
&         '(" -----------------------------------------------------------",/ &
&           "      Natural Population Analysis (Orbitals, Beta)",/ &
&           "    NAO   Atom   nlm    Type    Occupancy        Energy",/ &
&           " -----------------------------------------------------------")')
      endif
!
!   Transform beta Fock matrix to Natural Atomic Orbital basis
!
      call expand(dmtrxb,pnao,nao)
      call expand2(overlap,fnao,nao)
      call dsymm('L','U',nao,nao,one,pnao,nao,fnao,nao,zero,workb,nao)
      call dsymm('L','U',nao,nao,one,fnao,nao,workb,nao,zero,pnao,nao)
      call dsymm('L','U',nao,nao,one,pnao,nao,trans,nao,zero,workb,nao)
      call dgemm('T','N',nao,nao,nao,one,trans,nao,workb,nao,zero,pnao,nao)
      call expand(fockb,fnao,nao)
      call dsymm('L','U',nao,nao,one,fnao,nao,trans,nao,zero,workb,nao)
      call dgemm('T','N',nao,nao,nao,one,trans,nao,workb,nao,zero,fnao,nao)
!
!   Summarize beta Natural Atomic Orbital Energies
!
      locao= 0
      jao  = 0
      do iatom= 1,datamol%natom
        do iang= 0,maxang
          numshell= infobasis(1,iang,iatom)
!
          do ishell= 1,numshell
            faverage= zero
            do iao= 1,2*iang+1
              fdiag((ishell-1)*(2*iang+1)+iao)= fnao(locao+iao,locao+iao)
              pdiag((ishell-1)*(2*iang+1)+iao)= pnao(locao+iao,locao+iao)
              faverage= faverage+fnao(locao+iao,locao+iao)
            enddo
            fshell(ishell)= faverage/dble(2*iang+1)
            list(ishell)= ishell-1
            locao= locao+2*iang+1
          enddo
!
          do ii= 1,numshell-1
            do jj= ii+1, numshell
              if(fshell(ii) > fshell(jj)) then
                tmp= fshell(ii)
                fshell(ii)= fshell(jj)
                fshell(jj)= tmp
                itmp= list(ii)
                list(ii)= list(jj)
                list(jj)= itmp
              endif
            enddo
          enddo
!
          do ishell= 1,numshell
            if(ishell <= infonmb(1,iang,iatom)) then
              if(ishell <= infonmb(2,iang,iatom)) then
                imo= 1
              else
                imo= 2
              endif
            else
              imo= 3
            endif
            do iao= 1,2*iang+1
              jao= jao+1
              write(datacomp%iout,'(i6,i5,1x,a3,1x,i2,a3,4x,a3,f13.6,f16.6)') &
&                   jao, iatom, table(datamol%numatomic(iatom)), &
&                   ishell+infonmb(3,iang,iatom)+iang, anglabel(13*iang+iao), motype(imo),  &
&                   pdiag(list(ishell)*(2*iang+1)+iao), fdiag(list(ishell)*(2*iang+1)+iao)
              chargenpab(imo,iang,iatom)= chargenpab(imo,iang,iatom) &
&                                       +pdiag(list(ishell)*(2*iang+1)+iao)
              workb(jao)= pdiag(list(ishell)*(2*iang+1)+iao)
            enddo
          enddo
        enddo
      enddo
!
      write(datacomp%iout, &
&       '(" -----------------------------------------------------------",/, &
&         " ----------------------------------------------------------------------------",/, &
&         "      Natural Population Analysis (Atoms, Beta)",/, &
&         "     Atom       Charge     Pop(Total)     Core        Valence      Rydberg",/, &
&         " ----------------------------------------------------------------------------")')
      do iatom= 1,datamol%natom
        do imo= 1,3
          do iang= 0,maxang
            atomnpab(imo,iatom)= atomnpab(imo,iatom)+chargenpab(imo,iang,iatom)
          enddo
        enddo
        atomnpab(0,iatom)= atomnpab(1,iatom)+atomnpab(2,iatom)+atomnpab(3,iatom)
        atomnpab(0,0)= atomnpab(0,0)+atomnpab(0,iatom)
        atomnpab(1,0)= atomnpab(1,0)+atomnpab(1,iatom)
        atomnpab(2,0)= atomnpab(2,0)+atomnpab(2,iatom)
        atomnpab(3,0)= atomnpab(3,0)+atomnpab(3,iatom)
        totalchargeb= totalchargeb+datamol%znuc(iatom)*half-atomnpab(0,iatom)
        write(datacomp%iout,'(i6,1x,a3,5f13.6)') &
&             iatom, table(datamol%numatomic(iatom)), datamol%znuc(iatom)*half-atomnpab(0,iatom), &
&             atomnpab(0:3,iatom)
      enddo
      write(datacomp%iout, &
&       '(" ----------------------------------------------------------------------------",/, &
&         "    Total ",5f13.6,/, &
&         " ----------------------------------------------------------------------------",/)') &
&         totalchargeb, atomnpab(0:3,0)
!
! Alpha + Beta electrons
!
      if(datacomp%master) then
        write(datacomp%iout, &
&         '(" -----------------------------------------------------------",/ &
&           "      Natural Population Analysis (Orbitals, Alpha+Beta)",/ &
&           "    NAO   Atom   nlm    Type    Occupancy",/ &
&           " -----------------------------------------------------------")')
      endif
!
      jao  = 0
      do iatom= 1,datamol%natom
        do iang= 0,maxang
          numshell= infobasis(1,iang,iatom)
          do ishell= 1,numshell
            if(ishell <= infonmb(1,iang,iatom)) then
              if(ishell <= infonmb(2,iang,iatom)) then
                imo= 1
              else
                imo= 2
              endif
            else
              imo= 3
            endif
            do iao= 1,2*iang+1
              jao= jao+1
              write(datacomp%iout,'(i6,i5,1x,a3,1x,i2,a3,4x,a3,f13.6)') &
&                   jao, iatom, table(datamol%numatomic(iatom)), &
&                   ishell+infonmb(3,iang,iatom)+iang, anglabel(13*iang+iao), motype(imo),  &
&                   worka(jao)+workb(jao)
            enddo
          enddo
        enddo
      enddo
!
      write(datacomp%iout, &
&       '(" -----------------------------------------------------------",/, &
&         " ----------------------------------------------------------------------------",/, &
&         "      Natural Population Analysis (Atoms, Alpha+Beta)",/, &
&         "     Atom       Charge     Pop(Total)     Core        Valence      Rydberg",/, &
&         " ----------------------------------------------------------------------------")')
      do iatom= 1,datamol%natom
        write(datacomp%iout,'(i6,1x,a3,5f13.6)') &
&             iatom, table(datamol%numatomic(iatom)), &
&             datamol%znuc(iatom)-atomnpaa(0,iatom)-atomnpab(0,iatom), &
&             atomnpaa(0:3,iatom)+atomnpab(0:3,iatom)
      enddo
      write(datacomp%iout, &
&       '(" ----------------------------------------------------------------------------",/, &
&         "    Total ",5f13.6,/, &
&         " ----------------------------------------------------------------------------",/)') &
&         totalchargea+totalchargeb, atomnpaa(0:3,0)+atomnpab(0:3,0)
!
!
      return
end
