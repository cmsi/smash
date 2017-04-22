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
!-----------------------
  subroutine nucenergy
!-----------------------
!
! Calculate nuclear replusion energy
!
      use modparallel, only : master
      use modmolecule, only : natom, coord, znuc
      use modthresh, only : threshatom
      use modenergy, only : enuc
      use modwarn, only : nwarn
      implicit none
      integer :: iatom, jatom
      real(8),parameter :: zero=0.0D+00
      real(8) :: xyz(3), rr, chrgij
!
      enuc= zero

      do iatom= 2,natom
        do jatom= 1,iatom-1
          xyz(1)= coord(1,iatom)-coord(1,jatom)
          xyz(2)= coord(2,iatom)-coord(2,jatom)
          xyz(3)= coord(3,iatom)-coord(3,jatom)
          rr= sqrt(xyz(1)*xyz(1)+xyz(2)*xyz(2)+xyz(3)*xyz(3))
          chrgij= znuc(iatom)*znuc(jatom)
          if(rr /= zero) then
            enuc= enuc+chrgij/rr
            if((rr <= threshatom).and.master) then
              write(*,'("Warning! Distance of Atoms",i4," and",i4," is short!")') iatom, jatom
              nwarn= nwarn+1
            endif       
          else
            if((chrgij /= zero).and.master) then
              write(*,'("Error! Atoms",i4," and",i4," are the same position!")') iatom, jatom
              call iabort
            endif
          endif
        enddo
      enddo
      return
end


!---------------------------------------------
  subroutine nucgradient(egrad,nproc,myrank)
!---------------------------------------------
!
! Calculate gradinet of nuclear replusion energy
!
      use modmolecule, only : natom, coord, znuc
      implicit none
      integer,intent(in) :: nproc, myrank
      integer :: iatom, jatom, i
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: xyz(3), rr, chrgij
!
      do iatom= 2+myrank,natom,nproc
        do jatom= 1,iatom-1
          xyz(1)= coord(1,iatom)-coord(1,jatom)
          xyz(2)= coord(2,iatom)-coord(2,jatom)
          xyz(3)= coord(3,iatom)-coord(3,jatom)
          rr= xyz(1)*xyz(1)+xyz(2)*xyz(2)+xyz(3)*xyz(3)
          rr= rr*sqrt(rr)
          chrgij= znuc(iatom)*znuc(jatom)
          do i= 1,3
            egrad(i,iatom)= egrad(i,iatom)-xyz(i)*chrgij/rr
            egrad(i,jatom)= egrad(i,jatom)+xyz(i)*chrgij/rr
          enddo
        enddo
      enddo
      return
end


!-----------------------------------------------------------------------------------
  subroutine setredundantcoord(iredun,isizered,numbond,numangle,numtorsion,exceed)
!-----------------------------------------------------------------------------------
!
! Set redundant internal coordinate
! Covalent radii (H - Cn): P. Pyykko, M. Atsumi, Chem. Eur. J., 186 (2009) 15.
!
      use modparallel, only : master
      use modmolecule, only : coord, natom, numatomic
      use modunit, only : tobohr
      implicit none
      integer,parameter :: maxconnect=13
      integer,intent(in) :: isizered
      integer,intent(out) :: iredun(4,isizered/4), numbond, numangle, numtorsion
      integer :: numredun, maxsize, ijpair(natom,maxconnect), iatom, jatom, katom, icount
      integer :: ibond, kpair, iangle, jangle, numb, iblock(natom), ib, jb, ijatom(2) 
      real(8),parameter :: zero=0.0D+00
      real(8) :: radii(112), rrij, thresh, rrijmin
      logical,intent(out) :: exceed
      data radii/ &
&     0.32D+00, 0.46D+00, 1.33D+00, 1.02D+00, 0.85D+00, 0.75D+00, 0.71D+00, 0.63D+00, 0.64D+00, &
&     0.67D+00, 1.55D+00, 1.39D+00, 1.26D+00, 1.16D+00, 1.11D+00, 1.03D+00, 0.99D+00, 0.96D+00, &
&     1.96D+00, 1.71D+00, 1.48D+00, 1.36D+00, 1.34D+00, 1.22D+00, 1.19D+00, 1.16D+00, 1.11D+00, &
&     1.10D+00, 1.12D+00, 1.18D+00, 1.24D+00, 1.21D+00, 1.21D+00, 1.16D+00, 1.14D+00, 1.17D+00, &
&     2.10D+00, 1.85D+00, 1.63D+00, 1.54D+00, 1.47D+00, 1.38D+00, 1.28D+00, 1.25D+00, 1.25D+00, &
&     1.20D+00, 1.28D+00, 1.36D+00, 1.42D+00, 1.40D+00, 1.40D+00, 1.36D+00, 1.33D+00, 1.31D+00, &
&     2.32D+00, 1.96D+00, 1.80D+00, 1.63D+00, 1.76D+00, 1.74D+00, 1.73D+00, 1.72D+00, 1.68D+00, &
&     1.69D+00, 1.68D+00, 1.67D+00, 1.66D+00, 1.65D+00, 1.64D+00, 1.70D+00, 1.62D+00, 1.52D+00, &
&     1.46D+00, 1.37D+00, 1.31D+00, 1.29D+00, 1.22D+00, 1.23D+00, 1.24D+00, 1.33D+00, 1.44D+00, &
&     1.44D+00, 1.51D+00, 1.45D+00, 1.47D+00, 1.42D+00, 2.23D+00, 2.01D+00, 1.86D+00, 1.75D+00, &
&     1.69D+00, 1.70D+00, 1.71D+00, 1.72D+00, 1.66D+00, 1.66D+00, 1.68D+00, 1.68D+00, 1.65D+00, &
&     1.67D+00, 1.73D+00, 1.76D+00, 1.61D+00, 1.57D+00, 1.49D+00, 1.43D+00, 1.41D+00, 1.34D+00, &
&     1.29D+00, 1.28D+00, 1.21D+00, 1.22D+00/
!
      ijpair(:,:)= 0
      numredun= 0
      maxsize= isizered/4
      exceed=.false.
!
! Set bond strech
!
      numbond= 0
      do iatom= 1,natom
        if(numatomic(iatom) > 112) then
          if(master) &
&         write(*,'(" Error! This program supports up to Cn in Subroutine setredundantcoord.")')
          call iabort
        elseif(numatomic(iatom) < 1) then
          if(master) &
&         write(*,'(" Error! This program does not support dummy and ghost atoms ", &
&                   "in Subroutine setredundantcoord currently.")')
          call iabort
        endif
        icount= 0
        do jatom= 1,natom
          if(jatom == iatom) cycle
          thresh=(1.25D+00*tobohr*(radii(numatomic(iatom))+radii(numatomic(jatom))))**2
          rrij= (coord(1,jatom)-coord(1,iatom))**2+(coord(2,jatom)-coord(2,iatom))**2 &
&              +(coord(3,jatom)-coord(3,iatom))**2
          if(rrij <= thresh) then
            icount= icount+1
            if(icount > maxconnect) then
              if(master) write(*,'(" Error! There are too many atoms near Atom",i4,".")')iatom
              call iabort
            endif
            ijpair(iatom,icount)= jatom
            if(jatom > iatom) then
              numredun= numredun+1
              if(numredun > maxsize) then
                exceed=.true.
                return
              endif
              numbond= numbond+1
              iredun(1,numbond)= iatom
              iredun(2,numbond)= jatom
            endif
          endif
        enddo
      enddo
!
! Check the number of molecular fragments
!
      numb= 0
      iblock(:)= 0
 mblock:do iatom= 1,natom
          if(iblock(iatom) /= 0) cycle
          numb= numb+1
          iblock(iatom)= numb
          call checkbond(iatom,numb,iblock,ijpair,natom,maxconnect)
          do jatom= iatom+1,natom
            if(iblock(jatom) == 0) cycle mblock
          enddo
          exit mblock
        enddo mblock
!
     if(master.and.(numb /= 1)) then
       write(*,'(" There are",i3," moleclar blocks in redundant coordinate.")') numb
     endif
!
! Calculate length between molecular fragments
!
      do ib= 1, numb
        do jb= ib+1,numb
          rrijmin= zero
          do iatom= 1,natom
            if(iblock(iatom) /= ib) cycle
            do jatom= iatom+1,natom
              if(iblock(jatom) /= jb) cycle
                rrij= (coord(1,jatom)-coord(1,iatom))**2+(coord(2,jatom)-coord(2,iatom))**2 &
&                    +(coord(3,jatom)-coord(3,iatom))**2
              if((rrijmin == zero).or.(rrij < rrijmin)) then
                rrijmin= rrij
                ijatom(1)= min(iatom,jatom)
                ijatom(2)= max(iatom,jatom)
              endif
            enddo
          enddo
          do icount= 1,maxconnect
            if(ijpair(ijatom(1),icount) == 0) then
              ijpair(ijatom(1),icount)= ijatom(2)
              exit
            endif
            if(icount == maxconnect) then
              if(master) write(*,'(" Error! There are too many atoms near Atom",i4,".")')ijatom(1)
              call iabort
            endif
          enddo
          do icount= 1,maxconnect
            if(ijpair(ijatom(2),icount) == 0) then
              ijpair(ijatom(2),icount)= ijatom(1)
              exit
            endif
            if(icount == maxconnect) then
              if(master) write(*,'(" Error! There are too many atoms near Atom",i4,".")')ijatom(2)
              call iabort
            endif
          enddo
          numredun= numredun+1
          if(numredun > maxsize) then
            exceed=.true.
            return
          endif
          numbond= numbond+1
          iredun(1,numbond)= ijatom(1)
          iredun(2,numbond)= ijatom(2)
        enddo
      enddo
             
!
! Set bond angle
!
      numangle= 0
      do ibond= 1,numbond
        iatom= iredun(1,ibond)
        jatom= iredun(2,ibond)
        do kpair= 1,maxconnect
          katom= ijpair(jatom,kpair)
          if(katom == 0) exit
          if(katom > iatom) then
            numredun= numredun+1
            if(numredun > maxsize) then
              exceed=.true.
              return
            endif
            numangle= numangle+1
            iredun(1,numbond+numangle)= iatom
            iredun(2,numbond+numangle)= jatom
            iredun(3,numbond+numangle)= katom
          endif
        enddo
        iatom= iredun(2,ibond)
        jatom= iredun(1,ibond)
        do kpair= 1,maxconnect
          katom= ijpair(jatom,kpair)
          if(katom == 0) exit
          if(katom > iatom) then
            numredun= numredun+1
            if(numredun > maxsize) then
              exceed=.true.
              return
            endif
            numangle= numangle+1
            iredun(1,numbond+numangle)= iatom
            iredun(2,numbond+numangle)= jatom
            iredun(3,numbond+numangle)= katom
          endif
        enddo
      enddo
!
!
! Set dihedral angle
!
      numtorsion= 0
      do iangle= numbond+1,numbond+numangle
        iatom= iredun(1,iangle)
        jatom= iredun(2,iangle)
        katom= iredun(3,iangle)
        do jangle= iangle+1,numbond+numangle
          if(iredun(1,jangle) == jatom) then
            if(iredun(2,jangle) == katom) then
              if(iredun(3,jangle) /= iatom) then
                numredun= numredun+1
                if(numredun > maxsize) then
                  exceed=.true.
                  return
                endif
                numtorsion= numtorsion+1
                iredun(1,numbond+numangle+numtorsion)= iatom
                iredun(2,numbond+numangle+numtorsion)= jatom
                iredun(3,numbond+numangle+numtorsion)= katom
                iredun(4,numbond+numangle+numtorsion)= iredun(3,jangle)
              endif
            endif
          endif
          if(iredun(3,jangle) == jatom) then
            if(iredun(2,jangle) == katom) then
              if(iredun(1,jangle) /= iatom) then
                numredun= numredun+1
                if(numredun > maxsize) then
                  exceed=.true.
                  return
                endif
                numtorsion= numtorsion+1
                iredun(1,numbond+numangle+numtorsion)= iatom
                iredun(2,numbond+numangle+numtorsion)= jatom
                iredun(3,numbond+numangle+numtorsion)= katom
                iredun(4,numbond+numangle+numtorsion)= iredun(1,jangle)
              endif
            endif
          endif
        enddo
        iatom= iredun(3,iangle)
        jatom= iredun(2,iangle)
        katom= iredun(1,iangle)
        do jangle= iangle+1,numbond+numangle
          if(iredun(1,jangle) == jatom) then
            if(iredun(2,jangle) == katom) then
              if(iredun(3,jangle) /= iatom) then
                numredun= numredun+1
                if(numredun > maxsize) then
                  exceed=.true.
                  return
                endif
                numtorsion= numtorsion+1
                iredun(1,numbond+numangle+numtorsion)= iatom
                iredun(2,numbond+numangle+numtorsion)= jatom
                iredun(3,numbond+numangle+numtorsion)= katom
                iredun(4,numbond+numangle+numtorsion)= iredun(3,jangle)
              endif
            endif
          endif
          if(iredun(3,jangle) == jatom) then
            if(iredun(2,jangle) == katom) then
              if(iredun(1,jangle) /= iatom) then
                numredun= numredun+1
                if(numredun > maxsize) then
                  exceed=.true.
                  return
                endif
                numtorsion= numtorsion+1
                iredun(1,numbond+numangle+numtorsion)= iatom
                iredun(2,numbond+numangle+numtorsion)= jatom
                iredun(3,numbond+numangle+numtorsion)= katom
                iredun(4,numbond+numangle+numtorsion)= iredun(1,jangle)
              endif
            endif
          endif
        enddo
      enddo
!
      if((natom >= 4).and.(numtorsion == 0)) then
        do iangle= numbond+1,numbond+numangle
          iatom= iredun(1,iangle)
          jatom= iredun(2,iangle)
          katom= iredun(3,iangle)
          do jangle= iangle+1,numbond+numangle
            if(iredun(1,jangle) == katom) then
              if(iredun(2,jangle) == jatom) then
                numredun= numredun+1
                if(numredun > maxsize) then
                  exceed=.true.
                  return
                endif
                numtorsion= numtorsion+1
                iredun(1,numbond+numangle+numtorsion)= iatom
                iredun(2,numbond+numangle+numtorsion)= jatom
                iredun(3,numbond+numangle+numtorsion)= katom
                iredun(4,numbond+numangle+numtorsion)= iredun(3,jangle)
              endif
            endif
          enddo
        enddo
      endif
!
      return
end


!----------------------------------------------------------------------------
  recursive subroutine checkbond(iatom,numb,iblock,ijpair,natom,maxconnect)
!----------------------------------------------------------------------------
!
! Check the number of molecular fragments
!
      implicit none
      integer,intent(in) :: iatom, numb, natom, maxconnect, ijpair(natom,maxconnect)
      integer,intent(inout) :: iblock(natom)
      integer :: icount, jatom
!
      do icount= 1,maxconnect
        jatom= ijpair(iatom,icount)
        if(jatom == 0) exit
        if(iblock(jatom) /= 0) cycle
        iblock(jatom)= numb
        call checkbond(jatom,numb,iblock,ijpair,natom,maxconnect)
      enddo
!
      return
end


!------------------------------------------------------------------------------------
  subroutine calcnewcoord(coord,coordold,egrad,egradold,ehess,displc,natom3,iopt, &
&                         nproc,myrank,mpi_comm)
!------------------------------------------------------------------------------------
!
! Calculate new Cartesian coordinate with gradient and hessian
!
      use modparallel, only : master
      use modprint, only : iprint
      use modunit, only : toang
      use modmolecule, only : numatomic
      implicit none
      integer,intent(in) :: natom3, iopt, nproc, myrank, mpi_comm
      integer :: i, j, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, third=0.3333333333333333D+00
      real(8),intent(inout) :: egrad(natom3), egradold(natom3), ehess(natom3*(natom3+1)/2)
      real(8),intent(inout) :: coord(natom3), coordold(natom3), displc(natom3*3)
!
      real(8) :: work(natom3,natom3),eigen(natom3),work2(natom3,natom3),work3(natom3,natom3)
!
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
!
      if(iopt == 1) then
!
! Set initial hessian
!
        ehess(:)= zero
        do i= 1,natom3
          ehess(i*(i+1)/2)= third
        enddo 
      else
!
! Updata hessian matrix
!
        call hessianbfgs(ehess,coord,coordold,egrad,egradold,displc,natom3)
      endif
!
      do i=1,natom3
        ii= i*(i-1)/2
        do j=1,i
          work(j,i)=ehess(ii+j)
        enddo
      enddo
      call diag('V','U',natom3,work,natom3,eigen,nproc,myrank,mpi_comm)
!
      do i=1,natom3
        eigen(i)= one/eigen(i)
      enddo
      do i=1,natom3
        do j=1,natom3
          work2(j,i)=work(j,i)*eigen(i)
        enddo
      enddo
      call dgemm('N','T',natom3,natom3,natom3,one,work,natom3,work2,natom3,zero,work3,natom3)
!
! Copy old coordinate and gradient
!
      coordold(:)= coord(:)
      egradold(:)= egrad(:)
!
      do i=1,natom3
        do j=1,natom3
          coord(i)=coord(i)-work3(i,j)*egrad(j)
        enddo
      enddo
!
! Print delta xyz
!
      if(master.and.(iprint >= 2)) then
        write(*,'(" ----------------------------------------------------")')
        write(*,'("          Delta xyz (Angstrom)")')
        write(*,'("  Atom            X             Y             Z")')
        write(*,'(" ----------------------------------------------------")')
        do i= 1,natom3/3
          write(*,'(3x,a3,3x,3f14.7)')table(numatomic(i)), &
&              ((coord((i-1)*3+j)-coordold((i-1)*3+j))*toang,j=1,3)
        enddo
        write(*,'(" ----------------------------------------------------")')
      endif
!
      return
end


!---------------------------------------------------------------------------------------
  subroutine calcnewcoordred(coord,coordold,coordredun,egrad,egradredun,ehess,work1, &
&                            work2,work3,work4,workv,iopt,iredun,isizered, &
&                            maxredun,numbond,numangle,numtorsion,numredun, &
&                            nproc,myrank,mpi_comm)
!---------------------------------------------------------------------------------------
!
! Calculate new Cartesian coordinate with gradient, hessian and redundant coordinate
! using Rational Function Optimization (RFO) method
!
! Parameters for force constants (bond strech (up to 3rd period), angle bend, torsion):
!                          H. B. Schlegel, Theoret. Chim. Acta, 333 (1984) 66.
! Parameters for force constants (bond strech (4th - 6th period)):
!                          J. M. Wittbrodt, H. B. Schlegel, J. Mol. Strut., 398 (1997) 55.
! Covalent raddi (H - Kr): H. B. Schlegel, Theoret. Chim. Acta, 333 (1984) 66.
! Covalent radii (Rb- Cn): P. Pyykko, M. Atsumi, Chem. Eur. J., 186 (2009) 15.
!
      use modparallel, only : master
      use modprint, only : iprint
      use modunit, only : toang, tobohr
      use modmolecule, only : numatomic, natom
      implicit none
      integer,parameter :: maxiterdx=100, maxiterrfo=1000
      integer,intent(in) :: iopt, isizered, maxredun, iredun(4,isizered/4)
      integer,intent(in) :: numbond, numangle, numtorsion, numredun, nproc, myrank, mpi_comm
      integer :: irow(112)
      integer :: natom3, ii, jj, ij, kk, iatom, jatom, katom, iterrfo, iterdx
      integer :: numdim
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, convl=1.0D-08, convrms=1.0D-06
      real(8),parameter :: rad2deg=5.729577951308232D+01
      real(8),intent(inout) :: coord(natom*3), coordold(natom*3), coordredun(numredun,2)
      real(8),intent(inout) :: egrad(natom*3), egradredun(numredun,2)
      real(8),intent(inout) :: ehess(numredun*(numredun+1)/2), work1(maxredun,maxredun)
      real(8),intent(inout) :: work2(maxredun,maxredun), work3(maxredun,maxredun)
      real(8),intent(inout) :: work4(maxredun,maxredun), workv(maxredun,3) 
      real(8) :: parambond(7,7), rij, paramb
      real(8) :: rcov(112), rjk, suml, rlambda, rmsdx, rmsqx, ddot
      character(len=33) :: paramred
      character(len=5) :: chartmp(5)
      character(len=3) :: table(112)= &
&     (/'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ','Al ','Si ','P  ',&
&       'S  ','Cl ','Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',&
&       'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ',&
&       'Pd ','Ag ','Cd ','In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',&
&       'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ','Ta ','W  ','Re ',&
&       'Os ','Ir ','Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',&
&       'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ','Md ','No ','Lr ','Rf ','Db ',&
&       'Sg ','Bh ','Hs ','Mt ','Uun','Uuu','Uub'/)
      data irow/2*1, 8*2, 8*3, 18*4, 18*5, 32*6, 26*7/
      data parambond/ &
&     -0.2440D+00, 0.3520D+00, 0.6600D+00, 0.7126D+00, 0.8335D+00, 0.9491D+00, 1.0D+00, &
&      0.3520D+00, 1.0850D+00, 1.5220D+00, 1.4725D+00, 1.6549D+00, 1.7190D+00, 1.8D+00, &
&      0.6600D+00, 1.5220D+00, 2.0680D+00, 1.8238D+00, 2.1164D+00, 2.3185D+00, 2.5D+00, &
&      0.7126D+00, 1.4725D+00, 1.8238D+00, 2.0203D+00, 2.2137D+00, 2.5206D+00, 2.8D+00, &
&      0.8335D+00, 1.6549D+00, 2.1164D+00, 2.2137D+00, 2.3718D+00, 2.5110D+00, 2.7D+00, &
&      0.9491D+00, 1.7190D+00, 2.3185D+00, 2.5206D+00, 2.5110D+00, 2.5D+00,    2.5D+00, &
&      1.0D+00,    1.8D+00,    2.5D+00,    2.8D+00,    2.7D+00,    2.5D+00,    2.5D+00/
      data rcov/ &
&     0.32D+00, 0.60D+00, 1.20D+00, 1.05D+00, 0.81D+00, 0.77D+00, 0.74D+00, 0.74D+00, 0.72D+00, &
&     0.72D+00, 1.50D+00, 1.40D+00, 1.30D+00, 1.17D+00, 1.10D+00, 1.04D+00, 0.99D+00, 0.99D+00, &
&     1.80D+00, 1.60D+00, 1.40D+00, 1.40D+00, 1.40D+00, 1.40D+00, 1.40D+00, 1.40D+00, 1.40D+00, &
&     1.40D+00, 1.40D+00, 1.40D+00, 1.40D+00, 1.30D+00, 1.20D+00, 1.20D+00, 1.10D+00, 1.10D+00, &
&     2.10D+00, 1.85D+00, 1.63D+00, 1.54D+00, 1.47D+00, 1.38D+00, 1.28D+00, 1.25D+00, 1.25D+00, &
&     1.20D+00, 1.28D+00, 1.36D+00, 1.42D+00, 1.40D+00, 1.40D+00, 1.36D+00, 1.33D+00, 1.31D+00, &
&     2.32D+00, 1.96D+00, 1.80D+00, 1.63D+00, 1.76D+00, 1.74D+00, 1.73D+00, 1.72D+00, 1.68D+00, &
&     1.69D+00, 1.68D+00, 1.67D+00, 1.66D+00, 1.65D+00, 1.64D+00, 1.70D+00, 1.62D+00, 1.52D+00, &
&     1.46D+00, 1.37D+00, 1.31D+00, 1.29D+00, 1.22D+00, 1.23D+00, 1.24D+00, 1.33D+00, 1.44D+00, &
&     1.44D+00, 1.51D+00, 1.45D+00, 1.47D+00, 1.42D+00, 2.23D+00, 2.01D+00, 1.86D+00, 1.75D+00, &
&     1.69D+00, 1.70D+00, 1.71D+00, 1.72D+00, 1.66D+00, 1.66D+00, 1.68D+00, 1.68D+00, 1.65D+00, &
&     1.67D+00, 1.73D+00, 1.76D+00, 1.61D+00, 1.57D+00, 1.49D+00, 1.43D+00, 1.41D+00, 1.34D+00, &
&     1.29D+00, 1.28D+00, 1.21D+00, 1.22D+00/
!
      natom3= natom*3
!
! Calculate B-matrix
!
      call calcbmatrix(coordredun,coord,work1,iredun,numbond,numangle,numtorsion)
!
! Calculate G=B*Bt
!
      call dgemm('N','T',numredun,numredun,natom3,one,work1,numredun,work1,numredun, &
&                zero,work2,maxredun)
!
! Calculate G-inverse and K-matrix
!
      call diag('V','U',numredun,work2,maxredun,workv(1,3),nproc,myrank,mpi_comm)
      numdim=0
      do ii= 1,numredun
        if(abs(workv(ii,3)) >= 1.0D-5) then
          numdim= numdim+1
          workv(ii,3)= one/workv(ii,3)
        endif
      enddo
!$OMP parallel do
      do ii= numredun-numdim+1,numredun
        do jj= 1,numredun
          work3(jj,ii)= work2(jj,ii)*workv(ii,3)
        enddo
      enddo
!$OMP end parallel do
      call dgemm('N','T',numredun,numredun,numdim,one,work2(1,numredun-numdim+1),maxredun, &
&                work3(1,numredun-numdim+1),maxredun,zero,work4,maxredun)
!
! Calculate (G-inverse)*B
!
      call dgemm('N','N',numredun,natom3,numredun,one,work4,maxredun,work1,numredun, &
&                zero,work3,maxredun)
!
! Calculate Fred=(G-inverse)*B*Fcart
!
      call dgemv('N',numredun,natom3,one,work3,maxredun,egrad,1,zero,egradredun,1)
!
      if(iopt >= 2) then
!
! Update Hessian
!
        call hessianbfgsred(ehess,coordredun,coordredun(1,2),egradredun,egradredun(1,2), &
&                           workv,numredun,numbond,numangle,numtorsion)
!
        do ii= 1,numredun
          ij= ii*(ii-1)/2
          do jj= 1,ii
            work4(jj,ii)=ehess(ij+jj)
          enddo
        enddo
      else
!
! Set initial Hessian
!
        do ii= 1,numbond
          iatom= numatomic(iredun(1,ii))
          jatom= numatomic(iredun(2,ii))
          rij= coordredun(ii,1)
          paramb= parambond(irow(iatom),irow(jatom))
          work4(ii,ii)= 1.734D+00/((rij-paramb)*(rij-paramb)*(rij-paramb))
        enddo
        do ii= numbond+1,numbond+numangle
          iatom= numatomic(iredun(1,ii))
          katom= numatomic(iredun(3,ii))
          if((iatom == 1).or.(katom == 1)) then
            work4(ii,ii)= 0.16D+00
          else
            work4(ii,ii)= 0.25D+00
          endif
        enddo
        do ii= numbond+numangle+1,numredun
          jatom= numatomic(iredun(2,ii))
          katom= numatomic(iredun(3,ii))
          jj=(iredun(2,ii)-1)*3
          kk=(iredun(3,ii)-1)*3
          rjk=sqrt((coord(jj+1)-coord(kk+1))**2+(coord(jj+2)-coord(kk+2))**2 &
&                 +(coord(jj+3)-coord(kk+3))**2)-(rcov(jatom)+rcov(katom))*tobohr
          if(rjk > zero) rjk= zero
          work4(ii,ii)= 0.0023D+00-0.07D+00*rjk
        enddo
!
        do ii= 1,numredun
          ij= ii*(ii-1)/2
          do jj= 1,ii-1
            work4(jj,ii)= zero
            ehess(ij+jj)= zero
          enddo
          ehess(ij+ii)= work4(ii,ii)
        enddo
      endif
!
! Copy coordinates and gradients
!
      do ii= 1,natom3
        coordold(ii)= coord(ii)
      enddo
      do ii= 1,numredun
        coordredun(ii,2)= coordredun(ii,1)
        egradredun(ii,2)= egradredun(ii,1)
      enddo
!
! Calculate Kt*FC*K and diagonalize, and then obtain eigenvalues and orthogonal basis
!
      call dsymm('L','U',numredun,numdim,one,work4,maxredun,work2(1,numredun-numdim+1),maxredun, &
&                zero,work1,maxredun)
      call dgemm('T','N',numdim,numdim,numredun,one,work2(1,numredun-numdim+1),maxredun, &
&                work1,maxredun,zero,work4,maxredun)
      call diag('V','U',numdim,work4,maxredun,workv(1,3),nproc,myrank,mpi_comm)
!
! Calculate K*(orthogonal basis)
!
      call dgemm('N','N',numredun,numdim,numdim,one,work2(1,numredun-numdim+1),maxredun, &
&                work4,maxredun,zero,work1,maxredun)
!
! Calculate (K*(orthogonal basis))t*Fred
!
      call dgemv('T',numredun,numdim,one,work1,maxredun,egradredun,1,zero,workv,1)
!
! Rational Function Optimization (RFO) step
!
      rlambda= zero
      do iterrfo= 1,maxiterrfo
        suml= zero
        do ii= 1,numdim
          suml= suml+workv(ii,1)*workv(ii,1)/(rlambda-workv(ii,3))
        enddo
        if(master.and.(iprint >= 3)) then
          write(*,'(" Lambda iteration of RFO",i3,3x,"Lambda=",1p,d15.8,4x,"Sum=",1p,d15.8)') &
&               iterrfo,rlambda,suml
        endif
        if(abs(rlambda-suml) <= convl) exit
        rlambda= suml
        if(iterrfo == maxiterrfo) then
          if(master) write(*,'(" Error! RFO step in calcnewcoordred dit not converge.")')
          call iabort
        endif
      enddo
!
! Calculate displacement work(*,3) in coordinate where Hessian is diagonal
!
      do ii=1,numdim
        workv(ii,1)= workv(ii,1)/(rlambda-workv(ii,3))
      enddo
!
! Calculate displacement work(*,2) in redundant coordinate
!
      call dgemv('N',numredun,numdim,one,work1,maxredun,workv,1,zero,workv(1,2),1)
!
! Calculate displacement work(*,1) in Cartesian coordinate
!
      call dgemv('T',numredun,natom3,one,work3,maxredun,workv(1,2),1,zero,workv,1)
!
! Calculate next Cartesian coordinate
!
      do iterdx= 1,maxiterdx
!
! Update Catesian coordinate
!
        do ii= 1,natom3
          coord(ii)=coord(ii)+workv(ii,1)
        enddo
!
! Check convergence of displacement in Cartesian coordinate
!
        rmsdx= sqrt(ddot(natom3,workv,1,workv,1)/natom3)
        rmsqx= sqrt(ddot(numredun,workv(1,2),1,workv(1,2),1)/numredun)
        if(master.and.(iprint >= 3)) then
          write(*,'(" Displacement Iteration",i3,2x,"RMS(Cart)=",1p,d10.3,4x, &
&                   "RMS(Red)=",1p,d10.3)') iterdx, rmsdx, rmsqx
        endif
        if(rmsdx < convrms) exit
!
! Calculate B-matrix and new displacement in redundant coordinate
!
! delta-q workv(*,3)
!
        call calcbmatrix(workv,coord,work1,iredun,numbond,numangle,numtorsion)
        do ii= 1,numredun
          workv(ii,3)= workv(ii,1)-coordredun(ii,1)
        enddo
        call fixdtor(workv(1,3),numbond,numangle,numtorsion)
!
! delta(delta-q) workv(*,2)
!
        do ii= 1,numredun
          workv(ii,2)=workv(ii,2)-workv(ii,3)
        enddo
        call fixdtor(workv(1,2),numbond,numangle,numtorsion)
!
        do ii= 1,numredun
          coordredun(ii,1)= workv(ii,1)
        enddo
!
! Calculate G=B*Bt
!
        call dgemm('N','T',numredun,numredun,natom3,one,work1,numredun,work1,numredun, &
&                  zero,work2,maxredun)
!
!   Calculate G-inverse
!
        call diag('V','U',numredun,work2,maxredun,workv(1,1),nproc,myrank,mpi_comm)
        numdim=0
        do ii= 1,numredun
          if(abs(workv(ii,1)) >= 1.0D-5) then
            numdim= numdim+1
            workv(ii,1)= one/workv(ii,1)
          endif
        enddo
        do ii= numredun-numdim+1,numredun
          do jj= 1,numredun
            work3(jj,ii)= work2(jj,ii)*workv(ii,1)
          enddo
        enddo
        call dgemm('N','T',numredun,numredun,numdim,one,work2(1,numredun-numdim+1),maxredun, &
&                  work3(1,numredun-numdim+1),maxredun,zero,work4,maxredun)
!
! Calculate Bt*(G-inverse)
!
        call dgemm('T','N',natom3,numredun,numredun,one,work1,numredun,work4,maxredun, &
&                  zero,work3,maxredun)
!
! Calculate displacement work(*,1) in Cartesian coordinate
!
        call dgemv('N',natom3,numredun,one,work3,maxredun,workv(1,2),1,zero,workv,1)
!
        if(iterdx == maxiterdx) then
          if(master) then
            write(*,'(" Error! Transformation from redundant to Cartesian did not converge.")')
          endif
          call iabort
        endif
      enddo
      if(master) then
        if(iprint >= 3) write(*,*)
        write(*,'(" ---------------------------------------------------------------")')
        write(*,'("   Redundant coordinate parameters (Angstrom and Degree)")')
        write(*,'("                                        New           Old")')
        write(*,'(" ---------------------------------------------------------------")')
        do ii= 1,numbond
          write(chartmp(1:3),'(i5)')ii,iredun(1:2,ii)
          paramred= trim(trim("Bond"//adjustl(chartmp(1)) //"   ("//adjustl(chartmp(2)))//"," &
&                                  //adjustl(chartmp(3)))//")"
          write(*,'(3x,a33,f9.4,5x,f9.4)')paramred,coordredun(ii,1)*toang, &
&                                                  coordredun(ii,2)*toang
        enddo
        do ii= numbond+1,numbond+numangle
          write(chartmp(1:4),'(i5)')ii-numbond,iredun(1:3,ii)
          paramred= trim(trim(trim("Angle"//adjustl(chartmp(1)) //"  ("//adjustl(chartmp(2))) &
&                                    //","//adjustl(chartmp(3)))//","//adjustl(chartmp(4)))//")"
          write(*,'(3x,a33,f9.4,5x,f9.4)')paramred,coordredun(ii,1)*rad2deg, &
&                                                  coordredun(ii,2)*rad2deg
        enddo
        do ii= numbond+numangle+1,numbond+numangle+numtorsion
          write(chartmp(1:5),'(i5)')ii-numbond-numangle,iredun(1:4,ii)
          paramred= trim(trim(trim(trim("Torsion"//adjustl(chartmp(1)) //"("// &
&                   adjustl(chartmp(2)))//","//adjustl(chartmp(3)))//","// &
&                   adjustl(chartmp(4)))//","//adjustl(chartmp(5)))//")"
          write(*,'(3x,a33,f9.4,5x,f9.4)')paramred,coordredun(ii,1)*rad2deg, &
&                                                  coordredun(ii,2)*rad2deg
        enddo
        write(*,'(" ---------------------------------------------------------------")')
      endif

!
! Print delta xyz
!
      if(master.and.(iprint >= 2)) then
        write(*,'(" ----------------------------------------------------")')
        write(*,'("          Delta xyz (Angstrom)")')
        write(*,'("  Atom            X             Y             Z")')
        write(*,'(" ----------------------------------------------------")')
        do ii= 1,natom3/3
          write(*,'(3x,a3,3x,3f14.7)')table(numatomic(ii)), &
&              ((coord((ii-1)*3+jj)-coordold((ii-1)*3+jj))*toang,jj=1,3)
        enddo
        write(*,'(" ----------------------------------------------------")')
      endif
!
      return
end


!-----------------------------------------------------------------------------------
  subroutine calcbmatrix(coordredun,coord,bmat,iredun,numbond,numangle,numtorsion)
!-----------------------------------------------------------------------------------
!
! Calculate transformation matrix(B-matrix) from Cartesian to internal coordinate
!
      use modmolecule, only : natom
      implicit none
      integer,intent(in) :: numbond, numangle, numtorsion, iredun(4,numbond+numangle+numtorsion)
      integer :: iatom, jatom, katom, latom, ii, jj, kk, ll, mm, ibond, iangle, itorsion
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00, four=4.0D+00
      real(8),intent(in) :: coord(3,natom)
      real(8),intent(out) :: coordredun(numbond+numangle+numtorsion)
      real(8),intent(out) :: bmat(numbond+numangle+numtorsion,natom*3)
      real(8) :: rrij, rrji, rrjk, rji, rjk, rrkl, rij, rkl, pi, pi2
      real(8) :: dotj, dotk, dotjk, dotjk3, cp1(3), cp2(3), cp3(3), sinj, sink, b1, b2, b3
      real(8) :: xyzji(3), xyzjk(3), xyzij(3), xyzkl(3)
      real(8) :: unitji(3), unitjk(3), unitij(3), unitkl(3)
!
      pi= atan(one)*four
      pi2=pi*two
      bmat(:,:)= zero
!
! Bond strech
!
      do ibond= 1,numbond
        iatom= iredun(1,ibond)
        jatom= iredun(2,ibond)
        rrij= (coord(1,jatom)-coord(1,iatom))**2+(coord(2,jatom)-coord(2,iatom))**2 &
&            +(coord(3,jatom)-coord(3,iatom))**2
        rij= sqrt(rrij)
        coordredun(ibond)= rij
!
        ii=(iatom-1)*3
        jj=(jatom-1)*3
        do mm= 1,3
          bmat(ibond,ii+mm)=-(coord(mm,jatom)-coord(mm,iatom))/rij
          bmat(ibond,jj+mm)= (coord(mm,jatom)-coord(mm,iatom))/rij
        enddo
      enddo
!
! Bond angle
!
      do iangle= numbond+1,numbond+numangle
        iatom= iredun(1,iangle)
        jatom= iredun(2,iangle)
        katom= iredun(3,iangle)
        rrji= zero
        rrjk= zero
        do mm= 1,3
          xyzji(mm)= coord(mm,iatom)-coord(mm,jatom)
          xyzjk(mm)= coord(mm,katom)-coord(mm,jatom)
          rrji= rrji+xyzji(mm)*xyzji(mm)
          rrjk= rrjk+xyzjk(mm)*xyzjk(mm)
        enddo
        rji= sqrt(rrji)
        rjk= sqrt(rrjk)
        dotj= zero
        do mm= 1,3
          unitji(mm)= xyzji(mm)/rji
          unitjk(mm)= xyzjk(mm)/rjk
          dotj= dotj+unitji(mm)*unitjk(mm)
        enddo
        if(abs(dotj) >= one) then
          write(*,'(" Error! During calculation of bond angles in calcbmatrix.")')
          write(*,'(" Use Cartesian coordinate. The input is")')
          write(*,'("   opt cartesian=.true.",/)')
          call iabort
        endif
        coordredun(iangle)= acos(dotj)
!
        sinj= sqrt(one-dotj*dotj)
        ii=(iatom-1)*3
        jj=(jatom-1)*3
        kk=(katom-1)*3
        do mm= 1,3
          b1=(dotj*unitji(mm)-unitjk(mm))/(rji*sinj)
          b2=(dotj*unitjk(mm)-unitji(mm))/(rjk*sinj)
          bmat(iangle,ii+mm)= b1
          bmat(iangle,kk+mm)= b2
          bmat(iangle,jj+mm)=-(b1+b2)
        enddo
      enddo
!
! Set dihedral angle
!
      do itorsion= numbond+numangle+1,numbond+numangle+numtorsion
        iatom= iredun(1,itorsion)
        jatom= iredun(2,itorsion)
        katom= iredun(3,itorsion)
        latom= iredun(4,itorsion)
        rrij= zero
        rrjk= zero
        rrkl= zero
        do mm= 1,3
          xyzij(mm)= coord(mm,jatom)-coord(mm,iatom)
          xyzjk(mm)= coord(mm,katom)-coord(mm,jatom)
          xyzkl(mm)= coord(mm,latom)-coord(mm,katom)
          rrij= rrij+xyzij(mm)*xyzij(mm)
          rrjk= rrjk+xyzjk(mm)*xyzjk(mm)
          rrkl= rrkl+xyzkl(mm)*xyzkl(mm)
        enddo
        rij= sqrt(rrij)
        rjk= sqrt(rrjk)
        rkl= sqrt(rrkl)
        dotj= zero
        dotk= zero
        do mm= 1,3
          unitij(mm)= xyzij(mm)/rij
          unitjk(mm)= xyzjk(mm)/rjk
          unitkl(mm)= xyzkl(mm)/rkl
          dotj= dotj-unitij(mm)*unitjk(mm)
          dotk= dotk-unitjk(mm)*unitkl(mm)
        enddo
        cp1(1)= unitij(2)*unitjk(3)-unitij(3)*unitjk(2)
        cp1(2)= unitij(3)*unitjk(1)-unitij(1)*unitjk(3)
        cp1(3)= unitij(1)*unitjk(2)-unitij(2)*unitjk(1)
        cp2(1)= unitjk(2)*unitkl(3)-unitjk(3)*unitkl(2)
        cp2(2)= unitjk(3)*unitkl(1)-unitjk(1)*unitkl(3)
        cp2(3)= unitjk(1)*unitkl(2)-unitjk(2)*unitkl(1)
        cp3(1)= cp1(2)*cp2(3)-cp1(3)*cp2(2)
        cp3(2)= cp1(3)*cp2(1)-cp1(1)*cp2(3)
        cp3(3)= cp1(1)*cp2(2)-cp1(2)*cp2(1)
        if((abs(dotj) >= one).or.(abs(dotk) >= one)) then
          write(*,'(" Error! During calculation of torsion angles in calcbmatrix.")')
          call iabort
        endif
!
        sinj= sqrt(one-dotj*dotj)
        sink= sqrt(one-dotk*dotk)
        dotjk= zero
        dotjk3=zero
        do mm= 1,3
          dotjk= dotjk+cp1(mm)*cp2(mm)/(sinj*sink)
          dotjk3=dotjk3+unitjk(mm)*cp3(mm)
        enddo
        if(abs(dotjk) > one) dotjk= sign(one,dotjk)
        coordredun(itorsion)= acos(dotjk)
        if(dotjk3 < zero) coordredun(itorsion)=-coordredun(itorsion)
!
        if(coordredun(itorsion) > pi) coordredun(itorsion)= coordredun(itorsion)-pi2
        if(coordredun(itorsion) <-pi) coordredun(itorsion)= coordredun(itorsion)+pi2
!
        ii=(iatom-1)*3
        jj=(jatom-1)*3
        kk=(katom-1)*3
        ll=(latom-1)*3
        do mm= 1,3
          b1=-cp1(mm)/(rij*sinj*sinj)
          b2= cp1(mm)*(rjk-rij*dotj)/(rjk*rij*sinj*sinj)-dotk*cp2(mm)/(rjk*sink*sink)
          b3= cp2(mm)/(rkl*sink*sink)
          bmat(itorsion,ii+mm)= b1
          bmat(itorsion,jj+mm)= b2
          bmat(itorsion,ll+mm)= b3
          bmat(itorsion,kk+mm)=-(b1+b2+b3)
        enddo
      enddo
!
      return
end


!--------------------------------------------------------------
  subroutine fixdtor(dcoordredun,numbond,numangle,numtorsion)
!--------------------------------------------------------------
!
! Fix dihedral displacements
!
      implicit none
      integer,intent(in) :: numbond, numangle, numtorsion
      integer :: ii
      real(8),parameter :: pi2=6.283185307179586D+00
      real(8),intent(inout) :: dcoordredun(numbond+numangle+numtorsion)
!
!$OMP parallel do
      do ii= numbond+numangle+1,numbond+numangle+numtorsion
        if(abs(dcoordredun(ii)+pi2) < abs(dcoordredun(ii))) &
&         dcoordredun(ii)= dcoordredun(ii)+pi2
        if(abs(dcoordredun(ii)+pi2) < abs(dcoordredun(ii))) &
&         dcoordredun(ii)= dcoordredun(ii)+pi2
        if(abs(dcoordredun(ii)-pi2) < abs(dcoordredun(ii))) &
&         dcoordredun(ii)= dcoordredun(ii)-pi2
        if(abs(dcoordredun(ii)-pi2) < abs(dcoordredun(ii))) &
&         dcoordredun(ii)= dcoordredun(ii)-pi2
      enddo
!$OMP end parallel do
!
      return
end







