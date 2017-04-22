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
!---------------------------------------------------------------------------------------------
  subroutine formrfockexcor(fockdsum,fockd,energy,totalelec,cmo,atomvec,radpt,angpt, &
&                           rad,ptweight,vao,vmo,xyzpt,rsqrd,transcmo,work,idftex,idftcor, &
&                           nproc,myrank,mpi_comm)
!---------------------------------------------------------------------------------------------
!
! Driver of DFT Fock matrix formation from exchange-correlation functionals
!
! In  : cmo     (MO coefficients)
!       atomvec (atom vector and distance)
!       radpt   (radial point)
!       angpt   (angular point)
!       rad     (atom radius)
!       ptweight(weight of grid point)
! Out : fockdsum (Fock matrix of exchange-correlation terms)
!       fockd    (work space)
!       energy   (Exchange-correlation energy)
!       totalelec(Number of numerially integrated electrons)
!       vao,vmo,xyzpt,work (work space)
!
      use modmolecule, only : natom, neleca
      use moddft, only : nrad, nleb
      use modbasis, only : nao
      use modthresh, only : threshweight, threshrho, threshdfock, threshdftao
      implicit none
      integer,intent(in) :: idftex, idftcor, nproc, myrank, mpi_comm
      integer :: ngridatom, iatom, irad, ileb, icount, ilebstart, jatom, imo
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00
      real(8),intent(in) :: cmo(nao,neleca), atomvec(5,natom,natom), radpt(2,nrad)
      real(8),intent(in) :: angpt(4,nleb), rad(natom), ptweight(nleb,nrad,natom)
      real(8),intent(out) :: fockdsum(nao*(nao+1)/2), fockd(nao*(nao+1)/2), energy, totalelec
      real(8),intent(out) :: vao(nao,4), vmo(neleca,4), xyzpt(3,natom), rsqrd(natom)
      real(8),intent(out) :: transcmo(neleca,nao), work(nao)
      real(8) :: weight, rhoa, grhoa(3), excora(4)
      real(8) :: radpoint, tmp(2,2), wcutoff, rcutoff, fcutoff, aocutoff
!
      transcmo=transpose(cmo)
      fockd(:)= zero
!
      energy= zero
      totalelec= zero
      ngridatom= nrad*nleb
!
      wcutoff= threshweight/(natom*ngridatom)
      rcutoff= threshrho/(natom*ngridatom)
      fcutoff= threshdfock/(natom*ngridatom)
      aocutoff=threshdftao/(natom*ngridatom)
!
!$OMP parallel do collapse(2) schedule(dynamic,1) private(icount,ilebstart,xyzpt, &
!$OMP rsqrd,weight,rhoa,grhoa,excora,radpoint,vao,vmo,work) &
!$OMP reduction(+:energy,fockd,totalelec)
      do iatom= 1,natom
        do irad= 1,nrad
          icount=(iatom-1)*ngridatom+(irad-1)*nleb+1+myrank
          ilebstart=mod(icount,nproc)+1
          radpoint= rad(iatom)*radpt(1,irad)
          do ileb= ilebstart,nleb,nproc
!
            weight=ptweight(ileb,irad,iatom)
            if(weight < wcutoff) cycle
!
            do jatom= 1,natom
              xyzpt(1,jatom)= atomvec(1,iatom,jatom)+radpoint*angpt(1,ileb)
              xyzpt(2,jatom)= atomvec(2,iatom,jatom)+radpoint*angpt(2,ileb)
              xyzpt(3,jatom)= atomvec(3,iatom,jatom)+radpoint*angpt(3,ileb)
              rsqrd(jatom)= xyzpt(1,jatom)*xyzpt(1,jatom)+xyzpt(2,jatom)*xyzpt(2,jatom) &
&                          +xyzpt(3,jatom)*xyzpt(3,jatom)
            enddo
!
            call gridraomo(vao,vmo,transcmo,xyzpt,rsqrd,aocutoff)
!
            rhoa= zero
            grhoa(1:3)= zero
            do imo= 1,neleca
              rhoa=     rhoa    +vmo(imo,1)*vmo(imo,1)
              grhoa(1)= grhoa(1)+vmo(imo,2)*vmo(imo,1)
              grhoa(2)= grhoa(2)+vmo(imo,3)*vmo(imo,1)
              grhoa(3)= grhoa(3)+vmo(imo,4)*vmo(imo,1)
            enddo
            if((rhoa*two) < rcutoff)cycle
            grhoa(1)= grhoa(1)*two
            grhoa(2)= grhoa(2)*two
            grhoa(3)= grhoa(3)*two
            call calcexcor(excora,excora,energy,rhoa,rhoa,grhoa,grhoa,weight,idftex,idftcor,1)
            call fockexcor(fockd,excora,vao,vao(1,2),work,weight,nao,fcutoff)
            totalelec=totalelec+weight*rhoa*two
          enddo
        enddo
      enddo
!$OMP end parallel do
!
      tmp(1,1)= energy
      tmp(2,1)= totalelec
      call para_allreducer(fockd,fockdsum,nao*(nao+1)/2,mpi_comm)
      call para_allreducer(tmp(1,1),tmp(1,2),2,mpi_comm)
      energy    = tmp(1,2)
      totalelec = tmp(2,2)
      return
end


!----------------------------------------------------------------------------------------------
  subroutine calcexcor(excora,excorb,energy,rhoa,rhob,grhoa,grhob,weight,idftex,idftcor,iscf)
!----------------------------------------------------------------------------------------------
!
! Calculate exchange and correlation energy at a grid point
!
! In  : idft (1:B3LYP)
!       iscf (1:RHF, 2:UHF)
!
      implicit none
      integer,intent(in) :: idftex, idftcor, iscf
      real(8),parameter :: zero=0.0D+00, onethird=0.3333333333333333D+00, two=2.0D+00
      real(8),parameter :: four=4.0D+00
      real(8),intent(in) :: rhoa, rhob, grhoa(3), grhob(3), weight
      real(8),intent(out) :: excora(4), excorb(4)
      real(8),intent(inout) :: energy
      real(8) :: rhoa13, rhob13, csdlda, cb88, cvwn, clyp !, cpw91lda, cpw91
!     real(8) :: vrhoa, vrhob, zk, vsigmaaa, vsigmaab, vsigmabb, gradaa, gradab, gradbb
!
      excora(1:4)= zero
      excorb(1:4)= zero
      rhoa13= rhoa**onethird
      rhob13= rhob**onethird
!
! B3LYP
      select case(idftex)
        case(1)
          csdlda= 0.08D+00
          cb88=   0.72D+00
          call funcsdlda(excora,excorb,energy,rhoa,rhob,rhoa13,rhob13,weight,csdlda,iscf)
          call funcbecke88(excora,excorb,energy,rhoa,rhob,grhoa,grhob,rhoa13,rhob13, &
&                        weight,cb88,iscf)
      end select

      select case(idftcor)
        case(1)
          cvwn=   0.19D+00
          clyp=   0.81D+00
          call funcvwn3(excora,excorb,energy,rhoa,rhob,weight,cvwn,iscf)
          call funclyp(excora,excorb,energy,rhoa,rhob,grhoa,grhob,rhoa13,rhob13,weight,clyp,iscf)
        case(2)
          cvwn=   0.19D+00
          clyp=   0.81D+00
          call funcvwn5(excora,excorb,energy,rhoa,rhob,weight,cvwn,iscf)
          call funclyp(excora,excorb,energy,rhoa,rhob,grhoa,grhob,rhoa13,rhob13,weight,clyp,iscf)
      end select
      return
end


!----------------------------------------------------------------------
  subroutine fockexcor(fockd,excora,vao,vgao,work,weight,nao,fcutoff)
!----------------------------------------------------------------------
!
! Add exchange-correlation term to Fock matrix at a grid point
!
      implicit none
      integer,intent(in) :: nao
      integer :: ij, i, j
      real(8),parameter :: zero=0.0D+00, half=0.5D+00
      real(8),intent(in) :: excora(4), vao(nao), vgao(nao,3), weight, fcutoff
      real(8),intent(inout) :: fockd(nao*(nao+1)/2), work(nao)
      real(8) :: vaomax, gridmax, vaoabs, gridabs
!
      ij= 0
      vaomax= zero
      gridmax= zero
      do i= 1,nao
        work(i)=(excora(1)*vao(i)*half+excora(2)*vgao(i,1) &
&               +excora(3)*vgao(i,2)  +excora(4)*vgao(i,3))*weight
        vaoabs= abs(vao(i))
        gridabs= abs(work(i))
        vaomax= max(vaomax,vaoabs)
        gridmax= max(gridmax,gridabs)
        if((vaoabs*gridmax+gridabs*vaomax) >= fcutoff) then
          do j= 1,i
            ij= ij+1
            fockd(ij)= fockd(ij)+vao(i)*work(j)+work(i)*vao(j)
          enddo
        else
          ij= ij+i
        endif
      enddo
      return
end


!-------------------------------------------------------------------------------------------
  subroutine formufockexcor(fockd1,fockd2,fockd3,energy,totalelec,cmoa,cmob,atomvec, &
&                           radpt,angpt,rad,ptweight,vao,vmoa,vmob,xyzpt,rsqrd, &
&                           transcmoa,transcmob,work,idftex,idftcor,nproc,myrank,mpi_comm)
!-------------------------------------------------------------------------------------------
!
! Driver of unrestricted DFT Fock matrix formation from exchange-correlation functionals
!
! In  : cmoa  (Alpha MO coefficient)
!       cmob  (Beta MO coefficient)
! Out : fock1 (Alpha Fock matrix)
!       fock2 (Beta Fock matrix)
!       fock3 (Work space)
!
      use modmolecule, only : natom, neleca, nelecb
      use moddft, only : nrad, nleb
      use modbasis, only : nao
      use modthresh, only : threshweight, threshrho, threshdfock, threshdftao
      implicit none
      integer,intent(in) :: idftex, idftcor, nproc, myrank, mpi_comm
      integer :: ngridatom, iatom, irad, ileb, icount, ilebstart, jatom, imo
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00
      real(8),intent(in) :: cmoa(nao,neleca), cmob(nao,nelecb), atomvec(5,natom,natom)
      real(8),intent(in) :: radpt(2,nrad), angpt(4,nleb), rad(natom), ptweight(nleb,nrad,natom)
      real(8),intent(out) :: fockd1(nao*(nao+1)/2), fockd2(nao*(nao+1)/2)
      real(8),intent(out) :: fockd3(nao*(nao+1)/2), energy, totalelec, vao(nao,4)
      real(8),intent(out) :: vmoa(neleca,4), vmob(nelecb,4), xyzpt(3,natom), rsqrd(natom)
      real(8),intent(out) :: transcmoa(neleca,nao), transcmob(nelecb,nao), work(nao*2)
      real(8) :: weight, rhoa, rhob, grhoa(3), grhob(3), excora(4), excorb(4)
      real(8) :: radpoint, tmp(2,2), wcutoff, rcutoff, fcutoff, aocutoff
!
      transcmoa= transpose(cmoa)
      transcmob= transpose(cmob)
      fockd2(:)= zero
      fockd3(:)= zero
!
      energy= zero
      totalelec= zero
      ngridatom= nrad*nleb
!
      wcutoff= threshweight/(natom*ngridatom)
      rcutoff= threshrho/(natom*ngridatom)
      fcutoff= threshdfock/(natom*ngridatom)
      aocutoff=threshdftao/(natom*ngridatom)
!
!$OMP parallel do collapse(2) schedule(dynamic,1) private(icount,ilebstart,xyzpt, &
!$OMP rsqrd,weight,rhoa,rhob,grhoa,grhob,excora,excorb,radpoint,vao,vmoa,vmob,work) &
!$OMP reduction(+:energy,fockd2,fockd3,totalelec)
      do iatom= 1,natom
        do irad= 1,nrad
          icount=(iatom-1)*ngridatom+(irad-1)*nleb+1+myrank
          ilebstart=mod(icount,nproc)+1
          radpoint= rad(iatom)*radpt(1,irad)
          do ileb= ilebstart,nleb,nproc
!
            weight=ptweight(ileb,irad,iatom)
            if(weight < wcutoff) cycle
            do jatom= 1,natom
              xyzpt(1,jatom)= atomvec(1,iatom,jatom)+radpoint*angpt(1,ileb)
              xyzpt(2,jatom)= atomvec(2,iatom,jatom)+radpoint*angpt(2,ileb)
              xyzpt(3,jatom)= atomvec(3,iatom,jatom)+radpoint*angpt(3,ileb)
              rsqrd(jatom)= xyzpt(1,jatom)*xyzpt(1,jatom)+xyzpt(2,jatom)*xyzpt(2,jatom) &
&                          +xyzpt(3,jatom)*xyzpt(3,jatom)
            enddo
!
            call griduaomo(vao,vmoa,vmob,transcmoa,transcmob,xyzpt,rsqrd,aocutoff)
!
            rhoa= zero
            grhoa(1:3)= zero
            do imo= 1,neleca
              rhoa=     rhoa    +vmoa(imo,1)*vmoa(imo,1)
              grhoa(1)= grhoa(1)+vmoa(imo,2)*vmoa(imo,1)
              grhoa(2)= grhoa(2)+vmoa(imo,3)*vmoa(imo,1)
              grhoa(3)= grhoa(3)+vmoa(imo,4)*vmoa(imo,1)
            enddo
            rhob= zero
            grhob(1:3)= zero
            do imo= 1,nelecb
              rhob=     rhob    +vmob(imo,1)*vmob(imo,1)
              grhob(1)= grhob(1)+vmob(imo,2)*vmob(imo,1)
              grhob(2)= grhob(2)+vmob(imo,3)*vmob(imo,1)
              grhob(3)= grhob(3)+vmob(imo,4)*vmob(imo,1)
            enddo
            if((rhoa+rhob) < rcutoff)cycle
            grhoa(1)= grhoa(1)*two
            grhoa(2)= grhoa(2)*two
            grhoa(3)= grhoa(3)*two
            grhob(1)= grhob(1)*two
            grhob(2)= grhob(2)*two
            grhob(3)= grhob(3)*two
            call calcexcor(excora,excorb,energy,rhoa,rhob,grhoa,grhob,weight,idftex,idftcor,2)
            call ufockexcor(fockd2,fockd3,excora,excorb,vao,vao(1,2),work,weight,nao,fcutoff)
            totalelec=totalelec+weight*(rhoa+rhob)
          enddo
        enddo
      enddo
!$OMP end parallel do
!
      tmp(1,1)= energy
      tmp(2,1)= totalelec
      call para_allreducer(fockd2,fockd1,nao*(nao+1)/2,mpi_comm)
      call para_allreducer(fockd3,fockd2,nao*(nao+1)/2,mpi_comm)
      call para_allreducer(tmp(1,1),tmp(1,2),2,mpi_comm)
      energy    = tmp(1,2)
      totalelec = tmp(2,2)
      return
end


!--------------------------------------------------------------------------------------
  subroutine ufockexcor(fockda,fockdb,excora,excorb,vao,vgao,work,weight,nao,fcutoff)
!--------------------------------------------------------------------------------------
!
! Add exchange-correlation term to unrestricted Fock matrix at a grid point
!
      implicit none
      integer,intent(in) :: nao
      integer :: ij, i, j
      real(8),parameter :: zero=0.0D+00, half=0.5D+00
      real(8),intent(in) :: excora(4), excorb(4), vao(nao), vgao(nao,3), weight, fcutoff
      real(8),intent(inout) :: fockda(nao*(nao+1)/2), fockdb(nao*(nao+1)/2), work(nao,2)
      real(8) :: vaomax, gridamax, gridbmax, vaoabs, gridaabs, gridbabs
!
      ij= 0
      vaomax= zero
      gridamax= zero
      gridbmax= zero
      do i= 1,nao
        work(i,1)=(excora(1)*vao(i)*half+excora(2)*vgao(i,1) &
&                 +excora(3)*vgao(i,2)  +excora(4)*vgao(i,3))*weight
        work(i,2)=(excorb(1)*vao(i)*half+excorb(2)*vgao(i,1) &
&                 +excorb(3)*vgao(i,2)  +excorb(4)*vgao(i,3))*weight
        vaoabs= abs(vao(i))
        gridaabs= abs(work(i,1))
        gridbabs= abs(work(i,2))
        vaomax= max(vaomax,vaoabs)
        gridamax= max(gridamax,gridaabs)
        gridbmax= max(gridbmax,gridbabs)
        if((vaoabs*(gridamax+gridbmax)+(gridaabs+gridbabs)*vaomax) >= fcutoff) then
          do j= 1,i
            ij= ij+1
            fockda(ij)= fockda(ij)+vao(i)*work(j,1)+work(i,1)*vao(j)
            fockdb(ij)= fockdb(ij)+vao(i)*work(j,2)+work(i,2)*vao(j)
          enddo
        else
          ij= ij+i
        endif
      enddo
      return
end


!------------------------------------------
  subroutine calcatomvec(atomvec,surface)
!------------------------------------------
!
! Calculate atom vectors and surface shifting parameters
!
! Out : atomvec (atom vector and distance)
!       surface (surface shifting parameters)
!
      use modmolecule, only : natom, coord, numatomic
      use modatom, only : atomrad
      implicit none
      integer :: iatom, jatom, inum, jnum
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00
      real(8),intent(out) :: atomvec(5,natom,natom), surface(natom,natom)
      real(8) :: tmp1, tmp2
!
! Calculate atom vectors
!
!$OMP parallel do
      do iatom= 1,natom
        atomvec(1,iatom,iatom)= zero
        atomvec(2,iatom,iatom)= zero
        atomvec(3,iatom,iatom)= zero
        atomvec(4,iatom,iatom)= zero
        do jatom= 1,iatom-1
          atomvec(1,jatom,iatom)= coord(1,jatom)-coord(1,iatom)
          atomvec(2,jatom,iatom)= coord(2,jatom)-coord(2,iatom)
          atomvec(3,jatom,iatom)= coord(3,jatom)-coord(3,iatom)
          atomvec(4,jatom,iatom)= sqrt(atomvec(1,jatom,iatom)*atomvec(1,jatom,iatom) &
                                      +atomvec(2,jatom,iatom)*atomvec(2,jatom,iatom) &
                                      +atomvec(3,jatom,iatom)*atomvec(3,jatom,iatom))
          atomvec(5,jatom,iatom)= one/atomvec(4,jatom,iatom)
          atomvec(1,iatom,jatom)=-atomvec(1,jatom,iatom)
          atomvec(2,iatom,jatom)=-atomvec(2,jatom,iatom)
          atomvec(3,iatom,jatom)=-atomvec(3,jatom,iatom)
          atomvec(4,iatom,jatom)= atomvec(4,jatom,iatom)
          atomvec(5,iatom,jatom)= atomvec(5,jatom,iatom)
        enddo
      enddo
!$OMP end parallel do
!
! Calculate surface shifting parameters
!
!$OMP parallel do private(inum,jnum,tmp1,tmp2)
      do iatom= 1,natom
        surface(iatom,iatom)= zero
        inum= numatomic(iatom)
        if(inum == 0) then
          surface(1,iatom)= -one
          cycle
        endif
        do jatom= 1,iatom-1
          jnum= numatomic(jatom)
          if(jnum == 0) then
            surface(jatom,iatom)= one
            cycle
          endif
          tmp1= atomrad(inum)/atomrad(jnum)
          tmp2=(tmp1-one)/(tmp1+one)
          surface(jatom,iatom)= tmp2/(tmp2*tmp2-one)
          if(surface(jatom,iatom) > half) surface(jatom,iatom)= half
          if(surface(jatom,iatom) <-half) surface(jatom,iatom)=-half
        enddo
        do jatom= iatom+1,natom
          jnum= numatomic(jatom)
          if(jnum == 0) then
            surface(jatom,iatom)= one
            cycle
          endif
          tmp1= atomrad(inum)/atomrad(jnum)
          tmp2=(tmp1-one)/(tmp1+one)
          surface(jatom,iatom)= tmp2/(tmp2*tmp2-one)
          if(surface(jatom,iatom) > half) surface(jatom,iatom)= half
          if(surface(jatom,iatom) <-half) surface(jatom,iatom)=-half
        enddo
      enddo
!$OMP end parallel do
!
      return
end


!-----------------------------------
  subroutine calcradpt(radpt,nrad)
!-----------------------------------
!
! Calculate radial points and weights
!
! In  : nradpt (number of radial points)
! Out : radpt(1,*) (radial points)
!       radpt(2,*) (radial weights)
!
      implicit none
      integer,intent(in) :: nrad
      integer :: irad
      real(8),parameter :: one=1.0D+00, two=2.0D+00
      real(8),intent(out) :: radpt(2,nrad)
      real(8) :: xnrad1, xirad, xirad2, tmp
!
      xnrad1= nrad+one
!$OMP parallel do private(xirad,xirad2,tmp)
      do irad= 1,nrad
        xirad= xnrad1-irad
        xirad2= xirad*xirad
        radpt(1,irad)=(irad*irad)/xirad2
        tmp= xirad2*xirad2*xirad2*xirad
        radpt(2,irad)=two*xnrad1*irad*irad*irad*irad*irad/tmp
      enddo
!$OMP end parallel do
      return
end


!-----------------------------------
  subroutine calclebpt(angpt,nleb)
!-----------------------------------
!
! Calculate Levedev quadrature points and weights
!
! In  : nleb (number of Lebedev quadrature points)
! Out : angpt(1:3,*) (Levedev points)
!       angpt(4  ,*) (Levedev weights)
!
      implicit none
      integer,intent(in) :: nleb
      integer :: ileb
      real(8),parameter :: pi4=1.256637061435917D+01
      real(8),intent(out) :: angpt(4,nleb)
!
!
! Omit nleb= 74, 230, 266 because of negative weights
! Number of electrons and energy are not accurate.
!
      select case(nleb)
        case(6)
          call lebedev6(angpt)
        case(14)
          call lebedev14(angpt)
        case(26)
          call lebedev26(angpt)
        case(38)
          call lebedev38(angpt)
        case(50)
          call lebedev50(angpt)
!       case(74)
!         call lebedev74(angpt)
        case(86)
          call lebedev86(angpt)
        case(110)
          call lebedev110(angpt)
        case(146)
          call lebedev146(angpt)
        case(170)
          call lebedev170(angpt)
        case(194)
          call lebedev194(angpt)
!       case(230)
!         call lebedev230(angpt)
!       case(266)
!         call lebedev266(angpt)
        case(302)
          call lebedev302(angpt)
        case(350)
          call lebedev350(angpt)
        case(434)
          call lebedev434(angpt)
        case(590)
          call lebedev590(angpt)
        case(770)
          call lebedev770(angpt)
        case(974)
          call lebedev974(angpt)
        case(1202)
          call lebedev1202(angpt)
        case(1454)
          call lebedev1454(angpt)
        case default
          write(*,'(" Error! Nleb=",i4," is not supported. ")')nleb
          call iabort
      end select
!
      do ileb= 1,nleb
        angpt(4,ileb)= pi4*angpt(4,ileb)
      enddo
      return
end


!----------------------------------------------------------------------------------------------
  subroutine calcgridweight(ptweight,rad,radpt,angpt,atomvec,surface,xyzpt,work,nproc,myrank)
!----------------------------------------------------------------------------------------------
!
! Calculate weights of grid points
!
! In  : rad (atom radius)
!       radpt (radial point)
!       angpt (angular point)
!       atomvec (tom vector and distance)
!       surface (surface shifting parameter)
! Out : ptweight  (weight of grid point)
!       xyzpt, work (work space)
!
      use moddft, only : nrad, nleb
      use modmolecule, only : natom
      implicit none
      integer,intent(in) :: nproc, myrank
      integer :: katom, irad, ileb, iatom, jatom, i, icount, ilebstart, ngridatom
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, oneh=1.5D+00
      real(8),intent(in) :: rad(natom), radpt(2,nrad), angpt(4,nleb), atomvec(5,natom,natom)
      real(8),intent(in) :: surface(natom,natom)
      real(8),intent(out) :: ptweight(nleb,nrad,natom), xyzpt(3,natom), work(natom,2)
      real(8) :: xyzgrid(3), radpoint, radweight, wttot, cutij, cutji, xmuij, zmuij, f4, f2
!
      ngridatom= nrad*nleb
!
!$OMP parallel do collapse(2) private(icount,ilebstart,radpoint,radweight,wttot,xyzpt,&
!$OMP work,xyzgrid,cutij,cutji,zmuij,xmuij,f4,f2)
      do katom= 1,natom
        do irad= 1,nrad
          icount=(katom-1)*ngridatom+(irad-1)*nleb+1+myrank
          ilebstart=mod(icount,nproc)+1
          radpoint= rad(katom)*radpt(1,irad)
          radweight= rad(katom)*rad(katom)*rad(katom)*radpt(2,irad)
          do ileb= ilebstart,nleb,nproc
            wttot= zero
!
            xyzgrid(1)= radpoint*angpt(1,ileb)
            xyzgrid(2)= radpoint*angpt(2,ileb)
            xyzgrid(3)= radpoint*angpt(3,ileb)
            do iatom= 1,natom
              xyzpt(1,iatom)= atomvec(1,katom,iatom)+xyzgrid(1)
              xyzpt(2,iatom)= atomvec(2,katom,iatom)+xyzgrid(2)
              xyzpt(3,iatom)= atomvec(3,katom,iatom)+xyzgrid(3)
              work(iatom,1)= one
              work(iatom,2)= sqrt(xyzpt(1,iatom)*xyzpt(1,iatom)+xyzpt(2,iatom)*xyzpt(2,iatom) &
&                                +xyzpt(3,iatom)*xyzpt(3,iatom))
            enddo
!
            do iatom= 1,natom
              if(surface(1,iatom) == -one) then
                work(iatom,1)= zero
                cycle
              endif
              do jatom= 1,iatom-1
                if(surface(jatom,iatom) == one) cycle
                zmuij=(work(iatom,2)-work(jatom,2))*atomvec(5,iatom,jatom)
                xmuij= zmuij+surface(jatom,iatom)*(one-zmuij*zmuij)
                f4= xmuij
                do i= 1,4
                  f4=f4*(oneh-half*f4*f4)
                enddo
                f2= half*f4
                cutij= half-f2
                cutji= half+f2
                work(iatom,1)=work(iatom,1)*cutij
                work(jatom,1)=work(jatom,1)*cutji
              enddo
            enddo
            do iatom= 1,natom
              wttot= wttot+work(iatom,1)
            enddo
            ptweight(ileb,irad,katom)= work(katom,1)*radweight*angpt(4,ileb)/wttot
!
          enddo
        enddo
      enddo
      return
end


!--------------------------------------------------------------
  subroutine gridraomo(vao,vmo,transcmo,xyzpt,rsqrd,aocutoff)
!--------------------------------------------------------------
!
! Calculate closed-shell AO and MO values for a grid point
!
      use modmolecule, only : natom, neleca
      use modbasis, only : nao
      implicit none
      integer :: ii, jj
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: xyzpt(3,natom), rsqrd(natom), transcmo(neleca,nao), aocutoff
      real(8),intent(out) :: vao(nao,4), vmo(neleca,4)
!
! Calculate AO values for a grid point
!
      vao(:,:)= zero
      call gridao(vao,xyzpt,rsqrd)
!
! Calculate MO values for a grid point
!
      vmo(:,:)= zero
      do ii= 1,nao
        if(abs(vao(ii,1))+abs(vao(ii,2)) &
&         +abs(vao(ii,3))+abs(vao(ii,4)) > aocutoff) then
          do jj= 1,neleca
            vmo(jj,1)= vmo(jj,1)+vao(ii,1)*transcmo(jj,ii)
            vmo(jj,2)= vmo(jj,2)+vao(ii,2)*transcmo(jj,ii)
            vmo(jj,3)= vmo(jj,3)+vao(ii,3)*transcmo(jj,ii)
            vmo(jj,4)= vmo(jj,4)+vao(ii,4)*transcmo(jj,ii)
          enddo
        endif
      enddo
!
      return
end


!-------------------------------------------------------------------------------
  subroutine griduaomo(vao,vmoa,vmob,transcmoa,transcmob,xyzpt,rsqrd,aocutoff)
!-------------------------------------------------------------------------------
!
! Calculate open-shell AO and MO values for a grid point
!
      use modmolecule, only : natom, neleca, nelecb
      use modbasis, only : nao
      implicit none
      integer :: ii, jj
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: xyzpt(3,natom), rsqrd(natom), transcmoa(neleca,nao)
      real(8),intent(in) :: transcmob(nelecb,nao), aocutoff
      real(8),intent(out) :: vao(nao,4), vmoa(neleca,4), vmob(nelecb,4)
!
      vao(:,:)= zero
!
      call gridao(vao,xyzpt,rsqrd)
!
      vmoa(:,:)= zero
      vmob(:,:)= zero
      do ii= 1,nao
        if(abs(vao(ii,1))+abs(vao(ii,2)) &
&         +abs(vao(ii,3))+abs(vao(ii,4)) > aocutoff) then
          do jj= 1,neleca
            vmoa(jj,1)= vmoa(jj,1)+vao(ii,1)*transcmoa(jj,ii)
            vmoa(jj,2)= vmoa(jj,2)+vao(ii,2)*transcmoa(jj,ii)
            vmoa(jj,3)= vmoa(jj,3)+vao(ii,3)*transcmoa(jj,ii)
            vmoa(jj,4)= vmoa(jj,4)+vao(ii,4)*transcmoa(jj,ii)
          enddo
          do jj= 1,nelecb
            vmob(jj,1)= vmob(jj,1)+vao(ii,1)*transcmob(jj,ii)
            vmob(jj,2)= vmob(jj,2)+vao(ii,2)*transcmob(jj,ii)
            vmob(jj,3)= vmob(jj,3)+vao(ii,3)*transcmob(jj,ii)
            vmob(jj,4)= vmob(jj,4)+vao(ii,4)*transcmob(jj,ii)
          enddo
        endif
      enddo
!
      return
end


!-------------------------------------
  subroutine gridao(vao,xyzpt,rsqrd)
!-------------------------------------
!
! Calculate AO values for a grid point
!
      use modmolecule, only : natom
      use modbasis, only : ex, coeff, nshell, nao, locbf, locatom, &
&                          mprim, mbf, mtype
      use modthresh, only : threshex
      implicit none
      integer :: icount, ish, numprim, iatom, iprim, nang, nbf, ilocbf, ii, jj
      real(8),parameter :: half=0.5D+00, one=1.0D+00, two=2.0D+00
      real(8),parameter :: three=3.0D+00, four=4.0D+00, five=5.0D+00, six=6.0D+00
      real(8),parameter :: eight=8.0D+00, p9=9.0D+00, ten=10.0D+00, twelve=12.0D+00
      real(8),parameter :: p15=15.0D+00, p16=16.0D+00, p18=18.0D+00, p20=20.0D+00
      real(8),parameter :: p24=24.0D+00, p30=30.0D+00, p32=32.0D+00, p40=40.0D+00
      real(8),parameter :: p60=60.0D+00, p90=90.0D+00, p120=1.2D+02, p180=1.8D+02
      real(8),parameter :: eighth=0.125D+00, sixteenth=6.25D-02
      real(8),parameter :: sqrt3=1.73205080756888D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrt5=2.236067977499790D+00, sqrt15=3.872983346207417D+00
      real(8),parameter :: sqrt7=2.645751311064590D+00, sqrt35=5.916079783099616D+00
      real(8),parameter :: sqrt35third=3.415650255319866D+00
      real(8),parameter :: sqrt21=4.582575694955840D+00, sqrt63=7.937253933193772D+00
      real(8),parameter :: sqrt105=1.024695076595960D+01, sqrt11=3.316624790355400D+00
      real(8),parameter :: sqrt33=5.744562646538029D+00, sqrt99=9.949874371066200D+00
      real(8),parameter :: sqrt231=1.519868415357066D+01, sqrt231fifth=6.797058187186571D+00
      real(8),parameter :: sqrt385=1.962141687034858D+01
      real(8),parameter :: facf1=0.79056941504209483D+00 ! sqrt(5/2)/2
      real(8),parameter :: facf2=3.87298334620741688D+00 ! sqrt(15)
      real(8),parameter :: facf3=0.61237243569579452D+00 ! sqrt(3/2)/2
      real(8),parameter :: facf4=1.93649167310370844D+00 ! sqrt(15)/2
      real(8),parameter :: facg1=2.95803989154980802D+00 ! sqrt(35)/2
      real(8),parameter :: facg2=2.09165006633518887D+00 ! sqrt(35/2)/2
      real(8),parameter :: facg3=1.11803398874989484D+00 ! sqrt(5)/2
      real(8),parameter :: facg4=0.79056941504209483D+00 ! sqrt(5/2)/2
      real(8),parameter :: facg5=0.55901699437494742D+00 ! sqrt(5)/4
      real(8),parameter :: facg6=0.73950997288745200D+00 ! sqrt(35)/8
      real(8),parameter :: fach1=0.70156076002011400D+00 ! sqrt(63/2)/8
      real(8),parameter :: fach2=2.21852991866235601D+00 ! sqrt(315)/8
      real(8),parameter :: fach3=0.52291251658379721D+00 ! sqrt(35/2)/8
      real(8),parameter :: fach4=2.56173769148989959D+00 ! sqrt(105)/4
      real(8),parameter :: fach5=0.48412291827592711D+00 ! sqrt(15)/8
      real(8),parameter :: faci1=0.67169328938139615D+00 ! sqrt(231/2)/16
      real(8),parameter :: faci2=2.32681380862328561D+00 ! sqrt(693/2)/8
      real(8),parameter :: faci3=0.49607837082461073D+00 ! sqrt(63)/16
      real(8),parameter :: faci4=0.90571104663683991D+00 ! sqrt(105/2)/8
      real(8),parameter :: faci5=0.45285552331841995D+00 ! sqrt(105/2)/16
      real(8),parameter :: faci6=0.57282196186948000D+00 ! sqrt(21)/8
      real(8),intent(in) :: xyzpt(3,natom), rsqrd(natom)
      real(8),intent(inout) :: vao(nao,4)
      real(8) :: expval, fac, tmp(28,4), xx, yy, zz, xy, xz, yz
      real(8) :: xyz3(10), xyz4(15), xyz5(21), xyz6(28), xyz7(36)
!
      icount= 0
      do ish= 1,nshell
        nang= mtype(ish)
        numprim= mprim(ish)
        nbf = mbf(ish)
        ilocbf= locbf(ish)
        iatom= locatom(ish)
        select case(nang)
!
! S function
!
          case(0)
            do iprim= 1,numprim
              icount= icount+1
              if(ex(icount)*rsqrd(iatom) > threshex) cycle
              expval= exp(-ex(icount)*rsqrd(iatom))*coeff(icount)
              fac= ex(icount)*expval*two
              vao(ilocbf+1,1)= vao(ilocbf+1,1)+expval
              vao(ilocbf+1,2)= vao(ilocbf+1,2)-fac*xyzpt(1,iatom)
              vao(ilocbf+1,3)= vao(ilocbf+1,3)-fac*xyzpt(2,iatom)
              vao(ilocbf+1,4)= vao(ilocbf+1,4)-fac*xyzpt(3,iatom)
            enddo
!
! P function
!
          case(1)
            do iprim= 1,numprim
              icount= icount+1
              if(ex(icount)*rsqrd(iatom) > threshex) cycle
              expval= exp(-ex(icount)*rsqrd(iatom))*coeff(icount)
              fac= ex(icount)*expval*two
              vao(ilocbf+1,1)= vao(ilocbf+1,1)+expval*xyzpt(1,iatom)
              vao(ilocbf+2,1)= vao(ilocbf+2,1)+expval*xyzpt(2,iatom)
              vao(ilocbf+3,1)= vao(ilocbf+3,1)+expval*xyzpt(3,iatom)
              xx= xyzpt(1,iatom)*xyzpt(1,iatom)
              yy= xyzpt(2,iatom)*xyzpt(2,iatom)
              zz= xyzpt(3,iatom)*xyzpt(3,iatom)
              xy= xyzpt(1,iatom)*xyzpt(2,iatom)
              xz= xyzpt(1,iatom)*xyzpt(3,iatom)
              yz= xyzpt(2,iatom)*xyzpt(3,iatom)
              vao(ilocbf+1,2)= vao(ilocbf+1,2)+expval-fac*xx
              vao(ilocbf+2,2)= vao(ilocbf+2,2)       -fac*xy
              vao(ilocbf+3,2)= vao(ilocbf+3,2)       -fac*xz
              vao(ilocbf+1,3)= vao(ilocbf+1,3)       -fac*xy
              vao(ilocbf+2,3)= vao(ilocbf+2,3)+expval-fac*yy
              vao(ilocbf+3,3)= vao(ilocbf+3,3)       -fac*yz
              vao(ilocbf+1,4)= vao(ilocbf+1,4)       -fac*xz
              vao(ilocbf+2,4)= vao(ilocbf+2,4)       -fac*yz
              vao(ilocbf+3,4)= vao(ilocbf+3,4)+expval-fac*zz
            enddo
!
! D function
!
          case(2)
            do iprim= 1,numprim
              icount= icount+1
              if(ex(icount)*rsqrd(iatom) > threshex) cycle
              expval= exp(-ex(icount)*rsqrd(iatom))*coeff(icount)
              fac= ex(icount)*two
              tmp(1,1)= expval*xyzpt(1,iatom)*xyzpt(1,iatom)
              tmp(2,1)= expval*xyzpt(1,iatom)*xyzpt(2,iatom)*sqrt3
              tmp(3,1)= expval*xyzpt(1,iatom)*xyzpt(3,iatom)*sqrt3
              tmp(4,1)= expval*xyzpt(2,iatom)*xyzpt(2,iatom)
              tmp(5,1)= expval*xyzpt(2,iatom)*xyzpt(3,iatom)*sqrt3
              tmp(6,1)= expval*xyzpt(3,iatom)*xyzpt(3,iatom)
              xyz3( 1)=-xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(1,iatom)*fac
              xyz3( 2)=-xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(2,iatom)*fac
              xyz3( 3)=-xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(3,iatom)*fac
              xyz3( 4)=-xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)*fac
              xyz3( 5)=-xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)*fac
              xyz3( 6)=-xyzpt(1,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac
              xyz3( 7)=-xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)*fac
              xyz3( 8)=-xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)*fac
              xyz3( 9)=-xyzpt(2,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac
              xyz3(10)=-xyzpt(3,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac
              tmp(1,2)= expval*(xyzpt(1,iatom)*two+xyz3( 1))
              tmp(2,2)= expval*(xyzpt(2,iatom)    +xyz3( 2))*sqrt3
              tmp(3,2)= expval*(xyzpt(3,iatom)    +xyz3( 3))*sqrt3
              tmp(4,2)= expval*(                  +xyz3( 4))
              tmp(5,2)= expval*(                  +xyz3( 5))*sqrt3
              tmp(6,2)= expval*(                  +xyz3( 6))
              tmp(1,3)= expval*(                  +xyz3( 2))
              tmp(2,3)= expval*(xyzpt(1,iatom)    +xyz3( 4))*sqrt3
              tmp(3,3)= expval*(                  +xyz3( 5))*sqrt3
              tmp(4,3)= expval*(xyzpt(2,iatom)*two+xyz3( 7))
              tmp(5,3)= expval*(xyzpt(3,iatom)    +xyz3( 8))*sqrt3
              tmp(6,3)= expval*(                  +xyz3( 9))
              tmp(1,4)= expval*(                  +xyz3( 3))
              tmp(2,4)= expval*(                  +xyz3( 5))*sqrt3
              tmp(3,4)= expval*(xyzpt(1,iatom)    +xyz3( 6))*sqrt3
              tmp(4,4)= expval*(                  +xyz3( 8))
              tmp(5,4)= expval*(xyzpt(2,iatom)    +xyz3( 9))*sqrt3
              tmp(6,4)= expval*(xyzpt(3,iatom)*two+xyz3(10))
              if(nbf == 6) then
                do jj= 1,4
                  do ii= 1,6
                    vao(ilocbf+ii,jj)= vao(ilocbf+ii,jj)+tmp(ii,jj)
                  enddo
                enddo
              else
                do jj= 1,4
                  vao(ilocbf+1,jj)= vao(ilocbf+1,jj)+tmp(2,jj)
                  vao(ilocbf+2,jj)= vao(ilocbf+2,jj)+tmp(5,jj)
                  vao(ilocbf+3,jj)= vao(ilocbf+3,jj)+tmp(6,jj)-(tmp(1,jj)+tmp(4,jj))*half
                  vao(ilocbf+4,jj)= vao(ilocbf+4,jj)+tmp(3,jj)
                  vao(ilocbf+5,jj)= vao(ilocbf+5,jj)+(tmp(1,jj)-tmp(4,jj))*sqrt3h
                enddo
              endif
            enddo
!
! F function
!
          case(3)
            do iprim= 1,numprim
              icount= icount+1
              if(ex(icount)*rsqrd(iatom) > threshex) cycle
              expval= exp(-ex(icount)*rsqrd(iatom))*coeff(icount)
              fac= ex(icount)*two
              xx= xyzpt(1,iatom)*xyzpt(1,iatom)
              yy= xyzpt(2,iatom)*xyzpt(2,iatom)
              zz= xyzpt(3,iatom)*xyzpt(3,iatom)
              xy= xyzpt(1,iatom)*xyzpt(2,iatom)
              xz= xyzpt(1,iatom)*xyzpt(3,iatom)
              yz= xyzpt(2,iatom)*xyzpt(3,iatom)
              xyz4( 1)=-xx*xx*fac
              xyz4( 2)=-xx*xy*fac
              xyz4( 3)=-xx*xz*fac
              xyz4( 4)=-xx*yy*fac
              xyz4( 5)=-xx*yz*fac
              xyz4( 6)=-xx*zz*fac
              xyz4( 7)=-xy*yy*fac
              xyz4( 8)=-xy*yz*fac
              xyz4( 9)=-xy*zz*fac
              xyz4(10)=-xz*zz*fac
              xyz4(11)=-yy*yy*fac
              xyz4(12)=-yy*yz*fac
              xyz4(13)=-yy*zz*fac
              xyz4(14)=-yz*zz*fac
              xyz4(15)=-zz*zz*fac
              tmp( 1,1)= expval*xx*xyzpt(1,iatom)
              tmp( 2,1)= expval*xx*xyzpt(2,iatom)
              tmp( 3,1)= expval*xx*xyzpt(3,iatom)
              tmp( 4,1)= expval*xy*xyzpt(2,iatom)
              tmp( 5,1)= expval*xy*xyzpt(3,iatom)
              tmp( 6,1)= expval*xz*xyzpt(3,iatom)
              tmp( 7,1)= expval*yy*xyzpt(2,iatom)
              tmp( 8,1)= expval*yy*xyzpt(3,iatom)
              tmp( 9,1)= expval*yz*xyzpt(3,iatom)
              tmp(10,1)= expval*zz*xyzpt(3,iatom)
              tmp( 1,2)= expval*(xx*three+xyz4( 1))
              tmp( 2,2)= expval*(xy*two  +xyz4( 2))
              tmp( 3,2)= expval*(xz*two  +xyz4( 3))
              tmp( 4,2)= expval*(yy      +xyz4( 4))
              tmp( 5,2)= expval*(yz      +xyz4( 5))
              tmp( 6,2)= expval*(zz      +xyz4( 6))
              tmp( 7,2)= expval*(        +xyz4( 7))
              tmp( 8,2)= expval*(        +xyz4( 8))
              tmp( 9,2)= expval*(        +xyz4( 9))
              tmp(10,2)= expval*(        +xyz4(10))
              tmp( 1,3)= expval*(        +xyz4( 2))
              tmp( 2,3)= expval*(xx      +xyz4( 4))
              tmp( 3,3)= expval*(        +xyz4( 5))
              tmp( 4,3)= expval*(xy*two  +xyz4( 7))
              tmp( 5,3)= expval*(xz      +xyz4( 8))
              tmp( 6,3)= expval*(        +xyz4( 9))
              tmp( 7,3)= expval*(yy*three+xyz4(11))
              tmp( 8,3)= expval*(yz*two  +xyz4(12))
              tmp( 9,3)= expval*(zz      +xyz4(13))
              tmp(10,3)= expval*(        +xyz4(14))
              tmp( 1,4)= expval*(        +xyz4( 3))
              tmp( 2,4)= expval*(        +xyz4( 5))
              tmp( 3,4)= expval*(xx      +xyz4( 6))
              tmp( 4,4)= expval*(        +xyz4( 8))
              tmp( 5,4)= expval*(xy      +xyz4( 9))
              tmp( 6,4)= expval*(xz*two  +xyz4(10))
              tmp( 7,4)= expval*(        +xyz4(12))
              tmp( 8,4)= expval*(yy      +xyz4(13))
              tmp( 9,4)= expval*(yz*two  +xyz4(14))
              tmp(10,4)= expval*(zz*three+xyz4(15))
              if(nbf == 10) then
                do jj= 1,4
                  vao(ilocbf+ 1,jj)= vao(ilocbf+ 1,jj)+tmp( 1,jj)
                  vao(ilocbf+ 2,jj)= vao(ilocbf+ 2,jj)+tmp( 2,jj)*sqrt5
                  vao(ilocbf+ 3,jj)= vao(ilocbf+ 3,jj)+tmp( 3,jj)*sqrt5
                  vao(ilocbf+ 4,jj)= vao(ilocbf+ 4,jj)+tmp( 4,jj)*sqrt5
                  vao(ilocbf+ 5,jj)= vao(ilocbf+ 5,jj)+tmp( 5,jj)*sqrt15
                  vao(ilocbf+ 6,jj)= vao(ilocbf+ 6,jj)+tmp( 6,jj)*sqrt5
                  vao(ilocbf+ 7,jj)= vao(ilocbf+ 7,jj)+tmp( 7,jj)
                  vao(ilocbf+ 8,jj)= vao(ilocbf+ 8,jj)+tmp( 8,jj)*sqrt5
                  vao(ilocbf+ 9,jj)= vao(ilocbf+ 9,jj)+tmp( 9,jj)*sqrt5
                  vao(ilocbf+10,jj)= vao(ilocbf+10,jj)+tmp(10,jj)
                enddo
              else
                do jj= 1,4
                  vao(ilocbf+1,jj)= vao(ilocbf+1,jj)+(-tmp(7,jj)+three*tmp(2,jj)         )*facf1
                  vao(ilocbf+2,jj)= vao(ilocbf+2,jj)+  tmp(5,jj)                          *facf2
                  vao(ilocbf+3,jj)= vao(ilocbf+3,jj)+(-tmp(7,jj)-tmp(2,jj)+four*tmp(9,jj))*facf3
                  vao(ilocbf+4,jj)= vao(ilocbf+4,jj)+( two*tmp(10,jj)-three*tmp(3,jj) &
&                                                                    -three*tmp(8,jj)    )*half
                  vao(ilocbf+5,jj)= vao(ilocbf+5,jj)+(-tmp(1,jj)-tmp(4,jj)+four*tmp(6,jj))*facf3
                  vao(ilocbf+6,jj)= vao(ilocbf+6,jj)+( tmp(3,jj)-tmp(8,jj)               )*facf4
                  vao(ilocbf+7,jj)= vao(ilocbf+7,jj)+( tmp(1,jj)-three*tmp(4,jj)         )*facf1
                enddo
              endif
            enddo
!
! G function
!
          case(4)
            do iprim= 1,numprim
              icount= icount+1
              if(ex(icount)*rsqrd(iatom) > threshex) cycle
              expval= exp(-ex(icount)*rsqrd(iatom))*coeff(icount)
              fac= ex(icount)*two
              xx= xyzpt(1,iatom)*xyzpt(1,iatom)
              xy= xyzpt(1,iatom)*xyzpt(2,iatom)
              xz= xyzpt(1,iatom)*xyzpt(3,iatom)
              yy= xyzpt(2,iatom)*xyzpt(2,iatom)
              yz= xyzpt(2,iatom)*xyzpt(3,iatom)
              zz= xyzpt(3,iatom)*xyzpt(3,iatom)
              xyz3( 1)= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(1,iatom)
              xyz3( 2)= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(2,iatom)
              xyz3( 3)= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(3,iatom)
              xyz3( 4)= xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)
              xyz3( 5)= xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)
              xyz3( 6)= xyzpt(1,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)
              xyz3( 7)= xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)
              xyz3( 8)= xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)
              xyz3( 9)= xyzpt(2,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)
              xyz3(10)= xyzpt(3,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)
              xyz5( 1)=-xyz3( 1)*xx*fac
              xyz5( 2)=-xyz3( 1)*xy*fac
              xyz5( 3)=-xyz3( 1)*xz*fac
              xyz5( 4)=-xyz3( 1)*yy*fac
              xyz5( 5)=-xyz3( 1)*yz*fac
              xyz5( 6)=-xyz3( 1)*zz*fac
              xyz5( 7)=-xyz3( 2)*yy*fac
              xyz5( 8)=-xyz3( 2)*yz*fac
              xyz5( 9)=-xyz3( 2)*zz*fac
              xyz5(10)=-xyz3( 3)*zz*fac
              xyz5(11)=-xyz3( 4)*yy*fac
              xyz5(12)=-xyz3( 4)*yz*fac
              xyz5(13)=-xyz3( 4)*zz*fac
              xyz5(14)=-xyz3( 5)*zz*fac
              xyz5(15)=-xyz3( 6)*zz*fac
              xyz5(16)=-xyz3( 7)*yy*fac
              xyz5(17)=-xyz3( 7)*yz*fac
              xyz5(18)=-xyz3( 7)*zz*fac
              xyz5(19)=-xyz3( 8)*zz*fac
              xyz5(20)=-xyz3( 9)*zz*fac
              xyz5(21)=-xyz3(10)*zz*fac
              tmp( 1,1)= expval*xyz3( 1)*xyzpt(1,iatom)
              tmp( 2,1)= expval*xyz3( 1)*xyzpt(2,iatom)
              tmp( 3,1)= expval*xyz3( 1)*xyzpt(3,iatom)
              tmp( 4,1)= expval*xyz3( 2)*xyzpt(2,iatom)
              tmp( 5,1)= expval*xyz3( 2)*xyzpt(3,iatom)
              tmp( 6,1)= expval*xyz3( 3)*xyzpt(3,iatom)
              tmp( 7,1)= expval*xyz3( 4)*xyzpt(2,iatom)
              tmp( 8,1)= expval*xyz3( 4)*xyzpt(3,iatom)
              tmp( 9,1)= expval*xyz3( 5)*xyzpt(3,iatom)
              tmp(10,1)= expval*xyz3( 6)*xyzpt(3,iatom)
              tmp(11,1)= expval*xyz3( 7)*xyzpt(2,iatom)
              tmp(12,1)= expval*xyz3( 7)*xyzpt(3,iatom)
              tmp(13,1)= expval*xyz3( 8)*xyzpt(3,iatom)
              tmp(14,1)= expval*xyz3( 9)*xyzpt(3,iatom)
              tmp(15,1)= expval*xyz3(10)*xyzpt(3,iatom)
              tmp( 1,2)= expval*(xyz3( 1)*four +xyz5( 1))
              tmp( 2,2)= expval*(xyz3( 2)*three+xyz5( 2))
              tmp( 3,2)= expval*(xyz3( 3)*three+xyz5( 3))
              tmp( 4,2)= expval*(xyz3( 4)*two  +xyz5( 4))
              tmp( 5,2)= expval*(xyz3( 5)*two  +xyz5( 5))
              tmp( 6,2)= expval*(xyz3( 6)*two  +xyz5( 6))
              tmp( 7,2)= expval*(xyz3( 7)      +xyz5( 7))
              tmp( 8,2)= expval*(xyz3( 8)      +xyz5( 8))
              tmp( 9,2)= expval*(xyz3( 9)      +xyz5( 9))
              tmp(10,2)= expval*(xyz3(10)      +xyz5(10))
              tmp(11,2)= expval*(              +xyz5(11))
              tmp(12,2)= expval*(              +xyz5(12))
              tmp(13,2)= expval*(              +xyz5(13))
              tmp(14,2)= expval*(              +xyz5(14))
              tmp(15,2)= expval*(              +xyz5(15))
              tmp( 1,3)= expval*(              +xyz5( 2))
              tmp( 2,3)= expval*(xyz3( 1)      +xyz5( 4))
              tmp( 3,3)= expval*(              +xyz5( 5))
              tmp( 4,3)= expval*(xyz3( 2)*two  +xyz5( 7))
              tmp( 5,3)= expval*(xyz3( 3)      +xyz5( 8))
              tmp( 6,3)= expval*(              +xyz5( 9))
              tmp( 7,3)= expval*(xyz3( 4)*three+xyz5(11))
              tmp( 8,3)= expval*(xyz3( 5)*two  +xyz5(12))
              tmp( 9,3)= expval*(xyz3( 6)      +xyz5(13))
              tmp(10,3)= expval*(              +xyz5(14))
              tmp(11,3)= expval*(xyz3( 7)*four +xyz5(16))
              tmp(12,3)= expval*(xyz3( 8)*three+xyz5(17))
              tmp(13,3)= expval*(xyz3( 9)*two  +xyz5(18))
              tmp(14,3)= expval*(xyz3(10)      +xyz5(19))
              tmp(15,3)= expval*(              +xyz5(20))
              tmp( 1,4)= expval*(              +xyz5( 3))
              tmp( 2,4)= expval*(              +xyz5( 5))
              tmp( 3,4)= expval*(xyz3( 1)      +xyz5( 6))
              tmp( 4,4)= expval*(              +xyz5( 8))
              tmp( 5,4)= expval*(xyz3( 2)      +xyz5( 9))
              tmp( 6,4)= expval*(xyz3( 3)*two  +xyz5(10))
              tmp( 7,4)= expval*(              +xyz5(12))
              tmp( 8,4)= expval*(xyz3( 4)      +xyz5(13))
              tmp( 9,4)= expval*(xyz3( 5)*two  +xyz5(14))
              tmp(10,4)= expval*(xyz3( 6)*three+xyz5(15))
              tmp(11,4)= expval*(              +xyz5(17))
              tmp(12,4)= expval*(xyz3( 7)      +xyz5(18))
              tmp(13,4)= expval*(xyz3( 8)*two  +xyz5(19))
              tmp(14,4)= expval*(xyz3( 9)*three+xyz5(20))
              tmp(15,4)= expval*(xyz3(10)*four +xyz5(21))
              if(nbf == 15) then
                do jj= 1,4
                  vao(ilocbf+ 1,jj)= vao(ilocbf+ 1,jj)+tmp( 1,jj)
                  vao(ilocbf+ 2,jj)= vao(ilocbf+ 2,jj)+tmp( 2,jj)*sqrt7
                  vao(ilocbf+ 3,jj)= vao(ilocbf+ 3,jj)+tmp( 3,jj)*sqrt7
                  vao(ilocbf+ 4,jj)= vao(ilocbf+ 4,jj)+tmp( 4,jj)*sqrt35third
                  vao(ilocbf+ 5,jj)= vao(ilocbf+ 5,jj)+tmp( 5,jj)*sqrt35
                  vao(ilocbf+ 6,jj)= vao(ilocbf+ 6,jj)+tmp( 6,jj)*sqrt35third
                  vao(ilocbf+ 7,jj)= vao(ilocbf+ 7,jj)+tmp( 7,jj)*sqrt7
                  vao(ilocbf+ 8,jj)= vao(ilocbf+ 8,jj)+tmp( 8,jj)*sqrt35
                  vao(ilocbf+ 9,jj)= vao(ilocbf+ 9,jj)+tmp( 9,jj)*sqrt35
                  vao(ilocbf+10,jj)= vao(ilocbf+10,jj)+tmp(10,jj)*sqrt7
                  vao(ilocbf+11,jj)= vao(ilocbf+11,jj)+tmp(11,jj)
                  vao(ilocbf+12,jj)= vao(ilocbf+12,jj)+tmp(12,jj)*sqrt7
                  vao(ilocbf+13,jj)= vao(ilocbf+13,jj)+tmp(13,jj)*sqrt35third
                  vao(ilocbf+14,jj)= vao(ilocbf+14,jj)+tmp(14,jj)*sqrt7
                  vao(ilocbf+15,jj)= vao(ilocbf+15,jj)+tmp(15,jj)
                enddo
              else
                do jj= 1,4
                  vao(ilocbf+1,jj)= vao(ilocbf+1,jj)+(tmp(2,jj)-tmp(7,jj))*facg1
                  vao(ilocbf+2,jj)= vao(ilocbf+2,jj)+(-tmp(12,jj)+tmp(5,jj)*three)*facg2
                  vao(ilocbf+3,jj)= vao(ilocbf+3,jj)+(-tmp(2,jj)-tmp(7,jj)+tmp(9,jj)*six)*facg3
                  vao(ilocbf+4,jj)= vao(ilocbf+4,jj)+(-tmp(12,jj)*three+tmp(14,jj)*four &
&                                                     -tmp(5,jj)*three)*facg4
                  vao(ilocbf+5,jj)= vao(ilocbf+5,jj)+(tmp(1,jj)*three+tmp(11,jj)*three &
&                                                    +tmp(15,jj)*eight+tmp(4,jj)*six &
&                                                    -tmp(6,jj)*p24-tmp(13,jj)*p24)*eighth
                  vao(ilocbf+6,jj)= vao(ilocbf+6,jj)+(-tmp(3,jj)*three+tmp(10,jj)*four &
&                                                     -tmp(8,jj)*three)*facg4
                  vao(ilocbf+7,jj)= vao(ilocbf+7,jj)+(-tmp(1,jj)+tmp(11,jj)+tmp(6,jj)*six &
&                                                     -tmp(13,jj)*six)*facg5
                  vao(ilocbf+8,jj)= vao(ilocbf+8,jj)+(tmp(3,jj)-tmp(8,jj)*three)*facg2
                  vao(ilocbf+9,jj)= vao(ilocbf+9,jj)+(tmp(1,jj)+tmp(11,jj)-tmp(4,jj)*six)*facg6
                enddo
              endif
            enddo
!
! H function
!
          case(5)
            do iprim= 1,numprim
              icount= icount+1
              if(ex(icount)*rsqrd(iatom) > threshex) cycle
              expval= exp(-ex(icount)*rsqrd(iatom))*coeff(icount)
              fac= ex(icount)*two
              xx= xyzpt(1,iatom)*xyzpt(1,iatom)
              xy= xyzpt(1,iatom)*xyzpt(2,iatom)
              xz= xyzpt(1,iatom)*xyzpt(3,iatom)
              yy= xyzpt(2,iatom)*xyzpt(2,iatom)
              yz= xyzpt(2,iatom)*xyzpt(3,iatom)
              zz= xyzpt(3,iatom)*xyzpt(3,iatom)
              xyz4( 1)= xx*xx
              xyz4( 2)= xx*xy
              xyz4( 3)= xx*xz
              xyz4( 4)= xx*yy
              xyz4( 5)= xx*yz
              xyz4( 6)= xx*zz
              xyz4( 7)= xy*yy
              xyz4( 8)= xy*yz
              xyz4( 9)= xy*zz
              xyz4(10)= xz*zz
              xyz4(11)= yy*yy
              xyz4(12)= yy*yz
              xyz4(13)= yy*zz
              xyz4(14)= yz*zz
              xyz4(15)= zz*zz
              xyz6( 1)=-xyz4( 1)*xx*fac
              xyz6( 2)=-xyz4( 1)*xy*fac
              xyz6( 3)=-xyz4( 1)*xz*fac
              xyz6( 4)=-xyz4( 1)*yy*fac
              xyz6( 5)=-xyz4( 1)*yz*fac
              xyz6( 6)=-xyz4( 1)*zz*fac
              xyz6( 7)=-xyz4( 2)*yy*fac
              xyz6( 8)=-xyz4( 2)*yz*fac
              xyz6( 9)=-xyz4( 2)*zz*fac
              xyz6(10)=-xyz4( 3)*zz*fac
              xyz6(11)=-xyz4( 4)*yy*fac
              xyz6(12)=-xyz4( 4)*yz*fac
              xyz6(13)=-xyz4( 4)*zz*fac
              xyz6(14)=-xyz4( 5)*zz*fac
              xyz6(15)=-xyz4( 6)*zz*fac
              xyz6(16)=-xyz4( 7)*yy*fac
              xyz6(17)=-xyz4( 7)*yz*fac
              xyz6(18)=-xyz4( 7)*zz*fac
              xyz6(19)=-xyz4( 8)*zz*fac
              xyz6(20)=-xyz4( 9)*zz*fac
              xyz6(21)=-xyz4(10)*zz*fac
              xyz6(22)=-xyz4(11)*yy*fac
              xyz6(23)=-xyz4(11)*yz*fac
              xyz6(24)=-xyz4(11)*zz*fac
              xyz6(25)=-xyz4(12)*zz*fac
              xyz6(26)=-xyz4(13)*zz*fac
              xyz6(27)=-xyz4(14)*zz*fac
              xyz6(28)=-xyz4(15)*zz*fac
              tmp( 1,1)= expval*xyz4( 1)*xyzpt(1,iatom)
              tmp( 2,1)= expval*xyz4( 1)*xyzpt(2,iatom)
              tmp( 3,1)= expval*xyz4( 1)*xyzpt(3,iatom)
              tmp( 4,1)= expval*xyz4( 2)*xyzpt(2,iatom)
              tmp( 5,1)= expval*xyz4( 2)*xyzpt(3,iatom)
              tmp( 6,1)= expval*xyz4( 3)*xyzpt(3,iatom)
              tmp( 7,1)= expval*xyz4( 4)*xyzpt(2,iatom)
              tmp( 8,1)= expval*xyz4( 4)*xyzpt(3,iatom)
              tmp( 9,1)= expval*xyz4( 5)*xyzpt(3,iatom)
              tmp(10,1)= expval*xyz4( 6)*xyzpt(3,iatom)
              tmp(11,1)= expval*xyz4( 7)*xyzpt(2,iatom)
              tmp(12,1)= expval*xyz4( 7)*xyzpt(3,iatom)
              tmp(13,1)= expval*xyz4( 8)*xyzpt(3,iatom)
              tmp(14,1)= expval*xyz4( 9)*xyzpt(3,iatom)
              tmp(15,1)= expval*xyz4(10)*xyzpt(3,iatom)
              tmp(16,1)= expval*xyz4(11)*xyzpt(2,iatom)
              tmp(17,1)= expval*xyz4(11)*xyzpt(3,iatom)
              tmp(18,1)= expval*xyz4(12)*xyzpt(3,iatom)
              tmp(19,1)= expval*xyz4(13)*xyzpt(3,iatom)
              tmp(20,1)= expval*xyz4(14)*xyzpt(3,iatom)
              tmp(21,1)= expval*xyz4(15)*xyzpt(3,iatom)
              tmp( 1,2)= expval*(xyz4( 1)*five +xyz6( 1))
              tmp( 2,2)= expval*(xyz4( 2)*four +xyz6( 2))
              tmp( 3,2)= expval*(xyz4( 3)*four +xyz6( 3))
              tmp( 4,2)= expval*(xyz4( 4)*three+xyz6( 4))
              tmp( 5,2)= expval*(xyz4( 5)*three+xyz6( 5))
              tmp( 6,2)= expval*(xyz4( 6)*three+xyz6( 6))
              tmp( 7,2)= expval*(xyz4( 7)*two  +xyz6( 7))
              tmp( 8,2)= expval*(xyz4( 8)*two  +xyz6( 8))
              tmp( 9,2)= expval*(xyz4( 9)*two  +xyz6( 9))
              tmp(10,2)= expval*(xyz4(10)*two  +xyz6(10))
              tmp(11,2)= expval*(xyz4(11)      +xyz6(11))
              tmp(12,2)= expval*(xyz4(12)      +xyz6(12))
              tmp(13,2)= expval*(xyz4(13)      +xyz6(13))
              tmp(14,2)= expval*(xyz4(14)      +xyz6(14))
              tmp(15,2)= expval*(xyz4(15)      +xyz6(15))
              tmp(16,2)= expval*(              +xyz6(16))
              tmp(17,2)= expval*(              +xyz6(17))
              tmp(18,2)= expval*(              +xyz6(18))
              tmp(19,2)= expval*(              +xyz6(19))
              tmp(20,2)= expval*(              +xyz6(20))
              tmp(21,2)= expval*(              +xyz6(21))
              tmp( 1,3)= expval*(              +xyz6( 2))
              tmp( 2,3)= expval*(xyz4( 1)      +xyz6( 4))
              tmp( 3,3)= expval*(              +xyz6( 5))
              tmp( 4,3)= expval*(xyz4( 2)*two  +xyz6( 7))
              tmp( 5,3)= expval*(xyz4( 3)      +xyz6( 8))
              tmp( 6,3)= expval*(              +xyz6( 9))
              tmp( 7,3)= expval*(xyz4( 4)*three+xyz6(11))
              tmp( 8,3)= expval*(xyz4( 5)*two  +xyz6(12))
              tmp( 9,3)= expval*(xyz4( 6)      +xyz6(13))
              tmp(10,3)= expval*(              +xyz6(14))
              tmp(11,3)= expval*(xyz4( 7)*four +xyz6(16))
              tmp(12,3)= expval*(xyz4( 8)*three+xyz6(17))
              tmp(13,3)= expval*(xyz4( 9)*two  +xyz6(18))
              tmp(14,3)= expval*(xyz4(10)      +xyz6(19))
              tmp(15,3)= expval*(              +xyz6(20))
              tmp(16,3)= expval*(xyz4(11)*five +xyz6(22))
              tmp(17,3)= expval*(xyz4(12)*four +xyz6(23))
              tmp(18,3)= expval*(xyz4(13)*three+xyz6(24))
              tmp(19,3)= expval*(xyz4(14)*two  +xyz6(25))
              tmp(20,3)= expval*(xyz4(15)      +xyz6(26))
              tmp(21,3)= expval*(              +xyz6(27))
              tmp( 1,4)= expval*(              +xyz6( 3))
              tmp( 2,4)= expval*(              +xyz6( 5))
              tmp( 3,4)= expval*(xyz4( 1)      +xyz6( 6))
              tmp( 4,4)= expval*(              +xyz6( 8))
              tmp( 5,4)= expval*(xyz4( 2)      +xyz6( 9))
              tmp( 6,4)= expval*(xyz4( 3)*two  +xyz6(10))
              tmp( 7,4)= expval*(              +xyz6(12))
              tmp( 8,4)= expval*(xyz4( 4)      +xyz6(13))
              tmp( 9,4)= expval*(xyz4( 5)*two  +xyz6(14))
              tmp(10,4)= expval*(xyz4( 6)*three+xyz6(15))
              tmp(11,4)= expval*(              +xyz6(17))
              tmp(12,4)= expval*(xyz4( 7)      +xyz6(18))
              tmp(13,4)= expval*(xyz4( 8)*two  +xyz6(19))
              tmp(14,4)= expval*(xyz4( 9)*three+xyz6(20))
              tmp(15,4)= expval*(xyz4(10)*four +xyz6(21))
              tmp(16,4)= expval*(              +xyz6(23))
              tmp(17,4)= expval*(xyz4(11)      +xyz6(24))
              tmp(18,4)= expval*(xyz4(12)*two  +xyz6(25))
              tmp(19,4)= expval*(xyz4(13)*three+xyz6(26))
              tmp(20,4)= expval*(xyz4(14)*four +xyz6(27))
              tmp(21,4)= expval*(xyz4(15)*five +xyz6(28))
              if(nbf == 21) then
                do jj= 1,4
                  vao(ilocbf+ 1,jj)= vao(ilocbf+ 1,jj)+tmp( 1,jj)
                  vao(ilocbf+ 2,jj)= vao(ilocbf+ 2,jj)+tmp( 2,jj)*three
                  vao(ilocbf+ 3,jj)= vao(ilocbf+ 3,jj)+tmp( 3,jj)*three
                  vao(ilocbf+ 4,jj)= vao(ilocbf+ 4,jj)+tmp( 4,jj)*sqrt21
                  vao(ilocbf+ 5,jj)= vao(ilocbf+ 5,jj)+tmp( 5,jj)*sqrt63
                  vao(ilocbf+ 6,jj)= vao(ilocbf+ 6,jj)+tmp( 6,jj)*sqrt21
                  vao(ilocbf+ 7,jj)= vao(ilocbf+ 7,jj)+tmp( 7,jj)*sqrt21
                  vao(ilocbf+ 8,jj)= vao(ilocbf+ 8,jj)+tmp( 8,jj)*sqrt105
                  vao(ilocbf+ 9,jj)= vao(ilocbf+ 9,jj)+tmp( 9,jj)*sqrt105
                  vao(ilocbf+10,jj)= vao(ilocbf+10,jj)+tmp(10,jj)*sqrt21
                  vao(ilocbf+11,jj)= vao(ilocbf+11,jj)+tmp(11,jj)*three
                  vao(ilocbf+12,jj)= vao(ilocbf+12,jj)+tmp(12,jj)*sqrt63
                  vao(ilocbf+13,jj)= vao(ilocbf+13,jj)+tmp(13,jj)*sqrt105
                  vao(ilocbf+14,jj)= vao(ilocbf+14,jj)+tmp(14,jj)*sqrt63
                  vao(ilocbf+15,jj)= vao(ilocbf+15,jj)+tmp(15,jj)*three
                  vao(ilocbf+16,jj)= vao(ilocbf+16,jj)+tmp(16,jj)
                  vao(ilocbf+17,jj)= vao(ilocbf+17,jj)+tmp(17,jj)*three
                  vao(ilocbf+18,jj)= vao(ilocbf+18,jj)+tmp(18,jj)*sqrt21
                  vao(ilocbf+19,jj)= vao(ilocbf+19,jj)+tmp(19,jj)*sqrt21
                  vao(ilocbf+20,jj)= vao(ilocbf+20,jj)+tmp(20,jj)*three
                  vao(ilocbf+21,jj)= vao(ilocbf+21,jj)+tmp(21,jj)
                enddo
              else
                do jj= 1,4
                  vao(ilocbf+ 1,jj)= vao(ilocbf+ 1,jj)+(tmp(2,jj)*five-tmp(7,jj)*ten &
&                                                      +tmp(16,jj))*fach1
                  vao(ilocbf+ 2,jj)= vao(ilocbf+ 2,jj)+(tmp(5,jj)*four-tmp(12,jj)*four)*fach2
                  vao(ilocbf+ 3,jj)= vao(ilocbf+ 3,jj)+(-tmp(2,jj)*three-tmp(7,jj)*two &
&                                                      +tmp(9,jj)*p24+tmp(16,jj) &
&                                                      -tmp(18,jj)*eight)*fach3
                  vao(ilocbf+ 4,jj)= vao(ilocbf+ 4,jj)+(-tmp(5,jj)*two-tmp(12,jj)*two &
&                                                      +tmp(14,jj)*four)*fach4
                  vao(ilocbf+ 5,jj)= vao(ilocbf+ 5,jj)+(tmp(2,jj)+tmp(7,jj)*two-tmp(9,jj)*twelve &
&                                                      +tmp(16,jj)-tmp(18,jj)*twelve &
&                                                      +tmp(20,jj)*eight)*fach5
                  vao(ilocbf+ 6,jj)= vao(ilocbf+ 6,jj)+(tmp(3,jj)*p15+tmp(8,jj)*p30 &
&                                                      -tmp(10,jj)*p40+tmp(17,jj)*p15 &
&                                                      -tmp(19,jj)*p40 +tmp(21,jj)*eight)*eighth
                  vao(ilocbf+ 7,jj)= vao(ilocbf+ 7,jj)+(tmp(1,jj)+tmp(4,jj)*two-tmp(6,jj)*twelve &
&                                                      +tmp(11,jj)-tmp(13,jj)*twelve &
&                                                      +tmp(15,jj)*eight)*fach5
                  vao(ilocbf+ 8,jj)= vao(ilocbf+ 8,jj)+(-tmp(3,jj)+tmp(10,jj)*two+tmp(17,jj) &
&                                                      -tmp(19,jj)*two)*fach4
                  vao(ilocbf+ 9,jj)= vao(ilocbf+ 9,jj)+(-tmp(1,jj)+tmp(4,jj)*two+tmp(6,jj)*eight &
&                                                      +tmp(11,jj)*three-tmp(13,jj)*p24)*fach3
                  vao(ilocbf+10,jj)= vao(ilocbf+10,jj)+(tmp(3,jj)-tmp(8,jj)*six+tmp(17,jj))*fach2
                  vao(ilocbf+11,jj)= vao(ilocbf+11,jj)+(tmp(1,jj)-tmp(4,jj)*ten &
&                                                      +tmp(11,jj)*five)*fach1
                enddo
              endif
            enddo
!
! I function
!
          case(6)
            do iprim= 1,numprim
              icount= icount+1
              if(ex(icount)*rsqrd(iatom) > threshex) cycle
              expval= exp(-ex(icount)*rsqrd(iatom))*coeff(icount)
              fac= ex(icount)*two
              xx= xyzpt(1,iatom)*xyzpt(1,iatom)
              xy= xyzpt(1,iatom)*xyzpt(2,iatom)
              xz= xyzpt(1,iatom)*xyzpt(3,iatom)
              yy= xyzpt(2,iatom)*xyzpt(2,iatom)
              yz= xyzpt(2,iatom)*xyzpt(3,iatom)
              zz= xyzpt(3,iatom)*xyzpt(3,iatom)
              xyz3( 1)= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(1,iatom)
              xyz3( 2)= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(2,iatom)
              xyz3( 3)= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(3,iatom)
              xyz3( 4)= xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)
              xyz3( 5)= xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)
              xyz3( 6)= xyzpt(1,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)
              xyz3( 7)= xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)
              xyz3( 8)= xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)
              xyz3( 9)= xyzpt(2,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)
              xyz3(10)= xyzpt(3,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)
              xyz5( 1)= xyz3( 1)*xx
              xyz5( 2)= xyz3( 1)*xy
              xyz5( 3)= xyz3( 1)*xz
              xyz5( 4)= xyz3( 1)*yy
              xyz5( 5)= xyz3( 1)*yz
              xyz5( 6)= xyz3( 1)*zz
              xyz5( 7)= xyz3( 2)*yy
              xyz5( 8)= xyz3( 2)*yz
              xyz5( 9)= xyz3( 2)*zz
              xyz5(10)= xyz3( 3)*zz
              xyz5(11)= xyz3( 4)*yy
              xyz5(12)= xyz3( 4)*yz
              xyz5(13)= xyz3( 4)*zz
              xyz5(14)= xyz3( 5)*zz
              xyz5(15)= xyz3( 6)*zz
              xyz5(16)= xyz3( 7)*yy
              xyz5(17)= xyz3( 7)*yz
              xyz5(18)= xyz3( 7)*zz
              xyz5(19)= xyz3( 8)*zz
              xyz5(20)= xyz3( 9)*zz
              xyz5(21)= xyz3(10)*zz
              xyz7( 1)=-xyz5( 1)*xx*fac
              xyz7( 2)=-xyz5( 1)*xy*fac
              xyz7( 3)=-xyz5( 1)*xz*fac
              xyz7( 4)=-xyz5( 1)*yy*fac
              xyz7( 5)=-xyz5( 1)*yz*fac
              xyz7( 6)=-xyz5( 1)*zz*fac
              xyz7( 7)=-xyz5( 2)*yy*fac
              xyz7( 8)=-xyz5( 2)*yz*fac
              xyz7( 9)=-xyz5( 2)*zz*fac
              xyz7(10)=-xyz5( 3)*zz*fac
              xyz7(11)=-xyz5( 4)*yy*fac
              xyz7(12)=-xyz5( 4)*yz*fac
              xyz7(13)=-xyz5( 4)*zz*fac
              xyz7(14)=-xyz5( 6)*yz*fac
              xyz7(15)=-xyz5( 6)*zz*fac
              xyz7(16)=-xyz5( 7)*yy*fac
              xyz7(17)=-xyz5( 7)*yz*fac
              xyz7(18)=-xyz5( 7)*zz*fac
              xyz7(19)=-xyz5(10)*yy*fac
              xyz7(20)=-xyz5(10)*yz*fac
              xyz7(21)=-xyz5(10)*zz*fac
              xyz7(22)=-xyz5(11)*yy*fac
              xyz7(23)=-xyz5(11)*yz*fac
              xyz7(24)=-xyz5(11)*zz*fac
              xyz7(25)=-xyz5(14)*yy*fac
              xyz7(26)=-xyz5(14)*yz*fac
              xyz7(27)=-xyz5(14)*zz*fac
              xyz7(28)=-xyz5(15)*zz*fac
              xyz7(29)=-xyz5(16)*yy*fac
              xyz7(30)=-xyz5(16)*yz*fac
              xyz7(31)=-xyz5(16)*zz*fac
              xyz7(32)=-xyz5(19)*yy*fac
              xyz7(33)=-xyz5(19)*yz*fac
              xyz7(34)=-xyz5(19)*zz*fac
              xyz7(35)=-xyz5(21)*yz*fac
              xyz7(36)=-xyz5(21)*zz*fac
              tmp( 1,1)= expval*xyz5( 1)*xyzpt(1,iatom)
              tmp( 2,1)= expval*xyz5( 1)*xyzpt(2,iatom)
              tmp( 3,1)= expval*xyz5( 1)*xyzpt(3,iatom)
              tmp( 4,1)= expval*xyz5( 2)*xyzpt(2,iatom)
              tmp( 5,1)= expval*xyz5( 2)*xyzpt(3,iatom)
              tmp( 6,1)= expval*xyz5( 3)*xyzpt(3,iatom)
              tmp( 7,1)= expval*xyz5( 4)*xyzpt(2,iatom)
              tmp( 8,1)= expval*xyz5( 4)*xyzpt(3,iatom)
              tmp( 9,1)= expval*xyz5( 5)*xyzpt(3,iatom)
              tmp(10,1)= expval*xyz5( 6)*xyzpt(3,iatom)
              tmp(11,1)= expval*xyz5( 7)*xyzpt(2,iatom)
              tmp(12,1)= expval*xyz5( 7)*xyzpt(3,iatom)
              tmp(13,1)= expval*xyz5( 8)*xyzpt(3,iatom)
              tmp(14,1)= expval*xyz5( 9)*xyzpt(3,iatom)
              tmp(15,1)= expval*xyz5(10)*xyzpt(3,iatom)
              tmp(16,1)= expval*xyz5(11)*xyzpt(2,iatom)
              tmp(17,1)= expval*xyz5(11)*xyzpt(3,iatom)
              tmp(18,1)= expval*xyz5(12)*xyzpt(3,iatom)
              tmp(19,1)= expval*xyz5(13)*xyzpt(3,iatom)
              tmp(20,1)= expval*xyz5(14)*xyzpt(3,iatom)
              tmp(21,1)= expval*xyz5(15)*xyzpt(3,iatom)
              tmp(22,1)= expval*xyz5(16)*xyzpt(2,iatom)
              tmp(23,1)= expval*xyz5(16)*xyzpt(3,iatom)
              tmp(24,1)= expval*xyz5(17)*xyzpt(3,iatom)
              tmp(25,1)= expval*xyz5(18)*xyzpt(3,iatom)
              tmp(26,1)= expval*xyz5(19)*xyzpt(3,iatom)
              tmp(27,1)= expval*xyz5(20)*xyzpt(3,iatom)
              tmp(28,1)= expval*xyz5(21)*xyzpt(3,iatom)
              tmp( 1,2)= expval*(xyz5( 1)*six  +xyz7( 1))
              tmp( 2,2)= expval*(xyz5( 2)*five +xyz7( 2))
              tmp( 3,2)= expval*(xyz5( 3)*five +xyz7( 3))
              tmp( 4,2)= expval*(xyz5( 4)*four +xyz7( 4))
              tmp( 5,2)= expval*(xyz5( 5)*four +xyz7( 5))
              tmp( 6,2)= expval*(xyz5( 6)*four +xyz7( 6))
              tmp( 7,2)= expval*(xyz5( 7)*three+xyz7( 7))
              tmp( 8,2)= expval*(xyz5( 8)*three+xyz7( 8))
              tmp( 9,2)= expval*(xyz5( 9)*three+xyz7( 9))
              tmp(10,2)= expval*(xyz5(10)*three+xyz7(10))
              tmp(11,2)= expval*(xyz5(11)*two  +xyz7(11))
              tmp(12,2)= expval*(xyz5(12)*two  +xyz7(12))
              tmp(13,2)= expval*(xyz5(13)*two  +xyz7(13))
              tmp(14,2)= expval*(xyz5(14)*two  +xyz7(14))
              tmp(15,2)= expval*(xyz5(15)*two  +xyz7(15))
              tmp(16,2)= expval*(xyz5(16)      +xyz7(16))
              tmp(17,2)= expval*(xyz5(17)      +xyz7(17))
              tmp(18,2)= expval*(xyz5(18)      +xyz7(18))
              tmp(19,2)= expval*(xyz5(19)      +xyz7(19))
              tmp(20,2)= expval*(xyz5(20)      +xyz7(20))
              tmp(21,2)= expval*(xyz5(21)      +xyz7(21))
              tmp(22,2)= expval*(              +xyz7(22))
              tmp(23,2)= expval*(              +xyz7(23))
              tmp(24,2)= expval*(              +xyz7(24))
              tmp(25,2)= expval*(              +xyz7(25))
              tmp(26,2)= expval*(              +xyz7(26))
              tmp(27,2)= expval*(              +xyz7(27))
              tmp(28,2)= expval*(              +xyz7(28))
              tmp( 1,3)= expval*(              +xyz7( 2))
              tmp( 2,3)= expval*(xyz5( 1)      +xyz7( 4))
              tmp( 3,3)= expval*(              +xyz7( 5))
              tmp( 4,3)= expval*(xyz5( 2)*two  +xyz7( 7))
              tmp( 5,3)= expval*(xyz5( 3)      +xyz7( 8))
              tmp( 6,3)= expval*(              +xyz7( 9))
              tmp( 7,3)= expval*(xyz5( 4)*three+xyz7(11))
              tmp( 8,3)= expval*(xyz5( 5)*two  +xyz7(12))
              tmp( 9,3)= expval*(xyz5( 6)      +xyz7(13))
              tmp(10,3)= expval*(              +xyz7(14))
              tmp(11,3)= expval*(xyz5( 7)*four +xyz7(16))
              tmp(12,3)= expval*(xyz5( 8)*three+xyz7(17))
              tmp(13,3)= expval*(xyz5( 9)*two  +xyz7(18))
              tmp(14,3)= expval*(xyz5(10)      +xyz7(19))
              tmp(15,3)= expval*(              +xyz7(20))
              tmp(16,3)= expval*(xyz5(11)*five +xyz7(22))
              tmp(17,3)= expval*(xyz5(12)*four +xyz7(23))
              tmp(18,3)= expval*(xyz5(13)*three+xyz7(24))
              tmp(19,3)= expval*(xyz5(14)*two  +xyz7(25))
              tmp(20,3)= expval*(xyz5(15)      +xyz7(26))
              tmp(21,3)= expval*(              +xyz7(27))
              tmp(22,3)= expval*(xyz5(16)*six  +xyz7(29))
              tmp(23,3)= expval*(xyz5(17)*five +xyz7(30))
              tmp(24,3)= expval*(xyz5(18)*four +xyz7(31))
              tmp(25,3)= expval*(xyz5(19)*three+xyz7(32))
              tmp(26,3)= expval*(xyz5(20)*two  +xyz7(33))
              tmp(27,3)= expval*(xyz5(21)      +xyz7(34))
              tmp(28,3)= expval*(              +xyz7(35))
              tmp( 1,4)= expval*(              +xyz7( 3))
              tmp( 2,4)= expval*(              +xyz7( 5))
              tmp( 3,4)= expval*(xyz5( 1)      +xyz7( 6))
              tmp( 4,4)= expval*(              +xyz7( 8))
              tmp( 5,4)= expval*(xyz5( 2)      +xyz7( 9))
              tmp( 6,4)= expval*(xyz5( 3)*two  +xyz7(10))
              tmp( 7,4)= expval*(              +xyz7(12))
              tmp( 8,4)= expval*(xyz5( 4)      +xyz7(13))
              tmp( 9,4)= expval*(xyz5( 5)*two  +xyz7(14))
              tmp(10,4)= expval*(xyz5( 6)*three+xyz7(15))
              tmp(11,4)= expval*(              +xyz7(17))
              tmp(12,4)= expval*(xyz5( 7)      +xyz7(18))
              tmp(13,4)= expval*(xyz5( 8)*two  +xyz7(19))
              tmp(14,4)= expval*(xyz5( 9)*three+xyz7(20))
              tmp(15,4)= expval*(xyz5(10)*four +xyz7(21))
              tmp(16,4)= expval*(              +xyz7(23))
              tmp(17,4)= expval*(xyz5(11)      +xyz7(24))
              tmp(18,4)= expval*(xyz5(12)*two  +xyz7(25))
              tmp(19,4)= expval*(xyz5(13)*three+xyz7(26))
              tmp(20,4)= expval*(xyz5(14)*four +xyz7(27))
              tmp(21,4)= expval*(xyz5(15)*five +xyz7(28))
              tmp(22,4)= expval*(              +xyz7(30))
              tmp(23,4)= expval*(xyz5(16)      +xyz7(31))
              tmp(24,4)= expval*(xyz5(17)*two  +xyz7(32))
              tmp(25,4)= expval*(xyz5(18)*three+xyz7(33))
              tmp(26,4)= expval*(xyz5(19)*four +xyz7(34))
              tmp(27,4)= expval*(xyz5(20)*five +xyz7(35))
              tmp(28,4)= expval*(xyz5(21)*six  +xyz7(36))
              if(nbf == 28) then
                do jj= 1,4
                  vao(ilocbf+ 1,jj)= vao(ilocbf+ 1,jj)+tmp( 1,jj)
                  vao(ilocbf+ 2,jj)= vao(ilocbf+ 2,jj)+tmp( 2,jj)*sqrt11
                  vao(ilocbf+ 3,jj)= vao(ilocbf+ 3,jj)+tmp( 3,jj)*sqrt11
                  vao(ilocbf+ 4,jj)= vao(ilocbf+ 4,jj)+tmp( 4,jj)*sqrt33
                  vao(ilocbf+ 5,jj)= vao(ilocbf+ 5,jj)+tmp( 5,jj)*sqrt99
                  vao(ilocbf+ 6,jj)= vao(ilocbf+ 6,jj)+tmp( 6,jj)*sqrt33
                  vao(ilocbf+ 7,jj)= vao(ilocbf+ 7,jj)+tmp( 7,jj)*sqrt231fifth
                  vao(ilocbf+ 8,jj)= vao(ilocbf+ 8,jj)+tmp( 8,jj)*sqrt231
                  vao(ilocbf+ 9,jj)= vao(ilocbf+ 9,jj)+tmp( 9,jj)*sqrt231
                  vao(ilocbf+10,jj)= vao(ilocbf+10,jj)+tmp(10,jj)*sqrt231fifth
                  vao(ilocbf+11,jj)= vao(ilocbf+11,jj)+tmp(11,jj)*sqrt33
                  vao(ilocbf+12,jj)= vao(ilocbf+12,jj)+tmp(12,jj)*sqrt231
                  vao(ilocbf+13,jj)= vao(ilocbf+13,jj)+tmp(13,jj)*sqrt385
                  vao(ilocbf+14,jj)= vao(ilocbf+14,jj)+tmp(14,jj)*sqrt231
                  vao(ilocbf+15,jj)= vao(ilocbf+15,jj)+tmp(15,jj)*sqrt33
                  vao(ilocbf+16,jj)= vao(ilocbf+16,jj)+tmp(16,jj)*sqrt11
                  vao(ilocbf+17,jj)= vao(ilocbf+17,jj)+tmp(17,jj)*sqrt99
                  vao(ilocbf+18,jj)= vao(ilocbf+18,jj)+tmp(18,jj)*sqrt231
                  vao(ilocbf+19,jj)= vao(ilocbf+19,jj)+tmp(19,jj)*sqrt231
                  vao(ilocbf+20,jj)= vao(ilocbf+20,jj)+tmp(20,jj)*sqrt99
                  vao(ilocbf+21,jj)= vao(ilocbf+21,jj)+tmp(21,jj)*sqrt11
                  vao(ilocbf+22,jj)= vao(ilocbf+22,jj)+tmp(22,jj)
                  vao(ilocbf+23,jj)= vao(ilocbf+23,jj)+tmp(23,jj)*sqrt11
                  vao(ilocbf+24,jj)= vao(ilocbf+24,jj)+tmp(24,jj)*sqrt33
                  vao(ilocbf+25,jj)= vao(ilocbf+25,jj)+tmp(25,jj)*sqrt231fifth
                  vao(ilocbf+26,jj)= vao(ilocbf+26,jj)+tmp(26,jj)*sqrt33
                  vao(ilocbf+27,jj)= vao(ilocbf+27,jj)+tmp(27,jj)*sqrt11
                  vao(ilocbf+28,jj)= vao(ilocbf+28,jj)+tmp(28,jj)
                enddo
              else
                do jj= 1,4
                  vao(ilocbf+ 1,jj)= vao(ilocbf+ 1,jj) &
&                                  +(tmp(2,jj)*six-tmp(7,jj)*p20+tmp(16,jj)*six)*faci1
                  vao(ilocbf+ 2,jj)= vao(ilocbf+ 2,jj) &
&                                  +(tmp(5,jj)*five-tmp(12,jj)*ten+tmp(23,jj))*faci2
                  vao(ilocbf+ 3,jj)= vao(ilocbf+ 3,jj) &
&                                  +(-tmp(2,jj)*four+tmp(9,jj)*p40+tmp(16,jj)*four &
&                                   -tmp(18,jj)*p40)*faci3
                  vao(ilocbf+ 4,jj)= vao(ilocbf+ 4,jj) &
&                                  +(-tmp(5,jj)*p9-tmp(12,jj)*six+tmp(14,jj)*p24 &
&                                   +tmp(23,jj)*three-tmp(25,jj)*eight)*faci4
                  vao(ilocbf+ 5,jj)= vao(ilocbf+ 5,jj) &
&                                  +(tmp(2,jj)*two+tmp(7,jj)*four-tmp(9,jj)*p32 &
&                                   +tmp(16,jj)*two-tmp(18,jj)*p32+tmp(20,jj)*p32)*faci5
                  vao(ilocbf+ 6,jj)= vao(ilocbf+ 6,jj) &
&                                  +(tmp(5,jj)*five+tmp(12,jj)*ten-tmp(14,jj)*p20 &
&                                   +tmp(23,jj)*five-tmp(25,jj)*p20+tmp(27,jj)*eight)*faci6
                  vao(ilocbf+ 7,jj)= vao(ilocbf+ 7,jj) &
&                                  +(-tmp(1,jj)*five-tmp(4,jj)*p15+tmp(6,jj)*p90 &
&                                   -tmp(11,jj)*p15+tmp(13,jj)*p180-tmp(15,jj)*p120 &
&                                   -tmp(22,jj)*five+tmp(24,jj)*p90-tmp(26,jj)*p120 &
&                                   +tmp(28,jj)*p16)*sixteenth
                  vao(ilocbf+ 8,jj)= vao(ilocbf+ 8,jj) &
&                                  +(tmp(3,jj)*five+tmp(8,jj)*ten-tmp(10,jj)*p20 &
&                                   +tmp(17,jj)*five-tmp(19,jj)*p20+tmp(21,jj)*eight)*faci6
                  vao(ilocbf+ 9,jj)= vao(ilocbf+ 9,jj) &
&                                  +(tmp(1,jj)+tmp(4,jj)-tmp(6,jj)*p16-tmp(11,jj) &
&                                   +tmp(15,jj)*p16-tmp(22,jj)+tmp(24,jj)*p16 &
&                                   -tmp(26,jj)*p16)*faci5
                  vao(ilocbf+10,jj)= vao(ilocbf+10,jj) &
&                                  +(-tmp(3,jj)*three+tmp(8,jj)*six+tmp(10,jj)*eight &
&                                   +tmp(17,jj)*p9-tmp(19,jj)*p24)*faci4
                  vao(ilocbf+11,jj)= vao(ilocbf+11,jj) &
&                                  +(-tmp(1,jj)+tmp(4,jj)*five+tmp(6,jj)*ten+tmp(11,jj)*five &
&                                   -tmp(13,jj)*p60-tmp(22,jj)+tmp(24,jj)*ten)*faci3
                  vao(ilocbf+12,jj)= vao(ilocbf+12,jj) &
&                                  +(tmp(3,jj)-tmp(8,jj)*ten+tmp(17,jj)*five)*faci2
                  vao(ilocbf+13,jj)= vao(ilocbf+13,jj) &
&                                  +(tmp(1,jj)-tmp(4,jj)*p15+tmp(11,jj)*p15-tmp(22,jj))*faci1
                enddo
              endif
            enddo
          case default
            write(*,'(" Error! Subroutine Gridao supports up to i functions.")')
            call iabort
        end select
      enddo
!
      return
end


!-----------------------------------------
  subroutine gridd2ao(vg2ao,xyzpt,rsqrd)
!-----------------------------------------
!
! Calculate second derivatives of AO values for a grid point
!
      use modmolecule, only : natom
      use modbasis, only : ex, coeff, nshell, nao, locbf, locatom, mprim, mbf, mtype
      use modthresh, only : threshex
      implicit none
      integer :: icount, ish, numprim, iatom, iprim, nang, nbf, iloc, ii, jj
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, two=2.0D+00, three=3.0D+00
      real(8),parameter :: four=4.0D+00, five=5.0D+00, six=6.0D+00, seven=7.0D+00
      real(8),parameter :: eight=8.0D+00, p9=9.0D+00, ten=10.0D+00, eleven=11.0D+00
      real(8),parameter :: twelve=12.0D+00, p15=15.0D+00, p20=20.0D+00, p24=24.0D+00
      real(8),parameter :: p30=30.0D+00, p40=40.0D+00
      real(8),parameter :: eighth=0.125D+00
      real(8),parameter :: sqrt3=1.732050807568877D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrt5=2.236067977499790D+00, sqrt15=3.872983346207417D+00
      real(8),parameter :: sqrt7=2.645751311064590D+00, sqrt35=5.916079783099616D+00
      real(8),parameter :: sqrt35third=3.415650255319866D+00
      real(8),parameter :: sqrt21=4.582575694955840D+00, sqrt63=7.937253933193772D+00
      real(8),parameter :: sqrt105=1.024695076595960D+01, sqrt11=3.316624790355400D+00
      real(8),parameter :: facf1=0.79056941504209483D+00 ! sqrt(5/2)/2
      real(8),parameter :: facf2=3.87298334620741688D+00 ! sqrt(15)
      real(8),parameter :: facf3=0.61237243569579452D+00 ! sqrt(3/2)/2
      real(8),parameter :: facf4=1.93649167310370844D+00 ! sqrt(15)/2
      real(8),parameter :: facg1=2.95803989154980802D+00 ! sqrt(35)/2
      real(8),parameter :: facg2=2.09165006633518887D+00 ! sqrt(35/2)/2
      real(8),parameter :: facg3=1.11803398874989484D+00 ! sqrt(5)/2
      real(8),parameter :: facg4=0.79056941504209483D+00 ! sqrt(5/2)/2
      real(8),parameter :: facg5=0.55901699437494742D+00 ! sqrt(5)/4
      real(8),parameter :: facg6=0.73950997288745200D+00 ! sqrt(35)/8
      real(8),parameter :: fach1=0.70156076002011400D+00 ! sqrt(63/2)/8
      real(8),parameter :: fach2=2.21852991866235601D+00 ! sqrt(315)/8
      real(8),parameter :: fach3=0.52291251658379721D+00 ! sqrt(35/2)/8
      real(8),parameter :: fach4=2.56173769148989959D+00 ! sqrt(105)/4
      real(8),parameter :: fach5=0.48412291827592711D+00 ! sqrt(15)/8
      real(8),intent(in) :: xyzpt(3,natom), rsqrd(natom)
      real(8),intent(out) :: vg2ao(nao,6)
      real(8) :: fac, fac2, expval, tmp(21,6), xx, xy, xz, yy, yz, zz
      real(8) :: xyz3(10), xyz4(15), xyz5(21), xyz6(28), xyz7(36)
!
      vg2ao(:,:)= zero
      icount= 0
      do ish= 1,nshell
        nang= mtype(ish)
        numprim= mprim(ish)
        nbf = mbf(ish)
        iloc= locbf(ish)
        iatom= locatom(ish)
        select case(nang)
          case(0)
            do iprim= 1,numprim
              icount= icount+1
              if(ex(icount)*rsqrd(iatom) > threshex) cycle
              expval= exp(-ex(icount)*rsqrd(iatom))*coeff(icount)
              fac =-ex(icount)*two
              fac2= fac*fac
              xx= xyzpt(1,iatom)*xyzpt(1,iatom)
              xy= xyzpt(1,iatom)*xyzpt(2,iatom)
              xz= xyzpt(1,iatom)*xyzpt(3,iatom)
              yy= xyzpt(2,iatom)*xyzpt(2,iatom)
              yz= xyzpt(2,iatom)*xyzpt(3,iatom)
              zz= xyzpt(3,iatom)*xyzpt(3,iatom)
              vg2ao(iloc+1,1)= vg2ao(iloc+1,1)+expval*(fac+fac2*xx)
              vg2ao(iloc+1,2)= vg2ao(iloc+1,2)+expval*(    fac2*xy)
              vg2ao(iloc+1,3)= vg2ao(iloc+1,3)+expval*(    fac2*xz)
              vg2ao(iloc+1,4)= vg2ao(iloc+1,4)+expval*(fac+fac2*yy)
              vg2ao(iloc+1,5)= vg2ao(iloc+1,5)+expval*(    fac2*yz)
              vg2ao(iloc+1,6)= vg2ao(iloc+1,6)+expval*(fac+fac2*zz)
            enddo
          case(1)
            do iprim= 1,numprim
              icount= icount+1
              if(ex(icount)*rsqrd(iatom) > threshex) cycle
              expval= exp(-ex(icount)*rsqrd(iatom))*coeff(icount)
              fac =-ex(icount)*two
              fac2= fac*fac
              xyz3( 1)= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(1,iatom)*fac2
              xyz3( 2)= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(2,iatom)*fac2
              xyz3( 3)= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(3,iatom)*fac2
              xyz3( 4)= xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)*fac2
              xyz3( 5)= xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)*fac2
              xyz3( 6)= xyzpt(1,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac2
              xyz3( 7)= xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)*fac2
              xyz3( 8)= xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)*fac2
              xyz3( 9)= xyzpt(2,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac2
              xyz3(10)= xyzpt(3,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac2
              vg2ao(iloc+1,1)= vg2ao(iloc+1,1)+expval*(fac*xyzpt(1,iatom)*three+xyz3( 1))
              vg2ao(iloc+2,1)= vg2ao(iloc+2,1)+expval*(fac*xyzpt(2,iatom)      +xyz3( 2))
              vg2ao(iloc+3,1)= vg2ao(iloc+3,1)+expval*(fac*xyzpt(3,iatom)      +xyz3( 3))
              vg2ao(iloc+1,2)= vg2ao(iloc+1,2)+expval*(fac*xyzpt(2,iatom)      +xyz3( 2))
              vg2ao(iloc+2,2)= vg2ao(iloc+2,2)+expval*(fac*xyzpt(1,iatom)      +xyz3( 4))
              vg2ao(iloc+3,2)= vg2ao(iloc+3,2)+expval*(                        +xyz3( 5))
              vg2ao(iloc+1,3)= vg2ao(iloc+1,3)+expval*(fac*xyzpt(3,iatom)      +xyz3( 3))
              vg2ao(iloc+2,3)= vg2ao(iloc+2,3)+expval*(                        +xyz3( 5))
              vg2ao(iloc+3,3)= vg2ao(iloc+3,3)+expval*(fac*xyzpt(1,iatom)      +xyz3( 6))
              vg2ao(iloc+1,4)= vg2ao(iloc+1,4)+expval*(fac*xyzpt(1,iatom)      +xyz3( 4))
              vg2ao(iloc+2,4)= vg2ao(iloc+2,4)+expval*(fac*xyzpt(2,iatom)*three+xyz3( 7))
              vg2ao(iloc+3,4)= vg2ao(iloc+3,4)+expval*(fac*xyzpt(3,iatom)      +xyz3( 8))
              vg2ao(iloc+1,5)= vg2ao(iloc+1,5)+expval*(                        +xyz3( 5))
              vg2ao(iloc+2,5)= vg2ao(iloc+2,5)+expval*(fac*xyzpt(3,iatom)      +xyz3( 8))
              vg2ao(iloc+3,5)= vg2ao(iloc+3,5)+expval*(fac*xyzpt(2,iatom)      +xyz3( 9))
              vg2ao(iloc+1,6)= vg2ao(iloc+1,6)+expval*(fac*xyzpt(1,iatom)      +xyz3( 6))
              vg2ao(iloc+2,6)= vg2ao(iloc+2,6)+expval*(fac*xyzpt(2,iatom)      +xyz3( 9))
              vg2ao(iloc+3,6)= vg2ao(iloc+3,6)+expval*(fac*xyzpt(3,iatom)*three+xyz3(10))
            enddo
          case(2)
            do iprim= 1,numprim
              icount= icount+1
              if(ex(icount)*rsqrd(iatom) > threshex) cycle
              expval= exp(-ex(icount)*rsqrd(iatom))*coeff(icount)
              fac=-ex(icount)*two
              xx= xyzpt(1,iatom)*xyzpt(1,iatom)*fac
              yy= xyzpt(2,iatom)*xyzpt(2,iatom)*fac
              zz= xyzpt(3,iatom)*xyzpt(3,iatom)*fac
              xy= xyzpt(1,iatom)*xyzpt(2,iatom)*fac
              xz= xyzpt(1,iatom)*xyzpt(3,iatom)*fac
              yz= xyzpt(2,iatom)*xyzpt(3,iatom)*fac
              xyz4( 1)= xx*xx
              xyz4( 2)= xx*xy
              xyz4( 3)= xx*xz
              xyz4( 4)= xx*yy
              xyz4( 5)= xx*yz
              xyz4( 6)= xx*zz
              xyz4( 7)= xy*yy
              xyz4( 8)= xy*yz
              xyz4( 9)= xy*zz
              xyz4(10)= xz*zz
              xyz4(11)= yy*yy
              xyz4(12)= yy*yz
              xyz4(13)= yy*zz
              xyz4(14)= yz*zz
              xyz4(15)= zz*zz
              tmp(1,1)= expval*(two+xx*five +xyz4( 1))
              tmp(2,1)= expval*(   +xy*three+xyz4( 2))*sqrt3
              tmp(3,1)= expval*(   +xz*three+xyz4( 3))*sqrt3
              tmp(4,1)= expval*(   +yy      +xyz4( 4))
              tmp(5,1)= expval*(   +yz      +xyz4( 5))*sqrt3
              tmp(6,1)= expval*(   +zz      +xyz4( 6))
              tmp(1,2)= expval*(   +xy*two  +xyz4( 2))
              tmp(2,2)= expval*(one+xx+yy   +xyz4( 4))*sqrt3
              tmp(3,2)= expval*(   +yz      +xyz4( 5))*sqrt3
              tmp(4,2)= expval*(   +xy*two  +xyz4( 7))
              tmp(5,2)= expval*(   +xz      +xyz4( 8))*sqrt3
              tmp(6,2)= expval*(            +xyz4( 9))
              tmp(1,3)= expval*(   +xz*two  +xyz4( 3))
              tmp(2,3)= expval*(   +yz      +xyz4( 5))*sqrt3
              tmp(3,3)= expval*(one+xx+zz   +xyz4( 6))*sqrt3
              tmp(4,3)= expval*(            +xyz4( 8))
              tmp(5,3)= expval*(   +xy      +xyz4( 9))*sqrt3
              tmp(6,3)= expval*(   +xz*two  +xyz4(10))
              tmp(1,4)= expval*(   +xx      +xyz4( 4))
              tmp(2,4)= expval*(   +xy*three+xyz4( 7))*sqrt3
              tmp(3,4)= expval*(   +xz      +xyz4( 8))*sqrt3
              tmp(4,4)= expval*(two+yy*five +xyz4(11))
              tmp(5,4)= expval*(   +yz*three+xyz4(12))*sqrt3
              tmp(6,4)= expval*(   +zz      +xyz4(13))
              tmp(1,5)= expval*(            +xyz4( 5))
              tmp(2,5)= expval*(   +xz      +xyz4( 8))*sqrt3
              tmp(3,5)= expval*(   +xy      +xyz4( 9))*sqrt3
              tmp(4,5)= expval*(   +yz*two  +xyz4(12))
              tmp(5,5)= expval*(one+yy+zz   +xyz4(13))*sqrt3
              tmp(6,5)= expval*(   +yz*two  +xyz4(14))
              tmp(1,6)= expval*(   +xx      +xyz4( 6))
              tmp(2,6)= expval*(   +xy      +xyz4( 9))*sqrt3
              tmp(3,6)= expval*(   +xz*three+xyz4(10))*sqrt3
              tmp(4,6)= expval*(   +yy      +xyz4(13))
              tmp(5,6)= expval*(   +yz*three+xyz4(14))*sqrt3
              tmp(6,6)= expval*(two+zz*five +xyz4(15))
              if(nbf == 6) then
                do ii= 1,6
                  do jj= 1,6
                    vg2ao(iloc+jj,ii)= vg2ao(iloc+jj,ii)+tmp(jj,ii)
                  enddo
                enddo
              else
                do ii= 1,6
                  vg2ao(iloc+1,ii)= vg2ao(iloc+1,ii)+tmp(2,ii)
                  vg2ao(iloc+2,ii)= vg2ao(iloc+2,ii)+tmp(5,ii)
                  vg2ao(iloc+3,ii)= vg2ao(iloc+3,ii)+tmp(6,ii)-(tmp(1,ii)+tmp(4,ii))*half
                  vg2ao(iloc+4,ii)= vg2ao(iloc+4,ii)+tmp(3,ii)
                  vg2ao(iloc+5,ii)= vg2ao(iloc+5,ii)+(tmp(1,ii)-tmp(4,ii))*sqrt3h
                enddo
              endif
            enddo
          case(3)
            do iprim= 1,numprim
              icount= icount+1
              if(ex(icount)*rsqrd(iatom) > threshex) cycle
              expval= exp(-ex(icount)*rsqrd(iatom))*coeff(icount)
              fac=-ex(icount)*two
              xx= xyzpt(1,iatom)*xyzpt(1,iatom)*fac
              yy= xyzpt(2,iatom)*xyzpt(2,iatom)*fac
              zz= xyzpt(3,iatom)*xyzpt(3,iatom)*fac
              xy= xyzpt(1,iatom)*xyzpt(2,iatom)*fac
              xz= xyzpt(1,iatom)*xyzpt(3,iatom)*fac
              yz= xyzpt(2,iatom)*xyzpt(3,iatom)*fac
              xyz3( 1)= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(1,iatom)*fac
              xyz3( 2)= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(2,iatom)*fac
              xyz3( 3)= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(3,iatom)*fac
              xyz3( 4)= xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)*fac
              xyz3( 5)= xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)*fac
              xyz3( 6)= xyzpt(1,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac
              xyz3( 7)= xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)*fac
              xyz3( 8)= xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)*fac
              xyz3( 9)= xyzpt(2,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac
              xyz3(10)= xyzpt(3,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac
              xyz5( 1)= xyz3( 1)*xx
              xyz5( 2)= xyz3( 1)*xy
              xyz5( 3)= xyz3( 1)*xz
              xyz5( 4)= xyz3( 1)*yy
              xyz5( 5)= xyz3( 1)*yz
              xyz5( 6)= xyz3( 1)*zz
              xyz5( 7)= xyz3( 2)*yy
              xyz5( 8)= xyz3( 2)*yz
              xyz5( 9)= xyz3( 2)*zz
              xyz5(10)= xyz3( 3)*zz
              xyz5(11)= xyz3( 4)*yy
              xyz5(12)= xyz3( 4)*yz
              xyz5(13)= xyz3( 4)*zz
              xyz5(14)= xyz3( 5)*zz
              xyz5(15)= xyz3( 6)*zz
              xyz5(16)= xyz3( 7)*yy
              xyz5(17)= xyz3( 7)*yz
              xyz5(18)= xyz3( 7)*zz
              xyz5(19)= xyz3( 8)*zz
              xyz5(20)= xyz3( 9)*zz
              xyz5(21)= xyz3(10)*zz
              tmp( 1,1)= expval*(xyzpt(1,iatom)*six+xyz3( 1)*seven               +xyz5( 1))
              tmp( 2,1)= expval*(xyzpt(2,iatom)*two+xyz3( 2)*five                +xyz5( 2))
              tmp( 3,1)= expval*(xyzpt(3,iatom)*two+xyz3( 3)*five                +xyz5( 3))
              tmp( 4,1)= expval*(                  +xyz3( 4)*three               +xyz5( 4))
              tmp( 5,1)= expval*(                  +xyz3( 5)*three               +xyz5( 5))
              tmp( 6,1)= expval*(                  +xyz3( 6)*three               +xyz5( 6))
              tmp( 7,1)= expval*(                  +xyz3( 7)                     +xyz5( 7))
              tmp( 8,1)= expval*(                  +xyz3( 8)                     +xyz5( 8))
              tmp( 9,1)= expval*(                  +xyz3( 9)                     +xyz5( 9))
              tmp(10,1)= expval*(                  +xyz3(10)                     +xyz5(10))
              tmp( 1,2)= expval*(                  +xyz3( 2)*three               +xyz5( 2))
              tmp( 2,2)= expval*(xyzpt(1,iatom)*two+xyz3( 1)      +xyz3( 4)*two  +xyz5( 4))
              tmp( 3,2)= expval*(                                 +xyz3( 5)*two  +xyz5( 5))
              tmp( 4,2)= expval*(xyzpt(2,iatom)*two+xyz3( 2)*two  +xyz3( 7)      +xyz5( 7))
              tmp( 5,2)= expval*(xyzpt(3,iatom)    +xyz3( 3)      +xyz3( 8)      +xyz5( 8))
              tmp( 6,2)= expval*(                                 +xyz3( 9)      +xyz5( 9))
              tmp( 7,2)= expval*(                  +xyz3( 4)*three               +xyz5(11))
              tmp( 8,2)= expval*(                  +xyz3( 5)*two                 +xyz5(12))
              tmp( 9,2)= expval*(                  +xyz3( 6)                     +xyz5(13))
              tmp(10,2)= expval*(                                                +xyz5(14))
              tmp( 1,3)= expval*(                                 +xyz3( 3)*three+xyz5( 3))
              tmp( 2,3)= expval*(                                 +xyz3( 5)*two  +xyz5( 5))
              tmp( 3,3)= expval*(xyzpt(1,iatom)*two+xyz3( 1)      +xyz3( 6)*two  +xyz5( 6))
              tmp( 4,3)= expval*(                                 +xyz3( 8)      +xyz5( 8))
              tmp( 5,3)= expval*(xyzpt(2,iatom)    +xyz3( 2)      +xyz3( 9)      +xyz5( 9))
              tmp( 6,3)= expval*(xyzpt(3,iatom)*two+xyz3( 3)*two  +xyz3(10)      +xyz5(10))
              tmp( 7,3)= expval*(                                                +xyz5(12))
              tmp( 8,3)= expval*(                  +xyz3( 4)                     +xyz5(13))
              tmp( 9,3)= expval*(                  +xyz3( 5)*two                 +xyz5(14))
              tmp(10,3)= expval*(                  +xyz3( 6)*three               +xyz5(15))
              tmp( 1,4)= expval*(                  +xyz3( 1)                     +xyz5( 4))
              tmp( 2,4)= expval*(                  +xyz3( 2)*three               +xyz5( 7))
              tmp( 3,4)= expval*(                  +xyz3( 3)                     +xyz5( 8))
              tmp( 4,4)= expval*(xyzpt(1,iatom)*two+xyz3( 4)*five                +xyz5(11))
              tmp( 5,4)= expval*(                  +xyz3( 5)*three               +xyz5(12))
              tmp( 6,4)= expval*(                  +xyz3( 6)                     +xyz5(13))
              tmp( 7,4)= expval*(xyzpt(2,iatom)*six+xyz3( 7)*seven               +xyz5(16))
              tmp( 8,4)= expval*(xyzpt(3,iatom)*two+xyz3( 8)*five                +xyz5(17))
              tmp( 9,4)= expval*(                  +xyz3( 9)*three               +xyz5(18))
              tmp(10,4)= expval*(                  +xyz3(10)                     +xyz5(19))
              tmp( 1,5)= expval*(                                                +xyz5( 5))
              tmp( 2,5)= expval*(                                 +xyz3( 3)      +xyz5( 8))
              tmp( 3,5)= expval*(                  +xyz3( 2)                     +xyz5( 9))
              tmp( 4,5)= expval*(                                 +xyz3( 5)*two  +xyz5(12))
              tmp( 5,5)= expval*(xyzpt(1,iatom)    +xyz3( 4)      +xyz3( 6)      +xyz5(13))
              tmp( 6,5)= expval*(                  +xyz3( 5)*two                 +xyz5(14))
              tmp( 7,5)= expval*(                                 +xyz3( 8)*three+xyz5(17))
              tmp( 8,5)= expval*(xyzpt(2,iatom)*two+xyz3( 7)      +xyz3( 9)*two  +xyz5(18))
              tmp( 9,5)= expval*(xyzpt(3,iatom)*two+xyz3( 8)*two  +xyz3(10)      +xyz5(19))
              tmp(10,5)= expval*(                  +xyz3( 9)*three               +xyz5(20))
              tmp( 1,6)= expval*(                  +xyz3( 1)                     +xyz5( 6))
              tmp( 2,6)= expval*(                  +xyz3( 2)                     +xyz5( 9))
              tmp( 3,6)= expval*(                  +xyz3( 3)*three               +xyz5(10))
              tmp( 4,6)= expval*(                  +xyz3( 4)                     +xyz5(13))
              tmp( 5,6)= expval*(                  +xyz3( 5)*three               +xyz5(14))
              tmp( 6,6)= expval*(xyzpt(1,iatom)*two+xyz3( 6)*five                +xyz5(15))
              tmp( 7,6)= expval*(                  +xyz3( 7)                     +xyz5(18))
              tmp( 8,6)= expval*(                  +xyz3( 8)*three               +xyz5(19))
              tmp( 9,6)= expval*(xyzpt(2,iatom)*two+xyz3( 9)*five                +xyz5(20))
              tmp(10,6)= expval*(xyzpt(3,iatom)*six+xyz3(10)*seven               +xyz5(21))
              if(nbf == 10) then
                do ii= 1,6
                  vg2ao(iloc+ 1,ii)= vg2ao(iloc+ 1,ii)+tmp( 1,ii)
                  vg2ao(iloc+ 2,ii)= vg2ao(iloc+ 2,ii)+tmp( 2,ii)*sqrt5
                  vg2ao(iloc+ 3,ii)= vg2ao(iloc+ 3,ii)+tmp( 3,ii)*sqrt5
                  vg2ao(iloc+ 4,ii)= vg2ao(iloc+ 4,ii)+tmp( 4,ii)*sqrt5
                  vg2ao(iloc+ 5,ii)= vg2ao(iloc+ 5,ii)+tmp( 5,ii)*sqrt15
                  vg2ao(iloc+ 6,ii)= vg2ao(iloc+ 6,ii)+tmp( 6,ii)*sqrt5
                  vg2ao(iloc+ 7,ii)= vg2ao(iloc+ 7,ii)+tmp( 7,ii)
                  vg2ao(iloc+ 8,ii)= vg2ao(iloc+ 8,ii)+tmp( 8,ii)*sqrt5
                  vg2ao(iloc+ 9,ii)= vg2ao(iloc+ 9,ii)+tmp( 9,ii)*sqrt5
                  vg2ao(iloc+10,ii)= vg2ao(iloc+10,ii)+tmp(10,ii)
                enddo
              else
                do ii= 1,6
                  vg2ao(iloc+1,ii)= vg2ao(iloc+1,ii)+(-tmp(7,ii)+three*tmp(2,ii))*facf1
                  vg2ao(iloc+2,ii)= vg2ao(iloc+2,ii)+tmp(5,ii)*facf2
                  vg2ao(iloc+3,ii)= vg2ao(iloc+3,ii)+(-tmp(7,ii)-tmp(2,ii)+four*tmp(9,ii))*facf3
                  vg2ao(iloc+4,ii)= vg2ao(iloc+4,ii)+(two*tmp(10,ii) &
&                                                     -three*(tmp(3,ii)+tmp(8,ii)))*half
                  vg2ao(iloc+5,ii)= vg2ao(iloc+5,ii)+(-tmp(1,ii)-tmp(4,ii)+four*tmp(6,ii))*facf3
                  vg2ao(iloc+6,ii)= vg2ao(iloc+6,ii)+(tmp(3,ii)-tmp(8,ii))*facf4
                  vg2ao(iloc+7,ii)= vg2ao(iloc+7,ii)+(tmp(1,ii)-three*tmp(4,ii))*facf1
                enddo
              endif
            enddo
          case(4)
            do iprim= 1,numprim
              icount= icount+1
              if(ex(icount)*rsqrd(iatom) > threshex) cycle
              expval= exp(-ex(icount)*rsqrd(iatom))*coeff(icount)
              fac=-ex(icount)*two
              xx= xyzpt(1,iatom)*xyzpt(1,iatom)
              yy= xyzpt(2,iatom)*xyzpt(2,iatom)
              zz= xyzpt(3,iatom)*xyzpt(3,iatom)
              xy= xyzpt(1,iatom)*xyzpt(2,iatom)
              xz= xyzpt(1,iatom)*xyzpt(3,iatom)
              yz= xyzpt(2,iatom)*xyzpt(3,iatom)
              xyz4( 1)= xx*xx*fac
              xyz4( 2)= xx*xy*fac
              xyz4( 3)= xx*xz*fac
              xyz4( 4)= xx*yy*fac
              xyz4( 5)= xx*yz*fac
              xyz4( 6)= xx*zz*fac
              xyz4( 7)= xy*yy*fac
              xyz4( 8)= xy*yz*fac
              xyz4( 9)= xy*zz*fac
              xyz4(10)= xz*zz*fac
              xyz4(11)= yy*yy*fac
              xyz4(12)= yy*yz*fac
              xyz4(13)= yy*zz*fac
              xyz4(14)= yz*zz*fac
              xyz4(15)= zz*zz*fac
              xyz6( 1)= xyz4( 1)*xx*fac
              xyz6( 2)= xyz4( 1)*xy*fac
              xyz6( 3)= xyz4( 1)*xz*fac
              xyz6( 4)= xyz4( 1)*yy*fac
              xyz6( 5)= xyz4( 1)*yz*fac
              xyz6( 6)= xyz4( 1)*zz*fac
              xyz6( 7)= xyz4( 2)*yy*fac
              xyz6( 8)= xyz4( 2)*yz*fac
              xyz6( 9)= xyz4( 2)*zz*fac
              xyz6(10)= xyz4( 3)*zz*fac
              xyz6(11)= xyz4( 4)*yy*fac
              xyz6(12)= xyz4( 4)*yz*fac
              xyz6(13)= xyz4( 4)*zz*fac
              xyz6(14)= xyz4( 5)*zz*fac
              xyz6(15)= xyz4( 6)*zz*fac
              xyz6(16)= xyz4( 7)*yy*fac
              xyz6(17)= xyz4( 7)*yz*fac
              xyz6(18)= xyz4( 7)*zz*fac
              xyz6(19)= xyz4( 8)*zz*fac
              xyz6(20)= xyz4( 9)*zz*fac
              xyz6(21)= xyz4(10)*zz*fac
              xyz6(22)= xyz4(11)*yy*fac
              xyz6(23)= xyz4(11)*yz*fac
              xyz6(24)= xyz4(11)*zz*fac
              xyz6(25)= xyz4(12)*zz*fac
              xyz6(26)= xyz4(13)*zz*fac
              xyz6(27)= xyz4(14)*zz*fac
              xyz6(28)= xyz4(15)*zz*fac
              tmp( 1,1)= expval*(xx*twelve+xyz4( 1)*p9                  +xyz6( 1))
              tmp( 2,1)= expval*(xy*six   +xyz4( 2)*seven               +xyz6( 2))
              tmp( 3,1)= expval*(xz*six   +xyz4( 3)*seven               +xyz6( 3))
              tmp( 4,1)= expval*(yy*two   +xyz4( 4)*five                +xyz6( 4))
              tmp( 5,1)= expval*(yz*two   +xyz4( 5)*five                +xyz6( 5))
              tmp( 6,1)= expval*(zz*two   +xyz4( 6)*five                +xyz6( 6))
              tmp( 7,1)= expval*(          xyz4( 7)*three               +xyz6( 7))
              tmp( 8,1)= expval*(          xyz4( 8)*three               +xyz6( 8))
              tmp( 9,1)= expval*(          xyz4( 9)*three               +xyz6( 9))
              tmp(10,1)= expval*(          xyz4(10)*three               +xyz6(10))
              tmp(11,1)= expval*(          xyz4(11)                     +xyz6(11))
              tmp(12,1)= expval*(          xyz4(12)                     +xyz6(12))
              tmp(13,1)= expval*(          xyz4(13)                     +xyz6(13))
              tmp(14,1)= expval*(          xyz4(14)                     +xyz6(14))
              tmp(15,1)= expval*(          xyz4(15)                     +xyz6(15))
              tmp( 1,2)= expval*(          xyz4( 2)*four                +xyz6( 2))
              tmp( 2,2)= expval*(xx*three +xyz4( 1)      +xyz4( 4)*three+xyz6( 4))
              tmp( 3,2)= expval*(                         xyz4( 5)*three+xyz6( 5))
              tmp( 4,2)= expval*(xy*four  +xyz4( 2)*two  +xyz4( 7)*two  +xyz6( 7))
              tmp( 5,2)= expval*(xz*two   +xyz4( 3)      +xyz4( 8)*two  +xyz6( 8))
              tmp( 6,2)= expval*(                         xyz4( 9)*two  +xyz6( 9))
              tmp( 7,2)= expval*(yy*three +xyz4( 4)*three+xyz4(11)      +xyz6(11))
              tmp( 8,2)= expval*(yz*two   +xyz4( 5)*two  +xyz4(12)      +xyz6(12))
              tmp( 9,2)= expval*(zz       +xyz4( 6)      +xyz4(13)      +xyz6(13))
              tmp(10,2)= expval*(                         xyz4(14)      +xyz6(14))
              tmp(11,2)= expval*(          xyz4( 7)*four                +xyz6(16))
              tmp(12,2)= expval*(          xyz4( 8)*three               +xyz6(17))
              tmp(13,2)= expval*(          xyz4( 9)*two                 +xyz6(18))
              tmp(14,2)= expval*(          xyz4(10)                     +xyz6(19))
              tmp(15,2)= expval*(                                        xyz6(20))
              tmp( 1,3)= expval*(                         xyz4( 3)*four +xyz6( 3))
              tmp( 2,3)= expval*(                         xyz4( 5)*three+xyz6( 5))
              tmp( 3,3)= expval*(xx*three +xyz4( 1)      +xyz4( 6)*three+xyz6( 6))
              tmp( 4,3)= expval*(                         xyz4( 8)*two  +xyz6( 8))
              tmp( 5,3)= expval*(xy*two   +xyz4( 2)      +xyz4( 9)*two  +xyz6( 9))
              tmp( 6,3)= expval*(xz*four  +xyz4( 3)*two  +xyz4(10)*two  +xyz6(10))
              tmp( 7,3)= expval*(                         xyz4(12)      +xyz6(12))
              tmp( 8,3)= expval*(yy       +xyz4( 4)      +xyz4(13)      +xyz6(13))
              tmp( 9,3)= expval*(yz*two   +xyz4( 5)*two  +xyz4(14)      +xyz6(14))
              tmp(10,3)= expval*(zz*three +xyz4( 6)*three+xyz4(15)      +xyz6(15))
              tmp(11,3)= expval*(                                        xyz6(17))
              tmp(12,3)= expval*(          xyz4( 7)                     +xyz6(18))
              tmp(13,3)= expval*(          xyz4( 8)*two                 +xyz6(19))
              tmp(14,3)= expval*(          xyz4( 9)*three               +xyz6(20))
              tmp(15,3)= expval*(          xyz4(10)*four                +xyz6(21))
              tmp( 1,4)= expval*(          xyz4( 1)                     +xyz6( 4))
              tmp( 2,4)= expval*(          xyz4( 2)*three               +xyz6( 7))
              tmp( 3,4)= expval*(          xyz4( 3)                     +xyz6( 8))
              tmp( 4,4)= expval*(xx*two   +xyz4( 4)*five                +xyz6(11))
              tmp( 5,4)= expval*(          xyz4( 5)*three               +xyz6(12))
              tmp( 6,4)= expval*(          xyz4( 6)                     +xyz6(13))
              tmp( 7,4)= expval*(xy*six   +xyz4( 7)*seven               +xyz6(16))
              tmp( 8,4)= expval*(xz*two   +xyz4( 8)*five                +xyz6(17))
              tmp( 9,4)= expval*(          xyz4( 9)*three               +xyz6(18))
              tmp(10,4)= expval*(          xyz4(10)                     +xyz6(19))
              tmp(11,4)= expval*(yy*twelve+xyz4(11)*p9                  +xyz6(22))
              tmp(12,4)= expval*(yz*six   +xyz4(12)*seven               +xyz6(23))
              tmp(13,4)= expval*(zz*two   +xyz4(13)*five                +xyz6(24))
              tmp(14,4)= expval*(          xyz4(14)*three               +xyz6(25))
              tmp(15,4)= expval*(          xyz4(15)                     +xyz6(26))
              tmp( 1,5)= expval*(                                        xyz6( 5))
              tmp( 2,5)= expval*(                         xyz4( 3)      +xyz6( 8))
              tmp( 3,5)= expval*(          xyz4( 2)                     +xyz6( 9))
              tmp( 4,5)= expval*(                         xyz4( 5)*two  +xyz6(12))
              tmp( 5,5)= expval*(xx       +xyz4( 4)      +xyz4( 6)      +xyz6(13))
              tmp( 6,5)= expval*(          xyz4( 5)*two                 +xyz6(14))
              tmp( 7,5)= expval*(                         xyz4( 8)*three+xyz6(17))
              tmp( 8,5)= expval*(xy*two   +xyz4( 7)      +xyz4( 9)*two  +xyz6(18))
              tmp( 9,5)= expval*(xz*two   +xyz4( 8)*two  +xyz4(10)      +xyz6(19))
              tmp(10,5)= expval*(          xyz4( 9)*three               +xyz6(20))
              tmp(11,5)= expval*(                         xyz4(12)*four +xyz6(23))
              tmp(12,5)= expval*(yy*three +xyz4(11)      +xyz4(13)*three+xyz6(24))
              tmp(13,5)= expval*(yz*four  +xyz4(12)*two  +xyz4(14)*two  +xyz6(25))
              tmp(14,5)= expval*(zz*three +xyz4(13)*three+xyz4(15)      +xyz6(26))
              tmp(15,5)= expval*(          xyz4(14)*four                +xyz6(27))
              tmp( 1,6)= expval*(          xyz4( 1)                     +xyz6( 6))
              tmp( 2,6)= expval*(          xyz4( 2)                     +xyz6( 9))
              tmp( 3,6)= expval*(          xyz4( 3)*three               +xyz6(10))
              tmp( 4,6)= expval*(          xyz4( 4)                     +xyz6(13))
              tmp( 5,6)= expval*(          xyz4( 5)*three               +xyz6(14))
              tmp( 6,6)= expval*(xx*two   +xyz4( 6)*five                +xyz6(15))
              tmp( 7,6)= expval*(          xyz4( 7)                     +xyz6(18))
              tmp( 8,6)= expval*(          xyz4( 8)*three               +xyz6(19))
              tmp( 9,6)= expval*(xy*two   +xyz4( 9)*five                +xyz6(20))
              tmp(10,6)= expval*(xz*six   +xyz4(10)*seven               +xyz6(21))
              tmp(11,6)= expval*(          xyz4(11)                     +xyz6(24))
              tmp(12,6)= expval*(          xyz4(12)*three               +xyz6(25))
              tmp(13,6)= expval*(yy*two   +xyz4(13)*five                +xyz6(26))
              tmp(14,6)= expval*(yz*six   +xyz4(14)*seven               +xyz6(27))
              tmp(15,6)= expval*(zz*twelve+xyz4(15)*p9                  +xyz6(28))
              if(nbf == 15) then
                do ii= 1,6
                  vg2ao(iloc+ 1,ii)= vg2ao(iloc+ 1,ii)+tmp( 1,ii)
                  vg2ao(iloc+ 2,ii)= vg2ao(iloc+ 2,ii)+tmp( 2,ii)*sqrt7
                  vg2ao(iloc+ 3,ii)= vg2ao(iloc+ 3,ii)+tmp( 3,ii)*sqrt7
                  vg2ao(iloc+ 4,ii)= vg2ao(iloc+ 4,ii)+tmp( 4,ii)*sqrt35third
                  vg2ao(iloc+ 5,ii)= vg2ao(iloc+ 5,ii)+tmp( 5,ii)*sqrt35
                  vg2ao(iloc+ 6,ii)= vg2ao(iloc+ 6,ii)+tmp( 6,ii)*sqrt35third
                  vg2ao(iloc+ 7,ii)= vg2ao(iloc+ 7,ii)+tmp( 7,ii)*sqrt7
                  vg2ao(iloc+ 8,ii)= vg2ao(iloc+ 8,ii)+tmp( 8,ii)*sqrt35
                  vg2ao(iloc+ 9,ii)= vg2ao(iloc+ 9,ii)+tmp( 9,ii)*sqrt35
                  vg2ao(iloc+10,ii)= vg2ao(iloc+10,ii)+tmp(10,ii)*sqrt7
                  vg2ao(iloc+11,ii)= vg2ao(iloc+11,ii)+tmp(11,ii)
                  vg2ao(iloc+12,ii)= vg2ao(iloc+12,ii)+tmp(12,ii)*sqrt7
                  vg2ao(iloc+13,ii)= vg2ao(iloc+13,ii)+tmp(13,ii)*sqrt35third
                  vg2ao(iloc+14,ii)= vg2ao(iloc+14,ii)+tmp(14,ii)*sqrt7
                  vg2ao(iloc+15,ii)= vg2ao(iloc+15,ii)+tmp(15,ii)
                enddo
              else
                do ii= 1,6
                  vg2ao(iloc+1,ii)= vg2ao(iloc+1,ii)+(tmp(2,ii)-tmp(7,ii))*facg1
                  vg2ao(iloc+2,ii)= vg2ao(iloc+2,ii)+(-tmp(12,ii)+tmp(5,ii)*three)*facg2
                  vg2ao(iloc+3,ii)= vg2ao(iloc+3,ii)+(-tmp(2,ii)-tmp(7,ii)+tmp(9,ii)*six)*facg3
                  vg2ao(iloc+4,ii)= vg2ao(iloc+4,ii)+(-tmp(12,ii)*three+tmp(14,ii)*four &
&                                                     -tmp(5,ii)*three)*facg4
                  vg2ao(iloc+5,ii)= vg2ao(iloc+5,ii)+(tmp(1,ii)*three+tmp(11,ii)*three &
&                                                     +tmp(15,ii)*eight+tmp(4,ii)*six &
&                                                     -tmp(6,ii)*p24-tmp(13,ii)*p24)*eighth
                  vg2ao(iloc+6,ii)= vg2ao(iloc+6,ii)+(-tmp(3,ii)*three+tmp(10,ii)*four &
&                                                     -tmp(8,ii)*three)*facg4
                  vg2ao(iloc+7,ii)= vg2ao(iloc+7,ii)+(-tmp(1,ii)+tmp(11,ii)+tmp(6,ii)*six &
&                                                     -tmp(13,ii)*six)*facg5
                  vg2ao(iloc+8,ii)= vg2ao(iloc+8,ii)+(tmp(3,ii)-tmp(8,ii)*three)*facg2
                  vg2ao(iloc+9,ii)= vg2ao(iloc+9,ii)+(tmp(1,ii)+tmp(11,ii)-tmp(4,ii)*six)*facg6
                enddo
              endif
            enddo
          case(5)
            do iprim= 1,numprim
              icount= icount+1
              if(ex(icount)*rsqrd(iatom) > threshex) cycle
              expval= exp(-ex(icount)*rsqrd(iatom))*coeff(icount)
              fac=-ex(icount)*two
              xx= xyzpt(1,iatom)*xyzpt(1,iatom)*fac
              xy= xyzpt(1,iatom)*xyzpt(2,iatom)*fac
              xz= xyzpt(1,iatom)*xyzpt(3,iatom)*fac
              yy= xyzpt(2,iatom)*xyzpt(2,iatom)*fac
              yz= xyzpt(2,iatom)*xyzpt(3,iatom)*fac
              zz= xyzpt(3,iatom)*xyzpt(3,iatom)*fac
              xyz3( 1)= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(1,iatom)
              xyz3( 2)= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(2,iatom)
              xyz3( 3)= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(3,iatom)
              xyz3( 4)= xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)
              xyz3( 5)= xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)
              xyz3( 6)= xyzpt(1,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)
              xyz3( 7)= xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)
              xyz3( 8)= xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)
              xyz3( 9)= xyzpt(2,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)
              xyz3(10)= xyzpt(3,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)
              xyz5( 1)= xyz3( 1)*xx
              xyz5( 2)= xyz3( 1)*xy
              xyz5( 3)= xyz3( 1)*xz
              xyz5( 4)= xyz3( 1)*yy
              xyz5( 5)= xyz3( 1)*yz
              xyz5( 6)= xyz3( 1)*zz
              xyz5( 7)= xyz3( 2)*yy
              xyz5( 8)= xyz3( 2)*yz
              xyz5( 9)= xyz3( 2)*zz
              xyz5(10)= xyz3( 3)*zz
              xyz5(11)= xyz3( 4)*yy
              xyz5(12)= xyz3( 4)*yz
              xyz5(13)= xyz3( 4)*zz
              xyz5(14)= xyz3( 5)*zz
              xyz5(15)= xyz3( 6)*zz
              xyz5(16)= xyz3( 7)*yy
              xyz5(17)= xyz3( 7)*yz
              xyz5(18)= xyz3( 7)*zz
              xyz5(19)= xyz3( 8)*zz
              xyz5(20)= xyz3( 9)*zz
              xyz5(21)= xyz3(10)*zz
              xyz7( 1)= xyz5( 1)*xx
              xyz7( 2)= xyz5( 1)*xy
              xyz7( 3)= xyz5( 1)*xz
              xyz7( 4)= xyz5( 1)*yy
              xyz7( 5)= xyz5( 1)*yz
              xyz7( 6)= xyz5( 1)*zz
              xyz7( 7)= xyz5( 2)*yy
              xyz7( 8)= xyz5( 2)*yz
              xyz7( 9)= xyz5( 2)*zz
              xyz7(10)= xyz5( 3)*zz
              xyz7(11)= xyz5( 4)*yy
              xyz7(12)= xyz5( 4)*yz
              xyz7(13)= xyz5( 4)*zz
              xyz7(14)= xyz5( 6)*yz
              xyz7(15)= xyz5( 6)*zz
              xyz7(16)= xyz5( 7)*yy
              xyz7(17)= xyz5( 7)*yz
              xyz7(18)= xyz5( 7)*zz
              xyz7(19)= xyz5(10)*yy
              xyz7(20)= xyz5(10)*yz
              xyz7(21)= xyz5(10)*zz
              xyz7(22)= xyz5(11)*yy
              xyz7(23)= xyz5(11)*yz
              xyz7(24)= xyz5(11)*zz
              xyz7(25)= xyz5(14)*yy
              xyz7(26)= xyz5(14)*yz
              xyz7(27)= xyz5(14)*zz
              xyz7(28)= xyz5(15)*zz
              xyz7(29)= xyz5(16)*yy
              xyz7(30)= xyz5(16)*yz
              xyz7(31)= xyz5(16)*zz
              xyz7(32)= xyz5(19)*yy
              xyz7(33)= xyz5(19)*yz
              xyz7(34)= xyz5(19)*zz
              xyz7(35)= xyz5(21)*yz
              xyz7(36)= xyz5(21)*zz
              tmp( 1,1)= expval*(xyz3( 1)*p20   +xyz5( 1)*eleven              +xyz7( 1))
              tmp( 2,1)= expval*(xyz3( 2)*twelve+xyz5( 2)*p9                  +xyz7( 2))
              tmp( 3,1)= expval*(xyz3( 3)*twelve+xyz5( 3)*p9                  +xyz7( 3))
              tmp( 4,1)= expval*(xyz3( 4)*six   +xyz5( 4)*seven               +xyz7( 4))
              tmp( 5,1)= expval*(xyz3( 5)*six   +xyz5( 5)*seven               +xyz7( 5))
              tmp( 6,1)= expval*(xyz3( 6)*six   +xyz5( 6)*seven               +xyz7( 6))
              tmp( 7,1)= expval*(xyz3( 7)*two   +xyz5( 7)*five                +xyz7( 7))
              tmp( 8,1)= expval*(xyz3( 8)*two   +xyz5( 8)*five                +xyz7( 8))
              tmp( 9,1)= expval*(xyz3( 9)*two   +xyz5( 9)*five                +xyz7( 9))
              tmp(10,1)= expval*(xyz3(10)*two   +xyz5(10)*five                +xyz7(10))
              tmp(11,1)= expval*(                xyz5(11)*three               +xyz7(11))
              tmp(12,1)= expval*(                xyz5(12)*three               +xyz7(12))
              tmp(13,1)= expval*(                xyz5(13)*three               +xyz7(13))
              tmp(14,1)= expval*(                xyz5(14)*three               +xyz7(14))
              tmp(15,1)= expval*(                xyz5(15)*three               +xyz7(15))
              tmp(16,1)= expval*(                xyz5(16)                     +xyz7(16))
              tmp(17,1)= expval*(                xyz5(17)                     +xyz7(17))
              tmp(18,1)= expval*(                xyz5(18)                     +xyz7(18))
              tmp(19,1)= expval*(                xyz5(19)                     +xyz7(19))
              tmp(20,1)= expval*(                xyz5(20)                     +xyz7(20))
              tmp(21,1)= expval*(                xyz5(21)                     +xyz7(21))
              tmp( 1,2)= expval*(                xyz5( 2)*five                +xyz7( 2))
              tmp( 2,2)= expval*(xyz3( 1)*four  +xyz5( 1)      +xyz5( 4)*four +xyz7( 4))
              tmp( 3,2)= expval*(                               xyz5( 5)*four +xyz7( 5))
              tmp( 4,2)= expval*(xyz3( 2)*six   +xyz5( 2)*two  +xyz5( 7)*three+xyz7( 7))
              tmp( 5,2)= expval*(xyz3( 3)*three +xyz5( 3)      +xyz5( 8)*three+xyz7( 8))
              tmp( 6,2)= expval*(                               xyz5( 9)*three+xyz7( 9))
              tmp( 7,2)= expval*(xyz3( 4)*six   +xyz5( 4)*three+xyz5(11)*two  +xyz7(11))
              tmp( 8,2)= expval*(xyz3( 5)*four  +xyz5( 5)*two  +xyz5(12)*two  +xyz7(12))
              tmp( 9,2)= expval*(xyz3( 6)*two   +xyz5( 6)      +xyz5(13)*two  +xyz7(13))
              tmp(10,2)= expval*(                               xyz5(14)*two  +xyz7(14))
              tmp(11,2)= expval*(xyz3( 7)*four  +xyz5( 7)*four +xyz5(16)      +xyz7(16))
              tmp(12,2)= expval*(xyz3( 8)*three +xyz5( 8)*three+xyz5(17)      +xyz7(17))
              tmp(13,2)= expval*(xyz3( 9)*two   +xyz5( 9)*two  +xyz5(18)      +xyz7(18))
              tmp(14,2)= expval*(xyz3(10)       +xyz5(10)      +xyz5(19)      +xyz7(19))
              tmp(15,2)= expval*(                               xyz5(20)      +xyz7(20))
              tmp(16,2)= expval*(                xyz5(11)*five                +xyz7(22))
              tmp(17,2)= expval*(                xyz5(12)*four                +xyz7(23))
              tmp(18,2)= expval*(                xyz5(13)*three               +xyz7(24))
              tmp(19,2)= expval*(                xyz5(14)*two                 +xyz7(25))
              tmp(20,2)= expval*(                xyz5(15)                     +xyz7(26))
              tmp(21,2)= expval*(                                              xyz7(27))
              tmp( 1,3)= expval*(                               xyz5( 3)*five +xyz7( 3))
              tmp( 2,3)= expval*(                               xyz5( 5)*four +xyz7( 5))
              tmp( 3,3)= expval*(xyz3( 1)*four  +xyz5( 1)      +xyz5( 6)*four +xyz7( 6))
              tmp( 4,3)= expval*(                               xyz5( 8)*three+xyz7( 8))
              tmp( 5,3)= expval*(xyz3( 2)*three +xyz5( 2)      +xyz5( 9)*three+xyz7( 9))
              tmp( 6,3)= expval*(xyz3( 3)*six   +xyz5( 3)*two  +xyz5(10)*three+xyz7(10))
              tmp( 7,3)= expval*(                               xyz5(12)*two  +xyz7(12))
              tmp( 8,3)= expval*(xyz3( 4)*two   +xyz5( 4)      +xyz5(13)*two  +xyz7(13))
              tmp( 9,3)= expval*(xyz3( 5)*four  +xyz5( 5)*two  +xyz5(14)*two  +xyz7(14))
              tmp(10,3)= expval*(xyz3( 6)*six   +xyz5( 6)*three+xyz5(15)*two  +xyz7(15))
              tmp(11,3)= expval*(                               xyz5(17)      +xyz7(17))
              tmp(12,3)= expval*(xyz3( 7)       +xyz5( 7)      +xyz5(18)      +xyz7(18))
              tmp(13,3)= expval*(xyz3( 8)*two   +xyz5( 8)*two  +xyz5(19)      +xyz7(19))
              tmp(14,3)= expval*(xyz3( 9)*three +xyz5( 9)*three+xyz5(20)      +xyz7(20))
              tmp(15,3)= expval*(xyz3(10)*four  +xyz5(10)*four +xyz5(21)      +xyz7(21))
              tmp(16,3)= expval*(                                              xyz7(23))
              tmp(17,3)= expval*(                xyz5(11)                     +xyz7(24))
              tmp(18,3)= expval*(                xyz5(12)*two                 +xyz7(25))
              tmp(19,3)= expval*(                xyz5(13)*three               +xyz7(26))
              tmp(20,3)= expval*(                xyz5(14)*four                +xyz7(27))
              tmp(21,3)= expval*(                xyz5(15)*five                +xyz7(28))
              tmp( 1,4)= expval*(                xyz5( 1)                     +xyz7( 4))
              tmp( 2,4)= expval*(                xyz5( 2)*three               +xyz7( 7))
              tmp( 3,4)= expval*(                xyz5( 3)                     +xyz7( 8))
              tmp( 4,4)= expval*(xyz3( 1)*two   +xyz5( 4)*five                +xyz7(11))
              tmp( 5,4)= expval*(                xyz5( 5)*three               +xyz7(12))
              tmp( 6,4)= expval*(                xyz5( 6)                     +xyz7(13))
              tmp( 7,4)= expval*(xyz3( 2)*six   +xyz5( 7)*seven               +xyz7(16))
              tmp( 8,4)= expval*(xyz3( 3)*two   +xyz5( 8)*five                +xyz7(17))
              tmp( 9,4)= expval*(                xyz5( 9)*three               +xyz7(18))
              tmp(10,4)= expval*(                xyz5(10)                     +xyz7(19))
              tmp(11,4)= expval*(xyz3( 4)*twelve+xyz5(11)*p9                  +xyz7(22))
              tmp(12,4)= expval*(xyz3( 5)*six   +xyz5(12)*seven               +xyz7(23))
              tmp(13,4)= expval*(xyz3( 6)*two   +xyz5(13)*five                +xyz7(24))
              tmp(14,4)= expval*(                xyz5(14)*three               +xyz7(25))
              tmp(15,4)= expval*(                xyz5(15)                     +xyz7(26))
              tmp(16,4)= expval*(xyz3( 7)*p20   +xyz5(16)*eleven              +xyz7(29))
              tmp(17,4)= expval*(xyz3( 8)*twelve+xyz5(17)*p9                  +xyz7(30))
              tmp(18,4)= expval*(xyz3( 9)*six   +xyz5(18)*seven               +xyz7(31))
              tmp(19,4)= expval*(xyz3(10)*two   +xyz5(19)*five                +xyz7(32))
              tmp(20,4)= expval*(                xyz5(20)*three               +xyz7(33))
              tmp(21,4)= expval*(                xyz5(21)                     +xyz7(34))
              tmp( 1,5)= expval*(                                              xyz7( 5))
              tmp( 2,5)= expval*(                               xyz5( 3)      +xyz7( 8))
              tmp( 3,5)= expval*(                xyz5( 2)                     +xyz7( 9))
              tmp( 4,5)= expval*(                               xyz5( 5)*two  +xyz7(12))
              tmp( 5,5)= expval*(xyz3( 1)       +xyz5( 4)      +xyz5( 6)      +xyz7(13))
              tmp( 6,5)= expval*(                xyz5( 5)*two                 +xyz7(14))
              tmp( 7,5)= expval*(                               xyz5( 8)*three+xyz7(17))
              tmp( 8,5)= expval*(xyz3( 2)*two   +xyz5( 7)      +xyz5( 9)*two  +xyz7(18))
              tmp( 9,5)= expval*(xyz3( 3)*two   +xyz5( 8)*two  +xyz5(10)      +xyz7(19))
              tmp(10,5)= expval*(                xyz5( 9)*three               +xyz7(20))
              tmp(11,5)= expval*(                               xyz5(12)*four +xyz7(23))
              tmp(12,5)= expval*(xyz3( 4)*three +xyz5(11)      +xyz5(13)*three+xyz7(24))
              tmp(13,5)= expval*(xyz3( 5)*four  +xyz5(12)*two  +xyz5(14)*two  +xyz7(25))
              tmp(14,5)= expval*(xyz3( 6)*three +xyz5(13)*three+xyz5(15)      +xyz7(26))
              tmp(15,5)= expval*(                xyz5(14)*four                +xyz7(27))
              tmp(16,5)= expval*(                               xyz5(17)*five +xyz7(30))
              tmp(17,5)= expval*(xyz3( 7)*four  +xyz5(16)      +xyz5(18)*four +xyz7(31))
              tmp(18,5)= expval*(xyz3( 8)*six   +xyz5(17)*two  +xyz5(19)*three+xyz7(32))
              tmp(19,5)= expval*(xyz3( 9)*six   +xyz5(18)*three+xyz5(20)*two  +xyz7(33))
              tmp(20,5)= expval*(xyz3(10)*four  +xyz5(19)*four +xyz5(21)      +xyz7(34))
              tmp(21,5)= expval*(                xyz5(20)*five                +xyz7(35))
              tmp( 1,6)= expval*(                xyz5( 1)                     +xyz7( 6))
              tmp( 2,6)= expval*(                xyz5( 2)                     +xyz7( 9))
              tmp( 3,6)= expval*(                xyz5( 3)*three               +xyz7(10))
              tmp( 4,6)= expval*(                xyz5( 4)                     +xyz7(13))
              tmp( 5,6)= expval*(                xyz5( 5)*three               +xyz7(14))
              tmp( 6,6)= expval*(xyz3( 1)*two   +xyz5( 6)*five                +xyz7(15))
              tmp( 7,6)= expval*(                xyz5( 7)                     +xyz7(18))
              tmp( 8,6)= expval*(                xyz5( 8)*three               +xyz7(19))
              tmp( 9,6)= expval*(xyz3( 2)*two   +xyz5( 9)*five                +xyz7(20))
              tmp(10,6)= expval*(xyz3( 3)*six   +xyz5(10)*seven               +xyz7(21))
              tmp(11,6)= expval*(                xyz5(11)                     +xyz7(24))
              tmp(12,6)= expval*(                xyz5(12)*three               +xyz7(25))
              tmp(13,6)= expval*(xyz3( 4)*two   +xyz5(13)*five                +xyz7(26))
              tmp(14,6)= expval*(xyz3( 5)*six   +xyz5(14)*seven               +xyz7(27))
              tmp(15,6)= expval*(xyz3( 6)*twelve+xyz5(15)*p9                  +xyz7(28))
              tmp(16,6)= expval*(                xyz5(16)                     +xyz7(31))
              tmp(17,6)= expval*(                xyz5(17)*three               +xyz7(32))
              tmp(18,6)= expval*(xyz3( 7)*two   +xyz5(18)*five                +xyz7(33))
              tmp(19,6)= expval*(xyz3( 8)*six   +xyz5(19)*seven               +xyz7(34))
              tmp(20,6)= expval*(xyz3( 9)*twelve+xyz5(20)*p9                  +xyz7(35))
              tmp(21,6)= expval*(xyz3(10)*p20   +xyz5(21)*eleven              +xyz7(36))
              if(nbf == 21) then
                do ii= 1,6
                  vg2ao(iloc+ 1,ii)= vg2ao(iloc+ 1,ii)+tmp( 1,ii)
                  vg2ao(iloc+ 2,ii)= vg2ao(iloc+ 2,ii)+tmp( 2,ii)*three
                  vg2ao(iloc+ 3,ii)= vg2ao(iloc+ 3,ii)+tmp( 3,ii)*three
                  vg2ao(iloc+ 4,ii)= vg2ao(iloc+ 4,ii)+tmp( 4,ii)*sqrt21
                  vg2ao(iloc+ 5,ii)= vg2ao(iloc+ 5,ii)+tmp( 5,ii)*sqrt63
                  vg2ao(iloc+ 6,ii)= vg2ao(iloc+ 6,ii)+tmp( 6,ii)*sqrt21
                  vg2ao(iloc+ 7,ii)= vg2ao(iloc+ 7,ii)+tmp( 7,ii)*sqrt21
                  vg2ao(iloc+ 8,ii)= vg2ao(iloc+ 8,ii)+tmp( 8,ii)*sqrt105
                  vg2ao(iloc+ 9,ii)= vg2ao(iloc+ 9,ii)+tmp( 9,ii)*sqrt105
                  vg2ao(iloc+10,ii)= vg2ao(iloc+10,ii)+tmp(10,ii)*sqrt21
                  vg2ao(iloc+11,ii)= vg2ao(iloc+11,ii)+tmp(11,ii)*three
                  vg2ao(iloc+12,ii)= vg2ao(iloc+12,ii)+tmp(12,ii)*sqrt63
                  vg2ao(iloc+13,ii)= vg2ao(iloc+13,ii)+tmp(13,ii)*sqrt105
                  vg2ao(iloc+14,ii)= vg2ao(iloc+14,ii)+tmp(14,ii)*sqrt63
                  vg2ao(iloc+15,ii)= vg2ao(iloc+15,ii)+tmp(15,ii)*three
                  vg2ao(iloc+16,ii)= vg2ao(iloc+16,ii)+tmp(16,ii)
                  vg2ao(iloc+17,ii)= vg2ao(iloc+17,ii)+tmp(17,ii)*three
                  vg2ao(iloc+18,ii)= vg2ao(iloc+18,ii)+tmp(18,ii)*sqrt21
                  vg2ao(iloc+19,ii)= vg2ao(iloc+19,ii)+tmp(19,ii)*sqrt21
                  vg2ao(iloc+20,ii)= vg2ao(iloc+20,ii)+tmp(20,ii)*three
                  vg2ao(iloc+21,ii)= vg2ao(iloc+21,ii)+tmp(21,ii)
                enddo
              else
                do ii= 1,6
                  vg2ao(iloc+ 1,ii)= vg2ao(iloc+ 1,ii) &
&                                   +(tmp(2,ii)*five-tmp(7,ii)*ten+tmp(16,ii))*fach1
                  vg2ao(iloc+ 2,ii)= vg2ao(iloc+ 2,ii) &
&                                   +(tmp(5,ii)*four-tmp(12,ii)*four)*fach2
                  vg2ao(iloc+ 3,ii)= vg2ao(iloc+ 3,ii) &
&                                   +(-tmp(2,ii)*three-tmp(7,ii)*two+tmp(9,ii)*p24 &
&                                   +tmp(16,ii)-tmp(18,ii)*eight)*fach3
                  vg2ao(iloc+ 4,ii)= vg2ao(iloc+ 4,ii) &
&                                   +(-tmp(5,ii)*two-tmp(12,ii)*two+tmp(14,ii)*four)*fach4
                  vg2ao(iloc+ 5,ii)= vg2ao(iloc+ 5,ii) &
&                                   +(tmp(2,ii)+tmp(7,ii)*two-tmp(9,ii)*twelve &
&                                   +tmp(16,ii)-tmp(18,ii)*twelve+tmp(20,ii)*eight)*fach5
                  vg2ao(iloc+ 6,ii)= vg2ao(iloc+ 6,ii) &
&                                   +(tmp(3,ii)*p15+tmp(8,ii)*p30-tmp(10,ii)*p40 &
&                                   +tmp(17,ii)*p15-tmp(19,ii)*p40+tmp(21,ii)*eight)*eighth
                  vg2ao(iloc+ 7,ii)= vg2ao(iloc+ 7,ii) &
&                                   +(tmp(1,ii)+tmp(4,ii)*two-tmp(6,ii)*twelve+tmp(11,ii) &
&                                   -tmp(13,ii)*twelve+tmp(15,ii)*eight)*fach5
                  vg2ao(iloc+ 8,ii)= vg2ao(iloc+ 8,ii) &
&                                   +(-tmp(3,ii)+tmp(10,ii)*two+tmp(17,ii)-tmp(19,ii)*two)*fach4
                  vg2ao(iloc+ 9,ii)= vg2ao(iloc+ 9,ii) &
&                                   +(-tmp(1,ii)+tmp(4,ii)*two+tmp(6,ii)*eight &
&                                   +tmp(11,ii)*three-tmp(13,ii)*p24)*fach3
                  vg2ao(iloc+10,ii)= vg2ao(iloc+10,ii) &
&                                   +(tmp(3,ii)-tmp(8,ii)*six+tmp(17,ii))*fach2
                  vg2ao(iloc+11,ii)= vg2ao(iloc+11,ii) &
&                                   +(tmp(1,ii)-tmp(4,ii)*ten+tmp(11,ii)*five)*fach1
                enddo
              endif
            enddo
          case default
            write(*,'(" Error! Subroutine Gridgao supports up to h functions.")')
            call iabort
        end select
      enddo
!
      return
end


!-------------------------------------------------------------------------------------------------
  subroutine gradrexcor(egrad,edftgrad,cmo,fulldmtrx,atomvec,surface,radpt,angpt,rad,ptweight, &
&                       xyzpt,rsqrd,rr,uvec,vao,vmo,dweight,dpa,pa,transcmo,idftex,idftcor, &
&                       nproc,myrank)
!-------------------------------------------------------------------------------------------------
!
! Driver of derivatives for closed-shell exchange-correlation terms
!
      use modbasis, only : nao
      use modmolecule, only : natom, neleca
      use moddft, only : nrad, nleb
      use modthresh, only : threshweight, threshrho, threshdfock, threshdftao
      implicit none
      integer,intent(in) :: idftex, idftcor, nproc, myrank
      integer :: ngridatom, iatom, irad, ileb, icount, ilebstart, jatom, imo, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00
      real(8),intent(in) :: cmo(nao,neleca), fulldmtrx(nao,nao)
      real(8),intent(in) :: atomvec(5,natom,natom), surface(natom,natom), radpt(2,nrad)
      real(8),intent(in) :: angpt(4,nleb), rad(natom), ptweight(nleb,nrad,natom)
      real(8),intent(out) :: edftgrad(3*natom), xyzpt(3,natom), rsqrd(natom), rr(natom)
      real(8),intent(out) :: uvec(3,natom), vao(nao,10), vmo(neleca,4), dweight(3*natom)
      real(8),intent(out) :: dpa(3,natom,natom), pa(natom), transcmo(neleca,nao)
      real(8),intent(inout) :: egrad(3*natom)
      real(8) :: radpoint, radweight, weight, xgrid, ygrid, zgrid, tmp, rhoa, grhoa(3)
      real(8) :: sphweight, excora(4), ptenergy, wcutoff, rcutoff, fcutoff, aocutoff
!
      transcmo= transpose(cmo)
!
      edftgrad(:)= zero
      ngridatom= nrad*nleb
!
      wcutoff= threshweight/(natom*ngridatom)
      rcutoff= threshrho/(natom*ngridatom)
      fcutoff= threshdfock/(natom*ngridatom)
      aocutoff=threshdftao/(natom*ngridatom)
!
!$OMP parallel do collapse(2) schedule(dynamic,1) private(icount,ilebstart,radpoint,radweight, &
!$OMP weight,xgrid,ygrid,zgrid,xyzpt,rsqrd,rr,tmp,uvec,vao,vmo,rhoa,grhoa,sphweight, &
!$OMP dweight,dpa,pa,excora,ptenergy) reduction(+:edftgrad)
      do iatom= 1,natom
        do irad= 1,nrad
          icount=(iatom-1)*ngridatom+(irad-1)*nleb+1+myrank
          ilebstart=mod(icount,nproc)+1
          radpoint= rad(iatom)*radpt(1,irad)
          radweight= rad(iatom)*rad(iatom)*rad(iatom)*radpt(2,irad)
          do ileb= ilebstart,nleb,nproc
!
            weight=ptweight(ileb,irad,iatom)
            if(weight < wcutoff) cycle
!
            xgrid= radpoint*angpt(1,ileb)
            ygrid= radpoint*angpt(2,ileb)
            zgrid= radpoint*angpt(3,ileb)
            do jatom= 1,natom
              xyzpt(1,jatom)= atomvec(1,iatom,jatom)+xgrid
              xyzpt(2,jatom)= atomvec(2,iatom,jatom)+ygrid
              xyzpt(3,jatom)= atomvec(3,iatom,jatom)+zgrid
              rsqrd(jatom)= xyzpt(1,jatom)*xyzpt(1,jatom)+xyzpt(2,jatom)*xyzpt(2,jatom) &
&                          +xyzpt(3,jatom)*xyzpt(3,jatom)
              rr(jatom)= sqrt(rsqrd(jatom))
              tmp= one/rr(jatom)
              uvec(1,jatom)= xyzpt(1,jatom)*tmp
              uvec(2,jatom)= xyzpt(2,jatom)*tmp
              uvec(3,jatom)= xyzpt(3,jatom)*tmp
            enddo
!
            call gridraomo(vao,vmo,transcmo,xyzpt,rsqrd,aocutoff)
!
            rhoa= zero
            grhoa(1:3)= zero
            do imo= 1,neleca
              rhoa=     rhoa    +vmo(imo,1)*vmo(imo,1)
              grhoa(1)= grhoa(1)+vmo(imo,2)*vmo(imo,1)
              grhoa(2)= grhoa(2)+vmo(imo,3)*vmo(imo,1)
              grhoa(3)= grhoa(3)+vmo(imo,4)*vmo(imo,1)
            enddo
            if((rhoa*two) < rcutoff)cycle
            grhoa(1)= grhoa(1)*two
            grhoa(2)= grhoa(2)*two
            grhoa(3)= grhoa(3)*two
            call gridd2ao(vao(1,5),xyzpt,rsqrd)
!
! Calculate weight derivative
!
            sphweight=radweight*angpt(4,ileb)
            call calcdgridweight(dweight,dpa,pa,uvec,atomvec,surface,rr,sphweight,iatom)
!
            ptenergy= zero
            call calcexcor(excora,excora,ptenergy,rhoa,rhoa,grhoa,grhoa,one,idftex,idftcor,1)
!
            if(abs(excora(1))*two*weight < fcutoff)cycle
!
            call formgradexcor(edftgrad,fulldmtrx,fulldmtrx,vao,vao(1,2),vao(1,5), &
&                              excora,excora,weight,iatom,1)
!
! Add weight derivative contribution
!
            do ii= 1,3*natom
              edftgrad(ii)= edftgrad(ii)+ptenergy*dweight(ii)
            enddo
          enddo
        enddo
      enddo
!$OMP end parallel do
!
! Add derivatives of exchange-correlation terms
!
      do ii= 1,3*natom
        egrad(ii)= egrad(ii)+edftgrad(ii)
      enddo
!
      return
end


!------------------------------------------------------------------------------------------------
  subroutine graduexcor(egrad,edftgrad,cmoa,cmob,fulldmtrx1,fulldmtrx2,atomvec,surface,radpt, &
&                       angpt,rad,ptweight,xyzpt,rsqrd,rr,uvec,vao,vmoa,vmob,dweight, &
&                       dpa,pa,transcmoa,transcmob,idftex,idftcor,nproc,myrank)
!------------------------------------------------------------------------------------------------
!
! Driver of derivatives for open-shell exchange-correlation terms
!
      use modbasis, only : nao
      use modmolecule, only : natom, neleca, nelecb
      use moddft, only : nrad, nleb
      use modthresh, only : threshweight, threshrho, threshdfock, threshdftao
      implicit none
      integer,intent(in) :: idftex, idftcor, nproc, myrank
      integer :: ngridatom, iatom, irad, ileb, icount, ilebstart, jatom, imo, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00
      real(8),intent(in) :: cmoa(nao,neleca), cmob(nao,nelecb), fulldmtrx1(nao,nao)
      real(8),intent(in) :: fulldmtrx2(nao,nao), atomvec(5,natom,natom), surface(natom,natom)
      real(8),intent(in) :: radpt(2,nrad), angpt(4,nleb), rad(natom), ptweight(nleb,nrad,natom)
      real(8),intent(out) :: edftgrad(3*natom), xyzpt(3,natom), rsqrd(natom), rr(natom)
      real(8),intent(out) :: uvec(3,natom), vao(nao,10), vmoa(neleca,4)
      real(8),intent(out) :: vmob(nelecb,4), dweight(3*natom), dpa(3,natom,natom), pa(natom)
      real(8),intent(out) :: transcmoa(neleca,nao), transcmob(nelecb,nao)
      real(8),intent(inout) :: egrad(3*natom)
      real(8) :: radpoint, radweight, weight, xgrid, ygrid, zgrid, tmp, rhoa, rhob
      real(8) :: grhoa(3), grhob(3), sphweight, excora(4), excorb(4), ptenergy
      real(8) :: wcutoff, rcutoff, fcutoff, aocutoff
!
      transcmoa= transpose(cmoa)
      transcmob= transpose(cmob)
!
      edftgrad(:)= zero
      ngridatom= nrad*nleb
!
      wcutoff= threshweight/(natom*ngridatom)
      rcutoff= threshrho/(natom*ngridatom)
      fcutoff= threshdfock/(natom*ngridatom)
      aocutoff=threshdftao/(natom*ngridatom)
!
!$OMP parallel do collapse(2) schedule(dynamic,1) private(icount,ilebstart,radpoint,radweight, &
!$OMP weight,xgrid,ygrid,zgrid,xyzpt,rsqrd,rr,tmp,uvec,vao,vmoa,vmob,rhoa,rhob,grhoa, &
!$OMP grhob,sphweight,dweight,dpa,pa,excora,excorb,ptenergy) reduction(+:edftgrad)
      do iatom= 1,natom
        do irad= 1,nrad
          icount=(iatom-1)*ngridatom+(irad-1)*nleb+1+myrank
          ilebstart=mod(icount,nproc)+1
          radpoint= rad(iatom)*radpt(1,irad)
          radweight= rad(iatom)*rad(iatom)*rad(iatom)*radpt(2,irad)
          do ileb= ilebstart,nleb,nproc
!
            weight=ptweight(ileb,irad,iatom)
            if(weight < wcutoff) cycle
!
            xgrid= radpoint*angpt(1,ileb)
            ygrid= radpoint*angpt(2,ileb)
            zgrid= radpoint*angpt(3,ileb)
            do jatom= 1,natom
              xyzpt(1,jatom)= atomvec(1,iatom,jatom)+xgrid
              xyzpt(2,jatom)= atomvec(2,iatom,jatom)+ygrid
              xyzpt(3,jatom)= atomvec(3,iatom,jatom)+zgrid
              rsqrd(jatom)= xyzpt(1,jatom)*xyzpt(1,jatom)+xyzpt(2,jatom)*xyzpt(2,jatom) &
&                          +xyzpt(3,jatom)*xyzpt(3,jatom)
              rr(jatom)= sqrt(rsqrd(jatom))
              tmp= one/rr(jatom)
              uvec(1,jatom)= xyzpt(1,jatom)*tmp
              uvec(2,jatom)= xyzpt(2,jatom)*tmp
              uvec(3,jatom)= xyzpt(3,jatom)*tmp
            enddo
!
            call griduaomo(vao,vmoa,vmob,transcmoa,transcmob,xyzpt,rsqrd,aocutoff)
!
            rhoa= zero
            grhoa(1:3)= zero
            do imo= 1,neleca
              rhoa=     rhoa    +vmoa(imo,1)*vmoa(imo,1)
              grhoa(1)= grhoa(1)+vmoa(imo,2)*vmoa(imo,1)
              grhoa(2)= grhoa(2)+vmoa(imo,3)*vmoa(imo,1)
              grhoa(3)= grhoa(3)+vmoa(imo,4)*vmoa(imo,1)
            enddo
            rhob= zero
            grhob(1:3)= zero
            do imo= 1,nelecb
              rhob=     rhob    +vmob(imo,1)*vmob(imo,1)
              grhob(1)= grhob(1)+vmob(imo,2)*vmob(imo,1)
              grhob(2)= grhob(2)+vmob(imo,3)*vmob(imo,1)
              grhob(3)= grhob(3)+vmob(imo,4)*vmob(imo,1)
            enddo
            if((rhoa+rhob) < rcutoff)cycle
            grhoa(1:3)= grhoa(1:3)*two
            grhob(1:3)= grhob(1:3)*two
            call gridd2ao(vao(1,5),xyzpt,rsqrd)
!
! Calculate weight derivative
!
            sphweight=radweight*angpt(4,ileb)
            call calcdgridweight(dweight,dpa,pa,uvec,atomvec,surface,rr,sphweight,iatom)
!
            ptenergy= zero
            call calcexcor(excora,excorb,ptenergy,rhoa,rhob,grhoa,grhob,one,idftex,idftcor,2)
!
            if((abs(excora(1))+abs(excorb(1)))*weight < fcutoff)cycle
!
            call formgradexcor(edftgrad,fulldmtrx1,fulldmtrx2,vao,vao(1,2),vao(1,5), &
&                              excora,excorb,weight,iatom,2)
!
! Add weight derivative contribution
!
            do ii= 1,3*natom
              edftgrad(ii)= edftgrad(ii)+ptenergy*dweight(ii)
            enddo
          enddo
        enddo
      enddo
!$OMP end parallel do
!
! Add derivatives of exchange-correlation terms
!
      do ii= 1,3*natom
        egrad(ii)= egrad(ii)+edftgrad(ii)
      enddo
!
      return
end


!-------------------------------------------------------------------------------------
  subroutine calcdgridweight(dweight,dpa,pa,uvec,atomvec,surface,rr,sphweight,katom)
!-------------------------------------------------------------------------------------
!
! Calculate weight derivatives of grid points
!
! In  : uvec     (unit vector in the direction from grid point to atom)
!       atomvec  (tom vector and distance)
!       surface  (surface shifting parameter)
!       sphweight(rad_weight*ang_weight)
! Out : dweight  (derivative of weight)
!       dpa      (work space (gradP))
!       pa       (work space (P))
!
      use modparallel, only : master
      use modmolecule, only : natom
      implicit none
      integer,intent(in) :: katom
      integer :: iatom, jatom, ii
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, oneh=1.5D+00, two=2.0D+00
      real(8),parameter :: pdcoeff=-2.53125D+00 !-81/32
      real(8),parameter :: threshcut=1.0D-12, threshg4=1.0D-08
      real(8),intent(in) :: uvec(3,natom), atomvec(5,natom,natom), surface(natom,natom)
      real(8),intent(in) :: rr(natom), sphweight
      real(8),intent(out) :: dweight(3,natom), dpa(3,natom,natom), pa(natom)
      real(8) :: cutij, cutji, xmuij, zmuij, f4, g4, f2, tmp1, tmp2, dmuji(3), dcoeff
      real(8) :: dcutij, dcutji, weighta, zz, dzz(3)
!
      dpa(:,:,:)= zero
      dweight(:,:)= zero
      pa(:)= one
!
      do iatom= 1,natom
        if(surface(1,iatom) == -one) then
          pa(iatom)= zero
          cycle
        endif
        do jatom= 1,natom
          if(iatom == jatom) cycle
          if(surface(jatom,iatom) == one) cycle
          zmuij=(rr(iatom)-rr(jatom))*atomvec(5,iatom,jatom)
          xmuij= zmuij+surface(jatom,iatom)*(one-zmuij*zmuij)
          f4= xmuij
          g4= one
          do ii= 1,4
            g4= g4*(one-f4*f4)
            f4= f4*(oneh-half*f4*f4)
          enddo
          f2= half*f4
          cutij= half+f2
          cutji= half-f2
          pa(iatom)= pa(iatom)*cutji
!
          if(iatom == katom) cycle
!
          tmp1= one-two*surface(jatom,iatom)*zmuij
          tmp2= zmuij*atomvec(5,iatom,jatom)*atomvec(5,iatom,jatom)
          dmuji(1)= tmp1*(-uvec(1,iatom)*atomvec(5,iatom,jatom)-tmp2*atomvec(1,iatom,jatom))
          dmuji(2)= tmp1*(-uvec(2,iatom)*atomvec(5,iatom,jatom)-tmp2*atomvec(2,iatom,jatom))
          dmuji(3)= tmp1*(-uvec(3,iatom)*atomvec(5,iatom,jatom)-tmp2*atomvec(3,iatom,jatom))
          dcoeff= pdcoeff*g4
!
          if(abs(cutij) > threshcut) then
            dcutij= dcoeff/cutij
            dpa(1,jatom,iatom)= -dcutij*dmuji(1)
            dpa(2,jatom,iatom)= -dcutij*dmuji(2)
            dpa(3,jatom,iatom)= -dcutij*dmuji(3)
          else
            if(abs(g4) > threshg4) then
              if(master) write(*,'(" Error! The value G4 in calcdergridweight is large.")')
              call iabort
            endif
          endif
          if(abs(cutji) > threshcut) then
            dcutji= dcoeff/cutji
            dpa(1,iatom,iatom)= dpa(1,iatom,iatom)+dcutji*dmuji(1)
            dpa(2,iatom,iatom)= dpa(2,iatom,iatom)+dcutji*dmuji(2)
            dpa(3,iatom,iatom)= dpa(3,iatom,iatom)+dcutji*dmuji(3)
          else
            if(abs(g4) > threshg4) then
              if(master) write(*,'(" Error! The value G4 in calcdergridweight is large.")')
              call iabort
            endif
          endif
        enddo
      enddo
!
! Calculate Z(r)
!
      zz= zero
      do iatom= 1,natom
        zz= zz+pa(iatom)           
      enddo
!
! Calculate omega_A(r)
!
      zz= one/zz
      weighta= pa(katom)*zz*sphweight
!
! Calculate grad_B(omega_A)
!
      do iatom= 1,natom
        if(iatom == katom) cycle
        dzz(1:3)=zero
        do jatom= 1,natom
          dzz(1)= dzz(1)+dpa(1,jatom,iatom)*pa(jatom)
          dzz(2)= dzz(2)+dpa(2,jatom,iatom)*pa(jatom)
          dzz(3)= dzz(3)+dpa(3,jatom,iatom)*pa(jatom)
        enddo
!
        dweight(1,iatom)= weighta*(dpa(1,katom,iatom)-dzz(1)*zz)
        dweight(2,iatom)= weighta*(dpa(2,katom,iatom)-dzz(2)*zz)
        dweight(3,iatom)= weighta*(dpa(3,katom,iatom)-dzz(3)*zz)
        dweight(1,katom)= dweight(1,katom)-dweight(1,iatom)
        dweight(2,katom)= dweight(2,katom)-dweight(2,iatom)
        dweight(3,katom)= dweight(3,katom)-dweight(3,iatom)
      enddo
!
      return
end


!-----------------------------------------------------------------------------
  subroutine formgradexcor(edftgrad,fulldmtrxa,fulldmtrxb,vao,vgao,vg2ao, &
&                          excora,excorb,weight,katom,itype)
!-----------------------------------------------------------------------------
!
! Calculate DFT energy gradient terms
!
      use modbasis, only : locatom, mbf, locbf, nao, nshell
      use modmolecule, only : natom
      implicit none
      integer,intent(in) :: katom, itype
      integer :: ish, iatom, nbf, ilocbf, iao, jao
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, two=2.0D+00
      real(8),intent(in) :: fulldmtrxa(nao,nao), fulldmtrxb(nao,nao), vao(nao)
      real(8),intent(in) :: vgao(nao,3), vg2ao(nao,6), excora(4), excorb(4), weight
      real(8),intent(inout) :: edftgrad(3,natom)
      real(8) :: dra, drxa, drya, drza, gradx, grady, gradz, gwx, gwy, gwz
      real(8) :: drb, drxb, dryb, drzb, dmtrxa, dmtrxb
!
! Closed-shell
!
      if(itype == 1) then
        do ish= 1,nshell
          iatom= locatom(ish)
          if(iatom == katom) cycle
          nbf = mbf(ish)
          ilocbf= locbf(ish)
          gradx= zero
          grady= zero
          gradz= zero
          do iao= ilocbf+1,ilocbf+nbf
            dra = zero
            drxa= zero
            drya= zero
            drza= zero
            do jao= 1,nao
              dra = dra +fulldmtrxa(jao,iao)*vao(jao)
              drxa= drxa+fulldmtrxa(jao,iao)*vgao(jao,1)
              drya= drya+fulldmtrxa(jao,iao)*vgao(jao,2)
              drza= drza+fulldmtrxa(jao,iao)*vgao(jao,3)
            enddo
!
            gradx= gradx+vgao(iao,1)*dra*excora(1)
            grady= grady+vgao(iao,2)*dra*excora(1)
            gradz= gradz+vgao(iao,3)*dra*excora(1)
            gradx= gradx+(vg2ao(iao,1)*dra+vgao(iao,1)*drxa)*excora(2) &
&                       +(vg2ao(iao,2)*dra+vgao(iao,1)*drya)*excora(3) &
&                       +(vg2ao(iao,3)*dra+vgao(iao,1)*drza)*excora(4)
            grady= grady+(vg2ao(iao,2)*dra+vgao(iao,2)*drxa)*excora(2) &
&                       +(vg2ao(iao,4)*dra+vgao(iao,2)*drya)*excora(3) &
&                       +(vg2ao(iao,5)*dra+vgao(iao,2)*drza)*excora(4)
            gradz= gradz+(vg2ao(iao,3)*dra+vgao(iao,3)*drxa)*excora(2) &
&                       +(vg2ao(iao,5)*dra+vgao(iao,3)*drya)*excora(3) &
&                       +(vg2ao(iao,6)*dra+vgao(iao,3)*drza)*excora(4)
          enddo
          gwx= two*weight*gradx
          gwy= two*weight*grady
          gwz= two*weight*gradz
          edftgrad(1,iatom)= edftgrad(1,iatom)-gwx
          edftgrad(2,iatom)= edftgrad(2,iatom)-gwy
          edftgrad(3,iatom)= edftgrad(3,iatom)-gwz
          edftgrad(1,katom)= edftgrad(1,katom)+gwx
          edftgrad(2,katom)= edftgrad(2,katom)+gwy
          edftgrad(3,katom)= edftgrad(3,katom)+gwz
        enddo
!
! Open-shell
!
      elseif(itype == 2) then
        do ish= 1,nshell
          iatom= locatom(ish)
          if(iatom == katom) cycle
          nbf = mbf(ish)
          ilocbf= locbf(ish)
          gradx= zero
          grady= zero
          gradz= zero
          do iao= ilocbf+1,ilocbf+nbf
            dra = zero
            drb = zero
            drxa= zero
            drya= zero
            drza= zero
            drxb= zero
            dryb= zero
            drzb= zero
            do jao= 1,nao
              dmtrxa=(fulldmtrxa(jao,iao)+fulldmtrxb(jao,iao))*half
              dmtrxb=(fulldmtrxa(jao,iao)-fulldmtrxb(jao,iao))*half
              dra = dra +dmtrxa*vao(jao)
              drb = drb +dmtrxb*vao(jao)
              drxa= drxa+dmtrxa*vgao(jao,1)
              drya= drya+dmtrxa*vgao(jao,2)
              drza= drza+dmtrxa*vgao(jao,3)
              drxb= drxb+dmtrxb*vgao(jao,1)
              dryb= dryb+dmtrxb*vgao(jao,2)
              drzb= drzb+dmtrxb*vgao(jao,3)
            enddo
!
            gradx= gradx+vgao(iao,1)*dra*excora(1)+vgao(iao,1)*drb*excorb(1)
            grady= grady+vgao(iao,2)*dra*excora(1)+vgao(iao,2)*drb*excorb(1)
            gradz= gradz+vgao(iao,3)*dra*excora(1)+vgao(iao,3)*drb*excorb(1)
            gradx= gradx+(vg2ao(iao,1)*dra+vgao(iao,1)*drxa)*excora(2) &
&                       +(vg2ao(iao,2)*dra+vgao(iao,1)*drya)*excora(3) &
&                       +(vg2ao(iao,3)*dra+vgao(iao,1)*drza)*excora(4) &
&                       +(vg2ao(iao,1)*drb+vgao(iao,1)*drxb)*excorb(2) &
&                       +(vg2ao(iao,2)*drb+vgao(iao,1)*dryb)*excorb(3) &
&                       +(vg2ao(iao,3)*drb+vgao(iao,1)*drzb)*excorb(4)
            grady= grady+(vg2ao(iao,2)*dra+vgao(iao,2)*drxa)*excora(2) &
&                       +(vg2ao(iao,4)*dra+vgao(iao,2)*drya)*excora(3) &
&                       +(vg2ao(iao,5)*dra+vgao(iao,2)*drza)*excora(4) &
&                       +(vg2ao(iao,2)*drb+vgao(iao,2)*drxb)*excorb(2) &
&                       +(vg2ao(iao,4)*drb+vgao(iao,2)*dryb)*excorb(3) &
&                       +(vg2ao(iao,5)*drb+vgao(iao,2)*drzb)*excorb(4)
            gradz= gradz+(vg2ao(iao,3)*dra+vgao(iao,3)*drxa)*excora(2) &
&                       +(vg2ao(iao,5)*dra+vgao(iao,3)*drya)*excora(3) &
&                       +(vg2ao(iao,6)*dra+vgao(iao,3)*drza)*excora(4) &
&                       +(vg2ao(iao,3)*drb+vgao(iao,3)*drxb)*excorb(2) &
&                       +(vg2ao(iao,5)*drb+vgao(iao,3)*dryb)*excorb(3) &
&                       +(vg2ao(iao,6)*drb+vgao(iao,3)*drzb)*excorb(4)
          enddo
          gwx= two*weight*gradx
          gwy= two*weight*grady
          gwz= two*weight*gradz
          edftgrad(1,iatom)= edftgrad(1,iatom)-gwx
          edftgrad(2,iatom)= edftgrad(2,iatom)-gwy
          edftgrad(3,iatom)= edftgrad(3,iatom)-gwz
          edftgrad(1,katom)= edftgrad(1,katom)+gwx
          edftgrad(2,katom)= edftgrad(2,katom)+gwy
          edftgrad(3,katom)= edftgrad(3,katom)+gwz
        enddo
      endif
!
      return
end

