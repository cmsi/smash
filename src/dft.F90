! Copyright 2014  Kazuya Ishimura
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
!---------------------------------------------------------------------------------------
  subroutine formrfockexcor(fockdsum,fockd,energy,totalelec,cmo,atomvec,radpt,angpt, &
&                           rad,ptweight,vao,vmo,xyzpt,rsqrd,transcmo,work,idft, &
&                           nproc,myrank,mpi_comm)
!---------------------------------------------------------------------------------------
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
      use modmolecule, only : natom, numatomic, neleca
      use moddft, only : nrad, nleb
      use modbasis, only : nao
      use modthresh, only : threshweight, threshrho, threshdfock, threshdftao
      implicit none
      integer,intent(in) :: idft, nproc, myrank
      integer(4),intent(in) :: mpi_comm
      integer :: ngridatom, iatom, irad, ileb, icount, ilebstart, jatom, imo
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00
      real(8),intent(in) :: cmo(nao,neleca), atomvec(5,natom,natom), radpt(2,nrad)
      real(8),intent(in) :: angpt(4,nleb), rad(natom), ptweight(nleb,nrad,natom)
      real(8),intent(out) :: fockdsum(nao*(nao+1)/2), fockd(nao*(nao+1)/2), energy, totalelec
      real(8),intent(out) :: vao(nao,4), vmo(neleca,4), xyzpt(3,natom), rsqrd(natom)
      real(8),intent(out) :: transcmo(neleca,nao), work(nao)
      real(8) :: weight, rhoa, grhoa(3), excora(4)
      real(8) :: radpoint, tmp(2,2), wcutoff, rcutoff, fcutoff, aocutoff
!ishimura
    real(8) :: ee1,ee2,ee3
    common/kazuyadft/ee1,ee2,ee3
    ee1=0.0D+0
    ee2=0.0D+0
    ee3=0.0D+0
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
            call calcexcor(excora,excora,energy,rhoa,rhoa,grhoa,grhoa,weight,idft,1)
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


!------------------------------------------------------------------------------------
  subroutine calcexcor(excora,excorb,energy,rhoa,rhob,grhoa,grhob,weight,idft,iscf)
!------------------------------------------------------------------------------------
!
! Calculate exchange and correlation energy at a grid point
!
! In  : idft (1:B3LYP)
!       iscf (1:RHF, 2:UHF)
!
      implicit none
      integer,intent(in) :: idft, iscf
      real(8),parameter :: zero=0.0D+00, onethird=0.3333333333333333D+00, two=2.0D+00
      real(8),parameter :: four=4.0D+00
      real(8),intent(in) :: rhoa, rhob, grhoa(3), grhob(3), weight
      real(8),intent(out) :: excora(4), excorb(4)
      real(8),intent(inout) :: energy
      real(8) :: rhoa13, rhob13, csdlda, cb88, cvwn, clyp, cpw91lda, cpw91
      real(8) :: dummy, vrhoa, vrhob, zk, vsigmaaa, vsigmaab, vsigmabb, gradaa, gradab, gradbb
!ishimura
    real(8) :: ee1,ee2,ee3,eold
    common/kazuyadft/ee1,ee2,ee3
!
      excora(1:4)= zero
      excorb(1:4)= zero
      rhoa13= rhoa**onethird
      rhob13= rhob**onethird
!
! B3LYP
      select case(idft)
        case(1)
          csdlda= 0.08D+00
          cb88=   0.72D+00
          cvwn=   0.19D+00
          clyp=   0.81D+00
!
  eold=energy
          call funcsdlda(excora,excorb,energy,rhoa,rhob,rhoa13,rhob13,weight,csdlda,iscf)
  ee1=ee1+energy-eold
  eold=energy
          call funcbecke88(excora,excorb,energy,rhoa,rhob,grhoa,grhob,rhoa13,rhob13, &
&                        weight,cb88,iscf)
  ee2=ee2+energy-eold
  eold=energy
          call funcvwn5(excora,excorb,energy,rhoa,rhob,weight,cvwn,iscf)
          call funclyp(excora,excorb,energy,rhoa,rhob,grhoa,grhob,rhoa13,rhob13,weight,clyp,iscf)
  ee3=ee3+energy-eold
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


!---------------------------------------------------------------------------------------
  subroutine formufockexcor(fockd1,fockd2,fockd3,energy,totalelec,cmoa,cmob,atomvec, &
&                           radpt,angpt,rad,ptweight,vao,vmoa,vmob,xyzpt,rsqrd, &
&                           transcmoa,transcmob,work,idft,nproc,myrank,mpi_comm)
!---------------------------------------------------------------------------------------
!
! Driver of unrestricted DFT Fock matrix formation from exchange-correlation functionals
!
! In  : cmoa  (Alpha MO coefficient)
!       cmob  (Beta MO coefficient)
! Out : fock1 (Alpha Fock matrix)
!       fock2 (Beta Fock matrix)
!       fock3 (Work space)
!
      use modmolecule, only : natom, numatomic, neleca, nelecb
      use moddft, only : nrad, nleb
      use modbasis, only : nao
      use modthresh, only : threshweight, threshrho, threshdfock, threshdftao
      implicit none
      integer,intent(in) :: idft, nproc, myrank
      integer(4),intent(in) :: mpi_comm
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
!ishimura
    real(8) :: ee1,ee2,ee3
    common/kazuyadft/ee1,ee2,ee3
    ee1=0.0D+0
    ee2=0.0D+0
    ee3=0.0D+0
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
            call calcexcor(excora,excorb,energy,rhoa,rhob,grhoa,grhob,weight,idft,2)
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
        case(74)
          call lebedev74(angpt)
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
        case(230)
          call lebedev230(angpt)
        case(266)
          call lebedev266(angpt)
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
      use modbasis, only : ex, coeff, nshell, nao, locprim, locbf, locatom, &
&                          mprim, mbf, mtype
      use modthresh, only : threshex
      implicit none
      integer :: icount, ish, numprim, iatom, iprim, nang, nbf, ilocbf, ii, jj
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00
      real(8),parameter :: two=2.0D+00, three=3.0D+00, four=4.0D+00
      real(8),parameter :: sqrt3=1.73205080756888D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrt5=2.236067977499790D+00, sqrt15=3.872983346207417D+00
      real(8),parameter :: facf1=0.36969351199675831D+00  ! 1/sqrt(10-6/sqrt(5))
      real(8),parameter :: facf2=0.86602540378443865D+00  ! 1/sqrt(4/3)
      real(8),parameter :: facf3=0.28116020334310144D+00  ! 1/sqrt(46/3-6/sqrt(5))
      real(8),parameter :: facf4=0.24065403274177409D+00  ! 1/sqrt(28-24/sqrt(5))
      real(8),intent(in) :: xyzpt(3,natom), rsqrd(natom), transcmo(neleca,nao), aocutoff
      real(8),intent(out) :: vao(nao,4), vmo(neleca,4)
      real(8) :: expval, fac, tmp(10,4), xx, yy, zz, xy, xz, yz, xxx, xxy, xxz, xyy, xyz
      real(8) :: xzz, yyy, yyz, yzz, zzz, xxxx, xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz
      real(8) :: xyzz, xzzz, yyyy, yyyz, yyzz, yzzz, zzzz
!
      vao(:,:)= zero
      vmo(:,:)= zero
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
              vao(ilocbf+1,1)= vao(ilocbf+1,1)+expval
              fac= ex(icount)*expval*two
              vao(ilocbf+1,2)= vao(ilocbf+1,2)-fac*xyzpt(1,iatom)
              vao(ilocbf+1,3)= vao(ilocbf+1,3)-fac*xyzpt(2,iatom)
              vao(ilocbf+1,4)= vao(ilocbf+1,4)-fac*xyzpt(3,iatom)
            enddo
            if(abs(vao(ilocbf+1,1))+abs(vao(ilocbf+1,2)) &
&             +abs(vao(ilocbf+1,3))+abs(vao(ilocbf+1,4)) > aocutoff) then
              do jj= 1,neleca
                vmo(jj,1)= vmo(jj,1)+vao(ilocbf+1,1)*transcmo(jj,ilocbf+1)
                vmo(jj,2)= vmo(jj,2)+vao(ilocbf+1,2)*transcmo(jj,ilocbf+1)
                vmo(jj,3)= vmo(jj,3)+vao(ilocbf+1,3)*transcmo(jj,ilocbf+1)
                vmo(jj,4)= vmo(jj,4)+vao(ilocbf+1,4)*transcmo(jj,ilocbf+1)
              enddo
            endif
!
! P function
!
          case(1)
            do iprim= 1,numprim
              icount= icount+1
              if(ex(icount)*rsqrd(iatom) > threshex) cycle
              expval= exp(-ex(icount)*rsqrd(iatom))*coeff(icount)
              vao(ilocbf+1,1)= vao(ilocbf+1,1)+expval*xyzpt(1,iatom)
              vao(ilocbf+2,1)= vao(ilocbf+2,1)+expval*xyzpt(2,iatom)
              vao(ilocbf+3,1)= vao(ilocbf+3,1)+expval*xyzpt(3,iatom)
              fac= ex(icount)*expval*two
              xx= xyzpt(1,iatom)*xyzpt(1,iatom)
              yy= xyzpt(2,iatom)*xyzpt(2,iatom)
              zz= xyzpt(3,iatom)*xyzpt(3,iatom)
              xy= xyzpt(1,iatom)*xyzpt(2,iatom)
              xz= xyzpt(1,iatom)*xyzpt(3,iatom)
              yz= xyzpt(2,iatom)*xyzpt(3,iatom)
              vao(ilocbf+1,2)= vao(ilocbf+1,2)+expval-fac*xx
              vao(ilocbf+2,2)= vao(ilocbf+2,2)               -fac*xy
              vao(ilocbf+3,2)= vao(ilocbf+3,2)               -fac*xz
              vao(ilocbf+1,3)= vao(ilocbf+1,3)               -fac*xy
              vao(ilocbf+2,3)= vao(ilocbf+2,3)+expval-fac*yy
              vao(ilocbf+3,3)= vao(ilocbf+3,3)               -fac*yz
              vao(ilocbf+1,4)= vao(ilocbf+1,4)               -fac*xz
              vao(ilocbf+2,4)= vao(ilocbf+2,4)               -fac*yz
              vao(ilocbf+3,4)= vao(ilocbf+3,4)+expval-fac*zz
            enddo
            do ii= 1,3
              if(abs(vao(ilocbf+ii,1))+abs(vao(ilocbf+ii,2)) &
&               +abs(vao(ilocbf+ii,3))+abs(vao(ilocbf+ii,4)) > aocutoff) then
                do jj= 1,neleca
                  vmo(jj,1)= vmo(jj,1)+vao(ilocbf+ii,1)*transcmo(jj,ilocbf+ii)
                  vmo(jj,2)= vmo(jj,2)+vao(ilocbf+ii,2)*transcmo(jj,ilocbf+ii)
                  vmo(jj,3)= vmo(jj,3)+vao(ilocbf+ii,3)*transcmo(jj,ilocbf+ii)
                  vmo(jj,4)= vmo(jj,4)+vao(ilocbf+ii,4)*transcmo(jj,ilocbf+ii)
                enddo
              endif
            enddo
!
! D function
!
          case(2)
            do iprim= 1,numprim
              icount= icount+1
              if(ex(icount)*rsqrd(iatom) > threshex) cycle
              expval= exp(-ex(icount)*rsqrd(iatom))*coeff(icount)
              tmp(1,1)= expval*xyzpt(1,iatom)*xyzpt(1,iatom)
              tmp(2,1)= expval*xyzpt(2,iatom)*xyzpt(2,iatom)
              tmp(3,1)= expval*xyzpt(3,iatom)*xyzpt(3,iatom)
              tmp(4,1)= expval*xyzpt(1,iatom)*xyzpt(2,iatom)*sqrt3
              tmp(5,1)= expval*xyzpt(1,iatom)*xyzpt(3,iatom)*sqrt3
              tmp(6,1)= expval*xyzpt(2,iatom)*xyzpt(3,iatom)*sqrt3
              fac= ex(icount)*two
              xxx=-xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(1,iatom)*fac
              xxy=-xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(2,iatom)*fac
              xxz=-xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(3,iatom)*fac
              xyy=-xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)*fac
              xyz=-xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)*fac
              xzz=-xyzpt(1,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac
              yyy=-xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)*fac
              yyz=-xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)*fac
              yzz=-xyzpt(2,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac
              zzz=-xyzpt(3,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac
              tmp(1,2)= expval*(xyzpt(1,iatom)*two+xxx)
              tmp(2,2)= expval*(                  +xyy)
              tmp(3,2)= expval*(                  +xzz)
              tmp(4,2)= expval*(xyzpt(2,iatom)    +xxy)*sqrt3
              tmp(5,2)= expval*(xyzpt(3,iatom)    +xxz)*sqrt3
              tmp(6,2)= expval*(                  +xyz)*sqrt3
              tmp(1,3)= expval*(                  +xxy)
              tmp(2,3)= expval*(xyzpt(2,iatom)*two+yyy)
              tmp(3,3)= expval*(                  +yzz)
              tmp(4,3)= expval*(xyzpt(1,iatom)    +xyy)*sqrt3
              tmp(5,3)= expval*(                  +xyz)*sqrt3
              tmp(6,3)= expval*(xyzpt(3,iatom)    +yyz)*sqrt3
              tmp(1,4)= expval*(                  +xxz)
              tmp(2,4)= expval*(                  +yyz)
              tmp(3,4)= expval*(xyzpt(3,iatom)*two+zzz)
              tmp(4,4)= expval*(                  +xyz)*sqrt3
              tmp(5,4)= expval*(xyzpt(1,iatom)    +xzz)*sqrt3
              tmp(6,4)= expval*(xyzpt(2,iatom)    +yzz)*sqrt3
              if(nbf == 6) then
                do jj= 1,4
                  do ii= 1,6
                    vao(ilocbf+ii,jj)= vao(ilocbf+ii,jj)+tmp(ii,jj)
                  enddo
                enddo
              else
                do jj= 1,4
                  vao(ilocbf+1,jj)= vao(ilocbf+1,jj)+tmp(3,jj)-(tmp(1,jj)+tmp(2,jj))*half
                  vao(ilocbf+2,jj)= vao(ilocbf+2,jj)+tmp(5,jj)
                  vao(ilocbf+3,jj)= vao(ilocbf+3,jj)+tmp(6,jj)
                  vao(ilocbf+4,jj)= vao(ilocbf+4,jj)+(tmp(1,jj)-tmp(2,jj))*sqrt3h
                  vao(ilocbf+5,jj)= vao(ilocbf+5,jj)+tmp(4,jj)
                enddo
              endif
            enddo
            do ii= 1,nbf
              if(abs(vao(ilocbf+ii,1))+abs(vao(ilocbf+ii,2)) &
&               +abs(vao(ilocbf+ii,3))+abs(vao(ilocbf+ii,4)) > aocutoff) then
                do jj= 1,neleca
                  vmo(jj,1)= vmo(jj,1)+vao(ilocbf+ii,1)*transcmo(jj,ilocbf+ii)
                  vmo(jj,2)= vmo(jj,2)+vao(ilocbf+ii,2)*transcmo(jj,ilocbf+ii)
                  vmo(jj,3)= vmo(jj,3)+vao(ilocbf+ii,3)*transcmo(jj,ilocbf+ii)
                  vmo(jj,4)= vmo(jj,4)+vao(ilocbf+ii,4)*transcmo(jj,ilocbf+ii)
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
              xx= xyzpt(1,iatom)*xyzpt(1,iatom)
              yy= xyzpt(2,iatom)*xyzpt(2,iatom)
              zz= xyzpt(3,iatom)*xyzpt(3,iatom)
              xy= xyzpt(1,iatom)*xyzpt(2,iatom)
              xz= xyzpt(1,iatom)*xyzpt(3,iatom)
              yz= xyzpt(2,iatom)*xyzpt(3,iatom)
              fac= ex(icount)*two
              xxxx=-xx*xx*fac
              xxxy=-xx*xy*fac
              xxxz=-xx*xz*fac
              xxyy=-xx*yy*fac
              xxyz=-xx*yz*fac
              xxzz=-xx*zz*fac
              xyyy=-xy*yy*fac
              xyyz=-xy*yz*fac
              xyzz=-xy*zz*fac
              xzzz=-xz*zz*fac
              yyyy=-yy*yy*fac
              yyyz=-yy*yz*fac
              yyzz=-yy*zz*fac
              yzzz=-yz*zz*fac
              zzzz=-zz*zz*fac
              tmp( 1,1)= expval*xx*xyzpt(1,iatom)
              tmp( 2,1)= expval*yy*xyzpt(2,iatom)
              tmp( 3,1)= expval*zz*xyzpt(3,iatom)
              tmp( 4,1)= expval*xx*xyzpt(2,iatom)*sqrt5
              tmp( 5,1)= expval*xx*xyzpt(3,iatom)*sqrt5
              tmp( 6,1)= expval*xy*xyzpt(2,iatom)*sqrt5
              tmp( 7,1)= expval*yy*xyzpt(3,iatom)*sqrt5
              tmp( 8,1)= expval*xz*xyzpt(3,iatom)*sqrt5
              tmp( 9,1)= expval*yz*xyzpt(3,iatom)*sqrt5
              tmp(10,1)= expval*xy*xyzpt(3,iatom)*sqrt15
              tmp( 1,2)= expval*(xx*three+xxxx)
              tmp( 2,2)= expval*(        +xyyy)
              tmp( 3,2)= expval*(        +xzzz)
              tmp( 4,2)= expval*(xy*two  +xxxy)*sqrt5
              tmp( 5,2)= expval*(xz*two  +xxxz)*sqrt5
              tmp( 6,2)= expval*(yy      +xxyy)*sqrt5
              tmp( 7,2)= expval*(        +xyyz)*sqrt5
              tmp( 8,2)= expval*(zz      +xxzz)*sqrt5
              tmp( 9,2)= expval*(        +xyzz)*sqrt5
              tmp(10,2)= expval*(yz      +xxyz)*sqrt15
              tmp( 1,3)= expval*(        +xxxy)
              tmp( 2,3)= expval*(yy*three+yyyy)
              tmp( 3,3)= expval*(        +yzzz)
              tmp( 4,3)= expval*(xx      +xxyy)*sqrt5
              tmp( 5,3)= expval*(        +xxyz)*sqrt5
              tmp( 6,3)= expval*(xy*two  +xyyy)*sqrt5
              tmp( 7,3)= expval*(yz*two  +yyyz)*sqrt5
              tmp( 8,3)= expval*(        +xyzz)*sqrt5
              tmp( 9,3)= expval*(zz      +yyzz)*sqrt5
              tmp(10,3)= expval*(xz      +xyyz)*sqrt15
              tmp( 1,4)= expval*(        +xxxz)
              tmp( 2,4)= expval*(        +yyyz)
              tmp( 3,4)= expval*(zz*three+zzzz)
              tmp( 4,4)= expval*(        +xxyz)*sqrt5
              tmp( 5,4)= expval*(xx      +xxzz)*sqrt5
              tmp( 6,4)= expval*(        +xyyz)*sqrt5
              tmp( 7,4)= expval*(yy      +yyzz)*sqrt5
              tmp( 8,4)= expval*(xz*two  +xzzz)*sqrt5
              tmp( 9,4)= expval*(yz*two  +yzzz)*sqrt5
              tmp(10,4)= expval*(xy      +xyzz)*sqrt15
              if(nbf == 10) then
                do jj= 1,4
                  do ii= 1,10
                    vao(ilocbf+ii,jj)= vao(ilocbf+ii,jj)+tmp(ii,jj)
                  enddo
                enddo
              else
                do jj= 1,4
                  vao(ilocbf+1,jj)= vao(ilocbf+1,jj)+( two*tmp(3,jj)-three*tmp(5,jj) &
&                                                                   -three*tmp(7,jj)     )*facf4
                  vao(ilocbf+2,jj)= vao(ilocbf+2,jj)+(-tmp(1,jj)-tmp(6,jj)+four*tmp(8,jj))*facf3
                  vao(ilocbf+3,jj)= vao(ilocbf+3,jj)+(-tmp(2,jj)-tmp(4,jj)+four*tmp(9,jj))*facf3
                  vao(ilocbf+4,jj)= vao(ilocbf+4,jj)+( tmp(5,jj)-tmp(7,jj)               )*facf2
                  vao(ilocbf+5,jj)= vao(ilocbf+5,jj)+  tmp(10,jj)
                  vao(ilocbf+6,jj)= vao(ilocbf+6,jj)+( tmp(1,jj)-three*tmp(6,jj)         )*facf1
                  vao(ilocbf+7,jj)= vao(ilocbf+7,jj)+(-tmp(2,jj)+three*tmp(4,jj)         )*facf1
                enddo
              endif
            enddo
            do ii= 1,nbf
              if(abs(vao(ilocbf+ii,1))+abs(vao(ilocbf+ii,2)) &
&               +abs(vao(ilocbf+ii,3))+abs(vao(ilocbf+ii,4)) > aocutoff) then
                do jj= 1,neleca
                  vmo(jj,1)= vmo(jj,1)+vao(ilocbf+ii,1)*transcmo(jj,ilocbf+ii)
                  vmo(jj,2)= vmo(jj,2)+vao(ilocbf+ii,2)*transcmo(jj,ilocbf+ii)
                  vmo(jj,3)= vmo(jj,3)+vao(ilocbf+ii,3)*transcmo(jj,ilocbf+ii)
                  vmo(jj,4)= vmo(jj,4)+vao(ilocbf+ii,4)*transcmo(jj,ilocbf+ii)
                enddo
              endif
            enddo
          case default
            write(*,'(" Error! Subroutine Gridraomo supports up to f functions.")')
            call iabort
        end select
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
      use modbasis, only : ex, coeff, nshell, nao, locprim, locbf, locatom, &
&                          mprim, mbf, mtype
      use modthresh, only : threshex
      implicit none
      integer :: icount, ish, numprim, iatom, iprim, nang, nbf, ilocbf, ii, jj
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00
      real(8),parameter :: two=2.0D+00, three=3.0D+00, four=4.0D+00
      real(8),parameter :: sqrt3=1.73205080756888D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrt5=2.236067977499790D+00, sqrt15=3.872983346207417D+00
      real(8),parameter :: facf1=0.36969351199675831D+00  ! 1/sqrt(10-6/sqrt(5))
      real(8),parameter :: facf2=0.86602540378443865D+00  ! 1/sqrt(4/3)
      real(8),parameter :: facf3=0.28116020334310144D+00  ! 1/sqrt(46/3-6/sqrt(5))
      real(8),parameter :: facf4=0.24065403274177409D+00  ! 1/sqrt(28-24/sqrt(5))
      real(8),intent(in) :: xyzpt(3,natom), rsqrd(natom), transcmoa(neleca,nao)
      real(8),intent(in) :: transcmob(nelecb,nao), aocutoff
      real(8),intent(out) :: vao(nao,4), vmoa(neleca,4), vmob(nelecb,4)
      real(8) :: expval, fac, tmp(10,4), xx, yy, zz, xy, xz, yz, xxx, xxy, xxz, xyy, xyz
      real(8) :: xzz, yyy, yyz, yzz, zzz, xxxx, xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz
      real(8) :: xyzz, xzzz, yyyy, yyyz, yyzz, yzzz, zzzz
!
      vao(:,:)= zero
      vmoa(:,:)= zero
      vmob(:,:)= zero
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
              vao(ilocbf+1,1)= vao(ilocbf+1,1)+expval
              fac= ex(icount)*expval*two
              vao(ilocbf+1,2)= vao(ilocbf+1,2)-fac*xyzpt(1,iatom)
              vao(ilocbf+1,3)= vao(ilocbf+1,3)-fac*xyzpt(2,iatom)
              vao(ilocbf+1,4)= vao(ilocbf+1,4)-fac*xyzpt(3,iatom)
            enddo
            if(abs(vao(ilocbf+1,1))+abs(vao(ilocbf+1,2)) &
&             +abs(vao(ilocbf+1,3))+abs(vao(ilocbf+1,4)) > aocutoff) then
              do jj= 1,neleca
                vmoa(jj,1)= vmoa(jj,1)+vao(ilocbf+1,1)*transcmoa(jj,ilocbf+1)
                vmoa(jj,2)= vmoa(jj,2)+vao(ilocbf+1,2)*transcmoa(jj,ilocbf+1)
                vmoa(jj,3)= vmoa(jj,3)+vao(ilocbf+1,3)*transcmoa(jj,ilocbf+1)
                vmoa(jj,4)= vmoa(jj,4)+vao(ilocbf+1,4)*transcmoa(jj,ilocbf+1)
              enddo
              do jj= 1,nelecb
                vmob(jj,1)= vmob(jj,1)+vao(ilocbf+1,1)*transcmob(jj,ilocbf+1)
                vmob(jj,2)= vmob(jj,2)+vao(ilocbf+1,2)*transcmob(jj,ilocbf+1)
                vmob(jj,3)= vmob(jj,3)+vao(ilocbf+1,3)*transcmob(jj,ilocbf+1)
                vmob(jj,4)= vmob(jj,4)+vao(ilocbf+1,4)*transcmob(jj,ilocbf+1)
              enddo
            endif
!
! P function
!
          case(1)
            do iprim= 1,numprim
              icount= icount+1
              if(ex(icount)*rsqrd(iatom) > threshex) cycle
              expval= exp(-ex(icount)*rsqrd(iatom))*coeff(icount)
              vao(ilocbf+1,1)= vao(ilocbf+1,1)+expval*xyzpt(1,iatom)
              vao(ilocbf+2,1)= vao(ilocbf+2,1)+expval*xyzpt(2,iatom)
              vao(ilocbf+3,1)= vao(ilocbf+3,1)+expval*xyzpt(3,iatom)
              fac= ex(icount)*expval*two
              xx= xyzpt(1,iatom)*xyzpt(1,iatom)
              yy= xyzpt(2,iatom)*xyzpt(2,iatom)
              zz= xyzpt(3,iatom)*xyzpt(3,iatom)
              xy= xyzpt(1,iatom)*xyzpt(2,iatom)
              xz= xyzpt(1,iatom)*xyzpt(3,iatom)
              yz= xyzpt(2,iatom)*xyzpt(3,iatom)
              vao(ilocbf+1,2)= vao(ilocbf+1,2)+expval-fac*xx
              vao(ilocbf+2,2)= vao(ilocbf+2,2)               -fac*xy
              vao(ilocbf+3,2)= vao(ilocbf+3,2)               -fac*xz
              vao(ilocbf+1,3)= vao(ilocbf+1,3)               -fac*xy
              vao(ilocbf+2,3)= vao(ilocbf+2,3)+expval-fac*yy
              vao(ilocbf+3,3)= vao(ilocbf+3,3)               -fac*yz
              vao(ilocbf+1,4)= vao(ilocbf+1,4)               -fac*xz
              vao(ilocbf+2,4)= vao(ilocbf+2,4)               -fac*yz
              vao(ilocbf+3,4)= vao(ilocbf+3,4)+expval-fac*zz
            enddo
            do ii= 1,3
              if(abs(vao(ilocbf+ii,1))+abs(vao(ilocbf+ii,2)) &
&               +abs(vao(ilocbf+ii,3))+abs(vao(ilocbf+ii,4)) > aocutoff) then
                do jj= 1,neleca
                  vmoa(jj,1)= vmoa(jj,1)+vao(ilocbf+ii,1)*transcmoa(jj,ilocbf+ii)
                  vmoa(jj,2)= vmoa(jj,2)+vao(ilocbf+ii,2)*transcmoa(jj,ilocbf+ii)
                  vmoa(jj,3)= vmoa(jj,3)+vao(ilocbf+ii,3)*transcmoa(jj,ilocbf+ii)
                  vmoa(jj,4)= vmoa(jj,4)+vao(ilocbf+ii,4)*transcmoa(jj,ilocbf+ii)
                enddo
                do jj= 1,nelecb
                  vmob(jj,1)= vmob(jj,1)+vao(ilocbf+ii,1)*transcmob(jj,ilocbf+ii)
                  vmob(jj,2)= vmob(jj,2)+vao(ilocbf+ii,2)*transcmob(jj,ilocbf+ii)
                  vmob(jj,3)= vmob(jj,3)+vao(ilocbf+ii,3)*transcmob(jj,ilocbf+ii)
                  vmob(jj,4)= vmob(jj,4)+vao(ilocbf+ii,4)*transcmob(jj,ilocbf+ii)
                enddo
              endif
            enddo
!
! D function
!
          case(2)
            do iprim= 1,numprim
              icount= icount+1
              if(ex(icount)*rsqrd(iatom) > threshex) cycle
              expval= exp(-ex(icount)*rsqrd(iatom))*coeff(icount)
              tmp(1,1)= expval*xyzpt(1,iatom)*xyzpt(1,iatom)
              tmp(2,1)= expval*xyzpt(2,iatom)*xyzpt(2,iatom)
              tmp(3,1)= expval*xyzpt(3,iatom)*xyzpt(3,iatom)
              tmp(4,1)= expval*xyzpt(1,iatom)*xyzpt(2,iatom)*sqrt3
              tmp(5,1)= expval*xyzpt(1,iatom)*xyzpt(3,iatom)*sqrt3
              tmp(6,1)= expval*xyzpt(2,iatom)*xyzpt(3,iatom)*sqrt3
              fac= ex(icount)*two
              xxx=-xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(1,iatom)*fac
              xxy=-xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(2,iatom)*fac
              xxz=-xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(3,iatom)*fac
              xyy=-xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)*fac
              xyz=-xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)*fac
              xzz=-xyzpt(1,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac
              yyy=-xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)*fac
              yyz=-xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)*fac
              yzz=-xyzpt(2,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac
              zzz=-xyzpt(3,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac
              tmp(1,2)= expval*(xyzpt(1,iatom)*two+xxx)
              tmp(2,2)= expval*(                  +xyy)
              tmp(3,2)= expval*(                  +xzz)
              tmp(4,2)= expval*(xyzpt(2,iatom)    +xxy)*sqrt3
              tmp(5,2)= expval*(xyzpt(3,iatom)    +xxz)*sqrt3
              tmp(6,2)= expval*(                  +xyz)*sqrt3
              tmp(1,3)= expval*(                  +xxy)
              tmp(2,3)= expval*(xyzpt(2,iatom)*two+yyy)
              tmp(3,3)= expval*(                  +yzz)
              tmp(4,3)= expval*(xyzpt(1,iatom)    +xyy)*sqrt3
              tmp(5,3)= expval*(                  +xyz)*sqrt3
              tmp(6,3)= expval*(xyzpt(3,iatom)    +yyz)*sqrt3
              tmp(1,4)= expval*(                  +xxz)
              tmp(2,4)= expval*(                  +yyz)
              tmp(3,4)= expval*(xyzpt(3,iatom)*two+zzz)
              tmp(4,4)= expval*(                  +xyz)*sqrt3
              tmp(5,4)= expval*(xyzpt(1,iatom)    +xzz)*sqrt3
              tmp(6,4)= expval*(xyzpt(2,iatom)    +yzz)*sqrt3
              if(nbf == 6) then
                do jj= 1,4
                  do ii= 1,6
                    vao(ilocbf+ii,jj)= vao(ilocbf+ii,jj)+tmp(ii,jj)
                  enddo
                enddo
              else
                do jj= 1,4
                  vao(ilocbf+1,jj)= vao(ilocbf+1,jj)+tmp(3,jj)-(tmp(1,jj)+tmp(2,jj))*half
                  vao(ilocbf+2,jj)= vao(ilocbf+2,jj)+tmp(5,jj)
                  vao(ilocbf+3,jj)= vao(ilocbf+3,jj)+tmp(6,jj)
                  vao(ilocbf+4,jj)= vao(ilocbf+4,jj)+(tmp(1,jj)-tmp(2,jj))*sqrt3h
                  vao(ilocbf+5,jj)= vao(ilocbf+5,jj)+tmp(4,jj)
                enddo
              endif
            enddo
            do ii= 1,nbf
              if(abs(vao(ilocbf+ii,1))+abs(vao(ilocbf+ii,2)) &
&               +abs(vao(ilocbf+ii,3))+abs(vao(ilocbf+ii,4)) > aocutoff) then
                do jj= 1,neleca
                  vmoa(jj,1)= vmoa(jj,1)+vao(ilocbf+ii,1)*transcmoa(jj,ilocbf+ii)
                  vmoa(jj,2)= vmoa(jj,2)+vao(ilocbf+ii,2)*transcmoa(jj,ilocbf+ii)
                  vmoa(jj,3)= vmoa(jj,3)+vao(ilocbf+ii,3)*transcmoa(jj,ilocbf+ii)
                  vmoa(jj,4)= vmoa(jj,4)+vao(ilocbf+ii,4)*transcmoa(jj,ilocbf+ii)
                enddo
                do jj= 1,nelecb
                  vmob(jj,1)= vmob(jj,1)+vao(ilocbf+ii,1)*transcmob(jj,ilocbf+ii)
                  vmob(jj,2)= vmob(jj,2)+vao(ilocbf+ii,2)*transcmob(jj,ilocbf+ii)
                  vmob(jj,3)= vmob(jj,3)+vao(ilocbf+ii,3)*transcmob(jj,ilocbf+ii)
                  vmob(jj,4)= vmob(jj,4)+vao(ilocbf+ii,4)*transcmob(jj,ilocbf+ii)
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
              xx= xyzpt(1,iatom)*xyzpt(1,iatom)
              yy= xyzpt(2,iatom)*xyzpt(2,iatom)
              zz= xyzpt(3,iatom)*xyzpt(3,iatom)
              xy= xyzpt(1,iatom)*xyzpt(2,iatom)
              xz= xyzpt(1,iatom)*xyzpt(3,iatom)
              yz= xyzpt(2,iatom)*xyzpt(3,iatom)
              fac= ex(icount)*two
              xxxx=-xx*xx*fac
              xxxy=-xx*xy*fac
              xxxz=-xx*xz*fac
              xxyy=-xx*yy*fac
              xxyz=-xx*yz*fac
              xxzz=-xx*zz*fac
              xyyy=-xy*yy*fac
              xyyz=-xy*yz*fac
              xyzz=-xy*zz*fac
              xzzz=-xz*zz*fac
              yyyy=-yy*yy*fac
              yyyz=-yy*yz*fac
              yyzz=-yy*zz*fac
              yzzz=-yz*zz*fac
              zzzz=-zz*zz*fac
              tmp( 1,1)= expval*xx*xyzpt(1,iatom)
              tmp( 2,1)= expval*yy*xyzpt(2,iatom)
              tmp( 3,1)= expval*zz*xyzpt(3,iatom)
              tmp( 4,1)= expval*xx*xyzpt(2,iatom)*sqrt5
              tmp( 5,1)= expval*xx*xyzpt(3,iatom)*sqrt5
              tmp( 6,1)= expval*xy*xyzpt(2,iatom)*sqrt5
              tmp( 7,1)= expval*yy*xyzpt(3,iatom)*sqrt5
              tmp( 8,1)= expval*xz*xyzpt(3,iatom)*sqrt5
              tmp( 9,1)= expval*yz*xyzpt(3,iatom)*sqrt5
              tmp(10,1)= expval*xy*xyzpt(3,iatom)*sqrt15
              tmp( 1,2)= expval*(xx*three+xxxx)
              tmp( 2,2)= expval*(        +xyyy)
              tmp( 3,2)= expval*(        +xzzz)
              tmp( 4,2)= expval*(xy*two  +xxxy)*sqrt5
              tmp( 5,2)= expval*(xz*two  +xxxz)*sqrt5
              tmp( 6,2)= expval*(yy      +xxyy)*sqrt5
              tmp( 7,2)= expval*(        +xyyz)*sqrt5
              tmp( 8,2)= expval*(zz      +xxzz)*sqrt5
              tmp( 9,2)= expval*(        +xyzz)*sqrt5
              tmp(10,2)= expval*(yz      +xxyz)*sqrt15
              tmp( 1,3)= expval*(        +xxxy)
              tmp( 2,3)= expval*(yy*three+yyyy)
              tmp( 3,3)= expval*(        +yzzz)
              tmp( 4,3)= expval*(xx      +xxyy)*sqrt5
              tmp( 5,3)= expval*(        +xxyz)*sqrt5
              tmp( 6,3)= expval*(xy*two  +xyyy)*sqrt5
              tmp( 7,3)= expval*(yz*two  +yyyz)*sqrt5
              tmp( 8,3)= expval*(        +xyzz)*sqrt5
              tmp( 9,3)= expval*(zz      +yyzz)*sqrt5
              tmp(10,3)= expval*(xz      +xyyz)*sqrt15
              tmp( 1,4)= expval*(        +xxxz)
              tmp( 2,4)= expval*(        +yyyz)
              tmp( 3,4)= expval*(zz*three+zzzz)
              tmp( 4,4)= expval*(        +xxyz)*sqrt5
              tmp( 5,4)= expval*(xx      +xxzz)*sqrt5
              tmp( 6,4)= expval*(        +xyyz)*sqrt5
              tmp( 7,4)= expval*(yy      +yyzz)*sqrt5
              tmp( 8,4)= expval*(xz*two  +xzzz)*sqrt5
              tmp( 9,4)= expval*(yz*two  +yzzz)*sqrt5
              tmp(10,4)= expval*(xy      +xyzz)*sqrt15
              if(nbf == 10) then
                do jj= 1,4
                  do ii= 1,10
                    vao(ilocbf+ii,jj)= vao(ilocbf+ii,jj)+tmp(ii,jj)
                  enddo
                enddo
              else
                do jj= 1,4
                  vao(ilocbf+1,jj)= vao(ilocbf+1,jj)+( two*tmp(3,jj)-three*tmp(5,jj) &
&                                                                   -three*tmp(7,jj)     )*facf4
                  vao(ilocbf+2,jj)= vao(ilocbf+2,jj)+(-tmp(1,jj)-tmp(6,jj)+four*tmp(8,jj))*facf3
                  vao(ilocbf+3,jj)= vao(ilocbf+3,jj)+(-tmp(2,jj)-tmp(4,jj)+four*tmp(9,jj))*facf3
                  vao(ilocbf+4,jj)= vao(ilocbf+4,jj)+( tmp(5,jj)-tmp(7,jj)               )*facf2
                  vao(ilocbf+5,jj)= vao(ilocbf+5,jj)+  tmp(10,jj)
                  vao(ilocbf+6,jj)= vao(ilocbf+6,jj)+( tmp(1,jj)-three*tmp(6,jj)         )*facf1
                  vao(ilocbf+7,jj)= vao(ilocbf+7,jj)+(-tmp(2,jj)+three*tmp(4,jj)         )*facf1
                enddo
              endif
            enddo
            do ii= 1,nbf
              if(abs(vao(ilocbf+ii,1))+abs(vao(ilocbf+ii,2)) &
&               +abs(vao(ilocbf+ii,3))+abs(vao(ilocbf+ii,4)) > aocutoff) then
                do jj= 1,neleca
                  vmoa(jj,1)= vmoa(jj,1)+vao(ilocbf+ii,1)*transcmoa(jj,ilocbf+ii)
                  vmoa(jj,2)= vmoa(jj,2)+vao(ilocbf+ii,2)*transcmoa(jj,ilocbf+ii)
                  vmoa(jj,3)= vmoa(jj,3)+vao(ilocbf+ii,3)*transcmoa(jj,ilocbf+ii)
                  vmoa(jj,4)= vmoa(jj,4)+vao(ilocbf+ii,4)*transcmoa(jj,ilocbf+ii)
                enddo
                do jj= 1,nelecb
                  vmob(jj,1)= vmob(jj,1)+vao(ilocbf+ii,1)*transcmob(jj,ilocbf+ii)
                  vmob(jj,2)= vmob(jj,2)+vao(ilocbf+ii,2)*transcmob(jj,ilocbf+ii)
                  vmob(jj,3)= vmob(jj,3)+vao(ilocbf+ii,3)*transcmob(jj,ilocbf+ii)
                  vmob(jj,4)= vmob(jj,4)+vao(ilocbf+ii,4)*transcmob(jj,ilocbf+ii)
                enddo
              endif
            enddo
          case default
            write(*,'(" Error! Subroutine Gridraomo supports up to f functions.")')
            call iabort
        end select
      enddo
!
      return
end


!------------------------------------------------
  subroutine gridd2ao(vg2ao,xyzpt,rsqrd)
!------------------------------------------------
!
! Calculate second derivatives of AO values for a grid point
!
      use modmolecule, only : natom, neleca
      use modbasis, only : ex, coeff, nshell, nao, locprim, locbf, locatom, mprim, mbf, mtype
      use modthresh, only : threshex
      implicit none
      integer :: icount, ish, numprim, iatom, iprim, nang, nbf, iloc, i, j
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, two=2.0D+00, three=3.0D+00
      real(8),parameter :: four=4.0D+00, five=5.0D+00, six=6.0D+00, seven=7.0D+00
      real(8),parameter :: sqrt3=1.732050807568877D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrt5=2.236067977499790D+00, sqrt15=3.872983346207417D+00
      real(8),parameter :: facf1=0.36969351199675831D+00  ! 1/sqrt(10-6/sqrt(5))
      real(8),parameter :: facf2=0.86602540378443865D+00  ! 1/sqrt(4/3)
      real(8),parameter :: facf3=0.28116020334310144D+00  ! 1/sqrt(46/3-6/sqrt(5))
      real(8),parameter :: facf4=0.24065403274177409D+00  ! 1/sqrt(28-24/sqrt(5))
      real(8),intent(in) :: xyzpt(3,natom), rsqrd(natom)
      real(8),intent(out) :: vg2ao(nao,6)
      real(8) :: fac, fac2, tmp(10,6), xx, yy, zz, xy, xz, yz, xxx, xxy, xxz, xyy, xyz, xzz
      real(8) :: yyy, yyz, yzz, zzz, xxxx, xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz
      real(8) :: xyzz, xzzz, yyyy, yyyz, yyzz, yzzz, zzzz, xxxxx, yyyyy, zzzzz, xxxxy
      real(8) :: xxxxz, xyyyy, yyyyz, xzzzz, yzzzz, xxxyy, xxxyz, xxxzz, xxyyy, xyyyz
      real(8) :: yyyzz, xxzzz, xyzzz, yyzzz, xxyyz, xxyzz, xyyzz, expval
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
              yy= xyzpt(2,iatom)*xyzpt(2,iatom)
              zz= xyzpt(3,iatom)*xyzpt(3,iatom)
              xy= xyzpt(1,iatom)*xyzpt(2,iatom)
              xz= xyzpt(1,iatom)*xyzpt(3,iatom)
              yz= xyzpt(2,iatom)*xyzpt(3,iatom)
              vg2ao(iloc+1,1)= vg2ao(iloc+1,1)+expval*(fac+fac2*xx)
              vg2ao(iloc+1,2)= vg2ao(iloc+1,2)+expval*(fac+fac2*yy)
              vg2ao(iloc+1,3)= vg2ao(iloc+1,3)+expval*(fac+fac2*zz)
              vg2ao(iloc+1,4)= vg2ao(iloc+1,4)+expval*(    fac2*xy)
              vg2ao(iloc+1,5)= vg2ao(iloc+1,5)+expval*(    fac2*xz)
              vg2ao(iloc+1,6)= vg2ao(iloc+1,6)+expval*(    fac2*yz)
            enddo
          case(1)
            do iprim= 1,numprim
              icount= icount+1
              if(ex(icount)*rsqrd(iatom) > threshex) cycle
              expval= exp(-ex(icount)*rsqrd(iatom))*coeff(icount)
              fac =-ex(icount)*two
              fac2= fac*fac
              xxx= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(1,iatom)*fac2
              xxy= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(2,iatom)*fac2
              xxz= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(3,iatom)*fac2
              xyy= xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)*fac2
              xyz= xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)*fac2
              xzz= xyzpt(1,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac2
              yyy= xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)*fac2
              yyz= xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)*fac2
              yzz= xyzpt(2,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac2
              zzz= xyzpt(3,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac2
              vg2ao(iloc+1,1)= vg2ao(iloc+1,1)+expval*(fac*xyzpt(1,iatom)*three+xxx)
              vg2ao(iloc+2,1)= vg2ao(iloc+2,1)+expval*(fac*xyzpt(2,iatom)      +xxy)
              vg2ao(iloc+3,1)= vg2ao(iloc+3,1)+expval*(fac*xyzpt(3,iatom)      +xxz)
              vg2ao(iloc+1,2)= vg2ao(iloc+1,2)+expval*(fac*xyzpt(1,iatom)      +xyy)
              vg2ao(iloc+2,2)= vg2ao(iloc+2,2)+expval*(fac*xyzpt(2,iatom)*three+yyy)
              vg2ao(iloc+3,2)= vg2ao(iloc+3,2)+expval*(fac*xyzpt(3,iatom)      +yyz)
              vg2ao(iloc+1,3)= vg2ao(iloc+1,3)+expval*(fac*xyzpt(1,iatom)      +xzz)
              vg2ao(iloc+2,3)= vg2ao(iloc+2,3)+expval*(fac*xyzpt(2,iatom)      +yzz)
              vg2ao(iloc+3,3)= vg2ao(iloc+3,3)+expval*(fac*xyzpt(3,iatom)*three+zzz)
              vg2ao(iloc+1,4)= vg2ao(iloc+1,4)+expval*(fac*xyzpt(2,iatom)      +xxy)
              vg2ao(iloc+2,4)= vg2ao(iloc+2,4)+expval*(fac*xyzpt(1,iatom)      +xyy)
              vg2ao(iloc+3,4)= vg2ao(iloc+3,4)+expval*(                        +xyz)
              vg2ao(iloc+1,5)= vg2ao(iloc+1,5)+expval*(fac*xyzpt(3,iatom)      +xxz)
              vg2ao(iloc+2,5)= vg2ao(iloc+2,5)+expval*(                        +xyz)
              vg2ao(iloc+3,5)= vg2ao(iloc+3,5)+expval*(fac*xyzpt(1,iatom)      +xzz)
              vg2ao(iloc+1,6)= vg2ao(iloc+1,6)+expval*(                        +xyz)
              vg2ao(iloc+2,6)= vg2ao(iloc+2,6)+expval*(fac*xyzpt(3,iatom)      +yyz)
              vg2ao(iloc+3,6)= vg2ao(iloc+3,6)+expval*(fac*xyzpt(2,iatom)      +yzz)
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
              xxxx= xx*xx
              xxxy= xx*xy
              xxxz= xx*xz
              xxyy= xx*yy
              xxyz= xx*yz
              xxzz= xx*zz
              xyyy= xy*yy
              xyyz= xy*yz
              xyzz= xy*zz
              xzzz= xz*zz
              yyyy= yy*yy
              yyyz= yy*yz
              yyzz= yy*zz
              yzzz= yz*zz
              zzzz= zz*zz
              tmp(1,1)= expval*(two+xx*five +xxxx)
              tmp(2,1)= expval*(   +yy      +xxyy)
              tmp(3,1)= expval*(   +zz      +xxzz)
              tmp(4,1)= expval*(   +xy*three+xxxy)*sqrt3
              tmp(5,1)= expval*(   +xz*three+xxxz)*sqrt3
              tmp(6,1)= expval*(   +yz      +xxyz)*sqrt3
              tmp(1,2)= expval*(   +xx      +xxyy)
              tmp(2,2)= expval*(two+yy*five +yyyy)
              tmp(3,2)= expval*(   +zz      +yyzz)
              tmp(4,2)= expval*(   +xy*three+xyyy)*sqrt3
              tmp(5,2)= expval*(   +xz      +xyyz)*sqrt3
              tmp(6,2)= expval*(   +yz*three+yyyz)*sqrt3
              tmp(1,3)= expval*(   +xx      +xxzz)
              tmp(2,3)= expval*(   +yy      +yyzz)
              tmp(3,3)= expval*(two+zz*five +zzzz)
              tmp(4,3)= expval*(   +xy      +xyzz)*sqrt3
              tmp(5,3)= expval*(   +xz*three+xzzz)*sqrt3
              tmp(6,3)= expval*(   +yz*three+yzzz)*sqrt3
              tmp(1,4)= expval*(   +xy*two  +xxxy)
              tmp(2,4)= expval*(   +xy*two  +xyyy)
              tmp(3,4)= expval*(            +xyzz)
              tmp(4,4)= expval*(one+xx+yy   +xxyy)*sqrt3
              tmp(5,4)= expval*(   +yz      +xxyz)*sqrt3
              tmp(6,4)= expval*(   +xz      +xyyz)*sqrt3
              tmp(1,5)= expval*(   +xz*two  +xxxz)
              tmp(2,5)= expval*(            +xyyz)
              tmp(3,5)= expval*(   +xz*two  +xzzz)
              tmp(4,5)= expval*(   +yz      +xxyz)*sqrt3
              tmp(5,5)= expval*(one+xx+zz   +xxzz)*sqrt3
              tmp(6,5)= expval*(   +xy      +xyzz)*sqrt3
              tmp(1,6)= expval*(            +xxyz)
              tmp(2,6)= expval*(   +yz*two  +yyyz)
              tmp(3,6)= expval*(   +yz*two  +yzzz)
              tmp(4,6)= expval*(   +xz      +xyyz)*sqrt3
              tmp(5,6)= expval*(   +xy      +xyzz)*sqrt3
              tmp(6,6)= expval*(one+yy+zz   +yyzz)*sqrt3
              if(nbf == 6) then
                do i= 1,6
                  do j= 1,6
                    vg2ao(iloc+j,i)= vg2ao(iloc+j,i)+tmp(j,i)
                  enddo
                enddo
              else
                do i= 1,6
                  vg2ao(iloc+1,i)= vg2ao(iloc+1,i)+tmp(3,i)-(tmp(1,i)+tmp(2,i))*half
                  vg2ao(iloc+2,i)= vg2ao(iloc+2,i)+tmp(5,i)
                  vg2ao(iloc+3,i)= vg2ao(iloc+3,i)+tmp(6,i)
                  vg2ao(iloc+4,i)= vg2ao(iloc+4,i)+(tmp(1,i)-tmp(2,i))*sqrt3h
                  vg2ao(iloc+5,i)= vg2ao(iloc+5,i)+tmp(4,i)
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
              xxx= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(1,iatom)*fac
              xxy= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(2,iatom)*fac
              xxz= xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(3,iatom)*fac
              xyy= xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)*fac
              xyz= xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)*fac
              xzz= xyzpt(1,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac
              yyy= xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)*fac
              yyz= xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)*fac
              yzz= xyzpt(2,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac
              zzz= xyzpt(3,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*fac
              xxxxx= xxx*xx
              yyyyy= yyy*yy
              zzzzz= zzz*zz
              xxxxy= xxx*xy
              xxxxz= xxx*xz
              xyyyy= xyy*yy
              yyyyz= yyy*yz
              xzzzz= xzz*zz
              yzzzz= yzz*zz
              xxxyy= xxx*yy
              xxxyz= xxx*yz
              xxxzz= xxx*zz
              xxyyy= xxy*yy
              xyyyz= xyy*yz
              yyyzz= yyy*zz
              xxzzz= xxz*zz
              xyzzz= xyz*zz
              yyzzz= yyz*zz
              xxyyz= xxy*yz
              xxyzz= xxy*zz
              xyyzz= xyy*zz
              tmp( 1,1)= expval*(xyzpt(1,iatom)*six+xxx*seven+xxxxx)
              tmp( 2,1)= expval*(                  +yyy      +xxyyy)
              tmp( 3,1)= expval*(                  +zzz      +xxzzz)
              tmp( 4,1)= expval*(xyzpt(2,iatom)*two+xxy*five +xxxxy)*sqrt5
              tmp( 5,1)= expval*(xyzpt(3,iatom)*two+xxz*five +xxxxz)*sqrt5
              tmp( 6,1)= expval*(                  +xyy*three+xxxyy)*sqrt5
              tmp( 7,1)= expval*(                  +yyz      +xxyyz)*sqrt5
              tmp( 8,1)= expval*(                  +xzz*three+xxxzz)*sqrt5
              tmp( 9,1)= expval*(                  +yzz      +xxyzz)*sqrt5
              tmp(10,1)= expval*(                  +xyz*three+xxxyz)*sqrt15
              tmp( 1,2)= expval*(                  +xxx      +xxxyy)
              tmp( 2,2)= expval*(xyzpt(2,iatom)*six+yyy*seven+yyyyy)
              tmp( 3,2)= expval*(                  +zzz      +yyzzz)
              tmp( 4,2)= expval*(                  +xxy*three+xxyyy)*sqrt5
              tmp( 5,2)= expval*(                  +xxz      +xxyyz)*sqrt5
              tmp( 6,2)= expval*(xyzpt(1,iatom)*two+xyy*five +xyyyy)*sqrt5
              tmp( 7,2)= expval*(xyzpt(3,iatom)*two+yyz*five +yyyyz)*sqrt5
              tmp( 8,2)= expval*(                  +xzz      +xyyzz)*sqrt5
              tmp( 9,2)= expval*(                  +yzz*three+yyyzz)*sqrt5
              tmp(10,2)= expval*(                  +xyz*three+xyyyz)*sqrt15
              tmp( 1,3)= expval*(                  +xxx      +xxxzz)
              tmp( 2,3)= expval*(                  +yyy      +yyyzz)
              tmp( 3,3)= expval*(xyzpt(3,iatom)*six+zzz*seven+zzzzz)
              tmp( 4,3)= expval*(                  +xxy      +xxyzz)*sqrt5
              tmp( 5,3)= expval*(                  +xxz*three+xxzzz)*sqrt5
              tmp( 6,3)= expval*(                  +xyy      +xyyzz)*sqrt5
              tmp( 7,3)= expval*(                  +yyz*three+yyzzz)*sqrt5
              tmp( 8,3)= expval*(xyzpt(1,iatom)*two+xzz*five +xzzzz)*sqrt5
              tmp( 9,3)= expval*(xyzpt(2,iatom)*two+yzz*five +yzzzz)*sqrt5
              tmp(10,3)= expval*(                  +xyz*three+xyzzz)*sqrt15
              tmp( 1,4)= expval*(                  +xxy*three+xxxxy)
              tmp( 2,4)= expval*(                  +xyy*three+xyyyy)
              tmp( 3,4)= expval*(                            +xyzzz)
              tmp( 4,4)= expval*(xyzpt(1,iatom)*two+xxx+xyy*two+xxxyy)*sqrt5
              tmp( 5,4)= expval*(                  +xyz*two  +xxxyz)*sqrt5
              tmp( 6,4)= expval*(xyzpt(2,iatom)*two+yyy+xxy*two+xxyyy)*sqrt5
              tmp( 7,4)= expval*(                  +xyz*two  +xyyyz)*sqrt5
              tmp( 8,4)= expval*(                  +yzz      +xxyzz)*sqrt5
              tmp( 9,4)= expval*(                  +xzz      +xyyzz)*sqrt5
              tmp(10,4)= expval*(xyzpt(3,iatom)    +xxz+yyz  +xxyyz)*sqrt15
              tmp( 1,5)= expval*(                  +xxz*three+xxxxz)
              tmp( 2,5)= expval*(                            +xyyyz)
              tmp( 3,5)= expval*(                  +xzz*three+xzzzz)
              tmp( 4,5)= expval*(                  +xyz*two  +xxxyz)*sqrt5
              tmp( 5,5)= expval*(xyzpt(1,iatom)*two+xxx+xzz*two+xxxzz)*sqrt5
              tmp( 6,5)= expval*(                  +yyz      +xxyyz)*sqrt5
              tmp( 7,5)= expval*(                  +xyy      +xyyzz)*sqrt5
              tmp( 8,5)= expval*(xyzpt(3,iatom)*two+zzz+xxz*two+xxzzz)*sqrt5
              tmp( 9,5)= expval*(                  +xyz*two  +xyzzz)*sqrt5
              tmp(10,5)= expval*(xyzpt(2,iatom)    +xxy+yzz  +xxyzz)*sqrt15
              tmp( 1,6)= expval*(                            +xxxyz)
              tmp( 2,6)= expval*(                  +yyz*three+yyyyz)
              tmp( 3,6)= expval*(                  +yzz*three+yzzzz)
              tmp( 4,6)= expval*(                  +xxz      +xxyyz)*sqrt5
              tmp( 5,6)= expval*(                  +xxy      +xxyzz)*sqrt5
              tmp( 6,6)= expval*(                  +xyz*two  +xyyyz)*sqrt5
              tmp( 7,6)= expval*(xyzpt(2,iatom)*two+yyy+yzz*two+yyyzz)*sqrt5
              tmp( 8,6)= expval*(                  +xyz*two  +xyzzz)*sqrt5
              tmp( 9,6)= expval*(xyzpt(3,iatom)*two+zzz+yyz*two+yyzzz)*sqrt5
              tmp(10,6)= expval*(xyzpt(1,iatom)    +xyy+xzz  +xyyzz)*sqrt15
              if(nbf == 10) then
                do i= 1,6
                  do j= 1,10
                    vg2ao(iloc+j,i)= vg2ao(iloc+j,i)+tmp(j,i)
                  enddo
                enddo
              else
                do i= 1,6
                  vg2ao(iloc+1,i)= vg2ao(iloc+1,i)+(two*tmp(3,i)-three*(tmp(5,i)+tmp(7,i)))*facf4
                  vg2ao(iloc+2,i)= vg2ao(iloc+2,i)+(-tmp(1,i)-tmp(6,i)+four*tmp(8,i)      )*facf3
                  vg2ao(iloc+3,i)= vg2ao(iloc+3,i)+(-tmp(2,i)-tmp(4,i)+four*tmp(9,i)      )*facf3
                  vg2ao(iloc+4,i)= vg2ao(iloc+4,i)+( tmp(5,i)-tmp(7,i)                    )*facf2
                  vg2ao(iloc+5,i)= vg2ao(iloc+5,i)+  tmp(10,i)
                  vg2ao(iloc+6,i)= vg2ao(iloc+6,i)+( tmp(1,i)-three*tmp(6,i)              )*facf1
                  vg2ao(iloc+7,i)= vg2ao(iloc+7,i)+(-tmp(2,i)+three*tmp(4,i)              )*facf1
                enddo
              endif
            enddo
          case default
            write(*,'(" Error! Subroutine Gridgao supports up to f functions.")')
            call iabort
        end select
      enddo
!
      return
end


!-------------------------------------------------------------------------------------------------
  subroutine gradrexcor(egrad,edftgrad,cmo,fulldmtrx,atomvec,surface,radpt,angpt,rad,ptweight, &
&                       xyzpt,rsqrd,rr,uvec,vao,vmo,dweight,dpa,pa,transcmo,idft,nproc,myrank)
!-------------------------------------------------------------------------------------------------
!
! Driver of derivatives for closed-shell exchange-correlation terms
!
      use modbasis, only : nao
      use modmolecule, only : natom, neleca
      use moddft, only : nrad, nleb
      use modthresh, only : threshweight, threshrho, threshdfock, threshdftao
      implicit none
      integer,intent(in) :: idft, nproc, myrank
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
            call calcexcor(excora,excora,ptenergy,rhoa,rhoa,grhoa,grhoa,one,idft,1)
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
&                       dpa,pa,transcmoa,transcmob,idft,nproc,myrank)
!------------------------------------------------------------------------------------------------
!
! Driver of derivatives for open-shell exchange-correlation terms
!
      use modbasis, only : nao
      use modmolecule, only : natom, neleca, nelecb
      use moddft, only : nrad, nleb
      use modthresh, only : threshweight, threshrho, threshdfock, threshdftao
      implicit none
      integer,intent(in) :: idft, nproc, myrank
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
            call calcexcor(excora,excorb,ptenergy,rhoa,rhob,grhoa,grhob,one,idft,2)
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
      use moddft, only : nrad, nleb
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
&                       +(vg2ao(iao,4)*dra+vgao(iao,1)*drya)*excora(3) &
&                       +(vg2ao(iao,5)*dra+vgao(iao,1)*drza)*excora(4)
            grady= grady+(vg2ao(iao,4)*dra+vgao(iao,2)*drxa)*excora(2) &
&                       +(vg2ao(iao,2)*dra+vgao(iao,2)*drya)*excora(3) &
&                       +(vg2ao(iao,6)*dra+vgao(iao,2)*drza)*excora(4)
            gradz= gradz+(vg2ao(iao,5)*dra+vgao(iao,3)*drxa)*excora(2) &
&                       +(vg2ao(iao,6)*dra+vgao(iao,3)*drya)*excora(3) &
&                       +(vg2ao(iao,3)*dra+vgao(iao,3)*drza)*excora(4)
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
&                       +(vg2ao(iao,4)*dra+vgao(iao,1)*drya)*excora(3) &
&                       +(vg2ao(iao,5)*dra+vgao(iao,1)*drza)*excora(4) &
&                       +(vg2ao(iao,1)*drb+vgao(iao,1)*drxb)*excorb(2) &
&                       +(vg2ao(iao,4)*drb+vgao(iao,1)*dryb)*excorb(3) &
&                       +(vg2ao(iao,5)*drb+vgao(iao,1)*drzb)*excorb(4)
            grady= grady+(vg2ao(iao,4)*dra+vgao(iao,2)*drxa)*excora(2) &
&                       +(vg2ao(iao,2)*dra+vgao(iao,2)*drya)*excora(3) &
&                       +(vg2ao(iao,6)*dra+vgao(iao,2)*drza)*excora(4) &
&                       +(vg2ao(iao,4)*drb+vgao(iao,2)*drxb)*excorb(2) &
&                       +(vg2ao(iao,2)*drb+vgao(iao,2)*dryb)*excorb(3) &
&                       +(vg2ao(iao,6)*drb+vgao(iao,2)*drzb)*excorb(4)
            gradz= gradz+(vg2ao(iao,5)*dra+vgao(iao,3)*drxa)*excora(2) &
&                       +(vg2ao(iao,6)*dra+vgao(iao,3)*drya)*excora(3) &
&                       +(vg2ao(iao,3)*dra+vgao(iao,3)*drza)*excora(4) &
&                       +(vg2ao(iao,5)*drb+vgao(iao,3)*drxb)*excorb(2) &
&                       +(vg2ao(iao,6)*drb+vgao(iao,3)*dryb)*excorb(3) &
&                       +(vg2ao(iao,3)*drb+vgao(iao,3)*drzb)*excorb(4)
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

