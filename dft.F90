!-----------------------------------------------------------------------------------------
  subroutine formfockexcor(fockd,energy,totalelec,cmoa,atomvec,radpt,angpt,rad,ptweight)
!-----------------------------------------------------------------------------------------
!
! Driver of DFT Fock matrix formation from exchange-correlation functionals
!
      use procpar, only : nproc, myrank, MPI_SUM, MPI_COMM_WORLD
      use molecule, only : natom, numatomic, neleca
      use dft, only : nrad, nleb
      use basis, only : nao, nprim
      use thresh, only : threshrho
      implicit none
      integer :: ngridatom, iatom, irad, ileb, icount, ilebstart, jatom, iao, imo
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00
      real(8),intent(in) :: cmoa(nao,nao), atomvec(5,natom,natom), radpt(2,nrad)
      real(8),intent(in) :: angpt(4,nleb), rad(natom), ptweight(nleb,nrad,natom)
      real(8),intent(out) :: fockd(nao*(nao+1)/2), energy, totalelec
      real(8) :: weight, expval(nprim), vao(nao,4), vmo(neleca,4), work(nao)
      real(8) :: rhoa, rhob, grhoa(3), grhob(3), excora(4), xyzpt(3,natom), rsqrd(natom)
      real(8) :: radpoint, tmp(2)
!ishimura
    real(8) :: ee1,ee2,ee3
    common/kazuyadft/ee1,ee2,ee3
    ee1=0.0D+0
    ee2=0.0D+0
    ee3=0.0D+0
!
      energy= zero
      totalelec= zero
      ngridatom= nrad*nleb
!
!$OMP parallel do collapse(2) private(icount,ilebstart,xyzpt,rsqrd,weight, &
!$OMP rhoa,rhob,grhoa,grhob,excora,radpoint,vao,vmo,work,expval) &
!$OMP reduction(+:energy,fockd,totalelec)
      do iatom= 1,natom
        do irad= 1,nrad
          icount=(iatom-1)*ngridatom+(irad-1)*nleb+1+myrank
          ilebstart=mod(icount,nproc)+1
          radpoint= rad(iatom)*radpt(1,irad)
          do ileb= ilebstart,nleb,nproc
!
            do jatom= 1,natom
              xyzpt(1,jatom)= atomvec(1,iatom,jatom)+radpoint*angpt(1,ileb)
              xyzpt(2,jatom)= atomvec(2,iatom,jatom)+radpoint*angpt(2,ileb)
              xyzpt(3,jatom)= atomvec(3,iatom,jatom)+radpoint*angpt(3,ileb)
              rsqrd(jatom)= xyzpt(1,jatom)*xyzpt(1,jatom)+xyzpt(2,jatom)*xyzpt(2,jatom) &
&                          +xyzpt(3,jatom)*xyzpt(3,jatom)
            enddo
!
            weight=ptweight(ileb,irad,iatom)
            call gridao(cmoa,vao,expval,xyzpt,rsqrd)
            call gridgao(cmoa,vao(1,2),expval,xyzpt)
!
            do imo= 1,neleca
              vmo(imo,1)= zero
              vmo(imo,2)= zero
              vmo(imo,3)= zero
              vmo(imo,4)= zero
              do iao= 1,nao
                vmo(imo,1)= vmo(imo,1)+cmoa(iao,imo)*vao(iao,1)
                vmo(imo,2)= vmo(imo,2)+cmoa(iao,imo)*vao(iao,2)
                vmo(imo,3)= vmo(imo,3)+cmoa(iao,imo)*vao(iao,3)
                vmo(imo,4)= vmo(imo,4)+cmoa(iao,imo)*vao(iao,4)
              enddo
            enddo
!
            rhoa= zero
            grhoa(1:3)= zero
            do imo= 1,neleca
              rhoa=     rhoa    +vmo(imo,1)*vmo(imo,1)
              grhoa(1)= grhoa(1)+vmo(imo,2)*vmo(imo,1)
              grhoa(2)= grhoa(2)+vmo(imo,3)*vmo(imo,1)
              grhoa(3)= grhoa(3)+vmo(imo,4)*vmo(imo,1)
            enddo
            rhob= rhoa
            if((rhoa+rhob).le.threshrho)cycle
            grhoa(1)= grhoa(1)*two
            grhoa(2)= grhoa(2)*two
            grhoa(3)= grhoa(3)*two
            grhob(1:3)= grhoa(1:3)
            call calcexcor(excora,energy,rhoa,rhob,grhoa,grhob,weight)
            call fockexcor(fockd,excora,vao,vao(1,2),work,weight,nao)
            totalelec=totalelec+weight*rhoa*2.0D+0
          enddo
        enddo
      enddo
!$OMP end parallel do
!
      tmp(1)= energy
      tmp(2)= totalelec
      call para_allreduce(fockd,rhoa,nao*(nao+1)/2,'D',MPI_SUM,MPI_COMM_WORLD,1)
      call para_allreduce(tmp,rhoa,2,'D',MPI_SUM,MPI_COMM_WORLD,1)
      energy    = tmp(1)
      totalelec = tmp(2)
      return
end


!------------------------------------------------------------------
  subroutine calcexcor(excora,energy,rhoa,rhob,grhoa,grhob,weight)
!------------------------------------------------------------------
!
! Calculate exchange and correlation potential and energy at a grid point
!
      implicit none
      real(8),parameter :: zero=0.0D+00, onethird=0.3333333333333333D+00
      real(8),intent(in) :: rhoa, rhob, grhoa(3), grhob(3), weight
      real(8),intent(out) :: excora(4), energy
      real(8) :: rhoa13, rhob13, csdlda, cb88, cvwn, clyp
!ishimura
    real(8) :: ee1,ee2,ee3,eold
    common/kazuyadft/ee1,ee2,ee3
!
      excora(1:4)= zero
      rhoa13= rhoa**onethird
      rhob13= rhoa13
!
      csdlda= 0.08D+00
      cb88=   0.72D+00
      cvwn=   0.19D+00
      clyp=   0.81D+00
!
!ishimura
  eold=energy
      call funcsdlda(excora,energy,rhoa,rhob,rhoa13,rhob13,weight,csdlda)
  ee1=ee1+energy-eold
  eold=energy
      call funcbecke88(excora,energy,rhoa,rhob,grhoa,grhob,rhoa13,rhob13,weight,cb88)
  ee2=ee2+energy-eold
  eold=energy
      call funcvwn5(excora,energy,rhoa,rhob,rhoa13,rhob13,weight,cvwn)
      call funclyp(excora,energy,rhoa,rhob,grhoa,grhob,rhoa13,rhob13,weight,clyp)
  ee3=ee3+energy-eold
      return
end


!--------------------------------------------------------------
  subroutine fockexcor(fockd,excora,vao,vgao,work,weight,nao)
!--------------------------------------------------------------
!
! ADD exchange-correlation term to Fock matrix at a grid point
!
      use thresh, only : threshdfock
      implicit none
      integer,intent(in) :: nao
      integer :: ij, i, j
      real(8),parameter :: zero=0.0D+00, half=0.5D+00
      real(8),intent(in) :: excora(4), vao(nao), vgao(nao,3), weight
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
        if((vaoabs*gridmax+gridabs*vaomax) >= threshdfock) then
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


!------------------------------------------
  subroutine calcatomvec(atomvec,surface)
!------------------------------------------
!
! Calculate atom vectors and surface shifting parameters
!
! Out : atomvec (atom vector and distance)
!       surface (surface shifting parameters)
!
      use molecule, only : natom, coord, numatomic
      use atominfo, only : atomrad
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
      use iofile, only : iout
      implicit none
      integer,intent(in) :: nleb
      integer :: ileb
      real(8),parameter :: pi4=1.256637061435917D+01
      real(8),intent(out) :: angpt(4,nleb)
      real(8) :: work(nleb,4)
!
      select case(nleb)
        case(6)
          call ld0006(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(14)
          call ld0014(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(26)
          call ld0026(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(38)
          call ld0038(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(50)
          call ld0050(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(74)
          call ld0074(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(86)
          call ld0086(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(110)
          call ld0110(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(146)
          call ld0146(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(170)
          call ld0170(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(194)
          call ld0194(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(230)
          call ld0230(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(266)
          call ld0266(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(302)
          call ld0302(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(350)
          call ld0350(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(434)
          call ld0434(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(590)
          call ld0590(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(770)
          call ld0770(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(974)
          call ld0974(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(1202)
          call ld1202(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(1454)
          call ld1454(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(1730)
          call ld1730(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(2030)
          call ld2030(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(2354)
          call ld2354(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(2702)
          call ld2702(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(3074)
          call ld3074(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(3470)
          call ld3470(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(3890)
          call ld3890(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(4334)
          call ld4334(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(4802)
          call ld4802(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(5294)
          call ld5294(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case(5810)
          call ld5810(work(1,1),work(1,2),work(1,3),work(1,4),nleb)
        case default
          write(iout,'(" Error! Nleb=",i4," is not supported. ")')nleb
          call iabort
      end select
!
      do ileb= 1,nleb
        work(ileb,4)= pi4*work(ileb,4)
      enddo
      angpt=transpose(work)
      return
end


!----------------------------------------------------------------------
  subroutine calcgridweight(ptweight,rad,radpt,angpt,atomvec,surface)
!----------------------------------------------------------------------
!
! Calculate weights of grid points
!
! In  : rad (atom radius)
!       radpt (radial point)
!       angpt (angular point)
!       atomvec (tom vector and distance)
!       surface (surface shifting parameter)
! Out : ptweight  (weight of grid point)
!
      use procpar, only : nproc, myrank
      use dft, only : nrad, nleb
      use molecule, only : natom
      implicit none
      integer :: katom, irad, ileb, iatom, jatom, i, icount, ilebstart, ngridatom
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, oneh=1.5D+00
      real(8),intent(in) :: rad(natom), radpt(2,nrad), angpt(4,nleb), atomvec(5,natom,natom)
      real(8),intent(in) :: surface(natom,natom)
      real(8),intent(out) :: ptweight(nleb,nrad,natom)
      real(8) :: xyzpt(3,natom), work(natom,2)
      real(8) :: radpoint, radweight, wttot, cutij, cutji, xmuij, zmuij, f4, f2
!
      ngridatom= nrad*nleb
!
!$OMP parallel do collapse(2) private(icount,ilebstart,radpoint,radweight,wttot,xyzpt,&
!$OMP work,cutij,cutji,zmuij,xmuij,f4,f2)
      do katom= 1,natom
        do irad= 1,nrad
          icount=(katom-1)*ngridatom+(irad-1)*nleb+1+myrank
          ilebstart=mod(icount,nproc)+1
          radpoint= rad(katom)*radpt(1,irad)
          radweight= rad(katom)*rad(katom)*rad(katom)*radpt(2,irad)
          do ileb= ilebstart,nleb,nproc
            wttot= zero
!
            do iatom= 1,natom
              xyzpt(1,iatom)= atomvec(1,katom,iatom)+radpoint*angpt(1,ileb)
              xyzpt(2,iatom)= atomvec(2,katom,iatom)+radpoint*angpt(2,ileb)
              xyzpt(3,iatom)= atomvec(3,katom,iatom)+radpoint*angpt(3,ileb)
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


!-------------------------------------------------
  subroutine gridao(cmoa,vao,expval,xyzpt,rsqrd)
!-------------------------------------------------
!
! Calculate AO values for a grid point
!
      use molecule, only : natom, neleca
      use basis, only : ex, coeff, nshell, nao, nprim, locprim, locbf, locatom, mprim, mbf, mtype
      use iofile, only : iout
      implicit none
      integer :: icount, ish, numprim, iatom, iprim, nang, nbf, ilocbf, i
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00
      real(8),parameter :: two=2.0D+00, three=3.0D+00, four=4.0D+00
      real(8),parameter :: sqrt3=1.73205080756888D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrt5=2.236067977499790D+00, sqrt15=3.872983346207417D+00
      real(8),parameter :: facf1=0.36969351199675831D+00  ! 1/sqrt(10-6/sqrt(5))
      real(8),parameter :: facf2=0.86602540378443865D+00  ! 1/sqrt(4/3)
      real(8),parameter :: facf3=0.28116020334310144D+00  ! 1/sqrt(46/3-6/sqrt(5))
      real(8),parameter :: facf4=0.24065403274177409D+00  ! 1/sqrt(28-24/sqrt(5))
      real(8),intent(in) :: cmoa(nao,nao), xyzpt(3,natom), rsqrd(natom)
      real(8),intent(out) :: vao(nao), expval(nprim)
      real(8) :: tmp(10)

      icount= 0
      do ish= 1,nshell
        numprim= mprim(ish)
        iatom= locatom(ish)
        do iprim= 1,numprim
          icount= icount+1
          expval(icount)= exp(-ex(icount)*rsqrd(iatom))*coeff(icount)
        enddo
      enddo
!
      vao= zero
      icount= 0
      do ish= 1,nshell
        nang= mtype(ish)
        numprim= mprim(ish)
        nbf = mbf(ish)
        ilocbf= locbf(ish)
        iatom= locatom(ish)
        select case(nang)
          case(0)
            do iprim= 1,numprim
              icount= icount+1
              vao(ilocbf+1)= vao(ilocbf+1)+expval(icount)
            enddo
          case(1)
            do iprim= 1,numprim
              icount= icount+1
              vao(ilocbf+1)= vao(ilocbf+1)+expval(icount)*xyzpt(1,iatom)
              vao(ilocbf+2)= vao(ilocbf+2)+expval(icount)*xyzpt(2,iatom)
              vao(ilocbf+3)= vao(ilocbf+3)+expval(icount)*xyzpt(3,iatom)
            enddo
          case(2)
            do iprim= 1,numprim
              icount= icount+1
              tmp(1)= expval(icount)*xyzpt(1,iatom)*xyzpt(1,iatom)
              tmp(2)= expval(icount)*xyzpt(2,iatom)*xyzpt(2,iatom)
              tmp(3)= expval(icount)*xyzpt(3,iatom)*xyzpt(3,iatom)
              tmp(4)= expval(icount)*xyzpt(1,iatom)*xyzpt(2,iatom)*sqrt3
              tmp(5)= expval(icount)*xyzpt(1,iatom)*xyzpt(3,iatom)*sqrt3
              tmp(6)= expval(icount)*xyzpt(2,iatom)*xyzpt(3,iatom)*sqrt3
              if(nbf == 6) then
                do i= 1,6
                  vao(ilocbf+i)= vao(ilocbf+i)+tmp(i)
                enddo
              else
                vao(ilocbf+1)= vao(ilocbf+1)+tmp(3)-(tmp(1)+tmp(2))*half
                vao(ilocbf+2)= vao(ilocbf+2)+tmp(5)
                vao(ilocbf+3)= vao(ilocbf+3)+tmp(6)
                vao(ilocbf+4)= vao(ilocbf+4)+(tmp(1)-tmp(2))*sqrt3h
                vao(ilocbf+5)= vao(ilocbf+5)+tmp(4)
              endif
            enddo
          case(3)
            do iprim= 1,numprim
              icount= icount+1
              tmp(1)= expval(icount)*xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(1,iatom)
              tmp(2)= expval(icount)*xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)
              tmp(3)= expval(icount)*xyzpt(3,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)
              tmp(4)= expval(icount)*xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(2,iatom)*sqrt5
              tmp(5)= expval(icount)*xyzpt(1,iatom)*xyzpt(1,iatom)*xyzpt(3,iatom)*sqrt5
              tmp(6)= expval(icount)*xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(2,iatom)*sqrt5
              tmp(7)= expval(icount)*xyzpt(2,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)*sqrt5
              tmp(8)= expval(icount)*xyzpt(1,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*sqrt5
              tmp(9)= expval(icount)*xyzpt(2,iatom)*xyzpt(3,iatom)*xyzpt(3,iatom)*sqrt5
              tmp(10)=expval(icount)*xyzpt(1,iatom)*xyzpt(2,iatom)*xyzpt(3,iatom)*sqrt15
              if(nbf == 10) then
                do i= 1,10
                  vao(ilocbf+i)= vao(ilocbf+i)+tmp(i)
                enddo
              else
                vao(ilocbf+1)= vao(ilocbf+1)+( tmp(1)-three*tmp(6)                 )*facf1
                vao(ilocbf+2)= vao(ilocbf+2)+( tmp(5)-tmp(7)                       )*facf2
                vao(ilocbf+3)= vao(ilocbf+3)+(-tmp(1)-tmp(6)+four*tmp(8)           )*facf3
                vao(ilocbf+4)= vao(ilocbf+4)+( two*tmp(3)-three*tmp(5)-three*tmp(7))*facf4
                vao(ilocbf+5)= vao(ilocbf+5)+(-tmp(2)-tmp(4)+four*tmp(9)           )*facf3
                vao(ilocbf+6)= vao(ilocbf+6)+  tmp(10)
                vao(ilocbf+7)= vao(ilocbf+7)+(-tmp(2)+three*tmp(4)                 )*facf1
              endif
            enddo
          case default
            write(iout,'(" Error! Subroutine Gridao supports up to f functions.")')
            call iabort
        end select
      enddo
!
! Transform form AO values to MO values
!
!     call dgemv('T',nao,neleca,one,cmoa,nao,vao,1,zero,vmo,1)
!
      return
end


!---------------------------------------------
  subroutine gridgao(cmoa,vgao,expval,xyzpt)
!---------------------------------------------
!
! Calculate AO gradient values for a grid point
!
      use molecule, only : natom, neleca
      use basis, only : ex, nshell, nao, nprim, locprim, locbf, locatom, mprim, mbf, mtype
      use iofile, only : iout
      implicit none
      integer :: icount, ish, numprim, iatom, iprim, nang, nbf, ilocbf, i, j
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, two=2.0D+00
      real(8),parameter :: three=3.0D+00, four=4.0D+00
      real(8),parameter :: sqrt3=1.732050807568877D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrt5=2.236067977499790D+00, sqrt15=3.872983346207417D+00
      real(8),parameter :: facf1=0.36969351199675831D+00  ! 1/sqrt(10-6/sqrt(5))
      real(8),parameter :: facf2=0.86602540378443865D+00  ! 1/sqrt(4/3)
      real(8),parameter :: facf3=0.28116020334310144D+00  ! 1/sqrt(46/3-6/sqrt(5))
      real(8),parameter :: facf4=0.24065403274177409D+00  ! 1/sqrt(28-24/sqrt(5))
      real(8),intent(in) :: cmoa(nao*nao), expval(nprim), xyzpt(3,natom)
      real(8),intent(out) :: vgao(nao,3)
      real(8) :: fac, tmp(10,3), xx, yy, zz, xy, xz, yz, xxx, xxy, xxz, xyy, xyz, xzz
      real(8) :: yyy, yyz, yzz, zzz, xxxx, xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz
      real(8) :: xyzz, xzzz, yyyy, yyyz, yyzz, yzzz, zzzz
!
      vgao= zero
      icount= 0
      do ish= 1,nshell
        nang= mtype(ish)
        numprim= mprim(ish)
        nbf = mbf(ish)
        ilocbf= locbf(ish)
        iatom= locatom(ish)
        select case(nang)
          case(0)
            do iprim= 1,numprim
              icount= icount+1
              fac= ex(icount)*expval(icount)*two
              vgao(ilocbf+1,1)= vgao(ilocbf+1,1)-fac*xyzpt(1,iatom)
              vgao(ilocbf+1,2)= vgao(ilocbf+1,2)-fac*xyzpt(2,iatom)
              vgao(ilocbf+1,3)= vgao(ilocbf+1,3)-fac*xyzpt(3,iatom)
            enddo
          case(1)
            do iprim= 1,numprim
              icount= icount+1
              fac= ex(icount)*expval(icount)*two
              xx= xyzpt(1,iatom)*xyzpt(1,iatom)
              yy= xyzpt(2,iatom)*xyzpt(2,iatom)
              zz= xyzpt(3,iatom)*xyzpt(3,iatom)
              xy= xyzpt(1,iatom)*xyzpt(2,iatom)
              xz= xyzpt(1,iatom)*xyzpt(3,iatom)
              yz= xyzpt(2,iatom)*xyzpt(3,iatom)
              vgao(ilocbf+1,1)= vgao(ilocbf+1,1)+expval(icount)-fac*xx
              vgao(ilocbf+2,1)= vgao(ilocbf+2,1)               -fac*xy
              vgao(ilocbf+3,1)= vgao(ilocbf+3,1)               -fac*xz
              vgao(ilocbf+1,2)= vgao(ilocbf+1,2)               -fac*xy
              vgao(ilocbf+2,2)= vgao(ilocbf+2,2)+expval(icount)-fac*yy
              vgao(ilocbf+3,2)= vgao(ilocbf+3,2)               -fac*yz
              vgao(ilocbf+1,3)= vgao(ilocbf+1,3)               -fac*xz
              vgao(ilocbf+2,3)= vgao(ilocbf+2,3)               -fac*yz
              vgao(ilocbf+3,3)= vgao(ilocbf+3,3)+expval(icount)-fac*zz
            enddo
          case(2)
            do iprim= 1,numprim
              icount= icount+1
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
              tmp(1,1)= expval(icount)*(xyzpt(1,iatom)*two+xxx)
              tmp(2,1)= expval(icount)*(                  +xyy)
              tmp(3,1)= expval(icount)*(                  +xzz)
              tmp(4,1)= expval(icount)*(xyzpt(2,iatom)    +xxy)*sqrt3
              tmp(5,1)= expval(icount)*(xyzpt(3,iatom)    +xxz)*sqrt3
              tmp(6,1)= expval(icount)*(                  +xyz)*sqrt3
              tmp(1,2)= expval(icount)*(                  +xxy)
              tmp(2,2)= expval(icount)*(xyzpt(2,iatom)*two+yyy)
              tmp(3,2)= expval(icount)*(                  +yzz)
              tmp(4,2)= expval(icount)*(xyzpt(1,iatom)    +xyy)*sqrt3
              tmp(5,2)= expval(icount)*(                  +xyz)*sqrt3
              tmp(6,2)= expval(icount)*(xyzpt(3,iatom)    +yyz)*sqrt3
              tmp(1,3)= expval(icount)*(                  +xxz)
              tmp(2,3)= expval(icount)*(                  +yyz)
              tmp(3,3)= expval(icount)*(xyzpt(3,iatom)*two+zzz)
              tmp(4,3)= expval(icount)*(                  +xyz)*sqrt3
              tmp(5,3)= expval(icount)*(xyzpt(1,iatom)    +xzz)*sqrt3
              tmp(6,3)= expval(icount)*(xyzpt(2,iatom)    +yzz)*sqrt3
              if(nbf == 6) then
                do i= 1,3
                  do j= 1,6
                    vgao(ilocbf+j,i)= vgao(ilocbf+j,i)+tmp(j,i)
                  enddo
                enddo
              else
                do i= 1,3
                  vgao(ilocbf+1,i)= vgao(ilocbf+1,i)+tmp(3,i)-(tmp(1,i)+tmp(2,i))*half
                  vgao(ilocbf+2,i)= vgao(ilocbf+2,i)+tmp(5,i)
                  vgao(ilocbf+3,i)= vgao(ilocbf+3,i)+tmp(6,i)
                  vgao(ilocbf+4,i)= vgao(ilocbf+4,i)+(tmp(1,i)-tmp(2,i))*sqrt3h
                  vgao(ilocbf+5,i)= vgao(ilocbf+5,i)+tmp(4,i)
                enddo
              endif
            enddo
          case(3)
            do iprim= 1,numprim
              icount= icount+1
              fac= ex(icount)*two
              xx= xyzpt(1,iatom)*xyzpt(1,iatom)
              yy= xyzpt(2,iatom)*xyzpt(2,iatom)
              zz= xyzpt(3,iatom)*xyzpt(3,iatom)
              xy= xyzpt(1,iatom)*xyzpt(2,iatom)
              xz= xyzpt(1,iatom)*xyzpt(3,iatom)
              yz= xyzpt(2,iatom)*xyzpt(3,iatom)
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
              tmp( 1,1)= expval(icount)*(xx*three+xxxx)
              tmp( 2,1)= expval(icount)*(        +xyyy)
              tmp( 3,1)= expval(icount)*(        +xzzz)
              tmp( 4,1)= expval(icount)*(xy*two  +xxxy)*sqrt5
              tmp( 5,1)= expval(icount)*(xz*two  +xxxz)*sqrt5
              tmp( 6,1)= expval(icount)*(yy      +xxyy)*sqrt5
              tmp( 7,1)= expval(icount)*(        +xyyz)*sqrt5
              tmp( 8,1)= expval(icount)*(zz      +xxzz)*sqrt5
              tmp( 9,1)= expval(icount)*(        +xyzz)*sqrt5
              tmp(10,1)= expval(icount)*(yz      +xxyz)*sqrt15
              tmp( 1,2)= expval(icount)*(        +xxxy)
              tmp( 2,2)= expval(icount)*(yy*three+yyyy)
              tmp( 3,2)= expval(icount)*(        +yzzz)
              tmp( 4,2)= expval(icount)*(xx      +xxyy)*sqrt5
              tmp( 5,2)= expval(icount)*(        +xxyz)*sqrt5
              tmp( 6,2)= expval(icount)*(xy*two  +xyyy)*sqrt5
              tmp( 7,2)= expval(icount)*(yz*two  +yyyz)*sqrt5
              tmp( 8,2)= expval(icount)*(        +xyzz)*sqrt5
              tmp( 9,2)= expval(icount)*(zz      +yyzz)*sqrt5
              tmp(10,2)= expval(icount)*(xz      +xyyz)*sqrt15
              tmp( 1,3)= expval(icount)*(        +xxxz)
              tmp( 2,3)= expval(icount)*(        +yyyz)
              tmp( 3,3)= expval(icount)*(zz*three+zzzz)
              tmp( 4,3)= expval(icount)*(        +xxyz)*sqrt5
              tmp( 5,3)= expval(icount)*(xx      +xxzz)*sqrt5
              tmp( 6,3)= expval(icount)*(        +xyyz)*sqrt5
              tmp( 7,3)= expval(icount)*(yy      +yyzz)*sqrt5
              tmp( 8,3)= expval(icount)*(xz*two  +xzzz)*sqrt5
              tmp( 9,3)= expval(icount)*(yz*two  +yzzz)*sqrt5
              tmp(10,3)= expval(icount)*(xy      +xyzz)*sqrt15
              if(nbf == 10) then
                do i= 1,3
                  do j= 1,10
                    vgao(ilocbf+j,i)= vgao(ilocbf+j,i)+tmp(j,i)
                  enddo
                enddo
              else
                do i= 1,3
                  vgao(ilocbf+1,i)= vgao(ilocbf+1,i)+( tmp(1,i)-three*tmp(6,i)              )*facf1
                  vgao(ilocbf+2,i)= vgao(ilocbf+2,i)+( tmp(5,i)-tmp(7,i)                    )*facf2
                  vgao(ilocbf+3,i)= vgao(ilocbf+3,i)+(-tmp(1,i)-tmp(6,i)+four*tmp(8,i)      )*facf3
                  vgao(ilocbf+4,i)= vgao(ilocbf+4,i)+(two*tmp(3,i)-three*(tmp(5,i)+tmp(7,i)))*facf4
                  vgao(ilocbf+5,i)= vgao(ilocbf+5,i)+(-tmp(2,i)-tmp(4,i)+four*tmp(9,i)      )*facf3
                  vgao(ilocbf+6,i)= vgao(ilocbf+6,i)+  tmp(10,i)
                  vgao(ilocbf+7,i)= vgao(ilocbf+7,i)+(-tmp(2,i)+three*tmp(4,i)              )*facf1
                enddo
              endif
            enddo
          case default
            write(iout,'(" Error! Subroutine Gridgao supports up to f functions.")')
            call iabort
        end select
      enddo
!
! Transform form AO values to MO values
!
!     call dgemm('T','N',neleca,3,nao,one,cmoa,nao,vgao,nao,zero,vgmo,neleca)
      return
end





