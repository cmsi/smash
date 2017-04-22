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
!--------------------------------------------------------------------
  subroutine gradoneei(egrad,egrad1,fulldmtrx,ewdmtrx,nproc,myrank)
!--------------------------------------------------------------------
!
! Driver of derivatives for one-electron and overlap integrals
!
      use modparallel, only : master
      use modbasis, only : nshell, nao, mtype
      use modmolecule, only : natom
      use modecp, only : flagecp
      implicit none
      integer,intent(in) :: nproc, myrank
      integer :: ish, jsh, len1, i, maxfunc(0:6), maxbasis
      real(8),parameter :: zero=0.0D+00, two=2.0D+00
      real(8),intent(in) :: fulldmtrx(nao*nao), ewdmtrx(nao*(nao+1)/2)
      real(8),intent(out) :: egrad1(3*natom)
      real(8),intent(inout) :: egrad(3*natom)
      data maxfunc/1,3,6,10,15,21,28/
!
      maxbasis= maxval(mtype(1:nshell))
      if(maxbasis > 6) then
        if(master) write(*,'(" Error! This program supports up to h function in gradoneei")')
        call iabort
      endif
      len1= maxfunc(maxbasis+1)
      egrad1(:)= zero
!
!$OMP parallel reduction(+:egrad1)
      do ish= nshell-myrank,1,-nproc
!$OMP do
        do jsh= 1,ish
          call calcdoverlap(egrad1,ewdmtrx,ish,jsh)
          call calchelfey(egrad1,fulldmtrx,ish,jsh)
        enddo
!$OMP enddo
      enddo
!
      do ish= myrank+1,nshell,nproc
!$OMP do
        do jsh= 1,nshell
          call calcdkinetic(egrad1,fulldmtrx,ish,jsh)
          call calcdcoulomb(egrad1,fulldmtrx,ish,jsh,len1)
        enddo
!$OMP enddo
      enddo
!$OMP end parallel
!
      do i= 1,3*natom
        egrad(i)= egrad(i)+egrad1(i)*two
      enddo
!
! Add ECP derivative terms
!
      if(flagecp) then
        egrad1(:)= zero
        call gradoneeiecp(egrad1,fulldmtrx,nproc,myrank)
        do i= 1,3*natom
          egrad(i)= egrad(i)+egrad1(i)*two
        enddo
      endif
! 
      return
end


!--------------------------------------------------------------------
  subroutine calcewdmtrx(cmo,energymo,fulldmtrx,ewdmtrx,ndim,nelec)
!--------------------------------------------------------------------
!
! Calculate energy-weighted and normal density matrix
!
! In  : cmo      (MO coefficient matrix)
!       energymo (MO energy)
!       ndim     (Dimension of basis functions)
!       nelec    (Number of electrons)
! Out : fulldmtrx(Full density matrix)
!       ewdmtrx  (Upper-triangle energy-weighted density matrix)
!
      implicit none
      integer,intent(in) :: ndim, nelec
      integer :: i, j, k, ij
      real(8),parameter :: zero=0.0D+00, two=2.0D+00
      real(8),intent(in) :: cmo(ndim,ndim), energymo(ndim)
      real(8),intent(out) :: fulldmtrx(ndim,ndim), ewdmtrx(ndim*(ndim+1)/2)
!
      fulldmtrx= transpose(cmo)
!$OMP parallel do schedule(guided) private(ij)
      do i= ndim,1,-1
        ij= i*(i-1)/2
        do j= 1,i
          ewdmtrx(ij+j)= zero
          do k= 1,nelec
            ewdmtrx(ij+j)= ewdmtrx(ij+j)-fulldmtrx(k,i)*fulldmtrx(k,j)*energymo(k)
          enddo
          ewdmtrx(ij+j)= ewdmtrx(ij+j)*two
        enddo
      enddo
!$OMP end parallel do
!
      call dgemm('N','T',ndim,ndim,nelec,two,cmo,ndim,cmo,ndim,zero,fulldmtrx,ndim)
      return
end


!-----------------------------------------------------------------------------------------
  subroutine calcuewdmtrx(cmoa,cmob,energymoa,energymob,fulldmtrx1,fulldmtrx2,ewdmtrx, &
&                         ndim,neleca,nelecb)
!-----------------------------------------------------------------------------------------
!
! Calculate energy-weighted and normal density matrix for open-shell
!
! In  : cmoa     (Alpha MO coefficient matrix)
!       cmob     (Beta MO coefficient matrix)
!       energymoa(Alpha MO energy)
!       energymob(Beta MO energy)
!       ndim     (Dimension of basis functions)
!       neleca   (Number of alpha electrons)
!       nelecb   (Number of beta electrons)
! Out : fulldmtrx1(Full alpha+beta density matrix)
!       fulldmtrx2(Full alpha-beta density matrix)
!       ewdmtrx   (Upper-triangle energy-weighted density matrix)
!
      implicit none
      integer,intent(in) :: ndim, neleca, nelecb
      integer :: i, j, k, ij
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: cmoa(ndim,ndim), cmob(ndim,ndim), energymoa(ndim), energymob(ndim)
      real(8),intent(out) :: fulldmtrx1(ndim,ndim), fulldmtrx2(ndim,ndim)
      real(8),intent(out) :: ewdmtrx(ndim*(ndim+1)/2)
      real(8) :: dena, denb
!
      fulldmtrx1= transpose(cmoa)
      fulldmtrx2= transpose(cmob)
!$OMP parallel do schedule(guided) private(ij)
      do i= ndim,1,-1
        ij= i*(i-1)/2
        do j= 1,i
          ewdmtrx(ij+j)= zero
          do k= 1,neleca
            ewdmtrx(ij+j)= ewdmtrx(ij+j)-fulldmtrx1(k,i)*fulldmtrx1(k,j)*energymoa(k)
          enddo
          do k= 1,nelecb
            ewdmtrx(ij+j)= ewdmtrx(ij+j)-fulldmtrx2(k,i)*fulldmtrx2(k,j)*energymob(k)
          enddo
        enddo
      enddo
!$OMP end parallel do
!
      call dgemm('N','T',ndim,ndim,neleca,one,cmoa,ndim,cmoa,ndim,zero,fulldmtrx1,ndim)
      call dgemm('N','T',ndim,ndim,nelecb,one,cmob,ndim,cmob,ndim,zero,fulldmtrx2,ndim)
!
!$OMP parallel do private(dena,denb)
      do i= 1,ndim
        do j= 1,ndim
          dena= fulldmtrx1(j,i)+fulldmtrx2(j,i)
          denb= fulldmtrx1(j,i)-fulldmtrx2(j,i)
          fulldmtrx1(j,i)= dena
          fulldmtrx2(j,i)= denb
        enddo
      enddo
!$OMP end parallel do
!
      return
end


!-------------------------------------------------
  subroutine calcdoverlap(egrad,ewdmtrx,ish,jsh)
!-------------------------------------------------
!
! Driver of overlap derivative term
!
! In : ewdmtrx  (Energy-weighted density matrix)
!    : ish, jsh (Shell indices)
! Inout : egrad (Energy gradient value)
!
      use modparam, only : mxprsh
      use modbasis, only : locatom, locprim, locbf, mprim, mbf, mtype, ex, coeff, nao
      use modthresh, only : threshex
      use modmolecule, only : natom, coord
      implicit none
      integer,intent(in) :: ish, jsh
      integer :: iatom, jatom, iloc, jloc, ilocbf, jlocbf, nprimi, nprimj, nangi, nangj
      integer :: nbfi, nbfj, iprim, jprim, ncarti, ncartj, i, j, iang, jang, ii, ij
      integer :: ix, jx, iy, jy, iz, jz, ncart(0:6)
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00
      real(8),intent(in) :: ewdmtrx(nao*(nao+1)/2)
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: xyzij(3), rij, rij2, fac, exi, exj, exj2, ci, cj
      real(8) :: ex1, ex2, ex3, xyzpij(3,2), cij
      real(8) :: xyzint(3), sx(0:7,0:6,2), sy(0:7,0:6,2), sz(0:7,0:6,2)
      real(8) :: dsint(28,28,3)
      data ncart/1,3,6,10,15,21,28/
!
      iatom = locatom(ish)
      jatom = locatom(jsh)
      if(iatom == jatom) return
!
      iloc  = locprim(ish)
      ilocbf= locbf(ish)
      nprimi= mprim(ish)
      nangi = mtype(ish)
      nbfi  = mbf(ish)
      jloc  = locprim(jsh)
      jlocbf= locbf(jsh)
      nprimj= mprim(jsh)
      nangj = mtype(jsh)
      nbfj  = mbf(jsh)
!
      if((nangi > 6).or.(nangj > 6)) then
        write(*,'(" Error! This program supports up to h function in calcdoverlap")')
        call iabort
      endif
!
      do i= 1,3
        xyzij(i)= coord(i,iatom)-coord(i,jatom)
      enddo
      rij= xyzij(1)*xyzij(1)+xyzij(2)*xyzij(2)+xyzij(3)*xyzij(3)
      ncarti= ncart(nangi)
      ncartj= ncart(nangj)
!
      do i= 1,ncarti
        do j= 1,ncartj
          dsint(j,i,1)= zero
          dsint(j,i,2)= zero
          dsint(j,i,3)= zero
        enddo
      enddo
!
! Calculate overlap derivative for each primitive
!
      do iprim= 1,nprimi
        exi= ex(iloc+iprim)
        ci = coeff(iloc+iprim)
        do jprim= 1,nprimj
          exj= ex(jloc+jprim)
          ex1= exi+exj
          ex2= one/ex1
          rij2=rij*exi*exj*ex2
          if(rij2 > threshex) cycle
          exj2= exj*two
          ex3= sqrt(ex2)
          fac= exp(-rij2)   !*ex2*ex3
          do i= 1,3
            xyzpij(i,1)=-exj*xyzij(i)*ex2
            xyzpij(i,2)= exi*xyzij(i)*ex2
          enddo
          cj = coeff(jloc+jprim)*fac
!
          do iang= 0,nangi
            do jang= 0,nangj+1
              call ghquad(xyzint,ex3,xyzpij,iang,jang)
              sx(jang,iang,1)= xyzint(1)*ex3
              sy(jang,iang,1)= xyzint(2)*ex3
              sz(jang,iang,1)= xyzint(3)*ex3
            enddo
          enddo
          do iang= 0,nangi
            sx(0,iang,2)= sx(1,iang,1)*exj2
            sy(0,iang,2)= sy(1,iang,1)*exj2
            sz(0,iang,2)= sz(1,iang,1)*exj2
          enddo
          do iang= 0,nangi
            do jang= 1,nangj
              sx(jang,iang,2)= sx(jang+1,iang,1)*exj2-sx(jang-1,iang,1)*jang
              sy(jang,iang,2)= sy(jang+1,iang,1)*exj2-sy(jang-1,iang,1)*jang
              sz(jang,iang,2)= sz(jang+1,iang,1)*exj2-sz(jang-1,iang,1)*jang
            enddo
          enddo
          cij= ci*cj
          i= 0
          do ix= nangi,0,-1
            do iy= nangi-ix,0,-1
              iz= nangi-ix-iy
              i= i+1
              j= 0
              do jx= nangj,0,-1
                do jy= nangj-jx,0,-1
                  jz= nangj-jx-jy
                  j= j+1
                  dsint(j,i,1)= dsint(j,i,1)+cij*sx(jx,ix,2)*sy(jy,iy,1)*sz(jz,iz,1)
                  dsint(j,i,2)= dsint(j,i,2)+cij*sx(jx,ix,1)*sy(jy,iy,2)*sz(jz,iz,1)
                  dsint(j,i,3)= dsint(j,i,3)+cij*sx(jx,ix,1)*sy(jy,iy,1)*sz(jz,iz,2)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
!
      if((nbfi >= 5).or.(nbfj >= 5)) then
        call nrmlz1(dsint(1,1,1),nbfi,nbfj,ncarti)
        call nrmlz1(dsint(1,1,2),nbfi,nbfj,ncarti)
        call nrmlz1(dsint(1,1,3),nbfi,nbfj,ncarti)
      endif
!
      do i= 1,nbfi
        ii= ilocbf+i
        ij= ii*(ii-1)/2+jlocbf
        do j= 1,nbfj
          egrad(1,iatom)= egrad(1,iatom)-ewdmtrx(ij+j)*dsint(j,i,1)
          egrad(2,iatom)= egrad(2,iatom)-ewdmtrx(ij+j)*dsint(j,i,2)
          egrad(3,iatom)= egrad(3,iatom)-ewdmtrx(ij+j)*dsint(j,i,3)
          egrad(1,jatom)= egrad(1,jatom)+ewdmtrx(ij+j)*dsint(j,i,1)
          egrad(2,jatom)= egrad(2,jatom)+ewdmtrx(ij+j)*dsint(j,i,2)
          egrad(3,jatom)= egrad(3,jatom)+ewdmtrx(ij+j)*dsint(j,i,3)
        enddo
      enddo
!
      return
end


!---------------------------------------------------
  subroutine calcdkinetic(egrad,fulldmtrx,ish,jsh)
!---------------------------------------------------
!
! Driver of kinetic derivative term
!
! In : fulldmtrx (Density matrix)
!    : ish, jsh  (Shell indices)
! Inout : egrad  (Energy gradient value)
!
      use modparam, only : mxprsh
      use modbasis, only : locatom, locprim, locbf, mprim, mbf, mtype, ex, coeff, nao
      use modthresh, only : threshex
      use modmolecule, only : natom, coord
      implicit none
      integer,intent(in) :: ish, jsh
      integer :: iatom, jatom, iloc, jloc, ilocbf, jlocbf, nprimi, nprimj, nangi, nangj
      integer :: nbfi, nbfj, iprim, jprim, ncarti, ncartj, i, j, iang, jang, ii
      integer :: ix, jx, iy, jy, iz, jz, ncart(0:6)
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00, half=0.5D+00
      real(8),intent(in) :: fulldmtrx(nao,nao)
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: xyzij(3), rij, rij2, fac, exi, exi2, exj, exj2, ci, cj
      real(8) :: ex1, ex2, ex3, xyzpij(3,2), cij
      real(8) :: xyzint(3), sx(0:7,0:8,4), sy(0:7,0:8,4), sz(0:7,0:8,4)
      real(8) :: dtint(28,28,3)
      data ncart/1,3,6,10,15,21,28/
!
      iatom = locatom(ish)
      iloc  = locprim(ish)
      ilocbf= locbf(ish)
      nprimi= mprim(ish)
      nangi = mtype(ish)
      nbfi  = mbf(ish)
      jatom = locatom(jsh)
      jloc  = locprim(jsh)
      jlocbf= locbf(jsh)
      nprimj= mprim(jsh)
      nangj = mtype(jsh)
      nbfj  = mbf(jsh)
!
      if((nangi > 6).or.(nangj > 6))then
        write(*,'(" Error! This program supports up to h function in calcdkinetic.")')
        call iabort
      endif
!
      do i= 1,3
        xyzij(i)= coord(i,iatom)-coord(i,jatom)
      enddo
      rij= xyzij(1)*xyzij(1)+xyzij(2)*xyzij(2)+xyzij(3)*xyzij(3)
      ncarti= ncart(nangi)
      ncartj= ncart(nangj)
!
      do i= 1,ncarti
        do j= 1,ncartj
          dtint(j,i,1)= zero
          dtint(j,i,2)= zero
          dtint(j,i,3)= zero
        enddo
      enddo
!
! Calculate overlap derivative for each primitive
!
      do iprim= 1,nprimi
        exi= ex(iloc+iprim)
        ci = coeff(iloc+iprim)
        exi2=-two*exi*exi
        do jprim= 1,nprimj
          exj= ex(jloc+jprim)
          ex1= exi+exj
          ex2= one/ex1
          rij2=rij*exi*exj*ex2
          if(rij2 > threshex) cycle
          ex3= sqrt(ex2)
          exj2= exj*two
          fac= exp(-rij2)   !*ex2*ex3
          do i= 1,3
            xyzpij(i,1)=-exj*xyzij(i)*ex2
            xyzpij(i,2)= exi*xyzij(i)*ex2
          enddo
          cj = coeff(jloc+jprim)*fac
!
! overlap term
!
          do iang= 0,nangi+2
            do jang= 0,nangj+1
              call ghquad(xyzint,ex3,xyzpij,iang,jang)
              sx(jang,iang,1)= xyzint(1)*ex3
              sy(jang,iang,1)= xyzint(2)*ex3
              sz(jang,iang,1)= xyzint(3)*ex3
            enddo
          enddo
!
! kinetic term
!
          do iang= 0,nangi
            do jang= 0,nangj+1
              sx(jang,iang,2)= sx(jang,iang+2,1)*exi2+sx(jang,iang,1)*(two*iang+one)*exi
              sy(jang,iang,2)= sy(jang,iang+2,1)*exi2+sy(jang,iang,1)*(two*iang+one)*exi
              sz(jang,iang,2)= sz(jang,iang+2,1)*exi2+sz(jang,iang,1)*(two*iang+one)*exi
              if(iang >= 2) then
                sx(jang,iang,2)=sx(jang,iang,2)-sx(jang,iang-2,1)*half*iang*(iang-1)
                sy(jang,iang,2)=sy(jang,iang,2)-sy(jang,iang-2,1)*half*iang*(iang-1)
                sz(jang,iang,2)=sz(jang,iang,2)-sz(jang,iang-2,1)*half*iang*(iang-1)
              endif
            enddo
          enddo
!
! overlap derivative term
!
          do iang= 0,nangi
            sx(0,iang,3)= sx(1,iang,1)*exj2
            sy(0,iang,3)= sy(1,iang,1)*exj2
            sz(0,iang,3)= sz(1,iang,1)*exj2
          enddo
          do iang= 0,nangi
            do jang= 1,nangj
              sx(jang,iang,3)= sx(jang+1,iang,1)*exj2-sx(jang-1,iang,1)*jang
              sy(jang,iang,3)= sy(jang+1,iang,1)*exj2-sy(jang-1,iang,1)*jang
              sz(jang,iang,3)= sz(jang+1,iang,1)*exj2-sz(jang-1,iang,1)*jang
            enddo
          enddo
!
! kinetic derivative term
!
          do iang= 0,nangi
            sx(0,iang,4)= sx(1,iang,2)*exj2
            sy(0,iang,4)= sy(1,iang,2)*exj2
            sz(0,iang,4)= sz(1,iang,2)*exj2
          enddo
          do iang= 0,nangi
            do jang= 1,nangj
              sx(jang,iang,4)= sx(jang+1,iang,2)*exj2-sx(jang-1,iang,2)*jang
              sy(jang,iang,4)= sy(jang+1,iang,2)*exj2-sy(jang-1,iang,2)*jang
              sz(jang,iang,4)= sz(jang+1,iang,2)*exj2-sz(jang-1,iang,2)*jang
            enddo
          enddo
!
          cij= ci*cj
          i= 0
          do ix= nangi,0,-1
            do iy= nangi-ix,0,-1
              iz= nangi-ix-iy
              i= i+1
              j= 0
              do jx= nangj,0,-1
                do jy= nangj-jx,0,-1
                  jz= nangj-jx-jy
                  j= j+1
                  dtint(j,i,1)= dtint(j,i,1)+cij*(sx(jx,ix,4)*sy(jy,iy,1)*sz(jz,iz,1) &
&                                                +sx(jx,ix,3)*sy(jy,iy,2)*sz(jz,iz,1) &
&                                                +sx(jx,ix,3)*sy(jy,iy,1)*sz(jz,iz,2)) 
                  dtint(j,i,2)= dtint(j,i,2)+cij*(sx(jx,ix,2)*sy(jy,iy,3)*sz(jz,iz,1) &
&                                                +sx(jx,ix,1)*sy(jy,iy,4)*sz(jz,iz,1) &
&                                                +sx(jx,ix,1)*sy(jy,iy,3)*sz(jz,iz,2)) 
                  dtint(j,i,3)= dtint(j,i,3)+cij*(sx(jx,ix,2)*sy(jy,iy,1)*sz(jz,iz,3) &
&                                                +sx(jx,ix,1)*sy(jy,iy,2)*sz(jz,iz,3) &
&                                                +sx(jx,ix,1)*sy(jy,iy,1)*sz(jz,iz,4)) 
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
!
      if((nbfi >= 5).or.(nbfj >= 5)) then
        call nrmlz1(dtint(1,1,1),nbfi,nbfj,ncarti)
        call nrmlz1(dtint(1,1,2),nbfi,nbfj,ncarti)
        call nrmlz1(dtint(1,1,3),nbfi,nbfj,ncarti)
      endif
!
      do i= 1,nbfi
        ii= ilocbf+i
        do j= 1,nbfj
          egrad(1,jatom)= egrad(1,jatom)+fulldmtrx(jlocbf+j,ii)*dtint(j,i,1)
          egrad(2,jatom)= egrad(2,jatom)+fulldmtrx(jlocbf+j,ii)*dtint(j,i,2)
          egrad(3,jatom)= egrad(3,jatom)+fulldmtrx(jlocbf+j,ii)*dtint(j,i,3)
        enddo
      enddo
!
      return
end


!--------------------------------------------------------
  subroutine calcdcoulomb(egrad,fulldmtrx,ish,jsh,len1)
!--------------------------------------------------------
!
! Driver of Coulomb derivative term
!
! In : fulldmtrx (density matrix)
!    : ish, jsh  (shell indices)
! Inout : egrad  (energy gradient value)
!
      use modparam, only : mxprsh
      use modbasis, only : locatom, locprim, locbf, mprim, mbf, mtype, ex, coeff, nao
      use modthresh, only : threshex
      use modmolecule, only : natom, coord, znuc
      implicit none
      integer,intent(in) :: ish, jsh, len1
      integer :: nangij(2), nprimij(2), nbfij(2), iatom, jatom, iloc, jloc, ilocbf, jlocbf
      integer :: iprim, jprim, i, j, k, ii, ncart(0:6)
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00, half=0.5D+00, three=3.0D+00
      real(8),parameter :: four=4.0D+00, five=5.0D+00, six=6.0D+00, eight=8.0D+00, ten=10.0D+00
      real(8),parameter :: twelve=12.0D+00, p15=15.0D+00, p24=24.0D+00, p30=30.0D+00
      real(8),parameter :: p40=40.0D+00
      real(8),parameter :: third=3.333333333333333D-01, eighth=0.125D+00
      real(8),parameter :: sqrt3=1.732050807568877D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrtthird=0.5773502691896258D+00, sqrtfifth=0.4472135954999579D+00
      real(8),parameter :: sqrt3fifth=0.7745966692414834D+00, sqrtseventh=0.3779644730092272D+00
      real(8),parameter :: sqrtinv35=0.1690308509457033D+00, sqrt3inv35=0.2927700218845599D+00
      real(8),parameter :: sqrt5=2.236067977499790D+00, sqrt15=3.872983346207417D+00
      real(8),parameter :: sqrt7=2.645751311064591D+00, sqrt35=5.916079783099616D+00
      real(8),parameter :: sqrt35third=3.415650255319866D+00
      real(8),parameter :: sqrtinv15=2.581988897471611D-01, sqrtinv21=2.182178902359924D-01
      real(8),parameter :: sqrtinv63=1.259881576697424D-01, sqrtinv105=9.759000729485332D-02
      real(8),parameter :: sqrtinv11=3.015113445777636D-01, sqrtinv33=1.740776559556978D-01
      real(8),parameter :: sqrtinv99=1.005037815259212D-01, sqrtinv231=6.579516949597690D-02
      real(8),parameter :: sqrtinv385=5.096471914376255D-02, sqrt5inv231=1.471224715841249D-01
      real(8),parameter :: sqrt21=4.582575694955840D+00, sqrt63=7.937253933193772D+00
      real(8),parameter :: sqrt105=1.024695076595960D+01
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
      real(8),intent(in) :: fulldmtrx(nao,nao)
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: exij(mxprsh,2), cij(mxprsh,2), coordij(3,2)
      real(8) :: cint1(len1,len1), cint2(len1,len1), dcint(28,28,3), work(21)
      data ncart/1,3,6,10,15,21,28/
!
      nangij(1)= mtype(ish)
      nangij(2)= mtype(jsh)
      nprimij(1)= mprim(ish)
      nprimij(2)= mprim(jsh)
      nbfij(1)  = mbf(ish)
      nbfij(2)  = mbf(jsh)
      iatom = locatom(ish)
      iloc  = locprim(ish)
      ilocbf= locbf(ish)
      jatom = locatom(jsh)
      jloc  = locprim(jsh)
      jlocbf= locbf(jsh)
!
      do i= 1,3
        coordij(i,1)= coord(i,iatom)
        coordij(i,2)= coord(i,jatom)
      enddo
      do iprim= 1,nprimij(1)
        exij(iprim,1)= ex(iloc+iprim)
        cij(iprim,1) = coeff(iloc+iprim)
      enddo
!
      nangij(2)= nangij(2)+1
      nbfij(2)= ncart(nangij(2))
      do jprim= 1,nprimij(2)
        exij(jprim,2)= ex(jloc+jprim)
        cij(jprim,2) = coeff(jloc+jprim)*two*ex(jloc+jprim)
      enddo
!
      if((nangij(1) <= 2).and.(nangij(2) <= 2)) then
        call int1cmd(cint1,exij,cij,coordij,coord,znuc,natom, &
&                    nprimij,nangij,nbfij,len1,mxprsh,threshex)
      else
        if((nangij(1) > 6).or.(nangij(2) > 7))then
          write(*,'(" Error! This program supports up to h function in int1c")')
          call iabort
        endif
        call int1rys(cint1,exij,cij,coordij,coord,znuc,natom, &
&                    nprimij,nangij,nbfij,len1,mxprsh,threshex)
      endif
!
      if(mtype(jsh) >= 1) then
        nangij(2)= mtype(jsh)-1
        nbfij(2)= ncart(nangij(2))
        do jprim= 1,nprimij(2)
          cij(jprim,2) = coeff(jloc+jprim)
        enddo
        if((nangij(1) <= 2).and.(nangij(2) <= 2)) then
          call int1cmd(cint2,exij,cij,coordij,coord,znuc,natom, &
&                      nprimij,nangij,nbfij,len1,mxprsh,threshex)
        else
          if((nangij(1) > 6).or.(nangij(2) > 5))then
            write(*,'(" Error! This program supports up to h function in int1c")')
            call iabort
          endif
!
          call int1rys(cint2,exij,cij,coordij,coord,znuc,natom, &
&                      nprimij,nangij,nbfij,len1,mxprsh,threshex)
!
        endif
      else
        cint2(1:nbfij(2),1:nbfij(1))= zero
      endif
!
      select case(mtype(jsh))
        case (0)
          do i= 1,nbfij(1)
            dcint(1,i,1)= cint1(1,i)
            dcint(1,i,2)= cint1(2,i)
            dcint(1,i,3)= cint1(3,i)
          enddo
        case (1)
          do i= 1,nbfij(1)
            dcint(1,i,1)= cint1(1,i)          -cint2(1,i)
            dcint(2,i,1)= cint1(2,i)*sqrtthird
            dcint(3,i,1)= cint1(3,i)*sqrtthird
            dcint(1,i,2)= cint1(2,i)*sqrtthird
            dcint(2,i,2)= cint1(4,i)          -cint2(1,i)
            dcint(3,i,2)= cint1(5,i)*sqrtthird
            dcint(1,i,3)= cint1(3,i)*sqrtthird
            dcint(2,i,3)= cint1(5,i)*sqrtthird
            dcint(3,i,3)= cint1(6,i)          -cint2(1,i)
          enddo
        case (2)
          do i= 1,nbfij(1)
            dcint(1,i,1)= cint1( 1,i)           -cint2(1,i)*two
            dcint(2,i,1)= cint1( 2,i)*sqrt3fifth-cint2(2,i)*sqrt3
            dcint(3,i,1)= cint1( 3,i)*sqrt3fifth-cint2(3,i)*sqrt3
            dcint(4,i,1)= cint1( 4,i)*sqrtfifth
            dcint(5,i,1)= cint1( 5,i)*sqrtfifth
            dcint(6,i,1)= cint1( 6,i)*sqrtfifth
            dcint(1,i,2)= cint1( 2,i)*sqrtfifth
            dcint(2,i,2)= cint1( 4,i)*sqrt3fifth-cint2(1,i)*sqrt3
            dcint(3,i,2)= cint1( 5,i)*sqrtfifth
            dcint(4,i,2)= cint1( 7,i)           -cint2(2,i)*two
            dcint(5,i,2)= cint1( 8,i)*sqrt3fifth-cint2(3,i)*sqrt3
            dcint(6,i,2)= cint1( 9,i)*sqrtfifth
            dcint(1,i,3)= cint1( 3,i)*sqrtfifth
            dcint(2,i,3)= cint1( 5,i)*sqrtfifth
            dcint(3,i,3)= cint1( 6,i)*sqrt3fifth-cint2(1,i)*sqrt3
            dcint(4,i,3)= cint1( 8,i)*sqrtfifth
            dcint(5,i,3)= cint1( 9,i)*sqrt3fifth-cint2(2,i)*sqrt3
            dcint(6,i,3)= cint1(10,i)           -cint2(3,i)*two
          enddo
        case (3)
          do i= 1,nbfij(1)
            dcint( 1,i,1)= cint1( 1,i)            -cint2(1,i)*three
            dcint( 2,i,1)= cint1( 2,i)*sqrtseventh-cint2(2,i)*two*sqrtthird
            dcint( 3,i,1)= cint1( 3,i)*sqrtseventh-cint2(3,i)*two*sqrtthird
            dcint( 4,i,1)= cint1( 4,i)*sqrt3inv35 -cint2(4,i)
            dcint( 5,i,1)= cint1( 5,i)*sqrtinv35  -cint2(5,i)*sqrtthird
            dcint( 6,i,1)= cint1( 6,i)*sqrt3inv35 -cint2(6,i)
            dcint( 7,i,1)= cint1( 7,i)*sqrtseventh
            dcint( 8,i,1)= cint1( 8,i)*sqrtinv35
            dcint( 9,i,1)= cint1( 9,i)*sqrtinv35
            dcint(10,i,1)= cint1(10,i)*sqrtseventh
            dcint( 1,i,2)= cint1( 2,i)*sqrtseventh
            dcint( 2,i,2)= cint1( 4,i)*sqrt3inv35 -cint2(1,i)
            dcint( 3,i,2)= cint1( 5,i)*sqrtinv35
            dcint( 4,i,2)= cint1( 7,i)*sqrtseventh-cint2(2,i)*two*sqrtthird
            dcint( 5,i,2)= cint1( 8,i)*sqrtinv35  -cint2(3,i)*sqrtthird
            dcint( 6,i,2)= cint1( 9,i)*sqrtinv35
            dcint( 7,i,2)= cint1(11,i)            -cint2(4,i)*three
            dcint( 8,i,2)= cint1(12,i)*sqrtseventh-cint2(5,i)*two*sqrtthird
            dcint( 9,i,2)= cint1(13,i)*sqrt3inv35 -cint2(6,i)
            dcint(10,i,2)= cint1(14,i)*sqrtseventh
            dcint( 1,i,3)= cint1( 3,i)*sqrtseventh
            dcint( 2,i,3)= cint1( 5,i)*sqrtinv35
            dcint( 3,i,3)= cint1( 6,i)*sqrt3inv35 -cint2(1,i)
            dcint( 4,i,3)= cint1( 8,i)*sqrtinv35
            dcint( 5,i,3)= cint1( 9,i)*sqrtinv35  -cint2(2,i)*sqrtthird
            dcint( 6,i,3)= cint1(10,i)*sqrtseventh-cint2(3,i)*two*sqrtthird
            dcint( 7,i,3)= cint1(12,i)*sqrtseventh
            dcint( 8,i,3)= cint1(13,i)*sqrt3inv35 -cint2(4,i)
            dcint( 9,i,3)= cint1(14,i)*sqrtseventh-cint2(5,i)*two*sqrtthird
            dcint(10,i,3)= cint1(15,i)            -cint2(6,i)*three
          enddo
        case (4)
          do i= 1,nbfij(1)
            dcint( 1,i,1)= cint1( 1,i)           -cint2( 1,i)*four
            dcint( 2,i,1)= cint1( 2,i)*third     -cint2( 2,i)*three*sqrtfifth
            dcint( 3,i,1)= cint1( 3,i)*third     -cint2( 3,i)*three*sqrtfifth
            dcint( 4,i,1)= cint1( 4,i)*sqrtinv21 -cint2( 4,i)*two*sqrtfifth
            dcint( 5,i,1)= cint1( 5,i)*sqrtinv63 -cint2( 5,i)*two*sqrtinv15
            dcint( 6,i,1)= cint1( 6,i)*sqrtinv21 -cint2( 6,i)*two*sqrtfifth
            dcint( 7,i,1)= cint1( 7,i)*sqrtinv21 -cint2( 7,i)
            dcint( 8,i,1)= cint1( 8,i)*sqrtinv105-cint2( 8,i)*sqrtfifth
            dcint( 9,i,1)= cint1( 9,i)*sqrtinv105-cint2( 9,i)*sqrtfifth
            dcint(10,i,1)= cint1(10,i)*sqrtinv21 -cint2(10,i)
            dcint(11,i,1)= cint1(11,i)*third
            dcint(12,i,1)= cint1(12,i)*sqrtinv63
            dcint(13,i,1)= cint1(13,i)*sqrtinv105
            dcint(14,i,1)= cint1(14,i)*sqrtinv63
            dcint(15,i,1)= cint1(15,i)*third
            dcint( 1,i,2)= cint1( 2,i)*third
            dcint( 2,i,2)= cint1( 4,i)*sqrtinv21 -cint2( 1,i)
            dcint( 3,i,2)= cint1( 5,i)*sqrtinv63
            dcint( 4,i,2)= cint1( 7,i)*sqrtinv21 -cint2( 2,i)*two*sqrtfifth
            dcint( 5,i,2)= cint1( 8,i)*sqrtinv105-cint2( 3,i)*sqrtfifth
            dcint( 6,i,2)= cint1( 9,i)*sqrtinv105
            dcint( 7,i,2)= cint1(11,i)*third     -cint2( 4,i)*three*sqrtfifth
            dcint( 8,i,2)= cint1(12,i)*sqrtinv63 -cint2( 5,i)*two*sqrtinv15
            dcint( 9,i,2)= cint1(13,i)*sqrtinv105-cint2( 6,i)*sqrtfifth
            dcint(10,i,2)= cint1(14,i)*sqrtinv63
            dcint(11,i,2)= cint1(16,i)           -cint2( 7,i)*four
            dcint(12,i,2)= cint1(17,i)*third     -cint2( 8,i)*three*sqrtfifth
            dcint(13,i,2)= cint1(18,i)*sqrtinv21 -cint2( 9,i)*two*sqrtfifth
            dcint(14,i,2)= cint1(19,i)*sqrtinv21 -cint2(10,i)
            dcint(15,i,2)= cint1(20,i)*third
            dcint( 1,i,3)= cint1( 3,i)*third
            dcint( 2,i,3)= cint1( 5,i)*sqrtinv63
            dcint( 3,i,3)= cint1( 6,i)*sqrtinv21 -cint2( 1,i)
            dcint( 4,i,3)= cint1( 8,i)*sqrtinv105
            dcint( 5,i,3)= cint1( 9,i)*sqrtinv105-cint2( 2,i)*sqrtfifth
            dcint( 6,i,3)= cint1(10,i)*sqrtinv21 -cint2( 3,i)*two*sqrtfifth
            dcint( 7,i,3)= cint1(12,i)*sqrtinv63
            dcint( 8,i,3)= cint1(13,i)*sqrtinv105-cint2( 4,i)*sqrtfifth
            dcint( 9,i,3)= cint1(14,i)*sqrtinv63 -cint2( 5,i)*two*sqrtinv15
            dcint(10,i,3)= cint1(15,i)*third     -cint2( 6,i)*three*sqrtfifth
            dcint(11,i,3)= cint1(17,i)*third
            dcint(12,i,3)= cint1(18,i)*sqrtinv21 -cint2( 7,i)
            dcint(13,i,3)= cint1(19,i)*sqrtinv21 -cint2( 8,i)*two*sqrtfifth
            dcint(14,i,3)= cint1(20,i)*third     -cint2( 9,i)*three*sqrtfifth
            dcint(15,i,3)= cint1(21,i)           -cint2(10,i)*four
          enddo
        case (5)
          do i= 1,nbfij(1)
            dcint( 1,i,1)= cint1( 1,i)            -cint2( 1,i)*five
            dcint( 2,i,1)= cint1( 2,i)*sqrtinv11  -cint2( 2,i)*four*sqrtseventh
            dcint( 3,i,1)= cint1( 3,i)*sqrtinv11  -cint2( 3,i)*four*sqrtseventh
            dcint( 4,i,1)= cint1( 4,i)*sqrtinv33  -cint2( 4,i)*three*sqrt3inv35
            dcint( 5,i,1)= cint1( 5,i)*sqrtinv99  -cint2( 5,i)*three*sqrtinv35
            dcint( 6,i,1)= cint1( 6,i)*sqrtinv33  -cint2( 6,i)*three*sqrt3inv35
            dcint( 7,i,1)= cint1( 7,i)*sqrt5inv231-cint2( 7,i)*two*sqrtseventh
            dcint( 8,i,1)= cint1( 8,i)*sqrtinv231 -cint2( 8,i)*two*sqrtinv35
            dcint( 9,i,1)= cint1( 9,i)*sqrtinv231 -cint2( 9,i)*two*sqrtinv35
            dcint(10,i,1)= cint1(10,i)*sqrt5inv231-cint2(10,i)*two*sqrtseventh
            dcint(11,i,1)= cint1(11,i)*sqrtinv33  -cint2(11,i)
            dcint(12,i,1)= cint1(12,i)*sqrtinv231 -cint2(12,i)*sqrtseventh
            dcint(13,i,1)= cint1(13,i)*sqrtinv385 -cint2(13,i)*sqrt3inv35
            dcint(14,i,1)= cint1(14,i)*sqrtinv231 -cint2(14,i)*sqrtseventh
            dcint(15,i,1)= cint1(15,i)*sqrtinv33  -cint2(15,i)
            dcint(16,i,1)= cint1(16,i)*sqrtinv11
            dcint(17,i,1)= cint1(17,i)*sqrtinv99
            dcint(18,i,1)= cint1(18,i)*sqrtinv231
            dcint(19,i,1)= cint1(19,i)*sqrtinv231
            dcint(20,i,1)= cint1(20,i)*sqrtinv99
            dcint(21,i,1)= cint1(21,i)*sqrtinv11
            dcint( 1,i,2)= cint1( 2,i)*sqrtinv11
            dcint( 2,i,2)= cint1( 4,i)*sqrtinv33  -cint2( 1,i)
            dcint( 3,i,2)= cint1( 5,i)*sqrtinv99
            dcint( 4,i,2)= cint1( 7,i)*sqrt5inv231-cint2( 2,i)*two*sqrtseventh
            dcint( 5,i,2)= cint1( 8,i)*sqrtinv231 -cint2( 3,i)*sqrtseventh
            dcint( 6,i,2)= cint1( 9,i)*sqrtinv231
            dcint( 7,i,2)= cint1(11,i)*sqrtinv33  -cint2( 4,i)*three*sqrt3inv35
            dcint( 8,i,2)= cint1(12,i)*sqrtinv231 -cint2( 5,i)*two*sqrtinv35
            dcint( 9,i,2)= cint1(13,i)*sqrtinv385 -cint2( 6,i)*sqrt3inv35
            dcint(10,i,2)= cint1(14,i)*sqrtinv231
            dcint(11,i,2)= cint1(16,i)*sqrtinv11  -cint2( 7,i)*four*sqrtseventh
            dcint(12,i,2)= cint1(17,i)*sqrtinv99  -cint2( 8,i)*three*sqrtinv35
            dcint(13,i,2)= cint1(18,i)*sqrtinv231 -cint2( 9,i)*two*sqrtinv35
            dcint(14,i,2)= cint1(19,i)*sqrtinv231 -cint2(10,i)*sqrtseventh
            dcint(15,i,2)= cint1(20,i)*sqrtinv99
            dcint(16,i,2)= cint1(22,i)            -cint2(11,i)*five
            dcint(17,i,2)= cint1(23,i)*sqrtinv11  -cint2(12,i)*four*sqrtseventh
            dcint(18,i,2)= cint1(24,i)*sqrtinv33  -cint2(13,i)*three*sqrt3inv35
            dcint(19,i,2)= cint1(25,i)*sqrt5inv231-cint2(14,i)*two*sqrtseventh
            dcint(20,i,2)= cint1(26,i)*sqrtinv33  -cint2(15,i)
            dcint(21,i,2)= cint1(27,i)*sqrtinv11
            dcint( 1,i,3)= cint1( 3,i)*sqrtinv11
            dcint( 2,i,3)= cint1( 5,i)*sqrtinv99
            dcint( 3,i,3)= cint1( 6,i)*sqrtinv33  -cint2( 1,i)
            dcint( 4,i,3)= cint1( 8,i)*sqrtinv231
            dcint( 5,i,3)= cint1( 9,i)*sqrtinv231 -cint2( 2,i)*sqrtseventh
            dcint( 6,i,3)= cint1(10,i)*sqrt5inv231-cint2( 3,i)*two*sqrtseventh
            dcint( 7,i,3)= cint1(12,i)*sqrtinv231
            dcint( 8,i,3)= cint1(13,i)*sqrtinv385 -cint2( 4,i)*sqrt3inv35
            dcint( 9,i,3)= cint1(14,i)*sqrtinv231 -cint2( 5,i)*two*sqrtinv35
            dcint(10,i,3)= cint1(15,i)*sqrtinv33  -cint2( 6,i)*three*sqrt3inv35
            dcint(11,i,3)= cint1(17,i)*sqrtinv99
            dcint(12,i,3)= cint1(18,i)*sqrtinv231 -cint2( 7,i)*sqrtseventh
            dcint(13,i,3)= cint1(19,i)*sqrtinv231 -cint2( 8,i)*two*sqrtinv35
            dcint(14,i,3)= cint1(20,i)*sqrtinv99  -cint2( 9,i)*three*sqrtinv35
            dcint(15,i,3)= cint1(21,i)*sqrtinv11  -cint2(10,i)*four*sqrtseventh
            dcint(16,i,3)= cint1(23,i)*sqrtinv11
            dcint(17,i,3)= cint1(24,i)*sqrtinv33  -cint2(11,i)
            dcint(18,i,3)= cint1(25,i)*sqrt5inv231-cint2(12,i)*two*sqrtseventh
            dcint(19,i,3)= cint1(26,i)*sqrtinv33  -cint2(13,i)*three*sqrt3inv35
            dcint(20,i,3)= cint1(27,i)*sqrtinv11  -cint2(14,i)*four*sqrtseventh
            dcint(21,i,3)= cint1(28,i)            -cint2(15,i)*five
          enddo
      end select
!
      nbfij(2)  = mbf(jsh)
      select case(nbfij(2))
        case(5)
          do k= 1,3
            do i= 1,nbfij(1)
              do j= 1,6
                work(j)= dcint(j,i,k)
              enddo
              dcint(1,i,k)= work(2)
              dcint(2,i,k)= work(5)
              dcint(3,i,k)=(work(6)*two-work(1)-work(4))*half
              dcint(4,i,k)= work(3)
              dcint(5,i,k)=(work(1)-work(4))*sqrt3h
            enddo
          enddo
        case(7)
          do k= 1,3
            do i= 1,nbfij(1)
              do j= 1,10
                work(j)= dcint(j,i,k)
              enddo
              dcint(1,i,k)=(-work(7)+work(2)*three                   )*facf1
              dcint(2,i,k)=  work(5)                                  *facf2
              dcint(3,i,k)=(-work(7)-work(2)+work(9)*four            )*facf3
              dcint(4,i,k)=( work(10)*two-work(3)*three-work(8)*three)*half
              dcint(5,i,k)=(-work(1)-work(4)+work(6)*four            )*facf3
              dcint(6,i,k)=( work(3)-work(8)                         )*facf4
              dcint(7,i,k)=( work(1)-work(4)*three                   )*facf1
            enddo
          enddo
        case(10)
          do k= 1,3
            do i= 1,nbfij(1)
              dcint(2,i,k)= dcint(2,i,k)*sqrt5
              dcint(3,i,k)= dcint(3,i,k)*sqrt5
              dcint(4,i,k)= dcint(4,i,k)*sqrt5
              dcint(5,i,k)= dcint(5,i,k)*sqrt15
              dcint(6,i,k)= dcint(6,i,k)*sqrt5
              dcint(8,i,k)= dcint(8,i,k)*sqrt5
              dcint(9,i,k)= dcint(9,i,k)*sqrt5
            enddo
          enddo
        case(9)
          do k= 1,3
            do i= 1,nbfij(1)
              do j= 1,15
                work(j)= dcint(j,i,k)
              enddo
              dcint(1,i,k)=(work(2)-work(7))*facg1
              dcint(2,i,k)=(-work(12)+work(5)*three)*facg2
              dcint(3,i,k)=(-work(2)-work(7)+work(9)*six)*facg3
              dcint(4,i,k)=(-work(12)*three+work(14)*four-work(5)*three)*facg4
              dcint(5,i,k)=(work(1)*three+work(11)*three+work(15)*eight+work(4)*six &
&                          -work(6)*p24-work(13)*p24)*eighth
              dcint(6,i,k)=(-work(3)*three+work(10)*four-work(8)*three)*facg4
              dcint(7,i,k)=(-work(1)+work(11)+work(6)*six-work(13)*six)*facg5
              dcint(8,i,k)=(work(3)-work(8)*three)*facg2
              dcint(9,i,k)=(work(1)+work(11)-work(4)*six)*facg6
            enddo
          enddo
        case(15)
          do k= 1,3
            do i= 1,nbfij(1)
              dcint( 2,i,k)= dcint( 2,i,k)*sqrt7
              dcint( 3,i,k)= dcint( 3,i,k)*sqrt7
              dcint( 4,i,k)= dcint( 4,i,k)*sqrt35third
              dcint( 5,i,k)= dcint( 5,i,k)*sqrt35
              dcint( 6,i,k)= dcint( 6,i,k)*sqrt35third
              dcint( 7,i,k)= dcint( 7,i,k)*sqrt7
              dcint( 8,i,k)= dcint( 8,i,k)*sqrt35
              dcint( 9,i,k)= dcint( 9,i,k)*sqrt35
              dcint(10,i,k)= dcint(10,i,k)*sqrt7
              dcint(12,i,k)= dcint(12,i,k)*sqrt7
              dcint(13,i,k)= dcint(13,i,k)*sqrt35third
              dcint(14,i,k)= dcint(14,i,k)*sqrt7
            enddo
          enddo
        case(11)
          do k= 1,3
            do i= 1,nbfij(1)
              do j= 1,21
                work(j)= dcint(j,i,k)
              enddo
              dcint( 1,i,k)=(work(2)*five-work(7)*ten+work(16))*fach1
              dcint( 2,i,k)=(work(5)*four-work(12)*four)*fach2
              dcint( 3,i,k)=(-work(2)*three-work(7)*two+work(9)*p24+work(16) &
&                           -work(18)*eight)*fach3
              dcint( 4,i,k)=(-work(5)*two-work(12)*two+work(14)*four)*fach4
              dcint( 5,i,k)=(work(2)+work(7)*two-work(9)*twelve+work(16)-work(18)*twelve &
&                           +work(20)*eight)*fach5
              dcint( 6,i,k)=(work(3)*p15+work(8)*p30-work(10)*p40+work(17)*p15 &
&                           -work(19)*p40+work(21)*eight)*eighth
              dcint( 7,i,k)=(work(1)+work(4)*two-work(6)*twelve+work(11)-work(13)*twelve &
&                           +work(15)*eight)*fach5
              dcint( 8,i,k)=(-work(3)+work(10)*two+work(17)-work(19)*two)*fach4
              dcint( 9,i,k)=(-work(1)+work(4)*two+work(6)*eight+work(11)*three &
&                           -work(13)*p24)*fach3
              dcint(10,i,k)=(work(3)-work(8)*six+work(17))*fach2
              dcint(11,i,k)=(work(1)-work(4)*ten+work(11)*five)*fach1
            enddo
          enddo
        case(21)
          do k= 1,3
            do i= 1,nbfij(1)
              dcint( 2,i,k)= dcint( 2,i,k)*three
              dcint( 3,i,k)= dcint( 3,i,k)*three
              dcint( 4,i,k)= dcint( 4,i,k)*sqrt21
              dcint( 5,i,k)= dcint( 5,i,k)*sqrt63
              dcint( 6,i,k)= dcint( 6,i,k)*sqrt21
              dcint( 7,i,k)= dcint( 7,i,k)*sqrt21
              dcint( 8,i,k)= dcint( 8,i,k)*sqrt105
              dcint( 9,i,k)= dcint( 9,i,k)*sqrt105
              dcint(10,i,k)= dcint(10,i,k)*sqrt21
              dcint(11,i,k)= dcint(11,i,k)*three
              dcint(12,i,k)= dcint(12,i,k)*sqrt63
              dcint(13,i,k)= dcint(13,i,k)*sqrt105
              dcint(14,i,k)= dcint(14,i,k)*sqrt63
              dcint(15,i,k)= dcint(15,i,k)*three
              dcint(17,i,k)= dcint(17,i,k)*three
              dcint(18,i,k)= dcint(18,i,k)*sqrt21
              dcint(19,i,k)= dcint(19,i,k)*sqrt21
              dcint(20,i,k)= dcint(20,i,k)*three
            enddo
          enddo
      end select
!
      do i= 1,nbfij(1)
        ii= ilocbf+i
        do j= 1,nbfij(2)
          egrad(1,jatom)= egrad(1,jatom)+fulldmtrx(jlocbf+j,ii)*dcint(j,i,1)
          egrad(2,jatom)= egrad(2,jatom)+fulldmtrx(jlocbf+j,ii)*dcint(j,i,2)
          egrad(3,jatom)= egrad(3,jatom)+fulldmtrx(jlocbf+j,ii)*dcint(j,i,3)
        enddo
      enddo
      return
end


!-------------------------------------------------
  subroutine calchelfey(egrad,fulldmtrx,ish,jsh)
!-------------------------------------------------
!
! Driver of Helmann-Feynman gradient term
!
! In : fulldmtrx (density matrix)
!    : ish, jsh (shell indices)
! Inout : egrad (energy gradient value)
!
      use modparam, only : mxprsh
      use modmolecule, only : natom, coord, znuc
      use modbasis, only : locatom, locprim, locbf, mprim, mbf, mtype, ex, coeff, nao
      use modthresh, only : threshex
      implicit none
      integer,intent(in) :: ish, jsh
      integer :: nangij(2), nprimij(2), nbfij(2), locbfij(2), iatom, jatom
      integer :: iloc, jloc, iprim, jprim, i
      real(8),intent(in) :: fulldmtrx(nao,nao)
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: exij(mxprsh,2), cij(mxprsh,2), coordij(3,2)
      logical :: iandj
!
      iandj=(ish == jsh)
      nangij(1)= mtype(ish)
      nangij(2)= mtype(jsh)
      nprimij(1)= mprim(ish)
      nprimij(2)= mprim(jsh)
      nbfij(1)  = mbf(ish)
      nbfij(2)  = mbf(jsh)
      locbfij(1)= locbf(ish)
      locbfij(2)= locbf(jsh)
      iatom = locatom(ish)
      iloc  = locprim(ish)
      jatom = locatom(jsh)
      jloc  = locprim(jsh)
      do i= 1,3
        coordij(i,1)= coord(i,iatom)
        coordij(i,2)= coord(i,jatom)
      enddo
      do iprim= 1,nprimij(1)
        exij(iprim,1)= ex(iloc+iprim)
        cij(iprim,1) = coeff(iloc+iprim)
      enddo
      do jprim= 1,nprimij(2)
        exij(jprim,2)= ex(jloc+jprim)
        cij(jprim,2) = coeff(jloc+jprim)
      enddo
!
      if((nangij(1) <= 2).and.(nangij(2) <= 2)) then
        call int1cgmd(egrad,fulldmtrx,exij,cij,coordij,coord,znuc,natom,nao, &
&                     nprimij,nangij,nbfij,locbfij,mxprsh,threshex,iandj)
      elseif((nangij(1) <= 5).and.(nangij(2) <= 5)) then
        call int1grys(egrad,fulldmtrx,exij,cij,coordij,coord,znuc,natom,nao, &
&                     nprimij,nangij,nbfij,locbfij,mxprsh,threshex,iandj)
      else
        write(*,'(" Error! This program supports up to h function in helfey")')
        call iabort
      endif
!
      return
end


!-------------------------------------------------------------------------------
  subroutine int1cgmd(egrad,fulldmtrx,exij,cij,coordij,coord,znuc,natom,nao, &
&                     nprimij,nangij,nbfij,locbfij,mxprsh,threshex,iandj)
!-------------------------------------------------------------------------------
!
! Driver of first derivative of 1-electron Coulomb integrals (j|(Z/r)'|i) 
! using McMurchie-Davidson method
!
! In : exij     (exponents of basis functions)
!      coij     (coefficients of basis functions)
!      coordij  (x,y,z coordinates of basis functions)
!      coord    (x,y,z coordinates of atoms)
!      znuc     (charges of atoms)
!      natom    (number of atoms)
!      nprimij  (numbers of primitive functions)
!      nangij   (degrees of angular momentum)
!      nbfij    (numbers of basis functions)
!      locbfij  (starting addresses of basis functions)
!      mxprsh   (size of primitive fuction array)
!      threshex (threshold of exponential calculation)
!      iandj    (flag of ish == jsh)
! Out: egrad    (energy gradient value)
!
      implicit none
      integer,intent(in) :: nprimij(2), nangij(2), nbfij(2), locbfij(2), natom, nao, mxprsh
      integer,parameter :: mxprsh2=30
      integer :: inttyp, nij, iprim, jprim, i, ii, jj, nbfij2(2), locbfij2(2)
      real(8),parameter :: one=1.0D+00, pi2=6.283185307179586D+00
      real(8),intent(in) :: fulldmtrx(nao,nao), exij(mxprsh,2), cij(mxprsh,2), coordij(3,2)
      real(8),intent(in) :: coord(3,natom), znuc(natom), threshex
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: xyz(3), rij, exi, exj, ci, cj, ex12, ex21, ex2i, ex2j, rij2, pixyz(3)
      real(8) :: exfac(5,mxprsh2*mxprsh2), pijxyz(3,mxprsh2*mxprsh2)
      logical,intent(in) :: iandj
!
      if(mxprsh > mxprsh2) then
        write(*,'(" Error! Parameter mxprsh2 in int1cgmd is small!")')
        call exit
      endif
!
      inttyp=nangij(2)*3+nangij(1)+1
      if(nangij(2).ge.nangij(1)) then
        ii=1
        jj=2
        nbfij2(1)= nbfij(1)
        nbfij2(2)= nbfij(2)
        locbfij2(1)= locbfij(1)
        locbfij2(2)= locbfij(2)
      else
        ii=2
        jj=1
        nbfij2(1)= nbfij(2)
        nbfij2(2)= nbfij(1)
        locbfij2(1)= locbfij(2)
        locbfij2(2)= locbfij(1)
      endif
!
      do i= 1,3
        xyz(i)= coordij(i,ii)-coordij(i,jj)
      enddo
      rij= xyz(1)*xyz(1)+xyz(2)*xyz(2)+xyz(3)*xyz(3)
      nij= 0
      do iprim= 1,nprimij(ii)
        exi = exij(iprim,ii)
        ci  = cij(iprim,ii)*pi2
        do i= 1,3
          pixyz(i)= exi*coordij(i,ii)
        enddo
        do jprim= 1,nprimij(jj)
          exj = exij(jprim,jj)
          ex12= exi+exj
          ex21= one/ex12
          ex2i= ex21*exi
          ex2j= ex21*exj
          rij2= rij*ex2i*exj
          if(rij2 > threshex) cycle
          nij= nij+1
          cj = cij(jprim,jj)
          exfac(1,nij)= ex12
          exfac(2,nij)= ex21
          exfac(3,nij)= ex2i
          exfac(4,nij)= ex2j
          exfac(5,nij)= exp(-rij2)*ex21*ci*cj
          do i= 1,3
            pijxyz(i,nij)=(pixyz(i)+exj*coordij(i,jj))*ex21
          enddo
        enddo
      enddo
!
      select case(inttyp)
        case (1)
          call int1gcss(egrad,fulldmtrx,exfac,pijxyz,nij,coord,znuc, &
&                       natom,nao,mxprsh,locbfij2,iandj)
        case (2,4)
          call int1gcps(egrad,fulldmtrx,exfac,pijxyz,xyz,nij,coord,znuc, &
&                       natom,nao,mxprsh,locbfij2)
        case (5)
          call int1gcpp(egrad,fulldmtrx,exfac,pijxyz,xyz,nij,coord,znuc, &
&                       natom,nao,mxprsh,locbfij2,iandj)
        case (3,7)
          call int1gcds(egrad,fulldmtrx,exfac,pijxyz,xyz,nij,coord,znuc, &
&                       natom,nao,mxprsh,locbfij2,nbfij2)
        case (6,8)
          call int1gcdp(egrad,fulldmtrx,exfac,pijxyz,xyz,nij,coord,znuc, &
&                       natom,nao,mxprsh,locbfij2,nbfij2)
        case (9)
          call int1gcdd(egrad,fulldmtrx,exfac,pijxyz,xyz,nij,coord,znuc, &
&                       natom,nao,mxprsh,locbfij2,nbfij2,iandj)
      end select
!
      return
end


!-----------------------------------------------------------------------
  subroutine int1gcss(egrad,fulldmtrx,exfac,pijxyz,nij,coord,znuc, &
&                     natom,nao,mxprsh,locbfij,iandj)
!-----------------------------------------------------------------------
!
! Calculate Helmann-Feynman gradient term <s|V'|s>
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nij, natom, nao, mxprsh, locbfij(2)
      integer :: ijprim, i, k, iatom, igrid
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, two=2.0D+00
      real(8),parameter :: pi=3.141592653589793D+00
      real(8),intent(in) :: fulldmtrx(nao,nao), exfac(5,mxprsh*mxprsh)
      real(8),intent(in) :: pijxyz(3,mxprsh*mxprsh), coord(3,natom), znuc(natom)
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: rc, ex12, c12, pcxyz(3), tinv
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: ft(0:1), r1(3)
      logical,intent(in) :: iandj
!
      do iatom= 1,natom
        do i= 1,3
          r1(i)= zero
        enddo
        do ijprim= 1,nij
          ex12= exfac(1,ijprim)
          c12 = exfac(5,ijprim)
          do i= 1,3
            pcxyz(i)= pijxyz(i,ijprim)-coord(i,iatom)
          enddo
          rc= pcxyz(1)*pcxyz(1)+pcxyz(2)*pcxyz(2)+pcxyz(3)*pcxyz(3)
          tval= ex12*rc
          if(tval >= threshtval) then
            tinv = one/tval
            ft(0)= half*sqrt(pi*tinv)
            ft(1)= ft(0)*half*tinv
          else
            igrid= int(tval)
            tval2= tval *tval
            tval3= tval2*tval
            tval4= tval2*tval2
            tval5= tval2*tval3
            tval6= tval3*tval3
            tval7= tval4*tval3
            tval8= tval4*tval4
            tval9= tval4*tval5
            tval10=tval5*tval5
            ft(0)= fgrid(0,0,igrid)      +fgrid( 1,0,igrid)*tval  +fgrid( 2,0,igrid)*tval2 &
&                 +fgrid(3,0,igrid)*tval3+fgrid( 4,0,igrid)*tval4 +fgrid( 5,0,igrid)*tval5 &
&                 +fgrid(6,0,igrid)*tval6+fgrid( 7,0,igrid)*tval7 +fgrid( 8,0,igrid)*tval8 &
&                 +fgrid(9,0,igrid)*tval9+fgrid(10,0,igrid)*tval10
            ft(1)= fgrid(0,1,igrid)      +fgrid( 1,1,igrid)*tval  +fgrid( 2,1,igrid)*tval2 &
&                 +fgrid(3,1,igrid)*tval3+fgrid( 4,1,igrid)*tval4 +fgrid( 5,1,igrid)*tval5 &
&                 +fgrid(6,1,igrid)*tval6+fgrid( 7,1,igrid)*tval7 +fgrid( 8,1,igrid)*tval8 &
&                 +fgrid(9,1,igrid)*tval9+fgrid(10,1,igrid)*tval10
          endif
          ft(1)=-ft(1)*two*ex12*c12*znuc(iatom)
          do i= 1,3
            r1(i)= r1(i)+ft(1)*pcxyz(i)
          enddo
        enddo
        if(.not.iandj) then
          do k= 1,3
            egrad(k,iatom)= egrad(k,iatom)+fulldmtrx(locbfij(2)+1,locbfij(1)+1)*r1(k)
          enddo
        else
          do k= 1,3
            egrad(k,iatom)= egrad(k,iatom)+fulldmtrx(locbfij(2)+1,locbfij(1)+1)*r1(k)*half
          enddo
        endif
      enddo
!
      return
end


!---------------------------------------------------------------------------
  subroutine int1gcps(egrad,fulldmtrx,exfac,pijxyz,xyz,nij,coord,znuc, &
&                     natom,nao,mxprsh,locbfij)
!---------------------------------------------------------------------------
!
! Calculate Helmann-Feynman gradient term <p|V'|s>
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nij, natom, nao, mxprsh, locbfij(2)
      integer :: ijprim, i, j, k, iatom, igrid, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, two=2.0D+00
      real(8),parameter :: three=3.0D+00, p15=1.5D+00
      real(8),parameter :: pi=3.141592653589793D+00
      real(8),intent(in) :: fulldmtrx(nao,nao), exfac(5,mxprsh*mxprsh), pijxyz(3,mxprsh*mxprsh)
      real(8),intent(in) :: xyz(3), coord(3,natom), znuc(natom)
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: rc, ex12, ex21, ex2i, ex2j, c12, pcxyz(3), extwo
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: ft(0:2), tinv, r1(3), r2(6), cint1(3,3)
!
      do iatom= 1,natom
        do i= 1,3
          r1(i)= zero
        enddo
        do i= 1,6
          r2(i)= zero
        enddo
        do ijprim= 1,nij
          ex12= exfac(1,ijprim)
          ex21= exfac(2,ijprim)
          ex2i= exfac(3,ijprim)
          ex2j= exfac(4,ijprim)
          c12 = exfac(5,ijprim)
          do i= 1,3
            pcxyz(i)= pijxyz(i,ijprim)-coord(i,iatom)
          enddo
          rc= pcxyz(1)*pcxyz(1)+pcxyz(2)*pcxyz(2)+pcxyz(3)*pcxyz(3)
          tval= ex12*rc
          if(tval >= threshtval) then
            tinv = one/tval
            ft(0)= half*sqrt(pi*tinv)
            ft(1)= half*tinv*ft(0)
            ft(2)= p15 *tinv*ft(1)
          else
            igrid= int(tval)
            tval2= tval *tval
            tval3= tval2*tval
            tval4= tval2*tval2
            tval5= tval2*tval3
            tval6= tval3*tval3
            tval7= tval4*tval3
            tval8= tval4*tval4
            tval9= tval4*tval5
            tval10=tval5*tval5
            do ii= 0,2
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
          endif
          do i= 1,2
            ft(i)= ft(i)*znuc(iatom)*c12
          enddo
          extwo= ex12*two
          do i= 1,3
            r1(i)= r1(i)-ft(1)*pcxyz(i)*ex2i*extwo
          enddo
          r2(1)= r2(1)+ft(2)*pcxyz(1)*pcxyz(1)*extwo-ft(1)
          r2(2)= r2(2)+ft(2)*pcxyz(2)*pcxyz(2)*extwo-ft(1)
          r2(3)= r2(3)+ft(2)*pcxyz(3)*pcxyz(3)*extwo-ft(1)
          r2(4)= r2(4)+ft(2)*pcxyz(1)*pcxyz(2)*extwo
          r2(5)= r2(5)+ft(2)*pcxyz(1)*pcxyz(3)*extwo
          r2(6)= r2(6)+ft(2)*pcxyz(2)*pcxyz(3)*extwo
        enddo
        cint1(1,1)= r2(1)+r1(1)*xyz(1)
        cint1(2,1)= r2(4)+r1(1)*xyz(2)
        cint1(3,1)= r2(5)+r1(1)*xyz(3)
        cint1(1,2)= r2(4)+r1(2)*xyz(1)
        cint1(2,2)= r2(2)+r1(2)*xyz(2)
        cint1(3,2)= r2(6)+r1(2)*xyz(3)
        cint1(1,3)= r2(5)+r1(3)*xyz(1)
        cint1(2,3)= r2(6)+r1(3)*xyz(2)
        cint1(3,3)= r2(3)+r1(3)*xyz(3)
        do k= 1,3
          do j= 1,3
            egrad(k,iatom)= egrad(k,iatom)+fulldmtrx(locbfij(2)+j,locbfij(1)+1)*cint1(j,k)
          enddo
        enddo
      enddo
!
      return
end


!---------------------------------------------------------------------------
  subroutine int1gcpp(egrad,fulldmtrx,exfac,pijxyz,xyz,nij,coord,znuc, &
&                     natom,nao,mxprsh,locbfij,iandj)
!---------------------------------------------------------------------------
!
! Calculate Helmann-Feynman gradient term <p|V'|p>
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nij, natom, nao, mxprsh, locbfij(2)
      integer :: ijprim, i, j, k, iatom, igrid, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, two=2.0D+00
      real(8),parameter :: three=3.0D+00, five=5.0D+00, p15=1.5D+00, p25=2.5D+00
      real(8),parameter :: pi=3.141592653589793D+00
      real(8),intent(in) :: fulldmtrx(nao,nao), exfac(5,mxprsh*mxprsh), pijxyz(3,mxprsh*mxprsh)
      real(8),intent(in) :: xyz(3), coord(3,natom), znuc(natom)
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: rc, ex12, ex21, ex2i, ex2j, c12, pcxyz(3), extwo
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: ft(0:3), fttwo, fpc2(6)
      real(8) :: pxx, pyy, pzz, pxy, pxz, pyz
      real(8) :: tinv, r1(6), r2(12), r3(10), xx, yy, zz, xy, xz, yz, cint1(3,3,3)
      logical,intent(in) :: iandj
!
      xx= xyz(1)*xyz(1)
      yy= xyz(2)*xyz(2)
      zz= xyz(3)*xyz(3)
      xy= xyz(1)*xyz(2)
      xz= xyz(1)*xyz(3)
      yz= xyz(2)*xyz(3)
      do iatom= 1,natom
        do i= 1,6
          r1(i) = zero
        enddo
        do i= 1,12
          r2(i) = zero
        enddo
        do i= 1,10
          r3(i) = zero
        enddo
        do ijprim= 1,nij
          ex12= exfac(1,ijprim)
          ex21= exfac(2,ijprim)
          ex2i= exfac(3,ijprim)
          ex2j= exfac(4,ijprim)
          c12 = exfac(5,ijprim)
          do i= 1,3
            pcxyz(i)= pijxyz(i,ijprim)-coord(i,iatom)
          enddo
          pxx= pcxyz(1)*pcxyz(1)
          pyy= pcxyz(2)*pcxyz(2)
          pzz= pcxyz(3)*pcxyz(3)
          pxy= pcxyz(1)*pcxyz(2)
          pxz= pcxyz(1)*pcxyz(3)
          pyz= pcxyz(2)*pcxyz(3)
          rc= pxx+pyy+pzz
          tval= ex12*rc
          if(tval >= threshtval) then
            tinv = one/tval
            ft(0)= half*sqrt(pi*tinv)
            ft(1)= half*tinv*ft(0)
            ft(2)= p15 *tinv*ft(1)
            ft(3)= p25 *tinv*ft(2)
          else
            igrid= int(tval)
            tval2= tval *tval
            tval3= tval2*tval
            tval4= tval2*tval2
            tval5= tval2*tval3
            tval6= tval3*tval3
            tval7= tval4*tval3
            tval8= tval4*tval4
            tval9= tval4*tval5
            tval10=tval5*tval5
            do ii= 0,3
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
          endif
!
          do i= 1,3
            ft(i)= ft(i)*znuc(iatom)*c12
          enddo
          extwo= ex12*two
          fttwo= ft(2)*extwo
          ft(3)= ft(3)*extwo
          do i= 1,3
            r1(i)  = r1(i)  +ft(1)*pcxyz(i)*ex2j*ex2i*extwo
            r1(i+3)= r1(i+3)-ft(1)*pcxyz(i)
          enddo
          fpc2(1)= fttwo*pxx-ft(1)
          fpc2(2)= fttwo*pyy-ft(1)
          fpc2(3)= fttwo*pzz-ft(1)
          fpc2(4)= fttwo*pxy
          fpc2(5)= fttwo*pxz
          fpc2(6)= fttwo*pyz
          do i= 1,6
            r2(i)  = r2(i)  -fpc2(i)*ex2j
            r2(i+6)= r2(i+6)+fpc2(i)*ex2i
          enddo
          r3(1) = r3(1) -ft(3)*pxx*pcxyz(1)+ft(2)*pcxyz(1)*three
          r3(2) = r3(2) -ft(3)*pyy*pcxyz(2)+ft(2)*pcxyz(2)*three
          r3(3) = r3(3) -ft(3)*pzz*pcxyz(3)+ft(2)*pcxyz(3)*three
          r3(4) = r3(4) -ft(3)*pxx*pcxyz(2)+ft(2)*pcxyz(2)
          r3(5) = r3(5) -ft(3)*pxx*pcxyz(3)+ft(2)*pcxyz(3)
          r3(6) = r3(6) -ft(3)*pxy*pcxyz(2)+ft(2)*pcxyz(1)
          r3(7) = r3(7) -ft(3)*pyy*pcxyz(3)+ft(2)*pcxyz(3)
          r3(8) = r3(8) -ft(3)*pxz*pcxyz(3)+ft(2)*pcxyz(1)
          r3(9) = r3(9) -ft(3)*pyz*pcxyz(3)+ft(2)*pcxyz(2)
          r3(10)= r3(10)-ft(3)*pxy*pcxyz(3)
        enddo
        cint1(1,1,1)= r3( 1)+r2(1)*xyz(1)+r2( 7)*xyz(1)+r1(1)*xx+r1(4)
        cint1(2,1,1)= r3( 4)+r2(4)*xyz(1)+r2( 7)*xyz(2)+r1(1)*xy
        cint1(3,1,1)= r3( 5)+r2(5)*xyz(1)+r2( 7)*xyz(3)+r1(1)*xz
        cint1(1,2,1)= r3( 4)+r2(1)*xyz(2)+r2(10)*xyz(1)+r1(1)*xy
        cint1(2,2,1)= r3( 6)+r2(4)*xyz(2)+r2(10)*xyz(2)+r1(1)*yy+r1(4)
        cint1(3,2,1)= r3(10)+r2(5)*xyz(2)+r2(10)*xyz(3)+r1(1)*yz
        cint1(1,3,1)= r3( 5)+r2(1)*xyz(3)+r2(11)*xyz(1)+r1(1)*xz
        cint1(2,3,1)= r3(10)+r2(4)*xyz(3)+r2(11)*xyz(2)+r1(1)*yz
        cint1(3,3,1)= r3( 8)+r2(5)*xyz(3)+r2(11)*xyz(3)+r1(1)*zz+r1(4)
        cint1(1,1,2)= r3( 4)+r2(4)*xyz(1)+r2(10)*xyz(1)+r1(2)*xx+r1(5)
        cint1(2,1,2)= r3( 6)+r2(2)*xyz(1)+r2(10)*xyz(2)+r1(2)*xy
        cint1(3,1,2)= r3(10)+r2(6)*xyz(1)+r2(10)*xyz(3)+r1(2)*xz
        cint1(1,2,2)= r3( 6)+r2(4)*xyz(2)+r2( 8)*xyz(1)+r1(2)*xy
        cint1(2,2,2)= r3( 2)+r2(2)*xyz(2)+r2( 8)*xyz(2)+r1(2)*yy+r1(5)
        cint1(3,2,2)= r3( 7)+r2(6)*xyz(2)+r2( 8)*xyz(3)+r1(2)*yz
        cint1(1,3,2)= r3(10)+r2(4)*xyz(3)+r2(12)*xyz(1)+r1(2)*xz
        cint1(2,3,2)= r3( 7)+r2(2)*xyz(3)+r2(12)*xyz(2)+r1(2)*yz
        cint1(3,3,2)= r3( 9)+r2(6)*xyz(3)+r2(12)*xyz(3)+r1(2)*zz+r1(5)
        cint1(1,1,3)= r3( 5)+r2(5)*xyz(1)+r2(11)*xyz(1)+r1(3)*xx+r1(6)
        cint1(2,1,3)= r3(10)+r2(6)*xyz(1)+r2(11)*xyz(2)+r1(3)*xy
        cint1(3,1,3)= r3( 8)+r2(3)*xyz(1)+r2(11)*xyz(3)+r1(3)*xz
        cint1(1,2,3)= r3(10)+r2(5)*xyz(2)+r2(12)*xyz(1)+r1(3)*xy
        cint1(2,2,3)= r3( 7)+r2(6)*xyz(2)+r2(12)*xyz(2)+r1(3)*yy+r1(6)
        cint1(3,2,3)= r3( 9)+r2(3)*xyz(2)+r2(12)*xyz(3)+r1(3)*yz
        cint1(1,3,3)= r3( 8)+r2(5)*xyz(3)+r2( 9)*xyz(1)+r1(3)*xz
        cint1(2,3,3)= r3( 9)+r2(6)*xyz(3)+r2( 9)*xyz(2)+r1(3)*yz
        cint1(3,3,3)= r3( 3)+r2(3)*xyz(3)+r2( 9)*xyz(3)+r1(3)*zz+r1(6)
        if(.not.iandj) then
          do k= 1,3
            do i= 1,3
              do j= 1,3
                egrad(k,iatom)= egrad(k,iatom)+fulldmtrx(locbfij(2)+j,locbfij(1)+i)*cint1(j,i,k)
              enddo
            enddo
          enddo
        else
          do k= 1,3
            do i= 1,3
              do j= 1,3
                egrad(k,iatom)= egrad(k,iatom) &
&                              +fulldmtrx(locbfij(2)+j,locbfij(1)+i)*cint1(j,i,k)*half
              enddo
            enddo
          enddo
        endif
      enddo
!
      return
end


!---------------------------------------------------------------------------
  subroutine int1gcds(egrad,fulldmtrx,exfac,pijxyz,xyz,nij,coord,znuc, &
&                     natom,nao,mxprsh,locbfij,nbfij)
!---------------------------------------------------------------------------
!
! Calculate Helmann-Feynman gradient term <d|V'|s>
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nij, natom, nao, mxprsh, locbfij(2), nbfij(2)
      integer :: ijprim, i, j, k, iatom, igrid, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, two=2.0D+00
      real(8),parameter :: three=3.0D+00, five=5.0D+00, p15=1.5D+00, p25=2.5D+00
      real(8),parameter :: pi=3.141592653589793D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrt3=1.732050807568877D+00
      real(8),intent(in) :: fulldmtrx(nao,nao), exfac(5,mxprsh*mxprsh), pijxyz(3,mxprsh*mxprsh)
      real(8),intent(in) :: xyz(3), coord(3,natom), znuc(natom)
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: rc, ex12, ex21, ex2i, ex2j, c12, pcxyz(3), extwo
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: ft(0:3), fttwo
      real(8) :: pxx, pyy, pzz, pxy, pxz, pyz
      real(8) :: tinv, r1(6), r2(6), r3(10), xx, yy, zz, xy, xz, yz, cint1(6,3), work(6)
!
      xx= xyz(1)*xyz(1)
      yy= xyz(2)*xyz(2)
      zz= xyz(3)*xyz(3)
      xy= xyz(1)*xyz(2)
      xz= xyz(1)*xyz(3)
      yz= xyz(2)*xyz(3)
      do iatom= 1,natom
        do i= 1,6
          r1(i) = zero
        enddo
        do i= 1,6
          r2(i) = zero
        enddo
        do i= 1,10
          r3(i) = zero
        enddo
        do ijprim= 1,nij
          ex12= exfac(1,ijprim)
          ex21= exfac(2,ijprim)
          ex2i= exfac(3,ijprim)
          ex2j= exfac(4,ijprim)
          c12 = exfac(5,ijprim)
          do i= 1,3
            pcxyz(i)= pijxyz(i,ijprim)-coord(i,iatom)
          enddo
          pxx= pcxyz(1)*pcxyz(1)
          pyy= pcxyz(2)*pcxyz(2)
          pzz= pcxyz(3)*pcxyz(3)
          pxy= pcxyz(1)*pcxyz(2)
          pxz= pcxyz(1)*pcxyz(3)
          pyz= pcxyz(2)*pcxyz(3)
          rc= pxx+pyy+pzz
          tval= ex12*rc
          if(tval >= threshtval) then
            tinv = one/tval
            ft(0)= half*sqrt(pi*tinv)
            ft(1)= half*tinv*ft(0)
            ft(2)= p15 *tinv*ft(1)
            ft(3)= p25 *tinv*ft(2)
          else
            igrid= int(tval)
            tval2= tval *tval
            tval3= tval2*tval
            tval4= tval2*tval2
            tval5= tval2*tval3
            tval6= tval3*tval3
            tval7= tval4*tval3
            tval8= tval4*tval4
            tval9= tval4*tval5
            tval10=tval5*tval5
            do ii= 0,3
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
          endif
!
          do i= 1,3
            ft(i)= ft(i)*znuc(iatom)*c12
          enddo
          extwo= ex12*two
          fttwo= ft(2)*extwo
          ft(3)= ft(3)*extwo
          do i= 1,3
            r1(i)  = r1(i)  -ft(1)*pcxyz(i)*ex2i*ex2i*extwo
            r1(i+3)= r1(i+3)-ft(1)*pcxyz(i)
          enddo
          r2(1)= r2(1)+(fttwo*pxx-ft(1))*ex2i
          r2(2)= r2(2)+(fttwo*pyy-ft(1))*ex2i
          r2(3)= r2(3)+(fttwo*pzz-ft(1))*ex2i
          r2(4)= r2(4)+(fttwo*pxy      )*ex2i
          r2(5)= r2(5)+(fttwo*pxz      )*ex2i
          r2(6)= r2(6)+(fttwo*pyz      )*ex2i
          r3( 1)=r3( 1)-ft(3)*pxx*pcxyz(1)+ft(2)*pcxyz(1)*three
          r3( 2)=r3( 2)-ft(3)*pyy*pcxyz(2)+ft(2)*pcxyz(2)*three
          r3( 3)=r3( 3)-ft(3)*pzz*pcxyz(3)+ft(2)*pcxyz(3)*three
          r3( 4)=r3( 4)-ft(3)*pxx*pcxyz(2)+ft(2)*pcxyz(2)
          r3( 5)=r3( 5)-ft(3)*pxx*pcxyz(3)+ft(2)*pcxyz(3)
          r3( 6)=r3( 6)-ft(3)*pxy*pcxyz(2)+ft(2)*pcxyz(1)
          r3( 7)=r3( 7)-ft(3)*pyy*pcxyz(3)+ft(2)*pcxyz(3)
          r3( 8)=r3( 8)-ft(3)*pxz*pcxyz(3)+ft(2)*pcxyz(1)
          r3( 9)=r3( 9)-ft(3)*pyz*pcxyz(3)+ft(2)*pcxyz(2)
          r3(10)=r3(10)-ft(3)*pxy*pcxyz(3)
        enddo
        cint1(1,1)= r3( 1)+r2(1)*xyz(1)*two         +r1(1)*xx+r1(4)
        cint1(2,1)= r3( 4)+r2(4)*xyz(1)+r2(1)*xyz(2)+r1(1)*xy
        cint1(3,1)= r3( 5)+r2(5)*xyz(1)+r2(1)*xyz(3)+r1(1)*xz
        cint1(4,1)= r3( 6)+r2(4)*xyz(2)*two         +r1(1)*yy+r1(4)
        cint1(5,1)= r3(10)+r2(5)*xyz(2)+r2(4)*xyz(3)+r1(1)*yz
        cint1(6,1)= r3( 8)+r2(5)*xyz(3)*two         +r1(1)*zz+r1(4)
        cint1(1,2)= r3( 4)+r2(4)*xyz(1)*two         +r1(2)*xx+r1(5)
        cint1(2,2)= r3( 6)+r2(2)*xyz(1)+r2(4)*xyz(2)+r1(2)*xy
        cint1(3,2)= r3(10)+r2(6)*xyz(1)+r2(4)*xyz(3)+r1(2)*xz
        cint1(4,2)= r3( 2)+r2(2)*xyz(2)*two         +r1(2)*yy+r1(5)
        cint1(5,2)= r3( 7)+r2(6)*xyz(2)+r2(2)*xyz(3)+r1(2)*yz
        cint1(6,2)= r3( 9)+r2(6)*xyz(3)*two         +r1(2)*zz+r1(5)
        cint1(1,3)= r3( 5)+r2(5)*xyz(1)*two         +r1(3)*xx+r1(6)
        cint1(2,3)= r3(10)+r2(6)*xyz(1)+r2(5)*xyz(2)+r1(3)*xy
        cint1(3,3)= r3( 8)+r2(3)*xyz(1)+r2(5)*xyz(3)+r1(3)*xz
        cint1(4,3)= r3( 7)+r2(6)*xyz(2)*two         +r1(3)*yy+r1(6)
        cint1(5,3)= r3( 9)+r2(3)*xyz(2)+r2(6)*xyz(3)+r1(3)*yz
        cint1(6,3)= r3( 3)+r2(3)*xyz(3)*two         +r1(3)*zz+r1(6)
!
        if(nbfij(2) == 6) then
          do k= 1,3
            cint1(2,k)= cint1(2,k)*sqrt3
            cint1(3,k)= cint1(3,k)*sqrt3
            cint1(5,k)= cint1(5,k)*sqrt3
          enddo
        else
          do k= 1,3
            do j= 1,6
              work(j)= cint1(j,k)
            enddo
            cint1(1,k)= work(2)*sqrt3
            cint1(2,k)= work(5)*sqrt3
            cint1(3,k)=(work(6)*two-work(1)-work(4))*half
            cint1(4,k)= work(3)*sqrt3
            cint1(5,k)=(work(1)-work(4))*sqrt3h
          enddo
        endif
        do k= 1,3
          do j= 1,nbfij(2)
            egrad(k,iatom)= egrad(k,iatom)+fulldmtrx(locbfij(2)+j,locbfij(1)+1)*cint1(j,k)
          enddo
        enddo
      enddo
!
      return
end


!---------------------------------------------------------------------------
  subroutine int1gcdp(egrad,fulldmtrx,exfac,pijxyz,xyz,nij,coord,znuc, &
&                     natom,nao,mxprsh,locbfij,nbfij)
!---------------------------------------------------------------------------
!
! Calculate Helmann-Feynman gradient term <d|V'|p>
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nij, natom, nao, mxprsh, locbfij(2), nbfij(2)
      integer :: ijprim, i, j, k, iatom, igrid, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, two=2.0D+00
      real(8),parameter :: three=3.0D+00, five=5.0D+00, six=6.0D+00, seven=7.0D+00
      real(8),parameter :: p15=1.5D+00, p25=2.5D+00, p35=3.5D+00
      real(8),parameter :: pi=3.141592653589793D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrt3=1.732050807568877D+00
      real(8),intent(in) :: fulldmtrx(nao,nao), exfac(5,mxprsh*mxprsh), pijxyz(3,mxprsh*mxprsh)
      real(8),intent(in) :: xyz(3), coord(3,natom), znuc(natom)
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: rc, ex12, ex21, ex2i, ex2j, c12, pcxyz(3), extwo
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: ft(0:4), ft2two, ft3two, ft2inv
      real(8) :: pxx, pyy, pzz, pxy, pxz, pyz
      real(8) :: tinv, r1(9), r2(18), r3(20), r4(15), xx, yy, zz, xy, xz, yz
      real(8) :: fpc2(6), fpc3(10), cint1(6,3,3), work(6)
!
      xx= xyz(1)*xyz(1)
      yy= xyz(2)*xyz(2)
      zz= xyz(3)*xyz(3)
      xy= xyz(1)*xyz(2)
      xz= xyz(1)*xyz(3)
      yz= xyz(2)*xyz(3)
      do iatom= 1,natom
        do i= 1,9
          r1(i) = zero
        enddo
        do i= 1,18
          r2(i) = zero
        enddo
        do i= 1,20
          r3(i) = zero
        enddo
        do i= 1,15
          r4(i) = zero
        enddo
        do ijprim= 1,nij
          ex12= exfac(1,ijprim)
          ex21= exfac(2,ijprim)
          ex2i= exfac(3,ijprim)
          ex2j= exfac(4,ijprim)
          c12 = exfac(5,ijprim)
          do i= 1,3
            pcxyz(i)= pijxyz(i,ijprim)-coord(i,iatom)
          enddo
          pxx= pcxyz(1)*pcxyz(1)
          pyy= pcxyz(2)*pcxyz(2)
          pzz= pcxyz(3)*pcxyz(3)
          pxy= pcxyz(1)*pcxyz(2)
          pxz= pcxyz(1)*pcxyz(3)
          pyz= pcxyz(2)*pcxyz(3)
          rc= pxx+pyy+pzz
          tval= ex12*rc
          if(tval >= threshtval) then
            tinv = one/tval
            ft(0)= half*sqrt(pi*tinv)
            ft(1)= half*tinv*ft(0)
            ft(2)= p15 *tinv*ft(1)
            ft(3)= p25 *tinv*ft(2)
            ft(4)= p35 *tinv*ft(3)
          else
            igrid= int(tval)
            tval2= tval *tval
            tval3= tval2*tval
            tval4= tval2*tval2
            tval5= tval2*tval3
            tval6= tval3*tval3
            tval7= tval4*tval3
            tval8= tval4*tval4
            tval9= tval4*tval5
            tval10=tval5*tval5
            do ii= 0,4
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
          endif
!
          do i= 1,4
            ft(i)= ft(i)*znuc(iatom)*c12
          enddo
          extwo= ex12*two
          ex21 = ex21*half
          ft2two= ft(2)*extwo
          ft3two= ft(3)*extwo
          ft(4) = ft(4)*extwo
          ft2inv= ft(2)*ex21
          do i= 1,3
            r1(i)  = r1(i)  +ft(1)*pcxyz(i)*ex2j
            r1(i+3)= r1(i+3)-ft(1)*pcxyz(i)*ex2i
            r1(i+6)= r1(i+6)+ft(1)*pcxyz(i)*ex2j*ex2i*ex2i*extwo
          enddo
          fpc2(1)= ft2two*pxx-ft(1)
          fpc2(2)= ft2two*pyy-ft(1)
          fpc2(3)= ft2two*pzz-ft(1)
          fpc2(4)= ft2two*pxy
          fpc2(5)= ft2two*pxz
          fpc2(6)= ft2two*pyz
          do i= 1,6
            r2(i   )= r2(i   )-fpc2(i)*ex2j*ex2i
            r2(i+ 6)= r2(i+ 6)+fpc2(i)*ex2i*ex2i
            r2(i+12)= r2(i+12)+fpc2(i)*ex21
          enddo
          fpc3(1) =-ft3two*pxx*pcxyz(1)+ft(2)*pcxyz(1)*three
          fpc3(2) =-ft3two*pyy*pcxyz(2)+ft(2)*pcxyz(2)*three
          fpc3(3) =-ft3two*pzz*pcxyz(3)+ft(2)*pcxyz(3)*three
          fpc3(4) =-ft3two*pxx*pcxyz(2)+ft(2)*pcxyz(2)
          fpc3(5) =-ft3two*pxx*pcxyz(3)+ft(2)*pcxyz(3)
          fpc3(6) =-ft3two*pxy*pcxyz(2)+ft(2)*pcxyz(1)
          fpc3(7) =-ft3two*pyy*pcxyz(3)+ft(2)*pcxyz(3)
          fpc3(8) =-ft3two*pxz*pcxyz(3)+ft(2)*pcxyz(1)
          fpc3(9) =-ft3two*pyz*pcxyz(3)+ft(2)*pcxyz(2)
          fpc3(10)=-ft3two*pxy*pcxyz(3)
          do i= 1,10
            r3(i   )= r3(i   )-fpc3(i)*ex2j
            r3(i+10)= r3(i+10)+fpc3(i)*ex2i
          enddo
          r4( 1)= r4( 1)+ft(4)*pxx*pxx-ft(3)*pxx*six      +ft2inv*three
          r4( 2)= r4( 2)+ft(4)*pyy*pyy-ft(3)*pyy*six      +ft2inv*three
          r4( 3)= r4( 3)+ft(4)*pzz*pzz-ft(3)*pzz*six      +ft2inv*three
          r4( 4)= r4( 4)+ft(4)*pxx*pxy-ft(3)*pxy*three
          r4( 5)= r4( 5)+ft(4)*pxx*pxz-ft(3)*pxz*three
          r4( 6)= r4( 6)+ft(4)*pxy*pyy-ft(3)*pxy*three
          r4( 7)= r4( 7)+ft(4)*pyy*pyz-ft(3)*pyz*three
          r4( 8)= r4( 8)+ft(4)*pxz*pzz-ft(3)*pxz*three
          r4( 9)= r4( 9)+ft(4)*pyz*pzz-ft(3)*pyz*three
          r4(10)= r4(10)+ft(4)*pxx*pyy-ft(3)*pxx-ft(3)*pyy+ft2inv
          r4(11)= r4(11)+ft(4)*pxx*pzz-ft(3)*pxx-ft(3)*pzz+ft2inv
          r4(12)= r4(12)+ft(4)*pyy*pzz-ft(3)*pyy-ft(3)*pzz+ft2inv
          r4(13)= r4(13)+ft(4)*pxx*pyz-ft(3)*pyz
          r4(14)= r4(14)+ft(4)*pxy*pyz-ft(3)*pxz
          r4(15)= r4(15)+ft(4)*pxy*pzz-ft(3)*pxy
        enddo
        cint1(1,1,1)=r4( 1)+r3( 1)*xyz(1)+r3(11)*xyz(1)*two+r2( 1)*xx*two+r2( 7)*xx &
&                   +r2(13)*three+r1(1)*xyz(1)+r1(4)*xyz(1)*two+r1(7)*xx*xyz(1)
        cint1(2,1,1)=r4( 4)+r3( 4)*xyz(1)+r3(14)*xyz(1)+r3(11)*xyz(2)+r2( 4)*xx &
&                   +r2( 1)*xy+r2( 7)*xy+r2(16)+r1(4)*xyz(2)+r1(7)*xy*xyz(1)
        cint1(3,1,1)=r4( 5)+r3( 5)*xyz(1)+r3(15)*xyz(1)+r3(11)*xyz(3)+r2( 5)*xx &
&                   +r2( 1)*xz+r2( 7)*xz+r2(17)+r1(4)*xyz(3)+r1(7)*xz*xyz(1)
        cint1(4,1,1)=r4(10)+r3( 6)*xyz(1)+r3(14)*xyz(2)*two+r2( 4)*xy*two+r2( 7)*yy &
&                   +r2(13)+r1(1)*xyz(1)+r1(7)*yy*xyz(1)
        cint1(5,1,1)=r4(13)+r3(10)*xyz(1)+r3(15)*xyz(2)+r3(14)*xyz(3)+r2( 5)*xy &
&                   +r2( 4)*xz+r2( 7)*yz+r1(7)*yz*xyz(1)
        cint1(6,1,1)=r4(11)+r3( 8)*xyz(1)+r3(15)*xyz(3)*two+r2( 5)*xz*two+r2( 7)*zz &
&                   +r2(13)+r1(1)*xyz(1)+r1(7)*zz*xyz(1)
        cint1(1,2,1)=r4( 4)+r3( 1)*xyz(2)+r3(14)*xyz(1)*two+r2( 1)*xy*two+r2(10)*xx &
&                   +r2(16)+r1(1)*xyz(2)+r1(7)*xx*xyz(2)
        cint1(2,2,1)=r4(10)+r3( 4)*xyz(2)+r3(16)*xyz(1)+r3(14)*xyz(2)+r2( 4)*xy &
&                   +r2( 1)*yy+r2(10)*xy+r2(13)+r1(4)*xyz(1)+r1(7)*xy*xyz(2)
        cint1(3,2,1)=r4(13)+r3( 5)*xyz(2)+r3(20)*xyz(1)+r3(14)*xyz(3)+r2( 5)*xy &
&                   +r2( 1)*yz+r2(10)*xz+r1(7)*xz*xyz(2)
        cint1(4,2,1)=r4( 6)+r3( 6)*xyz(2)+r3(16)*xyz(2)*two+r2( 4)*yy*two+r2(10)*yy &
&                   +r2(16)*three+r1(1)*xyz(2)+r1(4)*xyz(2)*two+r1(7)*yy*xyz(2)
        cint1(5,2,1)=r4(14)+r3(10)*xyz(2)+r3(20)*xyz(2)+r3(16)*xyz(3)+r2( 5)*yy &
&                   +r2( 4)*yz+r2(10)*yz+r2(17)+r1(4)*xyz(3)+r1(7)*yz*xyz(2)
        cint1(6,2,1)=r4(15)+r3( 8)*xyz(2)+r3(20)*xyz(3)*two+r2( 5)*yz*two+r2(10)*zz &
&                   +r2(16)+r1(1)*xyz(2)+r1(7)*zz*xyz(2)
        cint1(1,3,1)=r4( 5)+r3( 1)*xyz(3)+r3(15)*xyz(1)*two+r2( 1)*xz*two+r2(11)*xx &
&                   +r2(17)+r1(1)*xyz(3)+r1(7)*xx*xyz(3)
        cint1(2,3,1)=r4(13)+r3( 4)*xyz(3)+r3(20)*xyz(1)+r3(15)*xyz(2)+r2( 4)*xz &
&                   +r2( 1)*yz+r2(11)*xy+r1(7)*xy*xyz(3)
        cint1(3,3,1)=r4(11)+r3( 5)*xyz(3)+r3(18)*xyz(1)+r3(15)*xyz(3)+r2( 5)*xz &
&                   +r2( 1)*zz+r2(11)*xz+r2(13)+r1(4)*xyz(1)+r1(7)*xz*xyz(3)
        cint1(4,3,1)=r4(14)+r3( 6)*xyz(3)+r3(20)*xyz(2)*two+r2( 4)*yz*two+r2(11)*yy &
&                   +r2(17)+r1(1)*xyz(3)+r1(7)*yy*xyz(3)
        cint1(5,3,1)=r4(15)+r3(10)*xyz(3)+r3(18)*xyz(2)+r3(20)*xyz(3)+r2( 5)*yz &
&                   +r2( 4)*zz+r2(11)*yz+r2(16)+r1(4)*xyz(2)+r1(7)*yz*xyz(3)
        cint1(6,3,1)=r4( 8)+r3( 8)*xyz(3)+r3(18)*xyz(3)*two+r2( 5)*zz*two+r2(11)*zz &
&                   +r2(17)*three+r1(1)*xyz(3)+r1(4)*xyz(3)*two+r1(7)*zz*xyz(3)
        cint1(1,1,2)=r4( 4)+r3( 4)*xyz(1)+r3(14)*xyz(1)*two+r2( 4)*xx*two+r2(10)*xx &
&                   +r2(16)*three+r1(2)*xyz(1)+r1(5)*xyz(1)*two+r1(8)*xx*xyz(1)
        cint1(2,1,2)=r4(10)+r3( 6)*xyz(1)+r3(16)*xyz(1)+r3(14)*xyz(2)+r2( 2)*xx &
&                   +r2( 4)*xy+r2(10)*xy+r2(14)+r1(5)*xyz(2)+r1(8)*xy*xyz(1)
        cint1(3,1,2)=r4(13)+r3(10)*xyz(1)+r3(20)*xyz(1)+r3(14)*xyz(3)+r2( 6)*xx &
&                   +r2( 4)*xz+r2(10)*xz+r2(18)+r1(5)*xyz(3)+r1(8)*xz*xyz(1)
        cint1(4,1,2)=r4( 6)+r3( 2)*xyz(1)+r3(16)*xyz(2)*two+r2( 2)*xy*two+r2(10)*yy &
&                   +r2(16)+r1(2)*xyz(1)+r1(8)*yy*xyz(1)
        cint1(5,1,2)=r4(14)+r3( 7)*xyz(1)+r3(20)*xyz(2)+r3(16)*xyz(3)+r2( 6)*xy &
&                   +r2( 2)*xz+r2(10)*yz+r1(8)*yz*xyz(1)
        cint1(6,1,2)=r4(15)+r3( 9)*xyz(1)+r3(20)*xyz(3)*two+r2( 6)*xz*two+r2(10)*zz &
&                   +r2(16)+r1(2)*xyz(1)+r1(8)*zz*xyz(1)
        cint1(1,2,2)=r4(10)+r3( 4)*xyz(2)+r3(16)*xyz(1)*two+r2( 4)*xy*two+r2( 8)*xx &
&                   +r2(14)+r1(2)*xyz(2)+r1(8)*xx*xyz(2)
        cint1(2,2,2)=r4( 6)+r3( 6)*xyz(2)+r3(12)*xyz(1)+r3(16)*xyz(2)+r2( 2)*xy &
&                   +r2( 4)*yy+r2( 8)*xy+r2(16)+r1(5)*xyz(1)+r1(8)*xy*xyz(2)
        cint1(3,2,2)=r4(14)+r3(10)*xyz(2)+r3(17)*xyz(1)+r3(16)*xyz(3)+r2( 6)*xy &
&                   +r2( 4)*yz+r2( 8)*xz+r1(8)*xz*xyz(2)
        cint1(4,2,2)=r4( 2)+r3( 2)*xyz(2)+r3(12)*xyz(2)*two+r2( 2)*yy*two+r2( 8)*yy &
&                   +r2(14)*three+r1(2)*xyz(2)+r1(5)*xyz(2)*two+r1(8)*yy*xyz(2)
        cint1(5,2,2)=r4( 7)+r3( 7)*xyz(2)+r3(17)*xyz(2)+r3(12)*xyz(3)+r2( 6)*yy &
&                   +r2( 2)*yz+r2( 8)*yz+r2(18)+r1(5)*xyz(3)+r1(8)*yz*xyz(2)
        cint1(6,2,2)=r4(12)+r3( 9)*xyz(2)+r3(17)*xyz(3)*two+r2( 6)*yz*two+r2( 8)*zz &
&                   +r2(14)+r1(2)*xyz(2)+r1(8)*zz*xyz(2)
        cint1(1,3,2)=r4(13)+r3( 4)*xyz(3)+r3(20)*xyz(1)*two+r2( 4)*xz*two+r2(12)*xx &
&                   +r2(18)+r1(2)*xyz(3)+r1(8)*xx*xyz(3)
        cint1(2,3,2)=r4(14)+r3( 6)*xyz(3)+r3(17)*xyz(1)+r3(20)*xyz(2)+r2( 2)*xz &
&                   +r2( 4)*yz+r2(12)*xy+r1(8)*xy*xyz(3)
        cint1(3,3,2)=r4(15)+r3(10)*xyz(3)+r3(19)*xyz(1)+r3(20)*xyz(3)+r2( 6)*xz &
&                   +r2( 4)*zz+r2(12)*xz+r2(16)+r1(5)*xyz(1)+r1(8)*xz*xyz(3)
        cint1(4,3,2)=r4( 7)+r3( 2)*xyz(3)+r3(17)*xyz(2)*two+r2( 2)*yz*two+r2(12)*yy &
&                   +r2(18)+r1(2)*xyz(3)+r1(8)*yy*xyz(3)
        cint1(5,3,2)=r4(12)+r3( 7)*xyz(3)+r3(19)*xyz(2)+r3(17)*xyz(3)+r2( 6)*yz &
&                   +r2( 2)*zz+r2(12)*yz+r2(14)+r1(5)*xyz(2)+r1(8)*yz*xyz(3)
        cint1(6,3,2)=r4( 9)+r3( 9)*xyz(3)+r3(19)*xyz(3)*two+r2( 6)*zz*two+r2(12)*zz &
&                   +r2(18)*three+r1(2)*xyz(3)+r1(5)*xyz(3)*two+r1(8)*zz*xyz(3)
        cint1(1,1,3)=r4( 5)+r3( 5)*xyz(1)+r3(15)*xyz(1)*two+r2( 5)*xx*two+r2(11)*xx &
&                   +r2(17)*three+r1(3)*xyz(1)+r1(6)*xyz(1)*two+r1(9)*xx*xyz(1)
        cint1(2,1,3)=r4(13)+r3(10)*xyz(1)+r3(20)*xyz(1)+r3(15)*xyz(2)+r2( 6)*xx &
&                   +r2( 5)*xy+r2(11)*xy+r2(18)+r1(6)*xyz(2)+r1(9)*xy*xyz(1)
        cint1(3,1,3)=r4(11)+r3( 8)*xyz(1)+r3(18)*xyz(1)+r3(15)*xyz(3)+r2( 3)*xx &
&                   +r2( 5)*xz+r2(11)*xz+r2(15)+r1(6)*xyz(3)+r1(9)*xz*xyz(1)
        cint1(4,1,3)=r4(14)+r3( 7)*xyz(1)+r3(20)*xyz(2)*two+r2( 6)*xy*two+r2(11)*yy &
&                   +r2(17)+r1(3)*xyz(1)+r1(9)*yy*xyz(1)
        cint1(5,1,3)=r4(15)+r3( 9)*xyz(1)+r3(18)*xyz(2)+r3(20)*xyz(3)+r2( 3)*xy &
&                   +r2( 6)*xz+r2(11)*yz+r1(9)*yz*xyz(1)
        cint1(6,1,3)=r4( 8)+r3( 3)*xyz(1)+r3(18)*xyz(3)*two+r2( 3)*xz*two+r2(11)*zz &
&                   +r2(17)+r1(3)*xyz(1)+r1(9)*zz*xyz(1)
        cint1(1,2,3)=r4(13)+r3( 5)*xyz(2)+r3(20)*xyz(1)*two+r2( 5)*xy*two+r2(12)*xx &
&                   +r2(18)+r1(3)*xyz(2)+r1(9)*xx*xyz(2)
        cint1(2,2,3)=r4(14)+r3(10)*xyz(2)+r3(17)*xyz(1)+r3(20)*xyz(2)+r2( 6)*xy &
&                   +r2( 5)*yy+r2(12)*xy+r2(17)+r1(6)*xyz(1)+r1(9)*xy*xyz(2)
        cint1(3,2,3)=r4(15)+r3( 8)*xyz(2)+r3(19)*xyz(1)+r3(20)*xyz(3)+r2( 3)*xy &
&                   +r2( 5)*yz+r2(12)*xz+r1(9)*xz*xyz(2)
        cint1(4,2,3)=r4( 7)+r3( 7)*xyz(2)+r3(17)*xyz(2)*two+r2( 6)*yy*two+r2(12)*yy &
&                   +r2(18)*three+r1(3)*xyz(2)+r1(6)*xyz(2)*two+r1(9)*yy*xyz(2)
        cint1(5,2,3)=r4(12)+r3( 9)*xyz(2)+r3(19)*xyz(2)+r3(17)*xyz(3)+r2( 3)*yy &
&                   +r2( 6)*yz+r2(12)*yz+r2(15)+r1(6)*xyz(3)+r1(9)*yz*xyz(2)
        cint1(6,2,3)=r4( 9)+r3( 3)*xyz(2)+r3(19)*xyz(3)*two+r2( 3)*yz*two+r2(12)*zz &
&                   +r2(18)+r1(3)*xyz(2)+r1(9)*zz*xyz(2)
        cint1(1,3,3)=r4(11)+r3( 5)*xyz(3)+r3(18)*xyz(1)*two+r2( 5)*xz*two+r2( 9)*xx &
&                   +r2(15)+r1(3)*xyz(3)+r1(9)*xx*xyz(3)
        cint1(2,3,3)=r4(15)+r3(10)*xyz(3)+r3(19)*xyz(1)+r3(18)*xyz(2)+r2( 6)*xz &
&                   +r2( 5)*yz+r2( 9)*xy+r1(9)*xy*xyz(3)
        cint1(3,3,3)=r4( 8)+r3( 8)*xyz(3)+r3(13)*xyz(1)+r3(18)*xyz(3)+r2( 3)*xz &
&                   +r2( 5)*zz+r2( 9)*xz+r2(17)+r1(6)*xyz(1)+r1(9)*xz*xyz(3)
        cint1(4,3,3)=r4(12)+r3( 7)*xyz(3)+r3(19)*xyz(2)*two+r2( 6)*yz*two+r2( 9)*yy &
&                   +r2(15)+r1(3)*xyz(3)+r1(9)*yy*xyz(3)
        cint1(5,3,3)=r4( 9)+r3( 9)*xyz(3)+r3(13)*xyz(2)+r3(19)*xyz(3)+r2( 3)*yz &
&                   +r2( 6)*zz+r2( 9)*yz+r2(18)+r1(6)*xyz(2)+r1(9)*yz*xyz(3)
        cint1(6,3,3)=r4( 3)+r3( 3)*xyz(3)+r3(13)*xyz(3)*two+r2( 3)*zz*two+r2( 9)*zz &
&                    +r2(15)*three+r1(3)*xyz(3)+r1(6)*xyz(3)*two+r1(9)*zz*xyz(3)
!
        if(nbfij(2) == 6) then
          do k= 1,3
            do i= 1,3
              cint1(2,i,k)= cint1(2,i,k)*sqrt3
              cint1(3,i,k)= cint1(3,i,k)*sqrt3
              cint1(5,i,k)= cint1(5,i,k)*sqrt3
            enddo
          enddo
        else
          do k= 1,3
            do i= 1,3
              do j= 1,6
                work(j)= cint1(j,i,k)
              enddo
              cint1(1,i,k)= work(2)*sqrt3
              cint1(2,i,k)= work(5)*sqrt3
              cint1(3,i,k)=(work(6)*two-work(1)-work(4))*half
              cint1(4,i,k)= work(3)*sqrt3
              cint1(5,i,k)=(work(1)-work(4))*sqrt3h
            enddo
          enddo
        endif
        do k= 1,3
          do i= 1,3
            do j= 1,nbfij(2)
              egrad(k,iatom)= egrad(k,iatom)+fulldmtrx(locbfij(2)+j,locbfij(1)+i)*cint1(j,i,k)
            enddo
          enddo
        enddo
      enddo
!
      return
end


!---------------------------------------------------------------------------
  subroutine int1gcdd(egrad,fulldmtrx,exfac,pijxyz,xyz,nij,coord,znuc, &
&                     natom,nao,mxprsh,locbfij,nbfij,iandj)
!---------------------------------------------------------------------------
!
! Calculate Helmann-Feynman gradient term <d|V'|d>
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nij, natom, nao, mxprsh, locbfij(2), nbfij(2)
      integer :: ijprim, i, j, k, iatom, igrid, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, two=2.0D+00, three=3.0D+00
      real(8),parameter :: four=4.0D+00, five=5.0D+00, six=6.0D+00, seven=7.0D+00, nine=9.0D+00
      real(8),parameter :: ten=1.0D+01, fift=1.5D+01
      real(8),parameter :: p15=1.5D+00, p25=2.5D+00, p35=3.5D+00, p45=4.5D+00
      real(8),parameter :: pi=3.141592653589793D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrt3=1.732050807568877D+00
      real(8),intent(in) :: fulldmtrx(nao,nao), exfac(5,mxprsh*mxprsh), pijxyz(3,mxprsh*mxprsh)
      real(8),intent(in) :: xyz(3), coord(3,natom), znuc(natom)
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: rc, ex12, ex21, ex2i, ex2j, c12, pcxyz(3), extwo, ex2ii, ex2ij, ex2jj
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: ft(0:5), ft2two, ft3two, ft4two
      real(8) :: ft2inv, ft3inv, pxx, pyy, pzz, pxy, pxz, pyz
      real(8) :: pxxx, pyyy, pzzz, pxxy, pxxz, pxyy, pyyz, pxzz, pyzz, pxyz
      real(8) :: tinv, r1(15), r2(24), r3(40), r4(30), r5(21), xx, yy, zz, xy, xz, yz
      real(8) :: xxx, yyy, zzz, xxy, xxz, xyy, yyz, xzz, yzz, xyz2
      real(8) :: fpc2(6), fpc3(10), fpc4(15), cint1(6,6,3), work(6)
      logical,intent(in) :: iandj
!
      xx= xyz(1)*xyz(1)
      yy= xyz(2)*xyz(2)
      zz= xyz(3)*xyz(3)
      xy= xyz(1)*xyz(2)
      xz= xyz(1)*xyz(3)
      yz= xyz(2)*xyz(3)
      xxx= xx*xyz(1)
      yyy= yy*xyz(2)
      zzz= zz*xyz(3)
      xxy= xx*xyz(2)
      xxz= xx*xyz(3)
      xyy= yy*xyz(1)
      yyz= yy*xyz(3)
      xzz= zz*xyz(1)
      yzz= zz*xyz(2)
      xyz2= xy*xyz(3)
      do iatom= 1,natom
        do i= 1,15
          r1(i) = zero
        enddo
        do i= 1,24
          r2(i) = zero
        enddo
        do i= 1,40
          r3(i) = zero
        enddo
        do i= 1,30
          r4(i) = zero
        enddo
        do i= 1,21
          r5(i) = zero
        enddo
        do ijprim= 1,nij
          ex12= exfac(1,ijprim)
          ex21= exfac(2,ijprim)
          ex2i= exfac(3,ijprim)
          ex2j= exfac(4,ijprim)
          c12 = exfac(5,ijprim)
          do i= 1,3
            pcxyz(i)= pijxyz(i,ijprim)-coord(i,iatom)
          enddo
          pxx= pcxyz(1)*pcxyz(1)
          pyy= pcxyz(2)*pcxyz(2)
          pzz= pcxyz(3)*pcxyz(3)
          pxy= pcxyz(1)*pcxyz(2)
          pxz= pcxyz(1)*pcxyz(3)
          pyz= pcxyz(2)*pcxyz(3)
          pxxx= pxx*pcxyz(1)
          pyyy= pyy*pcxyz(2)
          pzzz= pzz*pcxyz(3)
          pxxy= pxx*pcxyz(2)
          pxxz= pxx*pcxyz(3)
          pxyy= pyy*pcxyz(1)
          pyyz= pyy*pcxyz(3)
          pxzz= pzz*pcxyz(1)
          pyzz= pzz*pcxyz(2)
          pxyz= pxy*pcxyz(3)
          rc= pxx+pyy+pzz
          tval= ex12*rc
          if(tval >= threshtval) then
            tinv = one/tval
            ft(0)= half*sqrt(pi*tinv)
            ft(1)= half*tinv*ft(0)
            ft(2)= p15 *tinv*ft(1)
            ft(3)= p25 *tinv*ft(2)
            ft(4)= p35 *tinv*ft(3)
            ft(5)= p45 *tinv*ft(4)
          else
            igrid= int(tval)
            tval2= tval *tval
            tval3= tval2*tval
            tval4= tval2*tval2
            tval5= tval2*tval3
            tval6= tval3*tval3
            tval7= tval4*tval3
            tval8= tval4*tval4
            tval9= tval4*tval5
            tval10=tval5*tval5
            do ii= 0,5
              ft(ii)= fgrid(0,ii,igrid)      +fgrid( 1,ii,igrid)*tval  +fgrid( 2,ii,igrid)*tval2 &
&                    +fgrid(3,ii,igrid)*tval3+fgrid( 4,ii,igrid)*tval4 +fgrid( 5,ii,igrid)*tval5 &
&                    +fgrid(6,ii,igrid)*tval6+fgrid( 7,ii,igrid)*tval7 +fgrid( 8,ii,igrid)*tval8 &
&                    +fgrid(9,ii,igrid)*tval9+fgrid(10,ii,igrid)*tval10
            enddo
          endif
!
          do i= 1,5
            ft(i)= ft(i)*znuc(iatom)*c12
          enddo
          extwo = ex12*two
          ex21  = ex21*half
          ft2two= ft(2)*extwo
          ft3two= ft(3)*extwo
          ft4two= ft(4)*extwo
          ft(5) = ft(5)*extwo
          ft2inv= ft(2)*ex21
          ft3inv= ft(3)*ex21
          ex2ii= ex2i*ex2i
          ex2ij= ex2i*ex2j
          ex2jj= ex2j*ex2j
          do i= 1,3
            r1(i   )= r1(i   )-ft(1)*pcxyz(i)*ex2jj*ex2ii*extwo
            r1(i+ 3)= r1(i+ 3)-ft(1)*pcxyz(i)*ex2jj
            r1(i+ 6)= r1(i+ 6)+ft(1)*pcxyz(i)*ex2ij
            r1(i+ 9)= r1(i+ 9)-ft(1)*pcxyz(i)*ex2ii
            r1(i+12)= r1(i+12)-ft(1)*pcxyz(i)*ex21
          enddo
          fpc2(1)= ft2two*pxx-ft(1)
          fpc2(2)= ft2two*pyy-ft(1)
          fpc2(3)= ft2two*pzz-ft(1)
          fpc2(4)= ft2two*pxy
          fpc2(5)= ft2two*pxz
          fpc2(6)= ft2two*pyz
          do i= 1,6
            r2(i   )= r2(i   )+fpc2(i)*ex2jj*ex2i
            r2(i+ 6)= r2(i+ 6)-fpc2(i)*ex2ii*ex2j
            r2(i+12)= r2(i+12)-fpc2(i)*ex2j*ex21
            r2(i+18)= r2(i+18)+fpc2(i)*ex2i*ex21
          enddo
          fpc3(1) =-ft3two*pxxx+ft(2)*pcxyz(1)*three
          fpc3(2) =-ft3two*pyyy+ft(2)*pcxyz(2)*three
          fpc3(3) =-ft3two*pzzz+ft(2)*pcxyz(3)*three
          fpc3(4) =-ft3two*pxxy+ft(2)*pcxyz(2)
          fpc3(5) =-ft3two*pxxz+ft(2)*pcxyz(3)
          fpc3(6) =-ft3two*pxyy+ft(2)*pcxyz(1)
          fpc3(7) =-ft3two*pyyz+ft(2)*pcxyz(3)
          fpc3(8) =-ft3two*pxzz+ft(2)*pcxyz(1)
          fpc3(9) =-ft3two*pyzz+ft(2)*pcxyz(2)
          fpc3(10)=-ft3two*pxyz
          do i= 1,10
            r3(i   )= r3(i   )+fpc3(i)*ex2jj
            r3(i+10)= r3(i+10)-fpc3(i)*ex2ij
            r3(i+20)= r3(i+20)+fpc3(i)*ex2ii
            r3(i+30)= r3(i+30)+fpc3(i)*ex21
          enddo
          fpc4( 1)= ft4two*pxx*pxx-ft(3)*pxx*six      +ft2inv*three
          fpc4( 2)= ft4two*pyy*pyy-ft(3)*pyy*six      +ft2inv*three
          fpc4( 3)= ft4two*pzz*pzz-ft(3)*pzz*six      +ft2inv*three
          fpc4( 4)= ft4two*pxx*pxy-ft(3)*pxy*three
          fpc4( 5)= ft4two*pxx*pxz-ft(3)*pxz*three
          fpc4( 6)= ft4two*pxy*pyy-ft(3)*pxy*three
          fpc4( 7)= ft4two*pyy*pyz-ft(3)*pyz*three
          fpc4( 8)= ft4two*pxz*pzz-ft(3)*pxz*three
          fpc4( 9)= ft4two*pyz*pzz-ft(3)*pyz*three
          fpc4(10)= ft4two*pxx*pyy-ft(3)*pxx-ft(3)*pyy+ft2inv
          fpc4(11)= ft4two*pxx*pzz-ft(3)*pxx-ft(3)*pzz+ft2inv
          fpc4(12)= ft4two*pyy*pzz-ft(3)*pyy-ft(3)*pzz+ft2inv
          fpc4(13)= ft4two*pxx*pyz-ft(3)*pyz
          fpc4(14)= ft4two*pxy*pyz-ft(3)*pxz
          fpc4(15)= ft4two*pxy*pzz-ft(3)*pxy
          do i= 1,15
            r4(i   )= r4(i   )-fpc4(i)*ex2j
            r4(i+15)= r4(i+15)+fpc4(i)*ex2i
          enddo
          r5( 1)= r5( 1)-ft(5)*pxxx*pxx+ft(4)*pxxx*ten             -ft3inv*pcxyz(1)*fift
          r5( 2)= r5( 2)-ft(5)*pyyy*pyy+ft(4)*pyyy*ten             -ft3inv*pcxyz(2)*fift
          r5( 3)= r5( 3)-ft(5)*pzzz*pzz+ft(4)*pzzz*ten             -ft3inv*pcxyz(3)*fift
          r5( 4)= r5( 4)-ft(5)*pxxx*pxy+ft(4)*pxxy*six             -ft3inv*pcxyz(2)*three
          r5( 5)= r5( 5)-ft(5)*pxxx*pxz+ft(4)*pxxz*six             -ft3inv*pcxyz(3)*three
          r5( 6)= r5( 6)-ft(5)*pxyy*pyy+ft(4)*pxyy*six             -ft3inv*pcxyz(1)*three
          r5( 7)= r5( 7)-ft(5)*pyyy*pyz+ft(4)*pyyz*six             -ft3inv*pcxyz(3)*three
          r5( 8)= r5( 8)-ft(5)*pxzz*pzz+ft(4)*pxzz*six             -ft3inv*pcxyz(1)*three
          r5( 9)= r5( 9)-ft(5)*pyzz*pzz+ft(4)*pyzz*six             -ft3inv*pcxyz(2)*three
          r5(10)= r5(10)-ft(5)*pxxx*pyy+ft(4)*pxyy*three+ft(4)*pxxx-ft3inv*pcxyz(1)*three
          r5(11)= r5(11)-ft(5)*pxxx*pyz+ft(4)*pxyz*three
          r5(12)= r5(12)-ft(5)*pxxx*pzz+ft(4)*pxzz*three+ft(4)*pxxx-ft3inv*pcxyz(1)*three
          r5(13)= r5(13)-ft(5)*pxxy*pyy+ft(4)*pxxy*three+ft(4)*pyyy-ft3inv*pcxyz(2)*three
          r5(14)= r5(14)-ft(5)*pxyy*pyz+ft(4)*pxyz*three
          r5(15)= r5(15)-ft(5)*pyyy*pzz+ft(4)*pyzz*three+ft(4)*pyyy-ft3inv*pcxyz(2)*three
          r5(16)= r5(16)-ft(5)*pxxz*pzz+ft(4)*pxxz*three+ft(4)*pzzz-ft3inv*pcxyz(3)*three
          r5(17)= r5(17)-ft(5)*pxyz*pzz+ft(4)*pxyz*three
          r5(18)= r5(18)-ft(5)*pyyz*pzz+ft(4)*pyyz*three+ft(4)*pzzz-ft3inv*pcxyz(3)*three
          r5(19)= r5(19)-ft(5)*pxxy*pyz+ft(4)*pxxz      +ft(4)*pyyz-ft3inv*pcxyz(3)
          r5(20)= r5(20)-ft(5)*pxxy*pzz+ft(4)*pxxy      +ft(4)*pyzz-ft3inv*pcxyz(2)
          r5(21)= r5(21)-ft(5)*pxyy*pzz+ft(4)*pxyy      +ft(4)*pxzz-ft3inv*pcxyz(1)
        enddo
        cint1(1,1,1)=r5( 1)+r4( 1)*xyz(1)*two+r4(16)*xyz(1)*two+r3( 1)*xx+r3(11)*xx*four &
&                   +r3(21)*xx+r3(31)*six+r2( 1)*xxx*two+r2( 7)*xxx*two+r2(13)*xyz(1)*six &
&                   +r2(19)*xyz(1)*six+r1( 1)*xx*xx+r1( 4)*xx+r1( 7)*xx*four+r1(10)*xx &
&                   +r1(13)*three
        cint1(2,1,1)=r5( 4)+r4( 4)*xyz(1)*two+r4(19)*xyz(1)+r4(16)*xyz(2)+r3( 4)*xx &
&                   +r3(14)*xx*two+r3(11)*xy*two+r3(21)*xy+r3(34)*three+r2( 4)*xxx &
&                   +r2( 1)*xxy+r2( 7)*xxy*two+r2(16)*xyz(1)*two+r2(22)*xyz(1) &
&                   +r2(19)*xyz(2)*three+r1( 1)*xx*xy+r1( 7)*xy*two+r1(10)*xy
        cint1(3,1,1)=r5( 5)+r4( 5)*xyz(1)*two+r4(20)*xyz(1)+r4(16)*xyz(3)+r3( 5)*xx &
&                   +r3(15)*xx*two+r3(11)*xz*two+r3(21)*xz+r3(35)*three+r2( 5)*xxx &
&                   +r2( 1)*xxz+r2( 7)*xxz*two+r2(17)*xyz(1)*two+r2(23)*xyz(1) &
&                   +r2(19)*xyz(3)*three+r1( 1)*xx*xz+r1( 7)*xz*two+r1(10)*xz
        cint1(4,1,1)=r5(10)+r4(10)*xyz(1)*two+r4(19)*xyz(2)*two+r3( 6)*xx+r3(14)*xy*four &
&                   +r3(21)*yy+r3(31)+r3(36)+r2( 4)*xxy*two+r2( 7)*xyy*two+r2(13)*xyz(1)*two &
&                   +r2(22)*xyz(2)*two+r1( 1)*xx*yy+r1( 4)*xx+r1(10)*yy+r1(13)
        cint1(5,1,1)=r5(11)+r4(13)*xyz(1)*two+r4(20)*xyz(2)+r4(19)*xyz(3)+r3(10)*xx &
&                   +r3(15)*xy*two+r3(14)*xz*two+r3(21)*yz+r3(40)+r2( 5)*xxy+r2( 4)*xxz &
&                   +r2( 7)*xyz2*two+r2(23)*xyz(2)+r2(22)*xyz(3)+r1( 1)*xx*yz+r1(10)*yz
        cint1(6,1,1)=r5(12)+r4(11)*xyz(1)*two+r4(20)*xyz(3)*two+r3( 8)*xx+r3(15)*xz*four &
&                   +r3(21)*zz+r3(31)+r3(38)+r2( 5)*xxz*two+r2( 7)*xzz*two+r2(13)*xyz(1)*two &
&                   +r2(23)*xyz(3)*two+r1( 1)*xx*zz+r1( 4)*xx+r1(10)*zz+r1(13)
        cint1(1,2,1)=r5( 4)+r4( 4)*xyz(1)+r4( 1)*xyz(2)+r4(19)*xyz(1)*two+r3( 1)*xy &
&                   +r3(14)*xx*two+r3(11)*xy*two+r3(24)*xx+r3(34)*three+r2( 1)*xxy*two &
&                   +r2(10)*xxx+r2( 7)*xxy+r2(16)*xyz(1)+r2(13)*xyz(2)*three &
&                   +r2(22)*xyz(1)*two+r1( 1)*xy*xx+r1( 4)*xy+r1( 7)*xy*two
        cint1(2,2,1)=r5(10)+r4(10)*xyz(1)+r4( 4)*xyz(2)+r4(25)*xyz(1)+r4(19)*xyz(2) &
&                   +r3( 4)*xy+r3(16)*xx+r3(14)*xy+r3(14)*xy+r3(11)*yy+r3(24)*xy+r3(31) &
&                   +r3(36)+r2( 4)*xxy+r2( 1)*xyy+r2(10)*xxy+r2( 7)*xyy+r2(13)*xyz(1) &
&                   +r2(16)*xyz(2)+r2(19)*xyz(1)+r2(22)*xyz(2)+r1( 1)*xy*xy+r1( 7)*xx &
&                   +r1( 7)*yy+r1(13)
        cint1(3,2,1)=r5(11)+r4(13)*xyz(1)+r4( 5)*xyz(2)+r4(28)*xyz(1)+r4(19)*xyz(3) &
&                   +r3( 5)*xy+r3(20)*xx+r3(14)*xz+r3(15)*xy+r3(11)*yz+r3(24)*xz+r3(40) &
&                   +r2( 5)*xxy+r2( 1)*xyz2+r2(10)*xxz+r2( 7)*xyz2+r2(17)*xyz(2) &
&                   +r2(22)*xyz(3)+r1( 1)*xy*xz+r1( 7)*yz
        cint1(4,2,1)=r5(13)+r4( 6)*xyz(1)+r4(10)*xyz(2)+r4(25)*xyz(2)*two+r3( 6)*xy &
&                   +r3(16)*xy*two+r3(14)*yy*two+r3(24)*yy+r3(34)*three+r2( 4)*xyy*two &
&                   +r2(10)*xyy+r2( 7)*yyy+r2(16)*xyz(1)*three+r2(13)*xyz(2) &
&                   +r2(19)*xyz(2)*two+r1( 1)*xy*yy+r1( 4)*xy+r1( 7)*xy*two
        cint1(5,2,1)=r5(19)+r4(14)*xyz(1)+r4(13)*xyz(2)+r4(28)*xyz(2)+r4(25)*xyz(3) &
&                   +r3(10)*xy+r3(20)*xy+r3(16)*xz+r3(15)*yy+r3(14)*yz+r3(24)*yz+r3(35) &
&                   +r2( 5)*xyy+r2( 4)*xyz2+r2(10)*xyz2+r2( 7)*yyz+r2(17)*xyz(1) &
&                   +r2(19)*xyz(3)+r1( 1)*xy*yz+r1( 7)*xz
        cint1(6,2,1)=r5(20)+r4(15)*xyz(1)+r4(11)*xyz(2)+r4(28)*xyz(3)*two+r3( 8)*xy &
&                   +r3(20)*xz*two+r3(15)*yz*two+r3(24)*zz+r3(34)+r2( 5)*xyz2*two+r2(10)*xzz &
&                   +r2( 7)*yzz+r2(16)*xyz(1)+r2(13)*xyz(2)+r1( 1)*xy*zz+r1( 4)*xy
        cint1(1,3,1)=r5( 5)+r4( 5)*xyz(1)+r4( 1)*xyz(3)+r4(20)*xyz(1)*two+r3( 1)*xz &
&                   +r3(15)*xx*two+r3(11)*xz*two+r3(25)*xx+r3(35)*three+r2( 1)*xxz*two &
&                   +r2(11)*xxx+r2( 7)*xxz+r2(17)*xyz(1)+r2(13)*xyz(3)*three &
&                   +r2(23)*xyz(1)*two+r1( 1)*xz*xx+r1( 4)*xz+r1( 7)*xz*two
        cint1(2,3,1)=r5(11)+r4(13)*xyz(1)+r4( 4)*xyz(3)+r4(28)*xyz(1)+r4(20)*xyz(2) &
&                   +r3( 4)*xz+r3(20)*xx+r3(15)*xy+r3(14)*xz+r3(11)*yz+r3(25)*xy+r3(40) &
&                   +r2( 4)*xxz+r2( 1)*xyz2+r2(11)*xxy+r2( 7)*xyz2+r2(16)*xyz(3) &
&                   +r2(23)*xyz(2)+r1( 1)*xz*xy+r1( 7)*yz
        cint1(3,3,1)=r5(12)+r4(11)*xyz(1)+r4( 5)*xyz(3)+r4(26)*xyz(1)+r4(20)*xyz(3) &
&                   +r3( 5)*xz+r3(18)*xx+r3(15)*xz+r3(15)*xz+r3(11)*zz+r3(25)*xz+r3(31) &
&                   +r3(38)+r2( 5)*xxz+r2( 1)*xzz+r2(11)*xxz+r2( 7)*xzz+r2(13)*xyz(1) &
&                   +r2(17)*xyz(3)+r2(19)*xyz(1)+r2(23)*xyz(3)+r1( 1)*xz*xz+r1( 7)*xx &
&                   +r1( 7)*zz+r1(13)
        cint1(4,3,1)=r5(19)+r4(14)*xyz(1)+r4(10)*xyz(3)+r4(28)*xyz(2)*two+r3( 6)*xz &
&                   +r3(20)*xy*two+r3(14)*yz*two+r3(25)*yy+r3(35)+r2( 4)*xyz2*two+r2(11)*xyy &
&                   +r2( 7)*yyz+r2(17)*xyz(1)+r2(13)*xyz(3)+r1( 1)*xz*yy+r1( 4)*xz
        cint1(5,3,1)=r5(20)+r4(15)*xyz(1)+r4(13)*xyz(3)+r4(26)*xyz(2)+r4(28)*xyz(3) &
&                   +r3(10)*xz+r3(18)*xy+r3(20)*xz+r3(15)*yz+r3(14)*zz+r3(25)*yz+r3(34) &
&                   +r2( 5)*xyz2+r2( 4)*xzz+r2(11)*xyz2+r2( 7)*yzz+r2(16)*xyz(1) &
&                   +r2(19)*xyz(2)+r1( 1)*xz*yz+r1( 7)*xy
        cint1(6,3,1)=r5(16)+r4( 8)*xyz(1)+r4(11)*xyz(3)+r4(26)*xyz(3)*two+r3( 8)*xz &
&                   +r3(18)*xz*two+r3(15)*zz*two+r3(25)*zz+r3(35)*three+r2( 5)*xzz*two &
&                   +r2(11)*xzz+r2( 7)*zzz+r2(17)*xyz(1)*three+r2(13)*xyz(3) &
&                   +r2(19)*xyz(3)*two+r1( 1)*xz*zz+r1( 4)*xz+r1( 7)*xz*two
        cint1(1,4,1)=r5(10)+r4( 4)*xyz(2)*two+r4(25)*xyz(1)*two+r3( 1)*yy+r3(14)*xy*four &
&                   +r3(26)*xx+r3(31)+r3(36)+r2( 1)*xyy*two+r2(10)*xxy*two+r2(16)*xyz(2)*two &
&                   +r2(19)*xyz(1)*two+r1( 1)*yy*xx+r1( 4)*yy+r1(10)*xx+r1(13)
        cint1(2,4,1)=r5(13)+r4(10)*xyz(2)*two+r4(21)*xyz(1)+r4(25)*xyz(2)+r3( 4)*yy &
&                   +r3(16)*xy*two+r3(14)*yy*two+r3(26)*xy+r3(34)*three+r2( 4)*xyy &
&                   +r2( 1)*yyy+r2(10)*xyy*two+r2(13)*xyz(2)*two+r2(22)*xyz(1)*three &
&                   +r2(19)*xyz(2)+r1( 1)*yy*xy+r1( 7)*xy*two+r1(10)*xy
        cint1(3,4,1)=r5(19)+r4(13)*xyz(2)*two+r4(29)*xyz(1)+r4(25)*xyz(3)+r3( 5)*yy &
&                   +r3(20)*xy*two+r3(14)*yz*two+r3(26)*xz+r3(35)+r2( 5)*xyy+r2( 1)*yyz &
&                   +r2(10)*xyz2*two+r2(23)*xyz(1)+r2(19)*xyz(3)+r1( 1)*yy*xz+r1(10)*xz
        cint1(4,4,1)=r5( 6)+r4( 6)*xyz(2)*two+r4(21)*xyz(2)*two+r3( 6)*yy+r3(16)*yy*four &
&                   +r3(26)*yy+r3(36)*six+r2( 4)*yyy*two+r2(10)*yyy*two+r2(16)*xyz(2)*six &
&                   +r2(22)*xyz(2)*six+r1( 1)*yy*yy+r1( 4)*yy+r1( 7)*yy*four+r1(10)*yy &
&                   +r1(13)*three
        cint1(5,4,1)=r5(14)+r4(14)*xyz(2)*two+r4(29)*xyz(2)+r4(21)*xyz(3)+r3(10)*yy &
&                   +r3(20)*yy*two+r3(16)*yz*two+r3(26)*yz+r3(40)*three+r2( 5)*yyy &
&                   +r2( 4)*yyz+r2(10)*yyz*two+r2(17)*xyz(2)*two+r2(23)*xyz(2) &
&                   +r2(22)*xyz(3)*three+r1( 1)*yy*yz+r1( 7)*yz*two+r1(10)*yz
        cint1(6,4,1)=r5(21)+r4(15)*xyz(2)*two+r4(29)*xyz(3)*two+r3( 8)*yy+r3(20)*yz*four &
&                   +r3(26)*zz+r3(36)+r3(38)+r2( 5)*yyz*two+r2(10)*yzz*two+r2(16)*xyz(2)*two &
&                   +r2(23)*xyz(3)*two+r1( 1)*yy*zz+r1( 4)*yy+r1(10)*zz+r1(13)
        cint1(1,5,1)=r5(11)+r4( 5)*xyz(2)+r4( 4)*xyz(3)+r4(28)*xyz(1)*two+r3( 1)*yz &
&                   +r3(15)*xy*two+r3(14)*xz*two+r3(30)*xx+r3(40)+r2( 1)*xyz2*two &
&                   +r2(11)*xxy+r2(10)*xxz+r2(17)*xyz(2)+r2(16)*xyz(3)+r1( 1)*yz*xx+r1( 4)*yz
        cint1(2,5,1)=r5(19)+r4(13)*xyz(2)+r4(10)*xyz(3)+r4(29)*xyz(1)+r4(28)*xyz(2) &
&                   +r3( 4)*yz+r3(20)*xy+r3(15)*yy+r3(16)*xz+r3(14)*yz+r3(30)*xy+r3(35) &
&                   +r2( 4)*xyz2+r2( 1)*yyz+r2(11)*xyy+r2(10)*xyz2+r2(13)*xyz(3) &
&                   +r2(23)*xyz(1)+r1( 1)*yz*xy+r1( 7)*xz
        cint1(3,5,1)=r5(20)+r4(11)*xyz(2)+r4(13)*xyz(3)+r4(30)*xyz(1)+r4(28)*xyz(3) &
&                   +r3( 5)*yz+r3(18)*xy+r3(15)*yz+r3(20)*xz+r3(14)*zz+r3(30)*xz+r3(34) &
&                   +r2( 5)*xyz2+r2( 1)*yzz+r2(11)*xyz2+r2(10)*xzz+r2(13)*xyz(2) &
&                   +r2(22)*xyz(1)+r1( 1)*yz*xz+r1( 7)*xy
        cint1(4,5,1)=r5(14)+r4(14)*xyz(2)+r4( 6)*xyz(3)+r4(29)*xyz(2)*two+r3( 6)*yz &
&                   +r3(20)*yy*two+r3(16)*yz*two+r3(30)*yy+r3(40)*three+r2( 4)*yyz*two &
&                   +r2(11)*yyy+r2(10)*yyz+r2(17)*xyz(2)+r2(16)*xyz(3)*three+ &
&                   r2(23)*xyz(2)*two+r1( 1)*yz*yy+r1( 4)*yz+r1( 7)*yz*two
        cint1(5,5,1)=r5(21)+r4(15)*xyz(2)+r4(14)*xyz(3)+r4(30)*xyz(2)+r4(29)*xyz(3) &
&                   +r3(10)*yz+r3(18)*yy+r3(20)*yz+r3(20)*yz+r3(16)*zz+r3(30)*yz+r3(36) &
&                   +r3(38)+r2( 5)*yyz+r2( 4)*yzz+r2(11)*yyz+r2(10)*yzz+r2(16)*xyz(2) &
&                   +r2(17)*xyz(3)+r2(22)*xyz(2)+r2(23)*xyz(3)+r1( 1)*yz*yz+r1( 7)*yy &
&                   +r1( 7)*zz+r1(13)
        cint1(6,5,1)=r5(17)+r4( 8)*xyz(2)+r4(15)*xyz(3)+r4(30)*xyz(3)*two+r3( 8)*yz &
&                   +r3(18)*yz*two+r3(20)*zz*two+r3(30)*zz+r3(40)*three+r2( 5)*yzz*two &
&                   +r2(11)*yzz+r2(10)*zzz+r2(17)*xyz(2)*three+r2(16)*xyz(3) &
&                   +r2(22)*xyz(3)*two+r1( 1)*yz*zz+r1( 4)*yz+r1( 7)*yz*two
        cint1(1,6,1)=r5(12)+r4( 5)*xyz(3)*two+r4(26)*xyz(1)*two+r3( 1)*zz+r3(15)*xz*four &
&                   +r3(28)*xx+r3(31)+r3(38)+r2( 1)*xzz*two+r2(11)*xxz*two+r2(17)*xyz(3)*two &
&                   +r2(19)*xyz(1)*two+r1( 1)*zz*xx+r1( 4)*zz+r1(10)*xx+r1(13)
        cint1(2,6,1)=r5(20)+r4(13)*xyz(3)*two+r4(30)*xyz(1)+r4(26)*xyz(2)+r3( 4)*zz &
&                   +r3(20)*xz*two+r3(15)*yz*two+r3(28)*xy+r3(34)+r2( 4)*xzz+r2( 1)*yzz &
&                   +r2(11)*xyz2*two+r2(22)*xyz(1)+r2(19)*xyz(2)+r1( 1)*zz*xy+r1(10)*xy
        cint1(3,6,1)=r5(16)+r4(11)*xyz(3)*two+r4(23)*xyz(1)+r4(26)*xyz(3)+r3( 5)*zz &
&                   +r3(18)*xz*two+r3(15)*zz*two+r3(28)*xz+r3(35)*three+r2( 5)*xzz &
&                   +r2( 1)*zzz+r2(11)*xzz*two+r2(13)*xyz(3)*two+r2(23)*xyz(1)*three &
&                   +r2(19)*xyz(3)+r1( 1)*zz*xz+r1( 7)*xz*two+r1(10)*xz
        cint1(4,6,1)=r5(21)+r4(14)*xyz(3)*two+r4(30)*xyz(2)*two+r3( 6)*zz+r3(20)*yz*four &
&                   +r3(28)*yy+r3(36)+r3(38)+r2( 4)*yzz*two+r2(11)*yyz*two+r2(17)*xyz(3)*two &
&                   +r2(22)*xyz(2)*two+r1( 1)*zz*yy+r1( 4)*zz+r1(10)*yy+r1(13)
        cint1(5,6,1)=r5(17)+r4(15)*xyz(3)*two+r4(23)*xyz(2)+r4(30)*xyz(3)+r3(10)*zz &
&                   +r3(18)*yz*two+r3(20)*zz*two+r3(28)*yz+r3(40)*three+r2( 5)*yzz &
&                   +r2( 4)*zzz+r2(11)*yzz*two+r2(16)*xyz(3)*two+r2(23)*xyz(2)*three &
&                   +r2(22)*xyz(3)+r1( 1)*zz*yz+r1( 7)*yz*two+r1(10)*yz
        cint1(6,6,1)=r5( 8)+r4( 8)*xyz(3)*two+r4(23)*xyz(3)*two+r3( 8)*zz+r3(18)*zz*four &
&                   +r3(28)*zz+r3(38)*six+r2( 5)*zzz*two+r2(11)*zzz*two+r2(17)*xyz(3)*six &
&                   +r2(23)*xyz(3)*six+r1( 1)*zz*zz+r1( 4)*zz+r1( 7)*zz*four+r1(10)*zz &
&                   +r1(13)*three
        cint1(1,1,2)=r5( 4)+r4( 4)*xyz(1)*two+r4(19)*xyz(1)*two+r3( 4)*xx+r3(14)*xx*four &
&                   +r3(24)*xx+r3(34)*six+r2( 4)*xxx*two+r2(10)*xxx*two+r2(16)*xyz(1)*six &
&                   +r2(22)*xyz(1)*six+r1( 2)*xx*xx+r1( 5)*xx+r1( 8)*xx*four+r1(11)*xx &
&                   +r1(14)*three
        cint1(2,1,2)=r5(10)+r4(10)*xyz(1)*two+r4(25)*xyz(1)+r4(19)*xyz(2)+r3( 6)*xx &
&                   +r3(16)*xx*two+r3(14)*xy*two+r3(24)*xy+r3(36)*three+r2( 2)*xxx &
&                   +r2( 4)*xxy+r2(10)*xxy*two+r2(14)*xyz(1)*two+r2(20)*xyz(1) &
&                   +r2(22)*xyz(2)*three+r1( 2)*xx*xy+r1( 8)*xy*two+r1(11)*xy
        cint1(3,1,2)=r5(11)+r4(13)*xyz(1)*two+r4(28)*xyz(1)+r4(19)*xyz(3)+r3(10)*xx &
&                   +r3(20)*xx*two+r3(14)*xz*two+r3(24)*xz+r3(40)*three+r2( 6)*xxx &
&                   +r2( 4)*xxz+r2(10)*xxz*two+r2(18)*xyz(1)*two+r2(24)*xyz(1) &
&                   +r2(22)*xyz(3)*three+r1( 2)*xx*xz+r1( 8)*xz*two+r1(11)*xz
        cint1(4,1,2)=r5(13)+r4( 6)*xyz(1)*two+r4(25)*xyz(2)*two+r3( 2)*xx+r3(16)*xy*four &
&                   +r3(24)*yy+r3(32)+r3(34)+r2( 2)*xxy*two+r2(10)*xyy*two+r2(16)*xyz(1)*two &
&                   +r2(20)*xyz(2)*two+r1( 2)*xx*yy+r1( 5)*xx+r1(11)*yy+r1(14)
        cint1(5,1,2)=r5(19)+r4(14)*xyz(1)*two+r4(28)*xyz(2)+r4(25)*xyz(3)+r3( 7)*xx &
&                   +r3(20)*xy*two+r3(16)*xz*two+r3(24)*yz+r3(37)+r2( 6)*xxy+r2( 2)*xxz &
&                   +r2(10)*xyz2*two+r2(24)*xyz(2)+r2(20)*xyz(3)+r1( 2)*xx*yz+r1(11)*yz
        cint1(6,1,2)=r5(20)+r4(15)*xyz(1)*two+r4(28)*xyz(3)*two+r3( 9)*xx+r3(20)*xz*four &
&                   +r3(24)*zz+r3(34)+r3(39)+r2( 6)*xxz*two+r2(10)*xzz*two+r2(16)*xyz(1)*two &
&                   +r2(24)*xyz(3)*two+r1( 2)*xx*zz+r1( 5)*xx+r1(11)*zz+r1(14)
        cint1(1,2,2)=r5(10)+r4(10)*xyz(1)+r4( 4)*xyz(2)+r4(25)*xyz(1)*two+r3( 4)*xy &
&                   +r3(16)*xx*two+r3(14)*xy*two+r3(26)*xx+r3(36)*three+r2( 4)*xxy*two &
&                   +r2( 8)*xxx+r2(10)*xxy+r2(14)*xyz(1)+r2(16)*xyz(2)*three &
&                   +r2(20)*xyz(1)*two+r1( 2)*xy*xx+r1( 5)*xy+r1( 8)*xy*two
        cint1(2,2,2)=r5(13)+r4( 6)*xyz(1)+r4(10)*xyz(2)+r4(21)*xyz(1)+r4(25)*xyz(2) &
&                   +r3( 6)*xy+r3(12)*xx+r3(16)*xy+r3(16)*xy+r3(14)*yy+r3(26)*xy+r3(32) &
&                   +r3(34)+r2( 2)*xxy+r2( 4)*xyy+r2( 8)*xxy+r2(10)*xyy+r2(16)*xyz(1) &
&                   +r2(14)*xyz(2)+r2(22)*xyz(1)+r2(20)*xyz(2)+r1( 2)*xy*xy+r1( 8)*xx &
&                   +r1( 8)*yy+r1(14)
        cint1(3,2,2)=r5(19)+r4(14)*xyz(1)+r4(13)*xyz(2)+r4(29)*xyz(1)+r4(25)*xyz(3) &
&                   +r3(10)*xy+r3(17)*xx+r3(16)*xz+r3(20)*xy+r3(14)*yz+r3(26)*xz+r3(37) &
&                   +r2( 6)*xxy+r2( 4)*xyz2+r2( 8)*xxz+r2(10)*xyz2+r2(18)*xyz(2)+r2(20)*xyz(3) &
&                   +r1( 2)*xy*xz+r1( 8)*yz
        cint1(4,2,2)=r5( 6)+r4( 2)*xyz(1)+r4( 6)*xyz(2)+r4(21)*xyz(2)*two+r3( 2)*xy &
&                   +r3(12)*xy*two+r3(16)*yy*two+r3(26)*yy+r3(36)*three+r2( 2)*xyy*two &
&                   +r2( 8)*xyy+r2(10)*yyy+r2(14)*xyz(1)*three+r2(16)*xyz(2) &
&                   +r2(22)*xyz(2)*two+r1( 2)*xy*yy+r1( 5)*xy+r1( 8)*xy*two
        cint1(5,2,2)=r5(14)+r4( 7)*xyz(1)+r4(14)*xyz(2)+r4(29)*xyz(2)+r4(21)*xyz(3) &
&                   +r3( 7)*xy+r3(17)*xy+r3(12)*xz+r3(20)*yy+r3(16)*yz+r3(26)*yz+r3(40) &
&                   +r2( 6)*xyy+r2( 2)*xyz2+r2( 8)*xyz2+r2(10)*yyz+r2(18)*xyz(1)+r2(22)*xyz(3) &
&                   +r1( 2)*xy*yz+r1( 8)*xz
        cint1(6,2,2)=r5(21)+r4(12)*xyz(1)+r4(15)*xyz(2)+r4(29)*xyz(3)*two+r3( 9)*xy &
&                   +r3(17)*xz*two+r3(20)*yz*two+r3(26)*zz+r3(36)+r2( 6)*xyz2*two+r2( 8)*xzz &
&                   +r2(10)*yzz+r2(14)*xyz(1)+r2(16)*xyz(2)+r1( 2)*xy*zz+r1( 5)*xy
        cint1(1,3,2)=r5(11)+r4(13)*xyz(1)+r4( 4)*xyz(3)+r4(28)*xyz(1)*two+r3( 4)*xz &
&                   +r3(20)*xx*two+r3(14)*xz*two+r3(30)*xx+r3(40)*three+r2( 4)*xxz*two &
&                   +r2(12)*xxx+r2(10)*xxz+r2(18)*xyz(1)+r2(16)*xyz(3)*three &
&                   +r2(24)*xyz(1)*two+r1( 2)*xz*xx+r1( 5)*xz+r1( 8)*xz*two
        cint1(2,3,2)=r5(19)+r4(14)*xyz(1)+r4(10)*xyz(3)+r4(29)*xyz(1)+r4(28)*xyz(2) &
&                   +r3( 6)*xz+r3(17)*xx+r3(20)*xy+r3(16)*xz+r3(14)*yz+r3(30)*xy+r3(37) &
&                   +r2( 2)*xxz+r2( 4)*xyz2+r2(12)*xxy+r2(10)*xyz2+r2(14)*xyz(3) &
&                   +r2(24)*xyz(2)+r1( 2)*xz*xy+r1( 8)*yz
        cint1(3,3,2)=r5(20)+r4(15)*xyz(1)+r4(13)*xyz(3)+r4(30)*xyz(1)+r4(28)*xyz(3) &
&                   +r3(10)*xz+r3(19)*xx+r3(20)*xz+r3(20)*xz+r3(14)*zz+r3(30)*xz+r3(34) &
&                   +r3(39)+r2( 6)*xxz+r2( 4)*xzz+r2(12)*xxz+r2(10)*xzz+r2(16)*xyz(1) &
&                   +r2(18)*xyz(3)+r2(22)*xyz(1)+r2(24)*xyz(3)+r1( 2)*xz*xz+r1( 8)*xx &
&                   +r1( 8)*zz+r1(14)
        cint1(4,3,2)=r5(14)+r4( 7)*xyz(1)+r4( 6)*xyz(3)+r4(29)*xyz(2)*two+r3( 2)*xz &
&                   +r3(17)*xy*two+r3(16)*yz*two+r3(30)*yy+r3(40)+r2( 2)*xyz2*two+r2(12)*xyy &
&                   +r2(10)*yyz+r2(18)*xyz(1)+r2(16)*xyz(3)+r1( 2)*xz*yy+r1( 5)*xz
        cint1(5,3,2)=r5(21)+r4(12)*xyz(1)+r4(14)*xyz(3)+r4(30)*xyz(2)+r4(29)*xyz(3) &
&                   +r3( 7)*xz+r3(19)*xy+r3(17)*xz+r3(20)*yz+r3(16)*zz+r3(30)*yz+r3(36) &
&                   +r2( 6)*xyz2+r2( 2)*xzz+r2(12)*xyz2+r2(10)*yzz+r2(14)*xyz(1) &
&                   +r2(22)*xyz(2)+r1( 2)*xz*yz+r1( 8)*xy
        cint1(6,3,2)=r5(17)+r4( 9)*xyz(1)+r4(15)*xyz(3)+r4(30)*xyz(3)*two+r3( 9)*xz &
&                   +r3(19)*xz*two+r3(20)*zz*two+r3(30)*zz+r3(40)*three+r2( 6)*xzz*two &
&                   +r2(12)*xzz+r2(10)*zzz+r2(18)*xyz(1)*three+r2(16)*xyz(3) &
&                   +r2(22)*xyz(3)*two+r1( 2)*xz*zz+r1( 5)*xz+r1( 8)*xz*two
        cint1(1,4,2)=r5(13)+r4(10)*xyz(2)*two+r4(21)*xyz(1)*two+r3( 4)*yy+r3(16)*xy*four &
&                   +r3(22)*xx+r3(32)+r3(34)+r2( 4)*xyy*two+r2( 8)*xxy*two+r2(14)*xyz(2)*two &
&                   +r2(22)*xyz(1)*two+r1( 2)*yy*xx+r1( 5)*yy+r1(11)*xx+r1(14)
        cint1(2,4,2)=r5( 6)+r4( 6)*xyz(2)*two+r4(17)*xyz(1)+r4(21)*xyz(2)+r3( 6)*yy &
&                   +r3(12)*xy*two+r3(16)*yy*two+r3(22)*xy+r3(36)*three+r2( 2)*xyy &
&                   +r2( 4)*yyy+r2( 8)*xyy*two+r2(16)*xyz(2)*two+r2(20)*xyz(1)*three &
&                   +r2(22)*xyz(2)+r1( 2)*yy*xy+r1( 8)*xy*two+r1(11)*xy
        cint1(3,4,2)=r5(14)+r4(14)*xyz(2)*two+r4(22)*xyz(1)+r4(21)*xyz(3)+r3(10)*yy &
&                   +r3(17)*xy*two+r3(16)*yz*two+r3(22)*xz+r3(40)+r2( 6)*xyy+r2( 4)*yyz &
&                   +r2( 8)*xyz2*two+r2(24)*xyz(1)+r2(22)*xyz(3)+r1( 2)*yy*xz+r1(11)*xz
        cint1(4,4,2)=r5( 2)+r4( 2)*xyz(2)*two+r4(17)*xyz(2)*two+r3( 2)*yy+r3(12)*yy*four &
&                   +r3(22)*yy+r3(32)*six+r2( 2)*yyy*two+r2( 8)*yyy*two+r2(14)*xyz(2)*six &
&                   +r2(20)*xyz(2)*six+r1( 2)*yy*yy+r1( 5)*yy+r1( 8)*yy*four+r1(11)*yy &
&                   +r1(14)*three
        cint1(5,4,2)=r5( 7)+r4( 7)*xyz(2)*two+r4(22)*xyz(2)+r4(17)*xyz(3)+r3( 7)*yy &
&                   +r3(17)*yy*two+r3(12)*yz*two+r3(22)*yz+r3(37)*three+r2( 6)*yyy &
&                   +r2( 2)*yyz+r2( 8)*yyz*two+r2(18)*xyz(2)*two+r2(24)*xyz(2) &
&                   +r2(20)*xyz(3)*three+r1( 2)*yy*yz+r1( 8)*yz*two+r1(11)*yz
        cint1(6,4,2)=r5(15)+r4(12)*xyz(2)*two+r4(22)*xyz(3)*two+r3( 9)*yy+r3(17)*yz*four &
&                   +r3(22)*zz+r3(32)+r3(39)+r2( 6)*yyz*two+r2( 8)*yzz*two+r2(14)*xyz(2)*two &
&                   +r2(24)*xyz(3)*two+r1( 2)*yy*zz+r1( 5)*yy+r1(11)*zz+r1(14)
        cint1(1,5,2)=r5(19)+r4(13)*xyz(2)+r4(10)*xyz(3)+r4(29)*xyz(1)*two+r3( 4)*yz &
&                   +r3(20)*xy*two+r3(16)*xz*two+r3(27)*xx+r3(37)+r2( 4)*xyz2*two+r2(12)*xxy &
&                   +r2( 8)*xxz+r2(18)*xyz(2)+r2(14)*xyz(3)+r1( 2)*yz*xx+r1( 5)*yz
        cint1(2,5,2)=r5(14)+r4(14)*xyz(2)+r4( 6)*xyz(3)+r4(22)*xyz(1)+r4(29)*xyz(2) &
&                   +r3( 6)*yz+r3(17)*xy+r3(20)*yy+r3(12)*xz+r3(16)*yz+r3(27)*xy+r3(40) &
&                   +r2( 2)*xyz2+r2( 4)*yyz+r2(12)*xyy+r2( 8)*xyz2+r2(16)*xyz(3) &
&                   +r2(24)*xyz(1)+r1( 2)*yz*xy+r1( 8)*xz
        cint1(3,5,2)=r5(21)+r4(15)*xyz(2)+r4(14)*xyz(3)+r4(27)*xyz(1)+r4(29)*xyz(3) &
&                   +r3(10)*yz+r3(19)*xy+r3(20)*yz+r3(17)*xz+r3(16)*zz+r3(27)*xz+r3(36) &
&                   +r2( 6)*xyz2+r2( 4)*yzz+r2(12)*xyz2+r2( 8)*xzz+r2(16)*xyz(2) &
&                   +r2(20)*xyz(1)+r1( 2)*yz*xz+r1( 8)*xy
        cint1(4,5,2)=r5( 7)+r4( 7)*xyz(2)+r4( 2)*xyz(3)+r4(22)*xyz(2)*two+r3( 2)*yz &
&                   +r3(17)*yy*two+r3(12)*yz*two+r3(27)*yy+r3(37)*three+r2( 2)*yyz*two &
&                   +r2(12)*yyy+r2( 8)*yyz+r2(18)*xyz(2)+r2(14)*xyz(3)*three &
&                   +r2(24)*xyz(2)*two+r1( 2)*yz*yy+r1( 5)*yz+r1( 8)*yz*two
        cint1(5,5,2)=r5(15)+r4(12)*xyz(2)+r4( 7)*xyz(3)+r4(27)*xyz(2)+r4(22)*xyz(3) &
&                   +r3( 7)*yz+r3(19)*yy+r3(17)*yz+r3(17)*yz+r3(12)*zz+r3(27)*yz+r3(32) &
&                   +r3(39)+r2( 6)*yyz+r2( 2)*yzz+r2(12)*yyz+r2( 8)*yzz+r2(14)*xyz(2) &
&                   +r2(18)*xyz(3)+r2(20)*xyz(2)+r2(24)*xyz(3)+r1( 2)*yz*yz+r1( 8)*yy &
&                   +r1( 8)*zz+r1(14)
        cint1(6,5,2)=r5(18)+r4( 9)*xyz(2)+r4(12)*xyz(3)+r4(27)*xyz(3)*two+r3( 9)*yz &
&                   +r3(19)*yz*two+r3(17)*zz*two+r3(27)*zz+r3(37)*three+r2( 6)*yzz*two &
&                   +r2(12)*yzz+r2( 8)*zzz+r2(18)*xyz(2)*three+r2(14)*xyz(3) &
&                   +r2(20)*xyz(3)*two+r1( 2)*yz*zz+r1( 5)*yz+r1( 8)*yz*two
        cint1(1,6,2)=r5(20)+r4(13)*xyz(3)*two+r4(30)*xyz(1)*two+r3( 4)*zz+r3(20)*xz*four &
&                   +r3(29)*xx+r3(34)+r3(39)+r2( 4)*xzz*two+r2(12)*xxz*two+r2(18)*xyz(3)*two &
&                   +r2(22)*xyz(1)*two+r1( 2)*zz*xx+r1( 5)*zz+r1(11)*xx+r1(14)
        cint1(2,6,2)=r5(21)+r4(14)*xyz(3)*two+r4(27)*xyz(1)+r4(30)*xyz(2)+r3( 6)*zz &
&                   +r3(17)*xz*two+r3(20)*yz*two+r3(29)*xy+r3(36)+r2( 2)*xzz+r2( 4)*yzz &
&                   +r2(12)*xyz2*two+r2(20)*xyz(1)+r2(22)*xyz(2)+r1( 2)*zz*xy+r1(11)*xy
        cint1(3,6,2)=r5(17)+r4(15)*xyz(3)*two+r4(24)*xyz(1)+r4(30)*xyz(3)+r3(10)*zz &
&                   +r3(19)*xz*two+r3(20)*zz*two+r3(29)*xz+r3(40)*three+r2( 6)*xzz &
&                   +r2( 4)*zzz+r2(12)*xzz*two+r2(16)*xyz(3)*two+r2(24)*xyz(1)*three &
&                   +r2(22)*xyz(3)+r1( 2)*zz*xz+r1( 8)*xz*two+r1(11)*xz
        cint1(4,6,2)=r5(15)+r4( 7)*xyz(3)*two+r4(27)*xyz(2)*two+r3( 2)*zz+r3(17)*yz*four &
&                   +r3(29)*yy+r3(32)+r3(39)+r2( 2)*yzz*two+r2(12)*yyz*two+r2(18)*xyz(3)*two &
&                   +r2(20)*xyz(2)*two+r1( 2)*zz*yy+r1( 5)*zz+r1(11)*yy+r1(14)
        cint1(5,6,2)=r5(18)+r4(12)*xyz(3)*two+r4(24)*xyz(2)+r4(27)*xyz(3)+r3( 7)*zz &
&                   +r3(19)*yz*two+r3(17)*zz*two+r3(29)*yz+r3(37)*three+r2( 6)*yzz &
&                   +r2( 2)*zzz+r2(12)*yzz*two+r2(14)*xyz(3)*two+r2(24)*xyz(2)*three &
&                   +r2(20)*xyz(3)+r1( 2)*zz*yz+r1( 8)*yz*two+r1(11)*yz
        cint1(6,6,2)=r5( 9)+r4( 9)*xyz(3)*two+r4(24)*xyz(3)*two+r3( 9)*zz+r3(19)*zz*four &
&                   +r3(29)*zz+r3(39)*six+r2( 6)*zzz*two+r2(12)*zzz*two+r2(18)*xyz(3)*six &
&                   +r2(24)*xyz(3)*six+r1( 2)*zz*zz+r1( 5)*zz+r1( 8)*zz*four+r1(11)*zz &
&                   +r1(14)*three
        cint1(1,1,3)=r5( 5)+r4( 5)*xyz(1)*two+r4(20)*xyz(1)*two+r3( 5)*xx+r3(15)*xx*four &
&                   +r3(25)*xx+r3(35)*six+r2( 5)*xxx*two+r2(11)*xxx*two+r2(17)*xyz(1)*six &
&                   +r2(23)*xyz(1)*six+r1( 3)*xx*xx+r1( 6)*xx+r1( 9)*xx*four+r1(12)*xx &
&                   +r1(15)*three
        cint1(2,1,3)=r5(11)+r4(13)*xyz(1)*two+r4(28)*xyz(1)+r4(20)*xyz(2)+r3(10)*xx &
&                   +r3(20)*xx*two+r3(15)*xy*two+r3(25)*xy+r3(40)*three+r2( 6)*xxx &
&                   +r2( 5)*xxy+r2(11)*xxy*two+r2(18)*xyz(1)*two+r2(24)*xyz(1) &
&                   +r2(23)*xyz(2)*three+r1( 3)*xx*xy+r1( 9)*xy*two+r1(12)*xy
        cint1(3,1,3)=r5(12)+r4(11)*xyz(1)*two+r4(26)*xyz(1)+r4(20)*xyz(3)+r3( 8)*xx &
&                   +r3(18)*xx*two+r3(15)*xz*two+r3(25)*xz+r3(38)*three+r2( 3)*xxx &
&                   +r2( 5)*xxz+r2(11)*xxz*two+r2(15)*xyz(1)*two+r2(21)*xyz(1) &
&                   +r2(23)*xyz(3)*three+r1( 3)*xx*xz+r1( 9)*xz*two+r1(12)*xz
        cint1(4,1,3)=r5(19)+r4(14)*xyz(1)*two+r4(28)*xyz(2)*two+r3( 7)*xx+r3(20)*xy*four &
&                   +r3(25)*yy+r3(35)+r3(37)+r2( 6)*xxy*two+r2(11)*xyy*two+r2(17)*xyz(1)*two &
&                   +r2(24)*xyz(2)*two+r1( 3)*xx*yy+r1( 6)*xx+r1(12)*yy+r1(15)
        cint1(5,1,3)=r5(20)+r4(15)*xyz(1)*two+r4(26)*xyz(2)+r4(28)*xyz(3)+r3( 9)*xx &
&                   +r3(18)*xy*two+r3(20)*xz*two+r3(25)*yz+r3(39)+r2( 3)*xxy+r2( 6)*xxz &
&                   +r2(11)*xyz2*two+r2(21)*xyz(2)+r2(24)*xyz(3)+r1( 3)*xx*yz+r1(12)*yz
        cint1(6,1,3)=r5(16)+r4( 8)*xyz(1)*two+r4(26)*xyz(3)*two+r3( 3)*xx+r3(18)*xz*four &
&                   +r3(25)*zz+r3(33)+r3(35)+r2( 3)*xxz*two+r2(11)*xzz*two+r2(17)*xyz(1)*two &
&                   +r2(21)*xyz(3)*two+r1( 3)*xx*zz+r1( 6)*xx+r1(12)*zz+r1(15)
        cint1(1,2,3)=r5(11)+r4(13)*xyz(1)+r4( 5)*xyz(2)+r4(28)*xyz(1)*two+r3( 5)*xy &
&                   +r3(20)*xx*two+r3(15)*xy*two+r3(30)*xx+r3(40)*three+r2( 5)*xxy*two &
&                   +r2(12)*xxx+r2(11)*xxy+r2(18)*xyz(1)+r2(17)*xyz(2)*three &
&                   +r2(24)*xyz(1)*two+r1( 3)*xy*xx+r1( 6)*xy+r1( 9)*xy*two
        cint1(2,2,3)=r5(19)+r4(14)*xyz(1)+r4(13)*xyz(2)+r4(29)*xyz(1)+r4(28)*xyz(2) &
&                   +r3(10)*xy+r3(17)*xx+r3(20)*xy+r3(20)*xy+r3(15)*yy+r3(30)*xy+r3(35) &
&                   +r3(37)+r2( 6)*xxy+r2( 5)*xyy+r2(12)*xxy+r2(11)*xyy+r2(17)*xyz(1) &
&                   +r2(18)*xyz(2)+r2(23)*xyz(1)+r2(24)*xyz(2)+r1( 3)*xy*xy+r1( 9)*xx &
&                   +r1( 9)*yy+r1(15)
        cint1(3,2,3)=r5(20)+r4(15)*xyz(1)+r4(11)*xyz(2)+r4(30)*xyz(1)+r4(28)*xyz(3) &
&                   +r3( 8)*xy+r3(19)*xx+r3(20)*xz+r3(18)*xy+r3(15)*yz+r3(30)*xz+r3(39) &
&                   +r2( 3)*xxy+r2( 5)*xyz2+r2(12)*xxz+r2(11)*xyz2+r2(15)*xyz(2) &
&                   +r2(24)*xyz(3)+r1( 3)*xy*xz+r1( 9)*yz
        cint1(4,2,3)=r5(14)+r4( 7)*xyz(1)+r4(14)*xyz(2)+r4(29)*xyz(2)*two+r3( 7)*xy &
&                   +r3(17)*xy*two+r3(20)*yy*two+r3(30)*yy+r3(40)*three+r2( 6)*xyy*two &
&                   +r2(12)*xyy+r2(11)*yyy+r2(18)*xyz(1)*three+r2(17)*xyz(2) &
&                   +r2(23)*xyz(2)*two+r1( 3)*xy*yy+r1( 6)*xy+r1( 9)*xy*two
        cint1(5,2,3)=r5(21)+r4(12)*xyz(1)+r4(15)*xyz(2)+r4(30)*xyz(2)+r4(29)*xyz(3) &
&                   +r3( 9)*xy+r3(19)*xy+r3(17)*xz+r3(18)*yy+r3(20)*yz+r3(30)*yz+r3(38) &
&                   +r2( 3)*xyy+r2( 6)*xyz2+r2(12)*xyz2+r2(11)*yyz+r2(15)*xyz(1) &
&                   +r2(23)*xyz(3)+r1( 3)*xy*yz+r1( 9)*xz
        cint1(6,2,3)=r5(17)+r4( 9)*xyz(1)+r4( 8)*xyz(2)+r4(30)*xyz(3)*two+r3( 3)*xy &
&                   +r3(19)*xz*two+r3(18)*yz*two+r3(30)*zz+r3(40)+r2( 3)*xyz2*two+r2(12)*xzz &
&                   +r2(11)*yzz+r2(18)*xyz(1)+r2(17)*xyz(2)+r1( 3)*xy*zz+r1( 6)*xy
        cint1(1,3,3)=r5(12)+r4(11)*xyz(1)+r4( 5)*xyz(3)+r4(26)*xyz(1)*two+r3( 5)*xz &
&                   +r3(18)*xx*two+r3(15)*xz*two+r3(28)*xx+r3(38)*three+r2( 5)*xxz*two &
&                   +r2( 9)*xxx+r2(11)*xxz+r2(15)*xyz(1)+r2(17)*xyz(3)*three &
&                   +r2(21)*xyz(1)*two+r1( 3)*xz*xx+r1( 6)*xz+r1( 9)*xz*two
        cint1(2,3,3)=r5(20)+r4(15)*xyz(1)+r4(13)*xyz(3)+r4(30)*xyz(1)+r4(26)*xyz(2) &
&                   +r3(10)*xz+r3(19)*xx+r3(18)*xy+r3(20)*xz+r3(15)*yz+r3(28)*xy+r3(39) &
&                   +r2( 6)*xxz+r2( 5)*xyz2+r2( 9)*xxy+r2(11)*xyz2+r2(18)*xyz(3) &
&                   +r2(21)*xyz(2)+r1( 3)*xz*xy+r1( 9)*yz
        cint1(3,3,3)=r5(16)+r4( 8)*xyz(1)+r4(11)*xyz(3)+r4(23)*xyz(1)+r4(26)*xyz(3) &
&                   +r3( 8)*xz+r3(13)*xx+r3(18)*xz+r3(18)*xz+r3(15)*zz+r3(28)*xz+r3(33) &
&                   +r3(35)+r2( 3)*xxz+r2( 5)*xzz+r2( 9)*xxz+r2(11)*xzz+r2(17)*xyz(1) &
&                   +r2(15)*xyz(3)+r2(23)*xyz(1)+r2(21)*xyz(3)+r1( 3)*xz*xz+r1( 9)*xx &
&                   +r1( 9)*zz+r1(15)
        cint1(4,3,3)=r5(21)+r4(12)*xyz(1)+r4(14)*xyz(3)+r4(30)*xyz(2)*two+r3( 7)*xz &
&                   +r3(19)*xy*two+r3(20)*yz*two+r3(28)*yy+r3(38)+r2( 6)*xyz2*two+r2( 9)*xyy &
&                   +r2(11)*yyz+r2(15)*xyz(1)+r2(17)*xyz(3)+r1( 3)*xz*yy+r1( 6)*xz
        cint1(5,3,3)=r5(17)+r4( 9)*xyz(1)+r4(15)*xyz(3)+r4(23)*xyz(2)+r4(30)*xyz(3) &
&                   +r3( 9)*xz+r3(13)*xy+r3(19)*xz+r3(18)*yz+r3(20)*zz+r3(28)*yz+r3(40) &
&                   +r2( 3)*xyz2+r2( 6)*xzz+r2( 9)*xyz2+r2(11)*yzz+r2(18)*xyz(1)+r2(23)*xyz(2) &
&                   +r1( 3)*xz*yz+r1( 9)*xy
        cint1(6,3,3)=r5( 8)+r4( 3)*xyz(1)+r4( 8)*xyz(3)+r4(23)*xyz(3)*two+r3( 3)*xz &
&                   +r3(13)*xz*two+r3(18)*zz*two+r3(28)*zz+r3(38)*three+r2( 3)*xzz*two &
&                   +r2( 9)*xzz+r2(11)*zzz+r2(15)*xyz(1)*three+r2(17)*xyz(3) &
&                   +r2(23)*xyz(3)*two+r1( 3)*xz*zz+r1( 6)*xz+r1( 9)*xz*two
        cint1(1,4,3)=r5(19)+r4(13)*xyz(2)*two+r4(29)*xyz(1)*two+r3( 5)*yy+r3(20)*xy*four &
&                   +r3(27)*xx+r3(35)+r3(37)+r2( 5)*xyy*two+r2(12)*xxy*two+r2(18)*xyz(2)*two &
&                   +r2(23)*xyz(1)*two+r1( 3)*yy*xx+r1( 6)*yy+r1(12)*xx+r1(15)
        cint1(2,4,3)=r5(14)+r4(14)*xyz(2)*two+r4(22)*xyz(1)+r4(29)*xyz(2)+r3(10)*yy &
&                   +r3(17)*xy*two+r3(20)*yy*two+r3(27)*xy+r3(40)*three+r2( 6)*xyy &
&                   +r2( 5)*yyy+r2(12)*xyy*two+r2(17)*xyz(2)*two+r2(24)*xyz(1)*three &
&                   +r2(23)*xyz(2)+r1( 3)*yy*xy+r1( 9)*xy*two+r1(12)*xy
        cint1(3,4,3)=r5(21)+r4(15)*xyz(2)*two+r4(27)*xyz(1)+r4(29)*xyz(3)+r3( 8)*yy &
&                   +r3(19)*xy*two+r3(20)*yz*two+r3(27)*xz+r3(38)+r2( 3)*xyy+r2( 5)*yyz &
&                   +r2(12)*xyz2*two+r2(21)*xyz(1)+r2(23)*xyz(3)+r1( 3)*yy*xz+r1(12)*xz
        cint1(4,4,3)=r5( 7)+r4( 7)*xyz(2)*two+r4(22)*xyz(2)*two+r3( 7)*yy+r3(17)*yy*four &
&                   +r3(27)*yy+r3(37)*six+r2( 6)*yyy*two+r2(12)*yyy*two+r2(18)*xyz(2)*six &
&                   +r2(24)*xyz(2)*six+r1( 3)*yy*yy+r1( 6)*yy+r1( 9)*yy*four+r1(12)*yy &
&                   +r1(15)*three
        cint1(5,4,3)=r5(15)+r4(12)*xyz(2)*two+r4(27)*xyz(2)+r4(22)*xyz(3)+r3( 9)*yy &
&                   +r3(19)*yy*two+r3(17)*yz*two+r3(27)*yz+r3(39)*three+r2( 3)*yyy &
&                   +r2( 6)*yyz+r2(12)*yyz*two+r2(15)*xyz(2)*two+r2(21)*xyz(2) &
&                   +r2(24)*xyz(3)*three+r1( 3)*yy*yz+r1( 9)*yz*two+r1(12)*yz
        cint1(6,4,3)=r5(18)+r4( 9)*xyz(2)*two+r4(27)*xyz(3)*two+r3( 3)*yy+r3(19)*yz*four &
&                   +r3(27)*zz+r3(33)+r3(37)+r2( 3)*yyz*two+r2(12)*yzz*two+r2(18)*xyz(2)*two &
&                   +r2(21)*xyz(3)*two+r1( 3)*yy*zz+r1( 6)*yy+r1(12)*zz+r1(15)
        cint1(1,5,3)=r5(20)+r4(11)*xyz(2)+r4(13)*xyz(3)+r4(30)*xyz(1)*two+r3( 5)*yz &
&                   +r3(18)*xy*two+r3(20)*xz*two+r3(29)*xx+r3(39)+r2( 5)*xyz2*two+r2( 9)*xxy &
&                   +r2(12)*xxz+r2(15)*xyz(2)+r2(18)*xyz(3)+r1( 3)*yz*xx+r1( 6)*yz
        cint1(2,5,3)=r5(21)+r4(15)*xyz(2)+r4(14)*xyz(3)+r4(27)*xyz(1)+r4(30)*xyz(2) &
&                   +r3(10)*yz+r3(19)*xy+r3(18)*yy+r3(17)*xz+r3(20)*yz+r3(29)*xy+r3(38) &
&                   +r2( 6)*xyz2+r2( 5)*yyz+r2( 9)*xyy+r2(12)*xyz2+r2(17)*xyz(3) &
&                   +r2(21)*xyz(1)+r1( 3)*yz*xy+r1( 9)*xz
        cint1(3,5,3)=r5(17)+r4( 8)*xyz(2)+r4(15)*xyz(3)+r4(24)*xyz(1)+r4(30)*xyz(3) &
&                   +r3( 8)*yz+r3(13)*xy+r3(18)*yz+r3(19)*xz+r3(20)*zz+r3(29)*xz+r3(40) &
&                   +r2( 3)*xyz2+r2( 5)*yzz+r2( 9)*xyz2+r2(12)*xzz+r2(17)*xyz(2) &
&                   +r2(24)*xyz(1)+r1( 3)*yz*xz+r1( 9)*xy
        cint1(4,5,3)=r5(15)+r4(12)*xyz(2)+r4( 7)*xyz(3)+r4(27)*xyz(2)*two+r3( 7)*yz &
&                   +r3(19)*yy*two+r3(17)*yz*two+r3(29)*yy+r3(39)*three+r2( 6)*yyz*two &
&                   +r2( 9)*yyy+r2(12)*yyz+r2(15)*xyz(2)+r2(18)*xyz(3)*three &
&                   +r2(21)*xyz(2)*two+r1( 3)*yz*yy+r1( 6)*yz+r1( 9)*yz*two
        cint1(5,5,3)=r5(18)+r4( 9)*xyz(2)+r4(12)*xyz(3)+r4(24)*xyz(2)+r4(27)*xyz(3) &
&                   +r3( 9)*yz+r3(13)*yy+r3(19)*yz+r3(19)*yz+r3(17)*zz+r3(29)*yz+r3(33) &
&                   +r3(37)+r2( 3)*yyz+r2( 6)*yzz+r2( 9)*yyz+r2(12)*yzz+r2(18)*xyz(2) &
&                   +r2(15)*xyz(3)+r2(24)*xyz(2)+r2(21)*xyz(3)+r1( 3)*yz*yz+r1( 9)*yy &
&                   +r1( 9)*zz+r1(15)
        cint1(6,5,3)=r5( 9)+r4( 3)*xyz(2)+r4( 9)*xyz(3)+r4(24)*xyz(3)*two+r3( 3)*yz &
&                   +r3(13)*yz*two+r3(19)*zz*two+r3(29)*zz+r3(39)*three+r2( 3)*yzz*two &
&                   +r2( 9)*yzz+r2(12)*zzz+r2(15)*xyz(2)*three+r2(18)*xyz(3) &
&                   +r2(24)*xyz(3)*two+r1( 3)*yz*zz+r1( 6)*yz+r1( 9)*yz*two
        cint1(1,6,3)=r5(16)+r4(11)*xyz(3)*two+r4(23)*xyz(1)*two+r3( 5)*zz+r3(18)*xz*four &
&                   +r3(23)*xx+r3(33)+r3(35)+r2( 5)*xzz*two+r2( 9)*xxz*two+r2(15)*xyz(3)*two &
&                   +r2(23)*xyz(1)*two+r1( 3)*zz*xx+r1( 6)*zz+r1(12)*xx+r1(15)
        cint1(2,6,3)=r5(17)+r4(15)*xyz(3)*two+r4(24)*xyz(1)+r4(23)*xyz(2)+r3(10)*zz &
&                   +r3(19)*xz*two+r3(18)*yz*two+r3(23)*xy+r3(40)+r2( 6)*xzz+r2( 5)*yzz &
&                   +r2( 9)*xyz2*two+r2(24)*xyz(1)+r2(23)*xyz(2)+r1( 3)*zz*xy+r1(12)*xy
        cint1(3,6,3)=r5( 8)+r4( 8)*xyz(3)*two+r4(18)*xyz(1)+r4(23)*xyz(3)+r3( 8)*zz &
&                   +r3(13)*xz*two+r3(18)*zz*two+r3(23)*xz+r3(38)*three+r2( 3)*xzz &
&                   +r2( 5)*zzz+r2( 9)*xzz*two+r2(17)*xyz(3)*two+r2(21)*xyz(1)*three &
&                   +r2(23)*xyz(3)+r1( 3)*zz*xz+r1( 9)*xz*two+r1(12)*xz
        cint1(4,6,3)=r5(18)+r4(12)*xyz(3)*two+r4(24)*xyz(2)*two+r3( 7)*zz+r3(19)*yz*four &
&                   +r3(23)*yy+r3(33)+r3(37)+r2( 6)*yzz*two+r2( 9)*yyz*two+r2(15)*xyz(3)*two &
&                   +r2(24)*xyz(2)*two+r1( 3)*zz*yy+r1( 6)*zz+r1(12)*yy+r1(15)
        cint1(5,6,3)=r5( 9)+r4( 9)*xyz(3)*two+r4(18)*xyz(2)+r4(24)*xyz(3)+r3( 9)*zz &
&                   +r3(13)*yz*two+r3(19)*zz*two+r3(23)*yz+r3(39)*three+r2( 3)*yzz &
&                   +r2( 6)*zzz+r2( 9)*yzz*two+r2(18)*xyz(3)*two+r2(21)*xyz(2)*three &
&                   +r2(24)*xyz(3)+r1( 3)*zz*yz+r1( 9)*yz*two+r1(12)*yz
        cint1(6,6,3)=r5( 3)+r4( 3)*xyz(3)*two+r4(18)*xyz(3)*two+r3( 3)*zz+r3(13)*zz*four &
&                   +r3(23)*zz+r3(33)*six+r2( 3)*zzz*two+r2( 9)*zzz*two+r2(15)*xyz(3)*six &
&                   +r2(21)*xyz(3)*six+r1( 3)*zz*zz+r1( 6)*zz+r1( 9)*zz*four+r1(12)*zz &
&                   +r1(15)*three
!
        if(nbfij(1) == 6) then
          do k= 1,3
            do j= 1,6
              cint1(j,2,k)= cint1(j,2,k)*sqrt3
              cint1(j,3,k)= cint1(j,3,k)*sqrt3
              cint1(j,5,k)= cint1(j,5,k)*sqrt3
            enddo
          enddo
        else
          do k= 1,3
            do j= 1,6
              do i= 1,6
                work(i)= cint1(j,i,k)
              enddo
              cint1(j,1,k)= work(2)*sqrt3
              cint1(j,2,k)= work(5)*sqrt3
              cint1(j,3,k)=(work(6)*two-work(1)-work(4))*half
              cint1(j,4,k)= work(3)*sqrt3
              cint1(j,5,k)=(work(1)-work(4))*sqrt3h
            enddo
          enddo
        endif
        if(nbfij(2) == 6) then
          do k= 1,3
            do i= 1,nbfij(1)
              cint1(2,i,k)= cint1(2,i,k)*sqrt3
              cint1(3,i,k)= cint1(3,i,k)*sqrt3
              cint1(5,i,k)= cint1(5,i,k)*sqrt3
            enddo
          enddo
        else
          do k= 1,3
            do i= 1,nbfij(1)
              do j= 1,6
                work(j)= cint1(j,i,k)
              enddo
              cint1(1,i,k)= work(2)*sqrt3
              cint1(2,i,k)= work(5)*sqrt3
              cint1(3,i,k)=(work(6)*two-work(1)-work(4))*half
              cint1(4,i,k)= work(3)*sqrt3
              cint1(5,i,k)=(work(1)-work(4))*sqrt3h
            enddo
          enddo
        endif
!
        if(.not.iandj) then
          do k= 1,3
            do i= 1,nbfij(1)
              do j= 1,nbfij(2)
                egrad(k,iatom)= egrad(k,iatom)+fulldmtrx(locbfij(2)+j,locbfij(1)+i)*cint1(j,i,k)
              enddo
            enddo
          enddo
        else
          do k= 1,3
            do i= 1,nbfij(1)
              do j= 1,nbfij(2)
                egrad(k,iatom)= egrad(k,iatom) &
&                              +fulldmtrx(locbfij(2)+j,locbfij(1)+i)*cint1(j,i,k)*half
              enddo
            enddo
          enddo
        endif
      enddo
!
      return
end


!--------------------------------------------------------------------------------
  subroutine int1grys(egrad,fulldmtrx,exij,coij,coordij,coord,znuc,natom,nao, &
&                     nprimij,nangij,nbfij,locbfij,mxprsh,threshex,iandj)
!--------------------------------------------------------------------------------
!
! Calculate derivative of 1-electron Coulomb integrals (j|Z/r|i) using Rys quadratures
!
      implicit none
      integer,intent(in) :: nprimij(2), nangij(2), nbfij(2), locbfij(2), natom, nao, mxprsh
      integer :: nroots, ncart(0:6), ncarti, ncartj, i, j, k, iprim, jprim, iatom, iroot
      integer :: iang, jang, ix, iy, iz, jx, jy, jz
      real(8),parameter :: sqrtpi2=1.128379167095513D+00 !2.0/sqrt(pi)
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, two=2.0D+00, three=3.0D+00
      real(8),parameter :: four=4.0D+00, six=6.0D+00, eight=8.0D+00, p24=24.0D+00, eighth=0.125D+00
      real(8),parameter :: sqrt3=1.732050807568877D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrt5=2.236067977499790D+00, sqrt15=3.872983346207417D+00
      real(8),parameter :: sqrt7=2.645751311064590D+00, sqrt35=5.916079783099616D+00
      real(8),parameter :: sqrt35third=3.415650255319866D+00
      real(8),parameter :: sqrt21=4.582575694955840D+00, sqrt63=7.937253933193772D+00
      real(8),parameter :: sqrt105=1.024695076595960D+01
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
      real(8),intent(in) :: fulldmtrx(nao,nao), exij(mxprsh,2), coij(mxprsh,2), coordij(3,2)
      real(8),intent(in) :: coord(3,natom), znuc(natom), threshex
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: xyz(3), rij, rij2, exi, exj, ci, cj, ex1, ex2, ex3, ex4, fac, rc, tval
      real(8) :: pixyz(3), pijxyz(3), pcxyz(3), xyzpijk(3,3), trys(13), wrys(13)
      real(8) :: cx(0:6,0:6,7,2), cy(0:6,0:6,7,2), cz(0:6,0:6,7,2), ww, xyzint(3)
      real(8) :: dcint1(28,28,3)
      logical,intent(in) :: iandj
      data ncart /1,3,6,10,15,21,28/
!
      nroots=(nangij(1)+nangij(2)+1)/2+1
      ncarti= ncart(nangij(1))
      ncartj= ncart(nangij(2))
!
      do i= 1,3
        xyz(i)= coordij(i,1)-coordij(i,2)
      enddo
      rij= xyz(1)*xyz(1)+xyz(2)*xyz(2)+xyz(3)*xyz(3)
!
      do iatom= 1,natom
        dcint1(1:ncartj,1:ncarti,1:3)= zero
        do iprim= 1,nprimij(1)
          exi= exij(iprim,1)
          ci = coij(iprim,1)*sqrtpi2
          do i= 1,3
            pixyz(i)= exi*coordij(i,1)
          enddo
          do jprim= 1,nprimij(2)
            exj= exij(jprim,2)
            ex1= exi+exj
            ex2= one/ex1
            ex3= exi*exj
            rij2=rij*ex2*ex3
            if(rij2 > threshex) cycle
            cj = coij(jprim,2)
            fac=exp(-rij2)*ex2*ci*cj
!
            do i= 1,3
              pijxyz(i)=(pixyz(i)+exj*coordij(i,2))*ex2
            enddo
            do i= 1,3
              pcxyz(i)= pijxyz(i)-coord(i,iatom)
            enddo
            rc= pcxyz(1)*pcxyz(1)+pcxyz(2)*pcxyz(2)+pcxyz(3)*pcxyz(3)
            tval= ex1*rc
            call rysquad(tval,trys,wrys,nroots)
            do iroot= 1,nroots
              ww=-wrys(iroot)*znuc(iatom)*two*ex1*trys(iroot)/(one-trys(iroot))
              ex4= sqrt(ex2*(one-trys(iroot)))
              do i= 1,3
                xyzpijk(i,1)=(one-trys(iroot))*pijxyz(i)+trys(iroot)*coord(i,iatom)-coordij(i,1)
                xyzpijk(i,2)=(one-trys(iroot))*pijxyz(i)+trys(iroot)*coord(i,iatom)-coordij(i,2)
                xyzpijk(i,3)=(one-trys(iroot))*pijxyz(i)+trys(iroot)*coord(i,iatom)-coord(i,iatom)
              enddo
              do iang= 0,nangij(1)
                do jang= 0,nangij(2)
                  call ghquad(xyzint,ex4,xyzpijk,iang,jang)
                  cx(jang,iang,iroot,1)= xyzint(1)
                  cy(jang,iang,iroot,1)= xyzint(2)
                  cz(jang,iang,iroot,1)= xyzint(3)*ww
                  call dghquad(xyzint,ex4,xyzpijk,iang,jang)
                  cx(jang,iang,iroot,2)= xyzint(1)
                  cy(jang,iang,iroot,2)= xyzint(2)
                  cz(jang,iang,iroot,2)= xyzint(3)*ww
                enddo
              enddo
            enddo
            i= 0
            do ix= nangij(1),0,-1
              do iy= nangij(1)-ix,0,-1
                iz= nangij(1)-ix-iy
                i= i+1
                j= 0
                do jx= nangij(2),0,-1
                  do jy= nangij(2)-jx,0,-1
                    jz= nangij(2)-jx-jy
                    j= j+1
                    do iroot= 1,nroots
                      dcint1(j,i,1)= dcint1(j,i,1) &
&                                   +fac*cx(jx,ix,iroot,2)*cy(jy,iy,iroot,1)*cz(jz,iz,iroot,1)
                      dcint1(j,i,2)= dcint1(j,i,2) &
&                                   +fac*cx(jx,ix,iroot,1)*cy(jy,iy,iroot,2)*cz(jz,iz,iroot,1)
                      dcint1(j,i,3)= dcint1(j,i,3) &
&                                   +fac*cx(jx,ix,iroot,1)*cy(jy,iy,iroot,1)*cz(jz,iz,iroot,2)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
!
        if((nbfij(1) >= 5).or.(nbfij(2) >= 5)) then
          call nrmlz1(dcint1(1,1,1),nbfij(1),nbfij(2),ncarti)
          call nrmlz1(dcint1(1,1,2),nbfij(1),nbfij(2),ncarti)
          call nrmlz1(dcint1(1,1,3),nbfij(1),nbfij(2),ncarti)
        endif
!
        if(.not.iandj) then
          do k= 1,3
            do i= 1,nbfij(1)
              do j= 1,nbfij(2)
                egrad(k,iatom)= egrad(k,iatom)+fulldmtrx(locbfij(2)+j,locbfij(1)+i)*dcint1(j,i,k)
              enddo
            enddo
          enddo
        else
          do k= 1,3
            do i= 1,nbfij(1)
              do j= 1,nbfij(2)
                egrad(k,iatom)= egrad(k,iatom) &
&                              +fulldmtrx(locbfij(2)+j,locbfij(1)+i)*dcint1(j,i,k)*half
              enddo
            enddo
          enddo
        endif
!
      enddo
!
      return
end


!-----------------------------------------------------
  subroutine dghquad(xyzint,expgh,xyzpijk,iang,jang)
!-----------------------------------------------------
!
! Calculate Gauss-Hermite quadrature
!
      implicit none
      integer,intent(in) :: iang, jang
      integer :: minh(12)= (/1,2,4, 7,11,16,22,29,37,46,56,67/)
      integer :: maxh(12)= (/1,3,6,10,15,21,28,36,45,55,66,78/)
      integer :: nroot, ij, i, j
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: expgh, xyzpijk(3,3)
      real(8),intent(out) :: xyzint(3)
      real(8) :: hnode(78), hweight(78)
      real(8) :: ghxyz(3), exnode, pxyz(3)
      data hnode/ &
        0.0000000000000000D+00,-0.7071067811865476D+00, 0.7071067811865476D+00, &
&      -0.1224744871391589D+01, 0.0000000000000000D+00, 0.1224744871391589D+01, &
&      -0.1650680123885785D+01,-0.5246476232752903D+00, 0.5246476232752903D+00, &
&       0.1650680123885785D+01,-0.2020182870456086D+01,-0.9585724646138185D+00, &
&       0.0000000000000000D+00, 0.9585724646138185D+00, 0.2020182870456086D+01, &
&      -0.2350604973674492D+01,-0.1335849074013697D+01,-0.4360774119276165D+00, &
&       0.4360774119276165D+00, 0.1335849074013697D+01, 0.2350604973674492D+01, &
&      -0.2651961356835233D+01,-0.1673551628767471D+01,-0.8162878828589647D+00, &
&       0.0000000000000000D+00, 0.8162878828589647D+00, 0.1673551628767471D+01, &
&       0.2651961356835233D+01,-0.2930637420257244D+01,-0.1981656756695843D+01, &
&      -0.1157193712446780D+01,-0.3811869902073221D+00, 0.3811869902073221D+00, &
&       0.1157193712446780D+01, 0.1981656756695843D+01, 0.2930637420257244D+01, &
&      -0.3190993201781528D+01,-0.2266580584531843D+01,-0.1468553289216668D+01, &
&      -0.7235510187528376D+00, 0.0000000000000000D+00, 0.7235510187528376D+00, &
&       0.1468553289216668D+01, 0.2266580584531843D+01, 0.3190993201781528D+01, &
&      -0.3436159118837737D+01,-0.2532731674232790D+01,-0.1756683649299882D+01, &
&      -0.1036610829789514D+01,-0.3429013272237046D+00, 0.3429013272237046D+00, &
&       0.1036610829789514D+01, 0.1756683649299882D+01, 0.2532731674232790D+01, &
&       0.3436159118837737D+01,-0.3668470846559583D+01,-0.2783290099781652D+01, &
&      -0.2025948015825755D+01,-0.1326557084494933D+01,-0.6568095668820998D+00, &
&       0.0000000000000000D+00, 0.6568095668820998D+00, 0.1326557084494933D+01, &
&       0.2025948015825755D+01, 0.2783290099781652D+01, 0.3668470846559583D+01, &
&      -0.3889724897869782D+01,-0.3020637025120890D+01,-0.2279507080501060D+01, &
&      -0.1597682635152605D+01,-0.9477883912401637D+00,-0.3142403762543591D+00, &
&       0.3142403762543591D+00, 0.9477883912401637D+00, 0.1597682635152605D+01, &
&       0.2279507080501060D+01, 0.3020637025120890D+01, 0.3889724897869782D+01/
      data hweight/ &
&       0.1772453850905516D+01, 0.8862269254527581D+00, 0.8862269254527581D+00, &
&       0.2954089751509194D+00, 0.1181635900603677D+01, 0.2954089751509194D+00, &
&       0.8131283544724517D-01, 0.8049140900055128D+00, 0.8049140900055128D+00, &
&       0.8131283544724517D-01, 0.1995324205904591D-01, 0.3936193231522412D+00, &
&       0.9453087204829419D+00, 0.3936193231522412D+00, 0.1995324205904591D-01, &
&       0.4530009905508846D-02, 0.1570673203228566D+00, 0.7246295952243925D+00, &
&       0.7246295952243925D+00, 0.1570673203228566D+00, 0.4530009905508846D-02, &
&       0.9717812450995191D-03, 0.5451558281912703D-01, 0.4256072526101278D+00, &
&       0.8102646175568073D+00, 0.4256072526101278D+00, 0.5451558281912703D-01, &
&       0.9717812450995191D-03, 0.1996040722113676D-03, 0.1707798300741347D-01, &
&       0.2078023258148919D+00, 0.6611470125582413D+00, 0.6611470125582413D+00, &
&       0.2078023258148919D+00, 0.1707798300741347D-01, 0.1996040722113676D-03, &
&       0.3960697726326439D-04, 0.4943624275536947D-02, 0.8847452739437657D-01, &
&       0.4326515590025558D+00, 0.7202352156060510D+00, 0.4326515590025558D+00, &
&       0.8847452739437657D-01, 0.4943624275536947D-02, 0.3960697726326439D-04, &
&       0.7640432855232621D-05, 0.1343645746781233D-02, 0.3387439445548106D-01, &
&       0.2401386110823147D+00, 0.6108626337353258D+00, 0.6108626337353258D+00, &
&       0.2401386110823147D+00, 0.3387439445548106D-01, 0.1343645746781233D-02, &
&       0.7640432855232621D-05, 0.1439560393714258D-05, 0.3468194663233455D-03, &
&       0.1191139544491153D-01, 0.1172278751677085D+00, 0.4293597523561250D+00, &
&       0.6547592869145918D+00, 0.4293597523561250D+00, 0.1172278751677085D+00, &
&       0.1191139544491153D-01, 0.3468194663233455D-03, 0.1439560393714258D-05, &
&       0.2658551684356306D-06, 0.8573687043587876D-04, 0.3905390584629067D-02, &
&       0.5160798561588394D-01, 0.2604923102641611D+00, 0.5701352362624796D+00, &
&       0.5701352362624796D+00, 0.2604923102641611D+00, 0.5160798561588394D-01, &
&       0.3905390584629067D-02, 0.8573687043587876D-04, 0.2658551684356306D-06/
!
      do i= 1,3
        xyzint(i)= zero
      enddo
!
      nroot=(iang+jang+1)/2+1
      do ij= minh(nroot),maxh(nroot)
        exnode= hnode(ij)*expgh
        do i= 1,3
          ghxyz(i)=(exnode+xyzpijk(i,3))*hweight(ij)
        enddo
        if(iang >= 1) then
          do i= 1,3
            pxyz(i)= exnode+xyzpijk(i,1)
          enddo
          do i= 1,iang
            do j= 1,3
              ghxyz(j)= ghxyz(j)*pxyz(j)
            enddo
          enddo
        endif
        if(jang >= 1) then
          do i= 1,3
            pxyz(i)= exnode+xyzpijk(i,2)
          enddo
          do i= 1,jang
            do j= 1,3
              ghxyz(j)= ghxyz(j)*pxyz(j)
            enddo
          enddo
        endif
        do i= 1,3
        xyzint(i)= xyzint(i)+ghxyz(i)
        enddo
      enddo
      return
end

