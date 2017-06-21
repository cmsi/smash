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
!--------------------------------------------------------------------------
  subroutine oneei(hstmat1,hstmat2,hstmat3,hstmat4,nproc,myrank,mpi_comm)
!--------------------------------------------------------------------------
!
! Driver of one-electron and overlap integrals
!
! Out : hstmat1 (One electron Hamiltonian matrix)
!       hstmat2 (Overlap integral matrix)
!       hstmat3 (Kinetic energy matrix)
!       hstmat4 (Work array)
!
      use modbasis, only : nao, nshell, mtype
      use modecp, only : flagecp
      implicit none
      integer,intent(in) :: nproc, myrank, mpi_comm
      integer :: ish, jsh, num, maxfunc(0:6), maxbasis, maxdim
      real(8),parameter :: zero=0.0D+00
      real(8),intent(out) :: hstmat1((nao*(nao+1))/2), hstmat2((nao*(nao+1))/2)
      real(8),intent(out) :: hstmat3((nao*(nao+1))/2), hstmat4((nao*(nao+1))/2)
      data maxfunc/1,3,6,10,15,21,28/
!
      maxbasis= maxval(mtype(1:nshell))
      maxdim= maxfunc(maxbasis)
!
      num=(nao*(nao+1))/2
      hstmat2(:)= zero
      hstmat3(:)= zero
      hstmat4(:)= zero
!
!$OMP parallel
      do ish= nshell-myrank,1,-nproc
!$OMP do
        do jsh= 1,ish
          call calcintst1c(hstmat2,hstmat3,hstmat4,ish,jsh,maxdim)
        enddo
!$OMP enddo
      enddo
!$OMP end parallel
!
! Calculate ECP integrals
!
      if(flagecp) call oneeiecp(hstmat2,nproc,myrank)
!
      call para_allreducer(hstmat2,hstmat1,num,mpi_comm)
      call para_allreducer(hstmat3,hstmat2,num,mpi_comm)
      call para_allreducer(hstmat4,hstmat3,num,mpi_comm)
!
      return
end


!------------------------------------------------------
  subroutine calcintst1c(hmat,smat,tmat,ish,jsh,len1)
!------------------------------------------------------
!
! Driver of overlap, kinetic, and 1-electron Coulomb integrals (j|Z/r|i)
!
      use modparam, only : mxprsh
      use modmolecule, only : natom, coord, znuc
      use modbasis, only : locatom, locprim, locbf, mprim, mbf, mtype, ex, coeff, nao
      use modthresh, only : threshex
      implicit none
      integer,intent(in) :: ish, jsh, len1
      integer :: nangij(2), nprimij(2), nbfij(2), iatom, jatom
      integer :: iloc, jloc, ilocbf, jlocbf, iprim, jprim, i, j, ii, ij, maxj
      real(8),intent(inout) :: hmat((nao*(nao+1))/2), smat((nao*(nao+1))/2)
      real(8),intent(inout) :: tmat((nao*(nao+1))/2)
      real(8) :: exij(mxprsh,2), coij(mxprsh,2), coordij(3,2)
      real(8) :: sint(len1,len1), tint(len1,len1), cint(len1,len1)
      logical :: iandj
!
      iandj=(ish == jsh)
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
      do i= 1,3
        coordij(i,1)= coord(i,iatom)
        coordij(i,2)= coord(i,jatom)
      enddo
      do iprim= 1,nprimij(1)
        exij(iprim,1)= ex(iloc+iprim)
        coij(iprim,1)= coeff(iloc+iprim)
      enddo
      do jprim= 1,nprimij(2)
        exij(jprim,2)= ex(jloc+jprim)
        coij(jprim,2)= coeff(jloc+jprim)
      enddo
!
      if((nangij(1) > 6).or.(nangij(2) > 6))then
        write(*,'(" Error! This program supports up to i function in calcintst1c.")')
        call exit
      endif
!
! Overlap and kinetic integrals
!
      call intst(sint,tint,exij,coij,coordij, &
&                nprimij,nangij,nbfij,len1,mxprsh,threshex)
!
! 1-electron Coulomb integrals
!
      if((nangij(1) <= 2).and.(nangij(2) <= 2)) then
        call int1cmd(cint,exij,coij,coordij,coord,znuc,natom, &
&                    nprimij,nangij,nbfij,len1,mxprsh,threshex)
      else
        call int1rys(cint,exij,coij,coordij,coord,znuc,natom, &
&                    nprimij,nangij,nbfij,len1,mxprsh,threshex)
      endif
!
      maxj= nbfij(2)
      do i= 1,nbfij(1)
        if(iandj) maxj= i
        ii= ilocbf+i
        ij= ii*(ii-1)/2+jlocbf
        do j= 1,maxj
          hmat(ij+j)= tint(j,i)+cint(j,i)
          smat(ij+j)= sint(j,i)
          tmat(ij+j)= tint(j,i)
        enddo
      enddo
!
      return
end


!--------------------------------------------------------------
  subroutine intst(sint,tint,exij,coij,coordij, &
&                  nprimij,nangij,nbfij,len1,mxprsh,threshex)
!--------------------------------------------------------------
!
! Calculate overlap and kinetic integrals
!
! In : exij     (Exponents of basis functions)
!      coij     (Coefficients of basis functions)
!      coordij  (x,y,z coordinates of basis functions)
!      nprimij  (Numbers of primitive functions)
!      nangij   (Degrees of angular momentum)
!      nbfij    (Numbers of basis functions)
!      len1     (Dimension of one-electron integral array)
!      mxprsh   (Size of primitive fuction array)
!      threshex (Threshold of exponential calculation)
! Out: sint     (Overlap integrals)
!      tint     (Kinetic integrals)
!
      implicit none
      integer,intent(in) :: nprimij(2), nangij(2), nbfij(2), len1, mxprsh
      integer :: ncarti, ncartj, iprim, jprim, ii, jj, iang, jang
      integer :: ncart(0:6), ix, iy, iz, jx, jy, jz
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00, half=0.5D+00
      real(8),intent(in) :: exij(mxprsh,2), coij(mxprsh,2), coordij(3,2), threshex
      real(8),intent(out) :: sint(len1,len1), tint(len1,len1)
      real(8) :: xyzij(3), rij, rij2, fac, exi, exi2, exj, ci, cj
      real(8) :: ex1, ex2, ex3, xyzpij(3,2), cij, sxyz, sxyz1, sxyz2, sxyz3
      real(8) :: xyzint(3), sx(0:6,0:8,2), sy(0:6,0:8,2), sz(0:6,0:8,2)
      real(8) :: sintmp(28,28), tintmp(28,28)
      data ncart /1,3,6,10,15,21,28/
!
      ncarti= ncart(nangij(1))
      ncartj= ncart(nangij(2))
      sintmp(1:ncartj,1:ncarti)= zero
      tintmp(1:ncartj,1:ncarti)= zero
!
      if((nangij(1) > 6).or.(nangij(2) > 6))then
        write(*,'(" Error! This program supports up to i function in intst.")')
        call abort
      endif
!
      do ii= 1,3
        xyzij(ii)= coordij(ii,1)-coordij(ii,2)
      enddo
      rij= xyzij(1)*xyzij(1)+xyzij(2)*xyzij(2)+xyzij(3)*xyzij(3)
!
! Calculate overlap and kinetic integrals for each primitive
!
      do iprim= 1,nprimij(1)
        exi= exij(iprim,1)
        ci = coij(iprim,1)
        exi2=-two*exi*exi
        do jprim= 1,nprimij(2)
          exj= exij(jprim,2)
          ex1= exi+exj
          ex2= one/ex1
          rij2=rij*exi*exj*ex2
          if(rij2 > threshex) cycle
          ex3= sqrt(ex2)
          fac= exp(-rij2)
          do ii= 1,3
            xyzpij(ii,1)=-exj*xyzij(ii)*ex2
            xyzpij(ii,2)= exi*xyzij(ii)*ex2
          enddo
          cj = coij(jprim,2)*fac
!
          do iang= 0,nangij(1)+2
            do jang= 0,nangij(2)
              call ghquad(xyzint,ex3,xyzpij,iang,jang)
              sx(jang,iang,1)= xyzint(1)*ex3
              sy(jang,iang,1)= xyzint(2)*ex3
              sz(jang,iang,1)= xyzint(3)*ex3
            enddo
          enddo
          do iang= 0,nangij(1)
            do jang= 0,nangij(2)
              sx(jang,iang,2)= sx(jang,iang+2,1)*exi2
              sy(jang,iang,2)= sy(jang,iang+2,1)*exi2
              sz(jang,iang,2)= sz(jang,iang+2,1)*exi2
              if(iang >= 2) then
                sx(jang,iang,2)=sx(jang,iang,2)-sx(jang,iang-2,1)*half*iang*(iang-1)
                sy(jang,iang,2)=sy(jang,iang,2)-sy(jang,iang-2,1)*half*iang*(iang-1)
                sz(jang,iang,2)=sz(jang,iang,2)-sz(jang,iang-2,1)*half*iang*(iang-1)
              endif
            enddo
          enddo
          cij= ci*cj
          ii= 0
          do ix= nangij(1),0,-1
            do iy= nangij(1)-ix,0,-1
              iz= nangij(1)-ix-iy
              ii= ii+1
              jj= 0
              do jx= nangij(2),0,-1
                do jy= nangij(2)-jx,0,-1
                  jz= nangij(2)-jx-jy
                  jj= jj+1
                  sxyz = sx(jx,ix,1)*sy(jy,iy,1)*sz(jz,iz,1)
                  sxyz1= sx(jx,ix,2)*sy(jy,iy,1)*sz(jz,iz,1)
                  sxyz2= sy(jy,iy,2)*sx(jx,ix,1)*sz(jz,iz,1)
                  sxyz3= sz(jz,iz,2)*sx(jx,ix,1)*sy(jy,iy,1)
                  sintmp(jj,ii)= sintmp(jj,ii)+cij*sxyz
                  tintmp(jj,ii)= tintmp(jj,ii)+cij*(sxyz1+sxyz2+sxyz3+sxyz*exi*(2*(ix+iy+iz)+3))
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
!
      if((nbfij(1) >= 5).or.(nbfij(2) >= 5)) then
        call nrmlz1(sintmp,nbfij(1),nbfij(2),ncarti)
        call nrmlz1(tintmp,nbfij(1),nbfij(2),ncarti)
      endif
      do ii= 1,nbfij(1)
        do jj= 1,nbfij(2)
          sint(jj,ii)= sintmp(jj,ii)
          tint(jj,ii)= tintmp(jj,ii)
        enddo
      enddo
!
      return
end


!---------------------------------------------------
  subroutine ghquad(xyzint,expgh,xyzpij,iang,jang)
!---------------------------------------------------
!
! Calculate Gauss-Hermite quadrature
!
      implicit none
      integer,intent(in) :: iang, jang
      integer :: minh(12)= (/1,2,4, 7,11,16,22,29,37,46,56,67/)
      integer :: maxh(12)= (/1,3,6,10,15,21,28,36,45,55,66,78/)
      integer :: nroot, ij, i, j
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: expgh, xyzpij(3,2)
      real(8),intent(out) :: xyzint(3)
      real(8) :: hnode(78), hweight(78)
      real(8) :: ghxyz(3), exnode, pxyz(3)
      data hnode/ &
&       0.0000000000000000D+00,-0.7071067811865476D+00, 0.7071067811865476D+00, &
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
      nroot=(iang+jang)/2+1
      do ij= minh(nroot),maxh(nroot)
        do i= 1,3
          ghxyz(i)= hweight(ij)
        enddo
        exnode= hnode(ij)*expgh
        if(iang >= 1) then
          do i= 1,3
            pxyz(i)= exnode+xyzpij(i,1)
          enddo
          do i= 1,iang
            do j= 1,3
              ghxyz(j)= ghxyz(j)*pxyz(j)
            enddo
          enddo
        endif
        if(jang >= 1) then
          do i= 1,3
            pxyz(i)= exnode+xyzpij(i,2)
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


!----------------------------------------------------------
  subroutine ghquadd(xyzint,expgh,xyzpijk,iang,jang,kang)
!----------------------------------------------------------
!
! Calculate Gauss-Hermite quadrature for dipole moment
!
      implicit none
      integer,intent(in) :: iang, jang, kang
      integer :: minh(12)= (/1,2,4, 7,11,16,22,29,37,46,56,67/)
      integer :: maxh(12)= (/1,3,6,10,15,21,28,36,45,55,66,78/)
      integer :: nroot, ijk, i, j
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: expgh, xyzpijk(3,3)
      real(8),intent(out) :: xyzint(3)
      real(8) :: hnode(78), hweight(78)
      real(8) :: ghxyz(3), exnode, pxyz(3)
      data hnode/ &
&       0.0000000000000000D+00,-0.7071067811865476D+00, 0.7071067811865476D+00, &
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
      nroot=(iang+jang+kang)/2+1
      do ijk= minh(nroot),maxh(nroot)
        do i= 1,3
          ghxyz(i)= hweight(ijk)
        enddo
        exnode= hnode(ijk)*expgh
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
        if(kang >= 1) then
          do i= 1,3
            pxyz(i)= exnode+xyzpijk(i,3)
          enddo
          do i= 1,kang
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



!------------------------------------------
  subroutine calcint1c(hmat,ish,jsh,len1)
!------------------------------------------
!
! Driver of 1-electron Coulomb integrals
!   (j|Z/r|i)
!
      use modparam, only : mxprsh
      use modmolecule, only : natom, coord, znuc
      use modbasis, only : locatom, locprim, locbf, mprim, mbf, mtype, ex, coeff, nao
      use modthresh, only : threshex
      implicit none
      integer,intent(in) :: ish, jsh, len1
      integer :: nangij(2), nprimij(2), nbfij(2), iatom, jatom
      integer :: iloc, jloc, ilocbf, jlocbf, iprim, jprim, i, j, ii, ij, maxj
      real(8),intent(inout) :: hmat((nao*(nao+1))/2)
      real(8) :: exij(mxprsh,2), coij(mxprsh,2), coordij(3,2), cint(len1,len1)
      logical :: iandj
!
      iandj=(ish == jsh)
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
      do i= 1,3
        coordij(i,1)= coord(i,iatom)
        coordij(i,2)= coord(i,jatom)
      enddo
      do iprim= 1,nprimij(1)
        exij(iprim,1)= ex(iloc+iprim)
        coij(iprim,1)= coeff(iloc+iprim)
      enddo
      do jprim= 1,nprimij(2)
        exij(jprim,2)= ex(jloc+jprim)
        coij(jprim,2)= coeff(jloc+jprim)
      enddo
!
      if((nangij(1) <= 2).and.(nangij(2) <= 2)) then
        call int1cmd(cint,exij,coij,coordij,coord,znuc,natom, &
&                    nprimij,nangij,nbfij,len1,mxprsh,threshex)
      else
        if((nangij(1) > 6).or.(nangij(2) > 6))then
          write(*,'(" Error! This program supports up to i function in calcint1c.")')
          call exit
        endif
!
        call int1rys(cint,exij,coij,coordij,coord,znuc,natom, &
&                    nprimij,nangij,nbfij,len1,mxprsh,threshex)
!
      endif
!
      maxj= nbfij(2)
      do i= 1,nbfij(1)
        if(iandj) maxj= i
        ii= ilocbf+i
        ij= ii*(ii-1)/2+jlocbf
        do j= 1,maxj
          hmat(ij+j)= hmat(ij+j)+cint(j,i)
        enddo
      enddo
!
      return
end


!-----------------------------------------------------------------
  subroutine int1cmd(cint,exij,coij,coordij,coord,znuc,natom, &
&                    nprimij,nangij,nbfij,len1,mxprsh,threshex)
!-----------------------------------------------------------------
!
! Driver of 1-electron Coulomb integrals (j|Z/r|i) using McMurchie-Davidson method
!
! In : exij     (Exponents of basis functions)
!      coij     (Coefficients of basis functions)
!      coordij  (x,y,z coordinates of basis functions)
!      coord    (x,y,z coordinates of atoms)
!      znuc     (Charges of atoms)
!      natom    (Number of atoms)
!      nprimij  (Numbers of primitive functions)
!      nangij   (Degrees of angular momentum)
!      nbfij    (Numbers of basis functions)
!      len1     (Dimension of one-electron integral array)
!      mxprsh   (Size of primitive fuction array)
!      threshex (Threshold of exponential calculation)
! Out: cint     (One-electron Coulomb integrals)
!
      implicit none
      integer,intent(in) :: nprimij(2), nangij(2), nbfij(2), len1, natom, mxprsh
      integer,parameter :: mxprsh2=30
      integer :: inttyp, nij, iprim, jprim, i, j, ii, jj, nbfij2(2)
      real(8),parameter :: one=1.0D+00, pi2=6.283185307179586D+00
      real(8),intent(in) :: exij(mxprsh,2), coij(mxprsh,2), coordij(3,2)
      real(8),intent(in) :: coord(3,natom), znuc(natom), threshex
      real(8),intent(out) :: cint(len1,len1)
      real(8) :: xyz(3), rij, exi, exj, ci, cj, ex12, ex21, ex2i, ex2j, rij2, pixyz(3)
      real(8) :: exfac(5,mxprsh2*mxprsh2), pijxyz(3,mxprsh2*mxprsh2), cint1(6,6)
!
      if(mxprsh > mxprsh2) then
        write(*,'(" Error! Parameter mxprsh2 in int1cmd is small!")')
        call abort
      endif
!
      inttyp=nangij(2)*3+nangij(1)+1
      if(nangij(2).ge.nangij(1)) then
        ii= 1
        jj= 2
        nbfij2(1)= nbfij(1)
        nbfij2(2)= nbfij(2)
      else
        ii= 2
        jj= 1
        nbfij2(1)= nbfij(2)
        nbfij2(2)= nbfij(1)
      endif
! 
      do i= 1,3
        xyz(i)= coordij(i,ii)-coordij(i,jj)
      enddo
      rij= xyz(1)*xyz(1)+xyz(2)*xyz(2)+xyz(3)*xyz(3)
      nij= 0
      do iprim= 1,nprimij(ii)
        exi = exij(iprim,ii)
        ci  = coij(iprim,ii)*pi2
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
          cj = coij(jprim,jj)
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
          call int1css(cint1,exfac,pijxyz,nij,coord,znuc,natom,mxprsh)
        case (2,4)
          call int1cps(cint1,exfac,pijxyz,xyz,nij,coord,znuc,natom,mxprsh)
        case (5)
          call int1cpp(cint1,exfac,pijxyz,xyz,nij,coord,znuc,natom,mxprsh)
        case (3,7)
          call int1cds(cint1,exfac,pijxyz,xyz,nij,coord,znuc,natom,mxprsh,nbfij2)
        case (6,8)
          call int1cdp(cint1,exfac,pijxyz,xyz,nij,coord,znuc,natom,mxprsh,nbfij2)
        case (9)
          call int1cdd(cint1,exfac,pijxyz,xyz,nij,coord,znuc,natom,mxprsh,nbfij2)
      end select
!
      if(ii == 1) then
        do i= 1,nbfij(1)
          do j= 1,nbfij(2)
            cint(j,i)= cint1(j,i)
          enddo
        enddo
      else
        do i= 1,nbfij(2)
          do j= 1,nbfij(1)
            cint(i,j)= cint1(j,i)
          enddo
        enddo
      endif
!
      return
end


!----------------------------------------------------------------
  subroutine int1rys(cint,exij,coij,coordij,coord,znuc,natom, &
&                    nprimij,nangij,nbfij,len1,mxprsh,threshex)
!----------------------------------------------------------------
!
! Calculate 1-electron Coulomb integrals (j|Z/r|i) using Rys quadratures
!
      implicit none
      integer,intent(in) :: nprimij(2), nangij(2), nbfij(2), len1, natom, mxprsh
      integer :: nroots, ncart(0:6), ncarti, ncartj, ii, jj, iprim, jprim, iatom, iroot
      integer :: iang, jang, ix, iy, iz, jx, jy, jz
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00
      real(8),parameter :: sqrtpi2=1.128379167095513D+00 !2.0/sqrt(pi)
      real(8),intent(in) :: exij(mxprsh,2), coij(mxprsh,2), coordij(3,2)
      real(8),intent(in) :: coord(3,natom), znuc(natom), threshex
      real(8),intent(out) :: cint(len1,len1)
      real(8) :: xyz(3), rij, rij2, exi, exj, ci, cj, ex1, ex2, ex3, ex4, fac, rc, tval
      real(8) :: pixyz(3), pijxyz(3), pcxyz(3), xyzpij(3,2), trys(13), wrys(13)
      real(8) :: cx(0:6,0:6,7), cy(0:6,0:6,7), cz(0:6,0:6,7), ww, xyzint(3)
      real(8) :: cintmp(28,28)
      data ncart /1,3,6,10,15,21,28/
!
      nroots=(nangij(1)+nangij(2))/2+1
      ncarti= ncart(nangij(1))
      ncartj= ncart(nangij(2))
!
      cintmp(1:ncartj,1:ncarti)= zero
!
      do ii= 1,3
        xyz(ii)= coordij(ii,1)-coordij(ii,2)
      enddo
      rij= xyz(1)*xyz(1)+xyz(2)*xyz(2)+xyz(3)*xyz(3)
!
      do iprim= 1,nprimij(1)
        exi= exij(iprim,1)
        ci = coij(iprim,1)*sqrtpi2
        do ii= 1,3
          pixyz(ii)= exi*coordij(ii,1)
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
          do ii= 1,3
            pijxyz(ii)=(pixyz(ii)+exj*coordij(ii,2))*ex2
          enddo
          do iatom= 1,natom
            do ii= 1,3
              pcxyz(ii)= pijxyz(ii)-coord(ii,iatom)
            enddo
            rc= pcxyz(1)*pcxyz(1)+pcxyz(2)*pcxyz(2)+pcxyz(3)*pcxyz(3)
            tval= ex1*rc
            call rysquad(tval,trys,wrys,nroots)
            do iroot= 1,nroots
              ww=-wrys(iroot)*znuc(iatom)
              ex4= sqrt(ex2*(one-trys(iroot)))
              do ii= 1,3
                xyzpij(ii,1)=(one-trys(iroot))*pijxyz(ii)+trys(iroot)*coord(ii,iatom)-coordij(ii,1)
                xyzpij(ii,2)=(one-trys(iroot))*pijxyz(ii)+trys(iroot)*coord(ii,iatom)-coordij(ii,2)
              enddo
              do iang= 0,nangij(1)
                do jang= 0,nangij(2)
                  call ghquad(xyzint,ex4,xyzpij,iang,jang)
                  cx(jang,iang,iroot)= xyzint(1)
                  cy(jang,iang,iroot)= xyzint(2)
                  cz(jang,iang,iroot)= xyzint(3)*ww
                enddo
              enddo
            enddo
            ii= 0
            do ix= nangij(1),0,-1
              do iy= nangij(1)-ix,0,-1
                iz= nangij(1)-ix-iy
                ii= ii+1
                jj= 0
                do jx= nangij(2),0,-1
                  do jy= nangij(2)-jx,0,-1
                    jz= nangij(2)-jx-jy
                    jj= jj+1
                    do iroot= 1,nroots
                      cintmp(jj,ii)= cintmp(jj,ii) &
&                                   +fac*cx(jx,ix,iroot)*cy(jy,iy,iroot)*cz(jz,iz,iroot)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
!
      if((nbfij(1) >= 5).or.(nbfij(2) >= 5)) then
        call nrmlz1(cintmp,nbfij(1),nbfij(2),ncarti)
      endif
      do ii= 1,nbfij(1)
        do jj= 1,nbfij(2)
          cint(jj,ii)= cintmp(jj,ii)
        enddo
      enddo
!
      return
end


!-----------------------------------------------------------------------
  subroutine int1css(cint1,exfac,pijxyz,nij,coord,znuc,natom,mxprsh)
!-----------------------------------------------------------------------
!
! Calculate 1-electron Coulomb integrals of (s|Z/r|s)
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nij, natom, mxprsh
      integer :: ijprim, i, iatom, igrid
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00
      real(8),parameter :: pi=3.141592653589793D+00
      real(8),intent(in) :: exfac(5,mxprsh*mxprsh), pijxyz(3,mxprsh*mxprsh)
      real(8),intent(in) :: coord(3,natom), znuc(natom)
      real(8),intent(out) :: cint1(6,6)
      real(8) :: rc, ex12, ex21, ex2i, ex2j, c12, pcxyz(3)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: ft0, fpc
!
      cint1(1,1)= zero
!
      do ijprim= 1,nij
        ex12= exfac(1,ijprim)
        ex21= exfac(2,ijprim)
        ex2i= exfac(3,ijprim)
        ex2j= exfac(4,ijprim)
        c12 = exfac(5,ijprim)
        fpc= zero
        do iatom= 1,natom
          do i= 1,3
            pcxyz(i)= pijxyz(i,ijprim)-coord(i,iatom)
          enddo
          rc= pcxyz(1)*pcxyz(1)+pcxyz(2)*pcxyz(2)+pcxyz(3)*pcxyz(3)
          tval= ex12*rc
          if(tval >= threshtval) then
            fpc= fpc+half*sqrt(pi/tval)*znuc(iatom)
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
            ft0= fgrid(0,0,igrid)      +fgrid( 1,0,igrid)*tval  +fgrid( 2,0,igrid)*tval2 &
&               +fgrid(3,0,igrid)*tval3+fgrid( 4,0,igrid)*tval4 +fgrid( 5,0,igrid)*tval5 &
&               +fgrid(6,0,igrid)*tval6+fgrid( 7,0,igrid)*tval7 +fgrid( 8,0,igrid)*tval8 &
&               +fgrid(9,0,igrid)*tval9+fgrid(10,0,igrid)*tval10
            fpc= fpc+ft0*znuc(iatom)
          endif
        enddo
        cint1(1,1)= cint1(1,1)-c12*fpc
      enddo
!
      return
end


!---------------------------------------------------------------------------
  subroutine int1cps(cint1,exfac,pijxyz,xyz,nij,coord,znuc,natom,mxprsh)
!---------------------------------------------------------------------------
!
! Calculate 1-electron Coulomb integrals of (p|Z/r|s)
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nij, natom, mxprsh
      integer :: ijprim, i, iatom, igrid
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, two=2.0D+00
      real(8),parameter :: pi=3.141592653589793D+00
      real(8),intent(in) :: exfac(5,mxprsh*mxprsh), pijxyz(3,mxprsh*mxprsh), xyz(3)
      real(8),intent(in) :: coord(3,natom), znuc(natom)
      real(8),intent(out) :: cint1(6,6)
      real(8) :: rc, ex12, ex21, ex2i, ex2j, c12, pcxyz(3)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: ft(0:1), fpc0, fpc1(3), tinv
      real(8) :: r0, r1(3)
!
      r0= zero
      do i= 1,3
        r1(i) = zero
      enddo
!
      do ijprim= 1,nij
        ex12= exfac(1,ijprim)
        ex21= exfac(2,ijprim)
        ex2i= exfac(3,ijprim)
        ex2j= exfac(4,ijprim)
        c12 = exfac(5,ijprim)
        fpc0= zero
        do i=1,3
          fpc1(i)= zero
        enddo
        do iatom= 1,natom
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
          ft(0)= ft(0)*znuc(iatom)
          ft(1)= ft(1)*znuc(iatom)
          fpc0= fpc0+ft(0)
          do i= 1,3
            fpc1(i)= fpc1(i)-ft(1)*pcxyz(i)
          enddo
        enddo
        r0= r0+c12*fpc0*ex2i
        do i= 1,3
          r1(i)= r1(i)+c12*fpc1(i)
        enddo
      enddo
      do i= 1,3
        cint1(i,1)=-(r1(i)+r0*xyz(i))
      enddo
!
      return
end


!----------------------------------------------------------------------------
  subroutine int1cpp(cint1,exfac,pijxyz,xyz,nij,coord,znuc,natom,mxprsh)
!----------------------------------------------------------------------------
!
! Calculate 1-electron Coulomb integrals of (p|Z/r|p)
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nij, natom, mxprsh
      integer :: ijprim, i, iatom, igrid, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, two=2.0D+00
      real(8),parameter :: three=3.0D+00, p15=1.5D+00
      real(8),parameter :: pi=3.141592653589793D+00
      real(8),intent(in) :: exfac(5,mxprsh*mxprsh), pijxyz(3,mxprsh*mxprsh), xyz(3)
      real(8),intent(in) :: coord(3,natom), znuc(natom)
      real(8),intent(out) :: cint1(6,6)
      real(8) :: rc, ex12, ex21, ex2i, ex2j, c12, pcxyz(3)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: ft(0:2), fpc0, fpc1(4), fpc2(6)
      real(8) :: tinv, r0(2), r1(7), r2(6)
!
      do i= 1,2
        r0(i) = zero
      enddo
      do i= 1,7
        r1(i) = zero
      enddo
      do i= 1,6
        r2(i) = zero
      enddo
!
      do ijprim= 1,nij
        ex12= exfac(1,ijprim)
        ex21= exfac(2,ijprim)
        ex2i= exfac(3,ijprim)
        ex2j= exfac(4,ijprim)
        c12 = exfac(5,ijprim)
        fpc0= zero
        do i= 1,4
          fpc1(i)= zero
        enddo
        do i= 1,6
          fpc2(i)= zero
        enddo
        do iatom= 1,natom
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
          do i= 0,2
            ft(i)= ft(i)*znuc(iatom) 
          enddo
          fpc0= fpc0+ft(0)
          do i= 1,3
            fpc1(i)= fpc1(i)+ft(1)*pcxyz(i) 
          enddo
          fpc1(4)= fpc1(4)+ft(1)
          fpc2(1)= fpc2(1)+ft(2)*pcxyz(1)*pcxyz(1)
          fpc2(2)= fpc2(2)+ft(2)*pcxyz(2)*pcxyz(2)
          fpc2(3)= fpc2(3)+ft(2)*pcxyz(3)*pcxyz(3)
          fpc2(4)= fpc2(4)+ft(2)*pcxyz(1)*pcxyz(2)
          fpc2(5)= fpc2(5)+ft(2)*pcxyz(1)*pcxyz(3)
          fpc2(6)= fpc2(6)+ft(2)*pcxyz(2)*pcxyz(3)
        enddo
        fpc0= fpc0*c12
        r0(1)= r0(1)-ex2i*ex2j*fpc0
        r0(2)= r0(2)+half*ex21*fpc0
        do i= 1,3
          r1(i)  = r1(i)  +c12*ex2j*fpc1(i)
          r1(i+3)= r1(i+3)-c12*ex2i*fpc1(i)
        enddo
        r1(7)= r1(7)-half*c12*ex21*fpc1(4)
        do i= 1,6
          r2(i)= r2(i)+c12*fpc2(i)
        enddo
      enddo
!
      cint1(1,1)=-(r2(1)+r1(1)*xyz(1)+r1(4)*xyz(1)+r1(7)+r0(1)*xyz(1)*xyz(1)+r0(2))
      cint1(2,1)=-(r2(4)+r1(2)*xyz(1)+r1(4)*xyz(2)      +r0(1)*xyz(1)*xyz(2)      )
      cint1(3,1)=-(r2(5)+r1(3)*xyz(1)+r1(4)*xyz(3)      +r0(1)*xyz(1)*xyz(3)      )
      cint1(1,2)=-(r2(4)+r1(1)*xyz(2)+r1(5)*xyz(1)      +r0(1)*xyz(1)*xyz(2)      )
      cint1(2,2)=-(r2(2)+r1(2)*xyz(2)+r1(5)*xyz(2)+r1(7)+r0(1)*xyz(2)*xyz(2)+r0(2))
      cint1(3,2)=-(r2(6)+r1(3)*xyz(2)+r1(5)*xyz(3)      +r0(1)*xyz(2)*xyz(3)      )
      cint1(1,3)=-(r2(5)+r1(1)*xyz(3)+r1(6)*xyz(1)      +r0(1)*xyz(1)*xyz(3)      )
      cint1(2,3)=-(r2(6)+r1(2)*xyz(3)+r1(6)*xyz(2)      +r0(1)*xyz(2)*xyz(3)      )
      cint1(3,3)=-(r2(3)+r1(3)*xyz(3)+r1(6)*xyz(3)+r1(7)+r0(1)*xyz(3)*xyz(3)+r0(2))
!
      return
end


!---------------------------------------------------------------------------------
  subroutine int1cds(cint1,exfac,pijxyz,xyz,nij,coord,znuc,natom,mxprsh,nbfij)
!---------------------------------------------------------------------------------
!
! Calculate 1-electron Coulomb integrals of (d|Z/r|s)
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nij, natom, mxprsh, nbfij(2)
      integer :: ijprim, i, iatom, igrid, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, two=2.0D+00
      real(8),parameter :: three=3.0D+00, p15=1.5D+00, pi=3.141592653589793D+00
      real(8),parameter :: sqrt3=1.732050807568877D+00, sqrt3h=8.660254037844386D-01
      real(8),intent(in) :: exfac(5,mxprsh*mxprsh), pijxyz(3,mxprsh*mxprsh), xyz(3)
      real(8),intent(in) :: coord(3,natom), znuc(natom)
      real(8),intent(out) :: cint1(6,6)
      real(8) :: rc, ex12, ex21, ex2i, ex2j, c12, pcxyz(3)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: ft(0:2), fpc0, fpc1(4), fpc2(6)
      real(8) :: tinv, r0(2), r1(4), r2(6), ctmp(6)
!
      do i= 1,2
        r0(i) = zero
      enddo
      do i= 1,4
        r1(i) = zero
      enddo
      do i= 1,6
        r2(i) = zero
      enddo
!
      do ijprim= 1,nij
        ex12= exfac(1,ijprim)
        ex21= exfac(2,ijprim)
        ex2i= exfac(3,ijprim)
        ex2j= exfac(4,ijprim)
        c12 = exfac(5,ijprim)
        fpc0= zero
        do i= 1,4
          fpc1(i)= zero
        enddo
        do i= 1,6
          fpc2(i)= zero
        enddo
        do iatom= 1,natom
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
          do i= 0,2
            ft(i)= ft(i)*znuc(iatom) 
          enddo
          fpc0= fpc0+ft(0)
          do i= 1,3
            fpc1(i)= fpc1(i)+ft(1)*pcxyz(i) 
          enddo
          fpc1(4)= fpc1(4)+ft(1)
          fpc2(1)= fpc2(1)+ft(2)*pcxyz(1)*pcxyz(1)
          fpc2(2)= fpc2(2)+ft(2)*pcxyz(2)*pcxyz(2)
          fpc2(3)= fpc2(3)+ft(2)*pcxyz(3)*pcxyz(3)
          fpc2(4)= fpc2(4)+ft(2)*pcxyz(1)*pcxyz(2)
          fpc2(5)= fpc2(5)+ft(2)*pcxyz(1)*pcxyz(3)
          fpc2(6)= fpc2(6)+ft(2)*pcxyz(2)*pcxyz(3)
        enddo
        fpc0= fpc0*c12
        r0(1)= r0(1)+ex2i*ex2i*fpc0
        r0(2)= r0(2)+half*ex21*fpc0
        do i= 1,3
          r1(i)= r1(i)-c12*ex2i*fpc1(i)
        enddo
        r1(4)= r1(4)-half*c12*ex21*fpc1(4)
        do i= 1,6
          r2(i)= r2(i)+c12*fpc2(i)
        enddo
      enddo
!
      ctmp(1)=-(r2(1)+r1(1)*xyz(1)*two+r1(4)   +r0(1)*xyz(1)*xyz(1)+r0(2))
      ctmp(2)=-(r2(2)+r1(2)*xyz(2)*two+r1(4)   +r0(1)*xyz(2)*xyz(2)+r0(2))
      ctmp(3)=-(r2(3)+r1(3)*xyz(3)*two+r1(4)   +r0(1)*xyz(3)*xyz(3)+r0(2))
      ctmp(4)=-(r2(4)+r1(1)*xyz(2)+r1(2)*xyz(1)+r0(1)*xyz(1)*xyz(2)      )
      ctmp(5)=-(r2(5)+r1(1)*xyz(3)+r1(3)*xyz(1)+r0(1)*xyz(1)*xyz(3)      )
      ctmp(6)=-(r2(6)+r1(2)*xyz(3)+r1(3)*xyz(2)+r0(1)*xyz(2)*xyz(3)      )
!
      if(nbfij(2) == 6) then
        cint1(1,1)= ctmp(1)
        cint1(2,1)= ctmp(4)*sqrt3
        cint1(3,1)= ctmp(5)*sqrt3
        cint1(4,1)= ctmp(2)
        cint1(5,1)= ctmp(6)*sqrt3
        cint1(6,1)= ctmp(3)
      else
        cint1(1,1)= ctmp(4)*sqrt3
        cint1(2,1)= ctmp(6)*sqrt3
        cint1(3,1)= ctmp(3)-(ctmp(1)+ctmp(2))*half
        cint1(4,1)= ctmp(5)*sqrt3
        cint1(5,1)=(ctmp(1)-ctmp(2))*sqrt3h
      endif
!
      return
end


!---------------------------------------------------------------------------------
  subroutine int1cdp(cint1,exfac,pijxyz,xyz,nij,coord,znuc,natom,mxprsh,nbfij)
!---------------------------------------------------------------------------------
!
! Calculate 1-electron Coulomb integrals of (d|Z/r|p)
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nij, natom, mxprsh, nbfij(2)
      integer :: ijprim, i, iatom, igrid, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, two=2.0D+00
      real(8),parameter :: three=3.0D+00, five=5.0D+00, p15=1.5D+00, p25=2.5D+00
      real(8),parameter :: pi=3.141592653589793D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrt3=1.732050807568877D+00
      real(8),intent(in) :: exfac(5,mxprsh*mxprsh), pijxyz(3,mxprsh*mxprsh), xyz(3)
      real(8),intent(in) :: coord(3,natom), znuc(natom)
      real(8),intent(out) :: cint1(6,6)
      real(8) :: rc, ex12, ex21, ex2i, ex2j, c12, pcxyz(3)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: ft(0:3), fpc0, fpc1(4), fpc2(9)
      real(8) :: fpc3(10), pxx, pyy, pzz, pxy, pxz, pyz
      real(8) :: tinv, r0(3), r1(11), r2(15), r3(10), xx, yy, zz, xy, xz, yz, ctmp(6,3)
!
      do i= 1,3
        r0(i) = zero
      enddo
      do i= 1,11
        r1(i) = zero
      enddo
      do i= 1,15
        r2(i) = zero
      enddo
      do i= 1,10
        r3(i) = zero
      enddo
      xx= xyz(1)*xyz(1)
      yy= xyz(2)*xyz(2)
      zz= xyz(3)*xyz(3)
      xy= xyz(1)*xyz(2)
      xz= xyz(1)*xyz(3)
      yz= xyz(2)*xyz(3)
!
      do ijprim= 1,nij
        ex12= exfac(1,ijprim)
        ex21= exfac(2,ijprim)
        ex2i= exfac(3,ijprim)
        ex2j= exfac(4,ijprim)
        c12 = exfac(5,ijprim)
        fpc0= zero
        do i= 1,4
          fpc1(i)= zero
        enddo
        do i= 1,9
          fpc2(i)= zero
        enddo
        do i= 1,10
          fpc3(i)= zero
        enddo
!
        do iatom= 1,natom
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
          do i= 0,3
            ft(i)= ft(i)*znuc(iatom) 
          enddo
          fpc0= fpc0+ft(0)
          do i= 1,3
            fpc1(i)= fpc1(i)+ft(1)*pcxyz(i) 
          enddo
          fpc1(4)= fpc1(4)+ft(1)
          fpc2(1)= fpc2(1)+ft(2)*pxx
          fpc2(2)= fpc2(2)+ft(2)*pyy
          fpc2(3)= fpc2(3)+ft(2)*pzz
          fpc2(4)= fpc2(4)+ft(2)*pxy
          fpc2(5)= fpc2(5)+ft(2)*pxz
          fpc2(6)= fpc2(6)+ft(2)*pyz
          fpc2(7)= fpc2(7)+ft(2)*pcxyz(1)
          fpc2(8)= fpc2(8)+ft(2)*pcxyz(2)
          fpc2(9)= fpc2(9)+ft(2)*pcxyz(3)
          fpc3(1) = fpc3(1) +ft(3)*pxx*pcxyz(1)
          fpc3(2) = fpc3(2) +ft(3)*pyy*pcxyz(2)
          fpc3(3) = fpc3(3) +ft(3)*pzz*pcxyz(3)
          fpc3(4) = fpc3(4) +ft(3)*pxx*pcxyz(2)
          fpc3(5) = fpc3(5) +ft(3)*pxx*pcxyz(3)
          fpc3(6) = fpc3(6) +ft(3)*pxy*pcxyz(2)
          fpc3(7) = fpc3(7) +ft(3)*pyy*pcxyz(3)
          fpc3(8) = fpc3(8) +ft(3)*pxz*pcxyz(3)
          fpc3(9) = fpc3(9) +ft(3)*pyz*pcxyz(3)
          fpc3(10)= fpc3(10)+ft(3)*pxy*pcxyz(3)
        enddo
!
        fpc0= fpc0*c12
        r0(1)= r0(1)-ex2i*ex2i*ex2j*fpc0
        r0(2)= r0(2)-half*ex2j*ex21*fpc0
        r0(3)= r0(3)+half*ex2i*ex21*fpc0
        do i= 1,3
          r1(i)  = r1(i)  +c12*ex2i*ex2j*fpc1(i)
          r1(i+3)= r1(i+3)-c12*ex2i*ex2i*fpc1(i)
          r1(i+6)= r1(i+6)-half*c12*ex21*fpc1(i)
        enddo
        r1(10)= r1(10)+half*c12*ex2j*ex21*fpc1(4)
        r1(11)= r1(11)-half*c12*ex2i*ex21*fpc1(4)
        do i= 1,6
          r2(i)=   r2(i)  -c12*ex2j*fpc2(i)
          r2(i+6)= r2(i+6)+c12*ex2i*fpc2(i)
        enddo
        do i= 1,3
          r2(i+12)= r2(i+12)+half*c12*ex21*fpc2(i+6)
        enddo
        do i= 1,10
          r3(i)= r3(i)-c12*fpc3(i)
        enddo
      enddo
!
      ctmp(1,1)=-(r3(1)+r2(1)*xyz(1)+r2(7)*xyz(1)*two+r2(13)*three+r1(1)*xx*two+r1(4)*xx &
&                +r1(7)*three+r1(10)*xyz(1)+r1(11)*xyz(1)*two+r0(1)*xyz(1)*xx &
&                +r0(2)*xyz(1)+r0(3)*xyz(1)*two)
      ctmp(2,1)=-(r3(6)+r2(2)*xyz(1)+r2(10)*xyz(2)*two+r2(13)+r1(2)*xy*two+r1(4)*yy &
&                +r1(7)+r1(10)*xyz(1)+r0(1)*xyz(1)*yy+r0(2)*xyz(1))
      ctmp(3,1)=-(r3(8)+r2(3)*xyz(1)+r2(11)*xyz(3)*two+r2(13)+r1(3)*xz*two+r1(4)*zz &
&                +r1(7)+r1(10)*xyz(1)+r0(1)*xyz(1)*zz+r0(2)*xyz(1))
      ctmp(4,1)=-(r3(4)+r2(4)*xyz(1)+r2(7)*xyz(2)+r2(10)*xyz(1)+r2(14)+r1(1)*xy+r1(2)*xx &
&                +r1(4)*xy+r1(8)+r1(11)*xyz(2)+r0(1)*xyz(1)*xy+r0(3)*xyz(2))
      ctmp(5,1)=-(r3(5)+r2(5)*xyz(1)+r2(7)*xyz(3)+r2(11)*xyz(1)+r2(15)+r1(1)*xz+r1(3)*xx &
&                +r1(4)*xz+r1(9)+r1(11)*xyz(3)+r0(1)*xyz(1)*xz+r0(3)*xyz(3))
      ctmp(6,1)=-(r3(10)+r2(6)*xyz(1)+r2(10)*xyz(3)+r2(11)*xyz(2)+r1(2)*xz+r1(3)*xy &
&                +r1(4)*yz+r0(1)*xyz(1)*yz)
      ctmp(1,2)=-(r3(4)+r2(1)*xyz(2)+r2(10)*xyz(1)*two+r2(14)+r1(1)*xy*two+r1(5)*xx &
&                +r1(8)+r1(10)*xyz(2)+r0(1)*xyz(2)*xx+r0(2)*xyz(2))
      ctmp(2,2)=-(r3(2)+r2(2)*xyz(2)+r2(8)*xyz(2)*two+r2(14)*three+r1(2)*yy*two+r1(5)*yy &
&                +r1(8)*three+r1(10)*xyz(2)+r1(11)*xyz(2)*two+r0(1)*xyz(2)*yy &
&                +r0(2)*xyz(2)+r0(3)*xyz(2)*two)
      ctmp(3,2)=-(r3(9)+r2(3)*xyz(2)+r2(12)*xyz(3)*two+r2(14)+r1(3)*yz*two+r1(5)*zz &
&                +r1(8)+r1(10)*xyz(2)+r0(1)*xyz(2)*zz+r0(2)*xyz(2))
      ctmp(4,2)=-(r3(6)+r2(4)*xyz(2)+r2(8)*xyz(1)+r2(10)*xyz(2)+r2(13)+r1(1)*yy+r1(2)*xy &
&                +r1(5)*xy+r1(7)+r1(11)*xyz(1)+r0(1)*xyz(2)*xy+r0(3)*xyz(1))
      ctmp(5,2)=-(r3(10)+r2(5)*xyz(2)+r2(10)*xyz(3)+r2(12)*xyz(1)+r1(1)*yz+r1(3)*xy &
&                +r1(5)*xz+r0(1)*xyz(2)*xz)
      ctmp(6,2)=-(r3(7)+r2(6)*xyz(2)+r2(8)*xyz(3)+r2(12)*xyz(2)+r2(15)+r1(2)*yz+r1(3)*yy &
&                +r1(5)*yz+r1(9)+r1(11)*xyz(3)+r0(1)*xyz(2)*yz+r0(3)*xyz(3))
      ctmp(1,3)=-(r3(5)+r2(1)*xyz(3)+r2(11)*xyz(1)*two+r2(15)+r1(1)*xz*two+r1(6)*xx &
&                +r1(9)+r1(10)*xyz(3)+r0(1)*xyz(3)*xx+r0(2)*xyz(3))
      ctmp(2,3)=-(r3(7)+r2(2)*xyz(3)+r2(12)*xyz(2)*two+r2(15)+r1(2)*yz*two+r1(6)*yy &
&                +r1(9)+r1(10)*xyz(3)+r0(1)*xyz(3)*yy+r0(2)*xyz(3))
      ctmp(3,3)=-(r3(3)+r2(3)*xyz(3)+r2(9)*xyz(3)*two+r2(15)*three+r1(3)*zz*two+r1(6)*zz &
&                +r1(9)*three+r1(10)*xyz(3)+r1(11)*xyz(3)*two+r0(1)*xyz(3)*zz &
&                +r0(2)*xyz(3)+r0(3)*xyz(3)*two)
      ctmp(4,3)=-(r3(10)+r2(4)*xyz(3)+r2(11)*xyz(2)+r2(12)*xyz(1)+r1(1)*yz+r1(2)*xz &
&                +r1(6)*xy+r0(1)*xyz(3)*xy)
      ctmp(5,3)=-(r3(8)+r2(5)*xyz(3)+r2(9)*xyz(1)+r2(11)*xyz(3)+r2(13)+r1(1)*zz+r1(3)*xz &
&                +r1(6)*xz+r1(7)+r1(11)*xyz(1)+r0(1)*xyz(3)*xz+r0(3)*xyz(1))
      ctmp(6,3)=-(r3(9)+r2(6)*xyz(3)+r2(9)*xyz(2)+r2(12)*xyz(3)+r2(14)+r1(2)*zz+r1(3)*yz &
&                +r1(6)*yz+r1(8)+r1(11)*xyz(2)+r0(1)*xyz(3)*yz+r0(3)*xyz(2))
!
      if(nbfij(2) == 6) then
        do i= 1,3
          cint1(1,i)= ctmp(1,i)
          cint1(2,i)= ctmp(4,i)*sqrt3
          cint1(3,i)= ctmp(5,i)*sqrt3
          cint1(4,i)= ctmp(2,i)
          cint1(5,i)= ctmp(6,i)*sqrt3
          cint1(6,i)= ctmp(3,i)
        enddo
      else
        do i= 1,3
          cint1(1,i)= ctmp(4,i)*sqrt3
          cint1(2,i)= ctmp(6,i)*sqrt3
          cint1(3,i)= ctmp(3,i)-(ctmp(1,i)+ctmp(2,i))*half
          cint1(4,i)= ctmp(5,i)*sqrt3
          cint1(5,i)=(ctmp(1,i)-ctmp(2,i))*sqrt3h
        enddo
      endif
!
      return
end


!---------------------------------------------------------------------------------
  subroutine int1cdd(cint1,exfac,pijxyz,xyz,nij,coord,znuc,natom,mxprsh,nbfij)
!---------------------------------------------------------------------------------
!
! Calculate 1-electron Coulomb integrals of (d|Z/r|d)
!
      use fmtgrid, only : fgrid, threshtval
      implicit none
      integer,intent(in) :: nij, natom, mxprsh, nbfij(2)
      integer :: ijprim, i, j, iatom, igrid, ii
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, half=0.5D+00, two=2.0D+00
      real(8),parameter :: three=3.0D+00, four=4.0D+00, five=5.0D+00, six=6.0D+00
      real(8),parameter :: seven=7.0D+00, p15=1.5D+00, p25=2.5D+00, p35=3.5D+00
      real(8),parameter :: pi=3.141592653589793D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrt3=1.732050807568877D+00
      real(8),intent(in) :: exfac(5,mxprsh*mxprsh), pijxyz(3,mxprsh*mxprsh), xyz(3)
      real(8),intent(in) :: coord(3,natom), znuc(natom)
      real(8),intent(out) :: cint1(6,6)
      real(8) :: rc, ex12, ex21, ex2i, ex2j, c12, pcxyz(3)
      real(8) :: tval, tval2, tval3, tval4, tval5, tval6, tval7, tval8, tval9, tval10
      real(8) :: ft(0:4), fpc0, fpc1(4), fpc2(10)
      real(8) :: fpc3(16), fpc4(15), pxx, pyy, pzz, pxy, pxz, pyz
      real(8) :: tinv, r0(5), r1(16), r2(31), r3(26), r4(15), xx, yy, zz, xy, xz, yz
      real(8) :: ctmp(6,6), ctmp2(6,6)
!
      do i= 1,5
        r0(i) = zero
      enddo
      do i= 1,16
        r1(i) = zero
      enddo
      do i= 1,31
        r2(i) = zero
      enddo
      do i= 1,26
        r3(i) = zero
      enddo
      do i= 1,15
        r4(i) = zero
      enddo
      xx= xyz(1)*xyz(1)
      yy= xyz(2)*xyz(2)
      zz= xyz(3)*xyz(3)
      xy= xyz(1)*xyz(2)
      xz= xyz(1)*xyz(3)
      yz= xyz(2)*xyz(3)
!
      do ijprim= 1,nij
        ex12= exfac(1,ijprim)
        ex21= exfac(2,ijprim)
        ex2i= exfac(3,ijprim)
        ex2j= exfac(4,ijprim)
        c12 = exfac(5,ijprim)
        fpc0= zero
        do i= 1,4
          fpc1(i)= zero
        enddo
        do i= 1,10
          fpc2(i)= zero
        enddo
        do i= 1,16
          fpc3(i)= zero
        enddo
        do i= 1,15
          fpc4(i)= zero
        enddo
!
        do iatom= 1,natom
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
          do i= 0,4
            ft(i)= ft(i)*znuc(iatom) 
          enddo
          fpc0= fpc0+ft(0)
          do i= 1,3
            fpc1(i)= fpc1(i)+ft(1)*pcxyz(i) 
          enddo
          fpc1(4)= fpc1(4)+ft(1)
          fpc2(1) = fpc2(1) +ft(2)*pxx
          fpc2(2) = fpc2(2) +ft(2)*pyy
          fpc2(3) = fpc2(3) +ft(2)*pzz
          fpc2(4) = fpc2(4) +ft(2)*pxy
          fpc2(5) = fpc2(5) +ft(2)*pxz
          fpc2(6) = fpc2(6) +ft(2)*pyz
          fpc2(7) = fpc2(7) +ft(2)*pcxyz(1)
          fpc2(8) = fpc2(8) +ft(2)*pcxyz(2)
          fpc2(9) = fpc2(9) +ft(2)*pcxyz(3)
          fpc2(10)= fpc2(10)+ft(2)
          fpc3(1) = fpc3(1) +ft(3)*pxx*pcxyz(1)
          fpc3(2) = fpc3(2) +ft(3)*pyy*pcxyz(2)
          fpc3(3) = fpc3(3) +ft(3)*pzz*pcxyz(3)
          fpc3(4) = fpc3(4) +ft(3)*pxx*pcxyz(2)
          fpc3(5) = fpc3(5) +ft(3)*pxx*pcxyz(3)
          fpc3(6) = fpc3(6) +ft(3)*pxy*pcxyz(2)
          fpc3(7) = fpc3(7) +ft(3)*pyy*pcxyz(3)
          fpc3(8) = fpc3(8) +ft(3)*pxz*pcxyz(3)
          fpc3(9) = fpc3(9) +ft(3)*pyz*pcxyz(3)
          fpc3(10)= fpc3(10)+ft(3)*pxy*pcxyz(3)
          fpc3(11)= fpc3(11)+ft(3)*pxx
          fpc3(12)= fpc3(12)+ft(3)*pyy
          fpc3(13)= fpc3(13)+ft(3)*pzz
          fpc3(14)= fpc3(14)+ft(3)*pxy
          fpc3(15)= fpc3(15)+ft(3)*pxz
          fpc3(16)= fpc3(16)+ft(3)*pyz
          fpc4(1) = fpc4(1) +ft(4)*pxx*pxx
          fpc4(2) = fpc4(2) +ft(4)*pyy*pyy
          fpc4(3) = fpc4(3) +ft(4)*pzz*pzz
          fpc4(4) = fpc4(4) +ft(4)*pxx*pxy
          fpc4(5) = fpc4(5) +ft(4)*pxx*pxz
          fpc4(6) = fpc4(6) +ft(4)*pxy*pyy
          fpc4(7) = fpc4(7) +ft(4)*pyy*pyz
          fpc4(8) = fpc4(8) +ft(4)*pxz*pzz
          fpc4(9) = fpc4(9) +ft(4)*pyz*pzz
          fpc4(10)= fpc4(10)+ft(4)*pxx*pyy
          fpc4(11)= fpc4(11)+ft(4)*pxx*pzz
          fpc4(12)= fpc4(12)+ft(4)*pyy*pzz
          fpc4(13)= fpc4(13)+ft(4)*pxx*pyz
          fpc4(14)= fpc4(14)+ft(4)*pxy*pyz
          fpc4(15)= fpc4(15)+ft(4)*pxy*pzz
        enddo
!
        fpc0= fpc0*c12
        r0(1)= r0(1)+ex2i*ex2i*ex2j*ex2j*fpc0
        r0(2)= r0(2)+half*ex2j*ex2j*ex21*fpc0
        r0(3)= r0(3)-half*ex2i*ex2j*ex21*fpc0
        r0(4)= r0(4)+half*ex2i*ex2i*ex21*fpc0
        r0(5)= r0(5)+half*half*ex21*ex21*fpc0
        do i= 1,4
          fpc1(i)= fpc1(i)*c12
        enddo
        do i= 1,3
          r1(i  )= r1(i  )-ex2i*ex2j*ex2j*fpc1(i)
          r1(i+3)= r1(i+3)+ex2i*ex2i*ex2j*fpc1(i)
          r1(i+6)= r1(i+6)+half*ex2j*ex21*fpc1(i)
          r1(i+9)= r1(i+9)-half*ex2i*ex21*fpc1(i)
        enddo
        r1(13)= r1(13)-half*ex2j*ex2j*ex21*fpc1(4)
        r1(14)= r1(14)+half*ex2i*ex2j*ex21*fpc1(4)
        r1(15)= r1(15)-half*ex2i*ex2i*ex21*fpc1(4)
        r1(16)= r1(16)-half*half*ex21*ex21*fpc1(4)
        do i= 1,10
          fpc2(i)= fpc2(i)*c12
        enddo
        do i= 1,6
          r2(i   )= r2(i   )+ex2j*ex2j*fpc2(i)
          r2(i+ 6)= r2(i+ 6)-ex2i*ex2j*fpc2(i)
          r2(i+12)= r2(i+12)+ex2i*ex2i*fpc2(i)
          r2(i+18)= r2(i+18)+half*ex21*fpc2(i)
        enddo
        do i= 1,3
          r2(i+24)= r2(i+24)-half*ex2j*ex21*fpc2(i+6)
          r2(i+27)= r2(i+27)+half*ex2i*ex21*fpc2(i+6)
        enddo
        r2(31)= r2(31)+half*half*ex21*ex21*fpc2(10)
        do i= 1,16
          fpc3(i)= fpc3(i)*c12
        enddo
        do i= 1,10
          r3(i   )= r3(i   )+ex2j*fpc3(i)
          r3(i+10)= r3(i+10)-ex2i*fpc3(i)
        enddo
        do i= 1,6
          r3(i+20)= r3(i+20)-half*ex21*fpc3(i+10)
        enddo
        do i= 1,15
          r4(i)= r4(i)+fpc4(i)*c12
        enddo
      enddo
!
      ctmp(1,1)=-(r4(1)+r3(1)*xyz(1)*two+r3(11)*xyz(1)*two+r3(21)*six+r2(1)*xx &
&                +r2(7)*xx*four+r2(13)*xx+r2(19)*six+r2(25)*xyz(1)*six+r2(28)*xyz(1)*six &
&                +r2(31)*three+r1(1)*xx*xyz(1)*two+r1(4)*xx*xyz(1)*two+r1(7)*xyz(1)*six &
&                +r1(10)*xyz(1)*six+r1(13)*xx+r1(14)*xx*four+r1(15)*xx+r1(16)*six &
&                +r0(1)*xx*xx+r0(2)*xx+r0(3)*xx*four+r0(4)*xx+r0(5)*three)
      ctmp(2,1)=-(r4(10)+r3(6)*xyz(1)*two+r3(14)*xyz(2)*two+r3(21)+r3(22)+r2(2)*xx &
&                +r2(10)*xy*four+r2(13)*yy+r2(19)+r2(20)+r2(25)*xyz(1)*two &
&                +r2(29)*xyz(2)*two+r2(31)+r1(2)*xx*xyz(2)*two+r1(4)*yy*xyz(1)*two &
&                +r1(7)*xyz(1)*two+r1(11)*xyz(2)*two+r1(13)*xx+r1(15)*yy+r1(16)*two &
&                +r0(1)*xx*yy+r0(2)*xx+r0(4)*yy+r0(5))
      ctmp(3,1)=-(r4(11)+r3(8)*xyz(1)*two+r3(15)*xyz(3)*two+r3(21)+r3(23)+r2(3)*xx &
&                +r2(11)*xz*four+r2(13)*zz+r2(19)+r2(21)+r2(25)*xyz(1)*two &
&                +r2(30)*xyz(3)*two+r2(31)+r1(3)*xx*xyz(3)*two+r1(4)*zz*xyz(1)*two &
&                +r1(7)*xyz(1)*two+r1(12)*xyz(3)*two+r1(13)*xx+r1(15)*zz+r1(16)*two &
&                +r0(1)*xx*zz+r0(2)*xx+r0(4)*zz+r0(5))
      ctmp(4,1)=-(r4(4)+r3(4)*xyz(1)*two+r3(11)*xyz(2)+r3(14)*xyz(1)+r3(24)*three &
&                +r2(4)*xx+r2(7)*xy*two+r2(10)*xx*two+r2(13)*xy+r2(22)*three &
&                +r2(26)*xyz(1)*two+r2(28)*xyz(2)*three+r2(29)*xyz(1)+r1(1)*xx*xyz(2) &
&                +r1(2)*xx*xyz(1)+r1(4)*xy*xyz(1)*two+r1(8)*xyz(1)*two+r1(10)*xyz(2)*three &
&                +r1(11)*xyz(1)+r1(14)*xy*two+r1(15)*xy+r0(1)*xx*xy+r0(3)*xy*two+r0(4)*xy)
      ctmp(5,1)=-(r4(5)+r3(5)*xyz(1)*two+r3(11)*xyz(3)+r3(15)*xyz(1)+r3(25)*three+r2(5)*xx &
&                +r2(7)*xz*two+r2(11)*xx*two+r2(13)*xz+r2(23)*three+r2(27)*xyz(1)*two &
&                +r2(28)*xyz(3)*three+r2(30)*xyz(1)+r1(1)*xx*xyz(3)+r1(3)*xx*xyz(1) &
&                +r1(4)*xz*xyz(1)*two+r1(9)*xyz(1)*two+r1(10)*xyz(3)*three+r1(12)*xyz(1) &
&                +r1(14)*xz*two+r1(15)*xz+r0(1)*xx*xz+r0(3)*xz*two+r0(4)*xz)
      ctmp(6,1)=-(r4(13)+r3(10)*xyz(1)*two+r3(14)*xyz(3)+r3(15)*xyz(2)+r3(26)+r2(6)*xx &
&                +r2(10)*xz*two+r2(11)*xy*two+r2(13)*yz+r2(24)+r2(29)*xyz(3)+r2(30)*xyz(2) &
&                +r1(2)*xx*xyz(3)+r1(3)*xx*xyz(2)+r1(4)*yz*xyz(1)*two+r1(11)*xyz(3) &
&                +r1(12)*xyz(2)+r1(15)*yz+r0(1)*xx*yz+r0(4)*yz)
      ctmp(1,2)=-(r4(10)+r3(4)*xyz(2)*two+r3(16)*xyz(1)*two+r3(21)+r3(22)+r2(1)*yy &
&                +r2(10)*xy*four+r2(14)*xx+r2(19)+r2(20)+r2(26)*xyz(2)*two &
&                +r2(28)*xyz(1)*two+r2(31)+r1(1)*yy*xyz(1)*two+r1(5)*xx*xyz(2)*two &
&                +r1(8)*xyz(2)*two+r1(10)*xyz(1)*two+r1(13)*yy+r1(15)*xx+r1(16)*two &
&                +r0(1)*yy*xx+r0(2)*yy+r0(4)*xx+r0(5))
      ctmp(2,2)=-(r4(2)+r3(2)*xyz(2)*two+r3(12)*xyz(2)*two+r3(22)*six+r2(2)*yy &
&                +r2(8)*yy*four+r2(14)*yy+r2(20)*six+r2(26)*xyz(2)*six+r2(29)*xyz(2)*six &
&                +r2(31)*three+r1(2)*yy*xyz(2)*two+r1(5)*yy*xyz(2)*two+r1(8)*xyz(2)*six &
&                +r1(11)*xyz(2)*six+r1(13)*yy+r1(14)*yy*four+r1(15)*yy+r1(16)*six &
&                +r0(1)*yy*yy+r0(2)*yy+r0(3)*yy*four+r0(4)*yy+r0(5)*three)
      ctmp(3,2)=-(r4(12)+r3(9)*xyz(2)*two+r3(17)*xyz(3)*two+r3(22)+r3(23)+r2(3)*yy &
&                +r2(12)*yz*four+r2(14)*zz+r2(20)+r2(21)+r2(26)*xyz(2)*two &
&                +r2(30)*xyz(3)*two+r2(31)+r1(3)*yy*xyz(3)*two+r1(5)*zz*xyz(2)*two &
&                +r1(8)*xyz(2)*two+r1(12)*xyz(3)*two+r1(13)*yy+r1(15)*zz+r1(16)*two &
&                +r0(1)*yy*zz+r0(2)*yy+r0(4)*zz+r0(5))
      ctmp(4,2)=-(r4(6)+r3(6)*xyz(2)*two+r3(12)*xyz(1)+r3(16)*xyz(2)+r3(24)*three+r2(4)*yy &
&                +r2(8)*xy*two+r2(10)*yy*two+r2(14)*xy+r2(22)*three+r2(25)*xyz(2)*two &
&                +r2(28)*xyz(2)+r2(29)*xyz(1)*three+r1(1)*yy*xyz(2)+r1(2)*yy*xyz(1) &
&                +r1(5)*xy*xyz(2)*two+r1(7)*xyz(2)*two+r1(10)*xyz(2)+r1(11)*xyz(1)*three &
&                +r1(14)*xy*two+r1(15)*xy+r0(1)*yy*xy+r0(3)*xy*two+r0(4)*xy)
      ctmp(5,2)=-(r4(14)+r3(10)*xyz(2)*two+r3(16)*xyz(3)+r3(17)*xyz(1)+r3(25)+r2(5)*yy &
&                +r2(10)*yz*two+r2(12)*xy*two+r2(14)*xz+r2(23)+r2(28)*xyz(3)+r2(30)*xyz(1) &
&                +r1(1)*yy*xyz(3)+r1(3)*yy*xyz(1)+r1(5)*xz*xyz(2)*two+r1(10)*xyz(3) &
&                +r1(12)*xyz(1)+r1(15)*xz+r0(1)*yy*xz+r0(4)*xz)
      ctmp(6,2)=-(r4(7)+r3(7)*xyz(2)*two+r3(12)*xyz(3)+r3(17)*xyz(2)+r3(26)*three+r2(6)*yy &
&                +r2(8)*yz*two+r2(12)*yy*two+r2(14)*yz+r2(24)*three+r2(27)*xyz(2)*two &
&                +r2(29)*xyz(3)*three+r2(30)*xyz(2)+r1(2)*yy*xyz(3)+r1(3)*yy*xyz(2) &
&                +r1(5)*yz*xyz(2)*two+r1(9)*xyz(2)*two+r1(11)*xyz(3)*three+r1(12)*xyz(2) &
&                +r1(14)*yz*two+r1(15)*yz+r0(1)*yy*yz+r0(3)*yz*two+r0(4)*yz)
      ctmp(1,3)=-(r4(11)+r3(5)*xyz(3)*two+r3(18)*xyz(1)*two+r3(21)+r3(23)+r2(1)*zz &
&                +r2(11)*xz*four+r2(15)*xx+r2(19)+r2(21)+r2(27)*xyz(3)*two &
&                +r2(28)*xyz(1)*two+r2(31)+r1(1)*zz*xyz(1)*two+r1(6)*xx*xyz(3)*two &
&                +r1(9)*xyz(3)*two+r1(10)*xyz(1)*two+r1(13)*zz+r1(15)*xx+r1(16)*two &
&                +r0(1)*zz*xx+r0(2)*zz+r0(4)*xx+r0(5))
      ctmp(2,3)=-(r4(12)+r3(7)*xyz(3)*two+r3(19)*xyz(2)*two+r3(22)+r3(23)+r2(2)*zz &
&                +r2(12)*yz*four+r2(15)*yy+r2(20)+r2(21)+r2(27)*xyz(3)*two &
&                +r2(29)*xyz(2)*two+r2(31)+r1(2)*zz*xyz(2)*two+r1(6)*yy*xyz(3)*two &
&                +r1(9)*xyz(3)*two+r1(11)*xyz(2)*two+r1(13)*zz+r1(15)*yy+r1(16)*two &
&                +r0(1)*zz*yy+r0(2)*zz+r0(4)*yy+r0(5))
      ctmp(3,3)=-(r4(3)+r3(3)*xyz(3)*two+r3(13)*xyz(3)*two+r3(23)*six+r2(3)*zz &
&                +r2(9)*zz*four+r2(15)*zz+r2(21)*six+r2(27)*xyz(3)*six+r2(30)*xyz(3)*six &
&                +r2(31)*three+r1(3)*zz*xyz(3)*two+r1(6)*zz*xyz(3)*two+r1(9)*xyz(3)*six &
&                +r1(12)*xyz(3)*six+r1(13)*zz+r1(14)*zz*four+r1(15)*zz+r1(16)*six &
&                +r0(1)*zz*zz+r0(2)*zz+r0(3)*zz*four+r0(4)*zz+r0(5)*three)
      ctmp(4,3)=-(r4(15)+r3(10)*xyz(3)*two+r3(18)*xyz(2)+r3(19)*xyz(1)+r3(24)+r2(4)*zz &
&                +r2(11)*yz*two+r2(12)*xz*two+r2(15)*xy+r2(22)+r2(28)*xyz(2) &
&                +r2(29)*xyz(1)+r1(1)*zz*xyz(2)+r1(2)*zz*xyz(1)+r1(6)*xy*xyz(3)*two &
&                +r1(10)*xyz(2)+r1(11)*xyz(1)+r1(15)*xy+r0(1)*zz*xy+r0(4)*xy)
      ctmp(5,3)=-(r4(8)+r3(8)*xyz(3)*two+r3(13)*xyz(1)+r3(18)*xyz(3)+r3(25)*three+r2(5)*zz &
&                +r2(9)*xz*two+r2(11)*zz*two+r2(15)*xz+r2(23)*three+r2(25)*xyz(3)*two &
&                +r2(28)*xyz(3)+r2(30)*xyz(1)*three+r1(1)*zz*xyz(3)+r1(3)*zz*xyz(1) &
&                +r1(6)*xz*xyz(3)*two+r1(7)*xyz(3)*two+r1(10)*xyz(3)+r1(12)*xyz(1)*three &
&                +r1(14)*xz*two+r1(15)*xz+r0(1)*zz*xz+r0(3)*xz*two+r0(4)*xz)
      ctmp(6,3)=-(r4(9)+r3(9)*xyz(3)*two+r3(13)*xyz(2)+r3(19)*xyz(3)+r3(26)*three+r2(6)*zz &
&                +r2(9)*yz*two+r2(12)*zz*two+r2(15)*yz+r2(24)*three+r2(26)*xyz(3)*two &
&                +r2(29)*xyz(3)+r2(30)*xyz(2)*three+r1(2)*zz*xyz(3)+r1(3)*zz*xyz(2) &
&                +r1(6)*yz*xyz(3)*two+r1(8)*xyz(3)*two+r1(11)*xyz(3)+r1(12)*xyz(2)*three &
&                +r1(14)*yz*two+r1(15)*yz+r0(1)*zz*yz+r0(3)*yz*two+r0(4)*yz)
      ctmp(1,4)=-(r4(4)+r3(1)*xyz(2)+r3(4)*xyz(1)+r3(14)*xyz(1)*two+r3(24)*three+r2(1)*xy &
&                +r2(7)*xy*two+r2(10)*xx*two+r2(16)*xx+r2(22)*three+r2(25)*xyz(2)*three &
&                +r2(26)*xyz(1)+r2(29)*xyz(1)*two+r1(1)*xy*xyz(1)*two+r1(4)*xx*xyz(2) &
&                +r1(5)*xx*xyz(1)+r1(7)*xyz(2)*three+r1(8)*xyz(1)+r1(11)*xyz(1)*two &
&                +r1(13)*xy+r1(14)*xy*two+r0(1)*xy*xx+r0(2)*xy+r0(3)*xy*two)
      ctmp(2,4)=-(r4(6)+r3(2)*xyz(1)+r3(6)*xyz(2)+r3(16)*xyz(2)*two+r3(24)*three+r2(2)*xy &
&                +r2(8)*xy*two+r2(10)*yy*two+r2(16)*yy+r2(22)*three+r2(25)*xyz(2) &
&                +r2(26)*xyz(1)*three+r2(28)*xyz(2)*two+r1(2)*xy*xyz(2)*two+r1(4)*yy*xyz(2) &
&                +r1(5)*yy*xyz(1)+r1(7)*xyz(2)+r1(8)*xyz(1)*three+r1(10)*xyz(2)*two &
&                +r1(13)*xy+r1(14)*xy*two+r0(1)*xy*yy+r0(2)*xy+r0(3)*xy*two)
      ctmp(3,4)=-(r4(15)+r3(8)*xyz(2)+r3(9)*xyz(1)+r3(20)*xyz(3)*two+r3(24)+r2(3)*xy &
&                +r2(11)*yz*two+r2(12)*xz*two+r2(16)*zz+r2(22)+r2(25)*xyz(2)+r2(26)*xyz(1) &
&                +r1(3)*xy*xyz(3)*two+r1(4)*zz*xyz(2)+r1(5)*zz*xyz(1)+r1(7)*xyz(2) &
&                +r1(8)*xyz(1)+r1(13)*xy+r0(1)*xy*zz+r0(2)*xy)
      ctmp(4,4)=-(r4(10)+r3(4)*xyz(2)+r3(6)*xyz(1)+r3(14)*xyz(2)+r3(16)*xyz(1)+r3(21) &
&                +r3(22)+r2(4)*xy+r2(7)*yy+r2(8)*xx+r2(10)*xy+r2(10)*xy+r2(16)*xy+r2(19) &
&                +r2(20)+r2(25)*xyz(1)+r2(26)*xyz(2)+r2(28)*xyz(1)+r2(29)*xyz(2)+r2(31) &
&                +r1(1)*xy*xyz(2)+r1(2)*xy*xyz(1)+r1(4)*xy*xyz(2)+r1(5)*xy*xyz(1) &
&                +r1(7)*xyz(1)+r1(8)*xyz(2)+r1(10)*xyz(1)+r1(11)*xyz(2)+r1(14)*xx+r1(14)*yy &
&               +r1(16)*two+r0(1)*xy*xy+r0(3)*xx+r0(3)*yy+r0(5))
      ctmp(5,4)=-(r4(13)+r3(5)*xyz(2)+r3(10)*xyz(1)+r3(14)*xyz(3)+r3(20)*xyz(1)+r3(26) &
&                +r2(5)*xy+r2(7)*yz+r2(10)*xz+r2(11)*xy+r2(12)*xx+r2(16)*xz+r2(24) &
&                +r2(27)*xyz(2)+r2(29)*xyz(3)+r1(1)*xy*xyz(3)+r1(3)*xy*xyz(1) &
&                +r1(4)*xz*xyz(2)+r1(5)*xz*xyz(1)+r1(9)*xyz(2)+r1(11)*xyz(3)+r1(14)*yz &
&                +r0(1)*xy*xz+r0(3)*yz)
      ctmp(6,4)=-(r4(14)+r3(7)*xyz(1)+r3(10)*xyz(2)+r3(16)*xyz(3)+r3(20)*xyz(2)+r3(25) &
&                +r2(6)*xy+r2(8)*xz+r2(10)*yz+r2(11)*yy+r2(12)*xy+r2(16)*yz+r2(23) &
&                +r2(27)*xyz(1)+r2(28)*xyz(3)+r1(2)*xy*xyz(3)+r1(3)*xy*xyz(2) &
&                +r1(4)*yz*xyz(2)+r1(5)*yz*xyz(1)+r1(9)*xyz(1)+r1(10)*xyz(3)+r1(14)*xz &
&                +r0(1)*xy*yz+r0(3)*xz)
      ctmp(1,5)=-(r4(5)+r3(1)*xyz(3)+r3(5)*xyz(1)+r3(15)*xyz(1)*two+r3(25)*three+r2(1)*xz &
&                +r2(7)*xz*two+r2(11)*xx*two+r2(17)*xx+r2(23)*three+r2(25)*xyz(3)*three &
&                +r2(27)*xyz(1)+r2(30)*xyz(1)*two+r1(1)*xz*xyz(1)*two+r1(4)*xx*xyz(3) &
&                +r1(6)*xx*xyz(1)+r1(7)*xyz(3)*three+r1(9)*xyz(1)+r1(12)*xyz(1)*two &
&                +r1(13)*xz+r1(14)*xz*two+r0(1)*xz*xx+r0(2)*xz+r0(3)*xz*two)
      ctmp(2,5)=-(r4(14)+r3(6)*xyz(3)+r3(7)*xyz(1)+r3(20)*xyz(2)*two+r3(25)+r2(2)*xz &
&                +r2(10)*yz*two+r2(12)*xy*two+r2(17)*yy+r2(23)+r2(25)*xyz(3)+r2(27)*xyz(1) &
&                +r1(2)*xz*xyz(2)*two+r1(4)*yy*xyz(3)+r1(6)*yy*xyz(1)+r1(7)*xyz(3) &
&                +r1(9)*xyz(1)+r1(13)*xz+r0(1)*xz*yy+r0(2)*xz)
      ctmp(3,5)=-(r4(8)+r3(3)*xyz(1)+r3(8)*xyz(3)+r3(18)*xyz(3)*two+r3(25)*three+r2(3)*xz &
&                +r2(9)*xz*two+r2(11)*zz*two+r2(17)*zz+r2(23)*three+r2(25)*xyz(3) &
&                +r2(27)*xyz(1)*three+r2(28)*xyz(3)*two+r1(3)*xz*xyz(3)*two+r1(4)*zz*xyz(3) &
&                +r1(6)*zz*xyz(1)+r1(7)*xyz(3)+r1(9)*xyz(1)*three+r1(10)*xyz(3)*two &
&                +r1(13)*xz+r1(14)*xz*two+r0(1)*xz*zz+r0(2)*xz+r0(3)*xz*two)
      ctmp(4,5)=-(r4(13)+r3(4)*xyz(3)+r3(10)*xyz(1)+r3(15)*xyz(2)+r3(20)*xyz(1)+r3(26) &
&                +r2(4)*xz+r2(7)*yz+r2(10)*xz+r2(11)*xy+r2(12)*xx+r2(17)*xy+r2(24) &
&                +r2(26)*xyz(3)+r2(30)*xyz(2)+r1(1)*xz*xyz(2)+r1(2)*xz*xyz(1) &
&                +r1(4)*xy*xyz(3)+r1(6)*xy*xyz(1)+r1(8)*xyz(3)+r1(12)*xyz(2)+r1(14)*yz &
&                +r0(1)*xz*xy+r0(3)*yz)
      ctmp(5,5)=-(r4(11)+r3(5)*xyz(3)+r3(8)*xyz(1)+r3(15)*xyz(3)+r3(18)*xyz(1)+r3(21) &
&                +r3(23)+r2(5)*xz+r2(7)*zz+r2(9)*xx+r2(11)*xz+r2(11)*xz+r2(17)*xz+r2(19) &
&                +r2(21)+r2(25)*xyz(1)+r2(27)*xyz(3)+r2(28)*xyz(1)+r2(30)*xyz(3)+r2(31) &
&                +r1(1)*xz*xyz(3)+r1(3)*xz*xyz(1)+r1(4)*xz*xyz(3)+r1(6)*xz*xyz(1) &
&                +r1(7)*xyz(1)+r1(9)*xyz(3)+r1(10)*xyz(1)+r1(12)*xyz(3)+r1(14)*xx+r1(14)*zz &
&                +r1(16)*two+r0(1)*xz*xz+r0(3)*xx+r0(3)*zz+r0(5))
      ctmp(6,5)=-(r4(15)+r3(9)*xyz(1)+r3(10)*xyz(3)+r3(18)*xyz(2)+r3(20)*xyz(3)+r3(24) &
&                +r2(6)*xz+r2(9)*xy+r2(10)*zz+r2(11)*yz+r2(12)*xz+r2(17)*yz+r2(22) &
&                +r2(26)*xyz(1)+r2(28)*xyz(2)+r1(2)*xz*xyz(3)+r1(3)*xz*xyz(2)+r1(4)*yz*xyz(3) &
&                +r1(6)*yz*xyz(1)+r1(8)*xyz(1)+r1(10)*xyz(2)+r1(14)*xy+r0(1)*xz*yz+r0(3)*xy)
      ctmp(1,6)=-(r4(13)+r3(4)*xyz(3)+r3(5)*xyz(2)+r3(20)*xyz(1)*two+r3(26)+r2(1)*yz &
&                +r2(10)*xz*two+r2(11)*xy*two+r2(18)*xx+r2(24)+r2(26)*xyz(3)+r2(27)*xyz(2) &
&                +r1(1)*yz*xyz(1)*two+r1(5)*xx*xyz(3)+r1(6)*xx*xyz(2)+r1(8)*xyz(3) &
&                +r1(9)*xyz(2)+r1(13)*yz+r0(1)*yz*xx+r0(2)*yz)
      ctmp(2,6)=-(r4(7)+r3(2)*xyz(3)+r3(7)*xyz(2)+r3(17)*xyz(2)*two+r3(26)*three+r2(2)*yz &
&               +r2(8)*yz*two+r2(12)*yy*two+r2(18)*yy+r2(24)*three+r2(26)*xyz(3)*three &
&               +r2(27)*xyz(2)+r2(30)*xyz(2)*two+r1(2)*yz*xyz(2)*two+r1(5)*yy*xyz(3) &
&               +r1(6)*yy*xyz(2)+r1(8)*xyz(3)*three+r1(9)*xyz(2)+r1(12)*xyz(2)*two &
&               +r1(13)*yz+r1(14)*yz*two+r0(1)*yz*yy+r0(2)*yz+r0(3)*yz*two)
      ctmp(3,6)=-(r4(9)+r3(3)*xyz(2)+r3(9)*xyz(3)+r3(19)*xyz(3)*two+r3(26)*three+r2(3)*yz &
&                +r2(9)*yz*two+r2(12)*zz*two+r2(18)*zz+r2(24)*three+r2(26)*xyz(3) &
&                +r2(27)*xyz(2)*three+r2(29)*xyz(3)*two+r1(3)*yz*xyz(3)*two+r1(5)*zz*xyz(3) &
&                +r1(6)*zz*xyz(2)+r1(8)*xyz(3)+r1(9)*xyz(2)*three+r1(11)*xyz(3)*two &
&                +r1(13)*yz+r1(14)*yz*two+r0(1)*yz*zz+r0(2)*yz+r0(3)*yz*two)
      ctmp(4,6)=-(r4(14)+r3(6)*xyz(3)+r3(10)*xyz(2)+r3(17)*xyz(1)+r3(20)*xyz(2)+r3(25) &
&                +r2(4)*yz+r2(8)*xz+r2(10)*yz+r2(11)*yy+r2(12)*xy+r2(18)*xy+r2(23) &
&                +r2(25)*xyz(3)+r2(30)*xyz(1)+r1(1)*yz*xyz(2)+r1(2)*yz*xyz(1)+r1(5)*xy*xyz(3) &
&                +r1(6)*xy*xyz(2)+r1(7)*xyz(3)+r1(12)*xyz(1)+r1(14)*xz+r0(1)*yz*xy+r0(3)*xz)
      ctmp(5,6)=-(r4(15)+r3(8)*xyz(2)+r3(10)*xyz(3)+r3(19)*xyz(1)+r3(20)*xyz(3)+r3(24) &
&                +r2(5)*yz+r2(9)*xy+r2(10)*zz+r2(11)*yz+r2(12)*xz+r2(18)*xz+r2(22) &
&                +r2(25)*xyz(2)+r2(29)*xyz(1)+r1(1)*yz*xyz(3)+r1(3)*yz*xyz(1)+r1(5)*xz*xyz(3) &
&                +r1(6)*xz*xyz(2)+r1(7)*xyz(2)+r1(11)*xyz(1)+r1(14)*xy+r0(1)*yz*xz+r0(3)*xy)
      ctmp(6,6)=-(r4(12)+r3(7)*xyz(3)+r3(9)*xyz(2)+r3(17)*xyz(3)+r3(19)*xyz(2)+r3(22) &
&                +r3(23)+r2(6)*yz+r2(8)*zz+r2(9)*yy+r2(12)*yz+r2(12)*yz+r2(18)*yz+r2(20) &
&                +r2(21)+r2(26)*xyz(2)+r2(27)*xyz(3)+r2(29)*xyz(2)+r2(30)*xyz(3)+r2(31) &
&                +r1(2)*yz*xyz(3)+r1(3)*yz*xyz(2)+r1(5)*yz*xyz(3)+r1(6)*yz*xyz(2) &
&                +r1(8)*xyz(2)+r1(9)*xyz(3)+r1(11)*xyz(2)+r1(12)*xyz(3)+r1(14)*yy+r1(14)*zz &
&                +r1(16)*two+r0(1)*yz*yz+r0(3)*yy+r0(3)*zz+r0(5))
!
      if(nbfij(1) == 6) then
        do j= 1,6
          ctmp2(j,1)= ctmp(j,1)
          ctmp2(j,2)= ctmp(j,4)*sqrt3
          ctmp2(j,3)= ctmp(j,5)*sqrt3
          ctmp2(j,4)= ctmp(j,2)
          ctmp2(j,5)= ctmp(j,6)*sqrt3
          ctmp2(j,6)= ctmp(j,3)
        enddo
      else
        do j= 1,6
          ctmp2(j,1)= ctmp(j,4)*sqrt3
          ctmp2(j,2)= ctmp(j,6)*sqrt3
          ctmp2(j,3)= ctmp(j,3)-(ctmp(j,1)+ctmp(j,2))*half
          ctmp2(j,4)= ctmp(j,5)*sqrt3
          ctmp2(j,5)=(ctmp(j,1)-ctmp(j,2))*sqrt3h
        enddo
      endif
      if(nbfij(2) == 6) then
        do i= 1,nbfij(1)
          cint1(1,i)= ctmp2(1,i)
          cint1(2,i)= ctmp2(4,i)*sqrt3
          cint1(3,i)= ctmp2(5,i)*sqrt3
          cint1(4,i)= ctmp2(2,i)
          cint1(5,i)= ctmp2(6,i)*sqrt3
          cint1(6,i)= ctmp2(3,i)
        enddo
      else
        do i= 1,nbfij(1)
          cint1(1,i)= ctmp2(4,i)*sqrt3
          cint1(2,i)= ctmp2(6,i)*sqrt3
          cint1(3,i)= ctmp2(3,i)-(ctmp2(1,i)+ctmp2(2,i))*half
          cint1(4,i)= ctmp2(5,i)*sqrt3
          cint1(5,i)=(ctmp2(1,i)-ctmp2(2,i))*sqrt3h
        enddo
      endif
!
      return
end


!-------------------------------------------------------------------------------
  subroutine ints(sint,exij,coij,coordij,nprimij,nangij,nbfij,mxprsh,threshex)
!-------------------------------------------------------------------------------
!
! Calculate overlap integrals
!
! In : exij     (Exponents of basis functions)
!      coij     (Coefficients of basis functions)
!      coordij  (x,y,z coordinates of basis functions)
!      nprimij  (Numbers of primitive functions)
!      nangij   (Degrees of angular momentum)
!      nbfij    (Numbers of basis functions)
!      mxprsh   (Size of primitive fuction array)
!      threshex (Threshold of exponential calculation)
! Out: sint(28,28) (Overlap integrals)
!
      implicit none
      integer,intent(in) :: nprimij(2), nangij(2), nbfij(2), mxprsh
      integer :: ncarti, ncartj, iprim, jprim, ii, jj, iang, jang
      integer :: ncart(0:6), ix, iy, iz, jx, jy, jz
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00, half=0.5D+00
      real(8),intent(in) :: exij(mxprsh,2), coij(mxprsh,2), coordij(3,2), threshex
      real(8),intent(out) :: sint(28,28)
      real(8) :: xyzij(3), rij, rij2, fac, exi, exi2, exj, ci, cj
      real(8) :: ex1, ex2, ex3, xyzpij(3,2), cij
      real(8) :: xyzint(3), sx(0:6,0:6), sy(0:6,0:6), sz(0:6,0:6)
      data ncart /1,3,6,10,15,21,28/
!
      ncarti= ncart(nangij(1))
      ncartj= ncart(nangij(2))
      sint(1:ncartj,1:ncarti)= zero
!
      if((nangij(1) > 6).or.(nangij(2) > 6))then
        write(*,'(" Error! This program supports up to i function in ints.")')
        call abort
      endif
!
      do ii= 1,3
        xyzij(ii)= coordij(ii,1)-coordij(ii,2)
      enddo
      rij= xyzij(1)*xyzij(1)+xyzij(2)*xyzij(2)+xyzij(3)*xyzij(3)
!
! Calculate overlap integrals for each primitive
!
      do iprim= 1,nprimij(1)
        exi= exij(iprim,1)
        ci = coij(iprim,1)
        exi2=-two*exi*exi
        do jprim= 1,nprimij(2)
          exj= exij(jprim,2)
          ex1= exi+exj
          ex2= one/ex1
          rij2=rij*exi*exj*ex2
          if(rij2 > threshex) cycle
          ex3= sqrt(ex2)
          fac= exp(-rij2)
          do ii= 1,3
            xyzpij(ii,1)=-exj*xyzij(ii)*ex2
            xyzpij(ii,2)= exi*xyzij(ii)*ex2
          enddo
          cj = coij(jprim,2)*fac
!
          do iang= 0,nangij(1)
            do jang= 0,nangij(2)
              call ghquad(xyzint,ex3,xyzpij,iang,jang)
              sx(jang,iang)= xyzint(1)*ex3
              sy(jang,iang)= xyzint(2)*ex3
              sz(jang,iang)= xyzint(3)*ex3
            enddo
          enddo
          cij= ci*cj
          ii= 0
          do ix= nangij(1),0,-1
            do iy= nangij(1)-ix,0,-1
              iz= nangij(1)-ix-iy
              ii= ii+1
              jj= 0
              do jx= nangij(2),0,-1
                do jy= nangij(2)-jx,0,-1
                  jz= nangij(2)-jx-jy
                  jj= jj+1
                  sint(jj,ii)= sint(jj,ii)+cij*sx(jx,ix)*sy(jy,iy)*sz(jz,iz)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
!
      if((nbfij(1) >= 5).or.(nbfij(2) >= 5)) then
        call nrmlz1(sint,nbfij(1),nbfij(2),ncarti)
      endif
!
      return
end


!------------------------------------------------------------------------
  subroutine calcmatdipole(dipmat,work,dipcenter,nproc,myrank,mpi_comm)
!------------------------------------------------------------------------
!
! Driver of dipole moment matrix calculation
!
! In  : dipcenter (Dipole moment center)
! Out : dipmat    (One electron Hamiltonian matrix)
!       work      (Overlap integral matrix)
!
      use modbasis, only : nao, nshell, mtype
      implicit none
      integer,intent(in) :: nproc, myrank, mpi_comm
      integer :: ish, jsh, num, maxfunc(0:6), maxbasis, maxdim
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: dipcenter(3)
      real(8),intent(out) :: dipmat((nao*(nao+1))/2*3), work((nao*(nao+1))/2*3)
      data maxfunc/1,3,6,10,15,21,28/
!
      maxbasis= maxval(mtype(1:nshell))
      maxdim= maxfunc(maxbasis)
!
      num=(nao*(nao+1))/2*3
      work(:)= zero
!
!$OMP parallel
      do ish= nshell-myrank,1,-nproc
!$OMP do
        do jsh= 1,ish
          call calcintdipole(work,dipcenter,ish,jsh,maxdim)
        enddo
!$OMP enddo
      enddo
!$OMP end parallel
!
      call para_allreducer(work,dipmat,num,mpi_comm)
!
      return
end


!----------------------------------------------------------
  subroutine calcintdipole(dipmat,dipcenter,ish,jsh,len1)
!----------------------------------------------------------
!
! Driver of dipole moment integrals (j|r|i)
!
      use modparam, only : mxprsh
      use modmolecule, only : coord
      use modbasis, only : locatom, locprim, locbf, mprim, mbf, mtype, ex, coeff, nao
      use modthresh, only : threshex
      implicit none
      integer,intent(in) :: ish, jsh, len1
      integer :: nangij(2), nprimij(2), nbfij(2), iatom, jatom
      integer :: iloc, jloc, ilocbf, jlocbf, iprim, jprim, i, j, ii, ij, maxj
      real(8),intent(in) :: dipcenter(3)
      real(8),intent(inout) :: dipmat((nao*(nao+1))/2,3)
      real(8) :: exij(mxprsh,2), coij(mxprsh,2), coordijk(3,3)
      real(8) :: dipint(len1,len1,3)
      logical :: iandj
!
      iandj=(ish == jsh)
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
      do i= 1,3
        coordijk(i,1)= coord(i,iatom)
        coordijk(i,2)= coord(i,jatom)
        coordijk(i,3)= dipcenter(i)
      enddo
      do iprim= 1,nprimij(1)
        exij(iprim,1)= ex(iloc+iprim)
        coij(iprim,1)= coeff(iloc+iprim)
      enddo
      do jprim= 1,nprimij(2)
        exij(jprim,2)= ex(jloc+jprim)
        coij(jprim,2)= coeff(jloc+jprim)
      enddo
!
      if((nangij(1) > 6).or.(nangij(2) > 6))then
        write(*,'(" Error! This program supports up to i function in calcintdipole.")')
        call exit
      endif
!
! Dipole moment integrals
!
      call intdipole(dipint,exij,coij,coordijk,nprimij,nangij,nbfij,len1,mxprsh,threshex)
!
      maxj= nbfij(2)
      do i= 1,nbfij(1)
        if(iandj) maxj= i
        ii= ilocbf+i
        ij= ii*(ii-1)/2+jlocbf
        do j= 1,maxj
          dipmat(ij+j,1)= dipint(j,i,1)
          dipmat(ij+j,2)= dipint(j,i,2)
          dipmat(ij+j,3)= dipint(j,i,3)
        enddo
      enddo
!
      return
end


!--------------------------------------------------------------------------------------------
  subroutine intdipole(dipint,exij,coij,coordijk,nprimij,nangij,nbfij,len1,mxprsh,threshex)
!--------------------------------------------------------------------------------------------
!
! Calculate dipole integrals
!
! In : exij     (Exponents of basis functions)
!      coij     (Coefficients of basis functions)
!      coordijk (x,y,z coordinates of basis functions and dipole center)
!      nprimij  (Numbers of primitive functions)
!      nangij   (Degrees of angular momentum)
!      nbfij    (Numbers of basis functions)
!      len1     (Dimension of dipint array)
!      mxprsh   (Size of primitive fuction array)
!      threshex (Threshold of exponential calculation)
! Out: dipint(len1,len1,3) (Dipole integrals)
!
      implicit none
      integer,intent(in) :: nprimij(2), nangij(2), nbfij(2), mxprsh, len1
      integer :: ncarti, ncartj, iprim, jprim, ii, jj, iang, jang, kang
      integer :: ncart(0:6), ix, iy, iz, jx, jy, jz
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00, half=0.5D+00
      real(8),intent(in) :: exij(mxprsh,2), coij(mxprsh,2), coordijk(3,3), threshex
      real(8),intent(out) :: dipint(len1,len1,3)
      real(8) :: xyzij(3), rij, rij2, fac, exi, exj, ci, cj
      real(8) :: ex1, ex2, ex3, xyzpijk(3,3), cij
      real(8) :: xyzint(3), dipx(0:6,0:6,0:1), dipy(0:6,0:6,0:1), dipz(0:6,0:6,0:1)
      real(8) :: diptmp(28,28,3)
      data ncart /1,3,6,10,15,21,28/
!
      ncarti= ncart(nangij(1))
      ncartj= ncart(nangij(2))
      diptmp(1:ncartj,1:ncarti,1:3)= zero
!
      if((nangij(1) > 6).or.(nangij(2) > 6))then
        write(*,'(" Error! This program supports up to i function in intdipole.")')
        call abort
      endif
!
      do ii= 1,3
        xyzij(ii)= coordijk(ii,1)-coordijk(ii,2)
      enddo
      rij= xyzij(1)*xyzij(1)+xyzij(2)*xyzij(2)+xyzij(3)*xyzij(3)
!
! Calculate dipole integrals for each primitive
!
      do iprim= 1,nprimij(1)
        exi= exij(iprim,1)
        ci = coij(iprim,1)
        do jprim= 1,nprimij(2)
          exj= exij(jprim,2)
          ex1= exi+exj
          ex2= one/ex1
          rij2=rij*exi*exj*ex2
          if(rij2 > threshex) cycle
          ex3= sqrt(ex2)
          fac= exp(-rij2)
          do ii= 1,3
            xyzpijk(ii,1)=-exj*xyzij(ii)*ex2
            xyzpijk(ii,2)= exi*xyzij(ii)*ex2
            xyzpijk(ii,3)=(exi*coordijk(ii,1)+exj*coordijk(ii,2))*ex2-coordijk(ii,3)
          enddo
          cj = coij(jprim,2)*fac
!
          do iang= 0,nangij(1)
            do jang= 0,nangij(2)
              do kang= 0,1
                call ghquadd(xyzint,ex3,xyzpijk,iang,jang,kang)
                dipx(jang,iang,kang)= xyzint(1)*ex3
                dipy(jang,iang,kang)= xyzint(2)*ex3
                dipz(jang,iang,kang)= xyzint(3)*ex3
              enddo
            enddo
          enddo
          cij= ci*cj
          ii= 0
          do ix= nangij(1),0,-1
            do iy= nangij(1)-ix,0,-1
              iz= nangij(1)-ix-iy
              ii= ii+1
              jj= 0
              do jx= nangij(2),0,-1
                do jy= nangij(2)-jx,0,-1
                  jz= nangij(2)-jx-jy
                  jj= jj+1
                  diptmp(jj,ii,1)= diptmp(jj,ii,1)+cij*dipx(jx,ix,1)*dipy(jy,iy,0)*dipz(jz,iz,0)
                  diptmp(jj,ii,2)= diptmp(jj,ii,2)+cij*dipx(jx,ix,0)*dipy(jy,iy,1)*dipz(jz,iz,0)
                  diptmp(jj,ii,3)= diptmp(jj,ii,3)+cij*dipx(jx,ix,0)*dipy(jy,iy,0)*dipz(jz,iz,1)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
!
      if((nbfij(1) >= 5).or.(nbfij(2) >= 5)) then
        call nrmlz1(diptmp(1,1,1),nbfij(1),nbfij(2),ncarti)
        call nrmlz1(diptmp(1,1,2),nbfij(1),nbfij(2),ncarti)
        call nrmlz1(diptmp(1,1,3),nbfij(1),nbfij(2),ncarti)
      endif
!
      do ii= 1,nbfij(1)
        do jj= 1,nbfij(2)
          dipint(jj,ii,1)= diptmp(jj,ii,1)
          dipint(jj,ii,2)= diptmp(jj,ii,2)
          dipint(jj,ii,3)= diptmp(jj,ii,3)
        enddo
      enddo
!
      return
end


!-------------------------------------------------------------------------------------------
  subroutine calcmatoctupole(dipmat,quadpmat,octpmat,work,dipcenter,nproc,myrank,mpi_comm)
!-------------------------------------------------------------------------------------------
!
! Driver of dipole, quadrupole, and octupole moment matrix calculation
!
! In  : dipcenter (Dipole moment center)
! Out : dipmat    (Dipole moment matrix)
!       quadpmat  (Quadrupole moment matrix)
!       octpmat   (Octupole moment matrix)
!       work      (Working matrix)
!
      use modbasis, only : nao, nshell, mtype
      implicit none
      integer,intent(in) :: nproc, myrank, mpi_comm
      integer :: ish, jsh, maxfunc(0:6), maxbasis, maxdim
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: dipcenter(3)
      real(8),intent(out) :: dipmat((nao*(nao+1))/2*3), quadpmat((nao*(nao+1))/2*6)
      real(8),intent(out) :: octpmat((nao*(nao+1))/2*10), work((nao*(nao+1))/2*10)
      data maxfunc/1,3,6,10,15,21,28/
!
      maxbasis= maxval(mtype(1:nshell))
      maxdim= maxfunc(maxbasis)
!
      quadpmat(:)= zero
      octpmat(:)= zero
      work(:)= zero
!
!$OMP parallel
      do ish= nshell-myrank,1,-nproc
!$OMP do
        do jsh= 1,ish
          call calcintoctupole(quadpmat,octpmat,work,dipcenter,ish,jsh,maxdim)
        enddo
!$OMP enddo
      enddo
!$OMP end parallel
!
      call para_allreducer(quadpmat,dipmat,(nao*(nao+1))/2*3,mpi_comm)
      call para_allreducer(octpmat,quadpmat,(nao*(nao+1))/2*6,mpi_comm)
      call para_allreducer(work,octpmat,(nao*(nao+1))/2*10,mpi_comm)
!
      return
end


!-----------------------------------------------------------------------------
  subroutine calcintoctupole(dipmat,quadpmat,octpmat,dipcenter,ish,jsh,len1)
!-----------------------------------------------------------------------------
!
! Driver of dipole (j|r|i), quadrupole (j|r^2|i), and octupole (j|r^3|i) 
! moment integrals
!
      use modparam, only : mxprsh
      use modmolecule, only : coord
      use modbasis, only : locatom, locprim, locbf, mprim, mbf, mtype, ex, coeff, nao
      use modthresh, only : threshex
      implicit none
      integer,intent(in) :: ish, jsh, len1
      integer :: nangij(2), nprimij(2), nbfij(2), iatom, jatom
      integer :: iloc, jloc, ilocbf, jlocbf, iprim, jprim, i, j, ii, ij, maxj, kk
      real(8),intent(in) :: dipcenter(3)
      real(8),intent(out) :: dipmat((nao*(nao+1))/2,3), quadpmat((nao*(nao+1))/2,6)
      real(8),intent(out) :: octpmat((nao*(nao+1))/2,10) 
      real(8) :: exij(mxprsh,2), coij(mxprsh,2), coordijk(3,3)
      real(8) :: dipint(len1,len1,3), quadpint(len1,len1,6), octpint(len1,len1,10)
      logical :: iandj
!
      iandj=(ish == jsh)
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
      do i= 1,3
        coordijk(i,1)= coord(i,iatom)
        coordijk(i,2)= coord(i,jatom)
        coordijk(i,3)= dipcenter(i)
      enddo
      do iprim= 1,nprimij(1)
        exij(iprim,1)= ex(iloc+iprim)
        coij(iprim,1)= coeff(iloc+iprim)
      enddo
      do jprim= 1,nprimij(2)
        exij(jprim,2)= ex(jloc+jprim)
        coij(jprim,2)= coeff(jloc+jprim)
      enddo
!
      if((nangij(1) > 6).or.(nangij(2) > 6))then
        write(*,'(" Error! This program supports up to i function in calcintoctupole.")')
        call exit
      endif
!
! Dipole moment integrals
!
      call intoctupole(dipint,quadpint,octpint,exij,coij,coordijk, &
&                      nprimij,nangij,nbfij,len1,mxprsh,threshex)
!
      maxj= nbfij(2)
      do i= 1,nbfij(1)
        if(iandj) maxj= i
        ii= ilocbf+i
        ij= ii*(ii-1)/2+jlocbf
        do j= 1,maxj
          dipmat(ij+j,1)= dipint(j,i,1)
          dipmat(ij+j,2)= dipint(j,i,2)
          dipmat(ij+j,3)= dipint(j,i,3)
        enddo
      enddo
      do kk= 1,6
        maxj= nbfij(2)
        do i= 1,nbfij(1)
          if(iandj) maxj= i
          ii= ilocbf+i
          ij= ii*(ii-1)/2+jlocbf
          do j= 1,maxj
            quadpmat(ij+j,kk)= quadpint(j,i,kk)
          enddo
        enddo
      enddo
      do kk= 1,10
        maxj= nbfij(2)
        do i= 1,nbfij(1)
          if(iandj) maxj= i
          ii= ilocbf+i
          ij= ii*(ii-1)/2+jlocbf
          do j= 1,maxj
            octpmat(ij+j,kk)= octpint(j,i,kk)
          enddo
        enddo
      enddo
!
      return
end


!-----------------------------------------------------------------------
  subroutine intoctupole(dipint,quadpint,octpint,exij,coij,coordijk, &
&                        nprimij,nangij,nbfij,len1,mxprsh,threshex)
!-----------------------------------------------------------------------
!
! Calculate dipole, quadrupole, and octupole moment integrals
!
! In : exij     (Exponents of basis functions)
!      coij     (Coefficients of basis functions)
!      coordijk (x,y,z coordinates of basis functions and dipole center)
!      nprimij  (Numbers of primitive functions)
!      nangij   (Degrees of angular momentum)
!      nbfij    (Numbers of basis functions)
!      len1     (Dimension of dipint array)
!      mxprsh   (Size of primitive fuction array)
!      threshex (Threshold of exponential calculation)
! Out: dipint(len1,len1,3) (Dipole integrals)
!      quadpint(len1,len1,6) (Quadrupole integrals)
!      octpint(len1,len1,10) (Octupole integrals)
!
      implicit none
      integer,intent(in) :: nprimij(2), nangij(2), nbfij(2), mxprsh, len1
      integer :: ncarti, ncartj, iprim, jprim, ii, jj, iang, jang, kang, kk
      integer :: ncart(0:6), ix, iy, iz, jx, jy, jz
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00, half=0.5D+00
      real(8),intent(in) :: exij(mxprsh,2), coij(mxprsh,2), coordijk(3,3), threshex
      real(8),intent(out) :: dipint(len1,len1,3), quadpint(len1,len1,6), octpint(len1,len1,10)
      real(8) :: xyzij(3), rij, rij2, fac, exi, exj, ci, cj
      real(8) :: ex1, ex2, ex3, xyzpijk(3,3), cij
      real(8) :: xyzint(3), octx(0:6,0:6,0:3), octy(0:6,0:6,0:3), octz(0:6,0:6,0:3)
      real(8) :: diptmp(28,28,3), quadptmp(28,28,6), octptmp(28,28,10)
      data ncart /1,3,6,10,15,21,28/
!
      ncarti= ncart(nangij(1))
      ncartj= ncart(nangij(2))
      diptmp(1:ncartj,1:ncarti,1:3)= zero
      quadptmp(1:ncartj,1:ncarti,1:6)= zero
      octptmp(1:ncartj,1:ncarti,1:10)= zero
!
      if((nangij(1) > 6).or.(nangij(2) > 6))then
        write(*,'(" Error! This program supports up to i function in intoctupole.")')
        call abort
      endif
!
      do ii= 1,3
        xyzij(ii)= coordijk(ii,1)-coordijk(ii,2)
      enddo
      rij= xyzij(1)*xyzij(1)+xyzij(2)*xyzij(2)+xyzij(3)*xyzij(3)
!
! Calculate dipole, quadrupole, and octupole integrals for each primitive
!
      do iprim= 1,nprimij(1)
        exi= exij(iprim,1)
        ci = coij(iprim,1)
        do jprim= 1,nprimij(2)
          exj= exij(jprim,2)
          ex1= exi+exj
          ex2= one/ex1
          rij2=rij*exi*exj*ex2
          if(rij2 > threshex) cycle
          ex3= sqrt(ex2)
          fac= exp(-rij2)
          do ii= 1,3
            xyzpijk(ii,1)=-exj*xyzij(ii)*ex2
            xyzpijk(ii,2)= exi*xyzij(ii)*ex2
            xyzpijk(ii,3)=(exi*coordijk(ii,1)+exj*coordijk(ii,2))*ex2-coordijk(ii,3)
          enddo
          cj = coij(jprim,2)*fac
!
          do iang= 0,nangij(1)
            do jang= 0,nangij(2)
              do kang= 0,3
                call ghquadd(xyzint,ex3,xyzpijk,iang,jang,kang)
                octx(jang,iang,kang)= xyzint(1)*ex3
                octy(jang,iang,kang)= xyzint(2)*ex3
                octz(jang,iang,kang)= xyzint(3)*ex3
              enddo
            enddo
          enddo
          cij= ci*cj
          ii= 0
          do ix= nangij(1),0,-1
            do iy= nangij(1)-ix,0,-1
              iz= nangij(1)-ix-iy
              ii= ii+1
              jj= 0
              do jx= nangij(2),0,-1
                do jy= nangij(2)-jx,0,-1
                  jz= nangij(2)-jx-jy
                  jj= jj+1
                  diptmp(jj,ii,1)= diptmp(jj,ii,1)+cij*octx(jx,ix,1)*octy(jy,iy,0)*octz(jz,iz,0)
                  diptmp(jj,ii,2)= diptmp(jj,ii,2)+cij*octx(jx,ix,0)*octy(jy,iy,1)*octz(jz,iz,0)
                  diptmp(jj,ii,3)= diptmp(jj,ii,3)+cij*octx(jx,ix,0)*octy(jy,iy,0)*octz(jz,iz,1)
                  quadptmp(jj,ii,1)= quadptmp(jj,ii,1)+cij*octx(jx,ix,2)*octy(jy,iy,0)*octz(jz,iz,0)
                  quadptmp(jj,ii,2)= quadptmp(jj,ii,2)+cij*octx(jx,ix,1)*octy(jy,iy,1)*octz(jz,iz,0)
                  quadptmp(jj,ii,3)= quadptmp(jj,ii,3)+cij*octx(jx,ix,1)*octy(jy,iy,0)*octz(jz,iz,1)
                  quadptmp(jj,ii,4)= quadptmp(jj,ii,4)+cij*octx(jx,ix,0)*octy(jy,iy,2)*octz(jz,iz,0)
                  quadptmp(jj,ii,5)= quadptmp(jj,ii,5)+cij*octx(jx,ix,0)*octy(jy,iy,1)*octz(jz,iz,1)
                  quadptmp(jj,ii,6)= quadptmp(jj,ii,6)+cij*octx(jx,ix,0)*octy(jy,iy,0)*octz(jz,iz,2)
                  octptmp(jj,ii, 1)= octptmp(jj,ii, 1)+cij*octx(jx,ix,3)*octy(jy,iy,0)*octz(jz,iz,0)
                  octptmp(jj,ii, 2)= octptmp(jj,ii, 2)+cij*octx(jx,ix,2)*octy(jy,iy,1)*octz(jz,iz,0)
                  octptmp(jj,ii, 3)= octptmp(jj,ii, 3)+cij*octx(jx,ix,2)*octy(jy,iy,0)*octz(jz,iz,1)
                  octptmp(jj,ii, 4)= octptmp(jj,ii, 4)+cij*octx(jx,ix,1)*octy(jy,iy,2)*octz(jz,iz,0)
                  octptmp(jj,ii, 5)= octptmp(jj,ii, 5)+cij*octx(jx,ix,1)*octy(jy,iy,1)*octz(jz,iz,1)
                  octptmp(jj,ii, 6)= octptmp(jj,ii, 6)+cij*octx(jx,ix,1)*octy(jy,iy,0)*octz(jz,iz,2)
                  octptmp(jj,ii, 7)= octptmp(jj,ii, 7)+cij*octx(jx,ix,0)*octy(jy,iy,3)*octz(jz,iz,0)
                  octptmp(jj,ii, 8)= octptmp(jj,ii, 8)+cij*octx(jx,ix,0)*octy(jy,iy,2)*octz(jz,iz,1)
                  octptmp(jj,ii, 9)= octptmp(jj,ii, 9)+cij*octx(jx,ix,0)*octy(jy,iy,1)*octz(jz,iz,2)
                  octptmp(jj,ii,10)= octptmp(jj,ii,10)+cij*octx(jx,ix,0)*octy(jy,iy,0)*octz(jz,iz,3)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
!
      if((nbfij(1) >= 5).or.(nbfij(2) >= 5)) then
        do kk= 1,3
          call nrmlz1(diptmp(1,1,kk),nbfij(1),nbfij(2),ncarti)
        enddo
        do kk= 1,6
          call nrmlz1(quadptmp(1,1,kk),nbfij(1),nbfij(2),ncarti)
        enddo
        do kk= 1,10
          call nrmlz1(octptmp(1,1,kk),nbfij(1),nbfij(2),ncarti)
        enddo
      endif
!
      do kk= 1,3
        do ii= 1,nbfij(1)
          do jj= 1,nbfij(2)
            dipint(jj,ii,kk)= diptmp(jj,ii,kk)
          enddo
        enddo
      enddo
      do kk= 1,6
        do ii= 1,nbfij(1)
          do jj= 1,nbfij(2)
            quadpint(jj,ii,kk)= quadptmp(jj,ii,kk)
          enddo
        enddo
      enddo
      do kk= 1,10
        do ii= 1,nbfij(1)
          do jj= 1,nbfij(2)
            octpint(jj,ii,kk)= octptmp(jj,ii,kk)
          enddo
        enddo
      enddo
!
      return
end


!-------------------------------------------
  subroutine nrmlz1(onei,nbfi,nbfj,ncarti)
!-------------------------------------------
!
! Normalize one-electron and overlap integrals
!
      implicit none
      integer,intent(in) :: nbfi, nbfj, ncarti
      integer :: i, j
      real(8),parameter :: half=0.5D+00, two=2.0D+00, three=3.0D+00, four=4.0D+00
      real(8),parameter :: five=5.0D+00, six=6.0D+00, eight=8.0D+00, p9=9.0D+00, ten=10.0D+00
      real(8),parameter :: twelve=12.0D+00, p15=15.0D+00, p16=16.0D+00
      real(8),parameter :: p20=20.0D+00, p24=24.0D+00, p30=30.0D+00, p32=32.0D+00
      real(8),parameter :: p40=40.0D+00, p60=60.0D+00, p90=90.0D+00, p120=1.2D+02
      real(8),parameter :: p180=1.8D+02, eighth=0.125D+00, sixteenth=6.25D-02
      real(8),parameter :: sqrt3=1.732050807568877D+00, sqrt3h=8.660254037844386D-01
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
      real(8),intent(inout) :: onei(28,28)
      real(8) :: work(28)
!
! Bra part
!
      select case(nbfj)
! D function
        case(5)
          do i= 1,ncarti
            do j= 1,6
              work(j)= onei(j,i)
            enddo
            onei(1,i)= work(2)*sqrt3
            onei(2,i)= work(5)*sqrt3
            onei(3,i)=(work(6)*two-work(1)-work(4))*half
            onei(4,i)= work(3)*sqrt3
            onei(5,i)=(work(1)-work(4))*sqrt3h
          enddo
        case(6)
          do i= 1,ncarti
            onei(2,i)= onei(2,i)*sqrt3
            onei(3,i)= onei(3,i)*sqrt3
            onei(5,i)= onei(5,i)*sqrt3
          enddo
! F function
        case(7)
          do i= 1,ncarti
            do j= 1,10
              work(j)= onei(j,i)
            enddo
            onei(1,i)=(-work(7)+three*work(2)                   )*facf1
            onei(2,i)=  work(5)                                  *facf2
            onei(3,i)=(-work(7)-work(2)+four*work(9)            )*facf3
            onei(4,i)=( two*work(10)-three*work(3)-three*work(8))*half
            onei(5,i)=(-work(1)-work(4)+four*work(6)            )*facf3
            onei(6,i)=( work(3)-work(8)                         )*facf4
            onei(7,i)=( work(1)-three*work(4)                   )*facf1
          enddo
        case(10)
          do i= 1,ncarti
            onei(2,i)= onei(2,i)*sqrt5
            onei(3,i)= onei(3,i)*sqrt5
            onei(4,i)= onei(4,i)*sqrt5
            onei(5,i)= onei(5,i)*sqrt15
            onei(6,i)= onei(6,i)*sqrt5
            onei(8,i)= onei(8,i)*sqrt5
            onei(9,i)= onei(9,i)*sqrt5
          enddo
! G function
        case(9)
          do i= 1,ncarti
            do j= 1,15
              work(j)= onei(j,i)
            enddo
            onei(1,i)=(work(2)-work(7))*facg1
            onei(2,i)=(-work(12)+work(5)*three)*facg2
            onei(3,i)=(-work(2)-work(7)+work(9)*six)*facg3
            onei(4,i)=(-work(12)*three+work(14)*four-work(5)*three)*facg4
            onei(5,i)=(work(1)*three+work(11)*three+work(15)*eight+work(4)*six &
&                     -work(6)*p24-work(13)*p24)*eighth
            onei(6,i)=(-work(3)*three+work(10)*four-work(8)*three)*facg4
            onei(7,i)=(-work(1)+work(11)+work(6)*six-work(13)*six)*facg5
            onei(8,i)=(work(3)-work(8)*three)*facg2
            onei(9,i)=(work(1)+work(11)-work(4)*six)*facg6
          enddo
        case(15)
          do i= 1,ncarti
            onei( 2,i)= onei( 2,i)*sqrt7
            onei( 3,i)= onei( 3,i)*sqrt7
            onei( 4,i)= onei( 4,i)*sqrt35third
            onei( 5,i)= onei( 5,i)*sqrt35
            onei( 6,i)= onei( 6,i)*sqrt35third
            onei( 7,i)= onei( 7,i)*sqrt7
            onei( 8,i)= onei( 8,i)*sqrt35
            onei( 9,i)= onei( 9,i)*sqrt35
            onei(10,i)= onei(10,i)*sqrt7
            onei(12,i)= onei(12,i)*sqrt7
            onei(13,i)= onei(13,i)*sqrt35third
            onei(14,i)= onei(14,i)*sqrt7
          enddo
! H function
        case(11)
          do i= 1,ncarti
            do j= 1,21
              work(j)= onei(j,i)
            enddo
            onei( 1,i)=(work(2)*five-work(7)*ten+work(16))*fach1
            onei( 2,i)=(work(5)*four-work(12)*four)*fach2
            onei( 3,i)=(-work(2)*three-work(7)*two+work(9)*p24+work(16)-work(18)*eight)*fach3
            onei( 4,i)=(-work(5)*two-work(12)*two+work(14)*four)*fach4
            onei( 5,i)=(work(2)+work(7)*two-work(9)*twelve+work(16)-work(18)*twelve &
&                      +work(20)*eight)*fach5
            onei( 6,i)=(work(3)*p15+work(8)*p30-work(10)*p40+work(17)*p15-work(19)*p40 &
&                      +work(21)*eight)*eighth
            onei( 7,i)=(work(1)+work(4)*two-work(6)*twelve+work(11)-work(13)*twelve &
&                      +work(15)*eight)*fach5
            onei( 8,i)=(-work(3)+work(10)*two+work(17)-work(19)*two)*fach4
            onei( 9,i)=(-work(1)+work(4)*two+work(6)*eight+work(11)*three-work(13)*p24)*fach3
            onei(10,i)=(work(3)-work(8)*six+work(17))*fach2
            onei(11,i)=(work(1)-work(4)*ten+work(11)*five)*fach1
          enddo
        case(21)
          do i= 1,ncarti
            onei( 2,i)= onei( 2,i)*three
            onei( 3,i)= onei( 3,i)*three
            onei( 4,i)= onei( 4,i)*sqrt21
            onei( 5,i)= onei( 5,i)*sqrt63
            onei( 6,i)= onei( 6,i)*sqrt21
            onei( 7,i)= onei( 7,i)*sqrt21
            onei( 8,i)= onei( 8,i)*sqrt105
            onei( 9,i)= onei( 9,i)*sqrt105
            onei(10,i)= onei(10,i)*sqrt21
            onei(11,i)= onei(11,i)*three
            onei(12,i)= onei(12,i)*sqrt63
            onei(13,i)= onei(13,i)*sqrt105
            onei(14,i)= onei(14,i)*sqrt63
            onei(15,i)= onei(15,i)*three
            onei(17,i)= onei(17,i)*three
            onei(18,i)= onei(18,i)*sqrt21
            onei(19,i)= onei(19,i)*sqrt21
            onei(20,i)= onei(20,i)*three
          enddo
! I function
        case(13)
          do i= 1,ncarti
            do j= 1,28
              work(j)= onei(j,i)
            enddo
            onei( 1,i)=(work(2)*six-work(7)*p20+work(16)*six)*faci1
            onei( 2,i)=(work(5)*five-work(12)*ten+work(23))*faci2
            onei( 3,i)=(-work(2)*four+work(9)*p40+work(16)*four-work(18)*p40)*faci3
            onei( 4,i)=(-work(5)*p9-work(12)*six+work(14)*p24+work(23)*three &
&                      -work(25)*eight)*faci4
            onei( 5,i)=(work(2)*two+work(7)*four-work(9)*p32+work(16)*two &
&                      -work(18)*p32+work(20)*p32)*faci5
            onei( 6,i)=(work(5)*five+work(12)*ten-work(14)*p20+work(23)*five &
&                      -work(25)*p20+work(27)*eight)*faci6
            onei( 7,i)=(-work(1)*five-work(4)*p15+work(6)*p90-work(11)*p15 &
&                      +work(13)*p180-work(15)*p120-work(22)*five+work(24)*p90 &
&                      -work(26)*p120+work(28)*p16)*sixteenth
            onei( 8,i)=(work(3)*five+work(8)*ten-work(10)*p20+work(17)*five &
&                      -work(19)*p20+work(21)*eight)*faci6
            onei( 9,i)=(work(1)+work(4)-work(6)*p16-work(11)+work(15)*p16 &
&                      -work(22)+work(24)*p16-work(26)*p16)*faci5
            onei(10,i)=(-work(3)*three+work(8)*six+work(10)*eight+work(17)*p9 &
&                      -work(19)*p24)*faci4
            onei(11,i)=(-work(1)+work(4)*five+work(6)*ten+work(11)*five &
&                      -work(13)*p60-work(22)+work(24)*ten)*faci3
            onei(12,i)=(work(3)-work(8)*ten+work(17)*five)*faci2
            onei(13,i)=(work(1)-work(4)*p15+work(11)*p15-work(22))*faci1
          enddo
        case(28)
          do i= 1,ncarti
            onei( 2,i)= onei( 2,i)*sqrt11
            onei( 3,i)= onei( 3,i)*sqrt11
            onei( 4,i)= onei( 4,i)*sqrt33
            onei( 5,i)= onei( 5,i)*sqrt99
            onei( 6,i)= onei( 6,i)*sqrt33
            onei( 7,i)= onei( 7,i)*sqrt231fifth
            onei( 8,i)= onei( 8,i)*sqrt231
            onei( 9,i)= onei( 9,i)*sqrt231
            onei(10,i)= onei(10,i)*sqrt231fifth
            onei(11,i)= onei(11,i)*sqrt33
            onei(12,i)= onei(12,i)*sqrt231
            onei(13,i)= onei(13,i)*sqrt385
            onei(14,i)= onei(14,i)*sqrt231
            onei(15,i)= onei(15,i)*sqrt33
            onei(16,i)= onei(16,i)*sqrt11
            onei(17,i)= onei(17,i)*sqrt99
            onei(18,i)= onei(18,i)*sqrt231
            onei(19,i)= onei(19,i)*sqrt231
            onei(20,i)= onei(20,i)*sqrt99
            onei(21,i)= onei(21,i)*sqrt11
            onei(23,i)= onei(23,i)*sqrt11
            onei(24,i)= onei(24,i)*sqrt33
            onei(25,i)= onei(25,i)*sqrt231fifth
            onei(26,i)= onei(26,i)*sqrt33
            onei(27,i)= onei(27,i)*sqrt11
          enddo
      end select
!
! Ket part
!
      select case(nbfi)
! D function
        case(5)
          do j= 1,nbfj
            do i= 1,6
              work(i)= onei(j,i)
            enddo
            onei(j,1)= work(2)*sqrt3
            onei(j,2)= work(5)*sqrt3
            onei(j,3)=(work(6)*two-work(1)-work(4))*half
            onei(j,4)= work(3)*sqrt3
            onei(j,5)=(work(1)-work(4))*sqrt3h
          enddo
        case(6)
          do j= 1,nbfj
            onei(j,2)= onei(j,2)*sqrt3
            onei(j,3)= onei(j,3)*sqrt3
            onei(j,5)= onei(j,5)*sqrt3
          enddo
! F function
        case(7)
          do j= 1,nbfj
            do i= 1,10
              work(i)= onei(j,i)
            enddo
            onei(j,1)=(-work(7)+three*work(2)                   )*facf1
            onei(j,2)=  work(5)                                  *facf2
            onei(j,3)=(-work(7)-work(2)+four*work(9)            )*facf3
            onei(j,4)=( two*work(10)-three*work(3)-three*work(8))*half
            onei(j,5)=(-work(1)-work(4)+four*work(6)            )*facf3
            onei(j,6)=( work(3)-work(8)                         )*facf4
            onei(j,7)=( work(1)-three*work(4)                   )*facf1
          enddo
        case(10)
          do j= 1,nbfj
            onei(j,2)= onei(j,2)*sqrt5
            onei(j,3)= onei(j,3)*sqrt5
            onei(j,4)= onei(j,4)*sqrt5
            onei(j,5)= onei(j,5)*sqrt15
            onei(j,6)= onei(j,6)*sqrt5
            onei(j,8)= onei(j,8)*sqrt5
            onei(j,9)= onei(j,9)*sqrt5
          enddo
! G function
        case(9)
          do j= 1,nbfj
            do i= 1,15
              work(i)= onei(j,i)
            enddo
            onei(j,1)=(work(2)-work(7))*facg1
            onei(j,2)=(-work(12)+work(5)*three)*facg2
            onei(j,3)=(-work(2)-work(7)+work(9)*six)*facg3
            onei(j,4)=(-work(12)*three+work(14)*four-work(5)*three)*facg4
            onei(j,5)=(work(1)*three+work(11)*three+work(15)*eight+work(4)*six &
&                     -work(6)*p24-work(13)*p24)*eighth
            onei(j,6)=(-work(3)*three+work(10)*four-work(8)*three)*facg4
            onei(j,7)=(-work(1)+work(11)+work(6)*six-work(13)*six)*facg5
            onei(j,8)=(work(3)-work(8)*three)*facg2
            onei(j,9)=(work(1)+work(11)-work(4)*six)*facg6
          enddo
        case(15)
          do j= 1,nbfj
            onei(j, 2)= onei(j, 2)*sqrt7
            onei(j, 3)= onei(j, 3)*sqrt7
            onei(j, 4)= onei(j, 4)*sqrt35third
            onei(j, 5)= onei(j, 5)*sqrt35
            onei(j, 6)= onei(j, 6)*sqrt35third
            onei(j, 7)= onei(j, 7)*sqrt7
            onei(j, 8)= onei(j, 8)*sqrt35
            onei(j, 9)= onei(j, 9)*sqrt35
            onei(j,10)= onei(j,10)*sqrt7
            onei(j,12)= onei(j,12)*sqrt7
            onei(j,13)= onei(j,13)*sqrt35third
            onei(j,14)= onei(j,14)*sqrt7
          enddo
! H function
        case(11)
          do j= 1,nbfj
            do i= 1,21
              work(i)= onei(j,i)
            enddo
            onei(j, 1)=(work(2)*five-work(7)*ten+work(16))*fach1
            onei(j, 2)=(work(5)*four-work(12)*four)*fach2
            onei(j, 3)=(-work(2)*three-work(7)*two+work(9)*p24+work(16)-work(18)*eight)*fach3
            onei(j, 4)=(-work(5)*two-work(12)*two+work(14)*four)*fach4
            onei(j, 5)=(work(2)+work(7)*two-work(9)*twelve+work(16)-work(18)*twelve &
&                      +work(20)*eight)*fach5
            onei(j, 6)=(work(3)*p15+work(8)*p30-work(10)*p40+work(17)*p15-work(19)*p40 &
&                      +work(21)*eight)*eighth
            onei(j, 7)=(work(1)+work(4)*two-work(6)*twelve+work(11)-work(13)*twelve &
&                      +work(15)*eight)*fach5
            onei(j, 8)=(-work(3)+work(10)*two+work(17)-work(19)*two)*fach4
            onei(j, 9)=(-work(1)+work(4)*two+work(6)*eight+work(11)*three-work(13)*p24)*fach3
            onei(j,10)=(work(3)-work(8)*six+work(17))*fach2
            onei(j,11)=(work(1)-work(4)*ten+work(11)*five)*fach1
          enddo
        case(21)
          do j= 1,nbfj
            onei(j, 2)= onei(j, 2)*three
            onei(j, 3)= onei(j, 3)*three
            onei(j, 4)= onei(j, 4)*sqrt21
            onei(j, 5)= onei(j, 5)*sqrt63
            onei(j, 6)= onei(j, 6)*sqrt21
            onei(j, 7)= onei(j, 7)*sqrt21
            onei(j, 8)= onei(j, 8)*sqrt105
            onei(j, 9)= onei(j, 9)*sqrt105
            onei(j,10)= onei(j,10)*sqrt21
            onei(j,11)= onei(j,11)*three
            onei(j,12)= onei(j,12)*sqrt63
            onei(j,13)= onei(j,13)*sqrt105
            onei(j,14)= onei(j,14)*sqrt63
            onei(j,15)= onei(j,15)*three
            onei(j,17)= onei(j,17)*three
            onei(j,18)= onei(j,18)*sqrt21
            onei(j,19)= onei(j,19)*sqrt21
            onei(j,20)= onei(j,20)*three
          enddo
! I function
        case(13)
          do j= 1,nbfj
            do i= 1,28
              work(i)= onei(j,i)
            enddo
            onei(j, 1)=(work(2)*six-work(7)*p20+work(16)*six)*faci1
            onei(j, 2)=(work(5)*five-work(12)*ten+work(23))*faci2
            onei(j, 3)=(-work(2)*four+work(9)*p40+work(16)*four-work(18)*p40)*faci3
            onei(j, 4)=(-work(5)*p9-work(12)*six+work(14)*p24+work(23)*three &
&                      -work(25)*eight)*faci4
            onei(j, 5)=(work(2)*two+work(7)*four-work(9)*p32+work(16)*two &
&                      -work(18)*p32+work(20)*p32)*faci5
            onei(j, 6)=(work(5)*five+work(12)*ten-work(14)*p20+work(23)*five &
&                      -work(25)*p20+work(27)*eight)*faci6
            onei(j, 7)=(-work(1)*five-work(4)*p15+work(6)*p90-work(11)*p15 &
&                      +work(13)*p180-work(15)*p120-work(22)*five+work(24)*p90 &
&                      -work(26)*p120+work(28)*p16)*sixteenth
            onei(j, 8)=(work(3)*five+work(8)*ten-work(10)*p20+work(17)*five &
&                      -work(19)*p20+work(21)*eight)*faci6
            onei(j, 9)=(work(1)+work(4)-work(6)*p16-work(11)+work(15)*p16 &
&                      -work(22)+work(24)*p16-work(26)*p16)*faci5
            onei(j,10)=(-work(3)*three+work(8)*six+work(10)*eight+work(17)*p9 &
&                      -work(19)*p24)*faci4
            onei(j,11)=(-work(1)+work(4)*five+work(6)*ten+work(11)*five &
&                      -work(13)*p60-work(22)+work(24)*ten)*faci3
            onei(j,12)=(work(3)-work(8)*ten+work(17)*five)*faci2
            onei(j,13)=(work(1)-work(4)*p15+work(11)*p15-work(22))*faci1
          enddo
        case(28)
          do j= 1,nbfj
            onei(j, 2)= onei(j, 2)*sqrt11
            onei(j, 3)= onei(j, 3)*sqrt11
            onei(j, 4)= onei(j, 4)*sqrt33
            onei(j, 5)= onei(j, 5)*sqrt99
            onei(j, 6)= onei(j, 6)*sqrt33
            onei(j, 7)= onei(j, 7)*sqrt231fifth
            onei(j, 8)= onei(j, 8)*sqrt231
            onei(j, 9)= onei(j, 9)*sqrt231
            onei(j,10)= onei(j,10)*sqrt231fifth
            onei(j,11)= onei(j,11)*sqrt33
            onei(j,12)= onei(j,12)*sqrt231
            onei(j,13)= onei(j,13)*sqrt385
            onei(j,14)= onei(j,14)*sqrt231
            onei(j,15)= onei(j,15)*sqrt33
            onei(j,16)= onei(j,16)*sqrt11
            onei(j,17)= onei(j,17)*sqrt99
            onei(j,18)= onei(j,18)*sqrt231
            onei(j,19)= onei(j,19)*sqrt231
            onei(j,20)= onei(j,20)*sqrt99
            onei(j,21)= onei(j,21)*sqrt11
            onei(j,23)= onei(j,23)*sqrt11
            onei(j,24)= onei(j,24)*sqrt33
            onei(j,25)= onei(j,25)*sqrt231fifth
            onei(j,26)= onei(j,26)*sqrt33
            onei(j,27)= onei(j,27)*sqrt11
          enddo
      end select
      return
end
