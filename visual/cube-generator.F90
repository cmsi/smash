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
!
! Cube file generator for SMASH
!
!------------------
  module modparam
!------------------
!
! maximum sizes for AOs (basis functions), MOs, atoms, shells
!
      implicit none
      integer,parameter :: mxatom=1000
      integer,parameter :: mxshell=5000
      integer,parameter :: mxprim=20000
      integer,parameter :: mxgridz=2000
end module


!------------------
  module modbasis
!------------------
!
      use modparam
      integer :: nshell, nao, nprim
      integer :: locprim(mxshell+1), locbf(mxshell+1), locatom(mxshell)
      integer :: mprim(mxshell),  mbf(mxshell), mtype(mxshell)
      real(8) :: ex(mxprim), coeff(mxprim)
end module


!------------------------
  program cube_generator
!------------------------
      use modparam
      use modbasis
      implicit none
      integer,parameter :: numperang=5
      integer :: ii, jj, kk, ilen, monum, maxxyz(3), minxyz(3), lenxyz(3), numxyz(3)
      integer :: natom, nmo, neleca, nelecb, multi, iatom, numatomic(mxatom)
      real(8),parameter :: space=0.20D+00, zero=0.0D+00, thresh=1.0D-15
      real(8),parameter :: toang= 0.5291772108D+00
      real(8) :: xyzmax(3), xyzmin(3), xx, yy, zz
      real(8) :: coord(3,mxatom), znuc(mxatom), charge, orbital(mxgridz)
      real(8),allocatable :: cmoa(:,:), cmob(:,:), dmtrxa(:), dmtrxb(:)
      character(len=64) :: filecheck, filecube, viewtype, moarg
      character(len=16) :: checkversion, scftype, method, runtype
!
      call getarg(1,filecheck)
      call getarg(2,filecube)
      call getarg(3,viewtype)
      ilen= len_trim(viewtype)
      call low2up(viewtype,ilen)

!
      open(unit=10,file=filecheck,form='unformatted',err=9999)
      open(unit=20,file=filecube,form='formatted')
!
      read(10) checkversion
!
      if(checkversion(1:2) =='1.') then
        write(*,'(" Error! This cube_generator supports SMASH version 2.0 or later.")')
        stop
      endif
!
      read(10,err=9998) scftype, natom, nao, nmo, nshell, nprim, neleca, nelecb,  &
&                       method, runtype, charge, multi
!
      if((viewtype == 'MO').or.(viewtype == 'MOA')) then
        call getarg(4,moarg)
        ilen= len_trim(moarg)
        call low2up(moarg,ilen)
        if(moarg == 'HOMO') then
          monum= neleca
        elseif(moarg == 'LUMO') then
          monum= neleca+1
        else
          read(moarg,*,err=9997) monum
        endif
      elseif(viewtype == 'MOB') then
        call getarg(4,moarg)
        ilen= len_trim(moarg)
        call low2up(moarg,ilen)
        if(moarg == 'HOMO') then
          monum= nelecb
        elseif(moarg == 'LUMO') then
          monum= nelecb+1
        else
          read(moarg,*,err=9997) monum
        endif
      else
        write(*,'(" Error! Requested type ",a10,"is not supported now.")')
        stop 
      endif
!
      if(monum > nmo) then
        write(*,'(" Error! Requested MO number is larger than the number of MOs.")')
        stop 
      endif
!
      if(natom > mxatom) then
        write(*,'(" Error! Natom in checkpoint file exceeds mxatom.")')
        stop
      endif
      if(nshell > mxshell) then
        write(*,'(" Error! Nshell in checkpoint file exceeds mxshell.")')
        stop
      endif
      if(nprim > mxprim) then
        write(*,'(" Error! Nprim in checkpoint file exceeds mxprim.")')
        stop
      endif
!
      read(10)
      read(10) (numatomic(ii),ii=1,natom)
      read(10)
      read(10)((coord(jj,ii),jj=1,3),ii=1,natom)
      read(10)
      read(10) (znuc(ii),ii=1,natom)
      read(10)
      read(10) (ex(ii),ii=1,nprim)
      read(10)
      read(10) (coeff(ii),ii=1,nprim)
      read(10)
      read(10) (locprim(ii),ii=1,nshell)
      read(10)
      read(10) (locbf(ii),ii=1,nshell)
      read(10)
      read(10) (locatom(ii),ii=1,nshell)
      read(10)
      read(10) (mprim(ii),ii=1,nshell)
      read(10)
      read(10) (mbf(ii),ii=1,nshell)
      read(10)
      read(10) (mtype(ii),ii=1,nshell)
!
      if((viewtype == 'MOB').and.(scftype == 'RHF')) then
        write(*,'(" Error! Scftype in checkpoint file is RHF.")')
        stop
      endif
!
      if(scftype == 'RHF') then
        allocate(cmoa(nao,nmo),dmtrxa(nao*(nao+1)/2))
        read(10)
        read(10)((cmoa(jj,ii),jj=1,nao),ii=1,nmo) 
        read(10)
        read(10) (dmtrxa(ii),ii=1,nao*(nao+1)/2)
      elseif(scftype == 'UHF') then
        allocate(cmoa(nao,nmo),dmtrxa(nao*(nao+1)/2))
        allocate(cmob(nao,nmo),dmtrxb(nao*(nao+1)/2))
        read(10)
        read(10)((cmoa(jj,ii),jj=1,nao),ii=1,nmo)
        read(10)
        read(10)((cmob(jj,ii),jj=1,nao),ii=1,nmo)
        read(10)
        read(10) (dmtrxa(ii),ii=1,nao*(nao+1)/2)
        read(10)
        read(10) (dmtrxb(ii),ii=1,nao*(nao+1)/2)
      endif
!
! Orthonormalize basis functions
!
      call bsnrmlz
!
! Calculate grid point range
!
      do ii= 1,3
        xyzmax(ii)= maxval(coord(ii,1:natom))
        xyzmin(ii)= minval(coord(ii,1:natom))
        maxxyz(ii)= nint(xyzmax(ii))+6
        minxyz(ii)= nint(xyzmin(ii))-6
        lenxyz(ii)= maxxyz(ii)-minxyz(ii)
        numxyz(ii)= lenxyz(ii)*numperang+1
      enddo
!
      if(numxyz(3) > mxgridz) then
        write(*,'(" Error! The number of z grid points exceeds mxgridz.")')
        stop
      endif
!
! Write header of cube file
!
      write(20,'(" MO=",a16)') moarg
      write(20,'(" MO coefficients")')
      write(20,'(i5,3f12.6,i5)')-natom,(dble(minxyz(ii)),ii=1,3),1
      write(20,'(i5,3f12.6)')numxyz(1),space,zero,zero
      write(20,'(i5,3f12.6)')numxyz(2),zero,space,zero
      write(20,'(i5,3f12.6)')numxyz(3),zero,zero,space
      do iatom= 1,natom
        write(20,'(i5,4f12.6)')numatomic(iatom),znuc(iatom),(coord(ii,iatom),ii=1,3)
      enddo
      write(20,'(2i5)')1,monum
!
! Write grid data
!
      if((viewtype == 'MO').or.(viewtype == 'MOA')) then
        do ii= 1,numxyz(1)
          xx= minxyz(1)+(ii-1)*space
          do jj= 1,numxyz(2)
            yy= minxyz(2)+(jj-1)*space
            do kk= 1,numxyz(3)
              zz= minxyz(3)+(kk-1)*space
              call calcorbital2(orbital(kk),cmoa,coord,xx,yy,zz,monum)
              if(abs(orbital(kk)) < thresh) orbital(kk)= zero
            enddo
            write(20,'(1p6e13.5)')(orbital(kk),kk=1,numxyz(3))
          enddo
        enddo
      elseif(viewtype == 'MOB') then
        do ii= 1,numxyz(1)
          xx= minxyz(1)+(ii-1)*space
          do jj= 1,numxyz(2)
            yy= minxyz(2)+(jj-1)*space
            do kk= 1,numxyz(3)
              zz= minxyz(3)+(kk-1)*space
              call calcorbital2(orbital(kk),cmob,coord,xx,yy,zz,monum)
              if(abs(orbital(kk)) < thresh) orbital(kk)= zero
            enddo
            write(20,'(1p6e13.5)')(orbital(kk),kk=1,numxyz(3))
          enddo
        enddo
      endif

      close(10)
      close(20)
      goto 100
!
 9999 write(*,'(" Error! Opening checkpoint file.")')
      stop
 9998 write(*,'(" Error! Reading checkpoint file.")')
      write(*,'(" Please check the integer sizes and endianness of smash", &
&               " and cube-generator.")')
      stop
 9997 write(*,'(" Error! MO number from the command line is not correct.")')
      stop
!
  100 continue
end


!-------------------------------
  subroutine low2up(line,ilen)
!-------------------------------
!
! Convert lower charcters to upper
!
      implicit none
      integer :: ilen, ii, inum
      character(len=64),intent(inout) :: line
      character(len=64) :: linecopy
      character(len=26) :: upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      character(len=26) :: lower='abcdefghijklmnopqrstuvwxyz'
!
      linecopy(1:ilen)=line(1:ilen)
      do ii= 1,ilen
        inum=(index(lower,line(ii:ii)))
        if(inum > 0) line(ii:ii)= upper(inum:inum)
      enddo
!
      return
end


!---------------------
  subroutine bsnrmlz
!---------------------
!
! Normalize basis functions
!
      use modbasis, only : nshell, mtype, ex, coeff, locprim, mprim
      implicit none
      integer :: ishell
!
!$OMP parallel do
      do ishell= 1,nshell
        call bsnor(ishell,ex,coeff,locprim,mprim,mtype)
      enddo
!$OMP end parallel do
      return
end


!--------------------------------------------------------
  subroutine bsnor(ishell,ex,coeff,locprim,mprim,mtype)
!--------------------------------------------------------
!
! Normalize basis functions
!
      use modparam, only : mxprim, mxshell
      implicit none
      integer,intent(in) :: ishell, locprim(mxshell+1), mprim(mxshell), mtype(mxshell)
      integer :: iloc, iprim, jprim, nprimi, nangi
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00
      real(8),parameter :: pi32=5.568327996831708D+00   ! pi32=pi**(3/2)
      real(8),parameter :: pii32=0.1795871221251666D+00 ! pii32=pi**(-3/2)
      real(8),intent(in) :: ex(mxprim)
      real(8),intent(inout) :: coeff(mxprim)
      real(8) :: exij, fac1, fac2
      real(8) :: factor(0:6)=(/1.0D+0,1.0D+0,3.0D+0,15.0D+0,105.0D+0,945.0D+0,10395.0D+0/)
!
! Normalize primitive functions
!
      iloc= locprim(ishell)
      nprimi= mprim(ishell)
      nangi= mtype(ishell)
      do iprim= iloc+1,iloc+nprimi
        exij= ex(iprim)*two
        fac1= pii32*exij*sqrt(exij)*(two*exij)**nangi
        coeff(iprim)= coeff(iprim)*sqrt(fac1/factor(nangi))
      enddo
!
! Normalize contracted functions
!
      fac2= zero
      do iprim= iloc+1,iloc+nprimi
        do jprim= iloc+1,iprim-1
          exij= ex(iprim)+ex(jprim)
          fac1= exij*sqrt(exij)*(two*exij)**nangi
          fac2= fac2+coeff(iprim)*coeff(jprim)*two*factor(nangi)/fac1
        enddo
        exij= ex(iprim)*two
        fac1= exij*sqrt(exij)*(two*exij)**nangi
        fac2= fac2+coeff(iprim)*coeff(iprim)*factor(nangi)/fac1
      enddo
      fac2= one/sqrt(fac2*pi32)
      do iprim= iloc+1,iloc+nprimi
        coeff(iprim)= coeff(iprim)*fac2
      enddo
      return
end


!---------------------------------------------------------------------
  subroutine calcorbital2(orbital,cmo,coord,xbohr,ybohr,zbohr,monum)
!---------------------------------------------------------------------
      use modparam
      use modbasis
      implicit none
      integer,intent(in) :: monum
      integer :: ish, ii, iatom, numprim, numbf, ibf, iprim
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, two=2.0D+00, three=3.0D+00, four=4.0D+00
      real(8),parameter :: six=6.0D+00, eight=8.0D+00, p24=24.0D+00
      real(8),parameter :: sqrt3=1.73205080756888D+00, sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrt5=2.236067977499790D+00, sqrt15=3.872983346207417D+00
      real(8),parameter :: sqrt7=2.645751311064590D+00, sqrt35=5.916079783099616D+00
      real(8),parameter :: sqrt35third=3.415650255319866D+00
      real(8),parameter :: facf1=0.36969351199675831D+00 ! 1/sqrt(10-6/sqrt(5))
      real(8),parameter :: facf2=0.86602540378443865D+00 ! 1/sqrt(4/3)
      real(8),parameter :: facf3=0.28116020334310144D+00 ! 1/sqrt(46/3-6/sqrt(5))
      real(8),parameter :: facf4=0.24065403274177409D+00 ! 1/sqrt(28-24/sqrt(5))
      real(8),parameter :: facg1=0.19440164201192295D+00 ! 1/sqrt(1336/35-8sqrt(15/7))
      real(8),parameter :: facg2=0.36969351199675831D+00 ! 1/sqrt(10-6/sqrt(5)
      real(8),parameter :: facg3=0.15721262982485929D+00 ! 1/sqrt(1774/35-8sqrt(15/7)-8sqrt(3/35))
      real(8),parameter :: facg4=0.24313189758394717D+00 ! 1/sqrt(98/5-6/sqrt(5))
      real(8),parameter :: facg5=3.20603188768051639D-02
!                                                 ! 1/sqrt(51512/35-984sqrt(5/21)-102/sqrt(105))
      real(8),parameter :: facg6=0.18742611911532351D+00 ! 1/sqrt(196/5-24/sqrt(5))
      real(8),parameter :: facg7=1.11803398874989484D+00 ! 1/sqrt(4/5)
      real(8),intent(in) :: cmo(nao,nao), coord(3,mxatom)
      real(8),intent(in) :: xbohr, ybohr, zbohr
      real(8),intent(out) :: orbital
      real(8) :: rr, xval, yval, zval, valbasis(28), work(28), expval
!
      orbital= zero
!
      do ish= 1,nshell
        iatom= locatom(ish)
        numprim= mprim(ish)
        numbf  = mbf(ish)
        ibf    = locbf(ish)
        rr=(xbohr-coord(1,iatom))**2+(ybohr-coord(2,iatom))**2+(zbohr-coord(3,iatom))**2
        xval= xbohr-coord(1,iatom)
        yval= ybohr-coord(2,iatom)
        zval= zbohr-coord(3,iatom)
        select case(numbf)
          case(1)
            valbasis(1)= zero
            do iprim= locprim(ish)+1,locprim(ish)+numprim
              valbasis(1)= valbasis(1)+exp(-ex(iprim)*rr)*coeff(iprim)
            enddo
            orbital= orbital+valbasis(1)*cmo(ibf+1,monum)
          case(3)
            valbasis(1:3)= zero
            do iprim= locprim(ish)+1,locprim(ish)+numprim
              expval= exp(-ex(iprim)*rr)*coeff(iprim)
              valbasis(1)= valbasis(1)+xval*expval
              valbasis(2)= valbasis(2)+yval*expval
              valbasis(3)= valbasis(3)+zval*expval
            enddo
            orbital= orbital+valbasis(1)*cmo(ibf+1,monum)
            orbital= orbital+valbasis(2)*cmo(ibf+2,monum)
            orbital= orbital+valbasis(3)*cmo(ibf+3,monum)
          case(5)
            work(1:6)= zero
            do iprim= locprim(ish)+1,locprim(ish)+numprim
              expval= exp(-ex(iprim)*rr)*coeff(iprim)
              work(1)= work(1)+xval*xval*expval
              work(2)= work(2)+xval*yval*expval*sqrt3
              work(3)= work(3)+xval*zval*expval*sqrt3
              work(4)= work(4)+yval*yval*expval
              work(5)= work(5)+yval*zval*expval*sqrt3
              work(6)= work(6)+zval*zval*expval
            enddo
            valbasis(1)= work(2)
            valbasis(2)= work(5)
            valbasis(3)=(work(6)*two-work(1)-work(4))*half
            valbasis(4)= work(3)
            valbasis(5)=(work(1)-work(4))*sqrt3h
            do ii= 1,5
              orbital= orbital+valbasis(ii)*cmo(ibf+ii,monum)
            enddo
          case(6)
            valbasis(1:6)= zero
            do iprim= locprim(ish)+1,locprim(ish)+numprim
              expval= exp(-ex(iprim)*rr)*coeff(iprim)
              valbasis(1)= valbasis(1)+xval*xval*expval
              valbasis(2)= valbasis(2)+xval*yval*expval*sqrt3
              valbasis(3)= valbasis(3)+xval*zval*expval*sqrt3
              valbasis(4)= valbasis(4)+yval*yval*expval
              valbasis(5)= valbasis(5)+yval*zval*expval*sqrt3
              valbasis(6)= valbasis(6)+zval*zval*expval
            enddo
            do ii= 1,6
              orbital= orbital+valbasis(ii)*cmo(ibf+ii,monum)
            enddo
          case(7)
            work(1:10)= zero
            do iprim= locprim(ish)+1,locprim(ish)+numprim
              expval= exp(-ex(iprim)*rr)*coeff(iprim)
              work( 1)= work( 1)+xval*xval*xval*expval
              work( 2)= work( 2)+xval*xval*yval*expval*sqrt5
              work( 3)= work( 3)+xval*xval*zval*expval*sqrt5
              work( 4)= work( 4)+xval*yval*yval*expval*sqrt5
              work( 5)= work( 5)+xval*yval*zval*expval*sqrt15
              work( 6)= work( 6)+xval*zval*zval*expval*sqrt5
              work( 7)= work( 7)+yval*yval*yval*expval
              work( 8)= work( 8)+yval*yval*zval*expval*sqrt5
              work( 9)= work( 9)+yval*zval*zval*expval*sqrt5
              work(10)= work(10)+zval*zval*zval*expval
            enddo
            valbasis(1)=(-work(7)+work(2)*three                   )*facf1
            valbasis(2)=  work(5)
            valbasis(3)=(-work(7)-work(2)+work(9)*four            )*facf3
            valbasis(4)=( work(10)*two-work(3)*three-work(8)*three)*facf4
            valbasis(5)=(-work(1)-work(4)+work(6)*four            )*facf3
            valbasis(6)=( work(3)-work(8)                         )*facf2
            valbasis(7)=( work(1)-work(4)*three                   )*facf1
            do ii= 1,7
              orbital= orbital+valbasis(ii)*cmo(ibf+ii,monum)
            enddo
          case(10)
            valbasis(1:10)= zero
            do iprim= locprim(ish)+1,locprim(ish)+numprim
              expval= exp(-ex(iprim)*rr)*coeff(iprim)
              valbasis( 1)= valbasis( 1)+xval*xval*xval*expval
              valbasis( 2)= valbasis( 2)+xval*xval*yval*expval*sqrt5
              valbasis( 3)= valbasis( 3)+xval*xval*zval*expval*sqrt5
              valbasis( 4)= valbasis( 4)+xval*yval*yval*expval*sqrt5
              valbasis( 5)= valbasis( 5)+xval*yval*zval*expval*sqrt15
              valbasis( 6)= valbasis( 6)+xval*zval*zval*expval*sqrt5
              valbasis( 7)= valbasis( 7)+yval*yval*yval*expval
              valbasis( 8)= valbasis( 8)+yval*yval*zval*expval*sqrt5
              valbasis( 9)= valbasis( 9)+yval*zval*zval*expval*sqrt5
              valbasis(10)= valbasis(10)+zval*zval*zval*expval
            enddo
            do ii= 1,10
              orbital= orbital+valbasis(ii)*cmo(ibf+ii,monum)
            enddo
          case(9)
            work(1:15)= zero
            do iprim= locprim(ish)+1,locprim(ish)+numprim
              expval= exp(-ex(iprim)*rr)*coeff(iprim)
              work( 1)= work( 1)+xval*xval*xval*xval*expval
              work( 2)= work( 2)+xval*xval*xval*yval*expval*sqrt7
              work( 3)= work( 3)+xval*xval*xval*zval*expval*sqrt7
              work( 4)= work( 4)+xval*xval*yval*yval*expval*sqrt35third
              work( 5)= work( 5)+xval*xval*yval*zval*expval*sqrt35
              work( 6)= work( 6)+xval*xval*zval*zval*expval*sqrt35third
              work( 7)= work( 7)+xval*yval*yval*yval*expval*sqrt7
              work( 8)= work( 8)+xval*yval*yval*zval*expval*sqrt35
              work( 9)= work( 9)+xval*yval*zval*zval*expval*sqrt35
              work(10)= work(10)+xval*zval*zval*zval*expval*sqrt7
              work(11)= work(11)+yval*yval*yval*yval*expval
              work(12)= work(12)+yval*yval*yval*zval*expval*sqrt7
              work(13)= work(13)+yval*yval*zval*zval*expval*sqrt35third
              work(14)= work(14)+yval*zval*zval*zval*expval*sqrt7
              work(15)= work(15)+zval*zval*zval*zval*expval
              valbasis(1)=(work(2)-work(7))*facg7
              valbasis(2)=(-work(12)+work(5)*three)*facg2
              valbasis(3)=(-work(2)-work(7)+work(9)*six)*facg6
              valbasis(4)=(-work(12)*three+work(14)*four-work(5)*three)*facg4
              valbasis(5)=(work(1)*three+work(11)*three+work(15)*eight+work(4)*six &
&                         -work(6)*p24-work(13)*p24)*facg5
              valbasis(6)=(-work(3)*three+work(10)*four-work(8)*three)*facg4
              valbasis(7)=(-work(1)+work(11)+work(6)*six-work(13)*six)*facg3
              valbasis(8)=(work(3)-work(8)*three)*facg2
              valbasis(9)=(work(1)+work(11)-work(4)*six)*facg1
            enddo
            do ii= 1,9
              orbital= orbital+valbasis(ii)*cmo(ibf+ii,monum)
            enddo
          case(15)
            valbasis(1:15)= zero
            do iprim= locprim(ish)+1,locprim(ish)+numprim
              expval= exp(-ex(iprim)*rr)*coeff(iprim)
              valbasis( 1)= valbasis( 1)+xval*xval*xval*xval*expval
              valbasis( 2)= valbasis( 2)+xval*xval*xval*yval*expval*sqrt7
              valbasis( 3)= valbasis( 3)+xval*xval*xval*zval*expval*sqrt7
              valbasis( 4)= valbasis( 4)+xval*xval*yval*yval*expval*sqrt35third
              valbasis( 5)= valbasis( 5)+xval*xval*yval*zval*expval*sqrt35
              valbasis( 6)= valbasis( 6)+xval*xval*zval*zval*expval*sqrt35third
              valbasis( 7)= valbasis( 7)+xval*yval*yval*yval*expval*sqrt7
              valbasis( 8)= valbasis( 8)+xval*yval*yval*zval*expval*sqrt35
              valbasis( 9)= valbasis( 9)+xval*yval*zval*zval*expval*sqrt35
              valbasis(10)= valbasis(10)+xval*zval*zval*zval*expval*sqrt7
              valbasis(11)= valbasis(11)+yval*yval*yval*yval*expval
              valbasis(12)= valbasis(12)+yval*yval*yval*zval*expval*sqrt7
              valbasis(13)= valbasis(13)+yval*yval*zval*zval*expval*sqrt35third
              valbasis(14)= valbasis(14)+yval*zval*zval*zval*expval*sqrt7
              valbasis(15)= valbasis(15)+zval*zval*zval*zval*expval
            enddo
            do ii= 1,15
              orbital= orbital+valbasis(ii)*cmo(ibf+ii,monum)
            enddo
          case default
            write(*,'(" Error! This program supports up to g function.")')
            stop
        end select
      enddo
!
      return
end
