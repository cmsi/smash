!---------------------------------------------------
  subroutine grad2eri(egrad,fulldmtrx,xint,maxdim)
!---------------------------------------------------
!
! Main driver of derivatives for two-electron integrals
!
      use procpar, only : nproc, myrank, MPI_SUM, MPI_COMM_WORLD
      use basis, only : nshell, nao
      use thresh, only : cutint2
      use molecule, only : natom
      implicit none
      integer,intent(in) :: maxdim
      integer :: ish, jsh, ksh, lsh, ij, kl
      integer :: ii, kk, kstart
      integer(8) :: ncount, icount
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: fulldmtrx(nao,nao), xint(nshell*(nshell+1)/2)
      real(8),intent(inout) :: egrad(3*natom)
      real(8) :: xijkl, twoeri(maxdim**4), dtwoeri(maxdim**4,3), pdmtrx(maxdim**4), egrad2(3*natom)
!
      egrad2= zero
      ncount=(2*nshell**3+3*nshell**2+nshell)/6+myrank
!
!$OMP parallel do schedule(dynamic,1) &
!$OMP private(ish,jsh,ksh,lsh,ij,kl,xijkl,twoeri,dtwoeri,pdmtrx,ii,kk,icount,kstart) &
!$OMP reduction(+:egrad2)
      do ish= nshell,1,-1
        ii= ish*(ish-1)/2
        icount=ncount-(2*ish*ish*ish-3*ish*ish+ish)/6
        do jsh= 1,ish
          ij= ii+jsh
          kstart=mod(icount-ish*(jsh-1),nproc)+1
 kloop:   do ksh= kstart,ish,nproc
            kk= ksh*(ksh-1)/2
            do lsh= 1,ksh
              kl= kk+lsh
              if(kl.gt.ij) exit kloop
              xijkl=xint(ij)*xint(kl)
              if(xijkl.lt.cutint2) cycle
              call calcd2eri(egrad2,fulldmtrx,pdmtrx,twoeri,dtwoeri,ish,jsh,ksh,lsh,maxdim)
            enddo
          enddo kloop
        enddo
      enddo
!$OMP end parallel do
!
      do ii= 1,3*natom
        egrad(ii)= egrad(ii)+egrad2(ii)
      enddo
      return
end


!--------------------------------------------------------------------------------------
  subroutine calcd2eri(egrad2,fulldmtrx,pdmtrx,twoeri,dtwoeri,ish,jsh,ksh,lsh,maxdim)
!--------------------------------------------------------------------------------------
!
! Driver of derivatives for two-electron integrals
!
! In    : fulldmtrx (full density matrix)
!         maxdim  (maximum dimension of twoeri and dtwoeri)
!         ish,jsh,ksh,lsh (shell indices)
! Out   : pdmtrx  (products of density matrix)
!       : twoeri  (derivatives for two-electron repulsion integrals)
!         dtwoeri (derivatives for two-electron repulsion integrals)
! Inout : egrad2  (energy gradient values)
!
      use basis, only : nao
      use param, only : mxprsh
      use molecule, only : coord, natom
      use basis, only : locatom, locprim, mprim, mbf, mtype, ex, coeff, locbf
      use thresh, only : threshex
      implicit none
      integer,parameter :: nsum(0:6)=(/1,3,6,10,15,21,28/)
      integer,intent(in) :: ish, jsh, ksh, lsh, maxdim
      integer :: i, j, k, l, iatom, jatom, katom, latom, iloc, jloc, kloc,lloc
      integer :: ilocbf, jlocbf, klocbf, llocbf, ii, jj, kk, ll, ider
      integer :: nangijkl(4), nbfijkl(4), nprimijkl(4)
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, two=2.0D+00, four=4.0D+00
      real(8),parameter :: sqrtthird=0.5773502691896258D+00, sqrt3=1.732050807568877D+00
      real(8),parameter :: sqrt3h=8.660254037844386D-01
      real(8),parameter :: sqrtfifth=0.4472135954999579D+00, sqrt3fifth=0.7745966692414834D+00
      real(8),intent(in) :: fulldmtrx(nao,nao)
      real(8),intent(out) :: pdmtrx(maxdim,maxdim,maxdim,maxdim)
      real(8),intent(out) :: twoeri(maxdim,maxdim,maxdim,maxdim)
      real(8),intent(out) :: dtwoeri(maxdim,maxdim,maxdim,maxdim,3)
      real(8),intent(inout) :: egrad2(3,natom)
      real(8) :: gradtwo(3,4), xyzijkl(3,4), exijkl(mxprsh,4), coijkl(mxprsh,4), factor, work(6)
!
      iatom= locatom(ish)
      jatom= locatom(jsh)
      katom= locatom(ksh)
      latom= locatom(lsh)
!
      if((iatom == jatom).and.(iatom == katom).and.(iatom == latom)) return
      gradtwo= zero
!
      nangijkl(1)= mtype(ish)
      nangijkl(2)= mtype(jsh)
      nangijkl(3)= mtype(ksh)
      nangijkl(4)= mtype(lsh)
!
      iloc  = locprim(ish)
      jloc  = locprim(jsh)
      kloc  = locprim(ksh)
      lloc  = locprim(lsh)
      ilocbf= locbf(ish)
      jlocbf= locbf(jsh)
      klocbf= locbf(ksh)
      llocbf= locbf(lsh)
      nbfijkl(1)= mbf(ish)
      nbfijkl(2)= mbf(jsh)
      nbfijkl(3)= mbf(ksh)
      nbfijkl(4)= mbf(lsh)
      nprimijkl(1)= mprim(ish)
      nprimijkl(2)= mprim(jsh)
      nprimijkl(3)= mprim(ksh)
      nprimijkl(4)= mprim(lsh)
!
      do i= 1,3
        xyzijkl(i,1)= coord(i,iatom)
        xyzijkl(i,2)= coord(i,jatom)
        xyzijkl(i,3)= coord(i,katom)
        xyzijkl(i,4)= coord(i,latom)
      enddo
!
      do i= 1,nprimijkl(1)
        exijkl(i,1)=ex(iloc+i)
        coijkl(i,1)=coeff(iloc+i)
      enddo
      do j= 1,nprimijkl(2)
        exijkl(j,2)=ex(jloc+j)
        coijkl(j,2)=coeff(jloc+j)
      enddo
      do k= 1,nprimijkl(3)
        exijkl(k,3)=ex(kloc+k)
        coijkl(k,3)=coeff(kloc+k)
      enddo
      do l= 1,nprimijkl(4)
        exijkl(l,4)=ex(lloc+l)
        coijkl(l,4)=coeff(lloc+l)
      enddo
!
! Check ish, jsh, ksh, lsh
!
      factor= four
      factor= one
      factor= one/four
      if(ish == jsh) factor= factor*half
      if(ksh == lsh) factor= factor*half
      if((ish == ksh).and.(jsh == lsh)) factor= factor*half
!
! 4*Dij*Dkl-DilDjk-DikDjl
!
      do i= 1,nbfijkl(1)
        ii= ilocbf+i
        do j= 1,nbfijkl(2)
          jj= jlocbf+j
          do k= 1,nbfijkl(3)
            kk= klocbf+k
            do l= 1,nbfijkl(4)
              ll= llocbf+l
              pdmtrx(l,k,j,i)= factor*(four*fulldmtrx(jj,ii)*fulldmtrx(ll,kk) &
&                         -fulldmtrx(kk,ii)*fulldmtrx(ll,jj)-fulldmtrx(ll,ii)*fulldmtrx(kk,jj))
            enddo
          enddo
        enddo
      enddo
!
! Lsh derivative 
!
      nangijkl(4)= mtype(lsh)+1
      nbfijkl(4) = nsum(nangijkl(4))
      do l= 1,nprimijkl(4)
        coijkl(l,4)= two*ex(lloc+l)*coeff(lloc+l)
      enddo
!
      if((nangijkl(1)<=2).and.(nangijkl(2)<=2).and.(nangijkl(3)<=2).and.(nangijkl(4)<=2)) then
!
! Pople-Hehre and McMurchie-Davidson scheme
!
        call int2phmd(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                     mxprsh,threshex)
      else
!
! Rys quadrature scheme
!
        call int2rys(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                    mxprsh,threshex)
      endif
!
      select case(nangijkl(4))
        case(1)
          do i= 1,nbfijkl(1)
            do j= 1,nbfijkl(2)
              do k= 1,nbfijkl(3)
                dtwoeri(1,k,j,i,1)= twoeri(1,k,j,i)
                dtwoeri(1,k,j,i,2)= twoeri(2,k,j,i)
                dtwoeri(1,k,j,i,3)= twoeri(3,k,j,i)
              enddo
            enddo
          enddo
        case(2)
          do i= 1,nbfijkl(1)
            do j= 1,nbfijkl(2)
              do k= 1,nbfijkl(3)
                dtwoeri(1,k,j,i,1)= twoeri(1,k,j,i)
                dtwoeri(2,k,j,i,1)= twoeri(4,k,j,i)*sqrtthird
                dtwoeri(3,k,j,i,1)= twoeri(5,k,j,i)*sqrtthird
                dtwoeri(1,k,j,i,2)= twoeri(4,k,j,i)*sqrtthird
                dtwoeri(2,k,j,i,2)= twoeri(2,k,j,i)
                dtwoeri(3,k,j,i,2)= twoeri(6,k,j,i)*sqrtthird
                dtwoeri(1,k,j,i,3)= twoeri(5,k,j,i)*sqrtthird
                dtwoeri(2,k,j,i,3)= twoeri(6,k,j,i)*sqrtthird
                dtwoeri(3,k,j,i,3)= twoeri(3,k,j,i)
              enddo
            enddo
          enddo
        case(3)
          do i= 1,nbfijkl(1)
            do j= 1,nbfijkl(2)
              do k= 1,nbfijkl(3)
                dtwoeri(1,k,j,i,1)= twoeri( 1,k,j,i)
                dtwoeri(2,k,j,i,1)= twoeri( 6,k,j,i)*sqrtfifth
                dtwoeri(3,k,j,i,1)= twoeri( 8,k,j,i)*sqrtfifth
                dtwoeri(4,k,j,i,1)= twoeri( 4,k,j,i)*sqrt3fifth
                dtwoeri(5,k,j,i,1)= twoeri( 5,k,j,i)*sqrt3fifth
                dtwoeri(6,k,j,i,1)= twoeri(10,k,j,i)*sqrtfifth
                dtwoeri(1,k,j,i,2)= twoeri( 4,k,j,i)*sqrtfifth
                dtwoeri(2,k,j,i,2)= twoeri( 2,k,j,i)
                dtwoeri(3,k,j,i,2)= twoeri( 9,k,j,i)*sqrtfifth
                dtwoeri(4,k,j,i,2)= twoeri( 6,k,j,i)*sqrt3fifth
                dtwoeri(5,k,j,i,2)= twoeri(10,k,j,i)*sqrtfifth
                dtwoeri(6,k,j,i,2)= twoeri( 7,k,j,i)*sqrt3fifth
                dtwoeri(1,k,j,i,3)= twoeri( 5,k,j,i)*sqrtfifth
                dtwoeri(2,k,j,i,3)= twoeri( 7,k,j,i)*sqrtfifth
                dtwoeri(3,k,j,i,3)= twoeri( 3,k,j,i)
                dtwoeri(4,k,j,i,3)= twoeri(10,k,j,i)*sqrtfifth
                dtwoeri(5,k,j,i,3)= twoeri( 8,k,j,i)*sqrt3fifth
                dtwoeri(6,k,j,i,3)= twoeri( 9,k,j,i)*sqrt3fifth
              enddo
            enddo
          enddo
        case default
          write(*,'(" Error! This program supports up to d function in calcd2eri")')
          call iabort
      end select
!
      if(mtype(lsh) >= 1) then
        nangijkl(4)= mtype(lsh)-1
        nbfijkl(4) = nsum(nangijkl(4))
        do l= 1,nprimijkl(4)
          coijkl(l,4)= coeff(lloc+l)
        enddo
!
        if((nangijkl(1)<=2).and.(nangijkl(2)<=2).and.(nangijkl(3)<=2).and.(nangijkl(4)<=2)) then
          call int2phmd(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                       mxprsh,threshex)
        else
          call int2rys(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                      mxprsh,threshex)
        endif
!
        select case(nangijkl(4))
          case(0)
            do i= 1,nbfijkl(1)
              do j= 1,nbfijkl(2)
                do k= 1,nbfijkl(3)
                  dtwoeri(1,k,j,i,1)= dtwoeri(1,k,j,i,1)-twoeri(1,k,j,i)
                  dtwoeri(2,k,j,i,2)= dtwoeri(2,k,j,i,2)-twoeri(1,k,j,i)
                  dtwoeri(3,k,j,i,3)= dtwoeri(3,k,j,i,3)-twoeri(1,k,j,i)
                enddo
              enddo
            enddo
          case(1)
            do i= 1,nbfijkl(1)
              do j= 1,nbfijkl(2)
                do k= 1,nbfijkl(3)
                  dtwoeri(1,k,j,i,1)= dtwoeri(1,k,j,i,1)-twoeri(1,k,j,i)*two
                  dtwoeri(4,k,j,i,1)= dtwoeri(4,k,j,i,1)-twoeri(2,k,j,i)*sqrt3
                  dtwoeri(5,k,j,i,1)= dtwoeri(5,k,j,i,1)-twoeri(3,k,j,i)*sqrt3
                  dtwoeri(2,k,j,i,2)= dtwoeri(2,k,j,i,2)-twoeri(2,k,j,i)*two
                  dtwoeri(4,k,j,i,2)= dtwoeri(4,k,j,i,2)-twoeri(1,k,j,i)*sqrt3
                  dtwoeri(6,k,j,i,2)= dtwoeri(6,k,j,i,2)-twoeri(3,k,j,i)*sqrt3
                  dtwoeri(3,k,j,i,3)= dtwoeri(3,k,j,i,3)-twoeri(3,k,j,i)*two
                  dtwoeri(5,k,j,i,3)= dtwoeri(5,k,j,i,3)-twoeri(1,k,j,i)*sqrt3
                  dtwoeri(6,k,j,i,3)= dtwoeri(6,k,j,i,3)-twoeri(2,k,j,i)*sqrt3
                enddo
              enddo
            enddo
            if(mbf(lsh) == 5) then
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do j= 1,nbfijkl(2)
                    do k= 1,nbfijkl(3)
                      do l= 1,6
                        work(l)= dtwoeri(l,k,j,i,ider)
                      enddo
                      dtwoeri(1,k,j,i,ider)=(work(3)*two-work(1)-work(2))*half
                      dtwoeri(2,k,j,i,ider)= work(5)
                      dtwoeri(3,k,j,i,ider)= work(6)
                      dtwoeri(4,k,j,i,ider)=(work(1)-work(2))*sqrt3h
                      dtwoeri(5,k,j,i,ider)= work(4)
                    enddo
                  enddo
                enddo
              enddo
            endif
        end select
      endif
!
      nangijkl(4)= mtype(lsh)
      nbfijkl(4) = mbf(lsh)
      do l= 1,nprimijkl(4)
        coijkl(l,4)= coeff(lloc+l)
      enddo
      do i= 1,nbfijkl(1)
        do j= 1,nbfijkl(2)
          do k= 1,nbfijkl(3)
            do l= 1,nbfijkl(4)
              gradtwo(1,4)= gradtwo(1,4)+pdmtrx(l,k,j,i)*dtwoeri(l,k,j,i,1)
              gradtwo(2,4)= gradtwo(2,4)+pdmtrx(l,k,j,i)*dtwoeri(l,k,j,i,2)
              gradtwo(3,4)= gradtwo(3,4)+pdmtrx(l,k,j,i)*dtwoeri(l,k,j,i,3)
            enddo
          enddo
        enddo
      enddo
!
! Ksh derivative
!
      nangijkl(3)= mtype(ksh)+1
      nbfijkl(3) = nsum(nangijkl(3))
      do k= 1,nprimijkl(3)
        coijkl(k,3)= two*ex(kloc+k)*coeff(kloc+k)
      enddo
!
      if((nangijkl(1)<=2).and.(nangijkl(2)<=2).and.(nangijkl(3)<=2).and.(nangijkl(4)<=2)) then
        call int2phmd(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                     mxprsh,threshex)
      else
        call int2rys(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                    mxprsh,threshex)
      endif
!
      select case(nangijkl(3))
        case(1)
          do i= 1,nbfijkl(1)
            do j= 1,nbfijkl(2)
              do l= 1,nbfijkl(4)
                dtwoeri(l,1,j,i,1)= twoeri(l,1,j,i)
                dtwoeri(l,1,j,i,2)= twoeri(l,2,j,i)
                dtwoeri(l,1,j,i,3)= twoeri(l,3,j,i)
              enddo
            enddo
          enddo
        case(2)
          do i= 1,nbfijkl(1)
            do j= 1,nbfijkl(2)
              do l= 1,nbfijkl(4)
                dtwoeri(l,1,j,i,1)= twoeri(l,1,j,i)
                dtwoeri(l,2,j,i,1)= twoeri(l,4,j,i)*sqrtthird
                dtwoeri(l,3,j,i,1)= twoeri(l,5,j,i)*sqrtthird
                dtwoeri(l,1,j,i,2)= twoeri(l,4,j,i)*sqrtthird
                dtwoeri(l,2,j,i,2)= twoeri(l,2,j,i)
                dtwoeri(l,3,j,i,2)= twoeri(l,6,j,i)*sqrtthird
                dtwoeri(l,1,j,i,3)= twoeri(l,5,j,i)*sqrtthird
                dtwoeri(l,2,j,i,3)= twoeri(l,6,j,i)*sqrtthird
                dtwoeri(l,3,j,i,3)= twoeri(l,3,j,i)
              enddo
            enddo
          enddo
        case(3)
          do i= 1,nbfijkl(1)
            do j= 1,nbfijkl(2)
              do l= 1,nbfijkl(4)
                dtwoeri(l,1,j,i,1)= twoeri(l, 1,j,i)
                dtwoeri(l,2,j,i,1)= twoeri(l, 6,j,i)*sqrtfifth
                dtwoeri(l,3,j,i,1)= twoeri(l, 8,j,i)*sqrtfifth
                dtwoeri(l,4,j,i,1)= twoeri(l, 4,j,i)*sqrt3fifth
                dtwoeri(l,5,j,i,1)= twoeri(l, 5,j,i)*sqrt3fifth
                dtwoeri(l,6,j,i,1)= twoeri(l,10,j,i)*sqrtfifth
                dtwoeri(l,1,j,i,2)= twoeri(l, 4,j,i)*sqrtfifth
                dtwoeri(l,2,j,i,2)= twoeri(l, 2,j,i)
                dtwoeri(l,3,j,i,2)= twoeri(l, 9,j,i)*sqrtfifth
                dtwoeri(l,4,j,i,2)= twoeri(l, 6,j,i)*sqrt3fifth
                dtwoeri(l,5,j,i,2)= twoeri(l,10,j,i)*sqrtfifth
                dtwoeri(l,6,j,i,2)= twoeri(l, 7,j,i)*sqrt3fifth
                dtwoeri(l,1,j,i,3)= twoeri(l, 5,j,i)*sqrtfifth
                dtwoeri(l,2,j,i,3)= twoeri(l, 7,j,i)*sqrtfifth
                dtwoeri(l,3,j,i,3)= twoeri(l, 3,j,i)
                dtwoeri(l,4,j,i,3)= twoeri(l,10,j,i)*sqrtfifth
                dtwoeri(l,5,j,i,3)= twoeri(l, 8,j,i)*sqrt3fifth
                dtwoeri(l,6,j,i,3)= twoeri(l, 9,j,i)*sqrt3fifth
              enddo
            enddo
          enddo
        case default
          write(*,'(" Error! This program supports up to d function in calcd2eri")')
          call iabort
      end select
!
      if(mtype(ksh) >= 1) then
        nangijkl(3)= mtype(ksh)-1
        nbfijkl(3) = nsum(nangijkl(3))
        do k= 1,nprimijkl(3)
          coijkl(k,3)= coeff(kloc+k)
        enddo
!
        if((nangijkl(1)<=2).and.(nangijkl(2)<=2).and.(nangijkl(3)<=2).and.(nangijkl(4)<=2)) then
          call int2phmd(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                       mxprsh,threshex)
        else
          call int2rys(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                      mxprsh,threshex)
        endif
!
        select case(nangijkl(3))
          case(0)
            do i= 1,nbfijkl(1)
              do j= 1,nbfijkl(2)
                do l= 1,nbfijkl(4)
                  dtwoeri(l,1,j,i,1)= dtwoeri(l,1,j,i,1)-twoeri(l,1,j,i)
                  dtwoeri(l,2,j,i,2)= dtwoeri(l,2,j,i,2)-twoeri(l,1,j,i)
                  dtwoeri(l,3,j,i,3)= dtwoeri(l,3,j,i,3)-twoeri(l,1,j,i)
                enddo
              enddo
            enddo
          case(1)
            do i= 1,nbfijkl(1)
              do j= 1,nbfijkl(2)
                do l= 1,nbfijkl(4)
                  dtwoeri(l,1,j,i,1)= dtwoeri(l,1,j,i,1)-twoeri(l,1,j,i)*two
                  dtwoeri(l,4,j,i,1)= dtwoeri(l,4,j,i,1)-twoeri(l,2,j,i)*sqrt3
                  dtwoeri(l,5,j,i,1)= dtwoeri(l,5,j,i,1)-twoeri(l,3,j,i)*sqrt3
                  dtwoeri(l,2,j,i,2)= dtwoeri(l,2,j,i,2)-twoeri(l,2,j,i)*two
                  dtwoeri(l,4,j,i,2)= dtwoeri(l,4,j,i,2)-twoeri(l,1,j,i)*sqrt3
                  dtwoeri(l,6,j,i,2)= dtwoeri(l,6,j,i,2)-twoeri(l,3,j,i)*sqrt3
                  dtwoeri(l,3,j,i,3)= dtwoeri(l,3,j,i,3)-twoeri(l,3,j,i)*two
                  dtwoeri(l,5,j,i,3)= dtwoeri(l,5,j,i,3)-twoeri(l,1,j,i)*sqrt3
                  dtwoeri(l,6,j,i,3)= dtwoeri(l,6,j,i,3)-twoeri(l,2,j,i)*sqrt3
                enddo
              enddo
            enddo
            if(mbf(ksh) == 5) then
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do j= 1,nbfijkl(2)
                    do l= 1,nbfijkl(4)
                      do k= 1,6
                        work(k)= dtwoeri(l,k,j,i,ider)
                      enddo
                      dtwoeri(l,1,j,i,ider)=(work(3)*two-work(1)-work(2))*half
                      dtwoeri(l,2,j,i,ider)= work(5)
                      dtwoeri(l,3,j,i,ider)= work(6)
                      dtwoeri(l,4,j,i,ider)=(work(1)-work(2))*sqrt3h
                      dtwoeri(l,5,j,i,ider)= work(4)
                    enddo
                  enddo
                enddo
              enddo
            endif
        end select
      endif
!
      nangijkl(3)= mtype(ksh)
      nbfijkl(3) = mbf(ksh)
      do k= 1,nprimijkl(3)
        coijkl(k,3)= coeff(kloc+k)
      enddo
      do i= 1,nbfijkl(1)
        do j= 1,nbfijkl(2)
          do k= 1,nbfijkl(3)
            do l= 1,nbfijkl(4)
              gradtwo(1,3)= gradtwo(1,3)+pdmtrx(l,k,j,i)*dtwoeri(l,k,j,i,1)
              gradtwo(2,3)= gradtwo(2,3)+pdmtrx(l,k,j,i)*dtwoeri(l,k,j,i,2)
              gradtwo(3,3)= gradtwo(3,3)+pdmtrx(l,k,j,i)*dtwoeri(l,k,j,i,3)
            enddo
          enddo
        enddo
      enddo
!
! Jsh derivative
!
      nangijkl(2)= mtype(jsh)+1
      nbfijkl(2) = nsum(nangijkl(2))
      do j= 1,nprimijkl(2)
        coijkl(j,2)= two*ex(jloc+j)*coeff(jloc+j)
      enddo
!
      if((nangijkl(1)<=2).and.(nangijkl(2)<=2).and.(nangijkl(3)<=2).and.(nangijkl(4)<=2)) then
        call int2phmd(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                     mxprsh,threshex)
      else
        call int2rys(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                    mxprsh,threshex)
      endif
!
      select case(nangijkl(2))
        case(1)
          do i= 1,nbfijkl(1)
            do k= 1,nbfijkl(3)
              do l= 1,nbfijkl(4)
                dtwoeri(l,k,1,i,1)= twoeri(l,k,1,i)
                dtwoeri(l,k,1,i,2)= twoeri(l,k,2,i)
                dtwoeri(l,k,1,i,3)= twoeri(l,k,3,i)
              enddo
            enddo
          enddo
        case(2)
          do i= 1,nbfijkl(1)
            do k= 1,nbfijkl(3)
              do l= 1,nbfijkl(4)
                dtwoeri(l,k,1,i,1)= twoeri(l,k,1,i)
                dtwoeri(l,k,2,i,1)= twoeri(l,k,4,i)*sqrtthird
                dtwoeri(l,k,3,i,1)= twoeri(l,k,5,i)*sqrtthird
                dtwoeri(l,k,1,i,2)= twoeri(l,k,4,i)*sqrtthird
                dtwoeri(l,k,2,i,2)= twoeri(l,k,2,i)
                dtwoeri(l,k,3,i,2)= twoeri(l,k,6,i)*sqrtthird
                dtwoeri(l,k,1,i,3)= twoeri(l,k,5,i)*sqrtthird
                dtwoeri(l,k,2,i,3)= twoeri(l,k,6,i)*sqrtthird
                dtwoeri(l,k,3,i,3)= twoeri(l,k,3,i)
              enddo
            enddo
          enddo
        case(3)
          do i= 1,nbfijkl(1)
            do k= 1,nbfijkl(3)
              do l= 1,nbfijkl(4)
                dtwoeri(l,k,1,i,1)= twoeri(l,k, 1,i)
                dtwoeri(l,k,2,i,1)= twoeri(l,k, 6,i)*sqrtfifth
                dtwoeri(l,k,3,i,1)= twoeri(l,k, 8,i)*sqrtfifth
                dtwoeri(l,k,4,i,1)= twoeri(l,k, 4,i)*sqrt3fifth
                dtwoeri(l,k,5,i,1)= twoeri(l,k, 5,i)*sqrt3fifth
                dtwoeri(l,k,6,i,1)= twoeri(l,k,10,i)*sqrtfifth
                dtwoeri(l,k,1,i,2)= twoeri(l,k, 4,i)*sqrtfifth
                dtwoeri(l,k,2,i,2)= twoeri(l,k, 2,i)
                dtwoeri(l,k,3,i,2)= twoeri(l,k, 9,i)*sqrtfifth
                dtwoeri(l,k,4,i,2)= twoeri(l,k, 6,i)*sqrt3fifth
                dtwoeri(l,k,5,i,2)= twoeri(l,k,10,i)*sqrtfifth
                dtwoeri(l,k,6,i,2)= twoeri(l,k, 7,i)*sqrt3fifth
                dtwoeri(l,k,1,i,3)= twoeri(l,k, 5,i)*sqrtfifth
                dtwoeri(l,k,2,i,3)= twoeri(l,k, 7,i)*sqrtfifth
                dtwoeri(l,k,3,i,3)= twoeri(l,k, 3,i)
                dtwoeri(l,k,4,i,3)= twoeri(l,k,10,i)*sqrtfifth
                dtwoeri(l,k,5,i,3)= twoeri(l,k, 8,i)*sqrt3fifth
                dtwoeri(l,k,6,i,3)= twoeri(l,k, 9,i)*sqrt3fifth
              enddo
            enddo
          enddo
        case default
          write(*,'(" Error! This program supports up to d function in calcd2eri")')
          call iabort
      end select
!
      if(mtype(jsh) >= 1) then
        nangijkl(2)= mtype(jsh)-1
        nbfijkl(2) = nsum(nangijkl(2))
        do j= 1,nprimijkl(2)
          coijkl(j,2)= coeff(jloc+j)
        enddo
!
        if((nangijkl(1)<=2).and.(nangijkl(2)<=2).and.(nangijkl(3)<=2).and.(nangijkl(4)<=2)) then
          call int2phmd(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                       mxprsh,threshex)
        else
          call int2rys(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                      mxprsh,threshex)
        endif
!
        select case(nangijkl(2))
          case(0)
            do i= 1,nbfijkl(1)
              do k= 1,nbfijkl(3)
                do l= 1,nbfijkl(4)
                  dtwoeri(l,k,1,i,1)= dtwoeri(l,k,1,i,1)-twoeri(l,k,1,i)
                  dtwoeri(l,k,2,i,2)= dtwoeri(l,k,2,i,2)-twoeri(l,k,1,i)
                  dtwoeri(l,k,3,i,3)= dtwoeri(l,k,3,i,3)-twoeri(l,k,1,i)
                enddo
              enddo
            enddo
          case(1)
            do i= 1,nbfijkl(1)
              do k= 1,nbfijkl(3)
                do l= 1,nbfijkl(4)
                  dtwoeri(l,k,1,i,1)= dtwoeri(l,k,1,i,1)-twoeri(l,k,1,i)*two
                  dtwoeri(l,k,4,i,1)= dtwoeri(l,k,4,i,1)-twoeri(l,k,2,i)*sqrt3
                  dtwoeri(l,k,5,i,1)= dtwoeri(l,k,5,i,1)-twoeri(l,k,3,i)*sqrt3
                  dtwoeri(l,k,2,i,2)= dtwoeri(l,k,2,i,2)-twoeri(l,k,2,i)*two
                  dtwoeri(l,k,4,i,2)= dtwoeri(l,k,4,i,2)-twoeri(l,k,1,i)*sqrt3
                  dtwoeri(l,k,6,i,2)= dtwoeri(l,k,6,i,2)-twoeri(l,k,3,i)*sqrt3
                  dtwoeri(l,k,3,i,3)= dtwoeri(l,k,3,i,3)-twoeri(l,k,3,i)*two
                  dtwoeri(l,k,5,i,3)= dtwoeri(l,k,5,i,3)-twoeri(l,k,1,i)*sqrt3
                  dtwoeri(l,k,6,i,3)= dtwoeri(l,k,6,i,3)-twoeri(l,k,2,i)*sqrt3
                enddo
              enddo
            enddo
            if(mbf(jsh) == 5) then
              do ider= 1,3
                do i= 1,nbfijkl(1)
                  do k= 1,nbfijkl(3)
                    do l= 1,nbfijkl(4)
                      do j= 1,6
                        work(j)= dtwoeri(l,k,j,i,ider)
                      enddo
                      dtwoeri(l,k,1,i,ider)=(work(3)*two-work(1)-work(2))*half
                      dtwoeri(l,k,2,i,ider)= work(5)
                      dtwoeri(l,k,3,i,ider)= work(6)
                      dtwoeri(l,k,4,i,ider)=(work(1)-work(2))*sqrt3h
                      dtwoeri(l,k,5,i,ider)= work(4)
                    enddo
                  enddo
                enddo
              enddo
            endif
        end select
      endif
!
!     nangijkl(2)= mtype(jsh)
      nbfijkl(2) = mbf(jsh)
!     do j= 1,nprimijkl(2)
!       coijkl(j,2)= coeff(jloc+j)
!     enddo
      do i= 1,nbfijkl(1)
        do j= 1,nbfijkl(2)
          do k= 1,nbfijkl(3)
            do l= 1,nbfijkl(4)
              gradtwo(1,2)= gradtwo(1,2)+pdmtrx(l,k,j,i)*dtwoeri(l,k,j,i,1)
              gradtwo(2,2)= gradtwo(2,2)+pdmtrx(l,k,j,i)*dtwoeri(l,k,j,i,2)
              gradtwo(3,2)= gradtwo(3,2)+pdmtrx(l,k,j,i)*dtwoeri(l,k,j,i,3)
            enddo
          enddo
        enddo
      enddo
!
! Ish derivative
!
      gradtwo(1,1)=-gradtwo(1,2)-gradtwo(1,3)-gradtwo(1,4)
      gradtwo(2,1)=-gradtwo(2,2)-gradtwo(2,3)-gradtwo(2,4)
      gradtwo(3,1)=-gradtwo(3,2)-gradtwo(3,3)-gradtwo(3,4)
!
! Add derivative terms of 2-ERI
!
      do i= 1,3
        egrad2(i,iatom)= egrad2(i,iatom)+gradtwo(i,1)
        egrad2(i,jatom)= egrad2(i,jatom)+gradtwo(i,2)
        egrad2(i,katom)= egrad2(i,katom)+gradtwo(i,3)
        egrad2(i,latom)= egrad2(i,latom)+gradtwo(i,4)
      enddo
!
      return
end
