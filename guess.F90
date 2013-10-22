!-----------------------------------
  subroutine guessmo(cmoa,overinv)
!-----------------------------------
!
! Initial guess generation
!   iguess = 1 : extended huckel MOs
!   iguess = 2 : MOs from previous ones
!
! In  : overinv (overlap integral inverse matrix)
! Out : cmoa (initial guess orbitals)
!
      use guess, only : iguess
      use basis, only : nao
      use iofile, only : iout
      implicit none
      real(8),intent(inout):: cmoa(nao*nao), overinv(nao*nao)
!
      if(iguess == 1) then
        call huckelguess(cmoa,overinv)
      elseif(iguess == 2) then
        call updatemo(cmoa,overinv)
      else
        write(iout,'(" Error! This program supports only iguess= 1 or 2 in guessmo.")')
        call iabort
      endif
      return
end


!---------------------------------------
  subroutine huckelguess(cmoa,overinv)
!---------------------------------------
!
! Initial guess calculation
!
! In  : overinv (overlap integral inverse matrix)
! Out : cmoa (initial guess orbitals)
!
      use guess, only : nao_g, spher_g, coord_g
      use basis, only : nao
      use molecule, only : coord, natom
      implicit none
      integer :: i, j, nao_g2
      real(8),intent(in):: overinv(nao*nao)
      real(8),intent(out):: cmoa(nao*nao)
      real(8),allocatable :: hmo(:), overlap(:), work1(:), work2(:), eigen(:)
!
! Set basis functions
!
      spher_g=.true.
      call setbasis_g
!
! Set coordinate
!
      do i= 1,natom
        do j= 1,3
          coord_g(j,i)=coord(j,i)
        enddo
      enddo
!
! Set required arrays
!
      nao_g2= nao_g*nao_g
      call memset(3*nao_g2+nao*nao_g+nao_g)
      allocate(hmo(nao_g2),overlap(nao*nao_g),work1(nao_g2),work2(nao_g2),eigen(nao_g))
!
! Calculate Extended Huckel method
!
      call calchuckelg(hmo,work1,work2,eigen)
!
! Calculate overlap integrals between input basis and Huckel basis
!
      call calcover2(overlap)
!
! Project orbitals from Huckel to SCF
!
      call projectmo(cmoa,overinv,overlap,hmo,work1,work2,eigen)
      deallocate(hmo,overlap,work1,work2,eigen)
      call memunset(3*nao_g2+nao*nao_g+nao_g)
      return
end

!-------------------------------------------------
  subroutine calchuckelg(hmo,huckel,ortho,eigen)
!-------------------------------------------------
!
! Driver of extended Huckel calculation for guess generation
!
! Out : hmo (extended Huckel orbitals)
!
      use guess, only : nao_g, nmo_g
      implicit none
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(out) :: hmo(nao_g*nao_g), huckel(nao_g*nao_g)
      real(8),intent(out) :: ortho(nao_g*nao_g), eigen(nao_g)
!
! Calculate overlap integrals
! (guess basis)x(guess basis)
!
      call calcover1(hmo)
!
! Set ionization potentials
!
      call huckelip(eigen)
!
! Form extended Huckel matrix
!
      call formhuckel(huckel,hmo,eigen)
!
! Calculate canonicalization matrix
!
      call mtrxcanon(ortho,hmo,eigen,nao_g,nmo_g)
!
! Canonicalize extended Huckel matrix
!
!ishimura
!     call canonicalize(huckel,ortho,hmo,nao_g,nmo_g)
!
! Diagonalize canonicalized matrix
!
      call diag('V','U',nmo_g,huckel,nao_g,eigen)
!
! Backtransform to AO basis
!
!ishimura
call dcopy(nmo_g*nao_g,huckel,1,hmo,1)
!     call dgemm('N','N',nao_g,nmo_g,nmo_g,one,ortho,nao_g,huckel,nao_g,zero,hmo,nao_g)
      return
end


!----------------------------------------------------------------------
  subroutine projectmo(cmoa,overinv,overlap,hmo,work1,work2,eigen)
!----------------------------------------------------------------------
!
! Project orbitals from Huckel to SCF
!    C1= S11^-1 * S12 * C2 [C2t * S12t * S11^-1 * S12 * C2]^-1/2
!
! In  :  overinv (overlap integral inverse matrix of SCF basis set)
! Inout: overlap (overlap integral of guess and SCF basis sets)
!        hmo (extended Huckel orbitals)
! Out :  cmoa (initial guess orbitals)
!
      use guess, only : nao_g, nmo_g
      use basis, only : nao
      implicit none
      integer :: i, j
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: overinv(nao,nao)
      real(8),intent(inout) :: overlap(nao,nao_g), hmo(nao_g,nao_g)
      real(8),intent(out) :: cmoa(nao,nao), work1(nao_g,nao_g), work2(nao_g,nao_g)
      real(8),intent(out) :: eigen(nao_g)
      real(8) :: eigeninv
!
! Calculate S12*C2
!
      call dgemm('N','N',nao,nmo_g,nao_g,one,overlap,nao,hmo,nao_g,zero,cmoa,nao)
!
! Calculate S11-1*S12*C2
!
      call dsymm('L','U',nao,nmo_g,one,overinv,nao,cmoa,nao,zero,overlap,nao)
!
! Calculate C2t*S12t*S11-1*S12*C2
!
      call dgemm('T','N',nmo_g,nmo_g,nao,one,cmoa,nao,overlap,nao,zero,work1,nao_g)
!
! Calculate (C2t*S12t*S11-1*S12*C2)-1/2
!
      call diag('V','U',nmo_g,work1,nao_g,eigen)
!$OMP parallel do private(eigeninv)
      do i= 1,nmo_g
        eigeninv= one/sqrt(eigen(i))
        do j= 1,nmo_g
          work2(j,i)= work1(j,i)*eigeninv
        enddo
      enddo
!$OMP end parallel do
      call dgemm('N','T',nmo_g,nmo_g,nmo_g,one,work1,nao_g,work2,nao_g,zero,hmo,nao_g)
!
! Calculate C1
!
      call dgemm('N','N',nao,nmo_g,nmo_g,one,overlap,nao,hmo,nao_g,zero,cmoa,nao)
      return
end


!------------------------------
  subroutine huckelip(energy)
!------------------------------
!
! Set ionization potentials
!
! Out : energy (ionization potential)
!
      use molecule, only : natom, numatomic
      use guess, only : nao_g
      use iofile, only : iout
      implicit none
      integer :: iao, iatom, i
      real(8),intent(out) :: energy(nao_g)
      real(8) :: row1(2)=(/-5.000D-01,-9.180D-01/)
      real(8) :: row2c(3:10)
      real(8) :: row2v(2,3:10)
      real(8) :: row3c(3,11:18)
      real(8) :: row3v(2,11:18)
      data row2c/-2.48D+00,-4.73D+00,-7.70D+00,-1.13D+01,-1.56D+01,-2.07D+01,-2.64D+01,-3.28D+01/
      data row2v/-1.960D-01, 0.000D+00,-3.090D-01, 0.000D+00,-4.950D-01,-3.100D-01,&
&                -7.060D-01,-4.330D-01,-9.450D-01,-5.680D-01,-1.244D+00,-6.320D-01,&
&                -1.573D+00,-7.300D-01,-1.930D+00,-8.500D-01/
      data row3c/-4.05D+01,-2.80D+00,-1.52D+00,-4.90D+01,-3.77D+00,-2.28D+00,-5.85D+01,-4.91D+00,&
&                -3.22D+00,-6.88D+01,-6.16D+00,-4.26D+00,-8.00D+01,-7.51D+00,-5.40D+00,-9.20D+01,&
&                -9.00D+00,-6.68D+00,-1.04D+02,-1.06D+01,-8.07D+00,-1.186D+02,-1.23D+01,-9.57D+00/
      data row3v/-1.820D-01, 0.000D+00,-2.530D-01, 0.000D+00,-3.930D-01,-2.100D-01,&
&                -5.400D-01,-2.970D-01,-6.960D-01,-3.920D-01,-8.800D-01,-4.370D-01,&
&                -1.073D+00,-5.060D-01,-1.278D+00,-5.910D-01/
!
      iao= 0
!
! Set valence ionization potentials
!
      do iatom= 1,natom
        select case(numatomic(iatom))
! H - He
          case(1:2)
            iao= iao+1
            energy(iao)= row1(numatomic(iatom))
! Li - Ne
          case(3:10)
            iao= iao+1
            energy(iao)= row2v(1,numatomic(iatom))
            do i= 1,3
              iao= iao+1
              energy(iao)= row2v(2,numatomic(iatom))
            enddo
! Na - Ar
          case(11:18)
            iao= iao+1
            energy(iao)= row3v(1,numatomic(iatom))
            do i= 1,3
              iao= iao+1
              energy(iao)= row3v(2,numatomic(iatom))
            enddo
          case default
            write(iout,'(" Error! This program supports up to Ar in huckelip.")')
            call iabort
        end select
      enddo
!
! Set core ionization potentials
!
      do iatom= 1,natom
        select case(numatomic(iatom))
! Li - Ne
          case(1:2)
          case(3:10)
            iao= iao+1
            energy(iao)= row2c(numatomic(iatom))
! Na - Ar
          case(11:18)
            iao= iao+1
            energy(iao)= row3c(1,numatomic(iatom))
            iao= iao+1
            energy(iao)= row3c(2,numatomic(iatom))
            do i= 1,3
              iao= iao+1
              energy(iao)= row3c(3,numatomic(iatom))
            enddo
          case default
            write(iout,'(" Error! This program supports up to Ar in huckelip.")')
            call iabort
        end select
      enddo
      return
end


!-----------------------------------------------
  subroutine formhuckel(huckel,overlap,energy)
!-----------------------------------------------
!
! Form extended Huckel matrix
!
! In  : overlap (overlap integral of guess basis set)
!       energy (ionization potential)
! Out : huckel (extended Huckel Hamiltonian)
!
      use guess, only : nao_g, nao_v
      implicit none
      integer :: i, j
      real(8),parameter :: factor=0.875D+00  !(=1.75/2.0)
      real(8),parameter :: fdown=0.05D+00
      real(8),intent(in) :: overlap(nao_g,nao_g), energy(nao_g)
      real(8),intent(out) :: huckel(nao_g,nao_g)
!
! Generate Huckel matrix from overlap integrals and ionization potentials
!
!$OMP parallel do schedule(static,1) private(j)
      do i= 1,nao_g
        do j= 1,i-1
          huckel(j,i)=factor*overlap(j,i)*(energy(i)+energy(j))
        enddo
        huckel(i,i)= energy(i)
      enddo
!$OMP end parallel do
!
! Scale down core-core and core-valence elements
!
!$OMP parallel do schedule(static,1) private(j)
      do i= nao_v+1,nao_g
        do j= 1,i-1
          huckel(j,i)= fdown*huckel(j,i)
        enddo
      enddo
!$OMP end parallel do
      return
end


!--------------------------------
  subroutine calcover1(overlap)
!--------------------------------
!
! Driver of overlap integral calculation
! (guess basis)x(guess basis)
!
! Out : overlap (overlap integral of guess basis set)
!
      use guess, only : nshell_g, nao_g
      implicit none
      integer :: ish, jsh
      real(8),intent(out) :: overlap(nao_g*nao_g)
!
!$OMP parallel do private(jsh)
      do ish= nshell_g,1,-1
        do jsh= 1,ish
          call intover1(overlap,ish,jsh)
        enddo
      enddo
!$OMP end parallel do
      return
end


!--------------------------------
  subroutine calcover2(overlap)
!--------------------------------
!
! Driver of overlap integral calculation
! (input basis)x(guess basis)
!
! Out : overlap (overlap integral of guess and SCF basis sets)
!
      use procpar, only : nproc, myrank, MPI_SUM, MPI_COMM_WORLD
      use basis, only : nshell, nao
      use guess, only : nshell_g, nao_g
      implicit none
      integer :: ish, jsh
      real(8),intent(out) :: overlap(nao,nao_g)
      real(8) :: dum
!
      call zeroclr(overlap,nao*nao_g)
!$OMP parallel
      do ish= nshell_g-myrank,1,-nproc
!$OMP do
        do jsh= 1,nshell
          call intover2(overlap,ish,jsh)
        enddo
!$OMP enddo
      enddo
!$OMP end parallel
!
      call para_allreduce(overlap,dum,nao*nao_g,"D",MPI_SUM,MPI_COMM_WORLD,1)
      return
end


!-----------------------------------------
  subroutine intover1(overlap,ish,jsh)
!-----------------------------------------
!
! Overlap integral calculation
! (guess basis)x(guess basis)
!
! In  : ish, jsh (shell index)
! Out : overlap (overlap integral of guess basis set)
!
      use param, only : mxprsh
      use thresh, only : threshex
      use basis, only : locatom, locprim, locbf, mprim, mtype, ex, coeff
      use guess, only : locatom_g, locprim_g, locbf_g, mprim_g, mbf_g, mtype_g, &
&                       ex_g, coeff_g, nao_g, coord_g
      use hermite, only : ix, iy, iz
      use iofile, only : iout
      implicit none
      integer,intent(in) :: ish, jsh
      integer :: iatom, jatom, iloc, jloc, ilocbf, jlocbf, nprimi, nprimj, nangi, nangj 
      integer :: nbfi, nbfj, iprim, jprim, nsumi, nsumj, i, j, iang, jang
      integer :: isx, jsx, isy, jsy, isz, jsz, maxj
      integer :: nsum(0:5)=(/1,3,6,10,15,21/)
      real(8),parameter :: zero=0.0D+0, one=1.0D+0, sqrt3=1.732050807568877D+0
      real(8),parameter :: sqrt5=2.236067977499790D+0, sqrt15=3.872983346207417D+0
      real(8),intent(out) :: overlap(nao_g,nao_g)
      real(8) :: xyzij(3), rij, rij2, fac, exi, exj, ci, cj, ex1, ex2, ex3, xyzpij(3,2)
      real(8) :: xyzint(3), sx(0:4,0:4), sy(0:4,0:4), sz(0:4,0:4), sint(28,28)
      logical :: iandj
!
! Set parameters
!
      iandj=(ish == jsh)
      iatom= locatom_g(ish)
      iloc= locprim_g(ish)
      ilocbf= locbf_g(ish)
      nprimi= mprim_g(ish)
      nangi= mtype_g(ish)
      nbfi= mbf_g(ish)
      nsumi= nsum(nangi)
      jatom= locatom_g(jsh)
      jloc= locprim_g(jsh)
      jlocbf= locbf_g(jsh)
      nprimj= mprim_g(jsh)
      nangj= mtype_g(jsh)
      nbfj= mbf_g(jsh)
      nsumj= nsum(nangj)
      do i= 1,3
        xyzij(i)= coord_g(i,iatom)-coord_g(i,jatom)
      enddo
      rij= xyzij(1)*xyzij(1)+xyzij(2)*xyzij(2)+xyzij(3)*xyzij(3)
      do i= 1,nsumi
        do j= 1,nsumj
          sint(j,i)=zero
        enddo
      enddo
!
! Calculate overlap integrals for each primitive
!
      do iprim= 1,nprimi
        exi= ex_g(iloc+iprim)
        ci = coeff_g(iloc+iprim)
        do jprim= 1,nprimj
          exj= ex_g(jloc+jprim)
          ex1= exi+exj
          ex2= one/ex1
          rij2=rij*exi*exj*ex2
          if(rij2 > threshex) cycle
          ex3= sqrt(ex2)
          fac= exp(-rij2)   !*ex2*ex3
          do i= 1,3
            xyzpij(i,1)=-exj*xyzij(i)*ex2
            xyzpij(i,2)= exi*xyzij(i)*ex2
          enddo
          cj = coeff_g(jloc+jprim)*fac
          do iang= 0,nangi
            do jang= 0,nangj
              call ghquad(xyzint, ex3, xyzpij, iang, jang)
              sx(jang,iang)= xyzint(1)*ex3
              sy(jang,iang)= xyzint(2)*ex3
              sz(jang,iang)= xyzint(3)*ex3
            enddo
          enddo
          do i= 1,nsumi
            isx= ix(i,nangi)
            isy= iy(i,nangi)
            isz= iz(i,nangi)
            do j= 1,nsumj
              jsx= ix(j,nangj)
              jsy= iy(j,nangj)
              jsz= iz(j,nangj)
              sint(j,i)= sint(j,i)+ci*cj*sx(jsx,isx)*sy(jsy,isy)*sz(jsz,isz)
            enddo
          enddo
        enddo
      enddo

      if((nbfi >= 5).or.(nbfj >= 5)) then
        call nrmlz1(sint,nbfi,nbfj,nsumi)
      endif

      maxj= nbfj
      do i= 1,nbfi
        if(iandj) maxj= i
        do j= 1,maxj
          overlap(jlocbf+j,ilocbf+i)= sint(j,i)
        enddo
      enddo
      return
end


!-----------------------------------------
  subroutine intover2(overlap,ish,jsh)
!-----------------------------------------
!
! Overlap integral calculation
! (input basis)x(guess basis)
!
! In  : ish, jsh (shell index)
! Out : overlap (overlap integral of guess and SCF basis sets)
!
      use param, only : mxprsh
      use thresh, only : threshex
      use molecule, only : coord
      use basis, only : locatom, locprim, locbf, mprim, mbf, mtype, ex, coeff, nao
      use guess, only : locatom_g, locprim_g, locbf_g, mprim_g, mbf_g, mtype_g, &
&                       ex_g, coeff_g, nao_g, coord_g
      use hermite, only : ix, iy, iz
      use iofile, only : iout
      implicit none
      integer,intent(in) :: ish, jsh
      integer :: iatom, jatom, iloc, jloc, ilocbf, jlocbf, nprimi, nprimj, nangi, nangj 
      integer :: nbfi, nbfj, iprim, jprim, nsumi, nsumj, i, j, iang, jang
      integer :: isx, jsx, isy, jsy, isz, jsz
      integer :: nsum(0:5)=(/1,3,6,10,15,21/)
      real(8),parameter :: zero=0.0D+0, one=1.0D+0, sqrt3=1.732050807568877D+0
      real(8),parameter :: sqrt5=2.236067977499790D+0, sqrt15=3.872983346207417D+0
      real(8),intent(out) :: overlap(nao,nao_g)
      real(8) :: xyzij(3), rij, rij2, fac, exi, exj, ci, cj, ex1, ex2, ex3, xyzpij(3,2)
      real(8) :: xyzint(3), sx(0:5,0:5), sy(0:5,0:5), sz(0:5,0:5),sint(28,28)
!
! Set parameters
!
      iatom= locatom_g(ish)
      iloc= locprim_g(ish)
      ilocbf= locbf_g(ish)
      nprimi= mprim_g(ish)
      nangi= mtype_g(ish)
      nbfi= mbf_g(ish)
      nsumi= nsum(nangi)
      jatom= locatom(jsh)
      jloc= locprim(jsh)
      jlocbf= locbf(jsh)
      nprimj= mprim(jsh)
      nangj= mtype(jsh)
      nbfj= mbf(jsh)
      nsumj= nsum(nangj)
      do i= 1,3
        xyzij(i)= coord(i,iatom)-coord_g(i,jatom)
      enddo
      rij= xyzij(1)*xyzij(1)+xyzij(2)*xyzij(2)+xyzij(3)*xyzij(3)
      do i= 1,nsumi
        do j= 1,nsumj
          sint(j,i)=zero
        enddo
      enddo
!
! Calculate overlap integrals for each primitive
!
      do iprim= 1,nprimi
        exi= ex_g(iloc+iprim)
        ci = coeff_g(iloc+iprim)
        do jprim= 1,nprimj
          exj= ex(jloc+jprim)
          ex1= exi+exj
          ex2= one/ex1
          rij2=rij*exi*exj*ex2
          if(rij2 > threshex) cycle
          ex3= sqrt(ex2)
          fac= exp(-rij2)   !*ex2*ex3
          do i= 1,3
            xyzpij(i,1)=-exj*xyzij(i)*ex2
            xyzpij(i,2)= exi*xyzij(i)*ex2
          enddo
          cj = coeff(jloc+jprim)*fac
          do iang= 0,nangi
            do jang= 0,nangj
              call ghquad(xyzint, ex3, xyzpij, iang, jang)
              sx(jang,iang)= xyzint(1)*ex3
              sy(jang,iang)= xyzint(2)*ex3
              sz(jang,iang)= xyzint(3)*ex3
            enddo
          enddo
          do i= 1,nsumi
            isx= ix(i,nangi)
            isy= iy(i,nangi)
            isz= iz(i,nangi)
            do j= 1,nsumj
              jsx= ix(j,nangj)
              jsy= iy(j,nangj)
              jsz= iz(j,nangj)
              sint(j,i)= sint(j,i)+ci*cj*sx(jsx,isx)*sy(jsy,isy)*sz(jsz,isz)
            enddo
          enddo
        enddo
      enddo

      if((nbfi >= 5).or.(nbfj >= 5)) then
        call nrmlz1(sint,nbfi,nbfj,nsumi)
      endif

      do i= 1,nbfi
        do j= 1,nbfj
          overlap(jlocbf+j,ilocbf+i)= sint(j,i)
        enddo
      enddo
      return
end


!------------------------------------
  subroutine updatemo(cmoa,overinv)
!------------------------------------
!
! Read and project MOs
!
! Inout : overinv (overlap integral inverse matrix)
!         cmoa (initial guess orbitals)
!
      use guess, only : nao_g
      use basis, only : nao
      use molecule, only : neleca, nmo
      implicit none
      integer :: ndim
      real(8),intent(inout) :: cmoa(nao*nao), overinv(nao*nao)
      real(8),allocatable :: overlap(:), work1(:), work2(:), eigen(:)
!
      ndim= min(nmo,neleca+5)
!
! Set arrays
!
      call memset(nao*nao_g+ndim*ndim+nao*ndim+ndim)
      allocate(overlap(nao*nao_g),work1(ndim*ndim),work2(nao*ndim),eigen(ndim))
!
! Calculate overlap integrals between previous and present bases
!
      call calcover2(overlap)
!
! Project orbitals from Huckel to SCF
!
      call projectmo2(cmoa,overinv,overlap,work1,work2,eigen,ndim)
!
! Unset arrays
!
      deallocate(overlap,work1,work2,eigen)
      call memunset(nao*nao_g+ndim*ndim+nao*ndim+ndim)
!
      return
end


!----------------------------------------------------------------------
  subroutine projectmo2(cmoa,overinv,overlap,work1,work2,eigen,ndim)
!----------------------------------------------------------------------
!
! Project orbitals from Huckel to SCF
!    C1= S11^-1 * S12 * C2 [C2t * S12t * S11^-1 * S12 * C2]^-1/2
!
! In this routine, nao >= ndim is assumed.
!
! Inout : cmoa (previous and updated orbitals)
!         overinv (overlap integral inverse matrix of SCF basis set)
!         overlap (overlap integral of guess and SCF basis sets)
!
      use iofile, only : iout
      use guess, only : nao_g
      use basis, only : nao
      implicit none
      integer,intent(in) :: ndim
      integer :: i, j
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(inout) :: cmoa(nao,nao), overinv(nao,nao), overlap(nao,nao_g)
      real(8),intent(out) :: work1(ndim,ndim), work2(ndim,nao), eigen(ndim)
      real(8) :: eigeninv
!
      if(nao < ndim) then
        write(iout,'(" Error! Nao is less than ndim in projectmo2.")')
        call iabort
      endif
!
! Calculate S12*C2
!
      call dgemm('N','N',nao,ndim,nao_g,one,overlap,nao,cmoa,nao_g,zero,work2,nao)
!
! Calculate S11-1*S12*C2
!
      call dsymm('L','U',nao,ndim,one,overinv,nao,work2,nao,zero,overlap,nao)
!
! Calculate C2t*S12t*S11-1*S12*C2
!
      call dgemm('T','N',ndim,ndim,nao,one,work2,nao,overlap,nao,zero,work1,ndim)
!
! Calculate (C2t*S12t*S11-1*S12*C2)-1/2
!
      call diag('V','U',ndim,work1,ndim,eigen)
!$OMP parallel do private(eigeninv)
      do i= 1,ndim
        eigeninv= one/sqrt(eigen(i))
        do j= 1,ndim
          work2(j,i)= work1(j,i)*eigeninv
        enddo
      enddo
!$OMP end parallel do
      call dgemm('N','T',ndim,ndim,ndim,one,work1,ndim,work2,ndim,zero,overinv,ndim)
!
! Calculate C1
!
      call dgemm('N','N',nao,ndim,ndim,one,overlap,nao,overinv,ndim,zero,cmoa,nao)
      return
end


