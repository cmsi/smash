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
!-----------------------------------------------------
  subroutine calc2eri(twoeri,ish,jsh,ksh,lsh,maxdim)
!-----------------------------------------------------
!
! Driver of two-electron repulsion integrals 
! 
! In  : ish, jsh, ksh, lsh (Shell indices)
! Out : twoeri (Two-electron repulsion integrals)
!       idxeri (Indices of two-electron repulsion integrals)
!       numeri (Number of two-electron repulsion integrals)
!
      use modparam, only : mxprsh
      use modmolecule, only : coord
      use modbasis, only : locatom, locprim, mprim, mbf, mtype, ex, coeff
      use modthresh, only : threshex
      implicit none
      integer,intent(in) :: ish, jsh, ksh, lsh, maxdim
      integer :: nangijkl(4), nprimijkl(4), nbfijkl(4)
      integer :: iloc, jloc, kloc, lloc, i, j, k, l
      real(8),intent(out) :: twoeri(maxdim,maxdim,maxdim,maxdim)
      real(8) :: xyzijkl(3,4), exijkl(mxprsh,4), coijkl(mxprsh,4)
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
        xyzijkl(i,1)= coord(i,locatom(ish))
        xyzijkl(i,2)= coord(i,locatom(jsh))
        xyzijkl(i,3)= coord(i,locatom(ksh))
        xyzijkl(i,4)= coord(i,locatom(lsh))
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
      call int2elec(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                   mxprsh,threshex)
      return
end


!-----------------------------------------------------------------------
  subroutine calcschwarzeri(xint,xinttmp,maxdim,nproc,myrank,mpi_comm)
!-----------------------------------------------------------------------
!
! Driver of (ij|ij) integral calculation
!
! Out : xint    ((ij|ij) integrals for Schwarz screening)
!       xinttmp (Work array)
!
      use modbasis, only : nshell, mbf
      implicit none
      integer,intent(in) :: maxdim, nproc, myrank, mpi_comm
      integer :: ish, jsh, nbfi, nbfj, i, j, ii, ij
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, two=2.0D+00
      real(8),intent(out) :: xint(nshell*(nshell+1)/2)
      real(8),intent(out) :: xinttmp(nshell*(nshell+1)/2)
      real(8) :: twoeri(maxdim,maxdim,maxdim,maxdim), xintmax, val
!
      xinttmp(:)= zero
!
!$OMP parallel private(jsh,xintmax,twoeri,nbfi,nbfj,i,j,val,ij,ii)
      do ish= nshell-myrank,1,-nproc
        nbfi= mbf(ish)
        ii  = ish*(ish-1)/2
!$OMP do
        do jsh = 1,ish
          nbfj = mbf(jsh)
          call calc2eri(twoeri,ish,jsh,ish,jsh,maxdim)
          xintmax= zero
          do i= 1,nbfi
            do j= 1,nbfj
              val= twoeri(j,i,j,i)
              if(val.gt.xintmax) xintmax= val
            enddo
          enddo
          ij= ii+jsh
          xinttmp(ij)=sqrt(xintmax)
        enddo
!$OMP end do
      enddo
!$OMP end parallel
!
      call para_allreducer(xinttmp,xint,nshell*(nshell+1)/2,mpi_comm)
!
      return
end
