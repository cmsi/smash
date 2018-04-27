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
  subroutine calc2eri(twoeri,ish,jsh,ksh,lsh,maxdim,databasis)
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
      use modjob, only : threshex
      use modtype, only : typebasis
      implicit none
      type(typebasis),intent(inout) :: databasis
      integer,intent(in) :: ish, jsh, ksh, lsh, maxdim
      integer :: nangijkl(4), nprimijkl(4), nbfijkl(4)
      integer :: iloc, jloc, kloc, lloc, i, j, k, l
      real(8),intent(out) :: twoeri(maxdim,maxdim,maxdim,maxdim)
      real(8) :: xyzijkl(3,4), exijkl(mxprsh,4), coijkl(mxprsh,4)
!
      nangijkl(1)= databasis%mtype(ish)
      nangijkl(2)= databasis%mtype(jsh)
      nangijkl(3)= databasis%mtype(ksh)
      nangijkl(4)= databasis%mtype(lsh)
!
      iloc  = databasis%locprim(ish)
      jloc  = databasis%locprim(jsh)
      kloc  = databasis%locprim(ksh)
      lloc  = databasis%locprim(lsh)
      nbfijkl(1)= databasis%mbf(ish)
      nbfijkl(2)= databasis%mbf(jsh)
      nbfijkl(3)= databasis%mbf(ksh)
      nbfijkl(4)= databasis%mbf(lsh)
      nprimijkl(1)= databasis%mprim(ish)
      nprimijkl(2)= databasis%mprim(jsh)
      nprimijkl(3)= databasis%mprim(ksh)
      nprimijkl(4)= databasis%mprim(lsh)
!
      do i= 1,3
        xyzijkl(i,1)= coord(i,databasis%locatom(ish))
        xyzijkl(i,2)= coord(i,databasis%locatom(jsh))
        xyzijkl(i,3)= coord(i,databasis%locatom(ksh))
        xyzijkl(i,4)= coord(i,databasis%locatom(lsh))
      enddo
!
      do i= 1,nprimijkl(1)
        exijkl(i,1)=databasis%ex(iloc+i)
        coijkl(i,1)=databasis%coeff(iloc+i)
      enddo
      do j= 1,nprimijkl(2)
        exijkl(j,2)=databasis%ex(jloc+j)
        coijkl(j,2)=databasis%coeff(jloc+j)
      enddo
      do k= 1,nprimijkl(3)
        exijkl(k,3)=databasis%ex(kloc+k)
        coijkl(k,3)=databasis%coeff(kloc+k)
      enddo
      do l= 1,nprimijkl(4)
        exijkl(l,4)=databasis%ex(lloc+l)
        coijkl(l,4)=databasis%coeff(lloc+l)
      enddo
!
      call int2elec(twoeri,exijkl,coijkl,xyzijkl,nprimijkl,nangijkl,nbfijkl,maxdim, &
&                   mxprsh,threshex)
      return
end


!-----------------------------------------------------------------------
  subroutine calcschwarzeri(xint,xinttmp,maxdim,nproc,myrank,mpi_comm,databasis)
!-----------------------------------------------------------------------
!
! Driver of (ij|ij) integral calculation
!
! Out : xint    ((ij|ij) integrals for Schwarz screening)
!       xinttmp (Work array)
!
      use modtype, only : typebasis
      implicit none
      type(typebasis),intent(inout) :: databasis
      integer,intent(in) :: maxdim, nproc, myrank, mpi_comm
      integer :: ish, jsh, nbfi, nbfj, i, j, ii, ij
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, two=2.0D+00
      real(8),intent(out) :: xint(databasis%nshell*(databasis%nshell+1)/2)
      real(8),intent(out) :: xinttmp(databasis%nshell*(databasis%nshell+1)/2)
      real(8) :: twoeri(maxdim,maxdim,maxdim,maxdim), xintmax, val
!
      xinttmp(:)= zero
!
!$OMP parallel private(jsh,xintmax,twoeri,nbfi,nbfj,i,j,val,ij,ii)
      do ish= databasis%nshell-myrank,1,-nproc
        nbfi= databasis%mbf(ish)
        ii  = ish*(ish-1)/2
!$OMP do
        do jsh = 1,ish
          nbfj = databasis%mbf(jsh)
          call calc2eri(twoeri,ish,jsh,ish,jsh,maxdim,databasis)
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
      call para_allreducer(xinttmp,xint,databasis%nshell*(databasis%nshell+1)/2,mpi_comm)
!
      return
end
