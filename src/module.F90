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
!------------------
  module modparam
!------------------
!
! maximum sizes for AOs (basis functions), MOs, atoms, shells
!
      implicit none
      integer,parameter :: mxao=10000
      integer,parameter :: mxmo=10000
      integer,parameter :: mxatom=1000
      integer,parameter :: mxshell=5000
      integer,parameter :: mxprim=20000
      integer,parameter :: mxprsh=30
      integer,parameter :: mxang=7
end


!----------------
  module modjob
!----------------
      implicit none
      character(len=16) :: method, runtype, scftype
end


!-----------------
  module modunit
!-----------------
      implicit none
      real(8),parameter :: toang= 0.5291772108D+00, tobohr= 1.889726125D+00
      real(8),parameter :: todebye= 2.541746D+00
      logical :: bohr
end


!---------------------
  module modparallel
!---------------------
      implicit none
      integer :: nproc1, myrank1, nproc2, myrank2, mpi_comm1, mpi_comm2
      logical :: master
!
      interface checkintsize
        subroutine checkintsize4(isize)
          integer(selected_int_kind(9)),intent(out) :: isize
        end subroutine checkintsize4
        subroutine checkintsize8(isize)
          integer(selected_int_kind(18)),intent(out) :: isize
        end subroutine checkintsize8
      end interface checkintsize
!
end


!-------------------
  module modiofile
!-------------------
      implicit none
      integer,parameter :: input=10, icheck=20, maxline=100000
      character(len=16) :: version
      character(len=64) :: check
end


!-------------------
  module modtclock
!-------------------
      implicit none
      integer ::  iwall0, iwall1
      real(8) :: cpu0, cpu1
end


!---------------------
  module modmolecule
!---------------------
!
! natom     : number of atoms
! numatomic : atomic number
! coord     : Cartesian coordinate 
! znuc      : atomic charge
! neleca    : number of alpha electrons
! nelecb    : number of beta electrons
! nmo       : number of molecular orbitals
! ncore     : number of core orbitals
! multi     : spin multiplicity
! charge    : molecular charge
!
      use modparam, only : mxatom
      implicit none
      integer :: numatomic(mxatom), natom, neleca, nelecb, nmo, ncore, multi
      real(8) :: coord(3,mxatom), znuc(mxatom), coordold(3,mxatom), charge
end


!-------------------
  module modthresh
!-------------------
!
! threshex    : threshold for overlap of two basis functions
! threshover  : threshold for linear depencency of basis functions
! threshatm   : threshold for distance of atoms
! cutint2     : threshold for two-electron integrals
! threshsoscf : threshold for second-order SCF
! threshqc    : threshold for quadratically convergent SCF
! threshweight: threshold for weight at a grid point
! threshrho   : threshold for density at a grid point
! threshdfock : threshold for functional at a grid point
! threshmp2cphf : threshold for MP2-CPHF
!
      implicit none
      real(8),parameter :: threshex=30.0D+00
      real(8) :: threshover, threshatom, threshdiis
      real(8) :: cutint2, threshsoscf, threshqc
      real(8) :: threshweight, threshrho, threshdfock, threshdftao 
      real(8) :: threshmp2cphf
      character(len=16) :: precision
end


!-------------------
  module modmemory
!-------------------
      integer :: memmax, memused, memusedmax
      character(len=16) :: memory, mem
end


!------------------
  module modbasis
!------------------
!
! nshell : Number of shells
! nao    : Number of AOs (=contracted basis functions)
! nprim  : Number of primitive basis functions
! ex    : Exponents of basis functions
! coeff : Normalized coefficients of basis functions
! coeffinp : Input coefficients of basis functions
!
! locprim : Starting address of primitive basis functions for a shell
! locbf   : Starting address of basis functions for a shell
! locatom : Atom center for a shell
! mprim   : Number of basis primitive functions for a shell
! mbf     : Number of basis functions for a shell
! mtype   : Type of a basis shell (s=0, p=1, d=2, f=3,...)
! basis   : Name of basis functions
! atombasis : Type of basis set for each atom
! exgen : Exponents of basis functions from input file
! coeffgen: Coefficients of basis functions from input file
! locgenprim: Starting address of primitive basis functions for a shell from input file
! mgenprim : Number of basis primitive functions for a shell from input file
! mgentype : Type of a basis shell (s=0, p=1, d=2, f=3,...) from input file
! locgenshell :  Starting address of basis shells
! ngenshell : Number of basis shells for an atom
!
      use modparam, only : mxprim, mxshell
      implicit none
      integer :: nshell, nao, nprim
      integer :: locprim(mxshell+1), locbf(mxshell+1), locatom(mxshell)
      integer :: mprim(mxshell),  mbf(mxshell), mtype(mxshell)
      integer :: locgenprim(mxshell+1), mgenprim(mxshell), mgentype(mxshell)
      integer :: locgenshell(-9:112), ngenshell(-9:112)
      real(8) :: ex(mxprim), coeff(mxprim), coeffinp(mxprim)
      real(8) :: exgen(mxprim), coeffgen(mxprim)
      character(len=16) :: basis, atombasis(-9:112)
      logical :: spher
end


!-----------------
  module modwarn
!-----------------
!
! nwarn : number of warnings
!
      implicit none
      integer :: nwarn
end


!------------------
  module modguess
!------------------
!
! These variables are for guess calculations.
! (Extended Huckel)
!
! nshell_v : the number of valence shells
! nao_v    : the number of valence AOs (=contracted basis functions)
! nprim_v  : the number of valence primitive basis functions
! nmo_v    : the number of valence MOs
! nshell_g : the number of shells
! nao_g    : the number of AOs (=contracted basis functions)
! nprim_g  : the number of primitive basis functions
! nmo_g    : the number of MOs
! nao_c    : the number of core MOs
!
! ex_g    : exponents
! coeff_g : coefficients 
! coord_g : coordinate
!
! locprim_g :(starting address of primitive basis functions for a shell) - 1
! locbf_g   :(starting address of basis functions for a shell) - 1
! locatom_g : atom center for a shell
! mprim_g   : the number of basis primitive functions for a shell
! mbf_g     : the number of basis funcions for a shell
! mtype_g   : the type of a basis shell (s=0, p=1, d=2, f=3,...)
! func_g    : the name of basis functions
! guess     : type of initail guess
!
      use modparam, only : mxprim, mxshell, mxatom
      implicit none
      integer :: nshell_v, nao_v, nprim_v, nmo_v, nshell_g, nao_g, nprim_g, nmo_g, nao_c
      integer :: locprim_g(mxshell+1), locbf_g(mxshell+1), locatom_g(mxshell)
      integer :: mprim_g(mxshell), mbf_g(mxshell), mtype_g(mxshell)
      integer :: nshell_gcore, nao_gcore, nprim_gcore
      integer :: locprim_gcore(mxshell+1), locbf_gcore(mxshell+1), locatom_gcore(mxshell)
      integer :: mprim_gcore(mxshell), mbf_gcore(mxshell), mtype_gcore(mxshell)
      real(8) :: ex_g(mxprim), coeff_g(mxprim), coord_g(3,mxatom)
      real(8) :: ex_gcore(mxprim), coeff_gcore(mxprim), coord_gcore(3,mxatom)
      logical :: spher_g
      character(len=16) :: guess
end


!------------------
  module modprint
!------------------
!
! iprint : print option
!    = 0 : minimal output
!    = 1 : standard output
!    = 2 : verbose output
!
      integer :: iprint
end


!-------------------
  module modenergy
!-------------------
!  energies 
!
      real(8) :: enuc, eelec, escf, escfe, emp2, escsmp2
end


!----------------
  module modscf
!----------------
!
! These variables are for SCF calculations.
!
      implicit none
      integer :: maxiter, maxdiis, maxsoscf, maxqc, maxqcdiag, maxqcdiagsub
      real(8) :: dconv
      logical :: fdiff, extrap
      character(len=16) :: scfconv
end


!----------------
  module moddft
!----------------
      implicit none
      integer :: idftex, idftcor, nrad, nleb
      real(8) :: hfexchange, bqrad(9)
end


!----------------
  module modmp2
!----------------
      implicit none
      integer :: ncore, nvfz, maxmp2diis, maxmp2iter
end


!-----------------
  module modatom
!-----------------
!
! Covalent radii  H       : Bohr radius
!                 He - Cn : P. Pyykko, M. Atsumi, Chem. Eur. J., 186 (2009) 15.
!                 
!
      implicit none
      real(8) :: atomrad(-9:112)
      data atomrad/ &
&     1.06D+00, 1.06D+00, 1.06D+00, 1.06D+00, 1.06D+00, 1.06D+00, 1.06D+00, 1.06D+00, 1.06D+00, &
&     0.00D+00, &
&     0.53D+00, 0.46D+00, 1.33D+00, 1.02D+00, 0.85D+00, 0.75D+00, 0.71D+00, 0.63D+00, 0.64D+00, &
&     0.67D+00, 1.55D+00, 1.39D+00, 1.26D+00, 1.16D+00, 1.11D+00, 1.03D+00, 0.99D+00, 0.96D+00, &
&     1.96D+00, 1.71D+00, 1.48D+00, 1.36D+00, 1.34D+00, 1.22D+00, 1.19D+00, 1.16D+00, 1.11D+00, &
&     1.10D+00, 1.12D+00, 1.18D+00, 1.24D+00, 1.21D+00, 1.21D+00, 1.16D+00, 1.14D+00, 1.17D+00, &
&     2.10D+00, 1.85D+00, 1.63D+00, 1.54D+00, 1.47D+00, 1.38D+00, 1.28D+00, 1.25D+00, 1.25D+00, &
&     1.20D+00, 1.28D+00, 1.36D+00, 1.42D+00, 1.40D+00, 1.40D+00, 1.36D+00, 1.33D+00, 1.31D+00, &
&     2.32D+00, 1.96D+00, 1.80D+00, 1.63D+00, 1.76D+00, 1.74D+00, 1.73D+00, 1.72D+00, 1.68D+00, &
&     1.69D+00, 1.68D+00, 1.67D+00, 1.66D+00, 1.65D+00, 1.64D+00, 1.70D+00, 1.62D+00, 1.52D+00, &
&     1.46D+00, 1.37D+00, 1.31D+00, 1.29D+00, 1.22D+00, 1.23D+00, 1.24D+00, 1.33D+00, 1.44D+00, &
&     1.44D+00, 1.51D+00, 1.45D+00, 1.47D+00, 1.42D+00, 2.23D+00, 2.01D+00, 1.86D+00, 1.75D+00, &
&     1.69D+00, 1.70D+00, 1.71D+00, 1.72D+00, 1.66D+00, 1.66D+00, 1.68D+00, 1.68D+00, 1.65D+00, &
&     1.67D+00, 1.73D+00, 1.76D+00, 1.61D+00, 1.57D+00, 1.49D+00, 1.43D+00, 1.41D+00, 1.34D+00, &
&     1.29D+00, 1.28D+00, 1.21D+00, 1.22D+00/
end


!----------------
  module modopt
!----------------
      implicit none
      integer :: nopt
      real(8) :: optconv
      logical :: cartesian
end


!-----------------
  module modprop
!-----------------
      implicit none
      logical :: octupole
end


!----------------
  module modecp
!----------------
!
! izcore    : number of core electrons per atom for ECP calculation
!
      use modparam, only : mxprim, mxshell, mxatom
      implicit none
      integer,parameter :: lfunc=16, nterm1=625291, nterm2=26841
      integer :: maxangecp(mxatom), izcore(mxatom)=0, mtypeecp(mxprim)
      integer :: locecp(0:5,mxatom), mprimecp(0:5,mxatom)
      integer :: maxgenangecp(-9:112), izgencore(-9:112), mgentypeecp(mxprim)
      integer :: locgenecp(0:5,-9:112), mgenprimecp(0:5,-9:112)
      integer :: lmf(122), lmx(581), lmy(581), lmz(581), nbfirst(0:7), nx(84), ny(84), nz(84)
      real(8) :: execp(mxprim), coeffecp(mxprim), zlm(581)
      real(8) :: exgenecp(mxprim), coeffgenecp(mxprim)
      character(lfunc) :: ecp, atomecp(-9:112)
      logical :: flagecp
      data lmf/1, 2,3,4, 5,7,8,10,11, 12,14,16,18,20,22,23, 25,28,30,34,36,39,41,43,45,&
&              47,50,53,57,61,64,67,70,72,76,78, 81,85,88,94,98,104,107,111,114,117,121,125,128,&
&              131,135,139,145,151,157,163,167,171,175,178,184,188,194,197,&
&              201,206,210,218,224,233,239,247,251,256,260,264,270,276,282,288,292,&
&              296,301,306,314,322,331,340,348,356,361,366,371,375,383,389,398,404,412,416,&
&              421,427,432,442,450,462,471,483,491,501,506,512,517,522,530,538,547,556,564,&
&              572,577,582/
      data lmx/0, 1,0,0, 2,0,1,0,0,0,1, 3,1,2,0,1,1,0,0,0,0,1,2,0,&
&              4,2,0,3,1,2,0,0,2,1,1,0,0,0,0,0,1,1,2,0,3,1, 5,3,1,0,2,4,3,1,1,3,&
&              2,0,2,0,1,1,1,0,0,0,0,0,0,1,1,2,2,0,0,3,1,4,2,0, 6,4,2,0,5,3,1,4,&
&              2,0,4,2,0,3,1,1,3,2,0,2,0,2,0,1,1,1,0,0,0,0,0,0,0,1,1,1,2,0,2,0,3,&
&              1,3,1,4,2,0,5,3,1, 7,3,5,1,6,4,0,2,5,3,1,5,3,1,4,2,0,4,2,0,3,1,3,&
&              1,3,1,2,0,2,0,2,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,2,0,2,0,2,0,3,1,3,&
&              1,4,2,0,4,2,0,5,3,1,0,4,2,6, 8,6,4,2,0,7,5,3,1,6,4,2,0,6,4,2,0,5,&
&              3,1,5,3,1,4,0,2,4,0,2,4,0,2,3,1,3,1,3,1,2,0,2,0,2,0,2,0,1,1,1,1,0,&
&              0,0,0,0,0,0,0,0,1,1,1,1,2,0,2,0,2,0,3,1,3,1,3,1,4,2,0,4,2,0,5,3,1,&
&              5,3,1,6,4,2,0,7,5,3,1, 9,7,5,3,1,8,6,4,2,0,7,5,3,1,7,5,3,1,6,4,2,&
&              0,6,4,2,0,5,3,1,5,3,1,5,3,1,4,0,2,4,0,2,4,0,2,3,1,3,1,3,1,3,1,2,0,&
&              2,0,2,0,2,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,2,0,2,0,2,0,2,0,&
&              3,1,3,1,3,1,4,2,0,4,2,0,4,2,0,5,3,1,5,3,1,6,4,2,0,6,4,2,0,7,5,3,1,&
&              8,6,4,2,0, 10,8,6,4,2,0,9,7,5,3,1,8,6,4,2,0,8,6,4,2,0,7,5,3,1,7,5,&
&              3,1,6,4,2,0,6,4,2,0,6,4,2,0,5,3,1,5,3,1,5,3,1,4,0,2,4,0,2,4,0,2,4,&
&              0,2,3,1,3,1,3,1,3,1,2,0,2,0,2,0,2,0,2,0,1,1,1,1,1,0,0,0,0,0,0,0,&
&              0,0,0,0,1,1,1,1,1,2,0,2,0,2,0,2,0,3,1,3,1,3,1,3,1,0,2,4,0,2,4,0,2,&
&              4,5,3,1,5,3,1,5,3,1,6,4,2,0,6,4,2,0,7,5,3,1,7,5,3,1,8,6,4,2,0,9,7,5,3,1/
      data lmy/0, 0,0,1, 0,2,0,0,0,1,1, 0,2,0,2,0,0,0,0,1,1,1,1,3,&
&              0,2,4,0,2,0,2,2,0,0,0,0,0,0,1,1,1,1,1,3,1,3, 0,2,4,4,2,0,0,2,2,0,&
&              0,2,0,2,0,0,0,0,0,0,1,1,1,1,1,1,1,3,3,1,3,1,3,5, 0,2,4,6,0,2,4,0,&
&              2,4,0,2,4,0,2,2,0,0,2,0,2,0,2,0,0,0,0,0,0,0,1,1,1,1,1,1,1,3,1,3,1,&
&              3,1,3,1,3,5,1,3,5, 0,4,2,6,0,2,6,4,0,2,4,0,2,4,0,2,4,0,2,4,0,2,0,&
&              2,0,2,0,2,0,2,0,2,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,1,&
&              3,1,3,5,1,3,5,1,3,5,7,3,5,1, 0,2,4,6,8,0,2,4,6,0,2,4,6,0,2,4,6,0,&
&              2,4,0,2,4,0,4,2,0,4,2,0,4,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,0,0,0,0,&
&              0,0,0,0,1,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,1,3,1,3,1,3,5,1,3,5,1,3,5,&
&              1,3,5,1,3,5,7,1,3,5,7, 0,2,4,6,8,0,2,4,6,8,0,2,4,6,0,2,4,6,0,2,4,&
&              6,0,2,4,6,0,2,4,0,2,4,0,2,4,0,4,2,0,4,2,0,4,2,0,2,0,2,0,2,0,2,0,2,&
&              0,2,0,2,0,2,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,&
&              1,3,1,3,1,3,1,3,5,1,3,5,1,3,5,1,3,5,1,3,5,1,3,5,7,1,3,5,7,1,3,5,7,&
&              1,3,5,7,9, 0,2,4,6,8,10,0,2,4,6,8,0,2,4,6,8,0,2,4,6,8,0,2,4,6,0,2,&
&              4,6,0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,0,2,4,0,2,4,0,4,2,0,4,2,0,4,2,0,&
&              4,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,0,0,0,0,0,0,0,0,0,0,1,&
&              1,1,1,1,1,1,1,1,1,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,5,3,1,5,3,1,5,3,&
&              1,1,3,5,1,3,5,1,3,5,1,3,5,7,1,3,5,7,1,3,5,7,1,3,5,7,1,3,5,7,9,1,3,5,7,9/
      data lmz/0, 0,1,0, 0,0,1,2,0,1,0, 0,0,1,1,2,0,3,1,2,0,1,0,0,&
&              0,0,0,1,1,2,2,0,0,3,1,4,2,0,3,1,2,0,1,1,0,0, 0,0,0,1,1,1,2,2,0,0,&
&              3,3,1,1,4,2,0,5,3,1,4,2,0,3,1,2,0,2,0,1,1,0,0,0, 0,0,0,0,1,1,1,2,&
&              2,2,0,0,0,3,3,1,1,4,4,2,2,0,0,5,3,1,6,4,2,0,5,3,1,4,2,0,3,3,1,1,2,&
&              2,0,0,1,1,1,0,0,0, 0,0,0,0,1,1,1,1,2,2,2,0,0,0,3,3,3,1,1,1,4,4,2,&
&              2,0,0,5,5,3,3,1,1,6,4,2,0,7,5,3,1,6,4,2,0,5,3,1,4,4,2,2,0,0,3,3,1,&
&              1,2,2,2,0,0,0,1,1,1,0,0,0,0, 0,0,0,0,0,1,1,1,1,2,2,2,2,0,0,0,0,3,&
&              3,3,1,1,1,4,4,4,2,2,2,0,0,0,5,5,3,3,1,1,6,6,4,4,2,2,0,0,7,5,3,1,8,&
&              6,4,2,0,7,5,3,1,6,4,2,0,5,5,3,3,1,1,4,4,2,2,0,0,3,3,3,1,1,1,2,2,2,&
&              0,0,0,1,1,1,1,0,0,0,0, 0,0,0,0,0,1,1,1,1,1,2,2,2,2,0,0,0,0,3,3,3,&
&              3,1,1,1,1,4,4,4,2,2,2,0,0,0,5,5,5,3,3,3,1,1,1,6,6,4,4,2,2,0,0,7,7,&
&              5,5,3,3,1,1,8,6,4,2,0,9,7,5,3,1,8,6,4,2,0,7,5,3,1,6,6,4,4,2,2,0,0,&
&              5,5,3,3,1,1,4,4,4,2,2,2,0,0,0,3,3,3,1,1,1,2,2,2,2,0,0,0,0,1,1,1,1,&
&              0,0,0,0,0, 0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,0,0,0,0,0,3,3,3,3,1,1,&
&              1,1,4,4,4,4,2,2,2,2,0,0,0,0,5,5,5,3,3,3,1,1,1,6,6,6,4,4,4,2,2,2,0,&
&              0,0,7,7,5,5,3,3,1,1,8,8,6,6,4,4,2,2,0,0,9,7,5,3,1,10,8,6,4,2,0,9,&
&              7,5,3,1,8,6,4,2,0,7,7,5,5,3,3,1,1,6,6,4,4,2,2,0,0,5,5,5,3,3,3,1,1,&
&              1,4,4,4,2,2,2,0,0,0,3,3,3,3,1,1,1,1,2,2,2,2,0,0,0,0,1,1,1,1,1,0,0,0,0,0/
      data nbfirst/1,2,5,11,21,36,57,85/
      data nx/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,&
&             1,0,2,2,0,2,1,1, 5,0,0,4,4,1,0,1,0,3,3,2,0,2,0,3,1,1,2,2,1,&
&             6,0,0,5,5,1,0,1,0,4,4,2,0,2,0,4,1,1,3,3,0,3,3,2,1,2,1,2/
      data ny/0,0,1,0,0,2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,&
&             0,1,2,0,2,1,2,1, 0,5,0,1,0,4,4,0,1,2,0,3,3,0,2,1,3,1,2,1,2,&
&             0,6,0,1,0,5,5,0,1,2,0,4,4,0,2,1,4,1,3,0,3,2,1,3,3,1,2,2/
      data nz/0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,1,2,2,1,0,0,4,0,1,0,1,&
&             3,3,0,2,2,1,1,2, 0,0,5,0,1,0,1,4,4,0,2,0,2,3,3,1,1,3,1,2,2,&
&             0,0,6,0,1,0,1,5,5,0,2,0,2,4,4,1,1,4,0,3,3,1,2,1,2,3,3,2/
end

