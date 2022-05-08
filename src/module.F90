! Copyright 2014-2021  Kazuya Ishimura
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
! Maximum sizes for AOs (basis functions), MOs, atoms, shells
! Factors of unit conversion
!
! Fundamental Physical Constants: 2018 CODATA
!   Eite Tiesinga, Peter J. Mohr, David B. Newell, and Barry N. Taylor (2020),
!   "The 2018 CODATA Recommended Values of the Fundamental Physical Constants"
!   (Web Version 8.1). Database developed by J. Baker, M. Douma, and S.
!   Kotochigova. Available at http://physics.nist.gov/constants, National
!   Institute of Standards and Technology, Gaithersburg, MD 20899.
!
      implicit none
      integer,parameter :: mxao=10000
      integer,parameter :: mxmo=10000
      integer,parameter :: mxatom=1000
      integer,parameter :: mxshell=5000
      integer,parameter :: mxprim=20000
      integer,parameter :: mxprsh=30
      integer,parameter :: mxang=7
      integer,parameter :: maxline=100000
      real(8),parameter :: toang= 0.5291772109D+00, tobohr= 1.889726125D+00
      real(8),parameter :: todebye= 2.541746473D+00
end


!----------------
  module modecp
!----------------
!
! ECP data
!
      implicit none
      integer,parameter :: nterm1=625291, nterm2=26841
      integer :: lmf(122), lmx(581), lmy(581), lmz(581), nbfirst(0:7), nx(84), ny(84), nz(84)
      real(8) :: zlm(581)
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


!-----------------
  module modtype
!-----------------
!
      use modparam, only : mxprim, mxshell, mxatom
      implicit none
!
! Basis set and ECP
!
! nshell   : Number of shells
! nao      : Number of AOs (=contracted basis functions)
! nprim    : Number of primitive basis functions
! locprim  : Starting address of primitive basis function for each shell
! locbf    : Starting address of basis function for each shell
! locatom  : Atom center for each shell
! mprim    : Number of basis primitive function for each shell
! mbf      : Number of basis function for each shell
! mtype    : Type of each basis shell (s=0, p=1, d=2, f=3,...)
! locecp   : Starting address of ECP function for each shell
! mprimecp : Number of ECP functions for each shell
! maxangecp: Maximum ECP angular momentum for each atom
! mtype    : Type of each ECP shell (s=0, p=1, d=2, f=3,...)
! izcore   : number of core electrons per atom for ECP calculation
! ex       : Exponents of basis functions
! coeff    : Normalized coefficients of basis functions
! coeffinp : Input coefficients of basis functions
! basis    : Name of basis functions
! spher    : Spherical Harmonics (True) or Cartesian (False)
      type typebasis
        integer :: nshell, nao, nprim
        integer :: locprim(mxshell+1), locbf(mxshell+1), locatom(mxshell)
        integer :: mprim(mxshell),  mbf(mxshell), mtype(mxshell)
        integer :: locecp(0:5,mxatom), mprimecp(0:5,mxatom)
        integer :: maxangecp(mxatom), izcore(mxatom)=0, mtypeecp(mxprim)
        real(8) :: ex(mxprim), coeff(mxprim), coeffinp(mxprim)
        real(8) :: execp(mxprim), coeffecp(mxprim)
        character(len=32) :: basis='STO-3G', ecp=''
        logical :: spher=.true.
      end type
!
! Job control
!
! iprint       : Print control
! maxiter      : Maximum number of SCF iterations
! maxdiis      : Maximum number of DIIS cycles
! maxsoscf     : Maximum number of Second-order SCF cycles
! maxqc        : Maximum number of Quadratic convergent (QC) cycles
! maxqcdiag    : Maximum number of Davidson diagonalization cycles in QC
! maxqcdiagsub : Maximum number of trial vectors of Davidson diagonalization in QC
! idftex       : Exchange functional
! idftco       : Correlation functional
! nrad         : Number of radial points in the Euler-MacLaurin quadrature
! nleb         : Number of anglualr points in the Levedev grids
! ncore        : Number of core orbitals
! nvfz         : Number of frozen virtual orbitals
! maxmp2diis   : Maximum number of DIIS iterations for MP2 CPHF calculation
! maxmp2iter   : Maximum number of iterations for MP2 CPHF calculation
! nopt         : Maximum number of geometry optimization cycles
! dconv        : Threshold for SCF density convergence
! optconv      : Threshold for geometry optimization convergence
! fbond        : Scaling factor of atom radii for bond detection
! hfexchange   : Factor of Hartree-Fock Exchange in DFT
! bqrad        : Radii of ghost (Bq) atoms
! threshex     : Threshold for overlap of two basis functions
! threshover   : Threshold for linear depencency of basis functions
! threshatm    : Threshold for distance of atoms
! cutint2      : Threshold for two-electron integrals
! threshdiis   : Threshold for DIIS start
! threshsoscf  : Threshold for second-order SCF
! threshqc     : Threshold for quadratically convergent SCF
! threshweight : Threshold for weight at a grid point
! threshrho    : Threshold for density at a grid point
! threshdfock  : Threshold for functional at a grid point
! threshmp2cphf: Threshold for MP2-CPHF
! method       : Hamiltonian
! runtype      : Calculation type
! scftype      : RHF or UHF
! memory       : Memory size
! version      : SMASH version
! guess        : Guess type
! precision    : Computational precision (high, medium, low)
! scfconv      : SCF convergence method
! print        : output control
! pop          : population control
! multipole    : Multipole memoent control
! check        : Checkpoint file
! xyz          : Xyz file
! bohr         : Length unit (True: atomic unit, False: Angstrom)
! flagecp      : Flag for ECP calculation
! extrap       : Pople extrapolation of Fock matrix
! cartesian    : Type of coordinate system (True: Cartesian, False: redundant coordinate)
      type typejob
        integer :: iprint=0
        integer :: maxiter=150, maxdiis=20, maxsoscf=20, maxqc=15, maxqcdiag=100, maxqcdiagsub=10
        integer :: idftex=0, idftcor=0, nrad=0, nleb=0
        integer :: ncore=-1, nvfz=0, maxmp2diis=20, maxmp2iter=100
        integer :: nopt=100
        real(8) :: dconv=-1.0D+00, optconv=-1.0D+00, fbond=1.20D+00 
        real(8) :: hfexchange=1.0D+00, bqrad(9)=1.0D+00
        real(8) :: threshex=30.0D+00, threshover=1.0D-06, threshatom=2.0D-01
        real(8) :: cutint2=-1.0D+00, threshdiis=6.0D-01, threshsoscf=0.25D+00, threshqc=1.0D-05
        real(8) :: threshweight=-1.0D+00, threshrho=-1.0D+00, threshdfock=-1.0D+00
        real(8) :: threshdftao=-1.0D+00, threshmp2cphf=1.0D-10
        character(len=32) :: method='HF', runtype='ENERGY', scftype='RHF', memory=''
        character(len=32) :: version='3.0.1', guess='HUCKEL', precision='MEDIUM'
        character(len=32) :: scfconv='DIIS', print='', pop='MULLIKEN', multipole='DIPOLE'
        character(len=256) :: check='', xyz=''
        logical :: bohr=.false., flagecp=.false., extrap=.false.
        logical :: cartesian=.false.
!       logical :: fdiff=.true.
      end type
!
! Molecule Information
!
! numatomic  : Atomic number
! natom      : Number of atoms
! neleca     : Number of alpha electrons
! nelecb     : Number of beta electrons
! nmo        : Number of molecular orbitals
! multi      : Spin multiplicity
! ndummyatom : Number of dummy (X) atoms
! coord      : Cartesian coordinate 
! znuc       : Atomic charge
! coordold   : Cartesian coordinate of previous optimization cycle
! charge     : Molecular charge
! enuc       : Nuclear repulsion energy
! escf       : SCF energy
! escfe      : Electronic energy
! emp2       : MP2 energy
! escsmp2    : SCS-MP2 energy
! dipole     : Dipole moment (x, y, z, total in a.u.)
! atomrad    : Atom radii
!                H       : Bohr radius
!                He - Cn : P. Pyykko, M. Atsumi, Chem. Eur. J., 15, 186-197 (2009).
      type typemol
        integer :: numatomic(mxatom), natom, neleca, nelecb, nmo, multi=1, ndummyatom
        real(8) :: coord(3,mxatom), znuc(mxatom), coordold(3,mxatom), charge=0.0D+00
        real(8) :: enuc, escf, escfe, emp2, escsmp2, dipole(4)
        real(8) :: atomrad(-9:112)=(/ &
&       1.06D+00, 1.06D+00, 1.06D+00, 1.06D+00, 1.06D+00, 1.06D+00, 1.06D+00, 1.06D+00, 1.06D+00, &
&       0.00D+00, &
&       0.53D+00, 0.46D+00, 1.33D+00, 1.02D+00, 0.85D+00, 0.75D+00, 0.71D+00, 0.63D+00, 0.64D+00, &
&       0.67D+00, 1.55D+00, 1.39D+00, 1.26D+00, 1.16D+00, 1.11D+00, 1.03D+00, 0.99D+00, 0.96D+00, &
&       1.96D+00, 1.71D+00, 1.48D+00, 1.36D+00, 1.34D+00, 1.22D+00, 1.19D+00, 1.16D+00, 1.11D+00, &
&       1.10D+00, 1.12D+00, 1.18D+00, 1.24D+00, 1.21D+00, 1.21D+00, 1.16D+00, 1.14D+00, 1.17D+00, &
&       2.10D+00, 1.85D+00, 1.63D+00, 1.54D+00, 1.47D+00, 1.38D+00, 1.28D+00, 1.25D+00, 1.25D+00, &
&       1.20D+00, 1.28D+00, 1.36D+00, 1.42D+00, 1.40D+00, 1.40D+00, 1.36D+00, 1.33D+00, 1.31D+00, &
&       2.32D+00, 1.96D+00, 1.80D+00, 1.63D+00, 1.76D+00, 1.74D+00, 1.73D+00, 1.72D+00, 1.68D+00, &
&       1.69D+00, 1.68D+00, 1.67D+00, 1.66D+00, 1.65D+00, 1.64D+00, 1.70D+00, 1.62D+00, 1.52D+00, &
&       1.46D+00, 1.37D+00, 1.31D+00, 1.29D+00, 1.22D+00, 1.23D+00, 1.24D+00, 1.33D+00, 1.44D+00, &
&       1.44D+00, 1.51D+00, 1.45D+00, 1.47D+00, 1.42D+00, 2.23D+00, 2.01D+00, 1.86D+00, 1.75D+00, &
&       1.69D+00, 1.70D+00, 1.71D+00, 1.72D+00, 1.66D+00, 1.66D+00, 1.68D+00, 1.68D+00, 1.65D+00, &
&       1.67D+00, 1.73D+00, 1.76D+00, 1.61D+00, 1.57D+00, 1.49D+00, 1.43D+00, 1.41D+00, 1.34D+00, &
&       1.29D+00, 1.28D+00, 1.21D+00, 1.22D+00/)
      end type
!
! Computational Information
!
! memmax       : Maximum memory size in words
! memused      : Used memory size in words
! memusedmax   : Maximum used memory size in words
! nwarn        : Number of warnings
! inpstd       : Unit number of input (default: 5, setting by command argument: 7)
! iout         : Unit number of standard output (default: 6, setting by command argument: 8)
! inpcopy      : Unit number of copied input file (default: 10)
! icheck       : Unit number of checkpoint file (default: 11)
! ixyz         : Unit number of xyz file (default: 12)
! iwall0       : System_clock count at the beginning of calculation
! iwall1       : System_clock count at the beginning of step
! nproc1       : Number of MPI processes of MPI_COMM_WORLD
! myrank1      : Rank number of MPI_COMM_WORLD
! mpi_comm1    : MPI_COMM_WORLD
! nproc2       : Number of MPI processes of 
! myrank2      : Rank number of new communicator
! mpi_comm2    : New communicator
! cpu0         : Cpu_time count at the beginning of calculation
! cpu1         : Cpu_time count at the beginning of step
! master       : Flag of master process of MPI_COMM_WORLD
! convergedscf : Flag of scf convergence
! convergedgeom: Flag of geometry optimization convergence
      type typecomp
        integer :: memmax=1000000000, memused=0, memusedmax=0
        integer :: inpstd=5, iout=6, inpcopy=10, icheck=11, ixyz=12
        integer :: nwarn=0
        integer :: iwall0, iwall1
        integer :: nproc1, myrank1, mpi_comm1, nproc2, myrank2, mpi_comm2
        real(8) :: cpu0, cpu1
        logical :: master=.true., convergedscf=.true., convergedgeom=.true.
      end type
end module
