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
      integer,parameter :: mxprsh=30, mxprsh2=30*30
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
      real(8),parameter :: toang= 0.5291772108D+00,tobohr= 1.889726125D+00
      logical :: bohr
end


!---------------------
  module modparallel
!---------------------
#ifdef MPI
      use mpi
#endif
      implicit none
      integer :: nproc, myrank
#ifndef MPI
      integer :: MPI_COMM_WORLD, MPI_SUM
#endif
      logical :: master, parallel
!
      interface para_bcast
        subroutine para_bcastd1(buff,num,irank,comm)
          integer,intent(in) :: num, irank
          integer(4),intent(in) :: comm
          real(8) :: buff
        end subroutine para_bcastd1
        subroutine para_bcastd2(buff,num,irank,comm)
          integer,intent(in) :: num, irank
          integer(4),intent(in) :: comm
          real(8) :: buff(*)
        end subroutine para_bcastd2
        subroutine para_bcasti41(buff,num,irank,comm)
          integer,intent(in) :: num, irank
          integer(4),intent(in) :: comm
          integer(4) :: buff
        end subroutine para_bcasti41
        subroutine para_bcasti42(buff,num,irank,comm)
          integer,intent(in) :: num, irank
          integer(4),intent(in) :: comm
          integer(4) :: buff(*)
        end subroutine para_bcasti42
        subroutine para_bcasti81(buff,num,irank,comm)
          integer,intent(in) :: num, irank
          integer(4),intent(in) :: comm
          integer(8) :: buff
        end subroutine para_bcasti81
        subroutine para_bcasti82(buff,num,irank,comm)
          integer,intent(in) :: num, irank
          integer(4),intent(in) :: comm
          integer(8) :: buff(*)
        end subroutine para_bcasti82
        subroutine para_bcastc1(buff,num,irank,comm)
          integer,intent(in) :: num, irank
          integer(4),intent(in) :: comm
          character(*) :: buff
        end subroutine para_bcastc1
        subroutine para_bcastc2(buff,num,irank,comm)
          integer,intent(in) :: num, irank
          integer(4),intent(in) :: comm
          character(*) :: buff(*)
        end subroutine para_bcastc2
        subroutine para_bcastl1(buff,num,irank,comm)
          integer,intent(in) :: num, irank
          integer(4),intent(in) :: comm
          logical :: buff
        end subroutine para_bcastl1
        subroutine para_bcastl2(buff,num,irank,comm)
          integer,intent(in) :: num, irank
          integer(4),intent(in) :: comm
          logical :: buff(*)
        end subroutine para_bcastl2
      end interface para_bcast
!
      interface para_allreduce
        subroutine para_allreduced1(sbuff,rbuff,num,op,comm)
          integer,intent(in) :: num
          integer(4),intent(in) :: op, comm
          real(8),intent(in) :: sbuff
          real(8),intent(out) :: rbuff
        end subroutine para_allreduced1
        subroutine para_allreduced2(sbuff,rbuff,num,op,comm)
          integer,intent(in) :: num
          integer(4),intent(in) :: op, comm
          real(8),intent(in) :: sbuff(*)
          real(8),intent(out) :: rbuff(*)
        end subroutine para_allreduced2
        subroutine para_allreducei1(sbuff,rbuff,num,op,comm)
          integer,intent(in) :: num
          integer(4),intent(in) :: op, comm
          integer,intent(in) :: sbuff
          integer,intent(out) :: rbuff
        end subroutine para_allreducei1
        subroutine para_allreducei2(sbuff,rbuff,num,op,comm)
          integer,intent(in) :: num
          integer(4),intent(in) :: op, comm
          integer,intent(in) :: sbuff(*)
          integer,intent(out) :: rbuff(*)
        end subroutine para_allreducei2
      end interface para_allreduce
      interface para_sendrecv
        subroutine para_sendrecvd1(sbuff,nums,dest,ntags,rbuff,numr,source,ntagr,comm)
          integer,intent(in) :: nums, dest, ntags, numr, source, ntagr
          integer(4),intent(in) :: comm
          real(8),intent(in) :: sbuff
          real(8),intent(out) :: rbuff
        end subroutine para_sendrecvd1
        subroutine para_sendrecvd2(sbuff,nums,dest,ntags,rbuff,numr,source,ntagr,comm)
          integer,intent(in) :: nums, dest, ntags, numr, source, ntagr
          integer(4),intent(in) :: comm
          real(8),intent(in) :: sbuff(*)
          real(8),intent(out) :: rbuff(*)
        end subroutine para_sendrecvd2
      end interface para_sendrecv
end


!-------------------
  module modiofile
!-------------------
      implicit none
      integer,parameter :: in=10
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
! threshweight: threshold for weight at a grid point
! threshrho   : threshold for density at a grid point
! threshdfock : threshold for functional at a grid point
!
      implicit none
!ishimura
      real(8) :: threshex=30.0D+00, threshover=1.0D-06
!     real(8) :: threshex=28.0D+00, threshover=1.0D-06
      real(8) :: threshatm=2.0D-01, thresherr=0.6D+00
      real(8) :: cutint2, threshsoscf
      real(8) :: threshweight=1.0D-09, threshrho=1.0D-07, threshdfock=1.0D-06
      real(8) :: threshdftao=1.0D-04
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
      integer :: locgenshell(112), ngenshell(112)
      real(8) :: ex(mxprim), coeff(mxprim), coeffinp(mxprim)
      real(8) :: exgen(mxprim), coeffgen(mxprim)
      character(len=16) :: basis, atombasis(112)
      logical :: spher
end


!--------------------
  module modhermite
!--------------------
      implicit none
      integer :: minh(12)= (/1,2,4, 7,11,16,22,29,37,46,56,67/)
      integer :: maxh(12)= (/1,3,6,10,15,21,28,36,45,55,66,78/)
      integer :: ix(15,0:4), iy(15,0:4), iz(15,0:4)
      real(8) :: hnode(78), hweight(78)
      data ix/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
&             1,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
&             2,0,0,1,1,0,0,0,0,0,0,0,0,0,0, &
&             3,0,0,2,2,1,0,1,0,1,0,0,0,0,0, &
&             4,0,0,3,3,1,0,1,0,2,2,0,2,1,1/
      data iy/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
&             0,1,0,0,0,0,0,0,0,0,0,0,0,0,0, &
&             0,2,0,1,0,1,0,0,0,0,0,0,0,0,0, &
&             0,3,0,1,0,2,2,0,1,1,0,0,0,0,0, &
&             0,4,0,1,0,3,3,0,1,2,0,2,1,2,1/
      data iz/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
&             0,0,1,0,0,0,0,0,0,0,0,0,0,0,0, &
&             0,0,2,0,1,1,0,0,0,0,0,0,0,0,0, &
&             0,0,3,0,1,0,1,2,2,1,0,0,0,0,0, &
&             0,0,4,0,1,0,1,3,3,0,2,2,1,1,2/
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
! These valuables are for guess calculations.
! (Extended Huckel)
!
! iguess : type of initial guess
!    = 1 : calculate extended Huckel orbitals
!    = 2 : calculate from previous orbitals
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
      integer :: iguess
      integer :: nshell_v, nao_v, nprim_v, nmo_v, nshell_g, nao_g, nprim_g, nmo_g, nao_c
      integer :: locprim_g(mxshell+1), locbf_g(mxshell+1), locatom_g(mxshell)
      integer :: mprim_g(mxshell), mbf_g(mxshell), mtype_g(mxshell)
      real(8) :: ex_g(mxprim), coeff_g(mxprim), coord_g(3,mxatom)
      logical :: spshell_g, spdshell_g, spher_g
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
! These valuables are for SCF calculations.
!
      implicit none
      integer :: maxiter, maxdiis, maxsoscf
      real(8) :: dconv
      logical :: fdiff, diis, extrap
end


!----------------
  module moddft
!----------------
      implicit none
      integer :: nrad, nleb
      real(8) :: hfexchange
end


!-----------------
  module modatom
!-----------------
      implicit none
      real(8) :: atomrad(137)
      data atomrad/1.0000D+00,0.5882D+00,3.0769D+00,2.0513D+00,1.5385D+00,1.2308D+00,&
&                  1.0256D+00,0.8791D+00,0.7692D+00,0.6838D+00,4.0909D+00,3.1579D+00,&
&                  2.5714D+00,2.1687D+00,1.8750D+00,1.6514D+00,1.4754D+00,1.3333D+00,&
&                  4.1574D+00,3.4015D+00,3.0236D+00,2.6456D+00,2.5511D+00,2.6456D+00,&
&                  2.6456D+00,2.6456D+00,2.5511D+00,2.5511D+00,2.5511D+00,2.5511D+00,&
&                  2.4566D+00,2.3622D+00,2.1732D+00,2.1732D+00,2.1732D+00,1.6630D+00,&
&                  4.4409D+00,3.7795D+00,3.4015D+00,2.9291D+00,2.7401D+00,2.7401D+00,&
&                  2.5511D+00,2.4566D+00,2.5511D+00,2.6456D+00,3.0236D+00,2.9291D+00,&
&                  2.9291D+00,2.7401D+00,2.7401D+00,2.6456D+00,2.6456D+00,2.0409D+00,&
&                  4.9133D+00,4.0629D+00,3.6850D+00,3.4960D+00,3.4960D+00,3.4960D+00,&
&                  3.4960D+00,3.4960D+00,3.4960D+00,3.4015D+00,3.3070D+00,3.3070D+00,&
&                  3.3070D+00,3.3070D+00,3.3070D+00,3.3070D+00,3.3070D+00,2.9291D+00,&
&                  2.7401D+00,2.5511D+00,2.5511D+00,2.4566D+00,2.5511D+00,2.5511D+00,&
&                  2.5511D+00,2.8346D+00,3.5905D+00,3.4015D+00,3.0236D+00,3.5905D+00,&
&                  2.4000D+00,2.2677D+00,4.9133D+00,4.0629D+00,3.6850D+00,3.4015D+00,&
&                  3.3070D+00,3.3070D+00,3.3070D+00,3.3070D+00, 43*3.3D+00/
end


!----------------
  module modopt
!----------------
      implicit none
      integer :: nopt
      real(8) :: optconv
      logical :: cartesian
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
      integer :: maxgenangecp(112), izgencore(112), mgentypeecp(mxprim)
      integer :: locgenecp(0:5,112), mgenprimecp(0:5,112)
      integer :: lmf(122), lmx(581), lmy(581), lmz(581), nbfirst(0:7), nx(84), ny(84), nz(84)
      real(8) :: execp(mxprim), coeffecp(mxprim), zlm(581)
      real(8) :: exgenecp(mxprim), coeffgenecp(mxprim)
      character(lfunc) :: ecp, atomecp(112)
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

