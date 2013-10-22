!-----------------------
  subroutine readinput
!-----------------------
!
! Read input data
!
      use iofile, only : in, iout
      use procpar, only : master, parallel, MPI_COMM_WORLD
      use basis, only : func, lfunc
      use molecule, only : numatomic, natom, coord, charge, neleca, nelecb
      use control, only : method, lmethod, runtype, lruntype
      implicit none
      integer :: i, ilen
      character(254) :: line
!
      if(master) then
!ISHIMURA
        do i= 1,1000
          read(5,'(a254)',end=100) line
          ilen=len_trim(line)
          call low2up(line,ilen)
          if(ilen > 0) write(in,*)line(1:ilen)
        enddo
100     rewind(in)
        read(in,*)runtype
        read(in,*)method
        read(in,*)func
      endif
!ISHIMURA
! if(func == 'GEN') call readatomfunc
! para_bcast of ex and coeff is needed.
      call readatom
!
      neleca= 0
      do i= 1,natom
        neleca= neleca+numatomic(i)
      enddo
      if(mod(neleca,2).eq.0) then
        neleca= neleca/2
        nelecb= neleca
      else
        write(iout,'(" Number of elctrons is not even.")')
        call iabort
      endif
!
      if(runtype == "OPT") runtype="OPTIMIZE"
!
      if(parallel) then
        call para_bcast(runtype,lruntype,"C",0,MPI_COMM_WORLD)
        call para_bcast(method,lmethod,"C",0,MPI_COMM_WORLD)
        call para_bcast(func,lfunc,"C",0,MPI_COMM_WORLD)
        call para_bcast(natom,1,"I",0,MPI_COMM_WORLD)
        call para_bcast(numatomic,natom,"I",0,MPI_COMM_WORLD)
        call para_bcast(neleca,1,"I",0,MPI_COMM_WORLD)
        call para_bcast(nelecb,1,"I",0,MPI_COMM_WORLD)
        call para_bcast(coord,natom*3,"D",0,MPI_COMM_WORLD)
        call para_bcast(charge,natom,"D",0,MPI_COMM_WORLD)
      endif
      return
end


!-----------------------
  subroutine writegeom 
!-----------------------
!
! Write molecular geometry
!
      use iofile, only : iout
      use procpar, only : master
      use molecule, only : numatomic, natom, coord
      use units, only : toang
      implicit none
      integer :: i, j
      character(3) :: table(112)= &
&     (/'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ','Al ','Si ','P  ',&
&       'S  ','Cl ','Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',&
&       'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ',&
&       'Pd ','Ag ','Cd ','In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',&
&       'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ','Ta ','W  ','Re ',&
&       'Os ','Ir ','Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',&
&       'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ','Md ','No ','Lr ','Rf ','Db ',&
&       'Sg ','Bh ','Hs ','Mt ','Uun','Uuu','Uub'/)
!
      if(master) then
        write(iout,'(" ------------------------------------------------")')
        write(iout,'("          Molecular Geometry (Angstrom)")')
        write(iout,'("  Atom           X            Y            Z   ")')
        write(iout,'(" ------------------------------------------------")')
        do i= 1,natom
          write(iout,'(3x,a3,3x,3f13.6)')table(numatomic(i)),(coord(j,i)*toang,j=1,3)
        enddo
        write(iout,'(" ------------------------------------------------",/)')
      endif
      return
end


!------------------------
  subroutine writebasis
!------------------------
!
! Write basis functions
!
      use iofile, only : iout
      use procpar, only : master
      use molecule, only : numatomic, natom, coord, neleca, nelecb
      use basis, only : nshell, nao, nprim, ex, coeffinp, locprim, locatom, mprim, mtype
      use print, only : iprint
      implicit none
      integer :: iatom, ishell, iloc, iprim
      character(3) :: table(112)= &
&     (/'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ','Al ','Si ','P  ',&
&       'S  ','Cl ','Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',&
&       'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ',&
&       'Pd ','Ag ','Cd ','In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',&
&       'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ','Ta ','W  ','Re ',&
&       'Os ','Ir ','Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',&
&       'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ','Md ','No ','Lr ','Rf ','Db ',&
&       'Sg ','Bh ','Hs ','Mt ','Uun','Uuu','Uub'/)
!
! Write basis functions
!
      if(master) then
        if(iprint == 2) then
          iatom= 0
          write(iout,'(" -------------")')
          write(iout,'("   Basis set")')
          write(iout,'(" -------------")')
          do ishell= 1, nshell
            if(iatom /= locatom(ishell)) then
              write(iout,'(2x,a3,i4)')table(numatomic(locatom(ishell))),locatom(ishell)
              iatom= locatom(ishell)
            endif
            if(mtype(ishell) == 0) write(iout,'(4x,"S",i3)') mprim(ishell)
            if(mtype(ishell) == 1) write(iout,'(4x,"P",i3)') mprim(ishell)
            if(mtype(ishell) == 2) write(iout,'(4x,"D",i3)') mprim(ishell)
            if(mtype(ishell) == 3) write(iout,'(4x,"F",i3)') mprim(ishell)
            if(mtype(ishell) == 4) write(iout,'(4x,"G",i3)') mprim(ishell)
!
            if(mtype(ishell) >  4) then
              write(iout,'(" Error! The subroutine writebasis supports up to g functions.")')
              call iabort
            endif 
!
            iloc= locprim(ishell)
            do iprim= 1,mprim(ishell)
              if(mtype(ishell) == 0) write(iout,'(3x,2f15.7)') &
&                                           ex(iloc+iprim), coeffinp(iloc+iprim)
              if(mtype(ishell) == 1) write(iout,'(3x,2f15.7)') &
&                                           ex(iloc+iprim), coeffinp(iloc+iprim)
              if(mtype(ishell) == 2) write(iout,'(3x,2f15.7)') &
&                                           ex(iloc+iprim), coeffinp(iloc+iprim)
              if(mtype(ishell) == 3) write(iout,'(3x,2f15.7)') &
&                                           ex(iloc+iprim), coeffinp(iloc+iprim)
              if(mtype(ishell) == 4) write(iout,'(3x,2f15.7)') &
&                                           ex(iloc+iprim), coeffinp(iloc+iprim)
            enddo
          enddo
          write(iout,'("")')
        endif
      endif
      return
end


!---------------------------
  subroutine writecondition
!---------------------------
!
! Write computational conditions
!
      use iofile, only : iout
      use procpar, only : master
      use molecule, only : natom, neleca, nelecb
      use basis, only : nshell, nao, nprim
      implicit none
!
      if(master) then
        write(iout,'(" ------------------------------------------------")')
        write(iout,'("   Computational condition")')
        write(iout,'(" ------------------------------------------------")')
        write(iout,'("   Number of atoms                      =",i5)') natom
        write(iout,'("   Number of alpha electrons            =",i5)') neleca
        write(iout,'("   Number of beta electrons             =",i5)') nelecb
        write(iout,'("   Number of basis shells               =",i5)') nshell
        write(iout,'("   Number of basis contracted functions =",i5)') nao
        write(iout,'("   Number of basis primitive functions  =",i5)') nprim
        write(iout,'(" ------------------------------------------------",/)')
      endif
      return
end


!-------------------------------
  subroutine low2up(line,ilen)
!-------------------------------
!
! Convert lower charcters to upper 
!
      implicit none
      integer :: ilen, i, inum
      character(254) :: line
      character(26) :: upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      character(26) :: lower='abcdefghijklmnopqrstuvwxyz'
!
      do i= 1,ilen
        inum=(index(lower,line(i:i)))
        if(inum > 0) line(i:i)= upper(inum:inum)
      enddo
      return
end


!----------------------
  subroutine readatom
!----------------------
!
! Read atomic data
!
      use iofile, only : in, iout
      use procpar, only : master, parallel
      use molecule, only : numatomic, natom, coord, charge
      use units, only : tobohr, bohr
      use param, only : mxatom
      implicit none
      integer :: i, j
      character(3) :: atomin(mxatom)
      character(3) :: table1(112)= &
&     (/'H  ','HE ','LI ','BE ','B  ','C  ','N  ','O  ','F  ','NE ','NA ','MG ','AL ','SI ','P  ',&
&       'S  ','CL ','AR ','K  ','CA ','SC ','TI ','V  ','CR ','MN ','FE ','CO ','NI ','CU ','ZN ',&
&       'GA ','GE ','AS ','SE ','BR ','KR ','RB ','SR ','Y  ','ZR ','NB ','MO ','TC ','RU ','RH ',&
&       'PD ','AG ','CD ','IN ','SN ','SB ','TE ','I  ','XE ','CS ','BA ','LA ','CE ','PR ','ND ',&
&       'PM ','SM ','EU ','GD ','TB ','DY ','HO ','ER ','TM ','YB ','LU ','HF ','TA ','W  ','RE ',&
&       'OS ','IR ','PT ','AU ','HG ','TL ','PB ','BI ','PO ','AT ','RN ','FR ','RA ','AC ','TH ',&
&       'PA ','U  ','NP ','PU ','AM ','CM ','BK ','CF ','ES ','FM ','MD ','NO ','LR ','RF ','DB ',&
&       'SG ','BH ','HS ','MT ','UUN','UUU','UUB'/)
      character(3) :: table2(112)= &
&     (/'1  ','2  ','3  ','4  ','5  ','6  ','7  ','8  ','9  ','10 ','11 ','12 ','13 ','14 ','15 ',&
&       '16 ','17 ','18 ','19 ','20 ','21 ','22 ','23 ','24 ','25 ','26 ','27 ','28 ','29 ','30 ',&
&       '31 ','32 ','33 ','34 ','35 ','36 ','37 ','38 ','39 ','40 ','41 ','42 ','43 ','44 ','45 ',&
&       '46 ','47 ','48 ','49 ','50 ','51 ','52 ','53 ','54 ','55 ','56 ','57 ','58 ','59 ','60 ',&
&       '61 ','62 ','63 ','64 ','65 ','66 ','67 ','68 ','69 ','70 ','71 ','72 ','73 ','74 ','75 ',&
&       '76 ','77 ','78 ','79 ','80 ','81 ','82 ','83 ','84 ','85 ','86 ','87 ','88 ','89 ','90 ',&
&       '91 ','92 ','93 ','94 ','95 ','96 ','97 ','98 ','99 ','100','101','102','103','104','105',&
&       '106','107','108','109','110','111','112'/)
!
      if(master) then
        natom= 0
        do i= 1,mxatom
          read(in,*,end=100) atomin(i),(coord(j,i),j=1,3)
          natom= natom+1
        enddo
100     continue
        do i= 1,natom
          do j= 1,112
            if((atomin(i) == table1(j)).or.(atomin(i) == table2(j))) then
              numatomic(i)= j
              charge(i)= dble(j)
              exit
            endif
            if(j == 112) then
              write(iout,'(" Error! Reading input file.")')
              call iabort
            endif
          enddo
        enddo
!
! Change to atomic unit
!
        if(.not.bohr) then
          do i=1,natom
            do j= 1,3
              coord(j,i)= coord(j,i)*tobohr
            enddo
          enddo
        endif
      endif
!
      return
end
