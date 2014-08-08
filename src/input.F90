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
!---------------------------------
  subroutine readinput(mpi_comm)
!---------------------------------
!
! Read input data
!
      use modparallel, only : master, parallel
      use modiofile, only : input
      use modbasis, only : basis, spher
      use modmolecule, only : numatomic, natom, coord, znuc, neleca, nelecb, charge, multi
      use modjob, only : method, runtype, scftype
      use modmemory, only : memory
      use modthresh, only : cutint2
      use modguess, only : guess
      use modprint, only : iprint
      use modscf, only : diis, maxiter, dconv, maxdiis, maxsoscf
      use modopt, only : nopt, optconv, cartesian
      use modunit, only : bohr
      use moddft, only : nrad, nleb
      use modecp, only : ecp, flagecp
      implicit none
      integer(4),intent(in) :: mpi_comm
      integer :: ii, ilen, intarray(9), info
      real(8) :: realarray(4)
      character(len=254) :: line
      character(len=16) :: chararray(7), mem=''
      logical :: logarray(5)
      namelist /job/ method, runtype, basis, scftype, memory, mem, charge, multi, ecp
      namelist /control/ cutint2, spher, guess, iprint, bohr
      namelist /scf/ diis, maxiter, dconv, maxdiis, maxsoscf
      namelist /opt/ nopt, optconv, cartesian
      namelist /dft/ nrad, nleb
!
      if(master) then
        do ii= 1,100000
          read(5,'(a)',end=100) line
          line=adjustl(line)
          ilen=len_trim(line)
          call low2up(line,ilen)
          select case(line(1:3))
            case('JOB')
              line="&"//trim(line)//" /"
              call addapos(line,'RUNTYPE=',8)
              call addapos(line,'METHOD=',7)
              call addapos(line,'BASIS=',6)
              call addapos(line,'SCFTYPE=',8)
              call addapos(line,'MEMORY=',7)
              call addapos(line,'MEM=',4)
              call addapos(line,'ECP=',4)
            case('CON')
              line="&"//trim(line)//" /"
              call addapos(line,'GUESS=',6)
            case('SCF','OPT','DFT')
              line="&"//trim(line)//" /"
          end select
          ilen=len_trim(line)
          if(ilen > 0) then
            write(input,'(a)')line(1:ilen)
          else
            write(input,*)
          endif
          if(ii == 100000) then
            write(*,'(" Input file is too long in Subroutine readinput!")')
            call iabort
          endif
        enddo
!
100     rewind(input)
        read(input,nml=job,end=110,iostat=info)
110     if(info > 0) then
          write(*,'(" Error was found in job line of input file!")')
          call iabort
        endif
        if(len_trim(mem) /= 0) memory= mem
!
        rewind(input)
        read(input,nml=control,end=120,iostat=info)
120     if(info > 0) then
          write(*,'(" Error was found in control line of input file!")')
          call iabort
        endif
!
        rewind(input)
        read(input,nml=scf,end=130,iostat=info)
130     if(info > 0) then
          write(*,'(" Error was found in scf line of input file!")')
          call iabort
        endif
!
        rewind(input)
        read(input,nml=opt,end=140,iostat=info)
140     if(info > 0) then
          write(*,'(" Error was found in opt line of input file!")')
          call iabort
        endif
!
        rewind(input)
        read(input,nml=dft,end=150,iostat=info)
150     if(info > 0) then
          write(*,'(" Error was found in dft line of input file!")')
          call iabort
        endif
!
      endif
      call readatom
!
      if(runtype == 'OPT') runtype='OPTIMIZE'
      if(method == 'HF') method='HARTREE-FOCK'
      if(ecp /= '') flagecp=.true.
!
      if(parallel) then
        chararray(1)= method
        chararray(2)= runtype
        chararray(3)= basis
        chararray(4)= scftype
        chararray(5)= memory
        chararray(6)= guess
        chararray(7)= ecp
        realarray(1)= charge
        realarray(2)= cutint2
        realarray(3)= dconv
        realarray(4)= optconv
        intarray(1)= natom
        intarray(2)= multi
        intarray(3)= iprint
        intarray(4)= maxiter
        intarray(5)= maxdiis
        intarray(6)= maxsoscf
        intarray(7)= nopt
        intarray(8)= nrad
        intarray(9)= nleb
        logarray(1)= spher
        logarray(2)= diis
        logarray(3)= bohr
        logarray(4)= flagecp
        logarray(5)= cartesian
!
        call para_bcastc(chararray,16*7,0,mpi_comm)
        call para_bcastr(realarray,4,0,mpi_comm)
        call para_bcasti(intarray,9,0,mpi_comm)
        call para_bcastl(logarray,5,0,mpi_comm)
!
        natom= intarray(1)
!
        call para_bcasti(numatomic,natom,0,mpi_comm)
        call para_bcastr(coord(1,1),natom*3,0,mpi_comm)
        call para_bcastr(znuc,natom,0,mpi_comm)
!
        method  = chararray(1)
        runtype = chararray(2)
        basis   = chararray(3)
        scftype = chararray(4)
        memory  = chararray(5)
        guess   = chararray(6)
        ecp     = chararray(7)
        charge  = realarray(1)
        cutint2 = realarray(2)
        dconv   = realarray(3)
        optconv = realarray(4)
        multi   = intarray(2)
        iprint  = intarray(3)
        maxiter = intarray(4)
        maxdiis = intarray(5)
        maxsoscf= intarray(6)
        nopt    = intarray(7)
        nrad    = intarray(8)
        nleb    = intarray(9)
        spher   = logarray(1)
        diis    = logarray(2)
        bohr    = logarray(3)
        flagecp = logarray(4)
        cartesian=logarray(5)
      endif
!
      return
end


!-----------------------
  subroutine writegeom 
!-----------------------
!
! Write molecular geometry
!
      use modparallel, only : master
      use modmolecule, only : numatomic, natom, coord
      use modunit, only : toang
      implicit none
      integer :: i, j
      character(len=3) :: table(112)= &
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
        write(*,'(" ----------------------------------------------------")')
        write(*,'("          Molecular Geometry (Angstrom)")')
        write(*,'("  Atom            X             Y             Z")')
        write(*,'(" ----------------------------------------------------")')
        do i= 1,natom
          write(*,'(3x,a3,3x,3f14.7)')table(numatomic(i)),(coord(j,i)*toang,j=1,3)
        enddo
        write(*,'(" ----------------------------------------------------",/)')
      endif
      return
end


!------------------------
  subroutine writebasis
!------------------------
!
! Write basis functions
!
      use modparallel, only : master
      use modmolecule, only : numatomic
      use modbasis, only : nshell, nao, nprim, ex, coeffinp, locprim, locatom, mprim, mtype
      use modprint, only : iprint
      implicit none
      integer :: iatom, ishell, iloc, iprim, jatomcheck(112)=0
      logical :: second
      character(len=3) :: table(112)= &
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
      second=.false.
      if(master.and.(iprint >= 1)) then
        iatom= 0
        write(*,'(" -------------")')
        write(*,'("   Basis set")')
        write(*,'(" -------------")')
        do ishell= 1, nshell
          if(iatom /= locatom(ishell)) then
            if(jatomcheck(numatomic(locatom(ishell))) /= 0)cycle
            jatomcheck(numatomic(locatom(ishell)))= 1
            if(second) write(*,'("  ****")')
            second=.true.
            write(*,'(2x,a3)')table(numatomic(locatom(ishell)))
            iatom= locatom(ishell)
          endif
          if(mtype(ishell) == 0) write(*,'(4x,"S",i3)') mprim(ishell)
          if(mtype(ishell) == 1) write(*,'(4x,"P",i3)') mprim(ishell)
          if(mtype(ishell) == 2) write(*,'(4x,"D",i3)') mprim(ishell)
          if(mtype(ishell) == 3) write(*,'(4x,"F",i3)') mprim(ishell)
          if(mtype(ishell) == 4) write(*,'(4x,"G",i3)') mprim(ishell)
          if(mtype(ishell) == 5) write(*,'(4x,"H",i3)') mprim(ishell)
          if(mtype(ishell) == 6) write(*,'(4x,"I",i3)') mprim(ishell)
!
          if(mtype(ishell) >  6) then
            write(*,'(" Error! The subroutine writebasis supports up to h functions.")')
            call iabort
          endif 
!
          iloc= locprim(ishell)
          do iprim= 1,mprim(ishell)
             write(*,'(3x,f16.7,1x,f15.8)') ex(iloc+iprim), coeffinp(iloc+iprim)
          enddo
        enddo
        write(*,'("  ****")')
      endif
      return
end


!----------------------
  subroutine writeecp
!----------------------
!
! Write ECP functions
!
      use modparallel, only : master
      use modmolecule, only : numatomic, natom
      use modecp, only : maxangecp, mtypeecp, locecp, mprimecp, execp, coeffecp, izcore
      use modprint, only : iprint
      implicit none
      integer :: iatom, ll, jprim, jloc, k, nprim, jatomcheck(112)=0
      character(len=7) :: tblecp
      character(len=3) :: table(112)= &
&     (/'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ','Al ','Si ','P  ',&
&       'S  ','Cl ','Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',&
&       'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ',&
&       'Pd ','Ag ','Cd ','In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',&
&       'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ','Ta ','W  ','Re ',&
&       'Os ','Ir ','Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',&
&       'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ','Md ','No ','Lr ','Rf ','Db ',&
&       'Sg ','Bh ','Hs ','Mt ','Uun','Uuu','Uub'/)
!
! Write ECP functions
!
      if(master.and.(iprint >= 1)) then
        write(*,'(" ----------------")')
        write(*,'("   ECP function")')
        write(*,'(" ----------------")')
!
        do iatom= 1,natom
          if(maxangecp(iatom) /= -1)then
            if(jatomcheck(numatomic(iatom)) == 0) then
              jatomcheck(numatomic(iatom))= iatom
              write(*,'(2x,a3)')table(numatomic(iatom))
              tblecp=table(numatomic(iatom))
              ll= len_trim(tblecp)
              tblecp= tblecp(1:ll)//'-ECP'
              write(*,'(2x,a7,2x,i3,2x,i3)')tblecp, maxangecp(iatom), izcore(iatom)
              nprim= mprimecp(0,iatom)
              if(maxangecp(iatom) == 0) write(*,'(2x,"s-ul potential")')
              if(maxangecp(iatom) == 1) write(*,'(2x,"p-ul potential")')
              if(maxangecp(iatom) == 2) write(*,'(2x,"d-ul potential")')
              if(maxangecp(iatom) == 3) write(*,'(2x,"f-ul potential")')
              if(maxangecp(iatom) == 4) write(*,'(2x,"g-ul potential")')
              write(*,'(3x,i2)')nprim
              do jprim= 1,nprim
                jloc= jprim+locecp(0,iatom)
                write(*,'(2x,i1,2f16.7)')mtypeecp(jloc), execp(jloc), coeffecp(jloc)
              enddo
              do k= 1,maxangecp(iatom)
                nprim= mprimecp(k,iatom)
                if(k == 1) write(*,'(2x,"s-ul potential")')
                if(k == 2) write(*,'(2x,"p-ul potential")')
                if(k == 3) write(*,'(2x,"d-ul potential")')
                if(k == 4) write(*,'(2x,"f-ul potential")')
                write(*,'(3x,i2)')nprim
                do jprim= 1,nprim
                  jloc= jprim+locecp(k,iatom)
                  write(*,'(2x,i1,2f16.7)')mtypeecp(jloc), execp(jloc), coeffecp(jloc)
                enddo
              enddo
            endif
          endif
        enddo
!
        write(*,*)
      endif
!
      return
end


!---------------------------
  subroutine writecondition
!---------------------------
!
! Write computational conditions
!
      use modparallel, only : master
      use modmolecule, only : natom, neleca, nelecb, charge, multi
      use modbasis, only : nshell, nao, nprim, basis, spher
      use modmemory, only : memmax
      use modjob, only : method, runtype, scftype
      use modopt, only : nopt, optconv, cartesian
      use modunit, only : bohr
      implicit none
!
      if(master) then
        write(*,'(" --------------------------------------------------------------------------")')
        write(*,'("   Job infomation")')
        write(*,'(" --------------------------------------------------------------------------")')
        write(*,'("   Runtype = ",a12,",   Method  = ",a12,",   Basis    = ",a12)') &
&                  runtype, method, basis
        write(*,'("   Memory  =",i10,"MB ,   SCFtype = ",a12,",   Spher    = ",l1)') &
&                  memmax/125000, scftype, spher
        write(*,'("   Charge  = ",F11.1," ,   Multi   =  ",i10," ,   Bohr     = ",l1)') &
&                  charge, multi, bohr
        if(runtype == 'OPTIMIZE') &
&       write(*,'("   Nopt    =  ",i10," ,   Optconv = ",1p,D11.1," ,   Cartesian= ",l1)') &
&                  nopt, optconv, cartesian
        write(*,'(" --------------------------------------------------------------------------")')

        write(*,'(/," ------------------------------------------------")')
        write(*,'("   Computational condition")')
        write(*,'(" ------------------------------------------------")')
        write(*,'("   Number of atoms                      =",i5)') natom
        write(*,'("   Number of alpha electrons            =",i5)') neleca
        write(*,'("   Number of beta electrons             =",i5)') nelecb
        write(*,'("   Number of basis shells               =",i5)') nshell
        write(*,'("   Number of basis contracted functions =",i5)') nao
        write(*,'("   Number of basis primitive functions  =",i5)') nprim
        write(*,'(" ------------------------------------------------",/)')
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
      character(len=254) :: line
      character(len=26) :: upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      character(len=26) :: lower='abcdefghijklmnopqrstuvwxyz'
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
      use modparallel, only : master
      use modiofile, only : input
      use modmolecule, only : numatomic, natom, coord, znuc
      use modunit, only : tobohr, bohr
      use modparam, only : mxatom
      implicit none
      integer :: ii, jj
      character(len=254) :: line
      character(len=3) :: atomin(mxatom)
      character(len=3) :: table1(112)= &
&     (/'H  ','HE ','LI ','BE ','B  ','C  ','N  ','O  ','F  ','NE ','NA ','MG ','AL ','SI ','P  ',&
&       'S  ','CL ','AR ','K  ','CA ','SC ','TI ','V  ','CR ','MN ','FE ','CO ','NI ','CU ','ZN ',&
&       'GA ','GE ','AS ','SE ','BR ','KR ','RB ','SR ','Y  ','ZR ','NB ','MO ','TC ','RU ','RH ',&
&       'PD ','AG ','CD ','IN ','SN ','SB ','TE ','I  ','XE ','CS ','BA ','LA ','CE ','PR ','ND ',&
&       'PM ','SM ','EU ','GD ','TB ','DY ','HO ','ER ','TM ','YB ','LU ','HF ','TA ','W  ','RE ',&
&       'OS ','IR ','PT ','AU ','HG ','TL ','PB ','BI ','PO ','AT ','RN ','FR ','RA ','AC ','TH ',&
&       'PA ','U  ','NP ','PU ','AM ','CM ','BK ','CF ','ES ','FM ','MD ','NO ','LR ','RF ','DB ',&
&       'SG ','BH ','HS ','MT ','UUN','UUU','UUB'/)
      character(len=3) :: table2(112)= &
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
        rewind(input)
        do ii= 1,10000
          read(input,*,end=9999)line
          if(line(1:4) == "GEOM") exit
          if(ii == 10000) then
            write(*,'(" Error! Molecular geometry is not found.")')
            call iabort
          endif
        enddo
        natom= 0
        do ii= 1,mxatom
          read(input,'(a)',end=100) line
          read(line,*,end=100) atomin(ii),(coord(jj,ii),jj=1,3)
          natom= natom+1
        enddo
100     continue
        do ii= 1,natom
          do jj= 1,112
            if((atomin(ii) == table1(jj)).or.(atomin(ii) == table2(jj))) then
              numatomic(ii)= jj
              znuc(ii)= dble(jj)
              exit
            endif
            if(jj == 112) then
              write(*,'(" Error! Atom type ",a3," is not supported.")')atomin(ii)
              call iabort
            endif
          enddo
        enddo
!
! Change to atomic unit
!
        if(.not.bohr) then
          do ii= 1,natom
            do jj= 1,3
              coord(jj,ii)= coord(jj,ii)*tobohr
            enddo
          enddo
        endif
      endif
      return
!
9999  write(*,'(" Error! Keyword GEOM is not found.")')
      call iabort
      
end


!------------------------------------
  subroutine addapos(line,word,num)
!------------------------------------
!
! Add apostrophes to input lines
!
      implicit none
      integer,intent(in) :: num
      integer :: ii, next, ispace
      character(len=num),intent(in) :: word
      character(len=254),intent(inout) :: line
      character(len=254) :: line2
!
      ii= index(line,word(1:num))
      if(ii /= 0) then
        next= ii+num
        if((line(next:next) /='"').and.(line(next:next) /= "'")) then
          ispace= index(line(ii:),' ')
          line2= line(:ii+num-1)//"'"//line(ii+num:ii+ispace-2)//"'"//line(ii+ispace-1:)
          line= line2
        endif
      endif
!
      return
end


!-----------------------
  subroutine readbasis
!-----------------------
!
! Read basis set
!
      use modparallel, only : master
      use modparam, only : mxprim, mxshell
      use modiofile, only : input
      use modbasis, only : exgen, coeffgen, locgenprim, mgenprim, mgentype, locgenshell, &
&                          ngenshell, atombasis
      implicit none
      integer :: ii, jj, iprim, ishell, ll, ielem(112), nelem, kprim, numprim, natomshell
      character(len=3) :: element(112)
      character(len=100) :: line
      character(len=16) :: symbol
      character(len=3) :: table(112)= &
&     (/'H  ','HE ','LI ','BE ','B  ','C  ','N  ','O  ','F  ','NE ','NA ','MG ','AL ','SI ','P  ',&
&       'S  ','CL ','AR ','K  ','CA ','SC ','TI ','V  ','CR ','MN ','FE ','CO ','NI ','CU ','ZN ',&
&       'GA ','GE ','AS ','SE ','BR ','KR ','RB ','SR ','Y  ','ZR ','NB ','MO ','TC ','RU ','RH ',&
&       'PD ','AG ','CD ','IN ','SN ','SB ','TE ','I  ','XE ','CS ','BA ','LA ','CE ','PR ','ND ',&
&       'PM ','SM ','EU ','GD ','TB ','DY ','HO ','ER ','TM ','YB ','LU ','HF ','TA ','W  ','RE ',&
&       'OS ','IR ','PT ','AU ','HG ','TL ','PB ','BI ','PO ','AT ','RN ','FR ','RA ','AC ','TH ',&
&       'PA ','U  ','NP ','PU ','AM ','CM ','BK ','CF ','ES ','FM ','MD ','NO ','LR ','RF ','DB ',&
&       'SG ','BH ','HS ','MT ','UUN','UUU','UUB'/)
!
      atombasis(:)= ''
      ngenshell(:)= 0
      iprim= 0
      ishell= 0
!
      if(master) then
        rewind(input)
        do ii= 1,20000
          read(input,*,end=9999)line
          if(line(1:5) == "BASIS") exit
          if(ii == 20000) then
            write(*,'(" Error! Keyword BASIS is not found.")')
            call iabort
          endif
        enddo
!
        do ll= 1,112
          line=''
          read(input,'(a)',end=300)line
          if(len_trim(line) == 0) exit
          element(:)=''
!
! Read elements
!
          read(line,*,end=100)(element(ii),ii=1,112)
 100      continue
          nelem= 0
!
! Check elements
!
          do ii= 1,112
            if((element(ii) == '0').or.(element(ii) == '')) exit
            do jj= 1,112
              if(element(ii) == table(jj)) then
                ielem(ii)= jj
                nelem= nelem+1
                exit
              endif
              if(jj == 112) then
                write(*,'(" Error! This program does not support Atom",a3,".")')element(ii)
                call iabort
              endif
            enddo
          enddo
!
! Read basis functions
!
          natomshell= 0
          do ii= 1,nelem
            locgenshell(ielem(ii))= ishell
          enddo
          do jj= 1,400
            symbol= ''
            read(input,'(a)',err=200,end=200) line
            read(line,*,end=200,err=9998) symbol,numprim
            ishell= ishell+1
            natomshell= natomshell+1
            locgenprim(ishell)= iprim
            mgenprim(ishell)= numprim
            select case(symbol)
              case('S')
                mgentype(ishell)= 0
              case('P')
                mgentype(ishell)= 1
              case('D')
                mgentype(ishell)= 2
              case('F')
                mgentype(ishell)= 3
              case('G')
                mgentype(ishell)= 4
              case('H')
                mgentype(ishell)= 5
              case('SP')
                mgentype(ishell)  = 0
                mgentype(ishell+1)= 1
            end select
            if(symbol /= 'SP') then
              do kprim= 1,numprim 
                iprim= iprim+1
                read(input,*,err=9998) exgen(iprim), coeffgen(iprim)
              enddo
            else
              do kprim= 1,numprim 
                iprim= iprim+1
                read(input,*,err=9998) exgen(iprim), coeffgen(iprim), coeffgen(iprim+numprim)
                exgen(iprim+numprim)= exgen(iprim)
              enddo
              ishell= ishell+1
              locgenprim(ishell)= iprim
              mgenprim(ishell)= numprim
              iprim= iprim+numprim
              natomshell= natomshell+1
            endif
            cycle
!
 200        if(symbol(1:2) == '**') then
              do ii= 1,nelem
                ngenshell(ielem(ii))= natomshell
              enddo
              exit
            elseif(symbol == '') then
              write(*,'(" Error! End of basis functions is not found.")')
              call iabort
            endif
            do ii= 1,nelem
              atombasis(ielem(ii))= symbol
            enddo
          enddo
!
        enddo
 300    continue
      endif
      return
!
9999  write(*,'(" Error! Keyword BASIS is not found.")')
      call iabort
9998  write(*,'(" Error! Format of basis functions is incorrect.")')
      call iabort
end


!---------------------
  subroutine readecp
!---------------------
!
! Read basis set
!
      use modparallel, only : master
      use modparam, only : mxprim, mxshell
      use modiofile, only : input
      use modecp, only : exgenecp, coeffgenecp, maxgenangecp, izgencore, mgentypeecp, &
&                        locgenecp, mgenprimecp, atomecp
      implicit none
      integer :: ii, jj, iprim, ll, ielem(112), nelem, jprim, numprim, lmax, ielec, iang
      character(len=3) :: element(112)
      character(len=100) :: line
      character(len=16) :: symbol
      character(len=3) :: table(112)= &
&     (/'H  ','HE ','LI ','BE ','B  ','C  ','N  ','O  ','F  ','NE ','NA ','MG ','AL ','SI ','P  ',&
&       'S  ','CL ','AR ','K  ','CA ','SC ','TI ','V  ','CR ','MN ','FE ','CO ','NI ','CU ','ZN ',&
&       'GA ','GE ','AS ','SE ','BR ','KR ','RB ','SR ','Y  ','ZR ','NB ','MO ','TC ','RU ','RH ',&
&       'PD ','AG ','CD ','IN ','SN ','SB ','TE ','I  ','XE ','CS ','BA ','LA ','CE ','PR ','ND ',&
&       'PM ','SM ','EU ','GD ','TB ','DY ','HO ','ER ','TM ','YB ','LU ','HF ','TA ','W  ','RE ',&
&       'OS ','IR ','PT ','AU ','HG ','TL ','PB ','BI ','PO ','AT ','RN ','FR ','RA ','AC ','TH ',&
&       'PA ','U  ','NP ','PU ','AM ','CM ','BK ','CF ','ES ','FM ','MD ','NO ','LR ','RF ','DB ',&
&       'SG ','BH ','HS ','MT ','UUN','UUU','UUB'/)
!
      atomecp(:)= ''
      maxgenangecp(:)= -1
      iprim= 0
!
      if(master) then
        rewind(input)
        do ii= 1,20000
          read(input,*,end=9999)line
          if(line(1:3) == "ECP") exit
          if(ii == 20000) then
            write(*,'(" Error! Keyword ECP is not found.")')
            call iabort
          endif
        enddo
!
        do ll= 1,112
          line=''
          read(input,'(a)',end=300)line
          if(len_trim(line) == 0) exit
          element(:)=''
!
! Read elements
!
          read(line,*,end=100)(element(ii),ii=1,112)
 100      nelem= 0
!
! Check elements
!
          do ii= 1,112
            if((element(ii) == '0').or.(element(ii) == '')) exit
            do jj= 1,112
              if(element(ii) == table(jj)) then
                ielem(ii)= jj
                nelem= nelem+1
                exit
              endif
              if(jj == 112) then
                write(*,'(" Error! This program does not support Atom",a3,".")')element(ii)
                call iabort
              endif
            enddo
          enddo
          if(nelem == 0) then
            write(*,'(" Error! Atom type is not found during reading ECP functions.")')
            call iabort
          endif
!
! Read ECP functions
!
          symbol= ''
          read(input,'(a)',err=200,end=200) line
          read(line,*,end=200) symbol,lmax,ielec
          do ii= 1,nelem
            maxgenangecp(ielem(ii))= lmax
            izgencore(ielem(ii))= ielec
          enddo
          do iang= 0,lmax
            read(input,*)line
            read(input,*,err=9998,end=9998)numprim
            do ii= 1,nelem
              locgenecp(iang,ielem(ii))= iprim
              mgenprimecp(iang,ielem(ii))= numprim
            enddo
            do jprim=1,numprim
              iprim= iprim+1
              read(input,*,err=9998,end=9998)mgentypeecp(iprim),exgenecp(iprim),coeffgenecp(iprim)
            enddo
          enddo
          cycle
!
 200      if(symbol == 'LANL2DZ') then
            do ii= 1,nelem
              atomecp(ielem(ii))= symbol
            enddo
            cycle
          elseif(symbol == '') then
            write(*,'(" Error! ECP functions are not found.")')
            call iabort
          else
            write(*,'(" Error! This program does not support ECP function ",a10,".")')symbol
            call iabort
          endif
        enddo
!
 300    continue
      endif
      return
!
9998  write(*,'(" Error! During reading ECP functions from input file.")')
      call iabort
9999  write(*,'(" Error! Keyword ECP is not found.")')
      call iabort
end


