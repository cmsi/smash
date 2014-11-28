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
! Read input data and open checkpoint file if necessary
!
      use modparallel, only : master, parallel
      use modiofile, only : input, check, maxline
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
      integer,intent(in) :: mpi_comm
      integer :: ii, ilen, intarray(9), info
      real(8) :: realarray(4)
      character(len=254) :: line
      character(len=16) :: chararray(7), mem=''
      logical :: logarray(5)
      namelist /job/ method, runtype, basis, scftype, memory, mem, charge, multi, ecp
      namelist /control/ cutint2, spher, guess, iprint, bohr, check
      namelist /scf/ diis, maxiter, dconv, maxdiis, maxsoscf
      namelist /opt/ nopt, optconv, cartesian
      namelist /dft/ nrad, nleb
!
      if(master) then
        do ii= 1,maxline
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
              call addapos(line,'CHECK=',6)
            case('SCF','OPT','DFT')
              line="&"//trim(line)//" /"
          end select
          ilen=len_trim(line)
          if(ilen > 0) then
            write(input,'(a)')line(1:ilen)
          else
            write(input,*)
          endif
          if(ii == maxline) then
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
!
! Open checkpoint file
!
      if(check /= '') call opencheckfile
!
! Read geometry data
!
      call readatom
!
      if(runtype == 'OPT') runtype='OPTIMIZE'
      if(method == 'HF') method='HARTREE-FOCK'
      if(ecp =='NONE') ecp = ''
      if(ecp /= '') flagecp=.true.
!
      if(parallel) then
        if(master) then
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
        endif
!
        call para_bcastc(chararray,16*7,0,mpi_comm)
        call para_bcastc(check,64,0,mpi_comm)
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
      use modiofile, only : check
      use modmolecule, only : natom, neleca, nelecb, charge, multi
      use modbasis, only : nshell, nao, nprim, basis, spher
      use modmemory, only : memmax
      use modjob, only : method, runtype, scftype
      use modopt, only : nopt, optconv, cartesian
      use modunit, only : bohr
      use modguess, only : guess
      implicit none
!
      if(master) then
!
! Check the numbers of electrons and basis functions
!
        if(neleca > nao) then
          write(*,'(" Error! The number of basis functions is smaller than that of electrons.")')
          call iabort
        endif
!
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
        if(check /= '') &
        write(*,'("   Guess   = ",a12)') guess
        write(*,'("   Check   = ",a64)') check
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
      integer :: ilen, ii, inum, ispace
      character(len=254),intent(inout) :: line
      character(len=254) :: linecopy
      character(len=26) :: upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      character(len=26) :: lower='abcdefghijklmnopqrstuvwxyz'
!
      linecopy(1:ilen)=line(1:ilen)
      do ii= 1,ilen
        inum=(index(lower,line(ii:ii)))
        if(inum > 0) line(ii:ii)= upper(inum:inum)
      enddo
!
      inum= index(line(1:ilen),'CHECK=')
      if(inum > 0) then
        ispace= index(line(inum:),' ')
        line(inum+6:inum+ispace-2)= linecopy(inum+6:inum+ispace-2)
      endif
      return
end


!----------------------
  subroutine readatom
!----------------------
!
! Read atomic data
!
      use modparallel, only : master
      use modiofile, only : input, icheck, check, maxline
      use modmolecule, only : numatomic, natom, coord, znuc
      use modunit, only : tobohr, bohr
      use modparam, only : mxatom
      implicit none
      integer :: ii, jj
      character(len=16) :: cdummy
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
        do ii= 1,maxline
          read(input,'(a)',end=9999)line
          if(line(1:4) == 'GEOM') exit
          if(ii == maxline) then
            write(*,'(" Error! Molecular geometry is not found.")')
            call iabort
          endif
        enddo
!
! Read data from input file
!
        if(index(line,'CHECK') == 0) then
          natom= 0
          do ii= 1,mxatom
            read(input,'(a)',end=100) line
            read(line,*,end=100) atomin(ii),(coord(jj,ii),jj=1,3)
            natom= natom+1
          enddo
100       continue
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
!
! Read data from checkpoint file
!
        else
          rewind(icheck)
          read(icheck,err=9998)
          read(icheck) cdummy, natom
          read(icheck)
          read(icheck)(numatomic(ii),ii=1,natom)
          read(icheck)
          read(icheck)((coord(jj,ii),jj=1,3),ii=1,natom)
          do ii= 1,natom
            znuc(ii)= dble(numatomic(ii))
          enddo
          write(*,'(" Geometry is read from checkpoint file.")')
        endif
      endif
      return
!
9999  write(*,'(" Error! Keyword GEOM is not found.")')
      call iabort
9998  write(*,'(" Error! Geometry cannot be read from checkpoint file.")')
      call iabort
!
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
          read(input,'(a)',end=9999)line
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
          do jj= 1,1000
            symbol= ''
            read(input,'(a)',err=200,end=200) line
            read(line,*,end=200,err=9998) symbol, numprim
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
              case default
                write(*,'(" Error! The angular momentum ",a2," is not supported.")') symbol
                call iabort
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


!---------------------------
  subroutine setcheckbasis
!---------------------------
!
! Read basis set from checkpoint file
!
      use modiofile, only : icheck
      use modbasis, only : nshell, nao, nprim, ex, coeff, locprim, locbf, locatom, &
&                          mprim, mbf, mtype
      use modparam, only : mxao, mxshell, mxprim
      implicit none
      integer :: idummy, ii
      character(len=16) :: cdummy
!
      rewind(icheck)
      read(icheck,err=9999)
      read(icheck) cdummy, idummy, nao, idummy, nshell, nprim
!
      write(*,'(" Basis set is read from checkpoint file.")')
      if(nshell+1 > mxshell) then
        write(*,'(" Error! The number of basis shells exceeds mxshell",i6,".")')mxshell
        call iabort
      endif
      if(nprim > mxprim ) then
        write(*,'(" Error! The number of primitive basis functions exceeds mxprim",i6,".")')mxprim
        call iabort
      endif
      if(nao > mxao ) then
        write(*,'(" Error! The number of basis functions exceeds mxao",i6,".")')mxao
        call iabort
      endif
!
      read(icheck)
      read(icheck)
      read(icheck)
      read(icheck)
      read(icheck)
      read(icheck) (ex(ii),ii=1,nprim)
      read(icheck)
      read(icheck) (coeff(ii),ii=1,nprim)
      read(icheck)
      read(icheck) (locprim(ii),ii=1,nshell)
      read(icheck)
      read(icheck) (locbf(ii),ii=1,nshell)
      read(icheck)
      read(icheck) (locatom(ii),ii=1,nshell)
      read(icheck)
      read(icheck) (mprim(ii),ii=1,nshell)
      read(icheck)
      read(icheck) (mbf(ii),ii=1,nshell)
      read(icheck)
      read(icheck) (mtype(ii),ii=1,nshell)
!
      locbf(nshell+1)= nao
      locprim(nshell+1)= nprim
!
      return
9999  write(*,'(" Error! Basis set cannot be read from checkpoint file.")')
      call iabort
!
end


!------------------------------------------------------------------------------------
  subroutine readcheckinfo(scftype_g,charge_g,flagecp_g,neleca_g,nelecb_g,mpi_comm)
!------------------------------------------------------------------------------------
!
! Read checkpoint information
!
      use modparallel, only : master
      use modiofile, only : icheck
      use modguess, only : nao_g, nmo_g, nshell_g, nprim_g
      use modmolecule, only : natom
      implicit none
      integer,intent(in) :: mpi_comm
      integer,intent(out) :: neleca_g, nelecb_g
      integer :: intarray(6), natom_g, idummy
      real(8),intent(out) :: charge_g
      character(len=16),intent(out) :: scftype_g
      character(len=16) :: cdummy
      logical,intent(out) :: flagecp_g
!
      if(master) then
        rewind(icheck)
        read(icheck)
        read(icheck,err=9999) scftype_g, natom_g, nao_g, nmo_g, nshell_g, nprim_g, neleca_g, &
&                             nelecb_g, cdummy, cdummy, charge_g, idummy, flagecp_g
        if(natom_g /= natom) then
          write(*,'(" Error! The numbers of atoms in checkpoint and input files are different.")')
          call iabort
        endif
        intarray(1)= nshell_g
        intarray(2)= nmo_g
        intarray(3)= neleca_g
        intarray(4)= nelecb_g
        intarray(5)= nao_g
        intarray(6)= nprim_g
      endif
      call para_bcasti(intarray,6,0,mpi_comm)
      call para_bcastc(scftype_g,16,0,mpi_comm)
      call para_bcastl(flagecp_g,1,0,mpi_comm)
      call para_bcastr(charge_g,1,0,mpi_comm)
      nshell_g= intarray(1)
      nmo_g   = intarray(2)
      neleca_g= intarray(3)
      nelecb_g= intarray(4)
      nao_g   = intarray(5)
      nprim_g = intarray(6)
!
      return
!
 9999 write(*,'(" Error! Checkpoint file cannot be read in checkguess.")')
      call iabort
end


!--------------------------------------------------------------
  subroutine readcheckguess(cmoa_g,cmob_g,scftype_g,mpi_comm)
!--------------------------------------------------------------
!
! Read guess basis functions and MOs from checkpoint file
!
      use modparallel, only : master
      use modiofile, only : icheck
      use modguess, only : locatom_g, locprim_g, locbf_g, mprim_g, mbf_g, mtype_g, &
&                          ex_g, coeff_g, nao_g, coord_g, nmo_g, nshell_g, nprim_g
      use modmolecule, only : natom
      use modjob, only : scftype
      implicit none
      integer,intent(in) :: mpi_comm
      integer :: ii, jj
      real(8),intent(out) :: cmoa_g(nao_g,nao_g), cmob_g(nao_g,nao_g)
      character(len=16),intent(in) :: scftype_g
!
      if(master) then
        read(icheck)
        read(icheck)
        read(icheck)
        read(icheck) ((coord_g(jj,ii),jj=1,3),ii=1,natom)
        read(icheck)
        read(icheck) (ex_g(ii),ii=1,nprim_g)
        read(icheck)
        read(icheck) (coeff_g(ii),ii=1,nprim_g)
        read(icheck)
        read(icheck) (locprim_g(ii),ii=1,nshell_g)
        read(icheck)
        read(icheck) (locbf_g(ii),ii=1,nshell_g)
        read(icheck)
        read(icheck) (locatom_g(ii),ii=1,nshell_g)
        read(icheck)
        read(icheck) (mprim_g(ii),ii=1,nshell_g)
        read(icheck)
        read(icheck) (mbf_g(ii),ii=1,nshell_g)
        read(icheck)
        read(icheck) (mtype_g(ii),ii=1,nshell_g)
!
        read(icheck)
        read(icheck)((cmoa_g(jj,ii),jj=1,nao_g),ii=1,nmo_g)
        if((scftype == 'UHF').and.(scftype_g == 'UHF')) then
          read(icheck)
          read(icheck)((cmob_g(jj,ii),jj=1,nao_g),ii=1,nmo_g)
        endif
      endif
!
! Broadcast guess basis functions and MOs
!
      call para_bcasti(locprim_g,nshell_g,0,mpi_comm)
      call para_bcasti(locbf_g  ,nshell_g,0,mpi_comm)
      call para_bcasti(locatom_g,nshell_g,0,mpi_comm)
      call para_bcastr(ex_g   ,nprim_g,0,mpi_comm)
      call para_bcastr(coeff_g,nprim_g,0,mpi_comm)
      call para_bcastr(coord_g,natom*3,0,mpi_comm)
      call para_bcasti(mprim_g,nshell_g,0,mpi_comm)
      call para_bcasti(mbf_g  ,nshell_g,0,mpi_comm)
      call para_bcasti(mtype_g,nshell_g,0,mpi_comm)
      call para_bcastr(cmoa_g,nao_g*nmo_g,0,mpi_comm)
      if((scftype == 'UHF').and.(scftype_g == 'UHF')) then
        call para_bcastr(cmob_g,nao_g*nmo_g,0,mpi_comm)
      endif
!
      return
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


!--------------------------------------------------
  subroutine writeeigenvalue(eigena,eigenb,itype)
!--------------------------------------------------
!
! Write eigenvalues
!
      use modparallel, only : master
      use modmolecule, only : nmo, neleca, nelecb
      use modprint, only : iprint
      integer,intent(in) :: itype
      integer :: imo
      real(8),intent(in) :: eigena(nmo), eigenb(nmo)
!
! Closed-shell
!
      if(master.and.(iprint >= 1)) then
        if(itype == 1) then
          write(*,'(1x,80("-"))')
          write(*,'("   Eigenvalues (Hartree)")')
          write(*,'(1x,80("-"))')
          write(*,'("   Alpha Occupied: ",5f12.5)')(eigena(imo),imo=1,neleca)
          write(*,'("   Alpha Virtual : ",5f12.5)')(eigena(imo),imo=neleca+1,nmo)
          write(*,'(1x,80("-"))')
!
! Open-shell
!
        elseif(itype == 2) then
          write(*,'(1x,80("-"))')
          write(*,'("   Eigenvalues (Hartree)")')
          write(*,'(1x,80("-"))')
          write(*,'("   Alpha Occupied: ",5f12.5)')(eigena(imo),imo=1,neleca)
          write(*,'("   Alpha Virtual : ",5f12.5)')(eigena(imo),imo=neleca+1,nmo)
          write(*,'("   Beta  Occupied: ",5f12.5)')(eigenb(imo),imo=1,nelecb)
          write(*,'("   Beta  Virtual : ",5f12.5)')(eigenb(imo),imo=nelecb+1,nmo)
          write(*,'(1x,80("-"))')
        endif
      endif
!
      return
end


!-----------------------------------------
  subroutine writeeigenvector(cmo,eigen)
!-----------------------------------------
!
! Write eigenvalues
!
      use modparallel, only : master
      use modmolecule, only : nmo, neleca, numatomic
      use modbasis, only : nao, nshell, mtype, spher, locatom
      use modparam, only : mxao
      implicit none
      integer :: maxmo, imin, imax, ii, jj, kk, iao, iatom
      real(8),intent(in) :: cmo(nao,nao), eigen(nmo)
      character(len=8) :: atomlabel(mxao)
      character(len=5) :: bflabel(mxao)
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
      character(len=5) :: anglabel(56)= &
&     (/'S    ','Px   ','Py   ','Pz   ','Dxx  ','Dyy  ','Dzz  ','Dxy  ','Dxz  ','Dyz  ', &
&       'D 0  ','D+1  ','D-1  ','D+2  ','D-2  ','Fxxx ','Fyyy ','Fzzz ','Fxxy ','Fxxz ', &
&       'Fxyy ','Fyyz ','Fxzz ','Fyzz ','Fxyz ','F 0  ','F+1  ','F-1  ','F+2  ','F-2  ', &
&       'F+3  ','F-3  ','Gxxxx','Gyyyy','Gzzzz','Gxxxy','Gxxxz','Gxyyy','Gyyyz','Gxzzz', &
&       'Gyzzz','Gxxyy','Gxxzz','Gyyzz','Gxxyz','Gxyyz','Gxyzz','G 0  ','G+1  ','G-1  ', &
&       'G+2  ','G-2  ','G+3  ','G-3  ','G+4  ','G-4  '/)
!
      if(maxval(mtype(1:nshell)) >= 5) then
        if(master) write(*,'(" Sorry! This program can not display MOs of h functions now.")')
        return
      endif
!
      atomlabel(1:nao)= ''
      iao= 1
      iatom= 0
      do ii= 1,nshell
        select case(mtype(ii))
          case(0)
            bflabel(iao)= anglabel(1)
            if(locatom(ii) /= iatom) then
              iatom= locatom(ii)
              write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
            endif
            iao= iao+1
          case(1)
            bflabel(iao:iao+2)= anglabel(2:4)
            if(locatom(ii) /= iatom) then
              iatom= locatom(ii)
              write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
            endif
            iao= iao+3
          case(2)
            if(spher) then
              bflabel(iao:iao+4)= anglabel(11:15)
              if(locatom(ii) /= iatom) then
                iatom= locatom(ii)
                write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
              endif
              iao= iao+5
            else
              bflabel(iao:iao+5)= anglabel(5:10)
              if(locatom(ii) /= iatom) then
                iatom= locatom(ii)
                write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
              endif
              iao= iao+6
            endif
          case(3)
            if(spher) then
              bflabel(iao:iao+6)= anglabel(26:32)
              if(locatom(ii) /= iatom) then
                iatom= locatom(ii)
                write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
              endif
              iao= iao+7
            else
              bflabel(iao:iao+9)= anglabel(16:25)
              if(locatom(ii) /= iatom) then
                iatom= locatom(ii)
                write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
              endif
              iao= iao+10
            endif
          case(4)
            if(spher) then
              bflabel(iao:iao+8)= anglabel(48:56)
              if(locatom(ii) /= iatom) then
                iatom= locatom(ii)
                write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
              endif
              iao= iao+9
            else
              bflabel(iao:iao+14)= anglabel(33:47)
              if(locatom(ii) /= iatom) then
                iatom= locatom(ii)
                write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
              endif
              iao= iao+15
            endif
        end select
      enddo
!
      maxmo=min(nmo,neleca+20)
      imin= 1
      imax= 5
      if(master) then
        do ii= 1,(maxmo-1)/5+1
          if(imax > maxmo) imax= maxmo
          write(*,*)
          write(*,'(18x,5(5x,i4,2x))')(jj,jj=imin,imax)
          write(*,'(18x,5f11.4)')(eigen(jj),jj=imin,imax)
          do kk= 1,nao
            write(*,'(i5,a8,a5,5f11.4)')kk,atomlabel(kk),bflabel(kk),(cmo(kk,jj),jj=imin,imax)
          enddo
          imin= imin+5
          imax= imax+5
        enddo
        write(*,*)
      endif
!
      return
end


!---------------------------------------------------------------------
  subroutine writecheck(cmoa,cmob,dmtrxa,dmtrxb,energymoa,energymob)
!---------------------------------------------------------------------
!
! Write checkpoint file
!
      use modiofile, only : icheck, version
      use modparallel, only : master
      use modmolecule, only : numatomic, natom, coord, nmo, neleca, nelecb, charge, multi
      use modjob, only : method, runtype, scftype
      use modbasis, only : nshell, nao, nprim, ex, coeffinp, locprim, locbf, locatom, &
&                          mprim, mbf, mtype
      use modecp, only : flagecp
      implicit none
      integer :: ii, jj
      real(8),intent(in) :: cmoa(nao,nao), dmtrxa(nao*(nao+1)/2), energymoa(nao)
      real(8),intent(in) :: cmob(nao,nao), dmtrxb(nao*(nao+1)/2), energymob(nao)
      character(len=16) :: datatype
!
      rewind(icheck)
      write(icheck) version
      write(icheck) scftype, natom, nao, nmo, nshell, nprim, neleca, nelecb,  &
&                   method, runtype, charge, multi, flagecp
!
      datatype= 'numatomic'
      write(icheck) datatype
      write(icheck) (numatomic(ii),ii=1,natom)
!
      datatype= 'coord'
      write(icheck) datatype
      write(icheck)((coord(jj,ii),jj=1,3),ii=1,natom)
!
      datatype= 'ex'
      write(icheck) datatype
      write(icheck) (ex(ii),ii=1,nprim)
!
      datatype= 'coeffinp'
      write(icheck) datatype
      write(icheck) (coeffinp(ii),ii=1,nprim)
!
      datatype= 'locprim'
      write(icheck) datatype
      write(icheck) (locprim(ii),ii=1,nshell)
!
      datatype= 'locbf'
      write(icheck) datatype
      write(icheck) (locbf(ii),ii=1,nshell)
!
      datatype= 'locatom'
      write(icheck) datatype
      write(icheck) (locatom(ii),ii=1,nshell)
!
      datatype= 'mprim'
      write(icheck) datatype
      write(icheck) (mprim(ii),ii=1,nshell)
!
      datatype= 'mbf'
      write(icheck) datatype
      write(icheck) (mbf(ii),ii=1,nshell)
!
      datatype= 'mtype'
      write(icheck) datatype
      write(icheck) (mtype(ii),ii=1,nshell)
!
      if(scftype == 'RHF') then
        datatype= 'cmo'
        write(icheck) datatype
        write(icheck)((cmoa(jj,ii),jj=1,nao),ii=1,nmo)
!
        datatype= 'dmtrx'
        write(icheck) datatype
        write(icheck) (dmtrxa(ii),ii=1,nao*(nao+1)/2)
!
        datatype= 'energymo'
        write(icheck) datatype
        write(icheck) (energymoa(ii),ii=1,nmo)
      elseif(scftype == 'UHF') then
        datatype= 'cmoa'
        write(icheck) datatype
        write(icheck)((cmoa(jj,ii),jj=1,nao),ii=1,nmo)
!
        datatype= 'cmob'
        write(icheck) datatype
        write(icheck)((cmob(jj,ii),jj=1,nao),ii=1,nmo)
!
        datatype= 'dmtrxa'
        write(icheck) datatype
        write(icheck) (dmtrxa(ii),ii=1,nao*(nao+1)/2)
!
        datatype= 'dmtrxb'
        write(icheck) datatype
        write(icheck) (dmtrxb(ii),ii=1,nao*(nao+1)/2)
!
        datatype= 'energymoa'
        write(icheck) datatype
        write(icheck) (energymoa(ii),ii=1,nmo)
!
        datatype= 'energymob'
        write(icheck) datatype
        write(icheck) (energymob(ii),ii=1,nmo)
      endif
!
      return
end


!----------------------------------
  subroutine setdetails(mpi_comm)
!----------------------------------
!
! Read and write settings (currently, charge only)
!
      use modparallel, only : master
      use modiofile, only : input, maxline
      use modmolecule, only : znuc, natom
      use modparam, only : mxatom
      implicit none
      integer,intent(in) :: mpi_comm
      integer :: ii, jj, iatom
      real(8) :: znew
      character(len=254) :: line
!
      if(master) then
        rewind(input)
        do ii= 1,maxline
          read(input,'(a)',end=200)line
          if(line(1:6) == 'CHARGE') then
            write(*,'(/," -----------------")')
            write(*,'(  "   Atomic charge")')
            write(*,'(  " -----------------")')
            write(*,'(  "   Atomic charges are set manually.")')
            do jj= 1,mxatom
              read(input,'(a)',end=100) line
              read(line,*,end=100) iatom, znew
              znuc(iatom)= znew
              write(*,'("   Charge of Atom ",i5,"     ",f6.2)')iatom, znew
            enddo
 100        write(*,*)
          endif
        enddo
      endif
 200  continue
!
      call para_bcastr(znuc,natom,0,mpi_comm)
!
      return
end





