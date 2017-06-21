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
!---------------------------------
  subroutine readinput(mpi_comm)
!---------------------------------
!
! Read input data and open checkpoint file if necessary
!
      use modiofile, only : input, check, maxline
      use modbasis, only : basis, spher
      use modmolecule, only : numatomic, natom, coord, znuc, charge, multi
      use modjob, only : method, runtype, scftype
      use modmemory, only : memory
      use modthresh, only : precision, cutint2, threshdiis, threshsoscf, threshqc, threshweight, &
&                           threshrho, threshdfock, threshdftao, threshover, threshatom, &
&                           threshmp2cphf
      use modguess, only : guess
      use modprint, only : iprint
      use modscf, only : scfconv, maxiter, dconv, maxdiis, maxsoscf, maxqc, maxqcdiag, maxqcdiagsub
      use modopt, only : nopt, optconv, cartesian
      use modunit, only : bohr
      use moddft, only : nrad, nleb, bqrad
      use modecp, only : ecp, flagecp
      use modmp2, only : ncore, nvfz, maxmp2diis, maxmp2iter
      use modprop, only : octupole
      implicit none
      integer,intent(in) :: mpi_comm
      integer :: myrank, ii, llen, intarray(16), info
      real(8) :: realarray(23)
      character(len=254) :: line
      character(len=16) :: chararray(9), mem=''
      logical :: logarray(5)
      namelist /job/ method, runtype, basis, scftype, memory, mem, charge, multi, ecp
      namelist /control/ precision, cutint2, spher, guess, iprint, bohr, check, threshover, &
&                        threshatom, octupole
      namelist /scf/ scfconv, maxiter, dconv, maxdiis, maxsoscf, maxqc, maxqcdiag, maxqcdiagsub, &
&                    threshdiis, threshsoscf, threshqc
      namelist /opt/ nopt, optconv, cartesian
      namelist /dft/ nrad, nleb, threshweight, threshrho, threshdfock, threshdftao, bqrad
      namelist /mp2/ ncore, nvfz, maxmp2diis, maxmp2iter, threshmp2cphf
!
      call para_comm_rank(myrank,mpi_comm)

      if(myrank == 0) then
        do ii= 1,maxline
          read(*,'(a)',end=100) line
          line=adjustl(line)
          llen=len_trim(line)
          call low2up(line,llen)
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
              call addapos(line,'PRECISION=',10)
            case('SCF')
              line="&"//trim(line)//" /"
              call addapos(line,'SCFCONV=',8)
            case('OPT','DFT','MP2')
              line="&"//trim(line)//" /"
          end select
          llen=len_trim(line)
          if(llen > 0) then
            write(input,'(a)')line(1:llen)
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
        rewind(input)
        read(input,nml=mp2,end=160,iostat=info)
160     if(info > 0) then
          write(*,'(" Error was found in mp2 line of input file!")')
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
      if(myrank == 0) then
        chararray(1)= method
        chararray(2)= runtype
        chararray(3)= basis
        chararray(4)= scftype
        chararray(5)= memory
        chararray(6)= guess
        chararray(7)= ecp
        chararray(8)= scfconv
        chararray(9)= precision
        realarray( 1)= charge
        realarray( 2)= cutint2
        realarray( 3)= dconv
        realarray( 4)= optconv
        realarray( 5)= threshdiis
        realarray( 6)= threshsoscf
        realarray( 7)= threshqc
        realarray( 8)= bqrad(1)
        realarray( 9)= bqrad(2)
        realarray(10)= bqrad(3)
        realarray(11)= bqrad(4)
        realarray(12)= bqrad(5)
        realarray(13)= bqrad(6)
        realarray(14)= bqrad(7)
        realarray(15)= bqrad(8)
        realarray(16)= bqrad(9)
        realarray(17)= threshweight
        realarray(18)= threshrho
        realarray(19)= threshdfock
        realarray(20)= threshdftao
        realarray(21)= threshover
        realarray(22)= threshatom
        realarray(23)= threshmp2cphf
        intarray( 1)= natom
        intarray( 2)= multi
        intarray( 3)= iprint
        intarray( 4)= maxiter
        intarray( 5)= maxdiis
        intarray( 6)= maxsoscf
        intarray( 7)= maxqc
        intarray( 8)= maxqcdiag
        intarray( 9)= maxqcdiagsub
        intarray(10)= nopt
        intarray(11)= nrad
        intarray(12)= nleb
        intarray(13)= ncore
        intarray(14)= nvfz
        intarray(15)= maxmp2diis
        intarray(16)= maxmp2iter
        logarray(1)= spher
        logarray(2)= bohr
        logarray(3)= flagecp
        logarray(4)= cartesian
        logarray(5)= octupole
      endif
!
      call para_bcastc(chararray,16*9,0,mpi_comm)
      call para_bcastc(check,64,0,mpi_comm)
      call para_bcastr(realarray,23,0,mpi_comm)
      call para_bcasti(intarray,16,0,mpi_comm)
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
      scfconv = chararray(8)
      precision=chararray(9)
      charge  = realarray( 1)
      cutint2 = realarray( 2)
      dconv   = realarray( 3)
      optconv = realarray( 4)
      threshdiis = realarray( 5)
      threshsoscf= realarray( 6)
      threshqc= realarray( 7)
      bqrad(1)= realarray( 8)
      bqrad(2)= realarray( 9)
      bqrad(3)= realarray(10)
      bqrad(4)= realarray(11)
      bqrad(5)= realarray(12)
      bqrad(6)= realarray(13)
      bqrad(7)= realarray(14)
      bqrad(8)= realarray(15)
      bqrad(9)= realarray(16)
      threshweight= realarray(17)
      threshrho   = realarray(18)
      threshdfock = realarray(19)
      threshdftao = realarray(20)
      threshover  = realarray(21)
      threshatom  = realarray(22)
      threshmp2cphf=realarray(23)
      multi   = intarray( 2)
      iprint  = intarray( 3)
      maxiter = intarray( 4)
      maxdiis = intarray( 5)
      maxsoscf= intarray( 6)
      maxqc   = intarray( 7)
      maxqcdiag=intarray( 8)
      maxqcdiagsub=intarray( 9)
      nopt    = intarray(10)
      nrad    = intarray(11)
      nleb    = intarray(12)
      ncore   = intarray(13)
      nvfz    = intarray(14)
      maxmp2diis= intarray(15)
      maxmp2iter= intarray(16)
      spher   = logarray(1)
      bohr    = logarray(2)
      flagecp = logarray(3)
      cartesian=logarray(4)
      octupole =logarray(5)
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
      character(len=3) :: table(-9:112)= &
&     (/'Bq9','Bq8','Bq7','Bq6','Bq5','Bq4','Bq3','Bq2','Bq ','X  ',&
&       'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ','Al ','Si ','P  ',&
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
      use modbasis, only : nshell, ex, coeffinp, locprim, locatom, mprim, mtype
      use modprint, only : iprint
      implicit none
      integer :: iatom, ishell, iloc, iprim, jatomcheck(-9:112)=0
      logical :: second
      character(len=3) :: table(-9:112)= &
&     (/'Bq9','Bq8','Bq7','Bq6','Bq5','Bq4','Bq3','Bq2','Bq ','X  ',&
&       'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ','Al ','Si ','P  ',&
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
            write(*,'(" Error! The subroutine writebasis supports up to i functions.")')
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
      integer :: iatom, ll, jprim, jloc, k, nprim, jatomcheck(-9:112)=0
      character(len=7) :: tblecp
      character(len=3) :: table(-9:112)= &
&     (/'Bq9','Bq8','Bq7','Bq6','Bq5','Bq4','Bq3','Bq2','Bq ','X  ',&
&       'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ','Al ','Si ','P  ',&
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


!----------------------------
  subroutine writecondition
!----------------------------
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
      use modthresh, only : precision
      use modprop, only : octupole
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
        write(*,'(" -------------------------------------------------------------------------")')
        write(*,'("   Job infomation")')
        write(*,'(" -------------------------------------------------------------------------")')
        write(*,'("   Runtype = ",a12,  ",  Method  = ",a12,     " ,  Basis    = ",a12)') &
&                  runtype, method, basis
        write(*,'("   Memory  =",i10, "MB ,  SCFtype = ",a12,     " ,  Precision= ",a12)') &
&                  memmax/125000, scftype, precision
        write(*,'("   Charge  = ",F11.1," ,  Multi   = ",i12,     " ,  Spher    = ",l1)') &
&                  charge, multi, spher
        write(*,'("   Bohr    = ",l1,11x,",  Guess   = ",a12,     " ,  Octupole = ",l1)') &
&                  bohr, guess, octupole
        if(runtype == 'OPTIMIZE') &
&       write(*,'("   Nopt    =  ",i10," ,  Optconv = ",1p,D12.1," ,  Cartesian= ",l1)') &
&                  nopt, optconv, cartesian
        if(check /= '') &
        write(*,'("   Check   = ",a64)') check
        write(*,'(" -------------------------------------------------------------------------")')

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
  subroutine low2up(line,llen)
!-------------------------------
!
! Convert lower charcters to upper 
!
      implicit none
      integer :: llen, ii, inum, ispace
      character(len=254),intent(inout) :: line
      character(len=254) :: linecopy
      character(len=26) :: upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      character(len=26) :: lower='abcdefghijklmnopqrstuvwxyz'
!
      linecopy(1:llen)=line(1:llen)
      do ii= 1,llen
        inum=(index(lower,line(ii:ii)))
        if(inum > 0) line(ii:ii)= upper(inum:inum)
      enddo
!
      inum= index(line(1:llen),'CHECK=')
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
      use modiofile, only : input, icheck, maxline
      use modmolecule, only : numatomic, natom, coord, znuc
      use modunit, only : tobohr, bohr
      use modparam, only : mxatom
      use modopt, only : cartesian
      use modwarn, only : nwarn
      use modjob, only : runtype
      implicit none
      integer :: ii, jj, minatomic
      real(8),parameter :: zero=0.0D+00
      character(len=16) :: cdummy
      character(len=254) :: line
      character(len=3) :: atomin(mxatom)
      character(len=3) :: table1(-9:112)= &
&     (/'BQ9','BQ8','BQ7','BQ6','BQ5','BQ4','BQ3','BQ2','BQ ','X  ',&
&       'H  ','HE ','LI ','BE ','B  ','C  ','N  ','O  ','F  ','NE ','NA ','MG ','AL ','SI ','P  ',&
&       'S  ','CL ','AR ','K  ','CA ','SC ','TI ','V  ','CR ','MN ','FE ','CO ','NI ','CU ','ZN ',&
&       'GA ','GE ','AS ','SE ','BR ','KR ','RB ','SR ','Y  ','ZR ','NB ','MO ','TC ','RU ','RH ',&
&       'PD ','AG ','CD ','IN ','SN ','SB ','TE ','I  ','XE ','CS ','BA ','LA ','CE ','PR ','ND ',&
&       'PM ','SM ','EU ','GD ','TB ','DY ','HO ','ER ','TM ','YB ','LU ','HF ','TA ','W  ','RE ',&
&       'OS ','IR ','PT ','AU ','HG ','TL ','PB ','BI ','PO ','AT ','RN ','FR ','RA ','AC ','TH ',&
&       'PA ','U  ','NP ','PU ','AM ','CM ','BK ','CF ','ES ','FM ','MD ','NO ','LR ','RF ','DB ',&
&       'SG ','BH ','HS ','MT ','UUN','UUU','UUB'/)
      character(len=3) :: table2(-9:112)= &
&     (/'-9 ','-8 ','-7 ','-6 ','-5 ','-4 ','-3 ','-2 ','-1 ','0  ',&
&       '1  ','2  ','3  ','4  ','5  ','6  ','7  ','8  ','9  ','10 ','11 ','12 ','13 ','14 ','15 ',&
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
            if(len_trim(line) == 0) exit
            read(line,*,err=9997,end=9997) atomin(ii),(coord(jj,ii),jj=1,3)
            natom= natom+1
            if(ii == mxatom) then
              write(*,'(" Error! Number of atoms exceeds mxatom=",i6,".")')mxatom
              write(*,'(" Increase mxatom in src/module.F90 and make again.",/)')
              call iabort
            endif
          enddo
100       continue
          do ii= 1,natom
            if(atomin(ii) == 'BQ1') atomin(ii)= 'BQ'
            do jj= -9,112
              if((atomin(ii) == table1(jj)).or.(atomin(ii) == table2(jj))) then
                numatomic(ii)= jj
                if(jj > 0) then
                  znuc(ii)= dble(jj)
                else
                  znuc(ii)= zero
                endif
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
            if(numatomic(ii) > 0) then
              znuc(ii)= dble(numatomic(ii))
            else
              znuc(ii)= zero
            endif
          enddo
          write(*,'(" Geometry is read from checkpoint file.")')
        endif
!
! Check dummy and ghost atoms
!
        if(runtype=='OPTIMIZE') then
          minatomic= minval(numatomic(1:natom))
          if((minatomic <= 0).and.(.not.cartesian)) then
            cartesian=.true.
            write(*,'(" Warning! Cartesian coordinate is used during geometry optimization.")')
            nwarn= nwarn+1
          endif
        endif
      endif
      return
!
9999  write(*,'(" Error! Keyword GEOM is not found.")')
      call iabort
9998  write(*,'(" Error! Geometry cannot be read from checkpoint file.")')
      call iabort
9997  write(*,'(" Error! Format of molecular geometry is incorrect.")')
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
      use modiofile, only : input, maxline
      use modbasis, only : exgen, coeffgen, locgenprim, mgenprim, mgentype, locgenshell, &
&                          ngenshell, atombasis
      implicit none
      integer :: ii, jj, iprim, ishell, ll, ielem(-9:112), nelem, kprim, numprim, natomshell
      character(len=3) :: element(-9:112)
      character(len=100) :: line
      character(len=16) :: symbol
      character(len=3) :: table(-9:112)= &
&     (/'BQ9','BQ8','BQ7','BQ6','BQ5','BQ4','BQ3','BQ2','BQ ','X  ',&
&       'H  ','HE ','LI ','BE ','B  ','C  ','N  ','O  ','F  ','NE ','NA ','MG ','AL ','SI ','P  ',&
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
        do ii= 1,maxline
          read(input,'(a)',end=9999)line
          if(line(1:5) == "BASIS") exit
          if(ii == maxline) then
            write(*,'(" Error! Keyword BASIS is not found.")')
            call iabort
          endif
        enddo
!
        do ll= -9,112
          line=''
          read(input,'(a)',end=300)line
          if(len_trim(line) == 0) exit
          element(:)=''
!
! Read elements
!
          read(line,*,end=100)(element(ii),ii=-9,112)
 100      continue
          nelem= 0
!
! Check elements
!
          do ii= -9,112
            if((element(ii) == '0').or.(element(ii) == '')) exit
            if(element(ii) == 'BQ1') element(ii)= 'BQ'
            do jj= -9,112
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
          do ii= -9,nelem-10
            locgenshell(ielem(ii))= ishell
          enddo
          do jj= 1,maxline
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
              case('I')
                mgentype(ishell)= 6
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
                read(input,*,end=9998,err=9998) exgen(iprim), coeffgen(iprim)
              enddo
            else
              do kprim= 1,numprim 
                iprim= iprim+1
                read(input,*,end=9998,err=9998) exgen(iprim), coeffgen(iprim), &
&                                                             coeffgen(iprim+numprim)
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
              do ii= -9,nelem-10
                ngenshell(ielem(ii))= natomshell
              enddo
              exit
            elseif(symbol == '') then
              write(*,'(" Error! End of basis functions is not found.")')
              call iabort
            endif
            do ii= -9,nelem-10
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
      character(len=16) :: checkversion, cdummy
!
      rewind(icheck)
      read(icheck,err=9999) checkversion
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
      if(checkversion(1:2) /= "1.") then
        read(icheck)
        read(icheck)
      endif
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
        read(icheck,end=9999)
        read(icheck,err=9998) scftype_g, natom_g, nao_g, nmo_g, nshell_g, nprim_g, neleca_g, &
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
      write(*,'(" Check the checkpoint file name.",/)')
      call iabort
 9998 write(*,'(" Error! Checkpoint file cannot be read in checkguess.")')
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
      character(len=16) :: checkversion
!
      if(master) then
        rewind(icheck)
        read(icheck) checkversion
        read(icheck)
        read(icheck)
        read(icheck)
        read(icheck)
        read(icheck) ((coord_g(jj,ii),jj=1,3),ii=1,natom)
        if(checkversion(1:2) /= "1.") then
          read(icheck)
          read(icheck)
        endif
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
! Interchange the order of basis functions in cmoa_g and cmob_g
        if(checkversion(1:2) == "1.") then
          call gcheckreorder(cmoa_g)
          if((scftype == 'UHF').and.(scftype_g == 'UHF')) call gcheckreorder(cmob_g)
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
      use modiofile, only : input, maxline
      use modecp, only : exgenecp, coeffgenecp, maxgenangecp, izgencore, mgentypeecp, &
&                        locgenecp, mgenprimecp, atomecp
      implicit none
      integer :: ii, jj, iprim, ll, ielem(-9:112), nelem, jprim, numprim, lmax, ielec, iang
      character(len=3) :: element(-9:112)
      character(len=100) :: line
      character(len=16) :: symbol
      character(len=3) :: table(-9:112)= &
&     (/'BQ9','BQ8','BQ7','BQ6','BQ5','BQ4','BQ3','BQ2','BQ ','X  ',&
&       'H  ','HE ','LI ','BE ','B  ','C  ','N  ','O  ','F  ','NE ','NA ','MG ','AL ','SI ','P  ',&
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
        do ii= 1,maxline
          read(input,*,end=9999)line
          if(line(1:3) == "ECP") exit
          if(ii == maxline) then
            write(*,'(" Error! Keyword ECP is not found.")')
            call iabort
          endif
        enddo
!
        do ll= -9,112
          line=''
          read(input,'(a)',end=300)line
          if(len_trim(line) == 0) exit
          element(:)=''
!
! Read elements
!
          read(line,*,end=100)(element(ii),ii=-9,112)
 100      nelem= 0
!
! Check elements
!
          do ii= -9,112
            if((element(ii) == '0').or.(element(ii) == '')) exit
            if(element(ii) == 'BQ1') element(ii)= 'BQ'
            do jj= -9,112
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
          do ii= -9,nelem-10
            maxgenangecp(ielem(ii))= lmax
            izgencore(ielem(ii))= ielec
          enddo
          do iang= 0,lmax
            read(input,*)line
            read(input,*,err=9998,end=9998)numprim
            do ii= -9,nelem-10
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
 200      if((symbol == 'LANL2DZ').or.(symbol == 'LANL2MB')) then
            do ii= -9,nelem-10
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
      integer :: minmo, maxmo, imin, imax, ii, jj, kk, iao, iatom
      real(8),intent(in) :: cmo(nao,nao), eigen(nmo)
      character(len=8) :: atomlabel(mxao)
      character(len=7) :: bflabel(mxao)
      character(len=3) :: table(-9:112)= &
&     (/'Bq9','Bq8','Bq7','Bq6','Bq5','Bq4','Bq3','Bq2','Bq ','X  ',&
&       'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ','Na ','Mg ','Al ','Si ','P  ',&
&       'S  ','Cl ','Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',&
&       'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru ','Rh ',&
&       'Pd ','Ag ','Cd ','In ','Sn ','Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',&
&       'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ','Ta ','W  ','Re ',&
&       'Os ','Ir ','Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',&
&       'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ','Md ','No ','Lr ','Rf ','Db ',&
&       'Sg ','Bh ','Hs ','Mt ','Uun','Uuu','Uub'/)
!
      character(len=7) :: anglabel(129)= &
&     (/'S      ','Px     ','Py     ','Pz     ','Dxx    ','Dxy    ','Dxz    ','Dyy    ', &
&       'Dyz    ','Dzz    ','D-2    ','D-1    ','D0     ','D+1    ','D+2    ','Fxxx   ', &
&       'Fxxy   ','Fxxz   ','Fxyy   ','Fxyz   ','Fxzz   ','Fyyy   ','Fyyz   ','Fyzz   ', &
&       'Fzzz   ','F-3    ','F-2    ','F-1    ','F0     ','F+1    ','F+2    ','F+3    ', &
&       'Gxxxx  ','Gxxxy  ','Gxxxz  ','Gxxyy  ','Gxxyz  ','Gxxzz  ','Gxyyy  ','Gxyyz  ', &
&       'Gxyzz  ','Gxzzz  ','Gyyyy  ','Gyyyz  ','Gyyzz  ','Gyzzz  ','Gzzzz  ','G-4    ', &
&       'G-3    ','G-2    ','G-1    ','G0     ','G+1    ','G+2    ','G+3    ','G+4    ', &
&       'Hxxxxx ','Hxxxxy ','Hxxxxz ','Hxxxyy ','Hxxxyz ','Hxxxzz ','Hxxyyy ','Hxxyyz ', &
&       'Hxxyzz ','Hxxzzz ','Hxyyyy ','Hxyyyz ','Hxyyzz ','Hxyzzz ','Hxzzzz ','Hyyyyy ', &
&       'Hyyyyz ','Hyyyzz ','Hyyzzz ','Hyzzzz ','Hzzzzz ','H-5    ','H-4    ','H-3    ', &
&       'H-2    ','H-1    ','H0     ','H+1    ','H+2    ','H+3    ','H+4    ','H+5    ', &
&       'Ixxxxxx','Ixxxxxy','Ixxxxxz','Ixxxxyy','Ixxxxyz','Ixxxxzz','Ixxxyyy','Ixxxyyz', &
&       'Ixxxyzz','Ixxxzzz','Ixxyyyy','Ixxyyyz','Ixxyyzz','Ixxyzzz','Ixxzzzz','Ixyyyyy', &
&       'Ixyyyyz','Ixyyyzz','Ixyyzzz','Ixyzzzz','Ixzzzzz','Iyyyyyy','Iyyyyyz','Iyyyyzz', &
&       'Iyyyzzz','Iyyzzzz','Iyzzzzz','Izzzzzz','I-6    ','I-5    ','I-4    ','I-3    ', &
&       'I-2    ','I-1    ','I0     ','I+1    ','I+2    ','I+3    ','I+4    ','I+5    ', &
&       'I+6    '/)
!
      if(maxval(mtype(1:nshell)) > 6) then
        if(master) write(*,'(" Sorry! This program can display MOs of up to i functions.")')
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
          case(5)
            if(spher) then
              bflabel(iao:iao+10)= anglabel(78:88)
              if(locatom(ii) /= iatom) then
                iatom= locatom(ii)
                write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
              endif
              iao= iao+11
            else
              bflabel(iao:iao+20)= anglabel(57:77)
              if(locatom(ii) /= iatom) then
                iatom= locatom(ii)
                write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
              endif
              iao= iao+21
            endif
          case(6)
            if(spher) then
              bflabel(iao:iao+12)= anglabel(117:129)
              if(locatom(ii) /= iatom) then
                iatom= locatom(ii)
                write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
              endif
              iao= iao+13
            else
              bflabel(iao:iao+27)= anglabel(89:116)
              if(locatom(ii) /= iatom) then
                iatom= locatom(ii)
                write(atomlabel(iao),'(i4,x,a3)')iatom, table(numatomic(locatom(ii)))
              endif
              iao= iao+28
            endif
        end select
      enddo
!
      minmo=max(1,neleca-15)
      maxmo=min(nmo,neleca+15)
      imin= minmo
      imax= minmo+4
      if(master) then
        do ii= 1,(maxmo-minmo-1)/5+1
          if(imax > maxmo) imax= maxmo
          write(*,*)
          write(*,'(20x,5(5x,i4,2x))')(jj,jj=imin,imax)
          write(*,'(20x,5f11.4)')(eigen(jj),jj=imin,imax)
          do kk= 1,nao
            write(*,'(i5,a8,a7,5f11.4)')kk,atomlabel(kk),bflabel(kk),(cmo(kk,jj),jj=imin,imax)
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
      use modmolecule, only : numatomic, natom, coord, nmo, neleca, nelecb, charge, multi, znuc
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
      datatype= 'znuc'
      write(icheck) datatype
      write(icheck)(znuc(ii),ii=1,natom)
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


!---------------------------------
  subroutine setcharge(mpi_comm)
!---------------------------------
!
! Set atom charge
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
              write(*,'("   Charge of Atom ",i5,"     ",f7.3)')iatom, znew
            enddo
          endif
        enddo
 100    write(*,*)
      endif
 200  continue
!
      call para_bcastr(znuc,natom,0,mpi_comm)
!
      return
end





