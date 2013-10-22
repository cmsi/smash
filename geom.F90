!-----------------------
  subroutine nucenergy
!-----------------------
!
! Calculate nuclear replusion energy
!
      use molecule, only : natom, coord, charge
      use iofile, only : iout
      use thresh, only : threshatm
      use energy, only : enuc
      use warn, only : nwarn
      use procpar, only : master
      implicit none
      integer :: iatom, jatom
      real(8),parameter :: zero=0.0D+00
      real(8) :: xyz(3), rr, chrgij
!
      enuc= zero

      do iatom= 2,natom
        do jatom= 1,iatom-1
          xyz(1)= coord(1,iatom)-coord(1,jatom)
          xyz(2)= coord(2,iatom)-coord(2,jatom)
          xyz(3)= coord(3,iatom)-coord(3,jatom)
          rr= sqrt(xyz(1)*xyz(1)+xyz(2)*xyz(2)+xyz(3)*xyz(3))
          chrgij= charge(iatom)*charge(jatom)
          if(rr /= zero) then
            enuc= enuc+chrgij/rr
            if((rr <= threshatm).and.master) then
              write(iout,'("Warning! Distance of Atoms",i4," and",i4," is short!")')
              nwarn= nwarn+1
            endif       
          else
            if((chrgij /= zero).and.master) then
              write(iout,'("Error! Atoms",i4," and",i4," are the same position!")') jatom, iatom
              call iabort
            endif
          endif
        enddo
      enddo
      return
end


!--------------------------------
  subroutine nucgradient(egrad)
!--------------------------------
!
! Calculate gradinet of nuclear replusion energy
!
      use procpar
      use molecule, only : natom, coord, charge
      implicit none
      integer :: iatom, jatom, i
      real(8),intent(inout) :: egrad(3,natom)
      real(8) :: xyz(3), rr, chrgij
!
      do iatom= 2+myrank,natom,nproc
        do jatom= 1,iatom-1
          xyz(1)= coord(1,iatom)-coord(1,jatom)
          xyz(2)= coord(2,iatom)-coord(2,jatom)
          xyz(3)= coord(3,iatom)-coord(3,jatom)
          rr= xyz(1)*xyz(1)+xyz(2)*xyz(2)+xyz(3)*xyz(3)
          rr= rr*sqrt(rr)
          chrgij= charge(iatom)*charge(jatom)
          do i= 1,3
            egrad(i,iatom)= egrad(i,iatom)-xyz(i)*chrgij/rr
            egrad(i,jatom)= egrad(i,jatom)+xyz(i)*chrgij/rr
          enddo
        enddo
      enddo
      return
end
