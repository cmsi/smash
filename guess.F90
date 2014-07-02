!-----------------------------------
  subroutine guessmo(cmo,overinv)
!-----------------------------------
!
! Initial guess generation
!   iguess = 1 : extended huckel MOs
!   iguess = 2 : MOs from previous ones
!
! In  : overinv (overlap integral inverse matrix)
! Out : cmo     (initial guess orbitals)
!
      use modguess, only : iguess
      use modbasis, only : nao
      implicit none
      real(8),intent(inout):: cmo(nao*nao), overinv(nao*nao)
!
      if(iguess == 1) then
        call huckelguess(cmo,overinv)
      elseif(iguess == 2) then
        call updatemo(cmo,overinv)
      else
        write(*,'(" Error! This program supports only iguess= 1 or 2 in guessmo.")')
        call iabort
      endif
      return
end


!---------------------------------------
  subroutine huckelguess(cmo,overinv)
!---------------------------------------
!
! Initial guess calculation
!
! In  : overinv (overlap integral inverse matrix)
! Out : cmo     (initial guess orbitals)
!
      use modguess, only : nao_g, spher_g, coord_g, nmo_g
      use modbasis, only : nao
      use modmolecule, only : coord, natom
      implicit none
      integer :: i, j, nao_g2
      real(8),intent(in):: overinv(nao*nao)
      real(8),intent(out):: cmo(nao*nao)
      real(8),allocatable :: hmo(:), overlap(:,:), work1(:), work2(:), eigen(:)
!
! Set basis functions
!
      spher_g=.true.
      call setbasis_g
!
! Set coordinate
!
      do i= 1,natom
        do j= 1,3
          coord_g(j,i)=coord(j,i)
        enddo
      enddo
!
! Set required arrays
!
      nao_g2= nao_g*nao_g
      call memset(3*nao_g2+2*nao*nao_g+nao_g)
      allocate(hmo(nao_g2),overlap(nao*nao_g,2),work1(nao_g2),work2(nao_g2),eigen(nao_g))
!
! Calculate Extended Huckel method
!
      call calchuckelg(hmo,work1,work2,eigen)
!
! Calculate overlap integrals between input basis and Huckel basis
!
      call calcover2(overlap(1,1),overlap(1,2))
!
! Project orbitals from Huckel to SCF
!
      call projectmo(cmo,overinv,overlap,hmo,work1,work2,eigen)
      deallocate(hmo,overlap,work1,work2,eigen)
      call memunset(3*nao_g2+2*nao*nao_g+nao_g)
      return
end

!-------------------------------------------------
  subroutine calchuckelg(hmo,huckel,ortho,eigen)
!-------------------------------------------------
!
! Driver of extended Huckel calculation for guess generation
!
! Out : hmo (extended Huckel orbitals)
!
      use modguess, only : nao_g, nmo_g
      implicit none
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(out) :: hmo(nao_g*nao_g), huckel(nao_g*nao_g)
      real(8),intent(out) :: ortho(nao_g*nao_g), eigen(nao_g)
!
! Calculate overlap integrals
! (guess basis)x(guess basis)
!
      call calcover1(hmo)
!
! Set ionization potentials
!
      call huckelip(eigen)
!
! Form extended Huckel matrix
!
      call formhuckel(huckel,hmo,eigen)
!
! Calculate canonicalization matrix
!
      call mtrxcanon(ortho,hmo,eigen,nao_g,nmo_g)
!
! Canonicalize extended Huckel matrix
!
!ishimura
!     call canonicalize(huckel,ortho,hmo,nao_g,nmo_g)
!
! Diagonalize canonicalized matrix
!
      call diag('V','U',nmo_g,huckel,nao_g,eigen)
!
! Backtransform to AO basis
!
!ishimura
      call dcopy(nmo_g*nao_g,huckel,1,hmo,1)
!     call dgemm('N','N',nao_g,nmo_g,nmo_g,one,ortho,nao_g,huckel,nao_g,zero,hmo,nao_g)
      return
end


!----------------------------------------------------------------------
  subroutine projectmo(cmo,overinv,overlap,hmo,work1,work2,eigen)
!----------------------------------------------------------------------
!
! Project orbitals from Huckel to SCF
!    C1= S11^-1 * S12 * C2 [C2t * S12t * S11^-1 * S12 * C2]^-1/2
!
! In  :  overinv (overlap integral inverse matrix of SCF basis set)
! Inout: overlap (overlap integral of guess and SCF basis sets)
!        hmo (extended Huckel orbitals)
! Out :  cmo (initial guess orbitals)
!
      use modguess, only : nao_g, nmo_g
      use modbasis, only : nao
      use modmolecule, only : neleca, nmo
      implicit none
      integer :: i, j
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: overinv(nao,nao)
      real(8),intent(inout) :: overlap(nao,nao_g), hmo(nao_g,nao_g)
      real(8),intent(out) :: cmo(nao,nao), work1(nao_g,nao_g), work2(nao_g,nao_g)
      real(8),intent(out) :: eigen(nao_g)
      real(8) :: eigeninv
!
! Calculate S12*C2
!
      call dgemm('N','N',nao,nmo_g,nao_g,one,overlap,nao,hmo,nao_g,zero,cmo,nao)
!
! Calculate S11^-1*S12*C2
!
      call dsymm('L','U',nao,nmo_g,one,overinv,nao,cmo,nao,zero,overlap,nao)
!
! Calculate C2t*S12t*S11^-1*S12*C2
!
      call dgemm('T','N',nmo_g,nmo_g,nao,one,cmo,nao,overlap,nao,zero,work1,nao_g)
!
! Calculate (C2t*S12t*S11^-1*S12*C2)^-1/2
!
      call diag('V','U',nmo_g,work1,nao_g,eigen)
!$OMP parallel do private(eigeninv)
      do i= 1,nmo_g
        eigeninv= one/sqrt(eigen(i))
        do j= 1,nmo_g
          work2(j,i)= work1(j,i)*eigeninv
        enddo
      enddo
!$OMP end parallel do
      call dgemm('N','T',nmo_g,nmo_g,nmo_g,one,work1,nao_g,work2,nao_g,zero,hmo,nao_g)
!
! Calculate C1
!
      call dgemm('N','N',nao,nmo_g,nmo_g,one,overlap,nao,hmo,nao_g,zero,cmo,nao)
!
      return
end


!------------------------------
  subroutine huckelip(energy)
!------------------------------
!
! Set ionization potentials
!
! Out : energy (ionization potential)
!
      use modmolecule, only : natom, numatomic
      use modguess, only : nao_g, spher_g
      use modecp, only : flagecp, izcore
      implicit none
      integer :: iao, iatom, i
      real(8),parameter :: one=1.0D+00
      real(8),intent(out) :: energy(nao_g)
      real(8) :: row1(2)=(/-5.000D-01,-9.180D-01/)
      real(8) :: row2c(3:10)
      real(8) :: row2v(2,3:10)
      real(8) :: row3c(3,11:18)
      real(8) :: row3v(2,11:18)
      real(8) :: row4c(5,19:36)
      real(8) :: row4v(3,19:36)
      real(8) :: row5c(8,37:54)
      real(8) :: row5v(3,37:54)
!ishimura
      real(8) :: row6v(4,79:79)
      data row2c/-2.48D+00,-4.73D+00,-7.70D+00,-1.13D+01,-1.56D+01,-2.07D+01,-2.64D+01,-3.28D+01/
      data row2v/-1.960D-01, 0.000D+00,-3.090D-01, 0.000D+00,-4.950D-01,-3.100D-01,&
&                -7.060D-01,-4.330D-01,-9.450D-01,-5.680D-01,-1.244D+00,-6.320D-01,&
&                -1.573D+00,-7.300D-01,-1.930D+00,-8.500D-01/
      data row3c/-4.05D+01,-2.80D+00,-1.52D+00,-4.90D+01,-3.77D+00,-2.28D+00,-5.85D+01,-4.91D+00,&
&                -3.22D+00,-6.88D+01,-6.16D+00,-4.26D+00,-8.00D+01,-7.51D+00,-5.40D+00,-9.20D+01,&
&                -9.00D+00,-6.68D+00,-1.04D+02,-1.06D+01,-8.07D+00,-1.186D+02,-1.23D+01,-9.57D+00/
      data row3v/-1.820D-01, 0.000D+00,-2.530D-01, 0.000D+00,-3.930D-01,-2.100D-01,&
&                -5.400D-01,-2.970D-01,-6.960D-01,-3.920D-01,-8.800D-01,-4.370D-01,&
&                -1.073D+00,-5.060D-01,-1.278D+00,-5.910D-01/
      data row4c/-1.335D+02,-1.450D+01,-1.150D+01,-1.750D+00,-9.500D-01,&
                 -1.494D+02,-1.680D+01,-1.360D+01,-2.240D+00,-1.340D+00,&
                 -1.659D+02,-1.908D+01,-1.567D+01,-2.570D+00,-1.575D+00,&
                 -1.833D+02,-2.142D+01,-1.779D+01,-2.874D+00,-1.795D+00,&
                 -2.013D+02,-2.370D+01,-1.980D+01,-2.990D+00,-1.840D+00,&
                 -2.204D+02,-2.620D+01,-2.210D+01,-3.290D+00,-2.050D+00,&
                 -2.404D+02,-2.890D+01,-2.460D+01,-3.620D+00,-2.300D+00,&
                 -2.612D+02,-3.170D+01,-2.720D+01,-3.960D+00,-2.550D+00,&
                 -2.829D+02,-3.460D+01,-2.990D+01,-4.300D+00,-2.800D+00,&
                 -3.054D+02,-3.770D+01,-3.270D+01,-4.650D+00,-3.060D+00,&
                 -3.288D+02,-4.080D+01,-3.560D+01,-5.010D+00,-3.320D+00,&
                 -3.533D+02,-4.440D+01,-3.890D+01,-5.630D+00,-3.840D+00,&
                 -3.788D+02,-4.820D+01,-4.250D+01,-6.400D+00,-4.480D+00,&
                 -4.052D+02,-5.210D+01,-4.620D+01,-7.190D+00,-5.170D+00,&
                 -4.326D+02,-5.630D+01,-5.020D+01,-8.030D+00,-5.880D+00,&
                 -4.609D+02,-6.070D+01,-5.430D+01,-8.930D+00,-6.660D+00,&
                 -4.901D+02,-6.520D+01,-5.860D+01,-9.870D+00,-7.480D+00,&
                 -5.202D+02,-6.990D+01,-6.300D+01,-1.080D+01,-8.330D+00/
      data row4v/-1.470D-01, 0.000D+00, 0.000D+00, -1.960D-01, 0.000D+00, 0.000D+00,&
                 -4.240D-01,-2.080D-01,-1.193D+00, -2.100D-01,-1.000D-01,-3.430D-01,&
                 -2.200D-01,-1.000D-01,-4.410D-01, -2.140D-01,-1.000D-01,-3.210D-01,&
                 -2.220D-01,-1.000D-01,-3.730D-01, -2.270D-01,-1.000D-01,-3.830D-01,&
                 -2.300D-01,-1.000D-01,-4.060D-01, -2.330D-01,-1.000D-01,-4.340D-01,&
                 -2.360D-01,-1.000D-01,-4.570D-01, -2.380D-01,-1.000D-01,-4.910D-01,&
                 -2.930D-01,-1.000D-01,-7.830D-01, -5.530D-01,-2.870D-01,-1.635D+00,&
                 -6.860D-01,-3.690D-01,-2.113D+00, -8.380D-01,-4.030D-01,-2.650D+00,&
                 -9.930D-01,-4.570D-01,-3.220D+00, -1.153D+00,-5.240D-01,-3.825D+00/
! The order is 1S,2S,2P,3S,3P,4S,4P,3D
      data row5c/ &
&       -5.510D+02,-7.500D+01,-6.790D+01,-1.210D+01,-9.500D+00,-1.520D+00,-8.100D-01,-4.700D+00,&
&       -5.837D+02,-8.040D+01,-7.300D+01,-1.350D+01,-1.070D+01,-1.900D+00,-1.100D+00,-5.700D+00,&
&       -6.168D+02,-8.581D+01,-7.816D+01,-1.476D+01,-1.185D+01,-2.168D+00,-1.300D+00,-6.599D+00,&
&       -6.507D+02,-9.138D+01,-8.348D+01,-1.606D+01,-1.302D+01,-2.418D+00,-1.487D+00,-7.515D+00,&
&       -6.854D+02,-9.700D+01,-8.880D+01,-1.720D+01,-1.400D+01,-2.530D+00,-1.550D+00,-8.300D+00,&
&       -7.212D+02,-1.029D+02,-9.450D+01,-1.860D+01,-1.530D+01,-2.760D+00,-1.720D+00,-9.300D+00,&
&       -7.579D+02,-1.089D+02,-1.002D+02,-2.000D+01,-1.660D+01,-3.000D+00,-1.910D+00,-1.030D+01,&
&       -7.955D+02,-1.152D+02,-1.062D+02,-2.140D+01,-1.780D+01,-3.260D+00,-2.100D+00,-1.130D+01,&
&       -8.340D+02,-1.216D+02,-1.124D+02,-2.290D+01,-1.920D+01,-3.500D+00,-2.290D+00,-1.240D+01,&
&       -8.735D+02,-1.281D+02,-1.187D+02,-2.440D+01,-2.050D+01,-3.750D+00,-2.480D+00,-1.350D+01,&
&       -9.138D+02,-1.349D+02,-1.252D+02,-2.590D+01,-2.190D+01,-4.000D+00,-2.680D+00,-1.470D+01,&
&       -9.554D+02,-1.421D+02,-1.321D+02,-2.770D+01,-2.360D+01,-4.450D+00,-3.050D+00,-1.610D+01,&
&       -9.978D+02,-1.494D+02,-1.392D+02,-2.960D+01,-2.540D+01,-4.980D+00,-3.510D+00,-1.760D+01,&
&       -1.041D+03,-1.570D+02,-1.465D+02,-3.160D+01,-2.720D+01,-5.510D+00,-3.970D+00,-1.920D+01,&
&       -1.086D+03,-1.648D+02,-1.540D+02,-3.360D+01,-1.920D+01,-6.060D+00,-4.450D+00,-2.080D+01,&
&       -1.131D+03,-1.728D+02,-1.617D+02,-3.580D+01,-3.110D+01,-6.650D+00,-4.950D+00,-2.250D+01,&
&       -1.177D+03,-1.809D+02,-1.697D+02,-3.790D+01,-3.310D+01,-7.240D+00,-5.470D+00,-2.430D+01,&
&       -1.224D+03,-1.893D+02,-1.778D+02,-4.020D+01,-3.520D+01,-7.860D+00,-6.010D+00,-2.610D+01/
! The order is 5S,5P,4D
      data row5v/-1.380D-01, 0.000D+00, 0.000D+00,-1.780D-01, 0.000D+00, 0.000D+00,&
&                -1.958D-01,-1.000D-01,-2.499D-01,-2.070D-01,-1.000D-01,-3.365D-01,&
&                -2.140D-01,-1.000D-01,-2.990D-01,-2.220D-01,-1.000D-01,-3.570D-01,&
&                -2.220D-01,-1.000D-01,-3.770D-01,-2.220D-01,-1.000D-01,-4.120D-01,&
&                -2.200D-01,-1.000D-01,-4.510D-01,-2.200D-01,-1.000D-01,-4.880D-01,&
&                -2.200D-01,-1.000D-01,-5.370D-01,-2.650D-01,-1.000D-01,-7.630D-01,&
&                -3.720D-01,-1.970D-01,-1.063D+00,-4.760D-01,-2.650D-01,-1.369D+00,&
&                -5.820D-01,-3.350D-01,-1.688D+00,-7.010D-01,-3.600D-01,-2.038D+00,&
&                -8.210D-01,-4.030D-01,-2.401D+00,-9.440D-01,-4.570D-01,-2.778D+00/
!ishimura
      data row6v/-8.657D+00,-7.618D+00,-5.069D+00,-1.0420D+00/
!
      iao= 0
!
! Set valence ionization potentials
!
      do iatom= 1,natom
        select case(numatomic(iatom))
! H - He
          case(1:2)
            iao= iao+1
            energy(iao)= row1(numatomic(iatom))
! Li - Ne
          case(3:10)
            iao= iao+1
            energy(iao)= row2v(1,numatomic(iatom))
            do i= 1,3
              iao= iao+1
              energy(iao)= row2v(2,numatomic(iatom))
            enddo
! Na - Ar
          case(11:18)
            iao= iao+1
            energy(iao)= row3v(1,numatomic(iatom))
            do i= 1,3
              iao= iao+1
              energy(iao)= row3v(2,numatomic(iatom))
            enddo
! K  - Ca, Ga - Kr
          case(19:20,31:36)
            iao= iao+1
            energy(iao)= row4v(1,numatomic(iatom))
            do i= 1,3
              iao= iao+1
              energy(iao)= row4v(2,numatomic(iatom))
            enddo
! Sc - Zn
          case(21:30)
            iao= iao+1
            energy(iao)= row4v(1,numatomic(iatom))
            do i= 1,3
              iao= iao+1
              energy(iao)= row4v(2,numatomic(iatom))
            enddo
            do i= 1,5
              iao= iao+1
              energy(iao)= row4v(3,numatomic(iatom))
            enddo
! Rb - Sr, In - Xe
          case(37:38,49:54)
            iao= iao+1
            energy(iao)= row5v(1,numatomic(iatom))
            do i= 1,3
              iao= iao+1
              energy(iao)= row5v(2,numatomic(iatom))
            enddo
! Y  - Cd
          case(39:48)
            iao= iao+1
            energy(iao)= row5v(1,numatomic(iatom))
            do i= 1,3
              iao= iao+1
              energy(iao)= row5v(2,numatomic(iatom))
            enddo
            do i= 1,5
              iao= iao+1
              energy(iao)= row5v(3,numatomic(iatom))
            enddo
! Au
          case(79)
            do i= 1,5
              iao= iao+1
              energy(iao)= row6v(3,numatomic(iatom))
            enddo
            iao= iao+1
            energy(iao)= row6v(4,numatomic(iatom))
          case default
            write(*,'(" Error! This program supports up to Xe in huckelip.")')
            call iabort
        end select
      enddo
!
! Set core ionization potentials
!
      do iatom= 1,natom
        select case(numatomic(iatom))
! Li - Ne
          case(1:2)
          case(3:10)
! 1s
            if(flagecp) then
              if(izcore(iatom) < 2) then
                iao= iao+1
                energy(iao)= row2c(numatomic(iatom))
              endif
            else
              iao= iao+1
              energy(iao)= row2c(numatomic(iatom))
            endif
! Na - Ar
          case(11:18)
            if(flagecp) then
! 1s
              if(izcore(iatom) < 2) then
                iao= iao+1
                energy(iao)=row3c(1,numatomic(iatom))
              endif
! 2s+2p
              if(izcore(iatom) < 10) then
                iao= iao+1
                energy(iao)=row3c(2,numatomic(iatom))
                do i= 1,3
                  iao= iao+1
                  energy(iao)=row3c(3,numatomic(iatom))
                enddo
              endif
            else
! 1s
              iao= iao+1
              energy(iao)= row3c(1,numatomic(iatom))
! 2s+2p
              iao= iao+1
              energy(iao)= row3c(2,numatomic(iatom))
              do i= 1,3
                iao= iao+1
                energy(iao)= row3c(3,numatomic(iatom))
              enddo
            endif
! K  - Kr
          case(19:36)
            if(flagecp) then
! 1s
              if(izcore(iatom) < 2) then
                iao= iao+1
                energy(iao)= row4c(1,numatomic(iatom))
              endif
! 2s+2p
              if(izcore(iatom) < 10) then
                iao= iao+1
                energy(iao)= row4c(2,numatomic(iatom))
                do i= 1,3
                  iao= iao+1
                  energy(iao)= row4c(3,numatomic(iatom))
                enddo
              endif
! 3s+3p
              if(izcore(iatom) < 18) then
                iao= iao+1
                energy(iao)= row4c(4,numatomic(iatom))
                do i= 1,3
                  iao= iao+1
                  energy(iao)= row4c(5,numatomic(iatom))
                enddo
              endif
! 3d
              if(numatomic(iatom) >= 31) then
                if(izcore(iatom) < 28) then
                  do i= 1,5
                    iao= iao+1
                    energy(iao)= row4v(3,numatomic(iatom))
                  enddo
                endif
              endif
            else
! 1s
              iao= iao+1
              energy(iao)= row4c(1,numatomic(iatom))
! 2s+2p
              iao= iao+1
              energy(iao)= row4c(2,numatomic(iatom))
              do i= 1,3
                iao= iao+1
                energy(iao)= row4c(3,numatomic(iatom))
              enddo
! 3s+3p
              iao= iao+1
              energy(iao)= row4c(4,numatomic(iatom))
              do i= 1,3
                iao= iao+1
                energy(iao)= row4c(5,numatomic(iatom))
              enddo
! 3d
              if(numatomic(iatom) >= 31) then
                do i= 1,5
                  iao= iao+1
                  energy(iao)= row4v(3,numatomic(iatom))
                enddo
              endif
            endif
! Rb - Xe
          case(37:54)
            if(flagecp) then
! 1s
              if(izcore(iatom) < 2) then
                iao= iao+1
                energy(iao)= row5c(1,numatomic(iatom))
              endif
! 2s+2p
              if(izcore(iatom) < 10) then
                iao= iao+1
                energy(iao)= row5c(2,numatomic(iatom))
                do i= 1,3
                  iao= iao+1
                  energy(iao)= row5c(3,numatomic(iatom))
                enddo
              endif
! 3s+3p
              if(izcore(iatom) < 18) then
                iao= iao+1
                energy(iao)= row5c(4,numatomic(iatom))
                do i= 1,3
                  iao= iao+1
                  energy(iao)= row5c(5,numatomic(iatom))
                enddo
              endif
! 4s+4p
              if(izcore(iatom) < 36) then
                iao= iao+1
                energy(iao)= row5c(6,numatomic(iatom))
                do i= 1,3
                  iao= iao+1
                  energy(iao)= row5c(7,numatomic(iatom))
                enddo
              endif
! 3d
              if(izcore(iatom) < 28) then
                do i= 1,5
                  iao= iao+1
                  energy(iao)= row5c(8,numatomic(iatom))
                enddo
              endif
! 4d
              if(numatomic(iatom) >= 49) then
                if(izcore(iatom) < 46) then
                  do i= 1,5
                    iao= iao+1
                    energy(iao)= row5v(3,numatomic(iatom))
                  enddo
                endif
              endif
            else
! 1s
              iao= iao+1
              energy(iao)= row5c(1,numatomic(iatom))
! 2s+2p
              iao= iao+1
              energy(iao)= row5c(2,numatomic(iatom))
              do i= 1,3
                iao= iao+1
                energy(iao)= row5c(3,numatomic(iatom))
              enddo
! 3s+3p
              iao= iao+1
              energy(iao)= row5c(4,numatomic(iatom))
              do i= 1,3
                iao= iao+1
                energy(iao)= row5c(5,numatomic(iatom))
              enddo
! 4s+4p
              iao= iao+1
              energy(iao)= row5c(6,numatomic(iatom))
              do i= 1,3
                iao= iao+1
                energy(iao)= row5c(7,numatomic(iatom))
              enddo
! 3d
              do i= 1,5
                iao= iao+1
                energy(iao)= row5c(8,numatomic(iatom))
              enddo
! 4d
              if(numatomic(iatom) >= 49) then
                do i= 1,5
                  iao= iao+1
                  energy(iao)= row5v(3,numatomic(iatom))
                enddo
              endif
            endif
! Au
          case(79)
            iao= iao+1
            energy(iao)= row6v(1,numatomic(iatom))
            do i= 1,3
              iao= iao+1
              energy(iao)= row6v(2,numatomic(iatom))
            enddo
          case default
            write(*,'(" Error! This program supports up to Xe in huckelip.")')
            call iabort
        end select
      enddo
      return
end


!-----------------------------------------------
  subroutine formhuckel(huckel,overlap,energy)
!-----------------------------------------------
!
! Form extended Huckel matrix
!
! In  : overlap (overlap integral of guess basis set)
!       energy (ionization potential)
! Out : huckel (extended Huckel Hamiltonian)
!
      use modguess, only : nao_g, nao_v
      implicit none
      integer :: i, j
      real(8),parameter :: factor=0.875D+00  !(=1.75/2.0)
      real(8),parameter :: fdown=0.05D+00
      real(8),intent(in) :: overlap(nao_g,nao_g), energy(nao_g)
      real(8),intent(out) :: huckel(nao_g,nao_g)
!
! Generate Huckel matrix from overlap integrals and ionization potentials
!
!$OMP parallel do schedule(static,1) private(j)
      do i= 1,nao_g
        do j= 1,i-1
          huckel(j,i)=factor*overlap(j,i)*(energy(i)+energy(j))
        enddo
        huckel(i,i)= energy(i)
      enddo
!$OMP end parallel do
!
! Scale down core-core and core-valence elements
!
!$OMP parallel do schedule(static,1) private(j)
      do i= nao_v+1,nao_g
        do j= 1,i-1
          huckel(j,i)= fdown*huckel(j,i)
        enddo
      enddo
!$OMP end parallel do
      return
end


!--------------------------------
  subroutine calcover1(overlap)
!--------------------------------
!
! Driver of overlap integral calculation
! (guess basis)x(guess basis)
!
! Out : overlap (overlap integral of guess basis set)
!
      use modguess, only : nshell_g, nao_g
      implicit none
      integer :: ish, jsh
      real(8),intent(out) :: overlap(nao_g*nao_g)
!
!$OMP parallel do private(jsh)
      do ish= nshell_g,1,-1
        do jsh= 1,ish
          call intover1(overlap,ish,jsh)
        enddo
      enddo
!$OMP end parallel do
      return
end


!-------------------------------------
  subroutine calcover2(overlap,work)
!-------------------------------------
!
! Driver of overlap integral calculation
! (input basis)x(guess basis)
!
! Out : overlap (overlap integral of guess and SCF basis sets)
!       work    (work array)
!
      use modparallel
      use modbasis, only : nshell, nao
      use modguess, only : nshell_g, nao_g
      implicit none
      integer :: ish, jsh
      real(8),parameter :: zero=0.0D+00
      real(8),intent(out) :: overlap(nao*nao_g), work(nao*nao_g)
!
      work(:)= zero
!$OMP parallel
      do ish= nshell_g-myrank,1,-nproc
!$OMP do
        do jsh= 1,nshell
          call intover2(work,ish,jsh)
        enddo
!$OMP enddo
      enddo
!$OMP end parallel
!
      call para_allreduce(work,overlap,nao*nao_g,MPI_SUM,MPI_COMM_WORLD)
      return
end


!-----------------------------------------
  subroutine intover1(overlap,ish,jsh)
!-----------------------------------------
!
! Overlap integral calculation
! (guess basis)x(guess basis)
!
! In  : ish, jsh (shell index)
! Out : overlap (overlap integral of guess basis set)
!
      use modparam, only : mxprsh
      use modthresh, only : threshex
      use modbasis, only : locatom, locprim, locbf, mprim, mtype, ex, coeff
      use modguess, only : locatom_g, locprim_g, locbf_g, mprim_g, mbf_g, mtype_g, &
&                       ex_g, coeff_g, nao_g, coord_g
      use modhermite, only : ix, iy, iz
      implicit none
      integer,intent(in) :: ish, jsh
      integer :: iatom, jatom, iloc, jloc, ilocbf, jlocbf, nprimi, nprimj, nangi, nangj 
      integer :: nbfi, nbfj, iprim, jprim, ncarti, ncartj, i, j, iang, jang
      integer :: isx, jsx, isy, jsy, isz, jsz, maxj
      integer :: ncart(0:5)=(/1,3,6,10,15,21/)
      real(8),parameter :: zero=0.0D+0, one=1.0D+0
      real(8),intent(out) :: overlap(nao_g,nao_g)
      real(8) :: xyzij(3), rij, rij2, fac, exi, exj, ci, cj, ex1, ex2, ex3, xyzpij(3,2)
      real(8) :: xyzint(3), sx(0:4,0:4), sy(0:4,0:4), sz(0:4,0:4), sint(28,28)
      logical :: iandj
!
! Set parameters
!
      iandj =(ish == jsh)
      iatom = locatom_g(ish)
      iloc  = locprim_g(ish)
      ilocbf= locbf_g(ish)
      nprimi= mprim_g(ish)
      nangi = mtype_g(ish)
      nbfi  = mbf_g(ish)
      ncarti= ncart(nangi)
      jatom = locatom_g(jsh)
      jloc  = locprim_g(jsh)
      jlocbf= locbf_g(jsh)
      nprimj= mprim_g(jsh)
      nangj = mtype_g(jsh)
      nbfj  = mbf_g(jsh)
      ncartj= ncart(nangj)
      do i= 1,3
        xyzij(i)= coord_g(i,iatom)-coord_g(i,jatom)
      enddo
      rij= xyzij(1)*xyzij(1)+xyzij(2)*xyzij(2)+xyzij(3)*xyzij(3)
      do i= 1,ncarti
        do j= 1,ncartj
          sint(j,i)=zero
        enddo
      enddo
!
! Calculate overlap integrals for each primitive
!
      do iprim= 1,nprimi
        exi= ex_g(iloc+iprim)
        ci = coeff_g(iloc+iprim)
        do jprim= 1,nprimj
          exj= ex_g(jloc+jprim)
          ex1= exi+exj
          ex2= one/ex1
          rij2=rij*exi*exj*ex2
          if(rij2 > threshex) cycle
          ex3= sqrt(ex2)
          fac= exp(-rij2)   !*ex2*ex3
          do i= 1,3
            xyzpij(i,1)=-exj*xyzij(i)*ex2
            xyzpij(i,2)= exi*xyzij(i)*ex2
          enddo
          cj = coeff_g(jloc+jprim)*fac
          do iang= 0,nangi
            do jang= 0,nangj
              call ghquad(xyzint,ex3,xyzpij,iang,jang)
              sx(jang,iang)= xyzint(1)*ex3
              sy(jang,iang)= xyzint(2)*ex3
              sz(jang,iang)= xyzint(3)*ex3
            enddo
          enddo
          do i= 1,ncarti
            isx= ix(i,nangi)
            isy= iy(i,nangi)
            isz= iz(i,nangi)
            do j= 1,ncartj
              jsx= ix(j,nangj)
              jsy= iy(j,nangj)
              jsz= iz(j,nangj)
              sint(j,i)= sint(j,i)+ci*cj*sx(jsx,isx)*sy(jsy,isy)*sz(jsz,isz)
            enddo
          enddo
        enddo
      enddo

      if((nbfi >= 5).or.(nbfj >= 5)) then
        call nrmlz1(sint,nbfi,nbfj,ncarti)
      endif

      maxj= nbfj
      do i= 1,nbfi
        if(iandj) maxj= i
        do j= 1,maxj
          overlap(jlocbf+j,ilocbf+i)= sint(j,i)
        enddo
      enddo
      return
end


!-----------------------------------------
  subroutine intover2(overlap,ish,jsh)
!-----------------------------------------
!
! Overlap integral calculation
! (input basis)x(guess basis)
!
! In  : ish, jsh (shell index)
! Out : overlap (overlap integral of guess and SCF basis sets)
!
      use modparam, only : mxprsh
      use modthresh, only : threshex
      use modmolecule, only : coord
      use modbasis, only : locatom, locprim, locbf, mprim, mbf, mtype, ex, coeff, nao
      use modguess, only : locatom_g, locprim_g, locbf_g, mprim_g, mbf_g, mtype_g, &
&                       ex_g, coeff_g, nao_g, coord_g
      use modhermite, only : ix, iy, iz
      implicit none
      integer,intent(in) :: ish, jsh
      integer :: iatom, jatom, iloc, jloc, ilocbf, jlocbf, nprimi, nprimj, nangi, nangj 
      integer :: nbfi, nbfj, iprim, jprim, ncarti, ncartj, i, j, iang, jang
      integer :: isx, jsx, isy, jsy, isz, jsz
      integer :: ncart(0:5)=(/1,3,6,10,15,21/)
      real(8),parameter :: zero=0.0D+0, one=1.0D+0
      real(8),intent(out) :: overlap(nao,nao_g)
      real(8) :: xyzij(3), rij, rij2, fac, exi, exj, ci, cj, ex1, ex2, ex3, xyzpij(3,2)
      real(8) :: xyzint(3), sx(0:5,0:5), sy(0:5,0:5), sz(0:5,0:5),sint(28,28)
!
! Set parameters
!
      iatom= locatom_g(ish)
      iloc= locprim_g(ish)
      ilocbf= locbf_g(ish)
      nprimi= mprim_g(ish)
      nangi= mtype_g(ish)
      nbfi= mbf_g(ish)
      ncarti= ncart(nangi)
      jatom= locatom(jsh)
      jloc= locprim(jsh)
      jlocbf= locbf(jsh)
      nprimj= mprim(jsh)
      nangj= mtype(jsh)
      nbfj= mbf(jsh)
      ncartj= ncart(nangj)
      do i= 1,3
!ishimura
!       xyzij(i)= coord(i,iatom)-coord_g(i,jatom)
        xyzij(i)= coord_g(i,iatom)-coord(i,jatom)
      enddo
      rij= xyzij(1)*xyzij(1)+xyzij(2)*xyzij(2)+xyzij(3)*xyzij(3)
      do i= 1,ncarti
        do j= 1,ncartj
          sint(j,i)=zero
        enddo
      enddo
!
! Calculate overlap integrals for each primitive
!
      do iprim= 1,nprimi
        exi= ex_g(iloc+iprim)
        ci = coeff_g(iloc+iprim)
        do jprim= 1,nprimj
          exj= ex(jloc+jprim)
          ex1= exi+exj
          ex2= one/ex1
          rij2=rij*exi*exj*ex2
          if(rij2 > threshex) cycle
          ex3= sqrt(ex2)
          fac= exp(-rij2)   !*ex2*ex3
          do i= 1,3
            xyzpij(i,1)=-exj*xyzij(i)*ex2
            xyzpij(i,2)= exi*xyzij(i)*ex2
          enddo
          cj = coeff(jloc+jprim)*fac
          do iang= 0,nangi
            do jang= 0,nangj
              call ghquad(xyzint,ex3,xyzpij,iang,jang)
              sx(jang,iang)= xyzint(1)*ex3
              sy(jang,iang)= xyzint(2)*ex3
              sz(jang,iang)= xyzint(3)*ex3
            enddo
          enddo
          do i= 1,ncarti
            isx= ix(i,nangi)
            isy= iy(i,nangi)
            isz= iz(i,nangi)
            do j= 1,ncartj
              jsx= ix(j,nangj)
              jsy= iy(j,nangj)
              jsz= iz(j,nangj)
              sint(j,i)= sint(j,i)+ci*cj*sx(jsx,isx)*sy(jsy,isy)*sz(jsz,isz)
            enddo
          enddo
        enddo
      enddo

      if((nbfi >= 5).or.(nbfj >= 5)) then
        call nrmlz1(sint,nbfi,nbfj,ncarti)
      endif

      do i= 1,nbfi
        do j= 1,nbfj
          overlap(jlocbf+j,ilocbf+i)= sint(j,i)
        enddo
      enddo
      return
end


!------------------------------------
  subroutine updatemo(cmo,overinv)
!------------------------------------
!
! Read and project MOs
!
! Inout : overinv (overlap integral inverse matrix)
!         cmo     (initial guess orbitals)
!
      use modguess, only : nao_g
      use modbasis, only : nao
      use modmolecule, only : neleca, nmo
      implicit none
      integer :: ndim
      real(8),intent(inout) :: cmo(nao*nao), overinv(nao*nao)
      real(8),allocatable :: overlap(:,:), work1(:), work2(:), eigen(:)
!
      ndim= min(nmo,neleca+5)
      ndim= nmo
!
! Set arrays
!
      call memset(2*nao*nao_g+ndim*ndim+nao*ndim+ndim)
      allocate(overlap(nao*nao_g,2),work1(ndim*ndim),work2(nao*ndim),eigen(ndim))
!
! Calculate overlap integrals between previous and present bases
!
      call calcover2(overlap(1,1),overlap(1,2))
!
! Project orbitals from previous basis to current basis
!
      call projectmo2(cmo,overinv,overlap,work1,work2,eigen,ndim)
!
! Unset arrays
!
      deallocate(overlap,work1,work2,eigen)
      call memunset(2*nao*nao_g+ndim*ndim+nao*ndim+ndim)
!
      return
end


!----------------------------------------------------------------------
  subroutine projectmo2(cmo,overinv,overlap,work1,work2,eigen,ndim)
!----------------------------------------------------------------------
!
! Project orbitals from previous basis to current basis
!    C1= S11^-1 * S12 * C2 [C2t * S12t * S11^-1 * S12 * C2]^-1/2
!
! In this routine, nao >= ndim is assumed.
!
! Inout : cmo     (previous and updated orbitals)
!         overinv (overlap integral inverse matrix of SCF basis set)
!         overlap (overlap integral of guess and SCF basis sets)
!
      use modguess, only : nao_g
      use modbasis, only : nao
      implicit none
      integer,intent(in) :: ndim
      integer :: i, j
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(inout) :: cmo(nao,nao), overinv(nao,nao), overlap(nao,nao_g)
      real(8),intent(out) :: work1(ndim,ndim), work2(ndim,nao), eigen(ndim)
      real(8) :: eigeninv
!
      if(nao < ndim) then
        write(*,'(" Error! Nao is less than ndim in projectmo2.")')
        call iabort
      endif
!
! Calculate S12*C2
!
      call dgemm('N','N',nao,ndim,nao_g,one,overlap,nao,cmo,nao_g,zero,work2,nao)
!
! Calculate S11^-1*S12*C2
!
      call dsymm('L','U',nao,ndim,one,overinv,nao,work2,nao,zero,overlap,nao)
!
! Calculate C2t*S12t*S11^-1*S12*C2
!
      call dgemm('T','N',ndim,ndim,nao,one,work2,nao,overlap,nao,zero,work1,ndim)
!
! Calculate (C2t*S12t*S11^-1*S12*C2)^-1/2
!
      call diag('V','U',ndim,work1,ndim,eigen)
!$OMP parallel do private(eigeninv)
      do i= 1,ndim
        eigeninv= one/sqrt(eigen(i))
        do j= 1,ndim
          work2(j,i)= work1(j,i)*eigeninv
        enddo
      enddo
!$OMP end parallel do
      call dgemm('N','T',ndim,ndim,ndim,one,work1,ndim,work2,ndim,zero,overinv,ndim)
!
! Calculate C1
!
      call dgemm('N','N',nao,ndim,ndim,one,overlap,nao,overinv,ndim,zero,cmo,nao)
      return
end


