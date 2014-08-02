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
      use modparallel
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
      call calcover2(overlap(1,1),overlap(1,2),nproc,myrank,MPI_COMM_WORLD)
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
      real(8) :: row2(3,3:10)
      real(8) :: row3(5,11:18)
      real(8) :: row4(8,19:36)
      real(8) :: row5(11,37:54)
      real(8) :: row6(15,55:86)
! The order is 1S,2S,2P
      data row2/-2.48D+00,-1.960D-01, 0.000D+00, -4.73D+00,-3.090D-01, 0.000D+00,&
&               -7.70D+00,-4.950D-01,-3.100D-01, -1.13D+01,-7.060D-01,-4.330D-01,&
&               -1.56D+01,-9.450D-01,-5.680D-01, -2.07D+01,-1.244D+00,-6.320D-01,&
&               -2.64D+01,-1.573D+00,-7.300D-01, -3.28D+01,-1.930D+00,-8.500D-01/
! The order is 1S,2S,2P,3S,3P
      data row3/-4.05D+01,-2.80D+00,-1.52D+00,-1.820D-01, 0.000D+00,&
&               -4.90D+01,-3.77D+00,-2.28D+00,-2.530D-01, 0.000D+00,&
&               -5.85D+01,-4.91D+00,-3.22D+00,-3.930D-01,-2.100D-01,&
&               -6.88D+01,-6.16D+00,-4.26D+00,-5.400D-01,-2.970D-01,&
&               -8.00D+01,-7.51D+00,-5.40D+00,-6.960D-01,-3.920D-01,&
&               -9.20D+01,-9.00D+00,-6.68D+00,-8.800D-01,-4.370D-01,&
&               -1.04D+02,-1.06D+01,-8.07D+00,-1.073D+00,-5.060D-01,&
&               -1.186D+02,-1.23D+01,-9.57D+00,-1.278D+00,-5.910D-01/
! The order is 1S,2S,2P,3S,3P,3D,4S,4P
      data row4/-1.335D+02,-1.450D+01,-1.150D+01,-1.750D+00,-9.500D-01,&
&                0.000D+00,-1.470D-01, 0.000D+00,&
&               -1.494D+02,-1.680D+01,-1.360D+01,-2.240D+00,-1.340D+00,&
&                0.000D+00,-1.960D-01, 0.000D+00,&
&               -1.659D+02,-1.908D+01,-1.567D+01,-2.570D+00,-1.575D+00,&
&               -1.193D+00,-4.240D-01,-2.080D-01,&
&               -1.833D+02,-2.142D+01,-1.779D+01,-2.874D+00,-1.795D+00,&
&               -3.430D-01,-2.100D-01,-1.000D-01,&
&               -2.013D+02,-2.370D+01,-1.980D+01,-2.990D+00,-1.840D+00,&
&               -4.410D-01,-2.200D-01,-1.000D-01,&
&               -2.204D+02,-2.620D+01,-2.210D+01,-3.290D+00,-2.050D+00,&
&               -3.210D-01,-2.140D-01,-1.000D-01,&
&               -2.404D+02,-2.890D+01,-2.460D+01,-3.620D+00,-2.300D+00,&
&               -3.730D-01,-2.220D-01,-1.000D-01,&
&               -2.612D+02,-3.170D+01,-2.720D+01,-3.960D+00,-2.550D+00,&
&               -3.830D-01,-2.270D-01,-1.000D-01,&
&               -2.829D+02,-3.460D+01,-2.990D+01,-4.300D+00,-2.800D+00,&
&               -4.060D-01,-2.300D-01,-1.000D-01,&
&               -3.054D+02,-3.770D+01,-3.270D+01,-4.650D+00,-3.060D+00,&
&               -4.340D-01,-2.330D-01,-1.000D-01,&
&               -3.288D+02,-4.080D+01,-3.560D+01,-5.010D+00,-3.320D+00,&
&               -4.570D-01,-2.360D-01,-1.000D-01,&
&               -3.533D+02,-4.440D+01,-3.890D+01,-5.630D+00,-3.840D+00,&
&               -4.910D-01,-2.380D-01,-1.000D-01,&
&               -3.788D+02,-4.820D+01,-4.250D+01,-6.400D+00,-4.480D+00,&
&               -7.830D-01,-2.930D-01,-1.000D-01,&
&               -4.052D+02,-5.210D+01,-4.620D+01,-7.190D+00,-5.170D+00,&
&               -1.635D+00,-5.530D-01,-2.870D-01,&
&               -4.326D+02,-5.630D+01,-5.020D+01,-8.030D+00,-5.880D+00,&
&               -2.113D+00,-6.860D-01,-3.690D-01,&
&               -4.609D+02,-6.070D+01,-5.430D+01,-8.930D+00,-6.660D+00,&
&               -2.650D+00,-8.380D-01,-4.030D-01,&
&               -4.901D+02,-6.520D+01,-5.860D+01,-9.870D+00,-7.480D+00,&
&               -3.220D+00,-9.930D-01,-4.570D-01,&
&               -5.202D+02,-6.990D+01,-6.300D+01,-1.080D+01,-8.330D+00,&
&               -3.825D+00,-1.153D+00,-5.240D-01/
! The order is 1S,2S,2P,3S,3P,3D,4S,4P,4D,5S,5P
      data row5/-5.510D+02,-7.500D+01,-6.790D+01,-1.210D+01,-9.500D+00,-4.700D+00,&
&               -1.520D+00,-8.100D-01, 0.000D+00,-1.380D-01, 0.000D+00,&
&               -5.837D+02,-8.040D+01,-7.300D+01,-1.350D+01,-1.070D+01,-5.700D+00,&
&               -1.900D+00,-1.100D+00, 0.000D+00,-1.780D-01, 0.000D+00,&
&               -6.168D+02,-8.581D+01,-7.816D+01,-1.476D+01,-1.185D+01,-6.599D+00,&
&               -2.168D+00,-1.300D+00,-2.499D-01,-1.958D-01,-1.000D-01,&
&               -6.507D+02,-9.138D+01,-8.348D+01,-1.606D+01,-1.302D+01,-7.515D+00,&
&               -2.418D+00,-1.487D+00,-3.365D-01,-2.070D-01,-1.000D-01,&
&               -6.854D+02,-9.700D+01,-8.880D+01,-1.720D+01,-1.400D+01,-8.300D+00,&
&               -2.530D+00,-1.550D+00,-2.990D-01,-2.140D-01,-1.000D-01,&
&               -7.212D+02,-1.029D+02,-9.450D+01,-1.860D+01,-1.530D+01,-9.300D+00,&
&               -2.760D+00,-1.720D+00,-3.570D-01,-2.220D-01,-1.000D-01,&
&               -7.579D+02,-1.089D+02,-1.002D+02,-2.000D+01,-1.660D+01,-1.030D+01,&
&               -3.000D+00,-1.910D+00,-3.770D-01,-2.220D-01,-1.000D-01,&
&               -7.955D+02,-1.152D+02,-1.062D+02,-2.140D+01,-1.780D+01,-1.130D+01,&
&               -3.260D+00,-2.100D+00,-4.120D-01,-2.220D-01,-1.000D-01,&
&               -8.340D+02,-1.216D+02,-1.124D+02,-2.290D+01,-1.920D+01,-1.240D+01,&
&               -3.500D+00,-2.290D+00,-4.510D-01,-2.200D-01,-1.000D-01,&
&               -8.735D+02,-1.281D+02,-1.187D+02,-2.440D+01,-2.050D+01,-1.350D+01,&
&               -3.750D+00,-2.480D+00,-4.880D-01,-2.200D-01,-1.000D-01,&
&               -9.138D+02,-1.349D+02,-1.252D+02,-2.590D+01,-2.190D+01,-1.470D+01,&
&               -4.000D+00,-2.680D+00,-5.370D-01,-2.200D-01,-1.000D-01,&
&               -9.554D+02,-1.421D+02,-1.321D+02,-2.770D+01,-2.360D+01,-1.610D+01,&
&               -4.450D+00,-3.050D+00,-7.630D-01,-2.650D-01,-1.000D-01,&
&               -9.978D+02,-1.494D+02,-1.392D+02,-2.960D+01,-2.540D+01,-1.760D+01,&
&               -4.980D+00,-3.510D+00,-1.063D+00,-3.720D-01,-1.970D-01,&
&               -1.041D+03,-1.570D+02,-1.465D+02,-3.160D+01,-2.720D+01,-1.920D+01,&
&               -5.510D+00,-3.970D+00,-1.369D+00,-4.760D-01,-2.650D-01,&
&               -1.086D+03,-1.648D+02,-1.540D+02,-3.360D+01,-1.920D+01,-2.080D+01,&
&               -6.060D+00,-4.450D+00,-1.688D+00,-5.820D-01,-3.350D-01,&
&               -1.131D+03,-1.728D+02,-1.617D+02,-3.580D+01,-3.110D+01,-2.250D+01,&
&               -6.650D+00,-4.950D+00,-2.038D+00,-7.010D-01,-3.600D-01,&
&               -1.177D+03,-1.809D+02,-1.697D+02,-3.790D+01,-3.310D+01,-2.430D+01,&
&               -7.240D+00,-5.470D+00,-2.401D+00,-8.210D-01,-4.030D-01,&
&               -1.224D+03,-1.893D+02,-1.778D+02,-4.020D+01,-3.520D+01,-2.610D+01,&
&               -7.860D+00,-6.010D+00,-2.778D+00,-9.440D-01,-4.570D-01/
! The order is 1S,2S,2P,3S,3P,3D,4S,4P,4D,4F,5S,5P,5D,6S,6P
      data row6/-2.5455D+03,-3.9630D+02,-3.7260D+02,-8.5390D+01,-7.5190D+01,-5.6450D+01,&
&               -1.7390D+01,-1.3540D+01,-6.7590D+00, 0.0000D+00,-2.4630D+00,-1.3670D+00,&
&                0.0000D+00,-2.4730D-01, 0.0000D+00,&
&               -2.6442D+03,-4.1430D+02,-3.9010D+02,-9.0560D+01,-8.0080D+01,-6.0800D+01,&
&               -1.9110D+01,-1.5100D+01,-8.0030D+00, 0.0000D+00,-3.0250D+00,-1.8080D+00,&
&                0.0000D+00,-3.1500D-01, 0.0000D+00,&
&               -2.7445D+03,-4.3260D+02,-4.0780D+02,-9.5690D+01,-8.4920D+01,-6.5110D+01,&
&               -2.0690D+01,-1.6520D+01,-9.1070D+00,-7.1800D-01,-3.4090D+00,-2.0990D+00,&
&                0.0000D+00,-3.4080D-01, 0.0000D+00,&
&               -2.8462D+03,-4.5050D+02,-4.2520D+02,-1.0020D+02,-8.9100D+01,-6.8770D+01,&
&               -2.1620D+01,-1.7310D+01,-9.6360D+00,-8.9300D-01,-3.5080D+00,-2.1560D+00,&
&                0.0000D+00,-3.4525D-01, 0.0000D+00,&
&               -2.9491D+03,-4.6820D+02,-4.4241D+02,-1.0400D+02,-9.2671D+01,-7.1830D+01,&
&               -2.1930D+01,-1.7492D+01,-9.6043D+00,-9.5104D-01,-3.3236D+00,-1.9759D+00,&
&                0.0000D+00,-3.2825D-01, 0.0000D+00,&
&               -3.0545D+03,-4.8680D+02,-4.6051D+02,-1.0859D+02,-9.6973D+01,-7.5608D+01,&
&               -2.2842D+01,-1.8266D+01,-1.0117D+01,-1.0277D+00,-3.4109D+00,-2.0226D+00,&
&                0.0000D+00,-3.3206D-01, 0.0000D+00,&
&               -3.1616D+03,-5.0580D+02,-4.7894D+02,-1.1324D+02,-1.0134D+02,-7.9449D+01,&
&               -2.3755D+01,-1.9038D+01,-1.0627D+01,-1.0952D+00,-3.4952D+00,-2.0667D+00,&
&                0.0000D+00,-3.3573D-01, 0.0000D+00,&
&               -3.2706D+03,-5.2510D+02,-4.9772D+02,-1.1796D+02,-1.0577D+02,-8.3355D+01,&
&               -2.4669D+01,-1.9811D+01,-1.1136D+01,-1.1553D+00,-3.5773D+00,-2.1091D+00,&
&                0.0000D+00,-3.3927D-01, 0.0000D+00,&
&               -3.3814D+03,-5.4470D+02,-5.1684D+02,-1.2275D+02,-1.1027D+02,-8.7329D+01,&
&               -2.5587D+01,-2.0586D+01,-1.1647D+01,-1.2092D+00,-3.6578D+00,-2.1500D+00,&
&                0.0000D+00,-3.4272D-01, 0.0000D+00,&
&               -3.4947D+03,-5.6537D+02,-5.3692D+02,-1.2831D+02,-1.1555D+02,-9.2059D+01,&
&               -2.7151D+01,-2.1995D+01,-1.2750D+01,-1.2170D+00,-4.0451D+00,-2.4489D+00,&
&                0.0000D+00,-3.6994D-01, 0.0000D+00,&
&               -3.6093D+03,-5.8572D+02,-5.5674D+02,-1.3326D+02,-1.2020D+02,-9.6182D+01,&
&               -2.8088D+01,-2.2787D+01,-1.9000D+00,-1.2010D+00,-4.1303D+00,-2.4936D+00,&
&                0.0000D+00,-3.7392D-01, 0.0000D+00,&
&               -3.7251D+03,-6.0580D+02,-5.7627D+02,-1.3756D+02,-1.2422D+02,-9.9678D+01,&
&               -2.8371D+01,-2.2939D+01,-1.3196D+01,-1.3411D+00,-3.8927D+00,-2.2663D+00,&
&                0.0000D+00,-3.5258D-01, 0.0000D+00,&
&               -3.8433D+03,-6.2681D+02,-5.9678D+02,-1.4265D+02,-1.2901D+02,-1.0394D+02,&
&               -2.9313D+01,-2.3736D+01,-1.3720D+01,-1.3769D+00,-3.9695D+00,-2.3035D+00,&
&                0.0000D+00,-3.5574D-01, 0.0000D+00,&
&               -3.9634D+03,-6.4820D+02,-6.1764D+02,-1.4781D+02,-1.3389D+02,-1.0830D+02,&
&               -3.0300D+01,-2.4539D+01,-1.4249D+01,-1.4093D+00,-4.0459D+00,-2.3401D+00,&
&                0.0000D+00,-3.5885D-01, 0.0000D+00,&
&               -4.0854D+03,-6.6990D+02,-6.3884D+02,-1.5306D+02,-1.3884D+02,-1.1260D+02,&
&               -3.1220D+01,-2.5349D+01,-1.4783D+01,-1.4385D+00,-4.1219D+00,-2.3761D+00,&
&                0.0000D+00,-3.6191D-01, 0.0000D+00,&
&               -4.2092D+03,-6.9200D+02,-6.6039D+02,-1.5838D+02,-1.4386D+02,-1.1720D+02,&
&               -3.2190D+01,-2.6167D+01,-1.5323D+01,-1.4648D+00,-4.1976D+00,-2.4117D+00,&
&                0.0000D+00,-3.6492D-01, 0.0000D+00,&
&               -4.3354D+03,-7.1510D+02,-6.8298D+02,-1.6454D+02,-1.4972D+02,-1.2250D+02,&
&               -3.3880D+01,-2.7690D+01,-1.6530D+01,-2.1540D+00,-4.6340D+00,-2.7520D+00,&
&               -4.8670D-01,-3.9770D-01, 0.0000D+00,&
&               -4.4636D+03,-7.3870D+02,-7.0600D+02,-1.7080D+02,-1.5570D+02,-1.2790D+02,&
&               -3.5610D+01,-2.9260D+01,-1.7770D+01,-2.8720D+00,-5.0500D+00,-3.0740D+00,&
&               -5.9830D-01,-4.2080D-01, 0.0000D+00,&
&               -4.5937D+03,-7.6260D+02,-7.2940D+02,-1.7730D+02,-1.6190D+02,-1.3350D+02,&
&               -3.7390D+01,-3.0870D+01,-1.9060D+01,-3.6300D+00,-5.4590D+00,-3.3930D+00,&
&               -7.0330D-01,-4.3950D-01, 0.0000D+00,&
&               -4.7256D+03,-7.8690D+02,-7.5310D+02,-1.8390D+02,-1.6810D+02,-1.3920D+02,&
&               -3.9220D+01,-3.2530D+01,-2.0390D+01,-4.4290D+00,-5.8680D+00,-3.7130D+00,&
&               -8.0580D-01,-4.5540D-01, 0.0000D+00,&
&               -4.8595D+03,-8.1170D+02,-7.7730D+02,-1.9060D+02,-1.7460D+02,-1.4510D+02,&
&               -4.1100D+01,-3.4240D+01,-2.1770D+01,-5.2670D+00,-6.2770D+00,-4.0340D+00,&
&               -9.0760D-01,-4.6930D-01, 0.0000D+00,&
&               -4.9953D+03,-8.3680D+02,-8.0190D+02,-1.9740D+02,-1.8110D+02,-1.5100D+02,&
&               -4.3030D+01,-3.5990D+01,-2.3180D+01,-6.1430D+00,-6.6880D+00,-4.3590D+00,&
&               -1.0096D+00,-4.8183D-01, 0.0000D+00,&
&               -5.1329D+03,-8.6240D+02,-8.2690D+02,-2.0450D+02,-1.8780D+02,-1.5720D+02,&
&               -4.4990D+01,-3.7780D+01,-2.4630D+01,-7.0580D+00,-7.1020D+00,-4.6880D+00,&
&               -1.1124D+00,-4.9317D-01, 0.0000D+00,&
&               -5.2721D+03,-8.8800D+02,-8.5200D+02,-2.1130D+02,-1.9430D+02,-1.6310D+02,&
&               -4.6680D+01,-3.9290D+01,-2.5800D+01,-7.6870D+00,-7.2130D+00,-4.7470D+00,&
&               -9.5290D-01,-4.3590D-01, 0.0000D+00,&
&               -5.4136D+03,-9.1440D+02,-8.7780D+02,-2.1850D+02,-2.0120D+02,-1.6940D+02,&
&               -4.8710D+01,-4.1140D+01,-2.7310D+01,-8.6570D+00,-7.6180D+00,-5.0690D+00,&
&               -1.0420D+00,-4.4155D-01, 0.0000D+00,&
&               -5.5573D+03,-9.4150D+02,-9.0440D+02,-2.2630D+02,-2.0870D+02,-1.7630D+02,&
&               -5.1150D+01,-4.3400D+01,-2.9220D+01,-1.0020D+01,-8.3640D+00,-5.7020D+00,&
&               -1.4284D+00,-5.2209D-01, 0.0000D+00,&
&               -5.7030D+03,-9.6910D+02,-9.3150D+02,-2.3430D+02,-2.1640D+02,-1.8340D+02,&
&               -5.3770D+01,-4.5840D+01,-3.1310D+01,-1.1570D+01,-9.2370D+00,-6.4630D+00,&
&               -1.9370D+00,-7.2220D-01,-3.8480D-01,&
&               -5.8507D+03,-9.9720D+02,-9.5900D+02,-2.4250D+02,-2.2430D+02,-1.9070D+02,&
&               -5.6450D+01,-4.8330D+01,-3.3450D+01,-1.3170D+01,-1.0120D+01,-7.2290D+00,&
&               -2.4490D+00,-9.1770D-01,-4.7970D-01,&
&               -6.0002D+03,-1.0257D+03,-9.8690D+02,-2.5080D+02,-2.3230D+02,-1.9810D+02,&
&               -5.9200D+01,-5.0900D+01,-3.5660D+01,-1.4840D+01,-1.1020D+01,-8.0100D+00,&
&               -2.9750D+00,-1.1160D+00,-5.7240D-01,&
&               -6.1517D+03,-1.0546D+03,-1.0153D+03,-2.5930D+02,-2.4050D+02,-2.0570D+02,&
&               -6.2010D+01,-5.3530D+01,-3.7930D+01,-1.6570D+01,-1.1930D+01,-8.8060D+00,&
&               -3.5170D+00,-1.3200D+00,-6.6540D-01,&
&               -6.3052D+03,-1.0840D+03,-1.0441D+03,-2.6800D+02,-2.4880D+02,-2.1350D+02,&
&               -6.4890D+01,-5.6220D+01,-4.0260D+01,-1.8360D+01,-1.2860D+01,-9.6200D+00,&
&               -4.0760D+00,-1.5310D+00,-7.5970D-01,&
&               -6.4605D+03,-1.1138D+03,-1.0734D+03,-2.7680D+02,-2.5730D+02,-2.2140D+02,&
&               -6.7840D+01,-5.8980D+01,-4.2660D+01,-2.0220D+01,-1.3810D+01,-1.0450D+01,&
&               -4.6530D+00,-1.7480D+00,-8.5600D-01/
! The order is 1S,2S,2P,3S,3P,3D,4S,4P,4D,4F,5S,5P,5D,5F,6S,6P,6D,7S
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
            energy(iao)= row2(2,numatomic(iatom))
            do i= 1,3
              iao= iao+1
              energy(iao)= row2(3,numatomic(iatom))
            enddo
! Na - Ar
          case(11:18)
            iao= iao+1
            energy(iao)= row3(4,numatomic(iatom))
            do i= 1,3
              iao= iao+1
              energy(iao)= row3(5,numatomic(iatom))
            enddo
! K  - Ca, Ga - Kr
          case(19:20,31:36)
            iao= iao+1
            energy(iao)= row4(7,numatomic(iatom))
            do i= 1,3
              iao= iao+1
              energy(iao)= row4(8,numatomic(iatom))
            enddo
! Sc - Zn
          case(21:30)
            iao= iao+1
            energy(iao)= row4(7,numatomic(iatom))
            do i= 1,3
              iao= iao+1
              energy(iao)= row4(8,numatomic(iatom))
            enddo
            do i= 1,5
              iao= iao+1
              energy(iao)= row4(6,numatomic(iatom))
            enddo
! Rb - Sr, In - Xe
          case(37:38,49:54)
            iao= iao+1
            energy(iao)= row5(10,numatomic(iatom))
            do i= 1,3
              iao= iao+1
              energy(iao)= row5(11,numatomic(iatom))
            enddo
! Y  - Cd
          case(39:48)
            iao= iao+1
            energy(iao)= row5(10,numatomic(iatom))
            do i= 1,3
              iao= iao+1
              energy(iao)= row5(11,numatomic(iatom))
            enddo
            do i= 1,5
              iao= iao+1
              energy(iao)= row5(9,numatomic(iatom))
            enddo
! Cs - Ba
          case(55:56)
            iao= iao+1
            energy(iao)= row6(14,numatomic(iatom))
! La - Yb
          case(57:70)
            iao= iao+1
            energy(iao)= row6(14,numatomic(iatom))
            do i= 1,7
              iao= iao+1
              energy(iao)= row6(9,numatomic(iatom))
            enddo
! Lu - Hg
          case(71:80)
            iao= iao+1
            energy(iao)= row6(14,numatomic(iatom))
            do i= 1,5
              iao= iao+1
              energy(iao)= row6(13,numatomic(iatom))
            enddo
! Tl - Rn
          case(81:86)
            iao= iao+1
            energy(iao)= row6(14,numatomic(iatom))
            do i= 1,3
              iao= iao+1
              energy(iao)= row6(15,numatomic(iatom))
            enddo
          case default
            write(*,'(" Error! This program supports up to Rn in huckelip.")')
            call iabort
        end select
      enddo
!
! Set core ionization potentials
!
      do iatom= 1,natom
        select case(numatomic(iatom))
! H  - He
          case(1:2)
! Li - Ne
          case(3:10)
!  1s
            if(.not.flagecp.or.(izcore(iatom) < 2)) then
              iao= iao+1
              energy(iao)= row2(1,numatomic(iatom))
            endif
! Na - Ar
          case(11:18)
!  1s
            if(.not.flagecp.or.(izcore(iatom) < 2)) then
              iao= iao+1
              energy(iao)=row3(1,numatomic(iatom))
            endif
!  2s+2p
            if(.not.flagecp.or.(izcore(iatom) < 10)) then
              iao= iao+1
              energy(iao)=row3(2,numatomic(iatom))
              do i= 1,3
                iao= iao+1
                energy(iao)=row3(3,numatomic(iatom))
              enddo
            endif
! K  - Kr
          case(19:36)
!  1s
            if(.not.flagecp.or.(izcore(iatom) < 2)) then
              iao= iao+1
              energy(iao)= row4(1,numatomic(iatom))
            endif
!  2s+2p
            if(.not.flagecp.or.(izcore(iatom) < 10)) then
              iao= iao+1
              energy(iao)= row4(2,numatomic(iatom))
              do i= 1,3
                iao= iao+1
                energy(iao)= row4(3,numatomic(iatom))
              enddo
            endif
!  3s+3p
            if(.not.flagecp.or.(izcore(iatom) < 18)) then
              iao= iao+1
              energy(iao)= row4(4,numatomic(iatom))
              do i= 1,3
                iao= iao+1
                energy(iao)= row4(5,numatomic(iatom))
              enddo
            endif
!  3d
            if(numatomic(iatom) >= 31) then
              if(.not.flagecp.or.(izcore(iatom) < 28)) then
                do i= 1,5
                  iao= iao+1
                  energy(iao)= row4(6,numatomic(iatom))
                enddo
              endif
            endif
! Rb - Xe
          case(37:54)
!  1s
            if(.not.flagecp.or.(izcore(iatom) < 2)) then
              iao= iao+1
              energy(iao)= row5(1,numatomic(iatom))
            endif
!  2s+2p
            if(.not.flagecp.or.(izcore(iatom) < 10)) then
              iao= iao+1
              energy(iao)= row5(2,numatomic(iatom))
              do i= 1,3
                iao= iao+1
                energy(iao)= row5(3,numatomic(iatom))
              enddo
            endif
!  3s+3p
            if(.not.flagecp.or.(izcore(iatom) < 18)) then
              iao= iao+1
              energy(iao)= row5(4,numatomic(iatom))
              do i= 1,3
                iao= iao+1
                energy(iao)= row5(5,numatomic(iatom))
              enddo
            endif
!  3d
            if(.not.flagecp.or.(izcore(iatom) < 28)) then
              do i= 1,5
                iao= iao+1
                energy(iao)= row5(6,numatomic(iatom))
              enddo
            endif
!  4s+4p
            if(.not.flagecp.or.(izcore(iatom) < 36)) then
              iao= iao+1
              energy(iao)= row5(7,numatomic(iatom))
              do i= 1,3
                iao= iao+1
                energy(iao)= row5(8,numatomic(iatom))
              enddo
            endif
!  4d
            if(numatomic(iatom) >= 49) then
              if(.not.flagecp.or.(izcore(iatom) < 46)) then
                do i= 1,5
                  iao= iao+1
                  energy(iao)= row5(9,numatomic(iatom))
                enddo
              endif
            endif
! Cs - Rn
          case(55:86)
!  1s
            if(.not.flagecp.or.(izcore(iatom) < 2)) then
              iao= iao+1
              energy(iao)= row6(1,numatomic(iatom))
            endif
!  2s+2p
            if(.not.flagecp.or.(izcore(iatom) < 10)) then
              iao= iao+1
              energy(iao)= row6(2,numatomic(iatom))
              do i= 1,3
                iao= iao+1
                energy(iao)= row6(3,numatomic(iatom))
              enddo
            endif
!  3s+3p
            if(.not.flagecp.or.(izcore(iatom) < 18)) then
              iao= iao+1
              energy(iao)= row6(4,numatomic(iatom))
              do i= 1,3
                iao= iao+1
                energy(iao)= row6(5,numatomic(iatom))
              enddo
            endif
!  3d
            if(.not.flagecp.or.(izcore(iatom) < 28)) then
              do i= 1,5
                iao= iao+1
                energy(iao)= row6(6,numatomic(iatom))
              enddo
            endif
!  4s+4p
            if(.not.flagecp.or.(izcore(iatom) < 36)) then
              iao= iao+1
              energy(iao)= row6(7,numatomic(iatom))
              do i= 1,3
                iao= iao+1
                energy(iao)= row6(8,numatomic(iatom))
              enddo
            endif
!  4d
            if(.not.flagecp.or.(izcore(iatom) < 46)) then
              do i= 1,5
                iao= iao+1
                energy(iao)= row6(9,numatomic(iatom))
              enddo
            endif
!  4f
            if(numatomic(iatom) >= 71) then
              if(.not.flagecp.or.(izcore(iatom) < 60)) then
                do i= 1,7
                  iao= iao+1
                  energy(iao)= row6(10,numatomic(iatom))
                enddo
              endif
            endif
!  5s+5p
            if(.not.flagecp.or.((izcore(iatom) /= 54).and.(izcore(iatom) < 68))) then
              iao= iao+1
              energy(iao)= row6(11,numatomic(iatom))
              do i= 1,3
                iao= iao+1
                energy(iao)= row6(12,numatomic(iatom))
              enddo
            endif
!  5d
            if(numatomic(iatom) >= 81) then
              if(.not.flagecp.or.(izcore(iatom) < 78)) then
                do i= 1,5
                  iao= iao+1
                  energy(iao)= row6(13,numatomic(iatom))
                enddo
              endif
            endif
          case default
            write(*,'(" Error! This program supports up to Rn in huckelip.")')
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


!-----------------------------------------------------------
  subroutine calcover2(overlap,work,nproc,myrank,mpi_comm)
!-----------------------------------------------------------
!
! Driver of overlap integral calculation
! (input basis)x(guess basis)
!
! Out : overlap (overlap integral of guess and SCF basis sets)
!       work    (work array)
!
      use modbasis, only : nshell, nao
      use modguess, only : nshell_g, nao_g
      implicit none
      integer,intent(in) :: nproc, myrank
      integer(4),intent(in) :: mpi_comm
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
      call para_allreducer(work,overlap,nao*nao_g,mpi_comm)
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
      use modguess, only : locatom_g, locprim_g, locbf_g, mprim_g, mbf_g, mtype_g, &
&                          ex_g, coeff_g, nao_g, coord_g
      implicit none
      integer,intent(in) :: ish, jsh
      integer :: iatom, jatom, iloc, jloc, ilocbf, jlocbf, iprim, jprim
      integer :: nbfij(2), nprimij(2), nangij(2), ii, jj, maxj
      real(8),intent(out) :: overlap(nao_g,nao_g)
      real(8) :: sint(28,28), exij(mxprsh,2), coij(mxprsh,2), coordij(3,2)
      logical :: iandj
!
! Set parameters
!
      iandj =(ish == jsh)
      nangij(1)= mtype_g(ish)
      nangij(2)= mtype_g(jsh)
      nprimij(1)= mprim_g(ish)
      nprimij(2)= mprim_g(jsh)
      nbfij(1)  = mbf_g(ish)
      nbfij(2)  = mbf_g(jsh)
      iatom = locatom_g(ish)
      iloc  = locprim_g(ish)
      ilocbf= locbf_g(ish)
      jatom = locatom_g(jsh)
      jloc  = locprim_g(jsh)
      jlocbf= locbf_g(jsh)
      do ii= 1,3
        coordij(ii,1)= coord_g(ii,iatom)
        coordij(ii,2)= coord_g(ii,jatom)
      enddo
      do iprim= 1,nprimij(1)
        exij(iprim,1)= ex_g(iloc+iprim)
        coij(iprim,1)= coeff_g(iloc+iprim)
      enddo
      do jprim= 1,nprimij(2)
        exij(jprim,2)= ex_g(jloc+jprim)
        coij(jprim,2)= coeff_g(jloc+jprim)
      enddo
!
! Calculate overlap integrals
!
      call ints(sint,exij,coij,coordij,nprimij,nangij,nbfij,mxprsh,threshex)
!
      maxj= nbfij(2)
      do ii= 1,nbfij(1)
        if(iandj) maxj= ii
        do jj= 1,maxj
          overlap(jlocbf+jj,ilocbf+ii)= sint(jj,ii)
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
&                          ex_g, coeff_g, nao_g, coord_g
      implicit none
      integer,intent(in) :: ish, jsh
      integer :: iatom, jatom, iloc, jloc, ilocbf, jlocbf, iprim, jprim
      integer :: nbfij(2), nprimij(2), nangij(2), ii, jj
      real(8),parameter :: zero=0.0D+0, one=1.0D+0
      real(8),intent(out) :: overlap(nao,nao_g)
      real(8) :: sint(28,28), exij(mxprsh,2), coij(mxprsh,2), coordij(3,2)
!
! Set parameters
!
      nangij(1)= mtype_g(ish)
      nangij(2)= mtype(jsh)
      nprimij(1)= mprim_g(ish)
      nprimij(2)= mprim(jsh)
      nbfij(1)  = mbf_g(ish)
      nbfij(2)  = mbf(jsh)
      iatom = locatom_g(ish)
      iloc  = locprim_g(ish)
      ilocbf= locbf_g(ish)
      jatom = locatom(jsh)
      jloc  = locprim(jsh)
      jlocbf= locbf(jsh)
      do ii= 1,3
        coordij(ii,1)= coord_g(ii,iatom)
        coordij(ii,2)= coord(ii,jatom)
      enddo
      do iprim= 1,nprimij(1)
        exij(iprim,1)= ex_g(iloc+iprim)
        coij(iprim,1)= coeff_g(iloc+iprim)
      enddo
      do jprim= 1,nprimij(2)
        exij(jprim,2)= ex(jloc+jprim)
        coij(jprim,2)= coeff(jloc+jprim)
      enddo
!
! Calculate overlap integrals
!
      call ints(sint,exij,coij,coordij,nprimij,nangij,nbfij,mxprsh,threshex)
!
      do ii= 1,nbfij(1)
        do jj= 1,nbfij(2)
          overlap(jlocbf+jj,ilocbf+ii)= sint(jj,ii)
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
      use modparallel
      use modguess, only : nao_g, nmo_g
      use modbasis, only : nao
      use modmolecule, only : neleca, nmo
      implicit none
      integer :: ndim
      real(8),intent(inout) :: cmo(nao*nao), overinv(nao*nao)
      real(8),allocatable :: overlap(:,:), work1(:), work2(:), eigen(:)
!
!     ndim= min(nmo,neleca+5)
      ndim= nmo_g
!
! Set arrays
!
      call memset(2*nao*nao_g+ndim*ndim+nao*ndim+ndim)
      allocate(overlap(nao*nao_g,2),work1(ndim*ndim),work2(nao*ndim),eigen(ndim))
!
! Calculate overlap integrals between previous and present bases
!
      call calcover2(overlap(1,1),overlap(1,2),nproc,myrank,MPI_COMM_WORLD)
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


