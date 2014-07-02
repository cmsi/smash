!--------------------------------------------
  subroutine rysquad(tval,trys,wrys,nroots)
!--------------------------------------------
!
! Calculate abscissas and weights of Rys quadrature
!
! In  : tval, nroots
! Out : trys, wrys
!
      use rysgrid1, only : t1, w1
      use rysgrid2, only : t2, w2
      use rysgrid3, only : t3, w3
      use rysgrid4, only : t4, w4
      use rysgrid5, only : t5, w5
      use rysgrid6, only : t6, w6
      use rysgrid7, only : t7, w7
      use rysgrid8, only : t8, w8
      use rysgrid9, only : t9, w9
      implicit none
      integer,intent(in) :: nroots
      integer :: j, igrid
      real(8),parameter :: half=0.5D+00, one=1.0D+00, two=2.0D+00, three=3.0D+00
      real(8),parameter :: third=0.3333333333333333D+00, fifth=0.2D+00, five=5.0D+00
      real(8),intent(in) :: tval
      real(8),intent(out) :: trys(13), wrys(13)
      real(8) :: tgh1, tgh2(2), tgh3(3), tgh4(4), tgh5(5), tgh6(6), tgh7(7), tgh8(8), tgh9(9)
      real(8) :: wgh1, wgh2(2), wgh3(3), wgh4(4), wgh5(5), wgh6(6), wgh7(7), wgh8(8), wgh9(9)
      real(8) :: tinv, tinvrt, tval2
      data tgh1/ 5.000000000000000D-01/
      data tgh2/ 2.752551286084110D-01, 2.724744871391589D+00/
      data tgh3/ 1.901635091934881D-01, 1.784492748543252D+00, 5.525343742263260D+00/
      data tgh4/ 1.453035215033171D-01, 1.339097288126361D+00, 3.926963501358287D+00, &
&                8.588635689012034D+00/
      data tgh5/ 1.180718948997174D+01, 6.414729733662030D+00, 3.085937443717550D+00, &
&                1.074562012436904D+00, 1.175813202117781D-01/
      data tgh6/ 9.874701406848118D-02, 8.983028345696177D-01, 2.552589802668171D+00, &
&                5.196152530054466D+00, 9.124248037531179D+00, 1.512995978110809D+01/
      data tgh7/ 8.511544299759403D-02, 7.721379200427770D-01, 2.180591888450459D+00, &
&                4.389792886731014D+00, 7.554091326101784D+00, 1.198999303982388D+01, &
&                1.852827749585249D+01/
      data tgh8/ 7.479188259681827D-02, 6.772490876492892D-01, 1.905113635031428D+00, &
&                3.809476361484907D+00, 6.483145428627170D+00, 1.009332367522134D+01, &
&                1.497262708842639D+01, 2.198427284096265D+01/
      data tgh9/ 6.670223095819440D-02, 6.032363570817487D-01, 1.692395079793179D+00, &
&                3.369176270243269D+00, 5.694423342957755D+00, 8.769756730268602D+00, &
&                1.277182535486919D+01, 1.804650546772898D+01, 2.548597916609908D+01/
      data wgh1/ 8.862269254527580D-01/
      data wgh2/ 8.049140900055128D-01, 8.131283544724518D-02/
      data wgh3/ 7.246295952243925D-01, 1.570673203228566D-01, 4.530009905508846D-03/
      data wgh4/ 6.611470125582413D-01, 2.078023258148919D-01, 1.707798300741348D-02, &
&                1.996040722113676D-04/
      data wgh5/ 7.640432855232621D-06, 1.343645746781233D-03, 3.387439445548106D-02, &
&                2.401386110823147D-01, 6.108626337353258D-01/
      data wgh6/ 5.701352362624796D-01, 2.604923102641611D-01, 5.160798561588393D-02, &
&                3.905390584629062D-03, 8.573687043587859D-05, 2.658551684356302D-07/
      data wgh7/ 5.364059097120901D-01, 2.731056090642466D-01, 6.850553422346521D-02, &
&                7.850054726457944D-03, 3.550926135519236D-04, 4.716484355018917D-06, &
&                8.628591168125158D-09/
      data wgh8/ 5.079294790166137D-01, 2.806474585285337D-01, 8.381004139898583D-02, &
&                1.288031153550997D-02, 9.322840086241805D-04, 2.711860092537882D-05, &
&                2.320980844865211D-07, 2.654807474011182D-10/
      data wgh9/ 4.834956947254556D-01, 2.848072856699796D-01, 9.730174764131543D-02, &
&                1.864004238754465D-02, 1.888522630268418D-03, 9.181126867929404D-05, &
&                1.810654481093430D-06, 1.046720579579208D-08, 7.828199772115891D-12/
!
      select case (nroots)
        case(1)
          if(tval < 35.0D+00)then
            igrid=nint(tval*half)
            tval2=tval-igrid*two
            trys(1)=((((((((((( &
&               t1( 1,igrid) *tval2+t1( 2,igrid))*tval2+t1( 3,igrid))*tval2+ &
&               t1( 4,igrid))*tval2+t1( 5,igrid))*tval2+t1( 6,igrid))*tval2+ &
&               t1( 7,igrid))*tval2+t1( 8,igrid))*tval2+t1( 9,igrid))*tval2+ &
&               t1(10,igrid))*tval2+t1(11,igrid))*tval2+t1(12,igrid))*tval2+ &
&               t1(13,igrid)
            wrys(1)=((((((((((( &
&               w1( 1,igrid) *tval2+w1( 2,igrid))*tval2+w1( 3,igrid))*tval2+ &
&               w1( 4,igrid))*tval2+w1( 5,igrid))*tval2+w1( 6,igrid))*tval2+ &
&               w1( 7,igrid))*tval2+w1( 8,igrid))*tval2+w1( 9,igrid))*tval2+ &
&               w1(10,igrid))*tval2+w1(11,igrid))*tval2+w1(12,igrid))*tval2+ &
&               w1(13,igrid)
          else
            tinvrt=one/sqrt(tval)
            trys(1)= tgh1*tinvrt*tinvrt
            wrys(1)= wgh1*tinvrt
          endif
        case (2)
          if(tval < 43.0D+00)then
            igrid=nint(tval*half)
            tval2=tval-igrid*two
            do j= 1,2
              trys(j)=((((((((((( &
&                 t2( 1,j,igrid) *tval2+t2( 2,j,igrid))*tval2+t2( 3,j,igrid))*tval2+ &
&                 t2( 4,j,igrid))*tval2+t2( 5,j,igrid))*tval2+t2( 6,j,igrid))*tval2+ &
&                 t2( 7,j,igrid))*tval2+t2( 8,j,igrid))*tval2+t2( 9,j,igrid))*tval2+ &
&                 t2(10,j,igrid))*tval2+t2(11,j,igrid))*tval2+t2(12,j,igrid))*tval2+ &
&                 t2(13,j,igrid)
            enddo
            do j= 1,2
              wrys(j)=((((((((((( &
&                 w2( 1,j,igrid) *tval2+w2( 2,j,igrid))*tval2+w2( 3,j,igrid))*tval2+ &
&                 w2( 4,j,igrid))*tval2+w2( 5,j,igrid))*tval2+w2( 6,j,igrid))*tval2+ &
&                 w2( 7,j,igrid))*tval2+w2( 8,j,igrid))*tval2+w2( 9,j,igrid))*tval2+ &
&                 w2(10,j,igrid))*tval2+w2(11,j,igrid))*tval2+w2(12,j,igrid))*tval2+ &
&                 w2(13,j,igrid)
            enddo
          else
            tinvrt=one/sqrt(tval)
            tinv=tinvrt*tinvrt
            do j= 1,2
              trys(j)= tgh2(j)*tinv
              wrys(j)= wgh2(j)*tinvrt
            enddo
          endif
        case (3)
          if(tval < 49.0D+00)then
            igrid=nint(tval*half)
            tval2=tval-igrid*two
            do j= 1,3
              trys(j)=((((((((((( &
&                 t3( 1,j,igrid) *tval2+t3( 2,j,igrid))*tval2+t3( 3,j,igrid))*tval2+ &
&                 t3( 4,j,igrid))*tval2+t3( 5,j,igrid))*tval2+t3( 6,j,igrid))*tval2+ &
&                 t3( 7,j,igrid))*tval2+t3( 8,j,igrid))*tval2+t3( 9,j,igrid))*tval2+ &
&                 t3(10,j,igrid))*tval2+t3(11,j,igrid))*tval2+t3(12,j,igrid))*tval2+ &
&                 t3(13,j,igrid)
            enddo
            do j= 1,3
              wrys(j)=((((((((((( &
&                 w3( 1,j,igrid) *tval2+w3( 2,j,igrid))*tval2+w3( 3,j,igrid))*tval2+ &
&                 w3( 4,j,igrid))*tval2+w3( 5,j,igrid))*tval2+w3( 6,j,igrid))*tval2+ &
&                 w3( 7,j,igrid))*tval2+w3( 8,j,igrid))*tval2+w3( 9,j,igrid))*tval2+ &
&                 w3(10,j,igrid))*tval2+w3(11,j,igrid))*tval2+w3(12,j,igrid))*tval2+ &
&                 w3(13,j,igrid)
            enddo
          else
            tinvrt=one/sqrt(tval)
            tinv=tinvrt*tinvrt
            do j= 1,3
              trys(j)= tgh3(j)*tinv
              wrys(j)= wgh3(j)*tinvrt
            enddo
          endif
        case (4)
          if(tval < 55.0D+00)then
            igrid=nint(tval*half)
            tval2=tval-igrid*two
            do j= 1,4
              trys(j)=((((((((((( &
&                 t4( 1,j,igrid) *tval2+t4( 2,j,igrid))*tval2+t4( 3,j,igrid))*tval2+ &
&                 t4( 4,j,igrid))*tval2+t4( 5,j,igrid))*tval2+t4( 6,j,igrid))*tval2+ &
&                 t4( 7,j,igrid))*tval2+t4( 8,j,igrid))*tval2+t4( 9,j,igrid))*tval2+ &
&                 t4(10,j,igrid))*tval2+t4(11,j,igrid))*tval2+t4(12,j,igrid))*tval2+ &
&                 t4(13,j,igrid)
            enddo
            do j= 1,4
              wrys(j)=((((((((((( &
&                 w4( 1,j,igrid) *tval2+w4( 2,j,igrid))*tval2+w4( 3,j,igrid))*tval2+ &
&                 w4( 4,j,igrid))*tval2+w4( 5,j,igrid))*tval2+w4( 6,j,igrid))*tval2+ &
&                 w4( 7,j,igrid))*tval2+w4( 8,j,igrid))*tval2+w4( 9,j,igrid))*tval2+ &
&                 w4(10,j,igrid))*tval2+w4(11,j,igrid))*tval2+w4(12,j,igrid))*tval2+ &
&                 w4(13,j,igrid)
            enddo
          else
            tinvrt=one/sqrt(tval)
            tinv=tinvrt*tinvrt
            do j= 1,4
              trys(j)= tgh4(j)*tinv
              wrys(j)= wgh4(j)*tinvrt
            enddo
          endif
        case (5)
          if(tval < 61.0D+00)then
!     call rys(tval,nroots,trys,wrys)
!     return
            igrid=nint(tval*half)
            tval2=tval-igrid*two
            do j= 1,5
              trys(j)=((((((((((( &
&                 t5( 1,j,igrid) *tval2+t5( 2,j,igrid))*tval2+t5( 3,j,igrid))*tval2+ &
&                 t5( 4,j,igrid))*tval2+t5( 5,j,igrid))*tval2+t5( 6,j,igrid))*tval2+ &
&                 t5( 7,j,igrid))*tval2+t5( 8,j,igrid))*tval2+t5( 9,j,igrid))*tval2+ &
&                 t5(10,j,igrid))*tval2+t5(11,j,igrid))*tval2+t5(12,j,igrid))*tval2+ &
&                 t5(13,j,igrid)
            enddo
            do j= 1,5
              wrys(j)=((((((((((( &
&                 w5( 1,j,igrid) *tval2+w5( 2,j,igrid))*tval2+w5( 3,j,igrid))*tval2+ &
&                 w5( 4,j,igrid))*tval2+w5( 5,j,igrid))*tval2+w5( 6,j,igrid))*tval2+ &
&                 w5( 7,j,igrid))*tval2+w5( 8,j,igrid))*tval2+w5( 9,j,igrid))*tval2+ &
&                 w5(10,j,igrid))*tval2+w5(11,j,igrid))*tval2+w5(12,j,igrid))*tval2+ &
&                 w5(13,j,igrid)
            enddo
          else
            tinvrt=one/sqrt(tval)
            tinv=tinvrt*tinvrt
            do j= 1,5
              trys(j)= tgh5(j)*tinv
              wrys(j)= wgh5(j)*tinvrt
            enddo
          endif
        case (6)
          if(tval < 67.5D+00)then
            igrid=nint(tval*third)
            tval2=tval-igrid*three
            do j= 1,6
              trys(j)=((((((((((((( &
&                 t6( 1,j,igrid) *tval2+t6( 2,j,igrid))*tval2+t6( 3,j,igrid))*tval2+ &
&                 t6( 4,j,igrid))*tval2+t6( 5,j,igrid))*tval2+t6( 6,j,igrid))*tval2+ &
&                 t6( 7,j,igrid))*tval2+t6( 8,j,igrid))*tval2+t6( 9,j,igrid))*tval2+ &
&                 t6(10,j,igrid))*tval2+t6(11,j,igrid))*tval2+t6(12,j,igrid))*tval2+ &
&                 t6(13,j,igrid))*tval2+t6(14,j,igrid))*tval2+t6(15,j,igrid)
            enddo
            do j= 1,6
              wrys(j)=((((((((((((( &
&                 w6( 1,j,igrid) *tval2+w6( 2,j,igrid))*tval2+w6( 3,j,igrid))*tval2+ &
&                 w6( 4,j,igrid))*tval2+w6( 5,j,igrid))*tval2+w6( 6,j,igrid))*tval2+ &
&                 w6( 7,j,igrid))*tval2+w6( 8,j,igrid))*tval2+w6( 9,j,igrid))*tval2+ &
&                 w6(10,j,igrid))*tval2+w6(11,j,igrid))*tval2+w6(12,j,igrid))*tval2+ &
&                 w6(13,j,igrid))*tval2+w6(14,j,igrid))*tval2+w6(15,j,igrid)
            enddo
          else
            tinvrt=one/sqrt(tval)
            tinv=tinvrt*tinvrt
            do j= 1,6
              trys(j)= tgh6(j)*tinv
              wrys(j)= wgh6(j)*tinvrt
            enddo
          endif
        case (7)
          if(tval < 73.5D+00)then
            igrid=nint(tval*third)
            tval2=tval-igrid*three
            do j= 1,7
              trys(j)=((((((((((((( &
&                 t7( 1,j,igrid) *tval2+t7( 2,j,igrid))*tval2+t7( 3,j,igrid))*tval2+ &
&                 t7( 4,j,igrid))*tval2+t7( 5,j,igrid))*tval2+t7( 6,j,igrid))*tval2+ &
&                 t7( 7,j,igrid))*tval2+t7( 8,j,igrid))*tval2+t7( 9,j,igrid))*tval2+ &
&                 t7(10,j,igrid))*tval2+t7(11,j,igrid))*tval2+t7(12,j,igrid))*tval2+ &
&                 t7(13,j,igrid))*tval2+t7(14,j,igrid))*tval2+t7(15,j,igrid)
            enddo
            do j= 1,7
              wrys(j)=((((((((((((( &
&                 w7( 1,j,igrid) *tval2+w7( 2,j,igrid))*tval2+w7( 3,j,igrid))*tval2+ &
&                 w7( 4,j,igrid))*tval2+w7( 5,j,igrid))*tval2+w7( 6,j,igrid))*tval2+ &
&                 w7( 7,j,igrid))*tval2+w7( 8,j,igrid))*tval2+w7( 9,j,igrid))*tval2+ &
&                 w7(10,j,igrid))*tval2+w7(11,j,igrid))*tval2+w7(12,j,igrid))*tval2+ &
&                 w7(13,j,igrid))*tval2+w7(14,j,igrid))*tval2+w7(15,j,igrid)
            enddo
          else
            tinvrt=one/sqrt(tval)
            tinv=tinvrt*tinvrt
            do j= 1,7
              trys(j)= tgh7(j)*tinv
              wrys(j)= wgh7(j)*tinvrt
            enddo
          endif
        case (8)
          if(tval < 76.5D+00)then
            igrid=nint(tval*third)
            tval2=tval-igrid*three
            do j= 1,8
              trys(j)=((((((((((((( &
&                 t8( 1,j,igrid) *tval2+t8( 2,j,igrid))*tval2+t8( 3,j,igrid))*tval2+ &
&                 t8( 4,j,igrid))*tval2+t8( 5,j,igrid))*tval2+t8( 6,j,igrid))*tval2+ &
&                 t8( 7,j,igrid))*tval2+t8( 8,j,igrid))*tval2+t8( 9,j,igrid))*tval2+ &
&                 t8(10,j,igrid))*tval2+t8(11,j,igrid))*tval2+t8(12,j,igrid))*tval2+ &
&                 t8(13,j,igrid))*tval2+t8(14,j,igrid))*tval2+t8(15,j,igrid)
            enddo
            do j= 1,8
              wrys(j)=((((((((((((( &
&                 w8( 1,j,igrid) *tval2+w8( 2,j,igrid))*tval2+w8( 3,j,igrid))*tval2+ &
&                 w8( 4,j,igrid))*tval2+w8( 5,j,igrid))*tval2+w8( 6,j,igrid))*tval2+ &
&                 w8( 7,j,igrid))*tval2+w8( 8,j,igrid))*tval2+w8( 9,j,igrid))*tval2+ &
&                 w8(10,j,igrid))*tval2+w8(11,j,igrid))*tval2+w8(12,j,igrid))*tval2+ &
&                 w8(13,j,igrid))*tval2+w8(14,j,igrid))*tval2+w8(15,j,igrid)
            enddo
          else
            tinvrt=one/sqrt(tval)
            tinv=tinvrt*tinvrt
            do j= 1,8
              trys(j)= tgh8(j)*tinv
              wrys(j)= wgh8(j)*tinvrt
            enddo
          endif
        case (9)
          if(tval < 82.5D+00)then
            igrid=nint(tval*fifth)
            tval2=tval-igrid*five
            do j= 1,9
              trys(j)=(((((((((((((((( &
&                 t9( 1,j,igrid) *tval2+t9( 2,j,igrid))*tval2+t9( 3,j,igrid))*tval2+ &
&                 t9( 4,j,igrid))*tval2+t9( 5,j,igrid))*tval2+t9( 6,j,igrid))*tval2+ &
&                 t9( 7,j,igrid))*tval2+t9( 8,j,igrid))*tval2+t9( 9,j,igrid))*tval2+ &
&                 t9(10,j,igrid))*tval2+t9(11,j,igrid))*tval2+t9(12,j,igrid))*tval2+ &
&                 t9(13,j,igrid))*tval2+t9(14,j,igrid))*tval2+t9(15,j,igrid))*tval2+ &
&                 t9(16,j,igrid))*tval2+t9(17,j,igrid))*tval2+t9(18,j,igrid)
            enddo
            do j= 1,9
              wrys(j)=(((((((((((((((( &
&                 w9( 1,j,igrid) *tval2+w9( 2,j,igrid))*tval2+w9( 3,j,igrid))*tval2+ &
&                 w9( 4,j,igrid))*tval2+w9( 5,j,igrid))*tval2+w9( 6,j,igrid))*tval2+ &
&                 w9( 7,j,igrid))*tval2+w9( 8,j,igrid))*tval2+w9( 9,j,igrid))*tval2+ &
&                 w9(10,j,igrid))*tval2+w9(11,j,igrid))*tval2+w9(12,j,igrid))*tval2+ &
&                 w9(13,j,igrid))*tval2+w9(14,j,igrid))*tval2+w9(15,j,igrid))*tval2+ &
&                 w9(16,j,igrid))*tval2+w9(17,j,igrid))*tval2+w9(18,j,igrid)
            enddo
          else
            tinvrt=one/sqrt(tval)
            tinv=tinvrt*tinvrt
            do j= 1,9
              trys(j)= tgh9(j)*tinv
              wrys(j)= wgh9(j)*tinvrt
            enddo
          endif
        case default
          write(*,'("Subroutine rysroot supports up to nroots=9")')
          call iabort
      end select
      return
end
