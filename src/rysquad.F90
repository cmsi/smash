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
      use rysgrid10, only : t10, w10
      use rysgrid11, only : t11, w11
      use rysgrid12, only : t12, w12
      use rysgrid13, only : t13, w13
      implicit none
      integer,intent(in) :: nroots
      integer :: j, igrid
      real(8),parameter :: half=0.5D+00, one=1.0D+00, two=2.0D+00, three=3.0D+00
      real(8),parameter :: third=0.3333333333333333D+00, fifth=0.2D+00, five=5.0D+00
      real(8),parameter :: twenty=20.0D+00
      real(8),intent(in) :: tval
      real(8),intent(out) :: trys(13), wrys(13)
      real(8) :: tgh1, tgh2(2), tgh3(3), tgh4(4), tgh5(5), tgh6(6), tgh7(7), tgh8(8), tgh9(9)
      real(8) :: tgh10(10), tgh11(11), tgh12(12), tgh13(13)
      real(8) :: wgh1, wgh2(2), wgh3(3), wgh4(4), wgh5(5), wgh6(6), wgh7(7), wgh8(8), wgh9(9)
      real(8) :: wgh10(10), wgh11(11), wgh12(12), wgh13(13)
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
      data tgh10/6.019206314958792D-02, 5.438675002946460D-01, 1.522944105404444D+00, &
&                3.022513376451574D+00, 5.084907750098524D+00, 7.777439231525445D+00, &
&                1.120813020434866D+01, 1.556116333218935D+01, 2.119389209630154D+01, &
&                2.902495034023623D+01/
      data tgh11/5.483986957881849D-02, 4.951741233503564D-01, 1.384655740084600D+00, &
&                2.741919940106702D+00, 4.597737700485711D+00, 6.999397469528836D+00, &
&                1.001890827595723D+01, 1.376930586610169D+01, 1.844111968097819D+01, &
&                2.440196124238704D+01, 3.259498009144082D+01/
      data tgh12/5.036188911729395D-02, 4.545066815637803D-01, 1.269589940103961D+00, &
&                2.509848097232128D+00, 4.198415644878413D+00, 6.369975388030635D+00, &
&                9.075434230961203D+00, 1.239044796380947D+01, 1.643219508767531D+01, &
&                2.139675593616611D+01, 2.766110877984609D+01, 3.619136036061560D+01/
      data tgh13/4.656008324502477D-02, 4.200274064012136D-01, 1.172310773277780D+00, &
&                2.314540864349434D+00, 3.864585038228159D+00, 5.848734811306343D+00, &
&                8.304553489985900D+00, 1.128575099351764D+01, 1.487096037752540D+01, &
&                1.918091948561046D+01, 2.441669233305652D+01, 3.096393827474680D+01, &
&                3.981042606874934D+01/
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
      data wgh10/4.622436696006101D-01, 2.866755053628341D-01, 1.090172060200233D-01, &
&                2.481052088746361D-02, 3.243773342237862D-03, 2.283386360163540D-04, &
&                7.802556478532064D-06, 1.086069370769282D-07, 4.399340992273181D-10, &
&                2.229393645534151D-13/
      data wgh11/4.435452264349593D-01, 2.869714332469071D-01, 1.191023609587825D-01, &
&                3.114037088442385D-02, 4.978399335051647D-03, 4.648850508842522D-04, &
&                2.365512855251046D-05, 5.884287563301006D-07, 5.966990986059653D-09, &
&                1.744339007547993D-11, 6.167183424404049D-15/
      data wgh12/4.269311638686992D-01, 2.861795353464430D-01, 1.277396217845592D-01, &
&                3.744547050323075D-02, 7.048355810072671D-03, 8.236924826884175D-04, &
&                5.688691636404380D-05, 2.158245704902334D-06, 4.018971174941430D-08, &
&                3.046254269987564D-10, 6.584620243078170D-13, 1.664368496489109D-16/
      data wgh13/4.120436505903693D-01, 2.846322411767845D-01, 1.351133279117879D-01, &
&                4.359822721725080D-02, 9.397901291159580D-03, 1.319064722323857D-03, &
&                1.162297016031096D-04, 6.103291717396013D-06, 1.770106337397356D-07, &
&                2.524494034490534D-09, 1.460999933981612D-11, 2.383148659372161D-14, &
&                4.396916094753902D-18/
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
        case (10)
          if(tval < 90.0D+00)then
            igrid=nint(tval/twenty)
            tval2=tval-igrid*twenty
            do j= 1,10
              trys(j)=((((((((((((((((((((((((((((((( &
&                 t10( 1,j,igrid) *tval2+t10( 2,j,igrid))*tval2+t10( 3,j,igrid))*tval2+ &
&                 t10( 4,j,igrid))*tval2+t10( 5,j,igrid))*tval2+t10( 6,j,igrid))*tval2+ &
&                 t10( 7,j,igrid))*tval2+t10( 8,j,igrid))*tval2+t10( 9,j,igrid))*tval2+ &
&                 t10(10,j,igrid))*tval2+t10(11,j,igrid))*tval2+t10(12,j,igrid))*tval2+ &
&                 t10(13,j,igrid))*tval2+t10(14,j,igrid))*tval2+t10(15,j,igrid))*tval2+ &
&                 t10(16,j,igrid))*tval2+t10(17,j,igrid))*tval2+t10(18,j,igrid))*tval2+ &
&                 t10(19,j,igrid))*tval2+t10(20,j,igrid))*tval2+t10(21,j,igrid))*tval2+ &
&                 t10(22,j,igrid))*tval2+t10(23,j,igrid))*tval2+t10(24,j,igrid))*tval2+ &
&                 t10(25,j,igrid))*tval2+t10(26,j,igrid))*tval2+t10(27,j,igrid))*tval2+ &
&                 t10(28,j,igrid))*tval2+t10(29,j,igrid))*tval2+t10(30,j,igrid))*tval2+ &
&                 t10(31,j,igrid))*tval2+t10(32,j,igrid))*tval2+t10(33,j,igrid)
            enddo
            do j= 1,10
              wrys(j)=((((((((((((((((((((((((((((((( &
&                 w10( 1,j,igrid) *tval2+w10( 2,j,igrid))*tval2+w10( 3,j,igrid))*tval2+ &
&                 w10( 4,j,igrid))*tval2+w10( 5,j,igrid))*tval2+w10( 6,j,igrid))*tval2+ &
&                 w10( 7,j,igrid))*tval2+w10( 8,j,igrid))*tval2+w10( 9,j,igrid))*tval2+ &
&                 w10(10,j,igrid))*tval2+w10(11,j,igrid))*tval2+w10(12,j,igrid))*tval2+ &
&                 w10(13,j,igrid))*tval2+w10(14,j,igrid))*tval2+w10(15,j,igrid))*tval2+ &
&                 w10(16,j,igrid))*tval2+w10(17,j,igrid))*tval2+w10(18,j,igrid))*tval2+ &
&                 w10(19,j,igrid))*tval2+w10(20,j,igrid))*tval2+w10(21,j,igrid))*tval2+ &
&                 w10(22,j,igrid))*tval2+w10(23,j,igrid))*tval2+w10(24,j,igrid))*tval2+ &
&                 w10(25,j,igrid))*tval2+w10(26,j,igrid))*tval2+w10(27,j,igrid))*tval2+ &
&                 w10(28,j,igrid))*tval2+w10(29,j,igrid))*tval2+w10(30,j,igrid))*tval2+ &
&                 w10(31,j,igrid))*tval2+w10(32,j,igrid))*tval2+w10(33,j,igrid)
            enddo
          else
            tinvrt=one/sqrt(tval)
            tinv=tinvrt*tinvrt
            do j= 1,10
              trys(j)= tgh10(j)*tinv
              wrys(j)= wgh10(j)*tinvrt
            enddo
          endif
        case (11)
          if(tval < 90.0D+00)then
            igrid=nint(tval/twenty)
            tval2=tval-igrid*twenty
            do j= 1,11
              trys(j)=((((((((((((((((((((((((((((((( &
&                 t11( 1,j,igrid) *tval2+t11( 2,j,igrid))*tval2+t11( 3,j,igrid))*tval2+ &
&                 t11( 4,j,igrid))*tval2+t11( 5,j,igrid))*tval2+t11( 6,j,igrid))*tval2+ &
&                 t11( 7,j,igrid))*tval2+t11( 8,j,igrid))*tval2+t11( 9,j,igrid))*tval2+ &
&                 t11(10,j,igrid))*tval2+t11(11,j,igrid))*tval2+t11(12,j,igrid))*tval2+ &
&                 t11(13,j,igrid))*tval2+t11(14,j,igrid))*tval2+t11(15,j,igrid))*tval2+ &
&                 t11(16,j,igrid))*tval2+t11(17,j,igrid))*tval2+t11(18,j,igrid))*tval2+ &
&                 t11(19,j,igrid))*tval2+t11(20,j,igrid))*tval2+t11(21,j,igrid))*tval2+ &
&                 t11(22,j,igrid))*tval2+t11(23,j,igrid))*tval2+t11(24,j,igrid))*tval2+ &
&                 t11(25,j,igrid))*tval2+t11(26,j,igrid))*tval2+t11(27,j,igrid))*tval2+ &
&                 t11(28,j,igrid))*tval2+t11(29,j,igrid))*tval2+t11(30,j,igrid))*tval2+ &
&                 t11(31,j,igrid))*tval2+t11(32,j,igrid))*tval2+t11(33,j,igrid)
            enddo
            do j= 1,11
              wrys(j)=((((((((((((((((((((((((((((((( &
&                 w11( 1,j,igrid) *tval2+w11( 2,j,igrid))*tval2+w11( 3,j,igrid))*tval2+ &
&                 w11( 4,j,igrid))*tval2+w11( 5,j,igrid))*tval2+w11( 6,j,igrid))*tval2+ &
&                 w11( 7,j,igrid))*tval2+w11( 8,j,igrid))*tval2+w11( 9,j,igrid))*tval2+ &
&                 w11(10,j,igrid))*tval2+w11(11,j,igrid))*tval2+w11(12,j,igrid))*tval2+ &
&                 w11(13,j,igrid))*tval2+w11(14,j,igrid))*tval2+w11(15,j,igrid))*tval2+ &
&                 w11(16,j,igrid))*tval2+w11(17,j,igrid))*tval2+w11(18,j,igrid))*tval2+ &
&                 w11(19,j,igrid))*tval2+w11(20,j,igrid))*tval2+w11(21,j,igrid))*tval2+ &
&                 w11(22,j,igrid))*tval2+w11(23,j,igrid))*tval2+w11(24,j,igrid))*tval2+ &
&                 w11(25,j,igrid))*tval2+w11(26,j,igrid))*tval2+w11(27,j,igrid))*tval2+ &
&                 w11(28,j,igrid))*tval2+w11(29,j,igrid))*tval2+w11(30,j,igrid))*tval2+ &
&                 w11(31,j,igrid))*tval2+w11(32,j,igrid))*tval2+w11(33,j,igrid)
            enddo
          else
            tinvrt=one/sqrt(tval)
            tinv=tinvrt*tinvrt
            do j= 1,11
              trys(j)= tgh11(j)*tinv
              wrys(j)= wgh11(j)*tinvrt
            enddo
          endif
        case (12)
          if(tval < 110.0D+00)then
            igrid=nint(tval/twenty)
            tval2=tval-igrid*twenty
            do j= 1,12
              trys(j)=((((((((((((((((((((((((((((((( &
&                 t12( 1,j,igrid) *tval2+t12( 2,j,igrid))*tval2+t12( 3,j,igrid))*tval2+ &
&                 t12( 4,j,igrid))*tval2+t12( 5,j,igrid))*tval2+t12( 6,j,igrid))*tval2+ &
&                 t12( 7,j,igrid))*tval2+t12( 8,j,igrid))*tval2+t12( 9,j,igrid))*tval2+ &
&                 t12(10,j,igrid))*tval2+t12(11,j,igrid))*tval2+t12(12,j,igrid))*tval2+ &
&                 t12(13,j,igrid))*tval2+t12(14,j,igrid))*tval2+t12(15,j,igrid))*tval2+ &
&                 t12(16,j,igrid))*tval2+t12(17,j,igrid))*tval2+t12(18,j,igrid))*tval2+ &
&                 t12(19,j,igrid))*tval2+t12(20,j,igrid))*tval2+t12(21,j,igrid))*tval2+ &
&                 t12(22,j,igrid))*tval2+t12(23,j,igrid))*tval2+t12(24,j,igrid))*tval2+ &
&                 t12(25,j,igrid))*tval2+t12(26,j,igrid))*tval2+t12(27,j,igrid))*tval2+ &
&                 t12(28,j,igrid))*tval2+t12(29,j,igrid))*tval2+t12(30,j,igrid))*tval2+ &
&                 t12(31,j,igrid))*tval2+t12(32,j,igrid))*tval2+t12(33,j,igrid)
            enddo
            do j= 1,12
              wrys(j)=((((((((((((((((((((((((((((((( &
&                 w12( 1,j,igrid) *tval2+w12( 2,j,igrid))*tval2+w12( 3,j,igrid))*tval2+ &
&                 w12( 4,j,igrid))*tval2+w12( 5,j,igrid))*tval2+w12( 6,j,igrid))*tval2+ &
&                 w12( 7,j,igrid))*tval2+w12( 8,j,igrid))*tval2+w12( 9,j,igrid))*tval2+ &
&                 w12(10,j,igrid))*tval2+w12(11,j,igrid))*tval2+w12(12,j,igrid))*tval2+ &
&                 w12(13,j,igrid))*tval2+w12(14,j,igrid))*tval2+w12(15,j,igrid))*tval2+ &
&                 w12(16,j,igrid))*tval2+w12(17,j,igrid))*tval2+w12(18,j,igrid))*tval2+ &
&                 w12(19,j,igrid))*tval2+w12(20,j,igrid))*tval2+w12(21,j,igrid))*tval2+ &
&                 w12(22,j,igrid))*tval2+w12(23,j,igrid))*tval2+w12(24,j,igrid))*tval2+ &
&                 w12(25,j,igrid))*tval2+w12(26,j,igrid))*tval2+w12(27,j,igrid))*tval2+ &
&                 w12(28,j,igrid))*tval2+w12(29,j,igrid))*tval2+w12(30,j,igrid))*tval2+ &
&                 w12(31,j,igrid))*tval2+w12(32,j,igrid))*tval2+w12(33,j,igrid)
            enddo
          else
            tinvrt=one/sqrt(tval)
            tinv=tinvrt*tinvrt
            do j= 1,12
              trys(j)= tgh12(j)*tinv
              wrys(j)= wgh12(j)*tinvrt
            enddo
          endif
        case (13)
          if(tval < 110.0D+00)then
            igrid=nint(tval/twenty)
            tval2=tval-igrid*twenty
            do j= 1,13
              trys(j)=((((((((((((((((((((((((((((((( &
&                 t13( 1,j,igrid) *tval2+t13( 2,j,igrid))*tval2+t13( 3,j,igrid))*tval2+ &
&                 t13( 4,j,igrid))*tval2+t13( 5,j,igrid))*tval2+t13( 6,j,igrid))*tval2+ &
&                 t13( 7,j,igrid))*tval2+t13( 8,j,igrid))*tval2+t13( 9,j,igrid))*tval2+ &
&                 t13(10,j,igrid))*tval2+t13(11,j,igrid))*tval2+t13(12,j,igrid))*tval2+ &
&                 t13(13,j,igrid))*tval2+t13(14,j,igrid))*tval2+t13(15,j,igrid))*tval2+ &
&                 t13(16,j,igrid))*tval2+t13(17,j,igrid))*tval2+t13(18,j,igrid))*tval2+ &
&                 t13(19,j,igrid))*tval2+t13(20,j,igrid))*tval2+t13(21,j,igrid))*tval2+ &
&                 t13(22,j,igrid))*tval2+t13(23,j,igrid))*tval2+t13(24,j,igrid))*tval2+ &
&                 t13(25,j,igrid))*tval2+t13(26,j,igrid))*tval2+t13(27,j,igrid))*tval2+ &
&                 t13(28,j,igrid))*tval2+t13(29,j,igrid))*tval2+t13(30,j,igrid))*tval2+ &
&                 t13(31,j,igrid))*tval2+t13(32,j,igrid))*tval2+t13(33,j,igrid)
            enddo
            do j= 1,13
              wrys(j)=((((((((((((((((((((((((((((((( &
&                 w13( 1,j,igrid) *tval2+w13( 2,j,igrid))*tval2+w13( 3,j,igrid))*tval2+ &
&                 w13( 4,j,igrid))*tval2+w13( 5,j,igrid))*tval2+w13( 6,j,igrid))*tval2+ &
&                 w13( 7,j,igrid))*tval2+w13( 8,j,igrid))*tval2+w13( 9,j,igrid))*tval2+ &
&                 w13(10,j,igrid))*tval2+w13(11,j,igrid))*tval2+w13(12,j,igrid))*tval2+ &
&                 w13(13,j,igrid))*tval2+w13(14,j,igrid))*tval2+w13(15,j,igrid))*tval2+ &
&                 w13(16,j,igrid))*tval2+w13(17,j,igrid))*tval2+w13(18,j,igrid))*tval2+ &
&                 w13(19,j,igrid))*tval2+w13(20,j,igrid))*tval2+w13(21,j,igrid))*tval2+ &
&                 w13(22,j,igrid))*tval2+w13(23,j,igrid))*tval2+w13(24,j,igrid))*tval2+ &
&                 w13(25,j,igrid))*tval2+w13(26,j,igrid))*tval2+w13(27,j,igrid))*tval2+ &
&                 w13(28,j,igrid))*tval2+w13(29,j,igrid))*tval2+w13(30,j,igrid))*tval2+ &
&                 w13(31,j,igrid))*tval2+w13(32,j,igrid))*tval2+w13(33,j,igrid)
            enddo
          else
            tinvrt=one/sqrt(tval)
            tinv=tinvrt*tinvrt
            do j= 1,13
              trys(j)= tgh13(j)*tinv
              wrys(j)= wgh13(j)*tinvrt
            enddo
          endif
        case default
          write(*,'("Subroutine rysroot supports up to nroots=13")')
          call exit
      end select
      return
end
