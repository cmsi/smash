!--------------------------------------------------------------------------------------
  subroutine setecp(term0ecp,term1ecp,term2ecp,label1ecp,num1ecp,label2ecp,num2ecp, &
&                   xyzintecp,maxecpdim,maxpangdim)
!--------------------------------------------------------------------------------------
!
! Prepare values for ECP calculation
!
      use ecp, only : nterm1, nterm2
      implicit none
      integer,parameter :: lenecp(0:6)=(/1,4,10,20,35,56,84/)
      integer,intent(in) :: maxecpdim, maxpangdim
      integer,intent(out) :: label1ecp(9*nterm1), label2ecp(6*nterm2)
      integer,intent(out) :: num1ecp(2,56,56), num2ecp(2,maxpangdim*maxpangdim*lenecp(maxecpdim))
      real(8),intent(out) :: term0ecp(maxpangdim*maxpangdim*lenecp(maxecpdim))
      real(8),intent(out) :: term1ecp(nterm1), term2ecp(nterm2), xyzintecp(25*25*25)
!
      call ecpzlm
      call ecpxyz(xyzintecp,maxecpdim)
      call ecpangint0(term0ecp,xyzintecp,maxecpdim,maxpangdim)
      call ecpangint1(term1ecp,label1ecp,num1ecp,xyzintecp,maxecpdim)
      call ecpangint2(term2ecp,label2ecp,num2ecp,xyzintecp,maxecpdim,maxpangdim)
!
      return
end


!-----------------------------------
! subroutine intecp(hmat,ish,jsh,maxdim)
!-----------------------------------
!
!
!





!--------------------
  subroutine ecpzlm
!--------------------
!
! Set up the coefficients for the real spherical harmonics
!
      use ecp, only : zlm
      implicit none
      integer :: i
      real(8),parameter :: one=1.0D+00, two=2.0D+00, three=3.0D+00, four=4.0D+00
      real(8),parameter :: five=5.0D+00, six=6.0D+00, seven=7.0D+00, eight=8.0D+00
      real(8),parameter :: nine=9.0D+00, ten=10.0D+00, eleven=11.0D+00, p13=13.0D+00
      real(8),parameter :: p14=14.0D+00, p15=15.0D+00, p17=17.0D+00, p19=19.0D+00
      real(8),parameter :: p20=20.0D+00, half=0.5D+00, pi4=1.256637061435917D+01
      real(8),parameter :: sr3=1.732050807568877D+00, sr5=2.236067977499790D+00
      real(8),parameter :: sr7=2.645751311064591D+00, sr11=3.316624790355400D+00
      real(8),parameter :: sr13=3.605551275463989D+00, sr15=3.872983346207417D+00
      real(8),parameter :: sr17=4.123105625617661D+00, sr19=4.358898943540674D+00
      real(8),parameter :: sr21=4.582575694955840D+00, sr35=5.916079783099616D+00
      real(8),parameter :: sr55=7.416198487095663D+00, sr77=8.774964387392122D+00
      real(8),parameter :: sr91=9.539392014169456D+00, sr119=1.090871211463571D+01
      real(8),parameter :: sr221=1.486606874731851D+01, sr247=1.571623364550171D+01
      real(8),parameter :: sr323=1.797220075561143D+01, sr385=1.962141687034858D+01
      real(8),parameter :: sr1001=3.163858403911275D+01, sr1309=3.618010503025109D+01
      real(8),parameter :: sr1365=3.694590640382233D+01, sr2431=4.930517214248420D+01
      real(8),parameter :: sr3003=5.479963503528103D+01, sr5005=7.074602462329597D+01
      real(8),parameter :: s13585=1.1655470818461174D+02
      real(8) :: fap(0:17), fat(0:17), tmp
!
      fap(0)= sqrt(one/pi4)
      fap(1)= fap(0)*sqrt(half)
      do i= 2,16,2
         fap(i  )= fap(i-2)*half
         fap(i+1)= fap(i-1)*half
      enddo
      do i= 0,17
        fat(i)= fap(i)*three
      enddo
!
! L= 0, ML= 0
!
      zlm(  1)= fap(0)
!
! L= 1, ML=+1,0,-1
!
      zlm(  2)= sr3*fap(0)
      zlm(  3)= zlm(  2)
      zlm(  4)= zlm(  2)
!
! L= 2, ML=+2...-2
!
      zlm(  5)= sr15*fap(2)
      zlm(  6)=-zlm(  5)
      zlm(  7)= zlm(  5)*two
      tmp     = sr5*fap(2)
      zlm(  8)= tmp*three
      zlm(  9)=-tmp
      zlm( 10)= zlm(  7)
      zlm( 11)= zlm(  7)
!
! L= 3, ML=+3...-3
!
      zlm( 12)= sr35*fap(3)
      zlm( 13)=-zlm( 12)*three
      zlm( 14)= sr21*sr5*fap(2)
      zlm( 15)=-zlm( 14)
      tmp     = sr21*fap(3)
      zlm( 16)= tmp*five
      zlm( 17)=-tmp
      tmp     = sr7*fap(2)
      zlm( 18)= tmp*five
      zlm( 19)=-tmp*three
      zlm( 20)= zlm( 16)
      zlm( 21)= zlm( 17)
      zlm( 22)= zlm( 14)*two
      zlm( 23)=-zlm( 13)
      zlm( 24)=-zlm( 12)
!
! L= 4, ML=+4...-4
!
      zlm( 25)= sr35*fat(6)
      zlm( 26)=-zlm( 25)*six
      zlm( 27)= zlm( 25)
      zlm( 28)= sr35*fat(3)
      zlm( 29)=-zlm( 28)*three
      tmp     = sr5*fat(4)
      zlm( 30)= tmp*seven
      zlm( 31)=-zlm( 30)
      zlm( 32)= tmp
      zlm( 33)=-tmp
      tmp     = sr5*fat(3)
      zlm( 34)= tmp*seven
      zlm( 35)=-tmp*three
      tmp     = fat(6)
      zlm( 36)= tmp*five*seven
      zlm( 37)=-tmp*five*six
      zlm( 38)= tmp*three
      zlm( 39)= zlm( 34)
      zlm( 40)= zlm( 35)
      tmp     = sr5*fat(2)
      zlm( 41)= tmp*seven
      zlm( 42)=-tmp
      zlm( 43)=-zlm( 29)
      zlm( 44)=-zlm( 28)
      zlm( 45)= sr35*fat(2)
      zlm( 46)=-zlm( 45)
!
! L= 5, ML=+5...-5
!
      zlm( 47)= sr77*fat(7)
      zlm( 48)=-zlm( 47)*ten
      zlm( 49)= zlm( 47)*five
      zlm( 50)= sr385*fat(6)
      zlm( 51)=-zlm( 50)*six
      zlm( 52)= zlm( 50)
      tmp     = sr385*fap(7)
      zlm( 53)= tmp*three*three
      zlm( 54)=-tmp*three*nine
      zlm( 55)= tmp*three
      zlm( 56)=-tmp
      tmp     = sr385*sr3*fap(4)
      zlm( 57)= tmp*three
      zlm( 58)=-zlm( 57)
      zlm( 59)=-tmp
      zlm( 60)= tmp
      tmp     = sr15*sr11*fap(6)
      zlm( 61)= tmp*2.1D+01
      zlm( 62)=-tmp*p14
      zlm( 63)= tmp
      tmp     = sr11*fap(6)
      zlm( 64)= tmp*6.3D+01
      zlm( 65)=-tmp*7.0D+01
      zlm( 66)= tmp*p15
      zlm( 67)= zlm( 61)
      zlm( 68)= zlm( 62)
      zlm( 69)= zlm( 63)
      tmp     = sr385*sr3*fap(2)
      zlm( 70)= tmp*three
      zlm( 71)=-tmp
      zlm( 72)=-zlm( 54)
      zlm( 73)=-zlm( 55)
      zlm( 74)=-zlm( 53)
      zlm( 75)=-zlm( 56)
      zlm( 76)= sr385*fat(2)
      zlm( 77)=-zlm( 76)
      zlm( 78)= zlm( 49)
      zlm( 79)= zlm( 48)
      zlm( 80)= zlm( 47)
!
! L= 6, ML=+6...-6
!
      zlm( 81)= sr3003*fap(9)
      zlm( 82)=-zlm( 81)*p15
      zlm( 83)=-zlm( 82)
      zlm( 84)=-zlm( 81)
      zlm( 85)= sr1001*fat(7)
      zlm( 86)=-zlm( 85)*ten
      zlm( 87)= zlm( 85)*five
      tmp     = sr91*fat(8)
      zlm( 88)= tmp*eleven
      zlm( 89)=-tmp*eleven*six
      zlm( 90)= zlm( 88)
      zlm( 91)=-tmp
      zlm( 92)= tmp*six
      zlm( 93)=-tmp
      tmp     = sr1365*fap(7)
      zlm( 94)= tmp*eleven
      zlm( 95)=-tmp*three*eleven
      zlm( 96)= tmp*three*three
      zlm( 97)=-tmp*three
      tmp     = sr1365*fap(9)
      zlm( 98)= tmp*three*eleven
      zlm( 99)=-zlm( 98)
      zlm(100)=-tmp*three*six
      zlm(101)= tmp*three*six
      zlm(102)= tmp
      zlm(103)=-tmp
      tmp     = sr21*sr13*fap(6)
      zlm(104)= tmp*three*eleven
      zlm(105)=-tmp*three*ten
      zlm(106)= tmp*five
      tmp     = sr13*fap(8)
      zlm(107)= tmp*2.31D+02
      zlm(108)=-tmp*3.15D+02
      zlm(109)= tmp*1.05D+02
      zlm(110)=-tmp*five
      zlm(111)= zlm(104)
      zlm(112)= zlm(105)
      zlm(113)= zlm(106)
      tmp     = sr1365*fap(7)
      zlm(114)= tmp*three*eleven
      zlm(115)=-tmp*three*six
      zlm(116)= tmp
      zlm(117)=-zlm( 95)
      zlm(118)=-zlm( 94)
      zlm(119)=-zlm( 96)
      zlm(120)=-zlm( 97)
      tmp     = sr91*fat(4)
      zlm(121)= tmp*eleven
      zlm(122)=-zlm(121)
      zlm(123)=-tmp
      zlm(124)= tmp
      zlm(125)= zlm( 87)
      zlm(126)= zlm( 86)
      zlm(127)= zlm( 85)
      tmp     = sr3003*fap(9)
      zlm(128)= tmp*six
      zlm(129)=-tmp*p20
      zlm(130)= zlm(128)
!
! L= 7, ML=-7
      zlm(131)= sr55*sr13*fat(10)
      zlm(132)= zlm(131)*seven*five
      zlm(133)=-zlm(131)*seven*three
      zlm(134)=-zlm(131)*seven
! L= 7, ML=-6
      zlm(135)= sr5005*fat(9)
      zlm(136)=-zlm(135)*p15
      zlm(137)=-zlm(135)
      zlm(138)=-zlm(136)
! L= 7, ML=-5
      tmp     = sr385*fat(10)
      zlm(139)= tmp*p13
      zlm(140)=-zlm(139)*ten
      zlm(141)= zlm(139)*five
      zlm(142)=-tmp
      zlm(143)=-zlm(142)*ten
      zlm(144)= zlm(142)*five
! L= 7, ML=-4
      tmp     = sr385*fat(8)
      zlm(145)= tmp*p13
      zlm(146)=-tmp*p13*six
      zlm(147)= zlm(145)
      zlm(148)=-tmp*three
      zlm(149)= tmp*three*six
      zlm(150)= zlm(148)
! L= 7, ML=-3
      tmp     = sr35*fat(10)
      zlm(151)= tmp*eleven*p13
      zlm(152)=-zlm(151)*three
      zlm(153)=-tmp*eleven*six
      zlm(154)=-zlm(153)*three
      zlm(155)= tmp*three
      zlm(156)=-tmp*nine
! L= 7, ML=-2
      tmp     = sr35*fat(9)
      zlm(157)= tmp*eleven*p13
      zlm(158)=-zlm(157)
      zlm(159)=-tmp*eleven*ten
      zlm(160)=-zlm(159)
      zlm(161)= tmp*p15
      zlm(162)=-zlm(161)
! L= 7, ML=-1
      tmp     = sr7*sr15*fap(10)
      zlm(163)= tmp*4.29D+02
      zlm(164)=-tmp*4.95D+02
      zlm(165)= tmp*1.35D+02
      zlm(166)=-tmp*five
! L= 7, ML= 0
      tmp     = sr15*fap(8)
      zlm(167)= tmp*4.29D+02
      zlm(168)=-tmp*6.93D+02
      zlm(169)= tmp*3.15D+02
      zlm(170)=-tmp*3.5D+01
! L= 7, ML= 1
      zlm(171)= zlm(163)
      zlm(172)= zlm(164)
      zlm(173)= zlm(165)
      zlm(174)= zlm(166)
! L= 7, ML= 2
      zlm(175)= zlm(157)*two
      zlm(176)= zlm(159)*two
      zlm(177)= zlm(161)*two
! L= 7, ML= 3
      zlm(178)=-zlm(152)
      zlm(179)=-zlm(151)
      zlm(180)=-zlm(154)
      zlm(181)=-zlm(153)
      zlm(182)=-zlm(156)
      zlm(183)=-zlm(155)
! L= 7, ML= 4
      zlm(184)= zlm(145)*four
      zlm(185)=-zlm(184)
      zlm(186)= zlm(148)*four
      zlm(187)=-zlm(186)
! L= 7, ML= 5
      zlm(188)= zlm(141)
      zlm(189)= zlm(140)
      zlm(190)= zlm(139)
      zlm(191)= zlm(144)
      zlm(192)= zlm(143)
      zlm(193)= zlm(142)
! L= 7, ML= 6
      tmp     = sr5005*fat(7)
      zlm(194)= tmp*three
      zlm(195)=-tmp*ten
      zlm(196)= zlm(194)
! L= 7, ML= 7
      zlm(197)=-zlm(131)
      zlm(198)=-zlm(132)
      zlm(199)=-zlm(133)
      zlm(200)=-zlm(134)
!
! L= 8, ML=-8
      zlm(201)= sr2431*sr5*fat(14)
      zlm(202)=-zlm(201)*2.8D+01
      zlm(203)= zlm(201)*7.0D+01
      zlm(204)= zlm(202)
      zlm(205)= zlm(201)
! L= 8, ML=-7
      zlm(206)= sr2431*sr5*fat(10)
      zlm(207)=-zlm(206)*seven*three
      zlm(208)= zlm(206)*seven*five
      zlm(209)=-zlm(206)*seven
! L= 8, ML=-6
      tmp     = sr2431*sr3*fap(11)
      zlm(210)= tmp*p15
      zlm(211)=-zlm(210)*p15
      zlm(212)=-zlm(211)
      zlm(213)=-zlm(210)
      zlm(214)=-tmp
      zlm(215)= zlm(210)
      zlm(216)=-zlm(215)
      zlm(217)= tmp
! L= 8, ML=-5
      tmp     = sr1001*sr17*fat(10)
      zlm(218)= tmp*five
      zlm(219)=-zlm(218)*ten
      zlm(220)= zlm(218)*five
      zlm(221)=-tmp
      zlm(222)= tmp*ten
      zlm(223)= zlm(221)*five
! L= 8, ML=-4
      tmp     = sr1309*fat(12)
      zlm(224)= tmp*p13*five
      zlm(225)= zlm(224)
      zlm(226)=-zlm(224)*six
      zlm(227)=-tmp*p13*two
      zlm(228)= zlm(227)
      zlm(229)=-zlm(227)*six
      zlm(230)= tmp
      zlm(231)= tmp
      zlm(232)=-tmp*six
! L= 8, ML=-3
      tmp     = sr1309*sr15*fap(10)
      zlm(233)= tmp*p13*three
      zlm(234)=-zlm(233)*three
      zlm(235)=-tmp*p13*two
      zlm(236)=-zlm(235)*three
      zlm(237)= tmp*three
      zlm(238)=-zlm(237)*three
! L= 8, ML=-2
      tmp     = sr119*sr5*fat(11)
      zlm(239)= tmp*eleven*p13
      zlm(240)=-zlm(239)
      zlm(241)= zlm(240)
      zlm(242)= zlm(239)
      zlm(243)= tmp*eleven*three
      zlm(244)=-zlm(243)
      zlm(245)=-tmp
      zlm(246)= tmp
! L= 8, ML=-1
      tmp     = sr17*fat(10)
      zlm(247)= tmp*7.15D+02
      zlm(248)=-tmp*1.001D+03
      zlm(249)= tmp*3.85D+02
      zlm(250)=-tmp*3.5D+01
! L= 8, ML= 0
      tmp     = sr17*fap(14)
      zlm(251)= tmp*6.435D+03
      zlm(252)=-tmp*1.2012D+04
      zlm(253)= tmp*6.930D+03
      zlm(254)=-tmp*1.260D+03
      zlm(255)= tmp*3.5D+01
! L= 8, ML= 1
      zlm(256)= zlm(247)
      zlm(257)= zlm(248)
      zlm(258)= zlm(249)
      zlm(259)= zlm(250)
! L= 8, ML= 2
      tmp     = sr119*sr5*fat(9)
      zlm(260)= tmp*eleven*p13
      zlm(261)=-zlm(260)
      zlm(262)= tmp*eleven*three
      zlm(263)=-tmp
! L= 8, ML= 3
      zlm(264)=-zlm(234)
      zlm(265)=-zlm(233)
      zlm(266)=-zlm(236)
      zlm(267)=-zlm(235)
      zlm(268)=-zlm(238)
      zlm(269)=-zlm(237)
! L= 8, ML= 4
      tmp     = sr1309*fat(8)
      zlm(270)= tmp*p13*five
      zlm(271)=-zlm(270)
      zlm(272)=-tmp*p13*two
      zlm(273)=-zlm(272)
      zlm(274)= tmp
      zlm(275)=-tmp
! L= 8, ML= 5
      tmp     = sr1001*sr17*fat(10)
      zlm(276)= tmp*five*five
      zlm(277)=-tmp*five*ten
      zlm(278)= tmp*five
      zlm(279)=-tmp*five
      zlm(280)= tmp*ten
      zlm(281)=-tmp
! L= 8, ML= 6
      tmp     = sr2431*sr3*fap(9)
      zlm(282)= tmp*p15*three
      zlm(283)=-tmp*p15*ten
      zlm(284)= zlm(282)
      zlm(285)=-tmp*three
      zlm(286)= tmp*ten
      zlm(287)= zlm(285)
! L= 8, ML= 7
      zlm(288)=-zlm(209)
      zlm(289)=-zlm(208)
      zlm(290)=-zlm(207)
      zlm(291)=-zlm(206)
! L= 8, ML= 8
      zlm(292)= zlm(201)*eight
      zlm(293)=-zlm(292)*seven
      zlm(294)=-zlm(293)
      zlm(295)=-zlm(292)
!
! L= 9, ML=-9
      zlm(296)= s13585*sr17*fap(15)
      zlm(297)=-zlm(296)*3.6D+01
      zlm(298)= zlm(296)*1.26D+02
      zlm(299)=-zlm(296)*8.4D+01
      zlm(300)= zlm(296)*nine
! L= 9, ML=-8
      zlm(301)= s13585*sr17*fat(14)
      zlm(302)=-zlm(301)*seven*four
      zlm(303)= zlm(301)*seven*ten
      zlm(304)= zlm(302)
      zlm(305)= zlm(301)
! L= 9, ML=-7
      tmp     = s13585*fat(15)
      zlm(306)= tmp*p17
      zlm(307)=-zlm(306)*seven*three
      zlm(308)= zlm(306)*seven*five
      zlm(309)=-zlm(306)*seven
      zlm(310)=-tmp
      zlm(311)= tmp*seven*three
      zlm(312)=-tmp*seven*five
      zlm(313)= tmp*seven
! L= 9, ML=-6
      tmp     = sr247*sr11*sr15*fap(11)
      zlm(314)= tmp*p17
      zlm(315)=-zlm(314)*p15
      zlm(316)=-zlm(315)
      zlm(317)=-zlm(314)
      zlm(318)=-tmp*three
      zlm(319)=-zlm(318)*p15
      zlm(320)=-zlm(319)
      zlm(321)=-zlm(318)
! L= 9, ML=-5
      tmp     = sr247*sr11*fat(13)
      zlm(322)= tmp*five*p17
      zlm(323)=-zlm(322)*ten
      zlm(324)= zlm(322)*five
      zlm(325)=-tmp*five*six
      zlm(326)=-zlm(325)*ten
      zlm(327)= zlm(325)*five
      zlm(328)= tmp
      zlm(329)=-tmp*ten
      zlm(330)= tmp*five
! L= 9, ML=-4
      tmp     = sr5005*sr19*fat(12)
      zlm(331)= tmp*p17
      zlm(332)= zlm(331)
      zlm(333)=-zlm(331)*six
      zlm(334)=-tmp*ten
      zlm(335)= zlm(334)
      zlm(336)=-zlm(334)*six
      zlm(337)= tmp
      zlm(338)= tmp
      zlm(339)=-tmp*six
! L= 9, ML=-3
      tmp     = sr15*sr77*sr19*fap(13)
      zlm(340)= tmp*2.21D+02
      zlm(341)=-zlm(340)*three
      zlm(342)=-tmp*p13*three*five
      zlm(343)=-zlm(342)*three
      zlm(344)= tmp*p13*three
      zlm(345)=-zlm(344)*three
      zlm(346)=-tmp
      zlm(347)= tmp*three
! L= 9, ML=-2
      tmp     = sr55*sr19*fat(11)
      zlm(348)= tmp*2.21D+02
      zlm(349)=-zlm(348)
      zlm(350)=-tmp*seven*p13*three
      zlm(351)=-zlm(350)
      zlm(352)= tmp*seven*p13
      zlm(353)=-zlm(352)
      zlm(354)=-tmp*seven
      zlm(355)=-zlm(354)
! L= 9, ML=-1
      tmp     = sr5*sr19*fat(14)
      zlm(356)= tmp*2.431D+03
      zlm(357)=-tmp*4.004D+03
      zlm(358)= tmp*2.002D+03
      zlm(359)=-tmp*3.08D+02
      zlm(360)= tmp*seven
! L= 9, ML= 0
      tmp     = sr19*fap(14)
      zlm(361)= tmp*1.2155D+04
      zlm(362)=-tmp*2.5740D+04
      zlm(363)= tmp*1.8018D+04
      zlm(364)=-tmp*4.620D+03
      zlm(365)= tmp*3.15D+02
! L= 9, ML= 1
      zlm(366)= zlm(356)
      zlm(367)= zlm(357)
      zlm(368)= zlm(358)
      zlm(369)= zlm(359)
      zlm(370)= zlm(360)
! L= 9, ML= 2
      tmp     = sr55*sr19*fat(9)
      zlm(371)= tmp*2.21D+02
      zlm(372)=-tmp*2.73D+02
      zlm(373)= tmp*9.1D+01
      zlm(374)=-tmp*seven
! L= 9, ML= 3
      zlm(375)=-zlm(341)
      zlm(376)=-zlm(340)
      zlm(377)=-zlm(343)
      zlm(378)=-zlm(342)
      zlm(379)=-zlm(345)
      zlm(380)=-zlm(344)
      zlm(381)=-zlm(347)
      zlm(382)=-zlm(346)
! L= 9, ML= 4
      tmp     = sr5005*sr19*fat(8)
      zlm(383)= tmp*p17
      zlm(384)=-zlm(383)
      zlm(385)=-tmp*ten
      zlm(386)=-zlm(385)
      zlm(387)= tmp
      zlm(388)=-zlm(387)
! L= 9, ML= 5
      zlm(389)= zlm(324)
      zlm(390)= zlm(323)
      zlm(391)= zlm(322)
      zlm(392)= zlm(327)
      zlm(393)= zlm(326)
      zlm(394)= zlm(325)
      zlm(395)= zlm(330)
      zlm(396)= zlm(329)
      zlm(397)= zlm(328)
! L= 9, ML= 6
      tmp     = sr247*sr15*sr11*fap(9)
      zlm(398)= tmp*p17*three
      zlm(399)=-tmp*p17*ten
      zlm(400)= zlm(398)
      zlm(401)=-tmp*three*three
      zlm(402)= tmp*three*ten
      zlm(403)= zlm(401)
! L= 9, ML= 7
      zlm(404)=-zlm(309)
      zlm(405)=-zlm(308)
      zlm(406)=-zlm(307)
      zlm(407)=-zlm(306)
      zlm(408)=-zlm(313)
      zlm(409)=-zlm(312)
      zlm(410)=-zlm(311)
      zlm(411)=-zlm(310)
! L= 9, ML= 8
      zlm(412)= zlm(301)*eight
      zlm(413)=-zlm(412)*seven
      zlm(414)=-zlm(413)
      zlm(415)=-zlm(412)
! L= 9, ML= 9
      zlm(416)= zlm(300)
      zlm(417)= zlm(299)
      zlm(418)= zlm(298)
      zlm(419)= zlm(297)
      zlm(420)= zlm(296)
!
! L=10, ML=-10
      zlm(421)= sr323*sr3003*fap(17)
      zlm(422)=-zlm(421)*4.5D+01
      zlm(423)= zlm(421)*2.1D+02
      zlm(424)=-zlm(423)
      zlm(425)=-zlm(422)
      zlm(426)=-zlm(421)
! L=10, ML=-9
      zlm(427)= sr323*sr3003*sr5*fap(15)
      zlm(428)=-zlm(427)*3.6D+01
      zlm(429)= zlm(427)*1.26D+02
      zlm(430)=-zlm(427)*8.4D+01
      zlm(431)= zlm(427)*nine
! L=10, ML=-8
      tmp     = sr5005*sr17*sr3*fap(16)
      zlm(432)= tmp*p19
      zlm(433)=-zlm(432)*seven*four
      zlm(434)= zlm(432)*seven*ten
      zlm(435)= zlm(433)
      zlm(436)= zlm(432)
      zlm(437)=-tmp
      zlm(438)=-zlm(437)*seven*four
      zlm(439)= zlm(437)*seven*ten
      zlm(440)= zlm(438)
      zlm(441)= zlm(437)
! L=10, ML=-7
      tmp     = sr5005*sr17*fat(15)
      zlm(442)= tmp*p19
      zlm(443)=-zlm(442)*seven*three
      zlm(444)= zlm(442)*seven*five
      zlm(445)=-zlm(442)*seven
      zlm(446)=-tmp*three
      zlm(447)=-zlm(446)*seven*three
      zlm(448)= zlm(446)*seven*five
      zlm(449)=-zlm(446)*seven
! L=10, ML=-6
      tmp     = sr5005*fat(17)
      zlm(450)= tmp*3.23D+02
      zlm(451)=-zlm(450)*p15
      zlm(452)=-zlm(451)
      zlm(453)=-zlm(450)
      zlm(454)=-tmp*1.02D+02
      zlm(455)=-zlm(454)*p15
      zlm(456)=-zlm(455)
      zlm(457)=-zlm(454)
      zlm(458)= tmp*three
      zlm(459)=-zlm(458)*p15
      zlm(460)=-zlm(459)
      zlm(461)=-zlm(458)
! L=10, ML=-5
      tmp     = sr1001*fat(13)
      zlm(462)= tmp*3.23D+02
      zlm(463)=-zlm(462)*ten
      zlm(464)= zlm(462)*five
      zlm(465)=-tmp*1.70D+02
      zlm(466)=-zlm(465)*ten
      zlm(467)= zlm(465)*five
      zlm(468)= tmp*p15
      zlm(469)=-zlm(468)*ten
      zlm(470)= zlm(468)*five
! L=10, ML=-4
      tmp     = sr5005*fat(14)
      zlm(471)= tmp*3.23D+02
      zlm(472)= zlm(471)
      zlm(473)=-zlm(471)*six
      zlm(474)=-tmp*2.55D+02
      zlm(475)= zlm(474)
      zlm(476)=-zlm(474)*six
      zlm(477)= tmp*4.5D+01
      zlm(478)= zlm(477)
      zlm(479)=-zlm(477)*six
      zlm(480)=-tmp
      zlm(481)= zlm(480)
      zlm(482)=-zlm(480)*six
! L=10, ML=-3
      tmp     = sr5005*fat(13)
      zlm(483)= tmp*3.23D+02
      zlm(484)=-zlm(483)*three
      zlm(485)=-tmp*3.57D+02
      zlm(486)=-zlm(485)*three
      zlm(487)= tmp*1.05D+02
      zlm(488)=-zlm(487)*three
      zlm(489)=-tmp*seven
      zlm(490)=-zlm(489)*three
! L=10, ML=-2
      tmp     = sr385*fat(16)
      zlm(491)= tmp*4.199D+03
      zlm(492)=-zlm(491)
      zlm(493)=-tmp*6.188D+03
      zlm(494)=-zlm(493)
      zlm(495)= tmp*2.730D+03
      zlm(496)=-zlm(495)
      zlm(497)=-tmp*3.64D+02
      zlm(498)=-zlm(497)
      zlm(499)= tmp*seven
      zlm(500)=-zlm(499)
! L=10, ML=-1
      tmp     = sr385*sr3*fap(14)
      zlm(501)= tmp*4.199D+03
      zlm(502)=-tmp*7.956D+03
      zlm(503)= tmp*4.914D+03
      zlm(504)=-tmp*1.092D+03
      zlm(505)= tmp*6.3D+01
! L=10, ML= 0
      tmp     = sr21*fap(16)
      zlm(506)= tmp*4.6189D+04
      zlm(507)=-tmp*1.09395D+05
      zlm(508)= tmp*9.0090D+04
      zlm(509)=-tmp*3.0030D+04
      zlm(510)= tmp*3.465D+03
      zlm(511)=-tmp*6.3D+01
! L=10, ML= 1
      zlm(512)= zlm(501)
      zlm(513)= zlm(502)
      zlm(514)= zlm(503)
      zlm(515)= zlm(504)
      zlm(516)= zlm(505)
! L=10, ML= 2
      zlm(517)= zlm(491)*two
      zlm(518)= zlm(493)*two
      zlm(519)= zlm(495)*two
      zlm(520)= zlm(497)*two
      zlm(521)= zlm(499)*two
! L=10, ML= 3
      zlm(522)=-zlm(484)
      zlm(523)=-zlm(483)
      zlm(524)=-zlm(486)
      zlm(525)=-zlm(485)
      zlm(526)=-zlm(488)
      zlm(527)=-zlm(487)
      zlm(528)=-zlm(490)
      zlm(529)=-zlm(489)
! L=10, ML= 4
      zlm(530)= zlm(471)*four
      zlm(531)=-zlm(530)
      zlm(532)= zlm(474)*four
      zlm(533)=-zlm(532)
      zlm(534)= zlm(477)*four
      zlm(535)=-zlm(534)
      zlm(536)= zlm(480)*four
      zlm(537)=-zlm(536)
! L=10, ML= 5
      zlm(538)= zlm(462)
      zlm(539)= zlm(463)
      zlm(540)= zlm(464)
      zlm(541)= zlm(465)
      zlm(542)= zlm(466)
      zlm(543)= zlm(467)
      zlm(544)= zlm(468)
      zlm(545)= zlm(469)
      zlm(546)= zlm(470)
! L=10, ML= 6
      zlm(547)= zlm(450)*six
      zlm(548)=-zlm(450)*p20
      zlm(549)= zlm(547)
      zlm(550)= zlm(454)*six
      zlm(551)=-zlm(454)*p20
      zlm(552)= zlm(550)
      zlm(553)= zlm(458)*six
      zlm(554)=-zlm(458)*p20
      zlm(555)= zlm(553)
! L=10, ML= 7
      zlm(556)=-zlm(445)
      zlm(557)=-zlm(444)
      zlm(558)=-zlm(443)
      zlm(559)=-zlm(442)
      zlm(560)=-zlm(449)
      zlm(561)=-zlm(448)
      zlm(562)=-zlm(447)
      zlm(563)=-zlm(446)
! L=10, ML= 8
      zlm(564)= zlm(432)*eight
      zlm(565)=-zlm(432)*5.6D+01
      zlm(566)=-zlm(565)
      zlm(567)=-zlm(564)
      zlm(568)= zlm(437)*eight
      zlm(569)=-zlm(437)*5.6D+01
      zlm(570)=-zlm(569)
      zlm(571)=-zlm(568)
! L=10, ML= 9
      zlm(572)= zlm(431)
      zlm(573)= zlm(430)
      zlm(574)= zlm(429)
      zlm(575)= zlm(428)
      zlm(576)= zlm(427)
! L=10, ML= 10
      zlm(577)= zlm(421)*ten
      zlm(578)=-zlm(421)*1.20D+02
      zlm(579)= zlm(421)*2.52D+02
      zlm(580)= zlm(578)
      zlm(581)= zlm(577)
!
      return
end


!-----------------------------------------
  subroutine ecpxyz(xyzintecp,maxecpdim)
!-----------------------------------------
!
! Set up angular parts for ECP calculations
!
      use iofile, only : iout
      implicit none
      integer,intent(in) :: maxecpdim
      integer :: maxdim, i, j, k
      real(8),parameter :: zero=0.0D+00, one=1.0D+00, two=2.0D+00, three=3.0D+00
      real(8),parameter :: pi4=1.256637061435917D+01
      real(8),intent(out) :: xyzintecp(0:24,0:24,0:24)
      real(8) :: factor(0:72), denom, xx, yy, zz
!
      maxdim= maxecpdim*4
!
      if(maxdim > 24) then
        write(iout,'(" Error! This program supports up to I function in Subroutine ecpxyz!")')
        call iabort
      endif
!
      denom= three
      do i= 0,3*maxdim-2,2
        factor(i)=one/denom
        denom= denom+two
      enddo
      xyzintecp= zero
      xyzintecp(0,0,0)= pi4

      xx=-one
      do i= 0,maxdim,2
        yy=-one
        do j= 0,maxdim,2
          zz=-one
          do k= 0,maxdim,2
            if(i > 0) then
              xyzintecp(k,j,i)= xx*xyzintecp(k,j,i-2)*factor(k+j+i-2)
            elseif(j > 0) then
              xyzintecp(k,j,i)= yy*xyzintecp(k,j-2,i)*factor(k+j-2+i)
            elseif(k > 0) then
              xyzintecp(k,j,i)= zz*xyzintecp(k-2,j,i)*factor(k-2+j+i)
            endif
            zz= zz+two
          enddo
          yy= yy+two
        enddo
        xx= xx+two
      enddo
!
      return
end


!--------------------------------------------------------------------------------
  subroutine ecpangint1(term1ecp,label1ecp,num1ecp,xyzintecp,maxecpdim)
!--------------------------------------------------------------------------------
!
! Generate angular parts of term1 for ECP calculations
!
      use ecp, only : nx, ny, nz, lmf, lmx, lmy, lmz, zlm, nterm1
      use iofile, only : iout
      implicit none
      integer,parameter :: itri(0:6)=(/0,1,3,6,10,15,21/), mbasis(0:7)=(/1,2,5,11,21,36,57,85/)
      integer,intent(in) :: maxecpdim
      integer,intent(out) :: label1ecp(9,nterm1), num1ecp(2,56,56)
      integer :: icount, i, j, ibasis, jbasis, mx1, my1, mz1, mx2, my2, mz2, lmax, lambda, mu
      integer :: lmindex, kx, ky, kz, kxp, kyp, kzp, ksum, lm, ijx, ijy, ijz
      real(8),parameter :: combi(0:27)= &
&     (/1.0D+00, 1.0D+00,1.0D+00,  1.0D+00,2.0D+00,1.0D+00, 1.0D+00,3.0D+00,3.0D+00,1.0D+00, &
&       1.0D+00,4.0D+00,6.0D+00,4.0D+00,1.0D+00, 1.0D+00,5.0D+00,10.0D+00,10.0D+00,5.0D+00, &
&       1.0D+00, 1.0D+00,6.0D+00,15.0D+00,20.0D+00,15.0D+00,6.0D+00,1.0D+00/)
      real(8),parameter :: zero=0.0D+00, tol=1.0D-15
      real(8),intent(in) :: xyzintecp(0:24,0:24,0:24)
      real(8),intent(out) :: term1ecp(nterm1)
      real(8) :: angint1
!
      icount= 0
!
      do i= 0,maxecpdim
        do ibasis= mbasis(i),mbasis(i+1)-1
          mx1= itri(nx(ibasis))
          my1= itri(ny(ibasis))
          mz1= itri(nz(ibasis))
          do j= 0,i
            do jbasis= mbasis(j),mbasis(j+1)-1
              num1ecp(1,jbasis,ibasis)=icount+1
              mx2= itri(nx(jbasis))
              my2= itri(ny(jbasis))
              mz2= itri(nz(jbasis))
              lmax= i+j
              do lambda= 0,lmax
                do mu=-lambda,lambda
                  lmindex= lambda*(lambda+1)-mu+1
                  do kx= 0,nx(ibasis)
                    do ky= 0,ny(ibasis)
                      do kz= 0,nz(ibasis)
                        do kxp= 0,nx(jbasis)
                          do kyp= 0,ny(jbasis)
                            do kzp= 0,nz(jbasis)
                              ksum=kx+ky+kz+kxp+kyp+kzp
                              if((mod(lambda+ksum,2) /= 1).and.(lambda <= ksum)) then
                                angint1= zero
                                do lm= lmf(lmindex),lmf(lmindex+1)-1
                                  ijx= kx+kxp+lmx(lm)
                                  ijy= ky+kyp+lmy(lm)
                                  ijz= kz+kzp+lmz(lm)
                                  angint1= angint1+zlm(lm)*xyzintecp(ijz,ijy,ijx)
                                enddo
                                if(abs(angint1) > tol)then
                                  icount= icount+1
                                  label1ecp(1,icount)= ksum
                                  label1ecp(2,icount)= lambda
                                  label1ecp(3,icount)= mu
                                  label1ecp(4,icount)= kx
                                  label1ecp(5,icount)= ky
                                  label1ecp(6,icount)= kz
                                  label1ecp(7,icount)= kxp
                                  label1ecp(8,icount)= kyp
                                  label1ecp(9,icount)= kzp
                                  term1ecp(icount)= combi(mx1+kx) *combi(my1+ky) *combi(mz1+kz) * &
&                                                   combi(mx2+kxp)*combi(my2+kyp)*combi(mz2+kzp)* &
&                                                   angint1
                                endif
                              endif
                            enddo
                          enddo
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              enddo
              num1ecp(2,jbasis,ibasis)= icount-num1ecp(1,jbasis,ibasis)+1
            enddo
          enddo
        enddo
      enddo
      if(icount > nterm1) then
        write(iout,'(" Error! Icount exceeds",i8," in Subroutine ecpangint1!")')nterm1
        call iabort
      endif
!
      return
end


!-------------------------------------------------------------------------------------------
  subroutine ecpangint2(term2ecp,label2ecp,num2ecp,xyzintecp,maxecpdim,maxpangdim)
!-------------------------------------------------------------------------------------------
!
! Generate angular parts of term2 for ECP calculations
!
      use ecp, only : nx, ny, nz, lmf, lmx, lmy, lmz, zlm, nterm2
      use iofile, only : iout
      implicit none
      integer,parameter :: itri(0:6)=(/0,1,3,6,10,15,21/), mbasis(0:7)=(/1,2,5,11,21,36,57,85/)
      integer,intent(in) :: maxecpdim, maxpangdim
      integer,intent(out) :: label2ecp(6,nterm2), num2ecp(2,*)
      integer :: icount, numbasis, lambda, mu, lmindex1, lmlabel, i, ltotal, ibasis
      integer :: mx, my, mz, ll, mm, lmindex2, kx, ky, kz, ksum, lm1, ijx, ijy, ijz, lm2
      real(8),parameter :: combi(0:27)= &
&     (/1.0D+00, 1.0D+00,1.0D+00,  1.0D+00,2.0D+00,1.0D+00, 1.0D+00,3.0D+00,3.0D+00,1.0D+00, &
&       1.0D+00,4.0D+00,6.0D+00,4.0D+00,1.0D+00, 1.0D+00,5.0D+00,10.0D+00,10.0D+00,5.0D+00, &
&       1.0D+00, 1.0D+00,6.0D+00,15.0D+00,20.0D+00,15.0D+00,6.0D+00,1.0D+00/)
      real(8),parameter :: zero=0.0D+00, tol=1.0D-15
      real(8),intent(in) :: xyzintecp(0:24,0:24,0:24)
      real(8),intent(out) :: term2ecp(nterm2)
      real(8) :: angint2, tmp
!
      icount= 0
      numbasis= mbasis(maxecpdim)
!
      do lambda= 0,maxpangdim-1
        do mu=-lambda,lambda
          lmindex1= lambda*(lambda+1)-mu+1
          lmlabel=(lambda*(lambda+1)-mu)*numbasis
          do i= 0,maxecpdim
            ltotal= lambda+i
            do ibasis= mbasis(i),mbasis(i+1)-1
              mx= itri(nx(ibasis))
              my= itri(ny(ibasis))
              mz= itri(nz(ibasis))
              num2ecp(1,lmlabel+ibasis)= icount+1
              do ll= 0,ltotal
                do mm=-ll,ll
                  lmindex2= ll*(ll+1)-mm+1
                  do kx= 0,nx(ibasis)
                    do ky= 0,ny(ibasis)
                      do kz= 0,nz(ibasis)
                      ksum= kx+ky+kz
                      angint2= zero
                      do lm1= lmf(lmindex1),lmf(lmindex1+1)-1
                        ijx= kx+lmx(lm1)
                        ijy= ky+lmy(lm1)
                        ijz= kz+lmz(lm1)
                        tmp= zero
                        do lm2= lmf(lmindex2),lmf(lmindex2+1)-1
                          tmp= tmp+zlm(lm2)*xyzintecp(ijz+lmz(lm2),ijy+lmy(lm2),ijx+lmx(lm2))
                        enddo
                        angint2= angint2+zlm(lm1)*tmp
                      enddo
                      if(abs(angint2) > tol) then
                        icount= icount+1
                        label2ecp(1,icount)= ll
                        label2ecp(2,icount)= ksum
                        label2ecp(3,icount)= mm
                        label2ecp(4,icount)= kx
                        label2ecp(5,icount)= ky
                        label2ecp(6,icount)= kz
                        term2ecp(icount)= combi(mx+kx)*combi(my+ky)*combi(mz+kz)*angint2
                      endif
                      enddo
                    enddo
                  enddo
                enddo
              enddo
              num2ecp(2,lmlabel+ibasis)= icount-num2ecp(1,lmlabel+ibasis)+1
            enddo
          enddo
        enddo
      enddo
      if(icount > nterm2) then
        write(iout,'(" Error! Icount exceeds",i8," in Subroutine ecpangint2!")')nterm2
        call iabort
      endif
!
      return
end


!-----------------------------------------------------------------
  subroutine ecpangint0(term0ecp,xyzintecp,maxecpdim,maxpangdim)
!-----------------------------------------------------------------
!
! Generate angular parts of term1 for ECP calculations with same center
!
      use ecp, only : nx, ny, nz, lmf, lmx, lmy, lmz, zlm
      implicit none
      integer,parameter :: mbasis(0:7)=(/1,2,5,11,21,36,57,85/)
      integer,intent(in) :: maxecpdim, maxpangdim
      integer :: numbasis, lambda, mu, lmindex, lmlabel, i, ibasis, lm, ijx, ijy, ijz
      real(8),parameter :: zero=0.0D+00
      real(8),intent(in) :: xyzintecp(0:24,0:24,0:24)
      real(8),intent(out) :: term0ecp(*)
      real(8) :: angint0
!
      numbasis= mbasis(maxecpdim)
!
      do lambda= 0,maxpangdim-1
        do mu=-lambda,lambda
          lmindex= lambda*(lambda+1)-mu+1
          lmlabel=(lambda*(lambda+1)-mu)*numbasis
          do i= 0,maxecpdim
            do ibasis= mbasis(i),mbasis(i+1)-1
              angint0= zero
              do lm= lmf(lmindex),lmf(lmindex+1)-1
                ijx= nx(ibasis)+lmx(lm)
                ijy= ny(ibasis)+lmy(lm)
                ijz= nz(ibasis)+lmz(lm)
                angint0= angint0+zlm(lm)*xyzintecp(ijz,ijy,ijx)
              enddo
              term0ecp(lmlabel+ibasis)= angint0
            enddo
          enddo
        enddo
      enddo
!
      return
end


!-----------------------------------------------------------------
  function dawson(val)
!-----------------------------------------------------------------
!
! Evaluate the Dawson function
!
      implicit none
      integer :: igrid, last, i
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, ten=1.0D+01
      real(8),parameter :: griddaw=5.0D+00
      real(8),intent(in) :: val
      real(8) :: dawson, val1, val2, val3, val4, val5, val6, val7, val8, val9
      real(8) :: val2inv, factor, total
      real(8) :: dawgrid(0:9,0:49)
      data dawgrid / &
&    6.051722509208561D-16,   9.999999999994005D-01,   9.777221521513968D-11, &
&   -6.666666728319638D-01,   1.965028896818846D-07,   2.666630955078337D-01, &
&    3.902825044190826D-05,  -7.644842301092071D-02,   9.777389039626225D-04, &
&    1.519711057206845D-02,   2.953748577202614D-08,   9.999989604795872D-01, &
&    1.632904172907173D-05,  -6.668173339622047D-01,   9.033269735325170D-04, &
&    2.629969013292320D-01,   1.018881107249708D-02,  -9.510163851001064D-02, &
&    2.186581860556165D-02,   4.206328735375168D-03,   5.022734342717940D-06, &
&    9.998964691434666D-01,   9.591681774814235D-04,  -6.719264553715743D-01, &
&    1.890290652354323D-02,   2.201777674225539D-01,   7.906931394260616D-02, &
&   -1.674426202438445D-01,   6.691761616978956D-02,  -8.476321032168930D-03, &
&    3.669932078830077D-05,   9.994086664526127D-01,   4.313802985030914D-03, &
&   -6.854492931066708D-01,   5.411705735859217D-02,   1.587484871726615D-01, &
&    1.508505021326560D-01,  -2.216150267576290D-01,   9.087352189345086D-02, &
&   -1.320491480520023D-02,  -9.240529554626878D-04,   1.009683221114372D+00, &
&   -4.461635662652316D-02,  -5.492383023037898D-01,  -1.901858785183352D-01, &
&    4.515560484687740D-01,  -8.369896234661265D-02,  -1.005130134305425D-01, &
&    5.429770655133835D-02,  -8.280903040450580D-03,  -1.029775943821545D-02, &
&    1.092399269442080D+00,  -3.697768628624183D-01,   1.981878620126198D-01, &
&   -1.297394353718835D+00,   1.547787220212697D+00,  -8.091518439150200D-01, &
&    2.089255992186689D-01,  -2.290165554070681D-02,   3.020315018984405D-04, &
&   -3.635700884765439D-02,   1.287744381657310D+00,  -1.021778458451531D+00, &
&    1.469941481761009D+00,  -2.894980361086176D+00,   2.888159254985652D+00, &
&   -1.560225529053189D+00,   4.799668635758824D-01,  -8.005976523151092D-02, &
&    5.668640259844556D-03,  -2.430672682332134D-02,   1.218214029581446D+00, &
&   -8.454866418703851D-01,   1.212966595174907D+00,  -2.658748167631678D+00, &
&    2.747139007450606D+00,  -1.506200145273536D+00,   4.674338456679980D-01, &
&   -7.853618452992746D-02,   5.604485954876450D-03,   2.086326626721770D-01, &
&   -7.659095371377669D-02,   2.356228131334686D+00,  -3.409613779527687D+00, &
&    1.635769789964874D+00,   8.474024580105899D-02,  -4.047464050401199D-01, &
&    1.742063022481715D-01,  -3.295397297470270D-02,   2.452076490858092D-03, &
&    7.819875729224917D-01,  -2.937321717501164D+00,   8.705241547587803D+00, &
&   -1.163605938718322D+01,   8.493710605657423D+00,  -3.729830062303937D+00, &
&    1.010943062539432D+00,  -1.638322120981495D-01,   1.416996484227287D-02, &
&   -4.700124848768664D-04,   1.304119931791932D+00,  -5.308704826333093D+00, &
&    1.349480033480972D+01,  -1.728228980794192D+01,   1.277511030971419D+01, &
&   -5.895380025040220D+00,   1.741581652446646D+00,  -3.223913333975149D-01, &
&    3.425315766198099D-02,  -1.601172584406407D-03,   8.453967279834436D-01, &
&   -3.473857303563860D+00,   1.023282373069887D+01,  -1.389938194430536D+01, &
&    1.051971338675129D+01,  -4.892911707033860D+00,   1.444531320532356D+00, &
&   -2.658060541288357D-01,   2.796549043503804D-02,  -1.290655862378755D-03, &
&   -1.168868023013557D+00,   4.051680384856491D+00,  -2.268990395019065D+00, &
&   -1.778797427179115D+00,   2.962044850454449D+00,  -1.749791207253024D+00, &
&    5.726770875040680D-01,  -1.102658329522576D-01,   1.177126126343049D-02, &
&   -5.409385300534173D-04,  -3.983678090018918D+00,   1.380643270044792D+01, &
&   -1.729954790908351D+01,   1.173647405757583D+01,  -4.853561872641671D+00, &
&    1.264482812164978D+00,  -2.026489282503666D-01,   1.798868732943707D-02, &
&   -6.095550168767980D-04,  -9.547785170957528D-06,  -5.814425354970969D+00, &
&    1.973235274865224D+01,  -2.582728197811491D+01,   1.889725061452376D+01, &
&   -8.720183404777921D+00,   2.656809882256313D+00,  -5.369883591179374D-01, &
&    6.961546456260244D-02,  -5.261176970280754D-03,   1.767788777026804D-04, &
&   -5.329689633295526D+00,   1.831861036236389D+01,  -2.399576707118949D+01, &
&    1.751397831213208D+01,  -8.049000441611372D+00,   2.439847689603503D+00, &
&   -4.902670710796377D-01,   6.315275734302311D-02,  -4.740151937244599D-03, &
&    1.581271836076713D-04,  -2.742874229129735D+00,   1.105990371185430D+01, &
&   -1.494089841574973D+01,   1.092324271927255D+01,  -4.964303598486863D+00, &
&    1.477099885778468D+00,  -2.898963332247264D-01,   3.633739182918682D-02, &
&   -2.646216460201094D-03,   8.543732483697165D-05,   5.334358295001977D-01, &
&    2.380065505895000D+00,  -4.718367967869862D+00,   3.898594616817992D+00, &
&   -1.860403980125197D+00,   5.625617963018695D-01,  -1.102134285396532D-01, &
&    1.363727777246505D-02,  -9.729401144281072D-04,   3.060634340168288D-05, &
&    3.080620069423487D+00,  -4.005453782311422D+00,   2.397694500496021D+00, &
&   -7.283035437468496D-01,   7.398635616053244D-02,   2.330434398922989D-02, &
&   -9.972367942442934D-03,   1.656144541521961D-03,  -1.374241748923812D-04, &
&    4.705353090582574D-06,   4.300353044663805D+00,  -6.909569738185462D+00, &
&    5.471347664874574D+00,  -2.626260134921639D+00,   8.275246738550115D-01, &
&   -1.761785998740552D-01,   2.523933144701089D-02,  -2.340119146346428D-03, &
&    1.271874623248163D-04,  -3.083104605409867D-06,   4.396530707222811D+00, &
&   -7.134196626551821D+00,   5.704257216451345D+00,  -2.766989300405958D+00, &
&    8.821357230714314D-01,  -1.902941758333474D-01,   2.766965259236612D-02, &
&   -2.608903835934439D-03,   1.445153102345324D-04,  -3.579244353753892D-06, &
&    3.896858160212054D+00,  -6.065851873388711D+00,   4.688902043479053D+00, &
&   -2.203994088320108D+00,   6.814255243986297D-01,  -1.425844799207229D-01, &
&    2.010800406935943D-02,  -1.838346719741425D-03,   9.870384880496008D-05, &
&   -2.368575908529228D-06,   3.246922643039968D+00,  -4.735897812566311D+00, &
&    3.479189852464805D+00,  -1.562037274177388D+00,   4.623945424588869D-01, &
&   -9.275624721235535D-02,   1.254985334572399D-02,  -1.101238261520193D-03, &
&    5.676415575057387D-05,  -1.307863413499134D-06,   2.675504458023278D+00, &
&   -3.616587250649823D+00,   2.504601603004815D+00,  -1.066969118221491D+00, &
&    3.007058498254446D-01,  -5.754679238803369D-02,   7.437700565744401D-03, &
&   -6.240180589474892D-04,   3.077414469741080D-05,  -6.786944578519707D-07, &
&    2.243099420626957D+00,  -2.804674807408287D+00,   1.826964650222504D+00, &
&   -7.370149564337481D-01,   1.974115492186115D-01,  -3.598626016422626D-02, &
&    4.437136966000534D-03,  -3.555380013897935D-04,   1.675926942372285D-05, &
&   -3.535061057294843D-07,   1.932029874176854D+00,  -2.243964191783463D+00, &
&    1.377717335005858D+00,  -5.270258066140214D-01,   1.343055779563454D-01, &
&   -2.334179536965425D-02,   2.747914965941448D-03,  -2.104488646833919D-04, &
&    9.489056293487565D-06,  -1.915773367624041D-07,   1.706708849350314D+00, &
&   -1.853496313508908D+00,   1.076949899974014D+00,  -3.918691452174971D-01, &
&    9.525721411066638D-02,  -1.582000633173914D-02,   1.781881255989836D-03, &
&   -1.306821199166410D-04,   5.646571822689784D-06,  -1.093031483485892D-07, &
&    1.537669196618148D+00,  -1.571461484279969D+00,   8.677913416794271D-01, &
&   -3.013779063895062D-01,   7.008664461576448D-02,  -1.115202231884632D-02, &
&    1.204695757535073D-03,  -8.479858763893169D-05,   3.518648920327665D-06, &
&   -6.543865306252688D-08,   1.405549318697144D+00,  -1.358932771279227D+00, &
&    7.158335591164556D-01,  -2.379933527355135D-01,   5.308868832732922D-02, &
&   -8.112838465920787D-03,   8.423992317003719D-04,  -5.703190720497197D-05, &
&    2.277172150589343D-06,  -4.076640715826195D-08,   1.298576035768416D+00, &
&   -1.192809905599154D+00,   6.011669931292463D-01,  -1.918194294121542D-01, &
&    4.113486133628721D-02,  -6.049542193356675D-03,   6.049555543509496D-04, &
&   -3.946441816132111D-05,   1.518924456429258D-06,  -2.621967464720331D-08, &
&    1.209577830183879D+00,  -1.059221212897871D+00,   5.120400290656999D-01, &
&   -1.571298707954724D-01,   3.245452949501371D-02,  -4.601384001133819D-03, &
&    4.438769850238281D-04,  -2.794561679474710D-05,   1.038387435358890D-06, &
&   -1.730927703431548D-08,   1.133997350075610D+00,  -9.494408737898824D-01, &
&    4.411658021914067D-01,  -1.304366957767376D-01,   2.599117585807237D-02, &
&   -3.557972489308219D-03,   3.315734790826391D-04,  -2.017462082284952D-05, &
&    7.246921372578375D-07,  -1.168083844809918D-08,   1.068780177926012D+00, &
&   -8.576790733538303D-01,   3.837794250426504D-01,  -1.095003008494338D-01, &
&    2.108052719188805D-02,  -2.790056105619201D-03,   2.515113598636368D-04, &
&   -1.480820888493274D-05,   5.148527578287036D-07,  -8.033830684925624D-09, &
&    1.011779626492335D+00,  -7.799125024116445D-01,   3.366217910576227D-01, &
&   -9.281801625720301D-02,   1.728649366216113D-02,  -2.214770949664605D-03, &
&    1.933544248618440D-04,  -1.102846916299806D-05,   3.715460199775798D-07, &
&   -5.618835277160696D-09,   9.614096479074705D-01,  -7.132165323271270D-01, &
&    2.973689581221718D-01,  -7.934126227493824D-02,   1.431181444060736D-02, &
&   -1.777017303985521D-03,   1.504053144171016D-04,  -8.319406490500501D-06, &
&    2.718617342711667D-07,  -3.988496802653806D-09,   9.165526309123725D-01, &
&   -6.555196224552808D-01,   2.643839356220485D-01,  -6.834051008923832D-02, &
&    1.195314538712647D-02,  -1.439850416169431D-03,   1.182719583461593D-04, &
&   -6.350585515470334D-06,   2.014901001305661D-07,  -2.870524681794359D-09, &
&    8.762449238894968D-01,  -6.051156577245481D-01,   2.363694886681162D-01, &
&   -5.925729795070593D-02,   1.005977937101519D-02,  -1.176725686645202D-03, &
&    9.389268779118348D-05,  -4.898415677199682D-06,   1.510293060240204D-07, &
&   -2.091176767916272D-09,   8.397973374736750D-01,  -5.607717546712178D-01, &
&    2.123900746196760D-01,  -5.169276657599643D-02,   8.525650051352768D-03, &
&   -9.692951718191823D-04,   7.519387477589534D-05,  -3.814758613123273D-06, &
&    1.143932740559919D-07,  -1.540667804508020D-09,   8.066453081016311D-01, &
&   -5.214997953854644D-01,   1.917127957710592D-01,  -4.534177013409804D-02, &
&    7.271568748050543D-03,  -8.041981986548870D-04,   6.070341686089757D-05, &
&   -2.997121855987458D-06,   8.747939440781839D-08,  -1.146908820242620D-09, &
&    7.763328700748650D-01,  -4.865131195188174D-01,   1.737645584369129D-01, &
&   -3.997049326403647D-02,   6.238170884938409D-03,  -6.716462489699137D-04, &
&    4.936811186568123D-05,  -2.373940361615341D-06,   6.749309265777579D-08, &
&   -8.620120298597762D-10,   7.484872168611639D-01,  -4.551777208971689D-01, &
&    1.580916582586361D-01,  -3.539751489571131D-02,   5.380378405930295D-03, &
&   -5.643726412581538D-04,   4.042412621397579D-05,  -1.894535342555393D-06, &
&    5.250281623647732D-08,  -6.536814794056658D-10,   7.228002002617220D-01, &
&   -4.269770225009569D-01,   1.443309279425442D-01,  -3.148047647059620D-02, &
&    4.663564236690544D-03,  -4.769181755117127D-04,   3.331061449736308D-05, &
&   -1.522556290172317D-06,   4.115564531122071D-08,  -4.998334049864280D-10, &
&    6.990146269097030D-01,  -4.014859925399112D-01,   1.321887743925377D-01, &
&   -2.810653190187287D-02,   4.060849619196801D-03,  -4.051368747339830D-04, &
&    2.761111481851752D-05,  -1.231622827617666D-06,   3.249232324267357D-08, &
&   -3.851739792444945D-10,   6.769139442201537D-01,  -3.783518223820482D-01, &
&    1.214257096016929D-01,  -2.518540056502214D-02,   3.551170329362119D-03, &
&   -3.458487212194843D-04,   2.301316857550481D-05,  -1.002382719404352D-06, &
&    2.582501734241042D-08,  -2.989865561703131D-10,   6.563143660892078D-01, &
&   -3.572792987091244D-01,   1.118447759141316D-01,  -2.264424363481021D-02, &
&    3.117873830445191D-03,  -2.965922633845162D-04,   1.928010426958492D-05, &
&   -8.204973417272264D-07,   2.065535195607614D-08,  -2.336795023974850D-10, &
&    6.370587839081798D-01,  -3.380195885781189D-01,   1.032827760925571D-01, &
&   -2.042383776876847D-02,   2.747688779923439D-03,  -2.554461321436570D-04, &
&    1.623106686867063D-05,  -6.752440992628909D-07,   1.661870682317888D-08, &
&   -1.838200346615760D-10,   6.190120033975020D-01,  -3.203615420730855D-01, &
&    9.560355253336551D-02,  -1.847568500303102D-02,   2.429959362224784D-03, &
&   -2.208987753263688D-04,   1.372671585753052D-05,  -5.585351186952662D-07, &
&    1.344589462985757D-08,  -1.454832089595699D-10,   6.020569792050626D-01, &
&   -3.041248766244393D-01,   8.869278306835553D-02,  -1.675980697559348D-02, &
&    2.156070393944811D-03,  -1.917523666799123D-04,   1.165887107556909D-05, &
&   -4.642205187486950D-07,   1.093649208022802D-08,  -1.158080645769749D-10, &
&    5.860918095088875D-01,  -2.891547845930740D-01,   8.245391253436112D-02, &
&   -1.524304526077679D-02,   1.919011474773674D-03,  -1.670512229576893D-04, &
&    9.942940784081151D-06,  -3.875885889643397D-07,   8.940082061173694D-09, &
&   -9.269166157280035D-11,   5.710273159428749D-01,  -2.753176292173221D-01, &
&    7.680494417517188D-02,  -1.389774006501874D-02,   1.713043581841007D-03, &
&   -1.460279888984735D-04,   8.512332244579305D-06,  -3.250037532305761D-07, &
&    7.342928218127729D-09,  -7.457599172953888D-11/
!
      val1= abs(val)
      if(val1 < ten) then
        igrid= aint(val1*griddaw)
        val2= val1*val1
        val3= val1*val1*val1
        val4= val2*val2
        val5= val2*val3
        val6= val3*val3
        val7= val3*val4
        val8= val4*val4
        val9= val4*val4*val1
        dawson= dawgrid(0,igrid)     +dawgrid(1,igrid)*val1+dawgrid(2,igrid)*val2 &
&              +dawgrid(3,igrid)*val3+dawgrid(4,igrid)*val4+dawgrid(5,igrid)*val5 &
&              +dawgrid(6,igrid)*val6+dawgrid(7,igrid)*val7+dawgrid(8,igrid)*val8 &
&              +dawgrid(9,igrid)*val9
        if(val < zero) dawson=-dawson
      else
        val2inv= one/(val1*val1)
        last= 9
        factor=(last+half)*val2inv
        total= factor
        do i= last,1,-1
          factor= factor-val2inv
          total = factor*(one+total)
        enddo
        dawson= val*factor*(one+total)
      endif
!
      return
end


!-----------------------------------------------------------------
  function dawerf(val)
!-----------------------------------------------------------------
!
! Evaluate the Dawson-Error function
!
      implicit none
      integer :: igrid, last, i
      real(8),parameter :: zero=0.0D+00, half=0.5D+00, one=1.0D+00, ten=1.0D+01
      real(8),parameter :: griddaw=5.0D+00
      real(8),intent(in) :: val
      real(8) :: dawerf, val1, val2, val3, val4, val5, val6, val7, val8, val9
      real(8) :: val2inv, factor, total
      real(8) :: dawerfgrid(0:9,0:49)
      data dawerfgrid / &
&   -1.464050227112128D-15,   1.467214183101818D-12,   5.641895833050387D-01, &
&    1.558305800769015D-08,  -3.761268974805824D-01,   9.539122392875117D-06, &
&    1.440726531464679D-01,   7.747133582156654D-04,  -4.274238191602858D-02, &
&    8.037433984828787D-03,  -9.951223039681468D-09,   3.665620582051934D-07, &
&    5.641835068621823D-01,   5.980782869308470D-05,  -3.765143811500408D-01, &
&    1.739802654227826D-03,   1.386911329865945D-01,   1.214475306151317D-02, &
&   -5.768347019627200D-02,   1.737017358803486D-02,   3.505730303251949D-06, &
&   -6.903611739405604D-05,   5.647934766198526D-01,  -3.073493000942437D-03, &
&   -3.661454290548012D-01,  -2.118740661715355D-02,   1.725643524039445D-01, &
&   -2.009535893707317D-02,  -3.975174841729590D-02,   1.293250194958627D-02, &
&    1.753521010416664D-04,  -2.509136076865635D-03,   5.802654485366705D-01, &
&   -6.059720330548457D-02,  -2.278842721571183D-01,  -2.440739289353691D-01, &
&    4.136483694069031D-01,  -1.888713464726069D-01,   2.966074254247206D-02, &
&    1.521467396299521D-04,   1.493399476360794D-03,  -1.710408420590437D-02, &
&    6.523669837322738D-01,  -2.691935415762747D-01,   1.616337380445307D-01, &
&   -7.309645117937897D-01,   8.210629072070301D-01,  -4.089382992885937D-01, &
&    9.928847617461442D-02,  -9.678922904382312D-03,   7.022349491783161D-04, &
&   -1.140029281594338D-02,   6.351255670624923D-01,  -2.418066058845062D-01, &
&    1.396558590557061D-01,  -7.280085092378150D-01,   8.317740595851261D-01, &
&   -4.189875032602776D-01,   1.031772638668831D-01,  -1.026592845497949D-02, &
&   -3.507828429386371D-02,   2.514441281612858D-01,  -2.243057271902791D-01, &
&    1.399931801042554D+00,  -1.879615077893551D+00,   9.304166624505790D-01, &
&   -7.777252579099405D-02,  -9.776713912647873D-02,   3.688751816961262D-02, &
&   -4.175226370687775D-03,  -1.635192995852709D-01,   1.072530607692978D+00, &
&   -2.560291425726306D+00,   5.281908470059739D+00,  -6.032401344817990D+00, &
&    3.896121221224385D+00,  -1.491676614240306D+00,   3.361654360188453D-01, &
&   -4.090475833958504D-02,   2.031481573663880D-03,  -2.969114851113720D-01, &
&    1.833603926731921D+00,  -4.491789807629176D+00,   8.143636663140860D+00, &
&   -8.760254633771366D+00,   5.630966936454882D+00,  -2.227785302457666D+00, &
&    5.371041861773432D-01,  -7.292479631053337D-02,   4.300877192829566D-03, &
&   -1.170572792872378D-03,   3.848954456804926D-01,  -1.336433123369666D+00, &
&    4.132965957640049D+00,  -5.481656626887961D+00,   3.843409548629058D+00, &
&   -1.577751169995361D+00,   3.850753397958487D-01,  -5.217397376886755D-02, &
&    3.041464194679571D-03,   1.174905657055873D+00,  -4.883509768159756D+00, &
&    9.159562518680080D+00,  -8.072949374706283D+00,   3.649353826728024D+00, &
&   -7.134423639453488D-01,  -6.066654191095896D-02,   6.016752490052443D-02, &
&   -1.155594302949185D-02,   7.831369932113179D-04,   2.807721647232346D+00, &
&   -1.157888576964923D+01,   2.136828102040118D+01,  -2.106630531361268D+01, &
&    1.254393798005916D+01,  -4.774847360791843D+00,   1.176340252893272D+00, &
&   -1.821681364947661D-01,   1.615253878701492D-02,  -6.256883756474562D-04, &
&    3.348333345889982D+00,  -1.366040981604734D+01,   2.492960224723117D+01, &
&   -2.462000868400332D+01,   1.482320891370940D+01,  -5.749289456829261D+00, &
&    1.454034752069000D+00,  -2.330356890766367D-01,   2.158734562461293D-02, &
&   -8.837375538626612D-04,   1.419529100885629D+00,  -7.035253702007297D+00, &
&    1.481226734600577D+01,  -1.560431657184057D+01,   9.656736117340066D+00, &
&   -3.774835259282677D+00,   9.508132419060235D-01,  -1.505577918487977D-01, &
&    1.369900313160904D-02,  -5.483047252040330D-04,  -2.544030075609960D+00, &
&    5.695392362927285D+00,  -3.367263672606511D+00,  -4.553991907628051D-01, &
&    1.538810907283111D+00,  -8.737041337394308D-01,   2.593827681045920D-01, &
&   -4.458505602238055D-02,   4.221228573272723D-03,  -1.714398979549848D-04, &
&   -6.348319107215338D+00,   1.714154675262704D+01,  -1.867786552351831D+01, &
&    1.149466787193640D+01,  -4.458975430819412D+00,   1.133760418179475D+00, &
&   -1.886815669411439D-01,   1.972420037565071D-02,  -1.164553455654676D-03, &
&    2.908398028937653D-05,  -7.854162583305704D+00,   2.142212197657538D+01, &
&   -2.408684341726629D+01,   1.548232235238074D+01,  -6.349167351893407D+00, &
&    1.731172160849376D+00,  -3.145799647601360D-01,   3.678306401855272D-02, &
&   -2.513096665959857D-03,   7.647153503181144D-05,  -6.587158759431202D+00, &
&    1.809711106668978D+01,  -2.020811667944246D+01,   1.284255197793339D+01, &
&   -5.194060318853190D+00,   1.394154820937142D+00,  -2.490172759971632D-01, &
&    2.858253201482559D-02,  -1.914670041028719D-03,   5.705978923133756D-05, &
&   -3.695796883540154D+00,   1.087367011914940D+01,  -1.218591014395560D+01, &
&    7.644375518853737D+00,  -3.028284844164388D+00,   7.924612082670810D-01, &
&   -1.375526789103503D-01,   1.530541438488127D-02,  -9.919323992952312D-04, &
&    2.855217532838284D-05,  -7.482906093884771D-01,   3.883066123961691D+00, &
&   -4.815793457346782D+00,   3.110883017076421D+00,  -1.235253709919397D+00, &
&    3.196005001244074D-01,  -5.440141038991109D-02,   5.903846745328878D-03, &
&   -3.717356719999930D-04,   1.036535141718575D-05,   1.294527932794106D+00, &
&   -7.257393636901485D-01,  -1.937132305125310D-01,   4.064505512011646D-01, &
&   -2.178302178553769D-01,   6.438426042200402D-02,  -1.171442818237512D-02, &
&    1.313255790973573D-03,  -8.371099516452187D-05,   2.332327334871007D-06, &
&    2.281922905825749D+00,  -2.850457119940149D+00,   1.838606956259343D+00, &
&   -7.276786094946810D-01,   1.890904180882148D-01,  -3.296386797751393D-02, &
&    3.813581584120174D-03,  -2.792448147132996D-04,   1.157321481202174D-05, &
&   -2.018677537492912D-07,   2.522613418072983D+00,  -3.347167236826221D+00, &
&    2.294223284946532D+00,  -9.714859435212977D-01,   2.729674006790809D-01, &
&   -5.220281515163948D-02,   6.755712070554076D-03,  -5.685059381042721D-04, &
&    2.816402860507764D-05,  -6.248236956132250D-07,   2.392463694043553D+00, &
&   -3.093955395832736D+00,   2.075255779445555D+00,  -8.610196537332251D-01, &
&    2.371386057955971D-01,  -4.445496071742221D-02,   5.638648897267226D-03, &
&   -4.649614499635758D-04,   2.256475450295823D-05,  -4.902398768025890D-07, &
&    2.145937647444013D+00,  -2.631801405379542D+00,   1.690150488548531D+00, &
&   -6.738047236875862D-01,   1.786235824007151D-01,  -3.226072843616258D-02, &
&    3.944304415481221D-03,  -3.136001516209016D-04,   1.467622326425910D-05, &
&   -3.074945512748113D-07,   1.902545939559848D+00,  -2.193377544986805D+00, &
&    1.339117421330826D+00,  -5.098338957524816D-01,   1.293802434109594D-01, &
&   -2.240054608830456D-02,   2.627931671017149D-03,  -2.006114668182249D-04, &
&    9.018321589031551D-06,  -1.815611881036022D-07,   1.698760592612476D+00, &
&   -1.840334358360287D+00,   1.067257847340245D+00,  -3.877038263953083D-01, &
&    9.410586515745797D-02,  -1.560774149743725D-02,   1.755780460356919D-03, &
&   -1.286180118934224D-04,   5.551312083652253D-06,  -1.073484549117828D-07, &
&    1.535758007637357D+00,  -1.568404091098504D+00,   8.656166142754264D-01, &
&   -3.004751789766771D-01,   6.984565529807552D-02,  -1.110911649012001D-02, &
&    1.199601203814210D-03,  -8.440957173131123D-05,   3.501314999817940D-06, &
&   -6.509526124657302D-08,   1.405138044663853D+00,  -1.358296585353194D+00, &
&    7.153960262808818D-01,  -2.378177606181386D-01,   5.304337147871756D-02, &
&   -8.105039014296582D-03,   8.415040453986343D-04,  -5.696583648760865D-05, &
&    2.274326714246108D-06,  -4.071192810443949D-08,   1.298496606911207D+00, &
&   -1.192690990251207D+00,   6.010878440391463D-01,  -1.917886898653514D-01, &
&    4.112718443148547D-02,  -6.048263688424006D-03,   6.048135693374317D-04, &
&   -3.945427884611453D-05,   1.518501982745022D-06,  -2.621184913044049D-08, &
&    1.209564029814801D+00,  -1.059201198185802D+00,   5.120271247188780D-01, &
&   -1.571250162707902D-01,   3.245335520362024D-02,  -4.601194585987844D-03, &
&    4.438566117780785D-04,  -2.794420778070469D-05,   1.038330578930851D-06, &
&   -1.730825714974198D-08,   1.133995188487756D+00,  -9.494378342324285D-01, &
&    4.411639021663465D-01,  -1.304360028035252D-01,   2.599101334929775D-02, &
&   -3.557947077695818D-03,   3.315708294828669D-04,  -2.017444318953867D-05, &
&    7.246851891792878D-07,  -1.168071763897105D-08,   1.068779872142251D+00, &
&   -8.576786561036843D-01,   3.837791719517366D-01,  -1.095002112816680D-01, &
&    2.108050681125562D-02,  -2.790053013426468D-03,   2.515110470422076D-04, &
&   -1.480818853741970D-05,   5.148519856610921D-07,  -8.033817659372651D-09, &
&    1.011779587362307D+00,  -7.799124505569585D-01,   3.366217605117409D-01, &
&   -9.281800575932100D-02,   1.728649134246527D-02,  -2.214770607895047D-03, &
&    1.933543912875351D-04,  -1.102846704240979D-05,   3.715459418360281D-07, &
&   -5.618833997239144D-09,   9.613901115608500D-01,  -7.131910271662479D-01, &
&    2.973541593751507D-01,  -7.933625347946173D-02,   1.431072463057332D-02, &
&   -1.776859225880105D-03,   1.503900283396192D-04,  -8.318456258752728D-06, &
&    2.718272775078469D-07,  -3.987941497667665D-09,   9.165526309123725D-01, &
&   -6.555196224552808D-01,   2.643839356220485D-01,  -6.834051008923832D-02, &
&    1.195314538712647D-02,  -1.439850416169431D-03,   1.182719583461593D-04, &
&   -6.350585515470334D-06,   2.014901001305661D-07,  -2.870524681794359D-09, &
&    8.762449238894968D-01,  -6.051156577245481D-01,   2.363694886681162D-01, &
&   -5.925729795070593D-02,   1.005977937101519D-02,  -1.176725686645202D-03, &
&    9.389268779118348D-05,  -4.898415677199682D-06,   1.510293060240204D-07, &
&   -2.091176767916272D-09,   8.397973374736750D-01,  -5.607717546712178D-01, &
&    2.123900746196760D-01,  -5.169276657599643D-02,   8.525650051352768D-03, &
&   -9.692951718191823D-04,   7.519387477589534D-05,  -3.814758613123273D-06, &
&    1.143932740559919D-07,  -1.540667804508020D-09,   8.066453081016311D-01, &
&   -5.214997953854644D-01,   1.917127957710592D-01,  -4.534177013409804D-02, &
&    7.271568748050543D-03,  -8.041981986548870D-04,   6.070341686089757D-05, &
&   -2.997121855987458D-06,   8.747939440781839D-08,  -1.146908820242620D-09, &
&    7.763328700748650D-01,  -4.865131195188174D-01,   1.737645584369129D-01, &
&   -3.997049326403647D-02,   6.238170884938409D-03,  -6.716462489699137D-04, &
&    4.936811186568123D-05,  -2.373940361615341D-06,   6.749309265777579D-08, &
&   -8.620120298597762D-10,   7.484872168611639D-01,  -4.551777208971689D-01, &
&    1.580916582586361D-01,  -3.539751489571131D-02,   5.380378405930295D-03, &
&   -5.643726412581538D-04,   4.042412621397579D-05,  -1.894535342555393D-06, &
&    5.250281623647732D-08,  -6.536814794056658D-10,   7.228002002617220D-01, &
&   -4.269770225009569D-01,   1.443309279425442D-01,  -3.148047647059620D-02, &
&    4.663564236690544D-03,  -4.769181755117127D-04,   3.331061449736308D-05, &
&   -1.522556290172317D-06,   4.115564531122071D-08,  -4.998334049864280D-10, &
&    6.990146269097030D-01,  -4.014859925399112D-01,   1.321887743925377D-01, &
&   -2.810653190187287D-02,   4.060849619196801D-03,  -4.051368747339830D-04, &
&    2.761111481851752D-05,  -1.231622827617666D-06,   3.249232324267357D-08, &
&   -3.851739792444945D-10,   6.769139442201537D-01,  -3.783518223820482D-01, &
&    1.214257096016929D-01,  -2.518540056502214D-02,   3.551170329362119D-03, &
&   -3.458487212194843D-04,   2.301316857550481D-05,  -1.002382719404352D-06, &
&    2.582501734241042D-08,  -2.989865561703131D-10,   6.563143660892078D-01, &
&   -3.572792987091244D-01,   1.118447759141316D-01,  -2.264424363481021D-02, &
&    3.117873830445191D-03,  -2.965922633845162D-04,   1.928010426958492D-05, &
&   -8.204973417272264D-07,   2.065535195607614D-08,  -2.336795023974850D-10, &
&    6.370587839081798D-01,  -3.380195885781189D-01,   1.032827760925571D-01, &
&   -2.042383776876847D-02,   2.747688779923439D-03,  -2.554461321436570D-04, &
&    1.623106686867063D-05,  -6.752440992628909D-07,   1.661870682317888D-08, &
&   -1.838200346615760D-10,   6.190120033975020D-01,  -3.203615420730855D-01, &
&    9.560355253336551D-02,  -1.847568500303102D-02,   2.429959362224784D-03, &
&   -2.208987753263688D-04,   1.372671585753052D-05,  -5.585351186952662D-07, &
&    1.344589462985757D-08,  -1.454832089595699D-10,   6.020569792050626D-01, &
&   -3.041248766244393D-01,   8.869278306835553D-02,  -1.675980697559348D-02, &
&    2.156070393944811D-03,  -1.917523666799123D-04,   1.165887107556909D-05, &
&   -4.642205187486950D-07,   1.093649208022802D-08,  -1.158080645769749D-10, &
&    5.860918095088875D-01,  -2.891547845930740D-01,   8.245391253436112D-02, &
&   -1.524304526077679D-02,   1.919011474773674D-03,  -1.670512229576893D-04, &
&    9.942940784081151D-06,  -3.875885889643397D-07,   8.940082061173694D-09, &
&   -9.269166157280035D-11,   5.710273159428749D-01,  -2.753176292173221D-01, &
&    7.680494417517188D-02,  -1.389774006501874D-02,   1.713043581841007D-03, &
&   -1.460279888984735D-04,   8.512332244579305D-06,  -3.250037532305761D-07, &
&    7.342928218127729D-09,  -7.457599172953888D-11/
!
      val1= abs(val)
      if(val1 < ten) then
        igrid= aint(val1*griddaw)
        val2= val1*val1
        val3= val1*val1*val1
        val4= val2*val2
        val5= val2*val3
        val6= val3*val3
        val7= val3*val4
        val8= val4*val4
        val9= val4*val4*val1
        dawerf= dawerfgrid(0,igrid)     +dawerfgrid(1,igrid)*val1+dawerfgrid(2,igrid)*val2 &
&              +dawerfgrid(3,igrid)*val3+dawerfgrid(4,igrid)*val4+dawerfgrid(5,igrid)*val5 &
&              +dawerfgrid(6,igrid)*val6+dawerfgrid(7,igrid)*val7+dawerfgrid(8,igrid)*val8 &
&              +dawerfgrid(9,igrid)*val9
      else
        val2inv= one/(val1*val1)
        last= 9
        factor=(last+half)*val2inv
        total= factor
        do i= last,1,-1
          factor= factor-val2inv
          total = factor*(one+total)
        enddo
        dawerf= val1*factor*(one+total)
      endif
!
      return
end








