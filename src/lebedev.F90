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
!
!
!--------------------------------------------
  subroutine lebedevquad1(angpt,weight,num)
!--------------------------------------------
!
!  Generate Lebedev grids for integration on a spher.
!  V.I. Lebedev, and D.N. Laikov, Doklady Mathematics, 59 (1999) 477-481.
!
      implicit none
      integer,intent(inout) :: num
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: weight
      real(8),intent(out) :: angpt(4,*)
!
      angpt(1,1)= one
      angpt(2,1)= zero
      angpt(3,1)= zero
      angpt(4,1)= weight
      angpt(1,2)=-one
      angpt(2,2)= zero
      angpt(3,2)= zero
      angpt(4,2)= weight
      angpt(1,3)= zero
      angpt(2,3)= one
      angpt(3,3)= zero
      angpt(4,3)= weight
      angpt(1,4)= zero
      angpt(2,4)=-one
      angpt(3,4)= zero
      angpt(4,4)= weight
      angpt(1,5)= zero
      angpt(2,5)= zero
      angpt(3,5)= one
      angpt(4,5)= weight
      angpt(1,6)= zero
      angpt(2,6)= zero
      angpt(3,6)=-one
      angpt(4,6)= weight
!
      num= num+6
!
      return
end


!--------------------------------------------
  subroutine lebedevquad2(angpt,weight,num)
!--------------------------------------------
!
!  Generate Lebedev grids for integration on a spher.
!  V.I. Lebedev, and D.N. Laikov, Doklady Mathematics, 59 (1999) 477-481.
!
      implicit none
      integer,intent(inout) :: num
      real(8),parameter :: zero=0.0D+00, half=0.5D+00
      real(8),intent(in) :: weight
      real(8),intent(out) :: angpt(4,*)
      real(8) :: a
!
      a= sqrt(half)
      angpt(1, 1)= zero
      angpt(2, 1)= a
      angpt(3, 1)= a
      angpt(4, 1)= weight
      angpt(1, 2)= zero
      angpt(2, 2)=-a
      angpt(3, 2)= a
      angpt(4, 2)= weight
      angpt(1, 3)= zero
      angpt(2, 3)= a
      angpt(3, 3)=-a
      angpt(4, 3)= weight
      angpt(1, 4)= zero
      angpt(2, 4)=-a
      angpt(3, 4)=-a
      angpt(4, 4)= weight
      angpt(1, 5)= a
      angpt(2, 5)= zero
      angpt(3, 5)= a
      angpt(4, 5)= weight
      angpt(1, 6)=-a
      angpt(2, 6)= zero
      angpt(3, 6)= a
      angpt(4, 6)= weight
      angpt(1, 7)= a
      angpt(2, 7)= zero
      angpt(3, 7)=-a
      angpt(4, 7)= weight
      angpt(1, 8)=-a
      angpt(2, 8)= zero
      angpt(3, 8)=-a
      angpt(4, 8)= weight
      angpt(1, 9)= a
      angpt(2, 9)= a
      angpt(3, 9)= zero
      angpt(4, 9)= weight
      angpt(1,10)=-a
      angpt(2,10)= a
      angpt(3,10)= zero
      angpt(4,10)= weight
      angpt(1,11)= a
      angpt(2,11)=-a
      angpt(3,11)= zero
      angpt(4,11)= weight
      angpt(1,12)=-a
      angpt(2,12)=-a
      angpt(3,12)= zero
      angpt(4,12)= weight
!
      num= num+12
!
      return
end


!--------------------------------------------
  subroutine lebedevquad3(angpt,weight,num)
!--------------------------------------------
!
!  Generate Lebedev grids for integration on a spher.
!  V.I. Lebedev, and D.N. Laikov, Doklady Mathematics, 59 (1999) 477-481.
!
      implicit none
      integer,intent(inout) :: num
      real(8),parameter :: one=1.0D+00, three=3.0D+00
      real(8),intent(in) :: weight
      real(8),intent(out) :: angpt(4,*)
      real(8) :: a
!
      a=sqrt(one/three)
      angpt(1,1)= a
      angpt(2,1)= a
      angpt(3,1)= a
      angpt(4,1)= weight
      angpt(1,2)=-a
      angpt(2,2)= a
      angpt(3,2)= a
      angpt(4,2)= weight
      angpt(1,3)= a
      angpt(2,3)=-a
      angpt(3,3)= a
      angpt(4,3)= weight
      angpt(1,4)=-a
      angpt(2,4)=-a
      angpt(3,4)= a
      angpt(4,4)= weight
      angpt(1,5)= a
      angpt(2,5)= a
      angpt(3,5)=-a
      angpt(4,5)= weight
      angpt(1,6)=-a
      angpt(2,6)= a
      angpt(3,6)=-a
      angpt(4,6)= weight
      angpt(1,7)= a
      angpt(2,7)=-a
      angpt(3,7)=-a
      angpt(4,7)= weight
      angpt(1,8)=-a
      angpt(2,8)=-a
      angpt(3,8)=-a
      angpt(4,8)= weight
!
      num= num+8
!
      return
end


!----------------------------------------------
  subroutine lebedevquad4(angpt,a,weight,num)
!----------------------------------------------
!
!  Generate Lebedev grids for integration on a spher.
!  V.I. Lebedev, and D.N. Laikov, Doklady Mathematics, 59 (1999) 477-481.
!
      implicit none
      integer,intent(inout) :: num
      real(8),parameter :: one=1.0D+00, two=2.0D+00
      real(8),intent(in) :: a, weight
      real(8),intent(out) :: angpt(4,*)
      real(8) :: b
!
      b= sqrt(one-two*a*a)
      angpt(1, 1) = a
      angpt(2, 1) = a
      angpt(3, 1) = b
      angpt(4, 1) = weight
      angpt(1, 2) =-a
      angpt(2, 2) = a
      angpt(3, 2) = b
      angpt(4, 2) = weight
      angpt(1, 3) = a
      angpt(2, 3) =-a
      angpt(3, 3) = b
      angpt(4, 3) = weight
      angpt(1, 4) =-a
      angpt(2, 4) =-a
      angpt(3, 4) = b
      angpt(4, 4) = weight
      angpt(1, 5) = a
      angpt(2, 5) = a
      angpt(3, 5) =-b
      angpt(4, 5) = weight
      angpt(1, 6) =-a
      angpt(2, 6) = a
      angpt(3, 6) =-b
      angpt(4, 6) = weight
      angpt(1, 7) = a
      angpt(2, 7) =-a
      angpt(3, 7) =-b
      angpt(4, 7) = weight
      angpt(1, 8) =-a
      angpt(2, 8) =-a
      angpt(3, 8) =-b
      angpt(4, 8) = weight
      angpt(1, 9) = a
      angpt(2, 9) = b
      angpt(3, 9) = a
      angpt(4, 9) = weight
      angpt(1,10) =-a
      angpt(2,10) = b
      angpt(3,10) = a
      angpt(4,10) = weight
      angpt(1,11) = a
      angpt(2,11) =-b
      angpt(3,11) = a
      angpt(4,11) = weight
      angpt(1,12) =-a
      angpt(2,12) =-b
      angpt(3,12) = a
      angpt(4,12) = weight
      angpt(1,13) = a
      angpt(2,13) = b
      angpt(3,13) =-a
      angpt(4,13) = weight
      angpt(1,14) =-a
      angpt(2,14) = b
      angpt(3,14) =-a
      angpt(4,14) = weight
      angpt(1,15) = a
      angpt(2,15) =-b
      angpt(3,15) =-a
      angpt(4,15) = weight
      angpt(1,16) =-a
      angpt(2,16) =-b
      angpt(3,16) =-a
      angpt(4,16) = weight
      angpt(1,17) = b
      angpt(2,17) = a
      angpt(3,17) = a
      angpt(4,17) = weight
      angpt(1,18) =-b
      angpt(2,18) = a
      angpt(3,18) = a
      angpt(4,18) = weight
      angpt(1,19) = b
      angpt(2,19) =-a
      angpt(3,19) = a
      angpt(4,19) = weight
      angpt(1,20) =-b
      angpt(2,20) =-a
      angpt(3,20) = a
      angpt(4,20) = weight
      angpt(1,21) = b
      angpt(2,21) = a
      angpt(3,21) =-a
      angpt(4,21) = weight
      angpt(1,22) =-b
      angpt(2,22) = a
      angpt(3,22) =-a
      angpt(4,22) = weight
      angpt(1,23) = b
      angpt(2,23) =-a
      angpt(3,23) =-a
      angpt(4,23) = weight
      angpt(1,24) =-b
      angpt(2,24) =-a
      angpt(3,24) =-a
      angpt(4,24) = weight
!
      num= num+24
!
      return
end


!----------------------------------------------
  subroutine lebedevquad5(angpt,a,weight,num)
!----------------------------------------------
!
!  Generate Lebedev grids for integration on a spher.
!  V.I. Lebedev, and D.N. Laikov, Doklady Mathematics, 59 (1999) 477-481.
!
      implicit none
      integer,intent(inout) :: num
      real(8),parameter :: zero=0.0D+00, one=1.0D+00
      real(8),intent(in) :: a, weight
      real(8),intent(out) :: angpt(4,*)
      real(8) :: b
!
      b= sqrt(one-a*a)
      angpt(1, 1)= a
      angpt(2, 1)= b
      angpt(3, 1)= zero
      angpt(4, 1)= weight
      angpt(1, 2)=-a
      angpt(2, 2)= b
      angpt(3, 2)= zero
      angpt(4, 2)= weight
      angpt(1, 3)= a
      angpt(2, 3)=-b
      angpt(3, 3)= zero
      angpt(4, 3)= weight
      angpt(1, 4)=-a
      angpt(2, 4)=-b
      angpt(3, 4)= zero
      angpt(4, 4)= weight
      angpt(1, 5)= b
      angpt(2, 5)= a
      angpt(3, 5)= zero
      angpt(4, 5)= weight
      angpt(1, 6)=-b
      angpt(2, 6)= a
      angpt(3, 6)= zero
      angpt(4, 6)= weight
      angpt(1, 7)= b
      angpt(2, 7)=-a
      angpt(3, 7)= zero
      angpt(4, 7)= weight
      angpt(1, 8)=-b
      angpt(2, 8)=-a
      angpt(3, 8)= zero
      angpt(4, 8)= weight
      angpt(1, 9)= a
      angpt(2, 9)= zero
      angpt(3, 9)= b
      angpt(4, 9)= weight
      angpt(1,10)=-a
      angpt(2,10)= zero
      angpt(3,10)= b
      angpt(4,10)= weight
      angpt(1,11)= a
      angpt(2,11)= zero
      angpt(3,11)=-b
      angpt(4,11)= weight
      angpt(1,12)=-a
      angpt(2,12)= zero
      angpt(3,12)=-b
      angpt(4,12)= weight
      angpt(1,13)= b
      angpt(2,13)= zero
      angpt(3,13)= a
      angpt(4,13)= weight
      angpt(1,14)=-b
      angpt(2,14)= zero
      angpt(3,14)= a
      angpt(4,14)= weight
      angpt(1,15)= b
      angpt(2,15)= zero
      angpt(3,15)=-a
      angpt(4,15)= weight
      angpt(1,16)=-b
      angpt(2,16)= zero
      angpt(3,16)=-a
      angpt(4,16)= weight
      angpt(1,17)= zero
      angpt(2,17)= a
      angpt(3,17)= b
      angpt(4,17)= weight
      angpt(1,18)= zero
      angpt(2,18)=-a
      angpt(3,18)= b
      angpt(4,18)= weight
      angpt(1,19)= zero
      angpt(2,19)= a
      angpt(3,19)=-b
      angpt(4,19)= weight
      angpt(1,20)= zero
      angpt(2,20)=-a
      angpt(3,20)=-b
      angpt(4,20)= weight
      angpt(1,21)= zero
      angpt(2,21)= b
      angpt(3,21)= a
      angpt(4,21)= weight
      angpt(1,22)= zero
      angpt(2,22)=-b
      angpt(3,22)= a
      angpt(4,22)= weight
      angpt(1,23)= zero
      angpt(2,23)= b
      angpt(3,23)=-a
      angpt(4,23)= weight
      angpt(1,24)= zero
      angpt(2,24)=-b
      angpt(3,24)=-a
      angpt(4,24)= weight
!
      num= num+24
!
      return
end

!------------------------------------------------
  subroutine lebedevquad6(angpt,a,b,weight,num)
!------------------------------------------------
!
!  Generate Lebedev grids for integration on a spher.
!  V.I. Lebedev, and D.N. Laikov, Doklady Mathematics, 59 (1999) 477-481.
!
      implicit none
      integer,intent(inout) :: num
      real(8),parameter :: one=1.0D+00
      real(8),intent(in) :: a, b, weight
      real(8),intent(out) :: angpt(4,*)
      real(8) :: c
!
      c= sqrt(one-a*a-b*b)
      angpt(1, 1)= a
      angpt(2, 1)= b
      angpt(3, 1)= c
      angpt(4, 1)= weight
      angpt(1, 2)=-a
      angpt(2, 2)= b
      angpt(3, 2)= c
      angpt(4, 2)= weight
      angpt(1, 3)= a
      angpt(2, 3)=-b
      angpt(3, 3)= c
      angpt(4, 3)= weight
      angpt(1, 4)=-a
      angpt(2, 4)=-b
      angpt(3, 4)= c
      angpt(4, 4)= weight
      angpt(1, 5)= a
      angpt(2, 5)= b
      angpt(3, 5)=-c
      angpt(4, 5)= weight
      angpt(1, 6)=-a
      angpt(2, 6)= b
      angpt(3, 6)=-c
      angpt(4, 6)= weight
      angpt(1, 7)= a
      angpt(2, 7)=-b
      angpt(3, 7)=-c
      angpt(4, 7)= weight
      angpt(1, 8)=-a
      angpt(2, 8)=-b
      angpt(3, 8)=-c
      angpt(4, 8)= weight
      angpt(1, 9)= a
      angpt(2, 9)= c
      angpt(3, 9)= b
      angpt(4, 9)= weight
      angpt(1,10)=-a
      angpt(2,10)= c
      angpt(3,10)= b
      angpt(4,10)= weight
      angpt(1,11)= a
      angpt(2,11)=-c
      angpt(3,11)= b
      angpt(4,11)= weight
      angpt(1,12)=-a
      angpt(2,12)=-c
      angpt(3,12)= b
      angpt(4,12)= weight
      angpt(1,13)= a
      angpt(2,13)= c
      angpt(3,13)=-b
      angpt(4,13)= weight
      angpt(1,14)=-a
      angpt(2,14)= c
      angpt(3,14)=-b
      angpt(4,14)= weight
      angpt(1,15)= a
      angpt(2,15)=-c
      angpt(3,15)=-b
      angpt(4,15)= weight
      angpt(1,16)=-a
      angpt(2,16)=-c
      angpt(3,16)=-b
      angpt(4,16)= weight
      angpt(1,17)= b
      angpt(2,17)= a
      angpt(3,17)= c
      angpt(4,17)= weight
      angpt(1,18)=-b
      angpt(2,18)= a
      angpt(3,18)= c
      angpt(4,18)= weight
      angpt(1,19)= b
      angpt(2,19)=-a
      angpt(3,19)= c
      angpt(4,19)= weight
      angpt(1,20)=-b
      angpt(2,20)=-a
      angpt(3,20)= c
      angpt(4,20)= weight
      angpt(1,21)= b
      angpt(2,21)= a
      angpt(3,21)=-c
      angpt(4,21)= weight
      angpt(1,22)=-b
      angpt(2,22)= a
      angpt(3,22)=-c
      angpt(4,22)= weight
      angpt(1,23)= b
      angpt(2,23)=-a
      angpt(3,23)=-c
      angpt(4,23)= weight
      angpt(1,24)=-b
      angpt(2,24)=-a
      angpt(3,24)=-c
      angpt(4,24)= weight
      angpt(1,25)= b
      angpt(2,25)= c
      angpt(3,25)= a
      angpt(4,25)= weight
      angpt(1,26)=-b
      angpt(2,26)= c
      angpt(3,26)= a
      angpt(4,26)= weight
      angpt(1,27)= b
      angpt(2,27)=-c
      angpt(3,27)= a
      angpt(4,27)= weight
      angpt(1,28)=-b
      angpt(2,28)=-c
      angpt(3,28)= a
      angpt(4,28)= weight
      angpt(1,29)= b
      angpt(2,29)= c
      angpt(3,29)=-a
      angpt(4,29)= weight
      angpt(1,30)=-b
      angpt(2,30)= c
      angpt(3,30)=-a
      angpt(4,30)= weight
      angpt(1,31)= b
      angpt(2,31)=-c
      angpt(3,31)=-a
      angpt(4,31)= weight
      angpt(1,32)=-b
      angpt(2,32)=-c
      angpt(3,32)=-a
      angpt(4,32)= weight
      angpt(1,33)= c
      angpt(2,33)= a
      angpt(3,33)= b
      angpt(4,33)= weight
      angpt(1,34)=-c
      angpt(2,34)= a
      angpt(3,34)= b
      angpt(4,34)= weight
      angpt(1,35)= c
      angpt(2,35)=-a
      angpt(3,35)= b
      angpt(4,35)= weight
      angpt(1,36)=-c
      angpt(2,36)=-a
      angpt(3,36)= b
      angpt(4,36)= weight
      angpt(1,37)= c
      angpt(2,37)= a
      angpt(3,37)=-b
      angpt(4,37)= weight
      angpt(1,38)=-c
      angpt(2,38)= a
      angpt(3,38)=-b
      angpt(4,38)= weight
      angpt(1,39)= c
      angpt(2,39)=-a
      angpt(3,39)=-b
      angpt(4,39)= weight
      angpt(1,40)=-c
      angpt(2,40)=-a
      angpt(3,40)=-b
      angpt(4,40)= weight
      angpt(1,41)= c
      angpt(2,41)= b
      angpt(3,41)= a
      angpt(4,41)= weight
      angpt(1,42)=-c
      angpt(2,42)= b
      angpt(3,42)= a
      angpt(4,42)= weight
      angpt(1,43)= c
      angpt(2,43)=-b
      angpt(3,43)= a
      angpt(4,43)= weight
      angpt(1,44)=-c
      angpt(2,44)=-b
      angpt(3,44)= a
      angpt(4,44)= weight
      angpt(1,45)= c
      angpt(2,45)= b
      angpt(3,45)=-a
      angpt(4,45)= weight
      angpt(1,46)=-c
      angpt(2,46)= b
      angpt(3,46)=-a
      angpt(4,46)= weight
      angpt(1,47)= c
      angpt(2,47)=-b
      angpt(3,47)=-a
      angpt(4,47)= weight
      angpt(1,48)=-c
      angpt(2,48)=-b
      angpt(3,48)=-a
      angpt(4,48)= weight
!
      num= num+48
!
      return
end


!-------------------------------
  subroutine lebedev6(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,6)
      real(8) :: weight
!
      num= 1
      weight= 0.1666666666666667D+00
      call lebedevquad1(angpt(1,num),weight,num)
!
      return
end


!-------------------------------
  subroutine lebedev14(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,14)
      real(8) :: weight
!
      num= 1
      weight= 0.6666666666666667D-01
      call lebedevquad1(angpt(1,num),weight,num)
      weight= 0.7500000000000000D-01
      call lebedevquad3(angpt(1,num),weight,num)
!
      return
end


!-------------------------------
  subroutine lebedev26(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,26)
      real(8) :: weight
!
      num= 1
      weight= 0.4761904761904762D-01
      call lebedevquad1(angpt(1,num),weight,num)
      weight= 0.3809523809523810D-01
      call lebedevquad2(angpt(1,num),weight,num)
      weight= 0.3214285714285714D-01
      call lebedevquad3(angpt(1,num),weight,num)
!
      return
end


!-------------------------------
  subroutine lebedev38(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,38)
      real(8) :: a, weight
!
      num= 1
      weight= 0.9523809523809524D-02
      call lebedevquad1(angpt(1,num),weight,num)
      weight= 0.3214285714285714D-01
      call lebedevquad3(angpt(1,num),weight,num)
      a= 0.4597008433809831D+00
      weight= 0.2857142857142857D-01
      call lebedevquad5(angpt(1,num),a,weight,num)
!
      return
end


!-------------------------------
  subroutine lebedev50(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,50)
      real(8) :: a, weight
!
      num= 1
      weight= 0.1269841269841270D-01
      call lebedevquad1(angpt(1,num),weight,num)
      weight= 0.2257495590828924D-01
      call lebedevquad2(angpt(1,num),weight,num)
      weight= 0.2109375000000000D-01
      call lebedevquad3(angpt(1,num),weight,num)
      a= 0.3015113445777636D+00
      weight= 0.2017333553791887D-01
      call lebedevquad4(angpt(1,num),a,weight,num)
!
      return
end


!-------------------------------
  subroutine lebedev74(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,74)
      real(8) :: a, weight
!
      num= 1
      weight= 0.5130671797338464D-03
      call lebedevquad1(angpt(1,num),weight,num)
      weight= 0.1660406956574204D-01
      call lebedevquad2(angpt(1,num),weight,num)
      weight=-0.2958603896103896D-01
      call lebedevquad3(angpt(1,num),weight,num)
      a= 0.4803844614152614D+00
      weight= 0.2657620708215946D-01
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.3207726489807764D+00
      weight= 0.1652217099371571D-01
      call lebedevquad5(angpt(1,num),a,weight,num)
!
      return
end


!-------------------------------
  subroutine lebedev86(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,86)
      real(8) :: a, weight
!
      num= 1
      weight= 0.1154401154401154D-01
      call lebedevquad1(angpt(1,num),weight,num)
      weight= 0.1194390908585628D-01
      call lebedevquad3(angpt(1,num),weight,num)
      a= 0.3696028464541502D+00
      weight= 0.1111055571060340D-01
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6943540066026664D+00
      weight= 0.1187650129453714D-01
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.3742430390903412D+00
      weight= 0.1181230374690448D-01
      call lebedevquad5(angpt(1,num),a,weight,num)
!
      return
end


!-------------------------------
  subroutine lebedev110(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,110)
      real(8) :: a, weight
!
      num= 1
      weight= 0.3828270494937162D-02
      call lebedevquad1(angpt(1,num),weight,num)
      weight= 0.9793737512487512D-02
      call lebedevquad3(angpt(1,num),weight,num)
      a= 0.1851156353447362D+00
      weight= 0.8211737283191111D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6904210483822922D+00
      weight= 0.9942814891178103D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.3956894730559419D+00
      weight= 0.9595471336070963D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.4783690288121502D+00
      weight= 0.9694996361663028D-02
      call lebedevquad5(angpt(1,num),a,weight,num)
!
      return
end


!-------------------------------
  subroutine lebedev146(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,146)
      real(8) :: a, b, weight
!
      num= 1
      weight= 0.5996313688621381D-03
      call lebedevquad1(angpt(1,num),weight,num)
      weight= 0.7372999718620756D-02
      call lebedevquad2(angpt(1,num),weight,num)
      weight= 0.7210515360144488D-02
      call lebedevquad3(angpt(1,num),weight,num)
      a= 0.6764410400114264D+00
      weight= 0.7116355493117555D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.4174961227965453D+00
      weight= 0.6753829486314477D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.1574676672039082D+00
      weight= 0.7574394159054034D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.1403553811713183D+00
      b= 0.4493328323269557D+00
      weight= 0.6991087353303262D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
!
      return
end


!-------------------------------
  subroutine lebedev170(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,170)
      real(8) :: a, b, weight
!
      num= 1
      weight= 0.5544842902037365D-02
      call lebedevquad1(angpt(1,num),weight,num)
      weight= 0.6071332770670752D-02
      call lebedevquad2(angpt(1,num),weight,num)
      weight= 0.6383674773515093D-02
      call lebedevquad3(angpt(1,num),weight,num)
      a= 0.2551252621114134D+00
      weight= 0.5183387587747790D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6743601460362766D+00
      weight= 0.6317929009813725D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.4318910696719410D+00
      weight= 0.6201670006589077D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.2613931360335988D+00
      weight= 0.5477143385137348D-02
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.4990453161796037D+00
      b= 0.1446630744325115D+00
      weight= 0.5968383987681156D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
!
      return
end


!-------------------------------
  subroutine lebedev194(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,194)
      real(8) :: a, b, weight
!
      num= 1
      weight= 0.1782340447244611D-02
      call lebedevquad1(angpt(1,num),weight,num)
      weight= 0.5716905949977102D-02
      call lebedevquad2(angpt(1,num),weight,num)
      weight= 0.5573383178848738D-02
      call lebedevquad3(angpt(1,num),weight,num)
      a= 0.6712973442695226D+00
      weight= 0.5608704082587997D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.2892465627575439D+00
      weight= 0.5158237711805383D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.4446933178717437D+00
      weight= 0.5518771467273614D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.1299335447650067D+00
      weight= 0.4106777028169394D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.3457702197611283D+00
      weight= 0.5051846064614808D-02
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.1590417105383530D+00
      b= 0.8360360154824589D+00
      weight= 0.5530248916233094D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
!
      return
end


!-------------------------------
  subroutine lebedev230(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,230)
      real(8) :: a, b, weight
!
      num= 1
      weight=-0.5522639919727325D-01
      call lebedevquad1(angpt(1,num),weight,num)
      weight= 0.4450274607445226D-02
      call lebedevquad3(angpt(1,num),weight,num)
      a= 0.4492044687397611D+00
      weight= 0.4496841067921404D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.2520419490210201D+00
      weight= 0.5049153450478750D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6981906658447242D+00
      weight= 0.3976408018051883D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6587405243460960D+00
      weight= 0.4401400650381014D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.4038544050097660D-01
      weight= 0.1724544350544401D-01
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.5823842309715585D+00
      weight= 0.4231083095357343D-02
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.3545877390518688D+00
      weight= 0.5198069864064399D-02
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.2272181808998187D+00
      b= 0.4864661535886647D+00
      weight= 0.4695720972568883D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
!
      return
end


!-------------------------------
  subroutine lebedev266(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,266)
      real(8) :: a, b, weight
!
      num= 1
      weight=-0.1313769127326952D-02
      call lebedevquad1(angpt(1,num),weight,num)
      weight=-0.2522728704859336D-02
      call lebedevquad2(angpt(1,num),weight,num)
      weight= 0.4186853881700583D-02
      call lebedevquad3(angpt(1,num),weight,num)
      a= 0.7039373391585475D+00
      weight= 0.5315167977810885D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.1012526248572414D+00
      weight= 0.4047142377086219D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.4647448726420539D+00
      weight= 0.4112482394406990D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.3277420654971629D+00
      weight= 0.3595584899758782D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6620338663699974D+00
      weight= 0.4256131351428158D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.8506508083520399D+00
      weight= 0.4229582700647240D-02
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.3233484542692899D+00
      b= 0.1153112011009701D+00
      weight= 0.4080914225780505D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.2314790158712601D+00
      b= 0.5244939240922365D+00
      weight= 0.4071467593830964D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
!
      return
end


!-------------------------------
  subroutine lebedev302(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,302)
      real(8) :: a, b, weight
!
      num= 1
      weight= 0.8545911725128148D-03
      call lebedevquad1(angpt(1,num),weight,num)
      weight= 0.3599119285025571D-02
      call lebedevquad3(angpt(1,num),weight,num)
      a= 0.3515640345570105D+00
      weight= 0.3449788424305883D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6566329410219612D+00
      weight= 0.3604822601419882D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.4729054132581005D+00
      weight= 0.3576729661743367D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.9618308522614784D-01
      weight= 0.2352101413689164D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.2219645236294178D+00
      weight= 0.3108953122413675D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.7011766416089545D+00
      weight= 0.3650045807677255D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.2644152887060663D+00
      weight= 0.2982344963171804D-02
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.5718955891878961D+00
      weight= 0.3600820932216460D-02
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.2510034751770465D+00
      b= 0.8000727494073952D+00
      weight= 0.3571540554273387D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.1233548532583327D+00
      b= 0.4127724083168531D+00
      weight= 0.3392312205006170D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
!
      return
end



!-------------------------------
  subroutine lebedev350(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,350)
      real(8) :: a, b, weight
!
      num= 1
      weight= 0.3006796749453936D-02
      call lebedevquad1(angpt(1,num),weight,num)
      weight= 0.3050627745650771D-02
      call lebedevquad3(angpt(1,num),weight,num)
      a= 0.7068965463912316D+00
      weight= 0.1621104600288991D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a=  0.4794682625712025D+00
      weight= 0.3005701484901752D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.1927533154878019D+00
      weight= 0.2990992529653774D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6930357961327123D+00
      weight= 0.2982170644107595D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.3608302115520091D+00
      weight= 0.2721564237310992D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6498486161496169D+00
      weight= 0.3033513795811141D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.1932945013230339D+00
      weight= 0.3007949555218533D-02
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.3800494919899303D+00
      weight= 0.2881964603055307D-02
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.2899558825499574D+00
      b= 0.7934537856582316D+00
      weight= 0.2958357626535696D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.9684121455103957D-01
      b= 0.8280801506686862D+00
      weight= 0.3036020026407088D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.1833434647041659D+00
      b= 0.9074658265305127D+00
      weight= 0.2832187403926303D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
!
      return
end


!-------------------------------
  subroutine lebedev434(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,434)
      real(8) :: a, b, weight
!
      num= 1
      weight= 0.5265897968224436D-03
      call lebedevquad1(angpt(1,num),weight,num)
      weight= 0.2548219972002607D-02
      call lebedevquad2(angpt(1,num),weight,num)
      weight= 0.2512317418927307D-02
      call lebedevquad3(angpt(1,num),weight,num)
      a= 0.6909346307509111D+00
      weight= 0.2530403801186355D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.1774836054609158D+00
      weight= 0.2014279020918528D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.4914342637784746D+00
      weight= 0.2501725168402936D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6456664707424256D+00
      weight= 0.2513267174597564D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.2861289010307638D+00
      weight= 0.2302694782227416D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.7568084367178018D-01
      weight= 0.1462495621594614D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.3927259763368002D+00
      weight= 0.2445373437312980D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.8818132877794288D+00
      weight= 0.2417442375638981D-02
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.9776428111182649D+00
      weight= 0.1910951282179532D-02
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.2054823696403044D+00
      b= 0.8689460322872412D+00
      weight= 0.2416930044324775D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.5905157048925271D+00
      b= 0.7999278543857286D+00
      weight= 0.2512236854563495D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.5550152361076807D+00
      b= 0.7717462626915901D+00
      weight= 0.2496644054553086D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.9371809858553722D+00
      b= 0.3344363145343455D+00
      weight= 0.2236607760437849D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
!
      return
end


!-------------------------------
  subroutine lebedev590(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,590)
      real(8) :: a, b, weight
!
      num= 1
      weight= 0.3095121295306187D-03
      call lebedevquad1(angpt(1,num),weight,num)
      weight= 0.1852379698597489D-02
      call lebedevquad3(angpt(1,num),weight,num)
      a= 0.7040954938227469D+00
      weight= 0.1871790639277744D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6807744066455243D+00
      weight= 0.1858812585438317D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6372546939258752D+00
      weight= 0.1852028828296213D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.5044419707800358D+00
      weight= 0.1846715956151242D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.4215761784010967D+00
      weight= 0.1818471778162769D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.3317920736472123D+00
      weight= 0.1749564657281154D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.2384736701421887D+00
      weight= 0.1617210647254411D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.1459036449157763D+00
      weight= 0.1384737234851692D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6095034115507196D-01
      weight= 0.9764331165051050D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6116843442009876D+00
      weight= 0.1857161196774078D-02
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.3964755348199858D+00
      weight= 0.1705153996395864D-02
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.1724782009907724D+00
      weight= 0.1300321685886048D-02
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.5610263808622060D+00
      b= 0.3518280927733519D+00
      weight= 0.1842866472905286D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.4742392842551980D+00
      b= 0.2634716655937950D+00
      weight= 0.1802658934377451D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.5984126497885380D+00
      b= 0.1816640840360209D+00
      weight= 0.1849830560443660D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.3791035407695563D+00
      b= 0.1720795225656878D+00
      weight= 0.1713904507106709D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.2778673190586244D+00
      b= 0.8213021581932511D-01
      weight= 0.1555213603396808D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.5033564271075117D+00
      b= 0.8999205842074875D-01
      weight= 0.1802239128008525D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
!
      return
end


!-------------------------------
  subroutine lebedev770(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,770)
      real(8) :: a, b, weight
!
      num= 1
      weight= 0.2192942088181184D-03
      call lebedevquad1(angpt(1,num),weight,num)
      weight= 0.1436433617319080D-02
      call lebedevquad2(angpt(1,num),weight,num)
      weight= 0.1421940344335877D-02
      call lebedevquad3(angpt(1,num),weight,num)
      a= 0.5087204410502360D-01
      weight= 0.6798123511050502D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.1228198790178831D+00
      weight= 0.9913184235294912D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.2026890814408786D+00
      weight= 0.1180207833238949D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.2847745156464294D+00
      weight= 0.1296599602080921D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.3656719078978026D+00
      weight= 0.1365871427428316D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.4428264886713469D+00
      weight= 0.1402988604775325D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.5140619627249735D+00
      weight= 0.1418645563595609D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6306401219166803D+00
      weight= 0.1421376741851662D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6716883332022612D+00
      weight= 0.1423996475490962D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6979792685336881D+00
      weight= 0.1431554042178567D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.1446865674195309D+00
      weight= 0.9254401499865368D-03
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.3390263475411216D+00
      weight= 0.1250239995053509D-02
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.5335804651263506D+00
      weight= 0.1394365843329230D-02
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.6944024393349413D-01
      b= 0.2355187894242326D+00
      weight= 0.1127089094671749D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.2269004109529460D+00
      b= 0.4102182474045730D+00
      weight= 0.1345753760910670D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.8025574607775339D-01
      b= 0.6214302417481605D+00
      weight= 0.1424957283316783D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.1467999527896572D+00
      b= 0.3245284345717394D+00
      weight= 0.1261523341237750D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.1571507769824727D+00
      b= 0.5224482189696630D+00
      weight= 0.1392547106052696D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.2365702993157246D+00
      b= 0.6017546634089558D+00
      weight= 0.1418761677877656D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.7714815866765732D-01
      b= 0.4346575516141163D+00
      weight= 0.1338366684479554D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.3062936666210730D+00
      b= 0.4908826589037616D+00
      weight= 0.1393700862676131D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.3822477379524787D+00
      b= 0.5648768149099500D+00
      weight= 0.1415914757466932D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
!
      return
end


!-------------------------------
  subroutine lebedev974(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,974)
      real(8) :: a, b, weight
!
      num= 1
      weight= 0.1438294190527431D-03
      call lebedevquad1(angpt(1,num),weight,num)
      weight= 0.1125772288287004D-02
      call lebedevquad3(angpt(1,num),weight,num)
      a= 0.4292963545341347D-01
      weight= 0.4948029341949241D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.1051426854086404D+00
      weight= 0.7357990109125470D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.1750024867623087D+00
      weight= 0.8889132771304384D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.2477653379650257D+00
      weight= 0.9888347838921435D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.3206567123955957D+00
      weight= 0.1053299681709471D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.3916520749849983D+00
      weight= 0.1092778807014578D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.4590825874187624D+00
      weight= 0.1114389394063227D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.5214563888415861D+00
      weight= 0.1123724788051555D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6253170244654199D+00
      weight= 0.1125239325243814D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6637926744523170D+00
      weight= 0.1126153271815905D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6910410398498301D+00
      weight= 0.1130286931123841D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.7052907007457760D+00
      weight= 0.1134986534363955D-02
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.1236686762657990D+00
      weight= 0.6823367927109931D-03
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.2940777114468387D+00
      weight= 0.9454158160447096D-03
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.4697753849207649D+00
      weight= 0.1074429975385679D-02
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.6334563241139567D+00
      weight= 0.1129300086569132D-02
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.5974048614181342D-01
      b= 0.2029128752777523D+00
      weight= 0.8436884500901954D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.1375760408473636D+00
      b= 0.4602621942484054D+00
      weight= 0.1075255720448885D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.3391016526336286D+00
      b= 0.5030673999662036D+00
      weight= 0.1108577236864462D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.1271675191439820D+00
      b= 0.2817606422442134D+00
      weight= 0.9566475323783357D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.2693120740413512D+00
      b= 0.4331561291720157D+00
      weight= 0.1080663250717391D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.1419786452601918D+00
      b= 0.6256167358580814D+00
      weight= 0.1126797131196295D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.6709284600738255D-01
      b= 0.3798395216859157D+00
      weight= 0.1022568715358061D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.7057738183256172D-01
      b= 0.5517505421423520D+00
      weight= 0.1108960267713108D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.2783888477882155D+00
      b= 0.6029619156159187D+00
      weight= 0.1122790653435766D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.1979578938917407D+00
      b= 0.3589606329589096D+00
      weight= 0.1032401847117460D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.2087307061103274D+00
      b= 0.5348666438135476D+00
      weight= 0.1107249382283854D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.4055122137872836D+00
      b= 0.5674997546074373D+00
      weight= 0.1121780048519972D-02
      call lebedevquad6(angpt(1,num),a,b,weight,num)
!
      return
end


!-------------------------------
  subroutine lebedev1202(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,1202)
      real(8) :: a, b, weight
!
      num= 1
      weight= 0.1105189233267572D-03
      call lebedevquad1(angpt(1,num),weight,num)
      weight= 0.9205232738090741D-03
      call lebedevquad2(angpt(1,num),weight,num)
      weight= 0.9133159786443561D-03
      call lebedevquad3(angpt(1,num),weight,num)
      a= 0.3712636449657089D-01
      weight= 0.3690421898017899D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.9140060412262223D-01
      weight= 0.5603990928680660D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.1531077852469906D+00
      weight= 0.6865297629282609D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.2180928891660612D+00
      weight= 0.7720338551145630D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.2839874532200175D+00
      weight= 0.8301545958894795D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.3491177600963764D+00
      weight= 0.8686692550179628D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.4121431461444309D+00
      weight= 0.8927076285846890D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.4718993627149127D+00
      weight= 0.9060820238568219D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.5273145452842337D+00
      weight= 0.9119777254940867D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6209475332444019D+00
      weight= 0.9128720138604181D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6569722711857291D+00
      weight= 0.9130714935691735D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6841788309070143D+00
      weight= 0.9152873784554116D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.7012604330123631D+00
      weight= 0.9187436274321654D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.1072382215478166D+00
      weight= 0.5176977312965694D-03
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.2582068959496968D+00
      weight= 0.7331143682101417D-03
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.4172752955306717D+00
      weight= 0.8463232836379928D-03
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.5700366911792503D+00
      weight= 0.9031122694253992D-03
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.9827986018263947D+00
      b= 0.1771774022615325D+00
      weight= 0.6485778453163257D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.9624249230326228D+00
      b= 0.2475716463426288D+00
      weight= 0.7435030910982369D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.9402007994128811D+00
      b= 0.3354616289066489D+00
      weight= 0.7998527891839054D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.9320822040143202D+00
      b= 0.3173615246611977D+00
      weight= 0.8101731497468018D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.9043674199393299D+00
      b= 0.4090268427085357D+00
      weight= 0.8483389574594331D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.8912407560074747D+00
      b= 0.3854291150669224D+00
      weight= 0.8556299257311812D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.8676435628462708D+00
      b= 0.4932221184851285D+00
      weight= 0.8803208679738260D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.8581979986041619D+00
      b= 0.4785320675922435D+00
      weight= 0.8811048182425720D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.8396753624049856D+00
      b= 0.4507422593157064D+00
      weight= 0.8850282341265444D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.8165288564022188D+00
      b= 0.5632123020762100D+00
      weight= 0.9021342299040653D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.8015469370783529D+00
      b= 0.5434303569693900D+00
      weight= 0.9010091677105086D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.7773563069070351D+00
      b= 0.5123518486419871D+00
      weight= 0.9022692938426915D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.7661621213900394D+00
      b= 0.6394279634749102D+00
      weight= 0.9158016174693465D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.7553584143533510D+00
      b= 0.6269805509024392D+00
      weight= 0.9131578003189435D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.7344305757559503D+00
      b= 0.6031161693096310D+00
      weight= 0.9107813579482705D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.7043837184021765D+00
      b= 0.5693702498468441D+00
      weight= 0.9105760258970126D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
!
      return
end


!-------------------------------
  subroutine lebedev1454(angpt)
!-------------------------------
      implicit none
      integer :: num
      real(8),intent(out) :: angpt(4,1454)
      real(8) :: a, b, weight
!
      num= 1
      weight= 0.7777160743261247D-04
      call lebedevquad1(angpt(1,num),weight,num)
      weight= 0.7557646413004701D-03
      call lebedevquad3(angpt(1,num),weight,num)
      a= 0.3229290663413854D-01
      weight= 0.2841633806090617D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.8036733271462222D-01
      weight= 0.4374419127053555D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.1354289960531653D+00
      weight= 0.5417174740872172D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.1938963861114426D+00
      weight= 0.6148000891358593D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.2537343715011275D+00
      weight= 0.6664394485800705D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.3135251434752570D+00
      weight= 0.7025039356923220D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.3721558339375338D+00
      weight= 0.7268511789249627D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.4286809575195696D+00
      weight= 0.7422637534208629D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.4822510128282994D+00
      weight= 0.7509545035841214D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.5320679333566263D+00
      weight= 0.7548535057718401D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6172998195394274D+00
      weight= 0.7554088969774001D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6510679849127481D+00
      weight= 0.7553147174442808D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6777315251687360D+00
      weight= 0.7564767653292297D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.6963109410648741D+00
      weight= 0.7587991808518730D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.7058935009831749D+00
      weight= 0.7608261832033027D-03
      call lebedevquad4(angpt(1,num),a,weight,num)
      a= 0.9955546194091857D+00
      weight= 0.4021680447874916D-03
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.9734115901794209D+00
      weight= 0.5804871793945964D-03
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.9275693732388626D+00
      weight= 0.6792151955945159D-03
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.8568022422795103D+00
      weight= 0.7336741211286294D-03
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.7623495553719372D+00
      weight= 0.7581866300989608D-03
      call lebedevquad5(angpt(1,num),a,weight,num)
      a= 0.5707522908892223D+00
      b= 0.4387028039889501D+00
      weight= 0.7538257859800743D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.5196463388403083D+00
      b= 0.3858908414762617D+00
      weight= 0.7483517247053123D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.4646337531215351D+00
      b= 0.3301937372343854D+00
      weight= 0.7371763661112059D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.4063901697557691D+00
      b= 0.2725423573563777D+00
      weight= 0.7183448895756934D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.3456329466643087D+00
      b= 0.2139510237495250D+00
      weight= 0.6895815529822191D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.2831395121050332D+00
      b= 0.1555922309786647D+00
      weight= 0.6480105801792886D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.2197682022925330D+00
      b= 0.9892878979686097D-01
      weight= 0.5897558896594636D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.1564696098650355D+00
      b= 0.4598642910675510D-01
      weight= 0.5095708849247346D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.6027356673721295D+00
      b= 0.3376625140173426D+00
      weight= 0.7536906428909755D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.5496032320255096D+00
      b= 0.2822301309727988D+00
      weight= 0.7472505965575118D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.4921707755234567D+00
      b= 0.2248632342592540D+00
      weight= 0.7343017132279698D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.4309422998598483D+00
      b= 0.1666224723456479D+00
      weight= 0.7130871582177445D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.3664108182313672D+00
      b= 0.1086964901822169D+00
      weight= 0.6817022032112776D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.2990189057758436D+00
      b= 0.5251989784120085D-01
      weight= 0.6380941145604121D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.6268724013144998D+00
      b= 0.2297523657550023D+00
      weight= 0.7550381377920310D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.5707324144834607D+00
      b= 0.1723080607093800D+00
      weight= 0.7478646640144802D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.5096360901960365D+00
      b= 0.1140238465390513D+00
      weight= 0.7335918720601220D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.4438729938312456D+00
      b= 0.5611522095882537D-01
      weight= 0.7110120527658118D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.6419978471082389D+00
      b= 0.1164174423140873D+00
      weight= 0.7571363978689501D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
      a= 0.5817218061802611D+00
      b= 0.5797589531445219D-01
      weight= 0.7489908329079234D-03
      call lebedevquad6(angpt(1,num),a,b,weight,num)
!
      return
end

