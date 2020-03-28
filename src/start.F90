! Copyright 2014-2020  Kazuya Ishimura
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
! implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!----------------
  program smash
!----------------
!
! Main program of Scalable Molecular Analysis Solver
! for High performance computing systems (SMASH).
!
      use modtype, only : typejob, typemol, typebasis, typecomp
      implicit none
      type(typejob) :: datajob
      type(typemol) :: datamol
      type(typebasis) :: databasis
      type(typecomp) :: datacomp
!
! Initialize and get MPI execution environment
!
      call setparallel(datacomp)
!
! Main driver of SMASH
!
      call smashmain(datajob,datamol,databasis,datacomp)
!
end program smash
