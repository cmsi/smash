! Copyright 2021  Kazuya Ishimura
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
      use modtype, only : typecomp
      implicit none
      type(typecomp) :: datacomp
!
! Initialize and get MPI execution environment
!
      call setparallel(datacomp)
!
! Open input and output files if necessary
!
      if((datacomp%master).and.(command_argument_count() >= 1)) then
        call openinputfile(datacomp)
        call openoutputfile(datacomp)
      endif
!
! Main driver of SMASH
!
      call smashmain(datacomp)
!
! Close input and output files if necessary
!
      if((datacomp%master).and.(command_argument_count() >= 1)) then
        call closeinputfile(datacomp)
        call closeoutputfile(datacomp)
      endif
!
! Finalize MPI
!
      call para_finalize
!
end program smash
