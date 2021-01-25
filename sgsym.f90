! 
!     Copyright (C) 2021  Chee Kwan Gan (ihpcganck@gmail.com)
!     Copyright (C) 2020  Chee Kwan Gan (ihpcganck@gmail.com)
! 
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
! 
! 
module sgsym
  use types_and_constants
  implicit none

  character(len=*),parameter :: sgtable='table-sg'
  character(len=*),parameter :: dsgtable='table-dsg'
  character(len=*),parameter :: pgtable='table-pg'
  character(len=*),parameter :: pgchartab='table-pg-char'

  real(double),parameter :: dis_tol_tryput = 1.0d-6

  real(double),parameter :: dev_tol_get_supercell_atom_map = 1.0d-8

  real(double),parameter :: FCCmat(1:3,1:3) = reshape( (/ zero,half,half, &
                                             half,zero,half,&
                                             half,half,zero/), (/3,3/) )

  real(double),parameter :: BCCmat(1:3,1:3) = reshape( (/ -half,half,half, &
                                             half,-half,half,&
                                             half,half,-half/), (/3,3/) )

  real(double),parameter :: Mat2Trigonal(1:3,1:3) = reshape(  (/ third,-third,third, &
                                                      third, two*third, third,&
                                                   -two*third,-third,third    /), (/3,3/)  )

  real(double),parameter :: Mat6Trigonal(1:3,1:3) = reshape( (/two*third,third,third,&
                                                 -third,third,third,&
                                                 -third,-two*third,third/),(/3,3/))

  real(double),parameter :: BCTmat(1:3,1:3) = reshape( (/ half,-half,half, &
                                             half,half,half,&
                                             -half,-half,half/), (/3,3/) )

  real(double),parameter :: OBaseCmat1(1:3,1:3) = reshape( (/ half,half,zero, &
                                             -half,half,zero,&
                                             zero,zero,one/), (/3,3/) )

  real(double),parameter :: OBaseCmat2(1:3,1:3) = reshape( (/ half,-half,zero, &
                                             half,half,zero,&
                                             zero,zero,one/), (/3,3/) )

  real(double),parameter :: OBaseAmat(1:3,1:3) = reshape( (/ one,zero,zero, &
                                             zero,half,-half,&
                                             zero,half,half/), (/3,3/) )

  real(double),parameter :: OFCmat(1:3,1:3) = reshape( (/ half,zero,half, &
                                             half,half,zero,&
                                             zero,half,half/), (/3,3/) )

  real(double),parameter :: OBCmat(1:3,1:3) = reshape( (/ half,half,half, &
                                             -half,half,half,&
                                             -half,-half,half/), (/3,3/) )

  real(double),parameter :: MonocBaxisBaseCentredCmat(1:3,1:3) = reshape( (/ half,-half,zero, &
                                             half,half,zero,&
                                             zero,zero,one/), (/3,3/) )

contains

  subroutine init_static_sg_database(qtablepath)
    character(len=*) :: qtablepath
    integer :: i
    integer :: SG, Y

    call assign_env_variable("QTABLEPATH",qtablepath)
    write(*,*) 'qtablepath is '//trim(qtablepath)
    if(len(trim(qtablepath)) == 0) then
      write(*,'(A)') 'Environment QTABLEPATH not set, the symmetry tables must be copied manually to the the current directory'
      qtablepath="."
    endif

    do i = 0, 230
      SGbase(1:maxSGvariant,i)%sglabel = ''
      SGbase(1:maxSGvariant,i)%reorder_IT_op_sequence = 0
      SGbase(1:maxSGvariant,i)%assigned=.FALSE.
    enddo


    SG = 191
    Y = 1
    SGbase(Y,SG)%sglabel    ='191-P6/m2/m2/m'
    SGbase(Y,SG)%PG         ='D6h'
    SGbase(Y,SG)%ncentering =1
    SGbase(Y,SG)%centering(1:3,1) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%PGNumOper  =24
    SGbase(Y,SG)%symmorphic =1
    SGbase(Y,SG)%nop        =24
    SGbase(Y,SG)%assigned   =.TRUE.

    SG = 194
    Y = 1
    SGbase(Y,SG)%sglabel    ='194-P6_3/m2/m2/c'
    SGbase(Y,SG)%PG         ='D6h'
    SGbase(Y,SG)%ncentering =1
    SGbase(Y,SG)%centering(1:3,1) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%PGNumOper  =24
    SGbase(Y,SG)%symmorphic =0
    SGbase(Y,SG)%nsy(1:3,1) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,2) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,3) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,4) = (/           0.0d0,           0.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,5) = (/           0.0d0,           0.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,6) = (/           0.0d0,           0.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,7) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,8) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,9) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,10) = (/           0.0d0,           0.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,11) = (/           0.0d0,           0.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,12) = (/           0.0d0,           0.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,13) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,14) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,15) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,16) = (/           0.0d0,           0.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,17) = (/           0.0d0,           0.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,18) = (/           0.0d0,           0.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,19) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,20) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,21) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,22) = (/           0.0d0,           0.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,23) = (/           0.0d0,           0.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,24) = (/           0.0d0,           0.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nop        =24
    SGbase(Y,SG)%assigned   =.TRUE.


    SG = 227
    Y = 1
    SGbase(Y,SG)%sglabel    ='227-F4_1/d-32/m-origin1'
    SGbase(Y,SG)%PG         ='Oh'
    SGbase(Y,SG)%ncentering =4
    SGbase(Y,SG)%centering(1:3,1) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%centering(1:3,2) = (/           0.0d0,     1.0d0/2.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%centering(1:3,3) = (/     1.0d0/2.0d0,           0.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%centering(1:3,4) = (/     1.0d0/2.0d0,     1.0d0/2.0d0,           0.0d0 /)
    SGbase(Y,SG)%PGNumOper  =48
    SGbase(Y,SG)%symmorphic =0
    SGbase(Y,SG)%nsy(1:3,1) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,2) = (/           0.0d0,     1.0d0/2.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,3) = (/     1.0d0/2.0d0,     1.0d0/2.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,4) = (/     1.0d0/2.0d0,           0.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,5) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,6) = (/     1.0d0/2.0d0,           0.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,7) = (/           0.0d0,     1.0d0/2.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,8) = (/     1.0d0/2.0d0,     1.0d0/2.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,9) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,10) = (/     1.0d0/2.0d0,     1.0d0/2.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,11) = (/     1.0d0/2.0d0,           0.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,12) = (/           0.0d0,     1.0d0/2.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,13) = (/     3.0d0/4.0d0,     1.0d0/4.0d0,     3.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,14) = (/     1.0d0/4.0d0,     1.0d0/4.0d0,     1.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,15) = (/     1.0d0/4.0d0,     3.0d0/4.0d0,     3.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,16) = (/     3.0d0/4.0d0,     3.0d0/4.0d0,     1.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,17) = (/     3.0d0/4.0d0,     1.0d0/4.0d0,     3.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,18) = (/     3.0d0/4.0d0,     3.0d0/4.0d0,     1.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,19) = (/     1.0d0/4.0d0,     1.0d0/4.0d0,     1.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,20) = (/     1.0d0/4.0d0,     3.0d0/4.0d0,     3.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,21) = (/     3.0d0/4.0d0,     1.0d0/4.0d0,     3.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,22) = (/     1.0d0/4.0d0,     3.0d0/4.0d0,     3.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,23) = (/     3.0d0/4.0d0,     3.0d0/4.0d0,     1.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,24) = (/     1.0d0/4.0d0,     1.0d0/4.0d0,     1.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,25) = (/     1.0d0/4.0d0,     1.0d0/4.0d0,     1.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,26) = (/     1.0d0/4.0d0,     3.0d0/4.0d0,     3.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,27) = (/     3.0d0/4.0d0,     3.0d0/4.0d0,     1.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,28) = (/     3.0d0/4.0d0,     1.0d0/4.0d0,     3.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,29) = (/     1.0d0/4.0d0,     1.0d0/4.0d0,     1.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,30) = (/     3.0d0/4.0d0,     1.0d0/4.0d0,     3.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,31) = (/     1.0d0/4.0d0,     3.0d0/4.0d0,     3.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,32) = (/     3.0d0/4.0d0,     3.0d0/4.0d0,     1.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,33) = (/     1.0d0/4.0d0,     1.0d0/4.0d0,     1.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,34) = (/     3.0d0/4.0d0,     3.0d0/4.0d0,     1.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,35) = (/     3.0d0/4.0d0,     1.0d0/4.0d0,     3.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,36) = (/     1.0d0/4.0d0,     3.0d0/4.0d0,     3.0d0/4.0d0 /)
    SGbase(Y,SG)%nsy(1:3,37) = (/     1.0d0/2.0d0,           0.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,38) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,39) = (/           0.0d0,     1.0d0/2.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,40) = (/     1.0d0/2.0d0,     1.0d0/2.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,41) = (/     1.0d0/2.0d0,           0.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,42) = (/     1.0d0/2.0d0,     1.0d0/2.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,43) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,44) = (/           0.0d0,     1.0d0/2.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,45) = (/     1.0d0/2.0d0,           0.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,46) = (/           0.0d0,     1.0d0/2.0d0,     1.0d0/2.0d0 /)
    SGbase(Y,SG)%nsy(1:3,47) = (/     1.0d0/2.0d0,     1.0d0/2.0d0,           0.0d0 /)
    SGbase(Y,SG)%nsy(1:3,48) = (/           0.0d0,           0.0d0,           0.0d0 /)
    SGbase(Y,SG)%reorder_IT_op_sequence = 1
    SGbase(Y,SG)%identity_map(1:48) = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48/)
    SGbase(Y,SG)%operation_map(1:48) = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,18,19,22,24,15,16,17,20,21,23,25,26,27,28,29,30,31,32,33,34,35,36,37,38,42,43,46,48,39,40,41,44,45,47/)
    SGbase(Y,SG)%nop        =192
    SGbase(Y,SG)%assigned   =.TRUE.


  end subroutine init_static_sg_database
end module sgsym
