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

!* Use the concept of the group of k

module extra
  use commod
  implicit none

  character(len=*),parameter :: ineqf='ineq.vasp'
  character(len=*),parameter :: ppf='p.vasp'
  character(len=*),parameter :: scpf='supercell.vasp'
  character(len=*),parameter :: cdlistname='list-of-fm-cd-ids'
  character(len=*),parameter :: fdlistname='list-of-fm-fd-ids'

  integer,parameter :: strlen=4

  integer,parameter :: controlu  =11
  integer,parameter :: iu        =12
  integer,parameter :: pgu       =13
  integer,parameter :: sgu       =14
  integer,parameter :: dsgu      =15
  integer,parameter :: vaspu     =16
  integer,parameter :: fmfu      =17
  integer,parameter :: partialfu =18
  integer,parameter :: pcontrolu =19
  integer,parameter :: debugu    =20

  character(len=*),parameter :: controlf='fm-forces.par'

  character(len=*),parameter :: fm_forces='fm-forces.dat'

  character(len=*),parameter :: proposed_controlf='proposed-fm-forces.par'
  character(len=*),parameter :: partial_forces='partial-forces.dat'

  real(double),parameter :: distance_thres=1.0d-8

  real(double) :: fm_delta
  real(double) :: xyz_delta

  type map_type
    integer :: fmdispn
    integer :: kop
  end type map_type

  character(len=*),parameter :: debugfile='DEBUG_FILE'
  LOGICAL :: debug_si

contains

  subroutine report_all_groups_of_k(ns0,karr,sg,pposcar)
    integer :: ns0,i,j,abs_op
    type(onekgroup) :: karr(ns0)
    type(onesg) :: sg
    type(supercell) :: pposcar
    real(double) :: neta(4),mat(3,3),tracem,determ,rotm(3,3)

    logical :: lprint=.false.

    write(*,*)
    write(*,*) 'Printing (in matrix form) group of k information (nature of each operation):'
    do i = 1, ns0
      write(*,*)
      write(*,*) '----------------------'
      write(*,*) 'ineq atom= '//trim(N2str(i))//', size of the group of k = '//trim(N2str(karr(i)%n))
      write(*,*)
      do j = 1, karr(i)%n

        abs_op = karr(i)%op(j)

        if(lprint) then
          write(*,*)
          write(*,'(A,I5,A,I0,I5,A,I0)') 'j/tot(site_symm),abs_op/tot(SG)=',j,'/',karr(i)%n,abs_op,'/',sg%nop
        endif
        mat(1:3,1:3) = sg%op(2:4,2:4,abs_op)

        determ = det3(mat(1,1))
        tracem = trace3(mat(1,1))

        karr(i)%det(j) = determ
        karr(i)%trace(j) = tracem

        call report_operation_type(nint(determ),nint(tracem))

        rotm = BAinvB(pposcar%a(1,1),mat(1,1))

        if(lprint) then
          write(*,*) 'rotm='
          write(*,'(3F10.4)') rotm(1,1:3)
          write(*,'(3F10.4)') rotm(2,1:3)
          write(*,'(3F10.4)') rotm(3,1:3)
        endif

        if(determ < zero) then
          rotm(1:3,1:3) = -one*rotm(1:3,1:3)
        endif

        neta = rotm_2_neta(rotm(1,1))

        write(*,'(A,4F12.5)') 'neta = ',neta(1:4)
      enddo
    enddo
  end subroutine report_all_groups_of_k

  subroutine equiv_classes_of_group_of_k(ns0,karr,ha,sg,site_symm,sudoku,op_label,indarr,ig,NC,pg_index)
    integer :: ns0
    integer :: ha
    integer :: NC,pg_index
    type(onekgroup) :: karr(ns0)
    real(double) :: site_symm(3,3,48)
    integer :: sudoku(48,48)
    integer :: indarr(48)
    integer :: ig(48+1)
    integer :: op_label(48)
    type(onesg) :: sg

    integer :: i,sg_op,j,k
    integer :: nop
    logical :: found
    real(double) :: targetsum,sumj,mat1(3,3),mat2(3,3),mat3(3,3),zdiff(3,3),prod(3,3)
    real(double) :: sortv(48)

    real(double) :: RG,determ,tracem
    integer :: pg_id,kop,classsize,idet,itrace

    write(*,'(A)') '/// Call equiv_classes_of_group_of_k... (note that the arrays in the argument may be used on the fly hence we have to call more than one time)///'

    nop = karr(ha)%n
    if(nop > 48) then
      write(*,'(A)') 'err: nop is too large.'
      stop 1
    endif

    do i = 1, nop
      sg_op = karr(ha)%op(i)
      site_symm(1:3,1:3,i) = sg%op(2:4,2:4,sg_op)
    enddo

    mat1(1:3,1:3) = zero
    do i = 1, 3
      mat1(i,i) = one
    enddo
    zdiff(1:3,1:3) = site_symm(1:3,1:3,1) - mat1(1:3,1:3)
    if(matmagn(3,zdiff) > 1.0d-5) then
      write(*,*) 'we want the first element to be identity.'
      write(*,'(A)') 'err: first element is not an identity operation'
      stop 1
    endif

    do i = 1, nop
      do j = 1, nop
        call matmatn(3,site_symm(1,1,i),site_symm(1,1,j),prod(1,1))
        found = .false.

        do k = 1, nop
          zdiff(1:3,1:3) = prod(1:3,1:3) - site_symm(1:3,1:3,k)
          if(matmagn(3,zdiff(1,1)) < 1.0d-5) then
            found = .true.
            sudoku(i,j) = k
          endif
        enddo
        if (.not. found) then
          write(*,'(A)') 'err: group operation is not closed'
          stop 1
        endif
      enddo
    enddo

    write(*,*) 'sudoku information'
      do i = 1, nop
      write(*,'(100I4)') sudoku(i,1:nop)
    enddo

    targetsum = nop*(nop+one)/two

    do i = 1, nop
      sumj = zero
      do j = 1, nop
        sumj = sumj + sudoku(i,j)
      enddo
      if(abs(sumj - targetsum) > 1.0d-5) then
        write(*,'(A)') 'err: sum failed 1'
        stop 1
      endif
    enddo

    do i = 1, nop
      sumj = zero
      do j = 1, nop
        sumj = sumj + sudoku(j,i)
      enddo
      if(abs(sumj - targetsum) > 1.0d-5) then
        write(*,'(A)') 'err: sum failed 2'
        stop 1
      endif
    enddo

    op_label(1:nop) = 0
    do i = 1, nop
      if(op_label(i) /= 0) then
        cycle
      endif

      mat1(1:3,1:3) = site_symm(1:3,1:3,i)

      do j = 1, nop

        mat2(1:3,1:3) = site_symm(1:3,1:3,j)

        mat3 = invBAB(mat2(1,1),mat1(1,1))

        found = .false.
        do k = 1, nop
          zdiff(1:3,1:3) = mat3(1:3,1:3) - site_symm(1:3,1:3,k)
          if(matmagn(3,zdiff(1,1)) < 1.0d-5) then
            found = .true.
            if(op_label(k) /= 0) then
              if(op_label(k) /= i) then
                write(*,*) 'i,op_label(k)= ',i,op_label(k)
                write(*,'(A)') 'err: inconsistency.'
                stop 1
              endif
            else
              op_label(k) = i
            endif
          endif
        enddo
        if( .not. found) then
          write(*,'(A)') 'err: not found.'
          stop 1
        endif
      enddo
    enddo

    do i = 1, nop
      sortv(i) = op_label(i)
      indarr(i) = i
    enddo
    do i = 1, nop

    enddo

    call dbl_sort(nop,sortv(1),indarr(1),1)

    do i = 1, nop

    enddo

    RG = sortv(1)
    NC = 1
    IG(1) = 1
    do k = 2, nop
      if(abs(sortv(k)-RG) > 1.0d-4) then
        RG = sortv(k)
        NC = NC + 1
        IG(NC) = k
      endif
    enddo
    ig(NC+1) = nop+1

    write(*,*)
    write(*,'(A)') 'Equivalence classes of the group of k:'

    pg_id = NC*PG_UUID_L5
    do i = 1, NC
      kop = indarr(ig(i))
      determ = karr(ha)%det(kop)
      tracem = karr(ha)%trace(kop)
      classsize=ig(i+1)-ig(i)
      write(*,'(A,I4,I4,2F8.4)') 'cluster, size of cluster, det, trace=',i,classsize,determ,tracem
      idet = nint(determ)
      itrace = nint(tracem)

      call report_operation_type(nint(determ),nint(tracem))
      pg_id = pg_id + pg_id_map(idet,itrace)
    enddo
    pg_index = get_pg_index(pg_id)
    karr(ha)%kname = trim(PG_NAME(pg_index))
    write(*,*) 'pg_index= '//trim(N2str(pg_index))//', site group (group of k) is '//trim(karr(ha)%kname)

  end subroutine equiv_classes_of_group_of_k

  subroutine find_oper_and_disp(fm_delta,xyz_delta,dir,pposcar,ns0,karr,sg,fm_vec,fm_map)
    real(double) :: fm_delta,xyz_delta
    type(supercell) :: pposcar
    integer :: ns0
    real(double) :: fm_vec(3,6,ns0)
    type(map_type) :: fm_map(3,2,ns0)
    type(onekgroup) :: karr(ns0)
    type(onesg) :: sg
    character(len=*) :: dir

    integer :: ha,i,itrace,idet
    real(double) :: tracem,determ
    integer :: nc
    integer :: kop
    integer :: pg_index
    real(double) :: r(3,3),rot_axis(3),rtp(3),neta(4)
    real(double) :: mat0(3,3),mat1(3,3),mat2(3,3),mat21(3,3),mat210(3,3),rotm(3,3)
    real(double) :: rotm2(3,3)
    real(double) :: basicv(3),dimen(6)
    integer :: kop2
    character(len=DSL) :: absproposed_controlf
    real(double) :: zerov(3),oper(3,3),perp(3)
    integer :: pm,fmdispn,negj,j
    real(double) :: rsum(3),nzvec(3),prodvec(3)
    logical :: done,foundneg
    real(double) :: testmat(3,3),tmpmat(3,3),invmat21(3,3)
    integer :: abs_op,idettrace
    logical :: C2done
    logical :: print_info
    character(len=DSL) :: abscontrolf
    real(double) :: determ2,tracem2
    integer :: idet2,itrace2,idettrace2
    real(double) :: oper2(3,3),tmpmat2(3,3)

    real(double) :: site_symm(3,3,48)
    integer :: sudoku(48,48)
    integer :: class_size(48)
    integer :: indarr(48)
    integer :: ig(48+1)
    real(double) :: cartdir(3),fracdir(3)
    logical :: three_disps_per_plusvol
    real(double) :: norm
    integer :: FDdispnum,CDdispnum
    integer :: det_choice

    print_info = .false.

    fm_vec(1:3,1:6,1:ns0) = 1.0d100

    write(*,*)
    write(*,'(A)') 'Check the groups of k derived from each head atom and suggest fundamental displacement patterns:'

    abscontrolf=trim(dir)//'/'//trim(controlf)
    open(controlu,file=trim(abscontrolf),status='old',action='read')
    read(controlu,*) fm_delta
    read(controlu,*) xyz_delta
    close(controlu)

    absproposed_controlf = trim(dir)//'/'//trim(proposed_controlf)
    open(unit=pcontrolu,file=trim(absproposed_controlf),status='replace')
    write(pcontrolu,'(F12.5,A)') fm_delta, ' # delta for displacement in DFT calculation'
    write(pcontrolu,'(F12.5,A)') xyz_delta, ' # delta in the qphonon, but forces are just scaled proportionally'
    write(pcontrolu,'(A)') trim(N2str(ns0))

    write(*,*)

    do ha = 1, ns0
      write(*,*)
      write(pcontrolu,'(A)') '# fm_vec(1:3,i), i=1,...,6. (0,0,0) for zero displacement. Next 6 lines of  "p q". Lines 1 to 3 (4 to 6) for are plus (minus) displacement. Each line with (p,q) relies on fm_vec(1:3,p) with q-th operation within the group of k'
      write(*,'(A)') 'Head atom= '//trim(N2str(ha))//' (of '//trim(N2str(ns0))//')'

      call equiv_classes_of_group_of_k(ns0,karr,ha,sg,site_symm(1,1,1),sudoku(1,1),class_size(1),indarr(1),ig(1),nc,pg_index)

      karr(ha)%inversion = .false.
      do j = 1, karr(ha)%n
        abs_op = karr(ha)%op(j)
        testmat(1:3,1:3) = sg%op(2:4,2:4,abs_op)
        if(abs(det3(  testmat(1,1))- minusone) < 1.0d-10 .and. abs( tracen(3,testmat(1,1)) - (-3.0d0)) < 1.0d-10) then

          write(*,*) 'There is an inversion operation.'
          karr(ha)%inversion = .true.
        endif
      enddo
      write(*,*)
      write(*,'(A)') 'Use the group of k (one of the 32 point groups) to suggest deduce displacement pattern.'

      zerov(1:3) = (/zero,zero,zero/)

      write(*,*)

      if(pg_index .ge. 1 .and. pg_index .le. 2) then
        write(*,*) 'Case 1 of 4 (PGs 1 to 2):'
        done = .true.
        r(1:3,1) = (/one,zero,zero/)
        fm_vec(1:3,1,ha) = (/one,zero,zero/)
        fm_vec(1:3,2,ha) = (/zero,one,zero/)
        fm_vec(1:3,3,ha) = (/zero,zero,one/)
        pm = 1
        fm_map(1,pm,ha)%fmdispn = 1
        fm_map(1,pm,ha)%kop = 1
        fm_map(2,pm,ha)%fmdispn = 2
        fm_map(2,pm,ha)%kop = 1
        fm_map(3,pm,ha)%fmdispn = 3
        fm_map(3,pm,ha)%kop = 1

      else if(pg_index .ge. 3 .and. pg_index .le. 5) then
        write(*,*) 'Case 2 of 4 (PGs 3 to 5)'
        done = .false.
        do i = 1, nc
          if(done) then
            cycle
          endif
          kop = indarr(ig(i))
          determ = karr(ha)%det(kop)
          tracem = karr(ha)%trace(kop)
          idet = nint(determ)
          itrace = nint(tracem)
          idettrace = idet*itrace

          if(idettrace == -1) then
            done = .true.
            oper(1:3,1:3) = site_symm(1:3,1:3,kop)
            if(idet == -1) then
              write(*,*) 'sigma_h: turn to positive determinant.'
              oper(1:3,1:3) = -oper(1:3,1:3)
            else
              write(*,*) 'C2: already in positive determinant'
            endif
            rotm = BAinvB(pposcar%a(1,1),oper(1,1))
            if(print_info) then
              write(*,*) 'rotm= '
              write(*,*) rotm(1,1:3)
              write(*,*) rotm(2,1:3)
              write(*,*) rotm(3,1:3)
            endif
            neta = rotm_2_neta(rotm(1,1))
            write(*,'(A,4F12.5)') 'neta = ',neta(1:4)
            rot_axis(1:3) = neta(1:3)
            rtp = rthetaphi(rot_axis(1))
            write(*,'(A,3F12.5)') 'rtp(1:3) = ',rtp(1:3)

            mat1 = rotmat_y(rtp(2)*rad2deg)
            mat2 = rotmat_z(rtp(3)*rad2deg)
            call matmatn(3,mat2,mat1,mat21)

            basicv(1:3) = (/zero,-one,one/)
            basicv(1:3) = basicv(1:3)/vecmagn(3,basicv(1))

            call matvecn(3,mat21(1,1),basicv(1),r(1,1))

            call matvecn(3,rotm(1,1),r(1:3,1),r(1:3,2))

            perp(1:3) = (/one,zero,zero/)
            call matvecn(3,mat21(1,1),perp(1),r(1,3))

            write(*,*) 'conventional displacement:'
            write(*,'(A,3F10.5)') 'basicv = ',basicv(1:3)
            write(*,*)
            write(*,*) 'Actual displacements (some are from sym operation):'
            write(*,'(A,3F10.5)') 'r1 = ',r(1:3,1)
            write(*,'(A,3F10.5)') 'r2 = ',r(1:3,2)
            write(*,'(A,3F10.5)') 'r3 = ',r(1:3,3)

            call getlenang(r(1,1),dimen(1))
            write(*,*) 'getlenang:'
            write(*,'(3F15.5)') dimen(1:3)
            write(*,'(3F15.5)') dimen(4:6)

            write(*,'(A,F10.5)') 'det3(r)= ',det3(r(1,1))
            if(abs(   abs(det3(r(1,1))) - one) > 1.0d-7) then
              write(*,'(A)') 'err: wrong determinant.'
              stop 1
            endif

            fm_vec(1:3,1,ha) = r(1:3,1)
            fm_vec(1:3,2,ha) = zerov(1:3)
            fm_vec(1:3,3,ha) = r(1:3,3)

            pm = 1
            fm_map(1,pm,ha)%fmdispn = 1
            fm_map(1,pm,ha)%kop = 1
            fm_map(2,pm,ha)%fmdispn = 1
            fm_map(2,pm,ha)%kop = kop
            fm_map(3,pm,ha)%fmdispn = 3
            fm_map(3,pm,ha)%kop = 1
          endif
        enddo

      else if(pg_index .ge. 6 .and. pg_index .le. 8) then
        write(*,*) 'Case 3 of 4 (PGs 6 to 8):'
        done = .false.
        do i = 1, nc
          if(done) then
            cycle
          endif
          kop = indarr(ig(i))
          determ = karr(ha)%det(kop)
          tracem = karr(ha)%trace(kop)
          idet = nint(determ)
          itrace = nint(tracem)
          idettrace = idet*itrace
          if(idettrace == -1) then

            write(*,*) 'Under pg_index: kop = ',kop
            done = .true.

            oper(1:3,1:3) = site_symm(1:3,1:3,kop)
            if(idet == -1) then
              write(*,*) 'turn to positive determinant.'
              oper(1:3,1:3) = -oper(1:3,1:3)
              write(*,'(A)') 'err: Why do we hit S2 rather than C2 first?'
              stop 1
            endif
            rotm = BAinvB(pposcar%a(1,1),oper(1,1))
            neta = rotm_2_neta(rotm(1,1))
            write(*,'(A,4F12.5)') 'neta = ',neta(1:4)
            rot_axis(1:3) = neta(1:3)
            rtp = rthetaphi(rot_axis(1))
            write(*,'(A,3F12.5)') 'rtp(1:3) = ',rtp(1:3)

            mat1 = rotmat_y(rtp(2)*rad2deg)
            mat2 = rotmat_z(rtp(3)*rad2deg)

            call matmatn(3,mat2,mat1,mat21)

            kop2 = indarr(ig(i)+1)
            write(*,*) 'kop2 is ',kop2
            determ2 = karr(ha)%det(kop2)
            tracem2 = karr(ha)%trace(kop2)
            idet2 = nint(determ2)
            itrace2 = nint(tracem2)
            idettrace2= idet2*itrace2
            if(idettrace2 /= -1) then
              write(*,'(A)') 'err: We must have C2 or S2'
              stop 1
            endif
            oper2(1:3,1:3) = site_symm(1:3,1:3,kop2)

            if(idet2 == -1) then
              oper2(1:3,1:3) = -oper2(1:3,1:3)
            endif

            invmat21 = inv3x3(mat21(1,1))
            call matmat3(invmat21(1,1),oper2(1,1),tmpmat(1,1))
            call matmat3(tmpmat(1,1),mat21(1,1),tmpmat2(1,1))

            rotm = BAinvB(pposcar%a(1,1),tmpmat2(1,1))
            neta = rotm_2_neta(rotm(1,1))
            write(*,'(A,4F12.5)') 'Critical: neta = ',neta(1:4)
            rot_axis(1:3) = neta(1:3)
            rtp = rthetaphi(rot_axis(1))
            write(*,'(A,3F12.5)') 'rtp(1:3) = ',rtp(1:3)
            if(abs(rtp(2) - pi/two) > 1.0d-6) then
              write(*,'(A)') 'angle theta must be 90 degrees since this C2 must be perpendicualr to the first C2'
              write(*,'(A)') 'err: not true.'
              stop 1
            endif

            write(*,'(A,F20.10)') 'To rotate the second C2 to x-axis, we need an angle (deg) of ',rtp(3)*rad2deg

            mat0 = rotmat_z(rtp(3)*rad2deg)

            basicv(1:3) = (/one,one,one/)
            basicv(1:3) = basicv(1:3)/vecmagn(3,basicv(1))

            call matmatn(3,mat21(1,1),mat0(1,1),mat210(1,1))
            call matvecn(3,mat210(1,1),basicv(1),r(1,1))

            rotm = BAinvB(pposcar%a(1,1),site_symm(1,1,kop))
            call matvecn(3,rotm(1,1),r(1,1),r(1,2))

            rotm = BAinvB(pposcar%a(1,1),site_symm(1,1,kop2))
            call matvecn(3,rotm(1,1),r(1:3,1),r(1:3,3))

            write(*,*) 'conventional displacement:'
            write(*,'(A,3F10.5)') 'basicv = ',basicv(1:3)
            write(*,*)
            write(*,*) 'Actual displacements (some are from sym operation):'
            write(*,'(A,3F10.5)') 'r1 = ',r(1:3,1)
            write(*,'(A,3F10.5)') 'r2 = ',r(1:3,2)
            write(*,'(A,3F10.5)') 'r3 = ',r(1:3,3)

            call getlenang(r(1,1),dimen(1))
            write(*,*) 'getlenang:'
            write(*,'(3F15.5)') dimen(1:3)
            write(*,'(3F15.5)') dimen(4:6)

            write(*,'(A,F10.5)') 'det3(r)= ',det3(r)
            if(abs(   abs(det3(r)) - 4.0d0/sqrt(27.0d0)) > 1.0d-7) then
              write(*,'(A)') 'err: wrong determinant.'
              stop 1
            endif

            fm_vec(1:3,1,ha) = r(1:3,1)
            fm_vec(1:3,2,ha) = zerov(1:3)
            fm_vec(1:3,3,ha) = zerov(1:3)
            pm = 1
            fm_map(1,pm,ha)%fmdispn = 1
            fm_map(1,pm,ha)%kop = 1
            fm_map(2,pm,ha)%fmdispn = 1
            fm_map(2,pm,ha)%kop = kop
            fm_map(3,pm,ha)%fmdispn = 1
            fm_map(3,pm,ha)%kop = kop2
          endif
        enddo

      else if(pg_index .ge. 9 .and. pg_index .le. 32)  then
        write(*,*) 'Case 4 (PGs 9 to 32):'
        write(*,'(/A)') 'Always find C3 (or 3-bar) first before C4 (or 4-bar)'

        done = .false.

        do i = 1, nc
          if(done) then
            cycle
          endif

          kop = indarr(ig(i))
          determ = karr(ha)%det(kop)
          tracem = karr(ha)%trace(kop)

          idet = nint(determ)
          itrace = nint(tracem)
          idettrace = idet*itrace

          if(idettrace == 1 .or. idettrace == 0) then

            write(*,'(/A,I4,A,I4,A,I4)') 'Enter idettrace loop: idet=',idet,', itrace=',itrace,', idettrace=',idettrace
            if(idet == 1 .and. itrace == 1) then
              write(*,'(A)') 'Operation is 4'
            else if(idet == -1 .and. itrace == -1) then
              write(*,'(A)') 'Operation is 4-bar'
            else if(idet == 1 .and. itrace == 0) then
              write(*,'(A)') 'Operation is 3'
            else if(idet == -1 .and. itrace == 0) then
              write(*,'(A)') 'Operation is 3-bar'
            endif
            write(*,'(/A,I4)') 'kop=',kop

            done = .true.

            oper(1:3,1:3) = site_symm(1:3,1:3,kop)
            if(idet == -1) then
              write(*,*) 'Negative determinant, so we turn it to become positive'
              oper(1:3,1:3) = -oper(1:3,1:3)
            endif

            rotm = BAinvB(pposcar%a(1,1),oper(1,1))

            neta = rotm_2_neta(rotm(1,1))
            write(*,'(A,4F12.5)') 'neta = ',neta(1:4)

            rot_axis(1:3) = neta(1:3)
            rtp = rthetaphi(rot_axis(1))
            write(*,'(A,3F12.5)') 'rtp(1:3) = ',rtp(1:3)

            mat1 = rotmat_y(rtp(2)*rad2deg)
            mat2 = rotmat_z(rtp(3)*rad2deg)
            call matmatn(3,mat2,mat1,mat21)

            write(*,'(//A)') 'Calling C2test...'
            call C2test(ha,kop,ns0,karr,sg,pposcar,basicv(1),C2done)

            if(.not. C2done) then
              write(*,*) 'No C2 can map d to -d.'

              basicv(1:3) = (/one,one,one/)
            endif

            write(*,*)
            write(*,'(A,3F12.5)') 'Canonical basicv used is ',basicv(1:3)

            basicv(1:3) = basicv(1:3)/vecmagn(3,basicv(1))

            call matvecn(3,mat21(1,1),basicv(1),r(1,1))

            call matvecn(3,rotm(1,1),r(1:3,1),r(1:3,2))

            kop2 = sudoku(kop,kop)
            rotm2 = BAinvB(pposcar%a(1,1),site_symm(1,1,kop2))
            call matvecn(3,rotm2(1,1),r(1:3,1),r(1:3,3))

            write(*,*)
            write(*,'(A,3F10.5)') 'r1 = ',r(1:3,1)
            write(*,'(A,3F10.5)') 'r2 = ',r(1:3,2)
            write(*,'(A,3F10.5)') 'r3 = ',r(1:3,3)

            call getlenang(r(1,1),dimen(1))
            write(*,*) 'getlenang:'
            write(*,'(3F15.5)') dimen(1:3)
            write(*,'(3F15.5)') dimen(4:6)

            write(*,'(A,F10.5)') 'det3(r)= ',det3(r)

            if(pg_index .ge. 9. .and. pg_index .le. 15) then
              if(abs(   abs(det3(r)) - 4.0d0/sqrt(27.0d0)) > 1.0d-7) then
                write(*,'(A)') 'err: wrong determinant.'
                stop 1
              endif
            endif

            if(pg_index .ge. 16. .and. pg_index .le. 32) then
              if(abs(   abs(det3(r)) - one) > 1.0d-7) then
                write(*,'(A)') 'err: wrong determinant.'
                stop 1
              endif
            endif
            fm_vec(1:3,1,ha) = r(1:3,1)
            fm_vec(1:3,2,ha) = zerov(1:3)
            fm_vec(1:3,3,ha) = zerov(1:3)
            pm = 1
            fm_map(1,pm,ha)%fmdispn = 1
            fm_map(1,pm,ha)%kop = 1
            fm_map(2,pm,ha)%fmdispn = 1
            fm_map(2,pm,ha)%kop = kop
            fm_map(3,pm,ha)%fmdispn = 1
            fm_map(3,pm,ha)%kop = kop2

            if(debug_si) then

              write(*,'(///A///)') 'WARNING WARNING WARNING WARNING: do si test.'
              write(*,'(A,3F20.12)') 'si_debug is turned on: r1 was: ',r(1:3,1)

              det_choice = 2

              write(*,*) 'det_choice = ',det_choice
              write(*,*) 'det = ',10**(-det_choice)
              if(det_choice == 0) then
                r(1:3,1) = (/  1.0000000000  ,     0.0000000000   ,    0.0000000000/)
              else if(det_choice == 1) then
                r(1:3,1) = (/  0.727700616604407d0,       0.485000934325680d0, 0.485000934325680d0/)
              else if(det_choice == 2) then
                r(1:3,1) = (/  0.626940602185246d0,      0.550883599924522d0, 0.550883599924522d0 /)
              else if(det_choice == 3) then
                r(1:3,1) = (/  0.593259410820008d0,       0.569228983570538d0, 0.569228983570538d0 /)
              else if(det_choice == 4) then
                r(1:3,1) = (/  0.582404777741018d0,     0.574806347766113d0, 0.574806347766113d0 /)
              else if(det_choice == 5) then
                r(1:3,1) = (/  0.578951033892801d0,     0.576548220166560d0, 0.576548220166560d0 /)
              else if(det_choice == 6) then
                r(1:3,1) = (/  0.577856715250969d0, 0.577096879492239d0, 0.577096879492239d0 /)
              else if(det_choice == 7) then
                r(1:3,1) = (/  0.57751044550768d0,  0.577270164363931d0, 0.577270164363931d0   /)
              else if(det_choice == 8) then
                r(1:3,1) = (/  0.57740092379098d0,  0.577324940222282d0, 0.577324940222282d0   /)
              else
                write(*,*) 'err: wrong det_choice = ',det_choice
                stop 1
              endif

              three_disps_per_plusvol = .false.
              if(three_disps_per_plusvol) then

                r(1:3,2) = (/  r(3,1),  r(1,1),    r(2,1)/)
                r(1:3,3) = (/  r(2,1),  r(3,1) ,   r(1,1)/)

                fm_vec(1:3,1,ha) = r(1:3,1)
                fm_vec(1:3,2,ha) = r(1:3,2)
                fm_vec(1:3,3,ha) = r(1:3,3)
                pm = 1

                fm_map(1,pm,ha)%fmdispn = 1
                fm_map(1,pm,ha)%kop = 1
                fm_map(2,pm,ha)%fmdispn = 2
                fm_map(2,pm,ha)%kop = 1
                fm_map(3,pm,ha)%fmdispn = 3
                fm_map(3,pm,ha)%kop = 1
              else

                fm_vec(1:3,1,ha) = r(1:3,1)
              endif
            endif
          endif
        enddo
      else
        write(*,'(A)') 'err: cannot be more than 32.'
        stop 1
      endif

      if(.not. done) then
        write(*,'(A)') 'err: cannot find form a nonzero volume for this site (or atom)?'
        stop 1
      endif

      write(*,*)
      write(*,'(A)') 'Entering the negative displacement part'
      write(*,*)

      do i = 1, 3
        write(*,*) 'Component = ',i
        fmdispn = fm_map(i,1,ha)%fmdispn
        kop = fm_map(i,1,ha)%kop
        nzvec(1:3) = fm_vec(1:3,fmdispn,ha)
        if(vecmag3(nzvec(1)) < 1.0d-8) then
          write(*,'(A)') 'err: must be nonzero vector.'
          stop 1
        endif
        foundneg = .false.
        negj = -1
        do j = 1, karr(ha)%n
          if(foundneg) then
            cycle
          endif

          rotm = BAinvB(pposcar%a(1,1),site_symm(1,1,j))
          call matvecn(3,rotm(1,1),nzvec(1),prodvec(1))

          rsum(1:3) = nzvec(1:3) + prodvec(1:3)
          if(vecmag3(rsum(1)) < 1.0d-8) then

            write(*,*) 'negj: rotm = '
            write(*,'(3F18.9)') rotm(1,1:3)
            write(*,'(3F18.9)') rotm(2,1:3)
            write(*,'(3F18.9)') rotm(3,1:3)
            foundneg = .true.
            negj = j
          endif
        enddo

        if(foundneg) then
          fm_vec(1:3,i+3,ha) = zerov(1:3)
          fm_map(i,2,ha)%fmdispn = fmdispn
          fm_map(i,2,ha)%kop = sudoku(kop,negj)
        endif
      enddo

      do i = 1, 3
        if(vecmagn(3,fm_vec(1:3,i+3,ha)) > 1.0d0) then
          fm_vec(1:3,i+3,ha) = -fm_vec(1:3,i,ha)
          fm_map(i,2,ha)%fmdispn = fm_map(i,1,ha)%fmdispn + 3
          fm_map(i,2,ha)%kop = fm_map(i,1,ha)%kop
        endif
      enddo

      write(pcontrolu,*) ha
      do i = 1, 6
        write(*,*)
        write(*,*) 'Cartesian to fractional coordinates: i = ',i

        cartdir(1:3) = fm_vec(1:3,i,ha)
        call abs2frac(pposcar%a,cartdir(1),fracdir(1))
        write(*,'(A,3F20.10)') 'Not normalized: fracdir(1:3) = ',fracdir(1:3)
        if(vecmag3(fracdir(1)) > zero) then
          norm = vecmag3(fracdir(1))
          write(*,'(A,3F20.10)') 'normalized: fracdir(1:3) =',fracdir(1:3)/norm
        else
          write(*,*) 'zero vector: skipped.'
        endif
        write(pcontrolu,'(3F20.10)') fm_vec(1:3,i,ha)
        write(*,'(A,3F20.10)') 'fm_vec = ',fm_vec(1:3,i,ha)
      enddo
      do pm = 1, 2
        do i = 1, 3
          write(pcontrolu,'(2I5)') fm_map(i,pm,ha)%fmdispn, fm_map(i,pm,ha)%kop
        enddo
      enddo

      FDdispnum = 0
      do i = 1, 3
        if(vecmagn(3,fm_vec(1:3,i,ha)) > 1.0d-5) then
          FDdispnum = FDdispnum + 1
        endif
      enddo
      CDdispnum = 0
      do i = 1, 6
        if(vecmagn(3,fm_vec(1:3,i,ha)) > 1.0d-5) then
          CDdispnum = CDdispnum + 1
        endif
      enddo
      write(*,'(A)') 'ineq atom '//trim(N2str(ha))//': group-of-k= '//trim(karr(ha)%kname)//', FDdispnum= '//trim(N2str(FDdispnum))//', CDdispnum= '//trim(N2str(CDdispnum))

    enddo

    close(pcontrolu)

  end subroutine find_oper_and_disp

  subroutine integer_translation_vecs(ns0,karr,sg,pposcar,m01)
    integer :: ns0,i,j,op
    integer :: m01(ns0)
    type(onekgroup) :: karr(ns0)
    type(onesg) :: sg
    type(supercell) :: pposcar
    integer :: ind1,nk,int1,int2,int3
    real(double) :: dev,fracA(3),fracB(3),avec(4),bvec(4),u(4,4),matv(4),rtrans(4)

    do i = 1, ns0
      ind1 = m01(i)
      fracA = pposcar%at(ind1)%f
      nk = karr(i)%n
      do j = 1, nk
        fracB = fracA
        avec(1:4) = (/one,fracA(1:3)/)
        bvec(1:4) = (/one,fracB(1:3)/)
        op = karr(i)%op(j)
        u(1:4,1:4) = sg%op(1:4,1:4,op)
        call matvecn(4,u(1,1),avec,matv)
        rtrans(1:4) = bvec(1:4) - matv(1:4)

        int1 = nint(rtrans(2))
        int2 = nint(rtrans(3))
        int3 = nint(rtrans(4))
        dev = abs(int1 - rtrans(2)) + abs(int2 - rtrans(3)) + abs(int3 - rtrans(4))
        if(dev > 1d-8) then
          write(*,*) 'dev = ',dev
          write(*,'(A)') 'err: Non integer solution.'
          stop 1
        endif

        karr(i)%trans(1:3,j) = (/int1,int2,int3/)
      enddo
    enddo
  end subroutine integer_translation_vecs

  subroutine C2test(ha,kop,ns0,karr,sg,pposcar,basicv,C2done)
    integer :: ha,kop,ns0
    type(onekgroup) :: karr(ns0)
    type(onesg) :: sg
    type(supercell) :: pposcar
    real(double) :: basicv(3)
    logical :: C2done
    integer :: i

    real(double) :: determ1,tracem1
    integer :: idet1,itrace1
    real(double) :: determ2,tracem2
    integer :: idet2,itrace2

    real(double) :: rotm(3,3),rotm2(3,3)
    real(double) :: rtp(3),vecr(3),site_symm(3,3),oper(3,3),neta(4)
    real(double) :: rneta(3),rotm3(3,3)
    integer :: abs_op
    real(double) :: Rz(3,3),Ry(3,3),neta2(4),diff
    real(double) :: newneta(4)
    real(double) :: nx,ny,nz,Aquad,Bquad,Cquad,DetQuad
    real(double) :: yconst,yminus,yplus,xminus,xplus,crossp(3),crosspm
    real(double) :: cangle,angle,xsq
    real(double) :: Tmat(3,3),invTmat(3,3)

    determ1 = karr(ha)%det(kop)
    tracem1 = karr(ha)%trace(kop)
    idet1 = nint(determ1)
    itrace1 = nint(tracem1)

    abs_op = karr(ha)%op(kop)
    site_symm(1:3,1:3) = sg%op(2:4,2:4,abs_op)
    oper(1:3,1:3) = site_symm(1:3,1:3)
    if(itrace1 < 0) then
      oper(1:3,1:3) = -oper(1:3,1:3)
    endif

    rotm = BAinvB(pposcar%a(1,1),oper(1,1))
    neta = rotm_2_neta(rotm(1,1))

    write(*,'(A,4F12.5)') 'C3: neta = ',neta(1:4)

    vecr(1:3) = neta(1:3)
    rtp = rthetaphi(vecr(1))
    write(*,'(A,3F12.5)') 'rtp(1:3) = ',rtp(1:3)

    Rz(1:3,1:3) = rotmat_z(rtp(3)*rad2deg)
    Ry(1:3,1:3) = rotmat_y(rtp(2)*rad2deg)

    call matmatn(3,Rz(1,1),Ry(1,1),Tmat(1,1))
    invTmat = inv3x3(Tmat(1,1))

    C2done = .false.

    do i = 1, karr(ha)%n
      if(C2done) then
        write(*,*)
        write(*,'(A)') 'C2 is done, Exit immediately'
        exit
      endif
      write(*,*)
      write(*,'(A,I4,A,I4,A)') 'Checking if i=',i,' (of ',karr(ha)%n,') can do the C2 job...'
      if(i == kop) then
        write(*,'(A)') 'Nope. same operation.'
        cycle
      endif

      abs_op = karr(ha)%op(i)
      site_symm(1:3,1:3) = sg%op(2:4,2:4,abs_op)
      determ2 = det3(site_symm(1,1))
      tracem2 = trace3(site_symm(1,1))
      idet2 = nint(determ2)
      itrace2 = nint(tracem2)
      write(*,'(A,I4,A,I4)') 'idet2=',idet2,', itrace2=',itrace2

      if(idet2 == 1 .and. itrace2 == -1) then

        write(*,*) 'A 2-fold operation is found'

        oper(1:3,1:3) = site_symm(1:3,1:3)
        rotm2 = BAinvB(pposcar%a(1,1),oper(1,1))
        neta2 = rotm_2_neta(rotm2(1,1))

        write(*,'(A,4F12.5)') 'neta2(1:4)= ',neta2(1:4)
        write(*,'(A,4F12.5)') 'compare neta(1:4)= ',neta(1:4)

        crossp = crossprod(neta2(1:3),neta(1:3))
        crosspm = vecmag3(crossp(1))
        write(*,'(A,F12.5)') 'crosspm = ',crosspm
        if( abs ( crosspm ) < 1.0d-8) then
          cycle
        endif

        angle = angle_between_vecs3(neta2(1:3),neta(1:3))
        cangle = atan2(one,sqrt2)*rad2deg
        write(*,*) 'angle(deg), cangle = ',angle,cangle
        if(angle < cangle .or. angle > 180.0d0-cangle) then
          write(*,*) 'skip this and test the next operation?'
          write(*,*) 'err: Quadratic equation has no solution.'
          stop 1
        endif

        call matvec3(invTmat(1,1),neta2(1),rneta(1))

        rotm3 = invBAB(Tmat(1,1),rotm2(1,1))

        newneta = rotm_2_neta(rotm3(1,1))
        diff = vecmag3( rneta(1:3) - newneta(1:3))

        if( abs(diff) < 1.0d-8 .or. abs(diff - 2.0d0) < 1.0d-8) then
          write(*,'(A,4F12.5)') 'rneta(1:3) = ',rneta(1:3)
          write(*,'(A,4F12.5)') 'newneta = ',newneta(1:4)
          write(*,*) 'diff',diff
          write(*,*) 'Good. diff should be either 0.0 or 2.0'
        else
           write(*,'(A)') 'err: too large difference.'
           stop 1
        endif

        vecr(1:3) = newneta(1:3)
        nx = vecr(1)
        ny = vecr(2)
        nz = vecr(3)
        write(*,'(A,3F20.10)') 'nx,ny,nz=',nx,ny,nz
        Aquad = nx**2  + ny**2
        Bquad = 2*ny*nz/sqrt3
        Cquad = nz**2/3.0d0 - 2.0d0*nx**2/3.0d0

        write(*,'(A,3F20.10)') 'Aquad,Bquad,Cquad=',Aquad,Bquad,Cquad
        Detquad = Bquad**2 - four*Aquad*Cquad
        if(Detquad < zero) then
          write(*,*) 'Detquad = ',Detquad
          write(*,'(A)') 'err: cannot take sqrt here.'
          stop 1
        endif
        if(abs(Aquad) < 1.0d-10) then
          write(*,'(A)') 'err: Aquad is too small.'
          stop 1
        endif

        C2done = .true.

        yminus = (-Bquad-sqrt(Detquad))/(two*Aquad)
        yplus = (-Bquad+sqrt(Detquad))/(two*Aquad)
        write(*,*) 'Aquad = ',Aquad
        write(*,*) 'solving quadratic equation: yminus,yplus=',yminus,yplus

        if(abs(nx) > 1.0d-8) then
          xminus = (-nz/sqrt3 - ny*yminus)/nx
          xplus = (-nz/sqrt3 - ny*yplus)/nx
        else
          yconst = ny*yminus+ nz/sqrt3
          if(abs(yconst) > 1.0d-8) then
            write(*,'(A)') 'err: impossible.'
            stop 1
          else

            xsq = two/three  - yminus**2
            if(xsq < zero) then
              xsq = zero
            endif
            xminus = -sqrt(xsq)
          endif

          yconst = ny*yplus+nz/sqrt3
          if(abs(yconst) > 1.0d-8) then
            write(*,'(A)') 'err: impossible.'
            stop 1
          else

            xsq = two/three  - yplus**2
            if(xsq < zero) then
              xsq = zero
            endif
            xplus = sqrt(xsq)
          endif
        endif

        write(*,'(A,2F20.10)') 'xminus,yminus=',xminus,yminus
        write(*,'(A,2F20.10)') 'xplus,yplus=',xplus,yplus

        yconst = xminus**2 + yminus**2 + third - one
        if(abs(yconst) > 1.0d-8) then
          write(*,'(A)') 'err: err1'
          stop 1
        endif
        yconst = nx*xminus + ny*yminus + nz/sqrt3
        if(abs(yconst) > 1.0d-8) then
          write(*,'(A)') 'err: err2'
          stop 1
        endif
        yconst = xplus**2 + yplus**2 + third - one
        if(abs(yconst) > 1.0d-8) then
          write(*,'(A)') 'err: err3'
          stop 1
        endif
        yconst = nx*xplus + ny*yplus + nz/sqrt3
        if(abs(yconst) > 1.0d-8) then
          write(*,'(A)') 'err: err4'
          stop 1
        endif

        basicv(1:3) = (/xminus,yminus,one/sqrt3/)
        write(*,'(A)') 'We have found a candidate. C2 test is done here.'
      else
        write(*,'(A)') 'Nope. not a C2.'
      endif
    enddo
  end subroutine C2test

  subroutine check_ndir(ns0,fm_vec,fm_map,ndir,sg,karr,pposcar)
    integer :: ns0
    real(double) :: fm_vec(3,6,ns0)
    type(map_type) :: fm_map(3,2,ns0)
    real(double) :: ndir(3,3,2,ns0)
    integer :: i
    integer :: dispn
    type(onekgroup) :: karr(ns0)
    type(onesg) :: sg
    type(supercell) :: pposcar

    integer :: kop
    integer :: fmdispn
    integer :: abs_op
    real(double) :: u(4,4)
    real(double) :: sgrotm(3,3),rotm(3,3),lattmat(3,3)
    real(double) :: dimen(6),box(3,3),invbox(3,3)
    integer :: pm
    real(double) :: norm

    write(*,'(/A)') 'Scale ndir(3,3,2) and make sure that we do not hit a singular matrix.'
    do i = 1, ns0
      write(*,*)
      write(*,*) 'Head atom= ',i
      do pm = 1, 2

        write(*,*)
        write(*,*) '---------------------------'
        write(*,*) 'pm= ',pm
        do dispn = 1, 3

          fmdispn = fm_map(dispn,pm,i)%fmdispn
          kop = fm_map(dispn,pm,i)%kop

          abs_op = karr(i)%op(kop)
          u(1:4,1:4) = sg%op(1:4,1:4,abs_op)
          u(2:4,1) = u(2:4,1) + karr(i)%trans(1:3,kop)

          lattmat(1:3,1:3) = pposcar%a(1:3,1:3)
          sgrotm(1:3,1:3) = u(2:4,2:4)

          rotm = BAinvB(lattMat(1,1),sgrotm(1,1))

          call matvecn(3,rotm(1,1),fm_vec(1,fmdispn,i),ndir(1,dispn,pm,i))

          norm = vecmag3(  ndir(1,dispn,pm,i) )
          if(norm < 1.0d-8) then
            write(*,'(A)') 'err: Too small norm after rotation?'
            stop 1
          endif
          ndir(1:3,dispn,pm,i) = fm_delta*ndir(1:3,dispn,pm,i)/norm
        enddo
        box(1:3,1:3) = ndir(1:3,1:3,pm,i)
        write(*,*) 'box = '
        write(*,'(A,3F20.10)') 'first  ndir vector:', box(1:3,1)
        write(*,'(A,3F20.10)') 'second ndir vector:', box(1:3,2)
        write(*,'(A,3F20.10)') 'third  ndir vector:', box(1:3,3)
        call getlenang(box,dimen(1))
        write(*,*) 'getlenang:'
        write(*,'(3F15.5)') dimen(1:3)
        write(*,'(3F15.5)') dimen(4:6)
        write(*,'(A,F20.10)') 'normalized volume is ',det3(box(1,1))/(dimen(1)*dimen(2)*dimen(3))

        invbox = inv3x3(box(1,1))
      enddo
    enddo
  end subroutine check_ndir

  subroutine output_all_struc_for_DFT_cal(ns0,ns1,dir,superc,m01,m12,fm_vec)

    integer :: ns0,ns1
    real(double) :: fm_vec(3,6,ns0)

    type(supercell) :: superc,dispp
    integer :: m01(ns0),m12(ns1)
    integer :: strucindex,ind1,ind2
    character(len=DSL) :: dir,poscaroutfile
    integer :: i,j,k
    real(double) :: norm,absr(3)
    integer :: FDu,CDu
    character(len=DSL) :: FDlist,CDlist

    strucindex = 0

    poscaroutfile = trim(dir)//'/'//'fmdisp-'//trim(num2str(strucindex,strlen))//'.vasp'
    superc%fracabs='frac'
    superc%sg_label = ''
    call supercell_2_vasp(vaspu,trim(poscaroutfile),superc%nsp,superc%zarr(1),superc)

    FDu = newunit()
    FDlist = trim(dir)//'/'//trim(fdlistname)
    write(*,*) 'FDlist is '//trim(FDlist)
    open(unit=FDu,file=trim(FDlist),status='replace')
    write(FDu,'(A)') trim(num2str(strucindex,strlen))

    CDu = newunit()
    CDlist = trim(dir)//'/'//trim(cdlistname)
    write(*,*) 'CDlist is '//trim(cdlist)
    open(unit=CDu,file=trim(CDlist),status='replace')
    write(CDu,'(A)') trim(num2str(strucindex,strlen))

    do i = 1, ns0
      ind1 = m01(i)
      ind2 = m12(ind1)
      do j = 1, 6
        norm = vecmag3(   fm_vec(1:3,j,i)  )
        if(norm > 1.0d-10) then
          strucindex = strucindex + 1
          if(j .le. 3) then
            write(FDu,'(A)') trim(num2str(strucindex,strlen))
          endif
          write(CDu,'(A)') trim(num2str(strucindex,strlen))

          dispp%n = superc%n
          dispp%a(1:3,1:3) = superc%a(1:3,1:3)
          call real2recip(dispp%a,dispp%b)
          dispp%nsp = superc%nsp
          dispp%zarr(1:dispp%nsp) = superc%zarr(1:superc%nsp)
          do k = 1, dispp%n
            dispp%at(k)%f(1:3) = superc%at(k)%f(1:3)
            dispp%at(k)%z = superc%at(k)%z
          enddo
          dispp%fracabs = 'frac'
          dispp%sg_label = ''

          call frac2abs(dispp%a,dispp%at(ind2)%f,absr(1))

          write(*,'(A,6F18.10)') 'absr,disp=',absr(1:3),fm_delta*fm_vec(1:3,j,i)/norm
          absr(1:3) = absr(1:3) + fm_delta*fm_vec(1:3,j,i)/norm
          call abs2frac(dispp%b,absr(1),dispp%at(ind2)%f)

          poscaroutfile = trim(dir)//'/'//'fmdisp-'//trim(num2str(strucindex,strlen))//'.vasp'

          dispp%fracabs='frac'
          call supercell_2_vasp(vaspu,trim(poscaroutfile),dispp%nsp,dispp%zarr(1),dispp)
        endif
      enddo
    enddo
    write(*,*)
    write(*,'(A)') 'The last poscaroutfile is '//trim(poscaroutfile)
    write(*,*)

    close(FDu)
    close(CDu)
  end subroutine output_all_struc_for_DFT_cal

  subroutine read_fm_forces(absfm_forces,fmfu,ns0,ns1,ns2,m01,m12,fm_vec,zf,nzf)
    integer :: fmfu,ns0,ns1,ns2,m01(ns0),m12(ns1)
    real(double) :: fm_vec(3,6,ns0)
    real(double) :: nzf(3,ns2,6,ns1),zf(3,ns2)

    integer :: strucindex,i,j,ii,ind1,ind2
    real(double) :: norm
    character(len=DSL) :: absfm_forces
    integer :: jj,error

    strucindex = 0
    open(unit=fmfu,file=trim(absfm_forces),status='old',action='read',iostat=error)
    if(error /= 0) then
      write(*,'(A)') 'err: Cannot find forces file='//trim(absfm_forces)//'. You may first carry out DFT force cals.'
      stop 1
    endif

    do i = 1, ns2
      read(fmfu,*) zf(1:3,i)
    enddo
    do i = 1, ns0
      ind1 = m01(i)
      ind2 = m12(ind1)
      do jj = 1, 6
        j = jj
        norm = vecmag3(  fm_vec(1:3,j,i ))
        if(norm > 1.0d-10) then
          strucindex = strucindex + 1

          do ii = 1, ns2
            read(fmfu,*) nzf(1:3,ii,j,ind1)
          enddo
        endif
      enddo
    enddo
    close(fmfu)
  end subroutine read_fm_forces

  subroutine forces_trans(ns0,ns1,ns2,m01,fm_map,ndir,sg,karr,pposcar,superc,nzf,targetDeltaMat,new_cartesian_nzf)
    real(double) :: targetDeltaMat(3,3)
    integer :: ns0,ns1,ns2
    integer :: m01(ns0)
    real(double) :: ndir(3,3,2,ns0)
    type(map_type) :: fm_map(3,2,ns0)
    real(double) :: nzf(3,ns2,6,ns1)
    real(double) :: new_cartesian_nzf(3,ns2,3,ns1,2)
    type(onekgroup) :: karr(ns0)
    type(onesg) :: sg
    type(supercell) :: pposcar,superc

    integer :: kop
    integer :: fmdispn
    integer :: op
    real(double) :: u(4,4)
    real(double) :: sgrotm(3,3),rotm(3,3),lattmat(3,3)
    integer :: ja,pm
    real(double) :: frcvec(3),new_frcvec(3)

    integer :: i
    integer :: dispn
    real(double) :: box(3,3),invbox(3,3),rotFmat(3,3),tmpmat(3,3),tmpmat2(3,3)
    integer :: ind1

    integer :: iprint

    integer,allocatable :: map(:)
    integer,allocatable :: inv_map(:)
    real(double),allocatable :: rot_nzf(:,:,:,:,:)
    real(double),allocatable :: final_rot_nzf(:,:,:,:,:)
    real(double) :: signedDeltaMat(3,3)

    allocate(rot_nzf(3,ns2,3,ns1,2))
    allocate(final_rot_nzf(3,ns2,3,ns1,2))

    allocate(map(ns2))
    allocate(inv_map(ns2))

    do i = 1, ns0
      ind1 = m01(i)

      do pm = 1, 2

        do dispn = 1, 3

          fmdispn = fm_map(dispn,pm,i)%fmdispn

          kop = fm_map(dispn,pm,i)%kop

          op = karr(i)%op(kop)
          u(1:4,1:4) = sg%op(1:4,1:4,op)
          u(2:4,1) = u(2:4,1) + karr(i)%trans(1:3,kop)

          lattmat(1:3,1:3) = pposcar%a(1:3,1:3)
          sgrotm(1:3,1:3) = u(2:4,2:4)

          rotm = BAinvB(lattMat(1,1),sgrotm(1,1))

          call get_supercell_atom_map(u(1,1),pposcar,superc,map(1))

          call inversemapping(ns2,map(1),inv_map(1))

          do ja = 1, ns2
            frcvec(1:3) = nzf(1:3,ja,fmdispn,ind1)
            call matvecn(3,rotm(1,1),frcvec(1),new_frcvec(1))
            rot_nzf(1:3,ja,dispn,ind1,pm) =  new_frcvec(1:3)
          enddo

          do ja = 1, ns2

            final_rot_nzf(1:3,ja,dispn,ind1,pm) = rot_nzf(1:3,inv_map(ja),dispn,ind1,pm)
          enddo
        enddo

        box(1:3,1:3) = ndir(1:3,1:3,pm,i)

        invbox = inv3x3(box(1,1))

        do ja = 1, ns2
          rotFmat(1:3,1:3) = final_rot_nzf(1:3,ja,1:3,ind1,pm)

          call matmatn(3,rotFmat(1,1),invbox(1,1),tmpmat(1,1))

          iprint = 1
          if(iprint == 1) then
            write(*,*) '--------------------'
            write(*,*) 'rotFmat is '
            write(*,'(3F20.10)') rotFmat(1,1:3)
            write(*,'(3F20.10)') rotFmat(2,1:3)
            write(*,'(3F20.10)') rotFmat(3,1:3)
            write(*,*) 'box is '
            write(*,'(3F20.10)') box(1,1:3)
            write(*,'(3F20.10)') box(2,1:3)
            write(*,'(3F20.10)') box(3,1:3)
            write(*,*) 'invbox is '
            write(*,'(3F20.10)') invbox(1,1:3)
            write(*,'(3F20.10)') invbox(2,1:3)
            write(*,'(3F20.10)') invbox(3,1:3)
            write(*,*) 'det rotFmat =',det3(rotFmat(1,1))
            write(*,*) 'det invbox =',det3(invbox(1,1))
            write(*,*) 'tmpmat is '
            write(*,'(3F20.10)') tmpmat(1,1:3)
            write(*,'(3F20.10)') tmpmat(2,1:3)
            write(*,'(3F20.10)') tmpmat(3,1:3)
          endif
          if(pm == 1) then
            signedDeltaMat(1:3,1:3) = targetDeltaMat(1:3,1:3)
          else
            signedDeltaMat(1:3,1:3) = -targetDeltaMat(1:3,1:3)
          endif

          call matmatn(3,tmpmat(1,1),signedDeltaMat(1,1),tmpmat2(1,1))
          if(iprint == 1) then
            write(*,*) 'signedDeltaMat is '
            write(*,'(3F20.10)') signedDeltaMat(1,1:3)
            write(*,'(3F20.10)') signedDeltaMat(2,1:3)
            write(*,'(3F20.10)') signedDeltaMat(3,1:3)
            write(*,*) 'tmpmat2 is '
            write(*,'(3F20.10)') tmpmat2(1,1:3)
            write(*,'(3F20.10)') tmpmat2(2,1:3)
            write(*,'(3F20.10)') tmpmat2(3,1:3)
          endif
          new_cartesian_nzf(1:3,ja,1:3,ind1,pm) = tmpmat2(1:3,1:3)
        enddo
      enddo
    enddo

    deallocate(map)
    deallocate(inv_map)
    deallocate(rot_nzf)
    deallocate(final_rot_nzf)
  end subroutine forces_trans

end module extra

program sample
  use extra
  implicit none

  character(len=DSL) :: dir
  character(len=DSL) :: absineqf
  character(len=DSL) :: absppf
  character(len=DSL) :: absscpf

  character(len=DSL) :: abspartial_forces
  character(len=DSL) :: absfm_forces

  character(len=DSL) :: sg_label,filestem,inputformat
  character(len=DSL) :: abscontrolf
  character(len=DSL) :: absfmdispposcardir

  character(len=DSL) :: cmd

  type(onesg) :: sg
  type(onekstar),allocatable :: kstar_arr(:)
  type(onekgroup),allocatable :: karr(:)

  type(supercell) :: ineqp,pposcar,superc

  integer :: nsp,zord(zordmax)
  integer :: ns0,ns1,ns2
  integer :: i
  integer,allocatable :: m01(:)
  integer,allocatable :: m12(:)
  real(double) :: targetDeltaMat(3,3)

  real(double),allocatable :: fm_vec(:,:,:)
  type(map_type),allocatable :: fm_map(:,:,:)

  real(double),allocatable :: ndir(:,:,:,:)

  real(double),allocatable :: nzf(:,:,:,:),zf(:,:)
  real(double),allocatable :: new_cartesian_nzf(:,:,:,:,:)
  real(double) :: tvec(3)
  logical :: fexist
  character(len=DSL) :: qtablepath

  call fargv(1,dir)
  if(trim(dir) == 'help') then
    write(*,'(A)') 'a.out (1) working-dir'
    write(*,'(A)') 'Prepare a control file '//trim(controlf)
    write(*,'(A)') 'Required input files: ineq.vasp, supercell.vasp.'
    write(*,'(A)') 'p.vasp is generated on the fly'
    stop 1
  endif
  write(*,*) '-----fm-forces starts --------------'

  inquire(file=trim(debugfile),exist=fexist)
  if(fexist) then
    write(*,*) trim(debugfile)//' is present. Hence we will do debugging'
    debug_si = .true.

  else
    write(*,*)
    write(*,*)
    write(*,*)
    write(*,*) trim(debugfile)//' is not present. Hence normal execution.'
    debug_si = .false.
  endif

  call init_static_sg_database(qtablepath)

  call PG_NAME_init()

  write(*,*)
  write(*,*)

  write(*,'(A)') '  ssssssssssssssssssssssssssssssss'
  call fargn(1)
  call fargv(1,dir)

  write(*,*) 'ineq atoms are in '//trim(ineqf)
  absineqf=trim(dir)//'/'//trim(ineqf)

  call read_struc(iu,trim(absineqf),filestem,inputformat,ineqp)

  if(trim('x')//trim(ineqp%sg_label) == 'x') then
    write(*,*) 'SG label is not found in ch.vasp.'
    write(*,'(A)') 'err: check the first line of ch.vasp .'
    stop 1
  else
    sg_label = trim(ineqp%sg_label)

  endif

  call read_sg(qtablepath,sg_label,sg,pgu,sgu,dsgu,sgtable,dsgtable)

  allocate(kstar_arr(ineqp%n))
  call create_primitive_cell(sg,ineqp,pposcar,kstar_arr)

  nsp = ineqp%nsp

  zord(1:nsp) = ineqp%zarr(1:nsp)

  absppf=trim(dir)//'/'//trim(ppf)
  write(*,'(A)') 'vasp output file is '//trim(absppf)
  pposcar%fracabs='frac'
  pposcar%sg_label = trim(ineqp%sg_label)
  call supercell_2_vasp(vaspu,trim(absppf),nsp,zord(1),pposcar)
  write(*,'(A)') trim(absppf)//' is produced'

  absscpf=trim(dir)//'/'//trim(scpf)
  write(*,'(A)')
  write(*,'(A)') 'To read supercell poscar:'
  call read_struc(iu,trim(absscpf),filestem,inputformat,superc)

  ns0 = ineqp%n
  allocate(m01(ns0))
  call assign_m01(m01(1),ns0,ineqp,pposcar,distance_thres)

  ns1 = pposcar%n
  ns2 = superc%n
  allocate(m12(ns1))

  tvec(1:3) = (/zero,zero,zero/)
  call assign_m12_tvec(tvec(1),m12(1),ns1,pposcar,superc,distance_thres)

  write(*,*) 'The indices of inequivalent atoms in ch.vasp, p.vasp and supercell.vasp:'
  do i = 1, ns0
    write(*,*) 'i,m01,m12=',i,m01(i),m12(m01(i))
  enddo

  allocate(karr(ineqp%n))
  call find_group_of_k(ns0,ineqp,pposcar,sg,m01(1),karr)
  write(*,*) 'Done with find_group_of_k'

  call report_all_groups_of_k(ns0,karr,sg,pposcar)

  write(*,*)
  write(*,*) 'Done with printing the group of k (in matrix form)'

  allocate(fm_vec(3,6,ns0))
  allocate(fm_map(3,2,ns0))

  write(*,*)

  write(*,'(A)') '/// calling find_oper_and_disp()...'

  call find_oper_and_disp(fm_delta,xyz_delta,dir,pposcar,ns0,karr,sg,fm_vec(1,1,1),fm_map(1,1,1))
  write(*,*) 'Done with minimal number of displacement needed'

  call integer_translation_vecs(ns0,karr,sg,pposcar,m01(1))

  cmd='rm -f ./pg-table ./sg-table ./dsg-table'
  call system(trim(cmd))

  allocate(zf(3,ns2))
  allocate(nzf(3,ns2,6,ns1))

  allocate(new_cartesian_nzf(3,ns2,3,ns1,2))

  new_cartesian_nzf(1:3,1:ns2,1:3,1:ns1,1:2) = zero

  abscontrolf=trim(dir)//'/'//trim(controlf)

  allocate(ndir(3,3,2,ns0))

  call check_ndir(ns0,fm_vec(1,1,1),fm_map(1,1,1),ndir(1,1,1,1),sg,karr,pposcar)

  write(*,*) 'Done with volume business for matrix inversion'
  write(*,*)

  targetDeltaMat(1:3,1:3) = reshape( (/xyz_delta, zero, zero, zero, xyz_delta,zero, zero,zero,xyz_delta/), (/3,3/))

  absfmdispposcardir=trim(dir)//'/fmdisp-poscar-dir'
  call system("rm -rf "//trim(absfmdispposcardir))
  call system("mkdir "//trim(absfmdispposcardir))

  call output_all_struc_for_DFT_cal(ns0,ns1,absfmdispposcardir,superc,m01(1),m12(1),fm_vec(1,1,1))

  write(*,*) 'ns1,ns2=',ns1,ns2

  zf(1:3,1:ns2) = zero
  nzf(1:3,1:ns2,1:6,1:ns1) = zero

  absfm_forces = trim(dir)//'/'//trim(fm_forces)
  write(*,'(A)') 'absfm_forces is '//trim(absfm_forces)
  call read_fm_forces(absfm_forces,fmfu,ns0,ns1,ns2,m01(1),m12(1),fm_vec(1,1,1),zf(1,1),nzf(1,1,1,1))

  write(*,*) 'Call forces_trans...'
  call forces_trans(ns0,ns1,ns2,m01(1),fm_map(1,1,1),ndir(1,1,1,1),sg,karr,pposcar,superc,nzf(1,1,1,1),targetDeltaMat(1,1),new_cartesian_nzf(1,1,1,1,1))

  abspartial_forces = trim(dir)//'/'//trim(partial_forces)
  write(*,*) 'Call write_forces...'
  call write_forces(partialfu,abspartial_forces,ns1,ns2,zf,new_cartesian_nzf(1,1,1,1,1))

  write(*,*)
  write(*,'(A)') 'input file= '//trim(absfm_forces)
  write(*,'(A)') 'output file= '//trim(abspartial_forces)

  write(*,*)
  write(*,'(A)') 'fullf or /opt1/full-forces . sym'

  deallocate(fm_map)
  deallocate(fm_vec)

end program sample
