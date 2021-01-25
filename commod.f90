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

module commod
  use sgsym
  use derived_constants
  implicit none

  interface fargv
    module procedure getarg_str,getarg_int,getarg_double
  end interface

contains

  subroutine getarg_str(n,str)
    integer :: n
    character(len=*) :: str
    call get_command_argument(n,str)
  end subroutine getarg_str

  subroutine getarg_int(n,k)
    integer :: n, k
    character(len=100) :: tmpstr
    call get_command_argument(n,tmpstr)
    read(unit=tmpstr,FMT='(I10)') k
  end subroutine getarg_int

  subroutine getarg_double(n,k)
    integer :: n
    real(double) :: k
    character(len=100) :: tmpstr
    call get_command_argument(n,tmpstr)
    read(unit=tmpstr,FMT='(D22.15)') k
  end subroutine getarg_double


  subroutine inversemapping(n,fdindex,bdindex)
    integer :: n,i,ind
    integer :: fdindex(n)
    integer :: bdindex(n)
    do i = 1, n
      ind = fdindex(i)
      bdindex(ind) = i
    enddo
  end subroutine inversemapping

  function num2str(v,nspaces) result (str)
    integer :: tmpv,i,v,nspaces,digit,ind
    character(len=nspaces) :: str
    tmpv = v
    if(v < 0) then
      write(*,'(A)') 'err: v is negative.'
      stop 1
    endif
    do i = 1, nspaces
      digit = modulo(tmpv,10)
      tmpv = tmpv/10
      ind = nspaces-i+1
      str(ind:ind) = DIGITARR(digit)
    enddo
  end function num2str

  function str2N(str) result(n)
    character(len=*) :: str
    integer :: n
    read(str,'(I10)') n
  end function str2N

  function N2str(i) result (str)
    integer :: i,j,tmpv,intarr(20),nd,digit
    character(len=20) :: str

    if(i < 0) then
      tmpv = -i
    else
      tmpv = I
    endif

    nd = 0
    do
      digit = modulo(tmpv,10)
      nd = nd + 1
      intarr(nd) = digit
      tmpv = tmpv/10
      if(tmpv == 0) exit
    enddo
    str = ' '
    do j = 1, nd
      str(j:j) = digitarr(intarr(nd-j+1))
    enddo
    if(i < 0) then
      nd = nd + 1
      do j = nd,2,-1
        str(j:j) = str(j-1:j-1)
      enddo
      str(1:1) = '-'
    endif
  end function N2str

  subroutine fargn(n)
    integer :: n
    integer :: m

    m = command_argument_count()
    if(n /= m) then
      write(*,'(A,I4,A)') 'we need ',n,' arguments'
      write(*,'(A)') 'err: check number of arguments'
      stop 1
    endif
  end subroutine fargn

  function num_of_strings(tmp)
    character(len=*) :: tmp
    integer :: num_of_strings
    character(len=DSL) :: line
    integer,parameter :: max_num_cols = 100

    character(len=100) :: test_array(max_num_cols)
    integer :: i,io

    line = trim(tmp)
    write(*,'(A)') 'line is '//trim(line)
    do i = 1, max_num_cols

      read(line,*,iostat=io) test_array(1:i)

      if(io /= 0) then
        exit
      endif
    enddo
    num_of_strings = i-1
    if(num_of_strings == max_num_cols) then

    endif
  end function num_of_strings

  function num_of_integers(tmp)
    character(len=*) :: tmp
    integer :: num_of_integers
    character(len=DSL) :: line
    integer,parameter :: max_num_cols = 100
    integer :: test_array(max_num_cols)
    integer :: i,io

    line = trim(tmp)
    write(*,'(A)') 'line is '//trim(line)
    do i = 1, max_num_cols

      read(line,*,iostat=io) test_array(1:i)
      if(io /= 0) then
        exit
      endif
    enddo
    num_of_integers = i-1
    if(num_of_integers == max_num_cols) then
      write(*,*) 'count_int_number= ',num_of_integers
      write(*,*) 'err: we may need to increase max_num_cols.'
      stop 1
    endif
  end function num_of_integers

  function newunit()
    integer :: newunit
    integer, parameter :: lower_unitn=10, upper_unitn=1000
    logical :: is_opened
    integer :: i

    newunit = -1
    do i = lower_unitn, upper_unitn
      inquire(unit=i,opened=is_opened)
      if (.not. is_opened) then
        newunit = i
        exit
      endif
    enddo
    if(newunit == -1) then
      write(*,*) 'lower_unitn,upper_unitn=',lower_unitn,upper_unitn
      write(*,'(A)') 'err: cannot find a unused unit number.'
      stop 1
    endif
  end function newunit

  function rthetaphi(p)
    real(double) :: p(3),rthetaphi(3),testv(3),diffv(3)
    real(double) :: rmag,phi,theta,rho,q(3),pmax,diff

    real(double),parameter :: mag_tol=1.0d-12

    real(double),parameter :: eps=1.0d-11

    q(1:3) = p(1:3)
    rmag = vecmag3(p)
    if(rmag < mag_tol) then
      write(*,*) 'p = ',p
      write(*,*) 'rmag = ',rmag
      write(*,*) 'mag_tol = ',mag_tol
      write(*,'(A)') 'err: The magnitude of p is too small (less than mag_tol)'
      stop 1
    else

      pmax = abs(p(1))
      if(abs(p(2)) > pmax) then
        pmax = abs(p(2))
      endif
      if(abs(p(3)) > pmax) then
        pmax = abs(p(3))
      endif

      q(1:3) = q(1:3)/pmax

    endif

    phi = atan2(q(2),q(1))

    if(sqrt(q(1)*q(1) + q(2)*q(2)) < 1.0d-12) then

      phi = zero

    endif

    if(phi < 0) then
      phi = 2*pi + phi
    endif

    rho = sqrt(q(1)*q(1)+q(2)*q(2))
    theta = atan2(rho,q(3))

    rthetaphi(1:3) = (/rmag,theta,phi/)

    testv(1:3) = (/rmag*sin(theta)*cos(phi),rmag*sin(theta)*sin(phi),rmag*cos(theta)/)
    diffv(1:3) = testv(1:3) - p(1:3)
    diff = vecmag3(diffv)
    if(diff > eps) then
      write(*,*) 'eps,diff= ',eps,diff
    endif

  end function rthetaphi

  function neta_2_rotm(n,eta) result (m)
    real(double) :: m(3,3),n(3),eta,a,b,c,ct,st,nmag

    ct = cos(eta*deg2rad)
    st = sin(eta*deg2rad)

    nmag = sqrt(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))
    if(abs(nmag) < 1.0d-5) then
      write(*,'(A)') 'err: in rotmat: mag u is too small.'
      stop 1
    endif
    a = n(1)/nmag; b = n(2)/nmag; c = n(3)/nmag
    m(1,1) = a*a*(1-ct)+ct; m(1,2) = a*b*(1-ct)-c*st; m(1,3) = a*c*(1-ct)+b*st
    m(2,1) = a*b*(1-ct)+c*st; m(2,2) = b*b*(1-ct)+ct; m(2,3) = b*c*(1-ct)-a*st
    m(3,1) = a*c*(1-ct)-b*st; m(3,2) = b*c*(1-ct)+a*st; m(3,3) = c*c*(1-ct)+ct
  end function neta_2_rotm

  function euler_2_rotm(alpha,beta,gamm) result(m)
    real(double) :: alpha,beta,gamm
    real(double) :: a,b,g
    real(double) :: m(3,3)
    real(double) :: sa,ca,sb,cb,sg,cg

    a = alpha*deg2rad
    b = beta*deg2rad
    g = gamm*deg2rad
    sa = sin(a)
    ca = cos(a)
    sb = sin(b)
    cb = cos(b)
    sg = sin(g)
    cg = cos(g)
    m(1,1:3) = (/ca*cb*cg-sa*sg,-ca*cb*sg-sa*cg,ca*sb/)
    m(2,1:3) = (/sa*cb*cg+ca*sg,-sa*cb*sg+ca*cg,sa*sb/)
    m(3,1:3) = (/-sb*cg,sb*sg,cb/)
  end function euler_2_rotm

  function euler_2_neta(EulerAngle) result (p)
    real(double) :: EulerAngle(3),p(4)
    real(double) :: heta,alpha,beta,gamm,cosheta,sinheta,eta,a,b,c

    alpha = EulerAngle(1)*deg2rad
    beta = EulerAngle(2)*deg2rad
    gamm = EulerAngle(3)*deg2rad
    cosheta = cos(beta/two)*cos((alpha+gamm)/two)

    sinheta = sqrt(one - cosheta**two)
    heta = atan2(sinheta,cosheta)
    eta = two*heta

    if(eta < 1.0d-12) then
      a = one; b = one; c = one
    else

      a = - sin(beta/two)*sin((alpha-gamm)/two)/sin(eta/two)
      b = sin(beta/two)*cos((alpha-gamm)/two)/sin(eta/two)
      c = cos(beta/two)*sin((alpha+gamm)/two)/sin(eta/two)
    endif

    p(1) = a
    p(2) = b
    p(3) = c
    p(4) = eta*rad2deg
  end function euler_2_neta

  function rotm_2_euler(m) result(EulerAngle)
    real(double) :: sinbeta,m(3,3),EulerAngle(3)
    real(double) :: alpha,beta,gamm,d1pd4,d2md3,sum_alpha_gamm
    real(double),parameter :: eps=1.0d-12
    real(double) :: diff,diffmat(3,3),checkmat(3,3)
    integer :: toprint

    sinbeta = sqrt(m(1,3)**2 + m(2,3)**2)

    gamm = atan2(m(3,2),-m(3,1))*rad2deg
    alpha = atan2(m(2,3),m(1,3))*rad2deg
    beta = atan2(sinbeta,m(3,3))*rad2deg

    toprint = 0

    if(toprint == 1) then
      write(*,*) 'sinbeta,s(3,3)= ',sinbeta,m(3,3)
      write(*,*) 'alpha = ',alpha
      write(*,*) 'beta = ',beta
      write(*,*) 'gamm = ', gamm

      if(abs(sinbeta) < eps) then
        write(*,*) 'sinbeta= ',sinbeta
        write(*,*) 'WARNING: sinbeta is too small. But we are okay if we do not handle near 2-fold rotations'
      endif
    endif

    if(abs(beta) > 1.0d-10) then

    else
      if(toprint == 1) then
        write(*,*) 'WARNING: beta is too close to zero, we symmetrize alpha and gamma'
        write(*,*) 'before special case handling...'
        write(*,*) 'alpha = ',alpha
        write(*,*) 'beta = ',beta
        write(*,*) 'gamm = ', gamm
      endif

      d1pd4 = m(1,1)+m(2,2)
      d2md3 = m(2,1)-m(1,2)
      sum_alpha_gamm = atan2(d2md3,d1pd4)
      alpha = (sum_alpha_gamm/2d0)*rad2deg
      gamm = (sum_alpha_gamm/2d0)*rad2deg

      toprint = 0
      if(toprint == 1) then
        write(*,*) 'after special case handling...'
        write(*,*) 'rotm_2_euler: alpha = ',alpha
        write(*,*) 'rotm_2_euler: beta = ',beta
        write(*,*) 'rotm_2_euler: gamm = ', gamm
      endif
    endif

    EulerAngle(1:3) = (/alpha,beta,gamm/)

    checkmat(1:3,1:3) = euler_2_rotm(EulerAngle(1),EulerAngle(2),EulerAngle(3))

    diffmat(1:3,1:3) = checkmat(1:3,1:3) - m(1:3,1:3)
    diff = matmagn(3,diffmat(1,1))

    if(diff > eps) then
      write(*,*) 'rotm_2_euler: eps,diff= ',eps,diff
    endif
  end function rotm_2_euler

  function rotm_2_neta(m) result(p)
    real(double) :: m(3,3),p(4),euler(3), checkm(3,3),diffm(3,3)
    real(double) :: eps
    integer :: info
    integer,parameter :: ldvl=3
    integer,parameter :: ldvr=3
    complex(double) :: VL(ldvl,3),VR(ldvr,3)
    integer,parameter :: LCWORK=6
    integer,parameter :: LRWORK=6
    complex(double) :: CWORK(LCWORK),cmat(3,3)
    real(double) :: RWORK(LRWORK)
    complex(double) :: eigval(3),sumeigval
    integer :: i
    real(double) :: rsum,costheta,sintheta,theta,xsq,tracem
    logical :: found,near2fold
    real(double) :: diff
    integer :: toprint

    eps = 1.0d-7

    call check_orthogonal(m(1,1),eps)

    tracem = m(1,1) + m(2,2) + m(3,3)

    cmat(1:3,1:3) = m(1:3,1:3)

    call zgeev('V','V',3,cmat(1,1),3,eigval(1),VL,LDVL,VR,LDVR,CWORK(1),LCWORK,RWORK(1),info)
    if(info /= 0) then
      write(*,'(A)') 'err: error in zgeev.'
      stop 1
    endif

    sumeigval = eigval(1)+eigval(2)+eigval(3)
    rsum = real(sumeigval)

    if(abs(rsum - tracem) > eps) then
      write(*,*) 'rsum,tracem= ',rsum,tracem
      write(*,'(A)') 'err: rsum and tracem must be the same.'
      stop 1
    endif

    if(abs ( abs(tracem) - 3.0d0)  < 1.0d-12) then
      write(*,*) 'Hit near identity: return a reasonabe n and eta'
      p(1:3) = one
      p(4) = zero
      return
    endif

    near2fold = .false.

    if(abs(rsum - (-one)) < eps) then

      found = .false.
      do i = 1, 3

        if(found) cycle

        if( abs (real(eigval(i) - one))  < eps) then
          found = .true.
          near2fold = .true.

          costheta = (rsum- one)/two
          xsq = one- costheta*costheta
          if(xsq < zero) then

            xsq = -xsq
          endif
          sintheta = sqrt(xsq)
          theta = atan2(sintheta,costheta)
          p(4) = theta*rad2deg
          p(1) = real(vr(1,i))
          p(2) = real(vr(2,i))
          p(3) = real(vr(3,i))
        endif
      enddo
    endif

    if(.not. near2fold) then

      euler = rotm_2_euler(m)

      p = euler_2_neta(euler)
    endif

    checkm = neta_2_rotm( (/p(1),p(2),p(3)/), p(4) )
    diffm = checkm-m
    diff = matmagn(3,diffm(1,1))

    toprint = 0
    if(toprint == 1) then
      write(*,*) 'm,checkm,diffm='
      write(*,'(3F8.4,A,3F8.4,A,3F8.4)') m(1,1:3), '  |  ', checkm(1,1:3),'  |  ',diffm(1,1:3)
      write(*,'(3F8.4,A,3F8.4,A,3F8.4)') m(2,1:3), '  |  ', checkm(2,1:3),'  |  ',diffm(2,1:3)
      write(*,'(3F8.4,A,3F8.4,A,3F8.4)') m(3,1:3), '  |  ', checkm(3,1:3),'  |  ',diffm(3,1:3)
      write(*,*) 'rotm_2_neta: eps,diff= ',eps,diff
    endif
    if(diff > eps) then
      write(*,'(A)') 'err: in rotm_2_neta: m and p not consistent.'
      stop 1
    endif
  end function rotm_2_neta

  subroutine getlenang(a,L)
    real(double) :: a(3,3),L(6)
    real(double) :: sintheta,costheta

    L(1) = vecmag3(a(1:3,1))
    L(2) = vecmag3(a(1:3,2))
    L(3) = vecmag3(a(1:3,3))

    costheta = dotprod3(a(1:3,2),a(1:3,3))/(L(2)*L(3))
    sintheta = sqrt((1d0-costheta)*(1+costheta))
    L(4) = atan2(sintheta,costheta)*rad2deg

    costheta = dotprod3(a(1:3,3),a(1:3,1))/(L(3)*L(1))
    sintheta = sqrt((1d0-costheta)*(1+costheta))
    L(5) = atan2(sintheta,costheta)*rad2deg

    costheta = dotprod3(a(1:3,1),a(1:3,2))/(L(1)*L(2))
    sintheta = sqrt((1d0-costheta)*(1+costheta))
    L(6) = atan2(sintheta,costheta)*rad2deg
  end subroutine getlenang

  subroutine real2recip(reallatt,reciplatt)
    real(double) :: reallatt(3,3),reciplatt(3,3)
    real(double) :: vol,tmp2(3,3)

    vol = det3(reallatt(1,1))
    if(vol < 0) then
      write(*,*) 'vol = ',vol
      write(*,'(A)') 'err: vol is negative. Check reallatt.'
      stop 1
    endif

    tmp2 = inv3x3(reallatt(1,1))

    reciplatt(1:3,1) = tmp2(1,1:3)
    reciplatt(1:3,2) = tmp2(2,1:3)
    reciplatt(1:3,3) = tmp2(3,1:3)
  end subroutine real2recip

  function angle_between_vecs3(A,B) result (C)
    real(double) :: A(3),B(3),C,maga,magb,AdotB,xsq,costheta,sintheta
    real(double),parameter :: eps=1.0d-8

    maga = vecmag3(A(1))
    magb = vecmag3(B(1))

    if(maga < eps .or. magb < eps) then
      write(*,'(A)') 'err: maga or magb too small.'
      stop 1
    endif
    adotb = dotprod3(A(1),B(1))
    costheta = adotb/(maga*magb)
    xsq = one-costheta*costheta
    if(xsq < zero) then
      xsq = -xsq
    endif
    sintheta = sqrt(xsq)
    c = atan2(sintheta,costheta)
    c = c*rad2deg
  end function angle_between_vecs3

  function crossprod(A,B)  result(C)
    real(double) :: A(1:3),B(1:3),C(1:3)
    C(1) = A(2)*B(3) - A(3)*B(2)
    C(2) = A(3)*B(1) - A(1)*B(3)
    C(3) = A(1)*B(2) - A(2)*B(1)
  end function crossprod

  function dotprod3(A,B)
    real(double) :: A(1:3),B(1:3),dotprod3
    dotprod3 = dotprodn(3,a(1),b(1))
  end function dotprod3

  function dotprodn(n,A,B)
    integer :: n,i
    real(double) :: A(n),B(n),dotprodn
    dotprodn = zero
    do i = 1, n
      dotprodn = dotprodn + A(i)*B(i)
    enddo
  end function dotprodn

  function tripprod(A)
    real(double) :: A(3,3),tripprod
    tripprod = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) - &
               A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) + &
               A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
  end function tripprod

  function det3(a)
    real(double) :: det3,a(3,3)
    det3 = tripprod(a)
  end function det3

  function tracen(n,a)
    integer :: n,i
    real(double) :: tracen,a(n,n)
    tracen = zero
    do i = 1, n
      tracen = tracen + a(i,i)
    enddo
  end function tracen

  function trace3(a)
    real(double) :: trace3,a(3,3)
    trace3 = tracen(3,a(1,1))
  end function trace3

  subroutine frac2abs(reallatt,FracC,AbsC)
    real(double) :: FracC(1:3),AbsC(1:3),reallatt(3,3)
    integer :: i,j
    do i = 1, 3
      AbsC(i) = 0
      do j = 1, 3
        AbsC(i) = AbsC(i) + reallatt(i,j)*FracC(j)
      enddo
    enddo
  end subroutine frac2abs

  subroutine abs2frac(reciplatt,AbsC,FracC)
    real(double) :: AbsC(3),FracC(3),reciplatt(3,3)
    integer :: i
    do i = 1, 3
      FracC(i) = dotprod3(AbsC(1:3),reciplatt(1:3,i))
    enddo
  end subroutine abs2frac

  function matmagn(n,m)
    integer :: n
    real(double) :: matmagn,m(n,n)
    integer :: i,j
    matmagn = zero
    do j = 1, n
      do i = 1, n
        matmagn = matmagn + m(i,j)**2
      enddo
    enddo
    matmagn = sqrt(matmagn)
  end function matmagn

  function vecmag3(x)
    real(double) :: vecmag3,x(3)
    vecmag3 = sqrt(dotprod3(x,x))
  end function vecmag3

  function vecmagn(n,x)
    integer :: n,i
    real(double) :: mag,vecmagn,x(n)
    mag = 0d0
    do i = 1, n
      mag = mag + x(i)**2
    enddo
    vecmagn = sqrt(mag)
  end function vecmagn

  function inv3x3(a) result(inv)
    real(double) :: a(3,3),inv(3,3)
    real(double) :: det,absdet

    inv(1,1) = a(2,2)*a(3,3)-a(2,3)*a(3,2)
    inv(2,1) = -(a(2,1)*a(3,3) - a(2,3)*a(3,1))
    inv(3,1) = a(2,1)*a(3,2) - a(2,2)*a(3,1)
    inv(1,2) = -(a(1,2)*a(3,3) - a(1,3)*a(3,2))
    inv(2,2) = a(1,1)*a(3,3) - a(1,3)*a(3,1)
    inv(3,2) = -(a(1,1)*a(3,2) - a(1,2)*a(3,1))
    inv(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
    inv(2,3) = -(a(1,1)*a(2,3) - a(1,3)*a(2,1))
    inv(3,3) = a(1,1)*a(2,2) - a(1,2)*a(2,1)

    det = a(1,1)*inv(1,1)+a(1,2)*inv(2,1)+a(1,3)*inv(3,1)
    absdet = abs(det)

    if(absdet < 1.0d-13) then
      write(*,*) 'absdet = ',absdet
      write(*,'(A)') 'err: det is too small in inv3x3.'
      stop 1
    endif
    inv = inv/det
  end function inv3x3

  subroutine matmat3(a,b,c)
    real(double):: a(3,3),b(3,3),c(3,3)
    call matmatn(3,a(1,1),b(1,1),c(1,1))
  end subroutine matmat3

  subroutine matmatn(n,a,b,c)
    integer :: i,j,k,n
    real(double) :: a(n,n),b(n,n),c(n,n),sumr
    do i = 1, n
      do j = 1, n
        sumr = zero
        do k = 1,n
          sumr = sumr + a(i,k)*b(k,j)
        enddo
        c(i,j) = sumr
      enddo
    enddo
  end subroutine matmatn

  subroutine matmatmnp(m,n,p,a,b,c)
    integer :: i,j,k,n,m,p
    real(double) :: a(m,n),b(n,p),c(m,p)

    do i = 1, m
      do j = 1, p
        c(i,j) = 0D0
        do k = 1,n
          c(i,j) = c(i,j) + a(i,k)*b(k,j)
        enddo
      enddo
    enddo
  end subroutine matmatmnp

  subroutine matvecn(n,a,b,c)
    integer :: j,k,n
    real(double) :: a(n,n),b(n),c(n)
    do j = 1, n
      c(j) = 0D0
      do k = 1,n
        c(j) = c(j) + a(j,k)*b(k)
      enddo
    enddo
  end subroutine matvecn

  subroutine matvec3(a,b,c)
    real(double) :: a(3,3),b(3),c(3)
    call matvecn(3,a(1,1),b(1),c(1))
  end subroutine matvec3

  function rotmat_x(t) result(m)
    real(double) :: m(3,3),t,u(3)
    u(1:3) = (/1,0,0/)
    m = neta_2_rotm(u,t)
  end function rotmat_x

  function rotmat_y(t) result(m)
    real(double) :: m(3,3),t,u(3)
    u(1:3) = (/0,1,0/)
    m = neta_2_rotm(u,t)
  end function rotmat_y

  function rotmat_z(t) result(m)
    real(double) :: m(3,3),t,u(3)
    u(1:3) = (/0,0,1/)
    m = neta_2_rotm(u,t)
  end function rotmat_z

  subroutine check_orthogonal(R,eps)
    real(double) :: diffabs,R(3,3),eps,tmp1(3,3),tmp2(3,3),tmp3(3,3),diffm(3,3)
    integer :: i,j

    do i = 1, 3
      do j = 1, 3
        tmp1(i,j) = R(j,i)
      enddo
    enddo
    call matmat3(R(1,1),tmp1(1,1),tmp2(1,1))
    tmp3 = zero
    tmp3(1,1) = one
    tmp3(2,2) = one
    tmp3(3,3) = one
    diffm = tmp2 - tmp3
    diffabs = matmagn(3,diffm(1,1))
    if(diffabs > eps) then
      write(*,*) ' R= '
      write(*,'(3F12.5)') R(1,1:3)
      write(*,'(3F12.5)') R(2,1:3)
      write(*,'(3F12.5)') R(3,1:3)
      write(*,*) 'diffabs = ',diffabs
      write(*,*) 'Are you sure that the matrix R is orthogonal ?'
      write(*,'(A)') 'err: R is not orthogonal within eps.'
      stop 1
    endif
  end subroutine check_orthogonal

  function invBAB(B,A) result (m)
    real(double) :: B(3,3),A(3,3),m(3,3),tmp(3,3)
    call matmatn(3,inv3x3(B),A,tmp)
    call matmatn(3,tmp,B,m)
  end function invBAB

  function BAinvB(B,A) result (m)
    real(double) :: B(3,3),A(3,3),m(3,3),tmp(3,3)
    call matmatn(3,B,A,tmp)
    call matmatn(3,tmp,inv3x3(B),m)
  end function BAinvB

  subroutine modulo1(x1)
    real(double) :: x1
    integer :: xint
    real(double) :: pluseps=1.0d-9

    if(x1 < 0.0d0) then

      xint = int(x1)
      x1 = x1 - xint + 1d0
      if(x1 == 1.0d0) x1 = 0d0
    else if(x1 >= 1.0d0) then

      xint = int(x1)
      x1 = x1 - xint
    endif

    if(x1 >= one .or. x1 < zero) then
      write(*,*) 'x1 is ',x1
      write(*,'(A)') 'err: must stop, x1 must be in [0,1) '
      stop 1
    endif

    if(abs(x1 - one) < pluseps) then

      x1 = zero
    endif
  end subroutine modulo1

  subroutine foldtofirstzone3(x)
    real(double) :: x(3)
    integer :: i
    do i = 1, 3
      call modulo1(x(i))
    enddo
  end subroutine foldtofirstzone3

  subroutine supercell_2_vasp(unitn,filename,atmtype,zord,s,writeforce,ATMSswitch)
    integer,optional :: writeforce
    integer,optional :: ATMSswitch
    integer :: unitn
    character(len=*) :: filename
    type(supercell) :: s
    integer :: atmtype,zord(atmtype)
    integer :: totat,i,ind,sizei,j,n
    integer :: groupsize(200)
    type(oneatom) :: outcell(maxnum)
    integer,allocatable :: fmap(:)
    character(len=DSL) :: commentline
    real(double),parameter :: maxlen=1000d0
    real(double) :: len6(6),frac(3),cartesian(3),num3(3)
    integer,allocatable :: localbmap(:)
    integer :: k

    n = s%n
    allocate(fmap(n))
    allocate(localbmap(n))

    ind = 0
    totat = 0
    do i = 1, atmtype
      sizei = 0
      do j = 1, n
        if(s%at(j)%z == zord(i)) then
          sizei = sizei + 1
          ind = ind + 1
          fmap(ind) = j
          outcell(ind)%z = zord(i)
          outcell(ind)%p(1:3) = s%at(j)%f(1:3)
          outcell(ind)%force(1:3) = s%at(j)%force(1:3)
        endif
      enddo
      groupsize(i) = sizei
      totat = totat + groupsize(i)
    enddo
    if(totat /= n) then
      write(*,*) 'totat,n=',totat,n
      write(*,'(A)') 'err: totat and n problem in supercell_2_vasp'
      stop 1
    endif
    if(ind /= n) then
      write(*,'(A,2I5)') 'ind,n = ',ind,n
      write(*,'(A)') 'err: new ind and n are not the same.'
      stop 1
    endif

    do i = 1, atmtype
      if(groupsize(i) < 1) then
        write(*,'(A,I4,A,I4)') 'i, groupsize(i)=',i,',',groupsize(i)
        write(*,'(A)') 'err: supercell_2_vasp: groupsize(i) is less than 1.'
        stop 1
      endif
    enddo
    do i = 1, n
      j = fmap(i)
      localbmap(j) = i
    enddo
    open(unit=unitn,file=trim(filename),status='replace')

    if(present(ATMSswitch)) then
      if(ATMSswitch == 1) then
        commentline = 'ATMS: '//trim(N2str(n))//' '//trim(N2str(atmtype))
        do i = 1, atmtype

          commentline = trim(commentline)//' '//trim(N2str(groupsize(i)))//' '//trim(Cats(zord(i)))
        enddo
        commentline = trim(commentline)//' SG '//trim(s%sg_label)
      endif
    else
      commentline = ''

      if(trim(s%sg_label) /= '' .and. trim(s%sg_label) /= 'undefined') then
        commentline = 'SG '//trim(s%sg_label)
      endif
    endif

    write(unitn,'(A)') trim(commentline)

    write(unitn,'(F6.3)') 1.0d0
    write(unitn,'(3F18.12)') s%a(1:3,1)
    write(unitn,'(3F18.12)') s%a(1:3,2)
    write(unitn,'(3F18.12)') s%a(1:3,3)

    commentline=''
    do i = 1, atmtype
      commentline = trim(commentline)//' '//trim(CATS(zord(i)))
    enddo
    write(unitn,'(A)') trim(commentline)
    write(unitn,'(10I4)') (groupsize(i),i=1,atmtype)
    if(trim(s%fracabs) == 'frac') then
      write(unitn,'(A)') 'Direct'
    else if(trim(s%fracabs) == 'abs') then
      write(unitn,'(A)') 'Cartesian'
    else
      write(*,'(A)') 'err: Must be either frac or abs.'
      stop 1
    endif

    ind = 0
    do i = 1, atmtype
      do j = 1, groupsize(i)
        ind = ind + 1
        frac(1:3) = outcell(ind)%p(1:3)

        call frac2abs(s%a,frac(1),cartesian(1))

        if(trim(s%fracabs) == 'frac') then
          num3(1:3) = frac(1:3)

          do k = 1, 3
            call modulo1(num3(k))
          enddo
        else if(trim(s%fracabs) == 'abs') then
          num3(1:3) = cartesian(1:3)
        endif
        if(.not. present(writeforce)) then
          write(unitn,'(3F18.12)') num3(1:3)
        else
          write(unitn,'(3F18.12,A,3F18.12)') num3(1:3),' | ',outcell(ind)%force(1:3)
        endif
      enddo
    enddo
    close(unitn)

    call getlenang(s%a,len6(1))
    deallocate(fmap)
    deallocate(localbmap)
  end subroutine supercell_2_vasp

  subroutine onespace(ns,comment1,comment)
    integer :: ns
    character(len=ns) :: comment1,comment
    integer :: i,j

    comment = 's'

    i = 1
    j = 1

    loop2: do
      loop1: do
        if(comment1(i:i) /= ' ') exit loop1
        i = i + 1
        if(i == ns) exit loop2
      enddo loop1

      loop3: do
        if(comment1(i:i) == ' ') exit loop3
        comment(j:j) = comment1(i:i)
        j = j + 1
        i = i + 1
      enddo loop3
      comment(j:j) = ' '
      j = j + 1
    enddo loop2
  end subroutine onespace

  subroutine assign_str_fr_struct_out(inu,s)
    integer :: inu
    type(supercell) :: s
    integer :: i

    do i = 1, 3
      read(inu,*) s%a(1:3,i)
    enddo
    read(inu,*) s%n
    do i = 1, s%n
      read(inu,*) s%at(i)%gr, s%at(i)%z,s%at(i)%f(1:3)
    enddo
    call real2recip(s%a,s%b)
    close(inu)
  end subroutine assign_str_fr_struct_out

  subroutine assign_str_fr_xsf(controlu,sc)
    integer :: i,controlu
    type(supercell) :: sc
    character(LEN=DSL) :: tmp

    read(controlu,*) tmp
    if(trim(tmp) /= 'CRYSTAL') then
      write(*,'(A)') 'err: Keyword CRYSTAL not found.'
      stop 1
    endif
    read(controlu,*) tmp
    if(trim(tmp) /= 'PRIMVEC') then
      write(*,'(A)') 'err: Keyword PRIMVEC not found.'
      stop 1
    endif

    read(controlu,*) sc%a(1:3,1)
    read(controlu,*) sc%a(1:3,2)
    read(controlu,*) sc%a(1:3,3)
    read(controlu,*) tmp
    if(trim(tmp) /= 'CONVVEC') then
      write(*,'(A)') 'err: Keyword CONVVEC not found.'
      stop 1
    endif

    read(controlu,*) sc%a(1:3,1)
    read(controlu,*) sc%a(1:3,2)
    read(controlu,*) sc%a(1:3,3)
    call real2recip(sc%a,sc%b)
    read(controlu,*) tmp
    if(trim(tmp) /= 'PRIMCOORD') then
      write(*,'(A)') 'err: Keyword PRIMCOORD not found.'
      stop 1
    endif
    read(controlu,*) sc%n
    do i = 1, sc%n
      read(controlu,*) sc%at(i)%z, sc%at(i)%ac
      call abs2frac(sc%b,sc%at(i)%ac,sc%at(i)%f)
    enddo
    sc%sg_label = '1-P1'
  end subroutine assign_str_fr_xsf

  subroutine assign_str_fr_poscar(controlu,sc)
    integer :: controlu
    type(supercell) :: sc
    real(double) :: totalscaling,vol
    integer :: i,j,s,nsp,v,z,ind
    character(LEN=DSL) :: comment1,comment2,tmp1,tmp2,valuestr,targetlabel
    real(double) :: pos(3)
    integer,allocatable :: numarr(:)
    integer :: totn,k,labellen,SGnumber,loc
    logical :: done,found,puredigits
    integer :: nvalid
    real(double) :: density

    read(controlu,DSF) comment1

    call onespace(DSL,comment1,comment2)

    sc%commentline = trim(comment2)
    write(*,'(A,A,A)') 'POSCAR comment is [',trim(comment2),'].'
    s = scan(comment2,' ')

    if(comment2(1:s-1) == 'ATMS:') then
      tmp1= comment2(s+1:DSL)

      read(tmp1,*) v
      sc%n = v
      if(sc%n < 1) then
        write(*,*) 'sc%n = ',sc%n
        write(*,'(A)') 'err: check sc%n in assign_str_fr_poscar'
        stop 1
      endif
      s = scan(tmp1,' ')
      tmp1 = tmp1(s+1:DSL)

      read(tmp1,*) v
      sc%nsp = v
      if(sc%nsp < 1) then
        write(*,*) 'sc%nsp = ',sc%nsp
        write(*,'(A)') 'err: check sc%nsp in assign_str_fr_poscar.'
        stop 1
      elseif(sc%nsp > zordmax) then
        write(*,*) 'sc%nsp,zordmax=',sc%nsp,zordmax
        write(*,'(A)') 'err: check sc%nsp, i.e., the number of species/elements'
        stop 1
      endif
      s = scan(tmp1,' ')
      tmp1 = tmp1(s+1:DSL)
      nsp = v
      if(sc%nsp > zordmax) then
        write(*,*) 'sc%nsp = ',sc%nsp
        write(*,*) 'zordmax = ',zordmax
        write(*,'(A)') 'err: out of bound.'
        stop 1
      endif

      do i = 1, nsp

        read(tmp1,*) v

        sc%n_per_species(i) = v

        s = scan(tmp1,' ')
        tmp1 = tmp1(s+1:DSL)

        s = scan(tmp1,' ')
        valuestr = tmp1(1:s-1)
        found = .false.
        do j = 1, zordmax
          if(found) exit
          if(trim(valuestr) == trim(Ats(j)) .or. trim(valuestr) == trim(CAts(j)) ) then
            z = j
            found = .true.
          endif
        enddo
        if(.not. found) then
          write(*,*) 'valuestr = ',trim(valuestr)
          write(*,'(A)') 'err: this valuestr is not found in the periodic table.'
          stop 1
        endif
        sc%zarr(i) = z
        tmp1 = tmp1(s+1:DSL)
      enddo

      ind = 0
      do j = 1, sc%nsp
        ind = ind + sc%n_per_species(j)
      enddo

      if(ind /= sc%n) then
        write(*,*) 'ind,sc%n=',ind,sc%n
        write(*,'(A)') 'err: missing atoms in the assign_str_fr_poscar.'
        stop 1
      endif

      sc%sg_label = ''
      if(trim(tmp1) /= '') then

        s = scan(tmp1,' ')
        valuestr = tmp1(1:s-1)
        if(trim(valuestr) /= 'SG') then
          write(*,*) 'valuestr is '//trim(valuestr)
          write(*,*) 'The first optional parameter must be SG.'
          write(*,'(A)') 'err: Check poscar.'
          stop 1
        endif
        tmp1 = tmp1(s+1:DSL)
        s = scan(tmp1,' ')
        valuestr = tmp1(1:s-1)
        sc%sg_label = trim(valuestr)
      else
        write(*,*) 'No SG label'
      endif

      read(controlu,*) totalscaling
      read(controlu,*) sc%a(1:3,1)
      read(controlu,*) sc%a(1:3,2)
      read(controlu,*) sc%a(1:3,3)

      allocate(numarr(sc%nsp))

      read(controlu,DSF) tmp1
      write(*,*) 'tmp1 is .'//trim(tmp1)//'.'
      call onespace(DSL,tmp1,tmp2)
      s = scan(tmp2,'0123456789')

      if(s == 0) then

        write(*,'(A)') 'We are going to process the species line |'//trim(tmp2)//'|'

        nsp = num_of_strings(tmp2)
        if(nsp /= sc%nsp) then
          write(*,*) 'err: nsp,sc%nsp = ',nsp,sc%nsp
          stop 1
        endif
        do i = 1, sc%nsp
          s = scan(tmp2,' ')
          valuestr = tmp2(1:s-1)
          write(*,*) 'valuestr is '//trim(valuestr)
          found = .false.
          do j = 1, zordmax
            if(found) exit
            if(trim(valuestr) == trim(CATS(j))) then
              z = j
              if(z /= sc%zarr(i)) then
                write(*,*) 'z, sc%zarr(i) = ',z,sc%zarr(i)
                write(*,'(A)') 'err: Inconsistency.'
                stop 1
              endif
              found = .true.
            endif
          enddo
          if( .not. found) then
            write(*,'(A)') 'err: not found.'
            stop 1
          endif
          tmp2 = tmp2(s+1:DSL)
        enddo

        read(controlu,DSF) tmp2
        nsp = num_of_integers(tmp2)
        if(nsp /= sc%nsp) then
          write(*,*) 'err: nsp,sc%nsp = ',nsp,sc%nsp
          stop 1
        endif

        read(tmp2,*) numarr(1:sc%nsp)
      else

        nsp = num_of_integers(tmp1)
        if(nsp /= sc%nsp) then
          write(*,*) 'err: nsp,sc%nsp = ',nsp,sc%nsp
          stop 1
        endif
        read(tmp1,*) numarr(1:sc%nsp)
      endif

      do i = 1, sc%nsp
        if(numarr(i) /= sc%n_per_species(i)) then
          write(*,*) 'species: i=',i
          write(*,*) 'numarr(i) = '//trim(N2str(numarr(i)))//', sc%n_per_species(i) = '//trim(N2str(sc%n_per_species(i)))
          write(*,'(A)') 'err: numarr(i) and sc%n_per_species(i) are not consistent'
          stop 1
        endif
      enddo
      deallocate(numarr)

    else
      if(comment2(1:s-1) == 'SG') then

        tmp2 = comment2(s+1:DSL)
        targetlabel = trim(tmp2)

        write(*,*) 'targetlabel = '//trim(targetlabel)

        puredigits = .true.
        labellen = len(trim(targetlabel))
        do i = 1, labellen
          if(.not. puredigits) then
            exit
          endif
         loc = scan(targetlabel(i:i),'0123456789')
         if(loc /= 1) then
           puredigits = .false.
         endif
        enddo

        if(puredigits)  then
          SGnumber = str2N(trim(targetlabel))
          write(*,*) 'SGnumber is ',SGnumber
          if(SGnumber >= 1 .and. SGnumber <= 230) then
            targetlabel = trim(SGbase(1,SGnumber)%sglabel)
            write(*,*) 'Fully expand the default SG label to the full SG label:'
            write(*,*) 'SGnumber= ',SGnumber
            write(*,*) 'Full targetlabel= ',trim(targetlabel)
          else
            write(*,*) 'err: must be between 1 and 230 inclusive. impossible. '
            stop 1
          endif
        endif
        sc%sg_label = trim(targetlabel)
        write(*,'(A)') 'Possibly adjusted sg_label is '//trim(sc%sg_label)
      else
        write(*,'(A)') 'WARNING WARNING: No ATMS:, no SG, hence we assume SG 1'
        sc%sg_label = '1-P1'
      endif

      read(controlu,*) totalscaling
      read(controlu,*) sc%a(1:3,1)
      read(controlu,*) sc%a(1:3,2)
      read(controlu,*) sc%a(1:3,3)

      read(controlu,DSF) tmp1
      call onespace(DSL,tmp1,tmp2)

      s = scan(tmp2,'0123456789')

      if(s == 0) then

        nsp = 0
        done = .false.
        do
          if(done) exit
          s = scan(tmp2,' ')
          valuestr = tmp2(1:s-1)

          found = .false.
          do j = 1, zordmax
            if(found) exit
            if(trim(valuestr) == trim(CATS(j))) then
              nsp = nsp + 1
              sc%zarr(nsp) = j
              found = .true.
            endif
          enddo
          if( .not. found) then
            write(*,'(A)') 'err: not found.'
            stop 1
          endif
          tmp2 = tmp2(s+1:DSL)
          if(len(trim(tmp2)) == 0) then
            write(*,'(A)') 'nsp = '//trim(N2str(nsp))
            sc%nsp = nsp
            done = .true.
          endif
        enddo

        read(controlu,DSF) tmp1

        nvalid = num_of_integers(tmp1)
        if(nsp /= nvalid) then
          write(*,*) 'nsp,nvalid=',nsp,nvalid
          write(*,*) 'err: inconsistency in species line and species line.'
          stop 1
        endif
        read(tmp1,*) sc%n_per_species(1:nsp)

      else
        write(*,'(A)') 'err: Since this is without ATMs, we must have species line.'
        stop 1
      endif
    endif

    totn = 0
    do j = 1, sc%nsp
      v = sc%n_per_species(j)
      do k = 1, v
        totn = totn + 1
        sc%at(totn)%z = sc%zarr(j)
        sc%at(totn)%mass = massofa(sc%zarr(j))
      enddo
    enddo
    write(*,'(A)') 'totn = '//trim(N2str(totn))
    sc%n = totn

    sc%a = sc%a*totalscaling

    call real2recip(sc%a,sc%b)

    write(*,'(A,/,3F15.10,A,F15.10,/,3F15.10,A,F15.10,/,3F15.10,A,F15.10)') &
      'Real latt:', sc%a(1:3,1),'    | ',vecmag3(sc%a(1:3,1)),sc%a(1:3,2),'    | ',vecmag3(sc%a(1:3,2)),sc%a(1:3,3),'    | ',vecmag3(sc%a(1:3,3))
    vol = det3(sc%a(1,1))
    if(vol .le. 0) then
      write(*,'(A)') 'err: vol is not positive.'
      stop 1
    endif

    sc%vol = vol

    write(*,'(A,/,3F15.10,A,F15.10,/,3F15.10,A,F15.10,/,3F15.10,A,F15.10)') &
      'Reciprocal latt:', sc%b(1:3,1),'    | ',vecmag3(sc%b(1:3,1)),sc%b(1:3,2),'    | ',vecmag3(sc%b(1:3,2)),sc%b(1:3,3),'    | ',vecmag3(sc%b(1:3,3))

    read(controlu,*) tmp1
    if(trim(tmp1) == 'Selective') then
      read(controlu,*) tmp1
    endif
    if(trim(tmp1) /= 'Direct' .and. trim(tmp1) /= 'direct' .and. trim(tmp1) /= 'Cartesian' ) then
      write(*,'(A)') 'err: dummy should be Direct, direct, or Cartesian'
      stop 1
    endif

    do i = 1, sc%n
      sc%at(i)%force(1:3) = (/zero,zero,zero/)
    enddo
    do i = 1, sc%n
      if(trim(tmp1) == 'Direct' .or. trim(tmp1) == 'direct') then
        read(controlu,DSF) comment1
        call onespace(DSL,comment1,comment2)
        read(comment2,*) sc%at(i)%f(1:3)

        call frac2abs(sc%a(1,1),sc%at(i)%f(1),sc%at(i)%ac(1))
        s = scan(comment2,'|')
        if(s > 0) then
          tmp2= comment2(s+1:DSL)

          read(tmp2,*) sc%at(i)%force(1:3)
          write(*,'(A,3F20.10)') 'sc%at(i)%force(1:3) = ',sc%at(i)%force(1:3)
        endif
      else if(trim(tmp1) == 'Cartesian') then

        read(controlu,DSF) comment1
        call onespace(DSL,comment1,comment2)
        read(comment2,*) pos(1:3)

        call abs2frac(sc%b(1,1),pos(1),sc%at(i)%f(1))

        call frac2abs(sc%a(1,1),sc%at(i)%f(1),sc%at(i)%ac(1))
        s = scan(comment2,'|')
        if(s > 0) then
          tmp2= comment2(s+1:DSL)

          read(tmp2,*) sc%at(i)%force(1:3)
          write(*,'(A,3F20.10)') 'sc%at(i)%force(1:3) = ',sc%at(i)%force(1:3)
        endif
      endif
    enddo
    call getlenang(sc%a,sc%la(1))
    write(*,'(A)') 'la(1:6)='
    do i = 1, 6
      write(*,'(3F15.10)') sc%la(i)
    enddo
    write(*,*) 'Volume of the supercell= ',vol
    density = crys_density(sc)
    sc%density = density
    write(*,*) 'Density is ',density,' g/cm^3'

  end subroutine assign_str_fr_poscar

  subroutine assign_str_fr_fdf(inu,s)
    integer :: inu,natom,nspecies,sp,i,j
    type(supercell) :: s
    character(len=DSL) :: str1,str2,str3
    real(double) :: x(3)
    character(len=DSL) :: icf

    read(inu,*) str1,str2
    if(trim(str1) /= 'NumberOfAtoms') then
      write(*,*) 'str1 = ',trim(str1)
      write(*,*) 'but str1 must be NumberofAtoms'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    else
      read(str2,'(I8)') natom
      write(*,*) 'natom = ',natom
      s%n = natom
    endif
    read(inu,*) str1,str2
    if(trim(str1) /= 'NumberOfSpecies') then
      write(*,*) 'str1 = ',trim(str1)
      write(*,*) 'but str1 must be NumberofSpecies'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    else
      read(str2,'(I8)') nspecies
      write(*,*) 'nspecies = ',nspecies
      s%nsp = nspecies
    endif
    read(inu,*) str1,str2
    if(trim(str2) /= 'ChemicalSpeciesLabel') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be ChemicalSpeciesLabel'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    do i = 1, nspecies
      read(inu,*) j,sp
      if(j/=i) then
        write(*,*) 'j,i=',j,i
        write(*,*) 'sp = ',sp
        write(*,'(A)') 'err: something wrong with the species numbering ?'
        stop 1
      else
        s%zarr(i) = sp
        write(*,*) 's%zarr(',i,')=',sp
      endif
    enddo
    read(inu,*) str1,str2
    if(trim(str2) /= 'ChemicalSpeciesLabel') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be ChemicalSpeciesLabel'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) str1,str2,str3
    if(trim(str1) /= 'LatticeConstant') then
      write(*,*) 'str1 = ',trim(str1)
      write(*,*) 'but str1 must be LatticeConstant'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    if(trim(str2) /= '1.0') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be 1.0'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    if(trim(str3) /= 'Ang') then
      write(*,*) 'str3 = ',trim(str3)
      write(*,*) 'but str3 must be Ang'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) str1,str2
    if(trim(str2) /= 'LatticeVectors') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be  LatticeVectors'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) s%a(1:3,1)
    read(inu,*) s%a(1:3,2)
    read(inu,*) s%a(1:3,3)
    write(*,*) 's%a is '
    write(*,*) s%a(1:3,1)
    write(*,*) s%a(1:3,2)
    write(*,*) s%a(1:3,3)
    call real2recip(s%a,s%b)
    read(inu,*) str1,str2
    if(trim(str2) /= 'LatticeVectors') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be LatticeVectors'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) str1,str2
    if(trim(str2) == 'Fractional') then
      icf = 'frac'
    else if(trim(str2) == 'NotScaledCartesianAng') then
      icf = 'abs'
    else
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be either Fractional or NotScaledCartesianAng'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) str1,str2
    if(trim(str2) /= 'AtomicCoordinatesAndAtomicSpecies') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be AtomicCoordinatesAndAtomicSpecies'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    do i = 1, natom
      read(inu,*) x(1:3),sp
      if(trim(icf) == 'frac') then
        s%at(i)%f(1:3) = x(1:3)
      else if(trim(icf) == 'abs') then
        call abs2frac(s%b,x,s%at(i)%f)
      else
        write(*,'(A)') 'err: frac and abs problem.'
        stop 1
      endif
      s%at(i)%z = s%zarr(sp)

    enddo
    read(inu,*) str1,str2
    if(trim(str2) /= 'AtomicCoordinatesAndAtomicSpecies') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be AtomicCoordinatesAndAtomicSpecies'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
  end subroutine assign_str_fr_fdf

  function massofa(z)
    integer :: z
    real(double) :: massofa
    if(z < 0 .or. z > zordmax) then
      write(*,*) 'z = ',z
      write(*,'(A)') 'massofa not defined yet.'
      write(*,*) 'err: massofa problem'
      stop 1
    endif
    massofa = MASSOFATOM(z)
  end function massofa

  function crys_density(s)
    type(supercell) :: s
    real(double) :: mass,crys_density
    integer :: i

    mass = 0.0d0
    do i = 1, s%n
      mass = mass + massofatom(s%at(i)%z)
    enddo
    crys_density = (mass*AMU/(det3(s%a(1:3,1:3))*1.d-30))*(1.0d3/1.0d6)
  end function crys_density

  subroutine dbl_sort(n,v,ind,kflag)
    integer,parameter :: maxrecur=50
    integer :: i,n,segment,tind,ind(n),eindex(maxrecur),bindex(maxrecur),bi,ei,kflag,rightindex,leftindex
    real(double) :: v(n),vref,t
    logical :: foundright,foundleft,cross,inseq

    if(kflag /= 1 .and. kflag /= -1) then
      write(*,*) 'kflag in dbl_sort = ',kflag
      write(*,'(A)') 'err: kflag should be 1 (ascending) or -1 (descending).'
      stop 1
    endif

    if(n < 1) then
      write(*,*) 'dbl_sort, n = ',n
      write(*,'(A)') 'err: wrong array size for sorting.'
      stop 1
    endif
    if(n == 1) then
      return
    endif

    inseq = .true.
    i = 2
    do
      if(i > n) then
        exit
      endif
      if(.not. inseq) then
        exit
      endif
      if(kflag == 1) then
        inseq = v(i) >= v(i-1)
      else if(kflag == -1) then
        inseq = v(i) <= v(i-1)
      endif
      i = i + 1
    enddo
    if(inseq) then
     return
    endif

    if(kflag == -1) then
      do i = 1, n
        v(i) = -v(i)
      enddo
    endif

    segment = 1
    bindex(1) = 1
    eindex(1) = n
    do

      if(segment == 0) exit

      bi = bindex(segment)
      ei = eindex(segment)
      vref = v(bi)
      rightindex = ei
      leftindex = bi+1
      cross = .false.
      do
        if(cross) then
          exit
        endif
        foundright = .false.
        do
          if(foundright) exit
          if(v(rightindex) >= vref) then
            rightindex = rightindex - 1
            if(rightindex < bi) then
              foundright = .true.
            endif
          else
            foundright = .true.
          endif
        enddo
        foundleft = .false.
        do
          if(foundleft) exit
          if(v(leftindex) <= vref) then
            leftindex = leftindex + 1
            if(leftindex > ei) then
              foundleft = .true.
            endif
          else
            foundleft = .true.
          endif
        enddo
        if(leftindex > rightindex) then
          cross = .true.

          if(rightindex < bi) then

            segment = segment - 1

            if(ei-bi > 1) then
              segment = segment + 1
              bindex(segment) = bi+1
              eindex(segment) = ei
            endif
          else

            if(v(rightindex) < vref) then
              v(bi) = v(rightindex)
              v(rightindex) = vref
              tind = ind(bi)
              ind(bi) = ind(rightindex)
              ind(rightindex) = tind
            endif

            segment = segment - 1

            if(rightindex - bi > 1) then
              segment = segment + 1
              bindex(segment) = bi
              eindex(segment) = rightindex - 1
            endif
            if(ei-rightindex > 1) then
              segment = segment + 1
              if(segment > maxrecur) then
                write(*,*) 'maxrecur = ',maxrecur
                write(*,'(A)') 'maxrecur is not large enough'
              endif
              bindex(segment) = rightindex + 1
              eindex(segment) = ei
            endif
          endif
        else
          t = v(leftindex)
          v(leftindex) = v(rightindex)
          v(rightindex) = t
          tind = ind(leftindex)
          ind(leftindex) = ind(rightindex)
          ind(rightindex) = tind
          leftindex = leftindex+1
          rightindex = rightindex-1
        endif
      enddo
    enddo
    do i = 2, n
      if(v(i) < v(i-1)) then
        write(*,*) 'v(i),v(i-1)=',v(i),v(i-1)
        write(*,'(A)') 'err: v(i) must be greater or equal to v(i-1).'
        stop 1
      endif
    enddo
    if(kflag == -1) then
      do i = 1, n
        v(i) = -v(i)
      enddo
    endif
  end subroutine dbl_sort

  subroutine tryput(z,x,s,newatom)
    integer :: z,i,j,k,m
    real(double) :: dis,x(3),y(3),r1(3),r2(3),f(3),q(3)
    type(supercell) :: s
    logical :: found,newatom

    call abs2frac(s%b,x,y)
    do i = 1, 3
      call modulo1(y(i))
    enddo

    call frac2abs(s%a,y,r1)
    found = .false.
    do i = 1, s%n
      if(found) exit
      f = s%at(i)%f

      do j = -1, 1
        do k = -1, 1
          do m = -1, 1
            q = (/j,k,m/) + f
            call frac2abs(s%a,q,r2)
            dis = vecmag3(r1-r2)
            if(dis < dis_tol_tryput) then
              found = .true.
              if(s%at(i)%z /= z) then
                write(*,*) 'i,s%at(i)%z=',i,s%at(i)%z
                write(*,*) 'z = ',z
                write(*,'(A)') 'err: problem. The identity of atoms are not the same.'
                stop 1
              endif
            endif
          enddo
        enddo
      enddo
    enddo
    if(.not. found) then
      s%n = s%n+1
      s%at(s%n)%f = y
      s%at(s%n)%z = z
    endif
    newatom = .not. found
  end subroutine tryput

  subroutine read_struc(inputu,inputfile,filestem,inputformat,s)
    integer :: inputu,length,ind
    character(len=*) :: inputfile,filestem,inputformat
    type(supercell) :: s

    length = len_trim(inputfile)
    ind = index(trim(inputfile),'.',BACK=.TRUE.)
    filestem=inputfile(1:ind-1)

    inputformat=inputfile(ind+1:length)

    open(inputu,file=trim(inputfile),status='old',action='read')

    if(trim(inputformat)=='vasp' .or. trim(inputformat) == 'VASP') then
      call assign_str_fr_poscar(inputu,s)

    else if(trim(inputformat)=='STRUCT_OUT') then
      call assign_str_fr_struct_out(inputu,s)
    else if(trim(inputformat)=='fdf') then
      call assign_str_fr_fdf(inputu,s)
    else if(trim(inputformat)=='xsf') then
      call assign_str_fr_xsf(inputu,s)
    else
      write(*,*)
      write(*,*) 'WARNING: accepted format are vasp,arc/car,STRUCT_OUT,gjf,fdf'
      write(*,*) 'but inputformat is ',trim(inputformat)
      write(*,*) 'unrecognized input file format.'
      write(*,'(A)') 'err: check input file format.'
      stop 1
    endif
    close(inputu)

  end subroutine read_struc

  subroutine read_pg_operations(qtablepath,pgu,sg)
    character(len=*) :: qtablepath
    type(onesg) :: sg
    integer :: pgu
    logical :: found
    character(len=DSL) :: tmp1,tmp2
    integer :: s,nop
    character(len=DSL) :: pgkeyword,numoperkeyword,pglabel
    integer :: i,j,iconfirm

    open(unit=pgu,file=trim(qtablepath)//'/table-pg',status='old',action='read')

    read(pgu,*)
    found = .false.
    do
      if(found) exit
      read(pgu,DSF) tmp1
      if(trim(tmp1) /= '') then
        s = scan(tmp1,' ')
        pgkeyword = tmp1(1:s-1)
        tmp2 = tmp1(s+1:DSL)
        s = scan(tmp2,' ')
        pglabel = tmp2(1:s-1)

      endif
      if(trim(pgkeyword)=='END') then
        write(*,*) 'PG under consideration is '//trim(sg%pg)
        write(*,'(A)') 'err: cannot find the PG label from the table'
        stop 1
      endif
      if(trim(pglabel)==trim(sg%pg)) then
        found=.true.
      else

      endif
    enddo
    if(.not. found) then
      write(*,*) 'trim(sg%pg) is '//trim(sg%pg)
      write(*,'(A)') 'err: not found.'
      stop 1
    else

    endif

    read(pgu,DSF) tmp1
    read(tmp1,*) numoperkeyword, nop
    if(trim(numoperkeyword) /= 'NumOper') then
      write(*,*) 'numoperkeyword is '//trim(numoperkeyword)
      write(*,'(A)') 'err: NumOper must be found.'
      stop 1
    endif
    write(*,'(A)') 'Point group '//trim(sg%pg)//' has '//trim(N2str(nop))//' operations'
    sg%PGNumOper = nop
    do i = 1, nop
      read(pgu,*) iconfirm
      if(iconfirm/= i) then
        write(*,*) 'iconfirm,i=',iconfirm,i
        write(*,'(A)') 'err: iconfirm is not equal to i'
        stop 1
      else

      endif
      sg%op(1,1:4,i) = (/one,zero,zero,zero/)
      do j = 2, 4
        read(pgu,*) sg%op(j,2:4,i)
      enddo
    enddo
    close(pgu)
  end subroutine read_pg_operations

  function pg_id_map(d,t)
    integer :: pg_id_map,d,t

    pg_id_map = d*pg_uuid_L4 + abs(t)*pg_uuid_L3 + t*pg_uuid_L2 + t*d
  end function pg_id_map

  subroutine pg_name_init()
    PG_NAME(1) =  trim("C1")
    PG_NAME(2) =  trim("Ci")
    PG_NAME(3) =  trim("C2")
    PG_NAME(4) =  trim("Cs")
    PG_NAME(5) =  trim("C2h")
    PG_NAME(6) =  trim("D2")
    PG_NAME(7) =  trim("C2v")
    PG_NAME(8) =  trim("D2h")
    PG_NAME(9) =  trim("C4")
    PG_NAME(10) = trim("S4")
    PG_NAME(11) = trim("C4h")
    PG_NAME(12) = trim("D4")
    PG_NAME(13) = trim("C4v")
    PG_NAME(14) = trim("D2d")
    PG_NAME(15) = trim("D4h")
    PG_NAME(16) = trim("C3")
    PG_NAME(17) = trim("S6")
    PG_NAME(18) = trim("D3")
    PG_NAME(19) = trim("C3v")
    PG_NAME(20) = trim("D3d")
    PG_NAME(21) = trim("C6")
    PG_NAME(22) = trim("C3h")
    PG_NAME(23) = trim("C6h")
    PG_NAME(24) = trim("D6")
    PG_NAME(25) = trim("C6v")
    PG_NAME(26) = trim("D3h-type1")
    PG_NAME(27) = trim("D6h")
    PG_NAME(28) = trim("T")
    PG_NAME(29) = trim("Th")
    PG_NAME(30) = trim("O")
    PG_NAME(31) = trim("Td")
    PG_NAME(32) = trim("Oh")
  end subroutine pg_name_init

  function get_pg_index(pg_id)
    integer :: pg_id,get_pg_index
    logical :: found
    integer :: i,ind

    found = .false.
    ind = -1
    do i = 1, 32
      if(.not. found .and. pg_id == PG_UUID(i)) then
        found = .true.
        ind = i
      endif
    enddo
    if(.not. found) then
      write(*,'(A)') 'err: get_pg_index: UUID not found.'
      stop 1
    endif
    get_pg_index = ind
  end function get_pg_index

  subroutine check_SG_closure(sg)
    type(onesg) :: sg
    integer :: nop
    real(double),allocatable :: op(:,:,:),origmat3(:,:,:),trace(:),pg_det3(:)
    integer :: i,j,k
    real(double) :: tmpmat(4,4),p(3),diffmat(4,4),diffmat3(3,3),mat1(3,3),mat2(3,3),mat3(3,3)
    real(double),parameter :: eps = 1.0d-8
    integer,allocatable :: sudoku(:,:),classes(:)
    logical :: found
    real(double) :: matmag,sumj, targetsum
    real(double),allocatable :: sortv(:)
    real(double) :: RG,mattrace,detm
    integer,allocatable :: ig(:),indarr(:)
    integer :: NG,abs_op,det_int,trace_int
    integer :: pg_id,pg_index
    logical :: print_op
    real(double) :: diffmag

    write(*,*)
    write(*,*) 'Checking space group information.'
    write(*,*) 'Name of SG is '//trim(sg%sglabel)
    write(*,*) 'Name of PG is '//trim(sg%pg)
    nop = sg%nop
    write(*,*) 'Number of operations= ',trim(N2str(nop))

    allocate(op(4,4,nop))

    print_op = .false.
    if(print_op) then
      write(*,'(A)') 'In check_SG_closure (for cubic SGs, the order might have'
      write(*,'(A)') 'been ordered according to conjugacy classes):'
    endif

    do i = 1, nop
      op(1:4,1:4,i) = sg%op(1:4,1:4,i)

      if(print_op) then
        write(*,*) 'i = ', i

        diffmag = abs(op(1,1,i)-one) + abs(op(1,2,i)) + abs(op(1,3,i)) + abs(op(1,4,i))
        if(diffmag > 1.0d-16) then
          write(*,'(A)') 'err: first row must be (1,0,0,0)'
          stop 1
        endif
        write(*,'(4F8.3)') op(2,1:4,i)
        write(*,'(4F8.3)') op(3,1:4,i)
        write(*,'(4F8.3)') op(4,1:4,i)
      endif
    enddo

    allocate(sudoku(nop,nop))

    do i = 1, nop
      do j = 1, nop
        call matmatn(4,op(1,1,i),op(1,1,j),tmpmat(1,1))

        p(1:3) = tmpmat(2:4,1)
        call foldtofirstzone3(p(1))
        do k = 1, 3
          if(p(k) > one - eps) then
            write(*,*) 'p(k) = ', p(k)
            p(k) = zero
            write(*,*) 'reset p(k), p(k) = ',p(k)
            write(*,'(A)') 'err: does this really happen?'
            stop 1
          endif
        enddo

        tmpmat(2:4,1) = p(1:3)

        found = .false.
        do k = 1, nop
          diffmat(1:4,1:4) = sg%op(1:4,1:4,k) - tmpmat(1:4,1:4)
          matmag = matmagn(4,diffmat(1,1))
          if(matmag < eps) then
            found = .true.
            sudoku(i,j) = k
          endif
        enddo
        if( .not. found) then
          write(*,*) 'i, j= ',i,j
          write(*,*) 'tmpmat(1:4,1:4)= '
          write(*,'(4F8.3)') tmpmat(1,1:4)
          write(*,'(4F8.3)') tmpmat(2,1:4)
          write(*,'(4F8.3)') tmpmat(3,1:4)
          write(*,'(4F8.3)') tmpmat(4,1:4)
          write(*,*) 'closure test failed.'
          write(*,'(A)') 'err: closure fail.'
          stop 1
        endif
      enddo
    enddo
    write(*,*)

    do i = 1, nop

    enddo
    targetsum = nop*(nop+one)/two

    do i = 1, nop
      sumj = zero
      do j = 1, nop
        sumj = sumj + sudoku(i,j)
      enddo
      if(abs(sumj - targetsum) > eps) then
        write(*,*) 'failed: sumj, targetsum= ',sumj, targetsum
        write(*,'(A)') 'err: matrix multiplication table sudoku property failed.'
        stop 1
      endif
    enddo

    do i = 1, nop
      sumj = zero
      do j = 1, nop
        sumj = sumj + sudoku(j,i)
      enddo
      if(abs(sumj - targetsum) > eps) then
        write(*,*) 'sumj, targetsum= ',sumj, targetsum
        write(*,'(A)') 'err: sudoku failed.'
        stop 1
      endif
    enddo
    write(*,'(A)') 'Passed closure property test'

    if(nop /= sg%PGNumOper) then

      write(*,*) 'WARNING: We will use the exact number of operations for the PG.'
      write(*,'(A,I4,/,A,I4)') 'We change the number of operations in the SG= ',nop,' to number of operations in the PG= ',sg%PGNumOper
    endif
    nop = sg%PGNumOper

    allocate(trace(nop))
    allocate(pg_det3(nop))

    do i = 1, nop
      sumj = zero
      do j = 2, 4
        sumj = sumj + op(j,j,i)
      enddo
      trace(i) = sumj

    enddo

    allocate(classes(nop))
    allocate(origmat3(3,3,nop))
    do i = 1, nop
      origmat3(1:3,1:3,i) = op(2:4,2:4,i)
      pg_det3(i) = det3(origmat3(1,1,i))
    enddo

    classes(1:nop) = 0
    do i = 1, nop
      if(classes(i) /= 0) cycle
      mat1(1:3,1:3) = origmat3(1:3,1:3,i)

      do j = 1, nop

        mat2 = origmat3(1:3,1:3,j)

        mat3 = BAinvB(mat2(1,1),mat1(1,1))

        found = .false.
        do k = 1, nop
          diffmat3(1:3,1:3) = mat3(1:3,1:3)-origmat3(1:3,1:3,k)

          if(.not. found .and. matmagn(3,diffmat3(1,1)) < eps) then
            found = .true.
            if(classes(k) /= 0) then
              if(classes(k) /= i) then
                write(*,*) 'i,classes(k)= ',i,classes(k)
                write(*,'(A)') 'err: inconsistency.'
                stop 1
               endif
            else
              classes(k) = i
            endif
          endif
        enddo
        if( .not. found) then
          write(*,'(A)') 'err: cannot find companion in equivalent classes.'
          stop 1
        endif
      enddo
    enddo
    if(nop < 1) then
      write(*,*) 'nop = ',nop
      write(*,'(A)') 'err: going to hit allocation problem since the array size is zero.'
      stop 1
    endif
    allocate(sortv(nop))
    allocate(indarr(nop))
    do i = 1, nop
      sortv(i) = classes(i)
      indarr(i) = i
    enddo
    call dbl_sort(nop,sortv(1),indarr(1),1)

    do i = 1, nop

    enddo
    do i = 1, nop
      if(i == indarr(i)) then

      else

        write(*,'(A)') ' We expect the operations have been fully ordered.'
        write(*,'(A)') 'err: Search the keyword OPERATION_MAP in commod.f90.'
        stop 1
      endif
    enddo
    write(*,'(A)') 'Passed the natural ordering test for equivalent classes.'

    allocate(ig(nop+1))
    RG = sortv(1)
    NG = 1
    IG(1) = 1
    do k = 2, nop
      if(abs(sortv(k)-RG) > 1.0d-4) then
        RG = sortv(k)
        NG = NG + 1
        IG(NG) = k
      endif
    enddo
    ig(ng+1) = nop+1

    write(*,*)
    write(*,'(A)') 'Equivalence classes of the point group: '//trim(sg%pg)

    pg_id = ng*PG_UUID_L5
    write(*,'(A)') '--------- classes begin -----'
    do i = 1, ng
      abs_op = indarr(ig(i))
      detm = det3( origmat3(1,1, abs_op ))
      mattrace = tracen(3, origmat3(1,1,abs_op))
      write(*,'(A,I4,I4,2F8.4)') 'cluster, size of cluster, det, trace=',i,ig(i+1)-ig(i),detm,mattrace
      det_int = nint(detm)
      trace_int = nint(mattrace)
      call report_operation_type(det_int,trace_int)
      pg_id  = pg_id + pg_id_map(det_int,trace_int)
    enddo
    write(*,'(A)') '--------- classes end -----'
    pg_index = get_pg_index(pg_id)
    deallocate(op)
    deallocate(sudoku)
    deallocate(trace)
    deallocate(pg_det3)
  end subroutine check_SG_closure

  subroutine read_sg(qtablepath,sg_label,sg,pgu,tu,dtu,sgtable,dsgtable)
    type(onesg) :: sg
    character(len=*) :: sg_label,qtablepath
    character(len=*) :: sgtable,dsgtable
    integer :: pgu,tu,dtu
    character(len=DSL) :: base
    character(len=DSL) :: lpart
    character(len=DSL) :: ext
    integer :: str_index
    integer :: strlen,ind,newnc,newnop
    integer :: readmode

    strlen = len(trim(sg_label))
    ind = scan(sg_label,':')

    if(ind > 0) then
      write(*,*) 'sg_label contains :'
      lpart = sg_label(1:ind-1)
      ext = sg_label(ind+1:strlen)
      write(*,*) 'lpart is '//trim(lpart)
      write(*,*) 'ext is '//trim(ext)

      str_index = index(trim(sg_label),'basic')
      if(str_index > 0) then
        write(*,*) 'found basic in the sglabel'
        base = sg_label(1:str_index-2)
        write(*,*) 'base is '//trim(base)
        call sg_table_read(qtablepath,pgu,tu,sgtable,base,sg,readmode)
        call check_SG_closure(sg)

        newnc = 1
        newnop = sg%PGNumOper
        write(*,'(/A)') 'With : case, since the label has "-basic", we set ncentering from '//trim(N2str(sg%ncentering))//' to '//trim(N2str(newnc))//', and nop from '//trim(N2str(sg%nop))//' to '//trim(N2str(newnop))
        sg%ncentering = newnc
        sg%nop = newnop

      else
        write(*,*) 'no basic in the sglabel'

        base = trim(lpart)
        call sg_table_read(qtablepath,pgu,tu,sgtable,base,sg,readmode)
        call check_SG_closure(sg)
      endif

      call update_derivsg_table(qtablepath,dtu,dsgtable,sg_label,sg,lpart,ext)

      call check_SG_closure(sg)
    else
      str_index = index(trim(sg_label),'basic')
      if(str_index > 0) then
        write(*,*) 'found basic in the sglabel'
        base = sg_label(1:str_index-2)
        write(*,*) 'base is '//trim(base)
        call sg_table_read(qtablepath,pgu,tu,sgtable,base,sg,readmode)
        if(readmode == 1 .or. readmode == 2) then

        else
          write(*,*) 'readmode =',readmode
          write(*,'(A)') 'err: impossible.'
          stop 1
        endif

        call check_SG_closure(sg)

        write(*,'(/A)') 'Without : case, since the label has "-basic", we set ncentering to '//trim(N2str(sg%ncentering))//' and nop to '//trim(N2str(sg%nop))
        sg%ncentering = 1
        sg%nop = sg%PGNumOper

      else

        base = trim(sg_label)
        call sg_table_read(qtablepath,pgu,tu,sgtable,base,sg,readmode)
        if(readmode == 1 .or. readmode == 2) then
          call check_SG_closure(sg)
        else if(readmode == 3) then
          write(*,*) 'readmode=3: we do not check the SG_closure.'
        else
          write(*,*) 'readmode =',readmode
          write(*,'(A)') 'err: readmode problem.'
          stop 1
        endif
      endif
    endif
  end subroutine read_sg

  subroutine sg_table_read(qtablepath,pgu,tu,table,targetlabel,sg,readmode)
    integer :: pgu,tu
    character(len=*) :: table,qtablepath
    character(len=*) :: targetlabel
    type(onesg) :: sg
    logical :: found
    integer :: i,j,k,iconfirm,nop
    character(len=DSL) :: sglabel
    character(len=DSL) :: sgkeyword,numoperkeyword
    character(len=DSL) :: tmp1,tmp2,tmp3,str,key
    integer :: s,v,ncentering,nopread
    character(len=DSL) :: string(3)
    real(double) :: f,nf,df
    integer :: pru
    real(double) :: shift(3)
    real(double) :: tmpMat(4,4,48)
    character(len=DSL) :: sg_f90_entry,master_sg_entry
    integer :: readmode
    integer :: labellen

    readmode = -10000

    write(*,*)
    labellen = len(trim(targetlabel))
    if(labellen == 0) then
      write(*,*) 'Empty targetlabel.'
      write(*,'(A)') 'err: Empty targetlabel.'
      stop 1
    endif
    write(*,'(A)') 'sg_table_read: targetlabel is '//trim(targetlabel)

    sg%assigned = .FALSE.

    do i = 0, 230
      do j = 1, maxSGvariant
        if(trim(SGbase(j,i)%sglabel)==trim(targetlabel)) then
          sg%sglabel = trim( SGbase(j,i)%sglabel)
          sg%PG = trim(SGbase(j,i)%PG)
          sg%ncentering = SGbase(j,i)%ncentering
          do k = 1, sg%ncentering
            sg%centering(1:3,k) = SGbase(j,i)%centering(1:3,k)
          enddo

          call read_pg_operations(qtablepath,pgu,sg)

          if(sg%PGNumOper /= SGbase(j,i)%PGNumOper) then
            write(*,*) 'sg%PGNumOper,SGbase(j,i)%PGNumOper=',sg%PGNumOper,SGbase(j,i)%PGNumOper
            write(*,'(A)') 'err: serious inconsistent problem.'
            stop 1
          endif
          sg%symmorphic = SGbase(j,i)%symmorphic
          if(sg%symmorphic == 1) then
            do k = 1, sg%PGNumOper
              sg%nsy(1:3,k) = (/zero,zero,zero/)
              sg%op(2:4,1,k) = (/zero,zero,zero/)
            enddo
          else
            do k = 1, sg%PGNumOper
              sg%nsy(1:3,k) = SGbase(j,i)%nsy(1:3,k)
              sg%op(2:4,1,k) = SGbase(j,i)%nsy(1:3,k)
            enddo
          endif
          sg%nop = SGbase(j,i)%nop
          sg%reorder_IT_op_sequence = SGbase(j,i)%reorder_IT_op_sequence
          if(sg%reorder_IT_op_sequence == 1) then
            sg%identity_map(1:sg%PGNumOper) = SGbase(j,i)%identity_map(1:sg%PGNumOper)
            sg%operation_map(1:sg%PGNumOper) = SGbase(j,i)%operation_map(1:sg%PGNumOper)
          endif
          sg%assigned = .TRUE.
          readmode = 1
        endif
      enddo
    enddo

    if(sg%assigned) then

      nop = sg%nop
      ncentering = nop/sg%PGNumOper
      if(ncentering /= sg%ncentering) then
        write(*,*) 'ncentering,sg%ncentering = ',ncentering,sg%ncentering
        write(*,'(A)') 'err: inconsistency.'
        stop 1
      endif

      if(sg%reorder_IT_op_sequence == 1) then
        do i = 1, sg%PGNumOper
          tmpMat(1:4,1:4,i) = sg%op(1:4,1:4,i)
        enddo
        do i = 1, sg%PGNumOper
          sg%op(1:4,1:4,i) = tmpMat(1:4,1:4,   sg%operation_map(i))
        enddo
        write(*,'(A)') 'In mode 1: Operations are reordered for O,Td,Oh since IT SGs equivalence classes are not contiguous'
      endif

      do i = 2, ncentering
        do j = 1, sg%PGNumOper
          sg%op(1:4,1:4,(i-1)*sg%PGNumOper+j) = sg%op(1:4,1:4,j)

          sg%op(2:4,1,(i-1)*sg%PGNumOper+j) = sg%op(2:4,1,(i-1)*sg%PGNumOper+j) + sg%centering(1:3,i)

          shift(1:3) = sg%op(2:4,1,(i-1)*sg%PGNumOper+j)
          call foldtofirstzone3(shift(1))

          sg%op(2:4,1,(i-1)*sg%PGNumOper+j) = shift(1:3)

        enddo
      enddo

    else
      write(*,'(A)') 'Do not use a SG entry in commod.f90, instead we use '//trim(table)

      write(*,'(A)') ' -------- To open SG table file '//trim(qtablepath)//'/'//trim(table)
      open(unit=tu,file=trim(qtablepath)//'/'//trim(table),status='old',action='read')

      read(tu,*)

      found=.false.

      do
        if(found) exit

        read(tu,DSF) tmp1

        if(trim(tmp1) /= '') then

          s = scan(tmp1,' ')
          sgkeyword = tmp1(1:s-1)
          tmp2 = tmp1(s+1:DSL)
          s = scan(tmp2,' ')
          sglabel = tmp2(1:s-1)

        endif

        if(trim(sgkeyword)=='END') then
          exit
        endif

        if(trim(sglabel)==trim(targetlabel)) then
          found=.true.
        else

        endif
      enddo

      if(found) then
        write(*,'(A)') 'Found SG label= '//trim(targetlabel)
      else
        write(*,'(A)') 'SG under consideration is '//trim(targetlabel)
        write(*,'(A)') 'err: in sg_table_read, cannot find SG label.'
        stop 1
      endif

      sg%PGNumOper=-1

      sg%ncentering=-1
      sg%reorder_IT_op_sequence = 0

      do
        read(tu,DSF) tmp1

        read(tmp1,*) str

        if(trim(str) == '#optional') then
          s = scan(str,' ')
          tmp2 = tmp1(s+1:DSL)

          s = scan(tmp2,' ')
          key=tmp2(1:s)

          tmp3 = tmp2(s+1:DSL)

          if(trim(key) == 'PG') then
            readmode = 2
            write(*,'(A)') 'Use a PG label to assign the matrix operations.'
            write(*,'(A)') 'PG: to process '//trim(tmp3)
            s = scan(tmp3,' ')
            sg%pg = tmp3(1:s)
            write(*,'(A)') 'Point group '//trim(sg%pg)

          else if(trim(key) == 'ncentering') then

            read(tmp3,'(I10)') v
            write(*,*) 'ncentering = ',v
            if(v < 1 .or. v > 4) then
              write(*,*) 'v = ',v
              write(*,'(A)') 'err: invalid v.'
              stop 1
            endif
            sg%ncentering = v
            do i = 1, v
              read(tu,*) tmp1,tmp2,sg%centering(1:3,i)
              iconfirm = str2N(tmp2)
              if(i /= iconfirm) then
                write(*,*) 'centering: i,iconfirm=',i,iconfirm
                write(*,'(A)') 'err: i and iconfirm must be the same.'
                stop 1
              endif
              write(*,'(A,3F10.5,A)') 'tmp1 is '//trim(tmp1)//', tmp2 is '//trim(tmp2)//', ncentering= |', sg%centering(1:3,i),'|'
            enddo

          else if(trim(key) == 'PGNumOper') then

            call read_pg_operations(qtablepath,pgu,sg)

            read(tmp3,'(I10)') v
            write(*,*) 'PGNumOper = ',v
            if(sg%PGNumOper /= v) then
              write(*,*) 'sg%PGNumOper = ',sg%PGNumOper
              write(*,*) 'v = ',v
              write(*,'(A)') 'err: serious inconsistency problem.'
              stop 1
            endif

            read(tu,*) tmp1,tmp2 !
            write(*,'(A,A)') 'tmp1 = '//trim(tmp1),', tmp2 = '//trim(tmp2)
            if(trim(tmp2) == 'symmorphic' .or. trim(tmp2) == 'nonsymmorphic') then

            else
              write(*,'(A)') 'err: the keyword must be symmorphic and nonsymmorphic.'
              stop 1
            endif

            if(trim(tmp2) == 'symmorphic') then
              sg%symmorphic = 1

              do i = 1, v
                sg%op(2:4,1,i) = (/zero,zero,zero/)
              enddo
            else
              sg%symmorphic = 0
              write(*,*) 'v = ',v
              do i = 1, v
                read(tu,*) tmp1,tmp2,sg%op(2:4,1,i)
                iconfirm = str2N(tmp2)
                if(i /= iconfirm) then
                  write(*,*) 'fractional translation: i,iconfirm=',i,iconfirm
                  write(*,'(A)') 'err: i and iconfirm must be the same.'
                  stop 1
                endif
                if(i > 48) then
                  write(*,'(A)') 'err: bound exceeded.'
                  stop 1
                endif
                sg%nsy(1:3,i) = sg%op(2:4,1,i)
                write(*,'(A,3F10.5)') 'tmp1 is '//trim(tmp1)//', tmp2 is '//trim(tmp2)//', translation vector= ', sg%op(2:4,1,i)
              enddo
            endif

          else if(trim(key) == 'reorder_IT_op_sequence') then
            if(sg%PGNumOper < 1) then
              write(*,'(A)') 'err: bad situation.'
              stop 1
            endif
            sg%reorder_IT_op_sequence = 1
            read(tu,*) tmp1, tmp2, sg%identity_map(1:sg%PGNumOper)
            if(trim(tmp1) /= '#optional' .or. trim(tmp2) /= 'identity_map') then
              write(*,*) 'tmp1= '//trim(tmp1),', tmp2= ',trim(tmp2)
              write(*,'(A)') 'err: keyword mismatched.'
              stop 1
            endif
            read(tu,*) tmp1, tmp2, sg%operation_map(1:sg%PGNumOper)
            if(trim(tmp1) /= '#optional' .or. trim(tmp2) /= 'operation_map') then
              write(*,*) 'tmp1= '//trim(tmp1),', tmp2= ',trim(tmp2)
              write(*,'(A)') 'err: keyword mismatched.'
              stop 1
            endif
          else
            write(*,*) 'keyword of |'//trim(key)//'| will not be processed.'
            write(*,'(A)') 'err: keyword not found.'
            stop 1
          endif
        else

          exit
        endif
      enddo

      if(readmode == 2) then
        if(sg%ncentering == -1 .or. sg%PGNumOper == -1) then
          write(*,*) 'sg%ncentering = ',sg%ncentering
          write(*,*) 'sg%PGNumOper = ',sg%PGNumOper
          write(*,'(A)') 'err: Serious inconsistency input.'
          stop 1
        endif
      endif

      if(sg%reorder_IT_op_sequence == 1) then
        do i = 1, sg%PGNumOper
          tmpMat(1:4,1:4,i) = sg%op(1:4,1:4,i)
        enddo
        do i = 1, sg%PGNumOper
          sg%op(1:4,1:4,i) = tmpMat(1:4,1:4,   sg%operation_map(i))
        enddo
        write(*,'(A)') 'In mode 2: Operations are reordered for O,Td,Oh since IT SGs equivalence classes are not contiguous'
      endif

      read(tmp1,*) numoperkeyword, nop

      if(trim(numoperkeyword) /= 'NumOper') then
        write(*,*) 'numoperkeyword is '//trim(numoperkeyword)
        write(*,'(A)') 'err: NumOper must be found.'
        stop 1
      endif
      write(*,'(A)') 'space group of '//trim(targetlabel)//' has '//trim(N2str(nop))//' operations'

      sg%sglabel = trim(targetlabel)
      sg%nop = nop

      if(readmode == 2) then

        if(mod(nop, sg%PGNumOper) /= 0) then
          write(*,'(A)') 'err: PGNumOper does not divide nop.'
          stop 1
        endif
        ncentering = nop/sg%PGNumOper
        if(ncentering /= sg%ncentering) then
          write(*,*) 'ncentering,sg%ncentering = ',ncentering,sg%ncentering
          write(*,'(A)') 'err: Inconsistency.'
          stop 1
        endif

        do i = 2, ncentering
          do j = 1, sg%PGNumOper
            sg%op(1:4,1:4,(i-1)*sg%PGNumOper+j) = sg%op(1:4,1:4,j)
            sg%op(2:4,1,(i-1)*sg%PGNumOper+j) = sg%op(2:4,1,(i-1)*sg%PGNumOper+j) + sg%centering(1:3,i)

            shift(1:3) = sg%op(2:4,1,(i-1)*sg%PGNumOper+j)
            call foldtofirstzone3(shift(1))
            sg%op(2:4,1,(i-1)*sg%PGNumOper+j) = shift(1:3)
          enddo
        enddo

      else
        readmode = 3
        write(*,*) 'This will be defaulted to reading operation by operation from the table'
        nopread = nop

        do i = 1, nopread
          read(tu,*) iconfirm
          if(iconfirm/= i) then
            write(*,*) 'iconfirm,i=',iconfirm,i
            write(*,'(A)') 'err: iconfirm is not equal to i'
            stop 1
          else

          endif
          sg%op(1,1:4,i) = (/one,zero,zero,zero/)
          do j = 2, 4
            read(tu,*) sg%op(j,1:4,i)
          enddo
        enddo

        return
      endif
      close(tu)
      write(*,'(A)') ' ------- Closing SG table file '//trim(table)
    endif

    if(readmode == 1 .or. readmode == 2) then
    else
      write(*,*) 'readmode = ',readmode
      write(*,'(A)') 'err: bad readmode.'
      stop 1
    endif

    write(*,*)
    sg_f90_entry='echo-sg-f90-entry.dat'

    pru = newunit()
    open(unit=pru,file=trim(sg_f90_entry),status='replace')

    write(pru,'(A)') '--------------------------------------'

    sglabel = trim(sg%sglabel)
    i = scan(sglabel,'-')
    write(pru,'(A)') '    SG = '//sglabel(1:i-1)
    write(pru,'(A)') '    Y = 1'
    write(pru,'(A)') '    SGbase(Y,SG)%sglabel    ='''//trim(sg%sglabel)//''''
    write(pru,'(A)') '    SGbase(Y,SG)%PG         ='''//trim(sg%pg)//''''
    write(pru,'(A)') '    SGbase(Y,SG)%ncentering ='//trim(N2str(sg%ncentering))
    do i = 1, sg%ncentering
      do j = 1, 3
        f = sg%centering(j,i)

        call frac2nd(f,nf,df)
        if(abs(df-one) < 1.0d-6) then

          string(j)=trim(N2str(nint(nf)))//'.0d0'
        else
          string(j)=trim(N2str(nint(nf)))//'.0d0/'//trim(N2str(nint(df)))//'.0d0'
        endif
      enddo
      write(pru,'(A,A16,A,A16,A,A16,A)') '    SGbase(Y,SG)%centering(1:3,'//trim(N2str(i))//') = (/',trim(string(1)),',',trim(string(2)),',',trim(string(3)),' /)'
    enddo
    write(pru,'(A)') '    SGbase(Y,SG)%PGNumOper  ='//trim(N2str(sg%PGNumOper))
    write(pru,'(A)') '    SGbase(Y,SG)%symmorphic ='//trim(N2str(sg%symmorphic))
    if(sg%symmorphic == 0) then
      do i = 1, sg%PGNumOper
        do j = 1, 3
          f = sg%nsy(j,i)

          call frac2nd(f,nf,df)
          if(abs(df-one) < 1.0d-6) then

            string(j)=trim(N2str(nint(nf)))//'.0d0'
          else
            string(j)=trim(N2str(nint(nf)))//'.0d0/'//trim(N2str(nint(df)))//'.0d0'
          endif
        enddo
        write(pru,'(A,A16,A,A16,A,A16,A)') '    SGbase(Y,SG)%nsy(1:3,'//trim(N2str(i))//') = (/',trim(string(1)),',',trim(string(2)),',',trim(string(3)),' /)'
      enddo
    endif

    if(sg%reorder_IT_op_sequence == 1) then
      write(pru,'(A)') '    SGbase(Y,SG)%reorder_IT_op_sequence = 1'
      tmp1=trim(N2str(sg%identity_map(1)))
      do i = 2, sg%PGNumOper
        tmp1=trim(tmp1)//','//trim(N2str(sg%identity_map(i)))
      enddo
      write(pru,'(A)') '    SGbase(Y,SG)%identity_map(1:'//trim(N2str(sg%PGNumOper))//') = (/'//trim(tmp1)//'/)'
      tmp1=trim(N2str(sg%operation_map(1)))
      do i = 2, sg%PGNumOper
        tmp1=trim(tmp1)//','//trim(N2str(sg%operation_map(i)))
      enddo
      write(pru,'(A)') '    SGbase(Y,SG)%operation_map(1:'//trim(N2str(sg%PGNumOper))//') = (/'//trim(tmp1)//'/)'
    else

    endif

    write(pru,'(A)') '    SGbase(Y,SG)%nop        ='//trim(N2str(sg%nop))
    write(pru,'(A)') '    SGbase(Y,SG)%assigned   =.TRUE.'
    write(pru,*)
    close(pru)

    master_sg_entry='echo-master-sg-entry.dat'

    pru = newunit()
    open(unit=pru,file=trim(master_sg_entry),status='replace')
    write(pru,'(A)') 'SG '//trim(sg%sglabel)
    write(pru,'(A)') '#optional PG '//trim(sg%pg)
    write(pru,'(A)') '#optional ncentering '//trim(N2str(sg%ncentering))
    do i = 1, sg%ncentering
      write(pru,'(A,3F20.15)') '#optional '//trim(N2str(i)),sg%centering(1:3,i)
    enddo
    write(pru,'(A)') '#optional PGNumOper '//trim(N2str(sg%PGNumOper))
    if(sg%symmorphic == 1) then
      write(pru,'(A)') '#optional symmorphic'
    else if(sg%symmorphic == 0) then
      write(pru,'(A)') '#optional nonsymmorphic'
      do i = 1, sg%PGNumOper
        write(pru,'(A,3F8.4)') '#optional '//trim(N2str(i)),sg%nsy(1:3,i)
      enddo
    else
      write(*,'(A)') 'err: problem in symmorphic.'
      stop 1
    endif
    if(sg%reorder_IT_op_sequence == 1) then
      write(pru,'(A)') '#optional reorder_IT_op_sequence'
      write(tmp2,'(100I3)') sg%identity_map(1:sg%PGNumOper)
      write(pru,'(A)') '#optional  identity_map '//trim(tmp2)
      write(tmp2,'(100I3)') sg%operation_map(1:sg%PGNumOper)
      write(pru,'(A)') '#optional operation_map '//trim(tmp2)
    else

    endif
    write(pru,'(A)') 'NumOper '//trim(N2str(sg%PGNumOper*sg%ncentering))

    close(pru)

  end subroutine sg_table_read

  subroutine update_derivsg_table(qtablepath,dtu,dtable,targetlabel,sg,lpart,ext)
    character(len=*) :: qtablepath
    integer :: dtu
    character(len=*) :: dtable
    character(len=*) :: targetlabel
    type(onesg) :: sg
    logical :: found
    character(len=DSL) :: line,tmp1,tmp2,sgkeyword,sglabel,ext,lpart
    integer :: s
    real(double) :: invT(3,3),blk4x4(4,4),newblk4x4(4,4),t1(3),R1(3,3),t2(3),R2(3,3)
    integer :: i,nop

    real(double) :: t3(3)
    real(double) :: origvec(3)
    integer :: j,pru
    real(double) :: f,nf,df
    character(len=DSL) :: string(3)
    character(len=DSL) :: filename
    integer :: ind,SGn
    character :: FI

    write(*,*)
    write(*,*)
    write(*,*) 'in update_derivsg_table.'

    if(len(trim(targetlabel)) == 0) then
      write(*,*) 'Empty targetlabel.'
      stop
    else
      write(*,*) 'targetlabel is '//trim(targetlabel)
    endif

    if(trim(ext) == 'stdbox' .or. trim(ext) == 'mat2' .or. trim(ext) == 'mat6') then
      write(*,*) 'lpart is '//trim(lpart)
      ind = scan(lpart,'-')
      tmp1 = lpart(1:ind-1)
      read(tmp1,'(I10)') SGn
      FI = lpart(ind+1:ind+1)
      write(*,*) 'SGn,FI = ',SGn,FI

      if(SGn >= 1 .and. SGn <= 2) then
        write(*,*) 'Triclinic case does not a chance to use stdbox'
        write(*,'(A)') 'triclinic'
      else if(SGn >= 3 .and. SGn <= 15) then
        if(trim(ext) /= 'stdbox') then
          write(*,*) 'For monoclinic, we must have stdbox case'
          write(*,'(A)') 'err: not stdbox.'
          stop 1
        endif

        s = index(lpart,'-')
        FI = lpart(s+1:s+1)
        if(FI == 'C') then
          sg%T(1:3,1:3) = MonocBaxisBaseCentredCmat(1:3,1:3)
        else
           write(*,'(A)') 'er: crash, P should call this subroutine.'
        endif
      else if(SGn >= 16 .and. SGn <= 74) then
        if(trim(ext) /= 'stdbox') then
          write(*,'(A)') 'err: not stdbox.'
          stop 1
        endif
        FI = lpart(4:4)
        if(FI == 'C') then
          sg%T(1:3,1:3) = OBaseCmat1(1:3,1:3)
        else if (FI == 'F') then
          sg%T(1:3,1:3) = OFCmat(1:3,1:3)
        else if (FI == 'I') then
          sg%T(1:3,1:3) = OBCmat(1:3,1:3)
        else if (FI == 'A') then
          sg%T(1:3,1:3) = OBaseAmat(1:3,1:3)
        else
          write(*,*) 'FI  is '//FI
          write(*,'(A)') 'err: It must be C, F, I, or A for orthorhombic'
          stop 1
        endif
      else if(SGn >= 75 .and. SGn <= 142) then
        if(trim(ext) /= 'stdbox') then
          write(*,'(A)') 'err: not stdbox.'
          stop 1
        endif
        if(FI == 'I') then
          sg%T(1:3,1:3) = BCTmat(1:3,1:3)
        else
          write(*,*) 'Bad SGn = ',SGn
          write(*,'(A)') 'err: crash.. FC should be I only.'
          stop 1
        endif
      else if(SGn >= 143 .and. SGn <= 167) then
        if(trim(ext) == 'stdbox') then
          write(*,'(A)') 'err: trigonal cannot be stdbox. It must be mat6 (Y) or mat2 (inverted Y)'
          stop 1
        endif
        write(*,*) 'triginal: SGn= ',SGn
        if(FI == 'R' .and. trim(ext) == 'mat2') then
          sg%T(1:3,1:3) = Mat2Trigonal(1:3,1:3)
        else if(FI == 'R' .and. trim(ext) == 'mat6') then
          sg%T(1:3,1:3) = Mat6Trigonal(1:3,1:3)
        else
          write(*,'(A)') 'err: crash, should be R only.'
          stop 1
        endif
      else if(SGn >= 168 .and. SGn <= 194) then
        write(*,*) 'Hexagonal case does not a chance to use stdbox'
        write(*,'(A)') 'Hexagonal'
      else if(SGn >= 195 .and. SGn <= 230) then
        write(*,*) 'cubic: SGn= ',SGn
        if(FI == 'F') then
          sg%T(1:3,1:3) = FCCmat(1:3,1:3)
        else if(FI == 'I') then
          sg%T(1:3,1:3) = BCCmat(1:3,1:3)
        else
          write(*,'(A)') 'err: crash.. should be F or I only.'
          stop 1
        endif
      else
         write(*,*) 'Bad SGn = ',SGn
         write(*,'(A)') 'err: SGn problem'
         stop 1
      endif
      sg%origvec(1:3) = zero
    else
      open(unit=dtu,file=trim(qtablepath)//'/'//trim(dtable),status='old',action='read')
      read(dtu,*)
      found=.false.
      do
        if(found) exit
        read(dtu,DSF) tmp1

        if(trim(tmp1) /= '') then
          s = scan(tmp1,' ')
          sgkeyword = tmp1(1:s-1)
          tmp2 = tmp1(s+1:DSL)
          s = scan(tmp2,' ')
          sglabel = tmp2(1:s-1)
        endif

        if(trim(sgkeyword) == 'END') then
          exit
        endif
        if(trim(sglabel) == trim(targetlabel)) then
          found = .true.
        else

        endif
      enddo

      if(.not. found) then
        write(*,'(A)') 'SG under consideration is '//trim(targetlabel)
        write(*,'(A)') 'err: in update_derivsg_table, cannot find SG label.'
        stop 1
      endif

      read(dtu,*) sg%T(1,1:3)
      read(dtu,*) sg%T(2,1:3)
      read(dtu,*) sg%T(3,1:3)
      write(*,*) 'To find keyword origvec:'
      read(dtu,DSF) line

      read(line,*) tmp1
      if(trim(tmp1) == '#optional') then
        read(line,*) tmp1,tmp2
        write(*,'(A)') 'tmp1,tmp2='//trim(tmp1),trim(tmp2)
        if(trim(tmp2) == 'origvec') then
          write(*,*) 'keyword origvec found.'
          read(dtu,*) sg%origvec(1:3)
        else
          write(*,*) 'I: keyword origvec not found'
          sg%origvec(1:3) = zero
        endif
      else
        write(*,*) 'II: keyword origvec not found'
        sg%origvec(1:3) = zero
      endif
      write(*,'(A,3F20.10)') 'sg%origvec(1:3)=', sg%origvec(1:3)
    endif

    sg%invT = inv3x3(sg%T(1,1))

    invT(1:3,1:3) = sg%invT(1:3,1:3)

    write(*,'(A)') 'Prepare to transform the translational and rotation parts: T='
    write(*,'(3F12.5)') sg%T(1,1:3)
    write(*,'(3F12.5)') sg%T(2,1:3)
    write(*,'(3F12.5)') sg%T(3,1:3)
    write(*,'(A,F12.5)') 'det3(sg%T) = ',det3(sg%T)
    write(*,'(A)') 'Prepare to transform the translational and rotation parts: invT=T^{-1}'
    write(*,'(3F12.5)') invT(1,1:3)
    write(*,'(3F12.5)') invT(2,1:3)
    write(*,'(3F12.5)') invT(3,1:3)

    nop = sg%nop
    write(*,*) 'nop = ',nop
    do i = 1, nop

      blk4x4(1:4,1:4) = sg%op(1:4,1:4,i)
      t1(1:3) = blk4x4(2:4,1)
      R1(1:3,1:3) = blk4x4(2:4,2:4)

      R2 = BAinvB(invT(1,1),R1(1,1))

      origvec(1:3) = sg%origvec(1:3)
      call matvecn(3,R1(1,1),origvec(1),t3(1))

      t1(1:3) = t1(1:3) + t3(1:3) - origvec(1:3)

      call matvecn(3,invT(1,1),t1(1),t2(1))

      call foldtofirstzone3(t2(1))

      newblk4x4(1:4,1:4) = zero
      newblk4x4(1,1) = one
      newblk4x4(2:4,1) = t2(1:3)
      newblk4x4(2:4,2:4) = R2(1:3,1:3)

      write(*,'(A,I5)') 'Pre || Post:  i =',i
      write(*,'(4F12.5,A,4F12.5)') blk4x4(2,1:4),' || ', newblk4x4(2,1:4)
      write(*,'(4F12.5,A,4F12.5)') blk4x4(3,1:4),' || ', newblk4x4(3,1:4)
      write(*,'(4F12.5,A,4F12.5)') blk4x4(4,1:4),' || ', newblk4x4(4,1:4)

      sg%op(1:4,1:4,i) = newblk4x4(1:4,1:4)
    enddo

    write(*,*)
    filename = 'echo-ext-operations.dat'
    write(*,'(A)') 'To produce '//trim(filename)

    pru = newunit()
    open(unit=pru,file=trim(filename),status='replace')

    write(pru,'(A)') '---------------------------------------'
    write(pru,'(A)') '    SG = '//trim(targetlabel)
    write(*,'(A)') '    SG = '//trim(targetlabel)
    write(pru,'(A)') '    Y = y'
    write(pru,'(A)') '    SGbase(Y,SG)%sglabel    ='''//trim(targetlabel)//''''
    write(pru,'(A)') '    SGbase(Y,SG)%PG         ='''//trim(sg%pg)//''''
    write(pru,'(A)') '    SGbase(Y,SG)%ncentering ='//trim(N2str(sg%ncentering))
    do i = 1, sg%ncentering
      do j = 1, 3
        f = sg%centering(j,i)

        call frac2nd(f,nf,df)
        if(abs(df-one) < 1.0d-6) then

          string(j)=trim(N2str(nint(nf)))//'.0d0'
        else
          string(j)=trim(N2str(nint(nf)))//'.0d0/'//trim(N2str(nint(df)))//'.0d0'
        endif
      enddo
      write(pru,'(A,A16,A,A16,A,A16,A)') '    SGbase(Y,SG)%centering(1:3,'//trim(N2str(i))//') = (/',trim(string(1)),',',trim(string(2)),',',trim(string(3)),' /)'
    enddo
    write(pru,'(A)') '    SGbase(Y,SG)%PGNumOper  ='//trim(N2str(sg%PGNumOper))
    write(pru,'(A)') '    SGbase(Y,SG)%symmorphic ='//trim(N2str(sg%symmorphic))
    if(sg%symmorphic == 0) then
      do i = 1, sg%PGNumOper
        do j = 1, 3
          f = sg%op(j+1,1,i)

          call frac2nd(f,nf,df)
          if(abs(df-one) < 1.0d-6) then

            string(j)=trim(N2str(nint(nf)))//'.0d0'
          else
            string(j)=trim(N2str(nint(nf)))//'.0d0/'//trim(N2str(nint(df)))//'.0d0'
          endif
        enddo
        write(pru,'(A,A16,A,A16,A,A16,A)') '    SGbase(Y,SG)%nsy(1:3,'//trim(N2str(i))//') = (/',trim(string(1)),',',trim(string(2)),',',trim(string(3)),' /)'
      enddo
    endif
    write(pru,'(A)') '    SGbase(Y,SG)%nop        ='//trim(N2str(sg%nop))
    write(pru,'(A)') '    SGbase(Y,SG)%assigned   =.TRUE.'
    write(pru,*)
    close(pru)
    write(*,'(A)') trim(filename)//' is generated.'
    write(*,*)

    write(*,*) 'original sg%sglabel = '//trim(sg%sglabel)
    sg%sglabel = trim(targetlabel)
    write(*,'(A)') 'Updated: sg%sglabel = '//trim(sg%sglabel)
  end subroutine update_derivsg_table

  subroutine assign_m01(m01,ns0,chPOSCAR,pPOSCAR,distance_thres)
    integer :: ns0,m01(ns0)
    type(supercell) :: chPOSCAR,pPOSCAR
    integer :: i,j,ii,jj,kk
    real(double) :: v1(3),v2(3),v3(3),newf(3),distance_thres
    m01(:) = -1
    do i = 1, ns0
      call frac2abs(chposcar%a,chposcar%at(i)%f,v1)
      do j = 1, pPOSCAR%n
        do kk = -1, 1
          do jj = -1, 1
            do ii = -1, 1
              newf = pPOSCAR%at(j)%f + (/ii,jj,kk/)
              call frac2abs(pPOSCAR%a,newf,v2)
              v3 = v2-v1
              if(vecmag3(v3) < distance_thres) then
                m01(i) = j
                exit
              endif
            enddo
          enddo
        enddo
      enddo
      if(m01(i) == -1) then
        write(*,*) 'atom', i, ' (in ch.poscar) cannot be found in pPOSCAR.'
        write(*,'(A)') 'err: serious mapping error in assign_m01'
        stop 1
      endif
    enddo
  end subroutine assign_m01

  subroutine assign_m12_tvec(tvec,m12,ns1,s1,s2,distance_thres)
    real(double) :: tvec(3)
    integer :: ns1,m12(ns1)
    type(supercell) :: s1,s2
    real(double) ::distance_thres
    integer :: i,j,ii,jj,kk
    real(double) :: v1(3),v2(3),v3(3),newf(3)
    real(double) :: mindis,dis
    m12(1:ns1) = -1

    do i = 1, ns1
      call frac2abs(s1%a(1,1),s1%at(i)%f(1),v1(1))

      v1(1:3) = v1(1:3) + tvec(1:3)

      mindis = 1d100
      do j = 1, s2%n

        do kk = -1, 1
          do jj = -1, 1
            do ii = -1, 1
              newf = s2%at(j)%f + (/ii,jj,kk/)
              call frac2abs(s2%a(1,1),newf(1),v2(1))
              v3 = v2-v1
              dis = vecmag3(v3(1))
              if(dis < mindis) then
                mindis = dis
              endif
              if(dis < distance_thres) then
                m12(i) = j
                exit
              endif
            enddo
          enddo
        enddo
      enddo

      if(m12(i) == -1) then
        write(*,*) 'mindis = ',mindis
        write(*,*) 'atm index=', i, ', cannot be found in supercell.poscar.'
        write(*,'(A)') 'err: serious mapping error in assign_m12'
        stop 1
      endif
    enddo

  end subroutine assign_m12_tvec

  subroutine create_primitive_cell(sp,s,t,kstar_arr)
    type(onesg) :: sp
    type(supercell) :: s
    type(supercell) :: t
    integer :: i,j,k
    real(double) :: vec(4,1),prod(4,1),p(3),absp(3)
    logical :: newatom
    integer :: indexu
    type(onekstar) :: kstar_arr(s%n)
    integer :: runindex,opind
    real(double) :: v1(4),v2(4),v3(3),matop(4,4),p1(3),p2(3),newf(3),rotm(3,3)

    t%n = 0
    t%a = s%a
    call real2recip(t%a,t%b)

    indexu = newunit()

    open(unit=indexu,file='ineq_atms_loc.dat',status='replace')

    t%nsp = s%nsp
    t%zarr(1:t%nsp) = s%zarr(1:s%nsp)

    do i = 1, s%n
      k = 0
      do j = 1, sp%nop

        vec(1,1) = 1d0
        vec(2:4,1) = s%at(i)%f

        call matmatmnp(4,4,1,sp%op(1:4,1:4,j),vec(1:4,1),prod(1:4,1))

        p = prod(2:4,1)
        call frac2abs(s%a,p,absp)
        call tryput(s%at(i)%z,absp,t,newatom)
        if(newatom) then
          k = k + 1
          if(k == 1) then

            write(indexu,'(I8)') t%n
          else
            write(indexu,'(A8)') 'equiv'
          endif
          kstar_arr(i)%op(k) = j
          kstar_arr(i)%ind1(k) = t%n
          kstar_arr(i)%trans(1:3,k) = (/0,0,0/)
        endif
      enddo
      kstar_arr(i)%n = k
    enddo
    close(unit=indexu)
    write(*,*)

    write(*,'(A,I5)') 'Full poscar has totaln = ',t%n

    write(*,*)

    runindex = 0
    do i = 1, s%n
      v1(1:4) = (/one,s%at(i)%f(1:3)/)
      do j = 1, kstar_arr(i)%n
        runindex = runindex + 1
        opind = kstar_arr(i)%op(j)
        matop(1:4,1:4) = sp%op(1:4,1:4,opind)
        call matvecn(4,matop(1,1),v1(1),v2(1))

        v3(1:3) = v2(2:4)
        call foldtofirstzone3(v3(1))

        call frac2abs(s%a,v3(1),p1(1))

        call foldtofirstzone3(t%at(runindex)%f)

        call frac2abs(t%a,t%at(runindex)%f,p2(1))

        if(vecmag3(p1-p2) > 1.0d-5) then
          write(*,'(A,I5,6F20.10,F20.10)') 'runindex,p1,p2=',runindex,p1(1:3),p2(1:3),vecmag3(p1-p2)
          write(*,'(A)') 'err: star of k is wrong?'
          stop 1
        endif

        rotm(1:3,1:3) = matop(2:4,2:4)
        call matvecn(3,rotm(1,1),s%at(i)%force(1),newf(1))
        t%at(runindex)%force(1:3) = newf(1:3)
      enddo
    enddo
  end subroutine create_primitive_cell

  subroutine find_group_of_k(ns0,chposcar,pposcar,sg,m01,kgroup_arr)
    integer :: ns0
    type(supercell) :: chposcar,pposcar
    type(onesg) :: sg
    integer :: m01(ns0)
    type(onekgroup) :: kgroup_arr(ns0)
    integer :: i,j,gs,ii,jj,kk
    integer :: ind1
    real(double) :: fracA(3),vec(4,1),prod(4,1),p(3),q(3),diffv(3),dis
    logical :: found

    integer :: absnrange=5

    ns0 = chposcar%n
    do i = 1, ns0
      ind1 = m01(i)
      gs = 0
      fracA = pposcar%at(ind1)%f
      do j = 1, sg%nop
        vec(1,1) = one
        vec(2:4,1) = fracA
        call matmatmnp(4,4,1,sg%op(1:4,1:4,j),vec(1:4,1),prod(1:4,1))
        p = prod(2:4,1)
        found = .false.

        do ii = -absnrange, absnrange
          do jj = -absnrange, absnrange
            do kk = -absnrange, absnrange

              q = p + (/ii,jj,kk/)
              diffv = q-fracA
              dis = vecmag3(diffv)
              if(dis < 1d-8) then
                found = .true.
                gs = gs + 1
                kgroup_arr(i)%op(gs) = j
              endif
            enddo
          enddo
        enddo
      enddo
      kgroup_arr(i)%n = gs
    enddo
  end subroutine find_group_of_k

  subroutine get_supercell_atom_map(u,pposcar,superc,map)

    real(double) :: u(4,4)
    type(supercell) :: pposcar,superc,tmp1
    integer :: map(superc%n)
    integer :: j,k
    real(double) :: jabsc(3),jfracc(3), prod(4,4),v(4,1)
    real(double) :: dev,testt(3),diffv(3),rundev
    integer :: n1scan,n2scan,n3scan
    real(double) :: summap,summap_analytic
    integer :: ns2
    integer,allocatable :: invmap(:)

    tmp1%a = superc%a
    call real2recip(tmp1%a,tmp1%b)
    tmp1%n = superc%n
    do j = 1, tmp1%n
      tmp1%at(j)%z = superc%at(j)%z
      tmp1%at(j)%f = superc%at(j)%f
    enddo

    do j = 1, tmp1%n

      call frac2abs(tmp1%a, tmp1%at(j)%f(1:3),jabsc)
      call abs2frac(pposcar%b,jabsc,jfracc)

      v(1,1) = one
      v(2:4,1) = jfracc(1:3)

      call matmatmnp(4,4,1,u(1:4,1:4),v(1:4,1),prod(1:4,1))

      jfracc(1:3) = prod(2:4,1)
      call frac2abs(pposcar%a,jfracc,jabsc)
      call abs2frac(tmp1%b,jabsc,tmp1%at(j)%f)

      call foldtofirstzone3(tmp1%at(j)%f)
    enddo

    allocate(invmap(superc%n))

    map(:) = 0
    invmap(:) = 0

    do j = 1, tmp1%n

      do k = 1, superc%n

        dev = 1d100
        do n1scan = -1, 1
          do n2scan = -1, 1
            do n3scan = -1, 1
              testt = tmp1%at(j)%f + (/n1scan,n2scan,n3scan/)
              diffv = testt - superc%at(k)%f
              rundev = vecmag3(diffv)
              if(rundev < dev) then
                dev = rundev
              endif
            enddo
          enddo
        enddo

        if(dev < dev_tol_get_supercell_atom_map) then
          if(map(j) /= 0) then
            write(*,*) 'Two atoms in the old list match with this atom j'
            write(*,'(A)') 'err: Check again.'
            stop 1
          endif
          map(j) = k
          if(invmap(k) /= 0) then
            write(*,*) 'j,k = ',j,k
            write(*,*) 'k is pre-mapped to a previous atom=',invmap(k)
            write(*,'(A)') 'err: Check again.'
            stop 1
          endif
          invmap(k) = j
        endif
      enddo
      if(map(j) == 0) then
        write(*,*) 'problematic j = ',j
        write(*,*) 'something is wrong. mapj must be nonzero'
        write(*,'(A)') 'err: stop here.'
        stop 1
      else

      endif
    enddo

    deallocate(invmap)

    ns2 = superc%n
    summap = 0.0d0
    do j = 1, ns2
      if(map(j) < 0 .or. map(j) > ns2) then
        write(*,*) 'j,map(j)=',j,map(j)
        write(*,'(A)') 'err: sanity check fails.'
        stop 1
      endif
      summap = summap + map(j)
    enddo
    summap_analytic=ns2*(one+ns2)/two
    if(abs(summap - summap_analytic) > 1.0d-10) then
      write(*,*) 'summap,summap_analytic=',summap,summap_analytic
      write(*,*) 'map array sanity check fails.'
      write(*,'(A)') 'err: check.'
      stop 1
    endif
  end subroutine get_supercell_atom_map

  subroutine write_forces(ofu,absfinalforcesfile,ns1,ns2,zf,nzf)
    integer :: ofu
    character(len=*) :: absfinalforcesfile
    integer :: ns1,ns2
    real(double) :: zf(3,ns2)
    real(double) :: nzf(3,ns2,3,ns1,2)
    integer :: i,j,k,m

    write(*,'(A)') 'write: To open file ='//trim(absfinalforcesfile)
    open(unit=ofu,file=trim(absfinalforcesfile),status='replace')
    do i = 1, ns2
      write(ofu,*) zf(1:3,i)

    enddo
    do i = 1, ns1
      do j = 1, 3
        do m = 1, 2
          do k = 1, ns2

            write(ofu,*) nzf(1:3,k,j,i,m)

          enddo
        enddo
      enddo
    enddo
    close(ofu)
    write(*,'(A)') 'file='//trim(absfinalforcesfile)//' is now closed'
  end subroutine write_forces

  subroutine report_operation_type(dm_int,tr_int)
    integer :: dm_int,tr_int

    if(dm_int == 1 .and. tr_int == -1) then
      write(*,*) 'operation is 2 (2-fold rotation)'
    else if(dm_int == 1 .and. tr_int == 0) then
      write(*,*) 'operation is 3 (3-fold rotation)'
    else if(dm_int == 1 .and. tr_int == 1) then
      write(*,*) 'operation is 4 (4-fold rotation)'
    else if(dm_int == 1 .and. tr_int == 2) then
      write(*,*) 'operation is 6 (6-fold rotation)'
    else if(dm_int == 1 .and. tr_int == 3) then
      write(*,*) 'operation is 1 (identity)'

    else if(dm_int == -1 .and. tr_int == 1) then
      write(*,*) 'operation is -2 (mirror)'
    else if(dm_int == -1 .and. tr_int == 0) then
      write(*,*) 'operation is -3'
    else if(dm_int == -1 .and. tr_int == -1) then
      write(*,*) 'operation is -4'
    else if(dm_int == -1 .and. tr_int == -2) then
      write(*,*) 'operation is -6'
    else if(dm_int == -1 .and. tr_int == -3) then
      write(*,*) 'operation is -1 (inversion)'
    else
      write(*,*) 'determinant,trace= ',dm_int,tr_int
      write(*,'(A)') 'err: bad combination of trace and dmerminant.'
      stop 1
    endif
  end subroutine report_operation_type

  subroutine frac2nd(f,n,d)
    real(double) :: f,n,d,absf,intf,ffrac,p,intp,pfrac
    integer :: i
    integer :: nscan=1000
    real(double),parameter :: eps=1.0d-6
    logical :: found

    if(abs(f) < eps) then
      n = zero
      d = one
      return
    endif

    if(f > zero) then
      absf = f
    else
      absf = -f
    endif

    intf = nint(absf)
    ffrac = absf - intf

    if(abs(ffrac) < eps) then

      n = f
      d = one
      return
    endif

    found = .false.
    do i = 1, nscan
      p = i*ffrac

      intp = nint(p)
      pfrac = p - intp
      if(abs(pfrac) < eps) then
        n = intf*i + intp
        d = i
        if(f < zero) then
          n = -n
        endif

        if(abs(f - n/d) > eps) then
          write(*,*) 'nscan = ',nscan
          write(*,*) 'f,n,d,n/d,f-n/d=',f,n,d,n/d,f-n/d
          write(*,'(A)') 'err: wrong.'
          stop 1
        endif
        found = .true.
        return
      endif
    enddo
    if(.not. found) then
      write(*,*) 'we cannot find n and d for f=',f
      write(*,*) 'eps =',eps
      write(*,*) 'nscan = ',nscan
      write(*,'(A)') 'err: check the value of f or decrease eps or increase nscan?'
      stop 1
    endif
  end subroutine frac2nd
end module commod
