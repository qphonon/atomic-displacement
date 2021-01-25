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

module types_and_constants
  implicit none

  integer,parameter::double=kind(1.0D0)
  integer,parameter::DSL=500
  character(len=*),parameter :: DSF='(A500)'
  real(double),parameter :: one=1.0D0,minusone=-one,zero=0.0D0,two=2.0D0,&
                            three=3.0D0,four=4.0D0,five=5.0D0,&
                            onethird=one/three,third=onethird,&
                            twothird=two/three,half=one/two,&
                            quarter=one/four,fifth=one/five,sixth=one/6.0D0
  real(double),parameter :: sqrt2=1.41421356237310D0
  real(double),parameter :: sqrt3=1.73205080756888D0
  real(double),parameter :: sqrt5=2.23606797749979D0

  real(double),parameter :: pi = 3.14159265358979323846D0
  real(double),parameter :: twopi = two*pi
  real(double),parameter :: efac =  2.718281828459045235D0
  real(double),parameter :: eulerg =0.5772156649015328606D0

  real(double),parameter :: m_pi=pi

  complex(double),parameter :: Imag=(zero,one)
  complex(double),parameter :: CZero=(zero,zero)
  complex(double),parameter :: COne=(one,zero)

  real(double),parameter :: CLight = 299792458.0D0
  real(double),parameter :: ECHG = 1.60217653D-19
  real(double),parameter :: hbar = 1.05457168D-34
  real(double),parameter :: EMss = 9.1093826D-31
  real(double),parameter :: kBoltz = 1.3806505D-23
  real(double),parameter :: NAvo = 6.0221415D23
  real(double),parameter :: Mu0 = 4.0D-7*pi
  real(double),parameter :: GGrav=6.6742D-11

  CHARACTER(LEN=1),DIMENSION(0:9),PARAMETER :: DIGITARR=(/'0','1','2','3','4','5','6','7','8','9'/)
  integer,parameter :: zordmax=104

  CHARACTER(LEN=2),DIMENSION(zordmax),PARAMETER:: ATS=(/'h ','he','li','be','b ','c ', &
     'n ','o ','f ','ne','na','mg','al','si','p ','s ','cl','ar','k ', &
     'ca','sc','ti','v ','cr','mn','fe','co','ni','cu','zn','ga','ge', &
     'as','se','br','kr','rb','sr','y ','zr','nb','mo','tc','ru','rh', &
     'pd','ag','cd','in','sn','sb','te','i ','xe','cs','ba','la','ce', &
     'pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu', &
     'hf','ta','w ','re','os','ir','pt','au','hg','tl','pb','bi','po', &
     'at','rn','fr','ra','ac','th','pa','u ','np','pu','am','cm','bk', &
     'cf','es','fm','md','no','lr','aa'/)

  CHARACTER(LEN=2),DIMENSION(zordmax),PARAMETER:: CATS=(/'H ','He','Li','Be','B ','C ', &
     'N ','O ','F ','Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar','K ', &
     'Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge', &
     'As','Se','Br','Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh', &
     'Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce', &
     'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu', &
     'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po', &
     'At','Rn','Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk', &
     'Cf','Es','Fm','Md','No','Lr','Aa'/)

  integer,parameter :: pg_uuid_L5=100000
  integer,parameter :: pg_uuid_L4=10000
  integer,parameter :: pg_uuid_L3=1000
  integer,parameter :: pg_uuid_L2=100

  integer,parameter :: PG_UUID(1:32) =(/ 113303,206006,224202,204402,408004,446000,406400,812000,&
                                         446404,406004,812008,557101,517501,517101,1014002,333303,606006,&
                                         334202,314402,608004,668606,608006,1216012,668202,628602,608202,&
                                         1216004,444202,808004,556202,516202,1012004 /)

  CHARACTER(DSL) :: PG_NAME(32)

  REAL(DOUBLE),DIMENSION(zordmax),PARAMETER :: MASSOFATOM=(/1.0079D0,4.0026D0,6.941D0,9.0122D0,10.81D0,12.01D0,&
     14.007D0,15.999D0,18.998D0,20.18D0,22.9898D0,24.305D0,26.982D0,28.086D0,30.974D0,32.064D0,&
     35.453D0,39.948D0,39.09D0,40.08D0,44.956D0,47.90D0,50.942D0,52.00D0,54.938D0,55.85D0,58.93D0,&
     58.71D0,63.55D0,65.38D0,69.72D0,72.59D0,74.922D0,78.96D0,79.91D0,83.80D0,85.47D0,87.62D0,&
     88.91D0,91.22D0,92.91D0,95.94D0,98.91D0,101.07D0,102.90D0,106.40D0,107.87D0,112.40D0,114.82D0,118.69D0,&
     121.75D0, 127.60D0, 126.90D0, 131.30D0, 132.91D0,137.34D0,138.91D0,140.12D0,140.91D0,144.24D0,145.0D0,&
     150.35D0,151.96D0,157.25D0,158.92D0,162.50D0,164.93D0,167.26D0,168.93D0,173.04D0,174.97D0,178.49D0,&
     180.95D0,183.85D0,186.2D0,190.2D0,192.22D0,195.09D0,196.97D0,200.59D0,204.37D0,207.19D0,208.98D0,&
     210.0D0,210.0D0,222.0D0,223.0D0,226.0D0,227.0D0,232.04D0,231.0D0,238.03D0,237.05D0,244.0D0,243.0D0,247.0D0,&
     247.0D0,251.0D0,254.0D0,257.0D0,256.0D0,254.0D0,257.0D0,55.85d0/)

  REAL(DOUBLE),DIMENSION(zordmax),PARAMETER :: COHBOFATOM=(/1.0079D99,4.0026D99,6.941D99,9.0122D99,10.81D99,12.01D99,&
     14.007D99,15.999D99,18.998D99,20.18D99,22.9898D99,24.305D99,26.982D99,28.086D99,30.974D99,32.064D99,&
     35.453D99,39.948D99,39.09D99,40.08D99,44.956D99,47.90D99,50.942D99,52.00D99,54.938D99,55.85D99,58.93D99,&
     58.71D99,63.55D99,65.38D99,69.72D99,72.59D99,74.922D99,78.96D99,79.91D99,83.80D99,85.47D99,87.62D99,&
     88.91D99,91.22D99,92.91D99,95.94D99,98.91D99,101.07D99,102.90D99,106.40D99,107.87D99,112.40D99,114.82D99,118.69D99,&
     5.57D0,5.80D0, 126.90D99, 131.30D99, 132.91D99,137.34D99,138.91D99,140.12D99,140.91D99,144.24D99,145.0D99,&
     150.35D99,151.96D99,157.25D99,158.92D99,162.50D99,164.93D99,167.26D99,168.93D99,173.04D99,174.97D99,178.49D99,&
     180.95D99,183.85D99,186.2D99,190.2D99,192.22D99,195.09D99,196.97D99,200.59D99,204.37D99,207.19D99,208.98D99,&
     210.0D99,210.0D99,222.0D99,223.0D99,226.0D99,227.0D99,232.04D99,231.0D99,238.03D99,237.05D99,244.0D99,243.0D99,247.0D99,&
     247.0D99,251.0D99,254.0D99,257.0D99,256.0D99,254.0D99,257.0D99,1.0d0/)

  CHARACTER(LEN=*),PARAMETER :: UPPER='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  CHARACTER(LEN=*),PARAMETER :: LOWER='abcdefghijklmnopqrstuvwxyz'
  CHARACTER(LEN=1),PARAMETER :: RTRN=CHAR(10)

  character(len=8) :: DATE
  character(len=10) :: TIME

  integer,parameter :: shellarr(21) = (/0,-1,1,-2,2,-3,3,-4,4,-5,5,&
                                       -6,6,-7,7,-8,8,-9,9,-10,10/)

  real(double),parameter :: deg2rad = pi/180.0D0
  real(double),parameter :: rad2deg = 180.0D0/pi

  real(double),parameter :: goldenratio=(sqrt5+one)/two
  real(double),parameter :: goldengamma=(sqrt5-one)/two

  type internalunitinsi
    real(double) :: M
    real(double) :: T
    real(double) :: L
    real(double) :: I
    real(double) :: Q
  end type internalunitinSI

  type internalfundamentalconstant
    real(double) :: internalmu0
    real(double) :: internalECharge
    real(double) :: internalmuB
    real(double) :: internalhbar
    real(double) :: internalkB
  end type internalfundamentalconstant

  type oneatom
    integer :: z
    real(double) :: p(3)
    real(double) :: force(3)
  end type oneatom

  type anatom
     integer :: z
     real(double) :: chg
     real(double) :: mass
     real(double) :: f(3),ac(3)
     real(double) :: force(3)
     integer :: gr
  end type anatom

  integer,parameter :: maxnum = 10000
  integer,parameter :: maxcs = 1000

  type supercell
     character(LEN=DSL) :: commentline
     character(LEN=DSL) :: fracabs
     real(double) :: la(6)
     real(double) :: a(3,3)
     real(double) :: b(3,3)
     real(double) :: vol
     integer :: n
     integer :: nsp
     integer :: n_per_species(zordmax)
     integer :: zarr(zordmax)
     type(anatom) :: at(maxnum)
     character(LEN=DSL) :: sg_label
     real(double) :: density
  end type supercell

  integer,parameter :: nopmax=300
  integer,parameter :: nclassesmax=100
  type onesg
    logical :: assigned
    character(len=DSL) :: sglabel
    character(len=DSL) :: pg
    integer :: nclass
    integer :: classsize(nclassesmax)
    integer :: classindex(nclassesmax)
    real(double) :: charac(nclassesmax)
    integer :: nop
    integer :: PGNumOper
    integer :: ncentering
    integer :: reorder_IT_op_sequence
    real(double) :: centering(3,4)
    real(double) :: op(4,4,nopmax)
    real(double) :: invT(3,3)
    real(double) :: T(3,3)
    real(double) :: origvec(3)
    integer :: symmorphic
    real(double) :: nsy(3,48)
    integer :: identity_map(48)
    integer :: operation_map(48)
  end type onesg

  integer,parameter :: maxSGvariant=10
  type(onesg) :: SGbase(maxSGvariant,0:230)

  type onepointgroup
    character(len=DSL) :: pgname
    integer :: nop
    integer :: nc
    integer :: classsize(nclassesmax)
    real(double) :: charac_mat(nclassesmax,nclassesmax)
    character(len=100) :: irrep(nclassesmax)
    real(double) :: m(3,3,100)
  end type onepointgroup

  type onekstar
    integer :: n
    integer :: op(maxcs)
    integer :: ind1(maxcs)
    integer :: trans(3,maxcs)
  end type onekstar

  type rmatrix33
    real(double) :: m(3,3)
  end type rmatrix33

  type rvec4
    real(double) :: v(4)
  end type rvec4

  integer,parameter :: maxops=1000
  type onekgroup
    character(len=DSL) :: kname
    integer :: n
    integer :: op(maxops)
    integer :: trans(3,maxops)
    real(double) :: det(maxops)
    real(double) :: trace(maxops)
    logical :: inversion
  end type onekgroup

  type memory_type
    integer :: n_int
    integer :: n_double
    integer :: n_dcomplex
    real(double) :: totmem
  end type memory_type

  type time_type
    real(double) :: init
    real(double) :: final
    real(double) :: total
  end type time_type

  integer,parameter :: WDMAX=500
  type RNGspec
    integer :: WDINTARR(WDMAX),WDLEN,RNGMODE
    integer :: SEEDSTR(WDMAX),SEEDLEN
  end type RNGspec
  type(RNGspec) :: globalRNG

  character(len=WDMAX) :: DATA3
  character(len=WDMAX) :: MYLIB
  character(len=WDMAX) :: MYPROG

  type rhombopara
    real(double) :: ah,ch
    real(double) :: ar,alphardeg
    real(double) :: alpha
    real(double) :: cosalpha
    real(double) :: caratio
    real(double) :: mat0(3,3)
    real(double) :: mat1(3,3)
    real(double) :: mat2(3,3)
    real(double) :: mat3(3,3)
    real(double) :: mat4(3,3)
    real(double) :: mat5(3,3)
    real(double) :: mat6(3,3)
  end type rhombopara

  integer,parameter :: DTYPEMAXLEN=20
  type dtype
    logical :: li
    logical :: lr
    logical :: ls
    integer :: i
    real(double) :: r
    character(len=DSL) :: s
    integer :: ivec(DTYPEMAXLEN)
    real(double) :: rvec(DTYPEMAXLEN)
    character(len=DSL) :: svec(DTYPEMAXLEN)
  end type dtype

  integer,parameter :: DEFAULTOPT=10
  integer,parameter :: NONDEFAULTOPT=20
  integer,parameter :: itype=1
  integer,parameter :: rtype=2
  integer,parameter :: stype=3

  type strarrtype
    integer :: n
    character(len=DSL) :: s(100)
  end type strarrtype

contains

  subroutine assign_env_variable(ename,strv)
    character(len=*) :: ename,strv
    call get_environment_variable(ename,strv)
    write(*,*) 'system variable '//trim(ename)//' is '//trim(strv)
  end subroutine assign_env_variable

end module types_and_constants

