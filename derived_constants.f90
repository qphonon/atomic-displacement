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

module derived_constants
  use types_and_constants
  implicit none

  real(double),parameter :: centimeter=1.0D-2
  real(double),parameter :: nm = 1.0D-9
  real(double),parameter :: Angstrom = 1.0D-10
  real(double),parameter :: RGas = kBoltz*NAvo

  real(double),parameter :: hPlanck=hbar*twopi
  real(double),parameter :: eV = EChg
  real(double),parameter :: Eps0 = one/(CLight*CLight*Mu0)
  real(double),parameter :: epsilon0 = eps0
  real(double),parameter :: Bohr = four*Pi*Eps0*hbar*hbar/(EMss*EChg*EChg)
  real(double),parameter :: FineConst = hbar/(EMss*CLight*Bohr)

  real(double),parameter :: tau0bohr=twopi*EMss*Bohr**2/hbar

  real(double),parameter :: Hartree = hbar*hbar/(EMss*Bohr*Bohr)
  real(double),parameter :: Hartree2 = EChg*EChg/(four*Pi*Eps0*Bohr)
  real(double),parameter :: Hartree3 = (ECHG**two/(four*Pi*Eps0))**2*EMss/hbar**two

  real(double),parameter :: Rydberg=half*Hartree
  real(double),parameter :: Stefan_Boltz =pi**two*kBoltz**four/hbar**three/60.0D0/CLight**two
  real(double),parameter :: Stefan_Boltz2=two*pi**five*kBoltz**four/hPlanck**three/15.0D0/CLight**two

  real(double),parameter :: BohrMagneton=EChg*hbar/(two*EMss)

  real(double),parameter :: QuantumConductanceG0=EChg*EChg/hPlanck
  real(double),parameter :: QuantizedHallResistance= one/QuantumConductanceG0

  real(double),parameter :: OneAUTime = EMss*Bohr*Bohr/hbar
  real(double),parameter :: OneAUTime1 =  (four*Pi*Eps0)**2*hbar**3/EMss/eV**four

  real(double),parameter :: unittime_HartreeBohr2_SI_1 = sqrt(EMss*Bohr**2/Hartree)
  real(double),parameter :: unittime_HartreeBohr2_SI_2 = (EMss*Bohr**2/Hartree)**half

  real(double),parameter :: Bohr2Ang = Bohr/Angstrom
  real(double),parameter :: BohrToAng = Bohr2Ang
  real(double),parameter :: Ang2Bohr = Angstrom/Bohr
  real(double),parameter :: AngToBohr = Ang2Bohr
  real(double),parameter :: Hartree2eV = Hartree/eV
  real(double),parameter :: Ryd2eV = Rydberg/eV
  real(double),parameter :: eV2Ryd = one/Ryd2eV
  real(double),parameter :: HartreeToeV = Hartree2eV
  real(double),parameter :: eV2Hartree = one/HartreeToeV
  real(double),parameter :: eVToHartree = eV2Hartree
  real(double),parameter :: Ang2nm = 0.1D0
  real(double),parameter :: nm2Ang = 10.0D0
  real(double),parameter :: J2eV=one/eV
  real(double),parameter :: DebyeRinBohr = 1.0D-18/(CLight*10.0D0)*1.0D-2/EChg/Bohr
  real(double),parameter :: Debye = DebyeRinBohr*Bohr*EChg
  real(double),parameter :: m2Ang = one/Angstrom
  real(double),parameter :: eVperAng = eV/Angstrom
  real(double),parameter :: eA2N = eVperAng
  real(double),parameter :: N2eA = one/eA2N

  real(double),parameter :: kBoltzDeV = kBoltz/eV
  real(double),parameter :: Hartree2K = Hartree/kBoltz

  real(double),parameter :: kg2g = 1.0D3
  real(double),parameter :: g2kg = one/kg2g
  real(double),parameter :: m2cm = 1.0D2
  real(double),parameter :: cm2m = one/m2cm

  real(double),parameter :: erg = g2kg*cm2m**2
  real(double),parameter :: J2erg = one/erg
  real(double),parameter :: dyn = g2kg*cm2m
  real(double),parameter :: CLight_cgs = CLight*m2cm
  real(double),parameter :: hbar_cgs = hbar*kg2g*m2cm**2
  real(double),parameter :: ECHG_CGS = ECHG*CLight*10.0D0
  real(double),parameter :: EMSS_cgs = EMSS*kg2g
  real(double),parameter :: bohr_cgs = hbar_cgs**2/(EMss_cgs*echg_cgs**2)
  real(double),parameter :: hartree_cgs = ECHG_CGS**2/bohr_cgs
  real(double),parameter :: Tesla2Gauss = 1.0D4
  real(double),parameter :: BohrMagneton_cgs = EChg_cgs*hbar_cgs/(two*emss_cgs*clight_cgs)

  real(double),parameter :: Oe2SI = 1000.0D0/(four*pi)
  real(double),parameter :: emucc2Adivm = 1.0D3
  real(double),parameter :: ergdivcm2Jdivm = erg/(cm2m)

  real(double),parameter :: Cal2J = 4.18400D0
  real(double),parameter :: J2Cal = one/Cal2J
  real(double),parameter :: pound2kg = 0.45359237D0
  real(double),parameter :: in2m=0.0254D0
  real(double),parameter :: m2in=one/in2m
  real(double),parameter :: foot=12.0D0*in2m
  real(double),parameter :: squarefeet=foot*foot
  real(double),parameter :: sm2sf=one/squarefeet

  real(double),parameter :: Hartree2kCalMol = Hartree*J2Cal*NAvo/1.0D3
  real(double),parameter :: eV2kCalMol = eV2Hartree*Hartree2kCalMol
  real(double),parameter :: kCalMol2eV = one/eV2kCalMol

  real(double),parameter :: Kilo=1.0D3,Mega=1.0D6,Giga=1.0D9,Tera=1.0D12
  real(double),parameter :: Mili=1.0D-3,Micro=1.0D-6,Nano=1.0D-9,Pico=1.0D-12,Femto=1.0D-15

  real(double),parameter :: AMU = 1.6605402D-27

  real(double),parameter :: AMU2 = 1.0D-3/NAvo

  real(double),parameter :: PM_U = 1.00727646688D0
  real(double),parameter :: PMss = PM_U*AMU2

  real(double),parameter :: THz=1.0D12
  real(double),parameter :: Hz2THz=1.0D-12
  real(double),parameter :: THz2Hz=one/Hz2THz
  real(double),parameter :: invcm=CLight/centimeter
  real(double),parameter :: Hz2invcm=one/invcm
  real(double),parameter :: invcm2THz=invcm*Hz2THz
  real(double),parameter :: THz2invcm=one/invcm2THz
  real(double),parameter :: invcm2eV=hPlanck*invcm/eV
  real(double),parameter :: Rydoverh=Rydberg/hPlanck
  real(double),parameter :: Rydoverh2invcm=Rydoverh/invcm
end module derived_constants
