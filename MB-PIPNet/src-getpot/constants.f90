module constants
  implicit none

  ! Define the mass of different atoms
  real,parameter::c_mass= 12.0000000  !21874.66
  real,parameter::d_mass=  2.0135532127
  real,parameter::h_mass=  1.0078250  !1837.15
  real,parameter::o_mass= 15.9949146  !29156.95
  real,parameter::ar_mass=39.9623831  !72846.97
  real,parameter::cl_mass=34.968853388
  !real,parameter::pi=acos(-1.0)

  ! Define constants
  real,parameter::raddeg=57.2957795
  real,parameter::auang=0.5291772083
  real,parameter::aucm=219474.6313710
  real,parameter::cmau=0.00000456
  real,parameter::aukcal=627.51
  real,parameter::cmhz=29979250000.0
  real,parameter::auhz=41341370000000000.0
  real,parameter::audb=0.393430307

  real,parameter::emass=1822.88848

   double precision, parameter ::     PI= 3.14159265358979323844d0
   double precision, parameter ::  TWOPI= 6.28318530717958647688d0
   double precision, parameter :: FOURPI=12.56637061435917295376d0
   double precision, parameter ::    TSP= 1.12837916709551257390d0 !  2/sqrt(pi)
   double precision, parameter ::   SQPI= 1.77245385090551602729d0 !  sqrt(pi)
   double precision, parameter :: DEG2RAD= .01745329251994329576d0   !degrees--> Radian
   !...
   double precision, parameter :: HBAR =0.0151704885d0 ! kcal/mol * ps
   double precision, parameter :: KCALKLVN=503.21659237841666799453d0
   double precision, parameter :: PRESCON=6.85695d+4  ! kcal/mole/Ang**3 to Atm
   double precision, parameter :: MASSCON=418.4d0!at.weig->(kcal/mol)*(ps/Ang)^2
   double precision, parameter :: DENSCON = 1.66053872801494672163d0*MASSCON
!   double precision, parameter :: CHARGECON = 18.22261720426243437986d0
   double precision, parameter :: CHARGECON = 18.22261534447435150793d0 !tinker
   double precision, parameter :: DEBYE  = 4.8033324d0/CHARGECON
   double precision, parameter :: FS2PS=0.001d0
   double precision, parameter :: CMIPSI=.18836515673088532773d0! cm-1 --> ps-1


end module constants
