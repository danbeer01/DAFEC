module m_Constants_Mod
!
! PURPOSE:
!   To define constants and functions that are used throughout the program
!
! Record of Revsions:
!   Date        Programer       Change
!   04/11/22    J.I.M           Original code
!

  Implicit None

  !! Private definition

  !! Precision Definitions
  Integer, Parameter        :: sp=selected_Real_Kind(4)   !! single precision
  Integer, Parameter        :: dp=selected_Real_Kind(8)  !! double precision
  Integer, Parameter        :: qp=selected_Real_Kind(16)  !! quad precision

  !! Constant Values
  Real(Kind=dp),  Parameter :: PI=4.D0*DATAN(1.0D0)
  Real(Kind=dp),  Parameter :: dp_EPSILON = 1.0E-12_dp
  Real(Kind=dp),  Parameter :: VSMALL_NUMBER = 1.0E-9_dp
  Real(Kind=dp),  Parameter :: SMALL_NUMBER = 1.0E-6_dp
  Real(Kind=dp),  Parameter :: LARGE_NUMBER = 1.0E+6_dp
  Real(Kind=dp),  Parameter :: VLARGE_NUMBER = 1.0E+9_dp
  INTEGER,PARAMETER         :: MAX_ITERATIONS = 10000
  real(dp),parameter        :: ADJUSTED_NUETRON_MASS = 1.04625e-8_dp
  
  !! Characters
  Character, Parameter          :: COMMENT_CHAR = '!'
  Character(len=2), Parameter   :: tab = "  "
  Character, Parameter          :: space = " "


  contains

  REAL(kind =dp) function interpolate_func(x1,x2,x3,val1,val2)
      real(KIND = dp) :: x1,x2,x3,val1,val2
      interpolate_func = ((val2-val1)*(x3-x1))/(x2-x1)+val1
  end function

  REAL(kind =dp) function twoD_interpolate_func(x1,x2,x3,y1,y2,y3,val11,val12,val21,val22)
      real(KIND = dp) :: x1,x2,x3,y1,y2,y3,val11,val12,val21,val22,av1,av2,av3,av4
      av1 = interpolate_func(x1,x2,x3,val11,val12)
      av2 = interpolate_func(x1,x2,x3,val21,val22)
      av3 = interpolate_func(y1,y2,y3,val11,val21)
      av4 = interpolate_func(y1,y2,y3,val12,val22)
      twoD_interpolate_func = (av1+av2+av3+av4)/4
  end function


end module
