C***********************************************************************
      subroutine chanco(V,UB,BV,VR,RI,g,xug,xgr,xri,xiz,xgi)
C---------------------------------------------------------------------
C     Transformation between the Johnson-Cousins UBVRI photometry 
C     system and the SDSS ugriz system.
C     Follow the eq. 1 to 8 from Jordi,Grebel & Ammon, 2006, A&A, 460
C---------------------------------------------------------------------
      implicit real (a-h,m,o-z)
      
      real a1,a2,a3,a5,a7,b1,b2,b3,b5,b7,c5
      real V,UB,BV,VR,RI,g,xug,xgr,xri,xiz,xgi

C     --- Used parameters (Table 3) ----
      parameter(a1=0.630, b1=-0.124)
      parameter(a2=1.007, b2=-0.236)
      parameter(a3=1.584, b3=-0.386)
      parameter(a5=0.750, b5=0.770, c5=0.720)
      parameter(a7=1.646, b7=-0.139)      

C     --- Performing transformations ----
      g=V+a1*BV+b1
      xug=a5*UB+b5*BV+c5
      xgr=a7*VR+b7
      xri=a2*RI+b2
      xiz=(a3-a2)*RI+(b3-b2)
      xgi=xgr+xri
      
      return
      end
C***********************************************************************
