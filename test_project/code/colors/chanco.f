C     TODO: rename to change_color
      subroutine chanco(V, UB, BV, VR, RI, g, xug, xgr, xri, xiz, xgi)
C         Transformation between the Johnson-Cousins UBVRI photometry 
C         system and the SDSS ugriz system. Follow the eq. 1 to 8 from 
C         Jordi,Grebel & Ammon, 2006, A&A,460
          implicit none
          
C         Used parameters (Table 3) 
          real, parameter :: a1 = 0.630, 
     &                       a2 = 1.007, 
     &                       a3 = 1.584, 
     &                       a5 = 0.750,  
     &                       a7 = 1.646, 
     &                       b1 = -0.124,
     &                       b2 = -0.236,
     &                       b3 = -0.386,
     &                       b5 = 0.770,
     &                       b7 = -0.139,
     &                       c5 = 0.720
C         TODO: give better names
          real :: V, UB, BV, VR, RI, g, xug, xgr, xri, xiz, xgi

          g = V + a1 * BV + b1
          xug = a5 * UB + b5 * BV + c5
          xgr = a7 * VR + b7
          xri = a2 * RI + b2
          xiz = (a3 - a2) * RI + (b3 - b2)
          xgi = xgr + xri
      end subroutine
