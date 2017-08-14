      subroutine errfot(x,xnew,i)
C=======================================================================
C
C     This soubrutine add a photometric error to the ugriz magnitudes
C     using a fifth polynomial fit.
C     See Alberto Rebassa Mansergas et al, 2016, DR12
C
C     Created 03/2016 (ER Cojocaru)
C     Revised 07/2017 (S Torres)
C
C-------------------------------------------------------------------
C
C     input parameters:
C
C     x = magnitude (u,g,r,i,z)
C     i = integer selector between (u,g,r,i,z), 1-5
C
C-------------------------------------------------------------------
C
C     output parameters:
C
C     xnew = magnitude (u,g,r,i,z) with the added error
C
C=======================================================================
      implicit none

      integer i, iseed
      real c0, c1, c2, c3, c4, c5
      real x, sigma, xnew, gasdev
      real temp

      common /ISEED/ iseed


C     coeff_U
      if(i .eq. 1) then
         c0 = -186.8250946999
         c1 = 51.1586639285
         c2 = -5.6026639938
         c3 = 0.3068921035
         c4 = -0.0084118287
         c5 = 0.0000923484
      end if

C     coeff_G
      if(i .eq. 2) then
         c0 = -9.2683243752
         c1 = 3.2245934010
         c2 = -0.4321658909
         c3 = 0.0282152630
         c4 = -0.0009030127
         c5 = 0.0000113837
      end if

C     coeff_R
      if(i .eq. 3) then
         c0 = -13.8539595604
         c1 = 4.6477341652
         c2 = -0.6098692715
         c3 = 0.0393348821
         c4 = -0.0012510254
         c5 = 0.0000157371
      end if

C     coeff_I
      if(i .eq. 4) then
         c0 = -53.9296059608
         c1 = 15.6034383774
         c2 = -1.8090582341
         c3 = 0.1050753752
         c4 = -0.0030573069
         c5 = 0.0000356513
      end if

C     coeff_Z
      if(i .eq. 5) then
         c0 = 196.1956806183
         c1 = -50.6647891998
         c2 = 5.1534942389
         c3 = -0.2570058424
         c4 = 0.0062429576
         c5 = -0.0000584448
      end if

      sigma = c0 + (c1 * x) + (c2 * x ** 2.0) + (c3 * x ** 3.0) 
     &        + (c4 * x ** 4.0) + (c5 * x ** 5.0)

      xnew = x + sigma * gasdev(iseed)

      return
      end
