C***********************************************************************
C     TODO: rewrite      
      subroutine velh(iseed,numberOfStarsInSample)
C=======================================================================
C
C     This subroutine calculates the heliocentrical velocities. 
C     From the height pattern it calculates the dispersions.
C
C     Revised in 22.09.07 by S. Torres
C
C--------------------------------------------------------------------
C     Input parameters:
C       NOTE: I should give iseed following descriptive name
C       iseed: random number generator parameter
C       numberOfStarsInSample
C=======================================================================
      implicit double precision (a-h,m,o-z)
      
      integer numberOfStars,iseed,numberOfStarsInSample,i,k,numberOfWDs
      double precision uo,vo,wo,a,b,solarGalactocentricDistance,uop,vop,
     &                 yy,gasdev,uom,vom
          
C     ---   Parameters   ---
C--------------------------------------------------------------------
C     a,b: Oort constants (values by Kerr and Lynden-Bell 1986)
C          a=14.4 Km/sKpc
C          b=-12.8 Km/sKpc
C     solarGalactocentricDistance=8.5 Kpc
C--------------------------------------------------------------------
      parameter (numberOfStars=6000000)
      parameter (a=14.4,b=-12.8,solarGalactocentricDistance=8.5d0)
      
C     ---   Dimensiones   ---
C--------------------------------------------------------------------
C     (uu,vv,ww): heliocentric velocities, B3 system
C     sigma(3): dispersion of velocities
C--------------------------------------------------------------------
      double precision uu(numberOfStars),vv(numberOfStars),
     &                 ww(numberOfStars)
      double precision sigma(3)
      double precision coordinate_R(numberOfStars),
     &                 coordinate_Theta(numberOfStars),
     &                 coordinate_Zcylindr(numberOfStars)
      double precision heightPattern(numberOfStars)
      double precision flagOfWD(numberOfStars)
      integer disk_belonging(numberOfStars)
      double precision zz(numberOfStars),zh(numberOfStars)
      
C     ---  Commons  ---
      common /vel/ uu,vv,ww
      common /patron/ heightPattern
      common /coorcil/ coordinate_R,coordinate_Theta,coordinate_Zcylindr
      common /index/ flagOfWD,numberOfWDs,disk_belonging

C     ---  Peculiar solar velocities  ---
C     ---  old values  ---
C      uo=-10.0
C      vo=-5.2
C      wo=-7.2
C     ---  values from Anguiano et al. 2017
      uo=-11.0
      vo=-12.0
      wo=-7.0

C     ---  Making the transfer of heightPattern, of z  ---
      do 2 i=1,numberOfStarsInSample
        zh(i)=heightPattern(i)
        zz(i)=coordinate_Zcylindr(i)
2     continue

      k=0      
      do 3 i=1,numberOfStarsInSample
        if (flagOfWD(i) > 0.95 .and. flagOfWD(i) < 1.05) then
          k=k+1
          heightPattern(k)=zh(i)
          coordinate_Zcylindr(k)=zz(i)
        else
      endif
3     continue

      do 1 i=1,numberOfWDs


C       ---  Calculating the dispersions  ---
C       ----------------------------------------------------------------      
C       ---  THIN DISK  ---
C       ---  model of sigmas depending on h variable
C         sigma(3)=dsqrt(heightPattern(i)/(6.25d-4))          
C         sigma(1)=(dsqrt(2.0d0))*sigma(3)
C         sigma(2)=(dsqrt(0.32+(1.67d-5)*sigma(1)*sigma(1)))*sigma(1)
C       ---  model of constant sigmas
        if (disk_belonging(i) == 1) then
            sigma(1)=32.4
            sigma(2)=23.0
            sigma(3)=18.1 
        else if (disk_belonging(i) == 2) then
            sigma(1)=50.0
            sigma(2)=56.0
            sigma(3)=34.0
        else
            write(6, *) "Error in velh:"
            write(6, *) "    Number of WD's =", numberOfWDs
            write(6, *) "    index =", i
            write(6, *) "    disk_belonging =", disk_belonging(i)
            stop
        end if
        
C       ---   Calling to the function of gaussian distribution  ---
        yy=gasdev(iseed)
        uop=uom(coordinate_R(i),coordinate_Theta(i),a,b,
     &      solarGalactocentricDistance,uo)
        uu(i)=sigma(1)*yy+uop
        yy=gasdev(iseed)
        vop=vom(coordinate_R(i),coordinate_Theta(i),a,b,
     &      solarGalactocentricDistance,vo)
        vv(i)=sigma(2)*yy+vop-sigma(1)*sigma(1)/120.0
        yy=gasdev(iseed)
        ww(i)=sigma(3)*yy+wo
        write(158,*) i,uu(i), vv(i), ww(i)
        write(400,*) uu(i)
        write(401,*) vv(i)
        write(402,*) ww(i)
1     continue

      return
      end
C***********************************************************************


C***********************************************************************
      function uom(r,th,a,b,solarGalactocentricDistance,uo)
C=======================================================================
C
C     Calculating uo taking into account the effect of galactic rotation
C
C=======================================================================
      implicit double precision (a-h,m,o-z)
      
      double precision r,th,a,b,solarGalactocentricDistance,uo,uom
      
      uom=uo+((3.0-(2.0*r)/solarGalactocentricDistance)*a-b)*r*sin(th)
      
      return
      end
C***********************************************************************



C***********************************************************************
      function vom(r,th,a,b,solarGalactocentricDistance,vo)
C=======================================================================
C
C     Calculating vo taking into account the effect of galactic rotation
C
C=======================================================================
      implicit double precision (a-h,m,o-z)
      
      double precision r,th,a,b,solarGalactocentricDistance,vo,vom

      vom = vo + ((3.0 - (2.0 * r) / solarGalactocentricDistance) * a 
     &            - b) * r * cos(th) - (a - b) 
     &                                 * solarGalactocentricDistance

      return
      end
C***********************************************************************
