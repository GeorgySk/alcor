C***********************************************************************
C     TODO: rewrite      
      subroutine coor(solarGalactocentricDistance)
C=======================================================================
C     This subroutine calculates the coordinates and proper motions
C     in the Galactic coordinate system.
C     Also it calculates the right ascension and the declination.
C     galactic coordinates calculated in radians
C     Distances in Kpc
C     Velocities in Km/s
C     Proper motions in  arcsec/yr
C
C     Revised in 22.09.07 by S. Torres
C-------------------------------------------------------------------
C     Input parameters:
C       solarGalactocentricDistance
C       numberOfWDs
C=======================================================================
      implicit real (a-h,m,o-z)
      
C     ---   Definition of variables  ---
      integer numberOfStars,i,numberOfWDs
      real pi,fi
      double precision :: solarGalactocentricDistance,
     &                    ros
      real k,alfag
      double precision :: theta,
     &                    deltag,
     &                    zsdg,
     &                    zcdg,
     &                    zcb,
     &                    zz,
     &                    zsb,
     &                    zzx
      real zsl,zcl,zkr,zkri
      real sinpsi,cospsi,xc,xs
      
C     ---   Parameters   ---
      parameter (numberOfStars=6000000)
      parameter (k=4.74d+3)
      parameter (alfag = 3.35, 
     &           deltag = 0.478d0, 
     &           theta = 2.147d0)
      
C     ---   Dimensions   ---
      double precision :: coordinate_R(numberOfStars),
     &                    coordinate_Theta(numberOfStars),
     &                    coordinate_Zcylindr(numberOfStars)
      double precision :: lgac(numberOfStars),
     &                    bgac(numberOfStars)
      real :: rgac(numberOfStars)
      real mpb(numberOfStars),mpl(numberOfStars),
     &                 vr(numberOfStars)
      real uu(numberOfStars),vv(numberOfStars),
     &                 ww(numberOfStars)
      double precision :: properMotion(numberOfStars),
     &                    rightAscension(numberOfStars),
     &                    declination(numberOfStars)
      real flagOfWD(numberOfStars)
      integer disk_belonging(numberOfStars)
       
C     ---   Commons   ---
      common /coorcil/ coordinate_R,coordinate_Theta,coordinate_Zcylindr
      common /vel/ uu,vv,ww
      common /paral/ rgac
      common /mad/ properMotion,rightAscension,declination
      common /index/ flagOfWD,numberOfWDs,disk_belonging
      common /mopro/ mpb,mpl,vr
      common /lb/ lgac,bgac
      
C     ---   Calculating some parameters  ---      
      pi = 4.0 * atan(1.0)
      fi = 180.0 / pi
      zsdg = dsin(deltag)
      zcdg = dcos(deltag)

C--------------------------------------------------------------------
C     ---   Calculating galactocentric coordinates (r,l,b)  ---
C--------------------------------------------------------------------
      i=1
      do 1 i=1,numberOfWDs
C        ---   Galactic coordinate r (Kpc) ---
        ros = solarGalactocentricDistance * solarGalactocentricDistance
     &        + coordinate_R(i) * coordinate_R(i) 
     &        - 2.d0 * coordinate_R(i) * solarGalactocentricDistance 
     &          * dcos(coordinate_Theta(i))
        rgac(i)=sqrt(ros+(coordinate_Zcylindr(i)*
     &                    coordinate_Zcylindr(i)))
C     ---   Galactic coordinate lgac ---
        ros = dsqrt(ros)
C       TODO: do smth about this asin/acos from 1+o(1)
C          if (abs(coordinate_R(i) * sin(coordinate_Theta(i)) / ros)
C     &        < 1.0) then
C              lgac(i) = dasin(coordinate_R(i) * sin(coordinate_Theta(i))
C     &                        / ros)
C          else
           if ((solarGalactocentricDistance ** 2 + ros**2
     &                         - coordinate_R(i) ** 2)
     &         /(2.d0 * solarGalactocentricDistance * ros) > 1.0) then
            lgac(i) = 0.0
           else
              lgac(i) = dacos((solarGalactocentricDistance ** 2 + ros**2
     &                         - coordinate_R(i) ** 2)
     &                  /(2.0d0 * solarGalactocentricDistance * ros))
           end if
        if (coordinate_R(i)*cos(coordinate_Theta(i)).gt.
     &     solarGalactocentricDistance) then 
          lgac(i)=pi-lgac(i)
        else
          if (sin(coordinate_Theta(i)).lt.0.0) then
            lgac(i)=2.0*pi+lgac(i)
          endif
          continue
        endif

        if (lgac(i) > 2.0 * pi) then
            lgac(i) = lgac(i) - 2.0 * pi
        end if

C        ---   galactic coordinate bgac ---
        zzx = dabs(coordinate_Zcylindr(i) / ros)
        bgac(i) = datan(zzx)
        if (coordinate_Zcylindr(i).lt.0.0) then
          bgac(i)=-bgac(i)
        else
         continue
        endif
C--------------------------------------------------------------------        
C       ---   Calculating the proper motions in galactic coordinates ---
C--------------------------------------------------------------------      
C       ---   Calculating some values that are going to be used later---
        zsl=sin(lgac(i))
        zcl=cos(lgac(i))
        zsb = dsin(bgac(i))
        zcb = dcos(bgac(i))
        zkr=k*rgac(i)
        zkri=(1.0/zkr)
C       ---   Calculating the components of the proper motion ---
        mpl(i)=(-zkri*(zsl/zcb)*uu(i))+(zkri*(zcl/zcb)*vv(i))
        mpb(i)=(-zkri*zcl*zsb*uu(i))+(-zkri*zsb*zsl*vv(i))+
     &         (zkri*zcb*ww(i))
        vr(i)=(zcb*zcl*uu(i))+(zcb*zsl*vv(i))+(zsb*ww(i)) 
        properMotion(i)=sqrt(mpl(i)*mpl(i)+mpb(i)*mpb(i))
C-------------------------------------------------------------------
C       ---   Calculating right ascension and the declination  ---
C-------------------------------------------------------------------      
C       ---   Calculating the declination   ---
        zz = zsdg * zsb + zcdg * zcb * dcos(theta - lgac(i))
        declination(i) = dasin(zz)
C       ---   Calculating the right ascension ---
        xs = ((zcb * dsin(theta - lgac(i))) / dcos(declination(i)))
        xc= ((zcdg*zsb-zsdg*zcb*cos(theta-lgac(i)))/
     &      cos(declination(i)))
C       --Looking at the sign that corresponds to the right ascension---
        if (xs .ge. 0.0) then
          if (xc .ge. 0.0) then 
            rightAscension(i) = asin(xs) + alfag
          else
            rightAscension(i) = acos(xc) + alfag
          endif    
        else
          if (xc .lt. 0.0) then
            rightAscension(i) = pi - asin(xs) + alfag
          else
            rightAscension(i) = 2 * pi + asin(xs) + alfag
          endif
        endif            
        if (rightAscension(i) * fi .gt. 360.0) then
          rightAscension(i) = rightAscension(i) - 2 * pi
        endif
C--------------------------------------------------------------------
C       ---  Calculating the proper motion in ecuatorial coordinates ---
C--------------------------------------------------------------------
        sinpsi=sin(theta)
        cospsi=cos(theta)    
C-------------------------------------------------------------------
C       ---   Calculating the proper motion  (arc sec/yr)   ---
C-------------------------------------------------------------------
  1   continue
           
      return
      end
C***********************************************************************
