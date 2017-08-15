      subroutine coor(solarGalactocentricDistance)
C     This subroutine calculates the coordinates and proper motions
C     in the Galactic coordinate system.
C     Also it calculates the right ascension and the declination.
C     Galactic coordinates in radians
C     Distances in Kpc
C     Velocities in Km/s
C     Proper motions in  arcsec/yr
      implicit none
      
C     TODO: this is numberOfStars, take it up as const for all functions
      integer, parameter :: MAX_STARS_COUNT = 6000000
C     TODO: find out the meaning of k(K_CONST) 
C     TODO: find out the meaning of alfag(ALPHA_G)
      real, parameter :: K_CONST = 4.74d+3,
     &                   ALPHA_G = 3.35,
     &                   PI = 4.0 * atan(1.0),
     &                   FI = 180.0 / PI
C     TODO: find out the meaning of deltag(DELTA_G) 
C     TODO: find out the meaning of theta(THETA)
      double precision, parameter :: DELTA_G = 0.478d0, 
     &                               THETA = 2.147d0,
     &                               SIN_DELTA_G = dsin(DELTA_G),
     &                               COS_DELTA_G = dcos(DELTA_G)

      integer :: wd_index,
     &           numberOfWDs
      real :: sin_lgac,
     &        cos_lgac,
     &        zkr,
     &        zkri,
     &        xc,
     &        xs
      double precision :: solarGalactocentricDistance,
     &                    ros,
     &                    cos_bgac,
     &                    zz,
     &                    sin_bgac,
     &                    zzx
      double precision :: coordinate_R(MAX_STARS_COUNT),
     &                    coordinate_Theta(MAX_STARS_COUNT),
     &                    coordinate_Zcylindr(MAX_STARS_COUNT),
     &                    lgac(MAX_STARS_COUNT),
     &                    bgac(MAX_STARS_COUNT),
     &                    properMotion(MAX_STARS_COUNT),
     &                    rightAscension(MAX_STARS_COUNT),
     &                    declination(MAX_STARS_COUNT)
      real :: rgac(MAX_STARS_COUNT),
     &        mpb(MAX_STARS_COUNT),
     &        mpl(MAX_STARS_COUNT),
     &        vr(MAX_STARS_COUNT),
     &        uu(MAX_STARS_COUNT),
     &        vv(MAX_STARS_COUNT),
     &        ww(MAX_STARS_COUNT),
     &        flagOfWD(MAX_STARS_COUNT)
      integer disk_belonging(MAX_STARS_COUNT)

      common /coorcil/ coordinate_R, 
     &                 coordinate_Theta,
     &                 coordinate_Zcylindr
      common /vel/ uu, vv, ww
      common /paral/ rgac
      common /mad/ properMotion,
     &             rightAscension,
     &             declination
      common /index/ flagOfWD, 
     &               numberOfWDs,
     &               disk_belonging
      common /mopro/ mpb, mpl, vr
      common /lb/ lgac, bgac

C     Calculating galactocentric coordinates (r,l,b)
      do wd_index = 1, numberOfWDs
C         Galactic coordinate r (Kpc)
C         TODO: find out the meaning of ros
          ros = solarGalactocentricDistance 
     &          * solarGalactocentricDistance
     &          + coordinate_R(wd_index) * coordinate_R(wd_index) 
     &          - 2.d0 * coordinate_R(wd_index) 
     &            * solarGalactocentricDistance 
     &            * dcos(coordinate_Theta(wd_index))
          rgac(wd_index) = real(
     &        sqrt(ros + coordinate_Zcylindr(wd_index)
     &                   * coordinate_Zcylindr(wd_index)))
C         Galactic coordinate lgac
          ros = dsqrt(ros)

C         TODO: check if precision can be lost here
          if ((solarGalactocentricDistance ** 2 + ros**2
     &         - coordinate_R(wd_index) ** 2)
     &        / (2.d0 * solarGalactocentricDistance * ros) > 1.0) then
              lgac(wd_index) = 0.0
          else
              lgac(wd_index) = dacos((solarGalactocentricDistance ** 2 
     &                         + ros**2
     &                         - coordinate_R(wd_index) ** 2)
     &                               / (2.0d0 * ros
     &                                  * solarGalactocentricDistance))
          end if

          if (coordinate_R(wd_index) * cos(coordinate_Theta(wd_index))
     &            > solarGalactocentricDistance) then 
              lgac(wd_index) = PI - lgac(wd_index)
          else
              if (sin(coordinate_Theta(wd_index)) < 0.0) then
                  lgac(wd_index)=2.0*PI+lgac(wd_index)
              end if
              continue
          end if

          if (lgac(wd_index) > 2.0 * PI) then
              lgac(wd_index) = lgac(wd_index) - 2.0 * PI
          end if

C         Galactic coordinate bgac ---
C         TODO: find out the meaning of zzx
          zzx = dabs(coordinate_Zcylindr(wd_index) / ros)
          bgac(wd_index) = datan(zzx)
          if (coordinate_Zcylindr(wd_index) < 0.0) then
              bgac(wd_index)=-bgac(wd_index)
          else
              continue
          end if
       
C         Calculating the proper motions in galactic coordinates
          sin_lgac = real(sin(lgac(wd_index)))
          cos_lgac = real(cos(lgac(wd_index)))
          sin_bgac = dsin(bgac(wd_index))
          cos_bgac = dcos(bgac(wd_index))
C         TODO: find out the meaning of zkr
          zkr = K_CONST * rgac(wd_index)
C         TODO: find out the meaning of zkr
          zkri = 1.0 / zkr

C         Calculating the components of the proper motion
          mpl(wd_index) = real(
     &        (-zkri * (sin_lgac / cos_bgac) * uu(wd_index))
     &        + (zkri * (cos_lgac / cos_bgac) * vv(wd_index)))
          mpb(wd_index) = real(
     &        (-zkri * cos_lgac * sin_bgac * uu(wd_index))
     &        + (-zkri * sin_bgac * sin_lgac * vv(wd_index))
     &        + (zkri * cos_bgac * ww(wd_index)))
          vr(wd_index) = real(
     &        (cos_bgac * cos_lgac * uu(wd_index))
     &        + (cos_bgac * sin_lgac * vv(wd_index))
     &        + (sin_bgac * ww(wd_index)))
          properMotion(wd_index) = sqrt(mpl(wd_index) * mpl(wd_index)
     &                                  + mpb(wd_index) * mpb(wd_index))

C         Calculating right ascension and the declination      
C         Calculating the declination
C         TODO: find out the meaning of zz
          zz = SIN_DELTA_G * sin_bgac 
     &         + COS_DELTA_G * cos_bgac * dcos(THETA - lgac(wd_index))
          declination(wd_index) = dasin(zz)

C         Calculating the right ascension
C         TODO: find out the meaning of xs
          xs = real((cos_bgac * dsin(THETA - lgac(wd_index))) 
     &              / dcos(declination(wd_index)))
C         TODO: find out the meaning of xc
          xc =  real((COS_DELTA_G * sin_bgac 
     &                - SIN_DELTA_G * cos_bgac 
     &                  * cos(THETA - lgac(wd_index)))
     &               / cos(declination(wd_index)))

C         Looking at the sign that corresponds to the right ascension
          if (xs >= 0.0) then
              if (xc >= 0.0) then 
                  rightAscension(wd_index) = asin(xs) + ALPHA_G
              else
                  rightAscension(wd_index) = acos(xc) + ALPHA_G
              end if    
          else
              if (xc < 0.0) then
                  rightAscension(wd_index) = PI - asin(xs) + ALPHA_G
              else
                  rightAscension(wd_index) = 2 * PI + asin(xs) + ALPHA_G
              end if
          end if            
          if (rightAscension(wd_index) * FI > 360.0) then
              rightAscension(wd_index) = rightAscension(wd_index) 
     &                                   - 2 * PI
          end if   
      end do
      end subroutine
