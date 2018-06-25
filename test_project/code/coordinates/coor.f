      subroutine coor(solarGalactocentricDistance)
C     This subroutine calculates the coordinates and proper motions
C     in the Galactic coordinate system.
C     Also it calculates the right ascension and the declination.
C     Galactic coordinates in radians
C     Distances in Kpc
C     Velocities in Km/s
C     Proper motions in  arcsec/yr
C     More info on conversions at: https://physics.stackexchange.com/questions/88663/converting-between-galactic-and-ecliptic-coordinates
      implicit none

      integer, parameter :: MAX_STARS_COUNT = 6000000
C     TODO: find out the meaning of alfag(ALPHA_G)
C     transverse_velocity = KAPPA * proper_motion * distance
C     4.74 km/sec is one au/year
      real, parameter :: KAPPA = 4.74,
     &                   ALPHA_G = 3.35,
     &                   PI = 4.0 * atan(1.0),
     &                   FI = 180.0 / PI,
     &                   PC_PER_KPC = 1e3
C     THETA - see BK angle at the link above
C     DELTA_G - declination of NGP 
C     TODO: find out why we use NGP declination for J1950 and not J2000 
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
     &        xs,
     &        velocity_by_prop_motion
      double precision :: solarGalactocentricDistance,
     &                    ros,
     &                    cos_bgac,
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
     &        latitude_proper_motion(MAX_STARS_COUNT),
     &        longitude_proper_motion(MAX_STARS_COUNT),
     &        radial_velocity(MAX_STARS_COUNT),
     &        uu(MAX_STARS_COUNT),
     &        vv(MAX_STARS_COUNT),
     &        ww(MAX_STARS_COUNT)
      integer ::  disk_belonging(MAX_STARS_COUNT),
     &            flagOfWD(MAX_STARS_COUNT)

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
      common /mopro/ latitude_proper_motion, 
     &               longitude_proper_motion, 
     &               radial_velocity
      common /lb/ lgac, bgac

      do wd_index = 1, numberOfWDs
C         TODO: rename ros to projection of distance from Sun to WD on
C               galactic plane 
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

C         Unfolding from 0-180 to 0-360
          if (coordinate_Theta(wd_index) > PI) then 
              lgac(wd_index) = 2.0 * PI - lgac(wd_index)
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

          zkr = KAPPA * rgac(wd_index) * PC_PER_KPC
          zkri = 1.0 / zkr
          velocity_by_prop_motion = 1.0 / (KAPPA * rgac(wd_index) 
     &                                     * PC_PER_KPC)

          longitude_proper_motion(wd_index) = real(
     &        velocity_by_prop_motion
     &        * (-uu(wd_index) * sin_lgac + vv(wd_index) * cos_lgac))
          latitude_proper_motion(wd_index) = real(
     &        velocity_by_prop_motion
     &        * (-uu(wd_index) * cos_lgac * sin_bgac
     &           -vv(wd_index) * sin_bgac * sin_lgac
     &           + ww(wd_index) * cos_bgac))
          radial_velocity(wd_index) = real(
     &        uu(wd_index) * cos_bgac * cos_lgac
     &        + vv(wd_index) * cos_bgac * sin_lgac
     &        + ww(wd_index) * sin_bgac)
          properMotion(wd_index) = sqrt(
     &        longitude_proper_motion(wd_index) 
     &        * longitude_proper_motion(wd_index)
     &        + latitude_proper_motion(wd_index) 
     &          * latitude_proper_motion(wd_index))

C         More info at the link above
          declination(wd_index) = dasin(
     &        SIN_DELTA_G * sin_bgac + COS_DELTA_G * cos_bgac 
     &                                 * dcos(THETA - lgac(wd_index)))

C         TODO: give better names for xs and xc
C         These vars are used for conversion from galactic to equatorial 
C         coordinates. More info at the link above
          xs = real((cos_bgac * dsin(THETA - lgac(wd_index))) 
     &              / dcos(declination(wd_index)))
          xc =  real((COS_DELTA_G * sin_bgac 
     &                - SIN_DELTA_G * cos_bgac 
     &                  * cos(THETA - lgac(wd_index)))
     &               / cos(declination(wd_index)))

C         Looking at the sign that corresponds to the right ascension
C         TODO: find out what is going on here
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
