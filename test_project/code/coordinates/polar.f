      subroutine polar(iseed,
     &                 minimum_sector_radius,
     &                 maximum_sector_radius,
     &                 angle_covering_sector,
     &                 sector_radius,
     &                 solar_galactocentric_distance,
     &                 scale_length)
C     minimum_sector_radius: minimum radius of the sector; in Kpc from 
C                            the Galactic Center (GC)
C     maximum_sector_radius: maximum radius
C     angle_covering_sector: angle covering the sector in degrees;
C                            from the GC
C     solar_galactocentric_distance: galactocentric distance of the Sun
      implicit none

      external ran       
      real :: ran

C     TODO: this is numberOfStars, take it up as const for all functions
      integer, parameter ::  MAX_STARS_COUNT = 6000000
      real, parameter :: PI = 4.0 * atan(1.0),
     &                   TAU = 2.0 * PI

      integer :: numberOfWDs,
     &           iseed,
     &           wd_index
      real :: minimum_sector_radius,
     &        maximum_sector_radius,
     &        angle_covering_sector,
     &        sector_radius,
     &        scale_length,
     &        angle_covering_sector_in_radians,
     &        squared_sector_radius,
     &        dist,
     &        squared_minimum_sector_radius,
     &        squared_maximum_sector_radius,
     &        random_valid_radius,
     &        zzr,
     &        zzy,
     &        zz,
     &        xx,
     &        xc,
     &        yc,
     &        squared_radii_difference,
     &        sector_diameter
      double precision :: solar_galactocentric_distance

      double precision :: coordinate_R(MAX_STARS_COUNT),
     &                    coordinate_Theta(MAX_STARS_COUNT),
     &                    coordinate_Zcylindr(MAX_STARS_COUNT)
      real :: x(MAX_STARS_COUNT), 
     &        y(MAX_STARS_COUNT)
      integer :: disk_belonging(MAX_STARS_COUNT),
     &           flagOfWD(MAX_STARS_COUNT)

      common /coorcil/ coordinate_R, 
     &                 coordinate_Theta,
     &                 coordinate_Zcylindr
C     TODO: find out what reference frame is used for x, y
      common /plano/ x, y
      common /index/ flagOfWD,
     &               numberOfWDs,
     &               disk_belonging

      squared_sector_radius = sector_radius * sector_radius
      
C     Calculating the angle in the sector
C     -angle_covering_sector / 2 and +angle_covering_sector / 2 degrees
C     and radius between minimum_sector_radius and maximum_sector_radius
      angle_covering_sector_in_radians = angle_covering_sector 
     &                                   * PI / 180.0
      squared_maximum_sector_radius = maximum_sector_radius 
     &                                * maximum_sector_radius
      squared_minimum_sector_radius = minimum_sector_radius 
     &                                * minimum_sector_radius
      squared_radii_difference = squared_maximum_sector_radius 
     &                           - squared_minimum_sector_radius
      sector_diameter = maximum_sector_radius - minimum_sector_radius
                
      do wd_index = 1, numberOfWDs
          do
              coordinate_Theta(wd_index) = (
     &            angle_covering_sector_in_radians
     &            * ran(iseed) - angle_covering_sector_in_radians / 2)
              
              if (coordinate_Theta(wd_index) < 0.0) then
                  coordinate_Theta(wd_index) = (
     &                coordinate_Theta(wd_index) + TAU)
              end if
              
              do 
                  random_valid_radius = minimum_sector_radius 
     &                                  + sector_diameter * ran(iseed)
C                 TODO: find out the meaning of zzy and 0.16
                  zzy = 0.16 * ran(iseed)
C                 TODO: find out the meaning of zzr
                  zzr = exp(-random_valid_radius / scale_length)
                  if (zzy <= zzr) then
                      exit
                  end if
              end do
    
C             TODO: give a good name for zz
              zz = (random_valid_radius - minimum_sector_radius) 
     &             / sector_diameter
C             TODO: find out the meaning of xx
              xx = squared_minimum_sector_radius 
     &             + squared_radii_difference 
     &               * zz
              coordinate_R(wd_index) = sqrt(xx)
              xc = real(coordinate_R(wd_index) 
     &                  * cos(coordinate_Theta(wd_index)))
              yc = real(coordinate_R(wd_index) 
     &                  * sin(coordinate_Theta(wd_index)))
C             TODO: find out the meanng of dist
              dist = real((xc - solar_galactocentric_distance)
     &                    * (xc - solar_galactocentric_distance) 
     &                    + yc * yc)

C             TODO: find out what this means       
C             Sol no hay más que uno
C             TODO: find out the meaning of 0.0000015 const
              if (dist <= squared_sector_radius 
     &                .and. dist >= 0.0000015) then 
                  exit
              end if
          end do

C         TODO: find out the meanng of x and y
          x(wd_index) = xc
          y(wd_index) = yc
      end do
      end subroutine
