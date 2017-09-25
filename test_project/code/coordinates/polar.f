      subroutine polar(iseed,
     &                 sector_radius,
     &                 solar_galactocentric_distance,
     &                 scale_length)
C     solar_galactocentric_distance: galactocentric distance of the Sun
      implicit none

      external ran       
      real :: ran

C     TODO: this is numberOfStars, take it up as const for all functions
      integer, parameter ::  MAX_STARS_COUNT = 6000000
      real, parameter :: PI = 4.0 * atan(1.0),
     &                   TAU = 2.0 * PI,
     &                   ALPHA_CENTAURI_DISTANCE = 1.5e-6,
C                        Value from "Simulating Gaia performances on
C                        white dwarfs" by Torres et al. 2005
     &                   HALO_CORE_RADIUS = 5.0

      integer :: numberOfWDs,
     &           iseed,
     &           wd_index
      real :: minimum_sector_radius,
     &        maximum_sector_radius,
     &        angle_covering_sector,
     &        sector_radius,
     &        scale_length,
     &        squared_sector_radius,
     &        dist,
     &        squared_minimum_sector_radius,
     &        squared_maximum_sector_radius,
     &        radius_try,
     &        radius_try_distrib,
     &        radial_distrib_random,
     &        random_value,
     &        xc,
     &        yc,
     &        squared_radii_difference,
     &        sector_diameter,
     &        radial_distrib_max,
     &        min_atan,
     &        max_atan,
     &        random_angle
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
      minimum_sector_radius = real(solar_galactocentric_distance) 
     &                        - sector_radius
      maximum_sector_radius = real(solar_galactocentric_distance) 
     &                        + sector_radius
      sector_diameter = maximum_sector_radius - minimum_sector_radius
      squared_maximum_sector_radius = maximum_sector_radius 
     &                                * maximum_sector_radius
      squared_minimum_sector_radius = minimum_sector_radius 
     &                                * minimum_sector_radius
      squared_radii_difference = squared_maximum_sector_radius 
     &                           - squared_minimum_sector_radius
      angle_covering_sector = 2. * asin(
     &    sector_radius / real(solar_galactocentric_distance))
      radial_distrib_max = exp(-minimum_sector_radius / scale_length)
                
      do wd_index = 1, numberOfWDs
          do
              coordinate_Theta(wd_index) = (
     &            angle_covering_sector * ran(iseed) 
     &            - angle_covering_sector / 2)
              
              if (coordinate_Theta(wd_index) < 0.0) then
                  coordinate_Theta(wd_index) = (
     &                coordinate_Theta(wd_index) + TAU)
              end if
              
              if (disk_belonging(wd_index) == 3) then
C                 Inverse transform sampling for halo distribution.
C                 See (4) at "Simulating Gaia performances on white 
C                 dwarfs" by Torres et al. 2005
                  min_atan = atan(minimum_sector_radius 
     &                            / HALO_CORE_RADIUS)
                  max_atan = atan(maximum_sector_radius 
     &                            / HALO_CORE_RADIUS)
                  radius_try = HALO_CORE_RADIUS * tan(
     &                min_atan + ran(iseed) * (max_atan - min_atan))
              else
C                 Accepting-rejecting method
                  do 
                      radius_try = minimum_sector_radius 
     &                             + sector_diameter * ran(iseed)
                      radius_try_distrib = exp(-radius_try 
     &                                         / scale_length)
                      radial_distrib_random = radial_distrib_max 
     &                                        * ran(iseed)
                      if (radial_distrib_random 
     &                        <= radius_try_distrib) then
                          exit
                      end if
                  end do
              end if 

C             Inverse transform sampling method for generating stars
C             uniformly in a circle sector in polar coordinates
              random_value = (radius_try - minimum_sector_radius) 
     &                       / sector_diameter
              coordinate_R(wd_index) = sqrt(
     &            squared_minimum_sector_radius 
     &            + squared_radii_difference * random_value)

              if (disk_belonging(wd_index) == 3) then
                  random_angle = (angle_covering_sector * ran(iseed) 
     &                            - angle_covering_sector / 2)
                  coordinate_Zcylindr(wd_index) = (
     &                coordinate_R(wd_index) * sin(random_angle))
              end if
              
              xc = real(coordinate_R(wd_index) 
     &                  * cos(coordinate_Theta(wd_index)))
              yc = real(coordinate_R(wd_index) 
     &                  * sin(coordinate_Theta(wd_index)))
C             TODO: find out the meanng of dist
              dist = real((xc - solar_galactocentric_distance)
     &                    * (xc - solar_galactocentric_distance) 
     &                    + yc * yc)

              if (dist <= squared_sector_radius 
     &                .and. dist >= ALPHA_CENTAURI_DISTANCE) then 
                  exit
              end if
          end do

C         TODO: find out the meanng of x and y
          x(wd_index) = xc
          y(wd_index) = yc
      end do
      end subroutine
