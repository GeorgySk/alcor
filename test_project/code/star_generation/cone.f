C     Generates stars population inside a square pyramid (pencil/cone) 
C     with its height direction set by longitude and latitude
      subroutine generate_cone_stars(cone_height_longitude,
     &                               cone_height_latitude,
     &                               numberOfStarsInSample,
     &                               iseed,
     &                               thick_disk_stars_fraction,
     &                               thin_disk_age,
     &                               massReductionFactor,
     &                               thick_disk_age,
     &                               thick_disk_sfr_param,
     &                               halo_age,
     &                               halo_stars_formation_time,
     &                               halo_stars_fraction)
          implicit none
          
          external ran
          real ran

          integer, parameter :: MAX_STARS_COUNT = 6000000
          real, parameter :: PI = 4.0 * atan(1.0),
     &                       FI = PI / 180.0,
     &                       DELTA_LATITUDE = 2.64 * FI,
     &                       NORMALIZATION_CONE_HEIGHT = 0.2,
     &                       CONE_HEIGHT = 2.0,
     &                       THIN_DISK_DENSITY = 0.095 * 1.0e9,
     &                       THIN_DISK_SCALEHEIGHT = 0.25,
     &                       THICK_DISK_SCALEHEIGHT = 1.5,
     &                       SOLAR_GALACTOCENTRIC_DISTANCE = 8.5,
     &                       HALO_CORE_RADIUS = 5.0

          integer :: iseed,
     &               numberOfStarsInSample,
     &               numberOfWDs,
     &               stars_count,
     &               disk_belonging(MAX_STARS_COUNT),
     &               flagOfWD(MAX_STARS_COUNT),
     &               thin_disk_stars_count

          real :: cone_height_longitude,
     &            cone_height_latitude,
     &            thin_disk_age,
     &            thick_disk_stars_fraction,
     &            thin_disk_stars_fraction,
     &            massReductionFactor,
     &            starBirthTime(MAX_STARS_COUNT),
     &            m(MAX_STARS_COUNT),
     &            delta_longitude,
     &            min_longitude,
     &            min_latitude,
     &            total_mass,
     &            longitude,
     &            latitude,
     &            density,
     &            random_valid_density,
     &            get_cone_mass,
     &            get_max_density,
     &            distance,
     &            get_density,
     &            generate_star_mass,
     &            thick_disk_max_sfr_relative_time,
     &            opposite_triangle_side,
     &            thick_disk_age,
     &            thick_disk_sfr_param,
     &            thin_disk_normalization_cone_mass,
     &            thin_disk_max_density,
     &            thick_disk_normalization_cone_mass,
     &            thick_disk_max_density,
     &            thick_disk_max_sfr,
     &            time_try,
     &            time_try_sfr,
     &            sfr_try,
     &            thick_disk_birth_init_time,
     &            max_age,
     &            thin_disk_birth_init_time,
     &            halo_age,
     &            halo_stars_formation_time,
     &            halo_age_birth_init_time,
     &            halo_stars_fraction,
     &            get_halo_density,
     &            halo_max_density,
     &            halo_normalization_cone_mass
          
          double precision :: coordinate_Theta(MAX_STARS_COUNT),
     &                        coordinate_R(MAX_STARS_COUNT),
     &                        coordinate_Zcylindr(MAX_STARS_COUNT)

          common /tm/ starBirthTime,
     &                m
          common /coorcil/ coordinate_R,
     &                     coordinate_Theta,
     &                     coordinate_Zcylindr
          common /index/ flagOfWD,numberOfWDs,disk_belonging

          max_age = max(thin_disk_age, thick_disk_age, halo_age)
          thin_disk_birth_init_time = max_age - thin_disk_age
          thick_disk_birth_init_time = max_age - thick_disk_age
          halo_age_birth_init_time = max_age - halo_age

C         This can be easily proved by taking derivative from
C         y = t * exp(-t / tau)
          thick_disk_max_sfr_relative_time = thick_disk_sfr_param

          thick_disk_max_sfr = (thick_disk_max_sfr_relative_time 
     &                          * exp(-thick_disk_max_sfr_relative_time 
     &                                / thick_disk_sfr_param))

          thin_disk_stars_fraction = 1.0 - thick_disk_stars_fraction 
     &                               - halo_stars_fraction

          ! TODO: figure out the case of b = 90º
          delta_longitude = DELTA_LATITUDE / cos(cone_height_latitude)
          min_longitude = cone_height_longitude - delta_longitude / 2.0
          min_latitude = cone_height_latitude  - DELTA_LATITUDE / 2.0

          total_mass = 0.0
          stars_count = 0

          thin_disk_normalization_cone_mass = (massReductionFactor 
     &        * get_cone_mass(
     &              cone_height_longitude, 
     &              cone_height_latitude, 
     &              THIN_DISK_DENSITY * thin_disk_stars_fraction,
     &              THIN_DISK_SCALEHEIGHT))
          thin_disk_max_density = get_max_density(
     &        cone_height_longitude,
     &        cone_height_latitude,
     &        THIN_DISK_SCALEHEIGHT)

          thick_disk_normalization_cone_mass = (massReductionFactor 
     &        * get_cone_mass(
     &              cone_height_longitude, 
     &              cone_height_latitude,
     &              THIN_DISK_DENSITY * thick_disk_stars_fraction,
     &              THICK_DISK_SCALEHEIGHT))
          thick_disk_max_density = get_max_density(
     &        cone_height_longitude,
     &        cone_height_latitude,
     &        THICK_DISK_SCALEHEIGHT)

          halo_normalization_cone_mass = (massReductionFactor 
     &        * get_cone_mass(
     &              cone_height_longitude, 
     &              cone_height_latitude,
     &              THIN_DISK_DENSITY * halo_stars_fraction,
     &              THICK_DISK_SCALEHEIGHT))
          halo_max_density = 1.0 + (SOLAR_GALACTOCENTRIC_DISTANCE 
     &                              / HALO_CORE_RADIUS) ** 2

          if (thin_disk_stars_fraction > 0.0) then
              outer_do1: do
                  stars_count = stars_count + 1          
                  disk_belonging(stars_count) = 1

                  longitude = min_longitude + delta_longitude 
     &                                        * ran(iseed)
                  latitude = min_latitude + DELTA_LATITUDE * ran(iseed)
          
C                 Accepting/rejecting method
                  inner_do1: do
                      distance = CONE_HEIGHT * ran(iseed)
          
                      density = get_density(distance,
     &                                      longitude,
     &                                      latitude,
     &                                      THIN_DISK_SCALEHEIGHT)
                      random_valid_density = thin_disk_max_density 
     &                                       * ran(iseed)

                      if (random_valid_density <= density) then
                          exit
                      end if
                  end do inner_do1

                  m(stars_count) = generate_star_mass(iseed)

                  coordinate_R(stars_count) = opposite_triangle_side(
     &                SOLAR_GALACTOCENTRIC_DISTANCE, 
     &                distance * abs(cos(latitude)), 
     &                longitude)
                  coordinate_Theta(stars_count) = asin(
     &                distance * abs(cos(latitude)) * sin(longitude)
     &                / coordinate_R(stars_count))
                  coordinate_Zcylindr(stars_count) = distance 
     &                                               * sin(latitude)
          
                  if (distance < NORMALIZATION_CONE_HEIGHT) then
                      total_mass = total_mass + m(stars_count)
                  end if
          
                  if (total_mass 
     &                    >= thin_disk_normalization_cone_mass) then
                      exit
                  end if

                  starBirthTime(stars_count) = thin_disk_birth_init_time
     &                                         + thin_disk_age 
     &                                           * ran(iseed)
              end do outer_do1
          end if

          thin_disk_stars_count = stars_count

          if (thick_disk_stars_fraction > 0.0) then
              outer_do2: do
                  stars_count = stars_count + 1
                  disk_belonging(stars_count) = 2
          
                  longitude = min_longitude + delta_longitude 
     &                                        * ran(iseed)
                  latitude = min_latitude + DELTA_LATITUDE * ran(iseed)
          
C                 Accepting/rejecting method
                  inner_do2: do
                      distance = CONE_HEIGHT * ran(iseed)
          
                      density = get_density(distance,
     &                                      longitude,
     &                                      latitude,
     &                                      THICK_DISK_SCALEHEIGHT)          
                      random_valid_density = thick_disk_max_density 
     &                                       * ran(iseed)
          
                      if (random_valid_density <= density) then
                          exit
                      end if
                  end do inner_do2
          
                  m(stars_count) = generate_star_mass(iseed)

                  coordinate_R(stars_count) = opposite_triangle_side(
     &                SOLAR_GALACTOCENTRIC_DISTANCE, 
     &                distance * abs(cos(latitude)), 
     &                longitude)
                  coordinate_Theta(stars_count) = asin(
     &                distance * abs(cos(latitude)) * sin(longitude)
     &                / coordinate_R(stars_count))
                  coordinate_Zcylindr(stars_count) = distance 
     &                                               * sin(latitude)

                  do
                      time_try = thick_disk_age * ran(iseed)
                      time_try_sfr = time_try 
     &                               * exp(-time_try 
     &                                      / thick_disk_sfr_param)
                      sfr_try = thick_disk_max_sfr * ran(iseed)

                      if (sfr_try <= time_try_sfr) then
                          starBirthTime(stars_count) = (
     &                        time_try + thick_disk_birth_init_time)
                          exit
                      end if
                  end do

                  if (stars_count > thin_disk_stars_count * 
     &                    (1.0 + thick_disk_stars_fraction)) then
                      exit
                  end if
              end do outer_do2
          end if

          total_mass = 0.0

          if (halo_stars_fraction > 0.0) then
              outer_do3: do
                  stars_count = stars_count + 1
                  disk_belonging(stars_count) = 3
          
                  longitude = min_longitude + delta_longitude 
     &                                        * ran(iseed)
                  latitude = min_latitude + DELTA_LATITUDE * ran(iseed)
          
C                 Accepting/rejecting method
                  inner_do3: do
                      distance = CONE_HEIGHT * ran(iseed)
          
                      density = get_halo_density(distance,
     &                                           longitude,
     &                                           latitude)          
                      random_valid_density = halo_max_density 
     &                                       * ran(iseed)
          
                      if (random_valid_density <= density) then
                          exit
                      end if
                  end do inner_do3
          
                  m(stars_count) = generate_star_mass(iseed)

                  coordinate_R(stars_count) = opposite_triangle_side(
     &                SOLAR_GALACTOCENTRIC_DISTANCE, 
     &                distance * abs(cos(latitude)), 
     &                longitude)
                  coordinate_Theta(stars_count) = asin(
     &                distance * abs(cos(latitude)) * sin(longitude)
     &                / coordinate_R(stars_count))
                  coordinate_Zcylindr(stars_count) = distance 
     &                                               * sin(latitude)
                  if (distance < NORMALIZATION_CONE_HEIGHT) then
                      total_mass = total_mass + m(stars_count)
                  end if

                  if (total_mass 
     &                    >= thick_disk_normalization_cone_mass) then 
                      exit
                  end if

                  do
                      time_try = thick_disk_age * ran(iseed)
                      time_try_sfr = time_try 
     &                               * exp(-time_try 
     &                                      / thick_disk_sfr_param)
                      sfr_try = thick_disk_max_sfr * ran(iseed)

                      if (sfr_try <= time_try_sfr) then
                          starBirthTime(stars_count) = (
     &                        time_try + thick_disk_birth_init_time)
                          exit
                      end if
                  end do
              end do outer_do3
          end if

          numberOfStarsInSample = stars_count
      end subroutine


      function get_cone_mass(cone_height_longitude, 
     &                       cone_height_latitude,
     &                       density,
     &                       scaleheight) result(mass)
          implicit none
          real, intent(in) :: cone_height_longitude, 
     &                                  cone_height_latitude
          real, parameter :: PI = 4.0 * atan(1.0),
     &                                   FI = PI / 180.0,
     &                                   DELTA_LATITUDE = 2.64 * FI,
     &                                THIN_DISK_DENSITY = 0.095 * 1.0e9,
     &                                   THIN_DISK_SCALEHEIGHT = 0.25
          real :: muted_cone_height_longitude, 
     &                        muted_cone_height_latitude, 
     &                        mass, 
     &                        delta_longitude, 
     &                        min_longitude, 
     &                        max_longitude, 
     &                        min_latitude, 
     &                        max_latitude,
     &                        get_kappa_integral,
     &                        get_iota_integral,
     &                        get_lambda_integral,
     &                        density,
     &                        scaleheight

          ! We need these vars as we don't want to change values of
          ! min_latitude and max_latitude outside this function
          muted_cone_height_longitude = cone_height_longitude
          muted_cone_height_latitude = cone_height_latitude

          ! NOTE: next 2 if-blocks are a workaround
          ! We do this as the picture is symmetrical
          if (muted_cone_height_latitude < 0.0) then
              muted_cone_height_latitude = -muted_cone_height_latitude
          end if
          ! This is done in order to avoid problems with spherical 
          ! coordinates near the pole
          ! Also we assume that angles are in [-90, 90] range
          if (muted_cone_height_latitude > 85.0 * FI) then
              muted_cone_height_latitude = 85.0 * FI
          end if

          delta_longitude = DELTA_LATITUDE 
     &                      / cos(muted_cone_height_latitude)
          min_longitude=muted_cone_height_longitude-delta_longitude/2.0
          max_longitude=muted_cone_height_longitude+delta_longitude/2.0
          min_latitude = muted_cone_height_latitude - DELTA_LATITUDE/2.0
          max_latitude = muted_cone_height_latitude + DELTA_LATITUDE/2.0   
  
          ! Next, we consider 3 cases:
          !   1) Both min_latitude and max_latitude are in [0; PI/2]
          !   2) min_latitude is negative and max_latitude is in [0; PI/2]
          !   3) min_latitude is in [0; PI/2] and max_latitude > PI/2
  
          ! TODO: add explanations about equations and integrals
          ! Case 1:
          if (min_latitude >= 0.0
     &        .and. max_latitude <= PI / 2.0) then
                  mass = density * delta_longitude 
     &                   * scaleheight 
     &                   * (get_kappa_integral(min_latitude,scaleheight) 
     &                      - get_kappa_integral(max_latitude,
     &                                           scaleheight))
          end if
  
          ! Case 2:
          if (min_latitude < 0.0) then
              mass = density * delta_longitude
     &               * (get_iota_integral(abs(min_latitude),
     &                                    scaleheight)
     &                 + get_iota_integral(max_latitude,scaleheight))
          end if
  
          ! Case 3:
          if (max_latitude > PI / 2.0) then
              mass = density * delta_longitude
     &               * (get_lambda_integral((PI - max_latitude),
     &                                      scaleheight)
     &                   + get_lambda_integral(min_latitude,
     &                                         scaleheight))
          end if
      end function


      function get_max_density(longitude,
     &                         latitude,
     &                         scaleheight) result(max_density)
          implicit none
          real, intent(in) :: longitude,
     &                        latitude
          real :: max_density,
     &            density,
     &            distance,
     &            get_density,
     &            scaleheight
          integer :: distance_index
          integer, parameter :: DISTANCE_ITER_COUNT = 1000
          real, parameter :: CONE_HEIGHT = 2.0,
     &                       DISTANCE_INCREMENT 
     &                           = CONE_HEIGHT
     &                             / float(DISTANCE_ITER_COUNT),
     &                       MONTE_CARLO_MAX_SHIFT_FACTOR = 1.1
  
          max_density = 0.0
          do distance_index = 1, DISTANCE_ITER_COUNT
              distance = DISTANCE_INCREMENT * float(distance_index)
  
              density = get_density(distance, 
     &                              longitude,
     &                              latitude,
     &                              scaleheight)
              max_density = max(max_density, density)
          end do
  
          max_density = MONTE_CARLO_MAX_SHIFT_FACTOR * max_density
      end function


      function get_density(distance,longitude,latitude,scaleheight)
     &                                           result(density)
          implicit none
          real, intent(in) :: distance, longitude, latitude
          real :: density, 
     &            pole_projection, 
     &            plane_projection, 
     &            galactocentric_distance,
     &            opposite_triangle_side,
     &            scaleheight
          real, parameter :: SOLAR_GALACTOCENTRIC_DISTANCE = 8.5,
     &                       THICK_DISK_SCALELENGTH = 3.0
  
          pole_projection = distance * abs(sin(latitude))
          plane_projection = distance * abs(cos(latitude))
          
          galactocentric_distance  = opposite_triangle_side(
     &                                   plane_projection, 
     &                                   SOLAR_GALACTOCENTRIC_DISTANCE, 
     &                                   longitude)
  
          density = distance * distance
     &              * exp(-abs(pole_projection) / scaleheight)
     &              * exp(-abs(galactocentric_distance)
     &                     / THICK_DISK_SCALELENGTH)
      end function


      function get_halo_density(distance,
     &                          longitude,
     &                          latitude) result(density)
          implicit none
          real, intent(in) :: distance, 
     &                        longitude,
     &                        latitude
          real :: density,
     &            galactocentric_distance
          real, parameter :: SOLAR_GALACTOCENTRIC_DISTANCE = 8.5,
     &                       THICK_DISK_SCALELENGTH = 3.0,
     &                       HALO_CORE_RADIUS = 5.0  ! kpc
          
          galactocentric_distance  = sqrt(
     &        (SOLAR_GALACTOCENTRIC_DISTANCE - distance * cos(latitude) 
     &                                         * cos(longitude)) ** 2
     &        + (distance * cos(latitude) * sin(latitude)) ** 2
     &        + (distance * sin(latitude)) ** 2)
  
          density = (HALO_CORE_RADIUS ** 2 
     &               + SOLAR_GALACTOCENTRIC_DISTANCE ** 2)
     &              / (HALO_CORE_RADIUS ** 2 + distance ** 2)
      end function


      function generate_star_mass(iseed) result(mass)
          implicit none

          external ran
          real ran

          integer, intent(inout) :: iseed
          real, parameter :: MIN_MASS = 0.4,
     &                            MAX_MASS = 50.0,
     &                            MASS_RANGE = MAX_MASS - MIN_MASS,
     &                            ALPHA_IMF = -2.35, 
     &                            YMAX = MIN_MASS ** ALPHA_IMF
          real :: zy,
     &                        zx, 
     &                        zyimf, 
     &                        mass
         do
              ! Launching(¿lanzamos?) a number (coord y) according to
              ! function of comparison cte. (?)
              zy = YMAX * ran(iseed)
              ! Launching (¿lanzamos?) another number (coord x)
              zx = MIN_MASS + MASS_RANGE * ran(iseed)
              ! Calculating the value of the IMF for m = zx
              zyimf = zx ** ALPHA_IMF
              ! Comparing the IMF with cte function
              if (zy <= zyimf) then
                  exit
              end if
          end do
  
          mass = zx
          
      end function



      function opposite_triangle_side(adjacent_1, 
     &                                adjacent_2, 
     &                                enclosed_angle) result(opposite)
          implicit none
          real, intent(in) :: adjacent_1,
     &                                    adjacent_2,
     &                                    enclosed_angle
          real :: opposite
  
          ! cosines law
          opposite = sqrt(adjacent_1 * adjacent_1
     &                      + adjacent_2 * adjacent_2
     &                      - 2.0 * adjacent_1 * adjacent_2
     &                        * cos(enclosed_angle))
      end function


      ! Calculates an integral that is needed to compute the mass of the
      ! spherical square pyramid. Name kappa is chosen randomly.
      function get_kappa_integral(latitude,scaleheight) result(kappa)
          implicit none
          real, intent(in) :: latitude
          real :: kappa, gamma, scaleheight
          real, parameter :: NORMALIZATION_CONE_HEIGHT=0.2
  
          ! We use this constant just to simplify the next expressions
          gamma = scaleheight / sin(latitude)
  
          kappa = gamma*(gamma-exp(-NORMALIZATION_CONE_HEIGHT / gamma)
     &                         * (gamma + NORMALIZATION_CONE_HEIGHT))
      end function
  
      
      ! Calculates an integral that is needed to compute the mass of the
      ! spherical square pyramid. Case 2. Name iota is chosen randomly.
      function get_iota_integral(latitude,scaleheight) result(iota)
          implicit none
          real, intent(in) :: latitude
          real :: iota, get_kappa_integral,scaleheight
          real, parameter :: THIN_DISK_SCALEHEIGHT = 0.25,
     &                                   NORMALIZATION_CONE_HEIGHT=0.2
  
          iota = (NORMALIZATION_CONE_HEIGHT * NORMALIZATION_CONE_HEIGHT 
     &             / 2.0
     &             - get_kappa_integral(latitude,scaleheight))
     &               * scaleheight
      end function
  
  
      ! Calculates an integral that is needed to compute the mass of the
      ! spherical square pyramid. Case 3. Name lambda is chosen randomly.
      function get_lambda_integral(latitude,scaleheight) result(lambda)
          implicit none
          real, intent(in) :: latitude
          real :: lambda, get_kappa_integral,scaleheight
          real, parameter :: THIN_DISK_SCALEHEIGHT = 0.25,
     &                                   NORMALIZATION_CONE_HEIGHT=0.2
          real, parameter :: PI = 4.0 * atan(1.0)
  
          lambda = -scaleheight
     &               * (get_kappa_integral((PI / 2.0),scaleheight)
     &                  - get_kappa_integral(latitude,scaleheight))
      end function 
