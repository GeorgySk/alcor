C     Generates stars population inside a square pyramid (pencil/cone) 
C     with its height direction set by longitude and latitude
      subroutine generate_cone_stars(cone_height_longitude,
     &                               cone_height_latitude,
     &                               numberOfStarsInSample,
     &                               iseed,
     &                               kinematicModel,
     &                               galacticDiskAge,
     &                               min_longitude,
     &                               max_longitude,
     &                               min_latitude,
     &                               max_latitude)
          implicit none
          
          external ran
          real ran

          real :: cone_height_longitude,
     &                        cone_height_latitude,
     &                        galacticDiskAge
          integer :: iseed,
     &               numberOfStarsInSample,
     &               kinematicModel

          integer, parameter :: numberOfStars = 6000000
          real :: starBirthTime(numberOfStars),
     &                        m(numberOfStars),
     &                        heightPattern(numberOfStars),
     &                        flagOfWD(numberOfStars)
          double precision :: coordinate_Theta(numberOfStars),
     &                        coordinate_R(numberOfStars),
     &                        coordinate_Zcylindr(numberOfStars)
          integer :: numberOfWDs,
     &               disk_belonging(numberOfStars)

          real, parameter :: PI = 4.0 * atan(1.0),
     &                                   FI = PI / 180.0,
     &                                   DELTA_LATITUDE = 2.64 * FI,
     &                                  NORMALIZATION_CONE_HEIGHT = 0.2,
     &                                   CONE_HEIGHT = 2.0,
     &                                THIN_DISK_DENSITY = 0.095 * 1.0e9,
     &                         THIN_DISK_STARS_COUNT_FRACTION = 0.92,
     &                         THICK_DISK_STARS_COUNT_FRACTION = 1.0 
     &                                 - THIN_DISK_STARS_COUNT_FRACTION,
     &                                  THIN_DISK_SCALEHEIGHT = 0.25,
     &                                  THICK_DISK_SCALEHEIGHT = 1.5,
     &                                  MASS_REDUCTION_FACTOR = 0.003,
     &                               SOLAR_GALACTOCENTRIC_DISTANCE = 8.5
          real :: delta_longitude,
     &                        min_longitude,
     &                        max_longitude,
     &                        longitude_range,
     &                        min_latitude,
     &                        max_latitude,
     &                        normalization_cone_mass,
     &                        max_density,
     &                        total_mass,
     &                        longitude,
     &                        latitude,
     &                        density,
     &                        random_valid_density,
     &                        get_cone_mass,
     &                        get_max_density,
     &                        distance,
     &                        get_density,
     &                        generate_star_mass,
     &                        tmax,
     &                        tmdisk = 10.0,
     &                        tau = 2.0,
     &                        ttry,
     &                        ttdisk = 12.0,
     &                        ft,
     &                        fz,
     &                        opposite_triangle_side
          integer :: stars_count

          common /tm/ starBirthTime,
     &                m
          common /coorcil/ coordinate_R,
     &                     coordinate_Theta,
     &                     coordinate_Zcylindr
          common /patron/ heightPattern
          common /index/ flagOfWD,numberOfWDs,disk_belonging

          ! NOTE: this can be infinity
          delta_longitude = DELTA_LATITUDE / cos(cone_height_latitude)
          min_longitude = cone_height_longitude  - delta_longitude / 2.0
          max_longitude = cone_height_longitude  + delta_longitude / 2.0
          longitude_range = max_longitude - min_longitude
          min_latitude = cone_height_latitude  - DELTA_LATITUDE / 2.0
          max_latitude = cone_height_latitude + DELTA_LATITUDE / 2.0

          total_mass = 0.0
          stars_count = 0

          if (kinematicModel == 1) then
            normalization_cone_mass=get_cone_mass(cone_height_longitude, 
     &                                            cone_height_latitude,
     &                                            THIN_DISK_DENSITY,
     &                                            THIN_DISK_SCALEHEIGHT)
            normalization_cone_mass = normalization_cone_mass 
     &                                * MASS_REDUCTION_FACTOR

            max_density = get_max_density(cone_height_longitude,
     &                                    cone_height_latitude,
     &                                    THIN_DISK_SCALEHEIGHT)

              outer_do1: do
                  stars_count = stars_count + 1
      
                  if (stars_count > numberOfStars) then
                      open(unit=725,file='SKIPPED_PLATES.txt',
     &                     access='append')
                      write(unit=725,fmt=*) cone_height_longitude,
     &                                      cone_height_latitude
                      stop
                  end if

                  disk_belonging(stars_count) = 1
                  heightPattern(stars_count) = THIN_DISK_SCALEHEIGHT
      
C                 assuming uniform distribution
                  longitude = min_longitude + delta_longitude*ran(iseed)
                  latitude = min_latitude + DELTA_LATITUDE * ran(iseed)
      
                  inner_do1: do
                      ! Distance between Sun and random point in the cone
                      distance = CONE_HEIGHT * ran(iseed)
      
                      density = get_density(distance,longitude,latitude,
     &                              THIN_DISK_SCALEHEIGHT)
      
                      ! Monte-Carlo accepting/rejecting method
                      ! QUESTION: can this be outside the inner loop?
                      random_valid_density = max_density * ran(iseed)
      
                      if (random_valid_density <= density) exit
                  end do inner_do1
      
C                 mass
                  m(stars_count) = generate_star_mass(iseed)
C                 converting from galactic to cyl.galactocentric
                  coordinate_R(stars_count) 
     &                = opposite_triangle_side(
     &                      SOLAR_GALACTOCENTRIC_DISTANCE, 
     &                      distance * abs(cos(latitude)), 
     &                      longitude)
                  coordinate_Theta(stars_count) 
     &               = asin(distance*abs(cos(latitude))
     &                 *sin(longitude)/coordinate_R(stars_count))
                  coordinate_Zcylindr(stars_count) = distance 
     &                                               * sin(latitude)
      
                  if(distance < NORMALIZATION_CONE_HEIGHT) then
                      total_mass = total_mass + m(stars_count)
                  endif
      
                  if (total_mass >= normalization_cone_mass) exit
      
                  ! Assuming constant star formation rate
                  starBirthTime(stars_count) =galacticDiskAge*ran(iseed)
              end do outer_do1

          else if (kinematicModel == 2) then

            normalization_cone_mass=get_cone_mass(cone_height_longitude, 
     &                                            cone_height_latitude,
     &               THIN_DISK_DENSITY * THIN_DISK_STARS_COUNT_FRACTION,
     &                                            THIN_DISK_SCALEHEIGHT)
            normalization_cone_mass = normalization_cone_mass 
     &                                * MASS_REDUCTION_FACTOR

            max_density = get_max_density(cone_height_longitude,
     &                                    cone_height_latitude,
     &                                    THIN_DISK_SCALEHEIGHT)
              !running for 92% of thin disk stars:
              outer_do2: do
                  stars_count = stars_count + 1
      
                  disk_belonging(stars_count) = 1
                  heightPattern(stars_count) = THIN_DISK_SCALEHEIGHT
      
C                 assuming uniform distribution
                  longitude = min_longitude + delta_longitude*ran(iseed)
                  latitude = min_latitude + DELTA_LATITUDE * ran(iseed)
      
                  inner_do2: do
                      ! Distance between Sun and random point in the cone
                      distance = CONE_HEIGHT * ran(iseed)
      
                      density = get_density(distance,
     &                                      longitude,
     &                                      latitude,
     &                                      THIN_DISK_SCALEHEIGHT)
      
                      ! Monte-Carlo accepting/rejecting method
                      ! QUESTION: can this be outside the inner loop?
                      random_valid_density = max_density * ran(iseed)
      
                      if (random_valid_density <= density) exit
                  end do inner_do2
      
C                 mass
                  m(stars_count) = generate_star_mass(iseed)
C                 converting from galactic to cyl.galactocentric
                  coordinate_R(stars_count) 
     &                = opposite_triangle_side(
     &                      SOLAR_GALACTOCENTRIC_DISTANCE, 
     &                      distance * abs(cos(latitude)), 
     &                      longitude)
                  coordinate_Theta(stars_count) 
     &                = asin(distance*abs(cos(latitude))
     &                  *sin(longitude)/coordinate_R(stars_count))
                  coordinate_Zcylindr(stars_count) = distance 
     &                                               * sin(latitude)
      
                  if(distance < NORMALIZATION_CONE_HEIGHT) then
                      total_mass = total_mass + m(stars_count)
                  endif
      
                  if (total_mass >= normalization_cone_mass) exit
      
                  ! Assuming constant star formation rate
                  starBirthTime(stars_count) =galacticDiskAge*ran(iseed)
              end do outer_do2

            total_mass = 0.0

            normalization_cone_mass=get_cone_mass(cone_height_longitude, 
     &                                            cone_height_latitude,
     &              THIN_DISK_DENSITY * THICK_DISK_STARS_COUNT_FRACTION,
     &                                          THICK_DISK_SCALEHEIGHT)
            normalization_cone_mass = normalization_cone_mass 
     &                                * MASS_REDUCTION_FACTOR

            max_density = get_max_density(cone_height_longitude,
     &                                    cone_height_latitude,
     &                                    THICK_DISK_SCALEHEIGHT)

              !running for 8% of thin disk stars:
              outer_do3: do
                  stars_count = stars_count + 1
      
                  disk_belonging(stars_count) = 2
                  heightPattern(stars_count) = THICK_DISK_SCALEHEIGHT
      
C                 assuming uniform distribution
                  longitude = min_longitude + delta_longitude*ran(iseed)
                  latitude = min_latitude + DELTA_LATITUDE * ran(iseed)
      
                  inner_do3: do
                      ! Distance between Sun and random point in the cone
                      distance = CONE_HEIGHT * ran(iseed)
      
                      density = get_density(distance,longitude,latitude,
     &                                      THICK_DISK_SCALEHEIGHT)
      
                      ! Monte-Carlo accepting/rejecting method
                      ! QUESTION: can this be outside the inner loop?
                      random_valid_density = max_density * ran(iseed)
      
                      if (random_valid_density <= density) exit
                  end do inner_do3
      
C                 mass
                  m(stars_count) = generate_star_mass(iseed)
C                 converting from galactic to cyl.galactocentric
                  coordinate_R(stars_count) 
     &                = opposite_triangle_side(
     &                      SOLAR_GALACTOCENTRIC_DISTANCE, 
     &                      distance * abs(cos(latitude)), 
     &                      longitude)
                  coordinate_Theta(stars_count) 
     &                = asin(distance*abs(cos(latitude))
     &                  *sin(longitude)/coordinate_R(stars_count))
                  coordinate_Zcylindr(stars_count) = distance 
     &                                               * sin(latitude)
      
                  if(distance < NORMALIZATION_CONE_HEIGHT) then
                      total_mass = total_mass + m(stars_count)
                  endif
      
                  if (total_mass >= normalization_cone_mass) exit
      
                  ! Assuming constant star formation rate
                  starBirthTime(stars_count) =galacticDiskAge*ran(iseed)
                  tmax = tmdisk * exp(-tmdisk/tau)
 33               ttry = ttdisk * ran(iseed)
                  ft = ttry * exp(-ttry/tau)
                  fz = tmax * ran(iseed)
                  if (fz .le. ft) then
                      starBirthTime(stars_count) = ttry
                  else
                      goto 33
                  end if
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
     &                        longitude_range, 
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
          longitude_range = max_longitude - min_longitude
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
                  mass = density * longitude_range 
     &                   * scaleheight 
     &                   * (get_kappa_integral(min_latitude,scaleheight) 
     &                      - get_kappa_integral(max_latitude,
     &                                           scaleheight))
          end if
  
          ! Case 2:
          if (min_latitude < 0.0) then
              mass = density * longitude_range
     &               * (get_iota_integral(abs(min_latitude),
     &                                    scaleheight)
     &                 + get_iota_integral(max_latitude,scaleheight))
          end if
  
          ! Case 3:
          if (max_latitude > PI / 2.0) then
              mass = density * longitude_range
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
     &                                    latitude
          real :: max_density,
     &                        density,
     &                        distance,
     &                        get_density,
     &                        scaleheight
          integer :: distance_index
          integer, parameter :: DISTANCE_ITER_COUNT = 1000
          real, parameter :: CONE_HEIGHT = 5.0,
     &                                   DISTANCE_INCREMENT=CONE_HEIGHT
     &                                   / float(DISTANCE_ITER_COUNT),
     &                               MONTE_CARLO_MAX_SHIFT_FACTOR = 1.1
  
          do distance_index = 1, DISTANCE_ITER_COUNT
              distance = DISTANCE_INCREMENT * float(distance_index)
  
              density = get_density(distance, 
     &                              longitude,
     &                              latitude,
     &                              scaleheight)
              max_density = max1(max_density, density)
          end do
  
          max_density = MONTE_CARLO_MAX_SHIFT_FACTOR * max_density
      end function


      function get_density(distance,longitude,latitude,scaleheight)
     &                                           result(density)
          implicit none
          real, intent(in) :: distance, longitude, latitude
          real :: density, 
     &                        pole_projection, 
     &                        plane_projection, 
     &                        galactocentric_distance,
     &                        opposite_triangle_side,
     &                        scaleheight
          real, parameter :: 
     &          SOLAR_GALACTOCENTRIC_DISTANCE = 8.5,
     &          THICK_DISK_SCALELENGTH = 3.0
  
          pole_projection = distance * abs(sin(latitude))
          plane_projection = distance * abs(cos(latitude))
          
          galactocentric_distance  = opposite_triangle_side(
     &                                  plane_projection, 
     &                                  SOLAR_GALACTOCENTRIC_DISTANCE, 
     &                                  longitude)
  
          density = distance * distance
     &            * exp(-abs(pole_projection) / scaleheight)
     &            * exp(-abs(galactocentric_distance)
     &                         / THICK_DISK_SCALELENGTH)
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
