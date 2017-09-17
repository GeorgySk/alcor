      subroutine gen(iseed,
     &               parameterOfSFR,
     &               radius,
     &               numberOfStarsInSample,
     &               thin_disk_age,
     &               burst_age,
     &               massReductionFactor,
     &               thick_disk_stars_fraction,
     &               halo_stars_fraction,
     &               thick_disk_age,
     &               thick_disk_sfr_param,
     &               halo_age,
     &               halo_stars_formation_time)
C     Divides the SFR in intervals of time. In each interval, the total 
C     mass of stars is distributed. The mass of each star follows the 
C     distribution law given by Initial Mass Function (IMF). The birth 
C     time is also calculated, from which the scale height and finally 
C     cylindrical z-coordinate are determined.
      implicit none

      external ran
      real ran

      integer, parameter :: MAX_STARS_COUNT = 6000000,
     &                      BINS_COUNT = 5000
      real, parameter :: THIN_DISK_SCALE_HEIGHT_KPC = 0.250,
     &                   THICK_DISK_SCALE_HEIGHT_KPC = 0.900,
     &                   BURST_FORMATION_FACTOR = 5.0,
     &                   HALO_STARS_BIRTH_END_TIME = 1.0,
     &                   HALO_DISTANCE_PARAM = 5.0,
     &                   PI = 4.0 * atan(1.0)
      integer ::  iseed,
     &            numberOfStarsInSample, 
     &            bin_index, 
     &            stars_count, 
     &            in, 
     &            numberOfWDs,
     &            flagOfWD(MAX_STARS_COUNT)
      real :: z_coordinate,                radius,
     &        time_try_sfr,                sfr_try,
     &        burst_age,                   burst_init_time, 
     &        current_bin_init_time,       random_value,
     &        time_try,                    thick_disk_stars_fraction,
     &        parameterOfSFR,              areaOfSector, 
     &        thin_disk_age,               normalization_const,
     &        time_increment,              massReductionFactor,
     &        psi,                         star_mass,
     &        total_generated_mass_in_bin, mrep,
     &        halo_stars_fraction,         burst_mrep,
     &        normal_mrep,                 thick_disk_max_sfr,
     &        thin_disk_birth_init_time,   thick_disk_age,
     &        thick_disk_max_sfr_relative_time,    thick_disk_sfr_param,
     &        halo_age,                    halo_stars_formation_time,
     &        thick_disk_birth_init_time,  max_age,
     &        halo_birth_init_time
      real :: mass_from_Salpeter_IMF, 
     &        get_normalization_const
      double precision :: coordinate_Theta(MAX_STARS_COUNT),
     &                    coordinate_R(MAX_STARS_COUNT),
     &                    coordinate_Zcylindr(MAX_STARS_COUNT)
C     TODO: rename this to main_sequence_mass? 
      real :: m(MAX_STARS_COUNT), 
     &        starBirthTime(MAX_STARS_COUNT),
     &        scale_height
      integer :: disk_belonging(MAX_STARS_COUNT)

      common /tm/ starBirthTime, 
     &            m
      common /coorcil/ coordinate_R,
     &                 coordinate_Theta,
     &                 coordinate_Zcylindr
      common /index/ flagOfWD,
     &               numberOfWDs,
     &               disk_belonging

C     TODO: find out the meaning of psi, mrep and 1.0e6
      normalization_const = get_normalization_const(parameterOfSFR, 
     &                                              thin_disk_age)
      time_increment = thin_disk_age / BINS_COUNT
      psi = normalization_const * time_increment
      areaOfSector = PI * radius ** 2
      mrep = psi * areaOfSector * 1.0e6
      mrep = mrep * massReductionFactor
      normal_mrep = mrep
      burst_mrep = mrep * BURST_FORMATION_FACTOR

      max_age = max(thin_disk_age, thick_disk_age, halo_age)
      burst_init_time = max_age - burst_age
      thin_disk_birth_init_time = max_age - thin_disk_age
      thick_disk_birth_init_time = max_age - thick_disk_age
      halo_birth_init_time = max_age - halo_age
C     This can be easily proved by taking derivative from
C     y = t * exp(-t / tau)
      thick_disk_max_sfr_relative_time = thick_disk_sfr_param
      thick_disk_max_sfr = (thick_disk_max_sfr_relative_time 
     &                      * exp(-thick_disk_max_sfr_relative_time 
     &                              / thick_disk_sfr_param))

      write(6,*) '      Mass reduction factor = ', massReductionFactor

      stars_count = 0
     
      do bin_index = 1, BINS_COUNT
          total_generated_mass_in_bin = 0.0

          current_bin_init_time = (thin_disk_birth_init_time 
     &                             + float(bin_index - 1) 
     &                               * time_increment)
    
          if (current_bin_init_time >= burst_init_time) then
              mrep = burst_mrep
          else
              mrep = normal_mrep
          end if

          do
              star_mass = mass_from_Salpeter_IMF(iseed)
              stars_count = stars_count + 1

              if (stars_count == MAX_STARS_COUNT) then 
                  write(6,*) '     ***  Dimension exceeded   ***'
                  write(6,*) '***  Increase the reduction factor   ***'
                  stop
              end if

              m(stars_count) = star_mass       

C             disk_belonging = 1 (thin disk), = 2 (thick disk)
              random_value = ran(iseed)
              if (random_value <= thick_disk_stars_fraction) then
                  disk_belonging(stars_count) = 2
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
              else if (random_value > thick_disk_stars_fraction 
     &                 .and. random_value <= thick_disk_stars_fraction 
     &                                       + halo_stars_fraction) then
                  disk_belonging(stars_count) = 3
                  starBirthTime(stars_count) = halo_birth_init_time
     &                + halo_stars_formation_time * ran(iseed)
              else
                  disk_belonging(stars_count) = 1
                  starBirthTime(stars_count) = thin_disk_birth_init_time 
     &                + float(bin_index - 1) * time_increment 
     &                + time_increment * ran(iseed)
                  total_generated_mass_in_bin = 
     &                total_generated_mass_in_bin + star_mass
              end if

C             TODO: add halo stars here
              if (disk_belonging(stars_count) == 1) then
                  scale_height = THIN_DISK_SCALE_HEIGHT_KPC
              else if (disk_belonging(stars_count) == 2) then
                  scale_height = THICK_DISK_SCALE_HEIGHT_KPC
              end if

C             Inverse transform sampling for y = exp(-z / H)
              do
                  random_value = ran(iseed)
                  if (random_value /= 0.0) then
                      exit
                  end if
              end do
              z_coordinate = -scale_height * log(random_value)
C             Assigning random sign
              in = int(2.0 * ran(iseed))
              coordinate_Zcylindr(stars_count) = z_coordinate 
     &                                           * dfloat(1 - 2 * in)

              if (total_generated_mass_in_bin >= mrep) then
                  exit
              end if
          end do
      end do
        
      numberOfStarsInSample = stars_count
      write(6,*) '     Number of stars in sample=',numberOfStarsInSample
      
      end subroutine


      function get_normalization_const(parameterOfSFR, 
     &                                 thin_disk_age)
C         Calculating the normalization constant of the SFR
          implicit none

C         TODO: find out what this sigma is
          real, parameter :: sigma = 51.0
          real :: parameterOfSFR,
     &            thin_disk_age,
     &            get_normalization_const

          get_normalization_const = sigma 
     &                              / (parameterOfSFR 
     &                                 * (1.0 - exp(-thin_disk_age 
     &                                               / parameterOfSFR)))
      end function


C     TODO: change this to inverse transform sampling
      function mass_from_Salpeter_IMF(iseed)
C         Calculating the mass following the IMF by Salpeter 55.
          implicit none
          
          external ran
          real ran

          real, parameter :: mmin = 0.4,
     &                       mmax = 50.0

          integer :: iseed,
     &               ISEED1,
     &               ISEED2
          real :: mass_from_Salpeter_IMF,
     &            fractionOfDB,
     &            thin_disk_age,
     &            parameterIMF,
     &            parameterIFMR,
     &            burst_age,
     &            ymax,
     &            zy,
     &            zx,
     &            zyimf

          common /param/ fractionOfDB,
     &                   thin_disk_age,
     &                   parameterIMF,
     &                   parameterIFMR,
     &                   burst_age
          common /RSEED/ ISEED1,
     &                   ISEED2
           
C         Releasing a number (y-coordinate) according to a comparison 
C         function cte.
          ymax = mmin ** parameterIMF

          do
              zy = ymax * ran(iseed)
C             Releasing another number (x-coordinate)
              zx = ((mmax - mmin) * ran(iseed)) + mmin
C             Calculating the value of the IMF for m = zx
              zyimf = zx ** parameterIMF
C             Comparing IMF with the function cte       
              if (zy <= zyimf) then
                  mass_from_Salpeter_IMF = zx
                  exit 
              end if
          end do
      end function
