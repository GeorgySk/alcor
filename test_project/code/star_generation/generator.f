      subroutine gen(iseed,
     &               parameterOfSFR,
     &               radius,
     &               numberOfStarsInSample,
     &               galacticDiskAge,
     &               timeOfBurst,
     &               massReductionFactor,
     &               thick_disk_stars_fraction,
     &               halo_stars_fraction)
C     Divides the SFR in intervals of time. In each interval, the total 
C     mass of stars is distributed. The mass of each star follows the 
C     distribution law given by Initial Mass Function (IMF). The birth 
C     time is also calculated, from which the scale height and finally 
C     cylindrical z-coordinate are determined.       
C           TODO: implement more modern RNG
C           iseed: seed for random generator
C           timeOfBurst: how many Gyr ago burst of star gener-n happened
      implicit none

C     Overloading intrinsic 'ran' function by our own RNG
      external ran
      real ran

C     TODO: this is numberOfStars, take it up as const
      integer, parameter :: MAX_STARS_COUNT = 6000000,
     &                      BINS_COUNT = 5000
      real, parameter :: Z_DISTRIBUTION_ZO = 1.0,
     &                   THIN_DISK_SCALE_HEIGHT_KPC = 0.250,
C                        TODO: this const is in cone generator, take up
     &                   THICK_DISK_SCALE_HEIGHT_KPC = 0.900,
C                        TODO: these are params of thick disk. exactly?
     &                   TTDISK = 12.0,
     &                   TMDISK = 10.0,
     &                   TAU = 2.0,
     &                   STARS_BIRTH_START_TIME = 0.0,
     &                   BURST_FORMATION_FACTOR = 5.0,
     &                   HALO_STARS_BIRTH_END_TIME = 1.0,
     &                   HALO_DISTANCE_PARAM = 5.0,
     &                   PI = 4.0*atan(1.0)
      integer ::  iseed,
     &            numberOfStarsInSample, 
     &            bin_index, 
     &            stars_count, 
     &            in, 
     &            numberOfWDs,
     &            flagOfWD(MAX_STARS_COUNT)
      real :: zz,                          tmax,
     &        ft,                          fz,
     &        timeOfBurst,                 burst_start_time, 
     &        current_bin_initial_time,    xx,
     &        ttry,                        thick_disk_stars_fraction,
     &        parameterOfSFR,              areaOfSector, 
     &        galacticDiskAge,             normalization_const,
     &        time_increment,              massReductionFactor,
     &        psi,                         star_mass,
     &        total_generated_mass_in_bin, mrep,
     &        halo_stars_fraction,         burst_mrep,
     &        normal_mrep,                 random_value,
     &        radius
      real :: mass_from_Salpeter_IMF, 
     &        get_normalization_const
      double precision :: coordinate_Theta(MAX_STARS_COUNT),
     &                    coordinate_R(MAX_STARS_COUNT),
     &                    coordinate_Zcylindr(MAX_STARS_COUNT)
C     TODO: find out what is m. massInMainSequence? 
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

      areaOfSector = PI * radius ** 2

C     Calculating the mass, time of birth and z-coordinate of every star
C     Calculating the normalization constant of the SFR
      normalization_const = get_normalization_const(parameterOfSFR, 
     &                                              galacticDiskAge)
      time_increment = (galacticDiskAge - STARS_BIRTH_START_TIME) 
     &                 / BINS_COUNT
C     TODO: find out the meaning
      psi = normalization_const * time_increment
C     TODO: find out the meaning of 1.0e6 and mrep
      mrep = psi * areaOfSector * 1.0e6
      mrep = mrep * massReductionFactor
      normal_mrep = mrep
      burst_mrep = mrep * BURST_FORMATION_FACTOR

      stars_count = 0

      write(6,*) '      Factor of mass reduction=', massReductionFactor

C     Calculating the mass to be distributed at each interval       
      do bin_index = 1, BINS_COUNT
          total_generated_mass_in_bin = 0.0

C         Recent burst
          burst_start_time = galacticDiskAge - timeOfBurst
          current_bin_initial_time = STARS_BIRTH_START_TIME 
     &                               + float(bin_index - 1) 
     &                                 * time_increment
    
          if (current_bin_initial_time >= burst_start_time 
     &            .and. current_bin_initial_time < galacticDiskAge) then
              mrep = burst_mrep
          else
              mrep = normal_mrep
          endif

          do
C             Calling to the IMF
              star_mass = mass_from_Salpeter_IMF(iseed)

C             We already have the mass
              stars_count = stars_count + 1

              if (stars_count == MAX_STARS_COUNT) then 
                  write(6,*) '     ***  Dimension exceeded   ***'
                  write(6,*) '***  Increase the reduction factor   ***'
                  stop
              end if

              m(stars_count) = star_mass       
         
C             Birth time from SFR constant 
C             disk_belonging = 1 (thin disk), = 2 (thick disk)
              random_value = ran(iseed)
              if (random_value <= thick_disk_stars_fraction 
     &                             + halo_stars_fraction
     &                .and. random_value > halo_stars_fraction) then
                  disk_belonging(stars_count) = 2
C                 TODO: find out what is going on here
                  tmax = TMDISK * exp(-TMDISK / TAU)
                  do
                      ttry = TTDISK * ran(iseed)
                      ft = ttry * exp(-ttry / TAU)
                      fz = tmax * ran(iseed)
                      if (fz <= ft) then
                          starBirthTime(stars_count) = ttry
                          exit
                      end if
                  end do
              else if (random_value <= halo_stars_fraction) then
                  disk_belonging(stars_count) = 3
                  starBirthTime(stars_count) = ran(iseed) 
     &                                       * HALO_STARS_BIRTH_END_TIME
              else
                  disk_belonging(stars_count) = 1
                  starBirthTime(stars_count) = STARS_BIRTH_START_TIME 
     &                               + float(bin_index - 1) 
     &                                 * time_increment 
     &                               + time_increment * ran(iseed)
                  total_generated_mass_in_bin = 
     &                total_generated_mass_in_bin + star_mass
              end if

C             Calculating z
              if (disk_belonging(stars_count) == 1) then
                  scale_height = THIN_DISK_SCALE_HEIGHT_KPC
              else if (disk_belonging(stars_count) == 2) then
                  scale_height = THICK_DISK_SCALE_HEIGHT_KPC
              end if

              if (disk_belonging(stars_count) /= 3) then
C               TODO: find out what is going on here
                do
                  xx = Z_DISTRIBUTION_ZO * scale_height * ran(iseed)
                  if (xx /= 0.0) then
                    exit
                  end if
                end do
                zz = scale_height * LOG(Z_DISTRIBUTION_ZO 
     &                                  * scale_height / xx)   
C               z-contstant       zz = 0.240 * ran(iseed)   
                in = int(2.0 * ran(iseed))
                coordinate_Zcylindr(stars_count) = zz * dfloat(1 - 2*in)
              end if

C             Checking if we have generated enough mass
              if (total_generated_mass_in_bin >= mrep) then
                  m(stars_count) = m(stars_count) 
     &                             - (total_generated_mass_in_bin 
     &                                - mrep)
                  total_generated_mass_in_bin = mrep
                  exit
              end if
          end do
      end do
        
      numberOfStarsInSample = stars_count
      write(6,*) '     Number of stars in sample=',numberOfStarsInSample
      
      end subroutine


      function get_normalization_const(parameterOfSFR, 
     &                                 galacticDiskAge)
C         Calculating the normalization constant of the SFR
          implicit none

C         TODO: find out what this sigma is
          real, parameter :: sigma = 51.0
          real :: parameterOfSFR,
     &            galacticDiskAge,
     &            get_normalization_const

          get_normalization_const = sigma 
     &                              / (parameterOfSFR 
     &                                 * (1.0 - exp(-galacticDiskAge 
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
     &            galacticDiskAge,
     &            parameterIMF,
     &            parameterIFMR,
     &            timeOfBurst,
     &            ymax,
     &            zy,
     &            zx,
     &            zyimf

          common /param/ fractionOfDB,
     &                   galacticDiskAge,
     &                   parameterIMF,
     &                   parameterIFMR,
     &                   timeOfBurst
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
