C----------------------------------------------------------------------
C The program monte.f simulates the population of white dwarfs in the 
C solar environment. It distributes randomly n points following a 
C uniform distribution in the galactic plane and an exponential 
C distribution in the perpendicular plane. It adds velocity distribution 
C at each point. From the cooling tables it interpolates brightness and 
C cooling times. Version B: In this version, the SFR function, within 
C each interval in which the time of the galaxy is divided, determines 
C the mass of stars to share. Each mass that is generated, following 
C the IMF, is assigned with a birth time according to the SFR. Finally  
C the maximum volume method is used to calculate the density function of 
C white dwarfs.
C----------------------------------------------------------------------

C     adding FileInfo type which carries all the info about
C     files of cooling and color tables: fort.xx links, numbers of rows
C     and columns, metallicities     
C     NOTE: modules are better than includes
      include 'code/external_types.f'

      program monte
      use external_types 
      implicit real (a-h,m,o-z)
   
C     'external' statement specifies that 'ran' function is no longer
C     intrinsic and must be defined in program
      external ran
      real ran

C     ---  Variables description ---
C     minimumSectorRadius - min radius of the sector of considered stars
C     maximumSectorRadius - max radius of the sector of considered stars
C     angleCoveringSector - angle in degrees, which covers the sector
C     radiusOfSector: radius (kpc) of the sector centered at Sun
C     galacticDiskAge (Gyr)
C     parameterOfSFR (taus): Y=exp(-t/taus)
C     solarGalactocentricDistance: distance from Sun to Galaxy center;
C     parameterIMF (alpha): M^{alpha}
C     Initial-to-Final Mass Relation (IFMR) : 
C         mfinal_new=parameterIFMR*mfinal_old
      integer numberOfStars
      real galacticDiskAge,parameterOfSFR,
     &                 minimumSectorRadius,
     &                 maximumSectorRadius,angleCoveringSector,
     &                 parameterIMF,radiusOfSector,scaleLength,
     &                 areaOfSector,pi
      double precision :: solarGalactocentricDistance
      parameter (numberOfStars=6000000)
      parameter (solarGalactocentricDistance=8.5d0)
      parameter (minimumSectorRadius=8.45)
      parameter (maximumSectorRadius=8.55)
      parameter (angleCoveringSector=0.674)
      parameter (radiusOfSector=0.050)
      parameter (parameterOfSFR=25.0)
      parameter (scaleLength=3.5)
      integer i,j,k,ISEED1,ISEED2,iseed,numberOfStarsInSample
      real randomNumber,fractionOfDB
      real parameterIFMR
C     For terminal:
      integer :: num_args
      character(len = 30) :: arg
      character(len = 30) :: temp_string
      real :: massReductionFactor
      integer :: kinematicModel
      character(len = 100) :: output_filename
      character(len = 6) :: geometry
      real :: cone_height_longitude,
     &                    cone_height_latitude
      real :: min_longitude, max_longitude, min_latitude,
     &                    max_latitude     

      TYPE(FileGroupInfo),DIMENSION(11) :: table

C     NOTE: use of commons is strongly discouraged!
      common /RSEED/ ISEED1,ISEED2
      common /param/ fractionOfDB,galacticDiskAge,parameterIMF,
     &               parameterIFMR,timeOfBurst

C     --- Filling info about groups of files (cooling, color tables) ---    
C ======================================================================
      call fillTable(table)
   

C     ---Reading parameters line from $temporary_files/grid_set_line.in
C ======================================================================
C       read (10,*) fractionOfDB,galacticDiskAge,parameterIMF,
C      &            parameterIFMR,timeOfBurst

C     Overwriting parameters
C     NOTE: this is not good
C       fractionOfDB=0.20 
C       galacticDiskAge=8.9
C       parameterIMF=-2.35
C       parameterIFMR=1.0
C       timeOfBurst=0.6     

C     Terminal reading:
      num_args = iargc()

      if (num_args .eq. 0) then
        fractionOfDB = 0.20 
        galacticDiskAge = 8.9
        parameterIMF = -2.35
        parameterIFMR = 1.0
        timeOfBurst = 0.6 
      else
        do i = 1, num_args
C           call get_command_argument(i, args(i))
          call getarg(i, arg)
          select case(arg)
            case ("-db")
              call getarg(i + 1, temp_string)
              read(temp_string, *) fractionOfDB 
            case ("-g")
              call getarg(i + 1, temp_string)
              read(temp_string, *) galacticDiskAge 
            case ("-mf")
              call getarg(i + 1, temp_string)
              read(temp_string, *) parameterIMF
            case ("-ifr")
              call getarg(i + 1, temp_string)
              read(temp_string, *) parameterIFMR
            case ("-bt")
              call getarg(i + 1, temp_string)
              read(temp_string, *) timeOfBurst
            case ("-mr")
              call getarg(i + 1, temp_string)
              read(temp_string, *) massReductionFactor
            case ("-km")
              call getarg(i + 1, temp_string)
              read(temp_string, *) kinematicModel
            case ("-o")
              call getarg(i + 1, output_filename)
            case ('-geom')
              call getarg(i + 1, geometry)
            case ('-cl')
              call getarg(i + 1, temp_string)
              read(temp_string, *) cone_height_longitude
            case ('-cb')
              call getarg(i + 1, temp_string)
              read(temp_string, *) cone_height_latitude
          end select
        end do
      end if

      write(6,*) '=========================================='
      write(6,*) ' '
      write(6,*) '            Programa monte.f'
      write(6,*) '          by S.Torres, 14.02.11 '
      write(6,*) ' '
      write(6,*) '            Used parameters:'
      write(6,*) 'numberOfStars=    ',numberOfStars
      write(6,*) 'SFR: parameterOfSFR=',parameterOfSFR,'Gyr'
      write(6,*) 'galacticDiskAge=    ',galacticDiskAge,'Gyr'
      write(6,*) 'minimumSectorRadius=',minimumSectorRadius,'kpc' 
      write(6,*) 'maximumSectorRadius=',maximumSectorRadius,'kpc'
      write(6,*) 'radiusOfSector=     ',radiusOfSector,'kpc'
      write(6,*) ' '
      write(6,*) '=========================================='
      write(6,*) ' '
      write(6,*) '          Start of calculations:'
      write(6,*) ' '
      write(6,*) 'Initializing random number generator and reading the s
     &eeds'
      write(6,*) ' '

C     NOTE: this just repeats inputs
C       write (157,157) fractionOfDB,galacticDiskAge,parameterIMF,
C      &                parameterIFMR,timeOfBurst
 157  format(5(f6.3,2x))

C     ---Reading seeds line from $temporary_files/seeds_line.in
C ======================================================================
      iseed=-9
C       read(72,100) iseed1,iseed2
C     This is done for cone geometry and 5000 plates as we don't want 
C     the sequence of RNs to repeat itself(init:805577 133547):
      open(unit=772,file='input_data/seeds_line.in')
      read(unit=772,fmt=*) iseed1, iseed2
      close(unit=772)
C       TODO: add if 
C       iseed1 = iseed1 - 1
C       iseed2 = iseed2 + 1
C       open(unit=772,file='input_data/seeds_line.in',status='replace')
C       write(unit=772,fmt=100) iseed1,iseed2
C       close(unit=772)
C     finished rewriting seeds here
100   format(I6,2x,I6)  
      write(6,*) 'iseed1=',iseed1
      write(6,*) 'iseed2=',iseed2
C     QUESTION: why do we need this part?      
      do i=1,10
        randomNumber=dble(ran(iseed))
        write (6,*) i,randomNumber
      end do  

C     ---  Calculating the area of the sector  ---
C ======================================================================
C     This style ensures maximum precision when assigning a value to PI.
      pi=4.0*atan(1.0)
      areaOfSector=pi*radiusOfSector**2

C     converting cone height parameters from deg to rad      
      if (geometry == 'cone') then
        cone_height_longitude = cone_height_longitude * pi / 180.0
        cone_height_latitude = cone_height_latitude * pi / 180.0
      end if

C     ---  Program itself  ---
C ======================================================================
      write(6,*) '1. Reading the cooling tables (1/9)'

      write(6,*) '   1.1 Tracks of CO DA WD Z=0.001;0.01;0.03;0.06'
C     Calling the function 'incoolda' for 4 metalicities that we have
      call incoolda(table(1))
      call incoolda(table(2))
      call incoolda(table(3))
      call incoolda(table(4))
      
      write(6,*) '   1.2 Tracks of CO non-DA (DB) WD'
C     Calling the function 'incooldb' for 3 metalicities that we have
      call incooldb(table(5))
      call incooldb(table(6))
      call incooldb(table(7))

      write(6,*) '   1.3 Tracks of ONe DA WD'
      call incoolone
      
      write(6,*) '   1.4 Reading the colors table of Rene(DAs) and Berge
     &ron(DBs)'
      call color(table(11))
      call colordb(table(10))

      write (6,*) '   1.5 Reading the tables of CO DA with G variable'      
      call incoolea

      write(6,*) '2. Calling the IMF and SFR (2/9)'
      if (geometry == 'sphere') then
        call gen(iseed,parameterOfSFR,areaOfSector,
     &           numberOfStarsInSample,galacticDiskAge,timeOfBurst,
     &           massReductionFactor,kinematicModel)
      else if (geometry == 'cone') then
        call generate_cone_stars(cone_height_longitude,
     &                           cone_height_latitude,
     &                           numberOfStarsInSample,iseed,
     &                           kinematicModel,galacticDiskAge,
     &                           min_longitude, max_longitude,
     &                           min_latitude, max_latitude)
      else 
        write(6,*) "wrong geometry name, use 'sphere' or 'cone'"
      end if
      write(6,*) "numberOfStarsInSample=", numberOfStarsInSample
      
      write(6,*) '3. Calculating luminosities (3/9)'
      call lumx(iseed,numberOfStarsInSample)      

      ! if cone then we already calculated them
      if (geometry == 'sphere') then
        write(6,*) '4. Calculating polar coordinates (4/9)'
        call polar(iseed,minimumSectorRadius,maximumSectorRadius,
     &         angleCoveringSector,radiusOfSector,
     &         solarGalactocentricDistance,scaleLength)
      end if

      write(6,*) '5. Generating heliocentric velocities (5/9)'
      call velh(iseed,numberOfStarsInSample,geometry)

C     QUESTION: why are we missing the next step?
      goto 7
C     QUESTION: what does this mean?
C     ---  Calculating the trajectories according to/along z-coordinate
      write(6,*) '6. Integrating trajectories (6/9)'
      call traject(galacticDiskAge)

7     write(6,*) '7. Calculating coordinates (7/9)'
      call coor(solarGalactocentricDistance)

      write(6,*) '8. Determinating visual magnitudes (8/9)'
      call magi(fractionOfDB,table) 

C     TODO: give a better description to this step
C     NOTE: This will be in a separate processing module
C       write(6,*) '9. Working with obtained sample (9/9)'
C       call volum_40pc
      call printForProcessing(output_filename, geometry, iseed,
     &         min_longitude, max_longitude, min_latitude, max_latitude,
     &         solarGalactocentricDistance,cone_height_longitude,
     &         cone_height_latitude)


      write (6,*) 'End'
 

      stop
      end
C***********************************************************************
      include 'code/star_generation/generator.f'

      include 'code/star_generation/cone.f'
     
      include 'code/cooling/DA/DA_cooling.f'

      include 'code/cooling/unknown/incoolea.f'

      include 'code/colors/DA/byRenedo.f'

      include 'code/luminosities/luminosities.f'

      include 'code/coordinates/polar.f'

      include 'code/velocities/velocities.f'      

      include 'code/math/random_number_generators.f'

      include 'code/trajectories/trajectories.f'
      
      include 'code/coordinates/coor.f'

      include 'code/magnitudes/magi.f'      

      include 'code/magnitudes/DA/interlumda.f'

      include 'code/DA_DB_fraction/dbd_fid.f'

      include 'code/magnitudes/DA/intermag.f'

      include 'code/magnitudes/DB/interlumdb.f'

      include 'code/colors/DB/byBergeron.f'

      include 'code/cooling/ONe/incoolone.f'

      include 'code/cooling/DB/incooldb.f'

      include 'code/magnitudes/ONe/interlumONe.f'

      include 'code/magnitudes/interp.f'

      include 'code/colors/chanco.f'

      include 'code/samples/volum_40pc.f'
      
      include 'code/velocities/vrado.f'

      include 'code/math/toSort.f'

      include 'code/tables_linking.f'

      subroutine printForProcessing(output_filename, geometry, iseed,
     &         min_longitude, max_longitude, min_latitude, max_latitude,
     &         solarGalactocentricDistance,cone_height_longitude,
     &         cone_height_latitude)
      implicit none
      external ran
      real ran
      character(len = 100), intent(in) :: output_filename
      character(len = 6) :: geometry
      integer i,j, iseed
      logical eleminationFlag
      integer numberOfStars,numberOfWDs
      integer eleminatedByParallax,eleminatedByDeclination,
     &        eleminatedByProperMotion,eleminatedByApparentMagn,
     &        eleminatedByReducedPropM
      real parameterIFMR
      real minimumProperMotion,declinationLimit,
     &                 minimumParallax
      real mbolmin,mbolinc,mbolmax      
      real errinfa,errsupa,mbol
      real fnora,fnor,pi,vvv,x,xx,xya
      real min_longitude, max_longitude, min_latitude,
     &                 max_latitude
      double precision :: solarGalactocentricDistance
      real prev_min_longitude, prev_max_longitude, 
     &                 prev_min_latitude, prev_max_latitude
      
      parameter (numberOfStars=6000000)
C     (Only northern hemisphere)
      parameter (declinationLimit=0.0)
C     Minimum parallax below which we discard results (0.025<=>40 pc)
      parameter (minimumParallax=0.025)
C     Binning of Luminosity Function 
      parameter (mbolmin=5.75,mbolmax=20.75,mbolinc=0.5)
C     Minimum proper motion
      parameter (minimumProperMotion=0.04)
C     Parameters of mass histograms

      double precision :: properMotion(numberOfStars),
     &                    rightAscension(numberOfStars),
     &                    declination(numberOfStars)
      real luminosityOfWD(numberOfStars),
     &                 massOfWD(numberOfStars),
     &                 metallicityOfWD(numberOfStars),
     &                 effTempOfWD(numberOfStars)
      real flagOfWD(numberOfStars)
C     rgac - galactocentric distance to WD TODO: give a better name
      real rgac(numberOfStars)
      double precision :: lgac(numberOfStars),
     &                    bgac(numberOfStars)
      real :: coolingTime(numberOfStars)
C     NOTE: this 70 comes from nowhere      
      integer numberOfWDsInBin(70),numberOfBins
      double precision :: coordinate_R(numberOfStars),
     &                    coordinate_Theta(numberOfStars),
     &                    coordinate_Zcylindr(numberOfStars)
      real parallax(numberOfStars)
      real tangenVelo(numberOfStars)
      real mpl(numberOfStars),mpb(numberOfStars),
     &                 vr(numberOfStars)
      real errora(70),ndfa(70)
C     ugriz-color system and V-band from Johnson system
      real v(numberOfStars)
      real massInBin(70)
      integer :: typeOfWD(numberOfStars)
C     values of LF in each bin. o-observational
      real xfl(19),xflo(19),xflcut(3),xflocut(3)
      real xflhot(11),xflohot(11)
C     bins for max-region: for synthetic and observational samples
      real xflMaxRegion(6), xfloMaxRegion(6)
C     number of WDs in bin of mass histogram
      real nbinmass(26)
C     WDs velocities. QUESTION: relative to what?
      real uu(numberOfStars), vv(numberOfStars), 
     &                 ww(numberOfStars)
C     sum of WDs velocities in specific bin, _u/_v/_w - components
      real sumOfWDVelocitiesInBin_u(70),
     &                 sumOfWDVelocitiesInBin_v(70),
     &                 sumOfWDVelocitiesInBin_w(70)
C     average velocity for WDs in specific bin      
      real averageWDVelocityInBin_u(70),
     &                 averageWDVelocityInBin_v(70),
     &                 averageWDVelocityInBin_w(70)  
C     this is used to calculate sigma (SD)
      real sumOfSquareDifferences_u,
     &                 sumOfSquareDifferences_v,
     &                 sumOfSquareDifferences_w
C     SD for velocities in each bin
      real standardDeviation_u(70),standardDeviation_v(70),
     &                 standardDeviation_w(70) 
C     2D-array of velocities (nº of bin; newly assigned to WD nº in bin)
C     needed to calculate Standart Deviation (SD) for velocities in each 
C     bin
C     TODO: make dynamic array or linked list
      real arrayOfVelocitiesForSD_u(25,50000)
      real arrayOfVelocitiesForSD_v(25,50000)
      real arrayOfVelocitiesForSD_w(25,50000)
C     2D-array of bolometric magnitudes for each WD; indexes are the 
C     same as for arrayOfVelocitiesForSD_u/v/w. (For cloud)
      real arrayOfMagnitudes(25,50000)
      real x_coordinate,y_coordinate,z_coordinate,
     & star_longitude, star_latitude
      integer disk_belonging(numberOfStars)
      logical :: overlapping_cone
      integer :: getNumberOfLines, overlapping_cones_count, lines_count
      real, allocatable :: overlap_min_longitudes(:)
      real, allocatable :: overlap_max_longitudes(:)
      real, allocatable :: overlap_min_latitudes(:)
      real, allocatable :: overlap_max_latitudes(:)
      real latitude, longitude,zzx,cone_height_longitude
      double precision :: ros
      real cone_height_latitude
C       TODO: uncomment when ubvri diagrams plotted
      real :: ugriz_ug(numberOfStars),
     &                    ugriz_gr(numberOfStars),
     &                    ugriz_ri(numberOfStars),
     &                    ugriz_iz(numberOfStars),
     &                    ugriz_g_apparent(numberOfStars)
      logical :: with_overlapping_checking = .true.
      integer processed_cones_count, overlappings, stars_counter,
     & eliminations_counter
      real DELTA_LATITUDE,delta_longitude,prev_cone_latitude
      real prev_cone_longitude
      logical longitude_overlapping, latitude_overlapping,
     &  star_in_intersection
      real :: UB(numberOfStars), BV(numberOfStars), 
     &                    VRR(numberOfStars), RI(numberOfStars)
C      &                    ugriz_g_apparent(numberOfStars)
      common /enanas/ luminosityOfWD,massOfWD,metallicityOfWD,
     &                effTempOfWD
      common /index/ flagOfWD,numberOfWDs,disk_belonging      
      common /mad/ properMotion,rightAscension,declination
      common /mopro/ mpb,mpl,vr
      common /paral/ rgac
      common /lb/ lgac,bgac
      common /coorcil/ coordinate_R,coordinate_Theta,coordinate_Zcylindr
      common /cool/ coolingTime
C       TODO: uncomment when ubvri diagrams plotted
      common /photo/ ugriz_ug, ugriz_gr, ugriz_ri, ugriz_iz, 
     &               ugriz_g_apparent
      common /indexdb/ typeOfWD
      common /johnson/ V
      common /vel/ uu,vv,ww
      common /ubvri/ UB, BV, VRR, RI
C      &, ugriz_g_apparent
      pi = 4.0 * atan(1.0)
      DELTA_LATITUDE = 2.64 * pi / 180.0
      overlappings = 0
      stars_counter = 0
      eliminations_counter = 0

      if (geometry == 'sphere') then

          open(421, file = output_filename)

          write(421, *) 'luminosity ',
     &                  'proper_motion ',
     &                  'proper_motion_component_b ',
     &                  'proper_motion_component_l ',
     &                  'proper_motion_component_vr ',
     &                  'right_ascension ',
     &                  'declination ',
     &                  'galactic_distance ',
     &                  'galactic_latitude ',
     &                  'galactic_longitude ',
     &                  'ugriz_g_apparent ',
     &                  'ugriz_ug ',    
     &                  'ugriz_gr ',    
     &                  'ugriz_ri ',
     &                  'ugriz_iz ',
     &                  'v_photometry ',
     &                  'velocity_u ',
     &                  'velocity_v ',
     &                  'velocity_w ',
     &                  'spectral_type'

          do i = 1, numberOfWDs
              write(421, *) luminosityOfWD(i),
     &                      properMotion(i),
     &                      mpb(i),
     &                      mpl(i),
     &                      vr(i),
     &                      rightAscension(i),
     &                      declination(i),
     &                      rgac(i),
     &                      bgac(i),
     &                      lgac(i),
     &                      ugriz_g_apparent(i),
C                           TODO: uncomment
     &                      ugriz_ug(i),
     &                      ugriz_gr(i),
     &                      ugriz_ri(i),
     &                      ugriz_iz(i),
     &                      V(i),
     &                      uu(i),
     &                      vv(i),
     &                      ww(i),
     &                      typeOfWD(i)
          end do
      else if (geometry == 'cone') then

         open(421,file='cone_stars_catalog.csv',access='append')
C          TODO: uncomment this when writing to Cassandra works well        
C          open(421,file=output_filename)

C          TODO: uncomment this when Cassandra is tuned  
C          write(421, *) 'velocity_u ',
C      &                 'velocity_v ',
C      &                 'velocity_w ',
C      &                 'galactic_distance ',
C      &                 'galactic_longitude ',
C      &                 'galactic_latitude ',
C      &                 'ugriz_g_apparent ',
C      &                 'ugriz_ug ',
C      &                 'ugriz_gr ',
C      &                 'ugriz_ri ',
C      &                 'ugriz_iz ',
C      &                 'disk_belonging ',
C      &                 'spectral_type'

C        TODO:this should be in processing. move it there when Cas is ok
         if (with_overlapping_checking .eqv. .true.) then
           write(6,*) 'Overlapping checking is on'
           delta_longitude = DELTA_LATITUDE / cos(cone_height_latitude)
           min_longitude = cone_height_longitude - delta_longitude/2.0
           max_longitude = cone_height_longitude + delta_longitude/2.0
           min_latitude = cone_height_latitude - DELTA_LATITUDE/2.0
           max_latitude = cone_height_latitude + DELTA_LATITUDE/2.0
           write(6,*) 'Current angles:'
           write(6,*) 'L min: ', min_longitude * 180.0 / pi
           write(6,*) 'L max: ', max_longitude * 180.0 / pi 
           write(6,*) 'B min: ', min_latitude * 180.0 / pi
           write(6,*) 'B max: ', max_latitude * 180.0 / pi

           processed_cones_count 
     &         = getNumberOfLines('processed_cones.txt')
           write(6,*) 'Processed cones: ', processed_cones_count 
           allocate(overlap_min_longitudes(processed_cones_count))
           allocate(overlap_max_longitudes(processed_cones_count))
           allocate(overlap_min_latitudes(processed_cones_count))
           allocate(overlap_max_latitudes(processed_cones_count))

           open(765, file='processed_cones.txt')
           do i = 1, processed_cones_count
            read(765,*) prev_cone_longitude, prev_cone_latitude
            prev_cone_longitude = prev_cone_longitude * pi / 180.0
            prev_cone_latitude = prev_cone_latitude * pi / 180.0

            delta_longitude=DELTA_LATITUDE/cos(prev_cone_latitude)
            prev_min_longitude=prev_cone_longitude-delta_longitude/2.0
            prev_max_longitude=prev_cone_longitude+delta_longitude/2.0
            prev_min_latitude=prev_cone_latitude - DELTA_LATITUDE/2.0
            prev_max_latitude=prev_cone_latitude + DELTA_LATITUDE/2.0

            if (min_longitude < (prev_max_longitude - 2.0*pi) 
     &          .and. min_longitude < 0.0) then
                   prev_min_longitude = prev_min_longitude - 2.0*pi
                   prev_max_longitude = prev_max_longitude - 2.0*pi
            end if
            if (prev_min_longitude < (max_longitude - 2.0*pi)
     &          .and. prev_min_longitude < 0.0) then
                   min_longitude = min_longitude - 2.0*pi
                   max_longitude = max_longitude - 2.0*pi
            end if

            longitude_overlapping = .false.
            latitude_overlapping = .false.
C            write(6,*) 'L overlap before check ', longitude_overlapping
C            write(6,*) 'B overlap before check ', latitude_overlapping
            if (((prev_min_longitude <= min_longitude) .and.
     &           (min_longitude <= prev_max_longitude)) 
     &          .or. ((prev_min_longitude <= max_longitude) .and.
     &                (max_longitude <= prev_max_longitude))
     &          .or. ((min_longitude <= prev_min_longitude) .and.
     &                (prev_min_longitude<= max_longitude))
     &          .or. ((min_longitude <= prev_max_longitude) .and.
     &                (prev_max_longitude<= max_longitude))) then
                    longitude_overlapping = .true.
C                     write(6,*) 'LONGITUDE OVERLAPPING'
            end if
            if (((prev_min_latitude <= min_latitude) .and.
     &            (min_latitude <= prev_max_latitude))
     &           .or. ((prev_min_latitude <= max_latitude) .and.
     &            (max_latitude <= prev_max_latitude))
     &           .or. ((min_latitude <= prev_min_latitude) .and.
     &            (prev_min_latitude <= max_latitude))
     &           .or. ((min_latitude <= prev_max_latitude) .and.
     &            (prev_max_latitude <= max_latitude))) then
                    latitude_overlapping = .true.
C                     write(6,*) 'LATITUDE OVERLAPPING'
            end if

            if ((longitude_overlapping .eqv. .true.) 
     &          .and. (latitude_overlapping .eqv. .true.)) then

C              write(6,*) 'L overlap after check', longitude_overlapping
C              write(6,*) 'B overlap after check', latitude_overlapping
              overlappings = overlappings + 1

              overlap_min_longitudes(overlappings) = prev_min_longitude
              overlap_max_longitudes(overlappings) = prev_max_longitude
              overlap_min_latitudes(overlappings) = prev_min_latitude
              overlap_max_latitudes(overlappings) = prev_max_latitude

              write(6,*) 'Overlapping with line: ', i
C              write(6,*) 'Prev Lmin: ', prev_min_longitude * 180.0 / pi
C              write(6,*) 'Prev Lmax: ', prev_max_longitude * 180.0 / pi
C              write(6,*) 'Prev Bmin: ', prev_min_latitude * 180.0 / pi
C              write(6,*) 'Prev Bmax: ', prev_max_latitude * 180.0 / pi
            end if
           end do

           write(6,*) 'Overlappings count:', overlappings

           do i = 1, numberOfWDs
            ros = solarGalactocentricDistance 
     &            * solarGalactocentricDistance + coordinate_R(i) 
     &                                            * coordinate_R(i)
     &            - 2.d0 * coordinate_R(i) 
     &              * solarGalactocentricDistance 
     &              * dcos(coordinate_Theta(i))
            ros = dsqrt(ros)
C           TODO: figurre out what to do with ill-conditioned cases            
            if ((solarGalactocentricDistance ** 2 + ros**2
     &           - coordinate_R(i) ** 2)
     &          /(2.0 * solarGalactocentricDistance * ros) > 1.0) then
              longitude = 0.0
            else
            longitude = dacos((solarGalactocentricDistance ** 2 + ros**2
     &                         - coordinate_R(i) ** 2)
     &                  /(2.0d0 * solarGalactocentricDistance * ros))
            end if
            zzx=coordinate_Zcylindr(i)/ros
            latitude = atan(zzx)
C           TODO: i am confused about how all this works now.
C           the reason for these conditions is that for angles > pi
C           longitude is simmetrically reflected on the top half-plane
            if (coordinate_Theta(i) < 0
     &          .and. cone_height_longitude > 3.0 * pi / 2) then
                longitude = 2.0 * pi - longitude
            end if
            if (coordinate_Theta(i) > 0 
     &          .and. cone_height_longitude > 3.0 * pi / 2) then
                longitude = 2.0 * pi + longitude
            end if
            if (coordinate_Theta(i) < 0 
     &          .and. cone_height_longitude < pi / 2) then
                longitude = -longitude
            end if
            if (coordinate_Theta(i) < 0 
     &          .and. cone_height_longitude < 3.0 * pi / 2
     &          .and. cone_height_longitude > pi / 2) then
                longitude = 2.0 * pi - longitude
            end if
C           if cone crosses 2pi, move it -2pi
            if (max_longitude > 2 * pi) then
                longitude = longitude - 2 * pi
            end if

            star_in_intersection = .false.
            if (overlappings > 0) then
              do j = 1, overlappings
                if (((overlap_min_longitudes(j) < longitude) .and.
     &                      (longitude < overlap_max_longitudes(j)))
     &                      .and. ((overlap_min_latitudes(j) 
     &                           < latitude) .and. 
     &                      (latitude < overlap_max_latitudes(j)))) then
                  star_in_intersection = .true.
                end if
              end do
            end if
            if (star_in_intersection .eqv. .false.) then
              stars_counter = stars_counter + 1
C               TODO: uncomment this after ubvri diag-s plotted
C               write(421,"(11(es17.8e3,x),i1,x,f3.1)") 
                write(421,"(es17.8e3,x,i1,x,i1,x,4(es17.8e3,x))") 
C      &                                              uu(i),
C      &                                              vv(i),
C      &                                              ww(i),
C      &                                              rgac(i),
C      &                                              longitude,
C      &                                              latitude,
     &                                              ugriz_g_apparent(i),
C      &                                              ugriz_ug(i),
C      &                                              ugriz_gr(i),
C      &                                              ugriz_ri(i),
C      &                                              ugriz_iz(i),
     &                                              disk_belonging(i),
     &                                              typeOfWD(i)
C                                            TODO: delete this after..
C      &                                       ,UB(i),BV(i),VRR(i),RI(i)
            else
                eliminations_counter = eliminations_counter + 1
            end if
           end do
           write(6, *) 'Recorded stars: ', stars_counter
           write(6, *) 'Eliminated stars: ', eliminations_counter
         else

         do i = 1, numberOfWDs
C              x_coordinate=8.5-coordinate_R(i) * cos(coordinate_Theta(i))
C              y_coordinate = coordinate_R(i) * sin(coordinate_Theta(i))
C              z_coordinate = coordinate_Zcylindr(i)

C              star_longitude = atan(y_coordinate / x_coordinate)
C              star_latitude =atan(z_coordinate/sqrt(x_coordinate ** 2
C    &                                             + y_coordinate ** 2))
C           TODO: i took this from coor.f
            ros = solarGalactocentricDistance 
     &            * solarGalactocentricDistance 
     &            + coordinate_R(i) * coordinate_R(i)
     &            -2.0 * coordinate_R(i) * solarGalactocentricDistance
     &             * cos(coordinate_Theta(i))
            ros = sqrt(ros)
C           TODO: figurre out what to do with ill-conditioned cases            
            if ((solarGalactocentricDistance ** 2 + ros**2
     &           - coordinate_R(i) ** 2)
     &          /(2.d0 * solarGalactocentricDistance * ros) > 1.0) then
              longitude = 0.0
            else
            longitude = dacos((solarGalactocentricDistance ** 2 + ros**2
     &                         - coordinate_R(i) ** 2)
     &                      /(2.d0 * solarGalactocentricDistance * ros))
            end if
            zzx=coordinate_Zcylindr(i)/ros
            latitude = atan(zzx)
C           TODO: i am confused about how all this works now.
C           the reason for these conditions is that for angles > pi
C           longitude is simmetrically reflected on the top half-plane
            if (coordinate_Theta(i) < 0
     &          .and. cone_height_longitude > 3.0 * pi / 2) then
                longitude = 2.0 * pi - longitude
            end if
            if (coordinate_Theta(i) > 0 
     &          .and. cone_height_longitude > 3.0 * pi / 2) then
                longitude = 2.0 * pi + longitude
            end if
            if (coordinate_Theta(i) < 0 
     &          .and. cone_height_longitude < pi / 2) then
                longitude = -longitude
            end if
            if (coordinate_Theta(i) < 0 
     &          .and. cone_height_longitude < 3.0 * pi / 2
     &          .and. cone_height_longitude > pi / 2) then
                longitude = 2.0 * pi - longitude
            end if
C           if cone crosses 2pi, move it -2pi
            if (max_longitude > 2 * pi) then
                longitude = longitude - 2 * pi
            end if
            write(421,"(11(es17.8e3,x),i1,x,f3.1)") uu(i),
     &                                              vv(i),
     &                                              ww(i),
     &                                              rgac(i),
     &                                              longitude,
     &                                              latitude,
     &                                              ugriz_g_apparent(i),
C      TODO: uncomment
C      &                                              ugriz_ug(i),
C      &                                              ugriz_gr(i),
C      &                                              ugriz_ri(i),
C      &                                              ugriz_iz(i),
     &                                              disk_belonging(i),
     &                                              typeOfWD(i)
         end do
         end if
      end if
      end subroutine


      function getNumberOfLines(filePath) result(n)
        character(len = *), intent(in) :: filePath
        integer :: n, ioStatus
        n = 0
        open(532, file = filePath)
        do
            read(532, *, iostat = ioStatus)
            if(is_iostat_end(ioStatus)) exit
            n = n + 1
        end do
        close(532)
      end function
