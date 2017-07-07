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
      implicit double precision (a-h,m,o-z)
   
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
      double precision galacticDiskAge,parameterOfSFR,
     &                 solarGalactocentricDistance,minimumSectorRadius,
     &                 maximumSectorRadius,angleCoveringSector,
     &                 parameterIMF,radiusOfSector,scaleLength,
     &                 areaOfSector,pi   
      parameter (numberOfStars=6000000)
      parameter (solarGalactocentricDistance=8.5)
      parameter (minimumSectorRadius=8.45)
      parameter (maximumSectorRadius=8.55)
      parameter (angleCoveringSector=0.674)
      parameter (radiusOfSector=0.050)
      parameter (parameterOfSFR=25.0)
      parameter (scaleLength=3.5)
      integer i,j,k,ISEED1,ISEED2,iseed,numberOfStarsInSample
      double precision randomNumber,fractionOfDB
      double precision parameterIFMR
C     For terminal:
      integer :: num_args
      character(len = 30) :: arg
      character(len = 30) :: temp_string
      double precision :: massReductionFactor
      integer :: kinematicModel
      character(len = 100) :: output_filename
      character(len = 6) :: geometry
      double precision :: cone_height_longitude,
     &                    cone_height_latitude
      double precision :: min_longitude, max_longitude, min_latitude,
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
      iseed1 = iseed1 - 1
      iseed2 = iseed2 + 1
      open(unit=772,file='input_data/seeds_line.in',status='replace')
      write(unit=772,fmt=100) iseed1,iseed2
      close(unit=772)
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
     &         solarGalactocentricDistance,cone_height_longitude)


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
     &         solarGalactocentricDistance,cone_height_longitude)
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
      double precision parameterIFMR
      double precision minimumProperMotion,declinationLimit,
     &                 minimumParallax
      double precision mbolmin,mbolinc,mbolmax      
      double precision errinfa,errsupa,mbol
      double precision fnora,fnor,pi,vvv,x,xx,xya
      double precision min_longitude, max_longitude, min_latitude,
     &                 max_latitude, solarGalactocentricDistance
      double precision prev_min_longitude, prev_max_longitude, 
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

      double precision properMotion(numberOfStars),
     &                 rightAscension(numberOfStars),
     &                 declination(numberOfStars)
      double precision luminosityOfWD(numberOfStars),
     &                 massOfWD(numberOfStars),
     &                 metallicityOfWD(numberOfStars),
     &                 effTempOfWD(numberOfStars)
      double precision flagOfWD(numberOfStars)
C     rgac - galactocentric distance to WD TODO: give a better name
      double precision rgac(numberOfStars)
      double precision lgac(numberOfStars)
      double precision bgac(numberOfStars)
      double precision coolingTime(numberOfStars)
C     NOTE: this 70 comes from nowhere      
      integer numberOfWDsInBin(70),numberOfBins
      double precision coordinate_R(numberOfStars),
     &                 coordinate_Theta(numberOfStars),
     &                 coordinate_Zcylindr(numberOfStars)
      double precision parallax(numberOfStars)
      double precision tangenVelo(numberOfStars)
      double precision mpl(numberOfStars),mpb(numberOfStars),
     &                 vr(numberOfStars)
      double precision errora(70),ndfa(70)
C     ugriz-color system and V-band from Johnson system
      double precision go(numberOfStars),gr(numberOfStars),
     &                 v(numberOfStars)
      double precision gi(numberOfStars),ur(numberOfStars),
     &                 rz(numberOfStars)
      double precision massInBin(70)
      integer typeOfWD(numberOfStars)
C     values of LF in each bin. o-observational
      double precision xfl(19),xflo(19),xflcut(3),xflocut(3)
      double precision xflhot(11),xflohot(11)
C     bins for max-region: for synthetic and observational samples
      double precision xflMaxRegion(6), xfloMaxRegion(6)
C     number of WDs in bin of mass histogram
      double precision nbinmass(26)
C     WDs velocities. QUESTION: relative to what?
      double precision uu(numberOfStars), vv(numberOfStars), 
     &                 ww(numberOfStars)
C     sum of WDs velocities in specific bin, _u/_v/_w - components
      double precision sumOfWDVelocitiesInBin_u(70),
     &                 sumOfWDVelocitiesInBin_v(70),
     &                 sumOfWDVelocitiesInBin_w(70)
C     average velocity for WDs in specific bin      
      double precision averageWDVelocityInBin_u(70),
     &                 averageWDVelocityInBin_v(70),
     &                 averageWDVelocityInBin_w(70)  
C     this is used to calculate sigma (SD)
      double precision sumOfSquareDifferences_u,
     &                 sumOfSquareDifferences_v,
     &                 sumOfSquareDifferences_w
C     SD for velocities in each bin
      double precision standardDeviation_u(70),standardDeviation_v(70),
     &                 standardDeviation_w(70) 
C     2D-array of velocities (nº of bin; newly assigned to WD nº in bin)
C     needed to calculate Standart Deviation (SD) for velocities in each 
C     bin
C     TODO: make dynamic array or linked list
      double precision arrayOfVelocitiesForSD_u(25,50000)
      double precision arrayOfVelocitiesForSD_v(25,50000)
      double precision arrayOfVelocitiesForSD_w(25,50000)
C     2D-array of bolometric magnitudes for each WD; indexes are the 
C     same as for arrayOfVelocitiesForSD_u/v/w. (For cloud)
      double precision arrayOfMagnitudes(25,50000)
      double precision x_coordinate,y_coordinate,z_coordinate,
     & star_longitude, star_latitude
      integer disk_belonging(numberOfStars)
      logical :: overlapping_cone
      integer :: getNumberOfLines, overlapping_cones_count, lines_count
      double precision, allocatable :: overlap_min_longitudes(:)
      double precision, allocatable :: overlap_max_longitudes(:)
      double precision, allocatable :: overlap_min_latitudes(:)
      double precision, allocatable :: overlap_max_latitudes(:)
      double precision latitude, longitude,ros,zzx,cone_height_longitude
      common /enanas/ luminosityOfWD,massOfWD,metallicityOfWD,
     &                effTempOfWD
      common /index/ flagOfWD,numberOfWDs,disk_belonging      
      common /mad/ properMotion,rightAscension,declination
      common /mopro/ mpb,mpl,vr
      common /paral/ rgac
      common /lb/ lgac,bgac
      common /coorcil/ coordinate_R,coordinate_Theta,coordinate_Zcylindr
      common /cool/ coolingTime
      common /photo/ go,gr,gi,ur,rz
      common /indexdb/ typeOfWD
      common /johnson/ V
      common /vel/ uu,vv,ww
      pi = 4.0 * atan(1.0)

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
     &                  'go_photometry ',
     &                  'gr_photometry ',
     &                  'rz_photometry ',
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
     &                      go(i),
     &                      gr(i),
     &                      rz(i),
     &                      V(i),
     &                      uu(i),
     &                      vv(i),
     &                      ww(i),
     &                      typeOfWD(i)
          end do
      else if (geometry == 'cone') then

         open(421,file=output_filename)

         write(421, *) 'velocity_u ',
     &                 'velocity_v ',
     &                 'velocity_w ',
     &                 'galactic_distance ',
     &                 'galactic_longitude ',
     &                 'galactic_latitude ',
     &                 'disk_belonging'

         do i = 1, numberOfWDs
C              x_coordinate=8.5-coordinate_R(i) * cos(coordinate_Theta(i))
C              y_coordinate = coordinate_R(i) * sin(coordinate_Theta(i))
C              z_coordinate = coordinate_Zcylindr(i)

C              star_longitude = atan(y_coordinate / x_coordinate)
C              star_latitude =atan(z_coordinate/sqrt(x_coordinate ** 2
C    &                                             + y_coordinate ** 2))
C           TODO: i took this from coor.f
            ros=solarGalactocentricDistance*solarGalactocentricDistance+
     &          coordinate_R(i)*coordinate_R(i)-2*coordinate_R(i)*
     &          solarGalactocentricDistance*dcos(coordinate_Theta(i))
            ros=dsqrt(ros)
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
            latitude=datan(zzx)
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
            write(421,"(6(es17.8e3,x),i1)") uu(i),
     &                                      vv(i),
     &                                      ww(i),
     &                                      rgac(i),
     &                                      longitude,
     &                                      latitude,
     &                                      disk_belonging(i)
         end do
      end if
      end subroutine
