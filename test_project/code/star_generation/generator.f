      subroutine gen(iseed,parameterOfSFR,areaOfSector,
     &           numberOfStarsInSample,galacticDiskAge,timeOfBurst,
     &           massReductionFactor,sfr_model)
C=======================================================================
C     Divides the SFR in intervals of time. In each interval, the total 
C     mass of stars is distributed. The mass of each star follows the 
C     distribution law given by IMF. The birth time is also calculated,
C     from which the height pattern and finally cylindrical 
C     z-coordinate are determined.
C-----------------------------------------------------------------------       
C           Input parameters:
C           iseed: seed for random generator
C           parameterOfSFR
C           areaOfSector: area of the considered sector in Kpc²
C           galacticDiskAge
C           timeOfBurst: how many Gyr ago burst of star gener-n happened
C-----------------------------------------------------------------------
C           Output parameters:
C           numberOfStarsInSample: number of participating stars
C=======================================================================
C     TODO: change to implicit none
      implicit real (a-h,m,o-z)

C     overloading intrinsic 'ran' function by our own RNG
      external ran
      real ran
      integer numberOfStars,
     &        nbins,
     &        iseed,
     &        sfr_model
      real zDistribution_zo,
     &                 deltaT_SFRThickDisk,
     &                 heightSFR_ThickDisk,
     &                 sheight, 
     &                 sheight_thick,
     &                 hDistr_zi,
     &                 hDistr_t,
     &                 hDistr_zf,
     &                 ttdisk,
     &                 tmdisk,
     &                 tmax,
     &                 tau,
     &                 ft,
     &                 fz

      parameter (numberOfStars = 6000000)
      parameter (zDistribution_zo = 1.0)
      parameter (deltaT_SFRThickDisk = 0.5)
      parameter (heightSFR_ThickDisk = 5.0)
C     Scale height of the thin disk?
      parameter (sheight = 0.250)
      parameter (hDistr_zi = 242.5)
      parameter (hDistr_t = 0.7)
      parameter (hDistr_zf = 0.250)
C     Scale height of the thick disk:
      parameter (sheight_thick=0.900)
C     Parameters of the thick disk
      parameter (ttdisk = 12.0)
      parameter (tmdisk = 10.0)
      parameter (tau = 2.0)
      
      integer numberOfStarsInSample, i, k, in
      real parameterOfSFR,areaOfSector,galacticDiskAge,cte,
     &                 deltat,massReductionFactor,psi
      real me,mgen,mrep,t,tf,to,xseed,xx,zmaxbin,zz
C     NOTE: ximf - mass from IMF by Salpeter, cnorm - const.of 
C           normalization of SFR. Both are from functions with same name
      real ximf,cnorm  
      double precision :: coordinate_Theta(numberOfStars),
     &                    coordinate_R(numberOfStars),
     &                    coordinate_Zcylindr(numberOfStars)
C     QUESTION: what is m? massInMainSequence? 
      real m(numberOfStars)
      real starBirthTime(numberOfStars)
      real heightPattern(numberOfStars)
      real flagOfWD(numberOfStars)
      integer disk_belonging(numberOfStars)

      common /tm/ starBirthTime,m
      common /coorcil/ coordinate_R,coordinate_Theta,coordinate_Zcylindr
      common /patron/ heightPattern
      common /index/ flagOfWD,numberOfWDs,disk_belonging

C------------------------------------------------------------------
C     ---   Calculating the mass, time of birth and z-coordinate of 
C           every star   ---
C------------------------------------------------------------------
C     --- Initialization of different variables ---

C     Variables referring to SFR:
C     Calculating the normalization constant of the SFR
      cte = cnorm(parameterOfSFR, galacticDiskAge)
C     QUESTION: is 'to' a time when stars begin to be born?      
      to = 0.0
C     QUESTION: why do we need tf if it's the same as galacticDiskAge?      
      tf = galacticDiskAge
      nbins = 5000
      deltat = (tf - to) / nbins  
C     NOTE: this k counter variable is used in GOTO-loop, which is bad
      k=0

C     ---  Reduction factor of the mass to be distributed ---
C     massReductionFactor=0.1 for 'fat' sample
C     massReductionFactor=3.0 individual samples

      write(6,*) '      Factor of mass reduction=',massReductionFactor

C     ---  Choosing properties of the thin/thick disk according to the 
C          corresponding line parameter ---

C     ---  Calculating the mass to be distributed at each interval ---        
      do 3 i = 1, nbins
          zmaxbin = 0.0
          psi = cte * deltat
          mrep = psi * areaOfSector * 1.06
          mrep = mrep / massReductionFactor
          mgen = 0.0

C         ---  Recent Burst   ---
C         NOTE: we don't need so many parameters
          tago = timeOfBurst
          tago = galacticDiskAge - tago
          tfin = galacticDiskAge

          tsfr = to + float(i-1) * deltat
      
          if(tsfr.ge.tago .and. tsfr.lt.tfin) then
              mrep = mrep * 5.0
          endif

C         ---  Calling to the IMF  ---
1         me = ximf(iseed)

C         QUESTION: what is it supposed to mean?
C         ---  Ya tenemos la masa  ---
          k=k+1

          if(k.eq.numberOfStars) then 
              write(6,*) '     ***  Dimension exceeded   ***'
              write(6,*) '***  Increase the reduction factor   ***'
              stop
          else
              continue
          endif
          m(k)=me
        
          mgen=mgen+me               
         
C         --- Birth time from SFR constant  --- 
C         --- Choosing with what model we generate stars ---
C         --- disk_belonging = 1 (thin disk), = 2 (thick disk)
          if (sfr_model == 2) then
              if (ran(iseed) .le. 0.08) then
                  disk_belonging(k) = 2
                  tmax = tmdisk * exp(-tmdisk/tau)
 33               ttry = ttdisk * ran(iseed)
                  ft = ttry * exp(-ttry/tau)
                  fz = tmax * ran(iseed)
                  if (fz .le. ft) then
                      starBirthTime(k) = ttry
                  else
                      goto 33
                  end if
              else
                  disk_belonging(k) = 1
                  xseed = deltat * ran(iseed)      
                  t = to + dfloat(i - 1)*deltat + xseed 
                  starBirthTime(k) = t
              end if
          else if (sfr_model == 1) then 
C             running ran once more to get same results as for model2 
              xseed = deltat * ran(iseed)
              xseed = deltat * ran(iseed)      
              t = to + dfloat(i - 1)*deltat + xseed 
              starBirthTime(k) = t
              disk_belonging(k) = 1
          else
              write(6,*) 'ERROR: wrong SFR model'
          end if
       
C         ---  Calculating the height pattern in kpc ---
C         ---  model of constant heightPattern  ---
          if (disk_belonging(k) == 1) then
              heightPattern(k) = sheight
          else if (disk_belonging(k) == 2) then
              heightPattern(k) = sheight_thick
          else
              write(6, *) "Error: wrong SFR model"
          end if
C         --- model of variable heightPattern  ---
C           heightPattern(k)=hDistr_zi*dexp(-starBirthTime(k)/hDistr_t)+
C      &                     hDistr_zf

C         --- Calculating z ---
C         TODO: delete this goto and put a loop here
2         xx=zDistribution_zo*heightPattern(k)*ran(iseed)
          if (xx.eq.0.0) goto 2 
          zz=heightPattern(k)*LOG(zDistribution_zo*heightPattern(k)/xx)

C         QUESTION: what is this?
C-------------------------------------------------------------------    
C       z-contstant       zz=0.240*ran(iseed)  
C-------------------------------------------------------------------    

          in=int(2.0*ran(iseed))

          coordinate_Zcylindr(k)=zz*dfloat(1-2*in)

          zmaxbin=dmax1(zmaxbin,coordinate_Zcylindr(k))

C         --- Checking if we have generated enough mass  ---
   
          if (mgen.lt.mrep) then
C             TODO: eleminate this goto
              goto 1
          else
              m(k)=m(k)-(mgen-mrep)
              mgen=mrep
          endif
3     continue
        
      numberOfStarsInSample=k
       
      write(6,*) '     Number of stars in sample=',numberOfStarsInSample

      return
      end


C***********************************************************************


C***********************************************************************
      function cnorm(parameterOfSFR, galacticDiskAge)
C=======================================================================
C     Calculating the normalization constant of the SFR
C=======================================================================
      implicit real (a-h,m,o-z)
      real sigma,parameterOfSFR,galacticDiskAge,cnorm
C     ---   Parameters   ---
      parameter (sigma=51.0)
C     ---   Calculating   ---
      cnorm=sigma/(parameterOfSFR*(1.0-
     &      exp(-galacticDiskAge/parameterOfSFR)))
      return
      end
C***********************************************************************


C***********************************************************************
      function ximf(iseed)
C=======================================================================
C
C     Calculating the mass following the IMF by Salpeter.
C    
C=======================================================================
      implicit real (a-h,m,o-z)
      external ran
      real ran

C     --- Parameter ---
C     QUESTION: what are these parameters?
      parameter(mmin=0.4)
      parameter(mmax=50.0)

      common /param/ fractionOfDB,galacticDiskAge,parameterIMF,
     &               parameterIFMR,timeOfBurst
      common /RSEED/ ISEED1,ISEED2
         
C     ---  IMF Salpeter 55  ---    

C     ---  Releasing a number (y-coordinate) according to a comparison 
C          function cte.   ---
      
      ymax=mmin**(parameterIMF)

 1    zy=ymax*ran(iseed)

C     ---  Releasing another number (x-coordinate)  ---
       
      zx=((mmax-mmin)*ran(iseed))+mmin

C     ---  Calculating the value of the IMF for m=zx ---

      zyimf=zx**(parameterIMF)
       
C     ---  Comparing IMF with the function cte  ---       
       
      if (zy.gt.zyimf) then
        goto 1
      else
        ximf=zx
      endif
       
      return
      end 
C***********************************************************************
