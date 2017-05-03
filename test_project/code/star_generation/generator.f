      subroutine gen(iseed,parameterOfSFR,areaOfSector,
     &           numberOfStarsInSample,galacticDiskAge,timeOfBurst,
     &           massReductionFactor,kinematicModel)
C=======================================================================
C
C     Divides the SFR in intervals of time. In each interval, the total 
C     mass of stars is distributed. The mass of each star follows the 
C     distribution law given by IMF. The birth time(t_b) is also
C     calculated from which the height pattern and finally cylindrical 
C     z-coordinate are determined.
C
C     NOTE: using version control system lets avoid lines like this
C     Revised by S. Torres 22/09/07
C     Modifications by E.R.Cojocaru 06/2013
C-----------------------------------------------------------------------       
C           Input parameters:
C           QUESTION: what does it mean?
C           iseed=entero generador de ran       
C           parameterOfSFR: parametro of the SFR        
C           areaOfSector: area of the considered sector in Kpc²
C-----------------------------------------------------------------------
C           Output parameters:
C           numberOfStarsInSample: number of participating stars
C=======================================================================
      implicit double precision (a-h,m,o-z)

C     ---   Declaration of variables  ---

C     QUESTION: what are the next 2 lines?
      external ran
      real ran
C     QUESTION: what is nbins and sheight?
      integer numberOfStars,nbins,iseed
      double precision zDistribution_zo,deltaT_SFRThickDisk,
     &                 heightSFR_ThickDisk,sheight
      double precision hDistr_zi,hDistr_t,hDistr_zf

      integer kinematicModel


C     ---   Parameters   ---


      parameter (numberOfStars=6000000)
      parameter (zDistribution_zo=1.0)
C     delta t SFR thick disk
      parameter (deltaT_SFRThickDisk=0.5)
C     height SFR thick disk
      parameter (heightSFR_ThickDisk=5.0)
      parameter (sheight=0.250)
      parameter (hDistr_zi=242.5)
      parameter (hDistr_t=0.7)
      parameter (hDistr_zf=0.250)

C     ---   Dimensions   ---
      
C     QUESTION: what are in, cte, deltat, psi?      
      integer numberOfStarsInSample,i,k,in
      double precision parameterOfSFR,areaOfSector,galacticDiskAge,cte,
     &                 deltat,massReductionFactor,psi
C     QUESTION: what are these variables?
      double precision me,mgen,mrep,t,tf,to,xseed,xx,zmaxbin,zz
C     NOTE; ximf - mass from IMF by Salpeter, cnorm - const.of 
C           normalization of SFR. Both are from functions with same name
      double precision ximf,cnorm  
      double precision coordinate_Theta(numberOfStars),
     &                 coordinate_R(numberOfStars),
     &                 coordinate_Zcylindr(numberOfStars)
C     QUESTION: what is m? massInMainSequence? 
      double precision m(numberOfStars)
      double precision starBirthTime(numberOfStars)
      double precision heightPattern(numberOfStars)

C     ---  Commons  ---
       
      common /tm/ starBirthTime,m
      common /coorcil/ coordinate_R,coordinate_Theta,coordinate_Zcylindr
      common /patron/ heightPattern


C------------------------------------------------------------------
C     ---   Calculating the mass, time of birth and z-coordinate of 
C           every star   ---
C------------------------------------------------------------------

C      --- Initialization de different variables ---

C     Variables referring to SFR:
C     Calculating the normalization constant of the SFR
      cte=cnorm(parameterOfSFR,galacticDiskAge)
C     QUESTION: is to - time when birth of stars begins?      
      to=0.0
C     QUESTION: why do we need tf if is the same as galacticDiskAge?      
      tf=galacticDiskAge
      nbins=5000
      deltat=(tf-to)/nbins  
C     NOTE: this k counter variable is used in GOTO-loop, which is bad
      k=0  

C      ---  Reduction factor of the mass to be distributed ---
C     massReductionFactor=0.1 for 'fat' sample
C     massReductionFactor=3.0 individual samples
C       massReductionFactor=0.5

      write(6,*) '      Factor of mass reduction=',massReductionFactor

C     ---  Choosing properties of the thin/thick disk according to the 
C          corresponding line parameter ---

C     ---  Calculating the mass to be distributed at each interval ---        



      do 3 i=1,nbins
C       QUESTION: what does all these parameters mean?       
        zmaxbin=0.0
        psi=cte*deltat
        mrep=psi*areaOfSector*1.0d6
        mrep=mrep/massReductionFactor
        mgen=0.0

C       ---  Recent Burst   ---
C       NOTE: we don't need so many parameters
        tago=timeOfBurst
        tago=galacticDiskAge-tago
        tfin=galacticDiskAge

        tsfr=to+dfloat(i-1)*deltat
      
        if(tsfr.ge.tago.and.tsfr.lt.tfin) then
          mrep=mrep*5.0
        endif

C       ---  Calling to the IMF  ---
1       me=ximf(iseed)

C       QUESTION: what is it supposed to mean?
C       ---  Ya tenemos la masa  ---
       
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
         
C       --- Birth time from SFR constant  --- 
        xseed=deltat*ran(iseed)      
        t=to+dfloat(i-1)*deltat+xseed 
        starBirthTime(k)=t 
       
C       ---  Calculating the height pattern in kpc ---
C       ---  modal of constant heightPattern  ---
      if (kinematicModel == 1) then
        heightPattern(k)=sheight
      else if (kinematicModel == 2) then
C       --- modal of variable heightPattern  ---
        heightPattern(k)=hDistr_zi*dexp(-starBirthTime(k)/hDistr_t)+
     &                    hDistr_zf
      else
        print*, "kinematic_model can be only 1 or 2"
      end if

C       --- Calculating z ---
C       TODO: delete this goto and put loop here
2       xx=zDistribution_zo*heightPattern(k)*ran(iseed)
        if (xx.eq.0.0) goto 2 
        zz=heightPattern(k)*DLOG(zDistribution_zo*heightPattern(k)/xx)

C       QUESTION: what is this?
C-------------------------------------------------------------------    
C       z-contstant       zz=0.240*ran(iseed)  
C-------------------------------------------------------------------    

        in=int(2.0*ran(iseed))

        coordinate_Zcylindr(k)=zz*dfloat(1-2*in)

        zmaxbin=dmax1(zmaxbin,coordinate_Zcylindr(k))

C       --- Checking if we have generated enough mass  ---
   
        if (mgen.lt.mrep) then
C         TODO: eleminate this goto
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
      function cnorm(parameterOfSFR,galacticDiskAge)
C=======================================================================
C
C     Calculating the normalization constant of the SFR
C    
C=======================================================================
      implicit double precision (a-h,m,o-z)
      double precision sigma,parameterOfSFR,galacticDiskAge,cnorm
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
      implicit double precision (a-h,m,o-z)
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