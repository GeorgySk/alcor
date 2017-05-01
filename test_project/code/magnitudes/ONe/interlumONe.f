C***********************************************************************
C     TODO: rewrite
      subroutine interlumONe(tcool,mass,lum,c1,c2,c3,c4,c5,teff,xlog)
C=======================================================================
C
C     This subroutine interpolates luminosity of the DA ONe WD's
C     together with Johnson colors (U,B,V,R,I)
C
C     Created by S. Torres
C     Modifications 10/2012 (ER Cojocaru) - subtracted subroutine 
C                                             wdcoolone,
C                                    interpolation for subroutine interp
C
C-------------------------------------------------------------------
C     Input parameters:
C       tcool: cooling time
C       mass: mass of the WD
C-------------------------------------------------------------------
C     Output parameters:
C       lum: luminosity
C       c1,c2,c3,c4,c5: Johnson colors (U,B,V,R,I)
C       teff: effective temperature [K]
C       log g: logarithm of the superficial gravity
C
C=======================================================================
      implicit double precision (a-h,m,o-z)

C     ---   Declaration of variables   ---
      integer ncol,ncol2,nrow,nrow2,i,model,modlog
      double precision mass,tcool,y
      double precision mv,cbv,cvi,cvr,cuv,zlte,zlr,G
      double precision lum,c1,c2,c3,c4,c5,teff,xlog

C     ---   Parameters   ---
      parameter (ncol=6)
      parameter (ncol2=5)
      parameter (nrow=300)
      parameter (nrow2=300)

C     ---   Dimensions   ---
      integer ndatsone(ncol),ndatsone2(ncol2)
      double precision mtabone(ncol),mtabone2(ncol2)
      double precision mvtabone(ncol,nrow),ltabone(ncol,nrow)
      double precision lgtabone(ncol,nrow),lgtetabone(ncol,nrow)
      double precision bvtabone(ncol,nrow),vitabone(ncol,nrow)
      double precision vrtabone(ncol,nrow),uvtabone(ncol,nrow)
      double precision lgrtabone(ncol2,nrow2),lgt2tabone(ncol2,nrow2)
      double precision tprewd1(ncol),tprewd2(ncol2)
 
C     ---   Commons   ---
      common /fredone/ lgtabone,ltabone,mvtabone,lgtetabone
      common /fredone2/ mtabone,ndatsone
      common /colorsone/ bvtabone,vitabone,vrtabone,uvtabone
      common /newone/ lgrtabone,lgt2tabone
      common /newone2/ mtabone2,ndatsone2
      
C     ---   Interpolation  ---
      model=1
      modlog=0
      y=dlog10(tcool)+9.0

      do i=1,ncol
        tprewd1(i)=0.0
      end do

      do i=1,ncol2
        tprewd2(i)=0.0
      end do

      call interp(model,modlog,y,mass,ncol,ndatsone,lgtabone,tprewd1,
     &     mtabone,ltabone,lum)
      call interp(model,modlog,y,mass,ncol,ndatsone,lgtabone,tprewd1,
     &     mtabone,mvtabone,mv)
      call interp(model,modlog,y,mass,ncol,ndatsone,lgtabone,tprewd1,
     &     mtabone,bvtabone,cbv)
      call interp(model,modlog,y,mass,ncol,ndatsone,lgtabone,tprewd1,
     &     mtabone,vitabone,cvi)
      call interp(model,modlog,y,mass,ncol,ndatsone,lgtabone,tprewd1,
     &     mtabone,vrtabone,cvr)
      call interp(model,modlog,y,mass,ncol,ndatsone,lgtabone,tprewd1,
     &     mtabone,uvtabone,cuv)
      call interp(model,modlog,y,mass,ncol,ndatsone,lgtabone,tprewd1,
     &     mtabone,lgtetabone,zlte)
      call interp(model,modlog,y,mass,ncol2,ndatsone2,lgt2tabone,
     &     tprewd2,mtabone2,lgrtabone,zlr)

      teff=10.0**zlte

C   --- G in cm/s² kg; M in kg and R in cm ---
      G=6.67d-5
      xlog=dlog10(G)+dlog10(mass*(1.989d30))-2.0*(zlr+dlog10(6.696d10))
      c1=cuv+mv
      c2=cbv+mv
      c3=mv
      c4=mv-cvr
      c5=mv-cvi
       
      return
      end
C***********************************************************************
