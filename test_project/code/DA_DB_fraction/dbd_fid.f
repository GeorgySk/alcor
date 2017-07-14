C***********************************************************************
      subroutine dbd_fid(iseed,fractionOfDB,in)
C=======================================================================
C
C     This subroutine determines if the WD is DA or non-DA (DB),
C     using the model A (20/80) by Torres (A&A, 2010) 
C    
C     in=0 DA
C     in=1 non-DA
C
C     Created by ER Cojocaru (11/2012)
C
C=======================================================================
      implicit real (a-h,m,o-z)
      
C     ---   Declaration of variables  ---
      external ran
      real ran
      integer iseed,in
      real x,fractionOfDB
      
C     --- Fiducial model: 20% DB, 80% DA
      x=ran(iseed)
      if(x.lt.fractionOfDB) then 
        in=1
      else
        in=0
      endif

      return
      end
C***********************************************************************
