C***********************************************************************
C     TODO: rewrite
      function gasdev(iseed)
C=======================================================================
C
C     Returns a normally distributed deviate with zero mean and unit
C     variance.
C
C=======================================================================
      implicit double precision (a-h,m,o-z)
      external ran
      real ran
C     QUESTION: what is save?
      save
      
      integer iseed,iset
      double precision v1,v2,r,fac,gset,gasdev
      
      data iset/0/
      if (iset.eq.0) then
  1     v1=2.0d0*(ran(iseed))-1.0d0
        v2=2.0d0*(ran(iseed))-1.0d0
        r=v1*v1+v2*v2
        if (r.ge.1.0d0.or.r.eq.0.0d0) goto 1
        fac=dsqrt(-2.0d0*dlog(r)/r)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif

      return
      end
C***********************************************************************      
