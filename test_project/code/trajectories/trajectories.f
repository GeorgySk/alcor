C     NOTE: the following subroutine is actually missed 
C***********************************************************************
C     TODO: rewrite
      subroutine traject(galacticDiskAge)
C=======================================================================
C     This subroutine calculates the trajectories of the WDs according 
C     to the z-axis. Using a 4th order Runge-Kuttta.
C
C     Revised in 22.09.07. by S. Torres
C
C-------------------------------------------------------------------
C     Input parameters
C       galacticDiskAge
C       numberOfStarsInSample
C=======================================================================
      implicit double precision (a-h,m,o-z)
      
      integer numberOfStars,numberOfWDs,njumps,n,i,NOK,NBAD
      double precision galacticDiskAge,wosun,hmin,eps,fcgys,tf,xcar,ycar
      double precision wo,zo,ecini,tb,tinc,htry,ti,ecinf
      double precision epotf,epoti,f

C     ---   Parameters  ---
      parameter (numberOfStars=6000000)
      parameter (wosun=-8.0d0)
      
C     ---   Dimensions   ---
      double precision uu(numberOfStars),vv(numberOfStars),
     &                 ww(numberOfStars)
      double precision coordinate_R(numberOfStars),
     &                 coordinate_Theta(numberOfStars),
     &                 coordinate_Zcylindr(numberOfStars)
      double precision starBirthTime(numberOfStars),m(numberOfStars)
      double precision flagOfWD(numberOfStars)
      double precision yscal(2),y(2),dydx(2)         
      double precision xpla(numberOfStars),ypla(numberOfStars)
      integer disk_belonging(numberOfStars)
        
C     ---   Commons   ---
      common /vel/ uu,vv,ww
      common /coorcil/ coordinate_R,coordinate_Theta,coordinate_Zcylindr
      common /tm/ starBirthTime,m
      common /index/ flagOfWD,numberOfWDs,disk_belonging         
C     NOTE: names are different in this block      
      common /plano/ xpla,ypla
      common /carte/ xcar,ycar

C     ---   Externals  ---
C     External - specifies procedures as external, and allows their 
C     symbolic names to be used as actual arguments.
      EXTERNAL DERIVS
      EXTERNAL RKQC
            
C     ---   Test   ---
      njumps=100
      hmin=0.0
      eps=1.0d-4
      n=2
      fcgys=(1.0d+9)*365.25*24.0*3600.0
      tf=galacticDiskAge*fcgys 

C     ---   Integrating trajectories   ---
      do 1 i=1,numberOfWDs
        xcar=xpla(i)
        ycar=ypla(i)
        wo=ww(i)+8.0
        zo=coordinate_Zcylindr(i)*(3.086d+16)
        ecini=0.5*wo*wo
        call epot(zo,epoti)
        tb=starBirthTime(i)
        tinc=(galacticDiskAge-tb)/dfloat(njumps)
C       ---  The time in seconds  ---
        htry=tinc*fcgys
C       ---  Initial conditions  ---
        y(1)=zo
        y(2)=wo
        dydx(1)=wo
        call fuerza(zo,f)
        dydx(2)=f
        ti=tb*fcgys 
C       ---  Calling to the Runge-Kutta integrator ---

        call ODEINT(y,n,ti,tf,eps,htry,hmin,NOK,NBAD,DERIVS,RKQC,yscal,
     &       y,dydx)      
        ecinf=0.5*y(2)*y(2)
        call epot(y(1),epotf)
        coordinate_Zcylindr(i)=y(1)/(3.086d+16)
        ww(i)=y(2)+wosun
 1    continue

      return
      end     
C***********************************************************************
     


C***********************************************************************      
C     TODO: rewrite      
      subroutine epot(z,e)
C=======================================================================
C
C     This function calculates the force along z-coordinate. 
C--------------------------------------------------------------------
C     Input parameters:
C       z: coordinate z (km)
C--------------------------------------------------------------------
C     Output parameters:
C       e: potential energy (km2/s2)
C=======================================================================
      implicit double precision (a-h,m,o-z)    
      
      double precision ro,vh,rc1,mc1,rc2,mc2,b,md1,md2,md3,a1,a2,a3,g
      double precision xcar,ycar,xpla,ypla,rpla,zsig,z,vh2,ro2
      double precision rc12,rc22,b2,rpla2,r2,poth,xa,xb,potc,xx
      double precision xd1,xd2,xd3,potd,potd1,potd2,potd3,pot,e

C     ---   Parameters   ---
      parameter(ro=8.5,vh=220.0)
      parameter(rc1=2.7,mc1=3.0d+09)
      parameter(rc2=0.42,mc2=1.6d+10)
      parameter(b=0.3)
      parameter(md1=6.6d+10,a1=5.81)
      parameter(md2=-2.9d+10,a2=17.43)
      parameter(md3=3.3d+09,a3=34.86)
      parameter(g=4.30026d-6)
      
C     ---   Common  
      common /carte/ xcar,ycar
                  
C     ---  Some calculations of interest   ---
      xpla=xcar
      ypla=ycar
      rpla=dsqrt(xpla*xpla+ypla*ypla)
      zsig=z
      z=dabs(z/(3.086d+16))                  
      vh2=0.5*vh*vh
      ro2=ro*ro
      rc12=rc1*rc1
      rc22=rc2*rc2
      b2=b*b
      rpla2=rpla*rpla
      r2=rpla2+z*z

C     ---   Calculating the potentials  ---           
C     ---   Dark halo  ---
      poth=vh2*dlog(r2+ro2)
C     ---   Central component  ---      
      xa=dsqrt(r2+rc12)
      xb=dsqrt(r2+rc22)
      potc=((-g*mc1)/xa)+((-g*mc2)/xb)
C     ---   Disk  ---
      xx=dsqrt(z*z+b2)
      xd1=rpla2+((a1+xx)*(a1+xx))
      xd2=rpla2+((a2+xx)*(a2+xx))
      xd3=rpla2+((a3+xx)*(a3+xx))
      potd1=(g*md1)/(dsqrt(xd1))
      potd2=(g*md2)/(dsqrt(xd2))
      potd3=(g*md3)/(dsqrt(xd3))
      potd=-potd1-potd2-potd3
C     ---   Total potential  ---
      pot=poth+potc+potd       
      e=pot
      z=zsig
      
      return
      end
C***********************************************************************



C     TODO: rewrite
C***********************************************************************
      subroutine fuerza(z,f)
C=======================================================================
C
C     This function calculates the force along z-coordinate.
C-------------------------------------------------------------------
C     Input parameters: z(km) xpla,ypla (kpc)
C-------------------------------------------------------------------
C     Output parameters: f
C=======================================================================
      implicit double precision (a-h,m,o-z)
      
      double precision ro,vh,rc1,mc1,rc2,mc2,b,md1,md2,md3,a1,a2,a3,g
      double precision xcar,ycar,xpla,ypla,rpla,zsig,z,vh2,ro2
      double precision rc12,rc22,b2,rpla2,r2,foh,foc1,foc2,foc,bzr
      double precision fod1,fod2,fod3,fod,ftot,fcv,f
      
C     ---   Parameters   ---
      parameter(ro=8.5,vh=220.0)
      parameter(rc1=2.7,mc1=3.0d+09)
      parameter(rc2=0.42,mc2=1.6d+10)
      parameter(b=0.3)
      parameter(md1=6.6d+10,a1=5.81)
      parameter(md2=-2.9d+10,a2=17.43)
      parameter(md3=3.3d+09,a3=34.86)
      parameter(g=4.30026d-6)
      
C     ---  Common  --- 
      common /carte/ xcar,ycar
                 
C     ---  Calculating some useful variables ---
      xpla=xcar
      ypla=ycar
      zsig=z
      z=dabs(z)/(3.086d+16)
      rpla2=xpla*xpla+ypla*ypla
      rpla=dsqrt(rpla2)
      vh2=vh*vh
      ro2=ro*ro
      rc12=rc1*rc1
      rc22=rc2*rc2
      b2=b*b
      r2=rpla*rpla+z*z
                   
C     ---   Calculating the forces  ---
C     ---   Dark halo  ---
      foh=vh2*z/(ro2+r2)
C     ---   Central component  ---      
      foc1=g*mc1*z/((rc12+r2)**(1.5d0))   
      foc2=g*mc2*z/((rc22+r2)**(1.5d0))      
      foc=foc1+foc2         
C     ---   Disk  ---
      bzr=dsqrt(b2+z*z)
      fod1=g*md1*z*(a1+bzr)/(bzr*(rpla2+(a1+bzr)*(a1+bzr))**(1.5d0))    
      fod2=g*md2*z*(a2+bzr)/(bzr*(rpla2+(a2+bzr)*(a2+bzr))**(1.5d0))     
      fod3=g*md3*z*(a3+bzr)/(bzr*(rpla2+(a3+bzr)*(a3+bzr))**(1.5d0))    
      fod=fod1+fod2+fod3
C     ---  Total force  ---
      ftot=foh+foc+fod 
C     ---  If we want the result in km/s²  ---
      fcv=1.0d0/(3.086d+16)
      ftot=fcv*dabs(ftot)
C     ---  The sign of z will be  ---
      f=-dsign(ftot,zsig)
      z=zsig
      
      return
      end
C***********************************************************************
