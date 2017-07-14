C***********************************************************************
C     TODO: rewrite
      subroutine intermag(mass,lumi,numberOfMassesWithColors,ntrk,mtrk,
     &           luminosity,color_U,color_B,color_V,color_R,color_I,c1,
     &           c2,c3,c4,c5)
C=======================================================================
C
C     This subroutine interpolates the luminosity of a DA WD
C     according to its mass and cooling time using the cooling sequence 
C     from input (corresponding to certain metallicity)
C
C     Modifications in 07/2012 (ER Cojocaru)
C
C-------------------------------------------------------------------
C     Input parameters:
C       mass: WD mass
C       lumi: luminosity
C       ltrx,color_U,color_B,color_V,color_R,color_I,ntrk,
C         numberOfMassesWithColors: information 
C           seq. colors (DA or DB)     
C
C-------------------------------------------------------------------
C     Output parameters:blanca
C       c1,c2,c3,c4,c5: Johnson colors
C
C=======================================================================
      implicit real (a-h,m,o-z)

C     ---   Declaraеion of variables   ---
      integer numberOfMassesWithColors,i,k,check1,check2,check3
      integer i1,i2,ntrk(numberOfMassesWithColors),ns1,ns2
      real mass,lumi,c1,c2,c3,c4,c5
      real c_1,c_2,a1,a2,b1,b2
      real mtrk(numberOfMassesWithColors)
      real luminosity(numberOfMassesWithColors,*),
     &                 color_U(numberOfMassesWithColors,*),
     &                 color_B(numberOfMassesWithColors,*)
      real color_V(numberOfMassesWithColors,*),
     &                 color_R(numberOfMassesWithColors,*),
     &                 color_I(numberOfMassesWithColors,*)

      check1=0
      check2=0
      check3=0

C     TODO: check if i mixed interpolation with extrapolation further
C     Smaller mass than known --> linear 2D extrapolation 
C     (using luminosity and mass)

      if(mass.le.mtrk(1)) then
        ns1=ntrk(1)
        ns2=ntrk(2)
C       QUESTION: what does it mean?
C       lum mayor que las conocidas --> linear 2D extrapolation
        if(lumi.gt.luminosity(1,1).OR.lumi.gt.luminosity(2,1)) then
          call extrap1(lumi,color_U(1,1),color_U(1,2),luminosity(1,1),
     &         luminosity(1,2),c_1)
          call extrap1(lumi,color_U(2,1),color_U(2,2),luminosity(2,1),
     &         luminosity(2,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c1)
          if(c1.lt.0.0) c1=0.0
          call extrap1(lumi,color_B(1,1),color_B(1,2),luminosity(1,1),
     &         luminosity(1,2),c_1)
          call extrap1(lumi,color_B(2,1),color_B(2,2),luminosity(2,1),
     &         luminosity(2,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c2)
          if(c2.lt.0.0) c2=0.0
          call extrap1(lumi,color_V(1,1),color_V(1,2),luminosity(1,1),
     &         luminosity(1,2),c_1)
          call extrap1(lumi,color_V(2,1),color_V(2,2),luminosity(2,1),
     &         luminosity(2,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c3)
          if(c3.lt.0.0) c3=0.0
          call extrap1(lumi,color_R(1,1),color_R(1,2),luminosity(1,1),
     &         luminosity(1,2),c_1)
          call extrap1(lumi,color_R(2,1),color_R(2,2),luminosity(2,1),
     &         luminosity(2,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c4)
          if(c4.lt.0.0) c4=0.0
          call extrap1(lumi,color_I(1,1),color_I(1,2),luminosity(1,1),
     &         luminosity(1,2),c_1)
          call extrap1(lumi,color_I(2,1),color_I(2,2),luminosity(2,1),
     &         luminosity(2,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c5)
          if(c5.lt.0.0) c5=0.0
          check3=1
          GOTO 45
C       lum menor que las conocidas --> linear 2D extrapolation
        elseif (lumi.lt.luminosity(1,ns1).OR.lumi.lt.luminosity(2,ns2)) 
     &  then
          call extrap1(lumi,color_U(1,ns1-1),color_U(1,ns1),
     &         luminosity(1,ns1-1),luminosity(1,ns1),c_1)
          call extrap1(lumi,color_U(2,ns2-1),color_U(2,ns2),
     &         luminosity(2,ns2-1),luminosity(2,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c1)
          if(c1.lt.0.0) c1=0.0
          call extrap1(lumi,color_B(1,ns1-1),color_B(1,ns1),
     &         luminosity(1,ns1-1),luminosity(1,ns1),c_1)
          call extrap1(lumi,color_B(2,ns2-1),color_B(2,ns2),
     &         luminosity(2,ns2-1),luminosity(2,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c2)
          if(c2.lt.0.0) c2=0.0
          call extrap1(lumi,color_V(1,ns1-1),color_V(1,ns1),
     &         luminosity(1,ns1-1),luminosity(1,ns1),c_1)
          call extrap1(lumi,color_V(2,ns2-1),color_V(2,ns2),
     &         luminosity(2,ns2-1),luminosity(2,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c3)
          if(c3.lt.0.0) c3=0.0
          call extrap1(lumi,color_R(1,ns1-1),color_R(1,ns1),
     &         luminosity(1,ns1-1),luminosity(1,ns1),c_1)
          call extrap1(lumi,color_R(2,ns2-1),color_R(2,ns2),
     &         luminosity(2,ns2-1),luminosity(2,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c4)
          if(c4.lt.0.0) c4=0.0      
          call extrap1(lumi,color_I(1,ns1-1),color_I(1,ns1),
     &         luminosity(1,ns1-1),luminosity(1,ns1),c_1)
          call extrap1(lumi,color_I(2,ns2-1),color_I(2,ns2),
     &         luminosity(2,ns2-1),luminosity(2,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c5)
          if(c5.lt.0.0) c5=0.0
          check3=1
          GOTO 45
C       lum contenida entre las conocidas
        else
          do i=1,ns1-1
            if(lumi.ge.luminosity(1,i+1).AND.lumi.le.luminosity(1,i)) 
     &      then
              i1=i
              check1=1
              GOTO 5
            end if
          end do
5         continue
          do i=1,ns2-1
            if(lumi.ge.luminosity(2,i+1).AND.lumi.le.luminosity(2,i)) 
     &      then
              i2=i
              check2=1
              GOTO 10
            end if
          end do
10        continue
          if(check1.eq.1.AND.check2.eq.1) then
            check3=1
            a1=lumi-luminosity(1,i1)
            a2=lumi-luminosity(2,i2)
            b1=luminosity(1,i1+1)-luminosity(1,i1)
            b2=luminosity(2,i2+1)-luminosity(2,i2)
            c_1=color_U(1,i1)+(color_U(1,i1+1)-color_U(1,i1))*a1/b1
            c_2=color_U(2,i2)+(color_U(2,i2+1)-color_U(2,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c1)
            c_1=color_B(1,i1)+(color_B(1,i1+1)-color_B(1,i1))*a1/b1
            c_2=color_B(2,i2)+(color_B(2,i2+1)-color_B(2,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c2)
            c_1=color_V(1,i1)+(color_V(1,i1+1)-color_V(1,i1))*a1/b1
            c_2=color_V(2,i2)+(color_V(2,i2+1)-color_V(2,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c3)
            c_1=color_R(1,i1)+(color_R(1,i1+1)-color_R(1,i1))*a1/b1
            c_2=color_R(2,i2)+(color_R(2,i2+1)-color_R(2,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c4)
            c_1=color_I(1,i1)+(color_I(1,i1+1)-color_I(1,i1))*a1/b1
            c_2=color_I(2,i2)+(color_I(2,i2+1)-color_I(2,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c5)
            GOTO 45
          end if
        end if
      end if
      
C     masa mayor que las conocidas --> linear interpolation
      if(mass.gt.mtrk(numberOfMassesWithColors)) then
        ns1=ntrk(numberOfMassesWithColors-1)
        ns2=ntrk(numberOfMassesWithColors)
C       lum mayor que las conocidas --> linear 2D extrapolation
        if(lumi.gt.luminosity(numberOfMassesWithColors-1,1).OR.lumi.gt.
     &  luminosity(numberOfMassesWithColors,1)) then
          call extrap1(lumi,color_U(numberOfMassesWithColors-1,1),
     &         color_U(numberOfMassesWithColors-1,2),
     &         luminosity(numberOfMassesWithColors-1,1),
     &         luminosity(numberOfMassesWithColors-1,2),c_1)
          call extrap1(lumi,color_U(numberOfMassesWithColors,1),
     &         color_U(numberOfMassesWithColors,2),
     &         luminosity(numberOfMassesWithColors,1),
     &         luminosity(numberOfMassesWithColors,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(numberOfMassesWithColors-1),
     &         mtrk(numberOfMassesWithColors),c1)
          if(c1.lt.0.0) c1=0.0    
          call extrap1(lumi,color_B(numberOfMassesWithColors-1,1),
     &         color_B(numberOfMassesWithColors-1,2),
     &         luminosity(numberOfMassesWithColors-1,1),
     &         luminosity(numberOfMassesWithColors-1,2),c_1)
          call extrap1(lumi,color_B(numberOfMassesWithColors,1),
     &         color_B(numberOfMassesWithColors,2),
     &         luminosity(numberOfMassesWithColors,1),
     &         luminosity(numberOfMassesWithColors,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(numberOfMassesWithColors-1),
     &         mtrk(numberOfMassesWithColors),c2)
          if(c2.lt.0.0) c2=0.0
          call extrap1(lumi,color_V(numberOfMassesWithColors-1,1),
     &         color_V(numberOfMassesWithColors-1,2),
     &         luminosity(numberOfMassesWithColors-1,1),
     &         luminosity(numberOfMassesWithColors-1,2),c_1)
          call extrap1(lumi,color_V(numberOfMassesWithColors,1),
     &         color_V(numberOfMassesWithColors,2),
     &         luminosity(numberOfMassesWithColors,1),
     &         luminosity(numberOfMassesWithColors,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(numberOfMassesWithColors-1),
     &         mtrk(numberOfMassesWithColors),c3)
          if(c3.lt.0.0) c3=0.0
          call extrap1(lumi,color_R(numberOfMassesWithColors-1,1),
     &         color_R(numberOfMassesWithColors-1,2),
     &         luminosity(numberOfMassesWithColors-1,1),
     &         luminosity(numberOfMassesWithColors-1,2),c_1)
          call extrap1(lumi,color_R(numberOfMassesWithColors,1),
     &         color_R(numberOfMassesWithColors,2),
     &         luminosity(numberOfMassesWithColors,1),
     &         luminosity(numberOfMassesWithColors,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(numberOfMassesWithColors-1),
     &         mtrk(numberOfMassesWithColors),c4)
          if(c4.lt.0.0) c4=0.0
          call extrap1(lumi,color_I(numberOfMassesWithColors-1,1),
     &         color_I(numberOfMassesWithColors-1,2),
     &         luminosity(numberOfMassesWithColors-1,1),
     &         luminosity(numberOfMassesWithColors-1,2),c_1)
          call extrap1(lumi,color_I(numberOfMassesWithColors,1),
     &         color_I(numberOfMassesWithColors,2),
     &         luminosity(numberOfMassesWithColors,1),
     &         luminosity(numberOfMassesWithColors,2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(numberOfMassesWithColors-1),
     &         mtrk(numberOfMassesWithColors),c5)
          if(c5.lt.0.0) c5=0.0
          check3=1
          GOTO 45
C       lum menor que las conocidas --> linear 2D extrapolation
        elseif (lumi.lt.luminosity(numberOfMassesWithColors-1,ns1).OR.
     &  lumi.lt.luminosity(numberOfMassesWithColors,ns1)) then
          call extrap1(lumi,color_U(numberOfMassesWithColors-1,ns1-1),
     &         color_U(numberOfMassesWithColors-1,ns1),
     &         luminosity(numberOfMassesWithColors-1,ns1-1),
     &         luminosity(numberOfMassesWithColors-1,ns1),c_1)
          call extrap1(lumi,color_U(numberOfMassesWithColors,ns2-1),
     &         color_U(numberOfMassesWithColors,ns2),
     &        luminosity(numberOfMassesWithColors,ns2-1),
     &        luminosity(numberOfMassesWithColors,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c1)
          call extrap1(lumi,color_B(numberOfMassesWithColors-1,ns1-1),
     &         color_B(numberOfMassesWithColors-1,ns1),
     &         luminosity(numberOfMassesWithColors-1,ns2-1),
     &         luminosity(numberOfMassesWithColors-1,ns1),c_1)
          call extrap1(lumi,color_B(numberOfMassesWithColors,ns2-1),
     &         color_B(numberOfMassesWithColors,ns2),
     &         luminosity(numberOfMassesWithColors,ns2-1),
     &         luminosity(numberOfMassesWithColors,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c2)
          call extrap1(lumi,color_V(numberOfMassesWithColors-1,ns1-1),
     &         color_V(numberOfMassesWithColors-1,ns1),
     &         luminosity(numberOfMassesWithColors-1,ns1-1),
     &         luminosity(numberOfMassesWithColors-1,ns1),c_1)
          call extrap1(lumi,color_V(numberOfMassesWithColors,ns2-1),
     &         color_V(numberOfMassesWithColors,ns2),
     &         luminosity(numberOfMassesWithColors,ns2-1),
     &         luminosity(numberOfMassesWithColors,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c3)
          call extrap1(lumi,color_R(numberOfMassesWithColors-1,ns1-1),
     &         color_R(numberOfMassesWithColors-1,ns1),
     &         luminosity(numberOfMassesWithColors-1,ns1-1),
     &         luminosity(numberOfMassesWithColors-1,ns1),c_1)
          call extrap1(lumi,color_R(numberOfMassesWithColors,ns2-1),
     &         color_R(numberOfMassesWithColors,ns2),
     &         luminosity(numberOfMassesWithColors,ns2-1),
     &         luminosity(numberOfMassesWithColors,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c4)
          call extrap1(lumi,color_I(numberOfMassesWithColors-1,ns1-1),
     &         color_I(numberOfMassesWithColors-1,ns1),
     &         luminosity(numberOfMassesWithColors-1,ns1-1),
     &         luminosity(numberOfMassesWithColors-1,ns1),c_1)
          call extrap1(lumi,color_I(numberOfMassesWithColors,ns2-1),
     &         color_I(numberOfMassesWithColors,ns2),
     &         luminosity(numberOfMassesWithColors,ns2-1),
     &         luminosity(numberOfMassesWithColors,ns2),c_2)
          call extrap1(mass,c_1,c_2,mtrk(1),mtrk(2),c5)      
          check3=1
          GOTO 45
C       lum contenida entre las conocidas
        else
          do i=1,ns1-1
            if(lumi.ge.luminosity(numberOfMassesWithColors-1,i+1).AND.
     &      lumi.le.luminosity(numberOfMassesWithColors-1,i)) then
              i1=i
              check1=1
              GOTO 15
            end if
          end do
15        continue
          do i=1,ns2-1
            if(lumi.ge.luminosity(numberOfMassesWithColors,i+1).AND.
     &      lumi.le.luminosity(numberOfMassesWithColors,i)) then
              i2=i
              check2=1
              GOTO 20
            end if
          end do
20        continue
          if(check1.eq.1.AND.check2.eq.1) then
            check3=1
            a1=lumi-luminosity(numberOfMassesWithColors-1,i1)
            a2=lumi-luminosity(numberOfMassesWithColors,i2)
            b1=luminosity(numberOfMassesWithColors-1,i1+1)-
     &         luminosity(numberOfMassesWithColors-1,i1)
            b2=luminosity(numberOfMassesWithColors,i2+1)-
     &         luminosity(numberOfMassesWithColors,i2)
            c_1=color_U(numberOfMassesWithColors-1,i1)+
     &          (color_U(numberOfMassesWithColors-1,i1+1)-
     &          color_U(numberOfMassesWithColors-1,i1))*a1/b1
            c_2=color_U(numberOfMassesWithColors,i2)+
     &          (color_U(numberOfMassesWithColors,i2+1)-
     &          color_U(numberOfMassesWithColors,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(1),
     &           mtrk(numberOfMassesWithColors),c1)
            c_1=color_B(numberOfMassesWithColors-1,i1)+
     &          (color_B(numberOfMassesWithColors-1,i1+1)-
     &          color_B(numberOfMassesWithColors-1,i1))*a1/b1
            c_2=color_B(numberOfMassesWithColors,i2)+
     &          (color_B(numberOfMassesWithColors,i2+1)-
     &          color_B(numberOfMassesWithColors,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(numberOfMassesWithColors-1),
     &           mtrk(numberOfMassesWithColors),c2)
            c_1=color_V(numberOfMassesWithColors-1,i1)+
     &          (color_V(numberOfMassesWithColors-1,i1+1)-
     &          color_V(numberOfMassesWithColors-1,i1))*a1/b1
            c_2=color_V(numberOfMassesWithColors,i2)+
     &          (color_V(numberOfMassesWithColors,i2+1)-
     &          color_V(numberOfMassesWithColors,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(1),
     &           mtrk(numberOfMassesWithColors),c3)
            c_1=color_R(numberOfMassesWithColors-1,i1)+
     &          (color_R(numberOfMassesWithColors-1,i1+1)-
     &          color_R(numberOfMassesWithColors-1,i1))*a1/b1
            c_2=color_R(numberOfMassesWithColors,i2)+
     &          (color_R(numberOfMassesWithColors,i2+1)-
     &          color_R(numberOfMassesWithColors,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(numberOfMassesWithColors-1),
     &           mtrk(numberOfMassesWithColors),c4)
            c_1=color_I(numberOfMassesWithColors-1,i1)+
     &          (color_I(numberOfMassesWithColors-1,i1+1)-
     &          color_I(numberOfMassesWithColors-1,i1))*a1/b1
            c_2=color_I(numberOfMassesWithColors,i2)+
     &          (color_I(numberOfMassesWithColors,i2+1)-
     &          color_I(numberOfMassesWithColors,i2))*a2/b2
            call extrap1(mass,c_1,c_2,mtrk(numberOfMassesWithColors-1),
     &           mtrk(numberOfMassesWithColors),c5)
            GOTO 45
          end if
        end if
      end if

C     masa contenida entre las conocidas --> linear interpolation
      do k=1,numberOfMassesWithColors-1
        if(mass.gt.mtrk(k).AND.mass.le.mtrk(k+1)) then
          ns1=ntrk(k)
          ns2=ntrk(k+1)
C         lum mayor que las conocidas --> linear 2D extrapolation
          if(lumi.gt.luminosity(k,1).OR.lumi.gt.luminosity(k+1,1)) then      
            call extrap1(lumi,color_U(k,1),color_U(k,2),luminosity(k,1),
     &           luminosity(k,2),c_1)
            call extrap1(lumi,color_U(k+1,1),color_U(k+1,2),
     &           luminosity(k+1,1),luminosity(k+1,2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c1)
            if(c1.lt.0.0) c1=0.0      
            call extrap1(lumi,color_B(k,1),color_B(k,2),luminosity(k,1),
     &           luminosity(k,2),c_1)
            call extrap1(lumi,color_B(k+1,1),color_B(k+1,2),
     &           luminosity(k+1,1),luminosity(k+1,2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c2)
            if(c2.lt.0.0) c2=0.0
            call extrap1(lumi,color_V(k,1),color_V(k,2),luminosity(k,1),
     &           luminosity(k,2),c_1)
            call extrap1(lumi,color_V(k+1,1),color_V(k+1,2),
     &           luminosity(k+1,1),luminosity(k+1,2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c3)
            if(c3.lt.0.0) c3=0.0
            call extrap1(lumi,color_R(k,1),color_R(k,2),luminosity(k,1),
     &           luminosity(k,2),c_1)
            call extrap1(lumi,color_R(k+1,1),color_R(k+1,2),
     &           luminosity(k+1,1),luminosity(k+1,2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c4)
            if(c4.lt.0.0) c4=0.0
            call extrap1(lumi,color_I(k,1),color_I(k,2),luminosity(k,1),
     &           luminosity(k,2),c_1)
            call extrap1(lumi,color_I(k+1,1),color_I(k+1,2),
     &           luminosity(k+1,1),luminosity(k+1,2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c5)
            if(c5.lt.0.0) c5=0.0
            check3=1      
            GOTO 45
C         lum menor que las conocidas --> linear 2D extrapolation
          elseif (lumi.lt.luminosity(k,ns1).OR.lumi.lt.
     &    luminosity(k+1,ns2)) then
            call extrap1(lumi,color_U(k,ns1-1),color_U(k,ns1),
     &           luminosity(k,ns1-1),luminosity(k,ns1),c_1)
            call extrap1(lumi,color_U(k+1,ns2-1),color_U(k+1,ns2),
     &           luminosity(k+1,ns2-1),luminosity(k+1,ns2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c1)
            call extrap1(lumi,color_B(k,ns1-1),color_B(k,ns1),
     &           luminosity(k,ns1-1),luminosity(k,ns1),c_1)
            call extrap1(lumi,color_B(k+1,ns2-1),color_B(k+1,ns2),
     &           luminosity(k+1,ns2-1),luminosity(k+1,ns2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c2)
            call extrap1(lumi,color_V(k,ns1-1),color_V(k,ns1),
     &          luminosity(k,ns1-1),luminosity(k,ns1),c_1)
            call extrap1(lumi,color_V(k+1,ns2-1),color_V(k+1,ns2),
     &           luminosity(k+1,ns2-1),luminosity(k+1,ns2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c3)
            call extrap1(lumi,color_R(k,ns1-1),color_R(k,ns1),
     &           luminosity(k,ns1-1),luminosity(k,ns1),c_1)
            call extrap1(lumi,color_R(k+1,ns2-1),color_R(k+1,ns2),
     &           luminosity(k+1,ns2-1),luminosity(k+1,ns2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c4)      
            call extrap1(lumi,color_I(k,ns1-1),color_I(k,ns1),
     &           luminosity(k,ns1-1),luminosity(k,ns1),c_1)
            call extrap1(lumi,color_I(k+1,ns2-1),color_I(k+1,ns2),
     &           luminosity(k+1,ns2-1),luminosity(k+1,ns2),c_2)
            call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c5)      
            check3=1
            GOTO 45
C         lum contenida entre las conocidas
          else
            do i=1,ns1-1
              if(lumi.ge.luminosity(k,i+1).AND.lumi.le.luminosity(k,i)) 
     &        then
                i1=i
                check1=1
                GOTO 25
              end if
            end do      
25          continue      
            do i=1,ns2-1
              if(lumi.ge.luminosity(k+1,i+1).AND.lumi.le.
     &        luminosity(k+1,i)) then
               i2=i
                check2=1
                GOTO 30
              end if
            end do
30          continue
            if(check1.eq.1.AND.check2.eq.1) then
              check3=1
              a1=lumi-luminosity(k,i1)
              a2=lumi-luminosity(k+1,i2)
              b1=luminosity(k,i1+1)-luminosity(k,i1)
              b2=luminosity(k+1,i2+1)-luminosity(k+1,i2)
              c_1=color_U(k,i1)+(color_U(k,i1+1)-color_U(k,i1))*a1/b1
              c_2=color_U(k+1,i2)+(color_U(k+1,i2+1)-color_U(k+1,i2))*
     &            a2/b2
              call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c1)
              c_1=color_B(k,i1)+(color_B(k,i1+1)-color_B(k,i1))*a1/b1
              c_2=color_B(k+1,i2)+(color_B(k+1,i2+1)-color_B(k+1,i2))*
     &            a2/b2
             call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c2)
              c_1=color_V(k,i1)+(color_V(k,i1+1)-color_V(k,i1))*a1/b1
              c_2=color_V(k+1,i2)+(color_V(k+1,i2+1)-color_V(k+1,i2))*
     &            a2/b2
              call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c3)
              c_1=color_R(k,i1)+(color_R(k,i1+1)-color_R(k,i1))*a1/b1
              c_2=color_R(k+1,i2)+(color_R(k+1,i2+1)-color_R(k+1,i2))*
     &            a2/b2
              call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c4)
              c_1=color_I(k,i1)+(color_I(k,i1+1)-color_I(k,i1))*a1/b1
              c_2=color_I(k+1,i2)+(color_I(k+1,i2+1)-color_I(k+1,i2))*
     &            a2/b2
              call extrap1(mass,c_1,c_2,mtrk(k),mtrk(k+1),c5)
              GOTO 45
            end if
          end if
        end if
      end do
45    continue

      if(check3.eq.0) write(1007,*) 'ERROR',lumi,mass

      return
      end
C***********************************************************************


C***********************************************************************
C     TODO:rewrite      
      subroutine extrap1(lumi,x1,x2,l1,l2,c)
      
      implicit real (a-h,m,o-z)
      real lumi,x1,x2,l1,l2,c,s,b
      
      s=(x2-x1)/(l2-l1)
      b=x1-s*l1
      c=s*lumi+b
         
      return
      end
C***********************************************************************
