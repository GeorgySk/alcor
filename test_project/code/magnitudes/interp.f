C***********************************************************************
C     TODO: rewrite      
      subroutine interp(model,modlog,y,m,ncol,ntrk,ttrkk,tprewd,mtrk,
     &           xtrk,xout)
C=======================================================================
C
C     This subroutine interpolates a certain related measure with the EB
C     DB or ONe, according to input data, considering all the cases of
C     interpolation and extrapolation needed.
C
C     Created by ER Cojocaru (10/2012)
C
C-------------------------------------------------------------------
C     Input parameters:
C       model: DA,DB (0) or ONe (1)
C       QUESTION: what does it mean?
C       modlog: hay que usar log de la medida (1 - teff) o no (0) (DA/DB)
C       y: cooling time
C       m: mass of the WD
C       ncol: number of sequences(columns) available
C       ntrk: number de rows for sequence
C       ttrkk: known cooling times in the sequence
C       mtrk: known masses in the sequence
C       xtrk: values of the measure of interest known in sequence
C-------------------------------------------------------------------
C     Output parameters:
C       xout: measure of interest
C=======================================================================
      implicit real (a-h,m,o-z)

C     ---   Declaration de variables   ---
      integer ncol,ii,jj,k,ns,j1,j2,case1,case2,model,modlog
C     TIME      
      real y,y1,y2,y3,y4,ym1,ym2
C     MASS
      real m,m1,m2
C     Measure of interest
      real x1,x2,x3,x4,xm1,xm2,xout
      
      real deltf,den,s,b,t

C     ---   Dimensions   ---
      integer ntrk(ncol)
      real mtrk(ncol),ttrkk(ncol,*)
      real tprewd(ncol),xtrk(ncol,*)
     
C     ---  Interpolating measure from tcool and mwd  ---

C     ---   Mass less than minimum   ---
      if (m .lt. mtrk(1)) then
        jj=1
        m1=mtrk(jj)
C       --- Checking if the time is less than the minimum---
        if(y.lt.ttrkk(jj,1)) then
C         making linear extrapolation in log-log
          if(model.eq.0) then
            if(modlog.eq.1) then
              deltf=log10((ttrkk(jj,2))/(ttrkk(jj,1)))
              s=log10(xtrk(jj,2)/xtrk(jj,1))/deltf
              b=log10(xtrk(jj,2))-s*log10(ttrkk(jj,2))      
              xm1=10.0**(s*log10(y)+b)
            else
              deltf=log10((ttrkk(jj,2))/(ttrkk(jj,1)))
              s=(xtrk(jj,2)-xtrk(jj,1))/deltf
              b=xtrk(jj,2)-s*log10(ttrkk(jj,2))      
              xm1=s*log10(y)+b
            endif
          else
            deltf=ttrkk(jj,2)-ttrkk(jj,1)
            s=(xtrk(jj,2)-xtrk(jj,1))/deltf
            b=xtrk(jj,2)-s*ttrkk(jj,2)
            xm1=s*y+b
          end if
          goto 20
        endif      
C       --- Checking if the time is greater than maximum ---
        ns=ntrk(jj)
        if(y.gt.ttrkk(jj,ns)) then
C         making linear extrapolation in log-log
          if(model.eq.0) then
            if(modlog.eq.1) then
              deltf=log10((ttrkk(jj,ns))/(ttrkk(jj,ns-1)))
              s=log10(xtrk(jj,ns)/xtrk(jj,ns-1))/deltf
              b=log10(xtrk(jj,ns))-s*log10(ttrkk(jj,ns))      
              xm1=10.0**(s*log10(y)+b)
            else
              deltf=log10((ttrkk(jj,ns))/(ttrkk(jj,ns-1)))
              s=(xtrk(jj,ns)-xtrk(jj,ns-1))/deltf
              b=xtrk(jj,ns)-s*log10(ttrkk(jj,ns))      
              xm1=s*log10(y)+b
            endif
          else
            deltf=ttrkk(jj,ns)-ttrkk(jj,ns-1)
            s=(xtrk(jj,ns)-xtrk(jj,ns-1))/deltf
            b=xtrk(jj,ns)-s*ttrkk(jj,ns)
            xm1=s*y+b
          end if
          goto 20
        endif
        do k=1,ns-1
          ii=k
C         ---   Interpolation of the tiempo   ---
          if (y .ge. ttrkk(jj,ii) .and. y .lt. ttrkk(jj,ii+1)) then            
            y1=ttrkk(jj,ii)
            y2=ttrkk(jj,ii+1)
            x1=xtrk(jj,ii)
            x2=xtrk(jj,ii+1)       
            deltf=(y-y1)/(y2-y1)
            xm1=x1+deltf*(x2-x1)   
            goto 20
          end if
        end do
20      CONTINUE
        jj=2
        m2=mtrk(jj)
C       --- Checking if the time is less than minimum---
        if(y.lt.ttrkk(jj,1)) then
C         hacemos un extrapolacion lineal en log-log
          if(model.eq.0) then
            if(modlog.eq.1) then
              deltf=log10((ttrkk(jj,2))/(ttrkk(jj,1)))
              s=log10(xtrk(jj,2)/xtrk(jj,1))/deltf
              b=log10(xtrk(jj,2))-s*log10(ttrkk(jj,2))      
              xm2=10.0**(s*log10(y)+b)
            else
              deltf=log10((ttrkk(jj,2))/(ttrkk(jj,1)))
              s=(xtrk(jj,2)-xtrk(jj,1))/deltf
              b=xtrk(jj,2)-s*log10(ttrkk(jj,2))      
              xm2=s*log10(y)+b
            endif
          else
            deltf=ttrkk(jj,2)-ttrkk(jj,1)
            s=(xtrk(jj,2)-xtrk(jj,1))/deltf
            b=xtrk(jj,2)-s*ttrkk(jj,2)
            xm2=s*y+b
          end if
          goto 30
        endif 
C       --- Checking if the time is greater than maximum---
        ns=ntrk(jj)
        if(y.gt.ttrkk(jj,ns)) then
C         making linear extrapolation in log-log
          if(model.eq.0) then
            if(modlog.eq.1) then
              deltf=log10((ttrkk(jj,ns))/(ttrkk(jj,ns-1)))
              s=log10(xtrk(jj,ns)/xtrk(jj,ns-1))/deltf
              b=log10(xtrk(jj,ns))-s*log10(ttrkk(jj,ns))      
              xm2=10.0**(s*log10(y)+b)
            else
              deltf=log10((ttrkk(jj,ns))/(ttrkk(jj,ns-1)))
              s=(xtrk(jj,ns)-xtrk(jj,ns-1))/deltf
              b=xtrk(jj,ns)-s*log10(ttrkk(jj,ns))      
              xm2=s*log10(y)+b
            endif
          else
            deltf=ttrkk(jj,ns)-ttrkk(jj,ns-1)
            s=(xtrk(jj,ns)-xtrk(jj,ns-1))/deltf
            b=xtrk(jj,ns)-s*ttrkk(jj,ns)
            xm2=s*y+b
          end if
          goto 30
        endif
        do k=1,ntrk(jj)-1
          ii=k
C         ---   Interpolation of the time   ---
          if (y .ge. ttrkk(jj,ii) .and. y .lt. ttrkk(jj,ii+1)) then
            y1=ttrkk(jj,ii)
            y2=ttrkk(jj,ii+1)
            x1=xtrk(jj,ii)
            x2=xtrk(jj,ii+1)  
            deltf=(y-y1)/(y2-y1)
            xm2=x1+deltf*(x2-x1) 
            goto 30
          end if
        end do
30      CONTINUE
        s=(xm2-xm1)/(m2-m1)
        t=xm2-s*m2
        xout=s*m+t
        goto 100
      end if
        
C     ---   Mass greater than maximum---
      if (m .ge. mtrk(ncol)) then
        jj=ncol-1
        m1=mtrk(jj)
C       --- Checking if time is less than minimum---
        if(y.lt.ttrkk(jj,1)) then
C         making linear extrapolation in log-log
          if(model.eq.0) then
            if(modlog.eq.1) then
              deltf=log10((ttrkk(jj,2)/(ttrkk(jj,1))))
              s=log10(xtrk(jj,2)/xtrk(jj,1))/deltf
              b=log10(xtrk(jj,2))-s*log10(ttrkk(jj,2))      
              xm1=10.0**(s*log10(y)+b)
            else
              deltf=log10((ttrkk(jj,2))/(ttrkk(jj,1)))
              s=(xtrk(jj,2)-xtrk(jj,1))/deltf
              b=xtrk(jj,2)-s*log10(ttrkk(jj,2))
              xm1=s*log10(y)+b
            endif
          else
            deltf=ttrkk(jj,2)-ttrkk(jj,1)
            s=(xtrk(jj,2)-xtrk(jj,1))/deltf
            b=xtrk(jj,2)-s*ttrkk(jj,2)
            xm1=s*y+b
          end if
          goto 200
        endif      
C       --- Checking if time is greater than maximum---
        ns=ntrk(jj)
        if(y.gt.ttrkk(jj,ns)) then
C         making linear extrapolation in log-log
          if(model.eq.0) then
            if(modlog.eq.1) then
              deltf=log10((ttrkk(jj,ns))/(ttrkk(jj,ns-1)))
              s=log10(xtrk(jj,ns)/xtrk(jj,ns-1))/deltf
              b=log10(xtrk(jj,ns))-s*log10(ttrkk(jj,ns))      
              xm1=10.0**(s*log10(y)+b)
            else
              deltf=log10((ttrkk(jj,ns))/(ttrkk(jj,ns-1)))
              s=(xtrk(jj,ns)-xtrk(jj,ns-1))/deltf
              b=xtrk(jj,ns)-s*log10(ttrkk(jj,ns))      
              xm1=s*log10(y)+b
            endif
          else
            deltf=ttrkk(jj,ns)-ttrkk(jj,ns-1)
            s=(xtrk(jj,ns)-xtrk(jj,ns-1))/deltf
            b=xtrk(jj,ns)-s*ttrkk(jj,ns)
            xm1=s*y+b
          end if
          goto 200
        endif
        do k=1,ns-1
          ii=k
C         ---   Interpolation of time   ---
          if (y .ge. ttrkk(jj,ii) .and. y .lt. ttrkk(jj,ii+1)) then
            y1=ttrkk(jj,ii)
            y2=ttrkk(jj,ii+1)
            x1=xtrk(jj,ii)
            x2=xtrk(jj,ii+1)       
            deltf=(y-y1)/(y2-y1)
            xm1=x1+deltf*(x2-x1)  
            goto 200
          end if
        end do
200     CONTINUE
        jj=ncol
        m2=mtrk(jj)
C       --- Checking if time is less than minimum ---
        if(y.lt.ttrkk(jj,1)) then
C       making linear extrapolation in log-log
          if(model.eq.0) then
            if(modlog.eq.1) then
              deltf=log10((ttrkk(jj,2))/(ttrkk(jj,1)))
              s=log10(xtrk(jj,2)/xtrk(jj,1))/deltf
              b=log10(xtrk(jj,2))-s*log10(ttrkk(jj,2))      
              xm2=10.0**(s*log10(y)+b)
            else
              deltf=log10((ttrkk(jj,2))/(ttrkk(jj,1)))
              s=(xtrk(jj,2)-xtrk(jj,1))/deltf
              b=xtrk(jj,2)-s*log10(ttrkk(jj,2))      
              xm2=s*log10(y)+b
            endif
          else
            deltf=ttrkk(jj,2)-ttrkk(jj,1)
            s=(xtrk(jj,2)-xtrk(jj,1))/deltf
            b=xtrk(jj,2)-s*ttrkk(jj,2)
            xm2=s*y+b
          end if      
          goto 300
        endif 
C       --- Checking if time is greater than maximum ---
        ns=ntrk(jj)
        if(y.gt.ttrkk(jj,ns)) then
C         making linear extrapolation in log-log
          if(model.eq.0) then
            if(modlog.eq.1) then
              deltf=log10((ttrkk(jj,ns))/(ttrkk(jj,ns-1)))
              s=log10(xtrk(jj,ns)/xtrk(jj,ns-1))/deltf
              b=log10(xtrk(jj,ns))-s*log10(ttrkk(jj,ns))      
              xm2=10.0**(s*log10(y)+b)
            else
              deltf=log10((ttrkk(jj,ns))/(ttrkk(jj,ns-1)))
              s=(xtrk(jj,ns)-xtrk(jj,ns-1))/deltf
              b=xtrk(jj,ns)-s*log10(ttrkk(jj,ns))      
              xm2=s*log10(y)+b
            endif
          else
            deltf=ttrkk(jj,ns)-ttrkk(jj,ns-1)
            s=(xtrk(jj,ns)-xtrk(jj,ns-1))/deltf
            b=xtrk(jj,ns)-s*ttrkk(jj,ns)
            xm2=s*y+b
          end if
          goto 300
        endif
        do k=1,ntrk(jj)-1
          ii=k
C         ---  Interpolation of time   ---
          if (y .ge. ttrkk(jj,ii) .and. y .lt. ttrkk(jj,ii+1)) then      
            y1=ttrkk(jj,ii)
            y2=ttrkk(jj,ii+1)
            x1=xtrk(jj,ii)
            x2=xtrk(jj,ii+1)  
            deltf=(y-y1)/(y2-y1)
            xm2=x1+deltf*(x2-x1)  
            goto 300
          end if
        end do
300     CONTINUE
        s=(xm2-xm1)/(m2-m1)
        t=xm2-s*m2
        xout=s*m+t
        goto 100
      end if
  
C     ---   Search for masses of interpolation ---
      do 3 k=1,ncol-1
        j1=k
        if (m .ge. mtrk(k) .and. m .lt. mtrk(k+1)) then
          goto 4
        end if
3     continue
4     j2=j1+1

C     ---   Search for the times of interpolation (1)   ---
C     ---   Times less than minimum  ---
      if (y. lt. ttrkk(j1,1)) then
C       making linear extrapolation in log-log
        if(model.eq.0) then
          if(modlog.eq.1) then
            deltf=log10((ttrkk(j1,2))/(ttrkk(j1,1)))
            s=log10(xtrk(j1,2)/xtrk(j1,1))/deltf
            b=log10(xtrk(j1,2))-s*log10(ttrkk(j1,2))
            x1=10.0**(s*log10(y)+b)
          else
            deltf=log10((ttrkk(j1,2))/(ttrkk(j1,1)))
            s=(xtrk(j1,2)-xtrk(j1,1))/deltf
            b=xtrk(j1,2)-s*log10(ttrkk(j1,2))
            x1=s*log10(y)+b
          end if
        else
          deltf=ttrkk(j1,2)-ttrkk(j1,1)
          s=(xtrk(j1,2)-xtrk(j1,1))/deltf
          b=xtrk(j1,2)-s*ttrkk(j1,2)
          x1=s*y+b
        end if
        case1=1
        goto 6      
      end if

C     ---   Times greater than maximum ---
      if(y .ge. ttrkk(j1,ntrk(j1))) then            
        ns=ntrk(j1)
C       making linear extrapolation in log-log
        if(model.eq.0) then
          if(modlog.eq.1) then
            deltf=log10((ttrkk(j1,ns))/(ttrkk(j1,ns-1)))
            s=log10(xtrk(j1,ns)/xtrk(j1,ns-1))/deltf
            b=log10(xtrk(j1,ns))-s*log10(ttrkk(j1,ns))
            x1=10.0**(s*log10(y)+b)
          else
            deltf=log10((ttrkk(j1,ns))/(ttrkk(j1,ns-1)))
            s=(xtrk(j1,ns)-xtrk(j1,ns-1))/deltf
            b=xtrk(j1,ns)-s*log10(ttrkk(j1,ns))
            x1=s*log10(y)+b
          end if
        else
          deltf=ttrkk(j1,ns)-ttrkk(j1,ns-1)
          s=(xtrk(j1,ns)-xtrk(j1,ns-1))/deltf
          b=xtrk(j1,ns)-s*ttrkk(j1,ns)
          x1=s*y+b
        end if
        case1=1
        goto 6
      end if

C     ---   Times between the minimum and maximum ---
      do 5 k=1,ntrk(j1)-1
        if (y .ge. ttrkk(j1,k) .and. y .lt. ttrkk(j1,k+1)) then
          y1=ttrkk(j1,k  )
          y2=ttrkk(j1,k+1)
          x1=xtrk(j1,k  )
          x2=xtrk(j1,k+1)
          case1=0
          goto 6
        end if
5     continue
 
            
C     ---   Search for interpolation times (2)   ---
C     ---   Times less than minimum---
6     if (y. lt. ttrkk(j2,1)) then
C       making linear extrapolation in log-log
        if(model.eq.0) then
          if(modlog.eq.1) then
            deltf=log10((ttrkk(j2,2))/(ttrkk(j2,1)))
            s=log10(xtrk(j2,2)/xtrk(j2,1))/deltf
            b=log10(xtrk(j2,2))-s*log10(ttrkk(j2,2))
            x3=10.0**(s*log10(y)+b)
          else
            deltf=log10((ttrkk(j2,2))/(ttrkk(j2,1)))
            s=(xtrk(j2,2)-xtrk(j2,1))/deltf
            b=xtrk(j2,2)-s*log10(ttrkk(j2,2))
            x3=s*log10(y)+b
          end if
        else
          deltf=ttrkk(j2,2)-ttrkk(j2,1)
            s=(xtrk(j2,2)-xtrk(j2,1))/deltf
            b=xtrk(j2,2)-s*ttrkk(j2,2)
            x3=s*y+b
        end if
        case2=1
        goto 8
      endif

C     ---   Times greater than maximum   ---
      if(y .ge. ttrkk(j2,ntrk(j2))) then  
        ns=ntrk(j2)
C       making linear extrapolation in log-log
        if(model.eq.0) then
          if(modlog.eq.1) then
            deltf=log10((ttrkk(j2,ns))/(ttrkk(j2,ns-1)))
            s=log10(xtrk(j2,ns)/xtrk(j2,ns-1))/deltf
            b=log10(xtrk(j2,ns))-s*log10(ttrkk(j2,ns))
            x3=10.0**(s*log10(y)+b)
          else
            deltf=log10((ttrkk(j2,ns))/(ttrkk(j2,ns-1)))
            s=(xtrk(j2,ns)-xtrk(j2,ns-1))/deltf
            b=xtrk(j2,ns)-s*log10(ttrkk(j2,ns))
            x3=s*log10(y)+b
          endif
        else
          deltf=ttrkk(j2,ns)-ttrkk(j2,ns-1)
          s=(xtrk(j2,ns)-xtrk(j2,ns-1))/deltf
          b=xtrk(j2,ns)-s*ttrkk(j2,ns)
          x3=s*y+b
        end if
        case2=1
        goto 8
      end if
      
C     ---   Times between the minimum and maximum ---
      do 7 k=1,ntrk(j2)-1
        if (y .ge. ttrkk(j2,k) .and. y .lt. ttrkk(j2,k+1)) then
          y3=ttrkk(j2,k  )
          y4=ttrkk(j2,k+1)
          x3=xtrk(j2,k  )
          x4=xtrk(j2,k+1)
          case2=0      
          goto 8
        end if
7     continue

C     ---   Bilinear interpolation ---
8     m1=mtrk(j1)
      m2=mtrk(j2)
      den=(m-m1)/(m2-m1)
      if(case1.eq.0.AND.case2.eq.0) then
        ym1=y1+(y3-y1)*den
        ym2=y2+(y4-y2)*den
        xm1=x1+(x3-x1)*den
        xm2=x2+(x4-x2)*den
        xout=xm1+(y-ym1)/(ym2-ym1)*(xm2-xm1)
      elseif(case1.eq.0.AND.case2.eq.1) then
        xm1=x1+(x2-x1)*(y-y1)/(y2-y1)
        xout=xm1+(x3-xm1)*den
      elseif(case1.eq.1.AND.case2.eq.0) then
        xm2=x3+(x4-x3)*(y-y3)/(y4-y3)
        xout=x1+(xm2-x1)*den
      else
        xout=x1+(x3-x1)*den
      end if
100   continue

      if(isNaN(xout)) then
        write(*,*) "Error interp(tcool,mass,model,log):",y,m,model,
     &             modlog
        stop
      end if

      return
      end
C***********************************************************************
