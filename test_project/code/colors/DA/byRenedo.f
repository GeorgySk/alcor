      subroutine color(table)
C=======================================================================
C     This subroutine reads the colors by Rene and interpolates 
C     according to the vector of reference time
C-----------------------------------------------------------------------
C     I/O parameters:
C       table: instance of FileGroupInfo
C=======================================================================
      use external_types
      implicit double precision (a-h,m,o-z)

      double precision vtanmin
      parameter (vtanmin=30)

      TYPE(FileGroupInfo) :: table
      integer i,k
      double precision a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14
      double precision a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26
      double precision a27,a28
C     TODO: size should be table%nrow      
      double precision Hg(650),ggii(650),Hbj(650),
     &                 bjr(650)
      double precision g,gr,bj,term

C     Curve to delimit RPMD, for Mwd=0.61Mo
C     NOTE: these common blocks are never used      
      common /Hgcurve/ Hg,ggii
      common /Hbjcurve/ Hbj,bjr

C     Read masses
      table%mass(1)=0.524
      table%mass(2)=0.570
      table%mass(3)=0.593
      table%mass(4)=0.610
      table%mass(5)=0.632
      table%mass(6)=0.659
      table%mass(7)=0.705
      table%mass(8)=0.767
      table%mass(9)=0.837
      table%mass(10)=0.877
      
C     read values from files
      do k=1,table%ncol
        do i=1,table%nrow
          read(table%initLink+k,*,end=15) a1,a2,a3,a4,a5,a6,a7,a8,a9,
     &        a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,
     &        a24,a25,a26,a27,a28
          table%luminosity(k,i)=a3
          table%color_U(k,i)=a24
          table%color_B(k,i)=a25
          table%color_V(k,i)=a26
          table%color_R(k,i)=a27
          table%color_I(k,i)=a28
        enddo
15      table%ntrk(k)=i-1
      end do
      
      term=5*log10(vtanmin)-3.379
      do i=1,table%ntrk(4)
C       SDSS
        g = table%color_V(4,i)+0.630*(table%color_B(4,i)-
     &      table%color_V(4,i))-0.124
C       Hg=Mg+5*log10(v_tan)-3.379
        Hg(i)=g+term
        ggii(i)=1.646*(table%color_V(4,i)-table%color_R(4,i))+
     &          1.007*(table%color_R(4,i)-table%color_I(4,i))-0.375
C       SuperCosmos
        gr=1.646*(table%color_V(4,i)-table%color_R(4,i))-0.139
        bj=0.15+0.13*gr+g
C       Hbj=Mbj+5*log10(v_tan)-3.379
        Hbj(i)=bj+term
        bjr(i)=0.28+1.13*gr
      end do
      
      return
      end
C***********************************************************************