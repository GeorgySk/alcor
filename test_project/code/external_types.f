C     TODO: delete unncecessary attributes
      MODULE external_types

      TYPE FileGroupInfo
        SEQUENCE
        CHARACTER(len=7) :: tableType
        CHARACTER(len=3) :: WDtype
        INTEGER :: flag
        REAL :: Z
        INTEGER :: initLink
        INTEGER :: link
        INTEGER :: ncol
        INTEGER :: nrow
C       should be ntrk(ncol) but compiler gets angry. 10 is max of ncol
C        INTEGER :: ntrk(10)
        INTEGER,allocatable :: ntrk(:)
C       should be mass(ncol)..        
c        DOUBLE PRECISION :: mass(10)
        DOUBLE PRECISION,allocatable :: mass(:)
C       should be coolingTime(ncol,nrow)..   
C        DOUBLE PRECISION :: coolingTime(10,650)
        DOUBLE PRECISION,allocatable :: coolingTime(:,:)
C       should be prevTime(ncol)..        
C        DOUBLE PRECISION :: prevTime(10)
        DOUBLE PRECISION,allocatable :: prevTime(:)
C       should be luminosity(ncol,nrow)..   
C        DOUBLE PRECISION :: luminosity(10,650)
        DOUBLE PRECISION,allocatable :: luminosity(:,:)
C       should be effTemp(ncol,nrow)..   
C        DOUBLE PRECISION :: effTemp(10,650)
        DOUBLE PRECISION ,allocatable:: effTemp(:,:)
C       should be gravAcc(ncol,nrow)..   
C        DOUBLE PRECISION :: gravAcc(10,650)
        DOUBLE PRECISION,allocatable :: gravAcc(:,:)
C       should be color_U(ncol,nrow)..   
C        DOUBLE PRECISION :: color_U(10,650)
        DOUBLE PRECISION,allocatable :: color_U(:,:)
C       should be color_B(ncol,nrow)..   
C        DOUBLE PRECISION :: color_B(10,650)
        DOUBLE PRECISION,allocatable :: color_B(:,:)
C       should be color_V(ncol,nrow)..   
C        DOUBLE PRECISION :: color_V(10,650)
        DOUBLE PRECISION,allocatable :: color_V(:,:)
C       should be color_R(ncol,nrow)..   
C        DOUBLE PRECISION :: color_R(10,650)
        DOUBLE PRECISION,allocatable :: color_R(:,:)
C       should be color_I(ncol,nrow)..   
C        DOUBLE PRECISION :: color_I(10,650)
        DOUBLE PRECISION,allocatable :: color_I(:,:)
      ENDTYPE FileGroupInfo

      END MODULE external_types