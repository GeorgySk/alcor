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
c        real :: mass(10)
        real,allocatable :: mass(:)
C       should be coolingTime(ncol,nrow)..   
C        real :: coolingTime(10,650)
        real,allocatable :: coolingTime(:,:)
C       should be prevTime(ncol)..        
C        real :: prevTime(10)
        real,allocatable :: prevTime(:)
C       should be luminosity(ncol,nrow)..   
C        real :: luminosity(10,650)
        real,allocatable :: luminosity(:,:)
C       should be effTemp(ncol,nrow)..   
C        real :: effTemp(10,650)
        real ,allocatable:: effTemp(:,:)
C       should be gravAcc(ncol,nrow)..   
C        real :: gravAcc(10,650)
        real,allocatable :: gravAcc(:,:)
C       should be color_U(ncol,nrow)..   
C        real :: color_U(10,650)
        real,allocatable :: color_U(:,:)
C       should be color_B(ncol,nrow)..   
C        real :: color_B(10,650)
        real,allocatable :: color_B(:,:)
C       should be color_V(ncol,nrow)..   
C        real :: color_V(10,650)
        real,allocatable :: color_V(:,:)
C       should be color_R(ncol,nrow)..   
C        real :: color_R(10,650)
        real,allocatable :: color_R(:,:)
C       should be color_I(ncol,nrow)..   
C        real :: color_I(10,650)
        real,allocatable :: color_I(:, :)
        real,allocatable :: color_J(:, :)
      ENDTYPE FileGroupInfo

      END MODULE external_types