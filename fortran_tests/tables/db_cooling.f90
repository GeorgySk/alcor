subroutine incooldb(flag, initLink, ncol, nrow, ntrk, mass, coolingTime, &
                    prevTime, luminosity, effTemp, gravAcc)
    implicit none
    
    integer :: i, j
    real :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, &
            a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, a26
    integer, intent(in) :: flag
    integer, intent(in) :: initLink
    integer, intent(in) :: ncol
    integer, intent(in) :: nrow
    integer, intent(inout) :: ntrk(ncol)
    real, intent(inout) :: mass(ncol)
    real, intent(inout) :: coolingTime(ncol, nrow)
    real, intent(inout) :: prevTime(ncol)
    real, intent(inout) :: luminosity(ncol, nrow)
    real, intent(inout) :: effTemp(ncol, nrow)
    real, intent(inout) :: gravAcc(ncol, nrow)
    
    if (flag == 1) then
        mass(1) = 0.5047
        mass(2) = 0.5527
        mass(3) = 0.59328
        mass(4) = 0.62738
        mass(5) = 0.6602
        mass(6) = 0.69289
        mass(7) = 0.8637
    else if (flag == 2) then
        mass(1) = 0.514
        mass(2) = 0.53
        mass(3) = 0.542
        mass(4) = 0.565
        mass(5) = 0.584
        mass(6) = 0.61
        mass(7) = 0.664
        mass(8) = 0.741
        mass(9) = 0.869
    else
        mass(1) = 0.524
        mass(2) = 0.570
        mass(3) = 0.593
        mass(4) = 0.61
        mass(5) = 0.632
        mass(6) = 0.659
        mass(7) = 0.70
        mass(8) = 0.76
        mass(9) = 0.87
    end if
    
    if (flag == 1) then
        prevTime(1) = 0.0
        prevTime(2) = 0.0
        prevTime(3) = 0.0
        prevTime(4) = 0.0
        prevTime(5) = 0.0
        prevTime(6) = 0.0
        prevTime(7) = 0.0
    else
        prevTime(1) = 11.117
        prevTime(2) = 2.7004
        prevTime(3) = 1.699
        prevTime(4) = 1.2114
        prevTime(5) = 0.9892
        prevTime(6) = 0.7422
        prevTime(7) = 0.4431
        prevTime(8) = 0.287
        prevTime(9) = 0.114
    end if
     
    if (flag == 1 .or. flag == 3) then
        do i = 1, ncol
            do j = 1, nrow
                read(initLink + i, *, end=2) a1, a2, a3, a4, a5, &
                        a6, a7, a8, a9, &
                        a10, a11, a12, a13, &
                        a14, a15, a16, a17, &
                        a18, a19, a20, a21, &
                        a22, a23, a24, a25, &
                        a26
                coolingTime(i, j) = 10.0 ** a9 / 1000.0
                coolingTime(i, j) = coolingTime(i, j) - prevTime(i)
                effTemp(i, j) = 10.0 ** a2
                gravAcc(i, j) = a23
                luminosity(i, j) = a1
            end do
2                ntrk(i) = j - 1
        end do
    else
        do i = 1, ncol
            do j = 1, nrow
                read (initLink + i, *, end=20) a1, a2, a3, a4
                coolingTime(i, j) = a3
                effTemp(i, j) = 10.0 ** a2
                gravAcc(i, j) = a4
                luminosity(i, j) = a1
            end do
20              ntrk(i) = j - 1
        end do
    end if
end subroutine
