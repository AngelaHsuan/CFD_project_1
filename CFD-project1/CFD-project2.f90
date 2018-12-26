! Written by Ting-Hsuan Hsu

! Start date 11 25 2017
! End date 11 22 2017

! Trying to use numerical method to solve the convection-diffusion mass transfer problem
!------------------------------------------------------------------------------------------------------------------------------
program Main
    implicit NONE
    
    integer :: N,G,printout, DCM
    real*8 :: Pe_m, del_x, del_y,conv,w,beta,divide_1, divide_2,Time_Exe,i,k,CVD_errMax
    integer Time_Start,Time_End,Rate
    real*8, allocatable, dimension(:) :: J_x, J_avg_x
    real*8, allocatable, dimension(:,:) :: CVD
    character :: again
    
    again = 'y'
    do while (again == 'y')
        call Initial(N, del_x, del_y, Pe_m, conv, w,beta,G,printout, DCM)
        allocate (CVD(N,N))
        call system_clock(Time_Start,Rate)
        call LineSORrow(N, del_x, del_y, conv, Pe_m, w, beta, G,CVD,divide_1, divide_2,CVD_errMax, DCM)
        call system_clock(Time_End,Rate)
        Time_Exe = (real(Time_End)-real(Time_Start))/real(Rate)
        write(*,*) '------------------------------------------------------------'
        Write(*,*)'The computational time(sec):', Time_Exe
        write(*,*) '------------------------------------------------------------'
        allocate (J_x(N), J_avg_x(N))
        if (CVD_errMax > 0) then
            call deposition_rate(N,divide_1, divide_2,del_x, del_y, CVD,J_x, J_avg_x)
            call output(del_x,del_y,N,CVD,printout,J_x,divide_1, divide_2, J_avg_x)
        end if
        write(*,*)'Do you want to do it again?(y/n)'
        read(*,*)again
        deallocate(CVD,J_x,J_avg_x)
    end do
    
end program Main