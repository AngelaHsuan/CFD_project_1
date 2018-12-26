! Written by Ting-Hsuan Hsu

! Start date 10 13 2017
! End date 10 25 2017

! Trying to use numerical method to solve the poisson eqqation in the channel flow
! Which only suffers from the pressure gradient
! A for the cross section aspect ratio, w for relaxation coefficient, iteration method is the initail condition to give

!------------------------------------------------------
program Main
    implicit NONE
    
    integer :: n,Ny,Nz,Pad,j,k
    real :: delta_y,delta_z,A,w,conv,Q
    real*8 :: Time_Exe 
    integer Time_Start,Time_End,Rate
    real*8, allocatable, dimension(:,:) :: u_old,u_new
    character :: SOR,printoutPLA
    
    write(*,*) 'Please enter the amount of A'
    write(*,*) 'for the rectengular channel of various aspect ratio.'
    write(*,*) 'A ='
    read(*,*)A
    
    if (A <= 0) then
        write(*,*)' A need to be positive!!!!'
        stop
    end if
    
    delta_y = 0.01
    delta_z = 0.01
    Nz = 2/delta_z+1       !zよ翴计
    Ny = (A/delta_y)+1	   !yよ翴计
    n = Ny*Nz              !羆翴计
    write(*,*)'Enter the relaxation coefficient'!肞Β玒计
    read(*,*)w
    conv = 1D-6     !Μ滥非玥
    
    allocate (u_old(Ny,Nz), u_new(Ny,Nz))
    
    write(*,*)'What method do you want to use?'
    write(*,*)'PointSOR type "P", LineSOR type "L", ADI type "A"'
    read(*,*)SOR
    
    if (SOR == 'P') then
        call system_clock(Time_Start,Rate)
        call PointSOR1(Ny,Nz,delta_y,delta_z,u_old,u_new,w,conv)
        printoutPLA = 'P'
        
        
    else if (SOR == 'L') then
        call system_clock(Time_Start,Rate)
        call LineSORrow(Ny,Nz,delta_y,delta_z,u_old,u_new,w,conv)
        printoutPLA = 'L'
        
        
    else if (SOR == 'A') then
        call system_clock(Time_Start,Rate)
        call ADI(Ny,Nz,delta_y,delta_z,u_old,u_new,w,conv)
        printoutPLA = 'A'
        
        
    else
        write(*,*)'Please enter the correct command'
        stop
                
    end if
    call printout(delta_y,delta_z,Ny,Nz,u_new,printoutPLA)
    call SYSTEM_CLOCK(Time_End,Rate)
    Time_Exe = (real(Time_End)-real(Time_Start))/real(Rate)
    write(*,*) '------------------------------------------------------------'
    Write(*,*)'The computational time(sec):', Time_Exe
    write(*,*) '------------------------------------------------------------'
    
! Computation of volumetric flow rate--------------------------
    Q = 0
    ! Center
    do k = 2, Nz-1, 1
        do j = 2, Ny-1, 1
            Q = Q + u_new(j,k)*delta_y*delta_z
        end do
    end do
    ! Bottom and Top
    do j = 2, Ny-1, 1
        Q = Q + u_new(j,1)*delta_y*delta_z/2 + u_new(j,Nz)*delta_y*delta_z/2
    end do
    ! Left and right
    do k = 2, Nz-1, 1
        Q = Q + u_new(1,k)*delta_y*delta_z/2 + u_new(Ny,k)*delta_y*delta_z/2
    end do
    ! Corner
    Q = Q + (u_new(1,1) + u_new(1,Nz) + u_new(Ny,1) + u_new(Ny,Nz))*delta_y*delta_z/4
    
    write(*,*) "Volumetric flow rate = ", Q
           
    
    write(*,*) '--------------------------Done------------------------------'
    
end program Main