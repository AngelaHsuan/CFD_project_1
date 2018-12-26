!------------------------------------------------------------
! Result output
subroutine output(del_x,del_y,N,CVD,printout,J_x,divide_1, divide_2, J_avg_x)
    implicit NONE
    
    integer :: i,j,N,Pad,printout
    real*8 :: del_x,del_y,divide_1, divide_2
    real*8, dimension(N) :: J_x, J_avg_x
    real*8, dimension(N,N) :: CVD
    
    Pad = 100
    ! Output CVD
    if (printout == 1) then
        open (unit=Pad, file="FOU.txt", status="UNKNOWN")
    else if (printout == 2) then
        open (unit=Pad, file="SOCD.txt", status="UNKNOWN")
    else if (printout == 3) then
        open (unit=Pad, file="QUICK.txt", status="UNKNOWN")
    else if (printout == 4) then
        open (unit=Pad, file="TOU.txt", status="UNKNOWN")
    else if (printout == 5) then
        open (unit=Pad, file="SOU.txt", status="UNKNOWN")
    else
        write(*,*)'there is error in printout'
        pause
    end if
        
    write(Pad,*) 'variables="x","y","CVD"'
    write(Pad,*) 'zone i=', N,'j=',N,'DATAPACKING=POINT'
    
    do j = 1,N,1
        do i = 1,N,1
            write(Pad,'(3F20.12)') 0.d0+(i-1)*del_x, 0.d0+(j-1)*del_y, CVD(i,j)
        end do
    end do
    close (Pad,status = 'Keep')
    
    ! Output J
    open (unit=Pad, file="deposition rate.txt", status="UNKNOWN")
    
    write(Pad,*) 'variables="x","J"'
    write(Pad,*) 'zone i=', divide_2-divide_1+1 ,'DATAPACKING=POINT'
    
    do i = divide_1, divide_2,1
        write(Pad,'(2F20.12)') 0.d0+i*del_x, J_x(i)
    end do
    close (Pad,status = 'Keep')
    
    ! Output J_avg
    open (unit=Pad, file="average deposition rate.txt", status="UNKNOWN")
    
    write(Pad,*) 'variables="x","J_bar"'
    write(Pad,*) 'zone i=', divide_2-divide_1+1 ,'DATAPACKING=POINT'
    
    do i = divide_1, divide_2,1
        write(Pad,'(2F20.12)') 0.d0+i*del_x, J_avg_x(i)
    end do
    close (Pad,status = 'Keep')
    
end subroutine output