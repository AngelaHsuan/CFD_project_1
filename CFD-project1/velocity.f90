!------------------------------------------------
subroutine velocity(N,Nv,u,v)
    implicit NONE
    
    integer :: i, j, N, Nv, Pad
    real*8 :: pi, del_x, del_y
    real*8, dimension(Nv,Nv) :: u, v
    
    del_x = 1.d0/(Nv-1)
    del_y = 1.d0/(Nv-1)
    
    pi = 4.d0*atan(1.d0)
    
    do i = 1,Nv,1
        do j = 1, Nv, 1
            u(i,j) = sin(pi*(i-1)*del_x)*cos(pi*(j-1)*del_y)
            v(i,j) = -cos(pi*(i-1)*del_x)*sin(pi*(j-1)*del_y)
        end do
    end do
    
    Pad = 100
    open(unit = Pad, file = "velocity.txt", status = "UNKNOWN")
    write(Pad, *)'variables="x","y","u","v"'
    write(Pad,*) 'zone i=', Nv,'j=',Nv,'DATAPACKING=POINT'
    
    do j = 1,Nv,1
        do i = 1,Nv,1
            write(Pad,'(4F20.12)') 0.d0+(i-1)*del_x, 0.d0+(j-1)*del_y, u(i,j), v(i,j)
        end do
    end do
    close (Pad,status = 'Keep')
    
end subroutine velocity