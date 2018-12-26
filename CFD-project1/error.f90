subroutine error(width,height,u_errMAX)
    implicit NONE
    
    integer :: width,height,m,n
    real :: err,u_errMAX
    real*8, dimension(width,height) :: u_old,u_new,u_err
        
    u_errMAX = 0
    do m = 1,width,1
        do n = 1,height,1
            if (u_errMAX < u_err(m,n)) then
            u_errMAX = u_err(m,n)
            end if
        end do
    end do
    
end subroutine error