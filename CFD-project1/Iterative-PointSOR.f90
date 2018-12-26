subroutine PointSOR(C_p,u_p,b,u_old,u_new,n,A,delta_y,delta_z,j,k)
    implicit NONE
    
    integer :: n,A,j,k,i
    real :: delta_y,delta_z
    real*8, dimension(n,n) :: C_p,u_old,u_new
    real*8, dimension(n) :: u_p,b
    
    u_old = 0
    
    !---Boundary condition------
    
    b(1:((1/0.01)+1)) = (-1/2)-(1/(delta_y)**2)*u_old(j+1,k)-(1/(2*(delta_z)**2))*u_old(j,k+1)-(1/(2*(delta_z)**2))*u_old(j,k-1)
    
    do i=1:101
        b(101*i) = 0
        b(101*i+1) = (-1/2)-(1/(2*(delta_y)**2))*u_old(j+1,k)-(1/(2*(delta_z)**2))*u_old(j-1,k)-(1/(delta_z)**2)*u_old(j,k-1)
        b((n-101):n) = 0
    end do
        
    b()
    C_p() = ((-1)/(delta_z)**2)-(1/(delta_y)**2)
    
end subroutine PointSOR