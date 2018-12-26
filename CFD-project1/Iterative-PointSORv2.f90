subroutine PointSOR(Ny,Nz,delta_y,delta_z,C_p,u_old,u_new)
    implicit NONE
    
    integer :: j,k,m,n,Ny,Nz,N_iter
    real :: delta_y,delta_z,C_p,b,err,u_errMAX
    real*8, dimension(Ny,Nz) :: u_old,u_new,u_err
    
    u_old = 0
    u_errMAX = 1.d0
    N_iter = 0
    do while (u_errMAX > 1D-6)
        do j = 1,Ny,1
            do k = 1,Nz,1
                if (j == 1 .and. k == Nz) then
                    b = -1/2 - (1/(delta_y)**2)*u_old(j+1,k) - (1/(delta_z)**2)*u_old(j,k-1)
                    C_p = -1/(delta_z)**2 - 1/(delta_y)**2
                    u_new(j,k) = b / C_p
                else if (j == 1 .and. k == 1) then
                    u_new(j,k) = 0
                else if (j == Ny .and. k == 1) then
                    u_new(j,k) = 0
                else if (j == Ny .and. k == Nz) then
                    u_new(j,k) = 0
                else if (j == 1) then
                    b = -1 - (2/(delta_y)**2)*u_old(j+1,k) - (1/(delta_z)**2)*u_old(j,k+1) - (1/(delta_z)**2)*u_old(j,k-1)
                    C_p = -2/(delta_z)**2 - 2/(delta_y)**2
                    u_new(j,k) = b / C_p
                else if (k == 1) then
                    u_new(j,k) = 0
                else if (k == Nz) then
                    b = -1 - (1/(delta_y)**2)*u_old(j+1,k) - (1/(delta_y)**2)*u_old(j-1,k) - (2/(delta_z)**2)*u_old(j,k-1)
                    C_p = -2/(delta_z)**2 - 2/(delta_y)**2
                    u_new(j,k) = b / C_p
                else if (j == Ny) then
                    u_new(j,k) = 0
                else
                    b = -1 - (1/(delta_y)**2)*u_old(j+1,k) - (1/(delta_y)**2)*u_old(j-1,k) - (1/(delta_z)**2)*u_old(j,k+1)- (1/(delta_z)**2)*u_old(j,k-1)
                    C_p = -2/(delta_z)**2 - 2/(delta_y)**2
                    u_new(j,k) = b / C_p
                end if           
            end do                
        end do
        
        ! ­pºâ»~®t
        u_err = u_new - u_old
        u_errMAX = 0
        do m = 1,Ny,1
            do n = 1,Nz,1
                if (u_errMAX < u_err(m,n)) then
                u_errMAX = u_err(m,n)
                end if
            end do
        end do
        write(*,*) N_iter, u_errMax
        
        u_old = u_new
        N_iter = N_iter + 1
        
    end do
    

    
end subroutine PointSOR