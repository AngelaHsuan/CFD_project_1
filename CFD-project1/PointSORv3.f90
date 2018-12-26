!-------------------------------------------------------
! Point SOR
subroutine PointSOR1(Ny,Nz,delta_y,delta_z,u_old,u_new,w,conv)
    implicit NONE
    
    integer :: j,k,Ny,Nz,N_iter
    real :: delta_y,delta_z,C_p,b,u_errMAX,C_N,C_S,C_E,C_W,w,conv
    real*8, dimension(Ny,Nz) :: u_old,u_new,u_err
    
    u_old = 0
    u_errMAX = 1.d0
    N_iter = 0
    
    do while (u_errMAX > conv)
        do j = 1,Ny,1
            do k = 1,Nz,1
                if (j == 1 .and. k == Nz) then  !Coner
                
				    C_E = (1/(delta_y)**2)
                    C_W = 0
                    C_N = 0
                    C_S = (1/(delta_z)**2)
                    
                    b = -1/2 - C_E*u_old(j+1,k) - C_S*u_old(j,k-1)
                    C_p = -1/(delta_z)**2 - 1/(delta_y)**2
                    u_new(j,k) = b / C_p
                else if (j == 1 .and. k == 1) then 
                    u_new(j,k) = 0
                else if (j == Ny .and. k == 1) then
                    u_new(j,k) = 0
                else if (j == Ny .and. k == Nz) then
                    u_new(j,k) = 0
                else if (j == 1) then
                
                    C_E = (2/(delta_y)**2)
                    C_W = 0
                    C_N = (1/(delta_z)**2)
                    C_S = (1/(delta_z)**2)
                    
                    b = -1 - C_E*u_old(j+1,k) - C_N*u_old(j,k+1) - C_S*u_old(j,k-1)
                    C_p = -2/(delta_z)**2 - 2/(delta_y)**2
                    u_new(j,k) = b / C_p
                else if (k == 1) then
                    u_new(j,k) = 0
                else if (k == Nz) then
                
                    C_E = (1/(delta_y)**2)
                    C_W = (1/(delta_y)**2)
                    C_N = 0
                    C_S = (2/(delta_z)**2)
                    
                    b = -1 - C_E*u_old(j+1,k) - C_W*u_old(j-1,k) - C_S*u_old(j,k-1)
                    C_p = -2/(delta_z)**2 - 2/(delta_y)**2
                    u_new(j,k) = b / C_p
                else if (j == Ny) then
                    u_new(j,k) = 0
                else
                
                    C_E = (1/(delta_y)**2)
                    C_W = (1/(delta_y)**2)
                    C_N = (1/(delta_z)**2)
                    C_S = (1/(delta_z)**2)
                    
                    b = -1 - C_E*u_old(j+1,k) - C_W*u_old(j-1,k) - C_N*u_old(j,k+1) - C_S*u_old(j,k-1)
                    C_p = -2/(delta_z)**2 - 2/(delta_y)**2
                    u_new(j,k) = b / C_p
                end if
                u_err(j,k) = abs(u_new(j,k) - u_old(j,k))
                u_old(j,k) = w*u_new(j,k) + (1-w)*u_old(j,k)
                
            end do                
        end do
    
    
    u_errMAX = 0
    do j = 1,Ny,1
        do k = 1,Nz,1
            if (u_errMAX < u_err(j,k)) then
            u_errMAX = u_err(j,k)
            end if
        end do
    end do
    write(*,*) 'number of iteration', N_iter,'    error    ', u_errMax 
    
    N_iter = N_iter + 1
    
    end do
    
end subroutine PointSOR1