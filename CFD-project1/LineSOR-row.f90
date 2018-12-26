!----------------------------------------------------------------
! Line SOR of row
subroutine LineSORrow(Ny,Nz,delta_y,delta_z,u_old,u_new,w,conv)
    implicit NONE
    
    integer :: j,k,m,n,Ny,Nz,N_iter,i
    real :: delta_y,delta_z,err,u_errMAX,w,conv
    real*8, dimension(Ny) :: b,C_p,C_N,C_S,C_E,C_W,u_fnew
    real*8, dimension(Ny,Nz) :: u_old,u_new,u_err
    
    u_old = 0.d0
    u_errMAX = 1.d0
    N_iter = 0.d0
    
    do while (u_errMAX > conv)
        do k = 1,Nz,1
        
            if (k == 1) then
                do i = 1,Ny,1
                    u_new(i,k) = 0.d0
                    
                end do
                
            else if (k == Nz) then
                do j = 1,Ny,1
                    if (j == 1) then        !Corner
                        
                        C_P(j) = -1/delta_z**2 - 1/delta_y**2
                        C_E(j) = 1/delta_y**2
                        C_S(j) = -1/delta_z**2
                        b(j) = -1/2 + C_S(j)*u_old(j,k-1)
                    
                    else if (j == Ny) then  !¥k¤W
                        C_W(j) = 0
                        C_P(j) = 1
                        b(j) = 0
                        u_new(j,k) = 0
                    else                    !Top
                        C_W(j) = 1/delta_y**2
                        C_P(j) = -2/delta_z**2 - 2/delta_y**2
                        C_E(j) = 1/delta_y**2
                        C_S(j) = -2/delta_z**2
                        b(j) = -1+C_S(j)*u_old(j,k-1)
                    
                    end if
                
                end do
                call TDMA(Ny,C_W,C_P,C_E,b,u_fnew)
                do i = 1,Ny,1
                    u_new(i,k) = u_fnew(i)
                    u_err(i,k) = abs(u_new(i,k)-u_old(i,k))
                    u_old(i,k) = w*u_new(i,k) + (1-w)*u_old(i,k)
                end do
                
            else
                do j = 1,Ny,1
                    if (j == 1) then        !Left
                        
                        C_P(j) = -2/delta_z**2 - 2/delta_y**2
                        C_E(j) = 2/delta_y**2
                        C_N(j) = -1/delta_z**2
                        C_S(j) = -1/delta_z**2
                        b(j) = -1 + C_S(j)*u_old(j,k-1) + C_N(j)*u_old(j,k+1)
                    
                    else if (j == Ny) then  !Right
                        C_W(j) = 0
                        C_P(j) = 1
                        b(j) = 0
                        u_new(j,k) = 0
                    else                    !Center
                        C_W(j) = 1/delta_y**2
                        C_P(j) = -2/delta_y**2 - 2/delta_z**2
                        C_E(j) = 1/delta_y**2
                        C_N(j) = -1/delta_z**2
                        C_S(j) = -1/delta_z**2
                        b(j) = -1 + C_S(j)*u_old(j,k-1) + C_N(j)*u_old(j,k+1)
                        
                    end if
                
                end do
                call TDMA(Ny,C_W,C_P,C_E,b,u_fnew)
                do i = 1,Ny,1
                    u_new(i,k) = u_fnew(i)
                    u_err(i,k) = abs(u_new(i,k)-u_old(i,k))
                    u_old(i,k) = w*u_new(i,k) + (1-w)*u_old(i,k)
                end do
                
            end if
        
        end do
        
        
        
    u_errMax = maxval(u_err)
    
    write(*,*) 'number of iteration', N_iter,'    error    ', u_errMax
    
    N_iter = N_iter + 1
          
    end do
    
    
end subroutine LineSORrow