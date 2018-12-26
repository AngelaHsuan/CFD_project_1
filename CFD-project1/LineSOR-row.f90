!----------------------------------------------------------------
! Line SOR of row
subroutine LineSORrow(N, del_x, del_y, conv, Pe_m, w, beta, G,CVD,divide_1, divide_2,CVD_errMax, DCM)
    implicit NONE
    
    integer :: i,j, N, Nv, N_iter,G, DCM
    real*8 :: del_x, del_y, u_errMAX, divide_1, divide_2, conv, CVD_errMAX, Pe_m,w,beta
    real*8, dimension(N) :: CVD_j,C_w_FOU,C_p_FOU,C_e_FOU,RHS,C_w_HO,C_p_HO,C_e_HO
    real*8, dimension(N,N) :: CVD, CVD_new, CVD_err,Res
    real*8, allocatable, dimension(:,:) :: u,v
    character(len=20) :: position
    
    Nv = N*2-1
    allocate(u(Nv,Nv), v(Nv,Nv))
    call velocity(N,Nv,u,v)
    
    Res = 1.d0
    CVD = 0.d0
    CVD_errMAX = 1.d0
    N_iter = 0
    divide_1 = (N-1)/5
    divide_2 = (N-1)/5*4
    
    do while (CVD_errMAX > conv)
        do j = 1,N,1
            if (j == 1) then
                do i = 1,N,1
                    if (i == 1) then    !(0,0)
                        C_w_FOU(i) = 0
                        C_p_FOU(i) = 1
                        C_e_FOU(i) = 0
                        RHS(i) = 1
                        !CVD(i,j) = 1
                    else if (i == N) then   !right-bottom
                        position = 'Right_bottom'
                        call coefficient(position,i,j,N,Pe_m,del_x,del_y, C_w_FOU, C_p_FOU, C_e_FOU,beta,G,RHS,CVD,u,v,Nv,C_w_HO,C_p_HO,C_e_HO, DCM)
                    else    !bottom
                        position = 'Bottom'
                        call coefficient(position,i,j,N,Pe_m,del_x,del_y, C_w_FOU, C_p_FOU, C_e_FOU,beta,G,RHS,CVD,u,v,Nv,C_w_HO,C_p_HO,C_e_HO, DCM)
                    end if
                end do
            else if (j == 2) then
                do i = 1,N,1
                    if (i == 1) then    !(1,0)
                        C_w_FOU(i) = 0
                        C_p_FOU(i) = 1
                        C_e_FOU(i) = 0
                        RHS(i) = 1
                        !CVD(i,j) = 1
                    else if (i == N) then   !right
                        position = 'Right'
                        call coefficient(position,i,j,N,Pe_m,del_x,del_y, C_w_FOU, C_p_FOU, C_e_FOU,beta,G,RHS,CVD,u,v,Nv,C_w_HO,C_p_HO,C_e_HO, DCM)
                    else    !bottom
                        position = 'Boundary_center'
                        call coefficient(position,i,j,N,Pe_m,del_x,del_y, C_w_FOU, C_p_FOU, C_e_FOU,beta,G,RHS,CVD,u,v,Nv,C_w_HO,C_p_HO,C_e_HO, DCM)
                    end if
                end do
            else if (j == N) then
                do i = 1,N,1
                    if (i == 1) then
                        C_w_FOU(i) = 0
                        C_p_FOU(i) = 1
                        C_e_FOU(i) = 0
                        RHS(i) = 1
                        !CVD(i,j) = 1
                    else if (i == N) then   !Right_top_coner
                        position = 'Right_top'
                        call coefficient(position,i,j,N,Pe_m,del_x,del_y, C_w_FOU, C_p_FOU, C_e_FOU,beta,G,RHS,CVD,u,v,Nv,C_w_HO,C_p_HO,C_e_HO, DCM)
                    else if (i < divide_1) then    !Left_top
                        position = 'Top'
                        call coefficient(position,i,j,N,Pe_m,del_x,del_y, C_w_FOU, C_p_FOU, C_e_FOU,beta,G,RHS,CVD,u,v,Nv,C_w_HO,C_p_HO,C_e_HO, DCM)
                    else if (i > divide_2) then   !Right_top
                        position = 'Top'
                        call coefficient(position,i,j,N,Pe_m,del_x,del_y, C_w_FOU, C_p_FOU, C_e_FOU,beta,G,RHS,CVD,u,v,Nv,C_w_HO,C_p_HO,C_e_HO, DCM)
                    else    !top
                        C_w_FOU(i) = 0
                        C_p_FOU(i) = 1
                        C_e_FOU(i) = 0
                        RHS(i) = 0
                        !CVD(i,j) = 0
                    end if
                end do
            else if (j == N-1) then
                do i = 1,N,1
                    if (i == 1) then    !(1,0)
                        C_w_FOU(i) = 0
                        C_p_FOU(i) = 1
                        C_e_FOU(i) = 0
                        RHS(i) = 1
                        !CVD(i,j) = 1
                    else if (i == N) then   !right
                        position = 'Right'
                        call coefficient(position,i,j,N,Pe_m,del_x,del_y, C_w_FOU, C_p_FOU, C_e_FOU,beta,G,RHS,CVD,u,v,Nv,C_w_HO,C_p_HO,C_e_HO, DCM)
                    else    !top_center
                        position = 'Boundary_center'
                        call coefficient(position,i,j,N,Pe_m,del_x,del_y, C_w_FOU, C_p_FOU, C_e_FOU,beta,G,RHS,CVD,u,v,Nv,C_w_HO,C_p_HO,C_e_HO, DCM)
                    end if                 
                end do
            else
                do i = 1,N,1
                    if (i == 1) then    !Left
                        C_w_FOU(i) = 0
                        C_p_FOU(i) = 1
                        C_e_FOU(i) = 0
                        RHS(i) = 1
                        !CVD(i,j) = 1
                    else if (i == 2) then   !Left-center
                        position = 'Boundary_center'    
                        call coefficient(position,i,j,N,Pe_m,del_x,del_y, C_w_FOU, C_p_FOU, C_e_FOU,beta,G,RHS,CVD,u,v,Nv,C_w_HO,C_p_HO,C_e_HO, DCM)
                    else if (i == N) then   !Right
                        position = 'Right'
                        call coefficient(position,i,j,N,Pe_m,del_x,del_y, C_w_FOU, C_p_FOU, C_e_FOU,beta,G,RHS,CVD,u,v,Nv,C_w_HO,C_p_HO,C_e_HO, DCM)
                    else if (i == N-1) then !Right-center
                        position = 'Boundary_center'
                        call coefficient(position,i,j,N,Pe_m,del_x,del_y, C_w_FOU, C_p_FOU, C_e_FOU,beta,G,RHS,CVD,u,v,Nv,C_w_HO,C_p_HO,C_e_HO, DCM)
                    else    !Center
                        position = 'Center'
                        call coefficient(position,i,j,N,Pe_m,del_x,del_y, C_w_FOU, C_p_FOU, C_e_FOU,beta,G,RHS,CVD,u,v,Nv,C_w_HO,C_p_HO,C_e_HO, DCM)
                        if (DCM == 1) then
                        else if (DCM == 0) then
                            C_w_FOU(i) = C_w_HO(i)
                            C_p_FOU(i) = C_p_HO(i)
                            C_e_FOU(i) = C_e_HO(i)
                        else
                            write(*,*)'Mistake of DCM in LSOR'
                        end if
                    end if
                end do
            end if
            call TDMA(N,C_w_FOU,C_p_FOU,C_e_FOU,RHS,CVD_j)
            do i = 1,N,1
                CVD_new(i,j) = CVD_j(i)
                CVD_err(i,j) = abs(CVD_new(i,j)-CVD(i,j))
                if (i == 1) then
                    Res(i,j) = C_p_FOU(i)*CVD_new(i,j) + C_e_FOU(i)*CVD_new(i+1,j) - RHS(i)
                else if (i == N) then
                    Res(i,j) = C_w_FOU(i)*CVD_new(i-1,j) + C_p_FOU(i)*CVD_new(i,j) - RHS(i)
                else
                    Res(i,j) = C_w_FOU(i)*CVD_new(i-1,j) + C_p_FOU(i)*CVD_new(i,j) + C_e_FOU(i)*CVD_new(i+1,j) - RHS(i)
                end if
                CVD(i,j) = w*CVD_new(i,j) + (1-w)*CVD(i,j)
            end do
        end do
        
        CVD_errMax = maxval(CVD_err)
        write(*,'(A22,I7,A11,1E12.5,A13,1E12.5)') 'number of iteration :', N_iter,'error=', CVD_errMax,'Residue=',maxval(Res)
        N_iter = N_iter + 1
        
        if (CVD_errMax>1D10) then
            write(*,*)'---------------------------------------------------------------'
            write(*,*)'The result divereged!!!!!!!!!!!'
            CVD_errMax = 0
        end if
    end do
    !write(*,'(A22,I7,A11,1E12.5)')'number of iteration :', N_iter,'error=', CVD_errMax
end subroutine LineSORrow