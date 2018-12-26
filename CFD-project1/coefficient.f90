!----------------------------------------------------------------------------
subroutine coefficient(position,i,j,N,Pe_m,del_x,del_y, C_w_FOU, C_p_FOU, C_e_FOU,beta,G,RHS,CVD,u,v,Nv,C_w_HO,C_p_HO,C_e_HO, DCM)
    implicit NONE
    
    integer :: i,j,N,Nv,G,DCM
    real*8 :: Pe_m,del_x,del_y,beta, u_EI_p, u_EI_n, u_WI_p, u_WI_n, v_NI_p, v_NI_n, v_SI_p, v_SI_n, u_p, v_p
    real*8, dimension(N) ::C_w_FOU,C_p_FOU,C_e_FOU,C_s_FOU,C_n_FOU,C_w_HO,C_p_HO,C_e_HO,C_s_HO,C_n_HO&
                           &,C_ww_HO,C_ee_HO,C_ss_HO,C_nn_HO,RHS
    real*8, dimension(N,N) :: CVD
    real*8, dimension(Nv,Nv) :: u,v
    character(len=20) :: position
        
    call uv(N,i,j,u_EI_p, u_EI_n, u_WI_p, u_WI_n, v_NI_p, v_NI_n, v_SI_p, v_SI_n,u_p,v_p,u,v,Nv)
    
    C_w_HO(i) = 0
    C_p_HO(i) = 0
    C_e_HO(i) = 0
    C_s_HO(i) = 0
    C_n_HO(i) = 0
    C_ww_HO(i) = 0
    C_ee_HO(i) = 0
    C_ss_HO(i) = 0
    C_nn_HO(i) = 0
    
    if (position == 'Boundary_center') then
        C_w_FOU(i) = -1.d0/del_x*u_WI_p - 1.d0/Pe_m/del_x**2
        C_p_FOU(i) = 1.d0/del_x*(u_EI_p - u_WI_n) + 1.d0/del_y*(v_NI_p - v_SI_n) + 1.d0/Pe_m*2.d0/del_x**2 + 1.d0/Pe_m*2.d0/del_y**2
        C_e_FOU(i) = 1.d0/del_x*u_EI_n - 1.d0/Pe_m/del_x**2
        C_s_FOU(i) = -1.d0/del_y*v_SI_p - 1.d0/Pe_m/del_y**2 
        C_n_FOU(i) = 1.d0/del_y*v_NI_n - 1.d0/Pe_m/del_y**2
        RHS(i) = -C_s_FOU(i)*CVD(i,j-1) - C_n_FOU(i)*CVD(i,j+1)
    else if (position == 'Center') then
        if (DCM == 1) then
            C_w_FOU(i) = -1.d0/del_x*u_WI_p - 1.d0/Pe_m/del_x**2
            C_p_FOU(i) = 1.d0/del_x*(u_EI_p - u_WI_n) + 1.d0/del_y*(v_NI_p - v_SI_n) + 1.d0/Pe_m*2.d0/del_x**2 + 1.d0/Pe_m*2.d0/del_y**2
            C_e_FOU(i) = 1.d0/del_x*u_EI_n - 1.d0/Pe_m/del_x**2
            C_s_FOU(i) = -1.d0/del_y*v_SI_p - 1.d0/Pe_m/del_y**2 
            C_n_FOU(i) = 1.d0/del_y*v_NI_n - 1.d0/Pe_m/del_y**2
            C_w_HO(i) = (1.d0/del_x*(-beta*u_EI_p - u_WI_p*(0.5d0 + 2*beta) - u_WI_n*(0.5d0 - beta)) - 1.d0/Pe_m/del_x**2) - C_w_FOU(i)
            C_p_HO(i) = (1.d0/del_x*(u_EI_p*(0.5d0 + 2*beta) + u_EI_n*(0.5d0 - beta) - u_WI_p*(0.5d0 - beta) - u_WI_n*(0.5d0 + 2*beta))&
                    & + 1.d0/del_y*(v_NI_p*(0.5d0 + 2*beta) + v_NI_n*(0.5d0 - beta) - v_SI_p*(0.5d0 - beta) - v_SI_n*(0.5d0 + 2*beta))&
                    & + 1.d0/Pe_m*2.d0/del_x**2 + 1.d0/Pe_m*2.d0/del_y**2) - C_p_FOU(i)
            C_e_HO(i) = (1.d0/del_x*(beta*u_WI_n + u_EI_p*(0.5d0 - beta) + u_EI_n*(0.5d0 + 2*beta)) - 1.d0/Pe_m/del_x**2) - C_e_FOU(i)
            C_s_HO(i) = (1.d0/del_y*(-beta*v_NI_p - v_SI_p*(0.5d0 + 2*beta) - v_SI_n*(0.5d0 - beta)) - 1.d0/Pe_m/del_y**2) - C_s_FOU(i)
            C_n_HO(i) = (1.d0/del_y*(beta*v_SI_n + v_NI_p*(0.5d0 - beta) + v_NI_n*(0.5d0 + 2*beta)) - 1.d0/Pe_m/del_y**2) - C_n_FOU(i)
            C_ww_HO(i) = 1.d0/del_x*(beta*u_WI_p)
            C_ee_HO(i) = -1.d0/del_x*(beta*u_EI_n)
            C_ss_HO(i) = 1.d0/del_y*(beta*v_SI_p)
            C_nn_HO(i) = -1.d0/del_y*(beta*v_NI_n)
            RHS(i) = -C_s_FOU(i)*CVD(i,j-1) - C_n_FOU(i)*CVD(i,j+1) - G*(C_w_HO(i)*CVD(i-1,j) + C_p_HO(i)*CVD(i,j) + C_e_HO(i)*CVD(i+1,j) &
                    &+ C_s_HO(i)*CVD(i,j-1) + C_n_HO(i)*CVD(i,j+1) + C_ww_HO(i)*CVD(i-2,j) + C_ee_HO(i)*CVD(i+2,j) + C_ss_HO(i)*CVD(i,j-2) + C_nn_HO(i)*CVD(i,j+2))
        else if (DCM == 0) then
            C_w_HO(i) = (1.d0/del_x*(-beta*u_EI_p - u_WI_p*(0.5d0 + 2*beta) - u_WI_n*(0.5d0 - beta)) - 1.d0/Pe_m/del_x**2)
            C_p_HO(i) = (1.d0/del_x*(u_EI_p*(0.5d0 + 2*beta) + u_EI_n*(0.5d0 - beta) - u_WI_p*(0.5d0 - beta) - u_WI_n*(0.5d0 + 2*beta))&
                    & + 1.d0/del_y*(v_NI_p*(0.5d0 + 2*beta) + v_NI_n*(0.5d0 - beta) - v_SI_p*(0.5d0 - beta) - v_SI_n*(0.5d0 + 2*beta))&
                    & + 1.d0/Pe_m*2.d0/del_x**2 + 1.d0/Pe_m*2.d0/del_y**2)
            C_e_HO(i) = (1.d0/del_x*(beta*u_WI_n + u_EI_p*(0.5d0 - beta) + u_EI_n*(0.5d0 + 2*beta)) - 1.d0/Pe_m/del_x**2)
            C_s_HO(i) = (1.d0/del_y*(-beta*v_NI_p - v_SI_p*(0.5d0 + 2*beta) - v_SI_n*(0.5d0 - beta)) - 1.d0/Pe_m/del_y**2)
            C_n_HO(i) = (1.d0/del_y*(beta*v_SI_n + v_NI_p*(0.5d0 - beta) + v_NI_n*(0.5d0 + 2*beta)) - 1.d0/Pe_m/del_y**2)
            C_ww_HO(i) = 1.d0/del_x*(beta*u_WI_p)
            C_ee_HO(i) = -1.d0/del_x*(beta*u_EI_n)
            C_ss_HO(i) = 1.d0/del_y*(beta*v_SI_p)
            C_nn_HO(i) = -1.d0/del_y*(beta*v_NI_n)
            RHS(i) = -C_s_HO(i)*CVD(i,j-1) - C_n_HO(i)*CVD(i,j+1) - C_ww_HO(i)*CVD(i-2,j) - C_ee_HO(i)*CVD(i+2,j) - C_ss_HO(i)*CVD(i,j-2) - C_nn_HO(i)*CVD(i,j+2)
        else
            write(*,*)'Mistake in DCM'
        end if
    else if (position == 'Right') then
        C_w_FOU(i) = -1.d0/del_x*u_WI_p - 1.d0/Pe_m/del_x**2
        C_p_FOU(i) = 1.d0/del_x*(u_p - u_WI_n) + 0.5d0/del_y*(v_NI_p - v_SI_n) + 1.d0/Pe_m/del_x**2 + 1.d0/Pe_m/del_y**2
        C_e_FOU(i) = 0
        C_s_FOU(i) = -0.5d0/del_y*v_SI_p - 0.5d0/Pe_m/del_y**2
        C_n_FOU(i) = 0.5d0/del_y*v_NI_n - 0.5d0/Pe_m/del_y**2
        RHS(i) = -C_s_FOU(i)*CVD(i,j-1) - C_n_FOU(i)*CVD(i,j+1)
    else if (position == 'Bottom') then
        C_w_FOU(i) = -0.5d0/del_x*u_WI_p - 0.5d0/Pe_m/del_x**2
        C_p_FOU(i) = 0.5d0/del_x*(u_EI_p - u_WI_n) + 1.d0/del_y*(v_NI_p - v_p) + 1.d0/Pe_m/del_x**2 + 1.d0/Pe_m/del_y**2
        C_e_FOU(i) = 0.5d0/del_x*u_EI_n - 0.5d0/Pe_m/del_x**2
        C_s_FOU(i) = 0
        C_n_FOU(i) = 1.d0/del_y*v_NI_n - 1.d0/Pe_m/del_y**2
        RHS(i) = - C_n_FOU(i)*CVD(i,j+1)
    else if (position == 'Top') then
        C_w_FOU(i) = -0.5d0/del_x*u_WI_p - 0.5d0/Pe_m/del_x**2
        C_p_FOU(i) = 0.5d0/del_x*(u_EI_p - u_WI_n) + 1.d0/del_y*(v_p - v_SI_n) + 1.d0/Pe_m/del_x**2 + 1.d0/Pe_m/del_y**2
        C_e_FOU(i) = 0.5d0/del_x*u_EI_n - 0.5d0/Pe_m/del_x**2
        C_s_FOU(i) = -1.d0/del_y*v_SI_p - 1.d0/Pe_m/del_y**2
        C_n_FOU(i) = 0
        RHS(i) = -C_s_FOU(i)*CVD(i,j-1)
    else if (position == 'Right_top') then
        C_w_FOU(i) = -0.5d0/del_x*u_WI_p - 0.5d0/Pe_m/del_x**2
        C_p_FOU(i) = 0.5d0/del_x*(u_p - u_WI_n) + 0.5d0/del_y*(v_p - v_SI_n) + 0.5d0/Pe_m/del_x**2 + 0.5d0/Pe_m/del_y**2
        C_e_FOU(i) = 0
        C_s_FOU(i) = -0.5d0/del_y*v_SI_p - 0.5d0/Pe_m/del_y**2
        C_n_FOU(i) = 0
        RHS(i) = -C_s_FOU(i)*CVD(i,j-1)
    else if (position == 'Right_bottom') then
        C_w_FOU(i) = -0.5d0/del_x*u_WI_p - 0.5d0/Pe_m/del_x**2
        C_p_FOU(i) = 0.5d0/del_x*(u_p - u_WI_n) + 0.5d0/del_y*(v_NI_p - v_p) + 0.5d0/Pe_m/del_x**2 + 0.5d0/Pe_m/del_y**2
        C_e_FOU(i) = 0
        C_s_FOU(i) = 0
        C_n_FOU(i) = 0.5d0/del_y*v_NI_n - 0.5d0/Pe_m/del_y**2
        RHS(i) = - C_n_FOU(i)*CVD(i,j+1)
    else
        write(*,*)'I type wrongQQ'
        pause
    end if
    
end subroutine coefficient