!-------------------------------------------------------------------
subroutine uv(N,i,j,u_EI_p, u_EI_n, u_WI_p, u_WI_n, v_NI_p, v_NI_n, v_SI_p, v_SI_n,u_p,v_p,u,v,Nv)
    implicit NONE
    
    integer :: i, j, Nv, N
    real*8 :: u_EI_p,u_EI_n,u_WI_p,u_WI_n,v_NI_p,v_NI_n,v_SI_p,v_SI_n,u_p,v_p
    real*8, dimension(Nv,Nv) :: u,v
    
    
    if (i.ne.N) then
        u_EI_p = (u(2*i,2*j-1)+abs(u(2*i,2*j-1)))/2
        u_EI_n = (u(2*i,2*j-1)-abs(u(2*i,2*j-1)))/2
    end if
    if (i.ne.1) then
        u_WI_p = (u(2*i-2,2*j-1)+abs(u(2*i-2,2*j-1)))/2
        u_WI_n = (u(2*i-2,2*j-1)-abs(u(2*i-2,2*j-1)))/2
    end if
    if (j.ne.N) then
        v_NI_p = (v(2*i-1,2*j)+abs(v(2*i-1,2*j)))/2
        v_NI_n = (v(2*i-1,2*j)-abs(v(2*i-1,2*j)))/2
    end if
    if (j.ne.1) then
        v_SI_p = (v(2*i-1,2*j-2)+abs(v(2*i-1,2*j-2)))/2
        v_SI_n = (v(2*i-1,2*j-2)-abs(v(2*i-1,2*j-2)))/2
    end if
    u_p = u(2*i-1,2*j-1)
    v_p = v(2*i-1,2*j-1)
    
end subroutine uv