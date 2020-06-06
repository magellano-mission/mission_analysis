function [lb, ub] = LBUB(XX0, X,  data, Bounds)
%X0 vector of propagated ODE45
    N = data.n_int;
    T_lb = 0;
    T_ub = data.Tmax;

    %limits definition
    ub          = zeros(1,10*N);
    lb          = zeros(1,10*N);
    
    
    alpha_lb    = Bounds.alpha_lb;
    alpha_ub    = Bounds.alpha_ub;
    beta_lb     = Bounds.beta_lb;
    beta_ub     = Bounds.beta_ub;
        
    %r
    r_ub        = Bounds.r_ub;
    r_lb        = Bounds.r_lb;

    %theta
    th_ub       = Bounds.th_ub;
    th_lb       = Bounds.th_lb;
    %z
    z_ub        =  Bounds.z_ub;
    z_lb        =  Bounds.z_lb;
    %vr
    vr_ub       = Bounds.vr_ub;
    vr_lb       = Bounds.vr_lb;
    %theta_dot
    thd_ub      = Bounds.thd_ub;
    thd_lb      = Bounds.thd_lb;
    %vz
    vz_ub       = Bounds.vz_ub;
    vz_lb       = Bounds.vz_lb;
    
    %
    m_ub        = X(1,7)*1.1;
    m_lb        = 1;
    %INITIAL CONDITION
    ub(1:10) = XX0(1:10); 
    ub(8) = T_ub; ub(9) = alpha_ub; ub(10) = beta_ub;
    ub(7) = m_ub;
    lb(1:10) = XX0(1:10); 
    lb(7) = m_lb;
    lb(8) = 0; lb(9) = alpha_lb; lb(10) = beta_lb;
    
    %PATH CONSTRAINTS
    ub(11:10:((N-1)*10)) = r_ub;
    ub(12:10:((N-1)*10)) = th_ub;
    ub(13:10:((N-1)*10)) = z_ub;
    ub(14:10:((N-1)*10)) = vr_ub;
    ub(15:10:((N-1)*10)) = thd_ub;
    ub(16:10:((N-1)*10)) = vz_ub;
    ub(17:10:((N-1)*10)) = m_ub;
    ub(18:10:((N-1)*10)) = T_ub;
    ub(19:10:((N-1)*10)) = alpha_ub;
    ub(20:10:((N-1)*10)) = beta_ub;

    lb(11:10:((N-1)*10)) = r_lb;
    lb(12:10:((N-1)*10)) = th_lb;
    lb(13:10:((N-1)*10)) = z_lb;
    lb(14:10:((N-1)*10)) = vr_lb;
    lb(15:10:((N-1)*10)) = thd_lb;
    lb(16:10:((N-1)*10)) = vz_lb;
    lb(17:10:((N-1)*10)) = m_lb;%data.Mdry;
    lb(18:10:((N-1)*10)) = T_lb;
    lb(19:10:((N-1)*10)) = alpha_lb;
    lb(20:10:((N-1)*10)) = beta_lb;

    %FINAL CONDITION
    ub(end-9:end) = XX0(end-9:end);
    ub(end-2)=T_ub; ub(end-1)=alpha_ub; ub(end)=beta_ub;
    ub(end-3) = m_lb;
    lb(end-9:end) = XX0(end-9:end);
    lb(end-2)=0; lb(end-1)=alpha_lb; lb(end)=beta_lb;
    lb(end-3) = m_lb;
end