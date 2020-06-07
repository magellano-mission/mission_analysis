function data = polar2cart(X, Tpolar, alpha, beta, r1vers, RCRRv)
    r  = X(:,1);    
    th = X(:,2);
    z  = X(:,3);    
    vr = X(:,4);
    vt = r.* X(:,5);   
    vz = X(:,6);
    
    rr = zeros(length(Tpolar), 3); vv = rr; T_cart = rr;

    for i = 1:length(Tpolar)
        [rr(i,:), vv(i,:)] = refplane2car( r(i), z(i),  vt(i), vr(i), vz(i), th(i), r1vers, RCRRv);
        T_cart(i,:) = thrust_cart2( Tpolar(i), alpha(i), beta(i), th(i), r1vers, RCRRv);
    end

    data.T_cart = T_cart;
    data.r_cart = rr;
    data.v_cart = vv;
    data.m = X(:,7);
    
end

