data = readtable('PHA_table');
% Calculate the deflection distance and Gain(%) through analytical method
% Write data into PHA_table.xlsx for specific column: BI to BL, and output
% to temporary result. 
% If the result is correct, please manually copy and paste the result into PHA_table.xlsx.

% Be noted true anomaly here is stored in variable named MOID, however it
% is actually true anomaly at the moment of Closest-Approach, this will
% increase the result instead of MOID
for ppp = 33:33
    if ~isnan(data.Best_Impact_Year(ppp))
        v_imp = [data.v_imp_X(ppp);data.v_imp_Y(ppp);data.v_imp_Z(ppp)];
        v_ast = [data.v_ast_X(ppp);data.v_ast_Y(ppp);data.v_imp_Z(ppp)];
        r_ast = [data.x_ast_impact(ppp);data.y_ast_impact(ppp);data.z_ast_impact(ppp)];
        delta_v_COG = [data.delta_v_t_COG(ppp);data.delta_v_n_COG(ppp);data.delta_v_h_COG(ppp)];
        delta_v_BIP = [data.delta_v_t_BIP(ppp);data.delta_v_n_BIP(ppp);data.delta_v_h_BIP(ppp)];
        a = data.a(ppp)/1000;
        e = data.e(ppp);
        i = data.i(ppp);
        omega = data.omega(ppp);
        w = data.w(ppp);
        theta_d = data.f(ppp); 
        theta_MOID = data.f_MOID(ppp); 
        Delta_t = data.delta_t(ppp);
        r_MOID = data.r_min_by_matlab(ppp)/1000;
        theta_MOID_star = theta_MOID + w;
        theta_d_star = theta_d + w;
        eta = sqrt(1-e^2);
        mu = 132712440041.279422464;
        % 换成km
        v_d = norm(v_ast)/1000;
        r_d = norm(r_ast)/1000;
        h = norm(cross(r_ast,v_ast))/1000000;
        b = a*sqrt(1-e^2);
        p = a * (1-e^2);
        
        % % Define the matrix A_MOID
        A_MOID = [
            r_MOID/a - 3*e*sind(theta_MOID)*sqrt(mu)/(2*eta*(a^(3/2)))*Delta_t , - 3*r_MOID/(2*eta^3)*(1 + e*cosd(theta_MOID))^2*sqrt(mu/(a^5))*Delta_t, 0;
            -a*cosd(theta_MOID), r_MOID*sind(theta_MOID)/eta^2*(2 + e*cosd(theta_MOID)), 0;
            0 , 0, r_MOID*sind(theta_MOID_star) ;
            0 ,  r_MOID*cosd(i) , -r_MOID*cosd(theta_MOID_star)*sind(i);
            0 , r_MOID , 0 ;
            a*e*sind(theta_MOID)/eta, r_MOID/(eta^3)*(1 + e*cosd(theta_MOID))^2, 0
        ]';
        % Define the matrix B
        G_d = [
            2*a^2*v_d/mu, 0, 0;
            2*(e + cosd(theta_d))/v_d, -r_d*sind(theta_d)/(a*v_d), 0;
            0, 0, r_d*cosd(theta_d_star)^2/h;
            0, 0, r_d*sind(theta_d_star)^2/(h*sind(i));
            2*sind(theta_d)/(e*v_d), (2*e+r_d/a*cosd(theta_d))/(e*v_d), -r_d*sind(theta_d_star)*cosd(i)/(h*sind(i));
            -b*2*(1 + e^2*r_d/p)*sind(theta_d)/(e*a*v_d), -b*r_d*cosd(theta_d)/(e*a^2*v_d), 0
        ];
        

        delta_r_COG = A_MOID*G_d*(delta_v_COG./1000);
        delta_r_BIP = A_MOID*G_d*(delta_v_BIP./1000);
        data.delta_r_COG_theory(ppp) = norm(delta_r_COG)*1000;
        data.delta_r_BIP_theory(ppp) = norm(delta_r_BIP)*1000;
        data.Gain_theory(ppp) = norm(delta_r_BIP) - norm(delta_r_COG);
        data.Gain_theory_percent(ppp) = data.Gain_theory(ppp) / norm(delta_r_COG)*100;


    end
end

writetable(data,'temporary_result.xlsx')