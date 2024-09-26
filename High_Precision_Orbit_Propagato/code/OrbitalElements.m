function [ h, mag_h, i, omega, e, mag_e, w, theta ] = OrbitalElements ( r,v,muo )
% This function is used to get OrbitalElements from position vector and its velocity vector
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% INPUTS:
% r        :  the position vector
% v       :  the position velocity
% muo :  Gravitational Parameter
%% OUTPUTS:
% h              : specific angular momentum vector
% mag_h     : specific angular momentum magnitude
% i               : inclination angle in degree
% omega     : right ascension of the ascending node in degree
% e              : eccentricity vector
% mag_e     : eccentricity magnitude
% w             : argument of perigee in degree
% theta       : true anomaly in degree
% ---------------------------------------------------------------------------------------------------------------------------------------------------------
mag_r=norm(r);
mag_v=norm(v);
vr=dot(r,v)/mag_r;
h=cross(r,v);
mag_h=norm(h);
i=acosd(h(3)/mag_h);
N=cross([0,0,1],h);
mag_N=norm(N);
if N(2) >= 0
    omega=acosd(N(1)/mag_N);
elseif N(2) < 0
    omega=360-acosd(N(1)/mag_N);
end
e=((mag_v^2-muo/mag_r)*r-mag_r*vr*v)/muo;
mag_e=norm(e);
if e(3) >= 0
    w=acosd(dot(N,e)/mag_N/mag_e);
elseif e(3) < 0
    w=360-acosd(dot(N,e)/mag_N/mag_e);
end
if vr >= 0
    theta=acosd(dot(e,r)/mag_r/mag_e);
elseif vr < 0
    theta=360-acosd(dot(e,r)/mag_r/mag_e);
end
end