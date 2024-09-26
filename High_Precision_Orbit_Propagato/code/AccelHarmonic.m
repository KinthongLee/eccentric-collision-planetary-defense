%--------------------------------------------------------------------------
%
% AccelHarmonic: Computes the acceleration due to the harmonic gravity
%                field of the central body
%
% Inputs:
%   r           Satellite position vector in the inertial system
%   E           Transformation matrix to body-fixed system
%   Cnm, Snm    Spherical harmonics coefficients (normalized)
%   n_max       Maximum degree
%   m_max       Maximum order (m_max<=n_max; m_max=0 for zonals, only)
%
% Output:
%   a           Acceleration (a=d^2r/dt^2)
%
% Last modified:   2022/06/16   Meysam Mahooti
% 
%--------------------------------------------------------------------------
function a = AccelHarmonic(r, E, n_max, m_max)

global Cnm Snm

gm    = 398600.4415e9; % [m^3/s^2]; GGM03C & GGM03S
r_ref = 6378.1363e3;   % Earth's radius [m]; GGM03C & GGM03S

% Body-fixed position 
r_bf = E * r;

% Auxiliary quantities
d = norm(r_bf);                     % distance
latgc = asin(r_bf(3)/d);
lon = atan2(r_bf(2),r_bf(1));

[pnm, dpnm] = Legendre(n_max,m_max,latgc);

dUdr = 0;
dUdlatgc = 0;
dUdlon = 0;
q3 = 0; q2 = q3; q1 = q2;
for n=0:n_max
    b1 = (-gm/d^2)*(r_ref/d)^n*(n+1);
    b2 =  (gm/d)*(r_ref/d)^n;
    b3 =  (gm/d)*(r_ref/d)^n;
    for m=0:m_max
        q1 = q1 + pnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
        q2 = q2 + dpnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
        q3 = q3 + m*pnm(n+1,m+1)*(Snm(n+1,m+1)*cos(m*lon)-Cnm(n+1,m+1)*sin(m*lon));
    end
    dUdr     = dUdr     + q1*b1;
    dUdlatgc = dUdlatgc + q2*b2;
    dUdlon   = dUdlon   + q3*b3;
    q3 = 0; q2 = q3; q1 = q2;
end

% Body-fixed acceleration
r2xy = r_bf(1)^2+r_bf(2)^2;

ax = (1/d*dUdr-r_bf(3)/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf(1)-(1/r2xy*dUdlon)*r_bf(2);
ay = (1/d*dUdr-r_bf(3)/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf(2)+(1/r2xy*dUdlon)*r_bf(1);
az =  1/d*dUdr*r_bf(3)+sqrt(r2xy)/d^2*dUdlatgc;

a_bf = [ax ay az]';

% Inertial acceleration 
a = E'*a_bf;

