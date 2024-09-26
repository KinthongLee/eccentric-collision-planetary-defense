%--------------------------------------------------------------------------
%
% ECEF2ECI: Transforms Earth Centered Earth Fixed (ECEF) coordinates to 
%           Earth Centered Inertial (ECI) coordinates
%
% Last modified:   2022/11/07   Meysam Mahooti
%
%--------------------------------------------------------------------------
function Y = ECEF2ECI(MJD_UTC, Y0)

global const eopdata

[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
MJD_UT1 = MJD_UTC + UT1_UTC/86400;
MJD_TT = MJD_UTC + TT_UTC/86400;

% Form bias-precession-nutation matrix
NPB = iauPnm06a(const.DJM0, MJD_TT);
% Form Earth rotation matrix
Theta = iauRz( iauGst06(const.DJM0, MJD_UT1, const.DJM0, MJD_TT, NPB),eye(3) );
% Polar motion matrix (TIRS->ITRS, IERS 2003)
Pi = iauPom00(x_pole, y_pole, iauSp00(const.DJM0, MJD_TT));

% ICRS to ITRS transformation matrix and derivative
S = zeros(3);
S(1,2) = 1;
S(2,1) = -1;              			 % Derivative of Earth rotation
% Omega  = 7292115.8553e-11+4.3e-15*( (MJD_UTC-const.MJD_J2000)/36525 ); % [rad/s]
Omega = const.omega_Earth-0.843994809*1e-9*LOD;  % IERS
dTheta = Omega*S*Theta;           				 % matrix [1/s]
U      = Pi*Theta*NPB;                			 % ICRS to ITRS transformation
dU     = Pi*dTheta*NPB;               			 % Derivative [1/s]

% Transformation from WGS to ICRS
r = U'*Y0(1:3)';
v = U'*Y0(4:6)' + dU'*Y0(1:3)';
Y = [r;v];

