%--------------------------------------------------------------------------
%
% Relativisty: Computes the perturbational acceleration due to relativistic
%              effects
%
% Inputs:
%   r           Satellite position vector
%   v           Satellite velocity vector
% 
% Output:
%   a    		Acceleration (a=d^2r/dt^2)
%
% Last modified:   2018/01/27   Meysam Mahooti
%
%--------------------------------------------------------------------------
function a = Relativity(r, v)


% Relative position vector of satellite w.r.t. point mass 
r_Sat = norm(r);
v_Sat = norm(v);

% Acceleration 
a = 132712440041279422464/(299792458^2*r_Sat^3)*((4*132712440041279422464/r_Sat-v_Sat^2)*r+4*dot(r,v)*v);

