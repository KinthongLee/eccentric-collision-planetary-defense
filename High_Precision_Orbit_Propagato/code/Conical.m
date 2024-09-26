%--------------------------------------------------------------------------
% 
% Conical
%
% Purpose:
%   Computes the fractional illumination of a spacecraft in the 
%   vicinity of the Earth assuming a conical shadow model
% 
% Inputs:
%   r         Spacecraft position vector [m]
%   r_Earth   Earth position vector (barycentric) [m]
%   r_Sun     Sun position vector (geocentric) [m]
%   
% Output:
%   nu        Illumination factor:
%             nu=0   Spacecraft in Earth shadow
%             0<nu<1 Spacecraft in Penumbra
%             nu=1   Spacecraft fully illuminated by the Sun
%
% Reference:
% Montenbruck O., and Gill E., "Satellite Orbits: Models, Methods, and 
% Applications," Springer Verlag, Heidelberg, Corrected 3rd Printing (2005).
% 
% Last modified:   2022/10/26   Meysam Mahooti
% 
%--------------------------------------------------------------------------
function nu = Conical(r, r_Earth, r_Sun)

global const

s = r-r_Earth;
smag = norm(s);
d     = r_Sun-r;
dmag  = norm(d);

a = asin(const.R_Sun/dmag);          % eq. 3.85
b = asin(const.R_Earth/smag);        % eq. 3.86
c = acos(-1.0*dot(s,d)/(smag*dmag)); % eq. 3.87

if( c >= (a+b) )            % in Sun light
    nu = 1.0;
elseif( c < abs(a-b) )      % in Umbra
    nu =  0.0;
else                        % in Penumbra 
    x = (c^2+a^2-b^2)/(2*c);                  % eq. 3.93
    y = sqrt(a^2-x^2);
    A = a^2*acos(x/a)+b^2*acos((c-x)/b)-c*y;  % eq. 3.92
    nu = 1.0 - A/(pi*a^2);                    % eq. 3.94
end

