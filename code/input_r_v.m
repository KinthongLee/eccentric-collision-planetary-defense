function [a,e,i,omega,w,f] = input_r_v(r,v,mu)
%   rx,ry,rz各别为位置矢量的X,Y,Z分量，单位为m（米）
%   vx,vy,vz各别为速度矢量的X,Y,Z分量，单位为m/s(米/秒)
%   mu:中心引力场引力系数,万有引力常数与中心天体质量的乘积,单位m^3/s^2,
%      默认为太阳引力系数132712440041279422464


if nargin<3
    mu=132712440041279422464;
    if nargin < 2
        error('NotEnoughInputs',...
              'Not enough input arguments'); 
    end
end


h = cross(r,v);
i = acosd(h(3)/norm(h));
omega = atan((-1)*h(1)/h(2))*180/pi;
a = (-mu/2)/(0.5*(norm(v))^2-mu/norm(r));
p = (norm(h))^2/mu;
vector_e = cross(v,h)/mu - r/norm(r);
e = norm(vector_e);
if(e<0.0)
    error('ORBIT:f2E:NotSuitableInputs',...
                     '偏心率不能为负.  See f2E.');
elseif(e>=1)
        error('ORBIT:f2E:NotSuitableInputs',...
                   '双曲线/抛物线没有定义的轨道根数.');
else
ON = [cos(omega); sin(omega); 0];
w = acosd((ON'*vector_e)/e);
f = acosd(vector_e'*r/(norm(r)*norm(e)));
end
end

