function [beta,gamma] = cal_eff(theta_)% theta[rad]
    num = length(theta_);
    beta = zeros(num,1);
    gamma = beta;
    % 数值结果用于插值
    angle = pi/6:pi/36:pi/2;
    data = [   2.211334939007717
               2.257745058522270
               2.305480842671130
               2.352808685962787
               2.397933307544245
               2.438912767230096
               2.473763347853509
               2.500859287024754
               2.519505366420655
               2.530348426957326
               2.535293768921257
               2.536840579999872
               2.537084725580740];
    A = 1.32;
    % 计算结果
    for i=1:num
        theta = theta_(i);
        if theta > pi/2
            theta = theta - pi/2;
        end
%         if theta<pi/6 || theta>pi/2
%             fprintf('Input out of range in function cal_beta!\n');
%             continue;
%         end
        beta(i) = interp1(angle,data,theta,'spline');% 样条插值
        phi_m = asin(A*sin(theta)-A);
        gamma(i) = (beta(i)-1)*tan(theta)*tan(phi_m)+1;
        disp(theta*180/pi)
        disp(beta(i))
        disp(gamma(i))
    end
end


