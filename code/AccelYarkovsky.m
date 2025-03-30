function a = AccelYarkovsky(x,  AU,v)
    % Define constants
    A2 = -2.899e-14;

    rr = x ;
    % Compute distance to the Sun
    r = norm(rr);
    
    % Compute the Yarkovsky effect force

    a2 = A2 * 1 / (r/AU)^2 * v / norm(v);
    a2 = a2 * AU / (86400)^2;
    
    % Return the acceleration due to non-gravitational forces
    a = a2;
end