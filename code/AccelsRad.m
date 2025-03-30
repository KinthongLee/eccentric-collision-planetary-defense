function a = AccelsRad(x,  AU)
    % Define constants
    A1 = 5e-13;

    rr = x ;
    % Compute distance to the Sun
    r = norm(rr);
    
    % Compute the solar radiation pressure force

    a1 = A1 * 1 / (r/AU)^2 * rr / r;
    a1 = a1 * AU / (86400)^2;
    
    % Return the acceleration due to non-gravitational forces
    a = a1;
end