function orbitalElements = calculateOrbitalElements(r, v)
    % Constants
    mu = 132712440041279422464; % Earth's gravitational parameter in km^3/s^2
    
    % Calculate the specific angular momentum
    h = cross(r, v);
    h_norm = norm(h);
    
    % Calculate the node vector
    K = [0 0 1]; % Unit vector along the Z-axis
    N = cross(K, h);
    N_norm = norm(N);
    
    % Calculate the eccentricity vector
    e_vec = (1/mu) * ((norm(v)^2 - mu/norm(r)) * r - dot(r, v) * v);
    e = norm(e_vec);
    
    % Calculate the inclination
    i = acos(h(3) / h_norm); % Inclination in radians
    
    % Calculate the longitude of the ascending node
    Omega = acos(N(1) / N_norm); % Longitude of ascending node in radians
    if N(2) < 0
        Omega = 2*pi - Omega;
    end
    
    % Calculate the argument of periapsis
    omega = acos(dot(N, e_vec) / (N_norm * e));
    if e_vec(3) < 0
        omega = 2*pi - omega;
    end
    
    % Calculate the true anomaly
    nu = acos(dot(e_vec, r) / (e * norm(r)));
    if dot(r, v) < 0
        nu = 2*pi - nu;
    end
    
    % Calculate the semi-major axis
    a = 1 / ((2/norm(r)) - (norm(v)^2 / mu));
    
    % Collect all orbital elements
    orbitalElements = struct('a', a, 'e', e, 'i', i, 'Omega', Omega, ...
                             'omega', omega, 'nu', nu);
end
