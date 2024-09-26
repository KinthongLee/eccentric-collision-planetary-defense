function intersects = rayTriangleIntersection(rayOrigin, rayDirection, triangleVertices)
    % Ray-triangle intersection using Möller–Trumbore algorithm
    
    v0 = triangleVertices(1, :);
    v1 = triangleVertices(2, :);
    v2 = triangleVertices(3, :);
    
    e1 = v1 - v0;
    e2 = v2 - v0;
    
    h = cross(rayDirection, e2);
    a = dot(e1, h);
    
    if a > -eps && a < eps
        intersects = false; % Ray is parallel to the triangle
        return;
    end
    
    f = 1 / a;
    s = rayOrigin - v0;
    u = f * dot(s, h);
    
    if u < 0 || u > 1
        intersects = false;
        return;
    end
    
    q = cross(s, e1);
    v = f * dot(rayDirection, q);
    
    if v < 0 || u + v > 1
        intersects = false;
        return;
    end
    
    t = f * dot(e2, q);
    
    intersects = t >= 0; % Check if t is greater than or equal to zero
end