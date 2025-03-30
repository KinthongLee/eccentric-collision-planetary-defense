function [best_delta_v,delta_v_COM,angles_output] = get_delta_v_from_momentum(v_imp,v_Apophis,mass_rocket,vertices,faces,beta)

m = mass_rocket; %kg
M = 2.7e10; %kg
v_r = v_imp-v_Apophis;
v_r_norm = v_r / norm(v_r);

% Define rotation angles
angles = [rand()*pi;rand()*pi;rand()*pi];

% Build rotation matrix
Rx = [cos(angles(1)),sin(angles(1)),0;-sin(angles(1)),cos(angles(1)),0;0,0,1];
Ry = [cos(angles(2)),0,sin(angles(2));0,1,0;-sin(angles(2)),0,cos(angles(2))];
Rz = [1,0,0;0,cos(angles(3)),sin(angles(3));0,-sin(angles(3)),cos(angles(3))];
R = Rx*Ry*Rz;

% Apply the rotation matrix to the vertices of the model
rotated_vertices = zeros(length(vertices),1);
for i = 1 : length(vertices)
    rotated_vertices(i,1:3) = (R * vertices(i,1:3)')';
end



% Compute the normal vectors of the every faces
face_normals = cross(rotated_vertices(faces(:,2),:)-rotated_vertices(faces(:,1),:), ...
                     rotated_vertices(faces(:,3),:)-rotated_vertices(faces(:,1),:));
face_normals = face_normals ./ vecnorm(face_normals, 2, 2);
irradiated_vector =  v_r_norm;
% Find all the surface which irradiated by the velocity vector
% Compute the dot product between the face normals and the light direction

dot_products = dot(face_normals, repmat(irradiated_vector', size(faces, 1), 1), 2);


% Find the indices of the irradiated faces
irradiated_face_indices = dot_products < 0;


% Calculate the tangential_vector_t & theta_impact of every faces
n = -face_normals(irradiated_face_indices,:);  %normal components of face
k = [0,0,1];
t = zeros(length(n),3);  % tangential components of face
theta_impact = zeros(length(n),1);
for i = 1 : length(n)
    t(i,:) = cross( cross(k,n(i,:)) , n(i,:) ) / norm( cross( k,n(i,:) ) );
    theta_impact(i) = acos(  dot( v_r_norm,t(i,:)' ) / norm(v_r_norm) / norm(t(i,:)')  ) ;
end



gamma = zeros(length(n),1);
for i = 1 : length(theta_impact)
    if theta_impact(i) > pi/2
        x = theta_impact(i) - pi/2;
    else
        x = theta_impact(i);
    end
    x = x*180/pi;
    gamma(i) =1 +  (beta-1)*tan(x*pi/180)*tan((-0.0052*x^2+1.3629*x-80.871)*pi/180);
end


% get delta_v
norm_delta_v_A = zeros(length(gamma),1);
delta_v = zeros(length(gamma),3);
% delta_v_mass_point = zeros(length(beta),3);
for i = 1 : length(gamma)
    delta_v(i,1:3) = m/M * (...
        beta*dot( v_r,n(i,:)' )*n(i,:)' + ...
        gamma(i)*dot(v_r,t(i,:)')*t(i,:)' ...
        );
    norm_delta_v_A(i,1) = abs( dot(delta_v(i,1:3),v_Apophis) / norm(v_Apophis) );
end
 [~,b] = max( norm_delta_v_A);
best_delta_v = delta_v(b,1:3);
angles_output(1,:) = angles(:,1);


% Through Center of Mass
% Calculate the faces that irradiated_vector pass though its center of mass
center_of_mass = mean(vertices);
ray_origin = center_of_mass;
ray_end = ray_origin - irradiated_vector' * 10;
effective_faces = faces(irradiated_face_indices, :);
vert1 = rotated_vertices(effective_faces(:,1),:);
vert2 = rotated_vertices(effective_faces(:,2),:);
vert3 = rotated_vertices(effective_faces(:,3),:);
% Check if the vector passes through the center of mass of the face
intersect = TriangleRayIntersection(ray_origin, ray_end, vert1, vert2, vert3);
face_through_COM_indices = intersect == 1;

% Calculate the tangential_vector_t & theta_impact of COM
n_COM = -face_normals(irradiated_face_indices,:);  %normal components of face
n_COM = n_COM(face_through_COM_indices,:);
k = [0,0,1];
t_COM = cross( cross(k,n_COM) , n_COM ) / norm( cross( k,n_COM ) );
theta_impact_COM = acos(  dot( v_r_norm,t_COM' ) / norm(v_r_norm) / norm(t_COM')  ) ;



    if theta_impact_COM > pi/2
        x_COM = theta_impact_COM - pi/2;
    else
        x_COM = theta_impact_COM;
    end
    x_COM = x_COM*180/pi;
    gamma_COM =1 +  (beta-1)*tan(x_COM*pi/180)*tan((-0.0052*x_COM^2+1.3629*x_COM-80.871)*pi/180);



% COM
delta_v_COM = m/M * (...
        beta*dot( v_r,n_COM' )*n_COM' + ...
        gamma_COM*dot(v_r,t_COM')*t_COM' ...
        );


end
