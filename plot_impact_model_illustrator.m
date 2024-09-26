clearvars
clc

% Display diagrams illustrating the two different strategies for targeting BIP and COG.

addpath('code\')
addpath('code\3D_model\')
% The state of impactor and Apophis at the moment of imapct at Apophis
% defense mission
% Modify if neccesary
% Impactor
v_imp = [5512.4087760978027290548197925091, -26735.322300538853596663102507591, -11345.535364285360628855414688587]';
% Asteroid (Apophis)
v_ast = [7026.8425407801;-24480.7553464794;-8922.30954540502];
vt = v_ast / norm(v_ast);
% Mass of impactor
m = 5955.12; %kg
% Mass of Asteroid
M = 2.7e10; %kg

% % --------------------------------------------------------------
% % Use to calculate different interception angle by rotating the vector v_imp through the axis v_imp x v_ast, mark it when not neccessory
% % Specify the desired angle theta2 (larger than theta1)
% % theta 1 = 11.450412891 deg
% % theta2 = 18.54958711*pi/180;  % Replace delta with the desired increase in angle
% theta2 = 30*pi/180;
% 
% % Rotate vector B to create vector C with the desired angle theta2
% axis_of_rotation = cross( v_ast,v_imp );  % Cross product gives the rotation axis
% axis_of_rotation = axis_of_rotation / norm(axis_of_rotation);  % Normalize the axis
% 
% % Create a rotation matrix for rotating around the specified axis
% R = axang2rotm([axis_of_rotation' theta2]);  % Convert axis-angle to rotation matrix
% 
% % Apply the rotation to vector B to get vector C
% C = R * v_imp;
% 
% % % Verify the angles between A and C
% % theta3 = acosd(dot(v_Apophis, C) / (norm(v_Apophis) * norm(C)));
% v_imp = C;
% % ---------------------------------------------------------------


% Read 3D models
[vertices, faces] = readObj('Apophis_Model.obj');

% Define rotation angles
% Here use rand() to generate random attitude of the asteroid
angles = [rand()*pi;rand()*pi;rand()*pi];
% Build rotation matrix
Rx = [cos(angles(1)),sin(angles(1)),0;-sin(angles(1)),cos(angles(1)),0;0,0,1];
Ry = [cos(angles(2)),0,sin(angles(2));0,1,0;-sin(angles(2)),0,cos(angles(2))];
Rz = [1,0,0;0,cos(angles(3)),sin(angles(3));0,-sin(angles(3)),cos(angles(3))];
R = Rx*Ry*Rz;

% Apply the rotation matrix to the vertices of the 3D model
rotated_vertices = zeros(length(vertices),1);
for i = 1 : length(vertices)
    rotated_vertices(i,1:3) = (R * vertices(i,1:3)')';
end

% Plot
% Plot the 3D model
figure(1)
p = patch('Faces',faces,'Vertices',rotated_vertices,'EdgeColor','black','FaceColor','blue','FaceAlpha', 0.1);
% Set plot properties
axis equal;
view(3);
xlabel('x')
ylabel('y')
zlabel('z')
hold on

% Calculate the Impact model
% Compute the normal vectors of the every faces
face_normals = cross(rotated_vertices(faces(:,2),:)-rotated_vertices(faces(:,1),:), ...
                     rotated_vertices(faces(:,3),:)-rotated_vertices(faces(:,1),:));
face_normals = face_normals ./ vecnorm(face_normals, 2, 2);


% Assume the direction vector of the relative velocity
v_r = v_imp-v_ast;
v_r_norm = v_r / norm(v_r);
irradiated_vector =  v_r_norm;
% Find all the surface which irradiated by the velocity vector
% Compute the dot product between the face normals and the light direction
dot_products = dot(face_normals, repmat(irradiated_vector', size(faces, 1), 1), 2);

% Find the vertices that are irradiated
irradiated_vertices_indices = unique(faces(dot_products < 0, :));

% Find the indices of the irradiated faces
irradiated_face_indices = find(dot_products < 0);


% Calculate the tangential_vector_t & theta_impact of every faces
n = -face_normals(irradiated_face_indices,:);  %normal components of face
k = [0,0,1];
t = zeros(length(n),3);  % tangential components of face
theta_impact = zeros(length(n),1);
for i = 1 : length(n)
    t(i,:) = cross( cross(k,n(i,:)) , n(i,:) ) / norm( cross( k,n(i,:) ) );
    theta_impact(i) = acos(  dot( v_r_norm,t(i,:)' ) / norm(v_r_norm) / norm(t(i,:)')  ) ;
end


effective_theta_impact= theta_impact;
effective_faces = faces(irradiated_face_indices, :);

% Beta coefficient, modify if necessary
beta = 3.61; 
% Gamma calculate through fitting equation as mentioned in article
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

% calculate delta_v
norm_delta_v= zeros(length(gamma),3);
delta_v_A_to_Apophis = zeros(length(gamma),1);
delta_v = zeros(length(gamma),3);
for i = 1 : length(gamma)
    delta_v(i,:) = m/M * (...
        beta*dot( v_r,n(i,:) )*n(i,:) + ...
        gamma(i)*dot(v_r,t(i,:))*t(i,:) ...
        );
    delta_v_A_to_Apophis(i,1) = norm( dot(delta_v(i,:)', vt )  );
end
 [~,b] = max(delta_v_A_to_Apophis);


patch('Faces', effective_faces,'Vertices', rotated_vertices, 'FaceColor','flat','FaceVertexCData',delta_v_A_to_Apophis)
ylabel(colorbar,'Velocity changes in direction of asteroid orbital velocity (m/s)')
hold on

% Calculate the surface hit by COG strategy
% Calculate the faces that irradiated_vector pass though its center of mass
center_of_mass = mean(vertices);
ray_origin = center_of_mass;
ray_end = ray_origin - irradiated_vector' * 10;
vert1 = rotated_vertices(effective_faces(:,1),:);
vert2 = rotated_vertices(effective_faces(:,2),:);
vert3 = rotated_vertices(effective_faces(:,3),:);
% Check if the vector passes through the center of mass of the face
intersect = TriangleRayIntersection(ray_origin, ray_end, vert1, vert2, vert3);
patch('Faces', effective_faces(intersect,:),'Vertices', rotated_vertices, 'FaceColor','red')
hold on
quiver3(0,0,0,irradiated_vector(1),irradiated_vector(2),irradiated_vector(3))



% Convert to (vt,vn,vh) coordinate：
r_ast = [-155716246356.708282470703125; -17671019879.613208770751953125; -10510360729.52490234375];
vh = cross(r_ast,vt) / norm(cross(r_ast,vt) );
vn = cross(vh,v_ast) / norm(cross(vh,v_ast));
delta_v_BIP = delta_v(b,:)';
delta_v_COG = delta_v(intersect,:);

delta_v_BIP = [dot(delta_v_BIP,vt);dot(delta_v_BIP,vn);dot(delta_v_BIP,vh)];
delta_v_COG = [dot(delta_v_COG,vt);dot(delta_v_COG,vn);dot(delta_v_COG,vh)];
disp("BIP:")
disp(vpa(delta_v_BIP))
disp('COG:')
disp(vpa(delta_v_COG))
disp('Gain of vt (%)')
disp( (delta_v_BIP(1,1)-delta_v_COG(1,1))/delta_v_COG(1,1)*100 )


