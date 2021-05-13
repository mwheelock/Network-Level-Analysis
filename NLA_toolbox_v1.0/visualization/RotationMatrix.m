function RotMat=RotationMatrix(direction,theta)

% This function generates a rotation matrix for the direction (x,y,z) given
% a angle (in radians!).

RotMat=zeros(3);

switch direction
    case 'x'
        RotMat=[1 0 0;...
            0 cos(theta) -sin(theta);...
            0 sin(theta) cos(theta)];
    case 'y'
        RotMat=[cos(theta) 0 sin(theta);...
            0 1 0;...
            -sin(theta) 0 cos(theta)];
    case 'z'
        RotMat=[cos(theta) -sin(theta) 0;...
            sin(theta) cos(theta) 0;...
            0 0 1];
end