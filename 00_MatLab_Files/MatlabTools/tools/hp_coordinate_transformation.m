function co = hp_coordinate_transformation(az, el, type)
% transform different coordinate conventions to matlab default
% if az containes values bigger than 2pi, it is assumed that the unit of az
% is degree. In this case it is converted to radians
%
% type
% 1: Matlab default
% 0 <= az < 360 (0=source in front 90 = source to the left of listener)
% 90 <= el <= -90 (90=north pole, 0=front, -90=south pole)
% positive x is at (0|90) (az|el)
% positive y is at (90|90)
% positive z is at (0|0);
%
% 2: Mathematical
% 0 <= az < 360 (0=source in front 90 = source to the left of listener)
% 0 <= el <= 180 (0=north pole, 90=front, 180=south pole)
% positive x is at (0|90) (az|el)
% positive y is at (90|90)
% positive z is at (0|0);

if max(max(az > 2*pi))
    az = deg2rad(az);
    el = deg2rad(el);
end

switch type
    case 1 % matlab default, see doc sph2cart
        % coordinate transformation
        co.az = az;
        co.el = el;
        % definition of axis normals for balloon plots
        co.x = [1 0 0]';
        co.y = [0 1 0]';
        co.z = [0 0 1]';
        % x and y axis for spherical planar plots
        co.xtick      = 1:45:361;
        co.xticklabel = 0:45:360;
        co.ytick      = 1:45:181;
        co.yticklabel = [90:-45:0 -45 -90];
    case 2 % mathematical correct
        
        % coordinate transformation
        co.az = az;
        id1     = el > pi/2;
        id2     = el <= pi/2;
        el(id1) = -(el(id1) - pi/2);
        el(id2) = abs(el(id2) - pi/2);
        co.el   = el;
        % definition of axis normals for balloon plots
        co.x = [1 0 0]';
        co.y = [0 1 0]';
        co.z = [0 0 1]';
        % x and y axis for spherical planar plots
        co.xtick      = 1:45:361;
        co.xticklabel = 0:45:360;
        co.ytick      = 1:45:181;
        co.yticklabel = 0:45:180;
end