function hp_draw_coordinates(co, type, marker_size, line_width)
% draws the coordinate axis in balloon plots
%
% type
% 1: lines are drawn, marker indicates positive x,y and z
%
% fabian.brinkmann@tu-berlin.de, Audio Communication Group TU Berlin,
% DFG research unit 'SEACEN', 7/2012

if ~exist('marker_size', 'var')
    % this is the size for the marker that indicates positive x,y and z
    marker_size = 10;
end

switch type
    case 1
        hold on
        % draw x axis and point for markin positive x
        plot3([-co.x(1)*co.L co.x(1)*co.L], ...
              [-co.x(2)*co.L co.x(2)*co.L], ...
              [-co.x(3)*co.L co.x(3)*co.L], 'r', 'LineWidth', line_width)
        plot3(co.x(1)*co.L, co.x(2)*co.L, co.x(3)*co.L, 'r.', 'Markersize', marker_size)
        % draw y axis and point for markin positive y
        plot3([-co.y(1)*co.L co.y(1)*co.L], ...
              [-co.y(2)*co.L co.y(2)*co.L], ...
              [-co.y(3)*co.L co.y(3)*co.L], 'g', 'LineWidth', line_width)
        plot3(co.y(1)*co.L, co.y(2)*co.L, co.y(3)*co.L, 'g.', 'Markersize', marker_size)
        % draw z axis and point for markin positive z
        plot3([-co.z(1)*co.L co.z(1)*co.L], ...
              [-co.z(2)*co.L co.z(2)*co.L], ...
              [-co.z(3)*co.L co.z(3)*co.L], 'b', 'LineWidth', line_width)
        plot3(co.z(1)*co.L, co.z(2)*co.L, co.z(3)*co.L, 'b.', 'Markersize', marker_size)
        hold off
end