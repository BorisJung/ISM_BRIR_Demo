%% Perceptual Mixing Time
%
% Calculates the perceptual mixing time tmp50 = 20*V/S+12[ms]
% from room volume and surface
%
%
function tmp50 = percMixT50(V,S)

tmp50 = 20 * V / S + 12;

