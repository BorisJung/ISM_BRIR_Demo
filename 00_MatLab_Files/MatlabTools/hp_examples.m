close all; clear; clc

% load exemplary dataset: Horizontal plane HRTFs
% (resolution 10 degree, 0=source in front, 90=source at left ear)
load('hp_example_data')

x.l = left;
x.r = right;

addpath('tools')

%% ----------------------------- check all 2D plots with default parameters
% (press enter for next plot)
% NOTE: nice_ticks only looks ugly on 'g2d' plot because of the MATLAB default labeling 3 x 10^4
clc
pt = {'ir' 'etc' 'edc' 'sr' 's' 'g' 'ph' 'pu' 'toa'};
for k = 1:length(pt)
    hp(left(:,1:9:28), [pt{k} '2D'])
    nice_ticks
    input('press any enter');
    close all
end

pt = {'itd' 'ild'};
for k = 1:length(pt)
    hp(x, [pt{k} '2D'])
    nice_ticks
    input('press any enter');
    close all
end
%%
hp(left, 's6', 'az', az, 'sph_f', 16000)


%% ---------------------------------- check optional arguments for 2D plots
figure
hp(left(:,1:9:28), 's2D', 'line_m', 's', 'line_s', '--', 'line_w', 1, 'c', [1 0 0; 0 0 1]); nice_ticks
%%
figure
hp(left(:,1:9:28), 's2D', 'frac', 3, 'c', 'cyc', 'dr', [-30 20], 'x', [1000 20000], 'xu', 'Hz'); nice_ticks
%%
figure
[data, h] = hp(x, 'itd2D', 'x', 0:10:350, 'du', 'n'); nice_ticks


%% ----------------------------- check all 3D plots with default parameters
%  (close figure and press return to go to next type)
pt = {'ir' 'etc' 'edc' 'sr' 's' 'g' 'ph' 'pu'};
for k = 1:length(pt)
    hp(left, [pt{k} '3D'])
    input('press any enter');
    close all
end


%% ---------------------------------- check optional arguments for 3D plots
hp(left, 's3D', 'x', [200 20000], 'y', 0:10:350, 'dr', [-20 20], 'cr', .5, 'cm', 'copper')

