% example code for generating an almost prefectly reconstructing filterbank
% using the auditory modeling toolbox.
close all; clear; clc
 
fs = 44100;                 % Sampling rate in Hz;
flow = 50;                  % Lowest center frequency in Hz;
basef = 1000;               % Base center frequency in Hz;
fhigh = 20000;              % Highest center frequency in Hz;
filters_per_ERBaud = 1.5;   % Filterband density on ERB scale;
filter_order = 4;           % Filter order;
bw_factor = 1.0;            % Bandwidth factor;
desired_delay = 0.004;      % Desired delay in seconds;
 
% Construct new analyzer object;
analyzer = gfb_analyzer_new(fs,flow,basef,fhigh,filters_per_ERBaud,filter_order,bw_factor);
% Build synthesizer for an analysis-synthesis delay of desired_delay in seconds.
synthesizer = gfb_synthesizer_new(analyzer, desired_delay);
% Impulse signal;
impulse = [1, zeros(1,8191)];
% Filter signal;
[analyzed_impulse, analyzer] = gfb_analyzer_process(analyzer, impulse);
% Resynthesize filtered impulse response from above.
[resynthesized_impulse, synthesizer] = gfb_synthesizer_process(synthesizer, analyzed_impulse);
 
 
subplot(2,1,1)
hp(real(analyzed_impulse)', 's2d', 'c', 'cyc', 'fs', fs)
 
subplot(2,1,2)
hp(resynthesized_impulse', 's2d', 'c', 'cyc', 'fs', fs)