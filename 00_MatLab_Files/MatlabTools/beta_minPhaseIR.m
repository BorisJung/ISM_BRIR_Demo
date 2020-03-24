% Example for generating a minimum phase impulse response from
% reflection coefficients given in octave frequencies.
%
% fabian.brinkmann@tu-berlin.de, Audio Communication Group, TU Berlin
close all; clear; clc

% beta and octave frequencies
fOct  = [125 250 500 1000 2000 4000 8000 16000]';
beta  = sqrt(1-[.2 .23 .34 .36   .5   .55   .7   .8]');

% length of impulse response, sampling frequency and frequency vector
N    = 2^6;
fs   = 44100;
f    = (0:fs/N:fs/2)';

% interpolate to desired frequencies of real, one-sided spectrum
beta_i = interp1(fOct, beta, f, 'nearest','extrap');
% octave smooth
beta_i = fract_oct_smooth(beta_i, 'welti', fs, 1);
% reconstruct both-sided spectrum
beta_i = single2both_sided_spectrum(beta_i, N);
% get zero phase impulse response
beta_i = ifft(beta_i, 'symmetric');
% get mininum phase impulse response
beta_i = phase_manipulation(beta_i, fs, 'min', 1, 0);

% plot result
subplot(1,2,1)
hp(beta_i, 'ir2d')
subplot(1,2,2)
hp(beta_i, 's2d', 'du', 'lin', 'dr', [0 1])