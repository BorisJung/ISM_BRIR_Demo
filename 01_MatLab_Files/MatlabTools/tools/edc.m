% y = edc(x, norm)
%
% calculates normalized (norm = 1, default) or unnormalized (norm = 0)
% energy decay curve in dB of input data x [samples x channels]
%
% F. Brinkmann, Audio Communicatin Group TU-Berlin
% fabian.brinkmann@tu-berlin.de, 04/2012

function y = edc(x, norm)

if ~exist('norm', 'var')
    norm = 1;
end

y = flipud(cumtrapz(flipud(x.^2)));

if norm
    y = y./repmat(trapz(x.^2), size(x,1),1);
end

y = 10*log10(y);