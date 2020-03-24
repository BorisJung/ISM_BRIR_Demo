function [single_sided is_even] = both2single_sided_spectrum(both_sided)
% single_sided = both2single_sided_spectrum(both_sided)
%
% INPUT
% both-sided spectrum, single- or multi-channel
% [1 column=1 channels x 1 rows=1 frequency bin]
%
% OUTPUT
% single-sided spectrum.
%
% F. Brinkmann, Audio Communication Group, TU Berlin, 04/2013
NFFT = size(both_sided,1);

if ~mod(NFFT,2)
    single_sided = both_sided(1:NFFT/2+1,:);
else
    single_sided = both_sided(1:ceil(NFFT/2),:);
end

is_even = 1-mod(size(both_sided,1),2);