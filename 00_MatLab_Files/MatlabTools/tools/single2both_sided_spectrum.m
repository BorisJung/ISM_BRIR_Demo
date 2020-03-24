function both_sided = single2both_sided_spectrum(single_sided, is_even)
% both_sided = single2both_sided_spectrum(single_sided, is_even)
%
% INPUT
% single-sided   - single sided spectrum, single- or multi-channel
%                  [1 column=1 channels x 1 rows=1 frequency bin]
% is_even        - 1 if both sided spectrum had even number of taps.
%                  if is_even>1, it denotes the number of samples of the
%                  both sided spectrum (default = 1)
%
% OUTPUT
% both-sided spectrum.
%
% F. Brinkmann, Audio Communication Group, TU Berlin, 04/2013
if ~exist('is_even', 'var')
    is_even = 1;
end

if is_even>1
    is_even = 1-mod(is_even,2);
end

N = size(single_sided,1);

if is_even
    both_sided = [single_sided; flipud(conj(single_sided(2:N-1,:)))];
else
    both_sided = [single_sided; flipud(conj(single_sided(2:N,:)))];
end