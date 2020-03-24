% ons = onset_detect(x, US, onset_threshold)
%
% INPUT:
% x               - data [samples x channels]
% US              - upsampling faktor US = {2,3,4...}
% onset_threshold - threshold for onset detection relative to maximum of
%                   applied separately for each column of data
%
% OUTPUT
% ons             - detected onsets in samples referring to orriginal
%                   samplerate before upsampling
function ons = onset_detect(x, US, onset_threshold)

ons = zeros(1, size(x, 2));

for r=1:size(x,2)
    % upsampling
    x_us=interp(x(:,r),US);
    % Maximum of x
    [MAX pos] = max(abs(x_us));
    % find position were where onset threshold is reached
    ons(r) = find(abs(x_us(1:pos)) >= MAX*10^(-abs(onset_threshold)/20), 1, 'first');
end

ons = ons/US + 1;