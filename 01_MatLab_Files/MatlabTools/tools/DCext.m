% Applies dc-correcting colva window on input
% [Müller S (1999) Digitale Signalverarbeitung für Lautsprecher]
%
% syntax: h = DCext(h,w)
%
% h         - multichannel impulse response to be corrected
% w         - window function to be used
% 
% Zora Schaerer, zora.schaerer@tu-berlin.de, 05.08.2008

function h = DCext(h,w)

N = length(h);
channels = size(h,2);

w_sum = sum(w);
h_sum = sum(h);
div = h_sum./w_sum;
mult = zeros(N,channels);
for i = 1:channels 
    mult(:,i) = div(i).*w;
    h(:,i) = h(:,i) - mult(:,i);
end