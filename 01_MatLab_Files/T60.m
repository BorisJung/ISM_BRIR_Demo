%% Reverbaration Time
%
% Calculates shoebox rooms reverberation time T60 from Volume and Surface
%
function t60 = T60(V,S) 
t60 = 0.161 * V / S;
