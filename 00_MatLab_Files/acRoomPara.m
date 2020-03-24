%% Acoustical Room Parameters
%
% Calculates acoustical room parameters from given wall lengths and
% reflexion coefficients.
%
% V: Room volume
% A: Total area
% S: Total surface (area * absorption coefficient)
% alpha: Absorption coefficients
% T60: Reverbaration time T60
% tmp50: Perceptual mixing time tmp50
%
%% 
function [V,A,S,alpha,delta,t60,tmp50] = acRoomPara(Lengths, beta)
%
% Input Check
if nargin < 2
    disp('Fehler! nargin');
end;
%
if size(Lengths) ~= [3,1]
    disp('Fehler! size(Lengths)');
end;
%
if size(beta) ~= [3,2]
    disp('Fehler! size(beta)');
end;
%
%% Fix Values
% Speed of sound
c = 343; %[m/s]
%
%% Volumen
V = Lengths(1) * Lengths(2) * Lengths(3);
%% Flächen
xyArea = Lengths(1) * Lengths(2);
yzArea = Lengths(2) * Lengths(3);
xzArea = Lengths(1) * Lengths(3);
% Gesamtfläche
A = 2*xyArea+2*yzArea+2*xzArea;
%
%% Absorptionskoeffizieten
% ß=sqrt(1-a)
% a = 1-b^2
% alpha = [alphax1 alphax2 ; alphay1 alphay2 ; alphaz1 alphaz2]
% 
for i = 1:size(beta,1)
    for j=1:size(beta,2)
        alpha(i,j) = 1-beta(i,j)^2;
    end;
end;
%
%% Absorptionsfläche
%
S = xyArea * alpha(1,1) + ...
    xyArea * alpha(1,2) + ...
    yzArea * alpha(2,1) + ...
    yzArea * alpha(2,2) + ...
    xzArea * alpha(3,1) + ...
    xzArea * alpha(3,2);
%
delta = (c*A) / (8*V);
%
%% Reverbaration Time
%
t60 = T60(V,S);
%
%% Perceptual mixing time
tmp50 = percMixT50(V,S);
%
