%% angleFix
%
% Fixes angle values resulting from ISM calculations to fit certain
% criteria matching stepsizes of the HRIR.
%
%
function [id_vektor_fixed] = angleFix(id_vektor,az,el)
%
if nargin ~= 3
    disp('Wrong number of input arguments');
    return;
end;

if size(id_vektor,1)~=2
    if size(id_vektor,2)~=2
        disp('Wrong input vector format');
        return;
    end;
    id_vektor = id_vektor';
end;
%
%
% Winkelwerte runden (max FEHLER +/- 0.5°)
id_vektor = round(id_vektor);
%
%
% Kontrolle ob Werte von Az und El vorhanden sind
for i=1:length(id_vektor)
id = az==id_vektor(1,i) & el==id_vektor(2,i);
if (isempty(find(id,1)))
  disp(['Wertepaar ' ,  num2str(i) , ' der Azimut/Elevationswinkel nicht vorhanden']);
end;
if (isempty(find(id,1)) == 0)
    disp(['Wertepaar ' , num2str(i) , ' vorhanden']);
end;
end;
%
%
%% Bearbeitung der Azi-Ele-Paare um eine entsprechende HRIR zu finden
%
% Winkelkorrektur
%
for i=1:length(id_vektor)
    if mod(id_vektor(2,i),2)~=0
        id_vektor(2,i) = id_vektor(2,i)+1;
    end;
    switch id_vektor(2,i)
        case 90
            id_vektor(1,i) = 0;
            disp('Bei 90° Elevation wird kein Azimutwinkel mehr berücksichtigt');
            
        case 88
            if mod(id_vektor(1,i),45) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/45) * 45;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;
        
        case {86 , 84 } 
            if mod(id_vektor(1,i),18) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/18) * 18;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;
        
        case {82 , 80}
            if mod(id_vektor(1,i),10) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/10) * 10;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;
            
        case {78}
            if mod(id_vektor(1,i),9) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/9) * 9;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;
            
        case {76 , 74 , 72}
            if mod(id_vektor(1,i),6) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/6) * 6;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
            
        case {70 , 68}
            if mod(id_vektor(1,i),5) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/5) * 5;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
            
        case {66 , 64 , 62 , 60 , 58 , 56 , 54 , 52 , 50}
            if mod(id_vektor(1,i),3) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/3) * 3;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
        
        case num2cell(48:-2:-48)
            
%             {48 , 46 , 44 , 42 , 40 , 38 , 36 , 34 , 32 , 30 , ...
%                 28 , 26 , 24 , 22 , 20 , 18 , 16 , 14 , 12 , 10 , ...
%                 8 , 6 , 4 , 2 , 0 , -2 , -4 , -6 , -8 , -10 , ...
%                 -12 , -14 , -16 , -18 , -20 , -22 , -24 , -26 , -28 , -30 , ...
%                 -32 , -34 , -36 , -38 , -40 , -42 , -44 , -46 , -48}
%             
            
            if mod(id_vektor(1,i),2) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/2) * 2;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
        
        case {-50 , -52 , -54 , -56 , -58 , -60 , -62 , -64}
            if mod(id_vektor(1,i),3) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/3) * 3;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
        
            
        case num2cell(-64:-2:-90)
            id_vektor(2,i) = -64;
            if mod(id_vektor(1,i),3) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/3) * 3;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
            
            
            
        otherwise
            disp(['Elevationswinkel Fehler bei Wertepaar ' , num2str(i)]);
            return;
    end;
    
end;
%
% Erneute Kontrolle ob Werte von Az und El vorhanden sind
for i=1:length(id_vektor)
id = az==id_vektor(1,i) & el==id_vektor(2,i);
if (isempty(find(id,1)))
  disp(['Wertepaar ' ,  num2str(i) , ' der Azimut/Elevationswinkel nicht vorhanden']);
  disp('Winkelkorrektur fehlgeschlagen!');
  return;
end;
if (isempty(find(id,1)) == 0)
    disp(['Wertepaar ' , num2str(i) , ' vorhanden']);
end;
end;
%
%
% write back value
id_vektor_fixed = id_vektor;
