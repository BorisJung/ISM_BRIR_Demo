%% angleRound
%
% converts angles to a certain grid
%
% angle_resolution:
%
%       'high': only changes 2 degree grid to 4 degree grid
%       'medium': 10 degree steps for most elevations
%       'low': 30 degree steps for most elevations
%
%
%

function [id_vektor_fixed] = angleRound(id_vektor,az,el,angle_resolution)
%
if nargin ~= 4
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
%% Bearbeitung der Azi-Ele-Paare 
%
% Winkelkorrektur
%
for i=1:length(id_vektor)
    if mod(id_vektor(2,i),2)~=0
        id_vektor(2,i) = id_vektor(2,i)+1;
    end;
    id_vektor_fixed = id_vektor;
    %% HIGH
    if strcmp(angle_resolution,'high');
        disp('angleRound high');
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
            if mod(id_vektor(1,i),4) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/4) * 4;
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
    
    %% MEDIUM
    if strcmp(angle_resolution,'medium');
        disp('angleRound medium');
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
            if mod(id_vektor(1,i),12) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/12) * 12;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
            
        case {70 , 68}
            if mod(id_vektor(1,i),10) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/10) * 10;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
            
        case {66 , 64 , 62 , 60 , 58 , 56 , 54 , 52 , 50}
            if mod(id_vektor(1,i),9) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/9) * 9;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
        
        case num2cell(48:-2:-48)
            if mod(id_vektor(1,i),10) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/10) * 10;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
        
        case {-50 , -52 , -54 , -56 , -58 , -60 , -62 , -64}
            if mod(id_vektor(1,i),9) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/9) * 9;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
        
            
        case num2cell(-64:-2:-90)
            id_vektor(2,i) = -64;
            if mod(id_vektor(1,i),9) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/9) * 9;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
            
            
            
        otherwise
            disp(['Elevationswinkel Fehler bei Wertepaar ' , num2str(i)]);
            return;
    end;
    end;

%% LOW
    if strcmp(angle_resolution,'low');
        disp('angleRound low');
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
            if mod(id_vektor(1,i),20) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/20) * 20;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;
            
        case {78}
            if mod(id_vektor(1,i),18) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/18) * 18;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;
            
        case {76 , 74 , 72}
            if mod(id_vektor(1,i),24) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/24) * 24;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
            
        case {70 , 68}
            if mod(id_vektor(1,i),20) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/20) * 20;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
            
        case {66 , 64 , 62 , 60 , 58 , 56 , 54 , 52 , 50}
            if mod(id_vektor(1,i),18) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/18) * 18;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
        
        case num2cell(48:-2:-48)
            if mod(id_vektor(1,i),20) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/20) * 20;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
        
        case {-50 , -52 , -54 , -56 , -58 , -60 , -62 , -64}
            if mod(id_vektor(1,i),18) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/18) * 18;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
        
            
        case num2cell(-64:-2:-90)
            id_vektor(2,i) = -64;
            if mod(id_vektor(1,i),18) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/18) * 18;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    

            
        otherwise
            disp(['Elevationswinkel Fehler bei Wertepaar ' , num2str(i)]);
            return;
    end;
    end;

    
    
    %% VERY LOW
    if strcmp(angle_resolution,'vlow');
        disp('angleRound high');
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
            if mod(id_vektor(1,i),36) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/36) * 36;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;
        
        case {82 , 80}
            if mod(id_vektor(1,i),30) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/30) * 30;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;
            
        case {78}
            if mod(id_vektor(1,i),36) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/36) * 36;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;
            
        case {76 , 74 , 72}
            if mod(id_vektor(1,i),30) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/30) * 30;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
            
        case {70 , 68}
            if mod(id_vektor(1,i),30) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/30) * 30;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
            
        case {66 , 64 , 62 , 60 , 58 , 56 , 54 , 52 , 50}
            if mod(id_vektor(1,i),30) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/30) * 30;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
        
        case num2cell(48:-2:-48)
            if mod(id_vektor(1,i),30) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/30) * 30;
            disp(['Azimuth geändert auf ' , num2str(id_vektor(1,i)) ]);
            end;    
        
        case {-50 , -52 , -54 , -56 , -58 , -60 , -62 , -64}
            if mod(id_vektor(1,i),30) ~= 0
            id_vektor(1,i) = round(id_vektor(1,i)/30) * 30;
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
    
    
    
    
    if id_vektor(1,i) == 360
        id_vektor(1,i) = 0;
        disp('0=    360');
    end;
    
    
    
end;

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
%
disp('angleRound finished successfully');
