% Example for reading the SOFA data:
%
% 1. Download the SOFA API for Matlan/Octave:
%    http://www.sofaconventions.org/mediawiki/index.php/Software_and_APIs
% 2. Go to the folder where you unpacked the data and run
%    SOFAstart.m
% 3. load
FABIAN = SOFAload('FABIANmeasuredHRIRs.sofa');

% 4. get left, and right ear hrirs
l = shiftdim(FABIAN.Data.IR(:,1,:), 2);
r = shiftdim(FABIAN.Data.IR(:,2,:), 2);

% 5. get source azimuth, and elevation
az = FABIAN.SourcePosition(:,1);
el = FABIAN.SourcePosition(:,2);

% Azi und Ele wählen
Azi = 41.5;
Ele = 12;


% 6. ... do whatever you feel like, for example plot the 90 degree HRTF on the
% horizontal plane
id = az==Azi & el==Ele;



% Kontrolle ob Werte von Az und El vorhanden sind
if(isempty(find(l(:,id))))
    disp('Azimuth / Elevation nicht kompatibel')
    break
end;


hp(l(:,id), 's2d', 'c', 'b')
hp(r(:,id), 's2d', 'c', 'r', 'dr', [-40 20])