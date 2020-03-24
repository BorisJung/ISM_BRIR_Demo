%% edit_version
%

%
close all;
clear all;
clc;
%
%% Spiegelquellen Modell
%
%
pause off;
diary('ISM_log1.txt');
tic
%
%% Konstanten
% Schallgeschwindigkeit
c = 343; %m/s
% Samplingfrequenz
fs = 44.1e3; %Hz
% HRIR von FABIAN mit 44.1kHz vorhanden!
%
%% Eingabe Parameter
%
% Berücksichtigungsradius (in Anzahl Quell-Empfänger-Entfernung)
dist_lim = 25;
% Limit für L (siehe Formel 11 Script)
N = 5
% Mix Methode für Late Reflections ('hardcut','linfade','logfade')
mixmethod = 'linfade';
% Raum Maße (in m)
x_length = 8.5;
y_length = 6.5;
z_length = 3.5;
%
% Quellposition
xs = 2.5;
ys = 1.5;
zs = 1.7;
% Empfängerposition
xr = 1.5;
yr = 2;
zr = 1.7;
% Reflektionskoeffizienten
bx1 = .8;
bx2 = .6;
by1 = .6;
by2 = .6;
bz1 = .7;
bz2 = .5;
%
%% Zusammenfassen der Eingabedaten
%
beta = [bx1 bx2; by1 by2; bz1 bz2];
wall_lengths = [x_length ; y_length ; z_length];
source_pos = [xs ; ys ; zs];
receiver_pos = [xr ; yr ; zr];
source_distance = sqrt((xs-xr)^2+(ys-yr)^2+(zs-zr)^2);
%
%
%% Berechnungen aus Raumparameter
[V,A,S,alpha,delta,T60,tmp50] = acRoomPara(wall_lengths, beta);
%
%% Zeitparameter
% Dauer
dur = 1.5; % s
% Zeitvektor
t = 0:1/fs:dur;
% tmp50 Sample Index
tmp50_sample_ind = find(floor(t.*1000-tmp50)==0,1,'first');
%
%% HRIR laden
%
run SOFAstart.m;
FABIAN = SOFAload('FABIANmeasuredHRIRs.sofa');
% get left, and right ear hrirs
l = shiftdim(FABIAN.Data.IR(:,1,:), 2);
r = shiftdim(FABIAN.Data.IR(:,2,:), 2);
% get source azimuth, and elevation
az = FABIAN.SourcePosition(:,1);
el = FABIAN.SourcePosition(:,2);
%
%% U und L Vektoren erstellen
% U Vektor erstellen
U = transpose(dec2bin(0:2^3-1)-'0');
% L initialisieren
L = zeros(3,2*N^3+1);  
zaehler=0;
for i=-N:N
    for k=-N:N
        for f=-N:N
            zaehler=zaehler+1;
            L(: , zaehler)=[i  k  f];
        end;
    end;
end;
%
%
%% Spiegelquell Positionen berechnen
%
% IS Pos init
IS_pos = [0;0;0];
IS_pos_rel = [0;0;0];
a = 0;
%
%
for n=1:length(L)
for i=1:length(U)
     % Spiegelquellenposition (Formel 9 Script)
     IS_pos_temp = ((ones(3,1) - 2 * U(:,i)) .* source_pos + (2 * ones(3,1) .* L(:,n)) .* wall_lengths);   
     IS_pos_rel_temp = IS_pos_temp - receiver_pos;
     %
     % zugehörige Amplitude (Formel 14 Script)
     a_temp = ( (bx1 ^ (abs(U(1,i)-L(1,n)))) * (bx2 ^ (abs(L(1,n)))) * ...
           (by1 ^ (abs(U(2,i)-L(2,n)))) * (by2 ^ (abs(L(2,n)))) * ...
           (bz1 ^ (abs(U(3,i)-L(3,n)))) * (bz2 ^ (abs(L(3,n)))) );
     %
     % Werte der aktuellen Spiegelquellen an Array anhängen
     IS_pos = [IS_pos IS_pos_temp];
     IS_pos_rel = [IS_pos_rel IS_pos_rel_temp];
     a = [a a_temp];
     %
end;
end;
%
%
% Initialisierungswerte entfernen
IS_pos(:,1)=[];
IS_pos_rel(:,1)=[];
a(:,1)=[];
%
%
% relative IS Pos in Kugelkoordinaten für Azi und Ele (WINKEL IN RAD!!)
[IS_az, IS_el, IS_r] = cart2sph(IS_pos_rel(1,:),IS_pos_rel(2,:),IS_pos_rel(3,:));
%IS_rel_SPH = [IS_az; IS_el; IS_r]; %in rad
IS_rel_SPH_deg = [rad2deg(IS_az); rad2deg(IS_el); IS_r]; %in Grad
%
tau = IS_r/c;
%
%
clc;
disp('Quellpositionen erfolgreich berechnet!');
disp('Entferne irrelevante Quellen...');
pause
%
%% Zu weit entfernte Quellen entfernen
%
% maximale Spiegelquellenentfernung
dist_limit = dist_lim * source_distance;
%
% Indizes der entsprechenden Quellen finden
dist_index=[find(IS_r>dist_limit) find(IS_r==0)];
% (IS_r=0: Spiegelquelle im Empfängerstandort!?)
%
% Quellen und zugehörige Werte entfernen
IS_rel_SPH_deg(:,dist_index) = [];
IS_az(:,dist_index) = []; IS_el(:,dist_index) = []; IS_r(:,dist_index) = [];
IS_pos(:,dist_index) = [];
IS_pos_rel(:,dist_index) = [];
tau(:,dist_index)=[];
a(:,dist_index)=[];
%
%
%% Wertebereich der genutzten HRIR:
% az:{0...360} ; el:{-64...90}
%
% => negative Azimut Winkel konvertieren (von -180bis180 auf 0bis360) 
for i=1:size(IS_rel_SPH_deg,2)
    if IS_rel_SPH_deg(1,i) < 0
        IS_rel_SPH_deg(1,i) = IS_rel_SPH_deg(1,i) + 360;
        %disp(['Winkel konvertiert auf ' , num2str(IS_rel_SPH_deg(1,i)) ]);
    end;
    if IS_rel_SPH_deg(1,i) == 360
        IS_rel_SPH_deg(1,i) = 0;
        disp('0=360');
    end;
end;
pause
%
%
% id Vektor zum Auswählen der HRIR erstellen
id_vektor_tmp = [IS_rel_SPH_deg(1,:);IS_rel_SPH_deg(2,:)];
%
clc;
disp('Winkel der Spiegelquellen dem Raster der HRIR anpassen');
pause;
%
%% Winkelgenauigkeit
%
% Winkel der Spiegelquellen dem Raster der HRIR anpassen
id_vektor = angleFix(id_vektor_tmp,az,el);
%
% Winkel runden
id_vektor_high = angleRound(round(id_vektor_tmp),az,el,'high');
id_vektor_med = angleRound(round(id_vektor_tmp),az,el,'medium');
id_vektor_low = angleRound(round(id_vektor_tmp),az,el,'low');
id_vektor_vlow = angleRound(round(id_vektor_tmp),az,el,'vlow');
%
% 0-Vektor (alle quellen 0°/0°)
id_vektor_0 = zeros(2,length(id_vektor));
% Az=0 Vektor
id_vektor_0az = [zeros(1,length(id_vektor)) ; id_vektor(2,:) ];
%
%
% Auswahl
%id_vektor = id_vektor_high;
%
%
%% Erneute Kontrolle
for i=1:length(id_vektor)
id = az==id_vektor(1,i) & el==id_vektor(2,i);
if (isempty(find(id,1)))
  disp(['Wertepaar ' ,  num2str(i) , ' der Azimut/Elevationswinkel nicht vorhanden']);
  disp(['Abbruch...']);
  return;
end;
if (isempty(find(id,1)) == 0)
    disp(['Wertepaar ' , num2str(i) , ' vorhanden']);
end;
end;
%
%
disp('Winkel Endkontrolle abgeschlossen...');
pause;
%
%% Impulse erstellen
% Zeitpunkte werden gerundet! (max Fehler +/- 1/fs )
%
gew_pulse = zeros(length(tau),length(t));
for i = 1:size(gew_pulse,1)
ind = knnsearch(t',tau(i));
gew_pulse(i,ind) = a(i);
end;
%    
%
clc;
disp('Impulse erstellt');
pause;
%
%% IS Pulse und HRIR falten
laenge = size(conv(gew_pulse(1,:)',l(:,id)),1);
%
L_BRIR = zeros(length(id_vektor),laenge);
R_BRIR = zeros(length(id_vektor),laenge);
%
% Puls jeder Quelle mit entsprechender HRIR falten
for i=1:size(gew_pulse,1)
id = az==id_vektor(1,i) & el==id_vektor(2,i);
L_BRIR(i,:) = conv(gew_pulse(i,:)',l(:,id));
R_BRIR(i,:) = conv(gew_pulse(i,:)',r(:,id));
end;
% Faltungsprodukte/Pulse aufsummieren
L_BRIR_SUM = sum(L_BRIR);
R_BRIR_SUM = sum(R_BRIR);
%
%
%clc;
disp('Impulse gefaltet und aufsummiert...');
pause;
%
%
%% Late Reflections
%
% Filter delay -    "the desired group delay of the analysis-synthesis
%                    system in seconds.  Greater delays result in better
%                    output signal quality.  Minimum delay is (1 /
%                    analyzer.fs)"
delay = 10/fs; %tmp50?
%
[x_binCoh,y_binCoh,h] = lateRefl(fs,t,delta,delay,false);
%
clc;
disp('Late Reflections berechnet...');
pause;
%
%% Decay
%L_BRIR_SUM = L_BRIR_SUM(1:length(t)).*h;
%R_BRIR_SUM = R_BRIR_SUM(1:length(t)).*h;
%
%
%% Normalisieren
%
normFac=max([abs(L_BRIR_SUM(:)) ; abs(R_BRIR_SUM(:))]);
L_BRIR_SUM_DEC = L_BRIR_SUM/normFac;
R_BRIR_SUM_DEC = R_BRIR_SUM/normFac;

normFac2=max([abs(x_binCoh(:)) ; abs(y_binCoh(:))]);
x_binCoh=x_binCoh/normFac2;
y_binCoh=y_binCoh/normFac2;
%
% Vorlauf
x_binCoh=[zeros(1,knnsearch(t',min(tau))) x_binCoh];
y_binCoh=[zeros(1,knnsearch(t',min(tau))) y_binCoh];
%
%
%% Mix ISM - Stochastic Reverb
%
if mixmethod == 'hardcut';
% Hard Cut
BRIR_L = [(L_BRIR_SUM_DEC(1:tmp50_sample_ind)+0.1*x_binCoh(1:tmp50_sample_ind)) ...
           x_binCoh(tmp50_sample_ind+1:length(t))]';
BRIR_R = [(R_BRIR_SUM_DEC(1:tmp50_sample_ind)+0.1*y_binCoh(1:tmp50_sample_ind)) ...
           y_binCoh(tmp50_sample_ind+1:length(t))]';
BRIR = [BRIR_L BRIR_R];
end;
%
%
% Fades
nsamples_fade = 2e2; % samples ; 1e3 samples = 22.7ms
fade_start = tmp50_sample_ind;
fade_end = fade_start + nsamples_fade -1;
early_factor = 0.05; % Amplitude factor for stochastic reverb in early part of BRIR
%
%
% Linear crossfade
if mixmethod == 'linfade';
fadefactor = linspace(0,1,nsamples_fade);
%
BRIR_L = [ L_BRIR_SUM_DEC(1:fade_start-1) + early_factor .* x_binCoh(1:fade_start-1) ...
           ... 
           (1-fadefactor) .* L_BRIR_SUM_DEC(fade_start:fade_end)...
           + fadefactor .* x_binCoh(fade_start:fade_end) ...
           ...
           x_binCoh(fade_end+1:length(t)) ]';
%
BRIR_R = [ R_BRIR_SUM_DEC(1:fade_start-1) + early_factor .* y_binCoh(1:fade_start-1) ...
            ...
           (1-fadefactor) .* R_BRIR_SUM_DEC(fade_start:fade_end)...
           + fadefactor .* y_binCoh(fade_start:fade_end) ...
            ...
           y_binCoh(fade_end+1:length(t)) ]';             
%       
BRIR = [BRIR_L BRIR_R];
end;
%
%
% Logarithmic crossfade
if mixmethod == 'logfade';
fadefactor = logspace(0,1,nsamples_fade);
%
BRIR_L = [ L_BRIR_SUM_DEC(1:fade_start-1)...
           ... 
           (1-fadefactor) .* L_BRIR_SUM_DEC(fade_start:fade_end)...
           + fadefactor .* x_binCoh(fade_start:fade_end) ...
           ...
           x_binCoh(fade_end+1:length(t)) ]';
%
BRIR_R = [ R_BRIR_SUM_DEC(1:fade_start-1)...
            ...
           (1-fadefactor) .* R_BRIR_SUM_DEC(fade_start:fade_end)...
           + fadefactor .* y_binCoh(fade_start:fade_end) ...
            ...
           y_binCoh(fade_end+1:length(t)) ]';       
%
BRIR = [BRIR_L BRIR_R];
end;
%
% backup
complexBRIR = BRIR;
%
BRIR = real(BRIR);
%
% Using the 2D vector BRIR and the function mono2binaural.m a binaural 
% can be created at this point.
% 
% Binaural_Signal =  =
% mono2binaural(filename,fs,BRIR,write_audio,play_audio);
%
%
clc;
disp('Signal erstellt!');
disp('Bereit zum plotten...');
pause;
%
%
%% ----------------------- PLOTS ------------------------------------------
%
%% Plot BRIRs roh
%
%
xlimit = 2.5*tmp50/1000;%500e-3; %ms
ylimit = 1.2*max(max(abs(L_BRIR_SUM),abs(R_BRIR_SUM)));
%
figure('NumberTitle','off','Name','BRIR - roh','units','normalized','position',[.1 .1 .8 .8])
hold all;
%
subplot(2,1,1);
hold all;
plot(t,L_BRIR_SUM(1:length(t)),'b');
title('Linkes Ohr');
xlim([0 xlimit]);
ylim([-ylimit ylimit]);
line([tmp50/1000 tmp50/1000],[-ylimit ylimit],'LineWidth',1,'LineStyle',':','Color',[0 0 0])
subplot(2,1,2);
hold all;
plot(t,R_BRIR_SUM(1:length(t)),'r');
title('Rechtes Ohr');
xlim([0.0 xlimit]);
ylim([-ylimit ylimit]);
line([tmp50/1000 tmp50/1000],[-ylimit ylimit],'LineWidth',1,'LineStyle',':','Color',[0 0 0])

%% BRIR & lateRefl getrennt
figure('NumberTitle','off','Name','BRIR mit Decay und stochastic binaural reverb','units','normalized','position',[.1 .1 .8 .8])
hold all;
subplot(2,1,1);
hold all;
plot(t,L_BRIR_SUM_DEC(1:length(t)),'b');
plot(t,real(x_binCoh(1:length(t))),'c');
title('Linkes Ohr');
xlim([0 xlimit]);
ylim([-ylimit ylimit]);
line([tmp50/1000 tmp50/1000],[-ylimit ylimit],'LineWidth',1,'LineStyle',':','Color',[0 0 0])
subplot(2,1,2);
hold all;
plot(t,R_BRIR_SUM_DEC(1:length(t)),'r');
plot(t,real(y_binCoh(1:length(t))),'m');
title('Rechtes Ohr');
xlim([0.0 xlimit]);
ylim([-ylimit ylimit]);
line([tmp50/1000 tmp50/1000],[-ylimit ylimit],'LineWidth',1,'LineStyle',':','Color',[0 0 0])
%
%
%
%% Plot kombinierte BRIR
%
figure('NumberTitle','off','Name','Final BRIR','units','normalized','position',[.1 .1 .8 .8])
subplot(2,1,1);
hold all;
plot(t,BRIR(:,1),'b');
line([tmp50/1000 tmp50/1000],[-ylimit ylimit],'LineWidth',1,'LineStyle',':','Color',[0 0 0])
xlim([0.0 xlimit]);
ylim([-ylimit ylimit]);
title('Links');
subplot(2,1,2);
hold all;
plot(t,BRIR(:,2),'r');
line([tmp50/1000 tmp50/1000],[-ylimit ylimit],'LineWidth',1,'LineStyle',':','Color',[0 0 0])
xlim([0.0 xlimit]);
ylim([-ylimit ylimit]);
title('Rechts');
%
%% Kontrollplot Impulse
figure('Name','Impulse','units','normalized','position',[.1 .1 .8 .8])
hold all;
plot(t,gew_pulse)
line([tmp50/1e3 tmp50/1e3],[-0.1 1.25],'LineStyle',':','Color',[0 0 0]);
xlim([-0.001 (tmp50+50)/1000])
ylim([-0.1 1.25])
%
%
%% Spiegelquellen Positionen grafisch
%
figure('Name','Spiegelquellen Modell','units','normalized','position',[.1 .1 .8 .8])
hold all;
xlim([1.25*(receiver_pos(1)-dist_limit) 1.25*(receiver_pos(1)+dist_limit)]);
axis equal;
%Spiegelräume erster Ordnung
rectangle('Position',[wall_lengths(1) 0 wall_lengths(1) wall_lengths(2)] , 'FaceColor',[.5 .5 .5]);
rectangle('Position',[-wall_lengths(1) 0 wall_lengths(1) wall_lengths(2)] , 'FaceColor',[.5 .5 .5]);
rectangle('Position',[0 wall_lengths(2) wall_lengths(1) wall_lengths(2)] , 'FaceColor',[.5 .5 .5]);
rectangle('Position',[0 -wall_lengths(2) wall_lengths(1) wall_lengths(2)] , 'FaceColor',[.5 .5 .5]);
%Spiegelräume zweiter Ordnung
rectangle('Position',[2*wall_lengths(1) 0 wall_lengths(1) wall_lengths(2)] , 'FaceColor',[.7 .7 .7]);
rectangle('Position',[0 2*wall_lengths(2) wall_lengths(1) wall_lengths(2)] , 'FaceColor',[.7 .7 .7]);
rectangle('Position',[1*wall_lengths(1) 1*wall_lengths(2) wall_lengths(1) wall_lengths(2)] , 'FaceColor',[.7 .7 .7]);
rectangle('Position',[1*wall_lengths(1) -1*wall_lengths(2) wall_lengths(1) wall_lengths(2)] , 'FaceColor',[.7 .7 .7]);
rectangle('Position',[0 -2*wall_lengths(2) wall_lengths(1) wall_lengths(2)] , 'FaceColor',[.7 .7 .7]);
rectangle('Position',[-1*wall_lengths(1) -1*wall_lengths(2) wall_lengths(1) wall_lengths(2)] , 'FaceColor',[.7 .7 .7]);
rectangle('Position',[-2*wall_lengths(1) -0*wall_lengths(2) wall_lengths(1) wall_lengths(2)] , 'FaceColor',[.7 .7 .7]);
rectangle('Position',[-1*wall_lengths(1) 1*wall_lengths(2) wall_lengths(1) wall_lengths(2)] , 'FaceColor',[.7 .7 .7]);
% Berücksichtungsradius einzeichnen
viscircles([receiver_pos(1) receiver_pos(2)],dist_limit,'LineWidth',1);%,'LineStyle',':');
% Symbol
plot(receiver_pos(1)+dist_limit,receiver_pos(2),'r','LineWidth',2);
%
% Wände/Räume zeichnen
%Raum
rectangle('Position',[0 0 wall_lengths(1) wall_lengths(2)],'FaceColor',[.9 .9 .9],'EdgeColor',[.5 0 0]);
%
% Quelle und Empfänger
plot(source_pos(1),source_pos(2),'pb');
plot(receiver_pos(1),receiver_pos(2),'>r');
%
for i=1:length(IS_pos)
plot(IS_pos(1,i),IS_pos(2,i),'+c');
end;
%
xlabel('m','FontSize',20);
ylabel('m','FontSize',20);
%title('Spiegelquellenmodell - Positionen','FontSize',26);
legend('Berücksichtigungsradius','Quelle','Empfaenger','Spiegelquellen');
%
% Ursprung einzeichen
plot(0,0,'xr');

toc
%
%
%% ======================= END OF SCRIPT ==================================



