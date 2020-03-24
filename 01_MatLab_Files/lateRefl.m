%% Late Reflections
%
% Creates a stochastic binaural reverb to use in ISM.
%
%%
function [x_binCoh,y_binCoh,h] = lateRefl(fs,t,delta,delay,plots)


%%
% Schallgeschwindigkeit
c = 343; %m/s
%
%% Decay of RIR 
%
% h(t) = h0 * e^-dt ; delta = (c*A)/(8*V)
h0 = 1;
h=h0 * exp(-delta*t);
%
%% Plot
if plots == true
figure ('Name','Room Impulse Response Decay','units','normalized','position',[.1 .1 .8 .8])
hold all;
subplot(2,1,1);
plot(t,h);
subplot(2,1,2);
plot(t,10*log10(h));
ylabel('dB');
end;
%
%% Rauschsignale
%
x = randn(length(t),1);
y = randn(length(t),1);
x = x/max(abs(x));
y = y/max(abs(y));
%
%% Kontrollplot
if plots == true
figure ('Name','Randn Signals','units','normalized','position',[.1 .1 .8 .8])
hold all;
subplot(2,1,1);
plot(t,x,'b');
title('Links');
subplot(2,1,2);
plot(t,y,'r');
title('Rechts');
end;
%
%
%% Filter
%
% Frequenzvektor
w = 2*pi*linspace(0,20e3,length(t));
%
w0 = 2*pi*550; %Hz
w1 = 2*pi*2700; %Hz
%
gamma = zeros(length(t),1);
for i=1:length(t)
gamma(i) = sinc(w(i)/w0) * max (0, 1-(w(i)/w1));
end;
gamma_dB = 20*log10(abs(gamma));
%
%% Kontrollplot
if plots == true
figure('NumberTitle','off','Name','Filterfunktion','units','normalized','position',[.1 .1 .8 .8])
plot(w,gamma_dB);
xlim([0 20e3]);
ylim([-80 0]);
end;
%
%
%% example code for generating an almost prefectly reconstructing filterbank
% using the auditory modeling toolbox.
flow = 50;                  % Lowest center frequency in Hz;
basef = 1000;               % Base center frequency in Hz;
fhigh = 20000;              % Highest center frequency in Hz;
filters_per_ERBaud = 1.5;   % Filterband density on ERB scale;
filter_order = 4;           % Filter order;
bw_factor = 1.0;            % Bandwidth factor;
desired_delay = delay;      % Desired delay in seconds;
% Construct new analyzer object;
analyzer = gfb_analyzer_new(fs,flow,basef,fhigh,filters_per_ERBaud,filter_order,bw_factor);
% Build synthesizer for an analysis-synthesis delay of desired_delay in seconds.
synthesizer = gfb_synthesizer_new(analyzer, desired_delay);
% Impulse signal;
impulse = x;
% Filter signal;
[analyzed_impulse, analyzer] = gfb_analyzer_process(analyzer, impulse);
% Resynthesize filtered impulse response from above.
[x_filtered, synthesizer] = gfb_synthesizer_process(synthesizer, analyzed_impulse);
%
impulse = [y];
% Filter signal;
[analyzed_impulse, analyzer] = gfb_analyzer_process(analyzer, impulse);
% Resynthesize filtered impulse response from above.
[y_filtered, synthesizer] = gfb_synthesizer_process(synthesizer, analyzed_impulse);
%
%% Kontrollplot
if plots == true
figure('Name','AMT Filter','units','normalized','position',[.1 .1 .8 .8])
hold all;
subplot(2,1,1);
plot(t,x,':b');
title('Links');
plot(t,x_filtered);

subplot(2,1,2);
plot(t,y_filtered);
title('Rechts');
end;
%
%
%% Binaural Coherence
%
H_beta=zeros(length(gamma),1);
for i=1:length(gamma)
    if gamma(i)>=0
        H_beta(i) = sqrt(0.5 * (1 - sqrt(1-(gamma(i))^2)));
    end;
    if gamma(i)<0
        H_beta(i) = - sqrt(0.5 * (1 - sqrt(1-(gamma(i))^2)));
    end;
end;
%
H_alpha = sqrt(1-H_beta.^2);
%
%
%% resynthesized impulse -> decay
%
x_f_dec = h .* x_filtered;
y_f_dec = h .* y_filtered;
%
%% -> F Domain 
%
L = length(t);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
X = fft(x_f_dec,NFFT)/L;
Y = fft(y_f_dec,NFFT)/L;
f = fs/2*linspace(0,1,NFFT/2+1);
%
% Plot single-sided amplitude spectrum.
if plots == true 
figure('Name','Spektrum','units','normalized','position',[.1 .1 .8 .8])
hold all;
subplot(2,1,1);
plot(f,2*abs(X(1:NFFT/2+1))) 
title('Links');
subplot(2,1,2);
plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Rechts')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
end;
%
%% -> * H_beta/alpha -> * TF
%
X_binCoh = X(1:NFFT/2+1) .* H_alpha(1:NFFT/2+1)' + Y(1:NFFT/2+1) .* H_beta(1:NFFT/2+1)';
Y_binCoh = Y(1:NFFT/2+1) .* H_alpha(1:NFFT/2+1)' + X(1:NFFT/2+1) .* H_beta(1:NFFT/2+1)';
%
% Kontrollplot
if plots == true
figure('Name','Spektrum Bin Coh','units','normalized','position',[.1 .1 .8 .8])
hold all;
subplot(2,1,1);
plot(w(1:NFFT/2+1),20*log10(abs(X_binCoh)));
title('Links');
subplot(2,1,2);
plot(w(1:NFFT/2+1),20*log10(abs(Y_binCoh)));
title('Rechts');
end;
%
%% Rücktransformation
x_binCoh = ifft(X_binCoh,NFFT)*L;
y_binCoh = ifft(Y_binCoh,NFFT)*L;
%
if plots == true 
figure('Name','Rücktransformierte Signale','units','normalized','position',[.1 .1 .8 .8])
subplot(2,1,1);
plot(t(1:NFFT/2+1),real(x_binCoh(1:NFFT/2+1)));
xlim([0 0.5]); %s
title('Links');
subplot(2,1,2);
plot(t(1:NFFT/2+1),real(y_binCoh(1:NFFT/2+1)));
xlim([0 0.5]); %s
title('Rechts');
end;
%
%% ======================= END OF SCRIPT ==================================
