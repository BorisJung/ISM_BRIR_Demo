% [y, dev_max] = phase_manipulation(x, fs, type, Nfft_double, verbose)
%
% manipulates the phase of input signal x, to be zero-, linear- or
% minimum-phase
%
% INPUT (default):
% x           - input time signal [rows=samples x columns=channels]
% fs          - sampling frequency in Hz
% type        - desired phase behaviour 'lin', 'min', 'zero' ('min_phase etc. does work as well')
% Nfft_double - if type='min', artifacts from the hilbert transform produce
% (1)           errors that usually get smaller if the in put is zero
%               padded. Thus, the length of x is changed by zero padding
%               according to Nfft_double. 1=doubled, 2=quadroupled, etc.
%               (Only affects min-phase processing)
% verbose (1) - Differences between original and min phase magnitude
%               spectra is printed to screen. It is given separately for
%               frequencies that are 0-60dB, 60-100dB, and 100-inf dB below
%               the maximum. Set 1 or 0.
%
%
% OUTPUT:
% y           - impulse responses showing the desired phase behaviour
% dev_max     - maximum spectral deviations between original and min-phase
%               magnitude response in dB. Values are given for three ranges
%               of the spectrum. First value is the maximal deviation in
%               the upper 60 dB of the magnitude response, the second value
%               in the -60 to -100 dB range, and the third value covers
%               everything below -100 dB
%
% fabian.brinmann@tu-berlin.de, alexander.lindau@tu-berlin.de
% Audio Communication Group, TU Berlin, 09/2013
function [y, dev_max] = phase_manipulation(x, fs, type, Nfft_double, verbose)

id = find(type == '_');
if id
    type = type(1:id-1);
end
if ~exist('Nfft_double', 'var')
    Nfft_double = 1;
end
if ~exist('verbose', 'var')
    verbose = 1;
end
if ~isreal(x)
    error('phase_manipulation:input', 'Input signal x must be real')
elseif isnan(x)
    error('phase_manipulation:input', 'Input signal x must not contain any NaN values')
end

type = lower(type);
N    = size(x,1);
Nfft = nextpow2(N)+Nfft_double;

switch type
    case 'zero'
        % results in zerophase impulse response
        % [Moeser (2008):"Theroetische Akustik." Skript zur Vorlesung S.25]
        y = ifft(abs(fft(x)), 'symmetric');
    case 'lin'
        % desired group delay in seconds
        group_delay = (N-1)/2 / fs;
        % spacing between frequency bins
        delta_omega = 2*pi* fs/N;
        % gradient of phase response
        delta_phi   = -group_delay * delta_omega;
        
        % get magnitude spectrum
        Y_abs = abs(fft(x));
        % generate phase spectrum
        Y_ang = linspace(0, (N-1)*delta_phi, N)';
        Y_ang = repmat(Y_ang, 1, size(Y_abs,2));
        
        % construct complex spectrum
        Y = Y_abs .* exp(1j*Y_ang);
        
        % obtain impulse response (matlab takes care of symmetry)
        y = ifft(Y, 'symmetric');
    case 'min'
        % pre-processing for removal of hilbert-artifacts: windowing, zeropadding
        y = ifft(abs(fft(x)));
        y = circshift(y, floor(N/2));
        
        % symmetrical zeropadding to avoid artifact from hilbert transform
        y = [zeros(floor(2^Nfft-N/2),size(y,2)); y; zeros(ceil(2^Nfft-N/2),size(y,2))];
        
        % create Cepstral/Hilbert minmumphase
        for n = 1:size(x, 2)
            [~,y(:, n)] = rceps(y(:, n));
        end
        
        % cut out "Hilbert-echo" .. dabei geht Sperrdaempfung in Arsch!
        y = y(1:N, :);
        
        
        % get difference between minimun-phase and original-phase
        % magnitude responses
        X = 20*log10(abs(fft(x)));
        deviation = 20*log10(abs(fft(y))) - X;
        % set deviation = 0, for bins where the magnitude response is more
        % than x dB below the maximum
        dev_max = zeros(2, size(x,2));
        for m = 1:2
            dB_val = [60 100];
            for n = 1:size(x,2)
                tmp    = abs(deviation(:,n));
                dB_max = max(X(:,n));
                dB_id  = X(:,n) > dB_max-dB_val(m);
                dev_max(m,n) = max(tmp(dB_id));
            end
        end
        dev_max(m+1,:) = max(abs(deviation));
        dev_max        = max(dev_max,[],2);
        
        if verbose
            disp('Difference between original and minimum-phase magnitude response:')
            disp([ num2str(round(dev_max(1)*100)/100) ' dB (0-60 dB dynamic)'])
            disp([ num2str(round(dev_max(2)*100)/100) ' dB (60-100 dB dynamic)'])
            disp([ num2str(round(dev_max(3)*100)/100) ' dB (100-inf dB dynamic)'])
        end
        
        if nargout == 1
            clear dev_max
        end
        
        
    otherwise
        error(['phase_type ' phase_type ' not known. Use ''zero'', ''min'' or ''lin''']);
end

