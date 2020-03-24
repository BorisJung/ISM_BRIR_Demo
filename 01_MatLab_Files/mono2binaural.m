%% mono2binaural
%
% Creates a binaural signal from a mono signal and a BRIR.
%
% function BRIR_signal = mono2binaural(filename,fs,BRIR,write_audio,play_audio)
%
% Output: Binaural signal for left and right ear in a 2D-vector.
% 
% Input:
%           filename:       Filename of mono signal
%           fs:             Sampling Frequency
%           BRIR:           2D-Vector containing left and right channel of
%           BRIR
%           write_audio:    Additionally outputs a wav-file. (Boolean)
%           play_audio:     Plays binaural signal after calculations.
%           (Boolean)
%
%
function BRIR_signal = mono2binaural(filename,fs,BRIR,write_audio,play_audio)
%
%
if nargin<5
    disp('Not enough input arguments');
    return;
end;
%
%
%
% Read audio file
mono_sig = audioread(filename);
%
% check vector format
if size(BRIR,2) ~= 2
        BRIR = BRIR';
end;
%
%
% convolve signal
BRIR_signal = [conv(mono_sig,BRIR(:,1)) conv(mono_sig,BRIR(:,2))];
%
% normalization
BRIR_signal = 0.99 * (real(BRIR_signal) / max(abs(real(BRIR_signal(:)))));
%
% write audio file
if write_audio == true
    % get current folder
    yourFolder = pwd;
    % get name of current directory
    [parentFolder deepestFolder] = fileparts(yourFolder);
    % name new subfolder
    newSubFolder = sprintf('%s/created_WAVs', yourFolder);
    % create Subfolder if it doesn't exist
    if ~exist(newSubFolder, 'dir')
        mkdir(newSubFolder);
    end
    audiowrite('./created_WAVs/Binaural_Signal.wav',BRIR_signal,fs);
end;
%
%
% play audio file
if play_audio == true
    sound(BRIR_signal,fs);
end;

%
%
%% ======================= END OF SCRIPT ==================================
%
%
%