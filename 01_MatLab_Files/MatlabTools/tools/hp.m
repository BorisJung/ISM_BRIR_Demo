function [data h] = hp(data, plot_type, varargin)
% data = hp(data, plot_type, varargin)
% Can be used for various plots. Only 'data' and 'plot_type' are required
% arguments, everything else is passed as attribute value pairs
% e.g. hp(data, 'ir', 'hp_norm', 0). hp is intended for quick generation of
% plots that look ok. If you want your plot to look good, use the the
% handels h for manual adjustment. hp means hyper plot - fresh!
%
% Examples:
% hp(data, 's2d') .  .  .  .  .  .  .  .  .  . plots spectrum 2D
% hp(data, 's3d', 'hp_view', 'top_h')  .  .  . plots spectrum 3D
% hp(data, 's4', 'az', az, 'el', el)   .  .  . plots spherical spectrum
% hp(data, 'itd1', , 'az', az, 'el', el)  .  . plots spherical ITD
%
%
% INPUT
% (default values) {possible plot looks}
% (default Value (M) means Matlab default is used)
%
% Required ----------------------------------------------------------------
% data      : impulse responses [samples=rows x channels=columns]
%             if plot type is 'itd' or 'ild' data must be a struct
%             containing data.l and data.r with left and right ear hrir's
%
% plot_type : string specifiying the plot followed by '2D', '3D' or a
%             number specifying a spherical plot. (See argument sph_type)
% 
%   'ir'  : time signal   .  .  .  .  .  .  .  .  . {2D, 3D}
%   'etc' : energy time curve   .  .  .  .  .  .  . {2D, 3D}
%   'edc' : energy decay curve  .  .  .  .  .  .  . {2D, 3D}
%   'sr'  : step response .  .  .  .  .  .  .  .  . {2D, 3D}
%   's'   : spectrum   .  .  .  .  .  .  .  .  .  . {2D, 3D, spherical}
%   'ph'  : phase wrapped .  .  .  .  .  .  .  .  . {2D, 3D, spherical}
%   'pu'  : phase unwrapped  .  .  .  .  .  .  .  . {2D, 3D, spherical}
%   'g'   : group delay   .  .  .  .  .  .  .  .  . {2D, 3D}
%   'toa' : time of arrival  .  .  .  .  .  .  .  . {2D, spherical}
%   'itd' : interaural time difference   .  .  .  . {2D, spherical}
%   'ild' : interaural level difference  .  .  .  . {2D, spherical}
%   [sc]  : spectogram, single channel   .  .  .  . {3D}
%   [w]   : waterfall diagram , single channel .  . {3D}
%   'x'   : plot data as it is (no fft, no nothing) {spherical}
%           NOTE: in this cas data MUST only contain
%           one row but multiple cloumns.
%
% Universal arguments -----------------------------------------------------
% dr            : restrict data range
%                 - 2D: limits for y-axis [ymin ymax] (M)
%                 - 3D: limits for z-axis [zmin zmax]
%                       ([min(data) max(data)])
% du            : desired unit of data for plotting
%                 - 'dB' or 'lin' for magnitude plots (dB is defualt)
%                   in spherical plots du will refer to the phase and
%                   magnitude will allways be dB.
%                 - 'rad', 'deg' for phase plots and for spherical plots if
%                   color according to phase (deg)
%                 - TOA, ITD, g , 'n', 's', 'ms', 'us' (us)
% xu            : unit of x-axis
%                 - time-axis , 'n' samples, 's', 'ms', 'us'
%                   ('us' for ITD, 'ms' for other time plots)
%                 - frequency-axis, 'Hz', 'kHz' (kHz)
% norm_d (off)  : normalize data to norm_d, unit of according to du
% fs (44100)    : sampling frequency in Hz
% overlay (on)  : for plotting multiple data in one plot, 'on', 'off'
% labeling (on) : for default axis labels and figure titles, 'on', 'off'
%
% Plot layout -------------------------------------------------------------
% f_size (M)    : font size used for labeling plots. Axis font size is (M).
%
% Frequency domain plots --------------------------------------------------
% frac (off)    : fractional octave smoothing (frac=3 means 1/3 oct. sm.)
%
% Arguments for 2D plots --------------------------------------------------
% x             : axis limits for time or frequency axis [xmin xmax],
%                 values according to xu. ([min(data) max(data)], x differs
%                 for toa, itd, ild, see below)
% c ('k')       : - color string or rgb vector for coloring all channels
%                   the same (see doc colorspec)
%                 - 'cyc' cycles through colors as defined by Matlab
%                 - rgb matrix [channel number x 3] specifies color for
%                   each channel
%                 - rgb matrix [2 x 3] channel wise interpolation between
%                   first and second column
%                 - single value denotes greyscale (c = [c c c])
% line_w (.5)   : linewidth in points
% line_s (-)    : linestyle as string (see doc linespec)
% line_m (none) : line_m style as strin (see doc linespec)
%
% Arguments for 3D plots --------------------------------------------------
% x             : see 2D plots
% hp_view(top_v): Sets view angle for 3D plots. 'side', 'top_h', 'top_v',
%                 or [az el] (see doc view)
% y             : y vector specifying the data, e.g. angles y=-90:90. y
%                 should be monotonic increasing/decreasing
%                 (1:size(data,2))
% cm (jet)      : colormap to be used see doc colormap
% cb            : location of colorbar. 0 for not showing the colorbar.
%                 ('EastOutside', see doc colorbar)
% cr (128 steps): resolution of colormap. Unit according to du. zmin
%                 will be adjusted if cr does not perfectly fit the range
%                 [zmin zmax], because abs(z(2)-z(1))/cr must be an integer
% c_freeze (1)  : freezes the colorbar and colormap to allow different
%                 setttings and different subplots
%                 
%
% Arguments for spherical plots -------------------------------------------
% - x, y and z axis are marked by a red, green and blue line
% - positive x,y and z are markerd by litlle circles
% - x,y and z axis are limted by default
% Types of spherical plot. Use together with plot_type e.g. 's1' for ballon
% plot of spectrum.
%                 1: Balloon, fixed radius, variable color
%                 2: Balloon, variable radius and color
%                 3: Ballon, radius is magnitude, color is phase
%                 4: Balloon, variable radius fixed color
%                    (c can be used to pass a color, only [r g b] mode)
%                 5: Planar plot
%                 6: Polar plot (In this case only specifiy az, see below)
% az            : specification of azimuth. This has to be either a
%                 - vector with one entry per channel of data. In this case
%                   data is interpolated to a 1x1 degree resolution for
%                   surf plots
%                 - a matrix with same size as data. In this case sph_proc
%                   is set to 'none'
%                 If az contains values > 2pi, it is assumed that its unit
%                 is degree (if not rad are assumed)
% el            : specification of elevation. Same as 'az'. For polar 
%                 plots, only specify az.
% sph_f (1000)  : Frequency in Hz to be plotted
% sph_proc      : Data display for balloon plots
%                 - 'interp' : for balloon and planar plots, data is 
%                              interpoalted/extrapolted to a full sphere
%                              with 1 degree resolution.
%                              for polar plots, data is interpolated to
%                              full circle if range(az) > 180 degree
%                              (default)
%                 - 'tri'    : data is displayed after triangularization
%                 - 'none'   : data is displayed as it is, this is
%                              autodetected, see 'az'
% coord (M)     : specify the used coordinate convention. This also defines
%                 the orientation of x,y and z axis. New coordinate systems
%                 can be implemented in hp_coordinate_transformation.m.
%                 Implemente conventions are described in the same file.
%                 1: Matlab default (default, see doc sph2cart.)
%                 2: Mathematical
% axis_l (.5)   : line width of axis
% axis_m (10)   : marker size for marking posititve x,y and z in spherical
%                 plots
% cb, cm, cr,   : see 3D plots
% c_freeze
%
% Plot type specific arguments --------------------------------------------
% ---- toa, itd, ild ----
% x             : x vector specifying the data e.g. angles [-90:90]
%                 (1:size(data,2))
% toa_us (10)   : upsampling factor for onset detection, (only toa, itd)
% toa_th (6)    : threshold for onset detection in dB relative to
%                 channel-wise abs(max(data)), (only toa, itd)
%
%
% OUTPUT:
% data          : processed input data. This is a struct containing the
%                   processed data and
%                 - time or freq axis (2D plots)
%                 - time or freq and y atrix (3D plots)
%                 - azimuth and elevation matrix (balloon and planar plots)
%                 - azimuth axis (polar plots)
% h             : handle to line series or surface
%
%
% fabian.brinkmann@tu-berlin.de, Audio Communication Group TU Berlin,
% DFG research unit 'SEACEN', 7/2012
% 
% Thanks for hints: Alexander Lindau, Frank Schultz, Vera Erbes, SOFiA
% Toolbox (Bennjamin Bernschuetz), ITA Toolbox RWTH Aachen.


% PROCESSING
% 0:  Categorization of plot type into
%     - time, freq, other
%     - 2D, 3D, spherical
% 1:  Parsing of varargin
% 2:  Calculation of auxiliary variables
% 3:  Transformation of data (e.g. fft), definition of title string and
%     string for data format (e.g. 'phase in degree')
% 4:  Fractional octave smoothing
% 5:  Normalization
% 6:  Autorange
% 7:  Process data for spherical plots
% 8:  Define x-axis and x-label
% 8:  Plotting
% 10: Delete output argument if not requested
% NOTE: The code within the blocks sometimes is a bit redundant, but that
% makes it easier to understand and maintain.
%
%
% Overload variables    Difference attribute
% x                     2D, 3D, toa, itd, ild
% WARNING process overloaded varibales only within switch or if cases to
% avoid misbehaviour of function.
%
% Hardcoded variables
% min_dB = -300; to avoid clipping of -inf dB
% set(gcf, 'color', [1 1 1])
%
% toDo (ordered according to importance)
% - output format of az and el in data struct is allways matlab default
% - add 'plop' feature from ita surf
% - spectogram, waterfall diagram
% - possibility to pass input as frequency data (single- and two-sided)
%   I don't need this - do you?


%% ----------------------------------------- 0. categorization of plot type
if ~isnan(str2double(plot_type(end)))
    plot_look = str2double(plot_type(end));
    plot_type = plot_type(1:end-1);
elseif length(plot_type) >= 3
    if strcmpi(plot_type(end-1:end), '2D')
        plot_look = '2D';
        plot_type = plot_type(1:end-2);
    elseif strcmpi(plot_type(end-1:end), '3D')
        plot_look = '3D';
        plot_type = plot_type(1:end-2);
    end
else
    error([plot_type ' is not a valid plot_type. See doc hp for a list of valid plot_types'])
end



if max(strcmpi(plot_type, {'ir' 'etc' 'edc' 'sr'}))
    hp_domain = 'time';
elseif max(strcmpi(plot_type, {'s' 'g' 'ph' 'pu', 'sc' 'w'}))
    hp_domain = 'freq';
elseif max(strcmpi(plot_type, {'toa' 'itd' 'ild' 'x'}))
    hp_domain = plot_type;
else
    error([plot_type ' is not a valid plot_type. See doc hp for a list of valid plot_types'])
end

%% ------------------------------------------ 1a. define default parameters
inputnames = {...
'fs', 'overlay', 'norm_d', 'labeling', ...  % universal parameters
'f_size', 'axis_m', 'axis_l', ...           % plot layout
'xu', 'frac', 'du', ...                     % frequency and time parameters
'hp_view', 'x', 'y', 'dr',...               % 2D and 3D plots
'c', 'line_w', 'line_s', 'line_m',...
'cb', 'cm', 'cr', 'c_freeze',...
'az', 'el', 'sph_f', 'sph_proc', 'coord',...% spherical plots parameters
'toa_us', 'toa_th'...                       % for plots toa, itd, g
};

% Define default values
def.fs       = 44100;
def.overlay  = 'on';
def.norm_d   = false;
def.labeling = 'on';

def.f_size  = get(0,'DefaultAxesFontSize');
def.axis_l  = .5;
def.axis_m  = 10;

if strcmpi(plot_type, 'itd')
    def.xu  = 'us';
elseif strcmpi(hp_domain, 'time')
    def.xu  = 'ms';
else
    def.xu  = 'kHz';
end
def.frac    = false;

if strcmpi(plot_type, 'itd')
    def.du  = 'us';
elseif max(strcmpi(plot_type, {'ph' 'pu'})) || ~ischar(plot_look)
    def.du  = 'deg';
elseif strcmpi(plot_type, 's')
    def.du  = 'dB';
else
    def.du  = 'us';
end

if strcmpi(plot_look, '3D')
    def.hp_view = 'top_v';
elseif ~ischar(plot_look) && plot_look <= 4
    def.hp_view = [135 15];
else
    def.hp_view = false;
end
def.x       = false;
def.y       = false;
def.dr      = false;

def.c       = 'k';
def.line_w  = .5;
def.line_s  = '-';
def.line_m  = 'none';

def.cb       = 'EastOutside';
def.cm       = 'jet';
def.cr       = false;
def.c_freeze = 1;

def.az      = false;
def.el      = false;
def.sph_f   = 1000;
if plot_look ~= 6;
    def.sph_proc= 'interp';
else
    def.sph_proc= 'interp';%'none';
end
def.coord   = 1;

def.toa_us  = 10;
def.toa_th  = 6;



% --------------------------------------------------------- 1b. parse input
% Check the number of arguments to be even
if mod(size(varargin,2),2)
    error('No even number of arguments. Check submitted attribute/value pairs.')
end

% Stock the values
for n=1:2:size(varargin,2)-1
    % checks for all possible attribute/value pairs
    is_parameter = 0;
    for m=1:size(inputnames,2) 
        % detect submitted attribute/value pairs
        if strcmp(inputnames{m},varargin{n})
            % create input-variables from submitted attribute-value-pairs (supress output)
            [~] = evalc([inputnames{m},'=','varargin{n+1}']);
            is_parameter = 1;
        end
    end
    if ~is_parameter
        error(['No such parameter: ' varargin{n}])
    end
end

% Create missing input-variables with default-values
for m=1:size(inputnames,2)
    default = ['def.' inputnames{m}];
    % if input-variable m hasn't already been defined...
    if ~exist(inputnames{m},'var')
        [~] = evalc([inputnames{m},'=',default]); % (supress output)
    end
end

% delete unnecassary variables to have a clean workspace inside the
% function (better for debugging)
clearvars def default in_args inputnames m n ninput trash varargin


%% ------------------------------------------------ 2.  auxiliary variables
% copy original input
data_cp = data;

% number of samples and channels
if ~isstruct(data)
    [N,C] = size(data);
else
    [N,C] = size(data.l);
end

% frequency vector
f     = (0:N-1)' * fs/N;

% minimum to be displayed in 3D plots if data has inf entries
min_dB = -300;

% check data format of az for spherical plots
if ~islogical(az)
    if size(az,1) == size(data,1) && size(az,2) == size(data,2)
        sph_proc = 'none';
        % disable processing of the data as done in next step
        plot_type_cp = plot_type;
        plot_type = 'spherical no transformation';
    else
        % transpose azimuth and elevation if necessary
        if size(az,1) < size(az,2)
            az = az';
        end
        if exist('el', 'var') && size(el,1) < size(el,2)
            el = el';
        end
    end
end


%% ----------------------------------- 3. data transformation and data type
% -------------------------------------------------------- ir (time signal)
if strcmpi(plot_type, 'ir')
    d_type = 'Amplitude';
    t_str  = 'Impulse response';
end
% ------------------------------------------------- etc (energy time curve)
if strcmpi(plot_type, 'etc')
    data   = 20*log10(abs(data));
    d_type = 'Amplitude in dB';
    t_str  = 'Energy Time Curve (ETC)';
end
% ------------------------------------------------ edc (energy decay curve)
if strcmpi(plot_type, 'edc')
    data   = edc(data);
    d_type = 'EDC in dB';
    t_str  = 'Energy Decay Curve (EDC)';
end
% ------------------------------------------------------ sr (step response)
if strcmpi(plot_type, 'sr')
    for k = 1:C
        data(:,k) = stepz(data(:,k), 1, N, fs);
    end
    clear k
    d_type = 'Amplitude';
    t_str  = 'Step response';
end
% -------------------------------------------------- s (magnitude spectrum)
if strcmpi(plot_type, 's')
    if strcmpi(du, 'lin')
        data = abs(fft(data));
        d_type = 'Amplitude';
        t_str  = 'Magnitude response';
    else
        data = 20*log10(abs(fft(data)));
        d_type = 'Amplitude in dB';
        t_str  = 'Magnitude response';   
    end
end
% --------------------------------------------------------- g (group delay)
if strcmpi(plot_type, 'g')
    for k = 1:C
        data(:,k)   = grpdelay(data(:,k), 1, f, fs);
    end
    clear k
    t_str  = 'Group delay';
    
    if strcmpi(du, 'n')
        d_type = 'group delay in samples';
    elseif strcmpi(du, 's')
        data   = data/fs;
        d_type = 'group delay in ms';
    elseif strcmpi(du, 'ms')
        data   = data/fs*10^3;
        d_type = 'group delay in ms';
    elseif strcmpi(du, 'us')
        data   = data/fs*10^6;
        d_type = 'group delay in us';
    else
        error([du ' is not a valid data unit. See doc hp for a list of valid d_units'])
    end
end
% ------------------------------------------- ph (phase response - wrapped)
if strcmpi(plot_type, 'ph')
    data   = angle(fft(data));
    t_str  = 'Phase response';
    
    if strcmpi(du, 'rad')
        d_type = 'phase in radians';
    elseif strcmpi(du, 'deg')
        data   = rad2deg(data);
        d_type = 'phase in degree';
    else
        error([du ' is not a valid data unit. See doc hp for a list of valid d_units'])
    end
end
% ----------------------------------------- ph (phase response - unwrapped)
if strcmpi(plot_type, 'pu')
    data   = unwrap(angle(fft(data)));
    t_str  = 'Phase response';
    
    if strcmpi(du, 'rad')
        d_type = 'phase in radians';
    elseif strcmpi(du, 'deg')
        data   = rad2deg(data);
        d_type = 'phase in degree';
    else
        error([du ' is not a valid data unit. See doc hp for a list of valid d_units'])
    end
end
% --------------------------------------------------- toa (time of arrival)
if strcmpi(plot_type, 'toa')
    data   = onset_detect(data, toa_us, -(abs(toa_th)));
    d_type = 'TOA in samples';
    t_str  = 'Time of Arrival (TOA)';
end
% ---------------------------------------- itd (interaural time difference)
if strcmpi(plot_type, 'itd')
    tmp.l  = onset_detect(data.l, toa_us, -(abs(toa_th)));
    tmp.r  = onset_detect(data.r, toa_us, -(abs(toa_th)));
    data   = tmp.l - tmp.r;
    d_type = 'ITD in samples';
    t_str  = 'Interaural Time Difference (ITD)';
    clear tmp
end

if max(strcmpi(plot_type, {'toa' 'itd'}))
    if strcmpi(du, 'n')
        d_type = [upper(plot_type) 'in samples'];
    elseif strcmpi(du, 's')
        data   = data/fs;
        d_type = [upper(plot_type) 'in s'];
    elseif strcmpi(du, 'ms')
        data   = data/fs*10^3;
        d_type = [upper(plot_type) 'in ms'];
    elseif strcmpi(du, 'us')
        data   = data/fs*10^6;
        d_type = [upper(plot_type) 'in us'];
    end
end
% --------------------------------------- ild (interaural level difference)
if strcmpi(plot_type, 'ild')
    % same calculation as in sound field synthesis toolbox, but sign is
    % switched
    data   = db(rms(data.l)) - db(rms(data.r));
    d_type = 'ILD in dB';
    t_str  = 'Interaural Level Difference (ILD)';
    clear tmp
end


if strcmpi(plot_type, 'spherical no transformation')
    plot_type = plot_type_cp;
    clear plot_type_cp
end
%% ----------------------------------------- 4. fractional octave smoothing
% (discard frequency bin at fs/2 if N is even)
if ~islogical(frac)
    data = fract_oct_smooth(data(1:ceil(N/2),:), 'welti', fs, frac);
    f    = f(1:ceil(N/2));
end


%% ------------------------------------------------------- 5. normalization
if ~islogical(norm_d)
    if strfind(upper(d_type), 'DB')
        data = data -max(max(data)) + norm_d;
    else
        data = data / max(max(abs(data))) * norm_d;
    end
end

%% ----------------------------------------------------------- 6. autorange
if islogical(dr)
    % find maximum and set visible plot range to next highest values that
    % is a divisor of 10
    tmp_max = max(max(max(data)));
    tmp_max = tmp_max - mod(tmp_max, 10)+10;
    if ceil(max(max(max(data)))) == tmp_max
        tmp_max = tmp_max + 10;
    end
    
    % set minimum value according to plot_type
    if strcmpi(plot_type, 'etc')
        dr = [tmp_max-100 tmp_max];
    elseif strcmpi(plot_type, 's')
        dr = [tmp_max-60 tmp_max];
    end
end

clear tmp_max

%% ------------------------------------ 7. process data for spherical plots
if ~ischar(plot_look) && max(strcmpi(plot_type, {'s' 'ph' 'pu' 'toa' 'itd' 'ild' 'x'}))
    if (strcmpi(sph_proc, 'none') && plot_look ~= 6) || ...
        strcmpi(plot_type, 'x')     
        data_c = data;
    else
    
        % get magnitude and phase spectrum or toa, itd, ild and assign to
        % new variables for further processing
        if strcmpi(plot_type, 's')
            data_abs = data;
            data_ang = angle(fft(data_cp));
        elseif max(strcmpi(plot_type, {'ph' 'pu'}))
            data_abs = 20*log10(abs(fft(data_cp)));
            data_ang = data;
        else max(strcmpi(plot_type, {'toa', 'itd', 'ild'}));
            if plot_look ~= 1
                error(['plot_type ' plot_type ' only supports ' plot_type '1'])
            end
            data_abs = data;
            data_ang = data;
        end
        
        % get data to be plotted (In ballon and polar plots this determines
        % the in planar plots it determines the color)
        if plot_look == 1
            data = ones(size(data));
        elseif strcmpi(plot_type, 's') || plot_look == 3 || ...
               max(strcmpi(plot_type, {'toa', 'itd', 'ild'}))
            data   = data_abs;
            d_type = 'amplitude in dB';
        else
            data = data_ang;
            if strcmpi(du, 'rad')
                d_type = 'phase in radians';
                if strcmpi(plot_type, 'pu')
                    dr = [min(min(data)) max(max(data))];
                else
                    dr = [-pi pi];
                end
            elseif strcmpi(du, 'deg')
                if strcmpi(plot_type, 's')
                    data   = rad2deg(data);
                end
                d_type = 'phase in degree';
                if strcmpi(plot_type, 'pu')
                    dr = [min(min(data)) max(max(data))];
                else
                    dr = [-180 180];
                end
            else
                error([du ' is not a valid data unit. See doc hp for a list of valid d_units'])
            end
        end
        
        % get data that determines the color to be plotted.
        % auto set data range for color, if the phase determines the color.
        if (max(plot_look == [1 2 4 5 6]) && strcmpi(plot_type, 's')) || ...
            max(strcmpi(plot_type, {'toa', 'itd', 'ild'}))
            data_c = data_abs;
        else
            data_c = data_ang;
            if strcmpi(du, 'rad')
                d_type = 'phase in radians';
                if strcmpi(plot_type, 'pu')
                    dr = [min(min(data_c)) max(max(data_c))];
                else
                    dr = [-pi pi];
                end
            elseif strcmpi(du, 'deg')
                if strcmpi(plot_type, 's')
                    data_c   = rad2deg(data_c);
                end
                d_type = 'phase in degree';
                if strcmpi(plot_type, 'pu')
                    dr = [min(min(data_c)) max(max(data_c))];
                else
                    dr = [-180 180];
                end
            else
                error([du ' is not a valid data unit. See doc hp for a list of valid d_units'])
            end 
        end
        
        % pick the frequency to be plotted
        if max(strcmpi(plot_type, {'s' 'ph' 'pu'}))
            tmp    = round(sph_f / (fs/N)) + 1;
            data   = data(tmp,:);
            data_c = data_c(tmp,:);
            sph_f  = f(tmp);
        else
            sph_f = upper(plot_type);
        end
        
        clear data_abs data_ang tmp
    end
    
    % coordinate conversion from user to matlab default format
    co = hp_coordinate_transformation(az, el, coord);
    
    % clip data to dr if its information is used as radius in balloon or
    % polar plots
    if ~islogical(dr) && max(plot_look == [2:4 6])
        data = max(data, dr(1));
        data = min(data, dr(2));
        if exist('data_c', 'var')
            data_c = max(data_c, dr(1));
            data_c = min(data_c, dr(2));
        end
    % clip data that is used for coloration to min_dB. Otherwise the
    % colorbar can not be constructed, if data_c contains -inf values.
    else
        data(isinf(data)) = min_dB;
        if exist('data_c', 'var')
            data_c(isinf(data_c)) = min_dB;
        end
    end
    % if data is as radius in balloon plots, negative values don't make any
    % sense. Everthing has to be shiftet to positive values
    if min(min(data)) < 0 && max(plot_look == 2:4)
        data = data - min(min(data)) + 1;
    end

    % processing for polar plots
    if plot_look == 6
        % data has to be sorted
        [co.az, id] = sort(co.az);
        data = data(id);
        % data has to be clipped for polar plots
        if ~islogical(dr)
            data = max(data, dr(1));
            data = min(data, dr(2));
        else
            dr = [min(data) max(data)];
        end
    end
    
    % interpolation and triangularization
    if plot_look == 6
    % polar plots
        if strcmpi(sph_proc, 'interp')
            if range(co.az) > pi
            % interpolate to full circle
                % copy smalles value to 0 and 360 degree to avoid
                % extrapolation which would result in a discontinuity at 0
                [tmp_min, tmp_id] = min(co.az);
                if tmp_min == 0 && max(co.az) < 2*pi
                    co.az = [co.az; 2*pi];
                    data  = [data data(tmp_id)];
                else
                    co.az = [0; co.az; 2*pi];
                    data  = [data(tmp_id) data data(tmp_id)];
                end
                % interpolate
                data = interp1(co.az, data, (0:360)/180*pi);
                co.az = (0:360)'/180*pi;
                
                clear tmp_min tmp_id
            else
            % interpolate to half circle
                % construct new azimuth
                tmp = min(co.az):pi/180:max(co.az);
                % check if the last value of co.az is included in the new
                % azimuth vector
                if co.az(end) ~= tmp(end)
                    tmp = [tmp co.az(end)];
                end
                % interpolate
                data = interp1(co.az, data, tmp);
                co.az = tmp;
            end   
        elseif strcmpi(sph_proc, 'tri') 
            error('Triangularization is not implemeted for polar plots')
        else
            error([sph_proc ' is not a valid shperical processing. See doc hp.'])
        end
    else
    % planar and balloon plots
        if strcmpi(sph_proc, 'interp')
            % this warning does not affect us
            warning('off', 'MATLAB:TriScatteredInterp:DupPtsAvValuesWarnId')
            
            % get radius data
            if plot_look == 1
                % ballon with radius 1
                [az_tmp, el_tmp] = meshgrid((0:360)*pi/180, (90:-1:-90)*pi/180);
                [X, Y, Z] = sph2cart(az_tmp, el_tmp, ones(181,361));
                clear az_tmp el_tmp
            else
                % radius of ballon
                [X Y Z data] = hp_data_interpolation(co, data);
            end
            
            % get color data
            [~,~,~, data_c] = hp_data_interpolation(co, data_c);
            
            warning('on', 'MATLAB:TriScatteredInterp:DupPtsAvValuesWarnId')
            
        elseif strcmpi(sph_proc, 'tri')
            % Balloon plots
            if plot_look ~= 5
                % create trangularization with unity radius
                [X Y Z] = sph2cart(co.az, co.el, ones(size(co.az)));
                tri = DelaunayTri([X Y Z]);
                ch = convexHull(tri);
                % obtain x, y, z coordinates of vertices with correct
                % radius for plotting
                [X Y Z] = sph2cart(co.az, co.el, data');
                
            % Planar plots
            else
                 % create 2D trangularization
                 tri = DelaunayTri(az, el);
            end
        elseif strcmpi(sph_proc, 'none')
            if plot_look == 1
                [X, Y, Z] = sph2cart(co.az, co.el, ones(size(data)));
            else
                [X, Y, Z] = sph2cart(co.az, co.el, data);
            end
        else
            error([sph_proc ' is not a valid argument. See doc hp'])
        end
    end    
    
elseif ~ischar(plot_look) && ~max(strcmpi(plot_type, {'s' 'ph' 'pu' 'toa' 'itd' 'ild', 'x'}))
    error(['spherical plots are not possible for plot_type ' plot_type])
end


%% ------------------------------------------------- 8. x-axis and x-labels
% (from samples to sec, etc.)
if strcmpi(xu, 's')
    n = (0:N-1) * 1/fs;
    x_label = 't in s';
elseif strcmpi(xu, 'ms')
    n =  (0:N-1) * 1/fs * 10^3;
    x_label = 't in ms';
elseif strcmpi(xu, 'us')
    n = (0:N-1) * 1/fs * 10^6;
    x_label = 't in us';
elseif strcmpi(xu, 'n')
    n = 1:N;
    x_label = 't in samples';
elseif strcmpi(xu, 'hz')
    n = f;
    x_label = 'f in Hz';
    % default xtick and label for frequency plots
    f_tick  = [1:10 20:10:100 200:100:1000 2000:1000:10000 20000:10000:40000]';
    f_label = {'1'     '' '' '' '' '' '' '' '' ...
               '10'    '' '' '' '' '' '' '' '' ...
               '100'   '' '' '' '' '' '' '' '' ...
               '1000'  '' '' '' '' '' '' '' '' ...
               '10000' '20000' '' '40000'}';
elseif strcmpi(xu, 'khz')
    n = f;
    x_label = 'f in kHz';
    % frequency axis   
    % default xtick and label for frequency plots
    f_tick  = [1:10 20:10:100 200:100:1000 2000:1000:10000 20000:10000:40000]';
    f_label = {'.001' '' '' '' '' '' '' '' '' ...
               '.01'  '' '' '' '' '' '' '' '' ...
               '.1'   '' '' '' '' '' '' '' '' ...
               '1'    '' '' '' '' '' '' '' '' ...
               '10'   '20' '' '40'}';
else
    error([xu ' is not a valid x axis unit. See doc hp for a list of valid x_units'])
end

% define default range and xlabels for toa, itd and ild
if max(strcmpi(hp_domain, {'toa' 'itd' 'ild'}))
    % range
    if islogical(x)
        n = 1:C;
    else
        n = x;
        x = [n(1) n(end)];
    end
    % xlabel is somehow a guess
    if strcmpi(hp_domain, {'toa'})
        x_label = 'data points';
    else
        x_label = 'source/head orientation';
    end
end



%% ----------------------------------------- 9. plot preparations and plots

warning('off', 'MATLAB:warn_r14_stucture_assignment');

% check x limits
if islogical(x)
    x = [min(n) max(n)];
    if strcmpi(hp_domain, 'freq')
        x(1) = n(2);
        x(2) = fs/2;
    end
end

% parse color and interpolate if necessary
if ~ischar(c) && size(c,1) == 2
    tmp = c;
    c      = linspace(tmp(1,1), tmp(2,1), C)';
    c(:,2) = linspace(tmp(1,2), tmp(2,2), C)';
    c(:,3) = linspace(tmp(1,3), tmp(2,3), C)';
elseif ~ischar(c) && numel(c) == 1
    c = [c c c];
end

% -------- 2D plots -------- %
if strcmpi(plot_look, '2D')
    % plot
    eval(['hold ' overlay])
    h = plot(n, data, 'LineWidth', line_w, 'LineStyle', line_s, 'Marker', line_m);
    hold off
    
    % color
    if ~strcmpi(c, 'cyc')
        if size(c,1) ~= C
            set(h, 'color', c)
        else
            for k = 1:C
                set(h(k), 'color', c(k,:))
            end
            clear k
        end
    end
    
    % logarithmic x axis for frequency plots
    if strcmpi(hp_domain, 'freq')
        set(gca, 'xscale', 'log')
    end
    
    box on
    
    % limits, grid and tick
    if length(n) > 1
        xlim([min(x) max(x)])
    end
    if ~islogical(dr)
        ylim([dr(1) dr(2)])
    end
    if strcmpi(hp_domain, 'freq')
        set(gca, 'xTick', f_tick(f_tick >= min(x) & f_tick <= max(x)),...
           'xtickLabel', f_label(f_tick >= min(x) & f_tick <= max(x)))
    end
    grid on
    
    % labeling
    if strcmpi(labeling, 'on')
        xlabel(x_label, 'fontsize', f_size)
        ylabel(d_type, 'fontsize', f_size)
        title(t_str, 'fontsize', f_size)
    end
    
    % set output struct
    data.data = data;
    data.n    = n;
    data.h    = h;
 
    
% -------- 3D plots -------- %
elseif strcmpi(plot_look, '3D')
    
    % default data range
    if islogical(y)
        y = 1:C;
    end
    
    % generate mesh
    [Y, N] = meshgrid(y, n);
    
    % clip data for plotting (do not clip the data that is returned by the function)
    data_p = data;
    if ~islogical(dr)
        data_p = max(dr(1), data_p);
        data_p = min(dr(2), data_p);
    end
    
    
    % surf with some hardcoded EdgeColor
    eval(['hold ' overlay])
    h = surf(N, Y, data_p, 'EdgeColor', 'none');
    hold off
    
    % set view
    if ischar(hp_view)
        if strcmpi(hp_view, 'top_h')
            if y(1)<y(2)
                view([0 -90])
            else
                view([0 90])
            end
        elseif strcmpi(hp_view, 'top_v')
            if y(1)<y(2)
                view([90 -90])
            else
                view([-90 90])
            end
        elseif strcmpi(hp_view, 'side')
            view([0 0])
        else
            error([hp_view 'is not a valid view type. See doc hp for valid view types'])
        end
    else
        view([hp_view(1) hp_view(2)])
    end
    
    % logarithmic x axis for frequency plots
    if strcmpi(hp_domain, 'freq')
        set(gca, 'xscale', 'log')
    end
    
    % limits, grid and tick
    grid on
    box on
    if strcmpi(hp_domain, 'freq')
        set(gca, 'xTick', f_tick(f_tick >= min(x) & f_tick <= max(x)),...
           'xtickLabel', f_label(f_tick >= min(x) & f_tick <= max(x)))
    end
    xlim([x(1) x(2)])
    ylim([min(y) max(y)])
    % clipping visualization of z axis (using 'dr') and defining colormap
    if islogical(dr)
        % if not specified, find min and max within given x-range (time, or
        % freqeuncy axis) and use that for clipping
        [~,a] = min(abs(x(1)-n));
        [~,b] = min(abs(x(2)-n));
        dr = [min(min(data(a:b,:))) max(max(data(a:b,:)))];
        % avoid zmin = inf
        dr(isinf(dr)) = min_dB;
        clear a b
    end
    
    % colorbar and colormap
    dr = hp_set_colormap_and_colorbar(cm, cb, cr, dr);
    set(gca, 'CLim', [dr(1) dr(2)], 'zlim', [dr(1) dr(2)]);
    
    if c_freeze
        if cb ~= 0
            cbfreeze
        end
        freezeColors
    end
    
    % labeling
    if strcmpi(labeling, 'on')
        xlabel(x_label, 'fontsize', f_size)
        ylabel('data points', 'fontsize', f_size)
        title(t_str, 'fontsize', f_size)
    end
    
    % set output struct
    data.data = data;
    data.N    = N;
    data.Y    = Y;
    data.h    = h;

    
% -------- spherical plots -------- %
elseif ~ischar(plot_look)
    
    % -------- balloon plots -------- %
    if plot_look <= 4
        
        % surf with some hardcoded EdgeColor
        eval(['hold ' overlay])
        
        % plot with or without color information
        if plot_look ~= 4
            if strcmpi(sph_proc, 'tri')
                h = trisurf(ch, X, Y, Z, data_c, 'EdgeColor', 'none');
            else
                h = surf(X, Y, Z, data_c, 'EdgeColor', 'none');
            end
        else
            if strcmpi(sph_proc, 'tri')
                if ~ischar(c)
                    h = trisurf(ch, X, Y, Z, data_c, 'EdgeColor', 'none', 'FaceColor', c);
                else
                    h = trisurf(ch, X, Y, Z, data_c, 'EdgeColor', 'none', 'FaceColor', [1 121/255 4/255]);
                end
            else
                if ~ischar(c)
                    h = surf(X, Y, Z, 'EdgeColor', 'none', 'FaceColor', c);
                else
                    h = surf(X, Y, Z, 'EdgeColor', 'none', 'FaceColor', [1 121/255 4/255]);
                end
            end
        end
        
        axis equal
        axis off
        box off
        hold off
        rotate3d on
                
        % colorbar and colormap or fixed color
        if plot_look ~= 4
            if islogical(dr)
                clear dr
                dr(1) = min(min(data_c));
                dr(2) = max(max(data_c));
            end
            dr = hp_set_colormap_and_colorbar(cm, cb, cr, dr);
            set(gca, 'CLim', [dr(1) dr(2)]);
        else
            light;
            lighting phong;
        end
        
        if c_freeze
            if cb ~= 0
                cbfreeze
            end
            freezeColors
        end
        
        % draw coordinate system and limit x,y,z axis
        co.L = 1.2 * max(max(sqrt(X.^2 + Y.^2 + Z.^2)));
        hp_draw_coordinates(co, 1, axis_m, axis_l)
        axis([-co.L co.L -co.L co.L -co.L co.L])
        
        % set view and use default if not specified
        if ischar(hp_view)
            if strcmpi(hp_view, 'top_h') || strcmpi(hp_view, 'top_v')
                view([0 90])
            elseif strcmpi(hp_view, 'side')
                view([0 0])
            else
                error([hp_view 'is not a valid view type. See doc hp for valid view types'])
            end
        elseif ~islogical(hp_view)
            view([hp_view(1) hp_view(2)])
        end
        
        % plot title (axis are not labeled because they are off)
        if sph_f < 10000 & ~max(strcmpi(plot_type, {'toa' 'itd' 'ild'}))
            title_str = ['f = ' num2str(round(sph_f)) ' Hz '];
        elseif sph_f >= 10000 & ~max(strcmpi(plot_type, {'toa' 'itd' 'ild'}))
            title_str = ['f = ' num2str(round(sph_f/10)*10/1000) ' kHz '];
        else
            title_str = '';
        end
        
        if strcmpi(plot_type, 's')
            tmp = 'abs';
        elseif strcmpi(plot_type, 'ph')
            tmp = 'phase';
        elseif strcmpi(plot_type, 'pu')
            tmp = 'phase';
        else
            tmp = upper(plot_type);
        end
        
        switch plot_look
            case 1; title([title_str '(radius=fix, color= ' tmp ')']);
            case 2; title([title_str '(radius=' tmp ', color=' tmp ')']);
            case 3; title([title_str '(radius=abs, color=phase)']);
            case 4; title([title_str '(radius=' tmp ', color=fix)']);
        end
        clear title_str tmp
        
        [az_tmp, el_tmp] = meshgrid((0:360)*pi/180, (90:-1:-90)*pi/180);
        
        % set output struct
        if strcmpi(sph_proc, 'interp')
            data.data = data_c;
            data.az   = az_tmp;
            data.el   = el_tmp;
            data.h    = h;
        else
            data.data = data_c;
            data.az   = co.az;
            data.el   = co.el;
            data.h    = h;
        end
        
        
    % -------- planar plots -------- %
    elseif plot_look == 5
        
        % plot data
        if strcmpi(sph_proc, 'tri')
            % get azimuth and elevation of triangles defined by tri
            az_tmp = az(tri(:,:))';
            el_tmp = el(tri(:,:))';
            
            % get corresponding color data (the mean of all three points
            % that constitute a triangle is used)
            c = zeros(3, size(az_tmp,2));
            for n = 1:size(az_tmp,2)
                for m = 1:3
                    c(m, n) = data(az == az_tmp(m,n) & el == el_tmp(m,n));
                end
            end
            
            % plot using patch
            patch(az_tmp, el_tmp, mean(c), 'EdgeColor', 'none')
            
            view([0 90])
            xlim([min(az) max(az)])
            ylim([min(el) max(el)])
        else
            imagesc(data)
        end
        
        % colormap and colorbar
        if islogical(dr)
            clear dr
            dr(1) = min(min(data));
            dr(2) = max(max(data));
        end
        dr = hp_set_colormap_and_colorbar(cm, cb, cr, dr);
        set(gca, 'CLim', [dr(1) dr(2)]);
        
        if c_freeze
            if cb ~= 0
                cbfreeze
            end
            freezeColors
        end
        
        % x and y axis labeling and ticks
        if strcmpi(sph_proc, 'interp')
            set(gca, 'xTick', co.xtick, 'xTickLabel', co.xticklabel,...
                     'yTick', co.ytick, 'yTickLabel', co.yticklabel)
        end
        
        if strcmpi(labeling, 'on')
            xlabel 'Azimuth'
            ylabel 'Elevation'
        end
        box on
        
        % set output struct
        data.data = data;
        data.az   = co.az;
        data.el   = co.el;
        
        
        
    % -------- polar plots -------- %
    elseif plot_look == 6
        
        eval(['hold ' overlay])
        if islogical(dr)
            h = mmpolar(co.az, data, 'Style', 'compass');
        else
            h = mmpolar(co.az, data, 'Style', 'compass', 'Rlimit', dr);
        end
        hold off
        
        % set line color and width
        set(h, 'color', c, 'lineWidth', line_w, 'LineStyle', line_s, 'marker', line_m)
    
        % set output struct
        data.data = data;
        data.az   = co.az;
        data.h    = h;
        
    % -------- wrong type -------- %
    else
        error([plot_look 'is not a valid spherical plot. See doc hp for valid plots'])
    end
end

set(gcf, 'color', [1 1 1])
warning('on', 'MATLAB:warn_r14_stucture_assignment');


%% ----------------------------- 10. delete output argument if not requested
if nargout == 0
    % no output on the screen if you foreget ';' in 'hp(data, 's');'
    clear data
end
