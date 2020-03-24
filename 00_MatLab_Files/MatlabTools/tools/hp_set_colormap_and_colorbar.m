function dr = hp_set_colormap_and_colorbar(cm, cb, cr, dr)

% set colormap
if islogical(cr)
    % fixed resolution of 128 steps
    eval(['colormap(' cm '(128))'])
    
else
    % resolution according to cr
    num_steps = dr(2):-(abs(cr)):dr(1);
    % adjust zmin, to make sure that cr fits into the data range
    if num_steps(end) > dr(1)
        dr(1) = num_steps(end) - abs(cr);
    end
    clear num_steps
    % still rounding is needed to ensure that an integer number is
    % passed to colormap (numerical inaccuracies if zmin and zmax
    % are not integer values)
    eval(['colormap(' cm '(round(abs(dr(2)-dr(1))/cr)))'])
end
if cb ~= 0
    colorbar('location', cb)
end