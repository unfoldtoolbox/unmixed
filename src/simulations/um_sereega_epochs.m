function EEG  =um_sereega_epochs(varargin)

cfg = finputcheck(varargin,...
    {'debug',   'boolean', [], 0;...
     'epochlength','real',[],0.5;... %s
     'noise_components','real',[],10;...
     'srate','real',[],100;... % Hz
     'n_epochs','real',[],100;... % two conditions, per condition
     'noise_orient','boolean',[],0; % add random orientation noise to sources
    },'mode','ignore');


%%
lf = lf_generate_fromnyhead('montage', 'S64');

if cfg.debug
    plot_headmodel(lf);
end
%% noise

noise = struct( ...
    'type', 'noise', ...
    'color', 'brown', ...
    'amplitude', 1);
noise = utl_check_class(noise);
% 
% noise_white = struct( ...
%     'type', 'noise', ...
%     'color', 'white', ...
%     'amplitude', cfg.noiseAmplitude);
% noise_white = utl_check_class(noise_white);


%%
n1_comp = [];

n1_comp.signal{1} = struct();
n1_comp.signal{1}.peakLatency = 170;      % in ms, starting at the start of the epoch
n1_comp.signal{1}.peakWidth = 100;        % in ms
n1_comp.signal{1}.peakAmplitude = -5;      % in microvolt
n1_comp.signal{1}.peakLatencyDv = 50;
n1_comp.signal{1}.peakAmplitudeDv = 3;

n1_comp.signal{1} = utl_check_class(n1_comp.signal{1}, 'type', 'erp');

if cfg.debug
    plot_signal_fromclass(n1_comp.signal{1}, epochs);
end
% Add the noise
% vis.signal{2} = noise;


% right visual cortex
n1_comp.source= lf_get_source_nearest(lf, [50 -40 -25]);
if cfg.debug
    plot_source_location(n1_comp.source, lf, 'mode', '3d');
end

% n1_comp.orientation = utl_get_orientation_pseudotangential(n1_comp.source,lf);
n1_comp.orientationDv = [0 0 0];

n1_comp.orientation = [0.5,-0.46,1];

if cfg.debug
    %%
    for k = 1:10
        n1_comp.orientation = utl_get_orientation_random(1);
        plot_source_projection(n1_comp.source, lf, 'orientation', n1_comp.orientation,'orientedonly',1);
        %     title(n1_comp.orientation)
    end
end

p1_comp = utl_get_component_fromtemplate('visual_p100_erp',lf);
if cfg.debug
    p1_comp(1).signal{1}.peakAmplitudeSlope = 0;
    p1_comp(2).signal{1}.peakAmplitudeSlope = 0;
    plot_signal_fromclass(p1_comp(1).signal{1}, epochs);
    
end

p3_comp = [];

p3_comp.signal{1} = struct();
p3_comp.signal{1}.peakLatency = 350;      % in ms, starting at the start of the epoch
p3_comp.signal{1}.peakWidth = 400;        % in ms
p3_comp.signal{1}.peakAmplitude = -3;      % in microvolt
p3_comp.signal{1}.peakLatencyDv = 100; % reaction time?
p3_comp.signal{1}.peakAmplitudeDv = 3;
p3_comp.signal{1} = utl_check_class(p3_comp.signal{1}, 'type', 'erp');
p3_comp.signal{2} = noise;

if cfg.debug
        plot_signal_fromclass(p3_comp.signal{1}, epochs);
end

% somewhere deep brain?
p3_comp.source= lf_get_source_nearest(lf, [0 -40 -25]);
if cfg.debug
    %%
    p3_comp.source= lf_get_source_nearest(lf, [0 -40 0]);
    plot_source_location(p3_comp.source, lf, 'mode', '3d');
end

p3_comp.orientation = utl_get_orientation_pseudotangential(p3_comp.source,lf);
p3_comp.orientationDv = [0 0 0];

p3_comp.orientation = [0.03,-0.16,1];
if cfg.debug
    %%
    for k = 1:10
        p3_comp.orientation = utl_get_orientation_random(1);
        plot_source_projection(p3_comp.source, lf, 'orientation', p3_comp.orientation,'orientedonly',1);
        %     title(vis.orientation)
    end
end


if cfg.noise_orient
    n1_comp.orientation = n1_comp.orientation + rand(1,3)*0.1;
    p3_comp.orientation = p3_comp.orientation + rand(1,3)*0.1;
    randorient =  rand(1,3)*0.1;
    p1_comp(1).orientation =  p1_comp(1).orientation+ randorient;
    p1_comp(2).orientation =  p1_comp(2).orientation-randorient;
end
%%
% noise sources
random_sources = lf_get_source_spaced(lf, cfg.noise_components, 25);
random_comp = utl_create_component(random_sources, noise, lf);

%%
epochs = struct();
epochs.n = cfg.n_epochs;             % the number of epochs to simulate
epochs.srate = cfg.srate;        % their sampling rate in Hz
epochs.length = 1000*cfg.epochlength;
EEG=[];
for fn = {'p1','random','p3','n1'}
    eval(sprintf('comp = %s_comp;',fn{1}))
    data = generate_scalpdata(comp, lf, epochs,'showprogress',0);
    EEG.(fn{1}) = utl_create_eeglabdataset(data, epochs, lf, 'marker', 'event1');
end

