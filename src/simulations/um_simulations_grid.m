function um_simulations_grid(varargin)
[cfg,residual_varargin] = finputcheck(varargin,...
    {'timelimits','real',[],[-.2 1.2]
    'noise','real',[],20
    'u_noise','real',[],20
    'srate','real',[],50
    'n_events','real',[],10
    'n_subjects','real',[],50
    'u_p1_item', 'real',[], [0]
    'optimizer','string',[],'bobyqa'
    'covariance','string',[],'diagonal'
    },'mode','ignore');



if 1 == 0
    %% donders specific Grid start
   addpath('/home/common/matlab/fieldtrip/qsub')
   cfg =[];
   cfg.timelimits = {[-.2 1.2],[-.1 2]};
   cfg.noise = {20 1 40};
   cfg.u_noise = {20 1};
   cfg.srate = {50 250};
   cfg.n_events = {10,300}; %in s
   cfg.n_subjects = {50 20 90};
   cfg.u_p1_item = {[0 10]};
   cfg.optimizer = {'bobyqa','quasinewton'};
   cfg.covariance = {'diagonal','fullCholesky'};
   for optim= cfg.optimizer
       for fn = fieldnames(cfg)'
           if fn == "optimizer"
               continue
           end
           % in order to not check out the complete grid of all
           % combinations, we set a default and vary the parameter in one
           % dimension only
           
           % set the first value as the default
           for fnDefault = fieldnames(cfg)'
               eval(sprintf('%s = cfg.%s{1};',fnDefault{1},fnDefault{1}))
           end
           
           % calculate all variations in the other dimensions
           for k = 1:length(cfg.(fn{1}))
               eval(sprintf('%s = cfg.%s{%i};',fn{1},fn{1},k))
               
%                qsubfeval(@um_simulations_grid,'timelimits',timelimits,'noise',noise,...
%                    'u_noise',u_noise,'srate',srate,'n_events',n_events,'n_subjects',n_subjects,...
%                    'optimizer',optim{1},'covariance',covariance,'u_p1_item',u_p1_item,'memreq',2*1024^3,'timreq',60*20*60)
           end
       end
   end
   %%
   if 1 == 0
       %% Load Results
       resultAll = struct();
       d = dir('cache/*.mat');
       for r = 1:length(d)
           tmp = load(fullfile('cache',d(r).name));
           result = tmp.result;
           if r == 1
               resultAll = result;
           end
           for fn = fieldnames(result)'
                    
                   tmp = result.(fn{1});
                   if length(tmp)>1
                       tmp = {tmp};
                   end
                   if r == 1
                       
                       resultAll.(fn{1}) = tmp;
                   else
                   resultAll.(fn{1}) = [resultAll.(fn{1}) tmp];
                   end
               
           end
       end
       
        %%
        resultAll.loglikHat = arrayfun(@(x)x.loglikHat,resultAll.model);
        resultAll.sizeXdc1 = cellfun(@(x)x(1),resultAll.sizeXdc);
        resultAll.sizeZdc1 = cellfun(@(x)x(1),resultAll.sizeZdc);
        %%
       figure
       g = gramm('x',resultAll.sizeZdc1+resultAll.sizeXdc1,'y',resultAll.timing/60,'color',resultAll.optimizer);g.geom_point();g.stat_smooth('geom','line' );g.draw()

        %% Plot Results
       %x = cellfun(@(x)x(2),resultAll.timelimits);
       %ix = resultAll.noise == 20& resultAll.srate == 20 & resultAll.nsubject == 10;
       
%         x = resultAll.noise;
%         ix = cellfun(@(x)x(2) == 0.5,resultAll.timelimits) & resultAll.srate == 20 & resultAll.nsubject == 10;
       
       %x = resultAll.nsubject;
       %ix = cellfun(@(x)x(2) == 0.5,resultAll.timelimits) & resultAll.srate == 20 & resultAll.noise == 20;
       
       
%       x = resultAll.srate;
%        ix = cellfun(@(x)x(2) == 0.5,resultAll.timelimits) & resultAll.noise == 20 & resultAll.nsubject == 10 ;
       
        ix = ix & resultAll.datalength == 60;
%        figure
%        g = gramm('x',x(ix),'y',resultAll.loglikHat(ix),'color',resultAll.optimizer(ix));g.geom_point();g.draw()
       figure
       g = gramm('x',x(ix),'y',resultAll.timing(ix)/60,'color',resultAll.optimizer(ix));g.geom_point();g.draw()
   end
 
    
end



rng(1)
input = [];
for k = 1:cfg.n_subjects
    cfg.randomItem = 1;
    cfg.simulationtype  = 'ideal_hanning';
    cfg.overlaptype = 'lognormal';
    cfg.b_p1_2x2=[5,5,0,0];
    cfg.u_p1_2x2=[10 1 0 0];
    cfg.b_p3_2x2=[0,0,0,0];
    cfg.u_p3_2x2=[0,0,0,0];
    cfg.b_n1_2x2=[0,0,0,0];
    cfg.u_n1_2x2=[0,0,0,0];
    cfg.u_p3_item=[0];
    cfg.u_n1_item=[0];
    cfg.noise_components = 5; % unused if ideal_hanning is used
    input{k} = simulate_data_lmm_v2(cfg);
end


EEG = um_designmat(input,'eventtypes','sim','formula','y~1+condA+condB+(1+condA|subject)');
EEG= um_timeexpandDesignmat(EEG,'timelimits',cfg.timelimits);


result = cfg;
result.events = length(EEG.event);
result.sizeXdc = size(EEG.unmixed.uf_fixef.Xdc);
result.sizeZdc = size(EEG.unmixed.uf_ranef{1}.Zdc);
tic
result.model = um_mmfit(EEG,input,'channel',1,'optimizer',cfg.optimizer);
result.timing = toc;
toc
save(fullfile('cache',DataHash(cfg)),'result','cfg')
