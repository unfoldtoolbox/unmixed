function um_simulations_grid(varargin)
if 1 == 0
    %% donders specific Grid start
   addpath('/home/common/matlab/fieldtrip/qsub')
   cfg =[];
   cfg.timelimits = {[-.1 0.5],[-.1 1.2]};
   cfg.noise = {20 10 5 1 40};
   cfg.srate = {20 50 150};
   cfg.datalength = {600};
   cfg.nsubject = {10 25 50};
   optimizer = {'fminunc','quasinewton','fminsearch'};
%    optimizer = {'bobyqa'};
   for optim= optimizer
       for fn = fieldnames(cfg)'
           for fnDefault = fieldnames(cfg)'
               eval(sprintf('%s = cfg.%s{1};',fnDefault{1},fnDefault{1}))
           end
           
           for k = 1:length(cfg.(fn{1}))
               eval(sprintf('%s = cfg.%s{%i};',fn{1},fn{1},k))
               
               qsubfeval(@um_simulations_grid,timelimits,noise,srate,datalength,nsubject,optim{1},'memreq',2*1024^3,'timreq',60*2*60)
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
cfg = [];
cfg.timelimits = varargin{1};%[-0.1,1.2];
cfg.noise = varargin{2};%varargin{1};%20;
cfg.srate= varargin{3};%10;
cfg.datalength = varargin{4};%60;
cfg.nsubject = varargin{5};%25;
cfg.optimizer = varargin{6};%'fminunc';%
cfg.covariance = varargin{7}; %FullCholesky,'Diagonal'

rng(1)
input = [];
for k = 1:cfg.nsubject
    input{k} = simulate_data_lmm('noise',cfg.noise,...
        'srate',cfg.srate,'datasamples',cfg.datalength*cfg.srate,...
        'basis','dirac');
end


EEG = um_designmat(input,'eventtypes','sim','formula','y~1+b+(1+b|subject)');
EEG= um_timeexpandDesignmat(EEG,'timelimits',cfg.timelimits);


result = cfg;
result.events = length(EEG.event);
result.sizeXdc = size(EEG.unmixed.uf_fixef.Xdc);
result.sizeZdc = size(EEG.unmixed.uf_ranef{1}.Zdc);
tic
result.model = um_mmfit(EEG,input,'channel',1,'optimizer',cfg.optimizer);
result.timing = toc;

save(fullfile('cache',DataHash(cfg)),'result')
