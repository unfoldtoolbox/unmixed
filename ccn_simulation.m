% if 1 == 0
    rng(1)
    input = [];
    for subj = 1:50
        %      input{k} = simulate_data_lmm_v2('noise',0.1,'noise_components',1,'srate',20,'n_events',64)
        input{subj} = simulate_data_lmm_v2('noise',20,'u_noise',10,'noise_components',5,'srate',20,'n_events',10,...
            'b_p1_2x2',[5,5,0,0], ... % P1: Intercept, MainA, MainB, Inter - effect coded beta
            'u_p1_2x2',[10,1,0,0], ... %   P1: Subject variability
            'b_p3_2x2',[0,0,0,0], ...
            'u_p3_2x2',[0,0,0,0], ...
            'b_n1_2x2',[0,0,0,0], ...
            'u_n1_2x2',[0,0,0,0], ...
            'u_p1_item',[0], ... % Item effect strength
            'u_p3_item',[0], ...
            'u_n1_item',[0], ...
            'simulationtype','ideal_hanning', ...
            'overlaptype','uniform',...
            'randomItem',1 ...
            );
        
    end
%     save('simulation_ccn_540events.mat','input')
% end
% tmp = load('simulation_ccn_540events.mat')
% input = tmp.input;
%%
% figure,
% for s = 1:30
%     subplot(5,6,s)
% plot(input{s}.times/1000,input{s}.data(1,:))
% title(sprintf('noiselevel: %.0f',input{s}.sim.noiselevel))
% box off
% xlabel('Time [s]')
% ylabel('Amplitude [µV]')
% end
%%
%     'formula','y~1+cat(condA)*cat(condB)+(1+condA*condB|subject)+(1|stimulus)','codingschema','effects');
%% Design
cfgDesign = struct();
cfgDesign.inputGroupingName='subject';
cfgDesign.eventtypes= 'sim';
% cfgDesign.formula= 'y~1+cat(condA)+cat(condB)+(1+cat(condA)+cat(condB)|subject)';
% cfgDesign.formula= 'y~1 + cat(condA) + (1+cat(condA)|subject)';
% cfgDesign.formula = 'y~1 + cat(condA)+ (1+cat(condA)|subject)+(1|stimulus)';
cfgDesign.formula = 'y~1+condA+(1+condA|subject)';
% cfgDesign.formula = 'y~1+condA+(1+condA|subject)+(1|stimulus)';
cfgDesign.codingschema = 'reference';

EEG = um_designmat(input,cfgDesign);

EEG= um_timeexpandDesignmat(EEG,'timelimits',[-.1,0.5]);

EEG = um_mmfit(EEG,input,'channel',1,'optimizer','bobyqa','covariance','Diagonal');

umresult = um_condense(EEG);
% umresult2 = um_condense(EEG2);

%%
figure
g =gramm('x',umresult.times,'y',umresult.fixef.estimate,'color',umresult.fixef.names,...
    'ymin',umresult.fixef.estimate+-2*umresult.fixef.se,...
    'ymax',umresult.fixef.estimate+2*umresult.fixef.se);
g.geom_interval();
g.set_names('x',EEG.unmixed.formula)
g.axe_property('YLim',[-5,15])
g.draw();
% title(EEG.unmixed.formula)

%% Comparison normal ERP
d2nd = [];
for subj = 1:length(input)
    EEGsing = uf_designmat(input{subj},'formula',char(EEG.unmixed.formulaFixef),'eventtypes','sim',...
        'codingschema','reference');
    EEGsing = uf_timeexpandDesignmat(EEGsing,'timelimits',[-.1,.5]);
    EEGsing = uf_glmfit(EEGsing,'channel',1);
    EEGsing = uf_epoch(EEGsing,'timelimits',[-.1,.5]);
    EEGsing = uf_glmfit_nodc(EEGsing);
    ufresult = uf_condense(EEGsing);
    if isempty(d2nd)
        d2nd = ufresult;
        d2nd.subject = subj;
    else
        d2nd.unfold(end+1) = ufresult.unfold;
        d2nd.beta(:,:,:,end+1) = ufresult.beta;
        d2nd.beta_nodc(:,:,:,end+1) = ufresult.beta_nodc;
        d2nd.subject(end+1) = subj;
    end
end

%% Two Stage Plot
figure
estimate = squeeze(d2nd.beta(1,:,:,:));
se = std(esti2mate,[],3)/sqrt(30);
m = mean(estimate,3);
g =gramm('x',d2nd.times','y',m',...
    'color',{ufresult.param.name},...
    'ymin',m'+-2.045*se',...
    'ymax',m'+2.045*se');
g.geom_interval();
g.set_names('x','two stage')
g.axe_property('YLim',[-5,15])
g.draw();
%% No Deconvolution 2 stage
figure
estimate = squeeze(d2nd.beta_nodc(1,:,:,:));
se = std(estimate,[],3)/sqrt(30);
m = mean(estimate,3);
g =gramm('x',d2nd.times','y',m',...
    'color',{ufresult.param.name},...
    'ymin',m'+-2.045*se',...
    'ymax',m'+2.045*se');
g.geom_interval();
g.set_names('x','no deconv two stage')
g.axe_property('YLim',[-5 15])
g.draw();

%% Complete Pooling Plot
mdl = fitlm(EEG.unmixed.modelfit.X,EEG.unmixed.modelfit.y);
%%
figure
m = reshape(mdl.Coefficients.Estimate(2:end),[],2);
se = reshape(mdl.Coefficients.SE(2:end),[],2);
g =gramm('x',d2nd.times','y',m',...
    'color',{ufresult.param.name},...
    'ymin',m'+-2.045*se',...
    'ymax',m'+2.045*se');
g.geom_interval();
g.set_names('x','complete pooling')
g.axe_property('YLim',[-5,15])
g.draw();
%% random effect
figure
plot(umresult.times,umresult.ranef(1).covmat(:,1,1))
hold all
try
plot(umresult.times,umresult.ranef(1).covmat(:,2,2))
catch
end
% 
% plot(umresult.times,umresult.ranef(2).covmat(:,1))

%% No Noise Plot

    rng(1)
    input = [];
    for subj = 1:50
        %      input{k} = simulate_data_lmm_v2('noise',0.1,'noise_components',1,'srate',20,'n_events',64)
        input{subj} = simulate_data_lmm_v2('noise',0,'u_noise',0,'noise_components',5,'srate',20,'n_events',10,...
            'b_p1_2x2',[5,5,0,0], ... % P1: Intercept, MainA, MainB, Inter - effect coded beta
            'u_p1_2x2',[10,1,0,0], ... %   P1: Subject variability
            'b_p3_2x2',[0,0,0,0], ...
            'u_p3_2x2',[0,0,0,0], ...
            'b_n1_2x2',[0,0,0,0], ...
            'u_n1_2x2',[0,0,0,0], ...
            'u_p1_item',[0], ... % Item effect strength
            'u_p3_item',[0], ...
            'u_n1_item',[0], ...
            'simulationtype','ideal_hanning', ...
            'overlaptype','uniform',...
            'randomItem',1 ...
            );
        
    end
%%
    d2nd = [];
for subj = 1:length(input)
    EEGsing = uf_designmat(input{subj},'formula',char(EEG.unmixed.formulaFixef),'eventtypes','sim',...
        'codingschema','effects');
    EEGsing = uf_timeexpandDesignmat(EEGsing,'timelimits',[-.1,.5]);
    EEGsing = uf_glmfit(EEGsing,'channel',1);
    EEGsing = uf_epoch(EEGsing,'timelimits',[-.1,.5]);
    EEGsing = uf_glmfit_nodc(EEGsing);
    ufresult = uf_condense(EEGsing);
    if isempty(d2nd)
        d2nd = ufresult;
        d2nd.subject = subj;
    else
        d2nd.unfold(end+1) = ufresult.unfold;
        d2nd.beta(:,:,:,end+1) = ufresult.beta;
        d2nd.beta_nodc(:,:,:,end+1) = ufresult.beta_nodc;
        d2nd.subject(end+1) = subj;
    end
end
%%
figure,
plot(d2nd.times,squeeze(d2nd.beta(:,:,1,:)))
hold on
plot(d2nd.times,mean(squeeze(d2nd.beta(:,:,2,:)),2),'-k','LineWidth',3)
