% % if 1 == 0
%     rng(1)
%     input = [];
%     for subj = 1:50
%         %      input{k} = simulate_data_lmm_v2('noise',0.1,'noise_components',1,'srate',20,'n_events',64)
%         input{subj} = simulate_data_lmm_v2('noise',5,'u_noise',10,'noise_components',5,'srate',50,'n_events',50,...
%             'b_p1_2x2',[10,5/2,0,0],... % P1: Intercept, MainA, MainB, Inter - effect coded beta
%             'u_p1_2x2',[3,3,0,0],... %   P1: Subject variability
%             'b_n1_2x2',[0,0,0,0],...
%             'u_n1_2x2',[0,0,0,0],...
%             'b_p3_2x2',[10,3/2,0,0],...
%             'u_p3_2x2',[5,0,0,0],...
%             'u_p1_item',[0], ... % Item effect strength
%             'u_p3_item',[0], ...
%             'u_n1_item',[0], ...
%             'simulationtype','realistic', ...
%             'overlaptype','lognormal',...
%             'randomItem',1 ...
%             );
%         
%     end
% save('simulation_ccn_poster.mat','input')
% % end
%%
% tmp = load('simulation_ccn_poster.mat')
% input = tmp.input;

% %%
% figure,
% for s = 1:10
%     subplot(5,6,s)
% plot(input{s}.times/1000,input{s}.data(30,:))
% title(sprintf('noiselevel: %.0f',input{s}.sim.noiselevel))
% box off
% xlabel('Time [s]')
% ylabel('Amplitude [uV]')
% end
%%
%     'formula','y~1+cat(condA)*cat(condB)+(1+condA*condB|subject)+(1|stimulus)','codingschema','effects');
%% Design
%%
% %%
% cfgDesign = struct();
% cfgDesign.inputGroupingName='subject';
% cfgDesign.eventtypes= 'sim';
% % cfgDesign.formula= 'y~1+cat(condA)+cat(condB)+(1+cat(condA)+cat(condB)|subject)';
% % cfgDesign.formula= 'y~1 + cat(condA) + (1+cat(condA)|subject)';
% % cfgDesign.formula = 'y~1 + cat(condA)+ (1+cat(condA)|subject)+(1|stimulus)';
% cfgDesign.formula = 'y~1+condA+(1+condA|subject)';
% % cfgDesign.formula = 'y~1+condA+(1+condA|subject)+(1|stimulus)';
% cfgDesign.codingschema = 'reference';
% 
% EEG = um_designmat(input,cfgDesign);
% 
% EEG= um_timeexpandDesignmat(EEG,'timelimits',[-.1,0.5]);
% 
% EEG = um_mmfit(EEG,input,'channel',52,'optimizer','bobyqa','covariance','Diagonal');
%  
% umresult = um_condense(EEG);
% % umresult2 = um_condense(EEG2);

%%
tmp = load('cache_realistic/e07d36a071ed2cb09ce84f69a47dbba5.mat')
input = tmp.input;
EEG = tmp.result.model;
umresult = um_condense(EEG)

%% Comparison normal ERP
d2nd = [];
for subj = 1:length(input)
    EEGsing = uf_designmat(input{subj},'formula',char(EEG.unmixed.formulaFixef),'eventtypes','sim',...
        'codingschema','reference');
    EEGsing = uf_timeexpandDesignmat(EEGsing,'timelimits',[-.1,.5]);
    EEGsing = uf_glmfit(EEGsing)%,'channel',21);
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
%% Mixed model Plot
plotYlim = [-0.75,1.75];
figure
g =gramm('x',umresult.times,'y',umresult.fixef.estimate,'color',umresult.fixef.names,...
    'ymin',umresult.fixef.estimate+-2*umresult.fixef.se,...
    'ymax',umresult.fixef.estimate+2*umresult.fixef.se);
g.geom_interval();
g.set_names('x',EEG.unmixed.formula)
g.axe_property('YLim',plotYlim)
g.draw();
% title(EEG.unmixed.formula)

%% Two Stage Plot
figure
estimate = squeeze(d2nd.beta(52,:,:,:));
se = std(estimate,[],3)/sqrt(length(input));
m = mean(estimate,3);
g =gramm('x',d2nd.times','y',m',...
    'color',{ufresult.param.name},...
    'ymin',m'+-2.045*se',...
    'ymax',m'+2.045*se');
g.geom_interval();
g.set_names('x','two stage')
g.axe_property('YLim',plotYlim)
g.draw();
%% No Deconvolution 2 stage
figure
estimate = squeeze(d2nd.beta_nodc(52,:,:,:));
se = std(estimate,[],3)/sqrt(50);
m = mean(estimate,3);
g =gramm('x',d2nd.times','y',m',...
    'color',{ufresult.param.name},...
    'ymin',m'+-2.045*se',...
    'ymax',m'+2.045*se');
g.geom_interval();
g.set_names('x','no deconv two stage')
g.axe_property('YLim',plotYlim)
g.draw();

%% Complete Pooling Plot
mdl = fitlm(EEG.unmixed.modelfit.X,EEG.unmixed.modelfit.y);
%
figure
m = reshape(mdl.Coefficients.Estimate(2:end),[],2);
se = reshape(mdl.Coefficients.SE(2:end),[],2);
g =gramm('x',d2nd.times','y',m',...
    'color',{ufresult.param.name},...
    'ymin',m'+-2.045*se',...
    'ymax',m'+2.045*se');
g.geom_interval();
g.set_names('x','complete pooling')
g.axe_property('YLim',plotYlim)
g.draw();
%% random effect
figure
plot(umresult.times,umresult.ranef(1).covmat(:,1,1))
hold all
try
plot(umresult.times,umresult.ranef(1).covmat(:,2,2))
catch
end
%%

estimate = squeeze(d2nd.beta(52,:,:,:));
two_se = std(estimate,[],3)/sqrt(length(input));
pool_se = reshape(mdl.Coefficients.SE(2:end),[],2);
mix_se = umresult.fixef.se';

figure,
ix = 1;
plot(pool_se(:,ix),mix_se(:,ix),'o'), hold on
plot(pool_se(:,ix),two_se(:,ix),'rx')

%%
figure,
ix = 1;
plot(d2nd.times,    pool_se(:,ix),'o-'), hold on, 
plot(d2nd.times,mix_se(:,ix),':x'),
plot(d2nd.times,two_se(:,ix),'-.x')
legend('pool','mixed','2stage')
box off
% 
% plot(umresult.times,umresult.ranef(2).covmat(:,1))
%% Save plots

FolderName = 'cache_realistic';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(get(FigHandle, 'Number'));
  set(0, 'CurrentFigure', FigHandle);
  export_fig(fullfile(FolderName, [FigName '.eps']),'-eps');
end
%% No Noise Plot
% 
%     rng(1)
%     input = [];
%     for subj = 1:50
%         %      input{k} = simulate_data_lmm_v2('noise',0.1,'noise_components',1,'srate',20,'n_events',64)
%         input{subj} = simulate_data_lmm_v2('noise',0,'u_noise',0,'noise_components',5,'srate',20,'n_events',10,...
%             'b_p1_2x2',[5,5,0,0], ... % P1: Intercept, MainA, MainB, Inter - effect coded beta
%             'u_p1_2x2',[10,1,0,0], ... %   P1: Subject variability
%             'b_p3_2x2',[0,0,0,0], ...
%             'u_p3_2x2',[0,0,0,0], ...
%             'b_n1_2x2',[0,0,0,0], ...
%             'u_n1_2x2',[0,0,0,0], ...
%             'u_p1_item',[0], ... % Item effect strength
%             'u_p3_item',[0], ...
%             'u_n1_item',[0], ...
%             'simulationtype','ideal_hanning', ...
%             'overlaptype','uniform',...
%             'randomItem',1 ...
%             );
%         
%     end
% %%
%     d2nd = [];
% for subj = 1:length(input)
%     EEGsing = uf_designmat(input{subj},'formula',char(EEG.unmixed.formulaFixef),'eventtypes','sim',...
%         'codingschema','reference');
%     EEGsing = uf_timeexpandDesignmat(EEGsing,'timelimits',[-.1,.5]);
%     EEGsing = uf_glmfit(EEGsing,'channel',52);
%     EEGsing = uf_epoch(EEGsing,'timelimits',[-.1,.5]);
%     EEGsing = uf_glmfit_nodc(EEGsing);
%     ufresult = uf_condense(EEGsing);
%     if isempty(d2nd)
%         d2nd = ufresult;
%         d2nd.subject = subj;
%     else
%         d2nd.unfold(end+1) = ufresult.unfold;
%         d2nd.beta(:,:,:,end+1) = ufresult.beta;
%         d2nd.beta_nodc(:,:,:,end+1) = ufresult.beta_nodc;
%         d2nd.subject(end+1) = subj;
%     end
% end
% %%
% figure,
% plot(d2nd.times,mean(d2nd.beta_nodc(52,:,1,:),4))
% hold on
% plot(d2nd.times,mean(sum(d2nd.beta_nodc(52,:,:,:),3),4),':')
% 
% plot(d2nd.times,mean(d2nd.beta(52,:,1,:),4))
% hold on
% plot(d2nd.times,mean(sum(d2nd.beta(52,:,:,:),3),4),':')