if 1 == 0
    rng(1)
    input = [];
    for k = 1:30
        %      input{k} = simulate_data_lmm_v2('noise',0.1,'noise_components',1,'srate',20,'n_events',64)
        input{k} = simulate_data_lmm_v2('noise',2,'noise_components',1,'srate',20,'n_events',64,...
            'b_p1_2x2',[10,5,0,0], ... % P1: Intercept, MainA, MainB, Inter - effect coded beta
            'u_p1_2x2',[5,5,0,0], ... %   P1: Subject variability
            'b_p3_2x2',[6,0,0,0], ...
            'u_p3_2x2',[3,0,0,0], ...
            'b_n1_2x2',[-6,0,0,0], ...
            'u_n1_2x2',[3,0,0,0], ...
            'u_p1_item',[5], ... % Item effect strength
            'u_p3_item',[5], ...
            'u_n1_item',[5], ...
            'simulationtype','ideal', ...
            'randomItem',1 ...
            );
        
    end
%     save('simulation_lmm.mat','input')
end
% load('simulation_lmm.mat')
%%

EEG = um_designmat(input,'inputGroupingName','subject',...
    'eventtypes','sim',...
    ... 'formula','y~1+(1|subject)+(1|stimulus)',...
    'formula','y~1+cat(condA)+(1+cat(condA)|subject)+(1|stimulus)',...
    'codingschema','effects');
%     'formula','y~1+cat(condA)*cat(condB)+(1+condA*condB|subject)+(1|stimulus)','codingschema','effects');

EEG= um_timeexpandDesignmat(EEG,'timelimits',[0,0.5]);
% EEG= um_timeexpandDesignmat(EEG,'timelimits',[0,0.03]);
%%
% Currently I recommend the bobyqa optimizer. Seems to be faster
tic
EEG = um_mmfit(EEG,input,'channel',1,'optimizer','bobyqa','covariance','FullCholesky'); % Todo: Directly read covariance from formula like lme4
toc

%%
umresult = um_condense(EEG);

%%
t =array2table(EEG.unmixed.uf_ranef{1}.Xdc);
t.image = EEG.unmixed.uf_ranef{2}.Xdc(:,2);
t.Properties.VariableNames = {'intercept','condA','subject','image'};
t(t.intercept == 0,:) = [];
t.y = EEG.unmixed.modelfit.y;
lmefit = fitlme(t,'y~1+condA+(1+condA|subject)+(1|image)');

assert(all(near(lmefit.Coefficients.Estimate,umresult.fixef.estimate,[],0.01)))
assert(all(near(lmefit.Coefficients.SE,umresult.fixef.se,[],0.01)))

% 
% figure,
% nTimeshifts = model.Psi.NumBlocks;
% subplot(2,nTimeshifts,1:nTimeshifts);
% plot(reshape(model.betaHat,nTimeshifts,[]),'o-')
% for k = 1:R
%     subplot(2,nTimeshifts,nTimeshifts + k);
%     imagesc(covmat{k})
%     caxis(prctile(tbl.Estimate,[1,99]))
% end
% 
% toc