function tests_LMM
% This function calculates an ideal spike response, analyses it for one
% elektrode and compares the results it to a "perfect" matlab response

% The bobyqa optimizer is tested here

rng(1)
input = [];
for k = 1:30
    %%
    %      input{k} = simulate_data_lmm_v2('noise',0.1,'noise_components',1,'srate',20,'n_events',64)
     input{k} =...
    simulate_data_lmm_v2('noise',20,'noise_components',5,'srate',20,'n_events',96,...
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

%%

EEG = um_designmat(input,'inputGroupingName','subject',...
    'eventtypes','sim',...
    ... 'formula','y~1+(1|subject)+(1|stimulus)',...
    'formula','y~1+cat(condA)+(1+cat(condA)|subject)+(1|stimulus)',...
    'codingschema','reference');
%     'formula','y~1+cat(condA)*cat(condB)+(1+condA*condB|subject)+(1|stimulus)','codingschema','effects');

% EEG= um_timeexpandDesignmat(EEG,'timelimits',[0,0.5]);
%% Single Time point
EEG= um_timeexpandDesignmat(EEG,'timelimits',[0,0.03]);
tic
EEG = um_mmfit(EEG,input,'channel',1,'optimizer','bobyqa','covariance','Diagonal'); % Todo: Directly read covariance from formula like lme4
toc

umresult = um_condense(EEG);

% Recalculate everything completly in matlab
t =array2table(EEG.unmixed.uf_ranef{1}.Xdc);
t.image = EEG.unmixed.uf_ranef{2}.Xdc(:,2);
t.Properties.VariableNames = {'intercept','condA','subject','image'};
t(t.intercept == 0,:) = [];
t.y = EEG.unmixed.modelfit.y;
lmefit = fitlme(t,'y~1+condA+(1+condA|subject)+(1|image)');

assert(all(near(lmefit.Coefficients.Estimate,umresult.fixef.estimate)))
assert(all(near(lmefit.Coefficients.SE,umresult.fixef.se)))
%% Multiple (two) Time Points
EEG= um_timeexpandDesignmat(EEG,'timelimits',[0,0.15]);
tic
EEG = um_mmfit(EEG,input,'channel',1,'optimizer','quasinewton','covariance','Diagonal'); % Todo: Directly read covariance from formula like lme4
toc

umresult = um_condense(EEG);
%%
% Recalculate everything completly in matlab
ntime = length(EEG.unmixed.uf_fixef.times);
t =    array2table(EEG.unmixed.uf_ranef{1}.Xdc);
t = [t array2table(EEG.unmixed.uf_ranef{2}.Xdc(:,1+end-ntime:end),'VariableNames',string(char(65+randsample(25,ntime))))];

tmp = repmat({'intercept','condA','subject','image'},ntime,1);
tmp = cellfun(@(x,y)[x,char(65+y)],vertcat(tmp(:))',repmat(num2cell(1:length(EEG.unmixed.uf_fixef.times)),1,4),'UniformOutput',0);
t.Properties.VariableNames = [regexprep(tmp,'\.','')];
t(all(t{:,1:ntime} == 0,2),:) = [];
t.y = EEG.unmixed.modelfit.y;

%fixed effects
allterms = t.Properties.VariableNames(1:2*ntime);
basestring = strjoin(allterms(1:ntime*2),'+');
form = ['y~-1+',basestring];
    %random effects
for k = 1:ntime
allterms = t.Properties.VariableNames(k:ntime:4*ntime);
basestring = strjoin(allterms(1:2),'+');

form = [form '+(-1+',basestring,'|',allterms{3} ')+(-1+',allterms{1},'|', allterms{4} ')'];
end
lmefit = fitlme(t,form,'CovariancePattern',repmat({'Diagonal'},ntime*2,1));
%%
assert(all(near(lmefit.Coefficients.Estimate,umresult.fixef.estimate)))
assert(all(near(lmefit.Coefficients.SE,umresult.fixef.se)))
