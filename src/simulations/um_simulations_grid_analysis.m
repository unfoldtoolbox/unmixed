% This script analyses the results of
% um_simulations_runParametricSearch.m


cfg.folder = 'cache';
% cfg.folder = 'cache_BobyqaAndQausinewton';
%% Load Results
resultAll = struct();
d = dir(fullfile(cfg.folder, '*.mat'));
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
        if fn{1} == "timelimits"
            if r == 1
                resultAll.timelimits_min = [];
                resultAll.timelimits_max = [];
            end
            resultAll.timelimits_min = [resultAll.timelimits_min tmp{1}(1)];
            resultAll.timelimits_max = [resultAll.timelimits_max tmp{1}(2)];
        end
    end
end


resultAll.loglikHat = arrayfun(@(x)x.unmixed.modelfit.loglikHat,resultAll.model);
resultAll.sizeXdc1 = cellfun(@(x)x(1),resultAll.sizeXdc);
resultAll.sizeZdc1 = cellfun(@(x)x(1),resultAll.sizeZdc);
%%
% convert to table
str = 't = table(';
for k = fieldnames(resultAll)'
    str = [str 'resultAll.' k{1} ''','];
end
str = [ str '''VariableNames'',fieldnames(resultAll));'];
eval(str); % t is now a nice table
% sort
t = removevars(t,{'timelimits','model'});
t = sortrows(t,{'optimizer','timelimits_min','covariance','n_events','noise','u_noise','srate','n_subjects','u_p1_item'});



if cfg.folder == "cache_BobyqaAndQuasinewton"
    %% How much (%) is bobyqa faster / slower?
    figure,
    quasi = t{size(t,1)/2+1:end,'timing'};
    boby = t{1:size(t,1)/2,'timing'};
    histogram((boby - quasi)./((boby +quasi)./2)*100,40)
    vline(0,'k:')
    xlabel('% change, neg => bobyqa faster')
    % %% Eval loglike of two optimizer
    % This should be pretty much 0.
    t{1:size(t,1)/2,'loglikHat'} -        t{size(t,1)/2+1:end,'loglikHat'}
end

%% Go over these result fields
for fn = {'noise','u_noise','srate','n_events','n_subjects','covariance','u_p1_item'}
    if fn{1} =="optimizer"
        continue
    end
    %%
    %subset loop
    % find the
    ix = ones(size(resultAll.events)); % start with all
    for fn2 = fieldnames(cfgSim)' % go over all fields
        display(fn2{1})
        if strcmp(fn{1},fn2{1}) || fn2{1} =="optimizer"
            continue
        end
        
        switch fn2{1}
            case 'covariance'
                ix = ix & (strcmp(cfgSim.covariance{1},resultAll.covariance));
            otherwise
                if iscell(resultAll.(fn2{1}))
                    ix = ix & (cellfun(@(x)all(x == cfgSim.(fn2{1}){1}),resultAll.(fn2{1})));
                else
                    ix = ix & (cfgSim.(fn2{1}){1}==resultAll.(fn2{1}));
                end
        end
        
        
        
    end
    %%
    figure
    g = gramm('x',resultAll.(fn{1}),'y',resultAll.timing/60,'color',resultAll.optimizer,'subset',ix);
    g.geom_point();
    g.set_names('x',fn{1},'y','time in s');
    %g.stat_smooth('geom','line' );
    g.draw()
end
%%
figure
g = gramm('x',resultAll.n_subjects,'y',resultAll.timing/60,'color',resultAll.n_events, 'subset');
g.geom_point();
%g.stat_smooth('geom','line' );
g.draw()

%% Loglikely difference
figure
g = gramm('x',resultAll.sizeXdc1,'y',resultAll.timing/60,'color',resultAll.optimizer);
g.geom_point()
g.draw()

