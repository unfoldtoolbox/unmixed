function grid_runParametricFunction(cfgSim,functionName,local)


defaults = {};
for fnDefault = fieldnames(cfgSim)'
    defaults(end+1,:) = {fnDefault{1},cfgSim.(fnDefault{1}){1}}; 
end

fields = fieldnames(cfgSim);
for fn = fields'
    % load defaults
    currentSettings = defaults;
    %find which to replace
    ix = strcmp(defaults(:,1),fn{1});
    for k = 1:length(cfgSim.(fn{1}))
        % proceed with the all default only for the first field
        if k == 1 && ~strcmp(fields{1},fn{1})
            continue
        end
        % replace by the k-th entry of the current active setting
        currentSettings(ix,2) = cfgSim.(fn{1})(k);
        
        % generate string
        cmd = [];
        for l = 1:size(currentSettings,1)
            cmd = [cmd sprintf('''%s'',%s,',currentSettings{l,1},mat2str(currentSettings{l,2}))];
        end
        cmd = cmd(1:end-1); % remove that last ','
        if local == 0 
            % start on grid
            cmd = ['qsubfeval(@' functionName ',' cmd ',''memreq'',2*1024^3,''timreq'',60*20*60);'];
        else
            cmd = [functionName '(' cmd ');'];
        end
        eval(cmd);
        
    end
end