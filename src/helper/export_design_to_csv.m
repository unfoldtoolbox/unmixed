function export_design_to_csv(input,EEG,varargin)
%asda
cfg = finputcheck(varargin,...
    {'out_fn','string',[],'';... # outfilename
    'channel','integer',[],[];... # which channel
    },'mode','error');

assert(isfield(EEG.unmixed,'uf_fixef'),'uf_fixef is missing')


assert(isfield(EEG.unmixed.uf_fixef,'Xdc'),'No time-expanded designmatrix found');



    function [t] = makeTableXdc(uf_struct)
        assert(length(size(input{1}.data))==2)
        assert(all(cellfun(@(x)size(x.data,1)>=cfg.channel,input)));
        data_y = cellfun(@(x)squeeze(x.data(cfg.channel,:)),input,'UniformOutput',0);
        
        y = double(cat(2,data_y{:})');
        t = table();
        tmp_times = cellfun(@(x)x.times,input,'UniformOutput',0)
        tmp_times = cat(2,tmp_times{:});
        for c = 1:length(uf_struct.colnames)
            ix = uf_struct.Xdc_terms2cols==c;
            Xdc_subset = full(uf_struct.Xdc(:,ix));
            t_subset = table();
            t_subset.Xdc =(Xdc_subset(:));
            
            t_subset.colnames = repmat(uf_struct.colnames(c),size(t_subset.Xdc));
            % have to temp it in order to rotate it
            tmp_tau = repmat(uf_struct.times,size(uf_struct.Xdc,1),1);
            t_subset.tau= tmp_tau(:);
            t_subset.times = repmat(tmp_times,1,size(Xdc_subset,2))';
            t = [t;t_subset];
        end
        
        % Possibly give them a reasonable name at some point
        % fixefNames = strcat(uf_fixef.colnames(uf_fixef.Xdc_terms2cols),'_',sprintfc('%g',repmat(uf_fixef.times,1,length(uf_fixef.colnames))));
        t.y = repmat(y,size(uf_struct.Xdc_terms2cols,2),1);
        t.subject = repmat(repelem(1:length(input),cellfun(@(x)size(x.data,2),input))',size(uf_struct.Xdc_terms2cols,2),1);
    end

    
    for r = 1:length(EEG.unmixed.uf_ranef)
        t_r = makeTableXdc(EEG.unmixed.uf_ranef{r});
        writetable(t_r,sprintf('cache/%s_ranef-%s.csv',cfg.out_fn,EEG.unmixed.uf_ranef{r}.ranefgrouping),'WriteVariableNames',1)
    end
t = makeTableXdc(EEG.unmixed.uf_fixef);

% end
writetable(t,sprintf('cache/%s.csv',cfg.out_fn),'WriteVariableNames',1)

%% Export Random Effects

end
% t2 = unstack(t(t.colnames=="(Intercept)",:),'Xdc',{'tau'});
% t3 = unstack(t(t.colnames=="condA",:),'Xdc',{'tau'});
% all(t3{:,4:end} == Xdc_subset)