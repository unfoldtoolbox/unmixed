resultAll = struct();
d = dir(fullfile('cache_realistic', '*.mat'));
for r = 1:length(d)
    tmp = load(fullfile('cache_realistic',d(r).name));
    result = tmp.result;
    umresult = um_condense(result.model);
    %     figure
    %     g =gramm('x',umresult.times,'y',umresult.fixef.estimate,'color',umresult.fixef.names,...
    %         'ymin',umresult.fixef.estimate+-2*umresult.fixef.se,...
    %         'ymax',umresult.fixef.estimate+2*umresult.fixef.se);
    %     g.geom_interval();
    %     % g.set_names('x',EEG.unmixed.formula)
    %     g.axe_property('YLim',[-2,7])
    %
    %     g.draw();
    %     title(sprintf('p1:%i,p3:%i,form:%s',tmp.cfg.u_p1_2x2(1),tmp.cfg.u_p3_2x2(1),tmp.cfg.formula))
    %     set(gcf,'name',d(r).name)
    
    figure
    plot(umresult.times,umresult.ranef(1).covmat(:,1,1))
    hold all
    try
        plot(umresult.times,umresult.ranef(1).covmat(:,2,2))
    catch
    end
    set(gcf,'name',d(r).name)
    title(sprintf('p1:%.1f,%.1f,p3:%.1f,%.1f,form:%s',tmp.cfg.u_p1_2x2(1),tmp.cfg.u_p1_2x2(2),tmp.cfg.u_p3_2x2(1),tmp.cfg.u_p3_2x2(2),tmp.cfg.formula))
end

