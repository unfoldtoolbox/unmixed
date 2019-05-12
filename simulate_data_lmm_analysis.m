%%
tic
EEG = simulate_data_lmm_v2('noise_components',1,'noise',0.01);
toc
EEG=uf_designmat(EEG,'formula','y~1+cat(condA)*cat(condB)','eventtypes','sim','codingschema','effects');
EEG = uf_timeexpandDesignmat(EEG,'timelimits',[0,simCFG.epochlength]);
EEG = uf_glmfit(EEG);
EEG_epoch = uf_epoch(EEG,'timelimits',[0,simCFG.epochlength]);
EEG_epoch = uf_glmfit_nodc(EEG_epoch);
ufresult = uf_condense(EEG_epoch);
%% Manual ERP
figure,
findevent= @(fn,val)[EEG.event.(fn)]==val;

ploterp = @(chan)plot(EEG_epoch.times,...
    [mean(squeeze(EEG_epoch.data(chan,:,findevent('condA',0))),2),...
    mean(squeeze(EEG_epoch.data(chan,:,findevent('condA',1))),2)]);
figure,ploterp(44),title(EEG_epoch.chanlocs(44).labels)
figure,ploterp(63),title(EEG_epoch.chanlocs(63).labels)
figure,ploterp(30),title(EEG_epoch.chanlocs(30).labels)
%% Visualize Deconv Beta

%% Within-Subject T-Test using LIMO
[m,dfe,ci,sd,n,t,p] = limo_ttest(2,EEG_epoch.data(:,:,findevent('condA',0)),EEG_epoch.data(:,:,findevent('condA',1)),0.05);
%% simple imagesc plots
figure,subplot(2,1,1),imagesc(EEG.times,1:64,mean(EEG.data,3));box off
subplot(2,1,2),imagesc(EEG.times,1:64,m),box off
% Show data
figure
imagesc(m)

figure
imagesc(t)
%% Prepare statistics


%%
ept_tfce_nb = ept_ChN2(EEG.chanlocs,0); % the 1 to plot
tfce_res = ept_TFCE(...
    permute(EEG_epoch.data(:,:,findevent('condA',0)),[3 1 2]),...
    permute(EEG_epoch.data(:,:,findevent('condA',1)),[3,1,2]),EEG.chanlocs,'nperm',800,'flag_save',0,'ChN',ept_tfce_nb);

%%
t_h0 = nan([size(EEG.data,1),size(EEG.data,2),800]);
for b = 1:800
    fprintf('%i/%i\n',b,800)
    r_perm      = randperm(epochs.n*2); % Consider using Shuffle mex here (50-85% faster)...
    
    nData       = EEG.data(:,:,r_perm);
    sData{1}    = nData(:,:,1:epochs.n);
    sData{2}    = nData(:,:,(epochs.n+1):(epochs.n*2));
    
    t_h0(:,:,b) = (mean(sData{1},3)-mean(sData{2},3))./sqrt(var(sData{1},[],3)/epochs.n+var(sData{2},[],3)/epochs.n);
    
    
end
p_h0 = tpdf(t_h0,199);
%% Generate Cluster Permutation
limo_nb = limo_neighbourdist(EEG,50);
[limomask,cluster_p,max_th] = limo_clustering(t.^2,p,t_h0.^2,p_h0,struct('data',struct('chanlocs',EEG.chanlocs','neighbouring_matrix',limo_nb)),2,0.05,0);
cluster_p(isnan(cluster_p(:))) = 1;
% F_tfce = limo_tfce(2,t.^2,limo_nb);

%% show significance methods
figure,
for k = 0:4
    switch k
        case 0
            plot_p = p;
            mask = plot_p<0.05;
        case 1
            plot_p = bonf_holm(double(p),0.05);
            mask = plot_p<0.05;
        case 2
            [~,~,~,plot_p] = fdr_bh(p,0.05,'pdep','no'); % see groppe 2011 why pdep is fine
            mask = plot_p<0.05;
            
        case 3
            plot_p = cluster_p;
            mask = plot_p<0.05;
            
        case 4
            plot_p = tfce_res.P_Values;
            mask = plot_p<0.05;
            
            
    end
    plot_data= plot_p;
    plot_data = log10(plot_data);
    subplot(5,1,(k)+1)
    eegvis_imagesc(m,t,'chanlocs',EEG.chanlocs,'times',EEG.times,...
        'contour',1,'clustermask',mask,'figure',0,'xlabel',0,'colorbar',k == 4)
%     caxis([-7,0])
%     if k ~=4
%         axis off
%     end
%     box off
    plot_data(~mask) = 1;
    box off
end

%%
hA = plot_topobutter(mean(EEG.data,3),EEG.times,EEG.chanlocs,'colormap',{'div'},'quality',40,'n_topos',12);

%%
figure
hA = plot_topobutter(cat(3,mean(EEG.data,3),m),EEG.times,EEG.chanlocs,'quality',32,'individualcolorscale','row','highlighted_channel',[44,63]);

%% Additional Blog Figures
plot_source_location([p100.source,vis.source,p3.source], lf, 'mode', '2d');




