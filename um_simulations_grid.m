cfg = [];
cfg.timelimits = [-0.1,1.2];
cfg.noise = 20;
cfg.srate= 50;
cfg.datalength = 60;
cfg.nsubject = 25;
cfg.optimizer = 'fminunc';%{'fminunc','quasinewton','fminsearch','bobyqa'};

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
result.model = um_mmfit(EEG,input,'channel',1,'optimizer',cfg.optim);
result.timing = toc;

save(DataHash(cfg),'result')
