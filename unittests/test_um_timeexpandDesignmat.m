function test_um_timeexpandDesignmat()

input = test_um_generateTestData('randomItem',1); % add a randomeffect
%%
EEG_um = um_designmat(input,'eventtypes','sim','formula','y~1+b+(b|subject)+(1|trialnum)');
%%
EEG= um_timeexpandDesignmat(EEG_um,'timelimits',[-0.1,0.2]);

assert(length(intersect(1:10,EEG.unmixed.uf_ranef{1}.Zdc_level2cols)) == 10)
assert(length(intersect(1:60,EEG.unmixed.uf_ranef{2}.Zdc_level2cols)) == 60)

% check that Z is the same dimension as Xdc in time
assert(size(EEG.unmixed.uf_ranef{1}.Zdc,1) == size(EEG.unmixed.uf_ranef{1}.Xdc,1))
assert(size(EEG.unmixed.uf_ranef{2}.Zdc,1) == size(EEG.unmixed.uf_ranef{2}.Xdc,1))
