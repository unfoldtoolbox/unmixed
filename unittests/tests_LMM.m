
EEG = simulate_data_lmm();
EEG = uf_designmat(EEG,'eventtypes','sim','formula','y~1+sub_slope');
EEG = uf_timeexpandDesignmat(EEG,'timelimits',[-0.5,1.5]);
EEG = uf_glmfit(EEG);
%%

