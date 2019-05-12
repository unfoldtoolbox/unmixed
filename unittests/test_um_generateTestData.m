function [input,cfgSim] = test_um_generateTestData(varargin)
cfgSim = [];
cfgSim.subject = 10;
cfgSim.srate = 10;
cfgSim.datasamples = cfgSim.srate*60;
input = {};
for k = 1:cfgSim.subject
    input{k} = simulate_data_lmm('noise',5,'srate',cfgSim.srate,'datasamples',cfgSim.datasamples,'basis','dirac',varargin{:});
end

