function test_um_designmat
[input,cfgSim]= test_um_generateTestData();
%%

EEG_um = um_designmat(input,'eventtypes','sim','formula','y~-1+b+c+(-1+b*c|subject)+(1|d)');

check_struct(EEG_um,2) % 2 random effect

% we simulate 1 event per second
% with 10 subjects, that should be 60*10 events
% if this fails, it might be worthwhile to check if the simulation messed
% up!
assert(size(EEG_um.unmixed.uf_fixef.X,1) == cfgSim.datasamples*cfgSim.subject/cfgSim.srate);
assert(size(EEG_um.unmixed.uf_ranef{1}.X,1) == cfgSim.datasamples*cfgSim.subject/cfgSim.srate);


assert(size(EEG_um.unmixed.uf_fixef.X,2) == 2);
assert(size(EEG_um.unmixed.uf_ranef{1}.X,2) == 4);
assert(size(EEG_um.unmixed.uf_ranef{2}.X,2) == 2);

function check_struct(EEG,nRanef)
assert(isfield(EEG,'unmixed'))
assert(isfield(EEG.unmixed,'uf_fixef'))
assert(isfield(EEG.unmixed,'uf_ranef'))
assert(length(EEG.unmixed.uf_ranef)==nRanef)
end


try
EEG_um = um_designmat(input,'eventtypes','sim','formula','y~1+subject+c+(b|subject)+(1|d)');
error('this should have raised an error because subject was a fixed effect and grouping variable')
catch
    
end

end
