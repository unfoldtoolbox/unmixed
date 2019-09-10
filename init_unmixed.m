try
    addpath('lib/unfold')
catch
    fprintf('could add unfold, assuming its already added\n')
end
addpath('src/helper')
addpath('src/simulations')
addpath('src/um_toolbox')
addpath('src')
addpath('lib/bobyqa')
addpath('lib/eeglab')
addpath('local')
addpath(genpath('lib/sereega/'))
init_unfold
