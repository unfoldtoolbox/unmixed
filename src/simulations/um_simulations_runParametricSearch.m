function um_simulations_runParametricSearch
% Parametric search of parameters
% Because i don't have good intuitions what drives runtime of the mixed
% model, I do some conditional search of the parameterspace.
% My main finding so far (maybe unsurprisingly) is that sampling rate
% drives quite a strong effect.
% I would like to test basis functions now, to strongly reduce the amount
% of columns.
addpath('/home/common/matlab/fieldtrip/qsub')
local = 0;
cfg.folder = 'cache_realistic';
% Config Files for some runs, the first value is alywazs the default, i.e.
% we are not running a grid (combination) of all parameters
if cfg.folder == "cache-BobyqaAndQuasinewton"
    % first set of simulations I ran to test bobyqa & quasinewton.
    % bobyqa usually wins, but not by super much (mostly tested in the
    % diagonal condition)
    cfgSim =[];
    cfgSim.timelimits = {[-.1 0.5],[-.2 1]};
    cfgSim.noise = {20 1 40};
    cfgSim.u_noise = {20 1};
    cfgSim.srate = {20 100};
    cfgSim.n_events = {10,300}; %in s
    cfgSim.n_subjects = {50 20 90};
    cfgSim.u_p1_item = {0 10};
    cfgSim.optimizer = {'bobyqa','quasinewton'};
    cfgSim.covariance = {'diagonal','fullCholesky'};
elseif cfg.folder == "cache_realistic"
    
    cfgSim =[];
    cfgSim.timelimits = {[-.1 0.5]};
    cfgSim.noise = {5};%{20 10 40};
    cfgSim.u_noise = {10,40};%{20 40};
    cfgSim.srate = {50};
    cfgSim.n_events = {50}; %in s
    cfgSim.n_subjects = {50};%{50 10 100};
    cfgSim.u_p1_item = {0 10};
    cfgSim.optimizer = {'bobyqa'};
    cfgSim.covariance = {'diagonal','FullCholesky'};
    cfgSim.simulationtype = {'realistic'};
    cfgSim.overlaptype = {'lognormal'};
    cfgSim.folder = {cfg.folder};
    cfgSim.formula = {'y~1+condA + (1+condA|subject)'};
    cfgSim.b_p1_2x2={[10,5/2,0,0]};
    cfgSim.u_p1_2x2={[5,3,0,0],[1,0.5,0,0],[3,1.5,0,0], [10,6,0,0]};
    cfgSim.b_n1_2x2={[0,0,0,0]};
    cfgSim.u_n1_2x2={[0,0,0,0]};
    cfgSim.b_p3_2x2={[10,3/2,0,0]};
    cfgSim.u_p3_2x2={[5,3,0,0],[1,0.5,0,0],[3,1.5,0,0], [10,6,0,0]};
    cfgSim.channel = {52};

else
    
    cfgSim =[];
    cfgSim.timelimits = {[-.1 0.5],[-.2 1], [-.1 .3], [-1 2]};
    cfgSim.noise = {20 1 5 10 80 40 0};
    cfgSim.u_noise = {20 40 10 5 1 0};
    cfgSim.srate = {20 30 40 60 100};
    cfgSim.n_events = {10,20,40,80,160,320,600}; %in s
    cfgSim.n_subjects = {50 10,20,70, 90};
    cfgSim.u_p1_item = {0 10,20,5};
    cfgSim.optimizer = {'bobyqa'};
    cfgSim.covariance = {'diagonal','fullCholesky'};
end

grid_runParametricFunction(cfgSim,'um_simulations_grid',local)


    