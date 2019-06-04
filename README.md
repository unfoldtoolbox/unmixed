# Disclaimer:
This extension is still being developed and needs optimization for speed and additional tests. **Use at your own risk!**. Probably best to contact me prior to using it to get to know the quirks :-)


# Unfold + Mixed Models = unmixed


Deconvolution, non linear modeling and Mixed Modeling toolbox

Specify ```y~1 + A + (1 + A|subject) + (1+C|image)``` Item and Subject Effects!

# Install
```
git clone https://github.com/behinger/unmixed
git submodule update --init --recursive
```

# Minimal Example:

```
rng(1)
input = [];
for k = 1:30
    % you can also add item effects by the adding the flag 'randomItem', 1
    % in simulate_data_lmm. The betas / thetas are currently hardcoded because I'm lazy
    input{k} = simulate_data_lmm('noise',10,...
        'srate',50,'datasamples',600*50,...
        'basis','dirac');
end


EEG = um_designmat(input,'eventtypes','fixation','formula','y~1+A+(1+A|subject)');

EEG= um_timeexpandDesignmat(EEG,'timelimits',[-0.1,0.5]);

% Currently I recommend the bobyqa optimizer. Seems to be faster
model_fv = um_mmfit(EEG,input,'channel',1,'optimizer','bobyqa','covariance','FullCholeksy'); % Todo: Directly read covariance from formula like lme4

```
