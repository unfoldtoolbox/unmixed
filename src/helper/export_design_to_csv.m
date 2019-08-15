input = load('simulation_ccn_540events.mat');
input = input.input;
%%

cfgDesign = struct();
cfgDesign.inputGroupingName='subject';
cfgDesign.eventtypes= 'sim';
% cfgDesign.formula= 'y~1+cat(condA)+cat(condB)+(1+cat(condA)+cat(condB)|subject)';
% cfgDesign.formula= 'y~1 + cat(condA) + (1+cat(condA)|subject)';
% cfgDesign.formula = 'y~1 + cat(condA)+ (1+cat(condA)|subject)+(1|stimulus)';
cfgDesign.formula = 'y~1+condA+(1+condA|subject)';
% cfgDesign.formula = 'y~1+condA+(1+condA|subject)+(1|stimulus)';
cfgDesign.codingschema = 'effects';

EEG = um_designmat(input,cfgDesign);

EEG= um_timeexpandDesignmat(EEG,'timelimits',[-.1,0.5]);

csvwrite('Xdc.csv',full(EEG.unmixed.uf_fixef.Xdc))

data_y = cellfun(@(x)squeeze(x.data(1,:,:)),input,'UniformOutput',0);

y = double(cat(2,data_y{:})');
t = array2table(full(EEG.unmixed.uf_fixef.Xdc));
% Possibly give them a reasonable name at some point
% t.Properties.VariableNames = 
t.y = y;
t.subject = repelem(1:length(input),cellfun(@(x)size(x.data,2),input))';
writetable(t,'cache/unfold_export.csv')
%%
tic
mres = fitlme(t,'y ~ 1 + Var1 + Var2 + Var3 + Var4 + Var5 + Var6 + Var7 + Var8 + Var9 + Var10 + Var11 + Var12 + Var13 + Var14 + Var15 + Var16 + Var17 + Var18 + Var19 + Var20 + Var21 + Var22 + Var23 + Var24 + (-1 + Var1| subject) + (-1 + Var2| subject) + (-1 + Var3| subject) + (-1 + Var4| subject) + (-1 + Var5| subject) + (-1 + Var6| subject) + (-1 + Var7| subject) + (-1 + Var8| subject) + (-1 + Var9| subject) + (-1 + Var1-1|subject) + (-1 + Var11| subject) + (-1 + Var12| subject) + (-1 + Var13| subject) + (-1 + Var14| subject) + (-1 + Var15| subject) + (-1 +Var16| subject) + (-1 + Var17| subject) + (-1 + Var18| subject) + (-1 + Var19| subject) + (-1 + Var2-1| subject) + (-1 + Var21| subject) + (-1 + Var22| subject) + (-1 + Var23| subject) + (-1 + Var24| subject)')
toc
%%
tic

EEG = um_mmfit(EEG,input,'channel',1,'covariance','Diagonal');
umresult = um_condense(EEG);

toc