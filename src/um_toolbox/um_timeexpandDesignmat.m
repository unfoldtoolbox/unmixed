function [EEG] = um_timeexpandDesignmat(EEG,varargin)
[cfg, uftimeexpandCFG] = finputcheck(varargin,...
    {'placeholder','boolean','',0},'mode','ignore');

cfg.nranef = length(EEG.unmixed.uf_ranef);

% First timeexpand fixef
EEGtmp = EEG;
EEGtmp.unfold = EEGtmp.unmixed.uf_fixef;
EEGtmp = uf_timeexpandDesignmat(EEGtmp,uftimeexpandCFG{:});
EEG.unmixed.uf_fixef = EEGtmp.unfold;
clear EEGtmp;

for k = 1:length(EEG.unmixed.uf_ranef)
    
    EEGtmp = EEG;
    EEGtmp.unfold = EEGtmp.unmixed.uf_ranef{k};
    EEGtmp = uf_timeexpandDesignmat(EEGtmp,uftimeexpandCFG{:});
    
    groupingvar = EEGtmp.unfold.ranefgrouping;
    nTimeshifts = length(EEGtmp.unfold.times);
    [un,ia,ic]= unique([EEGtmp.event.(groupingvar)]);
    nGroupingLevels =length(un);
    
    
    groupingvarCol = EEGtmp.unfold.Xdc_terms2cols == find(ismember(EEGtmp.unfold.variablenames,groupingvar));

    groupXdc= EEGtmp.unfold.Xdc(:,groupingvarCol);
    
    sz = size(EEGtmp.unfold.Xdc(:,~groupingvarCol));
    
    splitZdc = {};
    Zdc_terms2cols = [];Zdc_level2cols = [];
    for col = unique(groupXdc(groupXdc(:)~=0))'

        
        splitZdc{col} = spalloc(sz(1),sz(2),1);
        splitZdc{col}(any(groupXdc==col,2),:) = EEGtmp.unfold.Xdc(any(groupXdc==col,2),~groupingvarCol);
        Zdc_terms2cols = [Zdc_terms2cols EEGtmp.unfold.Xdc_terms2cols(~groupingvarCol)];
        Zdc_level2cols =  [Zdc_level2cols repmat(col,1,sum(~groupingvarCol))];
    end
        
    % We generated a block diagonal matrix.
%     Zdc =  blkdiag(splitZdc{:});
    Zdc = cat(2,splitZdc{:});
    clear splitZ
    % Now we need to shuffle the columns so that the organization is the
    % following
    % [ time lag (tl)   1 ] [ time lag  (tl) 2  ]
    % [ all N group-levels] [ all N group-levels]
    % [ [X1 X2 ...] [X1..]] [ [X1 X2 ...] [X1..]]
    
    % instead of 
    % 
    % [ N = 1                        ] [ N == 2                       ]
    % [[X1     ] [X2     ] [X3      ]] [[X1     ] [X2     ] [X3      ]]
    % [[tl1 tl2] [tl1 tl2] [tl1 tl2]]] [[tl1 tl2] [tl1 tl2] [tl1 tl2]]] 
    
    % We generate a permutation matrix for this
    P = 1:size(Zdc,2);
    P = reshape(P,nTimeshifts,[],nGroupingLevels);
    
    P = permute(P,[2 3 1]);
    
    
    % P2 = 1:size(P,2);
    % P2 = reshape(P2,[],nSubjects);
    % P2 = P2';
    % P = P(:,P2);
    % P = P';
    Zdc= Zdc(:,P(:));
    
    EEGtmp.unfold.Zdc_terms2cols = Zdc_terms2cols(P(:));
    EEGtmp.unfold.Zdc_level2cols = full(Zdc_level2cols(P(:)));
    EEGtmp.unfold.Zdc = Zdc;
    EEG.unmixed.uf_ranef{k} = EEGtmp.unfold;

    
end

