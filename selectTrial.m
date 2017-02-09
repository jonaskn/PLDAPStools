function [pa] = selectTrial(pa, trialLevelMatrix, iTrial, trialStartParams)

    if nargin<3
        trialStartParams=false;
    end
    
    if trialStartParams
        [~,structNames]=getAllStructs(pa);
        isData=cellfun(@(x) textscan(x,'data%d'),structNames);
        isData=~cellfun(@isempty,isData);
        isAnalysis=cellfun(@(x) textscan(x,'analysis%d'),structNames);
        isAnalysis=~cellfun(@isempty,isAnalysis);
        trialLevelMatrix(isData|isAnalysis,:)=false;
    end
    pa.setLevels(trialLevelMatrix(:,iTrial));