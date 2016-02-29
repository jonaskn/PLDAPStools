%%  recreateParams

%load a PDS
load test.PDS -mat
%recreate a valid params class with all information
[pa, trialLevelMatrix] = recreateParams(PDS);

%select the parameters that where valid after a trial:
iTrial=1;
pa.setLevels(trialLevelMatrix(:,iTrial));