function  [E2PTBfit, E2PTB,PTB2E, maxreconstructionerror, filenameE] = synchronizeEphysClock(PDS, filenameE)

openephys=false;
plexon=false;

if ~isstruct(PDS)
    filename=PDS;
    if nargin<2
        filenameE=fileparts(filename);
    end
    load(PDS,'-mat')
end

[~, ~, ext] = fileparts(filenameE);
if(isempty(ext))
    baseDir=filenameE;
    
    if isfield(PDS.initialParametersMerged, 'plexon') && PDS.initialParametersMerged.plexon.spikeserver.use
        plexon=true;

        if isfield(PDS.initialParametersMerged.plexon, 'filename')
            filenameE=PDS.initialParametersMerged.plexon.filename;
        else
            filenameE=cellfun(@(x) x.plexon.filename, PDS.data(cellfun(@(x) isfield(x.plexon, 'filename'), PDS.data)), 'UniformOutput', false);
            if isempty(filenameE)
                candidates=dir([baseDir filesep '*.pl*']);
                
%                 exp='([a-z A-Z _ 0-9 \- /]+)([/])([a-z A-Z])(?<date>[0-9]+)_e([0-9])([a-z A-Z _ 0-9 \- /]+)_t(?<time>[0-9]+).';
                exp='([a-z A-Z])(?<date>[0-9]+)_e([0-9])([a-z A-Z _ 0-9 \- /]+)_t(?<time>[0-9]+).';
                [finds,good]=regexp({candidates.name}, exp, 'names');
                candidates=candidates(~~[good{:}]);
                candidateTimes=cellfun(@(x)  datenum([x.date x.time], 'yyyymmddhhMM'), finds(~~[good{:}]));
                [candidateTimes, cti]=sort( candidateTimes);
                candidates=candidates(cti);
                
                choice=find(candidateTimes<PDS.initialParametersMerged.session.initTime,1,'last');
                
                candidates=candidates(candidateTimes==candidateTimes(choice));
                %multiple files with the same, create preference: pl2 over
                %plx, tdddd.pl2 over tddd-.pl2
                if length(candidates)>1
                    %%remove those that are longer than the shortest one
                    l=cellfun(@length,{candidates.name});
                    candidates=candidates(l==min(l));
                    
                    if length(candidates)>1
                         ispl2=cellfun(@(x) strcmpi(x(end-2:end), 'pl2'),{candidates.name});
                         candidates=candidates(ispl2);
                    end
                    
                end
                
                if ~isempty(candidates)
                    filenameE=[baseDir filesep candidates().name];
                else
                    error('there should be at least one file here')
                end
                
% %                 datevec([a.date a.time], 'yyyymmddhhMM')
%                 PDS.initialParametersMerged.session.initTime
% 
% 
% 
%                 [cfile, cpath] = uigetfile('*.pl*', 'choose ephys file', baseDir);
%     %             baseDir=fileparts(filename);
%     %             a=dir([baseDir filesep '*.pl*']);
%     %                         %N=datenum({a.date});
%     %                         [~,ii]=max([a.datenum]);
%     %                         %a(ii).name
%     %                         
%     %                         filename=[strrep(baseDir,'\','/') a(ii).name];
%     %                         WaitSecs(0.5);
%     %                         a2=dir(filename);

            elseif length(unique(filenameE))==1
                filenameE=filenameE{1};
            else
                error('implement how to choose if multiple files are present (e.g. from sorting...)')
            end

        end

    end

    if isfield(PDS.initialParametersMerged, 'openephys') && PDS.initialParametersMerged.openephys.use
        openephys=true;
    %     filenameE=
        PDS.initialParametersMerged.openephys.status.recordingPath
        PDS.initialParametersMerged.openephys.status.experimentNumber
        %correct way would be to check .status.recording, make sure thare is
        %only one recording period and check the filenames
%         filenameE=[baseDir 

    end
else
    if strcmp(ext,'.pl2')||strcmp(ext,'.plx')
        plexon=true;
    end
end
%parse p.trial.plexon.filename
if plexon
    [E2PTBfit,E2PTB,PTB2E, maxreconstructionerror] = ephys.syncPlexonClock(PDS, filenameE);
elseif openephys
    [E2PTBfit, E2PTB,PTB2E, maxreconstructionerror] = ephys.syncOpenEphysClock(PDS, filenameE);
end