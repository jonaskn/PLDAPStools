%% l = lazyload.lazyload(mat73filename)
%% loads information about a file saved in v7.3 mode without loading the
% data. 
%
%Will only the requested data from deep within a structure
% E.g. if you have data is PDS.data (a cell array of size 5) in a file 
% called test.mat and you'd like to have the data PDS.data{5}.timing.fliptimes(1,:) 
% 
% l=lazyload('test,mat')
% fliptimes=l.PDS.data{5}.timing.flipTimes(1,:);
%
%Lazyload allows you to get data from multiple cell fields at the ame time
% fliptimes=l.PDS.data{:}.timing.flipTimes(1,:);
% fliptimes=l.PDS.data{1:5}.timing.flipTimes(1,:);
% would return a cell array with the fliptime of each cell array
%
%Note that the first time a cell array is called is slow and slower the
%larger the cell array is. But next calls will be fast again.
% l=lazyload('test,mat')
% fliptimes=l.PDS.data{:}.timing.flipTimes(1,:); %a bit slow
% trstart=l.PDS.data{:}.trstart; %fast
%
%% written by Jonas Kn?ll 2015
classdef lazyload % <handle
    properties (Hidden = true)
        filename
        root=[];
        rootFields={};
        
        uberInfo
        
        referenceLoadMethod=2;%1:H5R.getname, 2:H5R.create, 3:guess

%         laziness = 2; % 1 (default): load full tree
%                       % 2 (don't dereference unless it has been requested)
%                       % 
%                       % 
                      
%TODO: find a way to better dereference cell arrays. On a nice and clean
%cell array, references are monotonically increasing.
%In that case we could just get the first elements reference and calculate
%the other ones from that
%but it seems that during savind some history of how cells got organized
%may stay preserved
                      
%TODO: test if we can combine multiple files into one lazyload

%TODO:
%even lazier: don't use the uberInfo, but get the initial list from
% tic;
% fid = H5F.open('testfile.mat');
% group_id = H5G.open(fid,'/');
% info = H5G.get_info(group_id);
% idx_type = 'H5_INDEX_NAME';
% order = 'H5_ITER_INC';
% lapl_id = 'H5P_DEFAULT';
% name=cell(1,info.nlinks);
% for iLink=1:info.nlinks
%     name{iLink} = H5L.get_name_by_idx(fid,'/',idx_type,order,iLink-1,lapl_id);
% end
%     
% H5G.close(group_id);
% H5F.close(fid);
% toc*1000
    end
    
	methods (Hidden = true)

        function l=lazyload(filename,location,info,uberInfo,root,referenceLoadMethod)           
            if nargin <2
                location='/';
            end
            
            if nargin<4
               l.uberInfo = lazyload.infos(h5info(filename));
            else
               l.uberInfo = uberInfo;
            end
                      
            if nargin<3 || isempty(info)
                if strcmp(location,'/')
                    info=getData(l.uberInfo);
                else
                    info=getSubInfo(l,uberInfo,location);
                end
            end

            if nargin >4
                l.root=root;
            end
            
            if nargin >5
                l.referenceLoadMethod=referenceLoadMethod;
            end
            
            if l.referenceLoadMethod==2
            %%later we will do this on initialization
               uI = getData(l.uberInfo);
               if ~isfield(uI,'refNames')
                   refGroup=uI.Groups(strcmp('/#refs#',{uI.Groups.Name}));
                   names=[];
                   if ~isempty(refGroup.Datasets)
                       names={refGroup.Datasets.Name};
                       names=cellfun(@(x) ['/#refs#/' x], names,'UniformOutput',false);
                   end
                   if ~isempty(refGroup.Groups)
                        names=[names {refGroup.Groups.Name}];
                   end
                    
                   fid=H5F.open(filename);
                   a=zeros(8,length(names));
                   for iName=1:length(names)
                       a(:,iName)=H5R.create(fid,names{iName},'H5R_OBJECT',-1);
                   end
                   H5F.close(fid);
                   l.uberInfo.refNames=names;
                   l.uberInfo.refData=double(a);
               end
            end
            
 
            if isfield(info,'Groups') && ~isempty(info.Groups)
                groupNames={info.Groups.Name};
                for iName=1:length(groupNames)
                    fin=(strfind(groupNames{iName}, '/'));fin=fin(end);
                    groupNames{iName}=groupNames{iName}(fin+1:end);
%                     splName=strsplit(groupNames{iName},'/');
%                     groupNames{iName}=splName{end};
                end
                 groupNames(~cellfun(@isempty,strfind(groupNames,'#')))=[];
            else
                groupNames={};
            end

            if isfield(info,'Datasets') && ~isempty(info.Datasets)
                datasetNames= {info.Datasets.Name};
            else
                datasetNames={};
            end
            
            if isfield(info,'Filters') && ~isempty(info.Filters)
                if ~isempty(datasetNames) || ~isempty(groupNames)
                    error('both cell root and struct/data???');
                end
                l.rootFields=info.Filters;
            else
                l.rootFields=cell2struct(cell(1,length(datasetNames)+length(groupNames)),sort([datasetNames groupNames]),2);
            end
            
            l.filename=filename;
        end
        
        function disp(l)
            builtin('disp',l);
            if length(l)==1
                builtin('disp',l.rootFields);
            else
                all_eq=true;
                for iLZ=1:length(l)-1
                   all_eq=all(strcmp(fieldnames(l(iLZ).rootFields),fieldnames(l(iLZ+1).rootFields)));
                   if ~all_eq
                       break;
                   end
                end
                if all_eq
                    builtin('disp',l(1).rootFields);
                end
            end
        end
        
          function is = isfield(l,s) 
             is = any(strcmp(fieldnames(l.rootFields),s));
        end
        
        function names = fieldnames(l) 
             names = sort(fieldnames(l.rootFields));
        end
        
        function names = properties(l)
            names = fieldnames(l);
        end
        
        %wrapper to handle multiple lazyloads
        function varargout = subsref(l,S)
            if length(l)>1 && strcmp(S(1).type,'{}')
                sz=size(l);
                testArray=ones(sz);
                testArray(:)=1:prod(sz);
                tmpS=S(1);
                tmpS.type='()';
                theseInds=builtin('subsref',testArray,tmpS);
                S(1)=[];
            else
                theseInds=1:length(l);
            end
            out = cell(1,length(l));
            for iLz = theseInds;%1:length(l)
                out{iLz} = subsrefSingle(l(iLz),S);
            end
                     
            if nargout==1 && length(out)>1
                varargout={[out{:}]};
            else
                varargout= out;
            end
        end
        
        function varargout = subsrefSingle(l,S)
            
            S=[l.root S];
            
            switch S(1).type
                case '.'
                numS=length(S);
                if strcmp(S(end).type,'()')
                    numS=numS-1;
                end
                    uI = getData(l.uberInfo);
                    locations={['/' S(1).subs]};
                    [infos{1}, subs{1}]=getSubInfo(l,uI, locations{1});  
                      roots={S(1)};

                    iSub=1;
                    while iSub<numS
                      iSub=iSub+1;
                        thisS=S(iSub);
                        nextLocations={};
                        nextInfos={};
                        nextSubs={};
                        nextRoot={};

                        switch thisS.type
                            case '{}' % index into cell array
                                for iLocation=1:length(locations)
                                    %ok, first we need to now the
                                    %dimensions
                                    sz=infos{iLocation}.Dataspace.Size;

                                    %I don't want to handle all 'end'
                                    %scenarios, will use matlabs
                                    %internal code for it:
                                    testArray=ones(sz);
                                    testArray(:)=1:prod(sz);
                                    tmpS=thisS;
                                    tmpS.type='()';
                                    theseInds=builtin('subsref',testArray,tmpS);

                                    %as far as I see a cell will always
                                    %use references only
                                    if ~strcmp(infos{iLocation}.Datatype.Class,'H5T_REFERENCE')
                                        error('non reference cell: not implemented yet...');
                                    end
                                    
                                    if isempty(infos{iLocation}.Filters)
                                        if l.referenceLoadMethod==2
                                            infos{iLocation}.Filters = getReferenceFast(l,l.filename,locations{iLocation},1:prod(sz));
                                            tmpS=subs{iLocation};
                                            tmpS(end+1).type='.';
                                            tmpS(end).subs='Filters';
                                            l.uberInfo = subsasgn(l.uberInfo,tmpS,infos{iLocation}.Filters);
                                        elseif l.referenceLoadMethod==3
                                            infos{iLocation}.Filters= guessReference(l,l.filename,locations{iLocation},1:prod(sz));
                                            tmpS=subs{iLocation};
                                            tmpS(end+1).type='.';
                                            tmpS(end).subs='Filters';
                                            l.uberInfo = subsasgn(l.uberInfo,tmpS,infos{iLocation}.Filters);
                                        else
                                            infos{iLocation}.Filters=cell(1,prod(sz));
                                        end
                                        
                                    end
                                    
                                    if l.referenceLoadMethod==1
                                        preLoaded=~cellfun(@isempty,infos{iLocation}.Filters);
                                        loadRefs=intersect(find(~preLoaded),theseInds);
                                        infos{iLocation}.Filters(loadRefs) = (getReference(l,l.filename,locations{iLocation},loadRefs));
                                        
                                        if ~isempty(loadRefs)
                                            tmpS=subs{iLocation};
                                            tmpS(end+1).type='.';
                                            tmpS(end).subs='Filters';
                                            l.uberInfo = subsasgn(l.uberInfo,tmpS,infos{iLocation}.Filters);
                                        end
                                    end
%                                     if ~isempty(loadRefs)
%                                         tic;
%                                         ref = guessReference(l,infos{iLocation}.Filters{loadRefs(1)},loadRefs(1),loadRefs);
%                                         toc
%                                         if ~all(cellfun(@(x,y) strcmp(x,y),ref, infos{iLocation}.Filters(loadRefs)' ))
%                                            warning('found not monotonically increasing references...') 
%                                         end
%                                         
%                                         tic;
%                                         ref = guessOrGetReference(l,l.filename,locations{iLocation},loadRefs);
%                                         toc
%                                         if ~all(cellfun(@(x,y) strcmp(x,y),ref, infos{iLocation}.Filters(loadRefs) ))
%                                            warning('guessOrGetReference failed.....') 
%                                         end
%                                     end


                                    nextLocations(end+1:end+length(theseInds)) = infos{iLocation}.Filters(theseInds);
                                    for iNextLocation=1:length(nextLocations)
                                        [nextInfos{end+1}, nextSubs{end+1}]=getSubInfo(l,uI,nextLocations{iLocation});
                                        
                                        nextRoot{end+1}=[roots{iLocation} thisS];
                                        nextRoot{end}(end).subs={theseInds(iNextLocation)};
                                    end
                                end
                            case '()' % index into aray
                                error('structs, classes not implemented yet')
                            case '.'  % get from a struct
                                %as far as I can see this is
                                %simple, just add the .subs to
                                %the location
                                nextLocations=strcat(locations,['/' thisS.subs]);
                                for iNextLocation=1:length(nextLocations)
%                                     [nextInfos{end+1}, nextSubs{end+1}]=getSubInfo(l,infos{iNextLocation},thisS.subs);
                                    [nextInfos{end+1}, nextSubs{end+1}]=getSubInfo(l,infos{iNextLocation},thisS.subs, subs{iNextLocation});
                                    nextRoot{end+1}=[roots{iNextLocation} thisS];
%                                     nextRoot{end+1}=[roots{:} thisS];
                                end
                        end

                        locations=nextLocations;
                        infos=nextInfos;
                        subs = nextSubs;
                        roots = nextRoot;
                        
                        %also condider resolving a final implicit
                        %{} in case the requested data ends on a cell
                        %i.e. treat a.b.c as a.b.c{:} if c is a cell
                        %If the last is a cell and the data inside is real
                        %data, we prob. wnant to add a {size}
                        if iSub==numS
                            %%assuming homogenous data for now that why we
                            %%use index 1
                           if isfield(infos{1},'Datatype') && strcmp(infos{1}.Datatype.Class,'H5T_REFERENCE') %if this is a cell array
                               sz=infos{1}.Dataspace.Size; 
%                                if ~isfield(infos{1}, 'Filters') || isempty(infos{1}.Filters) || all(cellfun(@isempty,infos{1}.Filters))
%                                    infos{1}.Filters=cell(1,prod(sz));
%                                    %get one
%                                    infos{1}.Filters(1) = (getReference(l,l.filename,locations{1},1));
% 
%                                    tmpS=subs{1};
%                                    tmpS(end+1).type='.';
%                                    tmpS(end).subs='Filters';
%                                    l.uberInfo = subsasgn(l.uberInfo,tmpS,infos{1}.Filters);
%                                end
%                                preLoaded=find(~cellfun(@isempty,infos{1}.Filters),1,'first');
                               %check if this leads to a data field
%                                testInfo=getSubInfo(l,uI,infos{1}.Filters{preLoaded});%getSubInfo(l,infos{1},thisS.subs, infos{1}.Filters(preLoaded));
                               %it's a data leaf
%                                if isfield(testInfo,'Datatype') && ~strcmp(testInfo.Datatype.Class,'H5T_REFERENCE')
                                   S(end+1).type='{}';
                                   S(end).subs=repmat({':'}, 1,length(sz));
                                   numS=numS+1;
%                                end
                           end
                           
                           
                           
                           
                        end
                    end %while iSub
                    
                    %ok, all but the last () checked. We could try
                    %to figure out how to reduce the amount read by
                    %only reading the requested data, but the input
                    %to h5read needs some carefull thinking for that
                    %for now I simply call subfsref on the read full
                    %data
                    
                    varargout=cell(1,length(locations));
                    isLZ = false(1,length(locations));
                    for iLocation=1:length(locations)
                        if isfield(infos{iLocation},'Datatype') && ~strcmp(infos{iLocation}.Datatype.Class,'H5T_REFERENCE')
                            varargout{iLocation} = h5read(l.filename,locations{iLocation});
                            mc_info=find(strcmp({infos{iLocation}.Attributes.Name},'MATLAB_class'));
                            if ~isempty(mc_info) && strcmp(infos{iLocation}.Attributes(mc_info).Value, 'char')
                                varargout{iLocation}=char(varargout{iLocation});
                            end
                        else                        
                            varargout{iLocation}=lazyload.lazyload(l.filename,locations{iLocation},infos{iLocation},l.uberInfo,roots{iLocation});
                            isLZ(iLocation) = true;
                        end
                    end
                    if numS~=length(S)
                        for iLocation=1:length(locations)
                            varargout{iLocation} = builtin('subsref',varargout{iLocation},S(end));
                        end
                    end
                    
                    if all(isLZ)
                       varargout={[varargout{:}]}; 
                    end
                    
                    if nargout==1 && length(varargout)>1
                        varargout= {varargout};
                    end
            end
        end
        
        function ref = getReferenceFast(l,filename,location,inds)
           plist = 'H5P_DEFAULT';
           space = 'H5S_ALL';
           fid=H5F.open(filename);
           
           %%later we will do this on initialization
           uI = getData(l.uberInfo);
           names=uI.refNames;
           a=uI.refData;
%            refGroup=uI.Groups(strcmp('/#refs#',{uI.Groups.Name}));
%            names={refGroup.Datasets.Name};
%            names=cellfun(@(x) ['/#refs#/' x], names,'UniformOutput',false);
%            names=[names {refGroup.Groups.Name}];
%            for iName=1:length(names)
%                a(:,iName)=H5R.create(fid,names{iName},'H5R_OBJECT',-1);
%            end
% %            
%            
%            for iName=1:length(names2)
%                b(:,iName)=H5R.create(fid,names2{iName},'H5R_OBJECT',-1);
%            end
%            
%            
%            a=[a b];
           
           multiplyers=256.^(0:7);
           
           did=H5D.open(fid,location);
%            ref=cell(1,length(inds));%cell(info.Dataspace.Size);
           ref_ind=ones(1,length(inds));%cell(info.Dataspace.Size);
           refdata = H5D.read(did,'H5T_STD_REF_OBJ',space,space,plist);
           for iRef=1:length(inds)
%                ref{iRef}=H5R.get_name(did,'H5R_OBJECT',refdata(:,inds(iRef)));
                ref_ind(iRef)=find( multiplyers*double(refdata(:,inds(iRef)))==multiplyers*double(a) );
           end
           
           ref=names(ref_ind);
           H5D.close(did);
           H5F.close(fid);
        end
        
        function ref = getReference(l,filename,location,inds)
           plist = 'H5P_DEFAULT';
           space = 'H5S_ALL';
           fid=H5F.open(filename);
           did=H5D.open(fid,location);
           ref=cell(1,length(inds));%cell(info.Dataspace.Size);
           refdata = H5D.read(did,'H5T_STD_REF_OBJ',space,space,plist);
           for iRef=1:length(inds)
               ref{iRef}=H5R.get_name(did,'H5R_OBJECT',refdata(:,inds(iRef)));
           end
           H5D.close(did);
           H5F.close(fid);
        end
        
        function ref = guessReference(l,reference,reference_ind,inds)
                %first get the reference number
                baseDictchar = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890';
%                 reference='cd';
                ref_num=1;
                reference=reference(9:end);
                for iN=1:length(reference)
                    ref_num=ref_num+((find(baseDictchar==reference(iN))-1)*62^(iN-1));
                end
                
                %now calculate references relative to that
                get_nums_10=ref_num-reference_ind+inds;
                %%go from number to 62 based char
                %'a' is zero, '0' is 61
                nNumbers=max(ceil(log(get_nums_10)/log(62)));
                get_nums_62=zeros(length(get_nums_10),nNumbers);
                get_nums_10=reshape(get_nums_10, length(get_nums_10),1);
                for iN=nNumbers:-1:1
                    get_nums_62(:,iN)=floor((get_nums_10-1)/62^(iN-1));
                    get_nums_10=get_nums_10-get_nums_62(:,iN)*62^(iN-1);
                end
                get_nums_62=get_nums_62+1;
                ref=reshape(baseDictchar(get_nums_62),size(get_nums_62));
                ref=num2cell(ref,2);
                %remove trailing zero (a)
                ref(all(get_nums_62(:,2:end)==1,2))=cellfun(@(x) x(1), ref(all(get_nums_62(:,2:end)==1,2)), 'UniformOutput', false);
                ref = cellfun(@(x) ['/#refs#/' x],ref, 'UniformOutput',false)';
        end

        function [info, S]=getSubInfo(r,uberInfo,location,S)
           if nargin < 4
            S=[];
           end
           start=strfind(location,'/');
           if isempty(start)
               split=location;
           else
               if start(1)==1
                    if length(start)==1
                        start(2)=length(location)+1;
                    end
                  start=start(2);
                  split=location(2:start-1);
               else
                   start=start(1);
                   split=location(1:start-1);
               end
           end
           if isempty(split)
               info=uberInfo;
               return
           end
           
           currentLocation=uberInfo.Name;
           if strcmp(currentLocation,'/')
               currentLocation='';
           end
           
           if isempty(uberInfo.Groups)
               groups=[];
           else
               groups=strcmp({uberInfo.Groups.Name},[currentLocation '/' split]);
           end
           if any(groups)
              S(end+1).type='.';
              S(end).subs=('Groups');
              S(end+1).type='()';
              S(end).subs={find(groups)};
              [info, S]=getSubInfo(r,uberInfo.Groups(groups),location(start+1:end),S);
           else
               datasets=strcmp({uberInfo.Datasets.Name},split);
               if any(datasets)
                   S(end+1).type='.';
                   S(end).subs=('Datasets');
                   S(end+1).type='()';
                   S(end).subs={find(datasets)};
                   [info, S]=getSubInfo(r,uberInfo.Datasets(datasets),location(start+1:end),S);
               else
                   error('location not found');
               end
           end
           
        end
    
	end
end