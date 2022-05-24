% runs k-means cluster analysis (k = n pops) on simulation results

% clusters all individuals by population at the transition to sympatry and
% the end of the simulations

%% cluster analysis
d = cd; %SET TO DIRECTORY CONTAINING SIMULATION RESULTS
    %d should contain a folder called 'simulation_results'
    %'simulation results' should contain folders for each simulated world
    %each of these folders should contain .mat and .tar files for all iterations

%temp folder used to extract tars from simulation results
if isfolder(fullfile(d,'tmp'))==0
    mkdir(fullfile(d,'tmp'))
end

%individual clustering results saved here
if isfolder(fullfile(d,'cluster_data'))==0
    mkdir(fullfile(d,'cluster_data'))
end

%get simulation
sims = dir(fullfile(d,'simulation_results/*')); 
sims = {sims.name}';
sims = sims(ismember(sims,{'.','..'})==0);

w = warning;
warning('off')

%set up data table for results
datClust = table('Size',[0 13], ...
                 'VariableTypes',[{'string'} repmat({'double'},[1 12])], ...
                 'VariableNames',{'simulation','iteration','generation','population','modeClust','assignedClust','nCorrect','nIncorrect','propCorrect','overallNCorrect','overallNIncorrect','overallPropCorrect','flag'});%,'feat1','feat2','feat3','feat4','feat5','feat6','feat7','feat8','feat9','feat10','feat11','feat12','feat13','feat14','feat15'}); %GAH

%run cluster analysis on each simulation
for s = 1:length(sims)
    disp(sims{s})
    for i = 1:28
        %determine transition to sympatry
        if contains(sims{s},'10pAllo')
            transGen = 2000;
        elseif contains(sims{s},'50pAllo')
            transGen = 10000;
        else
            transGen = 20000;
        end

        endGens = unique([transGen 20000]);
        
        %deteremine if already processed (if process was killed)
        processed = 0;
        if exist(fullfile(d,'cluster_data',[sims{s} '_itt' num2str(i,'%02.0f') '_gen' num2str(endGens(1)) '.mat']),'file') && ...
           exist(fullfile(d,'cluster_data',[sims{s} '_itt' num2str(i,'%02.0f') '_gen' num2str(endGens(end)) '.mat']),'file')
            processed = 1;
        else
            %extract simulation results
            untar(fullfile(d,'simulation_results',sims{s},['iteration' num2str(i,'%02.0f') '_generations.tar']), ...
                  fullfile(d,'tmp'));
        end
        
        %run cluster analysis for each relevant generation (transition to sympatry, end)
        for g = 1:length(endGens)
            if processed
                load(fullfile(d,'cluster_data',[sims{s} '_itt' num2str(i,'%02.0f') '_gen' num2str(endGens(g)) '.mat']))
            else
                %load data for relevant generation
                load(fullfile(d,'tmp',['gen' num2str(endGens(g)) '.mat']));
                
                %format data
                npop = size(featVals,3);
                
                f = featVals(1:15,:,1)';
                for p = 2:npop
                    f = [f; featVals(1:15,:,p)']; %#ok<AGROW> 
                end
                
                %run k means
                [idx,C] = kmeans(f,npop);  
                idx = reshape(idx,[],npop);
            
                %for each population number, determine the cluster with the most matches
                assignedClust = zeros(npop,1);
                for p = 1:npop
                    [~,maxIdx] = max(sum(idx==p));
                    assignedClust(maxIdx) = p;
                end
    
                %fix situations where some populations are unassigned
                n = 0;
                flag = 0;
                while sum(assignedClust==0)>0 
                    n = n + 1;
    
                    %find (first) unassigned pop
                    missing = find(~ismember(1:npop,assignedClust),1);
    
                    %find cluster it wants
                    [~,disputedClust] = max(sum(idx==missing)); 
    
                    %assign disputed cluster to the higher of the two options
                    if sum(idx(:,disputedClust)==missing) > sum(idx(:,disputedClust)==assignedClust(disputedClust))
                        assignedClust(find(assignedClust==0,1)) = assignedClust(disputedClust); %assign pop that currently has disputed clust to the first zero 
                        assignedClust(disputedClust) = missing;
                    else
                        assignedClust(find(assignedClust==0,1)) = missing; %assign missing to first zero 
                    end
    
                    %catch infinite loops
                    if n>npop
                        assignedClust = 1:npop;
                        flag = 3; %these should be re-done!
                    end
                end
                if n==1
                    flag = 1; %noting it; should be fine
                elseif n>2
                    flag = 2; %these should be checked!
                end
    
                %calculate correctly classified based on modeClust (will overwrite later when necessary)
                correct = idx==repmat(assignedClust',size(idx,1),1); 

                %save
                save(fullfile(d,'cluster_data',[sims{s} '_itt' num2str(i,'%02.0f') '_gen' num2str(endGens(g)) '.mat']), ... 
                     'idx','C','assignedClust','correct','n','flag','-v7.3')
            end

            npop = size(idx,2);
            nind = size(idx,1);
            
            %add to file
            rows = height(datClust)+1:height(datClust)+npop;

            datClust.simulation(rows) = sims{s};
            datClust.iteration(rows) = i;
            datClust.generation(rows) = endGens(g);
            datClust.population(rows) = 1:npop;
            datClust.modeClust(rows) = mode(idx)'; %most commonly identified cluster for each population
            datClust.assignedClust(rows) = assignedClust; %assigned cluster for each population
            datClust.nCorrect(rows) = sum(correct); %number correctly classified for each population
            datClust.nIncorrect(rows) = nind-datClust.nCorrect(rows); %number of incorrectly classified for each population
            datClust.propCorrect(rows) = datClust.nCorrect(rows)./nind; %proportion correct for each population
            datClust.overallNCorrect(rows) = sum(correct(:)); %number of correctly classified across all populations
            datClust.overallNIncorrect(rows) = (nind*npop)-datClust.overallNCorrect(rows); %number incorrectly classified across all populations
            datClust.overallPropCorrect(rows) = datClust.overallNCorrect(rows)./(nind*npop); %proportion correct across all populations
            datClust.flag(rows) = flag;
            
            %clean up
            delete([d 'tmp/*.mat'])
            clear featVals npop f p idx C assignedClust maxIdx n flag missing disputedClust correct rows
        end
        clear transGen endGens processed
    end
end

%save results
save(fullfile(d,'datClust.mat'),'datClust','-v7.3')
writetable(datClust,fullfile(d,'data_clusters.csv'))

warning(w);
