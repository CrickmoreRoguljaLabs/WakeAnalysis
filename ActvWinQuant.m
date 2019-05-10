
%4/30/2019 Bryan Song
% Deleted inputs to be used with Wakelax

%%

% Prime the cell to store data

dataOutput = cell(1,n_genos);

%% Find Start Time 
%Subtract 1 because the input is 1 hour higher than the CT value. 
startActivity4 = timesStart - 1;
%Multiply input by 60 and divide by 5 
startActivity3 = (startActivity4 * 60);
startActivity2 = (startActivity3 / 5);
startActivity = startActivity2 + 1;

%% Find End time 
%%Repeat for end activity
%Subtract 1 because the input is 1 hour higher than the CT value. 12 bins is 1 hour. 
endActivity3 = timesEnd - 1;
%Multiple input by 60 and divide by 5 
endActivity2 = (endActivity3 *60);
endActivity = (endActivity2 / 5);

%%
%Prime excel sheet for output
indexOutputCell = cell(sum([master_data_struct.num_processed_flies])+1, n_genos);
immediateChangeOutputCell = cell(sum([master_data_struct.num_processed_flies])+1, n_genos);
percentChangeOutputCell = cell(sum([master_data_struct.num_processed_flies])+1, n_genos);
deltaChangeOutputCell = cell(sum([master_data_struct.num_processed_flies])+1, n_genos);  
prePostOutputCell = cell(sum([master_data_struct.num_processed_flies])+1, (n_genos*2));
prePost60OutputCell = cell(sum([master_data_struct.num_processed_flies])+1, (n_genos*2));
immediateChange60OutputCell = cell(sum([master_data_struct.num_processed_flies])+1, n_genos); 
prePostXOutputCell = cell(sum([master_data_struct.num_processed_flies])+1, (n_genos*2));

%%

% Go through the genotypes
for j = 1 : n_genos

    % Obtain the number of flies
    
    n_flies = master_data_struct(j).num_alive_flies;
    
    alive_fly_indices = master_data_struct(j).alive_fly_indices > 0;
    
    if n_flies > 0

        % Rearrange the data to be 288 x days x flies
        data4d = reshape(master_data_struct(j).data(:,alive_fly_indices)...
            ,[5,288,n_days,n_flies]);

        data3d = squeeze(sum(data4d, 1));

        %%

        %%Make Prior Day Index (for future use)

        %%Sum activity for Day 1
        %Day1pre = squeeze(sum(data3d(startActivity:endActivity-1, theDay2analyze-1, :),1));
        %%Sum activity for Day 2
        %Day1post = squeeze(sum(data3d(startActivity:endActivity, theDay2analyze, :),1));

        %Get rid of divide by 0 errors
        %prePlus01 = Day1pre  + 0.0000001;
        %PostPlus01 = Day1post + 0.0000001;
        %preMinusPost = PostPlus01 - prePlus01;
        %prePlusPost = PostPlus01 + prePlus01;
        %performanceIndex = (preMinusPost ./ prePlusPost);


        %%
        %Make Immediate Change Index
        %%Sum activity for Light pulse day, X minutes before  the stimulus
        if n_days > 1
            preLightPulse = squeeze(sum(data3d(startActivity-percentChangeWindow:startActivity-1, theDay2analyze, :),1));
            % Sum activity for Light pulse day, X minutes after the stimulus
            postLightPulse = squeeze(sum(data3d(startActivity:startActivity+percentChangeWindow3, theDay2analyze, :),1));
        else
            preLightPulse = rot90(squeeze(sum(data3d(startActivity-percentChangeWindow:startActivity-1, :),1)));
            % Sum activity for Light pulse day, X minutes after the stimulus
            postLightPulse = rot90(squeeze(sum(data3d(startActivity:startActivity+percentChangeWindow3, :),1)));
        end

        preLightPulsePlus01 = preLightPulse + 0.0000001;
        postLightPulsePlus01 = postLightPulse + 0.0000001;

        immediateChange3 = postLightPulsePlus01 - preLightPulsePlus01;
        immediateChange2 = postLightPulsePlus01 + preLightPulsePlus01;
        immediateChange = (immediateChange3 ./ immediateChange2);

        %%
        %Make 60 min Immediate Change Index
        %%Sum activity for Light pulse day, 60 minutes before  the stimulus
        if n_days > 1
            pre60LightPulse = squeeze(sum(data3d(startActivity-12:startActivity-1, theDay2analyze, :),1));
            % Sum activity for Light pulse day, X minutes after the stimulus
            post60LightPulse = squeeze(sum(data3d(startActivity:startActivity+11, theDay2analyze, :),1));
        else 
            pre60LightPulse = rot90(squeeze(sum(data3d(startActivity-12:startActivity-1, :),1)));
            % Sum activity for Light pulse day, X minutes after the stimulus
            post60LightPulse = rot90(squeeze(sum(data3d(startActivity:startActivity+11, :),1)));
        end

        pre60LightPulsePlus01 = pre60LightPulse + 0.0000001;
        post60LightPulsePlus01 = post60LightPulse + 0.0000001;

        immediate60Change3 = post60LightPulsePlus01 - pre60LightPulsePlus01;
        immediate60Change2 = post60LightPulsePlus01 + pre60LightPulsePlus01;
        immediate60Change = (immediate60Change3 ./ immediate60Change2);

        %%
        %Make % change
        if n_days > 1
        %%Sum activity for Light pulse day, X minutes before  the stimulus
            prePercLightPulse = squeeze(sum(data3d(startActivity-percentChangeWindow:startActivity-1, theDay2analyze, :),1));
            % Sum activity for Light pulse day, X minutes after the stimulus
            postPercLightPulse = squeeze(sum(data3d(startActivity:startActivity+percentChangeWindow3, theDay2analyze, :),1));
        else
            prePercLightPulse = rot90(squeeze(sum(data3d(startActivity-percentChangeWindow:startActivity-1, :),1)));
            % Sum activity for Light pulse day, X minutes after the stimulus
            postPercLightPulse = rot90(squeeze(sum(data3d(startActivity:startActivity+percentChangeWindow3, :),1)));            
        end
            
        %Do Not try to get rid of divide by zero errors!! causes infinitely large 
        %percentage increases 

        percentChange3 = postPercLightPulse - prePercLightPulse;
        percentChange2 = (percentChange3 ./ prePercLightPulse);
        percentChange = percentChange2 * 100;

        %%
        %Make immediate delta change
        immediate30Delta = postLightPulse - preLightPulse;

        %%

        dataOutput{1,1} = 'Genotype';
        %dataOutput{2,1} = 'pre light activity';
        %dataOutput{3,1} = 'post light activity';
        %dataOutput{4,1} = 'performance index';
        dataOutput{5,1} = 'Copy labels to prism';
        dataOutput{6,1} = ['% Change ' num2str(percentChangeWindow2) ' min before and after'];
        dataOutput{7,1} = 'Delta Change';  
        dataOutput{8,1} = ['Immediate change index ' num2str(percentChangeWindow2) ' min before and after'];
        dataOutput{9,1} = ['pre' num2str(percentChangeWindow2) ' min immediately before and after'];
        dataOutput{10,1} = ['post' num2str(percentChangeWindow2) ' min immediately before and after'];
        dataOutput{11,1} = ['Immediate change index 60 min before and after'];
        dataOutput{12,1} = ['pre 60 min immediately before and after'];
        dataOutput{13,1} = ['post 60 min immediately before and after'];
        %dataOutput{X,1} = 'time analyzed';

        dataOutput{1,j+1} = master_data_struct(j).genotype;
        %dataOutput{2,j+1} = Day1pre;
        %dataOutput{3,j+1} = Day1post;
        %dataOutput{4,j+1} = performanceIndex;
        dataOutput{5,j+1} = strcat(';pre', '-', master_data_struct(j).genotype, ';post','-',master_data_struct(j).genotype); 
        dataOutput{6,j+1} = percentChange;
        dataOutput{7,j+1} = immediate30Delta;
        dataOutput{8,j+1} = immediateChange;
        dataOutput{9,j+1} = preLightPulse;
        dataOutput{10,j+1} = postLightPulse;
        dataOutput{11,j+1} = immediate60Change;
        dataOutput{12,j+1} = pre60LightPulse;
        dataOutput{13,j+1} = post60LightPulse;
        %dataOutput{X,2} = times2show(timesStart,3);

        %Throw Index data from dataOutput into cell for export
        %indexOutputCell(1,j) = {master_data_struct(j).genotype;};
        %for i = 1 : length(performanceIndex) 
        %indexOutputCell(i+2,j) = num2cell(performanceIndex(i,1));
        %end

        %Throw pre post data from dataOutput into cell for export
        %prePostOutputCell(1,(2*j-1)) ={strcat('pre', '-', master_data_struct(j).genotype)}; 
        %prePostOutputCell(1,(2*j)) = {strcat('post','-',master_data_struct(j).genotype)}; 
        %for i = 1 : length(Day1post)
        %prePostOutputCell(i+2,(2*j-1)) = num2cell(Day1pre(i,1));
        %prePostOutputCell(i+2,(2*j)) = num2cell(Day1post(i,1));
        %end

        %Throw Immediate Change Index data from dataOutput into cell for export
        immediateChangeOutputCell(1,j) = {master_data_struct(j).genotype;};
        for i = 1 : length(immediateChange) 
        immediateChangeOutputCell(i+2,j) = num2cell(immediateChange(i,1));
        end

        %Throw % Change data from dataOutput into cell for export
        percentChangeOutputCell(1,j) = {master_data_struct(j).genotype;};
        for i = 1 : length(percentChange) 
        percentChangeOutputCell(i+2,j) = num2cell(percentChange(i,1));
        end

        %Throw immediate30Delta data from dataOutput into cell for export
        deltaChangeOutputCell(1,j) = {master_data_struct(j).genotype;};
        for i = 1 : length(immediate30Delta) 
        deltaChangeOutputCell(i+2,j) = num2cell(immediate30Delta(i,1));
        end

        %Throw 60 min immediate change index data from dataOutput into cell for export
        immediateChange60OutputCell(1,j) = {master_data_struct(j).genotype;};
        for i = 1 : length(percentChange) 
        immediateChange60OutputCell(i+2,j) = num2cell(immediate60Change(i,1));
        end

        %Throw pre post 60 min data from dataOutput into cell for export
        prePost60OutputCell(1,(2*j-1)) ={strcat('pre', '-', master_data_struct(j).genotype)}; 
        prePost60OutputCell(1,(2*j)) ={strcat('post', '-', master_data_struct(j).genotype)}; 
        for i = 1 : length(pre60LightPulse)
        prePost60OutputCell(i+2,(2*j-1)) = num2cell(pre60LightPulse(i,1));
        prePost60OutputCell(i+2,(2*j)) = num2cell(post60LightPulse(i,1));
        end

        %Throw pre post X min data from dataOutput into cell for export
        prePostXOutputCell(1,(2*j-1)) ={strcat('pre', '-', master_data_struct(j).genotype)}; 
        prePostXOutputCell(1,(2*j)) = {strcat('post','-',master_data_struct(j).genotype)}; 
        for i = 1 : length(preLightPulse)
        prePostXOutputCell(i+2,(2*j-1)) = num2cell(preLightPulse(i,1));
        prePostXOutputCell(i+2,(2*j)) = num2cell(postLightPulse(i,1));

        end
    end
end

%Export cells to csv
cell2csv(fullfile(export_path,[filename_master(1:end-5),'_ULOSR.csv']),indexOutputCell);

%cell2csv(fullfile(export_path,[filename_master(1:end-5),'_Day1Day2prepostULOSR.csv']),prePostOutputCell);

cell2csv(fullfile(export_path,[filename_master(1:end-5),'_immediateIndexULOSR' num2str(percentChangeWindow2) 'minWindow.csv']),immediateChangeOutputCell);

%cell2csv(fullfile(export_path,[filename_master(1:end-5),'_percentChangeULOSR' num2str(percentChangeWindow2) 'minWindow.csv']),percentChangeOutputCell);

cell2csv(fullfile(export_path,[filename_master(1:end-5),'_immediate30delta' num2str(percentChangeWindow2) 'minWindow.csv']),deltaChangeOutputCell);

cell2csv(fullfile(export_path,[filename_master(1:end-5),'_immediatePrePostULOSR 60 min Window.csv']),prePost60OutputCell);

cell2csv(fullfile(export_path,[filename_master(1:end-5),'_immediateIndexULOSR 60 min Window.csv']),immediateChange60OutputCell);

cell2csv(fullfile(export_path,[filename_master(1:end-5),'_immediatePrePostULOSR' num2str(percentChangeWindow2) 'minWindow.csv']),prePostXOutputCell);


% Save the work space
save(fullfile(export_path,[filename_master(1:end-5),'_ULOSRworkspace.mat']));

