
%4/30/2019 Bryan Song
%Deleted inputs to be used with Wakelax


%%
%Prime color pallette for heatmap raster plots

%colorTable = flipud(cbrewer('seq','YlOrRd',100));
colorTable2 = cbrewer('seq','YlGnBu',100);
colorTable(1:70,:) = colorTable2(31:100,:);

%Make the first row white
colorTable(1,1:3) = ones;


%%


% Prime the cell to store data
dataOutput = cell(1,n_genos);

%Prime a variable to find maximum activity for y axis scale
maxAct = 0;

%%

%%Get 1 minute index values for the start input

%Subtract 1 because the input is 1 hour higher than the CT value. 60 bins is 1 hour. 
startActivity1minBins3 = baselineStart - 1;
%Multiply input by 60 and divide by 5 
startActivity1minBins2 = (startActivity1minBins3 *60);
startActivity1minBins = startActivity1minBins2 + 1;

startActivityVarBins2 = (startActivity1minBins2 / binLength);
startActivityVarBins = startActivityVarBins2 + 1;

%%

%%Repeat for end activity

%Subtract 1 because the input is 1 hour higher than the CT value. 12 bins is 1 hour. 
endActivity1minBins2 = timesEnd - 1;
%Multiple input by 60 and divide by 5 
endActivity1minBins = (endActivity1minBins2 *60);
endActivityVarBins = (endActivity1minBins / binLength);

%%
%Prime excel sheet for output
dataOutputCell = cell(1442, (n_genos*3));
dataOutputCellVar = cell(numOfBins+1, n_genos*3);

%%

% Go through the genotypes

for j = 1 : n_genos

    % Obtain the number of flies
    
    n_flies = master_data_struct(j).num_alive_flies;
    
    alive_fly_indices = master_data_struct(j).alive_fly_indices > 0;
    
    if n_flies > 0
    
        % Rearrange the data to be 1440 x days x flies (for 1 min bins)
        data3d = reshape(master_data_struct(j).data(:,alive_fly_indices)...
           ,[1440,n_days,n_flies]);
    
        %Make an emptry matrix to populate with activity. Arbitraily use the
        %first fly of the genotype to determine the length of activity string.
        %If Bug Replace data3d with data4d
        theVoid = zeros(length(data3d(startActivity1minBins:endActivity1minBins, theDay2analyze, 1)),n_flies);
    
        % Rearrange the data to be Number of vins x days x flies (for a user inputted bin)
        data4d = reshape(master_data_struct(j).data(:,alive_fly_indices)...
            ,[binLength,numOfBins,n_days,n_flies]);

        %Sum across the 5 minute bins
        sumData4d = sum(data4d,1);
        %Get rid of the extra dimension
        squeezeData4d = squeeze(sumData4d);
    
        %Make an emptry matrix to populate with activity. Arbitrarily the use the
        %first fly of the genotype to determine the length of activity string.
        theVoidVar = zeros(length(squeezeData4d(startActivityVarBins:endActivityVarBins, theDay2analyze, 1)),n_flies);
        
    else
        
        theVoid = zeros(length(data3d(startActivity1minBins:endActivity1minBins, theDay2analyze, 1)),n_flies);
        theVoidVar = zeros(length(squeezeData4d(startActivityVarBins:endActivityVarBins, theDay2analyze, 1)),n_flies);
        
    end
    
        for k = 1 : n_flies

            %Out activity values within defined window 
            activityString = data3d(startActivity1minBins:endActivity1minBins, theDay2analyze, k);
        
            if n_days > 1
              activityStringVar = squeezeData4d(startActivityVarBins:endActivityVarBins, theDay2analyze, k);
            else
                activityStringVar = squeezeData4d(startActivityVarBins:endActivityVarBins, k);
            end

            %Fill the empty matrix with each fly's activity string
            theVoid (:,k) = activityString;
            theVoidVar (:,k) = activityStringVar;
        
        end
    
 
        %Make rasterplots of individual flies
        %Make binary raster plots
    
        figure(j)
        %clf;
        tVec = 1:length(theVoid);
        %binarize data
        %binarytheVoid = theVoid > 0.1;
        binarytheVoid = logical(theVoid);
        %transpose the binary activity matrix such that time is on the x axis
        transbitheVoid = binarytheVoid';
        %hold all
        for flyCount=1:n_flies;
            hold on
            tickPos = tVec(transbitheVoid(flyCount, :));
            for spikeCount = 1:length(tickPos)
                plot([tickPos(spikeCount) tickPos(spikeCount)],...
                    [flyCount-0.4 flyCount+0.4], 'k');
            end
        end
        xlim([0 size(theVoid, 1)+1]);
        ylim([0 size(theVoid, 2)+1]);
        xlabel('Time (minutes)')
        ylabel('Individual Flies')
        hold off
    
        saveas(gcf, fullfile(export_path,[filename_master(1:end-5),'_',(master_data_struct(j).genotype),'_ULOSRaster.pdf']));
        close gcf
    
        figure(j)
        %Make heatmap raster plots
        %theVoid' is transposed theVoid
        %ULOSRheatMap = HeatMap(theVoid','Symmetric', false,'Colormap',colorTable);
        colormap(colorTable);
        clim = [0, 30];
        imagesc(theVoid',clim);
        colorbar;
        xlabel('time (mins)')
        ylabel('fly number')
        saveas(gcf, fullfile(export_path,[filename_master(1:end-5),'_',(master_data_struct(j).genotype),'_ULOSRheatPlot.fig']));

        close gcf
    
        %%
    
    
        %Estimate each genotypes' average. 
        averagedString = mean((theVoid),2);
        averagedStringVar = mean((theVoidVar),2);
    
        %Estimate each genotypes' std. Used x,1,2 instead of x,0,2 ie n vs n-1
        %because that's how the other programs estimate variability.
        stdString = std(theVoid ,0,2);
        semString = (stdString/ sqrt(n_flies));
        stdStringVar = std(theVoidVar ,0,2);
        semStringVar = (stdStringVar/ sqrt(n_flies));
    
        %Set up a test for largest y value in the array. Replace if larger than
        %existing one
        %maxActContender = max(averagedStringVar);
        %if maxActContender > maxAct;
        %    maxAct = maxActContender;
        %else
        %    maxAct = maxAct;
        %end
   
        %Generate a column of n flies for prism output
        nFliesPrism=zeros(length(activityString),1);
        nFliesPrism(:)=n_flies;
    
        nFliesPrismVar=zeros(length(activityStringVar),1);
        nFliesPrismVar(:)=n_flies;
    
        dataOutput{1,1} = 'Genotype';
        dataOutput{2,1} = 'Individual Activity Output 1 min Bins';
        dataOutput{3,1} = 'Averaged Activity Output 1 min Bins';
        dataOutput{4,1} = 'Individual Activity Output X min Bins';
        dataOutput{5,1} = 'Averaged Activity Output X min Bins';
        dataOutput{6,1} = 'First time point';
        dataOutput{7,1} = 'Last time point';
    
        dataOutput{1,j+1} = master_data_struct(j).genotype;
        dataOutput{2,j+1} = theVoid;
        dataOutput{3,j+1} = [averagedString, semString];
        dataOutput{4,j+1} = theVoidVar;
        dataOutput{5,j+1} = [averagedStringVar, semStringVar];
        %dataOutput{6,2} = times2show(timesStart,2);
        %dataOutput{6,3} = times2show(timesStart,3);
        %dataOutput{7,2} = times2show(timesEnd,2);
        %dataOutput{7,3} = times2show(timesEnd,3);

        %Throw data from dataOutput into cell for export
        dataOutputCell(1,(3*j-2)) = {strcat('Avg', '-', master_data_struct(j).genotype)};
        dataOutputCell(1,(3*j-1)) = {strcat('Sem', '-', master_data_struct(j).genotype)};
        dataOutputCell(1,(3*j)) = {'n flies'};

        dataOutputCell(3:length(averagedString)+2,3*j-2) = num2cell(averagedString);
        dataOutputCell(3:length(averagedString)+2,3*j-1) = num2cell(semString);
        dataOutputCell(3:length(averagedString)+2,3*j) = num2cell(nFliesPrism);

        %Throw data from dataOutput into cell for export
        dataOutputCellVar(1,(3*j-2)) = {strcat('Avg', '-', master_data_struct(j).genotype)};
        dataOutputCellVar(1,(3*j-1)) = {strcat('Sem', '-', master_data_struct(j).genotype)};
        dataOutputCellVar(1,(3*j)) = {'n flies'};

        dataOutputCellVar(3:length(averagedStringVar)+2,3*j-2) = num2cell(averagedStringVar);
        dataOutputCellVar(3:length(averagedStringVar)+2,3*j-1) = num2cell(semStringVar);
        dataOutputCellVar(3:length(averagedStringVar)+2,3*j) = num2cell(nFliesPrismVar);
    
    
     
end

%Export cells to csv
cell2csv(fullfile(export_path,[filename_master(1:end-5),'_ULOSRactivity1min.csv']),dataOutputCell);
cell2csv(fullfile(export_path,[filename_master(1:end-5),'_ULOSRactivity' num2str(binLength) 'mins.csv']),dataOutputCellVar);

%Set up a dummy array to make x axis
xAxisMat = zeros (1,length(averagedStringVar));
ticker = 0;


for m = 1 : length(averagedStringVar);
    xAxisMat(1,m) = ticker ;
    ticker = ticker + binLength;
end
    
 
% Construct a vector of rainbow groups corresponding to the genotypes
rainbowgroups_vector = cell2mat({master_data_struct.rainbowgroup})';

% Determine the unique rainbow groups (ignoring the NaNs) and their count
rainbowgroups_unique = unique(rainbowgroups_vector(rainbowgroups_vector>-99999),'stable');
rainbowgroups_unique = rainbowgroups_unique(rainbowgroups_unique~=0); %Group 0 reserved for universal controls
rainbowgroups_n = length(rainbowgroups_unique);

for j = 1:rainbowgroups_n
    
    % Find how many and which genotypes are of the current rainbow group
    geno_indices_of_the_current_rainbowgroup = [find(rainbowgroups_vector == rainbowgroups_unique(j));find(rainbowgroups_vector == 0)]; % Plot the current group and Group 0
    n_geno_of_the_current_rainbowgroup = length(geno_indices_of_the_current_rainbowgroup);

    % Prime the rainbow data matrix
    rainbow_mat = zeros(length(averagedStringVar),n_geno_of_the_current_rainbowgroup);
    rainbow_mat_sem = zeros(length(averagedStringVar),n_geno_of_the_current_rainbowgroup);

    % Prime the output rainbow cells
    rainbow_cell = cell(length(averagedStringVar)*2,n_geno_of_the_current_rainbowgroup);

    for i = 1:n_geno_of_the_current_rainbowgroup
        % Variables requested by people
        current_rainbow_geno = geno_indices_of_the_current_rainbowgroup(i);
        current_rainbow_alive_flies = master_data_struct(current_rainbow_geno).alive_fly_indices>0;

        % Construct rainbow data cell
        rainbow_cell{1,i} = genos{current_rainbow_geno};
        rainbow_cell{2,i} = master_data_struct(current_rainbow_geno).num_alive_flies;
        %rainbow_cell(3:50,i) = num2cell(temp_average_sleep_per_30_min);

        temp5minBins = cell2mat(dataOutput(5,(geno_indices_of_the_current_rainbowgroup(i))+1)); 
        longTemp5minBins = reshape (temp5minBins,[],1);
        rainbow_cell(3:(2*length(averagedStringVar))+2,i) = num2cell(longTemp5minBins);

        % Put the data in the rainbow matrices
        rainbow_mat(:,i) = temp5minBins(:,1);
        rainbow_mat_sem(:,i) = temp5minBins(:,2);

    end

    % Create the rainbox plots
    figure('Color', [1 1 1 ]);

    % setting the color scheme of the rainbow plot
    set(gcf,'Colormap',cbrewer('seq','PuBuGn',9));
    set(gcf,'DefaultAxesColorOrder',cbrewer('qual','Set2',8))
    set(gcf,'DefaultLineLineWidth',1.2)

    % Plotting
    %Use roundn to round to multiple of ten
    mseb(1:length(averagedStringVar),rainbow_mat',rainbow_mat_sem');
    %axis([1,length(averagedStringVar)+1,0,roundn(maxAct,1)]) 
    %set y axis here
    axis([1,length(averagedStringVar)+1,0,roundn(30,1)])
    set(gca,'XTick',1:1:length(averagedStringVar))
    set(gca,'XTickLabel',{xAxisMat})
    legend({master_data_struct(geno_indices_of_the_current_rainbowgroup).genotype},'Location', 'NorthEast')
    xlabel('Time')
    ylabel(['Activity / Beam crosses per' num2str(binLength) ' (min)'])

    % Save the fig and the data
    saveas(gcf, fullfile(export_path,[filename_master(1:end-5),'_',num2str(rainbowgroups_unique(j)),'_ULOSRactivityPlot.pdf']));
    savefig(fullfile(export_path,[filename_master(1:end-5),'_',num2str(rainbowgroups_unique(j)),'_ULOSRactivityPlot.fig']));
    cell2csv(fullfile(export_path,[filename_master(1:end-5),'_',num2str(rainbowgroups_unique(j)),'_ULOSRactivityPlot.csv']),rainbow_cell)
    close gcf
       
end

% Save the work space
save(fullfile(export_path,[filename_master(1:end-5),'_ULOSRactivityworkspace.mat']));
