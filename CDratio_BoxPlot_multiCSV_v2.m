% V2: for multi group plotting & statistic exam in sub functions

close all;

keyword = 'FA2_*.csv';    % apparatus name

Write2CSV =1;

listing = dir(keyword);
name = string(extractfield(listing, 'name'));
fprintf('Data in file:\n');
name_seq = [[1:length(name)]', name'];
disp(name_seq);
close 
for ii = 1:length(name)
    
    %% File input
    TableName = name(ii);
    T = readtable(TableName,'NumHeaderLines',1);
    %     T1 = readmatrix(TableName,'NumHeaderLines',1);
    T = table2array(T);
    
%     % Look into different FOV 
%     groupIndex = T(:,1);
%     Index = 3;
%     CD_Cy3 = T(groupIndex == Index,3);
%     CH_Cy3 = T(groupIndex == Index,4);
%     off_Cy3 = T(groupIndex == Index,5);
%     CD_Cy5 = T(groupIndex == Index,6);
%     CH_Cy5 = T(groupIndex == Index,7);
%     off_Cy5 = T(groupIndex == Index,8);
%     CD_allBac = T(groupIndex == Index,9);
%     CH_allBac = T(groupIndex == Index,10);
%     off_allBac = T(groupIndex == Index,11);
    
    CD_Cy3 = T(:,3);
    CH_Cy3 = T(:,4);
    off_Cy3 = T(:,5);
    
    CD_Cy5 = T(:,6);
    CH_Cy5 = T(:,7);
    off_Cy5 = T(:,8);
    
    CD_allBac = T(:,9);
    CH_allBac = T(:,10);
    off_allBac = T(:,11);
    
   % use neg ctrl sample 4 to calibrate the compensation coeff
    RescaleCoeff1 = 1.1;
    off_Cy3 = off_Cy3*RescaleCoeff1;
    CD_Cy3 = CD_Cy3 - off_Cy3;
    CH_Cy3 = CH_Cy3 - off_Cy3;
    RescaleCoeff2 = 1;
    CH_Cy3 = CH_Cy3*RescaleCoeff2;
    CH_low_Cy3 = 0.15;
    
    RescaleCoeff3 = 0.95;
    off_Cy5 = off_Cy5*RescaleCoeff3;
    CD_Cy5 = CD_Cy5 - off_Cy5;
    CH_Cy5 = CH_Cy5 - off_Cy5;
    RescaleCoeff4 = RescaleCoeff2;
    CH_Cy5 = CH_Cy5*RescaleCoeff4;
    CH_low_Cy5 = 0.2;
    
    RescaleCoeff5 = 1;
    off_allBac = off_allBac*RescaleCoeff5;
    CD_allBac = CD_allBac - off_allBac;
    CH_allBac = CH_allBac - off_allBac;
    RescaleCoeff6 = RescaleCoeff2;
    CH_allBac = CH_allBac*RescaleCoeff6;
    CH_low_allBac = 0.15;
    
    
    % normal
    CD_percent_Cy3 = StatisticBoxPlot(CD_Cy3, CH_Cy3, CH_low_Cy3,2);
    CD_percent_Cy5 = StatisticBoxPlot(CD_Cy5, CH_Cy5, CH_low_Cy5,2);
    CD_percent_allBac = StatisticBoxPlot(CD_allBac, CH_allBac, CH_low_allBac,2);
    
    pause(5);
    
    %% Save to original csv file
    if Write2CSV == 1
    Tupdate = array2table(T);
    Tupdate(:,12) = array2table(CD_percent_Cy3);
    Tupdate(:,13) = array2table(CD_percent_Cy5);
    Tupdate(:,14) = array2table(CD_percent_allBac);
    writetable(Tupdate, TableName);
    end
    
    %% Plotting default settings
    
    set(0, 'DefaultTextFontSize', 18);
    set(0, 'DefaultLineLineWidth', 2);
    set(0, 'DefaultAxesLIneWidth', 2);
    set(0, 'DefaultAxesFontSize', 14);
    set(0, 'DefaultAxesFontName', 'Arial');
    
end



function CD_percent_save = StatisticBoxPlot(CD, CH, CH_PrecisionLow, Range_const1, Range_const2)

    if nargin < 5
        Range_const2 = 2;
        if nargin < 4
            Range_const1 = 2;
            if nargin < 3
            CH_PrecisionLow = 0;
            end
        end
    end
        
    Length_vec = length(CD);
    Length_real = sum(~isnan(CD));
    
    CD = CD(find(~isnan(CD)));
    CH = CH(find(~isnan(CH)));
    
    % clean up NaN value when dealing with multiple kinds of objects
    
%% Statistics

%% Primary mask for rejecting extreme case
    CD_percent_raw  = (CD)./(CD+CH)*100;
    CDratio_LowLim0 = -10;
    CDratio_UpLim0 = 50;
    RatioMask_0 = (CD_percent_raw > CDratio_LowLim0).*(CD_percent_raw < CDratio_UpLim0);
    RatioMask_0 = RatioMask_0.*(CH > CH_PrecisionLow);
     
    CD_0  = CD.*RatioMask_0;
    CD_0 = CD_0(CD_0~=0);
    
    % off_0  = off.*CD_mask1.*off_mask1.*CH_mask1;
    % off_0 = off_0(off_0~=0);
    
    CH_0 = CH.*RatioMask_0;
    CH_0 = CH_0(CH_0~=0);
    
    Mean_CD1 = mean(CD_0);
    Std_CD1 = std(CD_0);
    % Mean_off2 = mean(off_0);
    % Std_off2 = std(off_0);
    Mean_CH1 = mean(CH_0);
    Std_CH1 = std(CH_0);
    
    %% Result masks 1 for abnormal results
%     Range_const1 = 2;
    
    % CD mask
    CD_LowLim1 = Mean_CD1 - Range_const1*Std_CD1;
    CD_UpLim1 = Mean_CD1 + Range_const1*Std_CD1;
    CD_mask1 = (CD > CD_LowLim1).*(CD < CD_UpLim1);
    
    % % off mask
    % off_LowLim1 = Mean_off1 - Range_const1*Std_off1;
    % off_UpLim1 = Mean_off1 + Range_const1*Std_off1;
    % off_mask1 = (off > off_LowLim1).*(off < off_UpLim1);
    
    % CH mask
    CH_LowLim1 = Mean_CH1 - Range_const1*Std_CH1;
    CH_UpLim1 = Mean_CH1 + Range_const1*Std_CH1;
    CH_mask1 = (CH > CH_PrecisionLow).*(CH > CH_LowLim1).*(CH < CH_UpLim1);
    
    %% Restatistics after picked out outliers
    CD_normal  = CD.*CD_mask1.*CH_mask1;
    CD_normal = CD_normal(CD_normal~=0);
    
    % off_normal  = off.*CD_mask1.*off_mask1.*CH_mask1;
    % off_normal = off_normal(off_normal~=0);
    
    CH_normal = CH.*CD_mask1.*CH_mask1;
    CH_normal = CH_normal(CH_normal~=0);
    
    Mean_CD2 = mean(CD_normal);
    Std_CD2 = std(CD_normal);
    % Mean_off2 = mean(off_normal);
    % Std_off2 = std(off_normal);
    Mean_CH2 = mean(CH_normal);
    Std_CH2 = std(CH_normal);
    
    
    %% Result masks 2 for outliers (did not precisely measured due to the sample shape/movement)
%     Range_const2 = 2;   % pure culture
    
    % CD mask
    CD_LowLim2 = Mean_CD2 - Range_const2*Std_CD2;
    CD_UpLim2 = Mean_CD2 + Range_const2*Std_CD2;
    CD_mask2 = (CD > CD_LowLim2).*(CD < CD_UpLim2);
    
    % % off mask
    % off_LowLim2 = Mean_off2 - Range_const2*Std_off2;
    % off_UpLim2 = Mean_off2 + Range_const2*Std_off2;
    % off_mask2 = (off > off_LowLim2).*(off < off_UpLim2);
    
    % CH mask
    CH_LowLim2 = Mean_CH2 - Range_const2*Std_CH2;
    CH_UpLim2 = Mean_CH2 + Range_const2*Std_CH2;
    CH_mask2 = (CH > CH_PrecisionLow).*(CH > CH_LowLim2).*(CH < CH_UpLim2);
    
    
    %% Masked CD ratio result
%     CD_percent_raw  = (CD)./(CD+CH)*100;
    CD_percent_0 = CD_percent_raw.*CD_mask1.*CH_mask1;
    CD_percent_0 = CD_percent_0(CD_percent_0~=0);
    CD_percent_masked = CD_percent_raw.*CD_mask2.*CH_mask2;
    CD_percent = CD_percent_masked(CD_percent_masked~=0);
    CD_percent_save = CD_percent_masked;
    CD_percent_save(CD_percent_save==0)= NaN;
    CD_percent_save = [CD_percent_save;NaN(Length_vec-Length_real,1)];
    
    %% Quick boxplot check
    figure;
    subplot(1,3,1);
    boxplot([CD_percent_raw]);
    hold on;
    scatter(ones(size(CD_percent_raw)).*(1+(rand(size(CD_percent_raw))-0.5)/10),CD_percent_raw,'r','filled')
    hold off;
%     ylim([-7 50]); 
    
    subplot(1,3,2);
    boxplot([CD_percent_0]);
    hold on;
    scatter(ones(size(CD_percent_0)).*(1+(rand(size(CD_percent_0))-0.5)/10),CD_percent_0,'r','filled')
    hold off;
    
    subplot(1,3,3);
    boxplot([CD_percent]);
    hold on;
    scatter(ones(size(CD_percent)).*(1+(rand(size(CD_percent))-0.5)/10),CD_percent,'r','filled')
    hold off;
    
    % print for quick quality check
    fprintf('\n');
    fprintf('CD percent std: %.3f \n',std(CD_percent));
    fprintf('CD percent mean: %.3f \n',mean(CD_percent));
    fprintf('Number of cell counted: %d \n', length(CD_percent));
    fprintf('Total cell number: %d \n\n', length(CD_percent_raw));
end