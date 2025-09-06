clc; 
close all; 
clear all;

records_dfc = dir("C:\Melbourne University\Study Four-DFC\draft\code and data\dynFC\"); 
win_len = 18; %window length
win_step = 3; %window step

dir_win = 'C:\Melbourne University\Study Four-DFC\draft\code and data\windows_time18_3\'; %files including the moving window time intervals
dir_onsets = 'C:\Melbourne University\Study Four-DFC\draft\code and data\processed_data\'; %window onsets

subs_CI = 1:37;
subs_CI([3,12,20,25,26,30]) = []; %CI users did not complete the experiment
subs_NH = [5,6,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28]; % NH subjects (Bad ones are removed)

T_chs = readmatrix("ch_config.csv");    %channles' ristricted rois
schs = [2,4,14,23,30,36,49,52];     %short channels
T_chs(schs,:) = [];    % remove short channels from the list
T_chs(find(T_chs(:,4)==6),4) = 5;  % relable channels to allocate the same labels for left and right channels in the same ROI
T_chs(find(T_chs(:,4)==7),4) = 6; % 1.L-STG 2.L-AG 3.R-STG 4.R-AG 5.OC 6. IFG
load("chsPos.mat")

%% reading dynFC networks & flex-in-time calculations
%NH
NH_rflex = []; %resting-state flexibility
NH_ralleg = []; %resting-state allegience
NH_tflex = [];  %task flexibility
NH_flexAveTask = []; %average task flexibility
NH_talleg = []; %task allegience
NH_rid = [];    %ID of resting-state recordings of normal hearing subjects
NH_tid = [];    %ID of task recordings of normal hearing subjects
NH_roiISIflex = []; %isi flexibility of different rois 
NH_chs_flex = [];
ctl_NH = 0;
%CI1
CI1_rflex = [];
CI1_ralleg = [];
CI1_tflex = [];
CI1_flexAveTask = [];
CI1_talleg = [];
CI1_rid = [];
CI1_tid = [];
CI1_roiISIflex = [];
CI1_chs_flex = [];
ctl_CI1 = 0;
%CI2
CI2_rflex = [];
CI2_ralleg = [];
CI2_tflex = [];
CI2_flexAveTask = [];
CI2_talleg = [];
CI2_rid = [];
CI2_tid = [];
CI2_roiISIflex = [];
CI2_chs_flex = [];
ctl_CI2 = 0;

%read dynFC files for subjects
for i = 3:numel(records_dfc)
    name_record = records_dfc(i).name;
    %NH
    if (strcmp(name_record(1:2), 'NH')) && (ismember(str2num(name_record(end-5:end-4)),subs_NH))
        load([records_dfc(i).folder '\' name_record]); %read dynFC files for the subject: F_: flexibility; Sopt: detected communities in each window; C_: supposed to include allegience but they donot. We calculate allegience here again
        [F_opt44, F_tempNull44, A_opt44, A_nodNull44, Sopt44 ] = resizeTo44(F_opt, F_tempNull, C_opt, C_nodNull, Sopt);  % resize mats for 44 channels to match the sizes    
        switch name_record(3)
            case 'r'
                NH_rflex = cat(2,NH_rflex,F_opt44);
                NH_ralleg = cat(3,NH_ralleg,A_opt44);
                NH_rid = cat(1,NH_rid,str2num(name_record(end-5:end-4)));
            case 't'
                ctl_NH = ctl_NH + 1;
                NH_tflex = cat(2,NH_tflex,F_opt44); 
                NH_tid = cat(1,NH_tid,str2num(name_record(end-5:end-4)));
                win_time = readmatrix([dir_win 'NH' name_record(end-5:end-4)]);
                onsets = readmatrix([dir_onsets 'NHonsets_' name_record(end-5:end-4)]);  
                [flex_ts, flex_ts44, flex_t, flex_l] = flexInTime(Sopt, F_opt(:,2), win_time, win_len, win_step, onsets); % flexibility-in-time calculation
                flex_tsNH( ctl_NH, : ) = mean(flex_ts(:,1:600));
                aveTrialFlexNH{ctl_NH} = ave_trial_flex(flex_l, mean(flex_ts));
                NH_ComSim{ctl_NH} = part_similarityMat(Sopt, flex_l);
                [audAve, visAve, ctlAve, isiAve, trnAve] = flexAveTask(flex_ts,flex_l);     % average flexibility in different conditions
                NH_flexAveTask = cat(1,NH_flexAveTask, [audAve, visAve, ctlAve, isiAve, trnAve]);
                [allChs_flex, LSTG_flex, LAG_flex, RSTG_flex, RAG_flex, OC_flex, IFG_flex] = ROIaveFlexISI(flex_ts, flex_l, T_chs(:,4), F_opt(:,2)); %average ROI flexibility
                NH_roiISIflex = cat(1, NH_roiISIflex, [allChs_flex, LSTG_flex, LAG_flex, RSTG_flex, RAG_flex, OC_flex, IFG_flex]);
                NH_talleg = cat(4, NH_talleg, alleg_cons(Sopt, flex_l, F_opt(:,2)));    %alleg_cons: allegience matrices during 1. whole recording, 2. auditory task, 3. visual task, 4. contro, 5. isi, 6. transition 
                % plot_flexInTime(flex_ts, flex_t, flex_l, onsets,name_record), pause, close all           % uncomment to draw flexibility in time 
                NH_chs_flex = [NH_chs_flex mean(flex_ts44,2)];
        end
    end
    %CI 1
    if ~strcmp(name_record(1:2), 'NH') && strcmp(name_record(end-6), '1') && ismember(str2num(name_record(end-5:end-4)),subs_CI)
        load([records_dfc(i).folder '\' name_record]);  
        [F_opt44, F_tempNull44, A_opt44, A_nodNull44, Sopt44 ] = resizeTo44(F_opt, F_tempNull, C_opt, C_nodNull, Sopt);     
        switch name_record(1)
            case 'r'                
                CI1_rflex = cat(2,CI1_rflex,F_opt44);
                CI1_ralleg = cat(3,CI1_ralleg,A_opt44);
                CI1_rid = cat(1,CI1_rid,str2num(name_record(end-5:end-4)));
            case 't'
                ctl_CI1 = ctl_CI1 + 1;
                CI1_tflex = cat(2,CI1_tflex,F_opt44);                
                CI1_tid = cat(1,CI1_tid,str2num(name_record(end-5:end-4)));
                win_time = readmatrix([dir_win 'CI1' name_record(end-5:end-4)]);
                onsets = readmatrix([dir_onsets 'onsets_1' name_record(end-5:end-4)]);
                [flex_ts, flex_ts44, flex_t, flex_l] = flexInTime(Sopt, F_opt(:,2), win_time, win_len, win_step, onsets);  % flexibility-in-time calculation
                flex_tsCI1( ctl_CI1, : ) = mean(flex_ts(:,1:600));
                aveTrialFlexCI1{ctl_CI1} = ave_trial_flex(flex_l, mean(flex_ts));
                CI1_ComSim{ctl_CI1} = part_similarityMat(Sopt, flex_l);
                [audAve, visAve, ctlAve, isiAve, trnAve] = flexAveTask(flex_ts,flex_l);
                CI1_flexAveTask = cat(1,CI1_flexAveTask, [audAve, visAve, ctlAve, isiAve, trnAve]);
                [allChs_flex, LSTG_flex, LAG_flex, RSTG_flex, RAG_flex, OC_flex, IFG_flex] = ROIaveFlexISI(flex_ts, flex_l, T_chs(:,4), F_opt(:,2));
                CI1_roiISIflex = cat(1, CI1_roiISIflex, [allChs_flex, LSTG_flex, LAG_flex, RSTG_flex, RAG_flex, OC_flex, IFG_flex]);
                CI1_talleg = cat(4, CI1_talleg, alleg_cons(Sopt, flex_l, F_opt(:,2)));
                % plot_flexInTime(flex_ts, flex_t, flex_l, onsets, name_record), pause, close all                %@@@@@@@@@@@ plot flex-in-time @@@@@@@@@@@@ 
                CI1_chs_flex = [CI1_chs_flex mean(flex_ts44,2)];
        end 
    end
    %CI 2
    if ~strcmp(name_record(1:2), 'NH') && strcmp(name_record(end-6), '2') && ismember(str2num(name_record(end-5:end-4)),subs_CI)
        load([records_dfc(i).folder '\' name_record]);
        [F_opt44, F_tempNull44, A_opt44, A_nodNull44, Sopt44 ] = resizeTo44(F_opt, F_tempNull, C_opt, C_nodNull, Sopt);       
        switch name_record(1)
            case 'r'
                CI2_rflex = cat(2,CI2_rflex,F_opt44);
                CI2_ralleg = cat(3,CI2_ralleg,A_opt44);
                CI2_rid = cat(1,CI2_rid,str2num(name_record(end-5:end-4)));
            case 't'
                ctl_CI2 = ctl_CI2 + 1;
                CI2_tflex = cat(2,CI2_tflex,F_opt44);
                CI2_tid = cat(1,CI2_tid,str2num(name_record(end-5:end-4)));
                win_time = readmatrix([dir_win 'CI2' name_record(end-5:end-4)]);
                onsets = readmatrix([dir_onsets 'onsets_2' name_record(end-5:end-4)]);
                [flex_ts, flex_ts44, flex_t, flex_l] = flexInTime(Sopt, F_opt(:,2), win_time, win_len, win_step, onsets);  % flexibility-in-time calculation
                flex_tsCI2( ctl_CI2, : ) = mean(flex_ts(:,1:600));
                aveTrialFlexCI2{ctl_CI2} = ave_trial_flex(flex_l, mean(flex_ts));
                CI2_ComSim{ctl_CI2} = part_similarityMat(Sopt, flex_l);
                [audAve, visAve, ctlAve, isiAve, trnAve] = flexAveTask(flex_ts,flex_l);
                CI2_flexAveTask = cat(1,CI2_flexAveTask, [audAve, visAve, ctlAve, isiAve, trnAve]);
                [allChs_flex, LSTG_flex, LAG_flex, RSTG_flex, RAG_flex, OC_flex, IFG_flex] = ROIaveFlexISI(flex_ts, flex_l, T_chs(:,4), F_opt(:,2));
                CI2_roiISIflex = cat(1, CI2_roiISIflex, [allChs_flex, LSTG_flex, LAG_flex, RSTG_flex, RAG_flex, OC_flex, IFG_flex]);
                CI2_talleg = cat(4, CI2_talleg, alleg_cons(Sopt, flex_l, F_opt(:,2)));
                % plot_flexInTime(flex_ts, flex_t, flex_l, onsets, name_record), pause, close all   close all      %@@@@@@@@@@@ plot flex-in-time @@@@@@@@@@@@
                CI2_chs_flex = [CI2_chs_flex mean(flex_ts44,2)];
        end
    end
end

%% Whole recording: ci flex & speech performance correlation (Alleginace is also processed here)
% resize the results to match the subject IDs(#37)
CI1_rflex37 = NaN(44,37);
CI1_rflex37(:,CI1_rid) = CI1_rflex;
CI1_tflex37 = NaN(44,37);
CI1_tflex37(:,CI1_tid) = CI1_tflex;
CI1_tflex37(:,13) = nan; % bad qualoty signal sub 13
CI1_flexAveTask37 = NaN(37,5);
CI1_flexAveTask37(CI1_tid,:) = CI1_flexAveTask;
CI1_flexAveTask37(13,:) = nan; % bad qualoty signal sub 13
CI1_roiISIflex37 = NaN(37,7);
CI1_roiISIflex37(CI1_tid,:) = CI1_roiISIflex;
CI1_roiISIflex37(13,:) = nan; % bad qualoty signal sub 13
CI2_rflex37 = NaN(44,37);
CI2_rflex37(:,CI2_rid) = CI2_rflex;
CI2_tflex37 = NaN(44,37);
CI2_tflex37(:,CI2_tid) = CI2_tflex;
CI2_flexAveTask37 = NaN(37,5);
CI2_flexAveTask37(CI2_tid,:) = CI2_flexAveTask;
CI2_roiISIflex37 = NaN(37,7);
CI2_roiISIflex37(CI2_tid,:) = CI2_roiISIflex;

% average flexibility & correlation with speech scores
NH_rflex_ave = mean(NH_rflex,"omitmissing");
CI1_rflex_ave = mean(CI1_rflex37,"omitmissing");
CI2_rflex_ave = mean(CI2_rflex37,"omitmissing");
NH_tflex_ave = mean(NH_tflex,"omitmissing");
CI1_tflex_ave = mean(CI1_tflex37,"omitmissing");
CI2_tflex_ave = mean(CI2_tflex37,"omitmissing");
CI1_tflex_ave(13) = nan; % bad qualoty signal sub 13
% correlation of average flexibility with speech scores
scores = load("Speech Scores_session two.txt");
scores(:,2:5) = rau_convert(scores(:,2:5)); %rau transormation on speech scores
scores37 = NaN(37,6);
scores37(scores(:,1),:) = scores;
scores37(2,2:5) = nan; %subject 2 did not complete speech test
%session one
figure('Name','correlations')
subplot(2,2,1)
lm = fitlm(CI1_tflex_ave ,scores37(:,3));
plot(lm)
xlim([0.1 0.2])
title(['fNIRS 1: p  = ' num2str(round(lm.ModelFitVsNullModel.Pvalue,3)) ', R=' num2str(round(sqrt(lm.Rsquared.Ordinary),2))])
ylabel('CNC phons test scores')
xlabel("average flexibility")
ylim([55 100])
subplot(2,2,2)
lm = fitlm(CI1_tflex_ave ,scores37(:,5));
plot(lm)
xlim([0.1 0.2])
title(['fNIRS 1: p  = ' num2str(round(lm.ModelFitVsNullModel.Pvalue,3)) ', R=' num2str(round(sqrt(lm.Rsquared.Ordinary),2))])
ylabel('BKB 10 test scores')
xlabel("average flexibility")
subplot(2,2,3)
lm = fitlm(CI2_tflex_ave ,scores37(:,3));
plot(lm)
xlim([0.1 0.2])
title(['fNIRS 2: p  = ' num2str(round(lm.ModelFitVsNullModel.Pvalue,3)) ', R=' num2str(round(sqrt(lm.Rsquared.Ordinary),2))])
ylabel('CNC phons test scores')
xlabel("average flexibility")
ylim([55 100])
subplot(2,2,4)
lm = fitlm(CI2_tflex_ave ,scores37(:,5));
plot(lm)
xlim([0.1 0.2])
title(['fNIRS 2: p  = ' num2str(round(lm.ModelFitVsNullModel.Pvalue,3)) ', R=' num2str(round(sqrt(lm.Rsquared.Ordinary),2))])
ylabel('BKB 10 test scores')
xlabel("average flexibility")

%% correlation of rois flexibilities with scores (whole recording)
CI1_ROIaveFlex = ROIaveFlex(CI1_tflex37, T_chs(:,4));
CI2_ROIaveFlex = ROIaveFlex(CI2_tflex37, T_chs(:,4));
CI1_ROIaveFlex = [CI1_tflex_ave' CI1_ROIaveFlex];
CI2_ROIaveFlex = [CI2_tflex_ave' CI2_ROIaveFlex];

ROIs = {'all montage', 'L_STG', 'L_AG', 'R_STG', 'R_AG', 'OC', 'IFG'};
for i  = 1:7
    figure('Name',['correlations for ' ROIs{i}])
    subplot(2,2,1)
    lm = fitlm(CI1_ROIaveFlex(:,i) ,scores37(:,3));
    plot(lm)
    title(['fNIRS 1: p  = ' num2str(round(lm.ModelFitVsNullModel.Pvalue,3)) ', R=' num2str(round(sqrt(lm.Rsquared.Ordinary),2))])
    ylabel('CNC phons test scores')
    xlabel("average recording flexibility")
    ylim([55 100])
    subplot(2,2,2)
    lm = fitlm(CI1_ROIaveFlex(:,i) ,scores37(:,5));
    plot(lm)
    title(['fNIRS 1: p  = ' num2str(round(lm.ModelFitVsNullModel.Pvalue,3)) ', R=' num2str(round(sqrt(lm.Rsquared.Ordinary),2))])
    ylabel('BKB 10 test scores')
    xlabel("average recording flexibility")
    subplot(2,2,3)
    lm = fitlm(CI2_ROIaveFlex(:,i) ,scores37(:,3));
    plot(lm)
    title(['fNIRS 2: p  = ' num2str(round(lm.ModelFitVsNullModel.Pvalue,3)) ', R=' num2str(round(sqrt(lm.Rsquared.Ordinary),2))])
    ylabel('CNC phons test scores')
    xlabel("average recording flexibility")
    ylim([55 100])
    subplot(2,2,4)
    lm = fitlm(CI2_ROIaveFlex(:,i) ,scores37(:,5));
    plot(lm)
    title(['fNIRS 2: p  = ' num2str(round(lm.ModelFitVsNullModel.Pvalue,3)) ', R=' num2str(round(sqrt(lm.Rsquared.Ordinary),2))])
    ylabel('BKB 10 test scores')
    xlabel("average recording flexibility")
end

%% channels' flexibility correlation with scores
for i = 1:44    
    lm = fitlm(CI1_tflex37(i,:) ,scores37(:,4));
    Rch11(i) = sqrt(lm.Rsquared.Ordinary);
    lm = fitlm(CI2_tflex37(i,:) ,scores37(:,4));
    Rch21(i) = sqrt(lm.Rsquared.Ordinary);
    lm = fitlm(CI1_tflex37(i,:) ,scores37(:,5));
    Rch12(i) = sqrt(lm.Rsquared.Ordinary);
    lm = fitlm(CI2_tflex37(i,:) ,scores37(:,5));
    Rch22(i) = sqrt(lm.Rsquared.Ordinary);
end

figure
subplot(2,2,1)
plot_chSignificance(Rch11, chsPos)
title('fNIRS1-CNC')
subplot(2,2,2)
plot_chSignificance(Rch12, chsPos)
title('fNIRS1-BKB')
subplot(2,2,3)
plot_chSignificance(Rch21, chsPos)
title('fNIRS2-CNC')
subplot(2,2,4)
plot_chSignificance(Rch22, chsPos)
title('fNIRS2-BKB')

%% Different conditions in the recording: ci flex & speech performance correlation
tasks = {'aud'; 'vis'; 'ctl'; 'isi';'trans'};
for i = 1:5
    figure('Name',['correlations for ' tasks{i,1}])
    subplot(2,2,1)
    lm = fitlm(CI1_flexAveTask37(:,i) ,scores37(:,3));
    plot(lm)
    title(['p  = ' num2str(round(lm.ModelFitVsNullModel.Pvalue,4)) ', R=' num2str(round(sqrt(lm.Rsquared.Ordinary),3))])
    ylabel('speech score')
    xlabel("average flexibility")
    ylim([55 100])
    xlim([0.1 0.2])
    subplot(2,2,2)
    lm = fitlm(CI1_flexAveTask37(:,i) ,scores37(:,5));
    plot(lm)
    xlim([0.1 0.2])
    title(['p  = ' num2str(round(lm.ModelFitVsNullModel.Pvalue,4)) ', R=' num2str(round(sqrt(lm.Rsquared.Ordinary),3))])
    ylabel('BKB 10 test scores')
    xlabel("average flexibility")
    subplot(2,2,3)
    lm = fitlm(CI2_flexAveTask37(:,i) ,scores37(:,3));
    plot(lm)
    title(['p  = ' num2str(round(lm.ModelFitVsNullModel.Pvalue,4)) ', R=' num2str(round(sqrt(lm.Rsquared.Ordinary),3))])
    ylabel('speech score')
    xlabel("average flexibility")
    ylim([55 100])
    xlim([0.1 0.2])
    subplot(2,2,4)
    lm = fitlm(CI2_flexAveTask37(:,i) ,scores37(:,5));
    plot(lm)
    xlim([0.1 0.2])
    title(['p  = ' num2str(round(lm.ModelFitVsNullModel.Pvalue,4)) ', R=' num2str(round(sqrt(lm.Rsquared.Ordinary),3))])
    ylabel('speech score')
    xlabel("average flexibility")
end

%% correlation between ROIs isi flexibility with scores
ROIs = {'all montage', 'L_STG', 'L_AG', 'R_STG', 'R_AG', 'OC', 'IFG'};
for i  = 1:7
    figure('Name',['correlations for ' ROIs{i}])
    subplot(2,2,1)
    lm = fitlm(CI1_roiISIflex37(:,i) ,scores37(:,3));
    plot(lm)
    title(['fNIRS 1: p  = ' num2str(lm.ModelFitVsNullModel.Pvalue) ', R=' num2str(sqrt(lm.Rsquared.Ordinary))])
    ylabel('CNC phons test scores')
    xlabel("average isi flexibility")
    ylim([55 100])
    subplot(2,2,2)
    lm = fitlm(CI1_roiISIflex37(:,i) ,scores37(:,5));
    plot(lm)
    title(['fNIRS 1: p  = ' num2str(lm.ModelFitVsNullModel.Pvalue) ', R=' num2str(sqrt(lm.Rsquared.Ordinary))])
    ylabel('BKB 10 test scores')
    xlabel("average isi flexibility")
    subplot(2,2,3)
    lm = fitlm(CI2_roiISIflex37(:,i) ,scores37(:,3));
    plot(lm)
    title(['fNIRS 2: p  = ' num2str(lm.ModelFitVsNullModel.Pvalue) ', R=' num2str(sqrt(lm.Rsquared.Ordinary))])
    ylabel('CNC phons test scores')
    xlabel("average isi flexibility")
    ylim([55 100])
    subplot(2,2,4)
    lm = fitlm(CI2_roiISIflex37(:,i) ,scores37(:,5));
    plot(lm)
    title(['fNIRS 2: p  = ' num2str(lm.ModelFitVsNullModel.Pvalue) ', R=' num2str(sqrt(lm.Rsquared.Ordinary))])
    ylabel('BKB 10 test scores')
    xlabel("average isi flexibility")
end

%% evaluate allegiance based on structural ROIs
% 1.L-STG 2.L-AG 3.R-STG 4.R-AG 5.OC 6.IFG
%ROIs integration
condition = 1; %whole recording
broad_ROIs = [1*ones(1,9) 2*ones(1,14) 3*ones(1,14), 4*ones(1,7)];
for condition = 1:6
for i = 1:23
    sub_integration = integration(NH_talleg(:,:,condition,i), broad_ROIs); 
    NH_IFGintegration(i, condition) = mean(sub_integration(1:9),"omitmissing");
    NH_STGintegration(i, condition) = mean(sub_integration(T_chs(:,4) == 2 | T_chs(:,4) == 4 ),"omitmissing");
    NH_AGintegration(i, condition) = mean(sub_integration(T_chs(:,4) == 3 | T_chs(:,4) == 5 ),"omitmissing");
    NH_VCintegration(i, condition) = mean(sub_integration(38:44),"omitmissing");
    NH_consensusCom(i, condition,:) = consensus_com(NH_talleg(:,:,condition,i));
end
for i = 1:31
    sub_integration = integration(CI1_talleg(:,:,condition,i), broad_ROIs);
    CI1_IFGintegration(i, condition) = mean(sub_integration(1:9),'omitmissing');
    CI1_STGintegration(i, condition) = mean(sub_integration(T_chs(:,4) == 2 | T_chs(:,4) == 4 ),"omitmissing");
    CI1_AGintegration(i, condition) = mean(sub_integration(T_chs(:,4) == 3 | T_chs(:,4) == 5 ),"omitmissing");    
    CI1_VCintegration(i, condition) = mean(sub_integration(38:44),"omitmissing");
    CI1_consensusCom(i, condition,:) = consensus_com(CI1_talleg(:,:,condition,i));
    sub_integration = integration(CI2_talleg(:,:,condition,i), broad_ROIs);
    CI2_IFGintegration(i, condition) = mean(sub_integration(1:9),'omitmissing');
    CI2_STGintegration(i, condition) = mean(sub_integration(T_chs(:,4) == 2 | T_chs(:,4) == 4 ),"omitmissing");
    CI2_AGintegration(i, condition) = mean(sub_integration(T_chs(:,4) == 3 | T_chs(:,4) == 5 ),"omitmissing");   
    CI2_VCintegration(i, condition) = mean(sub_integration(38:44),"omitmissing");
    CI2_consensusCom(i, condition,:) = consensus_com(CI2_talleg(:,:,condition,i));
end
end

%% consensus communities
NH_ci_aud = consensus_und( mean(NH_talleg(:,:,2,:), 4,"omitnan"), 0,50);
NH_ci_vis = consensus_und( mean(NH_talleg(:,:,3,:), 4,"omitnan"), 0,50);
NH_ci_ctl = consensus_und( mean(NH_talleg(:,:,4,:), 4,"omitnan"), 0,50);
NH_ci_isi = consensus_und( mean(NH_talleg(:,:,5,:), 4,"omitnan"), 0,50);

CI1_ci_aud = consensus_und( mean(CI1_talleg(:,:,2,:), 4,"omitnan"), 0,50);
CI1_ci_vis = consensus_und( mean(CI1_talleg(:,:,3,:), 4,"omitnan"), 0,50);
CI1_ci_ctl = consensus_und( mean(CI1_talleg(:,:,4,:), 4,"omitnan"), 0,50);
CI1_ci_isi = consensus_und( mean(CI1_talleg(:,:,5,:), 4,"omitnan"), 0,50);

CI2_ci_aud = consensus_und( mean(CI2_talleg(:,:,2,:), 4,"omitnan"), 0,50);
CI2_ci_vis = consensus_und( mean(CI2_talleg(:,:,3,:), 4,"omitnan"), 0,50);
CI2_ci_ctl = consensus_und( mean(CI2_talleg(:,:,4,:), 4,"omitnan"), 0,50);
CI2_ci_isi = consensus_und( mean(CI2_talleg(:,:,5,:), 4,"omitnan"), 0,50);



%% inter-hemesphere allegiance
for i = 1:23
    NH_ihemesphere_ac(i) = mean(NH_talleg(T_chs(:,4) == 1 | T_chs(:,4) == 2,T_chs(:,4) == 3 | T_chs(:,4) == 4, condition, i),"all","omitmissing");
    NH_ihemesphere_ag(i) = mean( NH_talleg( T_chs(:,4) == 2, T_chs(:,4) == 4 ,condition,i),"all","omitmissing");
    NH_ihemesphere_stg(i) = mean(NH_talleg(T_chs(:,4) == 1 ,T_chs(:,4) == 3 ,condition,i),"all","omitmissing");

end
for i = 1:31
    if ismember(i,[2,11]) % exclude subject 2 & 13
        CI1_ihemesphere_ac(i) = nan;
        CI1_ihemesphere_ag(i) = nan;
        CI1_ihemesphere_stg(i) = nan;
        CI2_ihemesphere_ac(i) = nan;
        CI2_ihemesphere_ag(i) = nan;
        CI2_ihemesphere_stg(i) = nan;
    else
        CI1_ihemesphere_ac(i) = mean(CI1_talleg(T_chs(:,4) == 1 | T_chs(:,4) == 2, T_chs(:,4) == 3 | T_chs(:,4) == 4, condition, i),"all","omitmissing");
        CI1_ihemesphere_ag(i) = mean( CI1_talleg( T_chs(:,4) == 2, T_chs(:,4) == 4 ,condition,i),"all","omitmissing");
        CI1_ihemesphere_stg(i) = mean(CI1_talleg(T_chs(:,4) == 1 ,T_chs(:,4) == 3 ,condition,i),"all","omitmissing");
        CI2_ihemesphere_ac(i) = mean(CI2_talleg(T_chs(:,4) == 1 | T_chs(:,4) == 2, T_chs(:,4) == 3 | T_chs(:,4) == 4, condition, i),"all","omitmissing");
        CI2_ihemesphere_ag(i) = mean( CI2_talleg( T_chs(:,4) == 2, T_chs(:,4) == 4 , condition, i),"all","omitmissing");
        CI2_ihemesphere_stg(i) = mean(CI2_talleg(T_chs(:,4) == 1 ,T_chs(:,4) == 3 , condition, i),"all","omitmissing");
    end
end
[p,~] = signrank(NH_ihemesphere_ag, NH_ihemesphere_stg)
[p,~] = signrank(CI1_ihemesphere_ag, CI1_ihemesphere_stg)
[p,~] = signrank(CI2_ihemesphere_ag, CI2_ihemesphere_stg)

% [p,~] = ranksum(NH_ihemesphere_ac, CI1_ihemesphere_ac)
[p,~] = ranksum(NH_ihemesphere_ag, CI1_ihemesphere_ag)
[p,~] = ranksum(NH_ihemesphere_ag, CI2_ihemesphere_ag)
[p,~] = signrank(CI2_ihemesphere_ag, CI1_ihemesphere_ag)

[p,~] = ranksum(NH_ihemesphere_stg, CI1_ihemesphere_stg)
[p,~] = ranksum(NH_ihemesphere_stg, CI2_ihemesphere_stg)
[p,~] = signrank(CI2_ihemesphere_stg, CI1_ihemesphere_stg)


CI1_ih = [CI1_ihemesphere_ag', CI1_ihemesphere_stg' ];
CI2_ih = [CI2_ihemesphere_ag', CI2_ihemesphere_stg' ];
NH_ih = NaN(size(CI1_ih));
NH_ih(1:23,:) = [NH_ihemesphere_ag', NH_ihemesphere_stg' ];
figure('Name','interhemispheric allegiance')
boxplot([NH_ih(:,1), CI1_ih(:,1), CI2_ih(:,1), NaN(31,1), NH_ih(:,2), CI1_ih(:,2), CI2_ih(:,2)], 'Labels',{'NH','CI1', 'CI2','','NH','CI1', 'CI2'}), hold on
dx1 = 0.12*randn(31,1);
dx2 = 0.12*randn(31,1);
dx3 = 0.12*randn(31,1);
plot(ones(31,1)+dx1,NH_ih(:,1),'o',"MarkerSize",3,"MarkerFaceColor",[0 0 1],"MarkerEdgeColor",[0.8 0.8 0.8])
plot(2*ones(31,1)+dx2,CI1_ih(:,1),'o',"MarkerSize",3,"MarkerFaceColor",[0 0 1],"MarkerEdgeColor",[0.8 0.8 0.8])
plot(3*ones(31,1)+dx3,CI2_ih(:,1),'o',"MarkerSize",3,"MarkerFaceColor",[0 0 1],"MarkerEdgeColor",[0.8 0.8 0.8])
for i = 1:31
    plot([2+dx2(i) 3+dx3(i)], [CI1_ih(i,1) CI2_ih(i,1)],"Color",[0.8 0.8 0.8])
end
plot([2,3], [mean(CI1_ih(:,1)) mean(CI2_ih(:,1))],'r')
plot(5*ones(31,1)+dx1,NH_ih(:,2),'o',"MarkerSize",3,"MarkerFaceColor",[0 0 1],"MarkerEdgeColor",[0.8 0.8 0.8])
plot(6*ones(31,1)+dx2,CI1_ih(:,2),'o',"MarkerSize",3,"MarkerFaceColor",[0 0 1],"MarkerEdgeColor",[0.8 0.8 0.8])
plot(7*ones(31,1)+dx3,CI2_ih(:,2),'o',"MarkerSize",3,"MarkerFaceColor",[0 0 1],"MarkerEdgeColor",[0.8 0.8 0.8])
for i = 1:31
    plot([6+dx2(i) 7+dx3(i)], [CI1_ih(i,2) CI2_ih(i,2)],"Color",[0.8 0.8 0.8])
end
plot([6,7], [nanmean(CI1_ih(:,2)) mean(CI2_ih(:,2))],'r'),
title('Interhemespheric Allegiance')
ylabel('Allegiance')
ylim([0.1 0.8])
box off
xline(4)

%% cross-modal allegiance
for i = 1:23
    NH_cm_ac(i) = mean(NH_talleg(T_chs(:,4) == 1 | T_chs(:,4) == 2 | T_chs(:,4) == 3 | T_chs(:,4) == 4, T_chs(:,4) == 5 ,condition,i),"all","omitmissing");
    NH_cm_ag(i) = mean( NH_talleg( T_chs(:,4) == 2 | T_chs(:,4) == 4 ,T_chs(:,4) == 5 ,condition,i),"all","omitmissing");
    NH_cm_stg(i) = mean(NH_talleg(T_chs(:,4) == 1 | T_chs(:,4) == 3,T_chs(:,4) == 5 ,condition,i),"all","omitmissing");
end
for i = 1:31
    if ismember(i,[2,11]) % exclude subject 2 & 13
        CI1_cm_ac(i) = nan;
        CI1_cm_ag(i) = nan;
        CI1_cm_stg(i) = nan;
        CI2_cm_ac(i) = nan;
        CI2_cm_ag(i) = nan;
        CI2_cm_stg(i) = nan;
    else
        CI1_cm_ac(i) = mean(CI1_talleg(T_chs(:,4) == 1 | T_chs(:,4) == 2 | T_chs(:,4) == 3 | T_chs(:,4) == 4, T_chs(:,4) == 5 ,condition,i),"all","omitmissing");
        CI1_cm_ag(i) = mean( CI1_talleg( T_chs(:,4) == 2 | T_chs(:,4) == 4 ,T_chs(:,4) == 5 ,condition,i),"all","omitmissing");
        CI1_cm_stg(i) = mean(CI1_talleg(T_chs(:,4) == 1 | T_chs(:,4) == 3,T_chs(:,4) == 5 ,condition,i),"all","omitmissing");
        CI2_cm_ac(i) = mean(CI2_talleg(T_chs(:,4) == 1 | T_chs(:,4) == 2 | T_chs(:,4) == 3 | T_chs(:,4) == 4, T_chs(:,4) == 5 ,condition,i),"all","omitmissing");
        CI2_cm_ag(i) = mean( CI2_talleg( T_chs(:,4) == 2 | T_chs(:,4) == 4 ,T_chs(:,4) == 5 ,condition,i),"all","omitmissing");
        CI2_cm_stg(i) = mean(CI2_talleg(T_chs(:,4) == 1 | T_chs(:,4) == 3,T_chs(:,4) == 5 ,condition,i),"all","omitmissing");
    end
end

[p,~] = signrank(NH_cm_ag, NH_cm_stg,'Tail','right')
[p,~] = signrank(CI1_cm_ag, CI1_cm_stg,'Tail','right')
[p,~] = signrank(CI2_cm_ag, CI2_cm_stg,'Tail','right')

% [p,h] = ranksum(NH_cm_ac, CI1_cm_ac,'Tail','right')
[p,~] = ranksum(NH_cm_ag, CI1_cm_ag)
[p,~] = ranksum(NH_cm_stg, CI1_cm_stg)

% [p,h] = ranksum(NH_cm_ac, CI2_cm_ac,'Tail','right')
[p,~] = ranksum(NH_cm_ag, CI2_cm_ag)
[p,~] = ranksum(NH_cm_stg, CI2_cm_stg)

% [p,h] = signrank(CI2_cm_ac, CI1_cm_ac,'Tail','right')
[p,~] = signrank(CI2_cm_ag, CI1_cm_ag)
[p,~] = signrank(CI2_cm_stg, CI1_cm_stg,'Tail','left')



CI1_cm = [CI1_cm_ag', CI1_cm_stg' ];
CI2_cm = [CI2_cm_ag', CI2_cm_stg' ];
NH_cm = NaN(size(CI1_cm));
NH_cm(1:23,:) = [NH_cm_ag', NH_cm_stg' ];
figure('Name','cross-modal allegiance')
boxplot([NH_cm(:,1), CI1_cm(:,1), CI2_cm(:,1), NaN(31,1), NH_cm(:,2), CI1_cm(:,2), CI2_cm(:,2)], 'Labels',{'NH','CI1', 'CI2','','NH','CI1', 'CI2'}), hold on
% dx1 = 0.12*randn(31,1);
% dx2 = 0.12*randn(31,1);
% dx3 = 0.12*randn(31,1);
plot(ones(31,1)+dx1,NH_cm(:,1),'o',"MarkerSize",3,"MarkerFaceColor",[0 0 1],"MarkerEdgeColor",[0.8 0.8 0.8])
plot(2*ones(31,1)+dx2,CI1_cm(:,1),'o',"MarkerSize",3,"MarkerFaceColor",[0 0 1],"MarkerEdgeColor",[0.8 0.8 0.8])
plot(3*ones(31,1)+dx3,CI2_cm(:,1),'o',"MarkerSize",3,"MarkerFaceColor",[0 0 1],"MarkerEdgeColor",[0.8 0.8 0.8])
for i = 1:31
    plot([2+dx2(i) 3+dx3(i)], [CI1_cm(i,1) CI2_cm(i,1)],"Color",[0.8 0.8 0.8])
end
plot([2,3], [mean(CI1_cm(:,1)) mean(CI2_cm(:,1))],'r')
plot(5*ones(31,1)+dx1,NH_cm(:,2),'o',"MarkerSize",3,"MarkerFaceColor",[0 0 1],"MarkerEdgeColor",[0.8 0.8 0.8])
plot(6*ones(31,1)+dx2,CI1_cm(:,2),'o',"MarkerSize",3,"MarkerFaceColor",[0 0 1],"MarkerEdgeColor",[0.8 0.8 0.8])
plot(7*ones(31,1)+dx3,CI2_cm(:,2),'o',"MarkerSize",3,"MarkerFaceColor",[0 0 1],"MarkerEdgeColor",[0.8 0.8 0.8])
for i = 1:31
    plot([6+dx2(i) 7+dx3(i)], [CI1_cm(i,2) CI2_cm(i,2)],"Color",[0.8 0.8 0.8])
end
plot([6,7], [mean(CI1_cm(:,2)) mean(CI2_cm(:,2))],'r')

title('Crossmodal Allegiance')
ylabel('Allegiance')
ylim([0.1 0.8])
box off
xline(4)


%% compare consensus partitions in different conditions
% par_distance_NH = parDistance_func(NH_talleg);
% par_distance_CI1 = parDistance_func(CI1_talleg);
% par_distance_CI2 = parDistance_func(CI2_talleg);
% 
% xx = for3Dboxplot(par_distance_NH(:,2:6,2:6));
% boxPlot3D(xx)
% xx = for3Dboxplot(par_distance_CI1(:,2:6,2:6));
% boxPlot3D(xx)
% xx = for3Dboxplot(par_distance_CI2(:,2:6,2:6));
% boxPlot3D(xx)
% 
% rep_ANOVA(par_distance_NH)
% rep_ANOVA(par_distance_CI1)
% rep_ANOVA(par_distance_CI2)
% 
% figure("Name","consensus partition distance-NH")
% plotParDistance(par_distance_NH)
% title("Compare consensus partitions: NH")
% figure("Name","consensus partition distance-CI1")
% plotParDistance(par_distance_CI1)
% title("Compare consensus partitions: CI1")
% figure("Name","consensus partition distance-CI2")
% plotParDistance(par_distance_CI2)
% title("Compare consensus partitions: CI2")

%% flexibility fluctuation
task_tran_pairsNH = [];
task_tran_aveNH = NaN(31,2);
for i = 1:size(aveTrialFlexNH,2)
    task_tran_pairsNH = cat(1, task_tran_pairsNH, flex_pairs(aveTrialFlexNH{i}));
    task_tran_aveNH(i,:) = mean(flex_pairs(aveTrialFlexNH{i}));
end
[~, p_flucNH] = ttest(task_tran_pairsNH(:,1), task_tran_pairsNH(:,2), "Tail","right")
% [~,p_flucAveNH] = ttest(task_tran_aveNH(:,1),task_tran_aveNH(:,2),'Tail','right')

task_tran_pairsCI1 = [];
for i = 1:size(aveTrialFlexCI1,2)
    task_tran_pairsCI1 = cat(1, task_tran_pairsCI1, flex_pairs(aveTrialFlexCI1{i}));
    task_tran_aveCI1(i,:) = mean(flex_pairs(aveTrialFlexCI1{i}));
end
[~, p_flucCI1] = ttest(task_tran_pairsCI1(:,1), task_tran_pairsCI1(:,2), "Tail","right")
%[~,p_flucAveCI1] = ttest(task_tran_aveCI1(:,1),task_tran_aveCI1(:,2),'Tail','right')


task_tran_pairsCI2 = [];
for i = 1:size(aveTrialFlexCI2,2)
    task_tran_pairsCI2 = cat(1, task_tran_pairsCI2, flex_pairs(aveTrialFlexCI2{i}));
    task_tran_aveCI2(i,:) = mean(flex_pairs(aveTrialFlexCI2{i}));
end
[~, p_flucCI2] = ttest(task_tran_pairsCI2(:,1), task_tran_pairsCI2(:,2), "Tail","right")
%[~,p_flucAveCI2] = ttest(task_tran_aveCI2(:,1),task_tran_aveCI2(:,2),'Tail','right')

fluc_data_ave(:, 1) = task_tran_aveNH(:,1) - task_tran_aveNH(:,2);
fluc_data_ave(:,2) = task_tran_aveCI1(:,1) - task_tran_aveCI1(:,2);
fluc_data_ave(:,3) = task_tran_aveCI2(:,1) - task_tran_aveCI2(:,2);

figure
boxplot(fluc_data_ave,'Labels',{'NH','CI1', 'CI2'}), hold on, yline(0)
plot(ones(31,1)+0.05*randn(31,1), fluc_data_ave(:,1),".b","MarkerSize",6)
plot(2*ones(31,1)+0.05*randn(31,1), fluc_data_ave(:,2),".b","MarkerSize",6)
plot(3*ones(31,1)+0.05*randn(31,1), fluc_data_ave(:,3),".b","MarkerSize",6)

%% Flexibility in time: group-level comparison
x = flex_t(1:600);
yNH = mean(flex_tsNH);
yCI1 = mean(flex_tsCI1);
yCI2 = mean(flex_tsCI2);
%comparing groups
lme = mixed_effects_modeling(flex_tsNH, flex_tsCI1, flex_tsCI2);
fixedEffectsTable = lme.Coefficients;
disp('Fixed effects coefficients (95% CIs):');
disp(fixedEffectsTable);
%effect of time
disp('Effect of time on individual groups:');
mixed_effects_modelingTime(flex_tsNH, flex_tsCI1, flex_tsCI2);
%%
win_len = 30;
for i = 1:600-win_len
    x_win(i) = x(i+ win_len/2) ;
    yNH_win(i) = mean( yNH(i:i+win_len) );
    yNH_error(i) = std( yNH(i:i+win_len) )/ sqrt(size(flex_tsNH,1));
    yCI1_win(i) = mean( yCI1(i:i+win_len) );
    yCI1_error(i) = std( yCI1(i:i+win_len) )/ sqrt(size(flex_tsCI1,1));
    yCI2_win(i) = mean( yCI2(i:i+win_len) );
    yCI2_error(i) = std( yCI2(i:i+win_len) )/ sqrt(size(flex_tsCI2,1));
end

figure("Name","flexibility in time (smoothed): group-level comparison")
shade(x_win, yNH_win - yNH_error,'w',x_win, yNH_win,'b', x_win, yNH_win + yNH_error,'w','FillType', [3,1], 'FillColor','b', 'FillAlpha',0.3), 
hold on, plot(x_win, yNH_win,'Color', "#0072BD")
shade(x_win, yCI1_win - yCI1_error,'w',x_win, yCI1_win, 'r', x_win, yCI1_win+ yCI1_error,'w','FillType', [3,1], 'FillColor','r', 'FillAlpha',0.3), 
plot(x_win, yCI1_win,'Color', "#D95319")
shade(x_win, yCI2_win - yCI2_error,'w',x_win, yCI2_win, 'g', x_win, yCI2_win+ yCI2_error,'w','FillType', [3,1], 'FillColor','g', 'FillAlpha',0.3), 
plot(x_win, yCI2_win,'Color', "#77AC30")
xlim([x_win(1), x_win(end)])
xlabel('time(seconds)')
ylabel('flexibility')

%% similarity in time: group-level comparison
%load("community_similarities.mat");
win_sim = 1:15; %number of neighbours
simInTimeNH = sim_in_time(NH_ComSim, win_sim);
simInTimeCI1 = sim_in_time(CI1_ComSim, win_sim);
simInTimeCI2 = sim_in_time(CI2_ComSim, win_sim);
win_len = 30;
for i = win_sim
    x = 1:size( simInTimeNH{i},2 );
    yNH = mean( simInTimeNH{i} );
    errNH = std( simInTimeNH{i} )/sqrt(size(simInTimeNH{i},1)); % Standard Error of the Mean
    yCI1 = mean( simInTimeCI1{i} );
    errCI1 = std( simInTimeCI1{i} )/sqrt(size(simInTimeCI1{i},1)); % Standard Error of the Mean
    yCI2 = mean( simInTimeCI2{i} );
    errCI2 = std( simInTimeCI2{i} )/sqrt(size(simInTimeCI2{i},1)); % Standard Error of the Mean
    %plot
    if i == 1
        figure("Name",num2str(i))
        x_ticks = linspace(x_win(1), x_win(end),length(x));
        errorbar( x_ticks, yNH, errNH,'Color',"#0072BD" ), 
        hold on, plot(x_ticks, yNH, "LineWidth",2,'Color',"#0072BD")
        errorbar( x_ticks, yCI1, errCI1, 'Color', "#D95319" ), 
        plot(x_ticks, yCI1, "LineWidth",2,'Color', "#D95319")
        errorbar( x_ticks, yCI2, errCI2, 'Color', "#77AC30" ), 
        plot(x_ticks, yCI2, "LineWidth",2,'Color', "#77AC30")
        xlim([min(x_ticks),max(x_ticks)])
        xlabel('time(seconds)')
        ylabel('community similarity (rand z-score)')
        title("Group-level comparison on community similarity for adjacent windows")
        legend({'NH', 'CI1','CI2'})
        % figure, plot(x_ticks, yNH, "LineWidth",2,'Color',"#0072BD"), hold on, plot(x_ticks, yNH, 'Color',"#0072BD","Marker","x","MarkerSize",8, "MarkerEdgeColor","#0072BD")
        lme = mixed_effects_modeling(simInTimeNH{i}, simInTimeCI1{i}, simInTimeCI2{i});
        fixedEffectsTable = lme.Coefficients;
        disp('Fixed effects coefficients (95% CIs):');
        disp(fixedEffectsTable);
        %effect of time
        disp('Effect of time on individual groups:');
        mixed_effects_modelingTime(simInTimeNH{i}, simInTimeCI1{i}, simInTimeCI2{i});

    end   
end

x = flex_t(1:599);
yNH = mean(simInTimeNH{1});
yCI1 = mean(simInTimeCI1{1});
yCI2 = mean(simInTimeCI2{1});
for i = 1:599-win_len
    x_winTSM(i) = x(i+ win_len/2) ;
    yNH_winTSM(i) = mean( yNH(i:i+win_len) );
    yNH_errorTSM(i) = std( yNH(i:i+win_len) )/ sqrt(size(flex_tsNH,1));
    yCI1_winTSM(i) = mean( yCI1(i:i+win_len) );
    yCI1_errorTSM(i) = std( yCI1(i:i+win_len) )/ sqrt(size(flex_tsCI1,1));
    yCI2_winTSM(i) = mean( yCI2(i:i+win_len) );
    yCI2_errorTSM(i) = std( yCI2(i:i+win_len) )/ sqrt(size(flex_tsCI2,1));
end
figure("Name","community similarity of adjacent windows (smoothed): group-level comparison")
shade(x_winTSM, yNH_winTSM - yNH_errorTSM,'w',x_winTSM, yNH_winTSM,'b', x_winTSM, yNH_winTSM + yNH_errorTSM,'w','FillType', [3,1], 'FillColor','b', 'FillAlpha',0.3), 
hold on, plot(x_winTSM, yNH_winTSM,'Color', "#0072BD")
shade(x_winTSM, yCI1_winTSM - yCI1_errorTSM,'w',x_winTSM, yCI1_winTSM, 'r', x_winTSM, yCI1_winTSM+ yCI1_errorTSM,'w','FillType', [3,1], 'FillColor','r', 'FillAlpha',0.3), 
plot(x_winTSM, yCI1_winTSM,'Color', "#D95319")
shade(x_winTSM, yCI2_winTSM - yCI2_errorTSM,'w',x_winTSM, yCI2_winTSM, 'g', x_winTSM, yCI2_winTSM+ yCI2_errorTSM,'w','FillType', [3,1], 'FillColor','g', 'FillAlpha',0.3), 
plot(x_winTSM, yCI2_winTSM,'Color', "#77AC30")
xlim([x_winTSM(1), x_winTSM(end)])
xlabel('time(seconds)')
ylabel('community similarity')



%% correlation between time-related flexibility and community similarity for adjacent windows
%NH
fNH = [];
sNH = [];
for i = 1:size(flex_tsNH,1)
    fNH = [fNH, flex_tsNH(i,2:600)];
    sNH = [sNH, simInTimeNH{1}(i,1:599)];
end
figure
subplot(3,1,1)
lm = fitlm(fNH, sNH);
plot(lm)
text(0.1, 0, ['p = ', num2str(sqrt(lm.ModelFitVsNullModel.Pvalue)),', R = ', num2str(sqrt(lm.Rsquared.Ordinary))])
title('NH')
xlabel('flexibility')
ylabel('community similarity')

%CI1
fCI1 = [];
sCI1 = [];
for i = 1:size(flex_tsCI1,1)
    fCI1 = [fCI1, flex_tsCI1(i,2:600)];
    sCI1 = [sCI1, simInTimeCI1{1}(i,1:599)];
end
subplot(3,1,2)
lm = fitlm(fCI1, sCI1);
plot(lm)
text(0.1, 0, ['p = ', num2str(sqrt(lm.ModelFitVsNullModel.Pvalue)),', R = ', num2str(sqrt(lm.Rsquared.Ordinary))])
title('CI1')
xlabel('flexibility')
ylabel('community similarity')

%CI1
fCI2 = [];
sCI2 = [];
for i = 1:size(flex_tsCI2,1)
    fCI2 = [fCI2, flex_tsCI2(i,2:600)];
    sCI2 = [sCI2, simInTimeCI2{1}(i,1:599)];
end
subplot(3,1,3)
lm = fitlm(fCI2, sCI2);
plot(lm)
text(0.1, 0, ['p = ', num2str(sqrt(lm.ModelFitVsNullModel.Pvalue)),', R = ', num2str(sqrt(lm.Rsquared.Ordinary))])
title('CI2')
xlabel('flexibility')
ylabel('community similarity')


%% functions
function [F_opt44, F_tempNull44, C_opt44, C_nodNull44, Sopt44 ] = resizeTo44(F_opt, F_tempNull, C_opt, C_nodNull, Sopt)
[F_opt44, F_tempNull44] = arrey44(F_opt, F_tempNull);
if size(F_opt,1)<31
    F_opt44 = NaN(44,1);
end
C_opt44 = matrix44(C_opt, F_opt(:,2));
C_nodNull44 = [];    
Sopt44 = NaN(44, size(Sopt,2), size(Sopt,3));
for i = 1:size(Sopt,3)
    Sopt44(F_opt(:,2),:,i) = Sopt(:,:,i);
end
end
%__________________________________________________________________________
function [F_opt44, F_tempNull44] = arrey44(F_opt, F_tempNull)
F_opt44 = NaN(44,1);
F_tempNull44 = NaN(44,50);
F_opt44(F_opt(:,2),1) = F_opt(:,1);
end
%__________________________________________________________________________
function C44 = matrix44(C,chs)
C44 = NaN(44,44);
for m = 1:44
    for n = 1:44
        if ismember(m,chs) && ismember(n,chs)
            C44(m,n) = C(find(chs==m), find(chs==n));
            if m == n
                C44(m,n) = 1;
            end
        end
    end
end
end
%__________________________________________________________________________
function scores_converted = rau_convert(scores)
for i = 1:size(scores,1)
    for j = 1: size(scores,2)
        scores_converted(i,j) = (146/pi)*(asin(sqrt(scores(i,j)/101)) + asin(sqrt((scores(i,j)+1)/101)))-23;
    end
end
end
%__________________________________________________________________________
function [flex_ts, flex_ts44, flex_t, flex_l] = flexInTime(Sopt, chs, win_time, win_len, win_step, onsets)
Fs = 3.9125;
rise_time = 0;
fall_time = 0;
flex_ts = NaN(size(Sopt,1), min(size(Sopt,2),length(win_time))-1);
flex_t = NaN(1,min(size(Sopt,2),length(win_time))-1);
flex_l = NaN(1,min(size(Sopt,2),length(win_time))-1);
t = 0:1/Fs:(max(win_time)+win_len/2);
trial_waves = zeros(3,length(t));
for i = 1:3
    for j = 1:18
        if onsets(i,j)~=0
            start_sample = find(abs(t - onsets(i,j)-rise_time) == min(abs(t - onsets(i,j)-rise_time)));
            end_sample = find(abs(t - onsets(i,j)-12.5-fall_time) == min(abs(t - onsets(i,j)-12.5-fall_time)));
            trial_waves(i,start_sample:end_sample) = 1;
        end
    end
end
trial_length  = end_sample - start_sample+1;
for i = 1:length(win_time)-1
    win_wave = zeros(1,length(t));
    start_edge = win_time(i) - win_len/2;
    end_edge = win_time(i+1) + win_len/2;
    start_edge_sample = find(abs(t-start_edge) == min(abs(t-start_edge)));
    end_edge_sample = find(abs(t-end_edge) == min(abs(t-end_edge)));
    win_wave(start_edge_sample:end_edge_sample) = 1;
    percent(i,1) = 100*sum(trial_waves(1,:).*win_wave)/trial_length;
    percent(i,2) = 100*sum(trial_waves(2,:).*win_wave)/trial_length;
    percent(i,3) = 100*sum(trial_waves(3,:).*win_wave)/trial_length;
end
for n = 1:size(flex_ts,2)
    for m = 1:size(flex_ts,1)
        flex_ts(m,n) = length(find(Sopt(m,n,:)-Sopt(m,n+1,:) ~= 0))/(size(Sopt,3));
    end
    flex_t(n) = 0.5*(win_time(n) + win_time(n+1));
end
flex_ts44 = NaN(44,size(flex_ts,2));
flex_ts44(chs,:) = flex_ts;

flex_l = NaN(1,min(size(Sopt,2),length(win_time))-1);
for i = 1:length(flex_ts)
    max_percent = find(percent(i,:)== max(percent(i,:)));
    switch max_percent(1)
        case 1            
            if percent(i,max_percent)>90    %audio
                flex_l(i) = 1;
            elseif percent(i,max_percent)<=5    %isi
                flex_l(i) = 4;
            else
                flex_l(i) = 5;      %trans
            end
        case 2
            if percent(i,max_percent)>90
                flex_l(i) = 2;
            elseif percent(i,max_percent)<=5
                flex_l(i) = 4;
            else
                flex_l(i) = 5;
            end
        case 3
                flex_l(i) = 3;
    end
end
end
%__________________________________________________________________________     
function plot_flexInTime(flex_ts, flex_t, flex_l, onsets, name_record)
figure
% stem(flex_t, mean(flex_ts),'LineStyle','-','Color',"#0072BD")
plot(flex_t, mean(flex_ts))
hold on
for i = 1:length(flex_t)
    if flex_l(i)==1
        plot(flex_t(i), mean(flex_ts(:,i)),"Color",[0.3010 0.7450 0.9330],"Marker","o","MarkerFaceColor",[0.3010 0.7450 0.9330], "MarkerSize",12 )     %blue for audio
    elseif flex_l(i)==2
        plot(flex_t(i), mean(flex_ts(:,i)),"Color",[0.8500 0.3250 0.0980],"Marker","o","MarkerFaceColor",[0.8500 0.3250 0.0980], "MarkerSize",12 )     %red for visual
    elseif flex_l(i)==3
        plot(flex_t(i), mean(flex_ts(:,i)),"Color",[0.6 0.6 0.6],"Marker","o","MarkerFaceColor",[0.6 0.6 0.6] , "MarkerSize",12 )     %green for ctl
    elseif flex_l(i)==4
        plot(flex_t(i), mean(flex_ts(:,i)),"Color",[0.9290 0.6940 0.1250],"Marker","o","MarkerFaceColor",[0.9290 0.6940 0.1250], "MarkerSize",12 )     %cyan for isi
    elseif flex_l(i) == 5
        plot(flex_t(i), mean(flex_ts(:,i)),"Color",'k',"Marker","o","MarkerFaceColor",'w', "MarkerSize",12  )       %yellow for transition
    end 
end
title(name_record)
% yline(mean(flex_ts(:,find(flex_l==1)),"all"),"Color",'b','LineStyle','--')
% yline(mean(flex_ts(:,find(flex_l==2)),"all"),"Color",'r','LineStyle','--')
% yline(mean(flex_ts(:,find(flex_l==3)),"all"),"Color",'g','LineStyle','--')
% yline(mean(flex_ts(:,find(flex_l==4)),"all"),"Color",'c','LineStyle','--')
% yline(mean(flex_ts(:,find(flex_l==5)),"all"),"Color",'y','LineStyle','--')

trialColors = [0.3010 0.7450 0.9330; [0.8500 0.3250 0.0980]; [0.6 0.6 0.6]];
for i = 1:3
    for j = 1:18
        if i==3 && j>10
            continue
        else
            xline(onsets(i,j),"Color",trialColors(i,:),'LineStyle','-.',"LineWidth",2)
        end
    end
end
xlabel('Time(Second)')
ylabel('Average Flexibility')
end
%__________________________________________________________________________     
function [audAve, visAve, ctlAve, isiAve, trnAve] = flexAveTask(flex_ts,flex_l)
audAve = mean(flex_ts(:,find(flex_l==1)),"all");
visAve = mean(flex_ts(:,find(flex_l==2)),"all");
ctlAve = mean(flex_ts(:,find(flex_l==3)),"all");
isiAve = mean(flex_ts(:,find(flex_l==4)),"all");
trnAve = mean(flex_ts(:,find(flex_l==5)),"all");
end
%__________________________________________________________________________     
function G_ROIaveFlex = ROIaveFlex(G_tflex, ROI_chs)
for i = 1:size(G_tflex, 2)
    for j = 1:6
        G_ROIaveFlex(i,j) = mean(G_tflex(find(ROI_chs==j),i),"omitmissing");
    end
end
end
%__________________________________________________________________________     
function [allChs_flex, LSTG_flex, LAG_flex, RSTG_flex, RAG_flex, OC_flex, IFG_flex] = ROIaveFlexISI(flex_ts,flex_l,roi_allchs, good_chs)
roi_chs = roi_allchs(good_chs);
allChs_flex = mean( flex_ts(:,find(flex_l==4)), "all" );
LSTG_flex = mean( flex_ts(find(roi_chs == 1), find(flex_l==4)), "all" );
LAG_flex = mean( flex_ts(find(roi_chs == 2), find(flex_l==4)), "all" );
RSTG_flex = mean( flex_ts(find(roi_chs == 3), find(flex_l==4)), "all" );
RAG_flex = mean( flex_ts(find(roi_chs == 4), find(flex_l==4)), "all" );
OC_flex = mean( flex_ts(find(roi_chs == 5), find(flex_l==4)), "all" );
IFG_flex = mean( flex_ts(find(roi_chs == 6), find(flex_l==4)), "all" );
end
%__________________________________________________________________________
function allegs44 = alleg_cons(Sopt, flex_l, chs)
allegs(:,:,1) = allegiance(Sopt);
allegs44(:,:,1) = matrix44(allegs(:,:,1),chs);
for i = 2:6
    Sopt_con = Sopt(:, find(flex_l==i-1), :);
    allegs(:,:,i) = allegiance(Sopt_con);
    allegs44(:,:,i) = matrix44(allegs(:,:,i),chs);
end
end
%__________________________________________________________________________
function mat_allegiance = allegiance(S)
for k = 1:size(S,3)
    for m = 1:size(S,1)
        for n = 1:size(S,1)
            mat3_allegiance(m,n,k) = length(find(S(n,:,k)==S(m,:,k)))/size(S,2);
        end
    end
end
mat_allegiance = mean(mat3_allegiance,3,"omitmissing");
end
%__________________________________________________________________________
function subParDis = sub_partition_distance( NHmean_talleg, CI1sub_talleg, CI2sub_talleg )
badChs1 = find( isnan(CI1sub_talleg(:,1,1))  ); 
badChs2 = find( isnan(CI2sub_talleg(:,1,1))  );
NHmean_talleg1 = removeBadChs(NHmean_talleg, badChs1);
NHmean_talleg2 = removeBadChs(NHmean_talleg, badChs2);
CI1sub_talleg = removeBadChs(CI1sub_talleg, badChs1);
CI2sub_talleg = removeBadChs(CI2sub_talleg, badChs2);
for i = 1:6
    [NH1_Ci(:,i) NH1_Q(i)] = modularity_und( NHmean_talleg1(:,:,i), 1.1 );
    [NH2_Ci(:,i) NH2_Q(i)] = modularity_und( NHmean_talleg2(:,:,i), 1.1 );
    [CI1_Ci(:,i) CI1_Q(i)] = modularity_und( CI1sub_talleg(:,:,i), 1.1 );
    [CI2_Ci(:,i) CI2_Q(i)] = modularity_und( CI2sub_talleg(:,:,i), 1.1 );
end
for i = 1:6
    [VIn, MIn] = partition_distance([NH1_Ci(:,i) CI1_Ci(:,i)]);
    subParDis(1,i) = VIn(1,2);
    [VIn, MIn] = partition_distance([NH2_Ci(:,i) CI2_Ci(:,i)]);
    subParDis(2,i) = VIn(1,2);
end
end
%__________________________________________________________________________
function mat = removeBadChs(mat44, badChs)
for i = 1:6
    m = mat44(:,:,i);
    m(badChs,:) = [];
    m(:,badChs) = [];
    mat(:,:,i) = m;
end
end
%__________________________________________________________________________
function plot_modules(Cis,chsPos)
figure
module_colors = {'b', 'r', 'g', 'y', 'k','c'};
titles = {'NH', 'CI1', 'CI2'};
for m = 1:3
    subplot(3,1,m)
    for n = 1:44
        plot3(chsPos(n,1), chsPos(n,2), chsPos(n,3), 'ok', 'MarkerFaceColor',module_colors{Cis(n,m)})
        hold on
    end
    title(titles{m})
end
end
%__________________________________________________________________________
function plot_chSignificance(R,chsPos)
for n = 1:44
    plot3(chsPos(n,1), chsPos(n,2), chsPos(n,3), 'ok', 'MarkerSize',2+15*R(n))
    hold on
end
end
%__________________________________________________________________________
function plot_allegiance(alleg, module_chs, module_labels)
module_chs{1} = reshape(module_chs{1},[1,length(module_chs{1})]);
module_chs{2} = reshape(module_chs{2},[1,length(module_chs{2})]);
module_chs{3} = reshape(module_chs{3},[1,length(module_chs{3})]);
module_chs{4} = reshape(module_chs{4},[1,length(module_chs{4})]);
imagesc(alleg( [module_chs{1},module_chs{2},module_chs{3},module_chs{4}],  [module_chs{1},module_chs{2},module_chs{3},module_chs{4}]),[0 1]), colorbar
xline(length(module_chs{1})+0.5,'w-.'), 
xline(length(module_chs{1})+length(module_chs{2})+0.5,'w-.'), 
xline(length(module_chs{1})+length(module_chs{2})+length(module_chs{3})+0.5,'w-.'), 
yline(length(module_chs{1})+0.5,'w-.'), 
yline(length(module_chs{1})+length(module_chs{2})+0.5,'w-.'), 
yline(length(module_chs{1})+length(module_chs{2})+length(module_chs{3})+0.5,'w-.'), 
xticks([length(module_chs{1})/2, ...
    length(module_chs{1})+length(module_chs{2})/2,...
    length(module_chs{1})+length(module_chs{2})+length(module_chs{3})/2,...
    length(module_chs{1})+length(module_chs{2})+length(module_chs{3})+length(module_chs{4})/2]),
    xticklabels(module_labels)
yticks([length(module_chs{1})/2, ...
    length(module_chs{1})+length(module_chs{2})/2,...
    length(module_chs{1})+length(module_chs{2})+length(module_chs{3})/2,...
    length(module_chs{1})+length(module_chs{2})+length(module_chs{3})+length(module_chs{4})/2]),
yticklabels(module_labels)
end
%__________________________________________________________________________
function integration = module_integration( talleg, Ci )
modules = unique(Ci);
for i = 1:size(talleg,4)
    sub_alleg = zeros(44,44);
    Ci_sub = Ci;
    sub_alleg(:,:) = talleg(:,:,1,i);
    %remove bad chs
    badChs = find( isnan(sub_alleg(:,1))  ); 
    if length(badChs)~=0
        sub_alleg(badChs,:) = [];
        sub_alleg(:,badChs) = [];
    end
    Ci_sub(badChs) = [];
    %re-order the channels
    chs_m1 = find(Ci_sub==2)';
    chs_m2 = find(Ci_sub==3)';
    chs_m3 = find(Ci_sub==4)';
    chs_m4 = find(Ci_sub==1)';
    sub_alleg = sub_alleg([chs_m1, chs_m2, chs_m3, chs_m4], [chs_m1, chs_m2, chs_m3, chs_m4]);
    integration(i,:) = [mean( sub_alleg(1:length(chs_m1), length(chs_m1)+1:end),"all"),...
        mean(sub_alleg(length(chs_m1)+1:length(chs_m1)+length(chs_m2), [1:length(chs_m1) length(chs_m1)+length(chs_m2)+1:end ]),"all"),...
        mean(sub_alleg(length(chs_m1)+length(chs_m2)+1:length(chs_m1)+length(chs_m2)+length(chs_m3), [1:length(chs_m1)+length(chs_m2) length(chs_m1)+length(chs_m2)+length(chs_m3)+1:end ]),"all"),...
        mean(sub_alleg(length(chs_m1)+length(chs_m2)+length(chs_m3)+1:end, 1:length(chs_m1)+length(chs_m2)+length(chs_m3)),"all")];
end
end
%__________________________________________________________________________
function par_distance = parDistance_func(talleg)
tau = 0;
reps = 30;
for m = 1:size(talleg,4)
    badChs = find( isnan(talleg(:,1,1,m))  ); 
    talleg_sub = removeBadChs(talleg(:,:,:,m), badChs);
    for n = 1:6
        Ci(:,n) = consensus_und( talleg_sub(:,:,n), tau, reps );        
    end
    for n = 1:6
        for k = 1:6
            par_distance(m,n,k) = partition_distance(Ci(:,n),Ci(:,k));
        end        
    end
    clear Ci
end
end
%__________________________________________________________________________
function rep_ANOVA(par_distance)
num_subs = size(par_distance,1);
data = [reshape(par_distance(:,2,3),[num_subs,1]) reshape(par_distance(:,2,4),[num_subs,1]) reshape(par_distance(:,2,5),[num_subs,1])...
    reshape(par_distance(:,3,4),[num_subs,1]) reshape(par_distance(:,3,5),[num_subs,1]) reshape(par_distance(:,4,5),[num_subs,1]) ...
    reshape(par_distance(:,2,6),[num_subs,1]) reshape(par_distance(:,3,6),[num_subs,1]) reshape(par_distance(:,4,6),[num_subs,1]) reshape(par_distance(:,5,6),[num_subs,1])];
T = array2table(data, 'VariableNames', {'dis1', 'dis2', 'dis3', 'dis4', 'dis5', 'dis6', 'dis7', 'dis8', 'dis9', 'dis10'});
rm = fitrm(T, 'dis1-dis10~1', 'WithinDesign', {'Time1', 'Time2', 'Time3', 'Time4', 'Time5', 'Time6', 'Time7', 'Time8', 'Time9', 'Time10'});
% Perform the repeated measures ANOVA
ranovatbl = ranova(rm);
% Display the results
disp(ranovatbl)
end
%__________________________________________________________________________
function xx = for3Dboxplot(par_distance)
xx = par_distance;
for m = 1:size(par_distance,3)
    for n = 1:size(par_distance,2)
        if m >= n
           xx(:,m,n) = nan;
        end
    end
end
end
%__________________________________________________________________________
function plotParDistance(par_distance)
num_subs = size(par_distance,1);
Xboxplot = {[reshape(par_distance(:,2,2),[num_subs,1]) , reshape(par_distance(:,3,2),[num_subs,1]), reshape(par_distance(:,4,2),[num_subs,1]) , reshape(par_distance(:,5,2),[num_subs,1]) ],...
    [reshape(par_distance(:,2,3),[num_subs,1]) , reshape(par_distance(:,3,3),[num_subs,1]), reshape(par_distance(:,4,3),[num_subs,1]) , reshape(par_distance(:,5,3),[num_subs,1]) ],...
    [reshape(par_distance(:,2,4),[num_subs,1]) , reshape(par_distance(:,3,4),[num_subs,1]), reshape(par_distance(:,4,4),[num_subs,1]) , reshape(par_distance(:,5,4),[num_subs,1]) ],...
    [reshape(par_distance(:,2,5),[num_subs,1]) , reshape(par_distance(:,3,5),[num_subs,1]), reshape(par_distance(:,4,5),[num_subs,1]) , reshape(par_distance(:,5,5),[num_subs,1]) ],...
    [reshape(par_distance(:,2,6),[num_subs,1]) , reshape(par_distance(:,3,6),[num_subs,1]), reshape(par_distance(:,4,6),[num_subs,1]) , reshape(par_distance(:,5,6),[num_subs,1]) ]};
boxplotGroup(Xboxplot,'Colors','kkkkr','GroupType','betweenGroups',...
    'PrimaryLabels', {'aud', 'vis', 'ctl', 'isi', 'trs'}, 'SecondaryLabels',{'aud', 'vis', 'ctl', 'isi'}, 'Notch','off')
ylabel('distance')
ylim([-0.05 0.6])
xline(6,'b'), xline(12,'b'), xline(18,'b'),
end
%__________________________________________________________________________
function ComSim = part_similarityMat(Sopt, flex_l)
ComSim = zeros(size(Sopt,2)-1, size(Sopt,2)-1);
for i = 1:size(Sopt,3)
    ComSim = ComSim + part_similarityItr(Sopt(:,2:end,i));
end
ComSim = ComSim/size(Sopt,3);
ComSim = cat(1,flex_l,ComSim);
end
%__________________________________________________________________________
function ComSim_itr = part_similarityItr(Sopt_itr)
for m = 1:size(Sopt_itr,2)
    for n = 1:size(Sopt_itr,2)
        if m<=n
            [ComSim_itr(m,n),~,~,~] = zrand(Sopt_itr(:,m),Sopt_itr(:,n));
            ComSim_itr(n,m) = ComSim_itr(m,n);
        end
    end
end
size(ComSim_itr), pause
end
%__________________________________________________________________________
function [comSim_mat, comSim_inTrials] = comSim_AVEs(ComSim)
for i = 1:numel(ComSim)
    SimMat = ComSim{1,i}(2:end, :);
    SimMat(eye(size(SimMat))==1) = nan; %remove elements on the main diagonal (high z-rand scores for comparing communities to themselves)
    win_lables = ComSim{1,i}(1,:); 
    for m = 1:4
        for n = 1:4
            if n~=m
                comSim_mat(i,m,n) = mean(SimMat(win_lables==m, win_lables==n),"all", "omitmissing");   
            else
                locs = find(win_lables==n);
                SimMat_task = SimMat(locs, locs);
                ave_sim = 0;
                cnt_sim = 0;
                comSim_inTrials(i,n) = 0;
                cnt_intrial = 0;
                for j = 1:size(SimMat_task,1)
                    for k = 1:size(SimMat_task,1)
                        if abs(j-k)>5
                            ave_sim = ave_sim + SimMat_task(j,k);
                            cnt_sim = cnt_sim + 1;
                        elseif abs(j-k)<7 && abs(j-k)~=0
                            comSim_inTrials(i,n) = comSim_inTrials(i,n) + SimMat_task(j,k);
                            cnt_intrial = cnt_intrial + 1;
                        end
                    end
                end
                comSim_mat(i,m,n) = ave_sim/cnt_sim;
                comSim_inTrials(i,n) = comSim_inTrials(i,n)/cnt_intrial;
            end            
        end
    end
end 
end
%__________________________________________________________________________
function p = comparing_communities(comSim_mat)
for i = 1:4
    for j = 1:4
        if i~=j
            if i<j
                [~,pp] = ttest(comSim_mat(:,i,i),comSim_mat(:,i,j),"Tail", "right");
            else
                [~,pp] = ttest(comSim_mat(:,i,i),comSim_mat(:,j,i),"Tail", "right");
            end
            p(i,j) = pp*12;
        else
            p(i,j) = nan;
        end
    end
end
end
%__________________________________________________________________________
function [comSim_mat, comSim_aveTask, comSim_aveNonTask] = comSimPerm_AVEs(ComSim)
for i = 1:numel(ComSim)
    SimMat = ComSim{1,i}(2:end, :);
    SimMat(eye(size(SimMat))==1) = nan; %remove elements on the main diagonal (high z-rand scores for comparing communities to themselves)
    win_lables = ComSim{1,i}(1,:);
    perm_lables = win_lables(randperm(length(win_lables)));
    comSim_aveTask(i) = mean(SimMat(perm_lables==1 | perm_lables==2 , perm_lables==1 | perm_lables==2 ),"all", "omitmissing");
    comSim_aveNonTask(i) = mean(SimMat(perm_lables==3 | perm_lables==4, perm_lables==3 | perm_lables==4),"all", "omitmissing");
    for m = 1:4
        for n = 1:4
            if n>=m
                comSim_mat(i,m,n) = mean(SimMat(perm_lables==m, perm_lables==n),"all", "omitmissing");   
            else
                comSim_mat(i,m,n) = nan;
            end            
        end
    end
end 
end
%__________________________________________________________________________
function p_perm = perm_test(comSimPerm_mat, comSim_mat, con1, con2)
if strcmp(con1,'aud')  
   x1 = 1; y1 = 1;
   switch con2
       case 'vis'
           x2 = 1; y2 = 2;
       case 'ctl'
           x2 = 1; y2 = 3;
       case 'isi'
           x2 = 1; y2 = 4;
   end
elseif strcmp(con1, 'vis')
   x1 = 2; y1 = 2;
   switch con2
       case 'aud'
           x2 = 1; y2 = 2;
       case 'ctl'
           x2 = 2; y2 = 3;
       case 'isi'
           x2 = 2; y2 = 4;
   end
elseif strcmp(con1,'ctl')
   x1 = 3; y1 = 3;
   switch con2
       case 'aud'
           x2 = 1; y2 = 3;
       case 'vis'
           x2 = 2; y2 = 3;
       case 'isi'
           x2 = 3; y2 = 4;
   end
elseif strcmp(con1,'isi')
   x1 = 4; y1 = 4;
   switch con2
       case 'aud'
           x2 = 1; y2 = 4;
       case 'vis'
           x2 = 2; y2 = 4;
       case 'ctl'
           x2 = 3; y2 = 4;
   end
end

for i = 1:size(comSimPerm_mat,1)
    ave_diff = 0;
    for j = 1:size(comSimPerm_mat,2)
        ave_diff = ave_diff + comSimPerm_mat(i,j,x1,y1) - comSimPerm_mat(i,j,x2,y2);
    end
    perm_val(i) = ave_diff/size(comSimPerm_mat,2);
end

real_val = 0;
for j = 1:size(comSim_mat,1)
    real_val = real_val + comSim_mat(j,x1,y1) - comSim_mat(j,x2,y2);
end
real_val = real_val/size(comSim_mat,1);

p_perm = length(find(perm_val > real_val))/1000;

p_perm = p_perm*12;

figure
histfit(perm_val,50), hold on
xline(real_val,'r')
end
%__________________________________________________________________________
function boxplot_comSims(comSim_mat)
labels = {'aud','vis','ctl','isi'};
data = [];
ctl = 0;
for i = 1:4
    for j = 1:4
        if i<=j
            ctl = ctl+1;
            data = cat(2,data, comSim_mat(:,i,j));
            tick_labels{ctl} = cat(2,labels{i},'-',labels{j}); 
        end
    end
end
data = data(:,[1,5,8,10,2,3,4,6,7,9]);
tick_labels = tick_labels([1,5,8,10,2,3,4,6,7,9]);
figure
boxplot(data)
xticklabels(tick_labels)
hold on
xline(4.5)
for i = 1:size(comSim_mat,1)
    hold on
    plot([1:10], data(i,:),"Color",[0.7,0.7,0.7])
    % plot([1:10], data(i,:),"Color",[0.7,0.7,0.7],"MarkerFaceColor",["auto"],'Marker','o')
end
end

%__________________________________________________________________________
function boxplot_comSims2(comSim_mat)
data = [];
ctl = 0;
for i = 1:4
    for j = 1:4
        if i<=j
            ctl = ctl+1;
            data = cat(2,data, comSim_mat(:,i,j));
        end
    end
end
data = data(:,[1,5,8,10,2,3,4,6,7,9]);
figure
subplot(2,2,1)
boxplot([data(:,1)-data(:,5) data(:,1)-data(:,6) data(:,1)-data(:,7)]),hold on
plot(ones(size(comSim_mat,1),1) + 0.1*randn(size(comSim_mat,1),1), data(:,1)-data(:,5),'o',"MarkerFaceColor",[0.7,0.7,0.7])
plot(2*ones(size(comSim_mat,1),1) + 0.1*randn(size(comSim_mat,1),1), data(:,1)-data(:,6),'o',"MarkerFaceColor",[0.7,0.7,0.7])
plot(3*ones(size(comSim_mat,1),1) + 0.1*randn(size(comSim_mat,1),1), data(:,1)-data(:,7),'o',"MarkerFaceColor",[0.7,0.7,0.7])
[~,p1] = ttest(data(:,1)-data(:,5),0,"Tail","right");
[~,p2] = ttest(data(:,1)-data(:,6),0,"Tail","right");
[~,p3] = ttest(data(:,1)-data(:,7),0,"Tail","right");
text(0.6,-0.3,num2str(12*p1))
text(1.6,-0.3,num2str(12*p2))
text(2.6,-0.3,num2str(12*p3))
title('AUD')
yline(0)
xticklabels({'vis','ctl','isi'})
subplot(2,2,2)
boxplot([data(:,2)-data(:,5) data(:,2)-data(:,8) data(:,2)-data(:,9)]), hold on
plot(ones(size(comSim_mat,1),1) + 0.1*randn(size(comSim_mat,1),1), data(:,2)-data(:,5),'o',"MarkerFaceColor",[0.7,0.7,0.7])
plot(2*ones(size(comSim_mat,1),1) + 0.1*randn(size(comSim_mat,1),1), data(:,2)-data(:,8),'o',"MarkerFaceColor",[0.7,0.7,0.7])
plot(3*ones(size(comSim_mat,1),1) + 0.1*randn(size(comSim_mat,1),1), data(:,2)-data(:,9),'o',"MarkerFaceColor",[0.7,0.7,0.7])
[~,p1] = ttest(data(:,2)-data(:,5),0,"Tail","right");
[~,p2] = ttest(data(:,2)-data(:,8),0,"Tail","right");
[~,p3] = ttest(data(:,2)-data(:,9),0,"Tail","right");
text(0.6,-0.3,num2str(12*p1))
text(1.6,-0.3,num2str(12*p2))
text(2.6,-0.3,num2str(12*p3))
title('VIS')
yline(0)
xticklabels({'aud','ctl','isi'})
subplot(2,2,3)
boxplot([data(:,3)-data(:,6) data(:,3)-data(:,8) data(:,3)-data(:,10)]), hold on
plot(ones(size(comSim_mat,1),1) + 0.1*randn(size(comSim_mat,1),1), data(:,3)-data(:,6),'o',"MarkerFaceColor",[0.7,0.7,0.7])
plot(2*ones(size(comSim_mat,1),1) + 0.1*randn(size(comSim_mat,1),1), data(:,3)-data(:,8),'o',"MarkerFaceColor",[0.7,0.7,0.7])
plot(3*ones(size(comSim_mat,1),1) + 0.1*randn(size(comSim_mat,1),1), data(:,3)-data(:,10),'o',"MarkerFaceColor",[0.7,0.7,0.7])
[~,p1] = ttest(data(:,3)-data(:,6),0,"Tail","right");
[~,p2] = ttest(data(:,3)-data(:,8),0,"Tail","right");
[~,p3] = ttest(data(:,3)-data(:,10),0,"Tail","right");
text(0.6,-0.3,num2str(12*p1))
text(1.6,-0.3,num2str(12*p2))
text(2.6,-0.3,num2str(12*p3))
title('CTL')
yline(0)
xticklabels({'aud','vis','isi'})
subplot(2,2,4)
boxplot([data(:,4)-data(:,7) data(:,4)-data(:,9) data(:,4)-data(:,10)]), hold on
plot(ones(size(comSim_mat,1),1) + 0.1*randn(size(comSim_mat,1),1), data(:,4)-data(:,7),'o',"MarkerFaceColor",[0.7,0.7,0.7])
plot(2*ones(size(comSim_mat,1),1) + 0.1*randn(size(comSim_mat,1),1), data(:,4)-data(:,9),'o',"MarkerFaceColor",[0.7,0.7,0.7])
plot(3*ones(size(comSim_mat,1),1) + 0.1*randn(size(comSim_mat,1),1), data(:,4)-data(:,10),'o',"MarkerFaceColor",[0.7,0.7,0.7])
[~,p1] = ttest(data(:,4)-data(:,7),0,"Tail","right");
[~,p2] = ttest(data(:,4)-data(:,9),0,"Tail","right");
[~,p3] = ttest(data(:,4)-data(:,10),0,"Tail","right");
text(0.6,-0.3,num2str(12*p1))
text(1.6,-0.3,num2str(12*p2))
text(2.6,-0.3,num2str(12*p3))
title('ISI')
yline(0)
xticklabels({'aud','vis','ctl'})
end

%__________________________________________________________________________
function comSim_dis = comSimDis(ComSim)
for i = 1:size(ComSim,2)
    sim_mat = ComSim{1,i}(2:end,:);
    % win_lables = ComSim{1,i}(1,:);
    cnt_dis = zeros(1, size(sim_mat,1));
    sim_dis = zeros(1, size(sim_mat,1));
    for m = 1:size(sim_mat,1)
        for n = 1:size(sim_mat,1)
            dis = abs(m-n);
            cnt_dis(dis+1) = cnt_dis(dis+1)+1;
            sim_dis(dis+1) = sim_dis(dis+1) + sim_mat(m,n);
        end
    end
    sim_dis = sim_dis./cnt_dis;
    comSim_dis{i} = sim_dis;
end
end
%__________________________________________________________________________
function m_x = plot_simDD(comSim_dis)
step = 5;
for i = 0:600/step-1
    sim_samples = [];
    for j = 1:size(comSim_dis,2)
        sim_samples = cat(2, sim_samples, comSim_dis{j}(i*step+1:(i+1)*step) );
    end
    m_x(i+1) = mean(sim_samples);
    sim_x(:,i+1) = sim_samples;
end
boxplot(sim_x)
end
%__________________________________________________________________________
function simInTime = sim_in_time(ComSim, win_sim)
for m = win_sim
    for n = 1:size(ComSim,2)
        TCSMat = ComSim{n}(2:end,:);
        TCSMat(eye(size(TCSMat))==1) = nan;
        cnt = 1;
        start_point = 1;
        end_point = 600;
        while start_point <= end_point-m
            simInTime{m}(n,cnt) = mean( TCSMat( start_point: start_point+m , start_point : start_point+m ), "all", "omitmissing" );
            start_point = start_point+m;
            cnt = cnt+1;
        end
    end
end
end
%__________________________________________________________________________
function aveTrialFlex = ave_trial_flex(labels, win_flex)
    % This function calculates the mean flexibility for consecutive identical labels
    % Inputs:
    %   labels - Array of labels for each window
    %   win_flex - Array of flexibility values corresponding to each window
    % Outputs:
    %   outLabels - Unique consecutive labels
    %   outMeanFlex - Mean flexibility values for each group of consecutive labels

    % Initialize output arrays
    outLabels = [];
    outMeanFlex = [];

    % Start processing
    startIdx = 1;
    while startIdx <= length(labels)
        % Get the current label
        currentLabel = labels(startIdx);

        % Find the end of the current consecutive segment
        endIdx = startIdx;
        while endIdx <= length(labels) && labels(endIdx) == currentLabel
            endIdx = endIdx + 1;
        end

        % Calculate mean flexibility for this segment
        meanFlex = mean(win_flex(startIdx:endIdx-1));

        % Store the results
        outLabels = [outLabels, currentLabel];
        outMeanFlex = [outMeanFlex, meanFlex];

        % Move to the next segment
        startIdx = endIdx;
    end

    aveTrialFlex = [outLabels', outMeanFlex'];
end
%__________________________________________________________________________
function task_tran_flex_pairs = flex_pairs(aveTrialFlex)
% task_tran_flex_pairs = [aveTrialFlex(find(aveTrialFlex(:,1)==1 | aveTrialFlex(:,1)==2),2), aveTrialFlex(find(aveTrialFlex(:,1)==1 | aveTrialFlex(:,1)==2)+1,2)];
% task_tran_flex_pairs = [aveTrialFlex(find(aveTrialFlex(:,1)==1) ,2), aveTrialFlex(find(aveTrialFlex(:,1)==2)+1,2)];
inx1 = find(aveTrialFlex(:,1)==1 | aveTrialFlex(:,1)==2);
if max(inx1)+1 <= length(aveTrialFlex)
    inx2 = inx1 + 1;
else
    inx1(end) = [];
    inx2 = inx1 + 1;
end
task_tran_flex_pairs = [aveTrialFlex(inx1,2), aveTrialFlex(inx2,2)];

end
%__________________________________________________________________________
function lme = mixed_effects_modeling(GroupNH, GroupCI1, GroupCI2)
    % Prepare the data for mixed-effects modeling
    numTimePoints = size(GroupNH, 2);
    GroupNH_flat = GroupNH(:); % Flatten NH group data
    GroupCI1_flat = GroupCI1(:); % Flatten CI1 group data
    GroupCI2_flat = GroupCI2(:); % Flatten CI2 group data

    % Combine all data
    allData = [GroupNH_flat; GroupCI1_flat; GroupCI2_flat];

    % Create a grouping variable for groups
    groupNH = repelem({'NH'}, numel(GroupNH_flat), 1);
    groupCI1 = repelem({'CI1'}, numel(GroupCI1_flat), 1);
    groupCI2 = repelem({'CI2'}, numel(GroupCI2_flat), 1);
    group = [groupNH; groupCI1; groupCI2];

    % Create a grouping variable for subjects
    subjectNH = repelem((1:size(GroupNH,1))', numTimePoints);
    subjectCI1 = repelem((1:size(GroupCI1,1))', numTimePoints);
    subjectCI2 = repelem((1:size(GroupCI2,1))', numTimePoints);
    subject = [subjectNH; subjectCI1; subjectCI2];

    % Create a time variable
    time = repmat((1:numTimePoints)', numel(group) / numTimePoints, 1);

    % Combine into a table
    tbl = table(allData, group, subject, time, ...
        'VariableNames', {'Values', 'Group', 'Subject', 'Time'});

    % Convert categorical variables and set GroupNH as the reference
    tbl.Group = categorical(tbl.Group, {'NH', 'CI1', 'CI2'}); % NH as reference
    tbl.Subject = categorical(tbl.Subject);

    % Fit a linear mixed-effects model
    lme = fitlme(tbl, 'Values ~ Group + Time + (1|Subject) + (1|Time)');
    
end

%__________________________________________________________________________
function [lme_NH lme_CI1 lme_CI2] = mixed_effects_modelingTime(GroupNH, GroupCI1, GroupCI2)
% Step 1: Compute average values for each group
ave_NH  = mean(GroupNH, 1)';  % Average curve for GroupNH
ave_CI1 = mean(GroupCI1, 1)'; % Average curve for GroupCI1
ave_CI2 = mean(GroupCI2, 1)'; % Average curve for GroupCI2

% Step 2: Create tables for each group
time = (1:length(ave_NH))'; % Time variable
tbl_NH  = table(ave_NH, time, 'VariableNames', {'Values', 'Time'});
tbl_CI1 = table(ave_CI1, time, 'VariableNames', {'Values', 'Time'});
tbl_CI2 = table(ave_CI2, time, 'VariableNames', {'Values', 'Time'});

% Step 3: Fit separate linear models for each group
lme_NH  = fitlme(tbl_NH,  'Values ~ Time');
lme_CI1 = fitlme(tbl_CI1, 'Values ~ Time');
lme_CI2 = fitlme(tbl_CI2, 'Values ~ Time');

% Step 4: Display only the main results (fixed effects table)
disp('NH Group:');
disp('Fixed effects coefficients (95% CIs):');
disp(lme_NH.Coefficients(:, {'Name', 'Estimate', 'SE', 'tStat', 'DF', 'pValue'}));

disp('CI1 Group:');
disp('Fixed effects coefficients (95% CIs):');
disp(lme_CI1.Coefficients(:, {'Name', 'Estimate', 'SE', 'tStat', 'DF', 'pValue'}));

disp('CI2 Group:');
disp('Fixed effects coefficients (95% CIs):');
disp(lme_CI2.Coefficients(:, {'Name', 'Estimate', 'SE', 'tStat', 'DF', 'pValue'}));

end

%__________________________________________________________________________
function [alleg_array, dis_array] = alleg_dis(alleg, dis_mat)
lvc_chs = [38, 41, 42];
rvc_chs = [40, 43, 44];
lac_chs = 24:37;
rac_chs = 10:23;
alleg_array = [];
dis_array = [];
for i = 1:size(alleg,4)
    alleg_mat = alleg(:,:,1,i);
    for j = 1:44
        for k = 1:44
            if j<k && ~isnan(alleg_mat(j,k))
               alleg_array = [alleg_array; alleg_mat(j,k)];
               dis_array = [dis_array; dis_mat(j,k)];
            end
        end
    end
end
end
%__________________________________________________________________________
function ciu44 = consensus_com(allegiance)
bad_chs = find(isnan(allegiance(:,1))); 
good_chs = 1:44;
good_chs(bad_chs) = [];
allegiance(bad_chs, :) = [];
allegiance(:, bad_chs) = [];
ciu = consensus_und(allegiance,0,50);
ciu44 = NaN(44,1);
ciu44(good_chs) = ciu;
end
%__________________________________________________________________________
function par_dis = par_disSub(consensusCom)
for m = 1:6
    for n = 1:6
        Cx = consensusCom(1,m,:);
        Cy = consensusCom(1,n,:);
        bad_chs = find(isnan(Cx));
        Cx(bad_chs) = []; Cy(bad_chs) = [];
        par_dis(m,n) = partition_distance(Cx,Cy);
    end
end
end


