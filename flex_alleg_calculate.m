% this code is to do:
% 1. dynamic modude detection (50 iteration) for each task/rest recording,
% 2. create 50 temporal & nodal null models for each task/rest recording and do
% dynamic module detection for each of them (30 iteration).
% 3. measure flexibility (for each node) in the networks.
% save variables for each recording (dynamic modules for the primary nets, node flexibilities and allegience matrices for primary & null nets)

clc, close all, clear all,

record_names = dir('C:\Melbourne University\Study Four-DFC\matlab code\data\flexibility_in_time\filteredData30'); %recording names

HB = 'hbo'; %choose hbo or hbr
%parameters
global win_len win_step gamma omega
sci_thre = 0.5;         %sci_thershold
gamma = 1;       %intralayer resolution parameter  
omega = 0.2;        %interlayer coupling strength
Fs = 3.9125;        %sampling frequencies 
win_len = round(Fs*18);       %window length
win_step = round(Fs*3);      %window step
cnt_records = 0;
for i = 238:numel(record_names)         %do the analysis for each recording
    if ~(strcmp(record_names(i,1).name(1:end-7),"NHonsets") || strcmp(record_names(i,1).name(1:end-7),"onsets_")) 
        cnt_records = cnt_records + 1;
        disp([num2str(cnt_records) '/119. ' record_names(i,1).name(1:end)])
        disp(datetime)
        recording = readmatrix([record_names(i,1).folder, '\', record_names(i,1).name(1:end-4) ]);         %read recording
        %remove bad channels
        sci_recording = recording(1,:);
        bad_chs = find(sci_recording<sci_thre);
        good_chs = find(sci_recording>=sci_thre);
        recording(:,bad_chs) = [];
        %daynamic module detection
        Sopt = dynamicModule( recording, HB );       %dynamic module detection of the recording: 30 trials
        StempNull = nan; %%%%%dynamicModule_tempNull( recording, HB );       %dynamic module detection for the temporal num models: 50 temporal null networks & 30 trials for each-->1500
        SnodNull = nan; %%%%%%%dynamicModule_nodalNull( recording, HB );         %dynamic module detection for the nodal num models: 50 nodal null networks & 30 trials for each-->1500
        %flexibility & coherence for nets    
        [F_opt, C_opt] = flex_coh(Sopt);
        F_opt = cat(2, F_opt, good_chs(1:length(good_chs)/2)'); %include good chs with F_opt variable
        F_tempNull = nan;
        C_tempNull = nan;
        F_nodNull = nan; 
        C_nodNull = nan;
        % [F_tempNull, C_tempNull] = flex_coh(StempNull);
        % [F_nodNull, C_nodNull] = flex_coh(SnodNull); C_nodNull] = flex_coh(SnodNull);
        %save results
        save(['C:\Melbourne University\Study Four-DFC\matlab code\data\parameter selection\dynFC_8\' record_names(i,1).name(1:end-4) '.mat'], "Sopt", "F_opt", "C_opt", "F_tempNull", "F_nodNull", "C_tempNull","C_nodNull")  
    end
end

%% functions
%%%%%%%%%%%%%%%
% dynamic module detection: main nets 
function Sopt = dynamicModule( recording, HB )
global win_len win_step gamma omega
switch HB %choose hbo or hbr
    case 'hbo'
        recording(:,size(recording,2)/2+1:end) = [];
    case 'hbr'
        recording(:,1:size(recording,2)/2) = [];
end
%windowing
start_n = 1;
end_n = win_len;
cnt_w = 1;
while end_n < size(recording,1)        
    A{cnt_w} = corr(recording(start_n:end_n,:));
    start_n = start_n + win_step;
    end_n = end_n + win_step;
    cnt_w = cnt_w + 1;
end
cnt_w = cnt_w - 1; %number of layers    
N=length(A{1});
T=length(A);    
k = 1;
while k <= 50   %dynamic module detection
      try
          [B,twom] = multiord(A,gamma,omega);
          PP = @(S) postprocess_ordinal_multilayer(S,T);
          [S,Q,n_it] = iterated_genlouvain(B,10000,0,1,'moverand',[], PP);         
          S = reshape(S, N, T);
          S = postprocess_ordinal_multilayer(S);
          Sopt(:,:,k) = multislice_pair_labeling(S);
          % num_community = length(unique(SOpt(:,:,k)))
          %imagesc(SOpt(:, :, k)), pause  
          k = k+1;
      catch
      end
end 
end

%%%%%%%%%%%%%%%%
% dynamic module detection: temporal null nets 
function StempNull = dynamicModule_tempNull( recording, HB )
global win_len win_step gamma omega
switch HB %choose hbo or hbr
    case 'hbo'
        recording(:,size(recording,2)/2+1:end) = [];
    case 'hbr'
        recording(:,1:size(recording,2)/2) = [];
end
%windowing
start_n = 1;
end_n = win_len;
cnt_w = 1;
while end_n < size(recording,1)        
    Amain{cnt_w} = corr(recording(start_n:end_n,:));
    start_n = start_n + win_step;
    end_n = end_n + win_step;
    cnt_w = cnt_w + 1;
end
cnt_w = cnt_w - 1;      %number of layers    
N=length(Amain{1});
T=length(Amain);    
cntNets = 0;
for i = 1:50        %50 null models
    A = permuteLayers(Amain);
    k = 1;
    while k <= 50   %dynamic module detection
          try
              [B,twom] = multiord(A,gamma,omega);
              PP = @(S) postprocess_ordinal_multilayer(S,T);
              [S,Q,n_it] = iterated_genlouvain(B,10000,0,1,'moverand',[], PP);               
              S = reshape(S, N, T);
              S = postprocess_ordinal_multilayer(S);
              cntNets = cntNets + 1;
              StempNull(:,:,cntNets) = multislice_pair_labeling(S);              
              k = k+1;
        catch
        end
    end
end
end

%%%%%%%%%%%%%%%%
% permute layers
function A = permuteLayers(Amain)
layer_order = randperm(numel(Amain));
for i = 1:numel(Amain)
    A{i} = Amain{layer_order(i)};
end
end

%%%%%%%%%%%%%%%
% dynamic module detection: nodal null nets 
function SnodNull = dynamicModule_nodalNull( recording, HB )
global win_len win_step gamma omega
switch HB %choose hbo or hbr
    case 'hbo'
        recording(:,size(recording,2)/2+1:end) = [];
    case 'hbr'
        recording(:,1:size(recording,2)/2) = [];
end
cntNets = 0;
for i = 1:50        %50 nodal null models
    %recordingNull = permuteNodes(recording);
    start_n = 1;
    end_n = win_len;
    cnt_w = 1;
    while end_n < size(recording,1)        
        A{cnt_w} = corr(permuteNodes(recording(start_n:end_n,:)));
        start_n = start_n + win_step;
        end_n = end_n + win_step;
        cnt_w = cnt_w + 1;
    end
    cnt_w = cnt_w - 1;      %number of layers    
    N=length(A{1});
    T=length(A);    
    k = 1;
    while k <= 50   %dynamic module detection
          try
              [B,twom] = multiord( A, gamma, omega );
              PP = @(S) postprocess_ordinal_multilayer(S,T);
              [ S,Q,n_it ] = iterated_genlouvain(B,10000,0,1,'moverand',[], PP);              
              S = reshape(S, N, T);
              S = postprocess_ordinal_multilayer(S);
              cntNets = cntNets + 1;
              SnodNull(:,:,cntNets) = multislice_pair_labeling(S);             
              k = k+1;
          catch
          end
    end
end
end

%%%%%%%%%%%%%%%%
% permute nodes
function recordingNull = permuteNodes(recording)
node_order = randperm(size(recording,2));
for i = 1:size(recording,2)
    recordingNull(:,i) = recording(:,node_order(i));
end
end

%%%%%%%%%%%%%%
function [F, C] = flex_coh(S)
for i = 1:size(S,3)
    Sw = S(:,:,i);
    node_flexibility(:,i) = flexibility(Sw');
    mat_allegiance(:,:,i) = allegiance(Sw);   
end
n = 0: 50: size(S,3);
for i = 1:length(n)-1
    C(:,:,i) = mean(mat_allegiance(:,:,n(i)+1:n(i+1)),3);
    F(:,i) = mean(node_flexibility(:,n(i)+1:n(i+1)),2);
end
end
%%%%%%%%%%%%
function mat_allegiance = allegiance(S)
for m = 1:size(S,1)
    for n = 1:size(S,1)
        mat_allegiance(m,n) = length(find(S(n,:)==S(m,:)))/size(S,2);
    end
end
end


