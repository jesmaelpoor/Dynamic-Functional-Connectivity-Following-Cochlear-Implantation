% this script includes preprocessing steps

pipeline = nirs.modules.OpticalDensity(); % conversion to optical density
pipeline = nirs.modules.TDDR(pipeline); % TDDR to remove motion artifacts
pipeline = nirs.modules.BeerLambertLaw(pipeline); % conversion to hbo and hbr
preProcessed_data = pipeline.run(raw);
% short-channel labeling
job = nirs.modules.LabelShortSeperation();
preProcessed_data = job.run(preProcessed_data);
% short-channels correction and band-pass filtering
pipeline = advanced.nirs.modules.ShortDistanceFilter();

preProcessed_data = pipeline.run(preProcessed_data);

% bandpass filtering
pipeline = nirsGUI.modules.LowPassFilter(); 
pipeline.FilterOrder = 8;
pipeline.CutOffFrequency = lpcf;
pipeline = nirsGUI.modules.HighPassFilter(pipeline);
pipeline.FilterOrder = 8;
pipeline.CutOffFrequency = hpcf;
preProcessed_data = pipeline.run( preProcessed_data );

saveExperiment( Filtered_ex,sci_ex, subs )
saveResting( Filtered_rs,sci_rs, subs )

%% functions
% _______save recording functions_______
% this function saves processed time series and trial onsets
function saveExperiment(exData,sci, subs)
dir = 'C:\Melbourne University\Study Four-DFC\matlab code\data\flexibility_in_time\filteredData33\';
schs = [2,4,14,23,30,36,49,52];
chs = (1:52);
chs = chs(~ismember(chs,schs));

for i = 1:length(subs)
    hbo = exData(i,1).data(:,1:2:end);
    hbo = hbo(:,chs);
    hbr = exData(i,1).data(:,2:2:end);
    hbr = hbr(:,chs);
    recordingsci = sci(chs,i);   
    data = [[recordingsci';hbo],[recordingsci';hbr]];
    tevents = zeros(3,18);
    tevents(1,:) = exData(i,1).stimulus.values{1,1}.onset;
    tevents(2,:) = exData(i,1).stimulus.values{1,2}.onset;
    tevents(3,1:10) = exData(i,1).stimulus.values{1,3}.onset;    
    if ~ismember(str2double(subs(i,:)),[1,2,3,4,7]) %bad subs
        writematrix(data,[dir 'NHtask_'  subs(i,:) '.txt'],'Delimiter','tab')  
        writematrix(tevents,[dir 'NHonsets_'  subs(i,:) '.txt'],'Delimiter','tab')
    end
    clear data
end
end

