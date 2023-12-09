% Here we make a function with all the different steps that need to be done in order to preproces the data. 

function processed_data = process_eeglab_data(dataset_params)

        % This tracks which version of EEGLAB is being used, you may ignore it
        EEG.etc.eeglabvers = '2023.1'; 
       
        % This command resamples the data. pop_resample(your EEG data structure, the sampling frequency you want) = step 1 downsampling
        EEG = pop_resample( EEG, 250);

        % Here we decide to save the dataset.
        EEG = pop_saveset( EEG, 'savemode','resave');

        
        
        
        % Here we filter the data. In this case a high-pass 1 Hz filter = step 2 high pass filter 1 Hz 
        EEG = pop_eegfiltnew(EEG, 'locutoff',1,'plotfreqz',1);

        % Save again 
        EEG = pop_saveset( EEG, 'savemode','resave');

       
       

       
        % Here we look up the names of the channels and lock them to the right place/sensors = step 3 import channel names
        EEG=pop_chanedit(EEG, 'lookup','\\\\client\\c$\\eeglab_current\\eeglab2023.1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc','lookup','\\\\client\\c$\\eeglab_current\\eeglab2023.1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc','eval','chans = pop_chancenter( chans, [],[]);','eval','chans = pop_chancenter( chans, [],[]);');

        % Save again 
        EEG = pop_saveset( EEG, 'savemode','resave');

        
        
        
        % Manual rejection of an artefact of a dataset = step 4 remove bad channels and indicate the reference channel
        
        pop_eegplot( EEG, 1, 1, 1);

        % Indicate which channels is a reference channel, normally this stays the same over all the datasets. 
        EEG=pop_chanedit(EEG, 'append',66,'changefield',{67,'labels','FCz'},'lookup','\\\\client\\c$\\eeglab_current\\eeglab2023.1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc','setref',{'1:67','FCz'});

        % Now we can remove the bad channels manually every time but this will take a lot of time. Here it is automated. 

        % Automatic bad channel rejection using clean_rawdata() function
        % FlatlineCriterion is used to see if the channel is a flat line or not, in this case it will be deleted if it has a sd of lower than 5
        % ChannelCriterion removes a channel if it has a high correlation with another channel (here a correlation of 0.8 or higher)
        % LineNoiseCriterion rejects channels that have to much line noise (in this case line noise power 4 times that of the normal power)
        % Highpass indicates if their is a high pass filter or not, in this case it is off
        % Burst cirterion indicate if bursts need to be removed or not. In this case it is off
        % Windowcriterion indicates data segments with extreme variance. Here it is turned off. 

        EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion', 5, 'ChannelCriterion', 0.8, 'LineNoiseCriterion', 4, 'Highpass', 'off', 'BurstCriterion', 'off', 'WindowCriterion', 'off');





        % Rereference the data, in this case average reference = Step 5 Average reference 
        EEG = pop_reref( EEG, [],'refloc',struct('labels',{'FCz'},'type',{''},'theta',{0.7867},'radius',{0.095376},'X',{27.39},'Y',{-0.3761},'Z',{88.668},'sph_theta',{-0.7867},'sph_phi',{72.8323},'sph_radius',{92.8028},'urchan',{67},'ref',{'FCz'},'datachan',{0}));

        % Now the previous line was only for one dataset, this is to automate it:
        % Apply average reference
        EEG = pop_reref(EEG, []);



        % Step 6: ICA -> use the automated way = IClabel way 

        % Run ICA through the runica funtion
        EEG = pop_runica(EEG, 'extended', 1, 'interupt', 'on');

        % label the components according to whether or not they are considered an artefact or noise
        EEG = iclabel(EEG);

        % deciding an (arbitrary) treshold to reject the ICA that are considered artefact or noise (here probablity >0.8 will be rejected)
        EEG = pop_icflag(EEG, [NaN NaN; 0.8 1; 0.8 1; 0 0; 0 0; 0 0], 'fields', 1);

        % Updae the new set
        EEG = pop_subcomp(EEG, find(EEG.reject.gcompreject));

        % save again 

        EEG = pop_saveset( EEG, 'savemode','resave');




        % Step 7: Interpolate bad channels
        
        EEG = pop_interp(EEG, badChannelIndices, 'spherical'); % You can use other interpolation methods as well

        %save again 
        EEG = pop_saveset( EEG, 'savemode','resave');



        processed_data = EEG; 
end


% Loop over all the datasets. 


% Common part of the dataset path
common_path = 'D:/UGent_gerelateerd/Masterproef/Data/RawData/';

% Variable parts for the last two blocks of the path
participants = {'sub-032/eeg/sub-032_task-memory_eeg', 'sub-033/eeg/sub-033_task-memory_eeg', 'sub-034/eeg/sub-034_task-memory_eeg','sub-035/eeg/sub-035_task-memory_eeg',
'sub-036/eeg/sub-036_task-memory_eeg','sub-038/eeg/sub-038_task-memory_eeg','sub-039/eeg/sub-039_task-memory_eeg','sub-040/eeg/sub-040_task-memory_eeg','sub-041/eeg/sub-041_task-memory_eeg',
'sub-042/eeg/sub-042_task-memory_eeg','sub-043/eeg/sub-043_task-memory_eeg','sub-044/eeg/sub-044_task-memory_eeg','sub-045/eeg/sub-045_task-memory_eeg','sub-046/eeg/sub-046_task-memory_eeg',
'sub-047/eeg/sub-047_task-memory_eeg','sub-048/eeg/sub-048_task-memory_eeg','sub-049/eeg/sub-049_task-memory_eeg','sub-050/eeg/sub-050_task-memory_eeg','sub-051/eeg/sub-051_task-memory_eeg',
'sub-052/eeg/sub-052_task-memory_eeg','sub-053/eeg/sub-053_task-memory_eeg','sub-054/eeg/sub-054_task-memory_eeg','sub-055/eeg/sub-055_task-memory_eeg','sub-056/eeg/sub-056_task-memory_eeg',
'sub-057/eeg/sub-057_task-memory_eeg','sub-058/eeg/sub-058_task-memory_eeg','sub-059/eeg/sub-059_task-memory_eeg','sub-060/eeg/sub-060_task-memory_eeg','sub-061/eeg/sub-061_task-memory_eeg',
'sub-062/eeg/sub-062_task-memory_eeg','sub-063/eeg/sub-063_task-memory_eeg','sub-064/eeg/sub-064_task-memory_eeg','sub-065/eeg/sub-065_task-memory_eeg','sub-067/eeg/sub-067_task-memory_eeg',
'sub-068/eeg/sub-068_task-memory_eeg','sub-069/eeg/sub-069_task-memory_eeg','sub-070/eeg/sub-070_task-memory_eeg','sub-071/eeg/sub-071_task-memory_eeg','sub-072/eeg/sub-072_task-memory_eeg',
'sub-073/eeg/sub-073_task-memory_eeg','sub-074/eeg/sub-074_task-memory_eeg','sub-075/eeg/sub-075_task-memory_eeg','sub-076/eeg/sub-076_task-memory_eeg','sub-077/eeg/sub-077_task-memory_eeg',
'sub-078/eeg/sub-078_task-memory_eeg','sub-079/eeg/sub-079_task-memory_eeg','sub-080/eeg/sub-080_task-memory_eeg','sub-081/eeg/sub-081_task-memory_eeg','sub-082/eeg/sub-082_task-memory_eeg',
'sub-083/eeg/sub-083_task-memory_eeg','sub-084/eeg/sub-084_task-memory_eeg','sub-085/eeg/sub-085_task-memory_eeg','sub-086/eeg/sub-086_task-memory_eeg','sub-087/eeg/sub-087_task-memory_eeg',
'sub-088/eeg/sub-088_task-memory_eeg','sub-089/eeg/sub-089_task-memory_eeg','sub-090/eeg/sub-090_task-memory_eeg','sub-091/eeg/sub-091_task-memory_eeg','sub-092/eeg/sub-092_task-memory_eeg',
'sub-093/eeg/sub-093_task-memory_eeg','sub-095/eeg/sub-095_task-memory_eeg','sub-096/eeg/sub-096_task-memory_eeg','sub-097/eeg/sub-097_task-memory_eeg','sub-098/eeg/sub-098_task-memory_eeg',};

% Loop over participants
for i = 1:length(participants)
   
    % Construct the full dataset path
    data_path = fullfile(common_path, participants{i});

    % Call the processing function
    processed_data = eeglab_processing_function(data_path);

    % Save the processed_data as needed
    save_results(processed_data, data_path);
end



