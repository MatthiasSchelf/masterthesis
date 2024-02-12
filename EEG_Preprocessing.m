function EEG_Preprocessing = EEG_Preprocessing(EEG)

       
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

        % Update the new set
        EEG = pop_subcomp(EEG, EEG.reject.gcompreject);

        % save again 

        EEG = pop_saveset( EEG, 'savemode','resave');

        % Step 7: Interpolate bad channels
        
        EEG = pop_interp(EEG, badChannelIndices, 'spherical'); % You can use other interpolation methods as well

        %save again 
        EEG = pop_saveset( EEG, 'savemode','resave');
        EEG_Preprocessing = EEG; 
end
