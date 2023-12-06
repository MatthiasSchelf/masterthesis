The steps indicate Makotos Preprocessing pipeline

% This tracks which version of EEGLAB is being used, you may ignore it
EEG.etc.eeglabvers = '2023.1'; 

% This command resamples the data. pop_resample(your EEG data structure, the sampling frequency you want) = step 1 downsampling
EEG = pop_resample( EEG, 250);

% Here we decide the name of the dataset, in this case sub 33
EEG.setname='sub-33-task-memory-pre-pro';

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

% Manual rejection of an artefact of a dataset = step 4 remove bad channels 
pop_eegplot( EEG, 1, 1, 1);

% Indicate which channels is a refernce channel
EEG=pop_chanedit(EEG, 'append',66,'changefield',{67,'labels','FCz'},'lookup','\\\\client\\c$\\eeglab_current\\eeglab2023.1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc','setref',{'1:67','FCz'});

% Rereference the data, in this case average reference = Step 5 Average reference 
EEG = pop_reref( EEG, [],'refloc',struct('labels',{'FCz'},'type',{''},'theta',{0.7867},'radius',{0.095376},'X',{27.39},'Y',{-0.3761},'Z',{88.668},'sph_theta',{-0.7867},'sph_phi',{72.8323},'sph_radius',{92.8028},'urchan',{67},'ref',{'FCz'},'datachan',{0}));

% This code was presented in the Makotos preprocessing pipeline to do average referencing. 
EEG.data = bsxfun( @minus, EEG.data, sum( EEG.data, 1 ) / ( EEG.nbchan + 1 ) );


% Step 6: ICA -> use the automated way 

% Step 7: Interpolate bad channels
