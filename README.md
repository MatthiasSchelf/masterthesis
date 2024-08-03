This code was used in the preprocessing and analysis of the Masterthesis titled "Evaluating different measures of working memory (overload) on the basis of Mutual Information (MI)".

This code consists of the preprocessing code and the analysis itself. 

The preprocessing file can be used as follows: 

  First, use the datareading R files to turn the files to .mat format. 
  
  Second, use HEPLAB as described in the paper for the ECG data.
  
  Third, use the controlcode to loop over all the data. 
  
  Fourth, check if all files have a sampling frequency of 120 Hz by using the resample files. 
  
  Fifth, the preprocessing files are referenced in the controlcode files and will be used automtically. 

The analysis file can be used as follows:

  First, the new_pupil_eeg_project can be used to get the Spearman correlation and MI for every participant separately between the EEG and pupillometry.
  
  Second, the new_pupil_eeg_ecg_project can be used to get the CMI for every participant separately between the EEG and pupillometry while contolling for the ECG.
  
  Third, the final_results_eeg_pupil can be used to combine all spearman correlations and MI together to get the final result, which is also reported in the paper       itself.
  
  Fourth, the final_results_eeg_pupil_ecg can be used to combine all CMI together to get the final result, which is also reported in the paper itself. 
  Fifth, all other codes are necessary since they support the main four codes. 

The code for the anlysis is based on the code of Ince (https://github.com/robince/gcmi), in which the Information Theoretic Framework and its functionalities are put into easy to use code. 
