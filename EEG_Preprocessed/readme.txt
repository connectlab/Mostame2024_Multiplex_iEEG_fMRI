You can download the data from here:


And here are some points that you should know about the data:
0- The data is in Fieldtrip format
1- There is no temporal cleaning, meaning that I did not remove any time points (not to mess up with the dynamics)
2- However, I have the vector that marks bad time points throughout the whole timecourse, just to make sure findings are not based on artifactual time points.
3- I have done filtering on the data. (It should be 0.5 to 90 Hz, with a notch at 50 Hz)
4- Jump artifacts (< 100ms) has been detected and interpolated appropriately.
5- I only cleaned P07 using ICA since it had significant gradient artifacts. But not for other subjects.
6- These electrodes have been excluded from the data:
  * Two versions of EEG for each subject are shared. The file name clarifies that either (NIZ-IZ2) or (NIZ-only) has been included in the data. I use (NIZ-only) meaning that I * exclude both IZ1 and IZ2 regions from my analyses.
  * Electrodes in white matter have been excluded, using the SPM gray matter mask.
  * Additional electrodes with excessive ictal activity or artifacts have been removed based on my visual inspections.
