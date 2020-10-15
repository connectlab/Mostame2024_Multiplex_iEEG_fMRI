% The function changes fMRI data to match to EEG electrodes, electrode order, and time interval.
% --------- subject: name of the subject as a string
% --------- fMRI: fMRI 2D matrix (ROI, time)
% --------- fMRI_time: time vector in seconds
function [fMRI, fMRI_time] = rematch_electrode_orders(subject, fMRI, fMRI_time)
if strcmp( subject, 'P05' )
    % rearrange electrodes
    temp = fMRI( 81:86 , :);
    fMRI( 81:86 , :) = [];
    fMRI = [fMRI; temp];
    % remove EEG-missing electrodes
    fMRI( [48 56 59:64 65 73 ], :) = [];
    % remove EEG-missing times
    T_End = 280; % doesn't have EEG after 280s
    fMRI_time( fMRI_time > T_End ) = [];
    fMRI = fMRI(:, 1 : fMRI_time );
    
elseif strcmp( subject, 'P07' )
    % remove EEG-missing electrodes
    fMRI( 63 : 64 , :) = [];
    
elseif strcmp( subject, 'P03' )
    % remove EEG-missing electrodes
    A = [ 25:88 17:24 1:16 ];
    A( ismember( A, [ 5 8 13:16 84 86 ] ) ) = [];
    fMRI  = fMRI(A, :);
    
elseif strcmp( subject, 'P02' )
    % remove EEG-missing electrodes
    A = [ 20 44 45:60 ];
    fMRI(A, :) = [];
    
elseif strcmp( subject, 'P04' )
    % remove EEG-missing electrodes
    A = [ 15:21 28 35:36 43 ];
    fMRI(A, :) = [];
    % rearrange electrodes
    temp = fMRI(9:14, :);
    fMRI(9:14, :) = [];
    fMRI = [ fMRI; temp ];
    
elseif strcmp( subject, 'P09' )
    % remove EEG-missing electrodes
    fMRI( [ 18 40 45 49 67 ] , :) = [];
    % fMRI-missing ROIs are filled with constant 1 timecourses here, but are
    % later excluded from both fMRI and EEG because of bad EEG.
    L = size(fMRI, 2);
    fMRI = [fMRI(1:16, :); ones(9, L); fMRI(17:40, :); ones(2, L); fMRI(41:45, :); ones(1, L); fMRI(46:56, :); ones(2, L); fMRI(57:end, :) ];
end
end

