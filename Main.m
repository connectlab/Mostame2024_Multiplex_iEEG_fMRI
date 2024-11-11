%% ____________________________________________ initialize
clear all; close all; clc;
% Fieldtrip directory
addpath( '..\fieldtrip-20160912' );
ft_defaults;
% Data directories
Path.path_main = '..';
Path.path_results = [Path.path_main '\Results_elife'];
Path.path_EEG = [Path.path_main '\EEG'];
Path.path_fMRI = [Path.path_main '\fMRI'];
Path.path_movement= [Path.path_main '\movement'];
Path.path_electrodes = [Path.path_main '\electrodes'];
% function directories
Path.path_func = '..\functions'; addpath(Path.path_func);





%% ____________________________________________ Load needed data of all modalities
mode = 'NIZ-only';
Reref = 'CmnAvg';
REMOVESPIKES = 0;

HRF = 0;

i_subs = [1:9];
freqs = 1:5;

for i_sub = i_subs
   Main_code( Path, i_sub, freqs, mode, Reref, HRF);
end





%% ____________________________________________ Main code
function [] = Main_code( Path, i_sub, freqs, mode, Reref, HRF)

cd( [ Path.path_main '\codes' ] ); load( 'Subjects.mat' );
clc; subject = Subjects{i_sub}


% --------- load iEEG data
cd( [Path.path_EEG '\' subject '\Preprocessed'] );
load( sprintf('%s_EEG_marked_%s_%s.mat', subject, mode, Reref) );
load( sprintf('%s_Bad_times_EEG.mat', subject) );


% --------- load fMRI data
TR = 3; 
cd( [Path.path_fMRI '\' subject] );
load(sprintf('%s_fMRI_regout_8mm_%s_forHRFanalysis.mat', subject, mode));
load( sprintf('%s_Bad_times_fMRI_8mm.mat', subject) );


% optionally truncate the first 6s of fMRI for hemodynamic delay
if HRF == 0
   fMRI_regout = fMRI_regout(:, floor(6/TR) + 1 : end);
end

% low-sample fMRI
fMRI_lowsample = fMRI_regout; clear fMRI_regout
fMRI_time_lowsample = [ 0 : (size(fMRI_lowsample, 2) - 1) ] * TR;

% high-sample fMRI
% fMRI_time = fMRI_time_lowsample(1) : STEP : ( fMRI_time_lowsample(end) ); 
% fMRI = interp1( fMRI_time_lowsample, fMRI_lowsample', fMRI_time, 'spline')';
% TR = STEP;

fMRI_time = fMRI_time_lowsample;
fMRI = fMRI_lowsample;




%% ____________________________________________ Extract good times for all modalities
EEG_time = EEG.time{1};
EEG_time = EEG_time(EEG_time <= fMRI_time(end) + TR);

Good_times_all_EEG = ones( size( EEG_time ) );
Good_times_all_fMRI = ones( size( fMRI_time ) );



%% ____________________________________________ clean both data
% clean EEG
EEG.trial{1} = EEG.trial{1}(:, Good_times_all_EEG == 1);
EEG.time{1} = EEG_time(1) : 1 / EEG.fsample : ( size( EEG.trial{1}, 2) - 1 )/ EEG.fsample;
EEG_time = EEG.time{1};
EEG.sampleinfo(2) = numel( EEG.time{1} );

% downsampled EEG time vector for eFC stimations
% EEG_time_ds = downsample(EEG_time, floor( STEP * EEG.fsample ) );
EEG_time_ds = fMRI_time;

% clean fMRI
fMRI = fMRI(:, Good_times_all_fMRI == 1);




%% ____________________________________________ Extract ROI information
% load electrode information
filename = 'electrodelabels.xlsx';
elec_info = load_elecinfo(filename, Path.path_electrodes, subject, mode);
locs = cell2mat( elec_info(:, 8:10) );

% extract type of electrode implantation
depth = cellfun(@(x) strcmp(x, 'depth'), elec_info(:, 5) );

% mark white matter electrodes to get removed
Gray = cell2mat( elec_info(:, 11) );

% EEG badelecs are already removed in the cleaning process. fMRI badelecs should be removed now
Badelecs = ~Gray | (cell2mat( elec_info(:, 7) ) == 1);
clear Gray
fMRI( Badelecs, : ) = [];
EEG = ft_selectdata(EEG, 'channel', elec_info(~Badelecs, 1));
elec_info( Badelecs , :) = [];

% extract number of electrodes
numelec = numel( EEG.label );


% ------------------------------ extract ROI distance (in cm) and implant type connection
dist_electrodes = nan(numelec, numelec);
ImplantType_Connections = dist_electrodes;
for i = 1 : size(dist_electrodes, 1)
   for j = 1 : size(dist_electrodes, 1)
      if i < j
         dist_electrodes(i, j) = 0.1 * sqrt( sum( abs( locs(i, :) - locs(j, :) ).^2 ) );
         dist_electrodes(j, i) = dist_electrodes(i, j);
         ImplantType_Connections(i, j) = ismember(i, depth) + ismember(j, depth) + 1;
         ImplantType_Connections(j, i) = ImplantType_Connections(i, j);
      end
   end
end





%% ____________________________________________ fMRI FC analyses (Betzel et al measure)
FC_fMRI_Package = [];

fMRI_zscore = zscore(fMRI')';

eFC = nan( size(fMRI_zscore, 1), size(fMRI_zscore, 1),  numel(EEG_time_ds) );
for i = 1 : size(fMRI_zscore, 1)
   for j = 1 : size(fMRI_zscore, 1)
      if i < j && dist_electrodes(i, j) >= 0.9
%          temp = interp1(fMRI_time, fMRI_zscore(i, :) .* fMRI_zscore(j, :), EEG_time_ds, 'spline');
         temp = fMRI_zscore(i, :) .* fMRI_zscore(j, :);
         eFC(i, j, :) = temp;
         eFC(j, i, :) = eFC(i, j, :);
      end
   end
end

FC_fMRI = eFC; clear eFC

% fMRI static FC
FC_fMRI_static = nanmean( FC_fMRI, 3 );
[FC_fMRI_static_regout,~,~,~] = Dist_Reg_Out(FC_fMRI_static, dist_electrodes);

% save in the fMRI Package
FC_fMRI_Package.FC_fMRI_MTD = FC_fMRI;
FC_fMRI_Package.FC_fMRI_static_MTD = FC_fMRI_static;
FC_fMRI_Package.FC_fMRI_static_regout_MTD = FC_fMRI_static_regout;

% Extract NAN FC values
NANS = isnan( FC_fMRI_Package.FC_fMRI_static_MTD );



%% ____________________________________________ iEEG FC analyses

if HRF == 0
   FC_EEG_AmpC = cell( 1, numel( freqs ) ); FC_EEG_AmpC_static = FC_EEG_AmpC; FC_EEG_AmpC_static_regout = FC_EEG_AmpC;
   FC_EEG_PhC = FC_EEG_AmpC; FC_EEG_PhC_static = FC_EEG_AmpC; FC_EEG_PhC_static_regout = FC_EEG_AmpC;
      
   for freq = freqs
      freq
      [ FC_EEG_AmpC{freq}, FC_EEG_AmpC_static{freq}, FC_EEG_AmpC_static_regout{freq} ] = ...
         AmpC_or_PhC_eFC('AmpC', [], EEG.trial{1}, EEG_time, EEG_time_ds, freq, TR,...
         dist_electrodes, []);
      
      [ FC_EEG_PhC{freq}, FC_EEG_PhC_static{freq}, FC_EEG_PhC_static_regout{freq} ] =  ...
         AmpC_or_PhC_eFC('PhC', [], EEG.trial{1}, EEG_time, EEG_time_ds, freq, TR,...
         dist_electrodes, []);

      [m, n] = find( ( NANS == 1 ) );
      while ~isempty( m )
         FC_EEG_AmpC{freq}( m(1), n(1), :) = nan(1, size(FC_EEG_AmpC{freq}, 3));
         FC_EEG_PhC{freq}( m(1), n(1), :) = nan(1, size(FC_EEG_PhC{freq}, 3));
         m(1) = []; n(1) = [];
      end
      
      FC_EEG_AmpC_static{freq}( NANS ) = nan;
      FC_EEG_AmpC_static_regout{freq}( NANS ) = nan;
      FC_EEG_PhC_static{freq}( NANS ) = nan;
      FC_EEG_PhC_static_regout{freq}( NANS ) = nan;
   end
   
   
   
elseif HRF == 1
   addpath('Y:\mostame2\spm12'); RT = diff(EEG_time); RT = RT(1); canonical_hrf = spm_hrf( RT ); clear RT
      
   FC_EEG_AmpC = cell( 1, numel( freqs ) ); FC_EEG_AmpC_static = FC_EEG_AmpC; FC_EEG_AmpC_static_regout = FC_EEG_AmpC;
      
   for freq = freqs
      freq
      [ FC_EEG_AmpC{freq}, FC_EEG_AmpC_static{freq}, FC_EEG_AmpC_static_regout{freq} ] = ...
         AmpC_or_PhC_eFC('AmpC', [], EEG.trial{1}, EEG_time, EEG_time_ds, freq, TR, ...
         dist_electrodes, canonical_hrf);
      
      [m, n] = find( ( NANS == 1 ) );
      while ~isempty( m )
         FC_EEG_AmpC{freq}( m(1), n(1), :) = nan(1, size(FC_EEG_AmpC{freq}, 3));
         m(1) = []; n(1) = [];
      end
      
      FC_EEG_AmpC_static{freq}( NANS ) = nan;
      FC_EEG_AmpC_static_regout{freq}( NANS ) = nan;
   end
end


%% ____________________________________________ Save required data/matrices
clear counter elec FC_fMRI FC_fMRI_zscored FC_fMRI_static FC_fMRI_static_regout freqs h i j k m n t
cd( Path.path_results );

if HRF == 1
   HRF_txt = '_HRF';
elseif HRF == 0
   HRF_txt = '';
end
txt = sprintf('%s_%s_eFC_TR%s.mat', subject, mode, HRF_txt);
save(txt, '-v7.3')


end
