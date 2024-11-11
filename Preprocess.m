%% initialize
% Fieldtrip directory
clear; clc
addpath( '..\fieldtrip-20160912' );
ft_defaults;
path_func = '..\functions'; addpath(path_func);
path_main = '..\fMRI-iEEG';
path_electrodes = [path_main '\electrodes'];
path_EEG = [path_main '\EEG'];
path_fMRI = [path_main '\fMRI'];
path_movement= [path_main '\movement'];
cd( [path_main '\codes'] )
load Subjects.mat
i_subs  = [2];
TR = 3;

%% clean the EEG data
mode = 'NIZ-only'
Reref = 'CmnAvg'

i_subs = 5;
for i_sub=i_subs
   clc
   subject = Subjects{i_sub}
   cd( [path_EEG '\' subject] );
   if ~exist('Preprocessed', 'dir')
      mkdir('Preprocessed');
   end
   path_data = [path_EEG '\' subject '\Preprocessed'];
   temp = dir();
   for k = 1:numel(temp)
      txt = temp(k).name;
      if endsWith(txt, '.set') || endsWith(txt, '.vhdr')
         break
      end
   end
   
   % load raw EEG (the data is already filtered)
   cfg = [];
   cfg.dataset = txt;
   cfg.channels = 'all';
   EEG_raw = ft_preprocessing(cfg);
   T_End = 600;
   if strcmp( subject, 'P05' )
      cfg1 = [];
      cfg1.latency = [0 70331/EEG_raw.fsample] ;
      cfg1.channel = ft_channelselection({'all', '-DA6'}, EEG_raw.label);
      EEG_raw = ft_selectdata( cfg1, EEG_raw);
      T_End = 280;
   end
   
   % identify bad electrodes from the excel file
   cd( path_electrodes )
   temp = EEG_raw.label;
   temp( find( cellfun( @(x) strcmp(x, 'ECG') || strcmp(x, 'ECG2'), temp ) ) ) = [];
   if strcmp( subject, 'P05' )
      temp = sprintf( 'G%d', numel(temp) + 2 );
   else
      temp = sprintf( 'G%d', numel(temp) + 1 );
   end
   [~, ~, elec_info] = xlsread('electrodelabels.xlsx', subject, [ 'A2:' temp ]);
   ictal = cell2mat( elec_info(:, 3:4) );
   switch mode
      case 'NIZ-only'
         ictal = or( ictal(:, 1), ictal(:, 2) );
      case 'NIZ-IZ2'
         ictal = ictal(:, 1);
   end
   Badelecs_EEG = find( cell2mat( elec_info(:, 6) ) == 1 );
   ictal(Badelecs_EEG) = 1;
   non_ictal_electrodes_EEG = elec_info(ictal == 0, 1); clear ictal
   
   % exclude ictal channels + ECG
   EEG_raw = ft_selectdata(EEG_raw, 'channel', non_ictal_electrodes_EEG);
   
   % identify extremely large spikes
   cd( [path_EEG '\' subject] );
   cfg= [];
   cfg.dataset = txt;
   cfg.trialdef.ntrials = 1;
   cfg.trialdef.triallength = inf;
   cfg = ft_definetrial(cfg);
   cfg.artfctdef.zvalue.channel = non_ictal_electrodes_EEG;
   cfg.artfctdef.zvalue.cutoff = 7;
   [~, artifact] = ft_artifact_zvalue(cfg);
   
   % interpret data of spike moments
   flag = zeros( size( artifact, 1 ), 1 );
   for k = size(artifact, 1): -1 : 1
      if diff( artifact(k, :) ) <= 100 && ( artifact(k, 2) < find( EEG_raw.time{1} > T_End, 1, 'first' ) )
         % To ensure the case of a equal to b is ok, we widen the window
         a = artifact(k, 1) - floor( 0.02 * EEG_raw.fsample ); a = max(a, 1);
         b = artifact(k, 2) + floor( 0.02 * EEG_raw.fsample );
         % get more data around the issue to interpolate
         L = 10 * EEG_raw.fsample - (b - a); L = floor( L / 2 );
         Start = a -  L;
         End = b + L;
         if Start < 1
            Start = 1;
         end
         if End > numel(EEG_raw.time{1})
            End = numel(EEG_raw.time{1});
         end
         % interpolate
         for elec = 1 : numel( EEG_raw.label )
            temp = EEG_raw.trial{1}(elec, [ Start:a b:End ] );
            temp_interp = interp1( EEG_raw.time{1}( [ Start:a b:End ] ), temp, EEG_raw.time{1}(a:b), 'spline');
            EEG_raw.trial{1}(elec, a:b ) = temp_interp;
         end
         flag(k) = 1;
      elseif artifact(k, 1) >= find( EEG_raw.time{1} > T_End, 1, 'first' )
         flag(k) = 1;
      end
   end
   % any unsolved issue remains for further removal of both EEG and fMRI data
   artifact(flag == 1, :) = [];
   
   % anything after 600s is removed
   temp = find( EEG_raw.time{1} > T_End, 1, 'first' );
   EEG_raw.trial{1}(:, temp:end) = [];
   EEG_raw.time{1}(:, temp:end) = [];
   EEG_raw.sampleinfo = [1, length(EEG_raw.time{1})];
   clear a b L Start End temp temp_interp flag k
   
   % filtering
   cfg = [];
   cfg.reref         = 'no';
   cfg.lpfilter    = 'yes'; cfg.lpfreq      = 90; cfg.lpfiltord     =  4;
   cfg.bsfilter    = 'yes'; cfg.bsfreq      =  [49,51]; cfg.bsfiltord     =  4;
   if strcmp(subject, 'P07')
      %                 cfg.bsfreq      =  [49,51; 24 27; 37 39];
   end
   cfg.hpfilter    = 'yes'; cfg.hpfreq      = 0.5; cfg.hpfiltord     =  6;
   cfg.demean='yes'; cfg.detrend='yes';
   EEG_raw = ft_preprocessing( cfg , EEG_raw);
   
   % visualize PSD
   figure('units', 'normalized', 'outerposition', [0.1 0.1 0.4 0.5]);
   pxx_total = [];
   for elec = 1 : numel( EEG_raw.label )
      %                 [pxx f] = cpsd(EEG.trial{1}(elec, :), EEG.trial{1}(elec, :), hann(100 * EEG_raw.fsample), 50 * EEG_raw.fsample, 512, EEG_raw.fsample);
      %                 [pxx, f] = pwelch(EEG.trial{1}(elec, :), hann(50 * EEG_raw.fsample), 25* EEG_raw.fsample, 512, EEG_raw.fsample);
      pxx = pwelch( EEG_raw.trial{1}(elec, :), hamming(0.8 * EEG_raw.time{1}(end) * EEG_raw.fsample), 20 * EEG_raw.fsample, 256, EEG_raw.fsample);
      pxx_total = [pxx_total; pxx'];
   end
   pxx_total = nanmean( pxx_total, 1);
   plot( linspace(0, 0.5 * EEG_raw.fsample, length( pxx_total )), 10 * log10( pxx_total ) );
   
   % visually inspect the data + mark to-be-rejected time ranges
   cfg = [];
   cfg.continuous = 'yes';
   cfg.channel = non_ictal_electrodes_EEG;
%    cfg.channel = 'all';
   cfg.viewmode = 'vertical';
   cfg.blocksize = 20;
   cfg.verticalpadding = 0.1;
   cfg.fontsize = 14;
   cfg = ft_databrowser(cfg, EEG_raw);
   Bad_times_EEG= [ cfg.artfctdef.visual.artifact; artifact ];
   
   % rereference data
   if strcmp(Reref, 'CmnAvg')
      cfg = [];
      cfg.reref         = 'yes';
      cfg.refchannel = 'all';
      EEG_raw = ft_preprocessing( cfg , EEG_raw);
   elseif strcmp(Reref, 'ElecTypeBased')
      EEG_raw = EEG_rereference_categories( EEG_raw );
   end
   % save the cleaned EEG
   EEG = EEG_raw; clear EEG_raw
   cd( path_data );
   save( sprintf('%s_EEG_marked_%s_%s.mat', subject, mode, Reref), 'EEG' );
%    load( sprintf('%s_Bad_times_EEG.mat', subject) )
   save( sprintf('%s_Bad_times_EEG.mat', subject), 'Bad_times_EEG' )
end

%% clean fMRI data (subject P05 8mm data is actually 4mm because 8mm was corrupted.)
mode = 'NIZ-only'
Reref = 'CmnAvg'
i_subs = [1:9];

regress_out_motion = 0;

HRF = 1;
if HRF == 1
   HRF_txt = '_forHRFanalysis';
elseif HRF == 0
   HRF_txt = '_shifted6s';
end

T_End = 600;

for i_sub = i_subs
   subject = Subjects{i_sub};
   % --------------- load fMRI + regress out motion
   %         cd( [path_fMRI '\' subject] );
   %         load(sprintf('%s_fMRI.mat', subject));
   
   cd( ['..' subject] );
   if i_sub == 5 || i_sub == 1
      load( 'RoiMinusCov.mat' );
      fMRI = RoiMinusCov;
   end
   load( 'RoiMinusCov.txt' );
   fMRI = RoiMinusCov;
   fMRI = fMRI';

   
   % shift fMRI for 6 seconds if not for the HRF analysis version
   if HRF == 0
      fMRI = fMRI(:, 3:end);
   end
   
   fMRI_time = [ 0:(size(fMRI, 2)) - 1 ] * TR;  
   
   % rematch electrode orders
   [fMRI, fMRI_time] = rematch_electrode_orders(subject, fMRI, fMRI_time);
   
   % load movement measurements
   if regress_out_motion
      cd( [path_movement '\' subject] )
      temp = dir(); temp = temp(end).name;
      movement = dlmread( temp);
      movement = movement'; movement = movement(:, 1 : numel( fMRI_time ) );
      % regress out movement from fMRI
      cd( [path_fMRI '\' subject] );
   else
      fMRI_regout = fMRI;
   end
   
   
   % load elec_info to mark bad electrodes in fMRI
   cd( [path_EEG '\' subject '\Preprocessed'] );
   load( sprintf('%s_EEG_marked_%s_%s.mat', subject, mode, Reref) );
   non_ictal_electrodes_EEG = EEG.label;
   
   
   % load original EEG to find out original number of electrodes in the excel file
   cd( [path_EEG '\' subject] )
   temp = dir();
   for k = 1:numel(temp)
      txt = temp(k).name;
      if endsWith(txt, '.set') || endsWith(txt, '.vhdr')
         break
      end
   end
   
   
   % load raw EEG (the data is already filtered)
   cfg = [];
   cfg.dataset = txt;
   cfg.channels = 'all';
   EEG_raw = ft_preprocessing(cfg);
   
   % load elec info
   cd( path_electrodes )
   temp = sprintf( 'G%d', numel(EEG_raw.label) + 1 );
   [~, ~, elec_info] = xlsread('electrodelabels.xlsx', subject, [ 'A2:' temp ]);
   temp = cellfun(@(x) find( strcmp(x, elec_info(:, 1)) ), non_ictal_electrodes_EEG);
   non_ictal_electrodes_fMRI = temp;
   if i_sub == 5
      non_ictal_electrodes_fMRI = non_ictal_electrodes_fMRI - ( non_ictal_electrodes_fMRI > 76);
   end
   fMRI_regout = fMRI_regout(non_ictal_electrodes_fMRI, :);
   clear temp;
   

   % save the regout fMRI file
   elec_info = elec_info( non_ictal_electrodes_fMRI, :);
   cd( [path_fMRI '\' subject] );
   
   save(sprintf('%s_fMRI_regout_8mm_%s%s.mat', subject, mode, HRF_txt), 'fMRI_regout');
%    save( sprintf('%s_Bad_times_fMRI_8mm%s.mat', subject, HRF_txt), 'Bad_times_fMRI' );
%    save( sprintf('%s_elec_info_8mm_%s.mat', subject, mode), 'elec_info' );
end

%% extract MNI locations
for i_sub=i_subs
   subject=Subjects{i_sub};
   cd( [path_electrodes '\coords\'] );
   MNI = textread( sprintf('roi_%s.txt', subject), '%s', 'delimiter', ' ' );
   for i = size(MNI, 1): -1 : 1
      if mod( i, 6 ) < 3 || mod( i, 6 ) >5
         MNI( i ) = [];
      end
   end
   MNI = reshape( MNI, [], 3 );
   MNI = cellfun(@(x) str2num(x), MNI );
end

%% Functions
function [fMRI_regout, R_before, R_after] = regout_movement(fMRI, movement)
% regress out movement from fMRI on every ROI
fMRI_regout = nan( size( fMRI ) );
R_before = zeros( size(fMRI, 1), size(movement, 1) );
R_after = R_before;
for ROI = 1:size(fMRI, 1)
   signal = fMRI(ROI, :);
   % identify
   for Direction = 1:size(movement, 1)
      [R_before(ROI, Direction), m(Direction), b(Direction)]=...
         regression(signal, movement(Direction, :));
   end
   % subtract
   for Direction = 1:size(movement, 1)
      signal = signal - ( m(Direction) .* movement(Direction, :) + b(Direction) * ones( size ( signal ) ) );
   end
   % estimate correlations with movement
   for Direction = 1:size(movement, 1)
      [R_after(ROI, Direction), ~, ~] = regression(signal, movement(Direction, :));
   end
   % assign to the new regressed out data
   fMRI_regout (ROI, :) = signal;
end
end
