%{
------------------------------
This function decomposes the input data (EEG_orig) into its ICs using AMICA
algorithm using EEGLAB. Fieltrip is needed to visualize the ICs.
timecourse and PSD of ICs are visualized, based on which the user can
choose bad ICs. The user will insert the label of bad components into a
dialogue box (as numbers, separated by space).

Output is the clean version of EEG_Orig.
------------------------------
%}
function [EEG_clean] = EEG_ICA_VisualCleaning(EEG_orig, EEGLAB_path, Fieldtrip_path)
% initialize
addpath( EEGLAB_path ); eeglab;
addpath( Fieldtrip_path ); ft_defaults;
data = EEG_orig.trial{1};
Fs = EEG_orig.fsample;

% run ICA
[ weights_blinkdat, sphere_blinkdat, ~ ] = runamica15(data, 'max_threads', 6);

% extract deweighting and weighting matrices (respectively)
A = weights_blinkdat * sphere_blinkdat;
A_inv = pinv(A);

% create IC data
act_blinkdata = A * data;
IC_labels = num2cell( 1 : size(act_blinkdata, 1) )';
IC_labels = cellfun(@(x) num2str(x), IC_labels, 'UniformOutput', 0);

% convert the IC data to the fieldtrip EEG format to visualize it
cfg = [];
cfg.channel = EEG_orig.label(1 : size(act_blinkdata, 1));
EEG_sample = ft_selectdata(cfg, EEG_orig);
EEG_sample.trial{1} = act_blinkdata;
EEG_sample.label = IC_labels;

% plot components
cfg = [];
cfg.continuous = 'yes';
cfg.channel = 'all';
cfg.viewmode = 'vertical';
cfg.blocksize = 20;
cfg.verticalpadding = 0.1;
cfg.fontsize = 14;
cfg.bandwidth = 6;
cfg = ft_databrowser(cfg, EEG_sample);

% plot PSD of components
figure('units', 'normalized', 'outerposition', [0 0 1 1])
set(gcf, 'color', 'w', 'inverthardcopy', 'off')
for IC_index = 1 : size(A,1)
   subplot(4, 5, IC_index);
   hold on;
   [pxx, w] = pwelch(squeeze( act_blinkdata(IC_index, :) ), 10 * Fs, 5 * Fs, 256);
   plot(w / (2*pi) * Fs, 10 * log10( pxx ));
   set(gca, 'fontsize', 6); xlabel('Hz'); grid minor
   line([10 10], ylim, 'color', 'r', 'linewidth', 1, 'linestyle', '-');
   xlim([0 Fs / 2])
   title(num2str(IC_index))
   ylabel('10*log10(Power)')
   xlabel('Freq (Hz)')
end

% identify bad ICs (by visual inspection)
Answer = inputdlg('Enter the number of bad ICs to be removed');
Answer = str2num( Answer{1} );
Comps_to_reject = zeros(size(A, 1), 1);
Comps_to_reject( Answer ) = 1;
Comps_to_reject = boolean( Comps_to_reject );

% remove bad ICs
data_out = A_inv(:, ~Comps_to_reject) * act_blinkdata(~Comps_to_reject, :);

% assign the clean data to the output
EEG_clean = EEG_orig;
EEG_clean.trial{1} = data_out;
end
% END OF FILE

