function [conn_Amp, conn_Amp_static, conn_Amp_static_RegOut] = ...
   AmpC_or_PhC_eFC(Measure, EEG_bp, EEG, EEG_time, EEG_time_ds, freq, T,...
   dist_electrodes, canonical_hrf)

HRF_flag = 0;
if ~isempty(canonical_hrf)
   HRF_flag = 1;
end


numelec = size(EEG, 1);
Fs = 1 / (EEG_time(2) - EEG_time(1) );


% bandpass data
switch freq
   case 1
      Freqrange=[1 4];
   case 2
      Freqrange=[5 7];
   case 3
      Freqrange=[8 13];
   case 4
      Freqrange=[14 30];
   case 5
      Freqrange=[31 60];
   case 6
      Freqrange=[61 110];
end



% band-pass filter
if isempty(EEG_bp)
   EEG_bp = nan( size( EEG ) );
   [B,A] = cheby2(4, 40, Freqrange / Fs*2, 'bandpass');
   
   for elec = 1 : numelec
      temp = squeeze( EEG(elec, :) );
      EEG_bp(elec, :) = filtfilt(B, A, temp); EEG_bp(elec, :) = EEG_bp(elec, :) / max( abs( EEG_bp(elec, :) ) ) * max( abs( EEG(elec, :) ) );
   end
end




% extract amplitude or phase of the band-limited signal
if strcmp(Measure, 'AmpC')
   EEG_bp_env = nan(size( EEG, 1), length( EEG_time ));
elseif strcmp(Measure, 'PhC')
   EEG_bp_env = nan(size( EEG, 1), length( EEG_time ));
end

for elec = 1 : numelec
   if strcmp(Measure, 'AmpC')
      EEG_bp_env(elec, :) = abs( hilbert( EEG_bp(elec, :) ) );
      
      if HRF_flag
         temp = conv(EEG_bp_env(elec, :), canonical_hrf, 'full');
         temp = temp(1: size(EEG_bp_env, 2));
         EEG_bp_env(elec, :) = temp;
      end

   elseif strcmp(Measure, 'PhC')
      Signal_Phase = unwrap( angle( hilbert( EEG_bp(elec, :) ) ) );
      EEG_bp_env(elec, :) = Signal_Phase;
   end
end


% estimate FC
conn_Amp = nan(numelec, numelec, size(EEG_time_ds, 2));
h = waitbar(0, sprintf('Calculating eFC (freq%d)...', freq) ); counter=0;
if strcmp(Measure, 'AmpC')
   EEG_bp_env_zscore = zscore( EEG_bp_env' )';
   
   for i = 1 : numelec
      for j = 1 : numelec
         if i < j
            counter=counter + 1;
            waitbar( counter / ( numelec * (numelec-1) / 2) );
            % estimate eFC
            AmpC_edge = nan( size( EEG_time_ds ) );
            AmpC_dot = EEG_bp_env_zscore(i, :) .* EEG_bp_env_zscore(j, :);
            for time  = 1 : numel( EEG_time_ds )
               a = EEG_time_ds(time); a = max(a, 0);
               b = EEG_time_ds(time) + T; b = min(b, EEG_time(end));
               AmpC_dot_interval = AmpC_dot( EEG_time >= a & EEG_time < b );
               AmpC_edge(time) = nanmean( (AmpC_dot_interval) );
            end
            conn_Amp(i,j,:) = AmpC_edge;
            conn_Amp(j,i,:) = AmpC_edge;
         end
      end
      conn_Amp(i,i,:) = nan;
   end
   
   
elseif strcmp(Measure, 'PhC')
   for i = 1 : numelec
      for j = 1 : numelec
         if i < j
            PhC_edge = nan( size( EEG_time_ds ) );
            counter=counter + 1;
            waitbar( counter / ( numelec * (numelec-1) / 2) );
            Phase_diff = EEG_bp_env(i, :) - EEG_bp_env(j, :);
            for time  = 1 : numel( EEG_time_ds )
               a = EEG_time_ds(time); a = max(a, 0);
               b = EEG_time_ds(time) + T; b = min(b, EEG_time(end));
               Phase_diff_interval = Phase_diff( EEG_time >= a & EEG_time < b );
               PhC_edge(time) = abs( nanmean( exp( 1i*(Phase_diff_interval) ) ) );
            end
            conn_Amp(i,j,:) = PhC_edge;
            conn_Amp(j,i,:) = PhC_edge;
         end
      end
      conn_Amp(i,i,:) = nan;
   end

end
close(h)

% static FC
conn_Amp_static = nanmean( conn_Amp, 3);

if ~isempty( dist_electrodes )
   [conn_Amp_static_RegOut,~,~,~] = Dist_Reg_Out(conn_Amp_static, dist_electrodes);
else
   conn_Amp_static_RegOut = [];
end
end

