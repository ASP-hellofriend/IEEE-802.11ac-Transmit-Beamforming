NumTxAnts = 4;  % Number of transmit antennas
NumSTS = 2;     % Number of space-time streams
NumRxAnts = 2;  % Number of receive antennas

cfgVHT = wlanVHTConfig;
cfgVHT.ChannelBandwidth = 'CBW20';
cfgVHT.APEPLength = 4000;
cfgVHT.NumTransmitAntennas = NumTxAnts;
cfgVHT.NumSpaceTimeStreams = NumSTS;
cfgVHT.MCS = 4; % 16-QAM, rate 3/4

tgacChannel = wlanTGacChannel;
tgacChannel.DelayProfile = 'Model-B';
tgacChannel.ChannelBandwidth = cfgVHT.ChannelBandwidth;
tgacChannel.SampleRate = wlanSampleRate(cfgVHT);
tgacChannel.NumReceiveAntennas = NumRxAnts;
tgacChannel.NumTransmitAntennas = NumTxAnts;
tgacChannel.TransmitReceiveDistance = 100; % Meters
tgacChannel.RandomStream = 'mt19937ar with seed';
tgacChannel.Seed = 70; % Seed to allow repeatability

noisePower = -37; % dBW

% Indices for extracting fields
ind = wlanFieldIndices(cfgVHT);

% AWGN channel to add noise with a specified noise power. The random
% process controlling noise generation is seeded to allow repeatability.
awgnChannel = comm.AWGNChannel;
awgnChannel.RandomStream = 'mt19937ar with seed';
awgnChannel.Seed = 5;
awgnChannel.NoiseMethod = 'Variance';
awgnChannel.Variance = 10^(noisePower/10);

% Calculate the expected noise variance after OFDM demodulation
noiseVar = vhtBeamformingNoiseVariance(noisePower,cfgVHT);

% Number of spatial streams
Nss = NumSTS/(cfgVHT.STBC+1);

% Get the number of occupied subcarriers in VHT fields
ofdmInfo = wlanVHTOFDMInfo('VHT-Data',cfgVHT);
Nst = ofdmInfo.NumTones;

% Generate a random PSDU which will be transmitted
rng(0); % Set random state for repeatability
psdu = randi([0 1],cfgVHT.PSDULength*8,1);

% Configure a spatial expansion transmission
vhtSE = cfgVHT;
vhtSE.SpatialMapping = 'Custom'; % Use custom spatial expansion matrix
vhtSE.SpatialMappingMatrix = helperSpatialExpansionMatrix(vhtSE);

% Generate waveform
tx = wlanWaveformGenerator(psdu,vhtSE);

% Pass waveform through a fading channel and add noise. Trailing zeros
% are added to allow for channel filter delay.
rx = tgacChannel([tx; zeros(15,NumTxAnts)]);
% Allow same channel realization to be used subsequently
reset(tgacChannel);
rx = awgnChannel(rx);
% Allow same noise realization to be used subsequently
reset(awgnChannel);

% Estimate symbol timing
tOff = wlanSymbolTimingEstimate(rx(ind.LSTF(1):ind.LSIG(2),:),vhtSE.ChannelBandwidth);

% Channel estimation
vhtltf = rx(tOff+(ind.VHTLTF(1):ind.VHTLTF(2)),:);
vhtltfDemod = wlanVHTLTFDemodulate(vhtltf,vhtSE);
chanEstSE = wlanVHTLTFChannelEstimate(vhtltfDemod,vhtSE);

% Demodulate and equalize the data
vhtdata = rx(tOff+(ind.VHTData(1):ind.VHTData(2)),:);
[~,~,symSE] = wlanVHTDataRecover(vhtdata,chanEstSE,noiseVar,vhtSE,...
    'PilotPhaseTracking','None');

refSym = wlanReferenceSymbols(cfgVHT); % Reference constellation
seConst = vhtBeamformingPlotConstellation(symSE,refSym, ...
    'Spatial Expansion Transmission Equalized Symbols');

disp('Mean received channel power per space-time stream with spatial expansion: ')
for i = 1:NumSTS
    fprintf('  Space-time stream %d: %2.2f W\n',i, ...
        sum(mean(chanEstSE(:,i,:).*conj(chanEstSE(:,i,:)),1),3))
end

% Configure a sounding packet
vhtSound = cfgVHT;
vhtSound.APEPLength = 0; % NDP so no data
vhtSound.NumSpaceTimeStreams = NumTxAnts;
vhtSound.SpatialMapping = 'Direct'; % Each TxAnt carries a STS

% Generate sounding waveform
soundingPSDU = [];
tx = wlanWaveformGenerator(soundingPSDU,vhtSound);

% Pass sounding waveform through the channel and add noise. Trailing zeros
% are added to allow for channel filter delay.
rx = tgacChannel([tx; zeros(15,NumTxAnts)]);
% Allow same channel realization to be used subsequently
reset(tgacChannel);
rx = awgnChannel(rx);
% Allow same noise realization to be used subsequently
reset(awgnChannel);

% Estimate symbol timing
tOff = wlanSymbolTimingEstimate(rx(ind.LSTF(1):ind.LSIG(2),:),vhtSound.ChannelBandwidth);

% Channel estimation
vhtLLTFInd = wlanFieldIndices(vhtSound,'VHT-LTF');
vhtltf = rx(tOff+(vhtLLTFInd(1):vhtLLTFInd(2)),:);
vhtltfDemod = wlanVHTLTFDemodulate(vhtltf,vhtSound);
chanEstSound = wlanVHTLTFChannelEstimate(vhtltfDemod,vhtSound);

% Remove impact of cyclic shift from channel estimate
chanEstSound = vhtBeamformingRemoveCSD(chanEstSound, ...
    vhtSound.ChannelBandwidth,vhtSound.NumSpaceTimeStreams);

chanEstPerm = permute(chanEstSound,[3 2 1]); % permute to Nr-by-Nt-by-Nst
V = zeros(Nst,NumTxAnts,NumRxAnts);
for i = 1:Nst
    [U,S,V(i,:,:)] = svd(chanEstPerm(:,:,i),'econ');
end
steeringMatrix = V(:,:,1:NumSTS); % Nst-by-Nt-by-Nsts

% Configure a transmission with beamforming
vhtBF = cfgVHT;
vhtBF.SpatialMapping = 'Custom';
% Permute steering matrix to Nst-by-Nsts-by-Nt
vhtBF.SpatialMappingMatrix = permute(steeringMatrix,[1 3 2]);

% Generate beamformed data transmission
tx = wlanWaveformGenerator(psdu,vhtBF);

% Pass through the channel and add noise. Trailing zeros
% are added to allow for channel filter delay.
rx = tgacChannel([tx; zeros(15,NumTxAnts)]);
rx = awgnChannel(rx);

% Estimate symbol timing
tOff = wlanSymbolTimingEstimate(rx(ind.LSTF(1):ind.LSIG(2),:),vhtBF.ChannelBandwidth);

% Channel estimation
vhtltf = rx(tOff+(ind.VHTLTF(1):ind.VHTLTF(2)),:);
vhtltfDemod = wlanVHTLTFDemodulate(vhtltf,vhtBF);
chanEstBF = wlanVHTLTFChannelEstimate(vhtltfDemod,vhtBF);

% Demodulate and equalize the data
vhtdata = rx(tOff+(ind.VHTData(1):ind.VHTData(2)),:);
[~,~,symBF] = wlanVHTDataRecover(vhtdata,chanEstBF,noiseVar,vhtBF,...
    'PilotPhaseTracking','None');

bfConst = vhtBeamformingPlotConstellation(symBF,refSym, ...
    'Beamformed Transmission Equalized Symbols');

disp('Mean received channel power per space-time stream with SVD transmit beamforming: ')
for i = 1:NumSTS
    fprintf('  Space-time stream %d: %2.2f W\n',i, ...
        sum(mean(chanEstBF(:,i,:).*conj(chanEstBF(:,i,:)),1),3))
end

str = sprintf('%dx%d',NumTxAnts,NumRxAnts);
compConst = vhtBeamformingPlotConstellation([symSE(:) symBF(:)],refSym, ...
    'Beamformed Transmission Equalized Symbols', ...
    {[str ' Spatial Expansion'],[str ' Transmit Beamforming']});

EVM = comm.EVM;
EVM.AveragingDimensions = [1 2]; % Average over all subcarriers and symbols
EVM.MaximumEVMOutputPort = true;
EVM.ReferenceSignalSource  = 'Estimated from reference constellation';
EVM.ReferenceConstellation = refSym;

[rmsEVMSE,maxEVMSE] = EVM(symSE); % EVM using spatial expansion
[rmsEVMBF,maxEVMBF] = EVM(symBF); % EVM using beamforming

for i = 1:Nss
    fprintf(['Spatial stream %d EVM:\n' ...
        '  Spatial expansion:    %2.1f%% RMS, %2.1f%% max\n' ...
        '  Transmit beamforming: %2.1f%% RMS, %2.1f%% max\n'], ...
        i,rmsEVMSE(i),maxEVMSE(i),rmsEVMBF(i),maxEVMBF(i));
end