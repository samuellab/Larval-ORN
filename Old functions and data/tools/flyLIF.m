function out=flyLIF_GS(params,connect,synType, inputs,inputVs)
% params, parameter about the LIF model and simulaiton.
% synType, synapse type, 1 means ACh, 2 means Glu.
% connect, connectivity matrix.
% inputs, input neuron index.
% inputVs, input neuron's time course membrane potential

numNeurons=size(connect,1); 
numInputs=length(inputs);

if numInputs ~=size(inputVs,2)
    error('input current count mismatch');
end
if length(synType) ~= numNeurons
    error('input current count mismatch');
end

if isfield(params,'time')
    time=params.time; else time=6.8; end                   % 0.5 (or 0.083s) seconds
if isfield(params,'dt')
    dt=params.dt; else dt=1e-4; end                        % 1e-4 s
if isfield(params,'V0')
    V0=params.V0; else V0=-52e-3; end                      % -52e-3 (or -40e-3, -60e-3) V
if isfield(params,'spikeMax')
    spikeMax=params.spikeMax; else spikeMax=20e-3; end     % 20e-3 (or -90e-3) V
if isfield(params,'spikeMin')
    spikeMin=params.spikeMin; else spikeMin=-72e-3; end    % -72e-3 (or -90e-3) V
if isfield(params,'spikeThr')
    spikeThr=params.spikeThr; else spikeThr=-45e-3; end    % -45e-3 (or -30e-3) V
if isfield(params,'Cmem')
    Cmem=params.Cmem; else Cmem=2e-10; end                  % 2e-8 F
if isfield(params,'Rmem')
    Rmem=params.Rmem; else Rmem=1e8; end                   % 1e6 Ohm
if isfield(params,'apDuration')
    apDuration=params.apDuration; else apDuration=2; end   % 2 (or 0.5, 8) ms
if isfield(params,'psgTimeFactor')
    psgTimeFactor=params.psgTimeFactor; else psgTimeFactor = 10; end % 10, 10 times of Tau
if isfield(params,'psgTauACh')
    psgTauACh=params.psgTauACh; else psgTauACh = 3e-3; end 	% 3ms, Acetylcholine synapse
if isfield(params,'gSynACh')
    gSynACh=params.gSynACh; else gSynACh = 3e-10; end      % S, Acetylcholine synapse
if isfield(params,'psgTauGlu')
    psgTauGlu=params.psgTauGlu; else psgTauGlu = 30e-3; end % 25ms, Glutamate
if isfield(params,'gSynGlu')
    gSynGlu=params.gSynGlu; else gSynGlu = 7.5e-11; end      %7.5e-11 S, Glutamate
if isfield(params,'VrevE')
    VrevE=params.VrevE; else VrevE = 0; end             % 0V,excitatory synaptic reversal potential
if isfield(params,'VrevI')
    VrevI=params.VrevI; else VrevI = -70e-3; end           % -70mV, inhibitory synaptic reversal potential

N=time/dt;

%% define the action potential
apRiseTime=round((apDuration*10)/2);
apRise=normpdf(-1:1/apRiseTime:0);
apRise=(apRise-min(apRise))/(max(apRise)-min(apRise));
apRise=apRise*(spikeMax-spikeThr)+spikeThr;
apFallTime=round((apDuration*9)/2);
apFall=sin(pi()/2:pi()/apFallTime:3*pi()/2);
apFall=(apFall-min(apFall))/(max(apFall)-min(apFall));
apFall=apFall*(spikeMax-spikeMin)+spikeMin-.0001;
AP=[apRise apFall]';
AP(1)=[];
AP(end+1)=AP(end)+0.0001;
AP(end+1)=AP(end)+0.0001;
apLength=length(AP);

%% define post synaptic current
PSGAChTime = (1 : psgTauACh * psgTimeFactor/dt)*dt ;
PSGACh = gSynACh * (PSGAChTime /psgTauACh) .* exp(1-PSGAChTime./psgTauACh);

PSGGluTime = (1 : psgTauGlu * psgTimeFactor/dt)*dt ;
PSGGlu = gSynGlu * (PSGGluTime /psgTauGlu) .* exp(1-PSGGluTime./psgTauGlu);

%% prep for the simulation 
tBuffACh = length(PSGACh);
tBuff=length(PSGGlu);
V=zeros(N+tBuff,numNeurons)+V0;
V(1,:)=V0;
V(1:N, inputs) = inputVs;

APMask=zeros(N+tBuff,numNeurons);
gSyn=zeros(N+tBuff,numNeurons);

% define the rev voltage, based on the synapse type
VerN = zeros(1, numNeurons);
for i =1:numNeurons
    if synType(i) == 1
        VerN(i) = VrevE;
    elseif synType(i) == 2
        VerN(i) = VrevI;
    end
end

%% calculate ODE
for t=2:N
    for i=1:numNeurons
        if APMask(t,i)==0
            if V(t-1,i)>spikeThr
                V(t-1,i)=spikeThr;
                V(t:t+apLength-1,i)=AP;
                APMask(t:t+apLength-1,i)=1;
                if synType(i) == 1 %ACh
                	gSyn(t:t+tBuffACh-1,i)=gSyn(t:t+tBuffACh-1,i)+PSGACh';
                elseif synType(i) == 2 %Glu
                    gSyn(t:t+tBuff-1,i)=gSyn(t:t+tBuff-1,i)+PSGGlu';
                end
            else
                if ~any(i == inputs) %in this version, only update non-input nodes
                    I_syn = - sum(connect(:,i) .* gSyn(t-1, :)' .* (V(t-1, i) - VerN)');
                    dV=(1/Cmem)*( (V0-V(t-1,i))/Rmem + I_syn );
                    V(t,i)=V(t-1,i)+dV*dt;
                end
            end
        end
    end
end

%% define output
out.V=V(1:N, :); %get rid of buff
out.gSyn=gSyn(1:N, :); %get rid of buff
out.inputPSCs=inputVs;
out.APMask=APMask;
out.connect=connect;
out.P=params;
out.P.time=time;
out.P.dt=dt;
out.P.V0=V0;
out.P.spikeMax=spikeMax;
out.P.spikeMin=spikeMin;
out.P.spikeThr=spikeThr;
out.P.Cmem=Cmem;
out.P.Rmem=Rmem;
out.P.apDuration=apDuration;
out.P.psgTauACh=psgTauACh;
out.P.psgTimeFactor=psgTimeFactor;
out.P.psgTauACh=psgTauACh;
out.P.gSynACh=gSynACh;
out.P.psgTauGlu=psgTauGlu;
out.P.gSynGlu=gSynGlu;
out.P.VrevE=VrevE;
out.P.VrevI=VrevI;