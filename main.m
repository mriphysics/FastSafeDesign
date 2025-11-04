%	Safe Spiral Out design
%	Edited by David Leitao 2025
clearvars; close all; clc

addpath(genpath('C:\Users\david\Documents\Projects\Spiral\SAFE_PNS'))

% Conversion rates:
% mT/m/ms --[ x1e1 ]--> G/m/ms --[ x1e-2 ]--> G/cm/ms --[ x1e3 ]--> G/cm/s
% mT/m/ms --------------------------[ x1e2 ]----------------------> G/cm/s
%
% mT/m   --[ x1e-2 ]--> mT/cm --[ x1e1 ]--> G/cm
% mT/m   ------------[ x1e-1 ]------------> G/cm#

%% Define SAFE model, forbidden bands and hardware limits

%%% Load SAFE model from gradient coil file (Siemens)
hw = safe_hw_from_asc('C:\Users\david\Documents\Projects\Spiral\SAFE_PNS\MP_GradSys_K2298_2250V_1250A_W60_SC72CD.asc',1,1);
% Check if hardware parameters are consistent
safe_hw_check(hw)
% Check if this hw is part of the library (validate hw)
safe_hw_verify(hw);

% %%% Alternartively use fictional SAFE model (from pulseq; NOT TESTED)
% hw.x.tau1 = 0.20; %[ms]
% hw.x.tau2 = 0.03; %[ms]
% hw.x.tau3 = 3.00; %[ms]
% hw.x.a1 = 0.40;
% hw.x.a2 = 0.10;
% hw.x.a3 = 0.50;
% hw.x.stim_limit = 30.0;  %[T/m/s]
% hw.x.stim_thresh = 24.0; %[T/m/s]
% hw.x.g_scale = 0.35;
% 
% hw.y.tau1 = 1.50; %[ms]
% hw.y.tau2 = 2.50; %[ms]
% hw.y.tau3 = 0.15; %[ms]
% hw.y.a1 = 0.55;
% hw.y.a2 = 0.15;
% hw.y.a3 = 0.30;
% hw.y.stim_limit = 15.0; %[T/m/s]
% hw.y.stim_thresh = 12.0; %[T/m/s]
% hw.y.g_scale = 0.31;
% 
% hw.z.tau1 = 0.00; %[ms]
% hw.z.tau2 = 0.00; %[ms]
% hw.z.tau3 = 0.00; %[ms]
% hw.z.a1 = 0.0;
% hw.z.a2 = 0.0;
% hw.z.a3 = 0.0;
% hw.z.stim_limit = 15.0; %[T/m/s]
% hw.z.stim_thresh = 12.0; %[T/m/s]
% hw.z.g_scale = 0.31;

% %%% Safe model can be made rotational invariant by taking the worst case scenario across the 3 axes
% sys.safeModel.RIV = true;
% sys.safeModel.tauW       = [1e-3*max([hw.x.tau1, hw.y.tau1, hw.z.tau1]), ...
%                             1e-3*max([hw.x.tau2, hw.y.tau2, hw.z.tau2]), ...
%                             1e-3*max([hw.x.tau3, hw.y.tau3, hw.z.tau3])];
% sys.safeModel.AW         = [max([hw.x.a1, hw.y.a1, hw.z.a1]), ...
%                             max([hw.x.a2, hw.y.a2, hw.z.a2]), ...
%                             max([hw.x.a3, hw.y.a3, hw.z.a3])];
% sys.safeModel.pnsScaling = max([hw.x.g_scale/hw.x.stim_limit, hw.y.g_scale/hw.y.stim_limit, hw.z.g_scale/hw.z.stim_limit]);

sys.safeModel.tauX       = 1e-3*[hw.x.tau1, hw.x.tau2, hw.x.tau3];
sys.safeModel.tauY       = 1e-3*[hw.y.tau1, hw.y.tau2, hw.y.tau3];
sys.safeModel.tauZ       = 1e-3*[hw.z.tau1, hw.z.tau2, hw.z.tau3];
sys.safeModel.AX         = [hw.x.a1, hw.x.a2, hw.x.a3];
sys.safeModel.AY         = [hw.y.a1, hw.y.a2, hw.y.a3];
sys.safeModel.AZ         = [hw.z.a1, hw.z.a2, hw.z.a3];
sys.safeModel.pnsScaling = [hw.x.g_scale/hw.x.stim_limit, hw.y.g_scale/hw.y.stim_limit, hw.z.g_scale/hw.z.stim_limit];


%%% Forbidden frequency bands (each row contains lower and upper limits of a band)
sys.resonFreq  = [1000  1200 
                  600   650
                  300   350];

%%% Hardware max gradient and slew rate (here depends on mode selected)
sys.gradientMode = 'max';
switch sys.gradientMode
    case 'max'
        sys.smax           = 250e2;	                 % maxium slew rate [G/cm/s]
        sys.gmax           = 135e-1;                 % maximum gradient amplitude [G/cm]
    case 'rapid'
        sys.smax           = 200e2;	                 % maxium slew rate [G/cm/s]
        sys.gmax           = 42e-1;                  % maximum gradient amplitude [G/cm]
    case 'normal'
        sys.smax           = 200e2;	                 % maxium slew rate [G/cm/s]
        sys.gmax           = 40e-1;                  % maximum gradient amplitude [G/cm]
    case 'slow'
        sys.smax           = 100e2;	                 % maxium slew rate [G/cm/s]
        sys.gmax           = 22e-1;                  % maximum gradient amplitude [G/cm]
end

Tgrad          = 10e-6;                  % gradient raster time [seconds]


%% Design safe spiral out

%%% Spiral design parameters
slewMargin     = 1.0;                    % design margin for slew rate (<=1)
Tadc           = 1e-6;                   % ADC sampling time [seconds]
N              = 1;		                 % Spiral interleaves
Fcoeff         = 24;                     % FOV coefficients (see formula: FOV(r) = Sum    Fcoeff(k)*(r/rmax)^(k-1)) 
res            = 2;                      % In-plane resolution [mm]
pnsDesignLimit = 0.88;                   % PNS limit used for spiral design
NrepDesign     = 1;                      % No. of times each design is repeated (computation time estimation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

krmax = 10/(2*res); % [cm^(-1)]
sys4des = sys;

% First estimate spiral duration with unsafe design
tmpSys = rmfield(sys4des, {'safeModel','resonFreq'});
for n=1:NrepDesign
    tic
    [traj,~,~,time,~,~,~,~] = SafeSpiralOut(tmpSys,Tadc,Tgrad,N,Fcoeff,krmax,slewMargin,pnsDesignLimit);
    tdesignUnsafe(n) = toc;
end
fprintf(1,'Unsafe spiral design took %.0fms %c %.0fms\n',1e3*mean(tdesignUnsafe),char(177),std(tdesignUnsafe))

% save(sprintf('output_unsafeSpiral_%s',sys.gradientMode),'traj','tdesign')

% Expand fordibben bands
Tspiral = time(end);
df = 1/Tspiral;
% Reduce lower resonance frequency of each band
if df>200 
    df = 200; %cap reduction to 200Hz (arbitrary choice)
end
sys4des.resonFreq(:,1) = sys4des.resonFreq(:,1) - df;

% Fast safe spiral design
for n=1:NrepDesign
    tic
    [k,g,s,time,r,theta,f,pnsSpiral] = SafeSpiralOut(sys4des,Tadc,Tgrad,N,Fcoeff,krmax,slewMargin,pnsDesignLimit);
    tdesignSafe(n) = toc;
end
fprintf(1,'Safe spiral design took %.0fms %c %.0fms\n',1e3*mean(tdesignSafe),char(177),std(tdesignSafe))

% save(sprintf('output_fastSafeSpiral_%s',sys.gradientMode),'k','g','s','time','tdesign','pnsDesignLimit')


% Convert to use mT, meter and seconds
g = 1e1*g; %[mT/m]
s = (g - [0,g(1:end-1)])/Tgrad/1e3;
gamma = 42577.478518; %[Hz/mT]
k = gamma*cumsum(g)*Tgrad;
gmax = 1e1*sys.gmax; %[mT/m]
smax = 1e1*sys.smax/1e3; %[mT/m/ms]

%% Plot results

% Plot k-space, gradient and slew
figure;
set(gcf,'position',[200 500 1200 325],'color','w')
hsp(1) = subplot(2,2,1);
plot(1e3*time,real(g),'linewidth',1.5); hold on
plot(1e3*time,imag(g),'linewidth',1.5)
plot(1e3*time,abs(g),'k','linewidth',1.5)
xlabel('time (ms)'); ylabel('Gradient (mT/m)'); legend('G_x','G_y','G_{max}','location','southeast')
grid on; box on; xlim([0 1e3*time(end)])
text(0.5,1.1,'Gradient','fontweight','bold','units','normalized','HorizontalAlignment','center','fontsize',14)

hsp(2) = subplot(2,2,2);
plot(1e3*time,real(s),'linewidth',1.5); hold on
plot(1e3*time,imag(s),'linewidth',1.5)
plot(1e3*time,abs(s),'k','linewidth',1.5)
xlabel('time (ms)'); ylabel('Slew rate (mT/m/ms)'); legend('S_x','S_y','S_{max}','location','southeast')
grid on; box on; xlim([0 1e3*time(end)])
text(0.5,1.1,'Slew','fontweight','bold','units','normalized','HorizontalAlignment','center','fontsize',14)

hsp(3) = subplot(2,2,3);
plot(1e3*time,real(k),'linewidth',1.5); hold on
plot(1e3*time,imag(k),'linewidth',1.5)
plot(1e3*time,abs(k),'k','linewidth',1.5)
xlabel('time (ms)'); ylabel('K-space (1/m)'); legend('K_x','K_y','K_r','location','southeast')
grid on; box on; xlim([0 1e3*time(end)])
text(0.5,1.1,'k-space','fontweight','bold','units','normalized','HorizontalAlignment','center','fontsize',14)
for n=1:3
    hsp(n).Position = [0.06+(n-1)*0.33 0.12 0.27 0.75];
end

% Predict PNS levels (concatenate spirals with 4ms gap to compute steady-state build up)
gwf = 1e-3*cat(2,real(g(:)),imag(g(:)),zeros(numel(g),1));
rf = zeros(numel(g),1);
dt = Tgrad;
nzp = round(4e-3/Tgrad);
gwf = cat(1,gwf,zeros(nzp,3));
gwf = repmat(gwf,[10 1]);
[pns,res] = safe_gwf_to_pns(gwf, rf, dt, hw, 1);

% Plot PNS
figure;
set(gcf,'position',[25 75 700 325],'color','w')
hsp(1) = subplot(1,2,1);
plot(1e3*(1:size(pnsSpiral,2))*Tgrad,100*pnsSpiral','linewidth',1.5); hold on;
plot(1e3*(1:size(pnsSpiral,2))*Tgrad,100*sqrt(sum(pnsSpiral'.^2,2)),'k','linewidth',1.5);
ylabel('PNS level (%)'); xlabel('time (ms)')
xlim([0 size(pnsSpiral,2)*Tgrad]*1e3)
legend('PNS_x','PNS_y','PNS_{rms}','location','southeast')
title('Predicted PNS during design'); grid on
hsp(2) = subplot(1,2,2);
idx = ((size(pns,1)-size(res.zp2,1)-size(pnsSpiral,2)):(size(pns,1)-size(res.zp2,1)))-nzp;
plot((1:numel(idx))*Tgrad*1e3,pns(idx,1:2),'linewidth',1.5); hold on;
plot((1:numel(idx))*Tgrad*1e3,sqrt(sum(pns(idx,:).^2,2)),'k','linewidth',1.5)
legend('PNS_x','PNS_y','PNS_{rms}','location','southeast')
title('Steady-state PNS'); grid on
ylabel('PNS level (%)'); xlabel('time (ms)')
xlim([0 numel(idx)*Tgrad]*1e3)
for n=1:2
    hsp(n).Position = [0.08+(n-1)*0.48 0.12 0.4 0.75];
end

% Plot gradient frequency and spectrum
figure;
set(gcf,'position',[750 75 700 325],'color','w')
hsp(1) = subplot(2,1,1);
[pxx,freq] = periodogram(real(g(:)), [], 1e5, 1/Tgrad);   % default: no window overlap
[pyy,freq] = periodogram(imag(g(:)), [], 1e5, 1/Tgrad);   % default: no window overlap
plot(freq,pxx,'linewidth',1.5); hold on;
plot(freq,pyy,'linewidth',1.5)
text(0.5,1.1,'Spectrum spiral','fontweight','bold','units','normalized','HorizontalAlignment','center','fontsize',14)
grid on; xlim([0 2000])
% hsp(1).Position = [0.06 0.12 0.42 0.75];
ylims{1} = ylim;
ha(1) = gca;
hold on;
for n=1:size(sys.resonFreq,1)
    fill([sys.resonFreq(n,1) sys.resonFreq(n,2) sys.resonFreq(n,2) sys.resonFreq(n,1)],...
        [ylims{1}(1) ylims{1}(1) ylims{1}(2) ylims{1}(2)],'r','facealpha',0.2,'edgealpha',0)
end
hln = findobj(gca,'type','line');
legend(hln([2 1]),'\epsilon{}_x','\epsilon{}_y','location','northeast')
xlabel('Frequency (Hz)'); ylabel('Power spectral density (mT/m)^2/Hz')

hsp(2) = subplot(2,1,2);
histogram(f,0:10:2000)
xlabel('Frequency (Hz)'); ylabel('No. of events')
text(0.5,1.1,'"Instantenous" frequency histogram','fontweight','bold','units','normalized','HorizontalAlignment','center','fontsize',14)
grid on; xlim([0 2000])
ylims{2} = ylim;
ha(2) = gca;
hold on;
for n=1:size(sys.resonFreq,1)
    fill([sys.resonFreq(n,1) sys.resonFreq(n,2) sys.resonFreq(n,2) sys.resonFreq(n,1)],...
        [ylims{2}(1) ylims{2}(1) ylims{2}(2) ylims{2}(2)],'r','facealpha',0.2,'edgealpha',0)
end
for n=1:2
    hsp(n).Position = [0.08+(n-1)*0.48 0.12 0.4 0.75];
end
ylim(ha(1),ylims{1})
ylim(ha(2),ylims{2})