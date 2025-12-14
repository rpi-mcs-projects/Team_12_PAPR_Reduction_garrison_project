function mid_term_check
% PAPR Reduction for OFDM under PA Nonlinearity (final, presentation-ready + extras)
% Adds:
%   1) Target-based view: min BO to hit EVM target at SNR=20/25 dB
%   2) Uncertainty bars: EVM error bars across seeds
%   3) Sensitivity: SLM U=2 vs 4 and PTS B=2 vs 4 at BO=0 dB, SNR=20 dB
%
% Outputs saved to ./results

clear; clc; close all; rng(42);
set(groot,'defaultFigureVisible','on');
set(0,'DefaultFigureWindowStyle','normal');

outdir = fullfile(pwd,'results'); if ~exist(outdir,'dir'), mkdir(outdir); end

cfg  = config();
meta = build_meta(cfg, cfg.Nsym);

% ---------- Frequency-domain frames (shared across methods) ----------
[framesFD, dataIdx, trIdx, ~, txSymIdxPerSym] = gen_ofdm_frames(cfg);

% ---------- Pre-PA waveforms ----------
x_BASE = ofdm_mod(framesFD, cfg, meta);

[frames_SLM, ~, ~] = slm_reduce(framesFD, cfg.SLM);
x_SLM = ofdm_mod(frames_SLM, cfg, meta);

[frames_PTS, ~] = pts_reduce(framesFD, cfg.PTS);
x_PTS = ofdm_mod(frames_PTS, cfg, meta);

[frames_TR, ~] = tone_reservation(framesFD, trIdx, cfg.TR, cfg);
x_TR = ofdm_mod(frames_TR, cfg, meta);

% ---------- Sanity (linear channel, high SNR) ----------
cfg_chk = cfg; cfg_chk.sim.SNRdB = 100; cfg_chk.PA.backoff_dB = 1e9;
[evm_chk, ber_chk] = link_metrics(x_BASE, framesFD, dataIdx, cfg_chk, meta, txSymIdxPerSym);
fprintf('Sanity: EVM=%.4f%%, BER=%.2e (expect ~0, ~0)\n', evm_chk, ber_chk);

% ---------- CCDF(PAPR) (pre-PA) ----------
f_ccdf = figure('Name','CCDF(PAPR)'); hold on; grid on; set(gcf,'Color','w');
ccdf_plot(x_BASE,'No PAPR reduction',meta);
ccdf_plot(x_SLM, sprintf('SLM (U=%d)', cfg.SLM.U),meta);
ccdf_plot(x_PTS, sprintf('PTS (B=%d)', cfg.PTS.B),meta);
ccdf_plot(x_TR,  sprintf('TR (%.0f%% tones)',100*numel(trIdx)/cfg.Nsub),meta);
xlabel('PAPR [dB]'); ylabel('CCDF (P\{PAPR>\gamma\})');
title('CCDF of PAPR (pre-PA)'); legend('show','Location','northeast');

% 0.1% CCDF callout
q = 0.001;
methodsL = {'BASE','SLM','PTS','TR'};
wavesL   = {x_BASE, x_SLM, x_PTS, x_TR};
vals = zeros(numel(wavesL),1);
for i = 1:numel(wavesL)
    vals(i) = papr_db_symbol(wavesL{i}, meta, (1-q)*100);
end
fprintf('PAPR @0.1%% CCDF [dB]: BASE=%.2f, SLM=%.2f, PTS=%.2f, TR=%.2f\n', vals);
mcs_save_fig(f_ccdf, fullfile(outdir,'ccdf_papr.png'));

% ---------- Metric sweeps across SNR and PA back-off ----------
methods = {'BASE','SLM','PTS','TR'};
waves   = {x_BASE, x_SLM, x_PTS, x_TR};
frames  = {framesFD, frames_SLM, frames_PTS, frames_TR};

SNRs = cfg.sim.SNRdB_list;
BOs  = cfg.PA.backoff_list_dB;

nAvg = 5;  % average over 5 noise draws
R = struct();

for m = 1:numel(methods)
    name = methods{m};
    xpre = waves{m};
    frms = frames{m};

    [~, preACLR] = measure_aclr(xpre, cfg);
    prePAPR = papr_db_symbol(xpre, meta, 99.9);

    for ib = 1:numel(BOs)
        cfg_i = cfg; cfg_i.PA.backoff_dB = BOs(ib);
        y = pa_nonlinearity(xpre, cfg_i.PA);
        [~, postACLR] = measure_aclr(y, cfg_i);

        for is = 1:numel(SNRs)
            cfg_i.sim.SNRdB = SNRs(is);
            evm_runs = zeros(nAvg,1); ber_runs = zeros(nAvg,1);

            for it = 1:nAvg
                rng(1000 + 100*ib + 10*is + it);
                [evm_runs(it), ber_runs(it)] = link_metrics(y, frms, dataIdx, cfg_i, meta, txSymIdxPerSym);
            end

            evm_mean = mean(evm_runs);
            evm_std  = std(evm_runs,0,1);
            ber_mean = mean(ber_runs);
            ber_std  = std(ber_runs,0,1);

            key = keyname(name, SNRs(is), BOs(ib));
            R.(key) = struct( ...
                'method',name,'SNRdB',SNRs(is),'BOdB',BOs(ib), ...
                'PAPR_pre_dB',prePAPR,'ACLR_pre_dB',preACLR,'ACLR_post_dB',postACLR, ...
                'EVM_pct',evm_mean,'EVM_std',evm_std,'BER',ber_mean,'BER_std',ber_std);
        end
    end
end

% ---- results table & saving ----
Tab = results_to_table(R);
mcs_save_results(Tab, cfg, meta, R, outdir);

% ---- plots: BER/EVM vs SNR at each back-off (with EVM error bars) ----
mk = {'-o','-s','-^','-d'};
for ib = 1:numel(BOs)
    bo = BOs(ib);

    f1 = figure('Name',sprintf('BER_vs_SNR_BO_%ddB',bo),'Visible','on','NumberTitle','off'); 
    hold on; grid on; set(gcf,'Color','w');
    f2 = figure('Name',sprintf('EVM_vs_SNR_BO_%ddB',bo),'Visible','on','NumberTitle','off'); 
    hold on; grid on; set(gcf,'Color','w');

    all_evm = [];

    for m = 1:numel(methods)
        ber = zeros(numel(SNRs),1); 
        evm = zeros(numel(SNRs),1);
        evm_sd = zeros(numel(SNRs),1);
        for is = 1:numel(SNRs)
            rec = R.(keyname(methods{m}, SNRs(is), bo));
            ber(is)   = rec.BER; 
            evm(is)   = rec.EVM_pct;
            evm_sd(is)= rec.EVM_std;
        end
        xj = SNRs + (m-2.5)*0.15;  % tiny jitter to prevent overlap
        figure(f1); semilogy(xj, max(ber,1e-7), mk{m}, 'DisplayName',methods{m}, 'LineWidth',1.4, 'MarkerSize',6);

        evm_plot = evm; evm_plot(~isfinite(evm_plot)) = NaN;
        evm_err  = evm_sd; evm_err(~isfinite(evm_err)) = 0;
        all_evm  = [all_evm; evm_plot(:)]; %#ok<AGROW>

        figure(f2); errorbar(xj, evm_plot, evm_err, mk{m}, 'DisplayName',methods{m}, ...
            'LineWidth',1.2, 'MarkerSize',6, 'CapSize',3);
    end

    figure(f1); xlabel('SNR [dB]'); ylabel('BER'); ylim([1e-5 1]); 
    title(sprintf('BER vs SNR (BO=%d dB)',bo)); legend('show','Location','southwest');

    all_evm = all_evm(isfinite(all_evm));
    if isempty(all_evm), yhi = 10; else, yhi = max(5, min(120, 1.10*max(all_evm))); end
    figure(f2); xlabel('SNR [dB]'); ylabel('EVM [%]'); ylim([0 yhi]);
    title(sprintf('EVM vs SNR (BO=%d dB)',bo)); legend('show','Location','northeast');

    mcs_save_fig(f1, fullfile(outdir,sprintf('ber_vs_snr_BO_%ddB.png',bo)));
    mcs_save_fig(f2, fullfile(outdir,sprintf('evm_vs_snr_BO_%ddB.png',bo)));
end

% ---- ACLR vs back-off ----
f_aclr = figure('Name','ACLR_vs_Backoff'); hold on; grid on; set(gcf,'Color','w');
for m = 1:numel(methods)
    aclr = zeros(numel(BOs),1);
    for ib = 1:numel(BOs)
        rec = R.(keyname(methods{m}, SNRs(1), BOs(ib)));
        aclr(ib) = rec.ACLR_post_dB;
    end
    plot(BOs, aclr,'-o','DisplayName',methods{m},'LineWidth',1.4,'MarkerSize',6);
end
xlabel('PA Output Back-Off [dB]'); ylabel('ACLR (post-PA) [dB]');
title('ACLR vs PA Output Back-Off'); legend('show','Location','southeast');
mcs_save_fig(f_aclr, fullfile(outdir,'aclr_vs_backoff.png'));

% ---------- Efficiency lens: EVM improvement vs BASE at BO = 0 dB ----------
target_BO    = 0;           % [dB] fixed back-off
target_SNRs  = [20 25];     % [dB]
methods_eff  = {'SLM','PTS','TR'};  % compare against BASE

nMeth = numel(methods_eff);
nTgt  = numel(target_SNRs);

evmGain = nan(nMeth, nTgt); % percentage-point EVM improvement vs BASE

for ti = 1:nTgt
    snr_t = target_SNRs(ti);

    % BASE reference EVM at this SNR and BO
    rec_base = R.(keyname('BASE', snr_t, target_BO));
    evm_base = rec_base.EVM_pct;

    for m = 1:nMeth
        rec_m = R.(keyname(methods_eff{m}, snr_t, target_BO));
        evm_m = rec_m.EVM_pct;

        % positive = better (lower EVM than BASE)
        evmGain(m, ti) = evm_base - evm_m;
    end
end

xCats = categorical(methods_eff);
xCats = reordercats(xCats, methods_eff);

f_tgt = figure('Name','EVM_Improvement_vs_BASE_BO0'); set(gcf,'Color','w');
hold on; grid on;

bh = bar(xCats, evmGain);  % grouped bars: one group per method

ylabel('EVM improvement vs BASE [percentage points]');
xlabel('Method');
legend(arrayfun(@(s) sprintf('SNR = %d dB', s), target_SNRs, ...
       'UniformOutput', false), 'Location','northwest');
title(sprintf('EVM Gain at BO = %d dB', target_BO), 'Interpreter','none');

gmax = max(evmGain(:), [], 'omitnan');
if ~isfinite(gmax) || gmax <= 0
    gmax = 1;   % safe default
end
ylim([0 gmax + 0.2]);

% optional: value labels above bars
for ti = 1:nTgt
    xEnd = bh(ti).XEndPoints;
    yVal = bh(ti).YData;
    for i = 1:numel(yVal)
        if ~isnan(yVal(i))
            text(xEnd(i), yVal(i) + 0.03, sprintf('%.2f', yVal(i)), ...
                'HorizontalAlignment','center','FontSize',8);
        end
    end
end

mcs_save_fig(f_tgt, fullfile(outdir,'evm_gain_vs_BASE_BO0.png'));


% ---------- Sensitivity plot (SLM U=2 vs 4, PTS B=2 vs 4) ----------
sens_SNR = 20; sens_BO = 0; % BO=0 dB, SNR=20 dB
f_sens = figure('Name','Sensitivity_SLM_PTS'); set(gcf,'Color','w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% SLM U sweep
nexttile; hold on; grid on;
Us = [2 4];
evmU = zeros(size(Us)); evmU_sd = zeros(size(Us));
for k = 1:numel(Us)
    SLM2 = cfg.SLM; SLM2.U = Us(k);
    [frames_SLM2, ~, ~] = slm_reduce(framesFD, SLM2);
    x_SLM2 = ofdm_mod(frames_SLM2, cfg, meta);
    cfg_i = cfg; cfg_i.sim.SNRdB = sens_SNR; cfg_i.PA.backoff_dB = sens_BO;
    y = pa_nonlinearity(x_SLM2, cfg_i.PA);

    evm_runs = zeros(nAvg,1);
    for it = 1:nAvg
        rng(900 + 10*k + it);
        [evm_runs(it), ~] = link_metrics(y, frames_SLM2, dataIdx, cfg_i, meta, txSymIdxPerSym);
    end
    evmU(k)   = mean(evm_runs);
    evmU_sd(k)= std(evm_runs,0,1);
end
errorbar(categorical("U="+string(Us)), evmU, evmU_sd, '-o', 'LineWidth',1.4, 'CapSize',4);
ylabel('EVM [%]'); title(sprintf('SLM Sensitivity @ BO=%d dB, SNR=%d dB', sens_BO, sens_SNR));

% PTS B sweep
nexttile; hold on; grid on;
Bs = [2 4];
evmB = zeros(size(Bs)); evmB_sd = zeros(size(Bs));
for k = 1:numel(Bs)
    PTS2 = cfg.PTS; PTS2.B = Bs(k);
    [frames_PTS2, ~] = pts_reduce(framesFD, PTS2);
    x_PTS2 = ofdm_mod(frames_PTS2, cfg, meta);
    cfg_i = cfg; cfg_i.sim.SNRdB = sens_SNR; cfg_i.PA.backoff_dB = sens_BO;
    y = pa_nonlinearity(x_PTS2, cfg_i.PA);

    evm_runs = zeros(nAvg,1);
    for it = 1:nAvg
        rng(1200 + 10*k + it);
        [evm_runs(it), ~] = link_metrics(y, frames_PTS2, dataIdx, cfg_i, meta, txSymIdxPerSym);
    end
    evmB(k)   = mean(evm_runs);
    evmB_sd(k)= std(evm_runs,0,1);
end
errorbar(categorical("B="+string(Bs)), evmB, evmB_sd, '-o', 'LineWidth',1.4, 'CapSize',4);
ylabel('EVM [%]'); title(sprintf('PTS Sensitivity @ BO=%d dB, SNR=%d dB', sens_BO, sens_SNR));

mcs_save_fig(f_sens, fullfile(outdir,'sensitivity_slm_pts.png'));

disp('Done. See results/ for figures, CSV, MAT, and added target/sensitivity plots.');
end

%% ========================= Local functions =========================
function cfg = config()
    % OFDM numerology
    cfg.Nsub     = 256;
    cfg.CP       = 16;
    cfg.oversamp = 4;
    cfg.Nsym     = 700;
    cfg.M        = 16;

    % Tone reservation fraction
    cfg.TR.reservePct = 0.10;

    % Sampling / occupied BW (only used for ACLR)
    cfg.fs = 15.36e6;
    cfg.BW = 3.84e6;

    % SNR sweep
    cfg.sim.SNRdB_list = 10:5:30;

    % PA (stronger nonlinearity + wider BO sweep)
    cfg.PA.Asat = 0.8;
    cfg.PA.p    = 5;
    cfg.PA.backoff_dB      = 6;
    cfg.PA.backoff_list_dB = [0 2 4 6 8];

    % SLM
    cfg.SLM.U        = 4;
    cfg.SLM.phaseSet = exp(1j*pi*[0 1]);   % {+1,-1}

    % PTS
    cfg.PTS.B        = 4;
    cfg.PTS.phaseSet = exp(1j*pi*[0 1]);   % {0,pi}
    cfg.PTS.search   = 'greedy';

    % Receiver equalization options (avoid over-correcting PA distortion)
    cfg.rx.use_equalizer      = false;   % set true for subset LS
    cfg.rx.eq_subset_frac     = 0.10;

    % Optional mild multipath (adds a bit of frequency selectivity)
    cfg.channel.use_multipath = false;
    cfg.channel.h             = [1 0 0.35*exp(1j*pi/4)];
end

function meta = build_meta(cfg, Ns)
    meta.Ns = Ns; meta.L = cfg.oversamp; meta.Nfft = cfg.Nsub*meta.L;
    meta.cp = cfg.CP*meta.L; meta.symLen = meta.Nfft + meta.cp;
end

function [framesFD, dataIdx, trIdx, bitsTx, txSymIdxPerSym] = gen_ofdm_frames(cfg) %#ok<INUSD>
    N = cfg.Nsub; M = cfg.M; bps = log2(M);
    % Reserve TR tones (uniformly inside band, avoid DC/end tones)
    trIdx = [];
    if cfg.TR.reservePct > 0
        k = max(1, round(N*cfg.TR.reservePct));
        trIdx = unique(round(linspace(2, N-1, k))).';
    end
    mask = true(N,1); mask(trIdx)=false; dataIdx = find(mask);
    numData = numel(dataIdx);

    % Bits for all symbols (column-major)
    bitsTx = randi([0 1], numData*bps*cfg.Nsym, 1);
    bitsMat = reshape(bitsTx, numData*bps, cfg.Nsym);

    % Frames and ground-truth symbol indices
    framesFD = zeros(N, cfg.Nsym);
    txSymIdxPerSym = zeros(numData, cfg.Nsym);
    for n = 1:cfg.Nsym
        bcol = bitsMat(:,n);
        s    = qammod(bcol, M, 'InputType','bit', 'UnitAveragePower', true); % numData x 1
        X    = zeros(N,1); X(dataIdx)=s; framesFD(:,n)=X;
        txSymIdxPerSym(:,n) = qamdemod(s, M, 'UnitAveragePower', true);
    end
end

function x = ofdm_mod(framesFD, cfg, meta)
    [N, Ns] = size(framesFD); Nfft = meta.Nfft;
    Xos = zeros(Nfft, Ns); half = floor(N/2);
    for n = 1:Ns
        X = framesFD(:,n);
        Xos([1:half, end-(N-half)+1:end], n) = [X(1:half); X(half+1:end)];
    end
    x_noCP = ifft(ifftshift(Xos,1), Nfft, 1) * sqrt(Nfft);
    x = [x_noCP(end-meta.cp+1:end,:); x_noCP];
    x = x(:);
end

function framesFD = ofdm_demod(y, cfg, meta)
    Nfft = meta.Nfft; cp=meta.cp; Ns=meta.Ns;
    Lsym  = Nfft + cp; Nneed = Lsym*Ns;
    if numel(y) < Nneed, y = [y(:); zeros(Nneed - numel(y),1)]; end
    y = y(1:Nneed);
    y = reshape(y, Lsym, Ns);
    Y = fftshift(fft(y(cp+1:end,:), Nfft, 1)/sqrt(Nfft),1);
    N = cfg.Nsub; half=floor(N/2);
    framesFD = zeros(N, Ns);
    for n=1:Ns
        Yo = Y(:,n);
        framesFD(:,n) = [Yo(1:half); Yo(end-(N-half)+1:end)];
    end
end

% ---------- PAPR helpers ----------
function p = papr_db_symbol(x, meta, prc)
    L  = meta.symLen;  Ns = floor(numel(x)/L);
    x  = x(1:Ns*L);    X  = reshape(x, L, Ns);
    Pav  = mean(abs(X).^2, 1);
    Ppk  = max(abs(X).^2, [], 1);
    papr = 10*log10(Ppk ./ max(Pav, eps));
    p    = prctile_local(papr, prc);
end

function ccdf_plot(x, labelstr, meta)
    symLen = meta.symLen;  Ns = floor(numel(x)/symLen);
    x = x(1:Ns*symLen);  x_mat = reshape(x, symLen, Ns);
    Pavg = mean(abs(x_mat).^2, 1);
    Ppk  = max(abs(x_mat).^2, [], 1);
    papr_db = 10*log10(Ppk ./ max(Pavg,eps));
    papr_db = sort(papr_db(:));
    N = numel(papr_db); thr = papr_db; 
    ccdf = 1 - ((1:N)' - 0.5)/N;
    semilogy(thr, max(ccdf,1e-6), 'DisplayName', labelstr, 'LineWidth',1.2); hold on;
    ylim([1e-6 1]); xlim([0 15]); grid on;
end

% ---------- Methods ----------
function [framesFD_opt, si, phases] = slm_reduce(framesFD, SLM)
    [N, Ns] = size(framesFD);
    U = SLM.U; K = numel(SLM.phaseSet);
    phases = cell(U,1);
    for u = 1:U
        idx = randi(K, N, 1); phu = SLM.phaseSet(idx);
        phases{u} = reshape(phu, [N,1]);
    end
    framesFD_opt = zeros(N, Ns); si = zeros(Ns,1);
    for n = 1:Ns
        X = framesFD(:,n);
        bestP = Inf; bestXu = X; bestu = 1;
        for u = 1:U
            Xu = X .* phases{u};
            x_tmp = time_proxy_from_fd(Xu);
            PAPR = 10*log10( max(abs(x_tmp).^2) / max(mean(abs(x_tmp).^2),eps) );
            if PAPR < bestP, bestP = PAPR; bestXu = Xu; bestu = u; end
        end
        framesFD_opt(:,n) = bestXu; si(n) = bestu;
    end
end

function [framesFD_opt, state] = pts_reduce(framesFD, PTS)
    [N, Ns] = size(framesFD);
    B = PTS.B; phaseSet = PTS.phaseSet;
    blkIdx = make_blocks(N, B);
    framesFD_opt = zeros(size(framesFD));
    state.blocks = blkIdx; state.phases = zeros(B, Ns);
    for n = 1:Ns
        X = framesFD(:,n); ph_sel = ones(B,1); bestX = X;
        for b = 1:B
            bestP = Inf; bestp = 1;
            for pi = 1:numel(phaseSet)
                ph_try = ph_sel; ph_try(b) = phaseSet(pi);
                Xb = apply_block_phases(X, blkIdx, ph_try);
                x_tmp = time_proxy_from_fd(Xb);
                PAPR = 10*log10( max(abs(x_tmp).^2) / max(mean(abs(x_tmp).^2),eps) );
                if PAPR < bestP, bestP = PAPR; bestp = pi; bestX = Xb; end
            end
            ph_sel(b) = phaseSet(bestp);
        end
        framesFD_opt(:,n) = bestX; state.phases(:,n) = ph_sel;
    end
end

function [framesFD_opt, c_all] = tone_reservation(framesFD, trIdx, TR, cfg)
    [~, Ns] = size(framesFD);
    if isempty(trIdx), framesFD_opt = framesFD; c_all = []; return; end
    K = numel(trIdx); lambda = 1e-3; framesFD_opt = framesFD; c_all = zeros(K,Ns);
    L = cfg.oversamp; Nfft = cfg.Nsub*L; half = floor(cfg.Nsub/2);
    T = zeros(Nfft, K);
    for k = 1:K
        unit = zeros(cfg.Nsub,1); unit(trIdx(k)) = 1;
        Xos = zeros(Nfft,1);
        Xos([1:half, end-(cfg.Nsub-half)+1:end]) = [unit(1:half); unit(half+1:end)];
        T(:,k) = ifft(ifftshift(Xos,1), Nfft, 1) * sqrt(Nfft);
    end
    for n = 1:Ns
        Xd = framesFD(:,n); 
        Xos_d = zeros(Nfft,1);
        Xos_d([1:half, end-(cfg.Nsub-half)+1:end]) = [Xd(1:half); Xd(half+1:end)];
        x_d = ifft(ifftshift(Xos_d,1), Nfft, 1) * sqrt(Nfft);
        c = -(T'*T + lambda*eye(K)) \ (T'*x_d);
        Xr = zeros(cfg.Nsub,1); Xr(trIdx) = c;
        framesFD_opt(:,n) = Xd + Xr; c_all(:,n) = c;
    end
end

% ---------- Metrics ----------
function y = pa_nonlinearity(x, PA)
    Ain = 10^(-PA.backoff_dB/20); xsc = Ain * x(:);
    A = PA.Asat; p = PA.p;
    mag = abs(xsc); ph = angle(xsc);
    gain = (1 + (mag./A).^(2*p)).^(-1/(2*p));
    y = gain .* mag .* exp(1j*ph);
end

function [Px, aclr_dB] = measure_aclr(x, cfg)
    % PSD via Welch
    nfft = 4096;
    [Px,F] = pwelch(x, hann(1024,'periodic'), 512, nfft, cfg.fs, 'centered');
    BW   = cfg.BW; fc = 0;      % baseband

    inband  = (abs(F-fc) <= BW/2);
    adjU    = (F >=  BW/2) & (F <=  3*BW/2);
    adjL    = (F <= -BW/2) & (F >= -3*BW/2);

    Pmain = trapz(F(inband), Px(inband));
    Padj  = 0.5*( trapz(F(adjU), Px(adjU)) + trapz(F(adjL), Px(adjL)) );

    aclr_dB = 10*log10(Pmain/max(Padj,eps));   % positive, larger is better
end

function [evm_percent, ber] = link_metrics(y_pa, framesFD_ref, dataIdx, cfg, meta, txSymIdxPerSym)
    % Optional channel coloring (mild multipath)
    y = y_pa;
    if isfield(cfg,'channel') && cfg.channel.use_multipath
        y = conv(y, cfg.channel.h(:), 'same');
    end
    % AWGN
    y = awgn(y, cfg.sim.SNRdB, 'measured');

    % OFDM demod
    framesFD_rx = ofdm_demod(y, cfg, meta);
    Srx  = framesFD_rx(dataIdx,:);       % data tones (Ndata x Nsym)
    Sref = framesFD_ref(dataIdx,:);      % reference constellation

    % Equalization policy
    if cfg.rx.use_equalizer
        frac = max(0,min(1,cfg.rx.eq_subset_frac));
        Nt = numel(Sref(:,1));
        K  = max(4, round(frac*Nt));
        alpha = zeros(1,size(Sref,2));
        for n = 1:size(Sref,2)
            idx = randperm(Nt, K);
            num = sum(conj(Sref(idx,n)).*Srx(idx,n));
            den = sum(abs(Sref(idx,n)).^2) + eps;
            alpha(n) = num/den;
        end
        Srx_eq = Srx ./ alpha;
    else
        Srx_eq = Srx;
    end

    % EVM
    err  = Srx_eq - Sref;
    evm_rms = sqrt( mean(abs(err(:)).^2) / max(mean(abs(Sref(:)).^2), eps) );
    evm_percent = 100 * evm_rms;

    % BER via symbol indices
    M   = cfg.M; bps = log2(M);
    rxSymIdx = qamdemod(Srx_eq(:), M, 'UnitAveragePower', true);
    txSymIdx = txSymIdxPerSym(:);
    rxBits = de2bi(rxSymIdx, bps, 'left-msb');
    txBits = de2bi(txSymIdx, bps, 'left-msb');
    ber = mean(rxBits(:) ~= txBits(:));
end

% ---------- Results & utilities ----------
function key = keyname(method, SNRdB, BOdB)
    key = sprintf('%s_SNR%02d_BO%02d', method, round(SNRdB), round(BOdB));
end

function Tab = results_to_table(R)
    K = fieldnames(R);
    Method = strings(numel(K),1); SNRdB = zeros(numel(K),1); BOdB = zeros(numel(K),1);
    PAPR_pre_dB = zeros(numel(K),1); ACLR_pre_dB = zeros(numel(K),1);
    ACLR_post_dB = zeros(numel(K),1); EVMpct = zeros(numel(K),1); EVMstd = zeros(numel(K),1); BER = zeros(numel(K),1); BERstd = zeros(numel(K),1);
    for i = 1:numel(K)
        rec = R.(K{i});
        Method(i) = string(rec.method);
        SNRdB(i)   = rec.SNRdB;
        BOdB(i)    = rec.BOdB;
        PAPR_pre_dB(i)  = rec.PAPR_pre_dB;
        ACLR_pre_dB(i)  = rec.ACLR_pre_dB;
        ACLR_post_dB(i) = rec.ACLR_post_dB;
        EVMpct(i)       = rec.EVM_pct;
        EVMstd(i)       = rec.EVM_std;
        BER(i)          = rec.BER;
        BERstd(i)       = rec.BER_std;
    end
    Tab = table(Method,SNRdB,BOdB,PAPR_pre_dB,ACLR_pre_dB,ACLR_post_dB,EVMpct,EVMstd,BER,BERstd);
end

function mcs_save_results(Tab, cfg, meta, R, outdir)
    ts = datestr(now,'yyyymmdd_HHMMSS');
    csv_main = fullfile(outdir,'summary_all.csv');
    mat_main = fullfile(outdir,'results_all.mat');
    try, writetable(Tab, csv_main);
    catch, warning('CSV locked. Writing timestamped copy instead.'); ...
          writetable(Tab, fullfile(outdir, ['summary_all_' ts '.csv']));
    end
    try, save(mat_main,'cfg','meta','R','Tab');
    catch, warning('MAT locked. Writing timestamped copy instead.'); ...
          save(fullfile(outdir, ['results_all_' ts '.mat']),'cfg','meta','R','Tab');
    end
end

function mcs_save_fig(h, fname)
    if ~ishandle(h), return; end
    try, set(h,'ToolBar','none','MenuBar','none'); catch, end
    ax = findall(h,'Type','axes');
    for k = 1:numel(ax)
        try, axtoolbar(ax(k),'Visible','off'); catch, end
    end
    drawnow;
    print(h, fname, '-dpng', '-r200');
end

% ---------- helpers ----------
function v = prctile_local(x, p)
    x = sort(x(:)); n = numel(x);
    if n == 0, v = NaN; return; end
    t = (p/100)*(n-1) + 1; lo = floor(t); hi = ceil(t);
    if lo == hi, v = x(lo); else, a=x(lo); b=x(hi); v = a + (t-lo)*(b-a); end
end

function s = lbl(x)
    if isnan(x), s = '—'; else, s = sprintf('%.0f', x); end
end

function x_proxy = time_proxy_from_fd(X)
    X = X(:); N = length(X); half = floor(N/2);
    Xn = zeros(N,1);
    Xn([1:half, end-(N-half)+1:end]) = [X(1:half); X(half+1:end)];
    x_proxy = ifft(ifftshift(Xn,1), N, 1) * sqrt(N);
end

function blkIdx = make_blocks(N, B)
    edges = round(linspace(0,N,B+1)); blkIdx = cell(B,1);
    for b=1:B, blkIdx{b} = (edges(b)+1):edges(b+1); end
end

function Xout = apply_block_phases(X, blkIdx, phSel)
    Xout = X;
    for b=1:numel(blkIdx)
        Xout(blkIdx{b}) = X(blkIdx{b}).*phSel(b);
    end
end
