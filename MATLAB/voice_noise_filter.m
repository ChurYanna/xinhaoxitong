%% voice_noise_filter_2024a_fix2.m
% 任务（第七点）：噪声抑制与滤波器设计（R2024a，严肃修复版）
% 更改点：
% - 表格写出时先构造列变量，再用 VariableNames 显式创建，避免 table 参数解析冲突
% - 修正分类打印变量名（isCover050）
% - 使用 Butterworth IIR 带通（300–3400 Hz）+ 可选 50 Hz IIR 陷波，filtfilt 零相位
% - 提供 FIR 线性相位备选（fir1），将 useIIR=false 可切换

clc; clear; close all;

%% 1) 路径与文件
baseDir = "C:\Users\Administrator\Desktop\学业备份\大三上\信号与系统\《信号与系统》期末大作业\WAV\";
file010 = baseDir + "voice_010.wav"; % 原唱-中文
file020 = baseDir + "voice_020.wav"; % 原唱-粤语
file030 = baseDir + "voice_030.wav"; % 安静翻唱
file040 = baseDir + "voice_040.wav"; % 嘈杂翻唱（待增强）
file050 = baseDir + "voice_050.wav"; % 增强输出

outDir = fullfile(pwd, "feature_outputs_task7");
if ~exist(outDir, "dir"), mkdir(outDir); end

%% 2) 环境检测
hasAudioTB      = license('test','Audio_Toolbox') && exist('pitch','file')==2;
hasDetectSpeech = hasAudioTB && exist('detectSpeech','file')==2;
hasSTFT         = exist('stft','file')==2;

%% 3) 读取并预处理
[x040, Fs040] = audioread(file040);
if size(x040,2)>1, x040 = mean(x040,2); end
x040 = x040 - mean(x040);
x040 = x040 / (max(abs(x040))+eps);

[x030, Fs030] = audioread(file030);
if size(x030,2)>1, x030 = mean(x030,2); end
x030 = x030 - mean(x030);
x030 = x030 / (max(abs(x030))+eps);

%% 4) 频谱分析（040，滤波前）
% STFT 参数
stftWinLen = max(256, 2^nextpow2(round(25e-3*Fs040)));
stftOverlapLen = round(0.75 * stftWinLen);
stftNFFT = max(stftWinLen, 2^nextpow2(stftWinLen));
wSTFT = hamming(stftWinLen,"periodic");

% STFT 可视化
if hasSTFT
    [S040, F040s, T040s] = stft(x040, Fs040, Window=wSTFT, OverlapLength=stftOverlapLen, FFTLength=stftNFFT);
else
    [S040, F040s, T040s] = spectrogram(x040, wSTFT, stftOverlapLen, stftNFFT, Fs040, "yaxis");
end
S040_dB = 20*log10(abs(S040)+eps);
figS040 = figure("Color","w","Position",[60,60,1200,500]);
imagesc(T040s, F040s, S040_dB); axis xy; colormap jet; colorbar;
xlabel("Time (s)"); ylabel("Frequency (Hz)");
title("Spectrogram (pre-filter): voice_040");
ylim([0, min(8000, Fs040/2)]);
set(gca,"CLim",[prctile(S040_dB(:),5), prctile(S040_dB(:),95)]);
exportgraphics(figS040, fullfile(outDir,"spectrogram_voice_040_pre.png"), "Resolution",160);
close(figS040);

% Welch PSD
welchWinLen = max(256, 2^nextpow2(round(0.046*Fs040)));
wWelch = hamming(welchWinLen,"periodic");
welchNFFT = max(4096, 2^nextpow2(welchWinLen));
welchOverlap = 0.5;
noWelch = round(welchOverlap*welchWinLen);
[Pxx040, f040] = pwelch(x040, wWelch, noWelch, welchNFFT, Fs040, "onesided");
Pxx040_dB = 10*log10(Pxx040+eps);
figP040 = figure("Color","w","Position",[80,80,950,420]);
plot(f040, Pxx040_dB, "LineWidth",1.2); grid on;
xlabel("Frequency (Hz)"); ylabel("PSD (dB/Hz)");
title("Welch PSD (pre-filter): voice_040");
xlim([0, min(8000, Fs040/2)]);
exportgraphics(figP040, fullfile(outDir,"psd_voice_040_pre.png"), "Resolution",160);
close(figP040);

% 频带能量（供报告）
bandsExplain = [0 300; 300 1000; 1000 3400; 3400 8000];
Eband = zeros(size(bandsExplain,1),1);
for i=1:size(bandsExplain,1)
    Eband(i) = bandpower(Pxx040, f040, bandsExplain(i,:), "psd");
end
fprintf("Pre-filter energy by bands (voice_040):\n");
fprintf("  0–300 Hz:     %.6f\n", Eband(1));
fprintf("  300–1000 Hz:  %.6f\n", Eband(2));
fprintf("  1–3.4 kHz:    %.6f\n", Eband(3));
fprintf("  3.4–8 kHz:    %.6f\n\n", Eband(4));

%% 5) 滤波器设计与应用
useNotch = true;  % 是否使用 50 Hz 陷波
useIIR   = true;  % true=Butterworth IIR；false=FIR(fir1)

xStep = x040;
if useNotch
    % 50 Hz IIR 陷波（Q≈35）
    wo = 50/(Fs040/2); bw = wo/35;
    [bN, aN] = iirnotch(wo, bw);
    xStep = filtfilt(bN, aN, xStep);
end

if useIIR
    % Butterworth IIR 带通，6阶，300–3400 Hz
    [bBP, aBP] = butter(6, [300 3400]/(Fs040/2), 'bandpass');
    x050 = filtfilt(bBP, aBP, xStep);
else
    % FIR 线性相位（fir1）
    firOrder = 1024;
    bFIR = fir1(firOrder, [300 3400]/(Fs040/2), 'bandpass', hamming(firOrder+1));
    x050 = filtfilt(bFIR, 1, xStep);
end
x050 = x050 / (max(abs(x050))+eps);
audiowrite(file050, x050, Fs040);

%% 6) 滤波后可视化
% STFT
if hasSTFT
    [S050, F050s, T050s] = stft(x050, Fs040, Window=wSTFT, OverlapLength=stftOverlapLen, FFTLength=stftNFFT);
else
    [S050, F050s, T050s] = spectrogram(x050, wSTFT, stftOverlapLen, stftNFFT, Fs040, "yaxis");
end
S050_dB = 20*log10(abs(S050)+eps);
figS050 = figure("Color","w","Position",[60,60,1200,500]);
imagesc(T050s, F050s, S050_dB); axis xy; colormap jet; colorbar;
xlabel("Time (s)"); ylabel("Frequency (Hz)");
title("Spectrogram (post-filter): voice_050");
ylim([0, min(8000, Fs040/2)]);
set(gca,"CLim",[prctile(S050_dB(:),5), prctile(S050_dB(:),95)]);
exportgraphics(figS050, fullfile(outDir,"spectrogram_voice_050_post.png"), "Resolution",160);
close(figS050);

% Welch PSD（post）
[Pxx050, f050] = pwelch(x050, wWelch, noWelch, welchNFFT, Fs040, "onesided");
Pxx050_dB = 10*log10(Pxx050+eps);
figP050 = figure("Color","w","Position",[80,80,950,420]);
plot(f050, Pxx050_dB,"LineWidth",1.2); grid on;
xlabel("Frequency (Hz)"); ylabel("PSD (dB/Hz)");
title("Welch PSD (post-filter): voice_050");
xlim([0, min(8000, Fs040/2)]);
exportgraphics(figP050, fullfile(outDir,"psd_voice_050_post.png"), "Resolution",160);
close(figP050);

% 叠加对比（pre vs post）
figOverlay = figure("Color","w","Position",[100,100,980,420]);
plot(f040, Pxx040_dB, "r-", "LineWidth",1.2); hold on;
plot(f050, Pxx050_dB, "b-", "LineWidth",1.2);
grid on; xlabel("Frequency (Hz)"); ylabel("PSD (dB/Hz)");
legend("voice_040 (pre)","voice_050 (post)","Location","best");
title("PSD Overlay: pre vs post (voice_040 -> voice_050)");
xlim([0, min(8000, Fs040/2)]);
exportgraphics(figOverlay, fullfile(outDir,"psd_overlay_040_pre_050_post.png"), "Resolution",160);
close(figOverlay);

%% 7) SNR 估计（前后）
winLenSNR = max(256, round(0.03*Fs040));
hopSNR = round(0.5*winLenSNR);
if hasDetectSpeech
    seg040 = detectSpeech(x040, Fs040);
    seg050 = detectSpeech(x050, Fs040);
else
    seg040 = vad_energy_segments(x040, Fs040, winLenSNR, hopSNR);
    seg050 = vad_energy_segments(x050, Fs040, winLenSNR, hopSNR);
end
SNR040 = estimate_snr_by_segments(x040, seg040, Fs040);
SNR050 = estimate_snr_by_segments(x050, seg050, Fs040);
fprintf("Estimated SNR (fullband): 040(pre)=%.2f dB, 050(post)=%.2f dB, improvement=%.2f dB\n\n", ...
    SNR040, SNR050, SNR050 - SNR040);

%% 8) 特征 Z 与分类（验证 050 是否归入“翻唱”）
feat010 = compute_features_for_files({file010}); feat010 = feat010(1);
feat020 = compute_features_for_files({file020}); feat020 = feat020(1);
feat030 = compute_features_for_files({file030}); feat030 = feat030(1);

thr_HF   = mid_mean([feat010.ratio_high, feat020.ratio_high], feat030.ratio_high);
thr_roll = mid_mean([feat010.rolloff95_Hz, feat020.rolloff95_Hz], feat030.rolloff95_Hz);

feat040 = compute_features_for_files({file040}); feat040 = feat040(1);
feat050 = compute_features_for_files({file050}); feat050 = feat050(1);

isCover050 = (feat050.ratio_high >= thr_HF) || (feat050.rolloff95_Hz >= thr_roll);

% 表格写出（修复：先构造列，再指定 VariableNames）
fileNames       = ["voice_040"; "voice_050"];
ratio_low       = [feat040.ratio_low;       feat050.ratio_low];
ratio_mid       = [feat040.ratio_mid;       feat050.ratio_mid];
ratio_high      = [feat040.ratio_high;      feat050.ratio_high];
rolloff95_Hz    = [feat040.rolloff95_Hz;    feat050.rolloff95_Hz];
SFM_0_2k        = [feat040.SFM_0_2k;        feat050.SFM_0_2k];
HNR_dB          = [feat040.HNR_dB;          feat050.HNR_dB];
SNR_est_dB      = [SNR040;                  SNR050];

varNames = ["file","ratio_low","ratio_mid","ratio_high","rolloff95_Hz","SFM_0_2k","HNR_dB","SNR_est_dB"];
T = table(fileNames, ratio_low, ratio_mid, ratio_high, rolloff95_Hz, SFM_0_2k, HNR_dB, SNR_est_dB, ...
          'VariableNames', varNames);
outCsv = fullfile(outDir, "results_task7_voice040_050_fix2.csv");
writetable(T, outCsv, "WriteMode","overwrite");
fprintf("Exported: %s\n", outCsv);

% 分类结论打印（修复变量名）
fprintf("Classification (Z-based) for voice_050: %s\n", ternary(isCover050,"cover","original"));
fprintf("  Thresholds: HF_ratio_thr=%.4f, rolloff95_thr=%.1f Hz\n", thr_HF, thr_roll);
fprintf("  voice_050 features: HF_ratio=%.4f, rolloff95=%.1f Hz, HNR=%.2f dB, SFM_0_2k=%.3f, SNR=%.2f dB\n\n", ...
    feat050.ratio_high, feat050.rolloff95_Hz, feat050.HNR_dB, feat050.SFM_0_2k, SNR050);

%% 9) 与 voice_030 的 PSD 叠加对比
welchWinLen30 = max(256, 2^nextpow2(round(0.046*Fs030)));
wWelch30 = hamming(welchWinLen30,"periodic");
noWelch30 = round(0.5*welchWinLen30);
welchNFFT30 = max(4096, 2^nextpow2(welchWinLen30));
[Pxx030, f030] = pwelch(x030, wWelch30, noWelch30, welchNFFT30, Fs030, "onesided");
Pxx030_dB = 10*log10(Pxx030+eps);

figComp = figure("Color","w","Position",[120,120,1000,430]);
plot(f030, Pxx030_dB, "k-", "LineWidth",1.2); hold on;
plot(f050, Pxx050_dB, "b-", "LineWidth",1.2);
grid on; xlabel("Frequency (Hz)"); ylabel("PSD (dB/Hz)");
legend("voice_030 (quiet cover)","voice_050 (filtered noisy cover)","Location","best");
title("PSD Comparison: 030 vs 050");
xlim([0, min([8000, Fs030/2, Fs040/2])]);
exportgraphics(figComp, fullfile(outDir,"psd_comparison_030_vs_050_fix2.png"), "Resolution",160);
close(figComp);

disp("Task7 (fix2) done. 请将生成的 PNG 与 CSV 发我，我会撰写报告的效果与局限性。");

%% ======= 本地函数 =======
function seg = vad_energy_segments(x, Fs, winLen, hop)
    n = length(x);
    idx = 1:hop:(n-winLen+1);
    E = zeros(numel(idx),1);
    for i=1:numel(idx)
        segWin = x(idx(i):idx(i)+winLen-1) .* hamming(winLen,"periodic");
        E(i) = sum(segWin.^2);
    end
    thr = median(E) + 0.5*(prctile(E,75)-median(E));
    mask = E >= thr;
    seg = [];
    if any(mask)
        starts = idx([true; diff(mask)>0] & mask);
        ends   = idx([diff(mask)<0; true] & mask) + winLen - 1;
        seg = [starts(:), ends(:)];
    end
end

function SNRdB = estimate_snr_by_segments(x, seg, Fs)
    if isempty(seg), SNRdB = NaN; return; end
    n = numel(x);
    voiceMask = false(n,1);
    for i=1:size(seg,1)
        voiceMask(seg(i,1):min(n,seg(i,2))) = true;
    end
    E_voice = sum(x(voiceMask).^2) / max(1, nnz(voiceMask));
    E_noise = sum(x(~voiceMask).^2) / max(1, nnz(~voiceMask));
    SNRdB = 10*log10((E_voice+eps)/(E_noise+eps));
end

function feats = compute_features_for_files(fileList)
    feats = struct([]);
    for k = 1:numel(fileList)
        [x, Fs] = audioread(fileList{k});
        if size(x,2)>1, x = mean(x,2); end
        x = x - mean(x);
        x = x / (max(abs(x))+eps);

        % Welch PSD
        winLen = max(256, 2^nextpow2(round(0.046*Fs)));
        wWelch = hamming(winLen,"periodic");
        noWelch = round(0.5*winLen);
        nfft = max(4096, 2^nextpow2(winLen));
        [Pxx, f] = pwelch(x, wWelch, noWelch, nfft, Fs, "onesided");

        % 频带能量比例
        E_low  = bandpower(Pxx, f, [0,    min(1000,Fs/2)], "psd");
        E_mid  = bandpower(Pxx, f, [1000, min(4000,Fs/2)], "psd");
        E_high = bandpower(Pxx, f, [4000, min(8000,Fs/2)], "psd");
        E_tot  = E_low + E_mid + E_high + eps;
        ratio_low  = E_low  / E_tot;
        ratio_mid  = E_mid  / E_tot;
        ratio_high = E_high / E_tot;

        % 滚降95%
        rolloff95 = spectral_rolloff_psd(Pxx, f, 0.95);
        % 谱平坦度 0–2k
        SFM = spectral_flatness_band(Pxx, f, [0, min(2000,Fs/2)]);
        % HNR（自相关）
        HNRdB = estimate_hnr_autocorr(x, Fs, [70, 400]);

        feats(k).ratio_low = ratio_low;
        feats(k).ratio_mid = ratio_mid;
        feats(k).ratio_high = ratio_high;
        feats(k).rolloff95_Hz = rolloff95;
        feats(k).SFM_0_2k = SFM;
        feats(k).HNR_dB = HNRdB;
    end
end

function fRoll = spectral_rolloff_psd(Pxx, f, p)
    P = Pxx(:); F = f(:);
    total = sum(P); if total<=0, fRoll = F(1); return; end
    csum = cumsum(P);
    idx = find(csum >= p*total, 1, "first");
    if isempty(idx), fRoll = F(end); else, fRoll = F(idx); end
end

function SFM = spectral_flatness_band(Pxx, f, bandHz)
    idx = (f>=bandHz(1)) & (f<bandHz(2));
    if ~any(idx), SFM = NaN; return; end
    Px = Pxx(idx);
    gm = exp(mean(log(Px+eps)));
    am = mean(Px+eps);
    SFM = gm / am;
end

function HNRdB = estimate_hnr_autocorr(x, Fs, f0RangeHz)
    winLen = max(256, 2^nextpow2(round(0.03*Fs)));
    hop = round(0.5*winLen);
    w = hamming(winLen,"periodic");
    nFrames = max(1, 1+floor((length(x)-winLen)/hop));
    rmaxs = zeros(nFrames,1);
    lagMin = floor(Fs/f0RangeHz(2));
    lagMax = ceil(Fs/f0RangeHz(1));
    ptr=1;
    for i=1:nFrames
        seg = x(ptr:ptr+winLen-1).*w;
        ptr = ptr + hop;
        r = xcorr(seg, seg, lagMax, "coeff");
        rpos = r(lagMax+1:end);
        search = rpos(lagMin:lagMax);
        rmaxs(i) = max(search);
    end
    rmax = median(rmaxs);
    HNRdB = 10*log10(rmax / max(1e-6, 1-rmax));
end

function mid = mid_mean(a,b), a=a(~isnan(a)); b=b(~isnan(b)); mid=(mean(a)+mean(b))/2; end
function out = ternary(cond,a,b), if cond, out=a; else, out=b; end, end