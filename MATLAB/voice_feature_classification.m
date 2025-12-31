%% voice_feature_classification.m
% 任务（第六点）：频谱特征提取与分类（R2024a 优化）
% 数据：voice_010（原唱中文）、voice_020（原唱粤语）、voice_030（安静翻唱）
%       可选：voice_040（嘈杂环境，参与特征但不参与“原唱 vs 翻唱”分类）
% 特征 Z：
% - ratio_high（4–8 kHz 高频能量占比）
% - rolloff95_Hz（95% 频谱能量滚降频率）
% - SFM_0_2k（0–2 kHz 谱平坦度）
% - HNR_dB（谐波噪声比，源周期性指标，若有 Audio Toolbox 也会尝试 pitch 轨迹）
%
% 输出：
% - STFT 语谱图（如可用叠加 pitch）PNG
% - 单文件 Welch PSD PNG
% - 频带能量比例堆叠柱状图 PNG
% - results_features.csv（包含各特征与分类结果）
%
% 说明：
% - 分类只对 010/020/030 进行（原唱 vs 翻唱）
% - 脚本会自动检测 Audio Toolbox（pitch/harmonicRatio），没有则回退到内置实现
%
% 作者：你（报告用），协作：Copilot（R2024a优化）

clc; clear; close all;

%% 1) 文件路径与文件名
baseDir = "C:\Users\Administrator\Desktop\学业备份\大三上\信号与系统\《信号与系统》期末大作业\WAV\";
files = {
    "voice_010.wav", ... % 原唱中文
    "voice_020.wav", ... % 原唱粤语
    "voice_030.wav", ... % 安静翻唱
    "voice_040.wav"  ... % 嘈杂环境（可选）
};
labels = {
    "voice_010 (Mandarin, original)", ...
    "voice_020 (Cantonese, original)", ...
    "voice_030 (quiet cover)", ...
    "voice_040 (noisy)"
};

% 分类只针对前三段（010/020/030）
classifyIdx = [1,2,3];     % 参与分类
originalIdx = [1,2];       % 原唱类
coverIdx    = [3];         % 翻唱类

outDir = fullfile(pwd, "feature_outputs");
if ~exist(outDir, "dir"), mkdir(outDir); end

%% 2) 工具箱与函数可用性检测
hasAudioTB = license('test','Audio_Toolbox') && exist('pitch','file')==2;
hasSTFT    = exist('stft','file')==2;  % R2024a Signal Processing Toolbox 有 stft

%% 3) 参数（R2024a优化）
% STFT：窗口 ~25 ms，75% 重叠，自适应 Fs
stftTargetMs = 25;
stftOverlap  = 0.75;

% Welch PSD（稳健、报告常用）
welchWinLenFactor = 0.046;   % ~46 ms 窗长（随 Fs 调整）
welchOverlap      = 0.5;
welchNFFTMin      = 4096;

% 频带划分（报告用）
bands = struct( ...
    "low",   [0, 1000], ...     % 0–1 kHz
    "mid",   [1000, 4000], ...  % 1–4 kHz
    "high",  [4000, 8000] ...   % 4–8 kHz
);

%% 4) 读取/预处理/特征提取
N = numel(files);
allFeats = table();

for k = 1:N
    % 读取
    fpath = baseDir + files{k};
    if ~isfile(fpath)
        warning("找不到文件：%s（跳过）", fpath);
        continue;
    end
    [x, Fs] = audioread(fpath);

    % 单声道、去直流、归一化
    if size(x,2) > 1, x = mean(x,2); end
    x = x - mean(x);
    x = x / (max(abs(x)) + eps);
    durSec = numel(x)/Fs;

    %% 4.1 STFT（优先 stft，回退 spectrogram）
    % 自适应窗口长度
    stftWinLen = max(256, 2^nextpow2(round(stftTargetMs*1e-3*Fs)));
    stftOverlapLen = round(stftOverlap * stftWinLen);
    stftNFFT = max(stftWinLen, 2^nextpow2(stftWinLen)); % 至少等于窗口长度
    w = hamming(stftWinLen, "periodic");

    if hasSTFT
        % R2024a stft
        [S, F, T] = stft(x, Fs, ...
            Window = w, ...
            OverlapLength = stftOverlapLen, ...
            FFTLength = stftNFFT);
        S_dB = 20*log10(abs(S)+eps);

        % 如可用 pitch，估计并叠加轨迹（便于展示基频）
        if hasAudioTB
            f0 = pitch(x, Fs, ...
                Method = "ACF", ...
                WindowLength = stftWinLen, ...
                OverlapLength = stftOverlapLen, ...
                Range = [70, 400]);       % 语音常见基频范围
            tPitch = linspace(0, durSec, numel(f0)).';  % 对齐时间
        else
            f0 = []; tPitch = [];
        end

        figSTFT = figure("Color","w","Position",[60,60,1200,500]);
        tiledlayout(figSTFT,1,1,"TileSpacing","compact","Padding","compact");
        ax = nexttile;
        imagesc(ax, T, F, S_dB); axis(ax,"xy"); colormap(ax,"jet"); colorbar(ax);
        xlabel(ax,"Time (s)"); ylabel(ax,"Frequency (Hz)");
        title(ax, sprintf("Spectrogram (STFT): %s", labels{k}));
        ylim(ax,[0, min(8000, Fs/2)]);

        % 动态范围自适应（更稳的可视化）
        climLow = prctile(S_dB(:), 5);
        climHigh= prctile(S_dB(:), 95);
        set(ax,"CLim",[climLow, climHigh]);

        % 叠加 pitch（若可用）
        hold(ax,"on");
        if ~isempty(f0)
            plot(ax, tPitch, f0, "w-", "LineWidth", 1.2);
            legend(ax, ["STFT", "Pitch"], "Location","northeast");
        end

        exportgraphics(figSTFT, fullfile(outDir, sprintf("spectrogram_%s.png", erase(files{k}, ".wav"))), "Resolution", 160);
        close(figSTFT);

    else
        % 回退 spectrogram（兼容所有版本）
        [S, F, T] = spectrogram(x, w, stftOverlapLen, stftNFFT, Fs, "yaxis");
        S_dB = 20*log10(abs(S)+eps);

        figSpec = figure("Color","w","Position",[60,60,1200,500]);
        imagesc(T, F, S_dB); axis xy; colormap jet; colorbar;
        xlabel("Time (s)"); ylabel("Frequency (Hz)");
        title(sprintf("Spectrogram (STFT fallback: spectrogram): %s", labels{k}));
        ylim([0, min(8000, Fs/2)]);
        climLow = prctile(S_dB(:), 5);
        climHigh= prctile(S_dB(:), 95);
        set(gca,"CLim",[climLow, climHigh]);
        exportgraphics(figSpec, fullfile(outDir, sprintf("spectrogram_%s.png", erase(files{k}, ".wav"))), "Resolution", 160);
        close(figSpec);
    end

    %% 4.2 Welch PSD（稳健能量估计）
    welchWinLen = max(256, 2^nextpow2(round(welchWinLenFactor*Fs)));
    welchNFFT   = max(welchNFFTMin, 2^nextpow2(welchWinLen));
    wWelch = hamming(welchWinLen, "periodic");
    noWelch = round(welchOverlap * welchWinLen);
    [Pxx, fPxx] = pwelch(x, wWelch, noWelch, welchNFFT, Fs, "onesided");
    Pxx_dB = 10*log10(Pxx + eps);

    % PSD 图保存
    figPSD = figure("Color","w","Position",[90,90,900,420]);
    plot(fPxx, Pxx_dB, "LineWidth",1.2);
    grid on; xlabel("Frequency (Hz)"); ylabel("PSD (dB/Hz)");
    title(sprintf("Welch PSD: %s", labels{k}));
    xlim([0, min(8000, Fs/2)]);
    exportgraphics(figPSD, fullfile(outDir, sprintf("psd_%s.png", erase(files{k}, ".wav"))), "Resolution", 160);
    close(figPSD);

    %% 4.3 频带能量比例（特征1：ratio_high 等）
    % 使用 bandpower(Pxx,f,'psd') 计算更稳的频带能量（线性能量）
    E_low  = bandpower(Pxx, fPxx, bands.low,  "psd");
    E_mid  = bandpower(Pxx, fPxx, bands.mid,  "psd");
    E_high = bandpower(Pxx, fPxx, bands.high, "psd");
    E_total = E_low + E_mid + E_high + eps;

    ratio_low  = E_low  / E_total;
    ratio_mid  = E_mid  / E_total;
    ratio_high = E_high / E_total;

    %% 4.4 高频滚降（特征2：rolloff95）
    rolloff95_Hz = spectral_rolloff_psd(Pxx, fPxx, 0.95);

    %% 4.5 谱平坦度（特征3：SFM_0_2k）
    sfmBand = [0, min(2000, Fs/2)];
    SFM_0_2k = spectral_flatness_band(Pxx, fPxx, sfmBand);

    %% 4.6 HNR 估计（特征4：HNR_dB）
    if hasAudioTB && exist('harmonicRatio','file')==2
        % Audio Toolbox 的 harmonicRatio（若可用）
        HNR_dB = 10*log10(harmonicRatio(x, Fs) + eps); % harmonicRatio 返回比值，取对数转 dB
    else
        % 回退到自相关法的 HNR
        HNR_dB = estimate_hnr_autocorr(x, Fs, [70, 400]);
    end

    %% 特征汇总
    row = table( ...
        string(files{k}), string(labels{k}), Fs, durSec, ...
        E_low, E_mid, E_high, E_total, ...
        ratio_low, ratio_mid, ratio_high, ...
        rolloff95_Hz, SFM_0_2k, HNR_dB, ...
        'VariableNames', ["file","label","Fs","duration_s","E_low","E_mid","E_high","E_total", ...
                          "ratio_low","ratio_mid","ratio_high","rolloff95_Hz","SFM_0_2k","HNR_dB"]);
    allFeats = [allFeats; row]; %#ok<AGROW>
end

%% 5) 二类分类：原唱(010/020) vs 翻唱(030)
origMask  = ismember(allFeats.file, string(files(originalIdx)));
coverMask = ismember(allFeats.file, string(files(coverIdx)));

% 自适应阈值：用原唱/翻唱组均值中点
thr_HF    = mid_mean(allFeats.ratio_high(origMask), allFeats.ratio_high(coverMask));
thr_roll  = mid_mean(allFeats.rolloff95_Hz(origMask), allFeats.rolloff95_Hz(coverMask));

% 规则：满足任一条件视为“原唱”（更稳健）
isOriginal = (allFeats.ratio_high >= thr_HF) | (allFeats.rolloff95_Hz >= thr_roll);

% 标注分类（仅 010/020/030）
classLabel = strings(height(allFeats),1);
for i = 1:height(allFeats)
    if ismember(allFeats.file(i), string(files(classifyIdx)))
        classLabel(i) = ternary(isOriginal(i), "original", "cover");
    else
        classLabel(i) = "n/a"; % 040 不参与该二类任务
    end
end

explain = strings(height(allFeats),1);
for i = 1:height(allFeats)
    explain(i) = sprintf("HF_ratio=%.4f (thr=%.4f), rolloff95=%.1f Hz (thr=%.1f)", ...
        allFeats.ratio_high(i), thr_HF, allFeats.rolloff95_Hz(i), thr_roll);
end

allFeats.result_class   = classLabel;
allFeats.result_explain = explain;

%% 6) 导出 CSV 与打印汇总
outCsv = fullfile(outDir, "results_features_2024a.csv");
writetable(allFeats, outCsv, "WriteMode","overwrite");
fprintf("\n已导出特征与分类结果 -> %s\n\n", outCsv);

disp("=== Classification (original vs cover) on 010/020/030 ===");
for i = 1:height(allFeats)
    if ismember(allFeats.file(i), string(files(classifyIdx)))
        fprintf("%s -> %s | %s\n", allFeats.file(i), allFeats.result_class(i), allFeats.result_explain(i));
    end
end

%% 7) 频带能量比例堆叠图（报告用）
figBar = figure("Color","w","Position",[120,120,1020,440]);
barData = [allFeats.ratio_low, allFeats.ratio_mid, allFeats.ratio_high];
bar(categorical(allFeats.file), barData, "stacked");
ylabel("Energy ratio"); ylim([0,1]); grid on;
legend("0–1k","1–4k","4–8k","Location","eastoutside");
title("Band Energy Ratios per file (R2024a)");
exportgraphics(figBar, fullfile(outDir, "band_energy_ratios_2024a.png"), "Resolution",160);
close(figBar);

%% ========= 本地函数（R2024a 兼容） =========
function fRoll = spectral_rolloff_psd(Pxx, f, p)
    % 频谱能量累计达到 p（0~1）时的频率
    P = Pxx(:); F = f(:);
    total = sum(P);
    if total<=0, fRoll = F(1); return; end
    csum = cumsum(P);
    idx = find(csum >= p*total, 1, "first");
    if isempty(idx), fRoll = F(end); else, fRoll = F(idx); end
end

function SFM = spectral_flatness_band(Pxx, f, bandHz)
    idx = (f >= bandHz(1)) & (f < bandHz(2));
    if ~any(idx), SFM = NaN; return; end
    Px = Pxx(idx);
    gm = exp(mean(log(Px + eps)));  % 几何均值
    am = mean(Px + eps);            % 算术均值
    SFM = gm / am;                  % 0~1，越小越“有谐波”
end

function HNRdB = estimate_hnr_autocorr(x, Fs, f0RangeHz)
    % 基于归一化自相关的 HNR 估计（无需 Audio Toolbox）
    winLen = max(256, 2^nextpow2(round(0.03*Fs))); % ~30 ms
    hop    = round(0.5*winLen);
    w      = hamming(winLen, "periodic");
    nFrames= max(1, 1 + floor((length(x)-winLen)/hop));
    rmaxs  = zeros(nFrames,1);
    lagMin = floor(Fs/f0RangeHz(2));
    lagMax = ceil(Fs/f0RangeHz(1));

    ptr = 1;
    for i = 1:nFrames
        seg = x(ptr:ptr+winLen-1) .* w;
        ptr = ptr + hop;
        r = xcorr(seg, seg, lagMax, "coeff"); % 归一化自相关
        rpos = r(lagMax+1:end);
        search = rpos(lagMin:lagMax);
        rmaxs(i) = max(search);
    end
    rmax = median(rmaxs);                     % 稳健统计
    HNRdB = 10*log10(rmax / max(1e-6, 1-rmax));
end

function mid = mid_mean(a, b)
    % 两组均值的中点（自适应阈值）
    a = a(~isnan(a)); b = b(~isnan(b));
    mid = (mean(a) + mean(b))/2;
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end