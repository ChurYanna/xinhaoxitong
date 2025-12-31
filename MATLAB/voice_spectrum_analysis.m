%% voice_spectrum_analysis.m
% 对四段语音wav做频谱分析并对比：单边幅度谱、Welch功率谱密度、语谱图等
clc; clear; close all;

%% 1) 文件路径与文件名
baseDir = "C:\Users\Administrator\Desktop\学业备份\大三上\信号与系统\《信号与系统》期末大作业\WAV\";

files = ["voice_010.wav", "voice_020.wav", "voice_030.wav", "voice_040.wav"];
labels = ["voice\_010", "voice\_020", "voice\_030", "voice\_040"];

nFiles = numel(files);

% 一些绘图颜色
colors = lines(nFiles);

%% 2) 读取音频并预处理（转单声道、去直流、归一化）
x = cell(nFiles,1);
Fs = zeros(nFiles,1);

for k = 1:nFiles
    fpath = baseDir + files(k);
    [sig, fs] = audioread(fpath);
    Fs(k) = fs;

    % 如果是双声道，转为单声道（取均值）
    if size(sig,2) > 1
        sig = mean(sig, 2);
    end

    % 去直流分量
    sig = sig - mean(sig);

    % 归一化（避免幅度差导致比较不直观）
    sig = sig / (max(abs(sig)) + eps);

    x{k} = sig;

    fprintf("[%s] Fs = %d Hz, length = %d samples, duration = %.3f s\n", ...
        labels(k), fs, length(sig), length(sig)/fs);
end

%% 3) 参数设置（FFT、Welch、语谱图）
NFFT = 2^nextpow2(max(cellfun(@length, x)));   % 统一NFFT，便于比较
% 也可以改成固定值，比如 NFFT = 8192;

% Welch参数（建议用于语音：更平滑稳定）
welchWinLen = 1024;
welchWin = hamming(welchWinLen, "periodic");
welchOverlap = round(0.5 * welchWinLen);
welchNFFT = 4096;

% 语谱图参数
specWinLen = 512;
specWin = hamming(specWinLen, "periodic");
specOverlap = round(0.75 * specWinLen);
specNFFT = 1024;

%% 4) 每个音频：时域 + 单边幅度谱 + PSD + 语谱图
for k = 1:nFiles
    sig = x{k};
    fs = Fs(k);
    t = (0:length(sig)-1)/fs;

    figure("Name", "Analysis - " + labels(k), "Color", "w", "Position", [100, 100, 1200, 800]);

    % ---- (a) 时域波形
    subplot(2,2,1);
    plot(t, sig, "Color", colors(k,:), "LineWidth", 1);
    grid on;
    xlabel("Time (s)");
    ylabel("Amplitude (norm.)");
    title("Time Domain Waveform: " + labels(k));

    % ---- (b) 单边幅度谱（Amplitude Spectrum, Single-sided）
    % FFT
    L = length(sig);
    X = fft(sig, NFFT);
    P2 = abs(X)/L;                 % 双边幅度谱（幅值）
    P1 = P2(1:NFFT/2+1);           % 单边
    P1(2:end-1) = 2*P1(2:end-1);   % 单边补偿
    f = fs*(0:(NFFT/2))/NFFT;

    subplot(2,2,2);
    plot(f, 20*log10(P1 + eps), "Color", colors(k,:), "LineWidth", 1);
    grid on;
    xlabel("Frequency (Hz)");
    ylabel("Magnitude (dB)");
    title("Single-sided Amplitude Spectrum (dB)");
    xlim([0, min(8000, fs/2)]);  % 语音常用观察0~8k（可按需要调整）

    % 主峰频率（简单找最大峰）
    [~, idxPeak] = max(P1);
    fPeak = f(idxPeak);
    hold on;
    plot(fPeak, 20*log10(P1(idxPeak)+eps), "ro", "MarkerSize", 6, "LineWidth", 1.2);
    text(fPeak, 20*log10(P1(idxPeak)+eps), sprintf("  Peak %.1f Hz", fPeak), ...
        "Color", "r", "FontSize", 10, "VerticalAlignment", "bottom");

    % ---- (c) Welch功率谱密度 PSD
    subplot(2,2,3);
    [Pxx, fPxx] = pwelch(sig, welchWin, welchOverlap, welchNFFT, fs, "onesided");
    plot(fPxx, 10*log10(Pxx + eps), "Color", colors(k,:), "LineWidth", 1);
    grid on;
    xlabel("Frequency (Hz)");
    ylabel("PSD (dB/Hz)");
    title("Power Spectral Density (Welch)");
    xlim([0, min(8000, fs/2)]);

    % ---- (d) 语谱图 Spectrogram
    subplot(2,2,4);
    [S,F,T] = spectrogram(sig, specWin, specOverlap, specNFFT, fs, "yaxis");
    imagesc(T, F, 20*log10(abs(S)+eps));
    axis xy;
    colormap jet;
    colorbar;
    xlabel("Time (s)");
    ylabel("Frequency (Hz)");
    title("Spectrogram (dB)");
    ylim([0, min(8000, fs/2)]);
end

%% 5) 四段音频的“叠加对比图”（单边幅度谱 + PSD）
% 如果采样率不同，会导致频轴上限不同；这里用各自的频轴，但统一显示到8k或fs/2
figure("Name","Comparison - Spectra Overlay","Color","w","Position",[120,120,1200,500]);

% ---- 叠加幅度谱
subplot(1,2,1); hold on; grid on;
for k = 1:nFiles
    sig = x{k}; fs = Fs(k);
    L = length(sig);
    X = fft(sig, NFFT);
    P2 = abs(X)/L;
    P1 = P2(1:NFFT/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(NFFT/2))/NFFT;

    plot(f, 20*log10(P1 + eps), "LineWidth", 1.2, "Color", colors(k,:));
end
xlabel("Frequency (Hz)");
ylabel("Magnitude (dB)");
title("Overlay: Single-sided Amplitude Spectrum");
legend(labels, "Interpreter","none", "Location","best");
xlim([0, min(8000, min(Fs)/2)]); % 用最小Fs保证大家都有该范围

% ---- 叠加PSD
subplot(1,2,2); hold on; grid on;
for k = 1:nFiles
    sig = x{k}; fs = Fs(k);
    [Pxx, fPxx] = pwelch(sig, welchWin, welchOverlap, welchNFFT, fs, "onesided");
    plot(fPxx, 10*log10(Pxx + eps), "LineWidth", 1.2, "Color", colors(k,:));
end
xlabel("Frequency (Hz)");
ylabel("PSD (dB/Hz)");
title("Overlay: PSD (Welch)");
legend(labels, "Interpreter","none", "Location","best");
xlim([0, min(8000, min(Fs)/2)]);

disp("Done. 已生成每段音频单独分析图 + 叠加对比图。");