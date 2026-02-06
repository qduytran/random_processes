clc; clear; close all;

%% Setup parameter for simulation
fs = 1000;              % frequency sampling
T = 2;                  % time range to plot
t = 0:1/fs:T-1/fs;      
N = length(t);  

s = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);

rng(42);
noise_power = 2.5;
v = sqrt(noise_power) * randn(size(t));

x = s + v;

%% Plot signal and noise
hFig1 = figure('Name', 'Analysis of Signal Components', 'Color', 'w');

t_zoom = t(t <= 0.2);
s_zoom = s(t <= 0.2);
v_zoom = v(t <= 0.2);
x_zoom = x(t <= 0.2);

% original signal (50 Hz and 120 Hz)
subplot(3,1,1);
plot(t_zoom, s_zoom, 'b', 'LineWidth', 1.5);
title('1. Original signal s[n] (50Hz + 120Hz)');
ylabel('Amplititude'); grid on;

% Gaussian Noise
subplot(3,1,2);
plot(t_zoom, v_zoom, 'r');
title(['2. Gaussian noise v[n] (Variance = ', num2str(noise_power), ')']);
ylabel('Amplititude'); grid on;

% Signal + noise
subplot(3,1,3);
plot(t_zoom, x_zoom, 'k');
hold on;
plot(t_zoom, s_zoom, 'b--', 'LineWidth', 1); % Vẽ đè tín hiệu gốc để so sánh
title('3. Observation signal x[n] = s[n] + v[n]');
xlabel('Time (s)'); ylabel('Amplititude');
legend('Real signal', 'Sin signal');
grid on;

exportgraphics(hFig1, 'Plot_signal.png', 'Resolution', 300);
exportgraphics(hFig1, 'Plot_signal.pdf', 'ContentType', 'vector');
%% Analysis 2 methods
[pxx_per, f_per] = periodogram(x, rectwin(N), N, fs);

window_len = 256;         % Window length
overlap = window_len/2;   % Overlap
nfft = 1024;              % number of FFT

[pxx_welch, f_welch] = pwelch(x, hamming(window_len), overlap, nfft, fs);

%% plot 
hFig2 = figure('Name', 'Periodogram va Welch', 'Color', 'w');

subplot(3,1,1);
plot(t(1:300), x(1:300), 'k'); 
title('Time-domain observation signal: x[n] = s[n] + v[n]');
xlabel('Time (s)'); ylabel('Amplititude');
grid on;

subplot(3,1,2);
plot(f_per, 10*log10(pxx_per), 'Color', [0.6 0.6 0.6]); 
title('Analysis 1: Periodogram (Inconsistent Estimation)');
xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
grid on; xlim([0 200]); ylim([-30 15]);
text(50, 5, '\leftarrow Strong signal interference', 'Color', 'r');

subplot(3,1,3);
plot(f_welch, 10*log10(pxx_welch), 'b', 'LineWidth', 1.5); % Màu xanh đậm
title(['Analysis 2: Welch Method (Hamming window, L=', num2str(window_len), ')']);
xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
grid on; xlim([0 200]); ylim([-30 15]);

[pks, locs] = findpeaks(10*log10(pxx_welch), f_welch, 'MinPeakHeight', -5);
text(locs, pks+2, num2str(locs, '%.1f Hz'), 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

exportgraphics(hFig2, 'Welch_analysis.png', 'Resolution', 300);
exportgraphics(hFig2, 'Welch_analysis.pdf', 'ContentType', 'vector');
%% calculate SD from 300 - 450 Hz
idx_noise = (f_per > 300 & f_per < 450);
std_per = std(10*log10(pxx_per(idx_noise)));

idx_noise_w = (f_welch > 300 & f_welch < 450);
std_welch = std(10*log10(pxx_welch(idx_noise_w)));

fprintf('Standard deviation of the Periodogram: %.2f dB\n', std_per);
fprintf('Standard deviation of the Welch: %.2f dB\n', std_welch);
fprintf('Welch helps reduce spectral noise %.1f times!\n', std_per/std_welch);

%% Analysis of stable
freq_signal = [50, 120];
freq_noise_range = [200 450]; 

idx_n_per = (f_per >= freq_noise_range(1) & f_per <= freq_noise_range(2));
idx_n_wel = (f_welch >= freq_noise_range(1) & f_welch <= freq_noise_range(2));

% Power in average noise floor
noise_floor_per = mean(pxx_per(idx_n_per));
noise_floor_wel = mean(pxx_welch(idx_n_wel));

% Power in peak Power (120Hz)
peak_per = max(pxx_per(f_per > 115 & f_per < 125));
peak_wel = max(pxx_welch(f_welch > 115 & f_welch < 125));

% SNR
snr_per = 10*log10(peak_per / noise_floor_per);
snr_wel = 10*log10(peak_wel / noise_floor_wel);

% Variance/Std of Noise
std_noise_per = std(10*log10(pxx_per(idx_n_per)));
std_noise_wel = std(10*log10(pxx_welch(idx_n_wel)));

fprintf('====================================================\n');
fprintf('   QUANTITATIVE ANALYSIS RESULTS (PERIODOGRAM vs WELCH)\n');
fprintf('====================================================\n');
fprintf('%-25s | %-12s | %-12s\n', 'Parameterization', 'Periodogram', 'Welch');
fprintf('----------------------------------------------------\n');
fprintf('%-25s | %-12.2f | %-12.2f\n', 'SNR in 120Hz (dB)', snr_per, snr_wel);
fprintf('%-25s | %-12.2f | %-12.2f\n', 'Noise standard deviation (dB)', std_noise_per, std_noise_wel);
fprintf('----------------------------------------------------\n');

if std_noise_per > std_noise_wel
    fprintf('=> The Welch method reduces spectral noise %.1f times more effectively than the Periodogram.\n', ...
        std_noise_per / std_noise_wel);
end

%% ACF
[acf, lags] = xcorr(x, 'biased'); 
lags_ms = lags / fs; 

hFig3 = figure('Name', 'Autocorrelation and Stationarity', 'Color', 'w');
tlo2 = tiledlayout(2, 1, 'Padding', 'normal', 'TileSpacing', 'compact');

nexttile;
plot(lags_ms, acf, 'g', 'LineWidth', 1.2);
xlim([-0.1 0.1]); 
title('1. Autocorrelation Function (ACF) of x[n]');
xlabel('Lag (s)'); ylabel('Amplitude');
grid on;
legend('ACF shows the periodicity of 50Hz & 120Hz');

%% check WSS
x1 = x(1:N/2);
x2 = x(N/2+1:end);

[p_part1, f_part] = pwelch(x1, hamming(window_len), overlap, nfft, fs);
[p_part2, ~]      = pwelch(x2, hamming(window_len), overlap, nfft, fs);

nexttile;
plot(f_part, 10*log10(p_part1), 'b', 'LineWidth', 1); hold on;
plot(f_part, 10*log10(p_part2), 'r--', 'LineWidth', 1);
xlim([0 200]); ylim([-30 15]);
title('2. Stationarity Check: Comparing PSD of two different time segments');
xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
legend('Segment 1 (0-1s)', 'Segment 2 (1-2s)');
grid on;

title(tlo2, 'Stochastic Process Analysis: ACF and Stationarity', 'FontSize', 12, 'FontWeight', 'bold');
exportgraphics(hFig3, 'Stochastic_Analysis.png', 'Resolution', 300);
exportgraphics(hFig3, 'Stochastic_Analysis.pdf', 'ContentType', 'vector');