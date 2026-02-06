clc; clear; close all;

%% simulation data
fs = 1000;
t = 0:1/fs:4-1/fs; 

s = sin(2*pi*150*t) + sin(2*pi*155*t); 
v = 1.5 * randn(size(t)); 
x = s + v;

%% Estimation 
nfft = 2048;

% Low Variance, High Bias
L1 = 64; 
[p1, f1] = pwelch(x, hamming(L1), L1/2, nfft, fs);

% Medium window
L2 = 256;
[p2, f2] = pwelch(x, hamming(L2), L2/2, nfft, fs);

% High Variance, Low Bias
L3 = 2048;
[p3, f3] = pwelch(x, hamming(L3), L3/2, nfft, fs);

%% plot
hFig = figure('Color', 'w', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
tlo = tiledlayout(3, 1, 'Padding', 'normal', 'TileSpacing', 'compact');

nexttile;
plot(f1, 10*log10(p1), 'r', 'LineWidth', 1.5);
grid on; xlim([0 500]); ylim([-50 25]);
title(['Short window (L=', num2str(L1), '): LOW variance but HIGH Bias'], 'FontSize', 11);
ylabel('PSD (dB/Hz)');

nexttile;
plot(f2, 10*log10(p2), 'b', 'LineWidth', 1.5);
grid on; xlim([0 500]); ylim([-50 25]);
title(['Medium window (L=', num2str(L2), '): Balanced Bias-Variance'], 'FontSize', 11);
ylabel('PSD (dB/Hz)');

nexttile;
plot(f3, 10*log10(p3), 'k', 'LineWidth', 0.5);
grid on; xlim([0 500]); ylim([-50 25]);
title(['Long window (L=', num2str(L3), '): LOW bias but HIGH variance'], 'FontSize', 11);
xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');

title(tlo, 'Bias-Variance Trade-off Analysis', 'FontSize', 14, 'FontWeight', 'bold');set(findall(hFig,'-property','FontSize'),'FontSize',10); 
hFig.PaperPositionMode = 'auto';
exportgraphics(hFig, 'Bias_Variance_Result.png', 'Resolution', 300);
exportgraphics(hFig, 'Bias_Variance_Result.pdf', 'ContentType', 'vector');