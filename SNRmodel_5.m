clc
%close 
echo off;

%%実行するsimulinkモデルの指定////////////////////////////////////////////
sim('Dnote3lv3rd');

%初期値設定//////////////////////////////////////////////////////////////////
% N = 128*48*50;
% N = length(simout);                     %plot数
% osr = 128;                                %OSR
% Fs = 48000*osr;                            %サンプリング周波数(正規化したい場合は1)
% Ts = 1/Fs;
% simtime = Ts*N;  
% w_in = 2*pi*1000;                         %入力周波数10kHzと設定
                                  
fBmax = ceil(20000/Fs*N);                  % 帯域幅を設定
fBmin = ceil(20/Fs*N);    % 帯域幅を設定
suso = 20;

simfft_norm = 0;
simfftdb_norm = 0;

%simlinkの出力を取り込みFFT///////////////////////////////////////////////
window = hanning(N);       % 使用する窓関数を定義するhanning(N), blackman(N)
sim_n = simout(Offset:N+Offset-1);     % 波形データを取り込む from simlink
simfft_n = fft(sim_n .* window);% 窓関数と乗算
%simfft_n = fft(sim_n);
simfft_norm = simfft_n/(N/4);            % Fsを１から-1に正規化
simfftdb_norm = 20*log10(abs(simfft_norm)); % 演算結果の絶対値を取りdB変換    

%横軸（周波数軸）の正規化
fre = (0:length(simfftdb_norm)-1)/length(simfftdb_norm)*Fs; %周波数軸を正規化

%シミュレーション結果をプロットする/////////////////////////////////////////////////////////
figure(1); %別ウインドウに対数スケールを表示 
semilogx(fre,simfftdb_norm,'r');                %正規化対数  b:青　r:赤　m:マゼンタ　g:緑 c:シアン y:黄 w:白 k:黒
axis tight;
%axis([100 Fs 0 200]);
ylabel('gain[dB]');
xlabel('frequency[Hz]');
grid on;
title('正規化したFFT出力結果');
hold on;

%SNRの算出/////////////////////////////////////////////////////////

bandfft = simfft_norm(fBmin:fBmax);  %FFTを帯域幅分だけ切り出す
band = simfftdb_norm(fBmin:fBmax);
[no_fin, fin] = max(band(:));       %帯域幅に存在する最大値のインデックスを入力中心周波数と定義

%入力中心周波数に裾野を加えた分が帯域幅を超えないよう処理
if fin-suso > fBmin
	finmin = fin-suso;
else
	finmin = fBmin;
end

if fin+suso < fBmax
	finmax = fin+suso;
else
	finmax = fBmax;
end

%信号電力を計算
signal = bandfft(finmin:finmax);  %入力中心周波数に裾野を加えたものを入力信号と定義
signalpower = norm(signal);

%雑音電力を計算
noise = vertcat(bandfft(fBmin:finmin-1), bandfft(finmax+1:fBmax)); %入力信号+裾野を抜いた分をノイズとする
noisepower = norm(noise);

%SNRを計算
snr = 20*log10(signalpower/noisepower);