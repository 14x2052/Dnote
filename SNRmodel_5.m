clc
%close 
echo off;

%%���s����simulink���f���̎w��////////////////////////////////////////////
sim('Dnote3lv3rd');

%�����l�ݒ�//////////////////////////////////////////////////////////////////
% N = 128*48*50;
% N = length(simout);                     %plot��
% osr = 128;                                %OSR
% Fs = 48000*osr;                            %�T���v�����O���g��(���K���������ꍇ��1)
% Ts = 1/Fs;
% simtime = Ts*N;  
% w_in = 2*pi*1000;                         %���͎��g��10kHz�Ɛݒ�
                                  
fBmax = ceil(20000/Fs*N);                  % �ш敝��ݒ�
fBmin = ceil(20/Fs*N);    % �ш敝��ݒ�
suso = 20;

simfft_norm = 0;
simfftdb_norm = 0;

%simlink�̏o�͂���荞��FFT///////////////////////////////////////////////
window = hanning(N);       % �g�p���鑋�֐����`����hanning(N), blackman(N)
sim_n = simout(Offset:N+Offset-1);     % �g�`�f�[�^����荞�� from simlink
simfft_n = fft(sim_n .* window);% ���֐��Ə�Z
%simfft_n = fft(sim_n);
simfft_norm = simfft_n/(N/4);            % Fs���P����-1�ɐ��K��
simfftdb_norm = 20*log10(abs(simfft_norm)); % ���Z���ʂ̐�Βl�����dB�ϊ�    

%�����i���g�����j�̐��K��
fre = (0:length(simfftdb_norm)-1)/length(simfftdb_norm)*Fs; %���g�����𐳋K��

%�V�~�����[�V�������ʂ��v���b�g����/////////////////////////////////////////////////////////
figure(1); %�ʃE�C���h�E�ɑΐ��X�P�[����\�� 
semilogx(fre,simfftdb_norm,'r');                %���K���ΐ�  b:�@r:�ԁ@m:�}�[���^�@g:�� c:�V�A�� y:�� w:�� k:��
axis tight;
%axis([100 Fs 0 200]);
ylabel('gain[dB]');
xlabel('frequency[Hz]');
grid on;
title('���K������FFT�o�͌���');
hold on;

%SNR�̎Z�o/////////////////////////////////////////////////////////

bandfft = simfft_norm(fBmin:fBmax);  %FFT��ш敝�������؂�o��
band = simfftdb_norm(fBmin:fBmax);
[no_fin, fin] = max(band(:));       %�ш敝�ɑ��݂���ő�l�̃C���f�b�N�X����͒��S���g���ƒ�`

%���͒��S���g���ɐ���������������ш敝�𒴂��Ȃ��悤����
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

%�M���d�͂��v�Z
signal = bandfft(finmin:finmax);  %���͒��S���g���ɐ�������������̂���͐M���ƒ�`
signalpower = norm(signal);

%�G���d�͂��v�Z
noise = vertcat(bandfft(fBmin:finmin-1), bandfft(finmax+1:fBmax)); %���͐M��+����𔲂��������m�C�Y�Ƃ���
noisepower = norm(noise);

%SNR���v�Z
snr = 20*log10(signalpower/noisepower);