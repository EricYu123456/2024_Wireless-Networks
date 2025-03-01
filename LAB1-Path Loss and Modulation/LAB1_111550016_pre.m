% Task 1: Calculate Rx power-----------------------------------------------
d = input('Enter the link distance in meters: '); % Get link distance from user

% Given parameters
freq = 2.4e9; % 2.4 GHz
c = 3e8; % Speed of light
Pt = 10; % Transmit power in dBm
Gt = 1; % Transmit antenna gain
Gr = 1; % Receive antenna gain

% Convert transmit power to Watts
Pt_watts = 10^((Pt - 30)/10); 

% Calculate wavelength
lambda = c / freq;

% Calculate path loss using Friss' free space model
Lp = (lambda / (4 * pi * d))^2; 

% Calculate received power in Watts
Prx_watts = Pt_watts * Gt * Gr * Lp;

% Convert received power to dBm
Prx_dBm = 10 * log10(Prx_watts / 0.001);

% Display results
% fprintf('Link distance (d) = %g meters\n', d);
fprintf('Received power (Prx) = %.4g (Watts) = %.2f (dBm)\n', Prx_watts, Prx_dBm);

% Task 2: Channel Generation-----------------------------------------------
% Generate random channel h = a + bi following Gaussian distribution
rng(0); % Fix random seed  
a = randn; % a ~ N(0,1)
b = randn; % b ~ N(0,1)  
h = a + b*1j; % h = a + bi

% Scale h to the receiving power Prx
h_scale = sqrt(Prx_watts)/sqrt(a^2 + b^2);
h = h_scale * h;

% fprintf('Random channel h = %.4g + %.4gi\n', real(h), imag(h));
% fprintf('|h|^2 = %.10f, Prx_watts = %.10f\n', abs(h)^2, Prx_watts);
% fprintf('|h|^2 = %.4g\n', abs(h)^2);

% Task 3: Modulation-------------------------------------------------------

% Generate 300,000 random bits
msg_bits = randi([0 1], 1, 300000); 

% BPSK Modulation
bpsk_samples = 1 - 2*msg_bits; % Map bits to +1, -1
bpsk_samples = bpsk_samples / sqrt(mean(abs(bpsk_samples).^2)); % Normalize power

% QPSK Modulation 
qpsk_bits = reshape(msg_bits, 2, 150000)'; % Reshape bits into QPSK symbols
qpsk_symbols = bi2de(qpsk_bits, 'left-msb'); % Convert bits to symbol indices
qpsk_table = [1+1j, 1-1j, -1+1j, -1-1j]; % QPSK constellation mapping table
qpsk_samples = qpsk_table(qpsk_symbols+1); % QPSK modulation
qpsk_samples = qpsk_samples / sqrt(mean(abs(qpsk_samples).^2)); % Normalize power

% 16QAM Modulation
qam16_bits = reshape(msg_bits, 4, 75000)'; % Reshape bits into 16QAM symbols  
qam16_symbols = bi2de(qam16_bits, 'left-msb'); % Convert bits to symbol indices
qam16_table = ...
    [3+3j,  3+1j,  3-1j,  3-3j, ...
     1+3j,  1+1j,  1-1j,  1-3j, ...
    -1+3j, -1+1j, -1-1j, -1-3j, ...
    -3+3j, -3+1j, -3-1j, -3-3j];
qam16_samples = qam16_table(qam16_symbols+1); % 16QAM modulation  
scale = sqrt(mean(abs(qam16_samples).^2));
qam16_samples = qam16_samples / scale; % Normalize power
qam16_table = qam16_table / scale;

% 64QAM Modulation
qam64_bits = reshape(msg_bits, 6, 50000)'; % Reshape bits into 64QAM symbols
qam64_symbols = bi2de(qam64_bits, 'left-msb'); % Convert bits to symbol indices
qam64_table = [...
      7+7j,  7+5j,  7+3j,  7+1j,  7-1j,  7-3j,  7-5j,  7-7j,...
      5+7j,  5+5j,  5+3j,  5+1j,  5-1j,  5-3j,  5-5j,  5-7j,...
      3+7j,  3+5j,  3+3j,  3+1j,  3-1j,  3-3j,  3-5j,  3-7j,...
      1+7j,  1+5j,  1+3j,  1+1j,  1-1j,  1-3j,  1-5j,  1-7j,...
     -1+7j, -1+5j, -1+3j, -1+1j, -1-1j, -1-3j, -1-5j, -1-7j,...
     -3+7j, -3+5j, -3+3j, -3+1j, -3-1j, -3-3j, -3-5j, -3-7j,...
     -5+7j, -5+5j, -5+3j, -5+1j, -5-1j, -5-3j, -5-5j, -5-7j,...
     -7+7j, -7+5j, -7+3j, -7+1j, -7-1j, -7-3j, -7-5j, -7-7j];
qam64_samples = qam64_table(qam64_symbols+1); % 64QAM modulation
scale = sqrt(mean(abs(qam64_samples).^2));
qam64_samples = qam64_samples / scale; % Normalize power
qam64_table = qam64_table / scale;


% Task 4: Transmit over the Air--------------------------------------------
N0_dBm = -90; % Noise power in dBm
N0_watts = 10^((N0_dBm - 30)/10); % Convert noise power to Watts

% Generate same noise vector for different modulations
noise_vec = randn(1, length(bpsk_samples)) + randn(1, length(bpsk_samples))* 1j;
noise_vec = noise_vec / sqrt(2) * sqrt(N0_watts);
fprintf('noise power = %.2e\n', mean(abs(noise_vec).^2));
% BPSK
bpsk_received = h * bpsk_samples + noise_vec(1:length(bpsk_samples));

% QPSK
qpsk_received = h * qpsk_samples + noise_vec(1:length(qpsk_samples));

% 16QAM
qam16_received = h * qam16_samples + noise_vec(1:length(qam16_samples)); 

% 64QAM
qam64_received = h * qam64_samples + noise_vec(1:length(qam64_samples));


% Task 5: Decode and Demodulation------------------------------------------

% BPSK Demodulation
bpsk_demod = sign(real(bpsk_received ./ h));
bpsk_decoded = (-bpsk_demod + 1)/2;

% QPSK Demodulation  
qpsk_decoded_pre = (-sign(real(qpsk_received ./ h)) + 1) + (-sign(imag(qpsk_received ./ h)) + 1)/2;
qpsk_decoded = zeros(1, 2*length(qpsk_decoded_pre));
for i = 2:2:2*length(qpsk_decoded_pre)
    qpsk_decoded(i - 1) = (-sign(real(qpsk_received(i/2) ./ h)) + 1)/2;
    qpsk_decoded(i) = (-sign(imag(qpsk_received(i/2) ./ h)) + 1)/2;
end

% 16QAM Demodulation
qam16_demod = qam16_received ./ h;
qam16_symbols_hat = zeros(size(qam16_demod));

for i = 1:length(qam16_demod)
    [~, qam16_symbols_hat(i)] = min(abs(qam16_demod(i) - qam16_table));
end
qam16_bits = de2bi(qam16_symbols_hat-1, 4, 'left-msb')';
qam16_decoded = qam16_bits(:)';

% 64QAM Demodulation  
qam64_demod = qam64_received ./ h;
qam64_symbols_hat = zeros(size(qam64_demod));

for i = 1:length(qam64_demod)
    [~, qam64_symbols_hat(i)] = min(abs(qam64_demod(i) - qam64_table));
end
qam64_bits = de2bi(qam64_symbols_hat-1, 6, 'left-msb')';
qam64_decoded = qam64_bits(:)';

% Task 6: Calculate SNR-----------------------------------------------------
% BPSK
bpsk_noise_power_watts = mean(abs((bpsk_received ./ h) - bpsk_samples).^2); % Empirical noise power
bpsk_noise_power_dBm = 10*log10(bpsk_noise_power_watts/0.001); % Noise power in dBm
bpsk_SNR = mean(abs(bpsk_samples).^2) / bpsk_noise_power_watts; % SNR in Watts
bpsk_SNR_dB = 10*log10(bpsk_SNR); % SNR in dB

% QPSK
qpsk_noise_power_watts = mean(abs((qpsk_received ./ h) - qpsk_samples).^2); % Empirical noise power
qpsk_noise_power_dBm = 10*log10(qpsk_noise_power_watts/0.001); % Noise power in dBm
qpsk_SNR = mean(abs(qpsk_samples).^2) / qpsk_noise_power_watts; % SNR in Watts
qpsk_SNR_dB = 10*log10(qpsk_SNR); % SNR in dB

% 16QAM
qam16_noise_power_watts = mean(abs((qam16_received ./ h) - qam16_samples).^2); % Empirical noise power
qam16_noise_power_dBm = 10*log10(qam16_noise_power_watts/0.001); % Noise power in dBm
qam16_SNR = mean(abs(qam16_samples).^2) / qam16_noise_power_watts; % SNR in Watts
qam16_SNR_dB = 10*log10(qam16_SNR); % SNR in dB

% 64QAM
qam64_noise_power_watts = mean(abs((qam64_received ./ h) - qam64_samples).^2); % Empirical noise power
qam64_noise_power_dBm = 10*log10(qam64_noise_power_watts/0.001); % Noise power in dBm
qam64_SNR = mean(abs(qam64_samples).^2) / qam64_noise_power_watts; % SNR in Watts
qam64_SNR_dB = 10*log10(qam64_SNR); % SNR in dB

theoretical_SNR = Prx_watts / N0_watts;
fprintf('SNR : BPSK = %.4f, QPSK = %.4f, 16QAM = %.4f, 64QAM = %.4f\n', ...
    bpsk_SNR, qpsk_SNR, qam16_SNR, qam64_SNR);
fprintf('SNR (dB): BPSK = %.2f, QPSK = %.2f, 16QAM = %.2f, 64QAM = %.2f\n', ...
    bpsk_SNR_dB, qpsk_SNR_dB, qam16_SNR_dB, qam64_SNR_dB);
fprintf('theoretical_SNR = %.2e\n', theoretical_SNR);

% Task 7: Calculate Throughput---------------------------------------------
sample_duration = 3.2e-6; % 3.2 microseconds
packet_size = 500; % 500 bytes
bit_num = packet_size * 8; % Number of bits per packet

% BPSK
bpsk_ber = sum(bpsk_decoded ~= msg_bits) / length(bpsk_decoded);
bpsk_pdr_theoretical = (1 - bpsk_ber)^bit_num;

bpsk_decoded_reshape = reshape(bpsk_decoded, bit_num, []);
msg_bits_reshape = reshape(msg_bits, bit_num, []);
bpsk_check = bpsk_decoded_reshape == msg_bits_reshape;
row_products = prod(bpsk_check, 1);
bpsk_num_packets = sum(row_products);

bpsk_pdr_empirical = bpsk_num_packets / (length(bpsk_decoded) / bit_num);
bpsk_throughput_theoretical = bpsk_pdr_theoretical / sample_duration;
bpsk_throughput_empirical = bpsk_pdr_empirical / sample_duration;

% QPSK  
qpsk_ber = sum(qpsk_decoded ~= msg_bits) / length(qpsk_decoded);
qpsk_pdr_theoretical = (1 - qpsk_ber)^bit_num;

qpsk_decoded_reshape = reshape(qpsk_decoded, bit_num, []);
qpsk_check = qpsk_decoded_reshape == msg_bits_reshape;
row_products = prod(qpsk_check, 1);
qpsk_num_packets = sum(row_products);

qpsk_pdr_empirical = qpsk_num_packets / (length(qpsk_decoded) / bit_num);
qpsk_throughput_theoretical = qpsk_pdr_theoretical / sample_duration* 2;
qpsk_throughput_empirical = qpsk_pdr_empirical / sample_duration* 2;

% 16QAM
qam16_ber = sum(qam16_decoded ~= msg_bits) / length(qam16_decoded);
qam16_pdr_theoretical = (1 - qam16_ber)^bit_num;

qam16_decoded_reshape = reshape(qam16_decoded, bit_num, []);
qam16_check = qam16_decoded_reshape == msg_bits_reshape;
colum_products = prod(qam16_check, 1);
qam16_num_packets = sum(colum_products);

qam16_pdr_empirical = qam16_num_packets / (length(qam16_decoded) / bit_num);
qam16_throughput_theoretical = qam16_pdr_theoretical / sample_duration* 4;
qam16_throughput_empirical = qam16_pdr_empirical / sample_duration* 4;
% 64QAM
qam64_bitcheck = sum(qam64_decoded ~= msg_bits);
qam64_ber = sum(qam64_decoded ~= msg_bits) / length(qam64_decoded);
qam64_pdr_theoretical = (1 - qam64_ber)^bit_num;
qam64_decoded_reshape = reshape(qam64_decoded, bit_num, []);
qam64_check = qam64_decoded_reshape == msg_bits_reshape;
row_products = prod(qam64_check, 1);
qam64_num_packets = sum(row_products);
% qam64_num_packets = sum(reshape(qam64_decoded, bit_num, []) == reshape(msg_bits, bit_num, []), 2);
qam64_pdr_empirical = qam64_num_packets / (length(qam64_decoded) / bit_num);
qam64_throughput_theoretical = qam64_pdr_theoretical / sample_duration * 6;
qam64_throughput_empirical = qam64_pdr_empirical / sample_duration * 6;

fprintf('Empirical BER: BPSK = %.4f, QPSK = %.4f, 16QAM = %.4f, 64QAM = %.4f\n', ...
    bpsk_ber, qpsk_ber, qam16_ber, qam64_ber);

fprintf('Theoretical Throughput (bit/s): BPSK = %.2e, QPSK = %.2e, 16QAM = %.2e, 64QAM = %.2e\n', ...
    bpsk_throughput_theoretical, qpsk_throughput_theoretical, qam16_throughput_theoretical, qam64_throughput_theoretical);

fprintf('Empirical Throughput (bit/s): BPSK = %.2e, QPSK = %.2e, 16QAM = %.2e, 64QAM = %.2e\n', ...
    bpsk_throughput_empirical, qpsk_throughput_empirical, qam16_throughput_empirical, qam64_throughput_empirical);