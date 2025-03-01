close all; clear; clc;
addpath ./tasks;

origin = [0, 0];
tx_location = origin;
P_tx_dBm = 10;          % Transmission power of Tx (dBm)
N0_dBm = -95;           % Assume noise power is -90 dBm

tx_node_number = 1;          % Number of Tx users
rx_node_number = 2;          % Number of Rx users
digital_antenna_number = 2;  % Number of Tx antennas of digital
rx_antenna_number = 1;       % Number of Rx antennas

% q1
distances = 50:50:500;
tx_nums = tx_node_number * digital_antenna_number;
rx_nums = rx_node_number * rx_antenna_number;
avg_snr_zf = zeros(length(distances), 2);
avg_snr_no_zf = zeros(length(distances), 2);

for d_idx = 1:length(distances)
    d = distances(d_idx);
    snr_zf_arr = zeros(10, 2);
    snr_no_zf_arr = zeros(10, 2);
    for topo_idx = 1:10
        numbers = 0:10:180;
        random_index1 = randi(length(numbers));
        random_index2 = randi(length(numbers));
        random_number1 = numbers(random_index1);
        random_number2 = numbers(random_index2);

        offset = -5 + 10*rand();
        rx1_x = d * cosd(random_number1 + offset);
        rx1_y = d * sind(random_number1 + offset);
        rx1_location = [rx1_x, rx1_y];
        
        offset = -5 + 10*rand();
        rx2_x = d * cosd(random_number2 + offset);
        rx2_y = d * sind(random_number2 + offset);
        rx2_location = [rx2_x, rx2_y];
        
        [snr1, snr2, snr1_no, snr2_no, ~, ~] = digital_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, tx_nums, rx_nums);
        snr_zf_arr(topo_idx, :) = [snr1, snr2];
        snr_no_zf_arr(topo_idx, :) = [snr1_no, snr2_no];
    end
    avg_snr_zf(d_idx, :) = mean(snr_zf_arr, 1);
    avg_snr_no_zf(d_idx, :) = mean(snr_no_zf_arr, 1);
end

figure;
plot(distances, avg_snr_zf(:,1), '-o', distances, avg_snr_no_zf(:,1), '--*');
hold on;
plot(distances, avg_snr_zf(:,2), '-s', distances, avg_snr_no_zf(:,2), '--d');
legend('R1 with ZFBF', 'R1 without ZFBF', 'R2 with ZFBF', 'R2 without ZFBF');
xlabel('Distance (m)');
ylabel('Average SNR (dBm)');
title('Average SNR with and without ZFBF');


% q2
d = 200;
tx_nums = tx_node_number * digital_antenna_number;
rx_nums = rx_node_number * rx_antenna_number;
heq_arr = zeros(1, 10);
error_arr = zeros(1, 10);

for topo_idx = 1:10
    numbers = 0:10:180;
    random_index1 = randi(length(numbers));
    random_index2 = randi(length(numbers));
    random_number1 = numbers(random_index1);
    random_number2 = numbers(random_index2);

    offset = -5 + 10*rand();
    rx1_x = d * cosd(random_number1 + offset);
    rx1_y = d * sind(random_number1 + offset);
    rx1_location = [rx1_x, rx1_y];
    
    offset = -5 + 10*rand();
    rx2_x = d * cosd(random_number2 + offset);
    rx2_y = d * sind(random_number2 + offset);
    rx2_location = [rx2_x, rx2_y];
    
    [~, ~, ~, ~, H_eq_1, x_err_power1_dBm] = digital_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, tx_nums, rx_nums);
    heq_arr(topo_idx) = H_eq_1;
    error_arr(topo_idx) = x_err_power1_dBm;
end

figure;
subplot(1,2,1);
plot(1:10, abs(heq_arr), '-o');
xlabel('Topology Index');
ylabel('|h_e_q|');
title('Magnitude of h_e_q for 10 Topologies, d=200m');

subplot(1,2,2);
plot(1:10, error_arr, '-o');
xlabel('Topology Index'); 
ylabel('Error (dBm)');
title('Error for 10 Topologies, d=200m');