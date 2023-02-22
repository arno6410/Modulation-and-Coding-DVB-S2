%% Modulation and Coding - Project
% DVB-S2 Simulation

%% 1) Communication chain sim
% simulation of the optimal communication chain over the ideal channel

%% a) Symbol mapping and demapping
close all; clc; clear;
test_vector = [0; 0; 0; 1; 1; 1; 1; 0];
Nbps = 2; %number of bits per symbol
result = mapping(test_vector, Nbps, 'pam'); %4-PAM ex
result2 = mapping(test_vector, Nbps, 'qam'); %4-QAM ex

figure
scatter(result,zeros(size(result)))

test = demapping(result, Nbps, 'pam')
test2 = demapping(result2, Nbps, 'qam')

test_vector2 = [0; 0; 1; 0; 1; 0; 1; 1; 1];
Nbps2 = 3;
result3 = mapping(test_vector2, Nbps2, 'qam');

% figure
% plot(result2,'x')
% xL = xlim;
% yL = ylim;
% line([0 0], yL);  %x-axis
% line(xL, [0 0]);  %y-axis