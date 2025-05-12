% Main Program to Execute Copula Functions
clear;
clc;

% Define common parameters
Hauptpfad = 'C:\Path\To\Output'; % Update this path to your desired output directory
Dim = [5, 10, 15, 20]; % Dimensions
alpha = [0.90, 0.925, 0.95, 0.975, 0.99, 0.999]; % Confidence levels
N = 1000; % Number of paths
Size = [1000, 2000, 5000, 10000, 15000, 20000, 25000]; % Sample sizes

% Ensure the output directory exists
if ~exist(Hauptpfad, 'dir')
    mkdir(Hauptpfad);
end

% Execute Frank Copula
disp('Executing Frank Copula...');
FrankCopula(Hauptpfad, Dim, alpha, N, Size);

% Execute Gaussian Copula
disp('Executing Gaussian Copula...');
GaussianCopula(Hauptpfad, Dim, alpha, N, Size);

% Execute Gumbel Copula
disp('Executing Gumbel Copula...');
GumbelCopula(Hauptpfad, Dim, alpha, N, Size);

% Execute t-Copula
disp('Executing t-Copula...');
TCopula(Hauptpfad, Dim, alpha, N, Size);

disp('All copula functions executed successfully. Results saved to output directory.');