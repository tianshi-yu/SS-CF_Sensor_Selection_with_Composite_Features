%% Main Demo Script for SS-CF
%% Add path
addpath(genpath('.\Data'))
addpath(genpath('.\Func'))

%% Load Data
close all
clear all
rng('default');
load('dataset.mat')

%% Selection settings
q = 1:9; % Select 1-9 sensors
L = Inf*ones(1,9); % No hard constraint on feature number

%% Optimizer settings
% --- PGD solver options --- %
pgdOption = {'gamma',0.7,'verbose',false,'step',1,'maxit',1000,'tol',5e-4};

% --- Optimal socring option --- %
alpha = 0.95; m = 2;  
% Lambda optimisation by line search
lambdaMax = 0.3; nPath = 3; frac = 0.1; 
lambdaMin = frac* lambdaMax;   
lambdaPath = logspace(log10(lambdaMax),log10(lambdaMin),nPath);
lambdaCV = 4;
% Options
tol = 1e-3; maxIte = 50; verbose = true; dispPlot = false;

%% Sensor  selection
% --- Run --- %
disp("Start SSCF Algorithm:")
[ssName,ss,fs,W,J] = runSSCF(X,Y,sensorCandidate,A,L,q,m,alpha,...
                                lambdaPath,lambdaCV,...
                                tol,maxIte,verbose,dispPlot,pgdOption);   

% --- Results --- %
disp("Complete. Sensor selection results:")
for i = 1:length(ssName)
    fprintf('\t - %d sensor(s): %s\n',i,strjoin(ssName{i,1}))
end
disp(" ")
disp("Check results:")
fprintf("\t - Check 'ssName' for the selection results of sensor names;\n")
fprintf("\t - Check 'ss' for the selection results of sensor indices;\n")
fprintf("\t - Check 'fs' for the selection results of feature indices;\n")
fprintf("\t - Check 'W' for the sparse projection matrix result;\n")
fprintf("\t - Check 'J' for the optimal socring objective function value result.\n")
