% 6 trusses are created and tested for an identical load case (a central
% load). The relative costs of each are compared

% Programwide Constants
clear all; clc; close all;
iSteps = 20000;
rLearningRate = 5e-6;

%% Truss 1 - Pratt
% Truss 1 - Trapezoidal Pratt
nodes = [...

    0     14;
    32.7  14;
    65.3  14;
    98    14;
    130.7 14;
    163.3 14;
    196   14;
    
    32.7  0;
    65.3  0;
    98    0;
    130.7 0;
    163.3 0;
    
    0    7;
    196  7]; 

nodeLock = ...
    [
    1 0;
    1 0;
    1 0;
    1 0;
    1 0;
    1 0;
    1 0;

    1 0;
    1 0;
    1 0;
    1 0;
    1 0;
    
    0 1;
    0 1];

XC1 = 0.3125;
members = [1  2 XC1; 
           2  3 XC1;
           3  4 XC1;
           4  5 XC1;
           5  6 XC1;
           6  7 XC1;
           
           1 13 XC1;
           2  8 XC1;
           3  9 XC1;
           4 10 XC1;
           5 11 XC1;
           6 12 XC1;
           7 14 XC1;
           
           13 8 XC1;
           8  9 XC1;
           9 10 XC1;
           10 11 XC1;
           11 12 XC1;
           12 14 XC1;
           
           1   8 XC1;
           2   9 XC1;
           3  10 XC1;
           10  5 XC1;
           11  6 XC1;
           12  7 XC1];
           

loads = [0 0 ...
         0 0 ...
         0 -300 ...
         0 -600 ...
         0 -300 ...
         0 0 ...
         0 0 ...
         ...
         0 0 ...
         0 0 ...
         0 0 ...
         0 0 ...
         0 0 ...
         0];

truss1 = Truss(nodes,nodeLock,members(:,1:2),members(:,3),loads);
truss1.SketchDeformedTruss();      

raParameters = truss1.GetParameters();
[raParameters, raCost] = Optimizer.Run(truss1, raParameters,iSteps, rLearningRate); 
truss1raCost = raCost;

%% Truss 2 - Howe

% Truss 2 - Trapezoidal Howe
nodes = [...

    0     14;
    32.7  14;
    65.3  14;
    98    14;
    130.7 14;
    163.3 14;
    196   14;
    
    32.7  0;
    65.3  0;
    98    0;
    130.7 0;
    163.3 0;
    
    0    7;
    196  7]; 

nodeLock = ...
    [
    1 0;
    1 0;
    1 0;
    1 0;
    1 0;
    1 0;
    1 0;

    1 0;
    1 0;
    1 0;
    1 0;
    1 0;
    
    0 1;
    0 1];

XC1 = 0.3125;
members = [1  2 XC1; 
           2  3 XC1;
           3  4 XC1;
           4  5 XC1;
           5  6 XC1;
           6  7 XC1;
           
           1 13 XC1;
           2  8 XC1;
           3  9 XC1;
           4 10 XC1;
           5 11 XC1;
           6 12 XC1;
           7 14 XC1;
           
           13 8 XC1;
           8  9 XC1;
           9 10 XC1;
           10 11 XC1;
           11 12 XC1;
           12 14 XC1;
           
           13  2 XC1;
           8   3 XC1;
           9   4 XC1;
           4  11 XC1;
           5  12 XC1;
           6  14 XC1];
           

loads = [0 0 ...
         0 0 ...
         0 -300 ...
         0 -600 ...
         0 -300 ...
         0 0 ...
         0 0 ...
         ...
         0 0 ...
         0 0 ...
         0 0 ...
         0 0 ...
         0 0 ...
         0];

truss2 = Truss(nodes,nodeLock,members(:,1:2),members(:,3),loads);

truss2.SketchDeformedTruss();      

raParameters = truss2.GetParameters();
[raParameters, raCost] = Optimizer.Run(truss2, raParameters,iSteps, rLearningRate); 
truss2raCost = raCost;


%% Truss 3 - Warren Truss
% Truss 3 - Trapezoidal Warren Up
rLearningRate = 3e-6;
nodes = [...

    0     14;
    32.7  14;
    65.3  14;
    98    14;
    130.7 14;
    163.3 14;
    196   14;
    
    32.7  0;
    65.3  0;
    98    0;
    130.7 0;
    163.3 0;
    
    0    7;
    196  7]; 

nodeLock = ...
    [
    1 0;
    1 0;
    1 0;
    1 0;
    1 0;
    1 0;
    1 0;

    1 0;
    1 0;
    1 0;
    1 0;
    1 0;
    
    0 1;
    0 1];

XC1 = 0.3125;
members = [1  2 XC1; 
           2  3 XC1;
           3  4 XC1;
           4  5 XC1;
           5  6 XC1;
           6  7 XC1;
           
           1 13 XC1;
           2  8 XC1;
           3  9 XC1;
           4 10 XC1;
           5 11 XC1;
           6 12 XC1;
           7 14 XC1;
           
           13 8 XC1;
           8  9 XC1;
           9 10 XC1;
           10 11 XC1;
           11 12 XC1;
           12 14 XC1;
           
           1   8 XC1;
           8   3 XC1;
           3  10 XC1;
           10  5 XC1;
           5  12 XC1;
           12  7 XC1];
           

loads = [0 0 ...
         0 0 ...
         0 -300 ...
         0 -600 ...
         0 -300 ...
         0 0 ...
         0 0 ...
         ...
         0 0 ...
         0 0 ...
         0 0 ...
         0 0 ...
         0 0 ...
         0];

truss3 = Truss(nodes,nodeLock,members(:,1:2),members(:,3),loads);
truss3.SketchDeformedTruss();      

raParameters = truss3.GetParameters();
[raParameters, raCost] = Optimizer.Run(truss3, raParameters,iSteps, rLearningRate); 
truss3raCost = raCost;

%% Truss 4 - Warren Down

% Truss 4 - Trapezoidal Warren Down
nodes = [...

    0     14;
    32.7  14;
    65.3  14;
    98    14;
    130.7 14;
    163.3 14;
    196   14;
    
    32.7  0;
    65.3  0;
    98    0;
    130.7 0;
    163.3 0;
    
    0    7;
    196  7]; 

nodeLock = ...
    [
    1 0;
    1 0;
    1 0;
    1 0;
    1 0;
    1 0;
    1 0;

    1 0;
    1 0;
    1 0;
    1 0;
    1 0;
    
    0 1;
    0 1];

XC1 = 0.3125;
members = [1  2 XC1; 
           2  3 XC1;
           3  4 XC1;
           4  5 XC1;
           5  6 XC1;
           6  7 XC1;
           
           1 13 XC1;
           2  8 XC1;
           3  9 XC1;
           4 10 XC1;
           5 11 XC1;
           6 12 XC1;
           7 14 XC1;
           
           13 8 XC1;
           8  9 XC1;
           9 10 XC1;
           10 11 XC1;
           11 12 XC1;
           12 14 XC1;
           
           13  2 XC1;
           2   9 XC1;
           9   4 XC1;
           4  11 XC1;
           11  6 XC1;
           6  14 XC1];
           

loads = [0 0 ...
         0 0 ...
         0 -300 ...
         0 -600 ...
         0 -300 ...
         0 0 ...
         0 0 ...
         ...
         0 0 ...
         0 0 ...
         0 0 ...
         0 0 ...
         0 0 ...
         0];

truss4 = Truss(nodes,nodeLock,members(:,1:2),members(:,3),loads);
truss4.SketchDeformedTruss();      

raParameters = truss4.GetParameters();
[raParameters, raCost] = Optimizer.Run(truss4, raParameters,iSteps, rLearningRate); 
truss4raCost = raCost;


%% Truss 5 - Offset
iSteps = 20000;
rLearningRate = 5e-6
% Truss 5 - Offset

nodes = [...

    0     14;
    32.7  14;
    65.3  14;
    98    14;
    130.7 14;
    163.3 14;
    196   14;
    
    16.3  0;
    49.0  0;
    81.6  0;
    114.3 0;
    147.0 0;
    179.6 0;
    
    0    7;
    196  7]; 

nodeLock = ...
    [
    1 0;
    1 0;
    1 0;
    1 0;
    1 0;
    1 0;
    1 0;

    1 0;
    1 0;
    1 0;
    1 0;
    1 0;
    1 0;
    
    0 1;
    0 1];

XC1 = 0.3125;
members = [1  2 XC1; 
           2  3 XC1;
           3  4 XC1;
           4  5 XC1;
           5  6 XC1;
           6  7 XC1;
           
           1 14 XC1;
           7 15 XC1;
           
           8   9 XC1;
           9  10 XC1;
           10 11 XC1;
           11 12 XC1;
           12 13 XC1;
           
           14  8 XC1;
           1   8 XC1;
           8   2 XC1;
           2   9 XC1;
           9   3 XC1;
           3  10 XC1;
           10  4 XC1;
           4  11 XC1;
           11  5 XC1;
           5  12 XC1;
           12  6 XC1;
           6  13 XC1;
           13  7 XC1;
           13 15 XC1;];
           

loads = [0 0 ...
         0 0 ...
         0 -300 ...
         0 -600 ...
         0 -300 ...
         0 0 ...
         0 0 ...
         ...
         0 0 ...
         0 0 ...
         0 0 ...
         0 0 ...
         0 0 ...
         0 0 ...
         0];

truss5 = Truss(nodes,nodeLock,members(:,1:2),members(:,3),loads);
truss5.SketchTruss();
truss5.SketchDeformedTruss();      

raParameters = truss5.GetParameters();
[raParameters, raCost] = Optimizer.Run(truss5, raParameters,iSteps, rLearningRate); 
truss5raCost = raCost;



%% Truss 6 - Bridget
nodes = [...

    0     14;
    32.7  14;
    65.3  14;
    98    14;
    130.7 14;
    163.3 14;
    196   14;
    
    32.7  0;
    65.3  0;
    98    0;
    130.7 0;
    163.3 0;
    
    16.3  7;
    49.0  7;
    81.6  7;
    114.3 7;
    147.0 7;
    179.6 7;
     
    0    7;
    196  7]; 

nodeLock = ...
    [
    1 0;
    1 0;
    1 0;
    1 0;
    1 0;
    1 0;
    1 0;

    1 0;
    1 0;
    1 0;
    1 0;
    1 0;
    
    1 1;
    1 1;
    1 1;
    1 1;
    1 1;
    1 1;
    
    0 1;
    0 1];

XC1 = 0.3125;
members = [1  2 XC1; 
           2  3 XC1;
           3  4 XC1;
           4  5 XC1;
           5  6 XC1;
           6  7 XC1;
           
           1 19 XC1;
%            2  8 XC1;
%            3  9 XC1;
%            4 10 XC1;
%            5 11 XC1;
%            6 12 XC1;
           7 20 XC1;
           
           19 8 XC1;
           8  9 XC1;
           9 10 XC1;
           10 11 XC1;
           11 12 XC1;
           12 20 XC1;
           
           1  13 XC1;
           19 13 XC1;
           8  13 XC1;
           2  13 XC1;
           2  14 XC1;
           3  14 XC1;
           8  14 XC1;
           9  14 XC1;
           3  15 XC1;
           4  15 XC1;
           9  15 XC1;
           10 15 XC1;
           4  16 XC1;
           5  16 XC1;
           10 16 XC1;
           11 16 XC1;
           5  17 XC1;
           6  17 XC1;
           11 17 XC1;
           12 17 XC1;
           6  18 XC1;
           7  18 XC1;
           12 18 XC1;
           20 18 XC1;
           
           13 14 XC1;
           14 15 XC1;
           15 16 XC1;
           16 17 XC1;
           17 18 XC1];
           

loads = [0 0 ...
         0 0 ...
         0 -300 ...
         0 -600 ...
         0 -300 ...
         0 0 ...
         0 0 ...
         ...
         0 0 ...
         0 0 ...
         0 0 ...
         0 0 ...
         0 0 ...
         ...
         0 0 ...
         0 0 ...
         0 0 ...
         0 0 ...
         0 0 ...
         0 0 ...
         0];

truss6 = Truss(nodes,nodeLock,members(:,1:2),members(:,3),loads);
truss6.SketchTruss();
truss6.SketchDeformedTruss();      

raParameters = truss6.GetParameters();
[raParameters, raCost] = Optimizer.Run(truss6, raParameters,iSteps, rLearningRate); 
truss6raCost = raCost;





%% Plots
% Yep



% Plot cost function over time
figure, semilogy(truss1raCost);
hold on;
semilogy(truss2raCost);
semilogy(truss3raCost);
semilogy(truss4raCost);
semilogy(truss5raCost);
semilogy(truss6raCost);
title('Cost Trajectory');
xlabel('Step');
ylabel('Cost');
legend('Pratt', 'Howe', 'Warren Up' ,'Warren Down', 'Offset', 'Bridget');
grid on;



