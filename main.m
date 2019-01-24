%% Runs shoreline dynamic model to calculate mean nourishment interval with option to conduct sensitivity analysis on parameter values
%
% Inputs:
%   location: case study location, choose 'VB', 'HI', or 'GI'
%   slr_scenario: sea level rise scenario, choose 0 for low, 1 for 
%       intermediate, or 2 for high
%   sensitivityAnalysis: set to 1 to conduct parameter senstivity analysis,
%       default = 0
%   iterations: number of iterations to calculate nourishment interval,
%       default = 1000
%
% Outputs:
%   CI: 90% confidence interval for nourishment interval
%   totalRange: minimum and maximum mean nourishment interval
%   aveInterval: mean nourishment interval
%   sensiPars: random sample of input parameter values
%   averageIntervalSampled: mean nourishment intervals calculated by Monte
%       Carlo simulation
%
% Saves sensitivity analysis output in "sensiOutputs" directory
% Saves nourishment interval output in "outputs" directory

%%
function [CI,totalRange,aveInterval,sensiPars,averageIntervalSampled,excluded,numberNourishments] = main(loc,slr_scenario,sensitivityAnalysis,iterations)

if nargin < 2
    error('Specify location and SLR scenario')
end
if nargin < 3
    sensitivityAnalysis = 0;
end
if nargin < 4
    iterations = 1000; 
end

%Calculate mean nourishment interval descriptive statistics for best guess
%parameter values
pars1 = parameters(loc,slr_scenario,0);
tic
[CI,totalRange,aveInterval,excluded,numberNourishments] = meanInterval(pars1,iterations);
toc

%Conduct sensitivity analysis
if sensitivityAnalysis
    
    averageIntervalSampled = [];
    sensiPars = [];
    included = []; 
    
    for i = 1:500
        % Randomly select input parameters
        pars = parameters(loc,slr_scenario,sensitivityAnalysis);
        
        % Calculate mean nourishment interval
        [~,~,aveInt,removed] = meanInterval(pars,iterations);
        averageIntervalSampled = cat(1,averageIntervalSampled,aveInt);
        included = cat(1,included,iterations-removed);
        
        % Store selected input parameters in sensiPars matrix
        senP = zeros(1,9);
        senP(1) = pars.x_nourished;
        senP(2) = pars.x_crit;
        senP(3) = pars.H;
        senP(4) = pars.theta;
        senP(5) = pars.r;
        senP(6) = pars.lambda;
        senP(7) = pars.k;
        senP(8) = pars.sigma;
        senP(9) = pars.m;
        sensiPars = cat(1,sensiPars,senP);
    end
    VAR = var(averageIntervalSampled);
    MEAN = mean(averageIntervalSampled);    
    varChange = 1;
    meanChange = 1;
    
    while varChange>0.005 && meanChange>0.01
        for i = 1:500
            % Randomly select input parameters
            pars = parameters(loc,slr_scenario,sensitivityAnalysis);
            
            % Calculate mean nourishment interval
            [~,~,aveInt,removed] = meanInterval(pars,iterations);
            averageIntervalSampled = cat(1,averageIntervalSampled,aveInt);
            included = cat(1,included,iterations-removed);
            
            % Store selected input parameters in sensiPars matrix
            senP = zeros(1,9);
            senP(1) = pars.x_nourished;
            senP(2) = pars.x_crit;
            senP(3) = pars.H;
            senP(4) = pars.theta;
            senP(5) = pars.r;
            senP(6) = pars.lambda;
            senP(7) = pars.k;
            senP(8) = pars.sigma;
            senP(9) = pars.m;
            sensiPars = cat(1,sensiPars,senP);
        end
        varChange = abs(var(averageIntervalSampled)-VAR);
        meanChange = abs(mean(averageIntervalSampled)-MEAN);
        VAR = var(averageIntervalSampled);
        MEAN = mean(averageIntervalSampled);
    end
    save(strcat('sensiOutputs/',loc,num2str(slr_scenario),'sensiOutput.mat'),'averageIntervalSampled','sensiPars','included');
end
save(strcat('outputs/',loc,num2str(slr_scenario),'output.mat'),'CI','totalRange','aveInterval','excluded','numberNourishments');

end
%% Function calls shorelineDynamics.m to calculate mean nourishment interval
% 
% Inputs:
%   pars: struct containing parameter values
%   iterations: number of iterations used to calculate mean interval
%
% Outputs:
%   interval90: 90% confidence interval on mean nourishment interval
%   Trange: range of mean nourishment interval
%   Tmean: mean nourishment interval

function [interval90,Trange,Tmean,excluded,numberNourishments] = meanInterval(pars,iterations)

T = {}; % cell array to hold nourishment times for each iteration

% Calculate nourishment times for each iteration
noRenourishments = []; % simulations in which renourishment does not occur
numNourish = zeros(iterations,1); % number of nourishments in simulation

for test_run = 1:iterations
    nourish_times = shorelineDynamics(pars);
    
    % Store number of times renourishment does not occur
    if length(nourish_times) > 1
        T = cat(1,T,nourish_times(2:end));
    else
        noRenourishments = cat(1,noRenourishments,test_run);
    end   
    
    % Store number of nourishments
    numNourish(test_run) = length(nourish_times);
end
numberNourishments = mean(numNourish);
excluded = length(noRenourishments);

if length(T)>1
    Tstar = zeros(length(T),1);
    
    % Calculate mean nourishment interval for each iteration
    for i = 1:length(T)
        Tstar(i) = mean(T{i});
    end
    
    Tmin = min(Tstar); % minimum mean nourishment interval
    Tmax = max(Tstar); % maximum mean nourishment interval
    Trange = [Tmin,Tmax]; % range of mean nourishment interval
    Tmean = mean(Tstar); % mean nourishment interval
    sd_interval = std(Tstar); % standard deviation of nourishment interval
    interval90 = [Tmean-norminv(0.05)*sd_interval/sqrt(length(Tstar)),Tmean+norminv(0.05)*sd_interval/sqrt(length(Tstar))]; % 90% confidence interval of nourishment interval
else
    disp('No renourishment for all runs')
end
end