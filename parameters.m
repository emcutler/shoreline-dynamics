%% Creates parameter struct object
%
% Inputs:
%   location: case study location, choose 'VB', 'HI', or 'GI'
%   slr_scenario: sea level rise scenario, choose 0 for low, 1 for 
%       intermediate, or 2 for high
%   sensitivityAnalysis: set to 1 to conduct parameter senstivity analysis,
%       default = 0
%
% Output:
%   pars: struct with the following fields:
%       x_init: initial beach width (m)
%       x_nourished: nourished beach width (m)
%       x_crit: beach width that triggers nourishment (m)
%       init_rate: historical shoreline change rate (m/yr)
%       theta: exponential erosion rate
%       a: historic sea level change rate (m/yr)
%       b: sea level rise acceleration (m/yr^2)
%       r: slope of the active profile
%       H: Bruun Rule correction term (m/yr)
%       lambda: storm frequency (1/yr)
%       k: GEV shape parameter
%       sigma: GEV scale parameter (m)
%       m: GEV location parameter (m)
%       c: exponent for subsequent storm impact
%       tau_init: time since last nourishment at start of simulation (yr)
%       start_year: year of simulation start
%       end_year: year of simulation end
%       mu: fraction of beach width that is nourished
%       sim_length: total length of simulation (yr)
%%
function pars = parameters(location,slr_scenario,sensitivityAnalysis)

if nargin < 2
    error('Specify location and SLR scenario')
end
if nargin == 2
    sensitivityAnalysis = 0;
end
if slr_scenario~=0 && slr_scenario~=1 && slr_scenario~=2
    error('SLR scenario incorrectly specified')
end

% Set best guess values
if strcmp(location,'VB')
    pars.x_init = 0;
    pars.x_nourished = 18.288;
    pars.x_crit = 0;
    pars.init_rate = -0.3587;
    pars.theta = 0.1;
    pars.a = 2.4e-3;
    pars.b = 0*(slr_scenario==0) + 0.0271e-3*(slr_scenario==1) + 0.113e-3*(slr_scenario==2);
    pars.r = 50;
    pars.H = pars.init_rate+pars.r*pars.a;
    pars.lambda = 0.41;
    pars.k = 0.277;
    pars.sigma = 4.24;
    pars.m = 1.68;
    pars.tau_init = Inf;
    pars.start_year = 2020;
    pars.end_year = 2070;
elseif strcmp(location,'HI')
    pars.x_init = 0;
    pars.x_nourished = 6.096;
    pars.x_crit = 0;
    pars.init_rate =   -0.0792; % Use -0.0122 for pre-2004
    pars.theta = 0.1;
    pars.a = 2.355e-3;
    pars.b = 0*(slr_scenario==0) + 0.0271e-3*(slr_scenario==1) + 0.113e-3*(slr_scenario==2);
    pars.r = 70.04;
    pars.H = pars.init_rate+pars.r*pars.a;
    pars.lambda = 0.35;
    pars.k = 0.277;
    pars.sigma = 4.24;
    pars.m = 1.68;
    pars.tau_init = 15;
    pars.start_year = 2020;
    pars.end_year = 2070;
elseif strcmp(location,'GI')
    pars.x_init = 42.17;
    pars.x_nourished = 74.02;
    pars.x_crit = 13.21;
    pars.init_rate = -0.2406;
    pars.theta = 0.1;
    pars.a = 2.4e-3;
    pars.b = 0*(slr_scenario==0) + 0.0271e-3*(slr_scenario==1) + 0.113e-3*(slr_scenario==2);
    pars.r = 53.33;
    pars.H = pars.init_rate+pars.r*pars.a;
    pars.lambda = 0.31;
    pars.k = 0.277;
    pars.sigma = 4.24;
    pars.m = 1.68;
    pars.tau_init = 3;
    pars.start_year = 2016;
    pars.end_year = 2056;
else
    error('Location incorrectly specified');
end

pars.c = 10;

% Modify values for sensitivity analysis
if sensitivityAnalysis
    if pars.x_crit == 0
        pars.x_crit = randomPar(pars.x_crit,0,0.4*pars.x_nourished);
    else
        pars.x_crit = randomPar(pars.x_crit);
    end
    pars.x_nourished = randomPar(pars.x_nourished);
    pars.H = randomPar(pars.H);
    pars.theta = randomPar(pars.theta);
    pars.r = randomPar(pars.r);
    pars.lambda = randomPar(pars.lambda);
    pars.k = randomPar(pars.k);
    pars.sigma = randomPar(pars.sigma);
    pars.m = randomPar(pars.m);
end

pars.mu = (pars.x_nourished-pars.x_crit)/pars.x_nourished;
pars.sim_length = pars.end_year-pars.start_year;

end

%% Selects random parameter value for sensitivity analysis
%
% Input:
%   par: selected parameter
%
% Output:
%   par: selected parameter

function par = randomPar(par,parMin,parMax)
if nargin<3 && par~=0
    sensiMin = 0.8;
    sensiMax = 1.2;
    parMin = par*sensiMin;
    parMax = par*sensiMax;
end
par = rand*(parMax-parMin)+parMin;
end
