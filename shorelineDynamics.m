%% Simulates shoreline change due to gradual erosion, storm impacts, and nourishment
%
% Inputs:
%   pars: parameter struct object. Must contain the following fields:
%       x_init: initial beach width (m)
%       x_nourished: nourished beach width (m)
%       x_crit: beach width that triggers nourishment (m)
%       theta: exponential erosion rate
%       a: historic sea level change rate (m/yr)
%       b: sea level rise acceleration (m/yr^2)
%       r: slope of the active profile
%       H: Bruun Rule correction term (m/yr)
%       lambda: storm frequency (1/yr)
%       k: GEV shape parameter
%       sigma: GEV scale parameter (m)
%       m: GEV location parameter (m)
%       tau_init: time since last nourishment at start of simulation (yr)
%       start_year: year of simulation start
%       end_year: year of simulation end
%       mu: fraction of beach width that is nourished
%       sim_length: total length of simulation (yr)
%   nourishing: boolean, set to 0 to run simulation with no nourishing,
%       default = 1
%
% Outputs:
%   nourish_times: vector containing years in which nourishment occurs
%   x: beach width at each simulation year (m)
%   E: storm induced shoreline change in each year (m/yr)
%   n_storms: number of storms in each year
%%
function [nourish_times,x,E,n_storms] = shorelineDynamics(pars,nourishing)

if nargin < 2
    nourishing = 1;
end
%% Poisson distribution for number of storms
lambda = pars.lambda*ones(pars.sim_length,1);
p0 = poisscdf(0,lambda);
p1 = poisscdf(1,lambda);
p2 = poisscdf(2,lambda);
p3 = poisscdf(3,lambda);
p4 = poisscdf(4,lambda);

%% Sea level rise and bruun rule
t_vec = 0:pars.sim_length;
RSL = pars.a*t_vec+pars.b*(t_vec.^2+2*t_vec*(pars.start_year-1992));
RSL_shift0 = RSL(1:end-1);
RSL_shift1 = RSL(2:end);
gamma_erosion = pars.r*(RSL_shift1-RSL_shift0) - pars.H;

%% Initialize variables
x = zeros(1,pars.sim_length); %beach width
x(1) = pars.x_init;
deltaX = zeros(pars.sim_length-1,1); %change in beach width
E = zeros(pars.sim_length,1); %storm induced shoreline change
n_storms = zeros(pars.sim_length,1); %number of storms
plan = zeros(pars.sim_length,1); %planning for nourishment in next year
nourish_times = []; %years in which nourishment occurs

tau = pars.tau_init;

%plan nourishment if beach width is narrower than critical width
if x(1) <= pars.x_crit
    plan(1) = 1*nourishing;
end

%calculate number of storms and storm induced shoreline change
rand_num = rand;
if rand_num < p0(1)
    n_storms(1) = 0;
elseif rand_num < p1(1)
    n_storms(1) = 1;
elseif rand_num < p2(1)
    n_storms(1) = 2;
elseif rand_num < p3(1)
    n_storms(1) = 3;
elseif rand_num < p4(1)
    n_storms(1) = 4;
else
    n_storms(1) = 5;
end
for i = 1:n_storms(1)
    E(1) = (E(1) + 0.1*gevrnd(pars.k,pars.sigma,pars.m)/i^pars.c);
end

%change in beach width and new beach width
deltaX(1) = -pars.mu*pars.x_nourished*pars.theta*exp(-pars.theta*tau)-gamma_erosion(1)-E(1);
x(2) = max(0,x(1) + deltaX(1));

%% Run simulation
for t = 2:pars.sim_length
   
    %Calculate number of storms in year t
    rand_num = rand;
    if rand_num < p0(t)
        n_storms(t) = 0;
    elseif rand_num < p1(t)
        n_storms(t) = 1;
    elseif rand_num < p2(t)
        n_storms(t) = 2;
    elseif rand_num < p3(t)
        n_storms(t) = 3;
    elseif rand_num < p4(t)
        n_storms(t) = 4;
    else
        n_storms(t) = 5;
    end
    
    %Storm erosion using GEV distribution
    for i = 1:n_storms(t)
        E(t) = (E(t) + 0.1*gevrnd(pars.k,pars.sigma,pars.m)/i);
    end
    
    %Plan for nourishment if beach is narrower than critical width
    nourish = 0;
    if plan(t-1) == 1
        if isempty(nourish_times)
            tau = t;
        end
        nourish = 1;
        nourish_times = cat(1,nourish_times,tau);
        tau = 0;
    elseif x(t) <= pars.x_crit
        plan(t) = 1*nourishing;
    end
    
    %Calculate change in beach width
    if nourish
        deltaX(t) = pars.x_nourished-x(t)-E(t);
    else
        deltaX(t) = -pars.mu*pars.x_nourished*pars.theta*exp(-pars.theta*tau)-gamma_erosion(t)-E(t);
    end
    x(t+1) = max(0,x(t) + deltaX(t));
    tau = tau+1; % Advance time since last nourishment
    
end

end