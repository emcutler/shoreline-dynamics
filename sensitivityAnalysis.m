%% Imports and analyzes output from sensitivity analysis created by main.m
%
% Creates scatter plots of model output vs. each input parameter for each
% location and sea level rise scenario
%
% Calculates the first order sensitivity measure using conditional
% variances for each parameter at each locaiton and sea level rise scenario
%
% Prints expected number of model runs included to calculate each data
% point for each location under each SLR scenario
%
% Imports sensitivity analysis outptut from "sensiOutputs" directory

%%
% Location and SLR scenario variables
locs = {'VB','HI','GI'};
locations = {'Vilano Beach','Hutchinson Island','Gasparilla Island'};
scenarios = {'Low SLR','Intermediate SLR','High SLR'};

% Number of data points included in each slice for conditional variance
sliceWidth = 25;

% Initialize variables for number of runs needed to sufficiently sample
% parameter space and first order sensitivity measures
Runs2Converge = zeros(1,9);
firstOrderS = zeros(9,9);

for loc = 1:3
    for scen = 0:2
        
        % Clear old and import new sensitivity output
        clear sensiPars averageIntervalSampled
        load(strcat('sensiOutputs/',locs{loc},num2str(scen),'sensiOutput.mat'));
        disp(strcat(locations{loc},',',{' '},scenarios{scen+1},':',{' '},'expected',{' '},num2str(mean(included)),{' '},'runs included'));
        
        % Save number of model runs needed to sample parameter space
        numRuns = length(averageIntervalSampled);
        Runs2Converge((scen+1)+(loc-1)*3) = numRuns;
        
        for j = 1:9
            
            % Sort model output by parameter j
            [par,I] = sort(sensiPars(:,j));
            intervalSort = averageIntervalSampled(I);
            
            % Initialize variables for scatter plots
            sliceCenter = zeros(numRuns/sliceWidth,1);
            meanVal = zeros(numRuns/sliceWidth,1);
            
            % Calculate mean output for each slice
            for slice = 1:numRuns/sliceWidth
                meanVal(slice) = mean(intervalSort(sliceWidth*(slice-1)+1:sliceWidth*slice));
                sliceCenter(slice) = par(sliceWidth*slice-floor(sliceWidth/2));
            end
            
            % Calculate first order sensitivity
            firstOrderS((scen+1)+(loc-1)*3,j) = var(meanVal)/var(averageIntervalSampled);
            
            % Create scatter plots
            inputPars = {'$$\hat{x}$$ (meters)','$$x_{crit}$$ (meters)','$$H$$ (meters/yr)','$$\theta$$','r','$$\lambda$$ (yr$$^{-1}$$)','$$k$$','$$\sigma$$ (meters)','$$m$$ (meters)'};
            subplot_labels = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};
            figure((scen+1)+(loc-1)*3)
            subplot(3,3,j);
            plot(par,intervalSort,'.','Color',[.5,.5,.5])
            hold on
            plot(sliceCenter,meanVal,'LineWidth',2,'Color','black')
            ylabel('Mean Interval (yr)','FontSize',18,'FontName','Times New Roman')
            xlabel(inputPars{j},'Interpreter','Latex','FontSize',18,'FontName','Times New Roman')
            set(gca, 'FontSize', 18,'FontName','Times New Roman')
            t = title(subplot_labels{j},'FontWeight','Normal','FontName','Times New Roman');
            set(t, 'horizontalAlignment', 'left')
            set(t, 'units', 'normalized')
            h1 = get(t, 'position');
            set(t, 'position', [0 h1(2) h1(3)])
            sgtitle(char(locations{loc} + ", " + scenarios{scen+1}))
        end
    end
end
