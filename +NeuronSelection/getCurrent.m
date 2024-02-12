classdef getCurrent < handle
    %GETCURRENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fs
        slope = 0.25
        ramp = 1
    end
    
    methods
        function obj = getCurrent(fs,slope,ramp)
            obj.fs = fs;
            if nargin == 2
                obj.slope = slope;
            elseif nargin == 3
                obj.ramp = ramp;
            end
        end
        
        function current = execute_get_current(obj,spikeTrains,plotDR)
            % You can define ramp = 1, slope = 0.25, fs = 2048
            % spikeTrains you get from getExperimentaSpiketrains()
            
            % Getting discharge pattern from MN pool
            [dischargeRatesPool,ISI_real,~] = obj.get_discharge_rate(spikeTrains);
            if plotDR == 1
                obj.plot_discharge_rate(dischargeRatesPool,ISI_real)
            end
            dischargeRatesPool = normalize(dischargeRatesPool,'zscore');
            
            % Estimating mean fire rate from every MN spiketrain reshaped into a single concatenated vector for calculating mean fire rate
            meanFiringRate = mean(1./diff(find(reshape(spikeTrains,[],1)==1)/obj.fs));
            stdFiringRate = std(1./diff(find(reshape(spikeTrains,[],1)==1)/obj.fs));
            
            % Setting discharge pattern from MN pool with mean firing rate from single MN
            MeanFiringRateSignal = meanFiringRate*ones(1,size(spikeTrains,1)) + stdFiringRate*(dischargeRatesPool)/max(dischargeRatesPool);
            
            % Setting input current to MN optimizer (to be implemented). For now, it's *0.9
            current = MeanFiringRateSignal*obj.slope;
            
            % To correct for constant offset in ramp tasks
            if obj.ramp ~= 0
                if abs(min(current)) ~= 0
                    current = current - min(current);
                end
            end
        end
        
        function [dischargeRatesFiltered,ISI_real,meanDischargeRate] = get_discharge_rate(obj,spikeTrains)
            if size(spikeTrains,2) > 1
                cumulativeSpikeTrain = sum(spikeTrains,2);
            else
                cumulativeSpikeTrain = spikeTrains;
            end
            
            spikeTimes = find(cumulativeSpikeTrain >= 1);
            rates = zeros(1,length(spikeTimes)-1);
            
            % Determining firing rate for each spike time
            for i = 1:length(spikeTimes)-1
                rates(i) = obj.fs/(spikeTimes(i+1) - spikeTimes(i)); % times 1000 to be in seconds
            end
            
            % Creating vector dischargesRates with zeroes (to put the firing rate to the corresponding spike time)
            dischargeRates = zeros(1,size(spikeTrains,1));
            
            % Assigning discharge rate to corresponding spike time
            for i = 1:length(spikeTimes)-1
                dischargeRates(spikeTimes(i)) = rates(i);
            end
            
            % If to calculate also discharge rate of single MNs
            if size(spikeTrains,2) > 1
                dischargeRatesFiltered = movmean(dischargeRates,500); % 50 for simulated, 500 experimental (Sartori et al)
            else
                dischargeRatesFiltered = dischargeRates;
            end
            
            % ISI real calculated as 1/dischargeRate times 1000 (to be in ms)
            ISI_real = 1000./dischargeRatesFiltered;
            
            % Getting mean discharge rate
            meanDischargeRate = mean(dischargeRatesFiltered);
        end
        
        function plot_discharge_rate(obj,dischargeRatesFiltered,ISI_real)
            % Ploting
            figure, clf
            subplot(2,1,1)
            plot(linspace(0,length(dischargeRatesFiltered)/obj.fs,length(dischargeRatesFiltered)),...
                dischargeRatesFiltered,'Color',[0,0,1,1]) %/100 to adjust time
            xlim([0,60])
            ylabel('Discharge Rates [s]')
            title('Discharge rate')
            subplot(2,1,2)
            plot(linspace(0,length(dischargeRatesFiltered)/obj.fs,length(dischargeRatesFiltered)),...
                ISI_real,'Color',[0,0,1,1])
            ylim([0,50])
            xlim([0,60])
            ylabel('ISI [ms]')
            title('Measured ISI')
        end
        
        function plot_current(~,tVec,current)
            figure, clf
            plot(tVec,current)
            xlabel('Time [s]')
            ylabel('Current')
            title('Current')
            set(gcf, 'units', 'normalized');
            set(gcf, 'Position', [0, 0, 1, 0.5]);
        end
    end
end

