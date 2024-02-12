classdef CoincidenceFactor
    %COINCIDENCEFACTOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Hidden, GetAccess=private)
        dt = 0.002; % tolerance = 2 ms, which is about the temporal width of an action potential in cortical neurons (Lynch and Houghton, 2015)
        tau = 0.02 % ms
        fs
    end
    
    methods
        function obj = CoincidenceFactor(fs)
            obj.fs = fs;
        end
        
        function [factor,factorCorrected,factorCorrectedFirstSpike] = execute_coincidence(obj,spikesExp,spikesSim,window)
            if window == 0
                [factor,factorCorrected,factorCorrectedFirstSpike] = obj.factor_no_window(spikesExp,spikesSim);
            else
                [factor,factorCorrected,factorCorrectedFirstSpike] = obj.factor_window(spikesExp,spikesSim,window);
            end
        end
        
        function [factor,frequencyFactor,firstSpikeTime,numberOfSpikes] = execute_moo(obj,spikesExp,spikesSim,window)
            if window == 0
                [factor,~,~,frequencyFactor,firstSpikeTime,numberOfSpikes] = obj.factor_no_window(spikesExp,spikesSim);
            else
                [factor,~,~,frequencyFactor,firstSpikeTime,numberOfSpikes] = obj.factor_window(spikesExp,spikesSim,window);
            end
        end
        
        function rms = execute_rms(~,spikesExp,spikesSim)
            N = length(spikesExp);
            rms = sqrt(1/N * sum((spikesExp - spikesSim).^2));
        end
        
        function VR = execute_van_rossum(obj,spikesExp,spikesSim)
            N = length(spikesExp);
            % Mapp experimental spike train
            mapExp = obj.mapping(spikesExp);
            % Mapp simulated spike train
            mapSim = obj.mapping(spikesSim);
            % Calculate difference with root mean square
            VR = sqrt(1/N * sum((mapExp - mapSim).^2));
        end
    end
    methods (Hidden, Access=private)
        function [factor,factorCorrected,factorCorrectedFirstSpike,frequencyFactor,firstSpikeTime,numberOfSpikes] = ...
                factor_no_window(obj,spikesExp,spikesSim)
            % Getting spike times of experimental and simulated spike trains
            spikeTimesExp = find(spikesExp==1);
            spikeTimesSim = find(spikesSim==1);
            
            % getting number fo experimental and simulated spikes
            numberOfSpikesExp = length(spikeTimesExp);
            numberOfSpikesSim = length(spikeTimesSim);
            
            % Calculating mean firing rate of experimental MNs
            meanFiringRateExp = mean(1./diff(spikeTimesExp/obj.fs));
            meanFiringRateSim = mean(1./diff(spikeTimesSim/obj.fs));
            
            % initializing number of coincidences
            numberOfCoincidences = 0;
            
            % Counting number of coincidences
            for i = 1:numberOfSpikesExp
                % Get all spike times that are above the lower threshold
                coincidences = sum((spikeTimesSim/obj.fs) >= ((spikeTimesExp(i)/obj.fs) - obj.dt) & ...
                    (spikeTimesSim/obj.fs) <= ((spikeTimesExp(i)/obj.fs) + obj.dt));
                % Add the coincidences for each spike in the experimental
                if coincidences >= 1
                    numberOfCoincidences = numberOfCoincidences + coincidences;
                end
            end
            
            % Calculating coincidence factor
            r0 = 2*obj.dt*meanFiringRateExp;
            r1 = 2/(1-r0);
            r2 = numberOfCoincidences - r0*numberOfSpikesExp;
            factor = 1 - r1*r2/(numberOfSpikesExp+numberOfSpikesSim);
            
            % Calcultating frecuency correction factor
            frequencyFactor = 2 * abs(meanFiringRateExp - meanFiringRateSim) / meanFiringRateExp;
            
            % Function to minimize accounting for correction factor
            factorCorrected = frequencyFactor + factor;
            
            % First spike time
            firstSpikeTime = 2 * abs(spikeTimesExp(1) - spikeTimesSim(1)) / spikeTimesExp(1);
            
            % Number of Spikes diffrence
            numberOfSpikes = abs(numberOfSpikesExp - numberOfSpikesSim)/numberOfSpikesExp;
            
            % Function including time to first spike correction
            factorCorrectedFirstSpike = factorCorrected + firstSpikeTime + numberOfSpikes;
        end
        
        function [factor,factorCorrected,factorCorrectedFirstSpike,frequencyFactor,firstSpikeTime,numberOfSpikes] = ...
                factor_window(obj,spikesExp,spikesSim,timeWindow)
            parts = round((length(spikesExp)/obj.fs)/timeWindow); % GEt number of segments for a given time window
            spikeSegment = round(length(spikesExp)/parts);
            
            % Set zeros to the variables
            factor = zeros(1,parts);
            frequencyFactor = zeros(1,parts);
            factorCorrected = zeros(1,parts);
            numberOfSpikes = zeros(1,parts);
            
            % Initialize firstSpike statement
            flag1 = 0;
            flag2 = 0;
            firstSpikes = NaN(2,2);
            
            for j = 0:parts-1
                try
                    spikeTrainsExp = spikesExp(spikeSegment*j+1:spikeSegment*(j+1));
                    spikeTrainsSim = spikesSim(spikeSegment*j+1:spikeSegment*(j+1));
                catch
                    spikeTrainsExp = spikesExp(spikeSegment*j+1:length(spikesExp));
                    spikeTrainsSim = spikesSim(spikeSegment*j+1:length(spikesSim));
                end
                
                % Getting spike times of experimental and simulated spike trains
                spikeTimesExp = find(spikeTrainsExp==1);
                spikeTimesSim = find(spikeTrainsSim==1);
                
                % Getting number fo experimental and simulated spikes
                numberOfSpikesExp = length(spikeTimesExp);
                numberOfSpikesSim = length(spikeTimesSim);
                
                % Calculating mean firing rate of experimental MNs
                meanFiringRateExp = mean(1./diff(spikeTimesExp/obj.fs));
                meanFiringRateSim = mean(1./diff(spikeTimesSim/obj.fs));
                
                % Initializing number of coincidences
                numberOfCoincidences = 0;
                
                % Counting number of coincidences
                for i = 1:numberOfSpikesExp
                    % Get all spike times that are above the lower threshold
                    coincidences = sum((spikeTimesSim/obj.fs) >= ((spikeTimesExp(i)/obj.fs) - obj.dt) & ...
                        (spikeTimesSim/obj.fs) <= ((spikeTimesExp(i)/obj.fs) + obj.dt));
                    % Add the coincidences for each spike in the experimental
                    if coincidences >= 1
                        numberOfCoincidences = numberOfCoincidences + 1;
                    end
                end
                
                % Calculating coincidence factor
                r0 = 2*obj.dt*meanFiringRateExp;
                r1 = 2/(1-r0);
                r2 = numberOfCoincidences - r0*numberOfSpikesExp;
                factor(j+1) = 1 - r1*r2/(numberOfSpikesExp+numberOfSpikesSim);
                
                % Calcultating frecuency correction factor
                frequencyFactor(j+1) = 2 * abs(meanFiringRateExp - meanFiringRateSim) / meanFiringRateExp;
                factorCorrected(j+1) = frequencyFactor(j+1) + factor(j+1);
                
                % First spike time
                if ~isempty(spikeTimesExp) && flag1 == 0
                    firstSpikes(1,1) = spikeTimesExp(1);
                    firstSpikes(1,2) = j;
                    flag1 = 1;
                elseif ~isempty(spikeTimesSim) && flag2 == 0
                    firstSpikes(2,1) = spikeTimesSim(1);
                    firstSpikes(2,2) = j;
                    flag2 = 1;
                end
                
                % Number of Spikes diffrence
                numberOfSpikes(j+1) = abs(numberOfSpikesExp - numberOfSpikesSim)/numberOfSpikesExp;
            end
            % First Spike
            dif = firstSpikes(2,2)-firstSpikes(2,1);
            firstSpikeTime = 2*abs(dif*spikeSegment+(firstSpikes(2,1)-firstSpikes(1,1)))/...
                (firstSpikes(1,2)*spikeSegment+firstSpikes(1,1));
            
            % Deleting NANs and getting mean factor from all segments
            factor = nanmean(factor);
            
            % Deleting NANs and getting mean frequencyFactor from all segments
            frequencyFactor = nanmean(frequencyFactor);
            
            % Deleting NANs and getting mean FactorCorrected from all segments
            factorCorrected = nanmean(factorCorrected);
            
            % Number of Spikes diffrence
            numberOfSpikes = nanmean(numberOfSpikes);
                        
            % Function including time to first spike correction
            factorCorrectedFirstSpike = factorCorrected + firstSpikeTime + numberOfSpikes;
        end
        
        function [map,t] = mapping(obj,spike)
            t0 = find(spike)/obj.fs;  % Esto es U :)
            t  = linspace(0,length(spike)/obj.fs,length(spike));
            map = zeros(size(t));
            vector = zeros(size(t0));
            spikeCounter = 1;
            run = 0;
            
            for i = 1:length(t)
                
                if t(i) >= t0(spikeCounter)
                    run = 1;
                    if spikeCounter < length(t0)
                        spikeCounter = spikeCounter+1;
                    end
                end
                
                if run == 1
                    vector(1:spikeCounter-1) = (2/obj.tau)* exp( -(t(i)-t0(1:spikeCounter-1)) / obj.tau);
                end
                map(i) = sum(vector);
            end
        end
    end
end

