classdef ResamplingSpikes
    %RESAMPLINGSPIKES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Hidden, GetAccess=private)
        fsExp
        tExp
        fsSim
        spikesExp_vec
    end
    
    methods
        function obj = ResamplingSpikes(spikesExp_vec,fsExp)
            %RESAMPLINGSPIKES Construct an instance of this class
            % Determining simulated sample frequency
            obj.fsExp = fsExp;
            obj.spikesExp_vec = spikesExp_vec;
            obj.tExp = length(spikesExp_vec)/obj.fsExp;
        end
        function [t,spikesResampled] = execute_resampling(obj,spikesSim)
            % Resampling spikes
            obj.fsSim = length(spikesSim)/obj.tExp;
            spikesResampled = obj.resampling(spikesSim);
            % Resampling time
            t = obj.resampling_t();
        end
    end
    methods (Hidden, Access=private)
        function spikesResampled = resampling(obj,spikesSim)
            % Creating matrix of zeros with m = experimental length and n = simulated trains to resample
            spikesResampled = zeros(length(obj.spikesExp_vec),1);
            
            % Getting index of original spike times in new sample frequency
            spikeTimesIndex = round(obj.fsExp*find(spikesSim == 1)/obj.fsSim); % Not very clear????
            spikesResampled(spikeTimesIndex) = 1;
        end
        function t = resampling_t(obj)
            % Generating time vector for new fs
            t = linspace(0,obj.tExp,obj.tExp*obj.fsExp)';
        end
    end
end

