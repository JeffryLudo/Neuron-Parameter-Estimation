classdef MakeFinalPlots
    %MAKEFINALPLOTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        tspan = 5
        fsExp
        paramsLabel
        lb
        ub
        current
        spikesExp
        n_neurons
        parametersFinal
        tVec
    end
    
    methods
        function obj = MakeFinalPlots(paramsLabel,lb,ub,current,spikesExp,tVec,parametersFinal,tspan)
            %MAKEFINALPLOTS Construct an instance of this class
            %   Detailed explanation goes here
            obj.fsExp = ceil(length(spikesExp)/obj.tspan);
            obj.paramsLabel = paramsLabel;
            obj.lb = lb;
            obj.ub = ub;
            obj.current = current;
            obj.spikesExp = spikesExp;
            obj.n_neurons = size(spikesExp,2);
            obj.tVec = tVec;
            obj.parametersFinal = parametersFinal;
            if nargin > 7
                obj.tspan = tspan;
            end
        end
        
        function spikesResampled = plot_spikes(obj)
            figure, clf
            for MNIndex = 1:obj.n_neurons
                subplot(obj.n_neurons+1,1,MNIndex), hold on;
                plot(obj.tVec,obj.spikesExp(:,MNIndex),'b','LineWidth',2)
                simulation = NeuronOptimization.Neuron(obj.paramsLabel,obj.lb,obj.ub,obj.current,obj.spikesExp(:,MNIndex),obj.tVec);
                [~,spikesResampled(:,MNIndex)] = simulation.simulate_resampling(obj.parametersFinal(MNIndex,:));
                plot(obj.tVec,spikesResampled(:,MNIndex)*0.8,'r','LineWidth',3)
            end
            hold off
            xlabel('Time: [s]')
            lg = legend({'Experiment','Model'});
            set(gcf,'units','normalized');
            set(gcf,'Position',[0, 0, 1, 1]);
            set(lg,'Position', [0.5, 0.07, 0.05, 0.05],'units','normalized');
            sgtitle('Spike Trains Comparison')
        end
        function h = plot_histogram(obj)
            h = cell(1,length(obj.paramsLabel));
            figure, clf
            for i = 1:length(obj.paramsLabel)
                subplot(length(obj.paramsLabel),1,i)
                h{1,i} = histfit(obj.parametersFinal(:,i),11,'kernel');
                xlim([min(obj.parametersFinal(:,i)),max(obj.parametersFinal(:,i))])
                xlabel(obj.paramsLabel{i})
            end
            set(gcf, 'units', 'normalized');
            set(gcf, 'Position', [0, 0, 0.6, 1]);
            sgtitle('Parameter Distribution')
        end
    end
end

