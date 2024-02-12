classdef getSpikeTrains < handle
    %GETSPIKETRAINS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Hidden, GetAccess=protected)
        color = {'#A2142F','#D95319','#EDB120','#4DBEEE','#0072BD'};
        linWidth = 0.5;
    end
    
    methods
        function [spikeTrains,CST_filtered,dischargeRates,averageDR,stdDR] = get_spike_trains(obj,selectedPulses,hannWindow,plotFlag,plotF)
            if ~isempty(selectedPulses)
                obj.MUPulses_vecReduced = obj.MUPulses_vec(selectedPulses);
            else
                obj.MUPulses_vecReduced = obj.MUPulses_vec;
            end
            % Initializing variables
            spikeTrains = zeros(length(obj.tVec),length(obj.MUPulses_vecReduced));
            dischargeRates = NaN(length(obj.tVec),length(obj.MUPulses_vecReduced));
            for muCount = 1:length(obj.MUPulses_vecReduced)
                spikeTrains(obj.MUPulses_vecReduced{muCount},muCount) = 1;
                dischargeRates(obj.MUPulses_vecReduced{muCount}(2:end),muCount) = obj.fs./diff(obj.MUPulses_vecReduced{muCount});
                % values below 3 hz may be due to big pauses (non physiological)
                dischargeRates(dischargeRates < 3 | dischargeRates > 100) = NaN;
            end
            averageDR = nanmean(dischargeRates);
            stdDR = nanstd(dischargeRates);
            % Getting cumulative spike train
            cumulativeSpikeTrain = sum(spikeTrains,2);
            % Making the last point zero (otherwise, the filter may go to infinite in the extremes)
            cumulativeSpikeTrain(length(cumulativeSpikeTrain)) = 0;
            % Low pass filtering with a Hann moving average window of 400ms
            w = hann(round(hannWindow*obj.fs));
            CST_filtered = normalize(filtfilt(w,1,cumulativeSpikeTrain),'range');
            % Plot
            if plotFlag
                if ~isempty(selectedPulses)
                    noselectedPulses = setdiff(1:length(obj.MUPulses_vec),selectedPulses);
                    obj.MUPulses_vecReduced = obj.MUPulses_vec;
                    for idx = noselectedPulses
                        obj.MUPulses_vecReduced{idx} = [];
                    end
                    dischargeRates2 = NaN(length(obj.tVec),length(obj.MUPulses_vec));
                    dischargeRates2(:,selectedPulses) = dischargeRates;
                    averageDR2 = zeros(1,length(obj.MUPulses_vec));
                    averageDR2(selectedPulses) = averageDR;
                    obj.plot_spikes(dischargeRates2,averageDR2,plotF)
                    obj.MUPulses_vecReduced = obj.MUPulses_vec(selectedPulses);
                else
                    obj.plot_spikes(dischargeRates,averageDR,plotF)
                end
            end
        end
        
        function plot_spikes(obj,dischargeRates,averageDR,plotF)
            colSeq = 1; % initializing color sequence
            
            figure, clf, hold on
            for muCount = 1:length(obj.MUPulses_vecReduced)
                % To plot spikes neatly, we need to give two values per instance
                % (sample), i.e. we need to go from 0 to 1 in a single sample,
                % hence the elements (samples) of t2 are repeated and the spike
                % locations are updated
                t2 = repelem(obj.tVec,2); % duplicating the samples of t vector
                spikePos = (obj.MUPulses_vecReduced{muCount})*2; % updating spike locations
                % Plot variable is full of NaNs to leave blank the samples without spikes
                plotVar = NaN(size(t2)); % creating plot variable
                if ~isempty(spikePos)
                    plotVar(spikePos) = muCount-0.4;
                    plotVar(spikePos+1) = muCount+0.4;
                end
                ax1 = subplot(1,2,1); hold on
                plot(t2,plotVar,'-','Color',obj.color{colSeq},'LineWidth',obj.linWidth);
                title('MU spike trains')
                %         plot(tVec(1:end-1),0.8*IPTs1(muCount,length(tVec(1:end-1)))/max(IPTs1(muCount,:))+ muCount-0.4,'Color', color{colSeq})
                ax2 = subplot(1,2,2); hold on
                plot(obj.tVec,dischargeRates(:,muCount)/30+muCount-0.5,...
                    '.','Color',obj.color{colSeq},'MarkerSize',10); % normalized to 30 pps
                title('MU discharge rates')
                if colSeq < 5
                    colSeq = colSeq + 1 ;
                else
                    colSeq = 1;
                end
            end
            % Add labels and more
            ax = [ax1,ax2];
            subplot(1,2,1)
            yyaxis left
            xlim([-1, max(obj.tVec)+1]);
            ylim([0 length(obj.MUPulses_vecReduced)+0.5]);
            yticks(0:length(obj.MUPulses_vecReduced))
            xlabel('Time (s)')
            ylabel('Motor Neuron (#)')
            yyaxis right
            ylim([0 length(obj.MUPulses_vecReduced)+0.5]);
            yticks(0:length(obj.MUPulses_vecReduced))
            B = num2str(obj.COV_ISI_vec(:),'%.2f');
            I = strfind(B(:)','NaN'); % Find NaNs
            B([I I+1 I+2]) = ' '; % Replace NaN with spaces
            yticklabels({'',B});
            ylabel(ax(1),'Coefficient of Variation (CoV)')
            ax(1).YAxis(1).Color = 'k';
            ax(1).YAxis(2).Color = 'k';
            subplot(1,2,2)
            yyaxis left
            xlim([-1, max(obj.tVec)+1]);
            ylim([0 length(obj.MUPulses_vecReduced)+0.5]);
            yticks(0:length(obj.MUPulses_vecReduced))
            yticklabels({'',num2str(averageDR(:),'%.1f')});
            xlabel('Time (s)')
            ylabel('Mean discharge rate (pps)')
            yyaxis right
            ylim([0 length(obj.MUPulses_vecReduced)+0.5]);
            yticks(0:length(obj.MUPulses_vecReduced))
            B2 = num2str(obj.PNR_vec(:),'%.1f');
            I2 = strfind(B(:)','NaN'); % Find NaNs
            B2([I2 I2+1 I2+2]) = ' '; % Replace NaN with spaces
            yticklabels({'',B2});
            ylabel(ax(2),'Pulse to Noise Ratio (PNR)')
            ax(2).YAxis(1).Color = 'k';
            ax(2).YAxis(2).Color = 'k';
            if ~isempty(obj.Force) && plotF
                for i = 1:2
                    %                     yyaxis(ax(i),'right')
                    %                     ax(i).YAxis(2).Color = 'k';
                    Force_res = resample(obj.Force,length(obj.tVec),length(obj.Force));
                    force_filtered = normalize(lowpass(Force_res,2,obj.fsForce),'range');
                    force_filtered = force_filtered*length(obj.MUPulses_vecReduced);
                    hold on
                    plot(ax(i),obj.tVec,force_filtered,'Color',[0.7,0.7,0.7,0.5]);
                    %                     ylabel(ax(i),'Force (N)')
                end
            end
            set(gcf, 'units', 'normalized');
            set(gcf, 'Position', [0, 0, 1, 1]);
            linkaxes(ax)
        end
    end
end

