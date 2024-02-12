classdef qualityControl < NeuronSelection.getSpikeTrains
    %QUALITYCONTROL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fs
        fsForce
        QCids
        MUPulses_vec
        MUPulses_vecReduced
        IPTs
        tVec
        Force
        ISI_vec
        COV_ISI_vec
        PNR_vec
        Pauses_vec
    end
    
    methods
        function obj = qualityControl(fs,MUPulses_vec,IPTs,tVec,Force,fsForce)
            obj.MUPulses_vec = MUPulses_vec;
            obj.IPTs = IPTs;
            obj.tVec = tVec;
            obj.fs = fs;
            if nargin >= 5
                obj.Force = Force;
                obj.fsForce = fsForce;
            end
        end
        
        function  execute_quality_control(obj,...
                PausesThres,CoV_ISIThres,PNRThres,ss,minNumPPS,ISI_lb,ISI_ub,withcorr)
            % Based on the quality criteria from Holobar et al.(CoV_IDI<0.3&PNR>30dB),
            % Negro et al.(SIL = 0.9), Martinez et al. (CoV_IDI=0.3 & SIL = 0.9), and
            % Laine et al. (pauses==0 & cov_ISI<thresholds(2) & PNR>thresholds(3))
            % PausesThres(1) is time period, PausesThres(2) is amount of pauses
            % ss --> time interval [start,end]
            
            % ----------------------- variables definition ------------------------
            % pause_threshold=[]; % Or, it can be a vector [timeOfPauseAllowed MaxNumOfPausesAllowed]
            % CoV_IDI_threshold=0.30;
            % PNR_threshold=20;% dB
            % SIL_threshold=0.9;
            % ---------------------------------------------------------------------
            if isempty(PausesThres) %~exist('PausesThres','var')
                maxNumPauses = 3;    % Pauses, COV_IDI, PNR, SIL
            else
                maxNumPauses = PausesThres(2);    % Pauses, COV_IDI, PNR, SIL
            end
            
            % Prelocate variables
            cov_ISI = zeros(1,length(obj.MUPulses_vec));
            pauses = zeros(1,length(obj.MUPulses_vec));
            PNR = zeros(1,length(obj.MUPulses_vec));
            DR = {zeros(1,length(obj.MUPulses_vec))};
            ISI = {zeros(1,length(obj.MUPulses_vec))};
            
            for j = 1:length(obj.MUPulses_vec)
                spikeTrains = zeros(length(obj.tVec),length(obj.MUPulses_vec{j}));
                spikeTrains(obj.MUPulses_vec{j}) = 1;
                limit = length(obj.tVec) + obj.fs;
                window = ceil(0.03*obj.fs);
                i = 0;
                initial = 1;
                
                % Minimum number of pulses
                while initial < limit
                    initial = window*(i)+1;
                    PulsesPerSecond(i+1) = sum(spikeTrains(initial:(window*(i)+obj.fs+1)));
                    i = i + 1;
                end
                
                % Negro decomposition
                if all(PulsesPerSecond < minNumPPS) % uniformity: if there is less than minNumPulses pulses
                    % (around 4 seconds of spikes firing continuosly) the MU is not
                    % considered as it doesn't provide useful information
                    cov_ISI(j) = NaN;
                    pauses(j) = NaN;
                    PNR(j) = NaN;
                else
                    % CoV ISI
                    for i = 2:length(obj.MUPulses_vec{j})
                        DR{j}(i-1) = obj.fs/(obj.MUPulses_vec{j}(i) - obj.MUPulses_vec{j}(i-1));
                        ISI{j}(i-1) = 1./DR{j}(i-1);    % Inter spike interval [s]
                    end
                    cov_ISI(j) = std(ISI{j}((ISI{j} > ISI_lb) & (ISI{j} < ISI_ub)))/...
                        mean(ISI{j}((ISI{j} > ISI_lb) & (ISI{j} < ISI_ub))); 
                    if ~isempty(PausesThres)
                        pauses(j) = sum(ISI{j} > PausesThres(1)); % gives the number for pauses above the threshold
                    else
                        pauses(j) = 0;      % if there is no value this condition does not affect
                    end
                    % PNR
                    t_hat = obj.IPTs(j,1:length(obj.tVec));
                    Discharge = obj.MUPulses_vec{j};
                    noDischarge = setdiff(1:length(t_hat),Discharge);
                    PNR(j) = 10*log10(mean(t_hat(Discharge).^2)/...
                        mean(t_hat(noDischarge).^2));
                end
            end
            obj.ISI_vec = ISI;
            obj.COV_ISI_vec = cov_ISI;
            obj.PNR_vec = PNR;
            obj.Pauses_vec = pauses;
            obj.QCids.Pauses = find(pauses <= maxNumPauses); % only MU with pauses < 0.5
            obj.QCids.COV = find(cov_ISI < CoV_ISIThres);  % only MU with CoV_IDI < 0.3
            obj.QCids.PNR = find(PNR > PNRThres);       % only MU with PNR < 20 dB
            
            if withcorr
                if ~isempty(obj.Force)
                    obj.QC_correlation(PNRThres,CoV_ISIThres)
                else
                    obj.QC_no_correlation(PNRThres,maxNumPauses,CoV_ISIThres)
                end
            else
                obj.QC_no_correlation(PNRThres,maxNumPauses,CoV_ISIThres)
            end
        end
        
        function QC_correlation(obj,PNRThres,CoV_ISIThres)
            % Quality control using correlation and the algorithm used
            % in the paper of Antonio
            Selected = 1:length(obj.MUPulses_vec);
            Selected = Selected(~isnan(obj.PNR_vec));
            
            for MN = 1:length(obj.MUPulses_vec)
                % STEP 1: PNR check
                %                 if obj.PNR_vec(MN) < PNRThres
                %                     Selected = Selected(Selected ~= MN);
                %                 end
                
                % Compute correlation z1
                [~,CST_filtered] = obj.get_spike_trains(Selected,0.4,0);
                Force_res = resample(obj.Force,length(CST_filtered),length(obj.Force));
                force_filtered = normalize(lowpass(Force_res,2,obj.fsForce),'range');
                z1 = corrcoef(CST_filtered,force_filtered);
                
                % STEP 2: CoVISI conditional check
                if obj.COV_ISI_vec(MN) > CoV_ISIThres
                    Selected_2 = Selected(Selected ~= MN);
                    % Compute z2
                    [~,CST_filtered] = obj.get_spike_trains(Selected_2,0.4,0);
                    Force_res = resample(obj.Force,length(CST_filtered),length(obj.Force));
                    force_filtered = normalize(lowpass(Force_res,2,obj.fsForce),'range');
                    z2 = corrcoef(CST_filtered,force_filtered);
                    
                    if z2(1,2) > z1(1,2)
                        Selected = Selected(Selected ~= MN);
                    end
                end
            end
            obj.QCids.Selected = Selected;
            if isempty(obj.QCids.Selected)
                obj.QCids.Selected = 1:length(obj.MUPulses_vec);
            end
        end
        
        function QC_no_correlation(obj,PNRThres,maxNumPauses,CoV_ISIThres)
            %             obj.QCids.Selected = find(obj.Pauses_vec <= maxNumPauses & obj.COV_ISI_vec < CoV_ISIThres & ...
            %                 obj.PNR_vec > PNRThres);
            obj.QCids.Selected = find(obj.Pauses_vec <= maxNumPauses & obj.COV_ISI_vec < CoV_ISIThres ...
                & obj.COV_ISI_vec ~= 0);
            if isempty(obj.QCids.Selected)
                obj.QCids.Selected = 1:length(obj.MUPulses_vec);
            end
        end
    end
end

