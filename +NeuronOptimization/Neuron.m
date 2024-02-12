classdef Neuron < handle
    %NEURONPOOL is a class that allow to optimize each motoneuron
    %independently and saves all the important informaation as properties.
    
    properties
        tspan
        window = 0;
        fsExp
        paramsLabel
        parameters0
        lb
        ub
        current
        current_simulation
        t_simulation
        spikesExp_vec
        tVec
        bestParams
        bestFactor
        Parameters
        OptOutput
    end
    properties (Hidden,GetAccess=private)
        errorFnc
        simulation
        resampling
    end
    
    methods
        function obj = Neuron(paramsLabel,lb,ub,current,spikesExp_vec,tVec,window)
            if nargin > 6
                obj.window = window;
            end
            obj.tspan = (tVec(end)-tVec(1));
            obj.fsExp = ceil(length(tVec)/obj.tspan);
            obj.paramsLabel = paramsLabel;
            obj.lb = lb;
            obj.ub = ub;
            obj.current = current;
            obj.spikesExp_vec = spikesExp_vec;
            obj.tVec = tVec;
            % Initialize instance of other classes
            obj.errorFnc = NeuronOptimization.ErrorFunctions.CoincidenceFactor(obj.fsExp);
            obj.simulation = NeuronOptimization.SimulateMN(obj.current,obj.tspan*1000);
            obj.resampling = NeuronOptimization.ResamplingSpikes(obj.spikesExp_vec,obj.fsExp);
            
            obj.Parameters = obj.simulation.Parameters;
            for i = 1:length(paramsLabel)
                obj.parameters0(i) = obj.Parameters(paramsLabel{i});
            end
            
            obj.current_simulation = obj.simulation.I_ext;
            obj.t_simulation = obj.simulation.t;
        end
        function optimize(obj,ErrFnc,OptAlg,max,parameters0)
            % Error function
            switch ErrFnc
                case 'coincidencefactor'
                    fun = @obj.err_coincidence_factor;
                case 'moo'
                    fun = @obj.err_moo;
                case 'rms'
                    fun = @obj.err_rms;
                case 'vanrossum'
                    fun = @obj.err_van_rossum;
                otherwise
                    error('No valid error function selected');
            end
            
            % Choose optimization algorithm
            switch OptAlg
                case 'ga'
                    obj.optimize_ga(fun,max)
                case 'gamultiobj'
                    obj.optimize_gamultiobj(fun,max)
                case 'simulanneal'
                    obj.optimize_simulanneal(fun,max,parameters0)
                case 'particleswarm'
                    obj.optimize_particleswarm(fun,max)
                case 'godlike'
                    obj.optimize_godlike(fun,max,contains(ErrFnc,'moo'))
                otherwise
                    error('No valid optimization algorithm selected');
            end
        end
    end
    %% Optimization Functions
    methods
        function optimize_ga(obj,fun,max)
            opts = optimoptions(@ga,'PlotFcn',@gaplotbestf,'MaxStallGenerations',max,'Display','diagnose','UseParallel',true,...
                'MigrationFraction',0.9);
            [x,fval,~,obj.OptOutput] = ga(fun,length(obj.paramsLabel),[],[],[],[],obj.lb,obj.ub,[],[],opts);
            obj.bestParams = x;
            obj.bestFactor = fval;
        end
        function optimize_gamultiobj(obj,fun,max)
            opts = optimoptions(@gamultiobj,'PlotFcn',@gaplotpareto,'MaxGenerations',max,...
                'Display','diagnose','UseParallel',true);
            [x,fval,~,obj.OptOutput] = gamultiobj(fun,length(obj.paramsLabel),[],[],[],[],obj.lb,obj.ub,opts);
            obj.bestParams = x;
            obj.bestFactor = fval;
        end
        function optimize_simulanneal(obj,fun,max,parameters0)
            opts = optimoptions(@simulannealbnd,'PlotFcn',@saplotbestf,'MaxIterations',max,'Display','diagnose');
            [x,fval,~,obj.OptOutput] = simulannealbnd(fun,parameters0,obj.lb,obj.ub,opts);
            obj.bestParams = x;
            obj.bestFactor = fval;
        end
        function optimize_particleswarm(obj,fun,max)
            opts = optimoptions(@particleswarm,'PlotFcn',@pswplotbestf,'MaxStallIterations',max,'Display','iter','UseParallel',true);
            [x,fval,~,obj.OptOutput] = particleswarm(fun,length(obj.paramsLabel),obj.lb,obj.ub,opts);
            obj.bestParams = x;
            obj.bestFactor = fval;
        end
        function optimize_godlike(obj,fun,max,ismoo)
            if ismoo
                opts = set_options('display','on','maxiters',max,'useparallel','on');
                [x,fval,~,~,~,obj.OptOutput] = GODLIKE(fun,obj.lb,obj.ub,[],opts);
            else
                opts = set_options('display','on','maxiters',max,'useparallel','on'); 
                [x,fval,~,obj.OptOutput] = GODLIKE(fun,obj.lb,obj.ub,[],opts);
            end
            obj.bestParams = x;
            obj.bestFactor = fval;
        end
    end
    %% Error functions
    methods
        function factor = err_coincidence_factor(obj,parameters)
            % Running simulation with default parameter values
            [~,spikesSimResampled] = obj.simulate_resampling(parameters);
            % Change factor position to change factor type
            [~,~,factor] = obj.errorFnc.execute_coincidence(obj.spikesExp_vec,spikesSimResampled,obj.window);
        end
        % Multi Objective Optimization
        function factor = err_moo(obj,parameters)
            % Running simulation with default parameter values
            [~,spikesSimResampled] = obj.simulate_resampling(parameters);
            [factor1,factor2,factor3,factor4] = obj.errorFnc.execute_moo(obj.spikesExp_vec,spikesSimResampled,obj.window);
            factor = [factor1,factor2,factor3,factor4];
        end
        function factor = err_rms(obj,parameters)
            % Running simulation with default parameter values
            [~,spikesSimResampled] = obj.simulate_resampling(parameters);
            factor = obj.errorFnc.execute_rms(obj.spikesExp_vec,spikesSimResampled);
        end
        function factor = err_van_rossum(obj,parameters)
            % Running simulation with default parameter values
            [~,spikesSimResampled] = obj.simulate_resampling(parameters);
            factor = obj.errorFnc.execute_van_rossum(obj.spikesExp_vec,spikesSimResampled);
        end
    end
    %% Plots
    methods
        function plot_spikes(obj,parameters)
            figure, clf, hold on
            plot(obj.tVec,obj.spikesExp_vec,'b','LineWidth',2)
            [~,spikesSimResampled] = obj.simulate_resampling(parameters);
            plot(obj.tVec,spikesSimResampled*0.8,'r','LineWidth',3)
            legend({'Experiment','Model'},'Location','NorthEastOutside');
            set(gcf, 'units', 'normalized');
            set(gcf, 'Position', [0, 0, 1, 0.1]);
            hold off
        end
        function plot_Voltage(obj,parameters)
            figure, clf, hold on
            [~,~,tSim,Vs] = obj.simulate_resampling(parameters);
            plot(tSim,Vs,'b')
            ylabel('Voltage (mV)')
            xlabel('time (s)')
            set(gcf, 'units', 'normalized');
            set(gcf, 'Position', [0, 0, 1, 0.5]);
            hold off
        end
        function plot_iIon(obj,parameters)
            figure, clf, hold on
            [~,~,tSim,~,iIon] = obj.simulate_resampling(parameters);
            plot(tSim,iIon,'b')
            ylabel('Current (mA)')
            xlabel('time (s)')
            set(gcf, 'units', 'normalized');
            set(gcf, 'Position', [0, 0, 1, 0.5]);
            hold off
        end
        function plot_mhnq(obj,parameters)
            figure, clf
            [~,~,tSim,~,~,m,h,n,q] = obj.simulate_resampling(parameters);
            subplot(2,2,1), hold on
            plot(tSim,m,'b')
            ylabel('Value (NaN)')
            xlabel('time (s)')
            title('Constant m')
            subplot(2,2,2), hold on
            plot(tSim,h,'b')
            ylabel('Value (NaN)')
            xlabel('time (s)')
            title('Constant h')
            subplot(2,2,3), hold on
            plot(tSim,n,'b')
            ylabel('Value (NaN)')
            xlabel('time (s)')
            title('Constant n')
            subplot(2,2,4), hold on
            plot(tSim,q,'b')
            ylabel('Value (NaN)')
            xlabel('time (s)')
            title('Constant q')
            set(gcf, 'units', 'normalized');
            set(gcf, 'Position', [0, 0, 1, 0.5]);
            hold off
        end
        function plot_current(obj)
            figure, clf, hold on
            plot(obj.tVec,obj.current,'b')
            ylabel('Current (mA)')
            xlabel('time (s)')
            set(gcf, 'units', 'normalized');
            set(gcf, 'Position', [0, 0, 1, 0.5]);
            hold off
        end
    end
    %% Simulate and Resample neuron
    methods
        function [t,spikesSimResampled,tSim,Vs,iIon,m,h,n,q] = simulate_resampling(obj,parameters)
            % Execute simulation
            [spikesSim,Vs,iIon,m,h,n,q] = obj.simulation.execute_simulation(parameters,obj.paramsLabel);
            tSim = obj.simulation.t;
            % Save Parameters Used
            keySet = {'gNa','gKf','gKs','dSoma','lSoma','resSoma','gLeakSoma','cSoma','threshold',...
                'betaQ','alphaQ','betaN','alphaN','betaH','alphaH','betaM','alphaM','fac'};
            valueSet = [obj.simulation.gNa,obj.simulation.gKf,obj.simulation.gKs,obj.simulation.dSoma,obj.simulation.lSoma,obj.simulation.resSoma,obj.simulation.gLeakSoma,obj.simulation.cSoma,...
                obj.simulation.threshold,obj.simulation.betaQ,obj.simulation.alphaQ,obj.simulation.betaN,obj.simulation.alphaN,obj.simulation.betaH,obj.simulation.alphaH,obj.simulation.betaM,obj.simulation.alphaM,obj.simulation.fac];
            obj.Parameters = containers.Map(keySet,valueSet);
            % Execute resampling
            [t,spikesSimResampled] = obj.resampling.execute_resampling(spikesSim);
        end
    end
end

