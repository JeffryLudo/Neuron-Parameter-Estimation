classdef SimulateMN < NeuronOptimization.SetInitialParameters
    %SIMULATEMN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Hidden, GetAccess=private)
        dt = 0.05   % Time step for forward euler method (ms) <- It is n = time[ms]*fs/1000. So dt = time[s]/fs *(1/1000) for ms
    end
    
    properties
        loop        % no. of iterations of euler method
        t           % Time vector [ms]
        I_ext
    end
    
    methods
        function obj = SimulateMN(I,tspan,dt)
            %SIMULATEMN Construct an instance of this class
            if nargin == 3
                obj.dt = dt;
                obj.loop = ceil(tspan/obj.dt);
                obj.t = (1:obj.loop)*obj.dt;
            else
                obj.loop = ceil(tspan/obj.dt);
                obj.t = (1:obj.loop)*obj.dt;
            end
            
            % Adjusting current sample size
            obj.I_ext = obj.adjust_size(I);
        end
        function [spikesSim,Vs,iIon,m,h,n,q] = execute_simulation(obj,parameters,paramsLabel)
            % Initial parameters (General and Conductance dependent)
            obj.set_parameters(parameters,paramsLabel);
            % Euler Method
            [Vs,iIon,m,h,n,q] = obj.euler_method;
            % Find Spikes
            spikesSim = obj.find_spikes(Vs);
        end
    end
    
    methods (Hidden, Access=private)
        function I_ext = adjust_size(obj,I)
            % Adjusting current sample size
            if length(I) > 1
                if  obj.loop ~= length(I)
                    I_ext_1 = resample(I(1:ceil(length(I)/2)),...
                        ceil(length(I(1:ceil(length(I)/2)))*obj.loop/length(I)),length(I(1:ceil(length(I)/2))));  
                    I_ext_2 = resample(I((ceil(length(I)/2)+1):end),...
                        ceil(length(I((ceil(length(I)/2)+1):end))*obj.loop/length(I))-1,length(I((ceil(length(I)/2)+1):end))); 
                    I_ext = [I_ext_1 I_ext_2]; % Resample current length to loop
                else
                    I_ext = I;
                end
            else
                I_ext = I*ones(obj.loop,1);    % time points for injected current  (*1e-6 to convert from nA to mA)
            end
        end
        function [Vs,iIon,m,h,n,q] = euler_method(obj)
            % Initializing variable vectors
            iIon = zeros(obj.loop,1);
            Vs = zeros(obj.loop,1);
            m = zeros(obj.loop,1);
            h = zeros(obj.loop,1);
            n = zeros(obj.loop,1);
            q = zeros(obj.loop,1);
            
            % Set initial values for the variables
            Vs(1) = 0;
            m(1) = 0*ones(1,1)+0.01*rand(1,1);
            h(1) = 1*ones(1,1)-0.01*rand(1,1);
            n(1) = 0*ones(1,1)+0.01*rand(1,1);
            q(1) = 0*ones(1,1)+0.01*rand(1,1);
            
            % State
            state = 0;
            t0 = 1;
            
            for i = 1:obj.loop-1
                iIon(i) = obj.gNa*(m(i)^3)*h(i)*(obj.eNa-(Vs(i))) + ...
                    obj.gKf*(n(i)^4)*(obj.eK-(Vs(i))) +  ...
                    obj.gKs*(q(i)^2)*(obj.eK-(Vs(i))); % 
                Vs(i+1) = Vs(i) + obj.dt*(obj.gLeakSoma*(obj.eL-(Vs(i))) + iIon(i) + obj.I_ext(i))/obj.cSoma;
                [m(i+1),h(i+1)] = obj.get_Na_conductance(m(t0),h(t0),obj.betaM,obj.alphaM,...
                    obj.betaH,obj.alphaH,state,obj.t(i),obj.t(t0));
                n(i+1) = obj.get_K_conductance(n(t0),obj.betaN,obj.alphaN,state,obj.t(i),obj.t(t0));
                q(i+1) = obj.get_K_conductance(q(t0),obj.betaQ,obj.alphaQ,state,obj.t(i),obj.t(t0));
                
                % Setting the onset time (t0) of the pulse (state,1) when Vs reaches threshold
                if Vs(i+1) >= obj.threshold
                    if state(1) == 0 % To ensure that the time t0 is only stored the first time that Vs > threshold
                        state(1) = 1; % Activate pulse
                        t0 = i; % Saving the index of time (t0) at which pulse starts
                    end
                end
                
                if state(1) == 1
                    state(2) = 0;
                    if obj.t(i)-obj.t(t0) > obj.t2peak % Turning of pulse (state,1) and switching activation variables (state,2) when 0.6s period is achieved
                        state(1) = 0;
                        state(2) = 1;
                        t0 = i; % Saving the index of time (t0) at which pulse ends
                    end
                end
            end
        end
        function spikesSim = find_spikes(obj,Vs)
            spikesSim = zeros(obj.loop,1);
            
            %Storing spike times
            [~,locs] = findpeaks(Vs,'MinPeakHeight',obj.threshold);
            spikesSim(locs) = 1;
        end
    end
    methods (Hidden, Access=private)
        % alpha and beta functions for the gating variables
        function [m,h] = get_Na_conductance(obj,m0,h0,betaM,alphaM,betaH,alphaH,state,t,t0)
            if state(1) == 1
                if state(2) == 1
                    m = obj.get_value_on(betaM,t,t0,m0);
                    h = obj.get_value_off(alphaH,t,t0,h0);
                else
                    m = obj.get_value_off(alphaM,t,t0,m0);
                    h = obj.get_value_on(betaH,t,t0,h0); 
                end
            else
                m = obj.get_value_on(betaM,t,t0,m0);
                h = obj.get_value_off(alphaH,t,t0,h0);
            end
        end
        
        function [kConstant] = get_K_conductance(obj,kConstant0,beta,alpha,state,t,t0)
            if state(1) == 1
                if state(2) == 1
                    kConstant = obj.get_value_on(beta,t,t0,kConstant0);
                else
                    kConstant = obj.get_value_off(alpha,t,t0,kConstant0);
                end
            else
                kConstant = obj.get_value_on(beta,t,t0,kConstant0);
            end
            
        end
        
        function value = get_value_on(~,beta,t,t0,v0)
            value = v0*exp(-beta*(t-t0));
        end
        
        function value = get_value_off(~,alpha,t,t0,v0)
            value = 1+(v0-1)*exp(-alpha*(t-t0));
        end
    end
end

