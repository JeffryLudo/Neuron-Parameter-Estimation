classdef SetInitialParameters < handle
    %SETINITIALPARAMETERS
    %   Detailed explanation goes here
    
    properties (Hidden, Constant, GetAccess=protected)
        eNa = 120  % Sodium reversal potential [120 mV]
        eK = -10   % Potassium reversal potential [-10 mV]
        eL = 0     % Reversal potential leak current [0 mV]
        
        C = 1      % Membrane capacitance [uF/cm^2]
        t2peak = 0.6 % Time to peak [ms]
    end
    properties (Hidden)
        %Initial default parameters
        threshold = 16 % soma firing threshold [mV] 
        gNa = 30   % Maximum gated sodium conductance [mS/cm^2]
        gKf = 4    % Maximum gated fast potassium conductance [mS/cm^2]
        gKs = 17   % Maximum gated slow potassium conductance [mS/cm^2]
        
        % Activation constants
        alphaM = 22  %Na (activation) [ms^-1]
        betaM = 13
        alphaH = 0.5 %Na (Deactivation)
        betaH = 4
        alphaN = 1.5 %Kf (activation)
        betaN = 0.1
        alphaQ = 1.5 %Kf (activation)
        betaQ = 0.025
        
        % Initial parameters (Soma/Dendrite compartment dependant)
        resSoma = 900   % Soma specific resistance [ohm*cm^2]
        dSoma = 77      % Soma diameter [um]
        lSoma = 77      % Soma length [um]
        
        fac = 0 % linear ration between threshold and diameter
        
        % Created here
        gLeakSoma
        cSoma
    end
    properties
        Parameters
    end
    
    methods
        function obj = SetInitialParameters()
            %SETINITIALPARAMETERS
            %   Set all the initial parameters ncessary for ...
            obj.gLeakSoma = 1e5*(pi*obj.dSoma*1e-4*obj.lSoma*1e-4)/obj.resSoma; %(*1e3 to convert from S to mS)
            obj.cSoma = 1e4*pi*obj.dSoma*1e-4*obj.lSoma*1e-4*obj.C; %(*1e8 to adjust diameter from cm back to um *1e-3 to convert microF to miliF = 1e5)
            keySet = {'gNa','gKf','gKs','dSoma','lSoma','resSoma','gLeakSoma','cSoma','threshold',...
                'betaQ','alphaQ','betaN','alphaN','betaH','alphaH','betaM','alphaM','fac'};
            valueSet = [obj.gNa,obj.gKf,obj.gKs,obj.dSoma,obj.lSoma,obj.resSoma,obj.gLeakSoma,obj.cSoma,...
                obj.threshold,obj.betaQ,obj.alphaQ,obj.betaN,obj.alphaN,obj.betaH,obj.alphaH,obj.betaM,obj.alphaM,obj.fac];
            obj.Parameters = containers.Map(keySet,valueSet);
        end
        function set_parameters(obj,parameters,paramsLabel)
            % Overwritting defined parameters
            opfac = false;
            for i = 1:length(paramsLabel)
                obj.Parameters(paramsLabel{i});     % Error if the key is not specified
                
                % Case for having just diameter or length
                if contains(paramsLabel{i},'dSoma')
                    obj.Parameters('dSoma') = parameters(i);
                    obj.Parameters('lSoma') = parameters(i);
                elseif contains(paramsLabel{i},'fac')
                    obj.Parameters(paramsLabel{i}) = parameters(i);
                    opfac = true;
                else
                    obj.Parameters(paramsLabel{i}) = parameters(i);
                end
            end
            
            %Initial default parameters      
            if opfac
                obj.threshold = 0.2 + obj.Parameters('dSoma')*obj.Parameters('fac');
            else 
                obj.threshold = obj.Parameters('threshold'); % soma firing threshold [mV]
            end
            obj.gNa = obj.Parameters('gNa');   % Maximum gated sodium conductance [30 mS/cm^2]
            obj.gKf = obj.Parameters('gKf');    % Maximum gated fast potassium conductance [4 mS/cm^2]
            obj.gKs = obj.Parameters('gKs');   % Maximum gated slow potassium conductance [17 mS/cm^2]
            
            % Activation constants
            obj.alphaM = obj.Parameters('alphaM');  %Na (activation)
            obj.betaM = obj.Parameters('betaM');
            obj.alphaH = obj.Parameters('alphaH'); %Na (Deactivation)
            obj.betaH = obj.Parameters('betaH');
            obj.alphaN = obj.Parameters('alphaN'); %Kf (activation)
            obj.betaN = obj.Parameters('betaN');
            obj.alphaQ = obj.Parameters('alphaQ'); %Kf (activation)
            obj.betaQ = obj.Parameters('betaQ');
            
            % Initial parameters (Soma/Dendrite compartment dependant)
            obj.resSoma = obj.Parameters('resSoma');   % Soma specific resistance [ohm*cm^2]
            obj.dSoma = obj.Parameters('dSoma');     % Soma diameter [um]
            obj.lSoma = obj.Parameters('lSoma');     % Soma length [um]
            
            % Created here
            obj.gLeakSoma = 1e3*(pi*obj.dSoma*1e-4*obj.lSoma*1e-4)/obj.resSoma;
            obj.cSoma = 1e4*pi*obj.dSoma*1e-4*obj.lSoma*1e-4*obj.C;
            
            %% Save new parameter set
%             writematrix([obj.betaQ,obj.threshold,obj.dSoma], strcat('path', 'fileName', '.csv'), 'writeMode', 'append')
        end
    end
end

