classdef Optimizer < handle
         
    % Optimizer- given an input vector of parameters raParameters 
    % corresponding to object oProblem with a cost function, optimizer uses
    % gradient descent to find values for raParameters to minimize the cost
    % function
 
    
    %-------------------------STATIC METHODS------------------------------%
    
    methods(Static)
        
        %% Run
        % Invoke this function to solve the optimization
        % problem using gradient descent
        function [raParameters, raCost] = Run(oProblem, raParameters, iSteps, rLearningRate)
                       
            % Create empty cost vector
            raCost = zeros(iSteps,1);
           
            % For each step...
            for iStep = 1:iSteps
                
                % Compute the gradient given the parameters
                [raGradient] = Optimizer.Gradient(oProblem, raParameters); 
                
                % Step downhill (use this order so cost is evaluated with final parameters) 
                raParameters = raParameters - raGradient*rLearningRate;
                
                % Compute the cost 
                [raCost(iStep)] = Cost(oProblem, raParameters);                
                %Supress Sketches to increase time efficiency
                if mod(iStep,100) == 1
                    SketchTruss(oProblem);
                end
                raCost(iStep);
                
                % Comment this block to supress early cutoff of optimization
                % when incremental improvement is sufficiently small
                if iStep > 100
                    if raCost(iStep)/raCost(iStep-100) > 0.9995
                        break; %break for
                    end
                end % end if
                
            end % for
            
        end % function
        
        %% Gradient
        % Compute the gradient numerically
        
        function [raGradient] = Gradient(oProblem, raParameters)
        
            % Use an epsion that is a small fraction of the parameter
            % values
            rEpsilon = max(abs(raParameters))/1e6;

            % Evaluate the cost with no perturbation of the parameters
            rCost0 = Cost(oProblem, raParameters);
            
            % Construct output vector
            raGradient = zeros(size(raParameters));
                
            % For each parameter...
            for iParameter = 1:length(raParameters)
                
                % Copy the parameter vector
                raParameters1 = raParameters;
                
                % Tweak just one component of the vector
                raParameters1(iParameter) = raParameters1(iParameter) + rEpsilon; 
                
                % Compute the cost with the tweaked paraemter
                rCost1 = Cost(oProblem, raParameters1);
                
                % Compute the slope 
                raGradient(iParameter) = (rCost1-rCost0)/rEpsilon;
                 
            end % for          
            
        end % function
          
    end % methods
    
end % classdef