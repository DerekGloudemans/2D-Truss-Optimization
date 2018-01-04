classdef Truss < matlab.mixin.Copyable %handle
    % Truss- a 2D simple truss on which stiffness matrix method analysis and
    % related functions can be performed
    
    % m- number of raNodes
    % n- number of raMembers
    
    % Variable Notation
    % r.. - real noninteger
    % n.. - integer
    % b.. - boolean
    % a.. - array
    % i,j,k- iterator variable

    
    %----------------------------PROPERTIES-------------------------------%
    
    properties
        %   raNodes- an mx2 matrix of xcoords and ycoords for each node
        %   raMembers- an nx2 matrix of start and end node for each member
        %   raMemberXC- an nx1 matrix of XC area for each member
        %   raNodeLock- an mx2 matrix which indicates whether the position of the
        %             node may be altered along each degree of freedom (1 = yes, 0 = no)
        %   raLoads- a 1x(2m-4) matrix of raLoads applied to each degree of freedom
        %          degrees of freedom for node 1 are 1 and 2, 3 and 4 for node
        %          2, etc.
        raNodes;
        raNodeLock;
        raMemberXC;
        raMembers;
        raLoads;
        
    end
    
    properties(Constant)
        % rAreaScaling- scales XC areas so that scale is similar to that of
            % XY node coordinates (for better gradient descent performance)
        % nNumReactions- set to 3 for pin-roller, or 4 for pin-pin truss
        rAreaScaling = 100;
        nNumReactions = 3;
        
            % (Properties for 1018 Steel - For analysis of many materials, 
            % input these variables in the constructor)
        % E- Elastic Modulus, psi
        % rDensity- lb/in^3
        % rYieldStrength- psi
        E = 2.9 * 10^7; 
        rDensity = 0.284; 
        rYieldStrength = 3.6*10^4; 
        
        % Weighting Coefficients (set to ASCE bridge competition specs)
        % a- controls the extent to which min XC area affects cost function
        % b- controls the extent to which min buckling FoS area affects cost function
        % y- controls the extent to which min yield FoS affects cost function
        % l- controls the extent to which max length affects cost function
        % d- controls the extent to which depth of truss affects cost function
        a = 10;
        b = 100000;
        y = 100000;
        l = 100000;
        d = 0;
        
    end
    
    
    %------------------------------METHODS--------------------------------%          
    
    methods
        
        % Returns a truss object with user-input raNodes, raNodeLock,
            % raMembers, raMemberXCs and raLoads
        function [obj] = Truss(n,nl,m,mxc,lo)
            obj.raNodes = n;
            obj.raNodeLock = nl;
            obj.raMemberXC = mxc;
            obj.raMembers = m;
            obj.raLoads = lo;
           
        end
        
        % Returns a nx1 matrix of member lengths
        function [ raMemLengths ] = MemberLengths (obj)
            raMemLengths = ((obj.raNodes(obj.raMembers(:,1),1)-obj.raNodes(obj.raMembers(:,2),1)).^2 + (obj.raNodes(obj.raMembers(:,1),2)-obj.raNodes(obj.raMembers(:,2),2)).^2).^.5;
        end
        
        % Returns a nx1 matrix of member weights
        function [ raMemWeights ] = MemberWeights(obj)
   
            raMemLengths = MemberLengths(obj);
            raMemWeights = raMemLengths.*obj.raMemberXC*obj.rDensity;

        end
        
        % Analyzes truss with matrix stiffness method
        % Returns node displacements, reaction forces and member forces    
        function [ raNodeDisplacements,raReactionForces,raMemberForces ] = StiffnessMethod(obj)
            %Number of raNodes
            m = size(obj.raNodes,1);
            %Number of raMembers
            n = size(obj.raMembers,1);
            % Member Lengths 
            raLengths = obj.MemberLengths();

            % Lambda x and y for each member (cos of angle between member (the order of nodes
            % matters) and x or y axis
            raMemXLambda = (obj.raNodes(obj.raMembers(:,2),1)-obj.raNodes(obj.raMembers(:,1),1))./raLengths;
            raMemYLambda = (obj.raNodes(obj.raMembers(:,2),2)-obj.raNodes(obj.raMembers(:,1),2))./raLengths;


            % Member K Matrix Calculation
            X = raMemXLambda * obj.E^.5  .* (obj.raMemberXC.^.5) ./ (raLengths.^.5);
            Y = raMemYLambda *  obj.E^.5  .* (obj.raMemberXC.^.5) ./ (raLengths.^.5);
            % row 1 for each member, then row 2, etc.
            raMemKUnsorted = [ X.^2  X.*Y  -X.^2  -X.*Y;
                             X.*Y  Y.^2  -X.*Y  -Y.^2;
                            -X.^2 -X.*Y   X.^2   X.*Y;
                            -X.*Y -Y.^2   X.*Y   Y.^2];


            % Overall Structure Stiffness Matrix
            raStructK = zeros(2*m,2*m);

            % defines indices of structK to assign values to
            raXCoords = [obj.raMembers(:,1)*2-1 obj.raMembers(:,1)*2 obj.raMembers(:,2)*2-1 obj.raMembers(:,2)*2];
            raXCoords = [raXCoords;raXCoords;raXCoords;raXCoords];
            raYCoords = [obj.raMembers(:,1)*2-1; obj.raMembers(:,1)*2; obj.raMembers(:,2)*2-1; obj.raMembers(:,2)*2];
            raYCoords = [raYCoords raYCoords raYCoords raYCoords];

            % assigns values
            for i = 1:16*n
                raStructK(raXCoords(i),raYCoords(i)) = raStructK(raXCoords(i),raYCoords(i)) + raMemKUnsorted(i);
            end

            % Node Displacements
            raUnrestrainedK = raStructK(1:2*m-obj.nNumReactions,1:2*m-obj.nNumReactions);
            raNodeDisplacements = (raUnrestrainedK)\ obj.raLoads';

            % Reaction Forces
            raRestrainedK = raStructK(2*m-obj.nNumReactions+1:2*m,1:2*m-obj.nNumReactions);
            raReactionForces = raRestrainedK*raNodeDisplacements;

            % Member Forces (positive= tension)
            if obj.nNumReactions == 3
                raNodeDisplacements = [raNodeDisplacements;0;0;0];
            else
                raNodeDisplacements = [raNodeDisplacements;0;0;0;0];
            end 
            raMemberForces = zeros(n,1);

            for i = 1:n
                raMemDelta = [raNodeDisplacements(obj.raMembers(i,1)*2-1); raNodeDisplacements(obj.raMembers(i,1)*2); raNodeDisplacements(obj.raMembers(i,2)*2-1); raNodeDisplacements(obj.raMembers(i,2)*2)];
                raForces = obj.raMemberXC(i)*obj.E/raLengths(i)*[1 -1;-1 1]*[raMemXLambda(i) raMemYLambda(i) 0 0; 0 0 raMemXLambda(i) raMemYLambda(i)]*raMemDelta;
                raMemberForces(i) = raForces(2);
            end
            
        end % function
        
        % Returns an nx1 matrix with member buckling factors of safety
        function [ raBucklingFoS ] = MemberBucklingCheck(obj,memberForces)

            raSmallMat = 1e-10 *ones(size(memberForces,1),1);
            raRemoveTensionForces = abs(min(raSmallMat,memberForces));
            raRadius = (obj.raMemberXC+1/256)*8/pi;
            raBucklingFailureLoads = obj.E.*(raRadius.^4*pi/4)./(obj.MemberLengths().^2);
            raBucklingFoS = raBucklingFailureLoads./raRemoveTensionForces;
        end
        
        % Returns an nx1 matrix with member yielding factors of safety
        function [ raYieldFoS ] = MemberYieldCheck(obj,memberForces)
      
            raMemStress = abs(memberForces)./obj.raMemberXC;
            raYieldFoS = obj.rYieldStrength./raMemStress;
        end
        
        % Sketches original and deformed truss
        function SketchDeformedTruss(obj)
                
            figure;
            hold on;
            axis equal;
            title("Deformed Truss: 100x Deflection");
            xlabel("inches");
            ylabel("inches");
            [nodeDisplacements,~,~] = obj.StiffnessMethod();
            raXYDisplacements = 100*reshape(nodeDisplacements, [2,size(obj.raNodes,1)])';
            raXYNew = obj.raNodes + raXYDisplacements;

            for i = 1:size(obj.raMembers,1)
                plot([obj.raNodes(obj.raMembers(i,1),1) obj.raNodes(obj.raMembers(i,2),1)],[obj.raNodes(obj.raMembers(i,1),2) obj.raNodes(obj.raMembers(i,2),2)],'k');
            end

            for i = 1:size(obj.raMembers,1)
                plot([raXYNew(obj.raMembers(i,1),1) raXYNew(obj.raMembers(i,2),1)],[raXYNew(obj.raMembers(i,1),2) raXYNew(obj.raMembers(i,2),2)],'b--');
            end
                
                
        end % function
        
        % Sketches truss, with differences in member XC
        % area exaggerated
        function SketchTruss(obj)
            
            
            hold off;
            for i = 1:size(obj.raMembers,1)
                plot([obj.raNodes(obj.raMembers(i,1),1) obj.raNodes(obj.raMembers(i,2),1)],[obj.raNodes(obj.raMembers(i,1),2) obj.raNodes(obj.raMembers(i,2),2)],'k', 'LineWidth', (obj.raMemberXC(i)-min(obj.raMemberXC())+.01)*30);
                hold on;
            end
            hold off;
            axis equal;
            title("Truss");
            xlabel("inches");
            ylabel("inches");
            drawnow;
        
        end
        
        % Returns the maximum negative vertical displacement 
        function [rMaxDisp] = MaxVertDisplacement(obj,displacements)
            raEvens = 2:2:size(obj.raNodes,1);
            rMaxDisp = -min(displacements(raEvens));
        end
        
        % Returns the total truss weight
        function [rTrussWeight] = TrussWeight(obj)
            rTrussWeight = sum(obj.MemberWeights());
        end
        
        % Returns the cost of the truss 
            % Contains additional costs in order to impose constraints on
            % minimum member XC area, maximum member length, minimum buckling
            % and yield factors of safety
        function [rCost] = InternalCostFunction(obj)
            [raDisplacements,~,raMemberForces] = obj.StiffnessMethod();

            % Find max - deflection
            rMaxDisp = obj.MaxVertDisplacement(raDisplacements);

            % Find member weights
            rTrussWeight = obj.TrussWeight();

            % Cost
            rCostRaw = rTrussWeight*5000 + rMaxDisp * 3000000;

            %Yield Check
            raYields = obj.y/min(obj.MemberYieldCheck(raMemberForces))^15;
            % Buckling Check
            raBuckles = obj.b/min(obj.MemberBucklingCheck(raMemberForces))^15;
            % Positive Area Check
            raArea = obj.a/min(obj.raMemberXC)^5;
            % Max Length Check
            raLengths = obj.l.*log(exp(10*(obj.MemberLengths()-35.7))+1);
            rLength= sum(raLengths); 
            
            rDepth = max(obj.raNodes(:,2)) - min (obj.raNodes(:,2));
            rDepthCost = obj.d*log(exp(4*(rDepth-14.25))+1);
                     
            rCost = rCostRaw + raYields + raBuckles + raArea +rLength + rDepthCost;
            
        end % function

        % Returns the cost of a truss after updating Truss parameters from
        % input vector 
        % (This method is designed to be used with Optimizer Class to 
        % iteratively update trusses)
        function [rCost] = Cost(obj,raParams)
            % overwrites parameters with input values
            % input values are in a ?x1 vector; for-loops reasign these 
            % values to correct parameters
            k = 1;
            for i = 1:size(obj.raNodeLock,1)
                if obj.raNodeLock(i,1) == 1
                    obj.raNodes(i,1) = raParams(k);
                    k = k+1;
                end
                if obj.raNodeLock(i,2) == 1
                    obj.raNodes(i,2) = raParams(k);
                    k = k+1;
                end
            end
            for j = 1:length(obj.raMemberXC)
                obj.raMemberXC(j) = raParams(k)/obj.rAreaScaling;
                k = k + 1;
            end
            
            rCost = obj.InternalCostFunction();
         
        end % function
        
        % Creates a parameter vector from existing truss configuration
        function [raParams] = GetParameters(obj)
            k = 1;
            for i = 1:length(obj.raNodeLock)
                if obj.raNodeLock(i,1) == 1
                    raParams(k) = obj.raNodes(i,1);
                    k = k + 1;
                end
                if obj.raNodeLock(i,2) == 1
                    raParams(k) = obj.raNodes(i,2);
                    k = k + 1;
                end
            end
            for j = 1:length(obj.raMemberXC)
                raParams(k)= obj.raMemberXC(j)*obj.rAreaScaling;
                k = k + 1;
            end
        end % function
            
    end % methods
    
    
    %------------------------STATIC METHODS-------------------------------%
    
    methods(Static)       
         % Method used to initially test the Truss class
         % Analyzes and optimizes given truss design
         % Final Truss is output
        function [truss1] = TrialDesign()
            
            iSteps = 20000;
            rLearningRate = 5e-6;
            
            raNodes = [...
                
                0   50;
                25  50;
                50  50;
                95  50;
                115 50;
                150 50;
                175 50;
                200 50;
                
               
                25  25;
                50  25;
                75  25;
                100 25;
                125 25;
                150 25;
                175 25;
                0   25;
                200 25]; 
                 
            raNodeLock = ...
                [
                1 1;
                1 1;
                1 1;
                1 0;
                0 1;
                1 1;
                1 1;
                1 1;
                0 0;
                0 0;
                0 0;
                0 0;
                0 0;
                0 0;
                0 0;
                0 0;
                0 0];

            XC1 = 0.3125;
            raMembers = [1  2 XC1;
                       2  3 XC1; 
                       3  4 XC1;
                       4  5 XC1;
                       5  6 XC1;
                       6  7 XC1;
                       7  8 XC1;
                       16 9 XC1;
                       9  10 XC1;
                       10 11 XC1;
                       11 12 XC1;
                       12 13 XC1;
                       13 14 XC1;
                       14 15 XC1;
                       15 17 XC1;
                       1  16 XC1;
                       2  9 XC1;
                       3  10 XC1;
                       4  11 XC1;
                       5  13 XC1;
                       6  14 XC1;
                       7  15 XC1;
                       8  17 XC1;
                       
                       1  9 XC1;
                       2  10 XC1;
                       3  11 XC1;
                       4  12 XC1;
                       5  12 XC1;
                       6  13 XC1;
                       7  14 XC1;
                       8  15 XC1];
                   
            raLoads = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -600 0 -1200 0 -600 0 0 0 0 0];

            truss1 = Truss(raNodes,raNodeLock,raMembers(:,1:2),raMembers(:,3),raLoads);
            truss1.SketchTruss();      
            
            raParameters = truss1.GetParameters();
            [raParameters, raCost] = Optimizer.Run(truss1, raParameters,iSteps, rLearningRate); 

            % Plot cost function over time
            figure, semilogy(raCost);
            title('Cost Trajectory');
            xlabel('Step');
            ylabel('Cost');
            grid on;
            SketchDeformedTruss(truss1);
          

        end % function      
               
   
    end % methods
    
end % classdef

