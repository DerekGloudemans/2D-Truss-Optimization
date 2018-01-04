classdef Truss < matlab.mixin.Copyable %handle
    % Truss a 2D simple truss on which stiffness matrix method analysis and
    % related functions can be performed
    
    properties
        %   nodes- an mx2 matrix of xcoords and ycoords for each node
        %   members- an nx2 matrix of start and end node for each member
        %   memberXC- an nx1 matrix of XC area for each member
        %   nodeLock- an mx2 matrix which indicates whether the position of the
        %             node may be altered along each degree of freedom (1 = yes, 0 = no)
        %   loads- a 1x(2m-4) matrix of loads applied to each degree of freedom
        %          degrees of freedom for node 1 are 1 and 2, 3 and 4 for node
        %          2, etc.
        nodes;
        nodeLock;
        memberXC;
        members;
        loads;
        
    end
    
    properties(Constant)
        % areaScaling- scales XC areas so that scale is similar to that of
            % XY node coordinates (for better gradient descent performance)
        % numReactions- set to 3 for pin-roller, or 4 for pin-pin truss
        areaScaling = 100;
        numReactions = 3;
        
        % 1018 Steel   
            % E- Elastic Modulus, psi
            % density- lb/in^3
            % yieldStrength- psi
        E = 2.9 * 10^7; 
        density = 0.284; 
        yieldStrength = 3.6*10^4; 
        
        % Weighting Coefficients
        % a- controls the extent to which min XC area affects cost function
        % b- controls the extent to which min buckling FoS area affects cost function
        % y- controls the extent to which min yield FoS affects cost function
        % l- controls the extent to which max length affects cost function
        % d- controls the extent to which max depth affects cost function
        a = 10;
        b = 100000;
        y = 100000;
        l = 100000;
        d = 0;%10000;
        
    end
              
    methods
        
        % Returns a truss object with user-input nodes, nodeLock,
            % members, memberXCs and loads
        function [obj] = Truss(n,nl,m,mxc,lo)
            obj.nodes = n;
            obj.nodeLock = nl;
            obj.memberXC = mxc;
            obj.members = m;
            obj.loads = lo;
           
        end
        
        % Returns a nx1 matrix of member lengths
        function [ memLengths ] = MemberLengths (obj)
            memLengths = ((obj.nodes(obj.members(:,1),1)-obj.nodes(obj.members(:,2),1)).^2 + (obj.nodes(obj.members(:,1),2)-obj.nodes(obj.members(:,2),2)).^2).^.5;
        end
        
        % Returns a nx1 matrix of member weights
        function [ memWeights ] = MemberWeights(obj)
   
            memLengths = MemberLengths(obj);
            memWeights = memLengths.*obj.memberXC*obj.density;

        end
        
        % Analyzes truss with matrix stiffness method
            % Returns node displacements, reaction forces and member forces           
        function [ nodeDisplacements,reactionForces,memberForces ] = StiffnessMethod(obj)
            %Number of nodes
            m = size(obj.nodes,1);
            %Number of members
            n = size(obj.members,1);
            % Member Lengths 
            memLength = obj.MemberLengths();

            % Lambda x and y for each member (cos of angle between member (the direction of nodes
            % matters) and x or y axis
            memXLambda = (obj.nodes(obj.members(:,2),1)-obj.nodes(obj.members(:,1),1))./memLength;
            memYLambda = (obj.nodes(obj.members(:,2),2)-obj.nodes(obj.members(:,1),2))./memLength;


            % Member K Matrix Calculation
            X = memXLambda * obj.E^.5  .* (obj.memberXC.^.5) ./ (memLength.^.5);
            Y = memYLambda *  obj.E^.5  .* (obj.memberXC.^.5) ./ (memLength.^.5);
            % row 1 for each member, then row 2, etc.
            memKUnsorted = [ X.^2  X.*Y  -X.^2  -X.*Y;
                             X.*Y  Y.^2  -X.*Y  -Y.^2;
                            -X.^2 -X.*Y   X.^2   X.*Y;
                            -X.*Y -Y.^2   X.*Y   Y.^2];


            % Overall Structure Stiffness Matrix
            structK = zeros(2*m,2*m);

            % defines indices of structK to assign values to
            xCoords = [obj.members(:,1)*2-1 obj.members(:,1)*2 obj.members(:,2)*2-1 obj.members(:,2)*2];
            xCoords = [xCoords;xCoords;xCoords;xCoords];
            yCoords = [obj.members(:,1)*2-1; obj.members(:,1)*2; obj.members(:,2)*2-1; obj.members(:,2)*2];
            yCoords = [yCoords yCoords yCoords yCoords];

            % assigns values
            for i = 1:16*n
                structK(xCoords(i),yCoords(i)) = structK(xCoords(i),yCoords(i)) + memKUnsorted(i);
            end

            % Node Displacements
            unrestrainedK = structK(1:2*m-obj.numReactions,1:2*m-obj.numReactions);
            nodeDisplacements = (unrestrainedK)\ obj.loads';

            % Reaction Forces
            restrainedK = structK(2*m-obj.numReactions+1:2*m,1:2*m-obj.numReactions);
            reactionForces = restrainedK*nodeDisplacements;

            % Member Forces (positive= tension)
            if obj.numReactions == 3
                nodeDisplacements = [nodeDisplacements;0;0;0];
            else
                nodeDisplacements = [nodeDisplacements;0;0;0;0];
            end 
            memberForces = zeros(n,1);

            for i = 1:n
                memDelta = [nodeDisplacements(obj.members(i,1)*2-1); nodeDisplacements(obj.members(i,1)*2); nodeDisplacements(obj.members(i,2)*2-1); nodeDisplacements(obj.members(i,2)*2)];
                forces = obj.memberXC(i)*obj.E/memLength(i)*[1 -1;-1 1]*[memXLambda(i) memYLambda(i) 0 0; 0 0 memXLambda(i) memYLambda(i)]*memDelta;
                memberForces(i) = forces(2);
            end
            
        end % function
        
        % Returns an nx1 matrix with member buckling factors of safety
        function [ bucklingFoS ] = MemberBucklingCheck(obj,memberForces)

            smallMat = 1e-10 *ones(size(memberForces,1),1);
            removeTensionForces = abs(min(smallMat,memberForces));
            radius = (obj.memberXC+1/256)*8/pi;
            bucklingFailureLoads = obj.E.*(radius.^4*pi/4)./(obj.MemberLengths().^2);
            bucklingFoS = bucklingFailureLoads./removeTensionForces;
        end
        
        % Returns an nx1 matrix with member yielding factors of safety
        function [ yieldFoS ] = MemberYieldCheck(obj,memberForces)
      
            memStress = abs(memberForces)./obj.memberXC;
            yieldFoS = obj.yieldStrength./memStress;
        end
        
        % Removes one member from the member Matrix
        function DeleteMembers(obj,rMemIndex)
            
            obj.members(rMemIndex,:) = [];
            obj.memberXC(rMemIndex,:) = [];
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
            xyDisplacements = 100*reshape(nodeDisplacements, [2,size(obj.nodes,1)])';
            xyNew = obj.nodes + xyDisplacements;

            for i = 1:size(obj.members,1)
                plot([obj.nodes(obj.members(i,1),1) obj.nodes(obj.members(i,2),1)],[obj.nodes(obj.members(i,1),2) obj.nodes(obj.members(i,2),2)],'k');
            end

            for i = 1:size(obj.members,1)
                plot([xyNew(obj.members(i,1),1) xyNew(obj.members(i,2),1)],[xyNew(obj.members(i,1),2) xyNew(obj.members(i,2),2)],'b--');
            end
                
                
        end % function
        
        % Sketches truss, with differences in member XC
        % area exaggerated
        function SketchTruss(obj)
            
            
            hold off;
            for i = 1:size(obj.members,1)
                plot([obj.nodes(obj.members(i,1),1) obj.nodes(obj.members(i,2),1)],[obj.nodes(obj.members(i,1),2) obj.nodes(obj.members(i,2),2)],'k', 'LineWidth', (obj.memberXC(i)-min(obj.memberXC())+.01)*30);
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
        function [maxDisp] = MaxVertDisplacement(obj,displacements)
            evens = 2:2:size(obj.nodes,1);
            maxDisp = -min(displacements(evens));
        end
        
        % Returns the total truss weight
        function [TrussWeight] = TrussWeight(obj)
            TrussWeight = sum(obj.MemberWeights());
        end
        
        % Returns the cost of the truss according to ASCE judging
            % Contains additional costs in order to impose constraints on
            % minimum member XC area, maximum member length, minimum buckling
            % and yield factors of safety
        function [cost] = InternalCostFunction(obj)
            [displacements,~,memberForces] = obj.StiffnessMethod();

            % Find max - deflection
            maxDisp = obj.MaxVertDisplacement(displacements);

            % Find member weights
            TrussWeight = obj.TrussWeight();

            % Cost
            costRaw = TrussWeight*5000 + maxDisp * 3000000;

  %%%%          %Yield Check
            yields = obj.y/min(obj.MemberYieldCheck(memberForces))^15;
            % Buckling Check
            buckles = obj.b/min(obj.MemberBucklingCheck(memberForces))^15;
            % Positive Area Check
            area = 0; %obj.a/min(obj.memberXC)^5;
            % Max Length Check
            lengths = obj.l.*log(exp(10*(obj.MemberLengths()-35.7))+1);
            length= sum(lengths); 
            
            depth = max(obj.nodes(:,2)) - min (obj.nodes(:,2));
            depthCost = obj.d*log(exp(4*(depth-14.25))+1);
                     
            cost = costRaw + yields + buckles + area +length + depthCost;
            
        end % function

        % Returns the cost of a truss given a vector of parameters to update
        function [cost] = Cost(obj,params)
            % overwrites parameters with input values
            k = 1;
            for i = 1:size(obj.nodeLock,1)
                if obj.nodeLock(i,1) == 1
                    obj.nodes(i,1) = params(k);
                    k = k+1;
                end
                if obj.nodeLock(i,2) == 1
                    obj.nodes(i,2) = params(k);
                    k = k+1;
                end
            end
            for j = 1:length(obj.memberXC)
                obj.memberXC(j) = params(k)/obj.areaScaling;
                k = k + 1;
            end
            
            cost = obj.InternalCostFunction();
         
        end % function
        
        % Creates a parameter vector from existing truss configuration
        function [params] = GetParameters(obj)
            k = 1;
            for i = 1:length(obj.nodeLock)
                if obj.nodeLock(i,1) == 1
                    params(k) = obj.nodes(i,1);
                    k = k + 1;
                end
                if obj.nodeLock(i,2) == 1
                    params(k) = obj.nodes(i,2);
                    k = k + 1;
                end
            end
            for j = 1:length(obj.memberXC)
                params(k)= obj.memberXC(j)*obj.areaScaling;
                k = k + 1;
            end
        end % function
            
    end % methods
    
    
    methods(Static)
        
        % Analyzes and optimizes given truss design
        function Design()
            
            iSteps = 20000;
            rLearningRate = 5e-6;
            
            nodes = [0   40;
                     40  0 ;
                     80  40;
                     80  0 ;
                     120 40;
                     0   0 ;
                     120 0 ]; 
                 
            nodeLock = ...
                [
                0;
                1;
                1;
                1;
                0;
                0;
                0];

            XC1 = 0.3125;
            members = [1 3 XC1;
                       3 5 XC1; 
                       1 6 XC1;
                       1 2 XC1;
                       2 3 XC1;
                       3 4 XC1;
                       3 7 XC1;
                       5 7 XC1;
                       6 2 XC1;
                       2 4 XC1;
                       4 7 XC1];
            loads = [0 -2000 0 0 0 0 0 0 -2000 1000];

            truss1 = Truss(nodes,nodeLock,members(:,1:2),members(:,3),loads);
            % l1 = MemberLengths(truss1);
            % w1 = MemberWeights(truss1);
            % [a,b,c] = truss1.StiffnessMethod();
            % buck = truss1.MemberBucklingCheck(c);
            % yield = truss1.MemberYieldCheck(c);
            % truss1.SketchTruss();
            % truss1.SketchDeformedTruss();
            % cost = truss1.Cost();

            
            raParameters = truss1.GetParameters();
            [raParameters, raCost] = Optimizer.Run(truss1, raParameters,iSteps, rLearningRate); 

            % Plot cost function over time
            figure, semilogy(raCost);
            title('Cost Trajectory');
            xlabel('Step');
            ylabel('Cost');
            grid on;
            raParameters 
            SketchDeformedTruss(truss1);
%             cost = truss1.Cost(params);

        end % function
        
        % Analyzes and optimizes given truss design
        function [truss1] = Design2()
            
            iSteps = 20000;
            rLearningRate = 5e-6;
            
            nodes = [...
                0   24;
                30  24;
                60  24;
                90 24;
                120 24;
                150 24;
                180 24;
                210 24;
                
                30  0;
                60  0;
                90  0;
                120 0;
                150 0;
                180 0;
                0   0;
                210 0]; 
                 
            nodeLock = ...
                [
                1;
                1;
                0;
                0;
                0;
                0;
                1;
                1;
                1;
                1;
                1;
                1;
                1;
                1;
                0;
                0];

            XC1 = 0.3125;
            members = [1  2 XC1;
                       2  3 XC1; 
                       3  4 XC1;
                       4  5 XC1;
                       5  6 XC1;
                       6  7 XC1;
                       7  8 XC1;
                       15 9 XC1;
                       9  10 XC1;
                       10 11 XC1;
                       11 12 XC1;
                       12 13 XC1;
                       13 14 XC1;
                       14 16 XC1;
                       1  15 XC1;
                       2  9 XC1;
                       3  10 XC1;
                       4  11 XC1;
                       5  12 XC1;
                       6  13 XC1;
                       7  14 XC1;
                       8  16 XC1;
                       1  9 XC1;
                       2  10 XC1;
                       3  11 XC1;
                       4  12 XC1;
                       5  11 XC1;
                       6  12 XC1;
                       7  13 XC1;
                       8  14 XC1];
            loads = [0 0 0 0 0 -600 0 -600 0 -600 0 -600 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

            truss1 = Truss(nodes,nodeLock,members(:,1:2),members(:,3),loads);
                        
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
        
        % Analyzes and optimizes given truss design
        function [truss1] = Design3()
            
            iSteps = 20000;
            rLearningRate = 5e-6;
            
            nodes = [...
                0   25;
                25  25;
                50  25;
                75  25;
                100 25;
                125 25;
                150 25;
                175 25;
                200 25;
                
                25  0;
                50  0;
                75  0;
                100 0;
                125 0;
                150 0;
                175 0;
                0   0;
                200 0]; 
                 
            nodeLock = ...
                [
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                1;
                1;
                1;
                1;
                1;
                1;
                1;
                0;
                0];

            XC1 = 0.3125;
            members = [1  2 XC1;
                       2  3 XC1; 
                       3  4 XC1;
                       4  5 XC1;
                       5  6 XC1;
                       6  7 XC1;
                       7  8 XC1;
                       8  9 XC1;
                       17 10 XC1;
                       10 11 XC1;
                       11 12 XC1;
                       12 13 XC1;
                       13 14 XC1;
                       14 15 XC1;
                       15 16 XC1;
                       16 18 XC1;
                       1  17 XC1;
                       2  10 XC1;
                       3  11 XC1;
                       4  12 XC1;
                       5  13 XC1;
                       6  14 XC1;
                       7  15 XC1;
                       8  16 XC1;
                       9  18 XC1;
                       1  10 XC1;
                       2  11 XC1;
                       3  12 XC1;
                       4  13 XC1;
                       6  13 XC1;
                       7  14 XC1;
                       8  15 XC1;
                       9  16 XC1];
                   
            loads = [0 0 0 0 0 0 0 -600 0 -1200 0 -600 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

            truss1 = Truss(nodes,nodeLock,members(:,1:2),members(:,3),loads);
                        
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
        
        % Analyzes and optimizes given truss design
        function [truss1] = Design4()
            
            iSteps = 20000;
            rLearningRate = 5e-6;
            
            nodes = [...
                0   25;
                25  25;
                50  25;
                75  25;
                100 25;
                125 25;
                150 25;
                175 25;
                200 25;
                
                25  0;
                50  0;
                75  0;
                100 0;
                125 0;
                150 0;
                175 0;
                0   0;
                200 0]; 
                 
            nodeLock = ...
                [
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                1;
                1;
                1;
                1;
                1;
                1;
                1;
                0;
                0];

            XC1 = 0.3125;
            members = [1  2 XC1;
                       2  3 XC1; 
                       3  4 XC1;
                       4  5 XC1;
                       5  6 XC1;
                       6  7 XC1;
                       7  8 XC1;
                       8  9 XC1;
                       17 10 XC1;
                       10 11 XC1;
                       11 12 XC1;
                       12 13 XC1;
                       13 14 XC1;
                       14 15 XC1;
                       15 16 XC1;
                       16 18 XC1;
                       1  17 XC1;
                       2  10 XC1;
                       3  11 XC1;
                       4  12 XC1;
                       5  13 XC1;
                       6  14 XC1;
                       7  15 XC1;
                       8  16 XC1;
                       9  18 XC1;
                       
                       2  17 XC1;
                       3  10 XC1;
                       4  11 XC1;
                       5  12 XC1;
                       5  13 XC1;
                       5  14 XC1;
                       6  15 XC1;
                       7  16 XC1;
                       8  18 XC1];
                   
            loads = [0 0 0 0 0 0 0 -600 0 -1200 0 -600 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ];

            truss1 = Truss(nodes,nodeLock,members(:,1:2),members(:,3),loads);
                        
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
        
         % Analyzes and optimizes given truss design
        function [truss1] = Design5()
            
            iSteps = 20000;
            rLearningRate = 5e-6;
            
            nodes = [...
                0   25;
                25  25;
                50  25;
                75  25;
                100 25;
                125 25;
                150 25;
                175 25;
                200 25;
                
                25  0;
                50  0;
                85  0;
                115 0;
                150 0;
                175 0;
                0   0;
                200 0]; 
                 
            nodeLock = ...
                [
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                0;
                1;
                1;
                1;
                1;
                1;
                1;
                0;
                0];

            XC1 = 0.3125;
            members = [1  2 XC1;
                       2  3 XC1; 
                       3  4 XC1;
                       4  5 XC1;
                       5  6 XC1;
                       6  7 XC1;
                       7  8 XC1;
                       8  9 XC1;
                       16 10 XC1;
                       10 11 XC1;
                       11 12 XC1;
                       12 13 XC1;
                       13 14 XC1;
                       14 15 XC1;
                       15 17 XC1;
                       1  16 XC1;
                       2  10 XC1;
                       3  11 XC1;
                       4  12 XC1;
                       6  13 XC1;
                       7  14 XC1;
                       8  15 XC1;
                       9  17 XC1;
                       
                       2  16 XC1;
                       3  10 XC1;
                       4  11 XC1;
                       5  12 XC1;
                       5  13 XC1;
                       6  14 XC1;
                       7  15 XC1;
                       8  17 XC1];
                   
            loads = [0 0 0 0 0 0 0 -600 0 -1200 0 -600 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

            truss1 = Truss(nodes,nodeLock,members(:,1:2),members(:,3),loads);
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
        
         % Analyzes and optimizes given truss design
        function [truss1] = Design6()
            
            iSteps = 20000;
            rLearningRate = 5e-6;
            
            nodes = [...
                
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
                 
            nodeLock = ...
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
            members = [1  2 XC1;
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
                   
            loads = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -600 0 -1200 0 -600 0 0 0 0 0];

            truss1 = Truss(nodes,nodeLock,members(:,1:2),members(:,3),loads);
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
        
         % Analyzes and optimizes given truss design
        function [truss1] = Design7()
            
            iSteps = 20000;
            rLearningRate = 5e-6;
            
            nodes = [...
                
                0   50;
                25  50;
                50  50;
                75  50;
                           
                37.5 0;
                0    0;
                75   0]; 
                 
            nodeLock = ...
                [
                0 0;
                0 0;
                0 0;
                0 1;
                1 0;
                0 0;
                0 0];

            XC1 = 0.3125;
            members = [1  2 XC1;
                       2  3 XC1; 
                       3  4 XC1;
                       1  6 XC1;
                       4  7 XC1;
                       5  6 XC1;
                       5  7 XC1;
                       1  5 XC1;
                       2  5 XC1;
                       3  5 XC1;
                       4  5 XC1];
                       
                   
            loads = [0 0 0 -500 0 -500 0 0 0 0 0];
            
            truss1 = Truss(nodes,nodeLock,members(:,1:2),members(:,3),loads);
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
        
    % Analyzes and optimizes given truss design
        function [truss1] = Design8()
            
            iSteps = 20000;
            rLearningRate = 5e-6;
            
            nodes = [...
                
                0   50;
                25  50;
                50  50;
                75  50;
                           
                37.5 0;
                0    0;
                75   0]; 
                 
            nodeLock = ...
                [
                0 0;
                0 0;
                0 0;
                0 0;
                1 1;
                0 0;
                0 0];

            XC1 = 0.3125;
            members = [1  2 XC1;
                       2  3 XC1; 
                       3  4 XC1;
                       1  6 XC1;
                       4  7 XC1;
                       5  6 XC1;
                       5  7 XC1;
                       1  5 XC1;
                       2  5 XC1;
                       3  5 XC1;
                       4  5 XC1];
                       
                   
            loads = [0 0 0 -500 0 -500 0 0 0 0 0];
            
            truss1 = Truss(nodes,nodeLock,members(:,1:2),members(:,3),loads);
            truss1.SketchTruss();      
            
            raParameters = truss1.GetParameters();
            [raParameters, raCost] = Optimizer.Run2(truss1, raParameters,iSteps, rLearningRate); 

            % Plot cost function over time
            figure, semilogy(raCost);
            title('Cost Trajectory');
            xlabel('Step');
            ylabel('Cost');
            grid on;
            SketchDeformedTruss(truss1);
          

        end % function
        
    % Analyzes and optimizes given truss design
        function [truss1] = Design9()
            
            iSteps = 20000;
            rLearningRate = 5e-6;
            
            nodes = [...
                
                0   25;
                25  25;
                50  25;
                25  0;
                           
                13 -25;
                37 -25;
                0   0;
                50  0]; 
                 
            nodeLock = ...
                [
                1 0;
                1 0;
                1 0;
                1 0;
                1 0;
                1 0;
                0 0;
                0 0];

            XC1 = 0.3125;
            members = [1  2 XC1;
                       2  3 XC1; 
                       3  4 XC1;
                       1  4 XC1;
                       2  7 XC1;
                       2  8 XC1;
                       1  7 XC1;
                       2  4 XC1;
                       3  8 XC1;
                       7  4 XC1;
                       4  8 XC1;
                       1  5 XC1;
                       2  5 XC1;
                       4  5 XC1;
                       4  6 XC1;
                       2  6 XC1;
                       3  6 XC1];
                       
                   
            loads = [0 0 0 -500 0 0 0 0 0 0 0 0 0];
            
            truss1 = Truss(nodes,nodeLock,members(:,1:2),members(:,3),loads);
            truss1.SketchTruss();      
            
            raParameters = truss1.GetParameters();
            [raParameters, raCost] = Optimizer.Run2(truss1, raParameters,iSteps, rLearningRate); 

            % Plot cost function over time
            figure, semilogy(raCost);
            title('Cost Trajectory');
            xlabel('Step');
            ylabel('Cost');
            grid on;
            SketchDeformedTruss(truss1);
          

        end % function
        
        
    % Analyzes and optimizes truss given input truss design. The truss is
    % constrained to a trapezoidal shape and not subjected to member length
    % limitations.
        function [truss1] = PreliminaryASCEDesign1()
            
            iSteps = 20000;
            rLearningRate = 5e-6;
            
            nodes = [...
                
                14  14;
                42  14;
                70  14;
                98  14;
                126 14;
                154 14;
                182 14;
                                                             
                28   0;
                56   0;
                84   0;
                112  0;
                140  0;
                168  0;
                
                0    7
                196  7]; 
                 
            nodeLock = ...
                [
                1 0;
                1 0;
                0 0;
                0 0;
                0 0;
                1 0;
                1 0;
               
                1 0;
                1 0;
                1 0;
                1 0;
                1 0;
                1 0;
                                
                0 1;
                0 1];

            XC1 = 0.3125;
            members = [14 1 XC1;
                       1  2 XC1; 
                       2  3 XC1;
                       3  4 XC1;
                       4  5 XC1;
                       5  6 XC1;
                       6  7 XC1;
                       7 15 XC1;
                       14 8 XC1;
                       8  9 XC1;
                       9 10 XC1;
                       10 11 XC1;
                       11 12 XC1;
                       12 13 XC1;
                       13 15 XC1;
                       1   8 XC1;
                       8   2 XC1;
                       2   9 XC1;
                       9   3 XC1;
                       3  10 XC1;
                       10  4 XC1;
                       4  11 XC1;
                       11  5 XC1;
                       5  12 XC1;
                       12  6 XC1;
                       6  13 XC1;
                       13  7 XC1];
                   
            loads = [0 0 ...
                     0 0 ...
                     0 -600 ...
                     0 -1200 ...
                     0 -600 ...
                     0 0 ...
                     0 0 ...
                     ...
                     0 0 ...
                     0 0 ...
                     0 0 ...
                     0 0 ...
                     0 0 ...
                     0 0 ...
                     0];

            truss1 = Truss(nodes,nodeLock,members(:,1:2),members(:,3),loads);

            truss1.SketchDeformedTruss();      

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

