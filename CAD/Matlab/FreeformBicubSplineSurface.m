classdef (Abstract) FreeformBicubSplineSurface < FreeformSurface & I3DEditable
    
    properties (Access = public)
        iMax;
    end
    
    properties (Access = {?Serializator, ?FreeformBicubSplineSurface})
        surface;
        optMatr;
        iMin;        
        isMirrored;
    end % properties
    
    properties (Dependent)
        Surface;
        IsMirrored;
    end % properties
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        function value = get.Surface(obj)
            value = obj.surface;
        end % function
        
        function value = get.IsMirrored(obj)
            value = obj.isMirrored;
        end
        function set.IsMirrored(obj, value)
            % get new value of isMirrored
            obj.isMirrored = value;
            
            % update lens
            [f, f_x, f_y, f_xy] = obj.surface.getUniquePart(obj.iMin, obj.iMax);
            obj.surface.updateByUniquePart(obj.iMin, obj.isMirrored, f, f_x, f_y, f_xy);
        end
                    
    end % methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Public %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    methods (Access = public)
        
        function UpdateProfileDerivs(obj, lineInd, isX, f_, f__)
            if ~isX && lineInd == 1
                obj.Surface.UpdateApexDerivs(f_, f__);
            else
                [f, f_x, f_y, f_xy] = obj.surface.getUniquePart(obj.iMin, obj.iMax);
                
                if isX
                    indBase = obj.getBaseAndDependentNodes([1 lineInd]);
                    for i = 1 : numel(f_)
                        f_(i) = obj.getRightAzimuthalValueForUniquePart(f_(i), [i lineInd]);
                        f__(i) = obj.getRightAzimuthalValueForUniquePart(f__(i), [i lineInd]);
                    end                    
                    f_x(:, indBase(2)) = f_;
                    f_xy(:, indBase(2)) = f__;
                else                    
                    f_y(lineInd, 1:obj.iMax) = f_;
                    f_xy(lineInd, 1:obj.iMax) = f__;
                end
                
                obj.surface.updateByUniquePart(obj.iMin, obj.isMirrored, f, f_x, f_y, f_xy);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% Surface %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function primitiveSet = getPrimitiveSet(obj, apprOpt)
            primitiveSet = TriangleSet(obj, apprOpt);        
        end % function
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% I3DEditable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [xBorder, yBorder, indsBorder] = getDefinitionArea(obj)
            xBorder = [obj.surface.X(1) obj.surface.X(obj.iMax)];
            yBorder = [obj.surface.Y(1) obj.surface.Y(end)];
            indsBorder = [numel(obj.surface.Y) obj.iMax];
        end
        
        function nodeInd = getNextXNodeInd(obj, nodeInd)
            if obj.isSurfaceAzimuthalClosed
                numProfiles = length(obj.surface.X) - 1;
            else
                numProfiles = length(obj.surface.X);
            end
            
            if nodeInd(2) == numProfiles
                if obj.isSurfaceAzimuthalClosed
                    nodeInd(2) = 1;
                else
                    nodeInd = NaN;
                end
            else
                nodeInd(2) = nodeInd(2) + 1;
            end
        end % function
                
        function nodeInd = getNextYNodeInd(obj, nodeInd)
            numPoints = length(obj.surface.Y);
            if nodeInd(1) ~= numPoints
                nodeInd(1) = nodeInd(1) + 1;
            else
                nodeInd = NaN;
            end
        end % function
        
        function nodeInd = getPrevXNodeInd(obj, nodeInd)
            if obj.isSurfaceAzimuthalClosed
                numProfiles = length(obj.surface.X) - 1;
            else
                numProfiles = length(obj.surface.X);
            end
            
            if nodeInd(2) == 1
                if obj.isSurfaceAzimuthalClosed
                    nodeInd(2) = numProfiles;
                else
                    nodeInd = NaN;
                end
            else
                nodeInd(2) = nodeInd(2) - 1;
            end
        end % function
        
        function nodeInd = getPrevYNodeInd(~, nodeInd)
            if nodeInd(1) ~= 1
                nodeInd(1) = nodeInd(1) - 1;
            else
                nodeInd = NaN;
            end
        end % function
        
        function res = hasNextXNode(obj, nodeInd)
            if obj.isSurfaceAzimuthalClosed
                res = true;
            else
                if nodeInd(2) == length(obj.surface.X)
                    res = false;
                else
                    res = true;
                end
            end
        end % function
        
        function res = hasNextYNode(obj, nodeInd)
            if nodeInd(1) == length(obj.surface.Y)
                res = false;
            else
                res = true;
            end
        end % function
        
        function res = hasPrevXNode(obj, nodeInd)
            if obj.isSurfaceAzimuthalClosed
                res = true;
            else
                if nodeInd(2) == 1
                    res = false;
                else
                    res = true;
                end
            end
        end % function
        
        function res = hasPrevYNode(~, nodeInd)
            if nodeInd(1) == 1
                res = false;
            else
                res = true;
            end
        end % function
        
        function res = canNodeBeRemoved(~, ~)
            res = false;
        end % function
       
        function res = isNodeXEditable(~, ~)
            res = false;
        end % function

        function res = isNodeYEditable(~, ~)
            res = false;
        end % function

        function res = isNodeFEditable(~, ~)
            res = true;
        end % function

        function res = isNodeF_XEditable(obj, nodeInd)
            indBase = obj.getBaseAndDependentNodes(nodeInd);
            
            
                if obj.iMin == obj.iMax
                    res = false;
                elseif obj.iMax == obj.surface.NumPatchesX
                    res = true;
                elseif indBase(2) == obj.iMin || indBase(2) == obj.iMax
                    res = false;
                else
                    res = true;
                end % if
        end % function

        function res = isNodeF_YEditable(obj, nodeInd)
            indBase = obj.getBaseAndDependentNodes(nodeInd);
            
%             if indBase(1) == 1 && obj.hasCentralPoint()
%                 res = false;
%             else
                res = true;
%             end
        end % function

        function res = isNodeF_XYEditable(obj, nodeInd)
            res = obj.isNodeF_XEditable(nodeInd);
        end % function

        function res = isNodeOptimOn(obj, nodeInd)
            isFOpt = obj.getNodeFOptim(nodeInd);
            isF_XOpt = obj.getNodeF_XOptim(nodeInd);
            isF_YOpt = obj.getNodeF_YOptim(nodeInd);
            isF_XYOpt = obj.getNodeF_XYOptim(nodeInd);
            isMainOpt = obj.getNodeMainOptim(nodeInd);
            if isMainOpt && (isFOpt || isF_XOpt || isF_YOpt || isF_XYOpt)
                res = true;
            else
                res = false;
            end % if
        end % function
        
        function res = areNodesDependent(obj, nodeInd1, nodeInd2)
            baseInd1 = obj.getBaseAndDependentNodes(nodeInd1);
            baseInd2 = obj.getBaseAndDependentNodes(nodeInd2);
            if isequal(baseInd1, baseInd2)
                res = true;
            else
                res = false;
            end % if
        end % function

        function res = getNodeFOptim(obj, nodeInd)
            baseInd = obj.getBaseAndDependentNodes(nodeInd);
            res = obj.optMatr(baseInd(1), baseInd(2), 1);
        end % function
        
        function res = getNodeF_XOptim(obj, nodeInd)
            baseInd = obj.getBaseAndDependentNodes(nodeInd);
            res = obj.optMatr(baseInd(1), baseInd(2), 2);
        end % function
        
        function res = getNodeF_YOptim(obj, nodeInd)
            baseInd = obj.getBaseAndDependentNodes(nodeInd);
            res = obj.optMatr(baseInd(1), baseInd(2), 3);
        end % function
        
        function res = getNodeF_XYOptim(obj, nodeInd)
            baseInd = obj.getBaseAndDependentNodes(nodeInd);
            res = obj.optMatr(baseInd(1), baseInd(2), 4);
        end % function
                
        function res = getNodeMainOptim(obj, nodeInd)
            baseInd = obj.getBaseAndDependentNodes(nodeInd);
            res = obj.optMatr(baseInd(1), baseInd(2), 5);
        end % function
        
        function setNodeFOptim(obj, nodeInd, value)
            baseInd = obj.getBaseAndDependentNodes(nodeInd);
            obj.optMatr(baseInd(1), baseInd(2), 1) = value;
        end % function
        
        function setNodeF_XOptim(obj, nodeInd, value)
            baseInd = obj.getBaseAndDependentNodes(nodeInd);
            obj.optMatr(baseInd(1), baseInd(2), 2) = value;
        end % function
        
        function setNodeF_YOptim(obj, nodeInd, value)
            baseInd = obj.getBaseAndDependentNodes(nodeInd);
            obj.optMatr(baseInd(1), baseInd(2), 3) = value;
        end % function
        
        function setNodeF_XYOptim(obj, nodeInd, value)
            baseInd = obj.getBaseAndDependentNodes(nodeInd);
            obj.optMatr(baseInd(1), baseInd(2), 4) = value;
        end % function
        
        function setNodeMainOptim(obj, nodeInd, value)
            baseInd = obj.getBaseAndDependentNodes(nodeInd);
            obj.optMatr(baseInd(1), baseInd(2), 5) = value;
            if value && sum(obj.optMatr(baseInd(1), baseInd(2), 1:4)) == 0
                if obj.isNodeFEditable(baseInd)
                    obj.optMatr(baseInd(1), baseInd(2), 1) = value;
                end
            end
        end % function
        
        function x = getNodeX(obj, nodeInd)
            x = obj.surface.X(nodeInd(2));
        end % function

        function y = getNodeY(obj, nodeInd)
            y = obj.surface.Y(nodeInd(1));
        end % function

        function f = getNodeF(obj, nodeInd)
            f = obj.surface.F(nodeInd(1), nodeInd(2));
        end % function
        
        function f_x = getNodeF_X(obj, nodeInd)
            f_x = obj.surface.F_X(nodeInd(1), nodeInd(2));
        end % function

        function f_y = getNodeF_Y(obj, nodeInd)
            f_y = obj.surface.F_Y(nodeInd(1), nodeInd(2));
        end % function

        function f_xy = getNodeF_XY(obj, nodeInd)
            f_xy = obj.surface.F_XY(nodeInd(1), nodeInd(2));
        end % function

        function setNodeF(obj, nodeInd, value)
            if obj.isNodeFEditable(nodeInd)
                % update lens
                indBase = obj.getBaseAndDependentNodes(nodeInd);
                [f, f_x, f_y, f_xy] = obj.surface.getUniquePart(obj.iMin, obj.iMax);
                if indBase(1) == 1 && obj.hasCentralPoint()
                    f(1, :) = value;
                else
                    f(indBase(1), indBase(2)) = value;
                end
                
                if obj.hasCentralPoint() && indBase(1) < 3
                    splineTmp = BicubSpline.makeByUniquePart(obj.isMirrored, f, f_x, f_y, f_xy, obj.surface.X, obj.surface.Y);
                    obj.updateCentralPoint(splineTmp);
                else
                    obj.surface.updateByUniquePart(obj.iMin, obj.isMirrored, f, f_x, f_y, f_xy);
                end
            end
        end % function
        
        function setNodeF_X(obj, nodeInd, value)
            if obj.isNodeF_XEditable(nodeInd)
                % update lens
                indBase = obj.getBaseAndDependentNodes(nodeInd);
                [f, f_x, f_y, f_xy] = obj.surface.getUniquePart(obj.iMin, obj.iMax);
                value = obj.getRightAzimuthalValueForUniquePart(value, nodeInd);
                f_x(indBase(1), indBase(2)) = value;
                obj.surface.updateByUniquePart(obj.iMin, obj.isMirrored, f, f_x, f_y, f_xy);
            end
        end % function
        
        function setNodeF_Y(obj, nodeInd, value)
            if obj.isNodeF_YEditable(nodeInd)
                % update lens
                indBase = obj.getBaseAndDependentNodes(nodeInd);
                [f, f_x, f_y, f_xy] = obj.surface.getUniquePart(obj.iMin, obj.iMax);
                f_y(indBase(1), indBase(2)) = value;
                obj.surface.updateByUniquePart(obj.iMin, obj.isMirrored, f, f_x, f_y, f_xy);
            end
        end % function
        
        function setNodeF_XY(obj, nodeInd, value)
            if obj.isNodeF_XYEditable(nodeInd)
                % update lens
                indBase = obj.getBaseAndDependentNodes(nodeInd);
                [f, f_x, f_y, f_xy] = obj.surface.getUniquePart(obj.iMin, obj.iMax);
                value = obj.getRightAzimuthalValueForUniquePart(value, nodeInd);
                f_xy(indBase(1), indBase(2)) = value;
                obj.surface.updateByUniquePart(obj.iMin, obj.isMirrored, f, f_x, f_y, f_xy);
            end
        end % function
                
        function setUniquePart(obj, symm)
            % get new iMax
            switch symm
                case Constants.SurfaceSymmetries.WHOLE
                    phiMax = 2*pi;
                case Constants.SurfaceSymmetries.HALF
                    phiMax = pi;
                case Constants.SurfaceSymmetries.QUARTER
                    phiMax = pi/2;
                case Constants.SurfaceSymmetries.EIGHTH
                    phiMax = pi/4;
                case Constants.SurfaceSymmetries.TWELTH
                    phiMax = pi/6;
            end % switch
            
            if phiMax == 2*pi
                newIMax = length(obj.surface.X) - 1;
            else
                newIMax = find(obj.surface.X == phiMax);
            end % if
            
            % get new part of optimization matrix
            numColumns = newIMax - obj.iMax;
            numRows = size(obj.optMatr,1);
            numLayers = size(obj.optMatr, 3);
            newPartOptMatr = ones(numRows, numColumns, numLayers);
            optMatrCur = [obj.optMatr newPartOptMatr];
            obj.iMax = newIMax;
            obj.optMatr = optMatrCur .* obj.getFullOptMatr();                        
            [f, f_x, f_y, f_xy] = obj.surface.getUniquePart(obj.iMin, obj.iMax);
            obj.surface.updateByUniquePart(obj.iMin, obj.isMirrored, f, f_x, f_y, f_xy);
        end % function
        
        function Scale(obj, factor, fFlag)
            if nargin == 2
                fFlag = false;
            end
            obj.surface.Scale(factor, fFlag);
        end
        
        function Rotate(obj)
            obj.iMax = numel(obj.surface.X);
            obj.optMatr = obj.getFullOptMatr();
            obj.surface.rotate90Sph();
        end
        
        function HalfSymmetry(obj)
            if obj.iMax < 3 || mod(obj.iMax, 2) == 0
                return;
            end
            obj.iMax = (obj.iMax - 1) / 2 + 1;
            oldSpline = obj.surface;
            
            obj.surface = BicubSpline.makeByUniquePart(obj.isMirrored, oldSpline.F(:,1:obj.iMax), oldSpline.F_X(:,1:obj.iMax), oldSpline.F_Y(:,1:obj.iMax), oldSpline.F_XY(:,1:obj.iMax), oldSpline.X, oldSpline.Y);
            
            obj.optMatr = obj.getFullOptMatr();            
        end
        
        function RemoveX(obj, xValue)
            ind = obj.surface.RemoveX(xValue, obj.iMax, obj.isMirrored);
            if ind > 0
                obj.iMax = obj.iMax - 1;
                obj.optMatr = [ obj.optMatr(:,1:ind-1, :) obj.optMatr(:,ind+1:end, :) ];
            end            
        end
        
        function RemoveY(obj, yValue)
            ind = obj.surface.RemoveY(yValue, obj.iMax, obj.isMirrored); 
            if ind > 0
                obj.optMatr = [ obj.optMatr(1:ind-1, :, :); obj.optMatr(ind+1:end, :, :) ];
            end
        end
        
        function AddX(obj, xValue)
            ind = obj.surface.AddX(xValue, obj.iMax, obj.isMirrored);
            if ind > 0
                obj.iMax = obj.iMax + 1;
                matSize = size(obj.optMatr);
                matSize(2) = 1;
                obj.optMatr = [ obj.optMatr(:,1:ind-1, :) zeros(matSize) obj.optMatr(:,ind:end, :) ];
            end            
        end
        
        function AddY(obj, yValue)
            ind = obj.surface.AddY(yValue, obj.iMax, obj.isMirrored);  
            if ind > 0
                matSize = size(obj.optMatr);
                matSize(1) = 1;
                obj.optMatr = [ obj.optMatr(1:ind-1, :, :); zeros(matSize); obj.optMatr(ind:end, :, :) ];
            end
        end
        
        function [symm, curSymm] = getAvailableSymmetries(obj)
            if obj.iMax == (length(obj.surface.X)-1)
                phiMax = obj.surface.X(obj.iMax+1);
            else
                phiMax = obj.surface.X(obj.iMax);
            end % if
            symmsStr = Constants.SurfaceSymmetries;
            
            switch phiMax
                case pi/6
                    symm(1) = symmsStr.WHOLE;
                    symm(2) = symmsStr.HALF;
                    symm(3) = symmsStr.QUARTER;
                    symm(4) = symmsStr.TWELTH;
                    curSymm = symmsStr.TWELTH;
                case pi/4
                    symm(1) = symmsStr.WHOLE;
                    symm(2) = symmsStr.HALF;
                    symm(3) = symmsStr.QUARTER;
                    symm(4) = symmsStr.EIGHTH;
                    curSymm = symmsStr.EIGHTH;
                case pi/2
                    symm(1) = symmsStr.WHOLE;
                    symm(2) = symmsStr.HALF;
                    symm(3) = symmsStr.QUARTER;
                    curSymm = symmsStr.QUARTER;
                case pi
                    symm(1) = symmsStr.WHOLE;
                    symm(2) = symmsStr.HALF;
                    curSymm = symmsStr.HALF;
                case 2*pi
                    symm(1) = symmsStr.WHOLE;
                    curSymm = symmsStr.WHOLE;
            end % switch         
        end % function                
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% IOptimizable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function res = isOptimizationAvailable(obj)
            res = ~isempty( find( sum(obj.optMatr(:,:,1:4), 3) .* obj.optMatr(:,:,5), 1) );
        end % function
        
        function [ values, packParams, LB, UB ] = getOptimParams(obj)
            % unpack & init
            optMatr = obj.optMatr;
            values = [];
            packParams = [];
                        
            % prepare parameters            
            opt_matrF = optMatr(:, :, 1) .* optMatr(:, :, 5);
            opt_matrF_X = optMatr(:, :, 2) .* optMatr(:, :, 5);
            opt_matrF_Y = optMatr(:, :, 3) .* optMatr(:, :, 5);
            opt_matrF_XY = optMatr(:, :, 4) .* optMatr(:, :, 5);
                        
            [f, f_x, f_y, f_xy] = obj.surface.getUniquePart(obj.iMin, obj.iMax);
            
            % pack f(y = 0)
            if obj.hasCentralPoint
                if sum(opt_matrF(1, :))
                    values = [values f(1, 1)];
                    packParams = [packParams; 1 1 1];
                end % if
                rStartInd = 2;
            else
                rStartInd = 1;
            end
            
            % other fs
            [j, i] = find(opt_matrF(rStartInd:end, :));
            j = j + rStartInd - 1;
            for k = 1 : length(j)
                values = [values f(j(k), i(k))];
                packParams = [packParams; 1, j(k), i(k)];
            end % for
            
            % f_xs
            [j, i] = find(opt_matrF_X(rStartInd : end, :));
            j = j + rStartInd - 1;
            for k = 1 : length(j)
                values = [values f_x(j(k), i(k))];
                packParams = [packParams; 2, j(k), i(k)];
            end % for
            
            % f_ys
            [j, i] = find(opt_matrF_Y);
            for k = 1 : length(j)
                values = [values f_y(j(k), i(k))];
                packParams = [packParams; 3, j(k), i(k)];
            end % for
            
            % other f_xys
            [j, i] = find(opt_matrF_XY);
            for k = 1 : length(j)
                values = [values f_xy(j(k), i(k))];
                packParams = [packParams; 4, j(k), i(k)];
            end % for

            % BORDERS!!!
            LB = -Inf(size(values));
            UB = Inf(size(values));
            ind = packParams(:,1) == 1;
            LB(ind) = 0;
            % r
            curDiff = max(max(f)) - min(min(f));
            if curDiff == 0; curDiff = inf; end
            ind = find(packParams(:,1) == 1);
            curMin = values(ind) - 0.5 * curDiff;
            curMin(curMin < 0) = 0;
            LB(ind) = curMin;
            UB(ind) = values(ind) + 0.5 * curDiff;
            % r' by phi
            curDiff = max(max(f_x)) - min(min(f_x));
            if curDiff == 0; curDiff = inf; end
            ind = find(packParams(:,1) == 2);
            LB(ind) = values(ind) - 10 * curDiff;
            UB(ind) = values(ind) + 10 * curDiff;
            % r' by psi
            curDiff = max(max(f_y)) - min(min(f_y));
            if curDiff == 0; curDiff = inf; end
            ind = find(packParams(:,1) == 3);
            LB(ind) = values(ind) - 10 * curDiff;
            UB(ind) = values(ind) + 10 * curDiff;
            % r' by phipsi
            curDiff = max(max(f_xy)) - min(min(f_xy));
            if curDiff == 0; curDiff = inf; end
            ind = find(packParams(:,1) == 4);
            LB(ind) = values(ind) - 10 * curDiff;
            UB(ind) = values(ind) + 10 * curDiff;
            
            
            %% CUSTOM CODE            
%             ind = find(packParams(:,1) == 1);
%             UB(ind) = 9;            
%             ind = find(packParams(:,1) == 3 & packParams(:,2) == size(opt_matrF, 1));
%             LB(ind) = 0;
        end % function
        
        function updateByOptimParams(obj, values, packParams)       
            [f, f_x, f_y, f_xy] = obj.surface.getUniquePart(obj.iMin, obj.iMax);
            
            % set values in f, f_x, f_y, f_xy
            for i = 1 : length(values)
                switch packParams(i, 1)
                    case 1
                        f(packParams(i, 2), packParams(i, 3)) = values(i);
                    case 2
                        f_x(packParams(i, 2), packParams(i, 3)) = values(i);
                    case 3
                        f_y(packParams(i, 2), packParams(i, 3)) = values(i);
                    case 4
                        f_xy(packParams(i, 2), packParams(i, 3)) = values(i);
                end % switch
            end % for
            
            if obj.hasCentralPoint()
                splineTmp = BicubSpline.makeByUniquePart(obj.isMirrored, f, f_x, f_y, f_xy, obj.surface.X, obj.surface.Y);
                obj.updateCentralPoint(splineTmp);
            else
                obj.surface.updateByUniquePart(obj.iMin, obj.isMirrored, f, f_x, f_y, f_xy);
            end
        end % function
        
        function updateSmoothlyByOptimParams(obj, values, packParams)       
            [f, f_x, f_y, f_xy] = obj.surface.getUniquePart(obj.iMin, obj.iMax);
            x = obj.surface.X(obj.iMin:obj.iMax)';
            y = obj.surface.Y;
            
            % set values in f, f_x, f_y, f_xy
            for i = 1 : length(values)
                switch packParams(i, 1)
                    case 1
                        f(packParams(i, 2), packParams(i, 3)) = values(i);
                end % switch
            end % for
            
            f_y(1, :) = (f(2,:) - f(1,:)) ./ repmat(y(2) - y(1), [1 size(f,2)]);
            f_y(end, :) = (f(end,:) - f(end-1,:)) ./ repmat(y(end) - y(end-1), [1 size(f,2)]);
            f_y(2:end-1, :) = (f(3:end,:) - f(1:end-2,:)) ./ repmat(y(3:end) - y(1:end-2), [1 size(f,2)]); 
            
            f_x(:, 1) = (f(:,2) - f(:,1)) ./ repmat(x(2) - x(1), [size(f,1) 1]);
            f_x(:, end) = (f(:, end) - f(:, end-1)) ./ repmat(x(end) - x(end-1), [size(f,1) 1]);
            f_x(:, 2:end-1) = (f(:, 3:end) - f(:, 1:end-2)) ./ repmat(x(3:end) - x(1:end-2), [size(f,1) 1]);
            
            f_xy(:, 1) = (f_y(:,2) - f_y(:,1)) ./ repmat(x(2) - x(1), [size(f,1) 1]);
            f_xy(:, end) = (f_y(:, end) - f_y(:, end-1)) ./ repmat(x(end) - x(end-1), [size(f,1) 1]);
            f_xy(:, 2:end-1) = (f_y(:, 3:end) - f_y(:, 1:end-2)) ./ repmat(x(3:end) - x(1:end-2), [size(f,1) 1]);
            
            if obj.hasCentralPoint()
                splineTmp = BicubSpline.makeByUniquePart(obj.isMirrored, f, f_x, f_y, f_xy, obj.surface.X, obj.surface.Y);
                obj.updateCentralPoint(splineTmp);
            else
                obj.surface.updateByUniquePart(obj.iMin, obj.isMirrored, f, f_x, f_y, f_xy);
            end
        end % function
        
        function turnOnFullOptimization(obj)
            obj.optMatr = obj.getFullOptMatr();
        end % function
        
        function turnOffOptimization(obj)
            obj.optMatr = zeros( obj.surface.NumPatchesY + 1, obj.iMax - obj.iMin + 1, 5 );
        end % function
        
    end % methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PROTECTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = protected)
        
        function value = getRightAzimuthalValueForUniquePart(obj, value, nodeNum)
            
            if obj.iMax < nodeNum(2)              
            t = abs(obj.iMax - obj.iMin);
                tmp = ceil(nodeNum(2) / t);
                if mod(tmp,2) == 0
                    value = - value;
                end % if
            end % if
            
        end % function  
        
        %
        function res = isSurfaceAzimuthalClosed(obj)
            if (obj.surface.X(1) ==  obj.surface.X(end)) || (obj.surface.X(1) ==  obj.surface.X(end) - 2*pi)
                res = true;
            else
                res = false;
            end
        end
        
        %
        function res = hasCentralPoint(obj)
            res = (obj.surface.Y(1) == 0);
        end
        
        %
        function [indBase, isReflectedNode, indDep] = getBaseAndDependentNodes(obj, nodeNum)
            % unpack
            iMinCur = obj.iMin;
            iMaxCur = obj.iMax;
            isMirroredCur = obj.isMirrored;
            numUniquePatchesPhi = iMaxCur - iMinCur;
            numPatchesPhi = obj.surface.NumPatchesX;
            numPatchesPsi = obj.surface.NumPatchesY;
            
            % get base node
            indBasePsi = nodeNum(1);
            if isMirroredCur
                periodSize = (iMaxCur - 1) * 2;
            else
                periodSize = iMaxCur - 1;
            end % if
            if periodSize == 0
                curPeriodNum = 1;
            elseif periodSize == 1
                curPeriodNum = nodeNum(2) - 1;
            else
                curPeriodNum = ceil((nodeNum(2) - 1) / periodSize);
                if curPeriodNum == 0
                    curPeriodNum = 1;
                end % if
            end % if
            indBasePhi = nodeNum(2) - (curPeriodNum - 1) * periodSize;
            if isMirroredCur && indBasePhi > iMaxCur
                indBasePhi = indBasePhi - 2 * (indBasePhi - iMaxCur);
                isReflectedNode = true;
            else
                isReflectedNode = false;
            end % if
            indBase = [indBasePsi indBasePhi];
            
            % get dependent nodes
            if nargout == 3
                % make base matrix
                flagMatr = zeros(numPatchesPsi + 1, numUniquePatchesPhi + 1);
                flagMatr(indBase(1), indBase(2)) = 1;
                
                % handle the case of the 1st column selected and non-mirrored surface
                if ~isMirroredCur   &&   indBase(2) == 1
                    flagMatr(indBase(1), iMaxCur) = 1;
                end % if
                % handle the case of last column selected and non-mirrored surface
                if ~isMirroredCur   &&   indBase(2) == iMaxCur
                    flagMatr(indBase(1), 1) = 1;
                end % if
                
                % extend it like in case of optMatrFull
                % extend matrices if isMirrored
                if isMirroredCur   &&   numUniquePatchesPhi <= numPatchesPhi / 2
                    flagMatr = [flagMatr   fliplr(flagMatr(:, 1:end-1))];
                    numUniquePatchesPhi = numUniquePatchesPhi * 2;
                end
                
                % check case when f is already full (but without last column)
                if numUniquePatchesPhi < numPatchesPhi - 1
                    periods = round(numPatchesPhi / numUniquePatchesPhi);
                    if periods > 1
                        flagMatr = [flagMatr(:,1)    repmat(flagMatr(:,2:end), [1 periods])];
                    end % if
                end % if
                
                % cut last column if it is neccessary
                if size(flagMatr, 2) == numPatchesPhi + 1
                    flagMatr = flagMatr(:, 1:end-1);
                end % if
                
                flagMatr(nodeNum(1), nodeNum(2)) = 0;
                
                [j, i] = find(flagMatr);
                indDep = [j i];
                
            end % if
        end % function
                
        %
        function splineChangedListener(obj, ~, eventData)
            bckp = Serializator.save(obj);
            if isa(eventData, 'GeometryEventData')
                bckp.surface = eventData.bckpStruct;
            end
            notify(obj,'SurfaceChanged', GeometryEventData(bckp));
        end
        
        % returns maximally filled optimization bool matrix
        function optMatr = getFullOptMatr(obj)
            optMatr = ones( obj.surface.NumPatchesY + 1, obj.iMax - obj.iMin + 1, 5 );
            for i = 1:size(optMatr, 1)
                for j = 1:size(optMatr, 2)
                    isF = obj.isNodeFEditable([i j]);
                    isF_X = obj.isNodeF_XEditable([i j]);
                    isF_Y = obj.isNodeF_YEditable([i j]);
                    isF_XY = obj.isNodeF_XYEditable([i j]);
                    optMatr(i, j, 1:4) = optMatr(i, j, 1:4) .* reshape([isF, isF_X, isF_Y, isF_XY], [1 1 4]);
                    optMatr(i, j, 5) = optMatr(i, j, 5) * (isF || isF_X || isF_Y || isF_XY);
                end
            end
        end % function
        
    end % methods
    
    methods (Abstract, Access = public)
        
        segLine = GetBottomContour(obj);
        
    end
    
    methods (Abstract, Access = protected)
        
        updateCentralPoint(obj, spline);
        
    end
    
end % classdef