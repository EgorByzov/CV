classdef (Sealed) TriangleSet < Surface & IPrimitiveSet
    
    properties (Access = {?Serializator})
        triangles;
        tree;
        origFunc;
        origFuncCoord;
    end % properties
    
    properties (Dependent)
        NumPrimitives;
        Triangles;
    end % properties
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        function value = get.NumPrimitives(obj)
            value = size(obj.triangles, 1);
        end % function
        function value = get.Triangles(obj)
            value = obj.triangles;
        end
        
    end % methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Public %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    methods (Access = public)
        
        % constructor
        function obj = TriangleSet(varArg1, varArg2)
			% varArg1 - surface
			if nargin == 1 && Serializator.isSerializedObj(varArg1)
				ValueChecker.checkDeserializationPossibility(obj, varArg1);
				obj = Serializator.load(obj, varArg1);
				obj.tree = [];
            elseif nargin >= 1 && size(varArg1, 2) == 12
                ValueChecker.checkNumerals(varArg1);
                rawTriangs = ValueChecker.checkAndUpdateSize(varArg1, [inf 12]);
                if nargin > 1
                    cs = varArg2;
                else
                    cs = CoordinateSystem.Global;
                end
                
                obj.triangles = TriangleSet.setDefaultTriangNormals(rawTriangs);
                obj.tree = [];
                obj.origFunc = [];
                obj.origFuncCoord = [];
                obj.opticalProperties = OpticalSurfaceProperties.Default;
                obj.cs = cs;
            else
                surf = varArg1;
                apprOpt = varArg2;
                descr = metaclass(surf);
                % get segmented line, original profile function and its type
                switch descr.Name
                    case 'AxisymCubSplinePolar'
                        origFunc = FactorFunc3D([], surf.Profile);
                        origFuncCoord = 'spherical';
                    case 'AxisymCubSplineCart'
                        origFunc = FactorFunc3D([], surf.Profile);
                        origFuncCoord = 'cylindrical';
                    case 'AxisymInterpLineCart'
                        prof = surf.Profile;
                        x = prof(:, 1);
                        z = prof(:, 2);
                        origFunc = FactorFunc3D([], Interp2D(x, z, 'cubic'));
                        origFuncCoord = 'cylindrical';
                    case 'AxisymSegLineCart'
                        origFunc = [];
                        origFuncCoord = [];
                    case 'AxisymSphereSegment'
                        origFunc = FactorFunc3D([], Constant(surf.Radius));
                        origFuncCoord = 'spherical';
                    case 'FreeformBicubSplineCyl'
                        origFunc = surf.Surface;
                        origFuncCoord = 'cylindrical';
                    case 'FreeformBicubSplineSpher'
                        origFunc = surf.Surface;
                        origFuncCoord = 'spherical';
                    case {'FreeformBicubSplineCart',...
                            'FreeformBicubSplineCuttedSpher'}
                        origFunc = [];
                        origFuncCoord = [];
%                         origFunc = surf.Surface;
%                         origFuncCoord = 'cartesian';
                    case 'FreeformCylinderSurface'
                        origFunc = FactorFunc3D(surf.Radius, []);
                        origFuncCoord = 'cylindrical2';
                    case {'Planar1ContourCubSplinePolar', ...
                            'Planar2ContourCubSplinePolar', ...
                            'PlanarQuadrangle'}
                        origFunc = [];
                        origFuncCoord = [];
                    case 'RingSet'
                        origFunc = []; % really, it not perfect solution but I don't care. This case will be never used, I hope.
                        origFuncCoord = [];
                    otherwise
                        throw(MException('TriangleSet:NotSupportedSurface', 'Such surface cannot be represented as RingSet.'));
                end % switch
                
                rawTriangs = surf.getTriangles(apprOpt);
                if size(rawTriangs, 2) == 12
                    if nnz(rawTriangs(1, 10:12)) == 0 % normals were not computed
                        rawTriangs = TriangleSet.setDefaultTriangNormals(rawTriangs);
                    end % if
                else
                    rawTriangs = ValueChecker.checkAndUpdateSize(rawTriangs, [inf 9]);
                    rawTriangs = TriangleSet.setDefaultTriangNormals(rawTriangs);
                end % if
                
                obj.tree = [];                
                obj.triangles = rawTriangs;
                obj.origFunc = origFunc;
                obj.origFuncCoord = origFuncCoord;
                obj.opticalProperties = surf.Properties.clone();
                obj.cs = surf.CS.clone();
			end % if
        end % function
        
        % destructor
        function delete(obj)
            if ~isempty(obj.tree)
                deleteTriangTreeMEX(obj.tree);
                obj.tree = [];
            end % if
        end % function
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% Surface %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function res = getPrimitiveSet(obj, ~)
            res = obj.clone();
        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% PrimitiveSet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [t, ind] = getInter(obj, ray)
            if isempty(obj.tree)
                obj.tree = TriangleSet.makeTriangleTree( obj.triangles );
            end % if
            ray = GeometryConvertor.getLocalRay(ray, obj.cs);
            [t, ind] = getIntersTriangMEX(obj.tree, ray);
        end % function
        
        function numInters = getNumInters(obj, ray)
            ray = GeometryConvertor.getLocalRay(ray, obj.cs);
            numInters = getNumIntersTriangMEX(ray, obj.triangles);
        end
         
        function res = getNormal(obj, ray, t, id)
            ValueChecker.checkConstraint(id, '>=', 1);
            ValueChecker.checkConstraint(id, '<=', obj.NumPrimitives);
            
            if isempty(obj.origFunc)
                % we have only triangles
                res = obj.triangles(id, 10:12);
            else
                ray = GeometryConvertor.getLocalRay(ray, obj.cs);
                % there is a possibility for exact computing normal
                intPoint = ray(:, 1:3) + ray(:, 4:6) .* repmat(t, [1 3]);
                switch obj.origFuncCoord
                    case 'spherical'
                        % compute spherical coordinates of intersection point
                        [phi, psi, r] = cart2sph(   intPoint(:, 1), ...
                                                    intPoint(:, 2), ...
                                                    intPoint(:, 3));
                        phi = mod(phi, 2*pi);
                        psi = pi/2 - psi;
                        % compute normal
                        [r, r_phi, r_psi] = obj.origFunc.eval(phi, psi);
                        sinPhi = sin(phi); sinPsi = sin(psi); cosPhi = cos(phi); cosPsi = cos(psi);
                        rVec_phi(:, 1) = (r_phi .* cosPhi - r .* sinPhi) .* sinPsi;
                        rVec_phi(:, 2) = (r_phi .* sinPhi + r .* cosPhi) .* sinPsi;
                        rVec_phi(:, 3) = r_phi .* cosPsi;
                        rVec_r(:, 1) = (r_psi .* sinPsi + r .* cosPsi) .* cosPhi;
                        rVec_r(:, 2) = (r_psi .* sinPsi + r .* cosPsi) .* sinPhi;
                        rVec_r(:, 3) = r_psi .* cosPsi - r .* sinPsi;
                        nVec = cross(rVec_r, rVec_phi, 2);
                        % norm normal
                        tmp = sqrt(dot(nVec, nVec, 2));
                        res = nVec ./ repmat(tmp, [1 3]);
                    case 'cylindrical'
                        % compute polar base
                        [phi, r] = cart2pol(    intPoint(:, 1), ...
                                                intPoint(:, 2));
                        phi = mod(phi, 2*pi);
                        % compute normal
                        [~, z_phi, z_r] = obj.origFunc.eval(phi, r);
                        sinPhi = sin(phi); cosPhi = cos (phi);
                        rVec_phi(:, 1) = - r .* sinPhi;
                        rVec_phi(:, 2) =   r .* cosPhi;
                        rVec_phi(:, 3) =   z_phi;
                        rVec_r(:, 1) = cosPhi;
                        rVec_r(:, 2) = sinPhi;
                        rVec_r(:, 3) = z_r;
                        nVec = cross(rVec_r, rVec_phi, 2);
                        % norm normal
                        tmp = sqrt(dot(nVec, nVec, 2));
                        res = nVec ./ repmat(tmp, [1 3]);
                    case 'cylindrical2'
                        % compute polar base
                        [phi, ~] = cart2pol(    intPoint(:, 1), ...
                                                intPoint(:, 2));
                        phi = mod(phi, 2*pi);
                        z = intPoint(:, 3);
                        % compute normal
                        [r, r_phi, r_z] = obj.origFunc.eval(phi, z);
                        sinPhi = sin(phi); cosPhi = cos (phi);
                        rVec_phi(:, 1) = r_phi .* cosPhi - r .* sinPhi;
                        rVec_phi(:, 2) = r_phi .* sinPhi + r .* cosPhi;
                        rVec_phi(:, 3) = 0;
                        rVec_r(:, 1) = r_z .* cosPhi;
                        rVec_r(:, 2) = r_z .* sinPhi;
                        rVec_r(:, 3) = 1;
                        nVec = cross(rVec_r, rVec_phi, 2);
                        % norm normal
                        tmp = sqrt(dot(nVec, nVec, 2));
                        res = nVec ./ repmat(tmp, [1 3]);
                    case 'cartesian'
                        % compute polar base
                        x = intPoint(:, 1);
                        y = intPoint(:, 2);
                        [~, z_x, z_y] = obj.origFunc.eval(x, y);
                        nVec = [z_x z_y -ones(size(z_y))];
                        % norm normal
                        tmp = sqrt(dot(nVec, nVec, 2));
                        res = nVec ./ repmat(tmp, [1 3]);
                end % switch
            end % if
            
            res = GeometryConvertor.getGlobalVector(res, obj.cs);
        end % function
                
        function ray = getOutgoingRay(obj)
%             triNum = ceil(size(obj.triangles, 1) / 2);
            meanX = mean([obj.triangles(:, 1);obj.triangles(:, 4);obj.triangles(:, 7)]);
            meanY = mean([obj.triangles(:, 2);obj.triangles(:, 5);obj.triangles(:, 8)]);
            meanZ = mean([obj.triangles(:, 3);obj.triangles(:, 6);obj.triangles(:, 9)]);
            d = (obj.triangles(:, 1)-meanX).^2 + (obj.triangles(:, 2)-meanY).^2 + (obj.triangles(:, 3)-meanZ).^2;
            triNum = find(d == min(d), 1);

            point = (obj.triangles(triNum, 1:3) + ...
                obj.triangles(triNum, 4:6) + ...
                obj.triangles(triNum, 7:9)) / 3;
            ray = [point 0 0 0 0];
            ray = GeometryConvertor.getGlobalRay(ray, obj.cs);
            ray(4:6) = obj.getNormal(ray, 0, triNum);            
        end % function
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% ICloneable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function objNew = clone(obj)
            s = struct('className__', 'TriangleSet', 'triangles', obj.triangles, 'tree', [], ...
                'origFunc', Serializator.save(obj.origFunc), ...
                'origFuncCoord', obj.origFuncCoord, 'cs', Serializator.save(obj.cs),...
                'opticalProperties', Serializator.save(obj.opticalProperties), ...
                'isRegistrator', obj.isRegistrator);
            objNew = TriangleSet(s);
        end % function
                    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% IDrawable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [pointsCloud, lines] = getPointsCloudAndLines(obj)
            lines = {};
            x = [obj.triangles(:, 1) obj.triangles(:, 4) obj.triangles(:, 7)];
            y = [obj.triangles(:, 2) obj.triangles(:, 5) obj.triangles(:, 8)];
            z = [obj.triangles(:, 3) obj.triangles(:, 6) obj.triangles(:, 9)];
            defSize = size(x);
            
            x = reshape(x, [numel(x), 1]);
            y = reshape(y, [numel(y), 1]);
            z = reshape(z, [numel(z), 1]);            
            points = GeometryConvertor.getGlobalPoint([x y z], obj.cs);
            
            pointsCloud(:,:,1) = reshape(points(:,1), defSize);
            pointsCloud(:,:,2) = reshape(points(:,2), defSize);
            pointsCloud(:,:,3) = reshape(points(:,3), defSize);            
        end
        
        function line2D = getProfile2D(~, ~, ~)
            line2D = [];
        end % function

        function line2D = getSection2D(~, ~, ~)
            line2D = [];
        end % function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% IRhinoExportable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function res = getRhinoCommandsString(~, ~, ~)
            % dummy
            res = [];
        end % function
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% ISTLExportable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function res = getTriangles(obj, ~)  
            res = obj.clone();
        end % function

    end % methods  
    
    methods (Static, Access = public)
        
        % creates 'triangles array' for extruded profile surface along Oy
        % profile defined by '2D segmented line'
        %   segLine         - '2D segmented line'
        %   apprOpt         - ApproximationOptions
        %   !!! normals are not computed
        function triangles = segLine2DProfile2ExtrudeTriangles(segLine, extrusionLength)            
            yMin = -extrusionLength/2;
            yMax = extrusionLength/2;
            x = segLine.Points(:,1);
            z = segLine.Points(:,2);
            
            % define maximal number of triangles and allocate memory
            d = sqrt((x(2:end)-x(1:end-1)).^2 + (z(2:end)-z(1:end-1)).^2);
            a = d/sin(pi/3);
            minA = min(a);
            maxTriangNum = ceil(length(x) * extrusionLength/minA * 2);
            triangles = zeros(maxTriangNum, 12);
            
            % build triangles
            curTriangNum = 1;
            
            numPoints1 = round(extrusionLength/a(1)) + 1;
            points1(:,2) = linspace(yMin, yMax,numPoints1);
            points1(:,1) = x(1);
            points1(:,3) = z(1);
            t1 = points1(:,2) - yMin;
                        
            for i = 2 : length(x)
                
                numPoints2 = round(extrusionLength/a(i-1)) + 1;
                points2 = zeros(numPoints2, 3);
                points2(:,2) = linspace(yMin, yMax,numPoints2);
                points2(:,1) = x(i);
                points2(:,3) = z(i);
                t2 = points2(:,2) - yMin;
                
                % pack points
                numInRing = numPoints1 + numPoints2 - 2;
                newTriangles = TriangleSet.twoMeasPointSets2RawTriangles(points1, t1, points2, t2);
                triangles(curTriangNum : curTriangNum + numInRing - 1, 1:9) = newTriangles;
                curTriangNum = curTriangNum + numInRing;
                
                % set new line of points to old
                points1 = points2;
                numPoints1 = numPoints2;
                t1 = t2;
                
            end % for
            triangles = triangles(1 : curTriangNum - 1, :);
            
            % Apply yOz symmetry
            symTriangs = triangles;
            symTriangs(:,1) = -symTriangs(:,1);
            symTriangs(:,4) = -symTriangs(:,4);
            symTriangs(:,7) = -symTriangs(:,7);
            tmp = symTriangs(:,1:3);
            symTriangs(:,1:3) = symTriangs(:,4:6);
            symTriangs(:,4:6) = tmp;
            triangles = [triangles; symTriangs];
        end % function
        
        % creates 'triangles array' for axisymmatrical surface with Oz
        % symmetry axis and profile defined by '2D segmented line'
        %   segLine         - '2D segmented line'
        %   apprOpt         - ApproximationOptions
        %   !!! normals are not computed
        function triangles = segLine2DProfile2RawTriangles(segLine, apprOpt, phiMin, phiMax)
            if nargin == 2
                phiMin = 0;
                phiMax = 2*pi;
            end
            ValueChecker.checkClass(segLine, 'SegLine2D');
            
            maxError = apprOpt.MaxError;
            R = segLine.Points(:,1);
            z = segLine.Points(:,2);
            
            % define maximal number of triangles and allocate memory
            dPhiMin = min(sqrt(8 * maxError ./ R));
            maxTriangNum = ceil(length(R) * (phiMax - phiMin)/dPhiMin * 2);
            triangles = zeros(maxTriangNum, 12);
            
            % build triangles
            curTriangNum = 1;
            [points1, phi1, numPoints1] = TriangleSet.circle2points(R(1), z(1), maxError, phiMin, phiMax);
            if length(phi1) == 2 && (diff(phi1) == 0 || diff(phi1) == 2*pi)
                points1 = points1(1,:);
                phi1 = phi1(1);
                numPoints1 = 1;
            end
            
            for i = 2 : length(R)
                
                % make new line of points
                [points2, phi2, numPoints2] = TriangleSet.circle2points(R(i), z(i), maxError, phiMin, phiMax);
                if length(phi2) == 2 && (diff(phi2) == 0 || diff(phi2) == 2*pi)
                    points2 = points2(1,:);
                    phi2 = phi2(1);
                    numPoints2 = 1;
                end
                
                % pack points
                numInRing = numPoints1 + numPoints2 - 2;
                newTriangles = TriangleSet.twoMeasPointSets2RawTriangles(points1, phi1, points2, phi2);
                triangles(curTriangNum : curTriangNum + numInRing - 1, 1:9) = newTriangles;
                curTriangNum = curTriangNum + numInRing;
                
                % set new line of points to old
                points1 = points2;
                numPoints1 = numPoints2;
                phi1 = phi2;
                
            end % for
            triangles = triangles(1 : curTriangNum - 1, :);
            
        end % function
        
        % creates truncated 'triangles array' for 2 measured sets of points
        % truncated means array with 9 columns only
        %   points1         - 1st 'points array'
        %   meas1           - 1st vector of measures (with the same length as points1)
        %   points2         - 2st 'points array'
        %   meas2           - 2st vector of measures (with the same length as points2)
        function triangles = twoMeasPointSets2RawTriangles(points1, meas1, points2, meas2)
            numPoints1 = size(points1, 1);
            numPoints2 = size(points2, 1);
            
            % allocate memory for pointsInds
            numInRing = numPoints1 + numPoints2 - 2;
            pointsInds = zeros(3*numInRing, 1);
            
            % go
            pointsInds = twoMeasPointSets2TrianglesMEX( pointsInds, meas1, meas2, numPoints1, numPoints2 );
            points = [points1; points2]';
            triangles = reshape(points(:, pointsInds), [9 numInRing])';
        end % function
        
        function triangles = uniformGrid2RawTriangles(x, y, z)
            numRows = size(x, 1);
            numCols = size(x, 2);
            
            % allocate memory for pointsInds
            numTri = (numRows-1) * (numCols - 1) * 2;
            pointsInds = zeros(numTri * 3, 1);
            
            % go
            pointsInds = uniformGrid2TrianglesMEX( pointsInds, numRows, numCols );            
            firstInds = pointsInds(:,1);
            secondInds = pointsInds(:,2);
            thirdInds = pointsInds(:,3);
            triangles = [x(firstInds) y(firstInds) z(firstInds) x(secondInds) y(secondInds) z(secondInds) x(thirdInds) y(thirdInds) z(thirdInds)];
        end % function
        
        % computes normals for triangles
        function tri = setDefaultTriangNormals(tri)
            v1 = tri(:, 7:9) - tri(:, 4:6);
            v2 = tri(:, 7:9) - tri(:, 1:3);
            n = cross(v1, v2);
            n = n ./ repmat(sqrt(dot(n, n, 2)), [1 3]);
            
            tri(:, 10:12) = n;
        end % function        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Private %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Static, Access = private)
        
        % returns C-created binary tree
        function tree = makeTriangleTree( triangles )
            numTriangs = size(triangles, 1);
            aabb = zeros(numTriangs, 6);            
            aabb(:, 1) = min([triangles(:, 1) triangles(:, 4) triangles(:, 7)], [], 2); % min x
            aabb(:, 2) = min([triangles(:, 2) triangles(:, 5) triangles(:, 8)], [], 2); % min y
            aabb(:, 3) = min([triangles(:, 3) triangles(:, 6) triangles(:, 9)], [], 2); % min z            
            aabb(:, 4) = max([triangles(:, 1) triangles(:, 4) triangles(:, 7)], [], 2); % max x
            aabb(:, 5) = max([triangles(:, 2) triangles(:, 5) triangles(:, 8)], [], 2); % max y
            aabb(:, 6) = max([triangles(:, 3) triangles(:, 6) triangles(:, 9)], [], 2); % max z
            rootBox = [min(aabb(:,1:3)) max(aabb(:,4:6))];
            
            maxPrimsInLeaf = 10;
            tree = buildTriangTreeMEX(triangles, aabb, rootBox, maxPrimsInLeaf);
        end % function
                
        function [points, phi, numPoints] = circle2points(R, z, maxError, phiMin, phiMax)
            if nargin == 3
                phiMin = 0;
                phiMax = 2*pi;
            end
            
            if R ~= 0
                sectorSize = (phiMax - phiMin);
                dPhi = sqrt(8 * maxError / R);
                numPhi = ceil(sectorSize / dPhi);
                dPhi = 2*pi / numPhi;
                phi = (phiMin : dPhi : phiMax)';
                points = [R * cos(phi),   R * sin(phi),   repmat(z, size(phi))];
                if sectorSize == 2*pi
                    points(end, :) = points(1, :);
%                 else
%                     points(end+1, :) = points(1, :);
%                     phi(end+1) = 2*pi;
                end
                numPoints = numel(phi);
            else
                points = [0 0 z];
                phi = NaN;
                numPoints = 1;
            end % if
            
        end % function
                
    end % methods
    
end % classdef