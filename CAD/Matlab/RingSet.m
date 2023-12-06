classdef (Sealed) RingSet < Surface & IPrimitiveSet
    
    properties (Access = {?Serializator})
        tree;
        r;
        z;
        k;
        b;
        nr;
        nz;
        origFunc;
        origFuncCoord;
    end % properties
    
    properties(Access = private, Transient)
        goodRings;
    end
    
    properties (Dependent)
        NumPrimitives;
    end % properties
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        function value = get.NumPrimitives(obj)
            value = length(obj.k);
        end % function
        
    end % methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Public %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = public)
        
        % constructor
        function obj = RingSet(varArg1, apprOpt)
			% varArg1 - rings
			if nargin == 1 && Serializator.isSerializedObj(varArg1)
				ValueChecker.checkDeserializationPossibility(obj, varArg1);
				obj = Serializator.load(obj, varArg1);
				obj.tree = [];
                obj.goodRings = zeros(obj.NumPrimitives, 4);
                for i = 1:obj.NumPrimitives
                    obj.goodRings(i, :) = [obj.r(i), obj.z(i), obj.r(i + 1), obj.z(i + 1)];
                end % for
            else
                surf = varArg1;
                descr = metaclass(surf);
                % get segmented line, original profile function and its type
                switch descr.Name
                    case 'AxisymCubSplinePolar'
                        origFunc = surf.Profile;
                        origFuncCoord = 'polar';
                        segLine = SegLine2D.makeFromPolarCubSpline(origFunc, apprOpt);
                    case 'AxisymCubSplineCart'
                        origFunc = surf.Profile;
                        origFuncCoord = 'cart';
                        segLine = SegLine2D.makeFromCartCubSpline(origFunc, apprOpt);
                    case 'AxisymInterpLineCart'
                        prof = surf.Profile;
                        x = prof(:, 1);
                        z = prof(:, 2);
                        origFunc = Interp2D(x, z, 'cubic');
                        origFuncCoord = 'cart';
                        segLine = SegLine2D.makeFromInterpLine(prof, apprOpt);
                    case 'AxisymSegLineCart'
                        origFunc = [];
                        origFuncCoord = [];
                        segLine = surf.Profile;
                    otherwise
                        throw(MException('RingSet:NotSupportedSurface', 'Such surface cannot be represented as RingSet.'));
                end % switch
                
				obj.tree = [];
                rings = RingSet.segLine2DProfile2Rings(segLine);
                obj.r = ValueChecker.checkAndUpdateSize(rings.R, [inf 1]);
				obj.z = ValueChecker.checkAndUpdateSize(rings.z, [inf 1]);
				obj.k = ValueChecker.checkAndUpdateSize(rings.k, [inf 1]);
				obj.b = ValueChecker.checkAndUpdateSize(rings.b, [inf 1]);
				obj.nr = ValueChecker.checkAndUpdateSize(rings.nr, [inf 1]);
				obj.nz = ValueChecker.checkAndUpdateSize(rings.nz, [inf 1]);
                obj.origFunc = origFunc;
                obj.origFuncCoord = origFuncCoord;
                obj.opticalProperties = surf.Properties.clone();
                obj.cs = surf.CS.clone();
                obj.goodRings = zeros(obj.NumPrimitives, 4);
                for i = 1:obj.NumPrimitives
                    obj.goodRings(i, :) = [obj.r(i), obj.z(i), obj.r(i + 1), obj.z(i + 1)];
                end % for
			end % if
        end % function
        
        % destructor
        function delete(obj)
            if ~isempty(obj.tree)
                deleteRingsTreeMEX(obj.tree);
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
        %%%%%%%%%% IPrimitiveSet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [t, primInd, numInters] = getInter(obj, ray)
            if isempty(obj.tree)
                obj.tree = RingSet.makeRingsTree(obj.r, obj.z, obj.NumPrimitives);
            end % if
            ray = GeometryConvertor.getLocalRay(ray, obj.cs);
            [t, primInd, numInters] = getIntersRingsMEX(obj.tree, ray);
        end % function
        
        function numInters = getNumInters(obj, ray)
            ray = GeometryConvertor.getLocalRay(ray, obj.cs);
            numInters = getNumIntersRingsMEX(ray, obj.goodRings);
        end
     
        function res = getNormal(obj, ray, t, id)
            ValueChecker.checkNumerals(ray);
            ray = ValueChecker.checkAndUpdateSize(ray, [inf 7]);
            ValueChecker.checkConstraint(id, '>=', 1);
            ValueChecker.checkConstraint(id, '<=', obj.NumPrimitives);
            ValueChecker.checkConstraint([length(id) size(ray,1)], '==', length(t));
            
            ray = GeometryConvertor.getLocalRay(ray, obj.cs);
            
            if isempty(obj.origFunc)
                % segmented profile
                nrCur = obj.nr(id);
                nzCur = obj.nz(id);
            else
                % there is a possibility for exact computing normal
                intPoint = ray(:, 1:3) + ray(:, 4:6) .* repmat(t, [1 3]);
                switch obj.origFuncCoord
                    case 'polar'
                        % compute spherical coordinates of intersection point
                        [~, psi, ~] = cart2sph(   intPoint(:, 1), ...
                                                  intPoint(:, 2), ...
                                                  intPoint(:, 3));
                        psi = pi/2 - psi;
                        [r, r_psi] = obj.origFunc.eval(psi);
                        N = sqrt(r .^ 2 + r_psi .^ 2);
                        sinPsi = sin(psi);
                        cosPsi = cos(psi);
                        nrCur = (r .* sinPsi - r_psi .* cosPsi) ./ N;
                        nzCur = (r .* cosPsi + r_psi .* sinPsi) ./ N;
                    otherwise % 'cart'
                        r = sqrt(intPoint(:, 1) .^ 2 + intPoint(:, 2) .^ 2);
                        [~, z_r] = obj.origFunc.eval(r);
                        N = sqrt(z_r .^ 2 + 1);
                        nrCur = z_r ./ N;
                        nzCur = -1 ./ N;
                end % switch
            end % if
            
            newStartXY = ray(:,1:2) + ray(:,4:5) .* [t t];
            newStartL = sqrt(dot(newStartXY, newStartXY, 2));
            res = [  nrCur .* newStartXY(:,1) ./ newStartL, ...
                nrCur .* newStartXY(:,2) ./ newStartL, ...
                nzCur   ];
            res = GeometryConvertor.getGlobalVector(res, obj.cs);
        end % function
        
        function ray = getOutgoingRay(obj)
            ringNum = ceil(size(obj.r, 1) / 2);
            x = (obj.r(ringNum) + obj.r(ringNum + 1)) / 2;
            z = (obj.z(ringNum) + obj.z(ringNum + 1)) / 2;
            ray = [x 0 z 0 0 0 0];
            ray = GeometryConvertor.getGlobalRay(ray, obj.cs);
            ray(4:6) = obj.getNormal(ray, 0, ringNum);            
        end % function
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% ICloneable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function objNew = clone(obj)
            s = struct('className__', 'RingSet', 'r', obj.r, 'z', obj.z, 'k', obj.k, 'b', obj.b,...
                'nr', obj.nr, 'nz', obj.nz, 'origFunc', Serializator.save(obj.origFunc), ...
                'origFuncCoord', obj.origFuncCoord, 'cs', Serializator.save(obj.cs),...
                'opticalProperties', Serializator.save(obj.opticalProperties), ...
                'isRegistrator', obj.isRegistrator);
            objNew = RingSet(s);
        end % function
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% IDrawable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [pointsCloud, lines] = getPointsCloudAndLines(obj)
            phi = 0 : 2*pi / 199 : 2*pi;
            phi = repmat(phi, [length(obj.z) 1]);
            z = repmat(obj.z, [1 200]);
            x = repmat(obj.r, [1 200]) .* cos(phi);
            y = repmat(obj.r, [1 200]) .* sin(phi);
            defSize = size(x);
            
            x = reshape(x, [numel(x), 1]);
            y = reshape(y, [numel(y), 1]);
            z = reshape(z, [numel(z), 1]);            
            points = GeometryConvertor.getGlobalPoint([x y z], obj.cs);
            
            pointsCloud(:,:,1) = reshape(points(:,1), defSize);
            pointsCloud(:,:,2) = reshape(points(:,2), defSize);
            pointsCloud(:,:,3) = reshape(points(:,3), defSize);
            lines = {};
        end % function
        
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

        function rawTriangs = getTriangles(obj, apprOpt)
            segLine = SegLine2D([obj.r obj.z]);
            rawTriangs = TriangleSet.segLine2DProfile2RawTriangles(segLine, apprOpt);
        end % function
        
    end % methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Private %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Static, Access = private)
        
        % returns C-created binary tree
        function tree = makeRingsTree(R, z, numRings)
            if numRings > 10
                numRingsInLeaf = 5;
            else
                numRingsInLeaf = numRings;
            end
            
            kdDepth = ceil(log2(numRings / numRingsInLeaf));
            if kdDepth == 0
                kdDepth = 1;
            end
            
            goodRings = zeros(numRings, 4);
            for i = 1:numRings
                goodRings(i, :) = [R(i), z(i), R(i + 1), z(i + 1)];
            end % for
            
            ringsCount = size(goodRings, 1);
            aabb = zeros(ringsCount, 4);
            
            aabb(:, 1) = min([goodRings(:, 1) goodRings(:, 3)], [], 2); % min r
            aabb(:, 2) = min([goodRings(:, 2) goodRings(:, 4)], [], 2); % min z
            
            aabb(:, 3) = max([goodRings(:, 1) goodRings(:, 3)], [], 2); % max r
            aabb(:, 4) = max([goodRings(:, 2) goodRings(:, 4)], [], 2); % max z
            
            box = zeros(2, 2);
            
            % TODO for k dimensions
            minR = min(aabb(:, 1));
            minZ = min(aabb(:, 2));
            
            maxR = max(aabb(:, 3));
            maxZ = max(aabb(:, 4));
            
            box(1, :) = [minR minZ]; % leftBottom
            box(2, :) = [maxR maxZ]; % rightTop
            
            tree = buildRingsTreeMEX(goodRings, aabb, box, kdDepth);
        end % function
        
        % creates 'rings structure' for axisymmatrical surface with Oz
        % symmetry axis and profile defined by '2D segmented line'
        %   segLine         - '2D segmented line'
        function rings = segLine2DProfile2Rings(segLine)
            
            ValueChecker.checkClass(segLine, 'SegLine2D');
            
            R = segLine.Points(:, 1);
            z = segLine.Points(:, 2);
            
            k = (R(2:end) - R(1:end-1)) ./ (z(2:end) - z(1:end-1));
            b = R(1:end-1) - k .* z(1:end-1);
            
            % normals
            l = sqrt( (z(2:end) - z(1:end-1)).^2 + (R(2:end) - R(1:end-1)).^2 );
            nr = (z(2:end) - z(1:end-1)) ./ l;
            nz = (R(1:end-1) - R(2:end)) ./ l;
            
            rings.R = R;
            rings.z = z;
            rings.k = k;
            rings.b = b;
            rings.nr = nr;
            rings.nz = nz;
  
        end % function
        
    end % methods
    
end % classdef